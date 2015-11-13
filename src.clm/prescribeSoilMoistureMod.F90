  module prescribeSoilMoistureMod

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: prescribeSoilMoistureMod
  !
  ! !DESCRIPTION:
  ! assemble all code neccesary to prescribe SM in clm4.0
  !
  ! USES

    use shr_kind_mod,     only : r8 => shr_kind_r8
    use clm_varpar, only       : nlevgrnd, nlevsoi
    use abortutils,       only : endrun
    use clm_varcon,       only : istsoil, watmin
    use clm_varctl,       only : iulog
    use spmdMod,          only : masterproc
    use clm_time_manager, only : get_curr_date, get_step_size, get_perp_date, is_perpetual
  ! !PUBLIC TYPES:
    implicit none
    save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
    public :: prescribeSoilMoisture  ! Sets Soil Moisture

  !

  ! !PRIVATE MEMBER FUNCTIONS:
    private :: interpSoilMoisture   ! interpolate monthly SM data
    private :: readSoilMoisture     ! read SM data from file
  !

  ! !PRIVATE TYPES:

    integer , private :: TimeStep_old = 0 ! time step at last call
    real(r8), private :: dtime            ! land model time step (sec)
    real(r8), private :: timwt_soil(2)    ! time weights for month 1 and month 2
    real(r8), private, allocatable :: mh2osoi_liq2t(:,:,:) !  liquid soil water for interpolation (2 months) read from input files
    real(r8), private, allocatable :: mh2osoi_ice2t(:,:,:) !  frozen soil water for interpolation (2 months) read from input files

    logical, private            :: monthly            ! if .true. -> monthly, else daily
    character(len=256), private :: pSMfile            ! file name mit SM data to prescribe
    integer, private            :: pSMtype = 1        ! how to prescribe SM [default = 1] 
    real(r8), private           :: nudge = 1._r8      ! nudging
    integer, private            :: levstart = 1       ! start level prescribing SM
    integer, private            :: levstop = nlevsoi  ! end level prescribing SM
    logical, private            :: use_qdrai = .true. ! use qdrai in pSMtype3

  ! !REVISION HISTORY:
  ! Created by 
  !
  !EOP
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: prescribeSoilMoisture
  !
  ! !INTERFACE:
    subroutine prescribeSoilMoisture(lbc,ubc, num_nolakec, filter_nolakec)
  !
  ! !DESCRIPTION:
  ! loads SM from a netCDF file and overwrites the SM state
  ! (SOILLIQ (mh2osoi_liq) and SOILICE (mh2osoi_ice)) of the
  ! model
  ! Needs a history file originally created by CLM4.0
  ! add the following to the clm namelist:
  ! hist_fincl2 = 'SOILLIQ','SOILICE'
  ! hist_fincl3 = 'SOILLIQ','SOILICE'
  ! hist_fincl4 = 'SOILLIQ','SOILICE'
  ! hist_nhtfrq = 0, -24, -24, 0
  ! hist_mfilt = 1, 365, 365, 1
  ! hist_dov2xy = .true., .true., .false., .false.
  !
  !
  ! !USES:
      use clmtype
      use shr_const_mod,  only : SHR_CONST_TKFRZ ! freezing temperature of water
      use decompMod,      only : get_proc_bounds
      use clm_varcon,     only : spval
      use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)

  ! !ARGUMENTS:
      implicit none
      integer, intent(in) :: lbc, ubc                    ! column bounds
      integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
      integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
      integer, pointer    :: ltype(:)                    ! landunit type index


  ! !CALLED FROM:
  ! clm_driver
  !
  ! !REVISION HISTORY:
  ! Author: Mathias Hauser
  !
  ! !LOCAL VARIABLES:
  !
  ! local pointers to implicit in arguments
  !
    integer , pointer :: clandunit(:) ! column's landunit
    real(r8), pointer :: t_soisno(:,:)  ! soil temperature (Kelvin)

  ! local pointers to implicit inout arguments
    real(r8), pointer :: qflx_drain(:)    ! sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_irrig(:)    ! irrigation flux (mm H2O /s)
    real(r8), pointer :: qflx_surf(:)     ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_runoff(:)   ! total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(r8), pointer :: wa(:)             !water in the unconfined aquifer (mm)
  !
  ! local pointers to implicit out arguments
  !
      real(r8), pointer :: h2osoi_liq(:,:)     ! liquid soil water content level
      real(r8), pointer :: h2osoi_ice(:,:)     ! frozen soil water content level


  !
  ! !OTHER LOCAL VARIABLES:
  !EOP
  !
      integer :: fc,j,l,c           ! indices
      integer :: begc,endc          ! beg and end local column index
      integer :: ier                ! error code
      logical :: FirstCall = .true. ! make sure initPrescribeSoilMoisture is only called once


      real(r8) :: frac        ! fraction
      real(r8) :: SM          ! total SM for a timestep
      real(r8) :: water_avail ! available runoff
      real(r8) :: SMdeficit   ! missing SM per column
      real(r8) :: SMassigned  ! missing SM per column

  !-----------------------------------------------------------------------


      ! implicit out arguments
      h2osoi_liq  => cws%h2osoi_liq
      h2osoi_ice  => cws%h2osoi_ice
      
      ! implicit inout arguments
      qflx_drain        => cwf%qflx_drain
      qflx_irrig        => cwf%qflx_irrig
      qflx_surf         => cwf%qflx_surf
      qflx_runoff       => cwf%qflx_runoff
      wa                => cws%wa

      ! implicit in arguments
      t_soisno    => ces%t_soisno
      clandunit   => col%landunit
      ltype       => lun%itype

      call get_proc_bounds(begc=begc,endc=endc)
      
      if (.not. allocated(mh2osoi_liq2t) .or. .not. allocated(mh2osoi_ice2t)) then
        allocate (mh2osoi_liq2t(begc:endc,1:nlevgrnd,2), &
                  mh2osoi_ice2t(begc:endc,1:nlevgrnd,2), stat=ier)

        if (ier /= 0) then
          write(iulog,*) 'prescribeSoilMoistureMod allocation error'
          call endrun
        endif
        
        mh2osoi_liq2t(:,:,:) = nan
        mh2osoi_ice2t(:,:,:) = nan
      endif

      ! get file name and 'monthly' or 'daily' from namelist (only once)
      if (FirstCall) then
        call initPrescribeSoilMoisture ()
        FirstCall = .false.

        ! make sure irrig is 0 (and not nan)
        qflx_irrig(:) = 0._r8
      end if ! FirstCall

      ! get time weight and possibly a new time step
      call interpSoilMoisture()

      ! get time step
      dtime = get_step_size()

      ! sets timwt_soil, mh2osoi_liq2t, mh2osoi_ice2t

      ! overwrite the current soil water and ice content
      ! only if SOILLIQ & SOILICE are >= 0

! ===========================================================================================================================================================================

      ! CLASSICAL: prescribe SOILLIQ and SOILICE
      ! NOW INCLUDES NUDGING
      if (pSMtype == 1) then

        do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ltype(l) == istsoil) then
            ! Assign SOILLIQ and SOILICE
            do j = levstart, levstop
              if (mh2osoi_liq2t(c,j,1) .ge. 0_r8 .and. mh2osoi_ice2t(c,j,1) .ge. 0_r8 .and. mh2osoi_liq2t(c,j,2) .ge. 0_r8 .and. mh2osoi_ice2t(c,j,2) .ge. 0_r8) then
                h2osoi_liq(c,j) = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2) * nudge + (1._r8 - nudge) * h2osoi_liq(c,j)
                h2osoi_ice(c,j) = timwt_soil(1)*mh2osoi_ice2t(c,j,1) + timwt_soil(2)*mh2osoi_ice2t(c,j,2) * nudge + (1._r8 - nudge) * h2osoi_ice(c,j)
              end if ! liq >= 0
            end do ! j = 1, nlevgrnd
          end if
        end do ! fc = 1, num_nolakec

! ===========================================================================================================================================================================

      ! PRESCRIBE SOILLIQ if SOILICE == 0
      else if (pSMtype == 2) then
       
          do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            l = clandunit(c)
            if (ltype(l) == istsoil) then
              ! Assign SOILLIQ
              do j = levstart, levstop
                ! only overwrite liq if no ice is present
                if (mh2osoi_liq2t(c,j,1) .ge. 0_r8 .and. mh2osoi_liq2t(c,j,2) .ge. 0_r8 .and. h2osoi_ice(c,j) == 0) then
                  h2osoi_liq(c,j) = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
                end if
              end do 
            end if
          end do

! ===========================================================================================================================================================================

      ! USE only water from qflx_surf to PRESCRIBE SM
      else if (pSMtype == 3) then

        do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ltype(l) == istsoil) then

            if (use_qdrai) then
              ! total available water (ignoring qflx_qrgwl)
              water_avail = (qflx_surf(c) + qflx_drain(c)) * dtime
            else
              ! total available water (ignoring qflx_qrgwl)
              water_avail = qflx_surf(c) * dtime
            end if

            ! water available?
            if (water_avail .gt. 0._r8) then

              ! for subtracting water from qflx_surf and qflx_drain
              frac = (qflx_surf(c) * dtime) / (water_avail)

              if (masterproc) then
                  write(iulog,*) '-----------------------------'
                  write(iulog,*) 'Start New Gridpoint'
                  write(iulog,*) 'Max water, water_avail: ', water_avail
              end if
                        
              do j = levstart, levstop
              
                ! check if soil is frozen 
                ! (if layer 3 is frozen, only accumulate deficit for layer 1 & 2)
                if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                  if (masterproc) then
                    write(iulog,*) 'there is ice on level: ', j
                  end if

                  exit ! leave do j = levstart, levstop
                end if

                ! desired SM state at this level
                SM = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
                
                ! TODO ADD SOILICE?

                ! SMdeficit
                SMdeficit = max(SM - h2osoi_liq(c,j), 0._r8)
                
                ! is there a deficit?
                if (SMdeficit .gt. 0._r8) then
                  
                  ! how much can we add?
                  SMassigned = min(water_avail, SMdeficit)

                  ! add water
                  h2osoi_liq(c,j) = h2osoi_liq(c,j) + SMassigned

                  ! subtract from runoff (protect from rounding errors)

                  ! according to the fraction available water
                  qflx_surf(c) = max(qflx_surf(c) - frac * (SMassigned / dtime), 0._r8)
                  ! if not use_qdrai: frac == 1
                  qflx_drain(c) = max(qflx_drain(c) - (1._r8 - frac) * (SMassigned / dtime), 0._r8)
                  
                  ! qflx_runoff = qflx_drain+qflx_surf+qflx_qrgwl
                  qflx_runoff(c) = max(qflx_runoff(c) - SMassigned / dtime, 0._r8)
                  
                  ! also decrease available water
                  water_avail = water_avail - SMassigned
                  
                  if (masterproc) then
                    write(iulog,*) 'Level: ', j, 'SMdeficit: ', SMdeficit, 'SMassigned: ', SMassigned
                  end if


                  ! check if water is available
                  if (water_avail .eq. 0._r8) then
                    if (masterproc) then
                      write(iulog,*) 'there is no more water on level: ', j
                    end if

                    exit
                  end if              
                end if ! (SMdeficit .gt. 0._r8)
              end do ! j = 1, nlevgrnd
              ! restrict the irrigation to the surface runoff of this time step

            end if ! not spval and > 0
          end if ! istsoil
        end do ! fc = 1, num_nolakec

! ===========================================================================================================================================================================

      ! ONLY PRESCRIBE IF LESS THAN 
      else if (pSMtype == 4) then

        do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ltype(l) == istsoil) then
            ! Assign SOILLIQ and SOILICE
            do j = levstart, levstop
              if (mh2osoi_liq2t(c,j,1) .ge. 0._r8 .and. mh2osoi_ice2t(c,j,1) .ge. 0._r8 .and. mh2osoi_liq2t(c,j,2) .ge. 0._r8 .and. mh2osoi_ice2t(c,j,2) .ge. 0._r8) then
                h2osoi_liq(c,j) = max(timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2), h2osoi_liq(c,j))
                h2osoi_ice(c,j) = max(timwt_soil(1)*mh2osoi_ice2t(c,j,1) + timwt_soil(2)*mh2osoi_ice2t(c,j,2), h2osoi_ice(c,j))
              end if ! liq >= 0
            end do ! j = 1, nlevgrnd
          end if
        end do ! fc = 1, num_nolakec

! ===========================================================================================================================================================================

      ! USE only water from qflx_surf and qflx_drain (sub- and surface runoff) 

      else if (pSMtype == 5) then
        call endrun(trim(subname)//' pSMtype 5 does not exist')



! ===========================================================================================================================================================================


      ! IRRIGATION 
      ! we use the runoff from this time step and add it to the prec flux
      ! (avoiding interception) at the ground in the next timestep

      else if (pSMtype == 6) then
        call endrun(trim(subname)//' something is wrong with pSMtype 6 (dont know what)')

        do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ltype(l) == istsoil) then

            ! make sure irrigation is not carried over
            qflx_irrig(c) = 0._r8
            if (qflx_surf(c) .gt. 0._r8 .and. qflx_surf(c) /= spval) then

              ! determine irrigation flux
              do j = levstart, levstop
                
                ! check if soil is frozen 
                ! (if layer 3 is frozen, only accumulate deficit for layer 1 & 2)
                if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then

                  if (masterproc) then
                    write(iulog,*) 'there is ice on level: ', j
                  end if

                  exit
                end if

                ! desired SM state at this level
                SM = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
                
                ! desired irrigation
                qflx_irrig(c) = qflx_irrig(c) + max(SM - h2osoi_liq(c,j), 0._r8) / dtime
                  

              end do ! j = 1, nlevgrnd
              ! restrict the irrigation to the surface runoff of this time step

              if (masterproc) then
                write(iulog,*) 'qflx_irrig(c) ',  qflx_irrig(c)
                write(iulog,*) 'qflx_surf(c) ',  qflx_surf(c)
              end if

              ! restrict irrigation to available water
              qflx_irrig(c) = max(qflx_surf(c), qflx_irrig(c))

              ! correct the two runoff terms
              qflx_surf(c) = qflx_surf(c) - qflx_irrig(c)
              qflx_runoff(c) = qflx_runoff(c) - qflx_irrig(c)

            end if ! gt 0 and not spval
          end if ! istsoil
        end do ! fc = 1, num_nolakec


      end if ! pSMtype


    end subroutine prescribeSoilMoisture




  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: interpMonthlySoil
  !
  ! !INTERFACE:
    subroutine initPrescribeSoilMoisture ()
  !
  ! !DESCRIPTION:
  ! Set monthly or daily SM, define file with SM data
  ! goal -> use namelist
  !
  ! !USES:
    use fileutils,      only : getfil, getavu, relavu
  ! !ARGUMENTS:
      implicit none
  !
  ! !REVISION HISTORY:
  ! Created by Mathias Hauser
  !
  !
  ! !LOCAL VARIABLES:
      integer :: unitn
      integer :: ierr 
      ! character(len=256) :: type_name    ! name of prescribe                   
      character(len=256) :: locfn        ! local file name
      character(len=256) :: subname      ! name of routine

  !EOP

      subname = 'initPrescribeSoilMoisture'

      ! Input datasets
      namelist /prescribe_sm/  &
           pSMfile, monthly, pSMtype, nudge, levstart, levstop, use_qdrai

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! CAREFUL: IF YOU USE A DAILY pSMfile BUT SET monthly = .true.  !
      ! THE MODEL WILL STILL RUN!                                     !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call getfil(trim('prescribe_SM_nl'), locfn, 0)
      
      if (masterproc) then
        write(iulog,*) 'local file name: ', trim(locfn)
      endif ! masterproc

      unitn = getavu()
     
      open(unitn, file=trim(locfn), status='old')
      ierr = 1
      do while (ierr /= 0)
          read(unitn, prescribe_sm, iostat=ierr)
          if (ierr < 0) then
              call endrun(trim(subname)//' encountered end-of-file on prescribe_SM_nl read')
          endif
      end do
      call relavu(unitn)


      if (masterproc) then
        write(iulog,*) 'Read pSMfile from namelist:', trim(pSMfile)
        write(iulog,*) 'Read monthly from namelist:', monthly
        write(iulog,*) 'Read pSMtype from namelist:', pSMtype
        write(iulog,*) 'Read nudge from namelist:', nudge
        write(iulog,*) 'Read levstart from namelist:', levstart
        write(iulog,*) 'Read levstop from namelist:', levstop
        write(iulog,*) 'Read use_qdrai from namelist:', use_qdrai
  
        if (nudge .lt. 0._r8 .or. nudge .gt. 1._r8) then
          call endrun(trim(subname)//'nudge must be in 0..1!')
        end if


        if (levstart .gt. levstop) then
          call endrun(trim(subname)//'levstop must be bigger than levend')
        end if

        if (max(levstart,  levstop) .gt. nlevsoi) then
          call endrun(trim(subname)//'levstop/ start must not exceed nlevsoi (10)')
        end if

        if (min(levstart,  levstop) .lt. 1) then
          call endrun(trim(subname)//'levstop/ start must not be smaller than 1')
        end if





        if (pSMtype == 1) then
          write(iulog,*) 'pSMtype: prescribe both (DEFAULT)'
          write(iulog,*) 'SOILICE AND SOILLIQ is prescribed, whatever the conditions'
          write(iulog,*) 'for any gridcell, level, timestep where LIQ !OR! ICE is < 0,'
          write(iulog,*) 'neither is prescribed (if only one is prescribed it yields unphysical values).'

        else if (pSMtype == 2) then
          write(iulog,*) 'pSMtype: SOILLIQ if no SOILICE'
          write(iulog,*) 'ONLY SOILLIQ is prescribed. SOILICE is calculated.'
          write(iulog,*) '-> can NOT prescribe SOILICE (is calculated interactively)'
          write(iulog,*) '-> if SOILICE is PRESENT -> SOILLIQ is INTERACTIVE'
        
        else if (pSMtype == 3) then
          ! write(iulog,*) 'pSMtype: FRACTION'
          ! write(iulog,*) 'At every timestep the sotal SM (SOILLIQ + SOILICE) is'
          ! write(iulog,*) 'prescribed. The fraction that becomes SOILLIQ and '
          ! write(iulog,*) 'SOILICE. Is determined by the current ICE and LIQ '
          ! write(iulog,*) 'fraction.'
          ! write(iulog,*) 'fr = SOILLIQ^n / (SOILLIQ^n + SOILICE^n), where n = TimeStep'
          ! write(iulog,*) 'SOILLIQ^n+1 = fr * (ICE + LIQ)'
          ! write(iulog,*) 'SOILICE^n+1 = (1 - fr) * (ICE + LIQ)'
          ! write(iulog,*) 'where LIQ and ICE are read from the file'

        else
          !call endrun(trim(subname)//'pSMtype namelist option must be one of 1, 2,3!')
        endif ! pSMtype


      endif ! masterproc

    end subroutine initPrescribeSoilMoisture



  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: interpSoilMoisture
  !
  ! !INTERFACE:
    subroutine interpSoilMoisture ()
  !
  ! !DESCRIPTION:
  ! Determine if 2 new months or days of data are to be read.
  ! if yes -> 
  ! !USES:
      use clm_varctl,       only : fsurdat
  !
  ! !ARGUMENTS:
      implicit none
  !
  ! !REVISION HISTORY:
  ! Created by Ruth Lorenz
  ! Mathias Hauser: added interpolation for daily SM input
  !
  !EOP
  !
  ! LOCAL VARIABLES:
      integer :: kyr         ! year (0, ...) for nstep+1
      integer :: kmo         ! month (1, ..., 12)
      integer :: kda         ! day of month (1, ..., 31)
      integer :: ksec        ! seconds into current date for nstep+1
      real(r8):: t           ! a fraction: kda/ndaypm
      integer :: it(2)       ! month 1 and month 2 (step 1)
      integer :: TimeStep(2) ! months to be interpolated (1 to 12)
      integer :: doy         ! day of year (1..365)    

      real(r8) :: nsec = 86400._r8 ! num of sec per day
      integer, dimension(12) :: ndaypm = &
           (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/) !days per month
      integer, dimension(12) :: cdaypm = &
          (/0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334/) ! cumulative
          ! number of days per month

  !-----------------------------------------------------------------------
      dtime = get_step_size()

      if ( is_perpetual() ) then
         call get_perp_date(kyr, kmo, kda, ksec, offset=int(dtime))
      else
         call get_curr_date(kyr, kmo, kda, ksec, offset=int(dtime))
      end if

      if (monthly) then ! interpolate monthly data

          t = (kda-0.5_r8) / ndaypm(kmo)
          it(1) = t + 0.5_r8
          it(2) = it(1) + 1
          TimeStep(1) = kmo + it(1) - 1
          TimeStep(2) = kmo + it(2) - 1
          if (TimeStep(1) <  1) TimeStep(1) = 12
          if (TimeStep(2) > 12) TimeStep(2) = 1
          timwt_soil(1) = (it(1)+0.5_r8) - t
          timwt_soil(2) = 1._r8-timwt_soil(1)

      else ! interpolate daily data
          doy = cdaypm(kmo) + kda

          ! fraction of day that has passed
          t = ksec / nsec

          ! if t < 0.5 we need 'doy - 1' and 'doy'
          ! else we need 'doy' and 'doy + 1'
          TimeStep(1) = doy + floor(t - 0.5)
          TimeStep(2) = doy + floor(t + 0.5)

          if (TimeStep(1) < 1) TimeStep(1) = 365
          if (TimeStep(2) > 365) TimeStep(2) = 1

          timwt_soil(1) = 1._r8 - abs(t - 0.5_r8)
          timwt_soil(2) = 1._r8 - timwt_soil(1)

      endif ! monthly

      
      if (TimeStep_old /= TimeStep(1)) then
        call readSoilMoisture (kmo, kda, TimeStep)
        TimeStep_old = TimeStep(1)
      end if

    end subroutine interpSoilMoisture





  !same as ReadMonthlyVegetation except for Soil Moisture
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: readSoilMoisture
  !
  ! !INTERFACE:
    subroutine readSoilMoisture (kmo, kda, months_soil)
  !
  ! !DESCRIPTION:
  ! Read monthly soil moisture data for two consec. months or days
  !
  ! !USES:
      use clmtype
      use ncdio_pio
      use netcdf
      use shr_kind_mod, only : r8 => shr_kind_r8
      use decompMod,    only : get_proc_bounds
      use clm_varpar,   only : nlevgrnd
      use clm_varcon,   only : istsoil
      use fileutils,    only : getfil

      
  ! !ARGUMENTS:
      implicit none

      integer, intent(in) :: kmo            ! month (1, ..., 12)
      integer, intent(in) :: kda            ! day of month (1, ..., 31)
      integer, intent(in) :: months_soil(2) ! months to be interpolated (1 to 12)
  !
  ! !REVISION HISTORY:
  ! Created by Ruth Lorenz
  ! Mathias Hauser: Rewritten for CLM 1.2.1
  !
  ! !LOCAL VARIABLES:
  !EOP
      integer :: n,j,k,l,m,c                ! indices
      character(len=256) :: locfn           ! local file name
      integer , pointer  :: clandunit(:)    ! column's landunit
      integer , pointer  :: ltype(:)        ! landunit type index
      type(file_desc_t)  :: ncid            ! netcdf id
      integer :: dimid,varid                ! input netCDF id's
      integer :: ntim                       ! number of input data time samples
      integer :: ncolumn_i                  ! number of columns
      integer :: begc,endc                  ! beg and end local c index
      integer :: ier,ret                    ! error code
      integer :: closelatidx,closelonidx
      real(r8):: closelat,closelon
      logical :: readvar

      real(r8), pointer :: mh2osoi_liq(:,:)  ! liquid soil water content read from input file
      real(r8), pointer :: mh2osoi_ice(:,:)  ! frozen soil water content read from input file
      character(len=32) :: subname = 'readSoilMoisture'
      character(len=32) :: cTimeStep ! string for the choosen time step (daily, monthly)
  !-----------------------------------------------------------------------


      ! Assign local pointers to derived subtypes components (landunit-level)
      ltype       => lun%itype
      ! Assign local pointers to derived subtypes components (column-level)
      clandunit   => col%landunit
      
      ! Determine necessary indices
      call get_proc_bounds(begc=begc, endc=endc)

      allocate(mh2osoi_liq(begc:endc, 1:nlevgrnd), &
               mh2osoi_ice(begc:endc, 1:nlevgrnd), stat=ier)

      if (ier /= 0) then
         write(iulog, *) subname, 'allocation big error '
         call endrun()
      end if

      if (monthly) then
        cTimeStep = 'month'
      else
        cTimeStep = 'day'
      endif

      ! ----------------------------------------------------------------------
      ! Open monthly soil moisture file
      ! Read data from column
      ! ----------------------------------------------------------------------
      
      if (masterproc) then
        write(iulog,*) 'Attempting to read ', trim(cTimeStep), 'ly soil moisture data...'
        write(iulog,*) 'month = ', kmo, ' day = ', kda
      endif ! masterproc


      do k=1,2  ! loop over months/ days and read SM data

          ! get file
          call getfil(trim(pSMfile), locfn, 0)
          call ncd_pio_openfile (ncid, trim(locfn), 0)


          if (masterproc) then
              write(iulog,*) 'Before read '
              write(iulog,*) trim(cTimeStep), '(k) ', months_soil(k)
          endif ! masterproc


          call ncd_io(ncid=ncid, varname='SOILLIQ', flag='read', data=mh2osoi_liq, dim1name=namec, &
            nt=months_soil(k), readvar=readvar)
          if (.not. readvar) call endrun(trim(subname) // ' ERROR: SOILLIQ NOT on pSMfile' // trim(pSMfile))
          

          call ncd_io(ncid=ncid, varname='SOILICE', flag='read', data=mh2osoi_ice, dim1name=namec, &
            nt=months_soil(k), readvar=readvar)
          if (.not. readvar) call endrun(trim(subname) // ' ERROR: SOILICE NOT on pSMfile ' // trim(pSMfile))     

        call ncd_pio_closefile(ncid)


        if (masterproc) then
            write(iulog,*) 'Successfully read ',  trim(cTimeStep), 'ly soil moisture data for'
            write(iulog,*) cTimeStep, ' ', months_soil(k)
        end if


         ! store data directly in clmtype structure
         do c = begc, endc
            l = clandunit(c)
            if (ltype(l) == istsoil) then 
               do j = 1, nlevsoi
                  mh2osoi_liq2t(c, j, k) = mh2osoi_liq(c, j)
                  mh2osoi_ice2t(c, j, k) = mh2osoi_ice(c, j)
               end do
            end if
         end do   ! end of loop over columns

      end do   ! end of loop over months
      
      if (masterproc) then
        write(iulog,*) 'Successfully read soil moisture data'
      endif ! masterproc

      deallocate(mh2osoi_liq, mh2osoi_ice)
    end subroutine readSoilMoisture




  end module prescribeSoilMoistureMod
