module prescribeSoilMoistureMod

  !---------------------------------------------------------------------
  ! BOP
  !
  ! !MODULE: prescribeSoilMoistureMod
  !
  ! !DESCRIPTION:
  ! assemble all code neccesary to prescribe SM in clm4.0
  !
  ! USES

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use clm_varpar,       only : nlevgrnd, nlevsoi
  use abortutils,       only : endrun
  use clm_varcon,       only : istsoil, watmin
  use clm_varctl,       only : iulog
  use spmdMod,          only : masterproc
  use clm_time_manager, only : get_curr_date, get_step_size, get_perp_date, is_perpetual

  !PUBLIC TYPES:
  implicit none
  save

  !PUBLIC MEMBER FUNCTIONS:
  public :: prescribeSoilMoisture  ! Sets Soil Moisture

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: interpSoilMoisture   ! interpolate monthly SM data
  private :: readSoilMoisture     ! read SM data from file

  ! PRIVATE TYPES:
  integer , private :: TimeStep_old = 0 ! time step at last call
  real(r8), private :: dtime            ! land model time step (sec)
  real(r8), private :: timwt_soil(2)    ! time weights for month 1 and month 2
  real(r8), private, allocatable :: reservoir(:)           ! reservoir
  real(r8), private, allocatable :: mh2osoi_liq2t(:, :, :) ! liquid soil water for interpolation (2 months) read from input files
  real(r8), private, allocatable :: mh2osoi_ice2t(:, :, :) ! frozen soil water for interpolation (2 months) read from input files

  logical, private            :: monthly                         ! if .true. -> monthly, else daily
  logical, private            :: interp_day = .true.             ! time interpolation when monthly == false
  character(len=256), private :: pSMfile                         ! file name mit SM data to prescribe
  integer, private            :: pSMtype = 1                     ! how to prescribe SM [default = 1] 
  real(r8), private           :: reservoir_capacity = 0._r8      ! capacity of the reservoir storing water for dry periods
  integer, private            :: levstart = 1                    ! start level prescribing SM
  integer, private            :: levstop = nlevsoi               ! end level prescribing SM
  logical, private            :: use_qdrai = .true.              ! use qdrai in pSMtype3
  logical, private            :: one_file_per_timestep = .true.  ! use one file for each timestep

  ! REVISION HISTORY:
  ! Created by 
  !
  ! EOP
  !---------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: prescribeSoilMoisture
  !
  ! !INTERFACE:
    subroutine prescribeSoilMoisture(lbc, ubc, num_nolakec, filter_nolakec)
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

      ! ARGUMENTS:
      implicit none
      integer, intent(in) :: lbc, ubc                  ! column bounds
      integer, intent(in) :: num_nolakec               ! number of column non-lake points in column filter
      integer, intent(in) :: filter_nolakec(ubc-lbc+1) ! column filter for non-lake points
      integer, pointer    :: ltype(:)                  ! landunit type index


      ! CALLED FROM:
      ! clm_driver
      !
      ! !REVISION HISTORY:
      ! Author: Mathias Hauser
      !
      ! LOCAL VARIABLES:

      ! local pointers to implicit in arguments
      integer , pointer :: clandunit(:)   ! column's landunit
      real(r8), pointer :: t_soisno(:, :) ! soil temperature (Kelvin)

      ! local pointers to implicit inout arguments
      real(r8), pointer :: qflx_drain(:)  ! sub-surface runoff (mm H2O /s)
      real(r8), pointer :: qflx_irrig(:)  ! irrigation flux (mm H2O /s)
      real(r8), pointer :: qflx_surf(:)   ! surface runoff (mm H2O /s)
      real(r8), pointer :: qflx_runoff(:) ! total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
      real(r8), pointer :: wa(:)          ! water in the unconfined aquifer (mm)
      
      real(r8), pointer :: reservoir(:)   ! water in the reservoir (mm)

      real(r8), pointer :: soilliq_prescribed(:, :) ! how much SL did we prescribe at this time step
      real(r8), pointer :: soilice_prescribed(:, :) ! how much SI did we prescribe at this time step



      ! local pointers to implicit out arguments
      real(r8), pointer :: h2osoi_liq(:,:) ! liquid soil water content level
      real(r8), pointer :: h2osoi_ice(:,:) ! frozen soil water content level


      !
      ! !OTHER LOCAL VARIABLES:
      !EOP
      !
      integer :: fc,j,l,c           ! indices
      integer :: begc,endc          ! beg and end local column index
      integer :: ier                ! error code
      logical :: FirstCall = .true. ! make sure initPrescribeSoilMoisture is only called once

      real(r8) :: frac        ! fraction
      real(r8) :: SL          ! total SOILLIQ for a timestep
      real(r8) :: SI          ! total SOILICE for a timestep
      real(r8) :: SM          ! total SOILMOISTURE (= SL + SI) for a timestep
      real(r8) :: water_avail ! available runoff
      real(r8) :: SMdeficit   ! missing SM per column
      real(r8) :: SMassigned  ! missing SM per column


      ! ---------------------------------------------------------------------


      ! implicit out arguments
      h2osoi_liq  => cws%h2osoi_liq
      h2osoi_ice  => cws%h2osoi_ice
      
      ! implicit inout arguments
      qflx_drain        => cwf%qflx_drain
      qflx_irrig        => cwf%qflx_irrig
      qflx_surf         => cwf%qflx_surf
      qflx_runoff       => cwf%qflx_runoff
      wa                => cws%wa



      reservoir          => cwf%reservoir

      soilliq_prescribed => cws%soilliq_prescribed
      soilice_prescribed => cws%soilice_prescribed
      
      ! implicit in arguments
      t_soisno    => ces%t_soisno
      clandunit   => col%landunit
      ltype       => lun%itype

      call get_proc_bounds(begc=begc,endc=endc)
      
      ! allocate mh2osoi_liq2t & mh2osoi_ice2t
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
      endif ! FirstCall

      ! get time weight and possibly a new time step
      call interpSoilMoisture()

      ! get time step
      dtime = get_step_size()

      ! sets timwt_soil, mh2osoi_liq2t, mh2osoi_ice2t

      ! overwrite the current soil water and ice content
      ! only if SOILLIQ & SOILICE are >= 0

! =======================================================================================================================

      ! FIND OUT HOW MUCH LIQ AND ICE IS PRESCRIBED
      do fc = 1, num_nolakec
        c = filter_nolakec(fc)
        l = clandunit(c)
        if (ltype(l) == istsoil) then
          do j = 1, nlevsoi
              ! save current state of h2osoi liq and ice
              soilliq_prescribed(c, j) = h2osoi_liq(c,j)
              soilice_prescribed(c, j) = h2osoi_ice(c,j)
          end do ! j = 1, nlevgrnd
        endif
      end do ! fc = 1, num_nolakec

! =======================================================================================================================

      ! CLASSICAL: prescribe SOILLIQ and SOILICE
      if (pSMtype == 1) then

        do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ltype(l) == istsoil) then
            ! Assign SOILLIQ and SOILICE
            do j = levstart, levstop
              ! only prescribe if all SL and SI values are >= 0
              if (min(mh2osoi_liq2t(c,j,1), mh2osoi_ice2t(c,j,1), mh2osoi_liq2t(c,j,2), mh2osoi_ice2t(c,j,2)) .ge. 0_r8) then

                ! obtain desired SL and SI at this timestep
                SL = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
                SI = timwt_soil(1)*mh2osoi_ice2t(c,j,1) + timwt_soil(2)*mh2osoi_ice2t(c,j,2)

                h2osoi_liq(c,j) = SL
                h2osoi_ice(c,j) = SI
              endif ! liq >= 0
            end do ! j = 1, nlevgrnd
          endif
        end do ! fc = 1, num_nolakec

! =======================================================================================================================

      ! PRESCRIBE SM if T_SOIL > 273.15
      else if (pSMtype == 2) then
       
        do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ltype(l) == istsoil) then
            do j = levstart, levstop

              ! check if soil is frozen 
              ! (if layer 3 is frozen, only accumulate deficit for layer 1 & 2)
              ! only overwrite liq if no ice is present
              if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                exit ! leave do j = levstart, levstop
              endif
              
              ! only prescribe if all SL and SI values are >= 0
              if (min(mh2osoi_liq2t(c,j,1), mh2osoi_ice2t(c,j,1), mh2osoi_liq2t(c,j,2), mh2osoi_ice2t(c,j,2)) .ge. 0_r8) then

                ! obtain desired SL and SI at this timestep
                SL = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
                SI = timwt_soil(1)*mh2osoi_ice2t(c,j,1) + timwt_soil(2)*mh2osoi_ice2t(c,j,2)

                ! prescribe total water (at places where there is no more soil ice)
                SM = SL + SI

                h2osoi_liq(c,j) = SM
              endif
            end do 
          endif
        end do

! =======================================================================================================================

      ! USE only water from qflx_surf and maybe qflx_drai to PRESCRIBE SM
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
              water_avail = (qflx_surf(c)) * dtime
            endif

            ! water available? (stop div0 errors)
            if (water_avail .gt. 0._r8) then
              ! for subtracting water from qflx_surf and qflx_drain
              frac = (qflx_surf(c) * dtime) / water_avail
            endif

            if (masterproc) then
                write(iulog,*) '-----------------------------'
                write(iulog,*) 'Start New Gridpoint'
                write(iulog,*) 'Max water, water_avail: ', water_avail
                write(iulog,*) 'Reservoir: ', reservoir(c)
            endif
                      
            do j = levstart, levstop
            
              ! check if soil is frozen 
              ! (if layer 3 is frozen, only 'correct' deficit for layer 1 & 2)
              if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                if (masterproc) then
                  write(iulog,*) 'there is ice on level: ', j
                endif

                exit ! leave do j = levstart, levstop
              endif

              if (min(mh2osoi_liq2t(c,j,1), mh2osoi_ice2t(c,j,1), mh2osoi_liq2t(c,j,2), mh2osoi_ice2t(c,j,2)) .ge. 0_r8) then

                ! obtain desired SL and SI
                SL = timwt_soil(1) * mh2osoi_liq2t(c, j, 1) + timwt_soil(2) * mh2osoi_liq2t(c, j, 2)
                SI = timwt_soil(1) * mh2osoi_ice2t(c, j, 1) + timwt_soil(2) * mh2osoi_ice2t(c, j, 2)

                ! prescribe total water (important at places where there is no more soil ice)
                SM = SL + SI

                ! SMdeficit
                SMdeficit = max(SM - h2osoi_liq(c, j), 0._r8)
                
                ! is there a deficit?
                if (SMdeficit .gt. 0._r8) then
                  
                  ! RUNOFF
                  if (water_avail .gt. 0._r8) then
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
                  endif ! (water_avail .gt. 0._r8)

                  ! RESERVOIR
                  if (reservoir(c) > 0._r8) then
                    ! get an update of the missing water
                    SMdeficit = max(0._r8, SMdeficit - SMassigned)

                    ! how much can we add?
                    SMassigned = min(reservoir(c), SMdeficit)

                    ! add water
                    h2osoi_liq(c,j) = h2osoi_liq(c,j) + SMassigned                

                    ! update the reservoir
                    reservoir(c) = max(0._r8, reservoir(c) - SMassigned)
                  endif ! reservoir(c) .gt. 0._r8

                  ! check if water is available
                  if (water_avail <= 0._r8 .and. reservoir(c) <= 0._r8) then
                    if (masterproc) then
                      write(iulog,*) 'there is no more water on level: ', j
                    endif

                    exit
                  endif              
                endif ! (SMdeficit .gt. 0._r8)
              endif ! (min(mh2osoi_liq2t(c,j,1), ...) .ge. 0_r8) then
            end do ! j = 1, nlevgrnd


            ! Fill the Reservoir if neccesary and possible
            SMdeficit = max(reservoir_capacity - reservoir(c),  0._r8)
            if (min(SMdeficit, water_avail) > 0._r8) then

              ! how much water is missing in the reservoir                  
              SMassigned = min(water_avail, SMdeficit)

              ! add water
              reservoir(c) = reservoir(c) + SMassigned

              ! subtract from runoff (protect from rounding errors)

              ! according to the fraction available water
              qflx_surf(c) = max(qflx_surf(c) - frac * (SMassigned / dtime), 0._r8)
              ! if not use_qdrai: frac == 1
              qflx_drain(c) = max(qflx_drain(c) - (1._r8 - frac) * (SMassigned / dtime), 0._r8)
              
              ! qflx_runoff = qflx_drain+qflx_surf+qflx_qrgwl
              qflx_runoff(c) = max(qflx_runoff(c) - SMassigned / dtime, 0._r8)

            endif

          endif ! istsoil
        end do ! fc = 1, num_nolakec

! =======================================================================================================================

      ! PRESCRIBE SOILLIQ if T_SOIL > 273.15
      ! ONLY PRESCRIBE IF LESS THAN
      else if (pSMtype == 4) then

        do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ltype(l) == istsoil) then
            do j = levstart, levstop

              ! check if soil is frozen 
              ! (if layer 3 is frozen, only accumulate deficit for layer 1 & 2)
              ! only overwrite liq if no ice is present
              if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                exit ! leave do j = levstart, levstop
              endif
              
              ! only prescribe if all SL and SI values are >= 0
              if (min(mh2osoi_liq2t(c,j,1), mh2osoi_ice2t(c,j,1), mh2osoi_liq2t(c,j,2), mh2osoi_ice2t(c,j,2)) .ge. 0_r8) then

                ! obtain desired SL and SI at this timestep
                SL = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
                SI = timwt_soil(1)*mh2osoi_ice2t(c,j,1) + timwt_soil(2)*mh2osoi_ice2t(c,j,2)

                ! prescribe total water (at places where there is no more soil ice)
                SM = SL + SI

                SM = SM
                ! ONLY PRESCRIBE IF LESS THAN
                h2osoi_liq(c,j) = max(h2osoi_liq(c,j), SM)
              endif
            end do 
          endif
        end do

! =======================================================================================================================

      else if (pSMtype == 5) then
      ! FRACTION of LIQ and ICE

        do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ltype(l) == istsoil) then
            do j = levstart, levstop
              
              ! only prescribe if all SL and SI values are >= 0
              if (min(mh2osoi_liq2t(c,j,1), mh2osoi_ice2t(c,j,1), mh2osoi_liq2t(c,j,2), mh2osoi_ice2t(c,j,2)) .ge. 0_r8) then

                ! obtain desired SL and SI at this timestep
                SL = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
                SI = timwt_soil(1)*mh2osoi_ice2t(c,j,1) + timwt_soil(2)*mh2osoi_ice2t(c,j,2)

                ! prescribe total water
                SM = SL + SI

                if (h2osoi_liq(c,j) + h2osoi_ice(c,j) .eq. 0) then
                  frac = 1._r8 ! all to water
                else
                  frac = h2osoi_liq(c,j) / (h2osoi_liq(c,j) + h2osoi_ice(c,j))
                endif

                frac = min(max(frac, 0._r8), 1._r8)

                h2osoi_liq(c,j) = frac * SM
                h2osoi_ice(c,j) = (1._r8 - frac) * SM

                if (masterproc) then
                  write(iulog,*) 'frac: ', frac, 'SM: ', SM
                endif


              endif
            end do 
          endif
        end do

! =======================================================================================================================


      ! IRRIGATION 
      ! we use the runoff from this time step and add it to the prec flux
      ! (avoiding interception) at the ground in the next timestep

      else if (pSMtype == 6) then

        call endrun('initPrescribeSoilMoisture something is wrong with pSMtype 6 (dont know what)')

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
                  endif

                  exit
                endif

                ! desired SM state at this level
                SM = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
                
                ! desired irrigation
                qflx_irrig(c) = qflx_irrig(c) + max(SM - h2osoi_liq(c,j), 0._r8) / dtime
                  

              end do ! j = 1, nlevgrnd
              ! restrict the irrigation to the surface runoff of this time step

              if (masterproc) then
                write(iulog,*) 'qflx_irrig(c) ',  qflx_irrig(c)
                write(iulog,*) 'qflx_surf(c) ',  qflx_surf(c)
              endif

              ! restrict irrigation to available water
              qflx_irrig(c) = max(qflx_surf(c), qflx_irrig(c))

              ! correct the two runoff terms
              qflx_surf(c) = qflx_surf(c) - qflx_irrig(c)
              qflx_runoff(c) = qflx_runoff(c) - qflx_irrig(c)

            endif ! gt 0 and not spval
          endif ! istsoil
        end do ! fc = 1, num_nolakec


      endif ! pSMtype


! =======================================================================================================================

      ! FIND OUT HOW MUCH LIQ AND ICE IS PRESCRIBED
      do fc = 1, num_nolakec
        c = filter_nolakec(fc)
        l = clandunit(c)
        if (ltype(l) == istsoil) then
          do j = 1, nlevsoi
              ! subtract the old from the new SL and SI state
              soilliq_prescribed(c, j) = h2osoi_liq(c,j) - soilliq_prescribed(c, j)
              soilice_prescribed(c, j) = h2osoi_ice(c,j) - soilice_prescribed(c, j)
          end do ! j = 1, nlevgrnd
        endif
      end do ! fc = 1, num_nolakec

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
           pSMfile, monthly, pSMtype, reservoir_capacity, levstart, levstop, use_qdrai, interp_day, one_file_per_timestep

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
        write(iulog,*) 'Read "pSMfile" from namelist:', trim(pSMfile)
        write(iulog,*) 'Read "monthly" from namelist:', monthly
        write(iulog,*) 'Read "interp_day" from namelist:', interp_day
        write(iulog,*) 'Read "pSMtype" from namelist:', pSMtype
        write(iulog,*) 'Read "reservoir_capacity" from namelist:', reservoir_capacity
        write(iulog,*) 'Read "levstart" from namelist:', levstart
        write(iulog,*) 'Read "levstop" from namelist:', levstop
        write(iulog,*) 'Read "use_qdrai" from namelist:', use_qdrai
        write(iulog,*) 'Read "one_file_per_timestep" from namelist:', one_file_per_timestep

        if (levstart .gt. levstop) then
          call endrun(trim(subname)//'levstop must be bigger than levstart')
        endif

        if (max(levstart,  levstop) .gt. nlevsoi) then
          call endrun(trim(subname)//'levstop/ start must not exceed nlevsoi (10)')
        endif

        if (min(levstart,  levstop) .lt. 1) then
          call endrun(trim(subname)//'levstop/ start must not be smaller than 1')
        endif

        if (reservoir_capacity < 0._r8) then
          call endrun(trim(subname)//'reservoir_capacity must be >= 0')
        endif

        if (one_file_per_timestep) then
          if ((monthly) .or. (interp_day)) then
            call endrun(trim(subname)//'"monthly" and "interp_day" must be .false. if "one_file_per_timestep" is .true.')
          endif
        endif

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
      integer :: TimeStep_current
      integer :: doy         ! day of year (1..365)    
      integer :: offset      ! how much must we shift time?

      real(r8) :: nsec = 86400._r8 ! num of sec per day
      integer, dimension(12) :: ndaypm = &
           (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/) !days per month
      integer, dimension(12) :: cdaypm = &
          (/0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334/) ! cumulative
          ! number of days per month

  !-----------------------------------------------------------------------
      dtime = get_step_size()

      if (monthly) then 
        ! as originally programmed
        offset = int(dtime)
      else
        offset = - int(dtime / 2._r8)
      endif


      if ( is_perpetual() ) then
         call get_perp_date(kyr, kmo, kda, ksec, offset=offset)
      else
         call get_curr_date(kyr, kmo, kda, ksec, offset=offset)
      endif

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
          TimeStep_current = TimeStep(1)
      else ! interpolate daily data
        doy = cdaypm(kmo) + kda
        if (interp_day) then

          ! fraction of day that has passed
          t = ksec / nsec

          ! if t < 0.5 we need 'doy - 1' and 'doy'
          ! else we need 'doy' and 'doy + 1'
          TimeStep(1) = doy + floor(t - 0.5_r8)
          TimeStep(2) = doy + floor(t + 0.5_r8)

          timwt_soil(1) = TimeStep(2) - doy - t + 0.5_r8
          timwt_soil(2) = 1._r8 - timwt_soil(1)

          if (TimeStep(1) < 1) TimeStep(1) = 365
          if (TimeStep(2) > 365) TimeStep(2) = 1
          TimeStep_current = TimeStep(1)
        else ! no interpolation during one day
          
          if (one_file_per_timestep) then
            TimeStep(1) = 1 ! always one
            TimeStep(2) = 1 ! constant / not used
          else
            TimeStep(1) = doy
            TimeStep(2) = 1 ! constant / not used
          endif

          ! so that we know if we need to read a new file
          TimeStep_current = doy

          timwt_soil(1) = 1._r8 ! all weight
          timwt_soil(2) = 0._r8 ! no weight
        
        endif ! interp_day
      endif ! monthly

      
      if (TimeStep_old /= TimeStep_current) then
        call readSoilMoisture (kmo, kda, TimeStep)
        TimeStep_old = TimeStep_current
      endif

    end subroutine interpSoilMoisture





  !same as ReadMonthlyVegetation except for Soil Moisture
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: readSoilMoisture
  !
  ! !INTERFACE:
    subroutine readSoilMoisture (kmo, kda, TimeStep)
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
      use filenames,    only : interpret_filename_spec
      
  ! !ARGUMENTS:
      implicit none

      integer, intent(in) :: kmo         ! month (1, ..., 12)
      integer, intent(in) :: kda         ! day of month (1, ..., 31)
      integer, intent(in) :: TimeStep(2) ! months to be interpolated (1 to 12)
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
      integer :: k_end
      logical :: readvar

      real(r8), pointer :: mh2osoi_liq(:,:)  ! liquid soil water content read from input file
      real(r8), pointer :: mh2osoi_ice(:,:)  ! frozen soil water content read from input file
      character(len=32) :: subname = 'readSoilMoisture'
      character(len=32) :: cTimeStep ! string for the choosen time step (daily, monthly)

      character(len=256) :: pSMfile_local ! file name mit SM data to prescribe

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
      endif

      if (monthly) then
        cTimeStep = 'monthly'
      else
        cTimeStep = 'daily'
      endif

      if (one_file_per_timestep) then
        ! parse file name according to date/time of model
        ! %y -> 2005 %m -> 05 %d -> 20
        pSMfile_local = interpret_filename_spec( pSMfile )

        ! ensure the string has changed
        if (trim(pSMfile_local) == trim(pSMfile)) then
            call endrun(trim(subname)//'set "one_file_per_timestep" to .false. when not providing one file per timestep')
        endif

      else
        pSMfile_local = pSMfile
      endif



      ! ----------------------------------------------------------------------
      ! Open monthly soil moisture file
      ! Read data from column
      ! ----------------------------------------------------------------------
      
      if (masterproc) then
        write(iulog,*) 'Attempting to read ', trim(cTimeStep), ' soil moisture data...'
        write(iulog,*) 'month = ', kmo, ' day = ', kda
        write(iulog,*) 'from file: ', pSMfile_local
      endif ! masterproc
      
      if ((interp_day) .or. (monthly)) then
        k_end = 2
      else
        k_end = 1
      endif

      do k=1,k_end  ! loop over months/ days and read SM data

          ! get file
          call getfil(trim(pSMfile_local), locfn, 0)
          call ncd_pio_openfile (ncid, trim(locfn), 0)


          if (masterproc) then
              write(iulog,*) 'Before read, timestep= ', TimeStep(k) 
          endif ! masterproc


          call ncd_io(ncid=ncid, varname='SOILLIQ', flag='read', data=mh2osoi_liq, dim1name=namec, &
            nt=TimeStep(k), readvar=readvar)
          if (.not. readvar) call endrun(trim(subname) // ' ERROR: SOILLIQ NOT on pSMfile' // trim(pSMfile_local))
          

          call ncd_io(ncid=ncid, varname='SOILICE', flag='read', data=mh2osoi_ice, dim1name=namec, &
            nt=TimeStep(k), readvar=readvar)
          if (.not. readvar) call endrun(trim(subname) // ' ERROR: SOILICE NOT on pSMfile ' // trim(pSMfile_local))     

        call ncd_pio_closefile(ncid)


        if (masterproc) then
            write(iulog,*) 'Successfully read soil moisture data for timestep ', TimeStep(k)
        endif


        ! store data directly in clmtype structure
        do c = begc, endc
          l = clandunit(c)
          if (ltype(l) == istsoil) then 
            do j = 1, nlevsoi
              mh2osoi_liq2t(c, j, k) = mh2osoi_liq(c, j)
              mh2osoi_ice2t(c, j, k) = mh2osoi_ice(c, j)
            end do
          endif
        end do   ! end of loop over columns

      end do   ! end of loop over months
      
      if (masterproc) then
        write(iulog,*) 'Successfully read soil moisture data'
      endif ! masterproc

      deallocate(mh2osoi_liq, mh2osoi_ice)
    end subroutine readSoilMoisture




  end module prescribeSoilMoistureMod
