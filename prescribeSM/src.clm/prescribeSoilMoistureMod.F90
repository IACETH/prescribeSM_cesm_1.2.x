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

  use shr_kind_mod,    only : r8 => shr_kind_r8
  use clm_varpar,      only : nlevgrnd
  use abortutils,      only : endrun
  use clm_varctl,      only : iulog
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: prescribeSoilMoisture  ! Sets Soil Moisture

!

! !PRIVATE MEMBER FUNCTIONS:
  private :: interpMonthlySoilMoisture   ! interpolate monthly SM data
  private :: readSoilMoisture     ! read SM data from file
!

! !PRIVATE TYPES:

  integer , private :: TimeStep_old = 0           ! time step at last call
  real(r8), private :: timwt_soil(2)              ! time weights for month 1 and month 2
  real(r8), private, allocatable :: mh2osoi_liq2t(:,:,:) !  liquid soil water for interpolation (2 months) read fro m input files
  real(r8), private, allocatable :: mh2osoi_ice2t(:,:,:) !  frozen soil water for interpolation (2 months) read fro m input files

  logical, private            :: monthly  ! if .true. -> monthly, else daily
  character(len=256), private :: pSMfile  ! file name mit SM data to prescribe

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
! 
! 
! 
!
! !USES:
    use clmtype
    use decompMod,  only : get_proc_bounds
    use clm_varpar, only : nlevgrnd
    use clm_varcon, only : istsoil
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
!    use pftvarcon, only : noveg, ncorn, nbrdlf_dcd_brl_shrub
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column bounds


    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, pointer    :: ltype(:)                    ! landunit type index

!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: 
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: clandunit(:) ! column's landunit

!
! local pointers to implicit out arguments
!

    real(r8), pointer :: h2osoi_liq(:,:)     ! liquid soil water content level
    real(r8), pointer :: h2osoi_ice(:,:)     ! frozen soil water content level


!
! !OTHER LOCAL VARIABLES:
!EOP
!

    integer :: fc,j,l,c         ! indices
    integer :: begc,endc        ! beg and end local c index
    integer :: ier              ! error code
!-----------------------------------------------------------------------


    ! implicit inout arguments
    h2osoi_liq  => cws%h2osoi_liq
    h2osoi_ice  => cws%h2osoi_ice

    ! implicit in arguments
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

    ! get file name and monthly or daily
    call initPrescribeSoilMoisture ()
    ! get time weight and possibly a new time step
    call interpMonthlySoilMoisture()


    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       if (ltype(l) == istsoil) then
        do j = 1, nlevgrnd
             h2osoi_liq(c,j) = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
             h2osoi_ice(c,j) = timwt_soil(1)*mh2osoi_ice2t(c,j,1) + timwt_soil(2)*mh2osoi_ice2t(c,j,2)
        end do
       end if
    end do

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
! goal -> use namlist
!
! !USES:
  use fileutils        , only : getfil, getavu, relavu
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
    character(len=256) :: locfn        ! local file name
    character(len=256) :: subname      ! name of routine
!EOP

	subname = 'initPrescribeSoilMoisture'

    ! Input datasets

    namelist /prescribe_sm/  &
         pSMfile, monthly




    ! CAREFUL: IF YOU USE A DAILY pSMfile BUT SET monthly = .true.
    ! THE MODEL WILL STILL RUN! 

    ! pSMfile = '/cluster/home03/uwis/mathause/data/SM_in_3D_daily.nc'
    ! monthly = .false.


    call getfil(trim('prescribe_SM_nl'), locfn, 0)
    
    write(iulog,*) 'local file name: ', trim(locfn)
    

    unitn = getavu()
    
    write(iulog,*) 'Read in prescribe_sm namelist from: prescribe_SM_nl'
    
    open( unitn, file=trim(locfn), status='old' )
    ierr = 1
    do while ( ierr /= 0 )
        read(unitn, prescribe_sm, iostat=ierr)
        if (ierr < 0) then
            call endrun( trim(subname)//' encountered end-of-file on clm_inparm read' )
        endif
    end do
    call relavu( unitn )


    write(iulog,*) 'Read pSMfile from namelist', trim(pSMfile)
    write(iulog,*) 'Read monthly from namelist', monthly




  end subroutine initPrescribeSoilMoisture



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: interpMonthlySoil
!
! !INTERFACE:
  subroutine interpMonthlySoilMoisture ()
!
! !DESCRIPTION:
! Determine if 2 new months of data are to be read.
!
! !USES:
    use shr_kind_mod,     only : r8 => shr_kind_r8
    use clm_varctl,       only : fsurdat
    use clm_time_manager, only : get_curr_date, get_step_size, get_perp_date, is_perpetual
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Ruth Lorenz
!
!EOP
!
! LOCAL VARIABLES:
    integer :: kyr         ! year (0, ...) for nstep+1
    integer :: kmo         ! month (1, ..., 12)
    integer :: kda         ! day of month (1, ..., 31)
    integer :: ksec        ! seconds into current date for nstep+1
    real(r8):: dtime       ! land model time step (sec)
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

        ! if t < 0.5 we need 'doy -1' and 'doy'
        ! else we need 'doy' and 'doy + 1'
        TimeStep(1) = doy + floor(t - 0.5)
        TimeStep(2) = doy + floor(t + 0.5)

        if (TimeStep(1) < 1) TimeStep(1) = 365
        if (TimeStep(2) > 365) TimeStep(2) = 1

        timwt_soil(1) = 1._r8 - abs(t - 0.5_r8)
        timwt_soil(2) = 1._r8 - timwt_soil(1)

    endif ! monthly

    
    if (TimeStep_old /= TimeStep(1)) then
        write(iulog,*) 'TimeStep_old ',  TimeStep_old
        write(iulog,*) 'timwt_soil(1) ',  timwt_soil(1)
       call readSoilMoisture (kmo, kda, TimeStep)
       TimeStep_old = TimeStep(1)
    end if

  end subroutine interpMonthlySoilMoisture





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
    use spmdMod,      only : masterproc
    
! !ARGUMENTS:
    implicit none

    integer, intent(in) :: kmo            ! month (1, ..., 12)
    integer, intent(in) :: kda            ! day of month (1, ..., 31)
    integer, intent(in) :: months_soil(2) ! months to be interpolated (1 to 12)
!
! !REVISION HISTORY:
! Created by Ruth Lorenz
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256) :: locfn           ! local file name
    integer :: g,n,i,j,k,l,m,c            ! indices
    integer , pointer :: clandunit(:)     ! column's landunit
    integer , pointer :: ltype(:)         ! landunit type index
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
    real(r8), pointer :: arrayl(:)  ! temp local array
    character(len=32) :: subname = 'readSoilMoisture'
    character(len=32) :: cTimeStep
!-----------------------------------------------------------------------



    ! Assign local pointers to derived subtypes components (landunit-level)
    ltype       => lun%itype
    ! Assign local pointers to derived subtypes components (column-level)
    !clandunit       => clm3%g%l%c%landunit
    clandunit   => col%landunit
    ! Determine necessary indices


    call get_proc_bounds(begc=begc,endc=endc)



    allocate(mh2osoi_liq(begc:endc,1:nlevgrnd), &
             mh2osoi_ice(begc:endc,1:nlevgrnd), stat=ier)

    if (ier /= 0) then
       write(iulog,*)subname, 'allocation big error '; call endrun()
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
          !call check_ret(nf90_close(ncid), subname)
          write(iulog,*) 'Successfully read ',  trim(cTimeStep), 'ly soil moisture data for'
          write(iulog,*) cTimeStep, ' ', months_soil(k)
          write(iulog,*)
      end if


       ! store data directly in clmtype structure
       do c = begc,endc
          l = clandunit(c)
          if (ltype(l) == istsoil) then 
             do j = 1, nlevgrnd
                mh2osoi_liq2t(c,j,k) = mh2osoi_liq(c,j)
                mh2osoi_ice2t(c,j,k) = mh2osoi_ice(c,j)
             end do
          end if
       end do   ! end of loop over columns

    end do   ! end of loop over months

    write(iulog,*) 'Successfully read soil moisture data'

    deallocate(mh2osoi_liq, mh2osoi_ice)

  end subroutine readSoilMoisture




end module prescribeSoilMoistureMod
