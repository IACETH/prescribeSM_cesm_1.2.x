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
  private :: readMonthlySoilMoisture     ! read SM data from file
!

! !PRIVATE TYPES:
  integer , private :: InterpMonths1            ! saved month index
  integer , private :: InterpMonths1_soil         ! saved month index
  real(r8), private :: timwt(2)                 ! time weights for month 1 and month 2
  real(r8), private :: timwt_soil(2)              ! time weights for month 1 and month 2
  real(r8), private, allocatable :: mlai2t(:,:) ! lai for interpolation (2 months)
  real(r8), private, allocatable :: msai2t(:,:) ! sai for interpolation (2 months)
  real(r8), private, allocatable :: mhvt2t(:,:) ! top vegetation height for interpolation (2 months)
  real(r8), private, allocatable :: mhvb2t(:,:) ! bottom vegetation height for interpolation(2 months)
  real(r8), private, allocatable :: mh2osoi_liq2t(:,:,:) !  liquid soil water for interpolation (2 months) read fro m input files
  real(r8), private, allocatable :: mh2osoi_ice2t(:,:,:) !  frozen soil water for interpolation (2 months) read fro m input files


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
  subroutine prescribeSoilMoisture(lbc,ubc,lbp,ubp, &
                                   num_nolakec,filter_nolakec)
!
! !DESCRIPTION:
! 
! 
! 
!
! !USES:
    use clmtype

    use clm_varpar, only : nlevsoi,nlevgrnd
    use clm_varcon, only : istsoil

!    use pftvarcon, only : noveg, ncorn, nbrdlf_dcd_brl_shrub
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: lbc, ubc                    ! column bounds


    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer , pointer :: ltype(:)                      ! landunit type index

!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: 
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)   ! column index associated with each pft
    integer , pointer :: clandunit(:) ! column's landunit
    real(r8), pointer :: snowdp(:)    ! snow height (m)
    integer , pointer :: ivt(:)       ! pft vegetation type

!
! local pointers to implicit out arguments
!

    real(r8), pointer :: h2osoi_liq(:,:)     ! liquid soil water content level
    real(r8), pointer :: h2osoi_ice(:,:)     ! frozen soil water content level


!
! !OTHER LOCAL VARIABLES:
!EOP
!
!    integer  :: fp,p,c   ! indices
    integer  :: fc,j,l,c      ! indices
 
!-----------------------------------------------------------------------


    !h2osoi_liq  => clm3%g%l%c%cws%h2osoi_liq
    !h2osoi_ice  => clm3%g%l%c%cws%h2osoi_ice

    h2osoi_liq  => cws%h2osoi_liq
    h2osoi_ice  => cws%h2osoi_ice

    !clandunit   => clm3%g%l%c%landunit
    clandunit   => col%landunit

    !ltype       => clm3%g%l%itype
    ltype       => lun%itype
    
    call interpMonthlySoilMoisture()

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       do j = 1, nlevsoi
          if (ltype(l) == istsoil) then
             h2osoi_liq(c,j) = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
             h2osoi_ice(c,j) = timwt_soil(1)*mh2osoi_ice2t(c,j,1) + timwt_soil(2)*mh2osoi_ice2t(c,j,2)
          endif
       end do
    end do
  end subroutine prescribeSoilMoisture









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
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clm_varctl, only : fsurdat
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
    integer :: months_soil(2) ! months to be interpolated (1 to 12)
    integer, dimension(12) :: ndaypm= &
         (/31,28,31,30,31,30,31,31,30,31,30,31/) !days per month
!-----------------------------------------------------------------------
    dtime = get_step_size()

    if ( is_perpetual() ) then
       call get_perp_date(kyr, kmo, kda, ksec, offset=int(dtime))
    else
       call get_curr_date(kyr, kmo, kda, ksec, offset=int(dtime))
    end if

    t = (kda-0.5_r8) / ndaypm(kmo)
    it(1) = t + 0.5_r8
    it(2) = it(1) + 1
    months_soil(1) = kmo + it(1) - 1
    months_soil(2) = kmo + it(2) - 1
    if (months_soil(1) <  1) months_soil(1) = 12
    if (months_soil(2) > 12) months_soil(2) = 1
    timwt_soil(1) = (it(1)+0.5_r8) - t
    timwt_soil(2) = 1._r8-timwt_soil(1)

    if (InterpMonths1_soil /= months_soil(1)) then
       call readMonthlySoilMoisture (kmo, kda, months_soil)
       InterpMonths1_soil = months_soil(1)
    end if

  end subroutine interpMonthlySoilMoisture





!same as ReadMonthlyVegetation except for Soil Moisture
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readMonthlySoilMoisture
!
! !INTERFACE:
  subroutine readMonthlySoilMoisture (kmo, kda, months_soil)
!
! !DESCRIPTION:
! Read monthly soil moisture data for two consec. months.
!
! !USES:
    use clmtype
    use shr_kind_mod, only : r8 => shr_kind_r8
    use decompMod   , only : get_proc_bounds, ldecomp, gsmap_lnd_gdc2glo
    use clm_varpar  , only : nlevgrnd
    use clm_varcon  , only : istsoil

    use fileutils   , only : getfil
    use spmdMod     , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
    use clm_time_manager, only : get_nstep
    use ncdio_pio
    use netcdf
    
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
    integer :: beg3d(3),len3d(3)          ! netCDF variable edges
    integer :: ntim                       ! number of input data time samples
    integer :: ncolumn_i                  ! number of columns
    integer :: nlev_i                     ! number of input data soil levels
    integer :: begc,endc                  ! beg and end local c index
    integer :: begp,endp
    integer :: ier,ret                    ! error code
    integer :: closelatidx,closelonidx
    real(r8):: closelat,closelon
    logical :: readvar

    real(r8), pointer :: mh2osoi_liq(:,:)  ! liquid soil water content read from input file
    real(r8), pointer :: mh2osoi_ice(:,:)  ! frozen soil water content read from input file
    real(r8), pointer :: arrayl(:)  ! temp local array
    character(len=32) :: subname = 'readMonthlySoilMoisture'
!-----------------------------------------------------------------------

    

  
    


    ! Assign local pointers to derived subtypes components (landunit-level)
    !ltype               => clm3%g%l%itype
    ltype       => lun%itype
    ! Assign local pointers to derived subtypes components (column-level)
    !clandunit       => clm3%g%l%c%landunit
    clandunit   => col%landunit
    ! Determine necessary indices


    call get_proc_bounds(begc=begc,endc=endc,begp=begp,endp=endp)
    write(iulog,*) 'begc, endc ',begc,endc
    write(iulog,*) 'begp, endp ',begp,endp

    allocate(mh2osoi_liq(begc:endc,1:nlevgrnd), &
             mh2osoi_ice(begc:endc,1:nlevgrnd), stat=ier)

    if (ier /= 0) then
       write(iulog,*)subname, 'allocation big error '; call endrun()
    end if

    ! ----------------------------------------------------------------------
    ! Open monthly soil moisture file
    ! Read data from column
    ! ----------------------------------------------------------------------
    
    if (masterproc) then

          write(iulog,*) 'Attempting to read monthly soil moisture data .....'
          write(iulog,*) 'nstep = ',get_nstep(),' month = ',kmo,' day = ',kda

          !call check_ret(nf90_open(locfn, 0, ncid), subname)
          !call check_ret(nf90_inq_dimid (ncid, 'column', dimid), subname)
          !call check_ret(nf90_inquire_dimension(ncid, dimid, len=ncolumn_i), subname)
          !call check_ret(nf90_inq_dimid(ncid, 'levgrnd', dimid), subname)
          !call check_ret(nf90_inquire_dimension(ncid, dimid, len=nlev_i), subname)
          !call check_ret(nf90_inq_dimid(ncid, 'time', dimid), subname)
          !call check_ret(nf90_inquire_dimension(ncid, dimid, len=ntim), subname)

       endif   ! masterproc

    call getfil('/cluster/home03/uwis/mathause/data/SM_test.nc', locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    write(iulog,*) 'size ',size(mh2osoi_liq)
    do k=1,2   !loop over months and read vegetated data

       !allocate(arrayl(begc:endc),stat=ier)
       !if (ier /= 0) then
       !   write(iulog,*)subname, 'allocation array error '; call endrun()
       !end if
       #if (masterproc) then
         write(iulog,*) 'Before read '
         write(iulog,*) size(mh2osoi_liq)
       #endif ! masterproc
        !do j = 1,nlevgrnd
          !beg3d(1) = j         ; len3d(1) = 1
          !beg3d(2) = 1         ; len3d(2) = ncolumn_i
          !beg3d(3) = months_soil(k) ; len3d(3) = 1

        call ncd_io(ncid=ncid, varname='SOILLIQ', flag='read', data=mh2osoi_liq, dim1name=grlnd, &
          nt=months_soil(k), readvar=readvar)
        if (.not. readvar) call endrun( trim(subname)//' ERROR: SOILLIQ NOT on pSM file' )

        call ncd_io(ncid=ncid, varname='SOILICE', flag='read', data=mh2osoi_liq, dim1name=grlnd, &
          nt=months_soil(k), readvar=readvar)
        if (.not. readvar) call endrun( trim(subname)//' ERROR: SOILICE NOT on pSM file' )

          !mh2osoi_liq(begc:endc,j) = arrayl(begc:endc)

          !call ncd_iolocal(ncid,'SOILLIQ','read',arrayl,namec,beg3d,len3d,status=ret)
          !write(iulog,*) 'After read '
          !if (ret /= 0) call endrun( trim(subname)//' ERROR: SOILLIQ NOT on clim file' )
          

          !call ncd_iolocal(ncid,'SOILICE','read',arrayl,namec,beg3d,len3d,status=ret)
          !if (ret /= 0) call endrun( trim(subname)//' ERROR: SOILICE NOT on clim file' )
          !mh2osoi_ice(begc:endc,j) = arrayl(begc:endc)

       !enddo

       !deallocate(arrayl)
      call ncd_pio_closefile(ncid)
      if (masterproc) then
          !call check_ret(nf90_close(ncid), subname)
          write(iulog,*) 'Successfully read monthly soil moisture data for'
          write(iulog,*) 'month ', months_soil(k)
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

  end subroutine readMonthlySoilMoisture




end module prescribeSoilMoistureMod
