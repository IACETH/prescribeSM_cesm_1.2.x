#include <misc.h>
#include <preproc.h>
module STATICEcosysdynMOD

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: STATICEcosysDynMod
!
! !DESCRIPTION:
! Static Ecosystem dynamics: phenology, vegetation. This is for the CLM Satelitte Phenology 
! model (CLMSP). Allow some subroutines to be used by the CLM Carbon Nitrogen model (CLMCN) 
! so that DryDeposition code can get estimates of LAI differences between months.
!
! !USES:
  use shr_kind_mod,    only : r8 => shr_kind_r8
  use abortutils,      only : endrun
  use clm_varctl,      only : scmlat,scmlon,single_column
  use clm_varctl,      only : iulog
  use spmdGathScatMod, only : scatter_data_from_master
!RL
  use clm_varpar,      only : nlevgrnd
!RL
  use perf_mod,        only : t_startf, t_stopf
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
#if (!defined CN) && (!defined CNDV)
  public :: EcosystemDyn         ! CLMSP Ecosystem dynamics: phenology, vegetation
#endif
  public :: EcosystemDynini      ! Dynamically allocate memory
  public :: interpMonthlyVeg     ! interpolate monthly vegetation data
!RL
  public :: interpMonthlySoil   ! interpolate monthly soil data
!RL
  public :: readAnnualVegetation ! Read in annual vegetation (needed for Dry-deposition)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: readMonthlyVegetation   ! read monthly vegetation data for two months
!RL
!  public :: readMonthlySoil   ! read monthly soil data
!DML
  private :: readMonthlySoil   ! read monthly soil data
!DML
!RL
!
! !PRIVATE TYPES:
  integer , private :: InterpMonths1            ! saved month index
!RL
  integer , private :: InterpMonths1_soil         ! saved month index
!RL
  real(r8), private :: timwt(2)                 ! time weights for month 1 and month 2
!RL
  real(r8), private :: timwt_soil(2)              ! time weights for month 1 and month 2
!RL
  real(r8), private, allocatable :: mlai2t(:,:) ! lai for interpolation (2 months)
  real(r8), private, allocatable :: msai2t(:,:) ! sai for interpolation (2 months)
  real(r8), private, allocatable :: mhvt2t(:,:) ! top vegetation height for interpolation (2 months)
  real(r8), private, allocatable :: mhvb2t(:,:) ! bottom vegetation height for interpolation(2 months)
!RL
  real(r8), private, allocatable :: mh2osoi_liq2t(:,:,:) !  liquid soil water for interpolation (2 months) read fro m input files
  real(r8), private, allocatable :: mh2osoi_ice2t(:,:,:) !  frozen soil water for interpolation (2 months) read fro m input files
!RL
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EcosystemDynini
!
! !INTERFACE:
  subroutine EcosystemDynini ()
!
! !DESCRIPTION:
! Dynamically allocate memory and set to signaling NaN.
!
! !USES:
    use nanMod
    use decompMod, only : get_proc_bounds
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ier    ! error code
    integer :: begp,endp  ! local beg and end p index
!RL
    integer :: begc,endc  ! local beg and end column index
!RL
!-----------------------------------------------------------------------

    InterpMonths1 = -999  ! saved month index
!RL
!    call get_proc_bounds(begp=begp,endp=endp)
    call get_proc_bounds(begp=begp,endp=endp,begc=begc,endc=endc)
!RL

    ier = 0
    if(.not.allocated(mlai2t))allocate (mlai2t(begp:endp,2), &
              msai2t(begp:endp,2), &
              mhvt2t(begp:endp,2), &
              mhvb2t(begp:endp,2), stat=ier)
!RL
    if(.not.allocated(mh2osoi_liq2t))allocate (mh2osoi_liq2t(begc:endc,1:nlevgrnd,2), &
              mh2osoi_ice2t(begc:endc,1:nlevgrnd,2), stat=ier)
!RL
    if (ier /= 0) then
       write(iulog,*) 'EcosystemDynini allocation error'
       call endrun
    end if

    mlai2t(:,:) = nan
    msai2t(:,:) = nan
    mhvt2t(:,:) = nan
    mhvb2t(:,:) = nan
!RL
    mh2osoi_liq2t(:,:,:) = nan
    mh2osoi_ice2t(:,:,:) = nan
!RL

  end subroutine EcosystemDynini

#if (!defined CN) && (!defined CNDV)
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EcosystemDyn
!
! !INTERFACE:
!DML
!  subroutine EcosystemDyn(lbp, ubp, num_nolakep, filter_nolakep, doalb)
  subroutine EcosystemDyn(lbc,ubc,num_nolakec,filter_nolakec,lbp,ubp,num_nolakep,filter_nolakep,doalb)
!DML
!
! !DESCRIPTION:
! Ecosystem dynamics: phenology, vegetation
! Calculates leaf areas (tlai, elai),  stem areas (tsai, esai) and
! height (htop).
!
! !USES:
    use clmtype
!RL
    use clm_varpar, only : nlevsoi,nlevgrnd
!RL
!DML
    use clm_varcon  , only : istsoil
!DML
    use pftvarcon, only : noveg, ncorn, nbrdlf_dcd_brl_shrub
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
!RL
    integer, intent(in) :: lbc, ubc                    ! column bounds
!RL
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)   ! pft filter for non-lake points
!RL
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
!RL
!DML
    integer , pointer :: ltype(:)         ! landunit type index
!DML
    logical, intent(in) :: doalb                       ! true = surface albedo calculation time step
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/1/02, Peter Thornton: Migrated to new data structure.
! 2/29/08, David Lawrence: revised snow burial fraction for short vegetation   
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)  ! column index associated with each pft
    integer , pointer :: clandunit(:)       ! column's landunit
    real(r8), pointer :: snowdp(:)   ! snow height (m)
    integer , pointer :: ivt(:)      ! pft vegetation type
#if (defined CASA)
    real(r8), pointer :: plai(:)     ! Prognostic lai
#endif
!
! local pointers to implicit out arguments
!
!RL
    real(r8), pointer :: h2osoi_liq(:,:)     ! liquid soil water content level
    real(r8), pointer :: h2osoi_ice(:,:)     ! frozen soil water content level
!RL
    real(r8), pointer :: tlai(:)     ! one-sided leaf area index, no burying by snow
    real(r8), pointer :: tsai(:)     ! one-sided stem area index, no burying by snow
    real(r8), pointer :: htop(:)     ! canopy top (m)
    real(r8), pointer :: hbot(:)     ! canopy bottom (m)
    real(r8), pointer :: elai(:)     ! one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)     ! one-sided stem area index with burying by snow
    integer , pointer :: frac_veg_nosno_alb(:) ! frac of vegetation not covered by snow [-]
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: fp,p,c   ! indices
!RL
    integer  :: fc,j,l      ! indices
!RL
    real(r8) :: ol       ! thickness of canopy layer covered by snow (m)
    real(r8) :: fb       ! fraction of canopy layer covered by snow
!-----------------------------------------------------------------------

!RL
       h2osoi_liq  => clm3%g%l%c%cws%h2osoi_liq
       h2osoi_ice  => clm3%g%l%c%cws%h2osoi_ice

       clandunit   => clm3%g%l%c%landunit
       ltype       => clm3%g%l%itype

!DML
!       do f = 1, num_nolakec
!          c = filter_nolakec(f)
!           do j = 1, nlevgrnd
!              h2osoi_liq(c,j) = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
!              h2osoi_ice(c,j) = timwt_soil(1)*mh2osoi_ice2t(c,j,1) + timwt_soil(2)*mh2osoi_ice2t(c,j,2)
!           end do
!       end do
!DML
!RL
    if (doalb) then

       ! Assign local pointers to derived type scalar members (column-level)

       snowdp  => clm3%g%l%c%cps%snowdp

       ! Assign local pointers to derived type scalar members (pftlevel)

       pcolumn => clm3%g%l%c%p%column
       tlai    => clm3%g%l%c%p%pps%tlai
       tsai    => clm3%g%l%c%p%pps%tsai
       elai    => clm3%g%l%c%p%pps%elai
       esai    => clm3%g%l%c%p%pps%esai
       htop    => clm3%g%l%c%p%pps%htop
       hbot    => clm3%g%l%c%p%pps%hbot
       frac_veg_nosno_alb => clm3%g%l%c%p%pps%frac_veg_nosno_alb
       ivt     => clm3%g%l%c%p%itype
#if (defined CASA)
       plai    => clm3%g%l%c%p%pps%plai
#endif

!dir$ concurrent
!cdir nodep
       do fp = 1, num_nolakep
          p = filter_nolakep(fp)
          c = pcolumn(p)

          ! need to update elai and esai only every albedo time step so do not
          ! have any inconsistency in lai and sai between SurfaceAlbedo calls (i.e.,
          ! if albedos are not done every time step).
          ! leaf phenology
          ! Set leaf and stem areas based on day of year
          ! Interpolate leaf area index, stem area index, and vegetation heights
          ! between two monthly
          ! The weights below (timwt(1) and timwt(2)) were obtained by a call to
          ! routine InterpMonthlyVeg in subroutine NCARlsm.
          !                 Field   Monthly Values
          !                -------------------------
          ! leaf area index LAI  <- mlai1 and mlai2
          ! leaf area index SAI  <- msai1 and msai2
          ! top height      HTOP <- mhvt1 and mhvt2
          ! bottom height   HBOT <- mhvb1 and mhvb2

          tlai(p) = timwt(1)*mlai2t(p,1) + timwt(2)*mlai2t(p,2)
          tsai(p) = timwt(1)*msai2t(p,1) + timwt(2)*msai2t(p,2)
          htop(p) = timwt(1)*mhvt2t(p,1) + timwt(2)*mhvt2t(p,2)
          hbot(p) = timwt(1)*mhvb2t(p,1) + timwt(2)*mhvb2t(p,2)

#if (defined CASA)
! use PLAI from CASA
!! can also choose to use default tlai instead of plai - need option
! 03/03/11: don't reset lai of crops - use default lsm tlai/elai

          if (ivt(p) > noveg .and. ivt(p) < ncorn ) tlai(p) = plai(p)
#endif

          ! adjust lai and sai for burying by snow. if exposed lai and sai
          ! are less than 0.05, set equal to zero to prevent numerical
          ! problems associated with very small lai and sai.

          ! snow burial fraction for short vegetation (e.g. grasses) as in
          ! Wang and Zeng, 2007. 

          if (ivt(p) > noveg .and. ivt(p) <= nbrdlf_dcd_brl_shrub ) then
             ol = min( max(snowdp(c)-hbot(p), 0._r8), htop(p)-hbot(p))
             fb = 1._r8 - ol / max(1.e-06_r8, htop(p)-hbot(p))
          else
             fb = 1._r8 - max(min(snowdp(c),0.2_r8),0._r8)/0.2_r8   ! 0.2m is assumed
                  !depth of snow required for complete burial of grasses
          endif

          elai(p) = max(tlai(p)*fb, 0.0_r8)
          esai(p) = max(tsai(p)*fb, 0.0_r8)
          if (elai(p) < 0.05_r8) elai(p) = 0._r8
          if (esai(p) < 0.05_r8) esai(p) = 0._r8

          ! Fraction of vegetation free of snow

          if ((elai(p) + esai(p)) >= 0.05_r8) then
             frac_veg_nosno_alb(p) = 1
          else
             frac_veg_nosno_alb(p) = 0
          end if

       end do ! end of pft loop

    end if  !end of if-doalb block

!DML
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       do j = 1, nlevsoi
          if (ltype(l) == istsoil) then
!             write(iulog,*) 'LIQ1 ',j,c,h2osoi_liq(c,j)
             h2osoi_liq(c,j) = timwt_soil(1)*mh2osoi_liq2t(c,j,1) + timwt_soil(2)*mh2osoi_liq2t(c,j,2)
             h2osoi_ice(c,j) = timwt_soil(1)*mh2osoi_ice2t(c,j,1) + timwt_soil(2)*mh2osoi_ice2t(c,j,2)
!             write(iulog,*) 'LIQ2 ',j,c,h2osoi_liq(c,j),mh2osoi_liq2t(c,j,1)
!             write(iulog,*) 'ICE2 ',j,c,h2osoi_ice(c,j),mh2osoi_ice2t(c,j,1)
          endif
       end do
    end do
!DML

  end subroutine EcosystemDyn

#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: interpMonthlyVeg
!
! !INTERFACE:
  subroutine interpMonthlyVeg ()
!
! !DESCRIPTION:
! Determine if 2 new months of data are to be read.
!
! !USES:
    use clm_varctl  , only : fsurdat
    use clm_time_manager, only : get_curr_date, get_step_size, get_perp_date, is_perpetual
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: kyr         ! year (0, ...) for nstep+1
    integer :: kmo         ! month (1, ..., 12)
    integer :: kda         ! day of month (1, ..., 31)
    integer :: ksec        ! seconds into current date for nstep+1
    real(r8):: dtime       ! land model time step (sec)
    real(r8):: t           ! a fraction: kda/ndaypm
    integer :: it(2)       ! month 1 and month 2 (step 1)
    integer :: months(2)   ! months to be interpolated (1 to 12)
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
    months(1) = kmo + it(1) - 1
    months(2) = kmo + it(2) - 1
    if (months(1) <  1) months(1) = 12
    if (months(2) > 12) months(2) = 1
    timwt(1) = (it(1)+0.5_r8) - t
    timwt(2) = 1._r8-timwt(1)

    if (InterpMonths1 /= months(1)) then
       call t_startf('readMonthlyVeg')
       call readMonthlyVegetation (fsurdat, kmo, kda, months)
       InterpMonths1 = months(1)
       call t_stopf('readMonthlyVeg')
    end if

  end subroutine interpMonthlyVeg

!RL
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: interpMonthlySoil
!
! !INTERFACE:
  subroutine interpMonthlySoil ()
!
! !DESCRIPTION:
! Determine if 2 new months of data are to be read.
!
! !USES:
    use clm_varctl  , only : fsurdat
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
       call readMonthlySoil (kmo, kda, months_soil)
       InterpMonths1_soil = months_soil(1)
    end if

  end subroutine interpMonthlySoil
!RL

!-----------------------------------------------------------------------
! read 12 months of veg data for dry deposition
!-----------------------------------------------------------------------
  subroutine readAnnualVegetation ( )

    use clmtype
    use clm_varpar  , only : lsmlon, lsmlat, numpft
    use pftvarcon   , only : noveg
    use decompMod   , only : get_proc_bounds, ldecomp
    use spmdMod     , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
    use ncdio       , only : check_ret,ncd_iolocal
    use fileutils   , only : getfil
    use clm_varctl  , only : fsurdat
    use netcdf

    implicit none

    ! local vars

    real(r8), pointer     :: annlai(:,:)  ! 12 months of monthly lai from input data set 
    real(r8), allocatable :: mlai(:,:)    ! lai read from input files
    real(r8), pointer :: arrayl(:)        ! temp local array
    character(len=32) :: subname = 'readAnnualVegetation'
    integer :: ier                        ! error code
    character(len=256) :: locfn           ! local file name
    integer :: g,k,l,m,n,p,ivt            ! indices
    integer :: ncid,dimid,varid           ! input netCDF id's
    integer :: beg4d(4),len4d(4)          ! netCDF variable edges
    integer :: ntim                       ! number of input data time samples
    integer :: nlon_i                     ! number of input data longitudes
    integer :: nlat_i                     ! number of input data latitudes
    integer :: npft_i                     ! number of input data pft types
    integer :: begp,endp                  ! beg and end local p index
    integer :: nps,npe              ! start/end numpft limits
    integer :: closelatidx,closelonidx
    real(r8):: closelat,closelon
    integer :: begg,endg                  ! beg and end local g index

    annlai    => clm3%g%l%c%p%pps%annlai

    ! Determine necessary indices

    call get_proc_bounds(begg=begg,endg=endg,begp=begp,endp=endp)

    nps = 0
    npe = numpft

    allocate(mlai(begg:endg,0:numpft), stat=ier)
    if (ier /= 0) then
       write(6,*)subname, 'allocation error '; call endrun()
    end if

    call get_proc_bounds(begp=begp,endp=endp)

    if (masterproc) then

       write (6,*) 'Attempting to read annual vegetation data .....'

       call getfil(fsurdat, locfn, 0)
       call check_ret(nf90_open(locfn, 0, ncid), subname)

       call check_ret(nf90_inq_dimid (ncid, 'lsmlon', dimid), subname)
       call check_ret(nf90_inquire_dimension(ncid, dimid, len=nlon_i), subname)
       if (.not.single_column) then
          if (nlon_i /= lsmlon) then
             write(6,*)subname,' parameter lsmlon= ',lsmlon,'does not equal input nlat_i= ',nlon_i
             call endrun()
          end if
       endif
       call check_ret(nf90_inq_dimid(ncid, 'lsmlat', dimid), subname)
       call check_ret(nf90_inquire_dimension(ncid, dimid, len=nlat_i), subname)

       if (.not.single_column ) then
          if (nlat_i /= lsmlat) then
             write(6,*)subname,' parameter lsmlat= ',lsmlat,'does not equal input nlat_i= ',nlat_i
             call endrun()
          end if
       endif
       call check_ret(nf90_inq_dimid(ncid, 'lsmpft', dimid), subname)
       call check_ret(nf90_inquire_dimension(ncid, dimid, len=npft_i), subname)

       if (.not. single_column ) then
          if (npft_i /= numpft+1) then
             write(6,*)subname,' parameter numpft+1 = ',numpft+1,'does not equal input npft_i= ',npft_i
             call endrun()
          end if
       endif
       call check_ret(nf90_inq_dimid(ncid, 'time', dimid), subname)
       call check_ret(nf90_inquire_dimension(ncid, dimid, len=ntim), subname)

       call check_ret(nf90_inq_varid(ncid, 'MONTHLY_LAI', varid), subname)

    endif ! masterproc

    call mpi_bcast (nlon_i, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nlat_i, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (npft_i, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (ntim  , 1, MPI_INTEGER, 0, mpicom, ier)

    if (single_column) then
       call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
    endif

    allocate(arrayl(begg:endg),stat=ier)
    if (ier /= 0) then
       write(6,*)subname, 'allocation array error '; call endrun()
    end if

    do k=1,12   !! loop over months and read vegetated data

       do n = nps,npe
          beg4d(1) = 1         ; len4d(1) = nlon_i
          beg4d(2) = 1         ; len4d(2) = nlat_i
          beg4d(3) = n         ; len4d(3) = 1
          beg4d(4) = k         ; len4d(4) = 1
          if (single_column) then
             beg4d(1) = closelonidx; len4d(1) = 1
             beg4d(2) = closelatidx; len4d(2) = 1
          else
             beg4d(1) = 1         ; len4d(1) = nlon_i
             beg4d(2) = 1         ; len4d(2) = nlat_i
          end if
          if (nps == 0) beg4d(3) = n+1     ! shift by "1" if index starts at zero

          call ncd_iolocal(ncid,'MONTHLY_LAI','read',arrayl,grlnd,beg4d,len4d)
          mlai(begg:endg,n) = arrayl(begg:endg)

          !! store data directly in clmtype structure
          !! only vegetated pfts have nonzero values

       enddo

       do p = begp,endp

          g = clm3%g%l%c%p%gridcell(p)

          !! Assign lai/sai/hgtt/hgtb to the top [maxpatch_pft] pfts
          !! as determined in subroutine surfrd

          ivt = clm3%g%l%c%p%itype(p)
          if (ivt /= noveg) then     !! vegetated pft
             do l = 0, numpft
                if (l == ivt) then

                   annlai(k,p) = mlai(g,l)

                end if
             end do
          else                       !! non-vegetated pft

             annlai(k,p) = 0._r8

          end if

       end do   ! end of loop over pfts  

    enddo ! months loop

    if (masterproc) then
       call check_ret(nf90_close(ncid), subname)
    endif

    deallocate(mlai)

  endsubroutine readAnnualVegetation

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readMonthlyVegetation
!
! !INTERFACE:
  subroutine readMonthlyVegetation (fveg, kmo, kda, months)
!
! !DESCRIPTION:
! Read monthly vegetation data for two consec. months.
!
! !USES:
    use clmtype
    use decompMod   , only : get_proc_bounds, ldecomp, gsmap_lnd_gdc2glo
    use clm_varpar  , only : lsmlon, lsmlat, numpft
    use pftvarcon   , only : noveg
    use fileutils   , only : getfil
    use spmdMod     , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
    use clm_time_manager, only : get_nstep
    use ncdio       , only : check_ret,ncd_iolocal
    use netcdf
!
! !ARGUMENTS:
    implicit none

    character(len=*), intent(in) :: fveg  ! file with monthly vegetation data
    integer, intent(in) :: kmo            ! month (1, ..., 12)
    integer, intent(in) :: kda            ! day of month (1, ..., 31)
    integer, intent(in) :: months(2)      ! months to be interpolated (1 to 12)
!
! !REVISION HISTORY:
! Created by Sam Levis
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256) :: locfn           ! local file name
    integer :: g,n,k,l,m,p,ivt            ! indices
    integer :: ncid,dimid,varid           ! input netCDF id's
    integer :: beg4d(4),len4d(4)          ! netCDF variable edges
    integer :: ntim                       ! number of input data time samples
    integer :: nlon_i                     ! number of input data longitudes
    integer :: nlat_i                     ! number of input data latitudes
    integer :: npft_i                     ! number of input data pft types
    integer :: begp,endp                  ! beg and end local p index
    integer :: begg,endg                  ! beg and end local g index
    integer :: ier,ret                    ! error code
    integer :: closelatidx,closelonidx
    real(r8):: closelat,closelon

    integer :: nps,npe              ! start/end numpft limits
    real(r8), pointer :: mlai(:,:)  ! lai read from input files
    real(r8), pointer :: msai(:,:)  ! sai read from input files
    real(r8), pointer :: mhgtt(:,:) ! top vegetation height
    real(r8), pointer :: mhgtb(:,:) ! bottom vegetation height
    real(r8), pointer :: arrayl(:)  ! temp local array
    real(r8), pointer :: mlaidiff(:)! difference between lai month one and month two
    character(len=32) :: subname = 'readMonthlyVegetation'
!-----------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg=begg,endg=endg,begp=begp,endp=endp)

    nps = 0
    npe = numpft
    allocate(mlai(begg:endg,0:numpft), &
             msai(begg:endg,0:numpft), &
             mhgtt(begg:endg,0:numpft), &
             mhgtb(begg:endg,0:numpft), stat=ier)

    if (ier /= 0) then
       write(iulog,*)subname, 'allocation big error '; call endrun()
    end if

    ! ----------------------------------------------------------------------
    ! Open monthly vegetation file
    ! Read data and convert from [lsmlon] x [lsmlat] grid to patch data
    ! ----------------------------------------------------------------------

    do k=1,2   !loop over months and read vegetated data

       if (masterproc) then

          write(iulog,*) 'Attempting to read monthly vegetation data .....'
          write(iulog,*) 'nstep = ',get_nstep(),' month = ',kmo,' day = ',kda

          call getfil(fveg, locfn, 0)
          call check_ret(nf90_open(locfn, 0, ncid), subname)

          call check_ret(nf90_inq_dimid (ncid, 'lsmlon', dimid), subname)
          call check_ret(nf90_inquire_dimension(ncid, dimid, len=nlon_i), subname)
          if (.not.single_column) then
             if (nlon_i /= lsmlon) then
                write(iulog,*)subname,' parameter lsmlon= ',lsmlon,'does not equal input nlat_i= ',nlon_i
                call endrun()
             end if
          endif
          call check_ret(nf90_inq_dimid(ncid, 'lsmlat', dimid), subname)
          call check_ret(nf90_inquire_dimension(ncid, dimid, len=nlat_i), subname)

          if (.not.single_column ) then
             if (nlat_i /= lsmlat) then
                write(iulog,*)subname,' parameter lsmlat= ',lsmlat,'does not equal input nlat_i= ',nlat_i
                call endrun()
             end if
          endif
          call check_ret(nf90_inq_dimid(ncid, 'lsmpft', dimid), subname)
          call check_ret(nf90_inquire_dimension(ncid, dimid, len=npft_i), subname)

          if (.not. single_column ) then
             if (npft_i /= numpft+1) then
                write(iulog,*)subname,' parameter numpft+1 = ',numpft+1,'does not equal input npft_i= ',npft_i
                call endrun()
             end if
          endif
          call check_ret(nf90_inq_dimid(ncid, 'time', dimid), subname)
          call check_ret(nf90_inquire_dimension(ncid, dimid, len=ntim), subname)
 
       else

          nlon_i = lsmlon
          nlat_i = lsmlon
          npft_i = numpft+1

       endif   ! masterproc

       call mpi_bcast (nlon_i, 1, MPI_INTEGER, 0, mpicom, ier)
       call mpi_bcast (nlat_i, 1, MPI_INTEGER, 0, mpicom, ier)
       call mpi_bcast (npft_i, 1, MPI_INTEGER, 0, mpicom, ier)
       call mpi_bcast (ntim  , 1, MPI_INTEGER, 0, mpicom, ier)

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
       endif

       allocate(arrayl(begg:endg),stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname, 'allocation array error '; call endrun()
       end if

       do n = nps,npe
          beg4d(1) = 1         ; len4d(1) = nlon_i
          beg4d(2) = 1         ; len4d(2) = nlat_i
          beg4d(3) = n         ; len4d(3) = 1
          beg4d(4) = months(k) ; len4d(4) = 1
          if (single_column) then
             beg4d(1) = closelonidx; len4d(1) = 1
             beg4d(2) = closelatidx; len4d(2) = 1
          else
             beg4d(1) = 1         ; len4d(1) = nlon_i
             beg4d(2) = 1         ; len4d(2) = nlat_i
          end if
          if (nps == 0) beg4d(3) = n+1     ! shift by "1" if index starts at zero

          call ncd_iolocal(ncid,'MONTHLY_LAI','read',arrayl,grlnd,beg4d,len4d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_LAI NOT on fveg file' )
          mlai(begg:endg,n) = arrayl(begg:endg)

          call ncd_iolocal(ncid,'MONTHLY_SAI','read',arrayl,grlnd,beg4d,len4d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_SAI NOT on fveg file' )
          msai(begg:endg,n) = arrayl(begg:endg)

          call ncd_iolocal(ncid,'MONTHLY_HEIGHT_TOP','read',arrayl,grlnd,beg4d,len4d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_HEIGHT_TOP NOT on fveg file' )
          mhgtt(begg:endg,n) = arrayl(begg:endg)

          call ncd_iolocal(ncid,'MONTHLY_HEIGHT_BOT','read',arrayl,grlnd,beg4d,len4d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_HEIGHT_BOT NOT on fveg file' )
          mhgtb(begg:endg,n) = arrayl(begg:endg)
       enddo

       deallocate(arrayl)

       if (masterproc) then
          call check_ret(nf90_close(ncid), subname)
          write(iulog,*) 'Successfully read monthly vegetation data for'
          write(iulog,*) 'month ', months(k)
          write(iulog,*)
       end if


       ! store data directly in clmtype structure
       ! only vegetated pfts have nonzero values

       do p = begp,endp

          g = clm3%g%l%c%p%gridcell(p)

          ! Assign lai/sai/hgtt/hgtb to the top [maxpatch_pft] pfts
          ! as determined in subroutine surfrd


          ivt = clm3%g%l%c%p%itype(p)
          if (ivt /= noveg) then     ! vegetated pft
             do l = 0, numpft
                if (l == ivt) then
                   mlai2t(p,k) = mlai(g,l)
                   msai2t(p,k) = msai(g,l)
                   mhvt2t(p,k) = mhgtt(g,l)
                   mhvb2t(p,k) = mhgtb(g,l)
                end if
             end do
          else                        ! non-vegetated pft
             mlai2t(p,k) = 0._r8
             msai2t(p,k) = 0._r8
             mhvt2t(p,k) = 0._r8
             mhvb2t(p,k) = 0._r8
          end if

       end do   ! end of loop over pfts

    end do   ! end of loop over months

    deallocate(mlai, msai, mhgtt, mhgtb)

    mlaidiff => clm3%g%l%c%p%pps%mlaidiff
    do p = begp,endp
       mlaidiff(p)=mlai2t(p,1)-mlai2t(p,2)
    enddo

  end subroutine readMonthlyVegetation

!DML
!RL same as ReadMonthlyVegetation except for Soil Moisture
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readMonthlySoil
!
! !INTERFACE:
  subroutine readMonthlySoil (kmo, kda, months_soil)
!
! !DESCRIPTION:
! Read monthly soil moisture data for two consec. months.
!
! !USES:
    use clmtype
    use decompMod   , only : get_proc_bounds, ldecomp, gsmap_lnd_gdc2glo
    use clm_varpar  , only : nlevgrnd
    use clm_varcon  , only : istsoil

    use fileutils   , only : getfil
    use spmdMod     , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
    use clm_time_manager, only : get_nstep
    use ncdio       , only : check_ret,ncd_iolocal
    use netcdf
!
! !ARGUMENTS:
    implicit none

    integer, intent(in) :: kmo            ! month (1, ..., 12)
    integer, intent(in) :: kda            ! day of month (1, ..., 31)
    integer, intent(in) :: months_soil(2)      ! months to be interpolated (1 to 12)
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
    integer :: ncid,dimid,varid           ! input netCDF id's
    integer :: beg3d(3),len3d(3)          ! netCDF variable edges
    integer :: ntim                       ! number of input data time samples
    integer :: ncolumn_i                  ! number of columns
    integer :: nlev_i                     ! number of input data soil levels
    integer :: begc,endc                  ! beg and end local c index
    integer :: ier,ret                    ! error code
    integer :: closelatidx,closelonidx
    real(r8):: closelat,closelon

    real(r8), pointer :: mh2osoi_liq(:,:)  ! liquid soil water content read from input file
    real(r8), pointer :: mh2osoi_ice(:,:)  ! frozen soil water content read from input file
    real(r8), pointer :: arrayl(:)  ! temp local array
    character(len=32) :: subname = 'readMonthlySoil'
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype               => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit       => clm3%g%l%c%landunit

    ! Determine necessary indices


    call get_proc_bounds(begc=begc,endc=endc)

    allocate(mh2osoi_liq(begc:endc,1:nlevgrnd), &
             mh2osoi_ice(begc:endc,1:nlevgrnd), stat=ier)

    if (ier /= 0) then
       write(iulog,*)subname, 'allocation big error '; call endrun()
    end if

    ! ----------------------------------------------------------------------
    ! Open monthly soil moisture file
    ! Read data from column
    ! ----------------------------------------------------------------------

    do k=1,2   !loop over months and read vegetated data

       if (masterproc) then

          write(iulog,*) 'Attempting to read monthly soil moisture data .....'
          write(iulog,*) 'nstep = ',get_nstep(),' month = ',kmo,' day = ',kda

          call getfil('/ptmp/dlawren/SM_climatol/sm_control.h1.clim.0001.nc', locfn, 0)

          call check_ret(nf90_open(locfn, 0, ncid), subname)

          call check_ret(nf90_inq_dimid (ncid, 'column', dimid), subname)
          call check_ret(nf90_inquire_dimension(ncid, dimid, len=ncolumn_i), subname)
          call check_ret(nf90_inq_dimid(ncid, 'levgrnd', dimid), subname)
          call check_ret(nf90_inquire_dimension(ncid, dimid, len=nlev_i), subname)

          call check_ret(nf90_inq_dimid(ncid, 'time', dimid), subname)
          call check_ret(nf90_inquire_dimension(ncid, dimid, len=ntim), subname)

       endif   ! masterproc

       call mpi_bcast (ncolumn_i, 1, MPI_INTEGER, 0, mpicom, ier)

       allocate(arrayl(begc:endc),stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname, 'allocation array error '; call endrun()
       end if

       write(iulog,*) 'Before read '
       do j = 1,nlevgrnd
          beg3d(1) = j         ; len3d(1) = 1
          beg3d(2) = 1         ; len3d(2) = ncolumn_i
          beg3d(3) = months_soil(k) ; len3d(3) = 1

          call ncd_iolocal(ncid,'SOILLIQ','read',arrayl,namec,beg3d,len3d,status=ret)
          write(iulog,*) 'After read '
          if (ret /= 0) call endrun( trim(subname)//' ERROR: SOILLIQ NOT on clim file' )
          mh2osoi_liq(begc:endc,j) = arrayl(begc:endc)

          call ncd_iolocal(ncid,'SOILICE','read',arrayl,namec,beg3d,len3d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: SOILICE NOT on clim file' )
          mh2osoi_ice(begc:endc,j) = arrayl(begc:endc)

       enddo

       deallocate(arrayl)

       if (masterproc) then
          call check_ret(nf90_close(ncid), subname)
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

  end subroutine readMonthlySoil
!DML


end module STATICEcosysDynMod
