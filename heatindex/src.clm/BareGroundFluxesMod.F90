module BareGroundFluxesMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: BareGroundFluxesMod
!
! !DESCRIPTION:
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
!KO
   use shr_const_mod      , only : SHR_CONST_TKFRZ
!KO
!
! !PUBLIC TYPES:
   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: BareGroundFluxes   ! Calculate sensible and latent heat fluxes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BareGroundFluxes
!
! !INTERFACE:
  subroutine BareGroundFluxes(lbp, ubp, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! !USES:
    use clmtype
    use clm_atmlnd         , only : clm_a2l
    use clm_varpar         , only : nlevgrnd
    use clm_varcon         , only : cpair, vkc, grav, denice, denh2o, istsoil
    use clm_varcon         , only : istcrop
    use shr_const_mod      , only : SHR_CONST_RGAS
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use QSatMod            , only : QSat

!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                     ! pft bounds
    integer, intent(in) :: num_nolakep                  ! number of pft non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)    ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 12/19/01, Peter Thornton
! This routine originally had a long list of parameters, and also a reference to
! the entire clm derived type.  For consistency, only the derived type reference
! is passed (now pointing to the current column and pft), and the other original
! parameters are initialized locally. Using t_grnd instead of tg (tg eliminated
! as redundant).
! 1/23/02, PET: Added pft reference as parameter. All outputs will be written
! to the pft data structures, and averaged to the column level outside of
! this routine.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)        ! pft's column index
    integer , pointer :: pgridcell(:)      ! pft's gridcell index
    integer , pointer :: plandunit(:)      ! pft's landunit index
    integer , pointer :: ltype(:)          ! landunit type
    integer , pointer :: frac_veg_nosno(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
    real(r8), pointer :: t_grnd(:)         ! ground surface temperature [K]
    real(r8), pointer :: thm(:)            ! intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
    real(r8), pointer :: qg(:)             ! specific humidity at ground surface [kg/kg]
    real(r8), pointer :: thv(:)            ! virtual potential temperature (kelvin)
    real(r8), pointer :: dqgdT(:)          ! temperature derivative of "qg"
    real(r8), pointer :: htvp(:)           ! latent heat of evaporation (/sublimation) [J/kg]
    real(r8), pointer :: beta(:)           ! coefficient of conective velocity [-]
    real(r8), pointer :: zii(:)            ! convective boundary height [m]
    real(r8), pointer :: forc_u(:)         ! atmospheric wind speed in east direction (m/s)
    real(r8), pointer :: forc_v(:)         ! atmospheric wind speed in north direction (m/s)
    real(r8), pointer :: forc_t(:)         ! atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_th(:)        ! atmospheric potential temperature (Kelvin)
    real(r8), pointer :: forc_q(:)         ! atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_rho(:)       ! density (kg/m**3)
    real(r8), pointer :: forc_pbot(:)      ! atmospheric pressure (Pa)
    real(r8), pointer :: forc_hgt_u_pft(:) ! observational height of wind at pft level [m]
    real(r8), pointer :: psnsun(:)         ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: psnsha(:)         ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: z0mg_col(:)       ! roughness length, momentum [m]
    real(r8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
    real(r8), pointer :: dz(:,:)           ! layer depth (m)
    real(r8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: frac_sno(:)       ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: soilbeta(:)       ! soil wetness relative to field capacity
!KO
    real(r8), pointer :: u10_clm(:)        ! 10-m wind (m/s)
    real(r8), pointer :: u10_clm_r(:)      ! rural 10-m wind (m/s)
!KO
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: z0hg_col(:)       ! roughness length, sensible heat [m]
    real(r8), pointer :: z0qg_col(:)       ! roughness length, latent heat [m]
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: dlrad(:)         ! downward longwave radiation below the canopy [W/m2]
    real(r8), pointer :: ulrad(:)         ! upward longwave radiation above the canopy [W/m2]
    real(r8), pointer :: cgrnds(:)        ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    real(r8), pointer :: cgrndl(:)        ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    real(r8), pointer :: cgrnd(:)         ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(r8), pointer :: taux(:)          ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8), pointer :: tauy(:)          ! wind (shear) stress: n-s (kg/m/s**2)
    real(r8), pointer :: eflx_sh_grnd(:)  ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot(:)   ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi(:) ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot(:) ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(r8), pointer :: t_ref2m(:)       ! 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: q_ref2m(:)       ! 2 m height surface specific humidity (kg/kg)
    real(r8), pointer :: t_ref2m_r(:)     ! Rural 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: rh_ref2m_r(:)    ! Rural 2 m height surface relative humidity (%)
    real(r8), pointer :: rh_ref2m(:)      ! 2 m height surface relative humidity (%)
    real(r8), pointer :: t_veg(:)         ! vegetation temperature (Kelvin)
    real(r8), pointer :: btran(:)         ! transpiration wetness factor (0 to 1)
    real(r8), pointer :: rssun(:)         ! sunlit stomatal resistance (s/m)
    real(r8), pointer :: rssha(:)         ! shaded stomatal resistance (s/m)
    real(r8), pointer :: ram1(:)          ! aerodynamical resistance (s/m)
    real(r8), pointer :: fpsn(:)          ! photosynthesis (umol CO2 /m**2 /s)
    real(r8), pointer :: rootr(:,:)       ! effective fraction of roots in each soil layer
    real(r8), pointer :: rresis(:,:)      ! root resistance by layer (0-1)  (nlevgrnd)	
!KO
    real(r8), pointer :: nws_heatindex_ref2m(:)   ! 2 m height NWS heat index (degrees F)
    real(r8), pointer :: appar_temp_ref2m(:)      ! 2 m height apparent temperature (degrees C)
    real(r8), pointer :: swbgt_ref2m(:)           ! 2 m height simplified wet bulb global temperature (degrees C)
    real(r8), pointer :: humidex_ref2m(:)         ! 2 m height humidex (degrees C)
    real(r8), pointer :: discomf_index_ref2m(:)   ! 2 m height discomfort index (degrees C)
    real(r8), pointer :: nws_heatindex_ref2m_r(:) ! Rural 2 m height NWS heat index (degrees F)
    real(r8), pointer :: appar_temp_ref2m_r(:)    ! Rural 2 m height apparent temperature (degrees C)
    real(r8), pointer :: swbgt_ref2m_r(:)         ! Rural 2 m height simplified wet bulb global temperature (degrees C)
    real(r8), pointer :: humidex_ref2m_r(:)       ! Rural 2 m height humidex (degrees C)
    real(r8), pointer :: discomf_index_ref2m_r(:) ! Rural 2 m height discomfort index (degrees C)
!KO
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer, parameter  :: niters = 3  ! maximum number of iterations for surface temperature
    integer  :: p,c,g,f,j,l            ! indices
    integer  :: filterp(ubp-lbp+1)     ! pft filter for vegetated pfts
    integer  :: fn                     ! number of values in local pft filter
    integer  :: fp                     ! lake filter pft index
    integer  :: iter                   ! iteration index
    real(r8) :: zldis(lbp:ubp)         ! reference height "minus" zero displacement height [m]
    real(r8) :: displa(lbp:ubp)        ! displacement height [m]
    real(r8) :: zeta                   ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: wc                     ! convective velocity [m/s]
    real(r8) :: dth(lbp:ubp)           ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                   ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh(lbp:ubp)           ! diff of humidity between ref. height and surface
    real(r8) :: obu(lbp:ubp)           ! Monin-Obukhov length (m)
    real(r8) :: ur(lbp:ubp)            ! wind speed at reference height [m/s]
    real(r8) :: um(lbp:ubp)            ! wind speed including the stablity effect [m/s]
    real(r8) :: temp1(lbp:ubp)         ! relation for potential temperature profile
    real(r8) :: temp12m(lbp:ubp)       ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(lbp:ubp)         ! relation for specific humidity profile
    real(r8) :: temp22m(lbp:ubp)       ! relation for specific humidity profile applied at 2-m
    real(r8) :: ustar(lbp:ubp)         ! friction velocity [m/s]
    real(r8) :: tstar                  ! temperature scaling parameter
    real(r8) :: qstar                  ! moisture scaling parameter
    real(r8) :: thvstar                ! virtual potential temperature scaling parameter
    real(r8) :: cf                     ! heat transfer coefficient from leaves [-]
    real(r8) :: ram                    ! aerodynamical resistance [s/m]
    real(r8) :: rah                    ! thermal resistance [s/m]
    real(r8) :: raw                    ! moisture resistance [s/m]
    real(r8) :: raih                   ! temporary variable [kg/m2/s]
    real(r8) :: raiw                   ! temporary variable [kg/m2/s]
    real(r8) :: fm(lbp:ubp)            ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: z0mg_pft(lbp:ubp)
    real(r8) :: z0hg_pft(lbp:ubp)
    real(r8) :: z0qg_pft(lbp:ubp)
    real(r8) :: e_ref2m                ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                 ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m             ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT              ! derivative of 2 m height surface saturated specific humidity on t_ref2m 
    real(r8) :: www                    ! surface soil wetness [-]
!KO
    real(r8) :: e                      ! vapor pressure for calculation of apparent temperature (hPa)
    real(r8) :: wbt                    ! wet-bulb temperature for calculation of discomfort index (degrees C)
    real(r8) :: t_ref2m_C              ! 2-m height surface air temperature (degrees C)
    real(r8) :: t_ref2m_F              ! 2-m height surface air temperature (degrees F)
    real(r8) :: rh                     ! 2-m height surface relative humidity for calculation of wet-bulb temperature (%)
    real(r8) :: rh_min                 ! minimum 2-m height surface relative humidity for calculation of wet-bulb temperature (%)
!KO
!------------------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    forc_u     => clm_a2l%forc_u
    forc_v     => clm_a2l%forc_v

    ! Assign local pointers to derived type members (landunit-level)

    ltype      => clm3%g%l%itype

    ! Assign local pointers to derived type members (column-level)

    forc_th    => clm3%g%l%c%ces%forc_th
    forc_t     => clm3%g%l%c%ces%forc_t
    forc_pbot  => clm3%g%l%c%cps%forc_pbot
    forc_rho   => clm3%g%l%c%cps%forc_rho
    forc_q     => clm3%g%l%c%cws%forc_q
    pcolumn    => clm3%g%l%c%p%column
    pgridcell  => clm3%g%l%c%p%gridcell
    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    dlrad  => clm3%g%l%c%p%pef%dlrad
    ulrad  => clm3%g%l%c%p%pef%ulrad
    t_grnd => clm3%g%l%c%ces%t_grnd
    qg     => clm3%g%l%c%cws%qg
    z0mg_col => clm3%g%l%c%cps%z0mg
    z0hg_col => clm3%g%l%c%cps%z0hg
    z0qg_col => clm3%g%l%c%cps%z0qg
    thv    => clm3%g%l%c%ces%thv
    beta   => clm3%g%l%c%cps%beta
    zii    => clm3%g%l%c%cps%zii
    ram1   => clm3%g%l%c%p%pps%ram1
    cgrnds => clm3%g%l%c%p%pef%cgrnds
    cgrndl => clm3%g%l%c%p%pef%cgrndl
    cgrnd  => clm3%g%l%c%p%pef%cgrnd
    dqgdT  => clm3%g%l%c%cws%dqgdT
    htvp   => clm3%g%l%c%cps%htvp
    watsat         => clm3%g%l%c%cps%watsat
    h2osoi_ice     => clm3%g%l%c%cws%h2osoi_ice
    dz             => clm3%g%l%c%cps%dz
    h2osoi_liq     => clm3%g%l%c%cws%h2osoi_liq
    frac_sno       => clm3%g%l%c%cps%frac_sno
    soilbeta       => clm3%g%l%c%cws%soilbeta

    ! Assign local pointers to derived type members (pft-level)

    taux => clm3%g%l%c%p%pmf%taux
    tauy => clm3%g%l%c%p%pmf%tauy
    eflx_sh_grnd => clm3%g%l%c%p%pef%eflx_sh_grnd
    eflx_sh_tot => clm3%g%l%c%p%pef%eflx_sh_tot
    qflx_evap_soi => clm3%g%l%c%p%pwf%qflx_evap_soi
    qflx_evap_tot => clm3%g%l%c%p%pwf%qflx_evap_tot
    t_ref2m => clm3%g%l%c%p%pes%t_ref2m
    q_ref2m => clm3%g%l%c%p%pes%q_ref2m
    t_ref2m_r => clm3%g%l%c%p%pes%t_ref2m_r
    rh_ref2m_r => clm3%g%l%c%p%pes%rh_ref2m_r
    plandunit => clm3%g%l%c%p%landunit
    rh_ref2m => clm3%g%l%c%p%pes%rh_ref2m
    t_veg => clm3%g%l%c%p%pes%t_veg
    thm => clm3%g%l%c%p%pes%thm
    btran => clm3%g%l%c%p%pps%btran
    rssun => clm3%g%l%c%p%pps%rssun
    rssha => clm3%g%l%c%p%pps%rssha
    rootr => clm3%g%l%c%p%pps%rootr
    rresis => clm3%g%l%c%p%pps%rresis
    psnsun => clm3%g%l%c%p%pcf%psnsun
    psnsha => clm3%g%l%c%p%pcf%psnsha
    fpsn => clm3%g%l%c%p%pcf%fpsn
    forc_hgt_u_pft => clm3%g%l%c%p%pps%forc_hgt_u_pft
!KO
    u10_clm               => clm3%g%l%c%p%pps%u10_clm
    nws_heatindex_ref2m   => clm3%g%l%c%p%pes%nws_heatindex_ref2m
    appar_temp_ref2m      => clm3%g%l%c%p%pes%appar_temp_ref2m
    swbgt_ref2m           => clm3%g%l%c%p%pes%swbgt_ref2m
    humidex_ref2m         => clm3%g%l%c%p%pes%humidex_ref2m
    discomf_index_ref2m   => clm3%g%l%c%p%pes%discomf_index_ref2m
    u10_clm_r             => clm3%g%l%c%p%pps%u10_clm_r
    nws_heatindex_ref2m_r => clm3%g%l%c%p%pes%nws_heatindex_ref2m_r
    appar_temp_ref2m_r    => clm3%g%l%c%p%pes%appar_temp_ref2m_r
    swbgt_ref2m_r         => clm3%g%l%c%p%pes%swbgt_ref2m_r
    humidex_ref2m_r       => clm3%g%l%c%p%pes%humidex_ref2m_r
    discomf_index_ref2m_r => clm3%g%l%c%p%pes%discomf_index_ref2m_r
!KO

    ! Filter pfts where frac_veg_nosno is zero

    fn = 0
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       if (frac_veg_nosno(p) == 0) then
          fn = fn + 1
          filterp(fn) = p
       end if
    end do

    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step

    do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       g = pgridcell(p)

       ! Initialization variables

       displa(p) = 0._r8
       dlrad(p)  = 0._r8
       ulrad(p)  = 0._r8

       ur(p) = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
       dth(p) = thm(p)-t_grnd(c)
       dqh(p) = forc_q(c) - qg(c)
       dthv = dth(p)*(1._r8+0.61_r8*forc_q(c))+0.61_r8*forc_th(c)*dqh(p)
       zldis(p) = forc_hgt_u_pft(p)

       ! Copy column roughness to local pft-level arrays

       z0mg_pft(p) = z0mg_col(c)
       z0hg_pft(p) = z0hg_col(c)
       z0qg_pft(p) = z0qg_col(c)

       ! Initialize Monin-Obukhov length and wind speed

       call MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg_pft(p), um(p), obu(p))

    end do

    ! Perform stability iteration
    ! Determine friction velocity, and potential temperature and humidity
    ! profiles of the surface boundary layer

    do iter = 1, niters

       call FrictionVelocity(lbp, ubp, fn, filterp, &
                             displa, z0mg_pft, z0hg_pft, z0qg_pft, &
                             obu, iter, ur, um, ustar, &
                             temp1, temp2, temp12m, temp22m, fm)

       do f = 1, fn
          p = filterp(f)
          c = pcolumn(p)
          g = pgridcell(p)

          tstar = temp1(p)*dth(p)
          qstar = temp2(p)*dqh(p)
          z0hg_pft(p) = z0mg_pft(p)/exp(0.13_r8 * (ustar(p)*z0mg_pft(p)/1.5e-5_r8)**0.45_r8)
          z0qg_pft(p) = z0hg_pft(p)
          thvstar = tstar*(1._r8+0.61_r8*forc_q(c)) + 0.61_r8*forc_th(c)*qstar
          zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

          if (zeta >= 0._r8) then                   !stable
             zeta = min(2._r8,max(zeta,0.01_r8))
             um(p) = max(ur(p),0.1_r8)
          else                                      !unstable
             zeta = max(-100._r8,min(zeta,-0.01_r8))
             wc = beta(c)*(-grav*ustar(p)*thvstar*zii(c)/thv(c))**0.333_r8
             um(p) = sqrt(ur(p)*ur(p) + wc*wc)
          end if
          obu(p) = zldis(p)/zeta
       end do

    end do ! end stability iteration

     do j = 1, nlevgrnd
       do f = 1, fn
          p = filterp(f)
          rootr(p,j) = 0._r8
          rresis(p,j) = 0._r8
        end do
     end do

    do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       g = pgridcell(p)
       l = plandunit(p)

       ! Determine aerodynamic resistances

       ram     = 1._r8/(ustar(p)*ustar(p)/um(p))
       rah     = 1._r8/(temp1(p)*ustar(p))
       raw     = 1._r8/(temp2(p)*ustar(p))
       raih    = forc_rho(c)*cpair/rah

       ! Soil evaporation resistance
       www     = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)/watsat(c,1)
       www     = min(max(www,0.0_r8),1._r8)

       !changed by K.Sakaguchi. Soilbeta is used for evaporation
       if (dqh(p) .gt. 0._r8) then   !dew  (beta is not applied, just like rsoil used to be)
          raiw    = forc_rho(c)/(raw)
       else
       ! Lee and Pielke 1992 beta is applied
          raiw    = soilbeta(c)*forc_rho(c)/(raw)
       end if

       ram1(p) = ram  !pass value to global variable

       ! Output to pft-level data structures
       ! Derivative of fluxes with respect to ground temperature

       cgrnds(p) = raih
       cgrndl(p) = raiw*dqgdT(c)
       cgrnd(p)  = cgrnds(p) + htvp(c)*cgrndl(p)

       ! Surface fluxes of momentum, sensible and latent heat
       ! using ground temperatures from previous time step

       taux(p)          = -forc_rho(c)*forc_u(g)/ram
       tauy(p)          = -forc_rho(c)*forc_v(g)/ram
       eflx_sh_grnd(p)  = -raih*dth(p)
       eflx_sh_tot(p)   = eflx_sh_grnd(p)
       qflx_evap_soi(p) = -raiw*dqh(p)
       qflx_evap_tot(p) = qflx_evap_soi(p)

       ! 2 m height air temperature

       t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

       ! 2 m height specific humidity

       q_ref2m(p) = forc_q(c) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

       ! 2 m height relative humidity
                                                                                
       call QSat(t_ref2m(p), forc_pbot(c), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)

       rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)
!KO
       ! National Weather Service Heat Index: Lans P. Rothfusz. "The heat index 'equation' (or
       ! more than you ever wanted to know about heat index)", Scientific Services Division (NWS
       ! Southern Region Headquarters), 1 July 1990
       ! Requires air temperature (F), relative humidity (%)
       ! Valid for air temperatures above 20C. If below this set heatindex to air temperature.
       t_ref2m_F = (t_ref2m(p) - SHR_CONST_TKFRZ) * 9._r8/5._r8 + 32._r8
       if ( (t_ref2m(p) - SHR_CONST_TKFRZ) .lt. 20._r8) then
         nws_heatindex_ref2m(p) = t_ref2m_F
       else
         nws_heatindex_ref2m(p) = -42.379_r8 + 2.04901523_r8*t_ref2m_F + 10.14333127_r8*rh_ref2m(p) + &
                                  (-0.22475541_r8*t_ref2m_F*rh_ref2m(p)) + &
                                  (-6.83783e-3_r8*t_ref2m_F**2._r8) + &
                                  (-5.481717e-2_r8*rh_ref2m(p)**2._r8) + &
                                  1.22874e-3_r8*(t_ref2m_F**2._r8)*rh_ref2m(p) + &
                                  8.5282e-4_r8*t_ref2m_F*rh_ref2m(p)**2._r8 + &
                                  (-1.99e-6_r8*(t_ref2m_F**2._r8)*(rh_ref2m(p)**2._r8))
       end if

       ! Apparent Temperature (Australian BOM): Steadman, R.G., 1994: Norms of apparent temperature
       ! in Australia, Aust. Met. Mag., 43, 1-16. Here we use equation 22 where AT is a function of
       ! air temperature (C), water vapor pressure (kPa), and 10-m wind speed (m/s).
       t_ref2m_C = t_ref2m(p) - SHR_CONST_TKFRZ
       ! e from Erich Fischer (consistent with CLM equations)
       e = (rh_ref2m(p)/100._r8) *e_ref2m   ! Pa
       appar_temp_ref2m(p) = t_ref2m_C + 3.30_r8*e/1000._r8 - 0.70_r8*u10_clm(p) - 4.0_r8

       ! Simplified Wet Bulb Globe Temperature: Willett, K.M., and S. Sherwood, 2010: Exceedance of heat
       ! index thresholds for 15 regions under a warming climate using the wet-bulb globe temperature,
       ! Int. J. Climatol., doi:10.1002/joc.2257
       ! Requires air temperature (C), water vapor pressure (hPa)
       swbgt_ref2m(p) = 0.567_r8*t_ref2m_C  + 0.393_r8*e/100._r8 + 3.94_r8

       ! Humidex: Masterson, J., and F. Richardson, 1979: Humidex, a method of quantifying human discomfort
       ! due to excessive heat and humidity, CLI 1-79, Environment Canada, Atmosheric Environment Service,
       ! Downsview, Ontario
       ! Requires air temperature (C), water vapor pressure (hPa)
       humidex_ref2m(p) = t_ref2m_C + ((5._r8/9._r8) * (e/100._r8 - 10._r8))

       ! Discomfort Index: Epstein, Y., and D.S. Moran, 2006: Thermal comfort and the heat stress indices,
       ! Ind. Health, 44, 388-398.
       ! The wet bulb temperature is from: Stull, R., 2011: Wet-bulb temperature from relative humidity
       ! and air temperature, J. Appl. Meteor. Climatol., doi:10.1175/JAMC-D-11-0143.1
       ! Requires air temperature (C), wet bulb temperature (C) (which in turn requires air temperature (C) and
       ! relative humidity (%))
       t_ref2m_C = min(t_ref2m_C,50._r8)
       rh = min(rh_ref2m(p),99._r8)
       rh = max(rh,5._r8)
       rh_min = t_ref2m_C*(-2.27_r8)+27.7_r8
       if (t_ref2m_C .lt. -20._r8 .or. rh .lt. rh_min) then
         discomf_index_ref2m(p) = t_ref2m_C
       else
         wbt = t_ref2m_C * atan(0.151977_r8*sqrt(rh + 8.313659_r8)) + &
               atan(t_ref2m_C+rh) - atan(rh-1.676331_r8) + &
               0.00391838_r8*rh**(3._r8/2._r8)*atan(0.023101_r8*rh) - &
               4.686035_r8
         discomf_index_ref2m(p) = 0.5_r8*wbt + 0.5_r8*t_ref2m_C
       end if
!KO

       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
         rh_ref2m_r(p) = rh_ref2m(p)
         t_ref2m_r(p) = t_ref2m(p)
!KO
         u10_clm_r(p) = u10_clm(p)
         nws_heatindex_ref2m_r(p) = nws_heatindex_ref2m(p)
         appar_temp_ref2m_r(p) = appar_temp_ref2m(p)
         swbgt_ref2m_r(p) = swbgt_ref2m(p)
         humidex_ref2m_r(p) = humidex_ref2m(p)
         discomf_index_ref2m_r(p) = discomf_index_ref2m(p)
!KO
       end if

       ! Variables needed by history tape

       t_veg(p) = forc_t(c)
       btran(p) = 0._r8
       cf = forc_pbot(c)/(SHR_CONST_RGAS*0.001_r8*thm(p))*1.e06_r8
       rssun(p) = 1._r8/1.e15_r8 * cf
       rssha(p) = 1._r8/1.e15_r8 * cf

       ! Add the following to avoid NaN

       psnsun(p) = 0._r8
       psnsha(p) = 0._r8
       fpsn(p) = 0._r8
       clm3%g%l%c%p%pps%lncsun(p) = 0._r8
       clm3%g%l%c%p%pps%lncsha(p) = 0._r8
       clm3%g%l%c%p%pps%vcmxsun(p) = 0._r8
       clm3%g%l%c%p%pps%vcmxsha(p) = 0._r8
       ! adding code for isotopes, 8/17/05, PET
       clm3%g%l%c%p%pps%cisun(p) = 0._r8
       clm3%g%l%c%p%pps%cisha(p) = 0._r8
#if (defined C13)
       clm3%g%l%c%p%pps%alphapsnsun(p) = 0._r8
       clm3%g%l%c%p%pps%alphapsnsha(p) = 0._r8
       clm3%g%l%c%p%pepv%rc13_canair(p) = 0._r8
       clm3%g%l%c%p%pepv%rc13_psnsun(p) = 0._r8
       clm3%g%l%c%p%pepv%rc13_psnsha(p) = 0._r8
       clm3%g%l%c%p%pc13f%psnsun(p) = 0._r8
       clm3%g%l%c%p%pc13f%psnsha(p) = 0._r8
       clm3%g%l%c%p%pc13f%fpsn(p) = 0._r8
#endif

    end do

  end subroutine BareGroundFluxes

end module BareGroundFluxesMod
