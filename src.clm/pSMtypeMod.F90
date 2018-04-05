module pSMtypeMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clmtype
!
! !DESCRIPTION: 
! Define derived type hierarchy for prescribing Soil Moisture.
! See clmtype.F90 and clmtypeInitMod.F90


  use shr_kind_mod,     only : r8 => shr_kind_r8
  use clm_varpar,       only : nlevgrnd, nlevsoi


! !PUBLIC TYPES:
  implicit none

  !PUBLIC MEMBER FUNCTIONS:
  public :: initPSMtype  ! Sets Soil Moisture

  public


  ! define prescribe_soilmoisture_type

  ! soilliq_prescribed and soilice_prescribed has the structure/ size of 
  ! column_wstate_type // cws

  ! reservoir has the structure/ size of
  ! column_wflux_type // cwf

  type, public :: prescribe_soilmoisture_type
     real(r8), pointer :: reservoir(:)             ! water reservoir for irrigation [mm]
     real(r8), pointer :: soilliq_prescribed(:, :) ! added SL at current timestep [mm]
     real(r8), pointer :: soilice_prescribed(:, :) ! added SI at current timestep [mm]
  end type prescribe_soilmoisture_type


type(prescribe_soilmoisture_type) :: psm




!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initClmtype
!
! !INTERFACE:

  subroutine initPSMtype()

    use decompMod,      only : get_proc_bounds

    integer :: begc, endc   ! per-proc beginning and ending column indices

    call get_proc_bounds(begc=begc,endc=endc)


    allocate(psm%reservoir(begc:endc))
    allocate(psm%soilliq_prescribed(begc:endc,1:nlevgrnd))
    allocate(psm%soilice_prescribed(begc:endc,1:nlevgrnd))

    psm%reservoir(begc:endc) = 0._r8
    psm%soilliq_prescribed(begc:endc,1:nlevgrnd) = 0._r8
    psm%soilice_prescribed(begc:endc,1:nlevgrnd) = 0._r8


  end subroutine initPSMtype

end module pSMtypeMod