module pSMtypeInitMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clmtype
!
! !DESCRIPTION: 
! Define derived type hierarchy for prescribing Soil Moisture.
! See clmtype.F90 and clmtypeInitMod.F90


  use shr_kind_mod,     only : r8 => shr_kind_r8
  use clmtype
  use pSMtypeMod
  use clm_varpar,       only : nlevgrnd, nlevsoi

! !PUBLIC TYPES:
  implicit none
  save

  !PUBLIC MEMBER FUNCTIONS:
  public :: initPSMtype  ! Sets Soil Moisture



  contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initClmtype
!
! !INTERFACE:

  subroutine initPSMtype()
    
    ! USES:
    use decompMod,      only : get_proc_bounds
    
    ! LOCAL VARAIBLES:
    integer :: begc, endc   ! per-proc beginning and ending column indices

    call get_proc_bounds(begc=begc,endc=endc)

    allocate(psm%reservoir(begc:endc))
    allocate(psm%soilliq_prescribed(begc:endc,1:nlevgrnd))
    allocate(psm%soilice_prescribed(begc:endc,1:nlevgrnd))

    psm%reservoir(begc:endc) = 0._r8
    psm%soilliq_prescribed(begc:endc,1:nlevgrnd) = 0._r8
    psm%soilice_prescribed(begc:endc,1:nlevgrnd) = 0._r8

    ! set irrigation flux to 0
    cwf%qflx_irrig(begc:endc)  = 0._r8

  end subroutine initPSMtype

end module pSMtypeInitMod
