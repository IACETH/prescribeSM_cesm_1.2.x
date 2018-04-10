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

! !PUBLIC TYPES:
  implicit none

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

end module pSMtypeMod
  