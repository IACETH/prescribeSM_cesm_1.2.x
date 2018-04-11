module pSMhistFldsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: pSMhistFldsMod
!
! !DESCRIPTION:
! Module containing initialization of clm history concerning pSM
! We use a new file because it was not possible to easily merge with
! Wim Thierrys energy balance adaptations of 
!
! !USES:
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public hist_init_pSM_Flds ! Build field list of all fields relevant for pSM

!
!EOP
!------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hist_init_pSM_Flds
!
! !INTERFACE:
  subroutine hist_init_pSM_Flds()
!
! !DESCRIPTION:
! Build master field list of all possible fields in a history file.
! Each field has associated with it a ``long\_name'' netcdf attribute that
! describes what the field is, and a ``units'' attribute. A subroutine is
! called to add each field to the masterlist.
!
! !USES:
    use pSMtypeMod
    use histFileMod, only : hist_addfld1d, hist_addfld2d

!
! !ARGUMENTS:
    implicit none

    call hist_addfld1d (fname='RESERVOIR', units='mm', &
         avgflag='A', long_name='water stored in the reservoir', &
         ptr_col=psm%reservoir, set_lake=0._r8)

    call hist_addfld2d (fname='SOILLIQ_PRES',  units='kg/m2', type2d='levgrnd', &
         avgflag='A', long_name='prescribed soil liquid water (vegetated landunits only)', &
         ptr_col=psm%soilliq_prescribed, l2g_scale_type='veg')

    call hist_addfld2d (fname='SOILICE_PRES',  units='kg/m2', type2d='levgrnd', &
         avgflag='A', long_name='prescribed soil ice (vegetated landunits only)', &
         ptr_col=psm%soilice_prescribed, l2g_scale_type='veg')

  end subroutine hist_init_pSM_Flds

end module pSMhistFldsMod