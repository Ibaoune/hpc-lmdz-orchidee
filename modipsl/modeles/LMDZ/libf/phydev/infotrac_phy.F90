
! $Id: $

MODULE infotrac_phy

! Information on tracers for physics;
! nqtot : total number of tracers and higher order of moment, water vapor and liquid included
  INTEGER, SAVE :: nqtot
!$OMP THREADPRIVATE(nqtot)

  CHARACTER(len=4),SAVE :: type_trac
!$OMP THREADPRIVATE(type_trac)

CONTAINS

  SUBROUTINE init_infotrac_phy(nqtot_,type_trac_)
  ! transfer information on tracers from dynamics to physics
  IMPLICIT NONE
    INTEGER,INTENT(IN) :: nqtot_
    CHARACTER(len=4),INTENT(IN) :: type_trac_

    nqtot=nqtot_
    type_trac=type_trac_
  
  END SUBROUTINE init_infotrac_phy

END MODULE infotrac_phy
