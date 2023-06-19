
! $Id: $

MODULE infotrac_phy

! Infotrac for physics; for now contains the same information as infotrac for
! the dynamics (could be further cleaned) and is initialized using values
! provided by the dynamics

! nqtot : total number of tracers and higher order of moment, water vapor and liquid included
  INTEGER, SAVE :: nqtot
!$OMP THREADPRIVATE(nqtot)
 
CONTAINS

  SUBROUTINE init_infotrac_phy(nqtot_)
  ! transfer information on tracers from dynamics to physics
  USE print_control_mod, ONLY: prt_level, lunout
  IMPLICIT NONE
    INTEGER,INTENT(IN) :: nqtot_

    CHARACTER(LEN=30) :: modname="init_infotrac_phy"

    nqtot=nqtot_
  
    IF(prt_level.ge.1) THEN
      write(lunout,*) TRIM(modname)//": nqtot",nqtot
    ENDIF
    
  END SUBROUTINE init_infotrac_phy

END MODULE infotrac_phy
