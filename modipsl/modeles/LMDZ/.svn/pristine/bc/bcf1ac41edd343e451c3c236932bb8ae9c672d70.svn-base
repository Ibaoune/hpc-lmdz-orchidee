MODULE advtrac_mod

  REAL,POINTER,SAVE :: finmasse(:,:)
  
CONTAINS

  SUBROUTINE advtrac_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE infotrac
  USE vlspltgen_mod
  IMPLICIT NONE
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  TYPE(distrib),POINTER :: d
    
    d=>distrib_vanleer
    CALL allocate_u(finmasse,llm,d)
    CALL vlspltgen_allocate
  END SUBROUTINE advtrac_allocate
  
  SUBROUTINE advtrac_switch_vanleer(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  USE vlspltgen_mod
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist
  
    CALL switch_u(finmasse,distrib_vanleer,dist)

    CALL vlspltgen_switch_vanleer(dist)

  END SUBROUTINE advtrac_switch_vanleer  
  
END MODULE advtrac_mod  
