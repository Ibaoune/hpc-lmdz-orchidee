MODULE vlspltgen_mod

  REAL,POINTER,SAVE :: qsat(:,:)
  REAL,POINTER,SAVE :: mu(:,:) ! CRisi: on ajoute une dimension
  REAL,POINTER,SAVE :: mv(:,:)
  REAL,POINTER,SAVE :: mw(:,:,:)
  REAL,POINTER,SAVE :: zm(:,:,:)
  REAL,POINTER,SAVE :: zq(:,:,:)
 
CONTAINS

  SUBROUTINE vlspltgen_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE infotrac
  USE vlz_mod,ONLY : vlz_allocate 
  IMPLICIT NONE
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  TYPE(distrib),POINTER :: d
    
    d=>distrib_vanleer
    CALL allocate_u(qsat,llm,d)
    CALL allocate_u(mu,llm,d)
    CALL allocate_v(mv,llm,d)
    CALL allocate_u(mw,llm+1,nqtot,d)
    CALL allocate_u(zm,llm,nqtot,d)
    CALL allocate_u(zq,llm,nqtot,d)

    CALL vlz_allocate

  END SUBROUTINE vlspltgen_allocate
  
  SUBROUTINE vlspltgen_switch_vanleer(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  USE vlz_mod,ONLY : vlz_switch_vanleer 
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist
  
    CALL switch_u(qsat,distrib_vanleer,dist)
    CALL switch_u(mu,distrib_vanleer,dist)
    CALL switch_u(mv,distrib_vanleer,dist)
    CALL switch_u(mw,distrib_vanleer,dist)
    CALL switch_u(zm,distrib_vanleer,dist)
    CALL switch_u(zq,distrib_vanleer,dist)

    CALL vlz_switch_vanleer(dist)

  END SUBROUTINE vlspltgen_switch_vanleer  
  
END MODULE vlspltgen_mod  
