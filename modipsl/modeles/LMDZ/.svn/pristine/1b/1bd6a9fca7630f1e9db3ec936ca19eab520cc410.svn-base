MODULE caladvtrac_mod

  REAL,POINTER,SAVE :: q_adv(:,:,:)
  REAL,POINTER,SAVE :: massem_adv(:,:)
  REAL,POINTER,SAVE :: wg_adv(:,:)
  REAL,POINTER,SAVE :: teta_adv(:,:)
  REAL,POINTER,SAVE :: p_adv(:,:)
  REAL,POINTER,SAVE :: pk_adv(:,:)
  REAL,POINTER,SAVE :: pbarug_adv(:,:)
  REAL,POINTER,SAVE :: pbarvg_adv(:,:)
  REAL,POINTER,SAVE :: pbaruc(:,:)
  REAL,POINTER,SAVE :: pbarvc(:,:)
  REAL,POINTER,SAVE :: pbarug(:,:)
  REAL,POINTER,SAVE :: pbarvg(:,:)
  REAL,POINTER,SAVE :: wg(:,:)

  REAL,POINTER,SAVE :: massem(:,:)
  
CONTAINS

  SUBROUTINE caladvtrac_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE infotrac
  USE advtrac_mod, ONLY : advtrac_allocate
  USE groupe_mod 
  IMPLICIT NONE
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  TYPE(distrib),POINTER :: d

    d=>distrib_vanleer
    CALL allocate_u(q_adv,llm,nqtot,d)
    CALL allocate_u(massem_adv,llm,d)
    CALL allocate_u(wg_adv,llm,d)
    CALL allocate_u(teta_adv,llm,d)
    CALL allocate_u(p_adv,llmp1,d)
    CALL allocate_u(pk_adv,llm,d)
    CALL allocate_u(pbarug_adv,llm,d)
    CALL allocate_v(pbarvg_adv,llm,d)

    d=>distrib_caldyn
    CALL allocate_u(massem,llm,d)
    CALL allocate_u(pbaruc,llm,d)
    pbaruc(:,:)=0
    CALL allocate_v(pbarvc,llm,d)
    pbarvc(:,:)=0
    CALL allocate_u(pbarug,llm,d)
    CALL allocate_v(pbarvg,llm,d)
    CALL allocate_u(wg,llm,d)

    CALL groupe_allocate
    CALL advtrac_allocate
    
  END SUBROUTINE caladvtrac_allocate
  
  SUBROUTINE caladvtrac_switch_caldyn(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  USE groupe_mod
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL switch_u(massem,distrib_caldyn,dist)
    CALL switch_u(pbaruc,distrib_caldyn,dist)
    CALL switch_v(pbarvc,distrib_caldyn,dist,up=1)
    CALL switch_u(pbarug,distrib_caldyn,dist)
    CALL switch_v(pbarvg,distrib_caldyn,dist)
    CALL switch_u(wg,distrib_caldyn,dist)
    
    CALL groupe_switch_caldyn(dist)

  END SUBROUTINE caladvtrac_switch_caldyn
  
  SUBROUTINE caladvtrac_switch_vanleer(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  USE advtrac_mod, ONLY : advtrac_switch_vanleer 
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist
  
    CALL switch_u(q_adv,distrib_vanleer,dist)
    CALL switch_u(massem_adv,distrib_vanleer,dist)
    CALL switch_u(wg_adv,distrib_vanleer,dist)
    CALL switch_u(teta_adv,distrib_vanleer,dist)
    CALL switch_u(p_adv,distrib_vanleer,dist)
    CALL switch_u(pk_adv,distrib_vanleer,dist)
    CALL switch_u(pbarug_adv,distrib_vanleer,dist)
    CALL switch_v(pbarvg_adv,distrib_vanleer,dist)

    CALL advtrac_switch_vanleer(dist)
    
  END SUBROUTINE caladvtrac_switch_vanleer  
  
END MODULE caladvtrac_mod  
