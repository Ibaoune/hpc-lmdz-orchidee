MODULE caldyn_mod

  REAL,POINTER,SAVE :: vcont(:,:)
  REAL,POINTER,SAVE :: ucont(:,:)
  REAL,POINTER,SAVE :: ang(:,:)
  REAL,POINTER,SAVE :: p(:,:)
  REAL,POINTER,SAVE :: massebx(:,:)
  REAL,POINTER,SAVE :: masseby(:,:)
  REAL,POINTER,SAVE :: psexbarxy(:,:)
  REAL,POINTER,SAVE :: vorpot(:,:)
  REAL,POINTER,SAVE :: ecin(:,:)
  REAL,POINTER,SAVE :: bern(:,:)
  REAL,POINTER,SAVE :: massebxy(:,:)
  REAL,POINTER,SAVE :: convm(:,:)


  
CONTAINS

  SUBROUTINE caldyn_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  USE advect_new_mod,ONLY : advect_new_allocate
  IMPLICIT NONE
  TYPE(distrib),POINTER :: d


    d=>distrib_caldyn
    CALL allocate_v(vcont,llm,d)
    CALL allocate_u(ucont,llm,d)
    CALL allocate_u(ang,llm,d)
    CALL allocate_u(p,llmp1,d)
    CALL allocate_u(massebx,llm,d)
    CALL allocate_v(masseby,llm,d)
    CALL allocate_v(psexbarxy,llm,d)
    CALL allocate_v(vorpot,llm,d)
    CALL allocate_u(ecin,llm,d)
    CALL allocate_u(bern,llm,d)
    CALL allocate_v(massebxy,llm,d)
    CALL allocate_u(convm,llm,d)
    
    CALL advect_new_allocate
    
  END SUBROUTINE caldyn_allocate
  
  SUBROUTINE caldyn_switch_caldyn(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  USE advect_new_mod,ONLY : advect_new_switch_caldyn
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL switch_v(vcont,distrib_caldyn,dist)
    CALL switch_u(ucont,distrib_caldyn,dist)
    CALL switch_u(ang,distrib_caldyn,dist)
    CALL switch_u(p,distrib_caldyn,dist)
    CALL switch_u(massebx,distrib_caldyn,dist)
    CALL switch_v(masseby,distrib_caldyn,dist)
    CALL switch_v(psexbarxy,distrib_caldyn,dist)
    CALL switch_v(vorpot,distrib_caldyn,dist)
    CALL switch_u(ecin,distrib_caldyn,dist)
    CALL switch_u(bern,distrib_caldyn,dist)
    CALL switch_v(massebxy,distrib_caldyn,dist)
    CALL switch_u(convm,distrib_caldyn,dist)
    
    CALL advect_new_switch_caldyn(dist)
    
  END SUBROUTINE caldyn_switch_caldyn
  

  
END MODULE caldyn_mod  
