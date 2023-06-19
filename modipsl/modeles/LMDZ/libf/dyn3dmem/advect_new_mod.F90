MODULE advect_new_mod

  REAL,POINTER,SAVE :: dv1(:,:)
  REAL,POINTER,SAVE :: du1(:,:)
  REAL,POINTER,SAVE :: dteta1(:,:)
  REAL,POINTER,SAVE :: dv2(:,:)
  REAL,POINTER,SAVE :: du2(:,:)
  REAL,POINTER,SAVE :: dteta2(:,:)
  REAL,POINTER,SAVE :: uav(:,:)
  REAL,POINTER,SAVE :: vav(:,:)

  
CONTAINS

  SUBROUTINE advect_new_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
  TYPE(distrib),POINTER :: d


    d=>distrib_caldyn
    CALL allocate_v(dv1,llm,d)
    CALL allocate_u(du1,llm,d)
    CALL allocate_u(dteta1,llm,d)
    CALL allocate_v(dv2,llm,d)
    CALL allocate_u(du2,llm,d)
    CALL allocate_u(dteta2,llm,d)
    CALL allocate_u(uav,llm,d)
    CALL allocate_v(vav,llm,d)
 
    
  END SUBROUTINE advect_new_allocate
  
  SUBROUTINE advect_new_switch_caldyn(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL switch_v(dv1,distrib_caldyn,dist)
    CALL switch_u(du1,distrib_caldyn,dist)
    CALL switch_u(dteta1,distrib_caldyn,dist)
    CALL switch_v(dv2,distrib_caldyn,dist)
    CALL switch_u(du2,distrib_caldyn,dist)
    CALL switch_u(dteta2,distrib_caldyn,dist)
    CALL switch_u(uav,distrib_caldyn,dist)
    CALL switch_v(vav,distrib_caldyn,dist)

  END SUBROUTINE advect_new_switch_caldyn
  
END MODULE advect_new_mod  
