MODULE integrd_mod

  REAL,POINTER,SAVE :: p(:,:)
  REAL,POINTER,SAVE :: deltap(:,:)
  REAL,POINTER,SAVE :: ps(:)


  
CONTAINS

  SUBROUTINE integrd_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  USE advect_new_mod,ONLY : advect_new_allocate
  IMPLICIT NONE
  TYPE(distrib),POINTER :: d


    d=>distrib_caldyn
    CALL allocate_u(p,llmp1,d)
    CALL allocate_u(deltap,llm,d)
    CALL allocate_u(ps,d)
    ps(:)=0

    
  END SUBROUTINE integrd_allocate
  
  SUBROUTINE integrd_switch_caldyn(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL switch_u(p,distrib_caldyn,dist)
    CALL switch_u(deltap,distrib_caldyn,dist)
    CALL switch_u(ps,distrib_caldyn,dist)

    
    
  END SUBROUTINE integrd_switch_caldyn
  

  
END MODULE integrd_mod  
