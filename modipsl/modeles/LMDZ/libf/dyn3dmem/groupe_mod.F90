MODULE groupe_mod

  REAL,POINTER,SAVE :: zconvm(:,:,:)
  REAL,POINTER,SAVE :: zconvmm(:,:,:)
  
CONTAINS

  SUBROUTINE groupe_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE infotrac
  USE advtrac_mod, ONLY : advtrac_allocate 
  IMPLICIT NONE
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  TYPE(distrib),POINTER :: d

    d=>distrib_caldyn
    CALL allocate2d_u(zconvm,llm,d)
    CALL allocate2d_u(zconvmm,llm,d)


  END SUBROUTINE groupe_allocate
  
  SUBROUTINE groupe_switch_caldyn(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL switch2d_u(zconvm,distrib_caldyn,dist)
    CALL switch2d_u(zconvmm,distrib_caldyn,dist)

  END SUBROUTINE groupe_switch_caldyn
  

  
END MODULE groupe_mod  
