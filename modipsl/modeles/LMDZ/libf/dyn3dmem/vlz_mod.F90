MODULE vlz_mod

  REAL,POINTER,SAVE :: wq(:,:,:)
  REAL,POINTER,SAVE :: dzq(:,:)
  REAL,POINTER,SAVE :: dzqw(:,:)
  REAL,POINTER,SAVE :: adzqw(:,:)
  ! CRisi: pour les traceurs:  
  !REAL,POINTER,SAVE :: masseq(:,:,:)
  REAL,POINTER,SAVE :: Ratio(:,:,:)
  
CONTAINS

  SUBROUTINE vlz_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE infotrac
  USE dimensions_mod
  IMPLICIT NONE
  TYPE(distrib),POINTER :: d
    
    d=>distrib_vanleer
    CALL allocate_u(wq,llm+1,nqtot,d)
    CALL allocate_u(dzq,llm,d)
    CALL allocate_u(dzqw,llm,d)
    CALL allocate_u(adzqw,llm,d)
    if (nqdesc_tot.gt.0) then
    !CALL allocate_u(masseq,llm,nqtot,d)
    CALL allocate_u(Ratio,llm,nqtot,d)
    endif !if (nqdesc_tot.gt.0) then

  END SUBROUTINE vlz_allocate
  
  SUBROUTINE vlz_switch_vanleer(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  USE infotrac
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist
  
    CALL switch_u(wq,distrib_vanleer,dist)
    CALL switch_u(dzq,distrib_vanleer,dist)
    CALL switch_u(dzqw,distrib_vanleer,dist)
    CALL switch_u(adzqw,distrib_vanleer,dist)
    ! CRisi:
    if (nqdesc_tot.gt.0) then    
    !CALL switch_u(masseq,distrib_vanleer,dist)
    CALL switch_u(Ratio,distrib_vanleer,dist)
    endif !if (nqdesc_tot.gt.0) then     

  END SUBROUTINE vlz_switch_vanleer  
  
END MODULE vlz_mod  
