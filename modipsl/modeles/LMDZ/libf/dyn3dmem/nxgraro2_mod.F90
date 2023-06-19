MODULE nxgraro2_mod

  REAL,POINTER,SAVE ::  grx( :,: )
  REAL,POINTER,SAVE ::  gry( :,: )
  REAL,POINTER,SAVE ::  rot( :,: )
  
CONTAINS

  SUBROUTINE nxgraro2_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
    TYPE(distrib),POINTER :: d
    d=>distrib_dissip

    CALL allocate_u(grx,llm,d)
    CALL allocate_v(gry,llm,d)
    CALL allocate_v(rot,llm,d)

    
  END SUBROUTINE nxgraro2_allocate
  
  SUBROUTINE nxgraro2_switch_dissip(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL switch_u(grx,distrib_dissip,dist)
    CALL switch_v(gry,distrib_dissip,dist)
    CALL switch_v(rot,distrib_dissip,dist)


  END SUBROUTINE nxgraro2_switch_dissip
  
END MODULE nxgraro2_mod  
