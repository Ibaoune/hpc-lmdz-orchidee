MODULE gradiv2_mod

  REAL,POINTER,SAVE ::  gdx( :,: )
  REAL,POINTER,SAVE ::  gdy( :,: )
  REAL,POINTER,SAVE ::  div( :,: )
  
CONTAINS

  SUBROUTINE gradiv2_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
    TYPE(distrib),POINTER :: d
    d=>distrib_dissip

    CALL allocate_u(gdx,llm,d)
    CALL allocate_v(gdy,llm,d)
    CALL allocate_u(div,llm,d)

    
  END SUBROUTINE gradiv2_allocate
  
  SUBROUTINE gradiv2_switch_dissip(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL switch_u(gdx,distrib_dissip,dist)
    CALL switch_v(gdy,distrib_dissip,dist)
    CALL switch_u(div,distrib_dissip,dist)


  END SUBROUTINE gradiv2_switch_dissip
  
END MODULE gradiv2_mod  
