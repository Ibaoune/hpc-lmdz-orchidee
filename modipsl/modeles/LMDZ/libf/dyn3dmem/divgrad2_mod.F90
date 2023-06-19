MODULE divgrad2_mod

  REAL,POINTER,SAVE ::  divgra( :,: )
  
CONTAINS

  SUBROUTINE divgrad2_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
    TYPE(distrib),POINTER :: d
    d=>distrib_dissip

    CALL allocate_u(divgra,llm,d)

    
  END SUBROUTINE divgrad2_allocate
  
  SUBROUTINE divgrad2_switch_dissip(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL switch_u(divgra,distrib_dissip,dist)


  END SUBROUTINE divgrad2_switch_dissip
  
END MODULE divgrad2_mod  
