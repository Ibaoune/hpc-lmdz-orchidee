MODULE dissip_mod


  
CONTAINS

  SUBROUTINE dissip_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  USE gradiv2_mod, ONLY : gradiv2_allocate
  USE nxgraro2_mod, ONLY : nxgraro2_allocate
  USE divgrad2_mod, ONLY : divgrad2_allocate
  IMPLICIT NONE

    CALL gradiv2_allocate
    CALL nxgraro2_allocate
    CALL divgrad2_allocate

    
  END SUBROUTINE dissip_allocate
  
  SUBROUTINE dissip_switch_dissip(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  USE gradiv2_mod,ONLY : gradiv2_switch_dissip
  USE nxgraro2_mod,ONLY : nxgraro2_switch_dissip
  USE divgrad2_mod,ONLY : divgrad2_switch_dissip
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL gradiv2_switch_dissip(dist)
    CALL nxgraro2_switch_dissip(dist)
    CALL divgrad2_switch_dissip(dist)
    
  END SUBROUTINE dissip_switch_dissip
  
END MODULE dissip_mod  
