!
! $Id: $
!
MODULE phyaqua_mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE iniaqua(nlon, iflag_phys)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Create an initial state (startphy.nc) for the physics
  !  Usefull for idealised cases (e.g. aquaplanets or testcases)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE phys_state_var_mod, ONLY: phys_state_var_init
  USE mod_phys_lmdz_para, ONLY: klon_omp
  IMPLICIT NONE
      
  INTEGER,INTENT(IN) :: nlon,iflag_phys

  CALL phys_state_var_init()


  ! Here you could create an initial condition for the physics
  ! ...
  ! ... fill in the fields...
  ! ...
  ! ... and create a "startphy.nc" file
      CALL phyredem ("startphy.nc")

  END SUBROUTINE iniaqua

END MODULE phyaqua_mod
