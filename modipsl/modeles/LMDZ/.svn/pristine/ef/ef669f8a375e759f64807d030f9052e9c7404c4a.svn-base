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

  USE phys_state_var_mod, ONLY: rlat, rlon, phys_state_var_init
  USE mod_phys_lmdz_para, ONLY: klon_omp
  USE geometry_mod, ONLY: longitude_deg, latitude_deg
  IMPLICIT NONE
      
  INTEGER,INTENT(IN) :: nlon,iflag_phys

  ! local variables
  REAL :: pi

  ! initializations:
  pi=2.*ASIN(1.)

  CALL phys_state_var_init()

  rlat(1:klon_omp)=latitude_deg(1:klon_omp)
  rlon(1:klon_omp)=longitude_deg(1:klon_omp)


  ! Here you could create an initial condition for the physics
  ! ...
  ! ... fill in the fields...
  ! ...
  ! ... and create a "startphy.nc" file
      CALL phyredem ("startphy.nc")

  END SUBROUTINE iniaqua

END MODULE phyaqua_mod
