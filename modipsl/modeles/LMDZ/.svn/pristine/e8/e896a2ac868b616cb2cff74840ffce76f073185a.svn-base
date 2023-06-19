!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/geo2atm.F90,v 1.1 2008-12-05 17:56:40 lsce Exp $
!
SUBROUTINE geo2atm(im, jm, px, py, pz, plon, plat, pu, pv, pr)
  USE dimphy
  USE mod_phys_lmdz_para
  USE mod_grid_phy_lmdz, only: grid_type, unstructured, regular_lonlat
  IMPLICIT NONE
  INCLUDE 'YOMCST.h'
  CHARACTER (len = 6)                :: clmodnam
  CHARACTER (len = 20)               :: modname = 'geo2atm'
  CHARACTER (len = 80)               :: abort_message

! Change wind coordinates from cartesian geocentric to local spherical
! NB! Fonctionne probablement uniquement en MPI seul (sans OpenMP)
!
! Geocentric :
  ! axe x is eastward : crosses (0N,90E) point.
  ! axe y  crosses (0N,180E) point.
  ! axe z is 'up' : crosses north pole.
!
! NB! Aux poles, fonctionne probablement uniquement en MPI seul (sans OpenMP)

  INTEGER, INTENT (IN)                 :: im, jm
  REAL, DIMENSION (im,jm), INTENT(IN)  :: px, py, pz
  REAL, DIMENSION (im,jm), INTENT(IN)  :: plon, plat
  REAL, DIMENSION (im,jm), INTENT(OUT) :: pu, pv, pr

  REAL :: rad,reps


  rad = rpi / 180.0E0
  reps = 1.0e-5

  pu(:,:) = &
       - px(:,:) * SIN(rad * plon(:,:)) &
       + py(:,:) * COS(rad * plon(:,:))

  pv(:,:) = &
       - px(:,:) * SIN(rad * plat(:,:)) * COS(rad * plon(:,:)) &
       - py(:,:) * SIN(rad * plat(:,:)) * SIN(rad * plon(:,:)) &
       + pz(:,:) * COS(rad * plat(:,:))  

  pr(:,:) = &
       + px(:,:) * COS(rad * plat(:,:)) * COS(rad * plon(:,:)) &
       + py(:,:) * COS(rad * plat(:,:)) * SIN(rad * plon(:,:)) &
       + pz(:,:) * SIN(rad * plat(:,:))

  IF (grid_type==regular_lonlat) THEN
    ! Value at North Pole
    IF (is_north_pole_dyn) THEN
       pu(:, 1) = -px (1,1)
       pv(:, 1) = -py (1,1)
       pr(:, 1) = 0.0
    ENDIF
  
    ! Value at South Pole     
    IF (is_south_pole_dyn) THEN
       pu(:,jm) = -px (1,jm)
       pv(:,jm) = -py (1,jm)
       pr(:,jm) = 0.0
    ENDIF

  ELSE IF (grid_type==unstructured) THEN
     ! Pole nord pour Dynamico
     WHERE ( plat(:,:) >= 90.0-reps )
        pu(:,:) =  py(:,:)
        pv(:,:) = -px(:,:)
        pr(:,:) = 0.0e0
     END WHERE

  ELSE
     abort_message='Problem: unknown grid type'
     CALL abort_physic(modname,abort_message,1)
  END IF

  
  
  
END SUBROUTINE geo2atm
