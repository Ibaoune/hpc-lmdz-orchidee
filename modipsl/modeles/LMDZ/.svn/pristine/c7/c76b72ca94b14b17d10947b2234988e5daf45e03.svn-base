SUBROUTINE alpale_wk ( dtime, cell_area, zoccur, sigmaw, wdens, fip ,  &
                       fip_cond)

! **************************************************************
!                                                              *
! ALPALE_WK                                                    *
!                                                              *
!                                                              *
! written by   : Jean-Yves Grandpeix, 07/08/2017               *
! modified by :                                                *
! **************************************************************

  USE dimphy, ONLY: klon
  USE ioipsl_getin_p_mod, ONLY : getin_p
  USE print_control_mod, ONLY: mydebug=>debug , lunout, prt_level
!
  IMPLICIT NONE

!================================================================
! Auteur(s)   : Jean-Yves Grandpeix, 07/08/2017
! Objet : Contribution of the wake scheme to Ale and Alp
!================================================================

! Input arguments
!----------------
  REAL, INTENT(IN)                                           :: dtime
  REAL, DIMENSION(klon),    INTENT(IN)                       :: cell_area
  INTEGER, DIMENSION(klon), INTENT (IN)                      :: zoccur
  REAL, DIMENSION(klon),    INTENT(IN)                       :: sigmaw
  REAL, DIMENSION(klon),    INTENT(IN)                       :: wdens
  REAL, DIMENSION(klon),    INTENT(IN)                       :: fip
! Output arguments
!-----------------
  REAL, DIMENSION(klon), INTENT(OUT)                         :: fip_cond


! Local variables
!----------------
  INTEGER                                                    :: i
  LOGICAL, SAVE                                              :: first = .TRUE.
  !$OMP THREADPRIVATE(first)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)                      :: cellrad
  !$OMP THREADPRIVATE(cellrad)
  REAL, DIMENSION(klon)                                      :: wkrad
  REAL, DIMENSION(klon)                                      :: proba_gf

  INCLUDE "YOMCST.h"   ! rpi

IF (first) THEN
  ALLOCATE (cellrad(klon))
!  Compute pseudo grid-cell radius cellrad, such that pi*cellrad^2=cell_area
  print *,'alpale_wk: cell_area(1) ',cell_area(1)
  cellrad(:)=sqrt(cell_area(:)/rpi)
  first = .FALSE.
ENDIF

!  Compute wake radius
!!  print *,'alpale_wk: sigmaw(1), wdens(1) ', sigmaw(1), wdens(1)
  DO i = 1,klon
    IF (zoccur(i) .GE. 1) THEN
      wkrad(i) = sqrt(sigmaw(i)/(rpi*wdens(i)))
    ELSE
      wkrad(i) = 0.
    ENDIF ! (zoccur(i) .GE. 1)
  ENDDO

!  Compute probability that the grid-cell is intersected by a gust front
!!  print *,'alpale_wk: wkrad(1), cellrad(1) ', wkrad(1), cellrad(1)
!!  proba_gf(:) = exp(-wdens(:)*rpi*max(wkrad(:)-cellrad(:),0.)**2) - &   ! Formules
!!                exp(-wdens(:)*rpi*(wkrad(:)+cellrad(:))**2)             ! fausses !
  proba_gf(:) = 1. - exp(-wdens(:)*rpi*((wkrad(:)+cellrad(:))**2 - &
                                        max(wkrad(:)-cellrad(:),0.)**2) )
!
  proba_gf(:) = max(proba_gf(:),1.e-3)
!  Compute Fip conditionned on the presence of some gust front within the 
!  grid-cell
!!  print *,'alpale_wk: proba_gf(1), fip(1), ', proba_gf(1), fip(1)
  fip_cond(:) = fip(:)/proba_gf(:)
!!    print *,'alpale_wk: wkrad(1), cellrad(1), proba_gf(1), fip(1), fip_cond(1) ', &
!!                        wkrad(1), cellrad(1), proba_gf(1), fip(1), fip_cond(1)

   RETURN
   END SUBROUTINE alpale_wk

