MODULE YOERRTO6

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO6* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 6
!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO6  = 16

REAL(KIND=JPRB) , DIMENSION(NO6) :: FRACREFAO

REAL(KIND=JPRB) , DIMENSION(NO6) :: CFC11ADJO
REAL(KIND=JPRB) , DIMENSION(NO6) :: CFC12O
REAL(KIND=JPRB) , DIMENSION(NO6) :: ABSCO2O

REAL(KIND=JPRB) :: KAO(5,13,NO6)
REAL(KIND=JPRB) :: SELFREFO(10,NO6)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
!$OMP THREADPRIVATE(absco2o,cfc11adjo,cfc12o,fracrefao,kao,selfrefo)
END MODULE YOERRTO6
