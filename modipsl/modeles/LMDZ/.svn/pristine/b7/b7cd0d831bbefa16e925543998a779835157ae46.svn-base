MODULE YOERRTO7

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO7* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 7
!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO7  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO7,9)

REAL(KIND=JPRB) , DIMENSION(NO7) :: FRACREFBO
REAL(KIND=JPRB) , DIMENSION(NO7) :: ABSCO2O
REAL(KIND=JPRB) :: KAO(9,5,13,NO7)
REAL(KIND=JPRB) :: KBO(5,13:59,NO7)
REAL(KIND=JPRB) :: SELFREFO(10,NO7)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL  
!     -----------------------------------------------------------------
!$OMP THREADPRIVATE(absco2o,fracrefao,fracrefbo,kao,kbo,selfrefo)
END MODULE YOERRTO7
