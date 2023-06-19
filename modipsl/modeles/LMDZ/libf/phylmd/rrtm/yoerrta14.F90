MODULE YOERRTA14

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA14* - RRTM COEFFICIENTS FOR INTERVAL 14
!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG14 = 2

REAL(KIND=JPRB) , DIMENSION(NG14) :: FRACREFA
REAL(KIND=JPRB) , DIMENSION(NG14) :: FRACREFB

REAL(KIND=JPRB) :: KA(5,13,NG14)   ,ABSA(65,NG14)
REAL(KIND=JPRB) :: KB(5,13:59,NG14),ABSB(235,NG14)
REAL(KIND=JPRB) :: SELFREF(10,NG14)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTA14

