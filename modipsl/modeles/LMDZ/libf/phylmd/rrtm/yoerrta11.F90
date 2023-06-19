MODULE YOERRTA11

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA11* - RRTM COEFFICIENTS FOR INTERVAL 11
!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG11 = 8

REAL(KIND=JPRB) , DIMENSION(NG11) :: FRACREFA
REAL(KIND=JPRB) , DIMENSION(NG11) :: FRACREFB

REAL(KIND=JPRB) :: KA(5,13,NG11)   , ABSA(65,NG11)
REAL(KIND=JPRB) :: KB(5,13:59,NG11), ABSB(235,NG11)
REAL(KIND=JPRB) :: SELFREF(10,NG11)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTA11
