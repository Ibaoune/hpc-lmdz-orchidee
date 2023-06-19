MODULE YOERRTA6

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA6* - RRTM COEFFICIENTS FOR INTERVAL 6
!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG6  = 8

REAL(KIND=JPRB) , DIMENSION(NG6) :: FRACREFA

REAL(KIND=JPRB) , DIMENSION(NG6) :: CFC11ADJ
REAL(KIND=JPRB) , DIMENSION(NG6) :: CFC12
REAL(KIND=JPRB) , DIMENSION(NG6) :: ABSCO2

REAL(KIND=JPRB) :: KA(5,13,NG6),ABSA(65,NG6)
REAL(KIND=JPRB) :: SELFREF(10,NG6)

EQUIVALENCE (KA(1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSCO2  : REAL 
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! CFC11ADJ: REAL
! CFC12   : REAL
! KA      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTA6
