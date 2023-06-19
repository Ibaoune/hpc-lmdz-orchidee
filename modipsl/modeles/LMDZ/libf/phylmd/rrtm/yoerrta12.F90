MODULE YOERRTA12

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA12* - RRTM COEFFICIENTS FOR INTERVAL 12
!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG12 = 8

REAL(KIND=JPRB) :: FRACREFA(NG12,9)
REAL(KIND=JPRB) :: KA(9,5,13,NG12) ,ABSA(585,NG12)
REAL(KIND=JPRB) :: SELFREF(10,NG12)

REAL(KIND=JPRB) :: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL
! STRRAT1 : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTA12
