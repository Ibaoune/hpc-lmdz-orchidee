MODULE YOERRTA5

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA5* - RRTM COEFFICIENTS FOR INTERVAL 5
!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG5  = 16

REAL(KIND=JPRB) :: FRACREFA(NG5,9) ,FRACREFB(NG5,5)

REAL(KIND=JPRB) , DIMENSION(NG5) :: CCL4

REAL(KIND=JPRB) :: KA(9,5,13,NG5)   ,ABSA(585,NG5)
REAL(KIND=JPRB) :: KB(5,5,13:59,NG5),ABSB(1175,NG5)
REAL(KIND=JPRB) :: SELFREF(10,NG5)
REAL(KIND=JPRB) :: STRRAT1
REAL(KIND=JPRB) :: STRRAT2

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! CCL4    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
! STRRAT1 : REAL
! STRRAT2 : REAL    
!     -----------------------------------------------------------------
END MODULE YOERRTA5
