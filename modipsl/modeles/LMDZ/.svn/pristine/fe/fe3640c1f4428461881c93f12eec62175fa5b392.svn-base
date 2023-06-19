MODULE YOERRTA3

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA3* - RRTM COEFFICIENTS FOR INTERVAL 3
!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG3  = 16

REAL(KIND=JPRB) :: FRACREFA(NG3,10) ,FRACREFB(NG3,5)

REAL(KIND=JPRB) , DIMENSION(16) :: FORREF
REAL(KIND=JPRB) , DIMENSION(16) :: ABSN2OA
REAL(KIND=JPRB) , DIMENSION(16) :: ABSN2OB
REAL(KIND=JPRB) , DIMENSION(10) :: ETAREF
REAL(KIND=JPRB) , DIMENSION(59) :: H2OREF
REAL(KIND=JPRB) , DIMENSION(59) :: N2OREF
REAL(KIND=JPRB) , DIMENSION(59) :: CO2REF

REAL(KIND=JPRB) :: KA(10,5,13,NG3)  ,ABSA(650,NG3)
REAL(KIND=JPRB) :: KB(5,5,13:59,NG3),ABSB(1175,NG3)
REAL(KIND=JPRB) :: SELFREF(10,NG3)
REAL(KIND=JPRB) :: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,13,1),ABSB(1,1))

!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! ABSN2OA : REAL
! ABSN2OB : REAL
! CO2REF  : REAL
! ETAREF  : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! H2OREF  : REAL
! KA      : REAL     
! KB      : REAL     
! N2OREF  : REAL
! SELFREF : REAL     
! STRRAT  : REAL
!     -----------------------------------------------------------------
END MODULE YOERRTA3
