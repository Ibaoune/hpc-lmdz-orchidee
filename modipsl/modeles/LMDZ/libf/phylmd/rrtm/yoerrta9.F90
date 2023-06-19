MODULE YOERRTA9

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA9* - RRTM COEFFICIENTS FOR INTERVAL 9
!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG9  = 12

REAL(KIND=JPRB) :: FRACREFA(NG9,9)

REAL(KIND=JPRB) , DIMENSION(NG9) :: FRACREFB
REAL(KIND=JPRB) , DIMENSION(13) :: N2OREF
REAL(KIND=JPRB) , DIMENSION(13) :: H2OREF
REAL(KIND=JPRB) , DIMENSION(13) :: CH4REF
REAL(KIND=JPRB) , DIMENSION(11) :: ETAREF
! 36 = 3*NG9      
REAL(KIND=JPRB) , DIMENSION(36) :: ABSN2O

REAL(KIND=JPRB) :: KA(11,5,13,NG9) ,ABSA(715,NG9)
REAL(KIND=JPRB) :: KB(5,13:59,NG9) ,ABSB(235,NG9)
REAL(KIND=JPRB) :: SELFREF(10,NG9)
REAL(KIND=JPRB) :: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! ABSN2O  : REAL    
! CH4REF  : REAL
! ETAREF  : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! H2OREF  : REAL
! N2OREF  : REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
! STRRAT  : REAL
!     -----------------------------------------------------------------
END MODULE YOERRTA9
