MODULE YOERRTA8

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA8* - RRTM COEFFICIENTS FOR INTERVAL 8
!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG8  = 8

REAL(KIND=JPRB) , DIMENSION(NG8) :: FRACREFA
REAL(KIND=JPRB) , DIMENSION(NG8) :: FRACREFB
REAL(KIND=JPRB) , DIMENSION(NG8) :: CFC12
REAL(KIND=JPRB) , DIMENSION(NG8) :: CFC22ADJ
REAL(KIND=JPRB) , DIMENSION(NG8) :: ABSCO2A
REAL(KIND=JPRB) , DIMENSION(NG8) :: ABSCO2B
REAL(KIND=JPRB) , DIMENSION(NG8) :: ABSN2OA
REAL(KIND=JPRB) , DIMENSION(NG8) :: ABSN2OB
REAL(KIND=JPRB) , DIMENSION(59)  :: H2OREF
REAL(KIND=JPRB) , DIMENSION(59)  :: N2OREF
REAL(KIND=JPRB) , DIMENSION(59)  :: O3REF

REAL(KIND=JPRB) :: KA(5,7,NG8)    ,ABSA(35,NG8)
REAL(KIND=JPRB) :: KB(5,7:59,NG8) ,ABSB(265,NG8)
REAL(KIND=JPRB) :: SELFREF(10,NG8)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,7,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! ABSCO2A : REAL     
! ABSCO2B : REAL     
! ABSN2OA : REAL     
! ABSN2OB : REAL 
! CFC12   : REAL     
! CFC22ADJ: REAL     
! FRACREFA: REAL    
! FRACREFB: REAL    
! H2OREF  : REAL    
! KA      : REAL     
! KB      : REAL     
! N2OREF  : REAL    
! O3REF   : REAL    
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTA8
