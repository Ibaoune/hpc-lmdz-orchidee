MODULE YOERRTO8

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO8* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 8
!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO8  = 16

REAL(KIND=JPRB) , DIMENSION(NO8) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO8) :: FRACREFBO
REAL(KIND=JPRB) , DIMENSION(NO8) :: CFC12O
REAL(KIND=JPRB) , DIMENSION(NO8) :: CFC22ADJO
REAL(KIND=JPRB) , DIMENSION(NO8) :: ABSCO2AO
REAL(KIND=JPRB) , DIMENSION(NO8) :: ABSCO2BO
REAL(KIND=JPRB) , DIMENSION(NO8) :: ABSN2OAO
REAL(KIND=JPRB) , DIMENSION(NO8) :: ABSN2OBO

REAL(KIND=JPRB) :: KAO(5,7,NO8)
REAL(KIND=JPRB) :: KBO(5,7:59,NO8)
REAL(KIND=JPRB) :: SELFREFO(10,NO8)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSCO2A : REAL     
! ABSCO2B : REAL     
! ABSN2OA : REAL     
! ABSN2OB : REAL 
! CFC12   : REAL     
! CFC22ADJ: REAL     
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
!$OMP THREADPRIVATE(absco2ao,absco2bo,absn2oao,absn2obo,cfc12o,cfc22adjo,fracrefao,fracrefbo,kao,kbo,selfrefo)
END MODULE YOERRTO8
