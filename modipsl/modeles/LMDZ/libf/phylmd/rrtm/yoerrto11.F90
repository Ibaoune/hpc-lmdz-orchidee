MODULE YOERRTO11

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO11* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 11
!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO11 = 16

REAL(KIND=JPRB) , DIMENSION(NO11) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO11) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(5,13,NO11)
REAL(KIND=JPRB) :: KBO(5,13:59,NO11)
REAL(KIND=JPRB) :: SELFREFO(10,NO11)

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

!$OMP THREADPRIVATE(fracrefao,fracrefbo,kao,kbo,selfrefo)
END MODULE YOERRTO11
