MODULE YOERRTO1

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO1* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 1
!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO1  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO1)  , FRACREFBO(NO1)
REAL(KIND=JPRB) :: KAO(5,13,NO1)
REAL(KIND=JPRB) :: KBO(5,13:59,NO1)
REAL(KIND=JPRB) :: SELFREFO(10,NO1), FORREFO(NO1)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
!FRACREFAO: REAL    
!FRACREFBO: REAL
! FORREFO : REAL
! KAO     : REAL     
! KBO     : REAL     
! SELFREFO: REAL     
!     -----------------------------------------------------------------
!$OMP THREADPRIVATE(forrefo,fracrefao,fracrefbo,kao,kbo,selfrefo)
END MODULE YOERRTO1
