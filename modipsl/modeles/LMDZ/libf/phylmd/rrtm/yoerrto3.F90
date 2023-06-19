MODULE YOERRTO3

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO3* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 3
!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO3  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO3,10) ,FRACREFBO(NO3,5)

REAL(KIND=JPRB) , DIMENSION(NO3) :: FORREFO
REAL(KIND=JPRB) , DIMENSION(NO3) :: ABSN2OAO
REAL(KIND=JPRB) , DIMENSION(NO3) :: ABSN2OBO

REAL(KIND=JPRB) :: KAO(10,5,13,NO3)
REAL(KIND=JPRB) :: KBO(5,5,13:59,NO3)
REAL(KIND=JPRB) :: SELFREFO(10,NO3)

!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSN2OAO: REAL
! ABSN2OBO: REAL
!FRACREFAO: REAL    
!FRACREFBO: REAL
! KAO     : REAL     
! KBO     : REAL     
! SELFREFO: REAL     
!     -----------------------------------------------------------------
!$OMP THREADPRIVATE(absn2oao,absn2obo,forrefo,fracrefao,fracrefbo,kao,kbo,selfrefo)
END MODULE YOERRTO3
