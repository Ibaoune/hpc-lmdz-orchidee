MODULE YOERRTO15

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO15* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 15
!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO15 = 16

REAL(KIND=JPRB) :: FRACREFAO(NO15,9)

REAL(KIND=JPRB) :: KAO(9,5,13,NO15)
REAL(KIND=JPRB) :: SELFREFO(10,NO15)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL 
!     -----------------------------------------------------------------
!$OMP THREADPRIVATE(fracrefao,kao,selfrefo)
END MODULE YOERRTO15
