MODULE YOECND

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
REAL(KIND=JPRB),ALLOCATABLE:: CEVAPCU(:)

REAL(KIND=JPRB) :: REPFLM
REAL(KIND=JPRB) :: REPFLS
REAL(KIND=JPRB) :: REPQMI

!     -----------------------------------------------------------------
!*    CONTROL PARAMETERS FOR MOIST PROCESSES

! CEVAPCU(NFLEVG):
! REPFLM :  Minimum flux to avoid zero division in ice proportion
!           computations
! REPFLS :  Square-root of previous flux
! REPQMI :  Minimum specific humidity (security within QNEGAT)
!     -----------------------------------------------------------------

!$OMP THREADPRIVATE(repflm,repfls,repqmi)

!$OMP THREADPRIVATE(cevapcu)

END MODULE YOECND
