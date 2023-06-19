MODULE YOEDBUG

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -------------------------------------------------
!*    ** *YOEDBUG* - CONTROL OPTIONS FOR DEBUGGING HELP
!     -------------------------------------------------
INTEGER(KIND=JPIM) :: KSTPDBG(3)
!     ------------------------------------------------------------------

!$OMP THREADPRIVATE(kstpdbg)

END MODULE YOEDBUG

