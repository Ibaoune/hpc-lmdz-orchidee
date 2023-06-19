MODULE YOMJFH

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------

INTEGER(KIND=JPIM) :: N_VMASS=0

!    -----------------------------------------------------------------

!$OMP THREADPRIVATE(n_vmass)
END MODULE YOMJFH
