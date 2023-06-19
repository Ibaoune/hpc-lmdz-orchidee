MODULE YOMCT3

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Control variables for the model  (changed at level 3)

! NSTEP       : current timestep of model
! LGPQINSP    : .TRUE. if grid-point Q initialized from spectral field
! LSPC_FROM_DI: T (resp. F) if SPCH is called from STEPO (resp. STEPOTL).
! LRECALL_SUHDF_IN_CNT4: T if SUHDF must be called again in CNT4.

INTEGER(KIND=JPIM) :: NSTEP
LOGICAL :: LGPQINSP
LOGICAL :: LSPC_FROM_DI
LOGICAL :: LRECALL_SUHDF_IN_CNT4

!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(lgpqinsp,lrecall_suhdf_in_cnt4,lspc_from_di,nstep)
END MODULE YOMCT3
