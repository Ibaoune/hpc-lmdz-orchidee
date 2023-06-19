MODULE YOMSLPHY

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! * Variables for split ECMWF physics.

! LSLPHY  : for split physics (one part at t-Dt, one part at t+Dt);
!           can be used only at ECMWF with the ECMWF package.
! RSLWX   : level of implicitness of semi-Lagrangian/physics.
! NVTEND  : third dimension of SAVTEND (number of 3D fields).
! SAVTEND : buffer to store the physical tendencies.

LOGICAL :: LSLPHY
REAL(KIND=JPRB)   , PARAMETER :: RSLWX=0.5_JPRB
INTEGER(KIND=JPIM) :: NVTEND
REAL(KIND=JPRB),ALLOCATABLE :: SAVTEND(:,:,:,:)
! Pointers for SAVTEND
INTEGER(KIND=JPIM) :: MU_SAVTEND,MV_SAVTEND,MT_SAVTEND,MSAT_SAVTEND
INTEGER(KIND=JPIM) :: MU_SAVTEND_S,MV_SAVTEND_S,MT_SAVTEND_S,MSAT_SAVTEND_S
INTEGER(KIND=JPIM) :: MSAVTEND_S

!$OMP THREADPRIVATE(lslphy,msat_savtend,msat_savtend_s,msavtend_s,mt_savtend,mt_savtend_s)
!$OMP THREADPRIVATE(mu_savtend,mu_savtend_s,mv_savtend,mv_savtend_s,nvtend)
!$OMP THREADPRIVATE(savtend)
END MODULE YOMSLPHY
