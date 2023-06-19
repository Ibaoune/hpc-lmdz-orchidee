MODULE YOE_TILE_PROP

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!       ----------------------------------------------------------------
!PSEUDO HISTORICAL TESSEL QUANTITIES (INSTANTANEOUS) FOR EACH TILE

!     M.Hamrud      E.C.M.W.F.    20 Apr 2005
!     Moved out these arrays from GPPBUF     
!       ----------------------------------------------------------------
REAL(KIND=JPRB),ALLOCATABLE :: RUSTRTI(:,:,:) ! E-W  SURFACE STRESS 
REAL(KIND=JPRB),ALLOCATABLE :: RVSTRTI(:,:,:) ! N-S  SURFACE STRESS 
REAL(KIND=JPRB),ALLOCATABLE :: RAHFSTI(:,:,:) ! SURFACE SENSIBLE HEAT FLUX
REAL(KIND=JPRB),ALLOCATABLE :: REVAPTI(:,:,:) ! EVAPORATION 
REAL(KIND=JPRB),ALLOCATABLE :: RTSKTI (:,:,:) ! SKIN TEMPERATURE 

!$OMP THREADPRIVATE(rahfsti,revapti,rtskti,rustrti,rvstrti)
END MODULE YOE_TILE_PROP
