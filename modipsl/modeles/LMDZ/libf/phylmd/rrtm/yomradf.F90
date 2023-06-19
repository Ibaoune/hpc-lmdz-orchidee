MODULE YOMRADF

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

! EMTD    - longwave emissivity
! TRSW    - shortwave absorptivity
! EMTC    - clear-sky longwave emissivity
! TRSC    - clear-sky shortwave absorptivity

! SRSWD   - downward SW radiation at the surface 
! SRLWD   - downward LW radiation at the surface
! SRSWDCS - clear-sky downward SW radiation at the surface
! SRLWDCS - clear-sky downward LW radiation at the surface 
! SRSWDV  - downward SW visible radiation at the surface
! SRSWDUV - downward SW ultraviolet/visible radiation at the surface
! SRSWPAR - downward SW PAR radiation at the surface
! SRSWUVB - downward UV-B radiation at the surface
! SRSWPARC- downward clear-sky SW PAR radiation at the surface
! SRSWTINC- TOA incident solar radiation 

REAL(KIND=JPRB),ALLOCATABLE :: EMTD(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: TRSW(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: EMTC(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: TRSC(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: EMTU(:,:,:)

REAL(KIND=JPRB),ALLOCATABLE :: SRSWD(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRLWD(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWDCS(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRLWDCS(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWDV(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWDUV(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: EDRO(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWPAR(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWUVB(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWPARC(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWTINC(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RMOON(:,:)


!$OMP THREADPRIVATE(edro,emtc,emtd,emtu,rmoon,srlwd,srlwdcs,srswd,srswdcs,srswduv)
!$OMP THREADPRIVATE(srswdv,srswpar,srswparc,srswtinc,srswuvb,trsc,trsw)
END MODULE YOMRADF
