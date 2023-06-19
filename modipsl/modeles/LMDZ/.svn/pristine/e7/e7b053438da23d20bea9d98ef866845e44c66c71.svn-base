MODULE YOMCVER

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

! -----------------------------------------------------------------------------

! * Variables related to vertical discretisation in finite elements:
! LVERTFE : .T./.F. Finite element/conventional vertical discretisation.
! NVSCH   : type of basis if the finite element vertical discretisation is used.
!           1: linear functions.
!           3: Hermite cubic functions.
! RINTE   : matricial operator for vertical integrations in the
!           finite element vertical discretisation.
! RDERI   : matricial operator for vertical derivatives in the
!           finite element vertical discretisation.
LOGICAL :: LVERTFE
INTEGER(KIND=JPIM) :: NVSCH
REAL(KIND=JPRB),ALLOCATABLE :: RINTE(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RDERI(:,:)

! -----------------------------------------------------------------------------

! * Variables related to use of spline cubic vertical interpolations
!   in the semi-Lagrangian scheme:
! LVSPLIP            : .T. if vertical spline cubic SL interpolations for O3.
! RVSPTRI,RVSPC      : are used to re-profile the field to be interpolated
!                      in routine VSPLTRANS.
! RFAA,RFBB,RFCC,RFDD: are used in the computation of the vertical weights.
! VRDETAR            : ratio (eta(lbar)-eta(lbar-1))/(eta(lbar-1)-eta(lbar-2)),
!                      is used in the interpolation routine LAITVSPCQM to
!                      ensure monotonicity and conservation properties.
LOGICAL :: LVSPLIP
LOGICAL :: LSVTSM                ! Stratospheric vertical trajectory smoothed
REAL(KIND=JPRB),ALLOCATABLE :: RVSPTRI(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RVSPC(:)
REAL(KIND=JPRB),ALLOCATABLE :: RFAA(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RFBB(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RFCC(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RFDD(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: VRDETAR(:)

! -----------------------------------------------------------------------------

!$OMP THREADPRIVATE(lsvtsm,lvertfe,lvsplip,nvsch)
!$OMP THREADPRIVATE(rderi,rfaa,rfbb,rfcc,rfdd,rinte,rvspc,rvsptri,vrdetar)
END MODULE YOMCVER
