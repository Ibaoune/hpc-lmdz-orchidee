MODULE YOMLEG

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    Description of Legendre polynomials

! RW          : weights of the gaussian quadrature
! RMU         : mu              sin(theta)
! R1MU2       : 1.-MU*MU        cos(theta)**2
! R1MUI       : 1./R1MU2      1/cos(theta)**2
! R1MUA       : 1./R1MU2/RA   1/(a*cos(theta)**2)
! RSQM2       : SQRT(R1MU2)     cos(theta)
! R1QM2       : 1./SQRT(R1MU2) 1/cos(theta)
! RACTHE      : 1./SQRT(R1MU2) 1/(a*cos(theta))
! RLATIG      : arcsin(mu)      theta  GLOBAL VIEW
! RLATI       : arcsin(mu)      theta
! RIPI0       : bi-cubic interpolation coefficients
! RIPI1       : bi-cubic interpolation coefficients
! RIPI2       : bi-cubic interpolation coefficients

REAL(KIND=JPRB),ALLOCATABLE:: RW(:)
REAL(KIND=JPRB),ALLOCATABLE:: RMU(:)
REAL(KIND=JPRB),ALLOCATABLE:: R1MU2(:)
REAL(KIND=JPRB),ALLOCATABLE:: R1MUI(:)
REAL(KIND=JPRB),ALLOCATABLE:: R1MUA(:)
REAL(KIND=JPRB),ALLOCATABLE:: RSQM2(:)
REAL(KIND=JPRB),ALLOCATABLE:: R1QM2(:)
REAL(KIND=JPRB),ALLOCATABLE:: RACTHE(:)
REAL(KIND=JPRB),ALLOCATABLE:: RLATIG(:)
REAL(KIND=JPRB),ALLOCATABLE:: RLATI(:)
REAL(KIND=JPRB),ALLOCATABLE:: RIPI0(:)
REAL(KIND=JPRB),ALLOCATABLE:: RIPI1(:)
REAL(KIND=JPRB),ALLOCATABLE:: RIPI2(:)

!$OMP THREADPRIVATE(r1mu2,r1mua,r1mui,r1qm2,racthe,ripi0,ripi1,ripi2,rlati,rlatig,rmu,rsqm2,rw)
END MODULE YOMLEG
