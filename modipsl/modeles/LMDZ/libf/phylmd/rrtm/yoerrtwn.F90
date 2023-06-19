MODULE YOERRTWN

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

INTEGER(KIND=JPIM) , DIMENSION(16) :: NG
INTEGER(KIND=JPIM) , DIMENSION(16) :: NSPA
INTEGER(KIND=JPIM) , DIMENSION(16) :: NSPB

REAL(KIND=JPRB) , DIMENSION(16) :: WAVENUM1
REAL(KIND=JPRB) , DIMENSION(16) :: WAVENUM2
REAL(KIND=JPRB) , DIMENSION(16) :: DELWAVE

REAL(KIND=JPRB) , DIMENSION(181,16) :: TOTPLNK
REAL(KIND=JPRB) , DIMENSION(181)    :: TOTPLK16

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----    : -------
!  NG     : INTEGER : Number of k-coefficients in spectral intervals
!  NSPA   : INTEGER :
!  NSPB   : INTEGER :
! WAVENUM1: REAL    : Lower wavenumber spectral limit
! WAVENUM2: REAL    : Higher wavenumber spectral limit
! DELWAVE : REAL    : Spectral interval width
! TOTPLNK : REAL    :
! TOTPLK16: REAL    :
!     -----------------------------------------------------------------
!$OMP THREADPRIVATE(delwave,ng,nspa,nspb,totplk16,totplnk,wavenum1,wavenum2)
END MODULE YOERRTWN
