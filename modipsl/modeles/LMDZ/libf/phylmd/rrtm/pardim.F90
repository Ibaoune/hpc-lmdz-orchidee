MODULE PARDIM

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!     JPMXLE : MAXIMUM NUMBER OF LEVELS
!     JPMXGL : MAXIMUM NUMBER OF GAUSSIAN LATITUDES 

!     JPFPPYX : Maximum number of fields in catalogue cilipdy

!     JPNULNAM : Fortran unit number used for namelist reading

INTEGER(KIND=JPIM), PARAMETER :: JPMXLE=200
INTEGER(KIND=JPIM), PARAMETER :: JPMXGL=2056
INTEGER(KIND=JPIM), PARAMETER :: JPFPPYX=17
INTEGER(KIND=JPIM), PARAMETER :: JPNULNAM=4

!     ------------------------------------------------------------------
END MODULE PARDIM
