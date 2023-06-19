MODULE PARFPOS

USE PARKIND1  ,ONLY : JPIM

USE PARDIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! === basic dimensions for Full POST-PROCESSING ===

!     JPOSDOM : Maximum number of horizontal (sub)domains
!     JPOSLEN : Maximum length of a (sub)domain name
!     JPOSLIS : Maximum number of groups of subdomains
!     JPOSDIR : Maximum length of the path (or prefix) for the output files
!     JPOSLE  : Maximum number of eta levels on the output subdomain
!     JPOSGL  : Maximum number of latitude rows of the output gaussian grid

!     JPOSSCVA: Maximum number of post-processable passive scalars
!     JPOSPHY : Maximum number of physical fields
!     JPOSSGP : Maximum number of surface gridpoint fields
!     JPOSCFU : Maximum number of cumulated fluxes
!     JPOSXFU : Maximum number of instantaneous fluxes
!     JPOS3P  : Maximum number of pp. pressure levels
!     JPOS3H  : Maximum number of pp. height (above orography) levels
!     JPOS3TH : Maximum number of pp. potential temperature levels
!     JPOS3PV : Maximum number of pp. potential vorticity levels
!     JPOS3S  : Maximum number of pp. eta levels
!     JPOSSCVA: Maximum number of passive scalars
!     JPOSVX2 : Maximum number of free gp/sp upper air fields or extra surface fields
!     JPOSFSU : Maximum number of free gp/sp surface fields 
!     JPOS3DF : Maximum number of 3D dynamic fields
!     JPOS2DF : Maximum number of 2D dynamic fields
!     JPOSDYN : Maximum number of dynamic fields 
!     JPOSGHG : Maximum number of post-processable greenhouse gases
!     JPOSTRAC : Maximum number of post-processable tracers (used for diagnostics only)
!     JPOSGRG : Maximum number of post-processable reactive gases
!     JPOSAERO: Maximum number of post-processable aerosols

INTEGER(KIND=JPIM), PARAMETER :: JPOSDOM=15
INTEGER(KIND=JPIM), PARAMETER :: JPOSLEN=10
INTEGER(KIND=JPIM), PARAMETER :: JPOSLIS=10
INTEGER(KIND=JPIM), PARAMETER :: JPOSDIR=180
INTEGER(KIND=JPIM), PARAMETER :: JPOSLE=JPMXLE
INTEGER(KIND=JPIM), PARAMETER :: JPOSGL=JPMXGL

INTEGER(KIND=JPIM), PARAMETER :: JPOSCFU=63
INTEGER(KIND=JPIM), PARAMETER :: JPOSXFU=63
INTEGER(KIND=JPIM), PARAMETER :: JPOS3P=75
INTEGER(KIND=JPIM), PARAMETER :: JPOS3H=127
INTEGER(KIND=JPIM), PARAMETER :: JPOS3TH=75
INTEGER(KIND=JPIM), PARAMETER :: JPOS3PV=75
INTEGER(KIND=JPIM), PARAMETER :: JPOS3S=JPOSLE
INTEGER(KIND=JPIM), PARAMETER :: JPOSSCVA=5
INTEGER(KIND=JPIM), PARAMETER :: JPOSGHG=3
INTEGER(KIND=JPIM), PARAMETER :: JPOSTRAC=2
INTEGER(KIND=JPIM), PARAMETER :: JPOSGRG=5
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERO=12
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERO2=2*JPOSAERO
INTEGER(KIND=JPIM), PARAMETER :: JPOSVX2=70
! 15 + 2 (Non-Hydro) + 4 (Modis-Albedo) + 3 (ocean T)
INTEGER(KIND=JPIM), PARAMETER :: JPOSFSU=24
INTEGER(KIND=JPIM), PARAMETER :: JPOSSGP=157+JPOSVX2+JPOSFSU
INTEGER(KIND=JPIM), PARAMETER :: JPOSPHY=JPOSSGP+JPOSCFU+JPOSXFU
INTEGER(KIND=JPIM), PARAMETER :: &
 & JPOS3DF=63+JPOSSCVA+JPOSVX2+JPOSAERO+JPOSGHG+JPOSGRG+JPOSTRAC
INTEGER(KIND=JPIM), PARAMETER :: JPOS2DF=63+JPOSFSU+JPOSAERO2
INTEGER(KIND=JPIM), PARAMETER :: JPOSDYN=JPOS3DF+JPOS2DF

!     ------------------------------------------------------------------
END MODULE PARFPOS
