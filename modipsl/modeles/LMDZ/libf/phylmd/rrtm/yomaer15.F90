MODULE YOMAER15

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOMAER15* - RADIATIVE CHARACTERISTICS OF THE AEROSOLS
!*                          FROZEN VERSION (CYCLE 15) OF YOEAER
!     -----------------------------------------------------------------

!        * E.C.M.W.F. PHYSICS PACKAGE *

!        96-11: Ph. Dandin. Meteo-France
!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14
!        04-06: F. Bouyssel. Initialisation done in suaer15
!        05-09: A. Alias     Sulfate aerosol (Hu Ron Ming)

!  NAME      TYPE     PURPOSE
!  ----   :  ----   : -------
!  TAUA15 :  REAL     S.W. NORMALIZED OPTICAL THICKNESS AT 0.55 MICRON
!  RPIZA15:  REAL     S.W. SINGLE SCATTERING ALBEDO
!  RCGA15 :  REAL     S.W. ASSYMETRY FACTOR
!  RAER15 :  REAL     L.W. ABSORPTION COEFFICIENTS
!     -----------------------------------------------------------------

REAL(KIND=JPRB) :: TAUA15 (2,6)
REAL(KIND=JPRB) :: RPIZA15(2,6)
REAL(KIND=JPRB) :: RCGA15 (2,6)
REAL(KIND=JPRB) :: RAER15 (5,6)

!$OMP THREADPRIVATE(raer15,rcga15,rpiza15,taua15)
END MODULE YOMAER15
