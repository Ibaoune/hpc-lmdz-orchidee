MODULE YOMCAPE

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!*
!     ------------------------------------------------------------------

!     VARIABLES TO CONTROL CAPE COMPUTATION IN FULLPOS:

!        NCAPEITER : NUMBER OF ITERATIONS IN THE NEWTON LOOPS. Default
!                    is the same as NBITER in YOMPHY.

!        NETAPES : NUMBER OF INTERMEDIATE-LAYERS USED FOR CALCULATION OF
!                  VERTICAL ASCENT BETWEEN TWO MODEL PRESSURE LEVELS.
!                  Default value is 2.

!        GCAPEPSD : DEPTH OF LAYER ABOVE THE GROUND IN WHICH MOST 
!                   UNSTABLE PARCEL IS SEARCHED FOR (UNIT: Pa)
!                   USED JUST FOR TYPE 2 OF CAPE. Default value is 30000.
!                   (CAPE Pressure Search Depth)

!        GCAPERET: FRACTION OF THE CONDENSATES WHICH IS RETAINED, 
!                  I.E. WHICH DOES NOT PRECIPITATE.
!            IF GCAPERET=1. ==> REVERSIBLE MOIST ASCENT.
!                      IT IS ASSUMED THAT ALL THE PARCEL'S CONDENSED
!                      WATER IS RETAINED, THUS CLOUD CONDENSATES 
!                      REDUCE THE BUOYANCY.
!            IF GCAPERET=0. ==> "IRREVERSIBLE" (PSEUDO-ADIABATIC) MOIST ASCENT.
!                       CLOUD CONDENSATES PRECIPITATE  INSTANTANEOUSLY
!                       AND THUS DO NOT AFFECT THE BUOYANCY.
!            GCAPERET CAN BE USED WITH VALUES BETWEEN 0. AND 1..
!            Default value is 0.
!-------------------------------------------------

INTEGER(KIND=JPIM) :: NCAPEITER
INTEGER(KIND=JPIM) :: NETAPES

REAL(KIND=JPRB) :: GCAPERET
REAL(KIND=JPRB) :: GCAPEPSD

!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(gcapepsd,gcaperet,ncapeiter,netapes)
END MODULE YOMCAPE
