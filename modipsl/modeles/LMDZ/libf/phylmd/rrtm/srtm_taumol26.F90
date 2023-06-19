SUBROUTINE SRTM_TAUMOL26 &
 & ( KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1     , P_ONEMINUS,&
 & P_COLH2O  , P_COLCO2 , P_COLMOL,&
 & K_LAYTROP , P_SELFFAC, P_SELFFRAC, K_INDSELF  , P_FORFAC, P_FORFRAC, K_INDFOR,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    &
 & )  

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARSRTM  , ONLY : JPLAY, JPG, NG26
USE YOESRTA26, ONLY : SFLUXREFC, RAYLC
IMPLICIT NONE

!-- Output
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)                  :: P_FAC00(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_FAC01(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_FAC10(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_FAC11(JPLAY) ! Argument NOT used
INTEGER(KIND=JPIM)               :: K_JP(JPLAY) ! Argument NOT used
INTEGER(KIND=JPIM)               :: K_JT(JPLAY) ! Argument NOT used
INTEGER(KIND=JPIM)               :: K_JT1(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_ONEMINUS ! Argument NOT used
REAL(KIND=JPRB)                  :: P_COLH2O(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_COLCO2(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(JPLAY) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP 
REAL(KIND=JPRB)                  :: P_SELFFAC(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_SELFFRAC(JPLAY) ! Argument NOT used
INTEGER(KIND=JPIM)               :: K_INDSELF(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_FORFAC(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_FORFRAC(JPLAY) ! Argument NOT used
INTEGER(KIND=JPIM)               :: K_INDFOR(JPLAY) ! Argument NOT used

REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SFLUXZEN(JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUG(JPLAY,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUR(JPLAY,JPG) 
!- from AER
!- from INTFAC      
!- from INTIND
!- from PRECISE             
!- from PROFDATA             
!- from SELF             
INTEGER(KIND=JPIM) :: IG, I_LAY, I_LAYSOLFR, I_NLAYERS

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL26',0,ZHOOK_HANDLE)
I_NLAYERS = KLEV

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  
I_LAYSOLFR = K_LAYTROP

DO I_LAY = 1, K_LAYTROP
!  DO IG = 1, NG(26)
  DO IG = 1 , NG26 
!    TAUG(LAY,IG) = COLMOL(LAY) * RAYLC(IG)
!    SSA(LAY,IG) = 1.0
    IF (I_LAY == I_LAYSOLFR) P_SFLUXZEN(IG) = SFLUXREFC(IG) 
    P_TAUG(I_LAY,IG) = 0.0_JPRB
    P_TAUR(I_LAY,IG) = P_COLMOL(I_LAY) * RAYLC(IG) 
  ENDDO
ENDDO

DO I_LAY = K_LAYTROP+1, I_NLAYERS
!  DO IG = 1, NG(26)
  DO IG = 1 , NG26
!    TAUG(LAY,IG) = COLMOL(LAY) * RAYLC(IG)
!    SSA(LAY,IG) = 1.0
    P_TAUG(I_LAY,IG) = 0.0_JPRB
    P_TAUR(I_LAY,IG) = P_COLMOL(I_LAY) * RAYLC(IG) 
  ENDDO
ENDDO

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL26',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_TAUMOL26

