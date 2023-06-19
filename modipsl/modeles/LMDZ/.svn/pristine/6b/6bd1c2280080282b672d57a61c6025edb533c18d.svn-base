SUBROUTINE RRTM_SETCOEF_140GP (KLEV,P_COLDRY,P_WKL,&
 & P_FAC00,P_FAC01,P_FAC10,P_FAC11,P_FORFAC,K_JP,K_JT,K_JT1,&
 & P_COLH2O,P_COLCO2,P_COLO3,P_COLN2O,P_COLCH4,P_COLO2,P_CO2MULT,&
 & K_LAYTROP,K_LAYSWTCH,K_LAYLOW,PAVEL,P_TAVEL,P_SELFFAC,P_SELFFRAC,K_INDSELF)  

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714

!     Purpose:  For a given atmosphere, calculate the indices and
!     fractions related to the pressure and temperature interpolations.
!     Also calculate the values of the integrated Planck functions 
!     for each band at the level and layer temperatures.

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARRRTM  , ONLY : JPLAY     ,JPINPX
USE YOERRTRF , ONLY :       PREFLOG   ,TREF

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLDRY(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_WKL(JPINPX,JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FAC00(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FAC01(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FAC10(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FAC11(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FORFAC(JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_JP(JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_JT(JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_JT1(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLH2O(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLCO2(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLO3(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLN2O(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLCH4(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLO2(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_CO2MULT(JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_LAYTROP 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_LAYSWTCH 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_LAYLOW 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAVEL(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAVEL(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SELFFAC(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SELFFRAC(JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_INDSELF(JPLAY) 
!- from INTFAC      
!- from INTIND
!- from PROFDATA             
!- from PROFILE             
!- from SELF             
INTEGER(KIND=JPIM) :: JP1, I_LAY

REAL(KIND=JPRB) :: Z_CO2REG, Z_COMPFP, Z_FACTOR, Z_FP, Z_FT, Z_FT1, Z_PLOG, Z_SCALEFAC, Z_STPFAC, Z_WATER
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!#include "yoeratm.h"    

IF (LHOOK) CALL DR_HOOK('RRTM_SETCOEF_140GP',0,ZHOOK_HANDLE)
Z_STPFAC = 296._JPRB/1013._JPRB

K_LAYTROP  = 0
K_LAYSWTCH = 0
K_LAYLOW   = 0
DO I_LAY = 1, KLEV
!        Find the two reference pressures on either side of the
!        layer pressure.  Store them in JP and JP1.  Store in FP the
!        fraction of the difference (in ln(pressure)) between these
!        two values that the layer pressure lies.
  Z_PLOG = LOG(PAVEL(I_LAY))
  K_JP(I_LAY) = INT(36._JPRB - 5*(Z_PLOG+0.04_JPRB))
  IF (K_JP(I_LAY)  <  1) THEN
    K_JP(I_LAY) = 1
  ELSEIF (K_JP(I_LAY)  >  58) THEN
    K_JP(I_LAY) = 58
  ENDIF
  JP1 = K_JP(I_LAY) + 1
  Z_FP = 5._JPRB * (PREFLOG(K_JP(I_LAY)) - Z_PLOG)

!        Determine, for each reference pressure (JP and JP1), which
!        reference temperature (these are different for each  
!        reference pressure) is nearest the layer temperature but does
!        not exceed it.  Store these indices in JT and JT1, resp.
!        Store in FT (resp. FT1) the fraction of the way between JT
!        (JT1) and the next highest reference temperature that the 
!        layer temperature falls.

  K_JT(I_LAY) = INT(3._JPRB + (P_TAVEL(I_LAY)-TREF(K_JP(I_LAY)))/15._JPRB)
  IF (K_JT(I_LAY)  <  1) THEN
    K_JT(I_LAY) = 1
  ELSEIF (K_JT(I_LAY)  >  4) THEN
    K_JT(I_LAY) = 4
  ENDIF
  Z_FT = ((P_TAVEL(I_LAY)-TREF(K_JP(I_LAY)))/15._JPRB) - REAL(K_JT(I_LAY)-3)
  K_JT1(I_LAY) = INT(3._JPRB + (P_TAVEL(I_LAY)-TREF(JP1))/15._JPRB)
  IF (K_JT1(I_LAY)  <  1) THEN
    K_JT1(I_LAY) = 1
  ELSEIF (K_JT1(I_LAY)  >  4) THEN
    K_JT1(I_LAY) = 4
  ENDIF
  Z_FT1 = ((P_TAVEL(I_LAY)-TREF(JP1))/15._JPRB) - REAL(K_JT1(I_LAY)-3)

  Z_WATER = P_WKL(1,I_LAY)/P_COLDRY(I_LAY)
  Z_SCALEFAC = PAVEL(I_LAY) * Z_STPFAC / P_TAVEL(I_LAY)

!        If the pressure is less than ~100mb, perform a different
!        set of species interpolations.
!         IF (PLOG .LE. 4.56) GO TO 5300
!--------------------------------------         
  IF (Z_PLOG  >  4.56_JPRB) THEN
    K_LAYTROP =  K_LAYTROP + 1
!        For one band, the "switch" occurs at ~300 mb. 
    IF (Z_PLOG  >=  5.76_JPRB) K_LAYSWTCH = K_LAYSWTCH + 1
    IF (Z_PLOG  >=  6.62_JPRB) K_LAYLOW = K_LAYLOW + 1

    P_FORFAC(I_LAY) = Z_SCALEFAC / (1.0_JPRB+Z_WATER)

!        Set up factors needed to separately include the water vapor
!        self-continuum in the calculation of absorption coefficient.
!C           SELFFAC(LAY) = WATER * SCALEFAC / (1.+WATER)
    P_SELFFAC(I_LAY) = Z_WATER * P_FORFAC(I_LAY)
    Z_FACTOR = (P_TAVEL(I_LAY)-188.0_JPRB)/7.2_JPRB
    K_INDSELF(I_LAY) = MIN(9, MAX(1, INT(Z_FACTOR)-7))
    P_SELFFRAC(I_LAY) = Z_FACTOR - REAL(K_INDSELF(I_LAY) + 7)

!        Calculate needed column amounts.
    P_COLH2O(I_LAY) = 1.E-20_JPRB * P_WKL(1,I_LAY)
    P_COLCO2(I_LAY) = 1.E-20_JPRB * P_WKL(2,I_LAY)
    P_COLO3(I_LAY)  = 1.E-20_JPRB * P_WKL(3,I_LAY)
    P_COLN2O(I_LAY) = 1.E-20_JPRB * P_WKL(4,I_LAY)
    P_COLCH4(I_LAY) = 1.E-20_JPRB * P_WKL(6,I_LAY)
    P_COLO2(I_LAY)  = 1.E-20_JPRB * P_WKL(7,I_LAY)
    IF (P_COLCO2(I_LAY)  ==  0.0_JPRB) P_COLCO2(I_LAY) = 1.E-32_JPRB * P_COLDRY(I_LAY)
    IF (P_COLN2O(I_LAY)  ==  0.0_JPRB) P_COLN2O(I_LAY) = 1.E-32_JPRB * P_COLDRY(I_LAY)
    IF (P_COLCH4(I_LAY)  ==  0.0_JPRB) P_COLCH4(I_LAY) = 1.E-32_JPRB * P_COLDRY(I_LAY)
!        Using E = 1334.2 cm-1.
    Z_CO2REG = 3.55E-24_JPRB * P_COLDRY(I_LAY)
    P_CO2MULT(I_LAY)= (P_COLCO2(I_LAY) - Z_CO2REG) *&
     & 272.63_JPRB*EXP(-1919.4_JPRB/P_TAVEL(I_LAY))/(8.7604E-4_JPRB*P_TAVEL(I_LAY))  
!         GO TO 5400
!------------------
  ELSE
!        Above LAYTROP.
! 5300    CONTINUE

!        Calculate needed column amounts.
    P_FORFAC(I_LAY) = Z_SCALEFAC / (1.0_JPRB+Z_WATER)

    P_COLH2O(I_LAY) = 1.E-20_JPRB * P_WKL(1,I_LAY)
    P_COLCO2(I_LAY) = 1.E-20_JPRB * P_WKL(2,I_LAY)
    P_COLO3(I_LAY)  = 1.E-20_JPRB * P_WKL(3,I_LAY)
    P_COLN2O(I_LAY) = 1.E-20_JPRB * P_WKL(4,I_LAY)
    P_COLCH4(I_LAY) = 1.E-20_JPRB * P_WKL(6,I_LAY)
    P_COLO2(I_LAY)  = 1.E-20_JPRB * P_WKL(7,I_LAY)
    IF (P_COLCO2(I_LAY)  ==  0.0_JPRB) P_COLCO2(I_LAY) = 1.E-32_JPRB * P_COLDRY(I_LAY)
    IF (P_COLN2O(I_LAY)  ==  0.0_JPRB) P_COLN2O(I_LAY) = 1.E-32_JPRB * P_COLDRY(I_LAY)
    IF (P_COLCH4(I_LAY)  ==  0.0_JPRB) P_COLCH4(I_LAY) = 1.E-32_JPRB * P_COLDRY(I_LAY)
    Z_CO2REG = 3.55E-24_JPRB * P_COLDRY(I_LAY)
    P_CO2MULT(I_LAY)= (P_COLCO2(I_LAY) - Z_CO2REG) *&
     & 272.63_JPRB*EXP(-1919.4_JPRB/P_TAVEL(I_LAY))/(8.7604E-4_JPRB*P_TAVEL(I_LAY))  
!----------------     
  ENDIF
! 5400    CONTINUE

!        We have now isolated the layer ln pressure and temperature,
!        between two reference pressures and two reference temperatures 
!        (for each reference pressure).  We multiply the pressure 
!        fraction FP with the appropriate temperature fractions to get 
!        the factors that will be needed for the interpolation that yields
!        the optical depths (performed in routines TAUGBn for band n).

  Z_COMPFP = 1.0_JPRB - Z_FP
  P_FAC10(I_LAY) = Z_COMPFP * Z_FT
  P_FAC00(I_LAY) = Z_COMPFP * (1.0_JPRB - Z_FT)
  P_FAC11(I_LAY) = Z_FP * Z_FT1
  P_FAC01(I_LAY) = Z_FP * (1.0_JPRB - Z_FT1)

ENDDO

! MT 981104 
!-- Set LAYLOW for profiles with surface pressure less than 750 hPa. 
IF (K_LAYLOW == 0) K_LAYLOW=1

IF (LHOOK) CALL DR_HOOK('RRTM_SETCOEF_140GP',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_SETCOEF_140GP
