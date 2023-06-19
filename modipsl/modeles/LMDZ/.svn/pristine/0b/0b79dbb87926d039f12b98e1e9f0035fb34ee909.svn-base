SUBROUTINE SRTM_SETCOEF &
 & ( KLEV    , KNMOL ,&
 &   PAVEL   , PTAVEL   , PZ      , PTZ     , PTBOUND  ,&
 &   PCOLDRY , PWKL     ,&
 &   KLAYTROP, KLAYSWTCH, KLAYLOW ,&
 &   PCO2MULT, PCOLCH4  , PCOLCO2 , PCOLH2O , PCOLMOL  , PCOLN2O  , PCOLO2 , PCOLO3 ,&
 &   PFORFAC , PFORFRAC , KINDFOR , PSELFFAC, PSELFFRAC, KINDSELF ,&
 &   PFAC00  , PFAC01   , PFAC10  , PFAC11  ,&
 &   KJP     , KJT      , KJT1    &
 & )  

!     J. Delamere, AER, Inc. (version 2.5, 02/04/01)

!     Modifications:
!     JJMorcrette 030224   rewritten / adapted to ECMWF F90 system
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     Purpose:  For a given atmosphere, calculate the indices and
!     fractions related to the pressure and temperature interpolations.

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARSRTM , ONLY : JPLAY
USE YOESRTWN, ONLY :  PREFLOG, TREF
!!  USE YOESWN  , ONLY : NDBUG
                   
IMPLICIT NONE

!-- Input arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM)               :: KNMOL ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAVEL(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAVEL(JPLAY) 
REAL(KIND=JPRB)                  :: PZ(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTZ(0:JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTBOUND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLDRY(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWKL(35,JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLAYTROP 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLAYSWTCH 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLAYLOW 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCO2MULT(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLCH4(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLCO2(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLH2O(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLMOL(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLN2O(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLO2(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLO3(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFORFAC(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFORFRAC(JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KINDFOR(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSELFFAC(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSELFFRAC(JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KINDSELF(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFAC00(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFAC01(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFAC10(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFAC11(JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KJP(JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KJT(JPLAY) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KJT1(JPLAY) 
!-- Output arguments

!-- local integers

INTEGER(KIND=JPIM) :: I_NLAYERS, INDBOUND, INDLEV0, JK
INTEGER(KIND=JPIM) :: JP1

!-- local reals

REAL(KIND=JPRB) :: Z_STPFAC, Z_TBNDFRAC, Z_T0FRAC, Z_PLOG, Z_FP, Z_FT, Z_FT1, Z_WATER, Z_SCALEFAC
REAL(KIND=JPRB) :: Z_FACTOR, Z_CO2REG, Z_COMPFP
REAL(KIND=JPRB) :: ZHOOK_HANDLE




IF (LHOOK) CALL DR_HOOK('SRTM_SETCOEF',0,ZHOOK_HANDLE)
I_NLAYERS = KLEV

Z_STPFAC = 296._JPRB/1013._JPRB

INDBOUND = PTBOUND - 159._JPRB
Z_TBNDFRAC = PTBOUND - INT(PTBOUND)
INDLEV0  = PTZ(0) - 159._JPRB
Z_T0FRAC   = PTZ(0) - INT(PTZ(0))

KLAYTROP  = 0
KLAYSWTCH = 0
KLAYLOW   = 0

!IF (NDBUG.LE.3) THEN
!  print *,'-------- Computed in SETCOEF --------'
!  print 8990
8990 format(18x,'  T     PFAC00,    01,    10,    11  PCO2MULT     MOL   &
 & CH4      CO2      H2O      N2O      O2      O3      SFAC  &
 & SFRAC    FFAC    FFRAC  ISLF IFOR')  
!END IF

DO JK = 1, I_NLAYERS
!        Find the two reference pressures on either side of the
!        layer pressure.  Store them in JP and JP1.  Store in FP the
!        fraction of the difference (in ln(pressure)) between these
!        two values that the layer pressure lies.

  Z_PLOG = LOG(PAVEL(JK))
  KJP(JK) = INT(36. - 5*(Z_PLOG+0.04))
  IF (KJP(JK) < 1) THEN
    KJP(JK) = 1
  ELSEIF (KJP(JK) > 58) THEN
    KJP(JK) = 58
  ENDIF
  JP1 = KJP(JK) + 1
  Z_FP = 5. * (PREFLOG(KJP(JK)) - Z_PLOG)

!        Determine, for each reference pressure (JP and JP1), which
!        reference temperature (these are different for each  
!        reference pressure) is nearest the layer temperature but does
!        not exceed it.  Store these indices in JT and JT1, resp.
!        Store in FT (resp. FT1) the fraction of the way between JT
!        (JT1) and the next highest reference temperature that the 
!        layer temperature falls.

  KJT(JK) = INT(3. + (PTAVEL(JK)-TREF(KJP(JK)))/15.)
  IF (KJT(JK) < 1) THEN
    KJT(JK) = 1
  ELSEIF (KJT(JK) > 4) THEN
    KJT(JK) = 4
  ENDIF
  Z_FT = ((PTAVEL(JK)-TREF(KJP(JK)))/15.) - REAL(KJT(JK)-3)
  KJT1(JK) = INT(3. + (PTAVEL(JK)-TREF(JP1))/15.)
  IF (KJT1(JK) < 1) THEN
    KJT1(JK) = 1
  ELSEIF (KJT1(JK) > 4) THEN
    KJT1(JK) = 4
  ENDIF
  Z_FT1 = ((PTAVEL(JK)-TREF(JP1))/15.) - REAL(KJT1(JK)-3)

  Z_WATER = PWKL(1,JK)/PCOLDRY(JK)
  Z_SCALEFAC = PAVEL(JK) * Z_STPFAC / PTAVEL(JK)

!        If the pressure is less than ~100mb, perform a different
!        set of species interpolations.

  IF (Z_PLOG <= 4.56) GO TO 5300
  KLAYTROP =  KLAYTROP + 1
  IF (Z_PLOG >= 6.62) KLAYLOW = KLAYLOW + 1

!        Set up factors needed to separately include the water vapor
!        foreign-continuum in the calculation of absorption coefficient.

  PFORFAC(JK) = Z_SCALEFAC / (1.+Z_WATER)
  Z_FACTOR = (332.0-PTAVEL(JK))/36.0
  KINDFOR(JK) = MIN(2, MAX(1, INT(Z_FACTOR)))
  PFORFRAC(JK) = Z_FACTOR - REAL(KINDFOR(JK))

!        Set up factors needed to separately include the water vapor
!        self-continuum in the calculation of absorption coefficient.

  PSELFFAC(JK) = Z_WATER * PFORFAC(JK)
  Z_FACTOR = (PTAVEL(JK)-188.0)/7.2
  KINDSELF(JK) = MIN(9, MAX(1, INT(Z_FACTOR)-7))
  PSELFFRAC(JK) = Z_FACTOR - REAL(KINDSELF(JK) + 7)

!        Calculate needed column amounts.

  PCOLH2O(JK) = 1.E-20 * PWKL(1,JK)
  PCOLCO2(JK) = 1.E-20 * PWKL(2,JK)
  PCOLO3(JK) = 1.E-20 * PWKL(3,JK)
!         COLO3(LAY) = 0.
!         COLO3(LAY) = colo3(lay)/1.16
  PCOLN2O(JK) = 1.E-20 * PWKL(4,JK)
  PCOLCH4(JK) = 1.E-20 * PWKL(6,JK)
  PCOLO2(JK) = 1.E-20 * PWKL(7,JK)
  PCOLMOL(JK) = 1.E-20 * PCOLDRY(JK) + PCOLH2O(JK)
!         colco2(lay) = 0.
!         colo3(lay) = 0.
!         coln2o(lay) = 0.
!         colch4(lay) = 0.
!         colo2(lay) = 0.
!         colmol(lay) = 0.
  IF (PCOLCO2(JK) == 0.) PCOLCO2(JK) = 1.E-32 * PCOLDRY(JK)
  IF (PCOLN2O(JK) == 0.) PCOLN2O(JK) = 1.E-32 * PCOLDRY(JK)
  IF (PCOLCH4(JK) == 0.) PCOLCH4(JK) = 1.E-32 * PCOLDRY(JK)
  IF (PCOLO2(JK) == 0.) PCOLO2(JK) = 1.E-32 * PCOLDRY(JK)
!        Using E = 1334.2 cm-1.
  Z_CO2REG = 3.55E-24 * PCOLDRY(JK)
  PCO2MULT(JK)= (PCOLCO2(JK) - Z_CO2REG) * &
   & 272.63*EXP(-1919.4/PTAVEL(JK))/(8.7604E-4*PTAVEL(JK))  
  GO TO 5400

!        Above LAYTROP.
  5300    CONTINUE

!        Set up factors needed to separately include the water vapor
!        foreign-continuum in the calculation of absorption coefficient.

  PFORFAC(JK) = Z_SCALEFAC / (1.+Z_WATER)
  Z_FACTOR = (PTAVEL(JK)-188.0)/36.0
  KINDFOR(JK) = 3
  PFORFRAC(JK) = Z_FACTOR - 1.0

!        Calculate needed column amounts.

  PCOLH2O(JK) = 1.E-20 * PWKL(1,JK)
  PCOLCO2(JK) = 1.E-20 * PWKL(2,JK)
  PCOLO3(JK)  = 1.E-20 * PWKL(3,JK)
  PCOLN2O(JK) = 1.E-20 * PWKL(4,JK)
  PCOLCH4(JK) = 1.E-20 * PWKL(6,JK)
  PCOLO2(JK)  = 1.E-20 * PWKL(7,JK)
  PCOLMOL(JK) = 1.E-20 * PCOLDRY(JK) + PCOLH2O(JK)
  IF (PCOLCO2(JK) == 0.) PCOLCO2(JK) = 1.E-32 * PCOLDRY(JK)
  IF (PCOLN2O(JK) == 0.) PCOLN2O(JK) = 1.E-32 * PCOLDRY(JK)
  IF (PCOLCH4(JK) == 0.) PCOLCH4(JK) = 1.E-32 * PCOLDRY(JK)
  IF (PCOLO2(JK) == 0.) PCOLO2(JK)  = 1.E-32 * PCOLDRY(JK)
  Z_CO2REG = 3.55E-24 * PCOLDRY(JK)
  PCO2MULT(JK)= (PCOLCO2(JK) - Z_CO2REG) * &
   & 272.63*EXP(-1919.4/PTAVEL(JK))/(8.7604E-4*PTAVEL(JK))  

  PSELFFAC(JK) =0.0_JPRB
  PSELFFRAC(JK)=0.0_JPRB
  KINDSELF(JK) = 0

  5400    CONTINUE

!        We have now isolated the layer ln pressure and temperature,
!        between two reference pressures and two reference temperatures 
!        (for each reference pressure).  We multiply the pressure 
!        fraction FP with the appropriate temperature fractions to get 
!        the factors that will be needed for the interpolation that yields
!        the optical depths (performed in routines TAUGBn for band n).

  Z_COMPFP = 1. - Z_FP
  PFAC10(JK) = Z_COMPFP * Z_FT
  PFAC00(JK) = Z_COMPFP * (1. - Z_FT)
  PFAC11(JK) = Z_FP * Z_FT1
  PFAC01(JK) = Z_FP * (1. - Z_FT1)

!  IF (NDBUG.LE.3) THEN
!    print 9000,LAY,LAYTROP,JP(LAY),JT(LAY),JT1(LAY),TAVEL(LAY) &
!      &,FAC00(LAY),FAC01(LAY),FAC10(LAY),FAC11(LAY) &
!      &,CO2MULT(LAY),COLMOL(LAY),COLCH4(LAY),COLCO2(LAY),COLH2O(LAY) &
!      &,COLN2O(LAY),COLO2(LAY),COLO3(LAY),SELFFAC(LAY),SELFFRAC(LAY) &
!      &,FORFAC(LAY),FORFRAC(LAY),INDSELF(LAY),INDFOR(LAY)
  9000 format(1x,2I3,3I4,F6.1,4F7.2,12E9.2,2I5)
!  END IF

ENDDO

!----------------------------------------------------------------------- 
IF (LHOOK) CALL DR_HOOK('SRTM_SETCOEF',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_SETCOEF

