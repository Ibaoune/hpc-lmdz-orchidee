#ifdef RS6K
@PROCESS HOT NOSTRICT
#endif
SUBROUTINE SRTM_REFTRA &
 & ( KLEV  , KMODTS, &
 &   LDRTCHK, &
 &   PGG   , PRMUZ, PTAU , PW, &
 &   PREF  , PREFD, PTRA , PTRAD &
 & )  
  
!**** *SRTM_REFTRA* - REFLECTIVITY AND TRANSMISSIVITY

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLEAR OR 
!     CLOUDY LAYER USING A CHOICE OF VARIOUS APPROXIMATIONS.

!**   INTERFACE.
!     ----------
!          *SRTM_REFTRA* IS CALLED BY *SRTM_SPCVRT*

!        EXPLICIT ARGUMENTS :
!        --------------------
! INPUTS
! ------ 
!      KMODTS  = 1 EDDINGTON (JOSEPH ET AL., 1976)
!              = 2 PIFM (ZDUNKOWSKI ET AL., 1980)
!              = 3 DISCRETE ORDINATES (LIOU, 1973)
!      LDRTCHK = .T. IF CLOUDY
!              = .F. IF CLEAR-SKY
!      PGG     = ASSYMETRY FACTOR
!      PRMUZ   = COSINE SOLAR ZENITH ANGLE
!      PTAU    = OPTICAL THICKNESS
!      PW      = SINGLE SCATTERING ALBEDO

! OUTPUTS
! -------
!      PREF    : COLLIMATED BEAM REFLECTIVITY
!      PREFD   : DIFFUSE BEAM REFLECTIVITY 
!      PTRA    : COLLIMATED BEAM TRANSMISSIVITY
!      PTRAD   : DIFFUSE BEAM TRANSMISSIVITY

!     METHOD.
!     -------
!          STANDARD DELTA-EDDINGTON, P.I.F.M., OR D.O.M. LAYER CALCULATIONS.

!     EXTERNALS.
!     ----------
!          NONE

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 03-02-27
!        M.Hamrud   01-Oct-2003      CY28 Cleaning
!        Mike Iacono, AER, Mar 2004: bug fix 

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARSRTM , ONLY : JPLAY
!USE YOESW   , ONLY : NDBUG
USE YOERDU  , ONLY : REPLOG

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KMODTS 
LOGICAL           ,INTENT(IN)    :: LDRTCHK(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGG(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMUZ 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PW(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PREF(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PREFD(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PTRA(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PTRAD(JPLAY) 
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK

REAL(KIND=JPRB) :: ZA, ZA1, ZA2
REAL(KIND=JPRB) :: ZBETA, ZDEND, ZDENR, ZDENT
REAL(KIND=JPRB) :: ZE1, ZE2, ZEM1, ZEM2, ZEMM, ZEP1, ZEP2
REAL(KIND=JPRB) :: ZG, ZG3, ZGAMMA1, ZGAMMA2, ZGAMMA3, ZGAMMA4, ZGT
REAL(KIND=JPRB) :: ZR1, ZR2, ZR3, ZR4, ZR5, ZRK, ZRK2, ZRKG, ZRM1, ZRP, ZRP1, ZRPP
REAL(KIND=JPRB) :: ZSR3, ZT1, ZT2, ZT3, ZT4, ZT5, ZTO1
REAL(KIND=JPRB) :: ZW, ZWCRIT, ZWO
REAL(KIND=JPRB) :: ZHOOK_HANDLE,EXP500,ZTEMP

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SRTM_REFTRA',0,ZHOOK_HANDLE)
EXP500=EXP(500.0_JPRB)
ZSR3=SQRT(3._JPRB)
ZWCRIT=0.9995_JPRB
KMODTS=2

!NDBUG=3

DO JK=1,KLEV
!  if (NDBUG < 2) then
!    print 9000,JK,LRTCHK(JK),PTAU(JK),PW(JK),PGG(JK),PRMUZ
  9000 format(1x,'SRTM_REFTRA:inputs:',I3,L8,4E13.6)
!  end if
  IF (.NOT.LDRTCHK(JK)) THEN
    PREF(JK) =0.0_JPRB
    PTRA(JK) =1.0_JPRB
    PREFD(JK)=0.0_JPRB
    PTRAD(JK)=1.0_JPRB
!    if (NDBUG < 2) then
!      print 9001,JL,JK,PREF(JK),PTRA(JK),PREFD(JK),PTRAD(JK)
    9001  format(1x,'SRTM_REFTRA:not.LRTCKH:',2I3,4F10.6)
!    end if
  ELSE
    ZTO1=PTAU(JK)
    ZW  =PW(JK)
    ZG  =PGG(JK)  

!-- GENERAL TWO-STREAM EXPRESSIONS

    ZG3= 3._JPRB * ZG
    IF (KMODTS == 1) THEN
      ZGAMMA1= (7._JPRB - ZW * (4._JPRB + ZG3)) * 0.25_JPRB
      ZGAMMA2=-(1._JPRB - ZW * (4._JPRB - ZG3)) * 0.25_JPRB
      ZGAMMA3= (2._JPRB - ZG3 * PRMUZ ) * 0.25_JPRB
    ELSEIF (KMODTS == 2) THEN  
      ZGAMMA1= (8._JPRB - ZW * (5._JPRB + ZG3)) * 0.25_JPRB
      ZGAMMA2=  3._JPRB *(ZW * (1._JPRB - ZG )) * 0.25_JPRB
      ZGAMMA3= (2._JPRB - ZG3 * PRMUZ ) * 0.25_JPRB
    ELSEIF (KMODTS == 3) THEN  
      ZGAMMA1= ZSR3 * (2._JPRB - ZW * (1._JPRB + ZG)) * 0.5_JPRB
      ZGAMMA2= ZSR3 * ZW * (1._JPRB - ZG ) * 0.5_JPRB
      ZGAMMA3= (1._JPRB - ZSR3 * ZG * PRMUZ ) * 0.5_JPRB
    ENDIF
    ZGAMMA4= 1._JPRB - ZGAMMA3
    
!-- RECOMPUTE ORIGINAL S.S.A. TO TEST FOR CONSERVATIVE SOLUTION
    ZWO= ZW / (1._JPRB - (1._JPRB - ZW) * (ZG / (1._JPRB - ZG))**2)
!   ZTEMP=(1._JPRB - ZG)**2
!   ZWO= ZW*ZTEMP/ (ZTEMP - (1._JPRB - ZW)*(ZG **2))
    
    IF (ZWO >= ZWCRIT) THEN
!!-- conservative scattering

      ZA  = ZGAMMA1 * PRMUZ 
      ZA1 = ZA - ZGAMMA3
      ZGT = ZGAMMA1 * ZTO1
        
!-- Homogeneous reflectance and transmittance

! collimated beam

      ZE1 = MIN ( ZTO1 / PRMUZ , 500._JPRB)
      ZE2 = EXP ( - ZE1 )
      ZTEMP=1.0_JPRB/(1._JPRB + ZGT)
      PREF(JK) = (ZGT - ZA1 * (1._JPRB - ZE2)) *ZTEMP
      PTRA(JK) = 1._JPRB - PREF(JK)

! isotropic incidence

      PREFD(JK) = ZGT *ZTEMP
      PTRAD(JK) = 1._JPRB - PREFD(JK)        

!      if (NDBUG < 2) then
!        print 9002,JL,JK,PREF(JK),PTRA(JK),PREFD(JK),PTRAD(JK)
      9002  format(1x,'SRTM_REFTRA: consrv: LDRTCHK:',2I3,4F10.6)
!      end if

    ELSE

!-- non-conservative scattering

      ZA1 = ZGAMMA1 * ZGAMMA4 + ZGAMMA2 * ZGAMMA3
      ZA2 = ZGAMMA1 * ZGAMMA3 + ZGAMMA2 * ZGAMMA4
!      ZRK = SQRT ( ZGAMMA1**2 - ZGAMMA2**2)
      ZRK = SQRT ( MAX ( REPLOG, ZGAMMA1**2 - ZGAMMA2**2) )
      ZRP = ZRK * PRMUZ               
      ZRP1 = 1._JPRB + ZRP
      ZRM1 = 1._JPRB - ZRP
      ZRK2 = 2._JPRB * ZRK
      ZRPP = 1._JPRB - ZRP*ZRP
      ZRKG = ZRK + ZGAMMA1
      ZR1  = ZRM1 * (ZA2 + ZRK * ZGAMMA3)
      ZR2  = ZRP1 * (ZA2 - ZRK * ZGAMMA3)
      ZR3  = ZRK2 * (ZGAMMA3 - ZA2 * PRMUZ )
      ZR4  = ZRPP * ZRKG
      ZR5  = ZRPP * (ZRK - ZGAMMA1)
      ZT1  = ZRP1 * (ZA1 + ZRK * ZGAMMA4)
      ZT2  = ZRM1 * (ZA1 - ZRK * ZGAMMA4)
      ZT3  = ZRK2 * (ZGAMMA4 + ZA1 * PRMUZ )
      ZT4  = ZR4
      ZT5  = ZR5
      ZBETA = - ZR5 / ZR4
        
!-- Homogeneous reflectance and transmittance

      IF(ZRK * ZTO1 > 500._JPRB)THEN
        ZEP1=EXP500
      ELSE
        ZEP1=EXP(ZRK * ZTO1)
      ENDIF
      ZEM1=1.0_JPRB/ZEP1
      IF(ZTO1 > 500._JPRB*PRMUZ)THEN
        ZEP2=EXP500
      ELSE
        ZEP2=EXP(ZTO1 / PRMUZ)
      ENDIF
      ZEM2=1.0_JPRB/ZEP2

!     ZE1 = MIN ( ZRK * ZTO1, 500._JPRB)
!     ZE2 = MIN ( ZTO1 / PRMUZ , 500._JPRB)

!     ZEP1 = EXP( ZE1 )
!      ZEM1 = EXP(-ZE1 )
!     ZEM1=1.0_JPRB/ZEP1

!     ZEP2 = EXP( ZE2 )
!      ZEM2 = EXP(-ZE2 )
!     ZEM2=1.0_JPRB/ZEP2

! collimated beam

      ZDENR = ZR4*ZEP1 + ZR5*ZEM1
!- bug noticed by Mike Iacono
!      PREF(JK) = ZWO * (ZR1*ZEP1 - ZR2*ZEM1 - ZR3*ZEM2) / ZDENR
      PREF(JK) = ZW  * (ZR1*ZEP1 - ZR2*ZEM1 - ZR3*ZEM2) / ZDENR

      ZDENT = ZT4*ZEP1 + ZT5*ZEM1
!- bug noticed by Mike Iacono
!      PTRA(JK) = ZEM2 * (1._JPRB - ZWO * (ZT1*ZEP1 - ZT2*ZEM1 - ZT3*ZEP2) / ZDENT)
      PTRA(JK) = ZEM2 * (1._JPRB - ZW  * (ZT1*ZEP1 - ZT2*ZEM1 - ZT3*ZEP2) / ZDENT)

! diffuse beam

      ZEMM = ZEM1*ZEM1
      ZDEND = 1._JPRB / ( (1._JPRB - ZBETA*ZEMM ) * ZRKG)
      PREFD(JK) =  ZGAMMA2 * (1._JPRB - ZEMM) * ZDEND
      PTRAD(JK) =  ZRK2*ZEM1*ZDEND

!      if (NDBUG < 2) then        
!        print 9003,JL,JK,PREF(JK),PTRA(JK),PREFD(JK),PTRAD(JK)
      9003  format(1x,'SRTM_REFTRA: OMG<1:  LDRTCHK:',2I3,4F10.6)
!      end if
    ENDIF

  ENDIF         

ENDDO    

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_REFTRA',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_REFTRA     
