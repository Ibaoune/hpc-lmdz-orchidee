SUBROUTINE SRTM_SRTM_224GP_MCICA &
 & ( KIDIA , KFDIA  , KLON  , KLEV  , KSW , KCOLS , KCLDLY ,&
 &   PAER  , PALBD  , PALBP , PAPH  , PAP , &
 &   PTS   , PTH    , PT    ,&
 &   PQ    , PCCO2  , POZN  , PRMU0 ,&
 &   PFRCL , PTAUC  , PASYC , POMGC ,&
 &   PFSUX , PFSUC &
 & )  

!-- interface to RRTM_SW
!     JJMorcrette 030225
!     JJMorcrette 20050110  McICA version

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARSRTM  , ONLY : JPLAY
!MPL/IM 20160915 on prend GES de phylmd USE YOERDI   , ONLY : RCH4   , RN2O   
USE YOERAD   , ONLY : NAER
USE YOESRTAER, ONLY : RSRTAUA, RSRPIZA, RSRASYA
USE YOMPHY3  , ONLY : RII0
USE YOMCST   , ONLY : RI0 

IMPLICIT NONE

!-- Input arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW  
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOLS
INTEGER(KIND=JPIM),INTENT(IN)    :: KCLDLY(KCOLS) 

REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV)    ! top to bottom
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KLON,KSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,KSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCO2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KLON)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRCL(KLON,KCOLS,KLEV) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUC(KLON,KCOLS,KLEV) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PASYC(KLON,KCOLS,KLEV) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMGC(KLON,KCOLS,KLEV) ! bottom to top

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSUX(KLON,2,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSUC(KLON,2,KLEV+1) 

!-- Output arguments

!-----------------------------------------------------------------------

!-- dummy integers

INTEGER(KIND=JPIM) :: ICLDATM, INFLAG, ICEFLAG, I_LIQFLAG, ITMOL, I_NSTR

INTEGER(KIND=JPIM) :: IK, IMOL, J1, J2, JAE, JL, JK, JSW

!-- dummy reals

REAL(KIND=JPRB) :: ZPZ(0:JPLAY)   , ZTZ(0:JPLAY)   , ZPAVEL(JPLAY)  , ZTAVEL(JPLAY)
REAL(KIND=JPRB) :: ZCOLDRY(JPLAY) , ZCOLMOL(JPLAY) , ZWKL(35,JPLAY)
REAL(KIND=JPRB) :: ZCO2MULT(JPLAY), ZCOLCH4(JPLAY) , ZCOLCO2(JPLAY) , ZCOLH2O(JPLAY)
REAL(KIND=JPRB) :: ZCOLN2O(JPLAY) , ZCOLO2(JPLAY)  , ZCOLO3(JPLAY)
REAL(KIND=JPRB) :: ZFORFAC(JPLAY) , ZFORFRAC(JPLAY), ZSELFFAC(JPLAY), ZSELFFRAC(JPLAY)
REAL(KIND=JPRB) :: ZFAC00(JPLAY)  , ZFAC01(JPLAY)  , ZFAC10(JPLAY)  , ZFAC11(JPLAY)
REAL(KIND=JPRB) :: ZTBOUND        , ZONEMINUS    , ZRMU0 , ZADJI0
REAL(KIND=JPRB) :: ZALBD(KSW)    , ZALBP(KSW)    

REAL(KIND=JPRB) :: ZFRCL(KCOLS,JPLAY), ZTAUC(JPLAY,KCOLS), ZASYC(JPLAY,KCOLS), ZOMGC(JPLAY,KCOLS)
REAL(KIND=JPRB) :: ZTAUA(JPLAY,KSW), ZASYA(JPLAY,KSW), ZOMGA(JPLAY,KSW)

REAL(KIND=JPRB) :: ZBBCD(JPLAY+1), ZBBCU(JPLAY+1), ZBBFD(JPLAY+1), ZBBFU(JPLAY+1)
!REAL(KIND=JPRB) :: ZUVCD(JPLAY+1), ZUVCU(JPLAY+1), ZUVFD(JPLAY+1), ZUVFU(JPLAY+1)
!REAL(KIND=JPRB) :: ZVSCD(JPLAY+1), ZVSCU(JPLAY+1), ZVSFD(JPLAY+1), ZVSFU(JPLAY+1)
!REAL(KIND=JPRB) :: ZNICD(JPLAY+1), ZNICU(JPLAY+1), ZNIFD(JPLAY+1), ZNIFU(JPLAY+1)

INTEGER(KIND=JPIM) :: ILAYTROP, ILAYSWTCH, ILAYLOW
INTEGER(KIND=JPIM) :: INDFOR(JPLAY), INDSELF(JPLAY)
INTEGER(KIND=JPIM) :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

REAL(KIND=JPRB) :: ZAMD                  ! Effective molecular weight of dry air (g/mol)
REAL(KIND=JPRB) :: ZAMW                  ! Molecular weight of water vapor (g/mol)
REAL(KIND=JPRB) :: ZAMCO2                ! Molecular weight of carbon dioxide (g/mol)
REAL(KIND=JPRB) :: ZAMO                  ! Molecular weight of ozone (g/mol)
REAL(KIND=JPRB) :: ZAMCH4                ! Molecular weight of methane (g/mol)
REAL(KIND=JPRB) :: ZAMN2O                ! Molecular weight of nitrous oxide (g/mol)
REAL(KIND=JPRB) :: ZAMC11                ! Molecular weight of CFC11 (g/mol) - CFCL3
REAL(KIND=JPRB) :: ZAMC12                ! Molecular weight of CFC12 (g/mol) - CF2CL2
REAL(KIND=JPRB) :: ZAVGDRO               ! Avogadro's number (molecules/mole)
REAL(KIND=JPRB) :: ZGRAVIT               ! Gravitational acceleration (cm/sec2)
REAL(KIND=JPRB) :: ZAMM

REAL(KIND=JPRB) :: RAMW                  ! Molecular weight of water vapor (g/mol)
REAL(KIND=JPRB) :: RAMCO2                ! Molecular weight of carbon dioxide (g/mol)
REAL(KIND=JPRB) :: RAMO                  ! Molecular weight of ozone (g/mol)
REAL(KIND=JPRB) :: RAMCH4                ! Molecular weight of methane (g/mol)
REAL(KIND=JPRB) :: RAMN2O                ! Molecular weight of nitrous oxide (g/mol)

! Atomic weights for conversion from mass to volume mixing ratios; these
!  are the same values used in ECRT to assure accurate conversion to vmr
data ZAMD   /  28.970_JPRB    /
data ZAMW   /  18.0154_JPRB   /
data ZAMCO2 /  44.011_JPRB    /
data ZAMO   /  47.9982_JPRB   /
data ZAMCH4 /  16.043_JPRB    /
data ZAMN2O /  44.013_JPRB    /
data ZAMC11 / 137.3686_JPRB   /
data ZAMC12 / 120.9140_JPRB   /
data ZAVGDRO/ 6.02214E23_JPRB /
data ZGRAVIT/ 9.80665E02_JPRB /
data RAMW   /  0.05550_JPRB   /
data RAMCO2 /  0.02272_JPRB   /
data RAMO   /  0.02083_JPRB   /
data RAMCH4 /  0.06233_JPRB    /
data RAMN2O /  0.02272_JPRB    /


REAL(KIND=JPRB) :: ZCLEAR, ZCLOUD, ZEPSEC, ZTOTCC

INTEGER(KIND=JPIM) :: IOVLP
REAL(KIND=JPRB) :: ZHOOK_HANDLE


#include "srtm_setcoef.intfb.h"
#include "srtm_spcvrt_mcica.intfb.h"
!MPL/IM 20160915 on prend GES de phylmd
#include "clesphys.h"

!-----------------------------------------------------------------------
!-- calculate information needed ny the radiative transfer routine 

IF (LHOOK) CALL DR_HOOK('SRTM_SRTM_224GP_MCICA',0,ZHOOK_HANDLE)
ZEPSEC  = 1.E-06_JPRB
ZONEMINUS=1.0_JPRB -  ZEPSEC
ZADJI0 = RII0 / RI0
!-- overlap: 1=max-ran, 2=maximum, 3=random
IOVLP=3

!print *,'Entering srtm_srtm_224gp_mcica'

ICLDATM  = 1
INFLAG   = 2
ICEFLAG  = 3
I_LIQFLAG= 1
ITMOL    = 6
I_NSTR   = 2

DO JL = KIDIA, KFDIA
  ZRMU0=PRMU0(JL)
  IF (ZRMU0 > 0.0_JPRB) THEN

!- coefficients related to the cloud optical properties (original RRTM_SW)

!  print *,'just before SRTM_CLDPROP'

!  DO JK=1,KLEV 
!    CLDFRAC(JK) = PFRCL (JL,JK)
!    CLDDAT1(JK) = PSCLA1(JL,JK)
!    CLDDAT2(JK) = PSCLA2(JL,JK)
!    CLDDAT3(JK) = PSCLA3(JL,JK)
!    CLDDAT4(JK) = PSCLA4(JL,JK)
!    DO JMOM=0,16
!      CLDDATMOM(JMOM,JK)=PSCLMOM(JL,JMOM,JK)
!    ENDDO
!    print 9101,JK,CLDFRAC(JK),CLDDAT1(JK),CLDDAT2(JK),CLDDAT3(JK)&
!    &,CLDDAT4(JK),(CLDDATMOM(JMOM,JK),JMOM=0,NSTR)
    9101 format(1x,'srtm_srtm_224gp Cld :',I3,f7.4,7E12.5)
!  ENDDO

!  CALL SRTM_CLDPROP &
!    &( KLEV, ICLDATM, INFLAG, ICEFLAG, LIQFLAG, NSTR &
!    &, CLDFRAC, CLDDAT1, CLDDAT2, CLDDAT3, CLDDAT4, CLDDATMOM &
!    &, TAUCLDORIG, TAUCLOUD, SSACLOUD, XMOM &
!    &)

!- coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients

    DO J1=1,35
      DO J2=1,KLEV
        ZWKL(J1,J2)=0.0_JPRB 
      ENDDO
    ENDDO

    ZTBOUND=PTS(JL)
    ZPZ(0) = paph(JL,klev+1)*0.01_JPRB
    ZTZ(0) = pth (JL,klev+1)

    ZCLEAR=1.0_JPRB
    ZCLOUD=0.0_JPRB
    ZTOTCC=0.0_JPRB
    DO JK = 1, KLEV
      ZPAVEL(JK) = pap(JL,KLEV-JK+1) *0.01_JPRB
      ZTAVEL(JK) = pt (JL,KLEV-JK+1)
      ZPZ(JK)    = paph(JL,KLEV-JK+1) *0.01_JPRB
      ZTZ(JK)    = pth (JL,KLEV-JK+1)
      ZWKL(1,JK) = pq(JL,KLEV-JK+1)  *ZAMD*RAMW
      ZWKL(2,JK) = pcco2             *ZAMD*RAMCO2
      ZWKL(3,JK) = pozn(JL,KLEV-JK+1)*ZAMD*RAMO
      ZWKL(4,JK) = rn2o              *ZAMD*RAMN2O
      ZWKL(6,JK) = rch4              *ZAMD*RAMCH4
      ZAMM = (1-ZWKL(1,JK))*ZAMD + ZWKL(1,JK)*ZAMW
      ZCOLDRY(JK) = (ZPZ(JK-1)-ZPZ(JK))*1.E3_JPRB*ZAVGDRO/(ZGRAVIT*ZAMM*(1+ZWKL(1,JK)))
!    print 9200,JK,PAVEL(JK),TAVEL(JK),(WKL(JA,JK),JA=1,4),WKL(6,JK),COLDRY(JK)
      9200 format(1x,'SRTM ',I3,2F7.1,6E13.5)



    ENDDO

!  print *,'ZTOTCC ZCLEAR : ',ZTOTCC,' ',ZCLEAR

    DO IMOL=1,ITMOL
      DO JK=1,KLEV
        ZWKL(IMOL,JK)=ZCOLDRY(JK)* ZWKL(IMOL,JK)
      ENDDO
    ENDDO

!  print *,'just before SRTM_SETCOEF'

    CALL SRTM_SETCOEF &
     & ( KLEV   , ITMOL,&
     & ZPAVEL  , ZTAVEL   , ZPZ     , ZTZ     , ZTBOUND,&
     & ZCOLDRY , ZWKL,&
     & ILAYTROP, ILAYSWTCH, ILAYLOW,&
     & ZCO2MULT, ZCOLCH4  , ZCOLCO2 , ZCOLH2O , ZCOLMOL  , ZCOLN2O  , ZCOLO2 , ZCOLO3,&
     & ZFORFAC , ZFORFRAC , INDFOR  , ZSELFFAC, ZSELFFRAC, INDSELF, &
     & ZFAC00  , ZFAC01   , ZFAC10  , ZFAC11,&
     & JP      , JT       , JT1     &
     & )  
  
!  print *,'just after SRTM_SETCOEF'

!- call the radiation transfer routine
  
    DO JSW=1,KSW
      ZALBD(JSW)=PALBD(JL,JSW)
      ZALBP(JSW)=PALBP(JL,JSW)
    ENDDO

    DO JSW=1,KCOLS
      DO JK=1,KLEV        
        ZFRCL(JSW,JK) = PFRCL(JL,JSW,JK)
        ZTAUC(JK,JSW) = PTAUC(JL,JSW,JK)
        ZASYC(JK,JSW) = PASYC(JL,JSW,JK)
        ZOMGC(JK,JSW) = POMGC(JL,JSW,JK)

!---- security: might have to be moved upstream to radlswr -------
!        IF(ZTAUC(JK,JSW) == 0._JPRB) ZFRCL(JSW,JK) = 0._JPRB
!-----------------------------------------------------------------


!       IF (ZFRCL(JSW,JK) /= 0._JPRB) THEN
!          print 9002,JSW,JK,ZFRCL(JSW,JK),ZTAUC(JK,JSW),ZASYC(JK,JSW),ZOMGC(JK,JSW)
9002      format(1x,'srtm_224gp_McICA ClOPropECmodel ',2I3,f8.4,3E12.5)
!        ENDIF
      ENDDO
    ENDDO

!- mixing of aerosols
 
!  print *,'Aerosol optical properties computations'
!  DO JSW=1,KSW
!    print 9012,JSW,(JAE,RSRTAUA(JSW,JAE),RSRPIZA(JSW,JAE),RSRASYA(JSW,JAE),JAE=1,6)
    9012 format(I3,(/,I3,3E13.5))
!  ENDDO

!  DO JK=1,KLEV
!    print 9013,JK,(PAER(JL,JAE,JK),JAE=1,6)
    9013 format(1x,I3,6E12.5)
!  ENDDO

    IF (NAER == 0) THEN
      DO JSW=1,KSW
        DO JK=1,KLEV
          ZTAUA(JK,JSW)= 0.0_JPRB
          ZASYA(JK,JSW)= 0.0_JPRB
          ZOMGA(JK,JSW)= 1.0_JPRB
        ENDDO
      ENDDO
    ELSE
      DO JSW=1,KSW
        DO JK=1,KLEV
          IK=KLEV+1-JK
          ZTAUA(JK,JSW)=0.0_JPRB
          ZASYA(JK,JSW)=0.0_JPRB
          ZOMGA(JK,JSW)=0.0_JPRB
          DO JAE=1,6
            ZTAUA(JK,JSW)=ZTAUA(JK,JSW)+RSRTAUA(JSW,JAE)*PAER(JL,JAE,IK)
            ZOMGA(JK,JSW)=ZOMGA(JK,JSW)+RSRTAUA(JSW,JAE)*PAER(JL,JAE,IK) &
             & *RSRPIZA(JSW,JAE)  
            ZASYA(JK,JSW)=ZASYA(JK,JSW)+RSRTAUA(JSW,JAE)*PAER(JL,JAE,IK) &
             & *RSRPIZA(JSW,JAE)*RSRASYA(JSW,JAE)  
          ENDDO
          IF (ZOMGA(JK,JSW) /= 0.0_JPRB) THEN
            ZASYA(JK,JSW)=ZASYA(JK,JSW)/ZOMGA(JK,JSW)
          ENDIF
          IF (ZTAUA(JK,JSW) /= 0.0_JPRB) THEN
            ZOMGA(JK,JSW)=ZOMGA(JK,JSW)/ZTAUA(JK,JSW)
          ENDIF
!      print 9003,JSW,JK,ZTAUA(JK,JSW),ZOMGA(JK,JSW),ZASYA(JK,JSW)
9003  format(1x,'Aerosols ',2I3,3F10.4)
        ENDDO
      ENDDO
    ENDIF

    DO JK=1,KLEV+1
      ZBBCU(JK)=0.0_JPRB
      ZBBCD(JK)=0.0_JPRB
      ZBBFU(JK)=0.0_JPRB
      ZBBFD(JK)=0.0_JPRB
!      ZUVCU(JK)=0.0_JPRB
!      ZUVCD(JK)=0.0_JPRB
!      ZUVFU(JK)=0.0_JPRB
!      ZUVFD(JK)=0.0_JPRB
!      ZVSCU(JK)=0.0_JPRB
!      ZVSCD(JK)=0.0_JPRB
!      ZVSFU(JK)=0.0_JPRB
!      ZVSFD(JK)=0.0_JPRB
!      ZNICU(JK)=0.0_JPRB
!      ZNICD(JK)=0.0_JPRB
!      ZNIFU(JK)=0.0_JPRB
!      ZNIFD(JK)=0.0_JPRB
    ENDDO

!    print *,'just before calling STRM_SPCVRT for JL=',JL,' and ZRMU0=',ZRMU0

    CALL SRTM_SPCVRT_MCICA &
     &( KLEV   , ITMOL    , KSW    , KCOLS  , ZONEMINUS,&
     & ZPAVEL  , ZTAVEL   , ZPZ    , ZTZ    , ZTBOUND , ZALBD   , ZALBP,&
     & ZFRCL   , ZTAUC    , ZASYC  , ZOMGC  , ZTAUA   , ZASYA   , ZOMGA , ZRMU0,&
     & ZCOLDRY , ZWKL     ,&
     & ILAYTROP, ILAYSWTCH, ILAYLOW,&
     & ZCO2MULT, ZCOLCH4  , ZCOLCO2, ZCOLH2O , ZCOLMOL  , ZCOLN2O, ZCOLO2 , ZCOLO3,&
     & ZFORFAC , ZFORFRAC , INDFOR , ZSELFFAC, ZSELFFRAC, INDSELF,&
     & ZFAC00  , ZFAC01   , ZFAC10 , ZFAC11  ,&
     & JP      , JT       , JT1    ,&
     & ZBBFD   , ZBBFU    , ZBBCD  , ZBBCU )
     
!     & ZBBFD   , ZBBFU    , ZUVFD  , ZUVFU  , ZVSFD   , ZVSFU   , ZNIFD , ZNIFU,&
!     & ZBBCD   , ZBBCU    , ZUVCD  , ZUVCU  , ZVSCD   , ZVSCU   , ZNICD , ZNICU &
!     & )  

!  print *,'SRTM_SRTM_224GP before potential scaling'
!    IF (IOVLP == 3) THEN
!      DO JK=1,KLEV+1
!!      print 9004,JK,ZBBCU(JK),ZBBCD(JK),ZBBFU(JK),ZBBFD(JK)
        9004 format(1x,'Clear-sky and total fluxes U & D ',I3,4F10.3)
!        PFSUC(JL,1,JK)=ZBBCU(JK)
!        PFSUC(JL,2,JK)=ZBBCD(JK)
!        PFSUX(JL,1,JK)=ZBBFU(JK)
!        PFSUX(JL,2,JK)=ZBBFD(JK)
!      ENDDO
!    ELSE
!    print *,'SRTM_SRTM_224GP after potential scaling'
      DO JK=1,KLEV+1
        PFSUC(JL,1,JK)=ZADJI0 * ZBBCU(JK)
        PFSUC(JL,2,JK)=ZADJI0 * ZBBCD(JK)
        PFSUX(JL,1,JK)=ZADJI0 * ( (1.0_JPRB-ZCLEAR)*ZBBFU(JK)+ZCLEAR*ZBBCU(JK) )
        PFSUX(JL,2,JK)=ZADJI0 * ( (1.0_JPRB-ZCLEAR)*ZBBFD(JK)+ZCLEAR*ZBBCD(JK) )
!-- for testing only
        PFSUC(JL,1,JK)=ZADJI0 * ZBBCU(JK)
        PFSUC(JL,2,JK)=ZADJI0 * ZBBCD(JK)
        PFSUX(JL,1,JK)=ZADJI0 * ZBBFU(JK)
        PFSUX(JL,2,JK)=ZADJI0 * ZBBFD(JK)
      ENDDO
!    ENDIF

!  DO JK=1,KLEV+1
!    print 9005,JK,PFSUC(JL,1,JK),PFSUC(JL,2,JK),PFSUX(JL,1,JK),PFSUX(JL,2,JK)
    9005 format(1x,'Clear-sky and total fluxes U & D ',I3,4F10.3)
!  ENDDO
  
  ELSE
    DO JK=1,KLEV+1
      PFSUC(JL,1,JK)=0.0_JPRB
      PFSUC(JL,2,JK)=0.0_JPRB
      PFSUX(JL,1,JK)=0.0_JPRB
      PFSUX(JL,2,JK)=0.0_JPRB
    ENDDO
  ENDIF
ENDDO

!PRINT *,'OUT OF SRTM_224GP_MCICA'

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_SRTM_224GP_MCICA',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_SRTM_224GP_MCICA
