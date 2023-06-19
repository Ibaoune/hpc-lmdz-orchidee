!OPTIONS XOPT(NOEVAL)
SUBROUTINE SUPHY1(KULOUT)

!**** *SUPHY1*   - Initialize common YOMPHY1 physics land surface
!                  constants

!     Purpose.
!     --------
!           Initialize YOMPHY1, the common that contains the parameters
!           for the land surface part of the physics of the model.

!**   Interface.
!     ----------
!        *CALL* *SUPHY1(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY1

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation ARPEGE

!     Author.
!     -------
!        J.-F. Geleyn .

!     Modifications.
!     --------------
!        Original : 90-9-1
!        Modified 91-02-28 by Michel Deque (Relaxation of deep soil values)
!        Modified by R. EL Khatib : 93-04-02 Set-up defaults controled by LECMWF
!        Modified by M. Deque     : 94-10-19 4-layer soil temperature
!        Modified by D. Giard     : 94-10-07 Constants for ACSOL
!                                            Control keys : LIMC, LIMW
!        Modified by H. Douville : 95-01-13 Snow parameterization
!        Modified by D. Giard     : 95-03-09 Loop for SODELX modified
!                                   95-09-08 Constants for LSNV
!        Modified by P. Mercier   : 97-03-24 Ozone difusion + deposition
!        Modified by E. Bazile    : 97-05-05 C1 option vapour phase
!        Modified by J.M. Piriou : 97-02-26 soil inertia.
!        Modified by M. Deque     : 97-04-10 Sea-ice parameters
!        Modified by J.M. Piriou : 97-04-17 soil inertia default values.
!        Modified by D. Giard    : 97-11-13 defaults for veget. features
!        Modified by E. Bazile   : 97-12-08 Soil freezing.
!        Modified by E. Bazile   : 00-12-12 Default value for soil freezing.
!        Modified by E. Bazile   : 02-10-29 Defaults for the snow scheme LVGSN.
!        Modified by E. Bazile   : 04-02-24 Introduce NCHSP.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN   , ONLY : NULNAM
USE YOMCST   , ONLY : RPI      ,RDAY     ,RG       ,RCS
USE SURFACE_FIELDS   , ONLY : YSD_VVD, YSP_SBD
USE YOMCT0B  , ONLY : LECMWF
USE YOMPHY1  , ONLY : GF3      ,GF4      ,TREF4    ,RCTVEG   ,&
 & RGL      ,SODELX   ,GCZ0H    ,ALBGLA   ,ALBMAX   ,&
 & ALBMER   ,ALBMED   ,ALBMIN   ,ALCRIN   ,ALRCN1   ,ALRCN2   ,&
 & EA       ,EC2REF   ,EMCRIN   ,EMMGLA   ,EMMMER   ,&
 & EWFC     ,EWWILT   ,GA       ,GC1      ,GC1S1    ,&
 & GC1S2    ,GC1S3    ,GC1S4    ,GC1Y1    ,GTSVAP   ,&
 & GVEGMX   ,GLAIMX   ,GNEIMX   ,GWPIMX   ,GCGEL    ,&
 & GC2      ,GC2REF   ,GC3      ,GC31     ,GC32     ,&
 & GCONV    ,GF1      ,GWFC     ,GWLEX    ,GWLMX    ,&
 & GWWILT   ,G1B      ,G1CGSAT  ,G1C1SAT  ,G1P      ,&
 & G1WSAT   ,G2B      ,G2CGSAT  ,G2C1SAT  ,G2P      ,&
 & G2WSAT   ,G3CGSAT  ,GSNC1    ,GSNC2    ,HSOL     ,&
 & HSOLIWR  ,HSOLIT0  ,OMTPRO   ,OMWPRO   ,RC1MAX   ,&
 & RCTGLA   ,RCGMAX   ,RD1      ,RD2GLA   ,RD2MER   ,&
 & RHOMAX   ,RHOMIN   ,RSMAX    ,RTINER   ,RZ0GLA   ,&
 & RZ0MER   ,RZHZ0G   ,RZHZ0M   ,RZHGLA   ,RZHMER   ,&
 & TMERGL   ,TOEXP    ,&
 & TOLIN    ,WCRIN    ,WCRINC   ,WCRING   ,WNEW     ,&
 & WPMX     ,WSMX     ,XCRINR   ,XCRINV   ,LIMC     ,&
 & LIMW     ,LC1VAP   ,NTVGLA   ,NTVMER   ,GCGELS   ,&
 & GVEGMXS  ,GLAIMXS  ,GNEIMXS  ,ALB1     ,ALB2     ,&
 & RLAIMX   ,RLAI     ,NCHSP  
USE YOMVDOZ  , ONLY : VDHJS    ,VDHJH    ,VDHNS    ,VDHNH    ,&
 & VDPJS    ,VDPJH    ,VDPNS    ,VDPNH    ,VDEJS    ,&
 & VDEJH    ,VDENS    ,VDENH    ,VDAJS    ,VDAJH    ,&
 & VDANS    ,VDANH    ,VDNJS    ,VDNJH    ,VDNNS    ,&
 & VDNNH    ,VOZNJ    ,VOZHS    ,LRDIFOZ  ,LRDEPOZ  

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM) :: IVEG, J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
#include "abor1.intfb.h"
#include "posnam.intfb.h"

#include "namphy1.h"
#include "namvdoz.h"
!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

!        1.1 Set implicit default values

!  Snow , Sea-Ice
IF (LHOOK) CALL DR_HOOK('SUPHY1',0,ZHOOK_HANDLE)
ALCRIN=0.7_JPRB
ALRCN1=1.E-03_JPRB
ALRCN2=2.5E-03_JPRB
EMCRIN=0.98_JPRB
WCRIN=10._JPRB
TMERGL=271.23_JPRB
NCHSP=0
!  Characteristics of ice and sea
NTVGLA=2
NTVMER=1
RD2GLA=8._JPRB
RD2MER=8._JPRB
ALBGLA=.65_JPRB
ALBMER=.07_JPRB
ALBMED=.06_JPRB
EMMGLA=.97_JPRB
EMMMER=.96_JPRB
RZ0GLA=.001_JPRB
RZ0MER=.001_JPRB
RZHGLA=.001_JPRB
RZHMER=.001_JPRB
RZHZ0G=1.0_JPRB
RZHZ0M=1.0_JPRB
!  Usual - soil
HSOL=1.0E-05_JPRB
HSOLIT0=0.35_JPRB
HSOLIWR=6._JPRB
RTINER=5._JPRB
WPMX=100._JPRB
WSMX=20._JPRB
SODELX(0)=1.0_JPRB/SQRT(1.0_JPRB+2.0_JPRB*RPI)
DO J=1,9
  SODELX(J)=SODELX(J-1)*2.0_JPRB*RPI
ENDDO
!  Relaxation
OMTPRO=0._JPRB
OMWPRO=0._JPRB
!  ISBA - soil
EA=-0.54_JPRB
GA=732.42E-3_JPRB
G1B=0.137_JPRB
G2B=3.5_JPRB
G1P=0.134_JPRB
G2P=3.4_JPRB
GC1=0.5_JPRB
GC2=10._JPRB
GC3=8._JPRB
GCONV=1.E3_JPRB

G1WSAT=-1.08E-3_JPRB
G2WSAT=494.31E-3_JPRB
EWFC=0.35_JPRB
GWFC=89.0467E-3_JPRB
EWWILT=0.5_JPRB
GWWILT=37.1342E-3_JPRB
EC2REF=-0.95_JPRB
GC2REF=13.82_JPRB
G1CGSAT=-1.5571E-8_JPRB
G2CGSAT=-1.441E-8_JPRB
G3CGSAT=4.70217E-6_JPRB
G1C1SAT=5.58E-3_JPRB
G2C1SAT=84.88E-3_JPRB
GC31=5.3275_JPRB
GC32=-1.043_JPRB

RD1=1.E-2_JPRB
RC1MAX=500._JPRB
RCTGLA=5.5E-6_JPRB
RCGMAX=0.8E-5_JPRB
LIMC=.TRUE.
LIMW=.TRUE.
!  ISBA - vegetation
GF1=0.55_JPRB
GWLEX=2.0_JPRB/3._JPRB
GWLMX=0.2_JPRB
RSMAX=5000._JPRB
DO J=1,18
  GF3(J)=0.0_JPRB
  GF4(J)=0.0016_JPRB
  RCTVEG(J)=0.8E-5_JPRB
  RGL(J)=100._JPRB
  TREF4(J)=298._JPRB
ENDDO
GF3(4)=40._JPRB
RGL(4)=30._JPRB
!  ISBA - roughness length
GCZ0H(0,1)=7.5_JPRB
GCZ0H(1,1)=2.39037_JPRB
GCZ0H(2,1)=-.28583_JPRB
GCZ0H(3,1)=.01074_JPRB
GCZ0H(0,2)=0.5_JPRB
GCZ0H(1,2)=-.07028_JPRB
GCZ0H(2,2)=.01023_JPRB
GCZ0H(3,2)=-.00067_JPRB
GCZ0H(0,3)=5.0_JPRB
GCZ0H(1,3)=4.51268_JPRB
GCZ0H(2,3)=.34012_JPRB
GCZ0H(3,3)=-.05330_JPRB
GCZ0H(0,4)=0.5_JPRB
GCZ0H(1,4)=-.09421_JPRB
GCZ0H(2,4)=.01463_JPRB
GCZ0H(3,4)=-.00099_JPRB
!  ISBA - snow
ALBMAX=0.85_JPRB
ALBMIN=0.50_JPRB
RHOMAX=0.3_JPRB
RHOMIN=0.1_JPRB
TOEXP=0.24_JPRB/86400._JPRB
TOLIN=0.008_JPRB/86400._JPRB
WCRINC=70._JPRB
WCRING=10._JPRB
WNEW=10._JPRB
XCRINR=1.0_JPRB/RG
XCRINV=10000._JPRB
GSNC1=RPI/(2.22_JPRB*RCS*RDAY*1000._JPRB)
GSNC2=2.885_JPRB
!   ISBA - C1 vapour phase
LC1VAP=.TRUE.
GTSVAP=0._JPRB
GC1S1= 1.19_JPRB
GC1S2=-5.09_JPRB
GC1S3=-1.464E+2_JPRB
GC1S4= 17.86E+2_JPRB
GC1Y1=10._JPRB
! Deep Soil freezing
GVEGMX=5._JPRB
GLAIMX=30._JPRB
GNEIMX=1.8_JPRB
GWPIMX=150._JPRB
GCGEL=3.E-5_JPRB
! Surface soil freezing
GCGELS=5.E-5_JPRB
GVEGMXS=5._JPRB
GLAIMXS=30._JPRB
GNEIMXS=1.8_JPRB
!  OZONE DIFFUSION AND DEPOSITION
LRDIFOZ=.FALSE.
LRDEPOZ=.FALSE.
VOZNJ=1._JPRB
VOZHS=1.0_JPRB/86400._JPRB
DO J=1,99
  VDHJS(J)=-999._JPRB
  VDHJH(J)=-999._JPRB
  VDHNS(J)=-999._JPRB
  VDHNH(J)=-999._JPRB
  VDPJS(J)=-999._JPRB
  VDPJH(J)=-999._JPRB
  VDPNS(J)=-999._JPRB
  VDPNH(J)=-999._JPRB
  VDEJS(J)=-999._JPRB
  VDEJH(J)=-999._JPRB
  VDENS(J)=-999._JPRB
  VDENH(J)=-999._JPRB
  VDAJS(J)=-999._JPRB
  VDAJH(J)=-999._JPRB
  VDANS(J)=-999._JPRB
  VDANH(J)=-999._JPRB
  VDNJS(J)=-999._JPRB
  VDNJH(J)=-999._JPRB
  VDNNS(J)=-999._JPRB
  VDNNH(J)=-999._JPRB
ENDDO
! New snow scheme (LVGSN)
ALB1=0.87_JPRB
ALB2=0.84_JPRB
RLAIMX=7._JPRB
RLAI=3._JPRB

!        1.2 Modify default values according to LECMWF

IF (LECMWF) THEN
ELSE
ENDIF

!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

! Ce qui concerne NAMPHY1 et NAMVDOZ commente par MPL le 15.04.09
!CALL POSNAM(NULNAM,'NAMPHY1')
!READ(NULNAM,NAMPHY1)
!CALL POSNAM(NULNAM,'NAMVDOZ')
!READ(NULNAM,NAMVDOZ)

!        2.5   Check consistency
!              -----------------
IF (GC1Y1 > 60._JPRB) THEN
  CALL ABOR1 ('GC1Y1 FOR C1-VAPOUR PHASE IS BIGGER THAN 60.')
ENDIF

IF ((YSD_VVD%NUMFLDS < 8).AND.(RZHZ0M /= 1.0_JPRB)) THEN
  CALL ABOR1('YSD_VVD%NUMFLDS<8 IMPLIES RZHZ0M=1.0_JPRB !...')
ENDIF

!     ------------------------------------------------------------------

!*       3.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMPHY1 '')')

WRITE(UNIT=KULOUT,FMT='(&
 & '' ALCRIN ='',E10.4,'' ALRCN1 ='',E10.4,'' ALRCN2 ='',E10.4 &
 & ,'' EMCRIN ='',E10.4,'' WCRIN ='',E10.4,'' TMERGL ='',E10.4)')&
 & ALCRIN,ALRCN1,ALRCN2,EMCRIN,WCRIN,TMERGL  

WRITE(UNIT=KULOUT,FMT='(&
 & '' NCHSP ='',I3,&
 & '' NTVMER ='',I3,'' NTVGLA ='',I3 &
 & ,'' RD2MER ='',E10.4,'' ALBMER='',E10.4,'' EMMMER ='',E10.4,/&
 & ,'' RD2GLA ='',E10.4,'' ALBGLA='',E10.4,'' EMMGLA ='',E10.4,/&
 & ,'' RZ0MER ='',E10.4,'' RZHMER='',E10.4 &
 & ,'' RZ0GLA ='',E10.4,'' RZHGLA='',E10.4)')&
 & NCHSP,NTVMER,NTVGLA,RD2MER,ALBMER,EMMMER &
 & ,RD2GLA,ALBGLA,EMMGLA,RZ0MER,RZHMER,RZ0GLA,RZHGLA  
WRITE(UNIT=KULOUT,FMT='(&
 & '' RTINER ='',E10.4,'' HSOL ='',E10.4,'' HSOLIT0 ='',E10.4 &
 & ,'' HSOLIWR ='',E10.4,'' WPMX ='',E10.4 &
 & ,'' WSMX ='',E10.4,'' OMTPRO ='',E10.4,'' OMWPRO ='',E10.4)')&
 & RTINER,HSOL,HSOLIT0,HSOLIWR,WPMX,WSMX,OMTPRO,OMWPRO  

WRITE(UNIT=KULOUT,FMT='(&
 & '' EA ='',E10.4,'' GA ='',E10.4,'' G1B ='',E10.4 &
 & ,'' G2B ='',E10.4,'' G1P ='',E10.4,'' G2P = '',E10.4,/&
 & ,'' GC1 ='',E10.4,'' GC2 ='',E10.4,'' GC3 ='',E10.4 &
 & ,'' GCONV ='',E10.4)')&
 & EA,GA,G1B,G2B,G1P,G2P,GC1,GC2,GC3,GCONV  

WRITE(UNIT=KULOUT,FMT='(&
 & '' G1WSAT ='',E10.4,'' G2WSAT ='',E10.4,'' EWFC ='',E10.4 &
 & ,'' GWFC ='',E10.4,/,'' EWWILT ='',E10.4,'' GWWILT ='',E10.4 &
 & ,'' EC2REF ='',E10.4,'' GC2REF ='',E10.4,/&
 & ,'' G1CGSAT ='',E10.4,'' G2CGSAT ='',E10.4,'' G3CGSAT =''&
 & ,E10.4,'' G1C1SAT ='',E10.4,'' G2C1SAT ='',E10.4)')&
 & G1WSAT,G2WSAT,EWFC,GWFC,EWWILT,GWWILT,EC2REF,GC2REF,&
 & G1CGSAT,G2CGSAT,G3CGSAT,G1C1SAT,G2C1SAT  

WRITE(UNIT=KULOUT,FMT='(&
 & '' RD1 ='',E10.4,'' RC1MAX ='',E10.4,'' RCTGLA ='',E10.4 &
 & ,'' RCGMAX ='',E10.4,''  LIMC ='',L2,'' LIMW ='',L2)')&
 & RD1,RC1MAX,RCTGLA,RCGMAX,LIMC,LIMW  

WRITE(UNIT=KULOUT,FMT='(&
 & '' GC1S1 ='',E10.4,'' GC1S2 ='',E10.4,'' GC1S3 ='',E10.4 &
 & ,'' GC1S4 ='',E10.4,''  GC1Y1 ='',E10.4,'' LC1VAP ='',L2 &
 & ,'' GTSVAP ='',E10.4)')&
 & GC1S1,GC1S2,GC1S3,GC1S4,GC1Y1,LC1VAP,GTSVAP  

WRITE(UNIT=KULOUT,FMT='(&
 & '' GCGEL ='',E10.4,'' GVEGMX ='',E10.4,'' GLAIMX ='',E10.4 &
 & ,'' GWPIMX ='',E10.4,'' GNEIMX  ='',E10.4)')&
 & GCGEL,GVEGMX,GLAIMX,GWPIMX,GNEIMX  

WRITE(UNIT=KULOUT,FMT='(&
 & '' GCGELS  ='',E10.4,'' GVEGMXS ='',E10.4 &
 & ,'' GLAIMXS ='',E10.4,'' GNEIMXS ='',E10.4)')&
 & GCGELS,GVEGMXS,GLAIMXS,GNEIMXS  

WRITE(UNIT=KULOUT,FMT='(&
 & '' ALB1  ='',E10.4,'' ALB2 ='',E10.4 &
 & ,'' RLAIMX ='',E10.4,'' RLAI ='',E10.4)')&
 & ALB1,ALB2,RLAIMX,RLAI  

WRITE(UNIT=KULOUT,FMT='('' GCZ0H ='',/,4(1X,4E11.4,/))')GCZ0H

WRITE(UNIT=KULOUT,FMT='(&
 & '' GF1 ='',E10.4,'' GWLEX ='',E10.4,'' GWLMX ='',E10.4 &
 & ,'' RSMAX ='',E10.4,/&
 & ,'' GF3 ='',/,2(1X,9E11.4,/),'' GF4 ='',/,2(1X,9E11.4,/)&
 & ,'' RCTVEG ='',/,2(1X,9E11.4,/),'' RGL ='',/,2(1X,9E11.4,/)&
 & ,'' TREF4 ='',/,2(1X,9E11.4,/))')&
 & GF1,GWLEX,GWLMX,RSMAX,GF3,GF4,RCTVEG,RGL,TREF4  

WRITE(UNIT=KULOUT,FMT='(&
 & '' ALBMAX = '',E10.4,'' ALBMIN = '',E10.4 &
 & ,'' RHOMAX = '',E10.4,'' RHOMIN = '',E10.4,/&
 & ,'' TOEXP = '',E10.4,'' TOLIN = '',E10.4 &
 & ,'' WCRINC = '',E10.4,'' WCRING = '',E10.4,/&
 & ,'' WNEW = '',E10.4 &
 & ,'' XCRINR = '',E10.4,'' XCRINV = '',E10.4)')&
 & ALBMAX,ALBMIN,RHOMAX,RHOMIN,TOEXP,TOLIN,WCRINC,WCRING,WNEW,&
 & XCRINR,XCRINV  

WRITE(UNIT=KULOUT,FMT='('' SODELX = ''/5E11.4/5E11.4)') SODELX
IF(YSP_SBD%NLEVS > 9)  CALL ABOR1(' TOO MANY SOIL LAYERS !')

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMVDOZ '')')

WRITE(UNIT=KULOUT,FMT='('' LRDIFOZ ='',L2,'' LRDEPOZ ='',L2)')LRDIFOZ,LRDEPOZ

IF (LRDIFOZ.AND.LRDEPOZ) THEN
  WRITE(UNIT=KULOUT,FMT='('' VOZNJ ='',F10.5,'' VOZHS ='',F10.5)')VOZNJ,VOZHS
  DO IVEG=1,99
    IF((VDHJS(IVEG) /= -999._JPRB).OR.(VDHJH(IVEG) /= -999._JPRB).OR.&
       & (VDHNS(IVEG) /= -999._JPRB).OR.(VDHNH(IVEG) /= -999._JPRB).OR.&
       & (VDPJS(IVEG) /= -999._JPRB).OR.(VDPJH(IVEG) /= -999._JPRB).OR.&
       & (VDPNS(IVEG) /= -999._JPRB).OR.(VDPNH(IVEG) /= -999._JPRB).OR.&
       & (VDEJS(IVEG) /= -999._JPRB).OR.(VDEJH(IVEG) /= -999._JPRB).OR.&
       & (VDENS(IVEG) /= -999._JPRB).OR.(VDENH(IVEG) /= -999._JPRB).OR.&
       & (VDAJS(IVEG) /= -999._JPRB).OR.(VDAJH(IVEG) /= -999._JPRB).OR.&
       & (VDANS(IVEG) /= -999._JPRB).OR.(VDANH(IVEG) /= -999._JPRB).OR.&
       & (VDNJS(IVEG) /= -999._JPRB).OR.(VDNJH(IVEG) /= -999._JPRB).OR.&
       & (VDNNS(IVEG) /= -999._JPRB).OR.(VDNNH(IVEG) /= -999._JPRB))THEN  
      WRITE(UNIT=KULOUT,FMT='('' IVEJ ='',I2 &
       & ,'' VDHJS ='',E9.3,'' VDHJH ='',E9.3,'' VDHNS ='',E9.3 &
       & ,'' VDHNH ='',E9.3)')&
       & IVEG,VDHJS(IVEG),VDHJH(IVEG),VDHNS(IVEG),VDHNH(IVEG)  
      WRITE(UNIT=KULOUT,FMT='('' IVEJ ='',I2 &
       & ,'' VDPJS ='',E9.3,'' VDPJH ='',E9.3,'' VDPNS ='',E9.3 &
       & ,'' VDPNH ='',E9.3)')&
       & IVEG,VDPJS(IVEG),VDPJH(IVEG),VDPNS(IVEG),VDPNH(IVEG)  
      WRITE(UNIT=KULOUT,FMT='('' IVEJ ='',I2 &
       & ,'' VDEJS ='',E9.3,'' VDEJH ='',E9.3,'' VDENS ='',E9.3 &
       & ,'' VDENH ='',E9.3)')&
       & IVEG,VDEJS(IVEG),VDEJH(IVEG),VDENS(IVEG),VDENH(IVEG)  
      WRITE(UNIT=KULOUT,FMT='('' IVEJ ='',I2 &
       & ,'' VDAJS ='',E9.3,'' VDAJH ='',E9.3,'' VDANS ='',E9.3 &
       & ,'' VDANH ='',E9.3)')&
       & IVEG,VDAJS(IVEG),VDAJH(IVEG),VDANS(IVEG),VDANH(IVEG)  
      WRITE(UNIT=KULOUT,FMT='('' IVEJ ='',I2 &
       & ,'' VDNJS ='',E9.3,'' VDNJH ='',E9.3,'' VDNNS ='',E9.3 &
       & ,'' VDNNH ='',E9.3)')&
       & IVEG,VDNJS(IVEG),VDNJH(IVEG),VDNNS(IVEG),VDNNH(IVEG)  
    ENDIF
  ENDDO
ENDIF
!*
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPHY1',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHY1
