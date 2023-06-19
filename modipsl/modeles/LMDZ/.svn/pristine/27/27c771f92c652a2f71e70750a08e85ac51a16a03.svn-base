
!OPTIONS XOPT(NOEVAL)
SUBROUTINE SUPHY0(KULOUT)

!**** *SUPHY0*   - Initialize common YOMPHY0 physics atmospheric
!                  constants

!     Purpose.
!     --------
!           Initialize YOMPHY0, the common that contains the parameters
!           for the atmospheric part of the physics of the model.

!**   Interface.
!     ----------
!        *CALL* *SUPHY0(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY0

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
!        Ajout de GWDAMP (J.-F. Geleyn) : 91-2-2
!        Modified by Michel Deque 91-04-01 (param. for convect. clouds)
!        Ajout de HOBST, NPCLO1/2, XNBMAX et REVGSL, remplacement de
!           GWDCOE par GWDSE (J.-F. Geleyn, M. Deque, L. Labbe) : 92-4-8
!        Ajout de VZ0CM (E. Bazile) : 92-3-27
!        Modified by C. Castejon and E. Gerard 92-02-28 (stat. clouds)
!        Modified by R. EL Khatib : 93-04-02 Set-up defaults controled by LECMWF
!        Modified by M. Deque : 95-03-20 (USURIC and USUPRC)
!        Modified by J.F. Geleyn C. Bossuet : 95-12-20 (GWDBC)
!        Modified by Luc Gerard: 97-02-26 (entrainement de la qte de mvt horizontale).
!        Modified by J.M. Piriou: 97-02-28 (schema de nebulosite ACNEBN).
!        Modified by J.M. Piriou: 97-04-17 (valeurs par defaut).
!        Modified by J.L. Ricard : seuil (SCO) sur les precip. conv.
!                                  modif de la turbulence residuelle
!        Modified by J.M. Piriou: 97-08-21 (introduce Xu-Randall cloudiness).
!        Modified by M. Deque   : 97-10-21 (introduce GWD lift).
!        Modified by M. Deque   : 98-02-21 (new mixing length profile).
!        Modified by V. Lorant  : 98-08-05 (new mixing length profile).
!        Modified by J.M. Piriou: 98-02-11 (introduce downdrafts tuning parameters).
!        Modified by J.M. Piriou: 98-03-10 (introduce GRCVPP)
!        Modified by L. Gerard  : 98-11-30 (TUDGP, TDDGP, GCOMOD)
!        Modified by R. EL Khatib :98-12-14 Remove LRDSPIL
!        Modified by J.M. Piriou: 99-01-04 (introduce GCVADS, GCVBETA)
!        Modified by J.M. Piriou: 99-06-18 (introduce GCVPSI, GCVALFA, USURICL AND USURICE.
!                                 Change default value for GDDEVA).
!        Modified by J.M. Piriou: 2000-08-23 (new use of the Richardson critical number (USURID)).
!        Modified by J.M. Piriou: 2000-08-23 (cloud core buoyancy as a fraction of an undilute plume (GCVNU)).
!        Modified by J.M. Piriou: 2000-10-06 (exponent USURIDE).
!        Modified by E. Bazile  : 2000-12-12 CYCORA's default value.
!        Modified by J.M. Piriou: 2001-04-05 (introduce GCVMLT).
!        Modified by Y. Bouteloup:2002-06-14 (introduce RCVEVAP).
!        Modified by F. Bouyssel: 2002-06-25 (introduce UTILGUST, RRGAMMA, RRSCALE).
!        Modified by J.M. Piriou: 2002-08-19 (introduce GPBLH*).
!        Modified by D. Banciu:   2002-12-09 (introduce GCVPSIE).
!        R. El Khatib : 2001-08-07 Pruning options
!        J.M. Piriou  : 2002-01-10 set default values to operational ones.
!        03-06, move rnlcurv rnegat into yomphy0 (F. Bouyssel, C. Fischer)
!        E. Bazile : 2004-02-24 (introduce EDK).
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        E. Bazile : 2004-06-30 (inroduces XKLM).
!        R. Brozkova : 2004-11 modifs for Xue-Randall cloud. scheme
!        P. Marquet & F.Bouyssel : 2004-08-18 (Lopez)
!        F.Bouyssel : 2005-07-18 (new defaults and new parameters for Lopez)
!        F.Bouyssel : 2006-01-25 (new defaults for Lopez!)
!        F. Vana    : 2006-01-30 tunables for pTKE.
!        R.Brozkova : 2006-03-03 : tuning constants for mixing lengths and Charnock formulae
!        A.Alias    : 2006-03-10 renaming KDN in RKDN
!        Modified by GMGEC/EAC  : 2006-05 list of modif.
!                    V. Lorant  : 99-01-05 (new mixing length profile).
!                    M. Deque   : 00-03-21 (new entrainment rate VL  ).
!                    P. Marquet : 2002-11-05 TRENTRVL=0. if not ACPBLH
!                    P. Marquet : 2004-05-27 TFVR and TFVS for ADVPRC.
!                    P. Marquet : 2004-10-13 RAUTEFR for snow (ACMICRO).
!                    P. Marquet : 2004-10-14 RAUTSBET for snow (ACMICRO).
!                    A. Alias   : 2005-06-23 (param. for ACCVIMPGDY)
!        E. Bazile & P. Marquet : 2006-04-11 AGRE1,AGRERICR,AJBUMIN,RCOFLM 
!                                  for LPBLE.
!        F.Bouyssel : 2006-10-30 RQLCV,RQICVMAX,RQICVMIN,RHEVAP
!        M. Bellus  : 03-Oct-2006 : ALARO-0 phasing (defaults for prognostic
!                                   convection, pTKE and PIL microphysics)
!        A.Alias     : 2006-07   FEVAPC added for ACCVIMPGY (JF Gueremy)
! ------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN   , ONLY : NULNAM
USE YOMCT0B  , ONLY : LECMWF
USE YOMSIMPHL, ONLY : LSIMPH   ,LSMOOTHD ,LSMOOTHA
USE YOMPHY   , ONLY : LCVCAS   ,LCVLIS   ,LCVRA    ,LMPHYS   ,LCONDWT  ,&
 & LPROCLD  ,LNEBCO   ,LNEBR    ,LECT     ,LNEBGR   ,LCVRAV3 ,LPIL
USE YOMPHY0  , ONLY : TYM      ,NPCLO1   ,NPCLO2   ,AERCS1   ,AERCS3   ,&
 & AERCS5   ,ALMAV    ,BEDIFV   ,ECMNP    ,EDB      ,EDC      ,&
 & EDD      ,EVAP     ,FONT     ,GALP     ,GCCSV    ,GCOMOD   ,&
 & GCVADS   ,GCVALFA  ,GCVBETA  ,GCVPSI   ,GDDEVA   ,GDDSDE   ,&
 & GRCVPP   ,GWDAMP   ,GWDBC    ,GWDCCO   ,GWDCD    ,GWDLT    ,&
 & AHCLPV   ,GWDSE    ,HOBST    ,HUCOE    ,HUCOE2   ,HUTIL    ,&
 & HUTIL1   ,HUTIL2   ,&
 & QSSC     ,QSMIN    ,QSNEBC   ,QSNEBS   ,QSSUSC   ,QSSUSS   ,&
 & QXRAL    ,QXRDEL   ,QXRR     ,REVGSL   ,RTCAPE   ,SCO      ,&
 & SENSL    ,SNNBCO   ,SPNBCO   ,SXNBCO   ,TCA      ,TCT      ,&
 & TCW      ,TENTR    ,TENTRX   ,TDDGP    ,TUDGP    ,TURB     ,&
 & TVF      ,UHDIFV   ,USDMLT   ,USUPRC   ,USURIC   ,USURID   ,&
 & USURIDE  ,USURICE  ,USURICL  ,VCHRNK   ,VKARMN   ,VZ0CM    ,&
 & XNBMAX   ,RICRLM   ,XBLM     ,XMINLM   ,XMAXLM   ,XWSALM   ,&
 & XWSBLM   ,GCVNU    ,GCVMLT   ,RCVEVAP  ,GPBLHK0  ,GPBLHRA  ,&
 & UTILGUST ,RRGAMMA  ,RRSCALE  ,GCVPSIE  ,QSUSXC   ,QSUSXS   ,&
 & ETACUT   ,RNEGAT   ,RNLCURV  ,QSSUSV   ,QXRHX    ,GCISMIN  ,&
 & GWDPROF  ,GWDVALI  ,RCIN     ,EDK      ,XKLM     ,RPHI0    ,&
 & RPHIR    ,QXRTGH   ,ADISE    ,ADISI    ,AECLS3   ,AECLS4   ,&
 & AKN      ,ALD      ,ALPHAE   ,ALPHAT   ,ECTMIN   ,UCWSTAR  ,&
 & UDECT    ,USHEARM  ,UPRETMIN ,UPRETMAX ,ARSCH    ,ARSCQ    ,&
 & ARSC1    ,ARSB2    ,ACBRPHIM ,ALMAVE   ,RICRET   ,STTBMIN  ,&
 & AGREKE   ,RDTFAC   ,RAUTEFS  ,RNINTR   ,RNINTS   ,RQLCR    ,&
 & RQICRMAX ,RQICRMIN ,RACCEF   ,RRIMEF   ,RHCRIT1  ,RHCRIT2  ,&
 & RETAMIN  ,TFVR     ,TFVS     ,RAUTEFR  ,RAUTSBET ,GRHCMOD  ,&
 & RQICRT1  ,RQICRT2  ,RQICRSN  ,RQCRNS   ,RFACNSM  ,RAGGEF   ,&
 & RQLCV    ,RQICVMAX ,RQICVMIN ,RHEVAP   ,&
 & AGRE1    ,AGRERICR ,AJBUMIN  ,RCOFLM   ,&
 & A0ML_AU  ,A0ML_AT  ,A0ML_BU  ,A0ML_BT  ,VZIUSTAR0,&
 & TENTRVL  ,TRENTRV  ,UETEPS   ,UPRECLP  ,&
 & ARSC2    ,ARSCT    ,AGRE2    ,AGREF    ,&
 & AJ1PEPS  ,AJ1MEPS  ,NAJITER  ,&
 & ALFX     ,TCTC     ,TVFC     ,GAMAP1   ,RKDN      ,&
 & VVN      ,VVX      ,FENTRT   ,HCMIN    ,FQLIC    ,FNEBC    ,&
 & NUPTKE   ,GAMTKE   ,RCOLL    ,RFALLL   ,TDDBU    ,&
 & TDDFR    ,TUDBU    ,TUDFR    ,GCVACHI  ,GCVALMX  ,GCVADMW  ,&
 & GCVBEE   ,GCVEEX   ,ECMNPI   ,GFRIC    ,GCVSQDN  ,GCVSQDR  ,&
 & GCVSQDCX ,GRRINTE  ,GRRMINA  ,GDDBETA  ,GDDEVF   ,GDDWPF   ,&
 & TENTRD   ,RDPHIC   ,GWBFAUT  ,RWBF1    ,RWBF2    ,RAUITN   ,&
 & RAUITX   ,RAUIUSTE ,RSMDNEBX ,RSMDTX   ,NSMTPA   ,NSMTPB   ,&
 & FEVAPC


IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "surhcri.intfb.h"

#include "namphy0.h"

!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

!        1.1 Set implicit default values

IF (LHOOK) CALL DR_HOOK('SUPHY0',0,ZHOOK_HANDLE)
A0ML_AU=4.5_JPRB
A0ML_AT=5.0_JPRB
A0ML_BU=3.0_JPRB
A0ML_BT=0.8_JPRB
AERCS1=0.2726_JPRB
AERCS3=-0.4239_JPRB
AERCS5=0.3595_JPRB
ALMAV=300._JPRB
BEDIFV=0.05_JPRB
ECMNP=3000._JPRB
EDB=5._JPRB
EDC=5._JPRB
EDD=5._JPRB
EDK=1.0_JPRB
ETACUT=1.0_JPRB
EVAP=4.8E+06_JPRB
FONT=2.4E+04_JPRB
GALP=0.53_JPRB
GCISMIN=6.7E-05_JPRB
GCCSV=0._JPRB
GCOMOD=1._JPRB
GCVADS=0.8_JPRB
GCVALFA=4.5E-05_JPRB
GCVBETA=0.2_JPRB
GCVMLT=0.00016_JPRB
GCVNU=2.5E-05_JPRB
GCVPSI=0.5_JPRB
GCVPSIE=0.0_JPRB
GDDEVA=0.25_JPRB
GDDSDE=0.5_JPRB
GPBLHK0=0.25_JPRB
GPBLHRA=4._JPRB
GRCVPP=1._JPRB
GWDAMP=0.6_JPRB
GWDBC=2._JPRB
GWDCCO=1._JPRB
GWDCD=6._JPRB
GWDLT=0._JPRB
GWDSE=3.5E-03_JPRB
GWDPROF=1._JPRB
GWDVALI=0._JPRB
HOBST=3._JPRB
HUCOE=2._JPRB
HUCOE2=0.4_JPRB
HUTIL=1.8_JPRB
HUTIL1=-0.6_JPRB
HUTIL2=1.1_JPRB
NPCLO1=1
NPCLO2=1
QSSC=1600._JPRB
QSMIN=1.E-4_JPRB
QSNEBC=26000._JPRB
QSNEBS=0.7_JPRB
QSSUSC=1._JPRB
QSSUSS=0.25_JPRB
QSSUSV=0._JPRB
QSUSXC=3.3E-05_JPRB
QSUSXS=3.3E-05_JPRB
RPHI0=0._JPRB
RPHIR=1750._JPRB
QXRAL=10000._JPRB
QXRDEL=0._JPRB
QXRHX=1._JPRB
QXRR=0.5_JPRB
QXRTGH=3.5_JPRB
RCIN=0._JPRB
RCVEVAP=0._JPRB
REVGSL=80._JPRB
RTCAPE=10800._JPRB
SCO=-20._JPRB
SENSL=1._JPRB
SNNBCO=0._JPRB
SPNBCO=3000._JPRB
SXNBCO=0.5_JPRB
TCA=1._JPRB
TCT=1.E-4_JPRB
TCW=8.E-4_JPRB
TENTR=2.5E-06_JPRB
TENTRX=8.E-05_JPRB
TDDGP=0.8_JPRB
TUDGP=0.8_JPRB
TURB=1._JPRB
TVF=1._JPRB
TYM(1)=0.92_JPRB
TYM(2)=0.74_JPRB
TYM(3)=16.6_JPRB
TYM(4)=10.1_JPRB
TYM(5)=0.08_JPRB
UHDIFV=8.E-04_JPRB
USDMLT=1.25E+04_JPRB
USUPRC=0.0_JPRB
USURIC=1.0_JPRB
USURID=0.035_JPRB
USURIDE=1.0_JPRB
USURICE=0.5_JPRB
USURICL=4._JPRB
UTILGUST=0.125_JPRB
VCHRNK=0.021_JPRB
VKARMN=0.4_JPRB
VZ0CM=1.5E-04_JPRB
VZIUSTAR0=0._JPRB
XNBMAX=1._JPRB
AHCLPV=1000._JPRB
RICRLM=0.5_JPRB
RRGAMMA=0.8_JPRB
RRSCALE=1.15E-4_JPRB
XBLM=6.5_JPRB
XKLM=1.0_JPRB
XMAXLM=3000._JPRB
XMINLM=500._JPRB
XWSALM=0.1_JPRB
XWSBLM=7.0_JPRB
RNEGAT = -7.E-05_JPRB
RNLCURV = 7.E+04_JPRB

!     - - - - - - - - - - - - - -
!     The old Convective scheme :
!     - - - - - - - - - - - - - -
TRENTRV=1._JPRB
TENTRVL=-1.0_JPRB

IF ( LCONDWT.AND.LPROCLD ) THEN
  RDTFAC=0.5_JPRB
ELSE
  RDTFAC=1.0_JPRB
ENDIF

!        - - - - - - - - - -
!        Lopez Microphysics :
!        - - - - - - - - - -
RAUTEFR=1.E-03_JPRB
RAUTEFS=1.E-03_JPRB
RAUTSBET=0.025_JPRB
RNINTR=8.E+06_JPRB
RNINTS=2.E+06_JPRB
RQLCR=2.E-04_JPRB
RQICRMAX=0.3E-04_JPRB
RQICRMIN=0.2E-06_JPRB
RQLCV=2.E-04_JPRB
RQICVMAX=0.3E-04_JPRB
RQICVMIN=0.2E-06_JPRB
RQICRT1=-80._JPRB
RQICRT2=30._JPRB
RQICRSN=0.5_JPRB
RQCRNS=0.03_JPRB
RACCEF=1._JPRB
RAGGEF=0.2_JPRB
RRIMEF=1._JPRB
RHEVAP=0.0_JPRB
RHCRIT1=0.5_JPRB
RHCRIT2=0.91_JPRB
RETAMIN=0.4_JPRB
RFACNSM=1.4_JPRB
TFVR=5.0_JPRB
TFVS=0.6_JPRB
GRHCMOD=0.3_JPRB


!--------------------------
! PIL MICROPHYSICS
!--------------------------
RDPHIC=10000._JPRB
GWBFAUT=15._JPRB
! two constants for ACPLUIE_PROG:
RWBF1=300._JPRB
RWBF2=4._JPRB
!
RAUITN=233.15_JPRB
RAUITX=263.15_JPRB
RAUIUSTE=0.025_JPRB
RSMDNEBX=0.2_JPRB
RSMDTX=1.0_JPRB
NSMTPA=2
NSMTPB=3
RCOLL=6.9E-03_JPRB
RFALLL=1.0_JPRB
!--------------------------
! Prognostic convection physical parameters
!--------------------------
TUDBU=0.5_JPRB
TDDBU=0.5_JPRB
TUDFR=0.0012_JPRB
TDDFR=0.0006_JPRB
GCVALMX=0.95_JPRB
! Pas d'activite historique: defaut tres haut,
! pour que KUO joue tout seul:
GCVACHI=1.E9_JPRB
GCVADMW =0
! Explicit entrainment:
GCVBEE=0._JPRB ! 0.2_JPRB
GCVEEX=1._JPRB
ECMNPI=3000._JPRB
GFRIC=-1.0_JPRB !1.E-3_JPRB
! Sqeezing:
GCVSQDN=0.01_JPRB
GCVSQDR=0.8_JPRB
GCVSQDCX=1.0_JPRB ! No squeezing
! PRECIPITATING AREA PARAMETERS (aplpar):
GRRINTE=2._JPRB
GRRMINA=1.E-5_JPRB
! DD explicit detrainment:
GDDBETA=0.2_JPRB
! ACMODO DD PARAMETERS:
GDDEVF=0.5_JPRB
GDDWPF=0._JPRB
! DD ENTRAINMENT RATE:
TENTRD=1.E-4_JPRB

!     - - - - - - - -
!     For TKE scheme :
!     - - - - - - - -
ADISE=-0.5_JPRB
ADISI=1.5_JPRB
AECLS3=3.75_JPRB
AECLS4=0.3_JPRB
AKN=0.2_JPRB
ALD=1.4_JPRB
ALPHAE=1.0_JPRB
ALPHAT=1.0_JPRB
ECTMIN=1.E-10_JPRB

UCWSTAR=1.0_JPRB/3._JPRB
UDECT=5._JPRB
USHEARM=1.E-04_JPRB
UPRETMIN=60000._JPRB
UPRETMAX=97500._JPRB
ARSCH=4._JPRB
ARSCQ=1.2_JPRB
ARSC1=2.0_JPRB/(3._JPRB*ARSCH*ARSCQ)
ARSB2=3._JPRB*ARSCH*ARSC1/2.0_JPRB
ACBRPHIM=2.2_JPRB
ALMAVE=0._JPRB
RICRET=0.195_JPRB

STTBMIN=SQRT(3._JPRB)
! A minimum value for the adimentional jump in boyancy : d(Theta)/Theta
AJBUMIN=0.005_JPRB
! The "Master Length" is equal to "RCOFLM*Z_PBL"
RCOFLM=0.085_JPRB


! TKE (P.Marquet)
UETEPS=1.0_JPRB 

!     - - - - - - - - - - - - - - - -
!     For dry conv. adjustment scheme :
!     - - - - - - - - - - - - - - - -
AJ1MEPS=0.99_JPRB
AJ1PEPS=10.0_JPRB
NAJITER=30

!     - - - - - - - - - - - - - - - - - - - - - - -
!     For Grenier (2000) top-PBL entrainment scheme :
!     - - - - - - - - - - - - - - - - - - - - - - -
AGRE1=0.16_JPRB
AGRE2=15._JPRB
AGREF=0.8_JPRB
AGREKE=5.0_JPRB
AGRERICR=50._JPRB

! Pseudo prognostic TKE scheme
NUPTKE=0.52_JPRB
GAMTKE=0.5_JPRB



!     - - - - - - -
!     For ACCVIMPGY :
!     - - - - - - -
ALFX=10.E-02_JPRB
TCTC=1.60E-04_JPRB
TVFC=1._JPRB
GAMAP1=1.5_JPRB
RKDN=30.E-06_JPRB
VVN=0.0_JPRB
VVX=-45._JPRB
FENTRT=2.5_JPRB
HCMIN=0.0_JPRB
FQLIC=2.5_JPRB
FNEBC=25.0_JPRB
FEVAPC=3.5_JPRB

!        1.2 Modify default values according to LECMWF

IF (LECMWF) THEN
ELSE
ENDIF

!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

! Ce qui concerne NAMPHY0 commente par MPL le 15.04.09
!CALL POSNAM(NULNAM,'NAMPHY0')
!READ(NULNAM,NAMPHY0)

IF(.NOT.(LNEBCO.AND.(LNEBR.OR.LNEBGR.OR.LECT)) .AND.LCVRAV3) THEN
  WRITE(UNIT=KULOUT,FMT='(A)') ' '
  WRITE(UNIT=KULOUT,FMT='(A)') ' !'
  WRITE(UNIT=KULOUT,FMT='(A)') ' ! TENTRVL is set to 0. in SUPHY0  !!'
  WRITE(UNIT=KULOUT,FMT='(A)') ' !'
  WRITE(UNIT=KULOUT,FMT='(A)') ' '
  TENTRVL = 0.0_JPRB
ENDIF
!     ------------------------------------------------------------------

!*       3.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMPHY0 '')')
WRITE(UNIT=KULOUT,FMT='('' AERCS1 = '',E11.4,'' AERCS3 = '',E11.4 &
 & ,'' AERCS5 = '',E11.4,'' ALMAV = '',E11.4,'' ECMNP = '',E11.4 &
 & ,'' EDB = '',E11.4,/,'' EDC = '',E11.4,'' EDD = '',E11.4 &
 & ,'' EDK = '',E11.4,'' ETACUT = '',E11.4 &
 & ,'' EVAP = '',E11.4,'' FONT = '',E11.4,'' GWDAMP = '',E11.4 &
 & ,'' GWDSE = '',E11.4,'' GWDBC = '',E11.4,/&
 & ,'' GWDCD = '',E11.4,/&
 & ,'' GWDPROF = '',E11.4,'' GWDVALI = '',E11.4,/&
 & ,'' HUCOE = '',E11.4,'' HUTIL = '',E11.4 ,'' HUTIL1 = '',E11.4 ,'' HUTIL2 = '',E11.4 &
 & ,'' VCHRNK = '',E11.4,'' VKARMN = '',E11.4 &
 & ,'' SNNBCO = '',E11.4,/,'' SPNBCO = '',E11.4 &
 & ,'' SXNBCO = '',E11.4,'' HOBST = '',E11.4,'' NPCLO1 = '',I4 &
 & ,'' NPCLO2 = '',I4,'' RCIN = '',E11.4,'' RCVEVAP = '',E11.4 &
 & ,'' REVGSL = '',E11.4,/&
 & ,'' RTCAPE= '',E11.4,'' GCOMOD= '',E11.4,/&
 & ,'' UHDIFV = '',E11.4,'' VZ0CM = '',E11.4,'' VZIUSTAR0 = '',E11.4 &
 & ,'' XNBMAX = '',E11.4,'' GALP = '',E11.4,'' SENSL = '',E11.4 &
 & ,'' TCA = '',E11.4,/,'' TCT = '',E11.4,'' TCW = '',E11.4 &
 & ,'' TURB = '',E11.4,'' TVF = '',E11.4,/&
 & ,'' TYM = '',5E11.4,'' QSMIN = '',E11.4,/&
 & ,'' QSSC = '',E11.4,'' RPHI0 = '',E11.4,'' RPHIR = '',E11.4 &
 & ,'' BEDIFV = '',E11.4,'' SCO = '',E11.4 &
 & ,'' USDMLT = '',E11.4,/&
 & ,'' GDDEVA = '',E11.4,'' GDDSDE = '',E11.4,/&
 & ,'' GWDCCO = '',E11.4 &
 & ,'' HUCOE2 = '',E11.4,/&
 & ,'' TENTR = '',E11.4,'' TENTRX = '',E11.4 &
 & ,'' TUDGP = '',E11.4,'' TDDGP = '',E11.4 &
 & ,'' GRCVPP = '',E11.4 &
 & )')&
 & AERCS1,AERCS3,AERCS5,ALMAV,ECMNP,EDB,EDC,EDD,EDK,ETACUT,EVAP,FONT,&
 & GWDAMP,GWDSE,GWDBC,GWDCD,GWDPROF,GWDVALI,&
 & HUCOE,HUTIL,HUTIL1,HUTIL2,VCHRNK,VKARMN,SNNBCO,SPNBCO,SXNBCO,HOBST,&
 & NPCLO1,NPCLO2,RCIN,RCVEVAP,REVGSL,RTCAPE,GCOMOD,&
 & UHDIFV,VZ0CM,VZIUSTAR0,XNBMAX,GALP,SENSL,TCA,TCT,TCW,&
 & TURB,TVF,TYM,QSMIN,QSSC,RPHI0,RPHIR,BEDIFV,SCO,USDMLT,&
 & GDDEVA,GDDSDE,GWDCCO,&
 & HUCOE2,&
 & TENTR,TENTRX,TUDGP,TDDGP,GRCVPP  
WRITE(UNIT=KULOUT,FMT='('' USUPRC = '',E11.4,'' USURIC = '',E11.4 &
 & ,'' QSNEBC = '',E11.4,'' QSNEBS = '',E11.4 &
 & ,'' QSSUSC = '',E11.4,'' QSSUSS = '',E11.4,'' QSSUSV = '',E11.4 &
 & ,'' QSUSXC = '',E11.4,'' QSUSXS = '',E11.4 &
 & ,'' GCCSV  = '',E11.4,/&
 & ,'' QXRAL  = '',E11.4,'' QXRDEL = '',E11.4 &
 & ,'' QXRHX  = '',E11.4,'' QXRR   = '',E11.4,'' QXRTGH = '',E11.4,/&
 & ,'' GWDLT  = '',E11.4,'' AHCLPV = '',E11.4,/&
 & ,'' GCVADS = '',E11.4,'' GCVBETA= '',E11.4 &
 & ,'' RICRLM = '',E11.4,'' XBLM   = '',E11.4,'' XKLM = '',E11.4 &
 & ,'' XMINLM = '',E11.4,'' XMAXLM = '',E11.4 &
 & )')&
 & USUPRC,USURIC,QSNEBC,QSNEBS,QSSUSC,QSSUSS,QSSUSV,QSUSXC,QSUSXS,GCCSV,QXRAL &
 & ,QXRDEL,QXRHX,QXRR,QXRTGH,GWDLT,AHCLPV,GCVADS,GCVBETA,RICRLM,XBLM,XKLM,XMINLM &
 & ,XMAXLM  
WRITE(UNIT=KULOUT,FMT='('' XWSALM = '',E11.4,'' XWSBLM = '',E11.4 &
 & ,'' GCVALFA= '',E11.4,'' GCVPSI = '',E11.4,'' GCVPSIE = '',E11.4 &
 & ,'' USURICL= '',E11.4,'' USURICE= '',E11.4,/ &
 & ,'' USURID= '',E11.4,'' USURIDE= '',E11.4,'' GCVNU= '',E11.4 &
 & ,'' GCVMLT= '',E11.4,'' GPBLHK0= '',E11.4,'' GPBLHRA = '',E11.4 &
 & ,'' UTILGUST= '',E11.4,'' RRGAMMA= '',E11.4,'' RRSCALE = '',E11.4 &
 & )')&
 & XWSALM,XWSBLM,GCVALFA,GCVPSI,GCVPSIE,USURICL,USURICE,USURID,USURIDE,GCVNU,GCVMLT,GPBLHK0 &
 & ,GPBLHRA,UTILGUST,RRGAMMA,RRSCALE  
WRITE(UNIT=KULOUT,FMT='('' RNEGAT == '',E11.4,'' RNLCURV = '',E11.4 &
 & )') RNEGAT,RNLCURV  
WRITE(UNIT=KULOUT,FMT='('' RDTFAC == '',E11.4 )') RDTFAC
WRITE(UNIT=KULOUT,FMT='(9(5(A,E11.4),/))') 'GCISMIN=',GCISMIN
WRITE(UNIT=KULOUT,FMT='('' A0ML_AU == '',E11.4,'' A0ML_AT == '',E11.4 &
 &,'' A0ML_BU == '',E11.4,'' A0ML_BT == '',E11.4 )') A0ML_AU,A0ML_AT&
 &,A0ML_BU,A0ML_BT

!     - - - - - - - -
!     For TKE scheme :
!     - - - - - - - -
! WRITE(UNIT=KULOUT,FMT='( '' ADISE   = '',E11.4,'' ADISI   = '',E11.4,/ &
!  &,'' AECLS3  = '',E11.4,'' AECLS4  = '',E11.4,'' AKN     = '',E11.4,/ &
!  &,'' ALD     = '',E11.4,'' ALPHAE  = '',E11.4,'' ALPHAT  = '',E11.4,/ &
!  &,'' ECTMIN  = '',E11.4,'' UCWSTAR = '',E11.4,'' UDECT   = '',E11.4,/ &
!  &,'' USHEARM = '',E11.4,/ &
!  &,'' UPRETMIN= '',E11.4,'' UPRETMAX= '',E11.4,'' ARSCH   = '',E11.4,/ &
!  &,'' ARSCQ   = '',E11.4,'' ARSC1   = '',E11.4,/ &
!  &,'' ARSB2   = '',E11.4,'' STTBMIN = '',E11.4,/ &
!  &,'' ACBRPHIM= '',E11.4,'' ALMAVE  = '',E11.4,'' RICRET  = '',E11.4,/ &
!  &,'' UETEPS  = '',E11.4,/ &
!  &,'' ABJUMIN = '',E11.4,'' RCOFLM  = '',E11.4,&
!      &)')&
   WRITE(UNIT=KULOUT,FMT=*) &
  &ADISE,ADISI,AECLS3,AECLS4,AKN,ALD,ALPHAE,ALPHAT,ECTMIN,UCWSTAR,UDECT,&
  &USHEARM,UPRETMIN,UPRETMAX,ARSCH,ARSCQ,ARSC1,ARSB2,STTBMIN,&
  &ACBRPHIM,ALMAVE,RICRET,UETEPS,AJBUMIN,RCOFLM


!       - - - - - - - - - -
!       Lopez Microphysics :
!       - - - - - - - - - -
WRITE(UNIT=KULOUT,FMT='('' - - - - - - - - - - - '')')
WRITE(UNIT=KULOUT,FMT='(''  Microphysics scheme  '')')
WRITE(UNIT=KULOUT,FMT='('' - - - - - - - - - - - '')')
WRITE(UNIT=KULOUT,FMT='( '' RAUTEFR = '',E11.4   &
 &,'' RAUTEFS = '',E11.4,'' RAUTSBET= '',E11.4,/ &
 &,'' RNINTR  = '',E11.4,'' RNINTS  = '',E11.4,/ &
 &,'' RQLCR   = '',E11.4,'' RQCRNS  = '',E11.4,/ &
 &,'' RQICRMIN= '',E11.4,'' RQICRMAX= '',E11.4,/ &
 &,'' RQLCV   = '',E11.4 &
 &,'' RQICVMIN= '',E11.4,'' RQICVMAX= '',E11.4,/ &
 &,'' RQICRT1 = '',E11.4,'' RQICRT2 = '',E11.4,/ &
 &,'' RQICRSN = '',E11.4,'' RACCEF  = '',F11.4,/ &
 &,'' RAGGEF  = '',F11.4,'' RRIMEF  = '',F11.4,/ &
 &,'' RHCRIT1 = '',F11.6,'' RHCRIT2 = '',F11.6,/ &
 &,'' RETAMIN = '',F11.6,'' RFACNSM = '',F11.6,/ &
 &,'' TFVR    = '',F11.6,'' TFVS    = '',F11.6,/ &
 &,'' GRHCMOD = '',F11.6,'' RHEVAP  = '',F11.6,/ &
 &    )') &
 &RAUTEFR,RAUTEFS,RAUTSBET,RNINTR,RNINTS,RQLCR,RQCRNS,RQICRMIN,RQICRMAX, &
 &RQLCV,RQICVMIN,RQICVMAX,RQICRT1,RQICRT2,RQICRSN,RACCEF,RAGGEF,RRIMEF, &
 &RHCRIT1,RHCRIT2,RETAMIN,RFACNSM,TFVR,TFVS,GRHCMOD,RHEVAP

!     - - - - - - - - - - -
!     For Grenier scheme :
!     - - - - - - - - - - -
!WRITE(UNIT=KULOUT,FMT='('' AGRE1   = '',E11.4,'' AGRE2   = '',E11.4,'' AGREF   = '',E11.4,/ &
! & ,'' AGRERICR= '',E11.4,'' AGREKE  = '',E11.4, &
! & )')&
WRITE(UNIT=KULOUT,FMT=*) &
 & AGRE1,    AGRE2,   AGREF,      AGRERICR, AGREKE

!     - - - - - - - - - - - - - - - -
!     For dry conv. adjustment scheme :                                -
!     - - - - - - - - - - - - - - - -
WRITE(UNIT=KULOUT,FMT='(&
 & '' AJ1PEPS = '',E11.4,'' AJ1MEPS = '',E11.4,'' NAJITER = '',I3,8X &
 & )')&
 & AJ1PEPS, AJ1MEPS, NAJITER

!     - - - - - - -
!     For ACCVIMPGY :
!     - - - - - - -
WRITE(UNIT=KULOUT,FMT='( '' ALFX    = '',E11.4  &
 &,'' TCTC    = '',E11.4,'' TVFC    = '',F11.6,/ &
 &,'' GAMAP1  = '',F11.6,'' RKDN    = '',E11.4,/ &
 &,'' VVN     = '',F11.6,'' VVX     = '',F11.6,/ &
 &,'' FENTRT  = '',F11.6,'' HCMIN   = '',F11.6  &
 &,'' FQLIC   = '',F11.6,'' FNEBC   = '',F11.6,/ &
 &,'' FEVAPC  = '',F11.6)')&
 & ALFX, TCTC, TVFC, GAMAP1, RKDN, VVN, VVX, FENTRT, &
 & HCMIN, FQLIC, FNEBC, FEVAPC


WRITE(UNIT=KULOUT,FMT='(''  Pseudo prognostic TKE scheme  '')')
WRITE(UNIT=KULOUT,FMT='( '' NUPTKE = '',E11.4  &
 &'' GAMTKE = '',E11.4)') NUPTKE,GAMTKE

WRITE(KULOUT,'(/'' PIL MICROPHYSICS : '')')
WRITE(UNIT=KULOUT,FMT='('' RAUIUSTE = '',E10.4  &
 &,'' RAUITN = '',F8.3,'' RAUITX = '',F8.3  &
 &,'' RDPHIC = '',F8.2   &
 &,'' GWBFAUT = '',F5.2 &
 &,'' RWBF1 = '',F6.2,'' RWBF2 = '',F5.2 &
 &,'' RSMDNEBX = '',F5.2 &
 &,'' RSMDTX = '',F5.2 &
 &,'' NSMTPA = '',I2 &
 &,'' NSMTPB = '',I2 &
 &,'' RCOLL = '',E11.4 &
 &,'' RFALLL = '',E11.4 &
 &    )') &
 &RAUIUSTE,RAUITN,RAUITX,RDPHIC,GWBFAUT,RWBF1,RWBF2,&
 &RSMDNEBX, RSMDTX, NSMTPA, NSMTPB, RCOLL, RFALLL
!----------------------------------------
WRITE(KULOUT,'(/'' PROGNOSTIC CONVECTION '')')
if (1==0) then     !!!!! A REVOIR (MPL)
WRITE(KULOUT,&
 &'('' TUDBU='',E16.6,'' TUDFR='',E16.6,'' TDDBU='',E16.6,'' TDDFR='', &
 & E16.6,&
 & '' GCVALMX='',G10.4,'' GCVACHI='',E16.6,&
 & '' GCVADMW='',G10.4,&
 & '' GCVEEX ='',G10.4,&
 & '' GCVBEE ='',G10.4,&
 & '' GCVSQDN ='',G10.4,&
 & '' GCVSQDR ='',G10.4,&
 & '' GCVSQDCX ='',G10.4,&
 & '' ECMNPI ='',G10.4,&
 & '' GFRIC ='',G10.4&
 & )')TUDBU,TUDFR,TDDBU,TDDFR,GCVALMX,GCVACHI,&
 &    GCVADMW,GCVEEX,GCVBEE, &
 &    GCVSQDN, GCVSQDR, GCVSQDCX, ECMNPI,GFRIC
else
  print*,'>>>>> TOTO1 ', TUDBU,TUDFR,TDDBU,TDDFR,GCVALMX,GCVACHI,&
 &    GCVADMW,GCVEEX,GCVBEE, &
 &    GCVSQDN, GCVSQDR, GCVSQDCX, ECMNPI,GFRIC
endif
WRITE(KULOUT,&
 &'('' GDDEVF='',G10.4,'' GDDWPF='',G10.4,'' GDDBETA='',G10.4)')&
 & GDDEVF, GDDWPF, GDDBETA
WRITE(KULOUT,&
 & '('' GRRINTE ='',G10.4,'' GRRMINA ='',G10.4&
 & )') GRRINTE, GRRMINA
!     ------------------------------------------------------------------

!*       4.    Check consistency between logical and real namelist inputs.
!              -------------------

IF (LMPHYS) THEN

  IF(LCVRA.AND.((TUDGP /= 0.0_JPRB.OR.TDDGP /= 0.0_JPRB).AND..NOT.LCVCAS)) THEN
    WRITE(UNIT=KULOUT,FMT='(A)') &
     & 'INCONSISTENCY BETWEEN TUDGP, TDDGP AND LCVCAS!...'  
    CALL ABOR1('TUDP<>0. OR TDDGP<>0. IMPLIES LCVCAS=T!...')
  ENDIF

  IF(.NOT.LCVLIS.AND.GCVPSI /= 0.0_JPRB) THEN
    WRITE(UNIT=KULOUT,FMT='(A)') 'INCONSISTENCY BETWEEN LCVLIS AND GCVPSI!...'
    CALL ABOR1('LCVLIS=F IMPLIES GCVPSI=0.!...')
  ENDIF

ENDIF

!*  Consistency check for simplified physics keys

IF ((LMPHYS.OR.LSIMPH).AND.(LSMOOTHD.OR.LSMOOTHA)) THEN

  IF(RNLCURV == 0.0_JPRB) THEN
    WRITE(UNIT=KULOUT,FMT='(A)') 'SMOOTHING IN SIM. PH. BUT RNLCURV = ZERO'
    CALL ABOR1('PHYSICS AND SMOOTHING IMPLY RNLCURV /= ZERO')
  ENDIF

ENDIF
IF (LPIL) THEN
! Setup RHCRI profile(s)
! (for LAM case... GAW set up by suecuv, called before suphy)
! (For Not LAM, SULEG called by sugem called before suphy)
  CALL SURHCRI(KULOUT)
! CALL SULOCST(KULOUT)
ENDIF

IF (LHOOK) CALL DR_HOOK('SUPHY0',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHY0
