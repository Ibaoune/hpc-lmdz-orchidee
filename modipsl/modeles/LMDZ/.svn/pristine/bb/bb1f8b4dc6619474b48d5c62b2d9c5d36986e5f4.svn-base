SUBROUTINE SU0PHY(KULOUT)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!USE YOMLUN   , ONLY : NULOUT   ,NULNAM  MPL 14.04.09
USE YOMLUN   , ONLY : NULOUT
! Ce qui concerne LSLPHY commente par MPL 24.11.08
!USE YOMSLPHY , ONLY : LSLPHY
USE YOMPHY   , ONLY : NBITER   ,NOIR, NDPSFI, NPHYREP ,LMPHYS   ,LREASUR  ,&
 & LCAPE    ,LCONDWT  ,LCVPP    ,LCVDD    ,LHUNEG   ,&
 & LNEIGE   ,LRNUMX   ,LCLSATUR,L2PHYS    ,LCVRA    ,LGWD     ,&
 & LGWDC    ,LHMTO    ,LNEBCO   ,LNEBN    ,LNEBR    ,LQXRTGH, LHUCN,&
 & LNEBT    ,LND2DIFF  ,LOZONE   ,LRAY     ,LRAYFM   ,LRAYFM15 ,&
 & LRRMES   ,LSFHYD   ,LSNV     ,LSOLV    ,LFGEL    ,&
 & LSRCON   ,LSRCONT  ,LSLC     ,LRRGUST  ,LRELAXW  ,&
 & LAEROSEA ,LAEROLAN ,LAEROSOO ,LAERODES ,LAEROVOL ,LAEROSUL ,LRELAXT  ,&
 & LO3ABC   ,LSTRA    ,LSTRAS   ,LTHERMO  ,LVDIF    ,&
 & LRAYLU   ,LREWS    ,LRPROX   ,LRMIX    ,LRSTAB   ,&
 & LRAUTOEV ,LRAYPL   ,LCVLIS   ,LCVCAS   ,LVGSN    ,&
 & LNEBNXR  ,LFPCOR   ,LNOIAS   ,CGMIXLEN ,LPRGML   ,LGLT     ,&
 & LNEWD    ,LRTPP    ,LRTDL    ,LDIFCONS ,LECT     ,&
 & LCVPGY   ,LPROCLD  ,LEVAPP   ,LCOLLEC  ,LPTKE    ,L3MT     ,&
 & LCVPRO   ,LCDDPRO  ,LSCMF    ,LVOIGT   ,LVFULL   ,&
 & LNSMLIS, LPHCDPI   ,&
 & NPHY     ,JPHYEC   ,JPHYMF   ,JPHYARO  , & 
 & LAJUCV   ,LPBLE    ,LNEBGR   ,LNEBGY   ,&
 & LBCCOND  ,LCVRAV3  ,LZ0HSREL ,LBLVAR   ,&
 & LADJCLD  ,LAUTONEB ,LSSD     ,LCVPPKF  ,LECTFL ,&
 & LPIL     ,LPHSPSH  ,LSMROT   ,LSMNIMBT ,&
 & LSMTPS   ,L1DRHCRI ,LGWRHCRI ,NSMTBOT   , &
 & NSMDNEB  ,NPRAG    ,NPRAC    ,NPRRI    ,LSTRAPRO, LNEWSTAT
USE YOMARPHY , ONLY : LMPA    ,LMICRO   ,LTURB    ,&
 & LMSE  ,LKFBCONV ,LKFBD    ,LKFBS    ,LUSECHEM ,&
 & LORILAM  ,LRDUST, LBUFLUX ,CCOUPLING 
USE YOEPHY   , ONLY : LEPHYS    ,&
 & LECOND   ,LECUMF   ,LEDCLD   ,LEEVAP   ,LEGWDG   ,&
 & LEOZOC   ,LEQNGT   ,LERADI   ,LERADS   ,&
 & LESHCV   ,LESICE   ,LESURF   ,LEVDIF   ,&
 & LAGPHY   ,LEPCLD   ,LECO2DIU ,&
 & LEO3CH   ,LBUD23   ,LEMETHOX ,LERA40   ,LECURR   ,LVDFTRAC ,&
 & LEOCWA   ,LEOCCO   ,LEOCSA   ,LMFTRAC  ,LERAIN   ,LE4ALB   ,&
 & RTHRFRTI ,NEPHYS_PCFULL, LEMWAVE
USE YOEWCOU  , ONLY : NSTPW    ,RSOUTW   ,RNORTW   ,RDEGREW  ,&
 & LWCOU    ,LWCOU2W
! Ce qui concerne RCLDTOPP commente par MPL 24.11.08
!USE YOECLDP  , ONLY : RCLDTOPP  
USE YOPHLC   , ONLY : ALPHA    ,AH0      ,USTARL   ,&
 & USTARS   ,ALANDZ0  ,ASEAZ0   ,LSPHLC   ,LVDFLC   ,&
 & LSDRLC   ,LCZDEB   ,LZMCON   ,LKEXP    ,LVDFDS   ,&
 & LSDRDS  
USE YOMCT0   , ONLY : NCONF    ,&
 & LRETCFOU ,LWRTCFOU, LAROME,LFPOS,LPC_FULL
USE YOMCT0B  , ONLY : LECMWF
! Tous les YOM* ci-dessous commentes par MPL 24.11.08
!USE YOMINI   , ONLY : NEINI
!USE YOMDFI   , ONLY : NEDFI
!USE YOMVRTL  , ONLY : L131TL
!USE YOPHNC   , ONLY : LETRAJP  ,LETRAJPT ,LERADI2  ,LERADS2   ,&
! & LERADSW2 ,LERADN2 ,LERADFL2 ,LEDCLD2  ,LENCLD2   ,&
! & LEVAPLS2 ,LEVDIF2 ,LEGWDG2  ,LECUMF2  ,LECUBM2   ,&
! & LECOND2  ,LEQNGT2 ,LESURF2  ,LEKPERT  ,LTRACLNPH
!USE YOEPHLI  , ONLY : LENOPERT
!USE YOMNCL   , ONLY : LNCLIN   ,LREGCL
USE YOMSIMPHL, ONLY : LSIMPH   ,LTRAJPS  ,LTRAJPST ,&
 & LSMOOTHD ,LSMOOTHA ,LSMOOTHB ,LCVRASP  ,&
 & LGWDSP   ,LRAYSP   ,LSTRASP  ,LVDIFSP  ,LRRMESSP ,LCLOUDS  
USE YOMRCOEF , ONLY : LRCOEF    ,LTLADDIA ,LGLOBRAD
!USE YOMIOP   , ONLY : NPCKFT95 ,NEXPBT95
!USE YOMDYNA  , ONLY : LGWADV
!USE YOMCOAPHY   , ONLY : NPHYINT

!**** *SU0PHY*   - Initialize common YOxPHY controlling physics

!     Purpose.
!     --------
!           Initialize YOxPHY, the common that includes the
!           basic switches for the physics of the model.

!**   Interface.
!     ----------
!        *CALL* *SU0PHY(KULOUT) from SU0YOMA

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!        or

!        Documentation ARPEGE (depending on which physics will be used)

!     Author.
!     -------
!        J.-J. Morcrette                    *ECMWF*
!        J.-F. Geleyn for the ARPEGE rewriting.

!     Modifications.
!     --------------
!        Original : 91-11-12
!        Modified 92-02-22 by M. Deque (tests of consistency with *YOMDPHY*)
!        Modified by R. EL Khatib : 93-04-02 Set-up defaults controled by LECMWF
!        Modified 94-02-28 by M.  Deque  : Shallow convection clouds
!        Modified 93-10-28 by Ph. Dandin : FMR scheme with MF physics
!        Modified 93-08-24 by D. Giard (test of consistency with digital filter)
!        Modified by M. Hamrud    : 93-06-05 Make use of LECMWF for ECMWF
!        Modified 95-11-27 by M. Deque (2nd call to APLPAR)
!        Modified 96-01-10 by M. Janiskova (logical switches for simpl.ph.par.)
!        Modified 97-02-28 by J.M. Piriou (cloudiness scheme switch LNEBN)
!        Modified 97-04-17 by J.M. Piriou (default values)
!        Modified by F. Rabier    : 96-09-25 Full physics set-up for 801 job
!        Modified by G. Hello     : 97-07-31 Full MF physics set-up for 801
!        Modified by M. Deque     : 97-05-25 Frozen FMR
!        Modified by E.Bazile     : 97-11-18 Soil freezing (LFGEL)
!        Modified by M. Deque     : 98-01-05 Cleaning of NDPSFI
!        Modified by M. Janiskova : 98-11-18 Full MF physics set-up for 131
!                                 : 99-02-21 Set-up for radiation coef.
!        Modified by C. Jakob     : 98-04    Methane oxidation
!        Modified by T. Bergot    : 98-08  Full MF physics set-up for 601
!        Modified by E.Bazile     : 99-02-12 Superficial soil freezing (LFGELS)
!        Modified by L. Gerard    : 98-12-07 LSRCON
!        Modified by J.M. Piriou  : 99-04-19 Moon radiation
!        Modified 99-02-17 by K. YESSAD: options LRETCFOU, LWRTCFOU.
!        Modified by J.M. Piriou  : 99-06-18 LCVLIS and LCVCAS
!        Modified by J.M. Piriou  : 99-07-07 Introduce reproductibility in physics (NPHYREP).
!        Modified by D. Giard     : 2000-10-25 LVGSN (snow and vegetation)
!        Modified by E. Bazile    : 2000-11-12 CYCORA's default value.      
!        Modified by F. Bouyssel  : 2001-03-03 LRRMESSP     
!        Modified by J.M. Piriou  : 2001-11-27 LSRCONT
!        Modified by J.M. Piriou  : 2002-11-15 LNEBNXR
!        Modified by F. Bouyssel  : 2002-06-25 LCVPPLIS and LRRGUST
!        R. El Khatib : 2001-08-07 Pruning options
!        J.M. Piriou  : 2002-01-10 set default values to operational ones.
!        Modified by A.Beljaars   : 2002-11-12 LECURR (Ocean current)
!        Modified by Y. Bouteloup : 2002-03-05 LO3ABC
!        Modified 08-2002 C. Smith : use "w" as prognostic variable in the
!         semi-lag advection of vertical divergence in the NH model. 
!        Modified by F. Bouyssel  : 2002-12-18 Cleaning of LCVPPLIS
!        Modified by E. Bazile      : 2003-02-13 LFPCOR
!        Modified by E. Bazile      : 2003-02-18 LNOIAS
!        Modified by M. Janiskova : 2003-05 set-up for ECMWF stat.cloud scheme
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P. Viterbo   ECMWF   03-12-2004  Include user-defined RTHRFRTI
!        Modified by Y Seity      : 2004-11-16 For AROME setup, default values 
!            for namarphy keys,read namarphy, switch off arp/ald physics keys
!        Modified by R. Brozkova : 2004-11 modifs for Xu-Randall cloud. scheme
!        P. Marquet and F. Bouyssel : 2004-08-18 (Lopez) 
!        G. Hello : 2005-04-25 Lground initialization (surfex and arome)
!        Y. Seity : 2005-09-25 LRDUST, LORILAM and LUSECHEM initialisation
!                 for AROME
!        T. Kovacic : 2006-03-17 LPHCDPI, NPHY, LBUFLUX; for DDH and BUDGET
!        Modified by F. Bouyssel  : 2005-01-25 Change default of LNSMLIS
!        D.Salmond     22-Nov-2005 Mods for coarser/finer physics
!        P. Lopez      14-02-2006  Added switch LTRACLNPH for including  
!                                  tracers in linearized physics
!        R.Brozkova : pre ALARO 0 modset: mixing lengths computation
!        Modified by GMGEC/EAC   : 2006-03 list of modif.
!                 P. Marquet  : 99-01-18 Dry Conv. Adj.     (LAJUCV)
!                 P. Marquet  : 02-02-14 YOMPHY
!                 P. Marquet  : 02-06-18 LBCCOND (Becht/Chab ACCOND)
!                 P. Marquet  : 02-08-30 LCVRAV3    (old ACCVIMP_V3)
!                     and new cloud model : 2006-03-03
!        Modified by E. Bazile    : 2006-04-11 Add LPBLE in case of LECT
!        Modified by E. Bazile    : 2006-04-20 Add LADJCLD,LCVPPKF, LECTFL
!                 JF. Gueremy : LZ0HSREL initialised (used in ACHMT)
!        Modified by F. Bouyssel  : 2006-10-30 Add LAUTONEB, LSSD
!        Modified by F. Vana     : 2006-01-30 LPTKE
!        Modifed by D. Banciu    : 2006-08-31 LND2DIFF
!        M. Bellus : 28-Sep-2006 ALARO-0 phasing: L3MT,LCVPRO,LCDDPRO,LSCMF,
!                    LVOIGT,LVFULL,LPIL,LPHSPSH,LSMROT,LSMNIMBT,LSMTPS,
!                    L1DRHCRI,LGWRHCRI,NSMTBOT,NSMDNEB,NPRAG,NPRAC,NPRRI,
!                    LSTRAPRO,LNEWSTAT ; removed duplicity in LPBLE
!        JJMorcrette 20060721 MODIS albedo
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IERR
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "suscm.intfb.h"

!     ------------------------------------------------------------------

#include "namphy.h"
!#include "namarphy.h"
#include "naephy.h"
#include "naphlc.h"
!#include "namtrajp.h"
!#include "namsimphl.h"
!#include "namrcoef.h"

!     ------------------------------------------------------------------
print*,'Dans SUOPHY ', KULOUT

IF (LHOOK) CALL DR_HOOK('SU0PHY',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

!*       1.1.1 Set default values for Meteo-France physics
!              -------------------------------------------

LMPHYS=.FALSE.
LREASUR=.TRUE.

CGMIXLEN='Z'
LPRGML=.FALSE.
LCAPE=.FALSE.
LCONDWT=.FALSE.
LCVCAS=.FALSE.
LCVLIS=.FALSE.
LCVPP=.FALSE.
LCVPPKF=.FALSE.
LCVDD=.FALSE.
LHUNEG=.TRUE.
LNEIGE=.TRUE.
LRNUMX=.FALSE.
LCLSATUR=.FALSE.
LVOIGT=.FALSE.
LVFULL=.FALSE.
LRRGUST=.FALSE.
L2PHYS=.FALSE.
LO3ABC=.FALSE.
LAEROSEA=.FALSE.
LAEROLAN=.FALSE.
LAEROSOO=.FALSE.
LAERODES=.FALSE.
LAEROVOL=.FALSE.
LAEROSUL=.FALSE.
LRELAXT=.FALSE.
LRELAXW=.FALSE.
LE4ALB=.FALSE.
LGLT=.FALSE.
LNEWD=.FALSE.
LDIFCONS=.FALSE.
LECT=.FALSE.
LPTKE=.FALSE.
LECTFL=.FALSE.
LPBLE=.FALSE.
LCVPGY=.FALSE.
L3MT=.FALSE.

LCVRA=.FALSE.
LGWDC=.FALSE.
LGWD=.FALSE.
LHMTO=.FALSE.
LNEBCO=.FALSE.
LNEBN=.FALSE.
LNEBNXR=.FALSE.
LNEBR=.FALSE.
LNEBT=.FALSE.
LND2DIFF=.FALSE.
LQXRTGH=.FALSE.
LHUCN=.FALSE.
LOZONE=.FALSE.
LRAY=.FALSE.
LRAYLU=.FALSE.
LREWS=.FALSE.
LRPROX=.FALSE.
LRMIX=.FALSE.
LRSTAB=.FALSE.
LRAUTOEV=.FALSE.
LRTPP=.FALSE.
LRTDL=.FALSE.
LRAYPL=.FALSE.
LRAYFM=.FALSE.
LRAYFM15=.FALSE.
LRRMES=.FALSE.
LSFHYD=.FALSE.
LSNV=.FALSE.
LSOLV=.FALSE.
LFGEL=.FALSE.
LSRCON=.FALSE.
LSRCONT=.FALSE.
LSLC=.FALSE.
LSTRA=.FALSE.
LSTRAS=.FALSE.
LTHERMO=.FALSE.
LVDIF=.FALSE.
LVGSN=.FALSE.
LFPCOR=.FALSE.
LNOIAS=.FALSE.
LPHCDPI=.FALSE.
LBLVAR=.FALSE.

LZ0HSREL=.FALSE.

NBITER=2
NDPSFI=0
NPHYREP=1
NOIR=0

! ---------------------------------------------------
! ALARO-0 (cloud)
! ---------------------------------------------------
LPIL=.FALSE.
LSTRAPRO=.FALSE.
LNEWSTAT=.TRUE.
LPHSPSH=.FALSE. ! Pseudo Historic Surface Precip Sensible Heat Flux
LSMROT=.FALSE.
LSMTPS=.FALSE.
LSMNIMBT=.FALSE.
L1DRHCRI=.FALSE.
LGWRHCRI=.FALSE.
NSMTBOT=0 ! interpolate
NSMDNEB=2 ! gradient limitation
NPRAG=1
NPRAC=1
NPRRI=1
! ---------------------------------------------------
! ALARO-0 (prognostic convection)
! ---------------------------------------------------
LCVPRO=.FALSE.
LCDDPRO=.FALSE. 
LSCMF=.FALSE.

! - - - - - - - - - - - - - - - - - - - - - - - - - -
! Cloud and precipitation prognostic scheme (Lopez) :
! - - - - - - - - - - - - - - - - - - - - - - - - - -

LPROCLD=.FALSE.
LEVAPP=.TRUE.
LCOLLEC=.TRUE.
LNSMLIS=.TRUE.
LADJCLD=.TRUE.
LAUTONEB=.FALSE.
LSSD=.FALSE.

!AROME physics
LMPA=.FALSE.
LMICRO=.FALSE.
LTURB=.FALSE.
LMSE=.FALSE.
LKFBCONV=.FALSE.
LKFBD=.FALSE.
LKFBS=.FALSE.
LUSECHEM=.FALSE.
LORILAM=.FALSE.
LRDUST=.FALSE.
LBUFLUX=.TRUE.
CCOUPLING='E'

! - - - - - - - - - - -
! Module YOMPHY :
! - - - - - - - - - - -

LAJUCV=.FALSE.
LNEBGR=.FALSE.
LNEBGY=.FALSE.
LBCCOND=.FALSE.

LCVRAV3=.FALSE.

!*    1.1.2  Set default values for simplified physical parametrization
!            of Meteo-France
!     -----------------------------------------------------------------

LSIMPH=.FALSE.
LTRAJPS=.FALSE.
LTRAJPST=.FALSE.
LSMOOTHD=.FALSE.
LSMOOTHA=.FALSE.
LSMOOTHB=.FALSE.
LCLOUDS=.FALSE.

LCVRASP=.FALSE.
LGWDSP=.FALSE.
LRAYSP=.FALSE.
LSTRASP=.FALSE.
LVDIFSP=.FALSE.
LRRMESSP=.FALSE.

LRCOEF=.FALSE.
LTLADDIA=.FALSE.
LGLOBRAD=.FALSE.

!*    1.2.1  Set default values for ECMWF physics package
!            --------------------------------------------

LEPHYS=.FALSE.
LAGPHY=.TRUE.

LECOND=.FALSE.
LEPCLD=.FALSE.
LECUMF=.FALSE.
LEDCLD=.FALSE.
LEEVAP=.TRUE.
LEGWDG=.FALSE.
LEOZOC=.FALSE.
LEQNGT=.FALSE.
LERADI=.FALSE.
LERADS=.FALSE.
LESHCV=.FALSE.
LESICE=.TRUE.
LESURF=.FALSE.
LEVDIF=.FALSE.
LEOCWA=.FALSE.
LEOCCO=.FALSE.
LEOCSA=.FALSE.
LEMETHOX=.FALSE.
LERA40=.FALSE.
LECURR=.FALSE.
LVDFTRAC=.TRUE.
LMFTRAC=.TRUE.
LERAIN=.FALSE.
LE4ALB=.FALSE.
RTHRFRTI=0.0_JPRB
!NPHYINT=0
LE4ALB=.FALSE.

!-------------------------------------------------------
! pressure above which cloud scheme is not called
! !!!WARNING!!! this has to be in the part of the domain 
! where the pure pressure level grid is used, otherwise
! the code is not bit-reproducible!!!
! Don't call the cloud scheme for pressures lower than 1hPa
!RCLDTOPP=100.0_JPRB

!--------------------------------------------------------

!*     1.2.2  Set-up linearized physical parametrization of ECMWF
!             ---------------------------------------------------

!LETRAJP = .FALSE.
!LETRAJPT= .FALSE.
!LERADI2 = .FALSE.
!LERADS2 = .FALSE.
!LERADSW2= .FALSE.
!LERADN2 = .FALSE.
!LERADFL2= .FALSE.
!LEDCLD2 = .FALSE.
!LENCLD2 = .FALSE.
!LEVAPLS2= .FALSE.
!LEVDIF2 = .FALSE.
!LEGWDG2 = .FALSE.
!LECUMF2 = .FALSE.
!LECUBM2 = .FALSE.
!LECOND2 = .FALSE.
!LEQNGT2 = .FALSE.
!LESURF2 = .FALSE.
!LEKPERT = .FALSE.
!LNCLIN  = .FALSE.
!LREGCL  = .FALSE.
!LTRACLNPH = .FALSE.

! No perturbation of surface arrays
!LENOPERT = .TRUE.

!*       Packing parameters
!        -------------------
!NPCKFT95 = 1
!NEXPBT95 = 6

!     LOGICAL FOR THE VERT DIFF SCHEME VDIFLCZ USED IN CONF 601

LSPHLC  = .FALSE.
LVDFLC  = .FALSE.
LVDFDS  = .TRUE.
LSDRLC  = .TRUE.
LSDRDS  = .FALSE.
LCZDEB  = .FALSE.
LZMCON  = .TRUE.
LKEXP   = .TRUE.
ALPHA   = 3._JPRB
AH0     = 1000.0_JPRB
USTARL  = 0.5_JPRB
USTARS  = 0.2_JPRB
ALANDZ0 = 0.05_JPRB
ASEAZ0  = 0.0005_JPRB

NSTPW=2
!     NSTPW=30

LWCOU=.FALSE.
LWCOU2W=.FALSE.

!  Setup for 3 degree resolution wave model

!     RNORTW= 72.0
!     RSOUTW=-63.0
!     RDEGREW=3.0

!  Setup for 1.5 degree resolution wave model

RNORTW= 81.0_JPRB
RSOUTW=-81.0_JPRB
RDEGREW=1.5_JPRB

!        1.3 Modify default values according to LECMWF

IF (LECMWF) THEN
  LEPCLD=.TRUE.
  LEVDIF=.TRUE.
  LEOCWA=.FALSE.
  LEOCCO=.FALSE.
  LEOCSA=.TRUE.
  LESURF=.TRUE.
  LECOND=.FALSE.
  LECUMF=.TRUE.
  LEEVAP=.TRUE.
  LEGWDG=.TRUE.
  LEOZOC=.TRUE.
  LEQNGT=.TRUE.
  LERADI=.TRUE.
  LERADS=.TRUE.
  LESICE=.TRUE.
  LEDCLD=.TRUE.
  LEO3CH=.FALSE.
  LECO2DIU=.FALSE.
  LEMETHOX=.TRUE.
  IF(NCONF == 1) THEN
    LEPHYS=.TRUE.
  ENDIF
  IF(NCONF == 131) THEN
!   LERADI2 = .TRUE.
!   LERADS2 = .TRUE.
!   LEVDIF2 = .TRUE.
!   LEGWDG2 = .TRUE.
!   LECUMF2 = .TRUE.
!   LECOND2 = .TRUE.
!    IF(L131TL) THEN
!      LEPHYS=.TRUE.
!      LECOND=.TRUE.
!      LEPCLD=.FALSE.
!    ENDIF
    LSPHLC=.TRUE.
    LVDFLC=.TRUE.
    LSDRLC=.TRUE.
    LZMCON=.TRUE.
    LKEXP =.TRUE.
  ENDIF
  IF (NCONF == 401 .OR. NCONF == 501) THEN
    LEPHYS = .TRUE.
    LECOND = .TRUE.
    LEPCLD = .FALSE.

!   LETRAJP = .TRUE.
!   LERADI2 = .TRUE.
!   LERADS2 = .TRUE.
!   LEVDIF2 = .TRUE.
!   LEGWDG2 = .TRUE.
!   LECUMF2 = .TRUE.
!   LECOND2 = .TRUE.
  ENDIF

  IF(NCONF == 601.OR.NCONF == 801) THEN

!      Full physics when computing the trajectory

    LEPHYS=.TRUE.

!      Simple scheme for TL and ADJ

    LSPHLC  = .TRUE.
    LVDFLC  = .TRUE.
    LVDFDS  = .TRUE.
    LSDRLC  = .TRUE.
    LSDRDS  = .TRUE.
    LCZDEB  = .FALSE.
    LZMCON  = .TRUE.
    LKEXP   = .TRUE.
    ALPHA   = 3._JPRB
    AH0     = 1000.0_JPRB
    USTARL  = 0.5_JPRB
    USTARS  = 0.2_JPRB
    ALANDZ0 = 0.05_JPRB
    ASEAZ0  = 0.0005_JPRB

  ENDIF
  NPHY = JPHYEC
ELSE
  IF ((NCONF == 1).OR.(NCONF == 131).OR.(NCONF == 401)&
     & .OR.(NCONF == 501).OR.(NCONF == 601)&
     & .OR.(NCONF == 701).OR.(NCONF == 801)) THEN  
    LAGPHY=.FALSE.
    LCVCAS=.TRUE.
    LCVDD=.TRUE.
    LCVLIS=.TRUE.
    LCVPP=.TRUE.
    LSRCONT=.FALSE.
    LCVRA=.TRUE.
    LFGEL=.TRUE.
    LGWD=.TRUE.
    LHMTO=.TRUE.
    LMPHYS=.TRUE.
    LNEBN=.TRUE.
    LNEBT=.FALSE.
    LNEBNXR=.FALSE.
    LNEIGE=.TRUE.
    LRAY=.TRUE.
    LRAYLU=.TRUE.
    LRAYPL=.TRUE.
    LRRMES=.TRUE.
    LSFHYD=.TRUE.
    LSOLV=.TRUE.
    LSRCON=.TRUE.
    LSTRA=.TRUE.
    LTHERMO=.TRUE.
    LVDIF=.TRUE.
    NPHY = JPHYMF
    IF(LAROME.AND..NOT.LFPOS) THEN
      !extinction of MF phsic's keys for arome's run
      !lfpos is added for 927 confs for example.
      LCVRA=.FALSE.
      LGWD=.FALSE.
      LGWDC=.FALSE.
      LNEBCO=.FALSE.
      LNEBN=.FALSE.
      LNEBR=.FALSE.
      LNEBT=.FALSE.
      LOZONE=.FALSE.
      LRAY=.FALSE.
      LRAYLU=.FALSE.
      LREWS=.FALSE.
      LRAYPL=.FALSE.
      LRAYFM=.FALSE.
      LRAYFM15=.FALSE.
      LRRMES=.FALSE.
      LSFHYD=.FALSE.
      LSNV=.FALSE.
      LSOLV=.FALSE.
      LFGEL=.FALSE.
      LSTRA=.FALSE.
      LSTRAS=.FALSE.
      LVDIF=.FALSE.
      !initialisation of AROME's physic ones
      LMPA=.TRUE.
      LMICRO=.FALSE.
      LTURB=.FALSE.
      LMSE=.FALSE.
      LKFBCONV=.FALSE.
      LUSECHEM=.FALSE.
      LKFBD=.FALSE.
      LKFBS=.FALSE.
      LUSECHEM=.FALSE.
      LORILAM=.FALSE.
      LRDUST=.FALSE.
      NPHY = JPHYARO 
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

! Ce qui concerne NAMPHY commente par MPL le 14.04.09
!CALL POSNAM(NULNAM,'NAMPHY')
!READ(NULNAM,NAMPHY)

! Ce qui concerne NAMARPHY commente par MPL le 24.11.08
!CALL POSNAM(NULNAM,'NAMARPHY')
!READ(NULNAM,NAMARPHY)
!IF(LMSE) THEN
!  LSOLV=.FALSE.
!  LFGEL=.FALSE.
!  LSNV=.FALSE.
!  LVGSN=.FALSE.
!  WRITE(NULOUT,'('' INFO -old isba-lsolv,lfgel,lsnv,lvgsn- reset to .F. '', &
!                 &'' when using lground i.e. surfex '')')
!ENDIF

NEPHYS_PCFULL=3

! Ce qui concerne NAEPHY commente par MPL le 14.04.09
!CALL POSNAM(NULNAM,'NAEPHY')
!READ(NULNAM,NAEPHY)

IF(LPC_FULL)THEN
  IF(NEPHYS_PCFULL < 2 .OR.NEPHYS_PCFULL > 3)THEN
    CALL ABOR1(' SU0PHY: NEPHYS_PCFULL')
  ENDIF
ELSE
  NEPHYS_PCFULL=3
ENDIF

IF(     NCONF == 201.OR.NCONF == 202 &
   & .OR.NCONF == 421.OR.NCONF == 422 &
   & .OR.NCONF == 521.OR.NCONF == 522  )THEN  
  LREASUR=.FALSE.
ENDIF
!IF(.NOT.(LEPHYS.OR.LMPHYS)) THEN
!  LAGPHY=.FALSE.
!  LSLPHY=.FALSE.
!  WRITE(NULOUT,'('' INFO - LSLPHY RESET TO .FALSE. '')')
!ENDIF
IF(LEPCLD) THEN
  LECOND=.FALSE.
ENDIF
IF(LAGPHY.AND.NDPSFI /= 0) THEN
  NDPSFI=0
  WRITE(NULOUT,'('' WARNING - NDPSFI RESET TO 0 '')')
ENDIF

! Commente par MPL 24.11.08
!CALL POSNAM(NULNAM,'NAPHLC')
!READ(NULNAM,NAPHLC)

! Commente par MPL 24.11.08
!CALL POSNAM (NULNAM,'NAMTRAJP')
!READ (NULNAM,NAMTRAJP)

! Commente par MPL 24.11.08
!CALL POSNAM(NULNAM,'NAMSIMPHL')
!READ(NULNAM,NAMSIMPHL)

! Commente par MPL 24.11.08
!CALL POSNAM(NULNAM,'NAMRCOEF')
!READ(NULNAM,NAMRCOEF)

!-------------------------------------------------
! Initialize profile extractions for the Single Column Model .
!-------------------------------------------------

CALL SUSCM(KULOUT)

!     ------------------------------------------------------------------

!*       3.    Do tests of consistency.
!              ------------------------

IERR=0

! * Test that adequate physics is activated when at least one of the
!   options LRETCFOU, LWRTCFOU, LRCOEF,is set to .TRUE.
IF( LWRTCFOU .AND. ((.NOT.LMPHYS).OR.LEPHYS.OR.LAGPHY) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SU0PHY: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,'(1X,A,A,A)') ' SU0PHY: options ',&
   & 'LWRTC, LWRTCFOU ',&
   & 'require the following options for physics:'  
  WRITE(KULOUT,'(1X,A)') '         - LMPHYS=.TRUE. '
  WRITE(KULOUT,'(1X,A)') '         - LEPHYS=.FALSE. '
  WRITE(KULOUT,'(1X,A)') '         - LAGPHY=.FALSE. '
  IF(.NOT.LRAYSP) THEN
    WRITE(KULOUT,'(1X,A,A,A)') ' SU0PHY: option ',&
     & ' LWRTCFOU ',&
     & 'require the following option for simplified physics:'  
    WRITE(KULOUT,'(1X,A)') '         - LRAYSP=.TRUE. '
  ENDIF
ENDIF
IF( LRETCFOU .AND. ((.NOT.(LMPHYS.OR.LSIMPH)).OR.LEPHYS.OR.LAGPHY) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SU0PHY: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,'(1X,A,A,A)') ' SU0PHY: option ',&
   & ' LRETCFOU ',&
   & 'require the following options for physics:'  
  WRITE(KULOUT,'(1X,A)') '         - LMPHYS=.TRUE. or LSIMPH=.TRUE.'
  WRITE(KULOUT,'(1X,A)') '         - LEPHYS=.FALSE. '
  WRITE(KULOUT,'(1X,A)') '         - LAGPHY=.FALSE. '
ENDIF
IF( LRCOEF .AND. ((.NOT.LMPHYS).OR.LEPHYS.OR.LAGPHY) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SU0PHY: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,'(1X,A,A,A)') ' SU0PHY: options ',&
   & 'LRCOEF ',&
   & 'require the following options for physics:'  
  WRITE(KULOUT,'(1X,A)') '         - LMPHYS=.TRUE. '
  WRITE(KULOUT,'(1X,A)') '         - LEPHYS=.FALSE. '
  WRITE(KULOUT,'(1X,A)') '         - LAGPHY=.FALSE. '
  IF(.NOT.LRAYSP) THEN
    WRITE(KULOUT,'(1X,A,A,A)') ' SU0PHY: options ',&
     & 'LRCOEF ',&
     & 'require the following option for simplified physics:'  
    WRITE(KULOUT,'(1X,A)') '         - LRAYSP=.TRUE. '
  ENDIF
ENDIF

! * Test that options LRETCFOU, LWRTCFOU are set to .FALSE.
!   when at least LRCOEF=.TRUE. 
IF( (LRCOEF) .AND. (LRETCFOU.OR.LWRTCFOU) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SU0PHY: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,'(1X,A,A,A)') ' SU0PHY: options ',&
   & 'LRCOEF ',&
   & 'require that:'  
  WRITE(KULOUT,'(1X,A)') '         - LRETCFOU=.FALSE. '
  WRITE(KULOUT,'(1X,A)') '         - LWRTCFOU=.FALSE. '
ENDIF

! * Test that option LSRCONT is set to .FALSE. when LSRCON=.FALSE. 
IF( (LSRCONT) .AND. (.NOT.LSRCON) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SU0PHY: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,'(1X,A,A,A)') ' SU0PHY: options ',&
   & 'LSRCONT ',&
   & 'require that:'  
  WRITE(KULOUT,'(1X,A)') '         - LSRCON=.TRUE. '
ENDIF

! * Use NDPSFI=1 when not available or inconsistent with some other options?
!IF(LGWADV .AND. (NDPSFI==1)) THEN
  ! * these two options actually can run together on an informatic
  !   point of view, but the current assumptions done with LGWADV=T
  !   are that a particle which is on the Earth surface remains
  !   on the Earth surface; this is equivalent to assume that "etadot_surf"
  !   is always zero. This condition is satisfied when NDPSFI=0
  !   but not when NDPSFI=1 which gives a non-zero value to "etadot_surf".
  !   So one forbids combination "LGWADV .AND. (NDPSFI==1)" which leads
  !   to inconsistencies in the model at the surface.
!  IERR=IERR+1
!  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SU0PHY: ERROR NR ',IERR,' !!!'
!  WRITE(KULOUT,'(1X,A,A,A)') ' SU0PHY: options ',&
!   & 'NDPSFI=1 ',&
!   & 'require that:'  
!  WRITE(KULOUT,'(1X,A)') '         - LGWADV=.FALSE. '
!ENDIF

! TESTS FOR ALARO-0 (cloud)
IF (LPIL.AND..NOT.LCONDWT) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SU0PHY: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,'(1X,A)') ' SU0PHY: LPIL requires LCONDWT '
ENDIF

! TESTS FOR ALARO-0 (prognostic convection)
IF (LCVPRO) THEN
  LCDDPRO=.TRUE.
  IF (NCONF /= 1)THEN
    WRITE(KULOUT,'(1X, A)') 'WARNING: LCVPRO requires NCONF=1'
    LCVPRO=.FALSE.
    LCDDPRO=.FALSE.
  ENDIF
ENDIF

IF (LPHSPSH .AND. .NOT.(L3MT.OR.LSTRAPRO)) THEN
  ! * For the time being, LPHSPSH requires L3MT.OR.LSTRAPRO (missing code
  !   in CPTEND in the other cases, to be coded later).
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SU0PHY: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,'(1X,A)') ' SU0PHY: LPHSPSH requires L3MT.OR.LSTRAPRO '
ENDIF

! * Use simplified physics with LMPA=T?
IF(LMPA.AND.LSIMPH) THEN
  ! * Simplified physics is not consistent with AROME physics
  !   concerning most of the physical parameterizations.
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SU0PHY: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,'(1X,A)') ' SU0PHY: LMPA and LSIMPH cannot be both to T.'
ENDIF


IF (IERR >= 1) THEN
  CALL FLUSH(KULOUT)
  CALL ABOR1(' SU0PHY: ABOR1 CALLED')
ENDIF

!     ------------------------------------------------------------------

!*       4.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMPHY '')')
WRITE(UNIT=KULOUT,FMT='('' LMPHYS= '',L5,'' LREASUR = '',L5)')LMPHYS,LREASUR
if (1==0) then       !!!!! A REVOIR (MPL)
WRITE(UNIT=KULOUT,FMT='('' LCONDWT='',L5,'' LCVPP = '',L5,&
 & '' LCVPPKF = '',L5,&
 & '' LNEIGE = '',L5,'' LRNUMX = '',L5,'' LCLSATUR = '',L5,'' LHUNEG = '',L5,&
 & '' LVOIGT = '',L5,'' LVFULL = '',L5,&
 & '' L2PHYS = '',L5,&
 & '' LCVDD = '',L5,'' LCAPE= '',L5,'' LSRCON = '',L5,&
 & '' LCVLIS = '',L5,'' LCVCAS= '',L5,'' LSRCONT= '',L5,'' LO3ABC= '',L5,  &
 & '' LAEROSEA = '',L5,'' LAEROLAN = '',L5,'' LAEROSOO = '',L5, &
 & '' LAERODES = '',L5,'' LAEROVOL = '',L5,'' LAEROSUL = '',L5, &
 & '' LRELAXT  = '',L5,'' LRELAXW  = '',L5,&
 & '' LSLC = '',L5,'' LRRGUST= '',L5, &
 & '' LNEWD = '',L5,'' LGLT= '',L5,'' LDIFCONS= '',L5, &
 & '' LECT = '',L5,'' LECTFL = '',L5,'' LPBLE = '',L5,'' LCVPGY= '',L5, &
 & '' LPTKE= '',L5,'' L3MT= '',L5 &
 & )')&
 & LCONDWT,LCVPP,LCVPPKF,LNEIGE,LRNUMX,LCLSATUR,LHUNEG,LVOIGT,LVFULL,L2PHYS,LCVDD,&
 & LCAPE,LSRCON,LCVLIS,LCVCAS,LSRCONT,LO3ABC,LAEROSEA,LAEROLAN,LAEROSOO, &
 & LAERODES,LAEROVOL,LAEROSUL,LRELAXT,LRELAXW,LSLC,LRRGUST,LNEWD,LGLT,&
 & LDIFCONS,LECT,LECTFL,LPBLE,LCVPGY,LPTKE,L3MT 
else
  print*,'OKL'
  print*,'LOGICS',&
 & LCONDWT,LCVPP,LCVPPKF,LNEIGE,LRNUMX,LCLSATUR,LHUNEG,LVOIGT,LVFULL,L2PHYS,LCVDD,&
 & LCAPE,LSRCON,LCVLIS,LCVCAS,LSRCONT,LO3ABC,LAEROSEA,LAEROLAN,LAEROSOO, &
 & LAERODES,LAEROVOL,LAEROSUL,LRELAXT,LRELAXW,LSLC,LRRGUST,LNEWD,LGLT,&
 & LDIFCONS,LECT,LECTFL,LPBLE,LCVPGY,LPTKE,L3MT 
endif
WRITE(UNIT=KULOUT,FMT='('' LCVRA = '',L5 ,'' LFPCOR = '',L5&
 & ,'' LNOIAS = '',L5 &
 & ,'' LGWD = '',L5,'' LHMTO = '',L5,'' LRAY = '',L5 &
 & ,'' LSFHYD = '',L5,'' LSTRA = '',L5,'' LTHERMO = '',L5,/&
 & ,'' LVDIF = '',L5,'' LNEBCO = '',L5,'' LNEBT = '',L5,/&
 & ,'' LRRMES = '',L5,'' LOZONE = '',L5,'' LNEBR = '',L5,/&
 & ,'' LSNV = '',L5,'' LSOLV = '',L5,'' LFGEL = '',L5,/&
 & ,'' LVGSN ='',L5,'' LND2DIFF ='',L5,/&
 & ,'' LNEBN = '',L5,'' LNEBNXR = '',L5,'' LQXRTGH = '',L5,'' LSTRAS = '',L5,/&
 & ,'' LHUCN = '',L5,'' LGWDC     = '',L5,'' LRAYFM = '',L5,'' LRAYFM15= '',L5,/&
 & ,'' LRAYLU  = '',L5,'' LREWS  = '',L5,'' LRPROX = '',L5 /&
 & ,'' LRMIX   = '',L5,'' LRSTAB = '',L5,''LRAUTOEV = '',L5 /&
 & ,'' LRTPP = '',L5,'' LRTDL = '',L5,'' LRAYPL = '',L5 &
 & )')&
 & LCVRA,LFPCOR,LNOIAS,LGWD,LHMTO,LRAY,LSFHYD,LSTRA,LTHERMO,LVDIF,LNEBCO,&
 & LNEBT,LRRMES,LOZONE,LNEBR,LSNV,LSOLV,LFGEL,LVGSN,LND2DIFF,&
 & LNEBN,LNEBNXR,LQXRTGH,LSTRAS,LHUCN,LGWDC,&
 & LRAYFM,LRAYFM15,LRAYLU,LREWS,&
 & LRPROX,LRMIX,LRSTAB,LRAUTOEV,LRTPP,LRTDL,LRAYPL
WRITE(UNIT=KULOUT,FMT='('' NBITER = '',I2,'' NDPSFI = '',I2,'' NPHYREP = '',I2)') &
 & NBITER,NDPSFI,NPHYREP  
WRITE(UNIT=KULOUT,FMT='('' NOIR = '',I2)')NOIR
!WRITE(UNIT=KULOUT,FMT='('' NPHYINT = '',I2)')NPHYINT

WRITE(UNIT=KULOUT,FMT='('' - - - - - - - - - -'')')
WRITE(UNIT=KULOUT,FMT='('' COMMON YOMPHY'')')
WRITE(UNIT=KULOUT,FMT='('' - - - - - - - - - -'')')
WRITE(UNIT=KULOUT,FMT='(&
 & '' LECT    = '',L5,'' LAJUCV  = '',L5,/, &
 & '' LCVRAV3 = '',L5, &
 & '' LNEBGR  = '',L5,'' LNEBGY  = '',L5,'' LBCCOND = '',L5 &
 & )') LECT,    LAJUCV,  &
 & LCVRAV3, LNEBGR, LNEBGY, LBCCOND

WRITE(UNIT=KULOUT,FMT='('' - - - - - - - - -'')')
WRITE(UNIT=KULOUT,FMT='('' MICROPHYSICS KEY '')')
WRITE(UNIT=KULOUT,FMT='('' - - - - - - - - -'')')
WRITE(UNIT=KULOUT,FMT='(  '' LPROCLD = '',L5 &
 &,'' LEVAPP  = '',L5,'' LCOLLEC = '',L5,'' LNSMLIS  = '',L5 &
 &,'' LADJCLD  = '',L5,'' LAUTONEB  = '',L5,'' LSSD = '',L5)')&
 &LPROCLD,LEVAPP,LCOLLEC,LNSMLIS,LADJCLD,LAUTONEB,LSSD

WRITE(UNIT=KULOUT,FMT='('' ALARO-0 cloud '')')
WRITE(UNIT=KULOUT,FMT='('' LPIL = '',L5,'' LSTRAPRO = '',L5 &
  &,'' LNEWSTAT = '',L5 &
  &,'' LPHSPSH = '',L5 ,'' LSMROT = '',L5 &
  &,'' LSMTPS = '',L5 ,'' LSMNIMBT = '',L5 &
  &,'' L1DRHCRI = '',L5 ,'' LGWRHCRI = '',L5 &
  &,'' NSMTBOT = '',I2,'' NSMDNEB = '',I2 &
  &,'' NPRAG = '',I2,'' NPRAC = '',I2 &
  &,'' NPRRI = '',I2 &
  &) ')LPIL,LSTRAPRO,LNEWSTAT,LPHSPSH,LSMROT,LSMTPS, LSMNIMBT, &
  &    L1DRHCRI, LGWRHCRI,&
  &    NSMTBOT,NSMDNEB, NPRAG,NPRAC,NPRRI

WRITE(UNIT=KULOUT,FMT='('' ALARO-0 prognostic convection  '')')
WRITE(UNIT=KULOUT,&
  &FMT='('' LCVPRO ='',L1,'' LCDDPRO ='',L1,'' LSCMF ='',L1)') &
  &    LCVPRO,LCDDPRO,LSCMF

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMSIMPHL '')')
WRITE(UNIT=KULOUT,FMT='('' LSIMPH= '',L5,'' LTRAJPS = '',L5 &
 & ,'' LTRAJPST = '',L5 &
 & ,'' LSMOOTHD = '',L5,'' LSMOOTHA = '',L5,'' LSMOOTHB = '',L5 &
 & ,'' LCLOUDS = '',L5 )')&
 & LSIMPH,LTRAJPS,LTRAJPST,LSMOOTHD,LSMOOTHA, &
 & LSMOOTHB,LCLOUDS  
WRITE(UNIT=KULOUT,FMT='('' LCVRASP = '',L5,'' LGWDSP = '',L5 &
 & ,'' LRAYSP = '',L5,'' LSTRASP = '',L5,'' LVDIFSP = '',L5 &
 & ,'' LRRMESSP = '',L5)')&
 & LCVRASP,LGWDSP,LRAYSP,LSTRASP,LVDIFSP,LRRMESSP  
WRITE(UNIT=KULOUT,FMT='('' COMMON YOMRCOEF '')')
WRITE(UNIT=KULOUT,FMT='('' LRCOEF= '',L5 &
 & ,'' LTLADDIA = '',L5,'' LGLOBRAD = '',L5)')&
 & LRCOEF,LTLADDIA,LGLOBRAD  

! TEST OF CONSISTENCY FOR RADIATION SCHEMES

IF(LMPHYS) THEN
  IF((LRAY.AND.LRAYFM).OR.(LRAY.AND.LRAYFM15).OR.(LRAYFM15.AND.LRAYFM)) THEN
    WRITE(NULOUT,'('' WARNING - 2 RADIATION SCHEMES... '')')
    CALL ABOR1('SU0PHY: ABOR1 CALLED')
  ENDIF
ENDIF

!     TESTS OF CONSISTENCY INSIDE YOMPHY

IF(LSTRAS.AND.LRNUMX)THEN
  WRITE(NULOUT,FMT='('' ACPLUIS AND LRNUMX ARE NOT COMPATIBLE''&
   & ,''      FOR THE MOMENT'')')  
  CALL ABOR1('SU0PHY: ABOR1 CALLED')
ENDIF

!     TESTS OF CONSISTENCY BETWEEN YOMPHY AND INITIALIZATION

!IF (.NOT.LMPHYS.AND..NOT.LEPHYS) THEN
! Digital Filter Initialisation
!  IF ((NEINI == 2.OR.NEINI == 4).AND.NEDFI >= 2) THEN
!    IF (NEDFI <= 5) THEN
!      WRITE(NULOUT,FMT='('' YOMPHY AND YEMDFI ARE NOT COMPATIBLE'')')
!      CALL ABOR1(' DIABATIC DFI WITHOUT PHYSICS !')
!    ELSE
!      WRITE(NULOUT,FMT='('' CAUTION : FORWARD DFI WITHOUT PHYSICS'')')
!    ENDIF
!  ENDIF
!ENDIF

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMARPHY '')')
WRITE(UNIT=KULOUT,FMT='('' LMPA = '',L5, &
 & '' LMICRO = '',L5,'' LTURB = '',L5, &
 & '' LMSE = '',L5,'' LKFBCONV = '',L5,&
 & '' LKFBD = '',L5,'' LKFBS = '',L5,&
 & '' LUSECHEM = '',L5,'' LORILAM = '',L5,'' LRDUST = '',L5)')&
 & LMPA,LMICRO,LTURB,LMSE,LKFBCONV,LKFBD,LKFBS,LUSECHEM,&
 & LORILAM,LRDUST

WRITE(UNIT=KULOUT,FMT='('' COMMON YOEPHY '')')
WRITE(UNIT=KULOUT,FMT='('' LEPHYS = '',L5, &
 & '' LECOND = '',L5,'' LECUMF = '',L5 &
 & )')&
 & LEPHYS,LECOND,LECUMF  
!WRITE(NULOUT,'("RCLDTOPP=",F10.2)') RCLDTOPP
WRITE(UNIT=KULOUT,FMT='('' LEDCLD = '',L5,'' LEGWDG = '',L5 &
 & ,'' LEOZOC = '',L5,'' LEQNGT = '',L5 &
 & ,'' LEO3CH = '',L5)')&
 & LEDCLD,LEGWDG,LEOZOC,LEQNGT,LEO3CH  
WRITE(UNIT=KULOUT,FMT='('' LERADI = '',L5 &
 & ,'' LESHCV = '',L5,'' LESURF = '',L5,'' LEVDIF = '',L5 &
 & ,'' LEOCWA = '',L5,'' LEOCCO = '',L5,'' LEOCSA = '',L5 &
 & ,'' LECURR = '',L5,'' RTHRFRTI = '',F7.2 &
 & )')&
 & LERADI,LESHCV,LESURF,LEVDIF,LEOCWA,LEOCCO,LEOCSA,LECURR,RTHRFRTI 
WRITE(UNIT=KULOUT,FMT='('' LEPCLD = '',L5,'' LEMETHOX= '',L5 &
 & ,'' LE4ALB = '',L5 &
 & )')&
 & LEPCLD,LEMETHOX,LE4ALB  

! Do not allow storage of trajectory in TL physics if tracers 
! are to be included in linearized physics (temporary) 
!IF (LETRAJP .AND. LETRAJPT .AND. LTRACLNPH) THEN
!  WRITE(UNIT=KULOUT,FMT='('' CAUTION : STORAGE OF TRAJECTORY IN LINEAR. PHYSICS NOT YET '')')
!  WRITE(UNIT=KULOUT,FMT='('' IMPLEMENTED FOR TRACERS: LETRAJPT WILL BE RESET TO .FALSE. '')')
!ENDIF

!IF (LETRAJP) THEN
!  WRITE(UNIT=KULOUT,FMT='('' ======================== '')')
!  WRITE(UNIT=KULOUT,FMT='('' LINEAR PHYSICS ACTIVATED '')')
!  WRITE(UNIT=KULOUT,FMT='('' ======================== '')')
!  WRITE(UNIT=KULOUT,FMT='('' LERADI2 = '',L4 &
!   & ,'' LERADS2 = '',L4 &
!   & ,'' LERADSW2= '',L4 &
!   & ,'' LERADN2 = '',L4 &
!   & ,'' LERADFL2= '',L4 &
!   & )')&
!   & LERADI2,LERADS2,LERADSW2,LERADN2,LERADFL2  
!  WRITE(UNIT=KULOUT,FMT='('' LEDCLD2 = '',L4 &
!   & ,'' LENCLD2 = '',L4 &
!   & ,'' LEVAPLS2= '',L4 &
!   & ,'' LREGCL  = '',L4 &
!   & ,'' LNCLIN  = '',L4 &
!   & )')&
!   & LEDCLD2,LENCLD2,LEVAPLS2,LREGCL,LNCLIN
!  WRITE(UNIT=KULOUT,FMT='('' LECUMF2 = '',L4 &
!   & ,'' LECUBM2 = '',L4 &
!   & ,'' LECOND2 = '',L4    &
!   & )')&
!   & LECUMF2,LECUBM2,LECOND2  
!  WRITE(UNIT=KULOUT,FMT='('' LEVDIF2 = '',L4 &
!   & ,'' LEGWDG2 = '',L4 &
!   & ,'' LEQNGT2 = '',L4 &
!   & ,'' LESURF2 = '',L4 &
!   & )')&
!   & LEVDIF2,LEGWDG2,LEQNGT2,LESURF2  
!  WRITE(UNIT=KULOUT,FMT='('' LEKPERT = '',L4,'' LENOPERT = '',L4)')&
!   & LEKPERT,LENOPERT
!  WRITE(UNIT=KULOUT,FMT='('' ======================== '')')
!  WRITE(UNIT=KULOUT,FMT='('' LETRAJPT = '',L4)')LETRAJPT
!  WRITE(UNIT=KULOUT,FMT='('' ======================== '')')
!  WRITE(UNIT=KULOUT,FMT='('' LTRACLNPH = '',L4)')LTRACLNPH
!  WRITE(UNIT=KULOUT,FMT='('' ======================== '')')
!ELSE
!  WRITE(UNIT=KULOUT,FMT='('' ============================ '')')
!  WRITE(UNIT=KULOUT,FMT='('' LINEAR PHYSICS NOT ACTIVATED '')')
!  WRITE(UNIT=KULOUT,FMT='('' ============================ '')')
!  WRITE(UNIT=KULOUT,FMT='('' LENCLD2 = '',L4)')LENCLD2
!  WRITE(UNIT=KULOUT,FMT='('' LEVAPLS2= '',L4)')LEVAPLS2
!ENDIF

WRITE(UNIT=KULOUT,FMT='('' LSPHLC = '',L5,'' LVDFLC = '',L5 &
 & ,'' LVDFDS = '',L5,'' LSDRLC = '',L5,'' LSDRDS = '',L5)')&
 & LSPHLC,LVDFLC,LVDFDS,LSDRLC,LSDRDS  
WRITE(UNIT=KULOUT,FMT='('' LCZDEB = '',L5,'' LZMCON = '',L5 &
 & ,'' LKEXP  = '',L5)') LCZDEB,LZMCON,LKEXP  
WRITE(UNIT=KULOUT,FMT='('' ALPHA  = '',E16.6)') ALPHA
WRITE(UNIT=KULOUT,FMT='('' AH0    = '',E16.6,''m'')') AH0
WRITE(UNIT=KULOUT,FMT='('' USTARL = '',E16.6,''m s-1'')') USTARL
WRITE(UNIT=KULOUT,FMT='('' USTARS = '',E16.6,''m s-1'')') USTARS
WRITE(UNIT=KULOUT,FMT='('' ALANDZ0= '',E16.6,''m'')') ALANDZ0
WRITE(UNIT=KULOUT,FMT='('' ASEAZ0 = '',E16.6,''m'')') ASEAZ0

IF (LAGPHY) THEN
  WRITE(UNIT=KULOUT,FMT='('' LAGPHY = '',L5,'' THE PHYSICS PAC&
   & KAGE IS CALLED AFTER THE DYNAMICS IN LAGGED MODE'')')LAGPHY  
ELSE
  WRITE(UNIT=KULOUT,FMT='('' LAGPHY = '',L5,'' THE PHYSICS PAC&
   & KAGE IS CALLED BEFORE THE DYNAMICS '')')LAGPHY  
ENDIF
WRITE(UNIT=KULOUT,FMT=*) 'CGMIXLEN=',TRIM(CGMIXLEN)
WRITE(UNIT=KULOUT,FMT='('' LPRGML = '',L4)') LPRGML

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU0PHY',1,ZHOOK_HANDLE)
END SUBROUTINE SU0PHY
