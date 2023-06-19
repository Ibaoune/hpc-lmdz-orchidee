SUBROUTINE SUPHEC(KULOUT)

!**** *SUPHEC - INITIALISES PHYSICAL CONSTANTS OF UNCERTAIN VALUE.
!               WITHIN THE E.C.M.W.F. PHYSICS PACKAGE

!     PURPOSE.
!     --------

!          THIS ROUTINE SETS THE VALUES FOR THE PHYSICAL CONSTANTS USED
!     IN THE PARAMETERIZATION ROUTINES WHENEVER THESE VALUES ARE NOT
!     KNOWN WELL ENOUGH TO FORBID ANY TUNING OR WHENEVER THEY ARE
!     SUBJECT TO AN ARBITRARY CHOICE OF THE MODELLER. THESE CONSTANTS
!     ARE DISTRIBUTED IN COMMON DECKS *YOEXXXX* WHERE XXXX CORRESPONDS
!     TO THE INDIVIDUAL PHYSICAL PARAMETRIZATION

!**   INTERFACE.
!     ----------

!          *SUPHEC* IS CALLED FROM *SUPHY*

!     METHOD.
!     -------

!          NONE.

!     EXTERNALS.
!     ----------

!          *SUECRAD*, *SUCUMF*, *SUCUMF2*,*SUVDFS*, *SUSURF*
!          *SUECRAD15*, *SUCLOP15*
!          *SUGWD*, *SUCLD*, *SUCOND*, *SUPHLI*, *SUMETHOX*

!     REFERENCE.
!     ----------

!          SEE PHYSICAL ROUTINES FOR AN EXACT DEFINITION OF THE
!     CONSTANTS.

!     AUTHOR.
!     -------
!          J.-J. MORCRETTE  E.C.M.W.F.    91/06/15  ADAPTATION TO I.F.S.

!     MODIFICATIONS
!     -------------
!          MAY 1997 : M. Deque  - Frozen FMR
!          APRIL 1998: C. JAKOB - ADD METHANE OXIDATION
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P.Viterbo     24-May-2004 surf library
!        P.Viterbo     03-Dec-2004 Include user-defined RTHRFRTI
!        M.Ko"hler     03-Dec-2004 cp,moist=cp,dry
!        P.Viterbo     10-Jun-2005 Externalise surf
!        R. El Khatib & J-F Estrade  20-Jan-2005 Default PRSUN for FMR15
!        D.Salmond     22-Nov-2005 Mods for coarser/finer physics
!        P. Lopez      21-Aug-2006 Added call to SUCUMF2 
!                                 (new linearized convec)
!        JJMorcrette   20060525    MODIS albedo
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMDPHY  , ONLY : NTILES 
USE SURFACE_FIELDS, ONLY : YSP_SBD
USE YOELW    , ONLY : NSIL     ,TSTAND   ,XP
USE YOESW    , ONLY : RSUN
USE YOMSW15  , ONLY : RSUN15
USE YOMDIM   , ONLY : NFLEVG   ,NSMAX, NGPBLKS, NPROMA
USE YOMGEM   , ONLY : VBH      ,VAH      ,VP00, VAF   , VBF
USE YOMCST   , ONLY : RD       ,RV       ,RCPD     ,&
 & RLVTT    ,RLSTT    ,RLMLT    ,RTT      ,RATM 
!USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
! & R4IES    ,R5LES    ,R5IES    ,RVTMP2   ,RHOH2O   ,&
! & R5ALVCP  ,R5ALSCP  ,RALVDCP  ,RALSDCP  ,RALFDCP  ,&
! & RTWAT    ,RTBER    ,RTBERCU  ,RTICE    ,RTICECU  ,&
! & RTWAT_RTICE_R      ,RTWAT_RTICECU_R    ,&
! & RKOOP1   ,RKOOP2 
USE YOMPHY   , ONLY : LRAYFM15
!USE YOERAD   , ONLY : NSW      ,NTSW     ,&
! NSW mis dans .def MPL 20140211
USE YOERAD   , ONLY : NTSW     ,&
 & LCCNL    ,LCCNO    ,&
 & RCCNSEA  ,RCCNLND
USE YOE_TILE_PROP, ONLY : RUSTRTI, RVSTRTI, RAHFSTI, REVAPTI, RTSKTI
USE YOEPHY   , ONLY : RTHRFRTI ,LEOCWA   ,LEOCCO   ,LEOCSA, LE4ALB
USE YOEVDF   , ONLY : NVTYPES
USE YOMCOAPHY   , ONLY : NPHYINT
USE YOM_PHYS_GRID ,ONLY : PHYS_GRID
USE YOMCT0  , ONLY  : LSCMEC   ,LROUGH   ,REXTZ0M  ,REXTZ0H
USE vertical_layers_mod, ONLY: ap,bp

IMPLICIT NONE
include "YOETHF.h"
include "clesphys.h"

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTERFACE
#include "susurf.h"
#include "surf_inq.h"
END INTERFACE

#include "gppre.intfb.h"
#include "sucld.intfb.h"
#include "sucldp.intfb.h"
#include "suclop.intfb.h"
#include "suclop15.intfb.h"
#include "sucond.intfb.h"
#include "sucumf.intfb.h"
#include "sucumf2.intfb.h"
#include "suecrad.intfb.h"
#include "suecrad15.intfb.h"
#include "sugwd.intfb.h"
#include "sumethox.intfb.h"
#include "suphli.intfb.h"
#include "suvdf.intfb.h"
#include "suvdfs.intfb.h"
#include "suwcou.intfb.h"

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZPRES(0:NFLEVG),ZPRESF(NFLEVG), ZETA(NFLEVG),ZETAH(0:NFLEVG)

INTEGER(KIND=JPIM) :: JK,ISMAX,JLEV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*         0.2    DEFINING DERIVED CONSTANTS FROM UNIVERSAL CONSTANTS
!                 ---------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPHEC',0,ZHOOK_HANDLE)
!
  IF (OK_BAD_ECMWF_THERMO) THEN
!
     ! Modify constants defined in suphel.F90 and set RVTMP2 to 0.
     ! CALL GSTATS(1811,0) ! MPL 28.11.08
     ! RVTMP2=RCPV/RCPD-1.0_JPRB   !use cp,moist
     RVTMP2=0.0_JPRB              !neglect cp,moist
     RHOH2O=RATM/100._JPRB
     R2ES=611.21_JPRB*RD/RV
     R3LES=17.502_JPRB
     R3IES=22.587_JPRB
     R4LES=32.19_JPRB
     R4IES=-0.7_JPRB
     R5LES=R3LES*(RTT-R4LES)
     R5IES=R3IES*(RTT-R4IES)
     R5ALVCP=R5LES*RLVTT/RCPD
     R5ALSCP=R5IES*RLSTT/RCPD
     RALVDCP=RLVTT/RCPD
     RALSDCP=RLSTT/RCPD
     RALFDCP=RLMLT/RCPD
     RTWAT=RTT
     RTBER=RTT-5._JPRB
     RTBERCU=RTT-5.0_JPRB
     RTICE=RTT-23._JPRB
     RTICECU=RTT-23._JPRB
     
     RTWAT_RTICE_R=1.0_JPRB/(RTWAT-RTICE)
     RTWAT_RTICECU_R=1.0_JPRB/(RTWAT-RTICECU)
     IF(NPHYINT == 0) THEN
       ISMAX=NSMAX
     ELSE
       ISMAX=PHYS_GRID%NSMAX
     ENDIF
     
     RKOOP1=2.583_JPRB
     RKOOP2=0.48116E-2_JPRB
     
  ELSE 
     ! Keep constants defined in suphel.F90
     RTICE=RTT-23._JPRB
!
  ENDIF  ! (OK_BAD_ECMWF_THERMO)

!     ------------------------------------------------------------------
!*         0.5    DEFINE STANDARD ATMOSPHERE VERTICAL CONFIGURATION
!                 -------------------------------------------------
!ALLOCATE(VBH    (0:MAX(JPMXLE,NFLEVG)))  from suallo.F90
!! 
ALLOCATE(VAH    (0:NFLEVG))  ! Ajout ALLOCATE MPL 200509
ALLOCATE(VBH    (0:NFLEVG))
ALLOCATE(VAF    (NFLEVG))
ALLOCATE(VBF    (NFLEVG))
! Commente par MPL 28.11.08, puis decommente le 19.05.09
VP00=101325.     !!!!! A REVOIR (MPL)
ZPRES(NFLEVG)=VP00
! on recupere ap et bp de dyn3d (vertical_layers_mod) MPL 19.05.09
! Attention, VAH et VBH sont inverses, comme les niveaux
! plev(l)=PAPRS(klon,nlayer+1-l) de 1 a nlayer (apllmd.F)
DO JLEV = 0, NFLEVG  
!  VAH(JLEV)=ap(JLEV+1)ap(JLEV+1)
!  VBH(JLEV)=bp(JLEV+1)
!  print *,'SUPHEC: jlev ap bp',JLEV,ap(JLEV+1),bp(JLEV+1)
   VAH(JLEV)=ap(NFLEVG+1-JLEV)
   VBH(JLEV)=bp(NFLEVG+1-JLEV)
ENDDO
! Calcul de VAF et VBF, analogues de VAH et VBH mais aux niveaux pleins
DO JLEV = 1, NFLEVG   
   VAF(JLEV)=(VAH(JLEV)+VAH(JLEV-1))/2.
   VBF(JLEV)=(VBH(JLEV)+VBH(JLEV-1))/2.
ENDDO

! Appel a GPPRE commente par MPL 28.11.08, puis decommente le 19.05.09
CALL GPPRE ( 1 ,1, 1, NFLEVG, VAH, VBH, ZPRES, ZPRESF )

DO JK=0,NFLEVG
  ZETAH(JK)= ZPRES(JK)/ZPRES(NFLEVG)
ENDDO
DO JK=1,NFLEVG
  ZETA(JK)= ZPRESF(JK)/ZPRES(NFLEVG)
ENDDO

!     ------------------------------------------------------------------
!*         1.     SETTING CONSTANTS FOR DIAGNOSTIC CLOUD SCHEME
!                 ---------------------------------------------

!CALL SUCLD ( NFLEVG , ZETA ) ! MPL 28.11.08

!     ------------------------------------------------------------------

!*         2.     SETTING CONSTANTS FOR LARGE-SCALE CONDENSATION SCHEME
!                 -----------------------------------------------------

!CALL SUCOND ( KULOUT , NFLEVG , ZETA ) ! MPL 28.11.08

!     ------------------------------------------------------------------

!*         3.     SETTING CONSTANTS FOR CONVECTION SCHEME
!                 ---------------------------------------

!CALL SUCUMF(ISMAX)     ! MPL 28.11.08

!     ------------------------------------------------------------------

!*         3.     SETTING CONSTANTS FOR NEW LINEARIZED CONVECTION SCHEME
!                 ------------------------------------------------------

!CALL SUCUMF2(ISMAX)  ! MPL 28.11.08

!     ------------------------------------------------------------------
!*         4.     SETTING CONSTANTS FOR GRAVITY WAVE DRAG SCHEME
!                 ----------------------------------------------

!CALL SUGWD (KULOUT, NFLEVG, VAH, VBH )   ! MPL 28.11.08

!     ------------------------------------------------------------------

!*         5.     SETTING CONSTANTS FOR VERTICAL DIFFUSION
!                 ----------------------------------------

!CALL SUVDFS     ! MPL 28.11.08

!CALL SUVDF      ! MPL 28.11.08

!cccc CALL SUVDFD ( NABLPFR, ABLPLL ) cccccccccccccccccccccccccccccccccc

!     ------------------------------------------------------------------

!*         6.     SETTING CONSTANTS FOR RADIATION SCHEME
!                 --------------------------------------

IF (LRAYFM15) THEN
  CALL SUECRAD15 (KULOUT, NFLEVG, ZETAH )
ELSE
  CALL SUECRAD (KULOUT, NFLEVG, ZETAH )
ENDIF

!     ------------------------------------------------------------------
!*         7.     SETTING CONSTANTS FOR SURFACE SCHEME
!                 ------------------------------------

!IF (LRAYFM15) THEN
!   CALL SUSURF(KSW=NSW,KCSS=YSP_SBD%NLEVS,KSIL=NSIL,KTILES=NTILES,KTSW=NTSW,&
!    & LD_LLCCNL=LCCNL,LD_LLCCNO=LCCNO,&
!    & LD_LEOCWA=LEOCWA,LD_LEOCCO=LEOCCO,LD_LEOCSA=LEOCSA,LD_LLE4ALB=LE4ALB,&
!    & LD_LSCMEC=LSCMEC,LD_LROUGH=LROUGH,PEXTZ0M=REXTZ0M,PEXTZ0H=REXTZ0H,&
!    & PTHRFRTI=RTHRFRTI,PTSTAND=TSTAND,PXP=XP,PRCCNSEA=RCCNSEA,PRCCNLND=RCCNLND,&
!    & PRSUN=RSUN15)
!ELSE
!   CALL SUSURF(KSW=NSW,KCSS=YSP_SBD%NLEVS,KSIL=NSIL,KTILES=NTILES,KTSW=NTSW,&
!    & LD_LLCCNL=LCCNL,LD_LLCCNO=LCCNO,&
!    & LD_LEOCWA=LEOCWA,LD_LEOCCO=LEOCCO,LD_LEOCSA=LEOCSA,LD_LLE4ALB=LE4ALB,&
!    & LD_LSCMEC=LSCMEC,LD_LROUGH=LROUGH,PEXTZ0M=REXTZ0M,PEXTZ0H=REXTZ0H,&
!    & PTHRFRTI=RTHRFRTI,PTSTAND=TSTAND,PXP=XP,PRCCNSEA=RCCNSEA,PRCCNLND=RCCNLND,&
!    & PRSUN=RSUN)
!ENDIF


!CALL SURF_INQ(KNVTYPES=NVTYPES)


!          7.1    Allocate working arrays
!ALLOCATE(RUSTRTI(NPROMA,NTILES,NGPBLKS))
!ALLOCATE(RVSTRTI(NPROMA,NTILES,NGPBLKS))
!ALLOCATE(RAHFSTI(NPROMA,NTILES,NGPBLKS))
!ALLOCATE(REVAPTI(NPROMA,NTILES,NGPBLKS))
!ALLOCATE(RTSKTI (NPROMA,NTILES,NGPBLKS))
!RUSTRTI(:,:,:) = 0.0_JPRB
!RVSTRTI(:,:,:) = 0.0_JPRB
!RAHFSTI(:,:,:) = 0.0_JPRB
!REVAPTI(:,:,:) = 0.0_JPRB
!RTSKTI (:,:,:) = 0.0_JPRB
!CALL GSTATS(1811,1)

!     ------------------------------------------------------------------

!*         8.     SETTING CONSTANTS FOR CLOUD OPTICAL PROPERTIES
!                 ----------------------------------------------

IF (LRAYFM15) THEN
  CALL SUCLOP15
ELSE
  CALL SUCLOP
ENDIF

!     ------------------------------------------------------------------

!*         9.     SETTING CONSTANTS FOR PROGNOSTIC CLOUD SCHEME
!                 ----------------------------------------------

!CALL SUCLDP

!     ------------------------------------------------------------------

!*        10.     SETTING CONSTANTS FOR WAVE COUPLING
!                 -----------------------------------

!CALL SUWCOU

!     ------------------------------------------------------------------
!*         11.   SETTING CONSTANTS FOR LINEARIZED PHYSICS
!                ----------------------------------------

!CALL SUPHLI

!     ------------------------------------------------------------------
!*         12.   SETTING CONSTANTS FOR METHANE OXIDATION
!                ---------------------------------------

!CALL SUMETHOX

!     ------------------------------------------------------------------

WRITE(UNIT=KULOUT,FMT='('' SUPHEC IS OVER '')')

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPHEC',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHEC
