SUBROUTINE GPXYB(KPROMA,KSTART,KPROF,KFLEV,PVDELB,PVC,&
 & PRES,PDELP,PRDELP,PLNPR,PALPH,PRTGR,&
 & PRPRES,PRPP)  

!**** *GPXYB* - Computes auxillary arrays

!     Purpose.
!     --------
!           COMPUTES AUXILLARY ARRAYS RELATED TO THE HYBRID COORDINATE

!**   Interface.
!     ----------
!        *CALL* *GPXYB(..)

!        Explicit arguments :
!        --------------------
!     KPROMA : dimensioning
!     KSTART : start of work
!     KPROF  : depth of work
!     KFLEV     : vert. dimensioning

!     PVDELB(KPROMA,0:KFLEV) : related to vert. coordinate        (input)
!     PVC   (KPROMA,0:KFLEV) :  "     "      "     "    "         (input)
!     PRES (KPROMA,0:KFLEV)  : HALF LEVEL PRESSURE                (input)
!     PDELP (KPROMA,KFLEV)   : PRESSURE DIFFERENCE ACROSS LAYERS  (output)
!     PRDELP(KPROMA,KFLEV)   : THEIR INVERSE                      (output)
!     PLNPR (KPROMA,KFLEV)   : LOGARITHM OF RATIO OF PRESSURE     (output)
!     PALPH (KPROMA,KFLEV)   : COEFFICIENTS OF THE HYDROSTATICS   (output)
!     PRTGR (KPROMA,KFLEV)   : FOR PRES. GRAD. TERM AND ENE. CONV.(output)
!                              ((rssavnabla prehyd/prehyd)_[layer]
!                              = prtgr_[layer] * (rssavnabla prehyds))
!     PRPRES(KPROMA,KFLEV)   : inverse of HALF LEVEL PRESSURE     (output)
!     PRPP  (KPROMA,KFLEV)   : inverse of PRES(J)*PRES(J-1)       (output)

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.      None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-02-04
!        Modified : 94-10-11 by Radmila Bubnova: correction in the case
!                            of the other approximation of d (ln p).
!        Modified : 99-06-04 Optimisation   D.SALMOND
!        Modified : 02-03-08 K. YESSAD: consistent discretisations of
!                    "alpha" (PALPH) and "prtgr" (PRTGR)
!                    for finite element vertical discretisation
!                    to allow model to run with MF-physics and DDH.
!        Modified : 03-07-07 J. Hague:  Replace divides with reciprocal
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Modified : 15-Feb-2005 by K. YESSAD: ZTOPPRES becomes TOPPRES
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMDYN   , ONLY : NDLNPR   ,RHYDR0
USE YOMCST   , ONLY : RD       ,RCVD
USE YOMGEM   , ONLY : VDELA    ,VAF      ,VBF      ,TOPPRES
USE YOMCVER  , ONLY : LVERTFE

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDELB(KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVC(KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDELP(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRDELP(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLNPR(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALPH(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRTGR(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRPRES(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRPP(KPROMA,KFLEV) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IFIRST, JLEV, JLON, JJ, JTEMP, JM

REAL(KIND=JPRB) :: ZPRESF
REAL(KIND=JPRB) :: ZRPRES(KPROMA,2)
REAL(KIND=JPRB) :: ZPRESFD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPXYB',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       0.    Level to begin normal computations
!              ----------------------------------

! This is introduced to allow the use of GPXYB without the implicit
! assumption that the top level input for pressure is 0 hPa. This
! is used in the surface observation operators where you do not want
! to compute geopotential at all model levels.
! The first block if is for economy (no do loop start up) and the second
! for safety.
!print *,'GPXYB: NDLNPR RHYDR0=',NDLNPR,RHYDR0
TOPPRES=0.1  !!!!! A REVOIR (MPL) 29042010 passe de 0 a 0.1 comme ARPEGE
IF(PRES(KSTART,0) <= TOPPRES)THEN
  IFIRST=2
ELSE
  IFIRST=1
  DO JLON=KSTART,KPROF
    IF(PRES(JLON,0) <= TOPPRES)then
      IFIRST=2
      EXIT
    ENDIF
  ENDDO
ENDIF
!     ------------------------------------------------------------------

!*       1.    COMPUTES EVERYTHING.
!              --------------------

!print *,'NDLNPR LVERTFE',NDLNPR,LVERTFE
IF(NDLNPR == 0) THEN

  IF(LVERTFE) THEN
    DO JLEV=1,KFLEV
      DO JLON=KSTART,KPROF
        PDELP(JLON,JLEV)=VDELA(JLEV) + PVDELB(JLEV)*PRES(JLON,KFLEV)
        PRDELP(JLON,JLEV)=1.0_JPRB/PDELP(JLON,JLEV)
        ZPRESF =VAF(JLEV) + VBF(JLEV)*PRES(JLON,KFLEV)
        ZPRESFD=1.0_JPRB/ZPRESF
        PLNPR(JLON,JLEV)=PDELP(JLON,JLEV)*ZPRESFD
        ! * PRTGR needed for DDH and option LVERCOR=T.
        !   for finite element vertical discretisation,
        !   "prtgr_[layer]" is simply B_[layer]/prehyd_[layer]
        PRTGR (JLON,JLEV)=VBF(JLEV)*ZPRESFD
        ! * PALPH needed for MF physics:
        PALPH(JLON,JLEV)=(PRES(JLON,JLEV)-ZPRESF)*ZPRESFD
      ENDDO
    ENDDO
  ELSE
    JJ=1
    JM=2
    DO JLON=KSTART,KPROF
      ZRPRES(JLON,JM)=1.0_JPRB/PRES(JLON,IFIRST-1)
    ENDDO
    DO JLEV=IFIRST,KFLEV
      DO JLON=KSTART,KPROF
        ZRPRES(JLON,JJ)=1.0_JPRB/PRES(JLON,JLEV)
        PDELP (JLON,JLEV)=PRES(JLON,JLEV)-PRES(JLON,JLEV-1)
        PRDELP(JLON,JLEV)=1.0_JPRB/PDELP(JLON,JLEV)
        PLNPR (JLON,JLEV)=LOG(PRES(JLON,JLEV)*ZRPRES(JLON,JM))
        PRPRES(JLON,JLEV)=ZRPRES(JLON,JJ)
        PALPH (JLON,JLEV)=1.0_JPRB-PRES(JLON,JLEV-1)*PRDELP(JLON,JLEV)&
         & *PLNPR(JLON,JLEV)  
        PRPP  (JLON,JLEV)=ZRPRES(JLON,JJ)*ZRPRES(JLON,JM)
        PRTGR (JLON,JLEV)=PRDELP(JLON,JLEV)&
         & *(PVDELB(JLEV)+PVC(JLEV)*PLNPR(JLON,JLEV)*PRDELP(JLON,&
         & JLEV))  
!       print *,'GPXYB JLEV JLON JJ PRES ZPRES PDELP ', JLEV,JLON,JJ,PRES(JLON,JLEV),ZRPRES(JLON,JJ),PDELP(JLON,JLEV)
!       print *,'GPXYB JLEV JLON JM PRDELP PLNPR ', JLEV,JLON,JM,PRDELP(JLON,JLEV),PLNPR (JLON,JLEV)
!       print *,'GPXYB JLEV JLON JJ PRPRES PALPH ', JLEV,JLON,JJ,PRPRES(JLON,JLEV),PALPH (JLON,JLEV)
!       print *,'GPXYB JLEV JLON PRPP PRTGR PVDELB PVC ', JLEV,JLON,PRPP  (JLON,JLEV),PRTGR (JLON,JLEV),PVDELB(JLEV),PVC(JLEV)
      ENDDO
      JTEMP=JM
      JM=JJ
      JJ=JTEMP
    ENDDO
    DO JLEV=1,IFIRST-1
      DO JLON=KSTART,KPROF
        PDELP (JLON,JLEV)=PRES(JLON,JLEV)-PRES(JLON,JLEV-1)
        PRDELP(JLON,JLEV)=1.0_JPRB/PDELP(JLON,JLEV)
        PLNPR (JLON,JLEV)=LOG(PRES(JLON,1)/TOPPRES)
        PRPRES(JLON,JLEV)=1.0_JPRB/PRES(JLON,1)
        PALPH (JLON,JLEV)=RHYDR0
        PRPP  (JLON,JLEV)=1.0_JPRB/(PRES(JLON,1)*TOPPRES)
        PRTGR (JLON,JLEV)=PRDELP(JLON,JLEV)*PVDELB(JLEV)
      ENDDO
    ENDDO
  ENDIF

ELSEIF(NDLNPR == 1) THEN
  IF(LVERTFE) THEN
    CALL ABOR1(' LVERTFE=.T. NOT COMPATIBLE WITH NDLNPR == 1')
  ENDIF

  DO JLEV=IFIRST,KFLEV
    DO JLON=KSTART,KPROF
      PDELP (JLON,JLEV)=PRES(JLON,JLEV)-PRES(JLON,JLEV-1)
      PRDELP(JLON,JLEV)=1.0_JPRB/PDELP(JLON,JLEV)
      PRPP  (JLON,JLEV)=1.0_JPRB/(PRES(JLON,JLEV)*PRES(JLON,JLEV-1))
      PLNPR (JLON,JLEV)=PDELP(JLON,JLEV)*SQRT(PRPP(JLON,JLEV))
      PALPH (JLON,JLEV)=1.0_JPRB-PRES(JLON,JLEV-1)*PRDELP(JLON,JLEV)&
       & *PLNPR(JLON,JLEV)  
      PRTGR (JLON,JLEV)=PRDELP(JLON,JLEV)&
       & *(PVDELB(JLEV)+PVC(JLEV)*PLNPR(JLON,JLEV)*PRDELP(JLON,&
       & JLEV))  
      PRPRES(JLON,JLEV)=1.0_JPRB/PRES(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLEV=1,IFIRST-1
    DO JLON=KSTART,KPROF
      PDELP (JLON,JLEV)=PRES(JLON,JLEV)
      PRDELP(JLON,JLEV)=1.0_JPRB/PDELP(JLON,JLEV)
      PLNPR (JLON,JLEV)=2.0_JPRB+RCVD/RD
      PALPH (JLON,JLEV)=1.0_JPRB
      PRTGR (JLON,JLEV)=PRDELP(JLON,JLEV)*PVDELB(JLEV)
      PRPRES(JLON,JLEV)=1.0_JPRB/PRES(JLON,1)
      PRPP  (JLON,JLEV)=(PLNPR(JLON,JLEV)*PRDELP(JLON,JLEV))**2
    ENDDO
  ENDDO

ENDIF

! (PLNPR(JLON,1) AND PRPP(JLON,1) ARE A PRIORI NOT USED LATER)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPXYB',1,ZHOOK_HANDLE)
END SUBROUTINE GPXYB
