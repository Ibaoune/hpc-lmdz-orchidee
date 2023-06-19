SUBROUTINE PPPMER(KPROMA,KSTART,KPROF,PRPRESS,POROG,PTSTAR,PT0,PMSLPPP)

!**** *PPPMER* - POST-PROCESS MSL PRESSURE.

!     PURPOSE.
!     --------
!           COMPUTES MSL PRESSURE.

!**   INTERFACE.
!     ----------

!        *CALL* *PPPMER(KPROMA,KSTART,KPROF,PRPRESS,POROG,PTSTAR,PT0,
!    S                  PMSLPPP)

!        EXPLICIT ARGUMENTS
!        --------------------


!        KPROMA                    - HORIZONTAL DIMENSION.             (INPUT)
!        KSTART                    - START OF WORK.                    (INPUT)
!        KPROF                     - DEPTH OF WORK.                    (INPUT)
!        PRPRESS(KPROMA)           - SURFACE PRESSURE                  (INPUT)
!        POROG(KPROMA)             - MODEL OROGRAPHY.                  (INPUT)
!        PTSTAR(KPROMA)            - SURFACE TEMPERATURE               (INPUT)
!        PT0(KPROMA)               - STANDARD SURFACE TEMPERATURE      (INPUT)
!        PMSLPPP(KPROMA)           - POST-PROCESSED MSL PRESSURE       (OUTPUT)
!        IMPLICIT ARGUMENTS :  CONSTANTS FROM YOMCST,YOMGEM,YOMSTA.
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.  NONE
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-01-26

!     E. A-son, J-F Geleyn 920409 Mod. T*, T0 and alpha below surface.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

! USE PARKIND1 
!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/xrd/module/parkind1.F90.php#parkind1>  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK 
!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/xrd/module/yomhook.F90.php#yomhook>   ,ONLY : LHOOK,   DR_HOOK

!USE YOMCST, ONLY : RG, RD
!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/arp/module/yomcst.F90.php#yomcst>   , ONLY : RG

!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/xrd/support/rg.F.php#rg>       ,RD
! USE YOMSTA 
!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/arp/module/yomsta.F90.php#yomsta>   , ONLY : RDTDZ1

  IMPLICIT NONE

include "YOMCST.h"
!IM INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
!IM INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART
!IM INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
 INTEGER,INTENT(IN)    :: KPROMA
 INTEGER,INTENT(IN)    :: KSTART
 INTEGER,INTENT(IN)    :: KPROF
!IM REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESS(KPROMA)
!IM REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KPROMA)
!IM REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTAR(KPROMA)
!IM REAL(KIND=JPRB)   ,INTENT(IN)    :: PT0(KPROMA)
!IM REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMSLPPP(KPROMA)
!IM REAL(KIND=JPRB) :: ZTSTAR(KPROMA)
!IM REAL(KIND=JPRB) :: ZALPHA(KPROMA)
 REAL,INTENT(IN)    :: PRPRESS(KPROMA)
 REAL,INTENT(IN)    :: POROG(KPROMA)
 REAL,INTENT(IN)    :: PTSTAR(KPROMA)
 REAL,INTENT(IN)    :: PT0(KPROMA)
 REAL,INTENT(OUT)   :: PMSLPPP(KPROMA)
 REAL :: ZTSTAR(KPROMA)
 REAL :: ZALPHA(KPROMA)

!IM INTEGER(KIND=JPIM) :: JL
 INTEGER :: JL

!IM REAL(KIND=JPRB) :: ZDTDZSG, ZOROG, ZT0, ZTX, ZTY, ZX, ZY, ZY2
!IM REAL(KIND=JPRB) :: ZHOOK_HANDLE
 REAL :: ZDTDZSG, ZOROG, ZT0, ZTX, ZTY, ZX, ZY, ZY2
 REAL :: ZHOOK_HANDLE
!IM beg
REAL, PARAMETER                  :: RDTDZ1=-0.0065 !or USE YOMSTA
!IM end

!     ------------------------------------------------------------------

!*       1.    POST-PROCESS MSL PRESSURE.
!              --------------------------

!*       1.1   COMPUTATION OF MODIFIED ALPHA AND TSTAR.

!IM IF (LHOOK) CALL DR_HOOK('PPPMER',0,ZHOOK_HANDLE)
!IM ZTX=290.5_JPRB
!IM ZTY=255.0_JPRB
 ZTX=290.5
 ZTY=255.0
 ZDTDZSG=-RDTDZ1/RG 
!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/xrd/support/rg.F.php#rg>
 DO JL=KSTART,KPROF

   IF(PTSTAR(JL) < ZTY) THEN
!IM  ZTSTAR(JL)=0.5_JPRB*(ZTY+PTSTAR(JL))
     ZTSTAR(JL)=0.5*(ZTY+PTSTAR(JL))
   ELSEIF(PTSTAR(JL) < ZTX) THEN
     ZTSTAR(JL)=PTSTAR(JL)
   ELSE
!IM    ZTSTAR(JL)=0.5_JPRB*(ZTX+PTSTAR(JL))
     ZTSTAR(JL)=0.5*(ZTX+PTSTAR(JL))
   ENDIF

   ZT0=ZTSTAR(JL)+ZDTDZSG*POROG(JL)
   IF(ZTX > ZTSTAR(JL) .AND. ZT0 > ZTX) THEN
     ZT0=ZTX
   ELSEIF(ZTX <= ZTSTAR(JL) .AND. ZT0 > ZTSTAR(JL)) THEN
     ZT0=ZTSTAR(JL)
   ELSE
     ZT0=PT0(JL)
   ENDIF

!IM  ZOROG=SIGN(MAX(1.0_JPRB,ABS(POROG(JL))),POROG(JL))
   ZOROG=SIGN(MAX(1.0,ABS(POROG(JL))),POROG(JL))
   ZALPHA(JL)=RD*(ZT0-ZTSTAR(JL))/ZOROG
 ENDDO

!*       1.2   COMPUTATION OF MSL PRESSURE.

 DO JL=KSTART,KPROF
!IM  IF (ABS(POROG(JL)) >= 0.001_JPRB) THEN
   IF (ABS(POROG(JL)) >= 0.001) THEN
     ZX=POROG(JL)/(RD*ZTSTAR(JL))
     ZY=ZALPHA(JL)*ZX
     ZY2=ZY*ZY

!IM    PMSLPPP(JL)=PRPRESS(JL)*EXP(ZX*(1.0_JPRB-0.5_JPRB*ZY+1.0_JPRB/3._JPRB*ZY2))
     PMSLPPP(JL)=PRPRESS(JL)*EXP(ZX*(1.0-0.5*ZY+1.0/3.*ZY2))
   ELSE
     PMSLPPP(JL)=PRPRESS(JL)
   ENDIF
 ENDDO


!     ------------------------------------------------------------------

!IM IF (LHOOK) CALL DR_HOOK('PPPMER',1,ZHOOK_HANDLE)
 END SUBROUTINE PPPMER
