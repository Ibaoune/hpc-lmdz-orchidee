SUBROUTINE GPPRE(KPROMA,KSTART,KPROF,KFLEV,PVAH,PVBH,PRESH,PRESF)

!**** *GPPRE* - Computes half and full level pressure

!     Purpose.
!     --------
!           Computes pressures at half and full model levels.

!**   Interface.
!     ----------
!        *CALL* *GPPRE(KPROMA,KSTART,KPROF,KFLEV,PVAH,PVBH,PRESH,PRESF)

!        Explicit arguments :
!        --------------------
!                              KPROMA :  dimensioning
!                              KSTART :  start of work
!                              KPROF  :  depth of work
!                              KFLEV     : vert. dimensioning
!                              PRESH(KPROMA,0:KFLEV) - HALF LEVEL PRESSURE
!                              PRESF(KPROMA,KFLEV)   - FULL LEVEL PRESSURE
!                              PVAH(KFLEV),PVBH(KFLEV)- vertical coordinate
!        Implicit arguments :  NONE.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!     Half level P:     PHk = Ak + Bk * Ps

!                                PHk*ln(PHk) - PHk-1*ln(PHk-1)
!     Full level P: ln(PFk) = [ ------------------------------- - 1. ]
!                                        PHk - PHk-1

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-02-04
!     Erik Andersson  920326: Altered computation of full level pressure
!     Erik Andersson  930225: Use GPPREH/F.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMGEM   , ONLY : VC       ,VDELB

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVAH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVBH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRESH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRESF(KPROMA,KFLEV) 
REAL(KIND=JPRB) ::    ZLNPR(KPROMA,KFLEV) ,  ZALPH (KPROMA,KFLEV)
REAL(KIND=JPRB) ::    ZDELP(KPROMA,KFLEV) ,  ZRDELP(KPROMA,KFLEV)
REAL(KIND=JPRB) ::    ZRTGR(KPROMA,KFLEV) ,  ZRPRES(KPROMA,KFLEV)
REAL(KIND=JPRB) ::    ZRPP (KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "gppref.intfb.h"
#include "gppreh.intfb.h"
#include "gpxyb.intfb.h"

!     ------------------------------------------------------------------

!*       1.    COMPUTES HALF AND FULL LEVEL PRESSURES
!              --------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPRE',0,ZHOOK_HANDLE)
CALL GPPREH(KPROMA,KSTART,KPROF,KFLEV,PVAH,PVBH,PRESH)
CALL GPXYB(KPROMA,KSTART,KPROF,KFLEV,VDELB,VC,PRESH,ZDELP,&
 & ZRDELP,ZLNPR,ZALPH,ZRTGR,ZRPRES,ZRPP)  
CALL GPPREF(KPROMA,KSTART,KPROF,KFLEV,PVAH,PVBH,ZALPH,PRESH,PRESF)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPRE',1,ZHOOK_HANDLE)
END SUBROUTINE GPPRE
