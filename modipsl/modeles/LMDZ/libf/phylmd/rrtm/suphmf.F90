!OPTIONS XOPT(NOEVAL)
SUBROUTINE SUPHMF(KULOUT)

!**** *SUPHMF*   - Calls initialization of commons controlling physics
!                  in the Meteo-France version.

!     Purpose.
!     --------
!           Organise the setup of physical constants for Meteo-France
!           physics package.

!**   Interface.
!     ----------
!        *CALL* *SUPHMF(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        None.

!     Method.
!     -------
!        Irrelevant.

!     Externals.
!     ----------

!     SUPHY0
!     SUPHY1
!     SUPHY2
!     SUPHY3
!     SUTOPH

!     Reference.
!     ----------

!     Author.
!     -------
!        J.-F. Geleyn.

!     Modifications.
!     --------------
!        Original : 91-06-15
!        Modified 91-06-10 by A. Lasserre-Bigorry (call to SUTOPH)
!        Modified 99-03-01 by D. Giard (call to VAL923 for 923 and 927)
!        Modified 01-04-02 R. El Khatib setup for CAPE diagnostic
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Modified 04-11-16 Y. Seity : call suphmnh for AROME physics
!        R. Zaaboul 28-Feb-2006: call suparar, suphmpa and suphmse (ex suphmnh)
!        Y. Seity 06-07-10: nfpsurfex and lfpart2 in call suphmse (prepsurfex)  
!     ------------------------------------------------------------------

!*       1.    Call routines for specific physics' commons setup.
!              --------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMPHY   , ONLY : LSOLV
USE YOMARPHY , ONLY : LMPA, LMSE
USE YOMCT0   , ONLY : LSFORC, LFPART2
USE YOMFPC   , ONLY : NFPSURFEX

! Ce qui concerne NULNAM et JPNULNAM commente par MPL le 15.04.09
!USE PARDIM   , ONLY : JPNULNAM
!USE YOMLUN   , ONLY : NULNAM

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "sucape.intfb.h"
#include "su0phy.intfb.h"
#include "suphy0.intfb.h"
#include "suphy1.intfb.h"
#include "suphy2.intfb.h"
#include "suphy3.intfb.h"
#include "sutoph.intfb.h"
#include "val923.intfb.h"
! commente par MPL 20.11.08
!#include "suparar.intfb.h"
!#include "suphmpa.intfb.h"
!#include "suphmse.intfb.h"

! Ce qui concerne MNAM commente par MPL le 15.04.09
!NULNAM = JPNULNAM
!OPEN(NULNAM,ACTION="READ")
!OPEN(NULNAM,FILE='MNAM',ACTION="READ")

IF (LHOOK) CALL DR_HOOK('SUPHMF',0,ZHOOK_HANDLE)
print *,'SUPHMF: avant SU0PHY'
CALL SU0PHY(KULOUT)
print *,'SUPHMF: avant SUPHY0'
CALL SUPHY0(KULOUT)
print *,'SUPHMF: avant SUPHY1'
CALL SUPHY1(KULOUT)
print *,'SUPHMF: avant SUPHY2'
CALL SUPHY2(KULOUT)
print *,'SUPHMF: avant SUPHY3'
CALL SUPHY3(KULOUT)
print *,'SUPHMF: avant SUTOPH'
CALL SUTOPH(KULOUT)
print *,'SUPHMF: avant VAL923'

CALL VAL923(LSOLV)

print *,'SUPHMF: avant SUCAPE'
CALL SUCAPE(KULOUT)

! setup for AROME physics and SURFEX
! commente par MPL 20.11.08
!CALL SUPARAR(KULOUT)
!IF (LMPA) CALL SUPHMPA(KULOUT)
!IF (LMSE.AND.NFPSURFEX==0.AND.(.NOT.LFPART2)) CALL SUPHMSE(KULOUT)

!     ------------------------------------------------------------------


! Ce qui concerne NULNAM commente par MPL le 15.04.09
!CLOSE(NULNAM)
IF (LHOOK) CALL DR_HOOK('SUPHMF',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHMF
