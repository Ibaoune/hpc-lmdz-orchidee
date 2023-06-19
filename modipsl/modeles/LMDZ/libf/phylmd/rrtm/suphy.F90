SUBROUTINE SUPHY(KULOUT)

USE YOMCT0   , ONLY : NCONF
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!**** *SUPHY*   - Calls physic specific set-up routines

!     Purpose.
!     --------
!           Calls set-up routines specific to the different physics
!           packages that can be used in the IFS/ARPEGE model

!**   Interface.
!     ----------
!        *CALL* *SUPHY(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY, YOEPHY

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     SUPHEC
!     SUPHMF
!     SUPHYFL
!     SUHLPH
!
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!        or

!        Documentation ARPEGE (depending on which physics will be used)

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*
!        J.-F. Geleyn for the ARPEGE rewriting.

!     Modifications.
!     --------------
!        Original : 87-10-15
!        ARPEGE extension 90-9-1
!        ARPEGE modification 90-11-24
!        Modified by M. Deque 91-02-28 (key for Ozone)
!        Modified by M. Deque 91-03-18 (key for Negative humidity)
!        Modified by M. Deque 91-04-01 (keys for cloudiness and wavedrag)
!        Modified by J.-F. Geleyn 91-06-15 (call to SUPHMF and SUPHEC)
!        Modified by J.-J. Morcrette 91-11-12
!        Modified by K. YESSAD (MAY 2000): remove call to EC physics setup
!         in a 2D model because some dimensionings are inconsistent and
!         generate aborts, and because 2D models are purely adiabatic ones.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        B.Sass        01-June-2006 (call setup for HIRLAM physics)
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
LOGICAL :: LLSHW, LLVEQ, LL2D
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "suphec.intfb.h"
#include "suphmf.intfb.h"
#include "suhlph.intfb.h"
#include "suphyfl.intfb.h"
#include "sumts.intfb.h"

!     ------------------------------------------------------------------

!*       0.    Set-up LL2D (2D model switch).
!              ------------------------------

IF (LHOOK) CALL DR_HOOK('SUPHY',0,ZHOOK_HANDLE)
LLSHW=NCONF == 201.OR.NCONF == 421.OR.NCONF == 521
LLVEQ=NCONF == 202.OR.NCONF == 422.OR.NCONF == 522
LL2D=LLSHW.OR.LLVEQ

!     ------------------------------------------------------------------

!*       1.    Call initialization of specific physics' commons.
!              -------------------------------------------------

!*       1.1   Meteo-France Physics
!              --------------------

print *,'---- SUPHY: avant SUPHMF'
CALL SUPHMF(KULOUT)
!
print *,'---- SUPHY: avant SUGFL'
!SUGFL: Set up unified_treatment grid-point fields
CALL SUGFL


!*       1.2   ECMWF Physics
!              -------------

! IF Commente par MPL 20.11.08
!IF (.NOT.LL2D) THEN
   CALL SUPHEC(KULOUT)
!ENDIF

!        1.3   Initialize HIRLAM physics
!              -------------------------
! Commente par MPL 20.11.08
!CALL SUHLPH(KULOUT)    
!
!     ------------------------------------------------------------------

!*       2.    Initialize physics' flags commons.
!              ----------------------------------

! Commente par MPL 20.11.08
!CALL SUPHYFL

!     ------------------------------------------------------------------

!*       3.    Initialize "model to satellite" RTTOV parameters.
!              ------------------------------------------------

! Commente par MPL 20.11.08
!CALL SUMTS

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPHY',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHY
