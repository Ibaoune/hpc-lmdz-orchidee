SUBROUTINE LW &
 & ( KIDIA, KFDIA , KLON  , KLEV  , KMODE,&
 & PCCO2, PCLDLD, PCLDLU,&
 & PDP  , PDT0  , PEMIS , PEMIW,&
 & PPMB , PQOF  , PTL,&
 & PAER , PTAVE , PVIEW , PWV,&
 & PEMIT, PFLUX , PFLUC &
 & )  

!**** *LW*   - ORGANIZES THE LONGWAVE CALCULATIONS

!     PURPOSE.
!     --------
!           COMPUTES LONGWAVE FLUXES 

!**   INTERFACE.
!     ----------

!        *LW* IS CALLED FROM *RADLSW*

!        EXPLICIT ARGUMENTS :
!        --------------------
! PAER   : (KLON,6,KLEV)     ; OPTICAL THICKNESS OF THE AEROSOLS
! PCCO2  :                   ; CONCENTRATION IN CO2 (PA/PA)
! PCLDLD : (KLON,KLEV)       ; DOWNWARD EFFECTIVE CLOUD FRACTION
! PCLDLU : (KLON,KLEV)       ; UPWARD EFFECTIVE CLOUD FRACTION
! PDP    : (KLON,KLEV)       ; LAYER PRESSURE THICKNESS
! PDT0   : (KLON)            ; SURFACE TEMPERATURE DISCONTINUITY  
! PEMIS  : (KLON)            ; SURFACE LW EMISSIVITY
! PEMIW  : (KLON)            ; SURFACE LW WINDOW EMISSIVITY
! PPMB   : (KLON,KLEV+1)     ; HALF LEVEL PRESSURE
! PQOF   : (KLON,KLEV)       ; CONCENTRATION IN OZONE (PA/PA)
! PTAVE  : (KLON,KLEV)       ; TEMPERATURE
! PTL    : (KLON,KLEV+1)     ; HALF LEVEL TEMPERATURE
! PVIEW  : (KLON)            ; COSECANT OF VIEWING ANGLE
! PWV    : (KLON,KLEV)       ; SPECIFIC HUMIDITY  (PA/PA)
!     ==== OUTPUTS ===
! PEMIT(KLON)                ; SURFACE TOTAL LW EMISSIVITY
! PFLUX(KLON,2,KLEV+1)       ; RADIATIVE FLUXES :
!                     1  ==>  UPWARD   FLUX TOTAL
!                     2  ==>  DOWNWARD FLUX TOTAL
! PFLUC(KLON,2,KLEV+1)       ; RADIATIVE FLUXES CLEAR SKY:
!                     1  ==>  UPWARD   FLUX TOTAL
!                     2  ==>  DOWNWARD FLUX TOTAL

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
!     ABSORBERS.
!          2. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
!     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
!          3. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
!     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
!     BOUNDARIES.
!          4. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
!          5. INTRODUCES THE EFFECTS OF THE CLOUDS ON THE FLUXES.

!     EXTERNALS.
!     ----------

!          *LWU*, *LWBV*, *LWC*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        99-05-25   JJMorcrette    Revised aerosols
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOELW    , ONLY : NUA
IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMODE 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCO2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDLD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT0(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPMB(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQOF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAVE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVIEW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEMIT(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUX(KLON,2,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFLUC(KLON,2,KLEV+1) 
!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-------------------------------------------------------------------------

!              ------------
REAL(KIND=JPRB) :: ZABCU(KLON,NUA,3*KLEV+1)&
 & ,  ZBINT(KLON,KLEV+1)        , ZBSUI(KLON)&
 & ,  ZCNTRB(KLON,KLEV+1,KLEV+1)  

REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "lwbv.intfb.h"
#include "lwc.intfb.h"
#include "lwu.intfb.h"

!     ------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!100  CONTINUE

!     ------------------------------------------------------------------

!*         1.1   COMPUTES ABSORBER AMOUNTS
!                -------------------------

IF (LHOOK) CALL DR_HOOK('LW',0,ZHOOK_HANDLE)
print *,'	LW: Avant LWU'
CALL LWU &
 & (  KIDIA, KFDIA, KLON, KLEV,&
 & PAER , PCCO2, PDP , PPMB, PQOF , PTAVE, PVIEW, PWV,&
 & ZABCU &
 & )  

!     ------------------------------------------------------------------

!*         2.    COMPUTES PLANCK FUNCTIONS
!                -------------------------
!                PERFORMS THE VERTICAL INTEGRATION
!                ---------------------------------

print *,'	LW: Avant LWBV'
CALL LWBV &
 & ( KIDIA, KFDIA, KLON , KLEV  , KMODE,&
 & PDT0 , PEMIS, PEMIW, PTL   , PTAVE,&
 & PEMIT, PFLUC,&
 & ZABCU, ZBINT, ZBSUI, ZCNTRB &
 & )  

!     ------------------------------------------------------------------

!*         4.    INTRODUCES THE EFFECTS OF CLOUDS
!                --------------------------------

print *,'	LW: Avant LWC'
CALL LWC &
 & ( KIDIA , KFDIA, KLON  , KLEV,&
 & ZBINT , ZBSUI, PCLDLD, PCLDLU,&
 & ZCNTRB, PEMIT, PFLUC,&
 & PFLUX    &
 & )  

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LW',1,ZHOOK_HANDLE)
END SUBROUTINE LW
