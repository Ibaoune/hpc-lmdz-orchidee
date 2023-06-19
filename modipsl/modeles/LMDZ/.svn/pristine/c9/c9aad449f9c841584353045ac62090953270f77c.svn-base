SUBROUTINE LWBV &
 & ( KIDIA, KFDIA, KLON , KLEV  , KMODE,&
 & PDT0 , PEMIS, PEMIW,&
 & PTL  , PTAVE,&
 & PEMIT, PFLUC,&
 & PABCU, PBINT, PBSUI, PCNTRB &
 & )  

!**** *LWBV*   - COMPUTE PLANCK FUNC., PERF. VERT. INTEGRATION

!     PURPOSE.
!     --------
!           TO COMPUTE THE PLANCK FUNCTION AND PERFORM THE
!           VERTICAL INTEGRATION. SPLIT OUT FROM LW FOR MEMORY
!           SAVING

!**   INTERFACE.
!     ----------

!        *LWVB* IS CALLED FROM *LW*

!        EXPLICIT ARGUMENTS :
!        --------------------
! PDT0   : (KLON)            ; SURFACE TEMPERATURE DISCONTINUITY
! PEMIS  : (KLON)            ; SURFACE LW EMISSIVITY
! PEMIW  : (KLON)            ; SURFACE LW WINDOW EMISSIVITY
! PTAVE  : (KLON,KLEV)       ; TEMPERATURE
! PTL    : (KLON,KLEV+1)     ; HALF LEVEL TEMPERATURE
!     ==== OUTPUTS ===
! PABCU  :
! PBINT  :
! PBSUI  :
! PCNTRB : 
! PCOLC  :
! PEMIT  :
! PFLUC  :

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
!     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
!          2. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
!     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
!     BOUNDARIES.
!          3. COMPUTES THE CLEAR-SKY COOLING RATES.

!     EXTERNALS.
!     ----------

!          *LWB*, *LWV*

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
!        MODIFICATION : 93-10-15 M.HAMRUD (SPLIT OUT FROM LW TO SAVE
!                                          MEMORY)
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOELW    , ONLY : NSIL     ,NIPD     ,NUA
USE YOERDU   , ONLY : NUAER    ,NTRAER

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMODE 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT0(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAVE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMIT(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUC(KLON,2,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PABCU(KLON,NUA,3*KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBINT(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBSUI(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCNTRB(KLON,KLEV+1,KLEV+1) 
!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-------------------------------------------------------------------------

!              ------------
REAL(KIND=JPRB) ::&
 & ZB(KLON,NSIL,KLEV+1), ZBSUR(KLON,NSIL)   , ZBTOP(KLON,NSIL)&
 & ,  ZDBSL(KLON,NSIL,KLEV*2)&
 & ,  ZGA(KLON,NIPD,2,KLEV)  , ZGB(KLON,NIPD,2,KLEV)&
 & ,  ZGASUR(KLON,NIPD,2)    , ZGBSUR(KLON,NIPD,2)&
 & ,  ZGATOP(KLON,NIPD,2)    , ZGBTOP(KLON,NIPD,2)  

INTEGER(KIND=JPIM) :: JL, JLW
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "lwb.intfb.h"
#include "lwv.intfb.h"

!     ------------------------------------------------------------------

!*         2.    COMPUTES PLANCK FUNCTIONS
!                -------------------------

IF (LHOOK) CALL DR_HOOK('LWBV',0,ZHOOK_HANDLE)
print *,'LWBV: avant LWB'
CALL LWB &
 & ( KIDIA, KFDIA, KLON  , KLEV  , KMODE,&
 & PDT0 , PTAVE, PTL,&
 & ZB   , PBINT, ZBSUR , ZBTOP , ZDBSL,&
 & ZGA  , ZGB  , ZGASUR, ZGBSUR, ZGATOP, ZGBTOP    &
 & )  

!     ------------------------------------------------------------------

!*         3.    PERFORMS THE VERTICAL INTEGRATION
!                ---------------------------------

CALL LWV &
 & ( KIDIA , KFDIA, KLON  , KLEV,&
 & NUAER , NTRAER,&
 & PABCU , ZB   , PBINT , ZBSUR , ZBTOP , ZDBSL,&
 & PEMIS , PEMIW,&
 & ZGA   , ZGB  , ZGASUR, ZGBSUR, ZGATOP, ZGBTOP,&
 & PCNTRB, PFLUC &
 & )  

DO JL=KIDIA,KFDIA
  PEMIT(JL)=0.0_JPRB
  PBSUI(JL)=0.0_JPRB
ENDDO
DO JLW=1,NSIL
  DO JL=KIDIA,KFDIA
    PBSUI(JL)=PBSUI(JL)+ZBSUR(JL,JLW)
    IF (JLW >= 3.AND. JLW <= 4) THEN
      PEMIT(JL)=PEMIT(JL)+ZBSUR(JL,JLW)*PEMIW(JL)
    ELSE
      PEMIT(JL)=PEMIT(JL)+ZBSUR(JL,JLW)*PEMIS(JL)
    ENDIF
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
  PEMIT(JL)=PEMIT(JL)/PBSUI(JL)
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LWBV',1,ZHOOK_HANDLE)
END SUBROUTINE LWBV
