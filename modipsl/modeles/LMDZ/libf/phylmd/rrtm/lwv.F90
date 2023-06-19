SUBROUTINE LWV &
 & ( KIDIA, KFDIA, KLON , KLEV , KUAER , KTRAER,&
 & PABCU, PB   , PBINT, PBSUR, PBTOP , PDBSL,&
 & PEMIS, PEMIW,&
 & PGA  , PGB  , PGASUR,PGBSUR,PGATOP, PGBTOP,&
 & PCNTRB,PFLUC &
 & )  

!**** *LWV*   - LONGWAVE RADIATION, VERTICAL INTEGRATION

!     PURPOSE.
!     --------
!           CARRIES OUT THE VERTICAL INTEGRATION TO GIVE LONGWAVE
!           FLUXES OR RADIANCES

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
! PB     : (KLON,NSIL,KLEV+1); SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
! PBINT  : (KLON,KLEV+1)     ; HALF-LEVEL PLANCK FUNCTIONS
! PBSUR  : (KLON,NSIL)       ; SURFACE SPECTRAL PLANCK FUNCTION
! PBTOP  : (KLON,NSIL)       ; T.O.A. SPECTRAL PLANCK FUNCTION
! PDBSL  : (KLON,KLEV*2)     ; SUB-LAYER PLANCK FUNCTION GRADIENT
! PEMIS  : (KLON)            ; SURFACE LW EMISSIVITY
! PEMIW  : (KLON)            ; SURFACE LW WINDOW EMISSIVITY
! PGA, PGB                   ; PADE APPROXIMANTS
! PGASUR, PGBSUR             ; SURFACE PADE APPROXIMANTS
! PGATOP, PGBTOP             ; T.O.A. PADE APPROXIMANTS
!     ==== OUTPUTS ===
! PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
! PFLUC(KLON,2,KLEV)           ; RADIATIVE FLUXES CLEAR-SKY

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING BETWEEN
!     CONTRIBUTIONS BY -  THE NEARBY LAYERS
!                      -  THE DISTANT LAYERS
!                      -  THE BOUNDARY TERMS
!          2. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.

!     EXTERNALS.
!     ----------

!          *LWVN*, *LWVD*, *LWVB*

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
!        JJ Morcrette 96-06-07 Surface LW window emissivity
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOELW    , ONLY : NSIL     ,NIPD     ,NUA

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KUAER 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAER 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PABCU(KLON,NUA,3*KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB(KLON,NSIL,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBINT(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBSUR(KLON,NSIL) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBTOP(KLON,NSIL) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDBSL(KLON,NSIL,KLEV*2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGA(KLON,NIPD,2,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGB(KLON,NIPD,2,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGASUR(KLON,NIPD,2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGBSUR(KLON,NIPD,2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGATOP(KLON,NIPD,2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGBTOP(KLON,NIPD,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCNTRB(KLON,KLEV+1,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUC(KLON,2,KLEV+1) 
!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-----------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZADJD(KLON,KLEV+1)  , ZADJU(KLON,KLEV+1)&
 & ,  ZDBDT(KLON,NSIL,KLEV)&
 & ,  ZDISD(KLON,KLEV+1)  , ZDISU(KLON,KLEV+1)&
 & ,  ZDWFSU(KLON,NSIL)  

INTEGER(KIND=JPIM) :: JA, JK, JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "lwvb.intfb.h"
#include "lwvd.intfb.h"
#include "lwvn.intfb.h"

!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!*         1.1     INITIALIZE LAYER CONTRIBUTIONS
!                  ------------------------------

IF (LHOOK) CALL DR_HOOK('LWV',0,ZHOOK_HANDLE)
DO JK=1,KLEV+1
  DO JL=KIDIA,KFDIA
    ZADJD(JL,JK)=0.0_JPRB
    ZADJU(JL,JK)=0.0_JPRB
    ZDISD(JL,JK)=0.0_JPRB
    ZDISU(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
DO JA=1,NSIL
  DO JL=KIDIA,KFDIA
    ZDWFSU(JL,JA)=0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         2.      VERTICAL INTEGRATION
!                  --------------------

!     ------------------------------------------------------------------

!*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
!                  ---------------------------------

CALL LWVN &
 & ( KIDIA, KFDIA, KLON  , KLEV , KUAER,&
 & PABCU, PDBSL, PGA   , PGB,&
 & ZADJD, ZADJU, PCNTRB, ZDBDT, ZDWFSU  &
 & )  

!     ------------------------------------------------------------------

!*         2.2     CONTRIBUTION FROM DISTANT LAYERS
!                  ---------------------------------

CALL LWVD &
 & ( KIDIA , KFDIA, KLON , KLEV  , KTRAER,&
 & PABCU , ZDBDT, PGA  , PGB,&
 & PCNTRB, ZDISD, ZDISU, ZDWFSU &
 & )  

!     ------------------------------------------------------------------

!*         2.3     EXCHANGE WITH THE BOUNDARIES
!                  ----------------------------

CALL LWVB &
 & ( KIDIA , KFDIA , KLON  , KLEV  , KUAER,&
 & PABCU , ZADJD , ZADJU,&
 & PB    , PBINT , PBSUR , PBTOP,&
 & ZDISD , ZDISU , PEMIS , PEMIW,&
 & PGASUR, PGBSUR, PGATOP, PGBTOP,&
 & ZDWFSU,PFLUC  &
 & )  

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LWV',1,ZHOOK_HANDLE)
END SUBROUTINE LWV
