SUBROUTINE CTSTAR(KPROMA,KSTART,KPROF,PTB,PRESBH,PRESBF,POROG,PTSTAR,PT0)

!**** *CTSTAR* - COMPUTES STANDARD SURFACE TEMPERATURE
!                              AND SURFACE TEMPERATURE.

!     PURPOSE.
!     --------

!           COMPUTES THE STANDARD SURFACE TEMPERATURE AND THE SURFACE
!           TEMPERATURE TO BE USED FOR EXTRAPOLATIONS OF TEMPERATURE
!           AND GEOPOTENTIEL.

!**   INTERFACE.
!     ----------
!        *CALL* *CTSTAR(..)*

!        EXPLICIT ARGUMENTS
!        --------------------

!        KPROMA         - HORIZONTAL DIMENSIONS.             (INPUT)
!        KSTART         - START OF WORK                      (INPUT)
!        KPROF          - DEPTH OF WORK                      (INPUT)

!        PTB(KPROMA)    - TEMPERATURE AT NFLEVG-1             (INPUT)
!        PRESBH(KPROMA) - LOWEST MODEL HALF LEVEL PRESSURES  (INPUT)

!        PRESBF(KPROMA) - PRESSURE AT NFLEVG-1                (INPUT)
!        POROG(KPROMA)  - MODEL ORGRAPHY                     (INPUT)


!        PTSTAR(KPROMA) - SURFACE TEMPERATURE                (OUTPUT)

!        PT0(KPROMA)    - STANDARD SURFACE TEMPERATURE       (OUTPUT)

!        IMPLICIT ARGUMENTS :    CONSTANTS FROM YOMSTA,YOMCST.
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.   NONE.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-05-02

!      Modification : 93-06-01 M.Hamrud (Comment only, now T from NFLEVG-1)
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

!USE PARKIND1 
!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/xrd/module/parkind1.F90.php#parkind1>  ,ONLY : JPIM     ,JPRB
!USE YOMHOOK 
!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/xrd/module/yomhook.F90.php#yomhook>   ,ONLY : LHOOK,   DR_HOOK

!USE YOMCST, ONLY : RG, RD 
!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/arp/module/yomcst.F90.php#yomcst>   , ONLY :  RG

!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/xrd/support/rg.F.php#rg>       ,RD
!USE YOMSTA 
!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/arp/module/yomsta.F90.php#yomsta>   , ONLY : RDTDZ1

IMPLICIT NONE

include "YOMCST.h"
!IM INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
!IM INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART
!IM INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
INTEGER,INTENT(IN)    :: KPROMA
INTEGER,INTENT(IN)    :: KSTART
INTEGER,INTENT(IN)    :: KPROF
!IM REAL(KIND=JPRB)   ,INTENT(IN)    :: PTB(KPROMA)
REAL   ,INTENT(IN)    :: PTB(KPROMA)
!IM REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESBH(KPROMA)
REAL   ,INTENT(IN)    :: PRESBH(KPROMA)
!IM REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESBF(KPROMA)
REAL   ,INTENT(IN)    :: PRESBF(KPROMA)
!IM REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KPROMA)
REAL   ,INTENT(IN)    :: POROG(KPROMA)
!IM REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSTAR(KPROMA)
REAL   ,INTENT(OUT)   :: PTSTAR(KPROMA)
!IM REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT0(KPROMA)
REAL   ,INTENT(OUT)   :: PT0(KPROMA)
!IM INTEGER(KIND=JPIM) :: JL
INTEGER :: JL

!IM REAL(KIND=JPRB) :: ZALPHA, ZDTDZSG
REAL :: ZALPHA, ZDTDZSG
!IM REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL :: ZHOOK_HANDLE
!IM beg
REAL, PARAMETER                  :: RDTDZ1=-0.0065 !or USE YOMSTA
!IM end

!     ------------------------------------------------------------------

!*       1.    COMPUTES SURFACE TEMPERATURE
!*             THEN STANDARD SURFACE TEMPERATURE.

!IF (LHOOK) CALL DR_HOOK('CTSTAR',0,ZHOOK_HANDLE)
ZDTDZSG=-RDTDZ1/RG 
!<http://intra.cnrm.meteo.fr/eac/ARPCLI5.2/doci/code/arpcli5.2/xrd/support/rg.F.php#rg>
ZALPHA=ZDTDZSG*RD
DO JL=KSTART,KPROF

   !IM PTSTAR(JL)=PTB(JL)*(1.0_JPRB+ZALPHA*(PRESBH(JL)/PRESBF(JL)-1.0_JPRB))
   PTSTAR(JL)=PTB(JL)*(1.0+ZALPHA*(PRESBH(JL)/PRESBF(JL)-1.0))
   PT0(JL)=PTSTAR(JL)+ZDTDZSG*POROG(JL)
!  print*,'cstar JL ptb zalpha PRESBH PRESBF ptstar' &
!  ,JL,PTB(JL),ZALPHA,PRESBH(JL),PRESBF(JL),PTSTAR(JL)
ENDDO


!     ------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('CTSTAR',1,ZHOOK_HANDLE)
END SUBROUTINE CTSTAR
