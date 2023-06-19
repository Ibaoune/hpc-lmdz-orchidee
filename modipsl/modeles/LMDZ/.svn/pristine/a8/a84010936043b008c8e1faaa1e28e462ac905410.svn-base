SUBROUTINE SUOVLP ( KLEV )
  
!***** *SUOVLP*   -INITIALIZE PROFILE OF ALPHA1

!**   INTERFACE.
!     ----------
!        CALL *SUOVLP* FROM *SUECRAD*
!              ------        -------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 01-02-16
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOERAD   , ONLY : RAOVLP, RBOVLP
USE YOEOVLP  , ONLY : RA1OVLP
USE YOMSTA   , ONLY : STZ

IMPLICIT NONE

!     ------------------------------------------------------------------

!     DUMMY PARAMETERS

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 

INTEGER(KIND=JPIM) :: JK
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUOVLP',0,ZHOOK_HANDLE)
DO JK=1,KLEV
  RA1OVLP(JK)=RAOVLP*STZ(JK)+RBOVLP
  print *,'SU_OVLP: JK RAOVLP STZ RBOVLP:',JK,RAOVLP,STZ(JK),RBOVLP
ENDDO  

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUOVLP',1,ZHOOK_HANDLE)
END SUBROUTINE SUOVLP
