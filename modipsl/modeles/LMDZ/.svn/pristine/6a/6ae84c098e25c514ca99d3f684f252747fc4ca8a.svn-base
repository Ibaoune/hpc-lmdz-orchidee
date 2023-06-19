SUBROUTINE SURAYOLMD
#ifdef DOC

!     **** *SURAYOLMD - initialization rayonnement LMD

!     Purpose.
!     --------

!     **   Interface.
!     ----------

!     Explicit arguments :    None.
!     --------------------

!     Implicit arguments :    None.
!     --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------

!     Modifications.
!     --------------
!     Original :08-04-03
!     MP Lefebvre modified suinit.F90 program.
!     ------------------------------------------------------------------
#endif

USE PARKIND1  ,ONLY : JPIM     ,JPRB
!#include "tsmbkind.h"

!USE YOMFORC  , ONLY : CFORC
USE YOMCT0B  , ONLY : LECMWF
!USE YOMLUN   , ONLY : NULOUT
!USE YOMRIP   , ONLY : NINDAT   ,NSSSSS
USE YOMDIM   , ONLY : NFLEVG
!USE YOMINI   , ONLY : NEINI
!USE YOMVRTL  , ONLY : L131TL
!USE YOMMCC   , ONLY : LMCC01, LMCC02, LMCC03, LMCCEC
!USE YOMNMIA  , ONLY : LACUMTND
!USE YOMSC2   , ONLY : NSTABUF
!USE YOMCFU   , ONLY : NFRRC


IMPLICIT NONE
LOGICAL LLTRACE, LLDEBUG

LLTRACE=.TRUE.
LLDEBUG=.TRUE.

!     ------------------------
!     *    Initialize Physics
!     ------------------------

!IF(NFLEVG >= 1) THEN
IF (LLTRACE)  WRITE(*,*) " coucou SURAYOLMD : avant SUPHY"
WRITE(*,FMT='('' ---------------- '')')
WRITE(*,FMT='(''     SUPHY: '')')
WRITE(*,FMT='('' ---------------- '')')
  CALL SUPHY(6)    !!!!! A REVOIR (MPL) argument KULOUT=6 "en dur"

!     ------------------------------------------------------------------

END SUBROUTINE SURAYOLMD
