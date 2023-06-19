SUBROUTINE GPPREH(KPROMA,KSTART,KPROF,KFLEV,PVAH,PVBH,PRESH)

!**** *GPPREH* - Computes half and full level pressure

!     Purpose.
!     --------
!           Computes pressures at half and full model levels.

!**   Interface.
!     ----------
!        *CALL* *GPPREH(KPROMA,KSTART,KPROF,KFLEV,PVAH,PVBH,PRESH)

!        Explicit arguments :
!        --------------------
!                              KPROMA :  dimensioning
!                              KSTART :  start of work
!                              KPROF  :  depth of work
!                              KFLEV     : vert. dimensioning
!                              PRESH(KPROMA,0:KFLEV) - HALF LEVEL PRESSURE
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

!     Author.
!     -------
!        Erik Andersson, Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 92-11-23
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVAH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVBH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRESH(KPROMA,0:KFLEV) 
INTEGER(KIND=JPIM) :: JLEV, JLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    COMPUTES HALF LEVEL PRESSURES.
!              ------------------------------

IF (LHOOK) CALL DR_HOOK('GPPREH',0,ZHOOK_HANDLE)
DO JLEV=0,KFLEV-1
  DO JLON=KSTART,KPROF
    PRESH(JLON,JLEV)=PVAH(JLEV)+PVBH(JLEV)*PRESH(JLON,KFLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPREH',1,ZHOOK_HANDLE)
END SUBROUTINE GPPREH
