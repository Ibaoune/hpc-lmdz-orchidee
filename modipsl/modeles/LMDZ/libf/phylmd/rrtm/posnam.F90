SUBROUTINE POSNAM(KULNAM,CDNAML)

!**** *POSNAM* - position namelist file for reading

!     Purpose.
!     --------
!     To position namelist file at correct place for reading
!     namelist CDNAML. Replaces use of Cray specific ability
!     to skip to the correct namelist.

!**   Interface.
!     ----------
!        *CALL* *POSNAM*(..)

!        Explicit arguments :     KULNAM - file unit number (input)
!        --------------------     CDNAML - namelist name    (input)

!        Implicit arguments :     None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   None
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 93-06-22
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      R. El Khatib 04-08-10 Apply norms + proper abort if namelist is missing
!     --------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULNAM 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDNAML 


CHARACTER (LEN = 40) ::  CLINE
CHARACTER (LEN =  1) ::  CLTEST

INTEGER(KIND=JPIM) :: ILEN, IND1, ISTATUS, ISCAN
REAL(KIND=JPRB)    :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!      -----------------------------------------------------------

!*       1.    POSITION FILE
!              -------------

IF (LHOOK) CALL DR_HOOK('POSNAM',0,ZHOOK_HANDLE)

CLINE='                                        '
REWIND(KULNAM)
ILEN=LEN(CDNAML)
ISTATUS=0
ISCAN=0
print *,'On cherche a lire:',CDNAML
DO WHILE (ISTATUS==0 .AND. ISCAN==0)
  READ(KULNAM,'(A)',IOSTAT=ISTATUS) CLINE
! print *,'CLINE,ISTATUS= ',CLINE,ISTATUS
  SELECT CASE (ISTATUS)
  CASE (:-1)
    CLINE='POSNAM:CANNOT LOCATE '//CDNAML//' '
    CALL ABOR1(CLINE)
  CASE (0)
    IF (INDEX(CLINE(1:10),'&') == 0) THEN
      ISCAN=0
    ELSE
      IND1=INDEX(CLINE,'&'//CDNAML)
      IF (IND1 == 0) THEN
        ISCAN=0
      ELSE
        CLTEST=CLINE(IND1+ILEN+1:IND1+ILEN+1)
        IF (   (LGE(CLTEST,'0').AND.LLE(CLTEST,'9')) &
         & .OR.(LGE(CLTEST,'A').AND.LLE(CLTEST,'Z')) ) THEN
          ISCAN=0
        ELSE
          ISCAN=1
        ENDIF
      ENDIF
    ENDIF
  CASE (1:)
    CLINE='POSNAM:READ ERROR IN NAMELIST FILE'
    CALL ABOR1(CLINE)
  END SELECT
ENDDO
BACKSPACE(KULNAM)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('POSNAM',1,ZHOOK_HANDLE)
END SUBROUTINE POSNAM
