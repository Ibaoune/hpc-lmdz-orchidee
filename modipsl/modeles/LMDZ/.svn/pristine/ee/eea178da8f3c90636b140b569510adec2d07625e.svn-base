SUBROUTINE SUCAPE(KULOUT)

!**** *SUCAPE  * - ROUTINE TO INITIALIZE THE VARIABLES FOR CAPE 
!                  AND CIN COMPUTATION

!     PURPOSE.
!     --------
!        SET DEFAULT VALUES, THEN READS NAMELIST NAMCAPE

!**   INTERFACE.
!     ----------
!        *CALL* *SUCAPE(KULOUT)*

!         EXPLICIT ARGUMENTS :  KULOUT
!         --------------------

!         IMPLICIT ARGUMENTS :
!         --------------------
!            COMMON  YOMCAPE
!            COMMON  YOMPHY
!            COMMON  YOMLUN

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!        N.Pristov 03/2001

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 03/2001
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCAPE  , ONLY :  NCAPEITER   ,NETAPES  ,GCAPERET   ,GCAPEPSD
USE YOMPHY   , ONLY :  NBITER 
! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN   , ONLY :  NULNAM

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "posnam.intfb.h"

#include "namcape.h"

!      ----------------------------------------------------------------
!*       1.    SET DEFAULT VALUES.
!              -------------------
IF (LHOOK) CALL DR_HOOK('SUCAPE',0,ZHOOK_HANDLE)
NCAPEITER=NBITER
NETAPES=2
GCAPERET=0._JPRB
GCAPEPSD=30000._JPRB

!      ----------------------------------------------------------------
!*       2.    MODIFIES DEFAULT VALUES. READ NAMELIST.
!              ------------------------

! Ce qui concerne NAMCAPE commente par MPL le 15.04.09
!CALL POSNAM(NULNAM,'NAMCAPE')
!READ(NULNAM,NAMCAPE)

IF (NCAPEITER <= 0 ) CALL ABOR1('SUCAPE:  INVALID VALUE FOR NCAPEITER')
IF (NETAPES <= 0 ) CALL ABOR1('SUCAPE:  INVALID VALUE FOR NETAPES')
IF ((GCAPERET < 0.0_JPRB ) .OR. (GCAPERET > 1.0_JPRB ))  &
 & CALL ABOR1('SUCAPE:  INVALID VALUE FOR GCAPERET')  
IF (GCAPEPSD <= 0.0_JPRB ) CALL ABOR1('SUCAPE:  INVALID VALUE FOR GCAPEPSD')

!      -----------------------------------------------------------
!*       3.    PRINT FINAL VALUES.
!              -------------------

WRITE(KULOUT,'('' NCAPEITER = '',I2,'' NETAPES = '',I2,'' GCAPERET = '',&
 & E13.6,'' GCAPEPSD = '',E13.6)') NCAPEITER, NETAPES, GCAPERET, GCAPEPSD  
 
IF (LHOOK) CALL DR_HOOK('SUCAPE',1,ZHOOK_HANDLE)
END SUBROUTINE SUCAPE
