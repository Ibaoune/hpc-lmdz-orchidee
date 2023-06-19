!OPTIONS XOPT(NOEVAL)
SUBROUTINE SUCLOP15

!**** *SUCLOP15*  - INITIALIZE COMMON YOMCLOP15
!****               FROZEN VERSION (CYCLE 15) OF SUCLOP

!     PURPOSE.
!     --------
!           INITIALIZE YOMCLOP15, WITH CLOUD OPTICAL PARAMETERS

!**   INTERFACE.
!     ----------
!        *CALL*  SUCLOP15
!     FROM *SUPHEC*

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOMCLOP15

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "INTEGRATED FORECASTING SYSTEM"

!     AUTHOR.
!     -------
!        96-11: Ph. Dandin. Meteo-France
!        ORIGINAL : J.-J. MORCRETTE         *ECMWF*

!     MODIFICATIONS.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P. Marquet    14-Feb-2006 REFFWIA15 + NAMCLOP15 introduced

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
 
! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN   , ONLY : NULNAM, NULOUT
USE YOMLUN   , ONLY : NULOUT
USE YOMCLOP15, ONLY : RYFWCA15 ,RYFWCB15 ,RYFWCC15 ,RYFWCD15 ,&
 & RYFWCE15 ,RYFWCF15 ,REBCUA15 ,REBCUB15 ,REBCUC15 ,&
 & REBCUD15 ,REBCUE15 ,REBCUF15 ,REBCUG15 ,REBCUH15 ,&
 & REFFIA15 ,REFFIB15 ,RTIW15   ,RRIW15   ,REFFWIA15 

IMPLICIT NONE

#include "posnam.intfb.h"
#include "namclop15.h"


!*       1.    Set default values.
!              -------------------

!* Ice cloud properties - crystal: adapted from Ebert and Curry, 1992

! SW : 2 spectral intervals

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SUCLOP15',0,ZHOOK_HANDLE)
REBCUA15(1)= 3.448E-03_JPRB
REBCUA15(2)= 3.448E-03_JPRB
REBCUB15(1)= 2.431_JPRB
REBCUB15(2)= 2.431_JPRB
REBCUC15(1)= 0.99999_JPRB
REBCUC15(2)= 0.975634_JPRB
REBCUD15(1)= 0._JPRB
REBCUD15(2)= 2.487E-04_JPRB
REBCUE15(1)= 0.7661_JPRB
REBCUE15(2)= 0.7866_JPRB
REBCUF15(1)= 5.851E-04_JPRB
REBCUF15(2)= 5.937E-04_JPRB

! LW : spectrally averaged with reference Planck function at 257 K

REBCUG15= 1.07677_JPRB
REBCUH15= 0.00267_JPRB

! Ice particle Effective Radius as a function of LWC

REFFIA15= 40._JPRB
REFFIB15=  0._JPRB

! Water and Ice particle Effective Radius in the RADLSW15 formulae

REFFWIA15= 10._JPRB

!* Water cloud properties - from Fouquart (1987)

! SW : 2 spectral intervals: parameters as a function of Reff

RYFWCA15(1)= 0._JPRB
RYFWCA15(2)= 0._JPRB
RYFWCB15(1)= 1.5_JPRB
RYFWCB15(2)= 1.5_JPRB
RYFWCC15(1)= 0.9999_JPRB
RYFWCC15(2)= 0.9988_JPRB
RYFWCD15(1)= 5.000E-04_JPRB
RYFWCD15(2)= 2.500E-03_JPRB
RYFWCE15(1)= 0.5_JPRB
RYFWCE15(2)= 0.05_JPRB
RYFWCF15(1)= 0.865_JPRB
RYFWCF15(2)= 0.910_JPRB

!* Liquid/Solid water transition

RTIW15= 263._JPRB
RRIW15= 20._JPRB


!*       2.    Modify default values.
!              ----------------------

! Ce qui concerne NAMCLOP15 commente par MPL le 15.04.09
!CALL POSNAM(NULNAM,'NAMCLOP15')
!READ       (NULNAM, NAMCLOP15)


!*       3.    Print final values.
!              -------------------

WRITE(UNIT=NULOUT,FMT='('' - - - - - - - - -'')')
WRITE(UNIT=NULOUT,FMT='('' COMMON YOMCLOP15 '')')
WRITE(UNIT=NULOUT,FMT='('' - - - - - - - - -'')')
WRITE(UNIT=NULOUT,FMT='( '' REFFWIA15 = '',E11.4 )') &
                          & REFFWIA15

IF (LHOOK) CALL DR_HOOK('SUCLOP15',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE SUCLOP15
