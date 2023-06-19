SUBROUTINE SUCLD ( KLEV , PETA )

!**** *SUCLD*   - INITIALIZE COMMON YOECLD CONTROLLING *CLOUD*

!     PURPOSE.
!     --------
!           INITIALIZE YOECLD

!**   INTERFACE.
!     ----------
!        CALL *SUCLD* FROM *SUPHEC*
!              -----        ------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOECLD

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
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-12-15
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOECLD   , ONLY : CETA     ,RANVA    ,RANVB    ,RANVH    ,&
 & RCCA     ,RCCB     ,RCCC     ,RCFCT    ,RCLWMR   ,&
 & RCSCAL   ,RETAHB   ,RETAMB   ,RLOIA    ,RLOIB    ,&
 & RLOIC    ,RLOID    ,RLONIA   ,RLONIB   ,RRHH     ,&
 & RRHL     ,RRHM     ,RGAMMAS  ,REPSCR   ,REPSEC   ,&
 & LOMEGA  

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PETA(KLEV) 
!      ----------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

!*       1.    SET VALUES
!              ----------

IF (LHOOK) CALL DR_HOOK('SUCLD',0,ZHOOK_HANDLE)
RANVA  = 2._JPRB
RANVB  = 0.3_JPRB
RANVH  = 0.4_JPRB
RCCA   = 0.125_JPRB
RCCB   = 1.5_JPRB
RCCC   = 0.8_JPRB
RCFCT  = 0.400_JPRB
RCSCAL = 1.0E+11_JPRB

RETAHB = 0.45_JPRB
RETAMB = 0.80_JPRB

RLOIA  = 1.0E+02_JPRB
RLOIB  =-10.00_JPRB
RLOIC  =-0.9_JPRB
RLOID  = 5.0_JPRB

RLONIA = -0.1_JPRB
RLONIB = -10.0_JPRB

RRHH   = 0.9_JPRB
RRHM   = 0.8_JPRB
RRHL   = 0.70_JPRB

RGAMMAS= 0.05_JPRB
RCLWMR = 1.E-04_JPRB
LOMEGA =.TRUE.

REPSEC = 1.0E-12_JPRB
REPSCR = 1.0E-12_JPRB

DO JK=1,KLEV
  CETA(JK)=PETA(JK)
ENDDO

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUCLD',1,ZHOOK_HANDLE)
END SUBROUTINE SUCLD
