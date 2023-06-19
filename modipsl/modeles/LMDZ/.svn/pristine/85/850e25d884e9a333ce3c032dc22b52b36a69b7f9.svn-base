SUBROUTINE SUMETHOX

!**** *SUCLDP*   - INITIALIZE MODULE YOEMETH CONTROLLING *METHOX*

!     PURPOSE.
!     --------
!           INITIALIZE YOEMETH

!**   INTERFACE.
!     ----------
!        CALL *SUMETHOX* FROM *SUPHEC*
!              --------        ------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        MODULE YOEMETH

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
!        C.JAKOB   *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 98-04-07
!        Modified : 02-01-29  A.Simmons: increase RQLIM from 3.75e-6
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Cugnet     24-Feb-2012 Adapted for LMDZ
!     ------------------------------------------------------------------

USE YOEMETH   , ONLY : RALPHA1 ,RALPHA2  ,RQLIM   ,&
 & RPBOTOX,  RPBOTPH ,RPTOPOX  ,RPTOPPH ,&
 & RALPHA3,  RLOGPPH  

!*       1.    SET VALUES
!              ----------

IMPLICIT NONE
RALPHA1=(19.*LOG(10.))/(LOG(20.)**4)
RALPHA2=LOG(1.0/3.+0.01)
RQLIM=4.25E-6
RPBOTOX=10000.
RPBOTPH=20.
RPTOPOX=50.
RPTOPPH=0.1
RALPHA3=0.5*(LOG(100.)+RALPHA2)
RLOGPPH=LOG(RPTOPPH/RPBOTPH)

!     -----------------------------------------------------------------

END SUBROUTINE SUMETHOX


