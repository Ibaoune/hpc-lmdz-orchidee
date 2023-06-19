SUBROUTINE METHOX(KIDIA,  KFDIA,  KLON,  KLEV,PQ,     PTENQ,  PAP )

!**** *METHOX*   - Calculate humidity tendencies from methane
!                  oxidation and photolysis

!**   INTERFACE.
!     ----------
!        CALL *METHOX* FROM *CALLPAR*
!              ------        -------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS

!     INPUT PARAMETERS (REAL):

!    *PAP*          PRESSURE                                      PA
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG

!     UPDATED PARAMETERS (REAL):

!    *PTENQ*        TENDENCY OF SPECIFIC HUMIDITY                 KG/(KG*S)

!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        MODULE YOEMETH
!        MODULE YOMCST

!     METHOD.
!     -------
!        SEE RD-MEMO R60.1/AJS/31

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        SEE RD-MEMO R60.1/AJS/31

!     AUTHOR.
!     -------
!        C.JAKOB   *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 98-04-07
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Cugnet     24-Feb-2012 Adapted for LMDZ
!     ------------------------------------------------------------------

USE YOEMETH   , ONLY : RALPHA1 ,RALPHA2  ,RQLIM   ,&
 & RPBOTOX,  RPBOTPH ,RPTOPOX  ,RPTOPPH ,&
 & RALPHA3,  RLOGPPH  

IMPLICIT NONE

#include "YOMCST.h"

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
REAL   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL   ,INTENT(INOUT) :: PTENQ(KLON,KLEV) 
REAL   ,INTENT(IN)    :: PAP(KLON,KLEV) 
LOGICAL :: LLOXID,         LLPHOTO

INTEGER :: JK, JL

REAL :: ZARG, ZPRATIO, ZTAU1, ZTAU2, ZTDAYS

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

    LLOXID=PAP(JL,JK) < RPBOTOX.AND.PQ(JL,JK) < RQLIM
    LLPHOTO=PAP(JL,JK) < RPBOTPH

!     METHANE OXIDATION

    IF(LLOXID) THEN
      IF(PAP(JL,JK) <= RPTOPOX) THEN
        ZTDAYS=100.
      ELSE
        ZPRATIO=(LOG(PAP(JL,JK)/RPTOPOX))**4./LOG(RPBOTOX/PAP(JL,JK))
        ZTDAYS=100.*(1+RALPHA1*ZPRATIO)
      ENDIF
      ZTAU1=86400.*ZTDAYS
      PTENQ(JL,JK)=PTENQ(JL,JK)+(RQLIM-PQ(JL,JK))/ZTAU1
    ENDIF

!     PHOTOLYSIS

    IF(LLPHOTO) THEN
      IF(PAP(JL,JK) <= RPTOPPH) THEN
        ZTDAYS=3.
      ELSE
        ZARG=RALPHA2-RALPHA3*(1+COS((RPI*LOG(PAP(JL,JK)/RPBOTPH))/RLOGPPH))
        ZTDAYS=1.0/(EXP(ZARG)-0.01)
      ENDIF
      ZTAU2=86400.*ZTDAYS
      PTENQ(JL,JK)=PTENQ(JL,JK)-PQ(JL,JK)/ZTAU2
    ENDIF
  ENDDO
ENDDO

END SUBROUTINE METHOX


