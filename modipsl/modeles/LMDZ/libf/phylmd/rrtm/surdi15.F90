!OPTIONS XOPT(NOEVAL)
SUBROUTINE SURDI15

!**** *SURDI15*   - INITIALIZE COMMON YOMRDI15 CONTROLLING RADINT
!****               FROZEN VERSION (CYCLE 15) OF SURDI

!     PURPOSE.
!     --------
!           INITIALIZE YOMRDI15, THE COMMON THAT CONTROLS THE
!           RADIATION INTERFACE

!**   INTERFACE.
!     ----------
!        CALL *SURDI15* FROM *SUECRAD*
!              -------        -------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOMRDI15

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "IN CORE MODEL"

!     AUTHOR.
!     -------
!        96-11: Ph. Dandin. Meteo-France
!        ORIGINAL : 88-12-15 BY JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A.Alias       05-12-2005 greenhouse gases variables (M.Deque)
!        A.Alias       13-06-2006 RI0 value can be changed via namscen.h

!     ------------------------------------------------------------------

USE PARKIND1        ,ONLY : JPIM     ,JPRB
USE YOMHOOK         ,ONLY : LHOOK,   DR_HOOK
! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN          , ONLY : NULNAM
USE YOMLUN_IFSAUX   , ONLY : NULOUT

USE YOMCST          , ONLY : RI0
USE YOMRDI15        , ONLY : RSDTSN15 ,RRAE15   ,RMU0015  ,RALBICE15,&
 & RALBSEA15,RALBSNM15,RALBSNO15,RCARDI15 ,REMISS15 ,&
 & RSNOWAL15,RVLBDC15 ,RCH415   ,RN2O15   ,RCFC1115 ,&
 & RCFC1215 ,REPALB15 ,REPCLC15 ,REPH2O15  
USE YOMRDU15        , ONLY : REPSEC15

IMPLICIT NONE

REAL(KIND=JPRB) :: XCARDI, XCFC11, XCFC12, XCH4, XN2O
REAL(KIND=JPRB) :: ZAIRMWG, ZC11MWG, ZC12MWG, ZCH4MWG, ZCO2MWG, ZN2OMWG, ZSUPSAT
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "namscen.h"
#include "posnam.intfb.h"
!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

IF (LHOOK) CALL DR_HOOK('SURDI15',0,ZHOOK_HANDLE)
RSDTSN15  = 5.0_JPRB

RRAE15 = 0.1277E-02_JPRB
RMU0015 = RRAE15/SQRT(RRAE15*(RRAE15+2.0_JPRB))

RALBICE15 = 0.55_JPRB
RALBSEA15 = 0.07_JPRB
RALBSNO15 = 0.80_JPRB
RALBSNM15 = 0.40_JPRB
RSNOWAL15 = 0.01_JPRB
!*  Concentration of the various trace gases (IPCC/SACC values for 1990)
!        CO2         CH4        N2O        CFC11       CFC12
!      353ppmv     1.72ppmv   310ppbv     280pptv     484pptv

XCARDI  = 353.E-06_JPRB
XCH4    = 1.72E-06_JPRB
XN2O    = 310.E-09_JPRB
XCFC11  = 280.E-12_JPRB
XCFC12  = 484.E-12_JPRB

ZAIRMWG = 28.970_JPRB
ZCO2MWG = 44.011_JPRB
ZCH4MWG = 16.043_JPRB
ZN2OMWG = 44.013_JPRB
ZC11MWG = 137.3686_JPRB
ZC12MWG = 120.9140_JPRB

! Ce qui concerne NAMSCEN commente par MPL le 15.04.09
!CALL POSNAM(NULNAM,'NAMSCEN')
!READ       (NULNAM, NAMSCEN)

WRITE(NULOUT,'( &
 & '' CO2   = '',E14.7,'' CH4   = '',E14.7 &
 & ,'' N2O   = '',E14.7,'' CFC11 = '',E14.7 &
 & ,'' CFC12 = '',E14.7,'' RI0 = '',E14.7 &
 & )') XCARDI,XCH4,XN2O,XCFC11,XCFC12,RI0  

RCARDI15  = XCARDI*ZCO2MWG/ZAIRMWG
RCH415    = XCH4*ZCH4MWG/ZAIRMWG
RN2O15    = XN2O*ZN2OMWG/ZAIRMWG
RCFC1115  = XCFC11*ZC11MWG/ZAIRMWG
RCFC1215  = XCFC12*ZC12MWG/ZAIRMWG
REMISS15  = 0.996_JPRB
!ZSUPSAT = 0.01_JPRB
RVLBDC15  = 0.5_JPRB

REPSEC15=1.E-12_JPRB
REPCLC15=1.E-12_JPRB
REPH2O15=1.E-12_JPRB
REPALB15=1.E-12_JPRB

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURDI15',1,ZHOOK_HANDLE)
END SUBROUTINE SURDI15
