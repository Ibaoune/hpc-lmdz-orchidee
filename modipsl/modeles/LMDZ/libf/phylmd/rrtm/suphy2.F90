!OPTIONS XOPT(NOEVAL)
SUBROUTINE SUPHY2(KULOUT)

!**** *SUPHY2*   - Initialize common YOMPHY2 physics controlling
!                  constants

!     Purpose.
!     --------
!           Initialize YOMPHY2, the common that contains the parameters
!           for the control part of the physics of the model.

!**   Interface.
!     ----------
!        *CALL* *SUPHY2(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY2

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation ARPEGE

!     Author.
!     -------
!        J.-F. Geleyn .
!        Original : 90-9-1

!     Modifications.
!     --------------
!        R. EL Khatib : 93-04-02 Set-up defaults controled by LECMWF
!        J.-F. Geleyn : 93-08-19 New cloudiness diagnostics.
!        J.-F. Geleyn : 95-04-10 Anti-fibril. Girard-Delage.
!        P. Marquet   : 97-02-18 Value of VETAF=VAH/VP00+VBH.
!        J.M. Piriou  : 97-04-17 XMULAF default value.
!        E. Bazile    : 98-03-10 Introduce XMUCVPP.
!        W. Owcarz    : 2000-03-27 Set a default value for TSPHY
!        R. EL Khatib : 2000-06-13 RIPBLC
!        R. EL Khatib : 2000-08-21 Turbulent gusts setup
!        J.M. Piriou  : 2002-01-10 set default values to operational ones.
!        Modified by R. EL Khatib : 02-03-29 Control XMULAF<0 ; add LMULAF
!        Modified by D. Banciu    : 02-12-09 Introduction of XDAMP
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN   , ONLY : NULNAM
USE YOMCT0B  , ONLY : LECMWF
! commente par MPL 25.11.08
!USE YOMGEM   , ONLY : VALH     ,VBH
USE YOMDIM   , ONLY : NFLEVG
USE YOMPHY2  , ONLY : NTSHM    ,NTSML    ,XMUCVPP  ,LMULAF   ,&
 & XMULAF   ,XDAMP    ,HCLP     ,HTCLS    ,&
 & RIPBLC   ,&
 & LRAFTUR  ,GZ0RAF   ,FACRAF   ,&
 & HVCLS    ,HTSHM    ,HTSML    ,&
 & TSPHY  

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) :: ZVETAF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "posnam.intfb.h"

#include "namphy2.h"
!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

!        1.1 Set implicit default values

IF (LHOOK) CALL DR_HOOK('SUPHY2',0,ZHOOK_HANDLE)
XMULAF=-1.75_JPRB
XMUCVPP=0._JPRB
XDAMP=0._JPRB
HCLP=1500._JPRB
HTCLS=2._JPRB
HVCLS=10._JPRB
HTSHM=0.450_JPRB
HTSML=0.785_JPRB
TSPHY=1._JPRB
RIPBLC=0.5_JPRB
LRAFTUR=.FALSE.
GZ0RAF=10.0_JPRB
FACRAF=15.0_JPRB
LMULAF=.FALSE.

!        1.2 Modify default values according to LECMWF

IF (LECMWF) THEN
ELSE
  LRAFTUR=.TRUE.
ENDIF

!     Remark : values for TSPHY, NTSHM/ML are calculated and not set up.

!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

! Ce qui concerne NAMPHY2 commente par MPL le 15.04.09
!CALL POSNAM(NULNAM,'NAMPHY2')
!READ(NULNAM,NAMPHY2)
!     ------------------------------------------------------------------

!*       3.    Compute cloud transition indexes.
!              ---------------------------------

NTSHM=0
NTSML=0
! commente par MPL 25.11.08
!DO JLEV=1,NFLEVG
!  ZVETAF=(VALH(JLEV)+VBH(JLEV)+VALH(JLEV-1)+VBH(JLEV-1))*0.5_JPRB
!  IF (ZVETAF <= HTSHM) THEN
!    NTSHM=JLEV
!  ENDIF
!  IF (ZVETAF <= HTSML) THEN
!    NTSML=JLEV
!  ENDIF
!ENDDO

!     ------------------------------------------------------------------

!*       4.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMPHY2 '')')
WRITE(UNIT=KULOUT,FMT='('' XMUCVPP = '',E10.4,'' XMULAF = '',E10.4 &
 & ,'' XDAMP = '',E10.4 &
 & ,'' LMULAF = '',L2,/,'' HTCLS = '',E10.4 &
 & ,'' HVCLS = '',E10.4,'' HCLP = '',E10.4,/&
 & ,'' RIPBLC  = '',F8.4 &
 & ,'' LRAFTUR = '',L2,'' GZ0RAF = '',E10.4,'' FACRAF = '',E10.4 &
 & ,'' HTSHM = '',F8.4,'' NTSHM = '',I3,'' HTSML = '',F8.4 &
 & ,'' NTSML = '',I3 &
 & )')&
 & XMUCVPP,XMULAF,XDAMP,LMULAF,&
 & HTCLS,HVCLS,HCLP,&
 & RIPBLC,&
 & LRAFTUR,GZ0RAF,FACRAF,&
 & HTSHM,NTSHM,HTSML,NTSML  

!*       5.    Control
!              -------

IF (XMULAF > 0.0_JPRB) THEN
  WRITE(KULOUT,*) 'XMULAF SHOULD BE NEGATIVE'
  CALL ABOR1('SUPHY2 : ABOR1 CALLED')
ENDIF

IF ((XDAMP /= 0.0_JPRB).AND.(XMUCVPP /= 0.0_JPRB)) THEN
  WRITE(UNIT=KULOUT,FMT='(A)') 'INCONSISTENCY BETWEEN XDAMP AND XMUCVPP !'
  CALL ABOR1('XDAMP/=0. IMPLIES XMUCVPP=0.!...')
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPHY2',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHY2
