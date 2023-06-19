!OPTIONS XOPT(NOEVAL)
SUBROUTINE SUTOPH(KULOUT)

!**** *SUTOPH*   - Initialize common YOMTOPH top parameterization

!     Purpose.
!     --------
!           Initialize YOMTOPH, the common that contains the top pressure
!           and the first level of parameterization
!           it also contains mesospheric drag vertical profil

!**   Interface.
!     ----------
!        *CALL* *SUTOPH(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMTOPH, YOMSTA

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
!        A. Lasserre-Bigorry

!     Modifications.
!     --------------
!        Original : 91-06-10
!        Modified 92-02-22 by M. Deque (test of consistency between phys. para.)
!        Modified by R. EL Khatib : 93-04-02 Set-up defaults controled by LECMWF
!        Modified 93-11-17 by Ph. Dandin : FMR scheme with MF physics
!        Modified 97-05-17 by M. Deque   : frozen FMR                 
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F.Bouyssel 04-11-22 : NTCOET,ETCOET
!        P. Marquet 05-09-12 : NTAJUC
!        M. Deque   05-09-12 : default RCLX
!        M. Deque   05-09-12 : default TPSCLIM
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMDIM   , ONLY : NFLEVG
! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN   , ONLY : NULNAM
USE YOMCT0B  , ONLY : LECMWF
USE YOMSTA   , ONLY : STPRE
USE YOMTOPH  , ONLY : RMESOU   ,RMESOT   ,NTQSAT   ,NTDIFU   ,&
 & NTCOEF   ,NTDRAG   ,NTCVIM   ,NTPLUI   ,NTRADI   ,&
 & NTNEBU   ,NTOZON   ,NTDRME   ,ETQSAT   ,ETDIFU   ,&
 & ETCOEF   ,ETDRAG   ,ETCVIM   ,ETPLUI   ,ETRADI   ,&
 & ETNEBU   ,ETOZON   ,ETDRME   ,XDRMUK   ,XDRMUX   ,XDRMUP   ,&
 & XDRMTK   ,XDRMTX   ,XDRMTP   ,NTCOET   ,ETCOET   ,&
 & RMESOQ   ,XDRMQK   ,XDRMQP   ,RFMESOQ  ,RCLX     ,&
 & NTAJUC   ,ETAJUC   ,TPSCLIM
USE YOMPHY   , ONLY : LRAY     ,LRAYFM   ,LRAYFM15 ,LRRMES
USE YOEPHY   , ONLY : LAGPHY

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) :: PAP, PAPX, ZMEST, ZMESU, ZMESQ

REAL(KIND=JPRB) :: PMESQF
REAL(KIND=JPRB) :: PMESTF
REAL(KIND=JPRB) :: PMESUF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "namtoph.h"

!     ------------------------------------------------------------------

!*    Mesospheric drag shape function

!     PMESUF(PAP,PAPX) = MAX( (PAPX-PAP)/MAX(PAP**1.5,1.E-10),0. )
PMESUF(PAP,PAPX) = MAX( (PAPX-PAP)/MAX(PAP,1.E-10_JPRB),0.0_JPRB )
PMESTF(PAP,PAPX) = MAX( (PAPX-PAP)/MAX(PAP,1.E-10_JPRB),0.0_JPRB )
PMESQF(PAP,PAPX) = MAX( (PAPX-PAP)/MAX(PAP,1.E-10_JPRB),0.0_JPRB )

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "seapre.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUTOPH',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

!        1.1 Set implicit default values

ETQSAT=0._JPRB
ETDIFU=0._JPRB
ETCOEF=0._JPRB
ETDRAG=0._JPRB
ETCVIM=0._JPRB
ETPLUI=0._JPRB
ETRADI=0._JPRB
ETNEBU=0._JPRB
ETOZON=0._JPRB
ETDRME=0._JPRB
ETCOET=0._JPRB
ETAJUC=0._JPRB
NTQSAT=1
NTDIFU=1
NTCOEF=1
NTDRAG=1
NTCVIM=1
NTPLUI=1
NTRADI=1
NTNEBU=1
NTOZON=1
NTDRME=1
NTCOET=1
NTAJUC=1

XDRMUK=0._JPRB
XDRMUX=0._JPRB
XDRMUP=0._JPRB
XDRMTK=0._JPRB
XDRMTX=0._JPRB
XDRMTP=0._JPRB
XDRMQK=0._JPRB
XDRMQP=0._JPRB

RFMESOQ=3.725E-06_JPRB
RCLX=0.0_JPRB
TPSCLIM=197._JPRB

!        1.2 Modify default values according to LECMWF

IF (LECMWF) THEN
ELSE
ENDIF

!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

! Ce qui concerne NAMTOPH commente par MPL le 15.04.09
!CALL POSNAM(NULNAM,'NAMTOPH')
!READ(NULNAM,NAMTOPH)

!*       2.1  Search corresponding level, to pressure in NAMTOPH
!             for each parameterization

IF(ETQSAT /= 0.0_JPRB) CALL SEAPRE (ETQSAT,NTQSAT,STPRE,NFLEVG)
IF(ETDIFU /= 0.0_JPRB) CALL SEAPRE (ETDIFU,NTDIFU,STPRE,NFLEVG)
IF(ETCOEF /= 0.0_JPRB) CALL SEAPRE (ETCOEF,NTCOEF,STPRE,NFLEVG)
IF(ETDRAG /= 0.0_JPRB) CALL SEAPRE (ETDRAG,NTDRAG,STPRE,NFLEVG)
IF(ETCVIM /= 0.0_JPRB) CALL SEAPRE (ETCVIM,NTCVIM,STPRE,NFLEVG)
IF(ETPLUI /= 0.0_JPRB) CALL SEAPRE (ETPLUI,NTPLUI,STPRE,NFLEVG)
IF(ETRADI /= 0.0_JPRB) THEN
  IF (LRAY) THEN
    CALL SEAPRE (ETRADI,NTRADI,STPRE,NFLEVG)
  ENDIF
  IF (LRAYFM.OR.LRAYFM15) THEN
    ETRADI=0._JPRB
    NTRADI=1
  ENDIF
ENDIF
IF(ETNEBU /= 0.0_JPRB) CALL SEAPRE (ETNEBU,NTNEBU,STPRE,NFLEVG)
IF(ETOZON /= 0.0_JPRB) CALL SEAPRE (ETOZON,NTOZON,STPRE,NFLEVG)
IF(ETDRME /= 0.0_JPRB) CALL SEAPRE (ETDRME,NTDRME,STPRE,NFLEVG)
IF(ETCOET /= 0.0_JPRB) CALL SEAPRE (ETCOET,NTCOET,STPRE,NFLEVG)
IF(ETAJUC /= 0.0_JPRB) CALL SEAPRE (ETAJUC,NTAJUC,STPRE,NFLEVG)
!     ------------------------------------------------------------------

!*       3.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMTOPH '')')
WRITE(UNIT=KULOUT,FMT='('' ETQSAT = '',E10.4,'' NTQSAT = '',I10 &
 & ,'' ETDIFU = '',E10.4,'' NTDIFU = '',I10 &
 & ,/,'' ETCOEF = '',E10.4,'' NTCOEF = '',I10 &
 & ,'' ETDRAG = '',E10.4,'' NTDRAG = '',I10 &
 & ,/,'' ETCVIM = '',E10.4,'' NTCVIM = '',I10 &
 & ,'' ETPLUI = '',E10.4,'' NTPLUI = '',I10 &
 & ,/,'' ETRADI = '',E10.4,'' NTRADI = '',I10 &
 & ,'' ETNEBU = '',E10.4,'' NTNEBU = '',I10 &
 & ,/,'' ETOZON = '',E10.4,'' NTOZON = '',I10 &
 & ,'' ETDRME = '',E10.4,'' NTDRME = '',I10 &
 & ,/,'' ETCOET = '',E10.4,'' NTCOET = '',I10 &
 & ,/,'' ETAJUC = '',E10.4,'' NTAJUC = '',I10 &
 & ,/,'' XDRMUK = '',E10.4,'' XDRMUP = '',E10.4 &
 & ,'' XDRMUX = '',E10.4,'' XDRMTK = '',E10.4 &
 & ,'' XDRMTP = '',E10.4,'' XDRMTX = '',E10.4 &
 & ,'' XDRMQK = '',E11.4,'' XDRMQP = '',E11.4 &
 & ,/,'' RFMESOQ= '',E11.4,'' RCLX   = '',E11.4 &
 & )')&
 & ETQSAT,NTQSAT,ETDIFU,NTDIFU &
 & ,ETCOEF,NTCOEF,ETDRAG,NTDRAG &
 & ,ETCVIM,NTCVIM,ETPLUI,NTPLUI &
 & ,ETRADI,NTRADI,ETNEBU,NTNEBU &
 & ,ETOZON,NTOZON,ETDRME,NTDRME &
 & ,ETCOET,NTCOET &
 & ,ETAJUC,NTAJUC &
 & ,XDRMUK,XDRMUP,XDRMUX,XDRMTK,XDRMTP,XDRMTX &
 & ,XDRMQK,XDRMQP,RFMESOQ,RCLX

!     VERIFICATION OF CONSISTENCY BETWEEN PHYSICAL PARAMETERIZATION

IF (ETCOEF > ETDIFU.OR.ETCOEF > ETDRAG)THEN
  WRITE(UNIT=KULOUT,FMT='('' ETCOEF TOO LOW '')')
  CALL ABOR1('SUTOPH')
ENDIF
IF (ETQSAT > ETNEBU.OR.ETQSAT > ETPLUI.OR.ETQSAT > ETCVIM)THEN
  WRITE(UNIT=KULOUT,FMT='('' ETQSAT TOO LOW '')')
  CALL ABOR1('SUTOPH')
ENDIF
IF (ETCVIM > ETNEBU)THEN
  WRITE(UNIT=KULOUT,FMT='('' ETCVIM TOO LOW '')')
  CALL ABOR1('SUTOPH')
ENDIF

!     ------------------------------------------------------------------

!*       4.    INITIALIZE MESOSPHERIC DRAG FOR U,V AND T
!              -----------------------------------------

IF (LRRMES.AND..NOT.LAGPHY) THEN
  WRITE (UNIT=KULOUT,FMT='('' PROFIL VERTICAL DE DRAG MESO'',/&
   & ,'' LEV'',T15,''VITESSE'',T45,''TEMPERATURE'' &
   & , T65, ''HUMIDITE'' )')  
  DO JLEV=1,NFLEVG
    RMESOU(JLEV)=XDRMUK*PMESUF(STPRE(JLEV),XDRMUP)
    RMESOT(JLEV)=XDRMTK*PMESTF(STPRE(JLEV),XDRMTP)
    RMESOQ(JLEV)=XDRMQK*PMESQF(STPRE(JLEV),XDRMQP)
    IF (XDRMUX /= 0.0_JPRB) RMESOU(JLEV)=MIN(RMESOU(JLEV),XDRMUX)
    IF (XDRMTX /= 0.0_JPRB) RMESOT(JLEV)=MIN(RMESOT(JLEV),XDRMTX)
    ZMESU=1.0_JPRB/MAX(1.E-8_JPRB,RMESOU(JLEV)*3600._JPRB*24._JPRB)
    ZMEST=1.0_JPRB/MAX(1.E-8_JPRB,RMESOT(JLEV)*3600._JPRB*24._JPRB)
    ZMESQ=1.0_JPRB/MAX(1.E-8_JPRB,RMESOQ(JLEV)*3600._JPRB*24._JPRB)
    WRITE (UNIT=KULOUT,FMT='(I3,T10,E9.3,T20,G9.3,T40,E9.3,T50 &
     & ,G9.3, T70,E9.3, T80,G9.3)') JLEV,RMESOU(JLEV),ZMESU, &
     & RMESOT(JLEV),ZMEST, &
     & RMESOQ(JLEV),ZMESQ 
  ENDDO
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUTOPH',1,ZHOOK_HANDLE)
END SUBROUTINE SUTOPH
