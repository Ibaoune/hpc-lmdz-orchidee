#ifdef RS6K
@PROCESS HOT NOSTRICT
#endif
SUBROUTINE SRTM_SPCVRT_MCICA &
 & ( KLEV    , KMOL    , KSW    , KCOLS  , PONEMINUS, &
 &   PAVEL   , PTAVEL  , PZ     , PTZ    , PTBOUND  , PALBD   , PALBP, &
 &   PFRCL   , PTAUC   , PASYC  , POMGC  , PTAUA    , PASYA   , POMGA , PRMU0, &
 &   PCOLDRY , PWKL, &
 &   KLAYTROP, KLAYSWTCH, KLAYLOW ,&
 &   PCO2MULT, PCOLCH4  , PCOLCO2 , PCOLH2O , PCOLMOL  , PCOLN2O  , PCOLO2 , PCOLO3 ,&
 &   PFORFAC , PFORFRAC , KINDFOR , PSELFFAC, PSELFFRAC, KINDSELF ,&
 &   PFAC00  , PFAC01   , PFAC10  , PFAC11 ,&
 &   KJP     , KJT      , KJT1 ,&
!-- output arrays 
 &   PBBFD, PBBFU, PBBCD, PBBCU )

! &   PBBFD, PBBFU, PUVFD, PUVFU, PVSFD, PVSFU , PNIFD , PNIFU ,&
! &   PBBCD, PBBCU, PUVCD, PUVCU, PVSCD, PVSCU , PNICD , PNICU &
! & )  

!**** *SRTM_SPCVRT* - SPECTRAL LOOP TO COMPUTE THE SHORTWAVE RADIATION FLUXES.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE TWO-STREAM METHOD OF BARKER

!**   INTERFACE.
!     ----------

!          *SRTM_SPCVRT_MCICA* IS CALLED FROM *SRTM_SRTM_224GP*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          *SWVRTQDR*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION
!     AUTHOR.
!     -------
!        from Howard Barker
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 03-02-27
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        JJMorcrette   20050110 McICA version
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARSRTM  , ONLY : JPLAY, JPB1, JPB2, JPGPT

USE YOESRTWN , ONLY : NGC
USE YOERDI   , ONLY : REPCLC

!USE YOERAD   , ONLY : NSW
!USE YOERDU   , ONLY : RCDAY
!USE YOESWN   , ONLY : NTBANDS, NBANDS, NGS, NUV, NVS, RWGT, NDBUG

IMPLICIT NONE

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM),INTENT(IN)    :: KSW
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOLS
 
INTEGER(KIND=JPIM)               :: KLEV ! UNDETERMINED INTENT
INTEGER(KIND=JPIM)               :: KMOL ! Argument NOT used
!INTEGER(KIND=JPIM)               :: KPT

REAL(KIND=JPRB)                  :: PONEMINUS ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PAVEL(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: PTAVEL(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: PZ(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: PTZ(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: PTBOUND ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRCL(KCOLS,JPLAY)  ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUC(JPLAY,KCOLS)  ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PASYC(JPLAY,KCOLS)  ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMGC(JPLAY,KCOLS)  ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUA(JPLAY,KSW)    ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PASYA(JPLAY,KSW)    ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMGA(JPLAY,KSW)    ! bottom to top
REAL(KIND=JPRB)                  :: PRMU0 ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PCOLDRY(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: PWKL(35,JPLAY) ! Argument NOT used
INTEGER(KIND=JPIM)               :: KLAYTROP ! UNDETERMINED INTENT
INTEGER(KIND=JPIM)               :: KLAYSWTCH ! Argument NOT used
INTEGER(KIND=JPIM)               :: KLAYLOW ! Argument NOT used
REAL(KIND=JPRB)                  :: PCO2MULT(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: PCOLCH4(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PCOLCO2(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PCOLH2O(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PCOLMOL(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PCOLN2O(JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: PCOLO2(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PCOLO3(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PFORFAC(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PFORFRAC(JPLAY) ! UNDETERMINED INTENT
INTEGER(KIND=JPIM)               :: KINDFOR(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PSELFFAC(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PSELFFRAC(JPLAY) ! UNDETERMINED INTENT
INTEGER(KIND=JPIM)               :: KINDSELF(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PFAC00(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PFAC01(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PFAC10(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)                  :: PFAC11(JPLAY) ! UNDETERMINED INTENT
INTEGER(KIND=JPIM)               :: KJP(JPLAY) ! UNDETERMINED INTENT
INTEGER(KIND=JPIM)               :: KJT(JPLAY) ! UNDETERMINED INTENT
INTEGER(KIND=JPIM)               :: KJT1(JPLAY) ! UNDETERMINED INTENT
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBBFD(JPLAY+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBBFU(JPLAY+1) 
!REAL(KIND=JPRB)                  :: PUVFD(JPLAY+1) ! Argument NOT used
!REAL(KIND=JPRB)                  :: PUVFU(JPLAY+1) ! Argument NOT used
!REAL(KIND=JPRB)                  :: PVSFD(JPLAY+1) ! Argument NOT used
!REAL(KIND=JPRB)                  :: PVSFU(JPLAY+1) ! Argument NOT used
!REAL(KIND=JPRB)                  :: PNIFD(JPLAY+1) ! Argument NOT used
!REAL(KIND=JPRB)                  :: PNIFU(JPLAY+1) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBBCD(JPLAY+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBBCU(JPLAY+1) 
!REAL(KIND=JPRB)                  :: PUVCD(JPLAY+1) ! Argument NOT used
!REAL(KIND=JPRB)                  :: PUVCU(JPLAY+1) ! Argument NOT used
!REAL(KIND=JPRB)                  :: PVSCD(JPLAY+1) ! Argument NOT used
!REAL(KIND=JPRB)                  :: PVSCU(JPLAY+1) ! Argument NOT used
!REAL(KIND=JPRB)                  :: PNICD(JPLAY+1) ! Argument NOT used
!REAL(KIND=JPRB)                  :: PNICU(JPLAY+1) ! Argument NOT used
!     ------------------------------------------------------------------

!              ------------

LOGICAL :: LLRTCHK(JPLAY)

REAL(KIND=JPRB) :: &
 & ZCLEAR      , ZCLOUD       &
 & , ZDBT(JPLAY+1) &
 & , ZGCC(JPLAY)   , ZGCO(JPLAY)     &
 & , ZOMCC(JPLAY)  , ZOMCO(JPLAY)    &
 & , ZRDND(JPLAY+1), ZRDNDC(JPLAY+1)&
 & , ZREF(JPLAY+1) , ZREFC(JPLAY+1) , ZREFO(JPLAY+1)  &
 & , ZREFD(JPLAY+1), ZREFDC(JPLAY+1), ZREFDO(JPLAY+1) &
 & , ZRUP(JPLAY+1) , ZRUPD(JPLAY+1) &
 & , ZRUPC(JPLAY+1), ZRUPDC(JPLAY+1)&
 & , ZTAUC(JPLAY)  , ZTAUO(JPLAY)    &
 & , ZTDBT(JPLAY+1) &
 & , ZTRA(JPLAY+1) , ZTRAC(JPLAY+1) , ZTRAO(JPLAY+1)  &
 & , ZTRAD(JPLAY+1), ZTRADC(JPLAY+1), ZTRADO(JPLAY+1)   
REAL(KIND=JPRB) :: &
 & ZDBTC(JPLAY+1), ZTDBTC(JPLAY+1), ZINCFLX(JPGPT)  &
 & ,  ZINCF14(14)   , ZINCTOT   
  
INTEGER(KIND=JPIM) :: IB1, IB2, IBM, IGT, IKL, IW, JB, JG, JK, I_KMODTS

REAL(KIND=JPRB) :: ZARG1, ZARG2, ZDBTMC, ZDBTMO, ZF, ZINCFLUX, ZWF

!-- Output of SRTM_TAUMOLn routines

REAL(KIND=JPRB) :: ZTAUG(JPLAY,16), ZTAUR(JPLAY,16), ZSFLXZEN(16)

!-- Output of SRTM_VRTQDR routine
REAL(KIND=JPRB) :: &
 & ZCD(JPLAY+1,JPGPT), ZCU(JPLAY+1,JPGPT) &
 & ,  ZFD(JPLAY+1,JPGPT), ZFU(JPLAY+1,JPGPT)  
REAL(KIND=JPRB) :: ZHOOK_HANDLE


#include "srtm_taumol16.intfb.h"
#include "srtm_taumol17.intfb.h"
#include "srtm_taumol18.intfb.h"
#include "srtm_taumol19.intfb.h"
#include "srtm_taumol20.intfb.h"
#include "srtm_taumol21.intfb.h"
#include "srtm_taumol22.intfb.h"
#include "srtm_taumol23.intfb.h"
#include "srtm_taumol24.intfb.h"
#include "srtm_taumol25.intfb.h"
#include "srtm_taumol26.intfb.h"
#include "srtm_taumol27.intfb.h"
#include "srtm_taumol28.intfb.h"
#include "srtm_taumol29.intfb.h"
#include "srtm_reftra.intfb.h"
#include "srtm_vrtqdr.intfb.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_SPCVRT_MCICA',0,ZHOOK_HANDLE)

!-- Two-stream model 1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
! KMODTS is set in SWREFTRA
!NDBUG=4

IB1=JPB1
IB2=JPB2
!print *,'IB1, IB2, KSW, KMOL, KLEV: ', IB1,IB2,KSW,KMOL,KLEV

IW=0
ZINCFLUX=0.0_JPRB
ZINCTOT=0.0_JPRB

JB=IB1-1
DO JB = IB1, IB2
  IBM = JB-15
  IGT = NGC(IBM)
  ZINCF14(IBM)=0.0_JPRB

!  print *,'=== spectral band === JB= ',JB,' ====== i.e. IBM= ',IBM,' with IGT= ',IGT
        
!-- for each band, computes the gaseous and Rayleigh optical thickness 
!  for all g-points within the band

  IF (JB == 16) THEN
    CALL SRTM_TAUMOL16 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01   , PFAC10   , PFAC11   ,&
     &   KJP     , KJT      , KJT1     , PONEMINUS,&
     &   PCOLH2O , PCOLCH4  , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC , PSELFFRAC, KINDSELF, PFORFAC  , PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG    , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL16'

  ELSEIF (JB == 17) THEN
    CALL SRTM_TAUMOL17 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLCO2 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL17'

  ELSEIF (JB == 18) THEN
    CALL SRTM_TAUMOL18 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLCH4 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL18'

  ELSEIF (JB == 19) THEN
    CALL SRTM_TAUMOL19 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLCO2 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL19'

  ELSEIF (JB == 20) THEN
    CALL SRTM_TAUMOL20 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLCH4 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL20'

  ELSEIF (JB == 21) THEN
    CALL SRTM_TAUMOL21 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLCO2 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL21'

  ELSEIF (JB == 22) THEN
    CALL SRTM_TAUMOL22 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLMOL , PCOLO2   ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL22'

  ELSEIF (JB == 23) THEN
    CALL SRTM_TAUMOL23 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLMOL ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL23'

  ELSEIF (JB == 24) THEN
    CALL SRTM_TAUMOL24 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLMOL , PCOLO2   , PCOLO3 ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL24'

  ELSEIF (JB == 25) THEN
!--- visible 16000-22650 cm-1   0.4415 - 0.6250 um
    CALL SRTM_TAUMOL25 &
     & ( KLEV     ,&
     &   PFAC00   , PFAC01  , PFAC10 , PFAC11 ,&
     &   KJP      , KJT     , KJT1   , PONEMINUS ,&
     &   PCOLH2O  , PCOLMOL , PCOLO3 ,&
     &   KLAYTROP ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR   &
     & )  
!    print *,'After  SRTM_TAUMOL25'

  ELSEIF (JB == 26) THEN
!--- UV-A 22650-29000 cm-1   0.3448 - 0.4415 um
    CALL SRTM_TAUMOL26 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLCO2 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL26'

  ELSEIF (JB == 27) THEN
!--- UV-B 29000-38000 cm-1   0.2632 - 0.3448 um
    CALL SRTM_TAUMOL27 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLMOL , PCOLO3 ,&
     &   KLAYTROP ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL27'

  ELSEIF (JB == 28) THEN
!--- UV-C 38000-50000 cm-1   0.2000 - 0.2632 um
    CALL SRTM_TAUMOL28 &
     & ( KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10 , PFAC11 ,&
     &   KJP     , KJT     , KJT1   , PONEMINUS ,&
     &   PCOLMOL , PCOLO2  , PCOLO3 ,&
     &   KLAYTROP ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR  &
     & )  
!    print *,'After  SRTM_TAUMOL28'

  ELSEIF (JB == 29) THEN
    CALL SRTM_TAUMOL29 &
     & ( KLEV     ,&
     &   PFAC00   , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP      , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O  , PCOLCO2 , PCOLMOL  ,&
     &   KLAYTROP , PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN , ZTAUG   , ZTAUR    &
     & )  
!    print *,'After  SRTM_TAUMOL29'

  ENDIF

!  IF (NDBUG.LE.3) THEN
!    print *,'Incident Solar Flux'
!    PRINT 9010,(ZSFLXZEN(JG),JG=1,16)
  9010 format(1x,'SolFlx ',16F8.4)
!    print *,'Optical thickness for molecular absorption for JB= ',JB 
!    DO JK=1,KLEV
!      PRINT 9011,JK,(ZTAUG(JK,JG),JG=1,16)
  9011  format(1x,'TauGas ',I3,16E9.2)
!    ENDDO
!    print *,'Optical thickness for Rayleigh scattering for JB= ',JB 
!    DO JK=1,KLEV
!      PRINT 9012,JK,(ZTAUR(JK,JG),JG=1,16)
  9012  format(1x,'TauRay ',I3,16E9.2)
!    ENDDO
!  ENDIF

  DO JG=1,IGT
    IW=IW+1

!    IF (NDBUG.LE.1) THEN
!      print *,' === JG= ',JG,' === for JB= ',JB,' with IW, IBM, JPLAY, KLEV=',IW,IBM,JPLAY,KLEV
!    ENDIF
!    IF (NDBUG.LE.3) THEN
!      print *,'Cloud optical properties for JB= ',JB
!      DO JK=1,KLEV
!        PRINT 9013,JK,PFRCL(IW,JK),PTAUC(JK,IW),POMGC(JK,IW),PASYC(JK,IW)
  9013   format(1x,'Cloud optprop ',I3,f8.4,f8.3,2f8.5)
!      ENDDO
!    ENDIF

    ZINCFLX(IW) =ZSFLXZEN(JG)*PRMU0
    ZINCFLUX    =ZINCFLUX+ZSFLXZEN(JG)*PRMU0           
    ZINCTOT     =ZINCTOT+ZSFLXZEN(JG)
    ZINCF14(IBM)=ZINCF14(IBM)+ZSFLXZEN(JG)

!-- CALL to compute layer reflectances and transmittances for direct 
!  and diffuse sources, first clear then cloudy.
!   Use direct/parallel albedo for direct radiation and diffuse albedo
!   otherwise.

! ZREFC(JK)  direct albedo for clear
! ZREFO(JK)  direct albedo for cloud
! ZREFDC(JK) diffuse albedo for clear
! ZREFDO(JK) diffuse albedo for cloud
! ZTRAC(JK)  direct transmittance for clear
! ZTRAO(JK)  direct transmittance for cloudy
! ZTRADC(JK) diffuse transmittance for clear
! ZTRADO(JK) diffuse transmittance for cloudy

! ZREF(JK)   direct reflectance
! ZREFD(JK)  diffuse reflectance
! ZTRA(JK)   direct transmittance
! ZTRAD(JK)  diffuse transmittance

! ZDBTC(JK)  clear direct beam transmittance
! ZDBTO(JK)  cloudy direct beam transmittance
! ZDBT(JK)   layer mean direct beam transmittance
! ZTDBT(JK)  total direct beam transmittance at levels

!-- clear-sky    
!----- TOA direct beam    
    ZTDBTC(1)=1._JPRB
!----- surface values
    ZDBTC(KLEV+1) =0.0_JPRB
    ZTRAC(KLEV+1) =0.0_JPRB
    ZTRADC(KLEV+1)=0.0_JPRB
    ZREFC(KLEV+1) =PALBP(IBM)
    ZREFDC(KLEV+1)=PALBD(IBM)
    ZRUPC(KLEV+1) =PALBP(IBM)
    ZRUPDC(KLEV+1)=PALBD(IBM)
           
!-- total sky    
!----- TOA direct beam    
    ZTDBT(1)=1._JPRB
!----- surface values
    ZDBT(KLEV+1) =0.0_JPRB
    ZTRA(KLEV+1) =0.0_JPRB
    ZTRAD(KLEV+1)=0.0_JPRB
    ZREF(KLEV+1) =PALBP(IBM)
    ZREFD(KLEV+1)=PALBD(IBM)
    ZRUP(KLEV+1) =PALBP(IBM)
    ZRUPD(KLEV+1)=PALBD(IBM)
!    if (NDBUG < 2) print *,'SWSPCTRL after 1 with JB,JG,IBM and IW= ',JB,JG,IBM,IW
    

!-- NB: a two-stream calculations from top to bottom, but RRTM_SW quantities 
!       are given bottom to top (argh!)
!       Inputs for clouds and aerosols are bottom to top as inputs

    DO JK=1,KLEV
      IKL=KLEV+1-JK

!-- clear-sky optical parameters      
      LLRTCHK(JK)=.TRUE.

!      print 9000,JK,JG,IKL,ZTAUR(IKL,JG),ZTAUG(IKL,JG),PTAUC(IKL,IW)
      9000  format(1x,'Cloud quantities ',3I4,3E12.5)

!-- original
!      ZTAUC(JK)=ZTAUR(IKL,JG)+ZTAUG(IKL,JG)
!      ZOMCC(JK)=ZTAUR(IKL,JG)/ZTAUC(JK)
!      ZGCC (JK)=0.0001_JPRB

!-- total sky optical parameters        
!      ZTAUO(JK)=ZTAUR(IKL,JG)+ZTAUG(IKL,JG)+PTAUC(IKL,IW)
!      ZOMCO(JK)=PTAUC(IKL,IW)*POMGC(IKL,IW)+ZTAUR(IKL,JG)
!      ZGCO (JK)=(PTAUC(IKL,IW)*POMGC(IKL,IW)*PASYC(IKL,IW) &
!        & +ZTAUR(IKL,JG)*0.0001_JPRB)/ZOMCO(JK)
!      ZOMCO(JK)=ZOMCO(JK)/ZTAUO(JK)

!-- clear-sky optical parameters including aerosols
      ZTAUC(JK) = ZTAUR(IKL,JG) + ZTAUG(IKL,JG) + PTAUA(IKL,IBM)
      ZOMCC(JK) = ZTAUR(IKL,JG)*1.0_JPRB + PTAUA(IKL,IBM)*POMGA(IKL,IBM)
      ZGCC (JK) = PASYA(IKL,IBM)*POMGA(IKL,IBM)*PTAUA(IKL,IBM) / ZOMCC(JK)
      ZOMCC(JK) = ZOMCC(JK) / ZTAUC(JK)

    ENDDO
    DO JK=1,KLEV
      IKL=KLEV+1-JK
!-- total sky optical parameters        
      ZTAUO(JK) = ZTAUR(IKL,JG) + ZTAUG(IKL,JG) + PTAUA(IKL,IBM) + PTAUC(IKL,IW)
      ZOMCO(JK) = PTAUA(IKL,IBM)*POMGA(IKL,IBM) + PTAUC(IKL,IW)*POMGC(IKL,IW) &
       & + ZTAUR(IKL,JG)*1.0_JPRB  
      ZGCO (JK) = (PTAUC(IKL,IW)*POMGC(IKL,IW)*PASYC(IKL,IW)  &
       & +  PTAUA(IKL,IBM)*POMGA(IKL,IBM)*PASYA(IKL,IBM)) &
       & /  ZOMCO(JK)  
      ZOMCO(JK) = ZOMCO(JK) / ZTAUO(JK)

!      if (NDBUG <2) THEN
!        print 9001,JK,JG,LRTCHK(JK),0.00,ZTAUC(JK),ZOMCC(JK),ZGCC(JK),ZTAUR(IKL,JG),ZTAUG(IKL,JG)
      9001    format(1x,'clear :',2I3,L4,7(1x,E13.6))
!        print 9002,JK,JG,LRTCHK(JK),PFRCL(IW,IKL),ZTAUO(JK),ZOMCO(JK),ZGCO(JK) &
!          &,PTAUC(IKL,IW),POMGC(IKL,IW),PASYC(IKL,IW)
      9002    format(1x,'total0:',2I3,L4,7(1x,E13.6))
!      end if
    ENDDO    
!    if (NDBUG < 2) print *,'SWSPCTRL after 2'

!-- Delta scaling for clear-sky / aerosol optical quantities
    DO  JK=1,KLEV
      ZF=ZGCC(JK)*ZGCC(JK)
      ZWF=ZOMCC(JK)*ZF
      ZTAUC(JK)=(1._JPRB-ZWF)*ZTAUC(JK)
      ZOMCC(JK)=(ZOMCC(JK)-ZWF)/(1.0_JPRB-ZWF)
      ZGCC (JK)=(ZGCC(JK)-ZF)/(1.0_JPRB-ZF)
    ENDDO
   
    CALL SRTM_REFTRA ( KLEV, I_KMODTS ,&
     &   LLRTCHK, ZGCC  , PRMU0, ZTAUC , ZOMCC ,&
     &   ZREFC  , ZREFDC, ZTRAC, ZTRADC )  
!    if (NDBUG < 2) print *,'SWSPCTR after SWREFTRA for clear-sky'
    
!-- Delta scaling for cloudy quantities
    DO JK=1,KLEV
      IKL=KLEV+1-JK
      LLRTCHK(JK)=.FALSE.
      ZF=ZGCO(JK)*ZGCO(JK)
      ZWF=ZOMCO(JK)*ZF
      ZTAUO(JK)=(1._JPRB-ZWF)*ZTAUO(JK)
      ZOMCO(JK)=(ZOMCO(JK)-ZWF)/(1._JPRB-ZWF)
      ZGCO (JK)=(ZGCO(JK)-ZF)/(1._JPRB-ZF)
      LLRTCHK(JK)=(PFRCL(IW,IKL) > REPCLC)

!      if (NDBUG < 2) THEN 
!        print 9003,JK,LRTCHK(JK),PFRCL(IW,IKL),ZTAUO(JK),ZOMCO(JK),ZGCO(JK) &
!          &,PTAUC(IKL,IW),POMGC(IKL,IW),PASYC(IKL,IW)
      9003    format(1x,'totalD:',I3,L4,7(1x,E13.6))
!      end if

    ENDDO
!    if (NDBUG < 2) print *,'SWSPCTR after Delta scaling'
    
    CALL SRTM_REFTRA ( KLEV, I_KMODTS ,&
     &   LLRTCHK, ZGCO  , PRMU0, ZTAUO , ZOMCO ,&
     &   ZREFO , ZREFDO, ZTRAO, ZTRADO )  
!    if (NDBUG < 2) print *,'SWSPCTR after SWREFTRA for cloudy'

    DO JK=1,KLEV

!-- combine clear and cloudy contributions for total sky

      IKL=KLEV+1-JK 
      ZCLEAR   = 1.0_JPRB - PFRCL(IW,IKL)
      ZCLOUD   = PFRCL(IW,IKL)

      ZREF(JK) = ZCLEAR*ZREFC(JK) + ZCLOUD*ZREFO(JK)
      ZREFD(JK)= ZCLEAR*ZREFDC(JK)+ ZCLOUD*ZREFDO(JK)
      ZTRA(JK) = ZCLEAR*ZTRAC(JK) + ZCLOUD*ZTRAO(JK)
      ZTRAD(JK)= ZCLEAR*ZTRADC(JK)+ ZCLOUD*ZTRADO(JK)

!-- direct beam transmittance        
      ZARG1      = MIN( 200._JPRB, ZTAUC(JK)/PRMU0 )
      ZARG2      = MIN( 200._JPRB, ZTAUO(JK)/PRMU0 )
!      if (PRMU0 <= 0.05_JPRB ) THEN
!        print 9198,JB,IW,JK,PRMU0,ZTAUC(JK),ZTAUO(JK),PTAUC(IKL,IW),ZARG1,ZARG2,ZCLEAR,ZCLOUD,ZTDBT(JK),PFRCL(IW,IKL)
9198    format(1x,'Dbg:',3I4,10E13.6)
!        print 9198,KPT,JB,IW,JK,ZTAUC(JK),ZTAUO(JK),ZARG1,ZARG2,ZCLEAR,ZCLOUD,ZTDBT(JK),PFRCL(IW,IKL)
!9198    format(1x,'Dbg:',4I4,9E13.6)
!      endif
      ZDBTMC     = EXP(-ZARG1 )
      ZDBTMO     = EXP(-ZARG2 )
      ZDBT(JK)   = ZCLEAR*ZDBTMC+ZCLOUD*ZDBTMO
      ZTDBT(JK+1)= ZDBT(JK)*ZTDBT(JK)
        
!-- clear-sky        
      ZDBTC(JK)   =ZDBTMC
      ZTDBTC(JK+1)=ZDBTC(JK)*ZTDBTC(JK)

      
      IF (PRMU0 <= 0.05_JPRB) THEN
!      if (NDBUG < 2) print 9200,JK,ZREFC(JK),ZREFDC(JK),ZTRAC(JK),ZTRADC(JK),ZDBTC(JK),ZTDBTC(JK+1)
!      if (NDBUG < 2) print 9199,JK,ZREF(JK),ZREFD(JK),ZTRA(JK),ZTRAD(JK),ZDBT(JK),ZTDBT(JK+1)
!        print 9200,JK,ZREFC(JK),ZREFDC(JK),ZTRAC(JK),ZTRADC(JK),ZDBTC(JK),ZTDBTC(JK+1),ZCLEAR,ZCLOUD,PRMU0
!        print 9199,JK,ZREF (JK),ZREFD (JK),ZTRA (JK),ZTRAD (JK),ZDBT (JK),ZTDBT (JK+1),ZTAUC(JK),ZTAUO(JK)
      ENDIF
      9199  format(1x,'Comb total:',I3,9E13.6)
      9200  format(1x,'Comb clear:',I3,9E13.6)

    ENDDO           
!    if (NDBUG < 2) print *,'SRTM_SPCVRT after combining clear and cloudy'
                 
!-- vertical quadrature producing clear-sky fluxes

!    print *,'SRTM_SPCVRT after 3 before SRTM_VRTQDR clear'
    
    CALL SRTM_VRTQDR ( KLEV, IW ,&
     &   ZREFC, ZREFDC, ZTRAC , ZTRADC ,&
     &   ZDBTC, ZRDNDC, ZRUPC , ZRUPDC, ZTDBTC ,&
     &   ZCD  , ZCU   )  
      
!    IF (NDBUG < 2) THEN
!      print *,'SRTM_SPCVRT out of SRTM_VRTQDR for clear IW=',IW  
!      DO JK=1,KLEV+1
!        print 9201,JK,ZCD(JK,IW),ZCU(JK,IW)
    9201    format(1x,'clear-sky contrib to fluxes',I3,2F12.4)
!      ENDDO      
!    ENDIF

!-- vertical quadrature producing cloudy fluxes

!    print *,'SRTM_SPCVRT after 4 before SRTM_VRTQDR cloudy'
    
    CALL SRTM_VRTQDR ( KLEV, IW ,&
     &   ZREF , ZREFD , ZTRA , ZTRAD ,&
     &   ZDBT , ZRDND , ZRUP , ZRUPD , ZTDBT ,&
     &   ZFD  , ZFU   )  
 
!    IF (NDBUG < 2) THEN     
!      print *,'SRTM_SPCVRT out of SRTM_VRTQDR for cloudy IW=',IW
!      DO JK=1,KLEV+1
!        print 9202,JK,ZFD(JK,IW),ZFU(JK,IW)
    9202    format(1x,'cloudy sky contrib to fluxes',I3,2F12.4)
!      ENDDO      
!    ENDIF

!-- up and down-welling fluxes at levels
    DO JK=1,KLEV+1
!-- accumulation of spectral fluxes          
      PBBFU(JK) = PBBFU(JK) + ZINCFLX(IW)*ZFU(JK,IW)
      PBBFD(JK) = PBBFD(JK) + ZINCFLX(IW)*ZFD(JK,IW)
      PBBCU(JK) = PBBCU(JK) + ZINCFLX(IW)*ZCU(JK,IW)
      PBBCD(JK) = PBBCD(JK) + ZINCFLX(IW)*ZCD(JK,IW)

! to get NIR, visible and UV quantities

!      PBBFU(JK)=PBBFU(JK)+RWGT(IW)*ZFU(JK,IW)
!      PBBFD(JK)=PBBFD(JK)+RWGT(IW)*ZFD(JK,IW)
!      PBBCU(JK)=PBBCU(JK)+RWGT(IW)*ZCU(JK,IW)
!      PBBCD(JK)=PBBCD(JK)+RWGT(IW)*ZCD(JK,IW)
!      IF (IW <= NUV) THEN
!        PUVFD(JK)=PUVFD(JK)+RWGT(IW)*ZFD(JK,IW)
!        PUVFU(JK)=PUVFU(JK)+RWGT(IW)*ZFU(JK,IW)
!        PUVCD(JK)=PUVCD(JK)+RWGT(IW)*ZCD(JK,IW)
!        PUVCU(JK)=PUVCU(JK)+RWGT(IW)*ZCU(JK,IW)
!      ELSE IF (IW == NUV+1 .AND. IW <= NVS) THEN  
!        PVSFD(JK)=PVSFD(JK)+RWGT(IW)*ZFD(JK,IW)
!        PVSFU(JK)=PVSFU(JK)+RWGT(IW)*ZFU(JK,IW)
!        PVSCD(JK)=PVSCD(JK)+RWGT(IW)*ZCD(JK,IW)
!        PVSCU(JK)=PVSCU(JK)+RWGT(IW)*ZCU(JK,IW)
!      ELSE IF (IW > NVS) THEN  
!        PNIFD(JK)=PNIFD(JK)+RWGT(IW)*ZFD(JK,IW)
!        PNIFU(JK)=PNIFU(JK)+RWGT(IW)*ZFU(JK,IW)
!        PNICD(JK)=PNICD(JK)+RWGT(IW)*ZCD(JK,IW)
!        PNICU(JK)=PNICU(JK)+RWGT(IW)*ZCU(JK,IW)
!      ENDIF  
!      if (NDBUG < 2) then
!!      if (JG.EQ.IGT) THEN 
!           print 9206,JB,JG,JK,IW,PBBCU(JK),PBBCD(JK),PBBFU(JK),PBBFD(JK)
      9206      format(1x,'fluxes up to:',3I3,I4,6E13.6)       
!      end if
    ENDDO

!    if (NDBUG < 2) print *,'SRTM_SPCVRT end of JG=',JG,' for JB=',JB,' i.e. IW=',IW
  ENDDO             
!-- end loop on JG

!  print *,' --- JB= ',JB,' with IB1, IB2= ',IB1,IB2
ENDDO                    
!-- end loop on JB
!if (NDBUG < 2) print *,'SRTM_SPCVRT about to come out'
!print *,'SRTM_SPCVRT about to come out'

!DO IBM=1,14
!  print 9301,IBM,ZINCF14(IBM), ZINCTOT, ZINCF14(IBM)/ZINCTOT
9301 format(1x,'Incident Spectral Flux: ',I3,2E15.8,F12.8)
!ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_SPCVRT_MCICA',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_SPCVRT_MCICA
