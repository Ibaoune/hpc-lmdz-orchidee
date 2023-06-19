SUBROUTINE SW &
 & ( KIDIA, KFDIA , KLON  , KLEV , KAER,&
 & PSCT , PCARDI, PPSOL , PALBD, PALBP , PWV, PQS,&
 & PRMU0, PCG   , PCLDSW, PDP  , POMEGA, POZ, PPMB,&
 & PTAU , PTAVE , PAER,&
 & PFDOWN, PFUP,&
 & PCDOWN, PCUP,&
 & PFDNN, PFDNV , PFUPN, PFUPV,&
 & PCDNN, PCDNV , PCUPN, PCUPV,&
 & PSUDU, PUVDF , PPARF, PPARCF, PDIFFS , PDIRFS, &
 & LRDUST, PPIZA_DST,PCGA_DST,PTAUREL_DST &
 & )


!**** *SW* - COMPUTES THE SHORTWAVE RADIATION FLUXES.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SW* IS CALLED FROM *RADLSW*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES ABSORBER AMOUNTS                 (SWU)
!          2. COMPUTES FLUXES IN U.V./VISIBLE  SPECTRAL INTERVAL (SW1S)
!          3. COMPUTES FLUXES IN NEAR-INFRARED SPECTRAL INTERVAL (SWNI)

!     EXTERNALS.
!     ----------

!          *SWU*, *SW1S*, *SWNI*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        95-01-01   J.-J. MORCRETTE  Direct/Diffuse Albedo
!        95-12-07   J.-J. MORCRETTE  Near-Infrared in nsw-1 Intervals
!        990128     JJMorcrette      sunshine duration
!        99-05-25   JJMorcrette      Revised aerosols
!        00-12-18   JJMorcrette      6 spectral intervals
!        02-09-01   JJMorcrette      UV and PAR
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Y.Seity  04-11-18 : add two arguments for AROME extern. surface
!        Y.Seity  05-10-10 : add add 3 optional arg. for dust SW properties
!        JJMorcrette 20060721 PP of clear-sky PAR
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!USE YOERAD   , ONLY : NSW
! NSW mis dans .def MPL 20140211
USE write_field_phy

IMPLICIT NONE

include "clesphys.h"

integer, save :: icount=0
!$OMP THREADPRIVATE(icount)
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAER 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCARDI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPSOL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KLON,NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCG(KLON,NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDSW(KLON,KLEV) 
REAL(KIND=JPRB)                  :: PDP(KLON,KLEV) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMEGA(KLON,NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPMB(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU(KLON,NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAVE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) 
!++MODIFCODE
LOGICAL           ,INTENT(IN)    :: LRDUST
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZA_DST(KLON,KLEV,NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGA_DST(KLON,KLEV,NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUREL_DST(KLON,KLEV,NSW)
!--MODIFCODE
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFDOWN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFUP(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCDOWN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCUP(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFDNN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFDNV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFUPN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFUPV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCDNN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCDNV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCUPN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCUPV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSUDU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUVDF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPARF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPARCF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFFS(KLON,NSW) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIRFS(KLON,NSW) 
!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!     ------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZAKI(KLON,2,NSW)&
 & ,  ZCLD(KLON,KLEV)    , ZCLEAR(KLON) &
 & ,  ZDSIG(KLON,KLEV)   , ZFACT(KLON)&
 & ,  ZFD(KLON,KLEV+1)   , ZCD(KLON,KLEV+1)&
 & ,  ZCDOWN(KLON,KLEV+1), ZCDNIR(KLON,KLEV+1), ZCDUVS(KLON,KLEV+1)&
 & ,  ZFDOWN(KLON,KLEV+1), ZFDNIR(KLON,KLEV+1), ZFDUVS(KLON,KLEV+1)&
 & ,  ZFU(KLON,KLEV+1)   , ZCU(KLON,KLEV+1)&
 & ,  ZCUP(KLON,KLEV+1)  , ZCUNIR(KLON,KLEV+1), ZCUUVS(KLON,KLEV+1)&
 & ,  ZFUP(KLON,KLEV+1)  , ZFUNIR(KLON,KLEV+1), ZFUUVS(KLON,KLEV+1)&
 & ,  ZRMU(KLON)         , ZSEC(KLON)         &
 & ,  ZSUDU1(KLON)       , ZSUDU2(KLON)       &
 & ,  ZSUDU1T(KLON)      , ZSUDU2T(KLON)      &
 & ,  ZUD(KLON,5,KLEV+1) ,ZDIFF(KLON,KLEV)   ,ZDIRF(KLON,KLEV)    &
 & ,  ZDIFF2(KLON,KLEV)  , ZDIRF2(KLON,KLEV)

INTEGER(KIND=JPIM) ::  JK, JL, JNU, INUVS, INUIR

REAL(KIND=JPRB) :: ZHOOK_HANDLE
LOGICAL         :: LLDEBUG
character*1 str1

#include "sw1s.intfb.h"
#include "swni.intfb.h"
#include "swu.intfb.h"

!     ------------------------------------------------------------------

!*         1.     ABSORBER AMOUNTS AND OTHER USEFUL QUANTITIES
!                 --------------------------------------------

IF (LHOOK) CALL DR_HOOK('SW',0,ZHOOK_HANDLE)
LLDEBUG=.FALSE.
CALL SWU ( KIDIA,KFDIA ,KLON  ,KLEV,&
 & PSCT ,PCARDI,PCLDSW,PPMB ,PPSOL,&
 & PRMU0,PTAVE ,PWV,&
 & ZAKI ,ZCLD  ,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD )  

!     ------------------------------------------------------------------
!*         2.     INTERVAL (0.185/0.25-0.68 MICRON): U.V. AND VISIBLE
!                 ---------------------------------------------------
IF (NSW <= 4) THEN
  INUVS=1
  INUIR=2
ELSEIF (NSW == 6) THEN
  INUVS=1
  INUIR=4
ENDIF     

DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    ZFD(JL,JK) =0.0_JPRB
    ZFU(JL,JK) =0.0_JPRB
    ZCD(JL,JK) =0.0_JPRB
    ZCU(JL,JK) =0.0_JPRB
  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZSUDU1T(JL)=0.0_JPRB
  PUVDF(JL)  =0.0_JPRB
  PPARF(JL)  =0.0_JPRB
  PPARCF(JL) =0.0_JPRB
ENDDO

IF(LLDEBUG) THEN
call writefield_phy('sw_zsec',ZSEC,1)
call writefield_phy('sw_zrmu',ZRMU,1)
call writefield_phy('sw_prmu0',PRMU0,1)
call writefield_phy('sw_zfact',ZFACT,1)
ENDIF

icount=icount+1
DO JNU = INUVS , INUIR-1
   !++MODIFCODE
     CALL SW1S &
           &( KIDIA , KFDIA, KLON , KLEV , KAER  , JNU &
           &,  PAER , PALBD , PALBP, PCG  , ZCLD , ZCLEAR &
           &,  ZDSIG, POMEGA, POZ  , ZRMU , ZSEC , PTAU  , ZUD  &
           &,  ZFDUVS,ZFUUVS, ZCDUVS,ZCUUVS, ZSUDU1, ZDIFF,ZDIRF &
           &,  LRDUST,PPIZA_DST(:,:,JNU) &       ! SSA for this wavelength
           &,  PCGA_DST(:,:,JNU)   &            ! GCA for this wavelengt
           &,  PTAUREL_DST(:,:,JNU) )           ! TAUREL for this wavelength
  !--MODIFCODE
IF(LLDEBUG) THEN
! Ecriture des champs avec un indicage du fichier par l'intervalle spectral
  write(str1,'(i1)') jnu
  call writefield_phy("sw_zcduvs"//str1,zcduvs,klev+1)
ENDIF


  DO JL=KIDIA,KFDIA
  PDIFFS(JL,JNU)=ZDIFF(JL,1)*ZFACT(JL)
  PDIRFS(JL,JNU)=ZDIRF(JL,1)*ZFACT(JL)
  ENDDO
  DO JK = 1 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZFD(JL,JK)=ZFD(JL,JK)+ZFDUVS(JL,JK)
      ZFU(JL,JK)=ZFU(JL,JK)+ZFUUVS(JL,JK)
      ZCD(JL,JK)=ZCD(JL,JK)+ZCDUVS(JL,JK)
      ZCU(JL,JK)=ZCU(JL,JK)+ZCUUVS(JL,JK)
    ENDDO
  ENDDO
  DO JL = KIDIA,KFDIA
    ZSUDU1T(JL)=ZSUDU1T(JL)+ZSUDU1(JL)
  ENDDO

  IF (NSW == 6) THEN
    IF (JNU <= 2) THEN
      DO JL = KIDIA,KFDIA
        PUVDF(JL)=PUVDF(JL)+ZFDUVS(JL,1)
      ENDDO
    ELSEIF (JNU == 3) THEN
      DO JL=KIDIA,KFDIA
        PPARF(JL)=PPARF(JL)+ZFDUVS(JL,1)
        PPARCF(JL)=PPARCF(JL)+ZCDUVS(JL,1)
      ENDDO
    ENDIF    
  ENDIF  
ENDDO

!if (icount==5) stop'on arrete dans sw.F90 au bout de 5 appels'
!     ------------------------------------------------------------------

!*         3.     INTERVAL (0.68-4.00 MICRON): NEAR-INFRARED
!                 ------------------------------------------

DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    ZFDOWN(JL,JK)=0.0_JPRB
    ZFUP  (JL,JK)=0.0_JPRB
    ZCDOWN(JL,JK)=0.0_JPRB
    ZCUP  (JL,JK)=0.0_JPRB
    ZSUDU2T(JL)  =0.0_JPRB
  ENDDO
ENDDO

DO JNU = INUIR , NSW
   !++MODIFCODE
      CALL SWNI &
           &(  KIDIA ,KFDIA , KLON , KLEV , KAER , JNU &
           &,  PAER  ,ZAKI  , PALBD, PALBP, PCG  , ZCLD, ZCLEAR &
           &,  ZDSIG ,POMEGA, POZ  , ZRMU , ZSEC , PTAU, ZUD      &
           &,  PWV   ,PQS &
           &,  ZFDNIR,ZFUNIR,ZCDNIR,ZCUNIR,ZSUDU2,ZDIFF2,ZDIRF2 &
           &,  LRDUST,PPIZA_DST(:,:,JNU)  &
           &,  PCGA_DST(:,:,JNU)    &
           &,  PTAUREL_DST(:,:,JNU) &
           &)
    !--MODIFCODE

IF(LLDEBUG) THEN
! Ecriture des champs avec un indicage du fichier par l'intervalle spectral
  write(str1,'(i1)') jnu
  call writefield_phy("sw_zcdnir"//str1,zcdnir,klev+1)
ENDIF

  DO JL=KIDIA,KFDIA
    PDIFFS(JL,JNU)=ZDIFF2(JL,1)*ZFACT(JL)
    PDIRFS(JL,JNU)=ZDIRF2(JL,1)*ZFACT(JL)
  ENDDO
  DO JK = 1 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZFDOWN(JL,JK)=ZFDOWN(JL,JK)+ZFDNIR(JL,JK)
      ZFUP  (JL,JK)=ZFUP  (JL,JK)+ZFUNIR(JL,JK)
      ZCDOWN(JL,JK)=ZCDOWN(JL,JK)+ZCDNIR(JL,JK)
      ZCUP  (JL,JK)=ZCUP  (JL,JK)+ZCUNIR(JL,JK)
    ENDDO
  ENDDO
  DO JL = KIDIA,KFDIA
    ZSUDU2T(JL)=ZSUDU2T(JL)+ZSUDU2(JL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         4.     FILL THE DIAGNOSTIC ARRAYS
!                 --------------------------

DO JL = KIDIA,KFDIA
  PFDNN(JL)=ZFDOWN(JL,1)*ZFACT(JL)
  PFDNV(JL)=ZFD(JL,1)*ZFACT(JL)
  PFUPN(JL)=ZFUP(JL,KLEV+1)*ZFACT(JL)
  PFUPV(JL)=ZFU(JL,KLEV+1)*ZFACT(JL)

  PCDNN(JL)=ZCDOWN(JL,1)*ZFACT(JL)
  PCDNV(JL)=ZCD(JL,1)*ZFACT(JL)
  PCUPN(JL)=ZCUP(JL,KLEV+1)*ZFACT(JL)
  PCUPV(JL)=ZCU(JL,KLEV+1)*ZFACT(JL)

  PSUDU(JL)=(ZSUDU1T(JL)+ZSUDU2T(JL))*ZFACT(JL)
  PUVDF(JL)=PUVDF(JL)*ZFACT(JL)
  PPARF(JL)=PPARF(JL)*ZFACT(JL)
  PPARCF(JL)=PPARCF(JL)*ZFACT(JL)
ENDDO

!WRITE(*,'("---> Dans SW:")')
!WRITE(*,'("ZFUP  ",10E12.5)') (ZFUP(1,JK),JK=1,KLEV+1)
!WRITE(*,'("ZFU   ",10E12.5)') (ZFU(1,JK),JK=1,KLEV+1)
!WRITE(*,'("ZFUNIR",10E12.5)') (ZFUNIR(1,JK),JK=1,KLEV+1)
!WRITE(*,'("ZFACT ",E12.5)') ZFACT(1)

DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PFUP(JL,JK)   = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
    PFDOWN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
    PCUP(JL,JK)   = (ZCUP(JL,JK)   + ZCU(JL,JK)) * ZFACT(JL)
    PCDOWN(JL,JK) = (ZCDOWN(JL,JK) + ZCD(JL,JK)) * ZFACT(JL)
  ENDDO
ENDDO
IF(LLDEBUG) THEN
call writefield_phy('sw_pcdown',PCDOWN,KLEV+1)
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SW',1,ZHOOK_HANDLE)
END SUBROUTINE SW
