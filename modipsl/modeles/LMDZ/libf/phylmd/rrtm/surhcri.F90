!-----------------------------------------------------------------------
      SUBROUTINE SURHCRI(KULOUT)
!-----------------------------------------------------------------------
!
!**** *SURHCRI * - COMPUTATION OF THE CRTITICAL RELATIVE HUMIDITY
!                  PROFILE FOR SMITH'S CONDENSATION SCHEME.
!
!**   Interface.
!     ----------
!        *CALL* *SURHCRI*
!
!-----------------------------------------------------------------------
!
! -   ARGUMENTS D'ENTREE./INPUT ARGUMENTS.
!     ------------------------------------
!
!-----------------------------------------------------------------------
!
! -   ARGUMENTS IMPLICITES.
!     ---------------------
!
! COMMON /YOMPHY/
! COMMON /YOMPHY0/
!
!*
!-----------------------------------------------------------------------
!
!     Auteur.
!     -------
!         05-03, Luc Gerard (from Ph. Lopez acrhcri)
!
!     Modifications.
!     --------------
!         06-10, nettoyage - R. Brozkova
!-----------------------------------------------------------------------
!
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!
USE YOMPHY0   , ONLY : RHCRIT1, RHCRIT2, RETAMIN, GRHCMOD, RHCRI ,NRHCRI 
USE YOMDIM    , ONLY : NFLEVG
USE YOMGC     , ONLY : GAW, GM
USE YOMGEM    , ONLY : VP00, TEQH
USE YOMPHY    , ONLY : LGWRHCRI
USE YOMSTA    , ONLY : STPRE 
 
IMPLICIT NONE
 
INTEGER (KIND=JPIM),INTENT(IN) :: KULOUT
 
! LOCAL
 
REAL(KIND=JPRB) :: ZVETAF
REAL(KIND=JPRB) :: ZGAW0,ZEGAW,ZFACT,ZFACT0,ZFACTA,ZFACTB,ZFACTC, &
        & ZA,ZB,ZC
INTEGER (KIND=JPIM) :: JLEV
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SURHCRI',0,ZHOOK_HANDLE)
 
!     ------------------------------------------------------------------
 
IF(NRHCRI == NFLEVG) THEN
 
!     ------------------------------------------------------------------
 
! CRITICAL RELATIVE HUMIDITY PROFILE
! DEPENDENT ON THE LOCAL RESOLUTION.
! ----------------------------------
 
   ZFACT0=(RETAMIN*(RETAMIN-1._JPRB))**2
   ZFACTA=(2._JPRB*RETAMIN-1._JPRB)
   ZFACTB=(1._JPRB-3._JPRB*RETAMIN**2)
   ZFACTC=(3._JPRB*RETAMIN-2._JPRB)*RETAMIN
   IF (LGWRHCRI) THEN
     ZGAW0=2.42E-05_JPRB
     ZEGAW=(1.0_JPRB-EXP(-GAW(1)/ZGAW0))/ZFACT0
   ELSE
     ZEGAW=(1.0_JPRB-GRHCMOD*EXP(-1.0_JPRB/(TEQH*GM(1))))/ZFACT0
   ENDIF
 
   ZFACT=(RHCRIT2-RHCRIT1)*ZEGAW
   ZA=ZFACTA*ZFACT
   ZB=ZFACTB*ZFACT
   ZC=ZFACTC*ZFACT
 
   DO JLEV=1,NFLEVG
     ZVETAF=STPRE (JLEV)/VP00
     RHCRI(JLEV) = ZA * ZVETAF**3 &
                     & + ZB * ZVETAF**2 &
                     & + ZC * ZVETAF &
                     & + RHCRIT2
   ENDDO   
 
   WRITE(KULOUT,'(A)') 'RHCRI: SMITH SCHEME CRITICAL RH PROFILE:'
   WRITE(KULOUT,'((10(G9.3)))') (RHCRI(JLEV),JLEV=1,NFLEVG)
ELSE
   WRITE(KULOUT,*) 'NRHCRI=',NRHCRI,': GLOBAL MODEL, NO RHCRI INITIALISATION IN SETUP'
ENDIF
IF (LHOOK) CALL DR_HOOK('SURHCRI',1,ZHOOK_HANDLE)
END SUBROUTINE SURHCRI
