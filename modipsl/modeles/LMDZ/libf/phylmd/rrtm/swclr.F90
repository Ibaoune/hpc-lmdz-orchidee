SUBROUTINE SWCLR &
 & ( KIDIA , KFDIA , KLON  , KLEV  , KAER  , KNU,&
 & PAER  , PALBP , PDSIG , PRAYL , PSEC,&
 & PCGAZ , PPIZAZ, PRAY1 , PRAY2 , PREFZ , PRJ,&
 & PRK   , PRMU0 , PTAUAZ, PTRA1 , PTRA2 , PTRCLR, &
!++MODIFCODE
  & LDDUST,PPIZA_DST, PCGA_DST, PTAU_DST )
!--MODIFCODE

!**** *SWCLR* - CLEAR-SKY COLUMN COMPUTATIONS

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
!     CLEAR-SKY COLUMN

!**   INTERFACE.
!     ----------

!          *SWCLR* IS CALLED EITHER FROM *SW1S*
!                                OR FROM *SWNI*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 94-11-15
!        Modified : 96-03-19 JJM-PhD (loop 107 in absence of aerosols)
!        JJMorcrette 990128 : sunshine duration
!        JJMorcrette 990128 : sunshine duration
!        99-05-25   JJMorcrette    Revised aerosols
!        JJMorcrette 001218 : 6 spectral intervals
!        03-10-10 Deborah Salmond and Marta Janiskova Optimisation
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A.Grini (Meteo-France: 2005-11-10) 
!        Y.Seity 05-10-10 : add add 3 optional arg. for dust SW properties
!        Y.Seity 06-09-09 : add modset from O.Thouron (MesoNH) under NOVLP tests
!        O.Boucher fev.2014: modification sur les aerosols pour utiliser les variables DST
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOESW    , ONLY : RTAUA    ,RPIZA    ,RCGA
!USE YOERAD   , ONLY : NOVLP    ,NSW
! NSW mis dans .def MPL 20140211
USE YOERAD   , ONLY : NOVLP    
USE YOERDI   , ONLY : REPCLC
USE YOERDU   , ONLY : REPSCT

IMPLICIT NONE
INCLUDE "clesphys.h"

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAER 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSIG(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAYL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEC(KLON) 
!++MODIFCODE
LOGICAL           ,INTENT(IN)    :: LDDUST                   ! flag for DUST
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZA_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGA_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU_DST(KLON,KLEV)
!--MODIFCODE
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCGAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPIZAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAY1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAY2(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREFZ(KLON,2,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRJ(KLON,6,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRK(KLON,6,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRMU0(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAUAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRA1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRA2(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRCLR(KLON) 
!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!     ------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZC0I(KLON,KLEV+1)&
 & ,  ZCLE0(KLON,KLEV), ZCLEAR(KLON) &
 & ,  ZR21(KLON)&
 & ,  ZR23(KLON) , ZSS0(KLON) , ZSCAT(KLON)&
 & ,  ZTR(KLON,2,KLEV+1)  

INTEGER(KIND=JPIM) :: IKL, JA, JAE, JAJ, JK, JKL, JKLP1, JKM1, JL, INU1

REAL(KIND=JPRB) :: ZBMU0, ZBMU1, ZCORAE, ZDEN, ZDEN1, ZFACOA,&
 & ZFF, ZGAP, ZGAR, ZMU1, ZMUE, ZRATIO, ZRE11, &
 & ZTO, ZTRAY, ZWW, ZDENB   
REAL(KIND=JPRB) :: ZRR,ZMU0,ZI2MU1,ZIMU1,ZIDEN,ZIDEN1
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!++MODIFCODE
REAL(KIND=JPRB) ::ZFACOA_NEW(KLON,KLEV)
!--MODIFCODE


!     ------------------------------------------------------------------

!*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
!                --------------------------------------------

IF (LHOOK) CALL DR_HOOK('SWCLR',0,ZHOOK_HANDLE)
DO JK = 1 , KLEV+1
  DO JA = 1 , 6
    DO JL = KIDIA,KFDIA
      PRJ(JL,JA,JK) = 0.0_JPRB
      PRK(JL,JA,JK) = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO

! ------   NB: 'PAER' AEROSOLS ARE ENTERED FROM TOP TO BOTTOM

DO JK = 1 , KLEV
  IKL=KLEV+1-JK
  DO JL = KIDIA,KFDIA
    PCGAZ(JL,JK) = 0.0_JPRB
    PPIZAZ(JL,JK) =  0.0_JPRB
    PTAUAZ(JL,JK) = 0.0_JPRB
    ZFACOA_NEW(JL,JK) = 0.0_JPRB
  ENDDO

!++MODIFCODE  
!--OB on fait passer les aerosols LMDZ dans la variable DST
  IF(NOVLP < 5)THEN !ECMWF VERSION
!  DO JAE=1,6
      DO JL = KIDIA,KFDIA
!        PTAUAZ(JL,JK)=PTAUAZ(JL,JK)+PAER(JL,JAE,IKL)*RTAUA(KNU,JAE)
        PTAUAZ(JL,JK)=PTAU_DST(JL,IKL)
!        PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL,JAE,IKL)&
!         & * RTAUA(KNU,JAE)*RPIZA(KNU,JAE)  
        PPIZAZ(JL,JK)=PTAU_DST(JL,IKL)*PPIZA_DST(JL,IKL)
!        PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL,JAE,IKL)&
!         & * RTAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)  
        PCGAZ(JL,JK)=PTAU_DST(JL,IKL)*PPIZA_DST(JL,IKL)*PCGA_DST(JL,IKL)
      ENDDO
!    ENDDO
  ELSE ! MESONH VERSION
!--OB on utilise directement les aerosols LMDZ
!     DO JAE=1,6
        DO JL = KIDIA,KFDIA
           !Special optical properties for dust
!           IF (LDDUST.AND.(JAE==3)) THEN
           !Ponderation of aerosol optical properties:first step 
           !ti
!            PTAUAZ(JL,JK)=PTAUAZ(JL,JK) + PAER(JL,JAE,IKL) * PTAUREL_DST(JL,IKL)
            PTAUAZ(JL,JK)= PTAU_DST(JL,IKL)
           !wi*ti
!             PPIZAZ(JL,JK)=PPIZAZ(JL,JK) + PAER(JL,JAE,IKL)  &
!                   & *PTAUREL_DST(JL,IKL)*PPIZA_DST(JL,IKL)
             PPIZAZ(JL,JK)=PTAU_DST(JL,IKL)*PPIZA_DST(JL,IKL)
           !wi*ti*gi
!             PCGAZ(JL,JK) = PCGAZ(JL,JK) + PAER(JL,JAE,IKL) &
!                &  *PTAUREL_DST(JL,IKL)*PPIZA_DST(JL,IKL)*PCGA_DST(JL,IKL)
             PCGAZ(JL,JK) = PTAU_DST(JL,IKL)*PPIZA_DST(JL,IKL)*PCGA_DST(JL,IKL)
           !wi*ti*(gi**2)
!             ZFACOA_NEW(JL,JK)= ZFACOA_NEW(JL,JK)+PAER(JL, JAE, IKL)&
!                & *PTAUREL_DST(JL,IKL) *PPIZA_DST(JL,IKL)*PCGA_DST(JL,IKL)*&
!                & PCGA_DST(JL,IKL)
             ZFACOA_NEW(JL,JK)= ZFACOA_NEW(JL,JK)+&
                & PTAU_DST(JL,IKL) *PPIZA_DST(JL,IKL)*PCGA_DST(JL,IKL)*&
                & PCGA_DST(JL,IKL)
!           ELSE
           !Ponderation of aerosol optical properties:first step 
           !ti
!             PTAUAZ(JL,JK)=PTAUAZ(JL,JK)+PAER(JL, JAE, IKL)*RTAUA(KNU,JAE)
           !wi*ti
!             PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL, JAE, IKL)&
!                &* RTAUA(KNU,JAE)*RPIZA(KNU,JAE)
           !wi*ti*gi
!             PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL, JAE, IKL)&
!                &* RTAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)
           !wi*ti*(gi**2)
!             ZFACOA_NEW(JL,JK)= ZFACOA_NEW(JL,JK)+PAER(JL, JAE, IKL)&
!                &* RTAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)*RCGA(KNU,JAE)
!           ENDIF
        ENDDO
!     ENDDO
  ENDIF
!--MODIFCODE  

!++MODIFCODE  
  IF (NOVLP < 5) then !ECMWF VERSION
   DO JL = KIDIA,KFDIA
    IF (KAER /= 0) THEN
      PCGAZ(JL,JK)=PCGAZ(JL,JK)/PPIZAZ(JL,JK)
      PPIZAZ(JL,JK)=PPIZAZ(JL,JK)/PTAUAZ(JL,JK)
!!!! wrong  ZRATIO = ZTRAY / (ZTRAY + PTAUAZ(JL,JK))
!--     
      ZGAR = PCGAZ(JL,JK)
      ZFF = ZGAR * ZGAR

!-- bug-fix: ZRATIO must be defined from the transformed value of optical thickness
! MPLFH : ZTRAY N'EST PAS INITIALISE !!!!! A REVOIR (MPL)
      ZTRAY= PRAYL(JL) * PDSIG(JL,JK)
!     print *,'>>>>>>> swclr: ZTRAY ',ZTRAY
      ZDENB = ZTRAY + PTAUAZ(JL,JK)*(1.0_JPRB-PPIZAZ(JL,JK)*ZFF)
      ZRATIO=ZTRAY/ZDENB
 !--     
      PTAUAZ(JL,JK)=ZTRAY+PTAUAZ(JL,JK)*(1.0_JPRB-PPIZAZ(JL,JK)*ZFF)
      PCGAZ(JL,JK) = ZGAR * (1.0_JPRB - ZRATIO) / (1.0_JPRB + ZGAR)
      PPIZAZ(JL,JK) =ZRATIO+(1.0_JPRB-ZRATIO)*PPIZAZ(JL,JK)*(1.0_JPRB-ZFF)&
       & / (1.0_JPRB - PPIZAZ(JL,JK) * ZFF)  
    ELSE
      ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
      PTAUAZ(JL,JK) = ZTRAY
      PCGAZ(JL,JK) = 0.0_JPRB
      PPIZAZ(JL,JK) = 1.0_JPRB-REPSCT
    ENDIF
  ENDDO
  ELSE !MESONH VERSION
   DO JL = KIDIA,KFDIA
    IF (KAER /= 0) THEN
      ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
      ZRATIO =PPIZAZ(JL,JK)+ZTRAY
      !Ponderation G**2
      ZFACOA_NEW(JL,JK)= ZFACOA_NEW(JL,JK)/ZRATIO
      !Ponderation w
      PPIZAZ(JL,JK)=ZRATIO/(PTAUAZ(JL,JK)+ZTRAY)
      !Ponderation g
      PCGAZ(JL,JK)=PCGAZ(JL,JK)/ZRATIO
      !Ponderation+delta-modified parameters tau
      PTAUAZ(JL,JK)=(ZTRAY+PTAUAZ(JL,JK))*&
       &  (1.0_JPRB-PPIZAZ(JL,JK)*ZFACOA_NEW(JL,JK))
      !delta-modified parameters w
      PPIZAZ(JL,JK)=PPIZAZ(JL,JK)*(1.0_JPRB-ZFACOA_NEW(JL,JK))/&
          & (1.0_JPRB-ZFACOA_NEW(JL,JK)*PPIZAZ(JL,JK))     
      !delta-modified parameters g
      PCGAZ(JL,JK)=PCGAZ(JL,JK)/(1.0_JPRB+PCGAZ(JL,JK))
      
    ELSE
      ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
      ZFACOA_NEW(JL,JK)= 0.0_JPRB
      PTAUAZ(JL,JK) = ZTRAY
      PCGAZ(JL,JK) = 0.0_JPRB
      PPIZAZ(JL,JK) = 1.0_JPRB-REPSCT
    ENDIF
   ENDDO    
  ENDIF
!--MODIFCODE  
  
ENDDO

!     ------------------------------------------------------------------

!*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
!                ----------------------------------------------

DO JL = KIDIA,KFDIA
  ZR23(JL) = 0.0_JPRB
  ZC0I(JL,KLEV+1) = 0.0_JPRB
  ZCLEAR(JL) = 1.0_JPRB
  ZSCAT(JL) = 0.0_JPRB
ENDDO

JK = 1
JKL = KLEV+1 - JK
JKLP1 = JKL + 1
DO JL = KIDIA,KFDIA
!++MODIFCODE
  IF (NOVLP >= 5) THEN
   ZFACOA = PTAUAZ(JL,JK)
   ZCORAE = ZFACOA *  PSEC(JL)
  ELSE
   ZFACOA = 1.0_JPRB - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
   ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
  ENDIF
!--MODIFCODE
  ZR21(JL) = EXP(-ZCORAE   )
  ZSS0(JL) = 1.0_JPRB-ZR21(JL)
  ZCLE0(JL,JKL) = ZSS0(JL)

  IF (NOVLP == 1 .OR. NOVLP == 4) THEN
!* maximum-random      
    ZCLEAR(JL) = ZCLEAR(JL)&
     & *(1.0_JPRB-MAX(ZSS0(JL),ZSCAT(JL)))&
     & /(1.0_JPRB-MIN(ZSCAT(JL),1.0_JPRB-REPCLC))  
    ZC0I(JL,JKL) = 1.0_JPRB - ZCLEAR(JL)
    ZSCAT(JL) = ZSS0(JL)
  ELSEIF (NOVLP == 2) THEN
!* maximum
    ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
    ZC0I(JL,JKL) = ZSCAT(JL)
!++MODIFCODE
  ELSEIF ((NOVLP == 3).OR.(NOVLP  >=  5)) THEN
!--MODIFCODE
!* random
    ZCLEAR(JL)=ZCLEAR(JL)*(1.0_JPRB-ZSS0(JL))
    ZSCAT(JL) = 1.0_JPRB - ZCLEAR(JL)
    ZC0I(JL,JKL) = ZSCAT(JL)
  ENDIF
ENDDO

DO JK = 2 , KLEV
  JKL = KLEV+1 - JK
  JKLP1 = JKL + 1
  DO JL = KIDIA,KFDIA
!++MODIFCODE
    IF (NOVLP >= 5) THEN
     ZFACOA = PTAUAZ(JL,JK)
     ZCORAE = ZFACOA *  PSEC(JL)
    ELSE
    ZFACOA = 1.0_JPRB - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
    ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
    ENDIF
!--MODIFCODE
    ZR21(JL) = EXP(-ZCORAE   )
    ZSS0(JL) = 1.0_JPRB-ZR21(JL)
    ZCLE0(JL,JKL) = ZSS0(JL)

    IF (NOVLP == 1 .OR. NOVLP == 4) THEN
!* maximum-random      
      ZCLEAR(JL) = ZCLEAR(JL)&
       & *(1.0_JPRB-MAX(ZSS0(JL),ZSCAT(JL)))&
       & /(1.0_JPRB-MIN(ZSCAT(JL),1.0_JPRB-REPCLC))  
      ZC0I(JL,JKL) = 1.0_JPRB - ZCLEAR(JL)
      ZSCAT(JL) = ZSS0(JL)
    ELSEIF (NOVLP == 2) THEN
!* maximum
      ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
      ZC0I(JL,JKL) = ZSCAT(JL)
!++MODIFCODE
    ELSEIF ((NOVLP == 3).OR.(NOVLP >= 5)) THEN
!--MODIFCODE
!* random
      ZCLEAR(JL)=ZCLEAR(JL)*(1.0_JPRB-ZSS0(JL))
      ZSCAT(JL) = 1.0_JPRB - ZCLEAR(JL)
      ZC0I(JL,JKL) = ZSCAT(JL)
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
!                -----------------------------------------------

DO JL = KIDIA,KFDIA
  PRAY1(JL,KLEV+1) = 0.0_JPRB
  PRAY2(JL,KLEV+1) = 0.0_JPRB
  PREFZ(JL,2,1) = PALBP(JL,KNU)
  PREFZ(JL,1,1) = PALBP(JL,KNU)
  PTRA1(JL,KLEV+1) = 1.0_JPRB
  PTRA2(JL,KLEV+1) = 1.0_JPRB
ENDDO

DO JK = 2 , KLEV+1
  JKM1 = JK-1
  DO JL = KIDIA,KFDIA

!     ------------------------------------------------------------------

!*         3.1  EQUIVALENT ZENITH ANGLE
!               -----------------------

    ZMUE = (1.0_JPRB-ZC0I(JL,JK)) * PSEC(JL)+ ZC0I(JL,JK) * 1.66_JPRB
    PRMU0(JL,JK) = 1.0_JPRB/ZMUE
    ZMU0=PRMU0(JL,JK)

!     ------------------------------------------------------------------

!*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
!               ----------------------------------------------------

    ZGAP = PCGAZ(JL,JKM1)
    ZBMU0 = 0.5_JPRB - 0.75_JPRB * ZGAP *ZMU0
    ZWW = PPIZAZ(JL,JKM1)
    ZTO = PTAUAZ(JL,JKM1)
    ZDEN = 1.0_JPRB + (1.0_JPRB - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE &
     & + (1-ZWW) * (1.0_JPRB - ZWW +2.0_JPRB*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE  
    ZIDEN=1.0_JPRB / ZDEN
    PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE * ZIDEN
    PTRA1(JL,JKM1) = ZIDEN

    ZMU1 = 0.5_JPRB
    ZIMU1=2.0_JPRB
    ZI2MU1=4.0_JPRB
    ZBMU1 = 0.5_JPRB - 0.75_JPRB * ZGAP * ZMU1
    ZDEN1= 1.0_JPRB + (1.0_JPRB - ZWW + ZBMU1 * ZWW) * ZTO * ZIMU1 &
     & + (1-ZWW) * (1.0_JPRB - ZWW +2.0_JPRB*ZBMU1*ZWW)*ZTO*ZTO*ZI2MU1  
    ZIDEN1=1.0_JPRB / ZDEN1
    PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO * ZIMU1 *ZIDEN1
    PTRA2(JL,JKM1) = ZIDEN1

    ZRR=1.0_JPRB/(1.0_JPRB-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1))
    PREFZ(JL,1,JK) = PRAY1(JL,JKM1)&
     & + PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)&
     & * PTRA2(JL,JKM1)&
     & *ZRR  

    ZTR(JL,1,JKM1) = PTRA1(JL,JKM1)&
     & *ZRR  

    PREFZ(JL,2,JK) = PRAY1(JL,JKM1)&
     & + PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)&
     & * PTRA2(JL,JKM1)   

    ZTR(JL,2,JKM1) = PTRA1(JL,JKM1)

  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZMUE = (1.0_JPRB-ZC0I(JL,1))*PSEC(JL)+ZC0I(JL,1)*1.66_JPRB
  PRMU0(JL,1)=1.0_JPRB/ZMUE
  PTRCLR(JL)=1.0_JPRB-ZC0I(JL,1)
ENDDO

!     ------------------------------------------------------------------

!*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!                 -------------------------------------------------

IF (NSW <= 4) THEN
  INU1=1
ELSEIF (NSW == 6) THEN
  INU1=3
ENDIF    

IF (KNU <= INU1) THEN
  JAJ = 2
  DO JL = KIDIA,KFDIA
    PRJ(JL,JAJ,KLEV+1) = 1.0_JPRB
    PRK(JL,JAJ,KLEV+1) = PREFZ(JL, 1,KLEV+1)
  ENDDO

  DO JK = 1 , KLEV
    JKL = KLEV+1 - JK
    JKLP1 = JKL + 1
    DO JL = KIDIA,KFDIA
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
    ENDDO
  ENDDO

ELSE

  DO JAJ = 1 , 2
    DO JL = KIDIA,KFDIA
      PRJ(JL,JAJ,KLEV+1) = 1.0_JPRB
      PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      JKL = KLEV+1 - JK
      JKLP1 = JKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
        PRJ(JL,JAJ,JKL) = ZRE11
        PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
      ENDDO
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SWCLR',1,ZHOOK_HANDLE)
END SUBROUTINE SWCLR
