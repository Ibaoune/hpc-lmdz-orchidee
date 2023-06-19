SUBROUTINE SW1S &
 & ( KIDIA , KFDIA , KLON , KLEV , KAER , KNU,&
 & PAER  , PALBD , PALBP, PCG  , PCLD , PCLEAR,&
 & PDSIG , POMEGA, POZ  , PRMU , PSEC , PTAU  , PUD,&
 & PFD   , PFU   , PCD  , PCU  , PSUDU1,PDIFF , PDIRF, &
!++MODIFCODE
 & LRDUST,PPIZA_DST,PCGA_DST,PTAUREL_DST  &
!--MODIFCODE
 &)

!**** *SW1S* - SHORTWAVE RADIATION, FIRST SPECTRAL INTERVAL

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SW1S* IS CALLED FROM *SW*.

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES QUANTITIES FOR THE CLEAR-SKY FRACTION OF THE
!     COLUMN
!          2. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
!     CONTINUUM SCATTERING
!          3. MULTIPLY BY OZONE TRANSMISSION FUNCTION

!     EXTERNALS.
!     ----------

!          *SWCLR*, *SWR*, *SWTT*, *SWUVO3*

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
!        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
!        96-01-15   J.-J. MORCRETTE    SW in nsw SPECTRAL INTERVALS 
!        990128     JJMorcrette        sunshine duration
!        99-05-25   JJMorcrette        Revised aerosols
!        00-12-18   JJMorcrette        6 spectral intervals
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Y.Seity  04-11-19 : add two arguments for AROME externalized surface
!        Y.Seity  05-10-10 : add 3 optional arg. for dust SW properties
!        Y.Seity 06-09-09 : add modset from O.Thouron (MesoNH) under NOVLP tests
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOESW    , ONLY : RRAY     ,RSUN
!USE YOERAD   , ONLY : NSW
! NSW mis dans .def MPL 20140211
USE write_field_phy

IMPLICIT NONE

include "clesphys.h"

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAER 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KLON,NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCG(KLON,NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLEAR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSIG(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMEGA(KLON,NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEC(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU(KLON,NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUD(KLON,5,KLEV+1) 
!++MODIFCODE
LOGICAL           ,INTENT(IN)    :: LRDUST          ! flag for DUST
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZA_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGA_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUREL_DST(KLON,KLEV)
!--MODIFCODE
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFD(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCD(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSUDU1(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIRF(KLON,KLEV) 
!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!     ------------------------------------------------------------------

!              ------------

INTEGER(KIND=JPIM) :: IIND(6)

REAL(KIND=JPRB) :: ZCGAZ(KLON,KLEV)&
 & ,  ZDIFF(KLON)        , ZDIRF(KLON)        &
 & ,  ZDIFT(KLON)        , ZDIRT(KLON)        &
 & ,  ZPIZAZ(KLON,KLEV)&
 & ,  ZRAYL(KLON), ZRAY1(KLON,KLEV+1), ZRAY2(KLON,KLEV+1)&
 & ,  ZREFZ(KLON,2,KLEV+1)&
 & ,  ZRJ(KLON,6,KLEV+1), ZRJ0(KLON,6,KLEV+1)&
 & ,  ZRK(KLON,6,KLEV+1), ZRK0(KLON,6,KLEV+1)&
 & ,  ZRMUE(KLON,KLEV+1), ZRMU0(KLON,KLEV+1)&
 & ,  ZR(KLON,6)&
 & ,  ZTAUAZ(KLON,KLEV)&
 & ,  ZTRA1(KLON,KLEV+1), ZTRA2(KLON,KLEV+1)&
 & ,  ZTRCLD(KLON)      , ZTRCLR(KLON)&
 & ,  ZW(KLON,6)        , ZO(KLON,2) ,ZT(KLON,2)   

INTEGER(KIND=JPIM) :: IKL, IKM1, JAJ, JK, JL , JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
LOGICAL         :: LLDEBUG

#include "swclr.intfb.h"
#include "swr.intfb.h"
#include "swtt1.intfb.h"
#include "swuvo3.intfb.h"

!     ------------------------------------------------------------------

!*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
!                 ----------------------- ------------------

!*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
!                 -----------------------------------------

IF (LHOOK) CALL DR_HOOK('SW1S',0,ZHOOK_HANDLE)
LLDEBUG=.FALSE.
DO JL = KIDIA,KFDIA
  ZRAYL(JL) =  RRAY(KNU,1) + PRMU(JL) * (RRAY(KNU,2) + PRMU(JL)&
   & * (RRAY(KNU,3) + PRMU(JL) * (RRAY(KNU,4) + PRMU(JL)&
   & * (RRAY(KNU,5) + PRMU(JL) *  RRAY(KNU,6)       ))))  
ENDDO
!     ------------------------------------------------------------------

!*         2.    CONTINUUM SCATTERING CALCULATIONS
!                ---------------------------------

!*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
!                --------------------------------

!++MODIFCODE
CALL SWCLR &
   &( KIDIA  , KFDIA , KLON  , KLEV , KAER , KNU &
   &, PAER   , PALBP , PDSIG , ZRAYL, PSEC &
   &, ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
   &, ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2, ZTRCLR &
   &, LRDUST , PPIZA_DST,PCGA_DST  &
   &, PTAUREL_DST )

!--MODIFCODE

!*         2.2   CLOUDY FRACTION OF THE COLUMN
!                -----------------------------

CALL SWR &
 & ( KIDIA ,KFDIA ,KLON  ,KLEV  , KNU,&
 & PALBD ,PCG   ,PCLD  ,POMEGA, PSEC , PTAU,&
 & ZCGAZ ,ZPIZAZ,ZRAY1 ,ZRAY2 , ZREFZ, ZRJ  ,ZRK , ZRMUE,&
 & ZTAUAZ,ZTRA1 ,ZTRA2 ,ZTRCLD &
 & )  

! DO JK = 1 , KLEV
!   IKL = KLEV+1-JK
!   DO JL = KIDIA,KFDIA
!   print *,'Apres SWCLR,SWR RMU0 RMUE ',ZRMU0(JL,IKL),ZRMUE(JL,IKL)
!   ENDDO
! ENDDO
!     ------------------------------------------------------------------

!*         3.    OZONE ABSORPTION
!                ----------------

IF (NSW <= 4) THEN

!*         3.1   TWO OR FOUR SPECTRAL INTERVALS
!                ------------------------------

  IIND(1)=1
  IIND(2)=2
  IIND(3)=3
  IIND(4)=1
  IIND(5)=2
  IIND(6)=3

!*         3.1.1  DOWNWARD FLUXES
!                 ---------------

  JAJ = 2

  DO JL = KIDIA,KFDIA
    ZW(JL,1)=0.0_JPRB
    ZW(JL,2)=0.0_JPRB
    ZW(JL,3)=0.0_JPRB
    ZW(JL,4)=0.0_JPRB
    ZW(JL,5)=0.0_JPRB
    ZW(JL,6)=0.0_JPRB
    PFD(JL,KLEV+1)=((1.0_JPRB-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)&
     & + PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)  
    PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)
  ENDDO
  DO JK = 1 , KLEV
    IKL = KLEV+1-JK
    DO JL = KIDIA,KFDIA
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
      ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKL)/ZRMUE(JL,IKL)
      ZW(JL,3)=ZW(JL,3)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZW(JL,4)=ZW(JL,4)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW(JL,5)=ZW(JL,5)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
      ZW(JL,6)=ZW(JL,6)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
    ENDDO
    
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, 6,&
     & IIND,&
     & ZW,&
     & ZR                          )  

    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZRJ(JL,JAJ,IKL)
      ZDIRF(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZRJ0(JL,JAJ,IKL)
      PDIFF(JL,IKL) = ZDIFF(JL) * RSUN(KNU)*(1.0_JPRB-PCLEAR(JL))
      PDIRF(JL,IKL) = ZDIRF(JL) * RSUN(KNU)*PCLEAR(JL)
      PFD(JL,IKL) = ((1.0_JPRB-PCLEAR(JL)) * ZDIFF(JL)&
       & +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)  
      PCD(JL,IKL) = ZDIRF(JL) * RSUN(KNU)
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    ZDIFT(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZTRCLD(JL)
    ZDIRT(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZTRCLR(JL)
    PSUDU1(JL) = ((1.0_JPRB-PCLEAR(JL)) * ZDIFT(JL)&
     & +PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)  
  ENDDO

!*         3.1.2  UPWARD FLUXES
!                 -------------

  DO JL = KIDIA,KFDIA
    PFU(JL,1) = ((1.0_JPRB-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)&
     & + PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))&
     & * RSUN(KNU)  
    PCU(JL,1) = ZDIRF(JL) * PALBP(JL,KNU) * RSUN(KNU)
  ENDDO

  DO JK = 2 , KLEV+1
    IKM1=JK-1
    DO JL = KIDIA,KFDIA
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW(JL,3)=ZW(JL,3)+POZ(JL,  IKM1)*1.66_JPRB
      ZW(JL,4)=ZW(JL,4)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW(JL,5)=ZW(JL,5)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW(JL,6)=ZW(JL,6)+POZ(JL,  IKM1)*1.66_JPRB
    ENDDO
    
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, 6,&
     & IIND,&
     & ZW,&
     & ZR                          )  
  
    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZRK(JL,JAJ,JK)
      ZDIRF(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((1.0_JPRB-PCLEAR(JL)) * ZDIFF(JL)&
       & +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)  
      PCU(JL,JK) = ZDIRF(JL) * RSUN(KNU)
    ENDDO
!WRITE(*,'("---> Dans SW1S:")')
!WRITE(*,'("PFU",10E12.5)') (PFU(1,JJ),JJ=1,KLEV+1)
!WRITE(*,'("PCLEAR",10E12.5)') (PCLEAR(1))
!WRITE(*,'("ZDIFF",10E12.5)') (ZDIFF(1))
!WRITE(*,'("ZDIRF",10E12.5)') (ZDIRF(1))
!WRITE(*,'("RSUN",10E12.5)') (RSUN(KNU))
  ENDDO

ELSEIF (NSW == 6) THEN
!print *,'... dans SW1S: NSW=',NSW

!*         3.2   SIX SPECTRAL INTERVALS
!                ----------------------

  IIND(1)=1
  IIND(2)=2
  IIND(3)=1
  IIND(4)=2

!*         3.2,1  DOWNWARD FLUXES
!                 ---------------

  JAJ = 2

  DO JL = KIDIA,KFDIA
    ZW(JL,1)=0.0_JPRB
    ZW(JL,2)=0.0_JPRB
    ZW(JL,3)=0.0_JPRB
    ZW(JL,4)=0.0_JPRB
  
    ZO(JL,1)=0.0_JPRB
    ZO(JL,2)=0.0_JPRB
    PFD(JL,KLEV+1)=((1.0_JPRB-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)&
     & + PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)  
    PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)
  ENDDO
  DO JK = 1 , KLEV
    IKL = KLEV+1-JK
    DO JL = KIDIA,KFDIA
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
      ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKL)/ZRMUE(JL,IKL)
      ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW(JL,4)=ZW(JL,4)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
    
      ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
    ENDDO
 
!   WRITE(*,'("---> Dans SW1S avant SWTT1:")')
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, 4,&
     & IIND,&
     & ZW,&
     & ZR  &
     & )  

!   WRITE(*,'("---> Dans SW1S avant SWUVO3 flux dwn:")')
    CALL SWUVO3 ( KIDIA, KFDIA, KLON, KNU, 2,&
     & ZO,&
     & ZT  &
     & )  

    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZRJ(JL,JAJ,IKL)
      ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZRJ0(JL,JAJ,IKL)
      PDIFF(JL,IKL) = ZDIFF(JL) * RSUN(KNU)*(1.0_JPRB-PCLEAR(JL)) 
      PDIRF(JL,IKL) = ZDIRF(JL) * RSUN(KNU)*PCLEAR(JL)
      PFD(JL,IKL) = ((1.0_JPRB-PCLEAR(JL)) * ZDIFF(JL)&
       & +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)  
      PCD(JL,IKL) = ZDIRF(JL) * RSUN(KNU)
    ENDDO
  ENDDO

  IF(LLDEBUG) THEN
  call writefield_phy('sw1s_pud1',PUD(:,1,:),klev)
  call writefield_phy('sw1s_pud2',PUD(:,2,:),klev)
  call writefield_phy('sw1s_psec',PSEC,1)
  call writefield_phy('sw1s_zrmue',ZRMUE,klev+1)
  call writefield_phy('sw1s_zrmu0',ZRMU0,klev+1)
  call writefield_phy('sw1s_pdirf',PDIRF,klev)
  call writefield_phy('sw1s_pdiff',PDIFF,klev)
  call writefield_phy('sw1s_pfd',PFD,klev)
  ENDIF
  DO JL=KIDIA,KFDIA
    ZDIFT(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZTRCLD(JL)
    ZDIRT(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZTRCLR(JL)
    PSUDU1(JL) = ((1.0_JPRB-PCLEAR(JL)) * ZDIFT(JL)&
     & +PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)  
  ENDDO

!*         3.2.2  UPWARD FLUXES
!                 -------------

  DO JL = KIDIA,KFDIA
    PFU(JL,1) = ((1.0_JPRB-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)&
     & + PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))&
     & * RSUN(KNU)  
    PCU(JL,1) = ZDIRF(JL) * PALBP(JL,KNU) * RSUN(KNU)
  ENDDO

  DO JK = 2 , KLEV+1
    IKM1=JK-1
    DO JL = KIDIA,KFDIA
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW(JL,4)=ZW(JL,4)+PUD(JL,2,IKM1)*1.66_JPRB
      
      ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKM1)*1.66_JPRB
      ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKM1)*1.66_JPRB
    ENDDO

!   WRITE(*,'("---> Dans SW1S avant SWTT1:")')
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, 4,&
     & IIND,&
     & ZW,&
     & ZR  &
     & )  

!   WRITE(*,'("---> Dans SW1S avant SWUVO3 flux up:")')
    CALL SWUVO3 ( KIDIA, KFDIA, KLON, KNU, 2,&
     & ZO,&
     & ZT  &
     & )  

    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZRK(JL,JAJ,JK)
      ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((1.0_JPRB-PCLEAR(JL)) * ZDIFF(JL)&
       & +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)  
      PCU(JL,JK) = ZDIRF(JL) * RSUN(KNU)
!WRITE(*,'("---> Dans SW1S:")')
!print *,'===JL= ',jl
!WRITE(*,'("ZR1",10E12.5)') (ZR(JL,1))
!WRITE(*,'("ZR2",10E12.5)') (ZR(JL,2))
!WRITE(*,'("ZR3",10E12.5)') (ZR(JL,3))
!WRITE(*,'("ZR4",10E12.5)') (ZR(JL,4))
!WRITE(*,'("ZT1",10E12.5)') (ZT(JL,1))
!WRITE(*,'("ZT2",10E12.5)') (ZT(JL,2))
    ENDDO
  ENDDO
  
!WRITE(*,'("---> Dans SW1S:")')
!WRITE(*,'("PFU",10E12.5)') (PFU(1,JJ),JJ=1,KLEV+1)
!WRITE(*,'("ZR",10E12.5)') (ZR(1,JJ),JJ=1,4)
!WRITE(*,'("PCLEAR",10E12.5)') (PCLEAR(1))
!WRITE(*,'("ZDIFF",10E12.5)') (ZDIFF(1))
!WRITE(*,'("ZDIRF",10E12.5)') (ZDIRF(1))
!WRITE(*,'("RSUN",10E12.5)') (RSUN(KNU))
ENDIF  

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SW1S',1,ZHOOK_HANDLE)
END SUBROUTINE SW1S
