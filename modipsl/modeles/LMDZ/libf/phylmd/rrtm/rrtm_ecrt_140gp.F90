!
! $Id: rrtm_ecrt_140gp.F90 2626 2016-09-15 14:20:56Z musat $
!
!****************** SUBROUTINE RRTM_ECRT_140GP **************************

SUBROUTINE RRTM_ECRT_140GP &
 & ( K_IPLON, klon , klev, kcld,&
 & paer , paph , pap,&
 & pts  , pth  , pt,&
 & P_ZEMIS, P_ZEMIW,&
 & pq   , pcco2, pozn, pcldf, ptaucld, ptclear,&
 & P_CLDFRAC,P_TAUCLD,&
 & PTAU_LW,&
 & P_COLDRY,P_WKL,P_WX,&
 & P_TAUAERL,PAVEL,P_TAVEL,PZ,P_TZ,P_TBOUND,K_NLAYERS,P_SEMISS,K_IREFLECT )  

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714

!     Read in atmospheric profile from ECMWF radiation code, and prepare it
!     for use in RRTM.  Set other RRTM input parameters.  Values are passed
!     back through existing RRTM arrays and commons.

!- Modifications

!     2000-05-15 Deborah Salmond  Speed-up

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARRRTM  , ONLY : JPBAND   ,JPXSEC   ,JPLAY   ,&
 & JPINPX  
USE YOERAD   , ONLY : NLW      ,NOVLP
!MPL/IM 20160915 on prend GES de phylmd USE YOERDI   , ONLY :    RCH4     ,RN2O    ,RCFC11  ,RCFC12
USE YOESW    , ONLY : RAER

!------------------------------Arguments--------------------------------

IMPLICIT NONE


INTEGER(KIND=JPIM),INTENT(IN)    :: KLON! Number of atmospheres (longitudes) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV! Number of atmospheric layers 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IPLON 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCLD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) ! Aerosol optical thickness
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) ! Interface pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) ! Layer pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) ! Surface temperature (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1) ! Interface temperatures (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) ! Layer temperature (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ZEMIS(KLON) ! Non-window surface emissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ZEMIW(KLON) ! Window surface emissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) ! H2O specific humidity (mmr)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCO2 ! CO2 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZN(KLON,KLEV) ! O3 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDF(KLON,KLEV) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUCLD(KLON,KLEV,JPBAND) ! Cloud optical depth
!--C.Kleinschmitt
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU_LW(KLON,KLEV,NLW) ! LW Optical depth of aerosols  
!--end
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTCLEAR 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_CLDFRAC(JPLAY) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUCLD(JPLAY,JPBAND) ! Spectral optical thickness
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLDRY(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_WKL(JPINPX,JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_WX(JPXSEC,JPLAY) ! Amount of trace gases
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUAERL(JPLAY,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAVEL(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAVEL(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ(0:JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TZ(0:JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TBOUND 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_NLAYERS 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SEMISS(JPBAND) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_IREFLECT 
!      real rch4                       ! CH4 mass mixing ratio
!      real rn2o                       ! N2O mass mixing ratio
!      real rcfc11                     ! CFC11 mass mixing ratio
!      real rcfc12                     ! CFC12 mass mixing ratio
!- from AER
!- from PROFILE             
!- from SURFACE             
REAL(KIND=JPRB) :: ztauaer(5)
REAL(KIND=JPRB) :: zc1j(0:klev)               ! total cloud from top and level k
REAL(KIND=JPRB) :: Z_AMD                  ! Effective molecular weight of dry air (g/mol)
REAL(KIND=JPRB) :: Z_AMW                  ! Molecular weight of water vapor (g/mol)
REAL(KIND=JPRB) :: Z_AMCO2                ! Molecular weight of carbon dioxide (g/mol)
REAL(KIND=JPRB) :: Z_AMO                  ! Molecular weight of ozone (g/mol)
REAL(KIND=JPRB) :: Z_AMCH4                ! Molecular weight of methane (g/mol)
REAL(KIND=JPRB) :: Z_AMN2O                ! Molecular weight of nitrous oxide (g/mol)
REAL(KIND=JPRB) :: Z_AMC11                ! Molecular weight of CFC11 (g/mol) - CFCL3
REAL(KIND=JPRB) :: Z_AMC12                ! Molecular weight of CFC12 (g/mol) - CF2CL2
REAL(KIND=JPRB) :: Z_AVGDRO               ! Avogadro's number (molecules/mole)
REAL(KIND=JPRB) :: Z_GRAVIT               ! Gravitational acceleration (cm/sec2)

! Atomic weights for conversion from mass to volume mixing ratios; these
!  are the same values used in ECRT to assure accurate conversion to vmr
data Z_AMD   /  28.970_JPRB    /
data Z_AMW   /  18.0154_JPRB   /
data Z_AMCO2 /  44.011_JPRB    /
data Z_AMO   /  47.9982_JPRB   /
data Z_AMCH4 /  16.043_JPRB    /
data Z_AMN2O /  44.013_JPRB    /
data Z_AMC11 / 137.3686_JPRB   /
data Z_AMC12 / 120.9140_JPRB   /
data Z_AVGDRO/ 6.02214E23_JPRB /
data Z_GRAVIT/ 9.80665E02_JPRB /

INTEGER(KIND=JPIM) :: IATM, IMOL, IXMAX, J1, J2, JAE, JB, JK, JL, I_L
INTEGER(KIND=JPIM) :: I_NMOL, I_NXMOL

REAL(KIND=JPRB) :: Z_AMM, ZCLDLY, ZCLEAR, ZCLOUD, ZEPSEC
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!MPL/IM 20160915 on prend GES de phylmd
#include "clesphys.h"
! ***

! *** mji
! Initialize all molecular amounts and aerosol optical depths to zero here, 
! then pass ECRT amounts into RRTM arrays below.

!      DATA ZWKL /MAXPRDW*0.0/
!      DATA ZWX  /MAXPROD*0.0/
!      DATA KREFLECT /0/

! Activate cross section molecules:
!     NXMOL     - number of cross-sections input by user
!     IXINDX(I) - index of cross-section molecule corresponding to Ith
!                 cross-section specified by user
!                 = 0 -- not allowed in RRTM
!                 = 1 -- CCL4
!                 = 2 -- CFC11
!                 = 3 -- CFC12
!                 = 4 -- CFC22
!      DATA KXMOL  /2/
!      DATA KXINDX /0,2,3,0,31*0/

!      IREFLECT=KREFLECT
!      NXMOL=KXMOL

IF (LHOOK) CALL DR_HOOK('RRTM_ECRT_140GP',0,ZHOOK_HANDLE)
K_IREFLECT=0
I_NXMOL=2

DO J1=1,35
! IXINDX(J1)=0
  DO J2=1,KLEV
    P_WKL(J1,J2)=0.0_JPRB 
  ENDDO
ENDDO
!IXINDX(2)=2
!IXINDX(3)=3

!     Set parameters needed for RRTM execution:
IATM    = 0
!      IXSECT  = 1
!      NUMANGS = 0
!      IOUT    = -1
IXMAX   = 4

!     Bands 6,7,8 are considered the 'window' and allowed to have a
!     different surface emissivity (as in ECMWF).  Eli wrote this part....
P_SEMISS(1)  = P_ZEMIS(K_IPLON)
P_SEMISS(2)  = P_ZEMIS(K_IPLON)
P_SEMISS(3)  = P_ZEMIS(K_IPLON)
P_SEMISS(4)  = P_ZEMIS(K_IPLON)
P_SEMISS(5)  = P_ZEMIS(K_IPLON)
P_SEMISS(6)  = P_ZEMIW(K_IPLON)
P_SEMISS(7)  = P_ZEMIW(K_IPLON)
P_SEMISS(8)  = P_ZEMIW(K_IPLON)
P_SEMISS(9)  = P_ZEMIS(K_IPLON)
P_SEMISS(10) = P_ZEMIS(K_IPLON)
P_SEMISS(11) = P_ZEMIS(K_IPLON)
P_SEMISS(12) = P_ZEMIS(K_IPLON)
P_SEMISS(13) = P_ZEMIS(K_IPLON)
P_SEMISS(14) = P_ZEMIS(K_IPLON)
P_SEMISS(15) = P_ZEMIS(K_IPLON)
P_SEMISS(16) = P_ZEMIS(K_IPLON)

!     Set surface temperature.  

P_TBOUND = pts(K_IPLON)

!     Install ECRT arrays into RRTM arrays for pressure, temperature,
!     and molecular amounts.  Pressures are converted from Pascals
!     (ECRT) to mb (RRTM).  H2O, CO2, O3 and trace gas amounts are 
!     converted from mass mixing ratio to volume mixing ratio.  CO2
!     converted with same dry air and CO2 molecular weights used in 
!     ECRT to assure correct conversion back to the proper CO2 vmr.
!     The dry air column COLDRY (in molec/cm2) is calculated from 
!     the level pressures PZ (in mb) based on the hydrostatic equation
!     and includes a correction to account for H2O in the layer.  The
!     molecular weight of moist air (amm) is calculated for each layer.
!     Note: RRTM levels count from bottom to top, while the ECRT input
!     variables count from the top down and must be reversed here.

K_NLAYERS = klev
I_NMOL = 6
PZ(0) = paph(K_IPLON,klev+1)/100._JPRB
P_TZ(0) = pth(K_IPLON,klev+1)
DO I_L = 1, KLEV
  PAVEL(I_L) = pap(K_IPLON,KLEV-I_L+1)/100._JPRB
  P_TAVEL(I_L) = pt(K_IPLON,KLEV-I_L+1)
  PZ(I_L) = paph(K_IPLON,KLEV-I_L+1)/100._JPRB
  P_TZ(I_L) = pth(K_IPLON,KLEV-I_L+1)
  P_WKL(1,I_L) = pq(K_IPLON,KLEV-I_L+1)*Z_AMD/Z_AMW
  P_WKL(2,I_L) = pcco2*Z_AMD/Z_AMCO2
  P_WKL(3,I_L) = pozn(K_IPLON,KLEV-I_L+1)*Z_AMD/Z_AMO
  P_WKL(4,I_L) = rn2o*Z_AMD/Z_AMN2O
  P_WKL(6,I_L) = rch4*Z_AMD/Z_AMCH4
  Z_AMM = (1-P_WKL(1,I_L))*Z_AMD + P_WKL(1,I_L)*Z_AMW
  P_COLDRY(I_L) = (PZ(I_L-1)-PZ(I_L))*1.E3_JPRB*Z_AVGDRO/(Z_GRAVIT*Z_AMM*(1+P_WKL(1,I_L)))
ENDDO

!- Fill RRTM aerosol arrays with operational ECMWF aerosols,
!  do the mixing and distribute over the 16 spectral intervals

DO I_L=1,KLEV
  JK=KLEV-I_L+1
!       DO JAE=1,5
  JAE=1
  ZTAUAER(JAE) =&
   & RAER(JAE,1)*PAER(K_IPLON,1,JK)+RAER(JAE,2)*PAER(K_IPLON,2,JK)&
   & +RAER(JAE,3)*PAER(K_IPLON,3,JK)+RAER(JAE,4)*PAER(K_IPLON,4,JK)&
   & +RAER(JAE,5)*PAER(K_IPLON,5,JK)+RAER(JAE,6)*PAER(K_IPLON,6,JK)  
  P_TAUAERL(I_L, 1)=ZTAUAER(1)
  P_TAUAERL(I_L, 2)=ZTAUAER(1)
  JAE=2
  ZTAUAER(JAE) =&
   & RAER(JAE,1)*PAER(K_IPLON,1,JK)+RAER(JAE,2)*PAER(K_IPLON,2,JK)&
   & +RAER(JAE,3)*PAER(K_IPLON,3,JK)+RAER(JAE,4)*PAER(K_IPLON,4,JK)&
   & +RAER(JAE,5)*PAER(K_IPLON,5,JK)+RAER(JAE,6)*PAER(K_IPLON,6,JK)  
  P_TAUAERL(I_L, 3)=ZTAUAER(2)
  P_TAUAERL(I_L, 4)=ZTAUAER(2)
  P_TAUAERL(I_L, 5)=ZTAUAER(2)
  JAE=3
  ZTAUAER(JAE) =&
   & RAER(JAE,1)*PAER(K_IPLON,1,JK)+RAER(JAE,2)*PAER(K_IPLON,2,JK)&
   & +RAER(JAE,3)*PAER(K_IPLON,3,JK)+RAER(JAE,4)*PAER(K_IPLON,4,JK)&
   & +RAER(JAE,5)*PAER(K_IPLON,5,JK)+RAER(JAE,6)*PAER(K_IPLON,6,JK)  
  P_TAUAERL(I_L, 6)=ZTAUAER(3)
  P_TAUAERL(I_L, 8)=ZTAUAER(3)
  P_TAUAERL(I_L, 9)=ZTAUAER(3)
  JAE=4
  ZTAUAER(JAE) =&
   & RAER(JAE,1)*PAER(K_IPLON,1,JK)+RAER(JAE,2)*PAER(K_IPLON,2,JK)&
   & +RAER(JAE,3)*PAER(K_IPLON,3,JK)+RAER(JAE,4)*PAER(K_IPLON,4,JK)&
   & +RAER(JAE,5)*PAER(K_IPLON,5,JK)+RAER(JAE,6)*PAER(K_IPLON,6,JK)  
  P_TAUAERL(I_L, 7)=ZTAUAER(4)
  JAE=5
  ZTAUAER(JAE) =&
   & RAER(JAE,1)*PAER(K_IPLON,1,JK)+RAER(JAE,2)*PAER(K_IPLON,2,JK)&
   & +RAER(JAE,3)*PAER(K_IPLON,3,JK)+RAER(JAE,4)*PAER(K_IPLON,4,JK)&
   & +RAER(JAE,5)*PAER(K_IPLON,5,JK)+RAER(JAE,6)*PAER(K_IPLON,6,JK)  
!       END DO
  P_TAUAERL(I_L,10)=ZTAUAER(5)
  P_TAUAERL(I_L,11)=ZTAUAER(5)
  P_TAUAERL(I_L,12)=ZTAUAER(5)
  P_TAUAERL(I_L,13)=ZTAUAER(5)
  P_TAUAERL(I_L,14)=ZTAUAER(5)
  P_TAUAERL(I_L,15)=ZTAUAER(5)
  P_TAUAERL(I_L,16)=ZTAUAER(5)
ENDDO
!--Use LW AOD from own Mie calculations (C. Kleinschmitt)
DO I_L=1,KLEV
  JK=KLEV-I_L+1
  DO JAE=1, NLW 
    P_TAUAERL(I_L,JAE) = MAX( PTAU_LW(K_IPLON, JK, JAE), 1e-30 )
  ENDDO
ENDDO
!--end C. Kleinschmitt

DO J2=1,KLEV
  DO J1=1,JPXSEC
    P_WX(J1,J2)=0.0_JPRB
  ENDDO
ENDDO

DO I_L = 1, KLEV
!- Set cross section molecule amounts from ECRT; convert to vmr
  P_WX(2,I_L) = rcfc11*Z_AMD/Z_AMC11
  P_WX(3,I_L) = rcfc12*Z_AMD/Z_AMC12
  P_WX(2,I_L) = P_COLDRY(I_L) * P_WX(2,I_L) * 1.E-20_JPRB
  P_WX(3,I_L) = P_COLDRY(I_L) * P_WX(3,I_L) * 1.E-20_JPRB

!- Here, all molecules in WKL and WX are in volume mixing ratio; convert to
!  molec/cm2 based on COLDRY for use in RRTM

  DO IMOL = 1, I_NMOL
    P_WKL(IMOL,I_L) = P_COLDRY(I_L) * P_WKL(IMOL,I_L)
  ENDDO  
  
! DO IX = 1,JPXSEC
! IF (IXINDX(IX)  /=  0) THEN
!     WX(IXINDX(IX),L) = COLDRY(L) * WX(IX,L) * 1.E-20_JPRB
! ENDIF
! END DO  

ENDDO

!- Approximate treatment for various cloud overlaps
ZCLEAR=1.0_JPRB
ZCLOUD=0.0_JPRB
ZC1J(0)=0.0_JPRB
ZEPSEC=1.E-03_JPRB
JL=K_IPLON

!++MODIFCODE
IF ((NOVLP == 1).OR.(NOVLP ==6).OR.(NOVLP ==8)) THEN
!--MODIFCODE

  DO JK=1,KLEV
    IF (pcldf(JL,JK) > ZEPSEC) THEN
      ZCLDLY=pcldf(JL,JK)
      ZCLEAR=ZCLEAR &
       & *(1.0_JPRB-MAX( ZCLDLY , ZCLOUD ))&
       & /(1.0_JPRB-MIN( ZCLOUD , 1.0_JPRB-ZEPSEC ))  
      ZCLOUD = ZCLDLY
      ZC1J(JK)= 1.0_JPRB - ZCLEAR
    ELSE
      ZCLDLY=0.0_JPRB
      ZCLEAR=ZCLEAR &
       & *(1.0_JPRB-MAX( ZCLDLY , ZCLOUD ))&
       & /(1.0_JPRB-MIN( ZCLOUD , 1.0_JPRB-ZEPSEC ))  
      ZCLOUD = ZCLDLY
      ZC1J(JK)= 1.0_JPRB - ZCLEAR
    ENDIF
  ENDDO

!++MODIFCODE
ELSEIF ((NOVLP == 2).OR.(NOVLP ==7)) THEN
!--MODIFCODE

  DO JK=1,KLEV
    IF (pcldf(JL,JK) > ZEPSEC) THEN
      ZCLDLY=pcldf(JL,JK)
      ZCLOUD = MAX( ZCLDLY , ZCLOUD )
      ZC1J(JK) = ZCLOUD
    ELSE
      ZCLDLY=0.0_JPRB
      ZCLOUD = MAX( ZCLDLY , ZCLOUD )
      ZC1J(JK) = ZCLOUD
    ENDIF
  ENDDO

!++MODIFCODE
ELSEIF ((NOVLP == 3).OR.(NOVLP ==5)) THEN
!--MODIFCODE

  DO JK=1,KLEV
    IF (pcldf(JL,JK) > ZEPSEC) THEN
      ZCLDLY=pcldf(JL,JK)
      ZCLEAR = ZCLEAR * (1.0_JPRB-ZCLDLY)
      ZCLOUD = 1.0_JPRB - ZCLEAR
      ZC1J(JK) = ZCLOUD
    ELSE
      ZCLDLY=0.0_JPRB
      ZCLEAR = ZCLEAR * (1.0_JPRB-ZCLDLY)
      ZCLOUD = 1.0_JPRB - ZCLEAR
      ZC1J(JK) = ZCLOUD
    ENDIF
  ENDDO

ELSEIF (NOVLP == 4) THEN

ENDIF
PTCLEAR=1.0_JPRB-ZC1J(KLEV)

! Transfer cloud fraction and cloud optical depth to RRTM arrays; 
! invert array index for pcldf to go from bottom to top for RRTM

!- clear-sky column
IF (PTCLEAR  >  1.0_JPRB-ZEPSEC) THEN
  KCLD=0
  DO I_L = 1, KLEV
    P_CLDFRAC(I_L) = 0.0_JPRB
  ENDDO
  DO JB=1,JPBAND
    DO I_L=1,KLEV
      P_TAUCLD(I_L,JB) = 0.0_JPRB
    ENDDO
  ENDDO

ELSE

!- cloudy column
!   The diffusivity factor (Savijarvi, 1997) on the cloud optical 
!   thickness TAUCLD has already been applied in RADLSW

  KCLD=1
  DO I_L=1,KLEV
    P_CLDFRAC(I_L) = pcldf(K_IPLON,I_L)
  ENDDO
  DO JB=1,JPBAND
    DO I_L=1,KLEV
      P_TAUCLD(I_L,JB) = ptaucld(K_IPLON,I_L,JB)
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------
 
IF (LHOOK) CALL DR_HOOK('RRTM_ECRT_140GP',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_ECRT_140GP
