!OPTIONS XOPT(NOEVAL)
SUBROUTINE SUPHY3(KULOUT)

!**** *SUPHY3*   - Initialize common YOMPHY3 physics radiative
!                  constants

!     Purpose.
!     --------
!           Initialize YOMPHY3, the common that contains the parameters
!           for the radiation part of the physics of the model.

!**   Interface.
!     ----------
!        *CALL* *SUPHY3(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY3

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

!     Modifications.
!     --------------
!        Original : 90-9-1
!        Ajout de GCE4 et
!         modification de toutes les constantes optiques : 91-2-1
!        Ajout des nuages de glace : 92-04, L. Labbe.
!        Modified by R. EL Khatib : 93-04-02 Set-up defaults controled by LECMWF
!        Modified 97-04-17 J.M. Piriou: radiative coefficients default values, and new QSUSX default.
!        Modified 2000-08 R. Randriamampianina, J.F. Geleyn et J.M. Piriou: change infra-red interaction between layers (GIREC*).
!        Modified 2000-12-12 by E. Bazile: CYCORA's default value.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R.Brozkova    24-Sep-2006 ALARO-0 (new cloud model)
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN   , ONLY : NULNAM
USE YOMCT0B  , ONLY : LECMWF
USE YOMPHY   , ONLY : LCLSATUR
USE YOMPHY3  , ONLY : GCA      ,GCB      ,GCC      ,VDP      ,&
 & VNP      ,BSFSA    ,BSFSI    ,BSFSN    ,BSFTA    ,&
 & BSFTI    ,BSFTN    ,EARRT    ,EOASA    ,EOASI    ,&
 & EOASN    ,EOATA    ,EOATI    ,EOATN    ,EODSA    ,&
 & EODSI    ,EODSN    ,EODTA    ,EODTI    ,EODTN    ,&
 & EORAY    ,GCD4     ,GCE4     ,QCO2     ,QLIMI    ,&
 & QLIP0    ,RII0     ,USAA     ,USAI     ,&
 & USAN     ,USBA     ,USBI     ,USBN     ,&
 & GIREC1   ,GIREC2   ,GIREC3   ,GIREC4   ,&  
 & FCM_DEL_A,FCM_DEL_D,FCM_MU_A ,FCM_MU_D ,FCM_N_I  ,&
 & FCM_N_L  ,FCM_P_AI ,FCM_P_AL ,FCM_P_DI ,FCM_P_DL ,&
 & FCM_P_GI ,FCM_P_GL ,FCM_Q_AI ,FCM_Q_AL ,FCM_Q_DI ,&
 & FCM_Q_DL ,FCM_Q_GI ,FCM_Q_GL ,N_SPBAND ,REXP_NEB

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "posnam.intfb.h"

#include "namphy3.h"
!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

!        1.1 Set implicit default values

IF (LHOOK) CALL DR_HOOK('SUPHY3',0,ZHOOK_HANDLE)
BSFSA=.3490_JPRB
BSFSI=.3238_JPRB
BSFSN=.3300_JPRB
BSFTA=.3471_JPRB
BSFTI=.3231_JPRB
BSFTN=.3584_JPRB
EARRT=0.001324_JPRB
EOASA=.6727E-01_JPRB
EOASI=.1692E-01_JPRB
EOASN=.3353E-01_JPRB
EOATA=.5037E-01_JPRB
EOATI=.3455E+01_JPRB
EOATN=.9077E+01_JPRB
EODSA=.3665E+00_JPRB
EODSI=.1406E+01_JPRB
EODSN=.6558E+01_JPRB
EODTA=.1766E-01_JPRB
EODTI=.8435E+00_JPRB
EODTN=.4319E+01_JPRB
EORAY=.8606E-06_JPRB
GCA(1)=.8041E-01_JPRB
GCA(2)=.1456E+00_JPRB
GCA(3)=.4787E+01_JPRB
GCA(4)=.2102E+04_JPRB
GCA(5)=.1334E+01_JPRB
GCA(6)=.1551E-04_JPRB
GCB(1)=.8968E+07_JPRB
GCB(2)=.2413E+10_JPRB
GCB(3)=.3548E+05_JPRB
GCB(4)=.6370E+10_JPRB
GCB(5)=.8499E+11_JPRB
GCB(6)=.1012E+06_JPRB
GCC(1)=.5925E-10_JPRB
GCC(2)=.1842E-10_JPRB
GCC(3)=.2532E-07_JPRB
GCC(4)=.1953E+07_JPRB
GCC(5)=.1734E-11_JPRB
GCC(6)=.1225E-16_JPRB
GCD4=.3608E-69_JPRB
GCE4=.7563E+04_JPRB
GIREC1=0.0_JPRB
GIREC2=0.0_JPRB
GIREC3=0.0_JPRB
GIREC4=0.0_JPRB
QCO2=.5366E-03_JPRB
QLIMI=15000._JPRB
QLIP0=8.E+06_JPRB
USAA=-.3020_JPRB
USAI=-.3524_JPRB
USAN=-.3400_JPRB
USBA=0._JPRB
USBI=0._JPRB
USBN=0._JPRB
VDP(1,1)=.21868E+02_JPRB
VDP(2,1)=.17453E+02_JPRB
VDP(3,1)=.68918E+00_JPRB
VDP(4,1)=.23456E-03_JPRB
VDP(5,1)=.22317E-09_JPRB
VDP(1,2)=.48401E+02_JPRB
VDP(2,2)=.18648E+02_JPRB
VDP(3,2)=.35199E-01_JPRB
VDP(4,2)=.63691E-06_JPRB
VDP(5,2)=.86395E-13_JPRB
VDP(1,3)=.76948E+02_JPRB
VDP(2,3)=.41056E+01_JPRB
VDP(3,3)=.42667E-03_JPRB
VDP(4,3)=0._JPRB
VDP(5,3)=0._JPRB
VDP(1,4)=.27241E+03_JPRB
VDP(2,4)=.57091E+04_JPRB
VDP(3,4)=.14393E+05_JPRB
VDP(4,4)=.29879E+04_JPRB
VDP(5,4)=.25382E+02_JPRB
VDP(1,5)=.86942E+02_JPRB
VDP(2,5)=.32186E+03_JPRB
VDP(3,5)=.10775E+03_JPRB
VDP(4,5)=.21261E+01_JPRB
VDP(5,5)=.40003E-02_JPRB
VDP(1,6)=.24408E+02_JPRB
VDP(2,6)=.81919E+01_JPRB
VDP(3,6)=.72193E+01_JPRB
VDP(4,6)=.56230E+00_JPRB
VDP(5,6)=.25384E-02_JPRB
VNP(1,1)=.69926E+01_JPRB
VNP(2,1)=.63915E+00_JPRB
VNP(3,1)=.28896E-03_JPRB
VNP(4,1)=.10050E-08_JPRB
VNP(5,1)=.99037E-16_JPRB
VNP(1,2)=.31105E+01_JPRB
VNP(2,2)=.14225E-01_JPRB
VNP(3,2)=.69355E-06_JPRB
VNP(4,2)=.36087E-12_JPRB
VNP(5,2)=.44113E-20_JPRB
VNP(1,3)=.23659E+01_JPRB
VNP(2,3)=.11139E-02_JPRB
VNP(3,3)=.10618E-05_JPRB
VNP(4,3)=0._JPRB
VNP(5,3)=0._JPRB
VNP(1,4)=.19810E+03_JPRB
VNP(2,4)=.46954E+04_JPRB
VNP(3,4)=.22512E+04_JPRB
VNP(4,4)=.52461E+02_JPRB
VNP(5,4)=.11645E+00_JPRB
VNP(1,5)=.46348E+02_JPRB
VNP(2,5)=.35630E+02_JPRB
VNP(3,5)=.33005E+01_JPRB
VNP(4,5)=.18045E-01_JPRB
VNP(5,5)=.88667E-05_JPRB
VNP(1,6)=.47413E+01_JPRB
VNP(2,6)=.16334E+01_JPRB
VNP(3,6)=.48164E+00_JPRB
VNP(4,6)=.56140E-02_JPRB
VNP(5,6)=.67790E-04_JPRB

! Default values for cloud model - first index denotes spectral band:
!   1 - solar
!   2 - thermal
FCM_DEL_A(1)  =  9.55496e-01_JPRB
FCM_DEL_A(2)  =  1.79152e+02_JPRB
FCM_DEL_D(1)  =  5.72959e+01_JPRB
FCM_DEL_D(2)  =  2.00159e+02_JPRB
FCM_MU_A(1)   =  8.66200e-01_JPRB
FCM_MU_A(2)   =  6.65532e-01_JPRB
FCM_MU_D(1)   =  1.00000e+00_JPRB
FCM_MU_D(2)   =  6.70517e-01_JPRB
FCM_N_I       =  1.00000e-01_JPRB
FCM_N_L       =  1.00000e-01_JPRB
FCM_P_AI(1,0) =  1.48920e+01_JPRB
FCM_P_AI(1,1) = -3.71379e+01_JPRB
FCM_P_AI(1,2) =  0.00000e+00_JPRB
FCM_P_AI(1,3) =  0.00000e+00_JPRB
FCM_P_AI(2,0) =  9.76920e+00_JPRB
FCM_P_AI(2,1) = -7.54840e+00_JPRB
FCM_P_AI(2,2) =  0.00000e+00_JPRB
FCM_P_AI(2,3) =  0.00000e+00_JPRB
FCM_P_AL(1,0) =  1.42062e+00_JPRB
FCM_P_AL(1,1) = -3.54037e+00_JPRB
FCM_P_AL(1,2) =  1.73535e+00_JPRB
FCM_P_AL(1,3) =  0.00000e+00_JPRB
FCM_P_AL(2,0) =  5.08522e+00_JPRB
FCM_P_AL(2,1) = -1.08144e+01_JPRB
FCM_P_AL(2,2) =  5.16331e+00_JPRB
FCM_P_AL(2,3) =  0.00000e+00_JPRB
FCM_P_DI(1,0) =  1.30487e+01_JPRB
FCM_P_DI(1,1) =  1.98670e+01_JPRB
FCM_P_DI(1,2) = -5.08578e+01_JPRB
FCM_P_DI(1,3) =  0.00000e+00_JPRB
FCM_P_DI(2,0) =  1.23173e+01_JPRB
FCM_P_DI(2,1) = -3.32808e+01_JPRB
FCM_P_DI(2,2) =  5.77295e+01_JPRB
FCM_P_DI(2,3) = -7.35728e+01_JPRB
FCM_P_DL(1,0) =  4.68646e+00_JPRB
FCM_P_DL(1,1) = -8.77652e+00_JPRB
FCM_P_DL(1,2) =  2.12789e+00_JPRB
FCM_P_DL(1,3) =  0.00000e+00_JPRB
FCM_P_DL(2,0) =  4.52292e+00_JPRB
FCM_P_DL(2,1) = -1.43327e+01_JPRB
FCM_P_DL(2,2) =  1.16650e+01_JPRB
FCM_P_DL(2,3) =  0.00000e+00_JPRB
FCM_P_GI(1,0) = -2.01472e+00_JPRB
FCM_P_GI(1,1) =  9.16284e+01_JPRB
FCM_P_GI(1,2) =  7.82275e+01_JPRB
FCM_P_GI(1,3) =  0.00000e+00_JPRB
FCM_P_GI(2,0) = -7.55124e-01_JPRB
FCM_P_GI(2,1) = -4.68837e+00_JPRB
FCM_P_GI(2,2) =  6.72877e+01_JPRB
FCM_P_GI(2,3) =  0.00000e+00_JPRB
FCM_P_GL(1,0) =  8.82686e-01_JPRB
FCM_P_GL(1,1) = -2.02445e+00_JPRB
FCM_P_GL(1,2) =  1.45259e+00_JPRB
FCM_P_GL(1,3) =  0.00000e+00_JPRB
FCM_P_GL(2,0) =  1.93784e-01_JPRB
FCM_P_GL(2,1) = -6.07676e-01_JPRB
FCM_P_GL(2,2) =  5.37734e-01_JPRB
FCM_P_GL(2,3) =  0.00000e+00_JPRB
FCM_Q_AI(1,1) =  1.86588e+01_JPRB
FCM_Q_AI(1,2) =  0.00000e+00_JPRB
FCM_Q_AI(1,3) =  0.00000e+00_JPRB
FCM_Q_AI(2,1) =  2.78063e+00_JPRB
FCM_Q_AI(2,2) =  0.00000e+00_JPRB
FCM_Q_AI(2,3) =  0.00000e+00_JPRB
FCM_Q_AL(1,1) = -2.81210e+00_JPRB
FCM_Q_AL(1,2) =  2.35758e+00_JPRB
FCM_Q_AL(1,3) =  0.00000e+00_JPRB
FCM_Q_AL(2,1) = -2.31252e+00_JPRB
FCM_Q_AL(2,2) =  1.56848e+00_JPRB
FCM_Q_AL(2,3) =  0.00000e+00_JPRB
FCM_Q_DI(1,1) =  1.04582e+01_JPRB
FCM_Q_DI(1,2) =  6.15641e+00_JPRB
FCM_Q_DI(1,3) =  0.00000e+00_JPRB
FCM_Q_DI(2,1) =  1.18370e+00_JPRB
FCM_Q_DI(2,2) =  2.23187e+00_JPRB
FCM_Q_DI(2,3) =  7.00424e+00_JPRB
FCM_Q_DL(1,1) = -2.18205e+00_JPRB
FCM_Q_DL(1,2) =  1.34831e+00_JPRB
FCM_Q_DL(1,3) =  0.00000e+00_JPRB
FCM_Q_DL(2,1) = -3.12310e+00_JPRB
FCM_Q_DL(2,2) =  2.57959e+00_JPRB
FCM_Q_DL(2,3) =  0.00000e+00_JPRB
FCM_Q_GI(1,1) =  7.89999e+01_JPRB
FCM_Q_GI(1,2) =  5.24929e+01_JPRB
FCM_Q_GI(1,3) =  0.00000e+00_JPRB
FCM_Q_GI(2,1) =  4.93801e+00_JPRB
FCM_Q_GI(2,2) =  2.38096e+01_JPRB
FCM_Q_GI(2,3) =  0.00000e+00_JPRB
FCM_Q_GL(1,1) = -2.20830e+00_JPRB
FCM_Q_GL(1,2) =  1.39832e+00_JPRB
FCM_Q_GL(1,3) =  0.00000e+00_JPRB
FCM_Q_GL(2,1) = -3.27770e+00_JPRB
FCM_Q_GL(2,2) =  2.78446e+00_JPRB
FCM_Q_GL(2,3) =  0.00000e+00_JPRB
REXP_NEB      = 8.0_JPRB

!        1.2 Modify default values according to LECMWF

IF (LECMWF) THEN
ELSE
ENDIF

!     Remark : the values or RII0 is calculated and not set up

!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

! Ce qui concerne NAMPHY3 commente par MPL le 15.04.09
!CALL POSNAM(NULNAM,'NAMPHY3')
!READ(NULNAM,NAMPHY3)
!     ------------------------------------------------------------------

!*       3.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMPHY3 '')')
WRITE(UNIT=KULOUT,FMT='('' BSFSA = '',E10.4,'' BSFSI = '',E10.4 &
 & ,'' BSFSN = '',E10.4,'' BSFTA = '',E10.4,'' BSFTI = '',E10.4,/&
 & ,'' BSFTN = '',E10.4,'' EARRT = '',E10.4,'' EOASA = ''&
 & ,E10.4,'' EOASI = '',E10.4,'' EOASN = '',E10.4 &
 & ,'' EOATA = '',E10.4,/,'' EOATI = '',E10.4,'' EOATN = '',E10.4 &
 & ,'' EODSA = '',E10.4,'' EODSI = '',E10.4,'' EODSN = '',E10.4 &
 & ,'' EODTA = '',E10.4,/,'' EODTI = '',E10.4 &
 & ,'' EODTN = '',E10.4,'' EORAY = '',E10.4,/&
 & ,'' GCA(6) = '',6E11.4,/,'' GCB(6) = '',6E11.4,/,'' GCC(6) = ''&
 & ,6E11.4,/,'' GCD4 = '',E11.4,'' GCE4 = '',E11.4,/&
 & ,'' QCO2 = '',E10.4,'' QLIP0 = '',E10.4 &
 & ,'' USAA = '',E10.4,'' USAI = '',E10.4,'' USAN = '',E10.4 &
 & ,'' USBA = '',E10.4,/,'' USBI = '',E10.4,'' USBN = '',E10.4,/&
 & ,'' VDP(5,6) = '',5E12.5,/,5(12X,5E12.5,/)&
 & ,'' VNP(5,6) = '',5E12.5,/,5(12X,5E12.5,/)&
 & ,'' QLIMI = '',E10.4 &
 & ,'' GIREC1 = '',E10.4,'' GIREC2 = '',E10.4,'' GIREC3 = '',E10.4 &
 & ,'' GIREC4 = '',E10.4 &
 & )')&
 & BSFSA,BSFSI,BSFSN,BSFTA,BSFTI,BSFTN,EARRT,EOASA,EOASI,EOASN,&
 & EOATA,EOATI,EOATN,EODSA,EODSI,EODSN,EODTA,EODTI,EODTN,EORAY,&
 & GCA,GCB,GCC,GCD4,GCE4,QCO2,QLIP0,USAA,USAI,USAN,USBA,USBI,USBN,&
 & VDP,VNP,QLIMI,GIREC1,GIREC2,GIREC3,GIREC4  

IF (LCLSATUR) THEN
  WRITE(KULOUT,'('' N_SPBAND = '',I2,'' REXP_NEB = '',F5.2)') &
   & N_SPBAND,REXP_NEB
  WRITE(KULOUT,'('' FCM_DEL_A(:)   ='',2(1X,ES12.5))') FCM_DEL_A(:)
  WRITE(KULOUT,'('' FCM_DEL_D(:)   ='',2(1X,ES12.5))') FCM_DEL_D(:)
  WRITE(KULOUT,'('' FCM_MU_A (:)   ='',2(1X,ES12.5))') FCM_MU_A (:)
  WRITE(KULOUT,'('' FCM_MU_D (:)   ='',2(1X,ES12.5))') FCM_MU_D (:)
  WRITE(KULOUT,'('' FCM_N_I        ='',1X,ES12.5)')    FCM_N_I
  WRITE(KULOUT,'('' FCM_N_L        ='',1X,ES12.5)')    FCM_N_L
  WRITE(KULOUT,'('' FCM_P_AI (:,0) ='',2(1X,ES12.5))') FCM_P_AI (:,0)
  WRITE(KULOUT,'('' FCM_P_AI (:,1) ='',2(1X,ES12.5))') FCM_P_AI (:,1)
  WRITE(KULOUT,'('' FCM_P_AI (:,2) ='',2(1X,ES12.5))') FCM_P_AI (:,2)
  WRITE(KULOUT,'('' FCM_P_AI (:,3) ='',2(1X,ES12.5))') FCM_P_AI (:,3)
  WRITE(KULOUT,'('' FCM_P_AL (:,0) ='',2(1X,ES12.5))') FCM_P_AL (:,0)
  WRITE(KULOUT,'('' FCM_P_AL (:,1) ='',2(1X,ES12.5))') FCM_P_AL (:,1)
  WRITE(KULOUT,'('' FCM_P_AL (:,2) ='',2(1X,ES12.5))') FCM_P_AL (:,2)
  WRITE(KULOUT,'('' FCM_P_AL (:,3) ='',2(1X,ES12.5))') FCM_P_AL (:,3)
  WRITE(KULOUT,'('' FCM_P_DI (:,0) ='',2(1X,ES12.5))') FCM_P_DI (:,0)
  WRITE(KULOUT,'('' FCM_P_DI (:,1) ='',2(1X,ES12.5))') FCM_P_DI (:,1)
  WRITE(KULOUT,'('' FCM_P_DI (:,2) ='',2(1X,ES12.5))') FCM_P_DI (:,2)
  WRITE(KULOUT,'('' FCM_P_DI (:,3) ='',2(1X,ES12.5))') FCM_P_DI (:,3)
  WRITE(KULOUT,'('' FCM_P_DL (:,0) ='',2(1X,ES12.5))') FCM_P_DL (:,0)
  WRITE(KULOUT,'('' FCM_P_DL (:,1) ='',2(1X,ES12.5))') FCM_P_DL (:,1)
  WRITE(KULOUT,'('' FCM_P_DL (:,2) ='',2(1X,ES12.5))') FCM_P_DL (:,2)
  WRITE(KULOUT,'('' FCM_P_DL (:,3) ='',2(1X,ES12.5))') FCM_P_DL (:,3)
  WRITE(KULOUT,'('' FCM_P_GI (:,0) ='',2(1X,ES12.5))') FCM_P_GI (:,0)
  WRITE(KULOUT,'('' FCM_P_GI (:,1) ='',2(1X,ES12.5))') FCM_P_GI (:,1)
  WRITE(KULOUT,'('' FCM_P_GI (:,2) ='',2(1X,ES12.5))') FCM_P_GI (:,2)
  WRITE(KULOUT,'('' FCM_P_GI (:,3) ='',2(1X,ES12.5))') FCM_P_GI (:,3)
  WRITE(KULOUT,'('' FCM_P_GL (:,0) ='',2(1X,ES12.5))') FCM_P_GL (:,0)
  WRITE(KULOUT,'('' FCM_P_GL (:,1) ='',2(1X,ES12.5))') FCM_P_GL (:,1)
  WRITE(KULOUT,'('' FCM_P_GL (:,2) ='',2(1X,ES12.5))') FCM_P_GL (:,2)
  WRITE(KULOUT,'('' FCM_P_GL (:,3) ='',2(1X,ES12.5))') FCM_P_GL (:,3)
  WRITE(KULOUT,'('' FCM_Q_AI (:,1) ='',2(1X,ES12.5))') FCM_Q_AI (:,1)
  WRITE(KULOUT,'('' FCM_Q_AI (:,2) ='',2(1X,ES12.5))') FCM_Q_AI (:,2)
  WRITE(KULOUT,'('' FCM_Q_AI (:,3) ='',2(1X,ES12.5))') FCM_Q_AI (:,3)
  WRITE(KULOUT,'('' FCM_Q_AL (:,1) ='',2(1X,ES12.5))') FCM_Q_AL (:,1)
  WRITE(KULOUT,'('' FCM_Q_AL (:,2) ='',2(1X,ES12.5))') FCM_Q_AL (:,2)
  WRITE(KULOUT,'('' FCM_Q_AL (:,3) ='',2(1X,ES12.5))') FCM_Q_AL (:,3)
  WRITE(KULOUT,'('' FCM_Q_DI (:,1) ='',2(1X,ES12.5))') FCM_Q_DI (:,1)
  WRITE(KULOUT,'('' FCM_Q_DI (:,2) ='',2(1X,ES12.5))') FCM_Q_DI (:,2)
  WRITE(KULOUT,'('' FCM_Q_DI (:,3) ='',2(1X,ES12.5))') FCM_Q_DI (:,3)
  WRITE(KULOUT,'('' FCM_Q_DL (:,1) ='',2(1X,ES12.5))') FCM_Q_DL (:,1)
  WRITE(KULOUT,'('' FCM_Q_DL (:,2) ='',2(1X,ES12.5))') FCM_Q_DL (:,2)
  WRITE(KULOUT,'('' FCM_Q_DL (:,3) ='',2(1X,ES12.5))') FCM_Q_DL (:,3)
  WRITE(KULOUT,'('' FCM_Q_GI (:,1) ='',2(1X,ES12.5))') FCM_Q_GI (:,1)
  WRITE(KULOUT,'('' FCM_Q_GI (:,2) ='',2(1X,ES12.5))') FCM_Q_GI (:,2)
  WRITE(KULOUT,'('' FCM_Q_GI (:,3) ='',2(1X,ES12.5))') FCM_Q_GI (:,3)
  WRITE(KULOUT,'('' FCM_Q_GL (:,1) ='',2(1X,ES12.5))') FCM_Q_GL (:,1)
  WRITE(KULOUT,'('' FCM_Q_GL (:,2) ='',2(1X,ES12.5))') FCM_Q_GL (:,2)
  WRITE(KULOUT,'('' FCM_Q_GL (:,3) ='',2(1X,ES12.5))') FCM_Q_GL (:,3)
ENDIF
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPHY3',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHY3
