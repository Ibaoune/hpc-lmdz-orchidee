!***************************************************************************
!                                                                          *
!                RRTM :  RAPID RADIATIVE TRANSFER MODEL                    *
!                                                                          *
!             ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                 *
!                        840 MEMORIAL DRIVE                                *
!                        CAMBRIDGE, MA 02139                               *
!                                                                          *
!                           ELI J. MLAWER                                  *
!                         STEVEN J. TAUBMAN~                               *
!                         SHEPARD A. CLOUGH                                *
!                                                                          *
!                        ~currently at GFDL                                *
!                                                                          *
!                       email:  mlawer@aer.com                             *
!                                                                          *
!        The authors wish to acknowledge the contributions of the          *
!        following people:  Patrick D. Brown, Michael J. Iacono,           *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                    *
!                                                                          *
!***************************************************************************
!     Reformatted for F90 by JJMorcrette, ECMWF, 980714                    * 
!                                                                          *
!***************************************************************************
! *** mji ***
! *** This version of RRTM has been altered to interface with either
!     the ECMWF numerical weather prediction model or the ECMWF column 
!     radiation model (ECRT) package. 

!     Revised, April, 1997;  Michael J. Iacono, AER, Inc.
!          - initial implementation of RRTM in ECRT code
!     Revised, June, 1999;  Michael J. Iacono and Eli J. Mlawer, AER, Inc.
!          - to implement generalized maximum/random cloud overlap

SUBROUTINE RRTM_RRTM_140GP &
 & ( KIDIA , KFDIA , KLON , KLEV,&
 & PAER  , PAPH  , PAP,&
 & PTS   , PTH   , PT,&
 & P_ZEMIS , P_ZEMIW,&
 & PQ    , PCCO2 , POZN,&
 & PCLDF , PTAUCLD,&
 & PTAU_LW,&
 & PEMIT , PFLUX , PFLUC, PTCLEAR &
 & )  

! *** This program is the driver for RRTM, the AER rapid model.  
!     For each atmosphere the user wishes to analyze, this routine
!     a) calls ECRTATM to read in the atmospheric profile 
!     b) calls SETCOEF to calculate various quantities needed for 
!        the radiative transfer algorithm
!     c) calls RTRN to do the radiative transfer calculation for
!        clear or cloudy sky
!     d) writes out the upward, downward, and net flux for each
!        level and the heating rate for each layer

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOERAD    ,ONLY : NLW
USE PARRRTM  , ONLY : JPBAND   ,JPXSEC   ,JPGPT    ,JPLAY    ,&
 & JPINPX  
!------------------------------Arguments--------------------------------

! Input arguments

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON! Number of atmospheres (longitudes) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV! Number of atmospheric layers 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA ! First atmosphere index
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA ! Last atmosphere index
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) ! Aerosol optical thickness
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) ! Interface pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) ! Layer pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) ! Surface temperature (I_K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1) ! Interface temperatures (I_K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) ! Layer temperature (I_K)
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
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMIT(KLON) ! Surface LW emissivity
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUX(KLON,2,KLEV+1) ! LW total sky flux (1=up, 2=down)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUC(KLON,2,KLEV+1) ! LW clear sky flux (1=up, 2=down)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTCLEAR(KLON) ! clear-sky fraction of column
INTEGER(KIND=JPIM) :: ICLDLYR(JPLAY)        ! Cloud indicator
REAL(KIND=JPRB) :: Z_CLDFRAC(JPLAY)           ! Cloud fraction
REAL(KIND=JPRB) :: Z_TAUCLD(JPLAY,JPBAND)     ! Spectral optical thickness

REAL(KIND=JPRB) :: Z_ABSS1 (JPGPT*JPLAY)
REAL(KIND=JPRB) :: Z_ATR1  (JPGPT,JPLAY)
EQUIVALENCE (Z_ABSS1(1),Z_ATR1(1,1))

REAL(KIND=JPRB) :: Z_OD    (JPGPT,JPLAY)

REAL(KIND=JPRB) :: Z_TAUSF1(JPGPT*JPLAY)
REAL(KIND=JPRB) :: Z_TF1   (JPGPT,JPLAY)
EQUIVALENCE (Z_TAUSF1(1),Z_TF1(1,1))

REAL(KIND=JPRB) :: Z_COLDRY(JPLAY)
REAL(KIND=JPRB) :: Z_WKL(JPINPX,JPLAY)

REAL(KIND=JPRB) :: Z_WX(JPXSEC,JPLAY)         ! Amount of trace gases

REAL(KIND=JPRB) :: Z_CLFNET  (0:JPLAY)
REAL(KIND=JPRB) :: Z_CLHTR   (0:JPLAY)
REAL(KIND=JPRB) :: Z_FNET    (0:JPLAY)
REAL(KIND=JPRB) :: Z_HTR     (0:JPLAY)
REAL(KIND=JPRB) :: Z_TOTDFLUC(0:JPLAY)
REAL(KIND=JPRB) :: Z_TOTDFLUX(0:JPLAY)
REAL(KIND=JPRB) :: Z_TOTUFLUC(0:JPLAY)
REAL(KIND=JPRB) :: Z_TOTUFLUX(0:JPLAY)

INTEGER(KIND=JPIM) :: i, icld, iplon, I_K
INTEGER(KIND=JPIM) :: ISTART
INTEGER(KIND=JPIM) :: IEND

REAL(KIND=JPRB) :: Z_FLUXFAC, Z_HEATFAC, Z_PI, ZEPSEC, ZTCLEAR

!- from AER
REAL(KIND=JPRB) :: Z_TAUAERL(JPLAY,JPBAND)

!- from INTFAC      
REAL(KIND=JPRB) :: Z_FAC00(JPLAY)
REAL(KIND=JPRB) :: Z_FAC01(JPLAY)
REAL(KIND=JPRB) :: Z_FAC10(JPLAY)
REAL(KIND=JPRB) :: Z_FAC11(JPLAY)
REAL(KIND=JPRB) :: Z_FORFAC(JPLAY)

!- from INTIND
INTEGER(KIND=JPIM) :: JP(JPLAY)
INTEGER(KIND=JPIM) :: JT(JPLAY)
INTEGER(KIND=JPIM) :: JT1(JPLAY)

!- from PRECISE             
REAL(KIND=JPRB) :: Z_ONEMINUS

!- from PROFDATA             
REAL(KIND=JPRB) :: Z_COLH2O(JPLAY)
REAL(KIND=JPRB) :: Z_COLCO2(JPLAY)
REAL(KIND=JPRB) :: Z_COLO3 (JPLAY)
REAL(KIND=JPRB) :: Z_COLN2O(JPLAY)
REAL(KIND=JPRB) :: Z_COLCH4(JPLAY)
REAL(KIND=JPRB) :: Z_COLO2 (JPLAY)
REAL(KIND=JPRB) :: Z_CO2MULT(JPLAY)
INTEGER(KIND=JPIM) :: I_LAYTROP
INTEGER(KIND=JPIM) :: I_LAYSWTCH
INTEGER(KIND=JPIM) :: I_LAYLOW

!- from PROFILE             
REAL(KIND=JPRB) :: Z_PAVEL(JPLAY)
REAL(KIND=JPRB) :: Z_TAVEL(JPLAY)
REAL(KIND=JPRB) :: Z_PZ(0:JPLAY)
REAL(KIND=JPRB) :: Z_TZ(0:JPLAY)
REAL(KIND=JPRB) :: Z_TBOUND
INTEGER(KIND=JPIM) :: I_NLAYERS

!- from SELF             
REAL(KIND=JPRB) :: Z_SELFFAC(JPLAY)
REAL(KIND=JPRB) :: Z_SELFFRAC(JPLAY)
INTEGER(KIND=JPIM) :: INDSELF(JPLAY)

!- from SP             
REAL(KIND=JPRB) :: Z_PFRAC(JPGPT,JPLAY)

!- from SURFACE             
REAL(KIND=JPRB) :: Z_SEMISS(JPBAND)
REAL(KIND=JPRB) :: Z_SEMISLW
INTEGER(KIND=JPIM) :: IREFLECT
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "rrtm_ecrt_140gp.intfb.h"
#include "rrtm_gasabs1a_140gp.intfb.h"
#include "rrtm_rtrn1a_140gp.intfb.h"
#include "rrtm_setcoef_140gp.intfb.h"

!     HEATFAC is the factor by which one must multiply delta-flux/ 
!     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
!     the heating rate in units of degrees/day.  It is equal to 
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(86400)(1e-5)/(1.004)

IF (LHOOK) CALL DR_HOOK('RRTM_RRTM_140GP',0,ZHOOK_HANDLE)
ZEPSEC = 1.E-06_JPRB
Z_ONEMINUS = 1.0_JPRB - ZEPSEC
Z_PI = 2.0_JPRB*ASIN(1.0_JPRB)
Z_FLUXFAC = Z_PI * 2.D4
Z_HEATFAC = 8.4391_JPRB

! *** mji ***
! For use with ECRT, this loop is over atmospheres (or longitudes)
DO iplon = kidia,kfdia

! *** mji ***
!- Prepare atmospheric profile from ECRT for use in RRTM, and define
!  other RRTM input parameters.  Arrays are passed back through the
!  existing RRTM commons and arrays.
  ZTCLEAR=1.0_JPRB

  CALL RRTM_ECRT_140GP &
   & ( iplon, klon , klev, icld,&
   & paer , paph , pap,&
   & pts  , pth  , pt,&
   & P_ZEMIS, P_ZEMIW,&
   & pq   , pcco2, pozn, pcldf, ptaucld, ztclear,&
   & Z_CLDFRAC,Z_TAUCLD,&
   & PTAU_LW,&
   & Z_COLDRY,Z_WKL,Z_WX,&
   & Z_TAUAERL,Z_PAVEL,Z_TAVEL,Z_PZ,Z_TZ,Z_TBOUND,I_NLAYERS,Z_SEMISS,IREFLECT)  

  PTCLEAR(iplon)=ztclear

  ISTART = 1
  IEND   = 16

!  Calculate information needed by the radiative transfer routine
!  that is specific to this atmosphere, especially some of the 
!  coefficients and indices needed to compute the optical depths
!  by interpolating data from stored reference atmospheres. 

  CALL RRTM_SETCOEF_140GP (KLEV,Z_COLDRY,Z_WKL,&
   & Z_FAC00,Z_FAC01,Z_FAC10,Z_FAC11,Z_FORFAC,JP,JT,JT1,&
   & Z_COLH2O,Z_COLCO2,Z_COLO3,Z_COLN2O,Z_COLCH4,Z_COLO2,Z_CO2MULT,&
   & I_LAYTROP,I_LAYSWTCH,I_LAYLOW,Z_PAVEL,Z_TAVEL,Z_SELFFAC,Z_SELFFRAC,INDSELF)  

  CALL RRTM_GASABS1A_140GP (KLEV,Z_ATR1,Z_OD,Z_TF1,Z_COLDRY,Z_WX,&
   & Z_TAUAERL,Z_FAC00,Z_FAC01,Z_FAC10,Z_FAC11,Z_FORFAC,JP,JT,JT1,Z_ONEMINUS,&
   & Z_COLH2O,Z_COLCO2,Z_COLO3,Z_COLN2O,Z_COLCH4,Z_COLO2,Z_CO2MULT,&
   & I_LAYTROP,I_LAYSWTCH,I_LAYLOW,Z_SELFFAC,Z_SELFFRAC,INDSELF,Z_PFRAC)  

!- Call the radiative transfer routine.

! *** mji ***
!  Check for cloud in column.  Use ECRT threshold set as flag icld in
!  routine ECRTATM.  If icld=1 then column is cloudy, otherwise it is
!  clear.  Also, set up flag array, icldlyr, for use in radiative
!  transfer.  Set icldlyr to one for each layer with non-zero cloud
!  fraction.

  DO I_K = 1, KLEV
    IF (ICLD == 1.AND.Z_CLDFRAC(I_K) > ZEPSEC) THEN
      ICLDLYR(I_K) = 1
    ELSE
      ICLDLYR(I_K) = 0
    ENDIF
  ENDDO

!  Clear and cloudy parts of column are treated together in RTRN.
!  Clear radiative transfer is done for clear layers and cloudy radiative
!  transfer is done for cloudy layers as identified by icldlyr.

  CALL RRTM_RTRN1A_140GP (KLEV,ISTART,IEND,ICLDLYR,Z_CLDFRAC,Z_TAUCLD,Z_ABSS1,&
   & Z_OD,Z_TAUSF1,Z_CLFNET,Z_CLHTR,Z_FNET,Z_HTR,Z_TOTDFLUC,Z_TOTDFLUX,Z_TOTUFLUC,Z_TOTUFLUX,&
   & Z_TAVEL,Z_PZ,Z_TZ,Z_TBOUND,Z_PFRAC,Z_SEMISS,Z_SEMISLW,IREFLECT)  

! ***   Pass clear sky and total sky up and down flux profiles to ECRT
!       output arrays (zflux, zfluc). Array indexing from bottom to top 
!       is preserved for ECRT.
!       Invert down flux arrays for consistency with ECRT sign conventions.

  pemit(iplon) = Z_SEMISLW
  DO i = 0, KLEV
    PFLUC(iplon,1,i+1) =  Z_TOTUFLUC(i)*Z_FLUXFAC
    PFLUC(iplon,2,i+1) = -Z_TOTDFLUC(i)*Z_FLUXFAC
    PFLUX(iplon,1,i+1) =  Z_TOTUFLUX(i)*Z_FLUXFAC
    PFLUX(iplon,2,i+1) = -Z_TOTDFLUX(i)*Z_FLUXFAC
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('RRTM_RRTM_140GP',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_RRTM_140GP
