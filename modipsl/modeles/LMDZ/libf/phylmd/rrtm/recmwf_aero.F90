!
! $Id: recmwf_aero.F90 3479 2019-04-17 15:43:34Z fairhead $
!
!OPTIONS XOPT(NOEVAL)
SUBROUTINE RECMWF_AERO (KST, KEND, KPROMA, KTDIA , KLEV,&
 & KMODE,&
 & PALBD , PALBP , PAPRS , PAPRSF , PCCO2 , PCLFR,&
 & PQO3  , PAER  , PDP   , PEMIS  , PMU0,&
 & PQ    , PQS   , PQIWP , PQLWP , PSLM   , PT    , PTS,&
 & PREF_LIQ, PREF_ICE,&
!--OB
 & PREF_LIQ_PI, PREF_ICE_PI,&
!--fin
 & PEMTD , PEMTU , PTRSO,&
 & PTH   , PCTRSO, PCEMTR, PTRSOD,&
 & PLWFC, PLWFT, PSWFC, PSWFT, PSFSWDIR, PSFSWDIF,&
 & PFSDNN, PFSDNV,&  
 & PPIZA_TOT,PCGA_TOT,PTAU_TOT, & 
!--OB
 & PPIZA_NAT,PCGA_NAT,PTAU_NAT, & 
!--fin OB 
!--C.Kleinschmitt
 & PTAU_LW_TOT, PTAU_LW_NAT, & 
!--end
 & PFLUX,PFLUC,&
 & PFSDN ,PFSUP , PFSCDN , PFSCUP, PFSCCDN, PFSCCUP, PFLCCDN, PFLCCUP,&
!--OB diagnostics
 & PTOPSWADAERO,PSOLSWADAERO,&
 & PTOPSWAD0AERO,PSOLSWAD0AERO,&
 & PTOPSWAIAERO,PSOLSWAIAERO,&
 & PTOPSWCFAERO,PSOLSWCFAERO,&
 & PSWADAERO,& !--NL
!--LW diagnostics CK
 & PTOPLWADAERO,PSOLLWADAERO,&
 & PTOPLWAD0AERO,PSOLLWAD0AERO,&
 & PTOPLWAIAERO,PSOLLWAIAERO,&
 & PLWADAERO,& !--NL
!..end
 & ok_ade, ok_aie, ok_volcan, flag_aerosol,flag_aerosol_strat,&
 & flag_aer_feedback)
!--fin

!**** *RECMWF* - METEO-FRANCE RADIATION INTERFACE TO ECMWF RADIATION SCHEME

!     PURPOSE.
!     --------
!           SIMPLE INTERFACE TO RADLSW (NO INTERPOLATION)

!**   INTERFACE.
!     ----------

!     EXPLICIT ARGUMENTS :
!        --------------------
! KST    : START INDEX OF DATA IN KPROMA-LONG VECTOR
! KEND   : END   INDEX OF DATA IN KPROMA-LONG VECTOR
! KPROMA : VECTOR LENGTH
! KTDIA  : INDEX OF TOP LEVEL FROM WHICH COMPUTATIONS ARE ACTIVE
! KLEV   : NUMBER OF LEVELS
! PAER   : (KPROMA,KLEV ,6)     ; OPTICAL THICKNESS OF THE AEROSOLS
! PALBD  : (KPROMA,NSW)         ; DIFFUSE ALBEDO IN THE 2 SW INTERVALS
! PALBP  : (KPROMA,NSW)         ; PARALLEL ALBEDO IN THE 2 SW INTERVALS
! PAPRS  : (KPROMA,KLEV+1)      ; HALF LEVEL PRESSURE
! PAPRSF : (KPROMA,KLEV )       ; FULL LEVEL PRESSURE
! PCCO2  :                      ; CONCENTRATION IN CO2 (PA/PA)
! PCLFR  : (KPROMA,KLEV )       ; CLOUD FRACTIONAL COVER
! PQO3   : (KPROMA,KLEV )       ; OZONE MIXING RATIO (MASS)
! PDP    : (KPROMA,KLEV)        ; LAYER PRESSURE THICKNESS
! PEMIS  : (KPROMA)             ; SURFACE EMISSIVITY
! PMU0   : (KPROMA)             ; SOLAR ANGLE
! PQ     : (KPROMA,KLEV )       ; SPECIFIC HUMIDITY PA/PA
! PQS    : (KPROMA,KLEV )       ; SATURATION SPECIFIC HUMIDITY PA/PA
! PQIWP  : (KPROMA,KLEV )       ; ICE    WATER KG/KG
! PQLWP  : (KPROMA,KLEV )       ; LIQUID WATER KG/KG
! PSLM   : (KPROMA)             ; LAND-SEA MASK
! PT     : (KPROMA,KLEV)        ; FULL LEVEL TEMPERATURE
! PTS    : (KPROMA)             ; SURFACE TEMPERATURE
! PPIZA_TOT  : (KPROMA,KLEV,NSW); Single scattering albedo of total aerosol
! PCGA_TOT   : (KPROMA,KLEV,NSW); Assymetry factor for total aerosol
! PTAU_TOT: (KPROMA,KLEV,NSW)   ; Optical depth of total aerosol
! PREF_LIQ (KPROMA,KLEV)        ; Liquid droplet radius (um) - present-day
! PREF_ICE (KPROMA,KLEV)        ; Ice crystal radius (um) - present-day
!--OB
! PREF_LIQ_PI (KPROMA,KLEV)     ; Liquid droplet radius (um) - pre-industrial 
! PREF_ICE_PI (KPROMA,KLEV)     ; Ice crystal radius (um) - pre-industrial 
! ok_ade---input-L- apply the Aerosol Direct Effect or not?
! ok_aie---input-L- apply the Aerosol Indirect Effect or not?
! ok_volcan-input-L- activate volcanic diags (SW heat & LW cool rate, SW & LW flux)
! flag_aerosol-input-I- aerosol flag from 0 to 7
! flag_aerosol_strat-input-I- use stratospheric aerosols flag (T/F)
! flag_aer_feedback-input-I- use aerosols radiative effect flag (T/F)
! PPIZA_NAT  : (KPROMA,KLEV,NSW); Single scattering albedo of natural aerosol
! PCGA_NAT   : (KPROMA,KLEV,NSW); Assymetry factor for natural aerosol
! PTAU_NAT: (KPROMA,KLEV,NSW)   ; Optical depth of natural aerosol
! PTAU_LW_TOT  (KPROMA,KLEV,NLW); LW Optical depth of total aerosols  
! PTAU_LW_NAT  (KPROMA,KLEV,NLW); LW Optical depth of natural aerosols 
!--fin OB

!     ==== OUTPUTS ===
! PEMTD (KPROMA,KLEV+1)         ; TOTAL DOWNWARD LONGWAVE EMISSIVITY
! PEMTU (KPROMA,KLEV+1)         ; TOTAL UPWARD   LONGWAVE EMISSIVITY
! PTRSO (KPROMA,KLEV+1)         ; TOTAL SHORTWAVE TRANSMISSIVITY
! PTH   (KPROMA,KLEV+1)         ; HALF LEVEL TEMPERATURE
! PCTRSO(KPROMA,2)              ; CLEAR-SKY SHORTWAVE TRANSMISSIVITY
! PCEMTR(KPROMA,2)              ; CLEAR-SKY NET LONGWAVE EMISSIVITY
! PTRSOD(KPROMA)                ; TOTAL-SKY SURFACE SW TRANSMISSITY
! PLWFC (KPROMA,2)              ; CLEAR-SKY LONGWAVE FLUXES
! PLWFT (KPROMA,KLEV+1)         ; TOTAL-SKY LONGWAVE FLUXES
! PSWFC (KPROMA,2)              ; CLEAR-SKY SHORTWAVE FLUXES
! PSWFT (KPROMA,KLEV+1)         ; TOTAL-SKY SHORTWAVE FLUXES
! Ajout flux LW et SW montants et descendants, et ciel clair (MPL 19.12.08)
! PFLUX (KPROMA,2,KLEV+1)       ; LW total sky flux (1=up, 2=down)
! PFLUC (KPROMA,2,KLEV+1)       ; LW clear sky flux (1=up, 2=down)
! PFSDN(KPROMA,KLEV+1)          ; SW total sky flux down
! PFSUP(KPROMA,KLEV+1)          ; SW total sky flux up
! PFSCDN(KPROMA,KLEV+1)         ; SW clear sky flux down
! PFSCUP(KPROMA,KLEV+1)         ; SW clear sky flux up
! PFSCCDN(KPROMA,KLEV+1)        ; SW clear sky clean (no aerosol) flux down
! PFSCCUP(KPROMA,KLEV+1)        ; SW clear sky clean (no aerosol) flux up
! PFLCCDN(KPROMA,KLEV+1)        ; LW clear sky clean (no aerosol) flux down
! PFLCCUP(KPROMA,KLEV+1)        ; LW clear sky clean (no aerosol) flux up


!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------
!     SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHORS.
!     --------
!     ORIGINAL BY  B. RITTER   *ECMWF*        83-10-13
!     REWRITING FOR IFS BY J.-J. MORCRETTE    94-11-15
!     96-11: Ph. Dandin. Meteo-France
!     REWRITING FOR DM  BY J.PH. PIEDELIEVRE   1998-07
!     Duplication of RFMR to use present (cy25) ECMWF radiation scheme : Y. Bouteloup 09-2003
!     Use of 6 aerosols & introduce NSW : F. Bouyssel 09-2004
!     04-11-18 : 4 New arguments for AROME : Y. Seity 
!     2005-10-10 Y. Seity : 3 optional arguments for dust optical properties 
!     JJMorcrette 20060721 PP of clear-sky PAR and TOA incident solar radiation (ECMWF)
!     Olivier Boucher: added LMD radiation diagnostics 2014-03

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOEAERD  , ONLY : RCAEROS
USE YOMCST   , ONLY :         RMD      ,RMO3
USE YOMPHY3  , ONLY : RII0
USE YOERAD   , ONLY : NLW, NAER, RCCNLND  ,RCCNSEA  
USE YOERAD   , ONLY : NAER, RCCNLND  ,RCCNSEA  
USE YOERDU   , ONLY : REPSCQ
USE YOMGEM   , ONLY : NGPTOT
USE YOERDI   , ONLY : RRAE   ,REPCLC    ,REPH2O
USE YOMARPHY , ONLY : LRDUST 
USE phys_output_mod, ONLY : swaerofree_diag, swaero_diag

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

IMPLICIT NONE
INCLUDE "clesphys.h"

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM)               :: KTDIA ! Argument NOT used
INTEGER(KIND=JPIM),INTENT(IN)    :: KMODE
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KPROMA,NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KPROMA,NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KPROMA,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCO2
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLFR(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQO3(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KPROMA,KLEV,6) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDP(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(KPROMA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQIWP(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQLWP(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZA_TOT(KPROMA,KLEV,NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGA_TOT(KPROMA,KLEV,NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU_TOT(KPROMA,KLEV,NSW)
!--OB
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZA_NAT(KPROMA,KLEV,NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGA_NAT(KPROMA,KLEV,NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU_NAT(KPROMA,KLEV,NSW)
REAL(KIND=JPRB)                  :: PPIZA_ZERO(KPROMA,KLEV,NSW)
REAL(KIND=JPRB)                  :: PCGA_ZERO(KPROMA,KLEV,NSW)
REAL(KIND=JPRB)                  :: PTAU_ZERO(KPROMA,KLEV,NSW)
!--fin
!--C.Kleinschmitt
REAL(KIND=JPRB)                  :: PTAU_LW_ZERO(KPROMA,KLEV,NLW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU_LW_TOT(KPROMA,KLEV,NLW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU_LW_NAT(KPROMA,KLEV,NLW)
!--end
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF_LIQ(KPROMA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF_ICE(KPROMA,KLEV)
!--OB
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF_LIQ_PI(KPROMA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF_ICE_PI(KPROMA,KLEV)
LOGICAL, INTENT(in)  :: ok_ade, ok_aie         ! switches whether to use aerosol direct (indirect) effects or not
LOGICAL, INTENT(in)  :: ok_volcan              ! produce volcanic diags (SW/LW heat flux and rate)
INTEGER, INTENT(in)  :: flag_aerosol           ! takes value 0 (no aerosol) or 1 to 6 (aerosols)
LOGICAL, INTENT(in)  :: flag_aerosol_strat     ! use stratospheric aerosols
LOGICAL, INTENT(in)  :: flag_aer_feedback      ! use aerosols radiative feedback
REAL(KIND=JPRB)   ,INTENT(out)   :: PTOPSWADAERO(KPROMA), PSOLSWADAERO(KPROMA)       ! Aerosol direct forcing at TOA and surface
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOPSWAD0AERO(KPROMA), PSOLSWAD0AERO(KPROMA)     ! Aerosol direct forcing at TOA and surface
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOPSWAIAERO(KPROMA), PSOLSWAIAERO(KPROMA)       ! ditto, indirect
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOPSWCFAERO(KPROMA,3), PSOLSWCFAERO(KPROMA,3) !--do we keep this ?
!--fin
!--NL
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSWADAERO(KPROMA, KLEV+1)                        ! SW Aerosol direct forcing
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLWADAERO(KPROMA, KLEV+1)                        ! LW Aerosol direct forcing
!--CK
REAL(KIND=JPRB)   ,INTENT(out)   :: PTOPLWADAERO(KPROMA), PSOLLWADAERO(KPROMA)       ! LW Aerosol direct forcing at TOA + surface
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOPLWAD0AERO(KPROMA), PSOLLWAD0AERO(KPROMA)     ! LW Aerosol direct forcing at TOA + surface
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOPLWAIAERO(KPROMA), PSOLLWAIAERO(KPROMA)       ! LW Aer. indirect forcing at TOA + surface
!--end
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMTD(KPROMA,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMTU(KPROMA,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRSO(KPROMA,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTH(KPROMA,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCTRSO(KPROMA,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCEMTR(KPROMA,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRSOD(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLWFC(KPROMA,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLWFT(KPROMA,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSWFC(KPROMA,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSWFT(KPROMA,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSFSWDIR(KPROMA,NSW)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSFSWDIF(KPROMA,NSW)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSDNN(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSDNV(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUX(KPROMA,2,KLEV+1) ! LW total sky flux (1=up, 2=down)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUC(KPROMA,2,KLEV+1) ! LW clear sky flux (1=up, 2=down)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSDN(KPROMA,KLEV+1)   ! SW total sky flux down
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSUP(KPROMA,KLEV+1)   ! SW total sky flux up
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSCDN(KPROMA,KLEV+1)  ! SW clear sky flux down
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSCUP(KPROMA,KLEV+1)  ! SW clear sky flux up
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSCCDN(KPROMA,KLEV+1) ! SW clear sky clean (no aerosol) flux down
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSCCUP(KPROMA,KLEV+1) ! SW clear sky clean (no aerosol) flux up
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLCCDN(KPROMA,KLEV+1) ! LW clear sky clean (no aerosol) flux down
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLCCUP(KPROMA,KLEV+1) ! LW clear sky clean (no aerosol) flux up

!     ==== COMPUTED IN RADITE ===
!     ------------------------------------------------------------------
!*       0.2   LOCAL ARRAYS.
!              -------------
REAL(KIND=JPRB) :: ZRAER  (KPROMA,6,KLEV)
REAL(KIND=JPRB) :: ZRCLC  (KPROMA,KLEV)
REAL(KIND=JPRB) :: ZRMU0  (KPROMA)
REAL(KIND=JPRB) :: ZRPR   (KPROMA,KLEV)
REAL(KIND=JPRB) :: ZRTI   (KPROMA,KLEV)
REAL(KIND=JPRB) :: ZQLWP (KPROMA,KLEV ) , ZQIWP (KPROMA,KLEV )

REAL(KIND=JPRB) :: ZPQO3 (KPROMA,KLEV)
REAL(KIND=JPRB) :: ZQOZ (NGPTOT,KLEV)
REAL(KIND=JPRB) :: ZQS    (KPROMA,KLEV)
REAL(KIND=JPRB) :: ZQ     (KPROMA,KLEV)
REAL(KIND=JPRB) :: ZEMTD  (KPROMA,KLEV+1)
REAL(KIND=JPRB) :: ZEMTU  (KPROMA,KLEV+1)
REAL(KIND=JPRB) :: ZTRSOC (KPROMA,2)
REAL(KIND=JPRB) :: ZEMTC  (KPROMA,2)

REAL(KIND=JPRB) :: ZNBAS  (KPROMA)
REAL(KIND=JPRB) :: ZNTOP  (KPROMA)
REAL(KIND=JPRB) :: ZQRAIN (KPROMA,KLEV)
REAL(KIND=JPRB) :: ZQRAINT(KPROMA,KLEV)
REAL(KIND=JPRB) :: ZCCNL  (KPROMA)
REAL(KIND=JPRB) :: ZCCNO  (KPROMA)

!  output of radlsw

REAL(KIND=JPRB) :: ZEMIT  (KPROMA)
REAL(KIND=JPRB) :: ZFCT   (KPROMA,KLEV+1)
REAL(KIND=JPRB) :: ZFLT   (KPROMA,KLEV+1)
REAL(KIND=JPRB) :: ZFCS   (KPROMA,KLEV+1)
REAL(KIND=JPRB) :: ZFLS   (KPROMA,KLEV+1)
REAL(KIND=JPRB) :: ZFRSOD (KPROMA),ZSUDU(KPROMA)
REAL(KIND=JPRB) :: ZPARF  (KPROMA),ZUVDF(KPROMA),ZPARCF(KPROMA),ZTINCF(KPROMA)

INTEGER(KIND=JPIM) :: IBEG, IEND, JK, JL

REAL(KIND=JPRB) :: ZCRAE, ZRII0, ZEMIW(KPROMA)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!---aerosol radiative diagnostics
! Key to define the aerosol effect acting on climate
! OB: AEROSOLFEEDBACK_ACTIVE is now a LOGICAL
! TRUE: fluxes use natural and/or anthropogenic aerosols according to ok_ade and ok_aie, DEFAULT
! FALSE: fluxes use no aerosols (case 1) 
! to be used only for maintaining bit reproducibility with aerosol diagnostics activated
 LOGICAL :: AEROSOLFEEDBACK_ACTIVE ! now externalized from .def files

!OB - Fluxes including aerosol effects
!              |        direct effect
!ind effect    | no aerosol  NATural  TOTal
!standard      |   5
!natural (PI)  |               1       3     
!total   (PD)  |               2       4   
! so we need which case when ?
! if flag_aerosol is on
! ok_ade and ok_aie         = 4-2, 4-3 and 4 to proceed
! ok_ade and not ok_aie     = 3-1 and 3 to proceed
! not ok_ade and ok_aie     = 2-1 and 2 to proceed
! not ok_ade and not ok_aie = 1 to proceed
! therefore the cases have the following corresponding switches
! 1 = not ok_ade and not ok_aie OR not ok_ade and ok_aie and swaero_diag OR ok_ade and not ok_aie and swaero_diag
! 2 = not ok_ade and ok_aie OR ok_aie and ok_ade and swaero_diag
! 3 = ok_ade and not ok_aie OR ok_aie and ok_ade and swaero_diag
! 4 = ok_ade and ok_aie
! 5 = no aerosol feedback wanted or no aerosol at all 
! if they are called in this order then the correct call is used to proceed

REAL(KIND=JPRB) ::  ZFSUP_AERO(KPROMA,KLEV+1,5)
REAL(KIND=JPRB) ::  ZFSDN_AERO(KPROMA,KLEV+1,5)
REAL(KIND=JPRB) ::  ZFSUP0_AERO(KPROMA,KLEV+1,5)
REAL(KIND=JPRB) ::  ZFSDN0_AERO(KPROMA,KLEV+1,5)
!--LW (CK):
REAL(KIND=JPRB) ::  LWUP_AERO(KPROMA,KLEV+1,5)
REAL(KIND=JPRB) ::  LWDN_AERO(KPROMA,KLEV+1,5)
REAL(KIND=JPRB) ::  LWUP0_AERO(KPROMA,KLEV+1,5)
REAL(KIND=JPRB) ::  LWDN0_AERO(KPROMA,KLEV+1,5)

#include "radlsw.intfb.h"

IF (LHOOK) CALL DR_HOOK('RECMWF_AERO',0,ZHOOK_HANDLE)
IBEG=KST
IEND=KEND

AEROSOLFEEDBACK_ACTIVE = flag_aer_feedback !NL: externalize aer feedback


!*       1.    PREPARATORY WORK
!              ----------------
!--OB
!        1.0    INITIALIZATIONS
!               --------------

ZFSUP_AERO (:,:,:)=0.
ZFSDN_AERO (:,:,:)=0.
ZFSUP0_AERO(:,:,:)=0.
ZFSDN0_AERO(:,:,:)=0.

LWUP_AERO (:,:,:)=0.
LWDN_AERO (:,:,:)=0.
LWUP0_AERO(:,:,:)=0.
LWDN0_AERO(:,:,:)=0.

PTAU_ZERO(:,:,:) =1.e-15
PPIZA_ZERO(:,:,:)=1.0
PCGA_ZERO(:,:,:) =0.0

PTAU_LW_ZERO(:,:,:) =1.e-15


!*       1.1    LOCAL CONSTANTS
!                ---------------

ZRII0=RII0
ZCRAE=RRAE*(RRAE+2.0_JPRB)

!*       2.1    FULL-LEVEL QUANTITIES

ZRPR =PAPRSF

DO JK=1,KLEV
  DO JL=IBEG,IEND
!   ZPQO3(JL,JK)=PQO3(JL,JK)*PDP(JL,JK)*RMD/RMO3
    ZPQO3(JL,JK)=PQO3(JL,JK)*PDP(JL,JK)
    ZRCLC(JL,JK)=MAX( 0.0_JPRB ,MIN( 1.0_JPRB ,PCLFR(JL,JK)))
    IF (ZRCLC(JL,JK) > REPCLC) THEN
      ZQLWP(JL,JK)=PQLWP(JL,JK)
      ZQIWP(JL,JK)=PQIWP(JL,JK)
    ELSE
      ZQLWP(JL,JK)=REPH2O*ZRCLC(JL,JK)
      ZQIWP(JL,JK)=REPH2O*ZRCLC(JL,JK)
    ENDIF
    ZQRAIN(JL,JK)=0.
    ZQRAINT(JL,JK)=0.
    ZRTI(JL,JK) =PT(JL,JK)
    ZQS (JL,JK)=MAX(2.0_JPRB*REPH2O,PQS(JL,JK))
    ZQ  (JL,JK)=MAX(REPH2O,MIN(PQ(JL,JK),ZQS(JL,JK)*(1.0_JPRB-REPH2O)))
    ZEMIW(JL)=PEMIS(JL)
  ENDDO
ENDDO

IF (NAER == 0) THEN
  ZRAER=RCAEROS
ELSE
  DO JK=1,KLEV
    DO JL=IBEG,IEND
      ZRAER(JL,1,JK)=PAER(JL,JK,1)
      ZRAER(JL,2,JK)=PAER(JL,JK,2)
      ZRAER(JL,3,JK)=PAER(JL,JK,3)
      ZRAER(JL,4,JK)=PAER(JL,JK,4)
      ZRAER(JL,5,JK)=RCAEROS
      ZRAER(JL,6,JK)=PAER(JL,JK,6)
    ENDDO
  ENDDO
ENDIF

!*       2.2    HALF-LEVEL QUANTITIES

DO JK=2,KLEV
  DO JL=IBEG,IEND
    PTH(JL,JK)=&
     & (PT(JL,JK-1)*PAPRSF(JL,JK-1)*(PAPRSF(JL,JK)-PAPRS(JL,JK))&
     & +PT(JL,JK)*PAPRSF(JL,JK)*(PAPRS(JL,JK)-PAPRSF(JL,JK-1)))&
     & *(1.0_JPRB/(PAPRS(JL,JK)*(PAPRSF(JL,JK)-PAPRSF(JL,JK-1))))  
  ENDDO
ENDDO

!*       2.3     QUANTITIES AT BOUNDARIES

DO JL=IBEG,IEND
  PTH(JL,KLEV+1)=PTS(JL)
  PTH(JL,1)=PT(JL,1)-PAPRSF(JL,1)*(PT(JL,1)-PTH(JL,2))&
   & /(PAPRSF(JL,1)-PAPRS(JL,2))  
  ZNBAS(JL)=1. 
  ZNTOP(JL)=1.
  ZCCNL(JL)=RCCNLND
  ZCCNO(JL)=RCCNSEA
ENDDO

!*       3.1     SOLAR ZENITH ANGLE IS EARTH'S CURVATURE
!                CORRECTED

! CCMVAL: on impose ZRMU0=PMU0 MPL 25032010
! 2eme essai en 3D MPL 20052010
!DO JL=IBEG,IEND
! ZRMU0(JL)=PMU0(JL)
!ENDDO
!!!!! A REVOIR MPL 20091201: enleve cette correction pour comparer a AR4
 DO JL=IBEG,IEND
   IF (PMU0(JL) > 1.E-10_JPRB) THEN
     ZRMU0(JL)=RRAE/(SQRT(PMU0(JL)**2+ZCRAE)-PMU0(JL))
   ELSE
     ZRMU0(JL)= RRAE/SQRT(ZCRAE)
   ENDIF   
 ENDDO   

!*         4.1     CALL TO ACTUAL RADIATION SCHEME
!
!----now we make multiple calls to the radiation according to which 
!----aerosol flags are on

IF (flag_aerosol .GT. 0 .OR. flag_aerosol_strat) THEN

!--Case 1
IF ( ( .not. ok_ade .AND. .not. ok_aie ) .OR.             & 
   & ( .not. ok_ade .AND. ok_aie .AND. swaero_diag ) .OR. & 
   & ( ok_ade .AND. .not. ok_aie .AND. swaero_diag ) ) THEN

! natural aerosols for direct and indirect effect 
! PI cloud optical properties
! use PREF_LIQ_PI and PREF_ICE_PI
! use NAT aerosol optical properties
! store fluxes in index 1

CALL RADLSW (&
 & IBEG  , IEND   , KPROMA  , KLEV  , KMODE , NAER,&
 & ZRII0 ,&
 & ZRAER , PALBD  , PALBP   , PAPRS , ZRPR  ,&
 & ZCCNL , ZCCNO  ,&
 & PCCO2 , ZRCLC  , PDP     , PEMIS , ZEMIW ,PSLM    , ZRMU0 , ZPQO3,&
 & ZQ    , ZQIWP  , ZQLWP   , ZQS   , ZQRAIN,ZQRAINT ,&
 & PTH   , ZRTI   , PTS     , ZNBAS , ZNTOP ,&
 & PREF_LIQ_PI, PREF_ICE_PI,&
 & ZEMIT , ZFCT   , ZFLT    , ZFCS    , ZFLS  ,&
 & ZFRSOD, ZSUDU  , ZUVDF   , ZPARF   , ZPARCF, ZTINCF, PSFSWDIR,&
 & PSFSWDIF,PFSDNN, PFSDNV  ,&  
 & LRDUST,PPIZA_NAT,PCGA_NAT,PTAU_NAT,PTAU_LW_NAT,PFLUX,PFLUC,&
 & PFSDN , PFSUP  , PFSCDN  , PFSCUP )

!* SAVE VARIABLES IN INTERIM VARIABLES A LA SW_AEROAR4
ZFSUP0_AERO(:,:,1) = PFSCUP(:,:)
ZFSDN0_AERO(:,:,1) = PFSCDN(:,:)

ZFSUP_AERO(:,:,1) =  PFSUP(:,:)
ZFSDN_AERO(:,:,1) =  PFSDN(:,:)

LWUP0_AERO(:,:,1) = PFLUC(:,1,:)
LWDN0_AERO(:,:,1) = PFLUC(:,2,:)

LWUP_AERO(:,:,1) = PFLUX(:,1,:)
LWDN_AERO(:,:,1) = PFLUX(:,2,:)

ENDIF 

!--Case 2
IF ( ( .not. ok_ade .AND. ok_aie ) .OR. & 
   & ( ok_ade .AND. ok_aie .AND. swaero_diag ) ) THEN

! natural aerosols for direct indirect effect 
! use NAT aerosol optical properties 
! PD cloud optical properties
! use PREF_LIQ and PREF_ICE
! store fluxes in index 2

CALL RADLSW (&
 & IBEG  , IEND   , KPROMA  , KLEV  , KMODE , NAER,&
 & ZRII0 ,&
 & ZRAER , PALBD  , PALBP   , PAPRS , ZRPR  ,&
 & ZCCNL , ZCCNO  ,&
 & PCCO2 , ZRCLC  , PDP     , PEMIS , ZEMIW ,PSLM    , ZRMU0 , ZPQO3,&
 & ZQ    , ZQIWP  , ZQLWP   , ZQS   , ZQRAIN,ZQRAINT ,&
 & PTH   , ZRTI   , PTS     , ZNBAS , ZNTOP ,&
 & PREF_LIQ, PREF_ICE,&
 & ZEMIT , ZFCT   , ZFLT    , ZFCS    , ZFLS  ,&
 & ZFRSOD, ZSUDU  , ZUVDF   , ZPARF   , ZPARCF, ZTINCF, PSFSWDIR,&
 & PSFSWDIF,PFSDNN, PFSDNV  ,&  
 & LRDUST,PPIZA_NAT,PCGA_NAT,PTAU_NAT,PTAU_LW_NAT,PFLUX,PFLUC,&
 & PFSDN , PFSUP  , PFSCDN  , PFSCUP )

!* SAVE VARIABLES IN INTERIM VARIABLES A LA SW_AEROAR4
ZFSUP0_AERO(:,:,2) = PFSCUP(:,:)
ZFSDN0_AERO(:,:,2) = PFSCDN(:,:)

ZFSUP_AERO(:,:,2) =  PFSUP(:,:)
ZFSDN_AERO(:,:,2) =  PFSDN(:,:)

LWUP0_AERO(:,:,2) = PFLUC(:,1,:)
LWDN0_AERO(:,:,2) = PFLUC(:,2,:)

LWUP_AERO(:,:,2) = PFLUX(:,1,:)
LWDN_AERO(:,:,2) = PFLUX(:,2,:)

ENDIF ! ok_aie      

!--Case 3
IF ( ( ok_ade .AND. .not. ok_aie ) .OR. &
   & ( ok_ade .AND. ok_aie .AND. swaero_diag ) ) THEN

! direct effect of total aerosol activated
! TOT aerosols for direct effect
! PI cloud optical properties
! use PREF_LIQ_PI and PREF_ICE_PI
! STORE fluxes in index 3
 
CALL RADLSW (&
 & IBEG  , IEND   , KPROMA  , KLEV  , KMODE , NAER,&
 & ZRII0 ,&
 & ZRAER , PALBD  , PALBP   , PAPRS , ZRPR  ,&
 & ZCCNL , ZCCNO  ,&
 & PCCO2 , ZRCLC  , PDP     , PEMIS , ZEMIW ,PSLM    , ZRMU0 , ZPQO3,&
 & ZQ    , ZQIWP  , ZQLWP   , ZQS   , ZQRAIN,ZQRAINT ,&
 & PTH   , ZRTI   , PTS     , ZNBAS , ZNTOP ,&
 & PREF_LIQ_PI, PREF_ICE_PI,&
 & ZEMIT , ZFCT   , ZFLT    , ZFCS    , ZFLS  ,&
 & ZFRSOD, ZSUDU  , ZUVDF   , ZPARF   , ZPARCF, ZTINCF, PSFSWDIR,&
 & PSFSWDIF,PFSDNN, PFSDNV  ,&  
 & LRDUST,PPIZA_TOT,PCGA_TOT,PTAU_TOT,PTAU_LW_TOT,PFLUX,PFLUC,&
 & PFSDN , PFSUP  , PFSCDN  , PFSCUP )

!* SAVE VARIABLES IN INTERIM VARIABLES A LA SW_AEROAR4
ZFSUP0_AERO(:,:,3) = PFSCUP(:,:)
ZFSDN0_AERO(:,:,3) = PFSCDN(:,:)

ZFSUP_AERO(:,:,3) =  PFSUP(:,:)
ZFSDN_AERO(:,:,3) =  PFSDN(:,:)

LWUP0_AERO(:,:,3) = PFLUC(:,1,:)
LWDN0_AERO(:,:,3) = PFLUC(:,2,:)

LWUP_AERO(:,:,3) = PFLUX(:,1,:)
LWDN_AERO(:,:,3) = PFLUX(:,2,:)

ENDIF !-end ok_ade

!--Case 4
IF (ok_ade .and. ok_aie) THEN

! total aerosols for direct indirect effect 
! use TOT aerosol optical properties 
! PD cloud optical properties
! use PREF_LIQ and PREF_ICE
! store fluxes in index 4 

CALL RADLSW (&
 & IBEG  , IEND   , KPROMA  , KLEV  , KMODE , NAER,&
 & ZRII0 ,&
 & ZRAER , PALBD  , PALBP   , PAPRS , ZRPR  ,&
 & ZCCNL , ZCCNO  ,&
 & PCCO2 , ZRCLC  , PDP     , PEMIS , ZEMIW ,PSLM    , ZRMU0 , ZPQO3,&
 & ZQ    , ZQIWP  , ZQLWP   , ZQS   , ZQRAIN,ZQRAINT ,&
 & PTH   , ZRTI   , PTS     , ZNBAS , ZNTOP ,&
 & PREF_LIQ, PREF_ICE,&
 & ZEMIT , ZFCT   , ZFLT    , ZFCS    , ZFLS  ,&
 & ZFRSOD, ZSUDU  , ZUVDF   , ZPARF   , ZPARCF, ZTINCF, PSFSWDIR,&
 & PSFSWDIF,PFSDNN, PFSDNV  ,&  
 & LRDUST,PPIZA_TOT,PCGA_TOT,PTAU_TOT,PTAU_LW_TOT,PFLUX,PFLUC,&
 & PFSDN , PFSUP  , PFSCDN  , PFSCUP )

!* SAVE VARIABLES IN INTERIM VARIABLES A LA SW_AEROAR4
ZFSUP0_AERO(:,:,4) = PFSCUP(:,:)
ZFSDN0_AERO(:,:,4) = PFSCDN(:,:)

ZFSUP_AERO(:,:,4) =  PFSUP(:,:)
ZFSDN_AERO(:,:,4) =  PFSDN(:,:)

LWUP0_AERO(:,:,4) = PFLUC(:,1,:)
LWDN0_AERO(:,:,4) = PFLUC(:,2,:)

LWUP_AERO(:,:,4) = PFLUX(:,1,:)
LWDN_AERO(:,:,4) = PFLUX(:,2,:)

ENDIF ! ok_ade .and. ok_aie

ENDIF !--if flag_aerosol GT 0 OR flag_aerosol_strat

! case with no aerosols at all is also computed IF ACTIVEFEEDBACK_ACTIVE is false 
IF (.not. AEROSOLFEEDBACK_ACTIVE .OR. flag_aerosol .EQ. 0 .OR. swaerofree_diag) THEN    

! ZERO aerosol effect
! ZERO aerosol optical depth
! STANDARD cloud optical properties
! STORE fluxes in index 5

CALL RADLSW (&
 & IBEG  , IEND   , KPROMA  , KLEV  , KMODE , NAER,&
 & ZRII0 ,&
 & ZRAER , PALBD  , PALBP   , PAPRS , ZRPR  ,&
 & ZCCNL , ZCCNO  ,&
 & PCCO2 , ZRCLC  , PDP     , PEMIS , ZEMIW ,PSLM    , ZRMU0 , ZPQO3,&
 & ZQ    , ZQIWP  , ZQLWP   , ZQS   , ZQRAIN,ZQRAINT ,&
 & PTH   , ZRTI   , PTS     , ZNBAS , ZNTOP ,&
!--this needs to be changed to fixed cloud optical properties
 & PREF_LIQ_PI, PREF_ICE_PI,&
 & ZEMIT , ZFCT   , ZFLT    , ZFCS    , ZFLS  ,&
 & ZFRSOD, ZSUDU  , ZUVDF   , ZPARF   , ZPARCF, ZTINCF, PSFSWDIR,&
 & PSFSWDIF,PFSDNN, PFSDNV  ,&  
 & LRDUST,PPIZA_ZERO,PCGA_ZERO,PTAU_ZERO, PTAU_LW_ZERO,PFLUX,PFLUC,&
 & PFSDN , PFSUP  , PFSCDN  , PFSCUP )

!* SAVE VARIABLES IN INTERIM VARIABLES A LA SW_AEROAR4
ZFSUP0_AERO(:,:,5) = PFSCUP(:,:)
ZFSDN0_AERO(:,:,5) = PFSCDN(:,:)

ZFSUP_AERO(:,:,5) =  PFSUP(:,:)
ZFSDN_AERO(:,:,5) =  PFSDN(:,:)

LWUP0_AERO(:,:,5) = PFLUC(:,1,:)
LWDN0_AERO(:,:,5) = PFLUC(:,2,:)

LWUP_AERO(:,:,5) = PFLUX(:,1,:)
LWDN_AERO(:,:,5) = PFLUX(:,2,:)

ENDIF ! .not. AEROSOLFEEDBACK_ACTIVE

!*         4.2     TRANSFORM FLUXES TO MODEL HISTORICAL VARIABLES

DO JK=1,KLEV+1
  DO JL=IBEG,IEND
    PSWFT(JL,JK)=ZFLS(JL,JK)/(ZRII0*ZRMU0(JL))
    PLWFT(JL,JK)=ZFLT(JL,JK)
  ENDDO
ENDDO

ZEMTD=PLWFT
ZEMTU=PLWFT

DO JL=IBEG,IEND
  ZTRSOC(JL, 1)=ZFCS(JL,     1)/(ZRII0*ZRMU0(JL))
  ZTRSOC(JL, 2)=ZFCS(JL,KLEV+1)/(ZRII0*ZRMU0(JL))
  ZEMTC (JL, 1)=ZFCT(JL,     1)
  ZEMTC (JL, 2)=ZFCT(JL,KLEV+1)
ENDDO

!                 ------------ -- ------- -- ---- -----
!*         5.1    STORAGE OF TRANSMISSIVITY AND EMISSIVITIES
!*                IN KPROMA-LONG ARRAYS

DO JK=1,KLEV+1
  DO JL=IBEG,IEND
    PEMTD(JL,JK)=ZEMTD(JL,JK)
    PEMTU(JL,JK)=ZEMTU(JL,JK)
    PTRSO(JL,JK)=MAX(0.0_JPRB,MIN(1.0_JPRB,PSWFT(JL,JK)))
  ENDDO
ENDDO
DO JK=1,2
  DO JL=IBEG,IEND
    PCEMTR(JL,JK)=ZEMTC (JL,JK)
    PCTRSO(JL,JK)=MAX( 0.0_JPRB,MIN(1.0_JPRB,ZTRSOC(JL,JK)))
  ENDDO
ENDDO
DO JL=IBEG,IEND
  PTRSOD(JL)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZFRSOD(JL)/(ZRII0*ZRMU0(JL))))
ENDDO

!*         7.3   RECONSTRUCT FLUXES FOR DIAGNOSTICS

DO JL=IBEG,IEND
  IF (PMU0(JL) < 1.E-10_JPRB) ZRMU0(JL)=0.0_JPRB
ENDDO
DO JK=1,KLEV+1
  DO JL=IBEG,IEND
    PLWFT(JL,JK)=PEMTD(JL,JK)
    PSWFT(JL,JK)=ZRMU0(JL)*ZRII0*PTRSO(JL,JK)
  ENDDO
ENDDO
DO JK=1,2
  DO JL=IBEG,IEND
    PSWFC(JL,JK)=ZRMU0(JL)*ZRII0*PCTRSO(JL,JK)
    PLWFC(JL,JK)=PCEMTR(JL,JK)
  ENDDO
ENDDO

!*  8.0 DIAGNOSTICS
!---Now we copy back the correct fields to proceed to the next timestep

IF  ( AEROSOLFEEDBACK_ACTIVE .AND. (flag_aerosol .GT. 0 .OR. flag_aerosol_strat) ) THEN

  IF ( ok_ade .and. ok_aie  ) THEN
    PFSUP(:,:) =    ZFSUP_AERO(:,:,4)
    PFSDN(:,:) =    ZFSDN_AERO(:,:,4)
    PFSCUP(:,:) =   ZFSUP0_AERO(:,:,4)
    PFSCDN(:,:) =   ZFSDN0_AERO(:,:,4)

    PFLUX(:,1,:) =  LWUP_AERO(:,:,4)
    PFLUX(:,2,:) =  LWDN_AERO(:,:,4)
    PFLUC(:,1,:) =  LWUP0_AERO(:,:,4)
    PFLUC(:,2,:) =  LWDN0_AERO(:,:,4)    
  ENDIF

  IF ( ok_ade .and. (.not. ok_aie) )  THEN
    PFSUP(:,:) =    ZFSUP_AERO(:,:,3)
    PFSDN(:,:) =    ZFSDN_AERO(:,:,3)
    PFSCUP(:,:) =   ZFSUP0_AERO(:,:,3)
    PFSCDN(:,:) =   ZFSDN0_AERO(:,:,3)

    PFLUX(:,1,:) =  LWUP_AERO(:,:,3)
    PFLUX(:,2,:) =  LWDN_AERO(:,:,3)
    PFLUC(:,1,:) =  LWUP0_AERO(:,:,3)
    PFLUC(:,2,:) =  LWDN0_AERO(:,:,3) 
  ENDIF

  IF ( (.not. ok_ade) .and. ok_aie  )  THEN
    PFSUP(:,:) =    ZFSUP_AERO(:,:,2)
    PFSDN(:,:) =    ZFSDN_AERO(:,:,2)
    PFSCUP(:,:) =   ZFSUP0_AERO(:,:,2)
    PFSCDN(:,:) =   ZFSDN0_AERO(:,:,2)

    PFLUX(:,1,:) =  LWUP_AERO(:,:,2)
    PFLUX(:,2,:) =  LWDN_AERO(:,:,2)
    PFLUC(:,1,:) =  LWUP0_AERO(:,:,2)
    PFLUC(:,2,:) =  LWDN0_AERO(:,:,2) 
  ENDiF

  IF ((.not. ok_ade) .and. (.not. ok_aie)) THEN
    PFSUP(:,:) =    ZFSUP_AERO(:,:,1)
    PFSDN(:,:) =    ZFSDN_AERO(:,:,1)
    PFSCUP(:,:) =   ZFSUP0_AERO(:,:,1)
    PFSCDN(:,:) =   ZFSDN0_AERO(:,:,1)

    PFLUX(:,1,:) =  LWUP_AERO(:,:,1)
    PFLUX(:,2,:) =  LWDN_AERO(:,:,1)
    PFLUC(:,1,:) =  LWUP0_AERO(:,:,1)
    PFLUC(:,2,:) =  LWDN0_AERO(:,:,1) 
  ENDIF

! The following allows to compute the forcing diagostics without
! letting the aerosol forcing act on the meteorology
! SEE logic above

ELSE  !--not AEROSOLFEEDBACK_ACTIVE

    PFSUP(:,:) =    ZFSUP_AERO(:,:,5)
    PFSDN(:,:) =    ZFSDN_AERO(:,:,5)
    PFSCUP(:,:) =   ZFSUP0_AERO(:,:,5)
    PFSCDN(:,:) =   ZFSDN0_AERO(:,:,5)

    PFLUX(:,1,:) =  LWUP_AERO(:,:,5)
    PFLUX(:,2,:) =  LWDN_AERO(:,:,5)
    PFLUC(:,1,:) =  LWUP0_AERO(:,:,5)
    PFLUC(:,2,:) =  LWDN0_AERO(:,:,5) 

ENDIF

IF (swaerofree_diag) THEN
! copy shortwave clear-sky clean (no aerosol) case
  PFSCCUP(:,:) =   ZFSUP0_AERO(:,:,5)
  PFSCCDN(:,:) =   ZFSDN0_AERO(:,:,5)
! copy longwave clear-sky clean (no aerosol) case
  PFLCCUP(:,:) =   LWUP0_AERO(:,:,5)
  PFLCCDN(:,:) =   LWDN0_AERO(:,:,5)
ENDIF

!OB- HERE CHECK WITH MP IF BOTTOM AND TOP INDICES ARE OK !!!!!!!!!!!!!!!!!!
! net anthropogenic forcing direct and 1st indirect effect diagnostics
! requires a natural aerosol field read and used 
! Difference of net fluxes from double call to radiation
! Will need to be extended to LW radiation -> done by CK (2014-05-23)

IF (flag_aerosol .GT. 0 .OR. flag_aerosol_strat) THEN

IF (ok_ade.AND.ok_aie) THEN

! direct anthropogenic forcing
     PSOLSWADAERO(:)  = (ZFSDN_AERO(:,1,4)      -ZFSUP_AERO(:,1,4))      -(ZFSDN_AERO(:,1,2)      -ZFSUP_AERO(:,1,2))
     PTOPSWADAERO(:)  = (ZFSDN_AERO(:,KLEV+1,4) -ZFSUP_AERO(:,KLEV+1,4)) -(ZFSDN_AERO(:,KLEV+1,2) -ZFSUP_AERO(:,KLEV+1,2))
     PSOLSWAD0AERO(:) = (ZFSDN0_AERO(:,1,4)     -ZFSUP0_AERO(:,1,4))     -(ZFSDN0_AERO(:,1,2)     -ZFSUP0_AERO(:,1,2))
     PTOPSWAD0AERO(:) = (ZFSDN0_AERO(:,KLEV+1,4)-ZFSUP0_AERO(:,KLEV+1,4))-(ZFSDN0_AERO(:,KLEV+1,2)-ZFSUP0_AERO(:,KLEV+1,2))
     IF(ok_volcan) THEN
        PSWADAERO(:,:)  = (ZFSDN_AERO(:,:,4) -ZFSUP_AERO(:,:,4)) -(ZFSDN_AERO(:,:,2) -ZFSUP_AERO(:,:,2)) !--NL
     ENDIF

! indirect anthropogenic forcing
     PSOLSWAIAERO(:) = (ZFSDN_AERO(:,1,4)     -ZFSUP_AERO(:,1,4))     -(ZFSDN_AERO(:,1,3)     -ZFSUP_AERO(:,1,3))
     PTOPSWAIAERO(:) = (ZFSDN_AERO(:,KLEV+1,4)-ZFSUP_AERO(:,KLEV+1,4))-(ZFSDN_AERO(:,KLEV+1,3)-ZFSUP_AERO(:,KLEV+1,3))

! Cloud radiative forcing with natural aerosol for direct effect 
     PSOLSWCFAERO(:,1) = (ZFSDN_AERO(:,1,2)     -ZFSUP_AERO(:,1,2))     -(ZFSDN0_AERO(:,1,2)     -ZFSUP0_AERO(:,1,2))
     PTOPSWCFAERO(:,1) = (ZFSDN_AERO(:,KLEV+1,2)-ZFSUP_AERO(:,KLEV+1,2))-(ZFSDN0_AERO(:,KLEV+1,2)-ZFSUP0_AERO(:,KLEV+1,2))
! Cloud radiative forcing with anthropogenic aerosol for direct effect 
     PSOLSWCFAERO(:,2) = (ZFSDN_AERO(:,1,4)     -ZFSUP_AERO(:,1,4))     -(ZFSDN0_AERO(:,1,4)     -ZFSUP0_AERO(:,1,4))
     PTOPSWCFAERO(:,2) = (ZFSDN_AERO(:,KLEV+1,4)-ZFSUP_AERO(:,KLEV+1,4))-(ZFSDN0_AERO(:,KLEV+1,4)-ZFSUP0_AERO(:,KLEV+1,4))
! Cloud radiative forcing with no direct effect at all 
     PSOLSWCFAERO(:,3) = 0.0
     PTOPSWCFAERO(:,3) = 0.0

! LW direct anthropogenic forcing
     PSOLLWADAERO(:)  = (-LWDN_AERO(:,1,4)      -LWUP_AERO(:,1,4))      -(-LWDN_AERO(:,1,2)      -LWUP_AERO(:,1,2))
     PTOPLWADAERO(:)  = (-LWDN_AERO(:,KLEV+1,4) -LWUP_AERO(:,KLEV+1,4)) -(-LWDN_AERO(:,KLEV+1,2) -LWUP_AERO(:,KLEV+1,2))
     PSOLLWAD0AERO(:) = (-LWDN0_AERO(:,1,4)     -LWUP0_AERO(:,1,4))     -(-LWDN0_AERO(:,1,2)     -LWUP0_AERO(:,1,2))
     PTOPLWAD0AERO(:) = (-LWDN0_AERO(:,KLEV+1,4)-LWUP0_AERO(:,KLEV+1,4))-(-LWDN0_AERO(:,KLEV+1,2)-LWUP0_AERO(:,KLEV+1,2))
     IF(ok_volcan) THEN
        PLWADAERO(:,:)  = (-LWDN_AERO(:,:,4) -LWUP_AERO(:,:,4)) -(-LWDN_AERO(:,:,2) -LWUP_AERO(:,:,2)) !--NL
     ENDIF

! LW indirect anthropogenic forcing
     PSOLLWAIAERO(:) = (-LWDN_AERO(:,1,4)     -LWUP_AERO(:,1,4))     -(-LWDN_AERO(:,1,3)     -LWUP_AERO(:,1,3))
     PTOPLWAIAERO(:) = (-LWDN_AERO(:,KLEV+1,4)-LWUP_AERO(:,KLEV+1,4))-(-LWDN_AERO(:,KLEV+1,3)-LWUP_AERO(:,KLEV+1,3))

ENDIF

IF (ok_ade.AND..NOT.ok_aie) THEN

! direct anthropogenic forcing
     PSOLSWADAERO(:)  = (ZFSDN_AERO(:,1,3)      -ZFSUP_AERO(:,1,3))      -(ZFSDN_AERO(:,1,1)      -ZFSUP_AERO(:,1,1))
     PTOPSWADAERO(:)  = (ZFSDN_AERO(:,KLEV+1,3) -ZFSUP_AERO(:,KLEV+1,3)) -(ZFSDN_AERO(:,KLEV+1,1) -ZFSUP_AERO(:,KLEV+1,1))
     PSOLSWAD0AERO(:) = (ZFSDN0_AERO(:,1,3)     -ZFSUP0_AERO(:,1,3))     -(ZFSDN0_AERO(:,1,1)     -ZFSUP0_AERO(:,1,1))
     PTOPSWAD0AERO(:) = (ZFSDN0_AERO(:,KLEV+1,3)-ZFSUP0_AERO(:,KLEV+1,3))-(ZFSDN0_AERO(:,KLEV+1,1)-ZFSUP0_AERO(:,KLEV+1,1))
     IF(ok_volcan) THEN
        PSWADAERO(:,:)  = (ZFSDN_AERO(:,:,3) -ZFSUP_AERO(:,:,3)) -(ZFSDN_AERO(:,:,1) -ZFSUP_AERO(:,:,1)) !--NL
     ENDIF

! indirect anthropogenic forcing
     PSOLSWAIAERO(:) = 0.0
     PTOPSWAIAERO(:) = 0.0 

! Cloud radiative forcing with natural aerosol for direct effect 
     PSOLSWCFAERO(:,1) = (ZFSDN_AERO(:,1,1)     -ZFSUP_AERO(:,1,1))     -(ZFSDN0_AERO(:,1,1)     -ZFSUP0_AERO(:,1,1))
     PTOPSWCFAERO(:,1) = (ZFSDN_AERO(:,KLEV+1,1)-ZFSUP_AERO(:,KLEV+1,1))-(ZFSDN0_AERO(:,KLEV+1,1)-ZFSUP0_AERO(:,KLEV+1,1))
! Cloud radiative forcing with anthropogenic aerosol for direct effect 
     PSOLSWCFAERO(:,2) = (ZFSDN_AERO(:,1,3)     -ZFSUP_AERO(:,1,3))     -(ZFSDN0_AERO(:,1,3)     -ZFSUP0_AERO(:,1,3))
     PTOPSWCFAERO(:,2) = (ZFSDN_AERO(:,KLEV+1,3)-ZFSUP_AERO(:,KLEV+1,3))-(ZFSDN0_AERO(:,KLEV+1,3)-ZFSUP0_AERO(:,KLEV+1,3))
! Cloud radiative forcing with no direct effect at all 
     PSOLSWCFAERO(:,3) = 0.0
     PTOPSWCFAERO(:,3) = 0.0

! LW direct anthropogenic forcing
     PSOLLWADAERO(:)  = (-LWDN_AERO(:,1,3)      -LWUP_AERO(:,1,3))      -(-LWDN_AERO(:,1,1)      -LWUP_AERO(:,1,1))
     PTOPLWADAERO(:)  = (-LWDN_AERO(:,KLEV+1,3) -LWUP_AERO(:,KLEV+1,3)) -(-LWDN_AERO(:,KLEV+1,1) -LWUP_AERO(:,KLEV+1,1))
     PSOLLWAD0AERO(:) = (-LWDN0_AERO(:,1,3)     -LWUP0_AERO(:,1,3))     -(-LWDN0_AERO(:,1,1)     -LWUP0_AERO(:,1,1))
     PTOPLWAD0AERO(:) = (-LWDN0_AERO(:,KLEV+1,3)-LWUP0_AERO(:,KLEV+1,3))-(-LWDN0_AERO(:,KLEV+1,1)-LWUP0_AERO(:,KLEV+1,1))
     IF(ok_volcan) THEN
        PLWADAERO(:,:)  = (-LWDN_AERO(:,:,3) -LWUP_AERO(:,:,3)) -(-LWDN_AERO(:,:,1) -LWUP_AERO(:,:,1)) !--NL
     ENDIF
     
! LW indirect anthropogenic forcing
     PSOLLWAIAERO(:) = 0.0
     PTOPLWAIAERO(:) = 0.0 

ENDIF

IF (.NOT.ok_ade.AND.ok_aie) THEN

! direct anthropogenic forcing
     PSOLSWADAERO(:)  = 0.0
     PTOPSWADAERO(:)  = 0.0
     PSOLSWAD0AERO(:) = 0.0
     PTOPSWAD0AERO(:) = 0.0
     IF(ok_volcan) THEN
        PSWADAERO(:,:)  = 0.0 !--NL
     ENDIF
     
! indirect anthropogenic forcing
     PSOLSWAIAERO(:) = (ZFSDN_AERO(:,1,2)     -ZFSUP_AERO(:,1,2))     -(ZFSDN_AERO(:,1,1)     -ZFSUP_AERO(:,1,1))
     PTOPSWAIAERO(:) = (ZFSDN_AERO(:,KLEV+1,2)-ZFSUP_AERO(:,KLEV+1,2))-(ZFSDN_AERO(:,KLEV+1,1)-ZFSUP_AERO(:,KLEV+1,1))

! Cloud radiative forcing with natural aerosol for direct effect 
     PSOLSWCFAERO(:,1) = (ZFSDN_AERO(:,1,2)     -ZFSUP_AERO(:,1,2))     -(ZFSDN0_AERO(:,1,2)     -ZFSUP0_AERO(:,1,2))
     PTOPSWCFAERO(:,1) = (ZFSDN_AERO(:,KLEV+1,2)-ZFSUP_AERO(:,KLEV+1,2))-(ZFSDN0_AERO(:,KLEV+1,2)-ZFSUP0_AERO(:,KLEV+1,2))
! Cloud radiative forcing with anthropogenic aerosol for direct effect 
     PSOLSWCFAERO(:,2) = 0.0
     PTOPSWCFAERO(:,2) = 0.0
! Cloud radiative forcing with no direct effect at all 
     PSOLSWCFAERO(:,3) = 0.0
     PTOPSWCFAERO(:,3) = 0.0

! LW direct anthropogenic forcing
     PSOLLWADAERO(:)  = 0.0
     PTOPLWADAERO(:)  = 0.0
     PSOLLWAD0AERO(:) = 0.0
     PTOPLWAD0AERO(:) = 0.0
     IF(ok_volcan) THEN
        PLWADAERO(:,:)  = 0.0 !--NL
     ENDIF
     
! LW indirect anthropogenic forcing
     PSOLLWAIAERO(:) = (-LWDN_AERO(:,1,2)     -LWUP_AERO(:,1,2))     -(-LWDN_AERO(:,1,1)     -LWUP_AERO(:,1,1))
     PTOPLWAIAERO(:) = (-LWDN_AERO(:,KLEV+1,2)-LWUP_AERO(:,KLEV+1,2))-(-LWDN_AERO(:,KLEV+1,1)-LWUP_AERO(:,KLEV+1,1))

ENDIF

IF (.NOT.ok_ade.AND..NOT.ok_aie) THEN

! direct anthropogenic forcing
     PSOLSWADAERO(:)  = 0.0
     PTOPSWADAERO(:)  = 0.0
     PSOLSWAD0AERO(:) = 0.0
     PTOPSWAD0AERO(:) = 0.0
     IF(ok_volcan) THEN
        PSWADAERO(:,:)  = 0.0 !--NL
     ENDIF
     
! indirect anthropogenic forcing
     PSOLSWAIAERO(:) = 0.0 
     PTOPSWAIAERO(:) = 0.0

! Cloud radiative forcing with natural aerosol for direct effect 
     PSOLSWCFAERO(:,1) = (ZFSDN_AERO(:,1,1)     -ZFSUP_AERO(:,1,1))     -(ZFSDN0_AERO(:,1,1)     -ZFSUP0_AERO(:,1,1))
     PTOPSWCFAERO(:,1) = (ZFSDN_AERO(:,KLEV+1,1)-ZFSUP_AERO(:,KLEV+1,1))-(ZFSDN0_AERO(:,KLEV+1,1)-ZFSUP0_AERO(:,KLEV+1,1))
! Cloud radiative forcing with anthropogenic aerosol for direct effect 
     PSOLSWCFAERO(:,2) = 0.0
     PTOPSWCFAERO(:,2) = 0.0
! Cloud radiative forcing with no direct effect at all 
     PSOLSWCFAERO(:,3) = 0.0
     PTOPSWCFAERO(:,3) = 0.0

! LW direct anthropogenic forcing
     PSOLLWADAERO(:)  = 0.0
     PTOPLWADAERO(:)  = 0.0
     PSOLLWAD0AERO(:) = 0.0
     PTOPLWAD0AERO(:) = 0.0
     IF(ok_volcan) THEN
        PLWADAERO(:,:)  = 0.0 !--NL
     ENDIF
     
! LW indirect anthropogenic forcing
     PSOLLWAIAERO(:) = 0.0 
     PTOPLWAIAERO(:) = 0.0

ENDIF

ENDIF

!IF (swaero_diag .OR. .NOT. AEROSOLFEEDBACK_ACTIVE) THEN
IF (.NOT. AEROSOLFEEDBACK_ACTIVE) THEN
! Cloudforcing without aerosol at all
     PSOLSWCFAERO(:,3) = (ZFSDN_AERO(:,1,5)     -ZFSUP_AERO(:,1,5))     -(ZFSDN0_AERO(:,1,5)     -ZFSUP0_AERO(:,1,5))
     PTOPSWCFAERO(:,3) = (ZFSDN_AERO(:,KLEV+1,5)-ZFSUP_AERO(:,KLEV+1,5))-(ZFSDN0_AERO(:,KLEV+1,5)-ZFSUP0_AERO(:,KLEV+1,5))

ENDIF

IF (LHOOK) CALL DR_HOOK('RECMWF_AERO',1,ZHOOK_HANDLE)
END SUBROUTINE RECMWF_AERO
