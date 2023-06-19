MODULE YOMFPC

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARFPOS, ONLY : JPOS2DF  ,JPOS3DF  ,JPOS3H   ,JPOS3P   ,JPOS3PV  , &
 & JPOS3S   ,JPOS3TH  ,JPOSCFU  ,JPOSDIR  ,JPOSLEN  ,JPOSDOM  , &
 & JPOSSGP  ,JPOSXFU  

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Scientific and technical variables

! === TECHNICAL VARIABLES ===

!     CFPDIR : prefix (path) for the output files
!     CFPIDEN: identificator of the output files
!     CFPFMT : format of the output files ('MODEL','GAUSS','LELAM' or 'LALON')
!     CFPDOM, C1FPDOM : names of the subdomains
!     CFP3DF, C1FP3DF : names 3D dynamics fields to compute
!     CFP3DFS,C1FP3DFS: names 3D dynamics fields to compute on model levels
!     CFP3DFP,C1FP3DFP: names 3D dynamics fields to compute on Pressure levels
!     CFP3DFH,C1FP3DFH: names 3D dynamics fields to compute on H levels
!     CFP3DFT,C1FP3DFT: names 3D dynamics fields to compute on THETA levels
!     CFP3DFV,C1FP3DFV: names 3D dynamics fields to compute on PV levels
!     CFP2DF, C1FP2DF : names 2D dynamics fields to compute
!     CFPPHY, C1FPPHY : names of physical fields to be post-processed
!     CFPCFU, C1FPCFU : names of cumulated fluxes fields to be post-processed
!     CFPXFU, C1FPXFU : names of instantaneous fluxes fields to be post-proc
!     Remark: C1FPXXX=CFPXXX(1)//CFPXXX(2)//...//CFPXXX(JPOSXXX)
!      for XXX=DOM,3DF,2DF,PHY,CFU,XFU.

!     Variables used for ECMWF 

!     MFP3DFS: Gribcodes of  3D dynamics fields to compute on model levels
!     MFP3DFH: Gribcodes of  3D dynamics fields to compute on H levels
!     MFP3DFT: Gribcodes of  3D dynamics fields to compute on THETA levels
!     MFP3DFV: Gribcodes of  3D dynamics fields to compute on PV levels
!     MFP3DFP: Gribcodes of  3D dynamics fields to compute on P levels
!     MFP2DF : Gribcodes of  2D dynamics fields to compute
!     MFPPHY : Gribcodes of  physical fields to processed 

!     RFP3P  : post-processing pressure levels
!     RFP3H  : post-processing height (above orography) levels
!     RFP3TH : post-processing potential temperature levels
!     RFP3PV : post-processing potential vorticity levels
!     NRFP3S : post-processing eta levels (CONF. 927 only)

!     LFPCNT : control varying output variables according to time step
!     LTRACEFP: trace for Full-POS
!     NFPGRIB : level of GRIB coding in output file ARPEGE/ALADIN.
!               0 : no packing at all
!               1 : standart GRIB encoding
!               2 : modified GRIB encoding 

!     NFPDOM : useful dimension of CFPDOM
!     NFP3DF : useful dimension of CFP3DF
!     NFP3DFS: useful dimension of CFP3DFS
!     NFP3DFH: useful dimension of CFP3DFH
!     NFP3DFT: useful dimension of CFP3DFT
!     NFP3DFV: useful dimension of CFP3DFV
!     NFP3DFP: useful dimension of CFP3DFP
!     NFP2DF : useful dimension of CFP2DF
!     NFPPHY : useful dimension of CFPPHY
!     NFPCFU : useful dimension of CFPCFU
!     NFPXFU : useful dimension of CFPXFU
!     NFP3P  : useful dimension of RFP3P
!     NFP3H  : useful dimension of RFP3H
!     NFP3TH : useful dimension of RFP3TH
!     NFP3PV : useful dimension of RFP3PV
!     NFP3S  : useful dimension of NRFP3S
!     NFPXLEV: max. number of pp. levels
!     MFP3DYN : max. number of 3D-dynamic fields neeeded for LFPART2
!     MFP2DYN : max. number of 2D-dynamic fields neeeded for LFPART2
!     NFPDPHY : maximum number of physical fields

!     The next keys force variables to be primitive for the post-processing :
!     LFPNHPD : .TRUE. if pp of NH pressure departure
!     LFPNHVD : .TRUE. if pp of NH vertical divergence
!     LFPNHVW : .TRUE. if pp of NH true vertical velocity

!     LFPLOSP : .TRUE. = Fill Ps array with Log(Ps)
!     LFPRH100 : .TRUE. to convert relative humidity in percent

! === SCIENTIFIC VARIABLES ===

!     LFPSPEC: =.T. if pp. dyn. fields are written out as spectral coefficients
!              =.F. if pp. dyn. fields are written out as grid point values.
!     LFITP  : =1 if pp. fields on P levels should be fited ; =0 otherwise
!     LFITT  : =1 if pp. fields on THETA levels should be fited ; =0 otherwise
!     LFITV  : =1 if pp. fields on PV levels should be fited ; =0 otherwise

!     NFPCLI : usage level for climatology
!              =0 no climatology
!              =1 orography and land-sea mask of output only
!              =2 all available climatological fields of the current month
!              =3 shifting mean from the climatological fields of the current
!                 month to the ones of the closest month

!     NFPINDYN : type of interpolations for dynamical fields
!     NFPINPHY : type of interpolations for  physical fields
!                4 for bilinear, 12 for quadratic

!     LFPQ   : =.TRUE. if specific humidity is interpolated
!            : =.FALSE. if relative humidity is interpolated
!     LASQ   : =.TRUE. if ASQ set to 80% saturation

!     WSXI   : max. surface moisture in input
!     WDXI   : max. deep soil moisture in input
!     WSXO   : max. surface moisture in output
!     WDXO   : max. deep soil moisture in output

!     FPBL   : Critical Thickness of PBL
!     RFPCORR: Critical orography difference for correcting surface
!              temperature through standart profile.
!     RFPCSAB: Critical sand percentage difference for computing relative 
!              soil moisture in ISBA
!     RFPCD2 : Critical soil depth difference for computing relative
!              soil moisture in ISBA
!     RFPVCAP: Minimum pressure of model level to provide an equatorial
!              cap in the computation of variables on constant PV surfaces
!     NFPLNPR : 0 => conventional formulation of (delta P) : ln(P(l)/P(l-1))
!               1 => formulation of (delta P) used in non hydrostatic model,
!                    i.e. (P(l)-P(l-1))/SQRT(P(l)*P(l-1))
!     LFPMOIS: Month allowed for climatology usage : 
!               .F. => month of the model (forecast)
!               .T. => month of the file
!     NFPINCR : Bogussing indicator ; =0 for no bogussing ; =1 for bogussing
!     NFPLAKE : To overwrite created lakes or islands by specific data :
!               0 => do not overwrite
!              -1 => overwrite with rhoughly interpolated data
!              +1 => overwrite with climatology
!     NFPCAPE : Kind of computation for CAPE & CIN : 
!               1 => from bottom model layer
!               2 => from the most unstable layer
!               3 => from mto standart height (2 meters) as recomputed values
!               4 => from mto standart height (2 meters) out of fluxes
!                    (used for analysis)
!     NFPSURFEX : Subcontract surface fields to SURFEX
!               0 => no subcontract
!               1 => transform native arp/ald surface fields to surfex fields
!                    and write out by surfex
!     NFPMASK   : number of masks for the interpolation of surface fields
!                 0 => no mask
!                 1 => land-sea mask
!                 2 => land mask, sea mask
!     LMOCONVAR  : if true mocon is computed from N-1 level parameters (for varpack use)

CHARACTER (LEN = JPOSDIR) ::  CFPDIR
CHARACTER (LEN = 5) ::  CFPFMT
CHARACTER (LEN = 8) ::  CFPIDEN
CHARACTER (LEN = JPOSLEN) ::  CFPDOM(JPOSDOM)
CHARACTER (LEN = 12) ::  CFP3DF(JPOS3DF)
CHARACTER (LEN = 12) ::  CFP3DFS(JPOS3DF)
CHARACTER (LEN = 12) ::  CFP3DFP(JPOS3DF)
CHARACTER (LEN = 12) ::  CFP3DFH(JPOS3DF)
CHARACTER (LEN = 12) ::  CFP3DFT(JPOS3DF)
CHARACTER (LEN = 12) ::  CFP3DFV(JPOS3DF)
CHARACTER (LEN = 16) ::  CFP2DF(JPOS2DF)
CHARACTER (LEN = 16) ::  CFPPHY(JPOSSGP)
CHARACTER (LEN = 16) ::  CFPCFU(JPOSCFU)
CHARACTER (LEN = 16) ::  CFPXFU(JPOSXFU)
CHARACTER (LEN = 12*JPOS3DF) :: C1FP3DF
CHARACTER (LEN = 12*JPOS3DF) :: C1FP3DFS
CHARACTER (LEN = 12*JPOS3DF) :: C1FP3DFP
CHARACTER (LEN = 12*JPOS3DF) :: C1FP3DFH
CHARACTER (LEN = 12*JPOS3DF) :: C1FP3DFT
CHARACTER (LEN = 12*JPOS3DF) :: C1FP3DFV
CHARACTER (LEN = 16*JPOS2DF) :: C1FP2DF
CHARACTER (LEN = 16*JPOSSGP) :: C1FPPHY
CHARACTER (LEN = 16*JPOSCFU) :: C1FPCFU
CHARACTER (LEN = 16*JPOSXFU) :: C1FPXFU
CHARACTER (LEN = JPOSLEN*JPOSDOM) :: C1FPDOM
INTEGER(KIND=JPIM) :: MFP3DFS(JPOS3DF)
INTEGER(KIND=JPIM) :: MFP3DFP(JPOS3DF)
INTEGER(KIND=JPIM) :: MFP3DFH(JPOS3DF)
INTEGER(KIND=JPIM) :: MFP3DFT(JPOS3DF)
INTEGER(KIND=JPIM) :: MFP3DFV(JPOS3DF)
INTEGER(KIND=JPIM) :: MFP2DF(JPOS2DF)
INTEGER(KIND=JPIM) :: MFPPHY(JPOSSGP)

REAL(KIND=JPRB) :: RFP3P(JPOS3P)
REAL(KIND=JPRB) :: RFP3H(JPOS3H)
REAL(KIND=JPRB) :: RFP3TH(JPOS3TH)
REAL(KIND=JPRB) :: RFP3PV(JPOS3PV)
INTEGER(KIND=JPIM) :: NRFP3S(JPOS3S)

INTEGER(KIND=JPIM) :: NFPDOM
INTEGER(KIND=JPIM) :: NFP3DF
INTEGER(KIND=JPIM) :: NFP3DFS
INTEGER(KIND=JPIM) :: NFP3DFP
INTEGER(KIND=JPIM) :: NFP3DFH
INTEGER(KIND=JPIM) :: NFP3DFT
INTEGER(KIND=JPIM) :: NFP3DFV
INTEGER(KIND=JPIM) :: NFP2DF
INTEGER(KIND=JPIM) :: NFPPHY
INTEGER(KIND=JPIM) :: NFPCFU
INTEGER(KIND=JPIM) :: NFPXFU
INTEGER(KIND=JPIM) :: NFP3P
INTEGER(KIND=JPIM) :: NFP3H
INTEGER(KIND=JPIM) :: NFP3TH
INTEGER(KIND=JPIM) :: NFP3PV
INTEGER(KIND=JPIM) :: NFP3S
INTEGER(KIND=JPIM) :: NFPGRIB
INTEGER(KIND=JPIM) :: NFPCLI
INTEGER(KIND=JPIM) :: MFP3DYN
INTEGER(KIND=JPIM) :: MFP2DYN
INTEGER(KIND=JPIM) :: NFPLNPR
INTEGER(KIND=JPIM) :: NFPINCR
INTEGER(KIND=JPIM) :: NFPINDYN
INTEGER(KIND=JPIM) :: NFPINPHY
INTEGER(KIND=JPIM) :: NFPCAPE
REAL(KIND=JPRB) :: WSXI
REAL(KIND=JPRB) :: WDXI
REAL(KIND=JPRB) :: WSXO
REAL(KIND=JPRB) :: WDXO
REAL(KIND=JPRB) :: FPBL
REAL(KIND=JPRB) :: RFPCORR
REAL(KIND=JPRB) :: RFPCSAB
REAL(KIND=JPRB) :: RFPCD2
REAL(KIND=JPRB) :: RFPVCAP
LOGICAL :: LFPSPEC
LOGICAL :: LFITP
LOGICAL :: LFITT
LOGICAL :: LFITV
LOGICAL :: LFPQ
LOGICAL :: LTRACEFP
LOGICAL :: LFPCNT
LOGICAL :: LASQ
LOGICAL :: LFPNHPD
LOGICAL :: LFPNHVD
LOGICAL :: LFPNHVW
LOGICAL :: LFPMOIS
LOGICAL :: LFPLOSP
LOGICAL :: LFPRH100
LOGICAL :: LMOCONVAR
INTEGER(KIND=JPIM) :: NFPLAKE
INTEGER(KIND=JPIM) :: NFPXLEV
INTEGER(KIND=JPIM) :: NFPDPHY
INTEGER(KIND=JPIM) :: NFPSURFEX
INTEGER(KIND=JPIM) :: NFPMASK
!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(c1fp2df,c1fp3df,c1fp3dfh,c1fp3dfp,c1fp3dfs,c1fp3dft,c1fp3dfv,c1fpcfu,c1fpdom,c1fpphy,c1fpxfu)
!$OMP THREADPRIVATE(cfp2df,cfp3df,cfp3dfh,cfp3dfp,cfp3dfs,cfp3dft,cfp3dfv,cfpcfu,cfpdir,cfpdom,cfpfmt,cfpiden)
!$OMP THREADPRIVATE(cfpphy,cfpxfu,fpbl,lasq,lfitp,lfitt,lfitv,lfpcnt,lfplosp,lfpmois,lfpnhpd,lfpnhvd,lfpnhvw)
!$OMP THREADPRIVATE(lfpq,lfprh100,lfpspec,lmoconvar,ltracefp,mfp2df,mfp2dyn,mfp3dfh,mfp3dfp,mfp3dfs,mfp3dft)
!$OMP THREADPRIVATE(mfp3dfv,mfp3dyn,mfpphy,nfp2df,nfp3df,nfp3dfh,nfp3dfp,nfp3dfs,nfp3dft,nfp3dfv,nfp3h,nfp3p)
!$OMP THREADPRIVATE(nfp3pv,nfp3s,nfp3th,nfpcape,nfpcfu,nfpcli,nfpdom,nfpdphy,nfpgrib,nfpincr,nfpindyn,nfpinphy)
!$OMP THREADPRIVATE(nfplake,nfplnpr,nfpmask,nfpphy,nfpsurfex,nfpxfu,nfpxlev,nrfp3s,rfp3h,rfp3p,rfp3pv,rfp3th)
!$OMP THREADPRIVATE(rfpcd2,rfpcorr,rfpcsab,rfpvcap,wdxi,wdxo,wsxi,wsxo)
END MODULE YOMFPC
