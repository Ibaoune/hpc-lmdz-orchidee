MODULE YOMPHY3

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!*
!     ------------------------------------------------------------------
!     CONSTANTES PHYSIQUES REGLABLES UTILISEES POUR LES CALCULS
!     RADIATIFS :
!       BSFSA   : "BACK-SCATTERED FRACTION" SOLAIRE POUR LES AEROSOLS.
!               : SOLAR "BACK-SCATTERED FRACTION" FOR AEROSOLS.
!       BSFSI   : "BACK-SCATTERED FRACTION" SOLAIRE POUR LA GLACE.
!               : SOLAR "BACK-SCATTERED FRACTION" FOR ICE CLOUDS.
!       BSFSN   : "BACK-SCATTERED FRACTION" SOLAIRE POUR LES NUAGES.
!               : SOLAR "BACK-SCATTERED FRACTION" FOR CLOUDS.
!       BSFTA   : "BACK-SCATTERED FRACTION" THERMIQUE POUR LES AEROSOLS.
!               : THERMAL "BACK-SCATTERED FRACTION" FOR AEROSOLS.
!       BSFTI   : "BACK-SCATTERED FRACTION" THERMIQUE POUR LA GLACE.
!               : THERMAL "BACK-SCATTERED FRACTION" FOR ICE CLOUDS.
!       BSFTN   : "BACK-SCATTERED FRACTION" THERMIQUE POUR LES NUAGES.
!               : THERMAL "BACK-SCATTERED FRACTION" FOR CLOUDS.
!       EARRT   : EPAISSEUR DE L'ATMOSPHERE  /  RAYON DE LA TERRE.
!               : RATIO "DEPTH OF THE ATMOSPHERE / EARTH'S RADIUS".
!       EOASA   : COEFFICIENT D'ABSORPTION SOLAIRE PAR LES AEROSOLS.
!               : SOLAR ABSORPTION COEFFICIENT FOR AEROSOLS.
!       EOASI   : COEFFICIENT D'ABSORPTION SOLAIRE PAR LA GLACE.
!               : SOLAR ABSORPTION COEFFICIENT FOR ICE CLOUDS.
!       EOASN   : COEFFICIENT D'ABSORPTION SOLAIRE PAR LES NUAGES.
!               : SOLAR ABSORPTION COEFFICIENT FOR CLOUDS.
!       EOATA   : COEFFICIENT D'ABSORPTION THERMIQUE PAR LES AEROSOLS.
!               : THERMAL ABSORPTION COEFFICIENT FOR AEROSOLS.
!       EOATI   : COEFFICIENT D'ABSORPTION THERMIQUE PAR LA GLACE.
!               : THERMAL ABSORPTION COEFFICIENT FOR ICE CLOUDS.
!       EOATN   : COEFFICIENT D'ABSORPTION THERMIQUE PAR LES NUAGES.
!               : THERMAL ABSORPTION COEFFICIENT FOR CLOUDS.
!       EODSA   : COEFFICIENT DE DIFFUSION SOLAIRE PAR LES AEROSOLS.
!               : SOLAR SCATTERING COEFFICIENT FOR AEROSOLS.
!       EODSI   : COEFFICIENT DE DIFFUSION SOLAIRE PAR LA GLACE.
!               : SOLAR SCATTERING COEFFICIENT FOR ICE CLOUDS.
!       EODSN   : COEFFICIENT DE DIFFUSION SOLAIRE PAR LES NUAGES.
!               : SOLAR SCATTERING COEFFICIENT FOR CLOUDS.
!       EODTA   : COEFFICIENT DE DIFFUSION THERMIQUE PAR LES AEROSOLS.
!               : THERMAL SCATTERING COEFFICIENT FOR AEROSOLS.
!       EODTI   : COEFFICIENT DE DIFFUSION THERMIQUE PAR LA GLACE.
!               : THERMAL SCATTERING COEFFICIENT FOR ICE CLOUDS.
!       EODTN   : COEFFICIENT DE DIFFUSION THERMIQUE PAR LES NUAGES.
!               : THERMAL SCATTERING COEFFICIENT FOR CLOUDS.
!       EORAY   : COEFFICIENT DE DIFFUSION RAYLEIGH.
!               : RAYLEIGH SCATTERING COEFFICIENT.
!       GCA(6)  : POUR LE CALCUL "WEAK LINE" DE LA LARGEUR EQUIVALENTE.
!               : FOR THE "WEAK LINE" PART OF THE EQUIVALENT WIDTH.
!       GCB(6)  : POUR LE CALCUL "STRONG LINE" DE LA LARGEUR EQUIVAL..
!               : FOR THE "STRONG LINE" PART OF THE EQUIVALENT WIDTH.
!       GCC(6)  : POUR LE CALCUL "CONTINUUM" DE LA LARGEUR EQUIVALENTE.
!               : FOR THE "CONTINUUM" PART OF THE EQUIVALENT WIDTH.
!       GCD4    : POUR LA CONTRIBUTION "E-TYPE" A GCC(4) (H2O THERM.).
!               : FOR THE E-TYPE CONTRIBUTION TO GCC(4) (H2O THERM.).
!       GCE4    : POUR LA DEPENDANCE EN TEMPERATURE DU "E-TYPE" (GCD4).
!               : FOR THE TEMPERATURE DEPENDENCY OF THE E-TYPE (GCD4).
!       GIREC*  : JEU DE COEFFICIENTS MODULANT L'INTERACTION INFRA-ROUGE ENTRE COUCHES.
!               : COEFFICIENTS SET TO TUNE THE INFRA-RED EXCHANGE BETWEEN LAYERS.
!       QCO2    : CONCENTRATION MASSIQUE DU CO2.
!               : SPECIFIC RATIO OF CO2.
!       QLIMI   : INVERSE DU QL+QI MAXIMUM POUR UNE NEBULOSITE DE UN.
!               : INVERSE OF THE MAXIMUM QL+QI FOR CLOUD COVER ONE.
!       QLIP0   : PRESSION DE REFERENCE POUR LE CALCUL DE PQLI ET PQICE.
!               : SCALING PRESSURE FOR COMPUTING PQLI AND PQICE.
!       RII0    : VALEUR INSTANTANNEE DE LA CONST. SOLAIRE (CYCLE ANN.).
!               : INSTANTANEOUS VALUE OF THE SOLAR CONST. (ANN. CYCLE).
!       USAA    : AU NUMERATEUR DE "L'UPSCATTERED FRACTION" CAS AEROS.
!               : AT THE UPPER CASE OF THE UPSCATTERED FRACTION, AEROS.
!       USAI    : AU NUMERATEUR DE "L'UPSCATTERED FRACTION" CAS GLACE.
!               : AT THE UPPER CASE OF THE UPSCATTERED FRACTION, ICE.
!       USAN    : AU NUMERATEUR DE "L'UPSCATTERED FRACTION" CAS NUAGES.
!               : AT THE UPPER CASE OF THE UPSCATTERED FRACTION, CLOUDS.
!       USBA    : AU DENOMINATEUR DE "L'UPSCATTERED FRACT." CAS AEROS.
!               : AT THE LOWER CASE OF THE UPSCATTERED FRACTION, AEROS.
!       USBI    : AU DENOMINATEUR DE "L'UPSCATTERED FRACT." CAS GLACE.
!               : AT THE LOWER CASE OF THE UPSCATTERED FRACTION, ICE.
!       USBN    : AU DENOMINATEUR DE "L'UPSCATTERED FRACT." CAS NUAGES.
!               : AT THE LOWER CASE OF THE UPSCATTERED FRACTION, CLOUDS.
!       VDP(5,6): AU DENOMINATEUR DES FONCTIONS DE PADE POUR LES GAZ.
!               : AT THE LOWER CASE OF PADE FUNCTIONS FOR GASES.
!       VNP(5,6): AU NUMERATEUR DES FONCTIONS DE PADE POUR LES GAZ.
!               : AT THE UPPER CASE OF PADE FUNCTIONS FOR GASES.
! Parameters for cloud model:
!
!   Notations:
!     g      - asymmetry factor            (unscaled)
!     k_abs  - mass absorption coefficient (delta-scaled)
!     k_scat - mass scattering coefficient (delta-scaled)
!     delta0 - unsaturated optical depth
!     c_abs  - saturation factor for k_abs
!     c_scat - saturation factor for k_scat
!     iwc    - ice water content
!     lwc    - liquid water content
!
!   First index of FCM arrays (FCM = Fitting parameters for Cloud Model)
!   denotes spectral band:
!     1      - solar
!     2      - thermal
!
!   FCM_DEL_A(2)    : Critical value of delta0 for computation of c_abs.
!   FCM_DEL_D(2)    : Critical value of delta0 for computation of c_scat.
!   FCM_MU_A(2)     : Exponent mu for computation of c_abs.
!   FCM_MU_D(2)     : Exponent mu for computation of c_scat.
!   FCM_N_I         : Scaling exponent for iwc.
!   FCM_N_L         : Scaling exponent for lwc.
!   FCM_P_AI(2,0:3) : Pade coefficients in numerator for k_abs, ice.
!   FCM_P_AL(2,0:3) : Pade coefficients in numerator for k_abs, liquid.
!   FCM_P_DI(2,0:3) : Pade coefficients in numerator for k_scat, ice.
!   FCM_P_DL(2,0:3) : Pade coefficients in numerator for k_scat, liquid.
!   FCM_P_GI(2,0:3) : Pade coefficients in numerator for g, ice.
!   FCM_P_GL(2,0:3) : Pade coefficients in numerator for g, liquid.
!   FCM_Q_AI(2,1:3) : Pade coefficients in denominator for k_abs, ice.
!   FCM_Q_AL(2,1:3) : Pade coefficients in denominator for k_abs, liquid.
!   FCM_Q_DI(2,1:3) : Pade coefficients in denominator for k_scat, ice.
!   FCM_Q_DL(2,1:3) : Pade coefficients in denominator for k_scat, liquid.
!   FCM_Q_GI(2,1:3) : Pade coefficients in denominator for g, ice.
!   FCM_Q_GL(2,1:3) : Pade coefficients in denominator for g, liquid.
!   N_SPBAND        : Number of spectral bands.
!   N_CLOUD_MODEL   : Version of cloud model:
!                       0 - no dependency on iwc/lwc, mean saturation
!                       1 - dependency on iwc/lwc, saturation based on
!                           effective delta0 approach
!   REXP_NEB        : Scaling exponent for cloud fraction in definition
!                     of effective delta0.
INTEGER(KIND=JPIM), PARAMETER :: N_SPBAND = 2

REAL(KIND=JPRB) :: GCA(6)
REAL(KIND=JPRB) :: GCB(6)
REAL(KIND=JPRB) :: GCC(6)
REAL(KIND=JPRB) :: VDP(5,6)
REAL(KIND=JPRB) :: VNP(5,6)
REAL(KIND=JPRB) :: BSFSA
REAL(KIND=JPRB) :: BSFSI
REAL(KIND=JPRB) :: BSFSN
REAL(KIND=JPRB) :: BSFTA
REAL(KIND=JPRB) :: BSFTI
REAL(KIND=JPRB) :: BSFTN
REAL(KIND=JPRB) :: EARRT
REAL(KIND=JPRB) :: EOASA
REAL(KIND=JPRB) :: EOASI
REAL(KIND=JPRB) :: EOASN
REAL(KIND=JPRB) :: EOATA
REAL(KIND=JPRB) :: EOATI
REAL(KIND=JPRB) :: EOATN
REAL(KIND=JPRB) :: EODSA
REAL(KIND=JPRB) :: EODSI
REAL(KIND=JPRB) :: EODSN
REAL(KIND=JPRB) :: EODTA
REAL(KIND=JPRB) :: EODTI
REAL(KIND=JPRB) :: EODTN
REAL(KIND=JPRB) :: EORAY
REAL(KIND=JPRB) :: GCD4
REAL(KIND=JPRB) :: GCE4
REAL(KIND=JPRB) :: QCO2
REAL(KIND=JPRB) :: QLIMI
REAL(KIND=JPRB) :: QLIP0
REAL(KIND=JPRB) :: RII0
REAL(KIND=JPRB) :: USAA
REAL(KIND=JPRB) :: USAI
REAL(KIND=JPRB) :: USAN
REAL(KIND=JPRB) :: USBA
REAL(KIND=JPRB) :: USBI
REAL(KIND=JPRB) :: USBN
REAL(KIND=JPRB) :: GIREC1
REAL(KIND=JPRB) :: GIREC2
REAL(KIND=JPRB) :: GIREC3
REAL(KIND=JPRB) :: GIREC4
REAL(KIND=JPRB) :: FCM_DEL_A(N_SPBAND)
REAL(KIND=JPRB) :: FCM_DEL_D(N_SPBAND)
REAL(KIND=JPRB) :: FCM_MU_A(N_SPBAND)
REAL(KIND=JPRB) :: FCM_MU_D(N_SPBAND)
REAL(KIND=JPRB) :: FCM_N_I
REAL(KIND=JPRB) :: FCM_N_L
REAL(KIND=JPRB) :: FCM_P_AI(N_SPBAND,0:3)
REAL(KIND=JPRB) :: FCM_P_AL(N_SPBAND,0:3)
REAL(KIND=JPRB) :: FCM_P_DI(N_SPBAND,0:3)
REAL(KIND=JPRB) :: FCM_P_DL(N_SPBAND,0:3)
REAL(KIND=JPRB) :: FCM_P_GI(N_SPBAND,0:3)
REAL(KIND=JPRB) :: FCM_P_GL(N_SPBAND,0:3)
REAL(KIND=JPRB) :: FCM_Q_AI(N_SPBAND,1:3)
REAL(KIND=JPRB) :: FCM_Q_AL(N_SPBAND,1:3)
REAL(KIND=JPRB) :: FCM_Q_DI(N_SPBAND,1:3)
REAL(KIND=JPRB) :: FCM_Q_DL(N_SPBAND,1:3)
REAL(KIND=JPRB) :: FCM_Q_GI(N_SPBAND,1:3)
REAL(KIND=JPRB) :: FCM_Q_GL(N_SPBAND,1:3)
REAL(KIND=JPRB) :: REXP_NEB
!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(bsfsa,bsfsi,bsfsn,bsfta,bsfti,bsftn,earrt,eoasa,eoasi,eoasn,eoata,eoati,eoatn,eodsa,eodsi)
!$OMP THREADPRIVATE(eodsn,eodta,eodti,eodtn,eoray,fcm_del_a,fcm_del_d,fcm_mu_a,fcm_mu_d,fcm_n_i,fcm_n_l,fcm_p_ai)
!$OMP THREADPRIVATE(fcm_p_al,fcm_p_di,fcm_p_dl,fcm_p_gi,fcm_p_gl,fcm_q_ai,fcm_q_al,fcm_q_di,fcm_q_dl,fcm_q_gi)
!$OMP THREADPRIVATE(fcm_q_gl,gca,gcb,gcc,gcd4,gce4,girec1,girec2,girec3,girec4,qco2,qlimi,qlip0,rexp_neb,rii0)
!$OMP THREADPRIVATE(usaa,usai,usan,usba,usbi,usbn,vdp,vnp)
END MODULE YOMPHY3
