! $Id: physiq.F 1565 2011-08-31 12:53:29Z jghattas $
!#define IO_DEBUG
MODULE physiq_mod

IMPLICIT NONE

CONTAINS

!========================================================================================================================
! Adaptation de physiq.F90 de phydev pour mettre en oeuvre une interface dynamique de LMDZ/physique de MAR.
! Martin Ménégoz (Février 2014) martin.menegoz@lgge.obs.ujf-grenoble.fr
! Contact LGGE: gallee@lgge.obs.ujf-grenoble.fr (modèle MAR) ; contact LMD: Frederic.hourdin@lmd.jussieu.fr (modèle LMDZ)
!========================================================================================================================

      SUBROUTINE physiq (nlon,nlev, &
     &            debut,lafin,jD_cur, jH_cur,pdtphys, &
     &            paprs,pplay,pphi,pphis,ppresnivs, &
     &            u,v,rot,t,qx, &
     &            flxmass_w, &
     &            d_u, d_v, d_t, d_qx, d_ps &
     &            , dudyn )

      USE dimphy, only : klon,klev,klevp1
      USE infotrac_phy, only : nqtot
      USE geometry_mod, only : latitude,longitude,cell_area
      !USE comcstphy, only : rg
      USE iophy, only : histbeg_phy,histwrite_phy
      USE ioipsl, only : getin,histvert,histdef,histend,ymds2ju,ju2ymds,histsync
      USE mod_phys_lmdz_para, only : jj_nb
      USE phys_state_var_mod, only : phys_state_var_init
      USE vertical_layers_mod, ONLY: ap, bp
! Modules LMDZ utilisés pour LMDZ-MAR:
      USE control_mod      

! Modules utilisé pour MAR dans DYN_to_PHY_MAR.f90
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_S0_grd
      use Mod_SISVAT_grd
      use Mod_SISVAT_dat
      use Mod_SISVAT_gpt ! surface heat flxes in outputs
      use Mod_PHY_DY_kkl
      use Mod_PHY_CM_kkl ! pour écrire les snowfall et rainfall dans les outputs
      use Mod_PHY_CP_kkl ! pour écrire les snowfall et rainfall dans les outputs

#ifdef CPP_XIOS
      USE wxios
#endif

      IMPLICIT none

#include "dimensions.h"
! Include LMDZ utilisés ici pour MAR-LMDZ:
#include "YOMCST.h"

!======================================================================
! Arguments du moniteur general de la physique dans le modele LMDZ
!======================================================================
!
!  Arguments:
!
! nlon----input-I-nombre de points horizontaux
! nlev----input-I-nombre de couches verticales, doit etre egale a klev
! debut---input-L-variable logique indiquant le premier passage
! lafin---input-L-variable logique indiquant le dernier passage
! jD_cur       -R-jour courant a l'appel de la physique (jour julien)
! jH_cur       -R-heure courante a l'appel de la physique (jour julien)
! pdtphys-input-R-pas d'integration pour la physique (seconde)
! paprs---input-R-pression pour chaque inter-couche (en Pa)
! pplay---input-R-pression pour le mileu de chaque couche (en Pa)
! pphi----input-R-geopotentiel de chaque couche (g z) (reference sol)
! pphis---input-R-geopotentiel du sol
! ppresnivs-input_R_pressions approximat. des milieux couches ( en PA)
! u-------input-R-vitesse dans la direction X (de O a E) en m/s
! v-------input-R-vitesse Y (de S a N) en m/s
! t-------input-R-temperature (K)
! qx------input-R-humidite specifique (kg/kg) et d'autres traceurs
! d_t_dyn-input-R-tendance dynamique pour "t" (K/s)
! d_q_dyn-input-R-tendance dynamique pour "q" (kg/kg/s)
! flxmass_w -input-R- flux de masse verticale
! d_u-----output-R-tendance physique de "u" (m/s/s)
! d_v-----output-R-tendance physique de "v" (m/s/s)
! d_t-----output-R-tendance physique de "t" (K/s)
! d_qx----output-R-tendance physique de "qx" (kg/kg/s)
! d_ps----output-R-tendance physique de la pression au sol
!IM

!======================================================================
! Variables ajoutés pour l'utilisation de la physique de MAR
!======================================================================

      integer,parameter :: jjmp1=jjm+1-1/jjm
      integer,parameter :: iip1=iim+1

! Indices des traceurs (qx)
! NB: Tous les traceurs de MAR (vapeur d'eau et hydrométéores, solides et liquides) sont transportés via la variables qx (tendances: d_qx)
! Il doivent être déclarés dans un fichier traceur.def avant de lancer la simulation

      INTEGER, PARAMETER :: ivap  = 1  ! indice de traceurs pour vapeur d'eau
      INTEGER, PARAMETER :: iliq  = 2  ! indice de traceurs pour eau liquide
      INTEGER, PARAMETER :: iice  = 3  ! indice de traceurs pour eau solide
      INTEGER, PARAMETER :: icin  = 4  ! indice de traceurs pour CIN
      INTEGER, PARAMETER :: isnow = 5  ! indice de traceurs pour flocons de neige
      INTEGER, PARAMETER :: irain = 6  ! indice de traceurs pour gouttes de pluie
      INTEGER, PARAMETER :: itke  = 7  ! indice de traceurs pour tke
      INTEGER, PARAMETER :: itked = 8  ! indice de traceurs pour tke dissipation
! #cw INTEGER, PARAMETER :: iccn  = 9  ! indice de traceurs pour CCN

!======================================================================
! Arguments de la routine physiq.F90:
!======================================================================
      integer,intent(in) :: nlon ! number of atmospheric colums
      integer,intent(in) :: nlev ! number of vertical levels (should be =klev)
      real,intent(in) :: jD_cur ! current day number (Julian day)
      real,intent(in) :: jH_cur ! current time of day (as fraction of day)
      logical,intent(in) :: debut ! signals first call to physics
      logical,intent(in) :: lafin ! signals last call to physics
      real,intent(in) :: pdtphys ! physics time step (s)
      real,intent(in) :: paprs(klon,klev+1) ! interlayer pressure (Pa)
      real,intent(in) :: pplay(klon,klev) ! mid-layer pressure (Pa)
      real,intent(in) :: pphi(klon,klev) ! geopotential at mid-layer
      real,intent(in) :: pphis(klon) ! surface geopotential
      real,intent(in) :: ppresnivs(klev) ! pseudo-pressure (Pa) of mid-layers
      real,intent(in) :: u(klon,klev) ! eastward zonal wind (m/s)
      real,intent(in) :: v(klon,klev) ! northward meridional wind (m/s)
      REAL, intent(in):: rot(klon, klev) ! relative vorticity,
                                         ! in s-1, needed for frontal waves
      real,intent(in) :: t(klon,klev) ! temperature (K)
      real,intent(in) :: qx(klon,klev,nqtot) ! tracers (.../kg_air)
      real,intent(in) :: flxmass_w(klon,klev) ! vertical mass flux
      real,intent(out) :: d_u(klon,klev) ! physics tendency on u (m/s/s)
      real,intent(out) :: d_v(klon,klev) ! physics tendency on v (m/s/s)
      real,intent(out) :: d_t(klon,klev) ! physics tendency on t (K/s)
      real,intent(out) :: d_qx(klon,klev,nqtot) ! physics tendency on tracers
      real,intent(out) :: d_ps(klon) ! physics tendency on surface pressure
      real,intent(in) :: dudyn(iim+1,jjmp1,klev) ! Not used

!======================================================================
! Variables locales temporelles
!======================================================================
integer,save :: itau=0 ! counter to count number of calls to physics
!$OMP THREADPRIVATE(itau)
!real :: temp_newton(klon,klev)
integer :: k
logical, save :: first=.true.
!$OMP THREADPRIVATE(first)

!======================================================================
! Inputs/Outputs
!======================================================================
integer :: itau0
real :: zjulian
real :: dtime
integer :: nhori ! horizontal coordinate ID
integer,save :: nid_hist ! output file ID
!$OMP THREADPRIVATE(nid_hist)
integer :: zvertid ! vertical coordinate ID
integer,save :: iwrite_phys ! output every iwrite_phys physics step
!$OMP THREADPRIVATE(iwrite_phys)
integer :: iwrite_phys_omp ! intermediate variable to read iwrite_phys
real :: t_ops ! frequency of the IOIPSL operations (eg average over...)
real :: t_wrt ! frequency of the IOIPSL outputs


!======================================================================
! Arguments d'entrée pour la physique de MAR
!======================================================================

! NB : Pour plus de visibilité pour les développeurs de LMDZ et MAR, on conserve les noms de variables propres à MAR et LMDZ.
! La correspondance entre les variables se fait dans cette interface. 

! Flags pour MAR (A ranger ultérieurement dans un fichier .def, selon le protocole LMDZ).
! Attention, ces Flags pilotent l'initialisation de nombreuses variables. Toutes les configurations n'ont pas été testées.
! Ici, on trouve les déclarations utilisées pour une experience aquaplanète.
! Contacter Hubert Gallée pour connaître la signification des Flags.


       logical                    ::  FlagSV     = .TRUE.               ! Surface
       logical                    ::  FlagSV_Veg = .TRUE.               ! Végétation
       logical                    ::  FlagSV_SNo = .TRUE.               ! Neige
       logical                    ::  FlagSV_BSn = .FALSE.              !
       logical                    ::  FlagSV_KzT = .TRUE.               !
       logical                    ::  FlagSV_SWD = .TRUE.               !  (T,F) == SW INPUT is (downward,absorbed)
       ! FlagSV_SBC est utilisé pour lire des options de conditions aux limites. Ici, on le met à FALSE pour l'aquaplanète
       logical                    ::  FlagSV_SBC = .FALSE.              ! Flag pour utiliser le PHY_SISVAT_INP propre à notre experience
       logical                    ::  FlagSV_UBC = .FALSE.              !
       logical                    ::  FlagAT     = .TRUE.               !
       character(len=1)           ::  TypeAT     = 'e'                  !  Turbulence Model (e, K, L or H)
       logical                    ::  FlagAT_TKE = .TRUE.               !
       logical                    ::  FlagCM     = .TRUE.               !
!Gilles : T > F car LMDz utilise les tendances des traceurs (sauf frac.nuage)
       logical                    ::  FlagCM_UpD = .FALSE.              ! Update de la microphysique immédiat
       logical                    ::  FlagCP     = .TRUE.               ! Avec ou sans convection
       logical                    ::  FlagRT     = .TRUE.               !
       logical                    ::  FlagS0_SLO = .FALSE.              !
       logical                    ::  FlagS0_MtM = .FALSE.              !
       logical                    ::  Flag_O     = .TRUE.               !
!Gilles: FlagVR  T > F pas de sortie interessante
       logical                    ::  FlagVR     = .FALSE.              !

!======================================================================
! Arguments pour la physique de MAR :
!======================================================================

! Pas de temps par défaut. Ils sont modifiés plus bas en fonction de pdtphys défini dans les fichiers .def.

       real                       ::  dt0DYn     =   60.                !  Time Step, Dynamics                       [s] recalculé en fonction de day_step
       real                       ::  dt0_SV     =   60.                !  Time Step, SISVAT                         [s] f(pdtphys)
       real                       ::  dt0_AT     =   60.                !  Time Step, Atm_AT                         [s] f(pdtphys)
       real                       ::  dt0_CM     =   60.                !  Time Step, CMiPhy                         [s] f(pdtphys)
       real                       ::  dt0_CP     =  600.                !  Time Step, CVAmnh                         [s] f(pdtphys)
       real                       ::  dt0_RT     =  900.                !  Time Step, radCEP                         [s] f(pdtphys)
       real                       ::  dx         = 4.0e+4               !  Grid Spacing, MAX authorised for CMiPhy   [m] Limite max tolérée pour la microphysique MAR (40)
       real(kind=real8), dimension(klon)      ::  dx___HOST             !  Grid Spacing max along x axis.
       real(kind=real8), dimension(klon)      ::  dy___HOST             !  Grid   Size                               [m]
       real                       ::  DD_AxX =   90.00                  !  x-Axis Direction                          [degree] Standard (90dg = Eastward) [dg] utilise juste pour les conditions aux limite et la position du soleil.
       real, dimension(klevp1)    ::  s_HOST                            ! MAR ne fonctionne que avec des coordonnées sigma: sigma=(p-Psommet)/(psurface -Psommet)
       real, dimension(klon)      ::  sh___HOST                         !  Surface Height                            [m] =0 en aquaplanète. A calculer ultérieurement en fonction de l'environnement LMDZ.
       real(kind=real8), dimension(klon)         ::  sh_a_HOST          !  Topography Anomaly                                        [m]
       real(kind=real8), dimension(klon)         ::  slopxHOST          !  Slope, x-direction                                        [-]
       real(kind=real8), dimension(klon)         ::  slopyHOST          !  Slope, y-direction                                        [-]
       real(kind=real8), dimension(klon)         ::  slopeHOST          !  Slope                                                     [-]
       real(kind=real8), allocatable             ::  MMaskHOST(:,:)     !  Mountain Mask                                             [-]
       real, dimension(klon)                     ::  lonh_HOST          !  Longitude                                              [hour]
       real, dimension(klon)                     ::  latr_HOST          !  Latitude                                             [radian]
       real                                      ::  ptop_HOST          !  Pressure Model Top                                      [kPa] 0 in LMDZ
       real, dimension(klon,klevp1)              ::  pkta_HOST          !  Reduced  Potential Temperature                         [KX/s] !  Temperature * p ** (R/Cp) [KX]
       real, dimension(klon)                     ::  psa__HOST          !  Pressure Thickness                                      [kPa] a recalculer en fonction du geopotentiel de surface: P*=Psurface-Ptop
       real, dimension(klon,klevp1)              ::  gZa__HOST          !  Geopotential Height                                   [m2/s2]  pphi dans LMDZ
       real, dimension(klon,klevp1)              ::  gZam_HOST          !  Geopotential Height, mid-level                        [m2/s2] On fait la moyenne entre couche du dessus et couche du dessous.
       real, dimension(klon,klev)                ::  Ua___HOST          !  Wind Speed, x-direction                                 [m/s]
       real, dimension(klon,klev)                ::  Va___HOST          !  Wind Speed, y-direction                                 [m/s]
       real, dimension(klon,klev)                ::  Wa___HOST          !  Wind Speed, z-direction                                 [m/s] Non echangé entre phys et dyn dans LMDZ. Necessaire pour MAR. On le recalcule ci-aprés !
       real, dimension(klon,klevp1)              ::  qv___HOST          !  Specific Humidity                                     [kg/kg] qx(:,:,ivap)
       real, dimension(klon,klev)                ::  qw___HOST          !  Cloud Droplets Concentration                          [kg/kg] qx(:,:,iliq) 
! #cw  real, dimension(klon,klev)                ::  CCN__HOST          !  CCN            Concentration                           [-/kg] qx(:,:,iccn)
       real, dimension(klon,klev)                ::  qi___HOST          !  Cloud Crystals Concentration                          [kg/kg] qx(:,:,iice)
       real, dimension(klon,klev)                ::  CIN__HOST          !  CIN            Concentration                           [-/kg] qx(:,:,icin)
       real, dimension(klon,klev)                ::  CF___HOST          !  Cloud Fraction                                         [-/kg] Dans le futur, la CF sera transportée dans la dynamique, mais pour l'instant, on n'en a pas besoin.

       real, dimension(klon,klev)                ::  qs___HOST          !  Snow Particles Concentration                          [kg/kg] qx(:,:,isnow)
       real, dimension(klon,klev)                ::  qr___HOST          !  Rain Drops     Concentration                          [kg/kg] qx(:,:,irain)
       real, dimension(klon,klev)                ::  TKE__HOST          !  Turbulent Kinetic Energy                              [m2/s2] qx(:,:,itke)
       real, dimension(klon,klev)                ::  eps__HOST          !  Turbulent Kinetic Energy Dissipation                  [m2/s3] qx(:,:,itked)
       real, dimension(klon,klev)                ::  dpkt___dt          !  Reduced  Potential Temperature TENDENCY, ALL Contribut.[KX/s] a recalculer en fonction de d_t de LMDZ
       real, dimension(klon,klev)                ::  dua____dt          !  Wind Speed       (x-direc.)    TENDENCY, ALL Contribut.[m/s2] output commune MAR-LMD
       real, dimension(klon,klev)                ::  dva____dt          !  Wind Speed       (y-direc.)    TENDENCY, ALL Contribut.[m/s2] output commune MAR-LMD
       real, dimension(klon,klev)                ::  dqv____dt          !  Specific           Humidity    TENDENCY, ALL Contr. [kg/kg/s] output commune MAR-LMD
       real, dimension(klon,klev)                ::  dqw____dt          !  Cloud Droplets Concentration   TENDENCY, ALL Contr. [kg/kg/s] dqx
! #cw real, dimension(klon,klev)                 ::  dCw____dt          !  CCN            Concentration   TENDENCY, ALL Contr.  [1/kg/s] dqx
       real, dimension(klon,klev)                ::  dqi____dt          !  Cloud Crystals Concentration   TENDENCY, ALL Contr. [kg/kg/s] dqx
       real, dimension(klon,klev)                ::  dCi____dt          !  CIN            Concentration   TENDENCY, ALL Contr.  [1/kg/s] dqx
       real, dimension(klon,klev)                ::  dCF____dt          !  Cloud Fraction                 TENDENCY, ALL Contr. [kg/kg/s] dqx
       real, dimension(klon,klev)                ::  dqs____dt          !  Snow Particles Concentration   TENDENCY, ALL Contr. [kg/kg/s] dqx
       real, dimension(klon,klev)                ::  dqr____dt          !  Rain Drops     Concentration   TENDENCY, ALL Contr. [kg/kg/s] dqx
! #dT& real, dimension(klon,klev)                ::  dpktSV_dt          !  Reduced  Potential Temperature Tendency, SISVAT        [KX/s] 
! #dT       real, dimension(klon,klev)           ::  dpktAT_dt          !  Reduced  Potential Temperature Tendency, Atm_AT        [KX/s] 
! #dT       real, dimension(klon,klev)           ::  dqv_AT_dt          !  Specific           Humidity    TENDENCY, Atm_AT     [kg/kg/s] dqx
! #dT       real, dimension(klon,klev)           ::  dqw_AT_dt          !  Cloud Droplets Concentration   TENDENCY, Atm_AT     [kg/kg/s] dqx
! #cw real, dimension(klon,klev)                 ::  dCw_AT_dt          !  CCN            Concentration   TENDENCY, Atm_AT         [1/s] dqx
! #dT       real, dimension(klon,klev)           ::  dqi_AT_dt          !  Cloud Crystals Concentration   TENDENCY, Atm_AT     [kg/kg/s] dqx
! #dT	   real, dimension(klon,klev)            ::  dCi_AT_dt          !  CIN            Concentration   TENDENCY, Atm_AT         [1/s] dqx
! #dT	   real, dimension(klon,klev)            ::  dqs_AT_dt          !  Snow Particles Concentration   TENDENCY, Atm_AT     [kg/kg/s] dqx
! #dT	   real, dimension(klon,klev)            ::  dqr_AT_dt          !  Rain Drops     Concentration   TENDENCY, Atm_AT     [kg/kg/s] dqx
! #dT	   real, dimension(klon,klev)            ::  dpktCM_dt          !  Reduced  Potential Temperature Tendency, CMiPhy        [KX/s] dqx
! #dT	   real, dimension(klon,klev)            ::  dqv_CM_dt          !  Specific           Humidity    TENDENCY, CMiPhy     [kg/kg/s] dqx
! #dT	   real, dimension(klon,klev)            ::  dqw_CM_dt          !  Cloud Droplets Concentration   TENDENCY, CMiPhy     [kg/kg/s] dqx
! #cw real, dimension(klon,klev)                 ::  dCw_CM_dt          !  CCN            Concentration   TENDENCY, CMiPhy         [1/s] dqx
! #dT       real, dimension(klon,klev)           ::  dqi_CM_dt          !  Cloud Crystals Concentration   TENDENCY, CMiPhy     [kg/kg/s] dqx
! #dT       real, dimension(klon,klev)           ::  dCi_CM_dt          !  CIN            Concentration   TENDENCY, CMiPhy         [1/s] dqx
! #dT       real, dimension(klon,klev)           ::  dCF_CM_dt          !  Cloud Fraction                 TENDENCY, CMiPhy         [1/s] dqx
! #dT       real, dimension(klon,klev)           ::  dqs_CM_dt          !  Snow Particles Concentration   TENDENCY, CMiPhy     [kg/kg/s] dqx
! #dT       real, dimension(klon,klev)           ::  dqr_CM_dt          !  Rain Drops     Concentration   TENDENCY, CMiPhy     [kg/kg/s] dqx
! #dT       real, dimension(klon,klev)           ::  dpktCP_dt          !  Reduced  Potential Temperature TENDENCY, CVAmnh        [KX/s] 
! #dT       real, dimension(klon,klev)           ::  dqv_CP_dt          !  Specific           Humidity    TENDENCY, CVAmnh     [kg/kg/s] dqx
! #dT       real, dimension(klon,klev)           ::  dqw_CP_dt          !  Cloud Droplets Concentration   TENDENCY, CVAmnh     [kg/kg/s] dqx
! #dT       real, dimension(klon,klev)           ::  dqi_CP_dt          !  Cloud Crystals Concentration   TENDENCY, CVAmnh     [kg/kg/s] dqx
! #dT       real, dimension(klon,klev)           ::  dpktRT_dt          !  Reduced  Potential Temperature TENDENCY, radCEP        [KX/s]

       integer,save                              ::  it0EXP             ! Incremente en fonction de debut et lafin. NB: MAR fait des restarts dans le Fortran, LMDZ avec des scripts shells. On suit ici la logique LMDZ. it0EXP et it0RUN sont donc identiques.
       integer,save                              ::  it0RUN             ! Incremente en fonction de debut et lafin.
       integer                                   ::  Year_H             !  Time                                                   [year]
       integer                                   ::  Mon__H             !  Time                                                  [month]
       integer                                   ::  Day__H             !  Time                                                    [Day]
       integer                                   ::  Hour_H             !  Time                                                   [hour]
       integer                                   ::  minu_H             !  Time                                                 [minute]
       integer                                   ::  sec__H             !  Time                                                      [s]
       integer                                   ::  ixq1,i0x0,mxqq     !  Domain  Dimension: x                                      [-]
       integer                                   ::  jyq1,j0y0,myqq     !  Domain  Dimension: y                                      [-]
       !integer                                   ::  klev              !  Domain  Dimension: z                                      [-] Commun à MAR et LMDZ
       !integer                                   ::  klevp1            !  Domain  Dimension: z                                      [-] Commun à MAR et LMDZ
       integer                                   ::  mwq                !  Domain  Dimension: mosaic                                 [-]
       !integer                                   ::  klon              !  Domain  Dimension: x * y                                  [-] Commun à MAR et LMDZ
       integer                                   ::  kcolw              !  Domain  Dimension: x * y * mosaic                         [-]
       integer                                   ::  m_azim             !  Mountain Mask, nb of directions taken into account        [-]
       integer,parameter                         ::  n0pt=1             !  nombre de points qui peuvent être choisis pour le fichier PHY____.OUT
       integer, dimension(n0pt)                  ::  IOi0SV             !
       integer, dimension(n0pt)                  ::  IOj0SV             !
       real, dimension(klon)                     ::  sst__HOST          !  Ocean     FORCING (SST)                                   [K]
! #IP real, dimension(klon)                       ::  sif__HOST          !  Ocean     FORCING (Sea-Ice Fraction )                     [-]
! #AO real, dimension(ixq1:mxqq,jyq1:myqq,mwq)    ::  s_T__HOST          !  A - O    COUPLING                      n=1: Open Ocean    [K]
! #AO real, dimension(ixq1:mxqq,jyq1:myqq,mwq)    ::  Alb__HOST          !  A - O    COUPLING (Surface Albedo   )  n=2: Sea  Ice      [-]
! #AO real, dimension(ixq1:mxqq,jyq1:myqq,mwq)    ::  dSdT2HOST          !  A - O    COUPLING ( d(SH Flux) / dT )                [W/m2/K]
! #AO real, dimension(ixq1:mxqq,jyq1:myqq,mwq)    ::  dLdT2HOST          !  A - O    COUPLING ( d(SH Flux) / dT )                [W/m2/K]
! #TC integer, parameter                          ::   ntrac  =  28      !
! #TC real, dimension(ixq1:mxqq,jyq1:myqq,klev,ntrac) ::    qxTC             !  Aerosols: Atmospheric  Contentration
! #TC real, dimension(ixq1:mxqq,jyq1:myqq    ,ntrac) ::    qsTC             !  Aerosols: Near Surface Contentration
! #TC real, dimension(ixq1:mxqq,jyq1:myqq    ,ntrac) ::    uqTC             !  Aerosols: Surf.Flux

!======================================================================
! Variables locales:
!======================================================================

        real :: secondes ! pour le calcul de la date
        integer :: i ! Pour les boucles
        integer :: daysec ! Nombre de secondes dans ue journée
        integer,save :: annee_ref ! Année du début de la simulation
        !integer :: count_year ! le nombre d'années depuis le début de la simulation
        !real,parameter :: epsilon = 1.e-6 ! pour le calcul du temps

! Constantes définies comme dans MAR et utilisées ici:
        real    ::  RdA        =  287.                 !  Perfect Gas Law Constant for Dry Air [J/kg/K]
        real    ::  CpA        = 1004.                 !  Specific   Heat Capacity for Dry Air [J/kg/K]
        real    ::  kap                                !  RdA/CpA                                   [-]

! Pour l'initialisation:
        real, dimension(klon) :: Ts___HOST ! Temperature de surface pour l'initialisation.
        real, dimension(klon) :: pkt_DY_surf_tmp ! pkt de surface pour l'initialisation.
        real, dimension(klon) :: qv_HOST_surf_tmp ! qv de surface pour l'initialisation.

! Pour ecriture en netcdf
       real, dimension(klon,klev) :: cldfra   ! CF___HOST vertically inverted

!======================================================================
!! Ecriture des sorties netcdf:
!======================================================================

! NB : La version de MAR utilisée ici ne permet pas d'écrire de sorties netcdf.
! On n'utilise pas encore XIOS
! En attendant, on utilise ici le classique ioipsl, ainsi que l'outil iotd de Frédéric Hourdin.

if (debut) then ! Things to do only for the first call to physics
! load initial conditions for physics (including the grid)
  call phys_state_var_init() ! some initializations, required before calling phyetat0
  call phyetat0("startphy.nc")

! Initialize outputs:
  itau0=0
!$OMP MASTER
  iwrite_phys_omp=1 !default: output every physics timestep
  ! NB: getin() is not threadsafe; only one thread should call it.
  call getin("iwrite_phys",iwrite_phys_omp)
!$OMP END MASTER
!$OMP BARRIER
  iwrite_phys=iwrite_phys_omp
  t_ops=pdtphys*iwrite_phys ! frequency of the IOIPSL operation
  t_wrt=pdtphys*iwrite_phys ! frequency of the outputs in the file
  ! compute zjulian for annee0=1979 and month=1 dayref=1 and hour=0.0
  !CALL ymds2ju(annee0, month, dayref, hour, zjulian)
  call getin('anneeref',anneeref)
  call ymds2ju(anneeref, 1, 1, 0.0, zjulian)
  annee_ref=anneeref ! On la sauve pour le reste de la simulation
  dtime=pdtphys

! Initialisation de iophy pour l'utilisation de iotd créé par Frédéric Hourdin
  call iophys_ini ! NB: cette initialisation concerne l'ecriture de fichiers phys.nc dans PHY_MAR (rechercher iophys_ecrit)

!Gilles : replace zjulian by jD_cur (if correct init of the year)
! call histbeg_phy("histins.nc",itau0,zjulian,dtime,nhori,nid_hist)
  call histbeg_phy("histins.nc",itau0,jD_cur,dtime,nhori,nid_hist)

!$OMP MASTER

  ! define vertical coordinate
  call histvert(nid_hist,"presnivs","Vertical levels","Pa",klev, &
                ppresnivs,zvertid,'down')
  ! define variables which will be written in "histins.nc" file
  call histdef(nid_hist,'temperature','Atmospheric temperature','K', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'pplay','atmospheric pressure','K', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'ps','Surface Pressure','Pa', &
               iim,jj_nb,nhori,1,1,1,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'qx_vap','specific humidity','kg/kg', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'qx_liq','atmospheric liquid water content','kg/kg', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'qx_ice','atmospheric ice content','kg/kg', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'u','Eastward Zonal Wind','m/s', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'v','Northward Meridional Wind','m/s', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'lonh','longitude','Hours', &
               iim,jj_nb,nhori,1,1,1,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'latr','latitude','radian', &
               iim,jj_nb,nhori,1,1,1,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'sst','Prescribed SST','K', &
               iim,jj_nb,nhori,1,1,1,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  ! Tendances physiques de MAR ; pour le developpement MAR-LMDZ:
  call histdef(nid_hist,'d_t','Atmospheric temperature trends after MAR physics','K/s', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'d_u','Eastward Zonal Wind trend','m/s/s', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'d_v','Northward Meridional Wind trend','m/s/s', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  ! Définition de variables utilisée dans des sous-routines MAR
  call histdef(nid_hist,'snow','snowfall convectif + stratiform','m', &
               iim,jj_nb,nhori,1,1,1,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'rain','rainfall convectif + stratiform','m', &
               iim,jj_nb,nhori,1,1,1,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'snowCP','snowfall convective','m', &
               iim,jj_nb,nhori,1,1,1,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'rainCP','rainfall convective','m', &
               iim,jj_nb,nhori,1,1,1,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'nebulosity','CF___HOST','-', &
               iim,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  ! end definition sequence
  call histend(nid_hist)

!XIOS
#ifdef CPP_XIOS
    ! Déclaration de l'axe vertical du fichier:    
    !CALL wxios_add_vaxis("presnivs", "histins", klev, ppresnivs)

    !Déclaration du pas de temps:
    !CALL wxios_set_timestep(dtime)

    !Finalisation du contexte:
    !CALL wxios_closedef()
#endif
!$OMP END MASTER
endif ! of if (debut)

!======================================================================
! Calcul de variables environnement nécessaires à la physique de MAR
!======================================================================

! NB : Calculs à faire avant l'initialisation de la physique

! Calcul de constantes

kap    = RdA  / CpA
rpi  = 4.*ATAN(1.)
daysec = 86400.

! increment counter itau
itau=itau+1

PRINT*, 'dans physiq.F90 :'
PRINT*, 'itau=',itau

! set all tendencies to zero
PRINT*, 'On initialise les tendances à zéro:'
d_u(1:klon,1:klev)=0.
d_v(1:klon,1:klev)=0.
d_t(1:klon,1:klev)=0.
d_qx(1:klon,1:klev,1:nqtot)=0.
d_ps(1:klon)=0.

! Allocation de tableaux:
! Pour l'instant :
m_azim=1 ! Azimuths pour le calcul des masques de montagne
ALLOCATE(MMaskHOST(klon,m_azim)) !  Mountain Mask                                             [-]

! Appeller iniaqua et phyetat0 pour l'état initial avant le transfert vers les variables MAR !
! A faire pour lorsque l'on prendra en compte les continents

! Timesteps pour MAR (définis ici à partir des pas de temps LMDZ). A modifier éventuellement (via les fichiers .def)
! dt0_CP & dt0_RT modifiees ici selon Hubert avec daystep=720 ie dtdyn=2mn
! pdtphys modifie via gcm.def avec iphysiq=1 & nsplit_phys=1
dt0DYn=daysec/REAL(day_step)   !  Time Step, Dynamics                       [s]
dt0_AT=pdtphys                 !  Time Step, SISVAT                         [s]
dt0_SV=pdtphys                 !  Time Step, Atm_AT                         [s]
dt0_CM=pdtphys                 !  Time Step, CMiPhy                         [s]
dt0_CP=pdtphys*10              !  Time Step, CVAmnh                         [s]
dt0_RT=pdtphys*20              !  Time Step, radCEP                         [s]

! Correspondances de grilles: inutilisé car on utilise directement les variables LMD

! Dans un premier temps, une seule mozaique:
mwq=1

kcolw=mwq*klon

DO i=1,klon
  IF (longitude(i) .LT. 0) THEN
    lonh_HOST(i)=longitude(i)*12./rpi+24. ! from radians to hours
   ELSE
    lonh_HOST(i)=longitude(i)*12./rpi ! from radians to hours
  ENDIF
ENDDO
latr_HOST(:)=latitude(:) ! from radians to radians

!PRINT*, 'lonh_HOST(:)=',lonh_HOST(:)
!PRINT*, 'latr_HOST(:)=',latr_HOST(:)

! Niveaux verticaux:
! Attention, il faut faire un test pour vérifier si les coordonnées verticaux définis sont adéquats avec ceux de MAR!
DO k = 1, klev
  i=klev+1-k
  IF (ap(k) .NE. 0) THEN
    CALL abort_gcm('','MAR physic can be applied only with sigma levels, check vertical resolution in disvert.F90',1)
   ELSE
    s_HOST(i)=0.5*(bp(k)+bp(k+1))
  ENDIF
END DO
s_HOST(klev+1)=1 ! Dans MAR, sigma est défini depuis tout en haut jusqu'à la surface. s_HOST doit être égal à 1 en surface.

! En attendant que les routines PHY_____INI.F90, PHY_MAR_DAT.F90 et PHY_SISVAT_INP.F90 soient écrites selon un vecteur et non pas selon un tableau 2D lon-lat.
! On utilise ici la grille physique de LMDZ, qui a la taille klon (en argument de la physique): 2 + (JJM-1)*IIM (soit 1682 pour une grille dynamique de 48*36 déclarée à la compilation).

mxqq=iim
myqq=jjm+1

! Pour MAR, on a besoin des indices du premier point de la grille. Avec une grille globale, c'esty toujours (1,1):
ixq1=1
jyq1=1


! Temps:
IF (debut) THEN
  it0EXP=1
 ELSE
  it0EXP=it0EXP+1
ENDIF

IF (debut) THEN
  it0RUN=1
 ELSE
  it0RUN=it0RUN+1
ENDIF

! Pour l'instant, on ne peut démarrer une simulation seulement le 1er jour de l'année de référence:
! Attention, ju2ymds ne donne pas l'année courante, mais compte les années de simulations à partir de 0 !!! A checker !!!
CALL ymds2ju(annee_ref, 1, 1, 0.0, zjulian)
!Gilles: zjulian already counted in jD_cur if year correctly initialised
!CALL ju2ymds(jD_cur+jH_cur+zjulian,Year_H,Mon__H,Day__H,secondes)
CALL ju2ymds(jD_cur+jH_cur,Year_H,Mon__H,Day__H,secondes)
Hour_H=INT(secondes)/3600
minu_H=(INT(secondes)-Hour_H*3600)/60
!minu_H=INT((secondes-float(Hour_H)*3600.+epsilon)/60.) ! attention aux arrondis machine
sec__H=INT(secondes)-(Hour_H*3600)-(minu_H*60)

Print*,'CONTROL du temps avant d appeler la physique'
Print*,'zjulian=',zjulian
Print*,'annee_ref=',annee_ref
Print*, 'jD_cur=',JD_cur
Print*, 'jH_cur=',jH_cur
Print*, 'Year_H=',Year_H
Print*, 'Mon__H',Mon__H
Print*, 'Day__H',Day__H
Print*, 'secondes=',secondes
Print*, 'Hour_H=',Hour_H
Print*, 'minu_H=',minu_H
Print*, 'sec__H=',sec__H
Print*, 'klev=',klev
Print*, 'klevp1=',klevp1
Print*, 'nlev=',nlev
Print*, 'nlon=',nlon
Print*, 'klon=',klon
PRINT*, 'RCp=',RCp

! Coordonnées du point (en indice) pour lequel on imprime les sorties en ASCII dans le fichier PHY___________.OUT
i0x0=23
j0y0=48
! Variable utilisée juste pour les outputs:
IOi0SV(:)=1
IOj0SV(:)=1

! Variables inutilisees:
dy___HOST=0.
!CF___HOST=0. ! Pour l'instant pas utilisé, il faut tout de même l'initialiser?
dCF____dt=0.


!======================================================================
! INITIALISATIONS MAR
!======================================================================
!
! NB : les routines MAR ci-aprés doivent encore être travaillées. Ici, elles ne fonctionnent que pour une aquaplanète !

       it_RUN = it0RUN
       it_EXP = it0EXP

       IF (it0RUN.LE.1)                                              THEN

         PRINT*,'it0RUN=0; On appelle les routines d initialisation de la physique de MAR'

! Initialization of Mod_PHY____grd
! --------------------------------

! Anciennes initialisations MAR

         dxHOST = dx
         DirAxX = DD_AxX
         ixp1   = ixq1
         i_x0   = i0x0
         mxpp   = mxqq
         jyp1   = jyq1
         j_y0   = j0y0
         mypp   = myqq
         mzpp   = klevp1
         kcolp  = klon
         kcolv  = kcolw

! Initialization of Mod_PHY_S0_grd
! --------------------------------
         ! Pour l'instant, en attendant du relief :
         n_azim = m_azim

! Initialization of Mod_SISVAT_grd (Mosaic Discretization only                                      )
! -------------------------------- (       Discretization from Mod_SISVAT_dim : see PHY_SISVAT_ALLOC)

          mwp    = mwq

! Initialization of Mod_SISVAT_dat
! --------------------------------
! pour aquaplanete: FixedSST=1
          VarSST = 0.

! Initialization of Physical Constants and Allocation of main Variables
! ---------------------------------------------------------------------

! CONTROL
PRINT*, 'Appel de PHY INI'
!              **************
         CALL PHY________INI
!              **************


! INPUT DATA surface et Calcul des dependances a des Points Horizontaux autres que le Point courant
! =================================================================================================
!

! Pour l'instant, PHY_MAR_DAT n'est pas vectorisé, donc on initilise tout en dur pour l'aquaplanète.
! Ultérieurement, il faudra lires les données topo de LMDZ pour renseigner ces variables à l'initialisation.

PRINT*,'Initialisation en dur des variables topo pour une aquaplanete'
sh___HOST(:)=0.00
sh_a_HOST(:)=0.00
dx___HOST(:)   = 1.e3
dy___HOST(:)   = 1.e3
slopxHOST(:)   = 0.00
slopyHOST(:)   = 0.00
slopeHOST(:)   = 0.00
MMaskHOST(:,:) = 0.00

! CONTROL
!PRINT*, 'Appel de PHY MAR DAT'
!          ***********
!      CALL PHY_MAR_DAT(                                                &
!!          ***********
!&                 sh___HOST                                       &!  Topographie                                 [m]
!&                ,sh_a_HOST                                       &!  Topographie, Anomalie                       [m]
!&                ,dx___HOST                                       &!  Discretisation,   x                         [m]
!&                ,dy___HOST                                       &!  Discretisation,   y                         [m]
!&                ,slopxHOST                                       &!  Pente selon l'Axe x                         [-]
!&                ,slopyHOST                                       &!  Pente selon l'Axe y                         [-]
!&                ,slopeHOST                                       &!  Pente                                       [-]
!&                ,MMaskHOST                                       &!  Masque de Montagnes                         [-]
!&                ,ixq1,mxqq                                       &!
!&                ,jyq1,myqq                                       &!
!&                ,m_azim                                          &!
!&                )

! FIRST Initialization of the Surface Variables
! =============================================


! Initial Surface Temperature ( initialement dans DYN_to_PHY_MAR)
! ---------------------------

! CONTROL
PRINT*, 'Calcul de Ts___HOST dans la physique'
PRINT*, 'Initialisation de la temperature de surface avec les SST aquaplanète'
          DO i = 1,klon
            Ts___HOST(i)=273.+27.*(1-sin(1.5*latitude(i))**2)
            IF ((latitude(i).GT.1.0471975).OR.(latitude(i).LT.-1.0471975)) THEN
              Ts___HOST(i)=273.
            ENDIF

             !PRINT*, 'Initialisation de la temperature de surface avec la temperature du premier niveau'
             !Ts___HOST(i) = pkta_HOST(i,klevp1)*(psa__HOST(i)+ptop_HOST)**RCp
             !Ts___HOST(i) = t(i,klev)
          ENDDO

! CONTROL
PRINT*, 'Avant PHY SISVAT INP'
PRINT*,'MAXVAL(Ts___HOST(:)=)',maxval(Ts___HOST(:))
PRINT*,'MINVAL(Ts___HOST(:)=)',minval(Ts___HOST(:))
PRINT*, 'Appel de PHY SISVAT INP'
!              **************
          CALL PHY_SISVAT_INP(FlagSV_SBC,Ts___HOST,klon)
!              **************

PRINT*, 'Apres PHY SISVAT INP'
PRINT*,'MAXVAL(Ts___HOST(:)=)',maxval(Ts___HOST(:))
PRINT*,'MINVAL(Ts___HOST(:)=)',minval(Ts___HOST(:))

       END IF ! (it0RUN.LE.1)


! Horizontal Interactions
! ======================


! Interface: DATA for Stand-Alone Version of the Physical Package
! ===============================================================

       IF  (.NOT. FlagRT .OR.                                           &
      &     .NOT. FlagAT .OR.                                           &
      &     .NOT. FlagCM)                                           THEN
!     &     .NOT. FlagCP)                                            THEN
! CONTROL
PRINT*, 'Appel de PHY SISVAT UBC'
!                **************
           CALL   PHY_SISVAT_UBC(FlagRT,FlagAT,FlagCM,FlagCP)
!                **************

       END IF

PRINT*, 'Apres PHY_SISVAT_UBC'
Print*, 'mzp=',mzp
Print*, 'mzpp=',mzpp

!======================================================================
! Calcul des inputs de MAR à partir de LMDZ
!======================================================================

! NB : Calculs à faire aprés l'initialisation de la physique du MAR

PRINT*, 'Calcul des inputs de MAR à partir de LMDZ'

! Correspondance ou transformation de variables:
! RKAPPA=R/Cp dans un module LMD ;
! Probleme, on n'a pas la temperature de surface dans les arguments !!! a voir avec Hubert comment on fait...
!  En attendant, on prend la température du premier niveau au lieu de celle de la surface (idem pour le haut de l'atmosphère):
! NB : Attention aux unités !

IF (debut) THEN
! On initialise avec la temperature du premier niveau au lieu de la temperature de surface à la première itération:
   pkt_DY_surf_tmp(:)=t(:,1) / (pplay(:,1)/1000) ** kap
   ELSE
   pkt_DY_surf_tmp(:)=pkt_DY(:,klev+1)
ENDIF

pkta_HOST(:,klev+1)=pkt_DY_surf_tmp(:)
DO k = 1, klev
  i=klev+1-k
  pkta_HOST(:,i) = t(:,k) / (pplay(:,k)/1000) ** kap
END DO

ptop_HOST=0.
psa__HOST(:)=(paprs(:,1)-paprs(:,klev+1))/1000. ! [kPa]

! Correspondance du géopotentiel entre LMDZ et MAR :

! Aux niveaux:
DO k = 1, klev
  i=klev+1-k
  gZa__HOST(:,i)=pphi(:,k)
END DO
gZa__HOST(:,klev+1)=pphis(:)

! Aux inter-niveaux:
DO k = 1, klev
  gZam_HOST(:,k+1)=0.5*(gZa__HOST(:,k+1)+gZa__HOST(:,k))
END DO
gZam_HOST(:,1)=1.5*gZa__HOST(:,1)-0.5*gZa__HOST(:,2) ! Extrapolation au sommet de l'atmosphère
DO k = 1, klev
  i=klev+1-k
  Ua___HOST(:,i)=u(:,k)
  Va___HOST(:,i)=v(:,k)
ENDDO

! Calcul de vitesse verticale a partir de flux de masse verticale :
DO k = 1, klev
  i=klev+1-k
  !omega(i,k) = RG*flxmass_w(i,k) / cell_area(i) ! omega en Pa/s
  Wa___HOST(:,i)=flxmass_w(:,k) / cell_area(:) * (gZam_HOST(:,i+1) - gZam_HOST(:,i))/(paprs(:,k+1)-paprs(:,k)) ! Equilibre hydrostatique
END DO
Wa___HOST(:,klev)=0 ! Vitesse nulle à la surface.
Wa___HOST(:,1)=0 ! Vitesse nulle au sommet.

! Tendances:
DO k = 1, klev
  i=klev+1-k
  dpkt___dt(:,i)=d_t(:,k) / (pplay(:,k)/1000) ** kap
  dua____dt(:,i)=d_u(:,k)
  dva____dt(:,i)=d_v(:,k)
ENDDO


DO k = 1, klev
  i=klev+1-k

! Traceurs:
  qv___HOST(:,i)=qx(:,k,ivap)          !  Specific Humidity                                     [kg/kg]
  qw___HOST(:,i)=qx(:,k,iliq)          !  Cloud Droplets Concentration                          [kg/kg] 
  qi___HOST(:,i)=qx(:,k,iice)          !  Cloud Crystals Concentration                          [kg/kg]
  CIN__HOST(:,i)=qx(:,k,icin)          !  CIN            Concentration                           [-/kg]
  qs___HOST(:,i)=qx(:,k,isnow)         !  Snow Particles Concentration                          [kg/kg]
  qr___HOST(:,i)=qx(:,k,irain)         !  Rain Drops     Concentration                          [kg/kg]
  TKE__HOST(:,i)=qx(:,k,itke)          !  Turbulent Kinetic Energy                              [m2/s2]
  eps__HOST(:,i)=qx(:,k,itked)         !  Turbulent Kinetic Energy Dissipation                  [m2/s3]
  ! #cw CCN__HOST(:,:)=qx(:,:,iccn)    !  CCN            Concentration                           [-/kg]

END DO

! Dans MAR, qv(klev+1)=qv(klev)
IF (debut) THEN
  ! On initialise avec le qv du premier niveau au lieu du qv de surface à la première itération:
    qv_HOST_surf_tmp(:)=qv___HOST(:,klev)
   ELSE
    qv_HOST_surf_tmp(:)=qv__DY(:,klev+1)
ENDIF
qv___HOST(:,klev+1)=qv_HOST_surf_tmp(:)  !  Specific Humidity                                     [kg/kg]

DO k = 1, klev
  i=klev+1-k

! Dérivée des traceurs:
dqv____dt(:,i)=d_qx(:,k,ivap)        !  Specific           Humidity    TENDENCY, ALL Contr. [kg/kg/s]
dqw____dt(:,i)=d_qx(:,k,iliq)        !  Cloud Droplets Concentration   TENDENCY, ALL Contr. [kg/kg/s]
dqi____dt(:,i)=d_qx(:,k,iice)        !  Cloud Crystals Concentration   TENDENCY, ALL Contr. [kg/kg/s]
dCi____dt(:,i)=d_qx(:,k,icin)        !  CIN            Concentration   TENDENCY, ALL Contr.  [1/kg/s]
dqs____dt(:,i)=d_qx(:,k,isnow)       !  Snow Particles Concentration   TENDENCY, ALL Contr. [kg/kg/s]
dqr____dt(:,i)=d_qx(:,k,irain)       !  Rain Drops     Concentration   TENDENCY, ALL Contr. [kg/kg/s]
! #cw dCw____dt                      !  CCN            Concentration   TENDENCY, ALL Contr.  [1/kg/s]

END DO

!======================================================================
! Test Aquaplanète:
!======================================================================

! CONTROL
PRINT*, 'Calcul des SST en aquaplanète pour MAR comme dans le cas n°101 de LMDZ' ! (cf iniaqua de LMDZ)
!IF (debut) THEN
  DO i=1,klon
    sst__HOST(i)=273.+27.*(1-sin(1.5*latitude(i))**2)
    IF ((latitude(i).GT.1.0471975).OR.(latitude(i).LT.-1.0471975)) THEN
      sst__HOST(i)=273.
    ENDIF
  ENDDO
!ENDIF

!======================================================================
! CONTROL
!======================================================================
PRINT*, 'Control variables LMDZ entree de la physique'
PRINT*,'MAXVAL(pplay(:,klev))=',MAXVAL(pplay(:,klev))
PRINT*,'MINVAL(pplay(:,klev))=',MINVAL(pplay(:,klev))
PRINT*,'MAXVAL(t(:,klev))=',MAXVAL(t(:,klev))
PRINT*,'MINVAL(t(:,klev))=',MINVAL(t(:,klev))
PRINT*,'MAXVAL(pplay(:,1))=',MAXVAL(pplay(:,1))
PRINT*,'MINVAL(pplay(:,1))=',MINVAL(pplay(:,1))
PRINT*,'MAXVAL(t(:,1))=',MAXVAL(t(:,1))
PRINT*,'MINVAL(t(:,1))=',MINVAL(t(:,1))
PRINT*,'MAXVAL(pplay)=',MAXVAL(pplay)
PRINT*,'MINVAL(pplay)=',MINVAL(pplay)
PRINT*,'MAXVAL(t)=',MAXVAL(t)
PRINT*,'MINVAL(t)=',MINVAL(t)
PRINT*,'MAXVAL(pphis)=',MAXVAL(pphis)
PRINT*,'MINVAL(pphis)=',MINVAL(pphis)
PRINT*,'kap=',kap 
PRINT*, 'CONTROL variables MAR entree de la physique'
PRINT*,'MAXVAL(pkta_HOST)=',MAXVAL(pkta_HOST)
PRINT*,'MAXVAL(psa__HOST)=',MAXVAL(psa__HOST)
PRINT*,'MAXVAL(gZa__HOST)=',MAXVAL(gZa__HOST)
PRINT*,'MAXVAL(gZam_HOST)=',MAXVAL(gZam_HOST)
PRINT*,'MAXVAL(Ua___HOST)=',MAXVAL(Ua___HOST)
PRINT*,'MAXVAL(Va___HOST)=',MAXVAL(Va___HOST)
PRINT*,'MAXVAL(qv___HOST)=',MAXVAL(qv___HOST)
PRINT*,'MINVAL(pkta_HOST)=',MINVAL(pkta_HOST)
PRINT*,'MINVAL(psa__HOST)=',MINVAL(psa__HOST)
PRINT*,'MINVAL(gZa__HOST)=',MINVAL(gZa__HOST)
PRINT*,'MINVAL(gZam_HOST)=',MINVAL(gZam_HOST)
PRINT*,'MINVAL(Ua___HOST)=',MINVAL(Ua___HOST)
PRINT*,'MINVAL(Va___HOST)=',MINVAL(Va___HOST)
PRINT*,'MINVAL(qv___HOST)=',MINVAL(qv___HOST)
PRINT*,'Control iterations:'
PRINT*,'it0RUN=',it0RUN
PRINT*,'it0EXP=',it0EXP
PRINT*,'it_RUN=',it_RUN
PRINT*,'it_EXP=',it_EXP

! On choisi l'endroit où on fait les diagnostics de verification MAR (fichier ASCII):
! Pour la phase de débug, on cherche des Nan !

DO i=1,klon
  DO k=1,klev
    IF (isnan(Va___HOST(i,k))) THEN
      ikl0 = i
      PRINT*,'Attention : NaN at'
      PRINT*,'longitude=',longitude(i)*180/rpi
      PRINT*,'latitude=',latitude(i)*180/rpi
    ENDIF
  ENDDO
ENDDO

!======================================================================
! Appel de PHY_MAR
!======================================================================

! NB: Logiquement, il faut conserver cet appel tel quel, car Hubert Gallée s'engage à ne pas le changer lors des mises à jour !

! CONTROL
PRINT*, 'Appel de PHY MAR'
!PRINT*, 'FlagSV     =',FlagSV
!PRINT*, 'FlagSV_Veg =',FlagSV_Veg
!PRINT*, 'FlagSV_SNo =',FlagSV_SNo
!PRINT*, 'FlagSV_BSn =',FlagSV_BSn
!PRINT*, 'FlagSV_KzT =',FlagSV_KzT
!PRINT*, 'FlagSV_SWD =',FlagSV_SWD
!PRINT*, 'FlagSV_SBC =',FlagSV_SBC
!PRINT*, 'FlagSV_UBC =',FlagSV_UBC
!PRINT*, 'FlagAT     =',FlagAT

      CALL  PHY_MAR(FlagSV                                        &   ! FLAG  for SISVAT: (T,F) =                    (active OR NOT)
&                  ,FlagSV_Veg                                    &   ! FLAG  for SISVAT: (T,F) =       (Variable Vegetation OR NOT)
&                  ,FlagSV_SNo                                    &   ! FLAG  for SISVAT: (T,F) =         (Snow Model active OR NOT)
&                  ,FlagSV_BSn                                    &   ! FLAG  for SISVAT: (T,F) = (Blowing Snow Model active OR NOT)
&                  ,FlagSV_KzT                                    &   ! FLAG  for SISVAT: (T,F) = (pkt Turb.Transfert active OR NOT in SISVAT)
&                  ,FlagSV_SWD                                    &   ! FLAG  for SISVAT: (T,F) = (Modify SW INPUT->downward OR NOT)      
&                  ,FlagSV_SBC                                    &   ! FLAG  for SISVAT: (T,F) = (INPUT of Soil & Vege DATA OR NOT in SISVAT)
&                  ,FlagSV_UBC                                    &   ! FLAG  for SISVAT: (T,F) = (pkt UpperBC is Von Neuman OR NOT in SISVAT)
&                  ,FlagAT                                        &   ! FLAG  for Atm_AT: (T,F) = (Turbulent Transfer active OR NOT)
&                  ,TypeAT                                        &   ! TYPE  of  Atm_AT: (e= Ee Duynkerke, K= Ee Kitada, L= EL, H= Ee Huan-R)
&                  ,FlagAT_TKE                                    &   ! FLAG  for genTKE: (T,F) = (TKE-e     Model    active OR NOT)
&                  ,FlagCM                                        &   ! FLAG  for CMiPhy: (T,F) = (Cloud Microphysics active OR NOT)
&                  ,FlagCM_UpD                                    &   ! FLAG  for CMiPhy: (T,F) = (qv & hydrometeors updated OR NOT IN CMiPhy)
&                  ,FlagCP                                        &   ! FLAG  for Convection Paramet.
&                  ,FlagRT                                        &   ! FLAG  for Radiative Transfer
&                  ,FlagS0_SLO                                    &   ! FLAG  for Insolation, Surfac.Slope                 Impact  included  NEW
&                  ,FlagS0_MtM                                    &   ! FLAG  for Insolation, Surfac.Slope & Mountain Mask Impacts included  NEW
&                  ,Flag_O                                        &   ! FLAG  for OUTPUT
&                  ,FlagVR                                        &   ! FLAG  for OUTPUT for VERIFICATION
&                  ,dt0DYn                                        &   ! Time STEP between 2 CALLs of PHY_MAR                     [s] I, fix
&                  ,dt0_SV                                        &   ! Time STEP between 2 CALLs of SISVAT                      [s] I, fix
&                  ,dt0_AT                                        &   ! Time STEP between 2 CALLs of Atm_AT                      [s] I, fix
&                  ,dt0_CM                                        &   ! Time STEP between 2 CALLs of CMiPhy                      [s] I, fix
&                  ,dt0_CP                                        &   ! Time STEP between 2 CALLs of CVamnh                      [s] I, fix
&                  ,dt0_RT                                        &   ! Time STEP between 2 CALLs of radCEP                      [s] I, fix
&                  ,dx                                            &   ! Grid  Mesh size (Horizontal)                             [m] I, fix
&                  ,DD_AxX                                        &   ! Grid  x-Axis Direction                              [degree] I, fix
&                  ,s_HOST                                        &   ! Grid (Vertical)   of HOST (NORMALIZED PRESSURE assumed)  [-] I, fix
&                  ,sh___HOST                                     &   ! Topography                                               [m] I, fix
&                  ,sh_a_HOST                                     &   ! Topography Anomaly                                       [m] I, fix  NEW
&                  ,slopxHOST                                     &   ! Slope, x-direction                                       [-] I, fix  NEW
&                  ,slopyHOST                                     &   ! Slope, y-direction                                       [-] I, fix  NEW
&                  ,slopeHOST                                     &   ! Slope                                                    [-] I, fix  NEW
&                  ,MMaskHOST                                     &   ! Mountain Mask                                            [-] I, fix  NEW
&                  ,lonh_HOST                                     &   ! Longitude                                             [hour] I, fix
&                  ,latr_HOST                                     &   ! Latitude                                            [radian] I, fix
&                  ,pkta_HOST                                     &   ! Reduced Potential Temperature                           [XK] I, O
&                  ,ptop_HOST                                     &   ! Pressure, Model Top                                    [kPa] I, fix
&                  ,psa__HOST                                     &   ! Pressure  Thickness                                    [kPa] I
&                  ,gZa__HOST                                     &   ! Geopotential Height                                  [m2/s2] I
&                  ,gZam_HOST                                     &   ! Geopotential Height, mid-level                       [m2/s2] I
&                  ,Ua___HOST                                     &   ! Wind ,  x-Direction                                    [m/s] I
&                  ,Va___HOST                                     &   ! Wind ,  y-Direction                                    [m/s] I
&                  ,Wa___HOST                                     &   ! Wind ,  z-Direction                                    [m/s] I
&                  ,qv___HOST                                     &   ! Specific  Humidity                                   [kg/kg] I, O 
&                  ,qw___HOST                                     &   ! Cloud Droplets Concentration                         [kg/kg] I, O 
! #cw&             ,CCN__HOST                                     &   ! CCN            Concentration                          [-/kg] 
&                  ,qi___HOST                                     &   ! Cloud Crystals Concentration                         [kg/kg] I, O
&                  ,CIN__HOST                                     &   ! CIN            Concentration                          [-/kg] I, O 
&                  ,CF___HOST                                     &   ! Cloud Fraction                                           [-] I, O 
&                  ,qs___HOST                                     &   ! Snow Particles Concentration                         [kg/kg] I, O 
&                  ,qr___HOST                                     &   ! Rain Drops     Concentration                         [kg/kg] I, O 
&                  ,TKE__HOST                                     &   ! Turbulent Kinetic Energy                             [m2/s2] I, O
&                  ,eps__HOST                                     &   ! Turbulent Kinetic Energy Dissipation                 [m2/s3] I, O
&                  ,dpkt___dt                                     &   ! Reduced Potential Temperature TENDENCY, ALL Contribut.[KX/s]    O
&                  ,dua____dt                                     &   ! Wind Speed       (x-direc.)   TENDENCY, ALL Contribut.[m/s2]    O
&                  ,dva____dt                                     &   ! Wind Speed       (y-direc.)   TENDENCY, ALL Contribut.[m/s2]    O
&                  ,dqv____dt                                     &   ! Specific          Humidity    TENDENCY, ALL Contr. [kg/kg/s]    O
&                  ,dqw____dt                                     &   ! Cloud Droplets Concentration  TENDENCY, ALL Contr. [kg/kg/s]    O
! #cw&                  ,dCw____dt                                     &   ! CCN            Concentration  TENDENCY, ALL Contr.     [1/s]
&                  ,dqi____dt                                     &   ! Cloud Crystals Concentration  TENDENCY, ALL Contr. [kg/kg/s]    O
&                  ,dCi____dt                                     &   ! CIN            Concentration  TENDENCY, ALL Contr.     [1/s]    O
&                  ,dCF____dt                                     &   ! Cloud Fraction                TENDENCY, ALL Contr.     [1/s]    O
&                  ,dqs____dt                                     &   ! Snow Particles Concentration  TENDENCY, ALL Contr. [kg/kg/s]    O
&                  ,dqr____dt                                     &   ! Rain Drops     Concentration  TENDENCY, ALL Contr. [kg/kg/s]   O
! #dT&                  ,dpktSV_dt                                     &   ! Reduced Potential Temperature TENDENCY, SISVAT        [KX/s]  (O)
! #dT&                  ,dpktAT_dt                                     &   ! Reduced Potential Temperature TENDENCY, Atm_AT        [KX/s]  (O)
! #dT&                  ,dqv_AT_dt                                     &   ! Specific          Humidity    TENDENCY, Atm_AT     [kg/kg/s]  (O)
! #dT&                  ,dqw_AT_dt                                     &   ! Cloud Droplets Concentration  TENDENCY, Atm_AT     [kg/kg/s]  (O)
! #dT&                  ,dqi_AT_dt                                     &   ! Cloud Crystals Concentration  TENDENCY, Atm_AT     [kg/kg/s]  (O)
! #dT&                  ,dqs_AT_dt                                     &   ! Snow Particles Concentration  TENDENCY, Atm_AT     [kg/kg/s]  (O)
! #dT&                  ,dqr_AT_dt                                     &   ! Rain Drops     Concentration  TENDENCY, Atm_AT     [kg/kg/s]  (O)
! #cw&                  ,dCw_AT_dt                                     &   ! CCN            Concentration  TENDENCY, Atm_AT         [1/s]  (O)
! #dT&                  ,dCi_AT_dt                                     &   ! CIN            Concentration  TENDENCY, Atm_AT         [1/s]  (O)
! #dT&                  ,dpktCM_dt                                     &   ! Reduced Potential Temperature TENDENCY, CMiPhy        [KX/s]  (O)
! #dT&                  ,dqv_CM_dt                                     &   ! Specific          Humidity    TENDENCY, CMiPhy     [kg/kg/s]  (O)
! #dT&                  ,dqw_CM_dt                                     &   ! Cloud Droplets Concentration  TENDENCY, CMiPhy     [kg/kg/s]  (O)
! #dT&                  ,dCF_CM_dt                                     &   ! Cloud Fraction                TENDENCY, CMiPhy         [1/s]  (O)
! #dT&                  ,dqi_CM_dt                                     &   ! Cloud Crystals Concentration  TENDENCY, CMiPhy     [kg/kg/s]  (O)
! #dT&                  ,dqs_CM_dt                                     &   ! Snow Particles Concentration  TENDENCY, CMiPhy     [kg/kg/s]  (O)
! #dT&                  ,dqr_CM_dt                                     &   ! Rain Drops     Concentration  TENDENCY, CMiPhy     [kg/kg/s]  (O)
! #cw&                  ,dCw_CM_dt                                     &   ! CCN            Concentration  TENDENCY, CMiPhy         [1/s]  (O)
! #dT&                  ,dCi_CM_dt                                     &   ! CIN            Concentration  TENDENCY, CMiPhy         [1/s]  (O)
! #dT&                  ,dpktCP_dt                                     &   ! Reduced Potential Temperature TENDENCY, CVAmnh        [KX/s]  (O)
! #dT&                  ,dqv_CP_dt                                     &   ! Specific          Humidity    TENDENCY, CVAmnh     [kg/kg/s]  (O)
! #dT&                  ,dqw_CP_dt                                     &   ! Cloud Droplets Concentration  TENDENCY, CVAmnh     [kg/kg/s]  (O)
! #dT&                  ,dqi_CP_dt                                     &   ! Cloud Crystals Concentration  TENDENCY, CVAmnh     [kg/kg/s]  (O)
! #dT&                  ,dpktRT_dt                                     &   ! Reduced Potential Temperature TENDENCY, radCEP        [KX/s]  (O)
&                  ,sst__HOST                                     &   ! Ocean     FORCING (SST)                                  [K] I
! #IP&                  ,sif__HOST                                     &   ! Ocean     FORCING (Sea-Ice Fraction )                    [-] I
! #AO&                  ,s_T__HOST                                     &   ! Ocean    COUPLING (Surface Temperat.)  n=1: Open Ocean   [-] I,NEMO
! #AO&                  ,Alb__HOST                                     &   ! Ocean    COUPLING (Surface Albedo   )  n=2: Sea  Ice     [-] I,NEMO
! #AO&                  ,dSdT2HOST                                     &   ! Ocean    COUPLING ( d(SH Flux) / dT )               [W/m2/K]   O
! #AO&                  ,dLdT2HOST                                     &   ! Ocean    COUPLING ( d(LH Flux) / dT )               [W/m2/K]   O
!dead&                  ,it0EXP,it0RUN                                 &   ! Iteration
&                  ,Year_H,Mon__H,Day__H,Hour_H,minu_H,sec__H     &   ! Time
&                  ,ixq1  ,i0x0  ,mxqq                            &   ! Domain  Dimension: x
&                  ,jyq1  ,j0y0  ,myqq                            &   ! Domain  Dimension: y
&                  ,klev   ,klevp1                                   &   ! Domain  Dimension: z
&                  ,mwq                                           &   ! Domain  Dimension: mosaic
&                  ,klon                                         &   ! Domain  Dimension: x * y                                             NEW
&                  ,kcolw                                         &   ! Domain  Dimension: x * y * mosaic                                    NEW
&                  ,m_azim                                        &   ! Mountain Mask, nb of directions taken into account       [-]         NEW
&                  ,IOi0SV,IOj0SV,n0pt)                               ! Indices of OUTPUT Grid Point
!                *******

!======================================================================
! On renvoie les variables nécessaires à la dynamique de LMDZ:
!======================================================================

PRINT*, 'On est sorti de PHY_MAR, on peut actualiser les variables pour LMDZ'

d_ps(:)=0.0           ! d_ps----output-R-tendance physique de la pression au sol ! égal à zéro sur la terre (car on néglige la perte de poid de la colonne lorsqu'il pleut).

DO k = 1, klev
  i=klev+1-k

  d_u(:,k)=dua____dt(:,i) ! d_u-----output-R-tendance physique de "u" (m/s/s)
  d_v(:,k)=dva____dt(:,i) ! d_v-----output-R-tendance physique de "v" (m/s/s)
  d_t(:,k)=dpkt___dt(:,i)*(pplay(:,k)/1000) ** kap !d_t(:,:)= ! d_t-----output-R-tendance physique de "t" (K/s). Note: la pression ne change pas au cours de la physique.
  ! d_qx----output-R-tendance physique de "qx" (kg/kg/s):
  d_qx(:,k,ivap)=dqv____dt(:,i)        !  Specific           Humidity    TENDENCY, ALL Contr. [kg/kg/s]
  d_qx(:,k,iliq)=dqw____dt(:,i)        !  Cloud Droplets Concentration   TENDENCY, ALL Contr. [kg/kg/s]
  d_qx(:,k,iice)=dqi____dt(:,i)        !  Cloud Crystals Concentration   TENDENCY, ALL Contr. [kg/kg/s]
  d_qx(:,k,icin)=dCi____dt(:,i)        !  CIN            Concentration   TENDENCY, ALL Contr.  [1/kg/s]
  d_qx(:,k,isnow)=dqs____dt(:,i)       !  Snow Particles Concentration   TENDENCY, ALL Contr. [kg/kg/s]
  d_qx(:,k,irain)=dqr____dt(:,i)       !  Rain Drops     Concentration   TENDENCY, ALL Contr. [kg/kg/s]
  ! #cw dCw____dt                      !  CCN            Concentration   TENDENCY, ALL Contr.  [1/kg/s]
  cldfra(:,k)=CF___HOST(:,i)        ! cloud fraction

END DO

! Heat flux at surface (diagnostic)
      HsenSV_gpt(:) = -uts_SV_gpt(:) *1.e-3 *1005.    ! sensible heat flx in W/m2
      HLatSV_gpt(:) = -uqs_SV_gpt(:) *1.e-3 *2.5e6    ! latent heat flx in W/m2

! Question à Hubert : tendences de TKE non nécessaires ? Réponse: la TKE est mise à jour dans la physique, on n'a pas besoin de transporter la dérivée de TKE dans la dynamique.

! On désalloue les tableaux:
DEALLOCATE(MMaskHOST) ! Est-ce qu'il faut le faire ?

PRINT*, 'Control des tendances aprés la physique de MAR'

PRINT*, 'maxval(d_u)=',maxval(d_u)
PRINT*, 'maxval(d_v)=',maxval(d_v)
PRINT*, 'maxval(d_t)=',maxval(d_t)
PRINT*, 'maxval(paprs(:,1))=',maxval(paprs(:,1))
PRINT*, 'minval(d_u)=',minval(d_u)
PRINT*, 'minval(d_v)=',minval(d_v)
PRINT*, 'minval(d_t)=',minval(d_t)
PRINT*, 'minval(paprs(:,1))=',minval(paprs(:,1))
PRINT*,'Test snowCM'
PRINT*,'MAXVAL(snowCM)=',MAXVAL(snowCM)
PRINT*,'MAXVAL(rainCM)=',MAXVAL(rainCM)
PRINT*,'maxval(cldfra)=',maxval(cldfra)
PRINT*,'minval(cldfra)=',minval(cldfra)

!======================================================================
! Ecriture des sorties netcdf
!======================================================================

PRINT*, 'On ecrit les variables dans un fichier netcdf instantané'
! write some outputs:
if (debut) then ! On imprime les variables de sortie intemporelles seulement au début
  !call iophys_ini
  call histwrite_phy(nid_hist,.false.,"lonh",itau,lonh_HOST)
  call histwrite_phy(nid_hist,.false.,"latr",itau,latr_HOST)
  call histwrite_phy(nid_hist,.false.,"sst",itau,sst__HOST)
endif
if (modulo(itau,iwrite_phys)==0) then
  call histwrite_phy(nid_hist,.false.,"temperature",itau,t)
  call histwrite_phy(nid_hist,.false.,"pplay",itau,pplay)
  call histwrite_phy(nid_hist,.false.,"ps",itau,paprs(:,1))
  call histwrite_phy(nid_hist,.false.,"qx_vap",itau,qx(:,:,ivap))
  call histwrite_phy(nid_hist,.false.,"qx_liq",itau,qx(:,:,iliq))
  call histwrite_phy(nid_hist,.false.,"qx_ice",itau,qx(:,:,iice))
  call histwrite_phy(nid_hist,.false.,"u",itau,u)
  call histwrite_phy(nid_hist,.false.,"v",itau,v)
  call histwrite_phy(nid_hist,.false.,"d_t",itau,d_t)
  call histwrite_phy(nid_hist,.false.,"d_u",itau,d_u)
  call histwrite_phy(nid_hist,.false.,"d_v",itau,d_v)
  call histwrite_phy(nid_hist,.false.,"snow",itau,SnowCM)
  call histwrite_phy(nid_hist,.false.,"rain",itau,RainCM)
  call histwrite_phy(nid_hist,.false.,"snowCP",itau,snowCP)
  call histwrite_phy(nid_hist,.false.,"rainCP",itau,rainCP)
  call histwrite_phy(nid_hist,.false.,"nebulosity",itau,cldfra)
endif

!XIOS
#ifdef CPP_XIOS
!$OMP MASTER
    !On incrémente le pas de temps XIOS
    !CALL wxios_update_calendar(itau)

    !Et on écrit, avec la routine histwrite dédiée:
    !CALL histwrite_phy("temperature",t)
    !CALL histwrite_phy("u",u)
    !CALL histwrite_phy("v",v)
    !CALL histwrite_phy("ps",paprs(:,1))
!$OMP END MASTER
#endif
!
! On synchronise la fin de l'écriture du fichier à chaque pas de temps pour que le fichier soit OK même si le code plante...
CALL histsync(nid_hist)

PRINT*, 'Fin de physiq.f90'

end subroutine physiq


END MODULE physiq_mod
