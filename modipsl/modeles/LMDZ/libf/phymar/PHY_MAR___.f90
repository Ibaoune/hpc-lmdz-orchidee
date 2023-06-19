      subroutine PHY_MAR                                               &

!------------------------------------------------------------------------------+
!                                                         Mon  1-Jul-2013  MAR |
!     subroutine PHY_MAR is the MAR PHYsics Driver                             |
!     interfaces HOST   variables                                              |
!            and SISVAT variables                                              |
!                                                                              |
!     Applied to: MARthusalem               (variables in MAR***.inc files)    |
!                                                                              |
!                                                                              |
! # OPTIONS: #dT  Distinction among  Tendencies of MAR  Physical Parametr.     |
! # ^^^^^^^^ #cw  Cloud Condensation Nuclei (CCNw) Microphysics Activation     |
!                                                                              |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,               Mon  1-Jul-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

     &                  (FlagSV                                        &   ! FLAG  for SISVAT: (T,F) =                    (active OR NOT)
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
     &                  ,FlagS0_SLO                                    &   ! FLAG  for Insolation, Surfac.Slope                 Impact  included NEW
     &                  ,FlagS0_MtM                                    &   ! FLAG  for Insolation, Surfac.Slope & Mountain Mask Impacts included NEW
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
     &                  ,sh_a_HOST                                     &   ! Topography Anomaly                                       [m] I, fix NEW
     &                  ,slopxHOST                                     &   ! Slope, x-direction                                       [-] I, fix NEW
     &                  ,slopyHOST                                     &   ! Slope, y-direction                                       [-] I, fix NEW
     &                  ,slopeHOST                                     &   ! Slope                                                    [-] I, fix NEW
     &                  ,MMaskHOST                                     &   ! Mountain Mask                                            [-] I, fix NEW
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
! #cw&                  ,CCN__HOST                                     &   ! CCN            Concentration                          [-/kg]
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
     &                  ,mzq   ,mzqq                                   &   ! Domain  Dimension: z
     &                  ,mwq                                           &   ! Domain  Dimension: mosaic
     &                  ,kcolq                                         &   ! Domain  Dimension: x * y 
     &                  ,kcolw                                         &   ! Domain  Dimension: x * y * mosaic
     &                  ,m_azim                                        &   ! Mountain Mask, nb of directions taken into account       [-]
     &                  ,IOi0SV,IOj0SV,n0pt)                               ! Indices of OUTPUT Grid Point

!------------------------------------------------------------------------------+
!                                                         Sat 29-Jun-2013  MAR |
!     subroutine PHY_MAR is the MAR PHYsics Driver                             |
!     interfaces HOST   variables                                              |
!            and SISVAT variables                                              |
!                                                                              |
!     Applied to: MARthusalem               (variables in MAR***.inc files)    |
!                                                                              |
!                                                                              |
! # OPTIONS: #dT  Distinction among  Tendencies of MAR  Physical Parametr.     |
! # ^^^^^^^^ #cw  Cloud Condensation Nuclei (CCNw) Microphysics Activation     |
!                                                                              |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,               Sat 29-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY____kkl
      use Mod_PHY_CM_ctr
      use Mod_PHY_S0_ctr
      use Mod_SISVAT_ctr
      use Mod_PHY_CM_dat
      use Mod_PHY_AT_grd
      use Mod_PHY_CM_grd
      use Mod_PHY_CP_grd
      use Mod_PHY_RT_grd
      use Mod_PHY_S0_grd
      use Mod_SISVAT_grd
      use Mod_PHY_DY_kkl
      use Mod_PHY_AT_kkl
      use Mod_PHY_CM_kkl
      use Mod_PHY_CP_kkl
      use Mod_PHY_RT_kkl
      use Mod_PHY_S0_kkl
      use Mod_SISVAT_kkl
      use Mod_SISVAT_gpt

      IMPLICIT NONE

      logical                                               ::  FlagSV             !  Flag         (SISVAT)
      logical                                               ::  FlagSV_Veg         !  Flag         (SISVAT: Vegetation)
      logical                                               ::  FlagSV_SNo         !  Flag         (SISVAT: Surface * )
      logical                                               ::  FlagSV_BSn         !  Flag         (SISVAT: Blowing * )
      logical                                               ::  FlagSV_KzT         !  Flag         (SISVAT: d(KdT/dz)/dz)
      logical                                               ::  FlagSV_SWD         !  Flag: T/F :  (SISVAT: SW=Down/Abs.)
      logical                                               ::  FlagSV_SBC         !  Flag: T/F :  (SISVAT: SBC=INP/FIX.)
      logical                                               ::  FlagSV_UBC         !  Flag: T/F :  (SISVAT: UBC=VonN/Dr.)
      logical                                               ::  FlagAT             !  Flag         (Turbulent Transfer)
      character(len=1)                                      ::  TypeAT             !  Type         (Turbulent Transfer)
      logical                                               ::  FlagAT_TKE         !  Flag         (Turbulent Transfer, TKE-e Model: ON / OFF)
      logical                                               ::  FlagCM             !  Flag         (Cloud Microphysics)
      logical                                               ::  FlagCM_UpD         !  Flag: T/F :  (Cloud Microphysics: Update in/out PHY_MAR)
      logical                                               ::  FlagCP             !  Flag         (Convection Param. )
      logical                                               ::  FlagRT             !  Flag         (Radiative Transfer)
      logical                                               ::  FlagS0_SLO         !  FLAG         (Insolation, Surfac.Slope                )
      logical                                               ::  FlagS0_MtM         !  FLAG         (Insolation, Surfac.Slope & Mountain Mask)

      logical                                               ::  Flag_O             !  Flag         (OUTPUT)
      logical                                               ::  FlagVR             !  Flag         (OUTPUT for VERIFICATION)

      real                                                  ::  dt0DYn             !  Time   Step  (DYnamics, the shortest)                     [s]
      real                                                  ::  dt0_SV             !  Time   Step  (SISVAT)                                     [s]
      real                                                  ::  dt0_AT             !  Time   Step  (Atmo Turb.)                                 [s]
      real                                                  ::  dt0_CM             !  Time   Step  (Cloud Mic.)                                 [s]
      real                                                  ::  dt0_CP             !  Time   Step  (Convection)                                 [s]
      real                                                  ::  dt0_RT             !  Time   Step  (Radiat.Tr.)                                 [s]
      real                                                  ::  dx                 !  Grid   Size                                               [m]
      real                                                  ::  DD_AxX             !  x-Axis Direction                                     [degree]
      real, dimension(mzqq)                                 ::  s_HOST             !  Vertical Coordinate                                       [-]
      real            , dimension(kcolq)                    ::  sh___HOST          !  Topography                                                [m]
      real(kind=real8), dimension(kcolq)                    ::  sh_a_HOST          !  Topography Anomaly                                        [m]
      real(kind=real8), dimension(kcolq)                    ::  slopxHOST          !  Slope, x-direction                                        [-]
      real(kind=real8), dimension(kcolq)                    ::  slopyHOST          !  Slope, y-direction                                        [-]
      real(kind=real8), dimension(kcolq)                    ::  slopeHOST          !  Slope                                                     [-]
      real(kind=real8), dimension(kcolq,m_azim)             ::  MMaskHOST          !  Mountain Mask                                             [-]
      real, dimension(kcolq)                                ::  lonh_HOST          !  Longitude                                              [hour]
      real, dimension(kcolq)                                ::  latr_HOST          !   Latitude                                            [radian]
      real                                                  ::  ptop_HOST          !  Pressure Model Top                                      [kPa]

      real, dimension(kcolq,mzqq)                           ::  pkta_HOST          !  Reduced  Potential Temperature                         [KX/s]
      real, dimension(kcolq)                                ::  psa__HOST          !  Pressure Thickness                                      [kPa]
      real, dimension(kcolq,mzqq)                           ::  gZa__HOST          !  Geopotential Height                                   [m2/s2]
      real, dimension(kcolq,mzqq)                           ::  gZam_HOST          !  Geopotential Height, mid-level                        [m2/s2]
      real, dimension(kcolq,mzq)                            ::  Ua___HOST          !  Wind Speed, x-direction                                 [m/s]
      real, dimension(kcolq,mzq)                            ::  Va___HOST          !  Wind Speed, y-direction                                 [m/s]
      real, dimension(kcolq,mzq) ,INTENT(IN)                ::  Wa___HOST          !  Wind Speed, z-direction                                 [m/s]
      real, dimension(kcolq,mzqq)                           ::  qv___HOST          !  Specific Humidity                                     [kg/kg]
      real, dimension(kcolq,mzq)                            ::  qw___HOST          !  Cloud Droplets Concentration                          [kg/kg]
      real, dimension(kcolq,mzq)                            ::  CCN__HOST          !  CCN            Concentration                           [-/kg]
      real, dimension(kcolq,mzq)                            ::  qi___HOST          !  Cloud Crystals Concentration                          [kg/kg]
      real, dimension(kcolq,mzq)                            ::  CIN__HOST          !  CIN            Concentration                           [-/kg]
      real, dimension(kcolq,mzq)                            ::  CF___HOST          !  Cloud Fraction                                         [-/kg]
      real, dimension(kcolq,mzq)                            ::  qs___HOST          !  Snow Particles Concentration                          [kg/kg]
      real, dimension(kcolq,mzq)                            ::  qr___HOST          !  Rain Drops     Concentration                          [kg/kg]
      real, dimension(kcolq,mzq)                            ::  TKE__HOST          !  Turbulent Kinetic Energy                              [m2/s2]
      real, dimension(kcolq,mzq)                            ::  eps__HOST          !  Turbulent Kinetic Energy Dissipation                  [m2/s3]

      real, dimension(kcolq,mzq)                            ::  dpkt___dt          !  Reduced  Potential Temperature TENDENCY, ALL Contribut.[KX/s]
      real, dimension(kcolq,mzq)                            ::  dua____dt          !  Wind Speed       (x-direc.)    TENDENCY, ALL Contribut.[m/s2]
      real, dimension(kcolq,mzq)                            ::  dva____dt          !  Wind Speed       (y-direc.)    TENDENCY, ALL Contribut.[m/s2]
      real, dimension(kcolq,mzq)                            ::  dqv____dt          !  Specific           Humidity    TENDENCY, ALL Contr. [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dqw____dt          !  Cloud Droplets Concentration   TENDENCY, ALL Contr. [kg/kg/s]
! #cw real, dimension(kcolq,mzq)                            ::  dCw____dt          !  CCN            Concentration   TENDENCY, ALL Contr.  [1/kg/s]
      real, dimension(kcolq,mzq)                            ::  dqi____dt          !  Cloud Crystals Concentration   TENDENCY, ALL Contr. [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dCi____dt          !  CIN            Concentration   TENDENCY, ALL Contr.  [1/kg/s]
      real, dimension(kcolq,mzq)                            ::  dCF____dt          !  Cloud Fraction                 TENDENCY, ALL Contr. [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dqs____dt          !  Snow Particles Concentration   TENDENCY, ALL Contr. [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dqr____dt          !  Rain Drops     Concentration   TENDENCY, ALL Contr. [kg/kg/s]

      real, dimension(kcolq,mzq)                            ::  dpktSV_dt          !  Reduced  Potential Temperature Tendency, SISVAT        [KX/s]
      real, dimension(kcolq,mzq)                            ::  dpktAT_dt          !  Reduced  Potential Temperature Tendency, Atm_AT        [KX/s]
      real, dimension(kcolq,mzq)                            ::  dqv_AT_dt          !  Specific           Humidity    TENDENCY, Atm_AT     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dqw_AT_dt          !  Cloud Droplets Concentration   TENDENCY, Atm_AT     [kg/kg/s]
! #cw real, dimension(kcolq,mzq)                            ::  dCw_AT_dt          !  CCN            Concentration   TENDENCY, Atm_AT         [1/s]
      real, dimension(kcolq,mzq)                            ::  dqi_AT_dt          !  Cloud Crystals Concentration   TENDENCY, Atm_AT     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dCi_AT_dt          !  CIN            Concentration   TENDENCY, Atm_AT         [1/s]
      real, dimension(kcolq,mzq)                            ::  dqs_AT_dt          !  Snow Particles Concentration   TENDENCY, Atm_AT     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dqr_AT_dt          !  Rain Drops     Concentration   TENDENCY, Atm_AT     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dpktCM_dt          !  Reduced  Potential Temperature Tendency, CMiPhy        [KX/s]
      real, dimension(kcolq,mzq)                            ::  dqv_CM_dt          !  Specific           Humidity    TENDENCY, CMiPhy     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dqw_CM_dt          !  Cloud Droplets Concentration   TENDENCY, CMiPhy     [kg/kg/s]
! #cw real, dimension(kcolq,mzq)                            ::  dCw_CM_dt          !  CCN            Concentration   TENDENCY, CMiPhy         [1/s]
      real, dimension(kcolq,mzq)                            ::  dqi_CM_dt          !  Cloud Crystals Concentration   TENDENCY, CMiPhy     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dCi_CM_dt          !  CIN            Concentration   TENDENCY, CMiPhy         [1/s]
      real, dimension(kcolq,mzq)                            ::  dCF_CM_dt          !  Cloud Fraction                 TENDENCY, CMiPhy         [1/s]
      real, dimension(kcolq,mzq)                            ::  dqs_CM_dt          !  Snow Particles Concentration   TENDENCY, CMiPhy     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dqr_CM_dt          !  Rain Drops     Concentration   TENDENCY, CMiPhy     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dpktCP_dt          !  Reduced  Potential Temperature TENDENCY, CVAmnh        [KX/s]
      real, dimension(kcolq,mzq)                            ::  dqv_CP_dt          !  Specific           Humidity    TENDENCY, CVAmnh     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dqw_CP_dt          !  Cloud Droplets Concentration   TENDENCY, CVAmnh     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dqi_CP_dt          !  Cloud Crystals Concentration   TENDENCY, CVAmnh     [kg/kg/s]
      real, dimension(kcolq,mzq)                            ::  dpktRT_dt          !  Reduced  Potential Temperature TENDENCY, radCEP        [KX/s]

!dead integer                                               ::  it0EXP             !
!dead integer                                               ::  it0RUN             !
      integer                                               ::  Year_H             !  Time                                                   [year]
      integer                                               ::  Mon__H             !  Time                                                  [month]
      integer                                               ::  Day__H             !  Time                                                    [Day]
      integer                                               ::  Hour_H             !  Time                                                   [hour]
      integer                                               ::  minu_H             !  Time                                                 [minute]
      integer                                               ::  sec__H             !  Time                                                      [s]
      integer                                               ::  ixq1,i0x0,mxqq     !  Domain  Dimension: x                                      [-]
      integer                                               ::  jyq1,j0y0,myqq     !  Domain  Dimension: y                                      [-]
      integer                                               ::  mzq                !  Domain  Dimension: z                                      [-]
      integer                                               ::  mzqq               !  Domain  Dimension: z                                      [-]
      integer                                               ::  mwq                !  Domain  Dimension: mosaic                                 [-]
      integer                                               ::  kcolq              !  Domain  Dimension: x * y                                  [-]
      integer                                               ::  kcolw              !  Domain  Dimension: x * y * mosaic                         [-]
      integer                                               ::  m_azim             !  Mountain Mask, nb of directions taken into account        [-]
      integer, dimension(n0pt)                              ::  IOi0SV             !
      integer, dimension(n0pt)                              ::  IOj0SV             !
      integer                                               ::  n0pt               !

      real, dimension(kcolq)                                ::  sst__HOST          !  Ocean     FORCING (SST)                                   [K]
! #IP real, dimension(kcolq)                                ::  sif__HOST          !  Ocean     FORCING (Sea-Ice Fraction )                     [-]
! #AO real, dimension(ixq1:mxqq,jyq1:myqq,mwq)              ::  s_T__HOST          !  A - O    COUPLING                      n=1: Open Ocean    [K]
! #AO real, dimension(ixq1:mxqq,jyq1:myqq,mwq)              ::  Alb__HOST          !  A - O    COUPLING (Surface Albedo   )  n=2: Sea  Ice      [-]
! #AO real, dimension(ixq1:mxqq,jyq1:myqq,mwq)              ::  dSdT2HOST          !  A - O    COUPLING ( d(SH Flux) / dT )                [W/m2/K]
! #AO real, dimension(ixq1:mxqq,jyq1:myqq,mwq)              ::  dLdT2HOST          !  A - O    COUPLING ( d(SH Flux) / dT )                [W/m2/K]

! #TC integer, parameter                                    ::   ntrac  =  28      !
! #TC real, dimension(ixq1:mxqq,jyq1:myqq,mzq,ntrac)        ::    qxTC             !  Aerosols: Atmospheric  Contentration
! #TC real, dimension(ixq1:mxqq,jyq1:myqq    ,ntrac)        ::    qsTC             !  Aerosols: Near Surface Contentration
! #TC real, dimension(ixq1:mxqq,jyq1:myqq    ,ntrac)        ::    uqTC             !  Aerosols: Surf.Flux
! ---------------------------------------------------------------------------------!




! LOCAL VARIABLES
! ===============

      real(kind=real8), dimension(kcolq,mzq)                 ::  Wind_HOST            !  Wind Speed, Horizontal                  [m/s]
      real(kind=real8)                                       ::  dTdz   = 0.0065      ! -d(T) / dz   Lapse Rate                  [K/m]
      real(kind=real8)                                       ::  dTimAT               !  d(Time) between 2 calls of Atmos.Turbul.  [s]

      integer                                                ::  i   ,j   ,ikl   ,ikp !
      integer                                                ::  k   ,mn              !
      integer                                                ::  n   ,ipt ,iwr   ,kk  !





!dead     it_RUN = it0RUN
!dead     it_EXP = it0EXP




! INITIALIZATION
! ==============

! Initialization of local Variables (used each Time Step)
! -------------------------------------------------------

      DO k = 1,mzq
      DO ikl=1,kcolq  
        Wind_HOST(ikl,k) = sqrt(Ua___HOST(ikl,k)*Ua___HOST(ikl,k)+Va___HOST(ikl,k)*Va___HOST(ikl,k))
      ENDDO
      ENDDO



! Initialization of the run
! -------------------------

! Martin Control
     PRINT*, 'Dans PHY_MAR:'
     PRINT*, 'it_RUN=',it_RUN
! Martin Control

      IF (it_RUN.LE.1)                                              THEN 



! Initialization of Mod_SISVAT_grd
! --------------------------------

              jt__SV =              max(1,int(dt0_SV /dt0DYn))
              dt__SV =                   real(jt__SV)*dt0DYn
          IF (dt__SV .NE. dt0_SV) write(6,61) dt__SV ,dt0_SV
 61         format    ('dt__SV =',f9.3,' differs from dt0_SV =',f9.3)


! Initialization of Mod_PHY_AT_grd
! --------------------------------

              jt__AT =              max(1,int(dt0_AT /dt0DYn))
              dt__AT =                   real(jt__AT)*dt0DYn
          IF (dt__AT .NE. dt0_AT) write(6,62) dt__AT ,dt0_AT
 62         format    ('dt__AT =',f9.3,' differs from dt0_AT =',f9.3)


! Initialization of Mod_PHY_CM_grd
! --------------------------------

              jt__CM =              max(1,int(dt0_CM /dt0DYn))
              dt__CM =                   real(jt__CM)*dt0DYn
          IF (dt__CM .NE. dt0_CM) write(6,63) dt__CM ,dt0_CM
 63         format    ('dt__CM =',f9.3,' differs from dt0_CM =',f9.3)


! Initialization of Mod_PHY_CP_grd
! --------------------------------

              jt__CP =              max(1,int(dt0_CP /dt0DYn))
              dt__CP =                   real(jt__CP)*dt0DYn
          IF (dt__CP .NE. dt0_CP) write(6,64) dt__CP ,dt0_CP
 64         format    ('dt__CP =',f9.3,' differs from dt0_CP =',f9.3)


! Initialization of Mod_PHY_RT_grd
! --------------------------------

              jt__RT =              max(1,int(dt0_RT /dt0DYn))
              dt__RT =                   real(jt__RT)*dt0DYn
          IF (dt__RT .NE. dt0_RT) write(6,65) dt__RT ,dt0_RT
 65         format    ('dt__RT =',f9.3,' differs from dt0_RT =',f9.3)


! Initialization
! --------------

! Initialization of 1-D Axes Variables: Vertical   Axis (Atmosphere)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              pt__DY     = ptop_HOST
          DO k=1,mzpp
              sigma(k)   = s_HOST(k)
          ENDDO

! Martin control
!PRINT*,'s_HOST=',s_HOST 
! Martin control

          DO k=1,mzp
                k1m(k)   = max(k-1, 1)
                k1p(k)   = min(k+1,mzp)
                k2m(k)   = max(k-2, 1)
             dsigma(k)   =  sigma(k+1) - sigma(k)
              sigmi(k+1) = (sigma(k+1) + sigma(k)) * 0.5
          END DO
              sigmi(1)   = 0.0
              sigmi(mzpp)= 1.0

          DO k=1,mzp
             dsigmi(k)   =  sigmi(k+1) - sigmi(k)

! Guess of sigma-levels Height
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             hsigma(k)   =-(288./.0065)*(sigma(k) ** (287.*.0065/9.81)-1.)
          END DO

             write(6,501) (sigma(k),k=1,mzpp)
  501          format(/,'  sigma: ',10f8.5                             &
     &             ,5(/,'         ',10f8.5))

             write(6,502)(hsigma(k),k=1,mzp)
  502          format(/,' hsigma: ',10f8.1                             &
     &             ,5(/,'         ',10f8.1))



! Initialization of Mod_PHY____dat & Mod_PHY____kkl
! -------------------------------------------------

! Topography
! ~~~~~~~~~~
             sh_MAX      = 0.
        DO ikl = 1,kcolp
             sh__AP(ikl) =            sh___HOST(ikl)
             sha_AP(ikl) =            sh_a_HOST(ikl)

! Surface Slope
! ~~~~~~~~~~~~~
             sloxAP(ikl) =            slopxHOST(ikl)
             sloyAP(ikl) =            slopyHOST(ikl)
             slopAP(ikl) =            slopeHOST(ikl)
             sh_MAX      = max(sh_MAX,sh__AP   (ikl))
        END DO
             dzaMIN      = hsigma(mzp) *(hsigma(1) -sh_MAX) /hsigma(1)


        DO ikl = 1,kcolp

! Geographic Coordinates
! ~~~~~~~~~~~~~~~~~~~~~~

! Martin rearrangement pour que RADACA ne plante plus:

!           lon__r(ikl)   =     lonh_HOST(ikl) * 2.0 * piNmbr / 24.0
           IF ((lonh_HOST(ikl) ) .LT. 0) THEN
             lon__r(ikl) = 360 + (lonh_HOST(ikl) * 2.0 * piNmbr / 24.0)
            ELSE
             lon__r(ikl) = (lonh_HOST(ikl) * 2.0 * piNmbr / 24.0)
           ENDIF
           lon__h(ikl)   =     lonh_HOST(ikl)
           lat__r(ikl)   =     latr_HOST(ikl)
           sinLat(ikl)   = sin(lat__r(ikl))
           cosLat(ikl)   = cos(lat__r(ikl))

        ENDDO



! -----------------------------------------------------------------------------!
! Initialization of  Atm_DY (Counterpart of dynamical variables in PHY_MAR)
! -------------------------

! Martin control
!PRINT*,'Avant PHY_Atm_DY_INI'
!PRINT*,'size(psa_DY)=',size(psa_DY)
! Martin control

!                            **************
                       CALL  PHY_Atm_DY_INI
!                            **************

! Initialization of  Atm_DY: needs plausible Atmospheric Conditions
! ~~~~~~~~~~~~~~~~~~~~~~~~~ (here  idealized for Temperature Vertical Gradient)
! Martin control
!PRINT*,'Apres PHY_Atm_DY_INI'
!PRINT*,'size(psa_DY)=',size(psa_DY)
!PRINT*,'Dans PHY_MAR, calcul de Ta__DY:'
!PRINT*,'mzpp=',mzpp
!PRINT*,'pt__DY=',pt__DY
!PRINT*,'Dtdz=',dTdz
!PRINT*,'RCp=',RCp
!PRINT*,'minval(pkta(:,mzpp))=',minval(pkta_HOST(:,mzpp))
!PRINT*,'minval(psa__HOST(:))=',minval(psa__HOST(:))
!PRINT*,'minval(Z___DY(:,:))=',minval(Z___DY(:,:))
!PRINT*,'minval(Ta__DY(:,mzpp))=',minval(Ta__DY(:,mzpp))
! Martin control

        DO ikl = 1,kcolp
           i                    = ii__AP   (ikl)
           j                    = jj__AP   (ikl)
           psa_DY    (ikl     ) = psa__HOST(ikl)
           Ta__DY    (ikl,mzpp) = pkta_HOST(ikl,mzpp)*(psa_DY(ikl)+pt__DY)**RCp
!gilles : init
           Z___DY    (ikl,mzpp) = 0.
           qv__DY    (ikl,mzpp) = 0.001
         DO k =   mzp,1,-1
           Z___DY    (ikl,k   ) =                                      &
     &                 Ta__DY(ikl,mzpp)                                &
     &      * (1.0 - ((psa_DY(ikl)*sigma (k)+pt__DY)                   &
     &               /(psa_DY(ikl)          +pt__DY))                  &
     &              **(R_DAir     * dTdz    /Grav_F))       / dTdz
           Ta__DY    (ikl,k) =                                         &
     &                 Ta__DY(ikl,mzpp)     -Z___DY(ikl,k)  * dTdz

           qv__DY    (ikl,k) =   0.001

           WindDY    (ikl,k) =   Wind_HOST(ikl,k)
           ua__DY    (ikl,k) =   Ua___HOST(ikl,k)
           va__DY    (ikl,k) =   va___HOST(ikl,k)
           wa__DY    (ikl,k) =   Wa___HOST(ikl,k)
         ENDDO
        ENDDO

! Martin Control
!PRINT*,'Avant PHY_Atm_S0_INI'
!PRINT*,'minval(Ta__DY(:,:))=',minval(Ta__DY(:,:))


! -----------------------------------------------------------------------------!
! Initialization of  Atm_S0 (Insolation and cos of Sun Zenithal Distance)
! -------------------------

!                            **************
                       CALL  PHY_Atm_S0_INI
!                            **************

! Initialization of Mod_PHY_S0_ctr and Mod_PHY_S0_kkl
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              FaceS0 = FlagS0_SLO
              MMskS0 = FlagS0_MtM
          IF (FlagS0_SLO .AND. FlagS0_MtM)                          THEN
              DO k   = 1,m_azim
              DO ikl = 1,kcolp
              cszkS0(ikl,k) = MMaskHOST(ikl,k)
              ENDDO
              ENDDO
          ENDIF


! -----------------------------------------------------------------------------!
! Initialization of  SISVAT
! -------------------------

! Initialization of Mod_SISVAT_ctr  (SISVAT Switches
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        & Time/Space control Variables)
        VegMod = FlagSV_Veg
        SnoMod = FlagSV_SNo
        BloMod = FlagSV_BSn
        InpSWD = FlagSV_SWD
        InpSBC = FlagSV_SBC
        SVaKzT = FlagSV_KzT
        SVaUBC = FlagSV_UBC


! Initialization of Mod_SISVAT_kkl  (from INPUT from HOST Model)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DO ikp = 1,kcolp
           i   = ii__AP(ikp)
           j   = jj__AP(ikp)
        DO mn  = 1,mwp
        DO k   = 1,mzp
           kk  =   mzpp - k
           WindSV(ikp,mn,k) = Wind_HOST(ikp,k)
           pkt0SV(ikp,mn,kk)= pkta_HOST(ikp,k)
        END DO
           Ua__SV(ikp,mn)   = Ua___HOST(ikp,mzp)
           Va__SV(ikp,mn)   = Va___HOST(ikp,mzp)
        END DO
        END DO


! Initialization of SISVAT Variables
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                            **************
                       CALL  PHY_SISVAT_INI
!                            **************


! Initialization of Mod_SISVAT_kkl  (SISVAT OUTPUT Grid Points)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           iwr = 0
        DO ipt = 1,NbPts
           IF (IOi0SV(ipt).EQ.0)                                    THEN
               IOi_SV(ipt)=i_x0
           ELSE
               IOi_SV(ipt)=IOi0SV(ipt)
           END IF
           IF (IOj0SV(ipt).EQ.0)                                    THEN
               IOj_SV(ipt)=j_y0
           ELSE
               IOj_SV(ipt)=IOj0SV(ipt)
           END IF
           write(6,*) 'ipt   , IOi_SV, IOj_SV              = '               &
     &                ,ipt   , IOi_SV(ipt), IOj_SV(ipt)
        DO n   = 1,mwp
           iwr = 1+iwr
          IF (iwr.LE.nbwri)                                         THEN
           no__SV(iwr) = 0
           i___SV(iwr) = IOi_SV(ipt)
           j___SV(iwr) = IOj_SV(ipt)
           n___SV(iwr) = n
           write(6,*) 'n     , i___SV, j___SV, n___SV, iwr = '               &
     &                ,n     , i___SV(iwr), j___SV(iwr), n___SV(iwr), iwr
          END IF
        END DO
        END DO



! -----------------------------------------------------------------------------!
! Initialization of  Atm_RT (Radiative Transfert through the Atmosphere)
! -------------------------

!gilles: PHY_Atm_RT_INI requires correct dates
          YearTU = Year_H
          Mon_TU = Mon__H
          Day_TU = Day__H
          HourTU = Hour_H
          minuTU = minu_H
          sec_TU = sec__H

!                            **************
                       CALL  PHY_Atm_RT_INI
!                            **************

! -----------------------------------------------------------------------------!
! Initialization of  Atm_AT (Turbulent Transfert through the Atmosphere)
! -------------------------

!                            **************
                       CALL  PHY_Atm_AT_INI(FlagAT_TKE,TypeAT)
!                            **************



! -----------------------------------------------------------------------------!
! Initialization of  Atm_CP (Convectiv Transfert through the Atmosphere)
! -------------------------

!                            **************
                       CALL  PHY_Atm_CP_INI(mzp,kcolp)
!                            **************



! -----------------------------------------------------------------------------!
! Initialization of  CMiPhy (Cloud Microphysical Scheme)
! -------------------------

! Initialization of Mod_PHY_CM_ctr  (CMiPhy Switches
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        & Time/Space control Variables)
        CM_UpD = FlagCM_UpD

!                            **************
                       CALL  PHY_Atm_CM_INI
!                            **************

      END IF


! Martin CONTROL
!PRINT*, 'Impressions control'
!call iophys_ecrit('TA__DY_surf',1,'surface temperature','K',Ta__DY(:,61))
!call iophys_ecrit('TA__DY_air',60,'air temperature','K',Ta__DY(:,1:60))

!PRINT*,'TA__DY(:,61)=',TA__DY(kcolq/2,61)

! Interface: From HOST Model Variables to MARp Physics Variables
! ==============================================================


! Time
! ----

          YearTU = Year_H
          Mon_TU = Mon__H
          Day_TU = Day__H
          HourTU = Hour_H
          minuTU = minu_H
          sec_TU = sec__H

          TimeTU = (float(351)+(float(YearTU) -float(1902)) *float(365)&! Nb Days before YearTU
     &                        +(float(YearTU) -float(1901)) /float(  4)&! Nb Leap Years
     &                 + float(njYear(Mon_TU))                         &! Nb Days before Mon_TU
     &                 + float(njLeap(Mon_TU))                         &! (including Leap Day)
     &       *max(zer0,un_1-mod(float(YearTU),float(4)))               &!
     &                        + float(Day_TU)-float(1)  )   *float( 24)&!
     &             +float(HourTU)                                      &!
     &           + (float(minuTU) *float(60) +float(sec_TU))/3600.      !



! -----------------------------------------------------------------------------!
! Assignation    of  Mod_PHY_AT_grd
! ---------------------------------
          IF (it_RUN .EQ. 1) THEN
            TimeAT = TimeTU-dt0_AT/3600. ! Initialisation à la première itération
          END IF
          IF (FlagAT  .AND.      mod(it_RUN-1,jt__AT).EQ.0)         THEN
!gilles: leapfrog scheme stops time integration regularly
!              dTimAT =              (TimeTU-TimeAT) * 3600.
              dTimAT =    dt0_AT
              TimeAT =               TimeTU
          END IF

! Martin CONTROL
PRINT*,'CONTROL PHY_MAR temps'
PRINT*,'  Year_H =', Year_H
PRINT*,'  Mon__H =', Mon__H
PRINT*,'  Day__H =', Day__H
PRINT*,'  Hour_H =', Hour_H
PRINT*,'  minu_H =', minu_H
PRINT*,'  sec__H =', sec__H
PRINT*,'jt__AT=',jt__AT
PRINT*,'TimeTU=',TimeTU
PRINT*,'TimeAT=',TimeAT
PRINT*,'dTimAT=',dTimAT
! Martin CONTROL


! -----------------------------------------------------------------------------!
! Assignation    of  Mod_PHY_DY_kkl
! ---------------------------------

          DO ikl=1,kcolp
             i  =  ii__AP(ikl)
             j  =  jj__AP(ikl)
             k  =  mzpp
            ExnrDY    (ikl,k) = exp(RCp *log(psa__HOST(ikl)*sigma(k)+pt__DY))
            pkt_DY    (ikl,k) =              pkta_HOST(ikl,k)
            Ta__DY    (ikl,k) =              pkta_HOST(ikl,k) * ExnrDY(ikl,k)
!  CAUTION: Tas_SV_xy is not allowed to be changed by data coming from outside this routine
            Z___DY    (ikl,k) =              gZa__HOST(ikl,k)         *Grav_I
            ZmidDY    (ikl,k) =              gZam_HOST(ikl,k)         *Grav_I
            qv__DY    (ikl,k) =              qv___HOST(ikl,k)

          DO  k               =   1,mzp
            ExnrDY    (ikl,k) = exp(RCp *log(psa__HOST(ikl)*sigma(k  )+pt__DY))
            pkt_DY    (ikl,k) =              pkta_HOST(ikl,k)
            Ta__DY    (ikl,k) =              pkta_HOST(ikl,k) * ExnrDY(ikl,k)
            roa_DY    (ikl,k) =             (psa__HOST(ikl)*sigma(k  )+pt__DY)     &
     &                                     /(Ta__DY   (ikl,k)         *R_DAir)
            roamDY    (ikl,k) =             (psa__HOST(ikl)*sigmi(k+1)+pt__DY)     &
     &                                     /(Ta__DY   (ikl,k)         *R_DAir)
            Z___DY    (ikl,k) =              gZa__HOST(ikl,k)         *Grav_I
            ZmidDY    (ikl,k) =              gZam_HOST(ikl,k)         *Grav_I
            qv__DY    (ikl,k) =              qv___HOST(ikl,k)



! -----------------------------------------------------------------------------!
! Assignation    of  Mod_PHY_AT_kkl
! ---------------------------------


          IF (FlagAT.AND.it_EXP.GT.1.AND.mod(it_RUN-1,jt__AT).EQ.0) THEN
            TKE_AT    (ikl,k) =              TKE__HOST(ikl,k)
            eps_AT    (ikl,k) =              eps__HOST(ikl,k)
            TrT_AT    (ikl,k) =            (TKE_AT    (ikl,k) - TrT_AT    (ikl,k)) &
     &                                    / dTimAT
          END IF



! -----------------------------------------------------------------------------!
! Assignation    of  Mod_PHY_CM_kkl
! ---------------------------------

          IF (FlagCM.AND.it_EXP.GT.1.AND.mod(it_RUN-1,jt__CM).EQ.0) THEN

!           IF (qw___HOST(ikl,k).LT.qh_MIN)                         THEN
!               qv__DY(ikl,k) =             qv__DY(ikl,k) + qw__CM(ikl,k)
!               qw__CM(ikl,k) =             0.
! #cw           CCNwCM(ikl,k) =             0.
!           ELSE
                qw__CM(ikl,k) =             qw___HOST(ikl,k)
! #cw           CCNwCM(ikl,k) =             CCN__HOST(ikl,k)
!           END IF

!           IF (qi___HOST(ikl,k).LT.qh_MIN)                         THEN
!               qv__DY(ikl,k) =             qv__DY(ikl,k) + qi__CM(ikl,k)
!               qi__CM(ikl,k) =             0.
!               CCNiCM(ikl,k) =             0.
!           ELSE
                qi__CM(ikl,k) =             qi___HOST(ikl,k)
                CCNiCM(ikl,k) =             CIN__HOST(ikl,k)
!           END IF

! Gilles: CF___HOST non sauve & CFraCM reinitialise ici
!               CFraCM(ikl,k) =             CF___HOST(ikl,k)

!           IF (qw__CM(ikl,k).LT.qh_MIN  .AND.                         &
!    &          qi__CM(ikl,k).LT.qh_MIN)                            THEN
!               CFraCM(ikl,k) =             0.
!           ELSE
!               CFraCM(ikl,k) =  max(CFrMIN,CF___HOST(ikl,k))
!           END IF

                qs__CM(ikl,k) =             qs___HOST(ikl,k)
                qr__CM(ikl,k) =             qr___HOST(ikl,k)
          END IF

          ENDDO



! -----------------------------------------------------------------------------!
! Assignation    of  Mod_SISVAT_gpt                  (A-O FORCING OR COUPLING) 
! ---------------------------------                 

            sst_SB    (ikl)   =              sst__HOST(ikl)  
! #IP       sif_SB    (ikl)   =              sif__HOST(ikl)  
! #AO     DO k = 1,mwp
! #AO       s_T_AO_xyn(i,j,k) =              s_T__HOST(ikl,k)
! #AO       Alb_AO_xyn(i,j,k) =              Alb__HOST(ikl,k)
! #AO     ENDDO

          ENDDO

! -----------------------------------------------------------------------------!
! Assignation    of  Mod_PHY_DY_kkl
! ---------------------------------
          DO ikl = 1,kcolp
             i   = ii__AP(ikl)
             j   = jj__AP(ikl)
             psa_DY(ikl)   = psa__HOST(ikl)
          DO k   = 1,mzp
             WindDY(ikl,k) = Wind_HOST(ikl,k)
             ua__DY(ikl,k) = Ua___HOST(ikl,k)
             va__DY(ikl,k) = va___HOST(ikl,k)
             wa__DY(ikl,k) = Wa___HOST(ikl,k)  +  sqrt(2.*max(eps6,TKE_AT(ikl,k))/3.)
          END DO
          END DO

!                            **************
                       CALL  PHY_Atm_DY_RUN                                    ! Assignation of MAR Dyn. Variables
!                            **************

           IF(FlagVR)                                               THEN
!                            **************
                       CALL  PHY________OUT('After  PHY_Atm_DY_RUN                             ')
!                            **************  12345678901234567890123456789012345678901234567890
!                                                     1         2         3         4         5
           END IF



! -----------------------------------------------------------------------------!
! Saturation Specific Humidity
! ----------------------------

          IF (FlagCM .OR. &! ***************
     &        FlagSV)  CALL  PHY_Atm_CM_QSat
!                            ***************



! -----------------------------------------------------------------------------!
! Execution      of  Atm_S0 (Insolation and cos of Sun Zenithal Distance)
! -------------------------

!                            **************
                       CALL  PHY_Atm_S0_RUN
!                            **************

           IF(FlagVR)                                               THEN
!                            **************
                       CALL  PHY________OUT('After  PHY_Atm_S0_RUN                             ')
!                            **************  12345678901234567890123456789012345678901234567890
!                                                     1         2         3         4         5
           END IF



! -----------------------------------------------------------------------------!
! Execution      of  Atm_RT (Radiative Transfer through the Atmosphere)
! -------------------------

          IF (FlagRT  .AND.      mod(it_RUN-1,jt__RT).EQ.0)         THEN

!                            **************
                       CALL  PHY_Atm_RT_RUN(kcolp,ikl0)
!                            **************

           IF(FlagVR)                                               THEN
!                            **************
                       CALL  PHY________OUT('After  PHY_Atm_RT_RUN                             ')
!                            **************  12345678901234567890123456789012345678901234567890
!                                                     1         2         3         4         5
           END IF
          END IF



! -----------------------------------------------------------------------------!
! Execution      of  SISVAT (Soil-Ice-Snow-Vegetation-Atmosphere-Transfer Scheme)
! -------------------------

          IF (FlagSV)                                               THEN

! Assignation    of Mod_SISVAT_kkl  (from INPUT from HOST Model)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            DO ikp = 1,kcolp
            DO mn  = 1,mwp
               pkPaSV(ikp,mn)   =              psa__HOST(ikp)           +pt__DY
               ExnrSV(ikp,mn)   = exp(RCp *log(psa__HOST(ikp)*sigma(mzp)+pt__DY))
            DO k   = 1,mzp
               kk  =   mzpp - k
               WindSV(ikp,mn,k) =              Wind_HOST(ikp,k)
               pkt0SV(ikp,mn,kk)=              pkta_HOST(ikp,k)
               qv__SV(ikp,mn,kk)=              qv___HOST(ikp,k)
               zza_SV(ikp,mn,kk)=             (gZa__HOST(ikp,k)                &  !
      &                                       -gZa__HOST(ikp,mzpp))     *Grav_I   !
               roa_SV(ikp,mn,kk)=              roa_DY   (ikp,k)         * 1.e3    ! [kg/m3]
            END DO
               Ua__SV(ikp,mn)   =              Ua___HOST(ikp,mzp)
               Va__SV(ikp,mn)   =              Va___HOST(ikp,mzp)
               TaT_SV(ikp,mn)   =              pkta_HOST(ikp,mzp) * ExnrSV(ikp,mn)
            END DO
            END DO


!                            **************
                       CALL  PHY_SISVAT_RUN(                           &
!                            **************
! -----------------------------------------------------------------------------!
! #TC&                                        qxTC,uqTC  ,qsTC  ,      &       ! Aerosols: Atm.Conc., Surf.Flux
! -----------------------------------------------------------------------------!
     &                                     )

           IF(FlagVR)                                               THEN
!                            **************
                       CALL  PHY________OUT('After  PHY_SISVAT_RUN                             ')
!                            **************  12345678901234567890123456789012345678901234567890
!                                                     1         2         3         4         5
           END IF
          END IF



! -----------------------------------------------------------------------------!
! Execution      of  Atm_AT (Atmospheric Turbulence Contribution)
! -------------------------

          IF (FlagAT  .AND.      mod(it_RUN-1,jt__AT).EQ.0)         THEN

!                            **************
                       CALL  PHY_Atm_AT_RUN(FlagSV_KzT,FlagCM)
!                            **************

           IF(FlagVR)                                               THEN
!                            **************
                       CALL  PHY________OUT('After  PHY_Atm_AT_RUN                             ')
!                            **************  12345678901234567890123456789012345678901234567890
!                                                     1         2         3         4         5
           END IF
          END IF



! -----------------------------------------------------------------------------!
! Execution      of  Atm_CP (Convection Parameterisation - Mass Flux Scheme)
! -------------------------

          IF (FlagCP  .AND.      mod(it_RUN-1,jt__CP).EQ.0)         THEN

!                            **************
                       CALL  PHY_Atm_CP_RUN(mzp,kcolp)
                       ! CALL  PHY_Atm_CP_RUN(mxp,kcolp) ! ancien bug
!                            **************

           IF(FlagVR)                                               THEN
!                            **************
                       CALL  PHY________OUT('After  PHY_Atm_CP_RUN                             ')
!                            **************  12345678901234567890123456789012345678901234567890
!                                                     1         2         3         4         5
           END IF
          END IF



! -----------------------------------------------------------------------------!
! Execution      of  Atm_CM (Cloud Microphysics)
! -------------------------

          IF (FlagCM  .AND.      mod(it_RUN-1,jt__CM).EQ.0)         THEN

!                            **************
                       CALL  PHY_Atm_CM_RUN
!                            **************


           IF(FlagVR)                                               THEN
!                            **************
                       CALL  PHY________OUT('After  PHY_Atm_CM_RUN                             ')
!                            **************  12345678901234567890123456789012345678901234567890
!                                                     1         2         3         4         5
           END IF
          END IF




! Assignation of the Tendencies to transfer outside the physical Parameterizations Package
! ========================================================================================

          DO ikl=1,kcolp
             i   = ii__AP(ikl)
             j   = jj__AP(ikl)



! -----------------------------------------------------------------------------!
! Reinitialization of the Tendencies
! ----------------------------------

          DO k=   1,mzp

              dpkt___dt(ikl,k) = 0.0                                    !
              dua____dt(ikl,k) = 0.0                                    !
              dva____dt(ikl,k) = 0.0                                    !
              dqv____dt(ikl,k) = 0.0                                    !
              dqw____dt(ikl,k) = 0.0                                    !
              dCF____dt(ikl,k) = 0.0                                    !
              dqi____dt(ikl,k) = 0.0                                    !
              dCi____dt(ikl,k) = 0.0                                    !
              dqs____dt(ikl,k) = 0.0                                    !
              dqr____dt(ikl,k) = 0.0                                    !



! -----------------------------------------------------------------------------!
! Tendencies from SISVAT
! ----------------------

          IF (FlagSV  .AND.      mod(it_RUN-1,jt__SV).EQ.0)         THEN!
              dpktSV_dt(ikl,k) =                    dpktSV_gpt(ikl,k)   ! Reduced Potential Temperature TENDENCY, SISVAT        [KX/s]
          END IF

! Update           of the tendencies
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (FlagSV)                                               THEN!
              dpkt___dt(ikl,k) = dpkt___dt(ikl,k) + dpktSV_gpt(ikl,k)   !
          END IF



! -----------------------------------------------------------------------------!
! Tendencies from Atm_AT
! ----------------------

          IF (FlagAT  .AND.      mod(it_RUN-1,jt__AT).EQ.0)         THEN!
              dua____dt(ikl,k) =                    dua_AT(ikl,k)       ! Wind Speed       (x-direc.)   TENDENCY, ALL Contribut.[m/s2]
              dva____dt(ikl,k) =                    dva_AT(ikl,k)       ! Wind Speed       (y-direc.)   TENDENCY, ALL Contribut.[m/s2]
              dpktAT_dt(ikl,k) =                    dpktAT(ikl,k)       ! Reduced Potential Temperature TENDENCY, Atm_AT        [KX/s]
              dqv_AT_dt(ikl,k) =                    dqv_AT(ikl,k)       ! Specific          Humidity    TENDENCY, Atm_AT     [kg/kg/s]
              dqw_AT_dt(ikl,k) =                    dqw_AT(ikl,k)       ! Cloud Droplets Concentration  TENDENCY, Atm_AT     [kg/kg/s]
              dqi_AT_dt(ikl,k) =                    dqi_AT(ikl,k)       ! Cloud Crystals Concentration  TENDENCY, Atm_AT     [kg/kg/s]
              dqs_AT_dt(ikl,k) =                    dqs_AT(ikl,k)       ! Snow Particles Concentration  TENDENCY, Atm_AT     [kg/kg/s]
              dqr_AT_dt(ikl,k) =                    dqr_AT(ikl,k)       ! Rain Drops     Concentration  TENDENCY, Atm_AT     [kg/kg/s]
! #cw         dCw_AT_dt(ikl,k) =                    dCW_AT(ikl,k)       ! CCN            Concentration  TENDENCY, Atm_AT         [1/s]
              dCi_AT_dt(ikl,k) =                    dCi_AT(ikl,k)       ! CIN            Concentration  TENDENCY, Atm_AT         [1/s]
 
! NO Tendencies on TKE and its dissipation ==> TKE__HOST and eps__HOST are updated
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              TKE__HOST(ikl,k) =                    TKE_AT(ikl,k)       !
              eps__HOST(ikl,k) =                    eps_AT(ikl,k)       !
              
! Re-Initialization of the TKE Transport Rate
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              TrT_AT   (ikl,k) =                    TKE_AT(ikl,k)       !
          END IF

! Update           of the tendencies
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (FlagAT)                                               THEN!
              dpkt___dt(ikl,k) = dpkt___dt(ikl,k) + dpktAT(ikl,k)       !
              dqv____dt(ikl,k) = dqv____dt(ikl,k) + dqv_AT(ikl,k)       !
              dqw____dt(ikl,k) = dqw____dt(ikl,k) + dqw_AT(ikl,k)       !
              dqi____dt(ikl,k) = dqi____dt(ikl,k) + dqi_AT(ikl,k)       !
              dqs____dt(ikl,k) = dqs____dt(ikl,k) + dqs_AT(ikl,k)       !
              dqr____dt(ikl,k) = dqr____dt(ikl,k) + dqr_AT(ikl,k)       !
! #cw         dCw____dt(ikl,k) = dCw____dt(ikl,k) + dCw_AT(ikl,k)       !
              dCi____dt(ikl,k) = dCi____dt(ikl,k) + dCi_AT(ikl,k)       !
          END IF



! -----------------------------------------------------------------------------!
! Tendencies from CMiPhy (dqw_CM, dCw_CM, dqi_CM, dCi_CM
! ----------------------  dqs_CM, dqr_CM  are consumed in CMiPhy then reset to 0
!                                                      if FlagCM_UpD  = .TRUE. )

          IF (FlagCM  .AND.      mod(it_RUN-1,jt__CM).EQ.0)         THEN!
              dpktCM_dt(ikl,k) =                    dpktCM(ikl,k)       ! Reduced Potential Temperature TENDENCY, CMiPhy        [KX/s]
              dqv_CM_dt(ikl,k) =                    dqv_CM(ikl,k)       ! Specific          Humidity    TENDENCY, CMiPhy     [kg/kg/s]
              dqw_CM_dt(ikl,k) =                    dqw_CM(ikl,k)       ! Cloud Droplets Concentration  TENDENCY, CMiPhy     [kg/kg/s]
              dCF_CM_dt(ikl,k) =                    dCF_CM(ikl,k)       ! Cloud Fraction                TENDENCY, CMiPhy     [kg/kg/s]
              dqi_CM_dt(ikl,k) =                    dqi_CM(ikl,k)       ! Cloud Crystals Concentration  TENDENCY, CMiPhy     [kg/kg/s]
              dqs_CM_dt(ikl,k) =                    dqs_CM(ikl,k)       ! Snow Particles Concentration  TENDENCY, CMiPhy     [kg/kg/s]
              dqr_CM_dt(ikl,k) =                    dqr_CM(ikl,k)       ! Rain Drops     Concentration  TENDENCY, CMiPhy     [kg/kg/s]
! #cw         dCw_CM_dt(ikl,k) =                    dCW_CM(ikl,k)       ! CCN            Concentration  TENDENCY, CMiPhy         [1/s]
              dCi_CM_dt(ikl,k) =                    dCi_CM(ikl,k)       ! CIN            Concentration  TENDENCY, CMiPhy         [1/s]
          END IF

! Update           of the tendencies
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (FlagCM)                                               THEN!
              dpkt___dt(ikl,k) = dpkt___dt(ikl,k) + dpktCM(ikl,k)       !
              dqv____dt(ikl,k) = dqv____dt(ikl,k) + dqv_CM(ikl,k)       !
              dqw____dt(ikl,k) = dqw____dt(ikl,k) + dqw_CM(ikl,k)       !
              dCF____dt(ikl,k) = dCF____dt(ikl,k) + dCF_CM(ikl,k)       !
              dqi____dt(ikl,k) = dqi____dt(ikl,k) + dqi_CM(ikl,k)       !
              dqs____dt(ikl,k) = dqs____dt(ikl,k) + dqs_CM(ikl,k)       !
              dqr____dt(ikl,k) = dqr____dt(ikl,k) + dqr_CM(ikl,k)       !
! #cw         dCw____dt(ikl,k) = dCw____dt(ikl,k) + dCw_CM(ikl,k)       !
              dCi____dt(ikl,k) = dCi____dt(ikl,k) + dCi_CM(ikl,k)       !
          END IF



! -----------------------------------------------------------------------------!
! Tendencies from CVAmnh
! ----------------------

          IF (FlagCP  .AND.      mod(it_RUN-1,jt__CP).EQ.0)         THEN!
              dpktCP_dt(ikl,k) =                    dpktCP(ikl,k)       ! Reduced Potential Temperature TENDENCY, CVAmnh        [KX/s]
              dqv_CP_dt(ikl,k) =                    dqv_CP(ikl,k)       ! Specific          Humidity    TENDENCY, CVAmnh     [kg/kg/s]
              dqw_CP_dt(ikl,k) =                    dqw_CP(ikl,k)       ! Cloud Droplets Concentration  TENDENCY, CVAmnh     [kg/kg/s]
              dqi_CP_dt(ikl,k) =                    dqi_CP(ikl,k)       ! Cloud Crystals Concentration  TENDENCY, CVAmnh     [kg/kg/s]
          END IF

! Update           of the tendencies
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (FlagCP)                                               THEN!
              dpkt___dt(ikl,k) = dpkt___dt(ikl,k) + dpktCP(ikl,k)       !
              dqv____dt(ikl,k) = dqv____dt(ikl,k) + dqv_CP(ikl,k)       !
              dqw____dt(ikl,k) = dqw____dt(ikl,k) + dqw_CP(ikl,k)       !
              dqi____dt(ikl,k) = dqi____dt(ikl,k) + dqi_CP(ikl,k)       !
          END IF



! -----------------------------------------------------------------------------!
! Tendencies from radCEP
! ----------------------

          IF (FlagRT  .AND.      mod(it_RUN-1,jt__RT).EQ.0)         THEN!
              dpktRT_dt(ikl,k) =                    dpktRT(ikl,k)       ! Reduced Potential Temperature TENDENCY, radCEP        [KX/s]
          END IF

! Update           of the tendencies
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (FlagRT)                                               THEN!
              dpkt___dt(ikl,k) = dpkt___dt(ikl,k) + dpktRT(ikl,k)       !
          END IF

          ENDDO





! Assignation of the Variables  to keep      inside the Physical Parameterizations Package
! ========================================================================================

!         ...  NOTHING up to now



! Assignation of the Variables  to transfer outside the Physical Parameterizations Package
! ========================================================================================

! -----------------------------------------------------------------------------!
! Update of  pkta_HOST, qv___HOST 
! -------------------------------

          DO k=1,mzpp
            pkta_HOST(ikl,k) = pkt_DY(ikl,k)                                   !  Always on k=mzpp, possible from 1 to mzp (with dpkt = 0.)
            qv___HOST(ikl,k) = qv__DY(ikl,k)                                   !
          ENDDO



! -----------------------------------------------------------------------------!
! Update of  qw___HOST, CF___HOST, qi___HOST, CIN__HOST, qs___HOST, qr___HOST
! ---------------------------------------------------------------------------

      IF (FlagCM)                                                   THEN
          DO k=1,mzp
            qw___HOST(ikl,k) = qw__CM(ikl,k)
! #cw       CCN__HOST(ikl,k) = CCNwCM(ikl,k)
            CF___HOST(ikl,k) = CFraCM(ikl,k)
            qi___HOST(ikl,k) = qi__CM(ikl,k)
            CIN__HOST(ikl,k) = CCNiCM(ikl,k)
            qs___HOST(ikl,k) = qs__CM(ikl,k)
            qr___HOST(ikl,k) = qr__CM(ikl,k)
          ENDDO
      END IF



! -----------------------------------------------------------------------------!
! Update of  d(S,LH)/dT (needed in NEMO)
! ---------------------

! #AO     DO k=1,mwpp
! #AO       dSdT2HOST(ikl,k) = dSdTAO_xyn(i,j,k)
! #AO       dLdT2HOST(ikl,k) = dLdTAO_xyn(i,j,k)
! #AO     ENDDO



          ENDDO




! OUTPUT
! ======

! OUTPUT of Tendencies
! --------------------

          IF (FlagVR .OR.                                              &
     &       (FLAG_O .AND. ((minuTU.EQ.0 .AND. sec_TU.EQ.0) .OR.       &
     &                       it_RUN.EQ.1                        ))) THEN
             ikl = ikl0
             i   = ii__AP(ikl)
             j   = jj__AP(ikl)

! pkt TENDENCIES
! ~~~~~~~~~~~~~~
                write(4,400)
 400            format(//,'   pkt TENDENCIES',/,'   **************'/,1x)
                write(4,403) Day_TU,Mon_TU,YearTU,HourTU,MinuTU,Sec_TU,it_EXP
 403            format(3x,2(i2,'-'),i4,4x,3(i2,'-'),'  Simulation Iteration No ',i6,/,1x)
                write(4,404)
                write(4,401)
 401            format('    |   SISVAT    |   Atm_AT    |   CMiPhy    |   CVAmnh    |   radCEP    |')
                write(4,402)
 402            format('    |      [K/d]  |      [K/d]  |      [K/d]  |      [K/d]  |      [K/d]  |')
                write(4,404)
 404            format(4('-'),'+',5(13('-'),'+'))
             DO k=1,mzp
                write(4,405)        k                                  &
     &                     , ExnrDY(ikl,k)*dpktSV_gpt(ikl,k)*86400.    &
     &                     , ExnrDY(ikl,k)*dpktAT    (ikl,k)*86400.    &
     &                     , ExnrDY(ikl,k)*dpktCM    (ikl,k)*86400.    &
     &                     , ExnrDY(ikl,k)*dpktCP    (ikl,k)*86400.    &
     &                     , ExnrDY(ikl,k)*dpktRT    (ikl,k)*86400.
 405            format(i3,' |',5(f12.6,' |'))

                IF (mod(k,20).EQ.0)                                 THEN
                write(4,404)
                write(4,401)
                write(4,402)
                write(4,404)
                END IF
             ENDDO
                write(4,404)

! pkt TENDENCIES
! ~~~~~~~~~~~~~~
                write(4,410)
 410            format(//,'   Qv  TENDENCIES',/,'   **************'/,1x)
                write(4,403) Day_TU,Mon_TU,YearTU,HourTU,MinuTU,Sec_TU,it_EXP
!403            format(3x,2(i2,'-'),i4,4x,3(i2,'-'),'  Simulation Iteration No ',i6,/,1x)
                write(4,404)
                write(4,411)
 411            format('    |   SISVAT    |   Atm_AT    |   CMiPhy    |   CVAmnh    |    TOTAL    |')
                write(4,412)
 412            format('    | [g/kg/min]  | [g/kg/min]  | [g/kg/min]  | [g/kg/min]  | [g/kg/min]  |')
                write(4,404)
!404            format(4('-'),'+',5(13('-'),'+'))
             DO k=1,mzp
                write(4,415)        k                                  &
     &                     ,               dqv_AT   (ikl,k)*60000.     &
     &                     ,               dqv_CM   (ikl,k)*60000.     &
     &                     ,               dqv_CP   (ikl,k)*60000.     &
     &                     ,               dqv____dt(ikl,k)*60000.
 415            format(i3,' |',12x,' |',4(f12.6,' |'))

                IF (mod(k,20).EQ.0)                                 THEN
                write(4,404)
                write(4,411)
                write(4,412)
                write(4,404)
                END IF
             ENDDO
                write(4,404)
          END IF

          IF (FLAG_O .AND. ((minuTU.EQ.0 .AND. sec_TU.EQ.0) .OR.       &
     &                       it_RUN.EQ.1                        ))  THEN
              
!                            **************
                       CALL  PHY________OUT('After  PHY_MAR                                    ')
!                            **************  12345678901234567890123456789012345678901234567890
!                                                     1         2         3         4         5
          END IF  



      return
      end
