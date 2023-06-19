      module Mod_SISVAT_kkl

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVAT_kkl contains the main (prognostic) variables of    |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  4-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE


! SISVAT INPUT        Variables
! -----------------------------

      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  LSmask  ! Land-Sea   Mask
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  isotSV  ! Soil       Type
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  iWaFSV  ! Soil       Drainage:(1,0)=(y,n)
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  ivgtSV  ! Vegetation Type

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  FracSV  ! Grid Cell Fraction  (Mosaic)                       [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  coszSV  ! Cosine of Sun zenithal Angle
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  sol_SV  ! Downward  Solar    Radiation
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  IRd_SV  ! Downward  Longwave Radiation

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  drr_SV  ! Rain  Intensity                              [kg/m2/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dsn_SV  ! Snow  Intensity                              [kg/m2/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dsnbSV  ! Idem, fraction, from Drift                         [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  esnbSV  ! Idem, fraction, from Drift                         [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dbs_SV  ! Drift Amount                                   [kg/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  BrosSV  ! Buffer Snow Layer Density
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  BG1sSV  ! Buffer Snow Layer Dendricity / Sphericity          [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  BG2sSV  ! Buffer Snow Layer Sphericity / Size                [-] [0.0001 m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dz0_SV  ! dz0(Sastrugi dh)                                   [m]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  cld_SV  ! Cloudiness (seen from SBL)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  za__SV  ! SBL Height
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Ua__SV  !(SBL Top)  Wind Velocity, x-Direction, t          [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Ua0_SV  !(SBL Top)  Wind Velocity, x-Direction, t -dt      [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Va__SV  !(SBL Top)  Wind Velocity, y-Direction, t          [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Va0_SV  !(SBL Top)  Wind Velocity, y-Direction, t -dt      [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  VV__SV  !(SBL Top)  Wind Velocity
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  VV10SV  ! 10-m      Wind Velocity
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  VVs_SV  !(Sastr,V)  Relevance
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  RRsxSV  !(Sastr,V)  Counter
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  DDsxSV  !(Sastr,V)  Angle
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  DDs_SV  !(Sastr,V)  Angle
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  TaT_SV  ! SBL Top   Temperature
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Ts__SV  ! Surface   Air Temperature     (Mosaic)             [K]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  pkPaSV  ! Surface   Pressure                               [kPa]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  WindSV  ! Wind      Speed                                  [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  zza_SV  ! Atmospheric  Levels HEIGHTS                        [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  roa_SV  ! Air       Volumic   Mass                        [T/m3]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  Kz__SV  ! Turbulent Diffusion Coefficients                [m2/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  pktaSV  ! Temperature / Exner Potential (Current Value in (PHY_)SISVAT)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  pkt0SV  ! Temperature / Exner Potential (INPUT         of (PHY_)SISVAT)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  ExnrSV  ! Surface       Exner Potential
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  qv__SV  ! Specific  Humidity                             [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  rhT_SV  ! SBL Top   Air  Density
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  QaT_SV  ! SBL Top   Specific Humidity
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dQa_SV  ! SBL Flux  Limitation of Qa
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  SHumSV  ! Surface   Specific Humidity
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dSdTSV  ! Sensible Heat Flux T Derivat.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dLdTSV  ! Latent   Heat Flux T Derivat.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  qsnoSV  ! SBL Mean  Snow       Content

      real(kind=real8), SAVE                              ::  zSBLSV  ! SBL Height (Initial Value)

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  z0__SV  ! Roughness Length Momentum    (Mosaic)              [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  LAI0SV  ! Nominal Leaf Area Index
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  glf0SV  ! Green   Leaf Fraction

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  alb0SV  ! Soil    Albedo
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  slopSV  ! Snow/Ice/Soil-Water Surf. Slope                    [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  slorSV  ! Snow/Ice/Soil-Water Surf. Slope               [radian]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  ROF_SV  !  Cumulative Run-Off          (Mosaic)        [mm w.e.]

      character(len=18)                             ::  daHost  ! Date Host Model


! SISVAT INPUT/OUTPUT Variables
! -----------------------------

      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  isnoSV  ! Nb of Ice/Snow Layers
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  ispiSV  ! Uppermost superimposed ice
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  iiceSV  ! Nb of Ice      Layers
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:,:)  ::  istoSV  ! Snow Layer     History

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  albcSV  ! Coupl. Surface Albedo (Surface-Canopy / Ocean)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  alb_SV  ! Surface-Canopy Albedo
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  emi_SV  ! Surface-Canopy Emissivity
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  IRs_SV  ! Soil           IR Flux
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  LMO_SV  ! Monin-Obukhov  Scale
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  us__SV  ! Friction       Velocity
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  uts_SV  ! Temperature  Turbulent Scale
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  cutsSV  ! Temperature  Turbulent Scale C.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  uqs_SV  ! Spec.Humid.  Turbulent Scale
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  ussbSV  ! Blowing Snow Erosion   Flux  (Buffer) 
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  uss_SV  ! Blowing Snow Turbulent Scale
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  ussxSV  ! Blowing Snow Turbulent Scale (modified)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  usthSV  ! Blowing Snow Erosion Thresh.
! #BD real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  uds_SV  ! Blowing Dust Turbulent Scale               [kg/kg m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  rCDmSV  ! Square  Root Contribut. Drag_m
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  rCDhSV  ! Square  Root Contribut. Drag_h
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0m_SV  ! Momentum     Roughness Length
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0mmSV  !  z0(Momentum,    Time Mean)                        [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0mnSV  !  z0(Momentum,    instanta.)                        [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0roSV  ! Subgrid Topo Roughness Length
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0SaSV  !  z0(Sastrugi  h)                                   [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0e_SV  !  z0(Snow eroded)                                   [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0emSV  !  z0(Snow eroded, Time Mean)                        [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0enSV  !  z0(Snow eroded, instanta.)                        [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0h_SV  ! Heat         Roughness Length
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0hmSV  !  z0(Heat,        Time Mean)                        [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Z0hnSV  !  z0(Heat,        instanta.)                        [m]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  snCaSV  ! Canopy  Snow   Thickness
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  rrCaSV  ! Canopy  Water  Content
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  psivSV  ! Leaf    Water  Potential 
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  TvegSV  ! Vegetation     Temperature

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  TsisSV  ! Snow/Ice/Soil-Water Temperature
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  ro__SV  ! Snow/Ice/Soil-Water VolumicMass
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  eta_SV  ! Snow/Ice/Soil     Water Content
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  G1snSV  ! Snow Dendricity/Sphericity
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  G2snSV  ! Snow Sphericity/Size
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  dzsnSV  ! Snow Layer  Thickness
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  agsnSV  ! Snow Age
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  BufsSV  ! Snow Buffer Layer
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  rusnSV  ! Surficial   Water
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  SWf_SV  ! Normalized  Decay
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  SWS_SV  ! Surficial Water Status
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  HFraSV  ! Frazil      Thickness

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  zWE_SV  ! Current   Snow Thickness                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  zWEcSV  ! Compacted Snow Thickness                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dwemSV  ! Only Melting over dt__SV                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dwerSV  ! Refreezing   over dt__SV                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dwesSV  ! Sublimation  over dt__SV                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  wem0SV  ! Only Melting      Budget                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  wem_SV  ! Only Melting      Budget                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  wer0SV  ! Refreezing        Budget                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  wer_SV  ! Refreezing        Budget                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  wes0SV  ! Sublimation       Budget                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  wes_SV  ! Sublimation       Budget                        [mmWE]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  wee_SV  ! EvapotranspirationBudget                        [mmWE]


! SISVAT OUTPUT       Variables
! -----------------------------

      integer, SAVE         ,ALLOCATABLE ,dimension(:)      ::  no__SV  ! OUTPUT file Unit Number
      integer, SAVE         ,ALLOCATABLE ,dimension(:)      ::  IOi_SV  ! OUTPUT point   i Coordinate (independant txt file)
      integer, SAVE         ,ALLOCATABLE ,dimension(:)      ::  IOj_SV  ! OUTPUT point   j Coordinate (independant txt file)
      integer, SAVE         ,ALLOCATABLE ,dimension(:)      ::  i___SV  ! OUTPUT point   i Coordinate
      integer, SAVE         ,ALLOCATABLE ,dimension(:)      ::  j___SV  ! OUTPUT point   j Coordinate
      integer, SAVE         ,ALLOCATABLE ,dimension(:)      ::  n___SV  ! OUTPUT point   n Coordinate
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  lwriSV  ! OUTPUT point vec Index

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  IRu_SV  ! UPward    IR Flux (effective)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  hSalSV  ! Saltating Layer Height
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  qSalSV  ! Saltating Snow  Concentration
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  RnofSV  ! RunOFF    Intensity

      end module Mod_SISVAT_kkl
