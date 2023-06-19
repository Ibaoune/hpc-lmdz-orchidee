      subroutine SISVAT(SnoMod,BloMod,jjtime)

!--------------------------------------------------------------------------+
!   MAR          SISVAT                               Mon 04-Apr-2011  MAR |
!     SubRoutine SISVAT contains the fortran 77 code of the                |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!--------------------------------------------------------------------------+
!     PARAMETERS:  klonv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   daHost   : Date Host Model                                  |
!     ^^^^^                                                                |
!                                                                          |
!     INPUT:   LSmask   : 1:          Land       MASK                      |
!     ^^^^^               0:          Sea        MASK                      |
!              ivgtSV   = 0,...,12:   Vegetation Type                      |
!              isotSV   = 0,...,12:   Soil       Type                      |
!                         0:          Water,          Liquid (Sea, Lake)   |
!                        12:          Water, Solid           (Ice)         |
!                                                                          |
!     INPUT:   coszSV   : Cosine of the Sun Zenithal Distance          [-] |
!     ^^^^^    sol_SV   : Surface Downward  Solar      Radiation    [W/m2] |
!              IRd_SV   : Surface Downward  Longwave   Radiation    [W/m2] |
!              drr_SV   : Rain  Intensity                        [kg/m2/s] |
!              dsn_SV   : Snow  Intensity                      [mm w.e./s] |
!              dsnbSV   : Snow  Intensity,  Drift Fraction             [-] |
!              dbs_SV   : Drift Amount                           [mm w.e.] |
!              za__SV   : Surface Boundary Layer (SBL) Height          [m] |
!              VV__SV   :(SBL Top)   Wind Velocity                   [m/s] |
!              TaT_SV   : SBL Top    Temperature                       [K] |
!              rhT_SV   : SBL Top    Air  Density                  [kg/m3] |
!              QaT_SV   : SBL Top    Specific  Humidity            [kg/kg] |
!              qsnoSV   : SBL Mean   Snow      Content             [kg/kg] |
!              LAI0SV   : Leaf Area  Index                             [-] |
!              glf0SV   : Green Leaf Fraction                          [-] |
!              alb0SV   : Soil Basic Albedo                            [-] |
!              slopSV   : Surface    Slope                             [-] |
!              dt__SV   : Time  Step                                   [s] |
!                                                                          |
!     INPUT /  isnoSV   = total Nb of Ice/Snow Layers                      |
!     OUTPUT:  ispiSV   = 0,...,nsno: Uppermost Superimposed Ice Layer     |
!     ^^^^^^   iiceSV   = total Nb of Ice      Layers                      |
!              istoSV   = 0,...,5 :   Snow     History (see istdSV data)   |
!                                                                          |
!     INPUT /  alb_SV   : Surface-Canopy Albedo                        [-] |
!     OUTPUT:  emi_SV   : Surface-Canopy Emissivity                    [-] |
!     ^^^^^^   IRs_SV   : Soil           IR Flux  (negative)        [W/m2] |
!              LMO_SV   : Monin-Obukhov               Scale            [m] |
!              us__SV   : Friction          Velocity                 [m/s] |
!              uts_SV   : Temperature       Turbulent Scale          [m/s] |
!              uqs_SV   : Specific Humidity Velocity                 [m/s] |
!              uss_SV   : Blowing Snow      Turbulent Scale          [m/s] |
!              usthSV   : Blowing Snow      Erosion   Threshold      [m/s] |
!              Z0m_SV   : Momentum     Roughness Length                [m] |
!              Z0mmSV   : Momentum     Roughness Length (time mean)    [m] |
!              Z0mnSV   : Momentum     Roughness Length (instantaneous)[m] |
!              Z0SaSV   : Sastrugi     Roughness Length                [m] |
!              Z0e_SV   : Erosion Snow Roughness Length                [m] |
!              Z0emSV   : Erosion Snow Roughness Length (time mean)    [m] |
!              Z0enSV   : Erosion Snow Roughness Length (instantaneous)[m] |
!              Z0roSV   : Subgrid Topo Roughness Length                [m] |
!              Z0h_SV   : Heat         Roughness Length                [m] |
!              snCaSV   : Canopy   Snow     Thickness            [mm w.e.] |
!              rrCaSV   : Canopy   Water    Content                [kg/m2] |
!              psivSV   : Leaf     Water    Potential                  [m] |
!              TvegSV   : Canopy   Temperature                         [K] |
!              TsisSV   : Soil/Ice Temperatures (layers -nsol,-nsol+1,..,0)|
!                       & Snow     Temperatures (layers  1,2,...,nsno) [K] |
!              ro__SV   : Soil/Snow Volumic Mass                   [kg/m3] |
!              eta_SV   : Soil/Snow Water   Content                [m3/m3] |
!              G1snSV   : snow dendricity/sphericity                       |
!              G2snSV   : snow sphericity/grain size                       |
!              dzsnSV   : Snow Layer        Thickness                  [m] |
!              agsnSV   : Snow       Age                             [day] |
!              BufsSV   : Snow Buffer Layer              [kg/m2] .OR. [mm] |
!              BrosSV   : Snow Buffer Layer Density      [kg/m3]           |
!              BG1sSV   : Snow Buffer Layer Dendricity / Sphericity    [-] |
!              BG2sSV   : Snow Buffer Layer Sphericity / Size [-] [0.1 mm] |
!              rusnSV   : Surficial   Water              [kg/m2] .OR. [mm] |
!                                                                          |
!     OUTPUT:  no__SV   : OUTPUT file Unit Number                      [-] |
!     ^^^^^^   i___SV   : OUTPUT point   i Coordinate                  [-] |
!              j___SV   : OUTPUT point   j Coordinate                  [-] |
!              n___SV   : OUTPUT point   n Coordinate                  [-] |
!              lwriSV   : OUTPUT point vec Index                       [-] |
!                                                                          |
!     OUTPUT:  IRu_SV   : Upward     IR Flux (+, upw., effective)      [K] |
!     ^^^^^^   hSalSV   : Saltating Layer Height                       [m] |
!              qSalSV   : Saltating Snow  Concentration            [kg/kg] |
!              RnofSV   : RunOFF Intensity                       [kg/m2/s] |
!                                                                          |
!     Internal Variables:                                                  |
!     ^^^^^^^^^^^^^^^^^^                                                   |
!              NLaysv   = New            Snow Layer Switch             [-] |
!              albisv   : Snow/Ice/Water/Soil Integrated Albedo        [-] |
!              SoCasv   : Absorbed Solar Radiation by Canopy (Normaliz)[-] |
!              SoSosv   : Absorbed Solar Radiation by Surfac.(Normaliz)[-] |
!              tau_sv   : Fraction of Radiation transmitted by Canopy  [-] |
!              TBr_sv   : Brightness Temperature                       [K] |
!              IRupsv   : Upward     IR Flux (-, upw.)              [W/m2] |
!              IRv_sv   : Vegetation IR Flux                        [W/m2] |
!              rrMxsv   : Canopy Maximum Intercepted Rain          [kg/m2] |
!              Sigmsv   : Canopy Ventilation Factor                    [-] |
!              ram_sv   : Aerodynamic Resistance for Momentum        [s/m] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!              HSv_sv   : Vegetation Sensible Heat Flux             [W/m2] |
!              HLv_sv   : Vegetation Latent   Heat Flux             [W/m2] |
!              Rootsv   : Root Water Pump                        [kg/m2/s] |
!              Evp_sv   : Evaporation                              [kg/m2] |
!              EvT_sv   : Evapotranspiration                       [kg/m2] |
!              HSs_sv   : Surface    Sensible Heat Flux + => absorb.[W/m2] |
!              HLs_sv   : Surface    Latent   Heat Flux + => absorb.[W/m2] |
!              Lx_H2O   : Latent Heat of Vaporization/Sublimation   [J/kg] |
!              Tsrfsv   : Surface    Temperature                       [K] |
!              LAI_sv   : Leaf Area  Index (snow included)             [-] |
!              LAIesv   : Leaf Area  Index (effective / transpiration) [-] |
!              glf_sv   : Green Leaf Fraction of NOT fallen Leaves     [-] |
!              sEX_sv   : Verticaly Integrated Extinction Coefficient  [-] |
!              LSdzsv   : Vertical   Discretization Factor             [-] |
!                       =    1. Soil                                       |
!                       = 1000. Ocean                                      |
!              z_snsv   : Snow Pack  Thickness                         [m] |
!              zzsnsv   : Snow Pack  Thickness                         [m] |
!              albssv   : Soil       Albedo                            [-] |
!              Evg_sv   : Soil+Vegetation Emissivity                   [-] |
!              Eso_sv   : Soil+Snow       Emissivity                   [-] |
!              psi_sv   : Soil       Water    Potential                [m] |
!              Khydsv   : Soil   Hydraulic    Conductivity           [m/s] |
!                                                                          |
!              ETVg_d   : VegetationEnergy Power         Forcing    [W/m2] |
!              ETSo_0   : Snow/Soil Energy Power, before Forcing    [W/m2] |
!              ETSo_1   : Snow/Soil Energy Power, after  Forcing    [W/m2] |
!              ETSo_d   : Snow/Soil Energy Power         Forcing    [W/m2] |
!              EqSn_0   : Snow      Energy, before Phase Change     [J/m2] |
!              EqSn_1   : Snow      Energy, after  Phase Change     [J/m2] |
!              EqSn_d   : Snow      Energy,       net    Forcing    [J/m2] |
!              Enrsvd   : SVAT      Energy Power         Forcing    [W/m2] |
!              Enrbal   : SVAT      Energy Balance                  [W/m2] |
!              Wats_0   : Soil Water,  before Forcing                 [mm] |
!              Wats_1   : Soil Water,  after  Forcing                 [mm] |
!              Wats_d   : Soil Water          Forcing                 [mm] |
!              SIWm_0   : Snow        initial Mass               [mm w.e.] |
!              SIWm_1   : Snow        final   Mass               [mm w.e.] |
!              SIWa_i   : Snow Atmos. initial Forcing            [mm w.e.] |
!              SIWa_f   : Snow Atmos. final   Forcing(noConsumed)[mm w.e.] |
!              SIWe_i   : SnowErosion initial Forcing            [mm w.e.] |
!              SIWe_f   : SnowErosion final   Forcing(noConsumed)[mm w.e.] |
!              SIsubl   : Snow sublimed/deposed  Mass            [mm w.e.] |
!              SImelt   : Snow Melted            Mass            [mm w.e.] |
!              SIrnof   : Surficial Water + Run OFF Change       [mm w.e.] |
!              SIvAcr   : Sea-Ice    vertical Acretion           [mm w.e.] |
!              Watsvd   : SVAT Water          Forcing                 [mm] |
!              Watbal   : SVAT Water  Balance                       [W/m2] |
!                                                                          |
!              dsn_Ca,snCa_n :     Snow Contribution to the Canopy[m w.e.] |
!              drr_Ca,rrCa_n,drip: Rain Contribution to the Canopy [kg/m2] |
!              vk2      : Square of Von Karman Constant                [-] |
!              sqrCm0   : Factor of   Neutral Drag Coeffic.Momentum  [s/m] |
!              sqrCh0   : Factor of   Neutral Drag Coeffic.Heat      [s/m] |
!              EmiVeg   : Vegetation    Emissivity                     [-] |
!              EmiSol   : Soil          Emissivity                     [-] |
!              EmiSno   : Snow          Emissivity                     [-] |
!              EmiWat   : Water         Emissivity                     [-] |
!              Z0mSea   :          Sea  Roughness Length               [m] |
!              Z0mLnd   :          Land Roughness Length               [m] |
!              sqrrZ0   : u*t/u*                                           |
!              f_eff    : Marticorena & B. 1995 JGR (20)                   |
!              A_Fact   : Fundamental * Roughness                          |
!              Z0mBSn   :         BSnow Roughness Length               [m] |
!              Z0mBS0   : Mimimum BSnow Roughness Length (blown* )     [m] |
!              Z0m_Sn   :          Snow Roughness Length (surface)     [m] |
!              Z0m_S0   : Mimimum  Snow Roughness Length               [m] |
!              Z0m_S1   : Maximum  Snow Roughness Length               [m] |
!              Z0_GIM   : Minimum GIMEX Roughness Length               [m] |
!              Z0_ICE   : Sea Ice ISW   Roughness Length               [m] |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Possibility                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                          |
!     #AE: TURBULENCE: Aerosols Erosion / Turbulent Diffusion Coeff.       |
!     #BD: TraCer   Aeolian Erosion  Submodel           is turned ON       |
!     #BS: Explicit Cloud MICROPHYSICS: Blow. *(Snow)         Model)       |
!     #SN: SNOW Model                               may be turned ON       |
!     #NP: SNOW Model: Snow Properties may be those of Polar Snow          |
!     #GL: SNOW Model: ETH-Camp & Greenland 3D simulations                 |
!     #MB: SNOW Model: Erosion Efficiency (Marticorena & Berga.1995)       |
!     #SI: SISVAT: Sea-Ice Fraction calculated from prescribed SST         |
!     #MT: SISVAT: Monin-Obukhov Theory is linearized (Garrat schem)       |
!     #SH: Soil /Vegetation Model: Hapex-Sahel   Vegetation     DATA       |
!     #OR: SBL: Orography Roughness included from SL_z0 in MARdom          |
!     #ZS: SBL: Mom.: Roughn.Length= F(u*) Wang  MWR 129     ,  Sea        |
!     #TZ: SBL: Mom.: Roughn.Length= Typical value in polar models         |
!     #SZ: SBL: Mom.: Roughn.Length= F(u*) Andreas &al.(1995)  Snow        |
!     #ZA: SBL: Mom.: Roughn.Length= F(u*) Andreas &al.(2004), Snow        |
!     #za: SBL: Mom.: Roughn.Length= F(u*) Andreas &al.(2004), Snow (native|
!     #RN: SBL: Heat: Roughn.Length= F(u*,z0)  Andreas (1987)  Snow        |
!     #ZM: SBL: M/H   Roughn.Length: Box Moving Average (in Time)          |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Col de Porte                         |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^                                      |
!     #CP: Col de Porte  Turbulence        Parameterization                |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #zs: SBL: Mom.: Roughn.Length= F(u*) Wang  MWR 129 bis ,  Sea        |
!     #ZN: SBL: Mom.: Roughn.Length= F(u*) Shao  & Lin (1999), Snow        |
!     #CM: SBL: Z0mL  Roughn.Length= F(glf)                                |
!     #FL: SISVAT: LAI Assignation and Fallen Leaves Correction            |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_iii_jjj_n     | #E0: OUTPUT on ASCII  File (SISVAT Variables) |
!   #                      |      Energy Budg. Verif.: Soil+(Sea-Ice)+Snow |
!   #                      |(#E0  MUST BE PREPROCESSED BEFORE #e1 & #e2 !) |
!   # SISVAT_iii_jjj_n     | #m0: OUTPUT/Verification: H2O    Conservation |
!                          |                                               |
!   # stdout               | #s0: OUTPUT of Snow Buffer Layer              |
!                          |      unit  6, SubRoutine  SISVAT     **ONLY** |
!   # stdout               | #s2: OUTPUT of SnowFall, Snow Buffer          |
!                          |      unit  6, SubRoutine  SISVAT_BSn, _qSn    |
!   # stdout               | #b0: OUTPUT of Snow Erosion                   |
!                          |      unit  6, SubRoutine  SISVAT_BSn **ONLY** |
!   # stdout               | #sf: OUTPUT of SnowFall, Z0 and Drag Coeff.   |
!                          |      unit  6, SubRoutines PHY_SISVAT, SISVAT  |
!   # stdout               | #sz: OUTPUT of Roughness Length & Drag Coeff. |
!                          |      unit  6, SubRoutine  SISVAT     **ONLY** |
!                                                                          |
!     SUGGESTIONS of MODIFICATIONS: see lines beginning with "C +!!!"      |
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                         |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV
      USE VARlSV
      USE VARdSV
      USE VAR0SV
      USE VARxSV
      USE VARdCP

      USE VARySV

      IMPLICIT NONE

      logical   SnoMod
      logical   BloMod
      integer   jjtime


! Internal Variables
! ==================

! Non Local
! ---------

      real      TBr_sv(klonv)                 ! Brightness Temperature
      real      IRdwsv(klonv)                 ! DOWNward   IR Flux
      real      IRupsv(klonv)                 ! UPward     IR Flux
      real      d_Bufs,Bufs_N                 ! Buffer Snow Layer Increment
      real      Buf_ro,Bros_N                 ! Buffer Snow Layer Density
c #NP real      BufPro,VV__10                 ! Buffer Snow Layer Density
      real      Buf_G1,BG1__N                 ! Buffer Snow Layer Dendr/Sphe[-]
      real      Buf_G2,BG2__N                 ! Buffer Snow Layer Spher/Size[-]
      real      Bdzssv(klonv)                 ! Buffer Snow Layer Thickness
      real      z_snsv(klonv)                 ! Snow-Ice, current Thickness

! Energy         Budget
! ~~~~~~~~~~~~~~~~~~~~~
! #e1 real      ETVg_d(klonv)                 ! VegetationPower, Forcing
! #e1 real      ETSo_0(klonv)                 ! Soil/Snow Power, before Forcing
! #e1 real      ETSo_1(klonv)                 ! Soil/Snow Power, after  Forcing
! #e1 real      ETSo_d(klonv)                 ! Soil/Snow Power, Forcing
! #e1 real      EqSn_0(klonv)                 ! Snow Energy, befor Phase Change
! #e1 real      EqSn_1(klonv)                 ! Snow Energy, after Phase Change
! #e1 real      EqSn_d(klonv)                 ! Energy in Excess

! OUTPUT/Verification: H2O    Conservation
! #m0 real      Wats_0(klonv)                 ! Soil Water,  before Forcing
! #m0 real      Wats_1(klonv)                 ! Soil Water,  after  Forcing
! #m0 real      Wats_d(klonv)                 ! Soil Water,         Forcing

! OUTPUT/Verification: * Mass Conservation
! #m1 real      SIsubl(klonv)                 ! Snow Sublimed/Deposed Mass
! #m1 real      SImelt(klonv)                 ! Snow Melted           Mass
! #m1 real      SIrnof(klonv)                 ! Local Surficial Water + Run OFF

! OUTPUT/Verification: SeaIce Conservation
! #m2 real      SIvAcr(klonv)                 ! Sea-Ice      Vertical Acretion


! Local
! -----

c #MT logical      Garrat                     ! SBL     Scheme    Switch
      character* 1 SepLab                     ! OUTPUT ASCII File Labels
      character* 6 FilLab                     !
      character*16 FilNam                     !
      common/SISVAT_loc_abc/SepLab,FilLab     !

      integer               noUNIT            ! OUTPUT File  Unit Number
      integer               nwUNIT            ! OUTPUT File  Unit Number (New)
      common/SISVAT_loc_num/nwUNIT            !

      integer   iwr
      integer   ikl   ,isn   ,isl   ,ist      !
      integer   ist__s,ist__w                 ! Soil/Water Body Identifier
      integer   growth                        ! Seasonal               Mask
      integer   LISmsk                        ! Land+Ice / Open    Sea Mask
      integer   LSnMsk                        ! Snow-Ice / No Snow-Ice Mask
      integer   IceMsk,IcIndx(klonv)          !      Ice / No      Ice Mask
      integer   SnoMsk                        ! Snow     / No Snow     Mask

      real      drr_Ca,rrCa_n,drip            ! Rain Contribution to the Canopy
      real      dsn_Ca,snCa_n,FallOK(klonv)   ! Snow Contribution to the Canopy
      real      roSMin,roSn_1,roSn_2,roSn_3   ! Fallen Snow Density (PAHAUT)
      real      Dendr1,Dendr2,Dendr3          ! Fallen Snow Dendric.(GIRAUD)
      real      Spher1,Spher2,Spher3,Spher4   ! Fallen Snow Spheric.(GIRAUD)
      real      Polair                        ! Polar  Snow Switch
      real      PorSno,Por_BS,Salt_f,PorRef   !
c #sw real      PorVol,rWater                 !
c #sw real      rusNEW,rdzNEW,etaNEW          !
      real      ro_new                        !
      real      TaPole                        ! Maximum     Polar Temperature
      real      T__Min                        ! Minimum realistic Temperature

      real      EmiVeg                        ! Emissivity of Vegetation 
      real      EmiSol                        ! Emissivity of       Soil
      real      EmiSno                        ! Emissivity of            Snow
      real      EmiWat                        ! Emissivity of a Water Area
      real      vk2                           ! Square of Von Karman Constant
      real      u2star                        !(u*)**2
      real      fallen                        ! Fallen   Leaves         Switch
      real      Z0mSea,Z0hSea                 !          Sea  Roughness Length
      real      Z0mLnd                        !          Land Roughness Length
c #ZN real      sqrrZ0                        ! u*t/u*
      real      f_eff                         ! Marticorena & B. 1995 JGR (20)
      real      A_Fact                        ! Fundamental * Roughness
      real      Z0m_nu                        ! Smooth R Snow Roughness Length
      real      Z0mBSn                        !         BSnow Roughness Length
      real      Z0mBS0                        ! Mimimum BSnow Roughness Length
      real      Z0m_S0                        ! Mimimum  Snow Roughness Length
      real      Z0m_S1                        ! Maximum  Snow Roughness Length
c #SZ real      Z0Sa_N                        ! Regime   Snow Roughness Length
c #SZ real      Z0SaSi                        ! 1.IF Rgm Snow Roughness Length 
c #GL real      Z0_GIM                        ! Mimimum GIMEX Roughness Length
      real      Z0_ICE                        ! Sea-Ice ISW   Roughness Length
      real      Z0m_Sn,Z0m_90                 ! Snow  Surface Roughness Length
      real      SnoWat                        ! Snow Layer    Switch
c #RN real      rstar,alors                   !
c #RN real      rstar0,rstar1,rstar2          !
      real      SameOK                        ! 1. => Same Type of Grains
      real      G1same                        ! Averaged G1,  same Grains
      real      G2same                        ! Averaged G2,  same Grains
      real      typ__1                        ! 1. => Lay1 Type: Dendritic
      real      zroNEW                        ! dz X ro, if fresh Snow
      real      G1_NEW                        ! G1,      if fresh Snow
      real      G2_NEW                        ! G2,      if fresh Snow
      real      zroOLD                        ! dz X ro, if old   Snow
      real      G1_OLD                        ! G1,      if old   Snow
      real      G2_OLD                        ! G2,      if old   Snow
      real      SizNEW                        ! Size,    if fresh Snow
      real      SphNEW                        ! Spheric.,if fresh Snow
      real      SizOLD                        ! Size,    if old   Snow
      real      SphOLD                        ! Spheric.,if old   Snow
      real      Siz_av                        ! Averaged    Grain Size
      real      Sph_av                        ! Averaged    Grain Spher.
      real      Den_av                        ! Averaged    Grain Dendr.
      real      DendOK                        ! 1. => Average is  Dendr.
      real      G1diff                        ! Averaged G1, diff. Grains
      real      G2diff                        ! Averaged G2, diff. Grains
      real      G1                            ! Averaged G1
      real      G2                            ! Averaged G2

! Energy       Budget
! ~~~~~~~~~~~~~~~~~~~
! #e1 integer             noEBal              ! Energy Imbalances Counter
! #e1 common/SISVAT__EBal/noEBal              !
! #e1 real        Enrsvd(klonv)               ! Soil+Vegetat  Power  Forcing
! #e1 real        EnsBal                      ! Soil+Snow   , Power  Balance
! #e1 real        EnvBal                      !      Vegetat, Power  Balance

! OUTPUT/Verification: H2O    Conservation
! #m0 integer             noWBal              ! Water  Imbalances Counter
! #m0 common/SISVAT__WBal/noWBal              !
! #m0 real        Watsv0(klonv)               ! Soil+Vegetat, before Forcing
! #m0 real        Watsvd(klonv)               ! Soil+Vegetat  Water  Forcing
! #m0 real        Watbal                      ! Soil+Vegetat, Water  Balance

! OUTPUT/Verification: * Mass Conservation
! #m1 integer             noSBal              ! Water  Imbalances Counter
! #m1 common/SISVAT__SBal/noSBal              !
! #m1 real        SIWm_0(klonv),SIWm_1(klonv) ! Snow Initial/Final        Mass
! #m1 real        SIWa_i(klonv),SIWa_f(klonv) ! Snow Initial/Final ATM Forcing
! #m1 real        SIWe_i(klonv),SIWe_f(klonv) ! Snow Initial/Final BLS Forcing
! #m1 real        SnoBal                      ! Snow Pack     Mass   Balance


! Internal DATA
! =============

c #MT data      Garrat /.true. /              ! SBL     Scheme    Switch
      data      T__Min / 200.00/              ! Minimum realistic Temperature
      data      TaPole / 263.15/              ! Maximum Polar     Temperature
      data      roSMin /  30.  /              ! Minimum Snow  Density
      data      roSn_1 / 109.  /              ! Fall.Sno.Density, Indep. Param.
      data      roSn_2 /   6.  /              ! Fall.Sno.Density, Temper.Param.
      data      roSn_3 /  26.  /              ! Fall.Sno.Density, Wind   Param.
      data      Dendr1 /  17.12/              ! Fall.Sno.Dendric.,Wind 1/Param.
      data      Dendr2 / 128.  /              ! Fall.Sno.Dendric.,Wind 2/Param.
      data      Dendr3 / -20.  /              ! Fall.Sno.Dendric.,Indep. Param.
      data      Spher1 /   7.87/              ! Fall.Sno.Spheric.,Wind 1/Param.
      data      Spher2 /  38.  /              ! Fall.Sno.Spheric.,Wind 2/Param.
      data      Spher3 /  50.  /              ! Fall.Sno.Spheric.,Wind 3/Param.
      data      Spher4 /  90.  /              ! Fall.Sno.Spheric.,Indep. Param.
      data      EmiSol /   0.99999999/        ! 0.94Emissivity of Soil
      data      EmiVeg /   0.99999999/        ! 0.98Emissivity of Vegetation
      data      EmiWat /   0.99999999/        ! Emissivity of a Water Area
      data      EmiSno /   0.99999999/        ! Emissivity of Snow
!     DATA      Emissivities                  ! Pielke, 1984, pp. 383,409

      data      fallen /   0.    /            ! Fallen  Leaves         Switch
      data      Z0mBS0 /   0.5e-6/            ! MINimum Snow Roughness Length 
                                              ! for Momentum if Blowing Snow
                                              ! Gall?e et al. 2001 BLM 99 (19)
      data      Z0m_S0/    0.00005/           ! MINimum Snow Roughness Length
c #MG data      Z0m_S0/    0.00200/           ! MINimum Snow Roughness Length
                                              ! MegaDunes    included
      data      Z0m_S1/    0.030  /           ! MAXimum Snow Roughness Length
                                              !        (Sastrugis)
c #GL data      Z0_GIM/    0.0013/            ! Ice Min Z0 = 0.0013 m (Broeke)
!                                             ! Old Ice Z0 = 0.0500 m (Bruce)
!                                             !              0.0500 m (Smeets)
!                                             !              0.1200 m (Broeke)
      data      Z0_ICE/    0.0010/            ! Sea-Ice Z0 = 0.0010 m (Andreas)
!                                             !    (Ice Station Weddel -- ISW)
      vk2    =  vonKrm  *  vonKrm             ! Square of Von Karman Constant
c #FL fallen =             1.                 ! Fallen  Leaves         Switch


! BEGIN.main.
! SISVAT Forcing VERIFICATION
! ===========================

      IF (.not.iniOUT)                                              THEN
               iniOUT = .true.
        IF (IRs_SV(1).gt.-eps6)
     .  write(6,600) 
 600    format(/,'### SISVAT ERROR, Soil IR Upward  not defined ###',
     .         /,'###               Initialize and Store IRs_SV ###')
        IF (IRs_SV(1).gt.-eps6) THEN
          write(*,*)'ikl',ikl,'IR',IRs_SV(ikl)

          IRs_SV(ikl)=-IRs_SV(ikl)
          isn= isnoSV(ikl)
          write(*,*) isnoSV(ikl),TsisSV(ikl,isn),dzsnSV(ikl,isn) 
        ENDIF

! OUTPUT
! ======

                    FilLab              ='SISVAT'
                    SepLab              ='_'
                    nwUNIT              = 51
      END IF

c #E0 DO ikl=1,knonv
c #E0   IF   (lwriSV(ikl).ne.0.AND.no__SV(lwriSV(ikl)).eq.0)        THEN 
c #E0                nwUNIT              = nwUNIT+1
c #E0                no__SV(lwriSV(ikl)) = nwUNIT
c #E0      write(FilNam,'(a6,a1,2(i3.3,a1),i1)')
c #E0.           FilLab,SepLab,i___SV(lwriSV(ikl)),
c #E0.                  SepLab,j___SV(lwriSV(ikl)),
c #E0.                  SepLab,n___SV(lwriSV(ikl))
c #E0      open(unit=nwUNIT,status='unknown',file=FilNam)
c #E0      rewind    nwUNIT
c #E0   END IF
c #E0 END DO

c #E0 DO ikl=1,knonv
c #E0   IF (lwriSV(ikl).ne.0)                                       THEN 
c #E0           noUNIT=no__SV(lwriSV(ikl)) 
c #E0     write(noUNIT,5000) daHost,i___SV(lwriSV(ikl)),
c #E0.                              j___SV(lwriSV(ikl)),
c #E0.                              n___SV(lwriSV(ikl)),
c #E0.                                     Z0m_SV(ikl) ,
c #E0.                                     albisv(ikl)
 5000     format(
     .       /,              a18,'|           Grid Point ',2i4,
     .                                           ' (',i2,')',
     .         '    | Z0m =',f12.6,' | Albedo = ',f6.3,' |',
     .       /,' -------+',7('---------+'),2('--------+'))
c #E0   END IF
c #E0 END DO


! "Soil" Humidity of Water Bodies
! ===============================

      DO ikl=1,knonv
          ist    =      isotSV(ikl)                       ! Soil Type
          ist__s =  min(ist, 1)                           ! 1 => Soil
          ist__w =  1 - ist__s                            ! 1 => Water Body
        DO isl=-nsol,0
          eta_SV(ikl,isl) = eta_SV(ikl,isl) * ist__s      ! Soil 
     .                    + etadSV(ist)     * ist__w      ! Water Body
        END DO


! Vertical Discretization Factor
! ==============================

          LSdzsv(ikl)     =                   ist__s      ! Soil
     .                    + OcndSV          * ist__w      ! Water Body
      END DO


! Vegetation Temperature Limits
! =============================

          DO ikl=1,knonv
            TvegSV(ikl) = max(TvegSV(ikl),T__Min)         ! T__Min = 200.K


! LAI Assignation and Fallen Leaves Correction (#FL)
! ==================================================

            LAI0SV(ikl) =     LAI0SV(ikl)*min(1,ivgtSV(ikl)) ! NO LAI if 
!                                                            ! no vegetation
            glf_sv(ikl) =     glf0SV(ikl)
c #FL       glf_sv(ikl) =     1.
            LAI_sv(ikl) =     LAI0SV(ikl)
c #FL.               *        glf0SV(ikl) 
          END DO


! LAI in Presence of Snow
! =======================

!         ASSUMPTION: LAI decreases   when Snow Thickness increases, 
!         ^^^^^^^^^^      becoming  0 when Snow Thickn. = Displac.Height
          DO ikl=1,knonv
            LAI_sv(ikl) =     LAI_sv(ikl)
     .               * (1.0 - zzsnsv(       ikl, isnoSV(ikl))
     .                      /(DH_dSV(ivgtSV(ikl))+eps6)      )
            LAI_sv(ikl) = max(LAI_sv(ikl),zer0)
            LAI_sv(ikl) = min(LAI_sv(ikl),ea_Max)
          END DO


! Interception of Rain by the Canopy
! ==================================

! OUTPUT/Verification: H2O    Conservation: Vegetation Forcing
! #m0     DO ikl=1,knonv
! #m0       Watsv0(ikl) =      rrCaSV(ikl)           ! Canopy Water Cont.
! #m0       Watsvd(ikl) =      drr_SV(ikl)           ! Precipitation
! #m0     END DO


! New Canopy Water Content
! ------------------------

          DO ikl=1,knonv
            rrMxsv(ikl) = 0.2*max( eps6,LAI_sv(ikl)) ! Precip. Max. Intercept.
            Sigmsv(ikl) = 1.0-exp(-half*LAI_sv(ikl)) ! Canopy Ventilation Coe.
!                                                    ! (DR97, eqn 3.6)
            drr_Ca      = drr_SV(ikl)  *Sigmsv(ikl)  ! Intercepted Rain
     .                                 *dt__SV       !
            rrCa_n      = rrCaSV(ikl)  +drr_Ca       ! New Canopy Water Contnt
                                                     ! (DR97, eqn 3.28)
            drip        = rrCa_n       -rrMxsv(ikl)  ! Water  Drip
            drip        =      max(zer0,drip)        !
            rrCa_n      = rrCa_n       -drip         !
            drr_SV(ikl) = drr_SV(ikl) +(rrCaSV(ikl)  ! Update Rain  Contribut.
     .                                 -rrCa_n     ) !
     .                                 /dt__SV       !
            rrCaSV(ikl) = rrCa_n                     ! Upd.Canopy Water Contnt


! Interception of Snow by the Canopy
! ==================================

            dsn_Ca      = dsn_SV(ikl)  *Sigmsv(ikl)  ! Intercepted Snow
     .                                 *dt__SV       !
            snCa_n      = snCaSV(ikl)  +dsn_Ca       ! New Canopy Snow Thickn.
            drip        = snCa_n       -rrMxsv(ikl)  !
            drip        =      max(zer0,drip)        !
            snCa_n      = snCa_n       -drip         !
            dsn_SV(ikl) = dsn_SV(ikl) +(snCaSV(ikl)  ! Update Snow  Contribut.
     .                                 -snCa_n     ) ! 
     .                                 /dt__SV       !
            snCaSV(ikl) = snCa_n                     ! Upd.Canopy Snow Thickn.
          END DO 


! Snow Fall from the Canopy
! =========================

!         ASSUMPTION: snow fall from the canopy, 
!         ^^^^^^^^^^  when the temperature of the vegetation is positive
!               (.OR. when snow over the canopy is saturated  with water)

          DO ikl=1,knonv
            FallOK(ikl) =  max(zer0,sign(un_1,TvegSV(ikl)-Tf_Sno+eps6))
     .                  *  max(zer0,sign(un_1,snCaSV(ikl)       -eps6))
            dsn_SV(ikl) =      dsn_SV(ikl)   +snCaSV(ikl)*FallOK(ikl) 
     .                                       /dt__SV
            snCaSV(ikl) =      snCaSV(ikl) * (1.         -FallOK(ikl))


! Blowing Particles Threshold Friction velocity
! =============================================

c #AE       usthSV(ikl) =                     1.0e+2
          END DO


! Contribution of Snow to the Surface Snow Pack
! =============================================

      IF (SnoMod)                                                 THEN


! OUTPUT/Verification: * Mass Conservation
! #m1   DO ikl=1,knonv
! #m1     SIWa_i(ikl) =(drr_SV(ikl) + dsn_SV(ikl))   *dt__SV         ![mm w.e.]
! #m1     SIWe_i(ikl) = dbs_SV(ikl)                                  !
! #m1     SIWm_0(ikl) = BufsSV(ikl) + HFraSV(ikl)    *rhoIce         !
! #m1   DO isn=1,nsno                                                !
! #m1     SIWm_0(ikl) = SIWm_0(ikl) + dzsnSV(ikl,isn)*ro__SV(ikl,isn)!
! #m1   END DO                                                       !
! #m1   END DO                                                       !


! Blowing Snow
! ------------

!                         **********
c #SN   IF (BloMod)  call SISVAT_BSn(BloMod)
!                         **********

!                         **********
! #ve                call SISVAT_wEq('_BSn  ',1)
!                         **********


! Sea Ice
! -------

!            **********
c #SI   call SISVAT_SIc
! #m2.                 (SIvAcr)
!            **********

!            **********
! #ve   call SISVAT_wEq('_SIc  ',0)
!            **********


! Buffer Layer
! ------------

          DO ikl=1,knonv
            BufsSV(ikl) =      BufsSV(ikl)              !     [mm w.e.]
            d_Bufs      =  max(dsn_SV(ikl) *dt__SV,0.)  ! i.e., [kg/m2]
            dsn_SV(ikl) =      0.                       !
            Bufs_N      =      BufsSV(ikl) +d_Bufs      !

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #s0.          nn__SV(ikl).EQ.nwr_SV)
! #s0.      write(6,6601)      BufsSV(ikl) ,d_Bufs,Bufs_N
 6601       format(/,'Buffer *: ',3e15.6)
 
! Snow Density
! ^^^^^^^^^^^^
            Polair      =      0.00
c #NP       Polair      =  max(zer0,                    !
c #NP.                         sign(un_1,TaPole         !
c #NP.                                  -TaT_SV(ikl)))  !
            Buf_ro      =  max( rosMin,                 ! Fallen Snow Density
     .      roSn_1+roSn_2*     (TaT_SV(ikl)-Tf_Sno)     ! [kg/m3]
     .            +roSn_3*sqrt( VV__SV(ikl)))           ! Pahaut    (CEN)
c #NP       VV__10      =       VV__SV(ikl)             !
c #NP.         *log(10.        /Z0m_SV(ikl))            !
c #NP.         /log(za__SV(ikl)/Z0m_SV(ikl))            !
c #NP       BufPro      =  max( rosMin,                 ! Fallen Snow Density
c #NP.         104. *sqrt( max( VV__10  -   6.0,0.0)))  ! Kotlyakov (1961)
            Bros_N      = (1. - Polair) *   Buf_ro      ! Temperate Snow
c #NP.                        + Polair  *   BufPro      ! Polar     Snow
 
! Instantaneous Density of deposited blown Snow (de Montmollin, 1978)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #BS       PorSno =      1.0d00     -  BSnoRo
c #BS.                               /  rhoIce
c #BS       Salt_f =      usthSV(ikl)/  max(eps6,   us__SV(ikl))
c #BS       Salt_f =  min(Salt_f     ,  un_1)
c #BS       PorRef =      PorSno     /  max(eps6,1.-PorSno)
c #BS.               +log(Salt_f)
c #BS       Por_BS =      PorRef     /          (1.+PorRef)
c #BS       ro_new =      rhoIce     *          (1.-Por_BS)
c #BS       ro_new =  max(ro_new     ,  BSnoRo)
c #BS       Bros_N      = Bros_N     * (1.0-dsnbSV(ikl))
c #BS.                  + ro_new     *      dsnbSV(ikl)
 
! Instantaneous Density IF deposited blown Snow (Melted* from Canopy)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Bros_N      = Bros_N     * (1.0-FallOK(ikl))!
     .                  + 300.       *      FallOK(ikl) !

! Time averaged Density of deposited blown Snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            BrosSV(ikl) =(Bros_N     *      d_Bufs      !
     .                   +BrosSV(ikl)*      BufsSV(ikl))!
     .                   /         max(eps6,Bufs_N)     !


! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #s0.          nn__SV(ikl).EQ.nwr_SV)
! #s0.      write(6,6602) Buf_ro,Bros_N,BrosSV(ikl),dsnbSV(ikl)
 6602       format('rho    *: ',3e15.6,'    dsnbSV: ',e15.6)
 
!  S.Falling Snow Properties (computed as in SISVAT_zAg)
!    ^^^^^^^^^^^^^^^^^^^^^^^
            Buf_G1      =  max(-G1_dSV,                 ! Temperate Snow
     .               min(Dendr1*VV__SV(ikl)-Dendr2,     !     Dendricity
     .                   Dendr3                   ))    !
            Buf_G2      =  min( Spher4,                 ! Temperate Snow
     .               max(Spher1*VV__SV(ikl)+Spher2,     !     Sphericity
     .                   Spher3                   ))    !
            Buf_G1      = (1. - Polair) *   Buf_G1      ! Temperate Snow
     .                        + Polair  *   G1_dSV      ! Polar     Snow
            Buf_G2      = (1. - Polair) *   Buf_G2      ! Temperate Snow
     .                        + Polair  *   ADSdSV      ! Polar     Snow
                G1      =                   Buf_G1      ! NO  Blown Snow
                G2      =                   Buf_G2      ! NO  Blown Snow

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #s0.          nn__SV(ikl).EQ.nwr_SV)
! #s0.      write(6,6603)       BG1sSV(ikl),BG2sSV(ikl)  
 6603       format('G1,G2  *: ',3e15.6)
 
! S.1. Meme  Type  de Neige  / same Grain Type
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c #BS       SameOK  =  max(zer0,
c #BS.                     sign(un_1,    Buf_G1             *G1_dSV
c #BS.                                 - eps_21                    ))
c #BS       G1same  = ((1.0-dsnbSV(ikl))*Buf_G1+dsnbSV(ikl) *G1_dSV)
c #BS       G2same  = ((1.0-dsnbSV(ikl))*Buf_G2+dsnbSV(ikl) *ADSdSV)
!           Blowing Snow Properties:                         G1_dSV, ADSdSV
 
! S.2. Types differents / differents Types
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c #BS       typ__1  =  max(zer0,sign(un_1,eps6-Buf_G1))   ! =1.=> Dendritic
c #BS       zroNEW  =     typ__1  *(1.0-dsnbSV(ikl))      ! fract.Dendr.Lay.
c #BS.              + (1.-typ__1) *     dsnbSV(ikl)       !
c #BS       G1_NEW  =     typ__1  *Buf_G1                 ! G1 of Dendr.Lay.
c #BS.              + (1.-typ__1) *G1_dSV                 !
c #BS       G2_NEW  =     typ__1  *Buf_G2                 ! G2 of Dendr.Lay.
c #BS.              + (1.-typ__1) *ADSdSV                 !
c #BS       zroOLD  = (1.-typ__1) *(1.0-dsnbSV(ikl))      ! fract.Spher.Lay.
c #BS.              +     typ__1  *     dsnbSV(ikl)       !
c #BS       G1_OLD  = (1.-typ__1) *Buf_G1                 ! G1 of Spher.Lay.
c #BS.              +     typ__1  *G1_dSV                 !
c #BS       G2_OLD  = (1.-typ__1) *Buf_G2                 ! G2 of Spher.Lay.
c #BS.              +     typ__1  *ADSdSV                 !
c #BS       SizNEW  =    -G1_NEW  *DDcdSV/G1_dSV          ! Size  Dendr.Lay.
c #BS.               +(1.+G1_NEW         /G1_dSV)         !
c #BS.                  *(G2_NEW  *DScdSV/G1_dSV          !
c #BS.               +(1.-G2_NEW         /G1_dSV)*DFcdSV) !
c #BS       SphNEW  =     G2_NEW         /G1_dSV          ! Spher.Dendr.Lay.
c #BS       SizOLD  =     G2_OLD                          ! Size  Spher.Lay.
c #BS       SphOLD  =     G1_OLD         /G1_dSV          ! Spher.Spher.Lay.
c #BS       Siz_av =     (zroNEW*SizNEW+zroOLD*SizOLD)    ! Averaged Size
c #BS       Sph_av = min( zroNEW*SphNEW+zroOLD*SphOLD     !
c #BS.                   ,  un_1)                         ! Averaged Sphericity
c #BS       Den_av = min((Siz_av -(    Sph_av *DScdSV     !
c #BS.                            +(1.-Sph_av)*DFcdSV))   !
c #BS.                 / (DDcdSV -(    Sph_av *DScdSV     !
c #BS.                            +(1.-Sph_av)*DFcdSV))   !
c #BS.                   ,  un_1)                         !
c #BS       DendOK  = max(zer0,                           !
c #BS.                    sign(un_1,     Sph_av *DScdSV   ! Small   Grains
c #BS.                              +(1.-Sph_av)*DFcdSV   ! Faceted Grains
c #BS.                              -    Siz_av        )) !
!           REMARQUE: le  type moyen (dendritique ou non) depend
!           ^^^^^^^^  de la  comparaison avec le diametre optique
!                     d'une neige recente de   dendricite nulle
!           REMARK:   the mean type  (dendritic   or not) depends
!           ^^^^^^    on the comparaison with the optical diameter
!                     of a recent snow    having zero dendricity
 
c #BS       G1diff  =(   -DendOK *Den_av
c #BS.               +(1.-DendOK)*Sph_av) *G1_dSV
c #BS       G2diff  =     DendOK *Sph_av  *G1_dSV
c #BS.               +(1.-DendOK)*Siz_av
c #BS       G1      =     SameOK *G1same
c #BS.               +(1.-SameOK)*G1diff
c #BS       G2      =     SameOK *G2same
c #BS.               +(1.-SameOK)*G2diff
 
            BG1__N      =((1. - FallOK(ikl))*   G1      !
     .                        + FallOK(ikl) *   99.)    ! Melted *  from Canopy
     .                  *       d_Bufs/max(eps6,d_Bufs) !
            BG2__N      =((1. - FallOK(ikl))*   G2      !
     .                        + FallOK(ikl) *   30.)    ! Melted *  from Canopy
     .                  *       d_Bufs/max(eps6,d_Bufs) !

!  S.Buffer  Snow Properties (computed as in SISVAT_zAg)
!    ^^^^^^^^^^^^^^^^^^^^^^^
            Buf_G1      =       BG1__N                  ! Falling   Snow
            Buf_G2      =       BG2__N                  ! Falling   Snow

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #s0.          nn__SV(ikl).EQ.nwr_SV)
! #s0.      write(6,6604)      Buf_G1      ,Buf_G2         ,FallOK(ikl)
! #s0.                                                     ,TvegSV(ikl)
 6604       format('G1,G2 F*: ',3e15.6,'    T__Veg: ',e15.6)
 
! S.1. Meme  Type  de Neige  / same Grain Type
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            SameOK  =  max(zer0,
     .                     sign(un_1,    Buf_G1 *BG1sSV(ikl)
     .                                 - eps_21                    ))
            G1same  = (d_Bufs*Buf_G1+BufsSV(ikl)*BG1sSV(ikl))
     .                     / max(eps6,Bufs_N)
            G2same  = (d_Bufs*Buf_G2+BufsSV(ikl)*BG2sSV(ikl))
     .                     /max(eps6,Bufs_N)
 
! S.2. Types differents / differents Types
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            typ__1  =  max(zer0,sign(un_1,eps6-Buf_G1))   ! =1.=> Dendritic
            zroNEW  =(    typ__1  *d_Bufs                 ! fract.Dendr.Lay.
     .              + (1.-typ__1) *BufsSV(ikl))           !
     .                   /max(eps6,Bufs_N)                !
            G1_NEW  =     typ__1  *Buf_G1                 ! G1 of Dendr.Lay.
     .              + (1.-typ__1) *BG1sSV(ikl)            !
            G2_NEW  =     typ__1  *Buf_G2                 ! G2 of Dendr.Lay.
     .              + (1.-typ__1) *BG2sSV(ikl)            !
            zroOLD  =((1.-typ__1) *d_Bufs                 ! fract.Spher.Lay.
     .              +     typ__1  *BufsSV(ikl))           !
     .                   /max(eps6,Bufs_N)                !
            G1_OLD  = (1.-typ__1) *Buf_G1                 ! G1 of Spher.Lay.
     .              +     typ__1  *BG1sSV(ikl)            !
            G2_OLD  = (1.-typ__1) *Buf_G2                 ! G2 of Spher.Lay.
     .              +     typ__1  *BG2sSV(ikl)            !
            SizNEW  =    -G1_NEW  *DDcdSV/G1_dSV          ! Size  Dendr.Lay.
     .               +(1.+G1_NEW         /G1_dSV)         !
     .                  *(G2_NEW  *DScdSV/G1_dSV          !
     .               +(1.-G2_NEW         /G1_dSV)*DFcdSV) !
            SphNEW  =     G2_NEW         /G1_dSV          ! Spher.Dendr.Lay.
            SizOLD  =     G2_OLD                          ! Size  Spher.Lay.
            SphOLD  =     G1_OLD         /G1_dSV          ! Spher.Spher.Lay.
            Siz_av  =   ( zroNEW  *SizNEW+zroOLD*SizOLD)  ! Averaged Size
            Sph_av = min( zroNEW  *SphNEW+zroOLD*SphOLD   !
     .                  ,   un_1                       )  ! Averaged Sphericity
            Den_av = min((Siz_av  - (    Sph_av *DScdSV   !
     .                              +(1.-Sph_av)*DFcdSV)) !
     .                 / (DDcdSV  - (    Sph_av *DScdSV   !
     .                              +(1.-Sph_av)*DFcdSV)) !
     .                  ,   un_1                         )!
            DendOK  = max(zer0,                           !
     .                    sign(un_1,     Sph_av *DScdSV   ! Small   Grains
     .                              +(1.-Sph_av)*DFcdSV   ! Faceted Grains
     .                              -    Siz_av        )) !
!           REMARQUE: le  type moyen (dendritique ou non) depend
!           ^^^^^^^^  de la  comparaison avec le diametre optique
!                     d'une neige recente de   dendricite nulle
!           REMARK:   the mean type  (dendritic   or not) depends
!           ^^^^^^    on the comparaison with the optical diameter
!                     of a recent snow    having zero dendricity
 
            G1diff  =(   -DendOK *Den_av
     .               +(1.-DendOK)*Sph_av) *G1_dSV
            G2diff  =     DendOK *Sph_av  *G1_dSV
     .               +(1.-DendOK)*Siz_av
            G1      =     SameOK *G1same
     .               +(1.-SameOK)*G1diff
            G2      =     SameOK *G2same
     .               +(1.-SameOK)*G2diff
 
            BG1sSV(ikl) =                       G1      !
     .                  *       Bufs_N/max(eps6,Bufs_N) ! 
            BG2sSV(ikl) =                       G2      !
     .                  *       Bufs_N/max(eps6,Bufs_N) ! 

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #s0.          nn__SV(ikl).EQ.nwr_SV)
! #s0.      write(6,6605) Buf_G1     ,typ__1
! #s0.                   ,DendOK     ,Den_av     ,Sph_av     ,Siz_av
! #s0.                   ,G1same     ,G1diff     ,G1
 6605       format('B1,Typ  : ',2e15.6,11x,'OK,Den,Sph,Siz: ',4e15.6
     .          ,/,'          ',30x   ,11x,'sam,dif,G1    : ',3e15.6)

! Update of Buffer Layer Content & Decision about creating a new snow layer
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            BufsSV(ikl) =       Bufs_N                  !     [mm w.e.]
            NLaysv(ikl) = min(un_1,                     !
     .                    max(zer0,                     ! Allows to create 
     .                        sign(un_1,BufsSV(ikl)     ! a new snow Layer
     .                                 -SMndSV     ))   ! if Buffer > SMndSV
     .                   *max(zer0,                     ! Except if * Erosion
     .                        sign(un_1,half            ! dominates
     .                                 -dsnbSV(ikl)))   !
     .                   +max(zer0,                     ! Allows to create
     .                        sign(un_1,BufsSV(ikl)     ! a new snow Layer
     .                                 -SMndSV*3.00)))  ! is Buffer > SMndSV*3

            Bdzssv(ikl) = 1.e-3*BufsSV(ikl)*rhoWat      ! [mm w.e.] -> [m w.e.]
     .                            /max(eps6,BrosSV(ikl))!& [m w.e.] -> [m]
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #s0.          nn__SV(ikl).EQ.nwr_SV)
! #s0.      write(6,6606) BG1sSV(ikl),BG2sSV(ikl)
! #s0.                   ,NLaysv(ikl),BdzsSV(ikl)  
 6606       format('G1,G2 N*: ',2e15.6,i15,e27.6)

          END DO


! Snow Pack Discretization
! ========================

!            **********
        call SISVAT_zSn
!            **********

!            **********
! #ve   call SISVAT_wEq('_zSn  ',0)
!            **********

! OUTPUT in SISVAT for ikl = 1 (preferably for Stand Alone Version)
! OUTPUT           for SnowFall and Snow Buffer
! #s2   IF          (isnoSV(1) .GT. 0)
! #s2.  write(6,6004)isnoSV(1),    dsn_SV(1) *dt__SV + BufsSV(1),
! #s2.              (dzsnSV(1,isn)*ro__SV(1,isn),isn=1,isnoSV(1))
 6004   format(i3,'  dsn+Buf=',f6.2,6x,'z dz *ro =',10f6.2,
     .                                       (/,35x,10f6.2))


! Add a new Snow Layer
! ====================

          DO ikl=1,knonv
            isnoSV(ikl)     = isnoSV(ikl)         +NLaysv(ikl)
            isn             = isnoSV(ikl)
            dzsnSV(ikl,isn) = dzsnSV(ikl,isn) * (1-NLaysv(ikl))
     .                      + Bdzssv(ikl)     *    NLaysv(ikl) 
            TsisSV(ikl,isn) = TsisSV(ikl,isn) * (1-NLaysv(ikl))
     .                  + min(TaT_SV(ikl),Tf_Sno) *NLaysv(ikl) 
            ro__SV(ikl,isn) = ro__SV(ikl,isn) * (1-NLaysv(ikl))
     .                      + Brossv(ikl)     *    NLaysv(ikl) 
            eta_SV(ikl,isn) = eta_SV(ikl,isn) * (1-NLaysv(ikl))   ! + 0. 
            agsnSV(ikl,isn) = agsnSV(ikl,isn) * (1-NLaysv(ikl))   ! + 0. 
            G1snSV(ikl,isn) = G1snSV(ikl,isn) * (1-NLaysv(ikl))
     .                      + BG1ssv(ikl)     *    NLaysv(ikl)
            G2snSV(ikl,isn) = G2snSV(ikl,isn) * (1-NLaysv(ikl))
     .                      + BG2ssv(ikl)     *    NLaysv(ikl)
            istoSV(ikl,isn) = istoSV(ikl,isn) * (1-NLaysv(ikl))
     .   + max(zer0,sign(un_1,TaT_SV(ikl)
     .                       -Tf_Sno-eps_21)) *    istdSV(2)
     .                                        *    NLaysv(ikl)
            BufsSV(ikl)     = BufsSV(ikl)     * (1-NLaysv(ikl))
            NLaysv(ikl)     = 0
          END DO


! Snow Pack Thickness
! -------------------

          DO ikl=1,knonv
            z_snsv(ikl)     = 0.0
          END DO
        DO   isn=1,nsno
          DO ikl=1,knonv
            z_snsv(ikl)     = z_snsv(ikl) + dzsnSV(ikl,isn)
            zzsnsv(ikl,isn) = z_snsv(ikl)
          END DO
        END DO


! Diffusion of Surficial Water in the Snow Pack
! ---------------------------------------------

c #sw     DO isn=1,nsno
c #sw     DO ikl=1,knonv
c #sw       PorVol      = 1.     - ro__SV(ikl,isn) /  rhoIce          !
c #sw       PorVol      =      max(PorVol            ,zer0  )         !
c #sw       rWater      = ws0dSV * PorVol     *rhoWat*dzsnSV(ikl,isn)
c #sw.                  * max(zer0,
c #sw.                   sign(un_1,rusnSV(ikl)/rhoWat-zzsnsv(ikl,isn)
c #sw.                                               +dzsnSV(ikl,isn)))
c #sw       rusNEW      =      max(rusnSV(ikl)-rWater,zer0  )
c #sw       rWater      =          rusnSV(ikl)-rusNEW
c #sw       rdzNEW          =      rWater
c #sw.                           + ro__SV(ikl,isn) * dzsnSV(ikl,isn)
c #sw       etaNEW          =      rWater / max(eps6,rdzNEW)
c #sw       rusnSV(ikl) =          rusNEW
c #sw       ro__SV(ikl,isn) =      rdzNEW / max(eps6,dzsnSV(ikl,isn))
c #sw       eta_SV(ikl,isn) =      eta_SV(ikl,isn)  +etaNEW
c #sw     ENDDO
c #sw     ENDDO

      END IF

! OUTPUT in SISVAT for ikl = 1 (preferably for Stand Alone Version)
! OUTPUT           for SnowFall and Snow Buffer
! #s2   IF (knonv>0) THEN
! #s2   IF          (isnoSV(1) .GT. 0)
! #s2.  write(6,6006)isnoSV(1),    dsn_SV(1) *dt__SV + BufsSV(1),
! #s2.              (dzsnSV(1,isn)*ro__SV(1,isn),isn=1,isnoSV(1))
! #s2   END IF
 6006   format(i3,'  dsn+Buf=',f6.2,6x,'* dz *ro =',10f6.2,
     .                                       (/,35x,10f6.2))


! Blowing Dust
! ============

c #BD   IF (BloMod)                                               THEN

!         ***************
c #BD     call SISVAT_BDu
!         ***************

c #BD   END IF



! Soil      Albedo: Soil Humidity Correction
! ==========================================

!         REFERENCE: McCumber and Pielke (1981), Pielke (1984)
!         ^^^^^^^^^
          DO ikl=1,knonv
            albssv(ikl) = 
     .      alb0SV(ikl) *(1.0-min(half,eta_SV(       ikl,0)
     .                                /etadSV(isotSV(ikl))))
!         REMARK:    Albedo of Water Surfaces (isotSV=0): 
!         ^^^^^^     alb0SV := 2  X  effective value, while 
!                    eta_SV :=          etadSV
          END DO


! Snow Pack Optical Properties
! ============================

      IF (SnoMod)                                                 THEN

!            ******
        call SnOptP(jjtime)
!            ******

      ELSE
        DO ikl=1,knonv
          sEX_sv(ikl,1) = 1.0
          sEX_sv(ikl,0) = 0.0
          albisv(ikl)   = albssv(ikl)
        END DO
      END IF

!            **********
! #ve   call SISVAT_wEq('SnOptP',0)
!            **********


! Solar Radiation Absorption and Effective Leaf Area Index
! ========================================================

!            ******
        call VgOptP
!            ******


! Surface-Canopy Emissivity
! =========================

        DO ikl=1,knonv
            LSnMsk     =     min( 1,isnoSV(ikl))
            tau_sv(ikl)=     exp(  -LAI_sv(ikl))           ! Veg Transmit.Frac.
            Evg_sv(ikl)=  EmiVeg*(1-LSnMsk)+EmiSno*LSnMsk  ! Veg+Sno Emissivity
            Eso_sv(ikl)=  EmiSol*(1-LSnMsk)+EmiSno*LSnMsk  ! Sol+Sno Emissivity
            emi_SV(ikl)=
     .   (((EmiSol*     tau_sv(ikl) 
     .     +EmiVeg*(1.0-tau_sv(ikl))) *LSmask(ikl))
     .    + EmiWat                 *(1-LSmask(ikl)))*(1-LSnMsk)
     .   +  EmiSno                                     *LSnMsk
        END DO


! Soil/Vegetation Forcing/ Upward IR (INPUT, from previous time step)
! ===================================================================

        DO ikl=1,knonv
! #e1     Enrsvd(ikl) =    - IRs_SV(ikl)
          IRupsv(ikl) =      IRs_SV(ikl) *     tau_sv(ikl) ! Upward   IR
        END DO


! Turbulence
! ==========

! Latent Heat of Vaporization/Sublimation
! ---------------------------------------

        DO ikl=1,knonv
          SnoWat      =                     min(isnoSV(ikl),0)
          Lx_H2O(ikl) =
     .    (1.-SnoWat) * LhvH2O 
     .  +     SnoWat  *(LhsH2O * (1.-eta_SV(ikl,isnoSV(ikl)))
     .                 +LhvH2O *     eta_SV(ikl,isnoSV(ikl)) )
        END DO


! Roughness Length for Momentum
! -----------------------------

! Land+Sea-Ice / Ice-free Sea Mask
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,knonv
          IcIndx(ikl) = 0
        ENDDO
        DO isn=1,nsno
        DO ikl=1,knonv
          IcIndx(ikl) = max(IcIndx(ikl),
     .                      isn*max(0,
     .                              sign(1,
     .                                   int(ro__SV(ikl,isn)-900.))))
        ENDDO
        ENDDO

        DO ikl=1,knonv
          LISmsk    =     min(iiceSV(ikl),1     )
          LISmsk    =     max(LSmask(ikl),LISmsk)
          IceMsk    =     max(0,sign(1   ,IcIndx(ikl)-1)  )
          SnoMsk    = max(min(isnoSV(ikl)-iiceSV(ikl),1),0)

! Sea  Roughness Length
! ^^^^^^^^^^^^^^^^^^^^^
          Z0mSea =       0.0002
          Z0hSea =       0.000049

c #zs     Z0mSea =       0.0185*us__SV(ikl)*us__SV(ikl) ! Doyle MWR 130
c #zs.                         *Grav_I                  ! p.3088 2e col 

c #ZS     Z0mSea =       0.016 *us__SV(ikl)*us__SV(ikl) ! Wang  MWR 129
c #ZS.                         *Grav_I                  ! p.1377 (21)   
c #ZS.           +       0.11  *A_MolV                  !
c #ZS.                         /  max(eps6 ,us__SV(ikl))!

c #zs     Z0mSea =       0.0185*us__SV(ikl)*us__SV(ikl) ! Wang  MWR 129
c #zs.                         *Grav_I                  ! p.1377 (21)   
c #zs.           +       0.135 *A_MolV                  !   (adapted)
c #zs.                         /  max(eps6 ,us__SV(ikl))!

c #ZS     Z0hSea =   max(0.000049,                      ! Wang  MWR 129
c #ZS.                   0.20  *A_MolV                  ! p.1377 (22)   
c #ZS.                         /  max(eps6 ,us__SV(ikl)))

c #ZS     Z0mSea =   max(Z0mSea,eps6)                   ! 

! Land Roughness Length, Snow Contribution excluded
! ^^^^^^^^^^^^^^^^^^^^^^ Ice  Contribution included
!                        ^^^^^^^^^^^^^^^^^^^^^^^^^^
! If vegetation Seasonal Cycle described by  LAI     :
          growth      =min(max(0,7-ivgtSV(ikl)),1)
          Z0mLnd      =     Z0mdSV(ivgtSV(ikl))*LAI_sv(ikl)*growth
     .                                         /LAIdSV
     .                +     Z0mdSV(ivgtSV(ikl))*         (1-growth)

! If vegetation Seasonal Cycle described by  GLF only:
c #CM     Z0mLnd      =     
c #CM.             fallen * Z0mLnd     
c #CM.        +(1.-fallen)* Z0mdSV(ivgtSV(ikl))*glf_sv(ikl)*growth
c #CM.                 +    Z0mdSV(ivgtSV(ikl))*         (1-growth)

! Land Roughness Length, Influence of the Masking by Snow
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0mLnd      =max( Z0mLnd   ,
     .                      Z0mdSV(0)*(1-IceMsk)
     .                     +Z0_ICE   *   IceMsk )
          Z0mLnd      =     Z0mLnd             
     .                    -(zzsnsv(ikl,    isnoSV(ikl))
     .                     -zzsnsv(ikl,max(IcIndx(ikl),0)))/7. 
          Z0mLnd      =max( Z0mLnd    ,    5.e-5  )  ! Min set := Z0 on *
!         Roughness  disappears under Snow 
!         Assumption Height/Roughness Length =  7 is used

! Z0 Smooth Regime over Snow (Andreas 1995, CRREL Report 95-16, p. 8)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0m_nu =       5.e-5 ! z0s~(10-d)*exp(-vonKrm/sqrt(1.1e-03))  

! Z0 Saltat.Regime over Snow (Gallee  et al., 2001, BLM 99 (19) p.11)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          u2star =       us__SV(ikl) *us__SV(ikl)
          Z0mBSn =       u2star      *0.536e-3   -  61.8e-6
          Z0mBSn =   max(Z0mBS0      ,Z0mBSn)

! Z0 Smooth + Saltat. Regime
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0enSV(ikl) =  Z0m_nu
     .                +  Z0mBSn

! Rough   Snow Surface Roughness Length (Typical Value)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0m_Sn =      25.e-5 ! Andreas 1995, CRREL Report 95-16, fig.1&p.2
                               ! z0r~(10-d)*exp(-vonKrm/sqrt(1.5e-03))-5.e-5
c #TZ     Z0m_Sn =       5.e-5 ! Typical Tuning in polar mesoscale models

! Rough   Snow Surface Roughness Length (Variable Sastrugi Height)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          A_Fact      =  1.0000        ! Andreas et al., 2004, p.4
                                       ! ams.confex.com/ams/pdfpapers/68601.pdf

c #SZ     Z0Sa_N =                   (us__SV(ikl) -0.2)*0.0050 ! 0050=0.003/0.6
c #SZ     Z0SaSi = max(zer0,sign(un_1,Z0Sa_N))                 ! 0100=TUNING
c #SZ     Z0Sa_N = max(zer0,          Z0Sa_N)
c #SZ     Z0SaSV(ikl) =
c #SZ.             max(Z0SaSV(ikl)   ,Z0SaSV(ikl) 
c #SZ.               + Z0SaSi*(Z0Sa_N-Z0SaSV(ikl))*exp(-dt__SV/43200.))
c #SZ.               -            min(dz0_SV(ikl) ,     Z0SaSV(ikl))

c #SZ     A_Fact      =               Z0SaSV(ikl) *  5.0/0.15  ! A=5 if h~10cm
!         CAUTION: The influence of the sastrugi direction is not yet included

c #SZ     Z0m_Sn =                    Z0SaSV(ikl)              !
c #SZ.                              - Z0m_nu                   !

! Z0                         (Shao & Lin, 1999, BLM 91 (46)  p.222)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
! Z0 Saltat.Regime over Snow (Shao & Lin, 1999, BLM 91 (46)  p.222)
c #ZN     sqrrZ0 =       usthSV(ikl)/max( us__SV(ikl),0.001)
c #ZN     sqrrZ0 =                   min( sqrrZ0     ,0.999)
c #ZN     Z0mBSn =       0.55 *0.55 *exp(-sqrrZ0     *sqrrZ0)
c #ZN.                  *us__SV(ikl)*     us__SV(ikl)*Grav_I*0.5

! Z0 Smooth + Saltat. Regime (Shao & Lin, 1999, BLM 91 (46)  p.222)
c #ZN     Z0enSV(ikl) = (Z0m_nu     **    sqrrZ0 )
c #ZN.                * (Z0mBSn     **(1.-sqrrZ0))
c #ZN     Z0enSV(ikl) =  max(Z0enSV(ikl), Z0m_nu)

! Z0                         (Andreas etAl., 2004
! ^^^^^^^^^^^^^^^^^^^^^^^^^^  ams.confex.com/ams/pdfpapers/68601.pdf)
! Z0 Smooth Regime over Snow (Andreas etAl., 2004
c #ZA     Z0m_nu = 0.135*A_MolV / max(us__SV(ikl) , eps6)

! Z0 Saltat.Regime over Snow (Andreas etAl., 2004
c #ZA     Z0mBSn = 0.035*u2star      *Grav_I

! Z0 Smooth + Saltat. Regime (Andreas etAl., 2004
c #ZA     Z0enSV(ikl) =  Z0m_nu
c #ZA.                +  Z0mBSn

! Z0 Rough  Regime over Snow (Andreas etAl., 2004
! (.NOT. used by Erosion)     ams.confex.com/ams/pdfpapers/68601.pdf)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
c #ZA     Z0m_90 =(10.-0.025*VVs_SV(ikl)/5.)
c #ZA.            *exp(-0.4/sqrt(.0020+.00001*max(0.,VVs_SV(ikl)-5.))) 
c #ZA     Z0m_Sn =           DDs_SV(ikl)* Z0m_90 / 45.
c #ZA.         - DDs_SV(ikl)*DDs_SV(ikl)* Z0m_90 /(90.*90.)
! #za     u2star =      (us__SV(ikl) -0.1800)     / 0.1
! #za     Z0m_Sn =A_Fact*Z0mBSn *exp(-u2star*u2star)

! Z0 Rough  Regime over Snow (Andreas etAl., 2004
c #ZA     u2star =      (us__SV(ikl) -0.1800)     / 0.1
c #ZA     Z0m_Sn =A_Fact*Z0mBSn *exp(-u2star*u2star)

! Z0 Smooth + Saltat. Regime + Rough  Regime over Snow (Andreas etAl., 2004)
c #ZA     Z0enSV(ikl) =  Z0enSV(ikl)
c #ZA.                +  Z0m_Sn

! Z0               over Snow (instantaneous or time average)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0e_SV(ikl) =  Z0enSV(ikl)
c #ZM     Z0e_SV(ikl) =  Z0emSV(ikl)

! Momentum  Roughness Length
! ^^^^^^^^^^^^^^^^^^^^^^^^^^                                 ! Contribution of
          Z0mnSV(ikl) =  Z0mLnd                              ! Vegetation Form
     .                + (Z0m_Sn                              ! Sastrugi   Form
     .                +  Z0enSV(ikl))   *SnoMsk              ! Snow    Erosion

! Mom. Roughness Length, Discrimination among Ice/Land  and Ice-Free Ocean
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0mnSV(ikl) =  Z0mnSV(ikl)    *LISmsk              ! Ice and  Land
     .                  +Z0mSea      *(1-LISmsk)             ! Ice-Free Ocean
c #OR.                  +Z0roSV(ikl)                         ! Subgrid  Topogr.

! GIS  Roughness Length
! ^^^^^^^^^^^^^^^^^^^^^
c #GL     Z0mnSV(ikl) = 
c #GL.      (1-LSmask(ikl)) *     Z0mnSV(ikl)
c #GL.    +    LSmask(ikl)  * max(Z0mnSV(ikl),max(Z0_GIM,
c #GL.                                            Z0_GIM+
c #GL.      (0.0032-Z0_GIM)*(ro__SV(ikl,isnoSV(ikl))-600.)   !
c #GL.                     /(920.00                 -600.))) !

! Mom. Roughness Length, Instantaneous OR Box Moving Average in Time
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0m_SV(ikl) =  Z0mnSV(ikl)                         ! Z0mnSV  instant.
c #ZM     Z0m_SV(ikl) =  Z0mmSV(ikl)                         ! Z0mnSV  Average

! Corrected Threshold Friction Velocity before Erosion    ! Marticorena and
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ! Bergametti 1995
c #BS     Z0e_SV(ikl) =   min(Z0m_SV(ikl),Z0e_SV(ikl))       !
c #MB     f_eff=    log(0.35*(0.1        /Z0e_SV(ikl))**0.8) ! JGR 100
c #MB     f_eff=1.-(log(      Z0m_SV(ikl)/Z0e_SV(ikl)      ))! (20) p. 16420
c #MB.            /(max(      f_eff      ,eps6             ))! p.16426 2nd ?
c #MB     f_eff=    max(      f_eff      ,eps6              )! CONTROL
! #mB     f_eff=2.0*max(      f_eff      ,eps6              )! TUNING
c #MB     f_eff=    min(      f_eff      ,un_1              )!
c #MB     usthSV(ikl) =       usthSV(ikl)/f_eff              !


! Roughness Length for Scalars
! ----------------------------

          Z0hnSV(ikl) =     Z0mnSV(ikl)/  7.4
c #SH     Z0hnSV(ikl) =     Z0mnSV(ikl)/100.0
!                           Z0h = Z0m  /100.0   over the Sahel
!                                              (Taylor & Clark, QJRMS 127,p864)

c #RN     rstar       =     Z0mnSV(ikl) * us__SV(ikl) / A_MolV
c #RN     rstar       = max(eps6,min(rstar,R_1000))
c #RN     alors       =          log(rstar)   
c #RN     rstar0      = 1.250e0 * max(zer0,sign(un_1,0.135e0 - rstar))
c #RN.                +(1.      - max(zer0,sign(un_1,0.135e0 - rstar)))
c #RN.                *(0.149e0 * max(zer0,sign(un_1,2.500e0 - rstar)) 
c #RN.                + 0.317e0 
c #RN.                *(1.      - max(zer0,sign(un_1,2.500e0 - rstar))))      
c #RN     rstar1      = 0.      * max(zer0,sign(un_1,0.135e0 - rstar))
c #RN.                +(1.      - max(zer0,sign(un_1,0.135e0 - rstar)))
c #RN.                *(-0.55e0 * max(zer0,sign(un_1,2.500e0 - rstar)) 
c #RN.                - 0.565
c #RN.                *(1.      - max(zer0,sign(un_1,2.500e0 - rstar))))
c #RN     rstar2      = 0.      * max(zer0,sign(un_1,0.135e0 - rstar))
c #RN.                +(1.      - max(zer0,sign(un_1,0.135e0 - rstar)))
c #RN.                *(0.      * max(zer0,sign(un_1,2.500e0 - rstar)) 
c #RN.                - 0.183
c #RN.                *(1.00    - max(zer0,sign(un_1,2.500e0 - rstar))))
c #RN     Z0hnSV(ikl) = max(zer0
c #RN.                , sign(un_1,zzsnsv(ikl,isnoSV(ikl))-eps6)) 
c #RN.                * exp(rstar0+rstar1*alors+rstar2*alors*alors)
c #RN.                * 0.001e0 + Z0hnSV(ikl) * ( 1. - max(zer0   
c #RN.                , sign(un_1,zzsnsv(ikl,isnoSV(ikl))-eps6)))

          Z0hnSV(ikl) =     Z0hSea             *(1-LISmsk) ! Ice-free Ocean
     .                +     Z0hnSV(ikl)        *   LISmsk  ! Ice and  Land

          Z0h_SV(ikl) =     Z0hnSV(ikl)
c #ZM     Z0h_SV(ikl) =     Z0hmSV(ikl)


! Contributions of the Roughness Lenghths to the neutral Drag Coefficient
! -----------------------------------------------------------------------

c #MT     Z0m_SV(ikl) = max(2.0e-6     ,Z0m_SV(ikl)) ! Min Z0_m (Garrat Scheme)
          sqrCm0(ikl) = log(za__SV(ikl)/Z0m_SV(ikl))
          sqrCh0(ikl) = log(za__SV(ikl)/Z0h_SV(ikl))

! OUTPUT of SnowFall, Roughness Length and Drag Coefficients
! #sf     IF (ikl.EQ.1) write(6,6661) dsn_SV(ikl),us__SV(ikl),Z0SaSi
! #sf.                        ,Z0Sa_N,Z0SaSV(ikl),Z0m_Sn,Z0m_SV(ikl)
 6661     format(20x,7f9.6)

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           of Roughness Length and Drag Coefficients
! #sz     IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #sz.        nn__SV(ikl).EQ.nwr_SV)
! #sz.    write(6,6600)  za__SV(ikl) , Z0m_SV(ikl)
! #sz.                  ,sqrCm0(ikl) , za__SV(ikl)/Z0m_SV(ikl)
! #sz.                  ,Z0SaSV(ikl) , Z0h_SV(ikl)
! #sz.                  ,sqrCh0(ikl) , za__SV(ikl)/Z0h_SV(ikl)
 6600     format(/,' ** SISVAT     *0  '
     .            ,'  za__SV  = ',e12.4,'  Z0m_SV  = ',e12.4
     .            ,'  sqrCm0  = ',e12.4,'  Za/Z0m  = ',e12.4
     .          ,/,'                   '
     .            ,'  Z0SaSV  = ',e12.4,'  Z0h_SV  = ',e12.4
     .            ,'  sqrCh0  = ',e12.4,'  Za/Z0h  = ',e12.4)


! Vertical Stability Correction
! -----------------------------

! Surface/Canopy Temperature
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          Tsrfsv(ikl) = Sigmsv(ikl) * TvegSV(ikl)
     .          + (1. - Sigmsv(ikl))* TsisSV(ikl,isnoSV(ikl))
        END DO

! Aerodynamic Resistance
! ^^^^^^^^^^^^^^^^^^^^^^
c #CP   IF            (SnoMod.AND.ColPrt)                           THEN

!                **********
c #CP       call ColPrt_SBL
!                **********

c #CP   ELSE
c #MT    IF           (Garrat)                                      THEN

!                **********
c #MT       call SISVAT_SBL
!                **********

c #MT    ELSE

!                **********
            call SISVATeSBL
!                **********

c #MT    END IF
c #CP   END IF


! Canopy Energy Balance
! =====================

!            **********
        call SISVAT_TVg
! #e1.                 (ETVg_d)
!            **********


! Surface/Canopy Temperature
! ==========================

        DO ikl=1,knonv
          Tsrfsv(ikl) = Sigmsv(ikl) * TvegSV(ikl)
     .          + (1. - Sigmsv(ikl))* TsisSV(ikl,isnoSV(ikl))
        END DO


! Soil   Energy Balance
! =====================

!            **********
        call SISVAT_TSo
! #e1.                 (ETSo_0,ETSo_1,ETSo_d)
!            **********

!            **********
! #ve   call SISVAT_wEq('_TSo  ',0)
!            **********


! Canopy Water  Balance
! =====================

! Soil Water     Potential
! ------------------------

      DO   isl=-nsol,0
        DO ikl=1,knonv
          ist             =     isotSV(ikl)        ! Soil Type
          psi_sv(ikl,isl) =     psidSV(ist)        ! DR97, Eqn.(3.34)
     .  *(etadSV(ist) /max(eps6,eta_SV(ikl,isl)))  !
     .  **bCHdSV(ist)                              !


! Soil Hydraulic Conductivity
! ---------------------------

          Khydsv(ikl,isl) =    s2__SV(ist)         ! DR97, Eqn.(3.35)
     .  *(eta_SV(ikl,isl)**(2.*bCHdSV(ist)+3.))    !  
        END DO
      END DO

!            **********
        call SISVAT_qVg
!            **********


! OUTPUT/Verification: H2O    Conservation: Vegetation Forcing
! #m0   DO ikl=1,knonv
! #m0     Watsvd(ikl) =     (Watsvd(ikl)           ! Canopy Precip. IN
! #m0.                      -drr_SV(ikl)           ! Canopy Precip. OUT
! #m0.                      -Evp_sv(ikl))* dt__SV  ! Canopy Water Evap.
! #m0   END DO


! Melting / Refreezing in the Snow Pack
! =====================================

      IF (SnoMod)                                                 THEN

!            **********
        call SISVAT_qSn
     .                 (
! #e1.                  EqSn_0,EqSn_1,EqSn_d
! #m1.                 ,SIsubl,SImelt,SIrnof
     .                 )
!            **********

!            **********
! #ve   call SISVAT_wEq('_qSn  ',0)
!            **********

! OUTPUT in SISVAT for ikl = 1 (preferably for Stand Alone Version)
! OUTPUT           for SnowFall and Snow Buffer
! #s2   IF (knonv>0) THEN
! #s2   IF          (isnoSV(1) .GT. 0)
! #s2.  write(6,6007)isnoSV(1),    dsn_SV(1) *dt__SV + BufsSV(1),
! #s2.              (dzsnSV(1,isn)*ro__SV(1,isn),isn=1,isnoSV(1))
! #s2   END IF
 6007   format(i3,'  dsn+Buf=',f6.2,6x,'q dz *ro =',10f6.2,
     .                                       (/,35x,10f6.2))


! Snow Pack Thickness
! -------------------

          DO ikl=1,knonv
            z_snsv(ikl)     = 0.0
          END DO
        DO   isn=1,nsno
          DO ikl=1,knonv
            z_snsv(ikl)     = z_snsv(ikl) + dzsnSV(ikl,isn)
            zzsnsv(ikl,isn) = z_snsv(ikl)
          END DO
        END DO


! Energy in Excess is added to the first Soil Layer
! -------------------------------------------------

        DO ikl=1,knonv
            z_snsv(ikl)   = max(zer0,
     .                          sign(un_1,eps6-z_snsv(ikl)))
            TsisSV(ikl,0) = TsisSV(ikl,0)    + EExcsv(ikl) 
     .                                       /(rocsSV(isotSV(ikl))
     .                                        +rcwdSV*eta_SV(ikl,0))
            EExcsv(ikl)   = 0.
        END DO


! OUTPUT/Verification: * Mass Conservation: Mass (below the Canopy) and Forcing
! #m1   DO ikl=1,knonv
! #m1     SIWa_f(ikl) =(drr_SV(ikl) + dsn_SV(ikl))   *dt__SV         ![mm w.e.]
! #m1     SIWe_f(ikl) = dbs_SV(ikl)                                  !
! #m1     SIWm_1(ikl) = BufsSV(ikl) + HFraSV(ikl)    *rhoIce         !
! #m1   DO isn=1,nsno                                                !
! #m1     SIWm_1(ikl) = SIWm_1(ikl) + dzsnSV(ikl,isn)*ro__SV(ikl,isn)!
! #m1   END DO                                                       !
! #m1   END DO                                                       !

      END IF


! Soil   Water  Balance
! =====================

!            **********
        call SISVAT_qSo
! #m0.                 (Wats_0,Wats_1,Wats_d)
!            **********


! Surface/Canopy Fluxes
! =====================

        DO ikl=1,knonv
          IRdwsv(ikl)=tau_sv(ikl) *IRd_SV(ikl)*Eso_sv(ikl) ! Downward IR
     .          +(1.0-tau_sv(ikl))*IRd_SV(ikl)*Evg_sv(ikl) ! 
          IRupsv(ikl) =      IRupsv(ikl)                   ! Upward   IR
     .                + 0.5 *IRv_sv(ikl) * (1.-tau_sv(ikl))!
          IRu_SV(ikl) =     -IRupsv(ikl)                   ! Upward   IR
     .                      +IRd_SV(ikl)                   ! (effective)
     .                      -IRdwsv(ikl)                   ! (positive)

          TBr_sv(ikl) =sqrt(sqrt(IRu_SV(ikl)/StefBo))      ! Brightness
!                                                          ! Temperature
          uts_SV(ikl) =     (HSv_sv(ikl) +HSs_sv(ikl))     ! u*T*
     .                     /(rhT_SV(ikl) *CpdAir)          !
          uqs_SV(ikl) =     (HLv_sv(ikl) +HLs_sv(ikl))     ! u*q*
     .                     /(rhT_SV(ikl) *LhvH2O)          !

! Surface/Canopy Temperature
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          Tsrfsv(ikl) = Sigmsv(ikl) * TvegSV(ikl)
     .          + (1. - Sigmsv(ikl))* TsisSV(ikl,isnoSV(ikl))
        END DO


! Snow Pack Properties (sphericity, dendricity, size)
! ===================================================

      IF (SnoMod)                                                 THEN

!            **********
        call SISVAT_GSn
!            **********

!            **********
! #ve   call SISVAT_wEq('_GSn  ',0)
!            **********


! Surficial Water Freezing, including that of a Water Surface (isotSV=0)
! ======================================================================


      END IF


! OUTPUT
! ======

c #E0   DO ikl=1,knonv
c #E0   IF (lwriSV(ikl).ne.0)                                     THEN
c #E0           noUNIT =  no__SV(lwriSV(ikl))
c #E0     write(noUNIT,5001)
c #E0.       (SoSosv(ikl)+SoCasv(ikl))*sol_SV(ikl),
c #E0.        IRdwsv(ikl),IRu_SV(ikl),
c #E0.        HSv_sv(ikl)+HSs_sv(ikl),
c #E0.        HLv_sv(ikl)+HLs_sv(ikl), TaT_SV(ikl),
c #E0.        dsn_SV(ikl)*3.6e3,       drr_SV(ikl)*3.6e3,
c #E0.        SoSosv(ikl)             *sol_SV(ikl),
c #E0.                    IRv_sv(ikl) *0.5,
c #E0.        HSv_sv(ikl),HLv_sv(ikl), TvegSV(ikl),
c #E0.                    SoCasv(ikl) *sol_SV(ikl),
c #E0.        HSs_sv(ikl),HLs_sv(ikl), TsisSV(ikl,isnoSV(ikl)) 
 5001     format(
     .         '        |Net Solar| IR Down | IR Up   | HS/Dwn=+|',
     .          ' HL/Dwn=+| Temper. |         |  Snow  |  Rain  |',
     .       /,'        | [W/m2]  | [W/m2]  | [W/m2]  | [W/m2]  |',
     .          ' [W/m2]  | [K]     |         | [mm/h] | [mm/h] |',
     .       /,' -------+',7('---------+'),2('--------+'),
     .       /,' SISVAT |',f8.1,' |',f8.1,' |',f8.1,' |',f8.1,' |',
     .                     f8.1,' |A',f7.2,' |', 8x ,' |',2(f7.2,' |'),
     .       /,' Canopy |',f8.1,' |', 8x ,' |',f8.1,' |',f8.1,' |',
     .                     f8.1,' |',f8.2,' |', 8x ,' |',2( 7x ,' |')
     .       /,' Soil   |',f8.1,' |', 8x ,' |', 8x ,' |',f8.1,' |',
     .                     f8.1,' |',f8.2,' |', 8x ,' |',2( 7x ,' |'))


! OUTPUT/Verification: Energy/Water Budget
! #e1     Enrsvd(ikl) = Enrsvd(ikl)                            ! Up Surf. IR
! #e1.                + IRs_SV(ikl)                            ! Offset
! #e1.        + (      (SoSosv(ikl)                            ! Net   Solar
! #e1.                 +SoCasv(ikl)) *sol_SV(ikl)              !            
! #e1.          +                     IRdwsv(ikl)              ! Downward IR
! #e1.          +                     IRupsv(ikl)              ! Upward   IR
! #e1.          +                     HSv_sv(ikl)+HSs_sv(ikl)  ! Sensible
! #e1.          +                     HLv_sv(ikl)+HLs_sv(ikl)) ! Latent

! #e1     write(noUNIT,5002)           Enrsvd(ikl),
! #e1.                   ETSo_0(ikl),  ETSo_d(ikl),
! #e1.                   ETSo_0(ikl)+  ETSo_d(ikl), ETSo_1(ikl),
! #e1.                   EqSn_0(ikl)                            /dt__SV,
! #e1.                   EqSn_d(ikl)                            /dt__SV,
! #e1.                  (EqSn_1(ikl)-  EqSn_0(ikl)- EqSn_d(ikl))/dt__SV,
! #e1.                   EqSn_1(ikl)                            /dt__SV
 5002     format(
     .           ' -----------------+-------------------+', !
     .            '-----------------+-+-----------------+', !
     .          '-------------------+',                     !
     .         /,' SOIL/SNOW/VEGET. |                   |', !
     .            ' Power,  Forcing |                   |', ! Enrsvd
     .          '                   |',                     ! 
! #el.         /,' -----------------+-------------------+', !
! #el.            '-----------------+-------------------+', !
! #el.          '-------------------+',                     !
     .         /,'                  |',   11x ,'        |', !
     .                f9.2,' [W/m2] |',   11x ,'        |', ! Enrsvd
     .                11x ,'        |',                     !
     .         /,' -----------------+-------------------+', !
     .            '-----------------+-------------------+', !
     .          '-------------------+',                     ! 
     .         /,' SOIL/SNOW  (TSo) | Energy/dt, Time 0 |', !        ETSo_0
     .            ' Power,  Forcing |   Sum Tim.0+Forc. |', ! ETSo_d/ETSo_0+d
     .          ' Energy/dt, Time 1 |',                     ! ETSo_1
! #el.         /,' -----------------+-------------------+', !
! #el.            '-----------------+-------------------+', !
! #el.          '-------------------+',                     !
     .         /,'                  |',  f11.2,' [W/m2] |', !        ETSo_0
     .                f9.2,' [W/m2] |',  f11.2,' [W/m2] |', ! ETSo_d/ETSo_0+d
     .               f11.2,' [W/m2] |',                     ! ETSo_1
     .         /,' -----------------+-------------------+', !
     .            '-----------------+-------------------+', !
     .          '-------------------+',                     !
     .         /,'      SNOW  (qSn) | Energy/dt, Time 0 |', ! EqSn_0/dt
     .            ' Power,  Excess  |   D(Tim.1-0-Forc.)|', ! EqSn_d/dt, 1-0-d
     .          ' Energy/dt, Time 1 |',                     ! EqSn_1/dt
! #el.         /,' -----------------+-------------------+', !
! #el.            '-----------------+-------------------+', !
! #el.          '-------------------+',                     !
     .         /,'                  |',  f12.2, '[W/m2] |', ! EqSn_0/dt
     .                f9.2,' [W/m2] |',  f11.2,' [W/m2] |', ! EqSn_d/dt, 1-0-d
     .               f12.2, '[W/m2] | ',                    ! EqSn_1/dt
     .         /,' -----------------+-------------------+', !
     .            '-----------------+-------------------+', !
     .          '-------------------+')                     !

! #e1             EnsBal = ETSo_1(ikl)-(ETSo_0(ikl)+Enrsvd(ikl))
! #e1             EnvBal = Enrsvd(ikl)- ETVg_d(ikl)
! #e1     IF (abs(EnsBal).gt.5.e-1 
! #e2.    .OR.lwriSV(ikl).eq.    2 
! #e1.                            )                                 THEN
! #e1       write(6,6001) daHost,i___SV(lwriSV(ikl)),
! #e1.                           j___SV(lwriSV(ikl)),
! #e1.                           n___SV(lwriSV(ikl)),
! #e1.                           ETSo_1(ikl),ETSo_0(ikl),ETSo_d(ikl),
! #e1.                           ETSo_1(ikl)-ETSo_0(ikl)-ETSo_d(ikl),
! #e1.                           Enrsvd(ikl),ETVg_d(ikl),ETSo_d(ikl),
! #e1.                           Enrsvd(ikl)-ETVg_d(ikl)-ETSo_d(ikl)
 6001       format(a18,3i4,' (EB1'           ,f15.6,
     .                ')  - [(EB0           ',f15.6,')',
     .               /,55x,'+(ATM->Snow/Soil',f15.6,')] ',
     .                     '= EBAL'          ,f15.6,' [W/m2]',
     .               /,55x,' (ATM->SISVAT'   ,f18.6,
     .               /,55x,'- Veg. ImBal.',   f18.6,')  ',
     .               /,55x,'- ATM->SnoSol',   f18.6,')  ',
     .                     '= ????'          ,f15.6,' [W/m2]')
! #e1           noEBal = noEBal + 1
! #e2           noEBal = noEBal - 1
! #e1       IF (noEBal.GE.       10) stop 'TOO MUCH ENERGY IMBALANCES'
! #e1     END IF


! OUTPUT/Verification: * Mass Conservation: Budget [mm w.e.]
! #m1     write(noUNIT,5010) 
! #m1.             SIWm_0(ikl),  SIWa_i(ikl)-SIWa_f(ikl)
! #m1.            ,SIWm_0(ikl)+  SIWa_i(ikl)-SIWa_f(ikl)
! #m1.                          +SIWe_i(ikl)-SIWe_f(ikl)
! #m1.                          +SIsubl(ikl)
! #m1.                          -SImelt(ikl)
! #m1.                          -SIrnof(ikl)
! #m2.                          +SIvAcr(ikl)    
! #m1.            ,SIWm_1(ikl),  SIWe_i(ikl)-SIWe_f(ikl)
! #m1.            ,              SIsubl(ikl)
! #m1.            ,             -SImelt(ikl)
! #m1.            ,             -SIrnof(ikl)
! #m2.            ,              SIvAcr(ikl)
 5010     format(' SNOW             |   Snow,   Time 0  |',
     .            ' Snow,   Forcing |           Sum     |',
     .          '   Snow,   Time 1  |',
! #el.         /,' -----------------+-------------------+',
! #el.            '-----------------+-------------------+',
! #el.          '-------------------+',
     .         /,'                  |',    f13.3,' [mm] |',
     .           ' A',  f9.3,' [mm] |',    f13.3,' [mm] |',
     .                 f13.3,' [mm] |',
     .         /,'                  |',     13x ,'      |',
     .           ' E',  f9.3,' [mm] |',     13x ,'      |',
     .                  13x ,'      |',
     .         /,'                  |',     13x ,'      |',
     .           ' S',  f9.3,' [mm] |',     13x ,'      |',
     .                  13x ,'      |',
     .         /,'                  |',     13x ,'      |',
     .           '(M',  f9.3,' [mm])|  (included in A)  |',
     .                  13x ,'      |',
     .         /,'                  |',     13x ,'      |',
     .           ' R',  f9.3,' [mm] |',     13x ,'      |',
     .                  13x ,'      |',
! #m2.         /,'                  |',     13x ,'      |',
! #m2.           ' O',  f9.3,' [mm] |',     13x ,'      |',
! #m2.                  13x ,'      |',
     .         /,' -----------------+-------------------+',
     .            '-----------------+-------------------+',
     .          '-------------------+')
! #m1             SnoBal = SIWm_1(ikl)-(SIWm_0(ikl)
! #m1.                                 +SIWa_i(ikl)-SIWa_f(ikl)
! #m1.                                 +SIWe_i(ikl)-SIWe_f(ikl))
! #m1.                                 -SIsubl(ikl)
! #m1.                                 +SIrnof(ikl)
! #m2.                                 -SIvAcr(ikl) 
! #m1     IF (abs(SnoBal).gt.eps6)                                  THEN
! #m1       write(6,6010) daHost,i___SV(lwriSV(ikl)),
! #m1.                           j___SV(lwriSV(ikl)),
! #m1.                           n___SV(lwriSV(ikl)),
! #m1.                           SIWm_1(ikl),SIWm_0(ikl),
! #m1.                           SIWa_i(ikl),SIWa_f(ikl),
! #m1.                           SIWe_i(ikl),SIWe_f(ikl),
! #m1.                           SIsubl(ikl),SImelt(ikl),
! #m2.                           SIrnof(ikl),SIvAcr(ikl),
! #m1.                           SnoBal
 6010       format(a18,3i4,' (MB1'        ,f12.6,
     .                 ') - [(MB0        ',f12.6,        15x,')',
     .               /,51x,'+(ATM Forcing',f12.6,' - ',f12.6,')',
     .               /,51x,'+(BLS Forcing',f12.6,' - ',f12.6,')',
     .               /,51x,'-(Depo/Sublim',f12.6,        15x,')',
     .               /,51x,' !Melting    ',f12.6,'  included in A!',
     .               /,51x,'+(Run  OFF   ',f12.6,        15x,')',
! #m2.               /,51x,'-(Sea-Ice Acr',f12.6,        15x,')',
     .               /,29x,'= *BAL'       ,f12.6,      ' [mm w.e.]')
! #m1           noSBal = noSBal + 1
! #m1       IF (noSBal.GE.       10) stop 'TOO MUCH SNOW MASS IMBALANCE'
! #m1     END IF


! OUTPUT/Verification: H2O    Conservation: Water  Budget
! #m0     Watsv0(ikl) =  Watsv0(ikl)                    ! Canopy Water Cont.
! #m0.                 + Wats_0(ikl)                    ! Soil   Water Cont.
! #m0     Watsvd(ikl) =  Watsvd(ikl)                    ! Canopy Forcing
! #m0.                 + Wats_d(ikl)                    ! Soil   Forcing

! #m0     write(noUNIT,5003) 
! #m0.                   Wats_0(ikl),  Wats_d(ikl),
! #m0.                   Wats_0(ikl)+  Wats_d(ikl),    Wats_1(ikl),
! #m0.                   Watsv0(ikl),  Watsvd(ikl), 
! #m0.                   Watsv0(ikl)+  Watsvd(ikl),    Wats_1(ikl)
! #m0.                                                +rrCaSV(ikl)
 5003     format(' SOIL/SNOW  (qSo) |   Water,  Time 0  |',
     .            ' Water,  Forcing |           Sum     |',
     .          '   Water,  Time 1  |',
! #el.         /,' -----------------+-------------------+',
! #el.            '-----------------+-------------------+',
! #el.          '-------------------+',
     .         /,'                  |',    f13.3,' [mm] |',
     .                 f11.3,' [mm] |',    f13.3,' [mm] |',
     .                 f13.3,' [mm] |',
     .         /,' -----------------+-------------------+',
     .            '-----------------+-------------------+',
     .          '-------------------+',
     .         /,' SOIL/SNOW/VEGET. |   Water,  Time 0  |',
     .            ' Water,  Forcing |           Sum     |',
     .          '   Water,  Time 1  |',
! #el.         /,' -----------------+-------------------+',
! #el.            '-----------------+-------------------+',
! #el.          '-------------------+',
     .         /,'                  |',    f13.3,' [mm] |',
     .                 f11.3,' [mm] |',    f13.3,' [mm] |',
     .                 f13.3,' [mm] |',
     .         /,' -----------------+-------------------+',
     .            '-----------------+-------------------+',
     .          '-------------------+')

! #m0             WatBal = Wats_1(ikl)+rrCaSV(ikl)
! #m0.                   -(Watsv0(ikl)+Watsvd(ikl))
! #m0     IF (abs(WatBal).gt.eps6)                                  THEN
! #m0       write(6,6002) daHost,i___SV(lwriSV(ikl)),
! #m0.                           j___SV(lwriSV(ikl)),
! #m0.                           n___SV(lwriSV(ikl)),
! #m0.                           Wats_1(ikl),rrCaSV(ikl),
! #m0.                           Watsv0(ikl),Watsvd(ikl),WatBal,
! #m0.                           Wats_1(ikl),
! #m0.                           Wats_0(ikl),Wats_d(ikl),
! #m0.               Wats_1(ikl)-Wats_0(ikl)-Wats_d(ikl)
 6002       format(30x,' NEW Soil Water',3x,' Canopy   Water',3x,
     .                 ' OLD SVAT Water',4x,' FRC SVAT Water',
     .           /,a18,3i4,f15.6,' + ' ,f15.6,' - ' ,f15.6,
     .                           ' -  ',f15.6,'    ', 15x ,'    ',
     .      /,31x,'= ',f12.6,' [mm] (Water Balance)',
     .           /,30x,' NEW Soil Water',3x,'               ',3x,
     .                 ' OLD Soil Water',4x,' FRC Soil Water',
     .           /,30x,f15.6,'   ' , 15x ,' - ' ,f15.6,
     .                       ' -  ',f15.6,'    ', 15x ,'    ',
     .      /,31x,'= ',f12.6,' [mm] (3 terms SUM)')
! #m0           noWBal = noWBal + 1
! #m0       IF (noWBal.GE.       10) stop 'TOO MUCH WATER  IMBALANCES'
! #m0     END IF


! Water/Temperature Profiles
! --------------------------

c #E0       write(noUNIT,5004)
 5004       format(' -----+--------+--+-----+--------+----+---+',
     .  '--------+----+---+--------+------+-+--------+--------+',
     .           /,'    n |     z  |     dz |     ro |    eta |',
     .  '     T  |     G1 |     G2 | Extinc |        | HISTORY|',
     .           /,'      |    [m] |    [m] | [kg/m3]| [m3/m3]|',
     .  '    [K] |    [-] |    [-] |    [-] |        |   [-]  |',
     .           /,' -----+--------+--------+--------+--------+',
     .  '--------+--------+--------+--------+--------+--------+')
c #E0       write(noUNIT,5005) rusnSV(ikl),albisv(ikl)
 5005       format('      |        |        |        |W',f6.3,' |',
     .  '        |        |        |A',f6.3,' |        |        |')
c #E0       write(noUNIT,5015)
c #E0.                    (isn,zzsnsv(ikl,isn),dzsnSV(ikl,isn),
c #E0.                         ro__SV(ikl,isn),eta_SV(ikl,isn),
c #E0.                         TsisSV(ikl,isn),
c #E0.                         G1snSV(ikl,isn),G2snSV(ikl,isn),
c #E0.                         sEX_sv(ikl,isn),istoSV(ikl,isn),
c #E0.                     isn=isnoSV(ikl),1,-1)
 5015       format((i5,' |',2(f7.3,' |'),               f7.1,' |',
     .           f7.3,' |' ,  f7.2,' |', 2(f7.1,' |'),  f7.3,' |',
     .            7x ,' |' ,  i5,'   |'                          ))
c #E0       write(noUNIT,5006)
 5006       format(' -----+--------+--------+--------+--------+',
     .  '--------+--------+--------+--------+--------+--------+')
c #E0       write(noUNIT,5007) TBr_sv(ikl),
c #E0.                         TvegSV(ikl),rrCaSV(ikl)*1.e3,
c #E0.                         EvT_sv(ikl)*86.4e3
 5007       format(' Brgh |',4(8x,'|'),  f7.2,' | [micm] |',4(8x,'|'),
     .           /,' VEGE |',4(8x,'|'),2(f7.2,' |'),        2(8x,'|'),
     .                                   f7.3,' |',           8x,'|' )
c #E0       write(noUNIT,5014)
 5014       format(' -----+--------+--------+--------+--------+',
     .  '--------+--------+--------+--------+--------+--------+',
     .           /,'    n |        |     dz |        |    eta |',
     .  '     T  |        |        |        | Root W.| W.Flow |',
     .           /,'      |        |    [m] |        | [m3/m3]|',
     .  '    [K] |        |        |        | [mm/d] | [mm/h] |',
     .           /,' -----+--------+--------+--------+--------+',
     .  '--------+--------+--------+--------+--------+--------+')

c #E0       write(noUNIT,5008)
c #E0.                    (isl,    LSdzsv(ikl)*dz_dSV(    isl),
c #E0.                                         eta_SV(ikl,isl),
c #E0.                                         TsisSV(ikl,isl),
c #E0.                                  86.4e3*Rootsv(ikl,isl),
c #E0.                                   3.6e3*Khydsv(ikl,isl),
c #E0.                     isl=0,-nsol,-1)
 5008       format((i5,' |',   7x ,' |' ,  f7.3,' |' ,   7x ,' |',
     .           f7.3,' |' ,  f7.2,' |', 2( 7x ,' |'),   7x ,' |',
     .                                     f7.3,' |' ,  f7.2,' |'))
c #E0       write(noUNIT,5006)
c #E0       write(noUNIT,5009) RnofSV(ikl)* 3.6e3
 5009       format('      |',9(8x,'|'),f7.3,' |')
c #E0       write(noUNIT,5006)
c #E0   END IF
c #E0   END DO

! END  .main. (SISVAT)
      return
      end


      subroutine SISVAT_BSn(BloMod)

!--------------------------------------------------------------------------+
!   MAR          SISVAT_BSn                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_BSn treats Snow Erosion and Deposition             |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Possibility                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                          |
!     #BS: Explicit Cloud MICROPHYSICS: Blow. *(Snow)         Model        |
!     #EM: Explicit Cloud MICROPHYSICS: de Montmollin Parameterizat.       |
!     #MA: SNOW Model: Increased polar B* Mobility (Mann et al.2000)       |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #BA: Budd et al.            1966, Ant.Res.Ser.9    u* BS Threshold   |
!     #BY: Budd et al. 1966, 2~m Averag Blow. *(Snow)    Properties        |
!     #AG: Snow Aging Col de Porte     (Brun et al.1991) discard BS at CdP |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # stdout               | #s2: OUTPUT of SnowFall, Snow Buffer          |
!                          |      unit  6, SubRoutine  SISVAT_BSn, _qSn    |
!   # stdout               | #b0: OUTPUT of Snow Erosion                   |
!                          |      unit  6, SubRoutine  SISVAT_BSn **ONLY** |
!--------------------------------------------------------------------------+


       


! General Variables
! =================

      USE PHY_SV

      USE VAR_SV
      USE VARdSV
      USE VARxSV
      USE VARySV

      IMPLICIT NONE

      logical  BloMod


! Local Variables
! ===============

      logical           BlowIn
      common/llocal_BSn/BlowIn

      real              FacSBS,FacUBS        !
      real              Por_BS               ! Snow       Porosity
      real              SheaBS               !
      real              rCd10n               ! GM97:   assumed neutral stabil.
      common/rlocal_BSn/FacSBS,FacUBS,       !
     .                  Por_BS,SheaBS,rCd10n !

      integer  ikl   ,isn   ,isnMAX,is2      !
      integer  Mobilm,Mobiln                 !
      integer  Mobile(klonv)                 !

      real     DendOK                        ! Dendricity Switch
      real     SaltOK                        ! Saltation  Switch
      real     SnowOK                        ! Pack Top   Switch
      real     SaltM1,SaltM2,SaltMo,SaltMx   ! Saltation  Parameters
      real     ShearX                        ! Arg. Max Shear Stress
      real     SaltSU,Salt_U                 !
      real     ArgFac,Fac_Mo                 !
      real     FacRBS,FacTBS                 !
      real     ArguSi                        !
      real     SaltSI(klonv,nsno)            ! Snow Drift Index
      real     hdrift                        ! Inverse erodibl.Snow Lay.Thickn.
      real     h_mmWE                        ! Eroded Snow Layer Min Thickness
      real     tfv_vk                        ! * Fall Veloc. / Von Karman Cst
      real     sdrift(klonv,nsno)            !
      real     xdrift(klonv)                 !
      real     zdrift(klonv)                 !
      real     tdepos(klonv)                 !
      real     zdepos(klonv,nsno)            !
      real     dbsaux(klonv)                 ! Drift Amount   (Dummy Variable)
      real     dzweqo,dzweqn,bsno_x          ! 
      real                   hsno_x          !
      real     PorSno,Salt_f,PorRef,ro_new   !
      real     MIN_Mo                        ! Minimum Mobility Fresh Fallen *
      real     AgBlow                        ! Snow Mobility    Time  Scale
      real     snofOK                        ! Threshd Snow Fall

      integer  isagr1(klonv)                 ! 1st     Layer History
      integer  isagr2(klonv)                 ! 2nd     Layer History

      real     WEagre(klonv)                 ! Snow Water Equivalent Thickness
      real     Agrege(klonv)                 ! 1. when Agregation constrained
      real     dzagr1(klonv)                 ! 1st     Layer Thickness
      real     dzagr2(klonv)                 ! 2nd     Layer Thickness
      real     T_agr1(klonv)                 ! 1st     Layer Temperature
      real     T_agr2(klonv)                 ! 2nd     Layer Temperature
      real     roagr1(klonv)                 ! 1st     Layer Density
      real     roagr2(klonv)                 ! 2nd     Layer Density
      real     etagr1(klonv)                 ! 1st     Layer Water Content
      real     etagr2(klonv)                 ! 2nd     Layer Water Content
      real     G1agr1(klonv)                 ! 1st     Layer Dendricity/Spher.
      real     G1agr2(klonv)                 ! 2nd     Layer Dendricity/Spher.
      real     G2agr1(klonv)                 ! 1st     Layer Sphericity/Size
      real     G2agr2(klonv)                 ! 2nd     Layer Sphericity/Size
      real     agagr1(klonv)                 ! 1st     Layer Age
      real     agagr2(klonv)                 ! 2nd     Layer Age

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0 real     Sno0WE,Sno1WE                 ! Snow Mass before/after Erosion
! #b0 real     SnodWE                        ! Snow Mass              Erosion


! DATA
! ====

      data     AgBlow  / 1.00     /          ! 1 Day (F.Domine, pers.communic.)
      data     SaltMx  /-5.83e-2  /          !
      data     FacRBS  / 2.868    /          !
      data     FacTBS  / 0.085    /          !
      data     hdrift  / 1.00e+1  /          ! Inverse erodibl.Snow Lay.Thickn.
      data     h_mmWE  / 0.01e00  /          ! Eroded Snow Layer Min Thickness
      data     tfv_vk  / 5.10e-1  /          ! tfv (Terminal Fall Veloc. =.216) 
                                             ! /vk (Von Karman Constant  =.4  )
                                             ! (Wamser & Lykosov,   1995 
                                             !  Contr.Atm.Phys. 68, p.90)

! Initialization
! ==============

      IF (.NOT.BlowIn)                                             THEN
               BlowIn = .true.
               FacSBS =  1.             / FacRBS
               FacUBS =  1.             / FacTBS
               Por_BS =  1.             - BSnoRo/      rhoIce
               SheaBS =                   Por_BS/(1.00-Por_BS)
!              SheaBS =  Arg(sqrt(shear = max shear stress in snow)):
!              shear  =  3.420d00 * exp(-(Por_BS      +Por_BS)
!    .                                  /(1.00        -Por_BS))
!              SheaBS :  see de Montmollin         (1978),
!                        These Univ. Sci. Medic. Grenoble, Fig. 1 p. 124

             DO ikl=1,knonv             ! Parameterization of u*th
               rCd10n      =  1./  26.5 ! was developed from observations made
             END DO                     ! during assumed neutral conditions
                                         
               write(6,5000)  1./  rCd10n
 5000          format(/,' Blowing Snow Model  Initialization     ',
     .                /,' Vt / u*t =',f8.2,' (Neutral Assumption)',
     .                /,'           ', 8x ,' (Budd assumes  26.5)',/)
      END IF


! Snow Age (Influence on Snow Erosion Threshold)
! ==============================================

c #BS DO isn=1,nsno
c #BS DO ikl=1,knonv
c #BS   agsnSV(ikl,isn) = agsnSV(ikl,isn) + dt__SV/86400.
c #BS END DO
c #BS END DO
c #BS DO ikl=1,knonv
c #BS   isn    = max(1 ,        isnoSV(ikl))
c #BS   snofOK = max(0.,sign(1.,dsn_SV(ikl)-eps6))      !  Threshold=1.e-6
c #BS   agsnSV(ikl,isn) =   (1.-snofOK) *agsnSV(ikl,isn)! ~0.1 mm w.e./day
c #BS END DO
      IF (.NOT.BloMod)                                     GO TO 1000
c #AG STOP '?!&~@|@[#@#] --- INCONSISTANT SNOW AGE --- EMERGENCY STOP'
 1000 CONTINUE


! EROSION
! =======

      DO isn = 1, nsno
      DO ikl = 1,knonv

! Below the high Snow Density Threshold  (ro__SV < BSnoRo)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DendOK   =  max(zer0,sign(un_1,eps6-G1snSV(ikl,isn)  ))  !
        SaltOK   =  min(1   , max(istdSV(2)-istoSV(ikl,isn),0))  !
        SnowOK   =  min(1   , max(isnoSV(ikl)      +1 -isn ,0))  ! Snow Switch

        G1snSV(ikl,isn) =      SnowOK *    G1snSV(ikl,isn)
     .                  + (1.- SnowOK)*min(G1snSV(ikl,isn),G1_dSV)
        G2snSV(ikl,isn) =      SnowOK *    G2snSV(ikl,isn)
     .                  + (1.- SnowOK)*min(G2snSV(ikl,isn),G1_dSV)

        SaltOK   =  SaltOK   * SnowOK
        SaltM1   = -0.750e-2 * G1snSV(ikl,isn) 
     .             -0.500e-2 * G2snSV(ikl,isn)+ 0.500e00
!       SaltM1   :  Guyomarc'h & Merindol, 1997, Ann. Glac.
!         CAUTION:  Guyomarc'h & Merindol Dendricity Sign is +
!         ^^^^^^^^                    MAR Dendricity Sign is -
        SaltM2   = -0.833d-2 * G1snSV(ikl,isn) 
     .             -0.583d-2 * G2snSV(ikl,isn)+ 0.833d00
        SaltMo   = (DendOK   * SaltM1 + (1.-DendOK) *     SaltM2       )

! Increased Mobility of Deposed (blown) Snow (Mann et al., 2000, JGR 105,
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Fig.2 p.24496 & text below)
        MIN_Mo   =  0.
c #MA   MIN_Mo   =  0.6 * exp(-agsnSV(ikl,isn)                  /AgBlow)
        SaltMo   =                                    max(SaltMo,MIN_Mo)

        SaltMo   =  SaltOK   * SaltMo + (1.-SaltOK) * min(SaltMo,SaltMx)
!       SaltMo   =  SaltOK   * SaltMo - (1.-SaltOK) *     0.9500 ! Tuning
        SaltMo   =         max(SaltMo ,  eps6-un_1)

        SaltSU   =     (1.00d0+SaltMo)     *FacSBS

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0   Salt_U   =        -log(SaltSU)     *FacUBS
! #b0   IF (ii__SV(ikl).EQ.iwr_SV.AND.jj__SV(ikl).EQ.jwr_SV        .AND.
! #b0.      nn__SV(ikl).EQ.nwr_SV.AND.isn        .EQ.isnoSV(ikl))
! #b0.    write(6,6010) isnoSV(ikl),G1snSV(ikl,isn)/G1_dSV
! #b0.                             ,G2snSV(ikl,isn)/G1_dSV
! #b0.                             ,ro__SV(ikl,isn),agsnSV(ikl,isn)
! #b0.                             ,SaltM1, SaltM2, SaltMo, Salt_U
! #b0.                             ,us__SV(ikl)   / rCd10n
 6010     format(/,'SISVAT_BSn',6x
     .           ,6x,i3,2x,'G1         =',f6.3,'   G2         =',f7.3
     .           ,      '   ro [kg/m3] =',f9.3,'   Age* [Day] =',f9.3
     .           ,   /,27x,'SaltM1     =',f6.3,'   SaltM2     =',f7.3
     .           ,      '   Mobility I.=',f9.3,'   Vt   [m/s] =',f9.3
     .           ,   /,27x,'            ', 6x ,'               ', 7x
     .           ,      '               ', 9x ,'   Vn10 [m/s] =',f9.3)

! Above the high Snow Density Threshold  (ro__SV > BSnoRo)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Por_BS      =  1.000       - ro__SV(ikl,isn)     /rhoIce
        ShearX      =                Por_BS/max(eps6,un_1-Por_BS) 
!       ShearX ==> Arg(sqrt(shear)) with shear = max shear stress in snow:
!       shear       =  3.420d00 * exp(-(Por_BS      +Por_BS)
!    .                                /max(eps6,un_1-Por_BS))
!                      see de Montmollin         (1978),
!                      These Univ. Sci. Medic. Grenoble, Fig. 1 p. 124

! Influence of Density on Shear Stress if ro__SV > BSnoRo
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ArgFac      =  max(zer0  ,SheaBS-ShearX)     !
!       Fac_Mo      =  exp(       ArgFac       )     ! ** NOT ** tuned
        Fac_Mo   =     exp(       ArgFac       )     ! = 1 if ro__SV < BSnoRo
                                                     ! < 1 if ro__SV > BSnoRo
! Snow Drift Index
! ~~~~~~~~~~~~~~~~
        SaltSU      =  max(eps6  ,    SaltSU)
        SaltSU      =  exp(Fac_Mo*log(SaltSU))
        ArguSi      =     -FacTBS              *us__SV(ikl)/rCd10n
        SaltSI(ikl,isn) = (SaltSU-exp(ArguSi)) *FacRBS
!       SaltSI          :  Generalization of the Snow Drift Index of                       
!                          Guyomarc'h & Merindol (1997, Ann.Glaciol.)

! Threshold Friction Velocity
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SnowOK   =  1 -min(1,iabs(isn-isnoSV(ikl)))
        Salt_U   =               -log(SaltSU)  *FacUBS
!       Salt_U   :  Guyomarc'h & Merindol, 1997, Ann. Glac.

        usthSV(ikl) =     SnowOK *   (Salt_U   *rCd10n)
     .              + (1.-SnowOK)*    usthSV(ikl)

c #BA   usthSV(ikl) =     SnowOK *   (Salt_U   /26.5)
c #BA.              + (1.-SnowOK)*    usthSV(ikl)
!       Us(U10)     :  Budd et al.            1966, Ant.Res.Ser.9
!                 (see Pomeroy & Gray 1995 NHRI Sci.Rep.7(30)p.62)

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0   IF (ii__SV(ikl).EQ.iwr_SV.AND.jj__SV(ikl).EQ.jwr_SV        .AND.
! #b0.      nn__SV(ikl).EQ.nwr_SV.AND.isn        .EQ.isnoSV(ikl))
! #b0.    write(6,6011)     Fac_Mo,Por_BS,SaltSI(ikl,isn),usthSV(ikl)
 6011     format(      27x,'Fac_Mo     =',f6.3,'   Por_BS     =',f7.3
     .           ,      '   Drift    I.=',f9.3,'   ut*_0[m/s] =',f9.3)
      END DO
      END DO


! Deepest Mobile Snow Layer
! -------------------------

      DO ikl = 1,knonv
        Mobile(ikl) = nsno+1
      END DO
      DO isn =   nsno ,1,-1
      DO ikl = 1,knonv
        isnMAX      =       max(   1,          isnoSV(ikl)             )
        isnMAX      =       min( isn,          isnMAX                  )
        Mobiln      = isn * max(zer0,sign(un_1,SaltSI(ikl     ,isnMAX)))
        Mobilm      =   1 - min(1   ,          Mobile(ikl) -1 -Mobiln)
!       Mobilm      =   1   ONLY IF   Mobiln = Mobile(ikl) -1 (0 otherwise)

        Mobile(ikl) =                 Mobilm * Mobiln 
     .              +              (1-Mobilm)* Mobile(ikl)
      END DO
      END DO


! Weighting the Amount of Snow to erode
! -------------------------------------

      DO ikl = 1,knonv
        zdrift(ikl) = 0.0 
        xdrift(ikl) = 0.0
        dbsaux(ikl) = dbs_SV(ikl)
      END DO

      DO isn = 1, nsno
      DO ikl = 1,knonv
        zdrift(ikl)     =  zdrift(ikl)
     .            + 0.50 * dzsnSV(ikl,isn) * (3.25  -SaltSI(ikl,isn))
        sdrift(ikl,isn) =  SaltSI(ikl,isn)
     .          *exp(  max(Ea_Min, -zdrift(ikl)     *hdrift     ))
     .          *min(1,max(0     ,  isn +1          -Mobile(ikl)))
     .          *min(1,max(0     ,  isnoSV(ikl)     -isn +1     ))
!                Last 2 Lines force sdrift = 0 outside mobile Snow Layers
     .          *      max(zer0, sign(un_1,         -dbs_SV(ikl)))
!                Erosion is allowed only if available Blowing Snow
        xdrift(ikl)     =           sdrift(ikl,isn) +xdrift(ikl)
        zdrift(ikl)     =  zdrift(ikl)
     .            + 0.50 * dzsnSV(ikl,isn) * (3.25  -SaltSI(ikl,isn))
      END DO
      END DO

! Normalization
! ~~~~~~~~~~~~~
      DO isn = 1, nsno
      DO ikl = 1,knonv
        sdrift(ikl,isn) =  sdrift(ikl,isn) /max(eps6,xdrift(ikl))
      END DO
      END DO


! Weighting the Amount of Snow to depose
! --------------------------------------

      DO ikl = 1,knonv
        zdrift(ikl) = 0.0 
        tdepos(ikl) = 0.0
      END DO

      DO isn = 1, nsno
      DO ikl = 1,knonv
        zdepos(ikl,isn) =      exp(-zdrift(ikl)   )
     .          *min(1,max(0     ,  isn +1          -Mobile(ikl)))
     .          *min(1,max(0     ,  isnoSV(ikl    ) -isn +1     ))
!                Last 2 Lines force zdepos = 0 outside mobile Snow Layers
        tdepos(ikl) = tdepos(ikl) + zdepos(ikl,isn)
        zdrift(ikl) = zdrift(ikl) + dzsnSV(ikl,isn) *ro__SV(ikl,isn)
     .                                              /rhoWat
      END DO
      END DO

! Normalization
! ~~~~~~~~~~~~~
      DO isn = 1, nsno
      DO ikl = 1,knonv
        zdepos(ikl,isn) = zdepos(ikl,isn) / max(eps6,tdepos(ikl))
      END DO
      END DO

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0 DO ikl = 1,knonv
! #b0   IF (ii__SV(ikl).EQ.iwr_SV.AND.jj__SV(ikl).EQ.jwr_SV        .AND.
! #b0.      nn__SV(ikl).EQ.nwr_SV                          )        THEN
! #b0     Sno0WE =   0.
! #b0   DO isn=1,nsno
! #b0     Sno0WE =   Sno0WE 
! #b0.           +   dzsnSV(ikl,isn) *ro__SV(ikl,isn)
! #b0   END DO
! #b0     write(6,6005)   Sno0WE                    ,dbs_SV(ikl)
 6005     format(
     .      18x,'MB0',6x,'Sno1WE [mm]=',f9.3,19x,'0  dbs_SV [mm]=',f9.6)
! #b0     SnodWE =   dbs_SV(ikl)
! #b0   END IF
! #b0 END DO


! Weighted  Erosion (Erosion amount is distributed       ! dbs_SV decreases
! -----------------  over the upper Snow Pack)           ! dzsnSV decreases

      DO isn = 1, nsno
      DO ikl = 1,knonv
        SnowOK      = min(1,max(isnoSV(ikl)+1-isn ,0))      ! Snow Switch
        dzweqo      = dzsnSV(ikl,isn) *ro__SV(ikl,isn)      ! [kg/m2, mm w.e.]
        bsno_x      = dbsaux(ikl)     *sdrift(ikl,isn)
        dzweqn      = dzweqo          +bsno_x
        dzweqn  = max(dzweqn,          h_mmWE *SnowOK)
        dbs_SV(ikl) = dbs_SV(ikl)    +(dzweqo -dzweqn)
        dzsnSV(ikl,isn) =              dzweqn
     .                       /max(eps6,ro__SV(ikl,isn))
      END DO
      END DO

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0 DO ikl = 1,knonv
! #b0   IF (ii__SV(ikl).EQ. 1 .AND.  jj__SV(ikl) .EQ. 1)            THEN
! #b0     SnodWE =   SnodWE         -dbs_SV(ikl)
! #b0     Sno1WE =   0.
! #b0   DO isn=1,nsno
! #b0     Sno1WE =   Sno1WE 
! #b0.           +   dzsnSV(ikl,isn)*ro__SV(ikl,isn)
! #b0   END DO
! #b0     write(6,6006)Sno1WE      , dbs_SV(ikl)
 6006     format(
     .      18x,'MB1',6x,'Sno1WE [mm]=',f9.3,19x,'1  dbs_SV [mm]=',f9.6)
! #b0     write(6,6007)Sno1WE    ,SnodWE   ,Sno0WE,
! #b0.                (Sno1WE    -SnodWE   -Sno0WE)     
 6007     format(
     .      18x,'MB ',5x,'(After  [mm]=',f6.0, ')-(Erosion[mm]=', f7.3,
     .                                         ')-(Before [mm]=', f9.3,
     .                                         ')= Budget [mm]=', f9.6)
! #b0   END IF
! #b0 END DO


! ACCUMULATION of BLOWN SNOW                             ! dsn_SV decreases
! --------------------------                             ! dzsnSV increases

        DO ikl = 1,knonv
          tdepos(ikl) = dsn_SV(ikl) * dsnbSV(ikl) * dt__SV
          WEagre(ikl) = 0.

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0     IF (ii__SV(ikl).EQ.iwr_SV.AND.jj__SV(ikl).EQ.jwr_SV      .AND.
! #b0.        nn__SV(ikl).EQ.nwr_SV.AND.0          .LT.isnoSV(ikl))
! #b0.        write(6,6003) tdepos(ikl)  ,Mobile(ikl)
 6003         format(/,41x,'tdepos [-] =',f6.3,40x,'Mobil',i3
     .              ,/,27x,'Salt.Index    sdrift'
     .              ,      '    zdepos  ro__snow  ro_bsnow  roN_snow'
     .              ,                '  dz__snow  dz_bsnow  dzN_snow'
     .              ,                '  d___snow'
     .              ,/,27x,'             [kg/m3]   [kg/m3]   [kg/m3]'
     .              ,                '       [m]       [m]       [m]'
     .              ,                '   [kg/m2]')
        END DO

      DO isn =     nsno,1,-1
        DO ikl = 1,knonv
          WEagre(ikl) = WEagre(ikl) + ro__SV(ikl,isn)*dzsnSV(ikl,isn)
          isagr1(ikl) = istoSV(ikl,isn)
          isagr2(ikl) = 0.

! Density of deposited blown Snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ro_new =                    BSnoRo

! Density of deposited blown Snow (de Montmollin, 1978)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #EM     PorSno =      1.0d00     -  ro__SV(ikl,isn)
c #EM.                             /  rhoIce
c #EM     Salt_f =      usthSV(ikl)/  max(eps6,   us__SV(ikl))
c #EM     Salt_f =  min(Salt_f     ,  un_1)
c #EM     PorRef =      PorSno     /  max(eps6,1.-PorSno)
c #EM.             +log(Salt_f)
c #EM     Por_BS =      PorRef     / (1.0d00 + PorRef)
c #EM     ro_new =      rhoIce     * (1.0d00 - Por_BS)
c #EM     ro_new =  max(ro_new     ,  BSnoRo)

          roagr1(ikl) = ro__SV(ikl,isn)
          roagr2(ikl) = ro_new
          hsno_x      = tdepos(ikl)*  zdepos(ikl,isn)

          dzagr1(ikl) = dzsnSV(ikl,isn)
          dzagr2(ikl) = hsno_x     /  ro_new
!         Conversion    [kg/m2, i.e., mm w.e.] -----> [mSnow]

          dsn_SV(ikl) = dsn_SV(ikl)-  hsno_x / dt__SV

! Other Snow Properties
! ~~~~~~~~~~~~~~~~~~~~~
          T_agr1(ikl) =    TsisSV(ikl,isn)
          T_agr2(ikl) =min(Tf_Sno,TaT_SV(ikl))
          etagr1(ikl) =    eta_SV(ikl,isn)
          etagr2(ikl) =    0.0
          G1agr1(ikl) =    G1snSV(ikl,isn)
          G1agr2(ikl) =    G1_dSV
          G2agr1(ikl) =    G2snSV(ikl,isn)
          G2agr2(ikl) =    ADSdSV
c #BY     G2agr2(ikl) =    0.87d0
!         Budd et al. 1966, 2~m Average /Table 5 p. 97

          agagr1(ikl) =    agsnSV(ikl,isn)
          agagr2(ikl) =    0.
          Agrege(ikl) =    1.
        END DO

! Agregation
! ~~~~~~~~~~
!            **********
        call SISVAT_zAg
     .                 (isagr1,isagr2,WEagre
     .                 ,dzagr1,dzagr2,T_agr1,T_agr2
     .                 ,roagr1,roagr2,etagr1,etagr2
     .                 ,G1agr1,G1agr2,G2agr1,G2agr2
     .                 ,agagr1,agagr2,Agrege
     .                 )
!            **********

        DO ikl = 1,knonv

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0     IF (ii__SV(ikl).EQ.iwr_SV.AND.jj__SV(ikl).EQ.jwr_SV      .AND.
! #b0.        nn__SV(ikl).EQ.nwr_SV.AND.isn        .LE.isnoSV(ikl))
! #b0.        write(6,6004)   isn          ,SaltSI(ikl,isn)
! #b0.                     ,sdrift(ikl,isn),zdepos(ikl,isn)
! #b0.                     ,ro__SV(ikl,isn),roagr2(ikl),roagr1(ikl)
! #b0.                     ,dzsnSV(ikl,isn),dzagr2(ikl),dzagr1(ikl)
! #b0.                     ,dsn_SV(ikl)
 6004         format((27x,i3,f7.2,2f10.6,3f10.3,4f10.6))

          istoSV(ikl,isn) = isagr1(ikl)
          dzsnSV(ikl,isn) = dzagr1(ikl)
          TsisSV(ikl,isn) = T_agr1(ikl)
          ro__SV(ikl,isn) = roagr1(ikl)
          eta_SV(ikl,isn) = etagr1(ikl)
          G1snSV(ikl,isn) = G1agr1(ikl)
          G2snSV(ikl,isn) = G2agr1(ikl)
          agsnSV(ikl,isn) = agagr1(ikl)

        END DO

      END DO

! OUTPUT in SISVAT for ikl = 1 (preferably for Stand Alone Version)
! OUTPUT           for SnowFall and Snow Buffer
! #s2   IF          (isnoSV(1) .GT. 0)
! #s2.  write(6,6008)isnoSV(1),    dsn_SV(1) *dt__SV + BufsSV(1),
! #s2.              (dzsnSV(1,isn)*ro__SV(1,isn),isn=1,isnoSV(1))
 6008   format(i3,'  dsn+Buf=',f6.2,6x,'A dz *ro =',10f6.2,
     .                                       (/,35x,10f6.2))

        DO ikl = 1,knonv
          hdrift      =  tdepos(ikl)/dt__SV
          esnbSV(ikl) = (dsnbSV(ikl)-1.00)*hdrift/max(dsn_SV(ikl),eps6)
     .                  +dsnbSV(ikl)
          dsnbSV(ikl) =          min(un_1,   max(zer0,esnbSV(ikl) )   )
!         dsnbSV is now the Blown Snow fraction of precipitating snow
!                will be used for characterizing the Buffer Layer
!               (see update of  Bros_N, G1same, G2same, zroOLD, zroNEW)
        END DO

      return
      END


      subroutine SISVAT_BDu

!--------------------------------------------------------------------------+
!   MAR          SISVAT_BDu                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_BDu treats Dust Erosion                            |
!--------------------------------------------------------------------------+
!                                                                          |
!     OUTPUT:  usthSV   : Blowing Snow Erosion   Threshold           [m/s] |
!     ^^^^^^                                                               |
!                                                                          |
!     REFER. : Fecan, F., B. Marticorena and G. Bergametti, 1999   (Fal99) |
!     ^^^^^^^^ Ann. Geophysicae 17, 149--157                               |
!              u* threshold: adapted from Fig. 4 p. 153                    |
!              Clay Content:         from Tab. 2 p. 155                    |
!                                                                          |
!--------------------------------------------------------------------------+


       


! General Variables
! =================

      USE PHY_SV

      USE VAR_SV
      USE VARdSV
      USE VARxSV

      IMPLICIT NONE

! Local   Variables
! =================

      integer              ikl   ,  isot
      real                 ustdmn(0:nsot)
      real                 claypc(0:nsot)
      real                 f__ust(0:nvgt)
      real                 etaust(0:nsot)
      common /SISVAT_BDu_r/etaust
      logical              logust
      common /SISVAT_BDu_L/logust
      real                 eta_Du,usthDu


! DATA
! ====

      data         (ustdmn(isot),
     .                        claypc(isot),
     .                                                   isot=0,nsot)
     .             /10.000,  0.0000, !  0 WATER           !
     .               0.300,  0.0000, !  1 SAND            !
     .               0.300,  0.0920, !  2 LOAMY      SAND ! Fal99, Table 2
     .               0.300,  0.1420, !  3 SANDY      LOAM ! Fal99, Table 2
     .               0.300,  0.1630, !  4 SILT       LOAM ! Guess (Interpol.)
     .               0.300,  0.1840, !  5            LOAM ! Fal99, Table 2
     .               0.300,  0.2280, !  6 SANDY CLAY LOAM ! Guess (Interpol.)
     .               0.300,  0.2720, !  7 SILTY CLAY LOAM ! Guess (Interpol.)
     .               0.300,  0.3160, !  8       CLAY LOAM ! Fal99, Table 2
     .               0.300,  0.3750, !  9 SANDY CLAY      ! Guess (Interpol.)
     .               0.300,  0.4340, ! 10 SILTY CLAY      ! Guess (Interpol.)
     .               0.300,  0.4920, ! 11       CLAY      ! Fal99, Table 2
     .              10.000,  0.0000/ ! 12       ICE       !
      data          (f__ust(isot),                       isot=0,nvgt)
     .              /1.00,           !  0 NO     VEGETATION
     .               1.20,           !  1 CROPS      LOW
     .               5.00,           !  2 CROPS      MEDIUM
     .              10.00,           !  3 CROPS      HIGH
     .               1.20,           !  4 GRASS      LOW
     .               5.00,           !  5 GRASS      MEDIUM
     .              10.00,           !  6 GRASS      HIGH
     .               5.00,           !  7 BROADLEAF  LOW
     .              10.00,           !  8 BROADLEAF  MEDIUM
     .              12.00,           !  9 BROADLEAF  HIGH
     .              10.00,           ! 10 NEEDLELEAF LOW
     .              12.00,           ! 11 NEEDLELEAF MEDIUM
     .              50.00          / ! 12 NEEDLELEAF HIGH

      IF (.NOT.logust)                                              THEN
        DO isot=1,nsot
               etaust(isot) = 0.0014 * claypc(isot) * claypc(isot) ! Fal99
     .                      + 0.17   * claypc(isot)                ! Eqn.(14)
        END DO                                                     !  p. 154
               logust = .true.
      END IF


! Soil Erodibility
! ----------------

      DO ikl = 1,knonv
        eta_Du      =  max(     eta_SV(ikl,0),etaust(isotSV(ikl))) ! Fal99
        eta_Du      =  max(eps6,eta_SV(ikl,0)-eta_Du             ) ! Eqn.(15)
        usthDu      = sqrt(un_1+1.21*exp(0.68*   log(eta_Du)    )) ! p.  155
     .              *                         ustdmn(isotSV(ikl))  !
     .              *                         f__ust(ivgtSV(ikl))  !
        usthSV(ikl) = 
     .                 usthSV(ikl)*(1-max(0,1-isnoSV(ikl))) +
     .                 usthDu     *   max(0,1-isnoSV(ikl))
      END DO

      return
      END


      subroutine SISVAT_SIc
! #m2.                     (SIvAcr)

!--------------------------------------------------------------------------+
!   MAR          SISVAT_SIc                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_SIc treats Sea-Ice and Ocean Latent Heat Exchanges |
!--------------------------------------------------------------------------+
!                                                                          |
!     INPUT:   TaT_SV   : SBL Top    Temperature                       [K] |
!     ^^^^^    isnoSV   : total Nb of Ice/Snow Layers                  [-] |
!              LSmask   : Land-Sea   Mask                              [-] |
!              dsn_SV   : Snow  Intensity                      [mm w.e./s] |
!                                                                          |
!     INPUT /  TsisSV   : Snow/Ice/Soil-Water Temperature              [K] |
!     OUTPUT:  eta_SV   : Soil/Snow Water   Content                [m3/m3] |
!     ^^^^^^   dzsnSV   : Snow Layer        Thickness                  [m] |
!                                                                          |
!     OUTPUT:  HFraSV   : Frazil            Thickness                  [m] |
!     ^^^^^^                                                               |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #IA: Sea-Ice Bottom   accretion  and  ocean cooling                  |
!                                                                          |
!--------------------------------------------------------------------------+


       


! General Variables
! =================

      USE PHY_SV

      USE VAR_SV
      USE VARdSV


! INPUT/OUTPUT
! ------------

      USE VARxSV


      IMPLICIT NONE

! OUTPUT/Verification: SeaIce Conservation
! #m2 real      SIvAcr(klonv)                 ! Sea-Ice      Vertical Acretion


! Local Variables
! ===============

      integer             ikl   ,isn
      real                OCN_OK
      real                SIceOK
      real                SIcFrz
      real                Twat_n
      real                Crodzw,Lro__I
      common/SISVAT_SIc_R/Crodzw,Lro__I
      logical             SIcINI
      common/SISVAT_SIc_L/SIcINI

      real                SalIce              ! Sea-Ice   Salinity
      real                SalWat              ! Sea-Water Salinity


! DATA
! ====

      data                SalIce /10./        ! Sea-Ice   Salinity
      data                SalWat /35./        ! Sea-Water Salinity
!                         Typical Salinities in Terra Nova Bay
!                         (Bromwich and Kurtz,   1984, JGR, p.3568;
!                          Cavalieri and Martin, 1985,      p. 248)


! Initialisation
! ==============

      IF (.NOT.SIcINI)                                              THEN
               SIcINI =  .true.
               Crodzw =  hC_Wat*rhoWat *          dz_dSV(0)   ! [J/m2/K]
               Lro__I =  LhfH2O*rhoIce *(1.-1.e-3*SalIce      ! [J/m3]
     .                 -(SalIce/SalWat)*(1.-1.e-3*SalWat) )   !

! OUTPUT/Verification: Energy/Water Budget
! #e1          Lro__I =  LhfH2O*rhoIce

      END IF


! Snow Fall cools Sea Water
! =========================

      DO ikl=1,knonv
        OCN_OK        =  (1   -LSmask(ikl)   )                ! Free Ocean
     .             *max(0,1   -isnoSV(ikl)   )                !
c #IA   TsisSV(ikl,0) =        TsisSV(ikl,0)                  ! [K]
c #IA. -OCN_OK*(Cn_dSV*(Tf_Sno-TaT_SV(ikl)   )                ! [J/kg]
c #IA.         +LhfH2O*(1.    -eta_SV(ikl,0)))                ! [J/kg]
c #IA.        * dsn_SV(ikl)   *dt__SV          / Crodzw       ! [kg/m2]


! Sea-Ice Formation
! =================

c #IA   Twat_n      =      max(TsisSV(ikl,0  )  ,Tf_Sea)      ! [K]
c #IA   SIcFrz      =  (Twat_n-TsisSV(ikl,0  ) )*Crodzw/Lro__I! [m]
c #IA.                              * 0.75
!  ***  Hibler (1984), Ocean Heat Flux: 25% of cooling (ANTARCTIC Ocean)
!      (Hansen and Takahashi Eds) 
!       Geophys. Monogr. 29, M. Ewing Vol. 5, AGU, p. 241


! Frazil  Formation
! -----------------

c #IA   HFraSV(ikl) =          SIcFrz           *OCN_OK


! Growth of the Sea-Ice First Ice Floe
! ------------------------------------

c #IA   SIceOK        =  (1   -LSmask(ikl)     )              ! Ice Cover.Ocean
c #IA.             *min(  1   ,isnoSV(ikl)     )              !
c #IA   dzsnSV(ikl,1) =        dzsnSV(ikl,1)                  ! Vertical Acret.
c #IA.                +        SIcFrz           *SIceOK       !


! OUTPUT/Verification: SeaIce Conservation: Diagnostic of Surface Mass Balance
! #m2   SIvAcr(ikl) = rhoIce*SIcFrz     *(OCN_OK+SIceOK)
! #m2.              - dt__SV*dsn_SV(ikl)* OCN_OK


! Water Fluxes Update
! -------------------

        RnofSV(ikl) =          RnofSV(ikl)
     .              +          dsn_SV(ikl) *     OCN_OK
        dsn_SV(ikl) =          dsn_SV(ikl) * (1.-OCN_OK)

      END DO

      return
      end
      subroutine SISVAT_zSn

!--------------------------------------------------------------------------+
!   MAR          SISVAT_zSn                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_zSn manages the Snow Pack vertical Discretization  |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT /  NLaysv   = New             Snow Layer  Switch               |
!     OUTPUT:  isnoSV   = total Nb of Ice/Snow Layers                      |
!     ^^^^^^   ispiSV   = 0,...,nsno: Uppermost Superimposed Ice Layer     |
!              iiceSV   = total Nb of Ice      Layers                      |
!              istoSV   = 0,...,5 :   Snow     History (see istdSV data)   |
!                                                                          |
!     INPUT /  TsisSV   : Soil/Ice Temperatures (layers -nsol,-nsol+1,..,0)|
!     OUTPUT:           & Snow     Temperatures (layers  1,2,...,nsno) [K] |
!     ^^^^^^   ro__SV   : Soil/Snow Volumic Mass                   [kg/m3] |
!              eta_SV   : Soil/Snow Water   Content                [m3/m3] |
!              dzsnSV   : Snow Layer        Thickness                  [m] |
!              G1snSV   : Dendricity (<0) or Sphericity (>0) of Snow Layer |
!              G2snSV   : Sphericity (>0) or Size            of Snow Layer |
!              agsnSV   : Snow       Age                             [day] |
!                                                                          |
!     METHOD:  1) Agregate the thinest Snow Layer                          |
!     ^^^^^^      if a new Snow Layer has been precipitated   (NLaysv = 1) |
!              2) Divide   a too thick Snow Layer except                   |
!                 if the maximum Number of Layer is reached                |
!                 in this case forces                          NLay_s = 1  |
!              3) Agregate the thinest Snow Layer                          |
!                 in order to divide a too thick Snow Layer                |
!                 at next Time Step when                       NLay_s = 1  |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_zSn.vz        | #vz: OUTPUT/Verification: Snow Layers Agrega. |
!                          |      unit 41, SubRoutine  SISVAT_zSn **ONLY** |
!   # SISVAT_GSn.vp        | #vp: OUTPUT/Verification: Snow   Properties   |
!                          |      unit 47, SubRoutines SISVAT_zSn, _GSn    |
!   # stdout               | #s1: OUTPUT of Snow Layers Agregation         |
!                          |      unit  6, SubRoutine  SISVAT_zSn, _zAg    |
!--------------------------------------------------------------------------+
      USE PHY_SV 
                                 
      USE VAR_SV
      USE VARdSV
                       
      USE VARxSV 
      USE VARySV
                                 
      USE VAR0SV                                                 
!      USE VARphy 

      IMPLICIT NONE


! Global Variables
! ================

!      include  'PHY_SV.h'
!
!      include  'MAR_SV.inc'
!      include  'MARdSV.inc'
!      include  'MAR0SV.inc'
!
!      include  'MARxSV.inc'
!

! Internal Variables
! ==================

      integer   ikl   ,isn   ,i               !

!      include  'MARySV.inc'
      integer   NLay_s(klonv)                 ! Split Snow Layer         Switch
      integer   isagr1(klonv)                 ! 1st     Layer History
      integer   isagr2(klonv)                 ! 2nd     Layer History
      integer   LstLay                        ! 0 ====> isnoSV = 1
      integer   isno_n                        ! Snow Normal.Profile
      integer   iice_n                        ! Ice  Normal.Profile
      integer   iiceOK                        ! Ice         Switch
      integer   icemix                        ! 0 ====> Agregated Snow+Ice=Snow
                                              ! 1                          Ice
      integer   isn1  (klonv)                 ! 1st layer to stagger
      real      staggr                        !              stagger  Switch

      real      WEagre(klonv)                 ! Snow Water Equivalent Thickness
      real      dzthin(klonv)                 ! Thickness of the thinest layer
      real      OKthin                        ! Swich ON  a  new thinest layer
      real      dz_dif                        ! difference from ideal discret.
      real      thickL                        ! Thick Layer          Indicator
      real      OK_ICE                        ! Swich ON   uppermost Ice Layer

      real      Agrege(klonv)                 ! 1. when Agregation constrained
      real      dzepsi                        ! Min Single Snw Layer Thickness
      real      dzxmin                        ! Min Acceptable Layer Thickness
      real      dz_min                        ! Min            Layer Thickness
      real      dz_max                        ! Max            Layer Thickness
      real      dzagr1(klonv)                 ! 1st     Layer Thickness
      real      dzagr2(klonv)                 ! 2nd     Layer Thickness
      real      T_agr1(klonv)                 ! 1st     Layer Temperature
      real      T_agr2(klonv)                 ! 2nd     Layer Temperature
      real      roagr1(klonv)                 ! 1st     Layer Density
      real      roagr2(klonv)                 ! 2nd     Layer Density
      real      etagr1(klonv)                 ! 1st     Layer Water Content
      real      etagr2(klonv)                 ! 2nd     Layer Water Content
      real      G1agr1(klonv)                 ! 1st     Layer Dendricity/Spher.
      real      G1agr2(klonv)                 ! 2nd     Layer Dendricity/Spher.
      real      G2agr1(klonv)                 ! 1st     Layer Sphericity/Size
      real      G2agr2(klonv)                 ! 2nd     Layer Sphericity/Size
      real      agagr1(klonv)                 ! 1st     Layer Age
      real      agagr2(klonv)                 ! 2nd     Layer Age

! OUTPUT/Verification: Snow Layers Agregation
! #vz logical         as_opn                  ! IO   Switch
! #vz common/SI_zSn_L/as_opn                  !
! #vz real      dz_ref(            nsno)      ! Snow Reference Discretization
! #vz real      dzwdif(            nsno)      !

! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp logical         VP_opn                  ! IO   Switch
! #vp common/SI_GSn_L/VP_opn                  !


! DATA
! ====

!      data      icemix /   0    /             ! 0 ====> Agregated Snow+Ice=Snow
!      data      dzepsi /   0.0015/            ! Min single Layer Thickness 
!      data      dzxmin /   0.0020/            ! Min accept.Layer Thickness 
!      data      dz_min /   0.0050/            ! Min Local  Layer Thickness 
!      data      dz_max /   0.0300/            ! Min Gener. Layer Thickness 
!hjp250711
      data      dzepsi /   0.0045/            ! Min single Layer Thickness 
      data      dzxmin /   0.0060/            ! Min accept.Layer Thickness 
      data      dz_min /   0.0150/            ! Min Local  Layer Thickness 
      data      dz_max /   0.0900/            ! Min Gener. Layer Thickness 
!     CAUTION:  dz_max > dz_min*2 is required ! Otherwise re-agregation is
                                              ! activated  after splitting


! OUTPUT/Verification: Snow Layers Agregation
! #vz IF (.NOT.as_opn)                                              THEN
! #vz          as_opn=.true.
! #vz     open(unit=41,status='unknown',file='SISVAT_zSn.vz')
! #vz     rewind    41
! #vz END IF

! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp IF (.NOT.VP_opn)                                              THEN
! #vp          VP_opn=.true.
! #vp     open(unit=47,status='unknown',file='SISVAT_GSn.vp')
! #vp     rewind    47
! #vp END IF


! Constrains Agregation         of too thin  Layers
! =================================================

! Search the thinest  non-zero Layer
! ----------------------------------

        DO ikl=1,klonv
          dzthin(ikl) = 0.                              ! Arbitrary unrealistic
        END DO                                          !       Layer Thickness 
      DO   isn=1,nsno
        DO ikl=1,klonv
          isno_n    =             isnoSV(ikl)-isn+1     ! Snow Normal.Profile
          iice_n    =             iiceSV(ikl)-isn       ! Ice  Normal.Profile
          iiceOK    = min(1,max(0,iice_n         +1))   ! Ice         Switch

! OUTPUT/Verification: Snow Layers Agregation
! #vz     dz_ref(isn) =                                 !
! #vz.          dz_min *((1-iiceOK)*isno_n*isno_n       ! Theoretical Profile
! #vz.                 +    iiceOK *    2**iice_n)      !
! #vz.               /max(1,isnoSV(ikl))                !

          dz_dif      = max(zer0,                       ! Actual      Profile
     .          dz_min                                  !
     .                 *((1-iiceOK)*isno_n*isno_n       ! Theoretical Profile
     .                 +    iiceOK *2.   **iice_n)      !
     .        - dzsnSV(ikl, isn)                    )   ! Actual      Profile

! OUTPUT/Verification: Snow Layers Agregation
! #vz     dzwdif(isn) =     dz_dif                      !

          OKthin      = max(zer0,                       !
     .                      sign(un_1,                  !
     .                           dz_dif-dzthin(ikl)))   ! 1.=> New thinest Lay.
     .                * max(0,                          ! 1 => .le. isnoSV
     .                  min(1,                          ! 1 => isn is in the
     .                      isnoSV(ikl)-isn +1 ))       !          Snow Pack
     .                * min(un_1,                                !
!
!                       1st additional Condition to accept OKthin
     .                  max(zer0,                                ! combination
     .                      sign(un_1,G1snSV(ikl,      isn  )    ! G1 with same
     .                               *G1snSV(ikl,max(1,isn-1)))) !  sign => OK
!
!                       2nd additional Condition to accept OKthin
     .                + max(zer0,                                ! G1>0
     .                      sign(un_1,G1snSV(ikl,      isn   ))) !  =>OK
!
!                       3rd additional Condition to accept OKthin
     .                + max(zer0,                                ! dz too small
     .                      sign(un_1,dzxmin                     !  =>OK
     .                               -dzsnSV(ikl,      isn   ))))!

          i_thin(ikl) =    (1. - OKthin)  * i_thin(ikl) ! Update   thinest Lay.
     .                         + OKthin   * isn         !                Index
          dzthin(ikl) =    (1. - OKthin)  * dzthin(ikl) !
     .                         + OKthin   * dz_dif      !
        END DO
      END DO

! OUTPUT/Verification: Snow Layers Agregation
! #vz     write(41,4150) daHost   ,n___SV(  lwriSV(1))
! #vz.                  ,i_thin(1),dzsnSV(1,i_thin(1))
 4150     format(/,'-',a18,i5,' ',70('-'),
     .           /,' Thinest ',i3,':',f9.3)

      DO   isn=1,nsno
        DO ikl=1,klonv
          OKthin =      max(zer0,                       !
     .                      sign(un_1,                  !
     .                           dzxmin                 !
     .                          -dzsnSV(ikl,isn)))      !
     .                * max(zer0,                       ! ON if dz > 0
     .                      sign(un_1,                  !
     .                           dzsnSV(ikl,isn)-eps6)) !
     .           *min(1,max(0,                          ! Multiple Snow    Lay.
     .                      min (1,                     ! Switch = 1
     .                           isnoSV(ikl)            !   if isno > iice + 1
     .                          -iiceSV(ikl)-1))        !
                                                        !
     .             +int(max(zer0,                       !
     .                      sign(un_1,                  !
     .                           dzepsi                 ! Minimum accepted for
     .                          -dzsnSV(ikl,isn))))     ! 1 Snow Layer over Ice
     .             *int(max(zer0,                       ! ON if dz > 0
     .                      sign(un_1,                  !
     .                           dzsnSV(ikl,isn)-eps6)))!
     .                 *(1 -min (abs(isnoSV(ikl)        ! Switch = 1
     .                              -iiceSV(ikl)-1),1)) !   if isno = iice + 1
                                                        !
     .                 +max(0,                          ! Ice
     .                      min (1,                     ! Switch
     .                           iiceSV(ikl)+1-isn)))   !
     .             *min(un_1,                                    !
     .                  max(zer0,                                ! combination
     .                      sign(un_1,G1snSV(ikl,      isn  )    ! G1>0 + G1<0
     .                               *G1snSV(ikl,max(1,isn-1)))) ! NO
     .                + max(zer0,                                !
     .                      sign(un_1,G1snSV(ikl,      isn   ))) !
     .                + max(zer0,                                !
     .                      sign(un_1,dzxmin                     !
     .                               -dzsnSV(ikl,      isn   ))))!
          i_thin(ikl) =    (1. - OKthin)  * i_thin(ikl) ! Update   thinest Lay.
     .                         + OKthin   * isn         !                Index
        END DO
      END DO

! OUTPUT/Verification: Snow Layers Agregation
! #vz     write(41,4151) i_thin(1),dzsnSV(1,i_thin(1))
! #vz.                  ,isnoSV(1),dzsnSV(1,isnoSV(1))
 4151     format(' Thinest ',i3,':',f9.3,'   Max   =',i3,f12.3)

! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp   write(47,470)(G1snSV(1,isn),isn=1,isnoSV(1))
 470    format('Before _zCr1: G1 = ',10f8.1,(/,19x,10f8.1))
! #vp   write(47,472)(G2snSV(1,isn),isn=1,isnoSV(1))
 472    format('              G2 = ',10f8.1,(/,19x,10f8.1))


! Index of the contiguous Layer to agregate
! -----------------------------------------

!          **********
      call SISVAT_zCr
!          **********


! Assign the 2 Layers to agregate
! -------------------------------

        DO ikl=1,klonv
          isn         =    i_thin(ikl)
          isagr1(ikl) =    istoSV(ikl,isn)
          isagr2(ikl) =    istoSV(ikl,isn+LIndsv(ikl))
          dzagr1(ikl) =    dzsnSV(ikl,isn)
          dzagr2(ikl) =    dzsnSV(ikl,isn+LIndsv(ikl))
          T_agr1(ikl) =    TsisSV(ikl,isn)
          T_agr2(ikl) =    TsisSV(ikl,isn+LIndsv(ikl))
          roagr1(ikl) =    ro__SV(ikl,isn)
          roagr2(ikl) =    ro__SV(ikl,isn+LIndsv(ikl))
          etagr1(ikl) =    eta_SV(ikl,isn)
          etagr2(ikl) =    eta_SV(ikl,isn+LIndsv(ikl))
          G1agr1(ikl) =    G1snSV(ikl,isn)
          G1agr2(ikl) =    G1snSV(ikl,isn+LIndsv(ikl))
          G2agr1(ikl) =    G2snSV(ikl,isn)
          G2agr2(ikl) =    G2snSV(ikl,isn+LIndsv(ikl))
          agagr1(ikl) =    agsnSV(ikl,isn)
          agagr2(ikl) =    agsnSV(ikl,isn+LIndsv(ikl))
          LstLay      = min(1,max(  0,isnoSV(ikl) -1))  ! 0  if single Layer
          isnoSV(ikl) =               isnoSV(ikl)       ! decrement   isnoSV
     .     -(1-LstLay)* max(zer0,                       ! if downmost  Layer
     .                      sign(un_1,eps_21            ! <  1.e-21 m
     .                               -dzsnSV(ikl,1)))   !
          isnoSV(ikl) = max(   0,     isnoSV(ikl)   )   !
          Agrege(ikl) = max(zer0,                       !
     .                      sign(un_1,dz_min            ! No Agregation 
     .                               -dzagr1(ikl)  ))   ! if too thick Layer
     .                               *LstLay            ! if  a single Layer
     .                * min( max(0   ,isnoSV(ikl)+1     ! if Agregation
     .                               -i_thin(ikl)       !    with    a Layer
     .                               -LIndsv(ikl)  ),1) !    above the Pack

          WEagre(ikl) = 0.
        END DO

        DO isn=1,nsno
        DO ikl=1,klonv
          WEagre(ikl) = WEagre(ikl) + ro__SV(ikl,isn)*dzsnSV(ikl,isn)
     .                                *min(1,max(0,i_thin(ikl)+1-isn))
        ENDDO
        ENDDO

! OUTPUT/Verification: Snow Layers Agregation
! #vz write(41,410)
 410  format(/,' Agregation of too THIN Layers')
! #vz write(41,411) (100.*dz_ref(  isn),isn=1,nsno)
! #vz write(41,412) (100.*dzwdif(  isn),isn=1,nsno)
! #vz write(41,413) (100.*dzsnSV(1,isn),isn=1,nsno)
! #vz write(41,414) (              isn ,isn=1,nsno)
 411  format(' dz_ref [cm]:',10f8.2   ,/,('             ',10f8.2) )
 412  format(' dz_dif [cm]:',10f8.2   ,/,('             ',10f8.2) )
 413  format(' dzsnSV [cm]:',10f8.2   ,/,('             ',10f8.2) )
 414  format('             ',10(i5,3x),/,('             ',10(i5,3x)))
! #vz write(41,4111)      isnoSV(1    )
! #vz write(41,4112)      i_thin(1    )
! #vz write(41,4113)      LIndsv(1    )
! #vz write(41,4114)      Agrege(1    )
! #vz write(41,4115) 1.e2*dzagr1(1    )
! #vz write(41,4116) 1.e2*dzagr2(1    )
 4111 format(' isnoSV     :',  i8  )
 4112 format(' i_thin     :',  i8  )
 4113 format(' LIndsv     :',  i8  )
 4114 format(' Agrege     :',  f8.2)
 4115 format(' dzagr1     :',  f8.2)
 4116 format(' dzagr2     :',  f8.2)

! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp   write(47,471)(G1snSV(1,isn),isn=1,isnoSV(1))
 471    format('Before _zAg1: G1 = ',10f8.1,(/,19x,10f8.1))
! #vp   write(47,472)(G2snSV(1,isn),isn=1,isnoSV(1))


! Agregates
! ---------

!            **********
        call SISVAT_zAg
     .                 (isagr1,isagr2,WEagre
     .                 ,dzagr1,dzagr2,T_agr1,T_agr2
     .                 ,roagr1,roagr2,etagr1,etagr2
     .                 ,G1agr1,G1agr2,G2agr1,G2agr2
     .                 ,agagr1,agagr2,Agrege
     .                 )
!            **********


! Rearranges the Layers
! ---------------------

! New (agregated) Snow layer
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,klonv
          isn     =             i_thin(ikl)
          isn     = min(isn,isn+LIndsv(ikl))
          isnoSV(ikl) =         isnoSV(ikl) -Agrege(ikl)
          iiceSV(ikl) =         iiceSV(ikl) 
     .            -max(0,sign(1,iiceSV(ikl) -isn +icemix))
     .                                      *Agrege(ikl)
     .            *max(0,sign(1,iiceSV(ikl) -1          ))
          istoSV(ikl,isn) = (1.-Agrege(ikl))*istoSV(ikl,isn)
     .                      +   Agrege(ikl) *isagr1(ikl)
          dzsnSV(ikl,isn) = (1.-Agrege(ikl))*dzsnSV(ikl,isn)
     .                      +   Agrege(ikl) *dzagr1(ikl)
          TsisSV(ikl,isn) = (1.-Agrege(ikl))*TsisSV(ikl,isn)
     .                      +   Agrege(ikl) *T_agr1(ikl)
          ro__SV(ikl,isn) = (1.-Agrege(ikl))*ro__SV(ikl,isn)
     .                      +   Agrege(ikl) *roagr1(ikl)
          eta_SV(ikl,isn) = (1.-Agrege(ikl))*eta_SV(ikl,isn)
     .                      +   Agrege(ikl) *etagr1(ikl)
          G1snSV(ikl,isn) = (1.-Agrege(ikl))*G1snSV(ikl,isn)
     .                      +   Agrege(ikl) *G1agr1(ikl)
          G2snSV(ikl,isn) = (1.-Agrege(ikl))*G2snSV(ikl,isn)
     .                      +   Agrege(ikl) *G2agr1(ikl)
          agsnSV(ikl,isn) = (1.-Agrege(ikl))*agsnSV(ikl,isn)
     .                      +   Agrege(ikl) *agagr1(ikl)
        END DO

! Above
! ^^^^^
        DO ikl=1,klonv
          isn1(ikl)=max(i_thin(ikl),i_thin(ikl)+LIndsv(ikl))
        END DO
        DO i=  1,nsno-1
        DO ikl=1,klonv
            staggr        =  min(1,max(0,i +1 -isn1(ikl)   ))
            istoSV(ikl,i) = (1.-staggr     )*istoSV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*istoSV(ikl,i  )
     .                      +   Agrege(ikl) *istoSV(ikl,i+1))
            dzsnSV(ikl,i) = (1.-staggr     )*dzsnSV(ikl,i  ) 
     .            + staggr*((1.-Agrege(ikl))*dzsnSV(ikl,i  )
     .                      +   Agrege(ikl) *dzsnSV(ikl,i+1))
            TsisSV(ikl,i) = (1.-staggr     )*TsisSV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*TsisSV(ikl,i  )
     .                      +   Agrege(ikl) *TsisSV(ikl,i+1))
            ro__SV(ikl,i) = (1.-staggr     )*ro__SV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*ro__SV(ikl,i  )
     .                      +   Agrege(ikl) *ro__SV(ikl,i+1))
            eta_SV(ikl,i) = (1.-staggr     )*eta_SV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*eta_SV(ikl,i  )
     .                      +   Agrege(ikl) *eta_SV(ikl,i+1))
            G1snSV(ikl,i) = (1.-staggr     )*G1snSV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*G1snSV(ikl,i  )
     .                      +   Agrege(ikl) *G1snSV(ikl,i+1))
            G2snSV(ikl,i) = (1.-staggr     )*G2snSV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*G2snSV(ikl,i  )
     .                      +   Agrege(ikl) *G2snSV(ikl,i+1))
            agsnSV(ikl,i) = (1.-staggr     )*agsnSV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*agsnSV(ikl,i  )
     .                      +   Agrege(ikl) *agsnSV(ikl,i+1))
        END DO
        END DO

        DO ikl=1,klonv
          isn             = min(isnoSV(ikl) +1,nsno)
          istoSV(ikl,isn) = (1.-Agrege(ikl))*istoSV(ikl,isn)
          dzsnSV(ikl,isn) = (1.-Agrege(ikl))*dzsnSV(ikl,isn)
          TsisSV(ikl,isn) = (1.-Agrege(ikl))*TsisSV(ikl,isn)
          ro__SV(ikl,isn) = (1.-Agrege(ikl))*ro__SV(ikl,isn)
          eta_SV(ikl,isn) = (1.-Agrege(ikl))*eta_SV(ikl,isn)
          G1snSV(ikl,isn) = (1.-Agrege(ikl))*G1snSV(ikl,isn)
          G2snSV(ikl,isn) = (1.-Agrege(ikl))*G2snSV(ikl,isn)
          agsnSV(ikl,isn) = (1.-Agrege(ikl))*agsnSV(ikl,isn)
        END DO

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #s1   IF (ii__SV(ikl).EQ.iwr_SV.AND.jj__SV(ikl).EQ.jwr_SV        .AND.
! #s1.      nn__SV(ikl).EQ.nwr_SV                          )        THEN
! #s1     write(6,5991) i_thin(ikl)
 5991     format(/,'First Agregation / Layer',i3,
     .           /,'  i',11x,'T',9x,'rho',10x,'dz',11x,'H')
! #s1     write(6,5995) (isn,TsisSV(ikl,isn),ro__SV(ikl,isn)
! #s1.                      ,dzsnSV(ikl,isn),istoSV(ikl,isn),
! #s1.                   isn=isnoSV(ikl),1,-1)
 5995     format(i3,3f12.3,i12)
! #s1   END IF


! Constrains Splitting          of too thick Layers
! =================================================


! Search the thickest non-zero Layer
! ----------------------------------

        DO ikl=1,klonv
          dzthin(ikl) =   0.                            ! Arbitrary unrealistic
        END DO                                          !       Layer Thickness 
      DO   isn=1,nsno
        DO ikl=1,klonv
          isno_n    =             isnoSV(ikl)-isn+1     ! Snow Normal.Profile
          iice_n    =             iiceSV(ikl)-isn       ! Ice  Normal.Profile
          iiceOK    = min(1,max(0,iice_n         +1))   ! Ice         Switch
          dz_dif    =(      dzsnSV(ikl,isn)             ! Actual      Profile
     .        - dz_max *((1-iiceOK)*isno_n*isno_n       ! Theoretical Profile
     .                 +    iiceOK *2.   **iice_n)  )   !
     .                 /max(dzsnSV(ikl,isn),eps6)       !
          OKthin      = max(zer0,                       !
     .                      sign(un_1,                  !
     .                           dz_dif-dzthin(ikl)))   ! 1.=>New thickest Lay.
     .                * max(0,                          ! 1 =>.le. isnoSV
     .                  min(1,                          !
     .                      isnoSV(ikl)-isn +1 ))       !
          i_thin(ikl) =    (1. - OKthin)  * i_thin(ikl) !  Update thickest Lay.
     .                         + OKthin   * isn         !                Index
          dzthin(ikl) =    (1. - OKthin)  * dzthin(ikl) !
     .                         + OKthin   * dz_dif      !
        END DO
      END DO

      DO   ikl=1,klonv
          ThickL      = max(zer0,                       ! 1. => a too   thick
     .                      sign(un_1,dzthin(ikl)       !         Layer exists
     .                               -eps6       ))     !
     .        * max(0,1-max(0   ,     isnoSV(ikl)       ! No spliting allowed
     .                               -nsno+3     ))     ! if isno > nsno - 3
          Agrege(ikl) =               ThickL            ! 1. => effective split
     .        * max(0,1-max(0   ,     NLaysv(ikl)       !
     .                               +isnoSV(ikl)       !
     .                               -nsno+1     ))     ! 
          NLay_s(ikl) =               ThickL            ! Agregation
     .        * max(0,1-max(0   ,     NLaysv(ikl)       ! to allow  Splitting 
     .                               +isnoSV(ikl)       !   at next Time Step
     .                               -nsno       ))     !
     .                               -Agrege(ikl)       ! 
          NLay_s(ikl) = max(0   ,     NLay_s(ikl))      ! Agregation effective
      END DO

! OUTPUT/Verification: Snow Layers Agregation
! #vz     write(41,4152) i_thin(1),dzthin(1),ThickL
 4152     format(/,' Thickest',i3,':',f9.3,'   Split =',f4.0)


! Rearranges the Layers
! ---------------------

      DO isn=nsno,2,-1
      DO ikl=1,klonv
        IF (Agrege(ikl).gt.0..AND.i_thin(ikl).lt.isnoSV(ikl))       THEN
          staggr          =  min(1,max(0,isn-i_thin(ikl)    -1))
     .                    *  min(1,max(0,    isnoSV(ikl)-isn+2))
          istoSV(ikl,isn) =        staggr  * istoSV(ikl ,isn-1)
     .                    + (1. -  staggr) * istoSV(ikl ,isn  )
          dzsnSV(ikl,isn) =        staggr  * dzsnSV(ikl ,isn-1)
     .                    + (1. -  staggr) * dzsnSV(ikl ,isn  )
          TsisSV(ikl,isn) =        staggr  * TsisSV(ikl ,isn-1)
     .                    + (1. -  staggr) * TsisSV(ikl ,isn  )
          ro__SV(ikl,isn) =        staggr  * ro__SV(ikl ,isn-1)
     .                    + (1. -  staggr) * ro__SV(ikl ,isn  )
          eta_SV(ikl,isn) =        staggr  * eta_SV(ikl ,isn-1)
     .                    + (1. -  staggr) * eta_SV(ikl ,isn  )
          G1snSV(ikl,isn) =        staggr  * G1snSV(ikl ,isn-1)
     .                    + (1. -  staggr) * G1snSV(ikl ,isn  )
          G2snSV(ikl,isn) =        staggr  * G2snSV(ikl ,isn-1)
     .                    + (1. -  staggr) * G2snSV(ikl ,isn  )
          agsnSV(ikl,isn) =        staggr  * agsnSV(ikl ,isn-1)
     .                    + (1. -  staggr) * agsnSV(ikl ,isn  )
        END IF
      END DO
      END DO

      DO  ikl=1,klonv
          isn             =     i_thin(ikl)
          dzsnSV(ikl,isn) = 0.5*Agrege(ikl) *dzsnSV(ikl,isn)
     .                    + (1.-Agrege(ikl))*dzsnSV(ikl,isn)

          isn             = min(i_thin(ikl) +1,nsno)
          istoSV(ikl,isn) =     Agrege(ikl) *istoSV(ikl,isn-1)
     .                    + (1.-Agrege(ikl))*istoSV(ikl,isn)
          dzsnSV(ikl,isn) =     Agrege(ikl) *dzsnSV(ikl,isn-1)
     .                    + (1.-Agrege(ikl))*dzsnSV(ikl,isn)
          TsisSV(ikl,isn) =     Agrege(ikl) *TsisSV(ikl,isn-1)
     .                    + (1.-Agrege(ikl))*TsisSV(ikl,isn)
          ro__SV(ikl,isn) =     Agrege(ikl) *ro__SV(ikl,isn-1)
     .                    + (1.-Agrege(ikl))*ro__SV(ikl,isn)
          eta_SV(ikl,isn) =     Agrege(ikl) *eta_SV(ikl,isn-1)
     .                    + (1.-Agrege(ikl))*eta_SV(ikl,isn)
          G1snSV(ikl,isn) =     Agrege(ikl) *G1snSV(ikl,isn-1)
     .                    + (1.-Agrege(ikl))*G1snSV(ikl,isn)
          G2snSV(ikl,isn) =     Agrege(ikl) *G2snSV(ikl,isn-1)
     .                    + (1.-Agrege(ikl))*G2snSV(ikl,isn)
          agsnSV(ikl,isn) =     Agrege(ikl) *agsnSV(ikl,isn-1)
     .                    + (1.-Agrege(ikl))*agsnSV(ikl,isn)
          isnoSV(ikl)     =     Agrege(ikl) +isnoSV(ikl)
          iiceSV(ikl)     =                  iiceSV(ikl)
     .                    +     Agrege(ikl) *max(0,sign(1,iiceSV(ikl) 
     .                                                   -isn +icemix))
     .                                      *max(0,sign(1,iiceSV(ikl) 
     .                                                   -1          ))
      END DO


! Constrains Agregation in case of too much  Layers
! =================================================

! Search the thinest   non-zero Layer
! -----------------------------------

! OUTPUT/Verification: Snow Thinest Layer  
! #sd     write( 6,*)   ' '
! #sd     write( 6,*)   'Agregation 2'
! #sd     write( 6,6000) NLaysv(1)
 6000     format(i3,6x,
     .          'dzsnSV      dz_min      dz_dif      ',
     .          'OKthin      dzthin   i_thin')

        DO ikl=1,klonv
          dzthin(ikl) =   0.                            ! Arbitrary unrealistic
        END DO                                          !       Layer Thickness 
      DO   isn=1,nsno
        DO ikl=1,klonv
          isno_n    =             isnoSV(ikl)-isn+1     ! Snow Normal.Profile
          iice_n    =             iiceSV(ikl)-isn       ! Ice  Normal.Profile
          iiceOK    = min(1,max(0,iice_n         +1))   ! Ice         Switch

! OUTPUT/Verification: Snow Layers Agregation
! #vz     dz_ref(isn) =                                 !
! #vz.          dz_min *((1-iiceOK)*isno_n*isno_n       ! Theoretical Profile
! #vz.                 +    iiceOK *    2**iice_n)      !
! #vz.               /max(1,isnoSV(ikl))                !

          dz_dif      =     dz_min                      ! Actual      Profile
     .                    - dzsnSV(ikl    ,isn)         !
     .        /max(eps6,((1-iiceOK)*isno_n*isno_n       ! Theoretical Profile
     .                 +    iiceOK *2.   **iice_n))     !

! OUTPUT/Verification: Snow Layers Agregation
! #vz     dzwdif(isn) =     dz_dif                      !

          OKthin      = max(zer0,                       !
     .                      sign(un_1,                  !
     .                           dz_dif  - dzthin(ikl)))! 1.=> New thinest Lay.
     .                * max(0,                          ! 1 => .le. isnoSV
     .                  min(1,                          !
     .                      isnoSV(ikl)-isn +1 ))       !
          i_thin(ikl) =    (1. - OKthin) * i_thin(ikl)  ! Update   thinest Lay.
     .                         + OKthin  * isn          !                Index
          dzthin(ikl) =    (1. - OKthin) * dzthin(ikl)  !
     .                         + OKthin  * dz_dif       !

! OUTPUT/Verification: Snow Thinest Layer  
! #sd   IF(isn.LE.isnoSV(1).AND.ikl.EQ.1)
! #sd.    write( 6,6001) isn,dzsnSV(ikl,isn),dz_min*isno_n*isno_n,dz_dif
! #sd.               ,OKthin,dzthin(ikl),    i_thin(ikl)
 6001     format(i3,5f12.6,i9)

        END DO
      END DO

! OUTPUT/Verification: Snow Thinest Layer  
! #sd     write( 6,*)   ' '

! OUTPUT/Verification: Snow Layers Agregation
! #vz     write(41,4153) i_thin(1),dzsnSV(1,i_thin(1))
 4153     format(/,' Thinest ',i3,':',f9.3)
! #vz     write(41,4151) i_thin(1),dzsnSV(1,i_thin(1))
! #vz.                  ,isnoSV(1),dzsnSV(1,isnoSV(1))

! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp   write(47,473)(G1snSV(1,isn),isn=1,isnoSV(1))
 473    format('Before _zCr2: G1 = ',10f8.1,(/,19x,10f8.1))
! #vp   write(47,472)(G2snSV(1,isn),isn=1,isnoSV(1))


! Index of the contiguous Layer to agregate
! -----------------------------------------

!          **********
      call SISVAT_zCr
!          **********


! Assign the 2 Layers to agregate
! -------------------------------

        DO ikl=1,klonv
          isn         =    i_thin(ikl)
          isagr1(ikl) =    istoSV(ikl,isn)
          isagr2(ikl) =    istoSV(ikl,isn+LIndsv(ikl))
          dzagr1(ikl) =    dzsnSV(ikl,isn)
          dzagr2(ikl) =    dzsnSV(ikl,isn+LIndsv(ikl))
          T_agr1(ikl) =    TsisSV(ikl,isn)
          T_agr2(ikl) =    TsisSV(ikl,isn+LIndsv(ikl))
          roagr1(ikl) =    ro__SV(ikl,isn)
          roagr2(ikl) =    ro__SV(ikl,isn+LIndsv(ikl))
          etagr1(ikl) =    eta_SV(ikl,isn)
          etagr2(ikl) =    eta_SV(ikl,isn+LIndsv(ikl))
          G1agr1(ikl) =    G1snSV(ikl,isn)
          G1agr2(ikl) =    G1snSV(ikl,isn+LIndsv(ikl))
          G2agr1(ikl) =    G2snSV(ikl,isn)
          G2agr2(ikl) =    G2snSV(ikl,isn+LIndsv(ikl))
          agagr1(ikl) =    agsnSV(ikl,isn)
          agagr2(ikl) =    agsnSV(ikl,isn+LIndsv(ikl))
          LstLay      = min(1,max(  0,    isnoSV(ikl)-1   ))
          Agrege(ikl) = min(1,
     .                  max(0,
     .                      NLaysv(ikl)   +isnoSV(ikl)-nsno
     .                     +NLay_s(ikl)                    )
     .                                    *LstLay           )
          isnoSV(ikl) =                    isnoSV(ikl)
     .     -(1-LstLay)*max(zer0,
     .                     sign(un_1,      eps_21
     .                                    -dzsnSV(ikl,1)   ))
          isnoSV(ikl) =max(   0,           isnoSV(ikl)      )

          WEagre(ikl) = 0.
        END DO

        DO isn=1,nsno
        DO ikl=1,klonv
          WEagre(ikl) = WEagre(ikl) + ro__SV(ikl,isn)*dzsnSV(ikl,isn)
     .                                *min(1,max(0,i_thin(ikl)+1-isn))
        ENDDO
        ENDDO

! OUTPUT/Verification: Snow Layers Agregation
! #vz write(41,4120)
 4120 format(' Agregation of too MUCH Layers')
! #vz write(41,411) (100.*dz_ref(  isn),isn=1,nsno)
! #vz write(41,412) (100.*dzwdif(  isn),isn=1,nsno)
! #vz write(41,413) (100.*dzsnSV(1,isn),isn=1,nsno)
! #vz write(41,414) (              isn ,isn=1,nsno)
! #vz write(41,4111)      isnoSV(1    )
! #vz write(41,4112)      i_thin(1    )
! #vz write(41,4113)      LIndsv(1    )
! #vz write(41,4114)      Agrege(1    )

! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp   write(47,474)(G1snSV(1,isn),isn=1,isnoSV(1))
 474    format('Before _zAg2: G1 = ',10f8.1,(/,19x,10f8.1))
! #vp   write(47,472)(G2snSV(1,isn),isn=1,isnoSV(1))


! Agregates
! ---------

!            **********
        call SISVAT_zAg
     .                 (isagr1,isagr2,WEagre
     .                 ,dzagr1,dzagr2,T_agr1,T_agr2
     .                 ,roagr1,roagr2,etagr1,etagr2
     .                 ,G1agr1,G1agr2,G2agr1,G2agr2
     .                 ,agagr1,agagr2,Agrege
     .                 )
!            **********


! Rearranges the Layers
! ---------------------

! New (agregated) Snow layer
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,klonv
          isn     =             i_thin(ikl)
          isn     = min(isn,isn+LIndsv(ikl))
          isnoSV(ikl) =         isnoSV(ikl) -Agrege(ikl)
          iiceSV(ikl) =         iiceSV(ikl) 
     .            -max(0,sign(1,iiceSV(ikl) -isn +icemix))
     .                                      *Agrege(ikl)
     .            *max(0,sign(1,iiceSV(ikl) -1          ))
          istoSV(ikl,isn) = (1.-Agrege(ikl))*istoSV(ikl,isn)
     .                      +   Agrege(ikl) *isagr1(ikl)
          dzsnSV(ikl,isn) = (1.-Agrege(ikl))*dzsnSV(ikl,isn)
     .                      +   Agrege(ikl) *dzagr1(ikl)
          TsisSV(ikl,isn) = (1.-Agrege(ikl))*TsisSV(ikl,isn)
     .                      +   Agrege(ikl) *T_agr1(ikl)
          ro__SV(ikl,isn) = (1.-Agrege(ikl))*ro__SV(ikl,isn)
     .                      +   Agrege(ikl) *roagr1(ikl)
          eta_SV(ikl,isn) = (1.-Agrege(ikl))*eta_SV(ikl,isn)
     .                      +   Agrege(ikl) *etagr1(ikl)
          G1snSV(ikl,isn) = (1.-Agrege(ikl))*G1snSV(ikl,isn)
     .                      +   Agrege(ikl) *G1agr1(ikl)
          G2snSV(ikl,isn) = (1.-Agrege(ikl))*G2snSV(ikl,isn)
     .                      +   Agrege(ikl) *G2agr1(ikl)
          agsnSV(ikl,isn) = (1.-Agrege(ikl))*agsnSV(ikl,isn)
     .                      +   Agrege(ikl) *agagr1(ikl)
        END DO

! Above
! ^^^^^
        DO ikl=1,klonv
          isn1(ikl)=max(i_thin(ikl),i_thin(ikl)+LIndsv(ikl))
        END DO
        DO i=  1,nsno-1
        DO ikl=1,klonv
            staggr        =  min(1,max(0,i +1 -isn1(ikl)   ))
            istoSV(ikl,i) = (1.-staggr     )*istoSV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*istoSV(ikl,i  )
     .                      +   Agrege(ikl) *istoSV(ikl,i+1))
            dzsnSV(ikl,i) = (1.-staggr     )*dzsnSV(ikl,i  ) 
     .            + staggr*((1.-Agrege(ikl))*dzsnSV(ikl,i  )
     .                      +   Agrege(ikl) *dzsnSV(ikl,i+1))
            TsisSV(ikl,i) = (1.-staggr     )*TsisSV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*TsisSV(ikl,i  )
     .                      +   Agrege(ikl) *TsisSV(ikl,i+1))
            ro__SV(ikl,i) = (1.-staggr     )*ro__SV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*ro__SV(ikl,i  )
     .                      +   Agrege(ikl) *ro__SV(ikl,i+1))
            eta_SV(ikl,i) = (1.-staggr     )*eta_SV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*eta_SV(ikl,i  )
     .                      +   Agrege(ikl) *eta_SV(ikl,i+1))
            G1snSV(ikl,i) = (1.-staggr     )*G1snSV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*G1snSV(ikl,i  )
     .                      +   Agrege(ikl) *G1snSV(ikl,i+1))
            G2snSV(ikl,i) = (1.-staggr     )*G2snSV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*G2snSV(ikl,i  )
     .                      +   Agrege(ikl) *G2snSV(ikl,i+1))
            agsnSV(ikl,i) = (1.-staggr     )*agsnSV(ikl,i  )
     .            + staggr*((1.-Agrege(ikl))*agsnSV(ikl,i  )
     .                      +   Agrege(ikl) *agsnSV(ikl,i+1))
        END DO
        END DO

        DO ikl=1,klonv
          isn             = min(isnoSV(ikl) +1,nsno)
          istoSV(ikl,isn) = (1.-Agrege(ikl))*istoSV(ikl,isn)
          dzsnSV(ikl,isn) = (1.-Agrege(ikl))*dzsnSV(ikl,isn)
          TsisSV(ikl,isn) = (1.-Agrege(ikl))*TsisSV(ikl,isn)
          ro__SV(ikl,isn) = (1.-Agrege(ikl))*ro__SV(ikl,isn)
          eta_SV(ikl,isn) = (1.-Agrege(ikl))*eta_SV(ikl,isn)
          G1snSV(ikl,isn) = (1.-Agrege(ikl))*G1snSV(ikl,isn)
          G2snSV(ikl,isn) = (1.-Agrege(ikl))*G2snSV(ikl,isn)
          agsnSV(ikl,isn) = (1.-Agrege(ikl))*agsnSV(ikl,isn)
        END DO

! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp   write(47,475)(G1snSV(1,isn),isn=1,isnoSV(1))
 475    format('At End _zSn : G1 = ',10f8.1,(/,19x,10f8.1))
! #vp   write(47,472)(G2snSV(1,isn),isn=1,isnoSV(1))


! Search new Ice/Snow Interface
! =============================

c #II   DO ikl=1,klonv
c #II     iiceSV(ikl) =  0
c #II   END DO

c #II   DO isn=1,nsno
c #II   DO ikl=1,klonv
c #II     OK_ICE      = max(zer0,sign(un_1,ro__SV(ikl,isn)-850.))
c #II.                * max(zer0,sign(un_1,dzsnSV(ikl,isn)-eps6))
c #II     iiceSV(ikl) = (1.-OK_ICE)       *iiceSV(ikl)
c #II.                +     OK_ICE        *isn
c #II   END DO
c #II   END DO

      return 
      end

                                                            
                                                            
      subroutine SISVAT_zCr                                                     
C +                                                                             
C +------------------------------------------------------------------------+    
C | MAR          SISVAT_zCr                                12-12-2002  MAR |    
C |   SubRoutine SISVAT_zCr determines criteria for Layers Agregation      |    
C |                                                                        |    
C +------------------------------------------------------------------------+    
C |                                                                        |    
C |   PARAMETERS:  klonv: Total Number of columns =                        |    
C |   ^^^^^^^^^^        = Total Number of continental     grid boxes       |    
C |                     X       Number of Mosaic Cell per grid box         |    
C |                                                                        |    
C |   INPUT /  isnoSV   = total Nb of Ice/Snow Layers                      |    
C |   OUTPUT:  iiceSV   = total Nb of Ice      Layers                      |    
C |   ^^^^^^   ispiSV   = 0,...,nsno: Uppermost Superimposed Ice Layer     |    
C |            istoSV   = 0,...,5 :   Snow     History (see istdSV data)   |    
C |                                                                        |    
C |   INPUT /  ro__SV   : Soil/Snow Volumic Mass                   [kg/m3] |    
C |   OUTPUT:           & Snow     Temperatures (layers  1,2,...,nsno) [K] |    
C |   ^^^^^^   G1snSV   : Dendricity (<0) or Sphericity (>0) of Snow Layer |    
C |            G2snSV   : Sphericity (>0) or Size            of Snow Layer |    
C |            agsnSV   : Snow       Age                             [day] |    
C |                                                                        |    
C |   OUTPUT:  LIndsv   : Relative Index of a contiguous Layer to agregate |    
C |   ^^^^^^                                                               |    
C +------------------------------------------------------------------------+    
C +                                                                             
C +                                                                             
      USE VAR_SV
      USE VARdSV
 
      USE VARySV                      
      USE VARxSV 
                                 
      USE VAR0SV                                                 
      USE VARphy         

      IMPLICIT NONE                                                       
C +                                                                             
C +                                                                             
C +--Global Variables                                                           
C +  ================                                                           
C +                                                                             
c      include  "LMDZphy.inc"                                                     
c      include  "LMDZ_SV.inc"                                                     
c      include  "LMDZdSV.inc"                                                     
c      include  "LMDZ0SV.inc"                                                     
C +                                                                             
c      include  "LMDZxSV.inc"                                                     
c      include  "LMDZySV.inc"                                                     
C +                                                                             
C +                                                                             
C +--Internal Variables                                                         
C +  ==================                                                         
C +                                                                             
      integer   ikl   ,isn   ,is0   ,is1                                        
      integer   isno_1                        ! Switch:  ! Snow Layer over Ice  
      real*8    Dtyp_0,Dtyp_1                 ! Snow Grains Difference Measure  
      real*8    DenSph                        ! 1. when contiguous spheric      
C +                                           !     and dendritic  Grains       
      real*8    DendOK                        ! 1. when dendritic  Grains       
      real*8    dTypMx                        ! Grain Type Differ.              
      real*8    dTypSp                        ! Sphericity Weight               
      real*8    dTypRo                        ! Density    Weight               
      real*8    dTypDi                        ! Grain Diam.Weight               
      real*8    dTypHi                        ! History    Weight               
                                                                                
                                                                                
C +--DATA                                                                       
C +  ====                                                                       
                                                                                
      data      dTypMx / 200.0  /             ! Grain Type Weight               
      data      dTypSp /   0.5  /             ! Sphericity Weight               
      data      dTypRo /   0.5  /             ! Density    Weight               
      data      dTypDi /  10.0  /             ! Grain Diam.Weight               
      data      dTypHi / 100.0  /             ! History    Weight               
                                                                                
                                                                                
C +--Agregation Criteria                                                        
C +  ===================                                                        
C +                                                                             
      DO  ikl=1,knonv                                                           
          i_thin(ikl) = min(i_thin(ikl),isnoSV(ikl))                            
          isn         = max(1          ,i_thin(ikl))                            
C +                                                                             
C +                                                                             
C +--Comparison with the downward Layer                                         
C +  ----------------------------------                                         
C +                                                                             
          is0    = max(1,        i_thin(ikl)-1 )        ! Downward Layer Index  
          DenSph = max(zero,                            ! isn/is1               
     .                 sign(unun,                       ! Dendricity/Sphericity 
     .                      epsi-G1snSV(ikl,isn)        !            Switch     
     .                          *G1snSV(ikl,is0)))      !                       
          DendOK = max(zero,                            ! Dendricity Switch     
     .                 sign(unun,                       !                       
     .                      epsi-G1snSV(ikl,isn)))      !                       
C +                                                                             
          Dtyp_0 =                                                              
     .         DenSph *      dTypMx                                             
     .    +(1.-DenSph)                                                          
     .    *    DendOK *((abs(G1snSV(ikl,isn)            ! Dendricity            
     .                      -G1snSV(ikl,is0))           !     Contribution      
     .                  +abs(G2snSV(ikl,isn)            ! Sphericity            
     .                      -G2snSV(ikl,is0))) *dTypSp  !     Contribution      
     .                  +abs(ro__SV(ikl,isn)            ! Density               
     .                      -ro__SV(ikl,is0))  *dTypRo) !     Contribution      
     .    +(1.-DenSph)                                  !                       
     .    *(1.-DendOK)*((abs(G1snSV(ikl,isn)            ! Sphericity            
     .                      -G1snSV(ikl,is0))           !     Contribution      
     .                  +abs(G2snSV(ikl,isn)            ! Size                  
     .                      -G2snSV(ikl,is0))) *dTypDi  !     Contribution      
     .                  +abs(ro__SV(ikl,isn)            ! Density               
     .                      -ro__SV(ikl,is0))  *dTypRo) !     Contribution      
          Dtyp_0 =                                      !                       
     .                   min(dTypMx,                    !                       
     .                       Dtyp_0                     !                       
     .                  +abs(istoSV(ikl,isn)            ! History               
     .                      -istoSV(ikl,is0))  *dTypHi) !     Contribution      
     .        +             (1 -abs(isn-is0))  * 1.e+6  !"Same Layer"Score      
     .        +  max(0,1-abs(iiceSV(ikl)                !"Ice /Snow             
     .                                 -is0))  * 1.e+6  ! Interface" Score      
C +                                                                             
C +                                                                             
C +--Comparison with the   upward Layer                                         
C +  ----------------------------------                                         
C +                                                                             
          is1    = min(          i_thin(ikl)+1,         ! Upward   Layer Index  
     .                 max(1,    isnoSV(ikl)  ))        !                       
          DenSph = max(zero,                            ! isn/is1               
     .                 sign(unun,                       ! Dendricity/Sphericity 
     .                      epsi-G1snSV(ikl,isn)        !            Switch     
     .                          *G1snSV(ikl,is1)))      !                       
          DendOK = max(zero,                            ! Dendricity Switch     
     .                 sign(unun,                       !                       
     .                      epsi-G1snSV(ikl,isn)))      !                       
C +                                                                             
          Dtyp_1 =                                                              
     .         DenSph *      dTypMx                                             
     .    +(1.-DenSph)                                                          
     .    *    DendOK *((abs(G1snSV(ikl,isn)            ! Dendricity            
     .                      -G1snSV(ikl,is1))           !     Contribution      
     .                  +abs(G2snSV(ikl,isn)            ! Sphericity            
     .                      -G2snSV(ikl,is1))) *dTypSp  !     Contribution      
     .                  +abs(ro__SV(ikl,isn)            ! Density               
     .                      -ro__SV(ikl,is1))  *dTypRo) !     Contribution      
     .    +(1.-DenSph)                                  !                       
     .    *(1.-DendOK)*((abs(G1snSV(ikl,isn)            ! Sphericity            
     .                      -G1snSV(ikl,is1))           !     Contribution      
     .                  +abs(G2snSV(ikl,isn)            ! Size                  
     .                      -G2snSV(ikl,is1))) *dTypDi  !     Contribution      
     .                  +abs(ro__SV(ikl,isn)            ! Density               
     .                      -ro__SV(ikl,is1))  *dTypRo) !     Contribution      
          Dtyp_1 =                                      !                       
     .                   min(dTypMx,                    !                       
     .                       Dtyp_1                     !                       
     .                  +abs(istoSV(ikl,isn)            ! History               
     .                      -istoSV(ikl,is1))  *dTypHi) !     Contribution      
     .        +             (1 -abs(isn-is1))  * 1.e+6  !"Same Layer"Score      
     .        +  max(0,1-abs(iiceSV(ikl)                !"Ice /Snow             
     .                                 -isn))  * 1.e+6  ! Interface" Score      
C +                                                                             
C +                                                                             
C +--Index of the Layer to agregate                                             
C +  ==============================                                             
C +                                                                             
          LIndsv(ikl) = sign(unun,Dtyp_0                                        
     .                           -Dtyp_1)                                       
          isno_1      = (1 -min (abs(isnoSV(ikl)        ! Switch = 1            
     .                              -iiceSV(ikl)-1),1)) !   if isno = iice +1   
     .                * (1 -min (abs(isnoSV(ikl)        ! Switch = 1            
     .                              -i_thin(ikl)  ),1)) !   if isno = i_ithin   
          LIndsv(ikl) = (1 -isno_1) *LIndsv(ikl)        ! Contiguous Layer is   
     .                     -isno_1                      ! downward for top L.   
          i_thin(ikl) =  max(1,   i_thin(ikl)   )                               
      END DO                                                                    
C +                                                                             
      return                                                                    
      end                                                                       
                                                                                
                                                                                
      subroutine SISVAT_zAg                                                     
     .                     (isagra,isagrb,WEagra                                
     .                     ,dzagra,dzagrb,T_agra,T_agrb                         
     .                     ,roagra,roagrb,etagra,etagrb                         
     .                     ,G1agra,G1agrb,G2agra,G2agrb                         
     .                     ,agagra,agagrb,Agreg1                                
     .                     )                                                    
                                                                                
C +------------------------------------------------------------------------+    
C | MAR SURFACE                                            17-06-2004  MAR |    
C |   SubRoutine SISVAT_zAg aggregates two contiguous snow layers          |    
C |                                                                        |    
C +------------------------------------------------------------------------+    
C |                                                                        |    
C |   PARAMETERS:  klonv: Total Number of columns =                        |    
C |   ^^^^^^^^^^        = Total Number of continental     grid boxes       |    
C |                     X       Number of Mosaic Cell per grid box         |    
C |                                                                        |    
C |   INPUT:   isagrb   : 2nd     Layer History                            |    
C |   ^^^^^                                                                |    
C |                                                                        |    
C |   INPUT:   dzagrb   : 2nd     Layer Thickness                          |    
C |   ^^^^^    T_agrb   : 2nd     Layer Temperature                        |    
C |            roagrb   : 2nd     Layer Density                            |    
C |            etagrb   : 2nd     Layer Water Content                      |    
C |            G1agrb   : 2nd     Layer Dendricity/Spher.                  |    
C |            G2agrb   : 2nd     Layer Sphericity/Size                    |    
C |            agagrb   : 2nd     Age                                      |    
C |            Agreg1   : 1. when Agregation constrained                   |    
C |                                                                        |    
C |   INPUT /  isagra   : 1st     Layer History                            |    
C |   OUTPUT:                                                              |    
C |   ^^^^^^                                                               |    
C |                                                                        |    
C |   INPUT /  dzagra   : 1st     Layer Thickness                          |    
C |   OUTPUT:  T_agra   : 1st     Layer Temperature                        |    
C |   ^^^^^^   roagra   : 1st     Layer Density                            |    
C |            etagra   : 1st     Layer Water Content                      |    
C |            G1agra   : 1st     Layer Dendricity/Spher.                  |    
C |            G2agra   : 1st     Layer Sphericity/Size                    |    
C |            agagra   : 1st     Age                                      |    
C |                                                                        |    
C +------------------------------------------------------------------------+    
            
      USE VAR_SV
      USE VARdSV                      
      USE VARxSV                 
      USE VAR0SV                                                 
      USE VARphy                                                                    
                                                                                
      IMPLICIT NONE                                                             
                                                                                
                                                                                
C +--Global Variables                                                           
C +  ================                                                           
                                                                                
c      include  "LMDZphy.inc"                                                     
c      include  "LMDZ_SV.inc"                                                     
c      include  "LMDZdSV.inc"                                                     
c      include  "LMDZ0SV.inc"                                                     
c      include  "LMDZxSV.inc"                                                     
                                                                                
                                                                                
C +--INPUT                                                                      
C +  -----                                                                      
                                                                                
      integer   isagrb(klonv)                 ! 2nd Layer History               
      real*8    dzagrb(klonv)                 ! 2nd Layer Thickness             
      real*8    T_agrb(klonv)                 ! 2nd Layer Temperature           
      real*8    roagrb(klonv)                 ! 2nd Layer Density               
      real*8    etagrb(klonv)                 ! 2nd Layer Water Content         
      real*8    G1agrb(klonv)                 ! 2nd Layer Dendricity/Spher.     
      real*8    G2agrb(klonv)                 ! 2nd Layer Sphericity/Size       
      real*8    agagrb(klonv)                 ! 2nd Layer Age                   
                                                                                
                                                                                
C +--INPUT/OUTPUT                                                               
C +  ------------                                                               
                                                                                
      integer   isagra(klonv)                 ! 1st Layer History               
      real*8    WEagra(klonv)                 ! 1st Layer Height    [mm w.e.]   
      real*8    Agreg1(klonv)                 ! 1. ===>   Agregates             
      real*8    dzagra(klonv)                 ! 1st Layer Thickness             
      real*8    T_agra(klonv)                 ! 1st Layer Temperature           
      real*8    roagra(klonv)                 ! 1st Layer Density               
      real*8    etagra(klonv)                 ! 1st Layer Water Content         
      real*8    G1agra(klonv)                 ! 1st Layer Dendricity/Spher.     
      real*8    G2agra(klonv)                 ! 1st Layer Sphericity/Size       
      real*8    agagra(klonv)                 ! 1st Layer Age                   
                                                                                
                                                                                
C +--Internal Variables                                                         
C +  ==================                                                         
                                                                                
      integer   ikl                                                             
      integer   nh                            ! Averaged    Snow History        
      integer   nh__OK                        ! 1=>Conserve Snow History        
      real*8    rh                            !                                 
      real*8    dz                            ! Thickness                       
      real*8    dzro_1                        ! Thickness X Density, Lay.1      
      real*8    dzro_2                        ! Thickness X Density, Lay.2      
      real*8    dzro                          ! Thickness X Density, Aver.      
      real*8    ro                            ! Averaged    Density             
      real*8    wn                            ! Averaged    Water Content       
      real*8    tn                            ! Averaged    Temperature         
      real*8    ag                            ! Averaged    Snow Age            
      real*8    SameOK                        ! 1. => Same Type of Grains       
      real*8    G1same                        ! Averaged G1,  same Grains       
      real*8    G2same                        ! Averaged G2,  same Grains       
      real*8    typ__1                        ! 1. => Lay1 Type: Dendritic      
      real*8    zroNEW                        ! dz X ro, if fresh Snow          
      real*8    G1_NEW                        ! G1,      if fresh Snow          
      real*8    G2_NEW                        ! G2,      if fresh Snow          
      real*8    zroOLD                        ! dz X ro, if old   Snow          
      real*8    G1_OLD                        ! G1,      if old   Snow          
      real*8    G2_OLD                        ! G2,      if old   Snow          
      real*8    SizNEW                        ! Size,    if fresh Snow          
      real*8    SphNEW                        ! Spheric.,if fresh Snow          
      real*8    SizOLD                        ! Size,    if old   Snow          
      real*8    SphOLD                        ! Spheric.,if old   Snow          
      real*8    Siz_av                        ! Averaged    Grain Size          
      real*8    Sph_av                        ! Averaged    Grain Spher.        
      real*8    Den_av                        ! Averaged    Grain Dendr.        
      real*8    DendOK                        ! 1. => Average is  Dendr.        
      real*8    G1diff                        ! Averaged G1, diff. Grains       
      real*8    G2diff                        ! Averaged G2, diff. Grains       
      real*8    G1                            ! Averaged G1                     
      real*8    G2                            ! Averaged G2                     
                                                                                
C +--OUTPUT of Snow Agregat.Statistics (see assignation in PHY_SISVAT)          
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                          
c #EV integer                iEVwri,jEVwri,nEVwri,kEVwri,lEVwri                 
c #EV common  /SISVAT_EV/    iEVwri,jEVwri,nEVwri,kEVwri,lEVwri                 
                                                                                
                                                                                
C +--Mean   Properties                                                          
C +  =================                                                          
                                                                                
C +-- 1 Densite, Contenu en Eau, Temperature /                                  
C +     Density, Water Content,  Temperature                                    
C +     ------------------------------------                                    
                                                                                
      DO ikl = 1,knonv                                                          
          dz      =  dzagra(ikl) + dzagrb(ikl)                                  
          dzro_1  =  roagra(ikl) * dzagra(ikl)                                  
          dzro_2  =  roagrb(ikl) * dzagrb(ikl)                                  
          dzro    =  dzro_1      + dzro_2                                       
          ro      =  dzro                                                       
     .     /max(epsi,dz)                                                        
          wn      = (dzro_1*etagra(ikl) + dzro_2*etagrb(ikl))                   
     .     /max(epsi,dzro)                                                      
          tn      = (dzro_1*T_agra(ikl) + dzro_2*T_agrb(ikl))                   
     .     /max(epsi,dzro)                                                      
          ag      = (dzro_1*agagra(ikl) + dzro_2*agagrb(ikl))                   
     .     /max(epsi,dzro)                                                      
                                                                                
          rh      =  max(zero,sign(unun,zWEcSV(ikl)-WEagra(ikl)))               
          nh__OK  =  rh                                                         
          nh      =  nh__OK *       max(isagra(ikl),isagrb(ikl))                
     .          + (1-nh__OK)*       min(isagra(ikl),isagrb(ikl))                
                                                                                
C +--OUTPUT of Snow Agregation Statistics                                       
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                       
c #EV   IF (ikl .EQ. kEVwri .AND. lEVwri .EQ. 3)                    THEN        
c #EV     write(6,5995) zWEcSV(ikl),WEagra(ikl)                                 
c #EV.                 ,isagra(ikl),isagrb(ikl)                                 
c #EV.                 ,nh__OK     ,nh                                          
 5995     format(' WE2,WEa =',2f9.1,'  nha,b =',2i2,'  nh__OK,nh =',2i2)        
c #EV   END IF                                                                  
                                                                                
                                                                                
C +-- 2 Nouveaux Types de Grains /  new Grain Types                             
C +     -------------------------------------------                             
                                                                                
C +-- 2.1. Meme  Type  de Neige  / same Grain Type                              
C +        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                              
          SameOK  =  max(zero,                                                  
     .                   sign(unun, G1agra(ikl) *G1agrb(ikl)  - eps_21))        
          G1same  = (dzro_1*G1agra(ikl) + dzro_2*G1agrb(ikl))                   
     .     /max(epsi,dzro)                                                      
          G2same  = (dzro_1*G2agra(ikl) + dzro_2*G2agrb(ikl))                   
     .     /max(epsi,dzro)                                                      
                                                                                
C +-- 2.2. Types differents / differents Types                                  
C +        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                  
          typ__1  =  max(zero,sign(unun,epsi-G1agra(ikl))) ! =1.=> Dendritic    
          zroNEW  =     typ__1  *dzro_1                    ! ro of Dendr.Lay.   
     .            + (1.-typ__1) *dzro_2                    !                    
          G1_NEW  =     typ__1  *G1agra(ikl)               ! G1 of Dendr.Lay.   
     .            + (1.-typ__1) *G1agrb(ikl)               !                    
          G2_NEW  =     typ__1  *G2agra(ikl)               ! G2 of Dendr.Lay.   
     .            + (1.-typ__1) *G2agrb(ikl)               !                    
          zroOLD  = (1.-typ__1) *dzro_1                    ! ro of Spher.Lay.   
     .            +     typ__1  *dzro_2                    !                    
          G1_OLD  = (1.-typ__1) *G1agra(ikl)               ! G1 of Spher.Lay.   
     .            +     typ__1  *G1agrb(ikl)               !                    
          G2_OLD  = (1.-typ__1) *G2agra(ikl)               ! G2 of Spher.Lay.   
     .            +     typ__1  *G2agrb(ikl)               !                    
          SizNEW  =    -G1_NEW  *DDcdSV/G1_dSV             ! Size  Dendr.Lay.   
     .             +(1.+G1_NEW         /G1_dSV)            !                    
     .                *(G2_NEW  *DScdSV/G1_dSV             !                    
     .             +(1.-G2_NEW         /G1_dSV)*DFcdSV)    !                    
          SphNEW  =     G2_NEW         /G1_dSV             ! Spher.Dendr.Lay.   
          SizOLD  =     G2_OLD                             ! Size  Spher.Lay.   
          SphOLD  =     G1_OLD         /G1_dSV             ! Spher.Spher.Lay.   
          Siz_av  = (zroNEW*SizNEW+zroOLD*SizOLD)      ! Averaged Size          
     .     /max(epsi,dzro)                             !                        
          Sph_av  = (zroNEW*SphNEW+zroOLD*SphOLD)      ! Averaged Sphericity    
     .     /max(epsi,dzro)                             !                        
          Den_av  = (Siz_av -(    Sph_av *DScdSV       !                        
     .                       +(1.-Sph_av)*DFcdSV))     !                        
     .            / (DDcdSV -(    Sph_av *DScdSV       !                        
     .                       +(1.-Sph_av)*DFcdSV))     !                        
          DendOK  = max(zero,                          !                        
     .                  sign(unun,     Sph_av *DScdSV  ! Small   Grains Contr.  
     .                            +(1.-Sph_av)*DFcdSV  ! Faceted Grains Contr.  
     .                            -    Siz_av        ))!                        
C +...    REMARQUE: le  type moyen (dendritique ou non) depend                  
C +       ^^^^^^^^  de la  comparaison avec le diametre optique                 
C +                 d'une neige recente de   dendricite nulle                   
C +...    REMARK:   the mean type  (dendritic   or not) depends                 
C +       ^^^^^^    on the comparaison with the optical diameter                
C +                 of a recent snow    having zero dendricity                  
                                                                                
          G1diff  =(   -DendOK *Den_av                                          
     .             +(1.-DendOK)*Sph_av) *G1_dSV                                 
          G2diff  =     DendOK *Sph_av  *G1_dSV                                 
     .             +(1.-DendOK)*Siz_av                                          
          G1      =     SameOK *G1same                                          
     .             +(1.-SameOK)*G1diff                                          
          G2      =     SameOK *G2same                                          
     .             +(1.-SameOK)*G2diff                                          
                                                                                
                                                                                
C +--Assignation to new Properties                                              
C +  =============================                                              
                                                                                
          isagra(ikl)   = Agreg1(ikl) *nh +(1.-Agreg1(ikl)) *isagra(ikl)        
          dzagra(ikl)   = Agreg1(ikl) *dz +(1.-Agreg1(ikl)) *dzagra(ikl)        
          T_agra(ikl)   = Agreg1(ikl) *tn +(1.-Agreg1(ikl)) *T_agra(ikl)        
          roagra(ikl)   = Agreg1(ikl) *ro +(1.-Agreg1(ikl)) *roagra(ikl)        
          etagra(ikl)   = Agreg1(ikl) *wn +(1.-Agreg1(ikl)) *etagra(ikl)        
          G1agra(ikl)   = Agreg1(ikl) *G1 +(1.-Agreg1(ikl)) *G1agra(ikl)        
          G2agra(ikl)   = Agreg1(ikl) *G2 +(1.-Agreg1(ikl)) *G2agra(ikl)        
          agagra(ikl)   = Agreg1(ikl) *ag +(1.-Agreg1(ikl)) *agagra(ikl)        
                                                                                
      END DO                                                                    
                                                                                
      return                                                                    
      end                                             



      subroutine SnOptP(jjtime)

!--------------------------------------------------------------------------+
!   MAR/SISVAT   SnOptP                               Sat 12-Feb-2012  MAR |
!     SubRoutine SnOptP computes the Snow Pack optical Properties          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     Grid Boxes       |
!                       X       Number of Mosaic Cell per Grid Box         |
!                                                                          |
!     INPUT:   isnoSV   = total Nb of Ice/Snow Layers                      |
!     ^^^^^    ispiSV   = 0,...,nsno: Uppermost Superimposed Ice Layer     |
!                                                                          |
!              ivgtSV   = 0,...,12:   Vegetation Type                      |
!                         0:          Water, Solid or Liquid               |
!                                                                          |
!     INPUT:   G1snSV   : Dendricity (<0) or Sphericity (>0) of Snow Layer |
!     ^^^^^    G2snSV   : Sphericity (>0) or Size            of Snow Layer |
!              agsnSV   : Snow       Age                             [day] |
!              ro__SV   : Snow/Soil  Volumic Mass                  [kg/m3] |
!              eta_SV   : Water      Content                       [m3/m3] |
!              rusnSV   : Surficial  Water   Thickness   [kg/m2] .OR. [mm] |
!              SWS_SV   : Surficial  Water   Status                        |
!              dzsnSV   : Snow       Layer   Thickness                 [m] |
!                                                                          |
!              albssv   : Soil       Albedo                            [-] |
!              zzsnsv   : Snow       Pack    Thickness                 [m] |
!                                                                          |
!     OUTPUT:  albisv   : Snow/Ice/Water/Soil Integrated Albedo        [-] |
!     ^^^^^^   sEX_sv   : Verticaly Integrated Extinction Coefficient      |
!                                                                          |
!     Internal Variables:                                                  |
!     ^^^^^^^^^^^^^^^^^^                                                   |
!              SnOpSV   : Snow Grain optical Size                      [m] |
!              EX1_sv   : Integrated Snow Extinction (0.3--0.8micr.m)      |
!              EX2_sv   : Integrated Snow Extinction (0.8--1.5micr.m)      |
!              EX3_sv   : Integrated Snow Extinction (1.5--2.8micr.m)      |
!                                                                          |
!     METHODE:    Calcul de la taille optique des grains ? partir de       |
!     ^^^^^^^    -leur type decrit par les deux variables descriptives     |
!                      continues sur la plage -99/+99 passees en appel.    |
!                -la taille optique (1/10mm) des etoiles,                  |
!                                            des grains fins et            |
!                                            des jeunes faces planes       |
!                                                                          |
!     METHOD:     Computation of the optical diameter of the grains        |
!     ^^^^^^      described with the CROCUS formalism G1snSV / G2snSV      |
!                                                                          |
!     REFERENCE: Brun et al.      1989, J. Glaciol 35 pp. 333--342         |
!     ^^^^^^^^^  Brun et al.      1992, J. Glaciol 38 pp.  13-- 22         |
!                Eric Martin Sept.1996                                     |
!                                                                          |
!     CAUTION: Vegetation is not taken into account in albedo computations |
!     ^^^^^^^  Suggestion: 1) Reduce the displacement height  and/or LAI   |
!              (when snow)    for radiative transfert through vegetation   |
!                          2) Adapt leaf optical parameters                |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Possibility                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                          |
!     #CZ: Albedo Correction (Zenith Angle) (Warren,      1982)            |
!     #cz: Albedo Correction (Zenith Angle) (Segal etAl., 1991) (obsolete) |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Col de Porte                         |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^                                      |
!     #cp: Col de Porte Integrated Snow/Ice Albedo                         |
!     #AG: Snow Aging Col de Porte     (Brun et al.1991)                   |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SnOptP____.va        | #va: OUTPUT/Verification: Albedo Parameteriz. |
!                          |      unit 46, SubRoutine  SnOptP     **ONLY** |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV
      USE VARdSV
      USE VARdCP

      USE VARxSV
      USE VARySV

      IMPLICIT NONE

! Internal Variables
! ==================

      real      coalb1(klonv)                      ! weighted Coalbedo, Vis.
      real      coalb2(klonv)                      ! weighted Coalbedo, nIR 1
      real      coalb3(klonv)                      ! weighted Coalbedo, nIR 2
      real      coalbm                             ! weighted Coalbedo, mean
      real      sExt_1(klonv)                      ! Extinction Coeff., Vis.
      real      sExt_2(klonv)                      ! Extinction Coeff., nIR 1
      real      sExt_3(klonv)                      ! Extinction Coeff., nIR 2
      real      SnOpSV(klonv,      nsno)           ! Snow Grain optical Size
c #AG real      agesno

      integer   isn   ,ikl   ,isn1  ,jjtime
      real      sbeta1,sbeta2,sbeta3,sbeta4,sbeta5
      real      AgeMax,AlbMin,HSnoSV,HIceSV,doptmx,SignG1,Sph_OK
      real      dalbed,dalbeS,dalbeW
c #CZ real      bsegal,czemax,csegal
      real      RoFrez,DiffRo,SignRo,SnowOK,OpSqrt
      real      albSn1,albIc1,a_SnI1,a_SII1!,alb1sv(klonv)
      real      albSn2,albIc2,a_SnI2,a_SII2!,alb2sv(klonv)
      real      albSn3,albIc3,a_SnI3,a_SII3!,alb3sv(klonv)
      real      albSno,albIce,albSnI,albSII,albWIc,albmax
      real      doptic,Snow_H,SIce_H,SnownH,SIcenH
      real      exarg1,exarg2,exarg3,sign_0,sExt_0
      real      albedo_old
      real      ro_ave,dz_ave

! OUTPUT/Verification: Albedo Parameteriz.
! #va logical         aw_opn                       ! IO       Switch
! #va common/SnOptP_L/aw_opn                       !


! Local   DATA
! ============

! For the computation of the solar irradiance extinction in snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      data      sbeta1/0.0192/,sbeta2/0.4000/,sbeta3/0.1098/
      data      sbeta4/1.0000/
      data      sbeta5/2.00e1/

! Snow Age Maximum (Taiga, e.g. Col de Porte)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      data      AgeMax  /60.0/
!               AgeMax:  Snow Age Maximum                              [day]

      data      AlbMin  /0.94/
!               AlbMin:  Albedo   Minimum / visible (0.3--0.8 micrometers)

      data      HSnoSV  /0.01/
!               HSnoSV:  Snow     Thickness through witch 
!                        Albedo is interpolated to Ice  Albedo
      data      HIceSV  /0.10/
!               HIceSV:  Snow/Ice Thickness through witch 
!                        Albedo is interpolated to Soil Albedo
      data      doptmx  /2.3e-3/
!               doptmx:  Maximum  optical Diameter    (pi * R**2)        [m]
!      
c #CZ data      czeMAX  /0.173648178/  ! 80.deg (Segal et al., 1991 JAS)  
c #CZ data      bsegal  /4.00       /  !
      data      albmax  /0.99       /  ! Albedo max


! Snow Grain optical Size
! =======================

      DO   isn=1,nsno
        DO ikl=1,knonv

          G2snSV(ikl,isn) =  max(eps6,G2snSV(ikl,isn))
!         Avoid non physical Values

          SignG1          = sign(un_1,G1snSV(ikl,isn))
          Sph_OK          =  max(zer0,SignG1)

          SnOpSV(ikl,isn) =   1.e-4 *
!         SI:           (from 1/10 mm to m)


! Contribution of Non Dendritic Snow
! ----------------------------------

     .    (    Sph_OK *(      G2snSV(ikl,isn)*G1snSV(ikl,isn)/G1_dSV
     .              +max(half*G2snSV(ikl,isn),DFcdSV)
     .                 *(1.00-G1snSV(ikl,isn)                /G1_dSV))


! Contribution of     Dendritic Snow
! ----------------------------------

     .    +(1.-Sph_OK)*(     -G1snSV(ikl,isn)*DDcdSV         /G1_dSV
     .                 +(1.00+G1snSV(ikl,isn)                /G1_dSV)
     .                  *    (G2snSV(ikl,isn)*DScdSV         /G1_dSV
     .                 +(1.00-G2snSV(ikl,isn)                /G1_dSV)
     .                                       *DFcdSV                 )))
          SnOpSV(ikl,isn) =  max(zer0,SnOpSV(ikl,isn))
        END DO
      END DO


! Snow/Ice Albedo
! ===============

! Snow Age (Influence on Albedo)
! ------------------------------

c #AG IF (iabs(mod(jjtime,86400)).lt.dt__SV)                        THEN
c #AG   DO isn=1,nsno
c #AG   DO ikl=1,knonv
c #AG     agsnSV(ikl,isn) = agsnSV(ikl,isn) + 1.
c #AG.           + max(zer0,DH_dSV(ivgtSV(ikl))-DH_dSV(4)) ! High Vegetation
!                                                          ! Impurities
!                           CAUTION: crude parameterization
!                           ^^^^^^^       
c #AG   END DO
c #AG   END DO
c #AG END IF


! Uppermost effective Snow Layer
! ------------------------------

        DO ikl=1,knonv

           isn   =  max(1,isnoSV(ikl))

          SignRo = sign(un_1, rocdSV - ro__SV(ikl,isn))
          SnowOK =  max(zer0,SignRo) ! Ice Density Threshold

          OpSqrt = sqrt(SnOpSV(ikl,isn))

          albSn1 =  0.96-1.580*OpSqrt
          albSn1 =  max(albSn1,AlbMin)

          albSn1 =  max(albSn1,zer0)
          albSn1 =  min(albSn1,un_1)

          albSn2 =  0.95-15.40*OpSqrt
          albSn2 =  max(albSn2,zer0)
          albSn2 =  min(albSn2,un_1)

          doptic =  min(SnOpSV(ikl,isn),doptmx)
          albSn3 =  346.3*doptic -32.31*OpSqrt +0.88
          albSn3 =  max(albSn3,zer0)
          albSn3 =  min(albSn3,un_1)

          albSno =  So1dSV*albSn1
     .           +  So2dSV*albSn2
     .           +  So3dSV*albSn3

          SnowOK =  SnowOK*max(zer0,sign(un_1,albSno-aI3dSV))
                 !  Minimum snow albedo is aI3dSV

          albSn1 =  SnowOK*albSn1+(1.0-SnowOK)*max(albSno,aI3dSV)
          albSn2 =  SnowOK*albSn2+(1.0-SnowOK)*max(albSno,aI3dSV)
          albSn3 =  SnowOK*albSn3+(1.0-SnowOK)*max(albSno,aI3dSV)


! Snow/Ice Pack Thickness
! -----------------------
 
          isn    =    max(min(isnoSV(ikl) ,ispiSV(ikl)),0)
          Snow_H =  zzsnsv(ikl,isnoSV(ikl))-zzsnsv(ikl,isn)
          SIce_H =  zzsnsv(ikl,isnoSV(ikl))
          SnownH =  Snow_H  /  HSnoSV
          SnownH =  min(un_1,  SnownH)
          SIcenH =  SIce_H  / (HIceSV
     .           +  max(zer0,Z0mdSV(ivgtSV(ikl))
     .           -           Z0mdSV(4)         ))
          SIcenH =  min(un_1,  SIcenH)
 
!      The value of SnownH is set to 1 in case of ice lenses above 
!      1m of dry snow (ro<700kg/m3) for using CROCUS albedo

          ro_ave =  0. 
          dz_ave =  0.
          SnowOK =  1.
       do isn    =  isnoSV(ikl),1,-1
          ro_ave =  ro_ave + ro__SV(ikl,isn) * dzsnSV(ikl,isn) * SnowOK
          dz_ave =  dz_ave +                   dzsnSV(ikl,isn) * SnowOK
          SnowOK =  max(zer0,sign(un_1,1.-dz_ave))
       enddo

          ro_ave =  ro_ave / max(dz_ave,eps6)
          SnowOK =  max(zer0,sign(un_1,700.-ro_ave))

          SnownH =  SnowOK + SnownH * (1. - SnowOK)   


! Integrated Snow/Ice Albedo: Case of Water on Bare Ice
! -----------------------------------------------------

          isn    =  max(min(isnoSV(ikl) ,ispiSV(ikl)),0)

          albWIc =  aI1dSV-(aI1dSV-aI2dSV)
     .           *  exp(-rusnSV(ikl)                       !
     .           *  (1. -SWS_SV(ikl)                       ! 0 <=> freezing
     .           *  (1  -min(1,iabs(isn-isnoSV(ikl)))))    ! 1 <=> isn=isnoSV
     .                  /ru_dSV)                           !

          SignRo = sign(un_1,rhoIce-1.-ro__SV(ikl,isn))    ! RoSN<920kg/m3
          SnowOK =  max(zer0,SignRo)

          albWIc = (1. - SnowOK) * albWIc + SnowOK   
     .           * (aI2dSV + (aI3dSV -aI2dSV) 
     .           * (ro__SV(ikl,isn)-rhoIce)/(rocdSV-rhoIce))
  
!    rocdSV < ro < rhoIce | aI2dSV< al >aI3dSV (fct of density)) 
!             ro > rhoIce | aI1dSV< al >aI2dSV (fct of superficial water content)s


! Integrated Snow/Ice      Albedo
! -------------------------------
 
          a_SII1      =     albWIc      +(albSn1-albWIc)     *SnownH
          a_SII1      = min(a_SII1       ,albSn1)
 
          a_SII2      =     albWIc      +(albSn2-albWIc)     *SnownH
          a_SII2      = min(a_SII2       ,albSn2)
 
          a_SII3      =     albWIc      +(albSn3-albWIc)     *SnownH
          a_SII3      = min(a_SII3       ,albSn3)

c #AG     agesno =      min(agsnSV(ikl,isn)          ,AgeMax)
c #AG     a_SII1      =     a_SII1      -0.175*agesno/AgeMax
!                                        Impurities: Col de Porte Parameter.


!    Zenith Angle Correction (Segal et al.,         1991, JAS 48, p.1025)
!    ----------------------- (Wiscombe & Warren, dec1980, JAS   , p.2723)
!                            (Warren,               1982,  RG   , p.  81)
!                            --------------------------------------------


          dalbed = 0.0
c #CZ     csegal = max(czemax                   ,coszSV(ikl))
c #cz     dalbeS =   ((bsegal+1.00)/(1.00+2.0*bsegal*csegal)
c #cz.                -       1.00                          )*0.32
c #cz.              /  bsegal
c #cz     dalbeS = max(dalbeS,zer0)
c #cz     dalbed =     dalbeS      *       min(1,isnoSV(ikl))

c #CZ     dalbeW =(0.64 - csegal  )*0.0625  ! Warren 1982, RevGeo, fig.12b
                                            ! 0.0625 = 5% * 1/0.8,   p.81
                                            ! 0.64   = cos(50)
c #CZ     dalbed =     dalbeW      *       min(1,isnoSV(ikl))

! Col de Porte Integrated Snow/Ice Albedo
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #cp     IF (ColPrt.AND.TotSol.gt.0.)                              THEN
c #cp       albSII =  (((Dr_1SN*a_SII1+Dr_2SN*a_SII2+Dr_3SN*a_SII3)
c #cp.                  +dalbed                                    )
c #cp.                  *DirSol
c #cp.                 +(Df_1SN*a_SII1+Df_2SN*a_SII2+Df_3SN*a_SII3)
c #cp.                  *DifSol*(1.   -cld_SV(ikl))
c #cp.                 +(Dfc1SN*a_SII1+Dfc2SN*a_SII2+Dfc3SN*a_SII3)
c #cp.                  *DifSol*       cld_SV(ikl)                  )
c #cp.                 / TotSol

! Elsewhere    Integrated Snow/Ice Albedo
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #cp     ELSE
            albSII =     So1dSV*a_SII1
     .                 + So2dSV*a_SII2
     .                 + So3dSV*a_SII3
c #cp     END IF


! Integrated Snow/Ice/Soil Albedo
! -------------------------------

            alb1sv(ikl) =     albssv(ikl) +(a_SII1-albssv(ikl))*SIcenH
            alb1sv(ikl) = min(alb1sv(ikl)  ,a_SII1)

            alb2sv(ikl) =     albssv(ikl) +(a_SII2-albssv(ikl))*SIcenH
            alb2sv(ikl) = min(alb2sv(ikl)  ,a_SII2)

            alb3sv(ikl) =     albssv(ikl) +(a_SII3-albssv(ikl))*SIcenH
            alb3sv(ikl) = min(alb3sv(ikl)  ,a_SII3)

            albisv(ikl) =     albssv(ikl) +(albSII-albssv(ikl))*SIcenH
            albisv(ikl) = min(albisv(ikl)  ,albSII)


! Integrated Snow/Ice/Soil Albedo: Clouds Correction! Greuell & all., 1994
! --------------------------------------------------! Glob.&t Planet.Change
                                                    ! (9):91-114
c #cp     IF (.NOT.ColPrt)                                          THEN
            alb1sv(ikl) = alb1sv(ikl) + 0.05 *(cld_SV(ikl)-0.5)*SIcenH
c #CZ.                  + dalbed      *    (1.-cld_SV(ikl))
            alb2sv(ikl) = alb2sv(ikl) + 0.05 *(cld_SV(ikl)-0.5)*SIcenH
c #CZ.                  + dalbed      *    (1.-cld_SV(ikl))
            alb3sv(ikl) = alb3sv(ikl) + 0.05 *(cld_SV(ikl)-0.5)*SIcenH
c #CZ.                  + dalbed      *    (1.-cld_SV(ikl))
            albisv(ikl) = albisv(ikl) + 0.05 *(cld_SV(ikl)-0.5)*SIcenH
c #CZ.                  + dalbed      *    (1.-cld_SV(ikl))
c #cp     END IF


! Integrated Snow/Ice/Soil Albedo: Minimum snow albedo = 40%
! ----------------------------------------------------------

            albedo_old  = albisv(ikl)

            albisv(ikl) = max(albisv(ikl),0.400    * SIcenH 
     .                  + albssv(ikl) *(1.0        - SIcenH))
            alb1sv(ikl) = alb1sv(ikl) - 1.0/3.0             ! 33 %
     .                  * (albedo_old-albisv(ikl)) / So1dSV
            alb2sv(ikl) = alb2sv(ikl) - 1.0/3.0             ! 33 %
     .                  * (albedo_old-albisv(ikl)) / So2dSV  
            alb3sv(ikl) = alb3sv(ikl) - 1.0/3.0             ! 33 %
     .                  * (albedo_old-albisv(ikl)) / So3dSV


! Integrated Snow/Ice/Soil Albedo: Maximum albedo = 99%
! -----------------------------------------------------

            albedo_old  = albisv(ikl)
            albisv(ikl) = min(albisv(ikl),0.99)
            alb1sv(ikl) = alb1sv(ikl) - 1.0/3.0             ! 33 %
     .                  * (albedo_old-albisv(ikl)) / So1dSV
            alb2sv(ikl) = alb2sv(ikl) - 1.0/3.0             ! 33 %
     .                  * (albedo_old-albisv(ikl)) / So2dSV  
            alb3sv(ikl) = alb3sv(ikl) - 1.0/3.0             ! 33 %
     .                  * (albedo_old-albisv(ikl)) / So3dSV

            alb1sv(ikl) = min(max(zer0,alb1sv(ikl)),albmax)  
            alb2sv(ikl) = min(max(zer0,alb2sv(ikl)),albmax)
            alb3sv(ikl) = min(max(zer0,alb3sv(ikl)),albmax)

        END DO


! Extinction Coefficient: Exponential Factor
! ==========================================

        DO ikl=1,knonv
          sExt_1(ikl)        = 1.
          sExt_2(ikl)        = 1.
          sExt_3(ikl)        = 1.
          sEX_sv(ikl,nsno+1) = 1.

          coalb1(ikl) = (1.          -alb1sv(ikl))*So1dSV
          coalb2(ikl) = (1.          -alb2sv(ikl))*So2dSV
          coalb3(ikl) = (1.          -alb3sv(ikl))*So3dSV
          coalbm      =  coalb1(ikl) +coalb2(ikl) +coalb3(ikl)
          coalb1(ikl) =  coalb1(ikl)              /coalbm
          coalb2(ikl) =  coalb2(ikl)              /coalbm
          coalb3(ikl) =  coalb3(ikl)              /coalbm
        END DO

      DO   isn=nsno,1,-1
        DO ikl=1,knonv

          SignRo = sign(un_1, rocdSV - ro__SV(ikl,isn))
          SnowOK =  max(zer0,SignRo) ! Ice Density Threshold

          RoFrez =  1.e-3      * ro__SV(ikl,isn) * (1.0-eta_SV(ikl,isn))

          OpSqrt = sqrt(max(eps6,SnOpSV(ikl,isn)))
          exarg1 =      SnowOK  *1.e2 *max(sbeta1*RoFrez/OpSqrt,sbeta2)
     .            +(1.0-SnowOK)           *sbeta5
          exarg2 =      SnowOK  *1.e2 *max(sbeta3*RoFrez/OpSqrt,sbeta4)
     .            +(1.0-SnowOK)           *sbeta5
          exarg3 =      SnowOK  *1.e2     *sbeta5
     .            +(1.0-SnowOK)           *sbeta5

! Col de Porte Snow Extinction Coefficient
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #cp     IF (ColPrt.AND.TotSol.gt.0.)                              THEN
c #cp       exarg1 = exarg1*(Dr_1SN*DirSol
c #cp.                      +Df_1SN*DifSol*(1.-cld_SV(ikl))
c #cp.                      +Dfc1SN*DifSol*    cld_SV(ikl) )
c #cp.                     /(Dr_1SN*TotSol)
c #cp       exarg2 = exarg2*(Dr_2SN*DirSol
c #cp.                      +Df_2SN*DifSol*(1.-cld_SV(ikl))
c #cp.                      +Dfc2SN*DifSol*    cld_SV(ikl) )
c #cp.                     /(Dr_2SN*TotSol)
c #cp       exarg3 = exarg3*(Dr_3SN*DirSol
c #cp.                      +Df_3SN*DifSol*(1.-cld_SV(ikl))
c #cp.                      +Dfc3SN*DifSol*    cld_SV(ikl) )
c #cp.                     /(Dr_3SN*TotSol)
c #cp     END IF


! Integrated Extinction of Solar Irradiance (Normalized Value)
! ============================================================

          sExt_1(ikl) = sExt_1(ikl)
     .                          * exp(min(0.0,-exarg1 *dzsnSV(ikl,isn)))
          sign_0      =              sign(un_1,epsn   -sExt_1(ikl))
          sExt_0      =               max(zer0,sign_0)*sExt_1(ikl)
          sExt_1(ikl) = sExt_1(ikl)                   -sExt_0

          sExt_2(ikl) = sExt_2(ikl)
     .                          * exp(min(0.0,-exarg2 *dzsnSV(ikl,isn)))
          sign_0      =              sign(un_1,epsn   -sExt_2(ikl))
          sExt_0      =               max(zer0,sign_0)*sExt_2(ikl)
          sExt_2(ikl) = sExt_2(ikl)                   -sExt_0

          sExt_3(ikl) = sExt_3(ikl)
     .                          * exp(min(0.0,-exarg3 *dzsnSV(ikl,isn)))
          sign_0      =              sign(un_1,epsn   -sExt_3(ikl))
          sExt_0      =               max(zer0,sign_0)*sExt_3(ikl)
          sExt_3(ikl) = sExt_3(ikl)                   -sExt_0

          sEX_sv(ikl,isn) = coalb1(ikl) * sExt_1(ikl)
     .                    + coalb2(ikl) * sExt_2(ikl)
     .                    + coalb3(ikl) * sExt_3(ikl)
        END DO
      END DO

      DO   isn=0,-nsol,-1
        DO ikl=1,knonv
          sEX_sv(ikl,isn) = 0.0
        END DO
      END DO


! Albedo: IO
! ==========

! #va IF (.NOT.aw_opn)                                              THEN
! #va          aw_opn = .true.
! #va          open(unit=46,status='unknown',file='SnOptP____.va')
! #va          rewind(   46)
! #va END IF

! #va     ikl=1
! #va     write(46,460)daHost
 460      format('---------------------------------+----+',
     .          '-------+-------+-------+-------+-------+-------+',
     .                                  '-------+-------+-------+',
     .         /,'Snow/Ice Pack ',a18,' |    |',
     .          ' z [m] |0.3/0.8|0.8/1.5|1.5/2.8| Full  |Opt[mm]|',
     .                                  ' G1    | G2    | ro    |',
     .         /,'---------------------------------+----+',
     .          '-------+-------+-------+-------+-------+-------+',
     .                                  '-------+-------+-------+')
!         ______________________________________________________________
! #va     write(46,461)            SIce_H,
! #va.                             alb1sv(ikl),alb2sv(ikl),alb3sv(ikl),
! #va.                             albisv(ikl)
 461      format('Integrated Snow/Ice/Soil  Albedo |',
     .            3x,' |',  f6.3,' |' ,4(f6.3,' |'), 6x ,' |',
     .                                            3( 6x ,' |'))
!         ______________________________________________________________
! #va     write(46,462)ispiSV(ikl),a_SII1,a_SII2,a_SII3,albSII
 462      format('Integrated Snow/Ice       Albedo |',
     .            i3,' |',   6x ,' |' ,4(f6.3,' |'), 6x ,' |',
     .                                            3( 6x ,' |'))
!         ______________________________________________________________
! #va     write(46,463)            rusnSV(ikl),         albWIc,
! #va.                             SWS_SV(ikl)
 463      format('Integrated Water/Bare Ice Albedo |',
     .            3x,' |',  f6.3,'w|' ,3( 6x, ' |'),  
     .                                   f6.3,' |' ,f6.3,' |',
     .                                            3( 6x ,' |'))
!         ______________________________________________________________
! #va     write(46,465)isn1       ,zzsnsv(ikl,isn1),
! #va.                             albIc1,albIc2,albIc3,albIce,
! #va.                        1.e3*SnOpSV(ikl,max(1,isnoSV(ikl)-1)),
! #va.                             G1snSV(ikl,max(1,isnoSV(ikl)-1)),
! #va.                             G2snSV(ikl,max(1,isnoSV(ikl)-1)),
! #va.                             ro__SV(ikl,max(1,isnoSV(ikl)-1))
! #va.                      *(1. - eta_SV(ikl,max(1,isnoSV(ikl)-1)))
 465      format('Surficial       Ice Lense        |',
     .            i3,' |', (f6.3,'i|'),4(f6.3,' |'),f6.3,' |',
     .                                            3(f6.1,' |'))
!         ______________________________________________________________
! #va     write(46,466)isnoSV(ikl),zzsnsv(ikl,isnoSV(ikl)),
! #va.                             albSn1,albSn2,albSn3,albSno,
! #va.                        1.e3*SnOpSV(ikl,isnoSV(ikl)),
! #va.                             G1snSV(ikl,isnoSV(ikl)),
! #va.                             G2snSV(ikl,isnoSV(ikl)),
! #va.                             ro__SV(ikl,isnoSV(ikl))
! #va.                      *(1. - eta_SV(ikl,isnoSV(ikl)))
 466      format('Uppermost  Effective Snow Layer  |',
     .            i3,' |', (f6.3,'*|'),4(f6.3,' |'),f6.3,' |',
     .                                            3(f6.1,' |'))

      return
      end


      subroutine VgOptP

!--------------------------------------------------------------------------+
!   MAR/SISVAT   VgOptP                               Sat 12-Feb-2012  MAR |
!     SubRoutine VgOptP computes the Canopy    optical Properties          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     Grid Boxes       |
!                       X       Number of Mosaic Cell per Grid Box         |
!                                                                          |
!     INPUT:   ivgtSV   = 0,...,12:   Vegetation Type                      |
!     ^^^^^               0:          Water, Solid or Liquid               |
!                                                                          |
!     INPUT:   coszSV   : Cosine of the Sun Zenithal Distance          [-] |
!     ^^^^^    sol_SV   : Surface Downward  Solar   Radiation       [W/m2] |
!              snCaSV   : Canopy     Snow      Thickness         [mm w.e.] |
!                                                                          |
!              LAI_sv   : Leaf Area  Index      (snow included)        [-] |
!              glf_sv   : Green Leaf Fraction of NOT fallen Leaves     [-] |
!              albisv   : Snow/Ice/Water/Soil Integrated Albedo        [-] |
!                                                                          |
!     OUTPUT:  alb_SV   : Surface-Canopy Albedo                        [-] |
!     ^^^^^^   SoCasv   : Absorbed Solar Radiation by Canopy (Normaliz)[-] |
!              SoSosv   : Absorbed Solar Radiation by Surfac (Normaliz)[-] |
!              LAIesv   : Effective Leaf Area  Index for Transpiration [-] |
!                                                                          |
!     Internal Variables: Normalized Values:                               |
!     ^^^^^^^^^^^^^^^^^^                                                   |
!              u0_Vis   : Upward   Visible Radiation at Top Canopy     [-] |
!              absg_V   : Absorbed Visible Radiation by the Ground     [-] |
!              absv_V   : Absorbed Visible Radiation by the Canopy     [-] |
!              u0_nIR   : Upward   Near IR Radiation at Top Canopy     [-] |
!              absgnI   : Absorbed Near IR Radiation by the Ground     [-] |
!              absv_V   : Absorbed Near IR Radiation by the Canopy     [-] |
!                                                                          |
!     REFERENCE:   De Ridder, 1997, unpublished thesis, chapter 2 (DR97,2) |
!     ^^^^^^^^^                                                            |
!                                                                          |
!     ASSUMPTIONS: Leaf Inclination Index chi_l (eqn2.49 DR97) set to zero |
!     ^^^^^^^^^^^                         for all vegetation types         |
!                  Radiation Fluxes are normalized                         |
!                        with respect to incoming solar radiation (=I0+D0) |
!                                                                          |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV

      USE VARxSV
      USE VARySV


      IMPLICIT NONE
 
! Internal Variables
! ==================

      integer   ikl   ,kri

      real      exdRad,k_drad,k___sv(klonv)
      real      e_prad,e1pRad
      real      zv_fac,zv1fac,deadLF
      real      T_Rad0,A_Rad0,A0__sv(klonv)
      real                    r0_Rad,t0_Rad,nu_Rad
      real      Tr_Rad,Re_Rad,r__Rad,t__Rad,t1_Rad
      real      arggam, gamma,gamasv(klonv),gammaL
      real      denSig,Sig__c,Sigcsv(klonv)
      real      DDifH1,DDifC1,C1__sv(klonv)
      real      DDifH2,DDifC2,C2__sv(klonv)
      real      denS_s,denS_a,den_c1,DDif_L
      real      u0_Vis,absg_V,absv_V
      real      u0_nIR,absgnI,absvnI
      real      argexg,argexk,criLAI(klonv)
      real      residu,d_DDif,dDDifs,dDDifa


! Internal DATA
! =============

      integer            nvgt
      parameter         (nvgt=12)
      real      reVisL(0:nvgt)      ! Reflectivity  / Visible / Live Leaves
      real      renIRL(0:nvgt)      ! Reflectivity  / Near IR / Live Leaves
      real      trVisL(0:nvgt)      ! Transmitivity / Visible / Live Leaves
      real      trnIRL(0:nvgt)      ! Transmitivity / Near IR / Live Leaves
      real      reVisD(0:nvgt)      ! Reflectivity  / Visible / Dead Leaves
      real      renIRD(0:nvgt)      ! Reflectivity  / Near IR / Dead Leaves
      real      trVisD(0:nvgt)      ! Transmitivity / Visible / Dead Leaves
      real      trnIRD(0:nvgt)      ! Transmitivity / Near IR / Dead Leaves

      real      reVisS              ! Reflectivity  / Visible / Canopy Snow
      real      renIRS              ! Reflectivity  / Near IR / Canopy Snow
      real      trVisS              ! Transmitivity / Visible / Canopy Snow
      real      trnIRS              ! Transmitivity / Near IR / Canopy Snow

      real      snCaMx              ! Canopy Snow Thickness for having Snow
                                    ! Snow Reflectivity and Transmitivity
      real      CriStR              ! Critical Radiation Stomatal Resistance

      integer            ivg

      DATA (reVisL(ivg),renIRL(ivg),trVisL(ivg),trnIRL(ivg),
     .      reVisD(ivg),renIRD(ivg),trVisD(ivg),trnIRD(ivg),ivg=0,nvgt)

!     reVisL renIRL trVisL trnIRL reVisD renIRD trVisD trnIRD  IGBP     CLASSES
!     ------ ------ ------ ------ ------ ------ ------ ------+ ----------------
     ./0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38,!  0 NO VEGETATION
     . 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38,!  1 CROPS LOW
     . 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38,!  2 CROPS MEDIUM
     . 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38,!  3 CROPS HIGH
     . 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38,!  4 GRASS LOW
     . 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38,!  5 GRASS MEDIUM
     . 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38,!  6 GRASS HIGH
     . 0.10,  0.45,  0.05,  0.25,  0.16,  0.39,  0.01,  0.01,!  7 BROADL LOW
     . 0.10,  0.45,  0.05,  0.25,  0.16,  0.39,  0.01,  0.01,!  8 BROADL MEDIUM
     . 0.10,  0.45,  0.05,  0.25,  0.16,  0.39,  0.01,  0.01,!  9 BROADL HIGH
     . 0.07,  0.35,  0.05,  0.10,  0.10,  0.39,  0.01,  0.01,! 10 NEEDLE LOW
     . 0.07,  0.35,  0.05,  0.10,  0.10,  0.39,  0.01,  0.01,! 11 NEEDLE MEDIUM
     . 0.07,  0.35,  0.05,  0.10,  0.10,  0.39,  0.01,  0.01/! 12 NEEDLE HIGH

      DATA 
     .reVisS,renIRS,trVisS,trnIRS
!     ------ ------ ------ ------+
     ./0.85,  0.85,  0.00,  0.00/!
!     REMARK: Possible Refinement by taking actual Surface Snow Reflectivities
!     ^^^^^^

      DATA      snCaMx /0.5/

      DATA      CriStR /25./


! General Parameters, Solar Radiation Absorption
! ==============================================

      DO ikl=1,knonv

            k_dRad = 0.5 /max(coszSV(ikl),eps6) ! absorbed irradiance fraction
            e_pRad = 2.5   *  coszSV(ikl)       ! exponential argument,   
                                                ! V/nIR radiation partitioning,
                                                ! DR97, 2, eqn (2.53) & (2.54)
            exdRad =    exp(-k_dRad*LAI_sv(ikl))! exponential, Irradi. Absorpt.
            e1pRad = 1.-exp(-e_pRad)            ! exponential, V/nIR Rad. Part.

            ivg    =                ivgtSV(ikl) ! Vegetation Type
            zv_fac =    min( snCaSV(ikl)/snCaMx ! Contribution of Snow to Leaf
     .                      ,  un_1)            ! Reflectivity and Transmissiv.
            zv1fac = 1.     -       zv_fac      !
            deadLF = 1.     -       glf_sv(ikl) ! Dead Leaf Fraction


! Visible Part of the Solar Radiation Spectrum (V,   0.4--0.7mi.m)
! ================================================================

            A_Rad0 =      0.25 + 0.697 * e1pRad ! Absorbed    Vis. Radiation
            T_Rad0 = 1. - A_Rad0                ! Transmitted Vis  Radiation

! Reflectivity, Transmissivity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            Re_Rad = glf_sv(ikl) *ReVisL(ivg) 
     .             + deadLF      *ReVisD(ivg)
            Tr_Rad = glf_sv(ikl) *TrVisL(ivg) 
     .             + deadLF      *TrVisD(ivg)

! Adaptation to Snow
! ^^^^^^^^^^^^^^^^^^
            Re_Rad = zv1fac      *Re_Rad      + zv_fac *reVisS
            Tr_Rad = zv1fac      *Tr_Rad      + zv_fac *trVisS

! Scattering /DR97, 2, eqn (2.26) and (2.27)           ! Diffuse  Radiation:
! ^^^^^^^^^^                                           ! ^^^^^^^^^^^^^^^^^^
            r__Rad = (2. *Re_Rad +     Tr_Rad) / 3.    ! Upw.  Scatter.Fract.
            t__Rad = (    Re_Rad + 2. *Tr_Rad) / 3.    ! Downw.Scatter.Fract.

            t1_Rad =  1. -t__Rad                       !
            arggam =      t1_Rad*t1_Rad-r__Rad*r__Rad  !
            arggam =  max(arggam,zer0)                 !
            gamma  = sqrt(arggam)                      ! eqn (2.39)
            gammaL =  min( gamma*LAI_sv(ikl),40.0)     !
            DDifH1 =  exp( gammaL           )          ! Downw.Diffus.Solut.1
            DDifH2 =  exp(-gammaL           )          ! Downw.Diffus.Solut.2
!           REMARK:  These 2 contributions are zero in case of 0 Reflectivity
!           ^^^^^^

! Scattering /DR97, 2, eqn (2.19) and (2.20)           ! Direct   Radiation:
! ^^^^^^^^^^                                           ! ^^^^^^^^^^^^^^^^^^
            r0_Rad = 0.5 *((Re_Rad+Tr_Rad) *k_dRad     ! Upw.  Scatter.Fract.
     .                    +(Re_Rad-Tr_Rad) /    3.)    ! 
            t0_Rad = 0.5 *((Re_Rad+Tr_Rad) *k_dRad     ! Downw.Scatter.Fract.
     .                    -(Re_Rad-Tr_Rad) /    3.)    !

            nu_Rad = t1_Rad-r__Rad*albisv(ikl)         ! nu coeff., eqn 2.43
            den_c1 =  gamma*(DDifH1+DDifH2)            ! eqn (2.43) Denomin. 
     .              +nu_Rad*(DDifH1-DDifH2)            !(Constant for DDifH1)

            denSig =  gamma*gamma - k_dRad*k_dRad      ! eqn (2.40) Denomin.
            denS_s = sign(un_1,denSig)                 !
            denS_a =  abs(     denSig)                 !
            denSig =  max(eps6,denS_a) * denS_s        !
            Sig__c = (r__Rad* r0_Rad                   ! sigma_c, eqn (2.40)
     .               +t0_Rad*(k_dRad+t1_Rad)) / denSig ! 

            DDifC1 = ((gamma-nu_Rad)*(T_Rad0-Sig__c*A_Rad0)*DDifH2
     .             +((k_dRad-nu_Rad)* Sig__c
     .               +t0_Rad+r__Rad * albisv(ikl)) *A_Rad0 *exdRad)
     .           /max(den_c1,eps6)
            DDifC2 =  T_Rad0        - DDifC1-Sig__c*A_Rad0

! Visible Diffuse Fluxes
! ^^^^^^^^^^^^^^^^^^^^^^
            DDif_L =  DDifC1*DDifH1 + DDifC2*DDifH2    ! DOWNward,
     .             +  Sig__c*A_Rad0 *exdRad            ! Canopy Basis
            u0_Vis = ((gamma+t1_Rad)*DDifC1            ! UPward
     .               -(gamma-t1_Rad)*DDifC2            ! Canopy Top
     .             -((k_dRad-t1_Rad)*Sig__c            !
     .               +t0_Rad               )*A_Rad0)   !
     .          / max(r__Rad,eps6)                     !
            u0_Vis = min(0.99,max(eps6,u0_Vis))        ! ERROR 
            absg_V = (1.-albisv(ikl))*(A_Rad0*exdRad   ! Ground Absorption
     .                                +DDif_L       )  !
            absv_V = (1.-u0_Vis     )- absg_V          ! Veget. Absorption

! Parameters for Computing Effective LAI for Transpiration
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            gamasv(ikl) = gamma
            C1__sv(ikl) = DDifC1
            C2__sv(ikl) = DDifC2
            Sigcsv(ikl) = Sig__c
            k___sv(ikl) = k_dRad
            A0__sv(ikl) = A_Rad0


! Near-IR Part of the Solar Radiation Spectrum (nIR, 0.7--2.8mi.m)
! ================================================================

            A_Rad0 =      0.80 + 0.185 * e1pRad ! Absorbed    nIR. Radiation
            T_Rad0 = 1. - A_Rad0                ! Transmitted nIR  Radiation

! Reflectivity, Transmissivity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            Re_Rad = glf_sv(ikl) *RenIRL(ivg) 
     .             + deadLF      *RenIRD(ivg)
            Tr_Rad = glf_sv(ikl) *TrnIRL(ivg) 
     .             + deadLF      *TrnIRD(ivg)

! Adaptation to Snow
! ^^^^^^^^^^^^^^^^^^
            Re_Rad = zv1fac      *Re_Rad      + zv_fac *renIRS
            Tr_Rad = zv1fac      *Tr_Rad      + zv_fac *trnIRS

! Scattering /DR97, 2, eqn (2.26) and (2.27)           ! Diffuse  Radiation:
! ^^^^^^^^^^                                           ! ^^^^^^^^^^^^^^^^^^
            r__Rad = (2. *Re_Rad +     Tr_Rad) / 3.    ! Upw.  Scatter.Fract.
            t__Rad = (    Re_Rad + 2. *Tr_Rad) / 3.    ! Downw.Scatter.Fract.

            t1_Rad =  1. -t__Rad                       !
            arggam =      t1_Rad*t1_Rad-r__Rad*r__Rad  !
            arggam =  max(arggam,zer0)                 !
            gamma  = sqrt(arggam)                      ! eqn (2.39)
            DDifH1 =  exp( gamma*LAI_sv(ikl))          ! Downw.Diffus.Solut.1
            DDifH2 =  exp(-gamma*LAI_sv(ikl))          ! Downw.Diffus.Solut.2
!           REMARK:  These 2 contributions are zero in case of 0 Reflectivity
!           ^^^^^^

! Scattering /DR97, 2, eqn (2.19) and (2.20)           ! Direct   Radiation:
! ^^^^^^^^^^                                           ! ^^^^^^^^^^^^^^^^^^
            r0_Rad = 0.5 *((Re_Rad+Tr_Rad) *k_dRad     ! Upw.  Scatter.Fract.
     .                    +(Re_Rad-Tr_Rad) /    3.)    ! 
            t0_Rad = 0.5 *((Re_Rad+Tr_Rad) *k_dRad     ! Downw.Scatter.Fract.
     .                    -(Re_Rad-Tr_Rad) /    3.)    !

            nu_Rad = t1_Rad-r__Rad*albisv(ikl)         ! nu coeff., eqn 2.43
            den_c1 =  gamma*(DDifH1+DDifH2)            ! eqn (2.43) Denomin. 
     .              +nu_Rad*(DDifH1-DDifH2)            !(Constant for DDifH1)

            denSig =  gamma*gamma - k_dRad*k_dRad      ! eqn (2.40) Denomin.
            denS_s = sign(un_1,denSig)                 !
            denS_a =  abs(     denSig)                 !
            denSig =  max(eps6,denS_a) * denS_s        !
            Sig__c = (r__Rad* r0_Rad                   ! sigma_c, eqn (2.40)
     .               +t0_Rad*(k_dRad+t1_Rad)) / denSig ! 

            DDifC1 = ((gamma-nu_Rad)*(T_Rad0-Sig__c*A_Rad0)*DDifH2
     .             +((k_dRad-nu_Rad)* Sig__c
     .               +t0_Rad+r__Rad * albisv(ikl)) *A_Rad0 *exdRad)
     .           /max(den_c1,eps6)
            DDifC2 =  T_Rad0        - DDifC1-Sig__c*A_Rad0

! Near IR Diffuse Fluxes
! ^^^^^^^^^^^^^^^^^^^^^^
            DDif_L =  DDifC1*DDifH1 + DDifC2*DDifH2    ! DOWNward,
     .             +  Sig__c*A_Rad0 *exdRad            ! Canopy Basis
            u0_nIR = ((gamma+t1_Rad)*DDifC1            ! UPward
     .               -(gamma-t1_Rad)*DDifC2            ! Canopy Top
     .             -((k_dRad-t1_Rad)*Sig__c            !
     .               +t0_Rad               )*A_Rad0)   !
     .          / max(r__Rad,eps6)                     !
            u0_nIR = min(0.99,max(eps6,u0_nIR))        ! ERROR
            absgnI = (1.-albisv(ikl))*(A_Rad0*exdRad   ! Ground Absorption
     .                                +DDif_L       )  !
            absvnI = (1.-u0_nIR     )- absgnI          ! Veget. Absorption


! Surface-Canopy Albedo and Normalized Solar Radiation Absorption
! ===============================================================

            alb_SV(ikl) = (u0_Vis+u0_nIR)*0.5d0
            SoCasv(ikl) = (absv_V+absvnI)*0.5d0
            SoSosv(ikl) = (absg_V+absgnI)*0.5d0

      END DO


! Effective LAI for Transpiration
! ===============================

        DO ikl=1,knonv
              criLAI(ikl) = 2.                  ! LAI for which D0_Vis > 20W/m2
                                                ! DR97, 2, eqn (2.57)
        END DO

      DO   kri=1,10
        DO ikl=1,knonv

              argexg      =  min(criLAI(ikl)*gamasv(ikl),     ea_Max)
              argexk      =  min(criLAI(ikl)*k___sv(ikl),     ea_Max)
              residu      =      C1__sv(ikl)            *exp( argexg)
     .                          +C2__sv(ikl)            *exp(-argexg)
     .                          +A0__sv(ikl)*gamasv(ikl)*exp(-argexk)
     .                          -CriStR /max(sol_SV(ikl),       eps6)

              d_DDif      =      C1__sv(ikl)*gamasv(ikl)*exp( argexg)
     .                          -C2__sv(ikl)*gamasv(ikl)*exp(-argexg)
     .                          -A0__sv(ikl)*k___sv(ikl)*exp(-argexk)
              dDDifs      = sign(un_1,d_DDif)
              dDDifa      =  abs(     d_DDif)
              d_DDif      =  max(eps6,dDDifa) * dDDifs

              criLAI(ikl) =      criLAI(ikl)-residu/d_DDif
              criLAI(ikl) =  max(criLAI(ikl),zer0       )
              criLAI(ikl) =  min(criLAI(ikl),LAI_sv(ikl))

        END DO
      END DO

        DO ikl=1,knonv
              LAIesv(ikl) = criLAI(ikl) +(exp(-k___sv(ikl)*criLAI(ikl))
     .                                   -exp(-k___sv(ikl)*LAI_sv(ikl)))
     .                                  /      k___sv(ikl)
        END DO

      return
      end


      subroutine ColPrt_SBL

!--------------------------------------------------------------------------+
!   MAR          ColPrt_SBL                           Sat 12-Feb-2012  MAR |
!     SubRoutine ColPrt_SBL generates Surface Boundary Layers Properties   |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns                          |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   za__SV   : Surface Boundary Layer (SBL) Height          [m] |
!     ^^^^^    VV__SV   :(SBL Top)   Wind Velocity                   [m/s] |
!              TaT_SV   : SBL Top    Temperature                       [K] |
!              rhT_SV   : SBL Top  Air  Density                    [kg/m3] |
!              uqs_SV   : Specific  Humidity  Turbulent Flux         [m/s] |
!              Tsrfsv   : Surface    Temperature                       [K] |
!                                                                          |
!     INPUT /  LMO_SV   : Monin-Obukhov       Scale                    [m] |
!     OUTPUT:  us__SV   : Friction  Velocity                         [m/s] |
!     ^^^^^^   uts_SV   : Temperature         Turbulent Flux       [K.m/s] |
!                                                                          |
!     OUTPUT:  ram_sv   : Aerodynamic Resistance for Momentum        [s/m] |
!     ^^^^^^   rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!                                                                          |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV
      USE VARdSV

      USE VARxSV
      USE VARySV

      IMPLICIT NONE


! Internal Variables
! ==================

      integer   ikl    ,ist    ,ist__s ,ist__w
      real      d_TaTs ,CD_m
      real      uustar ,thstar ,qqstar ,ssstar
      real      thstarv,thstars,thstara
      real      zeta   ,zeta_S ,zeta_A
      real      fCdCdP ,Cd_min ,cCdUns
      real      RapCm0


! Internal DATA
! =============

      data      fCdCdP/ 3.09/   ! Drag Coefficient Factor, Col de Porte
      data      Cd_min/ 1.05/   ! Drag Coefficient Minimum Col de Porte
      data      cCdUns/-5.00/   ! Drag Coefficient Correction for Unstability


! Aerodynamic Resistances
! =======================

      DO ikl=1,knonv

! Surface Type
! ~~~~~~~~~~~~
        ist    =      isotSV(ikl)                    ! Soil Type
        ist__s =  min(ist, 1)                        ! 1 => Soil
        ist__w =  1 - ist__s                         ! 1 => Water Body

! Drag and Aerodynamic Resistance
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        d_TaTs =      TaT_SV(ikl)-Tsrfsv(ikl)
        RapCm0 =  log(za__SV(ikl)/Z0mdSV(4          ))
     .         /  log(za__SV(ikl)/Z0mdSV(ivgtSV(ikl)))
        RapCm0 =      RapCm0     *RapCm0             ! Neutral Drag Coefficient
                                                     ! Vegetation   Correction
        CD_m   =  max(Cd_min*RapCm0,                 ! Actual  Drag Coefficient
     .                fCdCdP*RapCm0*VV__SV(ikl)  )   !         for  Momentum
     .          *(1.+max(min(d_TaTs,zer0),cCdUns)    ! Unstability  Correction
     .                                   /cCdUns )
     .          * 1.5
        ram_sv(ikl) = rhT_SV(ikl) *CpdAir/CD_m
        rah_sv(ikl) = ram_sv(ikl)


! Turbulent Scales
! ================

! Friction Velocity                   u*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        uustar      =      VV__SV(ikl) / ram_sv(ikl)
        us__SV(ikl) = sqrt(uustar)

! Real    Temperature Turbulent Scale theta*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        uts_SV(ikl) =    - d_TaTs      / rah_sv(ikl)
        thstar      =      uts_SV(ikl) / us__SV(ikl)

! Specific Humidity   Turbulent Scale qq*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        qqstar      =      uqs_SV(ikl) / us__SV(ikl) 

! Virtual Temperature Turbulent Scale thetav*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        thstarv     = thstar      + TaT_SV(ikl) *(0.608*qqstar)
        thstars     =     sign(un_1,thstarv)
        thstara     =      abs(     thstarv)
        thstarv     =      max(eps6,thstara)    *thstars

! Monin Obukhov Scale Height
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        LMO_SV(ikl) = TaT_SV(ikl) * uustar 
     .              /(vonKrm      * Grav_F     * thstarv)
        zeta        = za__SV(ikl) / LMO_SV(ikl)
        zeta_S      =          sign(un_1  ,zeta)
        zeta_A      =           abs(       zeta)
        zeta        = zeta_S  * max(eps6  ,zeta_A)
        LMO_SV(ikl) = za__SV(ikl) / zeta

      END DO

      return
      end

      subroutine SISVATeSBL

!--------------------------------------------------------------------------+
!   MAR          SISVATeSBL                           Tue 12-Apr-2011  MAR |
!     SubRoutine SISVATeSBL generates Surface Boundary Layers Properties   |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns                          |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   za__SV   : Surface Boundary Layer (SBL) Height          [m] |
!     ^^^^^    VV__SV   :(SBL Top)   Wind Velocity                   [m/s] |
!              TaT_SV   : SBL Top    Temperature                       [K] |
!              ExnrSV   : Exner      Potential                         [-] |
!              qsnoSV   : SBL Mean   Snow      Content             [kg/kg] |
!              uqs_SV   : Specific   Humidity  Turbulent Flux        [m/s] |
!              usthSV   : Blowing Snow Erosion   Threshold           [m/s] |
!              Z0m_SV   : Momentum     Roughness Length                [m] |
!              Z0h_SV   : Heat         Roughness Length                [m] |
!              Tsrfsv   : Surface    Temperature                       [K] |
!              sqrCm0   : Contribution of Z0m to Neutral Drag Coefficient  |
!              sqrCh0   : Contribution of Z0h to Neutral Drag Coefficient  |
!                                                                          |
!     INPUT /  LMO_SV   : Monin-Obukhov       Scale                    [m] |
!     OUTPUT:  us__SV   : Friction  Velocity                         [m/s] |
!     ^^^^^^   uts_SV   : Temperature         Turbulent Flux       [K.m/s] |
!              uss_SV   : Blowing Snow        Turbulent Flux         [m/s] |
!                                                                          |
!     OUTPUT:  hSalSV   : Saltating Layer Height                       [m] |
!     ^^^^^^   qSalSV   : Saltating Snow  Concentration            [kg/kg] |
!              ram_sv   : Aerodynamic Resistance for Momentum        [s/m] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Possibility                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                          |
!     #AE: TURBULENCE: Aerosols Erosion / Turbulent Diffusion Coeff.       |
!                                                                          |
!     #AW  TURBULENCE: Wind Time Mean (BOX Moving Average)                 |
!     #AH  TURBULENCE: Ta-T Time Mean (BOX Moving Average)                 |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: OBSOLETE                                      |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^                                      |
!     #AM  TURBULENCE: u*   Time Mean (BOX Moving Average)                 |
!     #AT  TURBULENCE: u*T* Time Mean (BOX Moving Average)                 |
!     #AS  TURBULENCE: u*s* Time Mean (BOX Moving Average)                 |
!                                                                          |
!     #ZX  TURBULENCE: Strong Stability Limit (King    et al. 1996)        |
!     #zx  TURBULENCE: Strong Stability Limit (Mahalov et al. 2004)        |
!     #IX  TURBULENCE: recurrence                                          |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # stdout               | #ss: OUTPUT of Blowing Snow Variables         |
!                          |      unit  6, SubRoutine  SISVATeSBL **ONLY** |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV
      USE VARdSV

      USE VARxSV
      USE VARySV
      USE VARtSV

      IMPLICIT NONE

! V,  dT(a-s)    Time Moving Averages
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cc #AA     integer                       ntaver,nt !
cc #AA      parameter                    (ntaver= 4)! ntaver defined in MAR_SL.inc
c #AW real             V__mem(klonv,ntaver)   !                   only
c #AW real             VVmmem(klonv)          !
c #AW common/SVeSBLmem/V__mem,VVmmem          !
c #AH real             T__mem(klonv,ntaver)   !
c #AH real             dTmmem(klonv)          !
c #AH common/STeSBLmem/T__mem,dTmmem          !

! u*, u*T*, u*s* Time Moving Averages
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #AM real             u__mem(klonv,ntaver)
c #AT real             uT_mem(klonv,ntaver)
c #AS real             us_mem(klonv,ntaver)
c #AM common/S_eSBLmem/u__mem
c #AT.                ,uT_mem
c #AS.                ,us_mem


! Internal Variables
! ==================

      integer   ikl   ,icount
c #AE integer   nit   ,iit
      real      VVaSBL(klonv),VVa_OK          ! effective SBL wind speed
      real      dTa_Ts(klonv)                 ! effective SBL Temperature diff.
      real      Theta0                        ! Potential Reference Temperature
      real      LMOmom(klonv)                 ! Monin-Obukhov Scale Momentum
      real      LMOsgn                        ! Monin-Obukhov Scale Sign
      real      LMOabs                        ! Monin-Obukhov Scale Abs.Value
      real      uustar,thstar,qqstar,ssstar,thstarv,thstars,thstara
      real      zetam ,zetah ,zeta_S,zeta_A,zeta0m ,zeta0h
      real      psim_s,xpsimi,psim_i,psim_z
      real      psis_s,psis_z,psis_0
      real      psih_s,xpsihi,psih_i,psih_z
      real      psim_0,psih_0
      real      CDm(klonv)                    ! Drag Coefficient, Momentum
      real      CDs(klonv),rCDs(klonv)        ! Drag Coefficient, Blown **
      real      CDh(klonv)                    ! Drag Coefficient, Scalar
      real      dustar,u0star,uTstar,usstar
      real      sss__F,sss__N,usuth0
c #AE real      dusuth,signus
c #AE real      sss__K,sss__G
c #AE real      us_127,us_227,us_327,us_427,us_527
      real      zetMAX
      real      coef_m,coef_h,stab_s
c #AE real      SblPom
      real      Richar(klonv)                 ! Richardson Number
      real      fac_Ri,vuzvun,Kz_vun

! OUTPUT of Snow Erosion Turbulence
! #b1 real      W_pLMO                        ! Pseudo Obukhov Length  (WRITE)
! #b1 real      W_psim                        ! Pseudo psim(z)         (WRITE)

! OUTPUT of Snow Erosion Turbulence (2)
! #b2 real      W_NUs1                        ! Contrib to U* numerat.1(WRITE)
! #b2 real      W_NUs2                        ! Contrib to U* numerat.2(WRITE)
! #b2 real      W_NUs3                        ! Contrib to U* numerat.3(WRITE)
! #b2 real      W_DUs1                        ! Contrib to U* denomin.1(WRITE)
! #b2 real      W_DUs2                        ! Contrib to U* denomin.2(WRITE)


! Internal DATA
! =============

      data      Theta0/288.0/   ! Potential Reference Temperature
c #ZX data      zetMAX/ 1.e6/   ! Strong Stability Limit
c #zx 
c2306-280611      data      zetMAX/ 1.e1/   ! Strong Stability Lim   
c2806-290611      data      zetMAX/ 1.e0/   ! Strong Stability Limit
                                !(Mahalov et al. 2004, GRL  31 2004GL021055)
chj290911      data      zetMAX/ 4.28/   ! Strong Stability Limit
c                                !(King    et al. 1996, JGR 101(7) p.19121)
      data      zetMAX/ 0.0428/   ! Strong Stability Limit
      data      coef_m/20.  /   ! Stabil.Funct.for Moment.: unstab.coef.
      data      coef_h/15.  /   ! Stabil.Funct.for Heat:    unstab.coef.
c #AE data      SblPom/ 1.27/   ! Lower Boundary Height Parameter 
                                ! for Suspension
                                ! Pommeroy, Gray and Landine 1993, 
                                ! J. Hydrology, 144(8) p.169
c #AE data      nit   / 5   /   ! us(is0,uth) recursivity: Nb Iterations 


! Effective SBL variables
! =======================

      DO ikl=1,knonv
        VVaSBL(ikl)   = VV__SV(ikl)
c #AW !hj060511
        VVaSBL(ikl)   = VVmmem(ikl)
        dTa_Ts(ikl)   = TaT_SV(ikl)-Tsrfsv(ikl)
c #AH!hj060511
        dTa_Ts(ikl)   = dTmmem(ikl)
      ENDDO


! Convergence Criterion
! =====================

      icount = 0

 1    CONTINUE
      icount = icount + 1
      dustar = 0.

        DO ikl=1,knonv

          u0star      = us__SV(ikl)

! u*, u*T*, u*s* Time Moving Averages
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #AM     u0star      = 0.0
c #AT     uTstar      = 0.0
c #AS     usstar      = 0.0
c #AM     DO nt=1,ntaver
c #AM     u0star      = u0star + u__mem(ikl,nt)
c #AT     uTstar      = uTstar + uT_mem(ikl,nt)
c #AS     usstar      = usstar + us_mem(ikl,nt)
c #AM     ENDDO
c #AM     u0star      = u0star / ntaver
c #AM     us__SV(ikl) = u0star
c #AT     uts_SV(ikl) = uTstar / ntaver
c #AS     uss_SV(ikl) = usstar / ntaver


! Turbulent Scales from previous Time Step
! ----------------------------------------

          u0star      =      max(eps6,u0star)      ! Friction Velocity     u*
          uustar      = u0star      * u0star       ! Friction Velocity^2  uu*
          thstar      = uts_SV(ikl) / u0star       ! Temperature       theta*
          qqstar      = uqs_SV(ikl) / u0star       ! Specific Humidity    qq*
          ssstar      = uss_SV(ikl) / u0star       ! Blown    Snow        ss*


! Monin-Obukhov Stability Parameter for Momentum
! ----------------------------------------------

! Pseudo Virtual Temperature Turbulent Scale thetav*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          thstarv     = thstar      + Theta0      *(0.608*qqstar)
     .                             /(1.+0.608*QaT_SV(ikl)-qsnoSV(ikl))
          thstars     =     sign(un_1,thstarv)
          thstara     =      abs(     thstarv)
          thstarv     =      max(eps6,thstara)*thstars

! Pseudo Obukhov Length Scale        (Gall?e et al., 2001 BLM 99, (A2) p.17)
! Full   Obukhov Length Scale        (when Blowing * is ##NOT## switched ON)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          LMO_SV(ikl) = Theta0      * max(eps6,uustar)
     .                /(vonKrm      * Grav_F  *thstarv)

! OUTPUT of Snow Erosion Turbulence
! #b1     W_pLMO      = LMO_SV(ikl)

          zetah       = za__SV(ikl) / LMO_SV(ikl)
          zetam       =           min(zetMAX,zetah)! Strong Stability Limit
                                                   !(Mahalov et al. 2004
                                                   ! GRL 31 2004GL021055)
          LMOmom(ikl) = za__SV(ikl) /(max(eps6,abs(zetam))
     .                              *sign(un_1,    zetam ))
          zeta0m      = Z0m_SV(ikl) / LMOmom(ikl)
          zeta0h      = Z0h_SV(ikl) / LMO_SV(ikl)

! Momentum Pseudo Stability Function (Gall?e et al. 2001, BLM 99, (11) p. 7)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          stab_s      =  max(zer0,sign(un_1,zetam))

          psim_s      =  -A_Stab *zetam
          xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zetam)))
          psim_i      =   2. *log(half*(un_1+xpsimi)) 
     .                       +log(half*(un_1+xpsimi*xpsimi))
     .                   -2.*atan(xpsimi)   +half*piNmbr
          psim_z      =    stab_s*psim_s+(1.-stab_s)*psim_i

! OUTPUT of Snow Erosion Turbulence
! #b1     W_psim      =           psim_z

          psim_s      =  -A_Stab *zeta0m
          xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zeta0m)))
          psim_i      =   2. *log(half*(un_1+xpsimi)) 
     .                       +log(half*(un_1+xpsimi*xpsimi))
     .                   -2.*atan(xpsimi)   +half*piNmbr
          psim_0      =    stab_s*psim_s+(1.-stab_s)*psim_i

! Virtual Temperature Turbulent Scale thetav*    (ss* impact included   )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ needed for new ss*)
c #AE     thstarv     = thstar      + Theta0      *(0.608*qqstar
c #AE.                                                   -ssstar
c #AE.                                             )
c #AE.                             /(1.+0.608*QaT_SV(ikl)-qsnoSV(ikl))
c #AE     thstars     =     sign(un_1,thstarv)
c #AE     thstara     =      abs(     thstarv)
c #AE     thstarv     =      max(eps6,thstara)    *thstars

! Full   Obukhov Length Scale        (Gall?e et al. 2001, BLM 99, (A1) p.16)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #AE     LMO_SV(ikl) = Theta0      * us__SV(ikl)* us__SV(ikl)
c #AE.                /(vonKrm      * Grav_F     * thstarv)

c #AE     zetah       = za__SV(ikl) / LMO_SV(ikl)
c #AE     zetam       =           min(zetMAX,zetah)! Strong Stability Limit
                                                   !(Mahalov et al. 2004
                                                   ! GRL 31 2004GL021055)
c #AE     LMOmom(ikl) = za__SV(ikl) /(max(eps6,abs(zetam))
c #AE.                              *sign(un_1,    zetam ))
c #AE     zeta0m      = Z0m_SV(ikl) / LMOmom(ikl)

! Snow Erosion    Stability Function (Gall?e et al. 2001, BLM 99, (11) p. 7)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #AE     stab_s      =  max(zer0,sign(un_1,zetam))

c #AE     psis_s      =  -AsStab *zetam
c #AE     xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zetam)))
c #AE     psim_i      =   2. *log(half*(un_1+xpsimi)) 
c #AE.                       +log(half*(un_1+xpsimi*xpsimi))
c #AE.                   -2.*atan(xpsimi)   +half*piNmbr
c #AE     psis_z      =    stab_s*psis_s+(1.-stab_s)*psim_i

c #AE     psis_s      =  -AsStab *zeta0m
c #AE     xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zeta0m)))
c #AE     psim_i      =   2. *log(half*(un_1+xpsimi)) 
c #AE.                       +log(half*(un_1+xpsimi*xpsimi))
c #AE.                   -2.*atan(xpsimi)   +half*piNmbr
c #AE     psis_0      =    stab_s*psis_s+(1.-stab_s)*psim_i

! Square Roots of the Drag Coefficient for Snow Erosion Turbulent Flux
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #AE     rCDmSV(ikl) = vonKrm/(sqrCm0(ikl)-psim_z+psim_0)

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #ss     IF (ii__SV(ikl).EQ.iwr_SV.AND.jj__SV(ikl).EQ.jwr_SV      .AND.
! #ss.        nn__SV(ikl).EQ.nwr_SV                          )
! #ss.    write(6,6600)  Z0m_SV(ikl) , psim_z
! #ss.                  ,LMO_SV(ikl) , uustar
! #ss.                  ,sqrCm0(ikl) , psim_0
! #ss.                  ,LMOmom(ikl) , thstarv
 6600     format(/,' ** SISVATeSBL *0  '
     .            ,'  Z0m_SV  = ',e12.4,'  psim_z  = ',e12.4
     .            ,'  LMO_SV  = ',e12.4,'  uustar  = ',e12.4
     .          ,/,'                   '
     .            ,'  sqrCm0  = ',e12.4,'  psim_0  = ',e12.4
     .            ,'  LMOmom  = ',e12.4,'  thstarv = ',e12.4)


! Momentum            Turbulent Scale  u*
! ---------------------------------------

! Momentum            Turbulent Scale  u*          in case of NO Blow. Snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VVa_OK      =  max(0.000001,       VVaSBL(ikl))
          sss__N      =  vonKrm      *       VVa_OK
          sss__F      = (sqrCm0(ikl) - psim_z + psim_0)
          usuth0      =  sss__N /sss__F                ! u* if NO Blow. Snow

! Momentum            Turbulent Scale  u*          in case of    Blow. Snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #AE     sss__G      =  0.27417     * Grav_F

! ______________               _____
! Newton-Raphson (! Iteration, BEGIN)
! ~~~~~~~~~~~~~~               ~~~~~
c #AE     DO iit=1,nit
c #AE     sss__K      =  Grav_F      * r_Stab * A_Stab *za__SV(ikl)
c #AE.                                     *rCDmSV(ikl)*rCDmSV(ikl)
c #AE.                           /(1.+0.608*QaT_SV(ikl)-qsnoSV(ikl))
c #AE     us_127      =  exp(    SblPom *log(us__SV(ikl)))
c #AE     us_227      =  us_127         *    us__SV(ikl)
c #AE     us_327      =  us_227         *    us__SV(ikl)
c #AE     us_427      =  us_327         *    us__SV(ikl)
c #AE     us_527      =  us_427         *    us__SV(ikl)

c #AE     us__SV(ikl) =  us__SV(ikl)
c #AE.    - (  us_527     *sss__F     /sss__N
c #AE.      -  us_427
c #AE.      -  us_227     *qsnoSV(ikl)*sss__K
c #AE.      + (us__SV(ikl)*us__SV(ikl)-usthSV(ikl)*usthSV(ikl))/sss__G)
c #AE.     /(  us_427*5.27*sss__F     /sss__N
c #AE.      -  us_327*4.27
c #AE.      -  us_127*2.27*qsnoSV(ikl)*sss__K
c #AE.      +  us__SV(ikl)*2.0                                 /sss__G)

c #AE     us__SV(ikl)= min(us__SV(ikl),usuth0)
c #AE     us__SV(ikl)= max(us__SV(ikl),eps6  )
c #AE     rCDmSV(ikl)=     us__SV(ikl)/VVa_OK
! #aE     sss__F     =     vonKrm     /rCDmSV(ikl)
c #AE     ENDDO
! ______________               ___  
! Newton-Raphson (! Iteration, END  )
! ~~~~~~~~~~~~~~               ~~~  

c #AE     us_127      =  exp(    SblPom *log(us__SV(ikl)))
c #AE     us_227      =  us_127         *    us__SV(ikl)

! Momentum            Turbulent Scale  u*: 0-Limit in case of no Blow. Snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #AE     dusuth      =  us__SV(ikl) - usthSV(ikl)       ! u* - uth*
c #AE     signus      =  max(sign(un_1,dusuth),zer0)     ! 1 <=> u* - uth* > 0
          us__SV(ikl) =                                  !
c #AE.                   us__SV(ikl)  *signus  +         ! u* (_BS)
     .                   usuth0                          ! u* (nBS)
c #AE.                            *(1.-signus)           !


! Blowing Snow        Turbulent Scale ss*
! ---------------------------------------

c #AE     hSalSV(ikl) = 8.436e-2  *exp(SblPom  *log(us__SV(ikl)))
c #AE     qSalSV(ikl) = (us__SV(ikl) * us__SV(ikl)
c #AE.                  -usthSV(ikl) * usthSV(ikl))*signus
c #AE.                / (sss__G      * us_227     )

c #ae     qSalSV(ikl) = (us__SV(ikl) * us__SV(ikl)
c #ae.                  -usthSV(ikl) * usthSV(ikl))
c #ae.                  *signus      * us__SV(ikl) *3.25
c #ae.                 /(hSalSV(ikl) * Grav_F           )

c #AE     ssstar      =  rCDmSV(ikl) *(qsnoSV(ikl) -qSalSV(ikl))
c #AE.                 * r_Stab

c #AE     uss_SV(ikl) =  min(zer0    , us__SV(ikl) *ssstar)
c #BS     uss_SV(ikl) =  max(-0.002  , uss_SV(ikl)        )

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #ss     IF (ii__SV(ikl).EQ.iwr_SV.AND.jj__SV(ikl).EQ.jwr_SV      .AND.
! #ss.        nn__SV(ikl).EQ.nwr_SV                          )      THEN
! #ss     write(6,6000)  daHost     ,     icount     ,
! #ss.                   us__SV(ikl),1.e3*hSalSV(ikl),
! #ss.              1.e3*Z0m_SV(ikl),
! #ss.              1.e3*qsnoSV(ikl),1.e3*qSalSV(ikl)
! #ss.                  ,usthSV(ikl),     us__SV(ikl)-usthSV(ikl),
! #ss.              1.e3*ssstar     ,1.e3*us__SV(ikl)*ssstar
 6000     format(a18,i3,6x,'u*   [m/s] =',f6.3,'   hSalt[mm]='  ,e9.3,
     .                  '   Z0m   [mm] =',f9.3,'   q   [g/kg] =',f9.3,
     .               /,91x,                    '   qSa [g/kg] =',f9.3,
     .               /,27x, 'ut*[m/s]='  ,e9.3,'   u*-ut*   ='  ,e9.3,
     .                  '   s*  [g/kg] =',f9.3,'   us* [mm/s] =',f9.3)
! #ss     END IF


! Virtual Temperature Turbulent Scale thetav*    (ss* impact included)
! --------------------------------------------------------------------

c #AE     thstarv     = thstar      + Theta0      *(0.608*qqstar
c #AE.                                                   -ssstar
c #AE.                                             )
c #AE.                             /(1.+0.608*QaT_SV(ikl)-qsnoSV(ikl))
c #AE     thstars     =     sign(un_1,thstarv)
c #AE     thstara     =      abs(     thstarv)
c #AE     thstarv     =      max(eps6,thstara)    *thstars


! Full   Obukhov Length Scale (Gall?e et al., 2001, BLM 99, (A1) p.16)
! --------------------------------------------------------------------

c #AE     LMO_SV(ikl) = Theta0      * us__SV(ikl)* us__SV(ikl) 
c #AE.                /(vonKrm      * Grav_F     * thstarv)

c #AE     zetah       = za__SV(ikl) / LMO_SV(ikl)
c #AE     zetam       =           min(zetMAX,zetah)! Strong Stability Limit
                                                   !(Mahalov et al. 2004
                                                   ! GRL 31 2004GL021055)
c #AE     LMOmom(ikl) = za__SV(ikl) /(max(eps6,abs(zetam))
c #AE.                              *sign(un_1,    zetam ))
c #AE     zeta0m      = Z0m_SV(ikl) / LMOmom(ikl)
c #AE     zeta0h      = Z0h_SV(ikl) / LMO_SV(ikl)

! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #ss     IF (ii__SV(ikl).EQ.iwr_SV.AND.jj__SV(ikl).EQ.jwr_SV      .AND.
! #ss.        nn__SV(ikl).EQ.nwr_SV                          )      THEN
! #ss     write(6,6001)  LMO_SV(ikl)    ,    zetah
 6001     format(18x,9x,'LMO  [m]=',f9.1,'   zetah[-] =',f9.3)
! #ss     END IF


! Turbulent Scales
! ----------------

! Momentum Stability Function (Gall?e et al., 2001, BLM 99, (11) p. 7)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          stab_s      =  max(zer0,sign(un_1,zetam))

          psim_s      =  -A_Stab *zetam
          xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zetam)))
          psim_i      =   2. *log(half*(un_1+xpsimi)) 
     .                       +log(half*(un_1+xpsimi*xpsimi))
     .                   -2.*atan(xpsimi)   +half*piNmbr
          psim_z      =    stab_s*psim_s+(1.-stab_s)*psim_i

          psim_s      =  -A_Stab *zeta0m
          xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zeta0m)))
          psim_i      =   2. *log(half*(un_1+xpsimi)) 
     .                       +log(half*(un_1+xpsimi*xpsimi))
     .                   -2.*atan(xpsimi)   +half*piNmbr
          psim_0      =    stab_s*psim_s+(1.-stab_s)*psim_i

! Heat     Stability Function (Gall?e et al., 2001, BLM 99, (11) p. 7)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          stab_s      =  max(zer0,sign(un_1,zetah))

          psih_s      =  -AhStab *zetah
          xpsihi      = sqrt(sqrt(un_1-coef_h*min(zer0,zetah)))
          psih_i      =   2. *log(half*(un_1+xpsihi)) 
          psih_z      =    stab_s*psih_s+(1.-stab_s)*psih_i

          psih_s      =  -AhStab *zeta0h
          xpsihi      = sqrt(sqrt(un_1-coef_h*min(zer0,zeta0h)))
          psih_i      =   2. *log(half*(un_1+xpsihi)) 
          psih_0      =    stab_s*psih_s+(1.-stab_s)*psih_i

! Square Roots of the Drag Coefficients
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          rCDhSV(ikl) = vonKrm*(ExnrSV(ikl)/p0_kap)
     .                        /(sqrCh0(ikl)-psih_z+psih_0)
          rCDmSV(ikl) = vonKrm/(sqrCm0(ikl)-psim_z+psim_0)

! Drag Coefficients
! ~~~~~~~~~~~~~~~~~
          CDh(ikl)    = rCDmSV(ikl) * rCDhSV(ikl)
          CDm(ikl)    = rCDmSV(ikl) * rCDmSV(ikl)

! Real    Temperature Turbulent Scale theta*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          thstar      = rCDhSV(ikl) * dTa_Ts(ikl)
     .                              *(p0_kap     /ExnrSV(ikl))
          uts_SV(ikl) = us__SV(ikl) * thstar


! Convergence Criterion
! =====================

          dustar      = max(dustar,abs(us__SV(ikl)-u0star))

! u*, u*T*, u*s* Time Moving Averages
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #AM     DO nt=   1,ntaver-1
c #AM     u__mem(ikl,nt    ) = u__mem(ikl,nt+1)
c #AT     uT_mem(ikl,nt    ) = uT_mem(ikl,nt+1)
c #AS     us_mem(ikl,nt    ) = us_mem(ikl,nt+1)
c #AM     ENDDO
c #AM     u__mem(ikl,ntaver) = us__SV(ikl)
c #AT     uT_mem(ikl,ntaver) = uts_SV(ikl)
c #AS     us_mem(ikl,ntaver) = uss_SV(ikl)

! OUTPUT of Snow Erosion Turbulence
! #b1    IF (icount     .EQ.1  )                                    THEN
! #b1     write(6,6004)
 6004     format(122('-'))
! #b1    IF (mod(VVaSBL(ikl),4.).LT.0.1)                            THEN
! #b1     write(6,6003)
 6003     format('   V  Ta-Ts  Z0      It'
     .   ,' du*     u*    sss__F   CD       Qss       Qs*     '
     .   ,' PseudOL Full-OL zetam   zetah   psim_z  psih_z')
! #b1     write(6,6004)
! #b1    END IF
! #b1    END IF
! #b1     write(6,6002) VVaSBL(ikl),dTa_Ts(ikl),Z0m_SV(ikl),icount
! #b1.                 ,dustar     ,us__SV(ikl),sss__F
! #b1.                 ,   CDm(ikl),qSalSV(ikl),ssstar
! #b1.                 ,W_pLMO     ,LMO_SV(ikl)
! #b1.                 ,zetam      ,zetah      ,W_psim     ,psih_z
 6002     format(2f6.1,f8.4,i3,f9.6,f6.3,f9.3,3f9.6,2f8.2,2f8.4,2f8.2)

! OUTPUT of Snow Erosion Turbulence (2): u*_AE
! #b2    IF (icount     .EQ.1  )                                    THEN
! #b2     write(6,6014)
 6014     format(100('-'))
! #b2    IF (mod(VVaSBL(ikl),4.).LT.0.1)                            THEN
! #b2     write(6,6013)
 6013     format('   V  Ta-Ts  Z0      It'
     .   ,' du*     u*    sss__F   W_NUs1   W_NUs2   W_NUs3      '
     .   ,' W_DUs1     W_DUs2 ')
! #b2     write(6,6014)
! #b2    END IF
! #b2    END IF
! #b2     write(6,6012) VVaSBL(ikl),dTa_Ts(ikl),Z0m_SV(ikl),icount
! #b2.                 ,dustar     ,us__SV(ikl),sss__F
! #b2.                 ,W_NUs1     ,W_NUs2     ,W_NUs3     
! #b2.                 ,W_DUs1     ,W_DUs2     
 6012     format(2f6.1,f8.4,i3,f9.6,f6.3,f9.3,3f9.3,2f12.3)

        END DO
c #IX IF (                     icount.lt. 3)                     GO TO 1
! hjp if parallel mode, use   IF (                     icount.lt. 3)
!     IF (dustar.gt.0.0001.AND.icount.lt. 6)                     GO TO 1


c #AM   DO ikl=1,knonv
c #AM     u0star      =      0.0
c #AT     uTstar      =      0.0
c #AS     usstar      =      0.0
c #AM     DO nt=1,ntaver
c #AM     u0star      = u0star + u__mem(ikl,nt)
c #AT     uTstar      = uTstar + uT_mem(ikl,nt)
c #AS     usstar      = usstar + us_mem(ikl,nt)
c #AM     ENDDO
c #AM     us__SV(ikl) = u0star / ntaver
c #AT     uts_SV(ikl) = uTstar / ntaver
c #AS     uss_SV(ikl) = usstar / ntaver
c #AM   END DO


! Aerodynamic Resistances
! -----------------------

        DO ikl=1,knonv
          ram_sv(ikl) = 1./(CDm(ikl)*max(VVaSBL(ikl),eps6))
          rah_sv(ikl) = 1./(CDh(ikl)*max(VVaSBL(ikl),eps6))
        END DO


      return
      end




      subroutine SISVAT_SBL
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_SBL                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_SBL generates Surface Boundary Layers Properties   |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns                          |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   za__SV   : Surface Boundary Layer (SBL) Height          [m] |
!     ^^^^^    VV__SV   :(SBL Top)   Wind Velocity                   [m/s] |
!              TaT_SV   : SBL Top    Temperature                       [K] |
!              ExnrSV   : Exner      Potential                         [-] |
!              uqs_SV   : Specific   Humidity  Turbulent Flux        [m/s] |
!              Z0m_SV   : Momentum   Roughness Length                  [m] |
!              Z0h_SV   : Heat       Roughness Length                  [m] |
!              Tsrfsv   : Surface    Temperature                       [K] |
!              sqrCm0   : Contribution of Z0m to Neutral Drag Coefficient  |
!              sqrCh0   : Contribution of Z0h to Neutral Drag Coefficient  |
!                                                                          |
!     INPUT /  LMO_SV   : Monin-Obukhov       Scale                    [m] |
!     OUTPUT:  us__SV   : Friction  Velocity                         [m/s] |
!     ^^^^^^   uts_SV   : Temperature         Turbulent Flux       [K.m/s] |
!                                                                          |
!     OUTPUT:  Fh__sv   : Stability Function                           [-] |
!     ^^^^^^   dFh_sv   : Stability Function (Derivative)              [-] |
!              ram_sv   : Aerodynamic Resistance for Momentum        [s/m] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!                                                                          |
!     WARNING: SISVAT_SBL blows up for too small z0m values & large z_SBL  |
!     ^^^^^^^                      (z0m = 1.8e-6 m for z_SBL = 20 m)       |
!                                                                          |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # stdout               | #sb: OUTPUT/Verification: SISVAT_SBL          |
!                          |      unit  6, SubRoutine  SISVAT_SBL **ONLY** |
!--------------------------------------------------------------------------+
 
 
       
 
 
! Global Variables
! ================
 
      USE  PHY_SV

      USE VAR_SV
      USE VARdSV
 
      USE VARxSV
      USE VARySV

      IMPLICIT NONE 
 
! Internal Variables
! ==================
 
      integer   ikl    ,ist    ,ist__s ,ist__w
      real      CD_m_0 ,CD_h_0 ,ram0   ,rah0  ,rahMIN
      real      d_TaTs ,RiB__D ,RiBulk
      real      bmstab ,Am1_FU ,Am2_FU ,Fm_Uns
      real      bhstab ,Ah1_FU ,Ah2_FU ,Fh_Uns,dFh_Un
      real      Aux_FS ,FStabl ,dFSdRi ,Stabil,Fm_loc
      real      uustar ,thstar ,qqstar ,ssstar
      real      thstarv,thstars,thstara
      real      zeta   ,zeta_S ,zeta_A
      real      zetMAX
 
 
! Internal DATA
! =============
 
      data      zetMAX/ 4.28/   ! Strong Stability Limit
!                               !(King et al. 1996, JGR 101(7) p.19121)
 
! Aerodynamic Resistances
! =======================
 
      DO ikl=1,knonv
 
! Surface Type
! ~~~~~~~~~~~~
        ist    =      isotSV(ikl)                    ! Soil Type
        ist__s =  min(ist, 1)                        ! 1 => Soil
        ist__w =  1 - ist__s                         ! 1 => Water Body
 
! Neutral Parameters
! ~~~~~~~~~~~~~~~~~~
        CD_m_0 =  0.16/    (sqrCm0(ikl)*sqrCm0(ikl)) ! Neutral Drag Coeff.Mom.
        CD_h_0 =  0.16/    (sqrCm0(ikl)*sqrCh0(ikl)) ! Neutral Drag Coeff.Heat
        ram0   =  1.0 /    (CD_m_0     *VV__SV(ikl)) ! Neutral Aero Resis.Mom.
        rah0   =  1.0 /    (CD_h_0     *VV__SV(ikl)) ! Neutral Aero Resis.Heat
 
! Bulk Richardson Number
! ~~~~~~~~~~~~~~~~~~~~~~
        RiB__D =          VV__SV(ikl) *VV__SV(ikl)
     .                   *TaT_SV(ikl)
        d_TaTs =         (TaT_SV(ikl)- Tsrfsv(ikl))
     .                   *p0_kap     / ExnrSV(ikl)
        RiBulk =  Grav_F *za__SV(ikl)* d_TaTs
     .          / RiB__D
 
! OUTPUT/Verification: SISVAT_SBL
! #sb       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #sb.          nn__SV(ikl).GE.nwr_SV)
! #sb.      write(6,6600) Tsrfsv(ikl),TaT_SV(ikl),VV__SV(ikl)
! #sb.                  , d_TaTs     ,RiBulk
 6600       format(/,'Tem(s,a), Wind  , d_TaTs, RiBulk = ',5e15.6)
 
! Unstable Case
! ~~~~~~~~~~~~~
        bmstab =  ist__s * (13.7 -0.34 /sqrt(CD_m_0))! Momentum
     .          + ist__w *   4.9                     !
        bmstab =  10.    *  bmstab   * CD_m_0        !
     .              *sqrt(za__SV(ikl)/ Z0m_SV(ikl))  !
        Am1_FU =  bmstab *    sqrt(abs(RiBulk))      !
        Am2_FU =  Am1_FU +1.0 +10.*abs(RiBulk)       !
        Fm_Uns = (Am1_FU +1.0)/        Am2_FU        !
 
! OUTPUT/Verification: SISVAT_SBL
! #sb       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #sb.          nn__SV(ikl).GE.nwr_SV)
! #sb.      write(6,6601) CD_m_0     ,Z0m_SV(ikl),bmstab     
! #sb.                  , ist__s     ,ist__w
 6601       format(/,'CD_m_0  , Z0m_SV, bmstab, ist/sw = ',3e15.6,2i15)
 
        bhstab =  ist__s * ( 6.3 -0.18 /sqrt(CD_h_0))! Heat
     .          + ist__w *   2.6                     !
        bhstab =  10.    *  bhstab   * CD_h_0        !
     .              *sqrt(za__SV(ikl)/ Z0h_SV(ikl))  !
        Ah1_FU =  bhstab *    sqrt(abs(RiBulk))      !
        Ah2_FU =  Ah1_FU +1.0 +10.*abs(RiBulk)       !
        Fh_Uns = (Ah1_FU +1.0)/        Ah2_FU        !
        dFh_Un =((Ah1_FU +2.0)/(Ah2_FU*Ah2_FU)) * 5. !
 
! Stable   Case
! ~~~~~~~~~~~~~
        Aux_FS =          1.0 + 5.*    RiBulk
        FStabl =                Aux_FS*Aux_FS
        dFSdRi =                Aux_FS          *10.
 
! Effective Stability Functions and Derivatives
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Stabil      =          sign(un_1,d_TaTs)
        Fm_loc      =  FStabl * max(zer0,Stabil)
     .               - Fm_Uns * min(zer0,Stabil)
        Fh__sv(ikl) =  FStabl * max(zer0,Stabil)
     .               - Fh_Uns * min(zer0,Stabil)
        dFh_sv(ikl) =  dFSdRi * max(zer0,Stabil)
     .               - dFh_Un * min(zer0,Stabil)
 
! OUTPUT/Verification: SISVAT_SBL
! #sb       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #sb.          nn__SV(ikl).GE.nwr_SV)
! #sb.      write(6,6602) FStabl     ,Stabil     
! #sb.                   ,Fm_Uns     ,Fm_loc     
 6602       format(/,'FStabl  , Stabil, Fm_Uns, Fm_loc = ',4e15.6)
 
! Aerodynamic Resistances
! ~~~~~~~~~~~~~~~~~~~~~~~
        ram_sv(ikl) =  ram0   *          Fm_loc
        rah_sv(ikl) =  rah0   *          Fh__sv(ikl)
        rahMIN   = max(rah_sv(ikl),  abs(d_TaTs)*60./za__SV(ikl))
                                               ! 60 for 30dgC within 1/2 hour
        dFh_sv(ikl) =  rah0   *          dFh_sv(ikl)
     .              *  rahMIN /          rah_sv(ikl)
        rah_sv(ikl) =  rahMIN
 
 
! Square Root Contributions to the Drag Coefficients
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rCDmSV(ikl) =  sqrt(ram_sv(ikl) *VV__SV(ikl))
        rCDmSV(ikl) =  1.     / max(eps6,rCDmSV(ikl))
        rCDhSV(ikl) =       rah_sv(ikl) *VV__SV(ikl)
     .                                  *rCDmSV(ikl)
        rCDhSV(ikl) = (1.     / max(eps6,rCDhSV(ikl)))
     .              * (     ExnrSV(ikl) /p0_kap      )
 
! OUTPUT/Verification: SISVAT_SBL
! #sb       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #sb.          nn__SV(ikl).GE.nwr_SV)
! #sb.      write(6,6603) ram_sv(ikl),rah_sv(ikl)
! #sb.                   ,rCDmSV(ikl),rCDhSV(ikl)
 6603       format(/,'AeR(m,h), rCD(m,h)               = ',4e15.6)
 
 
! Turbulent Scales
! ================
 
! Friction Velocity                   u*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        uustar      =      VV__SV(ikl) / ram_sv(ikl)
        us__SV(ikl) = sqrt(uustar)
 
! Real    Temperature Turbulent Scale theta*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        uts_SV(ikl) =      d_TaTs      / rah_sv(ikl)
        thstar      =      uts_SV(ikl) / us__SV(ikl)
 
! Specific Humidity   Turbulent Scale qq*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        qqstar      =      uqs_SV(ikl) / us__SV(ikl)
 
! Virtual Temperature Turbulent Scale thetav*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        thstarv     = thstar      + TaT_SV(ikl) *(0.608*qqstar
     .                                             )
        thstars     =     sign(un_1,thstarv)
        thstara     =      abs(     thstarv)
        thstarv     =      max(eps6,thstara)    *thstars
 
! Monin Obukhov Scale Height
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        LMO_SV(ikl) = TaT_SV(ikl) * uustar
     .              /(vonKrm      * Grav_F     * thstarv)
        zeta        = za__SV(ikl) / LMO_SV(ikl)
        zeta        =           min(zetMAX,zeta)  ! Strong Stability Limit
!                                                 ! King et al.   1996
!                                                 ! JGR 101(7) p.19121
        zeta_S      =          sign(un_1  ,zeta)
        zeta_A      =           abs(       zeta)
        zeta        = zeta_S  * max(eps6  ,zeta_A)
        LMO_SV(ikl) = za__SV(ikl) / zeta
 
! OUTPUT/Verification: SISVAT_SBL
! #sb       IF (ii__SV(ikl).EQ.iwr_SV .AND. jj__SV(ikl).EQ. jwr_SV .AND.
! #sb.          nn__SV(ikl).GE.nwr_SV)
! #sb.      write(6,6604) us__SV(ikl),uts_SV(ikl)
! #sb.                   ,LMO_SV(ikl),zeta        
 6604       format(/,'***(m,h), LMO   , zeta           = ',4e15.6)
 
      END DO
 
      return
      end


      subroutine SISVAT_TVg
! #e1.                     (ETVg_d)

!--------------------------------------------------------------------------+
!   MAR          SISVAT_TVg                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_TVg computes the Canopy Energy Balance             |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   ivgtSV   = 0,...,12:   Vegetation Type                      |
!     ^^^^^               0:          Water, Solid or Liquid               |
!              isnoSV   = total Nb of Ice/Snow Layers                      |
!                                                                          |
!     INPUT:   sol_SV   : Downward Solar Radiation                  [W/m2] |
!     ^^^^^    IRd_SV   : Surface  Downward Longwave Radiation      [W/m2] |
!              TaT_SV   : SBL Top  Temperature                         [K] |
!              rhT_SV   : SBL Top  Air  Density                    [kg/m3] |
!              QaT_SV   : SBL Top  Specific  Humidity              [kg/kg] |
!              psivSV   : Leaf     Water     Potential                 [m] |
!              IRs_SV   : Soil     IR Flux  (previous time step)    [W/m2] |
!              dt__SV   : Time     Step                                [s] |
!                                                                          |
!              SoCasv   : Absorbed Solar Radiation by Canopy (Normaliz)[-] |
!              tau_sv   : Fraction of Radiation transmitted by Canopy  [-] |
!              Evg_sv   : Soil+Vegetation Emissivity                   [-] |
!              Eso_sv   : Soil+Snow       Emissivity                   [-] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!              Sigmsv   : Canopy Ventilation Factor                    [-] |
!              LAI_sv   : Leaf Area  Index                             [-] |
!              LAIesv   : Leaf Area  Index (effective / transpiration) [-] |
!              glf_sv   : Green Leaf Fraction of NOT fallen Leaves     [-] |
!              rrMxsv   : Canopy Maximum Intercepted Rain          [kg/m2] |
!                                                                          |
!     INPUT /  TvegSV   : Canopy   Temperature                         [K] |
!     OUTPUT:  rrCaSV   : Canopy     Water     Content             [kg/m2] |
!     ^^^^^^                                                               |
!                                                                          |
!     OUTPUT:  IRv_sv   : Vegetation IR Flux                        [W/m2] |
!     ^^^^^^   HSv_sv   : Sensible Heat Flux                        [W/m2] |
!              HLv_sv   : Latent   Heat Flux                        [W/m2] |
!              Evp_sv   : Evaporation                              [kg/m2] |
!              EvT_sv   : Evapotranspiration                       [kg/m2] |
!              ETVg_d   : Vegetation  Energy Power Forcing          [W/m2] |
!                                                                          |
!     Internal Variables:                                                  |
!     ^^^^^^^^^^^^^^^^^^                                                   |
!                                                                          |
!     METHOD: The Newton-Raphson Scheme is preferable                      |
!     ^^^^^^  when computing over a long time step the heat content        |
!             of a medium having a very small or zero heat capacity.       |
!             This is to handle strong non linearities arising             |
!             in conjunction with rapid temperature variations.            |
!                                                                          |
!     REFERENCE: DR97: Koen de Ridder thesis, UCL, 1997                    |
!     ^^^^^^^^^                                                            |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #NN: Newton-Raphson Increment not added in last Iteration            |
!     #NC: OUTPUT Preparation for Stand Alone NetCDF File                  |
!                                                                          |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV
      USE VARdSV

      USE VARxSV
      USE VARySV


      IMPLICIT NONE

! OUTPUT
! ------

! OUTPUT/Verification: Energy/Water Budget
! #e1 real      ETVg_d(klonv)                 ! VegetationPower, Forcing

! OUTPUT for Stand Alone NetCDF File
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #NC real          SOsoKL(klonv)             ! Absorbed Solar Radiation
c #NC real          IRsoKL(klonv)             ! Absorbed IR    Radiation
c #NC real          HSsoKL(klonv)             ! Absorbed Sensible Heat Flux
c #NC real          HLsoKL(klonv)             ! Absorbed Latent   Heat Flux
c #NC real          HLs_KL(klonv)             ! Evaporation
c #NC real          HLv_KL(klonv)             ! Transpiration
c #NC common/DumpNC/SOsoKL,IRsoKL
c #NC.             ,HSsoKL,HLsoKL
c #NC.             ,HLs_KL,HLv_KL


! Internal Variables
! ==================

      integer   ikl                           ! Grid Point Index
      integer   nitmax,nit                    ! Iterations Counter
      real      d_Tveg                        ! Canopy Temperat. Increment
      real      dTvMAX                        ! Canopy Temperat. Increment MAX
      real      dHvdTv                        ! Derivativ.of Canopy Energ.Budg.
      real      Hv_Tv0                        ! Imbalance of Canopy Energ.Budg.
      real      Hv_MAX                        ! MAX Imbal.of Canopy Energ.Budg.
      real      Hv_MIN                        ! MIN Imbal.of Canopy Energ.Budg.
      real      Hswich                        ! Newton-Raphson         Switch
      real      Tveg_0(klonv)                 ! Canopy Temperature, Previous t
      real      tau_Ca                        ! Canopy IR Radiation Absorption
      real      IR_net                        ! InfraRed  NET(t)
      real      dIRdTv(klonv)                 ! InfraRed  NET(t), Derivative(t)
      real      dHSdTv(klonv)                 ! Sensible Heat FL. Derivative(t)
      real      dHLdTv(klonv)                 ! Latent   Heat FL. Derivative(t)
      real      EvFrac                        ! Condensat./Transpirat. Switch
      real      SnoMsk                        ! Canopy Snow            Switch
      real      den_qs,arg_qs,qsatvg          ! Canopy Saturat. Spec. Humidity
      real      dqs_dT                        ! d(qsatvg)/dTv
      real      FacEvp,FacEvT,Fac_Ev          ! Evapo(transpi)ration Factor
      real      dEvpdT(klonv),dEvTdT(klonv)   ! Evapo(transpi)ration Derivative
      real      F_Stom                        ! Funct.  (Leaf Water Potential)
      real      R0Stom                        ! Minimum Stomatal Resistance
      real      R_Stom                        !         Stomatal Resistance
      real      LAI_OK                        ! 1. ==>  Leaves   exist
      real      rrCaOK,snCaOK,dEvpOK          ! Positive Definiteness Correct.


! Internal DATA
! =============

      data      nitmax /   5   /              ! Maximum  Iterations    Number
      data      dTvMAX /   5.  /              ! Canopy Temperat. Increment MAX
      data      Hv_MIN /   0.1 /              ! MIN Imbal. of Surf.Energy Budg.
      data      SnoMsk /   0.0 /              ! Canopy Snow Switch    (Default)


! Newton-Raphson Scheme
! =====================

      nit    = 0
  101 CONTINUE
      nit    = nit + 1
      HV_MAX = 0.


! Temperature of the Previous Time Step
! -------------------------------------

        DO ikl=1,knonv
          Tveg_0(ikl) = TvegSV(ikl)


! IR    Radiation Absorption
! --------------------------

          tau_Ca = 1.  - tau_sv(ikl)                        ! Canopy Absorption
          IRv_sv(ikl) = -2.0   *Evg_sv(ikl)  *StefBo        !
     .                         *TvegSV(ikl)  *TvegSV(ikl)   ! Downward IR (OUT)
     .                         *TvegSV(ikl)  *TvegSV(ikl)   ! + Upward IR (OUT)
          dIRdTv(ikl) =
     .     -Evg_sv(ikl)*                                    !
     .                8.*StefBo*TvegSV(ikl)  *TvegSV(ikl)   ! Downward IR (OUT)
     .                         *TvegSV(ikl)                 ! + Upward IR (OUT)
          IR_net =       tau_Ca                             !
     .    *(Evg_sv(ikl)* IRd_SV(ikl)                        ! Downward IR (IN)
     .     -             IRs_SV(ikl)                        !   Upward IR (IN)
     .     +             IRv_sv(ikl))                       !          IR (OUT)


! Sensible Heat Flux
! ------------------

          dHSdTv(ikl) = rhT_SV(ikl)* Sigmsv(ikl) *CpdAir    ! Derivative, t(n)
     .                / rah_sv(ikl)                         ! 
          HSv_sv(ikl) = dHSdTv(ikl)                         ! Value,      t(n)
     .                *(TaT_SV(ikl)-TvegSV(ikl))            !


! Latent   Heat Flux
! ------------------

! Canopy Saturation Specific Humidity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          den_qs =         TvegSV(ikl)- 35.8
          arg_qs = 17.27 *(TvegSV(ikl)-273.16)/den_qs
          qsatvg = .0038 * exp(arg_qs)
          dqs_dT = qsatvg* 4099.2    /(den_qs *den_qs)

! Canopy Stomatal Resistance
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          R0Stom = min(         StodSV(ivgtSV(ikl))
     .                /max(eps6,glf_sv(       ikl)),StxdSV)  ! Min Stomatal R.
          F_Stom = pscdSV / max(pscdSV-psivSV(ikl) ,eps6)    ! F(Leaf Wat.Pot.)
                                                             ! DR97, eqn. 3.22
          R_Stom =(R0Stom / max(LAIesv(ikl), R0Stom/StxdSV)) ! Can.Stomatal R.
     .           * F_Stom                                    ! DR97, eqn. 3.21

! Evaporation / Evapotranspiration
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          SnoMsk = max(zer0, sign(un_1,snCaSV(ikl)-eps_21))  !
          EvFrac = max(zer0, sign(un_1,QaT_SV(ikl)-qsatvg))  ! Condensation/
          EvFrac = EvFrac                                    ! Transpiration
     .       + (1.-EvFrac)*((1-SnoMsk)*         rrCaSV(ikl)  !        Switch
     .                                         /rrMxsv(ikl)  !
     .                      +  SnoMsk *min(un_1,snCaSV(ikl)  !
     .                                         /rrMxsv(ikl)))!
          Fac_Ev = rhT_SV(ikl) *Sigmsv(ikl)                  ! Idem,  Factor
          FacEvp = Fac_Ev      *EvFrac     / rah_sv(ikl)     !
          Evp_sv(ikl) = FacEvp*(qsatvg     - QaT_SV(ikl))    ! Evaporation
          dEvpdT(ikl) = FacEvp* dqs_dT                       ! Evp Derivative
          FacEvt = Fac_Ev * (1.-EvFrac)    /(rah_sv(ikl)     !
     .                              +R_Stom *Sigmsv(ikl))    !
          EvT_sv(ikl) = FacEvt*(qsatvg     - QaT_SV(ikl))    ! EvapoTranspir.
          dEvTdT(ikl) = FacEvt* dqs_dT                       ! EvT Derivative
          HLv_sv(ikl) =-LhvH2O*(Evp_sv(ikl)+ EvT_sv(ikl))    ! Latent   Heat
     .                 -LhfH2O* Evp_sv(ikl)* SnoMsk          !(Subli.Contrib.)
          dHLdTv(ikl) = LhvH2O*(dEvpdT(ikl)+ dEvTdT(ikl))    !
     .                 +LhfH2O* dEvpdT(ikl)* SnoMsk          !


! Imbalance  of the Canopy  Energy Budget
! ---------------------------------------

          LAI_OK = max(zer0,                                 ! NO Budget if 
     .                 sign(un_1, LAI_sv(ikl)-eps_21))       ! no Leaves
          Hv_Tv0 =    (  SoCasv(ikl)         *sol_SV(ikl)    ! Absorbed Solar
     .                 + IR_net                              ! NET      IR
     .                 + HSv_sv(ikl)                         ! Sensible Heat
     .                 + HLv_sv(ikl)                         ! Latent   Heat
     .                ) *LAI_OK                              !

! OUTPUT/Verification: Energy/Water Budget
! #e1     ETVg_d(ikl) =  Hv_Tv0                              ! Veg.Energ.Bal.

          Hswich      =          1.00
c #NN     Hswich      = max(zer0,                            ! Newton-Raphson
c #NN.                      sign(un_1,    abs(Hv_Tv0     )   !         Switch
c #NN.                                       -Hv_MIN      )) !


! Derivative of the Canopy  Energy Budget
! ---------------------------------------

          dHvdTv    =   dIRdTv(ikl) * max(eps_21,tau_Ca)
     .                - dHSdTv(ikl)
     .                - dHLdTv(ikl)


! Update Canopy and Surface/Canopy Temperatures
! ---------------------------------------------

          d_Tveg      = Hv_Tv0       / dHvdTv               !
          d_Tveg      =      sign(un_1,d_Tveg)              ! Increment
     .                       *min( abs(d_Tveg)     ,dTvMAX) ! Limitor
          TvegSV(ikl) = TvegSV(ikl)  - Hswich      *d_Tveg  ! Newton-Raphson
          Hv_MAX      = max(Hv_MAX,abs(Hv_Tv0     ))        !


! Update Vegetation Fluxes
! ------------------------

c #NN     IRv_sv(ikl) = IRv_sv(ikl)-dIRdTv(ikl)    *d_Tveg  ! Emitted  IR
c #NN     HSv_sv(ikl) = HSv_sv(ikl)+dHSdTv(ikl)    *d_Tveg  ! Sensible Heat
c #NN     Evp_sv(ikl) = Evp_sv(ikl)-dEvpdT(ikl)    *d_Tveg  ! Evapotranspir.
c #NN     EvT_sv(ikl) = EvT_sv(ikl)-dEvTdT(ikl)    *d_Tveg  ! Evapotranspir.
c #NN     HLv_sv(ikl) = HLv_sv(ikl)+dHLdTv(ikl)    *d_Tveg  ! Latent   Heat

          IRv_sv(ikl) = IRv_sv(ikl)                *LAI_OK
          HSv_sv(ikl) = HSv_sv(ikl)                *LAI_OK
          Evp_sv(ikl) = Evp_sv(ikl)                *LAI_OK
          EvT_sv(ikl) = EvT_sv(ikl)                *LAI_OK
          HLv_sv(ikl) = HLv_sv(ikl)                *LAI_OK
        END DO

c #IX IF (                     nit.lt.nitmax)    GO TO 101
c hj150311 for parallel      IF (Hv_MAX.gt.Hv_MIN.and.nit.lt.nitmax)    GO TO 101
      IF (                     nit.lt.nitmax)    GO TO 101                     
        DO ikl=1,knonv
          IRv_sv(ikl) = IRv_sv(ikl)                         ! Emitted  IR
     .   +dIRdTv(ikl) *(TvegSV(ikl)-Tveg_0(ikl))            ! 
          HSv_sv(ikl) = HSv_sv(ikl)                         ! Sensible Heat
     .   -dHSdTv(ikl) *(TvegSV(ikl)-Tveg_0(ikl))            !
          Evp_sv(ikl) = Evp_sv(ikl)                         ! Evaporation
     .   +dEvpdT(ikl) *(TvegSV(ikl)-Tveg_0(ikl))            !
          EvT_sv(ikl) = EvT_sv(ikl)                         ! Transpiration
     .   +dEvTdT(ikl) *(TvegSV(ikl)-Tveg_0(ikl))            !
          HLv_sv(ikl) = HLv_sv(ikl)                         ! Latent   Heat
     .   -dHLdTv(ikl) *(TvegSV(ikl)-Tveg_0(ikl))            !

! OUTPUT for Stand Alone NetCDF File
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #NC     HLv_KL(ikl) = HLv_sv(ikl)


! Update Canopy Water Content
! ---------------------------

          rrCaSV(ikl) = rrCaSV(ikl)-(1.-SnoMsk)*Evp_sv(ikl)*dt__SV
          snCaSV(ikl) = snCaSV(ikl)-    SnoMsk *Evp_sv(ikl)*dt__SV

! Correction for Positive Definiteness (see WKarea/EvpVeg/EvpVeg.f)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          rrCaOK =  max(rrCaSV(ikl), 0.)
          snCaOK =  max(snCaSV(ikl), 0.)
          dEvpOK =     (rrCaOK-rrCaSV(ikl)
     .                 +snCaOK-snCaSV(ikl))/dt__SV

          Evp_sv(ikl) = Evp_sv(ikl)       - dEvpOK          ! Evaporation
          HLv_sv(ikl) = HLv_sv(ikl)                         ! Latent   Heat
     .    +(1.-SnoMsk)* LhvH2O            * dEvpOK          !
     .    +    SnoMsk *(LhvH2O+LhfH2O)    * dEvpOK          !

          rrCaSV(ikl) = rrCaOK
          snCaSV(ikl) = snCaOK

        END DO


      return
      end


      subroutine SISVAT_TSo
! #e1.                     (ETSo_0,ETSo_1,ETSo_d)

!--------------------------------------------------------------------------+
!   MAR          SISVAT_TSo                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_TSo computes the Soil/Snow Energy Balance          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   isotSV   = 0,...,11:   Soil       Type                      |
!     ^^^^^               0:          Water, Solid or Liquid               |
!              isnoSV   = total Nb of Ice/Snow Layers                      |
!              dQa_SV   = Limitation of  Water Vapor  Turbulent Flux       |
!                                                                          |
!     INPUT:   sol_SV   : Downward Solar Radiation                  [W/m2] |
!     ^^^^^    IRd_SV   : Surface Downward  Longwave   Radiation    [W/m2] |
!              za__SV   : SBL Top    Height                            [m] |
!              VV__SV   : SBL Top    Wind Speed                      [m/s] |
!              TaT_SV   : SBL Top    Temperature                       [K] |
!              rhT_SV   : SBL Top    Air  Density                  [kg/m3] |
!              QaT_SV   : SBL Top    Specific  Humidity            [kg/kg] |
!              LSdzsv   : Vertical   Discretization Factor             [-] |
!                       =    1. Soil                                       |
!                       = 1000. Ocean                                      |
!              dzsnSV   : Snow Layer Thickness                         [m] |
!              ro__SV   : Snow/Soil  Volumic Mass                  [kg/m3] |
!              eta_SV   : Soil Water Content                       [m3/m3] |
!              dt__SV   : Time Step                                    [s] |
!                                                                          |
!              SoSosv   : Absorbed Solar Radiation by Surfac.(Normaliz)[-] |
!              IRv_sv   : Vegetation  IR Radiation                  [W/m2] |
!              tau_sv   : Fraction of Radiation transmitted by Canopy  [-] |
!              Evg_sv   : Soil+Vegetation Emissivity                   [-] |
!              Eso_sv   : Soil+Snow       Emissivity                   [-] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!              Lx_H2O   : Latent Heat of Vaporization/Sublimation   [J/kg] |
!              Sigmsv   : Canopy Ventilation Factor                    [-] |
!              sEX_sv   : Verticaly Integrated Extinction Coefficient  [-] |
!                                                                          |
!     INPUT /  TsisSV   : Soil/Ice Temperatures (layers -nsol,-nsol+1,..,0)|
!     OUTPUT:           & Snow     Temperatures (layers  1,2,...,nsno) [K] |
!     ^^^^^^                                                               |
!                                                                          |
!     OUTPUT:  IRs_SV   : Soil      IR Radiation                    [W/m2] |
!     ^^^^^^   HSs_sv   : Sensible  Heat Flux                       [W/m2] |
!              HLs_sv   : Latent    Heat Flux                       [W/m2] |
!              ETSo_0   : Snow/Soil Energy Power, before Forcing    [W/m2] |
!              ETSo_1   : Snow/Soil Energy Power, after  Forcing    [W/m2] |
!              ETSo_d   : Snow/Soil Energy Power         Forcing    [W/m2] |
!                                                                          |
!     METHOD: NO   Skin Surface Temperature                                |
!     ^^^^^^  Semi-Implicit Crank Nicholson Scheme                         |
!                                                                          |
!     REFERENCE: DR97: Koen de Ridder thesis, UCL, 1997                    |
!     ^^^^^^^^^                                                            |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #VX: TURBULENCE: u*q* limited to SBL Saturat.Specif.Humid.           |
!     #RC: TURBULENCE: Richardson Number:    T Derivative is used          |
!     #DL: TURBULENCE: Latent     Heat Flux: T Derivative is used          |
!     #NC: OUTPUT Preparation for Stand Alone NetCDF File                  |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_iii_jjj_n     | #e1: OUTPUT/Verification: Energy Conservation |
!                          |                                               |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV
      USE VARdSV
      USE VAR0SV

      USE VARxSV
      USE VARySV
!hj220711
      USE VARtSV

      IMPLICIT NONE

! OUTPUT/Verification: Energy/Water Budget
! #e1 real      ETSo_0(klonv)                 ! Soil/Snow Power, before Forcing
! #e1 real      ETSo_1(klonv)                 ! Soil/Snow Power, after  Forcing
! #e1 real      ETSo_d(klonv)                 ! Soil/Snow Power, Forcing


! Internal Variables
! ==================

      integer   ikl   ,isl   ,jsl   ,ist      !
      integer   ist__s,ist__w                 ! Soil/Water  Body Identifier
      integer   islsgn                        ! Soil/Snow Surfac.Identifier
      real      eps__3                        ! Arbitrary    Low Number
      real      etaMid,psiMid                 ! Layer Interface's Humidity
      real      mu_eta                        !     Soil thermal Conductivity
      real      mu_exp                        ! arg Soil thermal Conductivity
      real      mu_min                        ! Min Soil thermal Conductivity
      real      mu_max                        ! Max Soil thermal Conductivity
      real      mu_sno(klonv),mu_aux          !     Snow thermal Conductivity
      real      mu__dz(klonv,-nsol:nsno+1)    ! mu_(eta,sno)   / dz
      real      dtC_sv(klonv,-nsol:nsno)      ! dt      / C
      real      IRs__D(klonv)                 ! UpwardIR Previous Iter.Contr.
      real      dIRsdT(klonv)                 ! UpwardIR           T Derivat.
      real      f_HSHL(klonv)                 ! Factor common to HS and HL
      real      dRidTs(klonv)                 ! d(Rib)/d(Ts)
      real      HS___D(klonv)                 ! Sensible Heat Flux Atm.Contr.
      real      f___HL(klonv)                 !
      real      HL___D(klonv)                 ! Latent   Heat Flux Atm.Contr.
      REAL      TSurf0(klonv),dTSurf          ! Previous Surface Temperature
      real      qsatsg(klonv),den_qs,arg_qs   ! Soil   Saturat. Spec. Humidity
      real      dqs_dT(klonv)                 ! d(qsatsg)/dTv
      real      Psi(   klonv)                 ! 1st Soil Layer Water Potential
      real      RHuSol(klonv)                 ! Soil Surface Relative Humidity
      real      etaSol                        ! Soil Surface          Humidity
      real      d__eta                        ! Soil Surface Humidity Increm.
      real      Elem_A,Elem_C                 !   Diagonal Coefficients
      real      Diag_A(klonv,-nsol:nsno)      ! A Diagonal
      real      Diag_B(klonv,-nsol:nsno)      ! B Diagonal
      real      Diag_C(klonv,-nsol:nsno)      ! C Diagonal
      real      Term_D(klonv,-nsol:nsno)      !   Independant Term
      real      Aux__P(klonv,-nsol:nsno)      ! P Auxiliary Variable
      real      Aux__Q(klonv,-nsol:nsno)      ! Q Auxiliary Variable
      real      Ts_Min,Ts_Max                 ! Temperature Limits

! OUTPUT/Verification: Energy/Water Budget
! #e1 real      Exist0                        ! Existing Layer Switch

      integer   nt_srf,it_srf,itEuBk          ! HL: Surface Scheme
      parameter(nt_srf=10)                    !
      real      agpsrf,xgpsrf,dt_srf,dt_ver   !
      real      etaBAK(klonv)                 !
      real      etaNEW(klonv)                 !
      real      etEuBk(klonv)                 !
      real      fac_dt(klonv),faceta(klonv)   !
      real      PsiArg(klonv),SHuSol(klonv)   !

! OUTPUT for Stand Alone NetCDF File
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #NC real          SOsoKL(klonv)             ! Absorbed Solar Radiation
c #NC real          IRsoKL(klonv)             ! Absorbed IR    Radiation
c #NC real          HSsoKL(klonv)             ! Absorbed Sensible Heat Flux
c #NC real          HLsoKL(klonv)             ! Absorbed Latent   Heat Flux
c #NC real          HLs_KL(klonv)             ! Evaporation
c #NC real          HLv_KL(klonv)             ! Transpiration
c #NC common/DumpNC/SOsoKL,IRsoKL
c #NC.             ,HSsoKL,HLsoKL
c #NC.             ,HLs_KL,HLv_KL


! Internal DATA
! =============

      data      eps__3 /   1.e-3   /          ! Arbitrary    Low Number
      data      mu_exp /  -0.4343  /          !     Soil Thermal Conductivity
      data      mu_min /   0.172   /          ! Min Soil Thermal Conductivity
      data      mu_max /   2.000   /          ! Max Soil Thermal Conductivity
      data      Ts_Min / 175.      /          ! Temperature            Minimum
      data      Ts_Max / 300.      /          ! Temperature Acceptable Maximum
                                              ! including   Snow Melt  Energy


! Heat Conduction Coefficient (zero in the Layers over the highest one)
! ===========================
!                               ---------------- isl    eta_SV, rho C (isl)
!   
! Soil                          ++++++++++++++++        etaMid,    mu (isl)
! ----
!                               ---------------- isl-1  eta_SV, rho C (isl-1)
           isl=-nsol
        DO ikl=1,knonv
          mu__dz(ikl,isl) = 0.

          dtC_sv(ikl,isl) = dtz_SV(isl)                   ! dt / (dz X rho C)
     .                   /((rocsSV(isotSV(ikl))           ! [s / (m.J/m3/K)]
     .                     +rcwdSV*eta_SV(ikl,isl))       !
     .                     *LSdzsv(ikl)            )      !
        END DO
      DO   isl=-nsol+1,0
        DO ikl=1,knonv
          ist    =      isotSV(ikl)                       ! Soil Type
          ist__s =  min(ist, 1)                           ! 1 => Soil
          ist__w =  1 - ist__s                            ! 1 => Water Body

          etaMid = 0.5*(dz_dSV(isl-1)*eta_SV(ikl,isl-1)   ! eta at layers
     .                 +dz_dSV(isl)  *eta_SV(ikl,isl)  )  !     interface
     .                 /dzmiSV(isl)                       ! LSdzsv implicit !
          etaMid =  max(etaMid,eps6)
          psiMid =      psidSV(ist)
     .                *(etadSV(ist)/etaMid)**bCHdSV(ist)
          mu_eta =      3.82      *(psiMid)**mu_exp       ! Soil Thermal
          mu_eta =  min(max(mu_eta, mu_min), mu_max)      ! Conductivity
                                                          ! DR97 eq.3.31
          mu_eta =  ist__s *mu_eta +ist__w * vK_dSV       ! Water Bodies
                                                          ! Correction
          mu__dz(ikl,isl) = mu_eta/(dzmiSV(isl)           !
     .                             *LSdzsv(ikl))          !

          dtC_sv(ikl,isl) = dtz_SV(isl)                   ! dt / (dz X rho C)
     .                   /((rocsSV(isotSV(ikl))           !
     .                     +rcwdSV*eta_SV(ikl,isl))       !
     .                     *LSdzsv(ikl)            )      !
        END DO
      END DO


! Soil/Snow Interface
! -------------------

! Soil Contribution
! ^^^^^^^^^^^^^^^^^
           isl=1
        DO ikl=1,knonv
          ist    =      isotSV(ikl)                       ! Soil Type
          ist__s =  min(ist, 1)                           ! 1 => Soil
          ist__w =  1 - ist__s                            ! 1 => Water Body
          psiMid =      psidSV(ist)                       ! Snow => Saturation
          mu_eta =      3.82      *(psiMid)**mu_exp       ! Soil Thermal
          mu_eta =  min(max(mu_eta, mu_min), mu_max)      ! Conductivity
                                                          ! DR97 eq.3.31
          mu_eta =  ist__s *mu_eta +ist__w * vK_dSV       ! Water Bodies

! Snow Contribution
! ^^^^^^^^^^^^^^^^^
          mu_sno(ikl) =  CdidSV                           !
     .                 *(ro__SV(ikl,isl) /rhoWat) ** 1.88 !
          mu_sno(ikl) =          max(eps6,mu_sno(ikl))    !
!         mu_sno :  Snow Heat Conductivity Coefficient [Wm/K]
!                   (Yen 1981, CRREL Rep., 81-10)

! Combined Heat Conductivity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          mu__dz(ikl,isl) = 2./(dzsnSV(ikl,isl  )         ! Combined Heat
     .                         /mu_sno(ikl)               ! Conductivity
     .                         +LSdzsv(ikl)               !
     .                         *dz_dSV(    isl-1)/mu_eta) ! Coefficient

! Inverted Heat Capacity
! ^^^^^^^^^^^^^^^^^^^^^^
          dtC_sv(ikl,isl) = dt__SV/max(eps6,              ! dt / (dz X rho C)
     .    dzsnSV(ikl,isl) * ro__SV(ikl,isl) *Cn_dSV)      !
        END DO


! Snow
! ----

      DO   isl=1,nsno
        DO ikl=1,knonv
          ro__SV(ikl,isl) =                               !
     .                   ro__SV(ikl ,isl)                 !
     .       * max(0,min(isnoSV(ikl)-isl+1,1))            !
        END DO
      END DO

      DO   isl=1,nsno
        DO ikl=1,knonv

! Combined Heat Conductivity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          mu_aux      =  CdidSV                           !
     .                 *(ro__SV(ikl,isl) /rhoWat) ** 1.88 !
          mu__dz(ikl,isl) =                               !
     .      2.                        *mu_aux*mu_sno(ikl) ! Combined Heat
     .     /max(eps6,dzsnSV(ikl,isl  )*mu_sno(ikl)        ! Conductivity
     .              +dzsnSV(ikl,isl-1)*mu_aux     )       ! For upper Layer
          mu_sno(ikl)     =            mu_aux             !

! Inverted Heat Capacity
! ^^^^^^^^^^^^^^^^^^^^^^
          dtC_sv(ikl,isl) = dt__SV/max(eps__3,            ! dt / (dz X rho C)
     .    dzsnSV(ikl,isl) * ro__SV(ikl,isl) *Cn_dSV)      !
        END DO
      END DO


! Uppermost Effective Layer: NO conduction
! ----------------------------------------

        DO ikl=1,knonv
          mu__dz(ikl,isnoSV(ikl)+1) = 0.0
        END DO


! OUTPUT/Verification: Energy/Water Budget: Energy Budget (IN) 
! #e1   DO ikl=1,knonv
! #e1     ETSo_0(ikl) = 0.
! #e1   END DO
! #e1 DO   isl= -nsol,nsno
! #e1   DO ikl=1,knonv
! #e1     Exist0      = isl -           isnoSV(ikl)
! #e1     Exist0      = 1.  - max(zer0,min(un_1,Exist0))
! #e1     ETSo_0(ikl) = ETSo_0(ikl)
! #e1.                +(TsisSV(ikl,isl)-Tf_Sno)*Exist0
! #e1.                                 /dtC_sv(ikl,isl)
! #e1   END DO
! #e1 END DO


! Tridiagonal Elimination: Set Up
! ===============================

! Soil/Snow Interior
! ^^^^^^^^^^^^^^^^^^
      DO   isl= -nsol+1,nsno-1
        DO ikl=1,knonv
          Elem_A          =  dtC_sv(ikl,isl)         *mu__dz(ikl,isl)
          Elem_C          =  dtC_sv(ikl,isl)         *mu__dz(ikl,isl+1)
          Diag_A(ikl,isl) = -Elem_A  *Implic
          Diag_C(ikl,isl) = -Elem_C  *Implic
          Diag_B(ikl,isl) =  1.0d+0  -Diag_A(ikl,isl)-Diag_C(ikl,isl)
          Term_D(ikl,isl) =  Explic *(Elem_A         *TsisSV(ikl,isl-1)
     .                               +Elem_C         *TsisSV(ikl,isl+1))
     .             +(1.0d+0 -Explic *(Elem_A+Elem_C))*TsisSV(ikl,isl)
     .  + dtC_sv(ikl,isl)           * sol_SV(ikl)    *SoSosv(ikl)
     .                                              *(sEX_sv(ikl,isl+1)
     .                                               -sEX_sv(ikl,isl  ))
        END DO
      END DO

! Soil  lowest Layer
! ^^^^^^^^^^^^^^^^^^
           isl= -nsol
        DO ikl=1,knonv
          Elem_A          =  0.
          Elem_C          =  dtC_sv(ikl,isl)         *mu__dz(ikl,isl+1)
          Diag_A(ikl,isl) =  0.
          Diag_C(ikl,isl) = -Elem_C  *Implic
          Diag_B(ikl,isl) =  1.0d+0  -Diag_A(ikl,isl)-Diag_C(ikl,isl)
          Term_D(ikl,isl) =  Explic * Elem_C         *TsisSV(ikl,isl+1)
     .             +(1.0d+0 -Explic * Elem_C)        *TsisSV(ikl,isl)
     .  + dtC_sv(ikl,isl)           * sol_SV(ikl)    *SoSosv(ikl)
     .                                              *(sEX_sv(ikl,isl+1)
     .                                               -sEX_sv(ikl,isl  ))
        END DO

! Snow highest Layer (dummy!)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^
           isl=  nsno
        DO ikl=1,knonv
          Elem_A          =  dtC_sv(ikl,isl)  *mu__dz(ikl,isl)
          Elem_C          =  0.
          Diag_A(ikl,isl) = -Elem_A  *Implic
          Diag_C(ikl,isl) =  0.
          Diag_B(ikl,isl) =  1.0d+0  -Diag_A(ikl,isl)
          Term_D(ikl,isl) =  Explic * Elem_A  *TsisSV(ikl,isl-1)
     .             +(1.0d+0 -Explic * Elem_A) *TsisSV(ikl,isl)
     .  + dtC_sv(ikl,isl) * (sol_SV(ikl)      *SoSosv(ikl)
     .                                       *(sEX_sv(ikl,isl+1)
     .                                        -sEX_sv(ikl,isl  )))
        END DO

! Surface: UPwardIR Heat Flux
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,knonv
          isl         = isnoSV(ikl)
          dIRsdT(ikl) = Eso_sv(ikl)* StefBo          * 4.      ! - d(IR)/d(T)
     .                             * TsisSV(ikl,isl)           !
     .                             * TsisSV(ikl,isl)           !
     .                             * TsisSV(ikl,isl)           !
          IRs__D(ikl) = dIRsdT(ikl)* TsisSV(ikl,isl) * 0.75    !

! Surface: Richardson Number:   T Derivative
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c #RC     dRidTs(ikl) =-Grav_F      *    za__SV(ikl)
c #RC.                              *(1.-Sigmsv(ikl))
c #RC.                /(TaT_SV(ikl) *    VV__SV(ikl)
c #RC.                              *    VV__SV(ikl))

! Surface: Turbulent Heat Flux: Factors
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          f_HSHL(ikl) = rhT_SV(ikl) *(1.-Sigmsv(ikl))          !#common factor
     .                              /    rah_sv(ikl)           ! to  HS, HL
          f___HL(ikl) = f_HSHL(ikl) *    Lx_H2O(ikl)

! Surface: Sensible  Heat Flux: T Derivative
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          dSdTSV(ikl) = f_HSHL(ikl) *    CpdAir                !#- d(HS)/d(T)
c #RC.         *(1.0  -(TsisSV(ikl,isl) -TaT_SV(ikl))          !#Richardson
c #RC.         * dRidTs(ikl)*dFh_sv(ikl)/rah_sv(ikl))          ! Nb. Correct.
          HS___D(ikl) = dSdTSV(ikl) *    TaT_SV(ikl)           !

! Surface: Latent    Heat Flux: Saturation Specific Humidity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          den_qs      =         TsisSV(ikl,isl)- 35.8          !
          arg_qs      = 17.27 *(TsisSV(ikl,isl)-273.16)        !
     .                                   / den_qs              !
          qsatsg(ikl) = .0038 *        exp(arg_qs)             !
          dqs_dT(ikl) = qsatsg(ikl)* 4099.2   /(den_qs *den_qs)!
          fac_dt(ikl) = f_HSHL(ikl)/(rhoWat   * dz_dSV(0))     !
        END DO

! Surface: Latent    Heat Flux: Surface    Relative Humidity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              xgpsrf       =   1.05                            !
              agpsrf       = dt__SV*(   1.0-xgpsrf        )    !
     .                             /(   1.0-xgpsrf**nt_srf)    !
              dt_srf       = agpsrf                            !
              dt_ver       = 0.                                !
            DO ikl=1,knonv
              isl          =          isnoSV(ikl)              !
              etaBAK(ikl)  = max(eps6,eta_SV(ikl ,isl))        !
              etaNEW(ikl)  =          etaBAK(ikl)              !
              etEuBk(ikl)  =          etaNEW(ikl)              !
            END DO                                             !
        DO it_srf=1,nt_srf                                     !
              dt_ver       = dt_ver     +dt_srf                !
            DO ikl=1,knonv                                     !
              faceta(ikl)  = fac_dt(ikl)*dt_srf                !
c #VX         faceta(ikl)  = faceta(ikl)                       !
c #VX.                  /(1.+faceta(ikl)*dQa_SV(ikl))          !    Limitation
                                                               ! by Atm.Conten
!    .        *max(0,sign(1.,qsatsg(ikl)-QaT_SV(ikl))))        ! NO Limitation
                                                               ! of Downw.Flux
            END DO                                             !
          DO itEuBk=1,2                                        !
            DO ikl=1,knonv
              ist    = max(0,isotSV(ikl)-100*isnoSV(ikl))      ! 0 if    H2O
                                                               !
              Psi(ikl) =                                       !
     .                psidSV(ist)                              ! DR97, Eqn 3.34
     .              *(etadSV(ist)                              !
     .           /max(etEuBk(ikl),eps6))                       !
     .              **bCHdSV(ist)                              !
              PsiArg(ikl) = 7.2E-5*Psi(ikl)                    !
              RHuSol(ikl) =   exp(-min(ea_Max,PsiArg(ikl)))    !
              SHuSol(ikl) =     qsatsg(ikl)  *RHuSol(ikl)      ! DR97, Eqn 3.15
              etEuBk(ikl) =                                    !
     .       (etaNEW(ikl) + faceta(ikl)*(QaT_SV(ikl)           !
     .                                  -SHuSol(ikl)           !
     .                    *(1.          -bCHdSV(ist)           !
     .                                  *PsiArg(ikl))       )) !
     .      /(1.          + faceta(ikl)* SHuSol(ikl)           !
     .                                  *bCHdSV(ist)           !
     .                                  *PsiArg(ikl)           !
     .                                  /etaNEW(ikl))          !
              etEuBk(ikl) = etEuBk(ikl) -Rootsv(ikl,0)         !
     .                                 /(rhoWat*dz_dSV(0))     !
            END DO                                             !
          END DO                                               !
            DO ikl=1,knonv                                     !
              etaNEW(ikl) =  max(etEuBk(ikl),eps6)             !
            END DO                                             !
              dt_srf      =      dt_srf         * xgpsrf       !
        END DO                                                 !
 
! Surface: Latent    Heat Flux: Soil/Water Surface Contributions
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,knonv                                         !
          isl        =  isnoSV(ikl)                            !
          ist   = max(0,isotSV(ikl)-100*isnoSV(ikl))           ! 0 if    H2O
          ist__s= min(1,ist)                                   ! 1 if no H2O
          ist__w=     1-ist__s                                 ! 1 if    H2O
          d__eta     =  eta_SV(ikl,isl)-etaNEW(ikl)            !
          HL___D(ikl)=( ist__s *rhoWat *dz_dSV(0)              ! Soil Contrib.
     .                *(etaNEW(ikl)    -etaBAK(ikl)) / dt__SV  !
     .                 +ist__w         *f_HSHL(ikl)            ! H2O  Contrib.
     .                *(QaT_SV(ikl)    -qsatsg(ikl))         ) !
     .                * Lx_H2O(ikl)                            ! common factor

c #DL     RHuSol(ikl) =(QaT_SV(ikl)                            !
c #DL.                 -HL___D(ikl)    / f___HL(ikl))          !
c #DL.                / qsatsg(ikl)                            !
 
! Surface: Latent    Heat Flux: T Derivative
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          dLdTSV(ikl) = 0.
c #DL     dLdTSV(ikl) = f___HL(ikl) * RHuSol(ikl) *dqs_dT(ikl) ! - d(HL)/d(T)
c #DL     HL___D(ikl) = HL___D(ikl)                            !
c #DL.                 +dLdTSV(ikl) * TsisSV(ikl,isl)          !
        END DO                                                 !

! Surface: Tridiagonal Matrix Set Up
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,knonv
          isl             =  isnoSV(ikl)
          TSurf0(ikl)     =  TsisSV(ikl,isl)
          Elem_A          =  dtC_sv(ikl,isl)*mu__dz(ikl,isl)
          Elem_C          =  0.
          Diag_A(ikl,isl) = -Elem_A *Implic
          Diag_C(ikl,isl) =  0.
          Diag_B(ikl,isl) =  1.0d+0 -Diag_A(ikl,isl)
          Diag_B(ikl,isl) =  Diag_B(ikl,isl)
     .  + dtC_sv(ikl,isl) * (dIRsdT(ikl)                       ! Upw. Sol IR
     .                      +dSdTSV(ikl)                       ! HS/Surf.Contr.
     .                      +dLdTSV(ikl))                      ! HL/Surf.Contr.
          Term_D(ikl,isl) =  Explic *Elem_A *TsisSV(ikl,isl-1)
     .             +(1.0d+0 -Explic *Elem_A)*TsisSV(ikl,isl)
          Term_D(ikl,isl) =  Term_D(ikl,isl)
     .  + dtC_sv(ikl,isl) * (sol_SV(ikl)    *SoSosv(ikl)       ! Absorbed
     .                                     *(sEX_sv(ikl,isl+1) ! Solar
     .                                      -sEX_sv(ikl,isl  ))!
     .        +      tau_sv(ikl)      *IRd_SV(ikl)*Eso_sv(ikl) ! Down Atm IR
     .         -(1.0-tau_sv(ikl)) *0.5*IRv_sv(ikl)             ! Down Veg IR
     .                                +IRs__D(ikl)             ! Upw. Sol IR
     .                                +HS___D(ikl)             ! HS/Atmo.Contr.
     .                                +HL___D(ikl)            )! HL/Atmo.Contr.
        END DO


! Tridiagonal Elimination
! =======================

! Forward  Sweep
! ^^^^^^^^^^^^^^
        DO ikl=  1,knonv
          Aux__P(ikl,-nsol) = Diag_B(ikl,-nsol)
          Aux__Q(ikl,-nsol) =-Diag_C(ikl,-nsol)/Aux__P(ikl,-nsol)
        END DO

      DO   isl=-nsol+1,nsno
        DO ikl=      1,knonv
          Aux__P(ikl,isl)   = Diag_A(ikl,isl)  *Aux__Q(ikl,isl-1)
     .                       +Diag_B(ikl,isl)
          Aux__Q(ikl,isl)   =-Diag_C(ikl,isl)  /Aux__P(ikl,isl)
        END DO
      END DO

        DO ikl=      1,knonv
          TsisSV(ikl,-nsol) = Term_D(ikl,-nsol)/Aux__P(ikl,-nsol)
        END DO

      DO   isl=-nsol+1,nsno
        DO ikl=      1,knonv
          TsisSV(ikl,isl)   =(Term_D(ikl,isl)
     .                       -Diag_A(ikl,isl)  *TsisSV(ikl,isl-1))
     .                                         /Aux__P(ikl,isl)
        END DO
      END DO

! Backward Sweep
! ^^^^^^^^^^^^^^
      DO   isl=nsno-1,-nsol,-1
        DO ikl=     1,knonv
          TsisSV(ikl,isl)   = Aux__Q(ikl,isl)  *TsisSV(ikl,isl+1)
     .                                         +TsisSV(ikl,isl)
        END DO
      END DO

! Temperature Limits (avoids problems in case of no Snow Layers)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=     1,knonv
           isl              = isnoSV(ikl)
          dTSurf            = TsisSV(ikl,isl) -     TSurf0(ikl)
          TsisSV(ikl,isl)   = TSurf0(ikl) + sign(1.,dTSurf) ! 180.0 dgC/hr
     .              * min(abs(dTSurf),5.e-2*dt__SV)         ! =0.05 dgC/s
          IF (abs(dTSurf) > 5.e-2*dt__SV) THEN
            write(*,*) 'abrupt',ikl,'dTs ',dTSurf,TsisSV(ikl,isl) 
          ENDIF
        END DO
      
      DO ikl=     1,knonv
        DO   isl=isnoSV(ikl),1      ,-1

          IF (Ts_Min >  TsisSV(ikl,isl)) THEN
             write(*,*) 'cold', ikl, 'couche',isl, TsisSV(ikl,isl)
          ENDIF 
          IF (Ts_Max <  TsisSV(ikl,isl)) THEN
             write(*,*) 'hot ', ikl, 'couche',isl, TsisSV(ikl,isl) 
          ENDIF
          TsisSV(ikl,isl)   = max(Ts_Min,       TsisSV(ikl,isl))
          TsisSV(ikl,isl)   = min(Ts_Max,       TsisSV(ikl,isl))
        END DO
      END DO


! Update Surface    Fluxes
! ========================

        DO ikl=      1,knonv
          isl         = isnoSV(ikl)
          IRs_SV(ikl) = IRs__D(ikl)                          !
     .                - dIRsdT(ikl) * TsisSV(ikl,isl)        !
          HSs_sv(ikl) = HS___D(ikl)                          ! Sensible Heat
     .                - dSdTSV(ikl) * TsisSV(ikl,isl)        ! Downward > 0
          HLs_sv(ikl) = HL___D(ikl)                          ! Latent   Heat
     .                - dLdTSV(ikl) * TsisSV(ikl,isl)        ! Downward > 0
!hj220711
          TsfnSV(ikl) = TsisSV(ikl,isl)

! OUTPUT for Stand Alone NetCDF File
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c #NC     SOsoKL(ikl) = sol_SV(ikl) * SoSosv(ikl)            ! Absorbed Sol.
c #NC     IRsoKL(ikl) =               IRs_SV(ikl)            ! Up Surf. IR
c #NC.        +     tau_sv(ikl)      *IRd_SV(ikl)*Eso_sv(ikl)! Down Atm IR
c #NC.        -(1.0-tau_sv(ikl)) *0.5*IRv_sv(ikl)            ! Down Veg IR
c #NC     HSsoKL(ikl) = HSs_sv(ikl)                          ! HS
c #NC     HLsoKL(ikl) = HLs_sv(ikl)                          ! HL
c #NC     HLs_KL(ikl) = HLs_sv(ikl) / LhvH2O                 ! mm w.e./sec
        END DO


! OUTPUT/Verification: Energy/Water Budget: Energy Budget (OUT)
! #e1   DO ikl=1,knonv
! #e1     ETSo_d(ikl) =
! #e1.         (     SoSosv(ikl)      *sol_SV(ikl)             ! Net   Solar
! #e1.         +                       IRs_SV(ikl)             ! Up Surf. IR
! #e1.         +     tau_sv(ikl)      *IRd_SV(ikl)*Eso_sv(ikl) ! Down Atm IR
! #e1.         -(1.0-tau_sv(ikl)) *0.5*IRv_sv(ikl)             ! Down Veg IR
! #e1.                                +HSs_sv(ikl)             ! Sensible
! #e1.                                +HLs_sv(ikl)            )! Latent
! #e1     ETSo_1(ikl) = 0.
! #e1   END DO
! #e1 DO   isl= -nsol,nsno
! #e1   DO ikl=1,knonv
! #e1     Exist0      = isl -           isnoSV(ikl)
! #e1     Exist0      = 1.  - max(zer0,min(un_1,Exist0))
! #e1     ETSo_1(ikl) = ETSo_1(ikl)
! #e1.                +(TsisSV(ikl,isl)-Tf_Sno)*Exist0
! #e1.                                 /dtC_sv(ikl,isl)
! #e1   END DO
! #e1 END DO


      return
      end
      subroutine SISVAT_TS2                                                     
c #ES.                     (ETSo_0,ETSo_1,ETSo_d)                               
                                                                                
C +------------------------------------------------------------------------+    
C | MAR          SISVAT_TS2                            Mon 16-08-2009  MAR |    
C |   SubRoutine SISVAT_TS2 computes the Soil/Snow temperature and fluxes  |
C |   using the same method as in LMDZ for consistency.                    |
C |   The corresponding LMDZ routines are soil (soil.F90) and calcul_fluxs |
C |   (calcul_fluxs_mod.F90).                                              |   
C +------------------------------------------------------------------------+    
C |                                                                        |    
C |                                                                        |    
C |   PARAMETERS:  klonv: Total Number of columns =                        |    
C |   ^^^^^^^^^^        = Total Number of grid boxes of surface type       |    
C |                       (land ice for now)                               |    
C |                                                                        |    
C |   INPUT:   isnoSV   = total Nb of Ice/Snow Layers                      |    
C |   ^^^^^    sol_SV   : Downward Solar Radiation                  [W/m2] |    
C |            IRd_SV   : Surface Downward  Longwave   Radiation    [W/m2] |       
C |            VV__SV   : SBL Top    Wind Speed                      [m/s] |    
C |            TaT_SV   : SBL Top    Temperature                       [K] |     
C |            QaT_SV   : SBL Top    Specific  Humidity            [kg/kg] |      
C |            dzsnSV   : Snow Layer Thickness                         [m] |        
C |            dt__SV   : Time Step                                    [s] |    
C |                                                                        |    
C |            SoSosv   : Absorbed Solar Radiation by Surfac.(Normaliz)[-] |    
C |   ?        IRv_sv   : Vegetation  IR Radiation                  [W/m2] |    
C |   ?        tau_sv   : Fraction of Radiation transmitted by Canopy  [-] |        
C |            Eso_sv   : Soil+Snow       Emissivity                   [-] |    
C |   ?        rah_sv   : Aerodynamic Resistance for Heat            [s/m] |    
C |                                                                        | 
C |            dz1_SV    : "inverse" layer thickness (centered)      [1/m] |
C |            dz2_SV    : layer thickness (layer above (?))           [m] |
C |            AcoHSV    : coefficient for Enthalpy evolution, from atm.   |
C |            AcoHSV    : coefficient for Enthalpy evolution, from atm.   |
C |            AcoQSV    : coefficient for Humidity evolution, from atm.   | 
C |            BcoQSV    : coefficient for Humidity evolution, from atm.   |
C |            ps__SV    : surface pressure                           [Pa] |
C |            p1l_SV    : 1st atmospheric layer pressure             [Pa] |     
C |            cdH_SV    : drag coeff Energy (?)                           |
C |            rsolSV    : Radiation balance surface                [W/m2] |
C |            lambSV    : Coefficient for soil layer geometry         [-] |
C |                                                                        |    
C |   INPUT /  TsisSV   : Soil/Ice Temperatures (layers -nsol,-nsol+1,..,0)|    
C |   OUTPUT:           & Snow     Temperatures (layers  1,2,...,nsno) [K] |    
C |   ^^^^^^   rsolSV   : Radiation balance surface                 [W/m2] |     
C |                                                                        |    
C |   OUTPUT:  IRs_SV   : Soil      IR Radiation                    [W/m2] |    
C |   ^^^^^^   HSs_sv   : Sensible  Heat Flux                       [W/m2] |    
C |            HLs_sv   : Latent    Heat Flux                       [W/m2] |  
C |            TsfnSV   : new surface temperature                      [K] |
C |            Evp_sv   : Evaporation                              [kg/m2] |
C |            dSdTSV   : Sensible Heat Flux temp. derivative     [W/m2/K] |
C |            dLdTSV   : Latent Heat Flux temp. derivative       [W/m2/K] |
C |                                                                        | 
C |   ?        ETSo_0   : Snow/Soil Energy Power, before Forcing    [W/m2] |    
C |   ?        ETSo_1   : Snow/Soil Energy Power, after  Forcing    [W/m2] |    
C |   ?        ETSo_d   : Snow/Soil Energy Power         Forcing    [W/m2] |    
C |                                                                        | 
C |________________________________________________________________________|

      USE VAR_SV
      USE VARdSV
 
      USE VARySV                      
      USE VARtSV                      
      USE VARxSV 
      USE VARphy
      USE YOMCST_SISVAT
      USE indice_sol_mod


      IMPLICIT NONE                                                                    
                                                                                
                                                                                
C +--Global Variables                                                           
C +  ================                                                          

      INCLUDE "YOETHF.h"
      INCLUDE "FCTTRE.h"
!      INCLUDE "indicesol.h"
      INCLUDE "comsoil.h"
!      include  "LMDZphy.inc"   
                                                                        
C +--OUTPUT for Stand Alone NetCDF File                                         
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                         
c #NC      real*8        SOsoKL(klonv)             ! Absorbed Solar Radiation        
c #NC      real*8        IRsoKL(klonv)             ! Absorbed IR    Radiation        
c #NC      real*8        HSsoKL(klonv)             ! Absorbed Sensible Heat Flux     
c #NC      real*8        HLsoKL(klonv)             ! Absorbed Latent   Heat Flux     
c #NC      real*8        HLs_KL(klonv)             ! Evaporation                     
c #NC      real*8        HLv_KL(klonv)             ! Transpiration                   
c #NC      common/DumpNC/SOsoKL,IRsoKL                                               
c #NC     .             ,HSsoKL,HLsoKL                                               
c #NC     .             ,HLs_KL,HLv_KL          

C +--Internal Variables                                                         
C +  ==================    

      integer ig,jk,isl
      real mu     
      real Tsrf(klonv)               ! surface temperature as extrapolated from soil
      real mug(klonv)                 !hj coef top layers 
      real ztherm_i(klonv),zdz2(klonv,-nsol:nsno),z1s
      real pfluxgrd(klonv), pcapcal(klonv), cal(klonv)
      real beta(klonv), dif_grnd(klonv)
      real C_coef(klonv,-nsol:nsno),D_coef(klonv,-nsol:nsno)

      REAL, DIMENSION(klonv)   :: zx_mh, zx_nh, zx_oh
      REAL, DIMENSION(klonv)   :: zx_mq, zx_nq, zx_oq
      REAL, DIMENSION(klonv)   :: zx_pkh, zx_dq_s_dt, zx_qsat, zx_coef
      REAL, DIMENSION(klonv)   :: zx_sl, zx_k1
      REAL, DIMENSION(klonv)   :: d_ts
      REAL                     :: zdelta, zcvm5, zx_qs, zcor, zx_dq_s_dh
      REAL                     :: qsat_new, q1_new
C      REAL, PARAMETER          :: t_grnd = 271.35, t_coup = 273.15
C      REAL, PARAMETER          :: max_eau_sol = 150.0
      REAL, DIMENSION(klonv)   :: IRs__D, dIRsdT


      REAL t_grnd                      ! not used
      parameter(t_grnd = 271.35)       !
      REAL t_coup                      ! distinguish evap/sublimation
      parameter(t_coup = 273.15)       !
      REAL max_eau_sol
      parameter(max_eau_sol = 150.0)


!        write(*,*)'T check'
!     
!        DO  ig = 1,knonv  
!            DO  jk = 1,isnoSV(ig) !nsno       
!              IF (TsisSV(ig,jk) <= 1.) THEN          !hj check
!                TsisSV(ig,jk) = TsisSV(ig,isnoSV(ig))
!              ENDIF  
!
!              IF (TsisSV(ig,jk) <= 1.) THEN          !hj check
!                TsisSV(ig,jk) = 273.15
!              ENDIF   
!            END DO
!        END DO     

C!=======================================================================
C! I. First part: corresponds to soil.F90 in LMDZ
C!=======================================================================

      DO ig = 1,knonv      
        DO jk =1,isnoSV(ig)
          dz2_SV(ig,jk)=dzsnSV(ig,jk)
C! use arithmetic center between layers to derive dz1 for snow layers for simplicity: 
          dz1_SV(ig,jk)=2./(dzsnSV(ig,jk)+dzsnSV(ig,jk-1)) 
        ENDDO
      ENDDO

      DO ig = 1,knonv
        ztherm_i(ig)   = inertie_lic
        IF (isnoSV(ig) > 0) ztherm_i(ig)   = inertie_sno
      ENDDO

C!-----------------------------------------------------------------------
C! 1)
C! Calculation of Cgrf and Dgrd coefficients using soil temperature from 
C! previous time step.
C!
C! These variables are recalculated on the local compressed grid instead 
C! of saved in restart file.
C!-----------------------------------------------------------------------
      DO ig=1,knonv
        DO jk=-nsol,nsno
          zdz2(ig,jk)=dz2_SV(ig,jk)/dt__SV                                !ptimestep
        ENDDO
      ENDDO
 
      DO ig=1,knonv
        z1s = zdz2(ig,-nsol)+dz1_SV(ig,-nsol+1)
        C_coef(ig,-nsol+1)=zdz2(ig,-nsol)*TsisSV(ig,-nsol)/z1s
        D_coef(ig,-nsol+1)=dz1_SV(ig,-nsol+1)/z1s
      ENDDO

      DO ig=1,knonv
        DO jk=-nsol+1,isnoSV(ig)-1,1
          z1s = 1./(zdz2(ig,jk)+dz1_SV(ig,jk+1)+dz1_SV(ig,jk)           &
     &          *(1.-D_coef(ig,jk)))
          C_coef(ig,jk+1)=                                              &
     &          (TsisSV(ig,jk)*zdz2(ig,jk)                              &
     &          +dz1_SV(ig,jk)*C_coef(ig,jk)) * z1s
          D_coef(ig,jk+1)=dz1_SV(ig,jk+1)*z1s
        ENDDO
      ENDDO

C!-----------------------------------------------------------------------
C! 2)
C! Computation of the soil temperatures using the Cgrd and Dgrd
C! coefficient computed above
C!
C!-----------------------------------------------------------------------
C! Extrapolate surface Temperature                   !hj check
      mu=1./((2.**1.5-1.)/(2.**(0.5)-1.)-1.)

!     IF (knonv>0) THEN 
!      DO ig=1,8 
!        write(*,*)ig,'sisvat: Tsis ',TsisSV(ig,isnoSV(ig))
!        write(*,*)'max-1            ',TsisSV(ig,isnoSV(ig)-1)
!        write(*,*)'max-2            ',TsisSV(ig,isnoSV(ig)-2)
!        write(*,*)'0                ',TsisSV(ig,0)
!!        write(*,*)min(max(isnoSV(ig),0),1),max(1-isnoSV(ig),0)
!      ENDDO
!     END IF

      DO ig=1,knonv    
        IF (isnoSV(ig).GT.0) THEN
          IF (isnoSV(ig).GT.1) THEN
           mug(ig)=1./(1.+dzsnSV(ig,isnoSV(ig)-1)/dzsnSV(ig,isnoSV(ig))) !mu
          ELSE
           mug(ig) = 1./(1.+dzsnSV(ig,isnoSV(ig)-1)/dz_dSV(0)) !mu
          ENDIF
        ELSE
          mug(ig) = lambSV
        ENDIF

        IF (mug(ig)  .LE. 0.05) THEN
         write(*,*)'Attention mu low', mug(ig)
        ENDIF
        IF (mug(ig)  .GE. 0.98) THEN
         write(*,*)'Attention mu high', mug(ig)
        ENDIF

        Tsrf(ig)=(1.5*TsisSV(ig,isnoSV(ig))-0.5*TsisSV(ig,isnoSV(ig)-1))&
     &           *min(max(isnoSV(ig),0),1)+                             &
     &           ((mug(ig)+1)*TsisSV(ig,0)-mug(ig)*TsisSV(ig,-1))       &
     &           *max(1-isnoSV(ig),0)  
      ENDDO

  

C! Surface temperature
      DO ig=1,knonv
      TsisSV(ig,isnoSV(ig))=(mug(ig)*C_coef(ig,isnoSV(ig))+Tsf_SV(ig))/ &
     &        (mug(ig)*(1.-D_coef(ig,isnoSV(ig)))+1.)
      ENDDO
  
C! Other temperatures
      DO ig=1,knonv
        DO jk=isnoSV(ig),-nsol+1,-1
          TsisSV(ig,jk-1)=C_coef(ig,jk)+D_coef(ig,jk)                   &
     &          *TsisSV(ig,jk)
        ENDDO
      ENDDO
C      write(*,*)ig,'Tsis',TsisSV(ig,0)

C      IF (indice == is_sic) THEN
C        DO ig = 1,knonv
C          TsisSV(ig,-nsol) = RTT - 1.8
C        END DO
C      ENDIF

CC      !hj new 11 03 2010
        DO ig=1,knonv                                                          
          isl         = isnoSV(ig)                                             
C          dIRsdT(ig) = Eso_sv(ig)* stefan   * 4.                        & ! - d(IR)/d(T)   
C     &                             * Tsf_SV(ig)                         & !T TsisSV(ig,isl)           !                
C     &                             * Tsf_SV(ig)                         & !TsisSV(ig,isl)           !                
C     &                             * Tsf_SV(ig) !TsisSV(ig,isl)           !                
C          IRs__D(ig) = dIRsdT(ig)* Tsf_SV(ig) * 0.75 !TsisSV(ig,isl) * 0.75    !: 
          dIRsdT(ig) = Eso_sv(ig)* stefan   * 4.                        & ! - d(IR)/d(T)   
     &                             * TsisSV(ig,isl)                     & !
     &                             * TsisSV(ig,isl)                     & !
     &                             * TsisSV(ig,isl)                     & !
          IRs__D(ig) = dIRsdT(ig)* TsisSV(ig,isl) * 0.75                  !: 
         END DO
       !hj
C!-----------------------------------------------------------------------
C! 3)
C! Calculate the Cgrd and Dgrd coefficient corresponding to actual soil 
C! temperature
C!-----------------------------------------------------------------------
      DO ig=1,knonv
        z1s = zdz2(ig,-nsol)+dz1_SV(ig,-nsol+1)
        C_coef(ig,-nsol+1) = zdz2(ig,-nsol)*TsisSV(ig,-nsol)/z1s
        D_coef(ig,-nsol+1) = dz1_SV(ig,-nsol+1)/z1s
      ENDDO

      DO ig=1,knonv  
        DO jk=-nsol+1,isnoSV(ig)-1,1
          z1s = 1./(zdz2(ig,jk)+dz1_SV(ig,jk+1)+dz1_SV(ig,jk)           &
     &          *(1.-D_coef(ig,jk)))
          C_coef(ig,jk+1) = (TsisSV(ig,jk)*zdz2(ig,jk)+                 &
     &          dz1_SV(ig,jk)*C_coef(ig,jk)) * z1s
          D_coef(ig,jk+1) = dz1_SV(ig,jk+1)*z1s
        ENDDO
      ENDDO

C!-----------------------------------------------------------------------
C! 4)
C! Computation of the surface diffusive flux from ground and
C! calorific capacity of the ground
C!-----------------------------------------------------------------------
      DO ig=1,knonv
C! (pfluxgrd)
        pfluxgrd(ig) = ztherm_i(ig)*dz1_SV(ig,isnoSV(ig))*              &
     &      (C_coef(ig,isnoSV(ig))+(D_coef(ig,isnoSV(ig))-1.)           &
     &      *TsisSV(ig,isnoSV(ig)))
C! (pcapcal)
        pcapcal(ig)  = ztherm_i(ig)*                                    &
     &      (dz2_SV(ig,isnoSV(ig))+dt__SV*(1.-D_coef(ig,isnoSV(ig)))    &
     &      *dz1_SV(ig,isnoSV(ig)))
        z1s = mug(ig)*(1.-D_coef(ig,isnoSV(ig)))+1.
        pcapcal(ig)  = pcapcal(ig)/z1s
        pfluxgrd(ig) = ( pfluxgrd(ig)                                   &
     &       + pcapcal(ig) * (TsisSV(ig,isnoSV(ig)) * z1s               &
     &       - mug(ig)* C_coef(ig,isnoSV(ig))                           &
     &       - Tsf_SV(ig))       /dt__SV )
      ENDDO
        
      
      cal(1:knonv) = RCPD / pcapcal(1:knonv)
      rsolSV(1:knonv)  = rsolSV(1:knonv) + pfluxgrd(1:knonv)
C!=======================================================================
C! II. Second part: corresponds to calcul_fluxs_mod.F90 in LMDZ
C!=======================================================================

      Evp_sv = 0.
c #NC HSsoKL=0.
c #NC HLsoKL=0.
      dSdTSV = 0.
      dLdTSV = 0.
	
      beta(:) = 1.0
      dif_grnd(:) = 0.0

C! zx_qs = qsat en kg/kg
C!**********************************************************************x***************
!      write(*,*)'RCPD',RCPD,'RLVTT',RLVTT,'RD',RD,'RVTMP2',RVTMP2
      DO ig = 1,knonv
        IF (ps__SV(ig).LT.1.) THEN
!          write(*,*)'ig',ig,'ps',ps__SV(ig)
          ps__SV(ig)=max(ps__SV(ig),1.e-8)
        ENDIF
        IF (p1l_SV(ig).LT.1.) THEN
!          write(*,*)'ig',ig,'p1l',p1l_SV(ig)
          p1l_SV(ig)=max(p1l_SV(ig),1.e-8)
        ENDIF
        IF (TaT_SV(ig).LT.180.) THEN
!          write(*,*)'ig',ig,'TaT',TaT_SV(ig)
          TaT_SV(ig)=max(TaT_SV(ig),180.)
        ENDIF
        IF (QaT_SV(ig).LT.1.e-8) THEN
!          write(*,*)'ig',ig,'QaT',QaT_SV(ig)
          QaT_SV(ig)=max(QaT_SV(ig),1.e-8)
        ENDIF
        IF (Tsf_SV(ig).LT.100.) THEN
!          write(*,*)'ig',ig,'Tsf',Tsf_SV(ig)
          Tsf_SV(ig)=max(Tsf_SV(ig),180.)
        ENDIF
        IF (Tsf_SV(ig).GT.500.) THEN
!          write(*,*)'ig',ig,'Tsf',Tsf_SV(ig)
          Tsf_SV(ig)=min(Tsf_SV(ig),400.)
        ENDIF
!         IF (Tsrf(ig).LT.1.) THEN
!!          write(*,*)'ig',ig,'Tsrf',Tsrf(ig)  
!           Tsrf(ig)=max(Tsrf(ig),TaT_SV(ig)-20.)
!         ENDIF
         IF (cdH_SV(ig).LT.1.e-10) THEN
!           IF (ig.le.3)   write(*,*)'ig',ig,'cdH',cdH_SV(ig)
           cdH_SV(ig)=.5
         ENDIF
      ENDDO  
!      DO ig=1,3
!        write(*,*)'isnoSV',isnoSV(ig),'TsisSV',TsisSV(ig,isnoSV(ig))
!      ENDDO

      !write(*,*)'retv',retv,'thermcep',thermcep,'r2es',r2es
      !IF (r2es.LT.1.e-8) THEN
      !  r2es=1.e-8
      !ENDIF
      !IF (retv.LT.1.e-8) THEN
      !  retv=1.e-8
      !ENDIF

      DO ig = 1,knonv
        zx_pkh(ig) = 1. ! (ps__SV(ig)/ps__SV(ig))**RKAPPA
        IF (thermcep) THEN
          zdelta=MAX(0.,SIGN(1.,rtt-Tsf_SV(ig)))
          zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
          zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*QaT_SV(ig))
          zx_qs= r2es * FOEEW(Tsf_SV(ig),zdelta)/ps__SV(ig)
          zx_qs=MIN(0.5,zx_qs)
          !write(*,*)'zcor',retv*zx_qs
          zcor=1./(1.-retv*zx_qs)
          zx_qs=zx_qs*zcor
          zx_dq_s_dh = FOEDE(Tsf_SV(ig),zdelta,zcvm5,zx_qs,zcor)        &
     &          /RLVTT / zx_pkh(ig)
        ELSE
          IF (Tsf_SV(ig).LT.t_coup) THEN
             zx_qs = qsats(Tsf_SV(ig)) / ps__SV(ig)
             zx_dq_s_dh = dqsats(Tsf_SV(ig),zx_qs)/RLVTT                &
     &             / zx_pkh(ig)
          ELSE
             zx_qs = qsatl(Tsf_SV(ig)) / ps__SV(ig)
             zx_dq_s_dh = dqsatl(Tsf_SV(ig),zx_qs)/RLVTT                &
     &             / zx_pkh(ig)
          ENDIF
        ENDIF
        zx_dq_s_dt(ig) = RCPD * zx_pkh(ig) * zx_dq_s_dh
        zx_qsat(ig) = zx_qs
C        zx_coef(ig) = cdH_SV(ig) *                                     &
C     &       (1.0+SQRT(u1lay(ig)**2+v1lay(ig)**2)) *                   &
C     &       p1l_SV(ig)/(RD*t1lay(ig))
        zx_coef(ig) = cdH_SV(ig) *                                      &
     &       (1.0+VV__SV(ig)) *                                         &
     &       p1l_SV(ig)/(RD*TaT_SV(ig))
       
      ENDDO


C! === Calcul de la temperature de surface ===
C! zx_sl = chaleur latente d'evaporation ou de sublimation
C!****************************************************************************************

      DO ig = 1,knonv
        zx_sl(ig) = RLVTT
        IF (Tsf_SV(ig) .LT. RTT) zx_sl(ig) = RLSTT
        zx_k1(ig) = zx_coef(ig)
      ENDDO
    

      DO ig = 1,knonv
C! Q
        zx_oq(ig) = 1. - (beta(ig) * zx_k1(ig) * BcoQSV(ig) * dt__SV)
        zx_mq(ig) = beta(ig) * zx_k1(ig) *                              &
     &       (AcoQSV(ig) - zx_qsat(ig) +                                &
     &       zx_dq_s_dt(ig) * Tsf_SV(ig))                               &
     &       / zx_oq(ig)
        zx_nq(ig) = beta(ig) * zx_k1(ig) * (-1. * zx_dq_s_dt(ig))       &
     &       / zx_oq(ig)
       
C! H
        zx_oh(ig) = 1. - (zx_k1(ig) * BcoHSV(ig) * dt__SV)
        zx_mh(ig) = zx_k1(ig) * AcoHSV(ig) / zx_oh(ig)
        zx_nh(ig) = - (zx_k1(ig) * RCPD * zx_pkh(ig))/ zx_oh(ig)
     
C! surface temperature
        TsfnSV(ig) = (Tsf_SV(ig) + cal(ig)/RCPD * zx_pkh(ig) * dt__SV * &
     &       (rsolSV(ig) + zx_mh(ig) + zx_sl(ig) * zx_mq(ig))           & 
     &       + dif_grnd(ig) * t_grnd * dt__SV)/                         &
     &       ( 1. - dt__SV * cal(ig)/(RCPD * zx_pkh(ig)) *              &
     &       (zx_nh(ig) + zx_sl(ig) * zx_nq(ig))                        &  
     &       + dt__SV * dif_grnd(ig))

!hj rajoute 22 11 2010 tuning...
        TsfnSV(ig) = min(RTT+0.02,TsfnSV(ig))
       
        d_ts(ig) = TsfnSV(ig) - Tsf_SV(ig)


C!== flux_q est le flux de vapeur d'eau: kg/(m**2 s)  positive vers bas
C!== flux_t est le flux de cpt (energie sensible): j/(m**2 s)
        Evp_sv(ig) = - zx_mq(ig) - zx_nq(ig) * TsfnSV(ig) 
        HLs_sv(ig) = - Evp_sv(ig) * zx_sl(ig)
        HSs_sv(ig) = zx_mh(ig) + zx_nh(ig) * TsfnSV(ig)
       
C! Derives des flux dF/dTs (W m-2 K-1):
        dSdTSV(ig) = zx_nh(ig)
        dLdTSV(ig) = zx_sl(ig) * zx_nq(ig)


!hj  new 11 03 2010                                      
        isl         = isnoSV(ig)
!        TsisSV(ig,isl) = TsfnSV(ig)                                         
        IRs_SV(ig) = IRs__D(ig)                                         &!     
     &                - dIRsdT(ig) * TsfnSV(ig) !TsisSV(ig,isl)?      !  

! hj
c #NC   SOsoKL(ig) = sol_SV(ig) * SoSosv(ig)              ! Absorbed Sol.    
c #NC   IRsoKL(ig) =               IRs_SV(ig)                           & !Up Surf. IR
c #NC&        +     tau_sv(ig)      *IRd_SV(ig)*Eso_sv(ig)              & !Down Atm IR
c #NC&        -(1.0-tau_sv(ig)) *0.5*IRv_sv(ig)            ! Down Veg IR  
c #NC   HLsoKL(ig) = HLs_sv(ig) 
c #NC   HSsoKL(ig) = HSs_sv(ig) 
c #NC   HLs_KL(ig) = Evp_sv(ig)

C! Nouvelle valeure de l'humidite au dessus du sol
        qsat_new=zx_qsat(ig) + zx_dq_s_dt(ig) * d_ts(ig)
        q1_new = AcoQSV(ig) - BcoQSV(ig)* Evp_sv(ig)*dt__SV
        QaT_SV(ig)=q1_new*(1.-beta(ig)) + beta(ig)*qsat_new 

      ENDDO
!      DO ig=1,3                                                                
!        write(*,*)'   lat  HF',HLs_sv(ig),'sens HF',HSs_sv(ig)
!        write(*,*)'   dlHF/dT',dLdTSV(ig),'dsHF/dT',dSdTSV(ig)
!      END DO  
!      write(*,*)'RCPD',RCPD,'dt',dt__SV,'t_grnd',t_grnd
      DO ig = 1,1
!        write(*,*)ig,'Tsfn: ',TsfnSV(ig),'Tsrf',Tsrf(ig)
!        write(*,*)'  cal',cal(ig),'radsol',rsolSV(ig)
!        write(*,*)'  mh',zx_mh(ig),'sl', zx_sl(ig),'pkh',zx_pkh(ig)
!        write(*,*)'  mq',zx_mq(ig),'difgrnd',dif_grnd(ig)
!        write(*,*)'  nh',zx_nh(ig),'sl', zx_sl(ig),'nq',zx_nq(ig)
!        write(*,*)'term1:',cal(ig)/RCPD*zx_pkh(ig)*dt__SV
!        write(*,*)'*',(rsolSV(ig) + zx_mh(ig) + zx_sl(ig) * zx_mq(ig)) 
!        write(*,*)'+',dif_grnd(ig) * t_grnd * dt__SV,' / '
!        write(*,*)'(1-',dt__SV * cal(ig)/(RCPD * zx_pkh(ig)),'*'
!        write(*,*)'*',(zx_nh(ig) + zx_sl(ig) * zx_nq(ig))
!        write(*,*)'+',dt__SV * dif_grnd(ig),')'
      ENDDO

      end ! subroutine SISVAT_TS2  


      subroutine SISVAT_qVg

!--------------------------------------------------------------------------+
!   MAR          SISVAT_qVg                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_qVg computes the Canopy Water  Balance             |
!                                    including  Root   Extraction          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   ivgtSV   = 0,...,12:   Vegetation Type                      |
!     ^^^^^               0:          Water, Solid or Liquid               |
!                                                                          |
!     INPUT:   rhT_SV   : SBL Top    Air  Density                  [kg/m3] |
!     ^^^^^    QaT_SV   : SBL Top    Specific  Humidity            [kg/kg] |
!                                                                          |
!              TvegSV   : Canopy   Temperature                         [K] |
!              rrCaSV   : Canopy     Water     Content             [kg/m2] |
!              rrMxsv   : Canopy Maximum Intercepted Rain          [kg/m2] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!              EvT_sv   : EvapoTranspiration                       [kg/m2] |
!              Sigmsv   : Canopy Ventilation Factor                    [-] |
!              glf_sv   : Green Leaf Fraction of NOT fallen Leaves     [-] |
!              LAIesv   : Leaf Area  Index (effective / transpiration) [-] |
!              psi_sv   : Soil       Water    Potential                [m] |
!              Khydsv   : Soil   Hydraulic    Conductivity           [m/s] |
!                                                                          |
!     INPUT /  psivSV   : Leaf       Water    Potential                [m] |
!     OUTPUT:                                                              |
!     ^^^^^^                                                               |
!                                                                          |
!     OUTPUT:  Rootsv   : Root Water Pump                        [kg/m2/s] |
!     ^^^^^^                                                               |
!                                                                          |
!     REMARK: Water Extraction by roots calibrated by Evapotranspiration   |
!     ^^^^^^  (computed in the Canopy Energy Balance)                      |
!                                                                          |
!     REFERENCE: DR97: Koen de Ridder thesis, UCL, 1997                    |
!     ^^^^^^^^^                                                            |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #KW: Root Water Flow slowed down by Soil Hydraulic Conductivity      |
!                                                                          |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV
      USE VARdSV
      USE VAR0SV

      USE VARxSV
      USE VARySV


      IMPLICIT NONE

! Internal Variables
! ==================

      integer   ikl   ,isl                    ! Grid Point, Layer Indices
      integer   nitmax,nit                    ! Iterations Counter
      real      PlantW(klonv)                 ! Plant  Water
      real      dPdPsi(klonv)                 ! Plant  Water psi Derivative
      real      psidif                        ! Soil-Canopy  Water Pot. Differ.
      real      Root_W                        ! Root   Water     Flow
      real      RootOK                        ! Roots take Water in Soil Layer
      real      d_psiv                        ! Canopy Water     Increment
      real      dpvMAX                        ! Canopy Water     Increment MAX
      real      BWater                        ! Imbalance of Canopy Water Budg.
      real      BW_MAX                        ! MAX Imbal.of Canopy Water Budg.
      real      BW_MIN                        ! MIN Imbal.of Canopy Water Budg.
      real      dBwdpv                        ! Derivativ.of Canopy Water Budg.
      real      Bswich                        ! Newton-Raphson         Switch
      real      psiv_0(klonv)                 ! Canopy Temperature,  Previous t
      real      EvFrac                        ! Condensat./Transpiration Switch
      real      den_qs,arg_qs,qsatvg          ! Canopy Saturat. Spec. Humidity
      real      EvTran                        ! EvapoTranspiration
      real      dEdpsi                        ! Evapotranspiration  Derivative
      real      Fac_Ev,FacEvT                 ! Evapotranspiration  Factors
      real      denomE                        ! Evapotranspiration  Denominator
      real      F_Stom                        ! Funct.  (Leaf Water Potential)
      real      dFdpsi                        ! Derivative  of F_Stom
      real      denomF                        ! Denominator of F_Stom
      real      F___OK                        ! (psi>psi_c) => F_Stom swich  ON
      real      R0Stom                        ! Minimum Stomatal Resistance
      real      R_Stom                        !         Stomatal Resistance
      real      dRdpsi                        ! Derivat.Stomatal Resistance
      real      numerR                        ! Numerat.Stomatal Resistance


! Internal DATA
! =============

      data      nitmax /   5     /            ! Maximum  Iterations    Number
      data      dpvMAX /  20.    /            ! Canopy   Water   Increment MAX
      data      BW_MIN /   4.e-8 /            ! MIN Imbal. of Surf.Energy Budg.


! Newton-Raphson Scheme
! =====================

      nit    = 0
  101 CONTINUE
      nit    = nit + 1
      BW_MAX = 0.


! W.Potential of the Previous Time Step
! -------------------------------------

        DO ikl=1,knonv
          psiv_0(ikl) = psivSV(ikl)


! Extraction of Soil Water through the Plant Roots
! ------------------------------------------------

          PlantW(ikl)     = 0.                              ! Plant Water
          dPdPsi(ikl)     = 0.                              ! Idem, Derivat.
        END DO
      DO   isl=-nsol,0
        DO ikl=1,knonv
          psidif = psivSV(ikl)-(DH_dSV(ivgtSV(ikl))         ! Soil-Canopy Water
     .                         +psi_sv(       ikl ,isl))    ! Potential  Diff.
          Root_W = rhoWat     * RF__SV(ivgtSV(ikl),isl)     ! If > 0, Contrib.
     .              /max(eps_21,PR_dSV(ivgtSV(ikl))         !     to Root Water
c #KW.                         +Khydsv(ikl   ,isl )*1.e-4   ! (DR97, eqn.3.20)
     .                                                   )  !
!         Pas de prise en compte de la resistance sol/racine dans proto-svat 
!         (DR97, eqn.3.20)
          RootOK =      max(zer0, sign(un_1,psidif))
          Rootsv(ikl,isl) = Root_W*max(zer0,psidif)         ! Root  Water
          PlantW(ikl)     = PlantW(ikl) + Rootsv(ikl ,isl)  ! Plant Water
          dPdPsi(ikl)     = dPdPsi(ikl) + RootOK*Root_W     ! idem, Derivat.
        END DO
      END DO


! Latent   Heat Flux
! ------------------

! Canopy Saturation Specific Humidity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,knonv
          den_qs =         TvegSV(ikl)- 35.8
          arg_qs = 17.27 *(TvegSV(ikl)-273.16)/den_qs
          qsatvg = .0038 * exp(arg_qs)

! Canopy Stomatal Resistance
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          R0Stom = min(         StodSV(ivgtSV(ikl))
     .                /max(eps6,glf_sv(       ikl)),StxdSV) ! Min Stomatal R.
          denomF =              pscdSV-psivSV(ikl)
          F___OK =     max(zer0,sign(un_1,denomF))
          denomF =     max(eps6,          denomF)           !
          F_Stom =              pscdSV  / denomF            ! F(Leaf Wat.Pot.)
          dFdpsi =             -F_Stom  / denomF            !
                                                            ! DR97, eqn. 3.22
          numerR = R0Stom / max(LAIesv(ikl), R0Stom/StxdSV) !
          R_Stom = numerR *                  F_Stom         ! Can.Stomatal R.
                                                            ! DR97, eqn. 3.21
          dRdpsi = R_Stom *                  dFdpsi         !

! Evaporation / Evapotranspiration
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          EvFrac = max(zer0, sign(un_1,QaT_SV(ikl)-qsatvg)) ! Condensation/
          EvFrac = EvFrac                                   ! Transpiration
     .       + (1.-EvFrac)     *rrCaSV(ikl)/ rrMxsv(ikl)    !        Switch
          Fac_Ev = rhT_SV(ikl) *Sigmsv(ikl)                 ! idem,  Factor
          denomE = rah_sv(ikl) +R_Stom     * Sigmsv(ikl)
          FacEvT = Fac_Ev * (1.-EvFrac)    / denomE         !
          EvTran = FacEvT     *(qsatvg     - QaT_SV(ikl))   ! EvapoTranspir.
          dEdpsi =(EvTran     / denomE)    * dRdpsi         ! EvT Derivative


! Imbalance  of the Canopy  Water  Budget
! ---------------------------------------

          BWater      =( PlantW(ikl)                        ! Available  Water
     .                 - EvTran     )*        F___OK        ! Transpired Water

          Bswich      = max(zer0,                           ! Newton-Raphson
     .                      sign(un_1,    abs(BWater)       !         Switch
     .                                       -BW_MIN))      !


! Derivative of the Canopy  Water  Budget
! ---------------------------------------

          dBwdpv      = dPdpsi(ikl)
     .                - dEdpsi
          dBwdpv      = sign(  un_1,    dBwdpv)             !
     .                *  max(eps_21,abs(dBwdpv))            ! 


! Update Canopy and Surface/Canopy Temperatures
! ---------------------------------------------

          d_psiv      = BWater       / dBwdpv               !
          d_psiv      =      sign(un_1,d_psiv)              ! Increment
     .                       *min( abs(d_psiv)     ,dpvMAX) ! Limitor
          psivSV(ikl) = psivSV(ikl)  - Bswich      *d_psiv  ! Newton-Raphson
          BW_MAX      = max(BW_MAX,abs(BWater))
        END DO


! Update Root Water Fluxes | := Evapotranspiration 
! ------------------------------------------------

      DO   isl=-nsol,0
        DO ikl=1,knonv
          Rootsv(ikl,isl) = Rootsv(ikl,isl)*EvT_SV(ikl)     ! Root  Water
     .                          /max(eps_21,PlantW(ikl))    !
        END DO
      END DO

c hj150311 for parallel      IF (BW_MAX.gt.BW_MIN.and.nit.lt.nitmax)    GO TO 101
      IF (                     nit.lt.nitmax)    GO TO 101                    
      return
      end


      subroutine SISVAT_qSn
     .                     (
! #e1.                      EqSn_0,EqSn_1,EqSn_d
! #m1.                     ,SIsubl,SImelt,SIrnof
     .                     )

!--------------------------------------------------------------------------+
!   MAR          SISVAT_qSn                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_qSn updates  the Snow Water Content                |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   isnoSV   = total Nb of Ice/Snow Layers                      |
!     ^^^^^                                                                |
!                                                                          |
!     INPUT:   TaT_SV   : SBL Top    Temperature                       [K] |
!     ^^^^^    dt__SV   : Time Step                                    [s] |
!                                                                          |
!     INPUT /  drr_SV   : Rain Intensity                         [kg/m2/s] |
!     OUTPUT:  dzsnSV   : Snow Layer Thickness                         [m] |
!     ^^^^^^   eta_SV   : Snow Water Content                       [m3/m3] |
!              ro__SV   : Snow/Soil  Volumic Mass                  [kg/m3] |
!              TsisSV   : Soil/Ice Temperatures (layers -nsol,-nsol+1,..,0)|
!                       & Snow     Temperatures (layers  1,2,...,nsno) [K] |
!                                                                          |
!     OUTPUT:  SWS_SV   : Surficial Water Status                           |
!     ^^^^^^                                                               |
!              EExcsv   : Snow Energy in Excess, initial Forcing    [J/m2] |
!              EqSn_d   : Snow Energy in Excess, remaining          [J/m2] |
!              EqSn_0   : Snow Energy, before Phase Change          [J/m2] |
!              EqSn_1   : Snow Energy, after  Phase Change          [J/m2] |
!              SIsubl   : Snow sublimed/deposed Mass             [mm w.e.] |
!              SImelt   : Snow Melted           Mass             [mm w.e.] |
!              SIrnof   : Surficial Water + Run OFF Change       [mm w.e.] |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Possibility                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                          |
!     #IB: OUTPUT: Ice-Sheet Surface Mass Balance  (on NetCDF File )       |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:           (PLEASE VERIFY before USE)          |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #SU: SLUSH : Alternative Parameterization                            |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_iii_jjj_n     | #E0: OUTPUT on ASCII  File (SISVAT Variables) |
!   #                      |(#E0  MUST BE PREPROCESSED BEFORE #e1 & #e2 !) |
!   # SISVAT_iii_jjj_n     | #e1: OUTPUT/Verification: Energy Conservation |
!   # SISVAT_iii_jjj_n     | #m1: OUTPUT/Verification: * Mass Conservation |
!                          |                                               |
!   # SISVAT_qSn.vm        | #vm: OUTPUT/Verification: Energy/Water Budget |
!                          |      unit 43, SubRoutine  SISVAT_qSn **ONLY** |
!   # SISVAT_qSn.vu        | #vu: OUTPUT/Verification: Slush  Parameteriz. |
!                          |      unit 44, SubRoutine  SISVAT_qSn **ONLY** |
!                          |                                               |
!   # stdout               | #s2: OUTPUT of SnowFall, Snow Buffer          |
!                          |      unit  6, SubRoutine  SISVAT_BSn, _qSn    |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV
      USE VARdSV
      USE VAR0SV

      USE VARxSV


! OUTPUT
! ------

      USE VARySV

      IMPLICIT NONE

! OUTPUT/Verification: Energy/Water Budget
! #e1 real     EqSn_d(klonv)                 ! Energy in Excess, initial
! #e1 real     EqSn_0(klonv)                 ! Snow Energy, befor Phase Change
! #vm real     EqSn01(klonv)                 ! Snow Energy, after Phase Change
! #vm real     EqSn02(klonv)                 ! Snow Energy, after Phase Change
                                             !              .AND. Last Melting 
! #e1 real     EqSn_1(klonv)                 ! Snow Energy, after Phase Change
                                             !              .AND. Mass Redistr.
! OUTPUT/Verification: * Mass Conservation
! #m1 real     SIsubl(klonv)                 ! Snow Deposed Mass
! #m1 real     SImelt(klonv)                 ! Snow Melted  Mass
! #m1 real     SIrnof(klonv)                 ! Local Surficial Water + Run OFF


! Internal Variables
! ==================

      integer  ikl   ,isn                    !
      integer  nh                            ! Non erodible Snow: up.lay.Index
      integer  LayrOK                        ! 1 (0)  if In(Above) Snow Pack
      integer  k_face                        ! 1 (0)  if Crystal(no) faceted
      integer  LastOK                        ! 1 ==>  1! Snow Layer
      integer  NOLayr                        ! 1     Layer  Update
      integer  noSnow(klonv)                 ! Nb of Layers Updater
      integer  kSlush                        ! Slush Switch
      real     dTSnow                        ! Temperature                  [C]
      real     EExdum(klonv)                 ! Energy in Excess when no Snow
      real     OKmelt                        ! 1 (0)  if        (no) Melting
      real     EnMelt                        ! Energy in excess, for Melting
      real     SnHLat                        ! Energy consumed   in  Melting
      real     AdEnrg,B_Enrg                 ! Additional Energy from  Vapor
      real     dzVap0,dzVap1                 ! Vaporized Thickness          [m]
      real     dzMelt(klonv)                 ! Melted    Thickness          [m]
      real     rosDry                        ! Snow volumic Mass if no Water in
      real     PorVol                        ! Pore volume
      real     PClose                        ! Pore Hole Close OFF Switch
      real     SGDiam                        !      Snow Grain Diameter
      real     SGDmax                        ! Max. Snow Grain Diameter 
      real     rWater                        ! Retained Water           [kg/m2]
      real     drrNEW                        ! New available Water      [kg/m2]
      real     rdzNEW                        ! Snow          Mass       [kg/m2]
      real     rdzsno                        ! Snow          Mass       [kg/m2]
      real     EnFrez                        ! Energy Release    in  Freezing
      real     WaFrez                        ! Water  consumed   in  Melting
      real     RapdOK                        ! 1. ==> Snow melts rapidly
      real     ThinOK                        ! 1. ==> Snow Layer is thin
      real     dzepsi                        ! Minim. Snow Layer Thickness (!)
      real     dz_Min                        ! Minim. Snow Layer Thickness
      real     z_Melt                        ! Last (thin) Layer Melting
      real     rusnew                        ! Surficial Water Thickness   [mm]
      real     zWater                        ! Max Slush Water Thickness   [mm]
      real     zSlush                        !     Slush Water Thickness   [mm]
      real     ro_new                        ! New Snow/ice  Density    [kg/m3]
      real     zc,zt                         ! Non erod.Snow Thickness[mm w.e.]

! OUTPUT of SISVAT Trace Statistics (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer             isnnew,isinew,isnUpD,isnitr

! OUTPUT/Verification: Energy/Water Budget
! #vm real     WqSn_0(klonv)                 ! Snow Water+Forcing  Initial
! #vm real     WqSn_1(klonv)                 ! Snow Water+Forcing, Final
! #vm logical         emopen                 ! IO   Switch
! #vm common/Se_qSn_L/emopen                 !
! #vm integer         no_err                 !
! #vm common/Se_qSn_I/no_err                 !
! #vm real     hourer,timeer                 !
! #vm common/Se_qSn_R/timeer                 !

! OUTPUT/Verification: Slush  Parameterization
! #vu logical         su_opn                 ! IO   Switch
! #vu common/SI_qSn_L/su_opn                 !


! DATA
! ====
!hj1907      data      dzepsi/0.005/     !hj180711           ! Minim. Snow Layer Thickness (!) 
!hjp230611  
      data      dzepsi/0.0001/                ! Minim. Snow Layer Thickness (!)
!hjp290611      data      dzepsi/0.0003/                ! Minim. Snow Layer Thickness (!)

!hjp230611
      data      dz_Min/1.e-4/                 ! Minim. Snow Layer Thickness
!     data      dz_Min/0.005/                 ! Minim. Snow Layer Thickness
c ... Warning: Too high for Col de Porte: precludes 1st snow (layer) apparition
!hj1907      data      dz_Min/0.005/    !hj180711              ! Minim. Snow Layer Thickness 


      data      SGDmax/0.003/                 ! Maxim. Snow Grain Diameter  [m]
                                              ! (Rowe et al. 1995, JGR p.16268)


! OUTPUT/Verification: Energy/Water Budget: Energy Budget (IN)
! #e1   DO ikl=1,knonv
! #e1     EqSn_0(ikl) = 0.
! #e1   END DO
! #e1 DO   isn=nsno,1,-1
! #e1   DO ikl=1,knonv
! #e1     EqSn_0(ikl) = EqSn_0(ikl) + ro__SV(ikl,isn) *dzsnSV(ikl,isn)
! #e1.                *(Cn_dSV      *(TsisSV(ikl,isn) -Tf_Sno         )
! #e1.                 -LhfH2O      *(1.              -eta_SV(ikl,isn)))
! #e1   END DO
! #e1 END DO


! OUTPUT/Verification: Energy/Water Budget: Water  Budget (IN) 
! #vm   DO ikl=1,knonv
! #vm     WqSn_0(ikl) = drr_SV(ikl) * dt__SV          
! #vm.                 +rusnSV(ikl)
! #vm   END DO
! #vm DO   isn=nsno,1,-1
! #vm   DO ikl=1,knonv
! #vm     WqSn_0(ikl) = WqSn_0(ikl) + ro__SV(ikl,isn) *dzsnSV(ikl,isn)
! #vm   END DO
! #vm END DO


! OUTPUT/Verification: * Mass Conservation
! #m1   DO ikl=1,knonv
! #m1     SImelt(ikl) = 0.
! #m1     SIrnof(ikl) = rusnSV(ikl) + RnofSV(ikl) * dt__SV
! #m1   END DO


! Initialization
! ==============

      DO ikl=1,knonv
        noSnow(ikl)   = 0                   ! Nb of Layers Updater
        ispiSV(ikl)   = 0                   ! Pore Hole Close OFF Index
                                            ! (assumed to be the Top of
                                            !  the surimposed Ice Layer)
c #IB   wem_SV(ikl)   = 0.
c #IB   wer_SV(ikl)   = 0.
      END DO


! Melting/Freezing Energy
! =======================

!  REMARK: Snow liquid Water Temperature assumed = Tf_Sno
!  ^^^^^^
        DO ikl=1,knonv
          EExdum(ikl) = drr_SV(ikl)     * hC_Wat *(TaT_SV(ikl)-Tf_Sno) 
     .                                  * dt__SV
          EExcsv(ikl) = EExdum(ikl)     *    min(1,isnoSV(ikl)) ! Snow exists
          EExdum(ikl) = EExdum(ikl)     -          EExcsv(ikl)  !

! OUTPUT/Verification: Energy/Water Budget
! #e1     EqSn_d(ikl) = EExcsv(ikl)                             !

        END DO


! Surficial Water Status
! ----------------------

        DO ikl=1,knonv
          SWS_SV(ikl) = max(zer0,sign(un_1,Tf_Sno
     .                                    -TsisSV(ikl,isnoSV(ikl))))
        END DO

      DO   isn=nsno,1,-1
        DO ikl=1,knonv

! Energy, store Previous Content
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          dTSnow      = TsisSV(ikl,isn) -          Tf_Sno
          EExcsv(ikl) = EExcsv(ikl) 
     .                + ro__SV(ikl,isn) * Cn_dSV * dTSnow
     .                                           * dzsnSV(ikl,isn)
          TsisSV(ikl,isn) =                        Tf_Sno

! Water,  store Previous Content
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          drr_SV(ikl) = drr_SV(ikl)
     .                + ro__SV(ikl,isn)          * eta_SV(ikl,isn) 
     .                                           * dzsnSV(ikl,isn)
     .                / dt__SV
          ro__SV(ikl,isn) =
     .                  ro__SV(ikl,isn) *(1.     - eta_SV(ikl,isn))
          eta_SV(ikl,isn) =  0.


! Melting  if EExcsv > 0
! ======================

          EnMelt      =    max(zer0,          EExcsv(ikl) )

! Energy Consumption
! ^^^^^^^^^^^^^^^^^^
          SnHLat      = ro__SV(ikl,isn) * LhfH2O
          dzMelt(ikl) = EnMelt      / max(SnHLat,    eps6 )
          noSnow(ikl) = noSnow(ikl) 
     .      + max(zer0  ,sign(un_1,dzMelt(ikl)                !
     .                            -dzsnSV(ikl ,isn)))         ! 1 if full Melt
     .       *min(1     , max(0 ,1+isnoSV(ikl)-isn))          ! 1 in the  Pack
          dzMelt(ikl) =    
     .              min(dzsnSV(ikl, isn),dzMelt(ikl))
          dzsnSV(ikl,isn) = 
     .                  dzsnSV(ikl,isn) -dzMelt(ikl)
          EExcsv(ikl) = EExcsv(ikl)     -dzMelt(ikl)*SnHLat
c #IB     wem_SV(ikl) = wem_SV(ikl)     -dzMelt(ikl)*ro__SV(ikl,isn)

! Water  Production
! ^^^^^^^^^^^^^^^^^
          drr_SV(ikl) = drr_SV(ikl)     
     .                + ro__SV(ikl,isn) * dzMelt(ikl)/dt__SV

! OUTPUT/Verification: * Mass Conservation
! #m1     SImelt(ikl) = SImelt(ikl)
! #m1.                + ro__SV(ikl,isn) * dzMelt(ikl)

          OKmelt      =max(zer0,sign(un_1,drr_SV(ikl)-eps6))

! Snow History
! ^^^^^^^^^^^^
          k_face          =       min(    istoSV(ikl,isn),istdSV(1)) ! = 1  if
     .                           *max(0,2-istoSV(ikl,isn)          ) ! faceted
          istoSV(ikl,isn) =                                          !
     .        (1.-OKmelt) *               istoSV(ikl,isn)            !
     .      +     OKmelt  *((1-k_face) *  istdSV(2)                  !
     .                     +   k_face  *  istdSV(3)      )           !


! Freezing if EExcsv < 0
! ======================

          rdzsno      =          ro__SV(ikl,isn) * dzsnSV(ikl ,isn)
          LayrOK      = min(   1, max(0          , isnoSV(ikl)-isn+1))
          EnFrez      = min(zer0,                  EExcsv(ikl))
          WaFrez      =   -(     EnFrez          * LayrOK / LhfH2O)
          drrNEW      = max(zer0,drr_SV(ikl)     - WaFrez / dt__SV)
          WaFrez      =    (     drr_SV(ikl)     - drrNEW)* dt__SV
          drr_SV(ikl) =          drrNEW
          EExcsv(ikl) =          EExcsv(ikl)     + WaFrez * LhfH2O
          EnFrez      = min(zer0,EExcsv(ikl))    * LayrOK
          rdzNEW      = WaFrez + rdzsno
          ro__SV(ikl,isn) =      rdzNEW /max(eps6, dzsnSV(ikl,isn))
          TsisSV(ikl,isn) =      Tf_Sno
     .                + EnFrez /(Cn_dSV *max(eps6, rdzNEW)        )
          EExcsv(ikl) =          EExcsv(ikl)     - EnFrez
c #IB     wer_SV(ikl) = WaFrez                             
c #IB.                + wer_SV(ikl)


! Snow Water Content
! ==================

! Pore   Volume [-]
! ^^^^^^^^^^^^^^^^^
          rosDry      =(1.     - eta_SV(ikl,isn))* ro__SV(ikl,isn) !
          PorVol      = 1.     - rosDry          / rhoIce          !
          PorVol      =      max(PorVol          , zer0  )         !

! Water  Retention
! ^^^^^^^^^^^^^^^^
          rWater      = ws0dSV * PorVol * rhoWat * dzsnSV(ikl,isn)
          drrNEW      = max(zer0,drr_SV(ikl)     - rWater /dt__SV)
          rWater      =    (     drr_SV(ikl)     - drrNEW)*dt__SV
          drr_SV(ikl)     =      drrNEW
          rdzNEW          =      rWater 
     .                         + rosDry          * dzsnSV(ikl,isn)
          eta_SV(ikl,isn) =      rWater / max(eps6,rdzNEW)
          ro__SV(ikl,isn) =      rdzNEW / max(eps6,dzsnSV(ikl,isn))

! Pore Hole Close OFF
! ^^^^^^^^^^^^^^^^^^^
          PClose = max(zer0,
     .                 sign(un_1,ro__SV(ikl,isn)
     .                          -roCdSV         ))
          ispiSV(ikl) =          ispiSV(ikl)      *(1.-PClose)
     .                +      max(ispiSV(ikl),isn)    * Pclose
          PClose = max(0   ,                         ! Water under SuPer.Ice
     .                 min (1   ,ispiSV(ikl)         ! contributes to
     .                          -isn            ))   ! Surficial   Water
          rusnSV(ikl) =          rusnSV(ikl) 
     .                +          drr_SV(ikl) *dt__SV * PClose
          drr_SV(ikl) =          drr_SV(ikl)      *(1.-PClose)

        END DO
      END DO


! Remove Zero-Thickness Layers
! ============================

 1000 CONTINUE
 !          isnitr =          0
      DO   ikl=1,knonv
           isnUpD =          0
           isinew =          0
        DO isn=1,nsno-1
           isnnew =(un_1-max(zer0  ,sign(un_1,dzsnSV(ikl,isn)-dzepsi)))
     .             *     max(0     , min(1   ,isnoSV(ikl) +1 -isn ))
           isnUpD =      max(isnUpD,          isnnew)
!           isnitr =      max(isnitr,          isnnew)
           isinew =      isn*isnUpD *max(0, 1-isinew)      ! LowerMost  0-Layer
     .                                       +isinew       ! Index
           dzsnSV(ikl,isn) =                  dzsnSV(ikl,isn+isnnew)
           ro__SV(ikl,isn) =                  ro__SV(ikl,isn+isnnew)
           TsisSV(ikl,isn) =                  TsisSV(ikl,isn+isnnew)
           eta_SV(ikl,isn) =                  eta_SV(ikl,isn+isnnew)
           G1snSV(ikl,isn) =                  G1snSV(ikl,isn+isnnew)
           G2snSV(ikl,isn) =                  G2snSV(ikl,isn+isnnew)
           dzsnSV(ikl,isn+isnnew) =(1-isnnew)*dzsnSV(ikl,isn+isnnew)
           ro__SV(ikl,isn+isnnew) =(1-isnnew)*ro__SV(ikl,isn+isnnew)
           eta_SV(ikl,isn+isnnew) =(1-isnnew)*eta_SV(ikl,isn+isnnew)
           G1snSV(ikl,isn+isnnew) =(1-isnnew)*G1snSV(ikl,isn+isnnew)
           G2snSV(ikl,isn+isnnew) =(1-isnnew)*G2snSV(ikl,isn+isnnew)
        END DO
           isnoSV(ikl)   =   isnoSV(ikl)-isnUpD            ! Nb of Snow   Layer
           ispiSV(ikl)   =   ispiSV(ikl)                   ! Nb of SuperI Layer
     .    -isnUpD *max(0,min(ispiSV(ikl)-isinew,1))        ! Update  if I=0

      END DO
!      IF  (isnitr.GT.0)                                       GO TO 1000


! New upper Limit of the non erodible Snow (istoSV .GT. 1)
! ========================================

      DO   ikl=1,knonv
           nh =     0
        DO isn=  nsno,1,-1
           nh =    nh + isn* min(istoSV(ikl,isn)-1,1)*max(0,1-nh)
        ENDDO
           zc =     0. 
           zt =     0. 
        DO isn=1,nsno
           zc =    zc +          dzsnSV(ikl,isn) *ro__SV(ikl,isn)
     .                     * max(0,min(1,nh+1-isn))
           zt =    zt +          dzsnSV(ikl,isn) *ro__SV(ikl,isn)
        END DO
           zWE_SV(ikl) =                 zt
           zWEcSV(ikl) = min(zWEcSV(ikl),zt)
           zWEcSV(ikl) = max(zWEcSV(ikl),zc)
      END DO


! OUTPUT/Verification: Energy/Water Budget: Energy Budget (OUT)
! #vm   DO ikl=1,knonv
! #vm     EqSn01(ikl) =-EqSn_0(ikl)
! #vm.                 -EExcsv(ikl)
! #vm   END DO
! #vm DO   isn=nsno,1,-1
! #vm   DO ikl=1,knonv
! #vm     EqSn01(ikl) = EqSn01(ikl) + ro__SV(ikl,isn) *dzsnSV(ikl,isn)
! #vm.                *(Cn_dSV      *(TsisSV(ikl,isn) -Tf_Sno         )
! #vm.                 -LhfH2O      *(1.              -eta_SV(ikl,isn)))
! #vm   END DO
! #vm END DO


! "Negative Heat" from supercooled rain
!  ------------------------------------

      DO ikl=1,knonv
          EExcsv(ikl) = EExcsv(ikl) + EExdum(ikl)


! Surficial Water Run OFF
! -----------------------

          rusnew      = rusnSV(ikl) * SWf_SV(ikl)
          RnofSV(ikl) = RnofSV(ikl)
     .                +(rusnSV(ikl) - rusnew     ) / dt__SV
          rusnSV(ikl) = rusnew
      END DO


! Percolation down the Continental Ice Pack
! -----------------------------------------

        DO ikl=1,knonv
          drr_SV(ikl) = drr_SV(ikl) + rusnSV(ikl)
     .                     * (1-min(1,ispiSV(ikl)))/ dt__SV 
          rusnSV(ikl) = rusnSV(ikl)
     .                     *    min(1,ispiSV(ikl))
        END DO


! Slush Formation (CAUTION: ADD RunOff Possibility before Activation)
! ---------------  ^^^^^^^  ^^^

! OUTPUT/Verification: Slush  Parameterization
! #vu IF (.NOT.su_opn)                                              THEN
! #vu          su_opn=.true.
! #vu     open(unit=44,status='unknown',file='SISVAT_qSn.vu')
! #vu     rewind    44
! #vu END IF
! #vu     write(44,440) daHost
 440      format('iSupI    i       dz       ro      eta',
     .            '   PorVol   zSlush     ro_n    eta_n',2x,a18)

c #SU DO   isn=1,nsno
c #SU   DO ikl=1,knonv
c #SU     kSlush = min(1,max(0,isn+1-ispiSV(ikl)))        ! Slush Switch

! Available Additional Pore   Volume [-]
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c #SU     PorVol = 1. - ro__SV(ikl,isn)                   ! [--]
c #SU.           *(1. - eta_SV(ikl,isn))/ rhoIce          !
c #SU.           -      eta_SV(ikl,isn)                   !
c #SU.                 *ro__SV(ikl,isn) / rhoWat          !
c #SU     PorVol =  max(PorVol          , zer0  )         !
c #SU     zWater =      dzsnSV(ikl,isn) * PorVol * 1000.  ! [mm] OR [kg/m2]
c #SU.           * (1. -SWS_SV(ikl)                       ! 0 <=> freezing
c #SU.                *(1 -min(1,iabs(isn-isnoSV(ikl))))) ! 1 <=> isn=isnoSV
c #SU     zSlush =  min(rusnSV(ikl)     , zWater)         ! [mm] OR [kg/m2]
c #SU     rusnSV(ikl) = rusnSV(ikl)     - zSlush          ! [mm] OR [kg/m2]
c #SU     ro_new      =(dzsnSV(ikl,isn) * ro__SV(ikl,isn) !
c #SU.                 +zSlush                           )!
c #SU.            / max(dzsnSV(ikl,isn) , eps6           )!

! OUTPUT/Verification: Slush  Parameterization
! #vu     rusnew          = eta_SV(ikl,isn)               !

c #SU     eta_SV(ikl,isn) =(ro_new - ro__SV(ikl,isn)      !
c #SU.                    *(1.     - eta_SV(ikl,isn)))    !
c #SU.               / max (ro_new , eps6            )    !

! OUTPUT/Verification: Slush  Parameterization
! #vu     IF    (isn.le.isnoSV(ikl))                      !
! #vu.    write(44,441) ispiSV(ikl),isn,dzsnSV(ikl,isn)   !
! #vu.                 ,ro__SV(ikl,isn),rusnew            !
! #vu.                 ,PorVol         ,zSlush            !
! #vu.                 ,ro_new         ,eta_SV(ikl,isn)   !
 441      format(2i5,f9.3,f9.1,f9.6,f9.3,f9.6,f9.1,f9.6)  !

c #SU     ro__SV(ikl,isn) =      ro_new                   !
c #SU   END DO
c #SU END DO


! Impact of the Sublimation/Deposition on the Surface Mass Balance
! ================================================================

      DO ikl=1,knonv
          isn                     = isnoSV(ikl)
          dzVap0                  =                   dt__SV
     .  * HLs_sv(ikl)         * min(isn             , 1   )
     .  /(Lx_H2O(ikl)         * max(ro__SV(ikl,isn) , eps6))
          NOLayr=min(zer0,sign(un_1,dzsnSV(ikl,isn) + dzVap0))
          dzVap1=min(zer0,          dzsnSV(ikl,isn) + dzVap0)


! Additional Energy (CAUTION: Verification is not performed)
! -----------------

! OUTPUT/Verification: Energy Consrv. (HLS)
! #e4     AdEnrg = dzVap0 * ro__SV(ikl,isnoSV(ikl))           ! Water   Vapor
! #e4.            *hC_Wat *(TsisSV(ikl,isnoSV(ikl)) -Tf_Sno)  ! Sensible Heat

! OUTPUT/Verification: Energy Consrv. (HL)
! #e3     B_Enrg =(Cn_dSV      *(TsisSV(ikl,isn) -Tf_Sno         )
! #e3.            -LhfH2O      *(1.              -eta_SV(ikl,isn)))
! #e3.           /(1.          + dzVap0 /max(eps6,dzsnSV(ikl,isn)))
! #e3     eta_SV(ikl,isn) =
! #e3.           max(zer0,un_1 +(B_Enrg 
! #e3.                         -(TsisSV(ikl,isn) -Tf_Sno)*Cn_dSV)
! #e3.                          /LhfH2O                          )
! #e3     TsisSV(ikl,isn) =    ( B_Enrg
! #e3.                         +(1.              -eta_SV(ikl,isn))
! #e3.                          *LhfH2O                          )
! #e3.                         / Cn_dSV
! #e3.                         + Tf_Sno

! OUTPUT/Verification: Energy Conservation
! #e1     STOP "PLEASE add Energy (#e3) from deposition/sublimation"


! Update of the upper Snow layer Thickness
! ----------------------------------------

          dzsnSV(ikl,isn) =
     .           max(zer0,  dzsnSV(ikl,isnoSV(ikl)) + dzVap0)
          isnoSV(ikl)     = isnoSV(ikl)             + NOLayr
          isn             = isnoSV(ikl)
          dzsnSV(ikl,isn) = dzsnSV(ikl,isn) + dzVap1
c #IB     wes_SV(ikl)     = ro__SV(ikl,isn) * dzVap0
      END DO


! OUTPUT/Verification: Energy/Water Budget: Energy Budget (OUT)
! #vm   DO ikl=1,knonv
! #vm     EqSn02(ikl) =-EqSn_0(ikl)
! #vm.                 -EExcsv(ikl)
! #vm   END DO
! #vm DO   isn=nsno,1,-1
! #vm   DO ikl=1,knonv
! #vm     EqSn02(ikl) = EqSn01(ikl) + ro__SV(ikl,isn) *dzsnSV(ikl,isn)
! #vm.                *(Cn_dSV      *(TsisSV(ikl,isn) -Tf_Sno         )
! #vm.                 -LhfH2O      *(1.              -eta_SV(ikl,isn)))
! #vm   END DO
! #vm END DO


! OUTPUT/Verification: * Mass Conservation
! #m1   DO ikl=1,knonv
! #m1     SIsubl(ikl) = dt__SV*HLs_sv(ikl)*min(isnoSV(ikl),1)
! #m1.                        /Lx_H2O(ikl)
! #m1     SIrnof(ikl) = rusnSV(ikl) + RnofSV(ikl) * dt__SV
! #m1.                - SIrnof(ikl)
! #m1   END DO


! Anticipated Disappearance of a rapidly Melting too thin Last Snow Layer
! =======================================================================

      DO ikl=1,knonv
        LastOK = min(1   , max(0   ,iiceSV(ikl)-isnoSV(ikl)+2)
     .                    *min(1   ,isnoSV(ikl)-iiceSV(ikl))
     .                    +min(1   ,isnoSV(ikl))              )
        RapdOK = max(zer0,sign(un_1,dzMelt(ikl)-eps6         ))
        ThinOK = max(zer0,sign(un_1,dz_Min     -dzsnSV(ikl,1)))
        z_Melt = LastOK     *RapdOK*ThinOK
        noSnow(ikl)   = noSnow(ikl)   + z_Melt
        z_Melt        =                 z_Melt *dzsnSV(ikl,1) 
        dzsnSV(ikl,1) = dzsnSV(ikl,1) - z_Melt
        EExcsv(ikl)   = EExcsv(ikl)   - z_Melt *ro__SV(ikl,1)
     .                                *(1.     -eta_SV(ikl,1))*LhfH2O

! Water  Production
! ^^^^^^^^^^^^^^^^^
        drr_SV(ikl)   = drr_SV(ikl)     
     .                + ro__SV(ikl,1) * z_Melt /dt__SV
      END DO


! Update Nb of Layers
! ===================

! OUTPUT in SISVAT for ikl = 1 (preferably for Stand Alone Version)
! OUTPUT           for SnowFall and Snow Buffer
! #s2   IF          (isnoSV(1) .GT. 0)
! #s2.  write(6,6005)noSnow(1)
 6005   format(i3,' (noSnow) ')

      DO ikl=1,knonv
        isnoSV(ikl)   = isnoSV(ikl) 
     .     * min(1,iabs(isnoSV(ikl)-noSnow(ikl)))
      END DO


! OUTPUT/Verification: Energy Conservation: Energy Budget (OUT)
! #e1   DO ikl=1,knonv
! #e1     EqSn_1(ikl) = 0.
! #e1   END DO
! #e1 DO   isn=nsno,1,-1
! #e1   DO ikl=1,knonv
! #e1     EqSn_1(ikl) = EqSn_1(ikl) + ro__SV(ikl,isn) *dzsnSV(ikl,isn)
! #e1.                *(Cn_dSV      *(TsisSV(ikl,isn) -Tf_Sno         )
! #e1.                 -LhfH2O      *(1.              -eta_SV(ikl,isn)))
! #e1   END DO
! #e1 END DO


! OUTPUT/Verification: Energy/Water Budget: Water  Budget (OUT)
! #vm   DO ikl=1,knonv
! #vm     WqSn_0(ikl) = WqSn_0(ikl)
! #vm.                + HLs_sv(ikl)    * dt__SV      
! #vm.             *min(isnoSV(ikl),1) / Lx_H2O(ikl)
! #vm     WqSn_1(ikl) = drr_SV(ikl)    * dt__SV
! #vm.                + rusnSV(ikl)
! #vm.                + RnofSV(ikl)    * dt__SV
! #vm   END DO
! #vm DO   isn=nsno,1,-1
! #vm   DO ikl=1,knonv
! #vm     WqSn_1(ikl) = WqSn_1(ikl)    
! #vm.                + ro__SV(ikl,isn)* dzsnSV(ikl,isn)
! #vm   END DO
! #vm END DO


! OUTPUT/Verification: Energy/Water Budget
! #vm IF (.NOT.emopen)                                              THEN
! #vm          emopen = .true.
! #vm          open(unit=43,status='unknown',file='SISVAT_qSn.vm')
! #vm          rewind 43
! #vm   write(43,43)
 43     format('SubRoutine SISVAT_qSn: Local Energy and Water Budgets',
     .       /,'=====================================================')
! #vm END IF
! #vm DO ikl=1,knonv
! #vm IF (EqSn01(ikl).gt.1.e-3) write(43,431) dahost,EqSn01(ikl)
 431  format(' WARNING (1) in _qSn,',         a18,
     .       ': Energy Unbalance in Phase Change = ',e15.6)
! #vm END DO
! #vm DO ikl=1,knonv
! #vm IF (EqSn02(ikl).gt.1.e-3) write(43,432) dahost,EqSn01(ikl)
 432  format(' WARNING (2) in _qSn,',         a18,
     .       ': Energy Unbalance in Phase Change = ',e15.6)
! #vm END DO
! #vm         timeer=timeer + dt__SV
! #vm         hourer=3600.0
! #vm IF (mod(no_err,11).eq.0)                                      THEN
! #vm         no_err=       1
! #vm   write(43,435)timeer/hourer
 435    format(11('-'),'----------+-',3('-'),'----------+-',
     .          3('-'),'----------+-',3('-'),'----------+-',
     .                 '----------------+----------------+',
     .       /,f8.2,3x,'EqSn_0(1) | ',3x,'EqSn_d(1) | ',
     .              3x,'EqSn_1(1) | ',3x,'EExcsv(1) | ',
     .                 'E_0+E_d-E_1-EE  |   Water Budget |',
     .       /,11('-'),'----------+-',3('-'),'----------+-',
     .          3('-'),'----------+-',3('-'),'----------+-',
     .                 '----------------+----------------+')
! #vm END IF
! #vm IF (abs(EqSn_0(1)+EqSn_d(1)-EqSn_1(1)-EExcsv(1)).gt.eps6.OR.
! #vm.    abs(WqSn_1(1)-WqSn_0(1))                    .gt.eps6    ) THEN
! #vm         no_err=no_err+1
! #vm   write(43,436)   EqSn_0(1),EqSn_d(1)
! #vm.                 ,EqSn_1(1),EExcsv(1)
! #vm.                 ,EqSn_0(1)+EqSn_d(1)-EqSn_1(1)-EExcsv(1)
! #vm.                 ,WqSn_1(1)-WqSn_0(1)
 436    format(8x,f12.0,' + ',f12.0,' - ',f12.0,' - ',f12.3,' = ',f12.3,
     .      '    | ',f15.9)
! #vm END IF

! OUTPUT/Verification: Energy Conservation
! #e1   DO ikl=1,knonv
! #e1     EqSn_d(ikl) = EqSn_d(ikl) - EExcsv(ikl)
! #e1   END DO

      return
      end


      subroutine SISVAT_GSn

!--------------------------------------------------------------------------+
!   MAR          SISVAT_GSn                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_GSn simulates SNOW Metamorphism                    |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT /  isnoSV   = total Nb of Ice/Snow Layers                      |
!     OUTPUT:  iiceSV   = total Nb of Ice      Layers                      |
!     ^^^^^^   istoSV   = 0,...,5 :   Snow     History (see istdSV data)   |
!                                                                          |
!     INPUT:   TsisSV   : Soil/Ice Temperatures (layers -nsol,-nsol+1,..,0)|
!     ^^^^^             & Snow     Temperatures (layers  1,2,...,nsno) [K] |
!              ro__SV   : Soil/Snow Volumic Mass                   [kg/m3] |
!              eta_SV   : Soil/Snow Water   Content                [m3/m3] |
!              slopSV   : Surface Slope                                [-] |
!              dzsnSV   : Snow Layer        Thickness                  [m] |
!              dt__SV   : Time  Step                                   [s] |
!                                                                          |
!     INPUT /  G1snSV   : Dendricity (<0) or Sphericity (>0) of Snow Layer |
!     OUTPUT:  G2snSV   : Sphericity (>0) or Size            of Snow Layer |
!     ^^^^^^                                                               |
!                                                                          |
!     Formalisme adopte pour la Representation des Grains:                 |
!     Formalism         for the Representation of  Grains:                 |
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                 |
!                                                                          |
!               1       - -1                 Neige Fraiche                 |
!              / \      |                    -------------                 |
!             /   \     |  Dendricite        decrite  par Dendricite       |
!            /     \    |  Dendricity                  et Sphericite       |
!           /       \   |                                                  |
!          2---------3  -  0                 described by Dendricity       |
!                                                     and Sphericity       |
!          |---------|                                                     |
!          0         1                                                     |
!          Sphericite                                                      |
!          Sphericity                                                      |
!                                                                          |
!          4---------5  -                                                  |
!          |         |  |                                                  |
!          |         |  |  Diametre (1/10eme de mm) (ou Taille)            |
!          |         |  |  Diameter (1/10th  of mm) (or Size  )            |
!          |         |  |                                                  |
!          |         |  |                    Neige non dendritique         |
!          6---------7  -                    ---------------------         |
!                                            decrite  par Sphericite       |
!                                                      et     Taille       |
!                                            described by Sphericity       |
!                                                     and       Size       |
!                                                                          |
!     Les Variables du Modele:                                             |
!     Model         Variables:                                             |
!     ^^^^^^^^^^^^^^^^^^^^^^^^                                             |
!       Cas Dendritique               Cas non Dendritique                  |
!                                                                          |
!       G1snSV        : Dendricite    G1snSV        : Sphericite           |
!       G2snSV        : Sphericite    G2snSV        : Taille (1/10e mm)    |
!                                                     Size                 |
!                                                                          |
!     Cas Dendritique/ Dendritic Case                                      |
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                      |
!     Dendricite(Dendricity) G1snSV                                        |
!              varie     de -G1_dSV (-99 par defaut / etoile)          a 0 |
!              division par -G1_dSV pour obtenir des valeurs entre 1  et 0 |
!              varies  from -G1_dSV (default -99    / fresh snow)     to 0 |
!              division  by -G1_dSV to obtain values       between 1 and 0 |
!                                                                          |
!     Sphericite(Sphericity) G2snSV                                        |
!              varie     de  0         (cas completement anguleux)         |
!                         a  G1_dSV (99 par defaut, cas spherique)         |
!              division par  G1_dSV pour obtenir des valeurs entre 0  et 1 |
!              varies  from  0      (full faceted)               to G1_dSV |
!                                                                          |
!     Cas non Dendritique / non Dendritic Case                             |
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                             |
!     Sphericite(Sphericity) G1snSV                                        |
!              varie     de  0         (cas completement anguleux)         |
!                         a  G1_dSV (99 par defaut, cas spherique)         |
!              division par  G1_dSV pour obtenir des valeurs entre 0  et 1 |
!              varies  from  0      (full faceted)               to G1_dSV |
!                                                                          |
!     Taille    (Size)       G2snSV                                        |
!              superieure a  ADSdSV (.4 mm) et ne fait que croitre         |
!              greater than  ADSdSV (.4 mm) always increases               |
!                                                                          |
!     Exemples: Points caracteristiques des Figures ci-dessus              |
!     ^^^^^^^^^                                                            |
!                                                                          |
!                 G1snSV    G2snSV     dendricite  sphericite  taille      |
!                                      dendricity  sphericity  size        |
!     ------------------------------------------------------------------   |
!                                                              [1/10 mm]   |
!       1        -G1_dSV    sph3SN          1           0.5                |
!       2           0         0             0           0                  |
!       3           0       G1_dSV          0           1                  |
!       4           0       ADSdSV                      0       4.         |
!       5         G1_dSV    ADSdSV-vsphe1               1       3.         |
!       6           0         --                        0       --         |
!       7         G1_dSV      --                        1       --         |
!                                                                          |
!     par defaut: G1_dSV=99.                                               |
!                           sph3SN=50.                                     |
!                           ADSdSV= 4.                                     |
!                                  vsphe1=1.                               |
!                                                                          |
!     Methode:                                                             |
!     ^^^^^^^^                                                             |
!     1. Evolution Types de Grains selon Lois de Brun et al. (1992):       |
!        Grain metamorphism according to         Brun et al. (1992):       |
!        Plusieurs Cas sont a distiguer  /  the different Cases are:       |
!         1.1 Metamorphose Neige humide  /  wet Snow                       |
!         1.2 Metamorphose Neige seche   /  dry Snow                       |
!           1.2.1 Gradient faible        /  low      Temperature Gradient  |
!           1.2.2 Gradient moyen         /  moderate Temperature Gradient  |
!           1.2.3 Gradient fort          /  high     Temperature Gradient  |
!        Dans chaque Cas on separe Neige Dendritique et non Dendritique    |
!                             le Passage Dendritique -> non Dendritique    |
!                             se fait lorsque  G1snSV devient > 0          |
!        the Case of Dentritic or non Dendritic Snow is treated separately |
!        the Limit   Dentritic -> non Dendritic is reached when G1snSV > 0 |
!                                                                          |
!     2. Tassement: Loi de Viscosite adaptee selon le Type de Grains       |
!        Snow Settling:    Viscosity depends on the   Grain Type           |
!                                                                          |
!     3. Update Variables historiques (cas non dendritique seulement)      |
!        nhSNow defaut                                                     |
!                  0    Cas normal                                         |
!        istdSV(1) 1    Grains anguleux / faceted cristal                  |
!        istdSV(2) 2    Grains ayant ete en presence d eau liquide         |
!                       mais n'ayant pas eu de caractere anguleux    /     |
!                       liquid water and no faceted cristals before        |
!        istdSV(3) 3    Grains ayant ete en presence d eau liquide         |
!                       ayant eu auparavant un caractere anguleux    /     |
!                       liquid water and    faceted cristals before        |
!                                                                          |
!     REFER. : Brun et al.      1989, J. Glaciol 35 pp. 333--342           |
!     ^^^^^^^^ Brun et al.      1992, J. Glaciol 38 pp.  13-- 22           |
!              (CROCUS Model, adapted to MAR at CEN by H.Gallee)           |
!                                                                          |
!     REFER. : Marbouty, D.     1980, J. Glaciol 26 pp. xxx--xxx           |
!     ^^^^^^^^ (CROCUS Model, adapted to MAR at CEN by H.Gallee)           |
!              (for angular shapes)                                        |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_GSn.vp        | #vp: OUTPUT/Verification: Snow   Properties   |
!                          |      unit 47, SubRoutines SISVAT_zSn, _GSn    |
!   # stdout               | #wp: OUTPUT/Verification: Snow   Properties   |
!                          |      unit  6, SubRoutine  SISVAT_GSn          |
!                                                                          |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV
      USE VARdSV
      USE VAR0SV


! INPUT/OUTPUT
! ------------

      USE VARxSV

      IMPLICIT NONE

! OUTPUT
! ------

      integer   k


! Local  Variables
! ================

      logical   vector                        !
      integer   ikl                           !
      integer   isn   ,isnp                   !
      integer   istoOK                        !
      real      G1_bak,G2_bak                 ! Old Values of G1, G2
      real      ro_dry(klonv,      nsno)      ! Dry Density            [g/cm3]
      real      etaSno(klonv,      nsno)      ! Liquid Water Content   [g/cm2]
      real      SnMass(klonv)                 ! Snow   Mass            [kg/m2]
      real      dTsndz                        ! Temperature Gradient
      real      sWater                        !        Water Content       [%]  
      real      exp1Wa                        !
      real      dDENDR                        ! Dendricity Increment
      real      DENDRn                        ! Normalized Dendricity
      real      SPHERn                        ! Normalized Sphericity
      real      Wet_OK                        ! Wet Metamorphism Switch
      real      OK__DE                        !
      real      OK__wd                        ! New G*, from wet Dendritic
      real      G1__wd                        ! New G1, from wet Dendritic
      real      G2__wd                        ! New G2, from wet Dendritic
      real      OKlowT                        !
      real      facVap                        !
      real      OK_ldd                        !
      real      G1_ldd                        !
      real      G2_ldd                        !
      real      DiamGx                        !
      real      DiamOK                        !
      real      No_Big                        !
      real      dSPHER                        !
      real      SPHER0                        !
      real      SPHbig                        !
      real      G1_lds                        !
      real      OK_mdT                        !
      real      OKmidT                        !
      real      OKhigT                        !
      real      OK_mdd                        !
      real      G1_mdd                        !
      real      G2_mdd                        !
      real      G1_mds                        !
      real      OK_hdd                        !
      real      G1_hdd                        !
      real      G2_hdd                        !
      real      OK_hds                        !
      real      G1_hds                        !
      real      T1__OK,T2__OK                 !
      real      T3_xOK,T3__OK,T3_nOK          !
      real      ro1_OK,ro2_OK                 !
      real      dT1_OK,dT2_OK,dT3xOK,dT3_OK   !
      real      dT4xOK,dT4_OK,dT4nOK,AngSno   !
      real      G2_hds,SphrOK,HISupd          !
      real      H1a_OK,H1b_OK,H1__OK          !
      real      H23aOK,H23bOK,H23_OK          !
      real      H2__OK,H3__OK                 !
      real      H45_OK,H4__OK,H5__OK          !
      real      ViscSn,OK_Liq,OK_Ang,OKxLiq   !
      real      dSnMas,dzsnew,rosnew,rosmax   !

      real      epsi5                         ! Alpha ev67 single precision	
      real      vdiam1                        ! Small Grains Min.Diam.[.0001m]
      real      vdiam2                        ! Spher.Variat.Max Diam.    [mm]
      real      vdiam3                        ! Min.Diam.|Limit Spher.    [mm]
      real      vdiam4                        ! Min.Diam.|Viscosity Change
      real      vsphe1                        ! Max Sphericity
      real      vsphe2                        ! Low T Metamorphism  Coeff.
      real      vsphe3                        ! Max.Sphericity (history=1)
      real      vsphe4                        ! Min.Sphericity=>history=1      
      real      vtang1,vtang2,vtang3,vtang4   ! Temperature Contribution
      real      vtang5,vtang6,vtang7,vtang8   !
      real      vtang9,vtanga,vtangb,vtangc   ! 
      real      vrang1,vrang2                 ! Density     Contribution
      real      vgang1,vgang2,vgang3,vgang4   ! Grad(T)     Contribution
      real      vgang5,vgang6,vgang7,vgang8   !
      real      vgang9,vganga,vgangb,vgangc   !
      real      vgran6                        ! Max.Sphericity for Settling
      real      vtelv1                        ! Threshold | history = 2, 3
      real       vvap1                        ! Vapor Pressure Coefficient
      real       vvap2                        ! Vapor Pressure Exponent
      real      vgrat1                        ! Boundary weak/mid   grad(T)
      real      vgrat2                        ! Boundary mid/strong grad(T)
      real         vfi                        ! PHI,         strong grad(T)
      real      vvisc1,vvisc2,vvisc3,vvisc4   ! Viscosity Coefficients
      real      vvisc5,vvisc6,vvisc7          ! id., wet Snow
      real      rovisc                        ! Wet Snow Density  Influence
      real        vdz3                        ! Maximum Layer Densification
      real      OK__ws                        ! New G2
      real      G1__ws                        ! New G1, from wet Spheric
      real      G2__ws                        ! New G2, from wet Spheric
      real      husi_0,husi_1,husi_2,husi_3   ! Constants for New G2
      real      vtail1,vtail2                 ! Constants for New G2
      real      frac_j                        ! Time Step            [Day]

      real       vdent1                       ! Wet Snow Metamorphism
      integer   nvdent1                       ! (Coefficients for
      integer   nvdent2                       !           Dendricity)

! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp real      G_curr(18),Gcases(18)
! #vp common   /GSnLOC/    Gcases

! OUTPUT/Verification: Snow   Properties
! #wp real                 D__MAX
! #wp common   /GSnMAX/    D__MAX


! DATA
! ====

      data       vector/.true./               ! Vectorization  Switch
      data       vdent1/ 2.e8/                ! Wet Snow Metamorphism
      data      nvdent1/ 3   /                ! (Coefficients for
      data      nvdent2/16   /                !           Dendricity)

      data       husi_0 /20.      /           !   10  * 2
      data       husi_1 / 0.23873 /           ! (3/4) /pi
      data       husi_2 / 4.18880 /           ! (4/3) *pi
      data       husi_3 / 0.33333 /           !  1/3
      data       vtail1 / 1.28e-08/           !  Wet Metamorphism
      data       vtail2 / 4.22e-10/           ! (NON Dendritic / Spheric)

      data       epsi5  / 1.0e-5  /           ! 

      data       vdiam1 / 4.0     /           ! Small Grains Min.Diameter

      data       vdiam2 / 0.5     /           ! Spher.Variat.Max Diam.[mm]
      data       vdiam3 / 3.0     /           ! Min.Diam.|Limit Spher.[mm]
      data       vdiam4 / 2.0     /           ! Min.Diam.|Viscosity Change

      data       vsphe1 / 1.0     /           ! Max Sphericity
      data       vsphe2 / 1.0e9   /           ! Low T Metamorphism  Coeff.
      data       vsphe3 / 0.5     /           ! Max.Sphericity (history=1)
      data       vsphe4 / 0.1     /           ! Min.Sphericity=>history=1      

      data       vgran6 / 51.     /           ! Max.Sphericity for Settling
      data       vtelv1 / 5.e-1   /           ! Threshold | history = 2, 3

      data        vvap1 /-6.e3    /           ! Vapor Pressure Coefficient
      data        vvap2 / 0.4     /           ! Vapor Pressure Exponent

      data       vgrat1 /0.05     /           ! Boundary weak/mid   grad(T)
      data       vgrat2 /0.15     /           ! Boundary mid/strong grad(T)
      data          vfi /0.09     /           ! PHI,         strong grad(T)

      data       vvisc1 / 0.70    /           ! Viscosity Coefficients
      data       vvisc2 / 1.11e5  /           !
      data       vvisc3 /23.00    /           !
      data       vvisc4 / 0.10    /           !
      data       vvisc5 / 1.00    /           ! id., wet Snow
      data       vvisc6 / 2.00    /           !
      data       vvisc7 /10.00    /           !
      data       rovisc / 0.25    /           ! Wet Snow Density  Influence
      data         vdz3 / 0.30    /           ! Maximum Layer Densification


! DATA (Coefficient Fonction fort Gradient Marbouty)
! --------------------------------------------------

      data       vtang1 /40.0/                ! Temperature Contribution
      data       vtang2 / 6.0/                !
      data       vtang3 /22.0/                !
      data       vtang4 / 0.7/                !
      data       vtang5 / 0.3/                !
      data       vtang6 / 6.0/                !
      data       vtang7 / 1.0/                !
      data       vtang8 / 0.8/                !
      data       vtang9 /16.0/                !
      data       vtanga / 0.2/                !
      data       vtangb / 0.2/                !
      data       vtangc /18.0/                !

      data       vrang1 / 0.40/               ! Density     Contribution
      data       vrang2 / 0.15/               !

      data       vgang1 / 0.70/               ! Grad(T)     Contribution
      data       vgang2 / 0.25/               !
      data       vgang3 / 0.40/               !
      data       vgang4 / 0.50/               !
      data       vgang5 / 0.10/               !
      data       vgang6 / 0.15/               !
      data       vgang7 / 0.10/               !
      data       vgang8 / 0.55/               !
      data       vgang9 / 0.65/               !
      data       vganga / 0.20/               !
      data       vgangb / 0.85/               !
      data       vgangc / 0.15/               !

! OUTPUT/Verification: Snow   Properties
! #wp data       D__MAX / 4.00/               !


!  1. Metamorphoses dans les Strates
!     Metamorphism
!     ==============================

      frac_j = dt__SV / 86400.                        ! Time Step [Day]


!  1.1 Initialisation: teneur en eau liquide et gradient de temperature 
!  ------------------  liquid water content and temperature gradient

      DO   isn=1,nsno
        DO ikl=1,knonv
          ro_dry(ikl,isn) = 1.e-3 *ro__SV(ikl,isn)    ! Dry Density
     .                    *(1.    -eta_SV(ikl,isn))   !         [g/cm3]
          etaSno(ikl,isn) = 1.e-1 *dzsnSV(ikl,isn)    ! Liquid Water
     .                    *        ro__SV(ikl,isn)    ! Content [g/cm2]
     .                    *        eta_SV(ikl,isn)    !
        END DO
      END DO

      DO   isn=1,nsno
        DO ikl=1,knonv

          isnp   = min(isn+1,isnoSV(ikl))

          dTsndz = abs( (TsisSV(ikl,isnp)-TsisSV(ikl,isn-1)) *2.e-2
     .            /max(((dzsnSV(ikl,isnp)+dzsnSV(ikl,isn)  )
     .                 *(           isnp -           isn)
     .                 +(dzsnSV(ikl,isn )+dzsnSV(ikl,isn-1))),eps6))
!         Factor 1.d-2 for Conversion K/m --> K/cm


!  1.2 Metamorphose humide 
!      Wet Snow Metamorphism
!      ---------------------

          Wet_OK = max(zer0,sign(un_1,eta_SV(ikl,isn)-eps6))


!      Vitesse de diminution de la dendricite
!      Rate of the dendricity decrease
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          sWater=1.d-1*ro__SV(ikl,isn)*eta_SV(ikl,isn)
     .       /max(eps6,ro_dry(ikl,isn))
! .    sWater:Water Content [%]  
!             1.d-1= 1.d2(1->%) * 1.d-3(ro__SV*eta_SV:kg/m3->g/cm3) 

          exp1Wa=   sWater**nvdent1
          dDENDR=max(exp1Wa/nvdent2,vdent1*exp(vvap1/Tf_Sno))

!  1.2.1 Cas dendritique/dendritic Case
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          OK__wd=max(zer0,                            !
     .               sign(un_1,-G1snSV(ikl,isn)       !
     .                         -eps6           ))     !

          DENDRn=-G1snSV(ikl,isn)/G1_dSV  ! Normalized Dendricity (+)
          SPHERn= G2snSV(ikl,isn)/G1_dSV  ! Normalized Sphericity
          DENDRn= DENDRn -dDENDR *frac_j  ! New        Dendricity (+)
          SPHERn= SPHERn +dDENDR *frac_j  ! New        Sphericity

          OK__DE=max(zer0,                            ! IF 1., 
     .               sign(un_1, DENDRn                ! NO change
     .                         -eps6           ))     ! Dendr. -> Spheric

          G1__wd=OK__DE *    (      -DENDRn*G1_dSV)   ! Dendritic
     .      +(1.-OK__DE)* min(G1_dSV,SPHERn*G1_dSV)   ! Dendr. -> Spheric
          G2__wd=OK__DE * min(G1_dSV,SPHERn*G1_dSV)   ! Spheric
     .      +(1.-OK__DE)*(ADSdSV-min(SPHERn,vsphe1))  ! Spher. -> Size

!  1.2.2 Cas non dendritique non completement spherique
!        Evolution de la Sphericite seulement.
!        Non dendritic and not completely spheric Case
!        Evolution of    Sphericity only (not size)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          OK__ws=max(zer0,                            !
     .               sign(un_1, G1_dSV                !
     .                         -epsi5                 !
     .                         -G1snSV(ikl,isn)))     !

          SPHERn= G1snSV(ikl,isn)/G1_dSV
          SPHERn= SPHERn +dDENDR *frac_j
          G1__ws=         min(G1_dSV,SPHERn*G1_dSV)

!  1.2.3 Cas non dendritique et spherique / non dendritic and spheric
!        Evolution de la Taille seulement / Evolution of Size only
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          G2__ws =  husi_0
     .           *( husi_1
     .            *(husi_2 *( G2snSV(ikl,isn)/husi_0)**3
     .                      +(vtail1 +vtail2 *exp1Wa    )*dt__SV))
     .           ** husi_3


!  1.3 Metamorposes seches / Dry Metamorphism
!      --------------------------------------


!  1.3.1 Calcul Metamorphose faible/low Gradient (0.00-0.05 deg/cm)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          OKlowT=max(zer0,                            !
     .               sign(un_1, vgrat1                !
     .                         -dTsndz         ))     !

          facVap=exp(vvap1/TsisSV(ikl,isn))

!  1.3.1.1 Cas dendritique / dendritic Case

          OK_ldd=max(zer0,                            !
     .               sign(un_1,-G1snSV(ikl,isn)       !
     .                         -eps6           ))     !

          DENDRn=-G1snSV(ikl,isn)     /G1_dSV
          SPHERn= G2snSV(ikl,isn)     /G1_dSV
          DENDRn= DENDRn-vdent1*facVap*frac_j
          SPHERn= SPHERn+vsphe2*facVap*frac_j

          OK__DE=max(zer0,                            ! IF 1., 
     .               sign(un_1, DENDRn                ! NO change
     .                         -eps6           ))     ! Dendr. -> Spheric

          G1_ldd= OK__DE *    (      -DENDRn*G1_dSV)  ! Dendritic
     .       +(1.-OK__DE)* min(G1_dSV,SPHERn*G1_dSV)  ! Dendr. -> Spheric
          G2_ldd= OK__DE * min(G1_dSV,SPHERn*G1_dSV)  ! Spheric
     .       +(1.-OK__DE)*(ADSdSV-min(SPHERn,vsphe1)) ! Spher. -> Size 

!  1.3.1.2 Cas non dendritique / non dendritic Case

          SPHERn=G1snSV(ikl,isn)/G1_dSV
          DiamGx=G2snSV(ikl,isn)*0.1

          istoOK=min( abs(istoSV(ikl,isn)-
     .                    istdSV(1)      ),1)         ! zero if istoSV = 1
          DiamOK=max(zer0,  sign(un_1,vdiam2-DiamGx))
          No_Big=    istoOK+DiamOK
          No_Big=min(No_Big,un_1)

          dSPHER=           vsphe2*facVap*frac_j      !
          SPHER0=    SPHERn+dSPHER                    ! small grains
          SPHbig=    SPHERn+dSPHER                    ! big   grains
     .        *exp(min(zer0,vdiam3-G2snSV(ikl,isn)))  ! (history = 2 or 3)
          SPHbig=       min(vsphe3,SPHbig)            ! limited sphericity
          SPHERn= No_Big *  SPHER0
     .      + (1.-No_Big)*  SPHbig

          G1_lds=       min(G1_dSV,SPHERn*G1_dSV)

!  1.3.2 Calcul Metamorphose Gradient Moyen/Moderate (0.05-0.15)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          OK_mdT=max(zer0,                            !
     .               sign(un_1, vgrat2                !
     .                         -dTsndz))              !
          OKmidT=               OK_mdT  *(1.-OKlowT)  !
          OKhigT=          (1. -OK_mdT) *(1.-OKlowT)  !

          facVap=vdent1*exp(vvap1/TsisSV(ikl,isn))
     .                 *   (1.e2 *dTsndz)**vvap2

!  1.3.2.1 cas dendritique / dendritic case.

          OK_mdd=max(zer0,                            !
     .               sign(un_1,-G1snSV(ikl,isn)       !
     .                         -eps6           ))     !

          DENDRn=-G1snSV(ikl,isn)/G1_dSV
          SPHERn= G2snSV(ikl,isn)/G1_dSV
          DENDRn= DENDRn - facVap*frac_j
          SPHERn= SPHERn - facVap*frac_j

          OK__DE=max(zer0,                            ! IF 1., 
     .               sign(un_1, DENDRn                ! NO change
     .                         -eps6           ))     ! Dendr. -> Spheric

          G1_mdd= OK__DE *    (      -DENDRn*G1_dSV)  ! Dendritic
     .       +(1.-OK__DE)* max(zer0  ,SPHERn*G1_dSV)  ! Dendr. -> Spheric
          G2_mdd= OK__DE * max(zer0  ,SPHERn*G1_dSV)  ! Spheric
     .       +(1.-OK__DE)*(ADSdSV-max(SPHERn,zer0  )) ! Spher. -> Size 

!  1.3.2.2 Cas non dendritique / non dendritic Case

          SPHERn=G1snSV(ikl,isn)/G1_dSV
          SPHERn=         SPHERn-facVap*frac_j
          G1_mds=max(zer0,SPHERn*G1_dSV)

!  1.3.3 Calcul Metamorphose fort / high Gradient
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          facVap=vdent1*exp(vvap1/TsisSV(ikl,isn))
     .                 *   (1.e2 *dTsndz)**vvap2

!  1.3.3.1 Cas dendritique / dendritic Case

          OK_hdd=max(zer0,                            !
     .               sign(un_1,-G1snSV(ikl,isn)       !
     .                         -eps6           ))     !

          DENDRn=-G1snSV(ikl,isn)/G1_dSV              !
          SPHERn= G2snSV(ikl,isn)/G1_dSV              !
          DENDRn= DENDRn - facVap*frac_j              !
          SPHERn= SPHERn - facVap*frac_j              ! Non dendritic
                                                      ! and angular
          OK__DE=max(zer0,                            ! IF 1., 
     .               sign(un_1, DENDRn                ! NO change
     .                         -eps6  ))              ! Dendr. -> Spheric

          G1_hdd= OK__DE *    (      -DENDRn*G1_dSV)  ! Dendritic
     .       +(1.-OK__DE)* max(zer0  ,SPHERn*G1_dSV)  ! Dendr. -> Spheric
          G2_hdd= OK__DE * max(zer0  ,SPHERn*G1_dSV)  ! Spheric
     .       +(1.-OK__DE)*(ADSdSV-max(SPHERn,zer0  )) ! Spher. -> Size 

!  1.3.3.2 Cas non dendritique non completement anguleux.
!          non dendritic and spericity gt. 0

          OK_hds=max(zer0,                            !
     .               sign(un_1, G1snSV(ikl,isn)       !
     .                         -eps6           ))     !

          SPHERn= G1snSV(ikl,isn)/G1_dSV
          SPHERn= SPHERn - facVap*frac_j
          G1_hds= max(zer0,SPHERn*G1_dSV)

!  1.3.3.3 Cas non dendritique et anguleux
!          dendritic and spericity = 0.

          T1__OK = max(zer0,sign(un_1,TsisSV(ikl,isn)-Tf_Sno+vtang1))
          T2__OK = max(zer0,sign(un_1,TsisSV(ikl,isn)-Tf_Sno+vtang2))
          T3_xOK = max(zer0,sign(un_1,TsisSV(ikl,isn)-Tf_Sno+vtang3))
          T3__OK =                    T3_xOK  * (1. - T2__OK)
          T3_nOK =              (1. - T3_xOK) * (1. - T2__OK)
          ro1_OK = max(zer0,sign(un_1,vrang1-ro_dry(ikl,isn)))
          ro2_OK = max(zer0,sign(un_1,ro_dry(ikl,isn)-vrang2))
          dT1_OK = max(zer0,sign(un_1,vgang1-dTsndz         ))
          dT2_OK = max(zer0,sign(un_1,vgang2-dTsndz         ))
          dT3xOK = max(zer0,sign(un_1,vgang3-dTsndz         ))
          dT3_OK =                    dT3xOK  * (1. - dT2_OK)
          dT4xOK = max(zer0,sign(un_1,vgang4-dTsndz         ))
          dT4_OK =                    dT4xOK  * (1. - dT3_OK)
     .                                        * (1. - dT2_OK)
          dT4nOK =              (1. - dT4xOK) * (1. - dT3_OK)
     .                                        * (1. - dT2_OK)

!  Influence de la Temperature /Temperature Influence
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          AngSno = 
     .      T1__OK                                                       ! 11
     .    *(T2__OK*(vtang4+vtang5*(Tf_Sno       -TsisSV(ikl,isn))        ! 12
     .                    /vtang6)                                       !
     .     +T3__OK*(vtang7-vtang8*(Tf_Sno-vtang2-TsisSV(ikl,isn))        ! 13
     .                    /vtang9)                                       !
     .     +T3_nOK*(vtanga-vtangb*(Tf_Sno-vtang3-TsisSV(ikl,isn))        ! 14
     .                    /vtangc))                                      !

!  Influence de la Masse Volumique /Density Influence
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     .    * ro1_OK
     .        *(   ro2_OK*(1. - (ro_dry(ikl,isn)-vrang2)                 !   
     .                                  /(vrang1-vrang2))                !
     .         +1.-ro2_OK                                )               !   

!  Influence du Gradient de Temperature /Temperature Gradient Influence
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     .        *(   dT1_OK*(dT2_OK*vgang5*(dTsndz-vgang6)                 ! 15
     .                                  /(vgang2-vgang6)                 !
     .                    +dT3_OK*vgang7                                 ! 16
     .                    +dT4_OK*vgang9                                 ! 17
     .                    +dT4nOK*vgangb                )                ! 18
     .         +1.-dT1_OK                                )               !
     .    + ro1_OK
     .        *    dT1_OK*(dT3_OK*vgang8*(dTsndz-vgang2)
     .                                  /(vgang3-vgang2)
     .                    +dT4_OK*vganga*(dTsndz-vgang3)
     .                                  /(vgang4-vgang3)
     .                    +dT4nOK*vgangc*(dTsndz-vgang4)
     .                                  /(vgang1-vgang4))

          G2_hds = G2snSV(ikl,isn) + 1.d2 *AngSno*vfi     *frac_j


! New Properties
! --------------

          G1_bak          = G1snSV(ikl,isn)
          G2_bak          = G2snSV(ikl,isn)

          G1snSV(ikl,isn) = Wet_OK * (    OK__wd             *G1__wd    !  1
     .                               +(1.-OK__wd)*    OK__ws *G1__ws    !  2
     .                               +(1.-OK__wd)*(1.-OK__ws)*G1_bak)   !  3
     .               +(1. - Wet_OK)                                     !
     .                *(    OKlowT  *(    OK_ldd             *G1_ldd    !  4
     .                               +(1.-OK_ldd)            *G1_lds)   !  5
     .                    + OKmidT  *(    OK_mdd             *G1_mdd    !  6
     .                               +(1.-OK_mdd)            *G1_mds)   !  7
     .                    + OKhigT  *(    OK_hdd             *G1_hdd    !  8
     .                               +(1.-OK_hdd)*    OK_hds *G1_hds    !  9
     .                               +(1.-OK_hdd)*(1.-OK_hds)*G1_bak))  ! 10

          G2snSV(ikl,isn) = Wet_OK * (    OK__wd             *G2__wd    !  1
     .                               +(1.-OK__wd)*    OK__ws *G2_bak    !  2
     .                               +(1.-OK__wd)*(1.-OK__ws)*G2__ws)   !  3
     .               +(1. - Wet_OK)                                     !
     .                *(    OKlowT  *(    OK_ldd             *G2_ldd    !  4
     .                               +(1.-OK_ldd)            *G2_bak)   !  5
     .                    + OKmidT  *(    OK_mdd             *G2_mdd    !  6
     .                               +(1.-OK_mdd)            *G2_bak)   !  7
     .                    + OKhigT  *(    OK_hdd             *G2_hdd    !  8
     .                               +(1.-OK_hdd)*    OK_hds *G2_bak    !  9
     .                               +(1.-OK_hdd)*(1.-OK_hds)*G2_hds))  ! 10

! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp     G_curr( 1) =     Wet_OK             *    OK__wd             
! #vp     G_curr( 2) =     Wet_OK             *(1.-OK__wd)*    OK__ws 
! #vp     G_curr( 3) =     Wet_OK             *(1.-OK__wd)*(1.-OK__ws)
! #vp     G_curr( 4) = (1.-Wet_OK)*    OKlowT *    OK_ldd             
! #vp     G_curr( 5) = (1.-Wet_OK)*    OKlowT *(1.-OK_ldd)            
! #vp     G_curr( 6) = (1.-Wet_OK)*    OKmidT *    OK_mdd             
! #vp     G_curr( 7) = (1.-Wet_OK)*    OKmidT *(1.-OK_mdd)            
! #vp     G_curr( 8) = (1.-Wet_OK)*    OKhigT *    OK_hdd             
! #vp     G_curr( 9) = (1.-Wet_OK)*    OKhigT *(1.-OK_hdd)*    OK_hds 
! #vp     G_curr(10) = (1.-Wet_OK)*    OKhigT *(1.-OK_hdd)*(1.-OK_hds)
! #vp     G_curr(11) =     T1__OK                         * G_curr(10)
! #vp     G_curr(12) =     T2__OK                         * G_curr(10)
! #vp     G_curr(13) =     T3__OK                         * G_curr(10)
! #vp     G_curr(14) =     T3_nOK                         * G_curr(10)
! #vp     G_curr(15) =     ro1_OK*     dT1_OK *    dT2_OK * G_curr(10)
! #vp     G_curr(16) =     ro1_OK*     dT1_OK *    dT3_OK * G_curr(10)
! #vp     G_curr(17) =     ro1_OK*     dT1_OK *    dT4_OK * G_curr(10)
! #vp     G_curr(18) =     ro1_OK*     dT1_OK *    dT4nOK * G_curr(10)

! #vp     Gcases( 1) = max(Gcases( 1),G_curr( 1))
! #vp     Gcases( 2) = max(Gcases( 2),G_curr( 2))
! #vp     Gcases( 3) = max(Gcases( 3),G_curr( 3))
! #vp     Gcases( 4) = max(Gcases( 4),G_curr( 4))
! #vp     Gcases( 5) = max(Gcases( 5),G_curr( 5))
! #vp     Gcases( 6) = max(Gcases( 6),G_curr( 6))
! #vp     Gcases( 7) = max(Gcases( 7),G_curr( 7))
! #vp     Gcases( 8) = max(Gcases( 8),G_curr( 8))
! #vp     Gcases( 9) = max(Gcases( 9),G_curr( 9))
! #vp     Gcases(10) = max(Gcases(10),G_curr(10))
! #vp     Gcases(11) = max(Gcases(11),G_curr(11))
! #vp     Gcases(12) = max(Gcases(12),G_curr(12))
! #vp     Gcases(13) = max(Gcases(13),G_curr(13))
! #vp     Gcases(14) = max(Gcases(14),G_curr(14))
! #vp     Gcases(15) = max(Gcases(15),G_curr(15))
! #vp     Gcases(16) = max(Gcases(16),G_curr(16))
! #vp     Gcases(17) = max(Gcases(17),G_curr(17))
! #vp     Gcases(18) = max(Gcases(18),G_curr(18))

! #vp     IF          (isn    .le.     isnoSV(ikl))
! #vp.    write(47,471)isn            ,isnoSV(ikl)                    ,
! #vp.                 TsisSV(ikl,isn),ro__SV(ikl,isn),eta_SV(ikl,isn),
! #vp.                 G1_bak         ,G2_bak         ,istoSV(ikl,isn),
! #vp.                 dTsndz,
! #vp.                (       k ,k=1,18),
! #vp.                (G_curr(k),k=1,18),
! #vp.                (Gcases(k),k=1,18),
! #vp.                 Wet_OK,OK__wd,G1__wd,G2__wd,
! #vp.                     1.-OK__wd,OK__ws,G1__ws,1.-OK__ws,G2__ws,
! #vp.              1.-Wet_OK,OKlowT,OK_ldd,G1_ldd,          G2_ldd,
! #vp.                            1.-OK_ldd,G1_lds,
! #vp.                        OKmidT,OK_mdd,G1_mdd,          G1_mdd,
! #vp.                            1.-OK_mdd,G1_mds,
! #vp.                        OKhigT,OK_hdd,G1_hdd,          G2_hdd,
! #vp.                            1.-OK_hdd,OK_hds,          G1_hds,
! #vp.                                             1.-OK_hds,G2_hds,
! #vp.                 G1snSV(ikl,isn),
! #vp.                 G2snSV(ikl,isn)
 471      format(
     .         /,' isn     =  ',i4,6x,'(MAX.:',i4,')',  
     .         /,' T       =  ',f8.3,
     .         /,' ro      =  ',f8.3,
     .         /,' eta     =  ',f8.3,
     .         /,' G1      =  ',f8.3,
     .         /,' G2      =  ',f8.3,
     .         /,' Histor. =  ',i4  ,
     .         /,' Grad(T) =  ',f8.3,'                   ' ,18i3  ,
     ./,         '                       Current    Case: ',18f3.0,
     ./,         '                       Cases performed: ',18f3.0,
     ./,' ------------------------------------------------------------',
     .             '-----------+------------------+------------------+',
     ./,' Status                                                      ',
     .             '           | G1               | G2               |',
     ./,' ------------------------------------------------------------',
     .             '-----------+------------------+------------------+',
     ./,'    Wet_OK: ',f8.3,'                  OK__wd: ',f8.3,'       ',
     .             '           | G1__wd: ',f8.3,' | G2__wd: ',f8.5,' |',
     ./,'                                   1.-OK__wd: ',f8.3,' OK__ws',
     .             ': ',f8.3,' | G1__ws: ',f8.3,' |                  |',
     ./,'                                                    1.-OK__ws',
     .             ': ',f8.3,' |                  | G2__ws: ',f8.5,' |',
     ./,' 1.-Wet_OK: ',f8.3,' OKlowT: ',f8.3,' OK_ldd: ',f8.3,'       ',
     .             '           | G1_ldd: ',f8.3,' | G2_ldd: ',f8.5,' |',
     ./,'                                   1.-OK_ldd: ',f8.3,'       ',
     .             '           | G1_lds: ',f8.3,' |                  |',
     ./,'                     OKmidT: ',f8.3,' OK_mdd: ',f8.3,'       ',
     .             '           | G1_mdd: ',f8.3,' | G2_mdd: ',f8.5,' |',
     ./,'                                   1.-OK_mdd: ',f8.3,'       ',
     .             '           | G1_mds: ',f8.3,' |                  |',
     ./,'                     OKhigT: ',f8.3,' OK_hdd: ',f8.3,'       ',
     .             '           | G1_hdd: ',f8.3,' | G2_hdd: ',f8.5,' |',
     ./,'                                   1.-OK_hdd: ',f8.3,' OK_hds',
     .             ': ',f8.3,' | G1_hds: ',f8.3,' |                  |',
     ./,'                                                    1.-OK_hds',
     .             ': ',f8.3,' |                  | G2_hds: ',f8.5,' |',
     ./,' ------------------------------------------------------------',
     .             '-----------+------------------+------------------+',
     ./,'                                                             ',
     .             '           |         ',f8.3,' |         ',f8.5,' |',
     ./,' ------------------------------------------------------------',
     .             '-----------+------------------+------------------+')
        END DO
      END DO


!  2. Mise a Jour Variables Historiques (Cas non dendritique)
!     Update of the historical Variables
!     =======================================================

      IF (vector)                                                   THEN
        DO isn=1,nsno
        DO ikl=1,knonv

          SphrOK = max(zer0,sign(un_1,       G1snSV(ikl,isn)))
          H1a_OK = max(zer0,sign(un_1,vsphe4-G1snSV(ikl,isn)))
          H1b_OK =     1   - min(1   ,       istoSV(ikl,isn))
          H1__OK =                    H1a_OK*H1b_OK
          H23aOK = max(zer0,sign(un_1,vsphe4-G1_dSV
     .                                      +G1snSV(ikl,isn)))
          H23bOK = max(zer0,sign(un_1,etaSno(ikl,isn)
     .                      /max(eps6,dzsnSV(ikl,isn))
     .                                      -vtelv1         ))
          H23_OK =                    H23aOK*H23bOK
          H2__OK =     1   - min(1   ,       istoSV(ikl,isn))
          H3__OK =     1   - min(1   ,   abs(istoSV(ikl,isn)-istdSV(1)))
          H45_OK = max(zer0,sign(un_1,Tf_Sno-TsisSV(ikl,isn)+eps6))
          H4__OK =     1   - min(1   ,   abs(istoSV(ikl,isn)-istdSV(2)))
          H5__OK =     1   - min(1   ,   abs(istoSV(ikl,isn)-istdSV(3)))

          HISupd          =
     .    SphrOK*(H1__OK                             *istdSV(1)
     .       +(1.-H1__OK)*    H23_OK         *(H2__OK*istdSV(2)
     .                                        +H3__OK*istdSV(3))
     .       +(1.-H1__OK)*(1.-H23_OK) *H45_OK*(H4__OK*istdSV(4)
     .                                        +H5__OK*istdSV(5)))
          istoSV(ikl,isn) =   HISupd  +
     .           (1.-min(un_1,HISupd))               *istoSV(ikl,isn)
        END DO
        END DO
      ELSE 


!  2. Mise a Jour Variables Historiques (Cas non dendritique)
!     Update of the historical Variables
!     =======================================================

        DO ikl=1,knonv
        DO isn=iiceSV(ikl),isnoSV(ikl)
          IF  (G1snSV(ikl,isn).ge.0.)                               THEN
            IF(G1snSV(ikl,isn).lt.vsphe4.and.istoSV(ikl,isn).eq.0)  THEN
                   istoSV(ikl,isn)=istdSV(1)
            ELSEIF(G1_dSV-G1snSV(ikl,isn)         .lt.vsphe4.and.
     .             etaSno(ikl,isn)/dzsnSV(ikl,isn).gt.vtelv1)       THEN
              IF  (istoSV(ikl,isn).eq.0)
     .             istoSV(ikl,isn)=   istdSV(2)
              IF  (istoSV(ikl,isn).eq.istdSV(1))  
     .             istoSV(ikl,isn)=   istdSV(3)
            ELSEIF(TsisSV(ikl,isn).lt.Tf_Sno)                       THEN
              IF  (istoSV(ikl,isn).eq.istdSV(2))  
     .             istoSV(ikl,isn)=   istdSV(4)
              IF  (istoSV(ikl,isn).eq.istdSV(3))  
     .             istoSV(ikl,isn)=   istdSV(5)
            END IF
          END IF
        END DO
        END DO
      END IF


!  3. Tassement mecanique /mechanical Settlement
!     ==========================================

        DO ikl=1,knonv
          SnMass(ikl) = 0.
        END DO
      DO isn=nsno,1,-1
        DO ikl=1,knonv
          dSnMas     = 100.*dzsnSV(ikl,isn)*ro_dry(ikl,isn)
          SnMass(ikl)=      SnMass(ikl)+0.5*dSnMas
          ViscSn     =      vvisc1         *vvisc2
     .               *exp(vvisc3           *ro_dry(ikl,isn)
     .                   +vvisc4*abs(Tf_Sno-TsisSV(ikl,isn)))
     .                                     *ro_dry(ikl,isn)/rovisc

!  Changement de Viscosite si Teneur en Eau liquide
!  Change of the Viscosity if liquid Water Content
!  ------------------------------------------------

          OK_Liq     =    max(zer0,sign(un_1,etaSno(ikl,isn)-eps6))
          OK_Ang     =    max(zer0,sign(un_1,vgran6-G1snSV(ikl,isn)))
     .                *(1-min(1   , abs(istoSV(ikl,isn)-istdSV(1))))

! OUTPUT/Verification: Snow   Properties
! #wp     IF (G1snSV(ikl,isn).gt.0..AND.G1snSV(ikl,isn).lt.vsphe4
! #wp.                             .AND.istoSV(ikl,isn).eq.     0)
! #wp.    THEN
! #wp       write(6,*) ikl,isn,' G1,G2,hist,OK_Ang  ',
! #wp.          G1snSV(ikl,isn), G2snSV(ikl,isn),istoSV(ikl,isn),OK_Ang
! #wp       stop "Grains anguleux mal d?finis"
! #wp     END IF

          OKxLiq     =    max(zer0,sign(un_1,vtelv1-etaSno(ikl,isn)
     .                                    /max(eps6,dzsnSV(ikl,isn))))
     .               *    max(0   ,sign(1   ,istoSV(ikl,isn)
     .                                      -istdSV(1)      ))
          ViscSn     =  
     .    ViscSn*(    OK_Liq/(vvisc5+vvisc6*etaSno(ikl,isn)
     .                            /max(eps6,dzsnSV(ikl,isn)))
     .           +(1.-OK_Liq)                               )
     .          *(    OK_Ang*exp(min(ADSdSV,G2snSV(ikl,isn)-vdiam4))
     .           +(1.-OK_Ang)                                       )
     .          *(    OKxLiq*        vvisc7
     .           +(1.-OKxLiq)              )


!  Calcul nouvelle Epaisseur / new Thickness
!  -----------------------------------------

          dzsnew         =
     .    dzsnSV(ikl,isn) 
     .     *max(vdz3,
     .         (un_1-dt__SV*max(SnMass(ikl)*cos(slopSV(ikl)),un_1)
     .                     /max(ViscSn                      ,eps6)))
          rosnew         = ro__SV(ikl,isn) *dzsnSV(ikl,isn)  
     .                            /max(eps6,dzsnew)
          rosmax         = 1.d0   /( (1.d0 -eta_SV(ikl,isn)) /rhoIce
     .                               +      eta_SV(ikl,isn)  /rhoWat)
          rosnew         =                        min(rosnew ,rosmax)
          dzsnSV(ikl,isn)= dzsnSV(ikl,isn) *ro__SV(ikl,isn)
     .                            /max(eps6,rosnew)
          ro__SV(ikl,isn)= rosnew
          ro_dry(ikl,isn)= ro__SV(ikl,isn)*(1.-eta_SV(ikl,isn))*1.e-3
!         ro_dry: Dry Density (g/cm3)

          SnMass(ikl)    = SnMass(ikl)+dSnMas*0.5
        END DO
      END DO


! OUTPUT/Verification: Snow   Properties
! #wp DO ikl = 1,knonv
! #wp DO isn = 1,isnoSV(ikl)
! #wp   IF (G1snSV(ikl,isn).gt.0. .AND. G2snSV(ikl,isn).gt.D__MAX) THEN
! #wp     write(6,6600) G1snSV(ikl,isn),G2snSV(ikl,isn),ikl,isn
 6600     format(/,'WARNING in _GSn: G1,G2 =',2f9.3,'  (ikl,isn) =',2i4)
! #wp     D__MAX =                      G2snSV(ikl,isn)
! #wp   END IF
! #wp   IF (                            G2snSV(ikl,isn).lt.0.    ) THEN
! #wp     write(6,6601) G1snSV(ikl,isn),G2snSV(ikl,isn),ikl,isn
 6601     format(/,'ERROR 1 in _GSn: G1,G2 =',2f9.3,'  (ikl,isn) =',2i4)
! #wp     STOP
! #wp   END IF
! #wp   IF (G1snSV(ikl,isn).gt.G1_dSV+eps6                       ) THEN
! #wp     write(6,6602) G1snSV(ikl,isn),G2snSV(ikl,isn),ikl,isn
 6602     format(/,'ERROR 2 in _GSn: G1,G2 =',2f9.3,'  (ikl,isn) =',2i4)
! #wp     STOP
! #wp   END IF
! #wp   IF (G1snSV(ikl,isn).lt.0.                           .AND.
! #wp.      G2snSV(ikl,isn).gt.G1_dSV+eps6                       ) THEN
! #wp     write(6,6603) G1snSV(ikl,isn),G2snSV(ikl,isn),ikl,isn
 6603     format(/,'ERROR 3 in _GSn: G1,G2 =',2f9.3,'  (ikl,isn) =',2i4)
! #wp     STOP
! #wp   END IF
! #wp END DO
! #wp END DO

      return
      end


      subroutine SISVAT_qSo
! #m0.                     (Wats_0,Wats_1,Wats_d)

!--------------------------------------------------------------------------+
!   MAR          SISVAT_qSo                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_qSo computes the Soil      Water  Balance          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  klonv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   isnoSV   = total Nb of Ice/Snow Layers                      |
!     ^^^^^    isotSV   = 0,...,11:   Soil       Type                      |
!                         0:          Water, Solid or Liquid               |
!                                                                          |
!     INPUT:   rhT_SV   : SBL Top    Air  Density                  [kg/m3] |
!     ^^^^^    drr_SV   : Rain   Intensity                       [kg/m2/s] |
!              LSdzsv   : Vertical   Discretization Factor             [-] |
!                       =    1. Soil                                       |
!                       = 1000. Ocean                                      |
!              dt__SV   : Time   Step                                  [s] |
!                                                                          |
!              Lx_H2O   : Latent Heat of Vaporization/Sublimation   [J/kg] |
!              HLs_sv   : Latent Heat  Flux                         [W/m2] |
!              Rootsv   : Root   Water Pump                      [kg/m2/s] |
!                                                                          |
!     INPUT /  eta_SV   : Water      Content                       [m3/m3] |
!     OUTPUT:  Khydsv   : Soil   Hydraulic    Conductivity           [m/s] |
!     ^^^^^^                                                               |
!                                                                          |
!     OUTPUT:  RnofSV   : RunOFF Intensity                       [kg/m2/s] |
!     ^^^^^^   Wats_0   : Soil Water,  before Forcing                 [mm] |
!              Wats_1   : Soil Water,  after  Forcing                 [mm] |
!              Wats_d   : Soil Water          Forcing                 [mm] |
!                                                                          |
!     Internal Variables:                                                  |
!     ^^^^^^^^^^^^^^^^^^                                                   |
!              z_Bump   : (Partly)Bumpy Layers Height                  [m] |
!              z0Bump   :         Bumpy Layers Height                  [m] |
!              dzBump   :  Lowest Bumpy Layer:                         [m] |
!              etBump   :         Bumps Layer Averaged Humidity    [m3/m3] |
!              etaMid   : Layer Interface's Humidity               [m3/m3] |
!              eta__f   : Layer             Humidity  (Water Front)[m3/m3] |
!              Dhyd_f   : Soil  Hydraulic Diffusivity (Water Front) [m2/s] |
!              Dhydif   : Soil  Hydraulic Diffusivity               [m2/s] |
!              WgFlow   : Water         gravitational     Flux   [kg/m2/s] |
!              Wg_MAX   : Water MAXIMUM gravitational     Flux   [kg/m2/s] |
!              SatRat   : Water         Saturation        Flux   [kg/m2/s] |
!              WExces   : Water         Saturation Excess Flux   [kg/m2/s] |
!              Dhydtz   : Dhydif * dt / dz                             [m] |
!              FreeDr   : Free Drainage Fraction                       [-] |
!              Elem_A   : A Diagonal Coefficient                           |
!              Elem_C   : C Diagonal Coefficient                           |
!              Diag_A   : A Diagonal                                       |
!              Diag_B   : B Diagonal                                       |
!              Diag_C   : C Diagonal                                       |
!              Term_D   :   Independant Term                               |
!              Aux__P   : P Auxiliary Variable                             |
!              Aux__Q   : Q Auxiliary Variable                             |
!                                                                          |
!     TUNING PARAMETER:                                                    |
!     ^^^^^^^^^^^^^^^^                                                     |
!              z0soil   : Soil Surface averaged Bumps Height           [m] |
!                                                                          |
!     METHOD: NO   Skin Surface Humidity                                   |
!     ^^^^^^  Semi-Implicit Crank Nicholson Scheme                         |
!             (Partial) free Drainage, Water Bodies excepted (Lakes, Sea)  |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #GF: Saturation Front                                                |
!     #GH: Saturation Front allows Horton Runoff                           |
!     #GA: Soil Humidity Geometric Average                                 |
!     #BP: Parameterization of Terrain Bumps                               |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_iii_jjj_n     | #m0: OUTPUT/Verification: H2O    Conservation |
!   # SISVAT_iii_jjj_n     | #m1: OUTPUT/Verification: * Mass Conservation |
!   # stdout               | #mw: OUTPUT/Verification: H2O    Conservation |
!                          |      unit  6, SubRoutine  SISVAT_qSo **ONLY** |
!   # SISVAT_qSo.vw        | #vw: OUTPUT/Verif+Detail: H2O    Conservation |
!                          |      unit 42, SubRoutine  SISVAT_qSo **ONLY** |
!   # stdout               | #sg: OUTPUT/Verification: Gravitational Front |
!                          |      unit  6, SubRoutine  SISVAT_qSo **ONLY** |
!                                                                          |
!     REMARQUE: Inclure possibilite de creer mare sur bedrock impermeable  |
!     ^^^^^^^^                                                             |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE  PHY_SV

      USE VAR_SV
      USE VARdSV
      USE VAR0SV

      USE VARxSV
      USE VARySV


      IMPLICIT NONE

! OUTPUT
! ------

! OUTPUT/Verification: H2O    Conservation
! #m0 real      Wats_0(klonv)                 ! Soil Water,  before forcing
! #m0 real      Wats_1(klonv)                 ! Soil Water,  after  forcing
! #m0 real      Wats_d(klonv)                 ! Soil Water          forcing


! Internal Variables
! ==================

      integer   isl   ,jsl   ,ist   ,ikl      !
      integer   ikm   ,ikp   ,ik0   ,ik1      !
      integer   ist__s,ist__w                 ! Soil/Water Body Identifier
c #BP real      z0soil                        ! Soil Surface Bumps Height  [m]
c #BP real      z_Bump                        !(Partly)Bumpy Layers Height [m] 
c #BP real      z0Bump                        !        Bumpy Layers Height [m] 
c #BP real      dzBump                        ! Lowest Bumpy Layer: 

c #BP real      etBump(klonv)                 ! Bumps Layer Averaged Humidity
      real      etaMid                        ! Layer Interface's Humidity
      real      Dhydif                        ! Hydraulic Diffusivity   [m2/s]
      real      eta__f                        ! Water Front Soil Water Content
      real      Khyd_f                        ! Water Front Hydraulic Conduct.
      real      Khydav                        ! Hydraulic Conductivity   [m/s]
      real      WgFlow                        ! Water gravitat. Flux [kg/m2/s]
      real      Wg_MAX                        ! Water MAX.grav. Flux [kg/m2/s]
      real      SatRat                        ! Saturation      Flux [kg/m2/s]
      real      WExces                        ! Saturat. Excess Flux [kg/m2/s]
      real      SoRnOF(klonv)                 ! Soil     Run    OFF
      real      Dhydtz(klonv,-nsol:0)         ! Dhydif * dt / dz           [m]
      real      Elem_A,Elem_B,Elem_C          !   Diagonal Coefficients
      real      Diag_A(klonv,-nsol:0)         ! A Diagonal
      real      Diag_B(klonv,-nsol:0)         ! B Diagonal
      real      Diag_C(klonv,-nsol:0)         ! C Diagonal
      real      Term_D(klonv,-nsol:0)         !   Independant Term
      real      Aux__P(klonv,-nsol:0)         ! P Auxiliary Variable
      real      Aux__Q(klonv,-nsol:0)         ! Q Auxiliary Variable
      real      etaaux(klonv,-nsol:-nsol+1)   ! Soil Water Content     [m3/m3]
      real      FreeDr                        ! Free Drainage Fraction (actual)
      real      FreeD0                        ! Free Drainage Fraction (1=Full)

! OUTPUT/Verification: H2O    Conservation
! #mw logical         mwopen                  ! IO   Switch
! #mw common/Sm_qSo_L/mwopen                  !
! #mw real     hourwr,timewr                  !
! #mw common/Sm_qSo_R/timewr                  !
! #mw real            Evapor(klonv)           !


! Internal DATA
! =============

c #BP data      z0soil/0.020/                 ! Soil Surface Bumps Height  [m]
      data      FreeD0/1.000/                 ! Free Drainage Fraction (1=Full)


! OUTPUT/Verification: H2O    Conservation: Water  Budget (IN)
! #m0   DO ikl=1,knonv
! #m0     Wats_0(ikl) = 0.                    ! OLD RunOFF Contrib.
! #m0     Wats_d(ikl) = drr_SV(ikl)           ! Water Surface Forc.
! #m0   END DO

! #m0      isl= -nsol
! #m0   DO ikl=1,knonv
! #m0     Wats_0(ikl) = Wats_0(ikl) 
! #m0.      + rhoWat *( eta_SV(ikl,isl)   *dz78SV(isl)
! #m0.                + eta_SV(ikl,isl+1) *dz_8SV(isl) ) * LSdzsv(ikl)
! #m0   END DO

! #m0 DO   isl= -nsol+1,-1
! #m0   DO ikl=1,knonv
! #m0     Wats_0(ikl) = Wats_0(ikl) 
! #m0.      + rhoWat *( eta_SV(ikl,isl)   *dz34SV(isl)
! #m0.                +(eta_SV(ikl,isl-1)
! #m0.                 +eta_SV(ikl,isl+1))*dz_8SV(isl) ) * LSdzsv(ikl)
! #m0   END DO
! #m0 END DO 

! #m0      isl=  0
! #m0   DO ikl=1,knonv
! #m0     Wats_0(ikl) = Wats_0(ikl) 
! #m0.      + rhoWat *( eta_SV(ikl,isl)   *dz78SV(isl)
! #m0.                + eta_SV(ikl,isl-1) *dz_8SV(isl) ) * LSdzsv(ikl)
! #m0   END DO


! Gravitational Flow
! ==================

! .    METHOD: Surface Water Flux saturates successively the soil layers 
!      ^^^^^^  from up to below, but is limited by infiltration capacity.
!              Hydraulic Conductivity again contributes after this step,
!              not redundantly because of a constant (saturated) profile.

! Flux  Limitor
! ^^^^^^^^^^^^^
           isl=0
        DO ikl=1,knonv
          ist    = isotSV(ikl)                     ! Soil Type
          ist__s = min(ist, 1)                     ! 1 => Soil
          ist__w = 1 - ist__s                      ! 1 => Water Body
          Dhydif = s1__SV(ist)
     .               *max(eps6,eta_SV(ikl,isl))    ! Hydraulic Diffusivity
     .                      **(bCHdSV(ist)+2.)     ! DR97, Eqn.(3.36)
          Dhydif = ist__s    * Dhydif              !
     .           + ist__w    * vK_dSV              ! Water Bodies

          Khydav = ist__s    * Ks_dSV(ist)         ! DR97  Assumption
     .           + ist__w    * vK_dSV              ! Water Bodies

          Wg_MAX = rhoWat     *Dhydif              ! MAXimum  Infiltration
     .           *(etadSV(ist)-eta_SV(ikl,isl))    !          Rate
     .           /(dzAvSV(isl)*LSdzsv(ikl)    )    !
     .          +  rhoWat     *Khydav              ! 

! Surface Horton RunOFF
! ^^^^^^^^^^^^^^^^^^^^^
          SoRnOF(ikl) =
     .                max(zer0,drr_SV(ikl)-Wg_MAX)
          drr_SV(ikl) =        drr_SV(ikl)-SoRnOF(ikl)
        END DO

c #GF DO   isl=0,-nsol,-1
c #GF   DO ikl=1,knonv
c #GF     ist    = isotSV(ikl)                     ! Soil Type
c #GF     ist__s = min(ist, 1)                     ! 1 => Soil
c #GF     ist__w = 1 - ist__s                      ! 1 => Water Body

! Water Diffusion
! ^^^^^^^^^^^^^^^
c #GF     Dhydif = s1__SV(ist)
c #GF.               *max(eps6,eta_SV(ikl,isl))    ! Hydraulic Diffusivity
c #GF.                      **(bCHdSV(ist)+2.)     ! DR97, Eqn.(3.36)
c #GF     Dhydif = ist__s    * Dhydif              !
c #GF.           + ist__w    * vK_dSV              ! Water Bodies

! Water Conduction (without Horton Runoff)
! ^^^^^^^^^^^^^^^^
c #GF     Khyd_f =             Ks_dSV(ist)
!         Uses saturated K ==> Horton Runoff ~0    !

! Water Conduction (with    Horton Runoff)
! ^^^^^^^^^^^^^^^^
c #GH     ik0    = nkhy       *eta_SV(ikl,isl)
c #GH.                        /etadSV(ist)
c #GH     eta__f         =            1.
c #GH.   -aKdtSV(ist,ik0)/(2. *dzAvSV(isl)
c #GH.                        *LSdzsv(ikl))
c #GH     eta__f         = max(eps_21,eta__f)
c #GH     eta__f         = min(etadSV(ist),
c #GH.                         eta_SV(ikl,isl) +
c #GH.   (aKdtSV(ist,ik0)     *eta_SV(ikl,isl)
c #GH.   +bKdtSV(ist,ik0))   /(dzAvSV(isl)
c #GH.                        *LSdzsv(ikl))
c #GH.                       / eta__f          )

c #GH     eta__f         = .5*(eta_SV(ikl,isl)
c #GH.                        +eta__f)
!         eta__f         =     eta_SV(ikl,isl)     ! Another Possibility

c #GH     ik0    = nkhy       *eta__f
c #GH.                        /etadSV(ist)
c #GH     Khyd_f =
c #GH.   (aKdtSV(ist,ik0)     *eta__f
c #GH.   +bKdtSV(ist,ik0))    /dt__SV

c #GF     Khydav = ist__s    * Khyd_f              ! DR97  Assumption
c #GF.           + ist__w    * vK_dSV              ! Water Bodies

! Gravitational Flow
! ^^^^^^^^^^^^^^^^^^
c #GF     Wg_MAX =                                 ! MAXimum  Infiltration
c #GF.             rhoWat     *Dhydif              !          Rate
c #GF.           *(etadSV(ist)-eta_SV(ikl,isl))    ! 
c #GF.           /(dzAvSV(isl)*LSdzsv(ikl)    )    !
c #GF.          +  rhoWat     *Khydav              ! 

! OUTPUT/Verification: Gravitational Front
! #sg     write(6,6001) isl,drr_SV(ikl)*3.6e3,Wg_MAX     *3.6e3
 6001     format(i3,'  vRain ,Wg_MAX ',2e12.3,' mm/hr')

c #GF     WgFlow =  min(Wg_MAX,drr_SV(ikl))        ! Infiltration Rate
c #GF     WExces =  max(zer0  ,drr_SV(ikl)-WgFlow) ! Water Excess => RunOff

! OUTPUT/Verification: Gravitational Front
! #sg     write(6,6002)     WgFlow     *3.6e3,WExces     *3.6e3
 6002     format(3x,'  WgFlow,WExces ',2e12.3,' mm/hr')

c #GF     SoRnOF(ikl) =        SoRnOF(ikl)+WExces  !
c #GF     drr_SV(ikl) =        WgFlow              !

! OUTPUT/Verification: Gravitational Front
! #sg     write(6,6003)     SoRnOF(ikl)*3.6e3,drr_SV(ikl)*3.6e3
 6003     format(3x,'  SoRnOF,drr_SV ',2e12.3,' mm/hr')

c #GF     SatRat =(etadSV(ist)-eta_SV(ikl,isl))    ! Saturation   Rate
c #GF.            *rhoWat     *dzAvSV(isl)         !
c #GF.                        *LSdzsv(ikl)/dt__SV  !
c #GF     SatRat =  min(SatRat,drr_SV(ikl))        !
c #GF     drr_SV(ikl) =        drr_SV(ikl)-SatRat  ! Water Flux for Below

! OUTPUT/Verification: Gravitational Front
! #sg     write(6,6004)     SatRat     *3.6e3,drr_SV(ikl)*3.6e3
 6004     format(3x,'  SatRat,drr_SV ',2e12.3,' mm/hr')
! #sg     write(6,6005)     eta_SV(ikl,isl)*1.e3

c #GF     eta_SV(ikl,isl) =    eta_SV(ikl,isl)     ! Saturation
c #GF.                        +SatRat*dt__SV       !
c #GF.                       /(rhoWat*dzAvSV(isl)  !
c #GF.                               *LSdzsv(ikl)) !

! OUTPUT/Verification: Gravitational Front
! #sg     write(6,6005)     eta_SV(ikl,isl)*1.e3

 6005     format(3x,'  eta_SV,       ',e12.3,' g/kg')
c #GF   END DO
c #GF END DO
c #GF   DO ikl=1,knonv
c #GF     SoRnOF(ikl)     =    SoRnOF(ikl)         ! RunOFF Intensity
c #GF.                    +    drr_SV(ikl)         ! [kg/m2/s]
!         Inclure la possibilite de creer une mare sur un bedrock impermeable
c #GF     drr_SV(ikl) = 0.
c #GF   END DO


! Temperature Correction due to a changed Soil Energy Content
! ===========================================================

! REMARQUE: Mettre en oeuvre le couplage humidite-energie
! ^^^^^^^^


! Full Resolution of the Richard's Equation
! =========================================

!         METHOD: Water content evolution results from water fluxes 
!         ^^^^^^  at the layer boundaries
!                 Conductivity is approximated by a piecewise linear profile. 
!                 Semi-Implicit Crank-Nicholson scheme is used. 
!                (Bruen, 1997, Sensitivity of hydrological processes 
!                              at the land-atmosphere interface. 
!                              Proc. Royal Irish Academy,  IGBP symposium 
!                              on global change and the Irish Environment. 
!                              Publ.: Maynooth)

!                        - - - - - - - -   isl+1/2   - -  ^   
!                                                         |                       
!     eta_SV(isl)        ---------------   isl     -----  +--dz_dSV(isl)  ^
!                                                         |               |
!     Dhydtz(isl) etaMid - - - - - - - -   isl-1/2   - -  v  dzmiSV(isl)--+
!                                                                         |
!     eta_SV(isl-1)      ---------------   isl-1   -----                  v

! Transfert       Coefficients
! ----------------------------

      DO   isl=-nsol+1,0
        DO ikl=1,knonv
          ist    =      isotSV(ikl)                       ! Soil Type
          ist__s =      min(ist, 1)                       ! 1 => Soil
          ist__w =      1 - ist__s                        ! 1 => Water Body
          etaMid =     (dz_dSV(isl)  *eta_SV(ikl,isl-1)   ! eta at layers
     .                 +dz_dSV(isl-1)*eta_SV(ikl,isl)  )  !     interface
     .           /(2.0* dzmiSV(isl))                      ! LSdzsv implicit !
c #GA     etaMid = sqrt(dz_dSV(isl)  *eta_SV(ikl,isl-1)   ! Idem, geometric
c #GA.                 *dz_dSV(isl-1)*eta_SV(ikl,isl)  )  !       average
c #GA.           /(2.0* dzmiSV(isl))                      ! (Vauclin&al.1979)
          Dhydif          =    s1__SV(ist)                ! Hydraul.Diffusi.
     .  *(etaMid         **(   bCHdSV(ist)+2.))           ! DR97, Eqn.(3.36)
          Dhydtz(ikl,isl) =    Dhydif*dt__SV              !
     .                              /(dzmiSV(isl)         !
     .                               *LSdzsv(ikl))        !
          Dhydtz(ikl,isl) =    Dhydtz(ikl,isl) * ist__s   ! Soil
     .        +0.5*dzmiSV(isl)*LSdzsv(ikl)     * ist__w   ! Water bodies

        END DO
      END DO
           isl=-nsol
        DO ikl=1,knonv
          Dhydtz(ikl,isl) =    0.0                        !
        END DO


! Tridiagonal Elimination: Set Up
! -------------------------------

! Soil/Snow Interior
! ^^^^^^^^^^^^^^^^^^
      DO   isl=-nsol,-nsol+1
        DO ikl=1,knonv
          etaaux(ikl,isl) =  eta_SV(ikl,isl)
        END DO
      END DO

      DO   isl=-nsol+1,-1
        DO ikl=1,knonv
          ist      =         isotSV(ikl)
          ikm      = nkhy *  eta_SV(ikl,isl-1) / etadSV(ist)
          ik0      = nkhy *  eta_SV(ikl,isl)   / etadSV(ist)
          ikp      = nkhy *  eta_SV(ikl,isl+1) / etadSV(ist)
          Elem_A   =         Dhydtz(ikl,isl)   
     .                    -  aKdtSV(ist,ikm)* dziiSV(isl)  *LSdzsv(ikl)
          Elem_B   =      - (Dhydtz(ikl,isl)
     .                      +Dhydtz(ikl,isl+1)
     .                      -aKdtSV(ist,ik0)*(dziiSV(isl+1)
     .                                       -dzi_SV(isl) )*LSdzsv(ikl))
          Elem_C   =         Dhydtz(ikl,isl+1) 
     .                    +  aKdtSV(ist,ikp)* dzi_SV(isl+1)*LSdzsv(ikl)
          Diag_A(ikl,isl) =  dz_8SV(isl)        *LSdzsv(ikl)
     .                      -Implic            * Elem_A
          Diag_B(ikl,isl) =  dz34SV(isl)        *LSdzsv(ikl)
     .                      -Implic            * Elem_B
          Diag_C(ikl,isl) =  dz_8SV(isl)        *LSdzsv(ikl)
     .                      -Implic            * Elem_C

          Term_D(ikl,isl) = (dz_8SV(isl) *LSdzsv(ikl)
     .                      +Explic      *Elem_A     )*eta_SV(ikl,isl-1)
     .                    + (dz34SV(isl) *LSdzsv(ikl)
     .                      +Explic      *Elem_B     )*eta_SV(ikl,isl)
     .                    + (dz_8SV(isl) *LSdzsv(ikl)
     .                      +Explic      *Elem_C     )*eta_SV(ikl,isl+1)
     .                    + (bKdtSV(ist,ikp)* dzi_SV(isl+1)
     .                      +bKdtSV(ist,ik0)*(dziiSV(isl+1)
     .                                       -dzi_SV(isl)  )
     .                      -bKdtSV(ist,ikm)* dziiSV(isl)   )
     .                                      * LSdzsv(ikl)
     .                    -  dt__SV         * Rootsv(ikl,isl)/rhoWat
        END DO
      END DO

           isl=-nsol
        DO ikl=1,knonv
          ist      =         isotSV(ikl)
          FreeDr   =         iWaFSV(ikl)       *  min(ist,1)
!         FreeDr   =         FreeD0            *  min(ist,1) ! Free Drainage
          ik0      = nkhy *  eta_SV(ikl,isl  ) / etadSV(ist)
          ikp      = nkhy *  eta_SV(ikl,isl+1) / etadSV(ist)
          Elem_A   =         0.
          Elem_B   =      - (Dhydtz(ikl,isl+1) 
     .                      -aKdtSV(ist,ik0)*(dziiSV(isl+1)*LSdzsv(ikl)
     .                                       -FreeDr                  ))
          Elem_C   =         Dhydtz(ikl,isl+1) 
     .                    +  aKdtSV(ist,ikp)* dzi_SV(isl+1)*LSdzsv(ikl)
          Diag_A(ikl,isl) =  0.
          Diag_B(ikl,isl) =  dz78SV(isl) *LSdzsv(ikl)
     .                      -Implic      *Elem_B
          Diag_C(ikl,isl) =  dz_8SV(isl) *LSdzsv(ikl)
     .                      -Implic      *Elem_C

          Term_D(ikl,isl) = (dz78SV(isl) *LSdzsv(ikl)
     .                      +Explic      *Elem_B     )*eta_SV(ikl,isl)
     .                    + (dz_8SV(isl) *LSdzsv(ikl)
     .                      +Explic      *Elem_C     )*eta_SV(ikl,isl+1)
     .                    + (bKdtSV(ist,ikp)* dzi_SV(isl+1)*LSdzsv(ikl)
     .                      +bKdtSV(ist,ik0)*(dziiSV(isl+1)*LSdzsv(ikl)
     .                                       -FreeDr                  ))
     .                    -  dt__SV         * Rootsv(ikl,isl)/rhoWat
        END DO

           isl=0
        DO ikl=1,knonv
          ist      =         isotSV(ikl)
          ikm      = nkhy *  eta_SV(ikl,isl-1) / etadSV(ist)
          ik0      = nkhy *  eta_SV(ikl,isl)   / etadSV(ist)
          Elem_A   =         Dhydtz(ikl,isl)   
     .                    -  aKdtSV(ist,ikm)* dziiSV(isl)*LSdzsv(ikl)
          Elem_B   =      - (Dhydtz(ikl,isl)
     .                      +aKdtSV(ist,ik0)* dzi_SV(isl)*LSdzsv(ikl))
          Elem_C   =         0.
          Diag_A(ikl,isl) =  dz_8SV(isl) *LSdzsv(ikl)
     .                    -  Implic      *Elem_A
          Diag_B(ikl,isl) =  dz78SV(isl) *LSdzsv(ikl)
     .                    -  Implic      *Elem_B
          Diag_C(ikl,isl) =  0.

          Term_D(ikl,isl) = (dz_8SV(isl) *LSdzsv(ikl)
     .                      +Explic      *Elem_A     )*eta_SV(ikl,isl-1)
     .                    + (dz78SV(isl) *LSdzsv(ikl)
     .                      +Explic      *Elem_B     )*eta_SV(ikl,isl)
     .                    - (bKdtSV(ist,ik0)* dzi_SV(isl)
     .                      +bKdtSV(ist,ikm)* dziiSV(isl))*LSdzsv(ikl)
     .            + dt__SV *(HLs_sv(ikl)    *     (1-min(1,isnoSV(ikl)))
     .                                                   / Lx_H2O(ikl)
     .                      +drr_SV(ikl)
     .                      -Rootsv(ikl,isl)             )/rhoWat 
        END DO


! Tridiagonal Elimination
! =======================

! Forward  Sweep
! ^^^^^^^^^^^^^^
        DO ikl=  1,knonv
          Aux__P(ikl,-nsol) = Diag_B(ikl,-nsol)
          Aux__Q(ikl,-nsol) =-Diag_C(ikl,-nsol)/Aux__P(ikl,-nsol)
        END DO

      DO   isl=-nsol+1,0   
        DO ikl=      1,knonv
          Aux__P(ikl,isl)   = Diag_A(ikl,isl)  *Aux__Q(ikl,isl-1)
     .                       +Diag_B(ikl,isl)
          Aux__Q(ikl,isl)   =-Diag_C(ikl,isl)  /Aux__P(ikl,isl)
        END DO
      END DO

        DO ikl=      1,knonv
          eta_SV(ikl,-nsol) = Term_D(ikl,-nsol)/Aux__P(ikl,-nsol)
        END DO

      DO   isl=-nsol+1,0   
        DO ikl=      1,knonv
          eta_SV(ikl,isl)   =(Term_D(ikl,isl)
     .                       -Diag_A(ikl,isl)  *eta_SV(ikl,isl-1))
     .                                         /Aux__P(ikl,isl)
        END DO
      END DO

! Backward Sweep
! ^^^^^^^^^^^^^^
      DO   isl=-1,-nsol,-1
        DO ikl= 1,knonv
          eta_SV(ikl,isl)   = Aux__Q(ikl,isl)  *eta_SV(ikl,isl+1)
     .                                         +eta_SV(ikl,isl)
        END DO
      END DO


! Horton RunOFF Intensity
! =======================

      DO   isl=0,-nsol,-1
        DO ikl=1,knonv
          ist    =   isotSV(ikl)                   ! Soil Type
          SatRat =  (eta_SV(ikl,isl)-etadSV(ist))  ! OverSaturation Rate
     .              *rhoWat         *dzAvSV(isl)   !
     .                              *LSdzsv(ikl)   !
     .                              /dt__SV        !
          SoRnOF(ikl)     =          SoRnOF(ikl)   !
     .                    + max(zer0,SatRat)       !
          eta_SV(ikl,isl) = max(eps6               !
     .                         ,eta_SV(ikl,isl))   !
          eta_SV(ikl,isl) = min(eta_SV(ikl,isl)    !
     .                         ,etadSV(ist)    )   !
        END DO
      END DO

! OUTPUT/Verification: Soil Vertic.Discret.
! #so     write(6,6010) 

 6010     format(/,1x)
      DO   isl= 0,-nsol,-1
        DO ikl= 1,knonv
          ist      =          isotSV(ikl)
          ikp      = nkhy  *  eta_SV(ikl,isl)  /etadSV(ist)
          Khydsv(ikl,isl)   =(aKdtSV(ist,ikp)  *eta_SV(ikl,isl)
     .                       +bKdtSV(ist,ikp)) *2.0/dt__SV
! OUTPUT/Verification: Soil Vertic.Discret.
! #so     write(6,6011) ikl,isl,eta_SV(ikl,isl)*1.e3,
! #so.                  ikp,    aKdtSV(ist,ikp),bKdtSV(ist,ikp),
! #so.                          Khydsv(ikl,isl)
 6011     format(2i3,f8.1,i3,3e12.3)

        END DO
      END DO


! Additional RunOFF Intensity
! ===========================

        DO ikl=1,knonv
          ist      =          isotSV(ikl)
          ik0      = nkhy  *  etaaux(ikl,-nsol  ) /etadSV(ist)
          FreeDr   =          iWaFSV(ikl)       *  min(ist,1)
!         FreeDr   =          FreeD0            *  min(ist,1) ! Free Drainage
          SoRnOF(ikl) =  SoRnOF(ikl)
     .                + (aKdtSV(ist,ik0)*(etaaux(ikl,-nsol)*Explic
     .                                   +eta_SV(ikl,-nsol)*Implic)
     .                 + bKdtSV(ist,ik0)                           )
     .                 * FreeDr          *rhoWat           /dt__SV

! Full Run OFF: Update
! ~~~~~~~~~~~~~~~~~~~~
          RnofSV(ikl) = RnofSV(ikl)     + SoRnOF(ikl)
        END DO


! Temperature Correction due to a changed Soil Energy Content
! ===========================================================

! REMARQUE: Mettre en oeuvre le couplage humidite-energie
! ^^^^^^^^


! Bumps/Asperites Treatment
! =========================

! Average over Bump Depth (z0soil)
! --------------------------------

c #BP       z_Bump      = 0.
c #BP     DO ikl=1,knonv
c #BP       etBump(ikl) = 0.
c #BP     END DO

c #BP DO     isl=0,-nsol,-1
c #BP       z0Bump      = z_Bump
c #BP       z_Bump      = z_Bump      +  dzAvSV(isl)
c #BP   IF (z_Bump.lt.z0soil)                                       THEN
c #BP     DO ikl=1,knonv
c #BP       etBump(ikl) = etBump(ikl) +  dzAvSV(isl)   *eta_SV(ikl,isl)
c #BP     END DO
c #BP   END IF
c #BP   IF (z_Bump.gt.z0soil.AND.z0Bump.lt.z0soil)                  THEN
c #BP     DO ikl=1,knonv
c #BP       etBump(ikl) = etBump(ikl) + (z0soil-z0Bump)*eta_SV(ikl,isl)
c #BP       etBump(ikl) = etBump(ikl) /  z0soil
c #BP     END DO
c #BP   END IF
c #BP END DO


! Correction
! ----------

c #BP       z_Bump      = 0.
c #BP DO     isl=0,-nsol,-1
c #BP       z0Bump =  z_Bump
c #BP       z_Bump =  z_Bump +dzAvSV(isl)
c #BP   IF (z_Bump.lt.z0soil)                                       THEN
c #BP     DO ikl=1,knonv
c #BP       eta_SV(ikl,isl) = etBump(ikl)
c #BP     END DO
c #BP   END IF
c #BP   IF (z_Bump.gt.z0soil.AND.z0Bump.lt.z0soil)                  THEN
c #BP       dzBump          =    z_Bump -  z0soil
c #BP     DO ikl=1,knonv
c #BP       eta_SV(ikl,isl) =(etBump(ikl)    *(dzAvSV(isl)-dzBump)
c #BP.                      + eta_SV(ikl,isl)*             dzBump)
c #BP.                      /                  dzAvSV(isl)
c #BP     END DO
c #BP   END IF
c #BP END DO


! Positive Definite
! =================

c #BP DO   isl= 0,-nsol,-1
c #BP   DO ikl= 1,knonv
c #BP     eta_SV(ikl,isl)   =          max(eps6,eta_SV(ikl,isl))
c #BP   END DO
c #BP END DO


! OUTPUT/Verification: H2O    Conservation: Water  Budget (OUT)
! #m0   DO ikl=1,knonv
! #m0     Wats_d(ikl) = Wats_d(ikl)                    ! 
! #m0.                + drr_SV(ikl)     *0.00          ! Precipitation is
!                                        \______________ already included
! #m0.                + HLs_sv(ikl)
! #m0.          *(1-min(isnoSV(ikl),1)) /Lx_H2O(ikl)   ! Evaporation
! #m0.                - SoRnOF(ikl)                    ! Soil RunOFF Contrib.
! #m0     Wats_1(ikl) = 0.                             !

! OUTPUT/Verification: H2O    Conservation
! #mw     Evapor(ikl) = HLs_sv(ikl)     *dt__SV        !
! #mw.          *(1-min(isnoSV(ikl),1)) /Lx_H2O(ikl)   !

! #m0   END DO

! #m0 DO   isl= -nsol,0
! #m0   DO ikl=1,knonv
! #m0     Wats_d(ikl) = Wats_d(ikl)                    ! 
! #m0.                - Rootsv(ikl,isl)                ! Root Extract.
! #m0   END DO
! #m0 END DO 
! #m0   DO ikl=1,knonv
! #m0     Wats_d(ikl) = Wats_d(ikl)     *dt__SV        ! 
! #m0   END DO

! #m0      isl= -nsol
! #m0   DO ikl=1,knonv
! #m0     Wats_1(ikl) = Wats_1(ikl) 
! #m0.      + rhoWat *( eta_SV(ikl,isl)   *dz78SV(isl)
! #m0.                + eta_SV(ikl,isl+1) *dz_8SV(isl) ) *LSdzsv(ikl)
! #m0   END DO

! #m0 DO   isl= -nsol+1,-1
! #m0   DO ikl=1,knonv
! #m0     Wats_1(ikl) = Wats_1(ikl) 
! #m0.      + rhoWat *( eta_SV(ikl,isl)   *dz34SV(isl)
! #m0.                +(eta_SV(ikl,isl-1)
! #m0.                 +eta_SV(ikl,isl+1))*dz_8SV(isl) ) *LSdzsv(ikl)
! #m0   END DO
! #m0 END DO 

! #m0      isl=  0
! #m0   DO ikl=1,knonv
! #m0     Wats_1(ikl) = Wats_1(ikl) 
! #m0.      + rhoWat *( eta_SV(ikl,isl)   *dz78SV(isl)
! #m0.                + eta_SV(ikl,isl-1) *dz_8SV(isl) ) *LSdzsv(ikl)
! #m0   END DO


! OUTPUT/Verification: H2O    Conservation
! #mw IF (.NOT.mwopen)                                              THEN
! #mw          mwopen = .true.
! #mw          open(unit=42,status='unknown',file='SISVAT_qSo.vw')
! #mw          rewind 42
! #mw   write(42,42)
 42     format('SubRoutine SISVAT_qSo: Local Water Budget',
     .       /,'=========================================')
! #mw END IF
! #mw         timewr=timewr + dt__SV
! #mw         hourwr=3600.0
! #mw IF (mod(timewr,hourwr).lt.eps6)
! #mw.  write(42,420)timewr/hourwr
 420    format(11('-'),'----------+--------------+-',
     .          3('-'),'----------+--------------+-',
     .                 '----------------+----------------+',
     .       /,f8.2,3x,'Wats_0(1) |    Wats_d(1) | ',
     .              3x,'Wats_1(1) | W_0+W_d-W_1  | ',
     .                 '   Soil Run OFF |   Soil Evapor. |',
     .       /,11('-'),'----------+--------------+-',
     .          3('-'),'----------+--------------+-',
     .                 '----------------+----------------+')
! #mw   write(42,421)   Wats_0(1),Wats_d(1)
! #mw.                 ,Wats_1(1) 
! #mw.                 ,Wats_0(1)+Wats_d(1)-Wats_1(1)
! #mw.                 ,SoRnOF(1),Evapor(1)
 421    format(8x,f12.6,' + ',f12.6,' - ',f12.6,' = ',f12.6,' | ',f12.6,
     .      '      ',f15.6)

      return
      end


      subroutine SISVAT_wEq( labWEq ,istart)

!--------------------------------------------------------------------------+
!   MAR          SISVAT_wEq                           Sat 12-Feb-2012  MAR |
!     SubRoutine SISVAT_wEq computes the Snow/Ice  Water  Equivalent       |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_wEq.ve        | #ve: OUTPUT/Verification: Snow/Ice Water Eqv. |
!                          |      unit 45, SubRoutine  SISVAT_wEq **ONLY** |
!--------------------------------------------------------------------------+


       


! Global Variables
! ================

      USE VAR_SV
      USE VARxSV

      IMPLICIT NONE

      character*6         labWEq
      integer             istart

      logical             logWEq
      common/SISVAT_wEq_L/logWEq


! Local  Variables
! ================

      integer  ikl   ,isn
      real     SnoWEQ,IceWEQ


! Switch Initialization
! =====================

      IF (.NOT.logWEq)                                              THEN
               logWEq = .true.
               open(unit=45,status='unknown',file='SISVAT_wEq.ve')
               rewind    45
      END IF


! Snow Water Equivalent
! =====================

           ikl   = 1
      IF          (isnoSV(ikl).gt.iiceSV(ikl))                      THEN

          SnoWEQ = 0.
        DO isn   = iiceSV(ikl)+1 ,isnoSV(ikl)
          SnoWEQ = SnoWEQ       + ro__SV(ikl,isn) * dzsnSV(ikl,isn) 
        END DO

      END IF


! Ice  Water Equivalent
! =====================

      IF        (iiceSV(1).gt.0)                                    THEN

          IceWEQ = 0.
        DO isn   =             1 ,iiceSV(ikl)
          IceWEQ = IceWEQ       + ro__SV(ikl,isn) * dzsnSV(ikl,isn) 
        END DO

      END IF


! OUTPUT
! ======

      IF (istart.eq.1)                                              THEN
        write(45,45)dahost,i___SV(lwriSV(1)),j___SV(lwriSV(1)),
     .              n___SV(lwriSV(1))
 45     format(a18,10('-'),'Pt.',3i4,60('-'))
      END IF

      write(45,450) labWEq,IceWEQ,iiceSV(ikl),SnoWEQ
     .                    ,IceWEQ+SnoWEQ,isnoSV(ikl)
     .                                  ,drr_SV(ikl)*dt__SV
     .                                  ,dsn_SV(ikl)*dt__SV
     .                                  ,BufsSV(ikl)
 450  format(a6,3x,'  I+S =',f11.4,'(',i2,') +',f11.4,' =',
     .                       f11.4,'(',i2,')',
     .             '  drr =', f7.4,
     .             '  dsn =', f7.4,
     .             '  Buf =', f7.4)

      return
      end
