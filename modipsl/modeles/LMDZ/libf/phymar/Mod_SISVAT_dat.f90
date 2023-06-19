      module Mod_SISVAT_dat


!--------------------------------------------------------------------------+
!                                                     Sun 30-Jun-2013  MAR |
!     module Mod_SISVAT_dat contains the constants of the                  |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!                                                                          |
!   CAUTION: Soil Hydraulic Parameters: Please VERIFY ICE DATA             |
!   ^^^^^^^                                                                |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #fb: Fractions of total Solar Irradiance (Feagle and Businger 1981)  |
!                                                                          |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Wed 30-Jan-2013      |
!                    modified by H. Gallee,           30-Jun-May-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


! Global   Variables
! ==================

      use Mod_Real
      use Mod_SISVAT_dim


      IMPLICIT NONE



! Internal Variables
! ==================

      integer, SAVE                 ::   ivg,iso



! SISVAT Global Variables
! =======================

      real(kind=real8), SAVE        ::   eps_21 = 1.e-21         ! Arbitrary  very small value

      real(kind=real8), SAVE        ::   slop1d = 0.00           ! Local Slope (1-d simulations)
      real(kind=real8), SAVE        ::   Adz0dt                  ! Decay of Angle(Wind,Sastrugi) Influence on z0

      real(kind=real8), SAVE        ::   WatIsv =  273.16        ! used in Evaluation of Saturation Specific Humidity over Ice (see Dudhia)
      real(kind=real8), SAVE        ::   ExpIsv = 6150.00        ! used in Evaluation of Saturation Specific Humidity over Ice (see Dudhia)



! Snow
! ----

      integer, SAVE, dimension(5)   ::   istdSV                  ! Snow History
      data                        (istdSV(iso),iso=1,5) /1,2,3,4,5/
                                                           ! 1:             faceted cristal
                                                           ! 2: liq.watr/no faceted cristal befor
                                                           ! 3: liq.watr/   faceted cristal befor

      real(kind=real8), SAVE   ::         Cn_dSV = 2105.         ! Snow Heat Capacity          [J/kg/K] 
                                                           ! Loth et al. 1993, JGR 98 D6
      real(kind=real8), SAVE   ::         SMndSV =    1.00       ! New Snow Layer Min.Thickn. [mm w.e.]
      real(kind=real8), SAVE   ::         G1_dSV =   99.00       ! Conversion 0/99-->0/1 
                                                           ! Sphericity/Dendricity

                                                           ! Optical Diameter of:
      real(kind=real8), SAVE   ::         DDcdSV =    1.00       ! Dendritic     Crystals    [0.0001 m]
      real(kind=real8), SAVE   ::         DFcdSV =    4.00       ! Young Faceted Crystals    [0.0001 m]
      real(kind=real8), SAVE   ::         DScdSV =    3.00       ! Small         Crystals    [0.0001 m]

                                                           ! Actual  Diameter of:
      real(kind=real8), SAVE   ::         ADSdSV =    4.00       ! Small         Crystals    [0.0001 m]


! Fractions of total Solar Irradiance  in 3 spectral Intervals
! (see Feagle and Businger 1981, Int.Geoph.Ser. 25, p.215-222)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #fb real(kind=real8), SAVE   ::         So1dSV =    0.606      !                 0.3--0.8mim Interval
! #fb real(kind=real8), SAVE   ::         So2dSV =    0.301      !                 0.8--1.5mim Interval
! #fb real(kind=real8), SAVE   ::         So3dSV =    0.093      !                 1.5--2.8mim Interval

      real(kind=real8), SAVE   ::         So1dSV =    0.580      ! Tuning ETH camp 0.3--0.8mim Interval
      real(kind=real8), SAVE   ::         So2dSV =    0.320      ! Tuning ETH camp 0.8--1.5mim Interval
      real(kind=real8), SAVE   ::         So3dSV =    0.100      ! Tuning ETH camp 1.5--2.8mim Interval

      real(kind=real8), SAVE   ::         aI1dSV =    0.40       ! Minimum bare ICE albedo          [-]
      real(kind=real8), SAVE   ::         aI2dSV =    0.45       ! Maximum bare ICE albedo          [-]
      real(kind=real8), SAVE   ::         aI3dSV =    0.65       ! ICE lense albedo at 800 kg/m3    [-]
                                                           ! and minimum pure snow albedo

! Water in Snow
! -------------
      real(kind=real8), SAVE   ::         ws0dSV =    0.07       ! Irreducible Water Saturation in Snow
                                                           ! Coleou et al., 1998, A.Gla.26, 64-68
      real(kind=real8), SAVE   ::         roCdSV =  800.00       ! Pore Hole Close OFF Density  [kg/m3]
                                                           ! Greuell & Konzelmann (1994)
                                                           ! Glob.Plan.Change 9, 4.5 p.100
      real(kind=real8), SAVE   ::         ru_dSV =  200.00       ! Surficial Water Scale Factor [kg/m2]


! Ice
! ---

      real(kind=real8), SAVE   ::         CdidSV =    2.10       ! Conductivity of pure  Ice    [W/m/K]



! Vegetation                (SVAT Classification)
! -----------------------------------------------

      real(kind=real8), SAVE   ::         DH_dSV(0:nvgt)         ! Displacement            Height   [m]
      real(kind=real8), SAVE   ::         Z0mdSV(0:nvgt)         ! Roughness  Length for Momentum   [m]
      real(kind=real8), SAVE   ::         StodSV(0:nvgt)         ! Minimum    Stomatal Resistance [s/m]
      real(kind=real8), SAVE   ::         PR_dSV(0:nvgt)         ! Internal Plant      Resistance   [s]
      real(kind=real8), SAVE   ::         rbtdSV(0:nvgt)         ! Root Fraction Beta Coefficient   [-]

      data     (DH_dSV(ivg),                              &! Displacement            Height   [m]
     &                 Z0mdSV(ivg),                       &! Roughness  Length for Momentum   [m]
     &                        StodSV(ivg),                &! Minimum    Stomatal Resistance [s/m]
     &                               PR_dSV(ivg),         &! Internal Plant      Resistance   [s]
     &                                      rbtdSV(ivg),  &! Root beta coeffient              [-]
     &                                              ivg=0,nvgt) &
     &         /0.00,  0.01,5000.,   0.0,   0.000,        &!  0 NO     VEGETATION
     &          0.07,  0.01,  50.,   0.5e9, 0.961,        &!  1 CROPS      LOW
     &          0.21,  0.03,  50.,   0.5e9, 0.961,        &!  2 CROPS      MEDIUM
     &          0.70,  0.10,  50.,   0.5e9, 0.961,        &!  3 CROPS      HIGH
     &          0.07,  0.01,  50.,   0.5e9, 0.943,        &!  4 GRASS      LOW
     &          0.21,  0.03,  50.,   0.5e9, 0.964,        &!  5 GRASS      MEDIUM
     &          0.70,  0.10,  50.,   0.5e9, 0.972,        &!  6 GRASS      HIGH
     &          1.40,  0.20,  10.,   1.0e9, 0.968,        &!  7 BROADLEAF  LOW
     &          5.60,  0.80,  10.,   1.0e9, 0.962,        &!  8 BROADLEAF  MEDIUM
     &         14.00,  2.00,  10.,   1.0e9, 0.962,        &!  9 BROADLEAF  HIGH
     &          1.40,  0.20,  10.,   1.0e9, 0.971,        &! 10 NEEDLELEAF LOW
     &          5.60,  0.80,  10.,   1.0e9, 0.976,        &! 11 NEEDLELEAF MEDIUM
     &         14.00,  2.00,  10.,   1.0e9, 0.976/         ! 12 NEEDLELEAF HIGH


      real(kind=real8), SAVE   ::         pscdSV =  250.         ! Critical Leaf Water Potential    [m]
      real(kind=real8), SAVE   ::         StxdSV = 5000.         ! maXimum  Stomatal   Resistance [s/m]
      real(kind=real8), SAVE   ::         LAIdSV =    4.         ! maximum  LAI


      real(kind=real8), SAVE   ::         f__ust(0:nvgt)         ! 

      data          (f__ust(ivg),      ivg=0,nvgt)        &!
     &              /1.00,                                &!  0 NO     VEGETATION
     &               1.20,                                &!  1 CROPS      LOW
     &               5.00,                                &!  2 CROPS      MEDIUM
     &              10.00,                                &!  3 CROPS      HIGH
     &               1.20,                                &!  4 GRASS      LOW
     &               5.00,                                &!  5 GRASS      MEDIUM
     &              10.00,                                &!  6 GRASS      HIGH
     &               5.00,                                &!  7 BROADLEAF  LOW
     &              10.00,                                &!  8 BROADLEAF  MEDIUM
     &              12.00,                                &!  9 BROADLEAF  HIGH
     &              10.00,                                &! 10 NEEDLELEAF LOW
     &              12.00,                                &! 11 NEEDLELEAF MEDIUM
     &              50.00                               /  ! 12 NEEDLELEAF HIGH


! Soil
! ----

      real(kind=real8), SAVE   ::         rcwdSV = 4.180e+6      ! Density * Water Specific Heat

! Soil Vertical Discretization
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real8), SAVE   ::         dz_dSV(-nsol:0)        ! Vertical  Discretization
      data                         (dz_dSV(iso),iso=-4,0) &!
     &                     /0.72,0.20,0.060,0.019,0.001/   ! Layer's Thickness

      real(kind=real8), SAVE   ::         zz_dSV                 ! Soil      Thickness

! Soil Hydraulic Parameters (USDA Classification)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer,       parameter ::         nsot   =   12

      real(kind=real8), SAVE   ::         etadSV(0:nsot)         ! Water Content at Saturation  [kg/kg]
      real(kind=real8), SAVE   ::         psidSV(0:nsot)         ! Water Succion at Saturation      [m]
      real(kind=real8), SAVE   ::         Ks_dSV(0:nsot)         ! Hydraulic Conductivity
                                                                 !               at Saturation    [m/s]
      real(kind=real8), SAVE   ::         bCHdSV(0:nsot)         ! Clapp-Hornberger Coefficient b   [-]

      data         (etadSV(iso),                                     &
     &                        psidSV(iso),                           &
     &                                  Ks_dSV(iso),                 &
     &                                            bCHdSV(iso),       &
     &                                                   iso=0,nsot) &
     &             / 1.000,    1.000,   0.0e00,     1.00, &!  0 WATER
     &               0.395,    0.121, 176.0e-6,     4.05, &!  1 SAND
     &               0.410,    0.090, 156.3e-6,     4.38, &!  2 LOAMY      SAND
     &               0.435,    0.218,  34.1e-6,     4.90, &!  3 SANDY      LOAM
     &               0.485,    0.786,   7.2e-6,     5.30, &!  4 SILT       LOAM
     &               0.451,    0.478,   7.0e-6,     5.39, &!  5            LOAM
     &               0.420,    0.299,   6.3e-6,     7.12, &!  6 SANDY CLAY LOAM
     &               0.477,    0.356,   1.7e-6,     7.75, &!  7 SILTY CLAY LOAM
     &               0.476,    0.630,   2.5e-6,     8.52, &!  8       CLAY LOAM
     &               0.426,    0.153,   2.2e-6,    10.40, &!  9 SANDY CLAY
     &               0.492,    0.490,   1.0e-6,    10.40, &! 10 SILTY CLAY
     &               0.482,    0.405,   1.3e-6,    11.40, &! 11       CLAY
     &               0.001,    0.001,   0.0e00,     0.02/  ! 12       ICE 


      real(kind=real8), SAVE  ::          ustdmn(0:nsot)
      real(kind=real8), SAVE  ::          claypc(0:nsot)

      data         (ustdmn(iso),                                     &
     &                        claypc(iso),                           &
     &                                                   iso=0,nsot) &
     &             /10.000,  0.0000, &!  0 WATER           !
     &               0.300,  0.0000, &!  1 SAND            !
     &               0.300,  0.0920, &!  2 LOAMY      SAND ! Fal99, Table 2
     &               0.300,  0.1420, &!  3 SANDY      LOAM ! Fal99, Table 2
     &               0.300,  0.1630, &!  4 SILT       LOAM ! Guess (Interpol.)
     &               0.300,  0.1840, &!  5            LOAM ! Fal99, Table 2
     &               0.300,  0.2280, &!  6 SANDY CLAY LOAM ! Guess (Interpol.)
     &               0.300,  0.2720, &!  7 SILTY CLAY LOAM ! Guess (Interpol.)
     &               0.300,  0.3160, &!  8       CLAY LOAM ! Fal99, Table 2
     &               0.300,  0.3750, &!  9 SANDY CLAY      ! Guess (Interpol.)
     &               0.300,  0.4340, &! 10 SILTY CLAY      ! Guess (Interpol.)
     &               0.300,  0.4920, &! 11       CLAY      ! Fal99, Table 2
     &              10.000,  0.0000/  ! 12       ICE       !



! Water Bodies
! ------------

      real(kind=real8), SAVE   ::         vK_dSV = 1000.         ! Diffusivity in Water          [m2/s]
      real(kind=real8), SAVE   ::         SIcMIN =    0.1        ! Sea-Ice Layer Min Thickness      [m]
      real(kind=real8), SAVE   ::         dzSIce(4)              ! Sea-Ice Vertical  Discretisation [m]
      data                          dzSIce /0.5,0.05,0.001,0.0/
      real(kind=real8), SAVE   ::         SIc_OK(2)              ! Sea-Ice Switch
      data                          SIc_OK /1.0,0.00/ 
      real(kind=real8), SAVE   ::         TSIdSV =    0.50       ! Sea-Ice Fraction: SST Scale      [K]
      real(kind=real8), SAVE   ::         OcnMin =    0.05       ! Open Water Min Fraction (S.Hem.) [-]
      real(kind=real8), SAVE   ::         TocnSI =  270.70       ! Ocn Temp. for Full Sea-Ice Cover [K]
      real(kind=real8), SAVE   ::         TOF_SV                 ! Ocn Grid Cell Freez.Temperature  [K]
      real(kind=real8), SAVE   ::         VarSST                 ! Variable (0.) / Fixed    (1.)
      real(kind=real8), SAVE   ::         FixSST                 ! Fixed    (1.) / Variable (0.) SST
      real(kind=real8), SAVE   ::         SSTnud                 ! SST Nudging Parameter


! Auxiliary Variables
! -------------------

      integer,       parameter ::  nkhy=50

      real(kind=real8), SAVE   ::  rocsSV( 0:nsot)               ! Soil Contribution to (ro c)_s
      real(kind=real8), SAVE   ::  etamSV( 0:nsot)               ! Soil Minimum Humidity
      real(kind=real8), SAVE   ::  s1__SV( 0:nsot)               ! ... X eta**( b+2), DR97(3.36)
      real(kind=real8), SAVE   ::  s2__SV( 0:nsot)               ! ... X eta**(2b+3), DR97(3.35)
      real(kind=real8), SAVE   ::  aKdtSV( 0:nsot, 0:nkhy)       ! Khyd=a*eta+b: a * dt
      real(kind=real8), SAVE   ::  bKdtSV( 0:nsot, 0:nkhy)       ! Khyd=a*eta+b: b * dt



      end module Mod_SISVAT_dat
