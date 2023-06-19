MODULE VARdSV

USE VAR_SV

IMPLICIT NONE


! +--SISVAT Global Variables
! +  =======================

      LOGICAL           ::   INI_SV=.false.          ! Initialisation Switch     
      REAL,PARAMETER    ::   eps_21=1.e-21           ! Arbitrary  very small value


! +--Snow
! +  ----

      INTEGER,PARAMETER,DIMENSION(5) ::   istdSV(1:5)=(/1,2,3,4,5/)  ! Snow History

      REAL,PARAMETER    ::   Cn_dSV= 2105.           ! Snow Heat Capacity          [J/kg/K]
      REAL,PARAMETER    ::   SMndSV= 1.00            ! Minimum Thickness of new Layers
      REAL,PARAMETER    ::   G1_dSV= 99.             ! Conversion 0/99-->0/1 
      REAL,PARAMETER    ::   DDcdSV= 1.,DFcdSV= 4.,DScdSV= 3.   
                                                     ! Snow Grains Optical Diameter [1e-4m]
      REAL,PARAMETER    ::   ADSdSV= 4.              ! Snow Grains Actual  Diameter [1e-4m]
      REAL,PARAMETER    ::   So1dSV= 0.580,So2dSV= 0.320,So3dSV= 0.100
                                                     ! Total Solar Irradiance Fractions [-]
                                                     ! Tuning ETH camp 0.3--0.8mim Interval
                                                     ! Tuning ETH camp 0.8--1.5mim Interval
                                                     ! Tuning ETH camp 1.5--2.8mim Interval
                                                     ! So1dSV=0.606,So2dSV=0.301,So3dSV=0.093
      REAL,PARAMETER    ::   aI1dSV= 0.40,aI2dSV= 0.45,aI3dSV= 0.65
                                                     ! Bare Ice Albedo                  [-]
                                                     ! Minimum/Maximum/ICE lense albedo at 
                                                     ! 800 kg/m3 and minimum pure snow albedo
      REAL,PARAMETER    ::   ws0dSV= 0.07            ! Irreducible Water Saturation in Snow
      REAL,PARAMETER    ::   roCdSV= 800.            ! Pore Hole Close OFF Density  [kg/m3]
      REAL,PARAMETER    ::   ru_dSV= 200.            ! Surficial Water Scale Factor [kg/m2]

!C +--Ice
!C +  ---

      REAL,PARAMETER    ::   CdidSV= 2.1

!C +--Vegetation
!C +  ----------

      INTEGER,PARAMETER ::   nvgt=12
      REAL,PARAMETER    ::   DH_dSV(0:nvgt) = (/ 0.00, 0.07, 0.21,      &
     &       0.70, 0.07, 0.21, 0.70, 1.40, 5.60,14.00, 1.40, 5.60,14.00/)
                                                  ! Displacement            Height   [m]

      REAL,PARAMETER    ::   Z0mdSV(0:nvgt) = (/ 0.01, 0.01, 0.03,      &
     &       0.10, 0.01, 0.03, 0.10, 0.20, 0.80, 2.00, 0.20, 0.80, 2.00/)
                                                  ! Roughness  Length for Momentum   [m]

      REAL,PARAMETER    ::   StodSV(0:nvgt) = (/5000.,  50.,  50.,      &
     &        50.,  50.,  50.,  50.,  10.,  10.,  10.,  10.,  10.,  10./)
                                                  ! Minimum    Stomatal Resistance [s/m]

      REAL,PARAMETER    ::   PR_dSV(0:nvgt) = (/  0.0,0.5e9,0.5e9,      & 
     &      0.5e9,0.5e9,0.5e9,0.5e9,1.0e9,1.0e9,1.0e9,1.0e9,1.0e9,1.0e9/)
                                                  ! Roots Fraction Beta Coefficient  [-]

      REAL,PARAMETER    ::   rbtdSV(0:nvgt) = (/0.000,0.961,0.961,      &
     &      0.961,0.943,0.964,0.972,0.968,0.962,0.962,0.971,0.976,0.976/)

!        /0.00,  0.01,5000.,   0.0,   0.000,      !  0 NO     VEGETATION
!         0.07,  0.01,  50.,   0.5e9, 0.961,      !  1 CROPS      LOW
!         0.21,  0.03,  50.,   0.5e9, 0.961,      !  2 CROPS      MEDIUM
!         0.70,  0.10,  50.,   0.5e9, 0.961,      !  3 CROPS      HIGH
!         0.07,  0.01,  50.,   0.5e9, 0.943,      !  4 GRASS      LOW
!         0.21,  0.03,  50.,   0.5e9, 0.964,      !  5 GRASS      MEDIUM
!         0.70,  0.10,  50.,   0.5e9, 0.972,      !  6 GRASS      HIGH
!         1.40,  0.20,  10.,   1.0e9, 0.968,      !  7 BROADLEAF  LOW
!         5.60,  0.80,  10.,   1.0e9, 0.962,      !  8 BROADLEAF  MEDIUM
!        14.00,  2.00,  10.,   1.0e9, 0.962,      !  9 BROADLEAF  HIGH
!         1.40,  0.20,  10.,   1.0e9, 0.971,      ! 10 NEEDLELEAF LOW
!         5.60,  0.80,  10.,   1.0e9, 0.976,      ! 11 NEEDLELEAF MEDIUM
!        14.00,  2.00,  10.,   1.0e9, 0.976/      ! 12 NEEDLELEAF HIGH

                                                  ! Internal Plant      Resistance   [s]
      REAL,PARAMETER    ::   pscdSV = 250.        ! Critical Leaf Water Potential    [m]
      REAL,PARAMETER    ::   StxdSV = 5000.       ! maXimum  Stomatal   Resistance [s/m]
      REAL,PARAMETER    ::   LAIdSV = 4.          ! maximum  LAI


!C +--Soil
!C +  ----

      REAL,PARAMETER    ::   rcwdSV = 4.180e+6    ! Density * Water Specific Heat
      REAL              ::   dz_dSV(-nsol:0)      ! Vertical  Discretization MARSV:
                                                  !/0.72,0.20,0.060,0.019,0.001/ 
                                                  ! Layer's Thickness
      REAL              ::   zz_dSV               ! Soil      Thickness

      INTEGER,PARAMETER ::   nsot=12
      REAL,PARAMETER    ::   etadSV(0:nsot) = (/ 1.000,0.395,0.410,     &
     &     0.435,0.485,0.451,0.420,0.477,0.476,0.426,0.492,0.482,0.001 /)      
                                                  ! Water Content at Saturation  [m3/m3]

      REAL,PARAMETER    ::   psidSV(0:nsot) = (/ 1.000,0.121,0.090,     &
     &     0.218,0.786,0.478,0.299,0.356,0.630,0.153,0.490,0.405,0.001 /)
                                                  ! Water Succion at Saturation      [m]

      REAL,PARAMETER    ::   Ks_dSV(0:nsot) = (/ 0.e00, 176.0e-6,       &
     &          156.3e-6,  34.1e-6,   7.2e-6,   7.0e-6,   6.3e-6,       &
     &            1.7e-6,   2.5e-6,   2.2e-6,   1.0e-6,   1.3e-6,0.0e0 /)
                                                  ! Hydraulic Conductivity
                                                  !               at Saturation    [m/s]
      REAL,PARAMETER    ::   bCHdSV(0:nsot) = (/ 1.00, 4.05, 4.38,      &
     &      4.90, 5.30, 5.39, 7.12, 7.75, 8.52,10.40,10.40,11.40, 0.02 /)
                                                  ! Clapp-Hornberger Coefficient b   [-]

 !     etadSV,   psidSV,   Ks_dSV    bCHdSV  
 !     /1.000,    1.000,   0.0e00,     1.00,      !  0 WATER
 !      0.395,    0.121, 176.0e-6,     4.05,      !  1 SAND
 !      0.410,    0.090, 156.3e-6,     4.38,      !  2 LOAMY      SAND
 !      0.435,    0.218,  34.1e-6,     4.90,      !  3 SANDY      LOAM
 !      0.485,    0.786,   7.2e-6,     5.30,      !  4 SILT       LOAM
 !      0.451,    0.478,   7.0e-6,     5.39,      !  5            LOAM
 !      0.420,    0.299,   6.3e-6,     7.12,      !  6 SANDY CLAY LOAM
 !      0.477,    0.356,   1.7e-6,     7.75,      !  7 SILTY CLAY LOAM 
 !      0.476,    0.630,   2.5e-6,     8.52,      !  8       CLAY LOAM
 !      0.426,    0.153,   2.2e-6,    10.40,      !  9 SANDY CLAY
 !      0.492,    0.490,   1.0e-6,    10.40,      ! 10 SILTY CLAY
 !      0.482,    0.405,   1.3e-6,    11.40,      ! 11       CLAY
 !      0.001,    0.001,   0.0e00,     0.02/      ! 12       ICE 


!C +--Water Bodies
!C +  ------------

      REAL,PARAMETER    ::   vK_dSV = 1000.       ! Diffusivity in Water          [m2/s]
      REAL,PARAMETER    ::   TSIdSV = 0.50        ! Sea-Ice Fraction: SST Scale      [K]

      INTEGER ::  ivg1,iso1  !rajout hjp for ini - check if ivg,iso is possible

!

END MODULE VARdSV
