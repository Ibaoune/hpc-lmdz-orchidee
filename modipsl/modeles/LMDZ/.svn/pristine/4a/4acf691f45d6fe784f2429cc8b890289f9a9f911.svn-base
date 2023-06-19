      module Mod_PHY_S0_dat

!--------------------------------------------------------------------------+
!                                                     Sat 27-Apr-2013  MAR |
!     module Mod_PHY_S0_dat contains the main data                   of    |
!                Insolation                   Scheme                       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Thu 25-Apr-2013      |
!           Last Modification by H. Gallee,           Sat 27-Apr-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



      integer, SAVE                                       ::  kBP    = 0      ! Time                                   [kyr B.P.]
      real(kind=real8), SAVE                              ::  rsunS0          ! Insolation, Top of the Atmosphere          [W/m2] 
      real(kind=real8), SAVE                              ::  cszEPS = 1.e-3  ! Minimum accepted cos(Sun zenith.Dist.)        [-] 
      real(kind=real8), SAVE                              ::  ecc             ! Earth Orbit Eccentricity                      [-]
      real(kind=real8), SAVE                              ::  perh            ! Earth Orbit Longitude of the Perihelion  [degree]
      real(kind=real8), SAVE                              ::  xob             ! Earth Orbit Obliquity                    [degree]
      real(kind=real8), SAVE                              ::  ecc_EO(-10:0)   ! Earth Orbit Eccentricity                      [-] 
      real(kind=real8), SAVE                              ::  perhEO(-10:0)   ! Earth Orbit Longitude of the Perihelion  [degree]
      real(kind=real8), SAVE                              ::  xob_EO(-10:0)   ! Earth Orbit Obliquity                    [degree]

      data ecc_EO(  0)                                 /  0.01673 /     ! Eccentricity 
      data perhEO(  0)                                 /102.4     /     ! Longitude of the Perihelion              [degree]
      data xob_EO(  0)                                 / 23.445   /     ! Obliquity                                [degree]

      data ecc_EO( -6)                                 /  0.018682/     ! Eccentricity
      data perhEO( -6)                                 /  0.87    /     ! Longitude of the Perihelion              [degree]
      data xob_EO( -6)                                 /  24.105  /     ! Obliquity                                [degree]

      data ecc_EO(-10)                                 /  0.019419/     ! Eccentricity
      data perhEO(-10)                                 /294.81    /     ! Longitude of the Perihelion              [degree]
      data xob_EO(-10)                                 / 24.226   /     ! Obliquity                                [degree]


      real(kind=real8), SAVE                              ::  pirr            ! 180 / pi / 3600               [degree/radian/sec]

      real(kind=real8), SAVE                              ::  om = 0.0172142  ! Radian advance in one day                [radian]
      real(kind=real8), SAVE                              ::  Tyear = 365.25  ! Length         of one year                  [day]
      real(kind=real8), SAVE                              ::  Step            ! Advance on Earth Orbit in 1 day      [degree/day]

      real(kind=real8), SAVE                              ::  xee             ! Square    of Eccentricity (Ecc)               [-]
      real(kind=real8), SAVE                              ::  xse             ! Square Root (1-Ecc**2)                        [-]
      real(kind=real8), SAVE                              ::  xe3             ! Square Root (1-Ecc**2)  X  Ecc                [-]
      real(kind=real8), SAVE                              ::  xl              ! Longitude of Aphelion                    [degree]
      real(kind=real8), SAVE                              ::  xllp            ! Longitude of Aphelion                    [radian]
      real(kind=real8), SAVE                              ::  so              ! sinus     of Obliquity                        [-]
      real(kind=real8), SAVE                              ::  xlam            ! true long. sun for mean long. = 0        [degree]



      end module Mod_PHY_S0_dat
