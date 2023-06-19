      module Mod_PHY_DY_kkl

!--------------------------------------------------------------------------+
!                                                     Tue  4-Jun-2013  MAR |
!     module Mod_PHY_DY_kkl contains the main (prognostic) variables of    |
!                MAR Dynamics Variabbles on MAR Physics Grid               |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,           Tue  4-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



! Atm_DY INPUT        Variables
! -----------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  psa_DY       !  Pressure    Thickness                                    [kPa]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  ExnrDY       !  Potential   Exner                        pa  **(R/Cp)    [xxx]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  Z___DY       !  Geopotential, level k    , i.e. =  gZ(k)         /  g  [m2/s2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  ZmidDY       !  Geopotential, level k-1/2, i.e. = (gZ(k)+gZ(k-1))/(2g) [m2/s2]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  pkt_DY       !  Potential   Temperature, divided by (100 kPa)**(R/Cp)  [K/xxx]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  TmidDY       !  Temperature , level k+1/2, i.e. = (Ta(k)+Ta(k+1))/ 2       [K]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  Ta__DY       !  Temperature , level k    , i.e. =  Ta(k)                   [K]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  windDY       !  Wind Speed, Horizontal                                   [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  ua__DY       !  Wind Speed  x-Direction                                  [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  va__DY       !  Wind Speed  y-Direction                                  [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  wa__DY       !  Wind Speed  z-Direction                                  [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  roa_DY       !  Air Volumic Mass, Layer k                              [Mg/m3]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  roamDY       !  Air Volumic Mass, Level k+1/2                          [Mg/m3]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  qv__DY       !  Specific    Humidity                                   [kg/kg]
! #LD real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  ld_H2O       !  Loading    (Humidity, Hydrometeors, Aerosols ...)          [-]



! Atm_DY INPUT/OUTPUT Variables
! -----------------------------



! Atm_DY OUTPUT       Variables
! -----------------------------


      end module Mod_PHY_DY_kkl
