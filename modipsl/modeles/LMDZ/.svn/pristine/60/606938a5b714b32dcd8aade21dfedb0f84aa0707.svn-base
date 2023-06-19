
      subroutine PHY_Atm_DY_ALLOC

!------------------------------------------------------------------------------+
!                                                         Fri  7-Jun-2013  MAR |
!                                                                              |
!     subroutine PHY_Atm_DY_ALLOC  allocates prognostic variables of           |
!                Atmospheric Turbulence Scheme used by MAR                     |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,               Fri  7-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


      use Mod_Real
      use Mod_PHY____grd
      use Mod_PHY_DY_kkl


      IMPLICIT NONE



! =================================
! ALLOCATION Mod_PHY_DY_kkl - BEGIN
! =================================

      allocate  ( psa_DY(kcolp     ) )                       !  Pressure    Thickness                                    [kPa]
      allocate  ( ExnrDY(kcolp,mzpp) )                       !  Potential   Exner                        pa  **(R/Cp)    [xxx]
      allocate  ( Z___DY(kcolp,mzpp) )                       !  Geopotential, level k    , i.e. =  gZ(k)         /  g  [m2/s2]
      allocate  ( ZmidDY(kcolp,mzpp) )                       !  Geopotential, level k-1/2, i.e. = (gZ(k)+gZ(k-1))/(2g) [m2/s2]

      allocate  ( TmidDY(kcolp,mzpp) )                       !  Temperature , level k+1/2, i.e. = (Ta(k)+Ta(k+1))/ 2       [K]
      allocate  ( Ta__DY(kcolp,mzpp) )                       !  Temperature , level k    , i.e. =  Ta(k)                   [K]
      allocate  ( pkt_DY(kcolp,mzpp) )                       !  Pseudo P.T. , level k                                      [K]
      allocate  ( windDY(kcolp,mzp ) )                       !  Wind Speed, Horizontal                                   [m/s]
      allocate  ( ua__DY(kcolp,mzp ) )                       !  Wind Speed, x-Direction                                  [m/s]
      allocate  ( va__DY(kcolp,mzp ) )                       !  Wind Speed, y-Direction                                  [m/s]
      allocate  ( wa__DY(kcolp,mzp ) )                       !  Wind Speed  z-Direction                                  [m/s]
      allocate  ( roa_DY(kcolp,mzp ) )                       !  Air Volumic Mass, Layer k                              [Mg/m3]
      allocate  ( roamDY(kcolp,mzp ) )                       !  Air Volumic Mass, Level k+1/2                          [Mg/m3]

      allocate  ( qv__DY(kcolp,mzpp) )                       !  Specific    Humidity                                   [kg/kg]
! #LD allocate  ( ld_H2O(kcolp,mzpp) )                       !  Loading    (Humidity, Hydrometeors, Aerosols ...)          [-]

! =================================
! ALLOCATION Mod_PHY_DY_kkl -   END
! =================================



      end subroutine PHY_Atm_DY_ALLOC
