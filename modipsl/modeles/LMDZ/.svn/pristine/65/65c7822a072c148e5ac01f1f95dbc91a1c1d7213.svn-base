      subroutine PHY_Atm_CP_ALLOC

!------------------------------------------------------------------------------+
!                                                         Mon 17-May-2013  MAR |
!                                                                              |
!     subroutine PHY_Atm_CP_ALLOC  allocates prognostic variables of           |
!                Atmospheric Turbulence Scheme used by MAR                     |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,               Mon 17-May-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


      use Mod_Real
      use Mod_PHY____grd
      use Mod_PHY_CP_kkl
!     use Mod_PHY_CP_xyz


      IMPLICIT NONE



! =================================
! ALLOCATION Mod_PHY_CP_kkl - BEGIN
! =================================

      allocate  ( timeCP(kcolp    ) )                        !  Convective Adjustment Time                                 [s]

      allocate  ( dpktCP(kcolp,mzp) )                        !  Pseudo Potential Temperature Tendency                  [K/X/s]
      allocate  ( dqv_CP(kcolp,mzp) )                        !  Specific         Humidity    Tendency                [kg/kg/s]
      allocate  ( dqw_CP(kcolp,mzp) )                        !  Cloud Droplets Concentration Tendency                [kg/kg/s]
      allocate  ( dqi_CP(kcolp,mzp) )                        !  Cloud Crystals Concentration Tendency                [kg/kg/s]

      allocate  ( drr_CP(kcolp    ) )                        !  Rain (convective)            Tendency                    [m/s]
      allocate  ( rainCP(kcolp    ) )                        !  Rain (convective)                                        [m/s]
      allocate  ( dss_CP(kcolp    ) )                        !  Snow (convective)            Tendency                    [m/s]
      allocate  ( snowCP(kcolp    ) )                        !  Snow (convective)                                        [m/s]
      allocate  ( CAPECP(kcolp    ) )                        !  Convective Available Potentential Energy               [m2/s2]

! =================================
! ALLOCATION Mod_PHY_CP_kkl -   END
! =================================





      end subroutine PHY_Atm_CP_ALLOC
