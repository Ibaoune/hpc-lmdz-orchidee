
      subroutine PHY_Atm_RT_ALLOC

!------------------------------------------------------------------------------+
!                                                         Sun 16-Jun-2013  MAR |
!                                                                              |
!     subroutine PHY_Atm_RT_ALLOC  allocates prognostic variables of           |
!                Radiative Transfer Scheme used by MAR                         |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Wed  6-Mar-2013      |
!           Last Modification by H. Gallee,               Sun 16-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


      use Mod_Real
      use Mod_PHY____grd
      use Mod_PHY_RT_grd
      use Mod_PHY_RT_kkl


      IMPLICIT NONE



! =================================
! ALLOCATION Mod_PHY_RT_kkl - BEGIN
! =================================

      allocate  ( O3__RT   (kcolp,mzp      ) )                 !  Ozone    Concentration                               [Pa/Pa]
      allocate  ( AersRT   (kcolp,mzp,naero) )                 !  Aerosol  Optical Thickness                               [-]
      allocate  ( ODAzRT   (kcolp,mzp      ) )                 !  Aerosols Optical Thickness (Layer z)                     [-]
      allocate  ( ODA_RT   (kcolp          ) )                 !  Aerosols Optical Thickness (vertically integrated)       [-]

      allocate  ( ODCzRT   (kcolp,mzp      ) )                 !  Clouds   Optical Thickness (Layer z)                     [-]
      allocate  ( ODC_RT   (kcolp          ) )                 !  Clouds   Optical Thickness (vertically integrated)       [-]

      allocate  ( FIRn_c   (kcolp,mzpp     ) )                 !  CLEAR-SKY         LW NET      FLUXES                  [W/m2]
      allocate  ( FIRn_t   (kcolp,mzpp     ) )                 !  TOTAL             LW NET      FLUXES                  [W/m2]
      allocate  ( FSOn_c   (kcolp,mzpp     ) )                 !  CLEAR-SKY         SW NET      FLUXES                  [W/m2]
      allocate  ( FSOn_t   (kcolp,mzpp     ) )                 !  TOTAL             SW NET      FLUXES                  [W/m2]
      allocate  ( FSOs_t   (kcolp          ) )                 !  TOTAL-SKY SURFACE SW DOWNWARD FLUX                    [W/m2]
      allocate  ( FSOdir   (kcolp          ) )                 !  SOLAR RADIANCE  IN SUN'S  DIRECTION                   [W/m2]
      allocate  ( FSOsUV   (kcolp          ) )                 !  SURFACE   DOWNWARD U.V.   RADIATION                   [W/m2]
      allocate  ( FSOeff   (kcolp          ) )                 !  PHOTOSYNTHETICALLY ACTIVE RADIATION                   [W/m2]

      allocate  ( SWDsRT   (kcolp          ) )                 ! Surface ShrtWave Heat Flux (+)  (Downward)             [W/m2]
      allocate  ( SWAsRT   (kcolp          ) )                 ! Surface ShrtWave Heat Flux (+)  (Absorbed)             [W/m2]
      allocate  ( LWDsRT   (kcolp          ) )                 ! Surface LongWave Heat Flux (+)  (Downward)             [W/m2]
      allocate  ( LWUsRT   (kcolp          ) )                 ! Surface LongWave Heat Flux (+)  (  Upward)             [W/m2]
      allocate  ( ClouRT   (kcolp          ) )                 ! Total Cloudiness above lowest Atmospheric Level           [-]

      allocate  ( OLR_RT   (kcolp          ) )                 !  OutgoingLongWave Radiation (+)  (  Upward)            [W/m2]

      allocate  ( SWdTRT   (kcolp,mzp      ) )                 !  Radiative Heating SW                                 [K/Day]
      allocate  ( LWdTRT   (kcolp,mzp      ) )                 !  Radiative Heating LW                                 [K/Day]
      allocate  ( dpktRT   (kcolp,mzp      ) )                 !  Radiative Heating SW + LW                              [K/s]


! =================================
! ALLOCATION Mod_PHY_RT_kkl -   END
! =================================


      end subroutine PHY_Atm_RT_ALLOC
