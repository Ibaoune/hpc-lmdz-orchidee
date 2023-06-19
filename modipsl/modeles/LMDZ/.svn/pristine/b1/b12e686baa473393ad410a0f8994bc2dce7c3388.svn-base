
      subroutine PHY_Atm_CM_ALLOC

!------------------------------------------------------------------------------+
!                                                         Sun  9-Jun-2013  MAR |
!                                                                              |
!     subroutine PHY_Atm_CM_ALLOC  allocates prognostic variables of           |
!                Cloud Microphysical    Scheme used by MAR                     |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 19-Mar-2013      |
!           Last Modification by H. Gallee,               Sun  9-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


      use Mod_Real
      use Mod_PHY____grd
      use Mod_PHY_CM_kkl


      IMPLICIT NONE



! =================================
! ALLOCATION Mod_PHY_CM_kkl - BEGIN
! =================================

      allocate  ( Ta__CM(kcolp,mzp ) )                       !  Air   Temperature                                          [K]
      allocate  ( qvswCM(kcolp,mzpp) )                       !  Saturation Specific Humidity     (over liquid water)   [kg/kg]
      allocate  ( qvsiCM(kcolp,mzpp) )                       !  Saturation Specific Humidity     (over ice)            [kg/kg]
      allocate  ( qw__CM(kcolp,mzpp) )                       !  Cloud Droplets      Concentration                      [kg/kg]
      allocate  ( qwd_CM(kcolp,mzp ) )                       !  Cloud Droplets      Concentration Variation            [kg/kg]
      allocate  ( CCNwCM(kcolp,mzp ) )                       !  Cloud Droplets      Number                              [-/m3]
      allocate  ( qi__CM(kcolp,mzpp) )                       !  Cloud Ice Particles Concentration                      [kg/kg]
      allocate  ( qid_CM(kcolp,mzp ) )                       !  Cloud Ice Particles Concentration Variation            [kg/kg]
      allocate  ( CCNiCM(kcolp,mzp ) )                       !  Cloud Ice Particles Number                              [-/m3]
      allocate  ( CFraCM(kcolp,mzp ) )                       !  Cloud               Fraction                            [-]   
      allocate  ( qs__CM(kcolp,mzpp) )                       !  Snow      Particles Concentration                      [kg/kg]
! #qg allocate  ( qg__CM(kcolp,mzpp) )                       !  Graupels            Concentration                      [kg/kg]
      allocate  ( qr__CM(kcolp,mzpp) )                       !  Rain  Drops         Concentration                      [kg/kg]
      allocate  ( HLatCM(kcolp,mzp ) )                       !  Latent Heat Release                                     [W/m2]

      allocate  ( uss_CM(kcolp) )                            !  Snow  Particles     Turbulent Surface Flux            [kg m/s]
      allocate  ( Ice0CM(kcolp) )                            !  Ice C.Accumulation (time t-dt)                        [m w.e.]
      allocate  ( ICE_CM(kcolp) )                            !  Ice C.Accumulation (time t   )                        [m w.e.]
      allocate  ( Sno0CM(kcolp) )                            !  Snow  Accumulation (time t-dt, before snow erosion)   [m w.e.]
      allocate  ( SnobCM(kcolp) )                            !  Snow  Accumulation (time t-dt, after  snow erosion)   [m w.e.]
      allocate  ( SnowCM(kcolp) )                            !  Snow  Accumulation (time t   )                        [m w.e.]
      allocate  ( Rai0CM(kcolp) )                            !  Rain  Accumulation (time t-dt)                        [m w.e.]
      allocate  ( RainCM(kcolp) )                            !  Rain  Accumulation (time t   )                        [m w.e.]

      allocate  ( dpktCM(kcolp,mzp ) )                       !  Reduced Potential Temperature TENDENCY                  [KX/s]
      allocate  ( dqv_CM(kcolp,mzp ) )                       !  Specific          Humidity    TENDENCY               [kg/kg/s]
      allocate  ( dqw_CM(kcolp,mzp ) )                       !  Cloud Droplets Concentration  TENDENCY               [kg/kg/s]
      allocate  ( dqi_CM(kcolp,mzp ) )                       !  Cloud Crystals Concentration  TENDENCY               [kg/kg/s]
      allocate  ( dqs_CM(kcolp,mzp ) )                       !  Snow Particles Concentration  TENDENCY               [kg/kg/s]
      allocate  ( dqr_CM(kcolp,mzp ) )                       !  Rain Drops     Concentration  TENDENCY               [kg/kg/s]
! #cw allocate  ( dCw_CM(kcolp,mzp ) )                       !  CCN            Concentration  TENDENCY                   [1/s]
      allocate  ( dCi_CM(kcolp,mzp ) )                       !  CIN            Concentration  TENDENCY                   [1/s]
      allocate  ( dCF_CM(kcolp,mzp ) )                       !  Cloud Fraction                TENDENCY                   [1/s]

      allocate  ( wat0EW(kcolp) )                            !  Total Precipitable  Water  in the  Air Column         [m w.e.]
      allocate  ( wat1EW(kcolp) )                            !  Total Precipitable  Water  in the  Air Column         [m w.e.]
      allocate  ( wat2EW(kcolp) )                            !  Total Precipitable  Water  in the  Air Column         [m w.e.]
      allocate  ( watfEW(kcolp) )                            !  Water Flux (Atm. --> Srf.) during 1 Time Step         [m w.e.]
      allocate  ( enr0EW(kcolp) )                            !  Total Energy (Sens. +Lat.) in the  Air Column         [m w.e.]
      allocate  ( enr1EW(kcolp) )                            !  Total Energy (Sens. +Lat.) in the  Air Column         [m w.e.]
      allocate  ( enr2EW(kcolp) )                            !  Total Energy (Sens. +Lat.) in the  Air Column         [m w.e.]
      allocate  ( mphyEW(kcolp) )                            !                                                        [m w.e.]


!  Isotopes Proxies
!  ~~~~~~~~~~~~~~~~
      allocate  ( Hcd_CM(kcolp) )                            ! latent heat release                                   [mm w.e.]
      allocate  ( Tcd_CM(kcolp) )                            ! latent heat release weighted Air Temperature                [K]
      allocate  ( Zcd_CM(kcolp) )                            ! latent heat release weighted Altitude                       [m]
      allocate  ( Hsb_CM(kcolp) )                            ! latent heat absorb.                                   [mm w.e.]
      allocate  ( Tsb_CM(kcolp) )                            ! latent heat absorb. weighted Air Temperature                [K]
      allocate  ( Zsb_CM(kcolp) )                            ! latent heat absorb. weighted Altitude                       [m]



! =================================
! ALLOCATION Mod_PHY_CM_kkl -   END
! =================================



      end subroutine PHY_Atm_CM_ALLOC
