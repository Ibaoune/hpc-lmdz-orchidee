      module Mod_PHY_CM_kkl

!--------------------------------------------------------------------------+
!                                                     Sun  9-Jun-2013  MAR |
!     module Mod_PHY_CM_kkl contains the main (prognostic) variables of    |
!                Cloud Microphsics            Scheme                       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue 19-Mar-2013      |
!           Last Modification by H. Gallee,           Sun  9-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE


! Cloud Microphysical Model Parameters
! ------------------------------------

      logical   ::   Frac__Clouds  = .TRUE.   ! Fractional Cloud Cover



! PHY_CM INPUT        Variables
! -----------------------------



! PHY_CM INPUT/OUTPUT Variables
! -----------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  Ta__CM       !  Air   Temperature                                          [K]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  qvswCM       !  Saturation Specific Humidity     (over liquid water)   [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  qvsiCM       !  Saturation Specific Humidity     (over ice)            [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  qw__CM       !  Cloud Droplets      Concentration                      [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  qwd_CM       !  Cloud Droplets      Concentration Variation            [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  CCNwCM       !  Cloud Condensation  Nuclei                              [-/m3]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  qi__CM       !  Cloud Ice Particles Concentration                      [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  qid_CM       !  Cloud Ice Particles Concentration Variation            [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  CCNiCM       !  Cloud Ice Particles Number                              [-/m3]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  CFraCM       !  Cloud               Fraction                            [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  qs__CM       !  Snow      Particles Concentration                      [kg/kg]
! #qg real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  qg__CM       !  Graupels            Concentration                      [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  qr__CM       !  Rain  Drops         Concentration                      [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  HLatCM       !  Latent Heat Release                                     [W/m2]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  uss_CM       !  Snow  Particles     Turbulent Surface Flux            [kg m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  Ice0CM       !  Ice C.Accumulation (time t-dtR)                       [m w.e.]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  Ice_CM       !  Ice C.Accumulation (time t   )                        [m w.e.]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  Sno0CM       !  Snow  Accumulation (time t-dt, before snow erosion)   [m w.e.]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  SnobCM       !  Snow  Accumulation (time t-dt, after  snow erosion)   [m w.e.]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  SnowCM       !  Snow  Accumulation (time t   )                        [m w.e.]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  Rai0CM       !  Rain  Accumulation (time t-dt)                        [m w.e.]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  RainCM       !  Rain  Accumulation (time t   )                        [m w.e.]



! PHY_CM OUTPUT       Variables
! -----------------------------

      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:,:) ::  dpktCM       !  Reduced Potential Temperature TENDENCY                  [KX/s]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:,:) ::  dqv_CM       !  Specific          Humidity    TENDENCY               [kg/kg/s]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:,:) ::  dqw_CM       !  Cloud Droplets Concentration  TENDENCY               [kg/kg/s]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:,:) ::  dqi_CM       !  Cloud Crystals Concentration  TENDENCY               [kg/kg/s]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:,:) ::  dqs_CM       !  Snow Particles Concentration  TENDENCY               [kg/kg/s]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:,:) ::  dqr_CM       !  Rain Drops     Concentration  TENDENCY               [kg/kg/s]
! #cw real(kind=real8), SAVE, ALLOCATABLE ,dimension(:,:) ::  dCw_CM       !  CCN            Concentration  TENDENCY                   [1/s]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:,:) ::  dCi_CM       !  CIN            Concentration  TENDENCY                   [1/s]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:,:) ::  dCF_CM       !  Cloud Draction                TENDENCY                   [1/s]

      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:)   ::  wat0EW       !  Total Precipitable  Water  in the  Air Column         [m w.e.]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:)   ::  wat1EW       !  Total Precipitable  Water  in the  Air Column         [m w.e.]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:)   ::  wat2EW       !  Total Precipitable  Water  in the  Air Column         [m w.e.]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:)   ::  watfEW       !  Water Flux (Atm. --> Srf.) during 1 Time Step         [m w.e.]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:)   ::  enr0EW       !  Total Energy (Sens. +Lat.) in the  Air Column         [m w.e.]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:)   ::  enr1EW       !  Total Energy (Sens. +Lat.) in the  Air Column         [m w.e.]
      real(kind=real8), SAVE, ALLOCATABLE ,dimension(:)   ::  enr2EW       !  Total Energy (Sens. +Lat.) in the  Air Column         [m w.e.]
      character(len=20),ALLOCATABLE ,dimension(:)   ::  mphyEW       ! 
      character(len=20)                             ::  mauxEW       ! 


!  Isotopes Proxies
!  ~~~~~~~~~~~~~~~~
      logical                                         :: write_Proxy ! 
      integer, SAVE                                         :: nHL_CM      ! Counter
      real(kind=real8), SAVE,  ALLOCATABLE ,dimension(:)    :: Hcd_CM      ! latent heat release                                   [mm w.e.]
      real(kind=real8), SAVE,  ALLOCATABLE ,dimension(:)    :: Tcd_CM      ! latent heat release weighted Air Temperature                [K]
      real(kind=real8), SAVE,  ALLOCATABLE ,dimension(:)    :: Zcd_CM      ! latent heat release weighted Altitude                       [m]
      real(kind=real8), SAVE,  ALLOCATABLE ,dimension(:)    :: Hsb_CM      ! latent heat absorb.                                   [mm w.e.]
      real(kind=real8), SAVE,  ALLOCATABLE ,dimension(:)    :: Tsb_CM      ! latent heat absorb. weighted Air Temperature                [K]
      real(kind=real8), SAVE,  ALLOCATABLE ,dimension(:)    :: Zsb_CM      ! latent heat absorb. weighted Altitude                       [m]


      end module Mod_PHY_CM_kkl
