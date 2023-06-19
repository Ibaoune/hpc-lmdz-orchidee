      module Mod_SISVAT_loc

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVAT_loc contains the         diagnostic  variables of  |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Sat  8-Feb-2013      |
!                    modified by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+

      use    Mod_Real


      IMPLICIT NONE


      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  NLaysv  ! New   Snow     Layer   Switch
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  i_thin  ! Index of the thinest Layer
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)    ::  LIndsv  ! Contiguous Layer relative Index

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  albisv  ! Integrated Surface Albedo
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  albssv  ! Soil               Albedo [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  SoCasv  ! Canopy  Absorbed Solar Radiat.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  SoSosv  ! Surface Absorbed Solar Radiat.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  IRv_sv  ! Vegetation IR Flux  [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Evg_sv  ! Emissivity of Vegetation+Snow
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Eso_sv  ! Emissivity of       Soil+Snow
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  tau_sv  ! Transmited Radiation Fraction
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  rrMxsv  ! Canopy Maximum Intercepted Rain
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  LAIesv  ! effective LAI for transpirati.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  LAI_sv  ! corrected LAI in case of snow
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  glf_sv  ! Green  Leaf Fraction
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Sigmsv  ! Canopy Ventilation  Factor
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  HSv_sv  ! Sensible Heat Flux  [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  HLv_sv  ! Latent   Heat Flux  [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  HSs_sv  ! Sensible Heat Flux (t)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  HLs_sv  ! Latent   Heat Flux (t)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  sqrCm0  ! in Neutral Drag Coef.Moment.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  sqrCh0  ! in Neutral Drag Coef.Heat
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Lx_H2O  ! Latent Heat of Vaporiz./Sublim.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  ram_sv  ! Aerodyn.Resistance (Moment.)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  rah_sv  ! Aerodyn.Resistance (Heat)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Fh__sv  ! Stability Function
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dFh_sv  ! Stability Function (Deriv.)
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Evp_sv  ! Evaporation        [kg/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  EvT_sv  ! Evapotranspiration [kg/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  LSdzsv  ! Land/Sea Vert. Discretiz. Fact.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  Tsrfsv  ! Surface    Temperature
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  sEX_sv  ! Verticaly Integr.Extinct.Coef.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  zzsnsv  ! Snow  Pack Thickness      [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  psi_sv  ! Soil   Water        Potential
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  Khydsv  ! Soil   Hydraulic    Conductiv.
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  Rootsv  ! Root Water Pump      [kg/m2/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  EExcsv  ! Energy in Excess, current


      end module Mod_SISVAT_loc
