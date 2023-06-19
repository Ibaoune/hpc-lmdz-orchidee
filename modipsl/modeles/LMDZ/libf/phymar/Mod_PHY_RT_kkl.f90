      module Mod_PHY_RT_kkl

!--------------------------------------------------------------------------+
!                                                     Sun 16-Jun-2013  MAR |
!     module Mod_PHY_RT contains the representation in vector frame of     |
!                       Variables of the   Radiative Transfer Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Wed  3-Apr-2013      |
!           Last Modification by H. Gallee,           Sun 16-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE


! Radiative Transfer INPUT        Variables
! -----------------------------------------

!     real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  LWUsRT       ! Surface LongWave Heat Flux (+)    (Upward)             [W/m2]


! Radiative Transfer INPUT/OUTPUT Variables
! -----------------------------------------



! Radiative Transfer OUTPUT       Variables
! -----------------------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  O3__RT       ! Ozone   Concentration                                 [Pa/Pa]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  AersRT       ! Aerosol Optical Thickness                                 [-]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  ODC_RT       ! Clouds   Optical Depth (vertically integrated)            [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  ODCzRT       ! Clouds   Optical Depth (Layer z)                          [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  ODA_RT       ! Aerosols Optical Depth (vertically integrated)            [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  ODAzRT       ! Aerosols Optical Depth (Layer z)                          [-]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  FIRn_c       ! CLEAR-SKY         LW NET      FLUXES                   [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  FIRn_t       ! TOTAL             LW NET      FLUXES                   [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  FSOn_c       ! CLEAR-SKY         SW NET      FLUXES                   [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  FSOn_t       ! TOTAL             SW NET      FLUXES                   [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  FSOs_t       ! TOTAL-SKY SURFACE SW DOWNWARD FLUX                     [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  FSOdir       ! SOLAR RADIANCE  IN SUN'S  DIRECTION                    [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  FSOsUV       ! SURFACE   DOWNWARD U.V.   RADIATION                    [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  FSOeff       ! PHOTOSYNTHETICALLY ACTIVE RADIATION                    [W/m2]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  SWDsRT       ! Surface ShrtWave Heat Flux (+)  (Downward)             [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  SWAsRT       ! Surface ShrtWave Heat Flux (+)  (Absorbed)             [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  LWDsRT       ! Surface LongWave Heat Flux (+)  (Downward)             [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  LWUsRT       ! Surface LongWave Heat Flux (+)  (  Upward)             [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  ClouRT       ! Total Cloudiness above lowest Atmospheric Level           [-]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  OLR_RT       ! OutgoingLongWave Radiation (+)  (  Upward)             [W/m2] 

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  SWdTRT       ! Radiative Heating SW                                  [K/Day]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  LWdTRT       ! Radiative Heating LW                                  [K/Day]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dpktRT       ! Radiative Heating SW + LW                               [K/s]


      end module Mod_PHY_RT_kkl
