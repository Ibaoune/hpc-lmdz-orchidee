      module Mod_SISVAT_gpt

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVAT_gpt contains the Grid Point           variables of |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Sat 22-Jun-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



! SISVAT INPUT        Variables
! -----------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  sst_SB      ! Ocean  FORCING (SST)                         [K]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  sif_SB      ! Ocean  FORCING (Sea-Ice Fraction)            [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  MaskSV_gpt  ! Land(1)-Sea(0) Mask         (Cell value)     [-]



! SISVAT INPUT/OUTPUT Variables
! -----------------------------




! SISVAT OUTPUT       Variables
! -----------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  Alb_SV_gpt  ! Surface  Albedo             (Cell average)   [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  EmisSV_gpt  ! LongWave Surface Emissivity (Cell average)   [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  Tas_SV_gpt  ! Surface  Air    Temperature (Cell average)   [K]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  HSenSV_gpt  ! Sensible Heat   Flux (+ => Downward)      [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  HLatSV_gpt  ! Latent   Heat   Flux (+ => Downward)      [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  LMO_SV_gpt  ! Obukhov  Length             (Cell average)   [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  us__SV_gpt  ! Friction Velocity                          [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  uts_SV_gpt  ! Sensible Heat   Flux Turbulent Scale     [K m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  uqs_SV_gpt  ! Latent   Heat   Flux Turbulent Scale [kg/kg m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  WE2aSV_gpt  ! Cumulative H2O  Flux from the Surface  [mm w.e.]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  hFraSV_gpt  ! Frazil   Thickness                           [m]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dpktSV_gpt  ! Reduced  Potential Temperature Tendency   [KX/s]



      end module Mod_SISVAT_gpt
