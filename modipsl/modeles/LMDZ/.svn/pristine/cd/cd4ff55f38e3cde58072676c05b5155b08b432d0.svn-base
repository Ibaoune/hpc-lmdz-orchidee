      module Mod_PHY_CP_kkl

!--------------------------------------------------------------------------+
!                                                     Thu  9-May-2013  MAR |
!     module Mod_PHY_CP contains the representation in vector frame of     |
!                       Variables of the Convective Mass Flux Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  9-Apr-2013      |
!           Last Modification by H. Gallee,           Thu  9-May-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE


! Convective Mass Flux Scheme INPUT        Variables
! --------------------------------------------------



! Convective Mass Flux Scheme INPUT/OUTPUT Variables
! --------------------------------------------------

      real(kind=real8), SAVE                                ::  rANA         ! Subgrid Mountain Breeze: Horizontal Divergence           [1/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  hANA         ! D("Subgrid Mountain" Height - "Resolved Mountain" Height)  [m]



! Convective Mass Flux Scheme OUTPUT       Variables
! --------------------------------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  timeCP       ! Convective Adjustment Time                                 [s]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dpktCP       ! Reduced  Potential Temperature Tendency                [K/X/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dqv_CP       ! Specific Humidity              Tendency              [kg/kg/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dqw_CP       ! Cloud Droplets Concentration   Tendency              [kg/kg/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  dqi_CP       ! Cloud Crystals Concentration   Tendency              [kg/kg/s]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  drr_CP       ! Rain (convective)              Tendency                  [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  rainCP       ! Rain (convective)                                          [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  dss_CP       ! Snow (convective)              Tendency                  [m/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  snowCP       ! Snow (convective)                                          [m]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  CAPECP       ! Convective Available Potentential Energy             [.......]


      end module Mod_PHY_CP_kkl
