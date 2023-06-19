      module Mod_SISVATLmmm

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVATLmmm contains local variables of SISVAT main        |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Mon 17-Jun-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


! Global  Variables
! =================

      use Mod_Real



! Internal Variables
! ==================

      IMPLICIT NONE


! Non Local
! ---------

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  TBr_sv          ! Brightness Temperature
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  IRdwsv          ! DOWNward   IR Flux
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  IRupsv          ! UPward     IR Flux
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  Bdzssv          ! Buffer Snow Layer Thickness
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  z_snsv          ! Snow-Ice, current Thickness

! Energy         Budget
! ~~~~~~~~~~~~~~~~~~~~~
! #e1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  ETVg_d          ! VegetationPower, Forcing
! #e1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  EqSn_0          ! Snow Energy, befor Phase Change
! #e1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  EqSn_1          ! Snow Energy, after Phase Change
! #e1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  EqSn_d          ! Energy in Excess

! OUTPUT/Verification: H2O    Conservation
! #m0 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  Wats_0          ! Soil Water,  before Forcing
! #m0 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  Wats_1          ! Soil Water,  after  Forcing
! #m0 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  Wats_d          ! Soil Water,         Forcing

! OUTPUT/Verification: * Mass Conservation
! #m1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  SIsubl          ! Snow Sublimed/Deposed Mass
! #m1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  SImelt          ! Snow Melted           Mass
! #m1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  SIrnof          ! Local Surficial Water + Run OFF

! OUTPUT/Verification: SeaIce Conservation
! #m2 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  SIvAcr          ! Sea-Ice      Vertical Acretion


! Local
! -----

      integer, SAVE         , ALLOCATABLE, dimension(:,:)  ::  IcIndx          ! No   Ice               Mask

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  FallOK          ! Snow Contribution to the Canopy

! Energy and Mass Budget
! ~~~~~~~~~~~~~~~~~~~~~~
! #e1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::    Enrsvd        ! Soil+Vegetat  Power  Forcing

                                                                         ! H2O    Conservation
! #m0 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::    Watsv0        ! Soil+Vegetat, before Forcing
! #m0 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::    Watsvd        ! Soil+Vegetat  Water  Forcing

                                                                         ! * Mass Conservation
! #m1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::    SIWm_0,SIWm_1 ! Snow Initial/Final        Mass
! #m1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::    SIWa_i,SIWa_f ! Snow Initial/Final ATM Forcing
! #m1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::    SIWe_i,SIWe_f ! Snow Initial/Final BLS Forcing


      end module Mod_SISVATLmmm
