      module Mod_SISVATLTVg

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVATLTVg contains local variables of SISVAT_TVg         |
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



! OUTPUT/Verification: Energy/Water Budget
! #e1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  ETVg_d          ! VegetationPower, Forcing

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  Tveg_0          ! Canopy Temperature, Previous t
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  dIRdTv          ! InfraRed  NET(t), Derivative(t)
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  dHSdTv          ! Sensible Heat FL. Derivative(t)
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  dHLdTv          ! Latent   Heat FL. Derivative(t)
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  dEvpdT          ! Evapo(transpi)ration Derivative
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  dEvTdT          ! Evapo(transpi)ration Derivative


      end module Mod_SISVATLTVg
