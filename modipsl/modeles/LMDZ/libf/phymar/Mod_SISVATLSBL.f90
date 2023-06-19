      module Mod_SISVATLSBL

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVATLSBL contains local variables of SISVATeSBL         |
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



      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  VVaSBL          ! effective SBL wind speed
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  dTa_Ts          ! effective SBL Temperature diff.
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  LMOmom          ! Monin-Obukhov Scale Momentum
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  CDm             ! Drag Coefficient, Momentum
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  CDs             ! Drag Coefficient, Blown **
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  :: rCDs             ! Drag Coefficient, Blown **
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  CDh             ! Drag Coefficient, Scalar
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  Richar          ! Richardson Number


      end module Mod_SISVATLSBL
