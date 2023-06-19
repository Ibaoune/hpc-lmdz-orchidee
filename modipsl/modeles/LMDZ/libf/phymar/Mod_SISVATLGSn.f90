      module Mod_SISVATLGSn

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVATLGSn contains local variables of SISVAT_GSn         |
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



      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  ro_dry        ! Dry Density            [g/cm3]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  etaSno        ! Liquid Water Content   [g/cm2]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  SnMass        ! Snow   Mass            [kg/m2]


      end module Mod_SISVATLGSn
