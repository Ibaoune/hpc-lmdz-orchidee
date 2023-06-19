     module Mod_SISVAT_BDu

!--------------------------------------------------------------------------+
!    module Mod_SISVAT_BDu                            Fri  2-Feb-2013  MAR |
!    module Mod_SISVAT_BDu  contains specific variables (and constants)    |
!           used by Soil/Ice Snow Vegetation Atmosphere Transfer Scheme    |
!                                                                          |
!--------------------------------------------------------------------------+



! General Variables
! =================

      use Mod_Real
      use Mod_SISVAT_dim
      use Mod_SISVAT_dat


      IMPLICIT NONE


! SISVAT_BDu    specific variables (and constants)
! ================================================

      logical           ::  logust = .false.

      real(kind=real8), SAVE  ::  etaust(0:nsot)



     end module Mod_SISVAT_BDu
