     module Mod_SISVAT_SIc

!--------------------------------------------------------------------------+
!    module Mod_SISVAT_SIc                            Fri  2-Feb-2013  MAR |
!    module Mod_SISVAT_SIc  contains specific variables (and constants)    |
!           used by Soil/Ice Snow Vegetation Atmosphere Transfer Scheme    |
!                                                                          |
!--------------------------------------------------------------------------+



! General Variables
! =================

      use Mod_Real
      use Mod_SISVAT_dim


      IMPLICIT NONE


! SISVAT_SIc    specific variables (and constants)
! ================================================

      logical           ::  SIcINI = .false.

      real(kind=real8), SAVE  ::  Crodzw
      real(kind=real8), SAVE  ::  Lro__I



     end module Mod_SISVAT_SIc
