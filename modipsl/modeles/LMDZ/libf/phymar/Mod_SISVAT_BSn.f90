     module Mod_SISVAT_BSn

!--------------------------------------------------------------------------+
!    module Mod_SISVAT_BSn                            Fri  2-Feb-2013  MAR |
!    module Mod_SISVAT_BSn  contains specific variables (and constants)    |
!           used by Soil/Ice Snow Vegetation Atmosphere Transfer Scheme    |
!                                                                          |
!--------------------------------------------------------------------------+



! General Variables
! =================

      use Mod_Real
      use Mod_SISVAT_dim


      IMPLICIT NONE


! SISVAT_BSn    specific variables (and constants)
! ================================================

      logical           ::  BlowIn = .false.

      real(kind=real8), SAVE  ::  FacSBS,FacUBS                 !
      real(kind=real8), SAVE  ::  Por_BS                        ! Snow       Porosity
      real(kind=real8), SAVE  ::  SheaBS                        !
      real(kind=real8), SAVE  ::  rCd10n                        ! GM97:   assumed neutral stabil.


     end module Mod_SISVAT_BSn
