      module Mod_PHY_S0_grd

!--------------------------------------------------------------------------+
!                                                     Fri 14-Jun-2013  MAR |
!     module Mod_PHY_S0_grd contains the main grid descriptor        of    |
!                Insolation                   Scheme                       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Thu 25-Apr-2013      |
!           Last Modification by H. Gallee,           Fri 14-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



      real(kind=real8), SAVE                              ::  dST_UA          ! Distance Soleil-Terre                      [u.a.] 

      integer, SAVE                                       ::  n_azim          ! Nb of azimuth to look for Mountain Mask       [-]
      real(kind=real8), SAVE                              ::  d_azim          ! azimuth interval for  one Mountain Mask  [radian]



      end module Mod_PHY_S0_grd
