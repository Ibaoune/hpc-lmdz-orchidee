      module Mod_PHY_S0_ctr

!--------------------------------------------------------------------------+
!                                                     Sat 15-Jun-2013  MAR |
!     module Mod_PHY_S0_ctr contains the main Switches               of    |
!                Insolation                   Scheme                       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Thu 25-Apr-2013      |
!           Last Modification by H. Gallee,           Sat 15-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



      logical                                       ::  FaceS0            !  Slope Influence on Insolation
      logical                                       ::  MMskS0            !  Mount MASKing   of Insolation
      logical                                       ::  LDayS0 = .TRUE.   !  Day   Length    at Domain Origin



      end module Mod_PHY_S0_ctr
