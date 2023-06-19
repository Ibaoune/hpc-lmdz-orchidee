      module Mod_PHY_CP_grd

!--------------------------------------------------------------------------+
!                                                     Tue 30-Apr-2013  MAR |
!     module Mod_PHY_CP_grd contains the main grid descriptor of           |
!                Convection Parameterization                               |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Mon  8-Apr-2013      |
!           Last Modification by H. Gallee,           Tue 30-Apr-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



      real(kind=real8), SAVE                              ::  dt__CP          ! Time Step of                           Convection Parameterization  [s]
      integer, SAVE                                       ::  jt__CP          ! Number of DYnamical Time Steps for one Convection Param.  Time Step [-]




      end module Mod_PHY_CP_grd
