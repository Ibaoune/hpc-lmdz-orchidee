      module Mod_PHY_RT_grd

!--------------------------------------------------------------------------+
!                                                     Tue 30-Apr-2013  MAR |
!     module Mod_PHY_RT_grd contains the main grid descriptor        of    |
!                Radiative Transfert          Scheme                       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  2-Apr-2013      |
!           Last Modification by H. Gallee,           Tue 30-Apr-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



      integer, SAVE                                       ::  naero           ! Nb of Aerosol Classes                                               [-]
      integer, SAVE                                       ::  jt__RT          ! Number of DYnamical Time Steps for one Radiative Transfer Time Step [-]
      real(kind=real8), SAVE                              ::  dt__RT          ! Time Step of Radiative Transfert                                    [s]



      end module Mod_PHY_RT_grd
