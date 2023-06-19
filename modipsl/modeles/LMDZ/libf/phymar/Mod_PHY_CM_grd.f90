      module Mod_PHY_CM_grd

!--------------------------------------------------------------------------+
!                                                     Sun 30-Jun-2013  MAR |
!     module Mod_PHY_CM_grd contains the main grid descriptor        of    |
!                Cloud MicroPhysical          Scheme                       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Fri 22-Mar-2013      |
!           Last Modification by H. Gallee,           Sun 30-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



      real(kind=real8), SAVE                        ::  dt__CM          ! Time Step of                           Cloud MicroPhysical Scheme   [s]
      integer, SAVE                                 ::  jt__CM          ! Number of DYnamical Time Steps for one Cloud Microphysics Time Step [-]

      integer, parameter                            ::  mz1_CM = 2      ! Top       of Clouds     Vertical Index



      end module Mod_PHY_CM_grd
