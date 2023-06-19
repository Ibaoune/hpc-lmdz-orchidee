      module Mod_PHY_AT_grd

!--------------------------------------------------------------------------+
!                                                     Tue 30-Apr-2013  MAR |
!     module Mod_PHY_AT_grd contains the main grid descriptor        of    |
!                Turbulent Vertical Diffusion Scheme                       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Sun 17-Mar-2013      |
!           Last Modification by H. Gallee,           Tue 30-Apr-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



      integer, SAVE                                       ::  jt__AT  ! Number of DYnamical Time Steps for one Vertical Diffusion Time Step
      real(kind=real8), SAVE                              ::  dt__AT  ! Time Step          of Atmospheric Turbulence
      real(kind=real8), SAVE                              ::  TimeAT  ! Time previous call of Atmospheric Turbulence  [HOURS since 1901-01-15 00:00:00]
      real(kind=real8), SAVE                              ::  alphAT  ! Explicitness
      real(kind=real8), SAVE                              ::  betaAT  ! Implicitness
      real(kind=real8), SAVE                              ::  a_b_AT  !
      character(len=1)                              ::  schmAT  ! m = momentum / h = scalar / e = ekman spiral
      logical                                       ::  NewAAT  ! Diagonal Time Depend. Coeff. Calculation Switch



      end module Mod_PHY_AT_grd
