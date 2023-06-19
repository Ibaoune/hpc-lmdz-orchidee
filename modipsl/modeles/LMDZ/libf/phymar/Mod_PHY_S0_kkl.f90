      module Mod_PHY_S0_kkl

!--------------------------------------------------------------------------+
!                                                     Fri 26-Apr-2013  MAR |
!     module Mod_PHY_S0 contains the representation in vector   frame of   |
!                       Variables of the Insolation Scheme                 |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Thu 25-Apr-2013      |
!           Last Modification by H. Gallee,           Fri 26-Apr-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE


! Insolation        OUTPUT        Variables
! -----------------------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  csz0S0       ! cosine (solar zenithal Distance)                          [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  csz_S0       ! cosine (solar zenithal Distance), including Slope  Effect [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  cszkS0       ! cosine (solar zenithal Distance), including Mountain Mask [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  slopS0       ! Cosine of Fall Line Angle                                 [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)      ::  omenS0       ! Fall Line Azimuth (Downslope Direction)                   [-]


      end module Mod_PHY_S0_kkl
