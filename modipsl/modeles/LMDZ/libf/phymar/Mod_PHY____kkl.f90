      module Mod_PHY____kkl

!--------------------------------------------------------------------------+
!                                                     Fri 14-Jun-2013  MAR |
!     module Mod_PHY____kkl contains the main              variables of    |
!            MAR PHYsical   Parameterizations                              |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,           Fri 14-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



! PHY    INPUT        Variables
! -----------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    :: sh__AP        !  Topography          [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    :: sha_AP        !  Topography Anomaly  [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    :: slopAP        !  Topography Slope    [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    :: sloxAP        !  Topography Slope, x [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    :: sloyAP        !  Topography Slope, y [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  :: MMskAP        !  Mountain   Mask     [-]


! PHY    INPUT/OUTPUT Variables
! -----------------------------



! PHY    OUTPUT       Variables
! -----------------------------


      end module Mod_PHY____kkl
