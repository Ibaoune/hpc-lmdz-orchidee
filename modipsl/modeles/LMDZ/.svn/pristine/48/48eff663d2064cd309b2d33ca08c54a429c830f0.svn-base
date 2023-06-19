      module Mod_PHY____grd


!--------------------------------------------------------------------------+
!                                                     Mon 17-Jun-2013  MAR |
!     module Mod_PHY____grd contains the characteristics of the grid for   |
!            MAR PHYsics                                                   |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Sat 16-Feb-2013      |
!           Last Modification by H. Gallee,           Mon 17-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



      logical :: FlagDALLOC = .FALSE.

      integer, SAVE :: YearTU   !
      integer, SAVE :: Mon_TU   !
      integer, SAVE :: Day_TU   !
      integer, SAVE :: HourTU   ! Hour, Universal Time
      integer, SAVE :: MinuTU   !
      integer, SAVE :: Sec_TU   !

      integer, SAVE :: it_EXP   ! Nb of iterations    since the beginning of the EXPeriment
      integer, SAVE :: it_RUN   ! Nb of iterations    since the beginning of the RUN (job)

      integer, SAVE :: mxp      ! Nb of interior      Grid Points, x-Direction
      integer, SAVE :: mxpp     ! Nb of interior      Grid Points, x-Direction + 1
      integer, SAVE :: myp      ! Nb of interior      Grid Points, y-Direction
      integer, SAVE :: mypp     ! Nb of interior      Grid Points, y-Direction + 1
      integer, SAVE :: ixp1     ! 1er pt en x de la grille dynamique utile dans grille physique
      integer, SAVE :: jyp1     ! 1er pt en y de la grille dynamique utile dans grille physique
      integer, SAVE :: kcolp    ! Nb of interior      Vertical Columns (mxp * myp)        
      integer, SAVE :: mzp      ! Nb of Atmospheric   Levels
      integer, SAVE :: mzpp     ! Nb of Atmospheric   Levels                   + 1

      integer, SAVE                                        :: i_x0
      integer, SAVE                                        :: j_y0
      integer, SAVE                                        :: ikl0

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)    :: lat__r      !     Latitude                    [radian]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)    :: sinLat      ! sin(Latitude)                        [-]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)    :: cosLat      ! cos(Latitude)                        [-]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)    :: lon__r      !     Longitude                   [radian]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)    :: lon__h      !     Longitude                     [hour]

      real(kind=real8), SAVE                               :: timeTU      ! Time   [HOURS since 1901-01-15 00:00:00]
      real(kind=real8), SAVE                               :: dxHOST      ! dx
      real(kind=real8), SAVE                               :: dx2inv      ! 1 / (2 dx)
      real(kind=real8), SAVE                               :: dy2inv      ! 1 / (2 dy)
      real(kind=real8), SAVE                               :: pt__DY      ! Model Pressure Top                 [kPa]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)    ::  sigma      ! Vertical Coord. (normalized Pressure)
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)    ::  sigmi      !(sigma(k-1  )+sigma(k    )) / 2 
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)    :: dsigma      ! sigma(k+1  )-sigma(k    )
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)    :: dsigmi      ! sigma(k+1/2)-sigma(k-1/2)
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)    :: hsigma      ! Height of atmospheric layers      [magl]

      integer, SAVE,          ALLOCATABLE, dimension(:)    :: k1m         ! k - 1
      integer, SAVE,          ALLOCATABLE, dimension(:)    :: k1p         ! k + 1
      integer, SAVE,          ALLOCATABLE, dimension(:)    :: k2m         ! k - 2

      integer, SAVE,          ALLOCATABLE, dimension(:)    :: ii__AP      ! WORK   point   i Coordinate
      integer, SAVE,          ALLOCATABLE, dimension(:)    :: jj__AP      ! WORK   point   i Coordinate
      integer, SAVE,          ALLOCATABLE, dimension(:,:)  :: ikl_AP      ! WORK   point vec Coordinate



      end module Mod_PHY____grd
