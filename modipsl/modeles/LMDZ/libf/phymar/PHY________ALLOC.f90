      subroutine PHY________ALLOC

!------------------------------------------------------------------------------+
!                                                         Mon 17-May-2013  MAR |
!                                                                              |
!     subroutine PHY________ALLOC  allocates prognostic variables of           |
!            MAT PHYsics                                                       |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 26-Feb-2013      |
!           Last Modification by H. Gallee,               Mon 17-May-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


      use Mod_Real
      use Mod_PHY____grd
      use Mod_PHY____kkl
      use Mod_PHY_S0_grd


      IMPLICIT NONE




! =================================
! ALLOCATION Mod_PHY____grd - BEGIN
! =================================

      allocate      ( lat__r(kcolp) )                     !     Latitude                    [radian]
      allocate      ( sinLat(kcolp) )                     ! sin(Latitude)                        [-]
      allocate      ( cosLat(kcolp) )                     ! cos(Latitude)                        [-]
      allocate      ( lon__r(kcolp) )                     !     Longitude                   [radian]
      allocate      ( lon__h(kcolp) )                     !     Longitude                     [hour]

      allocate      ( k1m(mzp) )                          ! k - 1
      allocate      ( k1p(mzp) )                          ! k + 1
      allocate      ( k2m(mzp) )                          ! k - 2

      allocate      (  sigma(mzpp) )                      ! Normalized Pressure (Vertical Coordinate)
      allocate      (  sigmi(mzpp) )                      ! (sigma(k    )+sigma(k-1  )) / 2
      allocate      ( dsigma(mzp) )                       !  sigma(k+1  )-sigma(k    )      
      allocate      ( dsigmi(mzp) )                       !  sigma(k+1/2)-sigma(k-1/2)
      allocate      ( hsigma(mzp) )                       ! Height of atmospheric layers      [magl]

      allocate      ( ii__AP(kcolp) )                     ! WORK   point   i Coordinate
      allocate      ( jj__AP(kcolp) )                     ! WORK   point   j Coordinate
      allocate      ( ikl_AP(ixp1:mxpp,jyp1:mypp) )       ! WORK   point vec Coordinate


! =================================
! ALLOCATION Mod_PHY____grd -   END
! =================================




! =================================
! ALLOCATION Mod_PHY____kkl - BEGIN
! =================================

      allocate      ( sh__AP(kcolp) )                     ! Topography                           [m]
      allocate      ( sha_AP(kcolp) )                     ! Topography Anomaly                   [m]
      allocate      ( slopAP(kcolp) )                     ! Topography Slope                     [-]
      allocate      ( sloxAP(kcolp) )                     ! Topography Slope, x-direction        [-]
      allocate      ( sloyAP(kcolp) )                     ! Topography Slope, y-direction        [-]
      allocate      ( MMskAP(kcolp,n_azim) )              ! Mountain   Mask                      [-]



! =================================
! ALLOCATION Mod_PHY____kkl -   END
! =================================




      return
      end subroutine PHY________ALLOC
