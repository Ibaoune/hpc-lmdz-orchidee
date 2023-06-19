
      subroutine PHY_Atm_S0_ALLOC

!------------------------------------------------------------------------------+
!                                                         Sat 15-Jun-2013  MAR |
!                                                                              |
!     subroutine PHY_Atm_S0_ALLOC  allocates prognostic variables of           |
!                Insolation Scheme used by MAR                                 |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Thu 25-Apr-2013      |
!           Last Modification by H. Gallee,               Sat 15-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


      use Mod_Real
      use Mod_PHY____grd
      use Mod_PHY_S0_grd
      use Mod_PHY_S0_kkl


      IMPLICIT NONE



! =================================
! ALLOCATION Mod_PHY_S0_kkl - BEGIN
! =================================

      allocate  ( csz0S0   (kcolp) )                         !  cosine (solar zenithal Distance)                          [-]
      allocate  ( csz_S0   (kcolp) )                         !  cosine (solar zenithal Distance), including Slope  Effect [-]
      allocate  ( cszkS0   (kcolp, n_azim) )                 !  cosine (solar zenithal Distance), including Mountain Mask [-]
      allocate  ( slopS0   (kcolp) )                         !  Cosine of Fall Line Angle                                 [-]
      allocate  ( omenS0   (kcolp) )                         !  Fall Line Azimuth (Downslope Direction)                   [-]

! =================================
! ALLOCATION Mod_PHY_S0_kkl -   END
! =================================


      end subroutine PHY_Atm_S0_ALLOC
