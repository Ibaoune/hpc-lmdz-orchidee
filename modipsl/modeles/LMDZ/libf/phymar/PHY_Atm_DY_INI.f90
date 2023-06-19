      subroutine PHY_Atm_DY_INI

!------------------------------------------------------------------------------+
!                                                         Fri  7-Jun-2013  MAR |
!   MAR          PHY_Atm_DY_INI                                                |
!     subroutine PHY_Atm_DY_INI intializes Turbulent Vertical Diffusion Scheme |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,               Fri  7-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
!     use Mod_PHY____grd


      IMPLICIT NONE



!     integer    ::    i,     j,     k,     ikl



! Initialization of Mod_PHY_DY_kkl
! ================================

!                ****************
          CALL   PHY_Atm_DY_ALLOC
!                ****************


      end subroutine PHY_Atm_DY_INI
