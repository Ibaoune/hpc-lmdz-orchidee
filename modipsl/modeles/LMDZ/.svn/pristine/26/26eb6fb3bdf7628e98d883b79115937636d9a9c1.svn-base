
      module Mod_PHY_AT_ctr


!--------------------------------------------------------------------------+
!                                                     Sat 25-May-2013  MAR |
!     module Mod_PHY_AT_ctr contains the logical variables used to set up  |
!            MAR PHYsics                                                   |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Mon  1-Apr-2013      |
!                             as Mod_PHY____ctr                            |
!           Last Modification by H. Gallee,           Sat 25-May-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE


! E-e and K-l models switches
! ---------------------------

      logical            ::   AT_TKE                  ! Set-Up of Turbulent Kinetic Energy Computation

      logical            ::   Ee_Duynkerke            ! Dunkerke           (1988) E-epsilon model of turbulence
      logical            ::   Ee_Kitada               ! Kitada             (1987) E-epsilon model of turbulence
      logical            ::   Ee_HuangRamn            ! Huang and Raman    (1991) E-epsilon model of turbulence
      logical            ::   Kl_TherryLac            ! Therry & Lacarrere (1983) K-l       model of turbulence


! Box weighted vertical moving Averages of TKE and eps
! ----------------------------------------------------

      logical            ::   AT_vav = .FALSE.        ! Box weighted vertical moving Averages of TKE & eps
      integer, SAVE            ::   mz__KA =  1             ! Level No  below which moving Averages of TKE and eps is effective
      real(kind=real8), SAVE   ::   zz__KA = 20.            ! Height    below which moving Averages of TKE and eps is effective


! Flux Limitor
! ------------

      logical            ::   AT_LIM = .FALSE.        ! pkt  Turbulent Flux Limitor
      real(kind=real8), SAVE   ::   DThMAX =  5.            ! Maximum allowed d(PT ) in 1 d(t)                                    [K/h]
      real(kind=real8), SAVE   ::   Dpkt_X                  ! Maximum allowed d(pkt) in 1 d(t) = DThMAX /100.**(R/Cp) /3600.      [K/s]


      end module Mod_PHY_AT_ctr
