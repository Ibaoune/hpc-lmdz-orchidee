      module Mod_SISVAT_TRV


!--------------------------------------------------------------------------+
!                                                     Thu 28-Feb-2013  MAR |
!     module Mod_SISVAT_TRV contains the constants of the                  |
!                radiative tranfer model through Vegetation                |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Wed 30-Jan-2013      |
!                    modified by H. Gallee,           Thu 28-Feb-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


! Global   Variables
! ==================

      use Mod_Real
      use Mod_SISVAT_dim


      IMPLICIT NONE



! Internal Variables
! ==================

      integer, SAVE  ::   ivg


      real(kind=real8), SAVE    ::  reVisL(0:nvgt)      ! Reflectivity  / Visible / Live Leaves
      real(kind=real8), SAVE    ::  renIRL(0:nvgt)      ! Reflectivity  / Near IR / Live Leaves
      real(kind=real8), SAVE    ::  trVisL(0:nvgt)      ! Transmitivity / Visible / Live Leaves
      real(kind=real8), SAVE    ::  trnIRL(0:nvgt)      ! Transmitivity / Near IR / Live Leaves
      real(kind=real8), SAVE    ::  reVisD(0:nvgt)      ! Reflectivity  / Visible / Dead Leaves
      real(kind=real8), SAVE    ::  renIRD(0:nvgt)      ! Reflectivity  / Near IR / Dead Leaves
      real(kind=real8), SAVE    ::  trVisD(0:nvgt)      ! Transmitivity / Visible / Dead Leaves
      real(kind=real8), SAVE    ::  trnIRD(0:nvgt)      ! Transmitivity / Near IR / Dead Leaves

      DATA (reVisL(ivg),renIRL(ivg),trVisL(ivg),trnIRL(ivg),           &
     &      reVisD(ivg),renIRD(ivg),trVisD(ivg),trnIRD(ivg),ivg=0,nvgt)&

!     reVisL renIRL trVisL trnIRL reVisD renIRD trVisD trnIRD    SVAT     CLASSES
!     ------ ------ ------ ------ ------ ------ ------ --------+ ----------------
     &/0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38, &!  0 NO VEGETATION
     & 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38, &!  1 CROPS LOW
     & 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38, &!  2 CROPS MEDIUM
     & 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38, &!  3 CROPS HIGH
     & 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38, &!  4 GRASS LOW
     & 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38, &!  5 GRASS MEDIUM
     & 0.11,  0.58,  0.07,  0.25,  0.36,  0.58,  0.22,  0.38, &!  6 GRASS HIGH
     & 0.10,  0.45,  0.05,  0.25,  0.16,  0.39,  0.01,  0.01, &!  7 BROADL LOW
     & 0.10,  0.45,  0.05,  0.25,  0.16,  0.39,  0.01,  0.01, &!  8 BROADL MEDIUM
     & 0.10,  0.45,  0.05,  0.25,  0.16,  0.39,  0.01,  0.01, &!  9 BROADL HIGH
     & 0.07,  0.35,  0.05,  0.10,  0.10,  0.39,  0.01,  0.01, &! 10 NEEDLE LOW
     & 0.07,  0.35,  0.05,  0.10,  0.10,  0.39,  0.01,  0.01, &! 11 NEEDLE MEDIUM
     & 0.07,  0.35,  0.05,  0.10,  0.10,  0.39,  0.01,  0.01/  ! 12 NEEDLE HIGH



      real(kind=real8), SAVE    ::  reVisS = 0.85       ! Reflectivity  / Visible / Canopy Snow
      real(kind=real8), SAVE    ::  renIRS = 0.85       ! Reflectivity  / Near IR / Canopy Snow
      real(kind=real8), SAVE    ::  trVisS = 0.00       ! Transmitivity / Visible / Canopy Snow
      real(kind=real8), SAVE    ::  trnIRS = 0.00       ! Transmitivity / Near IR / Canopy Snow
!     REMARK: Possible Refinement by taking actual Surface Snow Reflectivities
!     ^^^^^^

      real(kind=real8), SAVE    ::  snCaMx = 0.5        ! Canopy Snow Thickness for having Snow
                                                  ! Snow Reflectivity and Transmitivity
      real(kind=real8), SAVE    ::  CriStR = 25.        ! Critical Radiation Stomatal Resistance



      end module Mod_SISVAT_TRV
