
      subroutine PHY_Atm_AT_ALLOC

!------------------------------------------------------------------------------+
!                                                         Sun 16-Jun-2013  MAR |
!                                                                              |
!     subroutine PHY_Atm_AT_ALLOC  allocates prognostic variables of           |
!                Atmospheric Turbulence Scheme used by MAR                     |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,               Sun 16-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


      use Mod_Real
      use Mod_PHY____grd
      use Mod_PHY_AT_kkl


      IMPLICIT NONE



! =================================
! ALLOCATION Mod_PHY_AT_kkl - BEGIN
! =================================

      allocate  ( var_AT(kcolp,mzpp) )                       ! Dummy            to Diffuse                                [x]
      allocate  ( Ac0_AT(      mzp ) )                       ! Tridiagonal Matrix Coefficient A: Common Factor        [m2/s3]
      allocate  ( Cc0_AT(      mzp ) )                       ! Tridiagonal Matrix Coefficient C: Common Factor        [m2/s3]
      allocate  ( Kz0_AT(      mzp ) )                       ! Vertical  Turbulent Diffusion Coefficient (MINIMUM)     [m2/s]
      allocate  ( Ac__AT(kcolp,mzp ) )                       ! Tridiagonal Matrix Coefficient A: Common Factor (t)     [s/m2]
      allocate  ( Cc__AT(kcolp,mzp ) )                       ! Tridiagonal Matrix Coefficient C: Common Factor (t)     [s/m2]
      allocate  ( Kz__AT(kcolp,mzp ) )                       ! Vertical  Turbulent Diffusion Coefficient               [m2/s]
      allocate  ( Kzm_AT(kcolp,mzp ) )                       ! Vertical  Turbulent Diffusion Coefficient (Momentum)    [m2/s]
      allocate  ( Kzh_AT(kcolp,mzp ) )                       ! Vertical  Turbulent Diffusion Coefficient (Scalars)     [m2/s]
      allocate  ( Kzh0AT(kcolp,mzp ) )                       ! Vertical  Turbulent Diffusion Coefficient (Scalars)     [m2/s]
      allocate  ( A___AT(kcolp,mzp ) )                       ! Tridiagonal Matrix Coefficient A                           [-]
      allocate  ( B___AT(kcolp,mzp ) )                       ! Tridiagonal Matrix Coefficient B                           [-]
      allocate  ( C___AT(kcolp,mzp ) )                       ! Tridiagonal Matrix Coefficient C                           [-]
      allocate  ( D___AT(kcolp,mzp ) )                       ! Independant Term               D                           [x]
      allocate  ( P___AT(      mzp ) )                       ! Auxiliary   Term               P                           [-]
      allocate  ( Q___AT(      mzp ) )                       ! Auxiliary   Term               Q                           [-]
      allocate  ( X___AT(      mzp ) )                       ! Auxiliary   Unknown            X                           [x]

      allocate  ( LMO_AT(kcolp     ) )                       ! Monin-Obukhov     Length         (Grid Cell Average)       [m]
      allocate  ( zi__AT(kcolp     ) )                       ! Inversion         Height         (Grid Cell Average)[m a.g.l.]
      allocate  ( TKE_AT(kcolp,mzp ) )                       ! Turbulent Kinetic Energy                               [m2/s2]
      allocate  ( eps_AT(kcolp,mzp ) )                       ! Turbulent Kinetic Energy Dissipation                   [m2/s3]
      allocate  ( TrT_AT(kcolp,mzp ) )                       ! Turbulent Kinetic Energy Transport                     [m2/s3]

      allocate  ( dua_AT(kcolp,mzp ) )                       ! Wind Speed (x-direc.) Tendency                          [m/s2]
      allocate  ( dva_AT(kcolp,mzp ) )                       ! Wind Speed (y-direc.) Tendency                          [m/s2]
      allocate  ( dpktAT(kcolp,mzp ) )                       ! Potential Temperature Tendency, divided by p0**(R/Cp)   [x/s]
      allocate  ( dqv_AT(kcolp,mzp ) )                       ! Specific  Humidity    Tendency                      [kg/kg/s]
      allocate  ( dqw_AT(kcolp,mzp ) )                       ! Cloud Droplets Concen.Tendency                      [kg/kg/s]
      allocate  ( dqi_AT(kcolp,mzp ) )                       ! Cloud Crystals Concen.Tendency                      [kg/kg/s]
      allocate  ( dCi_AT(kcolp,mzp ) )                       ! CCNi           Concen.Tendency                          [1/s]
      allocate  ( dqs_AT(kcolp,mzp ) )                       ! Snow Particles Concen.Tendency                      [kg/kg/s]
      allocate  ( dqr_AT(kcolp,mzp ) )                       ! Rain Drops     Concen.Tendency                      [kg/kg/s]

! =================================
! ALLOCATION Mod_PHY_AT_kkl -   END
! =================================


      end subroutine PHY_Atm_AT_ALLOC
