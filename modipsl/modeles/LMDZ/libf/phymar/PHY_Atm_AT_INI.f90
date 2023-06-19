      subroutine PHY_Atm_AT_INI(FlagAT_TKE,TypeAT)

!------------------------------------------------------------------------------+
!                                                         Sat 29-Jun-2013  MAR |
!   MAR          PHY_Atm_AT_INI                                                |
!     subroutine PHY_Atm_AT_INI intializes Turbulent Vertical Diffusion Scheme |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,               Sat 29-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_AT_ctr
      use Mod_PHY_AT_grd
      use Mod_PHY_AT_kkl



      IMPLICIT NONE




! Arguments
! =========

      logical                                    ::  FlagAT_TKE         !  Flag         (Turbulent Transfer, TKE-e Model: ON / OFF)
      character(len=1)                           ::  TypeAT             !  Type         (Turbulent Transfer        Model: Choice  )




! Local Variables
! ===============

      integer    ::    k
!     integer    ::    i     ,j     ,ikl




! Initialization of Mod_PHY_AT_ctr  (Atm_AT Switches
! ================================        & Time/Space control Variables)

! Atm_AT Switch
! -------------

      AT_TKE = FlagAT_TKE



! E-e and K-l models parameters
! -----------------------------

        Ee_Duynkerke = .FALSE.                                          ! Dunkerke           (1988) E-epsilon model of turbulence
        Ee_Kitada    = .FALSE.                                          ! Kitada             (1987) E-epsilon model of turbulence
        Ee_HuangRamn = .FALSE.                                          ! Huang and Raman    (1991) E-epsilon model of turbulence
        Kl_TherryLac = .FALSE.                                          ! Therry & Lacarrere (1983) K-l       model of turbulence

      IF      (TypeAT.EQ.'e')                                       THEN!
        Ee_Duynkerke =  .TRUE.                                          ! Dunkerke           (1988) E-epsilon model of turbulence
      ELSE IF (TypeAT.EQ.'K')                                       THEN!
        Ee_Kitada    =  .TRUE.                                          ! Kitada             (1987) E-epsilon model of turbulence
      ELSE IF (TypeAT.EQ.'H')                                       THEN!
        Ee_HuangRamn =  .TRUE.                                          ! Huang and Raman    (1991) E-epsilon model of turbulence
      ELSE IF (TypeAT.EQ.'L')                                       THEN!
        Kl_TherryLac =  .TRUE.                                          ! Therry & Lacarrere (1983) K-l       model of turbulence
      ELSE
            write(6,*)
            write(6,*) ' TypeAT = ',TypeAT,' is inadequate'
            write(6,*) ' Please Choose e (Duynkerke)  ,    K (Kitada)                           '
            write(6,*) '               H (Huang Raman), OR L (Therry-Lacarrere) Turbulence Model'
            write(6,*) '                                                    '

!                **************
            call MAR________BUG
!                **************

      END IF




! Initialisation of the TKE Scheme
! ================================

!                **************
          CALL   PHY_genTKE_INI
!                **************




! Initialisation of the Tri-Diagonal Matrix
! =========================================

        alphAT    =   0.25
        betaAT    =   1.00 - alphAT
        a_b_AT    =          alphAT / betaAT

        Ac0_AT(mzp) = - dt__AT *Grav_F *Grav_F *betaAT    /(dsigmi(mzp) * dsigma(mzp))
      DO k=    mzp-1,1,-1
        Ac0_AT(k)   = - dt__AT *Grav_F *Grav_F *betaAT    /(dsigmi(k)   * dsigma(k))
        Cc0_AT(k+1) =   Ac0_AT(k)                          *dsigmi(k)   / dsigmi(k+1)
      ENDDO



      end subroutine PHY_Atm_AT_INI
