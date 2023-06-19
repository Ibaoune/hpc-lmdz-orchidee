      subroutine PHY_genTKE_INI

!------------------------------------------------------------------------------+
!                                                         Sun 16-Jun-2013  MAR |
!   MAR          PHY_genTKE_INI                                                |
!     subroutine PHY_genTKE_INI initialises TKE Equation                       |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Sat 16-Mar-2013      |
!           Last Modification by H. Gallee,               Sun 16-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_Phy____dat
      use Mod_Phy____grd
      use Mod_PHY____kkl
      use Mod_PHY_AT_ctr
      use Mod_PHY_AT_kkl
      use Mod_PHY_DY_kkl


      IMPLICIT NONE



! Global Variables
! ================

      integer            ::   ikl   ,k



! Local  Variables
! ================

      real(kind=real8)   ::   dz_sig                  ! sigma Layer      Thickness                                        [m]



! Parameters
! ==========

            c1ep =   c1epd

      IF                     (Ee_Duynkerke)                         THEN
             cmu =    cmud
          sqrcmu = sqrcmud
            c1ep =   c1epd
            c2ep =   c2epd
            sige =   siged
            sigk =   sigkd

      ELSE IF                (Ee_Kitada)                            THEN
             cmu =    cmuk
          sqrcmu = sqrcmuk
            c1ep =   c1epk
            c2ep =   c2epk
            sige =   sigek
            sigk =   sigkk

      ELSE IF                (Kl_TherryLac)                         THEN
!            cmu =    cmut
          sqrcmu = sqrcmut
            sige =   siget
            sigk =   sigkt

      ELSE                 ! (Bintanja, Huang & Raman)
             cmu =    cmub
          sqrcmu = sqrcmub
            c1ep =   c1epb
            c2ep =   c2epb
            sige =   sigeb
            sigk =   sigkb
      END IF

          vK_inv =    1./ vonKrm




! Level below which box weighted vertical moving averages of TKE and e are performed
! ==================================================================================

             mz__KA=mzp-1
      11 CONTINUE
         IF (hsigma(mz__KA).GT.zz__KA.OR.mz__KA.LE.1)           GO TO 10
             mz__KA=mz__KA-1
                                                                GO TO 11
      10 CONTINUE
             write(6,1000)               mz__KA
      1000   format(/,' TKE: Moving Average until mz -',i2,' Level',/)




! pkt Turbulent Flux Limitor
! ==========================

      Dpkt_X = DThMAX / p0_kap / 3600.                                  !  From PT [K/h] to pkt = PT / p00**(R/Cp) [K/h]




! Minimum Vertical Turbulent Diffusion Coefficient (ARPS.4.0 Users Guide,
! ================================================  fin para 6.3.4 p.143)

         DO k=1,mzp-1
           dz_sig        =               hsigma(k)-hsigma(k+1)
           Kz0_AT(    k) =                      A_MolV
! #ARPS    Kz0_AT(    k) = min(0.15,     epsi * dz_sig * dz_sig )
         END DO

            k=  mzp
           dz_sig         =               hsigma(k)
           Kz0_AT(    k) =                      A_MolV
! #ARPS    Kz0_AT(    k) = min(0.15,     epsi * dz_sig * dz_sig )


! ++++++++++++++
! INITIALIZATION
! ++++++++++++++

      IF (it_EXP.EQ.1)                                              THEN


! Initial E,e
! -----------

         DO ikl=1,kcolp
         DO k=1,mzp
           TKE_AT(ikl,k) = eps6
           eps_AT(ikl,k) = epsn
! ...      These initial values of TKE and epsilon correspond to K ~ Kmol


! Initial Transport of TKE
! ------------------------

           TrT_AT(ikl,k) = 0.


! Initial Vertical Diffusion Coefficients
! ---------------------------------------

           Kzm_AT(ikl,k) = A_MolV
           Kzh_AT(ikl,k) = A_MolV


! Initial Inversion Height
! ------------------------

         END DO
           zi__AT(ikl)   =     ZmidDY(ikl,mzp)-sh__AP(ikl)
         END DO

      END IF




      end subroutine PHY_genTKE_INI
