      subroutine PHY_genTKE_RUN

!------------------------------------------------------------------------------+
!                                                         Sat 22-Jun-2013  MAR |
!   MAR          PHY_genTKE_RUN                                                |
!     subroutine PHY_genTKE_RUN computes   Turbulent Kinetic Energy            |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Fri 15-Mar-2013      |
!           Last Modification by H. Gallee,               Sat 22-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+
!                                                                              |
!   METHOD: 1. `Standard' Closure of Turbulence:                               |
!   ^^^^^^^     E - epsilon  , Duynkerke,           JAS 45, 865--880, 1988     |
!           2.  E - epsilon  , Huang and Raman,     BLM 55, 381--407, 1991     |
!           3.  K - l        , Therry et Lacarrere, BLM 25,  63-- 88, 1983     |
!                                                                              |
!   INPUT  : itexpe: Nb of iterations                                          |
!   ^^^^^^^^ dt__AT: Local Time Step                                      (s)  |
!            explIO: Experiment Label                                     (s)  |
!                                                                              |
!   INPUT / OUTPUT: The Vertical Turbulent Fluxes are included for:            |
!   ^^^^^^^^^^^^^^                                                             |
!     a) Turbulent Kinetic Energy             TKE_AT(kcolq,mzp)       (m2/s2)  |
!     b) Turbulent Kinetic Energy Dissipation eps_AT(kcolq,mzp)       (m2/s3)  |
!                                                                              |
!   OUTPUT :  Kzm_AT(kcolq,mzp): vertic. diffusion coeffic. (momentum) (m2/s)  |
!   ^^^^^^^^  Kzh_AT(kcolq,mzp): vertic. diffusion coeffic. (scalars ) (m2/s)  |
!             zi__AT(kcolq)    : inversion height                      (m)     |
!                                                                              |
!   OPTIONS: #HY: Latent Heat Exchanges associated with Cloud Microphysics     |
!   ^^^^^^^^                            contribute   to Buoyancy               |
!            #KA: replaces TKE & e below a prescribed height above the surface |
!                          by their      box weighted vertical moving averages |
!            #KC: Modification of TKE near the Lower Boundary                  |
!            #LD: Buoyancy includes Loading by the Hydrometeors                |
!            #RI: Correction of the Prandtl    Nb (Kzm/Kzh)                    |
!                        by the The Richardson Nb                              |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_Phy____dat
      use Mod_Phy____grd
      use Mod_PHY_AT_ctr
      use Mod_PHY_AT_grd
      use Mod_PHY____kkl
      use Mod_PHY_AT_kkl
      use Mod_PHY_DY_kkl
      use Mod_PHY_CM_kkl
      use Mod_SISVAT_gpt



!  Local  Variables
!  ================

      use Mod_genTKE_RUN


      IMPLICIT NONE


      real(kind=real8)   ::   z__SBL                  ! SBL Thickness assumed to be the lowest Model Layer Height         [m]
      real(kind=real8)   ::   TKE_zi                  ! TKE at the inversion height                                   [m2/s2]
      integer            ::   kTKEzi                  ! level correponding to TKE_zi
      real(kind=real8)   ::   dTKEdk                  ! TKE difference across model Layers                            [m2/s2]
      real(kind=real8)   ::   kz_inv                  ! 1 / (k z)
      real(kind=real8)   ::   Buoy_F                  ! Contribution of  Buoyancy Force
      real(kind=real8)   ::   dTKE_p                  ! Positive Contribution  to TKE
      real(kind=real8)   ::   r_eTKE                  ! e /  TKE Ratio  (used  in the Product./Destruct. Scheme of TKE)
      real(kind=real8)   ::   TKE_ds                  ! Verticaly Integrated      TKE
      real(kind=real8)   ::   eps_ds                  ! Verticaly Integrated      TKE Dissipation
      real(kind=real8)   ::   edt_HR                  ! Minimum Energy Dissipation Time              
      real(kind=real8)   ::   TKEnew                  ! New Value              of TKE
      real(kind=real8)   ::   TKEsbc                  ! SBC:                      TKE                                 [m2/s2]
      real(kind=real8)   ::   epssbc                  ! SBC:                      TKE Dissipation                     [m2/s3]
      real(kind=real8)   ::   zeta                    ! SBL:      z / LMO                                                 [-]
      real(kind=real8)   ::   u_star                  ! SBL:      u*             (Friction Velocity)                    [m/s]
      real(kind=real8)   ::   kz_max,kz_mix,kz2mix    ! Kzh       Auxiliary       Variables
      real(kind=real8)   ::   KzhMAX                  ! max.vert. Turbulent Diffusion Coefficient (Scalars)            [m2/s]

! #KC real(kind=real8)   ::   se                      ! 1    if TKE(mzp-2) > TKE(mzp-1) and 0 otherwise
! #KC integer            ::   ke                      ! index of the level where TKE is max 

      real(kind=real8)   ::   relHum                  ! Relative       Humidity                                           [-]
      real(kind=real8)   ::   QS_mid                  ! Saturation Sp. Humidity % Liquid Water                        [kg/kg]
      real(kind=real8)   ::   qwsLRT                  ! Qws * L / (Rd * T)
      real(kind=real8)   ::   surSat                  ! Normalized sur-Saturation (> 0 if relHum > 1)
      real(kind=real8)   ::   C_thq                   ! (Duynkerke & Driedonks 1987 JAS 44(1), Table 1 p.47)
      real(kind=real8)   ::   C_q_w                   ! (Duynkerke & Driedonks 1987)
      real(kind=real8)   ::   CX__hi                  !
      real(kind=real8)   ::   Ce1_hi                  ! ... Ce1 / hi
      real(kind=real8)   ::   Ck1_hi                  ! ... Ck1 / hi
      real(kind=real8)   ::   Ck2_hi                  ! ... Ck2 / hi
      real(kind=real8)   ::   Th_Lac                  !  Therry & Lacarrere (1983)      Parameter

      real(kind=real8)   ::   sgnLMO                  ! sign(LMO)
      real(kind=real8)   ::   absLMO                  !     |LMO|
! #RI real(kind=real8)   ::   fac_Ri,vuzvun,Kz_vun    ! Sukorianski    Parameterization Parameters

      integer            ::   ikl
      integer            ::   i
      integer            ::   j
      integer            ::   k




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    !
! ALLOCATION
! ==========

      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                          THEN !

          allocate  ( dukkp1(mzpp) )                                 !      Difference (u(k) - u(k+1))                      [m/s]
          allocate  ( dvkkp1(mzpp) )                                 !      Difference (v(k) - v(k+1))                      [m/s]
          allocate  ( kkp1dz(mzpp) )                                 !  1 / Difference (Z(k) - Z(k+1))                      [1/m]
          allocate  ( zShear(mzp)  )                                 !      Wind Shear Contribution to TKE                [m2/s3]
          allocate  ( REq_PT(mzpp) )                                 !  Reduced (Equivalent) Potential Temperature            [K]
          allocate  ( c_Buoy(mzpp) )                                 !  Buoyancy Coefficient (g/theta) X (dtheta/dz)       [1/s2]
          allocate  ( Ri__Nb(mzp)  )                                 !  Richardson Number                                     [-]
          allocate  ( Prandtl(mzp) )                                 !  Prandtl    Number (Kzm/Kzh)                           [-]
          allocate  ( Ls_inv(mzp)  )                                 !  1 / Ls  (Therry & Lacarrere, 1983)                  [1/m]
          allocate  ( ML_inv(mzp)  )                                 !  1 / ML  (Mixing      Length, Therry & Lacarr, 1983) [1/m]
          allocate  ( DL_inv(mzp)  )                                 !  1 / DL  (Dissipation Length, Therry & Lacarr, 1983) [1/m]
          allocate  ( Dissip(mzp)  )                                 !           Dissipation                              [m2/s3]
          allocate  ( TKEvav(mzp)  )                                 !           TKE         Vertical moving Average      [m2/s2]
          allocate  ( epsvav(mzp)  )                                 !           Dissipation Vertical moving Average      [m2/s3]
          allocate  ( pkt   (mzpp) )                                 !           Reduced     Potential       Temperature      [X]

      END IF
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!     ++++++++++++++
      DO ikl=1,kcolp
!     ++++++++++++++

!        i = ii__AP(ikl)
!        j = jj__AP(ikl)




! Friction Velocity
! =================

         u_star = us__SV_gpt(ikl)




! Verification: TKE must be Positive Definite
! ===========================================

         DO k=1,mzp
           TKE_AT(ikl,k)=max(eps6,TKE_AT(ikl,k))
           eps_AT(ikl,k)=max(epsn,eps_AT(ikl,k))
         END DO




! Reduced Potential Temperature
! =============================

         DO k=1,mzpp
           pkt   (k)    =         pkt_DY(ikl,k)
         END DO




! Inversion Height
! ================

           TKE_zi     = 0.05*max(max(TKE_AT(ikl,mzp-1),                  &
     &                               TKE_AT(ikl,mzp  )),TKEmin)
           kTKEzi     = 1 


         DO k=1,mzp
           IF   (TKE_AT(ikl,k) .lt. TKE_zi    )                     THEN
                 kTKEzi      =  min(mzp-1,k)
           END IF
         ENDDO

             k = kTKEzi
         IF      (TKE_AT(ikl,k+1).lt.TKEmin)                        THEN
                  TKE_zi     =      ZmidDY(ikl,mzp)-sh__AP(ikl)
         ELSE
                  dTKEdk     =      TKE_AT(ikl,k)  -TKE_AT(ikl,k+1)
                  TKE_zi     =      ZmidDY(ikl,k+2)                     &
     &                            +(ZmidDY(ikl,k+1)-ZmidDY(ikl,k+2))    &
     &                            *(TKE_zi         -TKE_AT(ikl,k+1))    &
     &             /(sign(un_1 ,    dTKEdk)                             &
     &               *max(epsn ,abs(dTKEdk)))      -sh__AP(ikl)
         END IF

                 zi__AT(ikl) =  min(TKE_zi, ZmidDY(ikl,  1)-sh__AP(ikl))

                 zi__AT(ikl) =  max(        ZmidDY(ikl,mzp)-sh__AP(ikl) &
     &                             ,zi__AT(ikl))
                 TKE_zi      =  0.
                 kTKEzi      =  0 




! TKE Production/Destruction by the Vertical Wind Shear
! =====================================================

         DO k=1,mzp-1
           dukkp1(k)  =  ua__DY(ikl,k) - ua__DY(ikl,k+1)
           dvkkp1(k)  =  va__DY(ikl,k) - va__DY(ikl,k+1)
         END DO
            k=  mzp
           dukkp1(k)  =  ua__DY(ikl,k)
           dvkkp1(k)  =  va__DY(ikl,k)

         DO k=1,mzp
           kkp1dz(k)  =  Z___DY(ikl,k) - Z___DY(ikl,k+1)                !  dz(k+1/2)
         END DO

         DO k=1,mzp-1
           zShear(k)    =                                              &
     &     Kzm_AT(ikl,k)*(dukkp1(k) *dukkp1(k) + dvkkp1(k) *dvkkp1(k)) &
     &                  /(kkp1dz(k) *kkp1dz(k))
         END DO
           zShear(mzp)  = 0.0




! Buoyancy
! ========

! Reduced (Equivalent) Potential Temperature
! ------------------------------------------

! Control Martin
!PRINT*,'CpdAir=',CpdAir
!PRINT*,'minval(Ta__DY(ikl,:))=',minval(Ta__DY(ikl,:))
! Control Martin

         DO k=1,mzpp
           REq_PT(k)     =     pkt   (    k)                            &
! #LD&        * (1.0 + 0.715 * ld_H2O(ikl,k)   )                        &
     &        *   exp(LhvH2O * qv__DY(ikl,k)                            &
     &             / (CpdAir * Ta__DY(ikl,k))  )                        &
     &        +  0.0
         END DO



! Buoyancy Coefficients
! ---------------------

         DO k=1,mzp

           relHum     = 0.50 *(qv__DY(ikl,k  ) /qvswCM(ikl,k  )         &
     &                        +qv__DY(ikl,k+1) /qvswCM(ikl,k+1))
           QS_mid     = 0.50 *(qvswCM(ikl,k  ) +qvswCM(ikl,k+1))

           surSat     = max(zer0,sign(un_1,relHum +epsp-un_1))
           qwsLRT     = LhvH2O*QS_mid /(R_DAir *TmidDY(ikl,k))

! Vectorization of the unsaturated (H<1) and saturated cases (H=1.)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           C_thq =1.-surSat                                             & ! H<1.
     &              +surSat*(1.+qwsLRT)                                 & ! H=1.
     &                     /(1.+qwsLRT*.605*LhvH2O                      & !
     &                                    /(CpdAir*TmidDY(ikl,k)))        !
! ...      C_thq (Duynkerke & Driedonks 1987 JAS 44(1), Table 1 p.47)     !

           C_q_w=(1.-surSat)              *(LhvH2O                      & !
     &                                    /(CpdAir*TmidDY(ikl,k))       & ! H<1.
     &                                     - 0.605                    ) & !
     &              +surSat                                               ! H=1.
! ...      C_q_w (Duynkerke and Driedonks 1987)                           !
!                 with (1-Ra/Rv)/(Ra/Rv) =  0.605 [Ra=287.J/kg/K;         !
!                                                  Rv=461.J/kg/K]         !

! Unsaturated Case
! ~~~~~~~~~~~~~~~~
!         IF(relHum.lt.1.0)                                         THEN  !
!            C_thq = 1.0
!            C_q_w =                   LhvH2O/(CpdAir*TmidDY(ikl,k))    & !
!    &                                        - 0.605                     !

! Saturated   Case
! ~~~~~~~~~~~~~~~~
!         ELSE                                                            !
!            qwsLRT=      QS_mid      *LhvH2O/(RDryAi*TmidDY(ikl,k))      !
!            C_thq = (1.0+qwsLRT)                                       & !
!    &              /(1.0+qwsLRT*0.605*LhvH2O/(CpdAir*TmidDY(ikl,k)))     !
!            C_q_w =  1.0
!         END IF


! Buoyancy
! --------

          IF(k.EQ.mzp)      C_q_w = 0.0

             c_Buoy(k)    = Grav_F                                      &
     &      *kkp1dz(k) * ( (REq_PT(k)-REq_PT(k+1))                      &
     &                 *2./(REq_PT(k)+REq_PT(k+1))                      &
     &                 *    C_thq                                       &
     &                 -    C_q_w *(qv__DY(ikl,k)-qv__DY(ikl,k+1)       &
     &                             +qw__CM(ikl,k)-qw__CM(ikl,k+1)       &
     &                             +qr__CM(ikl,k)-qr__CM(ikl,k+1)       &
     &                             +qi__CM(ikl,k)-qi__CM(ikl,k+1)       &
     &                             +qs__CM(ikl,k)-qs__CM(ikl,k+1))      &
     &                                                             )
! ...       (g/theta)                X(dtheta/dz) : 
!            Buoyancy Parameter beta X grad.vert.temp.pot. en k+1/2

         END DO


! Dissipation & Length Scales Parameters (Therry and Lacarrere 1983 Model)
! ======================================

         IF (Kl_TherryLac)                                          THEN
           CX__hi      =  1.0 / zi__AT(ikl)
           Ce1_hi      = 15.0 * CX__hi       ! ...    Ce1/hi 
           Ck1_hi      =  5.0 * CX__hi       ! ...    Ck1/hi 
           Ck2_hi      = 11.0 * CX__hi       ! ...    Ck2/hi 

           sgnLMO      = sign(un_1,LMO_AT(ikl))
           absLMO      =  abs(     LMO_AT(ikl))
           LMO_AT(ikl) =           sgnLMO    * max(absLMO,epsp)
           Th_Lac      = -min(0.00,sgnLMO)                              &
     &             /(1.0 -min(0.00,LMO_AT(ikl))  / zi__AT(ikl))

! Replacement of:
!          IF (LMO_AT(ikl).lt.0.0)                                  THEN
!              LMO_AT(ikl) =        min(LMO_AT(ikl),-epsp)
!              Th_Lac =  1.0/(1.d0-LMO_AT(ikl)/zi__AT(ikl)) 
!          ELSE
!              LMO_AT(ikl) =        max(LMO_AT(ikl), epsp)
!              Th_Lac =  0.0
!          END IF
! ...      m2


          DO k=1,mzp
           Ls_inv(k)   = sqrt( max(zer0,c_Buoy(k)) /TKE_AT(ikl,k) )      ! ...  1/ls
          END DO


! Dissipation Length
! ------------------

          DO k=1,mzp
           kz_inv=vK_inv/(ZmidDY(ikl,k+1)-sh__AP(ikl))                   ! ...  1/kz(i,j,k+1/2)
!   
           DL_inv(k)=kz_inv +Ce1_hi                                     &! ...  DL_inv=1/Dissipation Length
     &      -(kz_inv+Ck1_hi)*Th_Lac/(1.+5.0e-3*zi__AT(ikl)*kz_inv)      &!      (Therry and Lacarrere, 1983 BLM 25 p.75)
     &       +1.5*Ls_inv(k)


! Mixing Length
! -------------

           ML_inv(k)=kz_inv +Ce1_hi                                     &! ...  ML_inv=1/Mixing Length 
     &      -(kz_inv+Ck2_hi)*Th_Lac/(1.+2.5e-3*zi__AT(ikl)*kz_inv)      &!      (Therry and Lacarrere, 1983 BLM 25 p.78)
     &       +3.0*Ls_inv(k)

           Dissip(k) = 0.125*DL_inv(k)*sqrt(TKE_AT(ikl,k))*TKE_AT(ikl,k)
          ENDDO




! Dissipation                            (E-epsilon                 Models)
! ===========                                    

         ELSE

          DO k=1,mzp
           Dissip(k) =       eps_AT(ikl,k)
          ENDDO

         END IF




! Contribution of Vertical Wind Shear + Buoyancy + Dissipation to TKE 
! ===================================================================

         DO k=1,mzp
           Buoy_F=-Kzh_AT(ikl,k) *       c_Buoy(k)    

           TKEnew= TKE_AT(ikl,k) *                                           &
     &            (TKE_AT(ikl,k)+dt__AT*(zShear(k)    +max(zer0,Buoy_F)))    &
     &           /(TKE_AT(ikl,k)+dt__AT*(             -min(zer0,Buoy_F)      &
     &                                                +         Dissip(k) ))
! ...      Numerical Scheme : cfr. E. Deleersnijder, 1992 (thesis) pp.59-61 




! Contribution of Vertical Wind Shear + Buoyancy to epsilon
! =========================================================

          IF                     (Ee_Duynkerke)                     THEN
           dTKE_p = zShear(k)    +max(zer0,Buoy_F)+max(TrT_AT(ikl,k),zer0)
          ELSE
           dTKE_p = zShear(k)    +max(zer0,Buoy_F)
! ...      based on standard values of Kitada, 1987, BLM 41, p.220
          END IF

! #BH      dTKE_p = zShear(k)    +max(zer0,Buoy_F) * 1.80               &
! #BH&                           -min(Buoy_F,zer0) * 1.15
! ...      based on          values of Betts et Haroutunian, 1983
!          can be used by replacing strings `c #KI' (except the previous one) 
!                                       and `c #BH' by blanks 
!                                (cfr. Kitada, 1987, BLM 41, p.220):
!          buoyancy > 0 (unstability) => (1-ce3) X buoyancy = 1.8  X buoyancy
!          buoyancy < 0 (  stability) => (1-ce3) X buoyancy =-1.15 X buoyancy

           r_eTKE = eps_AT(ikl,k) /TKE_AT(ikl,k)
           eps_AT(ikl,k) =                                              &
     &              eps_AT(ikl,k)                                       &
     &            *(eps_AT(ikl,k) +dt__AT *r_eTKE *c1ep *dTKE_p)        &
     &            /(eps_AT(ikl,k) +dt__AT *r_eTKE *c2ep *eps_AT(ikl,k)) 
! ...      Numerical Scheme : cfr. E. Deleersnijder, 1992 (thesis) pp.59-61 

          IF                     (Kl_TherryLac)                         &
     &     eps_AT(ikl,k)=Dissip(k) 




! New TKE Value
! =============

           TKE_AT(ikl,k)=TKEnew

         END DO



! Lower Boundary Conditions
! =========================

           TKEsbc          =  u_star          * u_star                    ! -> TKE SBC
           z__SBL          =  Z___DY(ikl,mzp) - sh__AP(ikl)               !    z_SBL

           epssbc          =  TKEsbc          * u_star                    ! -> e   SBC
           TKE_AT(ikl,mzp) =  TKEsbc          * sqrcmu                    !    TKE SBC

           zeta            =  z__SBL          / LMO_AT(ikl)               !    zeta
           sgnLMO          =  max(0.0,sign(un_1,LMO_AT(ikl)))             !

           eps_AT(ikl,mzp) = epssbc                                     & !
     &              *(      (sgnLMO     *(1.+A_Stab*    zeta)           & ! phim Stab.
     &                 +(1.0-sgnLMO    )/(1.-20.*min(0.,zeta)))         & ! phim Inst.
     &              *vK_inv /z__SBL                                     & !
     &              -vK_inv /LMO_AT(ikl))
! ...      Duynkerke, 1988, JAS (45), (19) p. 869

! #KI      eps_AT(ikl,mzp) = epssbc                                     &
! #KI&              *vK_inv /z__SBL     


! When TKE Closure is Used, TKE is Modified near the Lower Boundary
! -----------------------------------------------------------------

! #KC      se = max(0.,sign(un_1,TKE_AT(ikl,mzp-2)-TKE_AT(ikl,mzp-1)))
! #KC      ke =                                           mzp-1-se
! #KC      TKE_AT(ikl,mzp-1)=     TKE_AT(ikl,ke  )
! #KC      eps_AT(ikl,mzp-1)=     eps_AT(ikl,ke  )
! ...      Schayes and Thunis, 1990, Contrib. 60 Inst.Astr.Geoph. p.8, 1.4.4.

! #KC      TKE_AT(ikl,mzp-1)=     TKE_AT(ikl,mzp )
! #KC      eps_AT(ikl,mzp-1)=     eps_AT(ikl,mzp )


! Upper Boundary Conditions
! =========================

           TKE_AT(ikl,1) = TKE_AT(ikl,2)
           eps_AT(ikl,1) = eps_AT(ikl,2)
        


! TKE-e Vertical Moving Average
! =============================

           DO k=     1,mzp
             TKEvav(k)=                TKE_AT(ikl,    k )
             epsvav(k)=                eps_AT(ikl,    k )
           END DO
         IF (AT_vav)                                                THEN
           DO k=     2,mzp-1
             TKEvav(k)=(dsigmi(k1p(k))*TKE_AT(ikl,k1p(k))               &
     &                 +dsigmi(    k )*TKE_AT(ikl,    k )*2.            &
     &                 +dsigmi(k1m(k))*TKE_AT(ikl,k1m(k))   )           &
     &                /(dsigmi(k1p(k))                                  &
     &                 +dsigmi(    k )                   *2.            &
     &                 +dsigmi(k1m(k))                      )
             epsvav(k)=(dsigmi(k1p(k))*eps_AT(ikl,k1p(k))               &
     &                 +dsigmi(    k )*eps_AT(ikl,    k )*2.            &
     &                 +dsigmi(k1m(k))*eps_AT(ikl,k1m(k))   )           &
     &                /(dsigmi(k1p(k))                                  &
     &                 +dsigmi(    k )                   *2.            &
     &                 +dsigmi(k1m(k))                      )
           END DO
         END IF

! #KA    DO k=mz__KA,mzp-1
! #KA      TKE_AT(ikl,k)= TKEvav(k)
! #KA      eps_AT(ikl,k)= epsvav(k)
! #KA    END DO


! Verification: TKE must be Positive Definite
! ===========================================

         DO k=1,mzp
           TKE_AT(ikl,k)=max(eps6,TKE_AT(ikl,k))
           eps_AT(ikl,k)=max(epsn,eps_AT(ikl,k))
           TKEvav(k)    =max(eps6,TKEvav(k)    )
           epsvav(k)    =max(eps6,epsvav(k)    )
         END DO


! Minimum Energy Dissipation Time
! ===============================

         IF                      (Ee_HuangRamn)                     THEN
           TKE_ds      = 0.0
           eps_ds      = 0.0
          DO k=1,mzp
           TKE_ds      = TKE_ds + TKE_AT(ikl,k)*dsigma(k)
           eps_ds      = eps_ds + eps_AT(ikl,k)*dsigma(k)
          END DO

          IF (eps_ds.gt.0.0)                                        THEN
           edt_HR      = betahr * TKE_ds       /eps_ds 
          ELSE
           edt_HR      = epsn 
! ...      edt_HR set to an arbitrary small value
          END IF
         END IF


! Turbulent Diffusion Coefficients
! ================================

         IF (it_EXP.gt.1)                                           THEN


! Richardson Number (contributors)
! -----------------

           Prandtl(mzp)    =   1.

          DO k=2,mzp-1
           c_Buoy(k) =  0.0
! #RI      c_Buoy(k) =      (pkt   (    k)-pkt   (    k+1))*p0_kap      &
! #RI&                   /   TmidDY(ikl,k)
! #RI      zShear(k) =   max(dukkp1(k)**2 +dvkkp1(k) **2 ,  eps6)


! Richardson Number
! -----------------

           Ri__Nb(k) =  0.0
! #RI      Ri__Nb(k) = (Grav_F/kkp1dz(k))                               & ! g * dz     (k+1/2)
! #RI&                        *c_Buoy(k)                                & ! d(theta)/T (k+1/2)
! #RI&                        /zShear(k)                                  ! d|V|
          END DO


! Diffusion Coefficient for Heat
! ------------------------------

          DO k=2,mzp

          IF                 (Kl_TherryLac)                         THEN
           Kzh_AT(ikl,k)=  0.50 * sqrt(TKEvav(k    ))/ML_inv(k)     
! ...      nu_t =c_mu X ECT          X ECT          / eps

          ELSE IF            (Ee_HuangRamn)                         THEN
           Kzh_AT(ikl,k)=                                               &
     &            cmu * TKE_AT(ikl,k)* TKE_AT(ikl,k)/(eps_AT(ikl,k)     &
     &                                +TKE_AT(ikl,k)/ edt_HR       )
! ...      nu_t =c_mu X ECT          X ECT          / eps

          ELSE          !    (Ee_Duynkerke / Ee_Kitada)
           Kzh_AT(ikl,k)=                                               &
     &            cmu * TKEvav(k    )* TKEvav(k    )/(epsvav(k    ))
! ...      nu_t =c_mu X ECT          X ECT          / eps

          END IF

           kz_max= vonKrm *     Z___DY(ikl,k+1)-Z___DY(ikl,mzpp)
           kz_mix= kz_max /    (1.             +kz_max          *0.1)
           kz2mix= kz_mix *     kz_mix
           KzhMAX= max(   5000.        ,   100.                         &
     &   * kz2mix *abs((WindDY    (ikl,k)-WindDY    (ikl,k1p(k)))       &
     &                         *kkp1dz(k)                      ))
           Kzh_AT(ikl,k)=   min( KzhMAX , Kzh_AT(ikl,k))
           Kzh_AT(ikl,k)=   max( A_MolV , Kzh_AT(ikl,k))
           Kzh0AT(ikl,k)=                 Kzh_AT(ikl,k)
          END DO



! Flux Limitor (Limitation of pkt Flux)
! ------------

      IF  (AT_LIM .AND. mzp.GE.3)                                   THEN
        DO k=min(3,mzp),mzp
          IF                 (pkt_DY(ikl,k+1) .GT. pkt_DY(ikl,k  )) THEN
           Kz_MAX=(          (Z___DY(ikl,k  ) -    Z___DY(ikl,k+1))    &
     &                      /(pkt_DY(ikl,k+1) -    pkt_DY(ikl,k  )))   &
     &     *(Dpkt_X    *0.5 *(Z___DY(ikl,k-1) -    Z___DY(ikl,k+1))    &
     &      -Kzh_AT(ikl,k-1)*(pkt_DY(ikl,k-1) -    pkt_DY(ikl,k  ))    &
     &                      /(Z___DY(ikl,k-1) -    Z___DY(ikl,k  )))
          ELSE IF            (pkt_DY(ikl,k+1) .LT. pkt_DY(ikl,k  )) THEN
           Kz_MAX=(          (Z___DY(ikl,k  ) -    Z___DY(ikl,k+1))    &
     &                      /(pkt_DY(ikl,k  ) -    pkt_DY(ikl,k+1)))   &
     &     *(Dpkt_X    *0.5 *(Z___DY(ikl,k-1) -    Z___DY(ikl,k+1))    &
     &      +Kzh_AT(ikl,k-1)*(pkt_DY(ikl,k-1) -    pkt_DY(ikl,k  ))    &
     &                      /(Z___DY(ikl,k-1) -    Z___DY(ikl,k  )))
          END IF
           Kzh_AT(ikl,k) = min(Kzh_AT(ikl,k),Kz_MAX)
           Kzh_AT(ikl,k) = max(Kzh_AT(ikl,k),A_MolV)
        ENDDO
      ENDIF



! Prandtl Number (Sukoriansky et al., 2005,
! --------------  BLM 117: 231-257, Eq.15, 19, 20 & Fig.2)

          DO k=2,mzp-1
! #RI      fac_Ri= 5.0 *     max(Ri__Nb(i,j,k), eps6)
! #RI      vuzvun= 0.4 *(1.-(fac_Ri-1./fac_Ri)/(fac_Ri+1./fac_Ri)) + 0.2
! #RI      fac_Ri= 4.2 *     max(Ri__Nb(i,j,k), eps6)
! #RI      Kz_vun= 0.7 *(1.-(fac_Ri-1./fac_Ri)/(fac_Ri+1./fac_Ri))
           Prandtl(k)    =                             1.
! #RI      Prandtl(k)    =   max(0.     ,sign(1.  , Kzh_AT(ikl,k)-0.20)) &
! #RI&                   -   min(0.     ,sign(1.  , Kzh_AT(ikl,k)-0.20)) &
! #RI&                   *   min(vuzvun / max(eps6,Kz_vun),     20.00)
          END DO


! Diffusion Coefficient for Momentum
! ----------------------------------

          DO k=2,mzp
          IF                 (Kl_TherryLac)                         THEN
           Kzm_AT(ikl,k)=  0.7           * Kzh_AT(ikl,k)

          ELSE
           Kzm_AT(ikl,k)=                  Kzh_AT(ikl,k)
! ...      cfr Machiels, 1992, TFE (FSA/UCL) (3.21) p.21

! #RI      Kzm_AT(ikl,k)=  Prandtl(k)    * Kzh_AT(ikl,k)
          END IF

          END DO

         END IF




! Work Arrays Reset
! =================

      DO k=1,mzp
        TrT_AT(ikl,k) = 0.0
      END DO




!     +++++++++++++++++++
      ENDDO ! ikl=1,kcolp
!     +++++++++++++++++++




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    !
! DE-ALLOCATION
! =============

      IF (FlagDALLOC)                                           THEN !
          deallocate   ( dukkp1 )                                    !      Difference (u(k) - u(k+1))                      [m/s]
          deallocate   ( dvkkp1 )                                    !      Difference (v(k) - v(k+1))                      [m/s]
          deallocate   ( kkp1dz )                                    !  1 / Difference (Z(k) - Z(k+1))
          deallocate   ( zShear )                                    !      Wind Shear Contribution to TKE                [m2/s3]
          deallocate   ( REq_PT )                                    !  Reduced (Equivalent) Potential Temperature            [K]
          deallocate   ( c_Buoy )                                    !  Buoyancy Coefficient (g/theta) X (dtheta/dz)       [1/s2]
          deallocate   ( Ri__Nb )                                    !  Richardson Number                                     [-]
          deallocate   ( Prandtl)                                    !  Prandtl    Number (Kzm/Kzh)                           [-]
          deallocate   ( Ls_inv )                                    !  1 / Ls  (Therry & Lacarrere, 1983)                  [1/m]
          deallocate   ( ML_inv )                                    !  1 / ML  (Mixing      Length, Therry & Lacarr, 1983) [1/m]
          deallocate   ( DL_inv )                                    !  1 / DL  (Dissipation Length, Therry & Lacarr, 1983) [1/m]
          deallocate   ( Dissip )                                    !           Dissipation                              [m2/s3]
          deallocate   ( TKEvav )                                    !           TKE         Vertical moving Average      [m2/s2]
          deallocate   ( epsvav )                                    !           Dissipation Vertical moving Average      [m2/s3]
          deallocate   ( pkt    )                                    !           Reduced     Potential       Temperature      [X]
      END IF                                                         !
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      return

      end subroutine PHY_genTKE_RUN
