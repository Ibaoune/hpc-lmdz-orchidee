      subroutine PHY_Atm_CP_RUN(mzc,kcolc)

!------------------------------------------------------------------------------+
!                                                         Mon 17-Jun-2013  MAR |
!     subroutine PHY_Atm_CP_RUN    interfaces                                  |
!                Bechtold et al. (2001) Convective Parameterization with MAR   |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Mon  8-Apr-2013      |
!           Last Modification by H. Gallee,               Mon 17-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+
!                                                                              |
!     INPUT                                                                    |
!     ^^^^^        it_EXP             : Experiment Iteration Counter           |
!                  dt__CP             : Mass Flux  Scheme:   Time Step         |
!                  dxHOST             : grid spacing  of  HOST MODEL       [m] |
!                  Ta__DY(kcolp,mzpp) : air  temperature                   [K] |
!                  qs__CM(kcolp,mzpp) : air  snow  Particl. concentr.  [kg/kg] |
!                  qr__CM(kcolp,mzp ) : air  rain  drops    concentr.  [kg/kg] |
!                                                                              |
!     INPUT/OUTPUT                                                             |
!     ^^^^^^^^^^^^ qv__DY(kcolp,mzpp) : Specific Humidity              [kg/kg] |
!                  qw__CM(kcolp,mzp ) : air  cloud droplets concentr.  [kg/kg] |
!                  qi__CM(kcolp,mzp ) : air  cloud crystals concentr.  [kg/kg] |
!                  snowCP(kcolp     ) : Snow (convective)                  [m] |
!                  snowCM(kcolp     ) : Snow (convective + stratiform)     [m] |
!                  rainCP(kcolp     ) : Rain (convective)                  [m] |
!                  rainCM(kcolp     ) : Rain (convective + stratiform)     [m] |
!                                                                              |
!     OUTPUT       dpktCP(kcolp,mzp ) : Reduc. Pot.Temperat.Tendency   [K/X/s] |
!     ^^^^^^       dqv_CP(kcolp,mzp ) : Specific   Humidity Tendency [kg/kg/s] |
!                  dqw_CP(kcolp,mzp ) : cloud dropl.Concent.Tendency [kg/kg/s] |
!                  dqi_CP(kcolp,mzp ) : cloud cryst.Concent.Tendency [kg/kg/s] |
!                                                                              |
!                  dss_CP(kcolp     ) : Snow (convective)   Tendency     [m/s] |
!                  drr_CP(kcolp     ) : Rain (convective)   Tendency     [m/s] |
!                                                                              |
!                  pkt_DY(kcolp,mzp ) : Reduced  Potential Temperature   [K/X] |
!                                                                              |
!                  CAPECP(kcolp     ) : Convective Avail.Potent.Energy         |
!                                                                              |
!     REFER. :  MesoNH MASS FLUX CONVECTION Routine                            |
!     ^^^^^^^^         (Bechtold et al., 2001, QJRMS 127, pp 869-886)          |
!                                                                              |
!   # OPTION : #EW  Energy and Water  Conservation                             |
!   # ^^^^^^^^                                                                 |
!                                                                              |
!     MODIF. HGallee: 18-11-2004: Adaptation to CVAmnh.f90.laurent             |
!     ^^^^^^                      (Argument kensbl of CONVECTION removed)      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_CP_ctr
      use Mod_PHY_CP_dat
      use Mod_PHY_CP_grd
      use Mod_PHY_CP_kkl
      use Mod_PHY_DY_kkl
      use Mod_PHY_AT_kkl
      use Mod_PHY_CM_kkl

      use Mod_Atm_CP_RUN

      IMPLICIT NONE



! Global Variables
! ================

      integer                                       ::  mzc                      !  Nb of levels
      integer                                       ::  kcolc                    !  Nb of columns




! Local  Variables
! ================

      real(kind=real8)                              ::  bANA                     !
      real(kind=real8)                              ::  zANA                     !

      integer                                       ::  i                        !  x-Axis Index 
      integer                                       ::  j                        !  y-Axis Index
      integer                                       ::  k                        !  Level  Index (from top    to bottom)
      integer                                       ::  klc                      !  Level  Index (from bottom to top   )
      integer                                       ::  ikl                      !  Column Index

      REAL                                          ::  Pdxdy0(kcolc)
      REAL                                          ::  P_pa_0(kcolc,mzc)
      REAL                                          ::  P_za_0(kcolc,mzc)
      REAL                                          ::  P_Ta_0(kcolc,mzc)
      REAL                                          ::  P_Qa_0(kcolc,mzc)
      REAL                                          ::  P_Qw_0(kcolc,mzc)
      REAL                                          ::  P_Qi_0(kcolc,mzc)
      REAL                                          ::  P_Ua_0(kcolc,mzc)
      REAL                                          ::  P_Va_0(kcolc,mzc)
      REAL                                          ::  P_Wa_0(kcolc,mzc)

      integer                                       ::  Kstep1(kcolc)            ! convective counter
      integer                                       ::  K_CbT1(kcolc)            ! cloud top  level
      integer                                       ::  K_CbB1(kcolc)            ! cloud base level

      integer, parameter                            ::  KTCCH0=1                 !
      REAL                                          ::  P_CH_0(kcolc,mzc,KTCCH0)
      REAL                                          ::  PdCH_1(kcolc,mzc,KTCCH0)

      REAL                                          ::  PdTa_1(kcolc,mzc)
      REAL                                          ::  PdQa_1(kcolc,mzc)
      REAL                                          ::  PdQw_1(kcolc,mzc)
      REAL                                          ::  PdQi_1(kcolc,mzc)
      REAL                                          ::  Pdrr_1(kcolc)
      REAL                                          ::  Pdss_1(kcolc)
      REAL                                          ::  PuMF_1(kcolc,mzc)        !   Upward        Mass Flux
      REAL                                          ::  PdMF_1(kcolc,mzc)        ! Downward        Mass Flux
      REAL                                          ::  Pfrr_1(kcolc,mzc)        ! Liquid Precipitation Flux
      REAL                                          ::  Pfss_1(kcolc,mzc)        !  Solid Precipitation Flux
      REAL                                          ::  Pcape1(kcolc)            !   CAPE [J/kg]


!  Diagnostic Variables
!  --------------------

! #EW integer                                       ::  irmx  ,jrmx  ,iter_0     !
! #EW real(kind=real8)                              ::  rr_max,temp_r,energ0     !
! #EW real(kind=real8)                              ::  water0,waterb            !



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                !
!  ALLOCATION
!  ----------

      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                                      THEN !
          allocate  ( wa_ANA(kcolc,mzc) )                                        !  ANAbatic       Wind        Speed     [m/s]
      END IF
!                                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! Update Convective Mass Flux
! ===========================

      IF     ( mod(iitCV0,jjtCV0).EQ.0 )                            THEN

! Martin control
        PRINT*,'Dans PHY_Atm_CP_RUN'
        PRINT*,'size(ii__AP)',size(ii__AP)
        PRINT*,'size(jj__AP)',size(jj__AP)
        PRINT*,'size(wa_ANA)=',size(wa_ANA)
        PRINT*,'ii__AP(1)=',ii__AP(1)
        PRINT*,'jj__AP(1)=',jj__AP(1)
        PRINT*,'ii__AP(kcolp)=',ii__AP(kcolp)
        PRINT*,'jj__AP(kcolp)=',jj__AP(kcolp)
        PRINT*,'mzp=',mzp
! Martin control

        DO       ikl = 1,kcolp

!           PRINT*,'ikl=',ikl

           i   = ii__AP(ikl)
           j   = jj__AP(ikl)

!          PRINT*,'ii__AP(',ikl,')=',ii__AP(ikl)
!          PRINT*,'jj__AP(',ikl,')=',jj__AP(ikl)


!  Contribution from Subgrid Mountain Breeze
!  -----------------------------------------

            DO k=1,mzp
                  wa_ANA(ikl,k   ) =  0.0
            END DO
          IF(Lo_ANA)                                                THEN
            DO k=1,mzp
                  bANA             =  min(Z___DY(ikl,k  ),  zi__AT(ikl))
                  zANA             =        hANA(ikl)      +  bANA * 2.0
              IF (Z___DY(ikl,k   )   .LE.   zANA             .AND.      &
     &            Ta__DY(ikl,mzpp)   .GT. Ta__DY(ikl,mzp))          THEN
                  wa_ANA(ikl,k   ) =        rANA           *  bANA * 0.5 !  Half Integrated Horizontal Divergence

              END IF
            END DO
          END IF



!  Mass Flux convective Scheme: Set Up Vertical Profiles
!  -----------------------------------------------------

             Pdxdy0(ikl)      = dxHOST * dxHOST                          !  grid area                           [m2]
          DO klc= 1,mzp
             k =    mzpp-klc
             P_pa_0(ikl,klc) = (psa_DY(  ikl)*sigma(k)  + pt__DY) *1.e3  !  pressure     in layer               [Pa]
             P_za_0(ikl,klc) =  Z___DY(ikl,k)                            !  height of model layer                [m]
             P_Ta_0(ikl,klc) =  Ta__DY(ikl,k)                            !  grid scale T           at time t     [K]
             P_Qa_0(ikl,klc) =  qv__DY(ikl,k)                            !  grid scale water vapor at time t [kg/kg]
             P_Qw_0(ikl,klc) =  qw__CM(ikl,k) / (1.0-qw__CM(ikl,k))      !  grid scale Cloud drops at time t [kg/kg]
             P_Qi_0(ikl,klc) =  qi__CM(ikl,k) / (1.0-qi__CM(ikl,k))      !  grid scale Cloud ice   at time t [kg/kg]
             P_Ua_0(ikl,klc) =  ua__DY(ikl,k)                            !  grid scale hor. wind u at time t   [m/s]
             P_Va_0(ikl,klc) =  va__DY(ikl,k)                            !  grid scale hor. wind v at time t   [m/s]
             P_Wa_0(ikl,klc) =  wa__DY(ikl,k) +      wa_ANA(ikl,k)       !  grid scale vertic.wind at time t   [m/s]
          END DO

        END DO ! ikl = 1,kcolp


!  Mass Flux convective Scheme: Bechtold et al. (2001) Convective Parameterization
!  -------------------------------------------------------------------------------

!         ***************
          call CONVECTION(                                               &
     &             kcolc , mzc   , kidia0, kfdia0, kbdia0, ktdia0,       &
     &             pdtCV , Odeep , Oshal , Orset0, Odown0, kIce_0,       &
     &             OsetA0, PTdcv , PTscv ,                               &
     &             kensbl,                                               &
     &     P_pa_0, P_za_0, Pdxdy0,                                       &
     &     P_Ta_0, P_Qa_0, P_Qw_0, P_Qi_0, P_Ua_0, P_Va_0, P_Wa_0,       &
     &     Kstep1, PdTa_1, PdQa_1, PdQw_1, PdQi_1,                       &
     &                     Pdrr_1, Pdss_1,                               &
     &     PuMF_1, PdMF_1, Pfrr_1, Pfss_1, Pcape1, K_CbT1, K_CbB1,       &
     &     OCvTC0, KTCCH0, P_CH_0, PdCH_1)
!         ***************



!  Mass Flux convective Scheme: products
!  -------------------------------------

        DO       ikl = 1,kcolp
             CAPECP(ikl)   = Pcape1(ikl)
             timeCP(ikl)   = Kstep1(ikl)     *dt__CP
             drr_CP(ikl)   = Pdrr_1(ikl)
             dss_CP(ikl)   = Pdss_1(ikl)

          DO klc= 1,mzp
             k  =   mzpp - klc
             dpktCP(ikl,k) = PdTa_1(ikl,klc) /ExnrDY(ikl,k)
             dqv_CP(ikl,k) = PdQa_1(ikl,klc)
             dqw_CP(ikl,k) = PdQw_1(ikl,klc)
             dqi_CP(ikl,k) = PdQi_1(ikl,klc)

          END DO

        END DO


      END IF ! mod(iitCV0,jjtCV0).EQ.0   (UPDATE)                   THEN




! Vertical Integrated Energy and Water Content
! ============================================

! #EW     irmx   =       i_x0
! #EW     jrmx   =       j_y0
! #EW     rr_max =       0.0

        DO       ikl = 1,kcolp
           i   = ii__AP(ikl)
           j   = jj__AP(ikl)

! #EW     enr0EW(ikl) = 0.0
! #EW     wat0EW(ikl) = 0.0

! #EW    DO k=1,mzp
! #EW     temp_r      = pkt_DY(ikl,k)*ExnrDY(ikl,k)
! #EW     enr0EW(ikl) = enr0EW(ikl)                                      &
! #EW&                +(temp_r                                           &
! #EW&                -(qw__CM(ikl,k)+qr__CM(ikl,k)) * Lv_CPd            &
! #EW&                -(qi__CM(ikl,k)+qs__CM(ikl,k)) * Ls_CPd)*dsigmi(k)
! #EW     wat0EW(ikl) = wat0EW(ikl)                                      &
! #EW&                +(qv__DY(ikl,k)                                    &
! #EW&                + qw__CM(ikl,k)+qr__CM(ikl,k)                      &
! #EW&                + qi__CM(ikl,k)+qs__CM(ikl,k)          )*dsigmi(k)
! #EW    END DO

! #EW     enr0EW(ikl) = enr0EW(ikl) * psa_DY(  ikl)  * Grav_I
! #EW     wat0EW(ikl) = wat0EW(ikl) * psa_DY(  ikl)  * Grav_I
!  ..     wat0EW [m]    contains implicit factor 1.d3 [kPa-->Pa] /ro_Wat

! #EW     energ0      = energ0 - enr0EW(ikl)
! #EW     water0      = water0 - wat0EW(ikl)




! Update of Mass Flux convective Tendencies
! =========================================

!       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         IF                    (timeCP(ikl) .GT. 0.)                THEN


!  Temporal tendencies on pkt_DY, qv__DY and rainCM
!  ------------------------------------------------

          DO k=1,mzp
            dqv_CP(ikl,k) = min(dqv_CP(ikl,k),(qv__DY(ikl,k)-epsq)/dt__CP)
            dqw_CP(ikl,k) = min(dqw_CP(ikl,k),(qw__CM(ikl,k)-epsq)/dt__CP)
            dqi_CP(ikl,k) = min(dqi_CP(ikl,k),(qi__CM(ikl,k)-epsq)/dt__CP)
          ENDDO

            rainCM(ikl)   =     rainCM(ikl)  + drr_CP(ikl)        *dt__CP   
            rainCP(ikl)   =     rainCP(ikl)  + drr_CP(ikl)        *dt__CP   

            snowCM(ikl)   =     snowCM(ikl)  + dss_CP(ikl)        *dt__CP   
            snowCP(ikl)   =     snowCP(ikl)  + dss_CP(ikl)        *dt__CP   

         ELSE

          DO k=1,mzp
            dpktCP(ikl,k) =     0.00
            dqv_CP(ikl,k) =     0.00
            dqw_CP(ikl,k) =     0.00
            dqi_CP(ikl,k) =     0.00
          ENDDO


         ENDIF          !      {timeCP(ikl) .gt. 0}

            timeCP(ikl)   = max(timeCP(ikl) - dt__CP,zer0)
!  ..       ^^^^ remaining time before the end of convection
!       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




! Vertical Integrated Energy and Water Content
! ============================================

! #EW     enr1EW(ikl) = 0.0
! #EW     wat1EW(ikl) = 0.0
! #EW     watfEW(ikl) =-drr_CP(ikl)

! #EW    DO k=1,mz
! #EW     temp_r      = pkt_DY(ikl,k)*ExnrDY(ikl,k)
! #EW     enr1EW(ikl) = enr1EW(ikl)                                      &
! #EW&                +(temp_r                                           &
! #EW&                -(qw__CM(ikl,k)+qr__CM(ikl,k)) *Lv_CPd             &
! #EW&                -(qi__CM(ikl,k)+qs__CM(ikl,k)) *Ls_CPd)*dsigmi(k)
! #EW     wat1EW(ikl) = wat1EW(ikl)                                      &
! #EW&                +(qv__DY(ikl,k)                                    &
! #EW&                + qw__CM(ikl,k)+qr__CM(ikl,k)                      &
! #EW&                + qi__CM(ikl,k)+qs__CM(ikl,k)         )*dsigmi(k)
! #EW    END DO

! #EW     enr1EW(ikl) = enr1EW(ikl) * psa_DY(  ikl)  *Grav_I             &
! #EW&                - drr_CP(ikl)                  *Lv_CPd
! #EW     wat1EW(ikl) = wat1EW(ikl) * psa_DY(  ikl)  *Grav_I
!  ..     wat1EW [m]    contains implicit factor 1.d3 [kPa-->Pa] /rhoWat

! #EW     energ0      = energ0 + enr1EW(ikl)
! #EW     water0      = water0 + wat1EW(ikl)
! #EW     iter_0      = iter_0 + 1




! Vertical Integrated Energy and Water Content: OUTPUT
! ====================================================

! #EW   IF (drr_CP(ikl).gt.rr_max)                                THEN
! #EW       rr_max =       drr_CP(ikl)
! #EW       irmx   =              i
! #EW       jrmx   =                j
! #EW   END IF

        END DO

! #EW   waterb = wat1EW(irmx,jrmx)-wat0EW(irmx,jrmx)-watfEW(irmx,jrmx)
! #EW   write(6,606) it_EXP,enr0EW(irmx,jrmx),1.d3*wat0EW(irmx,jrmx),    &
! #EW&            irmx,jrmx,enr1EW(irmx,jrmx),1.d3*wat1EW(irmx,jrmx),    &
! #EW&                                        1.d3*watfEW(irmx,jrmx),    &
! #EW&                                        1.d3*waterb           ,    &
! #EW&                      energ0/iter_0    ,     water0/iter_0
 606    format(i9,'  Before CVAj:  E0 =',f12.6,'  W0 = ',f9.6,           &
     &    /,i5,i4,'  After  CVAj:  E1 =',f12.6,'  W1 = ',f9.6,           &
     &                                         '  W Flux =',f9.6,        &
     &                                         '  Div(W) =',e9.3,        &
     &       /,9x,'         Mean   dE =',f12.9,'  dW = ',e9.3)




! Incrementation Step Nb
! ======================

      iitCV0 = iitCV0 + 1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                !
!  DE-ALLOCATION
!  =============

      IF (FlagDALLOC)                                      THEN !
          deallocate  ( wa_ANA )                                                 !  ANAbatic       Wind        Speed     [m/s]
      END IF

! Bug corrigé par Martin:
!      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                                      THEN !
!          deallocate  ( wa_ANA )                                                 !  ANAbatic       Wind        Speed     [m/s]
!      END IF
!                                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      return
      end subroutine PHY_Atm_CP_RUN
