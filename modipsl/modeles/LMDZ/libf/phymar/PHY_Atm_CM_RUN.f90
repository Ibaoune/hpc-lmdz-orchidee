      subroutine PHY_Atm_CM_RUN

!------------------------------------------------------------------------------+
!                                                         Sun  9-Jun-2013  MAR |
!   MAR          PHY_Atm_CM_RUN                                                |
!     subroutine PHY_Atm_CM_RUN  drives  Cmoud Microphysical Scheme CMiPhy     |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Thu 21-Mar-2013      |
!           Last Modification by H. Gallee,               Sun  9-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+
!                                                                              |
!   INPUT:          psa_DY : Pressure       Thickness                    [kPa] |
!   ^^^^^^                                                                     |
!                   qwd_CM : Cloud Droplets Concentr.Var.              [kg/kg] |
!                   qwi_CM : Ice   Crystals Concentr.Var.              [kg/kg] |
!                                                                              |
!   INPUT / OUTPUT: qv__DY : Air   specific Humidity                   [kg/kg] |
!   ^^^^^^^^^^^^^^^ qw__CM : Cloud Droplets Concentration              [kg/kg] |
!                   qi__CM : Ice   Crystals Concentration              [kg/kg] |
!                   qs__CM : Snow  Particl. Concentration              [kg/kg] |
!                   qr__CM : Rain  Drops    Concentration              [kg/kg] |
!                   CCNiCM : ice   crystals number                     [Nb/m3] |
!                  (CCNwCM : Cloud Condens. Nuclei(if #cw)             [Nb/m3])|
!                                                                              |
!   OUTPUT (CMiPhy) RainCM : rain  Precipitation                           [m] |
!   ^^^^^^^^^^^^^^^ SnowCM : snow  Precipitation                       [m w.e] |
!                   Ice_CM : ice   Precipitation                       [m w.e] |
!                                                                              |
!   OUTPUT:         qwd_CM : Cloud Droplets Concentr.Var.              [kg/kg] |
!   ^^^^^^^         qid_CM : Ice   Crystals Concentr.Var.              [kg/kg] |
!                                                                              |
!                   HLatCM : Latent Heat Release                         [K/s] |
!                                                                              |
!   REFER. : 1) Ntezimana, unpubl.thes.LLN, 115 pp,     1993                   |
!   ^^^^^^^^ 2) Lin et al.       JCAM   22, 1065--1092, 1983                   |
!                 (very similar, except that graupels are represented)         |
!              3) Emde and Kahlig, An.Geo. 7,  405-- 414, 1989                 |
!                                                                              |
! # OPTIONS: #qg  Graupels                  (qg)   Microphysics Activation     |
! # ^^^^^^^^ #cw  Cloud Condensation Nuclei (CCNw) Microphysics Activation     |
!                                                                              |
!            #qd  Variations in Time of qi and qw are stored for OUTPUT        |
!                                                                              |
!            #WH  Additional Output (Each Process  is detailed)                |
!            #EW  Additional Output (Energy and Water Conservation)            |
!                                                                              |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_CM_ctr
      use Mod_PHY_CM_dat
      use Mod_PHY_CM_grd
      use Mod_PHY_CM_kkl
      use Mod_PHY_DY_kkl




!  Local Variables
!  ===============

      use Mod_Atm_CM_RUN


      IMPLICIT NONE


      integer                                      ::  i  ,j  ,k  ,ikl  !
      real(kind=real8)                             ::  Hcd,Hsb,Tcd,Tsb  !
      real(kind=real8)                             ::  Zcd,Zsb,facLHR   !





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!  ALLOCATION                                                           !
!  ==========

      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                             THEN !
          allocate ( qv__00(kcolp,mzp) )                                !  Air   Specific Humidity  before Cloud Microphys.  [kg/kg]
          allocate ( qw__00(kcolp,mzp) )                                !  Cloud Droplets Concentr. before Cloud Microphys.  [kg/kg]
          allocate ( qi__00(kcolp,mzp) )                                !  Cloud Crystals Concentr. before Cloud Microphys.  [kg/kg]
          allocate ( qs__00(kcolp,mzp) )                                !  Snow  Particl. Concentr. before Cloud Microphys.  [kg/kg]
! #qg     allocate ( qg__00(kcolp,mzp) )                                !  Graupels       Concentr. before Cloud Microphys.  [kg/kg]
          allocate ( qr__00(kcolp,mzp) )                                !  Rain  Drops    Concentr. before Cloud Microphys.  [kg/kg]
          allocate ( CFra00(kcolp,mzp) )                                !  Cloud Fraction           before Cloud Microphys.      [-]
! #cw     allocate ( CCNw00(kcolp,mzp) )                                !  Cloud Condens. Nuclei    before Cloud Microphys.      [-]
          allocate ( CCNi00(kcolp,mzp) )                                !  Cloud Ice      Nuclei    before Cloud Microphys.      [-]
      END IF                                                            !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!  Microphysical Variables: from 3D to 2D Arrays
!  =======================  ====================

        DO   ikl = 1,kcolp

             i = ii__AP(ikl)
             j = jj__AP(ikl)

          DO k = 1,mzp
             qv__00(ikl,k) = qv__DY(ikl,k)
             qv__DY(ikl,k) =                                           &
     &                   max(qv__DY(ikl,k) , epsq)
             Ta__CM(ikl,k) = pkt_DY(ikl,k) * ExnrDY(ikl,k)

          IF(qw__CM(ikl,k) .LT. qh_MIN)                             THEN
             qw__CM(ikl,k) =                    max(zer0,qw__CM(ikl,k))   ! Sinon BOUM (possible)
             Ta__CM(ikl,k) =    Ta__CM(ikl,k) - Lv_Cpd * qw__CM(ikl,k)
             qv__DY(ikl,k) =    qv__DY(ikl,k) +          qw__CM(ikl,k)
             qw__CM(ikl,k) =    0.
          END IF
             qw__00(ikl,k) =    qw__CM(ikl,k)

          IF(qi__CM(ikl,k) .LT. qh_MIN)                             THEN
             qi__CM(ikl,k) =                    max(zer0,qi__CM(ikl,k))   ! Sinon BOUM (possible)
             Ta__CM(ikl,k) =    Ta__CM(ikl,k) - Ls_Cpd * qi__CM(ikl,k)
             qv__DY(ikl,k) =    qv__DY(ikl,k) +          qi__CM(ikl,k)
             qi__CM(ikl,k) =    0.
             CCNiCM(ikl,k) =    0.
          END IF
             qi__00(ikl,k) =    qi__CM(ikl,k)

          IF(qw__CM(ikl,k) .LT. qh_MIN  .AND.                          &
     &       qi__CM(ikl,k) .LT. qh_MIN)                             THEN
             CFraCM(ikl,k) =    0.
          ELSE
             CFraCM(ikl,k) =max(CFrMIN,CFraCM(ikl,k))
          END IF
             CFra00(ikl,k) =    CFraCM(ikl,k)

          IF(qr__CM(ikl,k) .LT. qh_MIN)                             THEN
             qr__CM(ikl,k) =                    max(zer0,qr__CM(ikl,k))   ! Sinon BOUM (possible)
             Ta__CM(ikl,k) =    Ta__CM(ikl,k) - Lv_Cpd * qr__CM(ikl,k)
             qv__DY(ikl,k) =    qv__DY(ikl,k) +          qr__CM(ikl,k)
             qr__CM(ikl,k) =    0.
          END IF
             qr__00(ikl,k) = qr__CM(ikl,k)

          IF(qs__CM(ikl,k) .LT. qh_MIN)                             THEN
             qs__CM(ikl,k) =                    max(zer0,qs__CM(ikl,k))   ! Sinon BOUM (possible)
             Ta__CM(ikl,k) =    Ta__CM(ikl,k) - Ls_Cpd * qs__CM(ikl,k)
             qv__DY(ikl,k) =    qv__DY(ikl,k) +          qs__CM(ikl,k)
             qs__CM(ikl,k) =    0.
          END IF
             qs__00(ikl,k) =    qs__CM(ikl,k)

! #qg        qg__00(ikl,k) =    qg__CM(ikl,k)

! #cw        CCNw00(ikl,k) =    CCNwCM(ikl,k)
             CCNi00(ikl,k) =    CCNiCM(ikl,k)

             qid_CM(ikl,k) =    0.
             qwd_CM(ikl,k) =    0.
          END DO




!  Vertical Integrated Energy and Water Content
!  ============================================

! #EW        wat01D(ikl)  =wat0EW(ikl)
! #EW        wat11D(ikl)  =wat1EW(ikl)
! #EW        wat21D(ikl)  =wat2EW(ikl)
! #EW        watf1D(ikl)  =watfEW(ikl)
! #EW        enr01D(ikl)  =enr0EW(ikl)
! #EW        enr11D(ikl)  =enr1EW(ikl)
! #EW        enr21D(ikl)  =enr2EW(ikl)
! #EW        mphy2D(ikl)  =mphyEW(ikl)


        END DO
                                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!  Call Cloud Microphysics
!  =======================

!               ******
           call CMiPhy
!               ******




!  Microphysical Variables: from 3D to 2D Arrays
!  =======================  ====================

        DO   ikl = 1,kcolp

             i = ii__AP(ikl)
             j = jj__AP(ikl)

            DO k = 1,mzp
             dpktCM(ikl,k) = (Ta__CM(ikl,k) / ExnrDY(ikl,k)            &
     &                                      - pkt_DY(ikl,k)) / dt__CM
             dqv_CM(ikl,k) = (qv__DY(ikl,k) - qv__00(ikl,k)) / dt__CM
                                  qv__DY(ikl,k) = qv__00(ikl,k)
! #qd        qid_CM(ikl,k) =  qid_CM(ikl,k) * psa_DY(ikl)    * dsigmi(k)
! #qd        qwd_CM(ikl,k) =  qwd_CM(ikl,k) * psa_DY(ikl)    * dsigmi(k)
            END DO



! Update of Microphysical Variables:  INSIDE PHY_MAR
! --------------------------------------------------

          IF    (CM_UpD)                                           THEN
            DO k = 1,mzp

             dqw_CM(ikl,k) =  0.0000
             dqr_CM(ikl,k) =  0.0000
             dqi_CM(ikl,k) =  0.0000
             dqs_CM(ikl,k) =  0.0000
! #qg        dqg_CM(ikl,k) =  0.0000

             dCF_CM(ikl,k) =  0.0000
! #cw        dCw_CM(ikl,k) =  0.0000
             dCi_CM(ikl,k) =  0.0000
            END DO



! Update of Microphysical Variables: OUTSIDE PHY_MAR
! --------------------------------------------------

          ELSE
            DO k = 1,mzp
             dqw_CM(ikl,k) = (qw__CM(ikl,k) -qw__00(ikl,k)) / dt__CM
             dqr_CM(ikl,k) = (qr__CM(ikl,k) -qr__00(ikl,k)) / dt__CM
             dqi_CM(ikl,k) = (qi__CM(ikl,k) -qi__00(ikl,k)) / dt__CM
             dqs_CM(ikl,k) = (qs__CM(ikl,k) -qs__00(ikl,k)) / dt__CM
! #qg        dqg_CM(ikl,k) = (qg__CM(ikl,k) -qg__00(ikl,k)) / dt__CM
             dCF_CM(ikl,k) = (CFraCM(ikl,k) -CFra00(ikl,k)) / dt__CM
! #cw        dCw_CM(ikl,k) = (CCNwCM(ikl,k) -CCNw00(ikl,k)) / dt__CM
             dCi_CM(ikl,k) = (CCNiCM(ikl,k) -CCNi00(ikl,k)) / dt__CM

             qw__CM(ikl,k) =  qw__00(ikl,k)
             qr__CM(ikl,k) =  qr__00(ikl,k)
             qi__CM(ikl,k) =  qi__00(ikl,k)
             qs__CM(ikl,k) =  qs__00(ikl,k)
! #qg        qg__CM(ikl,k) =  qg__00(ikl,k)
! Gilles : dCF_CM pas utilise pour mettre a jour CFraCM
!            CFraCM(ikl,k) =  CFra00(ikl,k)
! #cw        CCNwCM(ikl,k) =  CCNw00(ikl,k)
             CCNiCM(ikl,k) =  CCNi00(ikl,k)

            END DO
          END IF





!  Isotopes Proxies: Diagnostics
!  =============================

          IF (nHL_CM.EQ.0)                                          THEN
              Hcd_CM(ikl) = 0.
              Hsb_CM(ikl) = 0.
              Tcd_CM(ikl) = 0.
              Tsb_CM(ikl) = 0.
              Zcd_CM(ikl) = 0.
              Zsb_CM(ikl) = 0.
          END IF
              nHL_CM      = nHL_CM     + 1

              Hcd =  0.0
              Hsb =  0.0
              Tcd =  0.0
              Tsb =  0.0
              Zcd =  0.0
              Zsb =  0.0
          DO k=2,mzp
             HLatCM(ikl,k) = (Ta__CM(ikl,k)                            &
     &                       -pkt_DY(ikl,k) *ExnrDY(ikl,k)) / dt__CM

              Hcd =  Hcd + dsigmi(k)*max(HLatCM(ikl,k),0.)
              Hsb =  Hsb - dsigmi(k)*min(HLatCM(ikl,k),0.)
              Tcd =  Tcd + dsigmi(k)*max(HLatCM(ikl,k),0.)*Ta__CM(ikl,k)
              Tsb =  Tsb - dsigmi(k)*min(HLatCM(ikl,k),0.)*Ta__CM(ikl,k)
              Zcd =  Zcd + dsigmi(k)*max(HLatCM(ikl,k),0.)*Z___DY(ikl,k)
              Zsb =  Zsb - dsigmi(k)*min(HLatCM(ikl,k),0.)*Z___DY(ikl,k)
          END DO

              facLHR = (CpdAir/LhsH2O)*psa_DY(ikl)   *1.e3*Grav_I*dt__CM

          IF (write_Proxy)                                          THEN
              Hcd_CM(ikl) =(Hcd_CM(ikl)  + Hcd * facLHR) / nHL_CM
              Hsb_CM(ikl) =(Hsb_CM(ikl)  + Hsb * facLHR) / nHL_CM
              Tcd_CM(ikl) =(Tcd_CM(ikl)  + Tcd * facLHR) / nHL_CM
              Tsb_CM(ikl) =(Tsb_CM(ikl)  + Tsb * facLHR) / nHL_CM
              Zcd_CM(ikl) =(Zcd_CM(ikl)  + Zcd * facLHR) / nHL_CM
              Zsb_CM(ikl) =(Zsb_CM(ikl)  + Zsb * facLHR) / nHL_CM
          ELSE
              Hcd_CM(ikl) = Hcd_CM(ikl)  + Hcd * facLHR    
              Hsb_CM(ikl) = Hsb_CM(ikl)  + Hsb * facLHR    
              Tcd_CM(ikl) = Tcd_CM(ikl)  + Tcd * facLHR     
              Tsb_CM(ikl) = Tsb_CM(ikl)  + Tsb * facLHR     
              Zcd_CM(ikl) = Zcd_CM(ikl)  + Zcd * facLHR     
              Zsb_CM(ikl) = Zsb_CM(ikl)  + Zsb * facLHR     
          END IF

        END DO

          IF (write_Proxy)  nHL_CM = 0




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!  DE-ALLOCATION                                                        !
!  =============                                                        !

      IF (FlagDALLOC)                                              THEN !
          deallocate ( qv__00 )                                         !  Air   Specific Humidity  before Cloud Microphys.  [kg/kg]
          deallocate ( qw__00 )                                         !  Cloud Droplets Concentr. before Cloud Microphys.  [kg/kg]
          deallocate ( qi__00 )                                         !  Cloud Crystals Concentr. before Cloud Microphys.  [kg/kg]
          deallocate ( qs__00 )                                         !  Snow  Particl. Concentr. before Cloud Microphys.  [kg/kg]
! #qg     deallocate ( qg__00 )                                         !  Graupels       Concentr. before Cloud Microphys.  [kg/kg]
          deallocate ( qr__00 )                                         !  Rain  Drops    Concentr. before Cloud Microphys.  [kg/kg]
          deallocate ( CFra00 )                                         !  Cloud Fraction           before Cloud Microphys.      [-]
! #cw     deallocate ( CCNw00 )                                         !  Cloud Condens. Nuclei    before Cloud Microphys.      [-]
          deallocate ( CCNi00 )                                         !  Cloud Ice      Nuclei    before Cloud Microphys.      [-]
      END IF                                                            !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      return

      end subroutine PHY_Atm_CM_RUN
