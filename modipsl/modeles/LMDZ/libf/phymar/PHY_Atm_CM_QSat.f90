      subroutine PHY_Atm_CM_QSat

!------------------------------------------------------------------------------+
!                                                         Mon 17-Jun-2013  MAR |
!                                                                              |
!     subroutine PHY_Atm_CM_QSat   computes specifics Humidities at Saturation |
!            for Cloud Microphysical    Scheme used by MAR                     |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Sat 23-Mar-2013      |
!           Last Modification by H. Gallee,               Mon 17-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+
!                                                                              |
!     INPUT :   Ta__DY:         Air Temperature                            [K] |
!     ^^^^^     psa_DY: Model Pressure Thickness                         [kPa] |
!                                                                              |
!     OUTPUT :  qvswCM     : Saturation Specific Humidity  over Water  [kg/kg] |
!     ^^^^^^    qvsiCM     : Saturation Specific Humidity  over Ice    [kg/kg] |
!                                                                              |
!     REFERENCE :  Dudhia (1989) JAS, (B1) and (B2) p.3103                     |
!     ^^^^^^^^^                                                                |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_DY_kkl
      use Mod_PHY_CM_kkl




! Local  Variables
! ================

      use Mod_Atm_CM_QSa


      IMPLICIT NONE


      integer i     ,j     ,ikl   ,k

      real(kind=real8)                             ::  WatIce =  273.16e0
      real(kind=real8)                             ::  ExpWat =    5.138e0
      real(kind=real8)                             ::  ExpWa2 = 6827.e0
      real(kind=real8)                             ::  ExpIce = 6150.e0
      real(kind=real8)                             ::  pr__75 =   75.00        !                                      75 [hPa]
      real(kind=real8)                             ::  qs__75 =    0.001       !                                 0.001 [kg/kg]

      real(kind=real8)                             ::  Wphase
      real(kind=real8)                             ::  qvSatI




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
! ALLOCATION
! ==========

      IF (it_RUN.EQ.1 .OR. FlagDALLOC)               THEN !
          allocate   ( pa_hPa(mzpp) )
          allocate   ( pr_b75(mzpp) )
          allocate   ( ei_sat(mzpp) )
          allocate   ( ew_sat(mzpp) )
      END IF
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!     ++++++++++++++++
      DO ikl = 1,kcolp
!     ++++++++++++++++



! Temperature (K) and Pressure (hPa)
! ==================================

        DO k=1,mzpp
          pa_hPa(k)   = (psa_DY(ikl) * sigma(k) + pt__DY) * 10.0d0


! Pressure Discriminator
! ----------------------

          pr_b75(k)   = max(zer0,sign(un_1,pa_hPa(k)-pr__75))




! Saturation Vapor Pressure over Ice
! ==================================

          ei_sat(k)    =                                               &
     &    6.1070d0 * exp (ExpIce *(un_1/WatIce-un_1/Ta__DY(ikl,k) ) )
!  ..     Dudhia (1989) JAS, (B1) and (B2) p.3103 

          qvSatI       = .622d0*ei_sat(k)                              &
     &  /(pa_hPa(k)    - .378d0*ei_sat(k) )




! Saturation Vapor Pressure over Water
! ====================================

          ew_sat(k)    =                                               &
     &    6.1078d0 * exp (ExpWat *  log(WatIce     /Ta__DY(ikl,k) ) )  &
     &             * exp (ExpWa2 *(un_1/WatIce-un_1/Ta__DY(ikl,k) ) )
!  ..     Dudhia (1989) JAS, (B1) and (B2) p.3103 
!         See also Pielke (1984), p.234 and Stull (1988), p.276 

          qvswCM(ikl,k) = max(epsn  ,  .622d0*ew_sat(k)                &
     &  /(pa_hPa(k) -                  .378d0*ew_sat(k)   ))
!  ..     Saturation Vapor Specific Concentration over Water
!         (even for temperatures less than freezing point)




! Water Phase Discriminator
! =========================

          Wphase       = max(zer0,sign(un_1,Ta__DY(ikl,k)-WatIce))
!  ..     Wphase       =     1    if        Tair     >    273.16
!                            0    if        Tair     <    273.16




! Saturation Vapor Specific Concentration over Ice
! ================================================

          qvsiCM(ikl,k) =                                              &
     &       max(epsn    , qvswCM(ikl,k) *      Wphase                 &
     &                   + qvSatI        *(un_1-Wphase))


          qvswCM    (ikl,k) = pr_b75(k) *qvswCM(ikl,k) +(1.0-pr_b75(k)) *qs__75
          qvsiCM    (ikl,k) = pr_b75(k) *qvsiCM(ikl,k) +(1.0-pr_b75(k)) *qs__75


        END DO

!     ++++++++++++++++++++++
      END DO ! ikl = 1,kcolp
!     ++++++++++++++++++++++



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
! DE-ALLOCATION
! =============

      IF (FlagDALLOC)                                THEN !
          deallocate ( pa_hPa )
          deallocate ( pr_b75 )
          deallocate ( ei_sat )
          deallocate ( ew_sat )
      END IF
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      return
      end subroutine PHY_Atm_CM_QSat
