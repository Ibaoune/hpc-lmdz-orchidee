      subroutine PHY_Atm_AT_RUN(FlagSV_KzT,FlagCM)

!------------------------------------------------------------------------------+
!                                                         Sat 22-Jun-2013  MAR |
!   MAR          PHY_Atm_AT_RUN                                                |
!     subroutine PHY_Atm_AT_RUN performs   Turbulent Vertical Diffusion Scheme |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Wed 13-Mar-2013      |
!           Last Modification by H. Gallee,               Sat 22-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_AT_ctr
      use Mod_PHY_AT_grd
      use Mod_PHY_AT_kkl
      use Mod_SISVAT_gpt


      IMPLICIT NONE


      logical    ::    FlagSV_KzT
      logical    ::    FlagCM

      integer    ::    i,     j,     k,     ikl




! Assignation    of Mod_PHY_AT_kkl
! ================================

        DO ikl = 1,kcolp

          i = ii__AP(ikl)
          j = jj__AP(ikl)

          LMO_AT(ikl)   = LMO_SV_gpt(ikl)

        ENDDO




! Turbulent Kinetic Energy
! ========================

      IF (AT_TKE)                                                   THEN

!              **************
        CALL   PHY_vdfTKE_RUN
        CALL   PHY_genTKE_RUN
!              **************

      ELSE
         DO ikl = 1,kcolp
         DO   k = 1,mzp
           TKE_AT(ikl,k)=    eps6
           eps_AT(ikl,k)=    epsn
         ENDDO
         ENDDO
      END IF




! Turbulent Vertical Diffusion
! ============================
 
          NewAAT =  .TRUE.
          SchmAT =  'u'    
!                    ******
          CALL       Atm_AT
!                    ******


          NewAAT =  .FALSE.
          SchmAT =  'v'    
!                    ******
          CALL       Atm_AT
!                    ******


      IF (.NOT.FlagSV_KzT)                                          THEN

          SchmAT =  'T'    
!                    ******
          CALL       Atm_AT
!                    ******
      ELSE

          DO k  =1,mzp
          DO ikl=1,kcolp
            dpktAT(ikl,k) = 0.0000
          ENDDO
          ENDDO

      END IF


          SchmAT =  'q'    
!                    ******
          CALL       Atm_AT
!                    ******


      IF (FlagCM)                                                   THEN

          SchmAT =  'w'    
!                    ******
          CALL       Atm_AT
!                    ******


          SchmAT =  'i'    
!                    ******
          CALL       Atm_AT
!                    ******


          SchmAT =  'C'    
!                    ******
          CALL       Atm_AT
!                    ******


          SchmAT =  's'    
!                    ******
          CALL       Atm_AT
!                    ******


          SchmAT =  'r'    
!                    ******
          CALL       Atm_AT
!                    ******


      END IF


      end subroutine PHY_Atm_AT_RUN




      subroutine Atm_AT

!------------------------------------------------------------------------------+
!                                                         Sat 22-Jun-2013  MAR |
!   MAR          Atm_AT                                                        |
!     subroutine Atm_AT solves Turbulent Vertical Diffusion Scheme             |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Wed 13-Mar-2013      |
!           Last Modification by H. Gallee,               Sat 22-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_AT_grd
      use Mod_PHY_AT_kkl
      use Mod_PHY_CM_kkl
      use Mod_PHY_DY_kkl
      use Mod_SISVAT_gpt



      IMPLICIT NONE



      integer           ::       i,     j,     k,     ikl
      real(kind=real8)  ::  psa_sq




!     ++++++++++++++++
      DO ikl = 1,kcolp
!     ++++++++++++++++

        i = ii__AP(ikl)
        j = jj__AP(ikl)



! Tri-Diagonal Matrix: Variable
! =============================

!  Wind Speed
!  ----------

        IF      (schmAT.EQ.'a')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = 1.0
            var_AT(ikl,k) = ua__DY(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = 0.

        ELSE IF (schmAT.EQ.'b')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = 1.0
            var_AT(ikl,k) = va__DY(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = 0.

        ELSE IF (schmAT.EQ.'u')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = Kzm_AT(ikl,k)
            var_AT(ikl,k) = ua__DY(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = 0.

        ELSE IF (schmAT.EQ.'v')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = Kzm_AT(ikl,k)
            var_AT(ikl,k) = va__DY(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = 0.


!  Temperature
!  -----------

        ELSE IF (schmAT.EQ.'T')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = Kzh_AT(ikl,k)
            var_AT(ikl,k) = pkt_DY(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = pkt_DY(ikl,k)


!  Water Species
!  -------------

        ELSE IF (schmAT.EQ.'q')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = Kzh_AT(ikl,k)
            var_AT(ikl,k) = qv__DY(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = qv__DY(ikl,mzp)                            &
     &                    -(Z___DY(ikl,mzp)-Z___DY(ikl,mzpp))          &
     &                    * uqs_SV_gpt(ikl)    /Kzh_AT(ikl,mzp )
            var_AT(ikl,k) =            max(epsq,var_AT(ikl,k   ))
            var_AT(ikl,k) =            min(     var_AT(ikl,k   ),qvsiCM(ikl,k))
            qv__DY(ikl,mzpp) =                  var_AT(ikl,k   )

        ELSE IF (schmAT.EQ.'w')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = Kzh_AT(ikl,k)
            var_AT(ikl,k) = qw__CM(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = qw__CM(ikl,mzp)

        ELSE IF (schmAT.EQ.'i')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = Kzh_AT(ikl,k)
            var_AT(ikl,k) = qi__CM(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = qi__CM(ikl,mzp)

        ELSE IF (schmAT.EQ.'C')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = Kzh_AT(ikl,k)
            var_AT(ikl,k) = CCNiCM(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = CCNiCM(ikl,mzp)

        ELSE IF (schmAT.EQ.'s')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = Kzh_AT(ikl,k)
            var_AT(ikl,k) = qs__CM(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = qs__CM(ikl,k)

        ELSE IF (schmAT.EQ.'r')                                     THEN
          DO k=1,mzp
            Kz__AT(ikl,k) = Kzh_AT(ikl,k)
            var_AT(ikl,k) = qr__CM(ikl,k)
          ENDDO
             k=  mzpp
            var_AT(ikl,k) = qr__CM(ikl,mzp)

        ELSE
            write(6,*) 
            write(6,*) ' schmAT = ',schmAT
            write(6,*) ' Inadequate..., NOR for Ekman Spiral, NOR for momentum, NOR for scalars'
            write(6,*) '                                                    '

!                **************
            call MAR________BUG
!                **************
            STOP
        END IF



! Tri-Diagonal Matrix Coefficients
! ================================

        IF (NewAAT)                                                 THEN
            psa_sq          = psa_DY(ikl)     *psa_DY(ikl)    
            Ac__AT(ikl,mzp) = Ac0_AT(mzp)                                   &
     &                      * roa_DY(ikl,mzp) *roamDY(ikl,mzp)              &
     &                      / psa_sq
          DO k=mzp-1,1,-1
            Ac__AT(ikl,k  ) = Ac0_AT(k  )                                   &
     &                      * roa_DY(ikl,k  ) *roamDY(ikl,k  )              &
     &                      / psa_sq
            Cc__AT(ikl,k+1) = Ac__AT(ikl,k)                                 &
     &                      * roa_DY(ikl,k+1) /roa_DY(ikl,k  )
          ENDDO
        ENDIF

            A___AT(ikl,mzp) =           Ac__AT(ikl,mzp) * Kz__AT(ikl,mzp)
          DO k=mzp-1,1,-1
            A___AT(ikl,k  ) =           Ac__AT(ikl,k  ) * Kz__AT(ikl,k  )
            C___AT(ikl,k+1) =           Cc__AT(ikl,k+1) * Kz__AT(ikl,k  )
          ENDDO
          DO k=    2,mzp
            B___AT(ikl,k  ) =  1.0000 - A___AT(ikl,k  ) - C___AT(ikl,k  )
          ENDDO
            B___AT(ikl,1  ) =  1.0000 - A___AT(ikl,1  )
            C___AT(ikl,1  ) =  0.0000

            D___AT(ikl,1  ) =                                               &
     &      A___AT(ikl,1  ) *  a_b_AT *(var_AT(ikl,1  ) - var_AT(ikl,2  ))

             k=      mzp
            D___AT(ikl,k  ) =                                               &
     &      A___AT(ikl,k  ) * (a_b_AT * var_AT(ikl,k  ) - var_AT(ikl,k+1)   &
     &                                                  / betaAT)           &
     &    - C___AT(ikl,k  ) *  a_b_AT *(var_AT(ikl,k-1) - var_AT(ikl,k  ))  &
     &    +                                               var_AT(ikl,k  )
            A___AT(ikl,k  ) =  0.0000
          DO k=    2,mzp-1
            D___AT(ikl,k  ) =                                               &
     &      A___AT(ikl,k  ) *  a_b_AT *(var_AT(ikl,k  ) - var_AT(ikl,k+1))  &
     &    - C___AT(ikl,k  ) *  a_b_AT *(var_AT(ikl,k-1) - var_AT(ikl,k  ))  &
     &    +                                               var_AT(ikl,k  )
          ENDDO
            D___AT(ikl,1  ) =                                               &
     &      A___AT(ikl,1  ) *  a_b_AT *(var_AT(ikl,1  ) - var_AT(ikl,2  ))  &
     &    +                                               var_AT(ikl,1  )



! Tri-Diagonal Matrix Resolution (Gaussian Elimination / Thomas Algorithm)
! ==============================

! Forward  Sweep
! --------------

          P___AT(1) =  B___AT(ikl,1)
          Q___AT(1) = -A___AT(ikl,1) / P___AT(1)
        DO k=2,mzp
          P___AT(k) =  C___AT(ikl,k) * Q___AT(k-1)   + B___AT(ikl,k)
          Q___AT(k) = -A___AT(ikl,k) / P___AT(k)
        ENDDO
          X___AT(1) =  D___AT(ikl,1) / P___AT(1)
        DO k=2,mzp
          X___AT(k) = (D___AT(ikl,k) - C___AT(ikl,k) * X___AT(k-1)) / P___AT(k)
        END DO


! Backward Sweep
! --------------

        DO k=  mzp-1,1,-1
          X___AT(k) =  Q___AT(k)     * X___AT(k+1)   + X___AT(k)
        END DO



! Elimination of too small values
! -------------------------------

        IF (SchmAT.EQ.'w'.OR.SchmAT.EQ.'i'.OR.SchmAT.EQ.'C'.OR.        &
     &      SchmAT.EQ.'s'.OR.SchmAT.EQ.'r')                         THEN
          DO k=  mzp  ,1,-1
            IF (X___AT(k).NE.0.0000 .AND. X___AT(k).LT.1.e-18)         &
     &          X___AT(k) =  0.0000
          END DO
        END IF




! Tendency
! ========

!  Wind Speed
!  ----------

        IF      (schmAT.EQ.'a')                                     THEN
          DO k=1,mzp
            dua_AT(ikl,k) =(X___AT(k) - ua__DY(ikl,k))  / dt__AT
          ENDDO
        ELSE IF (schmAT.EQ.'b')                                     THEN
          DO k=1,mzp
            dva_AT(ikl,k) =(X___AT(k) - va__DY(ikl,k))  / dt__AT
          ENDDO
        ELSE IF (schmAT.EQ.'u')                                     THEN
          DO k=1,mzp
            dua_AT(ikl,k) =(X___AT(k) - ua__DY(ikl,k))  / dt__AT
          ENDDO
        ELSE IF (schmAT.EQ.'v')                                     THEN
          DO k=1,mzp
            dva_AT(ikl,k) =(X___AT(k) - va__DY(ikl,k))  / dt__AT
          ENDDO


!  Temperature
!  -----------

        ELSE IF (schmAT.EQ.'T')                                     THEN
          DO k=1,mzp
            dpktAT(ikl,k) =(X___AT(k) - pkt_DY(ikl,k))  / dt__AT
          ENDDO


!  Water Species
!  -------------

        ELSE IF (schmAT.EQ.'q')                                     THEN
          DO k=1,mzp
            dqv_AT(ikl,k) =(X___AT(k) - qv__DY(ikl,k))  / dt__AT
          ENDDO

        ELSE IF (schmAT.EQ.'w')                                     THEN
          DO k=1,mzp
            dqw_AT(ikl,k) =(X___AT(k) - qw__CM(ikl,k))  / dt__AT
          ENDDO

        ELSE IF (schmAT.EQ.'i')                                     THEN
          DO k=1,mzp
            dqi_AT(ikl,k) =(X___AT(k) - qi__CM(ikl,k))  / dt__AT
          ENDDO

        ELSE IF (schmAT.EQ.'C')                                     THEN
          DO k=1,mzp
            dCi_AT(ikl,k) =(X___AT(k) - CCNiCM(ikl,k))  / dt__AT
          ENDDO

        ELSE IF (schmAT.EQ.'s')                                     THEN
          DO k=1,mzp
            dqs_AT(ikl,k) =(X___AT(k) - qs__CM(ikl,k))  / dt__AT
          ENDDO

        ELSE IF (schmAT.EQ.'r')                                     THEN
          DO k=1,mzp
            dqr_AT(ikl,k) =(X___AT(k) - qr__CM(ikl,k))  / dt__AT
          ENDDO

        ELSE
            write(6,*) 
            write(6,*) ' schmAT = ',schmAT
            write(6,*) ' Inadequate..., NOR for Ekman Spiral, NOR for momentum, NOR for scalars'
            write(6,*) '                                                    '

!                **************
            call MAR________BUG
!                **************
            STOP
        END IF


!     ++++++
      END DO
!     ++++++




      end subroutine Atm_AT
