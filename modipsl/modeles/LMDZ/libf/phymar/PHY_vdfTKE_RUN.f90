      subroutine PHY_vdfTKE_RUN

!------------------------------------------------------------------------------+
!                                                         Mon 17-Jun-2013  MAR |
!     SubRoutine PHY_vdfTKE_RUN  includes Vertical Diffusion of TKE            |
!                                                           and epsilon        |
!                                                                              |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 19-Mar-2013      |
!           Last Modification by H. Gallee,               Mon 17-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+
!                                                                              |
!    INPUT: Kzm_AT           Vertical Turbulent Coeffic.(momentum) [m2/s2]     |
!    ^^^^^^                                                                    |
!                                                                              |
!    INPUT / OUTPUT: The Vertical Turbulent Fluxes are included for:           |
!    ^^^^^^^^^^^^^^                                                            |
!         a) Turbulent Kinetic Energy             TKE_AT(_xyz)     [m2/s2]     |
!         b) Turbulent Kinetic Energy Dissipation eps_AT(_xyz)     [m2/s3]     |
!                                                                              |
!   #OPTIONS: #De: Dirichlet Type Top Boundary Condit. for TKE_AT (TKE    )    |
!   #^^^^^^^^                                            & eps_AT (epsilon)    |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_AT_grd
      use Mod_PHY_AT_kkl
      use Mod_PHY_DY_kkl



! Local  Variables
! ================

      use Mod_vdfTKE_RUN


      IMPLICIT NONE


      real(kind=real8)                             ::  S3DSBC
      real(kind=real8)                             ::  sige2k
      real(kind=real8)                             ::  psa_sq
! #De real(kind=real8)                             ::  TKEtop = 0.

      integer                                      ::  i     ,j     ,k
      integer                                      ::  k1           ,ikl




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
! ALLOCATION
! ==========

      IF (it_RUN.EQ.1 .OR. FlagDALLOC)               THEN !
          allocate ( S3D__1(mzp-1) )
          allocate ( S3D__2(mzp-1) )
          allocate ( S3D__3(mzp-1) )
          allocate ( S3D__4(mzp-1) )
          allocate ( S3D__5(mzp-1) )
          allocate ( S3D__6(mzp-1) )
          allocate ( S3D__7(mzp-1) )
      END IF                                              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! INITIALIZATION
! ==============

      sige2k = sige / sigk



!     =====================
      DO      ikl = 1,kcolp
!     =====================

        i = ii__AP(ikl)
        j = jj__AP(ikl)

        psa_sq          = psa_DY(ikl)         *psa_DY(ikl)    




! Vertical Diffusion of Turbulent Kinetic Energy
! ==============================================


! Tridiagonal Matrix Coefficients - TKE_AT 
! ----------------------------------------

           k=mzp-1
          S3DSBC        = -GravF2*0.5*(Kzm_AT(ikl,k)+Kzm_AT(ikl,k+1))  &!  SBC: TKE and epsilon, Atm SBL
     &             * sigk *betaAT    * roamDY(ikl,k)*roa_DY(ikl,k  )   &!
     &                   /(psa_sq    * dsigmi(k)    *dsigma(    k+1))
          S3D__1(mzp-1) =  S3DSBC                                       !  SBC: TKE and epsilon, Atm SBL

        DO k=mzp-2,1,-1
          S3D__1(k)      =-GravF2*0.5*(Kzm_AT(ikl,k)+Kzm_AT(ikl,k+1))  &
     &             * sigk *betaAT    * roamDY(ikl,k)*roa_DY(ikl,k  )   &
     &                   /(psa_sq    * dsigmi(k)    *dsigma(    k+1))

          S3D__3(k+1)    = S3D__1(k) * dsigmi(    k)/dsigmi(    k+1)   &
     &                               / roa_DY(ikl,k)*roa_DY(ikl,k+1) 
        END DO

          S3D__3(1)      = 0.0                                          !  UBC: Von Neuman     , Atm Top
        DO k=   1,mzp-1
          S3D__1(k) =      S3D__1(k) * dt__AT
          S3D__3(k) =      S3D__3(k) * dt__AT
          S3D__2(k) = 1.0 -S3D__3(k) -S3D__1(k) 
        END DO


! Second Member of the Tridiagonal System - TKE_AT
! ------------------------------------------------

          S3D__4(1) =                                                   &
     &    S3D__1(1) *a_b_AT*(TKE_AT(ikl,1)-TKE_AT(ikl,k1p(1)))
! #De     S3D__1(1) = 0.0
! #De     S3D__2(1) = 1.0
! #De     S3D__4(1) = TKEtop

        DO k=k1p(1),mzp-2
          S3D__4(k) =                                                   &
     &    S3D__1(k) *a_b_AT*(TKE_AT(ikl,    k )-TKE_AT(ikl,k1p(1)))     &
     &   -S3D__3(k) *a_b_AT*(TKE_AT(ikl,k1m(k))-TKE_AT(ikl,    k ))
        END DO

           k=       mzp-1
          S3D__4(k) =                                                   &
     &    S3D__1(k)*(a_b_AT* TKE_AT(ikl,    k )-TKE_AT(ikl,k1p(1))      &
     &                                         /betaAT            )     &
     &   -S3D__3(k) *a_b_AT*(TKE_AT(ikl,k1m(k))-TKE_AT(ikl,    k ))

          S3D__1(k) = 0.000

!         S3D__4(mzp-1)=-(alphAT* TKE_AT(ikl,mzp-1)-TKE_AT(ikl,mzp  ))  &
!    &          * GravF2 *0.5000*(Kzm_AT(ikl,mzp-1)+Kzm_AT(ikl,mzp-2))  &
!    &                          * roamDY(ikl,mzp-1)*roamDY(ikl,mzp-1)   &
!    &                 / (psa_sq* dsigmi(    mzp-1)*dsigma(    mzp-1))  &
!    &   -S3D__3(mzp-1)*  a_b_AT*(TKE_AT(ikl,mzp-2)-TKE_AT(ikl,mzp-1))


! Tridiagonal Matrix Inversion - TKE_AT 
! -------------------------------------

             k1= 1
! #De        k1= 2
        DO k=k1,mzp-1
          S3D__4(k) = S3D__4(k) +  TKE_AT(ikl,k)
        END DO

! Forward  Sweep
! ~~~~~~~~~~~~~~
          S3D__5(1) = S3D__2(1)
          S3D__6(1) =-S3D__1(1)   /S3D__5(1)
        DO   k=k1p(1),mzp-1
          S3D__5(k) = S3D__3(k)   *S3D__6(k-1)+S3D__2(k)
          S3D__6(k) =-S3D__1(k)   /S3D__5(k)
        END DO
          S3D__7(1) = S3D__4(1)   /S3D__5(1)
        DO   k=k1p(1),mzp-1
          S3D__7(k) =(S3D__4(k)   -S3D__3(k)                           &
     &               *S3D__7(k-1))/S3D__5(k)
        END DO

! Backward Sweep
! ~~~~~~~~~~~~~~
        DO k=k1m(mzp-1),1,-1
          S3D__7(k) = S3D__6(k)   *S3D__7(k+1)+S3D__7(k)
        END DO


        DO k=1,mzp-1
          TrT_AT(ikl,k) =  TrT_AT(ikl,k)                               &
     &              +(S3D__7(k)   -TKE_AT(ikl,k))  /dt__AT
          TKE_AT(ikl,k) =                      S3D__7(k)
        END DO


! Vertical Diffusion of Dissipation
! =================================


! Update Tridiagonal Matrix Coefficients - eps_AT 
! -----------------------------------------------

          S3D__1(mzp-1) =   S3DSBC    * dt__AT                          !  SBC: TKE and epsilon, Atm SBL
        DO k=1,mzp-1
          S3D__1(k) =       S3D__1(k) * sige2k
          S3D__3(k) =       S3D__3(k) * sige2k
          S3D__2(k) = 1.0 - S3D__3(k) - S3D__1(k) 
        END DO


! Second Member of the Tridiagonal System - eps_AT
! ------------------------------------------------

          S3D__4(1) =                                                   &
     &    S3D__1(1) *a_b_AT*(eps_AT(ikl,1)-eps_AT(ikl,k1p(1)))
! #De     S3D__1(1) = 0.0
! #De     S3D__2(1) = 1.0
! #De     S3D__4(1) = eps_DI(i,j)

        DO k=k1p(1),mzp-2
          S3D__4(k) =                                                   &
     &    S3D__1(k) *a_b_AT*(eps_AT(ikl,k)-eps_AT(ikl,k1p(k)))          &
     &   -S3D__3(k) *a_b_AT*(eps_AT(ikl,k1m(k))-eps_AT(ikl,k))
        END DO

           k=       mzp-1
          S3D__4(k) =                                                   &
     &    S3D__1(k)*(a_b_AT* eps_AT(ikl,    k )-eps_AT(ikl,k1p(1))      &
     &                                         /betaAT            )     &
     &   -S3D__3(k) *a_b_AT*(eps_AT(ikl,k1m(k))-eps_AT(ikl,    k ))

          S3D__1(k) = 0.000

!         S3D__4(mzp-1)=-(alphAT* eps_AT(ikl,mzp-1)-eps_AT(ikl,mzp))    &
!    &          * GravF2* 0.5000*(Kzm_AT(ikl,mzp-1)+Kzm_AT(ikl,mzp-2))  &
!    &                          * roamDY(ikl,mzp-1)*roamDY(ikl,mzp-1)   &
!    &                 / (psa_sq* dsigmi(    mzp-1)*dsigma(    mzp-1))  &
!    &   -S3D__3(mzp-1)*  a_b_AT*(eps_AT(ikl,mzp-2)-eps_AT(ikl,mzp-1))


! Tridiagonal Matrix Inversion - eps_AT 
! -------------------------------------

             k1= 1
! #De        k1= 2
        DO k=k1,mzp-1
          S3D__4(k)    = S3D__4(k) + eps_AT(ikl,k)
        END DO

! Forward  Sweep
! ~~~~~~~~~~~~~~
          S3D__5(1) = S3D__2(1)
          S3D__6(1) =-S3D__1(1)   /S3D__5(1)
        DO   k=k1p(1),mzp-1
          S3D__5(k) = S3D__3(k)   *S3D__6(k-1)+S3D__2(k)
          S3D__6(k) =-S3D__1(k)   /S3D__5(k)
        END DO
          S3D__7(1) = S3D__4(1)   /S3D__5(1)
        DO   k=k1p(1),mzp-1
          S3D__7(k) =(S3D__4(k)   -S3D__3(k)                            &
     &               *S3D__7(k-1))/S3D__5(k)
        END DO

! Backward Sweep
! ~~~~~~~~~~~~~~
        DO k=k1m(mzp-1),1,-1
          S3D__7(k) = S3D__6(k)   *S3D__7(k+1)+S3D__7(k)
        END DO

        DO k=1,mzp-1
          eps_AT(ikl,k) =                      S3D__7(k)
        END DO




!     =====================
      ENDDO ! ikl = 1,kcolp
!     =====================




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
! DE-ALLOCATION
! =============

      IF (FlagDALLOC)                                THEN !
          deallocate ( S3D__1 )
          deallocate ( S3D__2 )
          deallocate ( S3D__3 )
          deallocate ( S3D__4 )
          deallocate ( S3D__5 )
          deallocate ( S3D__6 )
          deallocate ( S3D__7 )
      ENDIF
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      return
      end subroutine PHY_vdfTKE_RUN
