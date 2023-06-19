      subroutine PHY________INI

!------------------------------------------------------------------------------+
!                                                         Sun 30-Jun-2013  MAR |
!   MAR          PHY________INI                                                |
!     subroutine PHY________INI intializes MAR PHYsical parameterizations      |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,               Sun 30-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_RT_grd
      use Mod_PHY____kkl


      IMPLICIT NONE


      logical                               ::  search_argexp = .FALSE.         !  Lapack used to compute MAX/MIN exponential arguments
      integer                               ::       i,     j,   ikl            !




!=============================================================================================!
!                                                                                             !
!     include 'MARphy.inc'                                        ! MARthusalem constants
!                                                                                             !
!=============================================================================================!


!=============================================================================================!
!                                                                                             !
! Modification   of Mod_PHY____dat (needed if constants slighly differ  in the HOST model)    !
! ================================ (Here  the chosen HOST model is MAR)                       !
!                                                                                             !
! Initialization of Mod_PHY____dat from MARphy.inc (MARphy.inc contains MARthusalem constants)!
! ------------------------------------------------                                            !
!                                                                                             !
!     zer0   =       0.0               !                                                 [-]  !
!     half   =       0.5               !                                                 [-]  !
!     un_1   =       1.0               !                                                 [-]  !
      piNmbr = acos(-1.0)              !                                                 [-]  !
!     eps6   = epsi                    !                                                 [-]  !
!     epsn   = eps9                    !                                                 [-]  !
!     A_MolV = 1.35e-5                 !   Air        Viscosity            1.35e-5    [m2/s]  !
!     rhoIce = ro_Ice                  !   Ice        Specific Mass      920.e0      [kg/m3]  !
!     BSnoRo = blsno                   !   Blown Snow Specific Mass        2.55e+2   [kg/m3]  !
!     LhfH2O = Lf_H2O                  !
!     LhvH2O = Lv_H2O                  !   Latent Heat Vaporisation, Water 2.5008e+6  [J/kg]  !
!     LhsH2O = Ls_H2O                  !
!     CpdAir = Cp                      !   Air   Heat Capacity  (p=C)   1004.708845 [J/kg/K]  !
!     R_DAir = RDryAi                  !   Dry Air Perfect Gas Law C     287.05967  [J/kg/K]  !
!     RCp    = cap                     !   RDryAi         / Cp                           [-]  !
!     p0_kap = pcap                    !
!     hC_Wat = C__Wat                  !   H2O   Heat Capacity          4186.00e0   [J/kg/K]  !
!     rhoWat = ro_Wat
!     Tf_Sno = TfSnow
!     Tf_Sea = tfrwat
!     StefBo = stefan                  !   Stefan-Bolstzman Constant       5.67e-8 [W/m2/K4]  !
!     Grav_F = gravit                  !   Gravity    Acceleration         9.81e0     [m/s2]  !
!     vonKrm = vonkar                  !   von Karman Constant             0.40e0        [-]  !
!     A_Stab = A_Turb                  !                                   5.8           [-]  !
!     AhStab = AhTurb
!     AsStab = AsTurb
!     r_Stab = r_Turb
!                                                                                             !
!=============================================================================================!



! Initialization of Mod_PHY____dat (auxiliary Constants)                                      !
! ------------------------------------------------------                                      !

      Grav_I = 1.      /  Grav_F       !                                              [s2/m]  !
      GravF2 = Grav_F  *  Grav_F       !                                             [m2/s4]  !
      RCp    = R_DAir  /  CpdAir       !                        Case Sensitive           [-]  !
      Lv_CPd = LhvH2O  /  CpdAir
      Ls_CPd = LhsH2O  /  CpdAir
      Lc_CPd = LhfH2O  /  CpdAir
      vonKrm = 0.35                    !   von Karman Constant, Case Sensitive                !

      IF (search_argexp)                                            THEN

!          ***************
      call PHY_CPU_NumPrec
!          ***************

      ELSE

      ea_MIN =-80.
      ea_MAX = 80.

      END IF




! Initialization of Mod_PHY____grd
! ================================


! Correspondance entre la grille 2D horizontale dynamique et
! -------------------- la grille 2D horizontale physique utile mxp,myp,mzp
!                      ---------------------------------------------------
       mxp   = mxpp-ixp1+1          ! 
       myp   = mypp-jyp1+1          ! 
       !kcolp = mxp * myp            ! Déja calculé en amont dans physiq.F90 Martin
       mzp   = mzpp-1
       PRINT*, 'mxpp=',mxpp
       PRINT*, 'mxp=',mxp
       PRINT*, 'ixp1=',ixp1
       PRINT*, 'mypp=',mypp
       PRINT*, 'jyp1=',jyp1
       PRINT*, 'myp=',myp


! Horizontal Cartesian Grid
! -------------------------

           write(6,*) ' '
           write(6,*) 'i_x0  , j_y0                = '                  &
     &                ,i_x0  , j_y0


                                    !        dxHOST is Model Grid Size
         dx2inv  = 0.5/dxHOST       ! 1 / (2 dxHOST)
         dy2inv  = 0.5/dxHOST       ! 1 / (2 dxHOST)




! ALLOCATION
! ==========

!                ****************
          CALL   PHY________ALLOC
!                ****************




! Initialization of Mod_PHY____grd
! ================================

! Initialization of the Correspondance between 2-D cartesian and vector Grid
! --------------------------------------------------------------------------

! Adapted for MAR/LMDZ coupling:

!     ii__AP(1)=ixp1
!     jj__AP(1)=jyp1
!     PRINT*,'jyp1=',jyp1
!     PRINT*,'jj__AP(1)=',jj__AP(1)
!     ii__AP(kcolp)=ixp1
!     jj__AP(kcolp)=mypp
!
!     DO i=ixp1,mxpp
!     DO j=jyp1+1,mypp-1
!
!               ikl    = (j-(jyp1+1)) *mxpp +1 + (i-ixp1+1) ! Tout est décalé de 1 à cause du point isolé au pole dans la grille physique LMD
!        !       ikl    = (j-jyp1) *mxpp + i-ixp1+1
!        !       ikl    = (j-jyp1) *(mxp-1) + i-ixp1+1
!        ii__AP(ikl)   =                 i
!        jj__AP(ikl)   =  j
!        ikl_AP(i,j)   =  ikl
!        PRINT*,'ii__AP(',ikl,')=',ii__AP(ikl)
!        PRINT*,'jj__AP(',ikl,')=',jj__AP(ikl)
!     ENDDO
!     ENDDO

! Modification Gilles Delaygue 2014/07/14 !
      ikl=1
      ii__AP(ikl)=ixp1
      jj__AP(ikl)=jyp1
      ikl_AP(:,jyp1) = ikl

      DO j=jyp1+1,mypp-1
      DO i=ixp1,mxpp
         ikl=ikl+1
         ii__AP(ikl)   =  i
         jj__AP(ikl)   =  j
         ikl_AP(i,j)   =  ikl
      ENDDO
      ENDDO

      ikl=ikl+1
      ii__AP(ikl)=ixp1
      jj__AP(ikl)=mypp
      ikl_AP(:,mypp) = ikl





      PRINT*,'Control dans PHY_INI:'
      PRINT*,'ii__AP(1)=',ii__AP(1)
      PRINT*,'ii__AP(kcolp)=',ii__AP(kcolp)
      PRINT*,'jj__AP(1)=',jj__AP(1)
      PRINT*,'jj__AP(kcolp)=',jj__AP(kcolp)
        
! Martin control tout sauf les poles:
     WRITE(6,600)(ii__AP(ikl),ikl=2,kcolp-1)
     600 FORMAT (48i2) 
     WRITE(6,601)(jj__AP(ikl),ikl=2,kcolp-1)
     601 FORMAT (48i2) 

!      DO i=ixp1,mxpp
!      DO j=jyp1,mypp
!
!                ikl    = (j-jyp1) *mxp + i-ixp1+1
!         ii__AP(ikl)   =                 i
!         jj__AP(ikl)   =  j
!         ikl_AP(i,j)   =  ikl
!
!      ENDDO
!      ENDDO

         ikl0          =  ikl_AP(i_x0,j_y0)




! Allocation of radiative transfert Variables
! ===========================================


! Initialization of Mod_PHY_RT_grd
! --------------------------------

         naero         =  6


!                            ****************
                       CALL  PHY_Atm_RT_ALLOC
!                            ****************




! Allocation of microphysical       Variables
! ===========================================

!                            ****************
                       CALL  PHY_Atm_CM_ALLOC
!                            ****************




! Allocation of Turbulence          Variables
! ===========================================

!                            ****************
                       CALL  PHY_Atm_AT_ALLOC
!                            ****************




! Allocation of Surface Variables
! ===============================

!                            ****************
                       CALL  PHY_SISVAT_ALLOC
!                            ****************




! OUTPUT
! ======

      OPEN(unit=4,status='unknown',file='PHY___________.OUT')
      REWIND    4



      end subroutine PHY________INI
