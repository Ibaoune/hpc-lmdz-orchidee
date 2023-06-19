      module Mod_PHY_AT_kkl

!--------------------------------------------------------------------------+
!                                                     Sun 12-May-2013  MAR |
!     module Mod_PHY_AT_kkl contains the main (prognostic) variables of    |
!                Turbulent Vertical Diffusion Scheme                       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,           Sun 12-May-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



! E-e and K-l models parameters
! -----------------------------

      real(kind=real8), SAVE   ::   TKEmin       = 0.0001   ! Minimum SBL    turbulent kinetic   energy

      real(kind=real8), SAVE   ::      cmub      = 0.0900   ! Ee Model Parameter (Bintanja , 2000, BLM (95),       mid p.355)
      real(kind=real8), SAVE   ::   sqrcmub      = 3.333    !
      real(kind=real8), SAVE   ::     c1epb      = 1.46     !
      real(kind=real8), SAVE   ::     c2epb      = 1.92     !
      real(kind=real8), SAVE   ::     sigeb      = 0.862    !
      real(kind=real8), SAVE   ::     sigkb      = 1.000    !

      real(kind=real8), SAVE   ::      cmud      = 0.0330   ! Ee Model Parameter (Duynkerke, 1988, JAS (45), (19), top p.868)
      real(kind=real8), SAVE   ::   sqrcmud      = 5.500    !                    (c_mu)^1/2=(0.033)^1/2=5.50           p.869
      real(kind=real8), SAVE   ::     c1epd      = 1.46     !                                                          p.868
      real(kind=real8), SAVE   ::     c2epd      = 1.83     !                                                          p.868
      real(kind=real8), SAVE   ::     siged      = 0.420    !
      real(kind=real8), SAVE   ::     sigkd      = 1.000    !

      real(kind=real8), SAVE   ::      cmuk      = 0.0900   ! Ee Model Parameter (Kitada   , 1987, BLM (41),       top p.220)
      real(kind=real8), SAVE   ::   sqrcmuk      = 3.333    !                    (c_mu)^1/2=(0.090)^1/2=3.333
      real(kind=real8), SAVE   ::     c1epk      = 1.44     !                                                      top p.220
      real(kind=real8), SAVE   ::     c2epk      = 1.92     !
      real(kind=real8), SAVE   ::     sigek      = 0.769    !
      real(kind=real8), SAVE   ::     sigkk      = 1.000    !

      real(kind=real8), SAVE   ::   sqrcmut      = 4.0000   ! Kl Model Parameter (Schayes & Thunis, 1990, Contrib. 60 Inst.Astr.Geoph. p.8)
      real(kind=real8), SAVE   ::     siget      = 0.420    !
      real(kind=real8), SAVE   ::     sigkt      = 1.200    !

      real(kind=real8), SAVE   ::    betahr      = 2.0000   ! Ee Model Parameter (Huang and Raman,  1991, BLM (55), p.386 and (A22) p.405)

      real(kind=real8), SAVE   ::      cmu                  !
      real(kind=real8), SAVE   ::   sqrcmu                  !
      real(kind=real8), SAVE   ::     c1ep                  !
      real(kind=real8), SAVE   ::     c2ep                  !
      real(kind=real8), SAVE   ::     sige                  !
      real(kind=real8), SAVE   ::     sigk                  !
      real(kind=real8), SAVE   ::   vK_inv                  ! Inverse      of  Von-Karman Constant

! #KA real(kind=real8), SAVE   ::   zz__KA       = 5.0000   ! Height below which use a vertical weighted box filter of TKE, e  [m.agl]
! #KA integer, SAVE            ::   mz__KA                  ! Level  below which use a vertical weighted box filter of TKE, e  
! #KA logical            ::   log_KA       = .true.   ! Swich  deciding use of a vertical weighted box filter of TKE, e  

      

! PHY_AT INPUT        Variables
! -----------------------------

      integer, SAVE         ,ALLOCATABLE ,dimension(:)    ::  ii__AT  ! WORK   point   i Coordinate
      integer, SAVE         ,ALLOCATABLE ,dimension(:)    ::  jj__AT  ! WORK   point   j Coordinate
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:)  ::  ikl_AT  ! WORK   point vec Coordinate

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  var_AT  ! Dummy    to Diffuse                                        [x]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  Ac0_AT  ! Tridiagonal Matrix Coefficient A: Common Factor        [m2/s3]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  Cc0_AT  ! Tridiagonal Matrix Coefficient C: Common Factor        [m2/s3]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  Ac__AT  ! Tridiagonal Matrix Coefficient A: Common Factor (t)     [s/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  Cc__AT  ! Tridiagonal Matrix Coefficient C: Common Factor (t)     [s/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  Kz0_AT  ! Vertical  Turbulent Diffusion Coefficient (MINIMUM)     [m2/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  Kz__AT  ! Vertical  Turbulent Diffusion Coefficient               [m2/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  Kzm_AT  ! Vertical  Turbulent Diffusion Coefficient (Momentum)    [m2/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  Kzh_AT  ! Vertical  Turbulent Diffusion Coefficient (Scalars)     [m2/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  Kzh0AT  ! Vertical  Turbulent Diffusion Coefficient (Scalars)     [m2/s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  A___AT  ! Tridiagonal Matrix Coefficient A                           [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  B___AT  ! Tridiagonal Matrix Coefficient B                           [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  C___AT  ! Tridiagonal Matrix Coefficient C                           [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  D___AT  ! Independant Term               D                           [x]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  P___AT  ! Auxiliary   Term               P                           [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  Q___AT  ! Auxiliary   Term               Q                           [-]
!      real*16         ,ALLOCATABLE ,dimension(:)    ::  X___AT  ! Auxiliary   Unknown            X                           [x]
      double precision,ALLOCATABLE ,dimension(:)    ::  X___AT  ! Auxiliary   Unknown            X                           [x]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  LMO_AT  ! Monin-Obukhov     Length         (Grid Cell Average)       [m]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:)    ::  zi__AT  ! Inversion         Height         (Grid Cell Average)[m a.g.l.]



! PHY_AT INPUT/OUTPUT Variables
! -----------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  TKE_AT  ! Turbulent Kinetic Energy                               [m2/s2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  TrT_AT  ! Turbulent Kinetic Energy Transport                     [m2/s2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  eps_AT  ! Turbulent Kinetic Energy Dissipation                   [m2/s3]



! PHY_AT OUTPUT       Variables
! -----------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  dua_AT  ! Wind Speed (x-direc.) Tendency                          [m/s2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  dva_AT  ! Wind Speed (y-direc.) Tendency                          [m/s2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  dpktAT  ! Potential Temperature Tendency, divided by p0**(R/Cp)   [x/s] 
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  dqv_AT  ! Specific  Humidity    Tendency                      [kg/kg/s] 
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  dqw_AT  ! Cloud Droplets Concen.Tendency                      [kg/kg/s] 
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  dqi_AT  ! Cloud Crystals Concen.Tendency                      [kg/kg/s] 
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  dCi_AT  ! CCNi           Concen.Tendency                          [1/s] 
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  dqs_AT  ! Snow  Particls Concen.Tendency                      [kg/kg/s] 
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)  ::  dqr_AT  ! Rain  Drops    Concen.Tendency                      [kg/kg/s] 



      end module Mod_PHY_AT_kkl
