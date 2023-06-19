     module Mod_PHY____dat

!------------------------------------------------------------------------------+
!                                                         Sat  8-Jun-2013  MAR |
!    module Mod_PHY____dat  contains specific constants for                    |
!           MAR PHYsical    Parameterizations                                  |
!                                                                              |
!                     These constants may be slightly different in HOST Model  |
!                                 and may be modified in  PHY________INI       |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,               Sat  8-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+



! General Variables
! =================

      use Mod_Real


      IMPLICIT NONE


! Physical Constants (may be changed in the interface with the Host Model)
! ==================

      real(kind=real8), SAVE  ::  zer0   =    0.00         !                                !    0.00
      real(kind=real8), SAVE  ::  half   =    0.50         !                                !    0.50
      real(kind=real8), SAVE  ::  un_1   =    1.00         !                                !    1.00
      real(kind=real8), SAVE  ::  piNmbr =    3.1416       ! pi                             !    3.1416
      real(kind=real8), SAVE  ::  Dg2Rad = 1745.33d-5      ! pi / 180                       !           
      real(kind=real8), SAVE  ::  epsq   =    1.e-6        ! Arbirary Small Value, H2O      !    1.e-6
      real(kind=real8), SAVE  ::  eps1   =    1.e-1        ! Arbirary Small Value           !    1.e-1
      real(kind=real8), SAVE  ::  eps6   =    1.e-6        ! Arbirary Small Value           !    1.e-6
      real(kind=real8), SAVE  ::  epsn   =    1.e-9        ! Arbirary Small Value           !    1.e-9
      real(kind=real8), SAVE  ::  epsp   =    1.e-12       ! Arbirary Small Value           !    1.e-12
      real(kind=real8), SAVE  ::  R_1by3 =    0.333333     ! 1 / 3                          !
      real(kind=real8), SAVE  ::  R_5by3 =    1.666666     ! 5 / 3                          !
      real(kind=real8), SAVE  ::  R_1000 =    1.e+3        !                                !    1.e+3
      real(kind=real8), SAVE  ::  ea_MAX =   50.           ! MAX allowed exponential Argum. ! computed by HOST
      real(kind=real8), SAVE  ::  ea_MIN =  -50.           ! MIN allowed exponential Argum. ! computed by HOST
      real(kind=real8), SAVE  ::  A_MolV =    1.35d-5      ! Air Molecular Viscosity        !    1.35d-5 m2/s 
      real(kind=real8), SAVE  ::  rhoIce =  920.00d+0      ! Ice         Specific Mass      !  920.00d+0 kg/m3
      real(kind=real8), SAVE  ::  BSnoRo =  255.00d+0      ! Blown Snow  Specific Mass      !  255.00d+0 kg/m3
      real(kind=real8), SAVE  ::  LhvH2O = 2508.00d+3      ! Latent Heat of Vapor. of Rain  ! 2500.00d+3 J/kg
      real(kind=real8), SAVE  ::  LhfH2O =    3.34d+5      ! Latent Heat of Fusion of Snow  !    3.34d+5 J/kg
      real(kind=real8), SAVE  ::  LhsH2O = 2833.60d+3      ! Latent Heat of Sublim.of Snow  ! 2833.60d+3 J/kg
      real(kind=real8), SAVE  ::  CpdAir = 1004.708845     ! dry air specific heat at cst p ! 1004.00    J/kg/K
      real(kind=real8), SAVE  ::  R_DAir =  287.05967      ! dry air perfect gas law  cst   !  287.      J/kg/K
      real(kind=real8), SAVE  ::  RCp    =    0.285857     ! R / Cp                         ! 287./1004. -
      real(kind=real8), SAVE  ::  p0_kap =    3.730037     !                                ! 100 kPa ** (R/Cp)
      real(kind=real8), SAVE  ::  Lv_CPd                   ! LhvH2O      /  CpdAir          !
      real(kind=real8), SAVE  ::  Ls_CPd                   ! LhsH2O      /  CpdAir          !
      real(kind=real8), SAVE  ::  Lc_CPd                   ! LhfH2O      /  CpdAir          !
      real(kind=real8), SAVE  ::  hC_Wat = 4186.00d+0      ! Water Heat Capacity            ! 4186.00d+0 J/kg/K
      real(kind=real8), SAVE  ::  rhoWat = 1000.00d+0      ! Water Specific Mass            ! 1000.00d+0 kg/m3
      real(kind=real8), SAVE  ::  Tf_Sno =  273.16         ! Snow  Melting  Point           !  273.16    K
      real(kind=real8), SAVE  ::  Tf_Sea =  271.2          ! Sea   Melting  Point           !  271.2     K
      real(kind=real8), SAVE  ::  StefBo =    5.67d-8      ! Stefan-Boltzman Constant       !    5.67d-8 W/m2/K4
      real(kind=real8), SAVE  ::  Grav_F =    9.81         !     Gravitational  Force       !    9.81    m/s2
!     real(kind=real8), SAVE  ::  Grav_I =    0.101937     ! 1 /(Gravitational  Force)      ! 1 /9.81    s2/m
      real(kind=real8), SAVE  ::  Grav_I                   ! 1 /(Gravitational  Force)      ! 1 /9.81    s2/m
      real(kind=real8), SAVE  ::  GravF2                   !    (Gravitational  Force) ** 2 !    9.81    m2/s4
      real(kind=real8), SAVE  ::  vonKrm =    0.4          ! von Karman Constant            !    0.4
      real(kind=real8), SAVE  ::  A_Stab =    5.8          ! Stability  Coefficient Moment  !    5.8
      real(kind=real8), SAVE  ::  AhStab =    5.4          ! Stability  Coefficient Heat    !    5.4
      real(kind=real8), SAVE  ::  AsStab =    4.0          ! Stability  Coefficient Blown * !    4.0
      real(kind=real8), SAVE  ::  r_Stab =    3.0          ! Turbul.Diffusivit.Ratio K*/Km  !        

      real(kind=real8), SAVE  ::  EarthR = 6371.229e3      ! Earth   Radius                 !            m
      real(kind=real8), SAVE  ::  DirAxX                   ! x-Axis  Direction              !   90       degrees is the most natural choice
      real(kind=real8), SAVE  ::  sh_MAX                   ! Highest Domain  Grid  Point    !        
      real(kind=real8), SAVE  ::  dzaMIN                   ! Thinest Atmosph.Layer Thickness!        

      character(len=3), dimension(0:12)  ::  LabMon  !
      data LabMon /'---','Jan','Feb','Mar','Apr','May','Jun'           &
     &                  ,'Jul','Aug','Sep','Oct','Nov','Dec'/

      integer, SAVE         , dimension(0:12)  ::  njYear  ! Nb of Days since Begin of Year !
      data njYear / 0   ,   0 ,  31 ,  59 ,  90 , 120 , 151            &
     &                  , 181 , 212 , 243 , 273 , 304 , 334/

      integer, SAVE         , dimension(0:12)  ::  njLeap  ! Nb of added Days for Leap Year !
      data njLeap / 0   ,   0 ,   0 ,   1 ,   1 ,   1 ,   1            &
     &                  ,   1 ,   1 ,   1 ,   1 ,   1 ,   1/


!     CAUTION: values in the 3rd column are purely indicative

      end module Mod_PHY____dat
