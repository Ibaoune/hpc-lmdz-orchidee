      module Mod_PHY_CM_dat

!--------------------------------------------------------------------------+
!                                                     Sun 30-Jun-2013  MAR |
!     module Mod_PHY_CM_dat contains the      data                   of    |
!                Cloud Microphsics            Scheme                       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Sat 23-Mar-2013      |
!           Last Modification by H. Gallee,           Sun 30-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE


      integer,       parameter               ::  npt_CM = 1        ! Nb           of     OUTPUT Points
      integer, SAVE                          ::  ipt_CM            ! Index        of     OUTPUT Points
      integer, SAVE, dimension(npt_CM)       ::  i0__CM            ! x-Coordinate of the OUTPUT Point
      integer, SAVE, dimension(npt_CM)       ::  j0__CM            ! y-Coordinate of the OUTPUT Point
      integer, SAVE, dimension(npt_CM)       ::  k0__CM            ! z-Coordinate of the OUTPUT Point
      integer, SAVE, dimension(npt_CM)       ::  ikl0CM            ! v-Coordinate of the OUTPUT Point
      data                                (i0__CM(ipt_CM), ipt_CM=1,npt_CM) /   1/
      data                                (j0__CM(ipt_CM), ipt_CM=1,npt_CM) /   1/
      data                                (k0__CM(ipt_CM), ipt_CM=1,npt_CM) /   1/

      real(kind=real8), SAVE                 ::  qv_MIN = 3.e-6    ! Minimum        Specific       Humidity     (Ch. Tricot) [kg/kg]
      !real(kind=real8), SAVE                 ::  qh_MIN = 1.e-18   ! Minimum        Hydrometeors   Concentration             [kg/kg]
      real(kind=real8), SAVE                 ::  qh_MIN = 1.e-9   ! Minimum        Hydrometeors   Concentration             [kg/kg]
      real(kind=real8), SAVE                 ::  WatIce =  273.16  ! Computation of Saturation pressure
      real(kind=real8), SAVE                 ::  ExpWat =    5.138 ! Ref.: Dudhia        1989, JAS
      real(kind=real8), SAVE                 ::  ExpWa2 = 6827.    !

      real(kind=real8), SAVE                 ::  RH_MAX =   1.0    ! Maximum        Relative       Humidity     (1 = 100%)       [-]
      real(kind=real8), SAVE                 ::  RHcrit =   1.0    !                Sursaturation Threshold     (1 = 100%)       [-]
      real(kind=real8), SAVE                 ::  CFrMIN =   1.e-6  ! Minimum        Cloud Fraction when qw OR qi > 0             [-]
      real(kind=real8), SAVE                 ::  SSImax = 101.0    ! Maximum        Sursaturation % ICE (101 ==> RH= 201%)       [%]
      real(kind=real8), SAVE                 ::  n0___s = 3.0e+6   ! intercept parameter / snow    gamma distribution
      real(kind=real8), SAVE                 ::  n0___r = 8.0e+6   ! intercept parameter / rain    gamma distribution
      real(kind=real8), SAVE                 ::  n0___g = 4.0e+4   ! intercept parameter / graupel gamma distribution
                                                                   ! (Lin et al.   1983, JCAM 22, p.1068: 1, 2 and 3) 
      real(kind=real8), SAVE                 ::  qi0_DC = 0.10e-3  ! Ice   Crystals Critical Mixing Ratio    (tuned Dome C)  [kg/kg]
!                                                qi0_DC = 0.30e-3  ! Ice   Crystals Critical Mixing Ratio        (standard)  [kg/kg]
      real(kind=real8), SAVE                 ::  qs__D0 = 2.00e-4  ! Smallest Diameter of Particles in the snow Class            [m]
!                                                                  ! Ref.: Levkov et al. 1992, Contr.Atm.Phys.65, p.41, para 1
      real(kind=real8), SAVE                 ::  T_NuId =  -5.     ! Deposition & Condensation-Freezing Nucleation Temperature [dgC]
      real(kind=real8), SAVE                 ::  a_NuId=   -0.639  ! Deposition & Condensation-Freezing Nucleation Parameters
      real(kind=real8), SAVE                 ::  b_NuId =   0.1296 ! Ref.: Meyers et al. 1992                     p.713
      real(kind=real8), SAVE                 ::  T_NuIc =  -2.     ! Contact Freezing                   Nucleation Temperature [dgC]
      real(kind=real8), SAVE                 ::  a_NuIc=   -2.80   ! Contact Freezing                   Nucleation Parameters
      real(kind=real8), SAVE                 ::  b_NuIc =   0.262  ! Ref.: Meyers et al. 1992                     p.713

      real(kind=real8), SAVE                 ::  Di_Hex = 1.1e+4   ! Diameter of Hexagonal Plates
                                                                   ! Ref.: Emde & Kalig  1989, Ann.Geophys.    7, p.408  (15): D

      real(kind=real8), SAVE                 ::  TmaxHM =  -3.     !                                                           [dgC]
      real(kind=real8), SAVE                 ::  TminHM =  -8.     !                                                           [dgC]
      real(kind=real8), SAVE                 ::  wa__HM =   1.     ! V. Wind > wa__HM  => Nucleation II (Hall-Mossop) may occur[m/s]

      real(kind=real8), SAVE                 ::  TqwFrz =-35.e0    !       T < TqwFrz  => Instantaneous Freezing               [dgC]
                                                                   ! Ref.: Levkov et al. 1992, Contr.Atm.Phys.65, p.39
      real(kind=real8), SAVE                 ::  qw_VOL = 18.e-15  ! Typical Cloud Droplet Volume [m3] (typ. diam.: 32.5 mim)
                                                                   !      OR Cloud Droplet Weight [T]

      real(kind=real8), SAVE                 ::  qisMAX = 0.0010   ! Ice   Crystals MAX Concentration before autoconversion  [kg/kg]
                                                                   ! Ref.: Lin et al.    1983, JCAM           22, p.1070 (21)
      real(kind=real8), SAVE                 ::  qigMAX = 0.0006   ! Ice   Crystals MAX Concentration before autoconversion  [kg/kg]
                                                                   ! Ref.: Lin et al.    1983, JCAM           22, p.1074 (37)
      real(kind=real8), SAVE                 ::  qwTURB = 0.27     ! qwTURB=1/3 ln(1/k), k=0.8 (droplets dispersion parameter)
                                                                   ! Ref.: Martin et al. 1994, JAS 51, p.1823
      real(kind=real8), SAVE                 ::  C1_EkM = 0.14e-3  ! Partial Condensation Scheme
      real(kind=real8), SAVE                 ::  C2_EkM = 9.75e+0  ! Ref.: Ek and Mahrt 1991, An.Geoph. 9, 716--724

      real(kind=real8), SAVE                 ::  qw_MAX = 0.0001   ! Cloud Droplets MAX Concentration before autoconversion  [kg/kg]
      real(kind=real8), SAVE                 ::  qw_MAXL= 0.0020   ! Cloud droplets MAX concentration before autoconversion
      real(kind=real8), SAVE                 ::  rwCrit = 10.0e-6  ! Droplets Autoconversion:    Critical Radius
                                                                   ! Ref.: Liou and Ou, 1989
      real(kind=real8), SAVE                 ::  c_Sund =  1.0e-4  ! Droplets Autoconversion: 1/(characteristic time scale)
                                                                   ! Ref.: Sundqvist,   1988

      real(kind=real8), SAVE                 ::  cc1 = 1.200e-04   ! Cloud droplets autoconversion parameter
      real(kind=real8), SAVE                 ::  cc2 = 1.569e-12   ! Cloud droplets autoconversion parameter
      real(kind=real8), SAVE                 ::  dd0 = 0.15e0      ! Cloud droplets autoconversion parameter
                                                                   ! Ref.: Lin et al.    1983, JCAM           22, p.1076 (50)




!     ======================================================================
!     Bergeron Process Data (given by Koenig, 1971, J.A.S. 28,p235) ========

      real(kind=real8), SAVE, dimension(31)  ::  aa1

      data aa1/0.7939e-07 , 0.7841e-06 , 0.3369e-05 , 0.4336e-05 ,     &
     &         0.5285e-05 , 0.3728e-05 , 0.1852e-05 , 0.2991e-06 ,     &
     &         0.4248e-06 , 0.7434e-06 , 0.1812e-05 , 0.4394e-05 ,     &
     &         0.9145e-05 , 0.1725e-06 , 0.3348e-04 , 0.1725e-04 ,     &
     &         0.9175e-05 , 0.4412e-05 , 0.2252e-05 , 0.9115e-06 ,     &
     &         0.4876e-06 , 0.3473e-06 , 0.4758e-06 , 0.6306e-06 ,     &
     &         0.8573e-06 , 0.7868e-06 , 0.7192e-06 , 0.6513e-06 ,     &
     &         0.5956e-06 , 0.5333e-06 , 0.4834e-06 /


      real(kind=real8), SAVE, dimension(31)  ::  aa2

      data aa2/0.4006e0,  0.4831e0,  0.5320e0,  0.5307e0,  0.5319e0,   &
     &         0.5249e0,  0.4888e0,  0.3894e0,  0.4047e0,  0.4318e0,   &
     &         0.4771e0,  0.5183e0,  0.5463e0,  0.5651e0,  0.5813e0,   &
     &         0.5655e0,  0.5478e0,  0.5203e0,  0.4906e0,  0.4447e0,   &
     &         0.4126e0,  0.3960e0,  0.4149e0,  0.4320e0,  0.4506e0,   &
     &         0.4483e0,  0.4460e0,  0.4433e0,  0.4413e0,  0.4382e0,   &
     &         0.4361e0/

!     Bergeron Process Data (given by Koenig, 1971, J.A.S. 28,p235) ========
!     ======================================================================

      

      end module Mod_PHY_CM_dat
