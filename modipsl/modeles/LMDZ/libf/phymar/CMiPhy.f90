      subroutine CMiPhy

!------------------------------------------------------------------------------+
!                                                         Mon 17-Jun-2013  MAR |
!   MAR          CMiPhy                                                        |
!     subroutine CMiPhy contains the MAR    Cloud Microphysical Scheme         |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Thu 21-Mar-2013      |
!           Last Modification by H. Gallee,               Mon 17-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+
!                                                                              |
!   INPUT / OUTPUT: qv__DY(kcolp,mzp) : air   specific humidity        [kg/kg] |
!   ^^^^^^^^^^^^^^  qw__CM(kcolp,mzp) : cloud drops                    [kg/kg] |
!                   qi__CM(kcolp,mzp) : ice   crystals Concentration   [kg/kg] |
!                   qs__CM(kcolp,mzp) : Snow  Particl. Concentration   [kg/kg] |
!                   qr__CM(kcolp,mzp) : Rain  Drops    Concentration   [kg/kg] |
!                                                                              |
!   (to be added:   qg__CM(kcolp,mzp) : Graupels       Concentration   [kg/kg])|
!                                                                              |
!                   CCNwCM(kcolp,mzp) : cloud droplets number          [Nb/m3] |
!                   CCNiCM(kcolp,mzp) : ice   crystals number          [Nb/m3] |
!                                                                              |
!                   CFraCM(kcolp,mzp) : cloud fraction                    [-]  |
!                                                                              |
!                   RainCM(kcolp    ) : rain  Precipitation           [m w.e.] |
!                   SnowCM(kcolp    ) : snow  Precipitation           [m w.e.] |
!                   Ice_CM(kcolp    ) : ice   Precipitation           [m w.e.] |
!                                                                              |
!                   qid_CM(kcolp,mzp) : Ice    Water Formation         [kg/kg] |
!                   qwd_CM(kcolp,mzp) : Liquid Water Formation         [kg/kg] |
!                                                                              |
!   INPUT :         wa__DY(kcolp,mzp) : Vertical Wind Speed              [m/s] |
!   ^^^^^           roa_DY(kcolp,mzp) : Air Density                     [T/m3] |
!                   qvsiCM(kcolp,mzp) : Satur.specific humidity (ICE)  [kg/kg] |
!                   qvswCM(kcolp,mzp) : Satur.specific humidity (LIQ.) [kg/kg] |
!                                                                              |
!   OUTPUT:                                                                    |
!   ^^^^^^                                                                     |
!                                                                              |
!                                                                              |
!   REFER. : 1) Ntezimana, unpubl.thes.LLN,          115 pp,     1993          |
!   ^^^^^    2) Lin et al.       JCAM            22, 1065--1092, 1983          |
!               (very similar, except that graupels are represented)           |
!            3) Emde and Kahlig, Annal.Geophys.   7,  405-- 414, 1989          |
!            4) Levkov et al.,   Contr.Atm.Phys. 65,   35--  57, 1992          |
!            5) Meyers et al.,   JAM             31,  708-- 731, 1992          |
!               (Primary Ice-Nucleation Parameterization)                      |
!            6) Delobbe and Gallee, BLM          89,   75-- 107  1998          |
!               (Partial Condensation Scheme)                                  |
!                                                                              |
! # OPTIONS: #hy  Additional IF/THEN/ELSE added precluding vectorization       |
! # ^^^^^^^  #qg  Graupel Conservation Equation         (to  verify & include) |
! #          #cn  Intercept Parameter / snow gamma distribution = fct(Ta)      |
! #          #kk  Limitation of SCu fraction                                   |
! #          #VW  Cloud Drop. Sediment.(Duynkerke .. 1995, JAS  52, p.2763     |
! #          #pp  Emde & Kahlig Ice Crystal Deposition  (not included)         |
! #          #wi  QSat modified by  qw/qi   Ratio       (not included)         |
!                                                                              |
! # DEBUG:   #WH  Additional Output (Each Process  is detailled)               |
! # ^^^^^    #WQ  FULL       Output (Each Process  is detailled)               |
! #          #wH  Additional Output                                            |
! #          #wh  Additional Output (Include write CMiPhy_Debug.h)             |
! #          #EW  Additional Output (Energy and Water Conservation)            |
! #          #ew  Additional Output (Energy and Water Conservation)            |
!                                                                              |
!   REMARK : the sign '~' indicates that reference must be verified            |
!   ^^^^^^^^                                                                   |
!   CAUTION:     Partial Condensation Scheme NOT validated                     |
!   ^^^^^^^      for SCu -- Cu Transition                                      |
!                erf fonction is erroneous on HP                               |
!                                                                              |
!------------------------------------------------------------------------------+



!  Global Variables
!  ================

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_CM_dat
      use Mod_PHY_CM_grd
      use Mod_PHY_CM_kkl
      use Mod_PHY_DY_kkl
      use Mod_PHY_AT_kkl



!  Local  Variables
!  ================

      use Mod_CMiPhy_loc


      IMPLICIT NONE


      logical  ::  Heter_Freezng = .FALSE.  !  .TRUE. => Levkov et al. (1992)    Heterogeneous  Freezing of Cloud Droplets
      logical  ::  Homog_Sublima = .FALSE.  !  .TRUE. => Levkov et al. (1992)    Homogeneous Sublimation of Cloud Ice Particles
      logical  ::  Meyers        = .TRUE.   !  .TRUE. => Meyers et al. (1992)    Ice Nucleation I
      logical  ::  AUTO_w_Sundqv = .TRUE.   !  .TRUE. => Sundqvist     (1988)    Autoconversion Cloud Droplets --> Rain
      logical  ::  AUTO_w_LiouOu = .FALSE.  !  .TRUE. => Liou and Ou   (1989)    Autoconversion Cloud Droplets --> Rain (Tropical SCu ONLY)
      logical  ::  AUTO_w_LinAll = .FALSE.  !  .TRUE. => Lin    et al. (1983)    Autoconversion Cloud Droplets --> Rain
      logical  ::  AUTO_i_Levkov = .TRUE.   !  .TRUE. => Levkov et al. (1992)    Autoconversion Cloud Ice      --> Snow
      logical  ::  AUTO_i_LevkXX = .TRUE.   !  .TRUE. => Levkov et al. (1992)    Autoconversion Cloud Ice      --> Snow (Deposition.Growth)
      logical  ::  AUTO_i_EmdeKa = .FALSE.  !  .TRUE. => Emde & Kahlig (1989)    Autoconversion Cloud Ice      --> Snow
      logical  ::  AUTO_i_Sundqv = .FALSE.  !  .TRUE. => Sundqvist               Autoconversion Cloud Ice      --> Snow
                                            ! .FALSE. => Emde & Kahlig (1989)    Bergeron-Findeisen Process
      logical  ::  fracSC        = .FALSE.  !  .TRUE. => Delobbe                 SCu Fractional Cloudiness may be set up if Frac__Clouds = .TRUE.
      logical  ::  fraCEP        = .FALSE.  !  .TRUE. => ECMWF                   SCu Fractional Cloudiness
      logical  ::  HalMos        = .TRUE.   !  .TRUE. => Levkov et al. (1992)    Ice Nucleation II (Hallet-Mossop Process)
      logical  ::  graupel_shape = .TRUE.   !  .TRUE. => Snow Particles Shape:   Graupellike Snow Flakes of Hexagonal Type
      logical  ::  planes__shape = .FALSE.  !  .TRUE. => Snow Particles Shape:   Unrimed Side Planes
      logical  ::  aggrega_shape = .FALSE.  !  .TRUE. => Snow Particles Shape:   Aggregates of unrimed radiating assemblages

      logical  ::  NO_Vec        = .TRUE.   !  .TRUE. => Preference of IF/THEN/ELSE to sign and max/min Functions

      real(kind=real8)                             :: rad_ww            ! Cloud Droplets Radius                                     [...]
      real(kind=real8)                             :: qvs_wi            ! Saturation Specific Humididy over Liquid Water          [kg/kg]
                                                                        ! Ref.: Emde & Kahlig 1989, Ann.Geophys.    7, p.407 (5)
      real(kind=real8)                             :: BNUCVI            ! Nucleation  I: Deposition & Condensation-Freez. Nucl.   [kg/kg]
      real(kind=real8)                             :: SSat_I            ! Nucleation  I: Sursaturation % ICE                          [%]

      real(kind=real8)                             :: Flag_T_NuId       ! Flag: T < T_NuId  => Nucleation I (Dep./Cond.)  may occur   [-]
      real(kind=real8)                             :: CCNiId            !                                                         [Nb/m3]

      real(kind=real8)                             :: Flag___NuIc       ! Flag: 1           => Nucleation I (Cont.Frz.)   may occur   [-]
      real(kind=real8)                             :: Flag_T_NuIc       ! Flag: T < T_NuIc  => Nucleation I (Cont.Frz.)   may occur   [-]
      real(kind=real8)                             :: qw__OK            ! Flag: 1           => Non-zero cloud droplets Concentration
      real(kind=real8)                             :: qi__OK            ! Flag: 1           => Non-zero cloud ice      Concentration
      real(kind=real8)                             :: CCNiOK            ! Flag: 1           => Non-zero cloud ice      Particles     
      real(kind=real8)                             :: CCNiIc            !                                                         [Nb/m3]

      real(kind=real8)                             :: Flag_TmaxHM       ! Flag: T < TmaxHM  => Nucleation II (Hall-Mossop) may occur  [-]
      real(kind=real8)                             :: Flag_TminHM       ! Flag: T > TminHM  => Nucleation II (Hall-Mossop) may occur  [-]
      real(kind=real8)                             :: Flag_wa__HM       ! Flag: w > wa__HM  => Nucleation II (Hall-Mossop) may occur  [-]
      real(kind=real8)                             :: SplinJ            ! Hallet-Mossop Theory (Levkov et al., 1992,
      real(kind=real8)                             :: SplinP            !                       Contr.Atm.Phy. 65, p.40)

      real(kind=real8)                             :: Flag_Ta_Neg       ! Flag: T < 0.0dgC                                            [-]
      real(kind=real8)                             :: Flag_TqwFrz       ! Flag: T < Temper. => Instantaneous Freezing may occur       [-]

      real(kind=real8)                             :: RHumid            ! Relative Humidity                                           [-]

      real(kind=real8)                             :: qi_Nu1,qi_Nu2     ! Nucleation ( 0 > Ta > -35 dgC ), Generations 1 and 2    [kg/kg]
      real(kind=real8)                             :: qi_Nuc            ! Nucleation ( 0 > Ta > -35 dgC ), Generation (effective) [kg/kg]
      real(kind=real8)                             :: NuIdOK            ! Nucleation  I:              Contact     -Freez. Nucl.       [-]
      real(kind=real8)                             :: BSPRWI            ! Nucleation  I:              Contact     -Freez. Nucl.   [kg/kg]
      real(kind=real8)                             :: BHAMWI            ! Nucleation II: Hallett-Mossop Ice-Multiplic. Process    [kg/kg]
      real(kind=real8)                             :: BNUFWI            ! Heterogeneous Freezing of Cloud Droplets                [kg/kg]
      real(kind=real8)                             :: BDEPVI            ! Ice   Crystals Sublim.                                  [kg/kg]
      real(kind=real8)                             :: Flag_SURSat       ! Flag = 1/0 if (SUR/sub)Saturation                  FLAG     [-]
      real(kind=real8)                             :: Flag_SUBSat       ! Flag = 1/0 if (SUB/sur)Saturation                  FLAG     [-]
      real(kind=real8)                             :: Flag_Sublim       ! Flag = 1/0 for Sublimation Occurence or not                 [-]
      real(kind=real8)                             :: RH_Ice            ! Relative Humidity   vs   Ice                                [-]
      real(kind=real8)                             :: RH_Liq            ! Relative Humidity   vs   Water                              [-]
      real(kind=real8)                             :: DenDi1,DenDi2     ! dqi/dt|sublimation:  terms of the denominator             [...]
      real(kind=real8)                             :: dqsiqv            ! SUBsaturation  qsi - qv                                 [kg/kg]
      real(kind=real8)                             :: dqvSUB            ! SURsaturation  qv  - qsiEFF                             [kg/kg]
      real(kind=real8)                             :: dqvDUM            ! Water Vapor    Variation, dummy variable                [kg/kg]
      real(kind=real8)                             :: dqiDUM            ! Ice   Crystals Variation, dummy variable                [kg/kg]

      real(kind=real8)                             :: qCloud            ! qw + qi                                                 [kg/kg]
      real(kind=real8)                             :: coefC2            ! Coefficient of Ek and Mahrt parameter C2_EkM            [../..]
      real(kind=real8)                             :: pa_hPa,es_hPa     ! Pressure,     Pressure of Vapor at Saturat., over Water [hPa]
      real(kind=real8)                             :: Qsat_L            ! Specific Concentration of Vapor at Saturat., over Water [kg/kg]
                                                                        !       (even for temperatures smaller than freezing pt)
      real(kind=real8)                             :: t_qvqw            ! qv+qw mixing ratio used in Partial Condensation Scheme  [kg/kg]
      real(kind=real8)                             :: d_qvqw            ! qv+qw mixing ratio variation                            [kg/kg]
      real(kind=real8)                             :: Kdqvqw            ! qv+qw vertical turbulent Flux                           [kg/kg m/s]
      real(kind=real8)                             :: ww_TKE            !       vertical Velocity  Variance                       [m2/s2]
      real(kind=real8)                             :: RH_TKE            !       Relative Humidity  Variance                       [../..]
      real(kind=real8)                             :: qt_TKE            ! qv+qw                    Variance                       [../..]
      real(kind=real8)                             :: TLiqid            ! Liquid Temperature                                          [K]
      real(kind=real8)                             :: CFr_t1,CFr_t2     !                                                             [-]
      real(kind=real8)                             :: CFrCoe            !                                                             [-]
      real(kind=real8)                             :: CFraOK            !                                                             [-]
      real(kind=real8)                             :: qwCFra            ! CloudAveraged Liquid Water Mixing Ratio                 [kg/kg]
      real(kind=real8)                             :: qwMesh            ! Mesh Averaged Liquid Water Mixing Ratio                 [kg/kg]
      real(kind=real8)                             :: dwMesh            ! Mesh Averaged Liquid Water Mixing Ratio Variation       [kg/kg]
      real(kind=real8)                             :: signdw            ! Sign       of Liquid Water Mixing Ratio Variation  (-1/1)   [-]
      real(kind=real8)                             :: Flag_dqwPos       ! Flag = 1/0 if Liquid Water Mixing Ratio Variation    >/< 0  [-]
      real(kind=real8)                             :: updatw            !                                                             [-]
      real(kind=real8)                             :: SCuLim            ! Fraction      Limit                                         [-]
      real(kind=real8)                             :: ARGerf            ! Argument of erf function (used in partial condensation Scheme)
      real(kind=real8)                             :: OUTerf            ! OUTPUT   of erf function
      real(kind=real8)                             ::    erf            !             erf function
      real(kind=real8)                             :: dwTUR4,dwTURi     !
      real(kind=real8)                             :: dwTUR3,dwTUR2     !
      real(kind=real8)                             :: dwTUR8,dwTURc     !
      real(kind=real8)                             :: dwTUR5,dwTUR1     !

      real(kind=real8)                             :: Di_Pri            ! Pristine Ice Diameter
      real(kind=real8)                             :: c1saut,cnsaut     ! Ice   Crystals AUToconv.  Parameters
      real(kind=real8)                             :: dtsaut            ! Ice   Crystals AUToconv.  Time Scale
      real(kind=real8)                             :: ps_AUT            ! Ice   Crystals AUToconv. (BDEPIS, BAGRIS,...)      Rate [kg/kg/s]
      real(kind=real8)                             :: qs_AUT            ! Ice   Crystals AUToconv. (BDEPIS, BAGRIS,...)           [kg/kg]

      real(kind=real8)                             :: Flag_Ta_Pos       ! Flag = 1/0 if  Ta >/< 273.15 K                     FLAG     [-]
      real(kind=real8)                             :: Flag_qiMELT       ! Flag = 1/0 for qi Melting or not                   FLAG     [-]
      real(kind=real8)                             :: qxMelt            ! Potential qi to Melt                                    [kg/kg]
      real(kind=real8)                             :: qiMELT            ! Effective qi to Melt                                    [kg/kg]
      real(kind=real8)                             :: CiMelt            ! Effective CCNi removed   by qi   Melting                [nb/m3]

      real(kind=real8)                             :: Flag_qsMELT       ! Flag = 1/0 for qs Melting or not                   FLAG     [-]
      real(kind=real8)                             :: xCoefM            !
      real(kind=real8)                             :: AcoefM,BcoefM     !
      real(kind=real8)                             :: dTMELT            ! Temperature Offset before   Snow Particles Melting          [K]
      real(kind=real8)                             :: qsMELT            ! Melt                     of Snow Particles              [kg/kg]

      real(kind=real8)                             :: qs_ACW            ! Accretion of Cloud Drop. by Snow Particl.               [kg/kg]
      real(kind=real8)                             :: Flag_qs_ACW       ! Accretion of Cloud Drop. by Snow Particl.          FLAG     [-]
      real(kind=real8)                             :: Flag_qs_ACI       ! Accretion of Cloud Ice   by Snow Particl.          FLAG     [-]
      real(kind=real8)                             :: effACI            ! Accretion of Cloud Ice   by Snow Particl.    Efficiency     [-]
      real(kind=real8)                             :: ps_ACI            ! Accretion of Cloud Ice   by Snow Particl.          Rate [kg/kg/s]
      real(kind=real8)                             :: qs_ACI            ! Accretion of Cloud Ice   by Snow Particl.               [kg/kg]
      real(kind=real8)                             :: CNsACI            ! Accretion of Cloud Ice   by Snow Particl. CCNi decrease [nb/m3]
      real(kind=real8)                             :: Flag_qs_ACR       ! Accretion of Snow        by Rain                   FLAG     [-]
      real(kind=real8)                             :: coeACR            ! Accretion of Snow        by Rain Sedimentation Coeffic.   [...]
      real(kind=real8)                             :: qs_ACR            ! Accretion of Snow        by Rain                        [kg/kg]
      real(kind=real8)                             :: qs_ACR_R          ! Accretion of Snow        by Rain          => Rain       [kg/kg]
      real(kind=real8)                             :: qs_ACR_S          !                                           => Graupels   [kg/kg]
      real(kind=real8)                             :: Flag_qr_ACS       ! Accretion of Rain        by Snow Particl. => Snow  FLAG     [-]
      real(kind=real8)                             :: coeACS            ! Accretion of Rain        by Snow Sedimentation Coeffic.   [...]
      real(kind=real8)                             :: pr_ACS            ! Accretion of Rain        by Snow Particl. => Snow  Rate [kg/kg/s]
      real(kind=real8)                             :: qr_ACS            ! Accretion of Rain        by Snow Particl. => Snow       [kg/kg]
      real(kind=real8)                             :: qr_ACS_S          ! Accretion of Rain        by Snow Particl. => Snow       [kg/kg]
      real(kind=real8)                             :: ps_SUB            ! Deposition/Sublim.    on/of Snow Particl.          Rate [kg/kg/s]
      real(kind=real8)                             :: qs_SUB            ! Deposition/Sublim.    on/of Snow Particl.               [kg/kg]
      real(kind=real8)                             :: ls_NUM            ! Sublimation of Snow: NUMerator   of Snow Distribut.Coef.  [...]
      real(kind=real8)                             :: ls_DEN            ! Sublimation of Snow: DENominator of Snow Distribut.Coef.  [...]

      real(kind=real8)                             :: Flag_Freeze       ! Freezing                 of Rain                   FLAG     [-]   
      real(kind=real8)                             :: ps_FRZ            ! Freezing                 of Rain                   Rate [kg/kg/s]
      real(kind=real8)                             :: qs_FRZ            ! Freezing                 of Rain                        [kg/kg]

      real(kind=real8)                             :: rwMEAN            ! Droplets Autoconversion:    Mean     Radius

      real(kind=real8)                             :: th_AUT            ! Cloud Droplets AUToconv. Threshold                          [-]
      real(kind=real8)                             :: pr_AUT            ! Cloud Droplets AUToconv.                           Rate [kg/kg/s]
      real(kind=real8)                             :: qr_AUT            ! Cloud Droplets AUToconv.                                [kg/kg]
      real(kind=real8)                             :: Flag_qr_ACW       ! Accretion of Cloud Drop. by Rain                   FLAG     [-]
      real(kind=real8)                             :: pr_ACW            ! Accretion of Cloud Drop. by Rain                   Rate [kg/kg/s]
      real(kind=real8)                             :: qr_ACW            ! Accretion of Cloud Drop. by Rain                        [kg/kg]
      real(kind=real8)                             :: Flag_qr_ACI       ! Accretion of Cloud Drop. by Rain          => Rain  FLAG     [-]
      real(kind=real8)                             :: pr_ACI            ! Accretion of Cloud Ice   by Rain          => Rain  Rate [kg/kg/s]
      real(kind=real8)                             :: qr_ACI            ! Accretion of Cloud Ice   by Rain          => Rain       [kg/kg]
      real(kind=real8)                             :: CNrACI            ! Accretion of Cloud Ice   by Rain          CCNi decrease [nb/m3]
      real(kind=real8)                             :: pi_ACR            ! Accretion of Cloud Ice   by Rain          => Snow  Rate [kg/kg/s]
      real(kind=real8)                             :: qi_ACR            ! Accretion of Cloud Ice   by Rain          => Snow       [kg/kg]
      real(kind=real8)                             :: Flag_DryAir       ! 1 => RH_Liq < 1                                    FLAG     [-]
      real(kind=real8)                             :: Flag_qr_EVP       ! Evaporation of Rain                                FLAG     [-]
      real(kind=real8)                             :: lr_NUM            ! Evaporation of Rain: NUMerator   of rain Distribut.Coef.  [...]
      real(kind=real8)                             :: lr_DEN            ! Evaporation of Rain: DENominator of rain Distribut.Coef.  [...]
      real(kind=real8)                             :: pr_EVP            ! Evaporation of Rain                                Rate [kg/kg/s]
      real(kind=real8)                             :: qr_EVP            ! Evaporation of Rain                                     [kg/kg]

      real(kind=real8)                             :: effACS            ! Accretion of Snow        by Graupels         Efficiency     [-]

      real(kind=real8)                             :: a_rodz            ! Air               Mass
      real(kind=real8)                             :: qwFlux            ! qw  Sedimentation Flux          Coefficient
      real(kind=real8)                             :: qwrodz            ! Droplets          Mass
      real(kind=real8)                             :: wRatio            ! Droplets          Mass          Ratio                       [-]
      real(kind=real8)                             :: qiFlux            ! qi  Sedimentation Flux          Coefficient
      real(kind=real8)                             :: qirodz            ! Crystals          Mass
      real(kind=real8)                             :: iRatio            ! Crystals          Mass          Ratio                       [-]
      real(kind=real8)                             :: qsFlux            ! qs  Sedimentation Flux          Coefficient
      real(kind=real8)                             :: qsrodz            ! Snow              Mass
      real(kind=real8)                             :: sRatio            ! Snow              Mass          Ratio                       [-]
      real(kind=real8)                             :: qrFlux            ! qr  Sedimentation Flux          Coefficient
      real(kind=real8)                             :: qrrodz            ! Rain              Mass
      real(kind=real8)                             :: rRatio            ! Rain              Mass          Ratio                       [-]

      real(kind=real8)                             :: Vw_MAX            ! MAX Sedimentation Velocity   of Droplets
      real(kind=real8)                             :: Flag_Fall_i       !     Sedimentation            of Ice                FLAG     [-]
      real(kind=real8)                             :: Vi_MAX            ! MAX Sedimentation Velocity   of Ice  Particles
      real(kind=real8)                             :: Vs_MAX,VsMMAX     ! MAX Sedimentation Velocity   of Snow Particles
      real(kind=real8)                             :: Vs__OK            !     Sedimentation Velocity   of Snow Particles
      real(kind=real8)                             :: Vr_MAX,VrMMAX     ! MAX Sedimentation Velocity   of Rain Drops
      real(kind=real8)                             :: Vr__OK            !     Sedimentation Velocity   of Rain Drops
      real(kind=real8)                             :: dtPrec            !     Sedimentation Time Step
      real(kind=real8)                             :: d_Snow            !     Sedimented                  Snow Part. over 1 time step [m]
      real(kind=real8)                             :: d_Rain            !     Sedimented                  Rain Drops over 1 time step [m]

      real(kind=real8)                             :: SatiOK            !
      real(kind=real8)                             :: FlagNu            ! 1 =>Nucleation for -35 dgC < Ta < 0 dgc

      real(kind=real8)                             :: Qw0_OK            ! HYDROMETEORS INPUT STATUS
      real(kind=real8)                             :: Qi0_OK            !
      real(kind=real8)                             :: Qi0qOK            !
      real(kind=real8)                             :: Ci0_OK            !
      real(kind=real8)                             :: Ci0cOK            !
      real(kind=real8)                             :: Qs0_OK            !
      real(kind=real8)                             :: Qr0_OK            !

      real(kind=real8)                             :: qiBerg            ! Bergeron-Findeisen Process: avaiklable qi
      real(kind=real8)                             :: qwBerg            ! Bergeron-Findeisen Process: avaiklable qw
      real(kind=real8)                             :: qxBerg            ! Bergeron-Findeisen Process: potential qi Generation
      real(kind=real8)                             :: a1Berg,a2Berg     ! Bergeron-Findeisen Process: parameters
      real(kind=real8)                             :: a0Berg,afBerg     ! Bergeron-Findeisen Process: integration constants

      real(kind=real8)                             :: WaterB            ! Vertically Integrated Water Budget

      real(kind=real8)                             :: argEXP            !

      integer                                      :: k     ,ikl        !
      integer                                      :: i_Berg            !
      integer                                      :: nItMAX,itFall     !
      integer                                      :: ikl_io,io__Pt     !


!  Debug Variables
!  ---------------

! #wH character(len=70)                            :: debugH            !
! #wH character(len=10)                            :: proc_1,proc_2     !
! #wH character(len=10)                            :: proc_3,proc_4     !
! #wH real(kind=real8)                             :: procv1,procv2     !
! #wH real(kind=real8)                             :: procv3,procv4     !
! #wH integer                                      :: kv    ,nl         !




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! ALLOCATION                                                            !
! ==========                                                            !

      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                             THEN !

          allocate ( qw_io0(      mzp) )                                ! Droplets   Concentration entering CMiPhy                  [kg/kg]
          allocate ( qi_io0(      mzp) )                                ! Ice  Part. Concentration entering CMiPhy                  [kg/kg]
          allocate ( qs___0(kcolp,mzp) )                                ! Snow Part. Concentration entering CMiPhy                  [kg/kg]
! #qg     allocate ( qg___0(kcolp,mzp) )                                ! Graupels   Concentration entering CMiPhy                  [kg/kg]
          allocate ( qr___0(kcolp,mzp) )                                ! Rain Drops Concentration entering CMiPhy                  [kg/kg]
          allocate ( Ta_dgC(kcolp,mzp) )                                ! Air   Temperature                                           [dgC]
          allocate ( sqrrro(kcolp,mzp) )                                ! sqrt(roa(mzp)/roa(k))                                         [-]
          allocate ( qsiEFF(kcolp,mzp) )                                ! EFFective Saturation Specific Humidity over Ice           [kg/kg]
          allocate ( Fletch(kcolp,mzp) )                                ! Monodisperse Nb of hexagonal Plates, Fletcher (1962)          [-]
          allocate ( lamdaS(kcolp,mzp) )                                ! Marshall-Palmer distribution parameter for Snow Particl.
! #qg     allocate ( lamdaG(kcolp,mzp) )                                ! Marshall-Palmer distribution parameter for Graupels
          allocate ( lamdaR(kcolp,mzp) )                                ! Marshall-Palmer distribution parameter for Rain Drops
          allocate ( ps_ACR(kcolp,mzp) )                                ! Accretion of Snow        by Rain                   Rate [kg/kg/s]
          allocate ( ps_ACW(kcolp,mzp) )                                ! Accretion of Cloud Drop. by Snow Particl.          Rate [kg/kg/s]
          allocate ( FallVw(kcolp,mzp) )                                ! Sedimentation Velocity   of Droplets
          allocate ( FallVi(kcolp,mzp) )                                ! Sedimentation Velocity   of Ice  Particles
          allocate ( FallVs(kcolp,mzp) )                                ! Sedimentation Velocity   of Snow Particles
! #qg     allocate ( FallVg(kcolp,mzp) )                                ! Sedimentation Velocity   of Snow Particles
          allocate ( FallVr(kcolp,mzp) )                                ! Sedimentation Velocity   of Rain Drops 
          allocate ( qwLoss(      mzp) )                                ! Mass Loss related to Sedimentation of Rain Droplets
          allocate ( qiLoss(      mzp) )                                ! Mass Loss related to Sedimentation of Ice  Crystals
          allocate ( qsLoss(      mzp) )                                ! Mass Loss related to Sedimentation of Snow Particles
          allocate ( qrLoss(      mzp) )                                ! Mass Loss related to Sedimentation of Rain Drops
! #wH     allocate ( debugV(      mzp,16) )                             ! Debug Variable (of 16 microphysical processes)

! #WH     allocate ( wihm1(mzp) )                                       ! Cloud Droplets Freezing
! #WH     allocate ( wihm2(mzp) )                                       ! Ice   Crystals Homogeneous Sublimation
! #WH     allocate ( wicnd(mzp) )                                       ! Ice   Crystals Nucleation              (Emde & Kahlig)
! #WH     allocate ( widep(mzp) )                                       ! Ice   Crystals Growth Bergeron Process (Emde & Kahlig)
! #WH     allocate ( wisub(mzp) )                                       ! Ice   Crystals             Sublimation (Levkov)
! #WH     allocate ( wimlt(mzp) )                                       ! Ice   Crystals Melting  
! #WH     allocate ( wwevp(mzp) )                                       ! Water Vapor Condensation / Evaporation (Fractional Cloudiness)
! #WH     allocate ( wraut(mzp) )                                       ! Cloud Droplets AUTO-Conversion
! #WH     allocate ( wsaut(mzp) )                                       ! Ice   Crystals AUTO-Conversion
! #WH     allocate ( wracw(mzp) )                                       ! Accretion of Cloud Droplets by Rain, Ta > 0, --> Rain
! #WH     allocate ( wsacw(mzp) )                                       ! Accretion of Cloud Droplets by Rain, Ta < 0, --> Snow
! #WH     allocate ( wsaci(mzp) )                                       ! Accretion of Ice   Crystals by Snow          --> Snow
! #WH     allocate ( wraci(mzp) )                                       ! Accretion of Ice   Crystals by Rain          --> Snow
! #WH     allocate ( wiacr(mzp) )                                       ! Accretion of Rain by Ice   Crystals          --> Snow
! #WH     allocate ( wsacr(mzp) )                                       ! Accretion of Rain by Snow                    --> Snow
! #WH     allocate ( wracs(mzp) )                                       ! Accretion of Snow by Rain                    --> Snow, Rain
! #WH     allocate ( wrevp(mzp) )                                       ! Rain  Drops     Evaporation  
! #WH     allocate ( wssub(mzp) )                                       ! Snow  Particles Sublimation
! #WH     allocate ( wsmlt(mzp) )                                       ! Snow  Particles Melting 
! #WH     allocate ( wsfre(mzp) )                                       ! Rain  Drops     Freezing

      END IF                                                            !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!     +++++++++++++++++++
      DO      ikl=1,kcolp
!     +++++++++++++++++++




!  Debug
!  ~~~~~
! #wH   DO k=mz1_CM,mzp
! #wH     debugH( 1:35) = 'CMiPhy: Debugged Variables: Initial'
! #wH     debugH(36:70) = '                                   '
! #wH     proc_1        = 'R.Hum W[%]'
! #wH     procv1        =  0.1*qv__DY(ikl,k)/(RHcrit * qvswCM(ikl,k))
! #wH     proc_2        = 'R.Hum I[%]'
! #wH     procv2        =  0.1*qv__DY(ikl,k)/(RHcrit * qvsiCM(ikl,k))
! #wH     proc_3        = '          '
! #wH     procv3        =  0.
! #wH     proc_4        = '          '
! #wH     procv4        =  0.

! #wh     include 'CMiPhy_Debug.h'

! #wH     DO kv=1,16
! #wH     debugV(k,kv)  =  0.
! #wH     ENDDO
! #wH   END DO




!  Vertical Integrated Energy and Water Content
!  ============================================

! #EW     enr0EW(ikl) = 0.0
! #EW     wat0EW(ikl) = 0.0

! #EW   DO k=1,mzp
! #EW     enr0EW(ikl) = enr0EW(  ikl)                                   &
! #EW&                +(Ta__CM(ikl,k)                                   &
! #EW&                -(qw__CM(ikl,k)+qr__CM(ikl,k))*Lv_Cpd             &
! #EW&                -(qi__CM(ikl,k)+qs__CM(ikl,k))*Ls_Cpd)            &
! #EW&                * dsigmi(k)
! #EW     wat0EW(ikl) = wat0EW(  ikl)                                   &
! #EW&                +(qv__DY(ikl,k)                                   &
! #EW&                + qw__CM(ikl,k)+qr__CM(ikl,k)                     &
! #EW&                + qi__CM(ikl,k)+qs__CM(ikl,k)        )            &
! #EW&                * dsigmi(k)
! #EW   END DO

! #EW     mphyEW(ikl) ='                    '
!  ..     mphy2D -->   '12345678901234567890'

! #ew     enr0EW(ikl) = enr0EW(ikl) * psa_DY(ikl) * Grav_I
! #EW     wat0EW(ikl) = wat0EW(ikl) * psa_DY(ikl) * Grav_I
!  ..     wat0EW [m]    contains an implicit factor 1.d3 [kPa-->Pa] /ro_Wat


! #WH   VsMMAX = 0.0
! #WH   VrMMAX = 0.0




!  Set lower limit on Hydrometeor Concentration
!  ============================================

! #hy   IF (NO_Vec)                                               THEN

! #hy     DO k=mz1_CM,mzp

! #hy       IF (qw__CM(ikl,k).lt.epsn)                              THEN
! #hy           qv__DY(ikl,k) = qv__DY(ikl,k)+qw__CM(ikl,k)
! #hy           Ta__CM(ikl,k) = Ta__CM(ikl,k)-qw__CM(ikl,k)*Lv_Cpd
! #hy           qwd_CM(ikl,k)=  qwd_CM(ikl,k)-qw__CM(ikl,k)
! #hy           qw__CM(ikl,k) = 0.0
! #hy       END IF

! #hy       IF (qr__CM(ikl,k).lt.epsn)                              THEN
! #hy           qv__DY(ikl,k) = qv__DY(ikl,k)+qr__CM(ikl,k)
! #hy           Ta__CM(ikl,k) = Ta__CM(ikl,k)-qr__CM(ikl,k)*Lv_Cpd
! #hy           qwd_CM(ikl,k) = qwd_CM(ikl,k)-qr__CM(ikl,k)
! #hy           qr__CM(ikl,k) = 0.0
! #hy       END IF

! #hy       IF (qi__CM(ikl,k).lt.epsn.or.CCNiCM(ikl,k).lt.un_1)     THEN
! #hy           qv__DY(ikl,k) = qv__DY(ikl,k)+qi__CM(ikl,k)
! #hy           Ta__CM(ikl,k) = Ta__CM(ikl,k)-qi__CM(ikl,k)*Ls_Cpd
! #hy           qid_CM(ikl,k) = qid_CM(ikl,k)-qi__CM(ikl,k)
! #hy           qi__CM(ikl,k) = 0.0
! #hy           CCNiCM(ikl,k) = 0.0
! #hy       END IF

! #hy       IF (qs__CM(ikl,k).lt.epsn)                              THEN
! #hy           qv__DY(ikl,k) = qv__DY(ikl,k)+qs__CM(ikl,k)
! #hy           Ta__CM(ikl,k) = Ta__CM(ikl,k)-qs__CM(ikl,k)*Ls_Cpd
! #hy           qid_CM(ikl,k) = qid_CM(ikl,k)-qs__CM(ikl,k)
! #hy           qs__CM(ikl,k) = 0.0
! #hy       END IF
! #hy     END DO

! #hy   ELSE

          DO k=mz1_CM,mzp

            Qw0_OK        = max(zer0,sign(un_1,epsn-qw__CM(ikl,k)))*qw__CM(ikl,k)
            qw__CM(ikl,k) =                         qw__CM(ikl,k)  -Qw0_OK
            qwd_CM(ikl,k) =                         qwd_CM(ikl,k)  -Qw0_OK
            qv__DY(ikl,k) =                         qv__DY(ikl,k)  +Qw0_OK
            Ta__CM(ikl,k) =                         Ta__CM(ikl,k)  -Qw0_OK*Lv_Cpd

            Qr0_OK        = max(zer0,sign(un_1,epsn-qr__CM(ikl,k)))*qr__CM(ikl,k)
            qr__CM(ikl,k) =                         qr__CM(ikl,k)  -Qr0_OK
            qwd_CM(ikl,k) =                         qwd_CM(ikl,k)  -Qr0_OK
            qv__DY(ikl,k) =                         qv__DY(ikl,k)  +Qr0_OK
            Ta__CM(ikl,k) =                         Ta__CM(ikl,k)  -Qr0_OK*Lv_Cpd

            Qi0qOK        = max(zer0,sign(un_1,epsn-qi__CM(ikl,k)))
            Ci0cOK        = max(zer0,sign(un_1,un_1-CCNiCM(ikl,k)))

            Ci0_OK        = max(Ci0cOK,Qi0qOK)
            Qi0_OK        =     Ci0_OK*qi__CM(ikl,k)

            CCNiCM(ikl,k) =     Ci0_OK*CCNiCM(ikl,k)
            qi__CM(ikl,k) =            qi__CM(ikl,k) - Qi0_OK
            qid_CM(ikl,k) =            qid_CM(ikl,k) - Qi0_OK
            qv__DY(ikl,k) =            qv__DY(ikl,k) + Qi0_OK
            Ta__CM(ikl,k) =            Ta__CM(ikl,k) - Qi0_OK*Ls_Cpd

            Qs0_OK        = max(zer0,sign(un_1,epsn-qs__CM(ikl,k)))*qs__CM(ikl,k)
            qs__CM(ikl,k) =                         qs__CM(ikl,k)  -Qs0_OK
            qid_CM(ikl,k) =                         qid_CM(ikl,k)  -Qs0_OK
            qv__DY(ikl,k) =                         qv__DY(ikl,k)  +Qs0_OK
            Ta__CM(ikl,k) =                         Ta__CM(ikl,k)  -Qs0_OK*Ls_Cpd

          END DO

! #hy   END IF




!  Initial Concentrations
!  ======================

        DO k=1,mzp
          Ta_dgC(ikl,k) = Ta__CM(ikl,k) -Tf_Sno
          Fletch(ikl,k) = 1.e-2*exp(-0.6*Ta_dgC(ikl,k))                  ! Ice Crystals Number (Fletcher, 1962)

          qr___0(ikl,k) = qr__CM(ikl,k)
          qs___0(ikl,k) = qs__CM(ikl,k)
! #qg     qg___0(ikl,k) = qg__CM(ikl,k)

! #WH     IF (ikl.eq.ikl0CM(1))                                     THEN
! #WH       qw_io0(k)   = qw__CM(ikl,k)
! #WH       qi_io0(k)   = qi__CM(ikl,k)
! #WH     END IF

        END DO




!  Saturation Specific Humidity
!  ============================

        DO k=1,mzp

          qsiEFF(ikl,k) = RHcrit * qvsiCM(ikl,k  )                      ! Saturation Specific Humidity over Ice

          sqrrro(ikl,k) =    sqrt((psa_DY(ikl    )+pt__DY)             &!
     &                           /(roa_DY(ikl,k  )*R_DAir              &!
     &                           * Ta__CM(ikl,mzp)))                    !




!  Hydrometeors   Fall Velocities
!  ==============================

!  Cloud Droplets Fall Velocity (Calcul de la Vitesse Terminale Moyenne)                            FALL VELOCITY
!  ----------------------------

! #VW     IF (qw__CM(ikl,k).ge.epsn)                                THEN

! #VW       CCNwCM(ikl,k) = 1.2d+8                                       ! ASTEX case (Duynkerke et al. 1995, JAS 52, p.2763)

! #VW       qwCFra        =  qw__CM(ikl,k) / max(CFrMIN ,CFraCM(ikl,k))
! #VW       dwTUR4        =   4.5               *qwTURB *qwTURB
! #VW       dwTUR1        =  12.5               *qwTURB *qwTURB
! #VW       dwTURi        =  qwCFra        *     roa_DY(ikl,k)           &
! #VW&                    *  6.0d+0/(piNmbr*CCNwCM(ikl,k)*exp(dwTUR4))
! #VW       dwTUR5        =                    exp(R_5by3*log(dwTURi))

! #VW       FallVw(ikl,k) =  1.19d8* piNmbr*CCNwCM(ikl,k)    *dwTUR5     &
! #VW&                 * exp(dwTUR1)/(24.0 *roa_DY(ikl,k)    *qwCFra)
! #VW     ELSE
! #VW       FallVw(ikl,k) =  0.00
! #VW     END IF


!  Rain           Fall Velocity                                                                     FALL VELOCITY
!  ----------------------------

            lamdaR(ikl,k) = exp(0.25*log((piNmbr*n0___r)               &! Marshall-Palmer Distribution Parameter for Rain
     &              / (roa_DY(ikl,k)*max( epsn  ,qr__CM(ikl,k)))))      ! Ref.: Emde and Kahlig 1989, Ann.Geoph.      7, p.407 (3)
                                                                        ! Note that a simplification occurs
                                                                        ! between the 1000. factor of rho, and rho_water=1000.

! #hy     IF                            (qr__CM(ikl,k).gt.epsn)     THEN

            Vr__OK = max(zer0,sign(un_1, qr__CM(ikl,k)  - epsn))        ! Vr__OK = 1.0 if qr__CM(ikl,k)  > epsn
                                                                        !        = 0.0 otherwise

            FallVr(ikl,k) = Vr__OK*392. *sqrrro(ikl,k)                 &! Terminal Fall Velocity for Rain
     &                    / exp(0.8 *log(lamdaR(ikl,k)))                ! 392  = a Gamma[4+b] / 6
                                                                        ! where  a = 842.  b = 0.8

! #hy     ELSE
! #hy       FallVr(ikl,k)=   0.0
! #hy     END IF


!  Snow Fall Velocity (c and d parameters: see Locatelli and Hobbs, 1974, JGR: table 1 p.2188)      FALL VELOCITY
!  ------------------

! #cn       n0___s        = min(2.e8                                   &
! #cn&                         ,2.e6*exp(-.12*min(0.,Ta_dgC(ikl,k))))

            lamdaS(ikl,k) = exp(0.25*log((0.50*piNmbr*n0___s)          &! Marshall-Palmer distribution parameter for Snow
     &                       / (roa_DY(ikl,k)*max(epsn,qs__CM(ikl,k)))))! Ref.: Emde and Kahlig 1989, Ann.Geoph.      7,  p.407 (3)
                                                                        !       Levkov et al.   1992, Cont.Atm.Phys. 65(1) p.37 (5) (rho_snow)
                                                                        ! Note that a partial simplification occurs
                                                                        ! between the 1000. factor of rho, and rho_snow=500.

! #hy     IF                             (qs__CM(ikl,k).gt.epsn)    THEN

            Vs__OK = max(zer0,sign(un_1,  qs__CM(ikl,k)  - epsn))       ! Vs__OK = 1.0 if qs__CM(ikl,k)  > epsn
                                                                        !        = 0.0 otherwise

           IF      (graupel_shape)                                  THEN
            FallVs(ikl,k) = Vs__OK*2.19  *sqrrro(ikl,k)                &! Terminal Fall Velocity for Graupellike Snow Flakes of Hexagonal Type
     &                     / exp(0.25*log(lamdaS(ikl,k)))               ! 2.19 = c   Gamma[4+d] / 6
                                                                        ! where  c = 4.836,  d =  0.25
                                                                        !          = 0.86 *1000.**0.25
           ELSE IF (planes__shape)                                  THEN
            FallVs(ikl,k) =  Vs__OK*2976.*sqrrro(ikl,k)                &! Terminal Fall Velocity for Unrimed Side Planes
     &                     / exp(0.99*log(lamdaS(ikl,k)))               ! 2976.= c   Gamma[4+d] / 6
                                                                        ! where  c = 755.9,  d  = 0.99
                                                                        !          = 0.81 *1000.**0.99

           ELSE IF (aggrega_shape)                                  THEN
            FallVs(ikl,k) =  Vs__OK*20.06*sqrrro(ikl,k)                &! Terminal Fall Velocity for Aggregates of unrimed radiating assemblages
     &                     / exp(0.41*log(lamdaS(ikl,k)))               ! 2976.= c   Gamma[4+d] / 6
                                                                        ! where  c = 755.9,  d =  0.41
                                                                        !          = 0.69 *1000.**0.41
           ELSE
            STOP   'Snow Particles Shape             is not defined'
           END IF

! #hy     ELSE
! #hy       FallVs(ikl,k)   =   0.0                                      !                          FALL VELOCITY
! #hy     END IF


!  Graupel Fall Velocity
!  ---------------------

! #qg     IF (qg__CM(ikl,k).ge.epsn)                                THEN
! #qg       lamdaG(ikl,k) =exp(0.250*log((piNmbr*n0___g)               &! Marshall-Palmer distribution parameter for Graupel
! #qg&                /(roa_DY(ikl,k)*max(epsn  ,qg__CM(ikl,k)))))      ! Note that a simplification occurs
                                                                        ! between the 1000. factor of rho, and rho_ice=1000.
! #qg       FallVg(ikl,k) = 25.1 *sqrrro(ikl,k)                        &! 25.1 = c Gamma[4+d] / 6
! #qg&             / exp(0.57*log(lamdaG(ikl,k)))                       ! where  c = 4.836 = 1.10 *1000.**0.57 and d = 0.57
!                                                                       ! Hexagonal Graupel, Locatelli and Hobbs, 1974, JGR: table 1 p.2188:

! #qg     ELSE
! #qg       FallVg(ikl,k) =   0.0
! #qg       lamdaG(ikl,k) =   0.0
! #qg     END IF

        END DO


!===============================================================================                    CLOUD ICE  PARTICLES
!                                                                                                   ++++++++++++++++++++
!  Microphysical Processes affecting non Precipitating Cloud Particles
!  ===================================================================

        DO k=mz1_CM,mzp


!  Homogeneous Nucleation by Cloud Dropplets Solidification  ! BFREWI
!  Reference: Emde and Kahlig 1989, Ann.Geoph. 7, p.407 (11) ! Levkov (24) p.40
!  ---------------------------------------------------------

! #wH       Flag_Ta_Neg= 0.
! #wH       qw__OK     = 0.

! #hy    IF                                (Ta_dgC(ikl,k).lt.TqwFrz)THEN

            Flag_TqwFrz=max(zer0,-sign(un_1,Ta_dgC(ikl,k) -  TqwFrz))   ! Flag_TqwFrz = 1.0  if Ta_dgC(ikl,k) < TqwFrz
                                                                        !             = 0.0  otherwise

! #EW       IF(Flag_TqwFrz.gt.eps6)                                THEN
! #EW          mauxEW        =  mphyEW(ikl)
! #EW          mauxEW(01:01) = 'i'
! #EW          mphyEW(ikl)   =  mauxEW
! #EW       END IF

            qw__OK        = qw__CM(ikl,k) *                Flag_TqwFrz
            qi__CM(ikl,k) = qi__CM(ikl,k) +                qw__OK
            CCNiCM(ikl,k) = CCNiCM(ikl,k) + roa_DY(ikl,k) *qw__OK/qw_VOL
            Ta__CM(ikl,k) = Ta__CM(ikl,k) + Lc_Cpd        *qw__OK

! #WQ       write(6,*) 'Qihm1',qw__CM(ikl,k),                          &
! #WQ&                ' Qi'   ,qi__CM(ikl,k),                          &
! #WQ&                ' CcnI' ,CCNiCM(ikl,k),it_EXP,ikl,k
! #WH       IF (ikl.eq.ikl0CM(1))  wihm1(k)  = qw__OK

            qw__CM(ikl,k) = qw__CM(ikl,k) - qw__OK

! #hy    END IF


!  Heterogeneous Freezing of Cloud Droplets                  ! BNUFWI                               CLOUD ICE  PARTICLES
!  Reference: Levkov et al., 1992 (21) p.40                  ! Levkov (21) p.40                     ++++++++++++++++++++
!  ----------------------------------------   

         IF (Heter_Freezng)                                         THEN
! #hy     IF                               (Ta_dgC(ikl,k).lt.0.0)   THEN

            Flag_Ta_Neg=max(zer0,-sign(un_1,Ta_dgC(ikl,k) -  0.0))      ! Flag_Ta_Neg = 1.0 if Ta_dgC(ikl,k) < 0.00dgC
                                                                        !             = 0.0 otherwise

            argEXP = min(max(ea_MIN ,    -Ta_dgC(ikl,k)) ,ea_MAX)
            BNUFWI =     Flag_Ta_Neg*(exp(argEXP)        -1.    )      &!
     &                        * 100.*     qw__CM(ikl,k)  *qw_VOL
            BNUFWI = min(    BNUFWI ,     qw__CM(ikl,k)         )
              
            qi__CM(ikl,k) = qi__CM(ikl,k) +               BNUFWI
            CCNiCM(ikl,k) = CCNiCM(ikl,k) + roa_DY(ikl,k)*BNUFWI/qw_VOL
            Ta__CM(ikl,k) = Ta__CM(ikl,k) + Lc_Cpd       *BNUFWI
            qw__CM(ikl,k) = qw__CM(ikl,k) -               BNUFWI


!  Debug
!  ~~~~~
! #wH             debugH( 1:35)   = 'Homo+Hetero Nucleation by Droplets '
! #wH             debugH(36:70)   = 'Solidification (BFREWI+BNUFWI)     '
! #wH             proc_1          = 'BFREWI    '
! #wH             procv1          =  Flag_TqwFrz
! #wH             proc_2          = 'BNUFWI    '
! #wH             procv2          =  BNUFWI
! #wH             proc_3          = '          '
! #wH             procv3          =  0.
! #wH             proc_4          = '          '
! #wH             procv4          =  0.
! #wh             include 'CMiPhy_Debug.h'
! #wH         IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1)) &
! #wH&            debugV(k,01)    =  qw__OK+BNUFWI

! #hy      END IF
         END IF


!===============================================================================                    CLOUD ICE  PARTICLES
!                                                                                                   ++++++++++++++++++++
!  Homogeneous Sublimation                                   ! XXXXXX
!  Reference: Emde and Kahlig 1989, Ann.Geoph. 7, p.407 (12) ! Levkov
!  ---------------------------------------------------------

! #EW      IF (Flag_TqwFrz.gt.eps6)                                 THEN
! #EW          mauxEW        =  mphyEW(ikl)
! #EW          mauxEW(02:02) = 'I'
! #EW          mphyEW(ikl)   =  mauxEW
! #EW      END IF

         IF   (Homog_Sublima)                                       THEN
               dqvSUB =  (qv__DY(ikl,k)-qsiEFF(ikl,k))                 &!
     &                  /(1.00 +1.733e7*qsiEFF(ikl,k)                  &! 1.733e7=Ls*Ls*0.622/Cpa/Ra with Ls = 2833600 J/kg
     &                  /(Ta__CM(ikl,k)*Ta__CM(ikl,k)))                 !

               dqvSUB =   Flag_TqwFrz*max(zer0,dqvSUB)
               dqvDUM =                        dqvSUB

               qi__CM(ikl,k) = qi__CM(ikl,k) +          dqvDUM
               qid_CM(ikl,k) = qid_CM(ikl,k) +          dqvDUM
!              CCNiCM(ikl,k) : NO VARIATION
               qv__DY(ikl,k) = qv__DY(ikl,k) -          dqvDUM
               Ta__CM(ikl,k) = Ta__CM(ikl,k) + Ls_Cpd * dqvDUM

!  Full Debug
!  ~~~~~~~~~~
! #WQ         write(6,*) 'Qihm2',dqvDUM,                               &
! #WQ&                  ' Qi'   ,qi__CM(ikl,k),                        &
! #WQ&                  ' CcnI' ,CCNiCM(ikl,k),it_EXP,ikl,k
! #WH         if (ikl.eq.ikl0CM(1))  wihm2(k) =   dqvDUM

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)  = 'Emde and Kahlig: Homogeneous Sublim'
! #wH         debugH(36:70)  = 'ation                              '
! #wH         proc_1         = 'dQv   g/kg'
! #wH         procv1         =  dqvDUM
! #wH         proc_2         = '          '
! #wH         procv2         =  0.
! #wH         proc_3         = '          '
! #wH         procv3         =  0.
! #wH         proc_4         = 'CCNI/1.e15'
! #wH         procv4         =  CCNiCM(ikl,k)*1.e-18
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,01)  =  dqvDUM + debugV(k,01)

         END IF
        END DO




!===============================================================================                    CLOUD ICE  PARTICLES, MEYERS
!                                                                                                   ++++++++++++++++++++
!  Nucleation  I: Deposition & Condensation-Freezing Nucleat.
!  Source       : Water Vapor                                ! BNUCVI
!  Reference: Meyers et al., 1992, JAM 31, (2.4) p.712       ! Levkov (20) p.40
!  -----------------------------------------------------------

        IF (Meyers)                                               THEN
          DO k=mz1_CM,mzp

! #wH          NuIdOK      =  0.
! #wH          CCNiId      =  0.
! #wH          BNUCVI      =  0.
! #wH          BSPRWI      =  0.
! #wH          BHAMWI      =  0.

! #hy       IF                                (Ta_dgC(ikl,k).lt.T_NuId) THEN
               Flag_T_NuId=max(zer0,-sign(un_1,Ta_dgC(ikl,k) -  T_NuId))! Flag_T_NuId = 1.0 if Ta_dgC(ikl,k) < T_NuId
!                                                                       !             = 0.0 otherwise

                  dqvDUM  =    qv__DY(ikl,k) - qsiEFF(ikl,k)            ! Sursaturat.

! #hy         IF (dqvDUM.gt.0.)                                     THEN
                  SatiOK  =  max(zer0,sign(un_1,dqvDUM))                ! SatiOK      = 1.0 if qv__DY(ikl,k) > qsiEFF
!                                                                       !             = 0.0 otherwise
                  dqvDUM  =  max(zer0,          dqvDUM)

                  NuIdOK  =  Flag_T_NuId      * SatiOK

                  SSat_I  =  1.e2*dqvDUM      / qsiEFF(ikl,k)           ! Sursaturat.%I
                  SSat_I  =          min(SSat_I,SSImax)                 !
                  CCNiId  =  1.0e3 * exp(a_NuId+b_NuId*SSat_I)          ! Meyers et al. 1992 JAM, 2.4
                  CCNiId  =          max(CCNiId-CCNiCM(ikl,k),zer0)    &!
     &                             *     NuIdOK                         !
                  CCNiCM(ikl,k) =        CCNiId+CCNiCM(ikl,k)           !
                  dqiDUM  =  1.e-15*     CCNiId/roa_DY(ikl,k)           ! 1.e-15  =  0.001 * Initial Ice Crystal Mass
                  dqiDUM        =          min(dqiDUM , dqvDUM)
                  qi__CM(ikl,k)  =              qi__CM(ikl,k) + dqiDUM
                  qid_CM(ikl,k) =               qid_CM(ikl,k) + dqiDUM
                  qv__DY(ikl,k) =               qv__DY(ikl,k) - dqiDUM
                  Ta__CM(ikl,k) =               Ta__CM(ikl,k) + dqiDUM*Ls_Cpd
                  BNUCVI        =                               dqiDUM

! #hy         END IF
! #hy       END IF


!  Nucleation  I:              Contact     -Freezing Nucleat.                                       CLOUD ICE  PARTICLES, MEYERS
!  Source       : Cloud Dropplets                            ! BSPRWI                               ++++++++++++++++++++
!  Reference: Meyers et al., 1992, JAM 31, (2.6) p.713       ! Levkov (20) p.40
!  -----------------------------------------------------------

! #wH             CCNiIc =  0.
! #wH             dqiDUM =  0.

! #hy       IF                               (qw__CM(ikl,k).gt.0.)  THEN
                  qw__OK = max(zer0,sign(un_1,qw__CM(ikl,k)))           ! qw__OK = 1.0 if qw__CM(ikl,k) > 0.
                                                                        !        = 0.0 otherwise

! #hy         IF (Ta_dgC(ikl,k).lt.T_NuIc)                          THEN
                  Flag_T_NuIc   =   max(zer0,-sign(un_1,Ta_dgC(ikl,k) - T_NuIc))
!                 Flag_T_NuIc   =   1.0 if              Ta_dgC(ikl,k) < T_NuIc
!                               =   0.0 otherwise

                  Flag___NuIc   =   Flag_T_NuIc  * qw__OK

                  CCNiIc =   1.e3 *     Flag___NuIc                    &! Contact-Freez
     &                                   * exp(a_NuIc                  &! Potent.Nuclei
     &                                        -b_NuIc                  &! Meyers et al.
     &                                        *Ta_dgC(ikl,k))           ! 1992 JAM, 2.6
                  rad_ww =  (1.e3     * roa_DY(ikl,k)                  &! Drop.  Radius
     &                                * qw__CM(ikl,k)                  &!
     &                                * .2e-11       ) ** 0.33          !
!                 .2e-11 =   1. / (1.2e+8         * 1.e3 * 4.19)
!                                  CCNwCM (ASTEX)   ro_w   4 pi /3
                  CCNiIc = 603.2e+3  *  CCNiIc * rad_ww                &! Levkov et al.
     &                               *  roa_DY(ikl,k)                   ! 1992 CAM,(23)
!                          603.2e3 = 4.0e-7 * 4 pi * 1.2e+8 * 1.e3
!                                    DFar            CCNwCM   fact(rolv)
                  CCNiCM(ikl,k) =       CCNiCM(ikl,k)                  &!
     &                               +  CCNiIc                          !
                  dqiDUM =   1.e-15  *  CCNiIc/roa_DY(ikl,k)
!                            1.e-15  =  1.0e-3 * Ice Crystal Mass
                  dqiDUM =         min( qw__CM(ikl,k) , dqiDUM)
                  qi__CM(ikl,k)  =      qi__CM(ikl,k) + dqiDUM 
                  qw__CM(ikl,k)  =      qw__CM(ikl,k) - dqiDUM
                  Ta__CM(ikl,k)  =      Ta__CM(ikl,k) + dqiDUM*Lc_Cpd
                  BSPRWI         =                      dqiDUM

! #hy         END IF
! #hy       END IF


!  Nucleation II: Hallett-Mossop Ice-Multiplication Process  ! BSPRWI                               CLOUD ICE  PARTICLES, MEYERS
!  Reference: Levkov et al., 1992, Contr.Atm.Ph.65,(25) p.40 ! Levkov (25) p.40                     ++++++++++++++++++++
!  -----------------------------------------------------------

            IF   (HalMos)                                           THEN
! #hy        IF   (Ta_dgC(ikl,k).lt.TmaxHM.AND.                        &
! #hy&             Ta_dgC(ikl,k).gt.TminHM.AND.                        &
! #hy&             wa__DY(ikl,k).gt.wa__HM    )                     THEN
              Flag_TmaxHM = max(zer0,-sign(un_1,Ta_dgC(ikl,k) - TmaxHM))
!             Flag_TmaxHM = 1.0 if              Ta_dgC(ikl,k) < TmaxHM
!                         = 0.0 otherwise

              Flag_TminHM = max(zer0, sign(un_1,Ta_dgC(ikl,k) - TminHM))
!             Flag_TminHM = 1.0 if              Ta_dgC(ikl,k) > TminHM
!                         = 0.0 otherwise

              Flag_wa__HM = max(zer0, sign(un_1,wa__DY(ikl,k) - wa__HM))
!             Flag_wa__HM = 1.0 if              wa__DY(ikl,k) > wa__HM
!                         = 0.0 otherwise

! #cn         n0___s  = min(2.e8,2.e6*exp(-.12*min(Ta_dgC(ikl,k),0.)))

              SplinJ = 1.358e12 *qw__CM(ikl,k)                         &!
     &                          *n0___s          /(lamdaS(ikl,k)**.33)
!                      1.358e12=pi   *Gamma(3.5) *g   *ro_s /(3 *Cd  *4.19e-12)
!                             [=3.14 *3.3233625  *9.81*0.1  /(3 *0.6 *4.19e-12)]
              SplinP = 0.003 * (1. - 0.05 *SplinJ) * Flag_TmaxHM       &!
     &                              * Flag_TminHM  * Flag_wa__HM        !
              SplinP =      max(zer0,      SplinP)

              dqiDUM =          1.e-15  *  SplinP  / roa_DY(ikl,k)      ! 1.e-15  =  1.0e-3 * Ice Crystal Mass
              SplinP = (min(1.0,qs__CM(ikl,k)/max(dqiDUM,epsn))) *SplinP
              CCNiCM(ikl,k) =              CCNiCM(ikl,k)         +SplinP
              dqiDUM =      min(qs__CM(ikl,k),  dqiDUM)
              qi__CM(ikl,k)  =  qi__CM(ikl,k) + dqiDUM 
              qid_CM(ikl,k)  =  qid_CM(ikl,k) + dqiDUM 
              qs__CM(ikl,k)  =  qs__CM(ikl,k) - dqiDUM
              BHAMWI         =                  dqiDUM
! #hy        END IF
            END IF


!  Debug
!  ~~~~~
! #wH           debugH( 1:35)   = 'Meyers: Nucl. I, Depot & Cond-Freez'
! #wH           debugH(36:70)   = 'Nucl. / Freez / Nucl. II / Bergeron'
! #wH           proc_1          = 'dQi1 Meyer'
! #wH           procv1          =  BNUCVI
! #wH           proc_2          = 'dQi2 Meyer'
! #wH           procv2          =  BSPRWI
! #wH           proc_3          = 'dQi Ha-Mos'
! #wH           procv3          =  BHAMWI
! #wH           proc_4          = '          '
! #wH           procv4          =  0.
! #wh           include 'CMiPhy_Debug.h'
! #wH       IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1)) &
! #wH&          debugV(k,02)   =  BNUCVI + BSPRWI + BHAMWI

          END DO




!===============================================================================

        ELSE

!===============================================================================                    CLOUD ICE  PARTICLES, EMDE & KAHLIG
!                                                                                                   ++++++++++++++++++++
!  Ice Crystals Nucleation Process between 0.C and -35.C
!  (each crystal has a mass equal or less than 10d-12 kg)
!  Reference: Emde and Kahlig 1989, Ann.Geoph. 7, p.408 (13)
!  ---------------------------------------------------------

          DO k=mz1_CM,mzp

! #wH        qi_Nu1        =  0.
! #wH        qi_Nu2        =  0.
! #wH        qi_Nuc         =  0.

! #hy       IF                                (Ta_dgC(ikl,k).gt.TqwFrz) THEN
             Flag_TqwFrz = max(zer0, sign(un_1,Ta_dgC(ikl,k)  - TqwFrz))
!            Flag_TqwFrz =     1.0 if          Ta_dgC(ikl,k)  > TqwFrz
!                        =     0.0 otherwise

! #hy       IF                                (Ta_dgC(ikl,k).lt.0.e0  ) THEN
             Flag_Ta_Neg = max(zer0,-sign(un_1,Ta_dgC(ikl,k)          ))
!            Flag_Ta_Neg =     1.0 if          Ta_dgC(ikl,k)  < 0.e0
!                        =     0.0 otherwise

! #hy       IF                                (qv__DY(ikl,k).gt.qsiEFF(ikl,k)) THEN
             SatiOK      = max(zer0, sign(un_1,qv__DY(ikl,k)  - qsiEFF(ikl,k)))
!            SatiOK        =   1.0 if          qv__DY(ikl,k)  > qsiEFF(ikl,k)
!                          =   0.0 otherwise

             FlagNu        =   Flag_TqwFrz * Flag_Ta_Neg * SatiOK

! #EW        IF(FlagNu.gt.eps6)                                     THEN
! #EW           mauxEW        =  mphyEW(ikl)
! #EW           mauxEW(03:03) = 'I'
! #EW           mphyEW(ikl)   =  mauxEW
! #EW        END IF 

             qi_Nu1 = FlagNu * 1.d-15 * Fletch(ikl,k) /roa_DY(ikl,k)
!            qi_Nu1 : amount of nucleated ice crystals (first  condition)

             qi_Nu1 = qi_Nu1*max(zer0,sign(un_1,qi_Nu1-qi__CM(ikl,k)))

             qi_Nu2 = (  qv__DY(ikl,k)-qsiEFF(ikl,k))                  &
     &                 /(1.0d0+1.733d7*qsiEFF(ikl,k)                   &
     &                 /(Ta__CM(ikl,k)*Ta__CM(ikl,k)))
             qi_Nu2 =    FlagNu *  max(zer0  ,qi_Nu2)                   ! amount of nucleated ice crystals (second condition)

             qi_Nuc =              min(qi_Nu1,qi_Nu2)

             qi__CM(ikl,k) = qi__CM(ikl,k) +                qi_Nuc
             qid_CM(ikl,k) = qid_CM(ikl,k) +                qi_Nuc
             CCNiCM(ikl,k) = CCNiCM(ikl,k) + roa_DY(ikl,k) *qi_Nuc     &
     &                                                      *1.e15
             qv__DY(ikl,k) = qv__DY(ikl,k) -                qi_Nuc
             Ta__CM(ikl,k) = Ta__CM(ikl,k) + Ls_Cpd        *qi_Nuc

!  Full Debug
!  ~~~~~~~~~~
! #WQ        write(6,*) 'QiCnd',qi_Nuc,                                &
! #WQ&                 ' Qi'   ,qi__CM(ikl,k),                         &
! #WQ&                 ' CcnI' ,CCNiCM(ikl,k),it_EXP,ikl,k
! #WH        IF (ikl.eq.ikl0CM(1))  wicnd(k) =   qi_Nuc

!  Debug
!  ~~~~~
! #wH            debugH( 1:35)   = 'Emde and Kahlig: Ice Crystals Nucle'
! #wH            debugH(36:70)   = 'ation Process between 0.C and -35.C'
! #wH            proc_1          = 'Qicnd1    '
! #wH            procv1          =  qi_Nu1
! #wH            proc_2          = 'Qicnd2    '
! #wH            procv2          =  qi_Nu2
! #wH            proc_3          = 'Qicnd g/kg'
! #wH            procv3          =  qi_Nuc
! #wH            proc_4          = '          '
! #wH            procv4          =  0.
! #wh            include 'CMiPhy_Debug.h'
! #wH        IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))&
! #wH&           debugV(k,02)    =  qi_Nuc


! #hy       END IF
! #hy       END IF
! #hy       END IF


          END DO

        END IF




!==============================================================================                     CLOUD PARTICLES, MIXED PHASE
!                                                                                                   +++++++++++++++
!  Bergeron Process (water vapor diffusion-deposition on ice crystals)
!  Reference: Koenig          1971, J.A.S.    28, p.235
!             Emde and Kahlig 1989, Ann.Geoph. 7, p.408 (14)
!  ---------------------------------------------------------

        IF (.NOT.AUTO_i_LevkXX)                                     THEN

          DO k=mz1_CM,mzp

! #wH        qi0_OK       =  0.
! #wH        qxBerg       =  0.
! #wH        qwBerg       =  0.

! #hy       IF                                 (qi__CM(ikl,k).gt.epsn
! #hy&         .AND.                            Ta_dgC(ikl,k).lt.0.e0) THEN

              qi0_OK      = max(zer0, sign(un_1,qi__CM(ikl,k)  - epsn))
!             qi0_OK      = 1.0 if              qi__CM(ikl,k)  > epsn
!                         = 0.0 otherwise

              Flag_Ta_Neg = max(zer0,-sign(un_1,Ta_dgC(ikl,k)        ))
!             Flag_Ta_Neg = 1.0 if              Ta_dgC(ikl,k)  < 0.e0
!                         = 0.0 otherwise

              qi0_OK      = Flag_Ta_Neg * qi0_OK

! #EW        IF(qi0_OK.gt.eps6)                                     THEN
! #EW           mauxEW        =  mphyEW(ikl)
! #EW           mauxEW(04:04) = 'i'
! #EW           mphyEW(ikl)   =  mauxEW
! #EW        END IF

              i_Berg = abs(Ta_dgC(ikl,k)-un_1)
              i_Berg = min(i_Berg,31)
              i_Berg = max(i_Berg, 1)
              a1Berg = aa1(i_Berg)
              a2Berg = aa2(i_Berg)

              a0Berg = 1.d+3*roa_DY(ikl,k)*qi__CM(ikl,k) / Fletch(ikl,k)
              afBerg =(a1Berg *(1.0-a2Berg) *  dt__CM                  &! analytical integration of (14) p.408
     &                +a0Berg**(1.0-a2Berg))**(1.0/(1.0-a2Berg))        ! Emde and Kahlig 1989, Ann.Geoph. 7
              qxBerg =(1.d-3*Fletch(ikl,k)/roa_DY(ikl,k))              &!
     &               *(afBerg-a0Berg)                                   !
              qxBerg =     max(zer0,qxBerg)

              qwBerg =     max(zer0,qw__CM(ikl,k))                      ! qwBerg :  to avoid the use of qwd_CM < 0.

              qxBerg = qi0_OK*min(qwBerg,qxBerg)
              qi__CM(ikl,k)=  qi__CM(ikl,k)         +qxBerg
!             CCNiCM(ikl,k):NO VARIATION

              qw__CM(ikl,k)=  qw__CM(ikl,k)         -qxBerg
              Ta__CM(ikl,k)=  Ta__CM(ikl,k)+Lc_Cpd  *qxBerg

!  Full Debug
!  ~~~~~~~~~~
! #WQ         write(6,*) 'QiDep',qxBerg,
! #WQ&                  ' Qi'   ,qi__CM(ikl,k),
! #WQ&                  ' CcnI' ,CCNiCM(ikl,k),it_EXP,ikl,k
! #WH         IF (ikl.eq.ikl0CM(1))  widep(k)= qxBerg

!  Debug
!  ~~~~~
! #wH             debugH( 1:35)   = 'Bergeron Process (water vapor diffu'
! #wH             debugH(36:70)   = 'sion-deposition on ice crystals)   '
! #wH             proc_1          = 'qi0_OK ICE'
! #wH             procv1          =  qi0_OK
! #wH             proc_2          = 'Qicnd g/kg'
! #wH             procv2          =  qwBerg
! #wH             proc_3          = 'Qidep g/kg'
! #wH             procv3          =  qxBerg
! #wH             proc_4          = '          '
! #wH             procv4          =  0.
! #wh             include 'CMiPhy_Debug.h'
! #wH         IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1)) &
! #wH&            debugV(k,02)    =  qxBerg + debugV(k,02)


! #hy       END IF

          END DO

        END IF




!===============================================================================                    CLOUD ICE PARTICLES

!  Ice Crystals Sublimation                                  ! BDEPVI
!  Reference: Emde and Kahlig, 1989 p.408 (15)               ! Levkov (27) p.40
!  -------------------------------------------

        DO k=mz1_CM,mzp

! #wH       BDEPVI      =  0.

! #hy     IF                       (qsiEFF(ikl,k).gt.qv__DY(ikl,k))     THEN

            dqsiqv      =           qsiEFF(ikl,k) -  qv__DY(ikl,k)
! #pp       Flag_SUBSat =  max(zer0,sign(un_1,dqsiqv))
!           Flag_SUBSat =  1.0 if   qsiEFF(ikl,k) >  qv__DY(ikl,k)
!                       =  0.0 otherwise

! #hy     IF                                   (qi__CM(ikl,k).gt.epsn)  THEN

            qi0_OK      =  max(zer0,sign(un_1,  qi__CM(ikl,k)  - epsn))
!           qi0_OK      =  1.0 if               qi__CM(ikl,k)  > epsn
!                       =  0.0 otherwise

            Flag_Sublim =            qi0_OK                            &
! #pp&                  *            Flag_SUBSat                       &
     &                  +            0.0

! #EW      IF(Flag_Sublim.gt.eps6)                                  THEN
! #EW         mauxEW        =  mphyEW(ikl)
! #EW         mauxEW(05:05) = 'V'
! #EW         mphyEW(ikl)   =  mauxEW
! #EW      END IF

            RH_Ice   = qv__DY(ikl,k) /     qsiEFF(ikl,k)
            DenDi1   = 6.959d+11     /    (Ta__CM(ikl,k)*Ta__CM(ikl,k))
!                      6.959e+11
!                    = [Ls=2833600J/kg] * Ls / [kT=0.025W/m/K] / [Rv=461.J/kg/K]
!                                               kT: Air thermal Conductivity
            DenDi2   = 1.0d0 / (1.875d-2*roa_DY(ikl,k)*qsiEFF(ikl,k))
!                               1.875d-5: Water Vapor Diffusivity in Air
            BDEPVI   = dt__CM *(1.-RH_Ice)*4.0*Di_Hex *Fletch(ikl,k)   &!
     &                     /(DenDi1+DenDi2)
            BDEPVI   = max(BDEPVI, -qv__DY(ikl,k))                      ! H2O deposit.limit = H2O content
            BDEPVI   = min(BDEPVI,  qi__CM(ikl,k))                      ! qi  sublim. limit = qi  content
            BDEPVI   = min(BDEPVI,  dqsiqv       )                     &! qi  sublim. limit = Saturation
     &                   * Flag_Sublim

            qi__CM(ikl,k) = qi__CM(ikl,k) -           BDEPVI
            qid_CM(ikl,k) = qid_CM(ikl,k) -           BDEPVI
            qv__DY(ikl,k) = qv__DY(ikl,k) +           BDEPVI
            Ta__CM(ikl,k) = Ta__CM(ikl,k) - Ls_Cpd   *BDEPVI

!  Full Debug
!  ~~~~~~~~~~
! #WQ       write(6,*) 'QiSub', BDEPVI,                                &
! #WQ&                ' Qi'   ,  qi__CM(ikl,k),                        &
! #WQ&                ' CcnI' ,CCNiCM(ikl,k),it_EXP,ikl,k
! #WH       IF (ikl.eq.ikl0CM(1)) wisub(k) = BDEPVI

! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Emde and Kahlig: Ice Crystals Subli'
! #wH         debugH(36:70)   = 'mation                             '
! #wH         proc_1          = 'Qisub g/kg'
! #wH         procv1          =  BDEPVI
! #wH         proc_2          = 'R.Hum I[%]'
! #wH         procv2          =  0.1 * RH_Ice
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1)) &
! #wH&        debugV(k,03)    = -BDEPVI

        END DO

        DO k=mz1_CM,mzp
          IF (qi__CM(ikl,k).le.0.e0)                                THEN
              qi__CM(ikl,k) =  0.e0
              CCNiCM(ikl,k) =  0.e0
          END IF
        END DO




!===============================================================================

!  Ice Crystals Instantaneous Melting
!  ----------------------------------

        DO k=mz1_CM,mzp

! #wH       qiMELT      =  0.
! #wH       CiMelt      =  0.

! #hy     IF                                 (Ta_dgC(ikl,k).gt.0.e0)THEN

            Flag_Ta_Pos = max(zer0, sign(un_1,Ta_dgC(ikl,k)        ))
!           Flag_Ta_Pos = 1.0 if              Ta_dgC(ikl,k) >  0.e0
!                       = 0.0 otherwise

! #hy     IF                                 (qi__CM(ikl,k).gt.epsn)THEN
            qi0_OK      = max(zer0, sign(un_1,qi__CM(ikl,k) -  epsn))
!           qi0_OK      = 1.0 if              qi__CM(ikl,k) >  epsn
!                       = 0.0 otherwise

            Flag_qiMELT = Flag_Ta_Pos * qi0_OK

! #EW      IF(Flag_qiMELT .gt.eps6)                                 THEN
! #EW         mauxEW        =  mphyEW(ikl)
! #EW         mauxEW(06:06) = 'w'
! #EW         mphyEW(ikl)   =  mauxEW
! #EW      END IF

            qxMelt =        Ta_dgC(ikl,k) / Lc_Cpd
            qiMELT =    min(qi__CM(ikl,k) ,         qxMelt)*Flag_qiMELT
            CiMelt =        CCNiCM(ikl,k) *         qiMELT             &
     &                 /max(qi__CM(ikl,k) , epsn)
            qi__CM(ikl,k) = qi__CM(ikl,k) -         qiMELT
            CCNiCM(ikl,k) = CCNiCM(ikl,k) -         CiMelt
            qw__CM(ikl,k) = qw__CM(ikl,k) +         qiMELT
            Ta__CM(ikl,k) = Ta__CM(ikl,k) - Lc_Cpd *qiMELT

!  Full Debug
!  ~~~~~~~~~~
! #WQ       write(6,*) 'QiMlt',qiMELT,                                 &
! #WQ&                ' Qi'   ,qi__CM(ikl,k),                          &
! #WQ&                ' CcnI' ,CCNiCM(ikl,k),it_EXP,ikl,k               
! #WH       IF (ikl.eq.ikl0CM(1))  wimlt(k) =   qiMELT

! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Emde and Kahlig: Ice Crystals Insta'
! #wH         debugH(36:70)   = 'ntaneous Melting                   '
! #wH         proc_1          = 'Qimlt g/kg'
! #wH         procv1          =  qiMELT
! #wH         proc_2          = 'CiMelt /e15'
! #wH         procv2          =  CiMelt*1.e-18
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,04)    = -qiMELT

        END DO




!===============================================================================                    CONDENSATION, Delobbe SCu
!                                                                                                   +++++++++++++++++++++++++
!  Water Vapor Condensation / Evaporation (Fractional Cloudiness)
!  Reference: Laurent Delobbe Thesis (Ek &Mahrt 1991)
!  --------------------------------------------------------------

          DO k=mz1_CM,mzp                                                ! Zeroing needed since
            CFraCM(ikl,k) =  0.0                                        ! a maximization process
          END DO

        IF (Frac__Clouds.AND.fracSC)                                THEN

          DO k=mz1_CM,mzp

! #wH       dwMesh      = 0.

! #hy      IF                               (Ta_dgC(ikl,k).ge.TqwFrz) THEN

            Flag_TqwFrz = max(zer0,sign(un_1,Ta_dgC(ikl,k) -  TqwFrz))
!           Flag_TqwFrz = 1.0 if             Ta_dgC(ikl,k) >  TqwFrz
!                       = 0.0 otherwise

            t_qvqw = qv__DY(ikl,k) +              qw__CM(ikl,k)         ! Total Water Mixing Ratio
            TLiqid = Ta__CM(ikl,k) -  Lv_Cpd    * qw__CM(ikl,k)         ! Liquid Temperature

!  Saturation specific humidity over water,
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ corresponding to liquid temperature
!                               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            pa_hPa =(psa_DY(ikl)  * sigma(k) + pt__DY) * 10.0d0         ! Dudhia (1989) JAS, (B1) and (B2) p.3103
            es_hPa = 6.1078d0 * exp (ExpWat*  log(WatIce     /TLiqid)) &! (see also Pielke (1984), p.234, and
     &                        * exp (ExpWa2*(un_1/WatIce-un_1/TLiqid))  !           Stull  (1988), p.276

            Qsat_L = .622d0*es_hPa /(pa_hPa - .378d0*es_hPa)            ! Saturation Vapor Specific Concentration over Water
                                                                        ! (even for temperatures less than freezing point)

!  Partial Condensation/Scheme
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            d_qvqw = qv__DY(ikl,MIN(k+1,mzp))-qv__DY(ikl,k)            &
     &             + qwd_CM(ikl,MIN(k+1,mzp))-qw__CM(ikl,k)
            Kdqvqw = Kzh_AT(ikl,k)*d_qvqw/(Z___DY(ikl,k+1)-Z___DY(ikl,k))

            ww_TKE  = 0.66d+0 * TKE_AT(ikl,k)                           ! Vertical Velocity Variance

            coefC2 = Kdqvqw/(sqrt(ww_TKE)*Qsat_L)
            RH_TKE = C1_EkM + C2_EkM * coefC2 * coefC2                  ! Relative Humidity Variance
                                                                        ! (Ek and Mahrt, 1991, An. Geoph., 9, 716--724)
 
            qt_TKE  =         sqrt(RH_TKE)*Qsat_L                       ! Total    Water    Variance

            ARGerf = (t_qvqw-Qsat_L)/(1.414d+0*qt_TKE)
            OUTerf = erf(ARGerf)

            CFraCM(ikl,k) = 0.5d+0 * (1.d+0 + OUTerf)                   ! Cloud Fraction

            CFrCoe = 1.d+0/(1.d+0+1.349d7*Qsat_L/(TLiqid*TLiqid))       !
            CFr_t1 = qt_TKE/sqrt(piNmbr+piNmbr)                        &!
     &                * exp(-min(ARGerf*ARGerf,ea_MAX))                 !
            CFr_t2 = CFraCM(ikl,k)    *(t_qvqw-Qsat_L)                  !

            CFraOK =  max(zer0,sign(un_1,CFraCM(ikl,k) - CFrMIN))       ! CFraOK = 1.0 if  CFraCM(ikl,k) > CFrMIN
                                                                        !        = 0.0 otherwise

            CFraCM(ikl,k) = CFraCM(ikl,k) * CFraOK   * Flag_TqwFrz      !
            qwMesh        = CFrCoe * (CFr_t1+CFr_t2) * CFraOK           ! Mesh Averaged Liquid Water Mixing Ratio
            dwMesh        =    qwMesh   -  qw__CM(ikl,k)

!  Vectorisation of the Atmospheric Water Update
!  ~~~~~~~~~~~~~+-------------------------------------------------+
!               |       if (dwMesh.gt.0.d0)             then      |
!               |           dwMesh = min(qv__DY(ikl,k), dwMesh)   |
!               |       else                                      |
!               |           dwMesh =-min(qw__CM(ikl,k),-dwMesh)   |
!               |       end if                                    |
!               +-------------------------------------------------+

            signdw        =    sign(un_1,dwMesh)
            Flag_dqwPos   =     max(zer0,signdw)
            updatw        =    Flag_dqwPos *       qv__DY(ikl,k)       &!
     &                + (1.d0 -Flag_dqwPos)*       qw__CM(ikl,k)        !
! #kk       SCuLim        =        exp(min(0.,300.-Ta__CM(ikl,k)))      ! SCu Lim.
            dwMesh        =    signdw *min(updatw, signdw*dwMesh)      &!
     &                        *Flag_TqwFrz                             &!
! #kk&                                           * SCuLim              &! SCu
     &                    +    0.0
! #kk       CFraCM(ikl,k) =    CFraCM(ikl,k)     * SCuLim               ! Limitor

!  Update of qv__DY, qw__CM and Ta__CM
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            qw__CM(ikl,k) = qw__CM(ikl,k) +             dwMesh
            qwd_CM(ikl,k) = qwd_CM(ikl,k) +             dwMesh
            qv__DY(ikl,k) = qv__DY(ikl,k) -             dwMesh
            Ta__CM(ikl,k) = Ta__CM(ikl,k) + Lv_Cpd    * dwMesh

!  Full Debug
!  ~~~~~~~~~~
! #WQ         write(6,*) 'QwEvp',dwMesh,it_EXP,ikl,k
! #WH         if (ikl.eq.ikl0CM(1)) wwevp(k)     = dwMesh

! #EW        IF(Ta_dgC(ikl,k).ge.TqwFrz)                            THEN
! #EW           mauxEW        =  mphyEW(ikl)
! #EW           mauxEW(07:07) = 'W'
! #EW           mphyEW(ikl)   =  mauxEW
! #EW        END IF

! #hy       END IF

!  Debug
!  ~~~~~
! #wH           debugH( 1:35)   = 'Delobbe: Condensation              '
! #wH           debugH(36:70)   = '                                   '
! #wH           proc_1          = 'dQw   g/kg'
! #wH           procv1          =  dwMesh
! #wH           proc_2          = '          '
! #wH           procv2          =  0.
! #wH           proc_3          = '          '
! #wH           procv3          =  0.
! #wH           proc_4          = '          '
! #wH           procv4          =  0.
! #wh           include 'CMiPhy_Debug.h'
! #wH       IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1)) &
! #wH&          debugV(k,05)   =  dwMesh

          END DO




!===============================================================================                    CONDENSATION, NO SCu
!                                                                                                   ++++++++++++++++++++
!  Water Vapor Condensation / Evaporation
!  Reference: Emde and Kahlig 1989, Ann.Geoph. 7, p.407 (7)
!  --------------------------------------------------------

        ELSE

          DO k=mz1_CM,mzp

! #wH        dwMesh      = 0.

! #hy       IF                                (Ta_dgC(ikl,k).ge.TqwFrz)  THEN
             Flag_TqwFrz = max(zer0, sign(un_1,Ta_dgC(ikl,k) -  TqwFrz))
!            Flag_TqwFrz = 1.0 if              Ta_dgC(ikl,k) >  TqwFrz
!                        = 0.0 otherwise

             dwMesh = (qv__DY(ikl,k)  -qvswCM(ikl,k)*RHcrit)           &
     &                / (1.0d0+1.349d7*qvswCM(ikl,k)                   &
     &                               /(Ta__CM(ikl,k)*Ta__CM(ikl,k)))
!                              1.349e7=Lv*Lv*0.622/Cpa/Ra with Lv = 2500000 J/kg

!  Vectorisation of the Atmospheric Water Update
!  ~~~~~~~~~~~~~+-------------------------------------------------+
!               |       if (dwMesh.gt.0.d0)             then      |
!               |           dwMesh = min(qv__DY(ikl,k), dwMesh)   |
!               |       else                                      |
!               |           dwMesh =-min(qw__CM(ikl,k),-dwMesh)   |
!               |       end if                                    |
!               +-------------------------------------------------+

             signdw      =    sign(un_1,dwMesh)
             Flag_dqwPos =     max(zer0,signdw)
             updatw      =        Flag_dqwPos *    qv__DY(ikl,k)       &
     &                   + (1.d0 -Flag_dqwPos)*    qw__CM(ikl,k)
             dwMesh      =        signdw *min(updatw,signdw*dwMesh)    &
     &                           *Flag_TqwFrz

!  Update of qv__DY, qw__CM and Ta__CM
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             qw__CM(ikl,k) = qw__CM(ikl,k) +             dwMesh
             qwd_CM(ikl,k) = qwd_CM(ikl,k) +             dwMesh
             qv__DY(ikl,k) = qv__DY(ikl,k) -             dwMesh
             Ta__CM(ikl,k) = Ta__CM(ikl,k) + Lv_Cpd    * dwMesh
!            [Ls=2500000J/kg]/[Cp=1004J/kg/K]=2490.04

! #EW        IF(Ta_dgC(ikl,k).ge.TqwFrz)                            THEN
! #EW            mauxEW        =  mphyEW(ikl)
! #EW            mauxEW(07:07) = 'W'
! #EW            mphyEW(ikl)   =  mauxEW
! #EW        END IF

!  Full Debug
!  ~~~~~~~~~~
! #WQ         write(6,*) 'QwEvp',dwMesh,it_EXP,ikl,k
! #WH         if (ikl.eq.ikl0CM(1)) wwevp(k) = dwMesh

! #hy       END IF

!  Debug
!  ~~~~~
! #wH           debugH( 1:35)   = 'Emde and Kahlig: Water Vapor Conden'
! #wH           debugH(36:70)   = 'sation / Evaporation               '
! #wH           proc_1          = 'dQw   g/kg'
! #wH           procv1          =  dwMesh
! #wH           proc_2          = '          '
! #wH           procv2          =  0.
! #wH           proc_3          = '          '
! #wH           procv3          =  0.
! #wH           proc_4          = '          '
! #wH           procv4          =  0.
! #wh           include 'CMiPhy_Debug.h'
! #wH       IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1)) &
! #wH&          debugV(k,05)   =  dwMesh

          END DO

        END IF




!===============================================================================                    CONDENSATION, SCu added AFTER
!                                                                                                   +++++++++++++++++++++++++++++
!  Fractional  Cloudiness ! Guess may be computed (Ek&Mahrt91 fracSC=.T.)
!  ====================== ! Final value  computed  below

! #sc   IF (Frac__Clouds.AND..NOT.fracSC)                           THEN
        IF (Frac__Clouds)                                           THEN
         IF(fraCEP) THEN ! ECMWF Large Scale Cloudiness
                         ! ----------------------------
          DO k=mz1_CM,mzp
              CFraCM(ikl,k) =           (qi__CM(ikl,k) + qw__CM(ikl,k) &!
     &                                  +qs__CM(ikl,k) *  0.33         &!
     &               * (1.-min(un_1,exp((Ta__CM(ikl,k) -258.15)*0.1))))&!
     &               / (0.02       *     qvswCM(ikl,k)                ) !
              CFraCM(ikl,k) =min(un_1  , CFraCM(ikl,k))
              CFraCM(ikl,k) =max(0.001 , CFraCM(ikl,k))                &!
     &             *  max(zer0,sign(un_1,qi__CM(ikl,k) + qw__CM(ikl,k) &!
     &                                  +qs__CM(ikl,k) -3.E-9         ))!
          END DO
         ELSE            ! XU and Randall  1996, JAS 21, p.3099 (4)
                         ! ----------------------------
          DO k=mz1_CM,mzp
              qvs_wi=                                        qvswCM(ikl,k)
! #wi         qvs_wi=max(epsn,((qi__CM(ikl,k)+qs__CM(ikl,k))*qvsiCM(ikl,k)  &
! #wi&                         +qw__CM(ikl,k)               *qvswCM(ikl,k)) &
! #wi&              /max(epsn,  qi__CM(ikl,k)+qs__CM(ikl,k) +qw__CM(ikl,k)))
              RHumid= min( RH_MAX,        max(qv__DY(ikl,k) ,qv_MIN)   &
     &                                      / qvs_wi)
              argEXP=  (  (RH_MAX  -RHumid) * qvs_wi)      **  0.49
              argEXP= min(100.*(qi__CM(ikl,k)+qw__CM(ikl,k)            &
     &                                       +qs__CM(ikl,k) *  0.33    &
     &                 * (1.-min(1.,exp((Ta__CM(ikl,k) -258.15)*0.1))))&
     &                             /max( epsn         , argEXP       ) &
     &                             ,ea_MAX                            )
 
              CFraCM(ikl,k) =      (     RHumid       ** 0.25         )&
     &                         *   (1.  -   exp(-argEXP)              )
          END DO
         END IF

        ELSE
! #sc   ELSE IF (      .NOT.Frac__Clouds)                           THEN
! #sc     IF               (fracSC) stop 'fracSC set up when Frac__Clouds NOT'
          DO k=mz1_CM,mzp
              qCloud        =      qi__CM(ikl,k) + qw__CM(ikl,k)

! #hy       IF                                    (qCloud     &gt.epsn)  THEN
              CFraCM(ikl,k) = max(zer0,sign(un_1,  qCloud       - epsn))
!             CFraCM(ikl,k) = 1.0 if               qCloud       > epsn
!                           = 0.0 otherwise

! #hy       END IF
          END DO

        END IF


!  Debug
!  ~~~~~
! #wH     DO k=mz1_CM,mzp
! #wH         debugH( 1:35) = 'Fractional Cloudiness (XU .OR. CEP)'
! #wH         debugH(36:70) = '                                   '
! #wH         proc_1        = '          '
! #wH         procv1        =  0.
! #wH         proc_2        = '          '
! #wH         procv2        =  0.
! #wH         proc_3        = '          '
! #wH         procv3        =  0.
! #wH         proc_4        = '          '
! #wH         procv4        =  0.     
! #wh         include 'CMiPhy_Debug.h'
! #wH     END DO




!===============================================================================                    AUTO-CONVERSION, LIQUID
!                                                                                                   +++++++++++++++++++++++
!  Autoconversion (i.e., generation of precipitating particles), liquid water
!  ==========================================================================

!  Cloud Droplets Autoconversion
!  Reference: Sundqvist       1988, Schlesinger, Reidel, p.  433)
!  Reference: Lin et al.      1983, JCAM      22, p.1076 (50)
!  ----------------------------------------------------------

        DO k=mz1_CM,mzp

! #wH       qr_AUT = 0.0

! #hy     IF                           (qw__CM(ikl,k).gt.epsn)      THEN
            qw__OK = max(zer0,sign(un_1,qw__CM(ikl,k)  - epsn))
!           qw__OK = 1.0 if             qw__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                           (CFraCM(ikl,k).gt.CFrMIN)    THEN
            CFraOK = max(zer0,sign(un_1,CFraCM(ikl,k)  - CFrMIN))
!           CFraOK = 1.0 if             CFraCM(ikl,k)  > CFrMIN
!                  = 0.0 otherwise

            qw__OK = qw__OK * CFraOK

! #EW      IF(qw__OK.gt.eps6)                                       THEN
! #EW         mauxEW        =  mphyEW(ikl)
! #EW         mauxEW(08:08) = 'r'
! #EW         mphyEW(ikl)   =  mauxEW
! #EW      END IF

!  Sundqvist      (1988, Schlesinger, Reidel, p.  433) Autoconversion Scheme
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           IF      (AUTO_w_Sundqv)                                  THEN
            dwMesh = qw__OK *qw__CM(ikl,k)/qw_MAX                      &!
     &                                /max(CFrMIN,CFraCM(ikl,k))        !
            pr_AUT = qw__OK *qw__CM(ikl,k)*c_Sund                      &!
     &                       *(1.-exp(-min(dwMesh*dwMesh,ea_MAX)))     &!
     &                                /max(CFrMIN,CFraCM(ikl,k))        !

!  Liou and Ou    (1989, JGR  94, p. 8599)             Autoconversion Scheme
!  Boucher et al. (1995, JGR 100, p.16395)             ~~~~~~~~~~~~~~~~~~~~~
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
           ELSE IF (AUTO_w_LiouOu)                                  THEN
            CCNwCM(ikl,k) = 1.2e+8 ! ASTEX (Duynkerke&al.1995, JAS 52, p.2763)
! #lo       CCNwCM(ikl,k) = 1.e+11 !       (polluted air, Rogers&Yau 89, p.90)
      
            qwCFra        = qw__CM(ikl,k) / CFraCM(ikl,k)
            dwTUR4        = 4.5           * qwTURB        *    qwTURB
            dwTURi        = qwCFra        * roa_DY(ikl,k)              &!
     &                    * 6.0 /piNmbr   / CCNwCM(ikl,k) /exp(dwTUR4)
            dwTUR3        =  exp(R_1by3*log(dwTURi))
            dwTUR2        =                 dwTUR3        *    dwTUR3
            dwTUR8        = 8.0           * qwTURB        *    qwTURB
            dwTURc        =   exp(dwTUR8) * dwTUR2        *    dwTUR2
            rwMEAN        = 0.5  *sqrt(sqrt(dwTURc))

            th_AUT        =   max(zer0,sign(un_1,  rwMEAN -rwCrit))     ! Heaviside Function

            pr_AUT        = qw__OK*CFraCM(ikl,k) *th_AUT*4.09d6*piNmbr &!
     &                            *CCNwCM(ikl,k) *dwTURc*qwCFra

!  Lin et al.(1983)                                    Autoconversion Scheme
!  ~~~~~~~~~~~~~~~~                                    ~~~~~~~~~~~~~~~~~~~~~
           ELSE IF (AUTO_w_LinAll)                                  THEN
            dwMesh     = qw__OK * (qw__CM(ikl,k)-qw_MAX)
            pr_AUT     = dwMesh *  dwMesh       *dwMesh/(cc1*dwMesh+1000.d0*cc2/dd0)
           ELSE
            STOP   'AutoConversion of Cloud droplets is not defined'
           END IF

            qr_AUT     = pr_AUT * dt__CM
            qr_AUT     = min(qr_AUT,qw__CM(ikl,k))
            qw__CM(ikl,k) = qw__CM(ikl,k) - qr_AUT
            qr__CM(ikl,k) = qr__CM(ikl,k) + qr_AUT

! #WQ       write(6,*) 'QrAut',qr_AUT,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1)) wraut(k) = qr_AUT

! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Lin et al.(1983) Autoconversion Sch'
! #wH         debugH(36:70)   = 'eme                                '
! #wH         proc_1          = 'Qraut g/kg'
! #wH         procv1          =  qr_AUT
! #wH         proc_2          = '          '
! #wH         procv2          =  0.
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,06)   =  qr_AUT

        END DO




!===============================================================================                    AUTO-CONVERSION, SOLID
!                                                                                                   ++++++++++++++++++++++
!  Autoconversion (i.e., generation of precipitating particles), Ice --> Snow
!  ==========================================================================

!  Conversion from Cloud Ice Crystals to Snow Flakes
!  Reference: Levkov et al.   1992, Contr.Atm.Phys. 65, p.41
!  ---------------------------------------------------------

        IF      (AUTO_i_Levkov)                                     THEN


!  Depositional Growth: Ice Crystals  => Snow Flakes     (BDEPIS)
!  Reference: Levkov et al.   1992, Contr.Atm.Phys. 65, p.41 (28)
!  --------------------------------------------------------------

         IF     (AUTO_i_LevkXX)                                     THEN

          DO k=mz1_CM,mzp

! #wH       qs_AUT = 0.0

! #hy      IF                           (qi__CM(ikl,k).gt.epsn)     THEN
            qi__OK = max(zer0, sign(un_1,qi__CM(ikl,k)  - epsn))
!           qi__OK = 1.0 if              qi__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy      IF                           (CCNiCM(ikl,k).gt.1.e0)     THEN
            CCNiOK = max(zer0, sign(un_1,CCNiCM(ikl,k)  - 1.e0))
!           CCNiOK = 1.0 if              CCNiCM(ikl,k)  > 1.e0
!                  = 0.0 otherwise

            qi__OK = qi__OK * CCNiOK   * qi__CM(ikl,k)

!  Pristine Ice Crystals Diameter
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Di_Pri = 0.156 *exp(R_1by3*log(R_1000*roa_DY(ikl,k)        &! Pristine Ice Crystals Diameter
     &                                 *max(epsn ,qi__CM(ikl,k))       &! Levkov et al. 1992, Contr.Atm.Phys. 65, (5) p.37
     &                                 /max(un_1 ,CCNiCM(ikl,k))))      ! where 6/(pi*ro_I)**1/3 ~ 0.156

!  Deposition Time Scale
!  ~~~~~~~~~~~~~~~~~~~~~
            RH_Ice = max(epsq, qv__DY(ikl,k))   / qsiEFF(ikl,k)

            dtsaut = 0.125   *(qs__D0*qs__D0-Di_Pri*Di_Pri)            &!
     &             *(0.702e12/(Ta__CM(ikl,k)*Ta__CM(ikl,k))            &! 0.702e12 ~ 0.701987755e12 = (2.8345e+6)**2/0.0248/461.5
                                                                        !                              Ls_H2O    **2/Ka    /Rw
     &              +1.0     /(2.36e-2      *roa_DY(ikl,k)             &! 2.36e-2                   =  2.36e-5             *1.e3
     &               *max(epsq,qv__DY(ikl,k))*RH_Ice))                  !                              Dv

!  Deposition
!  ~~~~~~~~~~
            qs_AUT =    dt__CM *qi__OK*(RH_Ice-1.)/dtsaut
            qs_AUT =   min( qi__CM(ikl,k)  , qs_AUT)
            qs_AUT =   max(-qs__CM(ikl,k)  , qs_AUT)
            qi__CM(ikl,k) = qi__CM(ikl,k)  - qs_AUT
            qs__CM(ikl,k) = qs__CM(ikl,k)  + qs_AUT

! #hy      END IF
! #hy      END IF

!  Debug
!  ~~~~~
! #wH          debugH( 1:35)   = 'Lin et al.(1983) Depositional Growt'
! #wH          debugH(36:70)   = 'h                                  '
! #wH          proc_1          = 'QsAUT g/kg'
! #wH          procv1          =  qs_AUT
! #wH          proc_2          = '          '
! #wH          procv2          =  0.
! #wH          proc_3          = '          '
! #wH          procv3          =  0.
! #wH          proc_4          = '          '
! #wH          procv4          =  0.
! #wh          include 'CMiPhy_Debug.h'
! #wH      IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))  &
! #wH&         debugV(k,07)    =  qs_AUT

          END DO

         END IF


!  Ice Crystals Aggregation           => Snow Flakes     (BAGRIS)                                   AUTO-CONVERSION, SOLID
!  Reference: Levkov et al.   1992, Contr.Atm.Phys. 65, p.41 (31)                                   ++++++++++++++++++++++
!  --------------------------------------------------------------

          DO k=mz1_CM,mzp

! #wH       qs_AUT = 0.0
! #wH       dtsaut = 0.0

! #hy      IF                          (qi__CM(ikl,k).gt.epsn)      THEN
            qi__OK = max(zer0,sign(un_1,qi__CM(ikl,k)  - epsn))
!           qi__OK = 1.0 if             qi__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy      IF                          (CCNiCM(ikl,k).gt.1.e0)      THEN
            CCNiOK = max(zer0,sign(un_1,CCNiCM(ikl,k)  - 1.e0))
!           CCNiOK = 1.0 if             CCNiCM(ikl,k)  > 1.e0
!                  = 0.0 otherwise

            qi__OK      =  qi__OK  * CCNiOK * qi__CM(ikl,k)

! #EW       IF(qi__OK.gt.eps6)                                      THEN
! #EW          mauxEW        =  mphyEW(ikl)
! #EW          mauxEW(09:09) = 's'
! #EW          mphyEW(ikl)   =  mauxEW
! #EW       END IF

!  Pristine Ice Crystals Diameter
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Di_Pri = 0.156 *exp(R_1by3*log(R_1000*roa_DY(ikl,k)        &! Pristine Ice Crystals Diameter
     &                                 *max(epsn, qi__CM(ikl,k))       &! Levkov et al. 1992, Contr. Atm. Phys. 65, (5) p.37
     &                                 /max(un_1, CCNiCM(ikl,k))))      ! where [6/(pi*ro_I)]**1/3 ~ 0.156

!  Time needed for Ice Crystals Diameter to reach Snow Diameter Threshold
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c1saut =      max(epsn,qi__OK)    *roa_DY(ikl,k) *35.0     &!
     &        *exp(R_1by3*log(roa_DY(ikl,mzp) /roa_DY(ikl,k)))          !

            dtsaut =-6.d0*log(Di_Pri/qs__D0)  /c1saut                   !
            dtsaut =      max(dt__CM,          dtsaut)                  ! qi fully used if dtsaut<dt__CM

!  Time needed for Ice Crystals Diameter to reach Snow Diameter Threshold
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(ALTERNATE PARAMETERIZATION)~
! #nt       dtsaut =-2.0 *(3.0*log(    Di_Pri                /qs__D0)  &!
! #nt&                    +    log(max(qi__CM(ikl,k),epsn))) /c1saut    !
! #nt       dtsaut = max(epsn,dtsaut)

!  Aggregation
!  ~~~~~~~~~~~
            qs_AUT = dt__CM*qi__OK        / dtsaut
            qs_AUT =   min( qi__CM(ikl,k) , qs_AUT)
            qs_AUT =   max(-qs__CM(ikl,k) , qs_AUT)
            qi__CM(ikl,k) = qi__CM(ikl,k) - qs_AUT
            qs__CM(ikl,k) = qs__CM(ikl,k) + qs_AUT


!  Decrease of Ice Crystals Number                       (BAGRII)
!  Reference: Levkov et al.   1992, Contr.Atm.Phys. 65, p.41 (34)
!  --------------------------------------------------------------

            CCNiCM(ikl,k) = CCNiCM(ikl,k) * exp(-0.5*c1saut*dt__CM)

! #WQ       write(6,*) 'QsAut', qs_AUT,                                &!
! #WQ&                 ' Qi'   ,  qi__CM(ikl,k),                       &!
! #WQ&                 ' CcnI' ,CCNiCM(ikl,k),it_EXP,ikl,k              !
! #WH       if (ikl.eq.ikl0CM(1))   wsaut(k) =   qs_AUT

! #hy      END IF
! #hy      END IF

!  Debug
!  ~~~~~
! #wH          debugH( 1:35)   = 'Lin et al.(1983) Ice Crystals Aggre'
! #wH          debugH(36:70)   = 'gation                             '
! #wH          proc_1          = 'dtsaut sec'
! #wH          procv1          =  dtsaut
! #wH          proc_2          = 'QsAUT g/kg'
! #wH          procv2          =  qs_AUT
! #wH          proc_3          = '          '
! #wH          procv3          =  0.
! #wH          proc_4          = '          '
! #wH          procv4          =  0.
! #wh          include 'CMiPhy_Debug.h'
! #wH      IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))  &
! #wH&         debugV(k,07)   =  qs_AUT + debugV(k,07)

          END DO


!  Ice Crystals Autoconversion => Snow Flakes                                                       AUTO-CONVERSION, SOLID
!  Reference: Lin et al.      1983, JCAM      22, p.1070 (21)                                       ++++++++++++++++++++++
!             Lin et al.      1983, JCAM      22, p.1074 (38)
!             Emde and Kahlig 1989, Ann.Geoph. 7, p. 408 (18)
!  ----------------------------------------------------------

        ELSE IF (AUTO_i_EmdeKa)                                     THEN

          DO k=mz1_CM,mzp

! #wH       qs_AUT      =  0.0
! #wH       cnsaut      =  0.0

! #hy      IF (qi__CM(ikl,k) .ge. qisMAX)                           THEN

            ps_AUT   =      0.001d0*(qi__CM(ikl,k)-qisMAX)             &!
     &                 *exp(0.025d0* Ta_dgC(ikl,k))
            qs_AUT   =     ps_AUT  * dt__CM
            qs_AUT   =  max(qs_AUT,  zer0         )
            qs_AUT   =  min(qs_AUT,  qi__CM(ikl,k))
            cnsaut   =      qs_AUT*  CCNiCM(ikl,k)                     &!
     &                 /max(qisMAX , qi__CM(ikl,k))
            CCNiCM(ikl,k) = CCNiCM(ikl,k) - cnsaut
            qi__CM(ikl,k) = qi__CM(ikl,k) - qs_AUT
            qs__CM(ikl,k) = qs__CM(ikl,k) + qs_AUT
! #WQ       write(6,*) 'QsAut',qs_AUT   ,it_EXP,ikl,k
! #WH       IF (ikl.eq.ikl0CM(1)) wsaut(k)= qs_AUT

! #EW       IF (qi__CM(ikl,k) .ge. qisMAX)                          THEN
! #EW           mauxEW        =  mphyEW(ikl)
! #EW           mauxEW(09:09) = 's'
! #EW           mphyEW(ikl)   =  mauxEW
! #EW       END IF

! #hy      END IF

!  Debug
!  ~~~~~
! #wH          debugH( 1:35)   = 'Emde and Kahlig  Ice Crystals Autoc'
! #wH          debugH(36:70)   = 'onversion                          '
! #wH          proc_1          = 'QsAUT g/kg'
! #wH          procv1          =  qs_AUT
! #wH          proc_2          = 'cnsaut/e15'
! #wH          procv2          =  cnsaut*1.e-18
! #wH          proc_3          = '          '
! #wH          procv3          =  0.
! #wH          proc_4          = '          '
! #wH          procv4          =  0.
! #wh          include 'CMiPhy_Debug.h'
! #wH      IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))  &
! #wH&         debugV(k,07)   =  qs_AUT

          END DO


!  Sundqvist      (1988, Schlesinger, Reidel, p.  433) Autoconversion Scheme                        AUTO-CONVERSION, SOLID
!  -------------------------------------------------------------------------                        ++++++++++++++++++++++

        ELSE IF (AUTO_i_Sundqv)                                     THEN

          DO k=mz1_CM,mzp

! #wH       qs_AUT = 0.0
! #wH       cnsaut = 0.0

! #hy      IF                          (qi__CM(ikl,k).gt.epsn)      THEN
            qi__OK = max(zer0,sign(un_1,qi__CM(ikl,k)  - epsn))
!           qi__OK = 1.0 if             qi__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

            dqiDUM = qi__OK *qi__CM(ikl,k)/qi0_DC                      &!
! #mf&         /max( CFrMIN ,CFraCM(ikl,k))                            &!
     &             + 0.                                                 !
            ps_AUT = qi__OK *qi__CM(ikl,k)*c_Sund                      &!
     &     *(1.-exp(-dqiDUM *dqiDUM))                                  &!
! #mf&         *max( CFrMIN ,CFraCM(ikl,k))                            &!
     &             + 0.                                                 !
            qs_AUT =                        ps_AUT * dt__CM
            qs_AUT =    min(qi__CM(ikl,k) , qs_AUT)
            qs_AUT =    max(zer0          , qs_AUT)
            cnsaut =        CCNiCM(ikl,k) * qs_AUT                     &!
     &                 /max(qi__CM(ikl,k) , epsn)
            CCNiCM(ikl,k) = CCNiCM(ikl,k) - cnsaut
            qi__CM(ikl,k) = qi__CM(ikl,k) - qs_AUT
            qs__CM(ikl,k) = qs__CM(ikl,k) + qs_AUT

! #hy      END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Sundqvist (1988) Ice Crystals Autoc'
! #wH         debugH(36:70)   = 'onversion                          '
! #wH         proc_1          = 'QsAUT g/kg'
! #wH         procv1          =  qs_AUT
! #wH         proc_2          = 'cnsaut/e15'
! #wH         procv2          =  cnsaut*1.e-18
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,07)   =  qs_AUT

          END DO
        ELSE
            STOP   'AutoConversion of Cloud crystals is not defined'
        END IF




!===============================================================================                    AUTO-CONVERSION, SOLID
!                                                                                                   ++++++++++++++++++++++
!  Autoconversion (i.e., generation of precipitating particles), Ice --> Graupels
!  ==============================================================================

! #qg   DO k=mz1_CM,mzp

! #qg     IF (qi__CM(ikl,k) .ge. qigMAX)                              THEN

! #qg       pgaut     = 0.001*(   qi__CM(ikl,k)-qigMAX)*exp(0.090*  Ta_dgC(ikl,k))
! #qg       qgaut     =     pgaut * dt__CM
! #qg       qgaut     = max(qgaut,zer0       )
! #qg       qgaut     = min(qgaut,qi__CM(ikl,k))
! #qg       qi__CM(ikl,k) = qi__CM(ikl,k) - qgaut
! #qg       qg__CM(ikl,k) = qg__CM(ikl,k) + qgaut

! #qg     END IF

! #qg   END DO




!===============================================================================                    ACCRETION
!                                                                                                   +++++++++
!  Accretion Processes (i.e. increase in size of precipitating particles
!  ====================      through a collision-coalescence process)===
!                      ==============================================

!  Accretion of Cloud Droplets by Rain                                                              ACCRETION,  o > .
!  Reference: Lin et al.      1983, JCAM      22, p.1076 (51)                                       +++++++++++++++++
!             Emde and Kahlig 1989, Ann.Geoph. 7, p. 407 (10)
!  ----------------------------------------------------------

        DO k=mz1_CM,mzp

! #wH       qr_ACW = 0.0

! #hy     IF                           (qw__CM(ikl,k).gt.epsn)      THEN
            qw__OK = max(zer0,sign(un_1,qw__CM(ikl,k)  - epsn))
!           qw__OK = 1.0 if             qw__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                           (qr___0(ikl,k).gt.epsn)      THEN
            qr0_OK = max(zer0,sign(un_1,qr___0(ikl,k)  - epsn))
!           qr0_OK = 1.0 if             qr___0(ikl,k)  > epsn
!                  = 0.0 otherwise

            Flag_qr_ACW = qw__OK * qr0_OK

! #EW      IF(Flag_qr_ACW.gt.eps6)                                  THEN
! #EW         mauxEW        =  mphyEW(ikl)
! #EW         mauxEW(10:10) = 'r'
! #EW         mphyEW(ikl)   =  mauxEW
! #EW      END IF

            pr_ACW = 3104.28d0  * n0___r * sqrrro(ikl,k)               &! 3104.28 = a pi Gamma[3+b] / 4
     &        *qw__CM(ikl,k)/exp(3.8d0*log(lamdaR(ikl,k)))              !   where   a = 842. and b  = 0.8
            qr_ACW =        pr_ACW*dt__CM *Flag_qr_ACW
            qr_ACW =    min(qr_ACW,qw__CM(ikl,k))

            qw__CM(ikl,k) = qw__CM(ikl,k) -qr_ACW
            qr__CM(ikl,k) = qr__CM(ikl,k) +qr_ACW

! #WQ       write(6,*) 'Qracw',qr_ACW,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1)) wracw(k) =    qr_ACW

! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Lin et al.(1983): Accretion of Clou'
! #wH         debugH(36:70)   = 'd Droplets by Rain                 '
! #wH         proc_1          = 'Qracw g/kg'
! #wH         procv1          =  qr_ACW
! #wH         proc_2          = '          '
! #wH         procv2          =  0.
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,08)   =  qr_ACW

        END DO


!  Accretion of Cloud Droplets by Snow Flakes                                                       ACCRETION, * > .
!  Reference: Lin et al.      1983, JCAM      22, p.1070 (24)                                       ++++++++++++++++
!  ----------------------------------------------------------

        DO k=mz1_CM,mzp

! #hy     IF                           (qw__CM(ikl,k).gt.epsn)      THEN
            qw__OK = max(zer0,sign(un_1,qw__CM(ikl,k)  - epsn))
!           qw__OK = 1.0 if             qw__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                           (qs___0(ikl,k).gt.epsn)      THEN
            qs0_OK = max(zer0,sign(un_1,qs___0(ikl,k)  - epsn))
!           qs0_OK = 1.0 if             qs___0(ikl,k)  > epsn
!                  = 0.0 otherwise

            Flag_qs_ACW = qw__OK * qs0_OK

! #EW      IF(Flag_qs_ACW.gt.eps6)                                  THEN
! #EW         mauxEW        =  mphyEW(ikl)
! #EW         mauxEW(11:11) = 's'
! #EW         mphyEW(ikl)   =  mauxEW
! #EW      END IF

! #cn       n0___s = min(2.e8,2.e6*exp(-.12*min(0.,Ta_dgC(ikl,k))))

! ps_ACW is taken into account in the snow melting process (if positive temperatures)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           IF      (graupel_shape)                                  THEN! Graupellike Snow Flakes of Hexagonal Type
            ps_ACW(ikl,k)= 9.682d0 * n0___s * sqrrro(ikl,k)            &! 9.682 = c pi  Gamma[3+d] / 4
     &     *qw__CM(ikl,k)     /exp(3.25d0*log(lamdaS(ikl,k)))           ! where   c = 4.836 and d = 0.25
                                                                        ! Ref.: Locatelli and Hobbs, 1974, JGR: table 1 p.2188

           ELSE IF (planes__shape)                                  THEN! Unrimed Side Plane
            ps_ACW(ikl,k)= 3517.   * n0___s * sqrrro(ikl,k)            &! 3517. = c pi  Gamma[3+d] / 4
     &     *qw__CM(ikl,k)     /exp(3.99d0*log(lamdaS(ikl,k)))           ! where   c = 755.9 and d = 0.99

           ELSE IF (aggrega_shape)                                  THEN! Aggregates of unrimed radiating assemblages
            ps_ACW(ikl,k)= 27.73   * n0___s * sqrrro(ikl,k)            &! 27.73 = c pi  Gamma[3+d] / 4
     &     *qw__CM(ikl,k)     /exp(3.41d0*log(lamdaS(ikl,k)))           ! where   c = 11.718and d = 0.41

           ELSE
            STOP   'Snow Particles Shape             is not defined'
           END IF

            qs_ACW =     dt__CM*ps_ACW(ikl,k)*Flag_qs_ACW
            qs_ACW = min(qs_ACW,qw__CM(ikl,k))

            Flag_Ta_Pos = max(zer0,sign(un_1,Ta__CM(ikl,k) - Tf_Sno))
!           Flag_Ta_Pos = 1.0 if             Ta__CM(ikl,k) > Tf_Sno
!                       = 0.0 otherwise

            qw__CM(ikl,k) = qw__CM(ikl,k) -                       qs_ACW
            qr__CM(ikl,k) = qr__CM(ikl,k) +        Flag_Ta_Pos  * qs_ACW
            Flag_qs_ACW   =                (1.d0 - Flag_Ta_Pos) * qs_ACW
            qs__CM(ikl,k) = qs__CM(ikl,k) +                  Flag_qs_ACW
            Ta__CM(ikl,k) = Ta__CM(ikl,k) + Lc_Cpd         * Flag_qs_ACW
!           Negative Temperatures => Latent Heat is released by Freezing

!  Full Debug
!  ~~~~~~~~~~
! #WQ       write(6,*) 'Qsacw',qs_ACW,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1)) wsacw(k) =    qs_ACW

! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Lin et al.(1983): Accretion of Clou'
! #wH         debugH(36:70)   = 'd Droplets by Snow Particles       '
! #wH         proc_1          = 'Qsacw g/kg'
! #wH         procv1          =  Flag_qs_ACW
! #wH         proc_3          = '          '
! #wH         procv2          =  0.
! #wH         proc_2          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,09)   =  Flag_qs_ACW

        END DO


!  Accretion of Cloud Droplets by Graupels (Dry Growth Mode)                                        ACCRETION, # > . | #
!  Reference: Lin et al.      1983, JCAM      22, p.1075 (40)                                       ++++++++++++++++++++
!             Emde and Kahlig 1989, Ann.Geoph. 7, p. 407 (~20)
!  -----------------------------------------------------------

! #qg   DO k=mz1_CM,mzp

! #qg     IF                            (qw__CM(ikl,k).gt.epsn)     THEN
! #qg       WbyG_w = max(zer0, sign(un_1,qw__CM(ikl,k)  - epsn))
!           WbyG_w = 1.0 if              qw__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #qg     IF                            (qg__CM(ikl,k).gt.epsn)     THEN
! #qg       WbyG_g = max(zer0, sign(un_1,qg__CM(ikl,k)  - epsn))
!           WbyG_g = 1.0 if              qg__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #qg       WbyGOK = WbyG_w * WbyG_g

! #qg     IF                            (Ta__CM(ikl,k).lt.Tf_Sno)   THEN
! #qg       Fact_G = max(zer0,-sign(un_1,Ta__CM(ikl,k)  - Tf_Sno))
!           Fact_G = 1.0 if              Ta__CM(ikl,k)  > Tf_Sno
!                  = 0.0 otherwise

! #qg       pgacw  = PATATRAS
! #qg       qgacw  =     pgacw * dt__CM * WbyGOK
! #qg       qgacw  = min(qgacw,qw__CM(ikl,k))

! #qg       qw__CM(ikl,k) = qw__CM(ikl,k) -       qgacw
! #qg       qg__CM(ikl,k) = qg__CM(ikl,k) +       qgacw
! #qg       Ta__CM(ikl,k) = Ta__CM(ikl,k) +Lc_Cpd  gacw

! #qg     END IF
! #qg     END IF
! #qg     END IF

! #qg   END DO




!  Accretion of Cloud Ice      by Snow Particles                                                    ACCRETION, * > /
!  Reference: Lin et al.      1983, JCAM      22, p.1070 (22)                                       ++++++++++++++++
!  ----------------------------------------------------------

        DO k=mz1_CM,mzp

! #wH       qs_ACI = 0.0 
! #wH       CNsACI = 0.0

! #hy     IF                           (qi__CM(ikl,k).gt.epsn)      THEN
            qi__OK = max(zer0,sign(un_1,qi__CM(ikl,k)  - epsn))
!           qi__OK = 1.0 if             qi__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                           (qs___0(ikl,k).gt.epsn)      THEN
            qs0_OK = max(zer0,sign(un_1,qs___0(ikl,k)  - epsn))
!           qs0_OK = 1.0 if             qs___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                                 (Ta__CM(ikl,k).lt.Tf_Sno) THEN
            Flag_Ta_Neg = max(zer0,-sign(un_1,Ta__CM(ikl,k)  - Tf_Sno))
!           Flag_Ta_Neg = 1.0 if              Ta__CM(ikl,k)  < Tf_Sno
!                       = 0.0 otherwise

            Flag_qs_ACI = qi__OK * qs0_OK * Flag_Ta_Neg

! #EW      IF(Flag_qs_ACI.gt.eps6)                                  THEN
! #EW       mauxEW        =  mphyEW(ikl)
! #EW       mauxEW(12:12) = 's'
! #EW       mphyEW(ikl)   =  mauxEW
! #EW      END IF

            effACI = exp(0.025d0*Ta_dgC(ikl,k))                         ! Collection Efficiency
                                                                        ! Lin et al. 1983 JCAM 22 p.1070 (23)

! #cn       n0___s = min(2.e8,2.e6*exp(-.12*min(0.,Ta_dgC(ikl,k))))

! ps_ACI
! ~~~~~~
           IF      (graupel_shape)                                  THEN
            ps_ACI = effACI * 9.682d0 * n0___s * sqrrro(ikl,k)         &
     &     *qi__CM(ikl,k)        /exp(3.25d0*log(lamdaS(ikl,k)))

           ELSE IF (planes__shape)                                  THEN
            ps_ACI = effACI * 3517.d0 * n0___s * sqrrro(ikl,k)         &
     &     *qi__CM(ikl,k)        /exp(3.99d0*log(lamdaS(ikl,k)))

           ELSE IF (aggrega_shape)                                  THEN
            ps_ACI = effACI * 27.73d0 * n0___s * sqrrro(ikl,k)         &
     &     *qi__CM(ikl,k)        /exp(3.41d0*log(lamdaS(ikl,k)))
           ELSE
            STOP   'Snow Particles Shape             is not defined'
           END IF

            qs_ACI =     ps_ACI * dt__CM * Flag_qs_ACI
            qs_ACI = min(qs_ACI,qi__CM(ikl,k))

            CNsACI        = CCNiCM(ikl,k) * qs_ACI                     &!
     &                 /max(qi__CM(ikl,k) , epsn)
            CCNiCM(ikl,k) = CCNiCM(ikl,k) - CNsACI
            qi__CM(ikl,k) = qi__CM(ikl,k) - qs_ACI
            qs__CM(ikl,k) = qs__CM(ikl,k) + qs_ACI

! #WQ       write(6,*) 'Qsaci',qs_ACI,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1))  wsaci(k) =   qs_ACI

! #hy     END IF
! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Lin et al.(1983): Accretion of Clou'
! #wH         debugH(36:70)   = 'd Ice by Snow Particles            '
! #wH         proc_1          = 'Qsaci g/kg'
! #wH         procv1          =  qs_ACI
! #wH         proc_2          = 'CNsaci/e15'
! #wH         procv2          =  CNsACI*1.e-18
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,10)   =  qs_ACI

        END DO




!  Accretion of Cloud Ice      by Graupel (Cloud Ice Sink)                                          ACCRETION, # > /
!  Reference: Lin et al.      1983, JCAM      22, p.1075 (41)                                       ++++++++++++++++
!             Emde and Kahlig 1989, Ann.Geoph. 7, p. 407 (~19)
!  -----------------------------------------------------------

! #qg   DO k=mz1_CM,mzp

! #qg     IF                            (qi__CM(ikl,k).gt.epsn)     THEN
! #qg       CbyG_c = max(zer0, sign(un_1,qi__CM(ikl,k)  - epsn))
!           CbyG_c = 1.0 if              qi__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #qg     IF                            (qg__CM(ikl,k).gt.epsn)     THEN
! #qg       CbyG_g = max(zer0, sign(un_1,qg__CM(ikl,k)  - epsn))
!           CbyG_g = 1.0 if              qg__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #qg     IF                            (Ta__CM(ikl,k).lt.Tf_Sno)   THEN
! #qg       Fact_G = max(zer0,-sign(un_1,Ta__CM(ikl,k)  - Tf_Sno))
!           Fact_G = 1.0 if              Ta__CM(ikl,k)  < Tf_Sno
!                  = 0.0 otherwise

! #qg       CbyGOK = CbyG_c * CbyG_g * Fact_G

! #qg       pgaci = PATATRAS
! #qg       qgaci =     pgaci *dt__CM *CbyGOK
! #qg       qgaci = min(qgaci,qi__CM(ikl,k))

! #qg       qi__CM(ikl,k) = qi__CM(ikl,k) - qgaci
! #qg       qg__CM(ikl,k) = qg__CM(ikl,k) + qgaci

! #qg     END IF
! #qg     END IF
! #qg     END IF

! #qg   END DO


!  Accretion of Cloud Ice      by Rain (Cloud Ice Sink)                                             ACCRETION, o > / | o
!  Reference: Lin et al.      1983, JCAM      22, p.1071 (25)                                       ++++++++++++++++++++
!  ----------------------------------------------------------

        DO k=mz1_CM,mzp

! #wH       qr_ACI = 0.0 
! #wH       qi_ACR = 0.0 

! #hy     IF                           (qi__CM(ikl,k).gt.epsn)      THEN
            qi__OK = max(zer0,sign(un_1,qi__CM(ikl,k)  - epsn))
!           qi__OK = 1.0 if             qi__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                           (qr___0(ikl,k).gt.epsn)      THEN
            qr0_OK = max(zer0,sign(un_1,qr___0(ikl,k)  - epsn))
!           qr0_OK = 1.0 if             qr___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                                 (Ta__CM(ikl,k).lt.Tf_Sno)  THEN
            Flag_Ta_Neg = max(zer0,-sign(un_1,Ta__CM(ikl,k)  - Tf_Sno))
!           Flag_Ta_Neg = 1.0 if              Ta__CM(ikl,k)  < Tf_Sno
!                       = 0.0 otherwise

            Flag_qr_ACI = qi__OK * qr0_OK * Flag_Ta_Neg

! #EW      IF(Flag_qr_ACI.gt.eps6)                                  THEN
! #EW            mauxEW        =  mphyEW(ikl)
! #EW        IF (mauxEW(13:13).eq.'s'.or.mauxEW(13:13).eq.'A')      THEN
! #EW            mauxEW(13:13) =  'A'
! #EW        ELSE
! #EW            mauxEW(13:13) =  'r'
! #EW        END IF
! #EW            mphyEW(ikl)  =  mauxEW
! #EW      END IF

            pr_ACI = 3104.28d0  * n0___r * sqrrro(ikl,k)               &!
     &     *qi__CM(ikl,k)   /exp(3.8d0*log(lamdaR(ikl,k)))
            qr_ACI =     pr_ACI*dt__CM   * Flag_qr_ACI
            qr_ACI = min(qr_ACI,qi__CM(ikl,k))
            CNrACI =            CCNiCM(ikl,k)* qr_ACI/max(qi__CM(ikl,k),epsn)
            CCNiCM(ikl,k) =     CCNiCM(ikl,k)- CNrACI
            qi__CM(ikl,k) =     qi__CM(ikl,k)- qr_ACI

! #qg      IF(qr__CM(ikl,k) .gt. 1.e-4 )                           THEN
! #qg       qg__CM(ikl,k) = qg__CM(ikl,k) + qr_ACI
!           CAUTION : Graupels Formation is not taken into account
!                     This could be a reasonable assumption for Antarctica

! #qg      ELSE
            qs__CM(ikl,k) = qs__CM(ikl,k) + qr_ACI
! #qg      END IF

! #WQ       write(6,*) 'Qraci',qr_ACI,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1))  wraci(k) =   qr_ACI


!  Accretion of Rain           by Cloud Ice (Rain Sink)                                             ACCRETION, / > o | *
!  Reference: Lin et al.      1983, JCAM      22, p.1071 (26)                                       ++++++++++++++++++++
!  ----------------------------------------------------------

! #EW      IF  (Flag_qr_ACI.gt.eps6)                                THEN
! #EW           mauxEW        =  mphyEW(ikl)
! #EW       IF (mauxEW(13:13).eq.'r'.or.mauxEW(13:13).eq.'A')       THEN
! #EW           mauxEW(13:13) =  'A'
! #EW       ELSE
! #EW           mauxEW(13:13) =  's'
! #EW       END IF
! #EW           mphyEW(ikl)   =  mauxEW
! #EW      END IF

            pi_ACR =     4.1d20 * n0___r * sqrrro(ikl,k)               &! 4.1e20 = a pi**2 rhow/mi Gamma[6+b] / 24
     &     *qi__CM(ikl,k)   /exp(6.8d0*log(lamdaR(ikl,k)))              ! where    a=842., rhow=1000, mi=4.19e-13
                                                                        !                                  b = 0.8
                                                                        ! Lin et al, 1983, JAM,p1071: mi:Ice Crystal Mass
            qi_ACR =        pi_ACR*dt__CM * Flag_qr_ACI
            qi_ACR =    min(qi_ACR,qr__CM(ikl,k))
            qr__CM(ikl,k) = qr__CM(ikl,k) -          qi_ACR
            Ta__CM(ikl,k) = Ta__CM(ikl,k) + Lc_Cpd  *qi_ACR

! #qg      IF (qr__CM(ikl,k) .gt. 1.e-4 )                           THEN
! #qg       qg__CM(ikl,k) = qg__CM(ikl,k) +          qi_ACR
!           CAUTION : Graupels Formation is not taken into account
!                     This could be a reasonable assumption for Antarctica

! #qg      ELSE
            qs__CM(ikl,k) = qs__CM(ikl,k) +          qi_ACR
! #qg      END IF

!  Full Debug
!  ~~~~~~~~~~
! #WQ         write(6,*) 'Qiacr',qi_ACR,it_EXP,ikl,k
! #WH         if (ikl.eq.ikl0CM(1)) wiacr(k) =    qi_ACR

! #hy     END IF
! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Lin et al.(1983): Accretion of Clou'
! #wH         debugH(36:70)   = 'd Ice by Rain                      '
! #wH         proc_1          = 'Qraci g/kg'
! #wH         procv1          =  qr_ACI
! #wH         proc_2          = 'qi_ACR g/kg'
! #wH         procv2          =  qi_ACR
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,11)   =  qi_ACR

        END DO




!  Accretion of Rain           by Snow Flakes                                                       ACCRETION o > *, * > o
!  Accretion of Snow Flakes    by Rain                                                              ++++++++++++++++++++++
!  Reference: Lin et al.      1983, JCAM      22, p.1071 (27)
!             Lin et al.      1983, JCAM      22, p.1071 (28)
!             Emde and Kahlig 1989, Ann.Geoph. 7, p. 408 (~21)
!  -----------------------------------------------------------

        DO k=mz1_CM,mzp

            ps_ACR(ikl,k) =   0.0
            qs_ACR        =   0.0
            qs_ACR_S      =   0.0 
            qs_ACR_R      =   0.0 

! #hy     IF                           (qr___0(ikl,k).gt.epsn)      THEN
            qr0_OK = max(zer0,sign(un_1,qr___0(ikl,k)  - epsn))
!           qr0_OK = 1.0 if             qr___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                           (qs___0(ikl,k).gt.epsn)      THEN
            qs0_OK = max(zer0,sign(un_1,qs___0(ikl,k)  - epsn))
!           qs0_OK = 1.0 if             qs___0(ikl,k)  > epsn
!                  = 0.0 otherwise

            Flag_qr_ACS = qr0_OK * qs0_OK

! #EW      IF(Flag_qr_ACI.gt.eps6)                                  THEN
! #EW         mauxEW        =  mphyEW(ikl)
! #EW         mauxEW(14:14) = 'A'
! #EW         mphyEW(ikl)   =  mauxEW
! #EW      END IF

!  Accretion of Rain by Snow --> Snow           | lamdaR : lambda_r
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~           | lamdaS : lambda_s
            coeACS=(5.0d0/(lamdaS(ikl,k)*lamdaS(ikl,k)*lamdaR(ikl,k))  &!
     &             +2.0d0/(lamdaS(ikl,k)*lamdaR(ikl,k)*lamdaR(ikl,k))  &!
     &             +0.5d0/(lamdaR(ikl,k)*lamdaR(ikl,k)*lamdaR(ikl,k))) &!
     &     /(lamdaS(ikl,k)*lamdaS(ikl,k)*lamdaS(ikl,k)*lamdaS(ikl,k))   !

! #cn       n0___s = min(2.e8,2.e6*exp(-.12 *min(0.,Ta_dgC(ikl,k))))

            pr_ACS = 986.96d-3*(n0___r*n0___s/roa_DY(ikl,k))           &!  986.96: pi**2 * rhos
     &                    * abs(FallVr(ikl,k)-FallVs(ikl,k))*coeACS     ! (snow density assumed equal to  100 kg/m3)
            qr_ACS =            pr_ACS       *dt__CM  * Flag_qr_ACS
            qr_ACS =        min(qr_ACS,       qr__CM(ikl,k))

! #WQ       write(6,*) 'Qracs',qr_ACS,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1))  wracs(k) =   qr_ACS

!  Accretion of Snow by Rain --> Rain
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            qr0_OK = max(zer0,sign(un_1,qr___0(ikl,k)  - 1.e-4))
!           qr0_OK = 1.0 if             qr___0(ikl,k)  > 1.e-4
!                  = 0.0 otherwise

            qs0_OK = max(zer0,sign(un_1,qs___0(ikl,k)  - 1.e-4))
!           qs0_OK = 1.0 if             qs___0(ikl,k)  > 1.e-4
!                  = 0.0 otherwise

            Flag_qs_ACR      =   max(qr0_OK,qs0_OK)

! #hy      IF (Flag_qs_ACR.gt.eps6)                                 THEN
            coeACR=(5.0d0/(lamdaR(ikl,k)*lamdaR(ikl,k)*lamdaS(ikl,k))  &
     &             +2.0d0/(lamdaR(ikl,k)*lamdaS(ikl,k)*lamdaS(ikl,k))  &
     &             +0.5d0/(lamdaS(ikl,k)*lamdaS(ikl,k)*lamdaS(ikl,k))) &
     &    /(lamdaR(ikl,k) *lamdaR(ikl,k)*lamdaR(ikl,k)*lamdaR(ikl,k))

            ps_ACR(ikl,k)=9869.6d-3*(n0___r*n0___s/roa_DY(ikl,k))      &!  9869.6: pi**2 * rhow
     &                         * abs(FallVr(ikl,k)-FallVs(ikl,k))*coeACR! (water   density assumed equal to 1000 kg/m3)
            qs_ACR =     ps_ACR(ikl,k)*dt__CM *Flag_qr_ACS  *Flag_qs_ACR
            qs_ACR = min(qs_ACR,qs__CM(ikl,k))

! #WQ       write(6,*) 'Qsacr',qs_ACR,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1)) wsacr(k) =    qs_ACR
! #hy      ELSE
! #hy       ps_ACR(ikl,k) =  0.d0
! #hy       qs_ACR        =  0.d0
! #hy      END IF

            Flag_Ta_Neg = max(zer0,-sign(un_1,Ta__CM(ikl,k)  - Tf_Sno))
!           Flag_Ta_Neg = 1.0 if              Ta__CM(ikl,k)  < Tf_Sno
!                       = 0.0 otherwise

            qr_ACS_S      =                 qr_ACS *      Flag_Ta_Neg
            qs_ACR_R      =                 qs_ACR *(1.d0-Flag_Ta_Neg)
            qr__CM(ikl,k) = qr__CM(ikl,k) - qr_ACS_S
! #qg       IF (qr___0(ikl,k).lt.1.e-4 .and. qs___0(ikl,k).lt.1.e-4)THEN
!              CAUTION  : Graupel Formation is not taken into Account
                qs__CM(ikl,k)  = qs__CM(ikl,k) + qr_ACS_S
! #qg       ELSE²
! #qg           qs__CM(ikl,k)  = qs__CM(ikl,k) - qr_ACS_S
! #qg           qg__CM(ikl,k)  = qg__CM(ikl,k) + qs_ACR_S + qr_ACS_S
! #qg       REND IF
                Ta__CM(ikl,k)  = Ta__CM(ikl,k) + qs_ACR_S * Lc_Cpd

                qr__CM(ikl,k)  = qr__CM(ikl,k) + qs_ACR_R
                qs__CM(ikl,k)  = qs__CM(ikl,k) - qs_ACR_R
                Ta__CM(ikl,k)  = Ta__CM(ikl,k) - qs_ACR_R * Lc_Cpd

! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Lin et al.(1983): Accretion of Snow'
! #wH         debugH(36:70)   = '(Rain) by Rain(Snow)               '
! #wH         proc_1          = 'Qracs g/kg'
! #wH         procv1          =  qs_ACR_S
! #wH         proc_2          = 'Qsacr g/kg'
! #wH         procv2          =  qs_ACR_R
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,12)   =  qs_ACR_S - qs_ACR_R

        END DO




!  Accretion of Snow           by Graupels                                                          ACCRETION, # > *
!  Reference: Lin et al.      1983, JCAM      22, p.1071 (29)                                       ++++++++++++++++
!  ----------------------------------------------------------

! #qg   DO k=mz1_CM,mzp

! #qg     IF                           (qg___0(ikl,k).gt.epsn)      THEN
! #qg       SbyG_g = max(zer0,sign(un_1,qg___0(ikl,k)  - epsn))
!           SbyG_g = 1.0 if             qg___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #qg     IF                           (qs___0(ikl,k).gt.epsn)      THEN
! #qg       SbyG_s = max(zer0,sign(un_1,qs___0(ikl,k)  - epsn))
!           SbyG_s = 1.0 if             qs___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #qg       SbyGOK = SbyG_g *  SbyG_s
! #qg       effACS = exp(0.090*Ta_dgC(ikl,k))                          &! Collection Efficiency
!                                                                       ! Lin et al. 1983 JCAM 22 p.1072 (30)

! #qg       flg=exp(-6.0d0*log(lamdaS(ikl,k))                          &!
! #qg&         *(5.0/lamdaG(ikl,k)                                     &!
! #qg&          +2.0*lamdaS(ikl,k)/(lamdaG(ikl,k)*lamdaG(ikl,k))       &!
! #qg&          +0.5*lamdaS(ikl,k)* lamdaS(ikl,k)                      &!
! #qg&               /exp(3.0d0*log(lamdaG(ikl,k))))                    !

! #cn       n0___s = min(2.e8,2.e6*exp(-.12*min(0.,Ta_dgC(ikl,k))))

! #qg       pgacs  = 986.96d-3*(n0___g*n0___s/roa_DY(ikl,k))           &! 986.96: pi**2 * rhog
! #qg&                    * abs(FallVg(ikl,k)-FallVs(ikl,k))*flg*effACS !(graupel densitity assumed equal to snow density)
! #qg       qgacs  =     pgacs*dt__CM      * SbyGOK
! #qg       qgacs  = min(qgacs,qs__CM(ikl,k))
! #qg       qg__CM(ikl,k)  = qg__CM(ikl,k) + qgacs
! #qg       qs__CM(ikl,k)  = qs__CM(ikl,k) - qgacs

! #qg     END IF
! #qg     END IF

! #qg   END DO


!  Accretion of Rain           by Graupels (Dry Growth Mode)                                        ACCRETION, # > o
!  Reference: Lin et al.      1983, JCAM      22, p.1075 (42)                                       ++++++++++++++++
!  ----------------------------------------------------------

! #qg   DO k=mz1_CM,mzp

! #qg     IF                           (qg___0(ikl,k).gt.epsn)      THEN
! #qg       RbyG_g = max(zer0,sign(un_1,qg___0(ikl,k)  - epsn))
!           RbyG_g = 1.0 if             qg___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #qg     IF                           (qr___0(ikl,k).gt.epsn)      THEN
! #qg       RbyG_r = max(zer0,sign(un_1,qr___0(ikl,k)  - epsn))
!           RbyG_r = 1.0 if             qr___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #qg     IF                            (Ta__CM(ikl,k).lt.Tf_Sno)   THEN
! #qg       Fact_G = max(zer0,-sign(un_1,Ta__CM(ikl,k)  - Tf_Sno))
!           Fact_G = 1.0 if              Ta__CM(ikl,k)  < Tf_Sno
!                  = 0.0 otherwise

! #qg       RbyGOK = RbyG_g * RbyG_s * Fact_G

! #qg       flg=exp(-6.0d0*log(lamdaS(ikl,k))                          &!
! #qg&         *(5.0/lamdaG(ikl,k)                                     &!
! #qg&          +2.0*lamdaS(ikl,k)/(lamdaG(ikl,k)*lamdaG(ikl,k))       &!
! #qg&          +0.5*lamdaS(ikl,k)* lamdaS(ikl,k)                      &!
! #qg&               /exp(3.0d0*log(lamdaG(ikl,k))))                    !

! #cn       n0___s = min(2.e8,2.e6*exp(-.12*min(0.,Ta_dgC(ikl,k))))

! #qg       pgacr  = 986.96d-3*(n0___g*n0___s/roa_DY(ikl,k))           &!
! #qg&                    * abs(FallVg(ikl,k)-FallVr(ikl,k))*flg
! #qg       qgacr  = pgacr    * dt__CM       *RbyGOK
! #qg       qgacr  = min(qgacr,qr__CM(ikl,k))
! #qg       qg__CM(ikl,k) = qg__CM(ikl,k) +        qgacr
! #qg       qr__CM(ikl,k) = qr__CM(ikl,k) -        qgacr
! #qg       Ta__CM(ikl,k) = Ta__CM(ikl,k) + Lc_Cpd*qgacr

! #qg     END IF
! #qg     END IF
! #qg     END IF

! #qg   END DO


!  Graupels Wet Growth Mode
!  Reference: Lin et al.      1983, JCAM      22, p.1075 (43)
!  ----------------------------------------------------------

! #qg   ! TO BE ADDED !




!  Microphysical Processes affecting     Precipitating Cloud Particles
!  ===================================================================


!  Rain Drops Evaporation                                                                           RAIN, EVAPORATION
!  Reference: Lin et al.      1983, JCAM      22, p.1077 (52)                                       +++++++++++++++++
!  ----------------------------------------------------------

        DO k=mz1_CM,mzp

! #wH       qr_EVP = 0.0

! #hy     IF                           (qr___0(ikl,k).gt.epsn)      THEN
            qr0_OK = max(zer0,sign(un_1,qr___0(ikl,k)  - epsn))
!           qr0_OK = 1.0 if             qr___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #EW      IF(qr0_OK.gt.eps6)                                       THEN
! #EW         mauxEW        =  mphyEW(ikl)
! #EW         mauxEW(15:15) = 'v'
! #EW         mphyEW(ikl)   =  mauxEW
! #EW      END IF

            RH_Liq = qv__DY(ikl,k)/(RHcrit*qvswCM(ikl,k))
!           RH_Liq : grid scale saturation humidity

! #hy     IF                                 (RH_Liq.lt.un_1)       THEN
            Flag_DryAir = max(zer0,-sign(un_1,RH_Liq  - un_1))
!           Flag_DryAir = 1.0 if              RH_Liq  < un_1
!                       = 0.0 otherwise

            Flag_qr_EVP = qr0_OK * Flag_DryAir

            lr_NUM = 0.78d0  /(lamdaR(ikl,k) *lamdaR(ikl,k))           &!
     &          + 3940.d0    *           sqrt(sqrrro(ikl,k))           &! 3940.: 0.31 Sc**(1/3) *(a/nu)**(1/2) * Gamma[(b+5)/2]
     &                       /exp(2.9d0  *log(lamdaR(ikl,k)))           ! where       Sc=0.8(Schm.) nu=1.5e-5 (Air Kinematic Viscosity)
            lr_DEN = 5.423d11/(Ta__CM(ikl,k) *Ta__CM(ikl,k))           &! 5.423e11 = [Lv=2500000J/kg] * Lv / [kT=0.025W/m/K] / [Rv=461.J/kg/K]
     &       + 1.d0/(1.875d-2* roa_DY(ikl,k) *qvswCM(ikl,k))            !                                     kT:  Air Thermal Conductivity

            pr_EVP = 2. *piNmbr*(1.d0 -RH_Liq)*n0___r*lr_NUM/lr_DEN
            qr_EVP =     pr_EVP*       dt__CM
            qr_EVP = min(qr_EVP,       qr__CM(ikl,k))

            qr_EVP = min(qr_EVP,RHcrit*qvswCM(ikl,k) -qv__DY(ikl,k))
!           supersaturation is not allowed to occur

            qr_EVP = max(qr_EVP,zer0)       *         Flag_qr_EVP
!           condensation    is not allowed to occur

            qr__CM(ikl,k) = qr__CM(ikl,k) -           qr_EVP
            qwd_CM(ikl,k) = qwd_CM(ikl,k) -           qr_EVP
            qv__DY(ikl,k) = qv__DY(ikl,k) +           qr_EVP
            Ta__CM(ikl,k) = Ta__CM(ikl,k) - Lv_Cpd   *qr_EVP

!  Full Debug
!  ~~~~~~~~~~
! #WQ       write(6,*) 'Qrevp',qr_EVP,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1)) wrevp(k) =    qr_EVP

! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Lin et al.(1983): Rain Drops Evapor'
! #wH         debugH(36:70)   = 'ation                              '
! #wH         proc_1          = 'Qrevp g/kg'
! #wH         procv1          =  qr_EVP
! #wH         proc_2          = 'R.Hum  [%]'
! #wH         procv2          =  RH_Liq*0.1
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,13)    = -qr_EVP

        END DO


!  (Deposition on) Snow Flakes (Sublimation)                                                        SNOW, (SUBLI)/(DEPOSI)TION
!   Reference: Lin et al.      1983, JCAM      22, p.1072 (31)                                      ++++++++++++++++++++++++++
!   ----------------------------------------------------------

        DO k=mz1_CM,mzp

! #wH       qs_SUB = 0.0

! #hy     IF                           (qs___0(ikl,k).gt.epsn)      THEN
            qs0_OK = max(zer0,sign(un_1,qs___0(ikl,k)  - epsn))
!           qs0_OK = 1.0 if             qs___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #EW      IF(qs0_OK.gt.eps6)                                       THEN
! #EW         mauxEW        =  mphyEW(ikl)
! #EW         mauxEW(16:16) = 'V'
! #EW         mphyEW(ikl)   =  mauxEW
! #EW      END IF

            RH_ICE =           qv__DY(ikl,k)/qsiEFF(ikl,k)

            ls_NUM = 0.78d0  /(lamdaS(ikl,k)*lamdaS(ikl,k))            &!
     &             + 238.d0  *          sqrt(sqrrro(ikl,k))            &! 238.: 0.31 Sc**(1/3) *(c/nu)**(1/2) * Gamma[(d+5)/2]
     &                      /exp(2.625d0*log(lamdaS(ikl,k)))            ! where      Sc=0.8(Schm.) nu=1.5e-5 (Air Kinematic Viscosity)
            ls_DEN = 6.959d11/(Ta__CM(ikl,k)*Ta__CM(ikl,k))            &! 6.959e11 = [Ls=2833600J/kg]*Ls /[kT=0.025W/m/K] /[Rv=461.J/kg/K]
     &        + 1.d0/(1.875d-2*roa_DY(ikl,k)*qsiEFF(ikl,k))             !                                  kT: Air Thermal   Conductivity

! #cn       n0___s = min(2.e8,2.e6*exp(-.12*min(0.,Ta_dgC(ikl,k))))

            ps_SUB = 2*piNmbr*(1.d0-RH_ICE)*n0___s*ls_NUM              &! 
     &                       /(1.d3*roa_DY(ikl,k)*ls_DEN)
            qs_SUB = ps_SUB  * dt__CM

            dqsiqv = qsiEFF(ikl,k) -qv__DY(ikl,k)

            Flag_SURSat = max(zer0,sign(un_1,RH_ICE - un_1))
!           Flag_SURSat = 1.0 if             RH_ICE > un_1
!                       = 0.0 otherwise

            qs_SUB = max(qs_SUB               ,dqsiqv)*    Flag_SURSat &! qs_SUB < 0 ... Deposition
     &         + min(min(qs_SUB,qs__CM(ikl,k)),dqsiqv)*(1.-Flag_SURSat) !        > 0 ... Sublimation

            qs_SUB =     qs_SUB * qs0_OK

            qs__CM(ikl,k) = qs__CM(ikl,k)-          qs_SUB
            qid_CM(ikl,k) = qid_CM(ikl,k)-          qs_SUB
            qv__DY(ikl,k) = qv__DY(ikl,k)+          qs_SUB
            Ta__CM(ikl,k) = Ta__CM(ikl,k)-Ls_Cpd   *qs_SUB

!  Full Debug
!  ~~~~~~~~~~
! #WQ       write(6,*) 'Qssub',qs_SUB,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1)) wssub(k) =   -qs_SUB

! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Lin et al.(1983): (Deposition on) S'
! #wH         debugH(36:70)   = 'now Particles (Sublimation)        '
! #wH         proc_1          = 'Qssub g/kg'
! #wH         procv1          =  qs_SUB
! #wH         proc_2          = '          '
! #wH         procv2          =  0.
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,14)    = -qs_SUB

        END DO


!  Graupels Sublimation                                                                             GRAUPEL, SUBLIMATION
!  Reference: Lin et al.      1983, JCAM      22, p.1076 (46)                                       ++++++++++++++++++++
!  ----------------------------------------------------------

! #qg   ! TO BE ADDED !


!  Snow Flakes Melting        PSMLT                                                                 SNOW, MELT
!  Reference: Lin et al.      1983, JCAM      22, p.1072 (32)                                       ++++++++++
!  ----------------------------------------------------------

        DO k=mz1_CM,mzp

! #wH       qsMELT = 0.0

! #hy     IF                           (qs___0(ikl,k).gt.epsn)      THEN
            qs0_OK = max(zer0,sign(un_1,qs___0(ikl,k)  - epsn))
!           qs0_OK = 1.0 if             qs___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                                (Ta_dgC(ikl,k).gt.0.)   THEN! Ta_dgC : old Celsius Temperature
            Flag_Ta_Pos = max(zer0,sign(un_1,Ta_dgC(ikl,k)  - 0.))
!           Flag_Ta_Pos = 1.0 if             Ta_dgC(ikl,k)  > 0.
!                       = 0.0 otherwise

            Flag_qsMELT = qs0_OK * Flag_Ta_Pos

! #EW      IF(Flag_qsMELT.gt.eps6)                                  THEN
! #EW         mauxEW        =  mphyEW(ikl)
! #EW         mauxEW(17:17) = 'r'
! #EW         mphyEW(ikl)   =  mauxEW
! #EW      END IF

            ls_NUM = 0.78 /  (lamdaS(ikl,k) *lamdaS(ikl,k))            &
     &           + 238.   *             sqrt(sqrrro(ikl,k))            &
     &                    / exp(2.625d0 *log(lamdaS(ikl,k)))

! #cn       n0___s = min(2.e8,2.e6*exp(-.12*min(0.,Ta_dgC(ikl,k))))

            xCoefM = 1.904d-8 *n0___s *ls_NUM *Lc_Cpd /roa_DY(ikl,k)    ! 1.904e-8: 2 pi / Lc /[1.e3=rho Factor]

            AcoefM = 0.025d00 *xCoefM                                  &!
     &       +(ps_ACW(ikl,k) + ps_ACR(ikl,k)) *Lc_Cpd /78.8d0           ! 78.8    :        Lc /[Cpw=4.187e3 J/kg/K]

            BcoefM = 62.34d+3 *roa_DY(ikl,k)                           &! 62.34   :        Ls *[psiv=2.200e-5 m2/s]
     &                       *(qv__DY(ikl,k)-qsiEFF(ikl,k))            &! 46.88   :        Lv *[psiv=1.875e-5 m2/s]
     &                        *xCoefM
            BcoefM = min(-epsn,BcoefM)

            dTMELT =    ( Ta__CM(ikl,k) -Tf_Sno -AcoefM/BcoefM)        &!
     &              *exp(-AcoefM*dt__CM)                                !
            qsMELT =    ( Ta__CM(ikl,k) -Tf_Sno -dTMELT       )/ Lc_Cpd !
            qsMELT = max( qsMELT,0.           ) *Flag_qsMELT            !
            qsMELT = min( qsMELT,qs__CM(ikl,k))                         !

            qs__CM(ikl,k) = qs__CM(ikl,k) -          qsMELT
            qr__CM(ikl,k) = qr__CM(ikl,k) +          qsMELT
            Ta__CM(ikl,k) = Ta__CM(ikl,k) - Lc_Cpd  *qsMELT

!  Full Debug
!  ~~~~~~~~~~
! #WQ       write(6,*) 'Qsmlt',qsMELT,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1)) wsmlt(k) =    qsMELT

! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Lin et al.(1983): Snow Particles Me'
! #wH         debugH(36:70)   = 'lting                              '
! #wH         proc_1          = 'Qsmlt g/kg'
! #wH         procv1          =  qsMELT
! #wH         proc_2          = '          '
! #wH         procv2          =  0.
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,15)   = -qsMELT

        END DO


!  Graupels    Melting                                                                              GRAUPELS, MELT
!  Reference: Lin et al.      1983, JCAM      22, p.1076 (47)                                       ++++++++++++++
!  ----------------------------------------------------------

! #qg   ! TO BE ADDED !


!  Rain Freezing                                                                                    RAIN, FREEZING
!  Reference: Lin et al.      1983, JCAM      22, p.1075 (45)                                       ++++++++++++++
!  ----------------------------------------------------------

!  **CAUTION**: Graupel Formation TO BE ADDED !

        DO k=1,mzp     

! #wH       qs_FRZ = 0.0

! #hy     IF                           (qr___0(ikl,k).gt.epsn)      THEN
            qr0_OK = max(zer0,sign(un_1,qr___0(ikl,k)  - epsn))
!           qr0_OK = 1.0 if             qr___0(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                                 (Ta_dgC(ikl,k).lt.0.e0)THEN! Ta_dgC : old Celsius Temperature
            Flag_Ta_Neg = max(zer0,-sign(un_1,Ta_dgC(ikl,k)  - 0.e0))
!           Flag_Ta_Neg = 1.0 if              Ta_dgC(ikl,k)  < 0.e0
!                       = 0.0 otherwise

            Flag_Freeze = qr0_OK * Flag_Ta_Neg

! #EW      IF(Flag_Freeze.gt.eps6)                                  THEN
! #EW         mauxEW        =  mphyEW(ikl)
! #EW         mauxEW(19:19) = 's'
! #EW         mphyEW(ikl)   =  mauxEW
! #EW      END IF

            ps_FRZ = 1.974d4 *n0___r                                   &
     &       /(roa_DY(ikl,k)*exp(7.d0 *log(lamdaR(ikl,k))))            &
     &                     *(exp(-0.66d0  *Ta_dgC(ikl,k))-1.d0)
            qs_FRZ =     ps_FRZ * dt__CM  *Flag_Freeze
            qs_FRZ = min(qs_FRZ,qr__CM(ikl,k))

            qr__CM(ikl,k) = qr__CM(ikl,k) -          qs_FRZ
            qs__CM(ikl,k) = qs__CM(ikl,k) +          qs_FRZ
!           CAUTION : graupel production is included into snow production
!                     proposed modification in line below.
! #qg       qg__CM(ikl,k) = qg__CM(ikl,k) +          qs_FRZ
            Ta__CM(ikl,k) = Ta__CM(ikl,k) + Lc_Cpd  *qs_FRZ

!  Full Debug
!  ~~~~~~~~~~
! #WQ       write(6,*) 'Qsfre',qs_FRZ,it_EXP,ikl,k
! #WH       if (ikl.eq.ikl0CM(1)) wsfre(kl) = qs_FRZ

! #hy     END IF
! #hy     END IF

!  Debug
!  ~~~~~
! #wH         debugH( 1:35)   = 'Lin et al.(1983): Rain Freezing    '
! #wH         debugH(36:70)   = '                                   '
! #wH         proc_1          = 'Qsfr g/kg'
! #wH         procv1          =  qs_FRZ
! #wH         proc_2          = '          '
! #wH         procv2          =  0.
! #wH         proc_3          = '          '
! #wH         procv3          =  0.
! #wH         proc_4          = '          '
! #wH         procv4          =  0.
! #wh         include 'CMiPhy_Debug.h'
! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))   &
! #wH&        debugV(k,16)   =  qs_FRZ

        END DO




!  Debug (Summary)
!  ===============

! #wH   DO k=mz1_CM,mzp 

! #wH     IF (ii__AP(ikl).EQ.i0__CM(1).AND.jj__AP(ikl).EQ.j0__CM(1))THEN
! #wH       IF (k .EQ.mz1_CM)                                       THEN
! #wH          write(6,6022)
 6022          format(/,'CMiPhy STATISTICS'                            &
     &                /,'=================')
! #wH          write(6,6026)
 6026          format(  '    T_Air Qv   Qw g/kg  Qi g/kg  CLOUDS % '   &
     &                 ,              ' Qs g/kg  Qr g/kg'              &
     &                 ,' Qi+ E.K.'                                    &
     &                 ,' Qi+ Mey.'                                    &
     &                 ,' Qi- Sub.'                                    &
     &                 ,' Qi- Mlt.'                                    &
     &                 ,' Qw+ Cds.'                                    &
     &                 ,' Qraut r+'                                    &
     &                 ,' QsAUT s+'                                    &
     &                 ,' Qracw r+')
! #wH       END IF
! #wH          write(6,6023)      k                                    &
! #wH&              ,      Ta__CM(ikl,k)-Tf_Sno                        &
! #wH&              ,1.e3* qv__DY(ikl,k)                               &
! #wH&              ,1.e3* qw__CM(ikl,k)                               &
! #wH&              ,1.e3* qi__CM(ikl,k)                               &
! #wH&              ,1.e2* CFraCM(ikl,k)                               &
! #wH&              ,1.e3* qs__CM(ikl,k)                               &
! #wH&              ,1.e3* qr__CM(ikl,k)                               &
! #wH&             ,(1.e3* debugV(k,kv),kv=1,08)
 6023          format(i3,f6.1,f5.2,2f9.6,f9.1,2f9.3,8f9.6)
! #wH       IF (k .EQ.mzp )                                         THEN
! #wH          write(6,6026)
! #wH          write(6,*)  ' '
! #wH          write(6,6024)
 6024          format(  8x,'Z [km]'                                    &
     &                 ,' RH.w.[%]'                                    &
     &                 ,' RH.i.[%]'     ,9x                            &
     &                 ,' Vss cm/s'                                    &
     &                 ,' Vrr cm/s'                                    &
     &                 ,' Qsacw s+'                                    &
     &                 ,' Qsaci s+'                                    &
     &                 ,' Qiacr r+'                                    &
     &                 ,' Qracs ds'                                    &
     &                 ,' Qrevp w-'                                    &
     &                 ,' Qssub s-'                                    &
     &                 ,' Qsmlt s-'                                    &
     &                 ,' Qsfr  s+')
! #wH         DO nl=mz1_CM,mzp 
! #wH          write(6,6025)   nl       ,zsigma(   nl)*1.e-3           &
! #wH&              ,1.e2*   qv__DY(ikl,nl)/qvswCM(ikl,nl)             &
! #wH&              ,1.e2*   qv__DY(ikl,nl)/qvsiCM(ikl,nl)             &
! #wH&              ,1.e2*   FallVs(ikl,nl)                            &
! #wH&              ,1.e2*   FallVr(ikl,nl)                            &
! #wH&             ,(1.e3*   debugV(nl,kv),kv=9,16)
 6025          format(i3,f11.3,    2f9.1,9x,  2f9.1,8f9.6)
! #wH         END DO
! #wH          write(6,6024)
! #wH          write(6,*)  ' '
! #wH       END IF

! #wH     END IF

! #wH   END DO




!  Vertical Integrated Energy and Water Content
!  ============================================

! #EW     enr1EW(ikl) = 0.0d00
! #EW     wat1EW(ikl) = 0.0d00

! #EW   DO k=1,mzp 
! #EW     enr1EW(ikl) = enr1EW(ikl )                                   &
! #EW&              +(Ta__CM(ikl,k)                                    &
! #EW&              -(qw__CM(ikl,k)+qr__CM(ikl,k)) *Lv_Cpd             &
! #EW&              -(qi__CM(ikl,k)+qs__CM(ikl,k)) *Ls_Cpd)*dsigmi(k)
! #EW     wat1EW(ikl) = wat1EW(ikl )                                   &
! #EW&              +(qv__DY(ikl,k)                                    &
! #EW&              + qw__CM(ikl,k)+qr__CM(ikl,k)                      &
! #EW&              + qi__CM(ikl,k)+qs__CM(ikl,k)         )*dsigmi(k)
! #EW   END DO

! #ew     enr1EW(ikl) = enr1EW(ikl ) * psa_DY(ikl) * Grav_I
! #EW     wat1EW(ikl) = wat1EW(ikl ) * psa_DY(ikl) * Grav_I
!  ..     wat1EW [m]   contains implicit factor 1.d3 [kPa-->Pa] /ro_Wat




!  Precipitation
!  =============

!  Hydrometeors Fall Velocity
!  --------------------------

!  Pristine Ice Crystals Diameter and Fall Velocity
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DO k=mz1_CM,mzp

! #hy     IF                           (qi__CM(ikl,k).gt.epsn)      THEN
            qi__OK = max(zer0,sign(un_1,qi__CM(ikl,k)  - epsn))
!           qi__OK = 1.0 if             qi__CM(ikl,k)  > epsn
!                  = 0.0 otherwise

! #hy     IF                           (CCNiCM(ikl,k).gt.1.e0)      THEN
            CCNiOK = max(zer0,sign(un_1,CCNiCM(ikl,k)  - 1.e0))
!           CCNiOK = 1.0 if             CCNiCM(ikl,k)  > 1.e0
!                  = 0.0 otherwise

            Flag_Fall_i = qi__OK  * CCNiOK

            Di_Pri = 0.16d0 *exp(R_1by3*log(R_1000*roa_DY(ikl,k)       &! Pristine Ice Crystals Diameter, where 6/(pi*ro_I)**1/3 ~ 0.16
     &         *max(epsn,qi__CM(ikl,k))/max(un_1  ,CCNiCM(ikl,k))))     ! REF.: Levkov et al. 1992, Contr. Atm. Phys. 65, (5) p.37

            FallVi(ikl,k) =      Flag_Fall_i * 7.d2*Di_Pri             &! Terminal Fall Velocity for Pristine Ice Crystals
     &         *exp( 0.35d0 *log(roa_DY(ikl,mzp)  / roa_DY(ikl,k)))     ! REF.: Levkov et al. 1992, Contr. Atm. Phys. 65, (4) p.37
! #hy     ELSE
! #hy       FallVi(ikl,k) =  0.0d00

! #hy     END IF
! #hy     END IF

        END DO

!  Set Up of the Numerical Scheme
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #EW     watfEW(ikl) = 0.d0                                            ! Water Flux (Atmosphere --> Surface)

!  Snow and Rain Fall Velocity (Correction)
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DO k=mz1_CM,mzp
          FallVi(ikl,k) = FallVi(ikl,k) *qi__CM(ikl,k)/max(qi__CM(ikl,k),epsn)
          FallVs(ikl,k) = FallVs(ikl,k) *qs__CM(ikl,k)/max(qs__CM(ikl,k),epsn)
! #VW     FallVw(ikl,k) = FallVw(ikl,k) *qw__CM(ikl,k)/max(qw__CM(ikl,k),epsn)
          FallVr(ikl,k) = FallVr(ikl,k) *qr__CM(ikl,k)/max(qr__CM(ikl,k),epsn)
        END DO


!  -----------------------------------------------------
!  Droplets              Precipitation (Implicit Scheme)
!  Pristine Ice Crystals Precipitation (Implicit Scheme)
!  Snow Particles        Precipitation (Implicit Scheme)
!  Rain Drops            Precipitation (Implicit Scheme)
!  -----------------------------------------------------

          qwLoss(mz1_CM-1) = 0.
          qiLoss(mz1_CM-1) = 0.
          qsLoss(mz1_CM-1) = 0.
          qrLoss(mz1_CM-1) = 0.

!  Precipitation Mass & Flux
!  ~~~~~~~~~~~~~~~~~~~~~~~~~
        DO k= mz1_CM,mzp
          qwLoss(k)       = 0.

          a_rodz          = Grav_I* psa_DY(  ikl) *dsigmi(k)            ! Air  Mass
! #VW     qwFlux          = dt__CM* FallVw(ikl,k) *roa_DY(ikl,k)        ! Flux Fact. (droplets)
          qiFlux          = dt__CM* FallVi(ikl,k) *roa_DY(ikl,k)        ! Flux Fact. (crystals)
          qsFlux          = dt__CM* FallVs(ikl,k) *roa_DY(ikl,k)        ! Flux Fact. (snow)
          qrFlux          = dt__CM* FallVr(ikl,k) *roa_DY(ikl,k)        ! Flux Fact. (rain)

! #VW     qwrodz          =         qw__CM(ikl,k) *a_rodz              &! Droplets Mass
! #VW&                    +    0.5 *qwLoss(k-1)                         ! From abov.
          qirodz          =         qi__CM(ikl,k) *a_rodz              &! Crystals Mass
     &                    +    0.5 *qiLoss(k-1)                         ! From abov.
          qsrodz          =         qs__CM(ikl,k) *a_rodz              &! Snow Mass
     &                    +    0.5 *qsLoss(k-1)                         ! From abov.
          qrrodz          =         qr__CM(ikl,k) *a_rodz              &! Rain Mass
     &                    +    0.5 *qrLoss(k-1)                         ! From abov.

! #VW     wRatio          =                                            &! Var. Fact.
! #VW&                       min(2.,qwFlux        /a_rodz       )       ! Flux Limi.
          iRatio          =                                            &! Var. Fact.
     &                       min(2.,qiFlux        /a_rodz       )       ! Flux Limi.
          sRatio          =                                            &! Var. Fact.
     &                       min(2.,qsFlux        /a_rodz       )       ! Flux Limi.
          rRatio          =                                            &! Var. Fact.
     &                       min(2.,qrFlux        /a_rodz       )       ! Flux Limi.

! #VW     qwLoss(k)       =         qwrodz        *wRatio              &! Mass Loss
! #VW&                         /(1.+wRatio        *0.5)                 !
          qiLoss(k)       =         qirodz        *iRatio              &! Mass Loss
     &                         /(1.+iRatio        *0.5)                 !
          qsLoss(k)       =         qsrodz        *sRatio              &! Mass Loss
     &                         /(1.+sRatio        *0.5)                 !
          qrLoss(k)       =         qrrodz        *rRatio              &! Mass Loss
     &                         /(1.+rRatio        *0.5)                 !

! #VW     qwrodz          =         qwrodz        -qwLoss(k)           &!
! #VW&                    +    0.5 *qwLoss(k-1)                         ! From abov.
          qirodz          =         qirodz        -qiLoss(k)           &!
     &                    +    0.5 *qiLoss(k-1)                         ! From abov.
          qsrodz          =         qsrodz        -qsLoss(k)           &!
     &                    +    0.5 *qsLoss(k-1)                         ! From abov.
          qrrodz          =         qrrodz        -qrLoss(k)           &!
     &                    +    0.5 *qrLoss(k-1)                         ! From abov.

!  Cooling from above precipitating flux
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Ta__CM(ikl,k)   =                                            &!
     &   (Ta__CM(ikl,k  ) * a_rodz                                     &!
     &   +Ta__CM(ikl,k-1) *(qwLoss(k-1)+qiLoss(k-1)+qsLoss(k-1)+qrLoss(k-1))) &!
     &  /(a_rodz           +qwLoss(k-1)+qiLoss(k-1)+qsLoss(k-1)+qrLoss(k-1))

! #VW     qw__CM(ikl,k)   = qwrodz                /a_rodz       
          qi__CM(ikl,k)   = qirodz                /a_rodz       
          qs__CM(ikl,k)   = qsrodz                /a_rodz       
          qr__CM(ikl,k)   = qrrodz                /a_rodz       
        ENDDO

!  Precipitation reaching the Surface
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          d_rain          = qrLoss(mzp) + qwLoss(mzp)                   ! d_rain contains implicit factor 1.e3[kPa->Pa]/ro_Wat[kg/m2->m w.e.]
          RainCM(ikl)     = RainCM(ikl) + d_rain                        ! RainCM:         rain precipitation height since start of run    [m]
          d_snow          =               qiLoss(mzp)                   ! d_snow contains implicit factor 1.e3[kPa->Pa]/ro_Wat[kg/m2->m w.e.]
          Ice_CM(ikl)     = Ice_CM(ikl) + d_snow                        ! Ice_CM:   ice        precipitation height since start of run    [m]
          d_snow          = qsLoss(mzp) + qiLoss(mzp)                   ! d_snow contains implicit factor 1.e3[kPa->Pa]/ro_Wat[kg/m2->m w.e.]
          SnowCM(ikl)     = SnowCM(ikl) + d_snow                        ! SnowCM:   ice + snow precipitation height since start of run    [m]

! #EW     watfEW(ikl )     = watfEW(ikl ) - d_rain - d_snow

          d_rain          = 0.0
          d_snow          = 0.0




!  Fractional  Cloudiness ! Guess may be computed (Ek&Mahrt91 fracSC=.T.)
!  ====================== ! Final value  computed  below

! #sc   IF (Frac__Clouds.AND..NOT.fracSC)                           THEN
        IF (Frac__Clouds)                                           THEN
         IF(fraCEP) THEN ! ECMWF Large Scale Cloudiness
                         ! ----------------------------
          DO k=mz1_CM,mzp
            CFraCM(ikl,k) =              (qi__CM(ikl,k) + qw__CM(ikl,k)&
     &                                   +qs__CM(ikl,k) * 0.33         &
     &               * (1.-min(1.,exp((Ta__CM(ikl,k) -258.15)*0.1))))  &
     &               / (0.02     *     qvswCM(ikl,k)                )
            CFraCM(ikl,k) =min(1.000 , CFraCM(ikl,k))
            CFraCM(ikl,k) =max(0.001 , CFraCM(ikl,k))                  &
     &                *max(zer0,sign(un_1,qi__CM(ikl,k) + qw__CM(ikl,k)&
     &                                   +qs__CM(ikl,k) -3.E-9       ))
          END DO
         ELSE            ! XU and Randall  1996, JAS 21, p.3099 (4)
                         ! ----------------------------
          DO k=mz1_CM,mzp
            qvs_wi=                                        qvswCM(ikl,k)
! #wi       qvs_wi=max(epsn,((qi__CM(ikl,k)+qs__CM(ikl,k))*qvsiCM(ikl,k)    &
! #wi&                       +qw__CM(ikl,k)               *qvswCM(ikl,k))   &
! #wi&              /max(epsn,qi__CM(ikl,k)+qs__CM(ikl,k) +qw__CM(ikl,k)))
            RHumid=min(RH_MAX,max(qv__DY(ikl,k),qv_MIN) / qvs_wi)
            argEXP=  ((RH_MAX                  -RHumid) * qvs_wi)**0.49
            argEXP=min(100.     *(qi__CM(ikl,k)+qw__CM(ikl,k)          &
     &                                         +qs__CM(ikl,k) *  0.33  &
     &          * (1.-min(1.,exp((Ta__CM(ikl,k)-258.15)*0.1))))        &
     &                       /max(epsn         ,argEXP),ea_MAX)
              
            CFraCM(ikl,k) = ( RHumid ** 0.25 ) * ( 1. - exp(-argEXP) )
          END DO
         END IF

        ELSE
! #sc   ELSE IF (.NOT.Frac__Clouds)                                 THEN
          DO k=mz1_CM,mzp
              qCloud        =     qi__CM(ikl,k)  + qw__CM(ikl,k)
! #hy       IF                                    (qCloud.gt.epsn)  THEN
              CFraCM(ikl,k) = max(zer0,sign(un_1,  qCloud  - epsn))
!             CFraCM(ikl,k) = 1.0 if               qCloud  > epsn
!                           = 0.0 otherwise

! #hy       END IF
          END DO

        END IF




!  Vertically Integrated Energy and Water Content
!  ==============================================

! #EW     enr2EW(ikl) = 0.0d00
! #EW     wat2EW(ikl) = 0.0d00
!  ..     Vertical Integrated Energy and Water Content

! #EW   DO k=1,mzp
! #EW     enr2EW(ikl) = enr2EW(ikl)                                    &
! #EW&            +  (Ta__CM(ikl,k)                                    &
! #EW&            -  (qw__CM(ikl,k)+qr__CM(ikl,k))*Lv_Cpd              &
! #EW&            -  (qi__CM(ikl,k)+qs__CM(ikl,k))*Ls_Cpd) *dsigmi(k)
! #EW     wat2EW(ikl) = wat2EW(ikl)                                    &
! #EW&            +  (qv__DY(ikl,k)                                    &
! #EW&            +   qw__CM(ikl,k)+qr__CM(ikl,k)                      &
! #EW&            +   qi__CM(ikl,k)+qs__CM(ikl,k)        ) *dsigmi(k)
! #EW   END DO

! #ew     enr2EW(ikl) = enr2EW(ikl) * psa_DY(ikl) * Grav_I
! #EW     wat2EW(ikl) = wat2EW(ikl) * psa_DY(ikl) * Grav_I
!  ..     wat2EW [m]   contains implicit factor 1.d3 [kPa-->Pa] /ro_Wat




!  Limits on Microphysical Variables
!  =================================

        DO k=1,mz1_CM
            qv__DY(ikl,k)=max(qv__DY(ikl,k),qv_MIN)
            qv__DY(ikl,k)=min(qv__DY(ikl,k),qvsiCM(ikl,k))
            qw__CM(ikl,k)=    zer0
            qi__CM(ikl,k)=    zer0
            CCNiCM(ikl,k)=    zer0
            qr__CM(ikl,k)=    zer0
            qs__CM(ikl,k)=    zer0
        END DO

        DO k=mz1_CM,mzp
            qw__CM(ikl,k)=max(zer0,  qw__CM(ikl,k))
            qi__CM(ikl,k)=max(zer0,  qi__CM(ikl,k))
            CCNiCM(ikl,k)=max(zer0,  CCNiCM(ikl,k))
            qr__CM(ikl,k)=max(zer0,  qr__CM(ikl,k))
            qs__CM(ikl,k)=max(zer0,  qs__CM(ikl,k))
        END DO




!     +++++++++++++++++++
      ENDDO ! ikl=1,kcolp
!     +++++++++++++++++++




!  OUTPUT
!  ======

! #WH IF (mod(MinuTU,6).eq.0.and.Sec_TU.eq.0.and.ikl0CM(1).gt.0)       THEN
! #WH   write(6,1030) HourTU,MinuTU,Sec_TU,it_EXP,i0__CM(1),j0__CM(1)
 1030   format(//,i4,'UT',i2,'m',i2,'s (iter.',i6,')  /  Pt.(',2i4,')' &
     &         ,/,'  ==========================================')
! #WH   write(6,1031)(k,     Z___DY(ikl0CM(1),k),qv__DY(ikl0CM(1),k),  &
! #WH&   1.d3*qi_io0(k),                    1.d3*qi__CM(ikl0CM(1),k),  &
! #WH&   1.d3* wihm1(k),1.d3* wihm2(k),     1.d3* wicnd(k),            &
! #WH&   1.d3* widep(k),1.d3* wisub(k),     1.d3* wimlt(k),k=mz1_CM,mzp)
 1031   format(/,                                                      &
     &     '            |  Water Vapor |  Cloud Ice, Time n & n+1',    &
     &     '   Cloud Ice Nucleation Processes    |',                   &
     &     '   Bergeron   Sublimation   Melting  ',                    &
     &   /,'  k    z[m] |  qv   [g/kg] |  qi_n [g/kg] qi_n+[g/kg]',    &
     &     ' QiHm1[g/kg] QiHm2[g/kg] QiCnd[g/kg] |',                   &
     &     '  QiDep[g/kg] QiSub[g/kg] QiMlt[q/kg]',                    &
     &   /,'------------+--------------+-------------------------',    &
     &     '-------------------------------------+',                   &
     &     '-------------------------------------',                    &
     &   /,(i3,f8.1,' | ',f12.6,' | ',2f12.6,3d12.4,' | ',3d12.4))

! #WH   write(6,1032)(k,Z___DY(ikl0CM(1),k)                            &
! #WH&            ,1.d3*qs___0(ikl0CM(1),k),1.d3*qs__CM(ikl0CM(1),k)   &
! #WH&            ,1.d3* wsaut(k),          1.d3* wsaci(k)             &
! #WH&            ,1.d3* wsacw(k),          1.d3* wiacr(k)             &
! #WH&            ,1.d3* wsacr(k),          1.d3* wssub(k)             &
! #WH&            ,     FallVs(k,ikl0CM(1)),k=mz1_CM,mzp)
 1032   format(/,                                                      &
     &     '            |  Snow Flakes, Time n&n+1 Autoconver. |',     &
     &     '  Accretion Processes ===> Snow Flakes            |',      &
     &     '  Sublimation | Term.F.Vel',                               &
     &   /,'  k    z[m] |  qs_n [g/kg] qs_n+[g/kg] QsAUT[g/kg] |',     &
     &     '  Qsaci[g/kg] Qsacw[g/kg] Qiacr[g/kg] Qsacr[g/kg] |',      &
     &     '  QsSub[g/kg] | vs   [m/s]',                               &
     &   /,'------------+--------------------------------------+',     &
     &     '--------------------------------------------------+',      &
     &     '--------------+-----------',                               &
     &   /,(i3,f8.1,' | ',2f12.6,e12.4,' | ',4d12.4,' | ',e12.4,       &
     &              ' | ',f10.6))

! #WH   write(6,1033)(k,Z___DY(ikl0CM(1),k),Ta__CM(ikl0CM(1),k)
! #WH&            ,1.d3*qw_io0(k)  ,1.d3*qw__CM(ikl0CM(1),k)
! #WH&            ,1.d3* wwevp(k)  ,1.d2*CFraCM(ikl0CM(1),k),k=mz1_CM,mzp)
 1033   format(/,                                                      &
     &   /,'            | Temperat.|  Cloud Water, Time n&n+1',        &
     &     ' Condens/Evp | Cloud ',                                    &
     &   /,'  k    z[m] | T    [K] |  qw_n [g/kg] qw_n+[g/kg]',        &
     &     ' QwEvp[g/kg] | Fract.',                                    &
     &   /,'------------+----------+-------------------------',        &
     &     '-------------+-------',                                    &
     &   /,(i3,f8.1,' | ',f8.3,' | ',2f12.6,e12.4,' | ',f5.1))

! #WH   write(6,1034)(k,Z___DY(ikl0CM(1),k),                           &
! #WH&            ,1.d3*qr___0(k,ikl0CM(1)),1.d3*qr__CM(ikl0CM(1),k),  &
! #WH&            ,1.d3* wraut(k)          ,1.d3* wracw(k)             &
! #WH&            ,1.d3* wraci(k)          ,1.d3* wracs(k)             &
! #WH&            ,1.d3* wrevp(k)          ,1.d3* wsfre(k)             &
! #WH&            ,     FallVr(k,ikl0CM(1)),               k=mz1_CM,mzp)
 1034   format(/,                                                      &
     &  /,'            | Rain Drops, Time n&n+1   Autoconver. |',      &
     &    '  Accretion Processes ===> Rain Drops |',                   &
     &    '  Evaporation  Freezing   | Term.F.Vel',                    &
     &  /,'  k    z[m] |  qr_n [g/kg] qr_n+[g/kg] Qraut[g/kg] |',      &
     &    '  Qracw[g/kg] Qraci[g/kg] Qracs[g/kg] |',                   &
     &    '  QrEvp[g/kg] QsFre[g/kg] | vr   [m/s]',                    &
     &  /,'------------+--------------------------------------+',      &
     &    '--------------------------------------+',                   &
     &    '--------------------------+-----------',                    &
     &  /,(i3,f8.1,' | ',2f12.6,e12.4,' | ',3d12.4,' | ',2d12.4,       &
     &             ' | ',f10.6))

! #WH   DO k=mz1_CM,mzp
! #WH     wihm1(k) = 0.d0
! #WH     wihm2(k) = 0.d0
! #WH     wicnd(k) = 0.d0
! #WH     widep(k) = 0.d0
! #WH     wisub(k) = 0.d0
! #WH     wimlt(k) = 0.d0
! #WH     wwevp(k) = 0.d0
! #WH     wraut(k) = 0.d0
! #WH     wsaut(k) = 0.d0
! #WH     wracw(k) = 0.d0
! #WH     wsacw(k) = 0.d0
! #WH     wsaci(k) = 0.d0
! #WH     wraci(k) = 0.d0
! #WH     wiacr(k) = 0.d0
! #WH     wsacr(k) = 0.d0
! #WH     wracs(k) = 0.d0
! #WH     wrevp(k) = 0.d0
! #WH     wssub(k) = 0.d0
! #WH     wsmlt(k) = 0.d0
! #WH     wsfre(k) = 0.d0
! #WH   END DO
! #WH END IF


!  Vertical Integrated Energy and Water Content: OUTPUT
!  ====================================================

! #EW IF (ikl0CM(1).gt.0)                                              THEN
! #EW   WaterB = wat2EW(ikl0CM(1))-wat1EW(ikl0CM(1))-watfEW(ikl0CM(1))
! #EW   write(6,606) it_EXP,                                           &
! #EW&                     enr0EW(ikl0CM(1)),1.d3*wat0EW(ikl0CM(1)),   &
! #EW&                     mphyEW(ikl0CM(1)),                          &
! #EW&                     enr1EW(ikl0CM(1)),1.d3*wat1EW(ikl0CM(1)),   &
! #EW&                     enr2EW(ikl0CM(1)),1.d3*wat2EW(ikl0CM(1)),   &
! #EW&                                       1.d3*watfEW(ikl0CM(1)),   &
! #EW&                                       1.d3*WaterB
 606    format(i9,'  Before mPhy:  E0 =',f12.6,'  W0 = ',f9.6,3x,a20   &
     &   ,3x,/,9x,'  Before Prec:  E1 =',f12.6,'  W1 = ',f9.6          &
     &   ,   /,9x,'  After  Prec:  E2 =',f12.6,'  W2 = ',f9.6          &
     &   ,                                     '  W Flux =',f9.6       &
     &   ,                                     '  Div(W) =',e9.3)
! #EW END IF

      IF (MinuTU.eq.0.and.Sec_TU.eq.0.and.mod(HourTU,3).eq.0)       THEN
        DO  ipt_CM=1,npt_CM  
            ikl=ikl0CM(ipt_CM)

              write(4,1037)               HourTU,MinuTU,               &
     &                      i0__CM(ipt_CM),j0__CM(ipt_CM)
 1037         format(/,' Ice-Crystal mPhy ',                           &
     &                   2x,' ',2x,1x,i2,'h',i2,'UT',                  &
     &                 ' -- Grid Point (',i5,',',i5,')',               &
     &  /,' =========================================================='&
     &   ,     /,'     |  z  [m] | T  [K] | qi[g/kg] |'                &
     &   ,            ' Ni [m-3] | Ni0[m-3] | vi [m/s] | qs[g/kg] |'   &
     &   ,     /,'-----+---------+--------+----------+'                &
     &   ,            '----------+----------+----------+----------+')
              write(4,1038)(k,Z___DY(ikl,k),Ta__CM(ikl,k)              &
     &                       ,qi__CM(ikl,k)*1.d3                       &
     &                       ,CCNiCM(ikl,k),Fletch(ikl,k)              &
     &                       ,FallVi(ikl,k),qs__CM(ikl,k)*1.d3         &
     &                     ,k=mz1_CM,mzp)
 1038         format((i4,' |' ,  f8.1,' |',f7.2,' |',f9.6,' |',        &
     &            2(d9.3,' |'),2(f9.6,' |')))

        END DO
      END IF




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! DE-ALLOCATION                                                         !
! =============                                                         !
!                                                                       !
      IF (FlagDALLOC)                                              THEN !

          deallocate ( qw_io0 )                                         ! Droplets   Concentration entering CMiPhy                  [kg/kg]
          deallocate ( qi_io0 )                                         ! Ice  Part. Concentration entering CMiPhy                  [kg/kg]
          deallocate ( qs___0 )                                         ! Snow Part. Concentration entering CMiPhy                  [kg/kg]
! #qg     deallocate ( qg___0 )                                         ! Graupels   Concentration entering CMiPhy                  [kg/kg]
          deallocate ( qr___0 )                                         ! Rain Drops Concentration entering CMiPhy                  [kg/kg]
          deallocate ( Ta_dgC )                                         ! Air   Temperature                                           [dgC]
          deallocate ( sqrrro )                                         ! sqrt(roa(mzp)/roa(k))                                         [-]
          deallocate ( qsiEFF )                                         ! EFFective Saturation Specific Humidity over Ice           [kg/kg]
          deallocate ( Fletch )                                         ! Monodisperse Nb of hexagonal Plates, Fletcher (1962)          [-]
          deallocate ( lamdaS )                                         ! Marshall-Palmer distribution parameter for Snow Particl.
! #qg     deallocate ( lamdaG )                                         ! Marshall-Palmer distribution parameter for Graupels
          deallocate ( lamdaR )                                         ! Marshall-Palmer distribution parameter for Rain Drops
          deallocate ( ps_ACR )                                         ! Accretion of Snow        by Rain                   Rate [kg/kg/s]
          deallocate ( ps_ACW )                                         ! Accretion of Cloud Drop. by Snow Particl.          Rate [kg/kg/s]
          deallocate ( FallVw )                                         ! Sedimentation Velocity   of Droplets
          deallocate ( FallVi )                                         ! Sedimentation Velocity   of Ice  Particles
          deallocate ( FallVs )                                         ! Sedimentation Velocity   of Snow Particles
! #qg     deallocate ( FallVg )                                         ! Sedimentation Velocity   of Snow Particles
          deallocate ( FallVr )                                         ! Sedimentation Velocity   of Rain Drops 
          deallocate ( qwLoss )                                         ! Mass Loss related to Sedimentation of Rain Droplets
          deallocate ( qiLoss )                                         ! Mass Loss related to Sedimentation of Ice  Crystals
          deallocate ( qsLoss )                                         ! Mass Loss related to Sedimentation of Snow Particles
          deallocate ( qrLoss )                                         ! Mass Loss related to Sedimentation of Rain Drops
! #wH     deallocate ( debugV )                                         ! Debug Variable (of 16 microphysical processes)

! #WH     deallocate ( wihm1 )                                          ! Cloud Droplets Freezing
! #WH     deallocate ( wihm2 )                                          ! Ice   Crystals Homogeneous Sublimation
! #WH     deallocate ( wicnd )                                          ! Ice   Crystals Nucleation              (Emde & Kahlig)
! #WH     deallocate ( widep )                                          ! Ice   Crystals Growth Bergeron Process (Emde & Kahlig)
! #WH     deallocate ( wisub )                                          ! Ice   Crystals             Sublimation (Levkov)
! #WH     deallocate ( wimlt )                                          ! Ice   Crystals Melting  
! #WH     deallocate ( wwevp )                                          ! Water Vapor Condensation / Evaporation (Fractional Cloudiness)
! #WH     deallocate ( wraut )                                          ! Cloud Droplets AUTO-Conversion
! #WH     deallocate ( wsaut )                                          ! Ice   Crystals AUTO-Conversion
! #WH     deallocate ( wracw )                                          ! Accretion of Cloud Droplets by Rain, Ta > 0, --> Rain
! #WH     deallocate ( wsacw )                                          ! Accretion of Cloud Droplets by Rain, Ta < 0, --> Snow
! #WH     deallocate ( wsaci )                                          ! Accretion of Ice   Crystals by Snow          --> Snow
! #WH     deallocate ( wraci )                                          ! Accretion of Ice   Crystals by Rain          --> Snow
! #WH     deallocate ( wiacr )                                          ! Accretion of Rain by Ice   Crystals          --> Snow
! #WH     deallocate ( wsacr )                                          ! Accretion of Rain by Snow                    --> Snow
! #WH     deallocate ( wracs )                                          ! Accretion of Snow by Rain                    --> Snow, Rain
! #WH     deallocate ( wrevp )                                          ! Rain  Drops     Evaporation  
! #WH     deallocate ( wssub )                                          ! Snow  Particles Sublimation
! #WH     deallocate ( wsmlt )                                          ! Snow  Particles Melting 
! #WH     deallocate ( wsfre )                                          ! Rain  Drops     Freezing
!                                                                       !
      END IF                                                            !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      return
      end subroutine CMiPhy
