MODULE wx_pbl_mod
!
! Planetary Boundary Layer and Surface module
!
! This module manage the calculation of turbulent diffusion in the boundary layer 
! and all interactions towards the differents sub-surfaces.
!
!
  USE dimphy

  IMPLICIT NONE

  REAL, ALLOCATABLE, DIMENSION(:), SAVE        :: Kech_Tp, Kech_T_xp, Kech_T_wp
  REAL, ALLOCATABLE, DIMENSION(:), SAVE        :: dd_KTp, KxKwTp, dd_AT, dd_BT
!$OMP THREADPRIVATE(Kech_Tp, Kech_T_xp, Kech_T_wp, dd_KTp, KxKwTp, dd_AT, dd_BT)
  REAL, ALLOCATABLE, DIMENSION(:), SAVE        :: Kech_Qp, Kech_Q_xp, Kech_Q_wp
  REAL, ALLOCATABLE, DIMENSION(:), SAVE        :: dd_KQp, KxKwQp, dd_AQ, dd_BQ
!$OMP THREADPRIVATE(Kech_Qp, Kech_Q_xp, Kech_Q_wp, dd_KQp, KxKwQp, dd_AQ, dd_BQ)
  REAL, ALLOCATABLE, DIMENSION(:), SAVE        :: Kech_Up, Kech_U_xp, Kech_U_wp
  REAL, ALLOCATABLE, DIMENSION(:), SAVE        :: dd_KUp, KxKwUp, dd_AU, dd_BU
!$OMP THREADPRIVATE(Kech_Up, Kech_U_xp, Kech_U_wp, dd_KUp, KxKwUp, dd_AU, dd_BU)
  REAL, ALLOCATABLE, DIMENSION(:), SAVE        :: Kech_Vp, Kech_V_xp, Kech_V_wp
  REAL, ALLOCATABLE, DIMENSION(:), SAVE        :: dd_KVp, KxKwVp, dd_AV, dd_BV
!$OMP THREADPRIVATE(Kech_Vp, Kech_V_xp, Kech_V_wp, dd_KVp, KxKwVp, dd_AV, dd_BV)

CONTAINS
!
!****************************************************************************************
!
SUBROUTINE wx_pbl_init

! Local variables
!****************************************************************************************
    INTEGER                       :: ierr
 

!****************************************************************************************
! Allocate module variables
!
!****************************************************************************************    

    ierr = 0

    ALLOCATE(Kech_Tp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(Kech_T_xp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(Kech_T_wp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_KTp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(KxKwTp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_AT(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_BT(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

!----------------------------------------------------------------------------
    ALLOCATE(Kech_Qp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(Kech_Q_xp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(Kech_Q_wp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_KQp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(KxKwQp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_AQ(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_BQ(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

!----------------------------------------------------------------------------
    ALLOCATE(Kech_Up(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(Kech_U_xp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(Kech_U_wp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_KUp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(KxKwUp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_AU(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_BU(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

!----------------------------------------------------------------------------
    ALLOCATE(Kech_Vp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(Kech_V_xp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(Kech_V_wp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_KVp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(KxKwVp(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_AV(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

    ALLOCATE(dd_BV(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('wx_pbl_init', 'pb in allocation',1)

!----------------------------------------------------------------------------

END SUBROUTINE wx_pbl_init

SUBROUTINE wx_pbl0_fuse(knon, dtime, ypplay, ywake_s, &
                                 yt_x, yt_w, yq_x, yq_w, &
                                 yu_x, yu_w, yv_x, yv_w, &
                                 ycdragh_x, ycdragh_w, ycdragm_x, ycdragm_w, &
                                 AcoefT_x, AcoefT_w, AcoefQ_x, AcoefQ_w, &
                                 AcoefU_x, AcoefU_w, AcoefV_x, AcoefV_w, &
                                 BcoefT_x, BcoefT_w, BcoefQ_x, BcoefQ_w, &
                                 BcoefU_x, BcoefU_w, BcoefV_x, BcoefV_w, &
                                 AcoefT, AcoefQ, AcoefU, AcoefV, &
                                 BcoefT, BcoefQ, BcoefU, BcoefV, &
                                 ycdragh, ycdragm, &
                                 yt1, yq1, yu1, yv1 &
                                 )
!
    USE print_control_mod, ONLY: prt_level,lunout
!
    INCLUDE "YOMCST.h"
!
    INTEGER,                      INTENT(IN)        :: knon    ! number of grid cells
    REAL,                         INTENT(IN)        :: dtime   ! time step size (s)
    REAL, DIMENSION(knon,klev),   INTENT(IN)        :: ypplay  ! mid-layer pressure (Pa)
    REAL, DIMENSION(knon),        INTENT(IN)        :: ywake_s ! cold pools fractional area
    REAL, DIMENSION(knon,klev),   INTENT(IN)        :: yt_x, yt_w, yq_x, yq_w
    REAL, DIMENSION(knon,klev),   INTENT(IN)        :: yu_x, yu_w, yv_x, yv_w
    REAL, DIMENSION(knon),        INTENT(IN)        :: ycdragh_x, ycdragh_w, ycdragm_x, ycdragm_w
    REAL, DIMENSION(knon),        INTENT(IN)        :: AcoefT_x, AcoefT_w, AcoefQ_x, AcoefQ_w
    REAL, DIMENSION(knon),        INTENT(IN)        :: AcoefU_x, AcoefU_w, AcoefV_x, AcoefV_w
    REAL, DIMENSION(knon),        INTENT(IN)        :: BcoefT_x, BcoefT_w, BcoefQ_x, BcoefQ_w
    REAL, DIMENSION(knon),        INTENT(IN)        :: BcoefU_x, BcoefU_w, BcoefV_x, BcoefV_w
    REAL, DIMENSION(knon),        INTENT(OUT)       :: AcoefT, AcoefQ, AcoefU, AcoefV
    REAL, DIMENSION(knon),        INTENT(OUT)       :: BcoefT, BcoefQ, BcoefU, BcoefV
    REAL, DIMENSION(knon),        INTENT(OUT)       :: ycdragh, ycdragm
    REAL, DIMENSION(knon),        INTENT(OUT)       :: yt1, yq1, yu1, yv1  ! Apparent T, q, u, v at first level, as
                                                                           !seen by surface modules
!
! Local variables
    INTEGER                    :: j
    REAL                       :: rho1
    REAL                       :: mod_wind_x
    REAL                       :: mod_wind_w   
    REAL                       :: dd_Cdragh
    REAL                       :: dd_Cdragm
    REAL                       :: dd_Kh
    REAL                       :: dd_Km
    REAL                       :: dd_u 
    REAL                       :: dd_v 
    REAL                       :: dd_t 
    REAL                       :: dd_q 
!
    REAL                       :: KCT, KCQ, KCU, KCV
!
    REAL                       :: BBT, BBQ, BBU, BBV
    REAL                       :: DDT, DDQ, DDU, DDV
    REAL                       :: LambdaT, LambdaQ, LambdaU, LambdaV 
    REAL                       :: LambdaTs, LambdaQs, LambdaUs, LambdaVs 
!
    REAL, DIMENSION(knon)      :: sigx       ! fractional area of (x) region

    REAL, DIMENSION(knon)      :: Kech_h    ! Energy exchange coefficient
    REAL, DIMENSION(knon)      :: Kech_h_x, Kech_h_w
    REAL, DIMENSION(knon)      :: Kech_m    ! Momentum exchange coefficient
    REAL, DIMENSION(knon)      :: Kech_m_x, Kech_m_w

!!!
!!! jyg le 09/04/2013 ; passage aux nouvelles expressions en differences

        sigx(:) = 1.-ywake_s(:)

        DO j=1,knon
!
! Calcul des coefficients d echange
         mod_wind_x = 1.0+SQRT(yu_x(j,1)**2+yv_x(j,1)**2)
         mod_wind_w = 1.0+SQRT(yu_w(j,1)**2+yv_w(j,1)**2)
!!         rho1 = ypplay(j,1)/(RD*yt(j,1))
         rho1 = ypplay(j,1)/(RD*(yt_x(j,1) + ywake_s(j)*(yt_w(j,1)-yt_x(j,1))))
         Kech_h_x(j) = ycdragh_x(j) * mod_wind_x * rho1
         Kech_h_w(j) = ycdragh_w(j) * mod_wind_w * rho1
         Kech_m_x(j) = ycdragm_x(j) * mod_wind_x * rho1
         Kech_m_w(j) = ycdragm_w(j) * mod_wind_w * rho1
!
         dd_Kh = Kech_h_w(j) - Kech_h_x(j)
         dd_Km = Kech_m_w(j) - Kech_m_x(j)
         IF (prt_level >=10) THEN
          print *,' mod_wind_x, mod_wind_w ', mod_wind_x, mod_wind_w
          print *,' rho1 ',rho1
          print *,' ycdragh_x(j),ycdragm_x(j) ',ycdragh_x(j),ycdragm_x(j)
          print *,' ycdragh_w(j),ycdragm_w(j) ',ycdragh_w(j),ycdragm_w(j)
          print *,' dd_Kh: ',dd_Kh
         ENDIF
!
         Kech_h(j) = Kech_h_x(j) + ywake_s(j)*dd_Kh
         Kech_m(j) = Kech_m_x(j) + ywake_s(j)*dd_Km
!
! Calcul des coefficients d echange corriges des retroactions
        Kech_T_xp(j) = Kech_h_x(j)/(1.-BcoefT_x(j)*Kech_h_x(j)*dtime)
        Kech_T_wp(j) = Kech_h_w(j)/(1.-BcoefT_w(j)*Kech_h_w(j)*dtime)
        Kech_Q_xp(j) = Kech_h_x(j)/(1.-BcoefQ_x(j)*Kech_h_x(j)*dtime)
        Kech_Q_wp(j) = Kech_h_w(j)/(1.-BcoefQ_w(j)*Kech_h_w(j)*dtime)
        Kech_U_xp(j) = Kech_m_x(j)/(1.-BcoefU_x(j)*Kech_m_x(j)*dtime)
        Kech_U_wp(j) = Kech_m_w(j)/(1.-BcoefU_w(j)*Kech_m_w(j)*dtime)
        Kech_V_xp(j) = Kech_m_x(j)/(1.-BcoefV_x(j)*Kech_m_x(j)*dtime)
        Kech_V_wp(j) = Kech_m_w(j)/(1.-BcoefV_w(j)*Kech_m_w(j)*dtime)
!
         dd_KTp(j) = Kech_T_wp(j) - Kech_T_xp(j)
         dd_KQp(j) = Kech_Q_wp(j) - Kech_Q_xp(j)
         dd_KUp(j) = Kech_U_wp(j) - Kech_U_xp(j)
         dd_KVp(j) = Kech_V_wp(j) - Kech_V_xp(j)
!
        Kech_Tp(j) = Kech_T_xp(j) + ywake_s(j)*dd_KTp(j)
        Kech_Qp(j) = Kech_Q_xp(j) + ywake_s(j)*dd_KQp(j)
        Kech_Up(j) = Kech_U_xp(j) + ywake_s(j)*dd_KUp(j)
        Kech_Vp(j) = Kech_V_xp(j) + ywake_s(j)*dd_KVp(j)
!
! Calcul des differences w-x
       dd_Cdragm = ycdragm_w(j) - ycdragm_x(j)
       dd_Cdragh = ycdragh_w(j) - ycdragh_x(j)
       dd_u = yu_w(j,1) - yu_x(j,1)
       dd_v = yv_w(j,1) - yv_x(j,1)
       dd_t = yt_w(j,1) - yt_x(j,1)
       dd_q = yq_w(j,1) - yq_x(j,1)
       dd_AT(j) = AcoefT_w(j) - AcoefT_x(j)
       dd_AQ(j) = AcoefQ_w(j) - AcoefQ_x(j)
       dd_AU(j) = AcoefU_w(j) - AcoefU_x(j)
       dd_AV(j) = AcoefV_w(j) - AcoefV_x(j)
       dd_BT(j) = BcoefT_w(j) - BcoefT_x(j)
       dd_BQ(j) = BcoefQ_w(j) - BcoefQ_x(j)
       dd_BU(j) = BcoefU_w(j) - BcoefU_x(j)
       dd_BV(j) = BcoefV_w(j) - BcoefV_x(j)
!
       KxKwTp(j) = Kech_T_xp(j)*Kech_T_wp(j)
       KxKwQp(j) = Kech_Q_xp(j)*Kech_Q_wp(j)
       KxKwUp(j) = Kech_U_xp(j)*Kech_U_wp(j)
       KxKwVp(j) = Kech_V_xp(j)*Kech_V_wp(j)
       BBT = (BcoefT_x(j) + sigx(j)*dd_BT(j))*dtime
       BBQ = (BcoefQ_x(j) + sigx(j)*dd_BQ(j))*dtime
       BBU = (BcoefU_x(j) + sigx(j)*dd_BU(j))*dtime
       BBV = (BcoefV_x(j) + sigx(j)*dd_BV(j))*dtime
       KCT = Kech_h(j)
       KCQ = Kech_h(j)
       KCU = Kech_m(j)
       KCV = Kech_m(j)
       DDT = Kech_Tp(j)
       DDQ = Kech_Qp(j)
       DDU = Kech_Up(j)
       DDV = Kech_Vp(j)
       LambdaT = dd_Kh/KCT
       LambdaQ = dd_Kh/KCQ
       LambdaU = dd_Km/KCU
       LambdaV = dd_Km/KCV
       LambdaTs = dd_KTp(j)/DDT
       LambdaQs = dd_KQp(j)/DDQ
       LambdaUs = dd_KUp(j)/DDU
       LambdaVs = dd_KVp(j)/DDV
!
       IF (prt_level >=10) THEN
          print *,'Variables pour la fusion : Kech_T_xp(j)' ,Kech_T_xp(j)
          print *,'Variables pour la fusion : Kech_T_wp(j)' ,Kech_T_wp(j)
          print *,'Variables pour la fusion : Kech_Tp(j)' ,Kech_Tp(j)
          print *,'Variables pour la fusion : Kech_h(j)' ,Kech_h(j)
       ENDIF
!
! Calcul des coef A, B \'equivalents dans la couche 1
!
       AcoefT(j) = AcoefT_x(j) + ywake_s(j)*dd_AT(j)*(1.+sigx(j)*LambdaTs)
       AcoefQ(j) = AcoefQ_x(j) + ywake_s(j)*dd_AQ(j)*(1.+sigx(j)*LambdaQs)
       AcoefU(j) = AcoefU_x(j) + ywake_s(j)*dd_AU(j)*(1.+sigx(j)*LambdaUs)
       AcoefV(j) = AcoefV_x(j) + ywake_s(j)*dd_AV(j)*(1.+sigx(j)*LambdaVs)
!                                           
       BcoefT(j) = BcoefT_x(j) + ywake_s(j)*BcoefT_x(j)*sigx(j)*LambdaT*LambdaTs &
                               + ywake_s(j)*dd_BT(j)*(1.+sigx(j)*LambdaT)*(1.+sigx(j)*LambdaTs)
                                           
       BcoefQ(j) = BcoefQ_x(j) + ywake_s(j)*BcoefQ_x(j)*sigx(j)*LambdaQ*LambdaQs &
                               + ywake_s(j)*dd_BQ(j)*(1.+sigx(j)*LambdaQ)*(1.+sigx(j)*LambdaQs)
                                           
       BcoefU(j) = BcoefU_x(j) + ywake_s(j)*BcoefU_x(j)*sigx(j)*LambdaU*LambdaUs &
                               + ywake_s(j)*dd_BU(j)*(1.+sigx(j)*LambdaU)*(1.+sigx(j)*LambdaUs)
                                           
       BcoefV(j) = BcoefV_x(j) + ywake_s(j)*BcoefV_x(j)*sigx(j)*LambdaV*LambdaVs &
                               + ywake_s(j)*dd_BV(j)*(1.+sigx(j)*LambdaV)*(1.+sigx(j)*LambdaVs)

!
! Calcul des cdrag \'equivalents dans la couche 
!
       ycdragm(j) = ycdragm_x(j) + ywake_s(j)*dd_Cdragm
       ycdragh(j) = ycdragh_x(j) + ywake_s(j)*dd_Cdragh
!
! Calcul de T, q, u et v \'equivalents dans la couche 1
!!       yt1(j) = yt_x(j,1) + ywake_s(j)*dd_t*(1.+sigx(j)*dd_Kh/KCT)
!!       yq1(j) = yq_x(j,1) + ywake_s(j)*dd_q*(1.+sigx(j)*dd_Kh/KCQ)
!!       yu1(j) = yu_x(j,1) + ywake_s(j)*dd_u*(1.+sigx(j)*dd_Km/KCU)
!!       yv1(j) = yv_x(j,1) + ywake_s(j)*dd_v*(1.+sigx(j)*dd_Km/KCV)
       yt1(j) = yt_x(j,1) + ywake_s(j)*dd_t
       yq1(j) = yq_x(j,1) + ywake_s(j)*dd_q
       yu1(j) = yu_x(j,1) + ywake_s(j)*dd_u
       yv1(j) = yv_x(j,1) + ywake_s(j)*dd_v


        ENDDO

        RETURN

END SUBROUTINE wx_pbl0_fuse

SUBROUTINE wx_pbl0_split(knon, dtime, ywake_s, &
                       y_flux_t1, y_flux_q1, y_flux_u1, y_flux_v1, &
                       y_flux_t1_x, y_flux_t1_w, &
                       y_flux_q1_x, y_flux_q1_w, &
                       y_flux_u1_x, y_flux_u1_w, &
                       y_flux_v1_x, y_flux_v1_w, &
                       yfluxlat_x, yfluxlat_w, &
                       y_delta_tsurf &
                       )
!
    USE print_control_mod, ONLY: prt_level,lunout
!
    INCLUDE "YOMCST.h"
!
    INTEGER,                      INTENT(IN)        :: knon    ! number of grid cells
    REAL,                         INTENT(IN)        :: dtime   ! time step size (s)
    REAL, DIMENSION(knon),        INTENT(IN)        :: ywake_s ! cold pools fractional area
    REAL, DIMENSION(knon),        INTENT(IN)        :: y_flux_t1, y_flux_q1, y_flux_u1, y_flux_v1
!
    REAL, DIMENSION(knon),        INTENT(OUT)       :: y_flux_t1_x, y_flux_t1_w
    REAL, DIMENSION(knon),        INTENT(OUT)       :: y_flux_q1_x, y_flux_q1_w
    REAL, DIMENSION(knon),        INTENT(OUT)       :: y_flux_u1_x, y_flux_u1_w
    REAL, DIMENSION(knon),        INTENT(OUT)       :: y_flux_v1_x, y_flux_v1_w
    REAL, DIMENSION(knon),        INTENT(OUT)       :: yfluxlat_x, yfluxlat_w
    REAL, DIMENSION(knon),        INTENT(OUT)       :: y_delta_tsurf
!
!! Local variables
    INTEGER                    :: j
    REAL, DIMENSION(knon)      :: y_delta_flux_t1, y_delta_flux_q1, y_delta_flux_u1, y_delta_flux_v1
!
    REAL                       :: DDT, DDQ, DDU, DDV
    REAL                       :: LambdaTs, LambdaQs, LambdaUs, LambdaVs 
!
    REAL, DIMENSION(knon)      :: sigx       ! fractional area of (x) region
!!
        sigx(:) = 1.-ywake_s(:)

        DO j=1,knon
!
       DDT = Kech_Tp(j)
       DDQ = Kech_Qp(j)
       DDU = Kech_Up(j)
       DDV = Kech_Vp(j)
!
       LambdaTs =  dd_KTp(j)/DDT
       LambdaQs =  dd_KQp(j)/DDQ
       LambdaUs =  dd_KUp(j)/DDU
       LambdaVs =  dd_KVp(j)/DDV
!
         y_delta_flux_t1(j) = y_flux_t1(j)*LambdaTs + dd_AT(j)*KxKwTp(j)/DDT
         y_delta_flux_q1(j) = y_flux_q1(j)*LambdaQs + dd_AQ(j)*KxKwQp(j)/DDQ
         y_delta_flux_u1(j) = y_flux_u1(j)*LambdaUs + dd_AU(j)*KxKwUp(j)/DDU
         y_delta_flux_v1(j) = y_flux_v1(j)*LambdaVs + dd_AV(j)*KxKwVp(j)/DDV
!
         y_flux_t1_x(j)=y_flux_t1(j) - ywake_s(j)*y_delta_flux_t1(j)
         y_flux_t1_w(j)=y_flux_t1(j) + (1.-ywake_s(j))*y_delta_flux_t1(j)
         y_flux_q1_x(j)=y_flux_q1(j) - ywake_s(j)*y_delta_flux_q1(j)
         y_flux_q1_w(j)=y_flux_q1(j) + (1.-ywake_s(j))*y_delta_flux_q1(j)
         y_flux_u1_x(j)=y_flux_u1(j) - ywake_s(j)*y_delta_flux_u1(j)
         y_flux_u1_w(j)=y_flux_u1(j) + (1.-ywake_s(j))*y_delta_flux_u1(j)
         y_flux_v1_x(j)=y_flux_v1(j) - ywake_s(j)*y_delta_flux_v1(j)
         y_flux_v1_w(j)=y_flux_v1(j) + (1.-ywake_s(j))*y_delta_flux_v1(j)
!
         yfluxlat_x(j)=y_flux_q1_x(j)*RLVTT
         yfluxlat_w(j)=y_flux_q1_w(j)*RLVTT
!
!       Delta_tsurf computation
!!         y_delta_tsurf(j) = (1./RCPD)*(ah(j)*dd_AT(j) + &
!!                                       ah(j)*y_flux_t1(j)*dd_BT(j)*dtime + &
!!                                       y_delta_flux_t1(j)*(ah(j)*BBT+bh(j)) )
!
           y_delta_tsurf(j) = 0.
!
        ENDDO
!
        RETURN

END SUBROUTINE wx_pbl0_split

SUBROUTINE wx_pbl_final
!
!****************************************************************************************
! Deallocate module variables
!
!****************************************************************************************    
!
    IF (ALLOCATED(Kech_Tp))        DEALLOCATE(Kech_Tp)
    IF (ALLOCATED(Kech_T_xp))      DEALLOCATE(Kech_T_xp)
    IF (ALLOCATED(Kech_T_wp))      DEALLOCATE(Kech_T_wp)
    IF (ALLOCATED(dd_KTp))         DEALLOCATE(dd_KTp)
    IF (ALLOCATED(KxKwTp))         DEALLOCATE(KxKwTp)
    IF (ALLOCATED(dd_AT))          DEALLOCATE(dd_AT)
    IF (ALLOCATED(dd_BT))          DEALLOCATE(dd_BT)
    IF (ALLOCATED(Kech_Qp))        DEALLOCATE(Kech_Qp)
    IF (ALLOCATED(Kech_Q_xp))      DEALLOCATE(Kech_Q_xp)
    IF (ALLOCATED(Kech_Q_wp))      DEALLOCATE(Kech_Q_wp)
    IF (ALLOCATED(dd_KQp))         DEALLOCATE(dd_KQp)
    IF (ALLOCATED(KxKwQp))         DEALLOCATE(KxKwQp)
    IF (ALLOCATED(dd_AQ))          DEALLOCATE(dd_AQ)
    IF (ALLOCATED(dd_BQ))          DEALLOCATE(dd_BQ)
    IF (ALLOCATED(Kech_Up))        DEALLOCATE(Kech_Up)
    IF (ALLOCATED(Kech_U_xp))      DEALLOCATE(Kech_U_xp)
    IF (ALLOCATED(Kech_U_wp))      DEALLOCATE(Kech_U_wp)
    IF (ALLOCATED(dd_KUp))         DEALLOCATE(dd_KUp)
    IF (ALLOCATED(KxKwUp))         DEALLOCATE(KxKwUp)
    IF (ALLOCATED(dd_AU))          DEALLOCATE(dd_AU)
    IF (ALLOCATED(dd_BU))          DEALLOCATE(dd_BU)
    IF (ALLOCATED(Kech_Vp))        DEALLOCATE(Kech_Vp)
    IF (ALLOCATED(Kech_V_xp))      DEALLOCATE(Kech_V_xp)
    IF (ALLOCATED(Kech_V_wp))      DEALLOCATE(Kech_V_wp)
    IF (ALLOCATED(KxKwVp))         DEALLOCATE(KxKwVp)
    IF (ALLOCATED(dd_KVp))         DEALLOCATE(dd_KVp)
    IF (ALLOCATED(dd_AV))          DEALLOCATE(dd_AV)
    IF (ALLOCATED(dd_BV))          DEALLOCATE(dd_BV)

END SUBROUTINE wx_pbl_final

END MODULE wx_pbl_mod

