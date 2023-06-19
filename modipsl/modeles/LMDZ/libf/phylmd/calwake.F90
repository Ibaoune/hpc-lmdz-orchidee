
! $Id: calwake.F90 3208 2018-02-16 11:42:18Z jyg $

SUBROUTINE calwake(iflag_wake_tend, paprs, pplay, dtime, &
    t, q, omgb, &
    dt_dwn, dq_dwn, m_dwn, m_up, dt_a, dq_a, wgen, &
    sigd, Cin, &
    wake_deltat, wake_deltaq, wake_s, awake_dens, wake_dens, &
    wake_dth, wake_h, &
    wake_pe, wake_fip, wake_gfl, &
    dt_wake, dq_wake, wake_k, t_x, q_x, wake_omgbdth, &
    wake_dp_omgb, &
    wake_dtke, wake_dqke, &
    wake_omg, wake_dp_deltomg, &
    wake_spread, wake_cstar, wake_d_deltat_gw, &
    wake_ddeltat, wake_ddeltaq, wake_ds, awake_ddens, wake_ddens)
  ! **************************************************************
  ! *
  ! CALWAKE                                                     *
  ! interface avec le schema de calcul de la poche    *
  ! froide                                            *
  ! *
  ! written by   : CHERUY Frederique, 13/03/2000, 10.31.05      *
  ! modified by :  ROEHRIG Romain,    01/30/2007                *
  ! **************************************************************

  USE dimphy
  USE phys_state_var_mod, ONLY: pctsrf
  USE indice_sol_mod, ONLY: is_oce
  USE print_control_mod, ONLY: mydebug=>debug , lunout, prt_level
  IMPLICIT NONE
  ! ======================================================================
  include "YOMCST.h"

  ! Arguments
  ! ----------
  ! Input
  ! ----
  INTEGER,                       INTENT (IN)         :: iflag_wake_tend
  REAL,                          INTENT (IN)         :: dtime
  REAL, DIMENSION(klon, klev),   INTENT (IN)         :: pplay
  REAL, DIMENSION(klon, klev+1), INTENT (IN)         :: paprs
  REAL, DIMENSION(klon, klev),   INTENT (IN)         :: t, q, omgb
  REAL, DIMENSION(klon, klev),   INTENT (IN)         :: dt_dwn, dq_dwn
  REAL, DIMENSION(klon, klev),   INTENT (IN)         :: m_up, m_dwn
  REAL, DIMENSION(klon, klev),   INTENT (IN)         :: dt_a, dq_a
  REAL, DIMENSION(klon),         INTENT (IN)         :: wgen
  REAL, DIMENSION(klon),         INTENT (IN)         :: sigd
  REAL, DIMENSION(klon),         INTENT (IN)         :: Cin
  ! Input/Output
  ! ------------
  REAL, DIMENSION(klon, klev),   INTENT (INOUT)      :: wake_deltat, wake_deltaq
  REAL, DIMENSION(klon),         INTENT (INOUT)      :: wake_s
  REAL, DIMENSION(klon),         INTENT (INOUT)      :: awake_dens, wake_dens
  ! Output
  ! ------
  REAL, DIMENSION(klon, klev),   INTENT (OUT)        :: dt_wake, dq_wake
!!jyg  REAL, DIMENSION(klon),         INTENT (OUT)        :: wake_k
  INTEGER, DIMENSION(klon),      INTENT (OUT)        :: wake_k
  REAL, DIMENSION(klon, klev),   INTENT (OUT)        :: wake_d_deltat_gw
  REAL, DIMENSION(klon),         INTENT (OUT)        :: wake_h
  REAL, DIMENSION(klon, klev),   INTENT (OUT)        :: wake_dth
  REAL, DIMENSION(klon),         INTENT (OUT)        :: wake_pe, wake_fip, wake_gfl
  REAL, DIMENSION(klon, klev),   INTENT (OUT)        :: t_x, q_x
  REAL, DIMENSION(klon, klev),   INTENT (OUT)        :: wake_omgbdth, wake_dp_omgb
  REAL, DIMENSION(klon, klev),   INTENT (OUT)        :: wake_dtke, wake_dqke
  REAL, DIMENSION(klon, klev),   INTENT (OUT)        :: wake_omg, wake_dp_deltomg
  REAL, DIMENSION(klon, klev),   INTENT (OUT)        :: wake_spread
  REAL, DIMENSION(klon),         INTENT (OUT)        :: wake_cstar
  REAL, DIMENSION(klon, klev),   INTENT (OUT)        :: wake_ddeltat, wake_ddeltaq
  REAL, DIMENSION(klon),         INTENT (OUT)        :: wake_ds, awake_ddens, wake_ddens


  ! Variable internes
  ! -----------------
  LOGICAL, SAVE                                      :: first = .TRUE.
  !$OMP THREADPRIVATE(first)
  INTEGER                                            :: i, l
  INTEGER, DIMENSION(klon)                           :: znatsurf    ! 0 if pctsrf(is_oce)>0.1; 1 else.
  REAL                                               :: aire
  REAL, DIMENSION(klon, klev)                        :: p,  pi
  REAL, DIMENSION(klon, klev+1)                      ::  ph
  REAL, DIMENSION(klon, klev)                        ::  omgbe
  REAL, DIMENSION(klon, klev)                        :: te, qe
  REAL, DIMENSION(klon, klev)                        :: dtdwn, dqdwn
  REAL, DIMENSION(klon, klev)                        :: dta, dqa
  REAL, DIMENSION(klon, klev)                        :: amdwn, amup
  REAL, DIMENSION(klon, klev)                        :: dtw, dqw, dth
  REAL, DIMENSION(klon, klev)                        :: dtls, dqls
  REAL, DIMENSION(klon, klev)                        :: tx, qx
  REAL, DIMENSION(klon)                              :: hw, wape, fip, gfl
  REAL, DIMENSION(klon)                              :: sigmaw, awdens, wdens
  REAL, DIMENSION(klon, klev)                        :: omgbdth
  REAL, DIMENSION(klon, klev)                        :: dp_omgb
  REAL, DIMENSION(klon, klev)                        :: dtke, dqke
  REAL, DIMENSION(klon, klev)                        :: omg
  REAL, DIMENSION(klon, klev)                        :: dp_deltomg, spread
  REAL, DIMENSION(klon)                              :: cstar
  REAL, DIMENSION(klon)                              :: sigd0
  INTEGER, DIMENSION(klon)                           :: ktopw
  REAL, DIMENSION(klon, klev)                        :: d_deltat_gw
  REAL, DIMENSION(klon, klev)                        :: d_deltatw, d_deltaqw
  REAL, DIMENSION(klon)                              :: d_sigmaw, d_awdens, d_wdens

  REAL                                               :: rdcp


  IF (prt_level >= 10) THEN
    print *, '-> calwake, wake_s, wgen input ', wake_s(1), wgen(1)
  ENDIF

  rdcp = 1./3.5

  znatsurf(:) = 0
  DO i = 1,klon
    IF (pctsrf(i,is_oce) < 0.1) znatsurf(i) = 1
  ENDDO


  ! -----------------------------------------------------------
  ! IM 290108     DO 999 i=1,klon   ! a vectoriser
  ! ----------------------------------------------------------


  DO l = 1, klev
    DO i = 1, klon
      p(i, l) = pplay(i, l)
      ph(i, l) = paprs(i, l)
      pi(i, l) = (pplay(i,l)/100000.)**rdcp

      te(i, l) = t(i, l)
      qe(i, l) = q(i, l)
      omgbe(i, l) = omgb(i, l)

      dtdwn(i, l) = dt_dwn(i, l)
      dqdwn(i, l) = dq_dwn(i, l)
      dta(i, l) = dt_a(i, l)
      dqa(i, l) = dq_a(i, l)
    END DO
  END DO

!----------------------------------------------------------------
!         Initialize tendencies to zero
!----------------------------------------------------------------
dtls(:,:) = 0.
dqls(:,:) = 0.
d_deltat_gw(:,:) = 0.
d_deltatw(:,:) = 0.
d_deltaqw(:,:) = 0.
d_sigmaw(:) = 0.
d_awdens(:) = 0.
d_wdens(:) = 0.
!

  DO i = 1, klon
    sigd0(i) = sigd(i)
  END DO
  ! print*, 'sigd0,sigd', sigd0, sigd(i)
  DO i = 1, klon
    ph(i, klev+1) = 0.
  END DO

!!jyg!  DO i = 1, klon                  
!!jyg!    ktopw(i) = NINT(wake_k(i))    
!!jyg!  END DO                          

  DO i = 1, klon
    hw(i) = wake_h(i)
  END DO
!
!    Make a copy of state variables
  DO l = 1, klev
    DO i = 1, klon
      dtw(i, l) = wake_deltat(i, l)
      dqw(i, l) = wake_deltaq(i, l)
    END DO
  END DO

  DO i = 1, klon
    sigmaw(i) = wake_s(i)
  END DO

  DO i = 1, klon
    awdens(i) = max(0., awake_dens(i))
    wdens(i) = max(0., wake_dens(i))
  END DO

  ! fkc les flux de masses sont evalues aux niveaux et valent 0 a la surface
  ! fkc  on veut le flux de masse au milieu des couches

  DO l = 1, klev - 1
    DO i = 1, klon
      amdwn(i, l) = 0.5*(m_dwn(i,l)+m_dwn(i,l+1))
      amdwn(i, l) = (m_dwn(i,l+1))
    END DO
  END DO

  ! au sommet le flux de masse est nul

  DO i = 1, klon
    amdwn(i, klev) = 0.5*m_dwn(i, klev)
  END DO

  DO l = 1, klev
    DO i = 1, klon
      amup(i, l) = m_up(i, l)
    END DO
  END DO

  CALL wake(znatsurf, p, ph, pi, dtime, &
    te, qe, omgbe, &
    dtdwn, dqdwn, amdwn, amup, dta, dqa, wgen, &
    sigd0, Cin, &
    dtw, dqw, sigmaw, awdens, wdens, &                                   ! state variables
    dth, hw, wape, fip, gfl, &
    dtls, dqls, ktopw, omgbdth, dp_omgb, tx, qx, &
    dtke, dqke, omg, dp_deltomg, spread, cstar, &
    d_deltat_gw, &
    d_deltatw, d_deltaqw, d_sigmaw, d_awdens, d_wdens)                     ! tendencies

!
  DO l = 1, klev
    DO i = 1, klon
      IF (ktopw(i)>0) THEN
        wake_d_deltat_gw(i, l) = d_deltat_gw(i, l)
        wake_omgbdth(i, l) = omgbdth(i, l)
        wake_dp_omgb(i, l) = dp_omgb(i, l)
        wake_dtke(i, l) = dtke(i, l)
        wake_dqke(i, l) = dqke(i, l)
        wake_omg(i, l) = omg(i, l)
        wake_dp_deltomg(i, l) = dp_deltomg(i, l)
        wake_spread(i, l) = spread(i, l)
        wake_dth(i, l) = dth(i, l)
        dt_wake(i, l) = dtls(i, l)*dtime         ! derivative -> tendency
        dq_wake(i, l) = dqls(i, l)*dtime         ! derivative -> tendency
        t_x(i, l) = tx(i, l)
        q_x(i, l) = qx(i, l)
      ELSE
        wake_d_deltat_gw(i, l) = 0.
        wake_omgbdth(i, l) = 0.
        wake_dp_omgb(i, l) = 0.
        wake_dtke(i, l) = 0.
        wake_dqke(i, l) = 0.
        wake_omg(i, l) = 0.
        wake_dp_deltomg(i, l) = 0.
        wake_spread(i, l) = 0.
        wake_dth(i, l) = 0.
        dt_wake(i, l) = 0.
        dq_wake(i, l) = 0.
        t_x(i, l) = te(i, l)
        q_x(i, l) = qe(i, l)
      END IF
    END DO
  END DO

  DO i = 1, klon
    wake_h(i) = hw(i)
    wake_pe(i) = wape(i)
    wake_fip(i) = fip(i)
    wake_gfl(i) = gfl(i)
    wake_k(i) = ktopw(i)
    wake_cstar(i) = cstar(i)
  END DO

!  Tendencies of state variables
  DO l = 1, klev
    DO i = 1, klon
      IF (ktopw(i)>0) THEN
        wake_ddeltat(i, l) = d_deltatw(i, l)*dtime
        wake_ddeltaq(i, l) = d_deltaqw(i, l)*dtime
      ELSE
        wake_ddeltat(i, l) = -wake_deltat(i, l)
        wake_ddeltaq(i, l) = -wake_deltaq(i, l)
      END IF
    END DO
  END DO
  DO i = 1, klon
    IF (ktopw(i)>0) THEN
      wake_ds(i) = d_sigmaw(i)*dtime
      awake_ddens(i) = d_awdens(i)*dtime
      wake_ddens(i) = d_wdens(i)*dtime
    ELSE
      wake_ds(i)   = -wake_s(i)
      wake_ddens(i)= -wake_dens(i)
    END IF
  END DO
!

!jyg<  
  IF (iflag_wake_tend .EQ. 0) THEN
!  Update State variables
    DO l = 1, klev
      DO i = 1, klon
        IF (ktopw(i)>0) THEN
          wake_deltat(i, l) = dtw(i, l)
          wake_deltaq(i, l) = dqw(i, l)
        ELSE
          wake_deltat(i, l) = 0.
          wake_deltaq(i, l) = 0.
        END IF
      END DO
    END DO
    DO i = 1, klon
      wake_s(i) = sigmaw(i)
      awake_dens(i) = awdens(i)
      wake_dens(i) = wdens(i)
    END DO
  ENDIF  ! (iflag_wake_tend .EQ. 0)
!
  IF (first) THEN
    DO i = 1,klon
      IF (wake_dens(i) < -1.) THEN
        wake_dens(i) = wdens(i)
      ENDIF
    ENDDO
    first=.false.
  ENDIF  ! (first)
!>jyg

  RETURN
END SUBROUTINE calwake

