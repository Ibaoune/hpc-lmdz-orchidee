
! $Id: cv_routines.F90 2311 2015-06-25 07:45:24Z emillour $

SUBROUTINE cv_param(nd)
  IMPLICIT NONE

  ! ------------------------------------------------------------
  ! Set parameters for convectL
  ! (includes microphysical parameters and parameters that
  ! control the rate of approach to quasi-equilibrium)
  ! ------------------------------------------------------------

  ! *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
  ! ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
  ! ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
  ! ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
  ! ***               BETWEEN 0 C AND TLCRIT)                        ***
  ! ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
  ! ***                       FORMULATION                            ***
  ! ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
  ! ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
  ! ***                        OF CLOUD                              ***
  ! ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
  ! ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
  ! ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
  ! ***                          OF RAIN                             ***
  ! ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
  ! ***                          OF SNOW                             ***
  ! ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
  ! ***                         TRANSPORT                            ***
  ! ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
  ! ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
  ! ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
  ! ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
  ! ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
  ! ***                   (DAMP MUST BE LESS THAN 1)                 ***

  include "cvparam.h"
  INTEGER nd
  CHARACTER (LEN=20) :: modname = 'cv_routines'
  CHARACTER (LEN=80) :: abort_message

  ! noff: integer limit for convection (nd-noff)
  ! minorig: First level of convection

  noff = 2
  minorig = 2

  nl = nd - noff
  nlp = nl + 1
  nlm = nl - 1

  elcrit = 0.0011
  tlcrit = -55.0
  entp = 1.5
  sigs = 0.12
  sigd = 0.05
  omtrain = 50.0
  omtsnow = 5.5
  coeffr = 1.0
  coeffs = 0.8
  dtmax = 0.9

  cu = 0.70

  betad = 10.0

  damp = 0.1
  alpha = 0.2

  delta = 0.01 ! cld

  RETURN
END SUBROUTINE cv_param

SUBROUTINE cv_prelim(len, nd, ndp1, t, q, p, ph, lv, cpn, tv, gz, h, hm)
  IMPLICIT NONE

  ! =====================================================================
  ! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
  ! =====================================================================

  ! inputs:
  INTEGER len, nd, ndp1
  REAL t(len, nd), q(len, nd), p(len, nd), ph(len, ndp1)

  ! outputs:
  REAL lv(len, nd), cpn(len, nd), tv(len, nd)
  REAL gz(len, nd), h(len, nd), hm(len, nd)

  ! local variables:
  INTEGER k, i
  REAL cpx(len, nd)

  include "cvthermo.h"
  include "cvparam.h"


  DO k = 1, nlp
    DO i = 1, len
      lv(i, k) = lv0 - clmcpv*(t(i,k)-t0)
      cpn(i, k) = cpd*(1.0-q(i,k)) + cpv*q(i, k)
      cpx(i, k) = cpd*(1.0-q(i,k)) + cl*q(i, k)
      tv(i, k) = t(i, k)*(1.0+q(i,k)*epsim1)
    END DO
  END DO

  ! gz = phi at the full levels (same as p).

  DO i = 1, len
    gz(i, 1) = 0.0
  END DO
  DO k = 2, nlp
    DO i = 1, len
      gz(i, k) = gz(i, k-1) + hrd*(tv(i,k-1)+tv(i,k))*(p(i,k-1)-p(i,k))/ph(i, &
        k)
    END DO
  END DO

  ! h  = phi + cpT (dry static energy).
  ! hm = phi + cp(T-Tbase)+Lq

  DO k = 1, nlp
    DO i = 1, len
      h(i, k) = gz(i, k) + cpn(i, k)*t(i, k)
      hm(i, k) = gz(i, k) + cpx(i, k)*(t(i,k)-t(i,1)) + lv(i, k)*q(i, k)
    END DO
  END DO

  RETURN
END SUBROUTINE cv_prelim

SUBROUTINE cv_feed(len, nd, t, q, qs, p, hm, gz, nk, icb, icbmax, iflag, tnk, &
    qnk, gznk, plcl)
  IMPLICIT NONE

  ! ================================================================
  ! Purpose: CONVECTIVE FEED
  ! ================================================================

  include "cvparam.h"

  ! inputs:
  INTEGER len, nd
  REAL t(len, nd), q(len, nd), qs(len, nd), p(len, nd)
  REAL hm(len, nd), gz(len, nd)

  ! outputs:
  INTEGER iflag(len), nk(len), icb(len), icbmax
  REAL tnk(len), qnk(len), gznk(len), plcl(len)

  ! local variables:
  INTEGER i, k
  INTEGER ihmin(len)
  REAL work(len)
  REAL pnk(len), qsnk(len), rh(len), chi(len)

  ! -------------------------------------------------------------------
  ! --- Find level of minimum moist static energy
  ! --- If level of minimum moist static energy coincides with
  ! --- or is lower than minimum allowable parcel origin level,
  ! --- set iflag to 6.
  ! -------------------------------------------------------------------

  DO i = 1, len
    work(i) = 1.0E12
    ihmin(i) = nl
  END DO
  DO k = 2, nlp
    DO i = 1, len
      IF ((hm(i,k)<work(i)) .AND. (hm(i,k)<hm(i,k-1))) THEN
        work(i) = hm(i, k)
        ihmin(i) = k
      END IF
    END DO
  END DO
  DO i = 1, len
    ihmin(i) = min(ihmin(i), nlm)
    IF (ihmin(i)<=minorig) THEN
      iflag(i) = 6
    END IF
  END DO

  ! -------------------------------------------------------------------
  ! --- Find that model level below the level of minimum moist static
  ! --- energy that has the maximum value of moist static energy
  ! -------------------------------------------------------------------

  DO i = 1, len
    work(i) = hm(i, minorig)
    nk(i) = minorig
  END DO
  DO k = minorig + 1, nl
    DO i = 1, len
      IF ((hm(i,k)>work(i)) .AND. (k<=ihmin(i))) THEN
        work(i) = hm(i, k)
        nk(i) = k
      END IF
    END DO
  END DO
  ! -------------------------------------------------------------------
  ! --- Check whether parcel level temperature and specific humidity
  ! --- are reasonable
  ! -------------------------------------------------------------------
  DO i = 1, len
    IF (((t(i,nk(i))<250.0) .OR. (q(i,nk(i))<=0.0) .OR. (p(i,ihmin(i))< &
      400.0)) .AND. (iflag(i)==0)) iflag(i) = 7
  END DO
  ! -------------------------------------------------------------------
  ! --- Calculate lifted condensation level of air at parcel origin level
  ! --- (Within 0.2% of formula of Bolton, MON. WEA. REV.,1980)
  ! -------------------------------------------------------------------
  DO i = 1, len
    tnk(i) = t(i, nk(i))
    qnk(i) = q(i, nk(i))
    gznk(i) = gz(i, nk(i))
    pnk(i) = p(i, nk(i))
    qsnk(i) = qs(i, nk(i))

    rh(i) = qnk(i)/qsnk(i)
    rh(i) = min(1.0, rh(i))
    chi(i) = tnk(i)/(1669.0-122.0*rh(i)-tnk(i))
    plcl(i) = pnk(i)*(rh(i)**chi(i))
    IF (((plcl(i)<200.0) .OR. (plcl(i)>=2000.0)) .AND. (iflag(i)==0)) iflag(i &
      ) = 8
  END DO
  ! -------------------------------------------------------------------
  ! --- Calculate first level above lcl (=icb)
  ! -------------------------------------------------------------------
  DO i = 1, len
    icb(i) = nlm
  END DO

  DO k = minorig, nl
    DO i = 1, len
      IF ((k>=(nk(i)+1)) .AND. (p(i,k)<plcl(i))) icb(i) = min(icb(i), k)
    END DO
  END DO

  DO i = 1, len
    IF ((icb(i)>=nlm) .AND. (iflag(i)==0)) iflag(i) = 9
  END DO

  ! Compute icbmax.

  icbmax = 2
  DO i = 1, len
    icbmax = max(icbmax, icb(i))
  END DO

  RETURN
END SUBROUTINE cv_feed

SUBROUTINE cv_undilute1(len, nd, t, q, qs, gz, p, nk, icb, icbmax, tp, tvp, &
    clw)
  IMPLICIT NONE

  include "cvthermo.h"
  include "cvparam.h"

  ! inputs:
  INTEGER len, nd
  INTEGER nk(len), icb(len), icbmax
  REAL t(len, nd), q(len, nd), qs(len, nd), gz(len, nd)
  REAL p(len, nd)

  ! outputs:
  REAL tp(len, nd), tvp(len, nd), clw(len, nd)

  ! local variables:
  INTEGER i, k
  REAL tg, qg, alv, s, ahg, tc, denom, es, rg
  REAL ah0(len), cpp(len)
  REAL tnk(len), qnk(len), gznk(len), ticb(len), gzicb(len)

  ! -------------------------------------------------------------------
  ! --- Calculates the lifted parcel virtual temperature at nk,
  ! --- the actual temperature, and the adiabatic
  ! --- liquid water content. The procedure is to solve the equation.
  ! cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.
  ! -------------------------------------------------------------------

  DO i = 1, len
    tnk(i) = t(i, nk(i))
    qnk(i) = q(i, nk(i))
    gznk(i) = gz(i, nk(i))
    ticb(i) = t(i, icb(i))
    gzicb(i) = gz(i, icb(i))
  END DO

  ! ***  Calculate certain parcel quantities, including static energy   ***

  DO i = 1, len
    ah0(i) = (cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i) + qnk(i)*(lv0-clmcpv*(tnk(i)- &
      273.15)) + gznk(i)
    cpp(i) = cpd*(1.-qnk(i)) + qnk(i)*cpv
  END DO

  ! ***   Calculate lifted parcel quantities below cloud base   ***

  DO k = minorig, icbmax - 1
    DO i = 1, len
      tp(i, k) = tnk(i) - (gz(i,k)-gznk(i))/cpp(i)
      tvp(i, k) = tp(i, k)*(1.+qnk(i)*epsi)
    END DO
  END DO

  ! ***  Find lifted parcel quantities above cloud base    ***

  DO i = 1, len
    tg = ticb(i)
    qg = qs(i, icb(i))
    alv = lv0 - clmcpv*(ticb(i)-t0)

    ! First iteration.

    s = cpd + alv*alv*qg/(rrv*ticb(i)*ticb(i))
    s = 1./s
    ahg = cpd*tg + (cl-cpd)*qnk(i)*ticb(i) + alv*qg + gzicb(i)
    tg = tg + s*(ah0(i)-ahg)
    tg = max(tg, 35.0)
    tc = tg - t0
    denom = 243.5 + tc
    IF (tc>=0.0) THEN
      es = 6.112*exp(17.67*tc/denom)
    ELSE
      es = exp(23.33086-6111.72784/tg+0.15215*log(tg))
    END IF
    qg = eps*es/(p(i,icb(i))-es*(1.-eps))

    ! Second iteration.

    s = cpd + alv*alv*qg/(rrv*ticb(i)*ticb(i))
    s = 1./s
    ahg = cpd*tg + (cl-cpd)*qnk(i)*ticb(i) + alv*qg + gzicb(i)
    tg = tg + s*(ah0(i)-ahg)
    tg = max(tg, 35.0)
    tc = tg - t0
    denom = 243.5 + tc
    IF (tc>=0.0) THEN
      es = 6.112*exp(17.67*tc/denom)
    ELSE
      es = exp(23.33086-6111.72784/tg+0.15215*log(tg))
    END IF
    qg = eps*es/(p(i,icb(i))-es*(1.-eps))

    alv = lv0 - clmcpv*(ticb(i)-273.15)
    tp(i, icb(i)) = (ah0(i)-(cl-cpd)*qnk(i)*ticb(i)-gz(i,icb(i))-alv*qg)/cpd
    clw(i, icb(i)) = qnk(i) - qg
    clw(i, icb(i)) = max(0.0, clw(i,icb(i)))
    rg = qg/(1.-qnk(i))
    tvp(i, icb(i)) = tp(i, icb(i))*(1.+rg*epsi)
  END DO

  DO k = minorig, icbmax
    DO i = 1, len
      tvp(i, k) = tvp(i, k) - tp(i, k)*qnk(i)
    END DO
  END DO

  RETURN
END SUBROUTINE cv_undilute1

SUBROUTINE cv_trigger(len, nd, icb, cbmf, tv, tvp, iflag)
  IMPLICIT NONE

  ! -------------------------------------------------------------------
  ! --- Test for instability.
  ! --- If there was no convection at last time step and parcel
  ! --- is stable at icb, then set iflag to 4.
  ! -------------------------------------------------------------------

  include "cvparam.h"

  ! inputs:
  INTEGER len, nd, icb(len)
  REAL cbmf(len), tv(len, nd), tvp(len, nd)

  ! outputs:
  INTEGER iflag(len) ! also an input

  ! local variables:
  INTEGER i


  DO i = 1, len
    IF ((cbmf(i)==0.0) .AND. (iflag(i)==0) .AND. (tvp(i, &
      icb(i))<=(tv(i,icb(i))-dtmax))) iflag(i) = 4
  END DO

  RETURN
END SUBROUTINE cv_trigger

SUBROUTINE cv_compress(len, nloc, ncum, nd, iflag1, nk1, icb1, cbmf1, plcl1, &
    tnk1, qnk1, gznk1, t1, q1, qs1, u1, v1, gz1, h1, lv1, cpn1, p1, ph1, tv1, &
    tp1, tvp1, clw1, iflag, nk, icb, cbmf, plcl, tnk, qnk, gznk, t, q, qs, u, &
    v, gz, h, lv, cpn, p, ph, tv, tp, tvp, clw, dph)
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE

  include "cvparam.h"

  ! inputs:
  INTEGER len, ncum, nd, nloc
  INTEGER iflag1(len), nk1(len), icb1(len)
  REAL cbmf1(len), plcl1(len), tnk1(len), qnk1(len), gznk1(len)
  REAL t1(len, nd), q1(len, nd), qs1(len, nd), u1(len, nd), v1(len, nd)
  REAL gz1(len, nd), h1(len, nd), lv1(len, nd), cpn1(len, nd)
  REAL p1(len, nd), ph1(len, nd+1), tv1(len, nd), tp1(len, nd)
  REAL tvp1(len, nd), clw1(len, nd)

  ! outputs:
  INTEGER iflag(nloc), nk(nloc), icb(nloc)
  REAL cbmf(nloc), plcl(nloc), tnk(nloc), qnk(nloc), gznk(nloc)
  REAL t(nloc, nd), q(nloc, nd), qs(nloc, nd), u(nloc, nd), v(nloc, nd)
  REAL gz(nloc, nd), h(nloc, nd), lv(nloc, nd), cpn(nloc, nd)
  REAL p(nloc, nd), ph(nloc, nd+1), tv(nloc, nd), tp(nloc, nd)
  REAL tvp(nloc, nd), clw(nloc, nd)
  REAL dph(nloc, nd)

  ! local variables:
  INTEGER i, k, nn
  CHARACTER (LEN=20) :: modname = 'cv_compress'
  CHARACTER (LEN=80) :: abort_message


  DO k = 1, nl + 1
    nn = 0
    DO i = 1, len
      IF (iflag1(i)==0) THEN
        nn = nn + 1
        t(nn, k) = t1(i, k)
        q(nn, k) = q1(i, k)
        qs(nn, k) = qs1(i, k)
        u(nn, k) = u1(i, k)
        v(nn, k) = v1(i, k)
        gz(nn, k) = gz1(i, k)
        h(nn, k) = h1(i, k)
        lv(nn, k) = lv1(i, k)
        cpn(nn, k) = cpn1(i, k)
        p(nn, k) = p1(i, k)
        ph(nn, k) = ph1(i, k)
        tv(nn, k) = tv1(i, k)
        tp(nn, k) = tp1(i, k)
        tvp(nn, k) = tvp1(i, k)
        clw(nn, k) = clw1(i, k)
      END IF
    END DO
  END DO

  IF (nn/=ncum) THEN
    WRITE (lunout, *) 'strange! nn not equal to ncum: ', nn, ncum
    abort_message = ''
    CALL abort_physic(modname, abort_message, 1)
  END IF

  nn = 0
  DO i = 1, len
    IF (iflag1(i)==0) THEN
      nn = nn + 1
      cbmf(nn) = cbmf1(i)
      plcl(nn) = plcl1(i)
      tnk(nn) = tnk1(i)
      qnk(nn) = qnk1(i)
      gznk(nn) = gznk1(i)
      nk(nn) = nk1(i)
      icb(nn) = icb1(i)
      iflag(nn) = iflag1(i)
    END IF
  END DO

  DO k = 1, nl
    DO i = 1, ncum
      dph(i, k) = ph(i, k) - ph(i, k+1)
    END DO
  END DO

  RETURN
END SUBROUTINE cv_compress

SUBROUTINE cv_undilute2(nloc, ncum, nd, icb, nk, tnk, qnk, gznk, t, q, qs, &
    gz, p, dph, h, tv, lv, inb, inb1, tp, tvp, clw, hp, ep, sigp, frac)
  IMPLICIT NONE

  ! ---------------------------------------------------------------------
  ! Purpose:
  ! FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
  ! &
  ! COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
  ! FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
  ! &
  ! FIND THE LEVEL OF NEUTRAL BUOYANCY
  ! ---------------------------------------------------------------------

  include "cvthermo.h"
  include "cvparam.h"

  ! inputs:
  INTEGER ncum, nd, nloc
  INTEGER icb(nloc), nk(nloc)
  REAL t(nloc, nd), q(nloc, nd), qs(nloc, nd), gz(nloc, nd)
  REAL p(nloc, nd), dph(nloc, nd)
  REAL tnk(nloc), qnk(nloc), gznk(nloc)
  REAL lv(nloc, nd), tv(nloc, nd), h(nloc, nd)

  ! outputs:
  INTEGER inb(nloc), inb1(nloc)
  REAL tp(nloc, nd), tvp(nloc, nd), clw(nloc, nd)
  REAL ep(nloc, nd), sigp(nloc, nd), hp(nloc, nd)
  REAL frac(nloc)

  ! local variables:
  INTEGER i, k
  REAL tg, qg, ahg, alv, s, tc, es, denom, rg, tca, elacrit
  REAL by, defrac
  REAL ah0(nloc), cape(nloc), capem(nloc), byp(nloc)
  LOGICAL lcape(nloc)

  ! =====================================================================
  ! --- SOME INITIALIZATIONS
  ! =====================================================================

  DO k = 1, nl
    DO i = 1, ncum
      ep(i, k) = 0.0
      sigp(i, k) = sigs
    END DO
  END DO

  ! =====================================================================
  ! --- FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
  ! =====================================================================

  ! ---       The procedure is to solve the equation.
  ! cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.

  ! ***  Calculate certain parcel quantities, including static energy   ***


  DO i = 1, ncum
    ah0(i) = (cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i) + qnk(i)*(lv0-clmcpv*(tnk(i)- &
      t0)) + gznk(i)
  END DO


  ! ***  Find lifted parcel quantities above cloud base    ***


  DO k = minorig + 1, nl
    DO i = 1, ncum
      IF (k>=(icb(i)+1)) THEN
        tg = t(i, k)
        qg = qs(i, k)
        alv = lv0 - clmcpv*(t(i,k)-t0)

        ! First iteration.

        s = cpd + alv*alv*qg/(rrv*t(i,k)*t(i,k))
        s = 1./s
        ahg = cpd*tg + (cl-cpd)*qnk(i)*t(i, k) + alv*qg + gz(i, k)
        tg = tg + s*(ah0(i)-ahg)
        tg = max(tg, 35.0)
        tc = tg - t0
        denom = 243.5 + tc
        IF (tc>=0.0) THEN
          es = 6.112*exp(17.67*tc/denom)
        ELSE
          es = exp(23.33086-6111.72784/tg+0.15215*log(tg))
        END IF
        qg = eps*es/(p(i,k)-es*(1.-eps))

        ! Second iteration.

        s = cpd + alv*alv*qg/(rrv*t(i,k)*t(i,k))
        s = 1./s
        ahg = cpd*tg + (cl-cpd)*qnk(i)*t(i, k) + alv*qg + gz(i, k)
        tg = tg + s*(ah0(i)-ahg)
        tg = max(tg, 35.0)
        tc = tg - t0
        denom = 243.5 + tc
        IF (tc>=0.0) THEN
          es = 6.112*exp(17.67*tc/denom)
        ELSE
          es = exp(23.33086-6111.72784/tg+0.15215*log(tg))
        END IF
        qg = eps*es/(p(i,k)-es*(1.-eps))

        alv = lv0 - clmcpv*(t(i,k)-t0)
        ! print*,'cpd dans convect2 ',cpd
        ! print*,'tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd'
        ! print*,tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd
        tp(i, k) = (ah0(i)-(cl-cpd)*qnk(i)*t(i,k)-gz(i,k)-alv*qg)/cpd
        ! if (.not.cpd.gt.1000.) then
        ! print*,'CPD=',cpd
        ! stop
        ! endif
        clw(i, k) = qnk(i) - qg
        clw(i, k) = max(0.0, clw(i,k))
        rg = qg/(1.-qnk(i))
        tvp(i, k) = tp(i, k)*(1.+rg*epsi)
      END IF
    END DO
  END DO

  ! =====================================================================
  ! --- SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF
  ! --- PRECIPITATION FALLING OUTSIDE OF CLOUD
  ! --- THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)
  ! =====================================================================

  DO k = minorig + 1, nl
    DO i = 1, ncum
      IF (k>=(nk(i)+1)) THEN
        tca = tp(i, k) - t0
        IF (tca>=0.0) THEN
          elacrit = elcrit
        ELSE
          elacrit = elcrit*(1.0-tca/tlcrit)
        END IF
        elacrit = max(elacrit, 0.0)
        ep(i, k) = 1.0 - elacrit/max(clw(i,k), 1.0E-8)
        ep(i, k) = max(ep(i,k), 0.0)
        ep(i, k) = min(ep(i,k), 1.0)
        sigp(i, k) = sigs
      END IF
    END DO
  END DO

  ! =====================================================================
  ! --- CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL
  ! --- VIRTUAL TEMPERATURE
  ! =====================================================================

  DO k = minorig + 1, nl
    DO i = 1, ncum
      IF (k>=(icb(i)+1)) THEN
        tvp(i, k) = tvp(i, k)*(1.0-qnk(i)+ep(i,k)*clw(i,k))
        ! print*,'i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)'
        ! print*, i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)
      END IF
    END DO
  END DO
  DO i = 1, ncum
    tvp(i, nlp) = tvp(i, nl) - (gz(i,nlp)-gz(i,nl))/cpd
  END DO

  ! =====================================================================
  ! --- FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S
  ! --- HIGHEST LEVEL OF NEUTRAL BUOYANCY
  ! --- AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)
  ! =====================================================================

  DO i = 1, ncum
    cape(i) = 0.0
    capem(i) = 0.0
    inb(i) = icb(i) + 1
    inb1(i) = inb(i)
  END DO

  ! Originial Code

  ! do 530 k=minorig+1,nl-1
  ! do 520 i=1,ncum
  ! if(k.ge.(icb(i)+1))then
  ! by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
  ! byp=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
  ! cape(i)=cape(i)+by
  ! if(by.ge.0.0)inb1(i)=k+1
  ! if(cape(i).gt.0.0)then
  ! inb(i)=k+1
  ! capem(i)=cape(i)
  ! endif
  ! endif
  ! 520    continue
  ! 530  continue
  ! do 540 i=1,ncum
  ! byp=(tvp(i,nl)-tv(i,nl))*dph(i,nl)/p(i,nl)
  ! cape(i)=capem(i)+byp
  ! defrac=capem(i)-cape(i)
  ! defrac=max(defrac,0.001)
  ! frac(i)=-cape(i)/defrac
  ! frac(i)=min(frac(i),1.0)
  ! frac(i)=max(frac(i),0.0)
  ! 540   continue

  ! K Emanuel fix

  ! call zilch(byp,ncum)
  ! do 530 k=minorig+1,nl-1
  ! do 520 i=1,ncum
  ! if(k.ge.(icb(i)+1))then
  ! by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
  ! cape(i)=cape(i)+by
  ! if(by.ge.0.0)inb1(i)=k+1
  ! if(cape(i).gt.0.0)then
  ! inb(i)=k+1
  ! capem(i)=cape(i)
  ! byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
  ! endif
  ! endif
  ! 520    continue
  ! 530  continue
  ! do 540 i=1,ncum
  ! inb(i)=max(inb(i),inb1(i))
  ! cape(i)=capem(i)+byp(i)
  ! defrac=capem(i)-cape(i)
  ! defrac=max(defrac,0.001)
  ! frac(i)=-cape(i)/defrac
  ! frac(i)=min(frac(i),1.0)
  ! frac(i)=max(frac(i),0.0)
  ! 540   continue

  ! J Teixeira fix

  CALL zilch(byp, ncum)
  DO i = 1, ncum
    lcape(i) = .TRUE.
  END DO
  DO k = minorig + 1, nl - 1
    DO i = 1, ncum
      IF (cape(i)<0.0) lcape(i) = .FALSE.
      IF ((k>=(icb(i)+1)) .AND. lcape(i)) THEN
        by = (tvp(i,k)-tv(i,k))*dph(i, k)/p(i, k)
        byp(i) = (tvp(i,k+1)-tv(i,k+1))*dph(i, k+1)/p(i, k+1)
        cape(i) = cape(i) + by
        IF (by>=0.0) inb1(i) = k + 1
        IF (cape(i)>0.0) THEN
          inb(i) = k + 1
          capem(i) = cape(i)
        END IF
      END IF
    END DO
  END DO
  DO i = 1, ncum
    cape(i) = capem(i) + byp(i)
    defrac = capem(i) - cape(i)
    defrac = max(defrac, 0.001)
    frac(i) = -cape(i)/defrac
    frac(i) = min(frac(i), 1.0)
    frac(i) = max(frac(i), 0.0)
  END DO

  ! =====================================================================
  ! ---   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL
  ! =====================================================================

  ! initialization:
  DO i = 1, ncum*nlp
    hp(i, 1) = h(i, 1)
  END DO

  DO k = minorig + 1, nl
    DO i = 1, ncum
      IF ((k>=icb(i)) .AND. (k<=inb(i))) THEN
        hp(i, k) = h(i, nk(i)) + (lv(i,k)+(cpd-cpv)*t(i,k))*ep(i, k)*clw(i, k &
          )
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE cv_undilute2

SUBROUTINE cv_closure(nloc, ncum, nd, nk, icb, tv, tvp, p, ph, dph, plcl, &
    cpn, iflag, cbmf)
  IMPLICIT NONE

  ! inputs:
  INTEGER ncum, nd, nloc
  INTEGER nk(nloc), icb(nloc)
  REAL tv(nloc, nd), tvp(nloc, nd), p(nloc, nd), dph(nloc, nd)
  REAL ph(nloc, nd+1) ! caution nd instead ndp1 to be consistent...
  REAL plcl(nloc), cpn(nloc, nd)

  ! outputs:
  INTEGER iflag(nloc)
  REAL cbmf(nloc) ! also an input

  ! local variables:
  INTEGER i, k, icbmax
  REAL dtpbl(nloc), dtmin(nloc), tvpplcl(nloc), tvaplcl(nloc)
  REAL work(nloc)

  include "cvthermo.h"
  include "cvparam.h"

  ! -------------------------------------------------------------------
  ! Compute icbmax.
  ! -------------------------------------------------------------------

  icbmax = 2
  DO i = 1, ncum
    icbmax = max(icbmax, icb(i))
  END DO

  ! =====================================================================
  ! ---  CALCULATE CLOUD BASE MASS FLUX
  ! =====================================================================

  ! tvpplcl = parcel temperature lifted adiabatically from level
  ! icb-1 to the LCL.
  ! tvaplcl = virtual temperature at the LCL.

  DO i = 1, ncum
    dtpbl(i) = 0.0
    tvpplcl(i) = tvp(i, icb(i)-1) - rrd*tvp(i, icb(i)-1)*(p(i,icb(i)-1)-plcl( &
      i))/(cpn(i,icb(i)-1)*p(i,icb(i)-1))
    tvaplcl(i) = tv(i, icb(i)) + (tvp(i,icb(i))-tvp(i,icb(i)+1))*(plcl(i)-p(i &
      ,icb(i)))/(p(i,icb(i))-p(i,icb(i)+1))
  END DO

  ! -------------------------------------------------------------------
  ! --- Interpolate difference between lifted parcel and
  ! --- environmental temperatures to lifted condensation level
  ! -------------------------------------------------------------------

  ! dtpbl = average of tvp-tv in the PBL (k=nk to icb-1).

  DO k = minorig, icbmax
    DO i = 1, ncum
      IF ((k>=nk(i)) .AND. (k<=(icb(i)-1))) THEN
        dtpbl(i) = dtpbl(i) + (tvp(i,k)-tv(i,k))*dph(i, k)
      END IF
    END DO
  END DO
  DO i = 1, ncum
    dtpbl(i) = dtpbl(i)/(ph(i,nk(i))-ph(i,icb(i)))
    dtmin(i) = tvpplcl(i) - tvaplcl(i) + dtmax + dtpbl(i)
  END DO

  ! -------------------------------------------------------------------
  ! --- Adjust cloud base mass flux
  ! -------------------------------------------------------------------

  DO i = 1, ncum
    work(i) = cbmf(i)
    cbmf(i) = max(0.0, (1.0-damp)*cbmf(i)+0.1*alpha*dtmin(i))
    IF ((work(i)==0.0) .AND. (cbmf(i)==0.0)) THEN
      iflag(i) = 3
    END IF
  END DO

  RETURN
END SUBROUTINE cv_closure

SUBROUTINE cv_mixing(nloc, ncum, nd, icb, nk, inb, inb1, ph, t, q, qs, u, v, &
    h, lv, qnk, hp, tv, tvp, ep, clw, cbmf, m, ment, qent, uent, vent, nent, &
    sij, elij)
  IMPLICIT NONE

  include "cvthermo.h"
  include "cvparam.h"

  ! inputs:
  INTEGER ncum, nd, nloc
  INTEGER icb(nloc), inb(nloc), inb1(nloc), nk(nloc)
  REAL cbmf(nloc), qnk(nloc)
  REAL ph(nloc, nd+1)
  REAL t(nloc, nd), q(nloc, nd), qs(nloc, nd), lv(nloc, nd)
  REAL u(nloc, nd), v(nloc, nd), h(nloc, nd), hp(nloc, nd)
  REAL tv(nloc, nd), tvp(nloc, nd), ep(nloc, nd), clw(nloc, nd)

  ! outputs:
  INTEGER nent(nloc, nd)
  REAL m(nloc, nd), ment(nloc, nd, nd), qent(nloc, nd, nd)
  REAL uent(nloc, nd, nd), vent(nloc, nd, nd)
  REAL sij(nloc, nd, nd), elij(nloc, nd, nd)

  ! local variables:
  INTEGER i, j, k, ij
  INTEGER num1, num2
  REAL dbo, qti, bf2, anum, denom, dei, altem, cwat, stemp
  REAL alt, qp1, smid, sjmin, sjmax, delp, delm
  REAL work(nloc), asij(nloc), smin(nloc), scrit(nloc)
  REAL bsum(nloc, nd)
  LOGICAL lwork(nloc)

  ! =====================================================================
  ! --- INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
  ! =====================================================================

  DO i = 1, ncum*nlp
    nent(i, 1) = 0
    m(i, 1) = 0.0
  END DO

  DO k = 1, nlp
    DO j = 1, nlp
      DO i = 1, ncum
        qent(i, k, j) = q(i, j)
        uent(i, k, j) = u(i, j)
        vent(i, k, j) = v(i, j)
        elij(i, k, j) = 0.0
        ment(i, k, j) = 0.0
        sij(i, k, j) = 0.0
      END DO
    END DO
  END DO

  ! -------------------------------------------------------------------
  ! --- Calculate rates of mixing,  m(i)
  ! -------------------------------------------------------------------

  CALL zilch(work, ncum)

  DO j = minorig + 1, nl
    DO i = 1, ncum
      IF ((j>=(icb(i)+1)) .AND. (j<=inb(i))) THEN
        k = min(j, inb1(i))
        dbo = abs(tv(i,k+1)-tvp(i,k+1)-tv(i,k-1)+tvp(i,k-1)) + &
          entp*0.04*(ph(i,k)-ph(i,k+1))
        work(i) = work(i) + dbo
        m(i, j) = cbmf(i)*dbo
      END IF
    END DO
  END DO
  DO k = minorig + 1, nl
    DO i = 1, ncum
      IF ((k>=(icb(i)+1)) .AND. (k<=inb(i))) THEN
        m(i, k) = m(i, k)/work(i)
      END IF
    END DO
  END DO


  ! =====================================================================
  ! --- CALCULATE ENTRAINED AIR MASS FLUX (ment), TOTAL WATER MIXING
  ! --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
  ! --- FRACTION (sij)
  ! =====================================================================


  DO i = minorig + 1, nl
    DO j = minorig + 1, nl
      DO ij = 1, ncum
        IF ((i>=(icb(ij)+1)) .AND. (j>=icb(ij)) .AND. (i<=inb(ij)) .AND. (j<= &
            inb(ij))) THEN
          qti = qnk(ij) - ep(ij, i)*clw(ij, i)
          bf2 = 1. + lv(ij, j)*lv(ij, j)*qs(ij, j)/(rrv*t(ij,j)*t(ij,j)*cpd)
          anum = h(ij, j) - hp(ij, i) + (cpv-cpd)*t(ij, j)*(qti-q(ij,j))
          denom = h(ij, i) - hp(ij, i) + (cpd-cpv)*(q(ij,i)-qti)*t(ij, j)
          dei = denom
          IF (abs(dei)<0.01) dei = 0.01
          sij(ij, i, j) = anum/dei
          sij(ij, i, i) = 1.0
          altem = sij(ij, i, j)*q(ij, i) + (1.-sij(ij,i,j))*qti - qs(ij, j)
          altem = altem/bf2
          cwat = clw(ij, j)*(1.-ep(ij,j))
          stemp = sij(ij, i, j)
          IF ((stemp<0.0 .OR. stemp>1.0 .OR. altem>cwat) .AND. j>i) THEN
            anum = anum - lv(ij, j)*(qti-qs(ij,j)-cwat*bf2)
            denom = denom + lv(ij, j)*(q(ij,i)-qti)
            IF (abs(denom)<0.01) denom = 0.01
            sij(ij, i, j) = anum/denom
            altem = sij(ij, i, j)*q(ij, i) + (1.-sij(ij,i,j))*qti - qs(ij, j)
            altem = altem - (bf2-1.)*cwat
          END IF
          IF (sij(ij,i,j)>0.0 .AND. sij(ij,i,j)<0.9) THEN
            qent(ij, i, j) = sij(ij, i, j)*q(ij, i) + (1.-sij(ij,i,j))*qti
            uent(ij, i, j) = sij(ij, i, j)*u(ij, i) + &
              (1.-sij(ij,i,j))*u(ij, nk(ij))
            vent(ij, i, j) = sij(ij, i, j)*v(ij, i) + &
              (1.-sij(ij,i,j))*v(ij, nk(ij))
            elij(ij, i, j) = altem
            elij(ij, i, j) = max(0.0, elij(ij,i,j))
            ment(ij, i, j) = m(ij, i)/(1.-sij(ij,i,j))
            nent(ij, i) = nent(ij, i) + 1
          END IF
          sij(ij, i, j) = max(0.0, sij(ij,i,j))
          sij(ij, i, j) = min(1.0, sij(ij,i,j))
        END IF
      END DO
    END DO

    ! ***   If no air can entrain at level i assume that updraft detrains
    ! ***
    ! ***   at that level and calculate detrained air flux and properties
    ! ***

    DO ij = 1, ncum
      IF ((i>=(icb(ij)+1)) .AND. (i<=inb(ij)) .AND. (nent(ij,i)==0)) THEN
        ment(ij, i, i) = m(ij, i)
        qent(ij, i, i) = q(ij, nk(ij)) - ep(ij, i)*clw(ij, i)
        uent(ij, i, i) = u(ij, nk(ij))
        vent(ij, i, i) = v(ij, nk(ij))
        elij(ij, i, i) = clw(ij, i)
        sij(ij, i, i) = 1.0
      END IF
    END DO
  END DO

  DO i = 1, ncum
    sij(i, inb(i), inb(i)) = 1.0
  END DO

  ! =====================================================================
  ! ---  NORMALIZE ENTRAINED AIR MASS FLUXES
  ! ---  TO REPRESENT EQUAL PROBABILITIES OF MIXING
  ! =====================================================================

  CALL zilch(bsum, ncum*nlp)
  DO ij = 1, ncum
    lwork(ij) = .FALSE.
  END DO
  DO i = minorig + 1, nl

    num1 = 0
    DO ij = 1, ncum
      IF ((i>=icb(ij)+1) .AND. (i<=inb(ij))) num1 = num1 + 1
    END DO
    IF (num1<=0) GO TO 789

    DO ij = 1, ncum
      IF ((i>=icb(ij)+1) .AND. (i<=inb(ij))) THEN
        lwork(ij) = (nent(ij,i)/=0)
        qp1 = q(ij, nk(ij)) - ep(ij, i)*clw(ij, i)
        anum = h(ij, i) - hp(ij, i) - lv(ij, i)*(qp1-qs(ij,i))
        denom = h(ij, i) - hp(ij, i) + lv(ij, i)*(q(ij,i)-qp1)
        IF (abs(denom)<0.01) denom = 0.01
        scrit(ij) = anum/denom
        alt = qp1 - qs(ij, i) + scrit(ij)*(q(ij,i)-qp1)
        IF (scrit(ij)<0.0 .OR. alt<0.0) scrit(ij) = 1.0
        asij(ij) = 0.0
        smin(ij) = 1.0
      END IF
    END DO
    DO j = minorig, nl

      num2 = 0
      DO ij = 1, ncum
        IF ((i>=icb(ij)+1) .AND. (i<=inb(ij)) .AND. (j>=icb( &
          ij)) .AND. (j<=inb(ij)) .AND. lwork(ij)) num2 = num2 + 1
      END DO
      IF (num2<=0) GO TO 783

      DO ij = 1, ncum
        IF ((i>=icb(ij)+1) .AND. (i<=inb(ij)) .AND. (j>=icb( &
            ij)) .AND. (j<=inb(ij)) .AND. lwork(ij)) THEN
          IF (sij(ij,i,j)>0.0 .AND. sij(ij,i,j)<0.9) THEN
            IF (j>i) THEN
              smid = min(sij(ij,i,j), scrit(ij))
              sjmax = smid
              sjmin = smid
              IF (smid<smin(ij) .AND. sij(ij,i,j+1)<smid) THEN
                smin(ij) = smid
                sjmax = min(sij(ij,i,j+1), sij(ij,i,j), scrit(ij))
                sjmin = max(sij(ij,i,j-1), sij(ij,i,j))
                sjmin = min(sjmin, scrit(ij))
              END IF
            ELSE
              sjmax = max(sij(ij,i,j+1), scrit(ij))
              smid = max(sij(ij,i,j), scrit(ij))
              sjmin = 0.0
              IF (j>1) sjmin = sij(ij, i, j-1)
              sjmin = max(sjmin, scrit(ij))
            END IF
            delp = abs(sjmax-smid)
            delm = abs(sjmin-smid)
            asij(ij) = asij(ij) + (delp+delm)*(ph(ij,j)-ph(ij,j+1))
            ment(ij, i, j) = ment(ij, i, j)*(delp+delm)*(ph(ij,j)-ph(ij,j+1))
          END IF
        END IF
      END DO
783 END DO
    DO ij = 1, ncum
      IF ((i>=icb(ij)+1) .AND. (i<=inb(ij)) .AND. lwork(ij)) THEN
        asij(ij) = max(1.0E-21, asij(ij))
        asij(ij) = 1.0/asij(ij)
        bsum(ij, i) = 0.0
      END IF
    END DO
    DO j = minorig, nl + 1
      DO ij = 1, ncum
        IF ((i>=icb(ij)+1) .AND. (i<=inb(ij)) .AND. (j>=icb( &
            ij)) .AND. (j<=inb(ij)) .AND. lwork(ij)) THEN
          ment(ij, i, j) = ment(ij, i, j)*asij(ij)
          bsum(ij, i) = bsum(ij, i) + ment(ij, i, j)
        END IF
      END DO
    END DO
    DO ij = 1, ncum
      IF ((i>=icb(ij)+1) .AND. (i<=inb(ij)) .AND. (bsum(ij, &
          i)<1.0E-18) .AND. lwork(ij)) THEN
        nent(ij, i) = 0
        ment(ij, i, i) = m(ij, i)
        qent(ij, i, i) = q(ij, nk(ij)) - ep(ij, i)*clw(ij, i)
        uent(ij, i, i) = u(ij, nk(ij))
        vent(ij, i, i) = v(ij, nk(ij))
        elij(ij, i, i) = clw(ij, i)
        sij(ij, i, i) = 1.0
      END IF
    END DO
789 END DO

  RETURN
END SUBROUTINE cv_mixing

SUBROUTINE cv_unsat(nloc, ncum, nd, inb, t, q, qs, gz, u, v, p, ph, h, lv, &
    ep, sigp, clw, m, ment, elij, iflag, mp, qp, up, vp, wt, water, evap)
  IMPLICIT NONE


  include "cvthermo.h"
  include "cvparam.h"

  ! inputs:
  INTEGER ncum, nd, nloc
  INTEGER inb(nloc)
  REAL t(nloc, nd), q(nloc, nd), qs(nloc, nd)
  REAL gz(nloc, nd), u(nloc, nd), v(nloc, nd)
  REAL p(nloc, nd), ph(nloc, nd+1), h(nloc, nd)
  REAL lv(nloc, nd), ep(nloc, nd), sigp(nloc, nd), clw(nloc, nd)
  REAL m(nloc, nd), ment(nloc, nd, nd), elij(nloc, nd, nd)

  ! outputs:
  INTEGER iflag(nloc) ! also an input
  REAL mp(nloc, nd), qp(nloc, nd), up(nloc, nd), vp(nloc, nd)
  REAL water(nloc, nd), evap(nloc, nd), wt(nloc, nd)

  ! local variables:
  INTEGER i, j, k, ij, num1
  INTEGER jtt(nloc)
  REAL awat, coeff, qsm, afac, sigt, b6, c6, revap
  REAL dhdp, fac, qstm, rat
  REAL wdtrain(nloc)
  LOGICAL lwork(nloc)

  ! =====================================================================
  ! --- PRECIPITATING DOWNDRAFT CALCULATION
  ! =====================================================================

  ! Initializations:

  DO i = 1, ncum
    DO k = 1, nl + 1
      wt(i, k) = omtsnow
      mp(i, k) = 0.0
      evap(i, k) = 0.0
      water(i, k) = 0.0
    END DO
  END DO

  DO i = 1, ncum
    qp(i, 1) = q(i, 1)
    up(i, 1) = u(i, 1)
    vp(i, 1) = v(i, 1)
  END DO

  DO k = 2, nl + 1
    DO i = 1, ncum
      qp(i, k) = q(i, k-1)
      up(i, k) = u(i, k-1)
      vp(i, k) = v(i, k-1)
    END DO
  END DO


  ! ***  Check whether ep(inb)=0, if so, skip precipitating    ***
  ! ***             downdraft calculation                      ***


  ! ***  Integrate liquid water equation to find condensed water   ***
  ! ***                and condensed water flux                    ***


  DO i = 1, ncum
    jtt(i) = 2
    IF (ep(i,inb(i))<=0.0001) iflag(i) = 2
    IF (iflag(i)==0) THEN
      lwork(i) = .TRUE.
    ELSE
      lwork(i) = .FALSE.
    END IF
  END DO

  ! ***                    Begin downdraft loop                    ***


  CALL zilch(wdtrain, ncum)
  DO i = nl + 1, 1, -1

    num1 = 0
    DO ij = 1, ncum
      IF ((i<=inb(ij)) .AND. lwork(ij)) num1 = num1 + 1
    END DO
    IF (num1<=0) GO TO 899


    ! ***        Calculate detrained precipitation             ***

    DO ij = 1, ncum
      IF ((i<=inb(ij)) .AND. (lwork(ij))) THEN
        wdtrain(ij) = g*ep(ij, i)*m(ij, i)*clw(ij, i)
      END IF
    END DO

    IF (i>1) THEN
      DO j = 1, i - 1
        DO ij = 1, ncum
          IF ((i<=inb(ij)) .AND. (lwork(ij))) THEN
            awat = elij(ij, j, i) - (1.-ep(ij,i))*clw(ij, i)
            awat = max(0.0, awat)
            wdtrain(ij) = wdtrain(ij) + g*awat*ment(ij, j, i)
          END IF
        END DO
      END DO
    END IF

    ! ***    Find rain water and evaporation using provisional   ***
    ! ***              estimates of qp(i)and qp(i-1)             ***


    ! ***  Value of terminal velocity and coeffecient of evaporation for snow
    ! ***

    DO ij = 1, ncum
      IF ((i<=inb(ij)) .AND. (lwork(ij))) THEN
        coeff = coeffs
        wt(ij, i) = omtsnow

        ! ***  Value of terminal velocity and coeffecient of evaporation for
        ! rain   ***

        IF (t(ij,i)>273.0) THEN
          coeff = coeffr
          wt(ij, i) = omtrain
        END IF
        qsm = 0.5*(q(ij,i)+qp(ij,i+1))
        afac = coeff*ph(ij, i)*(qs(ij,i)-qsm)/(1.0E4+2.0E3*ph(ij,i)*qs(ij,i))
        afac = max(afac, 0.0)
        sigt = sigp(ij, i)
        sigt = max(0.0, sigt)
        sigt = min(1.0, sigt)
        b6 = 100.*(ph(ij,i)-ph(ij,i+1))*sigt*afac/wt(ij, i)
        c6 = (water(ij,i+1)*wt(ij,i+1)+wdtrain(ij)/sigd)/wt(ij, i)
        revap = 0.5*(-b6+sqrt(b6*b6+4.*c6))
        evap(ij, i) = sigt*afac*revap
        water(ij, i) = revap*revap

        ! ***  Calculate precipitating downdraft mass flux under     ***
        ! ***              hydrostatic approximation                 ***

        IF (i>1) THEN
          dhdp = (h(ij,i)-h(ij,i-1))/(p(ij,i-1)-p(ij,i))
          dhdp = max(dhdp, 10.0)
          mp(ij, i) = 100.*ginv*lv(ij, i)*sigd*evap(ij, i)/dhdp
          mp(ij, i) = max(mp(ij,i), 0.0)

          ! ***   Add small amount of inertia to downdraft              ***

          fac = 20.0/(ph(ij,i-1)-ph(ij,i))
          mp(ij, i) = (fac*mp(ij,i+1)+mp(ij,i))/(1.+fac)

          ! ***      Force mp to decrease linearly to zero
          ! ***
          ! ***      between about 950 mb and the surface
          ! ***

          IF (p(ij,i)>(0.949*p(ij,1))) THEN
            jtt(ij) = max(jtt(ij), i)
            mp(ij, i) = mp(ij, jtt(ij))*(p(ij,1)-p(ij,i))/ &
              (p(ij,1)-p(ij,jtt(ij)))
          END IF
        END IF

        ! ***       Find mixing ratio of precipitating downdraft     ***

        IF (i/=inb(ij)) THEN
          IF (i==1) THEN
            qstm = qs(ij, 1)
          ELSE
            qstm = qs(ij, i-1)
          END IF
          IF (mp(ij,i)>mp(ij,i+1)) THEN
            rat = mp(ij, i+1)/mp(ij, i)
            qp(ij, i) = qp(ij, i+1)*rat + q(ij, i)*(1.0-rat) + &
              100.*ginv*sigd*(ph(ij,i)-ph(ij,i+1))*(evap(ij,i)/mp(ij,i))
            up(ij, i) = up(ij, i+1)*rat + u(ij, i)*(1.-rat)
            vp(ij, i) = vp(ij, i+1)*rat + v(ij, i)*(1.-rat)
          ELSE
            IF (mp(ij,i+1)>0.0) THEN
              qp(ij, i) = (gz(ij,i+1)-gz(ij,i)+qp(ij,i+1)*(lv(ij,i+1)+t(ij, &
                i+1)*(cl-cpd))+cpd*(t(ij,i+1)-t(ij, &
                i)))/(lv(ij,i)+t(ij,i)*(cl-cpd))
              up(ij, i) = up(ij, i+1)
              vp(ij, i) = vp(ij, i+1)
            END IF
          END IF
          qp(ij, i) = min(qp(ij,i), qstm)
          qp(ij, i) = max(qp(ij,i), 0.0)
        END IF
      END IF
    END DO
899 END DO

  RETURN
END SUBROUTINE cv_unsat

SUBROUTINE cv_yield(nloc, ncum, nd, nk, icb, inb, delt, t, q, u, v, gz, p, &
    ph, h, hp, lv, cpn, ep, clw, frac, m, mp, qp, up, vp, wt, water, evap, &
    ment, qent, uent, vent, nent, elij, tv, tvp, iflag, wd, qprime, tprime, &
    precip, cbmf, ft, fq, fu, fv, ma, qcondc)
  IMPLICIT NONE

  include "cvthermo.h"
  include "cvparam.h"

  ! inputs
  INTEGER ncum, nd, nloc
  INTEGER nk(nloc), icb(nloc), inb(nloc)
  INTEGER nent(nloc, nd)
  REAL delt
  REAL t(nloc, nd), q(nloc, nd), u(nloc, nd), v(nloc, nd)
  REAL gz(nloc, nd)
  REAL p(nloc, nd), ph(nloc, nd+1), h(nloc, nd)
  REAL hp(nloc, nd), lv(nloc, nd)
  REAL cpn(nloc, nd), ep(nloc, nd), clw(nloc, nd), frac(nloc)
  REAL m(nloc, nd), mp(nloc, nd), qp(nloc, nd)
  REAL up(nloc, nd), vp(nloc, nd)
  REAL wt(nloc, nd), water(nloc, nd), evap(nloc, nd)
  REAL ment(nloc, nd, nd), qent(nloc, nd, nd), elij(nloc, nd, nd)
  REAL uent(nloc, nd, nd), vent(nloc, nd, nd)
  REAL tv(nloc, nd), tvp(nloc, nd)

  ! outputs
  INTEGER iflag(nloc) ! also an input
  REAL cbmf(nloc) ! also an input
  REAL wd(nloc), tprime(nloc), qprime(nloc)
  REAL precip(nloc)
  REAL ft(nloc, nd), fq(nloc, nd), fu(nloc, nd), fv(nloc, nd)
  REAL ma(nloc, nd)
  REAL qcondc(nloc, nd)

  ! local variables
  INTEGER i, j, ij, k, num1
  REAL dpinv, cpinv, awat, fqold, ftold, fuold, fvold, delti
  REAL work(nloc), am(nloc), amp1(nloc), ad(nloc)
  REAL ents(nloc), uav(nloc), vav(nloc), lvcp(nloc, nd)
  REAL qcond(nloc, nd), nqcond(nloc, nd), wa(nloc, nd) ! cld
  REAL siga(nloc, nd), ax(nloc, nd), mac(nloc, nd) ! cld


  ! -- initializations:

  delti = 1.0/delt

  DO i = 1, ncum
    precip(i) = 0.0
    wd(i) = 0.0
    tprime(i) = 0.0
    qprime(i) = 0.0
    DO k = 1, nl + 1
      ft(i, k) = 0.0
      fu(i, k) = 0.0
      fv(i, k) = 0.0
      fq(i, k) = 0.0
      lvcp(i, k) = lv(i, k)/cpn(i, k)
      qcondc(i, k) = 0.0 ! cld
      qcond(i, k) = 0.0 ! cld
      nqcond(i, k) = 0.0 ! cld
    END DO
  END DO


  ! ***  Calculate surface precipitation in mm/day     ***

  DO i = 1, ncum
    IF (iflag(i)<=1) THEN
      ! c            precip(i)=precip(i)+wt(i,1)*sigd*water(i,1)*3600.*24000.
      ! c     &                /(rowl*g)
      ! c            precip(i)=precip(i)*delt/86400.
      precip(i) = wt(i, 1)*sigd*water(i, 1)*86400/g
    END IF
  END DO


  ! ***  Calculate downdraft velocity scale and surface temperature and  ***
  ! ***                    water vapor fluctuations                      ***

  DO i = 1, ncum
    wd(i) = betad*abs(mp(i,icb(i)))*0.01*rrd*t(i, icb(i))/(sigd*p(i,icb(i)))
    qprime(i) = 0.5*(qp(i,1)-q(i,1))
    tprime(i) = lv0*qprime(i)/cpd
  END DO

  ! ***  Calculate tendencies of lowest level potential temperature  ***
  ! ***                      and mixing ratio                        ***

  DO i = 1, ncum
    work(i) = 0.01/(ph(i,1)-ph(i,2))
    am(i) = 0.0
  END DO
  DO k = 2, nl
    DO i = 1, ncum
      IF ((nk(i)==1) .AND. (k<=inb(i)) .AND. (nk(i)==1)) THEN
        am(i) = am(i) + m(i, k)
      END IF
    END DO
  END DO
  DO i = 1, ncum
    IF ((g*work(i)*am(i))>=delti) iflag(i) = 1
    ft(i, 1) = ft(i, 1) + g*work(i)*am(i)*(t(i,2)-t(i,1)+(gz(i,2)-gz(i, &
      1))/cpn(i,1))
    ft(i, 1) = ft(i, 1) - lvcp(i, 1)*sigd*evap(i, 1)
    ft(i, 1) = ft(i, 1) + sigd*wt(i, 2)*(cl-cpd)*water(i, 2)*(t(i,2)-t(i,1))* &
      work(i)/cpn(i, 1)
    fq(i, 1) = fq(i, 1) + g*mp(i, 2)*(qp(i,2)-q(i,1))*work(i) + &
      sigd*evap(i, 1)
    fq(i, 1) = fq(i, 1) + g*am(i)*(q(i,2)-q(i,1))*work(i)
    fu(i, 1) = fu(i, 1) + g*work(i)*(mp(i,2)*(up(i,2)-u(i,1))+am(i)*(u(i, &
      2)-u(i,1)))
    fv(i, 1) = fv(i, 1) + g*work(i)*(mp(i,2)*(vp(i,2)-v(i,1))+am(i)*(v(i, &
      2)-v(i,1)))
  END DO
  DO j = 2, nl
    DO i = 1, ncum
      IF (j<=inb(i)) THEN
        fq(i, 1) = fq(i, 1) + g*work(i)*ment(i, j, 1)*(qent(i,j,1)-q(i,1))
        fu(i, 1) = fu(i, 1) + g*work(i)*ment(i, j, 1)*(uent(i,j,1)-u(i,1))
        fv(i, 1) = fv(i, 1) + g*work(i)*ment(i, j, 1)*(vent(i,j,1)-v(i,1))
      END IF
    END DO
  END DO

  ! ***  Calculate tendencies of potential temperature and mixing ratio  ***
  ! ***               at levels above the lowest level                   ***

  ! ***  First find the net saturated updraft and downdraft mass fluxes  ***
  ! ***                      through each level                          ***

  DO i = 2, nl + 1

    num1 = 0
    DO ij = 1, ncum
      IF (i<=inb(ij)) num1 = num1 + 1
    END DO
    IF (num1<=0) GO TO 1500

    CALL zilch(amp1, ncum)
    CALL zilch(ad, ncum)

    DO k = i + 1, nl + 1
      DO ij = 1, ncum
        IF ((i>=nk(ij)) .AND. (i<=inb(ij)) .AND. (k<=(inb(ij)+1))) THEN
          amp1(ij) = amp1(ij) + m(ij, k)
        END IF
      END DO
    END DO

    DO k = 1, i
      DO j = i + 1, nl + 1
        DO ij = 1, ncum
          IF ((j<=(inb(ij)+1)) .AND. (i<=inb(ij))) THEN
            amp1(ij) = amp1(ij) + ment(ij, k, j)
          END IF
        END DO
      END DO
    END DO
    DO k = 1, i - 1
      DO j = i, nl + 1
        DO ij = 1, ncum
          IF ((i<=inb(ij)) .AND. (j<=inb(ij))) THEN
            ad(ij) = ad(ij) + ment(ij, j, k)
          END IF
        END DO
      END DO
    END DO

    DO ij = 1, ncum
      IF (i<=inb(ij)) THEN
        dpinv = 0.01/(ph(ij,i)-ph(ij,i+1))
        cpinv = 1.0/cpn(ij, i)

        ft(ij, i) = ft(ij, i) + g*dpinv*(amp1(ij)*(t(ij,i+1)-t(ij, &
          i)+(gz(ij,i+1)-gz(ij,i))*cpinv)-ad(ij)*(t(ij,i)-t(ij, &
          i-1)+(gz(ij,i)-gz(ij,i-1))*cpinv)) - sigd*lvcp(ij, i)*evap(ij, i)
        ft(ij, i) = ft(ij, i) + g*dpinv*ment(ij, i, i)*(hp(ij,i)-h(ij,i)+t(ij &
          ,i)*(cpv-cpd)*(q(ij,i)-qent(ij,i,i)))*cpinv
        ft(ij, i) = ft(ij, i) + sigd*wt(ij, i+1)*(cl-cpd)*water(ij, i+1)*(t( &
          ij,i+1)-t(ij,i))*dpinv*cpinv
        fq(ij, i) = fq(ij, i) + g*dpinv*(amp1(ij)*(q(ij,i+1)-q(ij, &
          i))-ad(ij)*(q(ij,i)-q(ij,i-1)))
        fu(ij, i) = fu(ij, i) + g*dpinv*(amp1(ij)*(u(ij,i+1)-u(ij, &
          i))-ad(ij)*(u(ij,i)-u(ij,i-1)))
        fv(ij, i) = fv(ij, i) + g*dpinv*(amp1(ij)*(v(ij,i+1)-v(ij, &
          i))-ad(ij)*(v(ij,i)-v(ij,i-1)))
      END IF
    END DO
    DO k = 1, i - 1
      DO ij = 1, ncum
        IF (i<=inb(ij)) THEN
          awat = elij(ij, k, i) - (1.-ep(ij,i))*clw(ij, i)
          awat = max(awat, 0.0)
          fq(ij, i) = fq(ij, i) + g*dpinv*ment(ij, k, i)*(qent(ij,k,i)-awat-q &
            (ij,i))
          fu(ij, i) = fu(ij, i) + g*dpinv*ment(ij, k, i)*(uent(ij,k,i)-u(ij,i &
            ))
          fv(ij, i) = fv(ij, i) + g*dpinv*ment(ij, k, i)*(vent(ij,k,i)-v(ij,i &
            ))
          ! (saturated updrafts resulting from mixing)               ! cld
          qcond(ij, i) = qcond(ij, i) + (elij(ij,k,i)-awat) ! cld
          nqcond(ij, i) = nqcond(ij, i) + 1. ! cld
        END IF
      END DO
    END DO
    DO k = i, nl + 1
      DO ij = 1, ncum
        IF ((i<=inb(ij)) .AND. (k<=inb(ij))) THEN
          fq(ij, i) = fq(ij, i) + g*dpinv*ment(ij, k, i)*(qent(ij,k,i)-q(ij,i &
            ))
          fu(ij, i) = fu(ij, i) + g*dpinv*ment(ij, k, i)*(uent(ij,k,i)-u(ij,i &
            ))
          fv(ij, i) = fv(ij, i) + g*dpinv*ment(ij, k, i)*(vent(ij,k,i)-v(ij,i &
            ))
        END IF
      END DO
    END DO
    DO ij = 1, ncum
      IF (i<=inb(ij)) THEN
        fq(ij, i) = fq(ij, i) + sigd*evap(ij, i) + g*(mp(ij,i+1)*(qp(ij, &
          i+1)-q(ij,i))-mp(ij,i)*(qp(ij,i)-q(ij,i-1)))*dpinv
        fu(ij, i) = fu(ij, i) + g*(mp(ij,i+1)*(up(ij,i+1)-u(ij, &
          i))-mp(ij,i)*(up(ij,i)-u(ij,i-1)))*dpinv
        fv(ij, i) = fv(ij, i) + g*(mp(ij,i+1)*(vp(ij,i+1)-v(ij, &
          i))-mp(ij,i)*(vp(ij,i)-v(ij,i-1)))*dpinv
        ! (saturated downdrafts resulting from mixing)               ! cld
        DO k = i + 1, inb(ij) ! cld
          qcond(ij, i) = qcond(ij, i) + elij(ij, k, i) ! cld
          nqcond(ij, i) = nqcond(ij, i) + 1. ! cld
        END DO ! cld
        ! (particular case: no detraining level is found)            ! cld
        IF (nent(ij,i)==0) THEN ! cld
          qcond(ij, i) = qcond(ij, i) + (1.-ep(ij,i))*clw(ij, i) ! cld
          nqcond(ij, i) = nqcond(ij, i) + 1. ! cld
        END IF ! cld
        IF (nqcond(ij,i)/=0.) THEN ! cld
          qcond(ij, i) = qcond(ij, i)/nqcond(ij, i) ! cld
        END IF ! cld
      END IF
    END DO
1500 END DO

  ! *** Adjust tendencies at top of convection layer to reflect  ***
  ! ***       actual position of the level zero cape             ***

  DO ij = 1, ncum
    fqold = fq(ij, inb(ij))
    fq(ij, inb(ij)) = fq(ij, inb(ij))*(1.-frac(ij))
    fq(ij, inb(ij)-1) = fq(ij, inb(ij)-1) + frac(ij)*fqold*((ph(ij, &
      inb(ij))-ph(ij,inb(ij)+1))/(ph(ij,inb(ij)-1)-ph(ij, &
      inb(ij))))*lv(ij, inb(ij))/lv(ij, inb(ij)-1)
    ftold = ft(ij, inb(ij))
    ft(ij, inb(ij)) = ft(ij, inb(ij))*(1.-frac(ij))
    ft(ij, inb(ij)-1) = ft(ij, inb(ij)-1) + frac(ij)*ftold*((ph(ij, &
      inb(ij))-ph(ij,inb(ij)+1))/(ph(ij,inb(ij)-1)-ph(ij, &
      inb(ij))))*cpn(ij, inb(ij))/cpn(ij, inb(ij)-1)
    fuold = fu(ij, inb(ij))
    fu(ij, inb(ij)) = fu(ij, inb(ij))*(1.-frac(ij))
    fu(ij, inb(ij)-1) = fu(ij, inb(ij)-1) + frac(ij)*fuold*((ph(ij, &
      inb(ij))-ph(ij,inb(ij)+1))/(ph(ij,inb(ij)-1)-ph(ij,inb(ij))))
    fvold = fv(ij, inb(ij))
    fv(ij, inb(ij)) = fv(ij, inb(ij))*(1.-frac(ij))
    fv(ij, inb(ij)-1) = fv(ij, inb(ij)-1) + frac(ij)*fvold*((ph(ij, &
      inb(ij))-ph(ij,inb(ij)+1))/(ph(ij,inb(ij)-1)-ph(ij,inb(ij))))
  END DO

  ! ***   Very slightly adjust tendencies to force exact   ***
  ! ***     enthalpy, momentum and tracer conservation     ***

  DO ij = 1, ncum
    ents(ij) = 0.0
    uav(ij) = 0.0
    vav(ij) = 0.0
    DO i = 1, inb(ij)
      ents(ij) = ents(ij) + (cpn(ij,i)*ft(ij,i)+lv(ij,i)*fq(ij,i))*(ph(ij,i)- &
        ph(ij,i+1))
      uav(ij) = uav(ij) + fu(ij, i)*(ph(ij,i)-ph(ij,i+1))
      vav(ij) = vav(ij) + fv(ij, i)*(ph(ij,i)-ph(ij,i+1))
    END DO
  END DO
  DO ij = 1, ncum
    ents(ij) = ents(ij)/(ph(ij,1)-ph(ij,inb(ij)+1))
    uav(ij) = uav(ij)/(ph(ij,1)-ph(ij,inb(ij)+1))
    vav(ij) = vav(ij)/(ph(ij,1)-ph(ij,inb(ij)+1))
  END DO
  DO ij = 1, ncum
    DO i = 1, inb(ij)
      ft(ij, i) = ft(ij, i) - ents(ij)/cpn(ij, i)
      fu(ij, i) = (1.-cu)*(fu(ij,i)-uav(ij))
      fv(ij, i) = (1.-cu)*(fv(ij,i)-vav(ij))
    END DO
  END DO

  DO k = 1, nl + 1
    DO i = 1, ncum
      IF ((q(i,k)+delt*fq(i,k))<0.0) iflag(i) = 10
    END DO
  END DO


  DO i = 1, ncum
    IF (iflag(i)>2) THEN
      precip(i) = 0.0
      cbmf(i) = 0.0
    END IF
  END DO
  DO k = 1, nl
    DO i = 1, ncum
      IF (iflag(i)>2) THEN
        ft(i, k) = 0.0
        fq(i, k) = 0.0
        fu(i, k) = 0.0
        fv(i, k) = 0.0
        qcondc(i, k) = 0.0 ! cld
      END IF
    END DO
  END DO

  DO k = 1, nl + 1
    DO i = 1, ncum
      ma(i, k) = 0.
    END DO
  END DO
  DO k = nl, 1, -1
    DO i = 1, ncum
      ma(i, k) = ma(i, k+1) + m(i, k)
    END DO
  END DO


  ! *** diagnose the in-cloud mixing ratio   ***            ! cld
  ! ***           of condensed water         ***            ! cld
  ! ! cld
  DO ij = 1, ncum ! cld
    DO i = 1, nd ! cld
      mac(ij, i) = 0.0 ! cld
      wa(ij, i) = 0.0 ! cld
      siga(ij, i) = 0.0 ! cld
    END DO ! cld
    DO i = nk(ij), inb(ij) ! cld
      DO k = i + 1, inb(ij) + 1 ! cld
        mac(ij, i) = mac(ij, i) + m(ij, k) ! cld
      END DO ! cld
    END DO ! cld
    DO i = icb(ij), inb(ij) - 1 ! cld
      ax(ij, i) = 0. ! cld
      DO j = icb(ij), i ! cld
        ax(ij, i) = ax(ij, i) + rrd*(tvp(ij,j)-tv(ij,j)) & ! cld
          *(ph(ij,j)-ph(ij,j+1))/p(ij, j) ! cld
      END DO ! cld
      IF (ax(ij,i)>0.0) THEN ! cld
        wa(ij, i) = sqrt(2.*ax(ij,i)) ! cld
      END IF ! cld
    END DO ! cld
    DO i = 1, nl ! cld
      IF (wa(ij,i)>0.0) &          ! cld
        siga(ij, i) = mac(ij, i)/wa(ij, i) & ! cld
        *rrd*tvp(ij, i)/p(ij, i)/100./delta ! cld
      siga(ij, i) = min(siga(ij,i), 1.0) ! cld
      qcondc(ij, i) = siga(ij, i)*clw(ij, i)*(1.-ep(ij,i)) & ! cld
        +(1.-siga(ij,i))*qcond(ij, i) ! cld
    END DO ! cld
  END DO ! cld

  RETURN
END SUBROUTINE cv_yield

SUBROUTINE cv_uncompress(nloc, len, ncum, nd, idcum, iflag, precip, cbmf, ft, &
    fq, fu, fv, ma, qcondc, iflag1, precip1, cbmf1, ft1, fq1, fu1, fv1, ma1, &
    qcondc1)
  IMPLICIT NONE

  include "cvparam.h"

  ! inputs:
  INTEGER len, ncum, nd, nloc
  INTEGER idcum(nloc)
  INTEGER iflag(nloc)
  REAL precip(nloc), cbmf(nloc)
  REAL ft(nloc, nd), fq(nloc, nd), fu(nloc, nd), fv(nloc, nd)
  REAL ma(nloc, nd)
  REAL qcondc(nloc, nd) !cld

  ! outputs:
  INTEGER iflag1(len)
  REAL precip1(len), cbmf1(len)
  REAL ft1(len, nd), fq1(len, nd), fu1(len, nd), fv1(len, nd)
  REAL ma1(len, nd)
  REAL qcondc1(len, nd) !cld

  ! local variables:
  INTEGER i, k

  DO i = 1, ncum
    precip1(idcum(i)) = precip(i)
    cbmf1(idcum(i)) = cbmf(i)
    iflag1(idcum(i)) = iflag(i)
  END DO

  DO k = 1, nl
    DO i = 1, ncum
      ft1(idcum(i), k) = ft(i, k)
      fq1(idcum(i), k) = fq(i, k)
      fu1(idcum(i), k) = fu(i, k)
      fv1(idcum(i), k) = fv(i, k)
      ma1(idcum(i), k) = ma(i, k)
      qcondc1(idcum(i), k) = qcondc(i, k)
    END DO
  END DO

  RETURN
END SUBROUTINE cv_uncompress

