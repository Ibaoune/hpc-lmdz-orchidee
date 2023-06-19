
! $Header$

SUBROUTINE tlift43(p, t, q, qs, gz, icb, nk, tvp, tpk, clw, nd, nl, kk)
  IMPLICIT NONE
  REAL gz(nd), tpk(nd), clw(nd), p(nd)
  REAL t(nd), q(nd), qs(nd), tvp(nd), lv0
  INTEGER icb, nk, nd, nl, kk
  REAL cpd, cpv,  cl, g, rowl, gravity, cpvmcl, eps, epsi
  REAL ah0, cpp, cpinv, tg, qg, alv, s, ahg, tc, denom, es
  INTEGER i, nst, nsb, j 
  ! ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***

  ! -- sb:
  ! !      CPD=1005.7
  ! !      CPV=1870.0
  ! !      CL=4190.0
  ! !      RV=461.5
  ! !      RD=287.04
  ! !      LV0=2.501E6
  ! !      G=9.8
  ! !      ROWL=1000.0
  ! ajouts:
  include "YOMCST.h"
  cpd = rcpd
  cpv = rcpv
  cl = rcw
  lv0 = rlvtt
  g = rg
  rowl = ratm/100.
  gravity = rg !sb: Pr que gravite ne devienne pas humidite!
  ! sb --

  cpvmcl = cl - cpv
  eps = rd/rv
  epsi = 1./eps

  ! ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***

  ah0 = (cpd*(1.-q(nk))+cl*q(nk))*t(nk) + q(nk)*(lv0-cpvmcl*(t(nk)-273.15)) + &
    gz(nk)
  cpp = cpd*(1.-q(nk)) + q(nk)*cpv
  cpinv = 1./cpp

  IF (kk==1) THEN

    ! ***   CALCULATE LIFTED PARCEL QUANTITIES BELOW CLOUD BASE   ***

    DO i = 1, icb - 1
      clw(i) = 0.0
    END DO
    DO i = nk, icb - 1
      tpk(i) = t(nk) - (gz(i)-gz(nk))*cpinv
      tvp(i) = tpk(i)*(1.+q(nk)*epsi)
    END DO
  END IF

  ! ***  FIND LIFTED PARCEL QUANTITIES ABOVE CLOUD BASE    ***

  nst = icb
  nsb = icb
  IF (kk==2) THEN
    nst = nl
    nsb = icb + 1
  END IF
  DO i = nsb, nst
    tg = t(i)
    qg = qs(i)
    alv = lv0 - cpvmcl*(t(i)-273.15)
    DO j = 1, 2
      s = cpd + alv*alv*qg/(rv*t(i)*t(i))
      s = 1./s
      ahg = cpd*tg + (cl-cpd)*q(nk)*t(i) + alv*qg + gz(i)
      tg = tg + s*(ah0-ahg)
      tg = max(tg, 35.0)
      tc = tg - 273.15
      denom = 243.5 + tc
      IF (tc>=0.0) THEN
        es = 6.112*exp(17.67*tc/denom)
      ELSE
        es = exp(23.33086-6111.72784/tg+0.15215*log(tg))
      END IF
      qg = eps*es/(p(i)-es*(1.-eps))
    END DO
    alv = lv0 - cpvmcl*(t(i)-273.15)
    tpk(i) = (ah0-(cl-cpd)*q(nk)*t(i)-gz(i)-alv*qg)/cpd
    clw(i) = q(nk) - qg
    clw(i) = max(0.0, clw(i))
    rg = qg/(1.-q(nk))
    tvp(i) = tpk(i)*(1.+rg*epsi)
  END DO

  ! -- sb:
  rg = gravity ! RG redevient la gravite de YOMCST (sb)
  ! sb --

  RETURN
END SUBROUTINE tlift43

