
! $Id: cv30_routines.F90 3577 2019-10-07 13:30:19Z fhourdin $



SUBROUTINE cv30_param(nd, delt)
  IMPLICIT NONE

  ! ------------------------------------------------------------
  ! Set parameters for convectL for iflag_con = 3
  ! ------------------------------------------------------------


  ! ***  PBCRIT IS THE CRITICAL CLOUD DEPTH (MB) BENEATH WHICH THE ***
  ! ***      PRECIPITATION EFFICIENCY IS ASSUMED TO BE ZERO     ***
  ! ***  PTCRIT IS THE CLOUD DEPTH (MB) ABOVE WHICH THE PRECIP. ***
  ! ***            EFFICIENCY IS ASSUMED TO BE UNITY            ***
  ! ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
  ! ***  SPFAC IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE ***
  ! ***                        OF CLOUD                         ***

  ! [TAU: CHARACTERISTIC TIMESCALE USED TO COMPUTE ALPHA & BETA]
  ! ***    ALPHA AND BETA ARE PARAMETERS THAT CONTROL THE RATE OF ***
  ! ***                 APPROACH TO QUASI-EQUILIBRIUM           ***
  ! ***    (THEIR STANDARD VALUES ARE 1.0 AND 0.96, RESPECTIVELY) ***
  ! ***           (BETA MUST BE LESS THAN OR EQUAL TO 1)        ***

  ! ***    DTCRIT IS THE CRITICAL BUOYANCY (K) USED TO ADJUST THE ***
  ! ***                 APPROACH TO QUASI-EQUILIBRIUM           ***
  ! ***                     IT MUST BE LESS THAN 0              ***

  include "cv30param.h"
  include "conema3.h"

  INTEGER nd
  REAL delt ! timestep (seconds)

  ! noff: integer limit for convection (nd-noff)
  ! minorig: First level of convection

  ! -- limit levels for convection:

  noff = 1
  minorig = 1
  nl = nd - noff
  nlp = nl + 1
  nlm = nl - 1

  ! -- "microphysical" parameters:

  sigd = 0.01
  spfac = 0.15
  pbcrit = 150.0
  ptcrit = 500.0
  ! IM cf. FH     epmax  = 0.993

  omtrain = 45.0 ! used also for snow (no disctinction rain/snow)

  ! -- misc:

  dtovsh = -0.2 ! dT for overshoot
  dpbase = -40. ! definition cloud base (400m above LCL)
  dttrig = 5. ! (loose) condition for triggering

  ! -- rate of approach to quasi-equilibrium:

  dtcrit = -2.0
  tau = 8000.
  beta = 1.0 - delt/tau
  alpha = 1.5E-3*delt/tau
  ! increase alpha to compensate W decrease:
  alpha = alpha*1.5

  ! -- interface cloud parameterization:

  delta = 0.01 ! cld

  ! -- interface with boundary-layer (gust factor): (sb)

  betad = 10.0 ! original value (from convect 4.3)

  RETURN
END SUBROUTINE cv30_param

SUBROUTINE cv30_prelim(len, nd, ndp1, t, q, p, ph, lv, cpn, tv, gz, h, hm, &
    th)
  IMPLICIT NONE

  ! =====================================================================
  ! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
  ! "ori": from convect4.3 (vectorized)
  ! "convect3": to be exactly consistent with convect3
  ! =====================================================================

  ! inputs:
  INTEGER len, nd, ndp1
  REAL t(len, nd), q(len, nd), p(len, nd), ph(len, ndp1)

  ! outputs:
  REAL lv(len, nd), cpn(len, nd), tv(len, nd)
  REAL gz(len, nd), h(len, nd), hm(len, nd)
  REAL th(len, nd)

  ! local variables:
  INTEGER k, i
  REAL rdcp
  REAL tvx, tvy ! convect3
  REAL cpx(len, nd)

  include "cvthermo.h"
  include "cv30param.h"


  ! ori      do 110 k=1,nlp
  DO k = 1, nl ! convect3
    DO i = 1, len
      ! debug          lv(i,k)= lv0-clmcpv*(t(i,k)-t0)
      lv(i, k) = lv0 - clmcpv*(t(i,k)-273.15)
      cpn(i, k) = cpd*(1.0-q(i,k)) + cpv*q(i, k)
      cpx(i, k) = cpd*(1.0-q(i,k)) + cl*q(i, k)
      ! ori          tv(i,k)=t(i,k)*(1.0+q(i,k)*epsim1)
      tv(i, k) = t(i, k)*(1.0+q(i,k)/eps-q(i,k))
      rdcp = (rrd*(1.-q(i,k))+q(i,k)*rrv)/cpn(i, k)
      th(i, k) = t(i, k)*(1000.0/p(i,k))**rdcp
    END DO
  END DO

  ! gz = phi at the full levels (same as p).

  DO i = 1, len
    gz(i, 1) = 0.0
  END DO
  ! ori      do 140 k=2,nlp
  DO k = 2, nl ! convect3
    DO i = 1, len
      tvx = t(i, k)*(1.+q(i,k)/eps-q(i,k)) !convect3
      tvy = t(i, k-1)*(1.+q(i,k-1)/eps-q(i,k-1)) !convect3
      gz(i, k) = gz(i, k-1) + 0.5*rrd*(tvx+tvy) & !convect3
        *(p(i,k-1)-p(i,k))/ph(i, k) !convect3

      ! ori         gz(i,k)=gz(i,k-1)+hrd*(tv(i,k-1)+tv(i,k))
      ! ori    &         *(p(i,k-1)-p(i,k))/ph(i,k)
    END DO
  END DO

  ! h  = phi + cpT (dry static energy).
  ! hm = phi + cp(T-Tbase)+Lq

  ! ori      do 170 k=1,nlp
  DO k = 1, nl ! convect3
    DO i = 1, len
      h(i, k) = gz(i, k) + cpn(i, k)*t(i, k)
      hm(i, k) = gz(i, k) + cpx(i, k)*(t(i,k)-t(i,1)) + lv(i, k)*q(i, k)
    END DO
  END DO

  RETURN
END SUBROUTINE cv30_prelim

SUBROUTINE cv30_feed(len, nd, t, q, qs, p, ph, hm, gz, nk, icb, icbmax, &
    iflag, tnk, qnk, gznk, plcl)
  IMPLICIT NONE

  ! ================================================================
  ! Purpose: CONVECTIVE FEED

  ! Main differences with cv_feed:
  ! - ph added in input
  ! - here, nk(i)=minorig
  ! - icb defined differently (plcl compared with ph instead of p)

  ! Main differences with convect3:
  ! - we do not compute dplcldt and dplcldr of CLIFT anymore
  ! - values iflag different (but tests identical)
  ! - A,B explicitely defined (!...)
  ! ================================================================

  include "cv30param.h"

  ! inputs:
  INTEGER len, nd
  REAL t(len, nd), q(len, nd), qs(len, nd), p(len, nd)
  REAL hm(len, nd), gz(len, nd)
  REAL ph(len, nd+1)

  ! outputs:
  INTEGER iflag(len), nk(len), icb(len), icbmax
  REAL tnk(len), qnk(len), gznk(len), plcl(len)

  ! local variables:
  INTEGER i, k
  INTEGER ihmin(len)
  REAL work(len)
  REAL pnk(len), qsnk(len), rh(len), chi(len)
  REAL a, b ! convect3
  ! ym
  plcl = 0.0
  ! @ !-------------------------------------------------------------------
  ! @ ! --- Find level of minimum moist static energy
  ! @ ! --- If level of minimum moist static energy coincides with
  ! @ ! --- or is lower than minimum allowable parcel origin level,
  ! @ ! --- set iflag to 6.
  ! @ !-------------------------------------------------------------------
  ! @
  ! @       do 180 i=1,len
  ! @        work(i)=1.0e12
  ! @        ihmin(i)=nl
  ! @  180  continue
  ! @       do 200 k=2,nlp
  ! @         do 190 i=1,len
  ! @          if((hm(i,k).lt.work(i)).and.
  ! @      &      (hm(i,k).lt.hm(i,k-1)))then
  ! @            work(i)=hm(i,k)
  ! @            ihmin(i)=k
  ! @          endif
  ! @  190    continue
  ! @  200  continue
  ! @       do 210 i=1,len
  ! @         ihmin(i)=min(ihmin(i),nlm)
  ! @         if(ihmin(i).le.minorig)then
  ! @           iflag(i)=6
  ! @         endif
  ! @  210  continue
  ! @ c
  ! @ !-------------------------------------------------------------------
  ! @ ! --- Find that model level below the level of minimum moist static
  ! @ ! --- energy that has the maximum value of moist static energy
  ! @ !-------------------------------------------------------------------
  ! @
  ! @       do 220 i=1,len
  ! @        work(i)=hm(i,minorig)
  ! @        nk(i)=minorig
  ! @  220  continue
  ! @       do 240 k=minorig+1,nl
  ! @         do 230 i=1,len
  ! @          if((hm(i,k).gt.work(i)).and.(k.le.ihmin(i)))then
  ! @            work(i)=hm(i,k)
  ! @            nk(i)=k
  ! @          endif
  ! @  230     continue
  ! @  240  continue

  ! -------------------------------------------------------------------
  ! --- Origin level of ascending parcels for convect3:
  ! -------------------------------------------------------------------

  DO i = 1, len
    nk(i) = minorig
  END DO

  ! -------------------------------------------------------------------
  ! --- Check whether parcel level temperature and specific humidity
  ! --- are reasonable
  ! -------------------------------------------------------------------
  DO i = 1, len
    IF (((t(i,nk(i))<250.0) .OR. (q(i,nk(i))<=0.0)) & ! @      &       .or.(
                                                      ! p(i,ihmin(i)).lt.400.0
                                                      ! )  )
      .AND. (iflag(i)==0)) iflag(i) = 7
  END DO
  ! -------------------------------------------------------------------
  ! --- Calculate lifted condensation level of air at parcel origin level
  ! --- (Within 0.2% of formula of Bolton, MON. WEA. REV.,1980)
  ! -------------------------------------------------------------------

  a = 1669.0 ! convect3
  b = 122.0 ! convect3

  DO i = 1, len

    IF (iflag(i)/=7) THEN ! modif sb Jun7th 2002

      tnk(i) = t(i, nk(i))
      qnk(i) = q(i, nk(i))
      gznk(i) = gz(i, nk(i))
      pnk(i) = p(i, nk(i))
      qsnk(i) = qs(i, nk(i))

      rh(i) = qnk(i)/qsnk(i)
      ! ori        rh(i)=min(1.0,rh(i)) ! removed for convect3
      ! ori        chi(i)=tnk(i)/(1669.0-122.0*rh(i)-tnk(i))
      chi(i) = tnk(i)/(a-b*rh(i)-tnk(i)) ! convect3
      plcl(i) = pnk(i)*(rh(i)**chi(i))
      IF (((plcl(i)<200.0) .OR. (plcl(i)>=2000.0)) .AND. (iflag(i)==0)) iflag &
        (i) = 8

    END IF ! iflag=7

  END DO

  ! -------------------------------------------------------------------
  ! --- Calculate first level above lcl (=icb)
  ! -------------------------------------------------------------------

  ! @      do 270 i=1,len
  ! @       icb(i)=nlm
  ! @ 270  continue
  ! @c
  ! @      do 290 k=minorig,nl
  ! @        do 280 i=1,len
  ! @          if((k.ge.(nk(i)+1)).and.(p(i,k).lt.plcl(i)))
  ! @     &    icb(i)=min(icb(i),k)
  ! @ 280    continue
  ! @ 290  continue
  ! @c
  ! @      do 300 i=1,len
  ! @        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
  ! @ 300  continue

  DO i = 1, len
    icb(i) = nlm
  END DO

  ! la modification consiste a comparer plcl a ph et non a p:
  ! icb est defini par :  ph(icb)<plcl<ph(icb-1)
  ! @      do 290 k=minorig,nl
  DO k = 3, nl - 1 ! modif pour que icb soit sup/egal a 2
    DO i = 1, len
      IF (ph(i,k)<plcl(i)) icb(i) = min(icb(i), k)
    END DO
  END DO

  DO i = 1, len
    ! @        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
    IF ((icb(i)==nlm) .AND. (iflag(i)==0)) iflag(i) = 9
  END DO

  DO i = 1, len
    icb(i) = icb(i) - 1 ! icb sup ou egal a 2
  END DO

  ! Compute icbmax.

  icbmax = 2
  DO i = 1, len
    ! !        icbmax=max(icbmax,icb(i))
    IF (iflag(i)<7) icbmax = max(icbmax, icb(i)) ! sb Jun7th02
  END DO

  RETURN
END SUBROUTINE cv30_feed

SUBROUTINE cv30_undilute1(len, nd, t, q, qs, gz, plcl, p, nk, icb, tp, tvp, &
    clw, icbs)
  IMPLICIT NONE

  ! ----------------------------------------------------------------
  ! Equivalent de TLIFT entre NK et ICB+1 inclus

  ! Differences with convect4:
  ! - specify plcl in input
  ! - icbs is the first level above LCL (may differ from icb)
  ! - in the iterations, used x(icbs) instead x(icb)
  ! - many minor differences in the iterations
  ! - tvp is computed in only one time
  ! - icbs: first level above Plcl (IMIN de TLIFT) in output
  ! - if icbs=icb, compute also tp(icb+1),tvp(icb+1) & clw(icb+1)
  ! ----------------------------------------------------------------

  include "cvthermo.h"
  include "cv30param.h"

  ! inputs:
  INTEGER len, nd
  INTEGER nk(len), icb(len)
  REAL t(len, nd), q(len, nd), qs(len, nd), gz(len, nd)
  REAL p(len, nd)
  REAL plcl(len) ! convect3

  ! outputs:
  REAL tp(len, nd), tvp(len, nd), clw(len, nd)

  ! local variables:
  INTEGER i, k
  INTEGER icb1(len), icbs(len), icbsmax2 ! convect3
  REAL tg, qg, alv, s, ahg, tc, denom, es, rg
  REAL ah0(len), cpp(len)
  REAL tnk(len), qnk(len), gznk(len), ticb(len), gzicb(len)
  REAL qsicb(len) ! convect3
  REAL cpinv(len) ! convect3

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
    ! ori        ticb(i)=t(i,icb(i))
    ! ori        gzicb(i)=gz(i,icb(i))
  END DO

  ! ***  Calculate certain parcel quantities, including static energy   ***

  DO i = 1, len
    ah0(i) = (cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i) + qnk(i)*(lv0-clmcpv*(tnk(i)- &
      273.15)) + gznk(i)
    cpp(i) = cpd*(1.-qnk(i)) + qnk(i)*cpv
    cpinv(i) = 1./cpp(i)
  END DO

  ! ***   Calculate lifted parcel quantities below cloud base   ***

  DO i = 1, len !convect3
    icb1(i) = min(max(icb(i), 2), nl)
    ! if icb is below LCL, start loop at ICB+1:
    ! (icbs est le premier niveau au-dessus du LCL)
    icbs(i) = icb1(i) !convect3
    IF (plcl(i)<p(i,icb1(i))) THEN
      icbs(i) = min(icbs(i)+1, nl) !convect3
    END IF
  END DO !convect3

  DO i = 1, len !convect3
    ticb(i) = t(i, icbs(i)) !convect3
    gzicb(i) = gz(i, icbs(i)) !convect3
    qsicb(i) = qs(i, icbs(i)) !convect3
  END DO !convect3


  ! Re-compute icbsmax (icbsmax2):        !convect3
  ! !convect3
  icbsmax2 = 2 !convect3
  DO i = 1, len !convect3
    icbsmax2 = max(icbsmax2, icbs(i)) !convect3
  END DO !convect3

  ! initialization outputs:

  DO k = 1, icbsmax2 ! convect3
    DO i = 1, len ! convect3
      tp(i, k) = 0.0 ! convect3
      tvp(i, k) = 0.0 ! convect3
      clw(i, k) = 0.0 ! convect3
    END DO ! convect3
  END DO ! convect3

  ! tp and tvp below cloud base:

  DO k = minorig, icbsmax2 - 1
    DO i = 1, len
      tp(i, k) = tnk(i) - (gz(i,k)-gznk(i))*cpinv(i)
      tvp(i, k) = tp(i, k)*(1.+qnk(i)/eps-qnk(i)) !whole thing (convect3)
    END DO
  END DO

  ! ***  Find lifted parcel quantities above cloud base    ***

  DO i = 1, len
    tg = ticb(i)
    ! ori         qg=qs(i,icb(i))
    qg = qsicb(i) ! convect3
    ! debug         alv=lv0-clmcpv*(ticb(i)-t0)
    alv = lv0 - clmcpv*(ticb(i)-273.15)

    ! First iteration.

    ! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
    s = cpd*(1.-qnk(i)) + cl*qnk(i) & ! convect3
      +alv*alv*qg/(rrv*ticb(i)*ticb(i)) ! convect3
    s = 1./s
    ! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
    ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gzicb(i) ! convect3
    tg = tg + s*(ah0(i)-ahg)
    ! ori          tg=max(tg,35.0)
    ! debug          tc=tg-t0
    tc = tg - 273.15
    denom = 243.5 + tc
    denom = max(denom, 1.0) ! convect3
    ! ori          if(tc.ge.0.0)then
    es = 6.112*exp(17.67*tc/denom)
    ! ori          else
    ! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
    ! ori          endif
    ! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
    qg = eps*es/(p(i,icbs(i))-es*(1.-eps))

    ! Second iteration.


    ! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
    ! ori          s=1./s
    ! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
    ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gzicb(i) ! convect3
    tg = tg + s*(ah0(i)-ahg)
    ! ori          tg=max(tg,35.0)
    ! debug          tc=tg-t0
    tc = tg - 273.15
    denom = 243.5 + tc
    denom = max(denom, 1.0) ! convect3
    ! ori          if(tc.ge.0.0)then
    es = 6.112*exp(17.67*tc/denom)
    ! ori          else
    ! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
    ! ori          end if
    ! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
    qg = eps*es/(p(i,icbs(i))-es*(1.-eps))

    alv = lv0 - clmcpv*(ticb(i)-273.15)

    ! ori c approximation here:
    ! ori         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i)
    ! ori     &   -gz(i,icb(i))-alv*qg)/cpd

    ! convect3: no approximation:
    tp(i, icbs(i)) = (ah0(i)-gz(i,icbs(i))-alv*qg)/(cpd+(cl-cpd)*qnk(i))

    ! ori         clw(i,icb(i))=qnk(i)-qg
    ! ori         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
    clw(i, icbs(i)) = qnk(i) - qg
    clw(i, icbs(i)) = max(0.0, clw(i,icbs(i)))

    rg = qg/(1.-qnk(i))
    ! ori         tvp(i,icb(i))=tp(i,icb(i))*(1.+rg*epsi)
    ! convect3: (qg utilise au lieu du vrai mixing ratio rg)
    tvp(i, icbs(i)) = tp(i, icbs(i))*(1.+qg/eps-qnk(i)) !whole thing

  END DO

  ! ori      do 380 k=minorig,icbsmax2
  ! ori       do 370 i=1,len
  ! ori         tvp(i,k)=tvp(i,k)-tp(i,k)*qnk(i)
  ! ori 370   continue
  ! ori 380  continue


  ! -- The following is only for convect3:

  ! * icbs is the first level above the LCL:
  ! if plcl<p(icb), then icbs=icb+1
  ! if plcl>p(icb), then icbs=icb

  ! * the routine above computes tvp from minorig to icbs (included).

  ! * to compute buoybase (in cv3_trigger.F), both tvp(icb) and tvp(icb+1)
  ! must be known. This is the case if icbs=icb+1, but not if icbs=icb.

  ! * therefore, in the case icbs=icb, we compute tvp at level icb+1
  ! (tvp at other levels will be computed in cv3_undilute2.F)


  DO i = 1, len
    ticb(i) = t(i, icb(i)+1)
    gzicb(i) = gz(i, icb(i)+1)
    qsicb(i) = qs(i, icb(i)+1)
  END DO

  DO i = 1, len
    tg = ticb(i)
    qg = qsicb(i) ! convect3
    ! debug         alv=lv0-clmcpv*(ticb(i)-t0)
    alv = lv0 - clmcpv*(ticb(i)-273.15)

    ! First iteration.

    ! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
    s = cpd*(1.-qnk(i)) + cl*qnk(i) & ! convect3
      +alv*alv*qg/(rrv*ticb(i)*ticb(i)) ! convect3
    s = 1./s
    ! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
    ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gzicb(i) ! convect3
    tg = tg + s*(ah0(i)-ahg)
    ! ori          tg=max(tg,35.0)
    ! debug          tc=tg-t0
    tc = tg - 273.15
    denom = 243.5 + tc
    denom = max(denom, 1.0) ! convect3
    ! ori          if(tc.ge.0.0)then
    es = 6.112*exp(17.67*tc/denom)
    ! ori          else
    ! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
    ! ori          endif
    ! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
    qg = eps*es/(p(i,icb(i)+1)-es*(1.-eps))

    ! Second iteration.


    ! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
    ! ori          s=1./s
    ! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
    ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gzicb(i) ! convect3
    tg = tg + s*(ah0(i)-ahg)
    ! ori          tg=max(tg,35.0)
    ! debug          tc=tg-t0
    tc = tg - 273.15
    denom = 243.5 + tc
    denom = max(denom, 1.0) ! convect3
    ! ori          if(tc.ge.0.0)then
    es = 6.112*exp(17.67*tc/denom)
    ! ori          else
    ! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
    ! ori          end if
    ! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
    qg = eps*es/(p(i,icb(i)+1)-es*(1.-eps))

    alv = lv0 - clmcpv*(ticb(i)-273.15)

    ! ori c approximation here:
    ! ori         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i)
    ! ori     &   -gz(i,icb(i))-alv*qg)/cpd

    ! convect3: no approximation:
    tp(i, icb(i)+1) = (ah0(i)-gz(i,icb(i)+1)-alv*qg)/(cpd+(cl-cpd)*qnk(i))

    ! ori         clw(i,icb(i))=qnk(i)-qg
    ! ori         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
    clw(i, icb(i)+1) = qnk(i) - qg
    clw(i, icb(i)+1) = max(0.0, clw(i,icb(i)+1))

    rg = qg/(1.-qnk(i))
    ! ori         tvp(i,icb(i))=tp(i,icb(i))*(1.+rg*epsi)
    ! convect3: (qg utilise au lieu du vrai mixing ratio rg)
    tvp(i, icb(i)+1) = tp(i, icb(i)+1)*(1.+qg/eps-qnk(i)) !whole thing

  END DO

  RETURN
END SUBROUTINE cv30_undilute1

SUBROUTINE cv30_trigger(len, nd, icb, plcl, p, th, tv, tvp, pbase, buoybase, &
    iflag, sig, w0)
  IMPLICIT NONE

  ! -------------------------------------------------------------------
  ! --- TRIGGERING

  ! - computes the cloud base
  ! - triggering (crude in this version)
  ! - relaxation of sig and w0 when no convection

  ! Caution1: if no convection, we set iflag=4
  ! (it used to be 0 in convect3)

  ! Caution2: at this stage, tvp (and thus buoy) are know up
  ! through icb only!
  ! -> the buoyancy below cloud base not (yet) set to the cloud base buoyancy
  ! -------------------------------------------------------------------

  include "cv30param.h"

  ! input:
  INTEGER len, nd
  INTEGER icb(len)
  REAL plcl(len), p(len, nd)
  REAL th(len, nd), tv(len, nd), tvp(len, nd)

  ! output:
  REAL pbase(len), buoybase(len)

  ! input AND output:
  INTEGER iflag(len)
  REAL sig(len, nd), w0(len, nd)

  ! local variables:
  INTEGER i, k
  REAL tvpbase, tvbase, tdif, ath, ath1


  ! ***   set cloud base buoyancy at (plcl+dpbase) level buoyancy

  DO i = 1, len
    pbase(i) = plcl(i) + dpbase
    tvpbase = tvp(i, icb(i))*(pbase(i)-p(i,icb(i)+1))/ &
      (p(i,icb(i))-p(i,icb(i)+1)) + tvp(i, icb(i)+1)*(p(i,icb(i))-pbase(i))/( &
      p(i,icb(i))-p(i,icb(i)+1))
    tvbase = tv(i, icb(i))*(pbase(i)-p(i,icb(i)+1))/ &
      (p(i,icb(i))-p(i,icb(i)+1)) + tv(i, icb(i)+1)*(p(i,icb(i))-pbase(i))/(p &
      (i,icb(i))-p(i,icb(i)+1))
    buoybase(i) = tvpbase - tvbase
  END DO


  ! ***   make sure that column is dry adiabatic between the surface  ***
  ! ***    and cloud base, and that lifted air is positively buoyant  ***
  ! ***                         at cloud base                         ***
  ! ***       if not, return to calling program after resetting       ***
  ! ***                        sig(i) and w0(i)                       ***


  ! oct3      do 200 i=1,len
  ! oct3
  ! oct3       tdif = buoybase(i)
  ! oct3       ath1 = th(i,1)
  ! oct3       ath  = th(i,icb(i)-1) - dttrig
  ! oct3
  ! oct3       if (tdif.lt.dtcrit .or. ath.gt.ath1) then
  ! oct3         do 60 k=1,nl
  ! oct3            sig(i,k) = beta*sig(i,k) - 2.*alpha*tdif*tdif
  ! oct3            sig(i,k) = AMAX1(sig(i,k),0.0)
  ! oct3            w0(i,k)  = beta*w0(i,k)
  ! oct3   60    continue
  ! oct3         iflag(i)=4 ! pour version vectorisee
  ! oct3c convect3         iflag(i)=0
  ! oct3cccc         return
  ! oct3       endif
  ! oct3
  ! oct3200   continue

  ! -- oct3: on reecrit la boucle 200 (pour la vectorisation)

  DO k = 1, nl
    DO i = 1, len

      tdif = buoybase(i)
      ath1 = th(i, 1)
      ath = th(i, icb(i)-1) - dttrig

      IF (tdif<dtcrit .OR. ath>ath1) THEN
        sig(i, k) = beta*sig(i, k) - 2.*alpha*tdif*tdif
        sig(i, k) = amax1(sig(i,k), 0.0)
        w0(i, k) = beta*w0(i, k)
        iflag(i) = 4 ! pour version vectorisee
        ! convect3         iflag(i)=0
      END IF

    END DO
  END DO

  ! fin oct3 --

  RETURN
END SUBROUTINE cv30_trigger

SUBROUTINE cv30_compress(len, nloc, ncum, nd, ntra, iflag1, nk1, icb1, icbs1, &
    plcl1, tnk1, qnk1, gznk1, pbase1, buoybase1, t1, q1, qs1, u1, v1, gz1, &
    th1, tra1, h1, lv1, cpn1, p1, ph1, tv1, tp1, tvp1, clw1, sig1, w01, &
    iflag, nk, icb, icbs, plcl, tnk, qnk, gznk, pbase, buoybase, t, q, qs, u, &
    v, gz, th, tra, h, lv, cpn, p, ph, tv, tp, tvp, clw, sig, w0)
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE

  include "cv30param.h"

  ! inputs:
  INTEGER len, ncum, nd, ntra, nloc
  INTEGER iflag1(len), nk1(len), icb1(len), icbs1(len)
  REAL plcl1(len), tnk1(len), qnk1(len), gznk1(len)
  REAL pbase1(len), buoybase1(len)
  REAL t1(len, nd), q1(len, nd), qs1(len, nd), u1(len, nd), v1(len, nd)
  REAL gz1(len, nd), h1(len, nd), lv1(len, nd), cpn1(len, nd)
  REAL p1(len, nd), ph1(len, nd+1), tv1(len, nd), tp1(len, nd)
  REAL tvp1(len, nd), clw1(len, nd)
  REAL th1(len, nd)
  REAL sig1(len, nd), w01(len, nd)
  REAL tra1(len, nd, ntra)

  ! outputs:
  ! en fait, on a nloc=len pour l'instant (cf cv_driver)
  INTEGER iflag(nloc), nk(nloc), icb(nloc), icbs(nloc)
  REAL plcl(nloc), tnk(nloc), qnk(nloc), gznk(nloc)
  REAL pbase(nloc), buoybase(nloc)
  REAL t(nloc, nd), q(nloc, nd), qs(nloc, nd), u(nloc, nd), v(nloc, nd)
  REAL gz(nloc, nd), h(nloc, nd), lv(nloc, nd), cpn(nloc, nd)
  REAL p(nloc, nd), ph(nloc, nd+1), tv(nloc, nd), tp(nloc, nd)
  REAL tvp(nloc, nd), clw(nloc, nd)
  REAL th(nloc, nd)
  REAL sig(nloc, nd), w0(nloc, nd)
  REAL tra(nloc, nd, ntra)

  ! local variables:
  INTEGER i, k, nn, j

  CHARACTER (LEN=20) :: modname = 'cv30_compress'
  CHARACTER (LEN=80) :: abort_message


  DO k = 1, nl + 1
    nn = 0
    DO i = 1, len
      IF (iflag1(i)==0) THEN
        nn = nn + 1
        sig(nn, k) = sig1(i, k)
        w0(nn, k) = w01(i, k)
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
        th(nn, k) = th1(i, k)
      END IF
    END DO
  END DO

  ! do 121 j=1,ntra
  ! do 111 k=1,nd
  ! nn=0
  ! do 101 i=1,len
  ! if(iflag1(i).eq.0)then
  ! nn=nn+1
  ! tra(nn,k,j)=tra1(i,k,j)
  ! endif
  ! 101  continue
  ! 111  continue
  ! 121  continue

  IF (nn/=ncum) THEN
    WRITE (lunout, *) 'strange! nn not equal to ncum: ', nn, ncum
    abort_message = ''
    CALL abort_physic(modname, abort_message, 1)
  END IF

  nn = 0
  DO i = 1, len
    IF (iflag1(i)==0) THEN
      nn = nn + 1
      pbase(nn) = pbase1(i)
      buoybase(nn) = buoybase1(i)
      plcl(nn) = plcl1(i)
      tnk(nn) = tnk1(i)
      qnk(nn) = qnk1(i)
      gznk(nn) = gznk1(i)
      nk(nn) = nk1(i)
      icb(nn) = icb1(i)
      icbs(nn) = icbs1(i)
      iflag(nn) = iflag1(i)
    END IF
  END DO

  RETURN
END SUBROUTINE cv30_compress

SUBROUTINE cv30_undilute2(nloc, ncum, nd, icb, icbs, nk, tnk, qnk, gznk, t, &
    q, qs, gz, p, h, tv, lv, pbase, buoybase, plcl, inb, tp, tvp, clw, hp, &
    ep, sigp, buoy)
    ! epmax_cape: ajout arguments
  IMPLICIT NONE

  ! ---------------------------------------------------------------------
  ! Purpose:
  ! FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
  ! &
  ! COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
  ! FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
  ! &
  ! FIND THE LEVEL OF NEUTRAL BUOYANCY

  ! Main differences convect3/convect4:
  ! - icbs (input) is the first level above LCL (may differ from icb)
  ! - many minor differences in the iterations
  ! - condensed water not removed from tvp in convect3
  ! - vertical profile of buoyancy computed here (use of buoybase)
  ! - the determination of inb is different
  ! - no inb1, only inb in output
  ! ---------------------------------------------------------------------

  include "cvthermo.h"
  include "cv30param.h"
  include "conema3.h"

  ! inputs:
  INTEGER ncum, nd, nloc
  INTEGER icb(nloc), icbs(nloc), nk(nloc)
  REAL t(nloc, nd), q(nloc, nd), qs(nloc, nd), gz(nloc, nd)
  REAL p(nloc, nd)
  REAL tnk(nloc), qnk(nloc), gznk(nloc)
  REAL lv(nloc, nd), tv(nloc, nd), h(nloc, nd)
  REAL pbase(nloc), buoybase(nloc), plcl(nloc)

  ! outputs:
  INTEGER inb(nloc)
  REAL tp(nloc, nd), tvp(nloc, nd), clw(nloc, nd)
  REAL ep(nloc, nd), sigp(nloc, nd), hp(nloc, nd)
  REAL buoy(nloc, nd)

  ! local variables:
  INTEGER i, k
  REAL tg, qg, ahg, alv, s, tc, es, denom, rg, tca, elacrit
  REAL by, defrac, pden
  REAL ah0(nloc), cape(nloc), capem(nloc), byp(nloc)
  LOGICAL lcape(nloc)

  ! =====================================================================
  ! --- SOME INITIALIZATIONS
  ! =====================================================================

  DO k = 1, nl
    DO i = 1, ncum
      ep(i, k) = 0.0
      sigp(i, k) = spfac
    END DO
  END DO

  ! =====================================================================
  ! --- FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
  ! =====================================================================

  ! ---       The procedure is to solve the equation.
  ! cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.

  ! ***  Calculate certain parcel quantities, including static energy   ***


  DO i = 1, ncum
    ah0(i) = (cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i) & ! debug     &
                                                  ! +qnk(i)*(lv0-clmcpv*(tnk(i)-t0))+gznk(i)
      +qnk(i)*(lv0-clmcpv*(tnk(i)-273.15)) + gznk(i)
  END DO


  ! ***  Find lifted parcel quantities above cloud base    ***


  DO k = minorig + 1, nl
    DO i = 1, ncum
      ! ori	    if(k.ge.(icb(i)+1))then
      IF (k>=(icbs(i)+1)) THEN ! convect3
        tg = t(i, k)
        qg = qs(i, k)
        ! debug	      alv=lv0-clmcpv*(t(i,k)-t0)
        alv = lv0 - clmcpv*(t(i,k)-273.15)

        ! First iteration.

        ! ori	       s=cpd+alv*alv*qg/(rrv*t(i,k)*t(i,k))
        s = cpd*(1.-qnk(i)) + cl*qnk(i) & ! convect3
          +alv*alv*qg/(rrv*t(i,k)*t(i,k)) ! convect3
        s = 1./s
        ! ori	       ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
        ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gz(i, k) ! convect3
        tg = tg + s*(ah0(i)-ahg)
        ! ori	       tg=max(tg,35.0)
        ! debug	       tc=tg-t0
        tc = tg - 273.15
        denom = 243.5 + tc
        denom = max(denom, 1.0) ! convect3
        ! ori	       if(tc.ge.0.0)then
        es = 6.112*exp(17.67*tc/denom)
        ! ori	       else
        ! ori			es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
        ! ori	       endif
        qg = eps*es/(p(i,k)-es*(1.-eps))

        ! Second iteration.

        ! ori	       s=cpd+alv*alv*qg/(rrv*t(i,k)*t(i,k))
        ! ori	       s=1./s
        ! ori	       ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
        ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gz(i, k) ! convect3
        tg = tg + s*(ah0(i)-ahg)
        ! ori	       tg=max(tg,35.0)
        ! debug	       tc=tg-t0
        tc = tg - 273.15
        denom = 243.5 + tc
        denom = max(denom, 1.0) ! convect3
        ! ori	       if(tc.ge.0.0)then
        es = 6.112*exp(17.67*tc/denom)
        ! ori	       else
        ! ori			es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
        ! ori	       endif
        qg = eps*es/(p(i,k)-es*(1.-eps))

        ! debug	       alv=lv0-clmcpv*(t(i,k)-t0)
        alv = lv0 - clmcpv*(t(i,k)-273.15)
        ! print*,'cpd dans convect2 ',cpd
        ! print*,'tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd'
        ! print*,tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd

        ! ori c approximation here:
        ! ori
        ! tp(i,k)=(ah0(i)-(cl-cpd)*qnk(i)*t(i,k)-gz(i,k)-alv*qg)/cpd

        ! convect3: no approximation:
        tp(i, k) = (ah0(i)-gz(i,k)-alv*qg)/(cpd+(cl-cpd)*qnk(i))

        clw(i, k) = qnk(i) - qg
        clw(i, k) = max(0.0, clw(i,k))
        rg = qg/(1.-qnk(i))
        ! ori               tvp(i,k)=tp(i,k)*(1.+rg*epsi)
        ! convect3: (qg utilise au lieu du vrai mixing ratio rg):
        tvp(i, k) = tp(i, k)*(1.+qg/eps-qnk(i)) ! whole thing
      END IF
    END DO
  END DO

  ! =====================================================================
  ! --- SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF
  ! --- PRECIPITATION FALLING OUTSIDE OF CLOUD
  ! --- THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)
  ! =====================================================================

  ! ori      do 320 k=minorig+1,nl
  DO k = 1, nl ! convect3
    DO i = 1, ncum
      pden = ptcrit - pbcrit
      ep(i, k) = (plcl(i)-p(i,k)-pbcrit)/pden*epmax
      ep(i, k) = amax1(ep(i,k), 0.0)
      ep(i, k) = amin1(ep(i,k), epmax)
      sigp(i, k) = spfac
      ! ori          if(k.ge.(nk(i)+1))then
      ! ori            tca=tp(i,k)-t0
      ! ori            if(tca.ge.0.0)then
      ! ori              elacrit=elcrit
      ! ori            else
      ! ori              elacrit=elcrit*(1.0-tca/tlcrit)
      ! ori            endif
      ! ori            elacrit=max(elacrit,0.0)
      ! ori            ep(i,k)=1.0-elacrit/max(clw(i,k),1.0e-8)
      ! ori            ep(i,k)=max(ep(i,k),0.0 )
      ! ori            ep(i,k)=min(ep(i,k),1.0 )
      ! ori            sigp(i,k)=sigs
      ! ori          endif
    END DO
  END DO

  ! =====================================================================
  ! --- CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL
  ! --- VIRTUAL TEMPERATURE
  ! =====================================================================

  ! dans convect3, tvp est calcule en une seule fois, et sans retirer
  ! l'eau condensee (~> reversible CAPE)

  ! ori      do 340 k=minorig+1,nl
  ! ori        do 330 i=1,ncum
  ! ori        if(k.ge.(icb(i)+1))then
  ! ori          tvp(i,k)=tvp(i,k)*(1.0-qnk(i)+ep(i,k)*clw(i,k))
  ! oric         print*,'i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)'
  ! oric         print*, i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)
  ! ori        endif
  ! ori 330    continue
  ! ori 340  continue

  ! ori      do 350 i=1,ncum
  ! ori       tvp(i,nlp)=tvp(i,nl)-(gz(i,nlp)-gz(i,nl))/cpd
  ! ori 350  continue

  DO i = 1, ncum ! convect3
    tp(i, nlp) = tp(i, nl) ! convect3
  END DO ! convect3

  ! =====================================================================
  ! --- EFFECTIVE VERTICAL PROFILE OF BUOYANCY (convect3 only):
  ! =====================================================================

  ! -- this is for convect3 only:

  ! first estimate of buoyancy:

  DO i = 1, ncum
    DO k = 1, nl
      buoy(i, k) = tvp(i, k) - tv(i, k)
    END DO
  END DO

  ! set buoyancy=buoybase for all levels below base
  ! for safety, set buoy(icb)=buoybase

  DO i = 1, ncum
    DO k = 1, nl
      IF ((k>=icb(i)) .AND. (k<=nl) .AND. (p(i,k)>=pbase(i))) THEN
        buoy(i, k) = buoybase(i)
      END IF
    END DO
    ! IM cf. CRio/JYG 270807   buoy(icb(i),k)=buoybase(i)
    buoy(i, icb(i)) = buoybase(i)
  END DO

  ! -- end convect3

  ! =====================================================================
  ! --- FIND THE FIRST MODEL LEVEL (INB) ABOVE THE PARCEL'S
  ! --- LEVEL OF NEUTRAL BUOYANCY
  ! =====================================================================

  ! -- this is for convect3 only:

  DO i = 1, ncum
    inb(i) = nl - 1
  END DO

  DO i = 1, ncum
    DO k = 1, nl - 1
      IF ((k>=icb(i)) .AND. (buoy(i,k)<dtovsh)) THEN
        inb(i) = min(inb(i), k)
      END IF
    END DO
  END DO

  ! -- end convect3

  ! ori      do 510 i=1,ncum
  ! ori        cape(i)=0.0
  ! ori        capem(i)=0.0
  ! ori        inb(i)=icb(i)+1
  ! ori        inb1(i)=inb(i)
  ! ori 510  continue

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

  ! ori      call zilch(byp,ncum)
  ! ori      do 515 i=1,ncum
  ! ori        lcape(i)=.true.
  ! ori 515  continue
  ! ori      do 530 k=minorig+1,nl-1
  ! ori        do 520 i=1,ncum
  ! ori          if(cape(i).lt.0.0)lcape(i)=.false.
  ! ori          if((k.ge.(icb(i)+1)).and.lcape(i))then
  ! ori            by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
  ! ori            byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
  ! ori            cape(i)=cape(i)+by
  ! ori            if(by.ge.0.0)inb1(i)=k+1
  ! ori            if(cape(i).gt.0.0)then
  ! ori              inb(i)=k+1
  ! ori              capem(i)=cape(i)
  ! ori            endif
  ! ori          endif
  ! ori 520    continue
  ! ori 530  continue
  ! ori      do 540 i=1,ncum
  ! ori          cape(i)=capem(i)+byp(i)
  ! ori          defrac=capem(i)-cape(i)
  ! ori          defrac=max(defrac,0.001)
  ! ori          frac(i)=-cape(i)/defrac
  ! ori          frac(i)=min(frac(i),1.0)
  ! ori          frac(i)=max(frac(i),0.0)
  ! ori 540  continue

  ! =====================================================================
  ! ---   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL
  ! =====================================================================

  ! ym      do i=1,ncum*nlp
  ! ym       hp(i,1)=h(i,1)
  ! ym      enddo

  DO k = 1, nlp
    DO i = 1, ncum
      hp(i, k) = h(i, k)
    END DO
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
END SUBROUTINE cv30_undilute2

SUBROUTINE cv30_closure(nloc, ncum, nd, icb, inb, pbase, p, ph, tv, buoy, &
    sig, w0, cape, m)
  IMPLICIT NONE

  ! ===================================================================
  ! ---  CLOSURE OF CONVECT3

  ! vectorization: S. Bony
  ! ===================================================================

  include "cvthermo.h"
  include "cv30param.h"

  ! input:
  INTEGER ncum, nd, nloc
  INTEGER icb(nloc), inb(nloc)
  REAL pbase(nloc)
  REAL p(nloc, nd), ph(nloc, nd+1)
  REAL tv(nloc, nd), buoy(nloc, nd)

  ! input/output:
  REAL sig(nloc, nd), w0(nloc, nd)

  ! output:
  REAL cape(nloc)
  REAL m(nloc, nd)

  ! local variables:
  INTEGER i, j, k, icbmax
  REAL deltap, fac, w, amu
  REAL dtmin(nloc, nd), sigold(nloc, nd)

  ! -------------------------------------------------------
  ! -- Initialization
  ! -------------------------------------------------------

  DO k = 1, nl
    DO i = 1, ncum
      m(i, k) = 0.0
    END DO
  END DO

  ! -------------------------------------------------------
  ! -- Reset sig(i) and w0(i) for i>inb and i<icb
  ! -------------------------------------------------------

  ! update sig and w0 above LNB:

  DO k = 1, nl - 1
    DO i = 1, ncum
      IF ((inb(i)<(nl-1)) .AND. (k>=(inb(i)+1))) THEN
        sig(i, k) = beta*sig(i, k) + 2.*alpha*buoy(i, inb(i))*abs(buoy(i,inb( &
          i)))
        sig(i, k) = amax1(sig(i,k), 0.0)
        w0(i, k) = beta*w0(i, k)
      END IF
    END DO
  END DO

  ! compute icbmax:

  icbmax = 2
  DO i = 1, ncum
    icbmax = max(icbmax, icb(i))
  END DO

  ! update sig and w0 below cloud base:

  DO k = 1, icbmax
    DO i = 1, ncum
      IF (k<=icb(i)) THEN
        sig(i, k) = beta*sig(i, k) - 2.*alpha*buoy(i, icb(i))*buoy(i, icb(i))
        sig(i, k) = amax1(sig(i,k), 0.0)
        w0(i, k) = beta*w0(i, k)
      END IF
    END DO
  END DO

  ! !      if(inb.lt.(nl-1))then
  ! !         do 85 i=inb+1,nl-1
  ! !            sig(i)=beta*sig(i)+2.*alpha*buoy(inb)*
  ! !     1              abs(buoy(inb))
  ! !            sig(i)=amax1(sig(i),0.0)
  ! !            w0(i)=beta*w0(i)
  ! !   85    continue
  ! !      end if

  ! !      do 87 i=1,icb
  ! !         sig(i)=beta*sig(i)-2.*alpha*buoy(icb)*buoy(icb)
  ! !         sig(i)=amax1(sig(i),0.0)
  ! !         w0(i)=beta*w0(i)
  ! !   87 continue

  ! -------------------------------------------------------------
  ! -- Reset fractional areas of updrafts and w0 at initial time
  ! -- and after 10 time steps of no convection
  ! -------------------------------------------------------------

  DO k = 1, nl - 1
    DO i = 1, ncum
      IF (sig(i,nd)<1.5 .OR. sig(i,nd)>12.0) THEN
        sig(i, k) = 0.0
        w0(i, k) = 0.0
      END IF
    END DO
  END DO

  ! -------------------------------------------------------------
  ! -- Calculate convective available potential energy (cape),
  ! -- vertical velocity (w), fractional area covered by
  ! -- undilute updraft (sig), and updraft mass flux (m)
  ! -------------------------------------------------------------

  DO i = 1, ncum
    cape(i) = 0.0
  END DO

  ! compute dtmin (minimum buoyancy between ICB and given level k):

  DO i = 1, ncum
    DO k = 1, nl
      dtmin(i, k) = 100.0
    END DO
  END DO

  DO i = 1, ncum
    DO k = 1, nl
      DO j = minorig, nl
        IF ((k>=(icb(i)+1)) .AND. (k<=inb(i)) .AND. (j>=icb(i)) .AND. (j<=(k- &
            1))) THEN
          dtmin(i, k) = amin1(dtmin(i,k), buoy(i,j))
        END IF
      END DO
    END DO
  END DO

  ! the interval on which cape is computed starts at pbase :
  DO k = 1, nl
    DO i = 1, ncum

      IF ((k>=(icb(i)+1)) .AND. (k<=inb(i))) THEN

        deltap = min(pbase(i), ph(i,k-1)) - min(pbase(i), ph(i,k))
        cape(i) = cape(i) + rrd*buoy(i, k-1)*deltap/p(i, k-1)
        cape(i) = amax1(0.0, cape(i))
        sigold(i, k) = sig(i, k)

        ! dtmin(i,k)=100.0
        ! do 97 j=icb(i),k-1 ! mauvaise vectorisation
        ! dtmin(i,k)=AMIN1(dtmin(i,k),buoy(i,j))
        ! 97     continue

        sig(i, k) = beta*sig(i, k) + alpha*dtmin(i, k)*abs(dtmin(i,k))
        sig(i, k) = amax1(sig(i,k), 0.0)
        sig(i, k) = amin1(sig(i,k), 0.01)
        fac = amin1(((dtcrit-dtmin(i,k))/dtcrit), 1.0)
        w = (1.-beta)*fac*sqrt(cape(i)) + beta*w0(i, k)
        amu = 0.5*(sig(i,k)+sigold(i,k))*w
        m(i, k) = amu*0.007*p(i, k)*(ph(i,k)-ph(i,k+1))/tv(i, k)
        w0(i, k) = w
      END IF

    END DO
  END DO

  DO i = 1, ncum
    w0(i, icb(i)) = 0.5*w0(i, icb(i)+1)
    m(i, icb(i)) = 0.5*m(i, icb(i)+1)*(ph(i,icb(i))-ph(i,icb(i)+1))/ &
      (ph(i,icb(i)+1)-ph(i,icb(i)+2))
    sig(i, icb(i)) = sig(i, icb(i)+1)
    sig(i, icb(i)-1) = sig(i, icb(i))
  END DO


  ! !      cape=0.0
  ! !      do 98 i=icb+1,inb
  ! !         deltap = min(pbase,ph(i-1))-min(pbase,ph(i))
  ! !         cape=cape+rrd*buoy(i-1)*deltap/p(i-1)
  ! !         dcape=rrd*buoy(i-1)*deltap/p(i-1)
  ! !         dlnp=deltap/p(i-1)
  ! !         cape=amax1(0.0,cape)
  ! !         sigold=sig(i)

  ! !         dtmin=100.0
  ! !         do 97 j=icb,i-1
  ! !            dtmin=amin1(dtmin,buoy(j))
  ! !   97    continue

  ! !         sig(i)=beta*sig(i)+alpha*dtmin*abs(dtmin)
  ! !         sig(i)=amax1(sig(i),0.0)
  ! !         sig(i)=amin1(sig(i),0.01)
  ! !         fac=amin1(((dtcrit-dtmin)/dtcrit),1.0)
  ! !         w=(1.-beta)*fac*sqrt(cape)+beta*w0(i)
  ! !         amu=0.5*(sig(i)+sigold)*w
  ! !         m(i)=amu*0.007*p(i)*(ph(i)-ph(i+1))/tv(i)
  ! !         w0(i)=w
  ! !   98 continue
  ! !      w0(icb)=0.5*w0(icb+1)
  ! !      m(icb)=0.5*m(icb+1)*(ph(icb)-ph(icb+1))/(ph(icb+1)-ph(icb+2))
  ! !      sig(icb)=sig(icb+1)
  ! !      sig(icb-1)=sig(icb)

  RETURN
END SUBROUTINE cv30_closure

SUBROUTINE cv30_mixing(nloc, ncum, nd, na, ntra, icb, nk, inb, ph, t, rr, rs, &
    u, v, tra, h, lv, qnk, hp, tv, tvp, ep, clw, m, sig, ment, qent, uent, &
    vent, sij, elij, ments, qents, traent)
  IMPLICIT NONE

  ! ---------------------------------------------------------------------
  ! a faire:
  ! - changer rr(il,1) -> qnk(il)
  ! - vectorisation de la partie normalisation des flux (do 789...)
  ! ---------------------------------------------------------------------

  include "cvthermo.h"
  include "cv30param.h"

  ! inputs:
  INTEGER ncum, nd, na, ntra, nloc
  INTEGER icb(nloc), inb(nloc), nk(nloc)
  REAL sig(nloc, nd)
  REAL qnk(nloc)
  REAL ph(nloc, nd+1)
  REAL t(nloc, nd), rr(nloc, nd), rs(nloc, nd)
  REAL u(nloc, nd), v(nloc, nd)
  REAL tra(nloc, nd, ntra) ! input of convect3
  REAL lv(nloc, na), h(nloc, na), hp(nloc, na)
  REAL tv(nloc, na), tvp(nloc, na), ep(nloc, na), clw(nloc, na)
  REAL m(nloc, na) ! input of convect3

  ! outputs:
  REAL ment(nloc, na, na), qent(nloc, na, na)
  REAL uent(nloc, na, na), vent(nloc, na, na)
  REAL sij(nloc, na, na), elij(nloc, na, na)
  REAL traent(nloc, nd, nd, ntra)
  REAL ments(nloc, nd, nd), qents(nloc, nd, nd)
  REAL sigij(nloc, nd, nd)

  ! local variables:
  INTEGER i, j, k, il, im, jm
  INTEGER num1, num2
  INTEGER nent(nloc, na)
  REAL rti, bf2, anum, denom, dei, altem, cwat, stemp, qp
  REAL alt, smid, sjmin, sjmax, delp, delm
  REAL asij(nloc), smax(nloc), scrit(nloc)
  REAL asum(nloc, nd), bsum(nloc, nd), csum(nloc, nd)
  REAL wgh
  REAL zm(nloc, na)
  LOGICAL lwork(nloc)

  ! =====================================================================
  ! --- INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
  ! =====================================================================

  ! ori        do 360 i=1,ncum*nlp
  DO j = 1, nl
    DO i = 1, ncum
      nent(i, j) = 0
      ! in convect3, m is computed in cv3_closure
      ! ori          m(i,1)=0.0
    END DO
  END DO

  ! ori      do 400 k=1,nlp
  ! ori       do 390 j=1,nlp
  DO j = 1, nl
    DO k = 1, nl
      DO i = 1, ncum
        qent(i, k, j) = rr(i, j)
        uent(i, k, j) = u(i, j)
        vent(i, k, j) = v(i, j)
        elij(i, k, j) = 0.0
        ! ym            ment(i,k,j)=0.0
        ! ym            sij(i,k,j)=0.0
      END DO
    END DO
  END DO

  ! ym
  ment(1:ncum, 1:nd, 1:nd) = 0.0
  sij(1:ncum, 1:nd, 1:nd) = 0.0

  ! do k=1,ntra
  ! do j=1,nd  ! instead nlp
  ! do i=1,nd ! instead nlp
  ! do il=1,ncum
  ! traent(il,i,j,k)=tra(il,j,k)
  ! enddo
  ! enddo
  ! enddo
  ! enddo
  zm(:, :) = 0.

  ! =====================================================================
  ! --- CALCULATE ENTRAINED AIR MASS FLUX (ment), TOTAL WATER MIXING
  ! --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
  ! --- FRACTION (sij)
  ! =====================================================================

  DO i = minorig + 1, nl

    DO j = minorig, nl
      DO il = 1, ncum
        IF ((i>=icb(il)) .AND. (i<=inb(il)) .AND. (j>=(icb(il)- &
            1)) .AND. (j<=inb(il))) THEN

          rti = rr(il, 1) - ep(il, i)*clw(il, i)
          bf2 = 1. + lv(il, j)*lv(il, j)*rs(il, j)/(rrv*t(il,j)*t(il,j)*cpd)
          anum = h(il, j) - hp(il, i) + (cpv-cpd)*t(il, j)*(rti-rr(il,j))
          denom = h(il, i) - hp(il, i) + (cpd-cpv)*(rr(il,i)-rti)*t(il, j)
          dei = denom
          IF (abs(dei)<0.01) dei = 0.01
          sij(il, i, j) = anum/dei
          sij(il, i, i) = 1.0
          altem = sij(il, i, j)*rr(il, i) + (1.-sij(il,i,j))*rti - rs(il, j)
          altem = altem/bf2
          cwat = clw(il, j)*(1.-ep(il,j))
          stemp = sij(il, i, j)
          IF ((stemp<0.0 .OR. stemp>1.0 .OR. altem>cwat) .AND. j>i) THEN
            anum = anum - lv(il, j)*(rti-rs(il,j)-cwat*bf2)
            denom = denom + lv(il, j)*(rr(il,i)-rti)
            IF (abs(denom)<0.01) denom = 0.01
            sij(il, i, j) = anum/denom
            altem = sij(il, i, j)*rr(il, i) + (1.-sij(il,i,j))*rti - &
              rs(il, j)
            altem = altem - (bf2-1.)*cwat
          END IF
          IF (sij(il,i,j)>0.0 .AND. sij(il,i,j)<0.95) THEN
            qent(il, i, j) = sij(il, i, j)*rr(il, i) + (1.-sij(il,i,j))*rti
            uent(il, i, j) = sij(il, i, j)*u(il, i) + &
              (1.-sij(il,i,j))*u(il, nk(il))
            vent(il, i, j) = sij(il, i, j)*v(il, i) + &
              (1.-sij(il,i,j))*v(il, nk(il))
            ! !!!      do k=1,ntra
            ! !!!      traent(il,i,j,k)=sij(il,i,j)*tra(il,i,k)
            ! !!!     :      +(1.-sij(il,i,j))*tra(il,nk(il),k)
            ! !!!      end do
            elij(il, i, j) = altem
            elij(il, i, j) = amax1(0.0, elij(il,i,j))
            ment(il, i, j) = m(il, i)/(1.-sij(il,i,j))
            nent(il, i) = nent(il, i) + 1
          END IF
          sij(il, i, j) = amax1(0.0, sij(il,i,j))
          sij(il, i, j) = amin1(1.0, sij(il,i,j))
        END IF ! new
      END DO
    END DO

    ! do k=1,ntra
    ! do j=minorig,nl
    ! do il=1,ncum
    ! if( (i.ge.icb(il)).and.(i.le.inb(il)).and.
    ! :       (j.ge.(icb(il)-1)).and.(j.le.inb(il)))then
    ! traent(il,i,j,k)=sij(il,i,j)*tra(il,i,k)
    ! :            +(1.-sij(il,i,j))*tra(il,nk(il),k)
    ! endif
    ! enddo
    ! enddo
    ! enddo


    ! ***   if no air can entrain at level i assume that updraft detrains
    ! ***
    ! ***   at that level and calculate detrained air flux and properties
    ! ***


    ! @      do 170 i=icb(il),inb(il)

    DO il = 1, ncum
      IF ((i>=icb(il)) .AND. (i<=inb(il)) .AND. (nent(il,i)==0)) THEN
        ! @      if(nent(il,i).eq.0)then
        ment(il, i, i) = m(il, i)
        qent(il, i, i) = rr(il, nk(il)) - ep(il, i)*clw(il, i)
        uent(il, i, i) = u(il, nk(il))
        vent(il, i, i) = v(il, nk(il))
        elij(il, i, i) = clw(il, i)
        ! MAF      sij(il,i,i)=1.0
        sij(il, i, i) = 0.0
      END IF
    END DO
  END DO

  ! do j=1,ntra
  ! do i=minorig+1,nl
  ! do il=1,ncum
  ! if (i.ge.icb(il) .and. i.le.inb(il) .and. nent(il,i).eq.0) then
  ! traent(il,i,i,j)=tra(il,nk(il),j)
  ! endif
  ! enddo
  ! enddo
  ! enddo

  DO j = minorig, nl
    DO i = minorig, nl
      DO il = 1, ncum
        IF ((j>=(icb(il)-1)) .AND. (j<=inb(il)) .AND. (i>=icb(il)) .AND. (i<= &
            inb(il))) THEN
          sigij(il, i, j) = sij(il, i, j)
        END IF
      END DO
    END DO
  END DO
  ! @      enddo

  ! @170   continue

  ! =====================================================================
  ! ---  NORMALIZE ENTRAINED AIR MASS FLUXES
  ! ---  TO REPRESENT EQUAL PROBABILITIES OF MIXING
  ! =====================================================================

  ! ym      call zilch(asum,ncum*nd)
  ! ym      call zilch(bsum,ncum*nd)
  ! ym      call zilch(csum,ncum*nd)
  CALL zilch(asum, nloc*nd)
  CALL zilch(csum, nloc*nd)
  CALL zilch(csum, nloc*nd)

  DO il = 1, ncum
    lwork(il) = .FALSE.
  END DO

  DO i = minorig + 1, nl

    num1 = 0
    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il)) num1 = num1 + 1
    END DO
    IF (num1<=0) GO TO 789


    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il)) THEN
        lwork(il) = (nent(il,i)/=0)
        qp = rr(il, 1) - ep(il, i)*clw(il, i)
        anum = h(il, i) - hp(il, i) - lv(il, i)*(qp-rs(il,i)) + &
          (cpv-cpd)*t(il, i)*(qp-rr(il,i))
        denom = h(il, i) - hp(il, i) + lv(il, i)*(rr(il,i)-qp) + &
          (cpd-cpv)*t(il, i)*(rr(il,i)-qp)
        IF (abs(denom)<0.01) denom = 0.01
        scrit(il) = anum/denom
        alt = qp - rs(il, i) + scrit(il)*(rr(il,i)-qp)
        IF (scrit(il)<=0.0 .OR. alt<=0.0) scrit(il) = 1.0
        smax(il) = 0.0
        asij(il) = 0.0
      END IF
    END DO

    DO j = nl, minorig, -1

      num2 = 0
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. j>=(icb( &
          il)-1) .AND. j<=inb(il) .AND. lwork(il)) num2 = num2 + 1
      END DO
      IF (num2<=0) GO TO 175

      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. j>=(icb( &
            il)-1) .AND. j<=inb(il) .AND. lwork(il)) THEN

          IF (sij(il,i,j)>1.0E-16 .AND. sij(il,i,j)<0.95) THEN
            wgh = 1.0
            IF (j>i) THEN
              sjmax = amax1(sij(il,i,j+1), smax(il))
              sjmax = amin1(sjmax, scrit(il))
              smax(il) = amax1(sij(il,i,j), smax(il))
              sjmin = amax1(sij(il,i,j-1), smax(il))
              sjmin = amin1(sjmin, scrit(il))
              IF (sij(il,i,j)<(smax(il)-1.0E-16)) wgh = 0.0
              smid = amin1(sij(il,i,j), scrit(il))
            ELSE
              sjmax = amax1(sij(il,i,j+1), scrit(il))
              smid = amax1(sij(il,i,j), scrit(il))
              sjmin = 0.0
              IF (j>1) sjmin = sij(il, i, j-1)
              sjmin = amax1(sjmin, scrit(il))
            END IF
            delp = abs(sjmax-smid)
            delm = abs(sjmin-smid)
            asij(il) = asij(il) + wgh*(delp+delm)
            ment(il, i, j) = ment(il, i, j)*(delp+delm)*wgh
          END IF
        END IF
      END DO

175 END DO

    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il)) THEN
        asij(il) = amax1(1.0E-16, asij(il))
        asij(il) = 1.0/asij(il)
        asum(il, i) = 0.0
        bsum(il, i) = 0.0
        csum(il, i) = 0.0
      END IF
    END DO

    DO j = minorig, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. j>=(icb( &
            il)-1) .AND. j<=inb(il)) THEN
          ment(il, i, j) = ment(il, i, j)*asij(il)
        END IF
      END DO
    END DO

    DO j = minorig, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. j>=(icb( &
            il)-1) .AND. j<=inb(il)) THEN
          asum(il, i) = asum(il, i) + ment(il, i, j)
          ment(il, i, j) = ment(il, i, j)*sig(il, j)
          bsum(il, i) = bsum(il, i) + ment(il, i, j)
        END IF
      END DO
    END DO

    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il)) THEN
        bsum(il, i) = amax1(bsum(il,i), 1.0E-16)
        bsum(il, i) = 1.0/bsum(il, i)
      END IF
    END DO

    DO j = minorig, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. j>=(icb( &
            il)-1) .AND. j<=inb(il)) THEN
          ment(il, i, j) = ment(il, i, j)*asum(il, i)*bsum(il, i)
        END IF
      END DO
    END DO

    DO j = minorig, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. j>=(icb( &
            il)-1) .AND. j<=inb(il)) THEN
          csum(il, i) = csum(il, i) + ment(il, i, j)
        END IF
      END DO
    END DO

    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. &
          csum(il,i)<m(il,i)) THEN
        nent(il, i) = 0
        ment(il, i, i) = m(il, i)
        qent(il, i, i) = rr(il, 1) - ep(il, i)*clw(il, i)
        uent(il, i, i) = u(il, nk(il))
        vent(il, i, i) = v(il, nk(il))
        elij(il, i, i) = clw(il, i)
        ! MAF        sij(il,i,i)=1.0
        sij(il, i, i) = 0.0
      END IF
    END DO ! il

    ! do j=1,ntra
    ! do il=1,ncum
    ! if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
    ! :     .and. csum(il,i).lt.m(il,i) ) then
    ! traent(il,i,i,j)=tra(il,nk(il),j)
    ! endif
    ! enddo
    ! enddo
789 END DO

  ! MAF: renormalisation de MENT
  DO jm = 1, nd
    DO im = 1, nd
      DO il = 1, ncum
        zm(il, im) = zm(il, im) + (1.-sij(il,im,jm))*ment(il, im, jm)
      END DO
    END DO
  END DO

  DO jm = 1, nd
    DO im = 1, nd
      DO il = 1, ncum
        IF (zm(il,im)/=0.) THEN
          ment(il, im, jm) = ment(il, im, jm)*m(il, im)/zm(il, im)
        END IF
      END DO
    END DO
  END DO

  DO jm = 1, nd
    DO im = 1, nd
      DO il = 1, ncum
        qents(il, im, jm) = qent(il, im, jm)
        ments(il, im, jm) = ment(il, im, jm)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cv30_mixing


SUBROUTINE cv30_unsat(nloc, ncum, nd, na, ntra, icb, inb, t, rr, rs, gz, u, &
    v, tra, p, ph, th, tv, lv, cpn, ep, sigp, clw, m, ment, elij, delt, plcl, &
    mp, rp, up, vp, trap, wt, water, evap, b & ! RomP-jyg
    , wdtraina, wdtrainm) ! 26/08/10  RomP-jyg
  IMPLICIT NONE


  include "cvthermo.h"
  include "cv30param.h"
  include "cvflag.h"

  ! inputs:
  INTEGER ncum, nd, na, ntra, nloc
  INTEGER icb(nloc), inb(nloc)
  REAL delt, plcl(nloc)
  REAL t(nloc, nd), rr(nloc, nd), rs(nloc, nd)
  REAL u(nloc, nd), v(nloc, nd)
  REAL tra(nloc, nd, ntra)
  REAL p(nloc, nd), ph(nloc, nd+1)
  REAL th(nloc, na), gz(nloc, na)
  REAL lv(nloc, na), ep(nloc, na), sigp(nloc, na), clw(nloc, na)
  REAL cpn(nloc, na), tv(nloc, na)
  REAL m(nloc, na), ment(nloc, na, na), elij(nloc, na, na)

  ! outputs:
  REAL mp(nloc, na), rp(nloc, na), up(nloc, na), vp(nloc, na)
  REAL water(nloc, na), evap(nloc, na), wt(nloc, na)
  REAL trap(nloc, na, ntra)
  REAL b(nloc, na)
  ! 25/08/10 - RomP---- ajout des masses precipitantes ejectees
  ! lascendance adiabatique et des flux melanges Pa et Pm.
  ! Distinction des wdtrain
  ! Pa = wdtrainA     Pm = wdtrainM
  REAL wdtraina(nloc, na), wdtrainm(nloc, na)

  ! local variables
  INTEGER i, j, k, il, num1
  REAL tinv, delti
  REAL awat, afac, afac1, afac2, bfac
  REAL pr1, pr2, sigt, b6, c6, revap, tevap, delth
  REAL amfac, amp2, xf, tf, fac2, ur, sru, fac, d, af, bf
  REAL ampmax
  REAL lvcp(nloc, na)
  REAL wdtrain(nloc)
  LOGICAL lwork(nloc)


  ! ------------------------------------------------------

  delti = 1./delt
  tinv = 1./3.

  mp(:, :) = 0.

  DO i = 1, nl
    DO il = 1, ncum
      mp(il, i) = 0.0
      rp(il, i) = rr(il, i)
      up(il, i) = u(il, i)
      vp(il, i) = v(il, i)
      wt(il, i) = 0.001
      water(il, i) = 0.0
      evap(il, i) = 0.0
      b(il, i) = 0.0
      lvcp(il, i) = lv(il, i)/cpn(il, i)
    END DO
  END DO

  ! do k=1,ntra
  ! do i=1,nd
  ! do il=1,ncum
  ! trap(il,i,k)=tra(il,i,k)
  ! enddo
  ! enddo
  ! enddo
  ! ! RomP >>>
  DO i = 1, nd
    DO il = 1, ncum
      wdtraina(il, i) = 0.0
      wdtrainm(il, i) = 0.0
    END DO
  END DO
  ! ! RomP <<<

  ! ***  check whether ep(inb)=0, if so, skip precipitating    ***
  ! ***             downdraft calculation                      ***


  DO il = 1, ncum
    lwork(il) = .TRUE.
    IF (ep(il,inb(il))<0.0001) lwork(il) = .FALSE.
  END DO

  CALL zilch(wdtrain, ncum)

  DO i = nl + 1, 1, -1

    num1 = 0
    DO il = 1, ncum
      IF (i<=inb(il) .AND. lwork(il)) num1 = num1 + 1
    END DO
    IF (num1<=0) GO TO 400


    ! ***  integrate liquid water equation to find condensed water   ***
    ! ***                and condensed water flux                    ***



    ! ***                    begin downdraft loop                    ***



    ! ***              calculate detrained precipitation             ***

    DO il = 1, ncum
      IF (i<=inb(il) .AND. lwork(il)) THEN
        IF (cvflag_grav) THEN
          wdtrain(il) = grav*ep(il, i)*m(il, i)*clw(il, i)
          wdtraina(il, i) = wdtrain(il)/grav !   Pa  26/08/10   RomP
        ELSE
          wdtrain(il) = 10.0*ep(il, i)*m(il, i)*clw(il, i)
          wdtraina(il, i) = wdtrain(il)/10. !   Pa  26/08/10   RomP
        END IF
      END IF
    END DO

    IF (i>1) THEN

      DO j = 1, i - 1
        DO il = 1, ncum
          IF (i<=inb(il) .AND. lwork(il)) THEN
            awat = elij(il, j, i) - (1.-ep(il,i))*clw(il, i)
            awat = amax1(awat, 0.0)
            IF (cvflag_grav) THEN
              wdtrain(il) = wdtrain(il) + grav*awat*ment(il, j, i)
            ELSE
              wdtrain(il) = wdtrain(il) + 10.0*awat*ment(il, j, i)
            END IF
          END IF
        END DO
      END DO
      DO il = 1, ncum
        IF (cvflag_grav) THEN
          wdtrainm(il, i) = wdtrain(il)/grav - wdtraina(il, i) !   Pm  26/08/10   RomP
        ELSE
          wdtrainm(il, i) = wdtrain(il)/10. - wdtraina(il, i) !   Pm  26/08/10   RomP
        END IF
      END DO

    END IF


    ! ***    find rain water and evaporation using provisional   ***
    ! ***              estimates of rp(i)and rp(i-1)             ***


    DO il = 1, ncum

      IF (i<=inb(il) .AND. lwork(il)) THEN

        wt(il, i) = 45.0

        IF (i<inb(il)) THEN
          rp(il, i) = rp(il, i+1) + (cpd*(t(il,i+1)-t(il, &
            i))+gz(il,i+1)-gz(il,i))/lv(il, i)
          rp(il, i) = 0.5*(rp(il,i)+rr(il,i))
        END IF
        rp(il, i) = amax1(rp(il,i), 0.0)
        rp(il, i) = amin1(rp(il,i), rs(il,i))
        rp(il, inb(il)) = rr(il, inb(il))

        IF (i==1) THEN
          afac = p(il, 1)*(rs(il,1)-rp(il,1))/(1.0E4+2000.0*p(il,1)*rs(il,1))
        ELSE
          rp(il, i-1) = rp(il, i) + (cpd*(t(il,i)-t(il, &
            i-1))+gz(il,i)-gz(il,i-1))/lv(il, i)
          rp(il, i-1) = 0.5*(rp(il,i-1)+rr(il,i-1))
          rp(il, i-1) = amin1(rp(il,i-1), rs(il,i-1))
          rp(il, i-1) = amax1(rp(il,i-1), 0.0)
          afac1 = p(il, i)*(rs(il,i)-rp(il,i))/(1.0E4+2000.0*p(il,i)*rs(il,i) &
            )
          afac2 = p(il, i-1)*(rs(il,i-1)-rp(il,i-1))/ &
            (1.0E4+2000.0*p(il,i-1)*rs(il,i-1))
          afac = 0.5*(afac1+afac2)
        END IF
        IF (i==inb(il)) afac = 0.0
        afac = amax1(afac, 0.0)
        bfac = 1./(sigd*wt(il,i))

        ! jyg1
        ! cc        sigt=1.0
        ! cc        if(i.ge.icb)sigt=sigp(i)
        ! prise en compte de la variation progressive de sigt dans
        ! les couches icb et icb-1:
        ! pour plcl<ph(i+1), pr1=0 & pr2=1
        ! pour plcl>ph(i),   pr1=1 & pr2=0
        ! pour ph(i+1)<plcl<ph(i), pr1 est la proportion a cheval
        ! sur le nuage, et pr2 est la proportion sous la base du
        ! nuage.
        pr1 = (plcl(il)-ph(il,i+1))/(ph(il,i)-ph(il,i+1))
        pr1 = max(0., min(1.,pr1))
        pr2 = (ph(il,i)-plcl(il))/(ph(il,i)-ph(il,i+1))
        pr2 = max(0., min(1.,pr2))
        sigt = sigp(il, i)*pr1 + pr2
        ! jyg2

        b6 = bfac*50.*sigd*(ph(il,i)-ph(il,i+1))*sigt*afac
        c6 = water(il, i+1) + bfac*wdtrain(il) - 50.*sigd*bfac*(ph(il,i)-ph( &
          il,i+1))*evap(il, i+1)
        IF (c6>0.0) THEN
          revap = 0.5*(-b6+sqrt(b6*b6+4.*c6))
          evap(il, i) = sigt*afac*revap
          water(il, i) = revap*revap
        ELSE
          evap(il, i) = -evap(il, i+1) + 0.02*(wdtrain(il)+sigd*wt(il,i)* &
            water(il,i+1))/(sigd*(ph(il,i)-ph(il,i+1)))
        END IF

        ! ***  calculate precipitating downdraft mass flux under     ***
        ! ***              hydrostatic approximation                 ***

        IF (i/=1) THEN

          tevap = amax1(0.0, evap(il,i))
          delth = amax1(0.001, (th(il,i)-th(il,i-1)))
          IF (cvflag_grav) THEN
            mp(il, i) = 100.*ginv*lvcp(il, i)*sigd*tevap*(p(il,i-1)-p(il,i))/ &
              delth
          ELSE
            mp(il, i) = 10.*lvcp(il, i)*sigd*tevap*(p(il,i-1)-p(il,i))/delth
          END IF

          ! ***           if hydrostatic assumption fails,             ***
          ! ***   solve cubic difference equation for downdraft theta  ***
          ! ***  and mass flux from two simultaneous differential eqns ***

          amfac = sigd*sigd*70.0*ph(il, i)*(p(il,i-1)-p(il,i))* &
            (th(il,i)-th(il,i-1))/(tv(il,i)*th(il,i))
          amp2 = abs(mp(il,i+1)*mp(il,i+1)-mp(il,i)*mp(il,i))
          IF (amp2>(0.1*amfac)) THEN
            xf = 100.0*sigd*sigd*sigd*(ph(il,i)-ph(il,i+1))
            tf = b(il, i) - 5.0*(th(il,i)-th(il,i-1))*t(il, i)/(lvcp(il,i)* &
              sigd*th(il,i))
            af = xf*tf + mp(il, i+1)*mp(il, i+1)*tinv
            bf = 2.*(tinv*mp(il,i+1))**3 + tinv*mp(il, i+1)*xf*tf + &
              50.*(p(il,i-1)-p(il,i))*xf*tevap
            fac2 = 1.0
            IF (bf<0.0) fac2 = -1.0
            bf = abs(bf)
            ur = 0.25*bf*bf - af*af*af*tinv*tinv*tinv
            IF (ur>=0.0) THEN
              sru = sqrt(ur)
              fac = 1.0
              IF ((0.5*bf-sru)<0.0) fac = -1.0
              mp(il, i) = mp(il, i+1)*tinv + (0.5*bf+sru)**tinv + &
                fac*(abs(0.5*bf-sru))**tinv
            ELSE
              d = atan(2.*sqrt(-ur)/(bf+1.0E-28))
              IF (fac2<0.0) d = 3.14159 - d
              mp(il, i) = mp(il, i+1)*tinv + 2.*sqrt(af*tinv)*cos(d*tinv)
            END IF
            mp(il, i) = amax1(0.0, mp(il,i))

            IF (cvflag_grav) THEN
              ! jyg : il y a vraisemblablement une erreur dans la ligne 2
              ! suivante:
              ! il faut diviser par (mp(il,i)*sigd*grav) et non par
              ! (mp(il,i)+sigd*0.1).
              ! Et il faut bien revoir les facteurs 100.
              b(il, i-1) = b(il, i) + 100.0*(p(il,i-1)-p(il,i))*tevap/(mp(il, &
                i)+sigd*0.1) - 10.0*(th(il,i)-th(il,i-1))*t(il, i)/(lvcp(il,i &
                )*sigd*th(il,i))
            ELSE
              b(il, i-1) = b(il, i) + 100.0*(p(il,i-1)-p(il,i))*tevap/(mp(il, &
                i)+sigd*0.1) - 10.0*(th(il,i)-th(il,i-1))*t(il, i)/(lvcp(il,i &
                )*sigd*th(il,i))
            END IF
            b(il, i-1) = amax1(b(il,i-1), 0.0)
          END IF

          ! ***         limit magnitude of mp(i) to meet cfl condition
          ! ***

          ampmax = 2.0*(ph(il,i)-ph(il,i+1))*delti
          amp2 = 2.0*(ph(il,i-1)-ph(il,i))*delti
          ampmax = amin1(ampmax, amp2)
          mp(il, i) = amin1(mp(il,i), ampmax)

          ! ***      force mp to decrease linearly to zero
          ! ***
          ! ***       between cloud base and the surface
          ! ***

          IF (p(il,i)>p(il,icb(il))) THEN
            mp(il, i) = mp(il, icb(il))*(p(il,1)-p(il,i))/ &
              (p(il,1)-p(il,icb(il)))
          END IF

        END IF ! i.eq.1

        ! ***       find mixing ratio of precipitating downdraft     ***


        IF (i/=inb(il)) THEN

          rp(il, i) = rr(il, i)

          IF (mp(il,i)>mp(il,i+1)) THEN

            IF (cvflag_grav) THEN
              rp(il, i) = rp(il, i+1)*mp(il, i+1) + &
                rr(il, i)*(mp(il,i)-mp(il,i+1)) + 100.*ginv*0.5*sigd*(ph(il,i &
                )-ph(il,i+1))*(evap(il,i+1)+evap(il,i))
            ELSE
              rp(il, i) = rp(il, i+1)*mp(il, i+1) + &
                rr(il, i)*(mp(il,i)-mp(il,i+1)) + 5.*sigd*(ph(il,i)-ph(il,i+1 &
                ))*(evap(il,i+1)+evap(il,i))
            END IF
            rp(il, i) = rp(il, i)/mp(il, i)
            up(il, i) = up(il, i+1)*mp(il, i+1) + u(il, i)*(mp(il,i)-mp(il,i+ &
              1))
            up(il, i) = up(il, i)/mp(il, i)
            vp(il, i) = vp(il, i+1)*mp(il, i+1) + v(il, i)*(mp(il,i)-mp(il,i+ &
              1))
            vp(il, i) = vp(il, i)/mp(il, i)

            ! do j=1,ntra
            ! trap(il,i,j)=trap(il,i+1,j)*mp(il,i+1)
            ! testmaf     :            +trap(il,i,j)*(mp(il,i)-mp(il,i+1))
            ! :            +tra(il,i,j)*(mp(il,i)-mp(il,i+1))
            ! trap(il,i,j)=trap(il,i,j)/mp(il,i)
            ! end do

          ELSE

            IF (mp(il,i+1)>1.0E-16) THEN
              IF (cvflag_grav) THEN
                rp(il, i) = rp(il, i+1) + 100.*ginv*0.5*sigd*(ph(il,i)-ph(il, &
                  i+1))*(evap(il,i+1)+evap(il,i))/mp(il, i+1)
              ELSE
                rp(il, i) = rp(il, i+1) + 5.*sigd*(ph(il,i)-ph(il,i+1))*(evap &
                  (il,i+1)+evap(il,i))/mp(il, i+1)
              END IF
              up(il, i) = up(il, i+1)
              vp(il, i) = vp(il, i+1)

              ! do j=1,ntra
              ! trap(il,i,j)=trap(il,i+1,j)
              ! end do

            END IF
          END IF
          rp(il, i) = amin1(rp(il,i), rs(il,i))
          rp(il, i) = amax1(rp(il,i), 0.0)

        END IF
      END IF
    END DO

400 END DO

  RETURN
END SUBROUTINE cv30_unsat

SUBROUTINE cv30_yield(nloc, ncum, nd, na, ntra, icb, inb, delt, t, rr, u, v, &
    tra, gz, p, ph, h, hp, lv, cpn, th, ep, clw, m, tp, mp, rp, up, vp, trap, &
    wt, water, evap, b, ment, qent, uent, vent, nent, elij, traent, sig, tv, &
    tvp, iflag, precip, vprecip, ft, fr, fu, fv, ftra, upwd, dnwd, dnwd0, ma, &
    mike, tls, tps, qcondc, wd)
  IMPLICIT NONE

  include "cvthermo.h"
  include "cv30param.h"
  include "cvflag.h"
  include "conema3.h"

  ! inputs:
  INTEGER ncum, nd, na, ntra, nloc
  INTEGER icb(nloc), inb(nloc)
  REAL delt
  REAL t(nloc, nd), rr(nloc, nd), u(nloc, nd), v(nloc, nd)
  REAL tra(nloc, nd, ntra), sig(nloc, nd)
  REAL gz(nloc, na), ph(nloc, nd+1), h(nloc, na), hp(nloc, na)
  REAL th(nloc, na), p(nloc, nd), tp(nloc, na)
  REAL lv(nloc, na), cpn(nloc, na), ep(nloc, na), clw(nloc, na)
  REAL m(nloc, na), mp(nloc, na), rp(nloc, na), up(nloc, na)
  REAL vp(nloc, na), wt(nloc, nd), trap(nloc, nd, ntra)
  REAL water(nloc, na), evap(nloc, na), b(nloc, na)
  REAL ment(nloc, na, na), qent(nloc, na, na), uent(nloc, na, na)
  ! ym      real vent(nloc,na,na), nent(nloc,na), elij(nloc,na,na)
  REAL vent(nloc, na, na), elij(nloc, na, na)
  INTEGER nent(nloc, na)
  REAL traent(nloc, na, na, ntra)
  REAL tv(nloc, nd), tvp(nloc, nd)

  ! input/output:
  INTEGER iflag(nloc)

  ! outputs:
  REAL precip(nloc)
  REAL vprecip(nloc, nd+1)
  REAL ft(nloc, nd), fr(nloc, nd), fu(nloc, nd), fv(nloc, nd)
  REAL ftra(nloc, nd, ntra)
  REAL upwd(nloc, nd), dnwd(nloc, nd), ma(nloc, nd)
  REAL dnwd0(nloc, nd), mike(nloc, nd)
  REAL tls(nloc, nd), tps(nloc, nd)
  REAL qcondc(nloc, nd) ! cld
  REAL wd(nloc) ! gust

  ! local variables:
  INTEGER i, k, il, n, j, num1
  REAL rat, awat, delti
  REAL ax, bx, cx, dx, ex
  REAL cpinv, rdcp, dpinv
  REAL lvcp(nloc, na), mke(nloc, na)
  REAL am(nloc), work(nloc), ad(nloc), amp1(nloc)
  ! !!      real up1(nloc), dn1(nloc)
  REAL up1(nloc, nd, nd), dn1(nloc, nd, nd)
  REAL asum(nloc), bsum(nloc), csum(nloc), dsum(nloc)
  REAL qcond(nloc, nd), nqcond(nloc, nd), wa(nloc, nd) ! cld
  REAL siga(nloc, nd), sax(nloc, nd), mac(nloc, nd) ! cld


  ! -------------------------------------------------------------

  ! initialization:

  delti = 1.0/delt

  DO il = 1, ncum
    precip(il) = 0.0
    wd(il) = 0.0 ! gust
    vprecip(il, nd+1) = 0.
  END DO

  DO i = 1, nd
    DO il = 1, ncum
      vprecip(il, i) = 0.0
      ft(il, i) = 0.0
      fr(il, i) = 0.0
      fu(il, i) = 0.0
      fv(il, i) = 0.0
      qcondc(il, i) = 0.0 ! cld
      qcond(il, i) = 0.0 ! cld
      nqcond(il, i) = 0.0 ! cld
    END DO
  END DO

  ! do j=1,ntra
  ! do i=1,nd
  ! do il=1,ncum
  ! ftra(il,i,j)=0.0
  ! enddo
  ! enddo
  ! enddo

  DO i = 1, nl
    DO il = 1, ncum
      lvcp(il, i) = lv(il, i)/cpn(il, i)
    END DO
  END DO



  ! ***  calculate surface precipitation in mm/day     ***

  DO il = 1, ncum
    IF (ep(il,inb(il))>=0.0001) THEN
      IF (cvflag_grav) THEN
        precip(il) = wt(il, 1)*sigd*water(il, 1)*86400.*1000./(rowl*grav)
      ELSE
        precip(il) = wt(il, 1)*sigd*water(il, 1)*8640.
      END IF
    END IF
  END DO

  ! ***  CALCULATE VERTICAL PROFILE OF  PRECIPITATIONs IN kg/m2/s  ===

  ! MAF rajout pour lessivage
  DO k = 1, nl
    DO il = 1, ncum
      IF (k<=inb(il)) THEN
        IF (cvflag_grav) THEN
          vprecip(il, k) = wt(il, k)*sigd*water(il, k)/grav
        ELSE
          vprecip(il, k) = wt(il, k)*sigd*water(il, k)/10.
        END IF
      END IF
    END DO
  END DO


  ! ***  Calculate downdraft velocity scale    ***
  ! ***  NE PAS UTILISER POUR L'INSTANT ***

  ! !      do il=1,ncum
  ! !        wd(il)=betad*abs(mp(il,icb(il)))*0.01*rrd*t(il,icb(il))
  ! !     :                                  /(sigd*p(il,icb(il)))
  ! !      enddo


  ! ***  calculate tendencies of lowest level potential temperature  ***
  ! ***                      and mixing ratio                        ***

  DO il = 1, ncum
    work(il) = 1.0/(ph(il,1)-ph(il,2))
    am(il) = 0.0
  END DO

  DO k = 2, nl
    DO il = 1, ncum
      IF (k<=inb(il)) THEN
        am(il) = am(il) + m(il, k)
      END IF
    END DO
  END DO

  DO il = 1, ncum

    ! convect3      if((0.1*dpinv*am).ge.delti)iflag(il)=4
    IF (cvflag_grav) THEN
      IF ((0.01*grav*work(il)*am(il))>=delti) iflag(il) = 1 !consist vect
      ft(il, 1) = 0.01*grav*work(il)*am(il)*(t(il,2)-t(il,1)+(gz(il,2)-gz(il, &
        1))/cpn(il,1))
    ELSE
      IF ((0.1*work(il)*am(il))>=delti) iflag(il) = 1 !consistency vect
      ft(il, 1) = 0.1*work(il)*am(il)*(t(il,2)-t(il,1)+(gz(il,2)-gz(il, &
        1))/cpn(il,1))
    END IF

    ft(il, 1) = ft(il, 1) - 0.5*lvcp(il, 1)*sigd*(evap(il,1)+evap(il,2))

    IF (cvflag_grav) THEN
      ft(il, 1) = ft(il, 1) - 0.009*grav*sigd*mp(il, 2)*t(il, 1)*b(il, 1)* &
        work(il)
    ELSE
      ft(il, 1) = ft(il, 1) - 0.09*sigd*mp(il, 2)*t(il, 1)*b(il, 1)*work(il)
    END IF

    ft(il, 1) = ft(il, 1) + 0.01*sigd*wt(il, 1)*(cl-cpd)*water(il, 2)*(t(il,2 &
      )-t(il,1))*work(il)/cpn(il, 1)

    IF (cvflag_grav) THEN
      ! jyg1  Correction pour mieux conserver l'eau (conformite avec
      ! CONVECT4.3)
      ! (sb: pour l'instant, on ne fait que le chgt concernant grav, pas
      ! evap)
      fr(il, 1) = 0.01*grav*mp(il, 2)*(rp(il,2)-rr(il,1))*work(il) + &
        sigd*0.5*(evap(il,1)+evap(il,2))
      ! +tard     :          +sigd*evap(il,1)

      fr(il, 1) = fr(il, 1) + 0.01*grav*am(il)*(rr(il,2)-rr(il,1))*work(il)

      fu(il, 1) = fu(il, 1) + 0.01*grav*work(il)*(mp(il,2)*(up(il,2)-u(il, &
        1))+am(il)*(u(il,2)-u(il,1)))
      fv(il, 1) = fv(il, 1) + 0.01*grav*work(il)*(mp(il,2)*(vp(il,2)-v(il, &
        1))+am(il)*(v(il,2)-v(il,1)))
    ELSE ! cvflag_grav
      fr(il, 1) = 0.1*mp(il, 2)*(rp(il,2)-rr(il,1))*work(il) + &
        sigd*0.5*(evap(il,1)+evap(il,2))
      fr(il, 1) = fr(il, 1) + 0.1*am(il)*(rr(il,2)-rr(il,1))*work(il)
      fu(il, 1) = fu(il, 1) + 0.1*work(il)*(mp(il,2)*(up(il,2)-u(il, &
        1))+am(il)*(u(il,2)-u(il,1)))
      fv(il, 1) = fv(il, 1) + 0.1*work(il)*(mp(il,2)*(vp(il,2)-v(il, &
        1))+am(il)*(v(il,2)-v(il,1)))
    END IF ! cvflag_grav

  END DO ! il

  ! do j=1,ntra
  ! do il=1,ncum
  ! if (cvflag_grav) then
  ! ftra(il,1,j)=ftra(il,1,j)+0.01*grav*work(il)
  ! :                     *(mp(il,2)*(trap(il,2,j)-tra(il,1,j))
  ! :             +am(il)*(tra(il,2,j)-tra(il,1,j)))
  ! else
  ! ftra(il,1,j)=ftra(il,1,j)+0.1*work(il)
  ! :                     *(mp(il,2)*(trap(il,2,j)-tra(il,1,j))
  ! :             +am(il)*(tra(il,2,j)-tra(il,1,j)))
  ! endif
  ! enddo
  ! enddo

  DO j = 2, nl
    DO il = 1, ncum
      IF (j<=inb(il)) THEN
        IF (cvflag_grav) THEN
          fr(il, 1) = fr(il, 1) + 0.01*grav*work(il)*ment(il, j, 1)*(qent(il, &
            j,1)-rr(il,1))
          fu(il, 1) = fu(il, 1) + 0.01*grav*work(il)*ment(il, j, 1)*(uent(il, &
            j,1)-u(il,1))
          fv(il, 1) = fv(il, 1) + 0.01*grav*work(il)*ment(il, j, 1)*(vent(il, &
            j,1)-v(il,1))
        ELSE ! cvflag_grav
          fr(il, 1) = fr(il, 1) + 0.1*work(il)*ment(il, j, 1)*(qent(il,j,1)- &
            rr(il,1))
          fu(il, 1) = fu(il, 1) + 0.1*work(il)*ment(il, j, 1)*(uent(il,j,1)-u &
            (il,1))
          fv(il, 1) = fv(il, 1) + 0.1*work(il)*ment(il, j, 1)*(vent(il,j,1)-v &
            (il,1))
        END IF ! cvflag_grav
      END IF ! j
    END DO
  END DO

  ! do k=1,ntra
  ! do j=2,nl
  ! do il=1,ncum
  ! if (j.le.inb(il)) then

  ! if (cvflag_grav) then
  ! ftra(il,1,k)=ftra(il,1,k)+0.01*grav*work(il)*ment(il,j,1)
  ! :                *(traent(il,j,1,k)-tra(il,1,k))
  ! else
  ! ftra(il,1,k)=ftra(il,1,k)+0.1*work(il)*ment(il,j,1)
  ! :                *(traent(il,j,1,k)-tra(il,1,k))
  ! endif

  ! endif
  ! enddo
  ! enddo
  ! enddo


  ! ***  calculate tendencies of potential temperature and mixing ratio  ***
  ! ***               at levels above the lowest level                   ***

  ! ***  first find the net saturated updraft and downdraft mass fluxes  ***
  ! ***                      through each level                          ***


  DO i = 2, nl + 1 ! newvecto: mettre nl au lieu nl+1?

    num1 = 0
    DO il = 1, ncum
      IF (i<=inb(il)) num1 = num1 + 1
    END DO
    IF (num1<=0) GO TO 500

    CALL zilch(amp1, ncum)
    CALL zilch(ad, ncum)

    DO k = i + 1, nl + 1
      DO il = 1, ncum
        IF (i<=inb(il) .AND. k<=(inb(il)+1)) THEN
          amp1(il) = amp1(il) + m(il, k)
        END IF
      END DO
    END DO

    DO k = 1, i
      DO j = i + 1, nl + 1
        DO il = 1, ncum
          IF (i<=inb(il) .AND. j<=(inb(il)+1)) THEN
            amp1(il) = amp1(il) + ment(il, k, j)
          END IF
        END DO
      END DO
    END DO

    DO k = 1, i - 1
      DO j = i, nl + 1 ! newvecto: nl au lieu nl+1?
        DO il = 1, ncum
          IF (i<=inb(il) .AND. j<=inb(il)) THEN
            ad(il) = ad(il) + ment(il, j, k)
          END IF
        END DO
      END DO
    END DO

    DO il = 1, ncum
      IF (i<=inb(il)) THEN
        dpinv = 1.0/(ph(il,i)-ph(il,i+1))
        cpinv = 1.0/cpn(il, i)

        ! convect3      if((0.1*dpinv*amp1).ge.delti)iflag(il)=4
        IF (cvflag_grav) THEN
          IF ((0.01*grav*dpinv*amp1(il))>=delti) iflag(il) = 1 ! vecto
        ELSE
          IF ((0.1*dpinv*amp1(il))>=delti) iflag(il) = 1 ! vecto
        END IF

        IF (cvflag_grav) THEN
          ft(il, i) = 0.01*grav*dpinv*(amp1(il)*(t(il,i+1)-t(il, &
            i)+(gz(il,i+1)-gz(il,i))*cpinv)-ad(il)*(t(il,i)-t(il, &
            i-1)+(gz(il,i)-gz(il,i-1))*cpinv)) - 0.5*sigd*lvcp(il, i)*(evap( &
            il,i)+evap(il,i+1))
          rat = cpn(il, i-1)*cpinv
          ft(il, i) = ft(il, i) - 0.009*grav*sigd*(mp(il,i+1)*t(il,i)*b(il,i) &
            -mp(il,i)*t(il,i-1)*rat*b(il,i-1))*dpinv
          ft(il, i) = ft(il, i) + 0.01*grav*dpinv*ment(il, i, i)*(hp(il,i)-h( &
            il,i)+t(il,i)*(cpv-cpd)*(rr(il,i)-qent(il,i,i)))*cpinv
        ELSE ! cvflag_grav
          ft(il, i) = 0.1*dpinv*(amp1(il)*(t(il,i+1)-t(il, &
            i)+(gz(il,i+1)-gz(il,i))*cpinv)-ad(il)*(t(il,i)-t(il, &
            i-1)+(gz(il,i)-gz(il,i-1))*cpinv)) - 0.5*sigd*lvcp(il, i)*(evap( &
            il,i)+evap(il,i+1))
          rat = cpn(il, i-1)*cpinv
          ft(il, i) = ft(il, i) - 0.09*sigd*(mp(il,i+1)*t(il,i)*b(il,i)-mp(il &
            ,i)*t(il,i-1)*rat*b(il,i-1))*dpinv
          ft(il, i) = ft(il, i) + 0.1*dpinv*ment(il, i, i)*(hp(il,i)-h(il,i)+ &
            t(il,i)*(cpv-cpd)*(rr(il,i)-qent(il,i,i)))*cpinv
        END IF ! cvflag_grav


        ft(il, i) = ft(il, i) + 0.01*sigd*wt(il, i)*(cl-cpd)*water(il, i+1)*( &
          t(il,i+1)-t(il,i))*dpinv*cpinv

        IF (cvflag_grav) THEN
          fr(il, i) = 0.01*grav*dpinv*(amp1(il)*(rr(il,i+1)-rr(il, &
            i))-ad(il)*(rr(il,i)-rr(il,i-1)))
          fu(il, i) = fu(il, i) + 0.01*grav*dpinv*(amp1(il)*(u(il,i+1)-u(il, &
            i))-ad(il)*(u(il,i)-u(il,i-1)))
          fv(il, i) = fv(il, i) + 0.01*grav*dpinv*(amp1(il)*(v(il,i+1)-v(il, &
            i))-ad(il)*(v(il,i)-v(il,i-1)))
        ELSE ! cvflag_grav
          fr(il, i) = 0.1*dpinv*(amp1(il)*(rr(il,i+1)-rr(il, &
            i))-ad(il)*(rr(il,i)-rr(il,i-1)))
          fu(il, i) = fu(il, i) + 0.1*dpinv*(amp1(il)*(u(il,i+1)-u(il, &
            i))-ad(il)*(u(il,i)-u(il,i-1)))
          fv(il, i) = fv(il, i) + 0.1*dpinv*(amp1(il)*(v(il,i+1)-v(il, &
            i))-ad(il)*(v(il,i)-v(il,i-1)))
        END IF ! cvflag_grav

      END IF ! i
    END DO

    ! do k=1,ntra
    ! do il=1,ncum
    ! if (i.le.inb(il)) then
    ! dpinv=1.0/(ph(il,i)-ph(il,i+1))
    ! cpinv=1.0/cpn(il,i)
    ! if (cvflag_grav) then
    ! ftra(il,i,k)=ftra(il,i,k)+0.01*grav*dpinv
    ! :         *(amp1(il)*(tra(il,i+1,k)-tra(il,i,k))
    ! :           -ad(il)*(tra(il,i,k)-tra(il,i-1,k)))
    ! else
    ! ftra(il,i,k)=ftra(il,i,k)+0.1*dpinv
    ! :         *(amp1(il)*(tra(il,i+1,k)-tra(il,i,k))
    ! :           -ad(il)*(tra(il,i,k)-tra(il,i-1,k)))
    ! endif
    ! endif
    ! enddo
    ! enddo

    DO k = 1, i - 1
      DO il = 1, ncum
        IF (i<=inb(il)) THEN
          dpinv = 1.0/(ph(il,i)-ph(il,i+1))
          cpinv = 1.0/cpn(il, i)

          awat = elij(il, k, i) - (1.-ep(il,i))*clw(il, i)
          awat = amax1(awat, 0.0)

          IF (cvflag_grav) THEN
            fr(il, i) = fr(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(qent(il,k &
              ,i)-awat-rr(il,i))
            fu(il, i) = fu(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(uent(il,k &
              ,i)-u(il,i))
            fv(il, i) = fv(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(vent(il,k &
              ,i)-v(il,i))
          ELSE ! cvflag_grav
            fr(il, i) = fr(il, i) + 0.1*dpinv*ment(il, k, i)*(qent(il,k,i)- &
              awat-rr(il,i))
            fu(il, i) = fu(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(uent(il,k &
              ,i)-u(il,i))
            fv(il, i) = fv(il, i) + 0.1*dpinv*ment(il, k, i)*(vent(il,k,i)-v( &
              il,i))
          END IF ! cvflag_grav

          ! (saturated updrafts resulting from mixing)        ! cld
          qcond(il, i) = qcond(il, i) + (elij(il,k,i)-awat) ! cld
          nqcond(il, i) = nqcond(il, i) + 1. ! cld
        END IF ! i
      END DO
    END DO

    ! do j=1,ntra
    ! do k=1,i-1
    ! do il=1,ncum
    ! if (i.le.inb(il)) then
    ! dpinv=1.0/(ph(il,i)-ph(il,i+1))
    ! cpinv=1.0/cpn(il,i)
    ! if (cvflag_grav) then
    ! ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv*ment(il,k,i)
    ! :        *(traent(il,k,i,j)-tra(il,i,j))
    ! else
    ! ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv*ment(il,k,i)
    ! :        *(traent(il,k,i,j)-tra(il,i,j))
    ! endif
    ! endif
    ! enddo
    ! enddo
    ! enddo

    DO k = i, nl + 1
      DO il = 1, ncum
        IF (i<=inb(il) .AND. k<=inb(il)) THEN
          dpinv = 1.0/(ph(il,i)-ph(il,i+1))
          cpinv = 1.0/cpn(il, i)

          IF (cvflag_grav) THEN
            fr(il, i) = fr(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(qent(il,k &
              ,i)-rr(il,i))
            fu(il, i) = fu(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(uent(il,k &
              ,i)-u(il,i))
            fv(il, i) = fv(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(vent(il,k &
              ,i)-v(il,i))
          ELSE ! cvflag_grav
            fr(il, i) = fr(il, i) + 0.1*dpinv*ment(il, k, i)*(qent(il,k,i)-rr &
              (il,i))
            fu(il, i) = fu(il, i) + 0.1*dpinv*ment(il, k, i)*(uent(il,k,i)-u( &
              il,i))
            fv(il, i) = fv(il, i) + 0.1*dpinv*ment(il, k, i)*(vent(il,k,i)-v( &
              il,i))
          END IF ! cvflag_grav
        END IF ! i and k
      END DO
    END DO

    ! do j=1,ntra
    ! do k=i,nl+1
    ! do il=1,ncum
    ! if (i.le.inb(il) .and. k.le.inb(il)) then
    ! dpinv=1.0/(ph(il,i)-ph(il,i+1))
    ! cpinv=1.0/cpn(il,i)
    ! if (cvflag_grav) then
    ! ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv*ment(il,k,i)
    ! :         *(traent(il,k,i,j)-tra(il,i,j))
    ! else
    ! ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv*ment(il,k,i)
    ! :             *(traent(il,k,i,j)-tra(il,i,j))
    ! endif
    ! endif ! i and k
    ! enddo
    ! enddo
    ! enddo

    DO il = 1, ncum
      IF (i<=inb(il)) THEN
        dpinv = 1.0/(ph(il,i)-ph(il,i+1))
        cpinv = 1.0/cpn(il, i)

        IF (cvflag_grav) THEN
          ! sb: on ne fait pas encore la correction permettant de mieux
          ! conserver l'eau:
          fr(il, i) = fr(il, i) + 0.5*sigd*(evap(il,i)+evap(il,i+1)) + &
            0.01*grav*(mp(il,i+1)*(rp(il,i+1)-rr(il,i))-mp(il,i)*(rp(il, &
            i)-rr(il,i-1)))*dpinv

          fu(il, i) = fu(il, i) + 0.01*grav*(mp(il,i+1)*(up(il,i+1)-u(il, &
            i))-mp(il,i)*(up(il,i)-u(il,i-1)))*dpinv
          fv(il, i) = fv(il, i) + 0.01*grav*(mp(il,i+1)*(vp(il,i+1)-v(il, &
            i))-mp(il,i)*(vp(il,i)-v(il,i-1)))*dpinv
        ELSE ! cvflag_grav
          fr(il, i) = fr(il, i) + 0.5*sigd*(evap(il,i)+evap(il,i+1)) + &
            0.1*(mp(il,i+1)*(rp(il,i+1)-rr(il,i))-mp(il,i)*(rp(il,i)-rr(il, &
            i-1)))*dpinv
          fu(il, i) = fu(il, i) + 0.1*(mp(il,i+1)*(up(il,i+1)-u(il, &
            i))-mp(il,i)*(up(il,i)-u(il,i-1)))*dpinv
          fv(il, i) = fv(il, i) + 0.1*(mp(il,i+1)*(vp(il,i+1)-v(il, &
            i))-mp(il,i)*(vp(il,i)-v(il,i-1)))*dpinv
        END IF ! cvflag_grav

      END IF ! i
    END DO

    ! sb: interface with the cloud parameterization:          ! cld

    DO k = i + 1, nl
      DO il = 1, ncum
        IF (k<=inb(il) .AND. i<=inb(il)) THEN ! cld
          ! (saturated downdrafts resulting from mixing)            ! cld
          qcond(il, i) = qcond(il, i) + elij(il, k, i) ! cld
          nqcond(il, i) = nqcond(il, i) + 1. ! cld
        END IF ! cld
      END DO ! cld
    END DO ! cld

    ! (particular case: no detraining level is found)         ! cld
    DO il = 1, ncum ! cld
      IF (i<=inb(il) .AND. nent(il,i)==0) THEN ! cld
        qcond(il, i) = qcond(il, i) + (1.-ep(il,i))*clw(il, i) ! cld
        nqcond(il, i) = nqcond(il, i) + 1. ! cld
      END IF ! cld
    END DO ! cld

    DO il = 1, ncum ! cld
      IF (i<=inb(il) .AND. nqcond(il,i)/=0.) THEN ! cld
        qcond(il, i) = qcond(il, i)/nqcond(il, i) ! cld
      END IF ! cld
    END DO

    ! do j=1,ntra
    ! do il=1,ncum
    ! if (i.le.inb(il)) then
    ! dpinv=1.0/(ph(il,i)-ph(il,i+1))
    ! cpinv=1.0/cpn(il,i)

    ! if (cvflag_grav) then
    ! ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv
    ! :     *(mp(il,i+1)*(trap(il,i+1,j)-tra(il,i,j))
    ! :     -mp(il,i)*(trap(il,i,j)-tra(il,i-1,j)))
    ! else
    ! ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv
    ! :     *(mp(il,i+1)*(trap(il,i+1,j)-tra(il,i,j))
    ! :     -mp(il,i)*(trap(il,i,j)-tra(il,i-1,j)))
    ! endif
    ! endif ! i
    ! enddo
    ! enddo

500 END DO


  ! ***   move the detrainment at level inb down to level inb-1   ***
  ! ***        in such a way as to preserve the vertically        ***
  ! ***          integrated enthalpy and water tendencies         ***

  DO il = 1, ncum

    ax = 0.1*ment(il, inb(il), inb(il))*(hp(il,inb(il))-h(il,inb(il))+t(il, &
      inb(il))*(cpv-cpd)*(rr(il,inb(il))-qent(il,inb(il), &
      inb(il))))/(cpn(il,inb(il))*(ph(il,inb(il))-ph(il,inb(il)+1)))
    ft(il, inb(il)) = ft(il, inb(il)) - ax
    ft(il, inb(il)-1) = ft(il, inb(il)-1) + ax*cpn(il, inb(il))*(ph(il,inb(il &
      ))-ph(il,inb(il)+1))/(cpn(il,inb(il)-1)*(ph(il,inb(il)-1)-ph(il, &
      inb(il))))

    bx = 0.1*ment(il, inb(il), inb(il))*(qent(il,inb(il),inb(il))-rr(il,inb( &
      il)))/(ph(il,inb(il))-ph(il,inb(il)+1))
    fr(il, inb(il)) = fr(il, inb(il)) - bx
    fr(il, inb(il)-1) = fr(il, inb(il)-1) + bx*(ph(il,inb(il))-ph(il,inb(il)+ &
      1))/(ph(il,inb(il)-1)-ph(il,inb(il)))

    cx = 0.1*ment(il, inb(il), inb(il))*(uent(il,inb(il),inb(il))-u(il,inb(il &
      )))/(ph(il,inb(il))-ph(il,inb(il)+1))
    fu(il, inb(il)) = fu(il, inb(il)) - cx
    fu(il, inb(il)-1) = fu(il, inb(il)-1) + cx*(ph(il,inb(il))-ph(il,inb(il)+ &
      1))/(ph(il,inb(il)-1)-ph(il,inb(il)))

    dx = 0.1*ment(il, inb(il), inb(il))*(vent(il,inb(il),inb(il))-v(il,inb(il &
      )))/(ph(il,inb(il))-ph(il,inb(il)+1))
    fv(il, inb(il)) = fv(il, inb(il)) - dx
    fv(il, inb(il)-1) = fv(il, inb(il)-1) + dx*(ph(il,inb(il))-ph(il,inb(il)+ &
      1))/(ph(il,inb(il)-1)-ph(il,inb(il)))

  END DO

  ! do j=1,ntra
  ! do il=1,ncum
  ! ex=0.1*ment(il,inb(il),inb(il))
  ! :      *(traent(il,inb(il),inb(il),j)-tra(il,inb(il),j))
  ! :      /(ph(il,inb(il))-ph(il,inb(il)+1))
  ! ftra(il,inb(il),j)=ftra(il,inb(il),j)-ex
  ! ftra(il,inb(il)-1,j)=ftra(il,inb(il)-1,j)
  ! :       +ex*(ph(il,inb(il))-ph(il,inb(il)+1))
  ! :          /(ph(il,inb(il)-1)-ph(il,inb(il)))
  ! enddo
  ! enddo


  ! ***    homoginize tendencies below cloud base    ***


  DO il = 1, ncum
    asum(il) = 0.0
    bsum(il) = 0.0
    csum(il) = 0.0
    dsum(il) = 0.0
  END DO

  DO i = 1, nl
    DO il = 1, ncum
      IF (i<=(icb(il)-1)) THEN
        asum(il) = asum(il) + ft(il, i)*(ph(il,i)-ph(il,i+1))
        bsum(il) = bsum(il) + fr(il, i)*(lv(il,i)+(cl-cpd)*(t(il,i)-t(il, &
          1)))*(ph(il,i)-ph(il,i+1))
        csum(il) = csum(il) + (lv(il,i)+(cl-cpd)*(t(il,i)-t(il, &
          1)))*(ph(il,i)-ph(il,i+1))
        dsum(il) = dsum(il) + t(il, i)*(ph(il,i)-ph(il,i+1))/th(il, i)
      END IF
    END DO
  END DO

  ! !!!      do 700 i=1,icb(il)-1
  DO i = 1, nl
    DO il = 1, ncum
      IF (i<=(icb(il)-1)) THEN
        ft(il, i) = asum(il)*t(il, i)/(th(il,i)*dsum(il))
        fr(il, i) = bsum(il)/csum(il)
      END IF
    END DO
  END DO


  ! ***           reset counter and return           ***

  DO il = 1, ncum
    sig(il, nd) = 2.0
  END DO


  DO i = 1, nd
    DO il = 1, ncum
      upwd(il, i) = 0.0
      dnwd(il, i) = 0.0
    END DO
  END DO

  DO i = 1, nl
    DO il = 1, ncum
      dnwd0(il, i) = -mp(il, i)
    END DO
  END DO
  DO i = nl + 1, nd
    DO il = 1, ncum
      dnwd0(il, i) = 0.
    END DO
  END DO


  DO i = 1, nl
    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il)) THEN
        upwd(il, i) = 0.0
        dnwd(il, i) = 0.0
      END IF
    END DO
  END DO

  DO i = 1, nl
    DO k = 1, nl
      DO il = 1, ncum
        up1(il, k, i) = 0.0
        dn1(il, k, i) = 0.0
      END DO
    END DO
  END DO

  DO i = 1, nl
    DO k = i, nl
      DO n = 1, i - 1
        DO il = 1, ncum
          IF (i>=icb(il) .AND. i<=inb(il) .AND. k<=inb(il)) THEN
            up1(il, k, i) = up1(il, k, i) + ment(il, n, k)
            dn1(il, k, i) = dn1(il, k, i) - ment(il, k, n)
          END IF
        END DO
      END DO
    END DO
  END DO

  DO i = 2, nl
    DO k = i, nl
      DO il = 1, ncum
        ! test         if (i.ge.icb(il).and.i.le.inb(il).and.k.le.inb(il))
        ! then
        IF (i<=inb(il) .AND. k<=inb(il)) THEN
          upwd(il, i) = upwd(il, i) + m(il, k) + up1(il, k, i)
          dnwd(il, i) = dnwd(il, i) + dn1(il, k, i)
        END IF
      END DO
    END DO
  END DO


  ! !!!      DO il=1,ncum
  ! !!!      do i=icb(il),inb(il)
  ! !!!
  ! !!!      upwd(il,i)=0.0
  ! !!!      dnwd(il,i)=0.0
  ! !!!      do k=i,inb(il)
  ! !!!      up1=0.0
  ! !!!      dn1=0.0
  ! !!!      do n=1,i-1
  ! !!!      up1=up1+ment(il,n,k)
  ! !!!      dn1=dn1-ment(il,k,n)
  ! !!!      enddo
  ! !!!      upwd(il,i)=upwd(il,i)+m(il,k)+up1
  ! !!!      dnwd(il,i)=dnwd(il,i)+dn1
  ! !!!      enddo
  ! !!!      enddo
  ! !!!
  ! !!!      ENDDO

  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! determination de la variation de flux ascendant entre
  ! deux niveau non dilue mike
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  DO i = 1, nl
    DO il = 1, ncum
      mike(il, i) = m(il, i)
    END DO
  END DO

  DO i = nl + 1, nd
    DO il = 1, ncum
      mike(il, i) = 0.
    END DO
  END DO

  DO i = 1, nd
    DO il = 1, ncum
      ma(il, i) = 0
    END DO
  END DO

  DO i = 1, nl
    DO j = i, nl
      DO il = 1, ncum
        ma(il, i) = ma(il, i) + m(il, j)
      END DO
    END DO
  END DO

  DO i = nl + 1, nd
    DO il = 1, ncum
      ma(il, i) = 0.
    END DO
  END DO

  DO i = 1, nl
    DO il = 1, ncum
      IF (i<=(icb(il)-1)) THEN
        ma(il, i) = 0
      END IF
    END DO
  END DO

  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! icb represente de niveau ou se trouve la
  ! base du nuage , et inb le top du nuage
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  DO i = 1, nd
    DO il = 1, ncum
      mke(il, i) = upwd(il, i) + dnwd(il, i)
    END DO
  END DO

  DO i = 1, nd
    DO il = 1, ncum
      rdcp = (rrd*(1.-rr(il,i))-rr(il,i)*rrv)/(cpd*(1.-rr(il, &
        i))+rr(il,i)*cpv)
      tls(il, i) = t(il, i)*(1000.0/p(il,i))**rdcp
      tps(il, i) = tp(il, i)
    END DO
  END DO


  ! *** diagnose the in-cloud mixing ratio   ***            ! cld
  ! ***           of condensed water         ***            ! cld
  ! ! cld

  DO i = 1, nd ! cld
    DO il = 1, ncum ! cld
      mac(il, i) = 0.0 ! cld
      wa(il, i) = 0.0 ! cld
      siga(il, i) = 0.0 ! cld
      sax(il, i) = 0.0 ! cld
    END DO ! cld
  END DO ! cld

  DO i = minorig, nl ! cld
    DO k = i + 1, nl + 1 ! cld
      DO il = 1, ncum ! cld
        IF (i<=inb(il) .AND. k<=(inb(il)+1)) THEN ! cld
          mac(il, i) = mac(il, i) + m(il, k) ! cld
        END IF ! cld
      END DO ! cld
    END DO ! cld
  END DO ! cld

  DO i = 1, nl ! cld
    DO j = 1, i ! cld
      DO il = 1, ncum ! cld
        IF (i>=icb(il) .AND. i<=(inb(il)-1) & ! cld
            .AND. j>=icb(il)) THEN ! cld
          sax(il, i) = sax(il, i) + rrd*(tvp(il,j)-tv(il,j)) & ! cld
            *(ph(il,j)-ph(il,j+1))/p(il, j) ! cld
        END IF ! cld
      END DO ! cld
    END DO ! cld
  END DO ! cld

  DO i = 1, nl ! cld
    DO il = 1, ncum ! cld
      IF (i>=icb(il) .AND. i<=(inb(il)-1) & ! cld
          .AND. sax(il,i)>0.0) THEN ! cld
        wa(il, i) = sqrt(2.*sax(il,i)) ! cld
      END IF ! cld
    END DO ! cld
  END DO ! cld

  DO i = 1, nl ! cld
    DO il = 1, ncum ! cld
      IF (wa(il,i)>0.0) &          ! cld
        siga(il, i) = mac(il, i)/wa(il, i) & ! cld
        *rrd*tvp(il, i)/p(il, i)/100./delta ! cld
      siga(il, i) = min(siga(il,i), 1.0) ! cld
      ! IM cf. FH
      IF (iflag_clw==0) THEN
        qcondc(il, i) = siga(il, i)*clw(il, i)*(1.-ep(il,i)) & ! cld
          +(1.-siga(il,i))*qcond(il, i) ! cld
      ELSE IF (iflag_clw==1) THEN
        qcondc(il, i) = qcond(il, i) ! cld
      END IF

    END DO ! cld
  END DO ! cld

  RETURN
END SUBROUTINE cv30_yield

! !RomP >>>
SUBROUTINE cv30_tracer(nloc, len, ncum, nd, na, ment, sij, da, phi, phi2, &
    d1a, dam, ep, vprecip, elij, clw, epmlmmm, eplamm, icb, inb)
  IMPLICIT NONE

  include "cv30param.h"

  ! inputs:
  INTEGER ncum, nd, na, nloc, len
  REAL ment(nloc, na, na), sij(nloc, na, na)
  REAL clw(nloc, nd), elij(nloc, na, na)
  REAL ep(nloc, na)
  INTEGER icb(nloc), inb(nloc)
  REAL vprecip(nloc, nd+1)
  ! ouputs:
  REAL da(nloc, na), phi(nloc, na, na)
  REAL phi2(nloc, na, na)
  REAL d1a(nloc, na), dam(nloc, na)
  REAL epmlmmm(nloc, na, na), eplamm(nloc, na)
  ! variables pour tracer dans precip de l'AA et des mel
  ! local variables:
  INTEGER i, j, k, nam1
  REAL epm(nloc, na, na)

  nam1=na-1 ! Introduced because ep is not defined for j=na
  ! variables d'Emanuel : du second indice au troisieme
  ! --->    tab(i,k,j) -> de l origine k a l arrivee j
  ! ment, sij, elij
  ! variables personnelles : du troisieme au second indice
  ! --->    tab(i,j,k) -> de k a j
  ! phi, phi2

  ! initialisations
  DO j = 1, na
    DO i = 1, ncum
      da(i, j) = 0.
      d1a(i, j) = 0.
      dam(i, j) = 0.
      eplamm(i, j) = 0.
    END DO
  END DO
  DO k = 1, na
    DO j = 1, na
      DO i = 1, ncum
        epm(i, j, k) = 0.
        epmlmmm(i, j, k) = 0.
        phi(i, j, k) = 0.
        phi2(i, j, k) = 0.
      END DO
    END DO
  END DO

  ! fraction deau condensee dans les melanges convertie en precip : epm
  ! et eau condensée précipitée dans masse d'air saturé : l_m*dM_m/dzdz.dzdz
  DO j = 1, nam1
    DO k = 1, j - 1
      DO i = 1, ncum
        IF (k>=icb(i) .AND. k<=inb(i) .AND. j<=inb(i)) THEN
          ! !jyg             epm(i,j,k)=1.-(1.-ep(i,j))*clw(i,j)/elij(i,k,j)
          epm(i, j, k) = 1. - (1.-ep(i,j))*clw(i, j)/max(elij(i,k,j), 1.E-16)
          ! !
          epm(i, j, k) = max(epm(i,j,k), 0.0)
        END IF
      END DO
    END DO
  END DO

  DO j = 1, nam1
    DO k = 1, nam1
      DO i = 1, ncum
        IF (k>=icb(i) .AND. k<=inb(i)) THEN
          eplamm(i, j) = eplamm(i, j) + ep(i, j)*clw(i, j)*ment(i, j, k)*(1.- &
            sij(i,j,k))
        END IF
      END DO
    END DO
  END DO

  DO j = 1, nam1
    DO k = 1, j - 1
      DO i = 1, ncum
        IF (k>=icb(i) .AND. k<=inb(i) .AND. j<=inb(i)) THEN
          epmlmmm(i, j, k) = epm(i, j, k)*elij(i, k, j)*ment(i, k, j)
        END IF
      END DO
    END DO
  END DO

  ! matrices pour calculer la tendance des concentrations dans cvltr.F90
  DO j = 1, nam1
    DO k = 1, nam1
      DO i = 1, ncum
        da(i, j) = da(i, j) + (1.-sij(i,k,j))*ment(i, k, j)
        phi(i, j, k) = sij(i, k, j)*ment(i, k, j)
        d1a(i, j) = d1a(i, j) + ment(i, k, j)*ep(i, k)*(1.-sij(i,k,j))
      END DO
    END DO
  END DO

  DO j = 1, nam1
    DO k = 1, j - 1
      DO i = 1, ncum
        dam(i, j) = dam(i, j) + ment(i, k, j)*epm(i, j, k)*(1.-ep(i,k))*(1.- &
          sij(i,k,j))
        phi2(i, j, k) = phi(i, j, k)*epm(i, j, k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cv30_tracer
! RomP <<<

SUBROUTINE cv30_uncompress(nloc, len, ncum, nd, ntra, idcum, iflag, precip, &
    vprecip, evap, ep, sig, w0, ft, fq, fu, fv, ftra, inb, ma, upwd, dnwd, &
    dnwd0, qcondc, wd, cape, da, phi, mp, phi2, d1a, dam, sij, elij, clw, &
    epmlmmm, eplamm, wdtraina, wdtrainm,epmax_diag, iflag1, precip1, vprecip1, evap1, &
    ep1, sig1, w01, ft1, fq1, fu1, fv1, ftra1, inb1, ma1, upwd1, dnwd1, &
    dnwd01, qcondc1, wd1, cape1, da1, phi1, mp1, phi21, d1a1, dam1, sij1, &
    elij1, clw1, epmlmmm1, eplamm1, wdtraina1, wdtrainm1,epmax_diag1) ! epmax_cape
  IMPLICIT NONE

  include "cv30param.h"

  ! inputs:
  INTEGER len, ncum, nd, ntra, nloc
  INTEGER idcum(nloc)
  INTEGER iflag(nloc)
  INTEGER inb(nloc)
  REAL precip(nloc)
  REAL vprecip(nloc, nd+1), evap(nloc, nd)
  REAL ep(nloc, nd)
  REAL sig(nloc, nd), w0(nloc, nd)
  REAL ft(nloc, nd), fq(nloc, nd), fu(nloc, nd), fv(nloc, nd)
  REAL ftra(nloc, nd, ntra)
  REAL ma(nloc, nd)
  REAL upwd(nloc, nd), dnwd(nloc, nd), dnwd0(nloc, nd)
  REAL qcondc(nloc, nd)
  REAL wd(nloc), cape(nloc)
  REAL da(nloc, nd), phi(nloc, nd, nd), mp(nloc, nd)
  REAL epmax_diag(nloc) ! epmax_cape
  ! RomP >>>
  REAL phi2(nloc, nd, nd)
  REAL d1a(nloc, nd), dam(nloc, nd)
  REAL wdtraina(nloc, nd), wdtrainm(nloc, nd)
  REAL sij(nloc, nd, nd)
  REAL elij(nloc, nd, nd), clw(nloc, nd)
  REAL epmlmmm(nloc, nd, nd), eplamm(nloc, nd)
  ! RomP <<<

  ! outputs:
  INTEGER iflag1(len)
  INTEGER inb1(len)
  REAL precip1(len)
  REAL vprecip1(len, nd+1), evap1(len, nd) !<<< RomP
  REAL ep1(len, nd) !<<< RomP
  REAL sig1(len, nd), w01(len, nd)
  REAL ft1(len, nd), fq1(len, nd), fu1(len, nd), fv1(len, nd)
  REAL ftra1(len, nd, ntra)
  REAL ma1(len, nd)
  REAL upwd1(len, nd), dnwd1(len, nd), dnwd01(len, nd)
  REAL qcondc1(nloc, nd)
  REAL wd1(nloc), cape1(nloc)
  REAL da1(nloc, nd), phi1(nloc, nd, nd), mp1(nloc, nd)
  REAL epmax_diag1(len) ! epmax_cape
  ! RomP >>>
  REAL phi21(len, nd, nd)
  REAL d1a1(len, nd), dam1(len, nd)
  REAL wdtraina1(len, nd), wdtrainm1(len, nd)
  REAL sij1(len, nd, nd)
  REAL elij1(len, nd, nd), clw1(len, nd)
  REAL epmlmmm1(len, nd, nd), eplamm1(len, nd)
  ! RomP <<<

  ! local variables:
  INTEGER i, k, j

  DO i = 1, ncum
    precip1(idcum(i)) = precip(i)
    iflag1(idcum(i)) = iflag(i)
    wd1(idcum(i)) = wd(i)
    inb1(idcum(i)) = inb(i)
    cape1(idcum(i)) = cape(i)
    epmax_diag1(idcum(i))=epmax_diag(i) ! epmax_cape
  END DO

  DO k = 1, nl
    DO i = 1, ncum
      vprecip1(idcum(i), k) = vprecip(i, k)
      evap1(idcum(i), k) = evap(i, k) !<<< RomP
      sig1(idcum(i), k) = sig(i, k)
      w01(idcum(i), k) = w0(i, k)
      ft1(idcum(i), k) = ft(i, k)
      fq1(idcum(i), k) = fq(i, k)
      fu1(idcum(i), k) = fu(i, k)
      fv1(idcum(i), k) = fv(i, k)
      ma1(idcum(i), k) = ma(i, k)
      upwd1(idcum(i), k) = upwd(i, k)
      dnwd1(idcum(i), k) = dnwd(i, k)
      dnwd01(idcum(i), k) = dnwd0(i, k)
      qcondc1(idcum(i), k) = qcondc(i, k)
      da1(idcum(i), k) = da(i, k)
      mp1(idcum(i), k) = mp(i, k)
      ! RomP >>>
      ep1(idcum(i), k) = ep(i, k)
      d1a1(idcum(i), k) = d1a(i, k)
      dam1(idcum(i), k) = dam(i, k)
      clw1(idcum(i), k) = clw(i, k)
      eplamm1(idcum(i), k) = eplamm(i, k)
      wdtraina1(idcum(i), k) = wdtraina(i, k)
      wdtrainm1(idcum(i), k) = wdtrainm(i, k)
      ! RomP <<<
    END DO
  END DO

  DO i = 1, ncum
    sig1(idcum(i), nd) = sig(i, nd)
  END DO


  ! do 2100 j=1,ntra
  ! do 2110 k=1,nd ! oct3
  ! do 2120 i=1,ncum
  ! ftra1(idcum(i),k,j)=ftra(i,k,j)
  ! 2120     continue
  ! 2110    continue
  ! 2100   continue
  DO j = 1, nd
    DO k = 1, nd
      DO i = 1, ncum
        sij1(idcum(i), k, j) = sij(i, k, j)
        phi1(idcum(i), k, j) = phi(i, k, j)
        phi21(idcum(i), k, j) = phi2(i, k, j)
        elij1(idcum(i), k, j) = elij(i, k, j)
        epmlmmm1(idcum(i), k, j) = epmlmmm(i, k, j)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cv30_uncompress

        subroutine cv30_epmax_fn_cape(nloc,ncum,nd &
                ,cape,ep,hp,icb,inb,clw,nk,t,h,lv &
                ,epmax_diag)
        implicit none

        ! On fait varier epmax en fn de la cape
        ! Il faut donc recalculer ep, et hp qui a déjà été calculé et
        ! qui en dépend
        ! Toutes les autres variables fn de ep sont calculées plus bas.

#include "cvthermo.h"
#include "cv30param.h"
#include "conema3.h"

! inputs:
      integer ncum, nd, nloc
      integer icb(nloc), inb(nloc)
      real cape(nloc)
      real clw(nloc,nd),lv(nloc,nd),t(nloc,nd),h(nloc,nd)
      integer nk(nloc)
! inouts:
      real ep(nloc,nd)
      real hp(nloc,nd)
! outputs ou local
      real epmax_diag(nloc)
! locals
      integer i,k    
      real hp_bak(nloc,nd)
      CHARACTER (LEN=20) :: modname='cv30_epmax_fn_cape'
      CHARACTER (LEN=80) :: abort_message

        ! on recalcule ep et hp
        
        if (coef_epmax_cape.gt.1e-12) then
        do i=1,ncum
           epmax_diag(i)=epmax-coef_epmax_cape*sqrt(cape(i))
           do k=1,nl
                ep(i,k)=ep(i,k)/epmax*epmax_diag(i)
                ep(i,k)=amax1(ep(i,k),0.0)
                ep(i,k)=amin1(ep(i,k),epmax_diag(i))
           enddo
        enddo

! On recalcule hp:
      do k=1,nl
        do i=1,ncum
	  hp_bak(i,k)=hp(i,k)
	enddo
      enddo
      do k=1,nlp
        do i=1,ncum
	  hp(i,k)=h(i,k)
	enddo
      enddo
      do k=minorig+1,nl
       do i=1,ncum
        if((k.ge.icb(i)).and.(k.le.inb(i)))then
          hp(i,k)=h(i,nk(i))+(lv(i,k)+(cpd-cpv)*t(i,k))*ep(i,k)*clw(i,k)
        endif
       enddo
      enddo !do k=minorig+1,n
!     write(*,*) 'cv30_routines 6218: hp(1,20)=',hp(1,20)
      do i=1,ncum  
       do k=1,nl
        if (abs(hp_bak(i,k)-hp(i,k)).gt.0.01) then
           write(*,*) 'i,k=',i,k
           write(*,*) 'coef_epmax_cape=',coef_epmax_cape
           write(*,*) 'epmax_diag(i)=',epmax_diag(i)
           write(*,*) 'ep(i,k)=',ep(i,k)
           write(*,*) 'hp(i,k)=',hp(i,k)
           write(*,*) 'hp_bak(i,k)=',hp_bak(i,k)
           write(*,*) 'h(i,k)=',h(i,k)
           write(*,*) 'nk(i)=',nk(i)
           write(*,*) 'h(i,nk(i))=',h(i,nk(i))
           write(*,*) 'lv(i,k)=',lv(i,k)
           write(*,*) 't(i,k)=',t(i,k)
           write(*,*) 'clw(i,k)=',clw(i,k)
           write(*,*) 'cpd,cpv=',cpd,cpv
           CALL abort_physic(modname,abort_message,0)
        endif
       enddo !do k=1,nl
      enddo !do i=1,ncum  
      endif !if (coef_epmax_cape.gt.1e-12) then

      return
      end subroutine cv30_epmax_fn_cape


