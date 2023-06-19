
! $Header$

SUBROUTINE convect1(len, nd, ndp1, noff, minorig, t, q, qs, u, v, p, ph, &
    iflag, ft, fq, fu, fv, precip, cbmf, delt, ma)
  ! .............................START PROLOGUE............................

  ! SCCS IDENTIFICATION:  @(#)convect1.f	1.1 04/21/00
  ! 19:40:52 /h/cm/library/nogaps4/src/sub/fcst/convect1.f_v

  ! CONFIGURATION IDENTIFICATION:  None

  ! MODULE NAME:  convect1

  ! DESCRIPTION:

  ! convect1     The Emanuel Cumulus Convection Scheme

  ! CONTRACT NUMBER AND TITLE:  None

  ! REFERENCES: Programmers  K. Emanuel (MIT), Timothy F. Hogan, M. Peng
  ! (NRL)

  ! CLASSIFICATION:  Unclassified

  ! RESTRICTIONS: None

  ! COMPILER DEPENDENCIES: FORTRAN 77, FORTRAN 90

  ! COMPILE OPTIONS: Fortran 77: -Zu -Wf"-ei -o aggress"
  ! Fortran 90: -O vector3,scalar3,task1,aggress,overindex  -ei -r 2

  ! LIBRARIES OF RESIDENCE: /a/ops/lib/libfcst159.a

  ! USAGE: call convect1(len,nd,noff,minorig,
  ! &                   t,q,qs,u,v,
  ! &                   p,ph,iflag,ft,
  ! &                   fq,fu,fv,precip,cbmf,delt)

  ! PARAMETERS:
  ! Name            Type         Usage            Description
  ! ----------      ----------     -------  ----------------------------

  ! len           Integer        Input        first (i) dimension
  ! nd            Integer        Input        vertical (k) dimension
  ! ndp1          Integer        Input        nd + 1
  ! noff          Integer        Input        integer limit for convection
  ! (nd-noff)
  ! minorig       Integer        Input        First level of convection
  ! t             Real           Input        temperature
  ! q             Real           Input        specific hum
  ! qs            Real           Input        sat specific hum
  ! u             Real           Input        u-wind
  ! v             Real           Input        v-wind
  ! p             Real           Input        full level pressure
  ! ph            Real           Input        half level pressure
  ! iflag         Integer        Output       iflag on latitude strip
  ! ft            Real           Output       temp tend
  ! fq            Real           Output       spec hum tend
  ! fu            Real           Output       u-wind tend
  ! fv            Real           Output       v-wind tend
  ! cbmf          Real           In/Out       cumulus mass flux
  ! delt          Real           Input        time step
  ! iflag         Integer        Output       integer flag for Emanuel
  ! conditions

  ! COMMON BLOCKS:
  ! Block      Name     Type    Usage              Notes
  ! --------  --------   ----    ------   ------------------------

  ! FILES: None

  ! DATA BASES: None

  ! NON-FILE INPUT/OUTPUT: None

  ! ERROR CONDITIONS: None

  ! ADDITIONAL COMMENTS: None

  ! .................MAINTENANCE SECTION................................

  ! MODULES CALLED:
  ! Name           Description
  ! convect2        Emanuel cumulus convection tendency calculations
  ! -------     ----------------------
  ! LOCAL VARIABLES AND
  ! STRUCTURES:
  ! Name     Type    Description
  ! -------  ------  -----------
  ! See Comments Below

  ! i        Integer loop index
  ! k        Integer loop index

  ! METHOD:

  ! See Emanuel, K. and M. Zivkovic-Rothman, 2000: Development and evaluation
  ! of a
  ! convective scheme for use in climate models.

  ! FILES: None

  ! INCLUDE FILES: None

  ! MAKEFILE: /a/ops/met/nogaps/src/sub/fcst/fcst159lib.mak

  ! ..............................END PROLOGUE.............................


  USE dimphy
  IMPLICIT NONE

  INTEGER len
  INTEGER nd
  INTEGER ndp1
  INTEGER noff
  REAL t(len, nd)
  REAL q(len, nd)
  REAL qs(len, nd)
  REAL u(len, nd)
  REAL v(len, nd)
  REAL p(len, nd)
  REAL ph(len, ndp1)
  INTEGER iflag(len)
  REAL ft(len, nd)
  REAL fq(len, nd)
  REAL fu(len, nd)
  REAL fv(len, nd)
  REAL precip(len)
  REAL cbmf(len)
  REAL ma(len, nd)
  INTEGER minorig
  REAL delt, cpd, cpv, cl, rv, rd, lv0, g
  REAL sigs, sigd, elcrit, tlcrit, omtsnow, dtmax, damp
  REAL alpha, entp, coeffs, coeffr, omtrain, cu

  ! -------------------------------------------------------------------
  ! --- ARGUMENTS
  ! -------------------------------------------------------------------
  ! --- On input:

  ! t:   Array of absolute temperature (K) of dimension ND, with first
  ! index corresponding to lowest model level. Note that this array
  ! will be altered by the subroutine if dry convective adjustment
  ! occurs and if IPBL is not equal to 0.

  ! q:   Array of specific humidity (gm/gm) of dimension ND, with first
  ! index corresponding to lowest model level. Must be defined
  ! at same grid levels as T. Note that this array will be altered
  ! if dry convective adjustment occurs and if IPBL is not equal to 0.

  ! qs:  Array of saturation specific humidity of dimension ND, with first
  ! index corresponding to lowest model level. Must be defined
  ! at same grid levels as T. Note that this array will be altered
  ! if dry convective adjustment occurs and if IPBL is not equal to 0.

  ! u:   Array of zonal wind velocity (m/s) of dimension ND, witth first
  ! index corresponding with the lowest model level. Defined at
  ! same levels as T. Note that this array will be altered if
  ! dry convective adjustment occurs and if IPBL is not equal to 0.

  ! v:   Same as u but for meridional velocity.

  ! tra: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
  ! where NTRA is the number of different tracers. If no
  ! convective tracer transport is needed, define a dummy
  ! input array of dimension (ND,1). Tracers are defined at
  ! same vertical levels as T. Note that this array will be altered
  ! if dry convective adjustment occurs and if IPBL is not equal to 0.

  ! p:   Array of pressure (mb) of dimension ND, with first
  ! index corresponding to lowest model level. Must be defined
  ! at same grid levels as T.

  ! ph:  Array of pressure (mb) of dimension ND+1, with first index
  ! corresponding to lowest level. These pressures are defined at
  ! levels intermediate between those of P, T, Q and QS. The first
  ! value of PH should be greater than (i.e. at a lower level than)
  ! the first value of the array P.

  ! nl:  The maximum number of levels to which convection can penetrate, plus
  ! 1.
  ! NL MUST be less than or equal to ND-1.

  ! delt: The model time step (sec) between calls to CONVECT

  ! ----------------------------------------------------------------------------
  ! ---   On Output:

  ! iflag: An output integer whose value denotes the following:
  ! VALUE   INTERPRETATION
  ! -----   --------------
  ! 0     Moist convection occurs.
  ! 1     Moist convection occurs, but a CFL condition
  ! on the subsidence warming is violated. This
  ! does not cause the scheme to terminate.
  ! 2     Moist convection, but no precip because ep(inb) lt 0.0001
  ! 3     No moist convection because new cbmf is 0 and old cbmf is 0.
  ! 4     No moist convection; atmosphere is not
  ! unstable
  ! 6     No moist convection because ihmin le minorig.
  ! 7     No moist convection because unreasonable
  ! parcel level temperature or specific humidity.
  ! 8     No moist convection: lifted condensation
  ! level is above the 200 mb level.
  ! 9     No moist convection: cloud base is higher
  ! then the level NL-1.

  ! ft:   Array of temperature tendency (K/s) of dimension ND, defined at
  ! same
  ! grid levels as T, Q, QS and P.

  ! fq:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
  ! defined at same grid levels as T, Q, QS and P.

  ! fu:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
  ! defined at same grid levels as T.

  ! fv:   Same as FU, but for forcing of meridional velocity.

  ! ftra: Array of forcing of tracer content, in tracer mixing ratio per
  ! second, defined at same levels as T. Dimensioned (ND,NTRA).

  ! precip: Scalar convective precipitation rate (mm/day).

  ! wd:   A convective downdraft velocity scale. For use in surface
  ! flux parameterizations. See convect.ps file for details.

  ! tprime: A convective downdraft temperature perturbation scale (K).
  ! For use in surface flux parameterizations. See convect.ps
  ! file for details.

  ! qprime: A convective downdraft specific humidity
  ! perturbation scale (gm/gm).
  ! For use in surface flux parameterizations. See convect.ps
  ! file for details.

  ! cbmf: The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
  ! BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
  ! ITS NEXT CALL. That is, the value of CBMF must be "remembered"
  ! by the calling program between calls to CONVECT.

  ! det:   Array of detrainment mass flux of dimension ND.

  ! -------------------------------------------------------------------

  ! Local arrays

  INTEGER nl
  INTEGER nlp
  INTEGER nlm
  INTEGER i, k, n
  REAL delti
  REAL rowl
  REAL clmcpv
  REAL clmcpd
  REAL cpdmcp
  REAL cpvmcpd
  REAL eps
  REAL epsi
  REAL epsim1
  REAL ginv
  REAL hrd
  REAL prccon1
  INTEGER icbmax
  REAL lv(klon, klev)
  REAL cpn(klon, klev)
  REAL cpx(klon, klev)
  REAL tv(klon, klev)
  REAL gz(klon, klev)
  REAL hm(klon, klev)
  REAL h(klon, klev)
  REAL work(klon)
  INTEGER ihmin(klon)
  INTEGER nk(klon)
  REAL rh(klon)
  REAL chi(klon)
  REAL plcl(klon)
  INTEGER icb(klon)
  REAL tnk(klon)
  REAL qnk(klon)
  REAL gznk(klon)
  REAL pnk(klon)
  REAL qsnk(klon)
  REAL ticb(klon)
  REAL gzicb(klon)
  REAL tp(klon, klev)
  REAL tvp(klon, klev)
  REAL clw(klon, klev)

  REAL ah0(klon), cpp(klon)
  REAL tg, qg, s, alv, tc, ahg, denom, es, rg

  INTEGER ncum
  INTEGER idcum(klon)

  cpd = 1005.7
  cpv = 1870.0
  cl = 4190.0
  rv = 461.5
  rd = 287.04
  lv0 = 2.501E6
  g = 9.8

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

  sigs = 0.12
  sigd = 0.05
  elcrit = 0.0011
  tlcrit = -55.0
  omtsnow = 5.5
  dtmax = 0.9
  damp = 0.1
  alpha = 0.2
  entp = 1.5
  coeffs = 0.8
  coeffr = 1.0
  omtrain = 50.0

  cu = 0.70
  damp = 0.1


  ! Define nl, nlp, nlm, and delti

  nl = nd - noff
  nlp = nl + 1
  nlm = nl - 1
  delti = 1.0/delt

  ! -------------------------------------------------------------------
  ! --- SET CONSTANTS
  ! -------------------------------------------------------------------

  rowl = 1000.0
  clmcpv = cl - cpv
  clmcpd = cl - cpd
  cpdmcp = cpd - cpv
  cpvmcpd = cpv - cpd
  eps = rd/rv
  epsi = 1.0/eps
  epsim1 = epsi - 1.0
  ginv = 1.0/g
  hrd = 0.5*rd
  prccon1 = 86400.0*1000.0/(rowl*g)

  ! dtmax is the maximum negative temperature perturbation.

  ! =====================================================================
  ! --- INITIALIZE OUTPUT ARRAYS AND PARAMETERS
  ! =====================================================================

  DO k = 1, nd
    DO i = 1, len
      ft(i, k) = 0.0
      fq(i, k) = 0.0
      fu(i, k) = 0.0
      fv(i, k) = 0.0
      tvp(i, k) = 0.0
      tp(i, k) = 0.0
      clw(i, k) = 0.0
      gz(i, k) = 0.
    END DO
  END DO
  DO i = 1, len
    precip(i) = 0.0
    iflag(i) = 0
  END DO

  ! =====================================================================
  ! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
  ! =====================================================================
  DO k = 1, nl + 1
    DO i = 1, len
      lv(i, k) = lv0 - clmcpv*(t(i,k)-273.15)
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
    alv = lv0 - clmcpv*(ticb(i)-273.15)

    ! First iteration.

    s = cpd + alv*alv*qg/(rv*ticb(i)*ticb(i))
    s = 1./s
    ahg = cpd*tg + (cl-cpd)*qnk(i)*ticb(i) + alv*qg + gzicb(i)
    tg = tg + s*(ah0(i)-ahg)
    tg = max(tg, 35.0)
    tc = tg - 273.15
    denom = 243.5 + tc
    IF (tc>=0.0) THEN
      es = 6.112*exp(17.67*tc/denom)
    ELSE
      es = exp(23.33086-6111.72784/tg+0.15215*log(tg))
    END IF
    qg = eps*es/(p(i,icb(i))-es*(1.-eps))

    ! Second iteration.

    s = cpd + alv*alv*qg/(rv*ticb(i)*ticb(i))
    s = 1./s
    ahg = cpd*tg + (cl-cpd)*qnk(i)*ticb(i) + alv*qg + gzicb(i)
    tg = tg + s*(ah0(i)-ahg)
    tg = max(tg, 35.0)
    tc = tg - 273.15
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

  ! -------------------------------------------------------------------
  ! --- Test for instability.
  ! --- If there was no convection at last time step and parcel
  ! --- is stable at icb, then set iflag to 4.
  ! -------------------------------------------------------------------

  DO i = 1, len
    IF ((cbmf(i)==0.0) .AND. (iflag(i)==0) .AND. (tvp(i, &
      icb(i))<=(tv(i,icb(i))-dtmax))) iflag(i) = 4
  END DO

  ! =====================================================================
  ! --- IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY
  ! =====================================================================

  ncum = 0
  DO i = 1, len
    IF (iflag(i)==0) THEN
      ncum = ncum + 1
      idcum(ncum) = i
    END IF
  END DO

  ! Call convect2, which compresses the points and computes the heating,
  ! moistening, velocity mixing, and precipiation.

  ! print*,'cpd avant convect2 ',cpd
  IF (ncum>0) THEN
    CALL convect2(ncum, idcum, len, nd, ndp1, nl, minorig, nk, icb, t, q, qs, &
      u, v, gz, tv, tp, tvp, clw, h, lv, cpn, p, ph, ft, fq, fu, fv, tnk, &
      qnk, gznk, plcl, precip, cbmf, iflag, delt, cpd, cpv, cl, rv, rd, lv0, &
      g, sigs, sigd, elcrit, tlcrit, omtsnow, dtmax, damp, alpha, entp, &
      coeffs, coeffr, omtrain, cu, ma)
  END IF

  RETURN
END SUBROUTINE convect1
