
! $Header$

SUBROUTINE cv_driver(len, nd, ndp1, ntra, iflag_con, t1, q1, qs1, u1, v1, &
    tra1, p1, ph1, iflag1, ft1, fq1, fu1, fv1, ftra1, precip1, vprecip1, &
    cbmf1, sig1, w01, icb1, inb1, delt, ma1, upwd1, dnwd1, dnwd01, qcondc1, &
    wd1, cape1, da1, phi1, mp1, phi21, d1a1, dam1, sij1, clw1, elij1, & !
                                                                        ! RomP
    evap1, ep1, epmlmmm1, eplamm1, & ! RomP
    wdtraina1, wdtrainm1, & ! RomP
    epmax_diag1) ! epmax_cape

  USE dimphy
  IMPLICIT NONE

  ! .............................START PROLOGUE............................


  ! All argument names (except len,nd,ntra,nloc,delt and the flags) have a
  ! "1" appended.
  ! The "1" is removed for the corresponding compressed (local) variables.

  ! PARAMETERS:
  ! Name            Type         Usage            Description
  ! ----------      ----------     -------  ----------------------------

  ! len           Integer        Input        first (i) dimension
  ! nd            Integer        Input        vertical (k) dimension
  ! ndp1          Integer        Input        nd + 1
  ! ntra          Integer        Input        number of tracors
  ! iflag_con     Integer        Input        version of convect (3/4)
  ! t1            Real           Input        temperature
  ! q1            Real           Input        specific hum
  ! qs1           Real           Input        sat specific hum
  ! u1            Real           Input        u-wind
  ! v1            Real           Input        v-wind
  ! tra1          Real           Input        tracors
  ! p1            Real           Input        full level pressure
  ! ph1           Real           Input        half level pressure
  ! iflag1        Integer        Output       flag for Emanuel conditions
  ! ft1           Real           Output       temp tend
  ! fq1           Real           Output       spec hum tend
  ! fu1           Real           Output       u-wind tend
  ! fv1           Real           Output       v-wind tend
  ! ftra1         Real           Output       tracor tend
  ! precip1       Real           Output       precipitation
  ! VPrecip1      Real           Output       vertical profile of
  ! precipitations
  ! cbmf1         Real           Output       cloud base mass flux
  ! sig1          Real           In/Out       section adiabatic updraft
  ! w01           Real           In/Out       vertical velocity within adiab
  ! updraft
  ! delt          Real           Input        time step
  ! Ma1           Real           Output       mass flux adiabatic updraft
  ! upwd1         Real           Output       total upward mass flux
  ! (adiab+mixed)
  ! dnwd1         Real           Output       saturated downward mass flux
  ! (mixed)
  ! dnwd01        Real           Output       unsaturated downward mass flux
  ! qcondc1       Real           Output       in-cld mixing ratio of
  ! condensed water
  ! wd1           Real           Output       downdraft velocity scale for
  ! sfc fluxes
  ! cape1         Real           Output       CAPE

  ! wdtrainA1     Real           Output   precipitation detrained from
  ! adiabatic draught;
  ! used in tracer transport (cvltr)
  ! wdtrainM1     Real           Output   precipitation detrained from mixed
  ! draughts;
  ! used in tracer transport (cvltr)
  ! da1           Real           Output   used in tracer transport (cvltr)
  ! phi1          Real           Output   used in tracer transport (cvltr)
  ! mp1           Real           Output   used in tracer transport (cvltr)

  ! phi21         Real           Output   used in tracer transport (cvltr)

  ! d1a1          Real           Output   used in tracer transport (cvltr)
  ! dam1          Real           Output   used in tracer transport (cvltr)

  ! evap1         Real           Output
  ! ep1           Real           Output
  ! sij1        Real           Output
  ! elij1         Real           Output

  ! S. Bony, Mar 2002:
  ! * Several modules corresponding to different physical processes
  ! * Several versions of convect may be used:
  ! - iflag_con=3: version lmd  (previously named convect3)
  ! - iflag_con=4: version 4.3b (vect. version, previously convect1/2)
  ! + tard: 	- iflag_con=5: version lmd with ice (previously named convectg)
  ! S. Bony, Oct 2002:
  ! * Vectorization of convect3 (ie version lmd)

  ! ..............................END PROLOGUE.............................


  ! Input
  INTEGER len
  INTEGER nd
  INTEGER ndp1
  INTEGER noff
  INTEGER iflag_con
  INTEGER ntra
  REAL delt
  REAL t1(len, nd)
  REAL q1(len, nd)
  REAL qs1(len, nd)
  REAL u1(len, nd)
  REAL v1(len, nd)
  REAL tra1(len, nd, ntra)
  REAL p1(len, nd)
  REAL ph1(len, ndp1)

  ! Output
  INTEGER iflag1(len)
  REAL ft1(len, nd)
  REAL fq1(len, nd)
  REAL fu1(len, nd)
  REAL fv1(len, nd)
  REAL ftra1(len, nd, ntra)
  REAL precip1(len)
  REAL cbmf1(len)
  REAL sig1(klon, klev)
  REAL w01(klon, klev)
  REAL vprecip1(len, nd+1)
  REAL evap1(len, nd) !RomP
  REAL ep1(len, nd) !RomP
  REAL ma1(len, nd)
  REAL upwd1(len, nd)
  REAL dnwd1(len, nd)
  REAL dnwd01(len, nd)

  REAL qcondc1(len, nd) ! cld
  REAL wd1(len) ! gust
  REAL cape1(len)

  ! RomP >>>
  REAL wdtraina1(len, nd), wdtrainm1(len, nd)
  REAL sij1(len, nd, nd), elij1(len, nd, nd)
  REAL da1(len, nd), phi1(len, nd, nd), mp1(len, nd)

  REAL phi21(len, nd, nd)
  REAL d1a1(len, nd), dam1(len, nd)
  REAL epmlmmm1(len, nd, nd), eplamm1(len, nd)
  ! RomP <<<
  REAL epmax_diag1 (len) ! epmax_cape     

  ! -------------------------------------------------------------------
  ! Original Prologue by Kerry Emanuel.
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

  ! VPrecip: Vertical profile of convective precipitation (kg/m2/s).

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


  INTEGER i, k, n, il, j
  INTEGER icbmax
  INTEGER nk1(klon)
  INTEGER icb1(klon)
  INTEGER inb1(klon)
  INTEGER icbs1(klon)

  REAL plcl1(klon)
  REAL tnk1(klon)
  REAL qnk1(klon)
  REAL gznk1(klon)
  REAL pnk1(klon)
  REAL qsnk1(klon)
  REAL pbase1(klon)
  REAL buoybase1(klon)

  REAL lv1(klon, klev)
  REAL cpn1(klon, klev)
  REAL tv1(klon, klev)
  REAL gz1(klon, klev)
  REAL hm1(klon, klev)
  REAL h1(klon, klev)
  REAL tp1(klon, klev)
  REAL tvp1(klon, klev)
  REAL clw1(klon, klev)
  REAL th1(klon, klev)

  INTEGER ncum

  ! (local) compressed fields:

  ! ym      integer nloc
  ! ym      parameter (nloc=klon) ! pour l'instant
#define nloc klon
  INTEGER idcum(nloc)
  INTEGER iflag(nloc), nk(nloc), icb(nloc)
  INTEGER nent(nloc, klev)
  INTEGER icbs(nloc)
  INTEGER inb(nloc), inbis(nloc)

  REAL cbmf(nloc), plcl(nloc), tnk(nloc), qnk(nloc), gznk(nloc)
  REAL t(nloc, klev), q(nloc, klev), qs(nloc, klev)
  REAL u(nloc, klev), v(nloc, klev)
  REAL gz(nloc, klev), h(nloc, klev), lv(nloc, klev), cpn(nloc, klev)
  REAL p(nloc, klev), ph(nloc, klev+1), tv(nloc, klev), tp(nloc, klev)
  REAL clw(nloc, klev)
  REAL dph(nloc, klev)
  REAL pbase(nloc), buoybase(nloc), th(nloc, klev)
  REAL tvp(nloc, klev)
  REAL sig(nloc, klev), w0(nloc, klev)
  REAL hp(nloc, klev), ep(nloc, klev), sigp(nloc, klev)
  REAL frac(nloc), buoy(nloc, klev)
  REAL cape(nloc)
  REAL m(nloc, klev), ment(nloc, klev, klev), qent(nloc, klev, klev)
  REAL uent(nloc, klev, klev), vent(nloc, klev, klev)
  REAL ments(nloc, klev, klev), qents(nloc, klev, klev)
  REAL sij(nloc, klev, klev), elij(nloc, klev, klev)
  REAL qp(nloc, klev), up(nloc, klev), vp(nloc, klev)
  REAL wt(nloc, klev), water(nloc, klev), evap(nloc, klev)
  REAL b(nloc, klev), ft(nloc, klev), fq(nloc, klev)
  REAL fu(nloc, klev), fv(nloc, klev)
  REAL upwd(nloc, klev), dnwd(nloc, klev), dnwd0(nloc, klev)
  REAL ma(nloc, klev), mike(nloc, klev), tls(nloc, klev)
  REAL tps(nloc, klev), qprime(nloc), tprime(nloc)
  REAL precip(nloc)
  REAL vprecip(nloc, klev+1)
  REAL tra(nloc, klev, ntra), trap(nloc, klev, ntra)
  REAL ftra(nloc, klev, ntra), traent(nloc, klev, klev, ntra)
  REAL qcondc(nloc, klev) ! cld
  REAL wd(nloc) ! gust

  ! RomP >>>
  REAL da(nloc, klev), phi(nloc, klev, klev), mp(nloc, klev)
  REAL epmlmmm(nloc, klev, klev), eplamm(nloc, klev)
  REAL phi2(nloc, klev, klev)
  REAL d1a(nloc, klev), dam(nloc, klev)
  REAL wdtraina(nloc, klev), wdtrainm(nloc, klev)
  REAL sigd(nloc)
  ! RomP <<<
  REAL epmax_diag(nloc) ! epmax_cape

  nent(:, :) = 0
  ! -------------------------------------------------------------------
  ! --- SET CONSTANTS AND PARAMETERS
  ! -------------------------------------------------------------------
  ! print *, '-> cv_driver'      !jyg
  ! -- set simulation flags:
  ! (common cvflag)

  CALL cv_flag(0)

  ! -- set thermodynamical constants:
  ! (common cvthermo)

  CALL cv_thermo(iflag_con)

  ! -- set convect parameters

  ! includes microphysical parameters and parameters that
  ! control the rate of approach to quasi-equilibrium)
  ! (common cvparam)


  IF (iflag_con==30) THEN
    CALL cv30_param(nd, delt)
  END IF

  IF (iflag_con==4) THEN
    CALL cv_param(nd)
  END IF

  ! ---------------------------------------------------------------------
  ! --- INITIALIZE OUTPUT ARRAYS AND PARAMETERS
  ! ---------------------------------------------------------------------

  inb(:) = 0.0
  inb1(:) = 0.0
  icb1(:) = 0.0

  ft1(:, :) = 0.0
  fq1(:, :) = 0.0
  fu1(:, :) = 0.0
  fv1(:, :) = 0.0
  tvp1(:, :) = 0.0
  tp1(:, :) = 0.0
  clw1(:, :) = 0.0
  ! ym
  clw(:, :) = 0.0
  gz1(:, :) = 0.
  vprecip1(:, :) = 0.
  ma1(:, :) = 0.0
  upwd1(:, :) = 0.0
  dnwd1(:, :) = 0.0
  dnwd01(:, :) = 0.0
  qcondc1(:, :) = 0.0

  ftra1(:, :, :) = 0.0

  elij1(:, :, :) = 0.0
  sij1(:, :, :) = 0.0

  precip1(:) = 0.0
  iflag1(:) = 0
  wd1(:) = 0.0
  cape1(:) = 0.0
  epmax_diag1(:) = 0.0 ! epmax_cape


  IF (iflag_con==30) THEN
    DO il = 1, len
      sig1(il, nd) = sig1(il, nd) + 1.
      sig1(il, nd) = amin1(sig1(il,nd), 12.1)
    END DO
  END IF

  ! RomP >>>
  wdtraina1(:, :) = 0.
  wdtrainm1(:, :) = 0.
  da1(:, :) = 0.
  phi1(:, :, :) = 0.
  epmlmmm1(:, :, :) = 0.
  eplamm1(:, :) = 0.
  mp1(:, :) = 0.
  evap1(:, :) = 0.
  ep1(:, :) = 0.
  sij1(:, :, :) = 0.
  elij1(:, :, :) = 0.
  phi21(:, :, :) = 0.
  d1a1(:, :) = 0.
  dam1(:, :) = 0.
  ! RomP <<<

  ! --------------------------------------------------------------------
  ! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
  ! --------------------------------------------------------------------

  IF (iflag_con==30) THEN

    ! print*,'Emanuel version 30 '
    CALL cv30_prelim(len, nd, ndp1, t1, q1, p1, ph1 & ! nd->na
      , lv1, cpn1, tv1, gz1, h1, hm1, th1)
  END IF

  IF (iflag_con==4) THEN
    CALL cv_prelim(len, nd, ndp1, t1, q1, p1, ph1, lv1, cpn1, tv1, gz1, h1, &
      hm1)
  END IF

  ! --------------------------------------------------------------------
  ! --- CONVECTIVE FEED
  ! --------------------------------------------------------------------

  IF (iflag_con==30) THEN
    CALL cv30_feed(len, nd, t1, q1, qs1, p1, ph1, hm1, gz1 & !
                                                             ! nd->na
      , nk1, icb1, icbmax, iflag1, tnk1, qnk1, gznk1, plcl1)
  END IF

  IF (iflag_con==4) THEN
    CALL cv_feed(len, nd, t1, q1, qs1, p1, hm1, gz1, nk1, icb1, icbmax, &
      iflag1, tnk1, qnk1, gznk1, plcl1)
  END IF

  ! --------------------------------------------------------------------
  ! --- UNDILUTE (ADIABATIC) UPDRAFT / 1st part
  ! (up through ICB for convect4, up through ICB+1 for convect3)
  ! Calculates the lifted parcel virtual temperature at nk, the
  ! actual temperature, and the adiabatic liquid water content.
  ! --------------------------------------------------------------------

  IF (iflag_con==30) THEN
    CALL cv30_undilute1(len, nd, t1, q1, qs1, gz1, plcl1, p1, nk1, icb1 & ! nd->na
      , tp1, tvp1, clw1, icbs1)
  END IF

  IF (iflag_con==4) THEN
    CALL cv_undilute1(len, nd, t1, q1, qs1, gz1, p1, nk1, icb1, icbmax, tp1, &
      tvp1, clw1)
  END IF

  ! -------------------------------------------------------------------
  ! --- TRIGGERING
  ! -------------------------------------------------------------------

  IF (iflag_con==30) THEN
    CALL cv30_trigger(len, nd, icb1, plcl1, p1, th1, tv1, tvp1 & !
                                                                 ! nd->na
      , pbase1, buoybase1, iflag1, sig1, w01)
  END IF

  IF (iflag_con==4) THEN
    CALL cv_trigger(len, nd, icb1, cbmf1, tv1, tvp1, iflag1)
  END IF

  ! =====================================================================
  ! --- IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY
  ! =====================================================================

  ncum = 0
  DO i = 1, len
    IF (iflag1(i)==0) THEN
      ncum = ncum + 1
      idcum(ncum) = i
    END IF
  END DO

  ! print*,'cv_driver : klon, ncum = ',len,ncum

  IF (ncum>0) THEN

    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! --- COMPRESS THE FIELDS
    ! (-> vectorization over convective gridpoints)
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    IF (iflag_con==30) THEN
      CALL cv30_compress(len, nloc, ncum, nd, ntra, iflag1, nk1, icb1, icbs1, &
        plcl1, tnk1, qnk1, gznk1, pbase1, buoybase1, t1, q1, qs1, u1, v1, &
        gz1, th1, tra1, h1, lv1, cpn1, p1, ph1, tv1, tp1, tvp1, clw1, sig1, &
        w01, iflag, nk, icb, icbs, plcl, tnk, qnk, gznk, pbase, buoybase, t, &
        q, qs, u, v, gz, th, tra, h, lv, cpn, p, ph, tv, tp, tvp, clw, sig, &
        w0)
    END IF

    IF (iflag_con==4) THEN
      CALL cv_compress(len, nloc, ncum, nd, iflag1, nk1, icb1, cbmf1, plcl1, &
        tnk1, qnk1, gznk1, t1, q1, qs1, u1, v1, gz1, h1, lv1, cpn1, p1, ph1, &
        tv1, tp1, tvp1, clw1, iflag, nk, icb, cbmf, plcl, tnk, qnk, gznk, t, &
        q, qs, u, v, gz, h, lv, cpn, p, ph, tv, tp, tvp, clw, dph)
    END IF

    ! -------------------------------------------------------------------
    ! --- UNDILUTE (ADIABATIC) UPDRAFT / second part :
    ! ---   FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
    ! ---   &
    ! ---   COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
    ! ---   FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
    ! ---   &
    ! ---   FIND THE LEVEL OF NEUTRAL BUOYANCY
    ! -------------------------------------------------------------------

    IF (iflag_con==30) THEN
      CALL cv30_undilute2(nloc, ncum, nd, icb, icbs, nk & !na->nd
        , tnk, qnk, gznk, t, q, qs, gz, p, h, tv, lv, pbase, buoybase, plcl, &
        inb, tp, tvp, clw, hp, ep, sigp, buoy)
    END IF

    IF (iflag_con==4) THEN
      CALL cv_undilute2(nloc, ncum, nd, icb, nk, tnk, qnk, gznk, t, q, qs, &
        gz, p, dph, h, tv, lv, inb, inbis, tp, tvp, clw, hp, ep, sigp, frac)
    END IF

    ! -------------------------------------------------------------------
    ! --- CLOSURE
    ! -------------------------------------------------------------------

    IF (iflag_con==30) THEN
      CALL cv30_closure(nloc, ncum, nd, icb, inb & ! na->nd
        , pbase, p, ph, tv, buoy, sig, w0, cape, m)

      ! epmax_cape
      call cv30_epmax_fn_cape(nloc,ncum,nd &
                ,cape,ep,hp,icb,inb,clw,nk,t,h,lv &
                ,epmax_diag)
        ! on Ã©crase ep et recalcule hp
    END IF

    IF (iflag_con==4) THEN
      CALL cv_closure(nloc, ncum, nd, nk, icb, tv, tvp, p, ph, dph, plcl, &
        cpn, iflag, cbmf)
    END IF
    

    ! -------------------------------------------------------------------
    ! --- MIXING
    ! -------------------------------------------------------------------

    IF (iflag_con==30) THEN
      CALL cv30_mixing(nloc, ncum, nd, nd, ntra, icb, nk, inb & !
                                                                ! na->nd
        , ph, t, q, qs, u, v, tra, h, lv, qnk, hp, tv, tvp, ep, clw, m, sig, &
        ment, qent, uent, vent, sij, elij, ments, qents, traent)
    END IF

    IF (iflag_con==4) THEN
      CALL cv_mixing(nloc, ncum, nd, icb, nk, inb, inbis, ph, t, q, qs, u, v, &
        h, lv, qnk, hp, tv, tvp, ep, clw, cbmf, m, ment, qent, uent, vent, &
        nent, sij, elij)
    END IF

    ! -------------------------------------------------------------------
    ! --- UNSATURATED (PRECIPITATING) DOWNDRAFTS
    ! -------------------------------------------------------------------

    IF (iflag_con==30) THEN
      ! RomP >>>
      CALL cv30_unsat(nloc, ncum, nd, nd, ntra, icb, inb & ! na->nd
        , t, q, qs, gz, u, v, tra, p, ph, th, tv, lv, cpn, ep, sigp, clw, m, &
        ment, elij, delt, plcl, mp, qp, up, vp, trap, wt, water, evap, b, &
        wdtraina, wdtrainm)
      ! RomP <<<
    END IF

    IF (iflag_con==4) THEN
      CALL cv_unsat(nloc, ncum, nd, inb, t, q, qs, gz, u, v, p, ph, h, lv, &
        ep, sigp, clw, m, ment, elij, iflag, mp, qp, up, vp, wt, water, evap)
    END IF

    ! -------------------------------------------------------------------
    ! --- YIELD
    ! (tendencies, precipitation, variables of interface with other
    ! processes, etc)
    ! -------------------------------------------------------------------

    IF (iflag_con==30) THEN
      CALL cv30_yield(nloc, ncum, nd, nd, ntra & ! na->nd
        , icb, inb, delt, t, q, u, v, tra, gz, p, ph, h, hp, lv, cpn, th, ep, &
        clw, m, tp, mp, qp, up, vp, trap, wt, water, evap, b, ment, qent, &
        uent, vent, nent, elij, traent, sig, tv, tvp, iflag, precip, vprecip, &
        ft, fq, fu, fv, ftra, upwd, dnwd, dnwd0, ma, mike, tls, tps, qcondc, &
        wd)
    END IF

    IF (iflag_con==4) THEN
      CALL cv_yield(nloc, ncum, nd, nk, icb, inb, delt, t, q, u, v, gz, p, &
        ph, h, hp, lv, cpn, ep, clw, frac, m, mp, qp, up, vp, wt, water, &
        evap, ment, qent, uent, vent, nent, elij, tv, tvp, iflag, wd, qprime, &
        tprime, precip, cbmf, ft, fq, fu, fv, ma, qcondc)
    END IF

    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! --- passive tracers
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    IF (iflag_con==30) THEN
      ! RomP >>>
      CALL cv30_tracer(nloc, len, ncum, nd, nd, ment, sij, da, phi, phi2, &
        d1a, dam, ep, vprecip, elij, clw, epmlmmm, eplamm, icb, inb)
      ! RomP <<<
    END IF

    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! --- UNCOMPRESS THE FIELDS
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! set iflag1 =42 for non convective points
    DO i = 1, len
      iflag1(i) = 42
    END DO

    IF (iflag_con==30) THEN
      CALL cv30_uncompress(nloc, len, ncum, nd, ntra, idcum, iflag, precip, &
        vprecip, evap, ep, sig, w0 & !RomP
        , ft, fq, fu, fv, ftra, inb, ma, upwd, dnwd, dnwd0, qcondc, wd, cape, &
        da, phi, mp, phi2, d1a, dam, sij & !RomP
        , elij, clw, epmlmmm, eplamm & !RomP
        , wdtraina, wdtrainm,epmax_diag &     !RomP
        , iflag1, precip1, vprecip1, evap1, ep1, sig1, w01 & !RomP
        , ft1, fq1, fu1, fv1, ftra1, inb1, ma1, upwd1, dnwd1, dnwd01, &
        qcondc1, wd1, cape1, da1, phi1, mp1, phi21, d1a1, dam1, sij1 & !RomP
        , elij1, clw1, epmlmmm1, eplamm1 & !RomP
        , wdtraina1, wdtrainm1,epmax_diag1) !RomP
    END IF

    IF (iflag_con==4) THEN
      CALL cv_uncompress(nloc, len, ncum, nd, idcum, iflag, precip, cbmf, ft, &
        fq, fu, fv, ma, qcondc, iflag1, precip1, cbmf1, ft1, fq1, fu1, fv1, &
        ma1, qcondc1)
    END IF

  END IF ! ncum>0

  ! print *, 'fin cv_driver ->'      !jyg
  RETURN
END SUBROUTINE cv_driver

! ==================================================================
SUBROUTINE cv_flag(iflag_ice_thermo)

  USE ioipsl_getin_p_mod, ONLY : getin_p

  IMPLICIT NONE

  ! Argument : iflag_ice_thermo : ice thermodynamics is taken into account if
  ! iflag_ice_thermo >=1
  INTEGER iflag_ice_thermo

  include "cvflag.h"

  ! -- si .TRUE., on rend la gravite plus explicite et eventuellement
  ! differente de 10.0 dans convect3:
  cvflag_grav = .TRUE.
  cvflag_ice = iflag_ice_thermo >= 1
  !
! si icvflag_Tpa=0, alors la fraction de glace dans l'ascendance adiabatique est
  ! fonction de la temperature de l'environnement et la temperature de l'ascendance est 
  ! calculee en deux itérations, une en supposant qu'il n'y a pas de glace et l'autre 
  ! en ajoutant la glace (ancien schéma d'Arnaud Jam).
! si icvflag_Tpa=1, alors la fraction de glace dans l'ascendance adiabatique est
  ! fonction de la temperature de l'environnement et la temperature de l'ascendance est 
  ! calculee en une seule iteration.
! si icvflag_Tpa=2, alors la fraction de glace dans l'ascendance adiabatique est
  ! fonction de la temperature de l'ascendance et la temperature de l'ascendance est 
  ! calculee en une seule iteration.
  icvflag_Tpa=0  
  call getin_p('icvflag_Tpa', icvflag_Tpa)

  RETURN
END SUBROUTINE cv_flag

! ==================================================================
SUBROUTINE cv_thermo(iflag_con)
  IMPLICIT NONE

  ! -------------------------------------------------------------
  ! Set thermodynamical constants for convectL
  ! -------------------------------------------------------------

  include "YOMCST.h"
  include "cvthermo.h"

  INTEGER iflag_con


  ! original set from convect:
  IF (iflag_con==4) THEN
    cpd = 1005.7
    cpv = 1870.0
    cl = 4190.0
    rrv = 461.5
    rrd = 287.04
    lv0 = 2.501E6
    g = 9.8
    t0 = 273.15
    grav = g
  ELSE

    ! constants consistent with LMDZ:
    cpd = rcpd
    cpv = rcpv
    cl = rcw
    ci = rcs
    rrv = rv
    rrd = rd
    lv0 = rlvtt
    lf0 = rlstt - rlvtt
    g = rg ! not used in convect3
    ! ori      t0  = RTT
    t0 = 273.15 ! convect3 (RTT=273.16)
    ! maf       grav= 10.    ! implicitely or explicitely used in convect3
    grav = g ! implicitely or explicitely used in convect3
  END IF

  rowl = 1000.0 !(a quelle variable de YOMCST cela correspond-il?)

  clmcpv = cl - cpv
  clmcpd = cl - cpd
  clmci = cl - ci
  cpdmcp = cpd - cpv
  cpvmcpd = cpv - cpd
  cpvmcl = cl - cpv ! for convect3
  eps = rrd/rrv
  epsi = 1.0/eps
  epsim1 = epsi - 1.0
  ! ginv=1.0/g
  ginv = 1.0/grav
  hrd = 0.5*rrd

  RETURN
END SUBROUTINE cv_thermo
