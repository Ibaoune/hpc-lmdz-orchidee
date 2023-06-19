
! $Header$

SUBROUTINE conccm(dtime, paprs, pplay, t, q, conv_q, d_t, d_q, rain, snow, &
    kbascm, ktopcm)

  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: le 14 mars 1996
  ! Objet: Schema simple (avec flux de masse) pour la convection
  ! (schema standard du modele NCAR CCM2)
  ! ======================================================================
  include "YOMCST.h"
  include "YOETHF.h"

  ! Entree:
  REAL dtime ! pas d'integration
  REAL paprs(klon, klev+1) ! pression inter-couche (Pa)
  REAL pplay(klon, klev) ! pression au milieu de couche (Pa)
  REAL t(klon, klev) ! temperature (K)
  REAL q(klon, klev) ! humidite specifique (g/g)
  REAL conv_q(klon, klev) ! taux de convergence humidite (g/g/s)
  ! Sortie:
  REAL d_t(klon, klev) ! incrementation temperature
  REAL d_q(klon, klev) ! incrementation vapeur
  REAL rain(klon) ! pluie (mm/s)
  REAL snow(klon) ! neige (mm/s)
  INTEGER kbascm(klon) ! niveau du bas de convection
  INTEGER ktopcm(klon) ! niveau du haut de convection

  REAL pt(klon, klev)
  REAL pq(klon, klev)
  REAL pres(klon, klev)
  REAL dp(klon, klev)
  REAL zgeom(klon, klev)
  REAL cmfprs(klon)
  REAL cmfprt(klon)
  INTEGER ntop(klon)
  INTEGER nbas(klon)
  INTEGER i, k
  REAL zlvdcp, zlsdcp, zdelta, zz, za, zb

  LOGICAL usekuo ! utiliser convection profonde (schema Kuo)
  PARAMETER (usekuo=.TRUE.)

  REAL d_t_bis(klon, klev)
  REAL d_q_bis(klon, klev)
  REAL rain_bis(klon)
  REAL snow_bis(klon)
  INTEGER ibas_bis(klon)
  INTEGER itop_bis(klon)
  REAL d_ql_bis(klon, klev)
  REAL rneb_bis(klon, klev)

  ! initialiser les variables de sortie (pour securite)
  DO i = 1, klon
    rain(i) = 0.0
    snow(i) = 0.0
    kbascm(i) = 0
    ktopcm(i) = 0
  END DO
  DO k = 1, klev
    DO i = 1, klon
      d_t(i, k) = 0.0
      d_q(i, k) = 0.0
    END DO
  END DO

  ! preparer les variables d'entree (attention: l'ordre des niveaux
  ! verticaux augmente du haut vers le bas)
  DO k = 1, klev
    DO i = 1, klon
      pt(i, k) = t(i, klev-k+1)
      pq(i, k) = q(i, klev-k+1)
      pres(i, k) = pplay(i, klev-k+1)
      dp(i, k) = paprs(i, klev+1-k) - paprs(i, klev+1-k+1)
    END DO
  END DO
  DO i = 1, klon
    zgeom(i, klev) = rd*t(i, 1)/(0.5*(paprs(i,1)+pplay(i, &
      1)))*(paprs(i,1)-pplay(i,1))
  END DO
  DO k = 2, klev
    DO i = 1, klon
      zgeom(i, klev+1-k) = zgeom(i, klev+1-k+1) + rd*0.5*(t(i,k-1)+t(i,k))/ &
        paprs(i, k)*(pplay(i,k-1)-pplay(i,k))
    END DO
  END DO

  CALL cmfmca(dtime, pres, dp, zgeom, pt, pq, cmfprt, cmfprs, ntop, nbas)

  DO k = 1, klev
    DO i = 1, klon
      d_q(i, klev+1-k) = pq(i, k) - q(i, klev+1-k)
      d_t(i, klev+1-k) = pt(i, k) - t(i, klev+1-k)
    END DO
  END DO

  DO i = 1, klon
    rain(i) = cmfprt(i)*rhoh2o
    snow(i) = cmfprs(i)*rhoh2o
    kbascm(i) = klev + 1 - nbas(i)
    ktopcm(i) = klev + 1 - ntop(i)
  END DO

  IF (usekuo) THEN
    CALL conkuo(dtime, paprs, pplay, t, q, conv_q, d_t_bis, d_q_bis, &
      d_ql_bis, rneb_bis, rain_bis, snow_bis, ibas_bis, itop_bis)
    DO k = 1, klev
      DO i = 1, klon
        d_t(i, k) = d_t(i, k) + d_t_bis(i, k)
        d_q(i, k) = d_q(i, k) + d_q_bis(i, k)
      END DO
    END DO
    DO i = 1, klon
      rain(i) = rain(i) + rain_bis(i)
      snow(i) = snow(i) + snow_bis(i)
      kbascm(i) = min(kbascm(i), ibas_bis(i))
      ktopcm(i) = max(ktopcm(i), itop_bis(i))
    END DO
    DO k = 1, klev ! eau liquide convective est
      DO i = 1, klon ! dispersee dans l'air
        zlvdcp = rlvtt/rcpd/(1.0+rvtmp2*q(i,k))
        zlsdcp = rlstt/rcpd/(1.0+rvtmp2*q(i,k))
        zdelta = max(0., sign(1.,rtt-t(i,k)))
        zz = d_ql_bis(i, k) ! re-evap. de l'eau liquide
        zb = max(0.0, zz)
        za = -max(0.0, zz)*(zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
        d_t(i, k) = d_t(i, k) + za
        d_q(i, k) = d_q(i, k) + zb
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE conccm
SUBROUTINE cmfmca(deltat, p, dp, gz, tb, shb, cmfprt, cmfprs, cnt, cnb)
  USE dimphy
  IMPLICIT NONE
  ! -----------------------------------------------------------------------
  ! Moist convective mass flux procedure:
  ! If stratification is unstable to nonentraining parcel ascent,
  ! complete an adjustment making use of a simple cloud model

  ! Code generalized to allow specification of parcel ("updraft")
  ! properties, as well as convective transport of an arbitrary
  ! number of passive constituents (see cmrb array).
  ! ----------------------------Code History-------------------------------
  ! Original version:  J. J. Hack, March 22, 1990
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          J. Hack, G. Taylor, August 1992
  ! Adaptation au LMD: Z.X. Li, mars 1996 (reference: Hack 1994,
  ! J. Geophys. Res. vol 99, D3, 5551-5568). J'ai
  ! introduit les constantes et les fonctions thermo-
  ! dynamiques du Centre Europeen. J'ai elimine le
  ! re-indicage du code en esperant que cela pourra
  ! simplifier la lecture et la comprehension.
  ! -----------------------------------------------------------------------
  INTEGER pcnst ! nombre de traceurs passifs
  PARAMETER (pcnst=1)
  ! ------------------------------Arguments--------------------------------
  ! Input arguments

  REAL deltat ! time step (seconds)
  REAL p(klon, klev) ! pressure
  REAL dp(klon, klev) ! delta-p
  REAL gz(klon, klev) ! geopotential (a partir du sol)

  REAL thtap(klon) ! PBL perturbation theta
  REAL shp(klon) ! PBL perturbation specific humidity
  REAL pblh(klon) ! PBL height (provided by PBL routine)
  REAL cmrp(klon, pcnst) ! constituent perturbations in PBL

  ! Updated arguments:

  REAL tb(klon, klev) ! temperature (t bar)
  REAL shb(klon, klev) ! specific humidity (sh bar)
  REAL cmrb(klon, klev, pcnst) ! constituent mixing ratios (cmr bar)

  ! Output arguments

  REAL cmfdt(klon, klev) ! dT/dt due to moist convection
  REAL cmfdq(klon, klev) ! dq/dt due to moist convection
  REAL cmfmc(klon, klev) ! moist convection cloud mass flux
  REAL cmfdqr(klon, klev) ! dq/dt due to convective rainout
  REAL cmfsl(klon, klev) ! convective lw static energy flux
  REAL cmflq(klon, klev) ! convective total water flux
  REAL cmfprt(klon) ! convective precipitation rate
  REAL cmfprs(klon) ! convective snowfall rate
  REAL qc(klon, klev) ! dq/dt due to rainout terms
  INTEGER cnt(klon) ! top level of convective activity
  INTEGER cnb(klon) ! bottom level of convective activity
  ! ------------------------------Parameters-------------------------------
  REAL c0 ! rain water autoconversion coefficient
  PARAMETER (c0=1.0E-4)
  REAL dzmin ! minimum convective depth for precipitation
  PARAMETER (dzmin=0.0)
  REAL betamn ! minimum overshoot parameter
  PARAMETER (betamn=0.10)
  REAL cmftau ! characteristic adjustment time scale
  PARAMETER (cmftau=3600.)
  INTEGER limcnv ! top interface level limit for convection
  PARAMETER (limcnv=1)
  REAL tpmax ! maximum acceptable t perturbation (degrees C)
  PARAMETER (tpmax=1.50)
  REAL shpmax ! maximum acceptable q perturbation (g/g)
  PARAMETER (shpmax=1.50E-3)
  REAL tiny ! arbitrary small num used in transport estimates
  PARAMETER (tiny=1.0E-36)
  REAL eps ! convergence criteria (machine dependent)
  PARAMETER (eps=1.0E-13)
  REAL tmelt ! freezing point of water(req'd for rain vs snow)
  PARAMETER (tmelt=273.15)
  REAL ssfac ! supersaturation bound (detrained air)
  PARAMETER (ssfac=1.001)

  ! ---------------------------Local workspace-----------------------------
  REAL gam(klon, klev) ! L/cp (d(qsat)/dT)
  REAL sb(klon, klev) ! dry static energy (s bar)
  REAL hb(klon, klev) ! moist static energy (h bar)
  REAL shbs(klon, klev) ! sat. specific humidity (sh bar star)
  REAL hbs(klon, klev) ! sat. moist static energy (h bar star)
  REAL shbh(klon, klev+1) ! specific humidity on interfaces
  REAL sbh(klon, klev+1) ! s bar on interfaces
  REAL hbh(klon, klev+1) ! h bar on interfaces
  REAL cmrh(klon, klev+1) ! interface constituent mixing ratio
  REAL prec(klon) ! instantaneous total precipitation
  REAL dzcld(klon) ! depth of convective layer (m)
  REAL beta(klon) ! overshoot parameter (fraction)
  REAL betamx ! local maximum on overshoot
  REAL eta(klon) ! convective mass flux (kg/m^2 s)
  REAL etagdt ! eta*grav*deltat
  REAL cldwtr(klon) ! cloud water (mass)
  REAL rnwtr(klon) ! rain water  (mass)
  REAL sc(klon) ! dry static energy   ("in-cloud")
  REAL shc(klon) ! specific humidity   ("in-cloud")
  REAL hc(klon) ! moist static energy ("in-cloud")
  REAL cmrc(klon) ! constituent mix rat ("in-cloud")
  REAL dq1(klon) ! shb  convective change (lower lvl)
  REAL dq2(klon) ! shb  convective change (mid level)
  REAL dq3(klon) ! shb  convective change (upper lvl)
  REAL ds1(klon) ! sb   convective change (lower lvl)
  REAL ds2(klon) ! sb   convective change (mid level)
  REAL ds3(klon) ! sb   convective change (upper lvl)
  REAL dcmr1(klon) ! cmrb convective change (lower lvl)
  REAL dcmr2(klon) ! cmrb convective change (mid level)
  REAL dcmr3(klon) ! cmrb convective change (upper lvl)
  REAL flotab(klon) ! hc - hbs (mesure d'instabilite)
  LOGICAL ldcum(klon) ! .true. si la convection existe
  LOGICAL etagt0 ! true if eta > 0.0
  REAL dt ! current 2 delta-t (model time step)
  REAL cats ! modified characteristic adj. time
  REAL rdt ! 1./dt
  REAL qprime ! modified specific humidity pert.
  REAL tprime ! modified thermal perturbation
  REAL pblhgt ! bounded pbl height (max[pblh,1m])
  REAL fac1 ! intermediate scratch variable
  REAL shprme ! intermediate specific humidity pert.
  REAL qsattp ! saturation mixing ratio for
  ! !  thermally perturbed PBL parcels
  REAL dz ! local layer depth
  REAL b1 ! bouyancy measure in detrainment lvl
  REAL b2 ! bouyancy measure in condensation lvl
  REAL g ! bounded vertical gradient of hb
  REAL tmass ! total mass available for convective exchange
  REAL denom ! intermediate scratch variable
  REAL qtest1 ! used in negative q test (middle lvl)
  REAL qtest2 ! used in negative q test (lower lvl)
  REAL fslkp ! flux lw static energy (bot interface)
  REAL fslkm ! flux lw static energy (top interface)
  REAL fqlkp ! flux total water (bottom interface)
  REAL fqlkm ! flux total water (top interface)
  REAL botflx ! bottom constituent mixing ratio flux
  REAL topflx ! top constituent mixing ratio flux
  REAL efac1 ! ratio cmrb to convectively induced change (bot lvl)
  REAL efac2 ! ratio cmrb to convectively induced change (mid lvl)
  REAL efac3 ! ratio cmrb to convectively induced change (top lvl)

  INTEGER i, k ! indices horizontal et vertical
  INTEGER km1 ! k-1 (index offset)
  INTEGER kp1 ! k+1 (index offset)
  INTEGER m ! constituent index
  INTEGER ktp ! temporary index used to track top
  INTEGER is ! nombre de points a ajuster

  REAL tmp1, tmp2, tmp3, tmp4
  REAL zx_t, zx_p, zx_q, zx_qs, zx_gam
  REAL zcor, zdelta, zcvm5

  REAL qhalf, sh1, sh2, shbs1, shbs2
  include "YOMCST.h"
  include "YOETHF.h"
  include "FCTTRE.h"
  qhalf(sh1, sh2, shbs1, shbs2) = min(max(sh1,sh2), &
    (shbs2*sh1+shbs1*sh2)/(shbs1+shbs2))

  ! -----------------------------------------------------------------------
  ! pas de traceur pour l'instant
  DO m = 1, pcnst
    DO k = 1, klev
      DO i = 1, klon
        cmrb(i, k, m) = 0.0
      END DO
    END DO
  END DO

  ! Les perturbations de la couche limite sont zero pour l'instant

  DO m = 1, pcnst
    DO i = 1, klon
      cmrp(i, m) = 0.0
    END DO
  END DO
  DO i = 1, klon
    thtap(i) = 0.0
    shp(i) = 0.0
    pblh(i) = 1.0
  END DO

  ! Ensure that characteristic adjustment time scale (cmftau) assumed
  ! in estimate of eta isn't smaller than model time scale (deltat)

  dt = deltat
  cats = max(dt, cmftau)
  rdt = 1.0/dt

  ! Compute sb,hb,shbs,hbs

  DO k = 1, klev
    DO i = 1, klon
      zx_t = tb(i, k)
      zx_p = p(i, k)
      zx_q = shb(i, k)
      zdelta = max(0., sign(1.,rtt-zx_t))
      zcvm5 = r5les*rlvtt*(1.-zdelta) + r5ies*rlstt*zdelta
      zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*zx_q)
      zx_qs = r2es*foeew(zx_t, zdelta)/zx_p
      zx_qs = min(0.5, zx_qs)
      zcor = 1./(1.-retv*zx_qs)
      zx_qs = zx_qs*zcor
      zx_gam = foede(zx_t, zdelta, zcvm5, zx_qs, zcor)
      shbs(i, k) = zx_qs
      gam(i, k) = zx_gam
    END DO
  END DO

  DO k = limcnv, klev
    DO i = 1, klon
      sb(i, k) = rcpd*tb(i, k) + gz(i, k)
      hb(i, k) = sb(i, k) + rlvtt*shb(i, k)
      hbs(i, k) = sb(i, k) + rlvtt*shbs(i, k)
    END DO
  END DO

  ! Compute sbh, shbh

  DO k = limcnv + 1, klev
    km1 = k - 1
    DO i = 1, klon
      sbh(i, k) = 0.5*(sb(i,km1)+sb(i,k))
      shbh(i, k) = qhalf(shb(i,km1), shb(i,k), shbs(i,km1), shbs(i,k))
      hbh(i, k) = sbh(i, k) + rlvtt*shbh(i, k)
    END DO
  END DO

  ! Specify properties at top of model (not used, but filling anyway)

  DO i = 1, klon
    sbh(i, limcnv) = sb(i, limcnv)
    shbh(i, limcnv) = shb(i, limcnv)
    hbh(i, limcnv) = hb(i, limcnv)
  END DO

  ! Zero vertically independent control, tendency & diagnostic arrays

  DO i = 1, klon
    prec(i) = 0.0
    dzcld(i) = 0.0
    cnb(i) = 0
    cnt(i) = klev + 1
  END DO

  DO k = 1, klev
    DO i = 1, klon
      cmfdt(i, k) = 0.
      cmfdq(i, k) = 0.
      cmfdqr(i, k) = 0.
      cmfmc(i, k) = 0.
      cmfsl(i, k) = 0.
      cmflq(i, k) = 0.
    END DO
  END DO

  ! Begin moist convective mass flux adjustment procedure.
  ! Formalism ensures that negative cloud liquid water can never occur

  DO k = klev - 1, limcnv + 1, -1
    km1 = k - 1
    kp1 = k + 1
    DO i = 1, klon
      eta(i) = 0.0
      beta(i) = 0.0
      ds1(i) = 0.0
      ds2(i) = 0.0
      ds3(i) = 0.0
      dq1(i) = 0.0
      dq2(i) = 0.0
      dq3(i) = 0.0

      ! Specification of "cloud base" conditions

      qprime = 0.0
      tprime = 0.0

      ! Assign tprime within the PBL to be proportional to the quantity
      ! thtap (which will be bounded by tpmax), passed to this routine by
      ! the PBL routine.  Don't allow perturbation to produce a dry
      ! adiabatically unstable parcel.  Assign qprime within the PBL to be
      ! an appropriately modified value of the quantity shp (which will be
      ! bounded by shpmax) passed to this routine by the PBL routine.  The
      ! quantity qprime should be less than the local saturation value
      ! (qsattp=qsat[t+tprime,p]).  In both cases, thtap and shp are
      ! linearly reduced toward zero as the PBL top is approached.

      pblhgt = max(pblh(i), 1.0)
      IF (gz(i,kp1)/rg<=pblhgt .AND. dzcld(i)==0.0) THEN
        fac1 = max(0.0, 1.0-gz(i,kp1)/rg/pblhgt)
        tprime = min(thtap(i), tpmax)*fac1
        qsattp = shbs(i, kp1) + rcpd/rlvtt*gam(i, kp1)*tprime
        shprme = min(min(shp(i),shpmax)*fac1, max(qsattp-shb(i,kp1),0.0))
        qprime = max(qprime, shprme)
      ELSE
        tprime = 0.0
        qprime = 0.0
      END IF

      ! Specify "updraft" (in-cloud) thermodynamic properties

      sc(i) = sb(i, kp1) + rcpd*tprime
      shc(i) = shb(i, kp1) + qprime
      hc(i) = sc(i) + rlvtt*shc(i)
      flotab(i) = hc(i) - hbs(i, k)
      dz = dp(i, k)*rd*tb(i, k)/rg/p(i, k)
      IF (flotab(i)>0.0) THEN
        dzcld(i) = dzcld(i) + dz
      ELSE
        dzcld(i) = 0.0
      END IF
    END DO

    ! Check on moist convective instability

    is = 0
    DO i = 1, klon
      IF (flotab(i)>0.0) THEN
        ldcum(i) = .TRUE.
        is = is + 1
      ELSE
        ldcum(i) = .FALSE.
      END IF
    END DO

    IF (is==0) THEN
      DO i = 1, klon
        dzcld(i) = 0.0
      END DO
      GO TO 70
    END IF

    ! Current level just below top level => no overshoot

    IF (k<=limcnv+1) THEN
      DO i = 1, klon
        IF (ldcum(i)) THEN
          cldwtr(i) = sb(i, k) - sc(i) + flotab(i)/(1.0+gam(i,k))
          cldwtr(i) = max(0.0, cldwtr(i))
          beta(i) = 0.0
        END IF
      END DO
      GO TO 20
    END IF

    ! First guess at overshoot parameter using crude buoyancy closure
    ! 10% overshoot assumed as a minimum and 1-c0*dz maximum to start
    ! If pre-existing supersaturation in detrainment layer, beta=0
    ! cldwtr is temporarily equal to RLVTT*l (l=> liquid water)

    DO i = 1, klon
      IF (ldcum(i)) THEN
        cldwtr(i) = sb(i, k) - sc(i) + flotab(i)/(1.0+gam(i,k))
        cldwtr(i) = max(0.0, cldwtr(i))
        betamx = 1.0 - c0*max(0.0, (dzcld(i)-dzmin))
        b1 = (hc(i)-hbs(i,km1))*dp(i, km1)
        b2 = (hc(i)-hbs(i,k))*dp(i, k)
        beta(i) = max(betamn, min(betamx,1.0+b1/b2))
        IF (hbs(i,km1)<=hb(i,km1)) beta(i) = 0.0
      END IF
    END DO

    ! Bound maximum beta to ensure physically realistic solutions

    ! First check constrains beta so that eta remains positive
    ! (assuming that eta is already positive for beta equal zero)
    ! La premiere contrainte de beta est que le flux eta doit etre positif.

    DO i = 1, klon
      IF (ldcum(i)) THEN
        tmp1 = (1.0+gam(i,k))*(sc(i)-sbh(i,kp1)+cldwtr(i)) - &
          (hbh(i,kp1)-hc(i))*dp(i, k)/dp(i, kp1)
        tmp2 = (1.0+gam(i,k))*(sc(i)-sbh(i,k))
        IF ((beta(i)*tmp2-tmp1)>0.0) THEN
          betamx = 0.99*(tmp1/tmp2)
          beta(i) = max(0.0, min(betamx,beta(i)))
        END IF

        ! Second check involves supersaturation of "detrainment layer"
        ! small amount of supersaturation acceptable (by ssfac factor)
        ! La 2e contrainte est que la convection ne doit pas sursaturer
        ! la "detrainment layer", Neanmoins, une petite sursaturation
        ! est acceptee (facteur ssfac).

        IF (hb(i,km1)<hbs(i,km1)) THEN
          tmp1 = (1.0+gam(i,k))*(sc(i)-sbh(i,kp1)+cldwtr(i)) - &
            (hbh(i,kp1)-hc(i))*dp(i, k)/dp(i, kp1)
          tmp1 = tmp1/dp(i, k)
          tmp2 = gam(i, km1)*(sbh(i,k)-sc(i)+cldwtr(i)) - hbh(i, k) + hc(i) - &
            sc(i) + sbh(i, k)
          tmp3 = (1.0+gam(i,k))*(sc(i)-sbh(i,k))/dp(i, k)
          tmp4 = (dt/cats)*(hc(i)-hbs(i,k))*tmp2/(dp(i,km1)*(hbs(i,km1)-hb(i, &
            km1))) + tmp3
          IF ((beta(i)*tmp4-tmp1)>0.0) THEN
            betamx = ssfac*(tmp1/tmp4)
            beta(i) = max(0.0, min(betamx,beta(i)))
          END IF
        ELSE
          beta(i) = 0.0
        END IF

        ! Third check to avoid introducing 2 delta x thermodynamic
        ! noise in the vertical ... constrain adjusted h (or theta e)
        ! so that the adjustment doesn't contribute to "kinks" in h

        g = min(0.0, hb(i,k)-hb(i,km1))
        tmp3 = (hb(i,k)-hb(i,km1)-g)*(cats/dt)/(hc(i)-hbs(i,k))
        tmp1 = (1.0+gam(i,k))*(sc(i)-sbh(i,kp1)+cldwtr(i)) - &
          (hbh(i,kp1)-hc(i))*dp(i, k)/dp(i, kp1)
        tmp1 = tmp1/dp(i, k)
        tmp1 = tmp3*tmp1 + (hc(i)-hbh(i,kp1))/dp(i, k)
        tmp2 = tmp3*(1.0+gam(i,k))*(sc(i)-sbh(i,k))/dp(i, k) + &
          (hc(i)-hbh(i,k)-cldwtr(i))*(1.0/dp(i,k)+1.0/dp(i,kp1))
        IF ((beta(i)*tmp2-tmp1)>0.0) THEN
          betamx = 0.0
          IF (tmp2/=0.0) betamx = tmp1/tmp2
          beta(i) = max(0.0, min(betamx,beta(i)))
        END IF
      END IF
    END DO

    ! Calculate mass flux required for stabilization.

    ! Ensure that the convective mass flux, eta, is positive by
    ! setting negative values of eta to zero..
    ! Ensure that estimated mass flux cannot move more than the
    ! minimum of total mass contained in either layer k or layer k+1.
    ! Also test for other pathological cases that result in non-
    ! physical states and adjust eta accordingly.

20  CONTINUE
    DO i = 1, klon
      IF (ldcum(i)) THEN
        beta(i) = max(0.0, beta(i))
        tmp1 = hc(i) - hbs(i, k)
        tmp2 = ((1.0+gam(i,k))*(sc(i)-sbh(i,kp1)+cldwtr(i))-beta(i)*(1.0+gam( &
          i,k))*(sc(i)-sbh(i,k)))/dp(i, k) - (hbh(i,kp1)-hc(i))/dp(i, kp1)
        eta(i) = tmp1/(tmp2*rg*cats)
        tmass = min(dp(i,k), dp(i,kp1))/rg
        IF (eta(i)>tmass*rdt .OR. eta(i)<=0.0) eta(i) = 0.0

        ! Check on negative q in top layer (bound beta)

        IF (shc(i)-shbh(i,k)<0.0 .AND. beta(i)*eta(i)/=0.0) THEN
          denom = eta(i)*rg*dt*(shc(i)-shbh(i,k))/dp(i, km1)
          beta(i) = max(0.0, min(-0.999*shb(i,km1)/denom,beta(i)))
        END IF

        ! Check on negative q in middle layer (zero eta)

        qtest1 = shb(i, k) + eta(i)*rg*dt*((shc(i)-shbh(i, &
          kp1))-(1.0-beta(i))*cldwtr(i)/rlvtt-beta(i)*(shc(i)-shbh(i, &
          k)))/dp(i, k)
        IF (qtest1<=0.0) eta(i) = 0.0

        ! Check on negative q in lower layer (bound eta)

        fac1 = -(shbh(i,kp1)-shc(i))/dp(i, kp1)
        qtest2 = shb(i, kp1) - eta(i)*rg*dt*fac1
        IF (qtest2<0.0) THEN
          eta(i) = 0.99*shb(i, kp1)/(rg*dt*fac1)
        END IF
      END IF
    END DO


    ! Calculate cloud water, rain water, and thermodynamic changes

    DO i = 1, klon
      IF (ldcum(i)) THEN
        etagdt = eta(i)*rg*dt
        cldwtr(i) = etagdt*cldwtr(i)/rlvtt/rg
        rnwtr(i) = (1.0-beta(i))*cldwtr(i)
        ds1(i) = etagdt*(sbh(i,kp1)-sc(i))/dp(i, kp1)
        dq1(i) = etagdt*(shbh(i,kp1)-shc(i))/dp(i, kp1)
        ds2(i) = (etagdt*(sc(i)-sbh(i,kp1))+rlvtt*rg*cldwtr(i)-beta(i)*etagdt &
          *(sc(i)-sbh(i,k)))/dp(i, k)
        dq2(i) = (etagdt*(shc(i)-shbh(i,kp1))-rg*rnwtr(i)-beta(i)*etagdt*(shc &
          (i)-shbh(i,k)))/dp(i, k)
        ds3(i) = beta(i)*(etagdt*(sc(i)-sbh(i,k))-rlvtt*rg*cldwtr(i))/dp(i, &
          km1)
        dq3(i) = beta(i)*etagdt*(shc(i)-shbh(i,k))/dp(i, km1)

        ! Isolate convective fluxes for later diagnostics

        fslkp = eta(i)*(sc(i)-sbh(i,kp1))
        fslkm = beta(i)*(eta(i)*(sc(i)-sbh(i,k))-rlvtt*cldwtr(i)*rdt)
        fqlkp = eta(i)*(shc(i)-shbh(i,kp1))
        fqlkm = beta(i)*eta(i)*(shc(i)-shbh(i,k))


        ! Update thermodynamic profile (update sb, hb, & hbs later)

        tb(i, kp1) = tb(i, kp1) + ds1(i)/rcpd
        tb(i, k) = tb(i, k) + ds2(i)/rcpd
        tb(i, km1) = tb(i, km1) + ds3(i)/rcpd
        shb(i, kp1) = shb(i, kp1) + dq1(i)
        shb(i, k) = shb(i, k) + dq2(i)
        shb(i, km1) = shb(i, km1) + dq3(i)
        prec(i) = prec(i) + rnwtr(i)/rhoh2o

        ! Update diagnostic information for final budget
        ! Tracking temperature & specific humidity tendencies,
        ! rainout term, convective mass flux, convective liquid
        ! water static energy flux, and convective total water flux

        cmfdt(i, kp1) = cmfdt(i, kp1) + ds1(i)/rcpd*rdt
        cmfdt(i, k) = cmfdt(i, k) + ds2(i)/rcpd*rdt
        cmfdt(i, km1) = cmfdt(i, km1) + ds3(i)/rcpd*rdt
        cmfdq(i, kp1) = cmfdq(i, kp1) + dq1(i)*rdt
        cmfdq(i, k) = cmfdq(i, k) + dq2(i)*rdt
        cmfdq(i, km1) = cmfdq(i, km1) + dq3(i)*rdt
        cmfdqr(i, k) = cmfdqr(i, k) + (rg*rnwtr(i)/dp(i,k))*rdt
        cmfmc(i, kp1) = cmfmc(i, kp1) + eta(i)
        cmfmc(i, k) = cmfmc(i, k) + beta(i)*eta(i)
        cmfsl(i, kp1) = cmfsl(i, kp1) + fslkp
        cmfsl(i, k) = cmfsl(i, k) + fslkm
        cmflq(i, kp1) = cmflq(i, kp1) + rlvtt*fqlkp
        cmflq(i, k) = cmflq(i, k) + rlvtt*fqlkm
        qc(i, k) = (rg*rnwtr(i)/dp(i,k))*rdt
      END IF
    END DO

    ! Next, convectively modify passive constituents

    DO m = 1, pcnst
      DO i = 1, klon
        IF (ldcum(i)) THEN

          ! If any of the reported values of the constituent is negative in
          ! the three adjacent levels, nothing will be done to the profile

          IF ((cmrb(i,kp1,m)<0.0) .OR. (cmrb(i,k,m)<0.0) .OR. (cmrb(i,km1, &
            m)<0.0)) GO TO 40

          ! Specify constituent interface values (linear interpolation)

          cmrh(i, k) = 0.5*(cmrb(i,km1,m)+cmrb(i,k,m))
          cmrh(i, kp1) = 0.5*(cmrb(i,k,m)+cmrb(i,kp1,m))

          ! Specify perturbation properties of constituents in PBL

          pblhgt = max(pblh(i), 1.0)
          IF (gz(i,kp1)/rg<=pblhgt .AND. dzcld(i)==0.) THEN
            fac1 = max(0.0, 1.0-gz(i,kp1)/rg/pblhgt)
            cmrc(i) = cmrb(i, kp1, m) + cmrp(i, m)*fac1
          ELSE
            cmrc(i) = cmrb(i, kp1, m)
          END IF

          ! Determine fluxes, flux divergence => changes due to convection
          ! Logic must be included to avoid producing negative values. A bit
          ! messy since there are no a priori assumptions about profiles.
          ! Tendency is modified (reduced) when pending disaster detected.

          etagdt = eta(i)*rg*dt
          botflx = etagdt*(cmrc(i)-cmrh(i,kp1))
          topflx = beta(i)*etagdt*(cmrc(i)-cmrh(i,k))
          dcmr1(i) = -botflx/dp(i, kp1)
          efac1 = 1.0
          efac2 = 1.0
          efac3 = 1.0

          IF (cmrb(i,kp1,m)+dcmr1(i)<0.0) THEN
            efac1 = max(tiny, abs(cmrb(i,kp1,m)/dcmr1(i))-eps)
          END IF

          IF (efac1==tiny .OR. efac1>1.0) efac1 = 0.0
          dcmr1(i) = -efac1*botflx/dp(i, kp1)
          dcmr2(i) = (efac1*botflx-topflx)/dp(i, k)

          IF (cmrb(i,k,m)+dcmr2(i)<0.0) THEN
            efac2 = max(tiny, abs(cmrb(i,k,m)/dcmr2(i))-eps)
          END IF

          IF (efac2==tiny .OR. efac2>1.0) efac2 = 0.0
          dcmr2(i) = (efac1*botflx-efac2*topflx)/dp(i, k)
          dcmr3(i) = efac2*topflx/dp(i, km1)

          IF (cmrb(i,km1,m)+dcmr3(i)<0.0) THEN
            efac3 = max(tiny, abs(cmrb(i,km1,m)/dcmr3(i))-eps)
          END IF

          IF (efac3==tiny .OR. efac3>1.0) efac3 = 0.0
          efac3 = min(efac2, efac3)
          dcmr2(i) = (efac1*botflx-efac3*topflx)/dp(i, k)
          dcmr3(i) = efac3*topflx/dp(i, km1)

          cmrb(i, kp1, m) = cmrb(i, kp1, m) + dcmr1(i)
          cmrb(i, k, m) = cmrb(i, k, m) + dcmr2(i)
          cmrb(i, km1, m) = cmrb(i, km1, m) + dcmr3(i)
        END IF
40    END DO
    END DO ! end of m=1,pcnst loop

    IF (k==limcnv+1) GO TO 60 ! on ne pourra plus glisser

    ! Dans la procedure de glissage ascendant, les variables thermo-
    ! dynamiques des couches k et km1 servent au calcul des couches
    ! superieures. Elles ont donc besoin d'une mise-a-jour.

    DO i = 1, klon
      IF (ldcum(i)) THEN
        zx_t = tb(i, k)
        zx_p = p(i, k)
        zx_q = shb(i, k)
        zdelta = max(0., sign(1.,rtt-zx_t))
        zcvm5 = r5les*rlvtt*(1.-zdelta) + r5ies*rlstt*zdelta
        zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*zx_q)
        zx_qs = r2es*foeew(zx_t, zdelta)/zx_p
        zx_qs = min(0.5, zx_qs)
        zcor = 1./(1.-retv*zx_qs)
        zx_qs = zx_qs*zcor
        zx_gam = foede(zx_t, zdelta, zcvm5, zx_qs, zcor)
        shbs(i, k) = zx_qs
        gam(i, k) = zx_gam

        zx_t = tb(i, km1)
        zx_p = p(i, km1)
        zx_q = shb(i, km1)
        zdelta = max(0., sign(1.,rtt-zx_t))
        zcvm5 = r5les*rlvtt*(1.-zdelta) + r5ies*rlstt*zdelta
        zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*zx_q)
        zx_qs = r2es*foeew(zx_t, zdelta)/zx_p
        zx_qs = min(0.5, zx_qs)
        zcor = 1./(1.-retv*zx_qs)
        zx_qs = zx_qs*zcor
        zx_gam = foede(zx_t, zdelta, zcvm5, zx_qs, zcor)
        shbs(i, km1) = zx_qs
        gam(i, km1) = zx_gam

        sb(i, k) = sb(i, k) + ds2(i)
        sb(i, km1) = sb(i, km1) + ds3(i)
        hb(i, k) = sb(i, k) + rlvtt*shb(i, k)
        hb(i, km1) = sb(i, km1) + rlvtt*shb(i, km1)
        hbs(i, k) = sb(i, k) + rlvtt*shbs(i, k)
        hbs(i, km1) = sb(i, km1) + rlvtt*shbs(i, km1)

        sbh(i, k) = 0.5*(sb(i,k)+sb(i,km1))
        shbh(i, k) = qhalf(shb(i,km1), shb(i,k), shbs(i,km1), shbs(i,k))
        hbh(i, k) = sbh(i, k) + rlvtt*shbh(i, k)
        sbh(i, km1) = 0.5*(sb(i,km1)+sb(i,k-2))
        shbh(i, km1) = qhalf(shb(i,k-2), shb(i,km1), shbs(i,k-2), &
          shbs(i,km1))
        hbh(i, km1) = sbh(i, km1) + rlvtt*shbh(i, km1)
      END IF
    END DO

    ! Ensure that dzcld is reset if convective mass flux zero
    ! specify the current vertical extent of the convective activity
    ! top of convective layer determined by size of overshoot param.

60  CONTINUE
    DO i = 1, klon
      etagt0 = eta(i) > 0.0
      IF (.NOT. etagt0) dzcld(i) = 0.0
      IF (etagt0 .AND. beta(i)>betamn) THEN
        ktp = km1
      ELSE
        ktp = k
      END IF
      IF (etagt0) THEN
        cnt(i) = min(cnt(i), ktp)
        cnb(i) = max(cnb(i), k)
      END IF
    END DO
70 END DO ! end of k loop

  ! determine whether precipitation, prec, is frozen (snow) or not

  DO i = 1, klon
    IF (tb(i,klev)<tmelt .AND. tb(i,klev-1)<tmelt) THEN
      cmfprs(i) = prec(i)*rdt
    ELSE
      cmfprt(i) = prec(i)*rdt
    END IF
  END DO

  RETURN ! we're all done ... return to calling procedure
END SUBROUTINE cmfmca
