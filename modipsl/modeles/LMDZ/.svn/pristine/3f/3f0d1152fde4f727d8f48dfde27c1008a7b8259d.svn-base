SUBROUTINE drag_noro_strato(partdrag, nlon, nlev, dtime, paprs, pplay, pmea, pstd, &
    psig, pgam, pthe, ppic, pval, kgwd, kdx, ktest, t, u, v, pulow, pvlow, &
    pustr, pvstr, d_t, d_u, d_v)

  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): F.Lott (LMD/CNRS) date: 19950201
  ! Object: Mountain drag interface. Made necessary because:
  ! 1. in the LMD-GCM Layers are from bottom to top,
  ! contrary to most European GCM.
  ! 2. the altitude above ground of each model layers
  ! needs to be known (variable zgeom)
  ! ======================================================================
  ! Explicit Arguments:
  ! ==================
  ! partdrag-input-I-control which part of the drag we consider (total part or GW part)
  ! nlon----input-I-Total number of horizontal points that get into physics
  ! nlev----input-I-Number of vertical levels
  ! dtime---input-R-Time-step (s)
  ! paprs---input-R-Pressure in semi layers    (Pa)
  ! pplay---input-R-Pressure model-layers      (Pa)
  ! t-------input-R-temperature (K)
  ! u-------input-R-Horizontal wind (m/s)
  ! v-------input-R-Meridional wind (m/s)
  ! pmea----input-R-Mean Orography (m)
  ! pstd----input-R-SSO standard deviation (m)
  ! psig----input-R-SSO slope
  ! pgam----input-R-SSO Anisotropy
  ! pthe----input-R-SSO Angle
  ! ppic----input-R-SSO Peacks elevation (m)
  ! pval----input-R-SSO Valleys elevation (m)

  ! kgwd- -input-I: Total nb of points where the orography schemes are active
  ! ktest--input-I: Flags to indicate active points
  ! kdx----input-I: Locate the physical location of an active point.

  ! pulow, pvlow -output-R: Low-level wind
  ! pustr, pvstr -output-R: Surface stress due to SSO drag      (Pa)

  ! d_t-----output-R: T increment
  ! d_u-----output-R: U increment
  ! d_v-----output-R: V increment

  ! Implicit Arguments:
  ! ===================

  ! iim--common-I: Number of longitude intervals
  ! jjm--common-I: Number of latitude intervals
  ! klon-common-I: Number of points seen by the physics
  ! (iim+1)*(jjm+1) for instance
  ! klev-common-I: Number of vertical layers
  ! ======================================================================
  ! Local Variables:
  ! ================

  ! zgeom-----R: Altitude of layer above ground
  ! pt, pu, pv --R: t u v from top to bottom
  ! pdtdt, pdudt, pdvdt --R: t u v tendencies (from top to bottom)
  ! papmf: pressure at model layer (from top to bottom)
  ! papmh: pressure at model 1/2 layer (from top to bottom)

  ! ======================================================================
  include "YOMCST.h"
  include "YOEGWD.h"

  ! ARGUMENTS

  INTEGER partdrag,nlon, nlev
  REAL dtime
  REAL paprs(nlon, nlev+1)
  REAL pplay(nlon, nlev)
  REAL pmea(nlon), pstd(nlon), psig(nlon), pgam(nlon), pthe(nlon)
  REAL ppic(nlon), pval(nlon)
  REAL pulow(nlon), pvlow(nlon), pustr(nlon), pvstr(nlon)
  REAL t(nlon, nlev), u(nlon, nlev), v(nlon, nlev)
  REAL d_t(nlon, nlev), d_u(nlon, nlev), d_v(nlon, nlev)

  INTEGER i, k, kgwd, kdx(nlon), ktest(nlon)

  ! LOCAL VARIABLES:

  REAL zgeom(klon, klev)
  REAL pdtdt(klon, klev), pdudt(klon, klev), pdvdt(klon, klev)
  REAL pt(klon, klev), pu(klon, klev), pv(klon, klev)
  REAL papmf(klon, klev), papmh(klon, klev+1)
  CHARACTER (LEN=20) :: modname = 'orografi_strato'
  CHARACTER (LEN=80) :: abort_message

  ! INITIALIZE OUTPUT VARIABLES

  DO i = 1, klon
    pulow(i) = 0.0
    pvlow(i) = 0.0
    pustr(i) = 0.0
    pvstr(i) = 0.0
  END DO
  DO k = 1, klev
    DO i = 1, klon
      d_t(i, k) = 0.0
      d_u(i, k) = 0.0
      d_v(i, k) = 0.0
      pdudt(i, k) = 0.0
      pdvdt(i, k) = 0.0
      pdtdt(i, k) = 0.0
    END DO
  END DO

  ! PREPARE INPUT VARIABLES FOR ORODRAG (i.e., ORDERED FROM TOP TO BOTTOM)
  ! CALCULATE LAYERS HEIGHT ABOVE GROUND)

  DO k = 1, klev
    DO i = 1, klon
      pt(i, k) = t(i, klev-k+1)
      pu(i, k) = u(i, klev-k+1)
      pv(i, k) = v(i, klev-k+1)
      papmf(i, k) = pplay(i, klev-k+1)
    END DO
  END DO
  DO k = 1, klev + 1
    DO i = 1, klon
      papmh(i, k) = paprs(i, klev-k+2)
    END DO
  END DO
  DO i = 1, klon
    zgeom(i, klev) = rd*pt(i, klev)*log(papmh(i,klev+1)/papmf(i,klev))
  END DO
  DO k = klev - 1, 1, -1
    DO i = 1, klon
      zgeom(i, k) = zgeom(i, k+1) + rd*(pt(i,k)+pt(i,k+1))/2.0*log(papmf(i,k+ &
        1)/papmf(i,k))
    END DO
  END DO

  ! CALL SSO DRAG ROUTINES

  CALL orodrag_strato(partdrag,klon, klev, kgwd, kdx, ktest, dtime, papmh, papmf, &
    zgeom, pt, pu, pv, pmea, pstd, psig, pgam, pthe, ppic, pval, pulow, &
    pvlow, pdudt, pdvdt, pdtdt)

  ! COMPUTE INCREMENTS AND STRESS FROM TENDENCIES

  DO k = 1, klev
    DO i = 1, klon
      d_u(i, klev+1-k) = dtime*pdudt(i, k)
      d_v(i, klev+1-k) = dtime*pdvdt(i, k)
      d_t(i, klev+1-k) = dtime*pdtdt(i, k)
      pustr(i) = pustr(i) + pdudt(i, k)*(papmh(i,k+1)-papmh(i,k))/rg
      pvstr(i) = pvstr(i) + pdvdt(i, k)*(papmh(i,k+1)-papmh(i,k))/rg
    END DO
  END DO

  RETURN
END SUBROUTINE drag_noro_strato

SUBROUTINE orodrag_strato(partdrag,nlon, nlev, kgwd, kdx, ktest, ptsphy, paphm1, &
    papm1, pgeom1, ptm1, pum1, pvm1, pmea, pstd, psig, pgam, pthe, ppic, pval &
  ! outputs
    , pulow, pvlow, pvom, pvol, pte)

  USE dimphy
  IMPLICIT NONE


  ! **** *orodrag* - does the SSO drag  parametrization.

  ! purpose.
  ! --------

  ! this routine computes the physical tendencies of the
  ! prognostic variables u,v  and t due to  vertical transports by
  ! subgridscale orographically excited gravity waves, and to
  ! low level blocked flow drag.

  ! **   interface.
  ! ----------
  ! called from *drag_noro*.

  ! the routine takes its input from the long-term storage:
  ! u,v,t and p at t-1.

  ! explicit arguments :
  ! --------------------
  ! ==== inputs ===
  ! partdrag-input-I-control which part of the drag we consider (total part or GW part)
  ! nlon----input-I-Total number of horizontal points that get into physics
  ! nlev----input-I-Number of vertical levels

  ! kgwd- -input-I: Total nb of points where the orography schemes are active
  ! ktest--input-I: Flags to indicate active points
  ! kdx----input-I: Locate the physical location of an active point.
  ! ptsphy--input-R-Time-step (s)
  ! paphm1--input-R: pressure at model 1/2 layer
  ! papm1---input-R: pressure at model layer
  ! pgeom1--input-R: Altitude of layer above ground
  ! ptm1, pum1, pvm1--R-: t, u and v
  ! pmea----input-R-Mean Orography (m)
  ! pstd----input-R-SSO standard deviation (m)
  ! psig----input-R-SSO slope
  ! pgam----input-R-SSO Anisotropy
  ! pthe----input-R-SSO Angle
  ! ppic----input-R-SSO Peacks elevation (m)
  ! pval----input-R-SSO Valleys elevation (m)

  INTEGER  nlon, nlev, kgwd
  REAL ptsphy

  ! ==== outputs ===
  ! pulow, pvlow -output-R: Low-level wind

  ! pte -----output-R: T tendency
  ! pvom-----output-R: U tendency
  ! pvol-----output-R: V tendency


  ! Implicit Arguments:
  ! ===================

  ! klon-common-I: Number of points seen by the physics
  ! klev-common-I: Number of vertical layers

  ! method.
  ! -------

  ! externals.
  ! ----------
  INTEGER ismin, ismax
  EXTERNAL ismin, ismax

  ! reference.
  ! ----------

  ! author.
  ! -------
  ! m.miller + b.ritter   e.c.m.w.f.     15/06/86.

  ! f.lott + m. miller    e.c.m.w.f.     22/11/94
  ! -----------------------------------------------------------------------


  include "YOMCST.h"
  include "YOEGWD.h"

  ! -----------------------------------------------------------------------

  ! *       0.1   arguments
  ! ---------

  INTEGER partdrag
  REAL pte(nlon, nlev), pvol(nlon, nlev), pvom(nlon, nlev), pulow(nlon), &
    pvlow(nlon)
  REAL pum1(nlon, nlev), pvm1(nlon, nlev), ptm1(nlon, nlev), pmea(nlon), &
    pstd(nlon), psig(nlon), pgam(nlon), pthe(nlon), ppic(nlon), pval(nlon), &
    pgeom1(nlon, nlev), papm1(nlon, nlev), paphm1(nlon, nlev+1)

  INTEGER kdx(nlon), ktest(nlon)
  ! -----------------------------------------------------------------------

  ! *       0.2   local arrays
  ! ------------
  INTEGER isect(klon), icrit(klon), ikcrith(klon), ikenvh(klon), iknu(klon), &
    iknu2(klon), ikcrit(klon), ikhlim(klon)

  REAL ztau(klon, klev+1), zstab(klon, klev+1), zvph(klon, klev+1), &
    zrho(klon, klev+1), zri(klon, klev+1), zpsi(klon, klev+1), &
    zzdep(klon, klev)
  REAL zdudt(klon), zdvdt(klon), zdtdt(klon), zdedt(klon), zvidis(klon), &
    ztfr(klon), znu(klon), zd1(klon), zd2(klon), zdmod(klon)


  ! local quantities:

  INTEGER jl, jk, ji
  REAL ztmst, zdelp, ztemp, zforc, ztend, rover, facpart
  REAL zb, zc, zconb, zabsv, zzd1, ratio, zbet, zust, zvst, zdis

  ! ------------------------------------------------------------------

  ! *         1.    initialization
  ! --------------

  ! print *,' in orodrag'

  ! ------------------------------------------------------------------

  ! *         1.1   computational constants
  ! -----------------------


  ! ztmst=twodt
  ! if(nstep.eq.nstart) ztmst=0.5*twodt
  ztmst = ptsphy

  ! ------------------------------------------------------------------

  ! *         1.3   check whether row contains point for printing
  ! ---------------------------------------------


  ! ------------------------------------------------------------------

  ! *         2.     precompute basic state variables.
  ! *                ---------- ----- ----- ----------
  ! *                define low level wind, project winds in plane of
  ! *                low level wind, determine sector in which to take
  ! *                the variance and set indicator for critical levels.





  CALL orosetup_strato(nlon, nlev, ktest, ikcrit, ikcrith, icrit, isect, &
    ikhlim, ikenvh, iknu, iknu2, paphm1, papm1, pum1, pvm1, ptm1, pgeom1, &
    pstd, zrho, zri, zstab, ztau, zvph, zpsi, zzdep, pulow, pvlow, pthe, &
    pgam, pmea, ppic, pval, znu, zd1, zd2, zdmod)

  ! ***********************************************************


  ! *         3.      compute low level stresses using subcritical and
  ! *                 supercritical forms.computes anisotropy coefficient
  ! *                 as measure of orographic twodimensionality.


  CALL gwstress_strato(nlon, nlev, ikcrit, isect, ikhlim, ktest, ikcrith, &
    icrit, ikenvh, iknu, zrho, zstab, zvph, pstd, psig, pmea, ppic, pval, &
    ztfr, ztau, pgeom1, pgam, zd1, zd2, zdmod, znu)

  ! *         4.      compute stress profile including
  ! trapped waves, wave breaking,
  ! linear decay in stratosphere.




  CALL gwprofil_strato(nlon, nlev, kgwd, kdx, ktest, ikcrit, ikcrith, icrit, &
    ikenvh, iknu, iknu2, paphm1, zrho, zstab, ztfr, zvph, zri, ztau &
    , zdmod, znu, psig, pgam, pstd, ppic, pval)

  ! *         5.      Compute tendencies from waves stress profile.
  ! Compute low level blocked flow drag.
  ! *                 --------------------------------------------




  ! explicit solution at all levels for the gravity wave
  ! implicit solution for the blocked levels

  DO jl = kidia, kfdia
    zvidis(jl) = 0.0
    zdudt(jl) = 0.0
    zdvdt(jl) = 0.0
    zdtdt(jl) = 0.0
  END DO


  DO jk = 1, klev


    ! WAVE STRESS
    ! -------------


    DO ji = kidia, kfdia

      IF (ktest(ji)==1) THEN

        zdelp = paphm1(ji, jk+1) - paphm1(ji, jk)
        ztemp = -rg*(ztau(ji,jk+1)-ztau(ji,jk))/(zvph(ji,klev+1)*zdelp)

        zdudt(ji) = (pulow(ji)*zd1(ji)-pvlow(ji)*zd2(ji))*ztemp/zdmod(ji)
        zdvdt(ji) = (pvlow(ji)*zd1(ji)+pulow(ji)*zd2(ji))*ztemp/zdmod(ji)

        ! Control Overshoots


        IF (jk>=nstra) THEN
          rover = 0.10
          IF (abs(zdudt(ji))>rover*abs(pum1(ji,jk))/ztmst) zdudt(ji) = rover* &
            abs(pum1(ji,jk))/ztmst*zdudt(ji)/(abs(zdudt(ji))+1.E-10)
          IF (abs(zdvdt(ji))>rover*abs(pvm1(ji,jk))/ztmst) zdvdt(ji) = rover* &
            abs(pvm1(ji,jk))/ztmst*zdvdt(ji)/(abs(zdvdt(ji))+1.E-10)
        END IF

        rover = 0.25
        zforc = sqrt(zdudt(ji)**2+zdvdt(ji)**2)
        ztend = sqrt(pum1(ji,jk)**2+pvm1(ji,jk)**2)/ztmst

        IF (zforc>=rover*ztend) THEN
          zdudt(ji) = rover*ztend/zforc*zdudt(ji)
          zdvdt(ji) = rover*ztend/zforc*zdvdt(ji)
        END IF

        ! BLOCKED FLOW DRAG:
        ! -----------------

        IF (partdrag .GE. 2) THEN
        facpart=0.
        ELSE
        facpart=gkwake
        ENDIF


        IF (jk>ikenvh(ji)) THEN
          zb = 1.0 - 0.18*pgam(ji) - 0.04*pgam(ji)**2
          zc = 0.48*pgam(ji) + 0.3*pgam(ji)**2
          zconb = 2.*ztmst*facpart*psig(ji)/(4.*pstd(ji))
          zabsv = sqrt(pum1(ji,jk)**2+pvm1(ji,jk)**2)/2.
          zzd1 = zb*cos(zpsi(ji,jk))**2 + zc*sin(zpsi(ji,jk))**2
          ratio = (cos(zpsi(ji,jk))**2+pgam(ji)*sin(zpsi(ji, &
            jk))**2)/(pgam(ji)*cos(zpsi(ji,jk))**2+sin(zpsi(ji,jk))**2)
          zbet = max(0., 2.-1./ratio)*zconb*zzdep(ji, jk)*zzd1*zabsv

          ! OPPOSED TO THE WIND

          zdudt(ji) = -pum1(ji, jk)/ztmst
          zdvdt(ji) = -pvm1(ji, jk)/ztmst

          ! PERPENDICULAR TO THE SSO MAIN AXIS:

          ! mod     zdudt(ji)=-(pum1(ji,jk)*cos(pthe(ji)*rpi/180.)
          ! mod *              +pvm1(ji,jk)*sin(pthe(ji)*rpi/180.))
          ! mod *              *cos(pthe(ji)*rpi/180.)/ztmst
          ! mod     zdvdt(ji)=-(pum1(ji,jk)*cos(pthe(ji)*rpi/180.)
          ! mod *              +pvm1(ji,jk)*sin(pthe(ji)*rpi/180.))
          ! mod *              *sin(pthe(ji)*rpi/180.)/ztmst

          zdudt(ji) = zdudt(ji)*(zbet/(1.+zbet))
          zdvdt(ji) = zdvdt(ji)*(zbet/(1.+zbet))
        END IF
        pvom(ji, jk) = zdudt(ji)
        pvol(ji, jk) = zdvdt(ji)
        zust = pum1(ji, jk) + ztmst*zdudt(ji)
        zvst = pvm1(ji, jk) + ztmst*zdvdt(ji)
        zdis = 0.5*(pum1(ji,jk)**2+pvm1(ji,jk)**2-zust**2-zvst**2)
        zdedt(ji) = zdis/ztmst
        zvidis(ji) = zvidis(ji) + zdis*zdelp
        zdtdt(ji) = zdedt(ji)/rcpd

        ! NO TENDENCIES ON TEMPERATURE .....

        ! Instead of, pte(ji,jk)=zdtdt(ji), due to mechanical dissipation

        pte(ji, jk) = 0.0

      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE orodrag_strato
SUBROUTINE orosetup_strato(nlon, nlev, ktest, kkcrit, kkcrith, kcrit, ksect, &
    kkhlim, kkenvh, kknu, kknu2, paphm1, papm1, pum1, pvm1, ptm1, pgeom1, &
    pstd, prho, pri, pstab, ptau, pvph, ppsi, pzdep, pulow, pvlow, ptheta, &
    pgam, pmea, ppic, pval, pnu, pd1, pd2, pdmod)

  ! **** *gwsetup*

  ! purpose.
  ! --------
  ! SET-UP THE ESSENTIAL PARAMETERS OF THE SSO DRAG SCHEME:
  ! DEPTH OF LOW WBLOCKED LAYER, LOW-LEVEL FLOW, BACKGROUND
  ! STRATIFICATION.....

  ! **   interface.
  ! ----------
  ! from *orodrag*

  ! explicit arguments :
  ! --------------------
  ! ==== inputs ===

  ! nlon----input-I-Total number of horizontal points that get into physics
  ! nlev----input-I-Number of vertical levels
  ! ktest--input-I: Flags to indicate active points

  ! ptsphy--input-R-Time-step (s)
  ! paphm1--input-R: pressure at model 1/2 layer
  ! papm1---input-R: pressure at model layer
  ! pgeom1--input-R: Altitude of layer above ground
  ! ptm1, pum1, pvm1--R-: t, u and v
  ! pmea----input-R-Mean Orography (m)
  ! pstd----input-R-SSO standard deviation (m)
  ! psig----input-R-SSO slope
  ! pgam----input-R-SSO Anisotropy
  ! pthe----input-R-SSO Angle
  ! ppic----input-R-SSO Peacks elevation (m)
  ! pval----input-R-SSO Valleys elevation (m)

  ! ==== outputs ===
  ! pulow, pvlow -output-R: Low-level wind
  ! kkcrit----I-: Security value for top of low level flow
  ! kcrit-----I-: Critical level
  ! ksect-----I-: Not used
  ! kkhlim----I-: Not used
  ! kkenvh----I-: Top of blocked flow layer
  ! kknu------I-: Layer that sees mountain peacks
  ! kknu2-----I-: Layer that sees mountain peacks above mountain mean
  ! kknub-----I-: Layer that sees mountain mean above valleys
  ! prho------R-: Density at 1/2 layers
  ! pri-------R-: Background Richardson Number, Wind shear measured along GW
  ! stress
  ! pstab-----R-: Brunt-Vaisala freq. at 1/2 layers
  ! pvph------R-: Wind in  plan of GW stress, Half levels.
  ! ppsi------R-: Angle between low level wind and SS0 main axis.
  ! pd1-------R-| Compared the ratio of the stress
  ! pd2-------R-| that is along the wind to that Normal to it.
  ! pdi define the plane of low level stress
  ! compared to the low level wind.
  ! see p. 108 Lott & Miller (1997).
  ! pdmod-----R-: Norme of pdi

  ! === local arrays ===

  ! zvpf------R-: Wind projected in the plan of the low-level stress.

  ! ==== outputs ===

  ! implicit arguments :   none
  ! --------------------

  ! method.
  ! -------


  ! externals.
  ! ----------


  ! reference.
  ! ----------

  ! see ecmwf research department documentation of the "i.f.s."

  ! author.
  ! -------

  ! modifications.
  ! --------------
  ! f.lott  for the new-gwdrag scheme november 1993

  ! -----------------------------------------------------------------------
  USE dimphy
  IMPLICIT NONE


  include "YOMCST.h"
  include "YOEGWD.h"

  ! -----------------------------------------------------------------------

  ! *       0.1   arguments
  ! ---------

  INTEGER nlon, nlev
  INTEGER kkcrit(nlon), kkcrith(nlon), kcrit(nlon), ksect(nlon), &
    kkhlim(nlon), ktest(nlon), kkenvh(nlon)


  REAL paphm1(nlon, klev+1), papm1(nlon, klev), pum1(nlon, klev), &
    pvm1(nlon, klev), ptm1(nlon, klev), pgeom1(nlon, klev), &
    prho(nlon, klev+1), pri(nlon, klev+1), pstab(nlon, klev+1), &
    ptau(nlon, klev+1), pvph(nlon, klev+1), ppsi(nlon, klev+1), &
    pzdep(nlon, klev)
  REAL pulow(nlon), pvlow(nlon), ptheta(nlon), pgam(nlon), pnu(nlon), &
    pd1(nlon), pd2(nlon), pdmod(nlon)
  REAL pstd(nlon), pmea(nlon), ppic(nlon), pval(nlon)

  ! -----------------------------------------------------------------------

  ! *       0.2   local arrays
  ! ------------


  INTEGER ilevh, jl, jk
  REAL zcons1, zcons2, zhgeo, zu, zphi
  REAL zvt1, zvt2, zdwind, zwind, zdelp
  REAL zstabm, zstabp, zrhom, zrhop
  LOGICAL lo
  LOGICAL ll1(klon, klev+1)
  INTEGER kknu(klon), kknu2(klon), kknub(klon), kknul(klon), kentp(klon), &
    ncount(klon)

  REAL zhcrit(klon, klev), zvpf(klon, klev), zdp(klon, klev)
  REAL znorm(klon), zb(klon), zc(klon), zulow(klon), zvlow(klon), znup(klon), &
    znum(klon)

  ! ------------------------------------------------------------------

  ! *         1.    initialization
  ! --------------

  ! PRINT *,' in orosetup'

  ! ------------------------------------------------------------------

  ! *         1.1   computational constants
  ! -----------------------


  ilevh = klev/3

  zcons1 = 1./rd
  zcons2 = rg**2/rcpd

  ! ------------------------------------------------------------------

  ! *         2.
  ! --------------


  ! ------------------------------------------------------------------

  ! *         2.1     define low level wind, project winds in plane of
  ! *                 low level wind, determine sector in which to take
  ! *                 the variance and set indicator for critical levels.



  DO jl = kidia, kfdia
    kknu(jl) = klev
    kknu2(jl) = klev
    kknub(jl) = klev
    kknul(jl) = klev
    pgam(jl) = max(pgam(jl), gtsec)
    ll1(jl, klev+1) = .FALSE.
  END DO

  ! Ajouter une initialisation (L. Li, le 23fev99):

  DO jk = klev, ilevh, -1
    DO jl = kidia, kfdia
      ll1(jl, jk) = .FALSE.
    END DO
  END DO

  ! *      define top of low level flow
  ! ----------------------------
  DO jk = klev, ilevh, -1
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        lo = (paphm1(jl,jk)/paphm1(jl,klev+1)) >= gsigcr
        IF (lo) THEN
          kkcrit(jl) = jk
        END IF
        zhcrit(jl, jk) = ppic(jl) - pval(jl)
        zhgeo = pgeom1(jl, jk)/rg
        ll1(jl, jk) = (zhgeo>zhcrit(jl,jk))
        IF (ll1(jl,jk) .NEQV. ll1(jl,jk+1)) THEN
          kknu(jl) = jk
        END IF
        IF (.NOT. ll1(jl,ilevh)) kknu(jl) = ilevh
      END IF
    END DO
  END DO
  DO jk = klev, ilevh, -1
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        zhcrit(jl, jk) = ppic(jl) - pmea(jl)
        zhgeo = pgeom1(jl, jk)/rg
        ll1(jl, jk) = (zhgeo>zhcrit(jl,jk))
        IF (ll1(jl,jk) .NEQV. ll1(jl,jk+1)) THEN
          kknu2(jl) = jk
        END IF
        IF (.NOT. ll1(jl,ilevh)) kknu2(jl) = ilevh
      END IF
    END DO
  END DO
  DO jk = klev, ilevh, -1
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        zhcrit(jl, jk) = amin1(ppic(jl)-pmea(jl), pmea(jl)-pval(jl))
        zhgeo = pgeom1(jl, jk)/rg
        ll1(jl, jk) = (zhgeo>zhcrit(jl,jk))
        IF (ll1(jl,jk) .NEQV. ll1(jl,jk+1)) THEN
          kknub(jl) = jk
        END IF
        IF (.NOT. ll1(jl,ilevh)) kknub(jl) = ilevh
      END IF
    END DO
  END DO

  DO jl = kidia, kfdia
    IF (ktest(jl)==1) THEN
      kknu(jl) = min(kknu(jl), nktopg)
      kknu2(jl) = min(kknu2(jl), nktopg)
      kknub(jl) = min(kknub(jl), nktopg)
      kknul(jl) = klev
    END IF
  END DO

  ! c*     initialize various arrays

  DO jl = kidia, kfdia
    prho(jl, klev+1) = 0.0
    ! ym correction en attendant mieux
    prho(jl, 1) = 0.0
    pstab(jl, klev+1) = 0.0
    pstab(jl, 1) = 0.0
    pri(jl, klev+1) = 9999.0
    ppsi(jl, klev+1) = 0.0
    pri(jl, 1) = 0.0
    pvph(jl, 1) = 0.0
    pvph(jl, klev+1) = 0.0
    ! ym correction en attendant mieux
    ! ym      pvph(jl,klev)    =0.0
    pulow(jl) = 0.0
    pvlow(jl) = 0.0
    zulow(jl) = 0.0
    zvlow(jl) = 0.0
    kkcrith(jl) = klev
    kkenvh(jl) = klev
    kentp(jl) = klev
    kcrit(jl) = 1
    ncount(jl) = 0
    ll1(jl, klev+1) = .FALSE.
  END DO

  ! *     define flow density and stratification (rho and N2)
  ! at semi layers.
  ! -------------------------------------------------------

  DO jk = klev, 2, -1
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        zdp(jl, jk) = papm1(jl, jk) - papm1(jl, jk-1)
        prho(jl, jk) = 2.*paphm1(jl, jk)*zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
        pstab(jl, jk) = 2.*zcons2/(ptm1(jl,jk)+ptm1(jl,jk-1))* &
          (1.-rcpd*prho(jl,jk)*(ptm1(jl,jk)-ptm1(jl,jk-1))/zdp(jl,jk))
        pstab(jl, jk) = max(pstab(jl,jk), gssec)
      END IF
    END DO
  END DO

  ! ********************************************************************

  ! *     define Low level flow (between ground and peacks-valleys)
  ! ---------------------------------------------------------
  DO jk = klev, ilevh, -1
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        IF (jk>=kknu2(jl) .AND. jk<=kknul(jl)) THEN
          pulow(jl) = pulow(jl) + pum1(jl, jk)*(paphm1(jl,jk+1)-paphm1(jl,jk) &
            )
          pvlow(jl) = pvlow(jl) + pvm1(jl, jk)*(paphm1(jl,jk+1)-paphm1(jl,jk) &
            )
          pstab(jl, klev+1) = pstab(jl, klev+1) + pstab(jl, jk)*(paphm1(jl,jk &
            +1)-paphm1(jl,jk))
          prho(jl, klev+1) = prho(jl, klev+1) + prho(jl, jk)*(paphm1(jl,jk+1) &
            -paphm1(jl,jk))
        END IF
      END IF
    END DO
  END DO
  DO jl = kidia, kfdia
    IF (ktest(jl)==1) THEN
      pulow(jl) = pulow(jl)/(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknu2(jl)))
      pvlow(jl) = pvlow(jl)/(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknu2(jl)))
      znorm(jl) = max(sqrt(pulow(jl)**2+pvlow(jl)**2), gvsec)
      pvph(jl, klev+1) = znorm(jl)
      pstab(jl, klev+1) = pstab(jl, klev+1)/(paphm1(jl,kknul(jl)+1)-paphm1(jl &
        ,kknu2(jl)))
      prho(jl, klev+1) = prho(jl, klev+1)/(paphm1(jl,kknul(jl)+1)-paphm1(jl, &
        kknu2(jl)))
    END IF
  END DO


  ! *******  setup orography orientation relative to the low level
  ! wind and define parameters of the Anisotropic wave stress.

  DO jl = kidia, kfdia
    IF (ktest(jl)==1) THEN
      lo = (pulow(jl)<gvsec) .AND. (pulow(jl)>=-gvsec)
      IF (lo) THEN
        zu = pulow(jl) + 2.*gvsec
      ELSE
        zu = pulow(jl)
      END IF
      zphi = atan(pvlow(jl)/zu)
      ppsi(jl, klev+1) = ptheta(jl)*rpi/180. - zphi
      zb(jl) = 1. - 0.18*pgam(jl) - 0.04*pgam(jl)**2
      zc(jl) = 0.48*pgam(jl) + 0.3*pgam(jl)**2
      pd1(jl) = zb(jl) - (zb(jl)-zc(jl))*(sin(ppsi(jl,klev+1))**2)
      pd2(jl) = (zb(jl)-zc(jl))*sin(ppsi(jl,klev+1))*cos(ppsi(jl,klev+1))
      pdmod(jl) = sqrt(pd1(jl)**2+pd2(jl)**2)
    END IF
  END DO

  ! ************ projet flow in plane of lowlevel stress *************
  ! ************ Find critical levels...                 *************

  DO jk = 1, klev
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        zvt1 = pulow(jl)*pum1(jl, jk) + pvlow(jl)*pvm1(jl, jk)
        zvt2 = -pvlow(jl)*pum1(jl, jk) + pulow(jl)*pvm1(jl, jk)
        zvpf(jl, jk) = (zvt1*pd1(jl)+zvt2*pd2(jl))/(znorm(jl)*pdmod(jl))
      END IF
      ptau(jl, jk) = 0.0
      pzdep(jl, jk) = 0.0
      ppsi(jl, jk) = 0.0
      ll1(jl, jk) = .FALSE.
    END DO
  END DO
  DO jk = 2, klev
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        zdp(jl, jk) = papm1(jl, jk) - papm1(jl, jk-1)
        pvph(jl, jk) = ((paphm1(jl,jk)-papm1(jl,jk-1))*zvpf(jl,jk)+(papm1(jl, &
          jk)-paphm1(jl,jk))*zvpf(jl,jk-1))/zdp(jl, jk)
        IF (pvph(jl,jk)<gvsec) THEN
          pvph(jl, jk) = gvsec
          kcrit(jl) = jk
        END IF
      END IF
    END DO
  END DO

  ! *         2.3     mean flow richardson number.


  DO jk = 2, klev
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        zdwind = max(abs(zvpf(jl,jk)-zvpf(jl,jk-1)), gvsec)
        pri(jl, jk) = pstab(jl, jk)*(zdp(jl,jk)/(rg*prho(jl,jk)*zdwind))**2
        pri(jl, jk) = max(pri(jl,jk), grcrit)
      END IF
    END DO
  END DO



  ! *      define top of 'envelope' layer
  ! ----------------------------

  DO jl = kidia, kfdia
    pnu(jl) = 0.0
    znum(jl) = 0.0
  END DO

  DO jk = 2, klev - 1
    DO jl = kidia, kfdia

      IF (ktest(jl)==1) THEN

        IF (jk>=kknu2(jl)) THEN

          znum(jl) = pnu(jl)
          zwind = (pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/ &
            max(sqrt(pulow(jl)**2+pvlow(jl)**2), gvsec)
          zwind = max(sqrt(zwind**2), gvsec)
          zdelp = paphm1(jl, jk+1) - paphm1(jl, jk)
          zstabm = sqrt(max(pstab(jl,jk),gssec))
          zstabp = sqrt(max(pstab(jl,jk+1),gssec))
          zrhom = prho(jl, jk)
          zrhop = prho(jl, jk+1)
          pnu(jl) = pnu(jl) + (zdelp/rg)*((zstabp/zrhop+zstabm/zrhom)/2.)/ &
            zwind
          IF ((znum(jl)<=gfrcrit) .AND. (pnu(jl)>gfrcrit) .AND. (kkenvh( &
            jl)==klev)) kkenvh(jl) = jk

        END IF

      END IF

    END DO
  END DO

  ! calculation of a dynamical mixing height for when the waves
  ! BREAK AT LOW LEVEL: The drag will be repartited over
  ! a depths that depends on waves vertical wavelength,
  ! not just between two adjacent model layers.
  ! of gravity waves:

  DO jl = kidia, kfdia
    znup(jl) = 0.0
    znum(jl) = 0.0
  END DO

  DO jk = klev - 1, 2, -1
    DO jl = kidia, kfdia

      IF (ktest(jl)==1) THEN

        znum(jl) = znup(jl)
        zwind = (pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/ &
          max(sqrt(pulow(jl)**2+pvlow(jl)**2), gvsec)
        zwind = max(sqrt(zwind**2), gvsec)
        zdelp = paphm1(jl, jk+1) - paphm1(jl, jk)
        zstabm = sqrt(max(pstab(jl,jk),gssec))
        zstabp = sqrt(max(pstab(jl,jk+1),gssec))
        zrhom = prho(jl, jk)
        zrhop = prho(jl, jk+1)
        znup(jl) = znup(jl) + (zdelp/rg)*((zstabp/zrhop+zstabm/zrhom)/2.)/ &
          zwind
        IF ((znum(jl)<=rpi/4.) .AND. (znup(jl)>rpi/4.) .AND. (kkcrith( &
          jl)==klev)) kkcrith(jl) = jk

      END IF

    END DO
  END DO

  DO jl = kidia, kfdia
    IF (ktest(jl)==1) THEN
      kkcrith(jl) = max0(kkcrith(jl), ilevh*2)
      kkcrith(jl) = max0(kkcrith(jl), kknu(jl))
      IF (kcrit(jl)>=kkcrith(jl)) kcrit(jl) = 1
    END IF
  END DO

  ! directional info for flow blocking *************************

  DO jk = 1, klev
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        lo = (pum1(jl,jk)<gvsec) .AND. (pum1(jl,jk)>=-gvsec)
        IF (lo) THEN
          zu = pum1(jl, jk) + 2.*gvsec
        ELSE
          zu = pum1(jl, jk)
        END IF
        zphi = atan(pvm1(jl,jk)/zu)
        ppsi(jl, jk) = ptheta(jl)*rpi/180. - zphi
      END IF
    END DO
  END DO

  ! forms the vertical 'leakiness' **************************

  DO jk = ilevh, klev
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        pzdep(jl, jk) = 0
        IF (jk>=kkenvh(jl) .AND. kkenvh(jl)/=klev) THEN
          pzdep(jl, jk) = (pgeom1(jl,kkenvh(jl))-pgeom1(jl,jk))/ &
            (pgeom1(jl,kkenvh(jl))-pgeom1(jl,klev))
        END IF
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE orosetup_strato
SUBROUTINE gwstress_strato(nlon, nlev, kkcrit, ksect, kkhlim, ktest, kkcrith, &
    kcrit, kkenvh, kknu, prho, pstab, pvph, pstd, psig, pmea, ppic, pval, &
    ptfr, ptau, pgeom1, pgamma, pd1, pd2, pdmod, pnu)

  ! **** *gwstress*

  ! purpose.
  ! --------
  ! Compute the surface stress due to Gravity Waves, according
  ! to the Phillips (1979) theory of 3-D flow above
  ! anisotropic elliptic ridges.

  ! The stress is reduced two account for cut-off flow over
  ! hill.  The flow only see that part of the ridge located
  ! above the blocked layer (see zeff).

  ! **   interface.
  ! ----------
  ! call *gwstress*  from *gwdrag*

  ! explicit arguments :
  ! --------------------
  ! ==== inputs ===
  ! ==== outputs ===

  ! implicit arguments :   none
  ! --------------------

  ! method.
  ! -------


  ! externals.
  ! ----------


  ! reference.
  ! ----------

  ! LOTT and MILLER (1997)  &  LOTT (1999)

  ! author.
  ! -------

  ! modifications.
  ! --------------
  ! f. lott put the new gwd on ifs      22/11/93

  ! -----------------------------------------------------------------------
  USE dimphy
  IMPLICIT NONE

  include "YOMCST.h"
  include "YOEGWD.h"

  ! -----------------------------------------------------------------------

  ! *       0.1   arguments
  ! ---------

  INTEGER nlon, nlev
  INTEGER kkcrit(nlon), kkcrith(nlon), kcrit(nlon), ksect(nlon), &
    kkhlim(nlon), ktest(nlon), kkenvh(nlon), kknu(nlon)

  REAL prho(nlon, nlev+1), pstab(nlon, nlev+1), ptau(nlon, nlev+1), &
    pvph(nlon, nlev+1), ptfr(nlon), pgeom1(nlon, nlev), pstd(nlon)

  REAL pd1(nlon), pd2(nlon), pnu(nlon), psig(nlon), pgamma(nlon)
  REAL pmea(nlon), ppic(nlon), pval(nlon)
  REAL pdmod(nlon)

  ! -----------------------------------------------------------------------

  ! *       0.2   local arrays
  ! ------------
  ! zeff--real: effective height seen by the flow when there is blocking

  INTEGER jl
  REAL zeff

  ! -----------------------------------------------------------------------

  ! *       0.3   functions
  ! ---------
  ! ------------------------------------------------------------------

  ! *         1.    initialization
  ! --------------

  ! PRINT *,' in gwstress'

  ! *         3.1     gravity wave stress.



  DO jl = kidia, kfdia
    IF (ktest(jl)==1) THEN

      ! effective mountain height above the blocked flow

      zeff = ppic(jl) - pval(jl)
      IF (kkenvh(jl)<klev) THEN
        zeff = amin1(gfrcrit*pvph(jl,klev+1)/sqrt(pstab(jl,klev+1)), zeff)
      END IF


      ptau(jl, klev+1) = gkdrag*prho(jl, klev+1)*psig(jl)*pdmod(jl)/4./ &
        pstd(jl)*pvph(jl, klev+1)*sqrt(pstab(jl,klev+1))*zeff**2


      ! too small value of stress or  low level flow include critical level
      ! or low level flow:  gravity wave stress nul.

      ! lo=(ptau(jl,klev+1).lt.gtsec).or.(kcrit(jl).ge.kknu(jl))
      ! *      .or.(pvph(jl,klev+1).lt.gvcrit)
      ! if(lo) ptau(jl,klev+1)=0.0

      ! print *,jl,ptau(jl,klev+1)

    ELSE

      ptau(jl, klev+1) = 0.0

    END IF

  END DO

  ! write(21)(ptau(jl,klev+1),jl=kidia,kfdia)

  RETURN
END SUBROUTINE gwstress_strato

SUBROUTINE gwprofil_strato(nlon, nlev, kgwd, kdx, ktest, kkcrit, kkcrith, &
    kcrit, kkenvh, kknu, kknu2, paphm1, prho, pstab, ptfr, pvph, pri, ptau, &
    pdmod, pnu, psig, pgamma, pstd, ppic, pval)

  ! **** *gwprofil*

  ! purpose.
  ! --------

  ! **   interface.
  ! ----------
  ! from *gwdrag*

  ! explicit arguments :
  ! --------------------
  ! ==== inputs ===

  ! ==== outputs ===

  ! implicit arguments :   none
  ! --------------------

  ! method:
  ! -------
  ! the stress profile for gravity waves is computed as follows:
  ! it decreases linearly with heights from the ground
  ! to the low-level indicated by kkcrith,
  ! to simulates lee waves or
  ! low-level gravity wave breaking.
  ! above it is constant, except when the waves encounter a critical
  ! level (kcrit) or when they break.
  ! The stress is also uniformly distributed above the level
  ! nstra.

  USE dimphy
  IMPLICIT NONE

  include "YOMCST.h"
  include "YOEGWD.h"

  ! -----------------------------------------------------------------------

  ! *       0.1   ARGUMENTS
  ! ---------

  INTEGER nlon, nlev, kgwd
  INTEGER kkcrit(nlon), kkcrith(nlon), kcrit(nlon), kdx(nlon), ktest(nlon), &
    kkenvh(nlon), kknu(nlon), kknu2(nlon)

  REAL paphm1(nlon, nlev+1), pstab(nlon, nlev+1), prho(nlon, nlev+1), &
    pvph(nlon, nlev+1), pri(nlon, nlev+1), ptfr(nlon), ptau(nlon, nlev+1)

  REAL pdmod(nlon), pnu(nlon), psig(nlon), pgamma(nlon), pstd(nlon), &
    ppic(nlon), pval(nlon)

  ! -----------------------------------------------------------------------

  ! *       0.2   local arrays
  ! ------------

  INTEGER jl, jk
  REAL zsqr, zalfa, zriw, zdel, zb, zalpha, zdz2n, zdelp, zdelpt

  REAL zdz2(klon, klev), znorm(klon), zoro(klon)
  REAL ztau(klon, klev+1)

  ! -----------------------------------------------------------------------

  ! *         1.    INITIALIZATION
  ! --------------

  ! print *,' entree gwprofil'


  ! *    COMPUTATIONAL CONSTANTS.
  ! ------------- ----------

  DO jl = kidia, kfdia
    IF (ktest(jl)==1) THEN
      zoro(jl) = psig(jl)*pdmod(jl)/4./pstd(jl)
      ztau(jl, klev+1) = ptau(jl, klev+1)
      ! print *,jl,ptau(jl,klev+1)
      ztau(jl, kkcrith(jl)) = grahilo*ptau(jl, klev+1)
    END IF
  END DO


  DO jk = klev + 1, 1, -1
    ! *         4.1    constant shear stress until top of the
    ! low-level breaking/trapped layer

    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        IF (jk>kkcrith(jl)) THEN
          zdelp = paphm1(jl, jk) - paphm1(jl, klev+1)
          zdelpt = paphm1(jl, kkcrith(jl)) - paphm1(jl, klev+1)
          ptau(jl, jk) = ztau(jl, klev+1) + zdelp/zdelpt*(ztau(jl,kkcrith(jl) &
            )-ztau(jl,klev+1))
        ELSE
          ptau(jl, jk) = ztau(jl, kkcrith(jl))
        END IF
      END IF
    END DO

    ! *         4.15   constant shear stress until the top of the
    ! low level flow layer.


    ! *         4.2    wave displacement at next level.


  END DO


  ! *         4.4    wave richardson number, new wave displacement
  ! *                and stress:  breaking evaluation and critical
  ! level


  DO jk = klev, 1, -1

    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        znorm(jl) = prho(jl, jk)*sqrt(pstab(jl,jk))*pvph(jl, jk)
        zdz2(jl, jk) = ptau(jl, jk)/amax1(znorm(jl), gssec)/zoro(jl)
      END IF
    END DO

    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        IF (jk<kkcrith(jl)) THEN
          IF ((ptau(jl,jk+1)<gtsec) .OR. (jk<=kcrit(jl))) THEN
            ptau(jl, jk) = 0.0
          ELSE
            zsqr = sqrt(pri(jl,jk))
            zalfa = sqrt(pstab(jl,jk)*zdz2(jl,jk))/pvph(jl, jk)
            zriw = pri(jl, jk)*(1.-zalfa)/(1+zalfa*zsqr)**2
            IF (zriw<grcrit) THEN
              ! print *,' breaking!!!',ptau(jl,jk)
              zdel = 4./zsqr/grcrit + 1./grcrit**2 + 4./grcrit
              zb = 1./grcrit + 2./zsqr
              zalpha = 0.5*(-zb+sqrt(zdel))
              zdz2n = (pvph(jl,jk)*zalpha)**2/pstab(jl, jk)
              ptau(jl, jk) = znorm(jl)*zdz2n*zoro(jl)
            END IF

            ptau(jl, jk) = amin1(ptau(jl,jk), ptau(jl,jk+1))

          END IF
        END IF
      END IF
    END DO
  END DO

  ! REORGANISATION OF THE STRESS PROFILE AT LOW LEVEL

  DO jl = kidia, kfdia
    IF (ktest(jl)==1) THEN
      ztau(jl, kkcrith(jl)-1) = ptau(jl, kkcrith(jl)-1)
      ztau(jl, nstra) = ptau(jl, nstra)
    END IF
  END DO

  DO jk = 1, klev

    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN

        IF (jk>kkcrith(jl)-1) THEN

          zdelp = paphm1(jl, jk) - paphm1(jl, klev+1)
          zdelpt = paphm1(jl, kkcrith(jl)-1) - paphm1(jl, klev+1)
          ptau(jl, jk) = ztau(jl, klev+1) + (ztau(jl,kkcrith(jl)-1)-ztau(jl, &
            klev+1))*zdelp/zdelpt

        END IF
      END IF

    END DO

    ! REORGANISATION AT THE MODEL TOP....

    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN

        IF (jk<nstra) THEN

          zdelp = paphm1(jl, nstra)
          zdelpt = paphm1(jl, jk)
          ptau(jl, jk) = ztau(jl, nstra)*zdelpt/zdelp
          ! ptau(jl,jk)=ztau(jl,nstra)

        END IF

      END IF

    END DO


  END DO


123 FORMAT (I4, 1X, 20(F6.3,1X))


  RETURN
END SUBROUTINE gwprofil_strato
SUBROUTINE lift_noro_strato(nlon, nlev, dtime, paprs, pplay, plat, pmea, &
    pstd, psig, pgam, pthe, ppic, pval, kgwd, kdx, ktest, t, u, v, pulow, &
    pvlow, pustr, pvstr, d_t, d_u, d_v)

  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): F.Lott (LMD/CNRS) date: 19950201
  ! Object: Mountain lift interface (enhanced vortex stretching).
  ! Made necessary because:
  ! 1. in the LMD-GCM Layers are from bottom to top,
  ! contrary to most European GCM.
  ! 2. the altitude above ground of each model layers
  ! needs to be known (variable zgeom)
  ! ======================================================================
  ! Explicit Arguments:
  ! ==================
  ! nlon----input-I-Total number of horizontal points that get into physics
  ! nlev----input-I-Number of vertical levels
  ! dtime---input-R-Time-step (s)
  ! paprs---input-R-Pressure in semi layers    (Pa)
  ! pplay---input-R-Pressure model-layers      (Pa)
  ! t-------input-R-temperature (K)
  ! u-------input-R-Horizontal wind (m/s)
  ! v-------input-R-Meridional wind (m/s)
  ! pmea----input-R-Mean Orography (m)
  ! pstd----input-R-SSO standard deviation (m)
  ! psig----input-R-SSO slope
  ! pgam----input-R-SSO Anisotropy
  ! pthe----input-R-SSO Angle
  ! ppic----input-R-SSO Peacks elevation (m)
  ! pval----input-R-SSO Valleys elevation (m)

  ! kgwd- -input-I: Total nb of points where the orography schemes are active
  ! ktest--input-I: Flags to indicate active points
  ! kdx----input-I: Locate the physical location of an active point.

  ! pulow, pvlow -output-R: Low-level wind
  ! pustr, pvstr -output-R: Surface stress due to SSO drag      (Pa)

  ! d_t-----output-R: T increment
  ! d_u-----output-R: U increment
  ! d_v-----output-R: V increment

  ! Implicit Arguments:
  ! ===================

  ! iim--common-I: Number of longitude intervals
  ! jjm--common-I: Number of latitude intervals
  ! klon-common-I: Number of points seen by the physics
  ! (iim+1)*(jjm+1) for instance
  ! klev-common-I: Number of vertical layers
  ! ======================================================================
  ! Local Variables:
  ! ================

  ! zgeom-----R: Altitude of layer above ground
  ! pt, pu, pv --R: t u v from top to bottom
  ! pdtdt, pdudt, pdvdt --R: t u v tendencies (from top to bottom)
  ! papmf: pressure at model layer (from top to bottom)
  ! papmh: pressure at model 1/2 layer (from top to bottom)

  ! ======================================================================

  include "YOMCST.h"
  include "YOEGWD.h"

  ! ARGUMENTS

  INTEGER nlon, nlev
  REAL dtime
  REAL paprs(klon, klev+1)
  REAL pplay(klon, klev)
  REAL plat(nlon), pmea(nlon)
  REAL pstd(nlon), psig(nlon), pgam(nlon), pthe(nlon)
  REAL ppic(nlon), pval(nlon)
  REAL pulow(nlon), pvlow(nlon), pustr(nlon), pvstr(nlon)
  REAL t(nlon, nlev), u(nlon, nlev), v(nlon, nlev)
  REAL d_t(nlon, nlev), d_u(nlon, nlev), d_v(nlon, nlev)

  INTEGER i, k, kgwd, kdx(nlon), ktest(nlon)

  ! Variables locales:

  REAL zgeom(klon, klev)
  REAL pdtdt(klon, klev), pdudt(klon, klev), pdvdt(klon, klev)
  REAL pt(klon, klev), pu(klon, klev), pv(klon, klev)
  REAL papmf(klon, klev), papmh(klon, klev+1)

  ! initialiser les variables de sortie (pour securite)


  ! print *,'in lift_noro'
  DO i = 1, klon
    pulow(i) = 0.0
    pvlow(i) = 0.0
    pustr(i) = 0.0
    pvstr(i) = 0.0
  END DO
  DO k = 1, klev
    DO i = 1, klon
      d_t(i, k) = 0.0
      d_u(i, k) = 0.0
      d_v(i, k) = 0.0
      pdudt(i, k) = 0.0
      pdvdt(i, k) = 0.0
      pdtdt(i, k) = 0.0
    END DO
  END DO

  ! preparer les variables d'entree (attention: l'ordre des niveaux
  ! verticaux augmente du haut vers le bas)

  DO k = 1, klev
    DO i = 1, klon
      pt(i, k) = t(i, klev-k+1)
      pu(i, k) = u(i, klev-k+1)
      pv(i, k) = v(i, klev-k+1)
      papmf(i, k) = pplay(i, klev-k+1)
    END DO
  END DO
  DO k = 1, klev + 1
    DO i = 1, klon
      papmh(i, k) = paprs(i, klev-k+2)
    END DO
  END DO
  DO i = 1, klon
    zgeom(i, klev) = rd*pt(i, klev)*log(papmh(i,klev+1)/papmf(i,klev))
  END DO
  DO k = klev - 1, 1, -1
    DO i = 1, klon
      zgeom(i, k) = zgeom(i, k+1) + rd*(pt(i,k)+pt(i,k+1))/2.0*log(papmf(i,k+ &
        1)/papmf(i,k))
    END DO
  END DO

  ! appeler la routine principale


  CALL orolift_strato(klon, klev, kgwd, kdx, ktest, dtime, papmh, papmf, &
    zgeom, pt, pu, pv, plat, pmea, pstd, psig, pgam, pthe, ppic, pval, pulow, &
    pvlow, pdudt, pdvdt, pdtdt)

  DO k = 1, klev
    DO i = 1, klon
      d_u(i, klev+1-k) = dtime*pdudt(i, k)
      d_v(i, klev+1-k) = dtime*pdvdt(i, k)
      d_t(i, klev+1-k) = dtime*pdtdt(i, k)
      pustr(i) = pustr(i) + pdudt(i, k)*(papmh(i,k+1)-papmh(i,k))/rg
      pvstr(i) = pvstr(i) + pdvdt(i, k)*(papmh(i,k+1)-papmh(i,k))/rg
    END DO
  END DO

  ! print *,' out lift_noro'

  RETURN
END SUBROUTINE lift_noro_strato
SUBROUTINE orolift_strato(nlon, nlev, kgwd, kdx, ktest, ptsphy, paphm1, &
    papm1, pgeom1, ptm1, pum1, pvm1, plat, pmea, pstd, psig, pgam, pthe, &
    ppic, pval &                   ! OUTPUTS
    , pulow, pvlow, pvom, pvol, pte)


  ! **** *OROLIFT: SIMULATE THE GEOSTROPHIC LIFT.

  ! PURPOSE.
  ! --------
  ! this routine computes the physical tendencies of the
  ! prognostic variables u,v  when enhanced vortex stretching
  ! is needed.

  ! **   INTERFACE.
  ! ----------
  ! CALLED FROM *lift_noro
  ! explicit arguments :
  ! --------------------
  ! ==== inputs ===
  ! nlon----input-I-Total number of horizontal points that get into physics
  ! nlev----input-I-Number of vertical levels

  ! kgwd- -input-I: Total nb of points where the orography schemes are active
  ! ktest--input-I: Flags to indicate active points
  ! kdx----input-I: Locate the physical location of an active point.
  ! ptsphy--input-R-Time-step (s)
  ! paphm1--input-R: pressure at model 1/2 layer
  ! papm1---input-R: pressure at model layer
  ! pgeom1--input-R: Altitude of layer above ground
  ! ptm1, pum1, pvm1--R-: t, u and v
  ! pmea----input-R-Mean Orography (m)
  ! pstd----input-R-SSO standard deviation (m)
  ! psig----input-R-SSO slope
  ! pgam----input-R-SSO Anisotropy
  ! pthe----input-R-SSO Angle
  ! ppic----input-R-SSO Peacks elevation (m)
  ! pval----input-R-SSO Valleys elevation (m)
  ! plat----input-R-Latitude (degree)

  ! ==== outputs ===
  ! pulow, pvlow -output-R: Low-level wind

  ! pte -----output-R: T tendency
  ! pvom-----output-R: U tendency
  ! pvol-----output-R: V tendency


  ! Implicit Arguments:
  ! ===================

  ! klon-common-I: Number of points seen by the physics
  ! klev-common-I: Number of vertical layers


  ! ----------

  ! AUTHOR.
  ! -------
  ! F.LOTT  LMD 22/11/95

  USE dimphy
  IMPLICIT NONE


  include "YOMCST.h"
  include "YOEGWD.h"
  ! -----------------------------------------------------------------------

  ! *       0.1   ARGUMENTS
  ! ---------


  INTEGER nlon, nlev, kgwd
  REAL ptsphy
  REAL pte(nlon, nlev), pvol(nlon, nlev), pvom(nlon, nlev), pulow(nlon), &
    pvlow(nlon)
  REAL pum1(nlon, nlev), pvm1(nlon, nlev), ptm1(nlon, nlev), plat(nlon), &
    pmea(nlon), pstd(nlon), psig(nlon), pgam(nlon), pthe(nlon), ppic(nlon), &
    pval(nlon), pgeom1(nlon, nlev), papm1(nlon, nlev), paphm1(nlon, nlev+1)

  INTEGER kdx(nlon), ktest(nlon)
  ! -----------------------------------------------------------------------

  ! *       0.2   local arrays

  INTEGER jl, ilevh, jk
  REAL zhgeo, zdelp, zslow, zsqua, zscav, zbet
  ! ------------
  INTEGER iknub(klon), iknul(klon)
  LOGICAL ll1(klon, klev+1)

  REAL ztau(klon, klev+1), ztav(klon, klev+1), zrho(klon, klev+1)
  REAL zdudt(klon), zdvdt(klon)
  REAL zhcrit(klon, klev)

  LOGICAL lifthigh
  REAL zcons1, ztmst
  CHARACTER (LEN=20) :: modname = 'orolift_strato'
  CHARACTER (LEN=80) :: abort_message


  ! -----------------------------------------------------------------------

  ! *         1.1  initialisations
  ! ---------------

  lifthigh = .FALSE.

  IF (nlon/=klon .OR. nlev/=klev) THEN
    abort_message = 'pb dimension'
    CALL abort_physic(modname, abort_message, 1)
  END IF
  zcons1 = 1./rd
  ztmst = ptsphy

  DO jl = kidia, kfdia
    zrho(jl, klev+1) = 0.0
    pulow(jl) = 0.0
    pvlow(jl) = 0.0
    iknub(jl) = klev
    iknul(jl) = klev
    ilevh = klev/3
    ll1(jl, klev+1) = .FALSE.
    DO jk = 1, klev
      pvom(jl, jk) = 0.0
      pvol(jl, jk) = 0.0
      pte(jl, jk) = 0.0
    END DO
  END DO


  ! *         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
  ! *                 LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
  ! *                 THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.



  DO jk = klev, 1, -1
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        zhcrit(jl, jk) = amax1(ppic(jl)-pval(jl), 100.)
        zhgeo = pgeom1(jl, jk)/rg
        ll1(jl, jk) = (zhgeo>zhcrit(jl,jk))
        IF (ll1(jl,jk) .NEQV. ll1(jl,jk+1)) THEN
          iknub(jl) = jk
        END IF
      END IF
    END DO
  END DO


  DO jl = kidia, kfdia
    IF (ktest(jl)==1) THEN
      iknub(jl) = max(iknub(jl), klev/2)
      iknul(jl) = max(iknul(jl), 2*klev/3)
      IF (iknub(jl)>nktopg) iknub(jl) = nktopg
      IF (iknub(jl)==nktopg) iknul(jl) = klev
      IF (iknub(jl)==iknul(jl)) iknub(jl) = iknul(jl) - 1
    END IF
  END DO

  DO jk = klev, 2, -1
    DO jl = kidia, kfdia
      zrho(jl, jk) = 2.*paphm1(jl, jk)*zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
    END DO
  END DO
  ! print *,'  dans orolift: 223'

  ! ********************************************************************

  ! *     define low level flow
  ! -------------------
  DO jk = klev, 1, -1
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        IF (jk>=iknub(jl) .AND. jk<=iknul(jl)) THEN
          pulow(jl) = pulow(jl) + pum1(jl, jk)*(paphm1(jl,jk+1)-paphm1(jl,jk) &
            )
          pvlow(jl) = pvlow(jl) + pvm1(jl, jk)*(paphm1(jl,jk+1)-paphm1(jl,jk) &
            )
          zrho(jl, klev+1) = zrho(jl, klev+1) + zrho(jl, jk)*(paphm1(jl,jk+1) &
            -paphm1(jl,jk))
        END IF
      END IF
    END DO
  END DO
  DO jl = kidia, kfdia
    IF (ktest(jl)==1) THEN
      pulow(jl) = pulow(jl)/(paphm1(jl,iknul(jl)+1)-paphm1(jl,iknub(jl)))
      pvlow(jl) = pvlow(jl)/(paphm1(jl,iknul(jl)+1)-paphm1(jl,iknub(jl)))
      zrho(jl, klev+1) = zrho(jl, klev+1)/(paphm1(jl,iknul(jl)+1)-paphm1(jl, &
        iknub(jl)))
    END IF
  END DO

  ! ***********************************************************

  ! *         3.      COMPUTE MOUNTAIN LIFT


  DO jl = kidia, kfdia
    IF (ktest(jl)==1) THEN
      ztau(jl, klev+1) = -gklift*zrho(jl, klev+1)*2.*romega* & ! *
                                                               ! (2*pstd(jl)+pmea(jl))*
        2*pstd(jl)*sin(rpi/180.*plat(jl))*pvlow(jl)
      ztav(jl, klev+1) = gklift*zrho(jl, klev+1)*2.*romega* & ! *
                                                              ! (2*pstd(jl)+pmea(jl))*
        2*pstd(jl)*sin(rpi/180.*plat(jl))*pulow(jl)
    ELSE
      ztau(jl, klev+1) = 0.0
      ztav(jl, klev+1) = 0.0
    END IF
  END DO

  ! *         4.      COMPUTE LIFT PROFILE
  ! *                 --------------------



  DO jk = 1, klev
    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        ztau(jl, jk) = ztau(jl, klev+1)*paphm1(jl, jk)/paphm1(jl, klev+1)
        ztav(jl, jk) = ztav(jl, klev+1)*paphm1(jl, jk)/paphm1(jl, klev+1)
      ELSE
        ztau(jl, jk) = 0.0
        ztav(jl, jk) = 0.0
      END IF
    END DO
  END DO


  ! *         5.      COMPUTE TENDENCIES.
  ! *                 -------------------
  IF (lifthigh) THEN
    ! EXPLICIT SOLUTION AT ALL LEVELS

    DO jk = 1, klev
      DO jl = kidia, kfdia
        IF (ktest(jl)==1) THEN
          zdelp = paphm1(jl, jk+1) - paphm1(jl, jk)
          zdudt(jl) = -rg*(ztau(jl,jk+1)-ztau(jl,jk))/zdelp
          zdvdt(jl) = -rg*(ztav(jl,jk+1)-ztav(jl,jk))/zdelp
        END IF
      END DO
    END DO

    ! PROJECT PERPENDICULARLY TO U NOT TO DESTROY ENERGY

    DO jk = 1, klev
      DO jl = kidia, kfdia
        IF (ktest(jl)==1) THEN

          zslow = sqrt(pulow(jl)**2+pvlow(jl)**2)
          zsqua = amax1(sqrt(pum1(jl,jk)**2+pvm1(jl,jk)**2), gvsec)
          zscav = -zdudt(jl)*pvm1(jl, jk) + zdvdt(jl)*pum1(jl, jk)
          IF (zsqua>gvsec) THEN
            pvom(jl, jk) = -zscav*pvm1(jl, jk)/zsqua**2
            pvol(jl, jk) = zscav*pum1(jl, jk)/zsqua**2
          ELSE
            pvom(jl, jk) = 0.0
            pvol(jl, jk) = 0.0
          END IF
          zsqua = sqrt(pum1(jl,jk)**2+pum1(jl,jk)**2)
          IF (zsqua<zslow) THEN
            pvom(jl, jk) = zsqua/zslow*pvom(jl, jk)
            pvol(jl, jk) = zsqua/zslow*pvol(jl, jk)
          END IF

        END IF
      END DO
    END DO

    ! 6.  LOW LEVEL LIFT, SEMI IMPLICIT:
    ! ----------------------------------

  ELSE

    DO jl = kidia, kfdia
      IF (ktest(jl)==1) THEN
        DO jk = klev, iknub(jl), -1
          zbet = gklift*2.*romega*sin(rpi/180.*plat(jl))*ztmst* &
            (pgeom1(jl,iknub(jl)-1)-pgeom1(jl,jk))/ &
            (pgeom1(jl,iknub(jl)-1)-pgeom1(jl,klev))
          zdudt(jl) = -pum1(jl, jk)/ztmst/(1+zbet**2)
          zdvdt(jl) = -pvm1(jl, jk)/ztmst/(1+zbet**2)
          pvom(jl, jk) = zbet**2*zdudt(jl) - zbet*zdvdt(jl)
          pvol(jl, jk) = zbet*zdudt(jl) + zbet**2*zdvdt(jl)
        END DO
      END IF
    END DO

  END IF

  ! print *,' out orolift'

  RETURN
END SUBROUTINE orolift_strato
SUBROUTINE sugwd_strato(nlon, nlev, paprs, pplay)


  ! **** *SUGWD* INITIALIZE COMMON YOEGWD CONTROLLING GRAVITY WAVE DRAG

  ! PURPOSE.
  ! --------
  ! INITIALIZE YOEGWD, THE COMMON THAT CONTROLS THE
  ! GRAVITY WAVE DRAG PARAMETRIZATION.
  ! VERY IMPORTANT:
  ! ______________
  ! THIS ROUTINE SET_UP THE "TUNABLE PARAMETERS" OF THE
  ! VARIOUS SSO SCHEMES

  ! **   INTERFACE.
  ! ----------
  ! CALL *SUGWD* FROM *SUPHEC*
  ! -----        ------

  ! EXPLICIT ARGUMENTS :
  ! --------------------
  ! PAPRS,PPLAY : Pressure at semi and full model levels
  ! NLEV        : number of model levels
  ! NLON        : number of points treated in the physics

  ! IMPLICIT ARGUMENTS :
  ! --------------------
  ! COMMON YOEGWD
  ! -GFRCRIT-R:  Critical Non-dimensional mountain Height
  ! (HNC in (1),    LOTT 1999)
  ! -GKWAKE--R:  Bluff-body drag coefficient for low level wake
  ! (Cd in (2),     LOTT 1999)
  ! -GRCRIT--R:  Critical Richardson Number
  ! (Ric, End of first column p791 of LOTT 1999)
  ! -GKDRAG--R:  Gravity wave drag coefficient
  ! (G in (3),      LOTT 1999)
  ! -GKLIFT--R:  Mountain Lift coefficient
  ! (Cl in (4),     LOTT 1999)
  ! -GHMAX---R:  Not used
  ! -GRAHILO-R:  Set-up the trapped waves fraction
  ! (Beta , End of first column,  LOTT 1999)

  ! -GSIGCR--R:  Security value for blocked flow depth
  ! -NKTOPG--I:  Security value for blocked flow level
  ! -nstra----I:  An estimate to qualify the upper levels of
  ! the model where one wants to impose strees
  ! profiles
  ! -GSSECC--R:  Security min value for low-level B-V frequency
  ! -GTSEC---R:  Security min value for anisotropy and GW stress.
  ! -GVSEC---R:  Security min value for ulow


  ! METHOD.
  ! -------
  ! SEE DOCUMENTATION

  ! EXTERNALS.
  ! ----------
  ! NONE

  ! REFERENCE.
  ! ----------
  ! Lott, 1999: Alleviation of stationary biases in a GCM through...
  ! Monthly Weather Review, 127, pp 788-801.

  ! AUTHOR.
  ! -------
  ! FRANCOIS LOTT        *LMD*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 90-01-01 (MARTIN MILLER, ECMWF)
  ! LAST:  99-07-09     (FRANCOIS LOTT,LMD)
  ! ------------------------------------------------------------------
  USE dimphy
  USE mod_phys_lmdz_para
  USE mod_grid_phy_lmdz
  USE geometry_mod
  IMPLICIT NONE

  ! -----------------------------------------------------------------
  include "YOEGWD.h"
  ! ----------------------------------------------------------------

  ! ARGUMENTS
  INTEGER nlon, nlev
  REAL paprs(nlon, nlev+1)
  REAL pplay(nlon, nlev)

  INTEGER jk
  REAL zpr, ztop, zsigt, zpm1r
  INTEGER :: cell,ij,nstra_tmp,nktopg_tmp
  REAL :: current_dist, dist_min,dist_min_glo

  ! *       1.    SET THE VALUES OF THE PARAMETERS
  ! --------------------------------


  PRINT *, ' DANS SUGWD NLEV=', nlev
  ghmax = 10000.

  zpr = 100000.
  ZTOP=0.00005
  zsigt = 0.94
  ! old  ZPR=80000.
  ! old  ZSIGT=0.85


!ym Take the point at equator close to (0,0) coordinates.
  dist_min=360
  dist_min_glo=360.
  cell=-1
  DO ij=1,klon
    current_dist=sqrt(longitude_deg(ij)**2+latitude_deg(ij)**2)
    current_dist=current_dist*(1+(1e-10*ind_cell_glo(ij))/klon_glo) ! For point unicity
    IF (dist_min>current_dist) THEN
      dist_min=current_dist
      cell=ij    
    ENDIF  
  ENDDO
  
  !PRINT *, 'SUGWD distmin cell=', dist_min,cell
  CALL reduce_min(dist_min,dist_min_glo)
  CALL bcast(dist_min_glo)
  IF (dist_min/=dist_min_glo) cell=-1
!ym in future find the point at equator close to (0,0) coordinates.
  PRINT *, 'SUGWD distmin dist_min_glo cell=', dist_min,dist_min_glo,cell

  nktopg_tmp=nktopg
  nstra_tmp=nstra
  
  IF (cell/=-1) THEN

    !print*,'SUGWD shape ',shape(pplay),cell+1

    DO jk = 1, nlev
      !zpm1r = pplay(cell+1, jk)/paprs(cell+1, 1)
      zpm1r = pplay(cell, jk)/paprs(cell, 1)
      IF (zpm1r>=zsigt) THEN
        nktopg_tmp = jk
      END IF
      IF (zpm1r>=ztop) THEN
        nstra_tmp = jk
      END IF
    END DO
  ELSE
    nktopg_tmp=0
    nstra_tmp=0
  ENDIF
  
  CALL reduce_sum(nktopg_tmp,nktopg)
  CALL bcast(nktopg)
  CALL reduce_sum(nstra_tmp,nstra)
  CALL bcast(nstra)
  
  ! inversion car dans orodrag on compte les niveaux a l'envers
  nktopg = nlev - nktopg + 1
  nstra = nlev - nstra
  PRINT *, ' DANS SUGWD nktopg=', nktopg
  PRINT *, ' DANS SUGWD nstra=', nstra
  if (nstra == 0) call abort_physic("sugwd_strato", "no level in stratosphere", 1)

!  Valeurs lues dans les .def, ou attribues dans conf_phys
  !gkdrag = 0.2   
  !grahilo = 0.1
  !grcrit = 1.00
  !gfrcrit = 0.70
  !gkwake = 0.40
  !gklift = 0.25

  gsigcr = 0.80 ! Top of low level flow
  gvcrit = 0.1

  WRITE (UNIT=6, FMT='('' *** SSO essential constants ***'')')
  WRITE (UNIT=6, FMT='('' *** SPECIFIED IN SUGWD ***'')')
  WRITE (UNIT=6, FMT='('' Gravity wave ct '',E13.7,'' '')') gkdrag
  WRITE (UNIT=6, FMT='('' Trapped/total wave dag '',E13.7,'' '')') grahilo
  WRITE (UNIT=6, FMT='('' Critical Richardson   = '',E13.7,'' '')') grcrit
  WRITE (UNIT=6, FMT='('' Critical Froude'',e13.7)') gfrcrit
  WRITE (UNIT=6, FMT='('' Low level Wake bluff cte'',e13.7)') gkwake
  WRITE (UNIT=6, FMT='('' Low level lift  cte'',e13.7)') gklift

  ! ----------------------------------------------------------------

  ! *       2.    SET VALUES OF SECURITY PARAMETERS
  ! ---------------------------------


  gvsec = 0.10
  gssec = 0.0001

  gtsec = 0.00001

  RETURN
END SUBROUTINE sugwd_strato
