
! $Header$

SUBROUTINE diagphy(airephy, tit, iprt, tops, topl, sols, soll, sens, evap, &
    rain_fall, snow_fall, ts, d_etp_tot, d_qt_tot, d_ec_tot, fs_bound, &
    fq_bound)
  ! ======================================================================

  ! Purpose:
  ! Compute the thermal flux and the watter mass flux at the atmosphere
  ! boundaries. Print them and also the atmospheric enthalpy change and
  ! the  atmospheric mass change.

  ! Arguments:
  ! airephy-------input-R-  grid area
  ! tit---------input-A15- Comment to be added in PRINT (CHARACTER*15)
  ! iprt--------input-I-  PRINT level ( <=0 : no PRINT)
  ! tops(klon)--input-R-  SW rad. at TOA (W/m2), positive up.
  ! topl(klon)--input-R-  LW rad. at TOA (W/m2), positive down
  ! sols(klon)--input-R-  Net SW flux above surface (W/m2), positive up
  ! (i.e. -1 * flux absorbed by the surface)
  ! soll(klon)--input-R-  Net LW flux above surface (W/m2), positive up
  ! (i.e. flux emited - flux absorbed by the surface)
  ! sens(klon)--input-R-  Sensible Flux at surface  (W/m2), positive down
  ! evap(klon)--input-R-  Evaporation + sublimation watter vapour mass flux
  ! (kg/m2/s), positive up
  ! rain_fall(klon)
  ! --input-R- Liquid  watter mass flux (kg/m2/s), positive down
  ! snow_fall(klon)
  ! --input-R- Solid  watter mass flux (kg/m2/s), positive down
  ! ts(klon)----input-R- Surface temperature (K)
  ! d_etp_tot---input-R- Heat flux equivalent to atmospheric enthalpy
  ! change (W/m2)
  ! d_qt_tot----input-R- Mass flux equivalent to atmospheric watter mass
  ! change (kg/m2/s)
  ! d_ec_tot----input-R- Flux equivalent to atmospheric cinetic energy
  ! change (W/m2)

  ! fs_bound---output-R- Thermal flux at the atmosphere boundaries (W/m2)
  ! fq_bound---output-R- Watter mass flux at the atmosphere boundaries
  ! (kg/m2/s)

  ! J.L. Dufresne, July 2002
  ! Version prise sur
  ! ~rlmd833/LMDZOR_201102/modipsl/modeles/LMDZ.3.3/libf/phylmd
  ! le 25 Novembre 2002.
  ! ======================================================================

  USE dimphy
  IMPLICIT NONE

  include "YOMCST.h"
  include "YOETHF.h"

  ! Input variables
  REAL airephy(klon)
  CHARACTER *15 tit
  INTEGER iprt
  REAL tops(klon), topl(klon), sols(klon), soll(klon)
  REAL sens(klon), evap(klon), rain_fall(klon), snow_fall(klon)
  REAL ts(klon)
  REAL d_etp_tot, d_qt_tot, d_ec_tot
  ! Output variables
  REAL fs_bound, fq_bound

  ! Local variables
  REAL stops, stopl, ssols, ssoll
  REAL ssens, sfront, slat
  REAL airetot, zcpvap, zcwat, zcice
  REAL rain_fall_tot, snow_fall_tot, evap_tot

  INTEGER i

  INTEGER pas
  SAVE pas
  DATA pas/0/
  !$OMP THREADPRIVATE(pas)

  pas = pas + 1
  stops = 0.
  stopl = 0.
  ssols = 0.
  ssoll = 0.
  ssens = 0.
  sfront = 0.
  evap_tot = 0.
  rain_fall_tot = 0.
  snow_fall_tot = 0.
  airetot = 0.

  ! Pour les chaleur specifiques de la vapeur d'eau, de l'eau et de
  ! la glace, on travaille par difference a la chaleur specifique de l'
  ! air sec. En effet, comme on travaille a niveau de pression donne,
  ! toute variation de la masse d'un constituant est totalement
  ! compense par une variation de masse d'air.

  zcpvap = rcpv - rcpd
  zcwat = rcw - rcpd
  zcice = rcs - rcpd

  DO i = 1, klon
    stops = stops + tops(i)*airephy(i)
    stopl = stopl + topl(i)*airephy(i)
    ssols = ssols + sols(i)*airephy(i)
    ssoll = ssoll + soll(i)*airephy(i)
    ssens = ssens + sens(i)*airephy(i)
    sfront = sfront + (evap(i)*zcpvap-rain_fall(i)*zcwat-snow_fall(i)*zcice)* &
      ts(i)*airephy(i)
    evap_tot = evap_tot + evap(i)*airephy(i)
    rain_fall_tot = rain_fall_tot + rain_fall(i)*airephy(i)
    snow_fall_tot = snow_fall_tot + snow_fall(i)*airephy(i)
    airetot = airetot + airephy(i)
  END DO
  stops = stops/airetot
  stopl = stopl/airetot
  ssols = ssols/airetot
  ssoll = ssoll/airetot
  ssens = ssens/airetot
  sfront = sfront/airetot
  evap_tot = evap_tot/airetot
  rain_fall_tot = rain_fall_tot/airetot
  snow_fall_tot = snow_fall_tot/airetot

  slat = rlvtt*rain_fall_tot + rlstt*snow_fall_tot
  ! Heat flux at atm. boundaries
  fs_bound = stops - stopl - (ssols+ssoll) + ssens + sfront + slat
  ! Watter flux at atm. boundaries
  fq_bound = evap_tot - rain_fall_tot - snow_fall_tot

  IF (iprt>=1) WRITE (6, 6666) tit, pas, fs_bound, d_etp_tot, fq_bound, &
    d_qt_tot

  IF (iprt>=1) WRITE (6, 6668) tit, pas, d_etp_tot + d_ec_tot - fs_bound, &
    d_qt_tot - fq_bound

  IF (iprt>=2) WRITE (6, 6667) tit, pas, stops, stopl, ssols, ssoll, ssens, &
    slat, evap_tot, rain_fall_tot + snow_fall_tot

  RETURN

6666 FORMAT ('Phys. Flux Budget ', A15, 1I6, 2F8.2, 2(1PE13.5))
6667 FORMAT ('Phys. Boundary Flux ', A15, 1I6, 6F8.2, 2(1PE13.5))
6668 FORMAT ('Phys. Total Budget ', A15, 1I6, F8.2, 2(1PE13.5))

END SUBROUTINE diagphy

! ======================================================================
SUBROUTINE diagetpq(airephy, tit, iprt, idiag, idiag2, dtime, t, q, ql, qs, &
    u, v, paprs, pplay, d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
  ! ======================================================================

  ! Purpose:
  ! Calcul la difference d'enthalpie et de masse d'eau entre 2 appels,
  ! et calcul le flux de chaleur et le flux d'eau necessaire a ces
  ! changements. Ces valeurs sont moyennees sur la surface de tout
  ! le globe et sont exprime en W/2 et kg/s/m2
  ! Outil pour diagnostiquer la conservation de l'energie
  ! et de la masse dans la physique. Suppose que les niveau de
  ! pression entre couche ne varie pas entre 2 appels.

  ! Plusieurs de ces diagnostics peuvent etre fait en parallele: les
  ! bilans sont sauvegardes dans des tableaux indices. On parlera
  ! "d'indice de diagnostic"


  ! ======================================================================
  ! Arguments:
  ! airephy-------input-R-  grid area
  ! tit-----imput-A15- Comment added in PRINT (CHARACTER*15)
  ! iprt----input-I-  PRINT level ( <=1 : no PRINT)
  ! idiag---input-I- indice dans lequel sera range les nouveaux
  ! bilans d' entalpie et de masse
  ! idiag2--input-I-les nouveaux bilans d'entalpie et de masse
  ! sont compare au bilan de d'enthalpie de masse de
  ! l'indice numero idiag2
  ! Cas parriculier : si idiag2=0, pas de comparaison, on
  ! sort directement les bilans d'enthalpie et de masse
  ! dtime----input-R- time step (s)
  ! t--------input-R- temperature (K)
  ! q--------input-R- vapeur d'eau (kg/kg)
  ! ql-------input-R- liquid watter (kg/kg)
  ! qs-------input-R- solid watter (kg/kg)
  ! u--------input-R- vitesse u
  ! v--------input-R- vitesse v
  ! paprs----input-R- pression a intercouche (Pa)
  ! pplay----input-R- pression au milieu de couche (Pa)

  ! the following total value are computed by UNIT of earth surface

  ! d_h_vcol--output-R- Heat flux (W/m2) define as the Enthalpy
  ! change (J/m2) during one time step (dtime) for the whole
  ! atmosphere (air, watter vapour, liquid and solid)
  ! d_qt------output-R- total water mass flux (kg/m2/s) defined as the
  ! total watter (kg/m2) change during one time step (dtime),
  ! d_qw------output-R- same, for the watter vapour only (kg/m2/s)
  ! d_ql------output-R- same, for the liquid watter only (kg/m2/s)
  ! d_qs------output-R- same, for the solid watter only (kg/m2/s)
  ! d_ec------output-R- Cinetic Energy Budget (W/m2) for vertical air column

  ! other (COMMON...)
  ! RCPD, RCPV, ....

  ! J.L. Dufresne, July 2002
  ! ======================================================================

  USE dimphy
  IMPLICIT NONE

  include "YOMCST.h"
  include "YOETHF.h"

  ! Input variables
  REAL airephy(klon)
  CHARACTER *15 tit
  INTEGER iprt, idiag, idiag2
  REAL dtime
  REAL t(klon, klev), q(klon, klev), ql(klon, klev), qs(klon, klev)
  REAL u(klon, klev), v(klon, klev)
  REAL paprs(klon, klev+1), pplay(klon, klev)
  ! Output variables
  REAL d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec

  ! Local variables

  REAL h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot, h_qs_tot, qw_tot, ql_tot, &
    qs_tot, ec_tot
  ! h_vcol_tot--  total enthalpy of vertical air column
  ! (air with watter vapour, liquid and solid) (J/m2)
  ! h_dair_tot-- total enthalpy of dry air (J/m2)
  ! h_qw_tot----  total enthalpy of watter vapour (J/m2)
  ! h_ql_tot----  total enthalpy of liquid watter (J/m2)
  ! h_qs_tot----  total enthalpy of solid watter  (J/m2)
  ! qw_tot------  total mass of watter vapour (kg/m2)
  ! ql_tot------  total mass of liquid watter (kg/m2)
  ! qs_tot------  total mass of solid watter (kg/m2)
  ! ec_tot------  total cinetic energy (kg/m2)

  REAL zairm(klon, klev) ! layer air mass (kg/m2)
  REAL zqw_col(klon)
  REAL zql_col(klon)
  REAL zqs_col(klon)
  REAL zec_col(klon)
  REAL zh_dair_col(klon)
  REAL zh_qw_col(klon), zh_ql_col(klon), zh_qs_col(klon)

  REAL d_h_dair, d_h_qw, d_h_ql, d_h_qs

  REAL airetot, zcpvap, zcwat, zcice

  INTEGER i, k

  INTEGER ndiag ! max number of diagnostic in parallel
  PARAMETER (ndiag=10)
  INTEGER pas(ndiag)
  SAVE pas
  DATA pas/ndiag*0/
  !$OMP THREADPRIVATE(pas)

  REAL h_vcol_pre(ndiag), h_dair_pre(ndiag), h_qw_pre(ndiag), &
    h_ql_pre(ndiag), h_qs_pre(ndiag), qw_pre(ndiag), ql_pre(ndiag), &
    qs_pre(ndiag), ec_pre(ndiag)
  SAVE h_vcol_pre, h_dair_pre, h_qw_pre, h_ql_pre, h_qs_pre, qw_pre, ql_pre, &
    qs_pre, ec_pre
  !$OMP THREADPRIVATE(h_vcol_pre, h_dair_pre, h_qw_pre, h_ql_pre)
  !$OMP THREADPRIVATE(h_qs_pre, qw_pre, ql_pre, qs_pre , ec_pre)
  ! ======================================================================

  DO k = 1, klev
    DO i = 1, klon
      ! layer air mass
      zairm(i, k) = (paprs(i,k)-paprs(i,k+1))/rg
    END DO
  END DO

  ! Reset variables
  DO i = 1, klon
    zqw_col(i) = 0.
    zql_col(i) = 0.
    zqs_col(i) = 0.
    zec_col(i) = 0.
    zh_dair_col(i) = 0.
    zh_qw_col(i) = 0.
    zh_ql_col(i) = 0.
    zh_qs_col(i) = 0.
  END DO

  zcpvap = rcpv
  zcwat = rcw
  zcice = rcs

  ! Compute vertical sum for each atmospheric column
  ! ================================================
  DO k = 1, klev
    DO i = 1, klon
      ! Watter mass
      zqw_col(i) = zqw_col(i) + q(i, k)*zairm(i, k)
      zql_col(i) = zql_col(i) + ql(i, k)*zairm(i, k)
      zqs_col(i) = zqs_col(i) + qs(i, k)*zairm(i, k)
      ! Cinetic Energy
      zec_col(i) = zec_col(i) + 0.5*(u(i,k)**2+v(i,k)**2)*zairm(i, k)
      ! Air enthalpy
      zh_dair_col(i) = zh_dair_col(i) + rcpd*(1.-q(i,k)-ql(i,k)-qs(i,k))* &
        zairm(i, k)*t(i, k)
      zh_qw_col(i) = zh_qw_col(i) + zcpvap*q(i, k)*zairm(i, k)*t(i, k)
      zh_ql_col(i) = zh_ql_col(i) + zcwat*ql(i, k)*zairm(i, k)*t(i, k) - &
        rlvtt*ql(i, k)*zairm(i, k)
      zh_qs_col(i) = zh_qs_col(i) + zcice*qs(i, k)*zairm(i, k)*t(i, k) - &
        rlstt*qs(i, k)*zairm(i, k)

    END DO
  END DO

  ! Mean over the planete surface
  ! =============================
  qw_tot = 0.
  ql_tot = 0.
  qs_tot = 0.
  ec_tot = 0.
  h_vcol_tot = 0.
  h_dair_tot = 0.
  h_qw_tot = 0.
  h_ql_tot = 0.
  h_qs_tot = 0.
  airetot = 0.

  DO i = 1, klon
    qw_tot = qw_tot + zqw_col(i)*airephy(i)
    ql_tot = ql_tot + zql_col(i)*airephy(i)
    qs_tot = qs_tot + zqs_col(i)*airephy(i)
    ec_tot = ec_tot + zec_col(i)*airephy(i)
    h_dair_tot = h_dair_tot + zh_dair_col(i)*airephy(i)
    h_qw_tot = h_qw_tot + zh_qw_col(i)*airephy(i)
    h_ql_tot = h_ql_tot + zh_ql_col(i)*airephy(i)
    h_qs_tot = h_qs_tot + zh_qs_col(i)*airephy(i)
    airetot = airetot + airephy(i)
  END DO

  qw_tot = qw_tot/airetot
  ql_tot = ql_tot/airetot
  qs_tot = qs_tot/airetot
  ec_tot = ec_tot/airetot
  h_dair_tot = h_dair_tot/airetot
  h_qw_tot = h_qw_tot/airetot
  h_ql_tot = h_ql_tot/airetot
  h_qs_tot = h_qs_tot/airetot

  h_vcol_tot = h_dair_tot + h_qw_tot + h_ql_tot + h_qs_tot

  ! Compute the change of the atmospheric state compare to the one
  ! stored in "idiag2", and convert it in flux. THis computation
  ! is performed IF idiag2 /= 0 and IF it is not the first CALL
  ! for "idiag"
  ! ===================================

  IF ((idiag2>0) .AND. (pas(idiag2)/=0)) THEN
    d_h_vcol = (h_vcol_tot-h_vcol_pre(idiag2))/dtime
    d_h_dair = (h_dair_tot-h_dair_pre(idiag2))/dtime
    d_h_qw = (h_qw_tot-h_qw_pre(idiag2))/dtime
    d_h_ql = (h_ql_tot-h_ql_pre(idiag2))/dtime
    d_h_qs = (h_qs_tot-h_qs_pre(idiag2))/dtime
    d_qw = (qw_tot-qw_pre(idiag2))/dtime
    d_ql = (ql_tot-ql_pre(idiag2))/dtime
    d_qs = (qs_tot-qs_pre(idiag2))/dtime
    d_ec = (ec_tot-ec_pre(idiag2))/dtime
    d_qt = d_qw + d_ql + d_qs
  ELSE
    d_h_vcol = 0.
    d_h_dair = 0.
    d_h_qw = 0.
    d_h_ql = 0.
    d_h_qs = 0.
    d_qw = 0.
    d_ql = 0.
    d_qs = 0.
    d_ec = 0.
    d_qt = 0.
  END IF

  IF (iprt>=2) THEN
    WRITE (6, 9000) tit, pas(idiag), d_qt, d_qw, d_ql, d_qs
9000 FORMAT ('Phys. Watter Mass Budget (kg/m2/s)', A15, 1I6, 10(1PE14.6))
    WRITE (6, 9001) tit, pas(idiag), d_h_vcol
9001 FORMAT ('Phys. Enthalpy Budget (W/m2) ', A15, 1I6, 10(F8.2))
    WRITE (6, 9002) tit, pas(idiag), d_ec
9002 FORMAT ('Phys. Cinetic Energy Budget (W/m2) ', A15, 1I6, 10(F8.2))
  END IF

  ! Store the new atmospheric state in "idiag"

  pas(idiag) = pas(idiag) + 1
  h_vcol_pre(idiag) = h_vcol_tot
  h_dair_pre(idiag) = h_dair_tot
  h_qw_pre(idiag) = h_qw_tot
  h_ql_pre(idiag) = h_ql_tot
  h_qs_pre(idiag) = h_qs_tot
  qw_pre(idiag) = qw_tot
  ql_pre(idiag) = ql_tot
  qs_pre(idiag) = qs_tot
  ec_pre(idiag) = ec_tot

  RETURN
END SUBROUTINE diagetpq
