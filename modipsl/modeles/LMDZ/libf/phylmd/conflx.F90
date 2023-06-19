
! $Header$

SUBROUTINE conflx(dtime, pres_h, pres_f, t, q, con_t, con_q, pqhfl, w, d_t, &
    d_q, rain, snow, pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, kcbot, kctop, &
    kdtop, pmflxr, pmflxs)

  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19941014
  ! Objet: Schema flux de masse pour la convection
  ! (schema de Tiedtke avec qqs modifications mineures)
  ! Dec.97: Prise en compte des modifications introduites par
  ! Olivier Boucher et Alexandre Armengaud pour melange
  ! et lessivage des traceurs passifs.
  ! ======================================================================
  include "YOMCST.h"
  include "YOETHF.h"
  ! Entree:
  REAL dtime ! pas d'integration (s)
  REAL pres_h(klon, klev+1) ! pression half-level (Pa)
  REAL pres_f(klon, klev) ! pression full-level (Pa)
  REAL t(klon, klev) ! temperature (K)
  REAL q(klon, klev) ! humidite specifique (g/g)
  REAL w(klon, klev) ! vitesse verticale (Pa/s)
  REAL con_t(klon, klev) ! convergence de temperature (K/s)
  REAL con_q(klon, klev) ! convergence de l'eau vapeur (g/g/s)
  REAL pqhfl(klon) ! evaporation (negative vers haut) mm/s
  ! Sortie:
  REAL d_t(klon, klev) ! incrementation de temperature
  REAL d_q(klon, klev) ! incrementation d'humidite
  REAL pmfu(klon, klev) ! flux masse (kg/m2/s) panache ascendant
  REAL pmfd(klon, klev) ! flux masse (kg/m2/s) panache descendant
  REAL pen_u(klon, klev)
  REAL pen_d(klon, klev)
  REAL pde_u(klon, klev)
  REAL pde_d(klon, klev)
  REAL rain(klon) ! pluie (mm/s)
  REAL snow(klon) ! neige (mm/s)
  REAL pmflxr(klon, klev+1)
  REAL pmflxs(klon, klev+1)
  INTEGER kcbot(klon) ! niveau du bas de la convection
  INTEGER kctop(klon) ! niveau du haut de la convection
  INTEGER kdtop(klon) ! niveau du haut des downdrafts
  ! Local:
  REAL pt(klon, klev)
  REAL pq(klon, klev)
  REAL pqs(klon, klev)
  REAL pvervel(klon, klev)
  LOGICAL land(klon)

  REAL d_t_bis(klon, klev)
  REAL d_q_bis(klon, klev)
  REAL paprs(klon, klev+1)
  REAL paprsf(klon, klev)
  REAL zgeom(klon, klev)
  REAL zcvgq(klon, klev)
  REAL zcvgt(klon, klev)
  ! AA
  REAL zmfu(klon, klev)
  REAL zmfd(klon, klev)
  REAL zen_u(klon, klev)
  REAL zen_d(klon, klev)
  REAL zde_u(klon, klev)
  REAL zde_d(klon, klev)
  REAL zmflxr(klon, klev+1)
  REAL zmflxs(klon, klev+1)
  ! AA


  INTEGER i, k
  REAL zdelta, zqsat

  include "FCTTRE.h"

  ! initialiser les variables de sortie (pour securite)
  DO i = 1, klon
    rain(i) = 0.0
    snow(i) = 0.0
    kcbot(i) = 0
    kctop(i) = 0
    kdtop(i) = 0
  END DO
  DO k = 1, klev
    DO i = 1, klon
      d_t(i, k) = 0.0
      d_q(i, k) = 0.0
      pmfu(i, k) = 0.0
      pmfd(i, k) = 0.0
      pen_u(i, k) = 0.0
      pde_u(i, k) = 0.0
      pen_d(i, k) = 0.0
      pde_d(i, k) = 0.0
      zmfu(i, k) = 0.0
      zmfd(i, k) = 0.0
      zen_u(i, k) = 0.0
      zde_u(i, k) = 0.0
      zen_d(i, k) = 0.0
      zde_d(i, k) = 0.0
    END DO
  END DO
  DO k = 1, klev + 1
    DO i = 1, klon
      zmflxr(i, k) = 0.0
      zmflxs(i, k) = 0.0
    END DO
  END DO

  ! calculer la nature du sol (pour l'instant, ocean partout)
  DO i = 1, klon
    land(i) = .FALSE.
  END DO

  ! preparer les variables d'entree (attention: l'ordre des niveaux
  ! verticaux augmente du haut vers le bas)
  DO k = 1, klev
    DO i = 1, klon
      pt(i, k) = t(i, klev-k+1)
      pq(i, k) = q(i, klev-k+1)
      paprsf(i, k) = pres_f(i, klev-k+1)
      paprs(i, k) = pres_h(i, klev+1-k+1)
      pvervel(i, k) = w(i, klev+1-k)
      zcvgt(i, k) = con_t(i, klev-k+1)
      zcvgq(i, k) = con_q(i, klev-k+1)

      zdelta = max(0., sign(1.,rtt-pt(i,k)))
      zqsat = r2es*foeew(pt(i,k), zdelta)/paprsf(i, k)
      zqsat = min(0.5, zqsat)
      zqsat = zqsat/(1.-retv*zqsat)
      pqs(i, k) = zqsat
    END DO
  END DO
  DO i = 1, klon
    paprs(i, klev+1) = pres_h(i, 1)
    zgeom(i, klev) = rd*pt(i, klev)/(0.5*(paprs(i,klev+1)+paprsf(i, &
      klev)))*(paprs(i,klev+1)-paprsf(i,klev))
  END DO
  DO k = klev - 1, 1, -1
    DO i = 1, klon
      zgeom(i, k) = zgeom(i, k+1) + rd*0.5*(pt(i,k+1)+pt(i,k))/paprs(i, k+1)* &
        (paprsf(i,k+1)-paprsf(i,k))
    END DO
  END DO

  ! appeler la routine principale

  CALL flxmain(dtime, pt, pq, pqs, pqhfl, paprsf, paprs, zgeom, land, zcvgt, &
    zcvgq, pvervel, rain, snow, kcbot, kctop, kdtop, zmfu, zmfd, zen_u, &
    zde_u, zen_d, zde_d, d_t_bis, d_q_bis, zmflxr, zmflxs)

  ! AA--------------------------------------------------------
  ! AA rem : De la meme facon que l'on effectue le reindicage
  ! AA       pour la temperature t et le champ q
  ! AA       on reindice les flux necessaires a la convection
  ! AA       des traceurs
  ! AA--------------------------------------------------------
  DO k = 1, klev
    DO i = 1, klon
      d_q(i, klev+1-k) = dtime*d_q_bis(i, k)
      d_t(i, klev+1-k) = dtime*d_t_bis(i, k)
    END DO
  END DO

  DO i = 1, klon
    pmfu(i, 1) = 0.
    pmfd(i, 1) = 0.
    pen_d(i, 1) = 0.
    pde_d(i, 1) = 0.
  END DO

  DO k = 2, klev
    DO i = 1, klon
      pmfu(i, klev+2-k) = zmfu(i, k)
      pmfd(i, klev+2-k) = zmfd(i, k)
    END DO
  END DO

  DO k = 1, klev
    DO i = 1, klon
      pen_u(i, klev+1-k) = zen_u(i, k)
      pde_u(i, klev+1-k) = zde_u(i, k)
    END DO
  END DO

  DO k = 1, klev - 1
    DO i = 1, klon
      pen_d(i, klev+1-k) = -zen_d(i, k+1)
      pde_d(i, klev+1-k) = -zde_d(i, k+1)
    END DO
  END DO

  DO k = 1, klev + 1
    DO i = 1, klon
      pmflxr(i, klev+2-k) = zmflxr(i, k)
      pmflxs(i, klev+2-k) = zmflxs(i, k)
    END DO
  END DO

  RETURN
END SUBROUTINE conflx
! --------------------------------------------------------------------
SUBROUTINE flxmain(pdtime, pten, pqen, pqsen, pqhfl, pap, paph, pgeo, ldland, &
    ptte, pqte, pvervel, prsfc, pssfc, kcbot, kctop, kdtop, & ! *
                                                              ! ldcum, ktype,
    pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, dt_con, dq_con, pmflxr, pmflxs)
  USE dimphy
  IMPLICIT NONE
  ! ------------------------------------------------------------------
  include "YOMCST.h"
  include "YOETHF.h"
  include "YOECUMF.h"
  ! ----------------------------------------------------------------
  REAL pten(klon, klev), pqen(klon, klev), pqsen(klon, klev)
  REAL ptte(klon, klev)
  REAL pqte(klon, klev)
  REAL pvervel(klon, klev)
  REAL pgeo(klon, klev), pap(klon, klev), paph(klon, klev+1)
  REAL pqhfl(klon)

  REAL ptu(klon, klev), pqu(klon, klev), plu(klon, klev)
  REAL plude(klon, klev)
  REAL pmfu(klon, klev)
  REAL prsfc(klon), pssfc(klon)
  INTEGER kcbot(klon), kctop(klon), ktype(klon)
  LOGICAL ldland(klon), ldcum(klon)

  REAL ztenh(klon, klev), zqenh(klon, klev), zqsenh(klon, klev)
  REAL zgeoh(klon, klev)
  REAL zmfub(klon), zmfub1(klon)
  REAL zmfus(klon, klev), zmfuq(klon, klev), zmful(klon, klev)
  REAL zdmfup(klon, klev), zdpmel(klon, klev)
  REAL zentr(klon), zhcbase(klon)
  REAL zdqpbl(klon), zdqcv(klon), zdhpbl(klon)
  REAL zrfl(klon)
  REAL pmflxr(klon, klev+1)
  REAL pmflxs(klon, klev+1)
  INTEGER ilab(klon, klev), ictop0(klon)
  LOGICAL llo1
  REAL dt_con(klon, klev), dq_con(klon, klev)
  REAL zmfmax, zdh
  REAL pdtime, zqumqe, zdqmin, zalvdcp, zhsat, zzz
  REAL zhhat, zpbmpt, zgam, zeps, zfac
  INTEGER i, k, ikb, itopm2, kcum

  REAL pen_u(klon, klev), pde_u(klon, klev)
  REAL pen_d(klon, klev), pde_d(klon, klev)

  REAL ptd(klon, klev), pqd(klon, klev), pmfd(klon, klev)
  REAL zmfds(klon, klev), zmfdq(klon, klev), zdmfdp(klon, klev)
  INTEGER kdtop(klon)
  LOGICAL lddraf(klon)
  ! ---------------------------------------------------------------------
  LOGICAL firstcal
  SAVE firstcal
  DATA firstcal/.TRUE./
  !$OMP THREADPRIVATE(firstcal)
  ! ---------------------------------------------------------------------
  IF (firstcal) THEN
    CALL flxsetup
    firstcal = .FALSE.
  END IF
  ! ---------------------------------------------------------------------
  DO i = 1, klon
    ldcum(i) = .FALSE.
  END DO
  DO k = 1, klev
    DO i = 1, klon
      dt_con(i, k) = 0.0
      dq_con(i, k) = 0.0
    END DO
  END DO
  ! ----------------------------------------------------------------------
  ! initialiser les variables et faire l'interpolation verticale
  ! ----------------------------------------------------------------------
  CALL flxini(pten, pqen, pqsen, pgeo, paph, zgeoh, ztenh, zqenh, zqsenh, &
    ptu, pqu, ptd, pqd, pmfd, zmfds, zmfdq, zdmfdp, pmfu, zmfus, zmfuq, &
    zdmfup, zdpmel, plu, plude, ilab, pen_u, pde_u, pen_d, pde_d)
  ! ---------------------------------------------------------------------
  ! determiner les valeurs au niveau de base de la tour convective
  ! ---------------------------------------------------------------------
  CALL flxbase(ztenh, zqenh, zgeoh, paph, ptu, pqu, plu, ldcum, kcbot, ilab)
  ! ---------------------------------------------------------------------
  ! calculer la convergence totale de l'humidite et celle en provenance
  ! de la couche limite, plus precisement, la convergence integree entre
  ! le sol et la base de la convection. Cette derniere convergence est
  ! comparee avec l'evaporation obtenue dans la couche limite pour
  ! determiner le type de la convection
  ! ---------------------------------------------------------------------
  k = 1
  DO i = 1, klon
    zdqcv(i) = pqte(i, k)*(paph(i,k+1)-paph(i,k))
    zdhpbl(i) = 0.0
    zdqpbl(i) = 0.0
  END DO

  DO k = 2, klev
    DO i = 1, klon
      zdqcv(i) = zdqcv(i) + pqte(i, k)*(paph(i,k+1)-paph(i,k))
      IF (k>=kcbot(i)) THEN
        zdqpbl(i) = zdqpbl(i) + pqte(i, k)*(paph(i,k+1)-paph(i,k))
        zdhpbl(i) = zdhpbl(i) + (rcpd*ptte(i,k)+rlvtt*pqte(i,k))*(paph(i,k+1) &
          -paph(i,k))
      END IF
    END DO
  END DO

  DO i = 1, klon
    ktype(i) = 2
    IF (zdqcv(i)>max(0.,-1.5*pqhfl(i)*rg)) ktype(i) = 1
    ! cc         if (zdqcv(i).GT.MAX(0.,-1.1*pqhfl(i)*RG)) ktype(i) = 1
  END DO

  ! ---------------------------------------------------------------------
  ! determiner le flux de masse entrant a travers la base.
  ! on ignore, pour l'instant, l'effet du panache descendant
  ! ---------------------------------------------------------------------
  DO i = 1, klon
    ikb = kcbot(i)
    zqumqe = pqu(i, ikb) + plu(i, ikb) - zqenh(i, ikb)
    zdqmin = max(0.01*zqenh(i,ikb), 1.E-10)
    IF (zdqpbl(i)>0. .AND. zqumqe>zdqmin .AND. ldcum(i)) THEN
      zmfub(i) = zdqpbl(i)/(rg*max(zqumqe,zdqmin))
    ELSE
      zmfub(i) = 0.01
      ldcum(i) = .FALSE.
    END IF
    IF (ktype(i)==2) THEN
      zdh = rcpd*(ptu(i,ikb)-ztenh(i,ikb)) + rlvtt*zqumqe
      zdh = rg*max(zdh, 1.0E5*zdqmin)
      IF (zdhpbl(i)>0. .AND. ldcum(i)) zmfub(i) = zdhpbl(i)/zdh
    END IF
    zmfmax = (paph(i,ikb)-paph(i,ikb-1))/(rg*pdtime)
    zmfub(i) = min(zmfub(i), zmfmax)
    zentr(i) = entrscv
    IF (ktype(i)==1) zentr(i) = entrpen
  END DO
  ! -----------------------------------------------------------------------
  ! DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
  ! -----------------------------------------------------------------------
  ! (A) calculer d'abord la hauteur "theorique" de la tour convective sans
  ! considerer l'entrainement ni le detrainement du panache, sachant
  ! ces derniers peuvent abaisser la hauteur theorique.

  DO i = 1, klon
    ikb = kcbot(i)
    zhcbase(i) = rcpd*ptu(i, ikb) + zgeoh(i, ikb) + rlvtt*pqu(i, ikb)
    ictop0(i) = kcbot(i) - 1
  END DO

  zalvdcp = rlvtt/rcpd
  DO k = klev - 1, 3, -1
    DO i = 1, klon
      zhsat = rcpd*ztenh(i, k) + zgeoh(i, k) + rlvtt*zqsenh(i, k)
      zgam = r5les*zalvdcp*zqsenh(i, k)/((1.-retv*zqsenh(i,k))*(ztenh(i, &
        k)-r4les)**2)
      zzz = rcpd*ztenh(i, k)*0.608
      zhhat = zhsat - (zzz+zgam*zzz)/(1.+zgam*zzz/rlvtt)*max(zqsenh(i,k)- &
        zqenh(i,k), 0.)
      IF (k<ictop0(i) .AND. zhcbase(i)>zhhat) ictop0(i) = k
    END DO
  END DO

  ! (B) calculer le panache ascendant

  CALL flxasc(pdtime, ztenh, zqenh, pten, pqen, pqsen, pgeo, zgeoh, pap, &
    paph, pqte, pvervel, ldland, ldcum, ktype, ilab, ptu, pqu, plu, pmfu, &
    zmfub, zentr, zmfus, zmfuq, zmful, plude, zdmfup, kcbot, kctop, ictop0, &
    kcum, pen_u, pde_u)
  IF (kcum==0) GO TO 1000

  ! verifier l'epaisseur de la convection et changer eventuellement
  ! le taux d'entrainement/detrainement

  DO i = 1, klon
    zpbmpt = paph(i, kcbot(i)) - paph(i, kctop(i))
    IF (ldcum(i) .AND. ktype(i)==1 .AND. zpbmpt<2.E4) ktype(i) = 2
    IF (ldcum(i)) ictop0(i) = kctop(i)
    IF (ktype(i)==2) zentr(i) = entrscv
  END DO

  IF (lmfdd) THEN ! si l'on considere le panache descendant

    ! calculer la precipitation issue du panache ascendant pour
    ! determiner l'existence du panache descendant dans la convection
    DO i = 1, klon
      zrfl(i) = zdmfup(i, 1)
    END DO
    DO k = 2, klev
      DO i = 1, klon
        zrfl(i) = zrfl(i) + zdmfup(i, k)
      END DO
    END DO

    ! determiner le LFS (level of free sinking: niveau de plonge libre)
    CALL flxdlfs(ztenh, zqenh, zgeoh, paph, ptu, pqu, ldcum, kcbot, kctop, &
      zmfub, zrfl, ptd, pqd, pmfd, zmfds, zmfdq, zdmfdp, kdtop, lddraf)

    ! calculer le panache descendant
    CALL flxddraf(ztenh, zqenh, zgeoh, paph, zrfl, ptd, pqd, pmfd, zmfds, &
      zmfdq, zdmfdp, lddraf, pen_d, pde_d)

    ! calculer de nouveau le flux de masse entrant a travers la base
    ! de la convection, sachant qu'il a ete modifie par le panache
    ! descendant
    DO i = 1, klon
      IF (lddraf(i)) THEN
        ikb = kcbot(i)
        llo1 = pmfd(i, ikb) < 0.
        zeps = 0.
        IF (llo1) zeps = cmfdeps
        zqumqe = pqu(i, ikb) + plu(i, ikb) - zeps*pqd(i, ikb) - &
          (1.-zeps)*zqenh(i, ikb)
        zdqmin = max(0.01*zqenh(i,ikb), 1.E-10)
        zmfmax = (paph(i,ikb)-paph(i,ikb-1))/(rg*pdtime)
        IF (zdqpbl(i)>0. .AND. zqumqe>zdqmin .AND. ldcum(i) .AND. &
            zmfub(i)<zmfmax) THEN
          zmfub1(i) = zdqpbl(i)/(rg*max(zqumqe,zdqmin))
        ELSE
          zmfub1(i) = zmfub(i)
        END IF
        IF (ktype(i)==2) THEN
          zdh = rcpd*(ptu(i,ikb)-zeps*ptd(i,ikb)-(1.-zeps)*ztenh(i,ikb)) + &
            rlvtt*zqumqe
          zdh = rg*max(zdh, 1.0E5*zdqmin)
          IF (zdhpbl(i)>0. .AND. ldcum(i)) zmfub1(i) = zdhpbl(i)/zdh
        END IF
        IF (.NOT. ((ktype(i)==1 .OR. ktype(i)==2) .AND. abs(zmfub1(i)-zmfub(i &
          ))<0.2*zmfub(i))) zmfub1(i) = zmfub(i)
      END IF
    END DO
    DO k = 1, klev
      DO i = 1, klon
        IF (lddraf(i)) THEN
          zfac = zmfub1(i)/max(zmfub(i), 1.E-10)
          pmfd(i, k) = pmfd(i, k)*zfac
          zmfds(i, k) = zmfds(i, k)*zfac
          zmfdq(i, k) = zmfdq(i, k)*zfac
          zdmfdp(i, k) = zdmfdp(i, k)*zfac
          pen_d(i, k) = pen_d(i, k)*zfac
          pde_d(i, k) = pde_d(i, k)*zfac
        END IF
      END DO
    END DO
    DO i = 1, klon
      IF (lddraf(i)) zmfub(i) = zmfub1(i)
    END DO

  END IF ! fin de test sur lmfdd

  ! -----------------------------------------------------------------------
  ! calculer de nouveau le panache ascendant
  ! -----------------------------------------------------------------------
  CALL flxasc(pdtime, ztenh, zqenh, pten, pqen, pqsen, pgeo, zgeoh, pap, &
    paph, pqte, pvervel, ldland, ldcum, ktype, ilab, ptu, pqu, plu, pmfu, &
    zmfub, zentr, zmfus, zmfuq, zmful, plude, zdmfup, kcbot, kctop, ictop0, &
    kcum, pen_u, pde_u)

  ! -----------------------------------------------------------------------
  ! determiner les flux convectifs en forme finale, ainsi que
  ! la quantite des precipitations
  ! -----------------------------------------------------------------------
  CALL flxflux(pdtime, pqen, pqsen, ztenh, zqenh, pap, paph, ldland, zgeoh, &
    kcbot, kctop, lddraf, kdtop, ktype, ldcum, pmfu, pmfd, zmfus, zmfds, &
    zmfuq, zmfdq, zmful, plude, zdmfup, zdmfdp, pten, prsfc, pssfc, zdpmel, &
    itopm2, pmflxr, pmflxs)

  ! ----------------------------------------------------------------------
  ! calculer les tendances pour T et Q
  ! ----------------------------------------------------------------------
  CALL flxdtdq(pdtime, itopm2, paph, ldcum, pten, zmfus, zmfds, zmfuq, zmfdq, &
    zmful, zdmfup, zdmfdp, zdpmel, dt_con, dq_con)

1000 CONTINUE
  RETURN
END SUBROUTINE flxmain
SUBROUTINE flxini(pten, pqen, pqsen, pgeo, paph, pgeoh, ptenh, pqenh, pqsenh, &
    ptu, pqu, ptd, pqd, pmfd, pmfds, pmfdq, pdmfdp, pmfu, pmfus, pmfuq, &
    pdmfup, pdpmel, plu, plude, klab, pen_u, pde_u, pen_d, pde_d)
  USE dimphy
  IMPLICIT NONE
  ! ----------------------------------------------------------------------
  ! THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
  ! TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
  ! AND INITIALIZES VALUES FOR UPDRAFTS
  ! ----------------------------------------------------------------------
  include "YOMCST.h"
  include "YOETHF.h"

  REAL pten(klon, klev) ! temperature (environnement)
  REAL pqen(klon, klev) ! humidite (environnement)
  REAL pqsen(klon, klev) ! humidite saturante (environnement)
  REAL pgeo(klon, klev) ! geopotentiel (g * metre)
  REAL pgeoh(klon, klev) ! geopotentiel aux demi-niveaux
  REAL paph(klon, klev+1) ! pression aux demi-niveaux
  REAL ptenh(klon, klev) ! temperature aux demi-niveaux
  REAL pqenh(klon, klev) ! humidite aux demi-niveaux
  REAL pqsenh(klon, klev) ! humidite saturante aux demi-niveaux

  REAL ptu(klon, klev) ! temperature du panache ascendant (p-a)
  REAL pqu(klon, klev) ! humidite du p-a
  REAL plu(klon, klev) ! eau liquide du p-a
  REAL pmfu(klon, klev) ! flux de masse du p-a
  REAL pmfus(klon, klev) ! flux de l'energie seche dans le p-a
  REAL pmfuq(klon, klev) ! flux de l'humidite dans le p-a
  REAL pdmfup(klon, klev) ! quantite de l'eau precipitee dans p-a
  REAL plude(klon, klev) ! quantite de l'eau liquide jetee du
  ! p-a a l'environnement
  REAL pdpmel(klon, klev) ! quantite de neige fondue

  REAL ptd(klon, klev) ! temperature du panache descendant (p-d)
  REAL pqd(klon, klev) ! humidite du p-d
  REAL pmfd(klon, klev) ! flux de masse du p-d
  REAL pmfds(klon, klev) ! flux de l'energie seche dans le p-d
  REAL pmfdq(klon, klev) ! flux de l'humidite dans le p-d
  REAL pdmfdp(klon, klev) ! quantite de precipitation dans p-d

  REAL pen_u(klon, klev) ! quantite de masse entrainee pour p-a
  REAL pde_u(klon, klev) ! quantite de masse detrainee pour p-a
  REAL pen_d(klon, klev) ! quantite de masse entrainee pour p-d
  REAL pde_d(klon, klev) ! quantite de masse detrainee pour p-d

  INTEGER klab(klon, klev)
  LOGICAL llflag(klon)
  INTEGER k, i, icall
  REAL zzs
  ! ----------------------------------------------------------------------
  ! SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
  ! ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
  ! ----------------------------------------------------------------------
  DO k = 2, klev

    DO i = 1, klon
      pgeoh(i, k) = pgeo(i, k) + (pgeo(i,k-1)-pgeo(i,k))*0.5
      ptenh(i, k) = (max(rcpd*pten(i,k-1)+pgeo(i,k-1),rcpd*pten(i,k)+pgeo(i, &
        k))-pgeoh(i,k))/rcpd
      pqsenh(i, k) = pqsen(i, k-1)
      llflag(i) = .TRUE.
    END DO

    icall = 0
    CALL flxadjtq(paph(1,k), ptenh(1,k), pqsenh(1,k), llflag, icall)

    DO i = 1, klon
      pqenh(i, k) = min(pqen(i,k-1), pqsen(i,k-1)) + &
        (pqsenh(i,k)-pqsen(i,k-1))
      pqenh(i, k) = max(pqenh(i,k), 0.)
    END DO

  END DO

  DO i = 1, klon
    ptenh(i, klev) = (rcpd*pten(i,klev)+pgeo(i,klev)-pgeoh(i,klev))/rcpd
    pqenh(i, klev) = pqen(i, klev)
    ptenh(i, 1) = pten(i, 1)
    pqenh(i, 1) = pqen(i, 1)
    pgeoh(i, 1) = pgeo(i, 1)
  END DO

  DO k = klev - 1, 2, -1
    DO i = 1, klon
      zzs = max(rcpd*ptenh(i,k)+pgeoh(i,k), rcpd*ptenh(i,k+1)+pgeoh(i,k+1))
      ptenh(i, k) = (zzs-pgeoh(i,k))/rcpd
    END DO
  END DO

  ! -----------------------------------------------------------------------
  ! INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
  ! -----------------------------------------------------------------------
  DO k = 1, klev
    DO i = 1, klon
      ptu(i, k) = ptenh(i, k)
      pqu(i, k) = pqenh(i, k)
      plu(i, k) = 0.
      pmfu(i, k) = 0.
      pmfus(i, k) = 0.
      pmfuq(i, k) = 0.
      pdmfup(i, k) = 0.
      pdpmel(i, k) = 0.
      plude(i, k) = 0.

      klab(i, k) = 0

      ptd(i, k) = ptenh(i, k)
      pqd(i, k) = pqenh(i, k)
      pmfd(i, k) = 0.0
      pmfds(i, k) = 0.0
      pmfdq(i, k) = 0.0
      pdmfdp(i, k) = 0.0

      pen_u(i, k) = 0.0
      pde_u(i, k) = 0.0
      pen_d(i, k) = 0.0
      pde_d(i, k) = 0.0
    END DO
  END DO

  RETURN
END SUBROUTINE flxini
SUBROUTINE flxbase(ptenh, pqenh, pgeoh, paph, ptu, pqu, plu, ldcum, kcbot, &
    klab)
  USE dimphy
  IMPLICIT NONE
  ! ----------------------------------------------------------------------
  ! THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)

  ! INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
  ! IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
  ! klab=1 FOR SUBCLOUD LEVELS
  ! klab=2 FOR CONDENSATION LEVEL

  ! LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
  ! (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
  ! ----------------------------------------------------------------------
  include "YOMCST.h"
  include "YOETHF.h"
  ! ----------------------------------------------------------------
  REAL ptenh(klon, klev), pqenh(klon, klev)
  REAL pgeoh(klon, klev), paph(klon, klev+1)

  REAL ptu(klon, klev), pqu(klon, klev), plu(klon, klev)
  INTEGER klab(klon, klev), kcbot(klon)

  LOGICAL llflag(klon), ldcum(klon)
  INTEGER i, k, icall, is
  REAL zbuo, zqold(klon)
  ! ----------------------------------------------------------------------
  ! INITIALIZE VALUES AT LIFTING LEVEL
  ! ----------------------------------------------------------------------
  DO i = 1, klon
    klab(i, klev) = 1
    kcbot(i) = klev - 1
    ldcum(i) = .FALSE.
  END DO
  ! ----------------------------------------------------------------------
  ! DO ASCENT IN SUBCLOUD LAYER,
  ! CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
  ! ADJUST T,Q AND L ACCORDINGLY
  ! CHECK FOR BUOYANCY AND SET FLAGS
  ! ----------------------------------------------------------------------
  DO k = klev - 1, 2, -1

    is = 0
    DO i = 1, klon
      IF (klab(i,k+1)==1) is = is + 1
      llflag(i) = .FALSE.
      IF (klab(i,k+1)==1) llflag(i) = .TRUE.
    END DO
    IF (is==0) GO TO 290

    DO i = 1, klon
      IF (llflag(i)) THEN
        pqu(i, k) = pqu(i, k+1)
        ptu(i, k) = ptu(i, k+1) + (pgeoh(i,k+1)-pgeoh(i,k))/rcpd
        zbuo = ptu(i, k)*(1.+retv*pqu(i,k)) - ptenh(i, k)*(1.+retv*pqenh(i,k) &
          ) + 0.5
        IF (zbuo>0.) klab(i, k) = 1
        zqold(i) = pqu(i, k)
      END IF
    END DO

    icall = 1
    CALL flxadjtq(paph(1,k), ptu(1,k), pqu(1,k), llflag, icall)

    DO i = 1, klon
      IF (llflag(i) .AND. pqu(i,k)/=zqold(i)) THEN
        klab(i, k) = 2
        plu(i, k) = plu(i, k) + zqold(i) - pqu(i, k)
        zbuo = ptu(i, k)*(1.+retv*pqu(i,k)) - ptenh(i, k)*(1.+retv*pqenh(i,k) &
          ) + 0.5
        IF (zbuo>0.) kcbot(i) = k
        IF (zbuo>0.) ldcum(i) = .TRUE.
      END IF
    END DO

290 END DO

  RETURN
END SUBROUTINE flxbase
SUBROUTINE flxasc(pdtime, ptenh, pqenh, pten, pqen, pqsen, pgeo, pgeoh, pap, &
    paph, pqte, pvervel, ldland, ldcum, ktype, klab, ptu, pqu, plu, pmfu, &
    pmfub, pentr, pmfus, pmfuq, pmful, plude, pdmfup, kcbot, kctop, kctop0, &
    kcum, pen_u, pde_u)
  USE dimphy
  IMPLICIT NONE
  ! ----------------------------------------------------------------------
  ! THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
  ! FOR CUMULUS PARAMETERIZATION
  ! ----------------------------------------------------------------------
  include "YOMCST.h"
  include "YOETHF.h"
  include "YOECUMF.h"

  REAL pdtime
  REAL pten(klon, klev), ptenh(klon, klev)
  REAL pqen(klon, klev), pqenh(klon, klev), pqsen(klon, klev)
  REAL pgeo(klon, klev), pgeoh(klon, klev)
  REAL pap(klon, klev), paph(klon, klev+1)
  REAL pqte(klon, klev)
  REAL pvervel(klon, klev) ! vitesse verticale en Pa/s

  REAL pmfub(klon), pentr(klon)
  REAL ptu(klon, klev), pqu(klon, klev), plu(klon, klev)
  REAL plude(klon, klev)
  REAL pmfu(klon, klev), pmfus(klon, klev)
  REAL pmfuq(klon, klev), pmful(klon, klev)
  REAL pdmfup(klon, klev)
  INTEGER ktype(klon), klab(klon, klev), kcbot(klon), kctop(klon)
  INTEGER kctop0(klon)
  LOGICAL ldland(klon), ldcum(klon)

  REAL pen_u(klon, klev), pde_u(klon, klev)
  REAL zqold(klon)
  REAL zdland(klon)
  LOGICAL llflag(klon)
  INTEGER k, i, is, icall, kcum
  REAL ztglace, zdphi, zqeen, zseen, zscde, zqude
  REAL zmfusk, zmfuqk, zmfulk, zbuo, zdnoprc, zprcon, zlnew

  REAL zpbot(klon), zptop(klon), zrho(klon)
  REAL zdprho, zentr, zpmid, zmftest, zmfmax
  LOGICAL llo1, llo2

  REAL zwmax(klon), zzzmb
  INTEGER klwmin(klon) ! level of maximum vertical velocity
  REAL fact
  ! ----------------------------------------------------------------------
  ztglace = rtt - 13.

  ! Chercher le niveau ou la vitesse verticale est maximale:
  DO i = 1, klon
    klwmin(i) = klev
    zwmax(i) = 0.0
  END DO
  DO k = klev, 3, -1
    DO i = 1, klon
      IF (pvervel(i,k)<zwmax(i)) THEN
        zwmax(i) = pvervel(i, k)
        klwmin(i) = k
      END IF
    END DO
  END DO
  ! ----------------------------------------------------------------------
  ! SET DEFAULT VALUES
  ! ----------------------------------------------------------------------
  DO i = 1, klon
    IF (.NOT. ldcum(i)) ktype(i) = 0
  END DO

  DO k = 1, klev
    DO i = 1, klon
      plu(i, k) = 0.
      pmfu(i, k) = 0.
      pmfus(i, k) = 0.
      pmfuq(i, k) = 0.
      pmful(i, k) = 0.
      plude(i, k) = 0.
      pdmfup(i, k) = 0.
      IF (.NOT. ldcum(i) .OR. ktype(i)==3) klab(i, k) = 0
      IF (.NOT. ldcum(i) .AND. paph(i,k)<4.E4) kctop0(i) = k
    END DO
  END DO

  DO i = 1, klon
    IF (ldland(i)) THEN
      zdland(i) = 3.0E4
      zdphi = pgeoh(i, kctop0(i)) - pgeoh(i, kcbot(i))
      IF (ptu(i,kctop0(i))>=ztglace) zdland(i) = zdphi
      zdland(i) = max(3.0E4, zdland(i))
      zdland(i) = min(5.0E4, zdland(i))
    END IF
  END DO

  ! Initialiser les valeurs au niveau d'ascendance

  DO i = 1, klon
    kctop(i) = klev - 1
    IF (.NOT. ldcum(i)) THEN
      kcbot(i) = klev - 1
      pmfub(i) = 0.
      pqu(i, klev) = 0.
    END IF
    pmfu(i, klev) = pmfub(i)
    pmfus(i, klev) = pmfub(i)*(rcpd*ptu(i,klev)+pgeoh(i,klev))
    pmfuq(i, klev) = pmfub(i)*pqu(i, klev)
  END DO

  DO i = 1, klon
    ldcum(i) = .FALSE.
  END DO
  ! ----------------------------------------------------------------------
  ! DO ASCENT: SUBCLOUD LAYER (klab=1) ,CLOUDS (klab=2)
  ! BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
  ! BY ADJUSTING T,Q AND L ACCORDINGLY IN *flxadjtq*,
  ! THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
  ! ----------------------------------------------------------------------
  DO k = klev - 1, 3, -1

    IF (lmfmid .AND. k<klev-1) THEN
      DO i = 1, klon
        IF (.NOT. ldcum(i) .AND. klab(i,k+1)==0 .AND. &
            pqen(i,k)>0.9*pqsen(i,k) .AND. pap(i,k)/paph(i,klev+1)>0.4) THEN
          ptu(i, k+1) = pten(i, k) + (pgeo(i,k)-pgeoh(i,k+1))/rcpd
          pqu(i, k+1) = pqen(i, k)
          plu(i, k+1) = 0.0
          zzzmb = max(cmfcmin, -pvervel(i,k)/rg)
          zmfmax = (paph(i,k)-paph(i,k-1))/(rg*pdtime)
          pmfub(i) = min(zzzmb, zmfmax)
          pmfu(i, k+1) = pmfub(i)
          pmfus(i, k+1) = pmfub(i)*(rcpd*ptu(i,k+1)+pgeoh(i,k+1))
          pmfuq(i, k+1) = pmfub(i)*pqu(i, k+1)
          pmful(i, k+1) = 0.0
          pdmfup(i, k+1) = 0.0
          kcbot(i) = k
          klab(i, k+1) = 1
          ktype(i) = 3
          pentr(i) = entrmid
        END IF
      END DO
    END IF

    is = 0
    DO i = 1, klon
      is = is + klab(i, k+1)
      IF (klab(i,k+1)==0) klab(i, k) = 0
      llflag(i) = .FALSE.
      IF (klab(i,k+1)>0) llflag(i) = .TRUE.
    END DO
    IF (is==0) GO TO 480

    ! calculer le taux d'entrainement et de detrainement

    DO i = 1, klon
      pen_u(i, k) = 0.0
      pde_u(i, k) = 0.0
      zrho(i) = paph(i, k+1)/(rd*ptenh(i,k+1))
      zpbot(i) = paph(i, kcbot(i))
      zptop(i) = paph(i, kctop0(i))
    END DO

    DO i = 1, klon
      IF (ldcum(i)) THEN
        zdprho = (paph(i,k+1)-paph(i,k))/(rg*zrho(i))
        zentr = pentr(i)*pmfu(i, k+1)*zdprho
        llo1 = k < kcbot(i)
        IF (llo1) pde_u(i, k) = zentr
        zpmid = 0.5*(zpbot(i)+zptop(i))
        llo2 = llo1 .AND. ktype(i) == 2 .AND. (zpbot(i)-paph(i,k)<0.2E5 .OR. &
          paph(i,k)>zpmid)
        IF (llo2) pen_u(i, k) = zentr
        llo2 = llo1 .AND. (ktype(i)==1 .OR. ktype(i)==3) .AND. &
          (k>=max(klwmin(i),kctop0(i)+2) .OR. pap(i,k)>zpmid)
        IF (llo2) pen_u(i, k) = zentr
        llo1 = pen_u(i, k) > 0. .AND. (ktype(i)==1 .OR. ktype(i)==2)
        IF (llo1) THEN
          fact = 1. + 3.*(1.-min(1.,(zpbot(i)-pap(i,k))/1.5E4))
          zentr = zentr*fact
          pen_u(i, k) = pen_u(i, k)*fact
          pde_u(i, k) = pde_u(i, k)*fact
        END IF
        IF (llo2 .AND. pqenh(i,k+1)>1.E-5) pen_u(i, k) = zentr + &
          max(pqte(i,k), 0.)/pqenh(i, k+1)*zrho(i)*zdprho
      END IF
    END DO

    ! ----------------------------------------------------------------------
    ! DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
    ! ----------------------------------------------------------------------

    DO i = 1, klon
      IF (llflag(i)) THEN
        IF (k<kcbot(i)) THEN
          zmftest = pmfu(i, k+1) + pen_u(i, k) - pde_u(i, k)
          zmfmax = min(zmftest, (paph(i,k)-paph(i,k-1))/(rg*pdtime))
          pen_u(i, k) = max(pen_u(i,k)-max(0.0,zmftest-zmfmax), 0.0)
        END IF
        pde_u(i, k) = min(pde_u(i,k), 0.75*pmfu(i,k+1))
        ! calculer le flux de masse du niveau k a partir de celui du k+1
        pmfu(i, k) = pmfu(i, k+1) + pen_u(i, k) - pde_u(i, k)
        ! calculer les valeurs Su, Qu et l du niveau k dans le panache
        ! montant
        zqeen = pqenh(i, k+1)*pen_u(i, k)
        zseen = (rcpd*ptenh(i,k+1)+pgeoh(i,k+1))*pen_u(i, k)
        zscde = (rcpd*ptu(i,k+1)+pgeoh(i,k+1))*pde_u(i, k)
        zqude = pqu(i, k+1)*pde_u(i, k)
        plude(i, k) = plu(i, k+1)*pde_u(i, k)
        zmfusk = pmfus(i, k+1) + zseen - zscde
        zmfuqk = pmfuq(i, k+1) + zqeen - zqude
        zmfulk = pmful(i, k+1) - plude(i, k)
        plu(i, k) = zmfulk*(1./max(cmfcmin,pmfu(i,k)))
        pqu(i, k) = zmfuqk*(1./max(cmfcmin,pmfu(i,k)))
        ptu(i, k) = (zmfusk*(1./max(cmfcmin,pmfu(i,k)))-pgeoh(i,k))/rcpd
        ptu(i, k) = max(100., ptu(i,k))
        ptu(i, k) = min(400., ptu(i,k))
        zqold(i) = pqu(i, k)
      ELSE
        zqold(i) = 0.0
      END IF
    END DO

    ! ----------------------------------------------------------------------
    ! DO CORRECTIONS FOR MOIST ASCENT BY ADJUSTING T,Q AND L
    ! ----------------------------------------------------------------------

    icall = 1
    CALL flxadjtq(paph(1,k), ptu(1,k), pqu(1,k), llflag, icall)

    DO i = 1, klon
      IF (llflag(i) .AND. pqu(i,k)/=zqold(i)) THEN
        klab(i, k) = 2
        plu(i, k) = plu(i, k) + zqold(i) - pqu(i, k)
        zbuo = ptu(i, k)*(1.+retv*pqu(i,k)) - ptenh(i, k)*(1.+retv*pqenh(i,k) &
          )
        IF (klab(i,k+1)==1) zbuo = zbuo + 0.5
        IF (zbuo>0. .AND. pmfu(i,k)>=0.1*pmfub(i)) THEN
          kctop(i) = k
          ldcum(i) = .TRUE.
          zdnoprc = 1.5E4
          IF (ldland(i)) zdnoprc = zdland(i)
          zprcon = cprcon
          IF ((zpbot(i)-paph(i,k))<zdnoprc) zprcon = 0.0
          zlnew = plu(i, k)/(1.+zprcon*(pgeoh(i,k)-pgeoh(i,k+1)))
          pdmfup(i, k) = max(0., (plu(i,k)-zlnew)*pmfu(i,k))
          plu(i, k) = zlnew
        ELSE
          klab(i, k) = 0
          pmfu(i, k) = 0.
        END IF
      END IF
    END DO
    DO i = 1, klon
      IF (llflag(i)) THEN
        pmful(i, k) = plu(i, k)*pmfu(i, k)
        pmfus(i, k) = (rcpd*ptu(i,k)+pgeoh(i,k))*pmfu(i, k)
        pmfuq(i, k) = pqu(i, k)*pmfu(i, k)
      END IF
    END DO

480 END DO
  ! ----------------------------------------------------------------------
  ! DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
  ! (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
  ! AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
  ! FROM PREVIOUS CALCULATIONS ABOVE)
  ! ----------------------------------------------------------------------
  DO i = 1, klon
    IF (kctop(i)==klev-1) ldcum(i) = .FALSE.
    kcbot(i) = max(kcbot(i), kctop(i))
  END DO

  ldcum(1) = ldcum(1)

  is = 0
  DO i = 1, klon
    IF (ldcum(i)) is = is + 1
  END DO
  kcum = is
  IF (is==0) GO TO 800

  DO i = 1, klon
    IF (ldcum(i)) THEN
      k = kctop(i) - 1
      pde_u(i, k) = (1.-cmfctop)*pmfu(i, k+1)
      plude(i, k) = pde_u(i, k)*plu(i, k+1)
      pmfu(i, k) = pmfu(i, k+1) - pde_u(i, k)
      zlnew = plu(i, k)
      pdmfup(i, k) = max(0., (plu(i,k)-zlnew)*pmfu(i,k))
      plu(i, k) = zlnew
      pmfus(i, k) = (rcpd*ptu(i,k)+pgeoh(i,k))*pmfu(i, k)
      pmfuq(i, k) = pqu(i, k)*pmfu(i, k)
      pmful(i, k) = plu(i, k)*pmfu(i, k)
      plude(i, k-1) = pmful(i, k)
    END IF
  END DO

800 CONTINUE
  RETURN
END SUBROUTINE flxasc
SUBROUTINE flxflux(pdtime, pqen, pqsen, ptenh, pqenh, pap, paph, ldland, &
    pgeoh, kcbot, kctop, lddraf, kdtop, ktype, ldcum, pmfu, pmfd, pmfus, &
    pmfds, pmfuq, pmfdq, pmful, plude, pdmfup, pdmfdp, pten, prfl, psfl, &
    pdpmel, ktopm2, pmflxr, pmflxs)
  USE dimphy
  USE print_control_mod, ONLY: prt_level
  IMPLICIT NONE
  ! ----------------------------------------------------------------------
  ! THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
  ! FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
  ! ----------------------------------------------------------------------
  include "YOMCST.h"
  include "YOETHF.h"
  include "YOECUMF.h"

  REAL cevapcu(klon, klev)
  ! -----------------------------------------------------------------
  REAL pqen(klon, klev), pqenh(klon, klev), pqsen(klon, klev)
  REAL pten(klon, klev), ptenh(klon, klev)
  REAL paph(klon, klev+1), pgeoh(klon, klev)

  REAL pap(klon, klev)
  REAL ztmsmlt, zdelta, zqsat

  REAL pmfu(klon, klev), pmfus(klon, klev)
  REAL pmfd(klon, klev), pmfds(klon, klev)
  REAL pmfuq(klon, klev), pmful(klon, klev)
  REAL pmfdq(klon, klev)
  REAL plude(klon, klev)
  REAL pdmfup(klon, klev), pdpmel(klon, klev)
  ! jq The variable maxpdmfdp(klon) has been introduced by Olivier Boucher
  ! jq 14/11/00 to fix the problem with the negative precipitation.
  REAL pdmfdp(klon, klev), maxpdmfdp(klon, klev)
  REAL prfl(klon), psfl(klon)
  REAL pmflxr(klon, klev+1), pmflxs(klon, klev+1)
  INTEGER kcbot(klon), kctop(klon), ktype(klon)
  LOGICAL ldland(klon), ldcum(klon)
  INTEGER k, kp, i
  REAL zcons1, zcons2, zcucov, ztmelp2
  REAL pdtime, zdp, zzp, zfac, zsnmlt, zrfl, zrnew
  REAL zrmin, zrfln, zdrfl
  REAL zpds, zpdr, zdenom
  INTEGER ktopm2, itop, ikb

  LOGICAL lddraf(klon)
  INTEGER kdtop(klon)

  include "FCTTRE.h"

  DO k = 1, klev
    DO i = 1, klon
      cevapcu(i, k) = 1.93E-6*261.*sqrt(1.E3/(38.3*0.293)*sqrt(0.5*(paph(i,k) &
        +paph(i,k+1))/paph(i,klev+1)))*0.5/rg
    END DO
  END DO

  ! SPECIFY CONSTANTS

  zcons1 = rcpd/(rlmlt*rg*pdtime)
  zcons2 = 1./(rg*pdtime)
  zcucov = 0.05
  ztmelp2 = rtt + 2.

  ! DETERMINE FINAL CONVECTIVE FLUXES

  itop = klev
  DO i = 1, klon
    itop = min(itop, kctop(i))
    IF (.NOT. ldcum(i) .OR. kdtop(i)<kctop(i)) lddraf(i) = .FALSE.
    IF (.NOT. ldcum(i)) ktype(i) = 0
  END DO

  ktopm2 = itop - 2
  DO k = ktopm2, klev
    DO i = 1, klon
      IF (ldcum(i) .AND. k>=kctop(i)-1) THEN
        pmfus(i, k) = pmfus(i, k) - pmfu(i, k)*(rcpd*ptenh(i,k)+pgeoh(i,k))
        pmfuq(i, k) = pmfuq(i, k) - pmfu(i, k)*pqenh(i, k)
        zdp = 1.5E4
        IF (ldland(i)) zdp = 3.E4

        ! l'eau liquide detrainee est precipitee quand certaines
        ! conditions sont reunies (sinon, elle est consideree
        ! evaporee dans l'environnement)

        IF (paph(i,kcbot(i))-paph(i,kctop(i))>=zdp .AND. pqen(i,k-1)>0.8* &
          pqsen(i,k-1)) pdmfup(i, k-1) = pdmfup(i, k-1) + plude(i, k-1)

        IF (lddraf(i) .AND. k>=kdtop(i)) THEN
          pmfds(i, k) = pmfds(i, k) - pmfd(i, k)*(rcpd*ptenh(i,k)+pgeoh(i,k))
          pmfdq(i, k) = pmfdq(i, k) - pmfd(i, k)*pqenh(i, k)
        ELSE
          pmfd(i, k) = 0.
          pmfds(i, k) = 0.
          pmfdq(i, k) = 0.
          pdmfdp(i, k-1) = 0.
        END IF
      ELSE
        pmfu(i, k) = 0.
        pmfus(i, k) = 0.
        pmfuq(i, k) = 0.
        pmful(i, k) = 0.
        pdmfup(i, k-1) = 0.
        plude(i, k-1) = 0.
        pmfd(i, k) = 0.
        pmfds(i, k) = 0.
        pmfdq(i, k) = 0.
        pdmfdp(i, k-1) = 0.
      END IF
    END DO
  END DO

  DO k = ktopm2, klev
    DO i = 1, klon
      IF (ldcum(i) .AND. k>kcbot(i)) THEN
        ikb = kcbot(i)
        zzp = ((paph(i,klev+1)-paph(i,k))/(paph(i,klev+1)-paph(i,ikb)))
        IF (ktype(i)==3) zzp = zzp**2
        pmfu(i, k) = pmfu(i, ikb)*zzp
        pmfus(i, k) = pmfus(i, ikb)*zzp
        pmfuq(i, k) = pmfuq(i, ikb)*zzp
        pmful(i, k) = pmful(i, ikb)*zzp
      END IF
    END DO
  END DO

  ! CALCULATE RAIN/SNOW FALL RATES
  ! CALCULATE MELTING OF SNOW
  ! CALCULATE EVAPORATION OF PRECIP

  DO k = 1, klev + 1
    DO i = 1, klon
      pmflxr(i, k) = 0.0
      pmflxs(i, k) = 0.0
    END DO
  END DO
  DO k = ktopm2, klev
    DO i = 1, klon
      IF (ldcum(i)) THEN
        IF (pmflxs(i,k)>0.0 .AND. pten(i,k)>ztmelp2) THEN
          zfac = zcons1*(paph(i,k+1)-paph(i,k))
          zsnmlt = min(pmflxs(i,k), zfac*(pten(i,k)-ztmelp2))
          pdpmel(i, k) = zsnmlt
          ztmsmlt = pten(i, k) - zsnmlt/zfac
          zdelta = max(0., sign(1.,rtt-ztmsmlt))
          zqsat = r2es*foeew(ztmsmlt, zdelta)/pap(i, k)
          zqsat = min(0.5, zqsat)
          zqsat = zqsat/(1.-retv*zqsat)
          pqsen(i, k) = zqsat
        END IF
        IF (pten(i,k)>rtt) THEN
          pmflxr(i, k+1) = pmflxr(i, k) + pdmfup(i, k) + pdmfdp(i, k) + &
            pdpmel(i, k)
          pmflxs(i, k+1) = pmflxs(i, k) - pdpmel(i, k)
        ELSE
          pmflxs(i, k+1) = pmflxs(i, k) + pdmfup(i, k) + pdmfdp(i, k)
          pmflxr(i, k+1) = pmflxr(i, k)
        END IF
        ! si la precipitation est negative, on ajuste le plux du
        ! panache descendant pour eliminer la negativite
        IF ((pmflxr(i,k+1)+pmflxs(i,k+1))<0.0) THEN
          pdmfdp(i, k) = -pmflxr(i, k) - pmflxs(i, k) - pdmfup(i, k)
          pmflxr(i, k+1) = 0.0
          pmflxs(i, k+1) = 0.0
          pdpmel(i, k) = 0.0
        END IF
      END IF
    END DO
  END DO

  ! jq The new variable is initialized here.
  ! jq It contains the humidity which is fed to the downdraft
  ! jq by evaporation of precipitation in the column below the base
  ! jq of convection.
  ! jq
  ! jq In the former version, this term has been subtracted from precip
  ! jq as well as the evaporation.
  ! jq
  DO k = 1, klev
    DO i = 1, klon
      maxpdmfdp(i, k) = 0.0
    END DO
  END DO
  DO k = 1, klev
    DO kp = k, klev
      DO i = 1, klon
        maxpdmfdp(i, k) = maxpdmfdp(i, k) + pdmfdp(i, kp)
      END DO
    END DO
  END DO
  ! jq End of initialization

  DO k = ktopm2, klev
    DO i = 1, klon
      IF (ldcum(i) .AND. k>=kcbot(i)) THEN
        zrfl = pmflxr(i, k) + pmflxs(i, k)
        IF (zrfl>1.0E-20) THEN
          zrnew = (max(0.,sqrt(zrfl/zcucov)-cevapcu(i, &
            k)*(paph(i,k+1)-paph(i,k))*max(0.,pqsen(i,k)-pqen(i,k))))**2* &
            zcucov
          zrmin = zrfl - zcucov*max(0., 0.8*pqsen(i,k)-pqen(i,k))*zcons2*( &
            paph(i,k+1)-paph(i,k))
          zrnew = max(zrnew, zrmin)
          zrfln = max(zrnew, 0.)
          zdrfl = min(0., zrfln-zrfl)
          ! jq At least the amount of precipiation needed to feed the
          ! downdraft
          ! jq with humidity below the base of convection has to be left and
          ! can't
          ! jq be evaporated (surely the evaporation can't be positive):
          zdrfl = max(zdrfl, min(-pmflxr(i,k)-pmflxs(i,k)-maxpdmfdp(i, &
            k),0.0))
          ! jq End of insertion

          zdenom = 1.0/max(1.0E-20, pmflxr(i,k)+pmflxs(i,k))
          IF (pten(i,k)>rtt) THEN
            zpdr = pdmfdp(i, k)
            zpds = 0.0
          ELSE
            zpdr = 0.0
            zpds = pdmfdp(i, k)
          END IF
          pmflxr(i, k+1) = pmflxr(i, k) + zpdr + pdpmel(i, k) + &
            zdrfl*pmflxr(i, k)*zdenom
          pmflxs(i, k+1) = pmflxs(i, k) + zpds - pdpmel(i, k) + &
            zdrfl*pmflxs(i, k)*zdenom
          pdmfup(i, k) = pdmfup(i, k) + zdrfl
        ELSE
          pmflxr(i, k+1) = 0.0
          pmflxs(i, k+1) = 0.0
          pdmfdp(i, k) = 0.0
          pdpmel(i, k) = 0.0
        END IF
        IF (pmflxr(i,k)+pmflxs(i,k)<-1.E-26 .AND. prt_level>=1) WRITE (*, *) &
          'precip. < 1e-16 ', pmflxr(i, k) + pmflxs(i, k)
      END IF
    END DO
  END DO

  DO i = 1, klon
    prfl(i) = pmflxr(i, klev+1)
    psfl(i) = pmflxs(i, klev+1)
  END DO

  RETURN
END SUBROUTINE flxflux
SUBROUTINE flxdtdq(pdtime, ktopm2, paph, ldcum, pten, pmfus, pmfds, pmfuq, &
    pmfdq, pmful, pdmfup, pdmfdp, pdpmel, dt_con, dq_con)
  USE dimphy
  IMPLICIT NONE
  ! ----------------------------------------------------------------------
  ! calculer les tendances T et Q
  ! ----------------------------------------------------------------------
  include "YOMCST.h"
  include "YOETHF.h"
  include "YOECUMF.h"
  ! -----------------------------------------------------------------
  LOGICAL llo1

  REAL pten(klon, klev), paph(klon, klev+1)
  REAL pmfus(klon, klev), pmfuq(klon, klev), pmful(klon, klev)
  REAL pmfds(klon, klev), pmfdq(klon, klev)
  REAL pdmfup(klon, klev)
  REAL pdmfdp(klon, klev)
  REAL pdpmel(klon, klev)
  LOGICAL ldcum(klon)
  REAL dt_con(klon, klev), dq_con(klon, klev)

  INTEGER ktopm2
  REAL pdtime

  INTEGER i, k
  REAL zalv, zdtdt, zdqdt

  DO k = ktopm2, klev - 1
    DO i = 1, klon
      IF (ldcum(i)) THEN
        llo1 = (pten(i,k)-rtt) > 0.
        zalv = rlstt
        IF (llo1) zalv = rlvtt
        zdtdt = rg/(paph(i,k+1)-paph(i,k))/rcpd*(pmfus(i,k+1)-pmfus(i,k)+ &
          pmfds(i,k+1)-pmfds(i,k)-rlmlt*pdpmel(i,k)-zalv*(pmful(i, &
          k+1)-pmful(i,k)-pdmfup(i,k)-pdmfdp(i,k)))
        dt_con(i, k) = zdtdt
        zdqdt = rg/(paph(i,k+1)-paph(i,k))*(pmfuq(i,k+1)-pmfuq(i,k)+pmfdq(i,k &
          +1)-pmfdq(i,k)+pmful(i,k+1)-pmful(i,k)-pdmfup(i,k)-pdmfdp(i,k))
        dq_con(i, k) = zdqdt
      END IF
    END DO
  END DO

  k = klev
  DO i = 1, klon
    IF (ldcum(i)) THEN
      llo1 = (pten(i,k)-rtt) > 0.
      zalv = rlstt
      IF (llo1) zalv = rlvtt
      zdtdt = -rg/(paph(i,k+1)-paph(i,k))/rcpd*(pmfus(i,k)+pmfds(i,k)+rlmlt* &
        pdpmel(i,k)-zalv*(pmful(i,k)+pdmfup(i,k)+pdmfdp(i,k)))
      dt_con(i, k) = zdtdt
      zdqdt = -rg/(paph(i,k+1)-paph(i,k))*(pmfuq(i,k)+pmfdq(i,k)+pmful(i,k)+ &
        pdmfup(i,k)+pdmfdp(i,k))
      dq_con(i, k) = zdqdt
    END IF
  END DO

  RETURN
END SUBROUTINE flxdtdq
SUBROUTINE flxdlfs(ptenh, pqenh, pgeoh, paph, ptu, pqu, ldcum, kcbot, kctop, &
    pmfub, prfl, ptd, pqd, pmfd, pmfds, pmfdq, pdmfdp, kdtop, lddraf)
  USE dimphy
  IMPLICIT NONE

  ! ----------------------------------------------------------------------
  ! THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
  ! CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES

  ! TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
  ! FOR MASSFLUX CUMULUS PARAMETERIZATION

  ! INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
  ! AND UPDRAFT VALUES T,Q,U AND V AND ALSO
  ! CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
  ! IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.

  ! CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
  ! MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
  ! ----------------------------------------------------------------------
  include "YOMCST.h"
  include "YOETHF.h"
  include "YOECUMF.h"

  REAL ptenh(klon, klev)
  REAL pqenh(klon, klev)
  REAL pgeoh(klon, klev), paph(klon, klev+1)
  REAL ptu(klon, klev), pqu(klon, klev)
  REAL pmfub(klon)
  REAL prfl(klon)

  REAL ptd(klon, klev), pqd(klon, klev)
  REAL pmfd(klon, klev), pmfds(klon, klev), pmfdq(klon, klev)
  REAL pdmfdp(klon, klev)
  INTEGER kcbot(klon), kctop(klon), kdtop(klon)
  LOGICAL ldcum(klon), lddraf(klon)

  REAL ztenwb(klon, klev), zqenwb(klon, klev), zcond(klon)
  REAL zttest, zqtest, zbuo, zmftop
  LOGICAL llo2(klon)
  INTEGER i, k, is, icall
  ! ----------------------------------------------------------------------
  DO i = 1, klon
    lddraf(i) = .FALSE.
    kdtop(i) = klev + 1
  END DO

  ! ----------------------------------------------------------------------
  ! DETERMINE LEVEL OF FREE SINKING BY
  ! DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS

  ! FOR EVERY POINT AND PROCEED AS FOLLOWS:
  ! (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
  ! (2) DO MIXING WITH CUMULUS CLOUD AIR
  ! (3) CHECK FOR NEGATIVE BUOYANCY

  ! THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
  ! OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
  ! TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
  ! EVAPORATION OF RAIN AND CLOUD WATER)
  ! ----------------------------------------------------------------------

  DO k = 3, klev - 3

    is = 0
    DO i = 1, klon
      ztenwb(i, k) = ptenh(i, k)
      zqenwb(i, k) = pqenh(i, k)
      llo2(i) = ldcum(i) .AND. prfl(i) > 0. .AND. .NOT. lddraf(i) .AND. &
        (k<kcbot(i) .AND. k>kctop(i))
      IF (llo2(i)) is = is + 1
    END DO
    IF (is==0) GO TO 290

    icall = 2
    CALL flxadjtq(paph(1,k), ztenwb(1,k), zqenwb(1,k), llo2, icall)

    ! ----------------------------------------------------------------------
    ! DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
    ! AND CHECK FOR NEGATIVE BUOYANCY.
    ! THEN SET VALUES FOR DOWNDRAFT AT LFS.
    ! ----------------------------------------------------------------------
    DO i = 1, klon
      IF (llo2(i)) THEN
        zttest = 0.5*(ptu(i,k)+ztenwb(i,k))
        zqtest = 0.5*(pqu(i,k)+zqenwb(i,k))
        zbuo = zttest*(1.+retv*zqtest) - ptenh(i, k)*(1.+retv*pqenh(i,k))
        zcond(i) = pqenh(i, k) - zqenwb(i, k)
        zmftop = -cmfdeps*pmfub(i)
        IF (zbuo<0. .AND. prfl(i)>10.*zmftop*zcond(i)) THEN
          kdtop(i) = k
          lddraf(i) = .TRUE.
          ptd(i, k) = zttest
          pqd(i, k) = zqtest
          pmfd(i, k) = zmftop
          pmfds(i, k) = pmfd(i, k)*(rcpd*ptd(i,k)+pgeoh(i,k))
          pmfdq(i, k) = pmfd(i, k)*pqd(i, k)
          pdmfdp(i, k-1) = -0.5*pmfd(i, k)*zcond(i)
          prfl(i) = prfl(i) + pdmfdp(i, k-1)
        END IF
      END IF
    END DO

290 END DO

  RETURN
END SUBROUTINE flxdlfs
SUBROUTINE flxddraf(ptenh, pqenh, pgeoh, paph, prfl, ptd, pqd, pmfd, pmfds, &
    pmfdq, pdmfdp, lddraf, pen_d, pde_d)
  USE dimphy
  IMPLICIT NONE

  ! ----------------------------------------------------------------------
  ! THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT

  ! TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
  ! (I.E. T,Q,U AND V AND FLUXES)

  ! INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
  ! IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
  ! AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS

  ! CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
  ! A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
  ! B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.

  ! ----------------------------------------------------------------------
  include "YOMCST.h"
  include "YOETHF.h"
  include "YOECUMF.h"

  REAL ptenh(klon, klev), pqenh(klon, klev)
  REAL pgeoh(klon, klev), paph(klon, klev+1)

  REAL ptd(klon, klev), pqd(klon, klev)
  REAL pmfd(klon, klev), pmfds(klon, klev), pmfdq(klon, klev)
  REAL pdmfdp(klon, klev)
  REAL prfl(klon)
  LOGICAL lddraf(klon)

  REAL pen_d(klon, klev), pde_d(klon, klev), zcond(klon)
  LOGICAL llo2(klon), llo1
  INTEGER i, k, is, icall, itopde
  REAL zentr, zseen, zqeen, zsdde, zqdde, zmfdsk, zmfdqk, zdmfdp
  REAL zbuo
  ! ----------------------------------------------------------------------
  ! CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
  ! (A) CALCULATING ENTRAINMENT RATES, ASSUMING
  ! LINEAR DECREASE OF MASSFLUX IN PBL
  ! (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
  ! AND MOISTENING IS CALCULATED IN *flxadjtq*
  ! (C) CHECKING FOR NEGATIVE BUOYANCY AND
  ! SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES

  DO k = 3, klev

    is = 0
    DO i = 1, klon
      llo2(i) = lddraf(i) .AND. pmfd(i, k-1) < 0.
      IF (llo2(i)) is = is + 1
    END DO
    IF (is==0) GO TO 180

    DO i = 1, klon
      IF (llo2(i)) THEN
        zentr = entrdd*pmfd(i, k-1)*rd*ptenh(i, k-1)/(rg*paph(i,k-1))* &
          (paph(i,k)-paph(i,k-1))
        pen_d(i, k) = zentr
        pde_d(i, k) = zentr
      END IF
    END DO

    itopde = klev - 2
    IF (k>itopde) THEN
      DO i = 1, klon
        IF (llo2(i)) THEN
          pen_d(i, k) = 0.
          pde_d(i, k) = pmfd(i, itopde)*(paph(i,k)-paph(i,k-1))/ &
            (paph(i,klev+1)-paph(i,itopde))
        END IF
      END DO
    END IF

    DO i = 1, klon
      IF (llo2(i)) THEN
        pmfd(i, k) = pmfd(i, k-1) + pen_d(i, k) - pde_d(i, k)
        zseen = (rcpd*ptenh(i,k-1)+pgeoh(i,k-1))*pen_d(i, k)
        zqeen = pqenh(i, k-1)*pen_d(i, k)
        zsdde = (rcpd*ptd(i,k-1)+pgeoh(i,k-1))*pde_d(i, k)
        zqdde = pqd(i, k-1)*pde_d(i, k)
        zmfdsk = pmfds(i, k-1) + zseen - zsdde
        zmfdqk = pmfdq(i, k-1) + zqeen - zqdde
        pqd(i, k) = zmfdqk*(1./min(-cmfcmin,pmfd(i,k)))
        ptd(i, k) = (zmfdsk*(1./min(-cmfcmin,pmfd(i,k)))-pgeoh(i,k))/rcpd
        ptd(i, k) = min(400., ptd(i,k))
        ptd(i, k) = max(100., ptd(i,k))
        zcond(i) = pqd(i, k)
      END IF
    END DO

    icall = 2
    CALL flxadjtq(paph(1,k), ptd(1,k), pqd(1,k), llo2, icall)

    DO i = 1, klon
      IF (llo2(i)) THEN
        zcond(i) = zcond(i) - pqd(i, k)
        zbuo = ptd(i, k)*(1.+retv*pqd(i,k)) - ptenh(i, k)*(1.+retv*pqenh(i,k) &
          )
        llo1 = zbuo < 0. .AND. (prfl(i)-pmfd(i,k)*zcond(i)>0.)
        IF (.NOT. llo1) pmfd(i, k) = 0.0
        pmfds(i, k) = (rcpd*ptd(i,k)+pgeoh(i,k))*pmfd(i, k)
        pmfdq(i, k) = pqd(i, k)*pmfd(i, k)
        zdmfdp = -pmfd(i, k)*zcond(i)
        pdmfdp(i, k-1) = zdmfdp
        prfl(i) = prfl(i) + zdmfdp
      END IF
    END DO

180 END DO
  RETURN
END SUBROUTINE flxddraf
SUBROUTINE flxadjtq(pp, pt, pq, ldflag, kcall)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Objet: ajustement entre T et Q
  ! ======================================================================
  ! NOTE: INPUT PARAMETER kcall DEFINES CALCULATION AS
  ! kcall=0    ENV. T AND QS IN*CUINI*
  ! kcall=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
  ! kcall=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)

  include "YOMCST.h"

  REAL pt(klon), pq(klon), pp(klon)
  LOGICAL ldflag(klon)
  INTEGER kcall

  REAL zcond(klon), zcond1
  REAL z5alvcp, z5alscp, zalvdcp, zalsdcp
  REAL zdelta, zcvm5, zldcp, zqsat, zcor
  INTEGER is, i
  include "YOETHF.h"
  include "FCTTRE.h"

  z5alvcp = r5les*rlvtt/rcpd
  z5alscp = r5ies*rlstt/rcpd
  zalvdcp = rlvtt/rcpd
  zalsdcp = rlstt/rcpd


  DO i = 1, klon
    zcond(i) = 0.0
  END DO

  DO i = 1, klon
    IF (ldflag(i)) THEN
      zdelta = max(0., sign(1.,rtt-pt(i)))
      zcvm5 = z5alvcp*(1.-zdelta) + zdelta*z5alscp
      zldcp = zalvdcp*(1.-zdelta) + zdelta*zalsdcp
      zqsat = r2es*foeew(pt(i), zdelta)/pp(i)
      zqsat = min(0.5, zqsat)
      zcor = 1./(1.-retv*zqsat)
      zqsat = zqsat*zcor
      zcond(i) = (pq(i)-zqsat)/(1.+foede(pt(i),zdelta,zcvm5,zqsat,zcor))
      IF (kcall==1) zcond(i) = max(zcond(i), 0.)
      IF (kcall==2) zcond(i) = min(zcond(i), 0.)
      pt(i) = pt(i) + zldcp*zcond(i)
      pq(i) = pq(i) - zcond(i)
    END IF
  END DO

  is = 0
  DO i = 1, klon
    IF (zcond(i)/=0.) is = is + 1
  END DO
  IF (is==0) GO TO 230

  DO i = 1, klon
    IF (ldflag(i) .AND. zcond(i)/=0.) THEN
      zdelta = max(0., sign(1.,rtt-pt(i)))
      zcvm5 = z5alvcp*(1.-zdelta) + zdelta*z5alscp
      zldcp = zalvdcp*(1.-zdelta) + zdelta*zalsdcp
      zqsat = r2es*foeew(pt(i), zdelta)/pp(i)
      zqsat = min(0.5, zqsat)
      zcor = 1./(1.-retv*zqsat)
      zqsat = zqsat*zcor
      zcond1 = (pq(i)-zqsat)/(1.+foede(pt(i),zdelta,zcvm5,zqsat,zcor))
      pt(i) = pt(i) + zldcp*zcond1
      pq(i) = pq(i) - zcond1
    END IF
  END DO

230 CONTINUE
  RETURN
END SUBROUTINE flxadjtq
SUBROUTINE flxsetup
  IMPLICIT NONE

  ! THIS ROUTINE DEFINES DISPOSABLE PARAMETERS FOR MASSFLUX SCHEME

  include "YOECUMF.h"

  entrpen = 1.0E-4 ! ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
  entrscv = 3.0E-4 ! ENTRAINMENT RATE FOR SHALLOW CONVECTION
  entrmid = 1.0E-4 ! ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
  entrdd = 2.0E-4 ! ENTRAINMENT RATE FOR DOWNDRAFTS
  cmfctop = 0.33 ! RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUO LEVEL
  cmfcmax = 1.0 ! MAXIMUM MASSFLUX VALUE ALLOWED FOR UPDRAFTS ETC
  cmfcmin = 1.E-10 ! MINIMUM MASSFLUX VALUE (FOR SAFETY)
  cmfdeps = 0.3 ! FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
  cprcon = 2.0E-4 ! CONVERSION FROM CLOUD WATER TO RAIN
  rhcdd = 1. ! RELATIVE SATURATION IN DOWNDRAFRS (NO LONGER USED)
  ! (FORMULATION IMPLIES SATURATION)
  lmfpen = .TRUE.
  lmfscv = .TRUE.
  lmfmid = .TRUE.
  lmfdd = .TRUE.
  lmfdudv = .TRUE.

  RETURN
END SUBROUTINE flxsetup
