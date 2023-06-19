
! $Header$

SUBROUTINE conlmd(dtime, paprs, pplay, t, q, conv_q, d_t, d_q, rain, snow, &
    ibas, itop)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: Schema de convection utilis'e dans le modele du LMD
  ! Ajustement humide (Manabe) + Ajustement convectif (Kuo)
  ! ======================================================================
  include "YOMCST.h"
  include "YOETHF.h"

  ! Arguments:

  REAL dtime ! pas d'integration (s)
  REAL paprs(klon, klev+1) ! pression inter-couche (Pa)
  REAL pplay(klon, klev) ! pression au milieu de couche (Pa)
  REAL t(klon, klev) ! temperature (K)
  REAL q(klon, klev) ! humidite specifique (kg/kg)
  REAL conv_q(klon, klev) ! taux de convergence humidite (g/g/s)

  REAL d_t(klon, klev) ! incrementation temperature
  REAL d_q(klon, klev) ! incrementation humidite
  REAL rain(klon) ! pluies (mm/s)
  REAL snow(klon) ! neige (mm/s)
  INTEGER ibas(klon) ! niveau du bas
  INTEGER itop(klon) ! niveau du haut

  LOGICAL usekuo ! utiliser convection profonde (schema Kuo)
  PARAMETER (usekuo=.TRUE.)

  REAL d_t_bis(klon, klev)
  REAL d_q_bis(klon, klev)
  REAL rain_bis(klon)
  REAL snow_bis(klon)
  INTEGER ibas_bis(klon)
  INTEGER itop_bis(klon)
  REAL d_ql(klon, klev), d_ql_bis(klon, klev)
  REAL rneb(klon, klev), rneb_bis(klon, klev)

  INTEGER i, k
  REAL zlvdcp, zlsdcp, zdelta, zz, za, zb

  ! cc      CALL fiajh ! ancienne version de Convection Manabe
  CALL conman &                    ! nouvelle version de Convection
                                   ! Manabe
    (dtime, paprs, pplay, t, q, d_t, d_q, d_ql, rneb, rain, snow, ibas, itop)

  IF (usekuo) THEN
    ! cc      CALL fiajc ! ancienne version de Convection Kuo
    CALL conkuo &                  ! nouvelle version de Convection
                                   ! Kuo
      (dtime, paprs, pplay, t, q, conv_q, d_t_bis, d_q_bis, d_ql_bis, &
      rneb_bis, rain_bis, snow_bis, ibas_bis, itop_bis)
    DO k = 1, klev
      DO i = 1, klon
        d_t(i, k) = d_t(i, k) + d_t_bis(i, k)
        d_q(i, k) = d_q(i, k) + d_q_bis(i, k)
        d_ql(i, k) = d_ql(i, k) + d_ql_bis(i, k)
      END DO
    END DO
    DO i = 1, klon
      rain(i) = rain(i) + rain_bis(i)
      snow(i) = snow(i) + snow_bis(i)
      ibas(i) = min(ibas(i), ibas_bis(i))
      itop(i) = max(itop(i), itop_bis(i))
    END DO
  END IF

  ! L'eau liquide convective est dispersee dans l'air:

  DO k = 1, klev
    DO i = 1, klon
      zlvdcp = rlvtt/rcpd/(1.0+rvtmp2*q(i,k))
      zlsdcp = rlstt/rcpd/(1.0+rvtmp2*q(i,k))
      zdelta = max(0., sign(1.,rtt-t(i,k)))
      zz = d_ql(i, k) ! re-evap. de l'eau liquide
      zb = max(0.0, zz)
      za = -max(0.0, zz)*(zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
      d_t(i, k) = d_t(i, k) + za
      d_q(i, k) = d_q(i, k) + zb
    END DO
  END DO

  RETURN
END SUBROUTINE conlmd
SUBROUTINE conman(dtime, paprs, pplay, t, q, d_t, d_q, d_ql, rneb, rain, &
    snow, ibas, itop)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19970324
  ! Objet: ajustement humide convectif avec la possibilite de faire
  ! l'ajustement sur une fraction de la maille.
  ! Methode: On impose une distribution uniforme pour la vapeur d'eau
  ! au sein d'une maille. On applique la procedure d'ajustement
  ! successivement a la totalite, 75%, 50%, 25% et 5% de la maille
  ! jusqu'a ce que l'ajustement a lieu. J'espere que ceci augmente
  ! les activites convectives et corrige le biais "trop froid et sec"
  ! du modele.
  ! ======================================================================
  include "YOMCST.h"

  REAL dtime ! pas d'integration (s)
  REAL t(klon, klev) ! temperature (K)
  REAL q(klon, klev) ! humidite specifique (kg/kg)
  REAL paprs(klon, klev+1) ! pression inter-couche (Pa)
  REAL pplay(klon, klev) ! pression au milieu de couche (Pa)

  REAL d_t(klon, klev) ! incrementation temperature
  REAL d_q(klon, klev) ! incrementation humidite
  REAL d_ql(klon, klev) ! incrementation eau liquide
  REAL rneb(klon, klev) ! nebulosite
  REAL rain(klon) ! pluies (mm/s)
  REAL snow(klon) ! neige (mm/s)
  INTEGER ibas(klon) ! niveau du bas
  INTEGER itop(klon) ! niveau du haut

  LOGICAL afaire(klon) ! .TRUE. implique l'ajustement
  LOGICAL accompli(klon) ! .TRUE. si l'ajustement est effectif

  INTEGER nb ! nombre de sous-fractions a considere
  PARAMETER (nb=1)
  ! cc      PARAMETER (nb=3)

  REAL ratqs ! largeur de la distribution pour vapeur d'eau
  PARAMETER (ratqs=0.05)

  REAL w_q(klon, klev)
  REAL w_d_t(klon, klev), w_d_q(klon, klev), w_d_ql(klon, klev)
  REAL w_rneb(klon, klev)
  REAL w_rain(klon), w_snow(klon)
  INTEGER w_ibas(klon), w_itop(klon)
  REAL zq1, zq2
  INTEGER i, k, n

  REAL t_coup
  PARAMETER (t_coup=234.0)
  REAL zdp1, zdp2
  REAL zqs1, zqs2, zdqs1, zdqs2
  REAL zgamdz
  REAL zflo ! flotabilite
  REAL zsat ! sur-saturation
  REAL zdelta, zcor, zcvm5
  LOGICAL imprim

  INTEGER ncpt
  SAVE ncpt
  !$OMP THREADPRIVATE(ncpt)
  REAL frac(nb) ! valeur de la maille fractionnelle
  SAVE frac
  !$OMP THREADPRIVATE(frac)
  INTEGER opt_cld(nb) ! option pour le modele nuageux
  SAVE opt_cld
  !$OMP THREADPRIVATE(opt_cld)
  LOGICAL appel1er
  SAVE appel1er
  !$OMP THREADPRIVATE(appel1er)

  ! Fonctions thermodynamiques:

  include "YOETHF.h"
  include "FCTTRE.h"

  DATA frac/1.0/
  DATA opt_cld/4/
  ! cc      DATA frac    / 1.0, 0.50, 0.25/
  ! cc      DATA opt_cld / 4,   4,    4/

  DATA appel1er/.TRUE./
  DATA ncpt/0/

  IF (appel1er) THEN
    PRINT *, 'conman, nb:', nb
    PRINT *, 'conman, frac:', frac
    PRINT *, 'conman, opt_cld:', opt_cld
    appel1er = .FALSE.
  END IF

  ! Initialiser les sorties a zero:

  DO k = 1, klev
    DO i = 1, klon
      d_t(i, k) = 0.0
      d_q(i, k) = 0.0
      d_ql(i, k) = 0.0
      rneb(i, k) = 0.0
    END DO
  END DO
  DO i = 1, klon
    ibas(i) = klev
    itop(i) = 1
    rain(i) = 0.0
    snow(i) = 0.0
  END DO

  ! S'il n'y a pas d'instabilite conditionnelle,
  ! pas la penne de se fatiguer:

  DO i = 1, klon
    afaire(i) = .FALSE.
  END DO
  DO k = 1, klev - 1
    DO i = 1, klon
      IF (thermcep) THEN
        zdelta = max(0., sign(1.,rtt-t(i,k)))
        zcvm5 = r5les*rlvtt*(1.-zdelta) + zdelta*r5ies*rlstt
        zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*q(i,k))
        zqs1 = r2es*foeew(t(i,k), zdelta)/pplay(i, k)
        zqs1 = min(0.5, zqs1)
        zcor = 1./(1.-retv*zqs1)
        zqs1 = zqs1*zcor
        zdqs1 = foede(t(i,k), zdelta, zcvm5, zqs1, zcor)

        zdelta = max(0., sign(1.,rtt-t(i,k+1)))
        zcvm5 = r5les*rlvtt*(1.-zdelta) + zdelta*r5ies*rlstt
        zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*q(i,k+1))
        zqs2 = r2es*foeew(t(i,k+1), zdelta)/pplay(i, k+1)
        zqs2 = min(0.5, zqs2)
        zcor = 1./(1.-retv*zqs2)
        zqs2 = zqs2*zcor
        zdqs2 = foede(t(i,k+1), zdelta, zcvm5, zqs2, zcor)
      ELSE
        IF (t(i,k)<t_coup) THEN
          zqs1 = qsats(t(i,k))/pplay(i, k)
          zdqs1 = dqsats(t(i,k), zqs1)

          zqs2 = qsats(t(i,k+1))/pplay(i, k+1)
          zdqs2 = dqsats(t(i,k+1), zqs2)
        ELSE
          zqs1 = qsatl(t(i,k))/pplay(i, k)
          zdqs1 = dqsatl(t(i,k), zqs1)

          zqs2 = qsatl(t(i,k+1))/pplay(i, k+1)
          zdqs2 = dqsatl(t(i,k+1), zqs2)
        END IF
      END IF
      zdp1 = paprs(i, k) - paprs(i, k+1)
      zdp2 = paprs(i, k+1) - paprs(i, k+2)
      zgamdz = -(pplay(i,k)-pplay(i,k+1))/paprs(i, k+1)/rcpd*(rd*(t(i, &
        k)*zdp1+t(i,k+1)*zdp2)/(zdp1+zdp2)+rlvtt*(zqs1*zdp1+zqs2*zdp2)/(zdp1+ &
        zdp2))/(1.0+(zdqs1*zdp1+zdqs2*zdp2)/(zdp1+zdp2))
      zflo = t(i, k) + zgamdz - t(i, k+1)
      zsat = (q(i,k)-zqs1)*zdp1 + (q(i,k+1)-zqs2)*zdp2
      IF (zflo>0.0) afaire(i) = .TRUE.
      ! erreur         IF (zflo.GT.0.0 .AND. zsat.GT.0.0) afaire(i) = .TRUE.
    END DO
  END DO

  imprim = mod(ncpt, 48) == 0
  DO n = 1, nb

    DO k = 1, klev
      DO i = 1, klon
        IF (afaire(i)) THEN
          zq1 = q(i, k)*(1.0-ratqs)
          zq2 = q(i, k)*(1.0+ratqs)
          w_q(i, k) = zq2 - frac(n)/2.0*(zq2-zq1)
        END IF
      END DO
    END DO

    CALL conmanv(dtime, paprs, pplay, t, w_q, afaire, opt_cld(n), w_d_t, &
      w_d_q, w_d_ql, w_rneb, w_rain, w_snow, w_ibas, w_itop, accompli, &
      imprim)
    DO k = 1, klev
      DO i = 1, klon
        IF (afaire(i) .AND. accompli(i)) THEN
          d_t(i, k) = w_d_t(i, k)*frac(n)
          d_q(i, k) = w_d_q(i, k)*frac(n)
          d_ql(i, k) = w_d_ql(i, k)*frac(n)
          IF (nint(w_rneb(i,k))==1) rneb(i, k) = frac(n)
        END IF
      END DO
    END DO
    DO i = 1, klon
      IF (afaire(i) .AND. accompli(i)) THEN
        rain(i) = w_rain(i)*frac(n)
        snow(i) = w_snow(i)*frac(n)
        ibas(i) = min(ibas(i), w_ibas(i))
        itop(i) = max(itop(i), w_itop(i))
      END IF
    END DO
    DO i = 1, klon
      IF (afaire(i) .AND. accompli(i)) afaire(i) = .FALSE.
    END DO

  END DO

  ncpt = ncpt + 1

  RETURN
END SUBROUTINE conman
SUBROUTINE conmanv(dtime, paprs, pplay, t, q, afaire, opt_cld, d_t, d_q, &
    d_ql, rneb, rain, snow, ibas, itop, accompli, imprim)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: ajustement humide (convection proposee par Manabe).
  ! Pour une colonne verticale, il peut avoir plusieurs blocs
  ! necessitant l'ajustement. ibas est le bas du plus bas bloc
  ! et itop est le haut du plus haut bloc
  ! ======================================================================
  include "YOMCST.h"

  ! Arguments:

  REAL dtime ! pas d'integration (s)
  REAL t(klon, klev) ! temperature (K)
  REAL q(klon, klev) ! humidite specifique (kg/kg)
  REAL paprs(klon, klev+1) ! pression inter-couche (Pa)
  REAL pplay(klon, klev) ! pression au milieu de couche (Pa)
  INTEGER opt_cld ! comment traiter l'eau liquide
  LOGICAL afaire(klon) ! .TRUE. si le point est a faire (Input)
  LOGICAL imprim ! .T. pour imprimer quelques diagnostiques

  REAL d_t(klon, klev) ! incrementation temperature
  REAL d_q(klon, klev) ! incrementation humidite
  REAL d_ql(klon, klev) ! incrementation eau liquide
  REAL rneb(klon, klev) ! nebulosite
  REAL rain(klon) ! pluies (mm/s)
  REAL snow(klon) ! neige (mm/s)
  INTEGER ibas(klon) ! niveau du bas
  INTEGER itop(klon) ! niveau du haut
  LOGICAL accompli(klon) ! .TRUE. si l'ajustement a eu lieu (Output)

  ! Quelques options:

  LOGICAL new_top ! re-calculer sommet quand re-ajustement est fait
  PARAMETER (new_top=.FALSE.)
  LOGICAL evap_prec ! evaporation de pluie au-dessous de convection
  PARAMETER (evap_prec=.TRUE.)
  REAL coef_eva
  PARAMETER (coef_eva=1.0E-05)
  REAL t_coup
  PARAMETER (t_coup=234.0)
  REAL seuil_vap
  PARAMETER (seuil_vap=1.0E-10)
  LOGICAL old_tau ! implique precip nulle, si vrai.
  PARAMETER (old_tau=.FALSE.)
  REAL toliq(klon) ! rapport entre l'eau nuageuse et l'eau precipitante
  REAL dpmin, tomax !Epaisseur faible, rapport eau liquide plus grande
  PARAMETER (dpmin=0.15, tomax=0.97)
  REAL dpmax, tomin !Epaisseur grande, rapport eau liquide plus faible
  PARAMETER (dpmax=0.30, tomin=0.05)
  REAL deep_sig, deep_to ! au dela de deep_sig, utiliser deep_to
  PARAMETER (deep_sig=0.50, deep_to=0.05)
  LOGICAL exigent ! implique un calcul supplementaire pour Qs
  PARAMETER (exigent=.FALSE.)

  INTEGER kbase
  PARAMETER (kbase=0)

  ! Variables locales:

  INTEGER nexpo
  INTEGER i, k, k1min, k1max, k2min, k2max, is
  REAL zgamdz(klon, klev-1)
  REAL zt(klon, klev), zq(klon, klev)
  REAL zqs(klon, klev), zdqs(klon, klev)
  REAL zqmqsdp(klon, klev)
  REAL ztnew(klon, klev), zqnew(klon, klev)
  REAL zcond(klon), zvapo(klon), zrapp(klon)
  REAL zrfl(klon), zrfln, zqev, zqevt
  REAL zsat(klon) ! sur-saturation
  REAL zflo(klon) ! flotabilite
  REAL za(klon), zb(klon), zc(klon)
  INTEGER k1(klon), k2(klon)
  REAL zdelta, zcor, zcvm5
  REAL delp(klon, klev)
  LOGICAL possible(klon), todo(klon), etendre(klon)
  LOGICAL aller(klon), todobis(klon)
  REAL zalfa
  INTEGER nbtodo, nbdone

  ! Fonctions thermodynamiques:

  include "YOETHF.h"
  include "FCTTRE.h"

  DO k = 1, klev
    DO i = 1, klon
      delp(i, k) = paprs(i, k) - paprs(i, k+1)
    END DO
  END DO

  ! Initialiser les sorties a zero

  DO k = 1, klev
    DO i = 1, klon
      d_t(i, k) = 0.0
      d_q(i, k) = 0.0
      d_ql(i, k) = 0.0
      rneb(i, k) = 0.0
    END DO
  END DO
  DO i = 1, klon
    ibas(i) = klev
    itop(i) = 1
    rain(i) = 0.0
    snow(i) = 0.0
    accompli(i) = .FALSE.
  END DO

  ! Preparations

  DO k = 1, klev
    DO i = 1, klon
      IF (afaire(i)) THEN
        zt(i, k) = t(i, k)
        zq(i, k) = q(i, k)

        ! Calculer Qs et L/Cp*dQs/dT

        IF (thermcep) THEN
          zdelta = max(0., sign(1.,rtt-zt(i,k)))
          zcvm5 = r5les*rlvtt*(1.-zdelta) + zdelta*r5ies*rlstt
          zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*zq(i,k))
          zqs(i, k) = r2es*foeew(zt(i,k), zdelta)/pplay(i, k)
          zqs(i, k) = min(0.5, zqs(i,k))
          zcor = 1./(1.-retv*zqs(i,k))
          zqs(i, k) = zqs(i, k)*zcor
          zdqs(i, k) = foede(zt(i,k), zdelta, zcvm5, zqs(i,k), zcor)
        ELSE
          IF (zt(i,k)<t_coup) THEN
            zqs(i, k) = qsats(zt(i,k))/pplay(i, k)
            zdqs(i, k) = dqsats(zt(i,k), zqs(i,k))
          ELSE
            zqs(i, k) = qsatl(zt(i,k))/pplay(i, k)
            zdqs(i, k) = dqsatl(zt(i,k), zqs(i,k))
          END IF
        END IF

        ! Calculer (q-qs)*dp
        zqmqsdp(i, k) = (zq(i,k)-zqs(i,k))*delp(i, k)
      END IF
    END DO
  END DO

  ! -----zgama is the moist convective lapse rate (-dT/dz).
  ! -----zgamdz(*,k) est la difference minimale autorisee des temperatures
  ! -----entre deux couches (k et k+1), c.a.d. si T(k+1)-T(k) est inferieur
  ! -----a zgamdz(*,k), alors ces 2 couches sont instables conditionnellement

  DO k = 1, klev - 1
    DO i = 1, klon
      IF (afaire(i)) THEN
        zgamdz(i, k) = -(pplay(i,k)-pplay(i,k+1))/paprs(i, k+1)/rcpd*(rd*(zt( &
          i,k)*delp(i,k)+zt(i,k+1)*delp(i,k+1))/(delp(i,k)+delp(i, &
          k+1))+rlvtt*(zqs(i,k)*delp(i,k)+zqs(i,k+1)*delp(i,k+1))/(delp(i, &
          k)+delp(i,k+1)))/(1.0+(zdqs(i,k)*delp(i,k)+zdqs(i,k+1)*delp(i, &
          k+1))/(delp(i,k)+delp(i,k+1)))
      END IF
    END DO
  END DO

  ! On cherche la presence simultanee d'instabilite conditionnelle
  ! et de sur-saturation. Sinon, pas la penne de se fatiguer:

  DO i = 1, klon
    possible(i) = .FALSE.
  END DO
  DO k = 2, klev
    DO i = 1, klon
      IF (afaire(i)) THEN
        zflo(i) = zt(i, k-1) + zgamdz(i, k-1) - zt(i, k)
        zsat(i) = zqmqsdp(i, k) + zqmqsdp(i, k-1)
        IF (zflo(i)>0.0 .AND. zsat(i)>0.0) possible(i) = .TRUE.
      END IF
    END DO
  END DO

  DO i = 1, klon
    IF (possible(i)) THEN
      k1(i) = kbase
      k2(i) = k1(i) + 1
    END IF
  END DO

810 CONTINUE ! chercher le bas de la colonne a ajuster

  k2min = klev
  DO i = 1, klon
    todo(i) = .FALSE.
    aller(i) = .TRUE.
    IF (possible(i)) k2min = min(k2min, k2(i))
  END DO
  IF (k2min==klev) GO TO 860
  DO k = k2min, klev - 1
    DO i = 1, klon
      IF (possible(i) .AND. k>=k2(i) .AND. aller(i)) THEN
        zflo(i) = zt(i, k) + zgamdz(i, k) - zt(i, k+1)
        zsat(i) = zqmqsdp(i, k) + zqmqsdp(i, k+1)
        IF (zflo(i)>0.0 .AND. zsat(i)>0.0) THEN
          k1(i) = k
          k2(i) = k + 1
          todo(i) = .TRUE.
          aller(i) = .FALSE.
        END IF
      END IF
    END DO
  END DO
  DO i = 1, klon
    IF (possible(i) .AND. aller(i)) THEN
      todo(i) = .FALSE.
      k1(i) = klev
      k2(i) = klev
    END IF
  END DO

  ! CC      DO i = 1, klon
  ! CC      IF (possible(i)) THEN
  ! CC  811    k2(i) = k2(i) + 1
  ! CC         IF (k2(i) .GT. klev) THEN
  ! CC            todo(i) = .FALSE.
  ! CC            GOTO 812
  ! CC         ENDIF
  ! CC         k = k2(i)
  ! CC         zflo(i) = zt(i,k-1) + zgamdz(i,k-1) - zt(i,k)
  ! CC         zsat(i) = zqmqsdp(i,k) + zqmqsdp(i,k-1)
  ! CC         IF (zflo(i).LE.0.0 .OR. zsat(i).LE.0.0) GOTO 811
  ! CC         k1(i) = k2(i) - 1
  ! CC         todo(i) = .TRUE.
  ! CC      ENDIF
  ! CC  812 CONTINUE
  ! CC      ENDDO

820 CONTINUE ! chercher le haut de la colonne

  k2min = klev
  DO i = 1, klon
    aller(i) = .TRUE.
    IF (todo(i)) k2min = min(k2min, k2(i))
  END DO
  IF (k2min<klev) THEN
    DO k = k2min, klev
      DO i = 1, klon
        IF (todo(i) .AND. k>k2(i) .AND. aller(i)) THEN
          zsat(i) = zsat(i) + zqmqsdp(i, k)
          zflo(i) = zt(i, k-1) + zgamdz(i, k-1) - zt(i, k)
          IF (zflo(i)<=0.0 .OR. zsat(i)<=0.0) THEN
            aller(i) = .FALSE.
          ELSE
            k2(i) = k
          END IF
        END IF
      END DO
    END DO
    ! error      is = 0
    ! error      DO i = 1, klon
    ! error      IF(todo(i).AND.aller(i)) THEN
    ! error         is = is + 1
    ! error         todo(i) = .FALSE.
    ! error         k2(i) = klev
    ! error      ENDIF
    ! error      ENDDO
    ! error      IF (is.GT.0) THEN
    ! error         PRINT*, "Bizard. je pourrais continuer mais j arrete"
    ! error         CALL abort
    ! error      ENDIF
  END IF

  ! CC      DO i = 1, klon
  ! CC      IF (todo(i)) THEN
  ! CC  821    CONTINUE
  ! CC         IF (k2(i) .EQ. klev) GOTO 822
  ! CC         k = k2(i) + 1
  ! CC         zsat(i) = zsat(i) + zqmqsdp(i,k)
  ! CC         zflo(i) = zt(i,k-1) + zgamdz(i,k-1) - zt(i,k)
  ! CC         IF (zflo(i).LE.0.0 .OR. zsat(i).LE.0.0) GOTO 822
  ! CC         k2(i) = k
  ! CC         GOTO 821
  ! CC      ENDIF
  ! CC  822 CONTINUE
  ! CC      ENDDO

830 CONTINUE ! faire l'ajustement en sachant k1 et k2

  is = 0
  DO i = 1, klon
    IF (todo(i)) THEN
      IF (k2(i)<=k1(i)) is = is + 1
    END IF
  END DO
  IF (is>0) THEN
    PRINT *, 'Impossible: k1 trop grand ou k2 trop petit'
    PRINT *, 'is=', is
    CALL abort
  END IF

  k1min = klev
  k1max = 1
  k2max = 1
  DO i = 1, klon
    IF (todo(i)) THEN
      k1min = min(k1min, k1(i))
      k1max = max(k1max, k1(i))
      k2max = max(k2max, k2(i))
    END IF
  END DO

  DO i = 1, klon
    IF (todo(i)) THEN
      k = k1(i)
      za(i) = 0.
      zb(i) = (rcpd*(1.+zdqs(i,k))*(zt(i,k)-za(i))-rlvtt*(zqs(i,k)-zq(i, &
        k)))*delp(i, k)
      zc(i) = delp(i, k)*rcpd*(1.+zdqs(i,k))
    END IF
  END DO

  DO k = k1min, k2max
    DO i = 1, klon
      IF (todo(i) .AND. k>=(k1(i)+1) .AND. k<=k2(i)) THEN
        za(i) = za(i) + zgamdz(i, k-1)
        zb(i) = zb(i) + (rcpd*(1.+zdqs(i,k))*(zt(i,k)-za(i))-rlvtt*(zqs(i, &
          k)-zq(i,k)))*delp(i, k)
        zc(i) = zc(i) + delp(i, k)*rcpd*(1.+zdqs(i,k))
      END IF
    END DO
  END DO

  DO i = 1, klon
    IF (todo(i)) THEN
      k = k1(i)
      ztnew(i, k) = zb(i)/zc(i)
      zqnew(i, k) = zqs(i, k) + (ztnew(i,k)-zt(i,k))*rcpd/rlvtt*zdqs(i, k)
    END IF
  END DO

  DO k = k1min, k2max
    DO i = 1, klon
      IF (todo(i) .AND. k>=(k1(i)+1) .AND. k<=k2(i)) THEN
        ztnew(i, k) = ztnew(i, k-1) + zgamdz(i, k-1)
        zqnew(i, k) = zqs(i, k) + (ztnew(i,k)-zt(i,k))*rcpd/rlvtt*zdqs(i, k)
      END IF
    END DO
  END DO

  ! Quantite de condensation produite pendant l'ajustement:

  DO i = 1, klon
    zcond(i) = 0.0
  END DO
  DO k = k1min, k2max
    DO i = 1, klon
      IF (todo(i) .AND. k>=k1(i) .AND. k<=k2(i)) THEN
        rneb(i, k) = 1.0
        zcond(i) = zcond(i) + (zq(i,k)-zqnew(i,k))*delp(i, k)/rg
      END IF
    END DO
  END DO

  ! Si condensation negative, effort completement perdu:

  DO i = 1, klon
    IF (todo(i) .AND. zcond(i)<=0.) todo(i) = .FALSE.
  END DO

  ! L'ajustement a ete accompli, meme les calculs accessoires
  ! ne sont pas encore faits:

  DO i = 1, klon
    IF (todo(i)) accompli(i) = .TRUE.
  END DO

  ! =====
  ! Une fois que la condensation a lieu, on doit construire un
  ! "modele nuageux" pour partager la condensation entre l'eau
  ! liquide nuageuse et la precipitation (leur rapport toliq
  ! est calcule selon l'epaisseur nuageuse). Je suppose que
  ! toliq=tomax quand l'epaisseur nuageuse est inferieure a dpmin,
  ! et que toliq=tomin quand l'epaisseur depasse dpmax (interpolation
  ! lineaire entre dpmin et dpmax).
  ! =====
  DO i = 1, klon
    IF (todo(i)) THEN
      toliq(i) = tomax - ((paprs(i,k1(i))-paprs(i,k2(i)+1))/paprs(i,1)-dpmin) &
        *(tomax-tomin)/(dpmax-dpmin)
      toliq(i) = max(tomin, min(tomax,toliq(i)))
      IF (pplay(i,k2(i))/paprs(i,1)<=deep_sig) toliq(i) = deep_to
      IF (old_tau) toliq(i) = 1.0
    END IF
  END DO
  ! =====
  ! On doit aussi determiner la distribution verticale de
  ! l'eau nuageuse. Plusieurs options sont proposees:

  ! (0) La condensation precipite integralement (toliq ne sera
  ! pas utilise).
  ! (1) L'eau liquide est distribuee entre k1 et k2 et proportionnelle
  ! a la vapeur d'eau locale.
  ! (2) Elle est distribuee entre k1 et k2 avec une valeur constante.
  ! (3) Elle est seulement distribuee aux couches ou la vapeur d'eau
  ! est effectivement diminuee pendant le processus d'ajustement.
  ! (4) Elle est en fonction (lineaire ou exponentielle) de la
  ! distance (epaisseur en pression) avec le niveau k1 (la couche
  ! k1 n'aura donc pas d'eau liquide).
  ! =====

  IF (opt_cld==0) THEN

    DO i = 1, klon
      IF (todo(i)) zrfl(i) = zcond(i)/dtime
    END DO

  ELSE IF (opt_cld==1) THEN

    DO i = 1, klon
      IF (todo(i)) zvapo(i) = 0.0 ! quantite integrale de vapeur d'eau
    END DO
    DO k = k1min, k2max
      DO i = 1, klon
        IF (todo(i) .AND. k>=k1(i) .AND. k<=k2(i)) zvapo(i) = zvapo(i) + &
          zqnew(i, k)*delp(i, k)/rg
      END DO
    END DO
    DO i = 1, klon
      IF (todo(i)) THEN
        zrapp(i) = toliq(i)*zcond(i)/zvapo(i)
        zrapp(i) = max(0., min(1.,zrapp(i)))
        zrfl(i) = (1.0-toliq(i))*zcond(i)/dtime
      END IF
    END DO
    DO k = k1min, k2max
      DO i = 1, klon
        IF (todo(i) .AND. k>=k1(i) .AND. k<=k2(i)) THEN
          d_ql(i, k) = d_ql(i, k) + zrapp(i)*zqnew(i, k)
        END IF
      END DO
    END DO

  ELSE IF (opt_cld==2) THEN

    DO i = 1, klon
      IF (todo(i)) zvapo(i) = 0.0 ! quantite integrale de masse
    END DO
    DO k = k1min, k2max
      DO i = 1, klon
        IF (todo(i) .AND. k>=k1(i) .AND. k<=k2(i)) zvapo(i) = zvapo(i) + &
          delp(i, k)/rg
      END DO
    END DO
    DO k = k1min, k2max
      DO i = 1, klon
        IF (todo(i) .AND. k>=k1(i) .AND. k<=k2(i)) THEN
          d_ql(i, k) = d_ql(i, k) + toliq(i)*zcond(i)/zvapo(i)
        END IF
      END DO
    END DO
    DO i = 1, klon
      IF (todo(i)) zrfl(i) = (1.0-toliq(i))*zcond(i)/dtime
    END DO

  ELSE IF (opt_cld==3) THEN

    DO i = 1, klon
      IF (todo(i)) zvapo(i) = 0.0 ! quantite de l'eau strictement condensee
    END DO
    DO k = k1min, k2max
      DO i = 1, klon
        IF (todo(i) .AND. k>=k1(i) .AND. k<=k2(i)) zvapo(i) = zvapo(i) + &
          max(0.0, zq(i,k)-zqnew(i,k))*delp(i, k)/rg
      END DO
    END DO
    DO k = k1min, k2max
      DO i = 1, klon
        IF (todo(i) .AND. k>=k1(i) .AND. k<=k2(i) .AND. zvapo(i)>0.0) d_ql(i, &
          k) = d_ql(i, k) + toliq(i)*zcond(i)/zvapo(i)*max(0.0, zq(i,k)-zqnew &
          (i,k))
      END DO
    END DO
    DO i = 1, klon
      IF (todo(i)) zrfl(i) = (1.0-toliq(i))*zcond(i)/dtime
    END DO

  ELSE IF (opt_cld==4) THEN

    nexpo = 3
    ! cc         nexpo = 1 ! distribution lineaire

    DO i = 1, klon
      IF (todo(i)) zvapo(i) = 0.0 ! quantite integrale de masse
    END DO ! (avec ponderation)
    DO k = k1min, k2max
      DO i = 1, klon
        IF (todo(i) .AND. k>=(k1(i)+1) .AND. k<=k2(i)) zvapo(i) = zvapo(i) + &
          delp(i, k)/rg*(pplay(i,k1(i))-pplay(i,k))**nexpo
      END DO
    END DO
    DO k = k1min, k2max
      DO i = 1, klon
        IF (todo(i) .AND. k>=(k1(i)+1) .AND. k<=k2(i)) d_ql(i, k) = d_ql(i, &
          k) + toliq(i)*zcond(i)/zvapo(i)*(pplay(i,k1(i))-pplay(i,k))**nexpo
      END DO
    END DO
    DO i = 1, klon
      IF (todo(i)) zrfl(i) = (1.0-toliq(i))*zcond(i)/dtime
    END DO

  ELSE ! valeur non-prevue pour opt_cld

    PRINT *, 'opt_cld est faux:', opt_cld
    CALL abort

  END IF ! fin de opt_cld

  ! L'eau precipitante peut etre evaporee:

  zalfa = 0.05
  IF (evap_prec .AND. (k1max>=2)) THEN
    DO k = k1max - 1, 1, -1
      DO i = 1, klon
        IF (todo(i) .AND. k<k1(i) .AND. zrfl(i)>0.0) THEN
          zqev = max(0.0, (zqs(i,k)-zq(i,k))*zalfa)
          zqevt = coef_eva*(1.0-zq(i,k)/zqs(i,k))*sqrt(zrfl(i))*delp(i, k)/ &
            pplay(i, k)*zt(i, k)*rd/rg
          zqevt = max(0.0, min(zqevt,zrfl(i)))*rg*dtime/delp(i, k)
          zqev = min(zqev, zqevt)
          zrfln = zrfl(i) - zqev*(delp(i,k))/rg/dtime
          zq(i, k) = zq(i, k) - (zrfln-zrfl(i))*(rg/(delp(i,k)))*dtime
          zt(i, k) = zt(i, k) + (zrfln-zrfl(i))*(rg/(delp(i, &
            k)))*dtime*rlvtt/rcpd/(1.0+rvtmp2*zq(i,k))
          zrfl(i) = zrfln
        END IF
      END DO
    END DO
  END IF

  ! La temperature de la premiere couche determine la pluie ou la neige:

  DO i = 1, klon
    IF (todo(i)) THEN
      IF (zt(i,1)>rtt) THEN
        rain(i) = rain(i) + zrfl(i)
      ELSE
        snow(i) = snow(i) + zrfl(i)
      END IF
    END IF
  END DO

  ! Mise a jour de la temperature et de l'humidite

  DO k = k1min, k2max
    DO i = 1, klon
      IF (todo(i) .AND. k>=k1(i) .AND. k<=k2(i)) THEN
        zt(i, k) = ztnew(i, k)
        zq(i, k) = zqnew(i, k)
      END IF
    END DO
  END DO

  ! Re-calculer certaines variables pour etendre et re-ajuster la colonne

  IF (exigent) THEN
    DO k = 1, klev
      DO i = 1, klon
        IF (todo(i)) THEN
          IF (thermcep) THEN
            zdelta = max(0., sign(1.,rtt-zt(i,k)))
            zcvm5 = r5les*rlvtt*(1.-zdelta) + zdelta*r5ies*rlstt
            zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*zq(i,k))
            zqs(i, k) = r2es*foeew(zt(i,k), zdelta)/pplay(i, k)
            zqs(i, k) = min(0.5, zqs(i,k))
            zcor = 1./(1.-retv*zqs(i,k))
            zqs(i, k) = zqs(i, k)*zcor
            zdqs(i, k) = foede(zt(i,k), zdelta, zcvm5, zqs(i,k), zcor)
          ELSE
            IF (zt(i,k)<t_coup) THEN
              zqs(i, k) = qsats(zt(i,k))/pplay(i, k)
              zdqs(i, k) = dqsats(zt(i,k), zqs(i,k))
            ELSE
              zqs(i, k) = qsatl(zt(i,k))/pplay(i, k)
              zdqs(i, k) = dqsatl(zt(i,k), zqs(i,k))
            END IF
          END IF
        END IF
      END DO
    END DO
  END IF

  IF (exigent) THEN
    DO k = 1, klev - 1
      DO i = 1, klon
        IF (todo(i)) THEN
          zgamdz(i, k) = -(pplay(i,k)-pplay(i,k+1))/paprs(i, k+1)/rcpd*(rd*( &
            zt(i,k)*delp(i,k)+zt(i,k+1)*delp(i,k+1))/(delp(i,k)+delp(i, &
            k+1))+rlvtt*(zqs(i,k)*delp(i,k)+zqs(i,k+1)*delp(i,k+1))/(delp(i, &
            k)+delp(i,k+1)))/(1.0+(zdqs(i,k)*delp(i,k)+zdqs(i,k+1)*delp(i, &
            k+1))/(delp(i,k)+delp(i,k+1)))
        END IF
      END DO
    END DO
  END IF

  ! Puisque l'humidite a ete modifiee, on re-fait (q-qs)*dp

  DO k = 1, klev
    DO i = 1, klon
      IF (todo(i)) THEN
        zqmqsdp(i, k) = (zq(i,k)-zqs(i,k))*delp(i, k)
      END IF
    END DO
  END DO

  ! Verifier si l'on peut etendre le bas de la colonne

  DO i = 1, klon
    etendre(i) = .FALSE.
  END DO

  k1max = 1
  DO i = 1, klon
    IF (todo(i) .AND. k1(i)>(kbase+1)) THEN
      k = k1(i)
      zflo(i) = zt(i, k-1) + zgamdz(i, k-1) - zt(i, k)
      zsat(i) = zqmqsdp(i, k) + zqmqsdp(i, k-1)
      ! sc voici l'ancienne ligne:
      ! sc         IF (zflo(i).LE.0.0 .OR. zsat(i).LE.0.0) THEN
      ! sc sylvain: il faut RESPECTER les 2 criteres:
      IF (zflo(i)>0.0 .AND. zsat(i)>0.0) THEN
        etendre(i) = .TRUE.
        k1(i) = k1(i) - 1
        k1max = max(k1max, k1(i))
        aller(i) = .TRUE.
      END IF
    END IF
  END DO

  IF (k1max>(kbase+1)) THEN
    DO k = k1max, kbase + 1, -1
      DO i = 1, klon
        IF (etendre(i) .AND. k<k1(i) .AND. aller(i)) THEN
          zsat(i) = zsat(i) + zqmqsdp(i, k)
          zflo(i) = zt(i, k) + zgamdz(i, k) - zt(i, k+1)
          IF (zsat(i)<=0.0 .OR. zflo(i)<=0.0) THEN
            aller(i) = .FALSE.
          ELSE
            k1(i) = k
          END IF
        END IF
      END DO
    END DO
    DO i = 1, klon
      IF (etendre(i) .AND. aller(i)) THEN
        k1(i) = 1
      END IF
    END DO
  END IF

  ! CC      DO i = 1, klon
  ! CC      IF (etendre(i)) THEN
  ! CC  840    k = k1(i)
  ! CC         IF (k.GT.1) THEN
  ! CC            zsat(i) = zsat(i) + zqmqsdp(i,k-1)
  ! CC            zflo(i) = zt(i,k-1) + zgamdz(i,k-1) - zt(i,k)
  ! CC            IF (zflo(i).GT.0.0 .AND. zsat(i).GT.0.0) THEN
  ! CC               k1(i) = k - 1
  ! CC               GOTO 840
  ! CC            ENDIF
  ! CC         ENDIF
  ! CC      ENDIF
  ! CC      ENDDO

  DO i = 1, klon
    todobis(i) = todo(i)
    todo(i) = .FALSE.
  END DO
  is = 0
  DO i = 1, klon
    IF (etendre(i)) THEN
      todo(i) = .TRUE.
      is = is + 1
    END IF
  END DO
  IF (is>0) THEN
    IF (new_top) THEN
      GO TO 820 ! chercher de nouveau le sommet k2
    ELSE
      GO TO 830 ! supposer que le sommet est celui deja trouve
    END IF
  END IF

  DO i = 1, klon
    possible(i) = .FALSE.
  END DO
  is = 0
  DO i = 1, klon
    IF (todobis(i) .AND. k2(i)<klev) THEN
      is = is + 1
      possible(i) = .TRUE.
    END IF
  END DO
  IF (is>0) GO TO 810 !on cherche en haut d'autres blocks
  ! a ajuster a partir du sommet de la colonne precedente

860 CONTINUE ! Calculer les tendances et diagnostiques
  ! cc      print*, "Apres 860"

  DO k = 1, klev
    DO i = 1, klon
      IF (accompli(i)) THEN
        d_t(i, k) = zt(i, k) - t(i, k)
        zq(i, k) = max(zq(i,k), seuil_vap)
        d_q(i, k) = zq(i, k) - q(i, k)
      END IF
    END DO
  END DO

  DO i = 1, klon
    IF (accompli(i)) THEN
      DO k = 1, klev
        IF (rneb(i,k)>0.0) THEN
          ibas(i) = k
          GO TO 807
        END IF
      END DO
807   CONTINUE
      DO k = klev, 1, -1
        IF (rneb(i,k)>0.0) THEN
          itop(i) = k
          GO TO 808
        END IF
      END DO
808   CONTINUE
    END IF
  END DO

  IF (imprim) THEN
    nbtodo = 0
    nbdone = 0
    DO i = 1, klon
      IF (afaire(i)) nbtodo = nbtodo + 1
      IF (accompli(i)) nbdone = nbdone + 1
    END DO
    PRINT *, 'nbTodo, nbDone=', nbtodo, nbdone
  END IF

  RETURN
END SUBROUTINE conmanv
SUBROUTINE conkuo(dtime, paprs, pplay, t, q, conv_q, d_t, d_q, d_ql, rneb, &
    rain, snow, ibas, itop)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: Schema de convection de type Kuo (1965).
  ! Cette version du code peut calculer le niveau de depart
  ! N.B. version vectorielle (le 6 oct. 1997)
  ! ======================================================================
  include "YOMCST.h"

  ! Arguments:

  REAL dtime ! intervalle du temps (s)
  REAL paprs(klon, klev+1) ! pression a inter-couche (Pa)
  REAL pplay(klon, klev) ! pression au milieu de couche (Pa)
  REAL t(klon, klev) ! temperature (K)
  REAL q(klon, klev) ! humidite specifique
  REAL conv_q(klon, klev) ! taux de convergence humidite (g/g/s)

  REAL d_t(klon, klev) ! incrementation temperature
  REAL d_q(klon, klev) ! incrementation humidite
  REAL d_ql(klon, klev) ! incrementation eau liquide
  REAL rneb(klon, klev) ! nebulosite
  REAL rain(klon) ! pluies (mm/s)
  REAL snow(klon) ! neige (mm/s)
  INTEGER itop(klon) ! niveau du sommet
  INTEGER ibas(klon) ! niveau du bas

  LOGICAL ldcum(klon) ! convection existe
  LOGICAL todo(klon)

  ! Quelsques options:

  LOGICAL calcfcl ! calculer le niveau de convection libre
  PARAMETER (calcfcl=.TRUE.)
  INTEGER ldepar ! niveau fixe de convection libre
  PARAMETER (ldepar=4)
  INTEGER opt_cld ! comment traiter l'eau liquide
  PARAMETER (opt_cld=4) ! valeur possible: 0, 1, 2, 3 ou 4
  LOGICAL evap_prec ! evaporation de pluie au-dessous de convection
  PARAMETER (evap_prec=.TRUE.)
  REAL coef_eva
  PARAMETER (coef_eva=1.0E-05)
  LOGICAL new_deh ! nouvelle facon de calculer dH
  PARAMETER (new_deh=.FALSE.)
  REAL t_coup
  PARAMETER (t_coup=234.0)
  LOGICAL old_tau ! implique precipitation nulle
  PARAMETER (old_tau=.FALSE.)
  REAL toliq(klon) ! rapport entre l'eau nuageuse et l'eau precipitante
  REAL dpmin, tomax !Epaisseur faible, rapport eau liquide plus grande
  PARAMETER (dpmin=0.15, tomax=0.97)
  REAL dpmax, tomin !Epaisseur grande, rapport eau liquide plus faible
  PARAMETER (dpmax=0.30, tomin=0.05)
  REAL deep_sig, deep_to ! au dela de deep_sig, utiliser deep_to
  PARAMETER (deep_sig=0.50, deep_to=0.05)

  ! Variables locales:

  INTEGER nexpo
  LOGICAL nuage(klon)
  INTEGER i, k, kbmin, kbmax, khmax
  REAL ztotal(klon, klev), zdeh(klon, klev)
  REAL zgz(klon, klev)
  REAL zqs(klon, klev)
  REAL zdqs(klon, klev)
  REAL ztemp(klon, klev)
  REAL zpres(klon, klev)
  REAL zconv(klon) ! convergence d'humidite
  REAL zvirt(klon) ! convergence virtuelle d'humidite
  REAL zfrac(klon) ! fraction convective
  INTEGER kb(klon), kh(klon)

  REAL zcond(klon), zvapo(klon), zrapp(klon)
  REAL zrfl(klon), zrfln, zqev, zqevt
  REAL zdelta, zcvm5, zcor
  REAL zvar

  LOGICAL appel1er
  SAVE appel1er
  !$OMP THREADPRIVATE(appel1er)

  ! Fonctions thermodynamiques

  include "YOETHF.h"
  include "FCTTRE.h"

  DATA appel1er/.TRUE./

  IF (appel1er) THEN
    PRINT *, 'conkuo, calcfcl:', calcfcl
    IF (.NOT. calcfcl) PRINT *, 'conkuo, ldepar:', ldepar
    PRINT *, 'conkuo, opt_cld:', opt_cld
    PRINT *, 'conkuo, evap_prec:', evap_prec
    PRINT *, 'conkuo, new_deh:', new_deh
    appel1er = .FALSE.
  END IF

  ! Initialiser les sorties a zero

  DO k = 1, klev
    DO i = 1, klon
      d_q(i, k) = 0.0
      d_t(i, k) = 0.0
      d_ql(i, k) = 0.0
      rneb(i, k) = 0.0
    END DO
  END DO
  DO i = 1, klon
    rain(i) = 0.0
    snow(i) = 0.0
    ibas(i) = 0
    itop(i) = 0
  END DO

  ! Calculer la vapeur d'eau saturante Qs et sa derive L/Cp * dQs/dT

  DO k = 1, klev
    DO i = 1, klon
      IF (thermcep) THEN
        zdelta = max(0., sign(1.,rtt-t(i,k)))
        zcvm5 = r5les*rlvtt*(1.-zdelta) + zdelta*r5ies*rlstt
        zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*q(i,k))
        zqs(i, k) = r2es*foeew(t(i,k), zdelta)/pplay(i, k)
        zqs(i, k) = min(0.5, zqs(i,k))
        zcor = 1./(1.-retv*zqs(i,k))
        zqs(i, k) = zqs(i, k)*zcor
        zdqs(i, k) = foede(t(i,k), zdelta, zcvm5, zqs(i,k), zcor)
      ELSE
        IF (t(i,k)<t_coup) THEN
          zqs(i, k) = qsats(t(i,k))/pplay(i, k)
          zdqs(i, k) = dqsats(t(i,k), zqs(i,k))
        ELSE
          zqs(i, k) = qsatl(t(i,k))/pplay(i, k)
          zdqs(i, k) = dqsatl(t(i,k), zqs(i,k))
        END IF
      END IF
    END DO
  END DO

  ! Calculer gz (energie potentielle)

  DO i = 1, klon
    zgz(i, 1) = rd*t(i, 1)/(0.5*(paprs(i,1)+pplay(i, &
      1)))*(paprs(i,1)-pplay(i,1))
  END DO
  DO k = 2, klev
    DO i = 1, klon
      zgz(i, k) = zgz(i, k-1) + rd*0.5*(t(i,k-1)+t(i,k))/paprs(i, k)*(pplay(i &
        ,k-1)-pplay(i,k))
    END DO
  END DO

  ! Calculer l'energie statique humide saturee (Cp*T + gz + L*Qs)

  DO k = 1, klev
    DO i = 1, klon
      ztotal(i, k) = rcpd*t(i, k) + rlvtt*zqs(i, k) + zgz(i, k)
    END DO
  END DO

  ! Determiner le niveau de depart et calculer la difference de
  ! l'energie statique humide saturee (ztotal) entre la couche
  ! de depart et chaque couche au-dessus.

  IF (calcfcl) THEN
    DO k = 1, klev
      DO i = 1, klon
        zpres(i, k) = pplay(i, k)
        ztemp(i, k) = t(i, k)
      END DO
    END DO
    CALL kuofcl(ztemp, q, zgz, zpres, ldcum, kb)
    DO i = 1, klon
      IF (ldcum(i)) THEN
        k = kb(i)
        IF (new_deh) THEN
          zdeh(i, k) = ztotal(i, k-1) - ztotal(i, k)
        ELSE
          zdeh(i, k) = rcpd*(t(i,k-1)-t(i,k)) - rd*0.5*(t(i,k-1)+t(i,k))/ &
            paprs(i, k)*(pplay(i,k-1)-pplay(i,k)) + &
            rlvtt*(zqs(i,k-1)-zqs(i,k))
        END IF
        zdeh(i, k) = zdeh(i, k)*0.5
      END IF
    END DO
    DO k = 1, klev
      DO i = 1, klon
        IF (ldcum(i) .AND. k>=(kb(i)+1)) THEN
          IF (new_deh) THEN
            zdeh(i, k) = zdeh(i, k-1) + (ztotal(i,k-1)-ztotal(i,k))
          ELSE
            zdeh(i, k) = zdeh(i, k-1) + rcpd*(t(i,k-1)-t(i,k)) - &
              rd*0.5*(t(i,k-1)+t(i,k))/paprs(i, k)* &
              (pplay(i,k-1)-pplay(i,k)) + rlvtt*(zqs(i,k-1)-zqs(i,k))
          END IF
        END IF
      END DO
    END DO
  ELSE
    DO i = 1, klon
      k = ldepar
      kb(i) = ldepar
      ldcum(i) = .TRUE.
      IF (new_deh) THEN
        zdeh(i, k) = ztotal(i, k-1) - ztotal(i, k)
      ELSE
        zdeh(i, k) = rcpd*(t(i,k-1)-t(i,k)) - rd*0.5*(t(i,k-1)+t(i,k))/paprs( &
          i, k)*(pplay(i,k-1)-pplay(i,k)) + rlvtt*(zqs(i,k-1)-zqs(i,k))
      END IF
      zdeh(i, k) = zdeh(i, k)*0.5
    END DO
    DO k = ldepar + 1, klev
      DO i = 1, klon
        IF (new_deh) THEN
          zdeh(i, k) = zdeh(i, k-1) + (ztotal(i,k-1)-ztotal(i,k))
        ELSE
          zdeh(i, k) = zdeh(i, k-1) + rcpd*(t(i,k-1)-t(i,k)) - &
            rd*0.5*(t(i,k-1)+t(i,k))/paprs(i, k)*(pplay(i,k-1)-pplay(i,k)) + &
            rlvtt*(zqs(i,k-1)-zqs(i,k))
        END IF
      END DO
    END DO
  END IF

  ! -----Chercher le sommet du nuage
  ! -----Calculer la convergence de l'humidite (en kg/m**2 a un facteur
  ! -----psolpa/RG pres) du bas jusqu'au sommet du nuage.
  ! -----Calculer la convergence virtuelle pour que toute la maille
  ! -----deviennt nuageuse (du bas jusqu'au sommet du nuage)

  DO i = 1, klon
    nuage(i) = .TRUE.
    zconv(i) = 0.0
    zvirt(i) = 0.0
    kh(i) = -999
  END DO
  DO k = 1, klev
    DO i = 1, klon
      IF (k>=kb(i) .AND. ldcum(i)) THEN
        nuage(i) = nuage(i) .AND. zdeh(i, k) > 0.0
        IF (nuage(i)) THEN
          kh(i) = k
          zconv(i) = zconv(i) + conv_q(i, k)*dtime*(paprs(i,k)-paprs(i,k+1))
          zvirt(i) = zvirt(i) + (zdeh(i,k)/rlvtt+zqs(i,k)-q(i,k))*(paprs(i,k) &
            -paprs(i,k+1))
        END IF
      END IF
    END DO
  END DO

  DO i = 1, klon
    todo(i) = ldcum(i) .AND. kh(i) > kb(i) .AND. zconv(i) > 0.0
  END DO

  kbmin = klev
  kbmax = 0
  khmax = 0
  DO i = 1, klon
    IF (todo(i)) THEN
      kbmin = min(kbmin, kb(i))
      kbmax = max(kbmax, kb(i))
      khmax = max(khmax, kh(i))
    END IF
  END DO

  ! -----Calculer la surface couverte par le nuage

  DO i = 1, klon
    IF (todo(i)) THEN
      zfrac(i) = max(0.0, min(zconv(i)/zvirt(i),1.0))
    END IF
  END DO

  ! -----Calculs essentiels:

  DO i = 1, klon
    IF (todo(i)) THEN
      zcond(i) = 0.0
    END IF
  END DO
  DO k = kbmin, khmax
    DO i = 1, klon
      IF (todo(i) .AND. k>=kb(i) .AND. k<=kh(i)) THEN
        zvar = zdeh(i, k)/(1.+zdqs(i,k))
        d_t(i, k) = zvar*zfrac(i)/rcpd
        d_q(i, k) = (zvar*zdqs(i,k)/rlvtt+zqs(i,k)-q(i,k))*zfrac(i) - &
          conv_q(i, k)*dtime
        zcond(i) = zcond(i) - d_q(i, k)*(paprs(i,k)-paprs(i,k+1))/rg
        rneb(i, k) = zfrac(i)
      END IF
    END DO
  END DO

  DO i = 1, klon
    IF (todo(i) .AND. zcond(i)<0.0) THEN
      PRINT *, 'WARNING: cond. negative (Kuo) ', i, kb(i), kh(i), zcond(i)
      zcond(i) = 0.0
      DO k = kb(i), kh(i)
        d_t(i, k) = 0.0
        d_q(i, k) = 0.0
      END DO
      todo(i) = .FALSE. ! effort totalement perdu
    END IF
  END DO

  ! =====
  ! Une fois que la condensation a lieu, on doit construire un
  ! "modele nuageux" pour partager la condensation entre l'eau
  ! liquide nuageuse et la precipitation (leur rapport toliq
  ! est calcule selon l'epaisseur nuageuse). Je suppose que
  ! toliq=tomax quand l'epaisseur nuageuse est inferieure a dpmin,
  ! et que toliq=tomin quand l'epaisseur depasse dpmax (interpolation
  ! lineaire entre dpmin et dpmax).
  ! =====
  DO i = 1, klon
    IF (todo(i)) THEN
      toliq(i) = tomax - ((paprs(i,kb(i))-paprs(i,kh(i)+1))/paprs(i,1)-dpmin) &
        *(tomax-tomin)/(dpmax-dpmin)
      toliq(i) = max(tomin, min(tomax,toliq(i)))
      IF (pplay(i,kh(i))/paprs(i,1)<=deep_sig) toliq(i) = deep_to
      IF (old_tau) toliq(i) = 1.0
    END IF
  END DO
  ! =====
  ! On doit aussi determiner la distribution verticale de
  ! l'eau nuageuse. Plusieurs options sont proposees:

  ! (0) La condensation precipite integralement (toliq ne sera
  ! pas utilise).
  ! (1) L'eau liquide est distribuee entre k1 et k2 et proportionnelle
  ! a la vapeur d'eau locale.
  ! (2) Elle est distribuee entre k1 et k2 avec une valeur constante.
  ! (3) Elle est seulement distribuee aux couches ou la vapeur d'eau
  ! est effectivement diminuee pendant le processus d'ajustement.
  ! (4) Elle est en fonction (lineaire ou exponentielle) de la
  ! distance (epaisseur en pression) avec le niveau k1 (la couche
  ! k1 n'aura donc pas d'eau liquide).
  ! =====

  IF (opt_cld==0) THEN

    DO i = 1, klon
      IF (todo(i)) zrfl(i) = zcond(i)/dtime
    END DO

  ELSE IF (opt_cld==1) THEN

    DO i = 1, klon
      IF (todo(i)) zvapo(i) = 0.0 ! quantite integrale de vapeur d'eau
    END DO
    DO k = kbmin, khmax
      DO i = 1, klon
        IF (todo(i) .AND. k>=kb(i) .AND. k<=kh(i)) THEN
          zvapo(i) = zvapo(i) + (q(i,k)+d_q(i,k))*(paprs(i,k)-paprs(i,k+1))/ &
            rg
        END IF
      END DO
    END DO
    DO i = 1, klon
      IF (todo(i)) THEN
        zrapp(i) = toliq(i)*zcond(i)/zvapo(i)
        zrapp(i) = max(0., min(1.,zrapp(i)))
      END IF
    END DO
    DO k = kbmin, khmax
      DO i = 1, klon
        IF (todo(i) .AND. k>=kb(i) .AND. k<=kh(i)) THEN
          d_ql(i, k) = zrapp(i)*(q(i,k)+d_q(i,k))
        END IF
      END DO
    END DO
    DO i = 1, klon
      IF (todo(i)) THEN
        zrfl(i) = (1.0-toliq(i))*zcond(i)/dtime
      END IF
    END DO

  ELSE IF (opt_cld==2) THEN

    DO i = 1, klon
      IF (todo(i)) zvapo(i) = 0.0 ! quantite integrale de masse
    END DO
    DO k = kbmin, khmax
      DO i = 1, klon
        IF (todo(i) .AND. k>=kb(i) .AND. k<=kh(i)) THEN
          zvapo(i) = zvapo(i) + (paprs(i,k)-paprs(i,k+1))/rg
        END IF
      END DO
    END DO
    DO k = kbmin, khmax
      DO i = 1, klon
        IF (todo(i) .AND. k>=kb(i) .AND. k<=kh(i)) THEN
          d_ql(i, k) = toliq(i)*zcond(i)/zvapo(i)
        END IF
      END DO
    END DO
    DO i = 1, klon
      IF (todo(i)) THEN
        zrfl(i) = (1.0-toliq(i))*zcond(i)/dtime
      END IF
    END DO

  ELSE IF (opt_cld==3) THEN

    DO i = 1, klon
      IF (todo(i)) THEN
        zvapo(i) = 0.0 ! quantite de l'eau strictement condensee
      END IF
    END DO
    DO k = kbmin, khmax
      DO i = 1, klon
        IF (todo(i) .AND. k>=kb(i) .AND. k<=kh(i)) THEN
          zvapo(i) = zvapo(i) + max(0.0, -d_q(i,k))*(paprs(i,k)-paprs(i,k+1)) &
            /rg
        END IF
      END DO
    END DO
    DO k = kbmin, khmax
      DO i = 1, klon
        IF (todo(i) .AND. k>=kb(i) .AND. k<=kh(i) .AND. zvapo(i)>0.0) THEN
          d_ql(i, k) = d_ql(i, k) + toliq(i)*zcond(i)/zvapo(i)*max(0.0, -d_q( &
            i,k))
        END IF
      END DO
    END DO
    DO i = 1, klon
      IF (todo(i)) THEN
        zrfl(i) = (1.0-toliq(i))*zcond(i)/dtime
      END IF
    END DO

  ELSE IF (opt_cld==4) THEN

    nexpo = 3
    ! cc         nexpo = 1 ! distribution lineaire

    DO i = 1, klon
      IF (todo(i)) THEN
        zvapo(i) = 0.0 ! quantite integrale de masse (avec ponderation)
      END IF
    END DO
    DO k = kbmin, khmax
      DO i = 1, klon
        IF (todo(i) .AND. k>=(kb(i)+1) .AND. k<=kh(i)) THEN
          zvapo(i) = zvapo(i) + (paprs(i,k)-paprs(i,k+1))/rg*(pplay(i,kb(i))- &
            pplay(i,k))**nexpo
        END IF
      END DO
    END DO
    DO k = kbmin, khmax
      DO i = 1, klon
        IF (todo(i) .AND. k>=(kb(i)+1) .AND. k<=kh(i)) THEN
          d_ql(i, k) = d_ql(i, k) + toliq(i)*zcond(i)/zvapo(i)*(pplay(i,kb(i) &
            )-pplay(i,k))**nexpo
        END IF
      END DO
    END DO
    DO i = 1, klon
      IF (todo(i)) THEN
        zrfl(i) = (1.0-toliq(i))*zcond(i)/dtime
      END IF
    END DO

  ELSE ! valeur non-prevue pour opt_cld

    PRINT *, 'opt_cld est faux:', opt_cld
    CALL abort

  END IF ! fin de opt_cld

  ! L'eau precipitante peut etre re-evaporee:

  IF (evap_prec .AND. kbmax>=2) THEN
    DO k = kbmax, 1, -1
      DO i = 1, klon
        IF (todo(i) .AND. k<=(kb(i)-1) .AND. zrfl(i)>0.0) THEN
          zqev = max(0.0, (zqs(i,k)-q(i,k))*zfrac(i))
          zqevt = coef_eva*(1.0-q(i,k)/zqs(i,k))*sqrt(zrfl(i))* &
            (paprs(i,k)-paprs(i,k+1))/pplay(i, k)*t(i, k)*rd/rg
          zqevt = max(0.0, min(zqevt,zrfl(i)))*rg*dtime/ &
            (paprs(i,k)-paprs(i,k+1))
          zqev = min(zqev, zqevt)
          zrfln = zrfl(i) - zqev*(paprs(i,k)-paprs(i,k+1))/rg/dtime
          d_q(i, k) = -(zrfln-zrfl(i))*(rg/(paprs(i,k)-paprs(i,k+1)))*dtime
          d_t(i, k) = (zrfln-zrfl(i))*(rg/(paprs(i,k)-paprs(i, &
            k+1)))*dtime*rlvtt/rcpd
          zrfl(i) = zrfln
        END IF
      END DO
    END DO
  END IF

  ! La temperature de la premiere couche determine la pluie ou la neige:

  DO i = 1, klon
    IF (todo(i)) THEN
      IF (t(i,1)>rtt) THEN
        rain(i) = rain(i) + zrfl(i)
      ELSE
        snow(i) = snow(i) + zrfl(i)
      END IF
    END IF
  END DO

  RETURN
END SUBROUTINE conkuo
SUBROUTINE kuofcl(pt, pq, pg, pp, ldcum, kcbot)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19940927
  ! adaptation du code de Tiedtke du ECMWF
  ! Objet: calculer le niveau de convection libre
  ! (FCL: Free Convection Level)
  ! ======================================================================
  ! Arguments:
  ! pt---input-R- temperature (K)
  ! pq---input-R- vapeur d'eau (kg/kg)
  ! pg---input-R- geopotentiel (g*z ou z est en metre)
  ! pp---input-R- pression (Pa)

  ! LDCUM---output-L- Y-t-il la convection
  ! kcbot---output-I- Niveau du bas de la convection
  ! ======================================================================
  include "YOMCST.h"
  include "YOETHF.h"

  REAL pt(klon, klev), pq(klon, klev), pg(klon, klev), pp(klon, klev)
  INTEGER kcbot(klon)
  LOGICAL ldcum(klon)

  REAL ztu(klon, klev), zqu(klon, klev), zlu(klon, klev)
  REAL zqold(klon), zbuo
  INTEGER is, i, k

  ! klab=1: on est sous le nuage convectif
  ! klab=2: le bas du nuage convectif
  ! klab=0: autres couches
  INTEGER klab(klon, klev)

  ! quand lflag=.true., on est sous le nuage, il faut donc appliquer
  ! le processus d'elevation.
  LOGICAL lflag(klon)

  DO k = 1, klev
    DO i = 1, klon
      ztu(i, k) = pt(i, k)
      zqu(i, k) = pq(i, k)
      zlu(i, k) = 0.0
      klab(i, k) = 0
    END DO
  END DO
  ! ----------------------------------------------------------------------
  DO i = 1, klon
    klab(i, 1) = 1
    kcbot(i) = 2
    ldcum(i) = .FALSE.
  END DO

  DO k = 2, klev - 1

    is = 0
    DO i = 1, klon
      IF (klab(i,k-1)==1) is = is + 1
      lflag(i) = .FALSE.
      IF (klab(i,k-1)==1) lflag(i) = .TRUE.
    END DO
    IF (is==0) GO TO 290

    ! on eleve le parcel d'air selon l'adiabatique sec

    DO i = 1, klon
      IF (lflag(i)) THEN
        zqu(i, k) = zqu(i, k-1)
        ztu(i, k) = ztu(i, k-1) + (pg(i,k-1)-pg(i,k))/rcpd
        zbuo = ztu(i, k)*(1.+retv*zqu(i,k)) - pt(i, k)*(1.+retv*pq(i,k)) + &
          0.5
        IF (zbuo>0.) klab(i, k) = 1
        zqold(i) = zqu(i, k)
      END IF
    END DO

    ! on calcule la condensation eventuelle

    CALL adjtq(pp(1,k), ztu(1,k), zqu(1,k), lflag, 1)

    ! s'il y a la condensation et la "buoyancy" force est positive
    ! c'est bien le bas de la tour de convection

    DO i = 1, klon
      IF (lflag(i) .AND. zqu(i,k)/=zqold(i)) THEN
        klab(i, k) = 2
        zlu(i, k) = zlu(i, k) + zqold(i) - zqu(i, k)
        zbuo = ztu(i, k)*(1.+retv*zqu(i,k)) - pt(i, k)*(1.+retv*pq(i,k)) + &
          0.5
        IF (zbuo>0.) THEN
          kcbot(i) = k
          ldcum(i) = .TRUE.
        END IF
      END IF
    END DO

290 END DO

  RETURN
END SUBROUTINE kuofcl
SUBROUTINE adjtq(pp, pt, pq, ldflag, kcall)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19940927
  ! adaptation du code de Tiedtke du ECMWF
  ! Objet: ajustement entre T et Q
  ! ======================================================================
  ! Arguments:
  ! pp---input-R- pression (Pa)
  ! pt---input/output-R- temperature (K)
  ! pq---input/output-R- vapeur d'eau (kg/kg)
  ! ======================================================================
  ! TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT

  ! NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
  ! KCALL=0    ENV. T AND QS IN*CUINI*
  ! KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
  ! KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)

  include "YOMCST.h"

  REAL pt(klon), pq(klon), pp(klon)
  LOGICAL ldflag(klon)
  INTEGER kcall

  REAL t_coup
  PARAMETER (t_coup=234.0)

  REAL zcond(klon), zcond1
  REAL zdelta, zcvm5, zldcp, zqsat, zcor, zdqsat
  INTEGER is, i
  include "YOETHF.h"
  include "FCTTRE.h"

  DO i = 1, klon
    zcond(i) = 0.0
  END DO

  DO i = 1, klon
    IF (ldflag(i)) THEN
      zdelta = max(0., sign(1.,rtt-pt(i)))
      zldcp = rlvtt*(1.-zdelta) + zdelta*rlstt
      zldcp = zldcp/rcpd/(1.0+rvtmp2*pq(i))
      IF (thermcep) THEN
        zcvm5 = r5les*rlvtt*(1.-zdelta) + zdelta*r5ies*rlstt
        zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*pq(i))
        zqsat = r2es*foeew(pt(i), zdelta)/pp(i)
        zqsat = min(0.5, zqsat)
        zcor = 1./(1.-retv*zqsat)
        zqsat = zqsat*zcor
        zdqsat = foede(pt(i), zdelta, zcvm5, zqsat, zcor)
      ELSE
        IF (pt(i)<t_coup) THEN
          zqsat = qsats(pt(i))/pp(i)
          zdqsat = dqsats(pt(i), zqsat)
        ELSE
          zqsat = qsatl(pt(i))/pp(i)
          zdqsat = dqsatl(pt(i), zqsat)
        END IF
      END IF
      zcond(i) = (pq(i)-zqsat)/(1.+zdqsat)
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
      zldcp = rlvtt*(1.-zdelta) + zdelta*rlstt
      zldcp = zldcp/rcpd/(1.0+rvtmp2*pq(i))
      IF (thermcep) THEN
        zcvm5 = r5les*rlvtt*(1.-zdelta) + zdelta*r5ies*rlstt
        zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*pq(i))
        zqsat = r2es*foeew(pt(i), zdelta)/pp(i)
        zqsat = min(0.5, zqsat)
        zcor = 1./(1.-retv*zqsat)
        zqsat = zqsat*zcor
        zdqsat = foede(pt(i), zdelta, zcvm5, zqsat, zcor)
      ELSE
        IF (pt(i)<t_coup) THEN
          zqsat = qsats(pt(i))/pp(i)
          zdqsat = dqsats(pt(i), zqsat)
        ELSE
          zqsat = qsatl(pt(i))/pp(i)
          zdqsat = dqsatl(pt(i), zqsat)
        END IF
      END IF
      zcond1 = (pq(i)-zqsat)/(1.+zdqsat)
      pt(i) = pt(i) + zldcp*zcond1
      pq(i) = pq(i) - zcond1
    END IF
  END DO

230 CONTINUE
  RETURN
END SUBROUTINE adjtq
SUBROUTINE fiajh(dtime, paprs, pplay, t, q, d_t, d_q, d_ql, rneb, rain, snow, &
    ibas, itop)
  USE dimphy
  IMPLICIT NONE

  ! Ajustement humide (Schema de convection de Manabe)
  ! .
  include "YOMCST.h"

  ! Arguments:

  REAL dtime ! intervalle du temps (s)
  REAL t(klon, klev) ! temperature (K)
  REAL q(klon, klev) ! humidite specifique (kg/kg)
  REAL paprs(klon, klev+1) ! pression a inter-couche (Pa)
  REAL pplay(klon, klev) ! pression au milieu de couche (Pa)

  REAL d_t(klon, klev) ! incrementation pour la temperature
  REAL d_q(klon, klev) ! incrementation pour vapeur d'eau
  REAL d_ql(klon, klev) ! incrementation pour l'eau liquide
  REAL rneb(klon, klev) ! fraction nuageuse

  REAL rain(klon) ! variable non utilisee
  REAL snow(klon) ! variable non utilisee
  INTEGER ibas(klon) ! variable non utilisee
  INTEGER itop(klon) ! variable non utilisee

  REAL t_coup
  PARAMETER (t_coup=234.0)
  REAL seuil_vap
  PARAMETER (seuil_vap=1.0E-10)

  ! Variables locales:

  INTEGER i, k
  INTEGER k1, k1p, k2, k2p
  LOGICAL itest(klon)
  REAL delta_q(klon, klev)
  REAL cp_new_t(klev)
  REAL cp_delta_t(klev)
  REAL new_qb(klev)
  REAL v_cptj(klev), v_cptjk1, v_ssig
  REAL v_cptt(klon, klev), v_p, v_t
  REAL v_qs(klon, klev), v_qsd(klon, klev)
  REAL zq1(klon), zq2(klon)
  REAL gamcpdz(klon, 2:klev)
  REAL zdp, zdpm

  REAL zsat ! sur-saturation
  REAL zflo ! flotabilite

  REAL local_q(klon, klev), local_t(klon, klev)

  REAL zdelta, zcor, zcvm5

  include "YOETHF.h"
  include "FCTTRE.h"

  DO k = 1, klev
    DO i = 1, klon
      local_q(i, k) = q(i, k)
      local_t(i, k) = t(i, k)
      rneb(i, k) = 0.0
      d_ql(i, k) = 0.0
      d_t(i, k) = 0.0
      d_q(i, k) = 0.0
    END DO
  END DO
  DO i = 1, klon
    rain(i) = 0.0
    snow(i) = 0.0
    ibas(i) = 0
    itop(i) = 0
  END DO

  ! Calculer v_qs et v_qsd:

  DO k = 1, klev
    DO i = 1, klon
      v_cptt(i, k) = rcpd*local_t(i, k)
      v_t = local_t(i, k)
      v_p = pplay(i, k)

      IF (thermcep) THEN
        zdelta = max(0., sign(1.,rtt-v_t))
        zcvm5 = r5les*rlvtt*(1.-zdelta) + zdelta*r5ies*rlstt
        zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*local_q(i,k))
        v_qs(i, k) = r2es*foeew(v_t, zdelta)/v_p
        v_qs(i, k) = min(0.5, v_qs(i,k))
        zcor = 1./(1.-retv*v_qs(i,k))
        v_qs(i, k) = v_qs(i, k)*zcor
        v_qsd(i, k) = foede(v_t, zdelta, zcvm5, v_qs(i,k), zcor)
      ELSE
        IF (v_t<t_coup) THEN
          v_qs(i, k) = qsats(v_t)/v_p
          v_qsd(i, k) = dqsats(v_t, v_qs(i,k))
        ELSE
          v_qs(i, k) = qsatl(v_t)/v_p
          v_qsd(i, k) = dqsatl(v_t, v_qs(i,k))
        END IF
      END IF
    END DO
  END DO

  ! Calculer Gamma * Cp * dz: (gamm est le gradient critique)

  DO k = 2, klev
    DO i = 1, klon
      zdp = paprs(i, k) - paprs(i, k+1)
      zdpm = paprs(i, k-1) - paprs(i, k)
      gamcpdz(i, k) = ((rd/rcpd/(zdpm+zdp)*(v_cptt(i,k-1)*zdpm+ &
        v_cptt(i,k)*zdp)+rlvtt/(zdpm+zdp)*(v_qs(i,k-1)*zdpm+ &
        v_qs(i,k)*zdp))*(pplay(i,k-1)-pplay(i,k))/paprs(i,k))/(1.0+(v_qsd(i, &
        k-1)*zdpm+v_qsd(i,k)*zdp)/(zdpm+zdp))
    END DO
  END DO

  ! ------------------------------------ modification des profils instables
  DO i = 1, klon
    itest(i) = .FALSE.

    k1 = 0
    k2 = 1

810 CONTINUE ! chercher k1, le bas de la colonne
    k2 = k2 + 1
    IF (k2>klev) GO TO 9999
    zflo = v_cptt(i, k2-1) - v_cptt(i, k2) - gamcpdz(i, k2)
    zsat = (local_q(i,k2-1)-v_qs(i,k2-1))*(paprs(i,k2-1)-paprs(i,k2)) + &
      (local_q(i,k2)-v_qs(i,k2))*(paprs(i,k2)-paprs(i,k2+1))
    IF (zflo<=0.0 .OR. zsat<=0.0) GO TO 810
    k1 = k2 - 1
    itest(i) = .TRUE.

820 CONTINUE ! chercher k2, le haut de la colonne
    IF (k2==klev) GO TO 821
    k2p = k2 + 1
    zsat = zsat + (paprs(i,k2p)-paprs(i,k2p+1))*(local_q(i,k2p)-v_qs(i,k2p))
    zflo = v_cptt(i, k2p-1) - v_cptt(i, k2p) - gamcpdz(i, k2p)
    IF (zflo<=0.0 .OR. zsat<=0.0) GO TO 821
    k2 = k2p
    GO TO 820
821 CONTINUE

    ! ------------------------------------------------------ ajustement local
830 CONTINUE ! ajustement proprement dit
    v_cptj(k1) = 0.0
    zdp = paprs(i, k1) - paprs(i, k1+1)
    v_cptjk1 = ((1.0+v_qsd(i,k1))*(v_cptt(i,k1)+v_cptj(k1))+rlvtt*(local_q(i, &
      k1)-v_qs(i,k1)))*zdp
    v_ssig = zdp*(1.0+v_qsd(i,k1))

    k1p = k1 + 1
    DO k = k1p, k2
      zdp = paprs(i, k) - paprs(i, k+1)
      v_cptj(k) = v_cptj(k-1) + gamcpdz(i, k)
      v_cptjk1 = v_cptjk1 + zdp*((1.0+v_qsd(i,k))*(v_cptt(i, &
        k)+v_cptj(k))+rlvtt*(local_q(i,k)-v_qs(i,k)))
      v_ssig = v_ssig + zdp*(1.0+v_qsd(i,k))
    END DO

    DO k = k1, k2
      cp_new_t(k) = v_cptjk1/v_ssig - v_cptj(k)
      cp_delta_t(k) = cp_new_t(k) - v_cptt(i, k)
      new_qb(k) = v_qs(i, k) + v_qsd(i, k)*cp_delta_t(k)/rlvtt
      local_q(i, k) = new_qb(k)
      local_t(i, k) = cp_new_t(k)/rcpd
    END DO

    ! --------------------------------------------------- sondage vers le bas
    ! -- on redefinit les variables prognostiques dans
    ! -- la colonne qui vient d'etre ajustee

    DO k = k1, k2
      v_cptt(i, k) = rcpd*local_t(i, k)
      v_t = local_t(i, k)
      v_p = pplay(i, k)

      IF (thermcep) THEN
        zdelta = max(0., sign(1.,rtt-v_t))
        zcvm5 = r5les*rlvtt*(1.-zdelta) + zdelta*r5ies*rlstt
        zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*local_q(i,k))
        v_qs(i, k) = r2es*foeew(v_t, zdelta)/v_p
        v_qs(i, k) = min(0.5, v_qs(i,k))
        zcor = 1./(1.-retv*v_qs(i,k))
        v_qs(i, k) = v_qs(i, k)*zcor
        v_qsd(i, k) = foede(v_t, zdelta, zcvm5, v_qs(i,k), zcor)
      ELSE
        IF (v_t<t_coup) THEN
          v_qs(i, k) = qsats(v_t)/v_p
          v_qsd(i, k) = dqsats(v_t, v_qs(i,k))
        ELSE
          v_qs(i, k) = qsatl(v_t)/v_p
          v_qsd(i, k) = dqsatl(v_t, v_qs(i,k))
        END IF
      END IF
    END DO
    DO k = 2, klev
      zdpm = paprs(i, k-1) - paprs(i, k)
      zdp = paprs(i, k) - paprs(i, k+1)
      gamcpdz(i, k) = ((rd/rcpd/(zdpm+zdp)*(v_cptt(i,k-1)*zdpm+ &
        v_cptt(i,k)*zdp)+rlvtt/(zdpm+zdp)*(v_qs(i,k-1)*zdpm+ &
        v_qs(i,k)*zdp))*(pplay(i,k-1)-pplay(i,k))/paprs(i,k))/(1.0+(v_qsd(i, &
        k-1)*zdpm+v_qsd(i,k)*zdp)/(zdpm+zdp))
    END DO

    ! Verifier si l'on peut etendre la colonne vers le bas

    IF (k1==1) GO TO 841 ! extension echouee
    zflo = v_cptt(i, k1-1) - v_cptt(i, k1) - gamcpdz(i, k1)
    zsat = (local_q(i,k1-1)-v_qs(i,k1-1))*(paprs(i,k1-1)-paprs(i,k1)) + &
      (local_q(i,k1)-v_qs(i,k1))*(paprs(i,k1)-paprs(i,k1+1))
    IF (zflo<=0.0 .OR. zsat<=0.0) GO TO 841 ! extension echouee

840 CONTINUE
    k1 = k1 - 1
    IF (k1==1) GO TO 830 ! GOTO 820 (a tester, Z.X.Li, mars 1995)
    zsat = zsat + (local_q(i,k1-1)-v_qs(i,k1-1))*(paprs(i,k1-1)-paprs(i,k1))
    zflo = v_cptt(i, k1-1) - v_cptt(i, k1) - gamcpdz(i, k1)
    IF (zflo>0.0 .AND. zsat>0.0) THEN
      GO TO 840
    ELSE
      GO TO 830 ! GOTO 820 (a tester, Z.X.Li, mars 1995)
    END IF
841 CONTINUE

    GO TO 810 ! chercher d'autres blocks en haut

9999 END DO ! boucle sur tous les points
  ! -----------------------------------------------------------------------

  ! Determiner la fraction nuageuse (hypothese: la nebulosite a lieu
  ! a l'endroit ou la vapeur d'eau est diminuee par l'ajustement):

  DO k = 1, klev
    DO i = 1, klon
      IF (itest(i)) THEN
        delta_q(i, k) = local_q(i, k) - q(i, k)
        IF (delta_q(i,k)<0.) rneb(i, k) = 1.0
      END IF
    END DO
  END DO

  ! Distribuer l'eau condensee en eau liquide nuageuse (hypothese:
  ! l'eau liquide est distribuee aux endroits ou la vapeur d'eau
  ! diminue et d'une maniere proportionnelle a cet diminution):

  DO i = 1, klon
    IF (itest(i)) THEN
      zq1(i) = 0.0
      zq2(i) = 0.0
    END IF
  END DO
  DO k = 1, klev
    DO i = 1, klon
      IF (itest(i)) THEN
        zdp = paprs(i, k) - paprs(i, k+1)
        zq1(i) = zq1(i) - delta_q(i, k)*zdp
        zq2(i) = zq2(i) - min(0.0, delta_q(i,k))*zdp
      END IF
    END DO
  END DO
  DO k = 1, klev
    DO i = 1, klon
      IF (itest(i)) THEN
        IF (zq2(i)/=0.0) d_ql(i, k) = -min(0.0, delta_q(i,k))*zq1(i)/zq2(i)
      END IF
    END DO
  END DO

  DO k = 1, klev
    DO i = 1, klon
      local_q(i, k) = max(local_q(i,k), seuil_vap)
    END DO
  END DO

  DO k = 1, klev
    DO i = 1, klon
      d_t(i, k) = local_t(i, k) - t(i, k)
      d_q(i, k) = local_q(i, k) - q(i, k)
    END DO
  END DO

  RETURN
END SUBROUTINE fiajh
SUBROUTINE fiajc(dtime, paprs, pplay, t, q, conv_q, d_t, d_q, d_ql, rneb, &
    rain, snow, ibas, itop)
  USE dimphy
  IMPLICIT NONE

  include "YOMCST.h"

  ! Options:

  INTEGER plb ! niveau de depart pour la convection
  PARAMETER (plb=4)

  ! Mystere: cette option n'est pas innocente pour les resultats !
  ! Qui peut resoudre ce mystere ? (Z.X.Li mars 1995)
  LOGICAL vector ! calcul vectorise
  PARAMETER (vector=.FALSE.)

  REAL t_coup
  PARAMETER (t_coup=234.0)

  ! Arguments:

  REAL q(klon, klev) ! humidite specifique (kg/kg)
  REAL t(klon, klev) ! temperature (K)
  REAL paprs(klon, klev+1) ! pression a inter-couche (Pa)
  REAL pplay(klon, klev) ! pression au milieu de couche (Pa)
  REAL dtime ! intervalle du temps (s)
  REAL conv_q(klon, klev) ! taux de convergence de l'humidite
  REAL rneb(klon, klev) ! fraction nuageuse
  REAL d_q(klon, klev) ! incrementaion pour la vapeur d'eau
  REAL d_ql(klon, klev) ! incrementation pour l'eau liquide
  REAL d_t(klon, klev) ! incrementation pour la temperature
  REAL rain(klon) ! variable non-utilisee
  REAL snow(klon) ! variable non-utilisee
  INTEGER itop(klon) ! variable non-utilisee
  INTEGER ibas(klon) ! variable non-utilisee

  INTEGER kh(klon), i, k
  LOGICAL nuage(klon), test(klon, klev)
  REAL zconv(klon), zdeh(klon, klev), zvirt(klon)
  REAL zdqs(klon, klev), zqs(klon, klev)
  REAL ztt, zvar, zfrac(klon)
  REAL zq1(klon), zq2(klon)
  REAL zdelta, zcor, zcvm5

  include "YOETHF.h"
  include "FCTTRE.h"

  ! Initialiser les sorties:

  DO k = 1, klev
    DO i = 1, klon
      rneb(i, k) = 0.0
      d_ql(i, k) = 0.0
      d_t(i, k) = 0.0
      d_q(i, k) = 0.0
    END DO
  END DO
  DO i = 1, klon
    itop(i) = 0
    ibas(i) = 0
    rain(i) = 0.0
    snow(i) = 0.0
  END DO

  ! Calculer Qs et L/Cp * dQs/dT:

  DO k = 1, klev
    DO i = 1, klon
      ztt = t(i, k)
      IF (thermcep) THEN
        zdelta = max(0., sign(1.,rtt-ztt))
        zcvm5 = r5les*rlvtt*(1.-zdelta) + zdelta*r5ies*rlstt
        zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*q(i,k))
        zqs(i, k) = r2es*foeew(ztt, zdelta)/pplay(i, k)
        zqs(i, k) = min(0.5, zqs(i,k))
        zcor = 1./(1.-retv*zqs(i,k))
        zqs(i, k) = zqs(i, k)*zcor
        zdqs(i, k) = foede(ztt, zdelta, zcvm5, zqs(i,k), zcor)
      ELSE
        IF (ztt<t_coup) THEN
          zqs(i, k) = qsats(ztt)/pplay(i, k)
          zdqs(i, k) = dqsats(ztt, zqs(i,k))
        ELSE
          zqs(i, k) = qsatl(ztt)/pplay(i, k)
          zdqs(i, k) = dqsatl(ztt, zqs(i,k))
        END IF
      END IF
    END DO
  END DO

  ! Determiner la difference de l'energie totale saturee:

  DO i = 1, klon
    k = plb
    zdeh(i, k) = rcpd*(t(i,k-1)-t(i,k)) - rd*0.5*(t(i,k-1)+t(i,k))/paprs(i, k &
      )*(pplay(i,k-1)-pplay(i,k)) + rlvtt*(zqs(i,k-1)-zqs(i,k))
    zdeh(i, k) = zdeh(i, k)*0.5 ! on prend la moitie
  END DO
  DO k = plb + 1, klev
    DO i = 1, klon
      zdeh(i, k) = zdeh(i, k-1) + rcpd*(t(i,k-1)-t(i,k)) - &
        rd*0.5*(t(i,k-1)+t(i,k))/paprs(i, k)*(pplay(i,k-1)-pplay(i,k)) + &
        rlvtt*(zqs(i,k-1)-zqs(i,k))
    END DO
  END DO

  ! Determiner le sommet du nuage selon l'instabilite
  ! Calculer les convergences d'humidite (reelle et virtuelle)

  DO i = 1, klon
    nuage(i) = .TRUE.
    zconv(i) = 0.0
    zvirt(i) = 0.0
    kh(i) = -999
  END DO
  DO k = plb, klev
    DO i = 1, klon
      nuage(i) = nuage(i) .AND. zdeh(i, k) > 0.0
      IF (nuage(i)) THEN
        kh(i) = k
        zconv(i) = zconv(i) + conv_q(i, k)*dtime*(paprs(i,k)-paprs(i,k+1))
        zvirt(i) = zvirt(i) + (zdeh(i,k)/rlvtt+zqs(i,k)-q(i,k))*(paprs(i,k)- &
          paprs(i,k+1))
      END IF
    END DO
  END DO

  IF (vector) THEN


    DO k = plb, klev
      DO i = 1, klon
        IF (k<=kh(i) .AND. kh(i)>plb .AND. zconv(i)>0.0) THEN
          test(i, k) = .TRUE.
          zfrac(i) = max(0.0, min(zconv(i)/zvirt(i),1.0))
        ELSE
          test(i, k) = .FALSE.
        END IF
      END DO
    END DO

    DO k = plb, klev
      DO i = 1, klon
        IF (test(i,k)) THEN
          zvar = zdeh(i, k)/(1.0+zdqs(i,k))
          d_q(i, k) = (zvar*zdqs(i,k)/rlvtt+zqs(i,k)-q(i,k))*zfrac(i) - &
            conv_q(i, k)*dtime
          d_t(i, k) = zvar*zfrac(i)/rcpd
        END IF
      END DO
    END DO

    DO i = 1, klon
      zq1(i) = 0.0
      zq2(i) = 0.0
    END DO
    DO k = plb, klev
      DO i = 1, klon
        IF (test(i,k)) THEN
          IF (d_q(i,k)<0.0) rneb(i, k) = zfrac(i)
          zq1(i) = zq1(i) - d_q(i, k)*(paprs(i,k)-paprs(i,k+1))
          zq2(i) = zq2(i) - min(0.0, d_q(i,k))*(paprs(i,k)-paprs(i,k+1))
        END IF
      END DO
    END DO

    DO k = plb, klev
      DO i = 1, klon
        IF (test(i,k)) THEN
          IF (zq2(i)/=0.) d_ql(i, k) = -min(0.0, d_q(i,k))*zq1(i)/zq2(i)
        END IF
      END DO
    END DO

  ELSE ! (.NOT. vector)

    DO i = 1, klon
      IF (kh(i)>plb .AND. zconv(i)>0.0) THEN
        ! cc         IF (kh(i).LE.plb) GOTO 999 ! il n'y a pas d'instabilite
        ! cc         IF (zconv(i).LE.0.0) GOTO 999 ! convergence insuffisante
        zfrac(i) = max(0.0, min(zconv(i)/zvirt(i),1.0))
        DO k = plb, kh(i)
          zvar = zdeh(i, k)/(1.0+zdqs(i,k))
          d_q(i, k) = (zvar*zdqs(i,k)/rlvtt+zqs(i,k)-q(i,k))*zfrac(i) - &
            conv_q(i, k)*dtime
          d_t(i, k) = zvar*zfrac(i)/rcpd
        END DO

        zq1(i) = 0.0
        zq2(i) = 0.0
        DO k = plb, kh(i)
          IF (d_q(i,k)<0.0) rneb(i, k) = zfrac(i)
          zq1(i) = zq1(i) - d_q(i, k)*(paprs(i,k)-paprs(i,k+1))
          zq2(i) = zq2(i) - min(0.0, d_q(i,k))*(paprs(i,k)-paprs(i,k+1))
        END DO
        DO k = plb, kh(i)
          IF (zq2(i)/=0.) d_ql(i, k) = -min(0.0, d_q(i,k))*zq1(i)/zq2(i)
        END DO
      END IF
    END DO

  END IF ! fin de teste sur vector

  RETURN
END SUBROUTINE fiajc
