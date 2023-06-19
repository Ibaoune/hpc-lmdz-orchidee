
! $Header$

! ======================================================================
SUBROUTINE orbite(xjour, longi, dist)
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) (adapte du GCM du LMD) date: 19930818
  ! Objet: pour un jour donne, calculer la longitude vraie de la terre
  ! (par rapport au point vernal-21 mars) dans son orbite solaire
  ! calculer aussi la distance terre-soleil (unite astronomique)
  ! ======================================================================
  ! Arguments:
  ! xjour--INPUT--R- jour de l'annee a compter du 1er janvier
  ! longi--OUTPUT-R- longitude vraie en degres par rapport au point
  ! vernal (21 mars) en degres
  ! dist---OUTPUT-R- distance terre-soleil (par rapport a la moyenne)
  REAL xjour, longi, dist
  ! ======================================================================
  include "YOMCST.h"

  ! -- Variables dynamiques locales
  REAL pir, xl, xllp, xee, xse, xlam, dlamm, anm, ranm, anv, ranv

  pir = 4.0*atan(1.0)/180.0
  xl = r_peri + 180.0
  xllp = xl*pir
  xee = r_ecc*r_ecc
  xse = sqrt(1.0-xee)
  xlam = (r_ecc/2.0+r_ecc*xee/8.0)*(1.0+xse)*sin(xllp) - &
    xee/4.0*(0.5+xse)*sin(2.0*xllp) + r_ecc*xee/8.0*(1.0/3.0+xse)*sin(3.0* &
    xllp)
  xlam = 2.0*xlam/pir
  dlamm = xlam + (xjour-81.0)
  anm = dlamm - xl
  ranm = anm*pir
  xee = xee*r_ecc
  ranv = ranm + (2.0*r_ecc-xee/4.0)*sin(ranm) + 5.0/4.0*r_ecc*r_ecc*sin(2.0* &
    ranm) + 13.0/12.0*xee*sin(3.0*ranm)

  anv = ranv/pir
  longi = anv + xl

  dist = (1-r_ecc*r_ecc)/(1+r_ecc*cos(pir*(longi-(r_peri+180.0))))
  RETURN
END SUBROUTINE orbite
! ======================================================================
SUBROUTINE angle(longi, lati, frac, muzero)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: Calculer la duree d'ensoleillement pour un jour et la hauteur
  ! du soleil (cosinus de l'angle zinithal) moyenne sur la journee
  ! ======================================================================
  ! Arguments:
  ! longi----INPUT-R- la longitude vraie de la terre dans son plan
  ! solaire a partir de l'equinoxe de printemps (degre)
  ! lati-----INPUT-R- la latitude d'un point sur la terre (degre)
  ! frac-----OUTPUT-R la duree d'ensoleillement dans la journee divisee
  ! par 24 heures (unite en fraction de 0 a 1)
  ! muzero---OUTPUT-R la moyenne du cosinus de l'angle zinithal sur
  ! la journee (0 a 1)
  ! ======================================================================
  REAL longi
  REAL lati(klon), frac(klon), muzero(klon)
  include "YOMCST.h"
  REAL lat, omega, lon_sun, lat_sun
  REAL pi_local, incl
  INTEGER i

  pi_local = 4.0*atan(1.0)
  incl = r_incl*pi_local/180.

  lon_sun = longi*pi_local/180.0
  lat_sun = asin(sin(lon_sun)*sin(incl))

  DO i = 1, klon
    lat = lati(i)*pi_local/180.0

    IF (lat>=(pi_local/2.+lat_sun) .OR. lat<=(-pi_local/2.+lat_sun)) THEN
      omega = 0.0 ! nuit polaire
    ELSE IF (lat>=(pi_local/2.-lat_sun) .OR. lat<=(-pi_local/2.-lat_sun)) &
        THEN
      omega = pi_local ! journee polaire
    ELSE
      omega = -tan(lat)*tan(lat_sun)
      omega = acos(omega)
    END IF

    frac(i) = omega/pi_local

    IF (omega>0.0) THEN
      muzero(i) = sin(lat)*sin(lat_sun) + cos(lat)*cos(lat_sun)*sin(omega)/ &
        omega
    ELSE
      muzero(i) = 0.0
    END IF
  END DO

  RETURN
END SUBROUTINE angle
! ====================================================================
SUBROUTINE zenang(longi, gmtime, pdtrad1, pdtrad2, lat, long, pmu0, frac)
  USE dimphy
  IMPLICIT NONE
  ! =============================================================
  ! Auteur : O. Boucher (LMD/CNRS)
  ! d'apres les routines zenith et angle de Z.X. Li
  ! Objet  : calculer les valeurs moyennes du cos de l'angle zenithal
  ! et l'ensoleillement moyen entre gmtime1 et gmtime2
  ! connaissant la declinaison, la latitude et la longitude.
  ! Rque   : Different de la routine angle en ce sens que zenang
  ! fournit des moyennes de pmu0 et non des valeurs
  ! instantanees, du coup frac prend toutes les valeurs
  ! entre 0 et 1. La routine integre entre gmtime+pdtrad1 et 
  ! gmtime+pdtrad2 avec pdtrad1 et pdtrad2 exprimes en secondes.
  ! Date   : premiere version le 13 decembre 1994
  ! revu pour  GCM  le 30 septembre 1996
  ! revu le 3 septembre 2015 pour les bornes de l'integrale
  ! ===============================================================
  ! longi : la longitude vraie de la terre dans son plan
  ! solaire a partir de l'equinoxe de printemps (degre)
  ! gmtime : temps universel en fraction de jour
  ! pdtrad1 : borne inferieure du pas de temps du rayonnement (secondes)
  ! pdtrad2 : borne inferieure du pas de temps du rayonnement (secondes)
  ! pdtrad2-pdtrad1 correspond a pdtrad, le pas de temps du rayonnement (secondes)
  ! lat------INPUT : latitude en degres
  ! long-----INPUT : longitude en degres
  ! pmu0-----OUTPUT: angle zenithal moyen entre gmtime+pdtrad1 et gmtime+pdtrad2
  ! frac-----OUTPUT: ensoleillement moyen entre gmtime+pdtrad1 et gmtime+pdtrad2
  ! ================================================================
  include "YOMCST.h"
  ! ================================================================
  REAL, INTENT (IN) :: longi, gmtime, pdtrad1, pdtrad2
  REAL lat(klon), long(klon), pmu0(klon), frac(klon)
  ! ================================================================
  INTEGER i
  REAL gmtime1, gmtime2
  REAL pi_local, deux_pi_local, incl
  REAL omega1, omega2, omega
  ! omega1, omega2 : temps 1 et 2 exprime en radian avec 0 a midi.
  ! omega : heure en radian du coucher de soleil
  ! -omega est donc l'heure en radian de lever du soleil
  REAL omegadeb, omegafin
  REAL zfrac1, zfrac2, z1_mu, z2_mu
  REAL lat_sun ! declinaison en radian
  REAL lon_sun ! longitude solaire en radian
  REAL latr    ! latitude du pt de grille en radian
  ! ================================================================

  pi_local = 4.0*atan(1.0)
  deux_pi_local = 2.0*pi_local
  incl = r_incl*pi_local/180.

  lon_sun = longi*pi_local/180.0
  lat_sun = asin(sin(lon_sun)*sin(incl))

  gmtime1 = gmtime*86400. + pdtrad1
  gmtime2 = gmtime*86400. + pdtrad2

  DO i = 1, klon

    latr = lat(i)*pi_local/180.

    omega = 0.0 !--nuit polaire

    IF (latr>=(pi_local/2.-lat_sun) .OR. latr<=(-pi_local/2.-lat_sun)) THEN
      omega = pi_local ! journee polaire
    END IF

    IF (latr<(pi_local/2.+lat_sun) .AND. latr>(-pi_local/2.+lat_sun) .AND. &
        latr<(pi_local/2.-lat_sun) .AND. latr>(-pi_local/2.-lat_sun)) THEN
      omega = -tan(latr)*tan(lat_sun)
      omega = acos(omega)
    END IF

    omega1 = gmtime1 + long(i)*86400.0/360.0
    omega1 = omega1/86400.0*deux_pi_local
    omega1 = mod(omega1+deux_pi_local, deux_pi_local)
    omega1 = omega1 - pi_local

    omega2 = gmtime2 + long(i)*86400.0/360.0
    omega2 = omega2/86400.0*deux_pi_local
    omega2 = mod(omega2+deux_pi_local, deux_pi_local)
    omega2 = omega2 - pi_local

    IF (omega1<=omega2) THEN !--on est dans la meme journee locale

      IF (omega2<=-omega .OR. omega1>=omega .OR. omega<1E-5) THEN !--nuit
        frac(i) = 0.0
        pmu0(i) = 0.0
      ELSE !--jour+nuit/jour
        omegadeb = max(-omega, omega1)
        omegafin = min(omega, omega2)
        frac(i) = (omegafin-omegadeb)/(omega2-omega1)
        pmu0(i) = sin(latr)*sin(lat_sun) + cos(latr)*cos(lat_sun)*(sin( &
          omegafin)-sin(omegadeb))/(omegafin-omegadeb)
      END IF

    ELSE !---omega1 GT omega2 -- a cheval sur deux journees

      ! -------------------entre omega1 et pi
      IF (omega1>=omega) THEN !--nuit
        zfrac1 = 0.0
        z1_mu = 0.0
      ELSE !--jour+nuit
        omegadeb = max(-omega, omega1)
        omegafin = omega
        zfrac1 = omegafin - omegadeb
        z1_mu = sin(latr)*sin(lat_sun) + cos(latr)*cos(lat_sun)*(sin(omegafin &
          )-sin(omegadeb))/(omegafin-omegadeb)
      END IF
      ! ---------------------entre -pi et omega2
      IF (omega2<=-omega) THEN !--nuit
        zfrac2 = 0.0
        z2_mu = 0.0
      ELSE !--jour+nuit
        omegadeb = -omega
        omegafin = min(omega, omega2)
        zfrac2 = omegafin - omegadeb
        z2_mu = sin(latr)*sin(lat_sun) + cos(latr)*cos(lat_sun)*(sin(omegafin &
          )-sin(omegadeb))/(omegafin-omegadeb)

      END IF
      ! -----------------------moyenne
      frac(i) = (zfrac1+zfrac2)/(omega2+deux_pi_local-omega1)
      pmu0(i) = (zfrac1*z1_mu+zfrac2*z2_mu)/max(zfrac1+zfrac2, 1.E-10)

    END IF !---comparaison omega1 et omega2

  END DO

END SUBROUTINE zenang
! ===================================================================
SUBROUTINE zenith(longi, gmtime, lat, long, pmu0, fract)
  USE dimphy
  IMPLICIT NONE

  ! Auteur(s): Z.X. Li (LMD/ENS)

  ! Objet: calculer le cosinus de l'angle zenithal du soleil en
  ! connaissant la declinaison du soleil, la latitude et la
  ! longitude du point sur la terre, et le temps universel

  ! Arguments d'entree:
  ! longi  : declinaison du soleil (en degres)
  ! gmtime : temps universel en second qui varie entre 0 et 86400
  ! lat    : latitude en degres
  ! long   : longitude en degres
  ! Arguments de sortie:
  ! pmu0   : cosinus de l'angle zenithal

  ! ====================================================================
  include "YOMCST.h"
  ! ====================================================================
  REAL longi, gmtime
  REAL lat(klon), long(klon), pmu0(klon), fract(klon)
  ! =====================================================================
  INTEGER n
  REAL zpi, zpir, omega, zgmtime
  REAL incl, lat_sun, lon_sun
  ! ----------------------------------------------------------------------
  zpi = 4.0*atan(1.0)
  zpir = zpi/180.0
  zgmtime = gmtime*86400.

  incl = r_incl*zpir

  lon_sun = longi*zpir
  lat_sun = asin(sin(lon_sun)*sin(incl))

  ! --initialisation a la nuit

  DO n = 1, klon
    pmu0(n) = 0.
    fract(n) = 0.0
  END DO

  ! 1 degre en longitude = 240 secondes en temps

  DO n = 1, klon
    omega = zgmtime + long(n)*86400.0/360.0
    omega = omega/86400.0*2.0*zpi
    omega = mod(omega+2.0*zpi, 2.0*zpi)
    omega = omega - zpi
    pmu0(n) = sin(lat(n)*zpir)*sin(lat_sun) + cos(lat(n)*zpir)*cos(lat_sun)* &
      cos(omega)
    pmu0(n) = max(pmu0(n), 0.0)
    IF (pmu0(n)>1.E-6) fract(n) = 1.0
  END DO

  RETURN
END SUBROUTINE zenith
