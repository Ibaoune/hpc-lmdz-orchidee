
! $Header$

SUBROUTINE flxtr(pdtime, pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, pt, pplay, &
    paprs, kcbot, kctop, kdtop, x, dx)
  USE dimphy
  IMPLICIT NONE
  ! =====================================================================
  ! Objet : Melange convectif de traceurs a partir des flux de masse
  ! Date : 13/12/1996 -- 13/01/97
  ! Auteur: O. Boucher (LOA) sur inspiration de Z. X. Li (LMD),
  ! Brinkop et Sausen (1996) et Boucher et al. (1996).
  ! ATTENTION : meme si cette routine se veut la plus generale possible,
  ! elle a herite de certaines notations et conventions du
  ! schema de Tiedtke (1993).
  ! --En particulier, les couches sont numerotees de haut en bas !!!
  ! Ceci est valable pour les flux, kcbot, kctop et kdtop
  ! mais pas pour les entrees x, pplay, paprs !!!!
  ! --Un schema amont est choisi pour calculer les flux pour s'assurer
  ! de la positivite des valeurs de traceurs, cela implique des eqs
  ! differentes pour les flux de traceurs montants et descendants.
  ! --pmfu est positif, pmfd est negatif
  ! --Tous les flux d'entrainements et de detrainements sont positifs
  ! contrairement au schema de Tiedtke d'ou les changements de signe!!!!
  ! =====================================================================

  include "YOMCST.h"
  include "YOECUMF.h"

  REAL pdtime
  ! --les flux sont definis au 1/2 niveaux
  ! --pmfu(klev+1) et pmfd(klev+1) sont implicitement nuls
  REAL pmfu(klon, klev) ! flux de masse dans le panache montant
  REAL pmfd(klon, klev) ! flux de masse dans le panache descendant
  REAL pen_u(klon, klev) ! flux entraine dans le panache montant
  REAL pde_u(klon, klev) ! flux detraine dans le panache montant
  REAL pen_d(klon, klev) ! flux entraine dans le panache descendant
  REAL pde_d(klon, klev) ! flux detraine dans le panache descendant
  ! --idem mais en variables locales
  REAL zpen_u(klon, klev)
  REAL zpde_u(klon, klev)
  REAL zpen_d(klon, klev)
  REAL zpde_d(klon, klev)

  REAL pplay(klon, klev) ! pression aux couches (bas en haut)
  REAL pap(klon, klev) ! pression aux couches (haut en bas)
  REAL pt(klon, klev) ! temperature aux couches (bas en haut)
  REAL zt(klon, klev) ! temperature aux couches (haut en bas)
  REAL paprs(klon, klev+1) ! pression aux 1/2 couches (bas en haut)
  REAL paph(klon, klev+1) ! pression aux 1/2 couches (haut en bas)
  INTEGER kcbot(klon) ! niveau de base de la convection
  INTEGER kctop(klon) ! niveau du sommet de la convection +1
  INTEGER kdtop(klon) ! niveau de sommet du panache descendant
  REAL x(klon, klev) ! q de traceur (bas en haut)
  REAL zx(klon, klev) ! q de traceur (haut en bas)
  REAL dx(klon, klev) ! tendance de traceur  (bas en haut)

  ! --variables locales
  ! --les flux de x sont definis aux 1/2 niveaux
  ! --xu et xd sont definis aux niveaux complets
  REAL xu(klon, klev) ! q de traceurs dans le panache montant
  REAL xd(klon, klev) ! q de traceurs dans le panache descendant
  REAL xe(klon, klev) ! q de traceurs dans l'environnement
  REAL zmfux(klon, klev+1) ! flux de x dans le panache montant
  REAL zmfdx(klon, klev+1) ! flux de x dans le panache descendant
  REAL zmfex(klon, klev+1) ! flux de x dans l'environnement
  INTEGER i, k
  REAL zmfmin
  PARAMETER (zmfmin=1.E-10)

  ! On remet les taux d'entrainement et de detrainement dans le panache
  ! descendant a des valeurs positives.
  ! On ajuste les valeurs de pen_u, pen_d pde_u et pde_d pour que la
  ! conservation de la masse soit realisee a chaque niveau dans les 2
  ! panaches.
  DO k = 1, klev
    DO i = 1, klon
      zpen_u(i, k) = pen_u(i, k)
      zpde_u(i, k) = pde_u(i, k)
    END DO
  END DO

  DO k = 1, klev - 1
    DO i = 1, klon
      zpen_d(i, k) = -pen_d(i, k+1)
      zpde_d(i, k) = -pde_d(i, k+1)
    END DO
  END DO

  DO i = 1, klon
    zpen_d(i, klev) = 0.0
    zpde_d(i, klev) = -pmfd(i, klev)
    ! Correction 03 11 97
    ! zpen_d(i,kdtop(i)-1) = pmfd(i,kdtop(i)-1)-pmfd(i,kdtop(i))
    IF (kdtop(i)==klev+1) THEN
      zpen_d(i, kdtop(i)-1) = pmfd(i, kdtop(i)-1)
    ELSE
      zpen_d(i, kdtop(i)-1) = pmfd(i, kdtop(i)-1) - pmfd(i, kdtop(i))
    END IF

    zpde_u(i, kctop(i)-2) = pmfu(i, kctop(i)-1)
    zpen_u(i, klev) = pmfu(i, klev)
  END DO

  DO i = 1, klon
    DO k = kcbot(i), klev - 1
      zpen_u(i, k) = pmfu(i, k) - pmfu(i, k+1)
    END DO
  END DO

  ! conversion des sens de notations bas-haut et haut-bas

  DO k = 1, klev + 1
    DO i = 1, klon
      paph(i, klev+2-k) = paprs(i, k)
    END DO
  END DO

  DO i = 1, klon
    DO k = 1, klev
      pap(i, klev+1-k) = pplay(i, k)
      zt(i, klev+1-k) = pt(i, k)
      zx(i, klev+1-k) = x(i, k)
    END DO
  END DO

  ! --initialisations des flux de traceurs aux extremites de la colonne

  DO i = 1, klon
    zmfux(i, klev+1) = 0.0
    zmfdx(i, 1) = 0.0
    zmfex(i, 1) = 0.0
  END DO

  ! --calcul des flux dans le panache montant

  DO k = klev, 1, -1
    DO i = 1, klon
      IF (k>=kcbot(i)) THEN
        xu(i, k) = zx(i, k)
        zmfux(i, k) = pmfu(i, k)*xu(i, k)
      ELSE
        zmfux(i, k) = (zmfux(i,k+1)+zpen_u(i,k)*zx(i,k))/ &
          (1.+zpde_u(i,k)/max(zmfmin,pmfu(i,k)))
        xu(i, k) = zmfux(i, k)/max(zmfmin, pmfu(i,k))
      END IF
    END DO
  END DO

  ! --calcul des flux dans le panache descendant

  DO k = 1, klev - 1
    DO i = 1, klon
      IF (k<=kdtop(i)-1) THEN
        xd(i, k) = (zx(i,k)+xu(i,k))/2.
        zmfdx(i, k+1) = pmfd(i, k+1)*xd(i, k)
      ELSE
        zmfdx(i, k+1) = (zmfdx(i,k)-zpen_d(i,k)*zx(i,k))/ &
          (1.-zpde_d(i,k)/min(-zmfmin,pmfd(i,k+1)))
        xd(i, k) = zmfdx(i, k+1)/min(-zmfmin, pmfd(i,k+1))
      END IF
    END DO
  END DO
  DO i = 1, klon
    zmfdx(i, klev+1) = 0.0
    xd(i, klev) = (zpen_d(i,klev)*zx(i,klev)-zmfdx(i,klev))/ &
      max(zmfmin, zpde_d(i,klev))
  END DO

  ! --introduction du flux de retour dans l'environnement

  DO k = 1, klev - 1
    DO i = 1, klon
      IF (k<=kctop(i)-3) THEN
        xe(i, k) = zx(i, k)
        zmfex(i, k+1) = -(pmfu(i,k+1)+pmfd(i,k+1))*xe(i, k)
      ELSE
        zmfex(i, k+1) = (zmfex(i,k)-(zpde_u(i,k)*xu(i,k)+zpde_d(i,k)*xd(i, &
          k)))/(1.-(zpen_d(i,k)+zpen_u(i,k))/min(-zmfmin,-pmfu(i,k+1)-pmfd(i, &
          k+1)))
        xe(i, k) = zmfex(i, k+1)/min(-zmfmin, -pmfu(i,k+1)-pmfd(i,k+1))
      END IF
    END DO
  END DO
  DO i = 1, klon
    zmfex(i, klev+1) = 0.0
    xe(i, klev) = (zpde_u(i,klev)*xu(i,klev)+zpde_d(i,klev)*xd(i,klev)-zmfex( &
      i,klev))/max(zmfmin, zpen_u(i,klev)+zpen_d(i,klev))
  END DO

  ! --calcul final des tendances

  DO k = 1, klev
    DO i = 1, klon
      dx(i, klev+1-k) = rg/(paph(i,k+1)-paph(i,k))*pdtime* &
        (zmfux(i,k+1)-zmfux(i,k)+zmfdx(i,k+1)-zmfdx(i,k)+zmfex(i,k+1)- &
        zmfex(i,k))
    END DO
  END DO

  RETURN
END SUBROUTINE flxtr
