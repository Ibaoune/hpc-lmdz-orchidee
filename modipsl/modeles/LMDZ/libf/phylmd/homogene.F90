
! $Header$

SUBROUTINE homogene(paprs, q, dq, u, v, du, dv)
  USE dimphy
  IMPLICIT NONE
  ! ==============================================================
  ! Schema ad hoc du melange vertical pour les vitesses u et v,
  ! a appliquer apres le schema de convection (fiajc et fiajh).

  ! paprs:input, pression demi-couche (inter-couche)
  ! q:    input, vapeur d'eau (kg/kg)
  ! dq:   input, incrementation de vapeur d'eau (de la convection)
  ! u:    input, vitesse u
  ! v:    input, vitesse v

  ! du:   output, incrementation pour u
  ! dv:   output, incrementation pour v
  ! ==============================================================

  REAL paprs(klon, klev+1)
  REAL q(klon, klev), dq(klon, klev)
  REAL u(klon, klev), du(klon, klev)
  REAL v(klon, klev), dv(klon, klev)

  REAL zm_dq(klon) ! quantite totale de l'eau deplacee
  REAL zm_q(klon) ! quantite totale de la vapeur d'eau
  REAL zm_u(klon) ! moyenne de u (brassage parfait et total)
  REAL zm_v(klon) ! moyenne de v (brassage parfait et total)
  REAL z_frac(klon) ! fraction du brassage parfait et total
  REAL zm_dp(klon)

  REAL zx
  INTEGER i, k
  REAL frac_max
  PARAMETER (frac_max=0.1)
  REAL seuil
  PARAMETER (seuil=1.0E-10)
  LOGICAL faisrien
  PARAMETER (faisrien=.TRUE.)

  DO k = 1, klev
    DO i = 1, klon
      du(i, k) = 0.0
      dv(i, k) = 0.0
    END DO
  END DO

  IF (faisrien) RETURN

  DO i = 1, klon
    zm_dq(i) = 0.
    zm_q(i) = 0.
    zm_u(i) = 0.
    zm_v(i) = 0.
    zm_dp(i) = 0.
  END DO
  DO k = 1, klev
    DO i = 1, klon
      IF (abs(dq(i,k))>seuil) THEN
        zx = paprs(i, k) - paprs(i, k+1)
        zm_dq(i) = zm_dq(i) + abs(dq(i,k))*zx
        zm_q(i) = zm_q(i) + q(i, k)*zx
        zm_dp(i) = zm_dp(i) + zx
        zm_u(i) = zm_u(i) + u(i, k)*zx
        zm_v(i) = zm_v(i) + v(i, k)*zx
      END IF
    END DO
  END DO

  ! Hypothese principale: apres la convection, la vitesse de chaque
  ! couche est composee de deux parties: celle (1-z_frac) de la vitesse
  ! original et celle (z_frac) de la vitesse moyenne qui serait la
  ! vitesse de chaque couche si le brassage etait parfait et total.
  ! La fraction du brassage est calculee par le rapport entre la quantite
  ! totale de la vapeur d'eau deplacee (ou condensee) et la quantite
  ! totale de la vapeur d'eau. Et cette fraction est limitee a frac_max
  ! (Est-ce vraiment raisonnable ? Z.X. Li, le 07-09-1995).

  DO i = 1, klon
    IF (zm_dp(i)>=1.0E-15 .AND. zm_q(i)>=1.0E-15) THEN
      z_frac(i) = min(frac_max, zm_dq(i)/zm_q(i))
      zm_u(i) = zm_u(i)/zm_dp(i)
      zm_v(i) = zm_v(i)/zm_dp(i)
    END IF
  END DO
  DO k = 1, klev
    DO i = 1, klon
      IF (zm_dp(i)>=1.E-15 .AND. zm_q(i)>=1.E-15 .AND. abs(dq(i, &
          k))>seuil) THEN
        du(i, k) = u(i, k)*(1.-z_frac(i)) + zm_u(i)*z_frac(i) - u(i, k)
        dv(i, k) = v(i, k)*(1.-z_frac(i)) + zm_v(i)*z_frac(i) - v(i, k)
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE homogene
