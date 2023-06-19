
! $Header$

SUBROUTINE haut2bas(klon, klev, varb2h, varh2b)
  IMPLICIT NONE

  INTEGER klon, klev
  REAL varb2h(klon, klev), varh2b(klon, klev)
  INTEGER i, k, kinv

  DO k = 1, klev
    kinv = klev - k + 1
    DO i = 1, klon
      varh2b(i, k) = varb2h(i, kinv)
    END DO
  END DO

  RETURN
END SUBROUTINE haut2bas
