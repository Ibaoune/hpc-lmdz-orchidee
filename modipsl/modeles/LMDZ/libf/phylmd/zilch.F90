
! $Header$

SUBROUTINE zilch(x, m)

  ! Zero the real array x dimensioned m.

  IMPLICIT NONE

  INTEGER m, i
  REAL x(m)

  DO i = 1, m
    x(i) = 0.0
  END DO
  RETURN
END SUBROUTINE zilch
