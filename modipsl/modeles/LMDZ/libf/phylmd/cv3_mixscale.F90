SUBROUTINE cv3_mixscale(nloc, ncum, na, ment, m)
  ! **************************************************************
  ! *
  ! CV3_MIXSCALE                                                *
  ! *
  ! *
  ! written by   : Jean-Yves Grandpeix, 30/05/2003, 16.34.37    *
  ! modified by :                                               *
  ! **************************************************************

  IMPLICIT NONE

  include "cv3param.h"

!inputs:
  INTEGER, INTENT (IN)                               :: ncum, na, nloc
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: m
!input/outputs:
  REAL, DIMENSION (nloc, na, na), INTENT (INOUT)     :: ment

!local variables:
  INTEGER i, j, il

    DO j = 1, nl
      DO i = 1, nl
        DO il = 1, ncum
          ment(il, i, j) = m(il, i)*ment(il, i, j)
        END DO
      END DO
    END DO


  RETURN
END SUBROUTINE cv3_mixscale
