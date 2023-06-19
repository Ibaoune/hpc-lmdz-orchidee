SUBROUTINE cv3_crit(nloc, ncum, nd, icb, inb, p, ph, pzero, v, threshold, &
    kcrit, pcrit)
  ! **************************************************************
  ! *
  ! CV3_CRIT   Find pressure level where vertical profile of    *
  ! variable 'v' intersects 'threshold'              *
  ! *
  ! written by   : FROHWIRTH Julie, 13/08/2003, 21.55.12        *
  ! modified by :                                               *
  ! **************************************************************

  IMPLICIT NONE

  include "cv3param.h"

  ! input:
  INTEGER ncum, nd, nloc
  INTEGER icb(nloc), inb(nloc)
  REAL p(nloc, nd), ph(nloc, nd+1)
  REAL pzero(nloc)
  REAL v(nloc, nd), threshold

  ! output:
  INTEGER kcrit(nloc)
  REAL pcrit(nloc)

  ! local variables
  INTEGER i, j, k, il
  LOGICAL ok(nloc)

  DO il = 1, ncum
    ok(il) = .TRUE.
    pcrit(il) = -1.
    kcrit(il) = 0
  END DO

  DO i = 1, nl
    DO il = 1, ncum
      IF (i>icb(il) .AND. i<=inb(il)) THEN
        IF (p(il,i)<=pzero(il) .AND. ok(il)) THEN
          IF ((v(il,i)-threshold)*(v(il,i-1)-threshold)<0.) THEN
            pcrit(il) = ((threshold-v(il,i))*p(il,i-1)-(threshold-v(il, &
              i-1))*p(il,i))/(v(il,i-1)-v(il,i))
            IF (pcrit(il)>pzero(il)) THEN
              pcrit(il) = -1.
            ELSE
              ok(il) = .FALSE.
              kcrit(il) = i
              IF (pcrit(il)<ph(il,i)) kcrit(il) = kcrit(il) + 1
            END IF
          END IF ! end IF (v(i) ...
        END IF ! end IF (P(i) ...
      END IF ! end IF (icb+1 le i le inb)
    END DO
  END DO


  RETURN
END SUBROUTINE cv3_crit
