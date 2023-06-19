SUBROUTINE cv3_buoy(nloc, ncum, nd, icb, inb, pbase, plcl, p, ph, ale, cin, &
    tv, tvp, buoy)
  ! **************************************************************
  ! *
  ! CV3_BUOY                                                    *
  ! Buoyancy corrections to account for ALE             *
  ! *
  ! written by   : MOREAU Cecile, 07/08/2003, 15.55.48          *
  ! modified by :                                               *
  ! **************************************************************

  IMPLICIT NONE

  include "cvthermo.h"
  include "cv3param.h"
  include "YOMCST2.h"

  ! input:
  INTEGER ncum, nd, nloc
  INTEGER icb(nloc), inb(nloc)
  REAL pbase(nloc), plcl(nloc)
  REAL p(nloc, nd), ph(nloc, nd+1)
  REAL ale(nloc), cin(nloc)
  REAL tv(nloc, nd), tvp(nloc, nd)

  ! output:
  REAL buoy(nloc, nd)

  ! local variables:
  INTEGER il, k
  INTEGER kmx(nloc)
  REAL bll(nloc), bmx(nloc)
  REAL gamma(nloc)
  LOGICAL ok(nloc)

  REAL dgamma
  REAL buoymin
  PARAMETER (dgamma=2.E-03) !dgamma gamma
  PARAMETER (buoymin=2.)

  LOGICAL fixed_bll
  SAVE fixed_bll
  DATA fixed_bll/.TRUE./
  !$OMP THREADPRIVATE(fixed_bll)


  ! print *,' Ale+cin ',ale(1)+cin(1)
  ! --------------------------------------------------------------
  ! Recompute buoyancies
  ! --------------------------------------------------------------
  DO k = 1, nl
    DO il = 1, ncum
      buoy(il, k) = tvp(il, k) - tv(il, k)
    END DO
  END DO

  ! -------------------------------------------------------------
  ! -- Compute low level buoyancy ( function of Ale+Cin )
  ! -------------------------------------------------------------
  IF (fixed_bll) THEN

    DO il = 1, ncum
      bll(il) = 0.5
    END DO
  ELSE

    DO il = 1, ncum
      IF (ale(il)+cin(il)>0.) THEN
        gamma(il) = 4.*buoy(il, icb(il))**2 + 8.*dgamma*(ale(il)+cin(il))*tv( &
          il, icb(il))/grav
        gamma(il) = max(gamma(il), 1.E-10)
      END IF
    END DO

    DO il = 1, ncum
      IF (ale(il)+cin(il)>0.) THEN
        bll(il) = 4.*dgamma*(ale(il)+cin(il))*tv(il, icb(il))/ &
          (grav*(abs(buoy(il,icb(il))+0.5*sqrt(gamma(il)))))
      END IF
    END DO

    DO il = 1, ncum
      IF (ale(il)+cin(il)>0.) THEN
        bll(il) = min(bll(il), buoymin)
      END IF
    END DO

  END IF !(fixed_bll)


  ! -------------------------------------------------------------
  ! --Get highest buoyancy among levels below LCL-200hPa
  ! -------------------------------------------------------------

  DO il = 1, ncum
    bmx(il) = -1000.
    kmx(il) = icb(il)
    ok(il) = .TRUE.
  END DO

  DO k = 1, nl
    DO il = 1, ncum
      IF (ale(il)+cin(il)>0. .AND. ok(il)) THEN
        IF (k>icb(il) .AND. k<=inb(il)) THEN
          ! c         print *,'k,p(il,k),plcl(il)-200. ',
          ! k,p(il,k),plcl(il)-200.
          IF (p(il,k)>plcl(il)-200.) THEN
            IF (buoy(il,k)>bmx(il)) THEN
              bmx(il) = buoy(il, k)
              kmx(il) = k
              IF (bmx(il)>=bll(il)) ok(il) = .FALSE.
            END IF
          END IF
        END IF
      END IF
    END DO
  END DO

  ! print *,' ==cv3_buoy== bll(1),bmx(1),icb(1),kmx(1) '
  ! $       ,bll(1),bmx(1),icb(1),kmx(1)

  ! -------------------------------------------------------------
  ! --Calculate modified buoyancies
  ! -------------------------------------------------------------

  DO il = 1, ncum
    IF (ale(il)+cin(il)>0.) THEN
      bll(il) = min(bll(il), bmx(il))
    END IF
  END DO

  DO k = 1, nl
    DO il = 1, ncum
      IF (ale(il)+cin(il)>0.) THEN
        IF (k>=icb(il) .AND. k<=kmx(il)-1) THEN
          buoy(il, k) = bll(il)
        END IF
      END IF
    END DO
  END DO

!CR:Correction of buoy for what comes next
!keep flag or to modify in all cases?
  IF (iflag_mix_adiab.eq.1) THEN
  DO k = 1, nl 
    DO il = 1, ncum
       IF ((k>=kmx(il)) .AND. (k<=inb(il)) .AND. (buoy(il,k).lt.0.)) THEN
          buoy(il,k)=buoy(il,k-1)
       END IF
    ENDDO
  ENDDO
  ENDIF

  RETURN
END SUBROUTINE cv3_buoy
