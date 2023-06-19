
! $Header$

SUBROUTINE ecribins(unit, pz)
  USE dimphy
  IMPLICIT NONE
  ! -----------------------------------------------------------------------
  include "dimensions.h"
  ! ccc#include "dimphy.h"
  include "paramet.h"
  include "comgeom.h"

  ! arguments:
  ! ----------
  INTEGER unit
  REAL pz(klon)

  ! local:
  ! ------
  INTEGER i, j, ig
  REAL zz(iim+1, jjm+1)
  ! -----------------------------------------------------------------------
  ! passage a la grille dynamique:
  ! ------------------------------
  DO i = 1, iim + 1
    zz(i, 1) = pz(1)
    zz(i, jjm+1) = pz(klon)
  END DO
  ! traitement des point normaux
  DO j = 2, jjm
    ig = 2 + (j-2)*iim
    CALL scopy(iim, pz(ig), 1, zz(1,j), 1)
    zz(iim+1, j) = zz(1, j)
  END DO
  ! -----------------------------------------------------------------------
#ifdef VPP
  CALL ecriture(unit, zz, (iim+1)*(jjm+1))
#else
  WRITE (unit) zz
#endif


  RETURN
END SUBROUTINE ecribins
SUBROUTINE ecribina(unit, pz)
  USE dimphy
  IMPLICIT NONE
  ! -----------------------------------------------------------------------
  include "dimensions.h"
  ! ccc#include "dimphy.h"
  include "paramet.h"
  include "comgeom.h"

  ! arguments:
  ! ----------
  INTEGER unit
  REAL pz(klon, klev)

  ! local:
  ! ------
  INTEGER i, j, ilay, ig
  REAL zz(iim+1, jjm+1, llm)
  ! -----------------------------------------------------------------------
  ! passage a la grille dynamique:
  ! ------------------------------
  DO ilay = 1, llm
    ! traitement des poles
    DO i = 1, iim + 1
      zz(i, 1, ilay) = pz(1, ilay)
      zz(i, jjm+1, ilay) = pz(klon, ilay)
    END DO
    ! traitement des point normaux
    DO j = 2, jjm
      ig = 2 + (j-2)*iim
      CALL scopy(iim, pz(ig,ilay), 1, zz(1,j,ilay), 1)
      zz(iim+1, j, ilay) = zz(1, j, ilay)
    END DO
  END DO
  ! -----------------------------------------------------------------------
  DO ilay = 1, llm
#ifdef VPP
    CALL ecriture(unit, zz(1,1,ilay), (iim+1)*(jjm+1))
#else
    WRITE (unit)((zz(i,j,ilay),i=1,iim+1), j=1, jjm+1)
#endif
  END DO

  RETURN
END SUBROUTINE ecribina
#ifdef VPP
@OPTIONS NODOUBLE
SUBROUTINE ecriture(nunit, r8, n)
  INTEGER nunit, n, i
  REAL (KIND=8) r8(n)
  REAL r4(n)

  DO i = 1, n
    r4(i) = r8(i)
  END DO
  WRITE (nunit) r4
  RETURN
END SUBROUTINE ecriture
#endif
