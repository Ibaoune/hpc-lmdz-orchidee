
! $Header$

SUBROUTINE ecriregs(unit, pz)
  USE dimphy
  IMPLICIT NONE
  ! -----------------------------------------------------------------------
  include "dimensions.h"
  ! ccc#include "dimphy.h"
  include "paramet.h"
  include "comgeom.h"
  include "regdim.h"

  ! arguments:
  ! ----------
  INTEGER unit
  REAL pz(klon)

  ! local:
  ! ------
  INTEGER i, j, ig
  REAL zz(iim, jjm+1)
  INTEGER nleng
  PARAMETER (nleng=(i2_fin-i2_deb+1+i1_fin-i1_deb+1)*(j_fin-j_deb+1))
  REAL zzz(nleng)

  ! -----------------------------------------------------------------------
  ! passage a la grille dynamique:
  ! ------------------------------
  DO i = 1, iim
    zz(i, 1) = pz(1)
    zz(i, jjm+1) = pz(klon)
  END DO

  ! traitement des point normaux
  DO j = 2, jjm
    ig = 2 + (j-2)*iim
    CALL scopy(iim, pz(ig), 1, zz(1,j), 1)
  END DO
  ! -----------------------------------------------------------------------
  ig = 0
  DO j = j_deb, j_fin
    DO i = i1_deb, i1_fin
      ig = ig + 1
      zzz(ig) = zz(i, j)
    END DO
    DO i = i2_deb, i2_fin
      ig = ig + 1
      zzz(ig) = zz(i, j)
    END DO
  END DO
#ifdef VPP
  CALL ecriture(unit, zzz, nleng)
#else
  WRITE (unit) zzz
#endif
  RETURN
END SUBROUTINE ecriregs
SUBROUTINE ecrirega(unit, pz)
  USE dimphy
  IMPLICIT NONE
  ! -----------------------------------------------------------------------
  include "dimensions.h"
  ! ccc#include "dimphy.h"
  include "paramet.h"
  include "comgeom.h"
  include "regdim.h"

  ! arguments:
  ! ----------
  INTEGER unit
  REAL pz(klon, klev)

  ! local:
  ! ------
  INTEGER i, j, ilay, ig
  REAL zz(iim, jjm+1, llm)
  INTEGER nleng
  PARAMETER (nleng=(i2_fin-i2_deb+1+i1_fin-i1_deb+1)*(j_fin-j_deb+1))
  REAL zzz(nleng)
  ! -----------------------------------------------------------------------
  ! passage a la grille dynamique:
  ! ------------------------------
  DO ilay = 1, llm
    ! traitement des poles
    DO i = 1, iim
      zz(i, 1, ilay) = pz(1, ilay)
      zz(i, jjm+1, ilay) = pz(klon, ilay)
    END DO
    ! traitement des point normaux
    DO j = 2, jjm
      ig = 2 + (j-2)*iim
      CALL scopy(iim, pz(ig,ilay), 1, zz(1,j,ilay), 1)
    END DO
  END DO
  ! -----------------------------------------------------------------------
  DO ilay = 1, llm
    ig = 0
    DO j = j_deb, j_fin
      DO i = i1_deb, i1_fin
        ig = ig + 1
        zzz(ig) = zz(i, j, ilay)
      END DO
      DO i = i2_deb, i2_fin
        ig = ig + 1
        zzz(ig) = zz(i, j, ilay)
      END DO
    END DO
#ifdef VPP
    CALL ecriture(unit, zzz, nleng)
#else
    WRITE (unit) zzz
#endif
  END DO

  RETURN
END SUBROUTINE ecrirega
