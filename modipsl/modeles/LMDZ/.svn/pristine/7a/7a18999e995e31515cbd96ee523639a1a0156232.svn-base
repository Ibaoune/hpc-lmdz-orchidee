FUNCTION qcheck(klon,klev,paprs,q,ql,aire)
  IMPLICIT none
  !
  ! Calculer et imprimer l'eau totale. A utiliser pour verifier
  ! la conservation de l'eau
  !
  include "YOMCST.h"
  INTEGER,INTENT(IN) :: klon,klev
  REAL,INTENT(IN) :: paprs(klon,klev+1), q(klon,klev), ql(klon,klev)
  REAL,INTENT(IN) :: aire(klon)
  REAL qtotal, zx, qcheck
  INTEGER i, k
  !
  zx = 0.0
  DO i = 1, klon
     zx = zx + aire(i)
  ENDDO
  qtotal = 0.0
  DO k = 1, klev
     DO i = 1, klon
        qtotal = qtotal + (q(i,k)+ql(i,k)) * aire(i) &
             *(paprs(i,k)-paprs(i,k+1))/RG
     ENDDO
  ENDDO
  !
  qcheck = qtotal/zx
  !
END FUNCTION qcheck

