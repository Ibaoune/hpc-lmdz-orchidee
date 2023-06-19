      SUBROUTINE kg_to_cm3(pplay,t_seri,tr_seri)
c      
      USE dimphy
      USE infotrac
      IMPLICIT NONE
c
#include "dimensions.h"
#include "YOMCST.h"
c      
      REAL t_seri(klon,klev), pplay(klon,klev)
      REAL tr_seri(klon,klev)
      REAL zrho
      INTEGER i, k
c
      DO k = 1, klev
      DO i = 1, klon
        zrho=pplay(i,k)/t_seri(i,k)/RD
        tr_seri(i,k)=tr_seri(i,k)/1.e6*zrho
      ENDDO
      ENDDO
c 
      END
