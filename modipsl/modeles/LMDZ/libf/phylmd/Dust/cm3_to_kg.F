      SUBROUTINE cm3_to_kg(pplay,t_seri,tr_seri)

      USE dimphy
      USE infotrac
      USE indice_sol_mod

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
!JE20150707      RD = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 28.9644
      DO k = 1, klev
      DO i = 1, klon
        zrho=pplay(i,k)/t_seri(i,k)/RD
        tr_seri(i,k)=tr_seri(i,k)*1.e6/zrho
      ENDDO
      ENDDO
c 
      END
