      SUBROUTINE pression_loc( ngrid, ap, bp, ps, p )
      USE parallel_lmdz, ONLY: ij_begin, ij_end, ijb_u, ije_u,
     &                         pole_nord, pole_sud, omp_chunk
c

c      Auteurs : P. Le Van , Fr.Hourdin  .

c  ************************************************************************
c     Calcule la pression p(l) aux differents niveaux l = 1 ( niveau du
c     sol) a l = llm +1 ,ces niveaux correspondant aux interfaces des (llm) 
c     couches , avec  p(ij,llm +1) = 0.  et p(ij,1) = ps(ij)  .      
c  ************************************************************************
c
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
c
      INTEGER,INTENT(IN) :: ngrid ! not used
      INTEGER l,ij
 
      REAL,INTENT(IN) :: ap( llmp1 ), bp( llmp1 ), ps( ijb_u:ije_u )
      REAL,INTENT(OUT) :: p( ijb_u:ije_u,llmp1 ) 
      
      INTEGER ijb,ije

      
      ijb=ij_begin-iip1
      ije=ij_end+2*iip1
      
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO    l    = 1, llmp1
        DO  ij   = ijb, ije
         p(ij,l) = ap(l) + bp(l) * ps(ij)
        ENDDO
      ENDDO
c$OMP END DO NOWAIT   
      RETURN
      END
