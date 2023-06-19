      SUBROUTINE covcont_loc (klevel,ucov, vcov, ucont, vcont )
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c  *********************************************************************
c    calcul des compos. contravariantes a partir des comp.covariantes
c  ********************************************************************
c
c=======================================================================

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

      INTEGER klevel
      REAL ucov( ijb_u:ije_u,klevel ),  vcov( ijb_v:ije_v,klevel )
      REAL ucont( ijb_u:ije_u,klevel ), vcont( ijb_v:ije_v,klevel )
      INTEGER   l,ij
      INTEGER ijbu,ijbv,ijeu,ijev

      
      ijbu=ij_begin-iip1
      ijbv=ij_begin-iip1
      ijeu=ij_end+iip1
      ijev=ij_end+iip1
      
      if (pole_nord) then 
        ijbu=ij_begin+iip1
        ijbv=ij_begin
      endif
      
      if (pole_sud) then
        ijeu=ij_end-iip1
        ijev=ij_end-iip1
      endif

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
      DO 10 l = 1,klevel

      DO 2  ij = ijb_u,ije_u
      ucont( ij,l ) = ucov( ij,l ) * unscu2( ij )
   2  CONTINUE

      DO 4 ij = ijb_v,ije_v
      vcont( ij,l ) = vcov( ij,l ) * unscv2( ij )
   4  CONTINUE

  10  CONTINUE
c$OMP END DO NOWAIT
      RETURN
      END
