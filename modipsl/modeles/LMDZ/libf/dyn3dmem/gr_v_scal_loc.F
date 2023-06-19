!
! $Header$
!
      SUBROUTINE gr_v_scal_loc(nx,x_v,x_scal)
c%W%    %G%
c=======================================================================
c
c   Author:    Frederic Hourdin      original: 11/11/92
c   -------
c
c   Subject:
c   ------
c
c   Method:
c   --------
c
c   Interface:
c   ----------
c
c      Input:
c      ------
c
c      Output:
c      -------
c
c=======================================================================
      USE parallel_lmdz
      IMPLICIT NONE
c-----------------------------------------------------------------------
c   Declararations:
c   ---------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

c   Arguments:
c   ----------

      INTEGER nx
      REAL x_v(ijb_v:ije_v,nx),x_scal(ijb_v:ije_v,nx)

c   Local:
c   ------

      INTEGER l,ij
      INTEGER :: ijb,ije
c-----------------------------------------------------------------------
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l=1,nx
         DO ij=ijb,ije
            x_scal(ij,l)=
     s      (airev(ij-iip1)*x_v(ij-iip1,l)+airev(ij)*x_v(ij,l))
     s      /(airev(ij-iip1)+airev(ij))
         ENDDO
      ENDDO
!$OMP ENDDO NOWAIT
      
      if (pole_nord) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
        DO l=1,nx
           DO ij=1,iip1
              x_scal(ij,l)=0.
           ENDDO
        ENDDO
!$OMP ENDDO NOWAIT
      endif
    
      if (pole_sud) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
        DO l=1,nx
           DO ij=ip1jm+1,ip1jmp1
              x_scal(ij,l)=0.
           ENDDO
        ENDDO
!$OMP ENDDO NOWAIT
      endif

      RETURN
      END
