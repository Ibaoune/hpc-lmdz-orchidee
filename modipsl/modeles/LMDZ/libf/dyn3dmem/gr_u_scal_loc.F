!
! $Header$
!
      SUBROUTINE gr_u_scal_loc(nx,x_u,x_scal)
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
      REAL x_u(ijb_u:ije_u,nx),x_scal(ijb_u:ije_u,nx)

c   Local:
c   ------

      INTEGER l,ij
      INTEGER :: ijb,ije

c-----------------------------------------------------------------------
      ijb=ij_begin
      ije=ij_end

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l=1,nx
         DO ij=ijb+1,ije
            x_scal(ij,l)=
     s      (aireu(ij)*x_u(ij,l)+aireu(ij-1)*x_u(ij-1,l))
     s      /(aireu(ij)+aireu(ij-1))
         ENDDO
      ENDDO
!$OMP ENDDO NOWAIT

      ijb=ij_begin
      ije=ij_end

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l=1,nx
         DO ij=ijb,ije-iip1+1,iip1
	   x_scal(ij,l)=x_scal(ij+iip1-1,l)
	 ENDDO
      ENDDO
!$OMP ENDDO NOWAIT
      RETURN
      
      END
