!
! $Header$
!
      SUBROUTINE covnat_loc(klevel,ucov, vcov, unat, vnat )
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  F Hourdin Phu LeVan
c   -------
c
c   Objet:
c   ------
c
c  *********************************************************************
c    calcul des compos. naturelles a partir des comp.covariantes
c  ********************************************************************
c
c=======================================================================

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

      INTEGER klevel
      REAL ucov( ijb_u:ije_u,klevel ),  vcov( ijb_v:ije_v,klevel )
      REAL unat( ijb_u:ije_u,klevel ), vnat( ijb_v:ije_v,klevel )
      INTEGER   l,ij
      INTEGER :: ijb,ije
      
      
      ijb=ij_begin
      ije=ij_end
      
      if (pole_nord) then

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
        DO l = 1,klevel
           DO ij = 1, iip1
              unat (ij,l) =0.
           END DO
        ENDDO
!$OMP ENDDO NOWAIT
      endif

      if (pole_sud) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
        DO l = 1,klevel
           DO ij = ip1jm+1, ip1jmp1  
            unat (ij,l) =0.
           END DO
        ENDDO
!$OMP ENDDO NOWAIT
      endif

      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l = 1,klevel
         DO ij = ijb, ije
            unat( ij,l ) = ucov( ij,l ) / cu(ij)
         ENDDO
      END DO
!$OMP ENDDO NOWAIT

      ijb=ij_begin-iip1
      ije=ij_end
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end-iip1
     
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l = 1,klevel
         DO ij = ijb,ije
            vnat( ij,l ) = vcov( ij,l ) / cv(ij)
         ENDDO
      ENDDO
!$OMP ENDDO NOWAIT
      
      RETURN
      END
