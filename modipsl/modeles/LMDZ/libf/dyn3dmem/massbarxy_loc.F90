SUBROUTINE massbarxy_loc(masse,massebxy)
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr. Hourdin.
!-------------------------------------------------------------------------------
! Purpose: Compute air mass mean along X and Y in each cell.
! See iniconst for more details.
  USE parallel_lmdz
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: masse   (ijb_u:ije_u,llm)
  REAL, INTENT(OUT) :: massebxy(ijb_v:ije_v,llm)
!===============================================================================
! Local variables:
  INTEGER :: ij, l, ijb, ije
!===============================================================================
  ijb=ij_begin-iip1
  ije=ij_end
  IF(pole_nord) ijb=ijb+iip1
  IF(pole_sud)  ije=ije-iip1
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
  DO l=1,llm
    DO ij=ijb,ije-1
      massebxy(ij,l)=masse(ij     ,l)*alpha2(ij     ) + &
     +               masse(ij+1   ,l)*alpha3(ij+1   ) + &
     +               masse(ij+iip1,l)*alpha1(ij+iip1) + &
     +               masse(ij+iip2,l)*alpha4(ij+iip2)
    END DO
    DO ij=ijb+iip1-1,ije+iip1-1,iip1; massebxy(ij,l)=massebxy(ij-iim,l); END DO
  END DO
!$OMP END DO NOWAIT

END SUBROUTINE massbarxy_loc

