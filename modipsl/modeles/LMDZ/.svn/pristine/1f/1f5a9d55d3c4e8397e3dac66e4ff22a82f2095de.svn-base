SUBROUTINE massbar_loc(masse,massebx,masseby)
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
  REAL, INTENT(IN)  :: masse  (ijb_u:ije_u,llm)
  REAL, INTENT(OUT) :: massebx(ijb_u:ije_u,llm)
  REAL, INTENT(OUT) :: masseby(ijb_v:ije_v,llm)
!-------------------------------------------------------------------------------
! Method used. Each scalar point is associated to 4 area coefficients:
!    * alpha1(i,j) at point ( i+1/4,j-1/4 )
!    * alpha2(i,j) at point ( i+1/4,j+1/4 )
!    * alpha3(i,j) at point ( i-1/4,j+1/4 )
!    * alpha4(i,j) at point ( i-1/4,j-1/4 )
! where alpha1(i,j) = aire(i+1/4,j-1/4)/ aire(i,j)
!
!   alpha4 .         . alpha1    . alpha4
!    (i,j)             (i,j)       (i+1,j)
!
!             P .        U .          . P
!           (i,j)       (i,j)         (i+1,j)
!
!   alpha3 .         . alpha2    .alpha3 
!    (i,j)              (i,j)     (i+1,j)
!
!             V .        Z .          . V
!           (i,j)
!
!   alpha4 .         . alpha1    .alpha4
!   (i,j+1)            (i,j+1)   (i+1,j+1) 
!
!             P .        U .          . P
!          (i,j+1)                    (i+1,j+1)
!
!
!    massebx(i,j) = masse(i  ,j) * ( alpha1(i  ,j) + alpha2(i,j))   +
!                   masse(i+1,j) * ( alpha3(i+1,j) + alpha4(i+1,j) )
!     localized at point  ... U (i,j) ...
!
!    masseby(i,j) = masse(i,j  ) * ( alpha2(i,j  ) + alpha3(i,j  )  +
!                   masse(i,j+1) * ( alpha1(i,j+1) + alpha4(i,j+1)  
!     localized at point  ... V (i,j) ...
!===============================================================================
! Local variables:
  INTEGER :: ij, l, ijb, ije
!===============================================================================
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)  
  DO l=1,llm
    ijb=ij_begin
    ije=ij_end+iip1
    IF(pole_sud) ije=ije-iip1
    DO ij=ijb,ije-1
      massebx(ij,l)=masse(ij,l)*alpha1p2(ij)+masse(ij+1   ,l)*alpha3p4(ij+1)
    END DO
    DO ij=ijb+iim,ije+iim,iip1; massebx(ij,l)=massebx(ij-iim,l); END DO
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF(pole_nord) ijb=ij_begin
    IF(pole_sud) ije=ij_end-iip1
    DO ij=ijb,ije
      masseby(ij,l)=masse(ij,l)*alpha2p3(ij)+masse(ij+iip1,l)*alpha1p4(ij+iip1)
    END DO
  END DO
!$OMP END DO NOWAIT

END SUBROUTINE massbar_loc

