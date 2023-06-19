SUBROUTINE tourpot_loc ( vcov, ucov, massebxy, vorpot )
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van.
!-------------------------------------------------------------------------------
! Purpose: Compute potential vorticity.
  USE parallel_lmdz
  USE mod_filtreg_p
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: vcov    (ijb_v:ije_v,llm)
  REAL, INTENT(IN)  :: ucov    (ijb_u:ije_u,llm)
  REAL, INTENT(IN)  :: massebxy(ijb_v:ije_v,llm)
  REAL, INTENT(OUT) :: vorpot  (ijb_v:ije_v,llm)
!===============================================================================
! Method used:
!   vorpot = ( Filtre( d(vcov)/dx - d(ucov)/dy ) + fext ) /psbarxy
!===============================================================================
! Local variables:
  INTEGER :: l, ij, ije, ijb, jje, jjb
  REAL    :: rot(ijb_v:ije_v,llm)
!===============================================================================

  ijb=ij_begin-iip1
  ije=ij_end
  IF(pole_nord) ijb=ij_begin

!--- Wind vorticity ; correction: rot(iip1,j,l) = rot(1,j,l)
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    IF(pole_sud) ije=ij_end-iip1-1
    DO ij=ijb,ije
      rot(ij,l)=vcov(ij+1,l)-vcov(ij,l)+ucov(ij+iip1,l)-ucov(ij,l)
    END DO
    IF(pole_sud) ije=ij_end-iip1
    DO ij=ijb+iip1-1,ije,iip1; rot(ij,l)=rot(ij-iim,l); END DO
  END DO
!$OMP END DO NOWAIT

!--- Filter
  jjb=jj_begin-1
  jje=jj_end
  IF(pole_nord) jjb=jjb+1
  IF(pole_sud)  jje=jje-1
  CALL filtreg_p(rot,jjb_v,jje_v,jjb,jje,jjm,llm,2,1,.FALSE.,1)

!--- Potential vorticity ; correction: rot(iip1,j,l) = rot(1,j,l)
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    IF(pole_sud) ije=ij_end-iip1-1
    DO ij=ijb,ije
      vorpot(ij,l)=(rot(ij,l)+fext(ij))/massebxy(ij,l)
    END DO
    IF(pole_sud) ije=ij_end-iip1
    DO ij=ijb+iip1-1,ije,iip1; vorpot(ij,l)=vorpot(ij-iim,l); END DO
  END DO
!$OMP END DO NOWAIT

END SUBROUTINE tourpot_loc

