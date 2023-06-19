SUBROUTINE vitvert_loc(convm, w)
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr. Hourdin.
!-------------------------------------------------------------------------------
! Purpose: Compute vertical speed at sigma levels.
  USE parallel_lmdz
  USE comvert_mod, ONLY: bp
  
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: convm(ijb_u:ije_u,llm)
  REAL, INTENT(OUT) :: w    (ijb_u:ije_u,llm)
!===============================================================================
! Notes: Vertical speed is oriented from bottom to top.
!   * At ground - level sigma(1):     w(i,j,1) = 0.
!   * At top    - level sigma(llm+1): w(i,j,l) = 0. (not stored in w)
!===============================================================================
! Local variables:
  INTEGER :: l, ijb, ije
!===============================================================================
  ijb=ij_begin
  ije=ij_end+iip1
  IF(pole_sud) ije=ij_end
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
  DO l=1,llmm1
    w(ijb:ije,l+1)=convm(ijb:ije,l+1)-bp(l+1)*convm(ijb:ije,1)
  END DO
!$OMP END DO
!$OMP MASTER
  w(ijb:ije,1)=0.
!$OMP END MASTER
!$OMP BARRIER

END SUBROUTINE vitvert_loc

