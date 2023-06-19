SUBROUTINE convmas2_loc (convm)
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr. Hourdin.
!-------------------------------------------------------------------------------
! Purpose: Compute mass flux convergence at p levels.
!          Equivalent to convmas_loc if convmas1_loc is called before.
  USE parallel_lmdz
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(INOUT) :: convm(ijb_u:ije_u,llm)
!===============================================================================
! Method used:   Computation from top to bottom.
!   Mass convergence at level llm is equal to zero and is not stored in convm.
!===============================================================================
! Local variables:
  INTEGER :: l, ijb, ije
!===============================================================================

!$OMP MASTER
!--- Mass convergence is integrated from top to bottom
  ijb=ij_begin
  ije=ij_end+iip1
  IF(pole_sud) ije=ij_end
  DO l=llmm1,1,-1
    convm(ijb:ije,l) = convm(ijb:ije,l) + convm(ijb:ije,l+1)
  END DO
!$OMP END MASTER

END SUBROUTINE convmas2_loc

