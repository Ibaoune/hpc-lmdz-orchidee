SUBROUTINE convmas_loc (pbaru, pbarv, convm)
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr. Hourdin.
!-------------------------------------------------------------------------------
! Purpose: Compute mass flux convergence at p levels.
  USE parallel_lmdz
  USE mod_filtreg_p
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: pbaru(ijb_u:ije_u,llm)
  REAL, INTENT(IN)  :: pbarv(ijb_v:ije_v,llm)
  REAL, INTENT(OUT) :: convm(ijb_u:ije_u,llm)
!===============================================================================
! Method used:   Computation from top to bottom.
!   Mass convergence at level llm is equal to zero and is not stored in convm.
!===============================================================================
! Local variables:
  INTEGER :: l, ijb, ije, jjb, jje
!===============================================================================

!--- Computation of - (d(pbaru)/dx + d(pbarv)/dy )
  CALL convflu_loc( pbaru, pbarv, llm, convm )

!--- Filter
  jjb=jj_begin
  jje=jj_end+1
  IF(pole_sud) jje=jj_end
  CALL filtreg_p(convm,jjb_u,jje_u,jjb,jje,jjp1,llm,2,2,.TRUE.,1)

!--- Mass convergence is integrated from top to bottom
!$OMP BARRIER
!$OMP MASTER
  ijb=ij_begin
  ije=ij_end+iip1
  IF(pole_sud) ije=ij_end
  DO l=llmm1,1,-1
    convm(ijb:ije,l) = convm(ijb:ije,l) + convm(ijb:ije,l+1)
  END DO
!$OMP END MASTER
!$OMP BARRIER

END SUBROUTINE convmas_loc

