SUBROUTINE convmas1_loc (pbaru, pbarv, convm)
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr. Hourdin.
!-------------------------------------------------------------------------------
! Purpose: Compute mass flux convergence at p levels.
!          Equivalent to convmas_loc if convmas2_loc is called after.
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
  REAL, TARGET, INTENT(OUT) :: convm(ijb_u:ije_u,llm)
!===============================================================================
! Method used:   Computation from top to bottom.
!   Mass convergence at level llm is equal to zero and is not stored in convm.
!===============================================================================
! Local variables:
  INTEGER :: l, jjb, jje
!===============================================================================

!--- Computation of - (d(pbaru)/dx + d(pbarv)/dy )
  CALL convflu_loc( pbaru, pbarv, llm, convm )

!--- Filter
  jjb=jj_begin
  jje=jj_end+1
  IF(pole_sud) jje=jj_end
  CALL filtreg_p(convm,jjb_u,jje_u,jjb,jje,jjp1,llm,2,2,.TRUE.,1)

END SUBROUTINE convmas1_loc

