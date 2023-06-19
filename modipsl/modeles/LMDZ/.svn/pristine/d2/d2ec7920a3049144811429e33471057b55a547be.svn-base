SUBROUTINE caldyn0(itau,ucov,vcov,teta,ps,masse,pk,phis,phi,w,pbaru,pbarv,time)
!
!-------------------------------------------------------------------------------
! Author: P. Le Van ; modif. 04/93: F.Forget.
!-------------------------------------------------------------------------------
! Purpose: Compute dynamic tendencies.
!-------------------------------------------------------------------------------
  USE control_mod, ONLY: resetvarc 
  USE comvert_mod, ONLY: ap, bp
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  INTEGER, INTENT(IN)  :: itau                      !--- TIME STEP INDEX
  REAL,    INTENT(IN)  :: vcov (ip1jm    ,llm)      !--- V COVARIANT WIND
  REAL,    INTENT(IN)  :: ucov (ip1jmp1  ,llm)      !--- U COVARIANT WIND
  REAL,    INTENT(IN)  :: teta (ip1jmp1  ,llm)      !--- POTENTIAL TEMPERATURE
  REAL,    INTENT(IN)  :: ps   (ip1jmp1)            !--- GROUND PRESSURE
  REAL,    INTENT(OUT) :: masse(ip1jmp1  ,llm)      !--- MASS IN EACH CELL
  REAL,    INTENT(IN)  :: pk   (iip1,jjp1,llm)      !--- PRESSURE
  REAL,    INTENT(IN)  :: phis (ip1jmp1)            !--- GROUND GEOPOTENTIAL
  REAL,    INTENT(IN)  :: phi  (ip1jmp1  ,llm)      !--- 3D GEOPOTENTIAL
  REAL,    INTENT(OUT) :: w    (ip1jmp1  ,llm)      !--- VERTICAL WIND
  REAL,    INTENT(OUT) :: pbaru(ip1jmp1  ,llm)      !--- U MASS FLUX
  REAL,    INTENT(OUT) :: pbarv(ip1jm    ,llm)      !--- V MASS FLUX
  REAL,    INTENT(IN)  :: time                      !--- TIME
!===============================================================================
! Local variables:
  REAL, DIMENSION(ip1jmp1,llmp1) :: p
  REAL, DIMENSION(ip1jmp1,llm)   :: ucont, massebx, ang, ecin, convm, bern
  REAL, DIMENSION(ip1jmp1)       :: dp
  REAL, DIMENSION(ip1jm  ,llm)   :: vcont, masseby, massebxy, vorpot
  REAL, DIMENSION(ip1jm)         :: psexbarxy
  INTEGER                        :: ij, l
!===============================================================================
  CALL covcont  ( llm    , ucov    , vcov , ucont, vcont        )
  CALL pression ( ip1jmp1, ap      , bp   ,  ps  , p            )
  CALL psextbar (   ps   , psexbarxy                            )
  CALL massdair (    p   , masse                                )
  CALL massbar  (   masse, massebx , masseby                    )
  CALL massbarxy(   masse, massebxy                             )
  CALL flumass  ( massebx, masseby , vcont, ucont ,pbaru, pbarv )
  CALL convmas  (   pbaru, pbarv   , convm                      )
  CALL vitvert  ( convm  , w                                    )
  CALL tourpot  ( vcov   , ucov    , massebxy  , vorpot         )
  CALL enercin  ( vcov   , ucov    , vcont     , ucont  , ecin  )
  CALL bernoui  ( ip1jmp1, llm     , phi       , ecin   , bern  )
  DO l=1,llm; ang(:,l) = ucov(:,l) + constang(:); END DO
  resetvarc=.true. ! force a recomputation of initial values in sortvarc
  dp(:)=convm(:,1)/airesurg(:)
  CALL sortvarc( itau,ucov,teta,ps,masse,pk,phis,vorpot,phi,bern,dp,time,vcov )

END SUBROUTINE caldyn0
