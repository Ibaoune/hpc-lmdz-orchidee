!
! $Id: minmaxsimple.F90 3526 2019-05-28 13:00:44Z fairhead $
!
SUBROUTINE minmaxsimple(zq,qmin,qmax,comment)
  USE dimphy
  IMPLICIT NONE

! Entrees
  REAL,DIMENSION(klon,klev), INTENT(IN)   :: zq
  REAL,INTENT(IN)                         :: qmin,qmax
  CHARACTER(LEN=*),INTENT(IN)             :: comment

! Local  
  REAL zmin, zmax
  
  zmin=MINVAL(zq)
  zmax=MAXVAL(zq)
  PRINT *, "qmin qmax=", zmin, zmax, comment
!  IF (zmin.LT.qmin.OR.zmax.GT.qmax) THEN
!      WRITE(*,*) "qmin qmax=", zmin, zmax, comment
!  ENDIF
  
END SUBROUTINE minmaxsimple
