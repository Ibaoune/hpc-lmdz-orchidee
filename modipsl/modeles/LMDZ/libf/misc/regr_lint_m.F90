MODULE regr_lint_m

  USE assert_eq_m,   ONLY: assert_eq
  USE assert_m,      ONLY: assert
  USE interpolation, ONLY: hunt

  IMPLICIT NONE

! Purpose: Each procedure regrids by linear interpolation along dimension "ix"
!          the input field "vs" to the output field "vt".
! Remark:
!   * "vs" and "vt" have the same dimensions except Nr. ix (ns for vs, nt for vt)

! Argument                         Type         Description
!-------------------------------------------------------------------------------
!  INTEGER, INTENT(IN)  :: ix      Scalar       dimension regridded <=rank(vs)
!  REAL,    INTENT(IN)  :: vs(*)   Rank>=1      source grid field values
!  REAL,    INTENT(IN)  :: xs(:)   Vector(ns)   centers of source grid, asc. order
!  REAL,    INTENT(IN)  :: xt(:)   Vector(nt)   centers of target grid, asc. order
!  REAL,    INTENT(OUT) :: vt(*)   Rank>=1      regridded field

  INTERFACE regr_lint
    ! The procedures differ only from the rank of the input/output fields.
    MODULE PROCEDURE regr1_lint, regr2_lint, regr3_lint, &
                     regr4_lint, regr5_lint
  END INTERFACE

  PRIVATE
  PUBLIC :: regr_lint

CONTAINS


!-------------------------------------------------------------------------------
!
SUBROUTINE regr1_lint(ix, vs, xs, xt, vt)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN)  :: ix
  REAL,    INTENT(IN)  :: vs(:)
  REAL,    INTENT(IN)  :: xs(:)
  REAL,    INTENT(IN)  :: xt(:)
  REAL,    INTENT(OUT) :: vt(:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: is, it, ns, nt, isb ! "is" bound between 1 and "ns - 1"
  REAL    :: r
!-------------------------------------------------------------------------------
  CALL check_size(ix,SHAPE(vs),SHAPE(vt),SIZE(xs),SIZE(xt),ns,nt)
  is = -1 ! go immediately to bisection on first call to "hunt"
  DO it=1,SIZE(xt)
    CALL hunt(xs,xt(it),is); isb=MIN(MAX(is,1),ns-1)
    r=(xs(isb+1)-xt(it))/(xs(isb+1)-xs(isb))
    vt(it)=r*vs(isb)+(1.-r)*vs(isb+1)
  END DO

END SUBROUTINE regr1_lint
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE regr2_lint(ix, vs, xs, xt, vt)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN)  :: ix
  REAL,    INTENT(IN)  :: vs(:,:)
  REAL,    INTENT(IN)  :: xs(:)
  REAL,    INTENT(IN)  :: xt(:)
  REAL,    INTENT(OUT) :: vt(:,:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: is, it, ns, nt, isb ! "is" bound between 1 and "ns - 1"
  REAL    :: r
!-------------------------------------------------------------------------------
  CALL check_size(ix,SHAPE(vs),SHAPE(vt),SIZE(xs),SIZE(xt),ns,nt)
  is = -1 ! go immediately to bisection on first call to "hunt"
  DO it=1,SIZE(xt)
    CALL hunt(xs,xt(it),is); isb=MIN(MAX(is,1),ns-1)
    r=(xs(isb+1)-xt(it))/(xs(isb+1)-xs(isb))
    IF(ix==1) vt(it,:)=r*vs(isb,:)+(1.-r)*vs(isb+1,:)
    IF(ix==2) vt(:,it)=r*vs(:,isb)+(1.-r)*vs(:,isb+1)
  END DO

END SUBROUTINE regr2_lint
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE regr3_lint(ix, vs, xs, xt, vt)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN)  :: ix
  REAL,    INTENT(IN)  :: vs(:,:,:)
  REAL,    INTENT(IN)  :: xs(:)
  REAL,    INTENT(IN)  :: xt(:)
  REAL,    INTENT(OUT) :: vt(:,:,:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: is, it, ns, nt, isb ! "is" bound between 1 and "ns - 1"
  REAL    :: r
!-------------------------------------------------------------------------------
  CALL check_size(ix,SHAPE(vs),SHAPE(vt),SIZE(xs),SIZE(xt),ns,nt)
  is = -1 ! go immediately to bisection on first call to "hunt"
  DO it=1,SIZE(xt)
    CALL hunt(xs,xt(it),is); isb=MIN(MAX(is,1),ns-1)
    r=(xs(isb+1)-xt(it))/(xs(isb+1)-xs(isb))
    IF(ix==1) vt(it,:,:)=r*vs(isb,:,:)+(1.-r)*vs(isb+1,:,:)
    IF(ix==2) vt(:,it,:)=r*vs(:,isb,:)+(1.-r)*vs(:,isb+1,:)
    IF(ix==3) vt(:,:,it)=r*vs(:,:,isb)+(1.-r)*vs(:,:,isb+1)
  END DO

END SUBROUTINE regr3_lint
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE regr4_lint(ix, vs, xs, xt, vt)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN)  :: ix
  REAL,    INTENT(IN)  :: vs(:,:,:,:)
  REAL,    INTENT(IN)  :: xs(:)
  REAL,    INTENT(IN)  :: xt(:)
  REAL,    INTENT(OUT) :: vt(:,:,:,:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: is, it, ns, nt, isb ! "is" bound between 1 and "ns - 1"
  REAL    :: r
!-------------------------------------------------------------------------------
  CALL check_size(ix,SHAPE(vs),SHAPE(vt),SIZE(xs),SIZE(xt),ns,nt)
  is = -1 ! go immediately to bisection on first call to "hunt"
  DO it=1,SIZE(xt)
    CALL hunt(xs,xt(it),is); isb=MIN(MAX(is,1),ns-1)
    r=(xs(isb+1)-xt(it))/(xs(isb+1)-xs(isb))
    IF(ix==1) vt(it,:,:,:)=r*vs(isb,:,:,:)+(1.-r)*vs(isb+1,:,:,:)
    IF(ix==2) vt(:,it,:,:)=r*vs(:,isb,:,:)+(1.-r)*vs(:,isb+1,:,:)
    IF(ix==3) vt(:,:,it,:)=r*vs(:,:,isb,:)+(1.-r)*vs(:,:,isb+1,:)
    IF(ix==4) vt(:,:,:,it)=r*vs(:,:,:,isb)+(1.-r)*vs(:,:,:,isb+1)
  END DO

END SUBROUTINE regr4_lint
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE regr5_lint(ix, vs, xs, xt, vt)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN)  :: ix
  REAL,    INTENT(IN)  :: vs(:,:,:,:,:)
  REAL,    INTENT(IN)  :: xs(:)
  REAL,    INTENT(IN)  :: xt(:)
  REAL,    INTENT(OUT) :: vt(:,:,:,:,:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: is, it, ns, nt, isb ! "is" bound between 1 and "ns - 1"
  REAL    :: r
!-------------------------------------------------------------------------------
  CALL check_size(ix,SHAPE(vs),SHAPE(vt),SIZE(xs),SIZE(xt),ns,nt)
  is = -1 ! go immediately to bisection on first call to "hunt"
  DO it=1,SIZE(xt)
    CALL hunt(xs,xt(it),is); isb=MIN(MAX(is,1),ns-1)
    r=(xs(isb+1)-xt(it))/(xs(isb+1)-xs(isb))
    IF(ix==1) vt(it,:,:,:,:)=r*vs(isb,:,:,:,:)+(1.-r)*vs(isb+1,:,:,:,:)
    IF(ix==2) vt(:,it,:,:,:)=r*vs(:,isb,:,:,:)+(1.-r)*vs(:,isb+1,:,:,:)
    IF(ix==3) vt(:,:,it,:,:)=r*vs(:,:,isb,:,:)+(1.-r)*vs(:,:,isb+1,:,:)
    IF(ix==4) vt(:,:,:,it,:)=r*vs(:,:,:,isb,:)+(1.-r)*vs(:,:,:,isb+1,:)
    IF(ix==5) vt(:,:,:,:,it)=r*vs(:,:,:,:,isb)+(1.-r)*vs(:,:,:,:,isb+1)
  END DO

END SUBROUTINE regr5_lint
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE check_size(ix,svs,svt,nxs,nxt,ns,nt)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN)  :: ix, svs(:), svt(:), nxs, nxt
  INTEGER, INTENT(OUT) :: ns, nt
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: rk, is
  CHARACTER(LEN=80) :: sub, msg
!-------------------------------------------------------------------------------
  rk=SIZE(svs)
  WRITE(sub,'(a,2i0,a)')"regr",rk,ix,"_lint"
  CALL assert(ix>=1.AND.ix<=rk,TRIM(sub)//": ix exceeds fields rank")
  DO is=1,rk; IF(is==ix) CYCLE
    WRITE(msg,'(a,i1)')TRIM(sub)//" n",is
    CALL assert(svs(is)==svt(is),msg)
  END DO
  ns=assert_eq(svs(ix),nxs,TRIM(sub)//" ns")
  nt=assert_eq(svt(ix),nxt,TRIM(sub)//" nt")

END SUBROUTINE check_size
!
!-------------------------------------------------------------------------------

END MODULE regr_lint_m
!
!-------------------------------------------------------------------------------

