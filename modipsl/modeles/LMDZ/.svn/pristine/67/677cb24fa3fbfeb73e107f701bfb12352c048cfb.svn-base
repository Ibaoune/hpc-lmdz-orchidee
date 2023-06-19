MODULE regr_conserv_m

  USE assert_eq_m,   ONLY: assert_eq
  USE assert_m,      ONLY: assert
  USE interpolation, ONLY: locate

  IMPLICIT NONE

! Purpose: Each procedure regrids a piecewise linear function (not necessarily
!          continuous) by averaging it. This is a conservative regridding.
!          The regridding operation is done along dimension "ix" of the input
!          field "vs". Input are positions of cell edges.
! Remarks:
!   * The target grid should be included in the source grid:
!     no extrapolation is allowed.
!   * Field on target grid "vt" has same rank, slopes and averages as "vs".
!   * If optional argument "slope" is not given, 0 is assumed for slopes.
!     Then the regridding is first order instead of second.
!   * "vs" and "vt" have the same dimensions except Nr. ix (ns for vs, nt for vt)

! Argument                         Type         Description
!-------------------------------------------------------------------------------
!  INTEGER, INTENT(IN)  :: ix      Scalar       dimension regridded <=rank(vs)
!  REAL,    INTENT(IN)  :: vs(*)   Rank>=1      averages in source grid cells
!  REAL,    INTENT(IN)  :: xs(:)   Vector(ns+1) edges of source grid, asc. order
!  REAL,    INTENT(IN)  :: xt(:)   Vector(nt+1) edges of target grid, asc. order
!  REAL,    INTENT(OUT) :: vt(*)   Rank>=1      regridded field
! [REAL,    INTENT(IN)  :: slope]  Rank>=1      slopes in source grid cells

  INTERFACE regr_conserv
    ! The procedures differ only from the rank of the input/output fields.
    MODULE PROCEDURE regr1_conserv, regr2_conserv, regr3_conserv, &
                     regr4_conserv, regr5_conserv
  END INTERFACE

  PRIVATE
  PUBLIC :: regr_conserv

CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE regr1_conserv(ix, vs, xs, xt, vt, slope)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,        INTENT(IN)  :: ix
  REAL,           INTENT(IN)  :: vs(:)
  REAL,           INTENT(IN)  :: xs(:), xt(:)
  REAL,           INTENT(OUT) :: vt(:)
  REAL, OPTIONAL, INTENT(IN)  :: slope(:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: is, it, ns, nt
  REAL    :: co, idt
  LOGICAL :: lslope
!-------------------------------------------------------------------------------
  lslope=PRESENT(slope)
  CALL check_size(ix,SHAPE(vs),SHAPE(vt),xs,xt,ns,nt)
  is=locate(xs,xt(1))                    !--- 1<= is <= ns (no extrapolation)
  vt(:)=0.
  DO it=1, nt; idt=1./(xt(it+1)-xt(it))
    IF(xt(it+1)<=xs(is+1)) THEN          !--- xs(is)<=xt(it)<=xt(it+1)<=xs(is+1)
      CALL mean_lin(xt(it),xt(it+1))
    ELSE
      CALL mean_lin(xt(it),xs(is+1)); is=is+1
      DO WHILE(xs(is+1)<xt(it+1))        !--- 1<=is<=ns-1
      CALL mean_lin(xs(is),xs(is+1)); is=is+1
      END DO
      CALL mean_lin(xs(is),xt(it+1))     !--- 1<=is<=ns
    END IF
    IF(xs(is+1)==xt(it+1)) is=is+1       !--- 1<=is<=ns .OR. it==nt
  END DO

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE mean_lin(a,b)  ! mean [a,b] of the linear function in [xs(is),xs(is+1)]
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN)    :: a, b
!-------------------------------------------------------------------------------
  IF(lslope.AND.(a/=xs(is).OR.b/=xs(is+1))) THEN; co=0.
    IF(a==xt(it  )) co=co+xt(it  )-xs(is  )
    IF(b==xt(it+1)) co=co+xt(it+1)-xs(is+1)
    vt(it) = vt(it)+idt*(b-a)*(vs(is)+co*slope(is)/2.)
  ELSE
    vt(it) = vt(it)+idt*(b-a)* vs(is)
  END IF

END SUBROUTINE mean_lin
!-------------------------------------------------------------------------------

END SUBROUTINE regr1_conserv
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE regr2_conserv(ix, vs, xs, xt, vt, slope)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,        INTENT(IN)  :: ix
  REAL,           INTENT(IN)  :: vs(:,:)
  REAL,           INTENT(IN)  :: xs(:), xt(:)
  REAL,           INTENT(OUT) :: vt(:,:)
  REAL, OPTIONAL, INTENT(IN)  :: slope(:,:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: is, it, ns, nt
  REAL    :: co, idt
  LOGICAL :: lslope
!-------------------------------------------------------------------------------
  lslope=PRESENT(slope)
  CALL check_size(ix,SHAPE(vs),SHAPE(vt),xs,xt,ns,nt)
  is=locate(xs,xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
  vt(:,:)=0.
  DO it=1, nt; idt=1./(xt(it+1)-xt(it))
    IF(xt(it+1)<=xs(is+1)) THEN          !--- xs(is)<=xt(it)<=xt(it+1)<=xs(is+1)
      CALL mean_lin(xt(it),xt(it+1))
    ELSE
      CALL mean_lin(xt(it),xs(is+1)); is=is+1
      DO WHILE(xs(is+1)<xt(it+1))        !--- 1<=is<=ns-1
      CALL mean_lin(xs(is),xs(is+1)); is=is+1
      END DO
      CALL mean_lin(xs(is),xt(it+1))     !--- 1<=is<=ns
    END IF
    IF(xs(is+1)==xt(it+1)) is=is+1       !--- 1<=is<=ns .OR. it==nt
  END DO

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE mean_lin(a,b)  ! mean [a,b] of the linear function in [xs(is),xs(is+1)]
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN)    :: a, b
!-------------------------------------------------------------------------------
  IF(lslope.AND.(a/=xs(is).OR.b/=xs(is+1))) THEN; co=0.
    IF(a==xt(it  )) co=co+xt(it  )-xs(is  )
    IF(b==xt(it+1)) co=co+xt(it+1)-xs(is+1)
    IF(ix==1) vt(it,:) = vt(it,:)+idt*(b-a)*(vs(is,:)+co*slope(is,:)/2.)
    IF(ix==2) vt(:,it) = vt(:,it)+idt*(b-a)*(vs(:,is)+co*slope(:,is)/2.)
  ELSE
    IF(ix==1) vt(it,:) = vt(it,:)+idt*(b-a)* vs(is,:)
    IF(ix==2) vt(:,it) = vt(:,it)+idt*(b-a)* vs(:,is)
  END IF

END SUBROUTINE mean_lin
!-------------------------------------------------------------------------------

END SUBROUTINE regr2_conserv
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE regr3_conserv(ix, vs, xs, xt, vt, slope)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,        INTENT(IN)  :: ix
  REAL,           INTENT(IN)  :: vs(:,:,:)
  REAL,           INTENT(IN)  :: xs(:), xt(:)
  REAL,           INTENT(OUT) :: vt(:,:,:)
  REAL, OPTIONAL, INTENT(IN)  :: slope(:,:,:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: is, it, ns, nt
  REAL    :: co, idt
  LOGICAL :: lslope
!-------------------------------------------------------------------------------
  lslope=PRESENT(slope)
  CALL check_size(ix,SHAPE(vs),SHAPE(vt),xs,xt,ns,nt)
  is=locate(xs,xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
  vt(:,:,:)=0.
  DO it=1, nt; idt=1./(xt(it+1)-xt(it))
    IF(xt(it+1)<=xs(is+1)) THEN          !--- xs(is)<=xt(it)<=xt(it+1)<=xs(is+1)
      CALL mean_lin(xt(it),xt(it+1))
    ELSE
      CALL mean_lin(xt(it),xs(is+1)); is=is+1
      DO WHILE(xs(is+1)<xt(it+1))        !--- 1<=is<=ns-1
      CALL mean_lin(xs(is),xs(is+1)); is=is+1
      END DO
      CALL mean_lin(xs(is),xt(it+1))     !--- 1<=is<=ns
    END IF
    IF(xs(is+1)==xt(it+1)) is=is+1       !--- 1<=is<=ns .OR. it==nt
  END DO

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE mean_lin(a,b)  ! mean [a,b] of the linear function in [xs(is),xs(is+1)]
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN)    :: a, b
!-------------------------------------------------------------------------------
  IF(lslope.AND.(a/=xs(is).OR.b/=xs(is+1))) THEN; co=0.
    IF(a==xt(it  )) co=co+xt(it  )-xs(is  )
    IF(b==xt(it+1)) co=co+xt(it+1)-xs(is+1)
    IF(ix==1) vt(it,:,:) = vt(it,:,:)+idt*(b-a)*(vs(is,:,:)+co*slope(is,:,:)/2.)
    IF(ix==2) vt(:,it,:) = vt(:,it,:)+idt*(b-a)*(vs(:,is,:)+co*slope(:,is,:)/2.)
    IF(ix==3) vt(:,:,it) = vt(:,:,it)+idt*(b-a)*(vs(:,:,is)+co*slope(:,:,is)/2.)
  ELSE
    IF(ix==1) vt(it,:,:) = vt(it,:,:)+idt*(b-a)* vs(is,:,:)
    IF(ix==2) vt(:,it,:) = vt(:,it,:)+idt*(b-a)* vs(:,is,:)
    IF(ix==3) vt(:,:,it) = vt(:,:,it)+idt*(b-a)* vs(:,:,is)
  END IF

END SUBROUTINE mean_lin
!-------------------------------------------------------------------------------

END SUBROUTINE regr3_conserv
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE regr4_conserv(ix, vs, xs, xt, vt, slope)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,        INTENT(IN)  :: ix
  REAL,           INTENT(IN)  :: vs(:,:,:,:)
  REAL,           INTENT(IN)  :: xs(:), xt(:)
  REAL,           INTENT(OUT) :: vt(:,:,:,:)
  REAL, OPTIONAL, INTENT(IN)  :: slope(:,:,:,:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: is, it, ns, nt
  REAL    :: co, idt
  LOGICAL :: lslope
!-------------------------------------------------------------------------------
  lslope=PRESENT(slope)
  CALL check_size(ix,SHAPE(vs),SHAPE(vt),xs,xt,ns,nt)
  is=locate(xs,xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
  vt(:,:,:,:)=0.
  DO it=1, nt; idt=1./(xt(it+1)-xt(it))
    IF(xt(it+1)<=xs(is+1)) THEN          !--- xs(is)<=xt(it)<=xt(it+1)<=xs(is+1)
      CALL mean_lin(xt(it),xt(it+1))
    ELSE
      CALL mean_lin(xt(it),xs(is+1)); is=is+1
      DO WHILE(xs(is+1)<xt(it+1))        !--- 1<=is<=ns-1
      CALL mean_lin(xs(is),xs(is+1)); is=is+1
      END DO
      CALL mean_lin(xs(is),xt(it+1))     !--- 1<=is<=ns
    END IF
    IF(xs(is+1)==xt(it+1)) is=is+1       !--- 1<=is<=ns .OR. it==nt
  END DO

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE mean_lin(a,b)  ! mean [a,b] of the linear function in [xs(is),xs(is+1)]
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN) :: a, b
!-------------------------------------------------------------------------------
  IF(lslope.AND.(a/=xs(is).OR.b/=xs(is+1))) THEN; co=0.
    IF(a==xt(it  )) co=co+xt(it  )-xs(is  )
    IF(b==xt(it+1)) co=co+xt(it+1)-xs(is+1)
    IF(ix==1) vt(it,:,:,:) = vt(it,:,:,:)+idt*(b-a)*(vs(is,:,:,:)+co*slope(is,:,:,:)/2.)
    IF(ix==2) vt(:,it,:,:) = vt(:,it,:,:)+idt*(b-a)*(vs(:,is,:,:)+co*slope(:,is,:,:)/2.)
    IF(ix==3) vt(:,:,it,:) = vt(:,:,it,:)+idt*(b-a)*(vs(:,:,is,:)+co*slope(:,:,is,:)/2.)
    IF(ix==4) vt(:,:,:,it) = vt(:,:,:,it)+idt*(b-a)*(vs(:,:,:,is)+co*slope(:,:,:,is)/2.)
  ELSE
    IF(ix==1) vt(it,:,:,:) = vt(it,:,:,:)+idt*(b-a)* vs(is,:,:,:)
    IF(ix==2) vt(:,it,:,:) = vt(:,it,:,:)+idt*(b-a)* vs(:,is,:,:)
    IF(ix==3) vt(:,:,it,:) = vt(:,:,it,:)+idt*(b-a)* vs(:,:,is,:)
    IF(ix==4) vt(:,:,:,it) = vt(:,:,:,it)+idt*(b-a)* vs(:,:,:,is)
  END IF

END SUBROUTINE mean_lin
!-------------------------------------------------------------------------------

END SUBROUTINE regr4_conserv
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE regr5_conserv(ix, vs, xs, xt, vt, slope)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,        INTENT(IN)  :: ix
  REAL,           INTENT(IN)  :: vs(:,:,:,:,:)
  REAL,           INTENT(IN)  :: xs(:), xt(:)
  REAL,           INTENT(OUT) :: vt(:,:,:,:,:)
  REAL, OPTIONAL, INTENT(IN)  :: slope(:,:,:,:,:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: is, it, ns, nt
  REAL    :: co, idt
  LOGICAL :: lslope
!-------------------------------------------------------------------------------
  lslope=PRESENT(slope)
  CALL check_size(ix,SHAPE(vs),SHAPE(vt),xs,xt,ns,nt)
  is=locate(xs,xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
  vt(:,:,:,:,:)=0.
  DO it=1, nt; idt=1./(xt(it+1)-xt(it))
    IF(xt(it+1)<=xs(is+1)) THEN          !--- xs(is)<=xt(it)<=xt(it+1)<=xs(is+1)
      CALL mean_lin(xt(it),xt(it+1))
    ELSE
      CALL mean_lin(xt(it),xs(is+1)); is=is+1
      DO WHILE(xs(is+1)<xt(it+1))        !--- 1<=is<=ns-1
      CALL mean_lin(xs(is),xs(is+1)); is=is+1
      END DO
      CALL mean_lin(xs(is),xt(it+1))     !--- 1<=is<=ns
    END IF
    IF(xs(is+1)==xt(it+1)) is=is+1     !--- 1<=is<=ns .OR. it==nt
  END DO

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE mean_lin(a,b)  ! mean [a,b] of the linear function in [xs(is),xs(is+1)]
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN) :: a, b
!-------------------------------------------------------------------------------
  IF(lslope.AND.(a/=xs(is).OR.b/=xs(is+1))) THEN; co=0.
    IF(a==xt(it  )) co=co+xt(it  )-xs(is  )
    IF(b==xt(it+1)) co=co+xt(it+1)-xs(is+1)
    IF(ix==1) vt(it,:,:,:,:) = vt(it,:,:,:,:)+idt*(b-a)*(vs(is,:,:,:,:)+co*slope(is,:,:,:,:)/2.)
    IF(ix==2) vt(:,it,:,:,:) = vt(:,it,:,:,:)+idt*(b-a)*(vs(:,is,:,:,:)+co*slope(:,is,:,:,:)/2.)
    IF(ix==3) vt(:,:,it,:,:) = vt(:,:,it,:,:)+idt*(b-a)*(vs(:,:,is,:,:)+co*slope(:,:,is,:,:)/2.)
    IF(ix==4) vt(:,:,:,it,:) = vt(:,:,:,it,:)+idt*(b-a)*(vs(:,:,:,is,:)+co*slope(:,:,:,is,:)/2.)
    IF(ix==5) vt(:,:,:,:,it) = vt(:,:,:,:,it)+idt*(b-a)*(vs(:,:,:,:,is)+co*slope(:,:,:,:,is)/2.)
  ELSE
    IF(ix==1) vt(it,:,:,:,:) = vt(it,:,:,:,:)+idt*(b-a)* vs(is,:,:,:,:)
    IF(ix==2) vt(:,it,:,:,:) = vt(:,it,:,:,:)+idt*(b-a)* vs(:,is,:,:,:)
    IF(ix==3) vt(:,:,it,:,:) = vt(:,:,it,:,:)+idt*(b-a)* vs(:,:,is,:,:)
    IF(ix==4) vt(:,:,:,it,:) = vt(:,:,:,it,:)+idt*(b-a)* vs(:,:,:,is,:)
    IF(ix==5) vt(:,:,:,:,it) = vt(:,:,:,:,it)+idt*(b-a)* vs(:,:,:,:,is)
  END IF

END SUBROUTINE mean_lin
!-------------------------------------------------------------------------------

END SUBROUTINE regr5_conserv
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE check_size(ix,svs,svt,xs,xt,ns,nt)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN)  :: ix, svs(:), svt(:)
  REAL,    INTENT(IN)  :: xs(:), xt(:)
  INTEGER, INTENT(OUT) :: ns,    nt
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: rk, is, ii, n
  CHARACTER(LEN=80) :: sub, msg
!-------------------------------------------------------------------------------
  WRITE(sub,'(a,2i0,a)')"regr",SIZE(svs),ix,"_conserv"
  rk=assert_eq(SIZE(svs),SIZE(svt),TRIM(sub)//': inconsistent ranks')
  WRITE(msg,'(a,2(i0,a))')TRIM(sub)//": ix (",ix,") exceeds field rank (",rk,")"
  CALL assert(ix>=1.AND.ix<=rk,msg)
  ii=0
  DO is=1,rk; IF(is==ix) CYCLE
    WRITE(msg,'(a,i0)')TRIM(sub)//" n",is
    n=assert_eq(svs(is),svt(is),msg)
    ii=ii+1
  END DO
  ns=assert_eq(svs(ix),SIZE(xs)-1,TRIM(sub)//" ns")
  nt=assert_eq(svt(ix),SIZE(xt)-1,TRIM(sub)//" nt")

  !--- Quick check on sort order:
  CALL assert(xs(1)<xs(2), TRIM(sub)//": xs bad order")
  CALL assert(xt(1)<xt(2), TRIM(sub)//": xt bad order")
  CALL assert(xs(1)<=xt(1).AND.xt(nt+1)<=xs(ns+1), TRIM(sub)//": extrapolation")

END SUBROUTINE check_size
!
!-------------------------------------------------------------------------------


END MODULE regr_conserv_m
