MODULE slopes_m

  ! Author: Lionel GUEZ
  ! Extension / factorisation: David CUGNET

  IMPLICIT NONE

  ! Those generic function computes second order slopes with Van
  ! Leer slope-limiting, given cell averages. Reference: Dukowicz,
  ! 1987, SIAM Journal on Scientific and Statistical Computing, 8,
  ! 305.

  ! The only difference between the specific functions is the rank
  ! of the first argument and the equal rank of the result.

  ! slope(ix,...) acts on ix th dimension.

  ! real, intent(in), rank >= 1:: f ! (n, ...) cell averages, n must be >= 1
  ! real, intent(in):: x(:) ! (n + 1) cell edges
  ! real slopes, same shape as f ! (n, ...)
  INTERFACE slopes
     MODULE procedure slopes1, slopes2, slopes3, slopes4, slopes5
  END INTERFACE

  PRIVATE
  PUBLIC :: slopes

CONTAINS

!-------------------------------------------------------------------------------
!
PURE FUNCTION slopes1(ix, f, x)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN) :: ix
  REAL,    INTENT(IN) :: f(:)
  REAL,    INTENT(IN) :: x(:)
  REAL :: slopes1(SIZE(f,1))
!-------------------------------------------------------------------------------
! Local:
  INTEGER :: n, i, j, sta(2), sto(2)
  REAL :: xc(SIZE(f,1))                             ! (n) cell centers
  REAL :: h(2:SIZE(f,1)-1), delta_xc(2:SIZE(f,1)-1) ! (2:n-1)
  REAL :: fm, ff, fp, dx
!-------------------------------------------------------------------------------
  n=SIZE(f,ix)
  FORALL(i=1:n) xc(i)=(x(i)+x(i+1))/2.
  FORALL(i=2:n-1)
    h(i)=ABS(x(i+1)-xc(i)) ; delta_xc(i)=xc(i+1)-xc(i-1)
  END FORALL
  slopes1(:)=0.
  DO i=2,n-1
    ff=f(i); fm=f(i-1); fp=f(i+1)
    IF(ff>=MAX(fm,fp).OR.ff<=MIN(fm,fp)) THEN
      slopes1(i)=0.; CYCLE           !--- Local extremum
      !--- 2nd order slope ; (fm, ff, fp) strictly monotonous
      slopes1(i)=(fp-fm)/delta_xc(i)
      !--- Slope limitation
      slopes1(i) = SIGN(MIN(ABS(slopes1(i)), &
        ABS(fp-ff)/h(i),ABS(ff-fm)/h(i)),slopes1(i) )
     END IF
  END DO

END FUNCTION slopes1
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
PURE FUNCTION slopes2(ix, f, x)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN) :: ix
  REAL,    INTENT(IN) :: f(:, :)
  REAL,    INTENT(IN) :: x(:)
  REAL :: slopes2(SIZE(f,1),SIZE(f,2))
!-------------------------------------------------------------------------------
! Local:
  INTEGER :: n, i, j, sta(2), sto(2)
  REAL, ALLOCATABLE :: xc(:)                        ! (n) cell centers
  REAL, ALLOCATABLE :: h(:), delta_xc(:)            ! (2:n-1)
  REAL :: fm, ff, fp, dx
!-------------------------------------------------------------------------------
  n=SIZE(f,ix); ALLOCATE(xc(n),h(2:n-1),delta_xc(2:n-1))
  FORALL(i=1:n) xc(i)=(x(i)+x(i+1))/2.
  FORALL(i=2:n-1)
    h(i)=ABS(x(i+1)-xc(i)) ; delta_xc(i)=xc(i+1)-xc(i-1)
  END FORALL
  slopes2(:,:)=0.
  sta=[1,1]; sta(ix)=2
  sto=SHAPE(f); sto(ix)=n-1
  DO j=sta(2),sto(2); IF(ix==2) dx=delta_xc(j)
    DO i=sta(1),sto(1); IF(ix==1) dx=delta_xc(i)
      ff=f(i,j)
      SELECT CASE(ix)
        CASE(1); fm=f(i-1,j); fp=f(i+1,j)
        CASE(2); fm=f(i,j-1); fp=f(i,j+1)
      END SELECT
      IF(ff>=MAX(fm,fp).OR.ff<=MIN(fm,fp)) THEN
        slopes2(i,j)=0.; CYCLE           !--- Local extremum
        !--- 2nd order slope ; (fm, ff, fp) strictly monotonous
        slopes2(i,j)=(fp-fm)/dx
        !--- Slope limitation
        slopes2(i,j) = SIGN(MIN(ABS(slopes2(i,j)), &
          ABS(fp-ff)/h(i),ABS(ff-fm)/h(i)),slopes2(i,j) )
       END IF
    END DO
  END DO
  DEALLOCATE(xc,h,delta_xc)

END FUNCTION slopes2
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
PURE FUNCTION slopes3(ix, f, x)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN) :: ix
  REAL,    INTENT(IN) :: f(:, :, :)
  REAL,    INTENT(IN) :: x(:)
  REAL :: slopes3(SIZE(f,1),SIZE(f,2),SIZE(f,3))
!-------------------------------------------------------------------------------
! Local:
  INTEGER :: n, i, j, k, sta(3), sto(3)
  REAL, ALLOCATABLE :: xc(:)                        ! (n) cell centers
  REAL, ALLOCATABLE :: h(:), delta_xc(:)            ! (2:n-1)
  REAL :: fm, ff, fp, dx
!-------------------------------------------------------------------------------
  n=SIZE(f,ix); ALLOCATE(xc(n),h(2:n-1),delta_xc(2:n-1))
  FORALL(i=1:n) xc(i)=(x(i)+x(i+1))/2.
  FORALL(i=2:n-1)
    h(i)=ABS(x(i+1)-xc(i)) ; delta_xc(i)=xc(i+1)-xc(i-1)
  END FORALL
  slopes3(:,:,:)=0.
  sta=[1,1,1]; sta(ix)=2
  sto=SHAPE(f); sto(ix)=n-1
  DO k=sta(3),sto(3); IF(ix==3) dx=delta_xc(k)
    DO j=sta(2),sto(2); IF(ix==2) dx=delta_xc(j)
      DO i=sta(1),sto(1); IF(ix==1) dx=delta_xc(i)
        ff=f(i,j,k)
        SELECT CASE(ix)
          CASE(1); fm=f(i-1,j,k); fp=f(i+1,j,k)
          CASE(2); fm=f(i,j-1,k); fp=f(i,j+1,k)
          CASE(3); fm=f(i,j,k-1); fp=f(i,j,k+1)
        END SELECT
        IF(ff>=MAX(fm,fp).OR.ff<=MIN(fm,fp)) THEN
          slopes3(i,j,k)=0.; CYCLE           !--- Local extremum
          !--- 2nd order slope ; (fm, ff, fp) strictly monotonous
          slopes3(i,j,k)=(fp-fm)/dx
          !--- Slope limitation
          slopes3(i,j,k) = SIGN(MIN(ABS(slopes3(i,j,k)), &
            ABS(fp-ff)/h(i),ABS(ff-fm)/h(i)),slopes3(i,j,k) )
         END IF
      END DO
    END DO
  END DO
  DEALLOCATE(xc,h,delta_xc)

END FUNCTION slopes3
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
PURE FUNCTION slopes4(ix, f, x)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN) :: ix
  REAL,    INTENT(IN) :: f(:, :, :, :)
  REAL,    INTENT(IN) :: x(:)
  REAL :: slopes4(SIZE(f,1),SIZE(f,2),SIZE(f,3),SIZE(f,4))
!-------------------------------------------------------------------------------
! Local:
  INTEGER :: n, i, j, k, l, m, sta(4), sto(4)
  REAL, ALLOCATABLE :: xc(:)                        ! (n) cell centers
  REAL, ALLOCATABLE :: h(:), delta_xc(:)            ! (2:n-1)
  REAL :: fm, ff, fp, dx
!-------------------------------------------------------------------------------
  n=SIZE(f,ix); ALLOCATE(xc(n),h(2:n-1),delta_xc(2:n-1))
  FORALL(i=1:n) xc(i)=(x(i)+x(i+1))/2.
  FORALL(i=2:n-1)
    h(i)=ABS(x(i+1)-xc(i)) ; delta_xc(i)=xc(i+1)-xc(i-1)
  END FORALL
  slopes4(:,:,:,:)=0.
  sta=[1,1,1,1]; sta(ix)=2
  sto=SHAPE(f); sto(ix)=n-1
  DO l=sta(4),sto(4); IF(ix==4) dx=delta_xc(l)
    DO k=sta(3),sto(3); IF(ix==3) dx=delta_xc(k)
      DO j=sta(2),sto(2); IF(ix==2) dx=delta_xc(j)
        DO i=sta(1),sto(1); IF(ix==1) dx=delta_xc(i)
          ff=f(i,j,k,l)
          SELECT CASE(ix)
            CASE(1); fm=f(i-1,j,k,l); fp=f(i+1,j,k,l)
            CASE(2); fm=f(i,j-1,k,l); fp=f(i,j+1,k,l)
            CASE(3); fm=f(i,j,k-1,l); fp=f(i,j,k+1,l)
            CASE(4); fm=f(i,j,k,l-1); fp=f(i,j,k,l+1)
          END SELECT
          IF(ff>=MAX(fm,fp).OR.ff<=MIN(fm,fp)) THEN
            slopes4(i,j,k,l)=0.; CYCLE           !--- Local extremum
            !--- 2nd order slope ; (fm, ff, fp) strictly monotonous
            slopes4(i,j,k,l)=(fp-fm)/dx
            !--- Slope limitation
            slopes4(i,j,k,l) = SIGN(MIN(ABS(slopes4(i,j,k,l)), &
              ABS(fp-ff)/h(i),ABS(ff-fm)/h(i)),slopes4(i,j,k,l) )
           END IF
        END DO
      END DO
    END DO
  END DO
  DEALLOCATE(xc,h,delta_xc)

END FUNCTION slopes4
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
PURE FUNCTION slopes5(ix, f, x)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN) :: ix
  REAL,    INTENT(IN) :: f(:, :, :, :, :)
  REAL,    INTENT(IN) :: x(:)
  REAL :: slopes5(SIZE(f,1),SIZE(f,2),SIZE(f,3),SIZE(f,4),SIZE(f,5))
!-------------------------------------------------------------------------------
! Local:
  INTEGER :: n, i, j, k, l, m, sta(5), sto(5)
  REAL, ALLOCATABLE :: xc(:)                        ! (n) cell centers
  REAL, ALLOCATABLE :: h(:), delta_xc(:)            ! (2:n-1)
  REAL :: fm, ff, fp, dx
!-------------------------------------------------------------------------------
  n=SIZE(f,ix); ALLOCATE(xc(n),h(2:n-1),delta_xc(2:n-1))
  FORALL(i=1:n) xc(i)=(x(i)+x(i+1))/2.
  FORALL(i=2:n-1)
    h(i)=ABS(x(i+1)-xc(i)) ; delta_xc(i)=xc(i+1)-xc(i-1)
  END FORALL
  slopes5(:,:,:,:,:)=0.
  sta=[1,1,1,1,1]; sta(ix)=2
  sto=SHAPE(f);    sto(ix)=n-1
  DO m=sta(5),sto(5); IF(ix==5) dx=delta_xc(m)
    DO l=sta(4),sto(4); IF(ix==4) dx=delta_xc(l)
      DO k=sta(3),sto(3); IF(ix==3) dx=delta_xc(k)
        DO j=sta(2),sto(2); IF(ix==2) dx=delta_xc(j)
          DO i=sta(1),sto(1); IF(ix==1) dx=delta_xc(i)
            ff=f(i,j,k,l,m)
            SELECT CASE(ix)
              CASE(1); fm=f(i-1,j,k,l,m); fp=f(i+1,j,k,l,m)
              CASE(2); fm=f(i,j-1,k,l,m); fp=f(i,j+1,k,l,m)
              CASE(3); fm=f(i,j,k-1,l,m); fp=f(i,j,k+1,l,m)
              CASE(4); fm=f(i,j,k,l-1,m); fp=f(i,j,k,l+1,m)
              CASE(5); fm=f(i,j,k,l,m-1); fp=f(i,j,k,l,m+1)
            END SELECT
            IF(ff>=MAX(fm,fp).OR.ff<=MIN(fm,fp)) THEN
              slopes5(i,j,k,l,m)=0.; CYCLE           !--- Local extremum
              !--- 2nd order slope ; (fm, ff, fp) strictly monotonous
              slopes5(i,j,k,l,m)=(fp-fm)/dx
              !--- Slope limitation
              slopes5(i,j,k,l,m) = SIGN(MIN(ABS(slopes5(i,j,k,l,m)), &
                ABS(fp-ff)/h(i),ABS(ff-fm)/h(i)),slopes5(i,j,k,l,m) )
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
  DEALLOCATE(xc,h,delta_xc)

END FUNCTION slopes5
!
!-------------------------------------------------------------------------------

END MODULE slopes_m
