MODULE array_copy_m

  IMPLICIT NONE

  INTERFACE array_copy
     MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
  END INTERFACE

  private array_copy_r, array_copy_d, array_copy_i

CONTAINS

  SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
    REAL, DIMENSION(:), INTENT(IN) :: src
    REAL, DIMENSION(:), INTENT(OUT) :: dest
    INTEGER, INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_r
  !BL
  SUBROUTINE array_copy_d(src,dest,n_copied,n_not_copied)
    double precision, DIMENSION(:), INTENT(IN) :: src
    double precision, DIMENSION(:), INTENT(OUT) :: dest
    INTEGER, INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_d
  !BL
  SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
    INTEGER, DIMENSION(:), INTENT(IN) :: src
    INTEGER, DIMENSION(:), INTENT(OUT) :: dest
    INTEGER, INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_i

END MODULE array_copy_m
