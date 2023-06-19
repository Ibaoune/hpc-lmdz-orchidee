MODULE swap_m

  IMPLICIT NONE

  INTERFACE swap
     MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
          swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
          masked_swap_rs,masked_swap_rv,masked_swap_rm
  END INTERFACE

  private
  public swap

CONTAINS

  SUBROUTINE swap_i(a,b)
    INTEGER, INTENT(INOUT) :: a,b
    INTEGER :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i

  !************************************************

  SUBROUTINE swap_r(a,b)
    REAL, INTENT(INOUT) :: a,b
    REAL :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r

  !************************************************

  SUBROUTINE swap_rv(a,b)
    REAL, DIMENSION(:), INTENT(INOUT) :: a,b
    REAL, DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv

  !************************************************

  SUBROUTINE swap_c(a,b)
    COMPLEX, INTENT(INOUT) :: a,b
    COMPLEX :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_c

  !************************************************

  SUBROUTINE swap_cv(a,b)
    COMPLEX, DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX, DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cv

  !************************************************

  SUBROUTINE swap_cm(a,b)
    COMPLEX, DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX, DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cm

  !************************************************

  SUBROUTINE swap_z(a,b)
    COMPLEX(KIND(0D0)), INTENT(INOUT) :: a,b
    COMPLEX(KIND(0D0)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_z

  !************************************************

  SUBROUTINE swap_zv(a,b)
    COMPLEX(KIND(0D0)), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(KIND(0D0)), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zv

  !************************************************

  SUBROUTINE swap_zm(a,b)
    COMPLEX(KIND(0D0)), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(KIND(0D0)), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zm

  !************************************************

  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL, INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    REAL :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs

  !************************************************

  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL, DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    REAL, DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rv

  !************************************************

  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL, DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    REAL, DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rm

END MODULE swap_m
