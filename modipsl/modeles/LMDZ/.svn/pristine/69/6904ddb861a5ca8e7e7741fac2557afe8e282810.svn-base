MODULE cumsum_m

  IMPLICIT NONE

  INTEGER, PARAMETER, private :: NPAR_CUMSUM=16

  INTERFACE cumsum
     MODULE PROCEDURE cumsum_r,cumsum_i
  END INTERFACE

  private cumsum_r,cumsum_i

CONTAINS

  RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
    REAL, DIMENSION(:), INTENT(IN) :: arr
    REAL, OPTIONAL, INTENT(IN) :: seed
    REAL, DIMENSION(size(arr)) :: ans
    INTEGER :: n,j
    REAL :: sd
    n=size(arr)
    if (n == 0) RETURN
    sd=0.0
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum_r
  !BL
  RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
    INTEGER, DIMENSION(:), INTENT(IN) :: arr
    INTEGER, OPTIONAL, INTENT(IN) :: seed
    INTEGER, DIMENSION(size(arr)) :: ans
    INTEGER :: n,j,sd
    n=size(arr)
    if (n == 0) RETURN
    sd=0
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum_i

END MODULE cumsum_m
