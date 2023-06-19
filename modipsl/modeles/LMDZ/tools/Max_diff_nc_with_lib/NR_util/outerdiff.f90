MODULE outerdiff_m

  IMPLICIT NONE

  INTERFACE outerdiff
     MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
  END INTERFACE

  private outerdiff_r,outerdiff_d,outerdiff_i

CONTAINS

  FUNCTION outerdiff_r(a,b)
    REAL, DIMENSION(:), INTENT(IN) :: a,b
    REAL, DIMENSION(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_r
  !BL
  FUNCTION outerdiff_d(a,b)
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: a,b
    DOUBLE PRECISION, DIMENSION(size(a),size(b)) :: outerdiff_d
    outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_d
  !BL
  FUNCTION outerdiff_i(a,b)
    INTEGER, DIMENSION(:), INTENT(IN) :: a,b
    INTEGER, DIMENSION(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i

END MODULE outerdiff_m
