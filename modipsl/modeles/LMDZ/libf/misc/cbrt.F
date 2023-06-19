      FUNCTION cbrt(x)
      ! Return the cubic root for positive or negative x
      IMPLICIT NONE

      REAL x,cbrt

      cbrt=sign(1.,x)*(abs(x)**(1./3.))

      RETURN
      END

