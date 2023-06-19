      SUBROUTINE nightingale(u, v, u_10m, v_10m, paprs, pplay, 
     .                       cdragh, cdragm, t, q, ftsol, tsol, 
     .                       pctsrf, lmt_dmsconc, lmt_dms)
c
      USE dimphy 
      USE indice_sol_mod
      IMPLICIT NONE 
c
#include "dimensions.h"
#include "YOMCST.h"
c
      REAL u(klon,klev), v(klon,klev)
      REAL u_10m(klon), v_10m(klon)
      REAL ftsol(klon,nbsrf)
      REAL tsol(klon)
      REAL paprs(klon,klev+1), pplay(klon,klev)
      REAL t(klon,klev)
      REAL q(klon,klev)
      REAL cdragh(klon), cdragm(klon)
      REAL pctsrf(klon,nbsrf)
      REAL lmt_dmsconc(klon)  ! concentration oceanique DMS
      REAL lmt_dms(klon)      ! flux de DMS
c
      REAL ustar(klon), obklen(klon)
      REAL u10(klon), u10n(klon)
      REAL tvelocity, schmidt_corr
      REAL t1, t2, t3, t4, viscosity_kin, diffusivity, schmidt
      INTEGER i
c
      CALL bl_for_dms(u, v, paprs, pplay, cdragh, cdragm,
     .                t, q, tsol, ustar, obklen)
c
      DO i=1,klon
        u10(i)=SQRT(u_10m(i)**2+v_10m(i)**2)
      ENDDO
c
      CALL neutral(u10, ustar, obklen, u10n)
c
      DO i=1,klon
c
c	tvelocity - transfer velocity, also known as kw (cm/s)
c	schmidt_corr - Schmidt number correction factor (dimensionless)
c Reference:  Nightingale, P.D., G. Malin, C. S. Law, J. J. Watson, P.S. Liss
c  M. I. Liddicoat, J. Boutin, R.C. Upstill-Goddard. 'In situ evaluation 
c  of air-sea gas exchange parameterizations using conservative and 
c  volatile tracers.'  Glob. Biogeochem. Cycles, 14:373-387, 2000.
c compute transfer velocity using u10neutral     
c
      tvelocity = 0.222*u10n(i)*u10n(i) + 0.333*u10n(i)
c
c above expression gives tvelocity in cm/hr. convert to cm/s. 1hr =3600 sec

      tvelocity = tvelocity / 3600. 	

c compute the correction factor, which for Nightingale parameterization is 
c based on how different the schmidt number is from 600.   
c correction factor based on temperature in Kelvin. good
c only for t<=30 deg C.  for temperatures above that, set correction factor
c equal to value at 30 deg C.

      IF (ftsol(i,is_oce) .LE. 303.15) THEN
         t1 = ftsol(i,is_oce) 
      ELSE
         t1 = 303.15
      ENDIF        

      t2 = t1 * t1
      t3 = t2 * t1
      t4 = t3 * t1
      viscosity_kin = 3.0363e-9*t4 - 3.655198e-6*t3 + 1.65333e-3*t2 
     +      - 3.332083e-1*t1 + 25.26819
      diffusivity = 0.01922 * exp(-2177.1/t1)
      schmidt = viscosity_kin / diffusivity
      schmidt_corr = (schmidt/600.)**(-.5)          	 
c
      lmt_dms(i) = tvelocity  *  pctsrf(i,is_oce) 
     .        * lmt_dmsconc(i)/1.0e12 * schmidt_corr * RNAVO
c
      IF (lmt_dmsconc(i).LE.1.e-20) lmt_dms(i)=0.0
c
      ENDDO
c
      END
