c       This subroutine estimateis Sea Salt emission fluxes over 
c       Oceanic surfaces.
c 
      SUBROUTINE seasalt(v_10m, u_10m, pct_ocean, lmt_sea_salt)
    
      USE dimphy
      IMPLICIT NONE 
c
#include "dimensions.h"
#include "chem.h"
#include "chem_spla.h"
#include "YOMCST.h"
#include "YOECUMF.h"
c
      INTEGER i, bin                 !local variables
      REAL pct_ocean(klon)           !hfraction of Ocean in each grid
      REAL v_10m(klon), u_10m(klon)  !V&H components of wind @10 m 
      REAL w_speed_10m(klon)         !wind speed at 10m from surface
      REAL lmt_sea_salt(klon,ss_bins)!sea salt emission flux - mg/m2/s
      REAL sea_salt_flux(ss_bins)    !sea salt emission flux per unit wind speed

      REAL wind, ocean
c
c------Sea salt emission fluxes for each size bin calculated 
c------based on on parameterisation of Gong et al. (1997).
c------Fluxes of sea salt for each size bin are given in mg/m^2/sec
c------at wind speed of 1 m/s at 10m height (at 80% RH).
c------Fluxes at various wind speeds (@10 m from sea 
c------surfaces are estimated using relationship: F=flux*U_10^3.14
c
cnhl for size bin of 0.03-0.5 and 0.5-20
      DATA sea_salt_flux/4.5E-09,8.7E-7/

      DO i=1, klon
      w_speed_10m(i)= (v_10m(i)**2.0+u_10m(i)**2.0)**0.5
      ENDDO
c
      DO bin=1,ss_bins
      wind=0.0
      ocean=0.0
      DO i=1, klon
      lmt_sea_salt(i,bin)=sea_salt_flux(bin)*(w_speed_10m(i)**3.41)
     . *pct_ocean(i)*1.e-4*1.e-3                       !g/cm2/s
      wind=wind+w_speed_10m(i)
      ocean=ocean+pct_ocean(i)
      ENDDO
!      print *,'Sea Salt flux = ',sea_salt_flux(bin)
      ENDDO
!      print *,'SUM OF WIND = ',wind
!      print *,'SUM OF OCEAN SURFACE = ',ocean
      RETURN 	
      END
