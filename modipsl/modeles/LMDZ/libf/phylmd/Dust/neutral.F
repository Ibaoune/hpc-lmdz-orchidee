c*********************************************************************** 
	subroutine neutral(u10_mps,ustar_mps,obklen_m,
     +         u10n_mps	)
c-----------------------------------------------------------------------	
c subroutine to compute u10 neutral wind speed
c inputs
c	u10_mps - wind speed at 10 m (m/s)
c	ustar_mps - friction velocity (m/s)
c	obklen_m - monin-obukhov length scale (m)
c outputs
c	u10n_mps - wind speed at 10 m under neutral conditions (m/s)
c following code assumes reference height Z is 10m, consistent with use
c of u10 and u10_neutral.  If not, code
c should be changed so that constants of 50. and 160. in equations
c below are changed to -5 * Z and -16 * Z respectively.
c Reference:  G. L. Geernaert.  'Bulk parameterizations for the 
c wind stress and heat fluxes,' in Surface Waves and Fluxes, Vol. I,
c Current Theory, Geernaert and W.J. Plant, editors, Kluwer Academic
c Publishers, Boston, MA, 1990.
c subroutine written Feb 2001 by eg chapman
c adapted to LMD-ZT by E. Cosme 310801
c Following Will Shaw (PNL, Seattle) the theory applied for flux
c calculation with the scheme of Nightingale et al. (2000) does not
c hold anymore when -1<obklen<20. In this case, u10n is set to 0,
c so that the transfer velocity  computed in nightingale.F will also
c be 0. The flux is then set to 0.
c----------------------------------------------------------------------		
c
      USE dimphy
#include "dimensions.h"
c
        real u10_mps(klon),ustar_mps(klon),obklen_m(klon)
        real u10n_mps(klon)
        real pi,von_karman
c       parameter (pi = 3.141592653589793, von_karman = 0.4)	
c pour etre coherent avec vk de bl_for_dms.F
        parameter (pi = 3.141592653589793, von_karman = 0.35)	
c
        real phi, phi_inv, phi_inv_sq, f1, f2, f3, dum1, psi
        integer i


        psi = 0.
        do i=1,klon

	if (u10_mps(i) .lt. 0.) u10_mps(i) = 0.0
	
	if  (obklen_m(i) .lt. 0.) then 
		phi = (1. - 160./obklen_m(i))**(-0.25)
		phi_inv = 1./phi
		phi_inv_sq = 1./phi * 1./phi
		f1 = (1. + phi_inv) / 2.
		f2 = (1. + phi_inv_sq)/2.
c following to avoid numerical overruns. recall tan(90deg)=infinity
		dum1 = min (1.e24, phi_inv)
		f3 = atan(dum1)
		psi = 2.*log(f1) + log(f2) - 2.*f3 + pi/2.   
	else if (obklen_m(i) .gt. 0.) then
		psi = -50. / obklen_m(i)
	end if

	u10n_mps(i) = u10_mps(i) + (ustar_mps(i) * psi /von_karman )
c u10n set to 0. if -1 < obklen < 20
        if ((obklen_m(i).gt.-1.).and.(obklen_m(i).lt.20.)) then
            u10n_mps(i) = 0.
        endif
	if (u10n_mps(i) .lt. 0.) u10n_mps(i) = 0.0

        enddo
	return
	end
c***********************************************************************
