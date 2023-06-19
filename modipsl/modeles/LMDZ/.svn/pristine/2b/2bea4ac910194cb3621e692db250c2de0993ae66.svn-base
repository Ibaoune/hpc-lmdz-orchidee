!
! $Header$
!
MODULE icefrac_lsc_mod

IMPLICIT NONE

CONTAINS
!*******************************************************************

SUBROUTINE icefrac_lsc(np,temp, sig, icefrac)
  !
  ! Compute the ice fraction 1-xliq (see e.g.
  ! Doutriaux-Boucher & Quaas 2004, section 2.2.)
  !
  ! (JBM 3/14 8/14 5/16)
  
  USE print_control_mod, ONLY: lunout, prt_level
  INCLUDE "nuage.h"

  ! nuage.h contains:
  ! t_glace_min: if T < Tmin, the cloud is only made of water ice
  ! t_glace_max: if T > Tmax, the cloud is only made of liquid water
  ! exposant_glace: controls the sharpness of the transition
  INTEGER :: np
  REAL, DIMENSION(np), INTENT(IN) :: temp ! temperature
  REAL, DIMENSION(np), INTENT(IN) :: sig
  REAL, DIMENSION(np), INTENT(OUT) :: icefrac

  REAL :: sig0,www,tmin_tmp,liqfrac_tmp
  INTEGER :: ip

  sig0=0.8

  DO ip=1,np
     IF (iflag_t_glace.EQ.1) THEN
       ! Transition to ice close to surface for T<Tmax
       ! w=1 at the surface and 0 for sig < sig0
       www=(max(sig(ip)-sig0,0.))/(1.-sig0) 
     ELSEIF (iflag_t_glace.GE.2) THEN
       ! No convertion to ice close to surface
       www = 0.
     ENDIF
     tmin_tmp=www*t_glace_max+(1.-www)*t_glace_min
     liqfrac_tmp=  (temp(ip)-tmin_tmp) / (t_glace_max-tmin_tmp)
     liqfrac_tmp = MIN(MAX(liqfrac_tmp,0.0),1.0)
     IF (iflag_t_glace.GE.3) THEN
       icefrac(ip) = 1.0-liqfrac_tmp**exposant_glace
     ELSE
       icefrac(ip) = (1.0-liqfrac_tmp)**exposant_glace
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE icefrac_lsc

!*******************************************************************
!
END MODULE icefrac_lsc_mod
