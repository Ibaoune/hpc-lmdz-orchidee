
! $Header$

SUBROUTINE ustarhb(knon, u, v, cd_m, ustar)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Laurent Li (LMD/CNRS), le 30 septembre 1998
  ! Couche limite non-locale. Adaptation du code du CCM3.
  ! Code non teste, donc a ne pas utiliser.
  ! ======================================================================
  ! Nonlocal scheme that determines eddy diffusivities based on a
  ! diagnosed boundary layer height and a turbulent velocity scale.
  ! Also countergradient effects for heat and moisture are included.

  ! For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
  ! Local versus nonlocal boundary-layer diffusion in a global climate
  ! model. J. of Climate, vol. 6, 1825-1842.
  ! ======================================================================
  include "YOMCST.h"

  ! Arguments:

  INTEGER knon ! nombre de points a calculer
  REAL u(klon, klev) ! vitesse U (m/s)
  REAL v(klon, klev) ! vitesse V (m/s)
  REAL cd_m(klon) ! coefficient de friction au sol pour vitesse
  REAL ustar(klon)

  INTEGER i, k
  REAL zxt, zxq, zxu, zxv, zxmod, taux, tauy
  REAL zx_alf1, zx_alf2 ! parametres pour extrapolation
  LOGICAL unssrf(klon) ! unstb pbl w/lvls within srf pbl lyr
  LOGICAL unsout(klon) ! unstb pbl w/lvls in outer pbl lyr
  LOGICAL check(klon) ! True=>chk if Richardson no.>critcal

  include "YOETHF.h"
  include "FCTTRE.h"
  DO i = 1, knon
    zx_alf1 = 1.0
    zx_alf2 = 1.0 - zx_alf1
    zxu = u(i, 1)*zx_alf1 + u(i, 2)*zx_alf2
    zxv = v(i, 1)*zx_alf1 + v(i, 2)*zx_alf2
    zxmod = 1.0 + sqrt(zxu**2+zxv**2)
    taux = zxu*zxmod*cd_m(i)
    tauy = zxv*zxmod*cd_m(i)
    ustar(i) = sqrt(taux**2+tauy**2)
    ! print*,'Ust ',zxu,zxmod,taux,ustar(i)
  END DO

  RETURN
END SUBROUTINE ustarhb
