      SUBROUTINE bl_for_dms(u,v,paprs,pplay,cdragh,cdragm
     .                     ,t,q,tsol,ustar,obklen)
      USE dimphy
      IMPLICIT NONE
c
c===================================================================
c Auteur : E. Cosme 
c Calcul de la vitesse de friction (ustar) et de la longueur de 
c Monin-Obukhov (obklen), necessaires pour calculer les flux de DMS
c par la methode de Nightingale.
c Cette subroutine est plus que fortement inspiree de la subroutine
c 'nonlocal' dans clmain.F .
c reference :  Holtslag, A.A.M., and B.A. Boville, 1993:
c Local versus nonlocal boundary-layer diffusion in a global climate
c model. J. of Climate, vol. 6, 1825-1842. (a confirmer)
c 31 08 01
c===================================================================
c
#include "dimensions.h"
#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
c
c Arguments :
      REAL u(klon,klev)          ! vent zonal
      REAL v(klon,klev)          ! vent meridien
      REAL paprs(klon,klev+1)    ! niveaux de pression aux intercouches (Pa)
      REAL pplay(klon,klev)      ! niveaux de pression aux milieux... (Pa)
      REAL cdragh(klon)          ! coefficient de trainee pour la chaleur
      REAL cdragm(klon)          ! coefficient de trainee pour le vent
      REAL t(klon,klev)          ! temperature
      REAL q(klon,klev)          ! humidite kg/kg
      REAL tsol(klon)            ! temperature du sol
      REAL ustar(klon)           ! vitesse de friction
      REAL obklen(klon)          ! longueur de Monin-Obukhov
c
c Locales :
      REAL vk
      PARAMETER (vk=0.35)
      REAL beta  ! coefficient d'evaporation reelle (/evapotranspiration)
                 ! entre 0 et 1, mais 1 au-dessus de la mer
      PARAMETER (beta=1.)
      INTEGER i,k
      REAL zxt, zxu, zxv, zxq, zxqs, zxmod, taux, tauy
      REAL zcor, zdelta, zcvm5
      REAL z(klon,klev)
      REAL zx_alf1, zx_alf2 ! parametres pour extrapolation
      REAL khfs(klon)       ! surface kinematic heat flux [mK/s]
      REAL kqfs(klon)       ! sfc kinematic constituent flux [m/s]
      REAL heatv(klon)      ! surface virtual heat flux

      
c
c======================================================================
c
c Calculer les hauteurs de chaque couche
c
! JE20150707      r2es=611.14 *18.0153/28.9644
      DO i = 1, klon
         z(i,1) = RD * t(i,1) / (0.5*(paprs(i,1)+pplay(i,1)))
     .               * (paprs(i,1)-pplay(i,1)) / RG
      ENDDO
      DO k = 2, klev
      DO i = 1, klon
         z(i,k) = z(i,k-1)
     .              + RD * 0.5*(t(i,k-1)+t(i,k)) / paprs(i,k)
     .                   * (pplay(i,k-1)-pplay(i,k)) / RG
      ENDDO
      ENDDO

      DO i = 1, klon
c
        zdelta=MAX(0.,SIGN(1.,RTT-tsol(i)))
        zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
        zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*q(i,1))
        zxqs= r2es * FOEEW(tsol(i),zdelta)/paprs(i,1)
        zxqs=MIN(0.5,zxqs)
        zcor=1./(1.-retv*zxqs)
        zxqs=zxqs*zcor
c
        zx_alf1 = 1.0
        zx_alf2 = 1.0 - zx_alf1
        zxt = (t(i,1)+z(i,1)*RG/RCPD/(1.+RVTMP2*q(i,1)))
     .        *(1.+RETV*q(i,1))*zx_alf1
     .      + (t(i,2)+z(i,2)*RG/RCPD/(1.+RVTMP2*q(i,2)))
     .        *(1.+RETV*q(i,2))*zx_alf2
        zxu = u(i,1)*zx_alf1+u(i,2)*zx_alf2
        zxv = v(i,1)*zx_alf1+v(i,2)*zx_alf2
        zxq = q(i,1)*zx_alf1+q(i,2)*zx_alf2
        zxmod = 1.0+SQRT(zxu**2+zxv**2)
        khfs(i) = (tsol(i)*(1.+RETV*q(i,1))-zxt) *zxmod*cdragh(i)
        kqfs(i) = (zxqs-zxq) *zxmod*cdragh(i) * beta
        heatv(i) = khfs(i) + 0.61*zxt*kqfs(i)
        taux = zxu *zxmod*cdragm(i)
        tauy = zxv *zxmod*cdragm(i)
        ustar(i) = SQRT(taux**2+tauy**2)
        ustar(i) = MAX(SQRT(ustar(i)),0.01)
c
      ENDDO
c
      DO i = 1, klon
         obklen(i) = -t(i,1)*ustar(i)**3/(RG*vk*heatv(i))
      ENDDO
c
      END SUBROUTINE
