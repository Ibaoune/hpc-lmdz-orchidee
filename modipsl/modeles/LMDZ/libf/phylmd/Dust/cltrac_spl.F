      SUBROUTINE cltrac_spl(dtime,coef,yu1,yv1,t,tr,
     .                     flux,paprs,pplay,d_tr)

      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): O. Boucher (LOA/LMD) date: 19961127
c            inspire de clvent
c Objet: diffusion verticale de traceurs avec flux fixe a la surface
c        ou/et flux du type c-drag
c======================================================================
c Arguments:
c dtime----input-R- intervalle du temps (en second)
c coef-----input-R- le coefficient d'echange (m**2/s) l>1
c yu1------input-R- le vent dans le 1iere couche
c yv1------input-R- le vent dans le 1iere couche
c t--------input-R- temperature (K)
c tr-------input-R- la q. de traceurs
c flux-----input-R- le flux de traceurs a la surface
c paprs----input-R- pression a inter-couche (Pa)
c pplay----input-R- pression au milieu de couche (Pa)
c delp-----input-R- epaisseur de couche (Pa)
c cdrag----input-R- cdrag pour le flux de surface (non active)
c tr0------input-R- traceurs a la surface ou dans l'ocean (non active)
c d_tr-----output-R- le changement de tr
c flux_tr--output-R- flux de tr
c======================================================================
#include "dimensions.h"
      REAL dtime
      REAL coef(klon,klev)
      REAL yu1(klon), yv1(klon)
      REAL t(klon,klev), tr(klon,klev)
      REAL paprs(klon,klev+1), pplay(klon,klev), delp(klon,klev)
      REAL d_tr(klon,klev)
      REAL flux(klon), cdrag(klon), tr0(klon)
c      REAL flux_tr(klon,klev)
c======================================================================
#include "YOMCST.h"
c======================================================================
      INTEGER i, k
      REAL zx_ctr(klon,2:klev)
      REAL zx_dtr(klon,2:klev)
      REAL zx_buf(klon)
      REAL zx_coef(klon,klev)
      REAL local_tr(klon,klev)
      REAL zx_alf1(klon), zx_alf2(klon), zx_flux(klon)
c======================================================================
c CHECKING VALUES
!      print *,'CHECKING VALUES IN CLTRAC (INI)'
!      print *,'d_tr = ',sum(d_tr),MINVAL(d_tr),MAXVAL(d_tr)
!      print *,'flux = ',sum(flux),MINVAL(flux),MAXVAL(flux)
!      print *,'tr = ',sum(tr),MINVAL(tr),MAXVAL(tr)
c======================================================================
      DO k = 1, klev
      DO i = 1, klon
         local_tr(i,k) = tr(i,k)
         delp(i,k) = paprs(i,k)-paprs(i,k+1)
      ENDDO
      ENDDO
c======================================================================
      DO i = 1, klon
         zx_alf1(i) = (paprs(i,1)-pplay(i,2))/(pplay(i,1)-pplay(i,2))
         zx_alf2(i) = 1.0 - zx_alf1(i)
         zx_flux(i) =  -flux(i)*dtime*RG
c--pour le moment le flux est prescrit
         cdrag(i) = 0.0 
c         cdrag(i) =  coef(i,1) * (1.0+SQRT(yu1(i)**2+yv1(i)**2))
c     .                * pplay(i,1)/(RD*t(i,1))
         tr0(i) = 0.0
         zx_coef(i,1) = cdrag(i)*dtime*RG 
      ENDDO
c======================================================================
      DO k = 2, klev
      DO i = 1, klon
         zx_coef(i,k) = coef(i,k)*RG/(pplay(i,k-1)-pplay(i,k))
     .                  *(paprs(i,k)*2/(t(i,k)+t(i,k-1))/RD)**2
         zx_coef(i,k) = zx_coef(i,k)*dtime*RG
      ENDDO
      ENDDO
c======================================================================
      DO i = 1, klon
         zx_buf(i) = delp(i,1) + zx_coef(i,1)*zx_alf1(i) + zx_coef(i,2)
         zx_ctr(i,2) = (local_tr(i,1)*delp(i,1)+
     .                  zx_coef(i,1)*tr0(i)-zx_flux(i))/zx_buf(i)
         zx_dtr(i,2) = (zx_coef(i,2)-zx_alf2(i)*zx_coef(i,1)) / 
     .                  zx_buf(i)
      ENDDO
c
      DO k = 3, klev
      DO i = 1, klon
         zx_buf(i) = delp(i,k-1) + zx_coef(i,k)
     .                  + zx_coef(i,k-1)*(1.-zx_dtr(i,k-1))
         zx_ctr(i,k) = (local_tr(i,k-1)*delp(i,k-1)
     .                  +zx_coef(i,k-1)*zx_ctr(i,k-1) )/zx_buf(i)
         zx_dtr(i,k) = zx_coef(i,k)/zx_buf(i)
      ENDDO
      ENDDO
      DO i = 1, klon
         local_tr(i,klev) = ( local_tr(i,klev)*delp(i,klev)
     .                        +zx_coef(i,klev)*zx_ctr(i,klev) )
     .                   / ( delp(i,klev) + zx_coef(i,klev)
     .                       -zx_coef(i,klev)*zx_dtr(i,klev) )
      ENDDO
      DO k = klev-1, 1, -1
      DO i = 1, klon
         local_tr(i,k) = zx_ctr(i,k+1) + zx_dtr(i,k+1)*local_tr(i,k+1)
      ENDDO
      ENDDO
c======================================================================
!      print *,'CHECKING VALUES IN CLTRAC (FIN)'
!      print *,'local_tr = ',sum(local_tr),MINVAL(local_tr),
!     .                                    MAXVAL(local_tr)
!      print *,'zx_ctr = ',sum(zx_ctr),MINVAL(zx_ctr),MAXVAL(zx_ctr)
!      print *,'zx_dtr = ',sum(zx_dtr),MINVAL(zx_dtr),MAXVAL(zx_dtr)
!      print *,'tr = ',sum(tr),MINVAL(tr),MAXVAL(tr)
c======================================================================
c== flux_tr est le flux de traceur (positif vers bas)
c      DO i = 1, klon
c         flux_tr(i,1) = zx_coef(i,1)/(RG*dtime)
c      ENDDO
c      DO k = 2, klev
c      DO i = 1, klon
c         flux_tr(i,k) = zx_coef(i,k)/(RG*dtime)
c     .               * (local_tr(i,k)-local_tr(i,k-1))
c      ENDDO
c      ENDDO
c======================================================================
      DO k = 1, klev
      DO i = 1, klon
         d_tr(i,k) = local_tr(i,k) - tr(i,k)
      ENDDO
      ENDDO
!      print *,'CHECKING VALUES IN CLTRAC (END)'
!      print *,'d_tr = ',sum(d_tr),MINVAL(d_tr),MAXVAL(d_tr)
c
      RETURN
      END
