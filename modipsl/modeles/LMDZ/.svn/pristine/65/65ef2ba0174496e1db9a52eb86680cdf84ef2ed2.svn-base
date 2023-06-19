
SUBROUTINE coefkzmin(knon, ypaprs, ypplay, yu, yv, yt, yq, ycdragm, km, kn)

  USE dimphy
  IMPLICIT NONE

  include "YOMCST.h"

  ! .......................................................................
  ! Entrees modifies en attendant une version ou les zlev, et zlay soient
  ! disponibles.

  REAL ycdragm(klon)

  REAL yu(klon, klev), yv(klon, klev)
  REAL yt(klon, klev), yq(klon, klev)
  REAL ypaprs(klon, klev+1), ypplay(klon, klev)
  REAL yustar(klon)
  REAL yzlay(klon, klev), yzlev(klon, klev+1), yteta(klon, klev)

  INTEGER i

  ! .......................................................................

  ! En entree :
  ! -----------

  ! zlev : altitude a chaque niveau (interface inferieure de la couche
  ! de meme indice)
  ! ustar : u*

  ! teta : temperature potentielle au centre de chaque couche
  ! (en entree : la valeur au debut du pas de temps)

  ! en sortier :
  ! ------------

  ! km : diffusivite turbulente de quantite de mouvement (au bas de chaque
  ! couche)
  ! (en sortie : la valeur a la fin du pas de temps)
  ! kn : diffusivite turbulente des scalaires (au bas de chaque couche)
  ! (en sortie : la valeur a la fin du pas de temps)

  ! .......................................................................

  REAL ustar(klon)
  REAL kmin, qmin, pblhmin(klon), coriol(klon)
  REAL zlev(klon, klev+1)
  REAL teta(klon, klev)

  REAL km(klon, klev)
  REAL kn(klon, klev)
  INTEGER knon


  INTEGER nlay, nlev
  INTEGER ig, k

  REAL, PARAMETER :: kap = 0.4

  nlay = klev
  nlev = klev + 1
  ! .......................................................................
  ! en attendant une version ou les zlev, et zlay soient
  ! disponibles.
  ! Debut de la partie qui doit etre unclue a terme dans clmain.

  DO i = 1, knon
    yzlay(i, 1) = rd*yt(i, 1)/(0.5*(ypaprs(i,1)+ypplay(i, &
      1)))*(ypaprs(i,1)-ypplay(i,1))/rg
  END DO
  DO k = 2, klev
    DO i = 1, knon
      yzlay(i, k) = yzlay(i, k-1) + rd*0.5*(yt(i,k-1)+yt(i,k))/ypaprs(i, k)*( &
        ypplay(i,k-1)-ypplay(i,k))/rg
    END DO
  END DO
  DO k = 1, klev
    DO i = 1, knon
      ! ATTENTION:on passe la temperature potentielle virt. pour le calcul de
      ! K
      yteta(i, k) = yt(i, k)*(ypaprs(i,1)/ypplay(i,k))**rkappa* &
        (1.+0.61*yq(i,k))
    END DO
  END DO
  DO i = 1, knon
    yzlev(i, 1) = 0.
    yzlev(i, klev+1) = 2.*yzlay(i, klev) - yzlay(i, klev-1)
  END DO
  DO k = 2, klev
    DO i = 1, knon
      yzlev(i, k) = 0.5*(yzlay(i,k)+yzlay(i,k-1))
    END DO
  END DO

  yustar(1:knon) = sqrt(ycdragm(1:knon)*(yu(1:knon,1)*yu(1:knon,1)+yv(1:knon, &
    1)*yv(1:knon,1)))

  ! Fin de la partie qui doit etre unclue a terme dans clmain.

  ! ette routine est ecrite pour avoir en entree ustar, teta et zlev
  ! Ici, on a inclut le calcul de ces trois variables dans la routine
  ! coefkzmin en attendant une nouvelle version de la couche limite
  ! ou ces variables seront disponibles.

  ! Debut de la routine coefkzmin proprement dite.

  ustar = yustar
  teta = yteta
  zlev = yzlev

  DO ig = 1, knon
    coriol(ig) = 1.E-4
    pblhmin(ig) = 0.07*ustar(ig)/max(abs(coriol(ig)), 2.546E-5)
  END DO

  DO k = 2, klev
    DO ig = 1, knon
      IF (teta(ig,2)>teta(ig,1)) THEN
        qmin = ustar(ig)*(max(1.-zlev(ig,k)/pblhmin(ig),0.))**2
        kmin = kap*zlev(ig, k)*qmin
      ELSE
        kmin = 0. ! kmin n'est utilise que pour les SL stables.
      END IF
      kn(ig, k) = kmin
      km(ig, k) = kmin
    END DO
  END DO


  RETURN
END SUBROUTINE coefkzmin
