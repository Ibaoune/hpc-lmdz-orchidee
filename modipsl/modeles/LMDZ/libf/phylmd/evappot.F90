SUBROUTINE evappot(klon,nbsrf,ftsol,pplay,cdragh,  &
       &    t_seri,q_seri,u_seri,v_seri,evap_pot)

IMPLICIT NONE

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"


INTEGER :: klon, nbsrf
REAL, DIMENSION(klon,nbsrf) :: ftsol,evap_pot
REAL, DIMENSION(klon) :: pplay,t_seri,wind,q_seri,u_seri,v_seri,cdragh

INTEGER :: nsrf,i
REAL, DIMENSION(klon,nbsrf) :: qsat_ftsol
REAL, DIMENSION(klon) :: rhos, norme_u
REAL :: t_coup

      t_coup=234.   ! Quelle horreur !!!!!

DO nsrf = 1, nbsrf
   DO i = 1, klon
      IF (ftsol(i,nsrf).LT.t_coup) THEN
         qsat_ftsol(i,nsrf) = qsats(ftsol(i,nsrf))/pplay(i)
      ELSE
         qsat_ftsol(i,nsrf) = qsatl(ftsol(i,nsrf))/pplay(i)
      ENDIF
   ENDDO
ENDDO
! ========================================================== c
! Calcul de l'evaporation Potentielle


rhos(:) = pplay(:)/(RD*t_seri(:))
norme_u(:)=1.+sqrt(u_seri(:)*u_seri(:)+v_seri(:)*v_seri(:))
DO nsrf = 1, nbsrf
  evap_pot(:,nsrf)=rhos(:)*cdragh(:)*norme_u(:)*(qsat_ftsol(:,nsrf)-q_seri(:))
ENDDO
RETURN

END
