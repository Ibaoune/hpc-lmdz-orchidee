!
MODULE screenc_mod
!
! This module contains some procedures for calculation of the correction
! of temperature, specific humidity and wind at a reference level
!
  USE cdrag_mod
  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!
!r original routine svn3623
!
      SUBROUTINE screenc(klon, knon, nsrf, zxli, &
                         speed, temp, q_zref, zref, &
                         ts, qsurf, z0m, z0h, psol, &
                         ustar, testar, qstar, okri, ri1, &
                         pref, delu, delte, delq)
      IMPLICIT NONE
!-----------------------------------------------------------------------
! 
! Objet : calcul "correcteur" des anomalies du vent, de la temperature 
!         potentielle et de l'humidite relative au niveau de reference zref et 
!         par rapport au 1er niveau (pour u) ou a la surface (pour theta et q) 
!         a partir des equations de Louis.
!
! Reference : Hess, Colman et McAvaney (1995)
!
! I. Musat, 01.07.2002
!-----------------------------------------------------------------------
!
! klon----input-I- dimension de la grille physique (= nb_pts_latitude X nb_pts_longitude)
! knon----input-I- nombre de points pour un type de surface
! nsrf----input-I- indice pour le type de surface; voir indice_sol_mod.F90
! zxli----input-L- TRUE si calcul des cdrags selon Laurent Li
! speed---input-R- module du vent au 1er niveau du modele
! temp----input-R- temperature de l'air au 1er niveau du modele
! q_zref--input-R- humidite relative au 1er niveau du modele
! zref----input-R- altitude de reference
! ts------input-R- temperature de l'air a la surface
! qsurf---input-R- humidite relative a la surface
! z0m, z0h---input-R- rugosite
! psol----input-R- pression au sol
! ustar---input-R- facteur d'echelle pour le vent
! testar--input-R- facteur d'echelle pour la temperature potentielle
! qstar---input-R- facteur d'echelle pour l'humidite relative
! okri----input-L- TRUE si on veut tester le nb. Richardson entre la sfce 
!                  et zref par rapport au Ri entre la sfce et la 1ere couche
! ri1-----input-R- nb. Richardson entre la surface et la 1ere couche 
!
! pref----input-R- pression au niveau de reference
! delu----input-R- anomalie du vent par rapport au 1er niveau
! delte---input-R- anomalie de la temperature potentielle par rapport a la surface
! delq----input-R- anomalie de l'humidite relative par rapport a la surface
!
      INTEGER, intent(in) :: klon, knon, nsrf
      LOGICAL, intent(in) :: zxli, okri 
      REAL, dimension(klon), intent(in) :: speed, temp, q_zref
      REAL, intent(in) :: zref
      REAL, dimension(klon), intent(in) :: ts, qsurf, z0m, z0h, psol
      REAL, dimension(klon), intent(in) :: ustar, testar, qstar, ri1         
!
      REAL, dimension(klon), intent(out) :: pref, delu, delte, delq 
!-----------------------------------------------------------------------
      include "YOMCST.h"
      include "flux_arp.h"
!
! Variables locales  
      INTEGER :: i 
      REAL, dimension(klon) :: cdram, cdrah, cdran, zri1, gref,ycdragm
!
!------------------------------------------------------------------------- 
      DO i=1, knon
        gref(i) = zref*RG
      ENDDO 
!
! Richardson at reference level 
!
!      CALL coefcdrag (klon, knon, nsrf, zxli, &
!                    speed, temp, q_zref, gref, &
!                    psol, ts, qsurf, rugos, &
!                    okri, ri1, &
!                    cdram, cdrah, cdran, zri1, &
!                    pref)
! Fuxing WANG, 04/03/2015, replace the coefcdrag by the merged version: cdrag
      CALL cdrag (knon, nsrf, &
                    speed, temp, q_zref, gref, &
                    psol, ts, qsurf, z0m, z0h, &
                    cdram, cdrah, zri1, pref)
      DO i = 1, knon
        IF(ok_prescr_ust) THEN
! La aussi il faut forcer avec ust (FC + MPL 20160210)
        ycdragm(i) = ust*ust/(1.+speed(i))/speed(i)
        cdram=ycdragm
        delu(i) = ust/sqrt(cdram(i))
        ELSE
        delu(i) = ustar(i)/sqrt(cdram(i))
        ENDIF
        delte(i)= (testar(i)* sqrt(cdram(i)))/ &
                   cdrah(i)
        delq(i)= (qstar(i)* sqrt(cdram(i)))/ &
                  cdrah(i)
      ENDDO 
!
      RETURN 
      END SUBROUTINE screenc
!
      SUBROUTINE screencn(klon, knon, nsrf, zxli, &
                         speed, temp, q_zref, zref, &
                         ts, qsurf, z0m, z0h, psol, &
                         cdrm, cdrh,  okri, &
                         ri1, iri1, &
                         pref, delm, delh, zri1)
      IMPLICIT NONE
!-----------------------------------------------------------------------
! 
! Objet : calcul "correcteur" des anomalies du vent, de la temperature 
!         potentielle et de l'humidite relative au niveau de reference zref et 
!         par rapport au 1er niveau (pour u) ou a la surface (pour theta et q) 
!         a partir des equations de Louis.
!
! Reference : Hess, Colman et McAvaney (1995)
!
! I. Musat, 01.07.2002
!-----------------------------------------------------------------------
!
! klon----input-I- dimension de la grille physique (= nb_pts_latitude X nb_pts_longitude)
! knon----input-I- nombre de points pour un type de surface
! nsrf----input-I- indice pour le type de surface; voir indice_sol_mod.F90
! zxli----input-L- TRUE si calcul des cdrags selon Laurent Li
! speed---input-R- module du vent au 1er niveau du modele
! temp----input-R- temperature de l'air au 1er niveau du modele
! q_zref--input-R- humidite relative au 1er niveau du modele
! zref----input-R- altitude de reference
! ts------input-R- temperature de l'air a la surface
! qsurf---input-R- humidite relative a la surface
! z0m, z0h---input-R- rugosite
! psol----input-R- pression au sol
! ustar---input-R- facteur d'echelle pour le vent
! testar--input-R- facteur d'echelle pour la temperature potentielle
! qstar---input-R- facteur d'echelle pour l'humidite relative
! okri----input-L- TRUE si on veut tester le nb. Richardson entre la sfce 
!                  et zref par rapport au Ri entre la sfce et la 1ere couche
! ri1-----input-R- nb. Richardson entre la surface et la 1ere couche 
!
! pref----input-R- pression au niveau de reference
! delu----input-R- anomalie du vent par rapport au 1er niveau
! delte---input-R- anomalie de la temperature potentielle par rapport a la surface
! delq----input-R- anomalie de l'humidite relative par rapport a la surface
!
      INTEGER, intent(in) :: klon, knon, nsrf
      LOGICAL, intent(in) :: zxli, okri 
      REAL, dimension(klon), intent(in) :: speed, temp, q_zref
      REAL, intent(in) :: zref
      REAL, dimension(klon), intent(in) :: ts, qsurf, z0m, z0h, psol
      REAL, dimension(klon), intent(in) :: cdrm, cdrh, ri1         
      INTEGER, INTENT(IN)  :: iri1 ! Richardson de la 1ere couche
!
      REAL, dimension(klon), intent(out) :: pref, delm, delh, zri1
!171220     REAL, dimension(klon) :: cdram, cdrah, zri1 
      REAL, dimension(klon) :: cdram, cdrah
!-----------------------------------------------------------------------
      include "YOMCST.h"
      include "flux_arp.h"
!
! Variables locales  
      INTEGER :: i 
      REAL, dimension(klon) :: cdran, gref,ycdragm
!
!------------------------------------------------------------------------- 
      DO i=1, knon
        gref(i) = zref*RG
      ENDDO 
!
! Richardson at reference level 
!
      CALL cdragn_ri1 (knon, nsrf, &
                    speed, temp, q_zref, gref, &
                    psol, ts, qsurf, z0m, z0h, &
                    ri1, iri1, &
                    cdram, cdrah, zri1, pref)
      DO i = 1, knon
        delm(i) = sqrt(cdrm(i))/sqrt(cdram(i))
!verifier que temp est un temperat potentielle.
        delh(i)= (cdrh(i)* sqrt(cdram(i)))/ &
                   (cdrah(i)*sqrt(cdrm(i)))
      ENDDO 
!
      RETURN 
      END SUBROUTINE screencn
END MODULE screenc_mod
