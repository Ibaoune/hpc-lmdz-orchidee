SUBROUTINE diag_slp(nlon,t,pab,pal,pphis,tasfc,tastd,pmer)
 USE dimphy
 USE phys_output_write_mod
 USE phys_output_ctrlout_mod
 USE phys_local_var_mod
 IMPLICIT NONE
  !>======================================================================
  !!
  !! Auteur(s) I.Musat (LMD/CNRS) date: 20151106
  !!
  !! Objet: Calcul pression au niveau de la mer cf. Arpege-IFS 
  !! ctstar: calcule la temperature de l'air a la surface (tasfc) et
  !!                 la temperature de l'air standard a la surface (tastd)
  !! pppmer: calcule la slp a partir de tasfc, tastd, de la pression a la surface (pab1)
  !!                 et du geopotentiel de la surface
  !!======================================================================
  !! nlon--input-R-temperature au milieu de chaque couche (en K)
  !! t--input-R-temperature au milieu de chaque couche (en K)
  !! pab--input-R-pression pour chaque inter-couche (en Pa)
  !! pal---input-R-pression pour le mileu de chaque couche (en Pa)
  !! pphis---input-R-geopotentiel du sol (en m2/s2)
  !! tasfc---output-R-temperature air au sol (en K)
  !! tastd---output-R-temperature air 'standard' au sol (en K)
  !! pmer---output-R-pression au niveau de la mer (en Pa)
  !!======================================================================
  INTEGER nlon
  REAL t(nlon,klev)
  REAL pab(nlon,klev+1)
  REAL pal(nlon,klev)
  REAL pphis(nlon)
  REAL tasfc(nlon), tastd(nlon)
  REAL pmer(nlon)
!!!
!!! calcul tasfc et tastd
  tal1(:)=t(:,1)
  pab2(:)=pab(:,2)
  pal1(:)=pal(:,1)
  CALL ctstar(nlon,1,nlon,tal1,pab2,pal1,pphis,tasfc,tastd)
!!!
!!! calcul slp
  pab1(:)=pab(:,1)
  CALL pppmer(nlon,1,nlon,pab1,pphis,tasfc,tastd,pmer)
!
  RETURN
  END
