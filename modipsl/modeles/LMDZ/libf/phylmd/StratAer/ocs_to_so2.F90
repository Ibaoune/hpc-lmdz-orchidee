!
! $Id: ocs_to_so2.F90 3677 2020-05-06 15:18:32Z oboucher $
!
SUBROUTINE ocs_to_so2(pdtphys,tr_seri,t_seri,pplay,paprs,is_strato)

  USE dimphy, ONLY : klon,klev
  USE aerophys
  USE infotrac_phy
  USE YOMCST, ONLY : RG
  USE phys_local_var_mod, ONLY : OCS_lifetime, budg_3D_ocs_to_so2, budg_ocs_to_so2 

  IMPLICIT NONE

  !--------------------------------------------------------
  ! transfer variables when calling this routine
  REAL,INTENT(IN)                               :: pdtphys ! Pas d'integration pour la physique (seconde)
  REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT)  :: tr_seri ! Concentration Traceur [U/KgA]
  REAL,DIMENSION(klon,klev),INTENT(IN)          :: t_seri  ! Temperature
  REAL,DIMENSION(klon,klev),INTENT(IN)          :: pplay   ! pression pour le mileu de chaque couche (en Pa)
  REAL,DIMENSION(klon,klev+1),INTENT(IN)        :: paprs   ! pression pour chaque inter-couche (en Pa)
  LOGICAL,DIMENSION(klon,klev),INTENT(IN)       :: is_strato

  ! local variables
  INTEGER                                       :: i,j,k,nb,ilon,ilev

!--convert OCS to SO2
  budg_3D_ocs_to_so2(:,:)=0.0
  budg_ocs_to_so2(:)=0.0

  DO ilon=1, klon
  DO ilev=1, klev
  !only in the stratosphere
  IF (is_strato(ilon,ilev)) THEN
    IF (OCS_lifetime(ilon,ilev).GT.0.0) THEN
      budg_3D_ocs_to_so2(ilon,ilev)=tr_seri(ilon,ilev,id_OCS_strat)*(1.0-exp(-pdtphys/OCS_lifetime(ilon,ilev)))
    ENDIF
    tr_seri(ilon,ilev,id_OCS_strat)=tr_seri(ilon,ilev,id_OCS_strat) - budg_3D_ocs_to_so2(ilon,ilev)
    tr_seri(ilon,ilev,id_SO2_strat)=tr_seri(ilon,ilev,id_SO2_strat) + mSO2mol/mOCSmol*budg_3D_ocs_to_so2(ilon,ilev)
    !convert budget from kg(OCS)/kgA to kg(S)/m2/layer/s for saving as diagnostic
    budg_3D_ocs_to_so2(ilon,ilev)=budg_3D_ocs_to_so2(ilon,ilev)*mSatom/mOCSmol*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG/pdtphys
    budg_ocs_to_so2(ilon)=budg_ocs_to_so2(ilon)+budg_3D_ocs_to_so2(ilon,ilev) 
  ENDIF
  ENDDO
  ENDDO

END SUBROUTINE ocs_to_so2
