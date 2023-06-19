!
! $Id: so2_to_h2so4.F90 3677 2020-05-06 15:18:32Z oboucher $
!
SUBROUTINE SO2_TO_H2SO4(pdtphys,tr_seri,t_seri,pplay,paprs,is_strato)

  USE dimphy, ONLY : klon,klev
  USE aerophys
  USE infotrac_phy
  USE YOMCST, ONLY : RG
  USE phys_local_var_mod, ONLY : SO2_lifetime, budg_3D_so2_to_h2so4, budg_so2_to_h2so4

  IMPLICIT NONE

  !--------------------------------------------------------

  ! transfer variables when calling this routine
  REAL,INTENT(IN)                               :: pdtphys ! Pas d'integration pour la physique (seconde)
  REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT)  :: tr_seri ! Concentration Traceur [U/KgA]
  REAL,DIMENSION(klon,klev),INTENT(IN)          :: t_seri  ! Temperature
  REAL,DIMENSION(klon,klev),INTENT(IN)          :: pplay   ! pression pour le mileu de chaque couche (en Pa)
  REAL,DIMENSION(klon,klev+1),INTENT(IN)        :: paprs   ! pression pour chaque inter-couche (en Pa)
  LOGICAL,DIMENSION(klon,klev),INTENT(IN)       :: is_strato ! stratospheric flag

  ! local variables in coagulation routine
  INTEGER                                       :: i,j,k,nb,ilon,ilev

!--convert SO2 to H2SO4
  budg_3D_so2_to_h2so4(:,:)=0.0
  budg_so2_to_h2so4(:)=0.0

  DO ilon=1, klon
  DO ilev=1, klev
  !only in the stratosphere
  IF (is_strato(ilon,ilev)) THEN
    IF (SO2_lifetime(ilon,ilev).GT.0.0) THEN
      budg_3D_so2_to_h2so4(ilon,ilev)=tr_seri(ilon,ilev,id_SO2_strat)*(1.0-exp(-pdtphys/SO2_lifetime(ilon,ilev)))
    ENDIF
    tr_seri(ilon,ilev,id_SO2_strat)=tr_seri(ilon,ilev,id_SO2_strat) - budg_3D_so2_to_h2so4(ilon,ilev)
    tr_seri(ilon,ilev,id_H2SO4_strat)=tr_seri(ilon,ilev,id_H2SO4_strat) + mH2SO4mol/mSO2mol*budg_3D_so2_to_h2so4(ilon,ilev)
    !convert budget from kg(SO2)/kgA to kg(S)/m2/layer/s for saving as diagnostic
    budg_3D_so2_to_h2so4(ilon,ilev)=budg_3D_so2_to_h2so4(ilon,ilev)*mSatom/mSO2mol*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG/pdtphys
    budg_so2_to_h2so4(ilon)=budg_so2_to_h2so4(ilon)+budg_3D_so2_to_h2so4(ilon,ilev)
  ENDIF
  ENDDO
  ENDDO

END SUBROUTINE SO2_TO_H2SO4
