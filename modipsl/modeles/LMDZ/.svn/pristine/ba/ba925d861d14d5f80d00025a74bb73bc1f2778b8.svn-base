MODULE init_ssrf_m
!
!*******************************************************************************

  USE indice_sol_mod, ONLY: is_ter, is_oce, is_oce, is_lic, epsfra
  USE dimphy,             ONLY: klon, zmasq
  USE phys_state_var_mod, ONLY: pctsrf
  USE geometry_mod, ONLY : longitude_deg, latitude_deg
  USE grid_atob_m,        ONLY: grille_m
  USE ioipsl,             ONLY: flininfo, flinopen, flinget, flinclo
  USE ioipsl_getin_p_mod, ONLY: getin_p
  USE comconst_mod, ONLY: im, pi

  CHARACTER(LEN=256), PARAMETER :: icefname="landiceref.nc", icevar="landice"
  PRIVATE
  PUBLIC :: start_init_subsurf
  include "iniprint.h"
  include "dimensions.h"
  include "paramet.h"
  include "comgeom2.h"

CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE start_init_subsurf(known_mask)
!
!-------------------------------------------------------------------------------
! Purpose: Subsurfaces initialization.
!-------------------------------------------------------------------------------
! Comment: Called by etat0phys_netcdf ; also called by limit_netcdf in case
!          no starting states are required (ok_etat0==.FALSE.).
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  LOGICAL, INTENT(IN) :: known_mask
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER           :: iml_lic, jml_lic
  INTEGER           :: fid, llm_tmp, ttm_tmp, itaul(1), ji, j
  REAL, ALLOCATABLE :: dlon_lic(:), lon_lic(:,:), fraclic (:,:)
  REAL, ALLOCATABLE :: dlat_lic(:), lat_lic(:,:), flic_tmp(:,:), vtmp(:,:)
  REAL              :: date, lev(1), dt, deg2rad
  LOGICAL           :: no_ter_antartique   ! If true, no land points are allowed at Antartic
!-------------------------------------------------------------------------------
  deg2rad= pi/180.0

!--- Physical grid points coordinates
  DO j=2,jjm; latitude_deg((j-2)*iim+2:(j-1)*iim+1)=rlatu(j);    END DO
  DO j=2,jjm; longitude_deg((j-2)*iim+2:(j-1)*iim+1)=rlonv(1:im); END DO
  latitude_deg(1) = pi/2.; latitude_deg(klon) = - pi/2.
  latitude_deg(:)=latitude_deg(:)/deg2rad
  longitude_deg(1) = 0.0;   longitude_deg(klon) = 0.0;
  longitude_deg(:)=longitude_deg(:)/deg2rad

! Compute ground geopotential, sub-cells quantities and possibly the mask.
! Sub-surfaces initialization
!*******************************************************************************
!--- Read and interpolate on model T-grid soil fraction and soil ice fraction.
  CALL flininfo(icefname, iml_lic, jml_lic, llm_tmp, ttm_tmp, fid)
  ALLOCATE(lat_lic(iml_lic,jml_lic),lon_lic(iml_lic,jml_lic))
  ALLOCATE(fraclic(iml_lic,jml_lic))
  CALL flinopen(icefname, .FALSE., iml_lic, jml_lic, llm_tmp,  &
 &               lon_lic, lat_lic, lev, ttm_tmp, itaul, date, dt, fid)
  CALL flinget(fid, icevar, iml_lic, jml_lic, llm_tmp, ttm_tmp, 1,1, fraclic)
  CALL flinclo(fid)
  WRITE(lunout,*)'landice dimensions: iml_lic, jml_lic : ',iml_lic,jml_lic

  ALLOCATE(dlon_lic(iml_lic),dlat_lic(jml_lic))
  dlon_lic(:)=lon_lic(:,1); IF(MAXVAL(dlon_lic)>pi) dlon_lic=dlon_lic*pi/180.
  dlat_lic(:)=lat_lic(1,:); IF(MAXVAL(dlat_lic)>pi) dlat_lic=dlat_lic*pi/180.
  DEALLOCATE(lon_lic,lat_lic); ALLOCATE(flic_tmp(iip1,jjp1))
  CALL grille_m(dlon_lic,dlat_lic,fraclic,rlonv(1:iim),rlatu,flic_tmp(1:iim,:))
  flic_tmp(iip1,:)=flic_tmp(1,:)

!--- To the physical grid
  pctsrf(:,:) = 0.
  CALL gr_dyn_fi(1, iip1, jjp1, klon, flic_tmp, pctsrf(:,is_lic))
  DEALLOCATE(flic_tmp)

!--- Adequation with soil/sea mask
  WHERE(pctsrf(:,is_lic)<EPSFRA) pctsrf(:,is_lic)=0. 
  WHERE(zmasq(:)<EPSFRA)         pctsrf(:,is_lic)=0.
  pctsrf(:,is_ter)=zmasq(:)
  DO ji=1,klon
    IF(zmasq(ji)>EPSFRA) THEN 
      IF(pctsrf(ji,is_lic)>=zmasq(ji)) THEN
        pctsrf(ji,is_lic)=zmasq(ji)
        pctsrf(ji,is_ter)=0.
      ELSE
        pctsrf(ji,is_ter)=zmasq(ji)-pctsrf(ji,is_lic)
        IF(pctsrf(ji,is_ter)<EPSFRA) THEN
          pctsrf(ji,is_ter)=0.
          pctsrf(ji,is_lic)=zmasq(ji)
        END IF 
      END IF 
    END IF 
  END DO 


  !--- Option no_ter_antartique removes all land fractions souther than 60S.
  !--- Land ice is set instead of the land fractions on these latitudes.
  !--- The ocean and sea-ice fractions are not changed.
  no_ter_antartique=.FALSE.
  CALL getin_p('no_ter_antartique',no_ter_antartique)
  WRITE(lunout,*)"no_ter_antartique=",no_ter_antartique
  IF (no_ter_antartique) THEN
     ! Remove all land fractions souther than 60S and set land-ice instead
     WRITE(lunout,*) "Remove land fractions souther than 60deg south by increasing"
     WRITE(lunout,*) "the continental ice fractions. No land can now be found at Antartic."
     DO ji=1, klon
        IF (latitude_deg(ji)<-60.0) THEN
           pctsrf(ji,is_lic) = pctsrf(ji,is_lic) + pctsrf(ji,is_ter)
           pctsrf(ji,is_ter) = 0
        END IF
     END DO
  END IF


!--- Sub-surface ocean and sea ice (sea ice set to zero for start).
  pctsrf(:,is_oce)=(1.-zmasq(:))
  WHERE(pctsrf(:,is_oce)<EPSFRA) pctsrf(:,is_oce)=0.
  IF(known_mask) pctsrf(:,is_oce)=1-zmasq(:)

!--- It is checked that the sub-surfaces sum is equal to 1.
  ji=COUNT((ABS(SUM(pctsrf(:,:),dim=2))-1.0)>EPSFRA)
  IF(ji/=0) WRITE(lunout,*) 'Sub-cell distribution problem for ',ji,' points'

END SUBROUTINE start_init_subsurf
!
!-------------------------------------------------------------------------------

END MODULE init_ssrf_m
!
!*******************************************************************************
