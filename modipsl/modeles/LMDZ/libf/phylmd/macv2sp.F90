SUBROUTINE MACv2SP(pphis,pplay,paprs,xlon,xlat,tau_allaer,piz_allaer,cg_allaer)
  !
  !--routine to read the MACv2SP plume and compute optical properties
  !--requires flag_aerosol = 7
  !--feeds into aerosol optical properties and newmicro cloud droplet size if ok_cdnc activated
  !--for this one needs to feed natural (pre-industrial) aerosols twice for nat and 1980 files
  !--pre-ind aerosols (index=1) are not changed, present-day aerosols (index=2) are incremented
  !--uses model year so year_cur needs to be correct in the model simulation
  !
  !--aod_prof = AOD per layer 
  !--ssa_prof = SSA 
  !--asy_prof = asymetry parameter
  !--dNovrN   = enhancement factor for CDNC
  !
  USE mo_simple_plumes, ONLY: sp_aop_profile
  USE phys_cal_mod, ONLY : year_cur, day_cur, year_len
  USE dimphy
  USE aero_mod
  USE phys_local_var_mod, ONLY: t_seri, od443aer, od550aer, od865aer, ec550aer, dryod550aer, od550lt1aer, dNovrN
  !!USE YOMCST, ONLY : RD, RG
  !
  IMPLICIT NONE
  !
  include "YOMCST.h"
  ! 
  REAL,DIMENSION(klon),INTENT(IN)        :: pphis   ! Geopotentiel de surface
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
  REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs   ! pression pour les interfaces de chaque couche (en Pa)
  REAL,DIMENSION(klon),INTENT(IN)        :: xlat    ! latitudes pour chaque point
  REAL,DIMENSION(klon),INTENT(IN)        :: xlon    ! longitudes pour chaque point
  !
  REAL, DIMENSION(klon,klev,2,nbands_sw_rrtm), INTENT(OUT) :: tau_allaer !  epaisseur optique aerosol
  REAL, DIMENSION(klon,klev,2,nbands_sw_rrtm), INTENT(OUT) :: piz_allaer !  single scattering albedo aerosol
  REAL, DIMENSION(klon,klev,2,nbands_sw_rrtm), INTENT(OUT) :: cg_allaer  !  asymmetry parameter aerosol
  !
  REAL,DIMENSION(klon,klev) :: aod_prof, ssa_prof, asy_prof
  REAL,DIMENSION(klon,klev) :: z, dz
  REAL,DIMENSION(klon)      :: oro, zrho, zt
  !
  INTEGER, PARAMETER :: nmon = 12
  !
  REAL, PARAMETER    :: l443 = 443.0, l550 = 550.0, l865 = 865.0 !--wavelengths in nm
  !
  INTEGER, PARAMETER :: Nwvmax=25
  REAL, DIMENSION(0:Nwvmax), PARAMETER :: lambda=(/ 240.0, &  !--this one is for band 1
                  280.0,  300.0,  330.0,  360.0,  400.0,   &  !--these are bounds of Streamer bands
                  440.0,  480.0,  520.0,  570.0,  640.0,   &
                  690.0,  750.0,  780.0,  870.0, 1000.0,   &
                 1100.0, 1190.0, 1280.0, 1530.0, 1640.0,   &
                 2130.0, 2380.0, 2910.0, 3420.0, 4000.0   /)
  !
  REAL, DIMENSION(1:Nwvmax-1), PARAMETER :: weight =(/    &   !--and the weights to be given to the bands
                 0.01,  4.05,  9.51, 15.99, 26.07, 33.10, &   !--corresponding to a typical solar spectrum 
                33.07, 39.91, 52.67, 27.89, 43.60, 13.67, &
                42.22, 40.12, 32.70, 14.44, 19.48, 14.23, &
                13.43, 16.42,  8.33,  0.95,  0.65,  2.76  /)
  !
  REAL :: zlambda, zweight
  REAL :: year_fr
  !
  INTEGER band, i, k, Nwv
  !
  ! define the height and dheight arrays
  !
  oro(:)  = pphis(:)/RG                             ! surface height in m
  !
  DO k = 1, klev
    zrho(:) = pplay(:,k)/t_seri(:,k)/RD                         ! air density in kg/m3
    dz(:,k) = (paprs(:,k)-paprs(:,k+1))/zrho(:)/RG              ! layer thickness in m
    IF (k==1) THEN
       z(:,1) = oro(:) + (paprs(:,1)-pplay(:,1))/zrho(:)/RG     ! altitude middle of first layer in m
       zt(:)  = oro(:) + dz(:,1)                                ! altitude top of first layer in m
    ELSE
      z(:,k) = zt(:) + (paprs(:,k)-pplay(:,k))/zrho(:)/RG       ! altitude middle of layer k in m
      zt(:)  = zt(:) + dz(:,k)                                  ! altitude top of layer k in m
    ENDIF
  ENDDO
  !
  !--fractional year
  !
  year_fr = FLOAT(year_cur) + (FLOAT(day_cur)-0.5) / FLOAT(year_len)
  IF (year_fr.LT.1850.0.OR.year_fr.GE.2017.0) THEN
     CALL abort_physic ('macv2sp','year not supported by plume model',1)
  ENDIF
  !
  !--call to sp routine -- 443 nm
  !
  CALL sp_aop_profile                                    ( &
       klev     ,klon ,l443 ,oro    ,xlon     ,xlat      , &
       year_fr  ,z    ,dz   ,dNovrN ,aod_prof ,ssa_prof  , &
       asy_prof )
  !
  !--AOD calculations for diagnostics
  od443aer(:)= od443aer(:)+SUM(aod_prof(:,:),dim=2)
  !
  !--call to sp routine -- 550 nm
  !
  CALL sp_aop_profile                                    ( &
       klev     ,klon ,l550 ,oro    ,xlon     ,xlat      , &
       year_fr  ,z    ,dz   ,dNovrN ,aod_prof ,ssa_prof  , &
       asy_prof )
  !
  !--AOD calculations for diagnostics
  od550aer(:)=od550aer(:)+SUM(aod_prof(:,:),dim=2)
  !
  !--dry AOD calculation for diagnostics
  dryod550aer(:)=dryod550aer(:)+od550aer(:)
  !
  !--fine-mode AOD calculation for diagnostics
  od550lt1aer(:)=od550lt1aer(:)+od550aer(:)
  !
  !--extinction coefficient for diagnostic
  ec550aer(:,:)=ec550aer(:,:)+aod_prof(:,:)/dz(:,:)
  !
  !--call to sp routine -- 865 nm
  !
  CALL sp_aop_profile                                    ( &
       klev     ,klon ,l865 ,oro    ,xlon     ,xlat      , &
       year_fr  ,z    ,dz   ,dNovrN ,aod_prof ,ssa_prof  , &
       asy_prof )
  !
  !--AOD calculations for diagnostics
  od865aer(:)=od865aer(:)+SUM(aod_prof(:,:),dim=2)
  !
  !--re-weighting of piz and cg arrays before adding the anthropogenic aerosols
  !--index 2 = all natural + anthropogenic aerosols
  piz_allaer(:,:,2,:)=piz_allaer(:,:,2,:)*tau_allaer(:,:,2,:)
  cg_allaer(:,:,2,:) =cg_allaer(:,:,2,:)*piz_allaer(:,:,2,:)
  !
  !--now computing the same at many wavelengths to fill the model bands
  !
  DO Nwv=0,Nwvmax-1

    IF (Nwv.EQ.0) THEN          !--RRTM spectral band 1
      zlambda=lambda(Nwv)
      zweight=1.0
      band=1
    ELSEIF (Nwv.LE.5) THEN      !--RRTM spectral band 2
      zlambda=0.5*(lambda(Nwv)+lambda(Nwv+1))
      zweight=weight(Nwv)/SUM(weight(1:5))
      band=2
    ELSEIF (Nwv.LE.10) THEN     !--RRTM spectral band 3
      zlambda=0.5*(lambda(Nwv)+lambda(Nwv+1))
      zweight=weight(Nwv)/SUM(weight(6:10))
      band=3
    ELSEIF (Nwv.LE.16) THEN     !--RRTM spectral band 4
      zlambda=0.5*(lambda(Nwv)+lambda(Nwv+1))
      zweight=weight(Nwv)/SUM(weight(11:16))
      band=4
    ELSEIF (Nwv.LE.21) THEN     !--RRTM spectral band 5
      zlambda=0.5*(lambda(Nwv)+lambda(Nwv+1))
      zweight=weight(Nwv)/SUM(weight(17:21))
      band=5
    ELSE                        !--RRTM spectral band 6
      zlambda=0.5*(lambda(Nwv)+lambda(Nwv+1))
      zweight=weight(Nwv)/SUM(weight(22:Nwvmax-1))
      band=6
    ENDIF
    !
    CALL sp_aop_profile                                       ( &
         klev     ,klon ,zlambda ,oro    ,xlon     ,xlat      , &
         year_fr  ,z    ,dz      ,dNovrN ,aod_prof ,ssa_prof  , &
         asy_prof )
    !
    !--adding up the quantities tau, piz*tau and cg*piz*tau
    tau_allaer(:,:,2,band)=tau_allaer(:,:,2,band)+zweight*MAX(aod_prof(:,:),1.e-15)
    piz_allaer(:,:,2,band)=piz_allaer(:,:,2,band)+zweight*MAX(aod_prof(:,:),1.e-15)*ssa_prof(:,:)
    cg_allaer(:,:,2,band) =cg_allaer(:,:,2,band) +zweight*MAX(aod_prof(:,:),1.e-15)*ssa_prof(:,:)*asy_prof(:,:)
    !
  ENDDO
  !
  !--renpomalizing cg and piz now that MACv2SP increments have been added
  cg_allaer(:,:,2,:) =cg_allaer(:,:,2,:) /piz_allaer(:,:,2,:)
  piz_allaer(:,:,2,:)=piz_allaer(:,:,2,:)/tau_allaer(:,:,2,:)
  !
END SUBROUTINE MACv2SP
