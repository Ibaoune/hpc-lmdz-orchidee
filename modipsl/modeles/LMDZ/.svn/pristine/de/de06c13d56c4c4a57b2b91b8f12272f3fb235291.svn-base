!Completed
MODULE ocean_slab_mod
!
! This module is used for both surface ocean and sea-ice when using the slab ocean,
! "ocean=slab".
!

  USE dimphy
  USE indice_sol_mod
  USE surface_data
  USE mod_grid_phy_lmdz, ONLY: klon_glo
  USE mod_phys_lmdz_mpi_data, ONLY: is_mpi_root

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ocean_slab_init, ocean_slab_frac, ocean_slab_noice, ocean_slab_ice

!***********************************************************************************
! Global saved variables
!***********************************************************************************
  ! number of slab vertical layers
  INTEGER, PUBLIC, SAVE :: nslay
  !$OMP THREADPRIVATE(nslay)
  ! timestep for coupling (update slab temperature) in timesteps
  INTEGER, PRIVATE, SAVE                           :: cpl_pas
  !$OMP THREADPRIVATE(cpl_pas)
  ! cyang = 1/heat capacity of top layer (rho.c.H)
  REAL, PRIVATE, SAVE                              :: cyang
  !$OMP THREADPRIVATE(cyang)
  ! depth of slab layers (1 or 2)
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE   :: slabh
  !$OMP THREADPRIVATE(slabh)
  ! slab temperature
  REAL, ALLOCATABLE, DIMENSION(:,:), PUBLIC, SAVE   :: tslab
  !$OMP THREADPRIVATE(tslab)
  ! heat flux convergence due to Ekman
  REAL, ALLOCATABLE, DIMENSION(:,:), PUBLIC, SAVE   :: dt_ekman
  !$OMP THREADPRIVATE(dt_ekman)
  ! heat flux convergence due to horiz diffusion
  REAL, ALLOCATABLE, DIMENSION(:,:), PUBLIC, SAVE   :: dt_hdiff
  !$OMP THREADPRIVATE(dt_hdiff)
  ! heat flux convergence due to GM eddy advection
  REAL, ALLOCATABLE, DIMENSION(:,:), PUBLIC, SAVE   :: dt_gm
  !$OMP THREADPRIVATE(dt_gm)
  ! Heat Flux correction 
  REAL, ALLOCATABLE, DIMENSION(:,:), PUBLIC, SAVE   :: dt_qflux
  !$OMP THREADPRIVATE(dt_qflux)
  ! fraction of ocean covered by sea ice (sic / (oce+sic))
  REAL, ALLOCATABLE, DIMENSION(:), PUBLIC, SAVE  :: fsic
  !$OMP THREADPRIVATE(fsic)
  ! temperature of the sea ice
  REAL, ALLOCATABLE, DIMENSION(:), PUBLIC, SAVE  :: tice
  !$OMP THREADPRIVATE(tice)
  ! sea ice thickness, in kg/m2
  REAL, ALLOCATABLE, DIMENSION(:), PUBLIC, SAVE  :: seaice
  !$OMP THREADPRIVATE(seaice)
  ! net surface heat flux, weighted by open ocean fraction
  ! slab_bils accumulated over cpl_pas timesteps
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE  :: bils_cum
  !$OMP THREADPRIVATE(bils_cum)
  ! net heat flux into the ocean below the ice : conduction + solar radiation
  REAL, ALLOCATABLE, DIMENSION(:), PUBLIC, SAVE  :: slab_bilg
  !$OMP THREADPRIVATE(slab_bilg)
  ! slab_bilg over cpl_pas timesteps
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE  :: bilg_cum
  !$OMP THREADPRIVATE(bilg_cum)
  ! wind stress saved over cpl_pas timesteps
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE  :: taux_cum
  !$OMP THREADPRIVATE(taux_cum)
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE  :: tauy_cum
  !$OMP THREADPRIVATE(tauy_cum)

!***********************************************************************************
! Parameters (could be read in def file: move to slab_init)
!***********************************************************************************
! snow and ice physical characteristics:
    REAL, PARAMETER :: t_freeze=271.35 ! freezing sea water temp
    REAL, PARAMETER :: t_melt=273.15   ! melting ice temp
    REAL, PARAMETER :: sno_den=300. !mean snow density, kg/m3
    REAL, PARAMETER :: ice_den=917. ! ice density
    REAL, PARAMETER :: sea_den=1025. ! sea water density
    REAL, PARAMETER :: ice_cond=2.17*ice_den !conductivity of ice
    REAL, PARAMETER :: sno_cond=0.31*sno_den ! conductivity of snow
    REAL, PARAMETER :: ice_cap=2067.   ! specific heat capacity, snow and ice
    REAL, PARAMETER :: sea_cap=3995.   ! specific heat capacity, snow and ice
    REAL, PARAMETER :: ice_lat=334000. ! freeze /melt latent heat snow and ice

! control of snow and ice cover & freeze / melt (heights converted to kg/m2)
    REAL, PARAMETER :: snow_min=0.05*sno_den !critical snow height 5 cm
    REAL, PARAMETER :: snow_wfact=0.4 ! max fraction of falling snow blown into ocean
    REAL, PARAMETER :: ice_frac_min=0.001
    REAL, PARAMETER :: ice_frac_max=1. ! less than 1. if min leads fraction 
    REAL, PARAMETER :: h_ice_min=0.01*ice_den ! min ice thickness 
    REAL, PARAMETER :: h_ice_thin=0.15*ice_den ! thin ice thickness 
    ! below ice_thin, priority is melt lateral / grow height
    ! ice_thin is also height of new ice
    REAL, PARAMETER :: h_ice_thick=2.5*ice_den ! thin ice thickness 
    ! above ice_thick, priority is melt height / grow lateral
    REAL, PARAMETER :: h_ice_new=1.*ice_den ! max height of new open ocean ice 
    REAL, PARAMETER :: h_ice_max=10.*ice_den ! max ice height 

! albedo  and radiation parameters
    REAL, PARAMETER :: alb_sno_min=0.55 !min snow albedo
    REAL, PARAMETER :: alb_sno_del=0.3  !max snow albedo = min + del
    REAL, PARAMETER :: alb_ice_dry=0.75 !dry thick ice
    REAL, PARAMETER :: alb_ice_wet=0.66 !melting thick ice
    REAL, PARAMETER :: pen_frac=0.3 !fraction of shortwave penetrating into the
    ! ice (no snow)
    REAL, PARAMETER :: pen_ext=1.5 !extinction of penetrating shortwave (m-1)

! horizontal transport
   LOGICAL, PUBLIC, SAVE :: slab_hdiff
   !$OMP THREADPRIVATE(slab_hdiff)
   LOGICAL, PUBLIC, SAVE :: slab_gm
   !$OMP THREADPRIVATE(slab_gm)
   REAL, PRIVATE, SAVE    :: coef_hdiff ! coefficient for horizontal diffusion
   !$OMP THREADPRIVATE(coef_hdiff)
   INTEGER, PUBLIC, SAVE :: slab_ekman, slab_cadj ! Ekman, conv adjustment
   !$OMP THREADPRIVATE(slab_ekman)
   !$OMP THREADPRIVATE(slab_cadj)

!***********************************************************************************

CONTAINS
!
!***********************************************************************************
!
  SUBROUTINE ocean_slab_init(dtime, pctsrf_rst)
  !, seaice_rst etc

    USE ioipsl_getin_p_mod, ONLY : getin_p
    USE mod_phys_lmdz_transfert_para, ONLY : gather
    USE slab_heat_transp_mod, ONLY : ini_slab_transp

    ! Input variables
!***********************************************************************************
    REAL, INTENT(IN)                         :: dtime
! Variables read from restart file
    REAL, DIMENSION(klon, nbsrf), INTENT(IN) :: pctsrf_rst 
    ! surface fractions from start file

! Local variables
!************************************************************************************
    INTEGER                :: error
    REAL, DIMENSION(klon_glo) :: zmasq_glo
    CHARACTER (len = 80)   :: abort_message
    CHARACTER (len = 20)   :: modname = 'ocean_slab_intit'

!***********************************************************************************
! Define some parameters
!***********************************************************************************
! Number of slab layers
    nslay=2
    CALL getin_p('slab_layers',nslay)
    print *,'number of slab layers : ',nslay
! Layer thickness
    ALLOCATE(slabh(nslay), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation slabh'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    slabh(1)=50.
    CALL getin_p('slab_depth',slabh(1))
    IF (nslay.GT.1) THEN 
        slabh(2)=150.
    END IF

! cyang = 1/heat capacity of top layer (rho.c.H)
    cyang=1/(slabh(1)*sea_den*sea_cap)

! cpl_pas  coupling period (update of tslab and ice fraction)
! pour un calcul a chaque pas de temps, cpl_pas=1
    cpl_pas = NINT(86400./dtime * 1.0) ! une fois par jour
    CALL getin_p('cpl_pas',cpl_pas)
    print *,'cpl_pas',cpl_pas

! Horizontal diffusion
    slab_hdiff=.FALSE.
    CALL getin_p('slab_hdiff',slab_hdiff)
    coef_hdiff=25000.
    CALL getin_p('coef_hdiff',coef_hdiff)
! Ekman transport
    slab_ekman=0
    CALL getin_p('slab_ekman',slab_ekman)
! GM eddy advection (2-layers only)
    slab_gm=.FALSE.
    CALL getin_p('slab_gm',slab_gm)
    IF (slab_ekman.LT.2) THEN
       slab_gm=.FALSE.
    ENDIF 
! Convective adjustment
    IF (nslay.EQ.1) THEN
        slab_cadj=0
    ELSE
        slab_cadj=1
    END IF
    CALL getin_p('slab_cadj',slab_cadj)

!************************************************************************************
! Allocate surface fraction read from restart file
!************************************************************************************
    ALLOCATE(fsic(klon), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation tmp_pctsrf_slab'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    fsic(:)=0.
    !zmasq = continent fraction
    WHERE (1.-zmasq(:)>EPSFRA)
        fsic(:) = pctsrf_rst(:,is_sic)/(1.-zmasq(:))
    END WHERE

!************************************************************************************
! Allocate saved fields
!************************************************************************************
    ALLOCATE(tslab(klon,nslay), stat=error)
       IF (error /= 0) CALL abort_physic &
         (modname,'pb allocation tslab', 1)

    ALLOCATE(bils_cum(klon), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation slab_bils_cum'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    bils_cum(:) = 0.0   

    IF (version_ocean=='sicINT') THEN ! interactive sea ice
        ALLOCATE(slab_bilg(klon), stat = error)
        IF (error /= 0) THEN
           abort_message='Pb allocation slab_bilg'
           CALL abort_physic(modname,abort_message,1)
        ENDIF
        slab_bilg(:) = 0.0   
        ALLOCATE(bilg_cum(klon), stat = error)
        IF (error /= 0) THEN
           abort_message='Pb allocation slab_bilg_cum'
           CALL abort_physic(modname,abort_message,1)
        ENDIF
        bilg_cum(:) = 0.0   
        ALLOCATE(tice(klon), stat = error)
        IF (error /= 0) THEN
           abort_message='Pb allocation slab_tice'
           CALL abort_physic(modname,abort_message,1)
        ENDIF
        ALLOCATE(seaice(klon), stat = error)
        IF (error /= 0) THEN
           abort_message='Pb allocation slab_seaice'
           CALL abort_physic(modname,abort_message,1)
        ENDIF
    END IF

    IF (slab_hdiff) THEN !horizontal diffusion
        ALLOCATE(dt_hdiff(klon,nslay), stat = error)
        IF (error /= 0) THEN
           abort_message='Pb allocation dt_hdiff'
           CALL abort_physic(modname,abort_message,1)
        ENDIF
        dt_hdiff(:,:) = 0.0   
    ENDIF

    ALLOCATE(dt_qflux(klon,nslay), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation dt_qflux'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    dt_qflux(:,:) = 0.0   

    IF (slab_gm) THEN !GM advection
        ALLOCATE(dt_gm(klon,nslay), stat = error)
        IF (error /= 0) THEN
           abort_message='Pb allocation dt_gm'
           CALL abort_physic(modname,abort_message,1)
        ENDIF
        dt_gm(:,:) = 0.0   
    ENDIF

    IF (slab_ekman.GT.0) THEN ! ekman transport
        ALLOCATE(dt_ekman(klon,nslay), stat = error)
        IF (error /= 0) THEN
           abort_message='Pb allocation dt_ekman'
           CALL abort_physic(modname,abort_message,1)
        ENDIF
        dt_ekman(:,:) = 0.0   
        ALLOCATE(taux_cum(klon), stat = error)
        IF (error /= 0) THEN
           abort_message='Pb allocation taux_cum'
           CALL abort_physic(modname,abort_message,1)
        ENDIF
        taux_cum(:) = 0.0   
        ALLOCATE(tauy_cum(klon), stat = error)
        IF (error /= 0) THEN
           abort_message='Pb allocation tauy_cum'
           CALL abort_physic(modname,abort_message,1)
        ENDIF
        tauy_cum(:) = 0.0   
    ENDIF

! Initialize transport
    IF (slab_hdiff.OR.(slab_ekman.GT.0)) THEN
      CALL gather(zmasq,zmasq_glo)
! Master thread/process only
!$OMP MASTER  
      IF (is_mpi_root) THEN 
          CALL ini_slab_transp(zmasq_glo)
      END IF
!$OMP END MASTER
    END IF

 END SUBROUTINE ocean_slab_init
!
!***********************************************************************************
!
  SUBROUTINE ocean_slab_frac(itime, dtime, jour, pctsrf_chg, is_modified)

! this routine sends back the sea ice and ocean fraction to the main physics
! routine. Called only with interactive sea ice

! Arguments
!************************************************************************************
    INTEGER, INTENT(IN)                        :: itime   ! current timestep 
    INTEGER, INTENT(IN)                        :: jour    !  day in year (not 
    REAL   , INTENT(IN)                        :: dtime   ! physics timestep (s)
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT) :: pctsrf_chg  ! sub-surface fraction
    LOGICAL, INTENT(OUT)                       :: is_modified ! true if pctsrf is 
                                                         ! modified at this time step

       pctsrf_chg(:,is_oce)=(1.-fsic(:))*(1.-zmasq(:))
       pctsrf_chg(:,is_sic)=fsic(:)*(1.-zmasq(:))
       is_modified=.TRUE.

  END SUBROUTINE ocean_slab_frac
!
!************************************************************************************
!
  SUBROUTINE ocean_slab_noice( & 
       itime, dtime, jour, knon, knindex, &
       p1lay, cdragh, cdragq, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, gustiness, tsurf_in, &
       radsol, snow, &
       qsurf, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
       tsurf_new, dflux_s, dflux_l, slab_bils)
    
    USE calcul_fluxs_mod
    USE slab_heat_transp_mod, ONLY: divgrad_phy,slab_ekman1,slab_ekman2,slab_gmdiff
    USE mod_phys_lmdz_para

    INCLUDE "clesphys.h"

! This routine 
! (1) computes surface turbulent fluxes over points with some open ocean    
! (2) reads additional Q-flux (everywhere)
! (3) computes horizontal transport (diffusion & Ekman)
! (4) updates slab temperature every cpl_pas ; creates new ice if needed.

! Note : 
! klon total number of points
! knon number of points with open ocean (varies with time)
! knindex gives position of the knon points within klon.
! In general, local saved variables have klon values
! variables exchanged with PBL module have knon.

! Input arguments
!***********************************************************************************
    INTEGER, INTENT(IN)                  :: itime ! current timestep INTEGER,
    INTEGER, INTENT(IN)                  :: jour  ! day in year (for Q-Flux) 
    INTEGER, INTENT(IN)                  :: knon  ! number of points 
    INTEGER, DIMENSION(klon), INTENT(IN) :: knindex 
    REAL, INTENT(IN) :: dtime  ! timestep (s)
    REAL, DIMENSION(klon), INTENT(IN)    :: p1lay 
    REAL, DIMENSION(klon), INTENT(IN)    :: cdragh, cdragq, cdragm 
    ! drag coefficients
    REAL, DIMENSION(klon), INTENT(IN)    :: precip_rain, precip_snow 
    REAL, DIMENSION(klon), INTENT(IN)    :: temp_air, spechum ! near surface T, q
    REAL, DIMENSION(klon), INTENT(IN)    :: AcoefH, AcoefQ, BcoefH, BcoefQ 
    REAL, DIMENSION(klon), INTENT(IN)    :: AcoefU, AcoefV, BcoefU, BcoefV 
    ! exchange coefficients for boundary layer scheme
    REAL, DIMENSION(klon), INTENT(IN)    :: ps  ! surface pressure
    REAL, DIMENSION(klon), INTENT(IN)    :: u1, v1, gustiness ! surface wind
    REAL, DIMENSION(klon), INTENT(IN)    :: tsurf_in ! surface temperature
    REAL, DIMENSION(klon), INTENT(INOUT) :: radsol ! net surface radiative flux

! In/Output arguments
!************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT) :: snow ! in kg/m2
    
! Output arguments
!************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)   :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)   :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)   :: flux_u1, flux_v1
    REAL, DIMENSION(klon), INTENT(OUT)   :: tsurf_new ! new surface tempearture
    REAL, DIMENSION(klon), INTENT(OUT)   :: dflux_s, dflux_l      
    REAL, DIMENSION(klon), INTENT(OUT)   :: slab_bils

! Local variables
!************************************************************************************
    INTEGER               :: i,ki,k
    REAL                  :: t_cadj
    !  for surface heat fluxes
    REAL, DIMENSION(klon) :: cal, beta, dif_grnd
    ! for Q-Flux computation: d/dt SST, d/dt ice volume (kg/m2), surf fluxes 
    REAL, DIMENSION(klon) :: diff_sst, diff_siv
    REAL, DIMENSION(klon,nslay) :: lmt_bils
    ! for surface wind stress
    REAL, DIMENSION(klon) :: u0, v0
    REAL, DIMENSION(klon) :: u1_lay, v1_lay
    ! for new ice creation
    REAL                  :: e_freeze, h_new, dfsic
    ! horizontal diffusion and Ekman local vars
    ! dimension = global domain (klon_glo) instead of // subdomains 
    REAL, DIMENSION(klon_glo,nslay) :: dt_hdiff_glo,dt_ekman_glo,dt_gm_glo
    ! dt_ekman_glo saved for diagnostic, dt_ekman_tmp used for time loop
    REAL, DIMENSION(klon_glo,nslay) :: dt_hdiff_tmp, dt_ekman_tmp
    REAL, DIMENSION(klon_glo,nslay) :: tslab_glo
    REAL, DIMENSION(klon_glo) :: taux_glo,tauy_glo

!****************************************************************************************
! 1) Surface fluxes calculation
!    
!****************************************************************************************
    cal(:)      = 0. ! infinite thermal inertia
    beta(:)     = 1. ! wet surface
    dif_grnd(:) = 0. ! no diffusion into ground
    
! Suppose zero surface speed
    u0(:)=0.0
    v0(:)=0.0
    u1_lay(:) = u1(:) - u0(:)
    v1_lay(:) = v1(:) - v0(:)

! Compute latent & sensible fluxes
    CALL calcul_fluxs(knon, is_oce, dtime, &
         tsurf_in, p1lay, cal, beta, cdragh, cdragq, ps, &
         precip_rain, precip_snow, snow, qsurf,  &
         radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, gustiness, &
         f_qsat_oce,AcoefH, AcoefQ, BcoefH, BcoefQ, &
         tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

! save total cumulated heat fluxes locally
! radiative + turbulent + melt of falling snow 
    slab_bils(:)=0.
    DO i=1,knon
        ki=knindex(i)
        slab_bils(ki)=(1.-fsic(ki))*(fluxlat(i)+fluxsens(i)+radsol(i) &
                      -precip_snow(i)*ice_lat*(1.+snow_wfact*fsic(ki)))
        bils_cum(ki)=bils_cum(ki)+slab_bils(ki)
    END DO

!  Compute surface wind stress
    CALL calcul_flux_wind(knon, dtime, &
         u0, v0, u1, v1, gustiness, cdragm, &
         AcoefU, AcoefV, BcoefU, BcoefV, &
         p1lay, temp_air, &
         flux_u1, flux_v1)  

! save cumulated wind stress
    IF (slab_ekman.GT.0) THEN 
      DO i=1,knon
          ki=knindex(i)
          taux_cum(ki)=taux_cum(ki)+flux_u1(i)*(1.-fsic(ki))/cpl_pas
          tauy_cum(ki)=tauy_cum(ki)+flux_v1(i)*(1.-fsic(ki))/cpl_pas
      END DO
    ENDIF

!****************************************************************************************
! 2) Q-Flux : get global variables lmt_bils, diff_sst and diff_siv from file limit_slab.nc
!
!****************************************************************************************
    CALL limit_slab(itime, dtime, jour, lmt_bils, diff_sst, diff_siv) 
    ! lmt_bils and diff_sst,siv saved by limit_slab
    ! qflux = total QFlux correction (in W/m2)
    dt_qflux(:,1)=lmt_bils(:,1)+diff_sst(:)/cyang/86400.-diff_siv(:)*ice_den*ice_lat/86400.
    IF (nslay.GT.1) THEN
       dt_qflux(:,2:nslay)=lmt_bils(:,2:nslay)
    END IF

!****************************************************************************************
! 3) Recalculate new temperature (add Surf fluxes, Q-Flux, Ocean transport)
!    Bring to freezing temp and make sea ice if necessary
!  
!***********************************************o*****************************************
    tsurf_new=tsurf_in
    IF (MOD(itime,cpl_pas).EQ.0) THEN ! time to update tslab & fraction
! ***********************************
!  Horizontal transport 
! ***********************************
      IF (slab_ekman.GT.0) THEN 
          ! copy wind stress to global var
          CALL gather(taux_cum,taux_glo)
          CALL gather(tauy_cum,tauy_glo)
      END IF

      IF (slab_hdiff.OR.(slab_ekman.GT.0)) THEN 
        CALL gather(tslab,tslab_glo)
      ! Compute horiz transport on one process only
        IF (is_mpi_root .AND. is_omp_root) THEN ! Only master processus          
          IF (slab_hdiff) THEN 
              dt_hdiff_glo(:,:)=0.
          END IF
          IF (slab_ekman.GT.0) THEN 
              dt_ekman_glo(:,:)=0.
          END IF
          IF (slab_gm) THEN 
              dt_gm_glo(:,:)=0.
          END IF
          DO i=1,cpl_pas ! time splitting for numerical stability
            IF (slab_ekman.GT.0) THEN
              SELECT CASE (slab_ekman)
                CASE (1)
                  CALL slab_ekman1(taux_glo,tauy_glo,tslab_glo,dt_ekman_tmp)
                CASE (2)
                  CALL slab_ekman2(taux_glo,tauy_glo,tslab_glo,dt_ekman_tmp,dt_hdiff_tmp,slab_gm)
                CASE DEFAULT
                  dt_ekman_tmp(:,:)=0.
              END SELECT
              dt_ekman_glo(:,:)=dt_ekman_glo(:,:)+dt_ekman_tmp(:,:)
              ! convert dt_ekman from K.s-1.(kg.m-2) to K.s-1   
              DO k=1,nslay
                dt_ekman_tmp(:,k)=dt_ekman_tmp(:,k)/(slabh(k)*sea_den)
              ENDDO
              tslab_glo=tslab_glo+dt_ekman_tmp*dtime
              IF (slab_gm) THEN ! Gent-McWilliams eddy advection
                dt_gm_glo(:,:)=dt_gm_glo(:,:)+ dt_hdiff_tmp(:,:)
                ! convert dt from K.s-1.(kg.m-2) to K.s-1   
                DO k=1,nslay
                  dt_hdiff_tmp(:,k)=dt_hdiff_tmp(:,k)/(slabh(k)*sea_den) 
                END DO
                tslab_glo=tslab_glo+dt_hdiff_tmp*dtime
              END IF
            ENDIF
! GM included in Ekman_2
!            IF (slab_gm) THEN ! Gent-McWilliams eddy advection
!              CALL slab_gmdiff(tslab_glo,dt_hdiff_tmp)
!              ! convert dt_gm from K.m.s-1 to K.s-1   
!              DO k=1,nslay
!                dt_hdiff_tmp(:,k)=dt_hdiff_tmp(:,k)/slabh(k) 
!              END DO
!              tslab_glo=tslab_glo+dt_hdiff_tmp*dtime
!              dt_gm_glo(:,:)=dt_gm_glo(:,:)+ dt_hdiff_tmp(:,:)
!            END IF
            IF (slab_hdiff) THEN ! horizontal diffusion
              ! laplacian of slab T
              CALL divgrad_phy(nslay,tslab_glo,dt_hdiff_tmp)
              ! multiply by diff coef and normalize to 50m slab equivalent
              dt_hdiff_tmp=dt_hdiff_tmp*coef_hdiff*50./SUM(slabh)
              dt_hdiff_glo(:,:)=dt_hdiff_glo(:,:)+ dt_hdiff_tmp(:,:)
              tslab_glo=tslab_glo+dt_hdiff_tmp*dtime
            END IF
          END DO ! time splitting
          IF (slab_hdiff) THEN
            !dt_hdiff_glo saved in W/m2
            DO k=1,nslay
              dt_hdiff_glo(:,k)=dt_hdiff_glo(:,k)*slabh(k)*sea_den*sea_cap/cpl_pas
            END DO
          END IF
          IF (slab_gm) THEN
            !dt_hdiff_glo saved in W/m2
            dt_gm_glo(:,:)=dt_gm_glo(:,:)*sea_cap/cpl_pas
          END IF
          IF (slab_ekman.GT.0) THEN 
            ! dt_ekman_glo saved in W/m2
            dt_ekman_glo(:,:)=dt_ekman_glo(:,:)*sea_cap/cpl_pas
          END IF
        END IF ! master process
!$OMP BARRIER
        ! Send new fields back to all processes
        CALL Scatter(tslab_glo,tslab)
        IF (slab_hdiff) THEN
          CALL Scatter(dt_hdiff_glo,dt_hdiff)
        END IF
        IF (slab_gm) THEN
          CALL Scatter(dt_gm_glo,dt_gm)
        END IF
        IF (slab_ekman.GT.0) THEN 
          CALL Scatter(dt_ekman_glo,dt_ekman)
          ! clear wind stress
          taux_cum(:)=0.
          tauy_cum(:)=0.
        END IF
      ENDIF ! transport

! ***********************************
! Other heat fluxes
! ***********************************
      ! Add read QFlux
      DO k=1,nslay
        tslab(:,k)=tslab(:,k)+dt_qflux(:,k)*cyang*dtime*cpl_pas &
                   *slabh(1)/slabh(k) 
      END DO
      ! Add cumulated surface fluxes
      tslab(:,1)=tslab(:,1)+bils_cum(:)*cyang*dtime
      ! Convective adjustment if 2 layers 
      IF ((nslay.GT.1).AND.(slab_cadj.GT.0)) THEN
        DO i=1,klon
          IF (tslab(i,2).GT.tslab(i,1)) THEN
            ! mean (mass-weighted) temperature
            t_cadj=SUM(tslab(i,:)*slabh(:))/SUM(slabh(:))
            tslab(i,1)=t_cadj
            tslab(i,2)=t_cadj
          END IF
        END DO
      END IF
! ***********************************
! Update surface temperature and ice 
! ***********************************
      SELECT CASE(version_ocean)
      CASE('sicNO') ! no sea ice even below freezing !
          DO i=1,knon
              ki=knindex(i)
              tsurf_new(i)=tslab(ki,1)
          END DO
      CASE('sicOBS') ! "realistic" case, for prescribed sea ice
        ! tslab cannot be below freezing, or above it if there is sea ice
          DO i=1,knon
              ki=knindex(i)
              IF ((tslab(ki,1).LT.t_freeze).OR.(fsic(ki).GT.epsfra)) THEN
                  tslab(ki,1)=t_freeze
              END IF
              tsurf_new(i)=tslab(ki,1)
          END DO
      CASE('sicINT') ! interactive sea ice
          DO i=1,knon
              ki=knindex(i)
              IF (fsic(ki).LT.epsfra) THEN ! Free of ice
                  IF (tslab(ki,1).LT.t_freeze) THEN ! create new ice
                      ! quantity of new ice formed
                      e_freeze=(t_freeze-tslab(ki,1))/cyang/ice_lat
                      ! new ice
                      tice(ki)=t_freeze
                      fsic(ki)=MIN(ice_frac_max,e_freeze/h_ice_thin)
                      IF (fsic(ki).GT.ice_frac_min) THEN
                          seaice(ki)=MIN(e_freeze/fsic(ki),h_ice_max)
                          tslab(ki,1)=t_freeze
                      ELSE
                          fsic(ki)=0.
                      END IF
                      tsurf_new(i)=t_freeze
                  ELSE
                      tsurf_new(i)=tslab(ki,1)
                  END IF 
              ELSE ! ice present
                  tsurf_new(i)=t_freeze
                  IF (tslab(ki,1).LT.t_freeze) THEN ! create new ice
                      ! quantity of new ice formed over open ocean
                      e_freeze=(t_freeze-tslab(ki,1))/cyang*(1.-fsic(ki)) &
                               /(ice_lat+ice_cap/2.*(t_freeze-tice(ki)))
                      ! new ice height and fraction
                      h_new=MIN(h_ice_new,seaice(ki)) ! max new height ice_new
                      dfsic=MIN(ice_frac_max-fsic(ki),e_freeze/h_new)
                      h_new=MIN(e_freeze/dfsic,h_ice_max)
                      ! update tslab to freezing over open ocean only
                      tslab(ki,1)=tslab(ki,1)*fsic(ki)+t_freeze*(1.-fsic(ki))
                      ! update sea ice
                      seaice(ki)=(h_new*dfsic+seaice(ki)*fsic(ki)) &
                                 /(dfsic+fsic(ki))
                      fsic(ki)=fsic(ki)+dfsic
                      ! update snow?
                  END IF ! tslab below freezing
              END IF ! sea ice present
          END DO
      END SELECT
      bils_cum(:)=0.0! clear cumulated fluxes
    END IF ! coupling time
  END SUBROUTINE ocean_slab_noice
!
!****************************************************************************************

  SUBROUTINE ocean_slab_ice(   &
       itime, dtime, jour, knon, knindex, &
       tsurf_in, p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, gustiness, &
       radsol, snow, qsurf, qsol, agesno, &
       alb1_new, alb2_new, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
       tsurf_new, dflux_s, dflux_l, swnet)

   USE calcul_fluxs_mod

   INCLUDE "YOMCST.h"
   INCLUDE "clesphys.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                  :: itime, jour, knon
    INTEGER, DIMENSION(klon), INTENT(IN) :: knindex
    REAL, INTENT(IN)                     :: dtime
    REAL, DIMENSION(klon), INTENT(IN)    :: tsurf_in
    REAL, DIMENSION(klon), INTENT(IN)    :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)    :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)    :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)    :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)    :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)    :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)    :: ps
    REAL, DIMENSION(klon), INTENT(IN)    :: u1, v1, gustiness
    REAL, DIMENSION(klon), INTENT(IN)    :: swnet

! In/Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)          :: snow, qsol
    REAL, DIMENSION(klon), INTENT(INOUT)          :: agesno
    REAL, DIMENSION(klon), INTENT(INOUT)          :: radsol

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)            :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)            :: alb1_new  ! new albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(OUT)            :: alb2_new  ! new albedo in near IR interval
    REAL, DIMENSION(klon), INTENT(OUT)            :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)            :: flux_u1, flux_v1
    REAL, DIMENSION(klon), INTENT(OUT)            :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)            :: dflux_s, dflux_l

! Local variables
!****************************************************************************************
    INTEGER               :: i,ki
    REAL, DIMENSION(klon) :: cal, beta, dif_grnd
    REAL, DIMENSION(klon) :: u0, v0
    REAL, DIMENSION(klon) :: u1_lay, v1_lay
    ! intermediate heat fluxes:
    REAL                  :: f_cond, f_swpen
    ! for snow/ice albedo:
    REAL                  :: alb_snow, alb_ice, alb_pond
    REAL                  :: frac_snow, frac_ice, frac_pond
    ! for ice melt / freeze
    REAL                  :: e_melt, snow_evap, h_test
    ! dhsic, dfsic change in ice mass, fraction.
    REAL                  :: dhsic, dfsic, frac_mf

!****************************************************************************************
! 1) Flux calculation
!****************************************************************************************
! Suppose zero surface speed
    u0(:)=0.0
    v0(:)=0.0
    u1_lay(:) = u1(:) - u0(:)
    v1_lay(:) = v1(:) - v0(:)

! set beta, cal, compute conduction fluxes inside ice/snow 
    slab_bilg(:)=0.
    dif_grnd(:)=0.
    beta(:) = 1.
    DO i=1,knon
    ki=knindex(i)
        IF (snow(i).GT.snow_min) THEN
            ! snow-layer heat capacity
            cal(i)=2.*RCPD/(snow(i)*ice_cap)
            ! snow conductive flux
            f_cond=sno_cond*(tice(ki)-tsurf_in(i))/snow(i)
            ! all shortwave flux absorbed
            f_swpen=0.
            ! bottom flux (ice conduction)
            slab_bilg(ki)=ice_cond*(tice(ki)-t_freeze)/seaice(ki)
            ! update ice temperature
            tice(ki)=tice(ki)-2./ice_cap/(snow(i)+seaice(ki)) &
                     *(slab_bilg(ki)+f_cond)*dtime
       ELSE ! bare ice
            ! ice-layer heat capacity
            cal(i)=2.*RCPD/(seaice(ki)*ice_cap)
            ! conductive flux
            f_cond=ice_cond*(t_freeze-tice(ki))/seaice(ki)
            ! penetrative shortwave flux...
            f_swpen=swnet(i)*pen_frac*exp(-pen_ext*seaice(ki)/ice_den)
            slab_bilg(ki)=f_swpen-f_cond
        END IF
        radsol(i)=radsol(i)+f_cond-f_swpen
    END DO
    ! weight fluxes to ocean by sea ice fraction
    slab_bilg(:)=slab_bilg(:)*fsic(:)

! calcul_fluxs (sens, lat etc)
    CALL calcul_fluxs(knon, is_sic, dtime, &
        tsurf_in, p1lay, cal, beta, cdragh, cdragh, ps, &
        precip_rain, precip_snow, snow, qsurf,  &
        radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, gustiness, &
        f_qsat_oce,AcoefH, AcoefQ, BcoefH, BcoefQ, &
        tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)
    DO i=1,knon 
        IF (snow(i).LT.snow_min) tice(knindex(i))=tsurf_new(i)
    END DO

! calcul_flux_wind
    CALL calcul_flux_wind(knon, dtime, &
         u0, v0, u1, v1, gustiness, cdragm, &
         AcoefU, AcoefV, BcoefU, BcoefV, &
         p1lay, temp_air, &
         flux_u1, flux_v1)

!****************************************************************************************
! 2) Update snow and ice surface
!****************************************************************************************
! snow precip
    DO i=1,knon
        ki=knindex(i)
        IF (precip_snow(i) > 0.) THEN
            snow(i) = snow(i)+precip_snow(i)*dtime*(1.-snow_wfact*(1.-fsic(ki)))
        END IF
! snow and ice sublimation
        IF (evap(i) > 0.) THEN 
           snow_evap = MIN (snow(i) / dtime, evap(i))
           snow(i) = snow(i) - snow_evap * dtime
           snow(i) = MAX(0.0, snow(i))
           seaice(ki) = MAX(0.0,seaice(ki)-(evap(i)-snow_evap)*dtime)
        ENDIF
! Melt / Freeze snow from above if Tsurf>0
        IF (tsurf_new(i).GT.t_melt) THEN
            ! energy available for melting snow (in kg of melted snow /m2)
            e_melt = MIN(MAX(snow(i)*(tsurf_new(i)-t_melt)*ice_cap/2. &
               /(ice_lat+ice_cap/2.*(t_melt-tice(ki))),0.0),snow(i))
            ! remove snow
            IF (snow(i).GT.e_melt) THEN
                snow(i)=snow(i)-e_melt
                tsurf_new(i)=t_melt
            ELSE ! all snow is melted
                ! add remaining heat flux to ice
                e_melt=e_melt-snow(i)
                tice(ki)=tice(ki)+e_melt*ice_lat*2./(ice_cap*seaice(ki))
                tsurf_new(i)=tice(ki)
            END IF
        END IF
! melt ice from above if Tice>0
        IF (tice(ki).GT.t_melt) THEN
            ! quantity of ice melted (kg/m2)
            e_melt=MAX(seaice(ki)*(tice(ki)-t_melt)*ice_cap/2. & 
             /(ice_lat+ice_cap/2.*(t_melt-t_freeze)),0.0)
            ! melt from above, height only
            dhsic=MIN(seaice(ki)-h_ice_min,e_melt)
            e_melt=e_melt-dhsic
            IF (e_melt.GT.0) THEN
            ! lateral melt if ice too thin
            dfsic=MAX(fsic(ki)-ice_frac_min,e_melt/h_ice_min*fsic(ki))
            ! if all melted add remaining heat to ocean
            e_melt=MAX(0.,e_melt*fsic(ki)-dfsic*h_ice_min)
            slab_bilg(ki)=slab_bilg(ki)+ e_melt*ice_lat/dtime
            ! update height and fraction
            fsic(ki)=fsic(ki)-dfsic
            END IF
            seaice(ki)=seaice(ki)-dhsic
            ! surface temperature at melting point
            tice(ki)=t_melt
            tsurf_new(i)=t_melt
        END IF
        ! convert snow to ice if below floating line
        h_test=(seaice(ki)+snow(i))*ice_den-seaice(ki)*sea_den
        IF (h_test.GT.0.) THEN !snow under water
            ! extra snow converted to ice (with added frozen sea water)
            dhsic=h_test/(sea_den-ice_den+sno_den)
            seaice(ki)=seaice(ki)+dhsic
            snow(i)=snow(i)-dhsic*sno_den/ice_den
            ! available energy (freeze sea water + bring to tice)
            e_melt=dhsic*(1.-sno_den/ice_den)*(ice_lat+ &
                   ice_cap/2.*(t_freeze-tice(ki)))
            ! update ice temperature
            tice(ki)=tice(ki)+2.*e_melt/ice_cap/(snow(i)+seaice(ki))
        END IF
    END DO

! New albedo
    DO i=1,knon
        ki=knindex(i)
       ! snow albedo: update snow age
        IF (snow(i).GT.0.0001) THEN
             agesno(i)=(agesno(i) + (1.-agesno(i)/50.)*dtime/86400.)&
                         * EXP(-1.*MAX(0.0,precip_snow(i))*dtime/5.)
        ELSE
            agesno(i)=0.0
        END IF
        ! snow albedo
        alb_snow=alb_sno_min+alb_sno_del*EXP(-agesno(i)/50.)
        ! ice albedo (varies with ice tkickness and temp)
        alb_ice=MAX(0.0,0.13*LOG(100.*seaice(ki)/ice_den)+0.1)
        IF (tice(ki).GT.t_freeze-0.01) THEN
            alb_ice=MIN(alb_ice,alb_ice_wet)
        ELSE
            alb_ice=MIN(alb_ice,alb_ice_dry)
        END IF
        ! pond albedo
        alb_pond=0.36-0.1*(2.0+MIN(0.0,MAX(tice(ki)-t_melt,-2.0)))
        ! pond fraction
        frac_pond=0.2*(2.0+MIN(0.0,MAX(tice(ki)-t_melt,-2.0)))
        ! snow fraction
        frac_snow=MAX(0.0,MIN(1.0-frac_pond,snow(i)/snow_min))
        ! ice fraction
        frac_ice=MAX(0.0,1.-frac_pond-frac_snow)
        ! total albedo
        alb1_new(i)=alb_snow*frac_snow+alb_ice*frac_ice+alb_pond*frac_pond
    END DO
    alb2_new(:) = alb1_new(:)

!****************************************************************************************
! 3) Recalculate new ocean temperature (add fluxes below ice)
!    Melt / freeze from below
!***********************************************o*****************************************
    !cumul fluxes
    bilg_cum(:)=bilg_cum(:)+slab_bilg(:)
    IF (MOD(itime,cpl_pas).EQ.0) THEN ! time to update tslab & fraction
        ! Add cumulated surface fluxes
        tslab(:,1)=tslab(:,1)+bilg_cum(:)*cyang*dtime
        DO i=1,knon
            ki=knindex(i)
            ! split lateral/top melt-freeze
            frac_mf=MIN(1.,MAX(0.,(seaice(ki)-h_ice_thin)/(h_ice_thick-h_ice_thin)))
            IF (tslab(ki,1).LE.t_freeze) THEN 
               ! ****** Form new ice from below *******
               ! quantity of new ice
                e_melt=(t_freeze-tslab(ki,1))/cyang & 
                       /(ice_lat+ice_cap/2.*(t_freeze-tice(ki)))
               ! first increase height to h_thin
               dhsic=MAX(0.,MIN(h_ice_thin-seaice(ki),e_melt/fsic(ki)))
               seaice(ki)=dhsic+seaice(ki)
               e_melt=e_melt-fsic(ki)*dhsic
               IF (e_melt.GT.0.) THEN
               ! frac_mf fraction used for lateral increase
               dfsic=MIN(ice_frac_max-fsic(ki),e_melt*frac_mf/seaice(ki))
               fsic(ki)=fsic(ki)+dfsic
               e_melt=e_melt-dfsic*seaice(ki)
               ! rest used to increase height
               seaice(ki)=MIN(h_ice_max,seaice(ki)+e_melt/fsic(ki))
               END IF
               tslab(ki,1)=t_freeze
           ELSE ! slab temperature above freezing
               ! ****** melt ice from below *******
               ! quantity of melted ice
               e_melt=(tslab(ki,1)-t_freeze)/cyang & 
                       /(ice_lat+ice_cap/2.*(tice(ki)-t_freeze))
               ! first decrease height to h_thick
               dhsic=MAX(0.,MIN(seaice(ki)-h_ice_thick,e_melt/fsic(ki)))
               seaice(ki)=seaice(ki)-dhsic
               e_melt=e_melt-fsic(ki)*dhsic
               IF (e_melt.GT.0) THEN
               ! frac_mf fraction used for height decrease
               dhsic=MAX(0.,MIN(seaice(ki)-h_ice_min,e_melt*frac_mf/fsic(ki)))
               seaice(ki)=seaice(ki)-dhsic
               e_melt=e_melt-fsic(ki)*dhsic
               ! rest used to decrease fraction (up to 0!)
               dfsic=MIN(fsic(ki),e_melt/seaice(ki))
               ! keep remaining in ocean
               e_melt=e_melt-dfsic*seaice(ki)
               END IF
               tslab(ki,1)=t_freeze+e_melt*ice_lat*cyang
               fsic(ki)=fsic(ki)-dfsic
           END IF
        END DO
        bilg_cum(:)=0.
    END IF ! coupling time
    
    !tests ice fraction 
    WHERE (fsic.LT.ice_frac_min)
        tslab(:,1)=tslab(:,1)-fsic*seaice*ice_lat*cyang
        tice=t_melt
        fsic=0.
        seaice=0.
    END WHERE

  END SUBROUTINE ocean_slab_ice
!
!****************************************************************************************
!
  SUBROUTINE ocean_slab_final

!****************************************************************************************
! Deallocate module variables
!****************************************************************************************
    IF (ALLOCATED(tslab)) DEALLOCATE(tslab)
    IF (ALLOCATED(fsic)) DEALLOCATE(fsic)
    IF (ALLOCATED(tice)) DEALLOCATE(tice)
    IF (ALLOCATED(seaice)) DEALLOCATE(seaice)
    IF (ALLOCATED(slab_bilg)) DEALLOCATE(slab_bilg)
    IF (ALLOCATED(bilg_cum)) DEALLOCATE(bilg_cum)
    IF (ALLOCATED(bils_cum)) DEALLOCATE(bils_cum)
    IF (ALLOCATED(taux_cum)) DEALLOCATE(taux_cum)
    IF (ALLOCATED(tauy_cum)) DEALLOCATE(tauy_cum)
    IF (ALLOCATED(dt_ekman)) DEALLOCATE(dt_ekman)
    IF (ALLOCATED(dt_hdiff)) DEALLOCATE(dt_hdiff)
    IF (ALLOCATED(dt_gm)) DEALLOCATE(dt_gm)
    IF (ALLOCATED(dt_qflux)) DEALLOCATE(dt_qflux)

  END SUBROUTINE ocean_slab_final

END MODULE ocean_slab_mod
