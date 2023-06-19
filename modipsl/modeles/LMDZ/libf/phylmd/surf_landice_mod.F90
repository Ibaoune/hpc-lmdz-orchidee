!
MODULE surf_landice_mod
  
  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE surf_landice(itime, dtime, knon, knindex, &
       rlon, rlat, debut, lafin, &
       rmu0, lwdownm, albedo, pphi1, &
       swnet, lwnet, tsurf, p1lay, &
       cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, gustiness, rugoro, pctsrf, &
       snow, qsurf, qsol, agesno, &
       tsoil, z0m, z0h, SFRWL, alb_dir, alb_dif, evap, fluxsens, fluxlat, &
       tsurf_new, dflux_s, dflux_l, &
       slope, cloudf, &
       snowhgt, qsnow, to_ice, sissnow, &
       alb3, runoff, &
       flux_u1, flux_v1)

    USE dimphy
    USE surface_data,     ONLY : type_ocean, calice, calsno, ok_snow
    USE fonte_neige_mod,  ONLY : fonte_neige, run_off_lic
    USE cpl_mod,          ONLY : cpl_send_landice_fields
    USE calcul_fluxs_mod
    USE phys_output_var_mod
!FC
    USE ioipsl_getin_p_mod, ONLY : getin_p

#ifdef CPP_SISVAT
    USE surf_sisvat_mod,  ONLY : surf_sisvat
#endif
    USE indice_sol_mod

!    INCLUDE "indicesol.h"
    INCLUDE "dimsoil.h"
    INCLUDE "YOMCST.h"
    INCLUDE "clesphys.h"

! Input variables 
!****************************************************************************************
    INTEGER, INTENT(IN)                           :: itime, knon
    INTEGER, DIMENSION(klon), INTENT(in)          :: knindex
    REAL, INTENT(in)                              :: dtime
    REAL, DIMENSION(klon), INTENT(IN)             :: swnet ! net shortwave radiance
    REAL, DIMENSION(klon), INTENT(IN)             :: lwnet ! net longwave radiance
    REAL, DIMENSION(klon), INTENT(IN)             :: tsurf
    REAL, DIMENSION(klon), INTENT(IN)             :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)             :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)             :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)             :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)             :: AcoefH, AcoefQ
    REAL, DIMENSION(klon), INTENT(IN)             :: BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)             :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)             :: ps
    REAL, DIMENSION(klon), INTENT(IN)             :: u1, v1, gustiness
    REAL, DIMENSION(klon), INTENT(IN)             :: rugoro
    REAL, DIMENSION(klon,nbsrf), INTENT(IN)       :: pctsrf

    LOGICAL,  INTENT(IN)                          :: debut   !true if first step
    LOGICAL,  INTENT(IN)                          :: lafin   !true if last step
    REAL, DIMENSION(klon), INTENT(IN)             :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)             :: rmu0
    REAL, DIMENSION(klon), INTENT(IN)             :: lwdownm !ylwdown
    REAL, DIMENSION(klon), INTENT(IN)             :: albedo  !mean albedo
    REAL, DIMENSION(klon), INTENT(IN)             :: pphi1    
    REAL, DIMENSION(klon), INTENT(IN)             :: slope   !mean slope in grid box  
    REAL, DIMENSION(klon), INTENT(IN)             :: cloudf  !total cloud fraction

! In/Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)          :: snow, qsol
    REAL, DIMENSION(klon), INTENT(INOUT)          :: agesno
    REAL, DIMENSION(klon, nsoilmx), INTENT(INOUT) :: tsoil

! Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)            :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)            :: z0m, z0h
!albedo SB >>>
!    REAL, DIMENSION(klon), INTENT(OUT)            :: alb1  ! new albedo in visible SW interval
!    REAL, DIMENSION(klon), INTENT(OUT)            :: alb2  ! new albedo in near IR interval
    REAL, DIMENSION(6), INTENT(IN)              ::SFRWL
    REAL, DIMENSION(klon,nsw), INTENT(OUT)        ::alb_dir,alb_dif
!albedo SB <<<
    REAL, DIMENSION(klon), INTENT(OUT)            :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)            :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)            :: dflux_s, dflux_l      
    REAL, DIMENSION(klon), INTENT(OUT)            :: flux_u1, flux_v1

    REAL, DIMENSION(klon), INTENT(OUT)           :: alb3
    REAL, DIMENSION(klon), INTENT(OUT)           :: qsnow   !column water in snow [kg/m2]
    REAL, DIMENSION(klon), INTENT(OUT)           :: snowhgt !Snow height (m)
    REAL, DIMENSION(klon), INTENT(OUT)           :: to_ice
    REAL, DIMENSION(klon), INTENT(OUT)           :: sissnow
    REAL, DIMENSION(klon), INTENT(OUT)           :: runoff  !Land ice runoff
 

! Local variables
!****************************************************************************************
    REAL, DIMENSION(klon)    :: soilcap, soilflux
    REAL, DIMENSION(klon)    :: cal, beta, dif_grnd
    REAL, DIMENSION(klon)    :: zfra, alb_neig
    REAL, DIMENSION(klon)    :: radsol
    REAL, DIMENSION(klon)    :: u0, v0, u1_lay, v1_lay
    INTEGER                  :: i,j

    REAL, DIMENSION(klon)    :: emis_new                  !Emissivity
    REAL, DIMENSION(klon)    :: swdown,lwdown
    REAL, DIMENSION(klon)    :: precip_snow_adv, snow_adv !Snow Drift precip./advection
    REAL, DIMENSION(klon)    :: bl_height, wind_velo      !height boundary layer, wind spd
    REAL, DIMENSION(klon)    :: dens_air,  snow_cont_air  !air density; snow content air
    REAL, DIMENSION(klon)    :: alb_soil                  !albedo of underlying ice
    REAL, DIMENSION(klon)    :: pexner                    !Exner potential
    REAL                     :: pref
    REAL, DIMENSION(klon,nsoilmx) :: tsoil0 !modfi

    CHARACTER (len = 20)                      :: modname = 'surf_landice'
    CHARACTER (len = 80)                      :: abort_message


!albedo SB >>>
    real,dimension(klon) :: alb1,alb2
!albedo SB <<<

! End definition
!****************************************************************************************
!FC 
!FC
   REAL,SAVE :: alb_vis_sno_lic
  !$OMP THREADPRIVATE(alb_vis_sno_lic)
   REAL,SAVE :: alb_nir_sno_lic
  !$OMP THREADPRIVATE(alb_nir_sno_lic)
  LOGICAL, SAVE :: firstcall = .TRUE.
  !$OMP THREADPRIVATE(firstcall)
!FC


  IF (firstcall) THEN
  alb_vis_sno_lic=0.77
  CALL getin_p('alb_vis_sno_lic',alb_vis_sno_lic)
           PRINT*, 'alb_vis_sno_lic',alb_vis_sno_lic
  alb_nir_sno_lic=0.77
  CALL getin_p('alb_nir_sno_lic',alb_nir_sno_lic)
           PRINT*, 'alb_nir_sno_lic',alb_nir_sno_lic
  firstcall=.false.
  ENDIF
!
! Initialize output variables
    alb3(:) = 999999.
    alb2(:) = 999999.
    alb1(:) = 999999.
    
    runoff(:) = 0.
!****************************************************************************************
! Calculate total absorbed radiance at surface
!
!****************************************************************************************
    radsol(:) = 0.0
    radsol(1:knon) = swnet(1:knon) + lwnet(1:knon)

!****************************************************************************************
!   ok_snow = TRUE  : prepare and call SISVAT snow model
!   ok_snow = FALSE : soil_model, calcul_flux, fonte_neige, ... 
!
!****************************************************************************************
    IF (ok_snow) THEN
#ifdef CPP_SISVAT
       ! Prepare for calling SISVAT
       
       ! Calculate incoming flux for SW and LW interval: swdown, lwdown
       swdown(:)        = 0.0
       lwdown(:)        = 0.0
       DO i = 1, knon
          swdown(i)        = swnet(i)/(1-albedo(i))
          lwdown(i)        = lwdownm(i)
       END DO
       
       ! Set constants and compute some input for SISVAT
       snow_adv(:)      = 0.                          ! no snow blown in for now
       snow_cont_air(:) = 0.       
       alb_soil(:)      = albedo(:)
       pref             = 100000.                     ! = 1000 hPa
       DO i = 1, knon
          wind_velo(i)     = u1(i)**2 + v1(i)**2
          wind_velo(i)     = wind_velo(i)**0.5
          pexner(i)        = (p1lay(i)/pref)**(RD/RCPD)
          dens_air(i)      = p1lay(i)/RD/temp_air(i)  ! dry air density
          bl_height(i)     = pphi1(i)/RG              
       END DO

!****************************************************************************************
! CALL to SISVAT interface
!
!****************************************************************************************
       ! config: compute everything with SV but temperatures afterwards with soil/calculfluxs
       DO i = 1, knon
          tsoil0(i,:)=tsoil(i,:)
       END DO
           ! Martin
           PRINT*, 'on appelle surf_sisvat'
           ! Martin
       CALL surf_sisvat(knon, rlon, rlat, knindex, itime, dtime, debut, lafin, &
            rmu0, swdown, lwdown, pexner, ps, p1lay, &
            precip_rain, precip_snow, precip_snow_adv, snow_adv, &
            bl_height, wind_velo, temp_air, dens_air, spechum, tsurf, &
            rugoro, snow_cont_air, alb_soil, slope, cloudf, &
            radsol, qsol, tsoil0, snow, snowhgt, qsnow, to_ice,sissnow, agesno, &
            AcoefH, AcoefQ, BcoefH, BcoefQ, cdragh, &
            run_off_lic, evap, fluxsens, fluxlat, dflux_s, dflux_l, &        
            tsurf_new, alb1, alb2, alb3, &
            emis_new, z0m, qsurf)
       z0h(1:knon)=z0m(1:knon) ! en attendant mieux
       
       ! Suppose zero surface speed
       u0(:)            = 0.0
       v0(:)            = 0.0
       ! The calculation of heat/water fluxes, otherwise done by "CALL calcul_fluxs" is 
       ! integrated in SISVAT, using the same method. It can be found in "sisvat.f", in the
       ! subroutine "SISVAT_TS2". 
       ! u0, v0=0., dif_grnd=0. and beta=1 are assumed there! 
       
       CALL calcul_flux_wind(knon, dtime, &
            u0, v0, u1, v1, gustiness, cdragm, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            p1lay, temp_air, &
            flux_u1, flux_v1)
#else
       abort_message='Pb de coherence: ok_snow = .true. mais CPP_SISVAT = .false.'
       CALL abort_physic(modname,abort_message,1)
#endif
    ELSE ! ok_snow=FALSE

!****************************************************************************************
! Soil calculations
! 
!****************************************************************************************
    IF (soil_model) THEN 
       CALL soil(dtime, is_lic, knon, snow, tsurf, tsoil, soilcap, soilflux)
       cal(1:knon) = RCPD / soilcap(1:knon)
       radsol(1:knon)  = radsol(1:knon) + soilflux(1:knon)
    ELSE 
       cal = RCPD * calice
       WHERE (snow > 0.0) cal = RCPD * calsno
    ENDIF


!****************************************************************************************
! Calulate fluxes
!
!****************************************************************************************
    beta(:) = 1.0
    dif_grnd(:) = 0.0

! Suppose zero surface speed
    u0(:)=0.0
    v0(:)=0.0
    u1_lay(:) = u1(:) - u0(:)
    v1_lay(:) = v1(:) - v0(:)

    CALL calcul_fluxs(knon, is_lic, dtime, &
         tsurf, p1lay, cal, beta, cdragh, cdragh, ps, &
         precip_rain, precip_snow, snow, qsurf,  &
         radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, gustiness, &
         1.,AcoefH, AcoefQ, BcoefH, BcoefQ, &
         tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

    CALL calcul_flux_wind(knon, dtime, &
         u0, v0, u1, v1, gustiness, cdragm, &
         AcoefU, AcoefV, BcoefU, BcoefV, &
         p1lay, temp_air, &
         flux_u1, flux_v1)

!****************************************************************************************
! Calculate snow height, age, run-off,..
!    
!****************************************************************************************
    CALL fonte_neige( knon, is_lic, knindex, dtime, &
         tsurf, precip_rain, precip_snow, &
         snow, qsol, tsurf_new, evap)


!****************************************************************************************
! Calculate albedo
!
!****************************************************************************************
    CALL albsno(klon,knon,dtime,agesno(:),alb_neig(:), precip_snow(:))  
    WHERE (snow(1 : knon) .LT. 0.0001) agesno(1 : knon) = 0.
    zfra(1:knon) = MAX(0.0,MIN(1.0,snow(1:knon)/(snow(1:knon)+10.0)))
    alb1(1:knon) = alb_neig(1:knon)*zfra(1:knon) + &
         0.6 * (1.0-zfra(1:knon))
!
!IM: plusieurs choix/tests sur l'albedo des "glaciers continentaux"
!       alb1(1 : knon)  = 0.6 !IM cf FH/GK 
!       alb1(1 : knon)  = 0.82
!       alb1(1 : knon)  = 0.77 !211003 Ksta0.77
!       alb1(1 : knon)  = 0.8 !KstaTER0.8 & LMD_ARMIP5
!IM: KstaTER0.77 & LMD_ARMIP6    

! Attantion: alb1 and alb2 are the same!
    alb1(1:knon)  = alb_vis_sno_lic
    alb2(1:knon)  = alb_nir_sno_lic


!****************************************************************************************
! Rugosity
!
!****************************************************************************************
    z0m=1.e-3
    z0h = z0m
    z0m = SQRT(z0m**2+rugoro**2)

    END IF ! ok_snow


!****************************************************************************************
! Send run-off on land-ice to coupler if coupled ocean.
! run_off_lic has been calculated in fonte_neige or surf_sisvat
!
!****************************************************************************************
    IF (type_ocean=='couple') THEN
       CALL cpl_send_landice_fields(itime, knon, knindex, run_off_lic)
    ENDIF

 ! transfer runoff rate [kg/m2/s](!) to physiq for output
    runoff(1:knon)=run_off_lic(1:knon)/dtime

  
!****************************************************************************************
       snow_o=0.
       zfra_o = 0.
       DO j = 1, knon
           i = knindex(j)
           snow_o(i) = snow(j) 
           zfra_o(i) = zfra(j)
       ENDDO

!****************************************************************************************
       snow_o=0.
       zfra_o = 0.
       DO j = 1, knon
           i = knindex(j)
           snow_o(i) = snow(j)
           zfra_o(i) = zfra(j)
       ENDDO


!albedo SB >>>
     select case(NSW)
     case(2)
       alb_dir(1:knon,1)=alb1(1:knon)
       alb_dir(1:knon,2)=alb2(1:knon)
     case(4)
       alb_dir(1:knon,1)=alb1(1:knon)
       alb_dir(1:knon,2)=alb2(1:knon)
       alb_dir(1:knon,3)=alb2(1:knon)
       alb_dir(1:knon,4)=alb2(1:knon)
     case(6)
       alb_dir(1:knon,1)=alb1(1:knon)
       alb_dir(1:knon,2)=alb1(1:knon)
       alb_dir(1:knon,3)=alb1(1:knon)
       alb_dir(1:knon,4)=alb2(1:knon)
       alb_dir(1:knon,5)=alb2(1:knon)
       alb_dir(1:knon,6)=alb2(1:knon)
     end select
alb_dif=alb_dir
!albedo SB <<<




  END SUBROUTINE surf_landice
!
!****************************************************************************************
!
END MODULE surf_landice_mod



