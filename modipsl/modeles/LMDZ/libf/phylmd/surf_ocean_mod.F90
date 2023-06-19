!
! $Id: surf_ocean_mod.F90 3395 2018-09-25 15:22:13Z lguez $
!
MODULE surf_ocean_mod

  IMPLICIT NONE

CONTAINS
  !
  !******************************************************************************
  !
  SUBROUTINE surf_ocean(rlon, rlat, swnet, lwnet, alb1, &
       windsp, rmu0, fder, tsurf_in, &
       itime, dtime, jour, knon, knindex, &
       p1lay, z1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, gustiness, rugoro, pctsrf, &
       snow, qsurf, agesno, &
       z0m, z0h, SFRWL, alb_dir_new, alb_dif_new, evap, fluxsens, fluxlat, & 
       tsurf_new, dflux_s, dflux_l, lmt_bils, &
       flux_u1, flux_v1)

    use albedo, only: alboc, alboc_cd
    USE dimphy, ONLY: klon, zmasq
    USE surface_data, ONLY     : type_ocean
    USE ocean_forced_mod, ONLY : ocean_forced_noice
    USE ocean_slab_mod, ONLY   : ocean_slab_noice
    USE ocean_cpl_mod, ONLY    : ocean_cpl_noice
    USE indice_sol_mod, ONLY : nbsrf, is_oce
    USE limit_read_mod
    !
    ! This subroutine will make a call to ocean_XXX_noice according to the ocean mode (force, 
    ! slab or couple). The calculations of albedo and rugosity for the ocean surface are 
    ! done in here because they are identical for the different modes of ocean. 


    INCLUDE "YOMCST.h"

    include "clesphys.h"
    ! for cycle_diurne and for iflag_z0_oce==-1 (prescribed z0)

    ! Input variables
    !******************************************************************************
    INTEGER, INTENT(IN)                      :: itime, jour, knon
    INTEGER, DIMENSION(klon), INTENT(IN)     :: knindex
    REAL, INTENT(IN)                         :: dtime
    REAL, DIMENSION(klon), INTENT(IN)        :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)        :: swnet  ! net shortwave radiation at surface  
    REAL, DIMENSION(klon), INTENT(IN)        :: lwnet  ! net longwave radiation at surface  
    REAL, DIMENSION(klon), INTENT(IN)        :: alb1   ! albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(IN)        :: windsp
    REAL, DIMENSION(klon), INTENT(IN)        :: rmu0  
    REAL, DIMENSION(klon), INTENT(IN)        :: fder
    REAL, DIMENSION(klon), INTENT(IN)        :: tsurf_in
    REAL, DIMENSION(klon), INTENT(IN)        :: p1lay,z1lay ! pression (Pa) et altitude (m) du premier niveau
    REAL, DIMENSION(klon), INTENT(IN)        :: cdragh
    REAL, DIMENSION(klon), INTENT(IN)        :: cdragm
    REAL, DIMENSION(klon), INTENT(IN)        :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)        :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)        :: ps
    REAL, DIMENSION(klon), INTENT(IN)        :: u1, v1, gustiness
    REAL, DIMENSION(klon), INTENT(IN)        :: rugoro
    REAL, DIMENSION(klon,nbsrf), INTENT(IN)  :: pctsrf

    ! In/Output variables
    !******************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)     :: snow
    REAL, DIMENSION(klon), INTENT(INOUT)     :: qsurf
    REAL, DIMENSION(klon), INTENT(INOUT)     :: agesno
    REAL, DIMENSION(klon), INTENT(inOUT):: z0h

    ! Output variables
    !******************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: z0m
    !albedo SB >>>
    !    REAL, DIMENSION(klon), INTENT(OUT)       :: alb1_new  ! new albedo in visible SW interval
    !    REAL, DIMENSION(klon), INTENT(OUT)       :: alb2_new  ! new albedo in near IR interval
    REAL, DIMENSION(6), INTENT(IN)          :: SFRWL 
    REAL, DIMENSION(klon,nsw), INTENT(OUT)       :: alb_dir_new,alb_dif_new
    !albedo SB <<<     
    REAL, DIMENSION(klon), INTENT(OUT)       :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)       :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: dflux_s, dflux_l      
    REAL, DIMENSION(klon), INTENT(OUT)       :: lmt_bils
    REAL, DIMENSION(klon), INTENT(OUT)       :: flux_u1, flux_v1

    ! Local variables
    !******************************************************************************
    INTEGER               :: i, k
    REAL                  :: tmp
    REAL, PARAMETER       :: cepdu2=(0.1)**2
    REAL, DIMENSION(klon) :: alb_eau, z0_lim
    REAL, DIMENSION(klon) :: radsol
    REAL, DIMENSION(klon) :: cdragq ! Cdrag pour l'evaporation
    CHARACTER(len=20),PARAMETER :: modname="surf_ocean"

    ! End definition
    !******************************************************************************


    !******************************************************************************
    ! Calculate total net radiance at surface
    !
    !******************************************************************************
    radsol(1:klon) = 0.0 ! initialisation a priori inutile
    radsol(1:knon) = swnet(1:knon) + lwnet(1:knon)

    !******************************************************************************
    ! Cdragq computed from cdrag
    ! The difference comes only from a factor (f_z0qh_oce) on z0, so that
    ! it can be computed inside surf_ocean
    ! More complicated appraches may require the propagation through
    ! pbl_surface of an independant cdragq variable.
    !******************************************************************************

    IF ( f_z0qh_oce .ne. 1.) THEN
       ! Si on suit les formulations par exemple de Tessel, on 
       ! a z0h=0.4*nu/u*, z0q=0.62*nu/u*, d'ou f_z0qh_oce=0.62/0.4=1.55
       cdragq(1:knon)=cdragh(1:knon)*                                      &
            log(z1lay(1:knon)/z0h(1:knon))/log(z1lay(1:knon)/(f_z0qh_oce*z0h(1:knon)))
    ELSE
       cdragq(1:knon)=cdragh(1:knon)
    ENDIF


    !******************************************************************************
    ! Switch according to type of ocean (couple, slab or forced)
    !******************************************************************************
    SELECT CASE(type_ocean)
    CASE('couple')
       CALL ocean_cpl_noice( &
            swnet, lwnet, alb1, &
            windsp, fder, & 
            itime, dtime, knon, knindex, &
            p1lay, cdragh, cdragq, cdragm, precip_rain, precip_snow,temp_air,spechum,& 
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            ps, u1, v1, gustiness, &
            radsol, snow, agesno, &
            qsurf, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
            tsurf_new, dflux_s, dflux_l)

    CASE('slab')
       CALL ocean_slab_noice( &
            itime, dtime, jour, knon, knindex, &
            p1lay, cdragh, cdragq, cdragm, precip_rain, precip_snow, temp_air, spechum,&
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            ps, u1, v1, gustiness, tsurf_in, &
            radsol, snow, &
            qsurf, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
            tsurf_new, dflux_s, dflux_l, lmt_bils)

    CASE('force')
       CALL ocean_forced_noice( &
            itime, dtime, jour, knon, knindex, &
            p1lay, cdragh, cdragq, cdragm, precip_rain, precip_snow, &
            temp_air, spechum, &
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            ps, u1, v1, gustiness, &
            radsol, snow, agesno, &
            qsurf, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
            tsurf_new, dflux_s, dflux_l)
    END SELECT

    !******************************************************************************
    ! fcodron: compute lmt_bils  forced case (same as wfbils_oce / 1.-contfracatm)
    !******************************************************************************
    IF (type_ocean.NE.'slab') THEN
       lmt_bils(1:klon)=0.
       DO i=1,knon
          lmt_bils(knindex(i))=(swnet(i)+lwnet(i)+fluxsens(i)+fluxlat(i)) &
               *pctsrf(knindex(i),is_oce)/(1.-zmasq(knindex(i))) 
       END DO
    END IF

    !******************************************************************************
    ! Calculate ocean surface albedo
    !******************************************************************************
    !albedo SB >>>
    IF (iflag_albedo==0) THEN
       !--old parametrizations of ocean surface albedo
       !
       IF (iflag_cycle_diurne.GE.1) THEN
          !
          CALL alboc_cd(rmu0,alb_eau)
          !
          !--ad-hoc correction for model radiative balance tuning
          !--now outside alboc_cd routine
          alb_eau(1:klon) = fmagic*alb_eau(1:klon) + pmagic
          alb_eau(1:klon)=MIN(MAX(alb_eau(1:klon),0.0),1.0)
          !
       ELSE
          !
          CALL alboc(REAL(jour),rlat,alb_eau)
          !--ad-hoc correction for model radiative balance tuning
          !--now outside alboc routine
          alb_eau(1:klon) = fmagic*alb_eau(1:klon) + pmagic
          alb_eau(1:klon)=MIN(MAX(alb_eau(1:klon),0.04),0.60)
          !
       ENDIF
       !
       DO i =1, knon
          DO  k=1,nsw
             alb_dir_new(i,k) = alb_eau(knindex(i))
          ENDDO
       ENDDO
       !IM 09122015 next line corresponds to the old way of doing in LMDZ5A/IPSLCM5A versions 
       !albedo for diffuse radiation is taken the same as for direct radiation
       alb_dif_new(1:knon,:)=alb_dir_new(1:knon,:)
       !IM 09122015 end
       !
    ELSE IF (iflag_albedo==1) THEN
       !--new parametrization of ocean surface albedo by Sunghye Baek
       !--albedo for direct and diffuse radiation are different
       !
       CALL ocean_albedo(knon,rmu0,knindex,windsp,SFRWL,alb_dir_new,alb_dif_new)
       !
       !--ad-hoc correction for model radiative balance tuning
       alb_dir_new(1:knon,:) = fmagic*alb_dir_new(1:knon,:) + pmagic
       alb_dif_new(1:knon,:) = fmagic*alb_dif_new(1:knon,:) + pmagic
       alb_dir_new(1:knon,:)=MIN(MAX(alb_dir_new(1:knon,:),0.0),1.0)
       alb_dif_new(1:knon,:)=MIN(MAX(alb_dif_new(1:knon,:),0.0),1.0)
       !
    ELSE IF (iflag_albedo==2) THEN
       ! F. Codron albedo read from limit.nc
       CALL limit_read_rug_alb(itime, dtime, jour,&
            knon, knindex, z0_lim, alb_eau)
       DO i =1, knon
          DO  k=1,nsw
             alb_dir_new(i,k) = alb_eau(i)
          ENDDO
       ENDDO
       alb_dif_new=alb_dir_new
    ENDIF
    !albedo SB <<<

    !******************************************************************************
    ! Calculate the rugosity
    !******************************************************************************
    IF (iflag_z0_oce==0) THEN
       DO i = 1, knon
          tmp = MAX(cepdu2,gustiness(i)+u1(i)**2+v1(i)**2)
          z0m(i) = 0.018*cdragm(i) * (gustiness(i)+u1(i)**2+v1(i)**2)/RG  &
               +  0.11*14e-6 / SQRT(cdragm(i) * tmp)
          z0m(i) = MAX(1.5e-05,z0m(i))
       ENDDO
       z0h(1:knon)=z0m(1:knon) ! En attendant mieux

    ELSE IF (iflag_z0_oce==1) THEN
       DO i = 1, knon
          tmp = MAX(cepdu2,gustiness(i)+u1(i)**2+v1(i)**2)
          z0m(i) = 0.018*cdragm(i) * (gustiness(i)+u1(i)**2+v1(i)**2)/RG  &
               + 0.11*14e-6 / SQRT(cdragm(i) * tmp)
          z0m(i) = MAX(1.5e-05,z0m(i))
          z0h(i)=0.4*14e-6 / SQRT(cdragm(i) * tmp)
       ENDDO
    ELSE IF (iflag_z0_oce==-1) THEN
       DO i = 1, knon
          z0m(i) = z0min
          z0h(i) = z0min
       ENDDO
    ELSE
       CALL abort_physic(modname,'version non prevue',1)
    ENDIF
    !
    !******************************************************************************
  END SUBROUTINE surf_ocean
  !******************************************************************************
  !
END MODULE surf_ocean_mod
