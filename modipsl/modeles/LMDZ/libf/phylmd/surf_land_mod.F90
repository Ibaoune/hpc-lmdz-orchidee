!
MODULE surf_land_mod

  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!  
  SUBROUTINE surf_land(itime, dtime, date0, jour, knon, knindex, &
       rlon, rlat, yrmu0, &
       debut, lafin, zlev, ccanopy, swnet, lwnet, albedo, &
       tsurf, p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, & 
       pref, u1, v1, gustiness, rugoro, pctsrf, &
       lwdown_m, q2m, t2m, &
       snow, qsol, agesno, tsoil, &
       z0m, z0h, SFRWL, alb_dir_new, alb_dif_new, evap, fluxsens, fluxlat, &   
       qsurf, tsurf_new, dflux_s, dflux_l, &
       flux_u1, flux_v1 , & 
       veget,lai,height)

    USE dimphy
    USE surface_data, ONLY    : ok_veget
! >> PC
    USE carbon_cycle_mod 
! << PC

    ! See comments in each module surf_land_orchidee_xxx for compatiblity with ORCHIDEE
#ifdef ORCHIDEE_NOOPENMP
    ! Compilation with cpp key ORCHIDEE NOOPENMP
    USE surf_land_orchidee_noopenmp_mod
#else
#if ORCHIDEE_NOZ0H
    ! Compilation with cpp key ORCHIDEE NOZ0H
    USE surf_land_orchidee_noz0h_mod
#else
#if ORCHIDEE_NOFREIN
    ! Compilation with cpp key ORCHIDEE_NOFREIN
    USE surf_land_orchidee_nofrein_mod
#else
#if ORCHIDEE_NOUNSTRUCT
    ! Compilation with cpp key ORCHIDEE_NOUNSTRUCT
    USE surf_land_orchidee_nounstruct_mod
#else
    USE surf_land_orchidee_mod
#endif
#endif
#endif
#endif

    USE surf_land_bucket_mod
    USE calcul_fluxs_mod
    USE indice_sol_mod

! >> PC
    USE print_control_mod, ONLY: lunout
! << PC

    INCLUDE "dimsoil.h"
    INCLUDE "YOMCST.h"
    INCLUDE "clesphys.h"
    INCLUDE "dimpft.h"

! Input variables  
!****************************************************************************************
    INTEGER, INTENT(IN)                     :: itime, jour, knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex
    REAL, INTENT(IN)                        :: date0
    REAL, DIMENSION(klon), INTENT(IN)       :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)       :: yrmu0  ! cosine of solar zenith angle
    LOGICAL, INTENT(IN)                     :: debut, lafin
    REAL, INTENT(IN)                        :: dtime
    REAL, DIMENSION(klon), INTENT(IN)       :: zlev, ccanopy
    REAL, DIMENSION(klon), INTENT(IN)       :: swnet, lwnet
    REAL, DIMENSION(klon), INTENT(IN)       :: albedo  ! albedo for whole short-wave interval
    REAL, DIMENSION(klon), INTENT(IN)       :: tsurf
    REAL, DIMENSION(klon), INTENT(IN)       :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)       :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)       :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)       :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)       :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)       :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)       :: pref   ! pressure reference
    REAL, DIMENSION(klon), INTENT(IN)       :: u1, v1, gustiness
    REAL, DIMENSION(klon), INTENT(IN)       :: rugoro
    REAL, DIMENSION(klon,nbsrf), INTENT(IN) :: pctsrf
    REAL, DIMENSION(klon), INTENT(IN)       :: lwdown_m  ! downwelling longwave radiation at mean surface
                                                         ! corresponds to previous sollwdown
    REAL, DIMENSION(klon), INTENT(IN)       :: q2m, t2m

! In/Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)          :: snow, qsol
    REAL, DIMENSION(klon), INTENT(INOUT)          :: agesno
    REAL, DIMENSION(klon, nsoilmx), INTENT(INOUT) :: tsoil

! Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: z0m, z0h
!albedo SB >>>
!    REAL, DIMENSION(klon), INTENT(OUT)       :: alb1_new ! albdeo for shortwave interval 1(visible)
!    REAL, DIMENSION(klon), INTENT(OUT)       :: alb2_new ! albedo for shortwave interval 2(near infrared)
    REAL, DIMENSION(6), INTENT(IN) :: SFRWL
    REAL, DIMENSION(klon,nsw), INTENT(OUT)       :: alb_dir_new,alb_dif_new
!albedo SB <<<
    REAL, DIMENSION(klon), INTENT(OUT)       :: evap
    REAL, DIMENSION(klon), INTENT(OUT)       :: fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)       :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)       :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: dflux_s, dflux_l      
    REAL, DIMENSION(klon), INTENT(OUT)       :: flux_u1, flux_v1  ! flux for U and V at first model level
    REAL, DIMENSION(klon,nvm_lmdz), INTENT(OUT) :: veget,lai
    REAL, DIMENSION(klon,nvm_lmdz), INTENT(OUT) :: height

! Local variables
!****************************************************************************************
    REAL, DIMENSION(klon) :: p1lay_tmp
    REAL, DIMENSION(klon) :: pref_tmp
    REAL, DIMENSION(klon) :: swdown     ! downwelling shortwave radiation at land surface
    REAL, DIMENSION(klon) :: epot_air           ! potential air temperature
    REAL, DIMENSION(klon) :: tsol_rad, emis_new ! output from interfsol not used
    REAL, DIMENSION(klon) :: u0, v0     ! surface speed
    INTEGER               :: i

!albedo SB >>>
    REAL, DIMENSION(klon)      :: alb1_new,alb2_new
!albedo SB <<<

!**************************************************************************************** 
! Choice between call to vegetation model (ok_veget=true) or simple calculation below
!
!****************************************************************************************
   IF (ok_veget) THEN
!****************************************************************************************
!  Call model sechiba in model ORCHIDEE
!
!****************************************************************************************
       p1lay_tmp(:)      = 0.0
       pref_tmp(:)       = 0.0
       p1lay_tmp(1:knon) = p1lay(1:knon)/100.
       pref_tmp(1:knon)  = pref(1:knon)/100.
! 
!* Calculate incoming flux for SW and LW interval: swdown
!
       swdown(:) = 0.0
       DO i = 1, knon
          swdown(i) = swnet(i)/(1-albedo(i))
       END DO
!
!* Calculate potential air temperature
!
       epot_air(:) = 0.0
       DO i = 1, knon
          epot_air(i) = RCPD*temp_air(i)*(pref(i)/p1lay(i))**RKAPPA
       END DO

       ! temporary for keeping same results using lwdown_m instead of lwdown
       CALL surf_land_orchidee(itime, dtime, date0, knon, &
            knindex, rlon, rlat, yrmu0, pctsrf, &
            debut, lafin, &
            zlev,  u1, v1, gustiness, temp_air, spechum, epot_air, ccanopy, & 
            cdragh, AcoefH, AcoefQ, BcoefH, BcoefQ, &
            precip_rain, precip_snow, lwdown_m, swnet, swdown, &
            pref_tmp, q2m, t2m, &
            evap, fluxsens, fluxlat, &              
            tsol_rad, tsurf_new, alb1_new, alb2_new, &
            emis_new, z0m, z0h, qsurf, &
            veget, lai, height)       
!  
!* Add contribution of relief to surface roughness
!  
       DO i=1,knon
          z0m(i) = MAX(1.5e-05,SQRT(z0m(i)**2 + rugoro(i)**2))
       ENDDO

    ELSE  ! not ok_veget
!****************************************************************************************
! No extern vegetation model choosen, call simple bucket calculations instead.
!
!****************************************************************************************
       CALL surf_land_bucket(itime, jour, knon, knindex, debut, dtime,&
            tsurf, p1lay, cdragh, precip_rain, precip_snow, temp_air, &
            spechum, AcoefH, AcoefQ, BcoefH, BcoefQ, pref, &
            u1, v1, gustiness, rugoro, swnet, lwnet, &
            snow, qsol, agesno, tsoil, &
            qsurf, z0m, alb1_new, alb2_new, evap, &
            fluxsens, fluxlat, tsurf_new, dflux_s, dflux_l)
        z0h(1:knon)=z0m(1:knon) ! En attendant mieux

    ENDIF ! ok_veget

!****************************************************************************************
! Calculation for all land models
! - Flux calculation at first modele level for U and V
!****************************************************************************************
! Suppose zero surface speed
    u0(:)=0.0
    v0(:)=0.0
    CALL calcul_flux_wind(knon, dtime, &
         u0, v0, u1, v1, gustiness, cdragm, &
         AcoefU, AcoefV, BcoefU, BcoefV, &
         p1lay, temp_air, &
         flux_u1, flux_v1)

!albedo SB >>>
     SELECT CASE(NSW)
     CASE(2)
       alb_dir_new(1:knon,1)=alb1_new(1:knon)
       alb_dir_new(1:knon,2)=alb2_new(1:knon)
     CASE(4)
       alb_dir_new(1:knon,1)=alb1_new(1:knon)
       alb_dir_new(1:knon,2)=alb2_new(1:knon)
       alb_dir_new(1:knon,3)=alb2_new(1:knon)
       alb_dir_new(1:knon,4)=alb2_new(1:knon)
     CASE(6)
       alb_dir_new(1:knon,1)=alb1_new(1:knon)
       alb_dir_new(1:knon,2)=alb1_new(1:knon)
       alb_dir_new(1:knon,3)=alb1_new(1:knon)
       alb_dir_new(1:knon,4)=alb2_new(1:knon)
       alb_dir_new(1:knon,5)=alb2_new(1:knon)
       alb_dir_new(1:knon,6)=alb2_new(1:knon)
     END SELECT

     alb_dif_new=alb_dir_new
!albedo SB <<<
    
  END SUBROUTINE surf_land
!
!****************************************************************************************
!  
END MODULE surf_land_mod
!
!****************************************************************************************
!  
