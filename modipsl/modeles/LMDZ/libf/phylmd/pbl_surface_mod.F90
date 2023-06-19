!
! $Id: pbl_surface_mod.F90 3838 2021-02-11 09:39:04Z musat $
!
MODULE pbl_surface_mod
!
! Planetary Boundary Layer and Surface module
!
! This module manages the calculation of turbulent diffusion in the boundary layer 
! and all interactions towards the differents sub-surfaces.
!
!
  USE dimphy
  USE mod_phys_lmdz_para,  ONLY : mpi_size
  USE mod_grid_phy_lmdz,   ONLY : klon_glo
  USE ioipsl
  USE surface_data,        ONLY : type_ocean, ok_veget
  USE surf_land_mod,       ONLY : surf_land
  USE surf_landice_mod,    ONLY : surf_landice
  USE surf_ocean_mod,      ONLY : surf_ocean
  USE surf_seaice_mod,     ONLY : surf_seaice
  USE cpl_mod,             ONLY : gath2cpl
  USE climb_hq_mod,        ONLY : climb_hq_down, climb_hq_up
  USE climb_wind_mod,      ONLY : climb_wind_down, climb_wind_up
  USE coef_diff_turb_mod,  ONLY : coef_diff_turb
  USE ioipsl_getin_p_mod,  ONLY : getin_p
  USE cdrag_mod
  USE stdlevvar_mod
  USE wx_pbl_mod,          ONLY : wx_pbl_init, wx_pbl_final, &
!!                                  wx_pbl_fuse_no_dts, wx_pbl_split_no_dts, &
!!                                  wx_pbl_fuse, wx_pbl_split
                                  wx_pbl0_fuse, wx_pbl0_split

  IMPLICIT NONE

! Declaration of variables saved in restart file
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE     :: fder   ! flux drift
  !$OMP THREADPRIVATE(fder)
  REAL, ALLOCATABLE, DIMENSION(:,:), PUBLIC, SAVE   :: snow   ! snow at surface
  !$OMP THREADPRIVATE(snow)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE, SAVE   :: qsurf  ! humidity at surface
  !$OMP THREADPRIVATE(qsurf)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: ftsoil ! soil temperature
  !$OMP THREADPRIVATE(ftsoil)

  INTEGER, SAVE :: iflag_pbl_surface_t2m_bug
  !$OMP THREADPRIVATE(iflag_pbl_surface_t2m_bug)
  integer, save :: iflag_new_t2mq2m
  !$OMP THREADPRIVATE(iflag_new_t2mq2m)
!FC
!  integer, save :: iflag_frein
!  !$OMP THREADPRIVATE(iflag_frein)

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE pbl_surface_init(fder_rst, snow_rst, qsurf_rst, ftsoil_rst)

! This routine should be called after the restart file has been read.
! This routine initialize the restart variables and does some validation tests
! for the index of the different surfaces and tests the choice of type of ocean.

    USE indice_sol_mod
    USE print_control_mod, ONLY: lunout
    USE ioipsl_getin_p_mod, ONLY : getin_p
    IMPLICIT NONE

    INCLUDE "dimsoil.h"
 
! Input variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(IN)                 :: fder_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(IN)          :: snow_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(IN)          :: qsurf_rst
    REAL, DIMENSION(klon, nsoilmx, nbsrf), INTENT(IN) :: ftsoil_rst
  
! Local variables
!****************************************************************************************
    INTEGER                       :: ierr
    CHARACTER(len=80)             :: abort_message
    CHARACTER(len = 20)           :: modname = 'pbl_surface_init'
    
!****************************************************************************************
! Allocate and initialize module variables with fields read from restart file.
!
!****************************************************************************************    
    ALLOCATE(fder(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('pbl_surface_init', 'pb in allocation',1)

    ALLOCATE(snow(klon,nbsrf), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('pbl_surface_init', 'pb in allocation',1)

    ALLOCATE(qsurf(klon,nbsrf), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('pbl_surface_init', 'pb in allocation',1)

    ALLOCATE(ftsoil(klon,nsoilmx,nbsrf), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('pbl_surface_init', 'pb in allocation',1)

    fder(:)       = fder_rst(:)
    snow(:,:)     = snow_rst(:,:)
    qsurf(:,:)    = qsurf_rst(:,:)
    ftsoil(:,:,:) = ftsoil_rst(:,:,:)

!****************************************************************************************
! Test for sub-surface indices
!
!****************************************************************************************
    IF (is_ter /= 1) THEN 
      WRITE(lunout,*)" *** Warning ***"
      WRITE(lunout,*)" is_ter n'est pas le premier surface, is_ter = ",is_ter
      WRITE(lunout,*)"or on doit commencer par les surfaces continentales"
      abort_message="voir ci-dessus"
      CALL abort_physic(modname,abort_message,1)
    ENDIF

    IF ( is_oce > is_sic ) THEN
      WRITE(lunout,*)' *** Warning ***'
      WRITE(lunout,*)' Pour des raisons de sequencement dans le code'
      WRITE(lunout,*)' l''ocean doit etre traite avant la banquise'
      WRITE(lunout,*)' or is_oce = ',is_oce, '> is_sic = ',is_sic
      abort_message='voir ci-dessus'
      CALL abort_physic(modname,abort_message,1)
    ENDIF

    IF ( is_lic > is_sic ) THEN
      WRITE(lunout,*)' *** Warning ***'
      WRITE(lunout,*)' Pour des raisons de sequencement dans le code'
      WRITE(lunout,*)' la glace contineltalle doit etre traite avant la glace de mer'
      WRITE(lunout,*)' or is_lic = ',is_lic, '> is_sic = ',is_sic
      abort_message='voir ci-dessus'
      CALL abort_physic(modname,abort_message,1)
    ENDIF

!****************************************************************************************
! Validation of ocean mode
!
!****************************************************************************************

    IF (type_ocean /= 'slab  ' .AND. type_ocean /= 'force ' .AND. type_ocean /= 'couple') THEN
       WRITE(lunout,*)' *** Warning ***'
       WRITE(lunout,*)'Option couplage pour l''ocean = ', type_ocean
       abort_message='option pour l''ocean non valable'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    iflag_pbl_surface_t2m_bug=0
    CALL getin_p('iflag_pbl_surface_t2m_bug',iflag_pbl_surface_t2m_bug)
!FC
!    iflag_frein = 0
!    CALL getin_p('iflag_frein',iflag_frein)
!
!jyg<
!****************************************************************************************
! Allocate variables for pbl splitting
!
!****************************************************************************************

    CALL wx_pbl_init
!>jyg

  END SUBROUTINE pbl_surface_init
!  
!****************************************************************************************
!  

  SUBROUTINE pbl_surface( &
       dtime,     date0,     itap,     jour,          &
       debut,     lafin,                              &
       rlon,      rlat,      rugoro,   rmu0,          &
       zsig,      lwdown_m,  pphi,     cldt,          &
       rain_f,    snow_f,    solsw_m,  sollw_m,       &
       gustiness,                                     &
       t,         q,         u,        v,             &
!!! nrlmd+jyg le 02/05/2011 et le 20/02/2012
!!       t_x,       q_x,       t_w,      q_w,           &
       wake_dlt,             wake_dlq,                &
       wake_cstar,           wake_s,                  &
!!!
       pplay,     paprs,     pctsrf,                  &
       ts,SFRWL,   alb_dir, alb_dif,ustar, u10m, v10m,wstar, &
       cdragh,    cdragm,   zu1,    zv1,              &
       alb_dir_m,    alb_dif_m,  zxsens,   zxevap,    &
       alb3_lic,  runoff,    snowhgt,   qsnow,     to_ice,    sissnow,  &
       zxtsol,    zxfluxlat, zt2m,     qsat2m, zn2mout, &
       d_t,       d_q,       d_u,      d_v, d_t_diss, & 
!!! nrlmd+jyg le 02/05/2011 et le 20/02/2012
       d_t_w,     d_q_w,                              &
       d_t_x,     d_q_x,                              & 
!!       d_wake_dlt,d_wake_dlq,                         &
       zxsens_x,  zxfluxlat_x,zxsens_w,zxfluxlat_w,   &
!!!
!!! nrlmd le 13/06/2011
       delta_tsurf,wake_dens,cdragh_x,cdragh_w,       &
       cdragm_x,cdragm_w,kh,kh_x,kh_w,                &
!!!
       zcoefh,    zcoefm,    slab_wfbils,             &
       qsol,    zq2m,      s_pblh,   s_plcl,        &
!!!
!!! jyg le 08/02/2012
       s_pblh_x, s_plcl_x,   s_pblh_w, s_plcl_w,      &
!!!
       s_capCL,   s_oliqCL,  s_cteiCL, s_pblT,        &
       s_therm,   s_trmb1,   s_trmb2,  s_trmb3,       &
       zustar,zu10m,  zv10m,    fder_print,    &
       zxqsurf,   rh2m,      zxfluxu,  zxfluxv,       &
       z0m, z0h,   agesno,  sollw,    solsw,         &
       d_ts,      evap,    fluxlat,  t2m,           &
       wfbils,    wfbilo, wfevap, wfrain, wfsnow,   & 
       flux_t,   flux_u, flux_v,                    &
       dflux_t,   dflux_q,   zxsnow,                  &
!jyg<
!!       zxfluxt,   zxfluxq,   q2m,      flux_q, tke,   &
       zxfluxt,   zxfluxq,   q2m,      flux_q, tke_x,   &
!>jyg
!!! nrlmd+jyg le 02/05/2011 et le 20/02/2012
!!        tke_x,     tke_w                              &
       wake_dltke,                                     &
        treedrg                                   &
!FC
!!!
                        )
!****************************************************************************************
! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
! Objet: interface de "couche limite" (diffusion verticale)
!
!AA REM:
!AA-----
!AA Tout ce qui a trait au traceurs est dans phytrac maintenant
!AA pour l'instant le calcul de la couche limite pour les traceurs
!AA se fait avec cltrac et ne tient pas compte de la differentiation
!AA des sous-fraction de sol.
!AA REM bis :
!AA----------
!AA Pour pouvoir extraire les coefficient d'echanges et le vent 
!AA dans la premiere couche, 3 champs supplementaires ont ete crees
!AA zcoefh, zu1 et zv1. Pour l'instant nous avons moyenne les valeurs
!AA de ces trois champs sur les 4 subsurfaces du modele. Dans l'avenir 
!AA si les informations des subsurfaces doivent etre prises en compte
!AA il faudra sortir ces memes champs en leur ajoutant une dimension, 
!AA c'est a dire nbsrf (nbre de subsurface).
!
! Arguments:
!
! dtime----input-R- interval du temps (secondes)
! itap-----input-I- numero du pas de temps
! date0----input-R- jour initial
! t--------input-R- temperature (K)
! q--------input-R- vapeur d'eau (kg/kg)
! u--------input-R- vitesse u
! v--------input-R- vitesse v
! wake_dlt-input-R- temperatre difference between (w) and (x) (K)
! wake_dlq-input-R- humidity difference between (w) and (x) (kg/kg)
!wake_cstar-input-R- wake gust front speed (m/s)
! wake_s---input-R- wake fractionnal area
! ts-------input-R- temperature du sol (en Kelvin)
! paprs----input-R- pression a intercouche (Pa)
! pplay----input-R- pression au milieu de couche (Pa)
! rlat-----input-R- latitude en degree
! z0m, z0h ----input-R- longeur de rugosite (en m)
! Martin
! zsig-----input-R- slope
! cldt-----input-R- total cloud fraction
! pphi-----input-R- geopotentiel de chaque couche (g z) (reference sol)
! Martin
!
! d_t------output-R- le changement pour "t"
! d_q------output-R- le changement pour "q"
! d_u------output-R- le changement pour "u"
! d_v------output-R- le changement pour "v"
! d_ts-----output-R- le changement pour "ts"
! flux_t---output-R- flux de chaleur sensible (CpT) J/m**2/s (W/m**2)
!                    (orientation positive vers le bas)
! tke_x---input/output-R- tke in the (x) region (kg/m**2/s)
! wake_dltke-input/output-R- tke difference between (w) and (x) (kg/m**2/s)
! flux_q---output-R- flux de vapeur d'eau (kg/m**2/s)
! flux_u---output-R- tension du vent X: (kg m/s)/(m**2 s) ou Pascal
! flux_v---output-R- tension du vent Y: (kg m/s)/(m**2 s) ou Pascal
! dflux_t--output-R- derive du flux sensible
! dflux_q--output-R- derive du flux latent
! zu1------output-R- le vent dans la premiere couche
! zv1------output-R- le vent dans la premiere couche
! trmb1----output-R- deep_cape
! trmb2----output-R- inhibition 
! trmb3----output-R- Point Omega
! cteiCL---output-R- Critere d'instab d'entrainmt des nuages de CL
! plcl-----output-R- Niveau de condensation
! pblh-----output-R- HCL
! pblT-----output-R- T au nveau HCL
! treedrg--output-R- tree drag (m)               
!
    USE carbon_cycle_mod,   ONLY : carbon_cycle_cpl, carbon_cycle_tr, level_coupling_esm 
    USE carbon_cycle_mod,   ONLY : co2_send, nbcf_out, fields_out, yfields_out, cfname_out
    USE indice_sol_mod
    USE time_phylmdz_mod,   ONLY : day_ini,annee_ref,itau_phy
    USE mod_grid_phy_lmdz,  ONLY : nbp_lon, nbp_lat, grid1dto2d_glo
    USE print_control_mod,  ONLY : prt_level,lunout
    USE ioipsl_getin_p_mod, ONLY : getin_p

    IMPLICIT NONE

    INCLUDE "dimsoil.h"
    INCLUDE "YOMCST.h"
    INCLUDE "YOETHF.h"
    INCLUDE "FCTTRE.h"
    INCLUDE "clesphys.h"
    INCLUDE "compbl.h"
    INCLUDE "flux_arp.h"
!FC 
    INCLUDE "dimpft.h"

!****************************************************************************************
    REAL,                         INTENT(IN)        :: dtime   ! time interval (s)
    REAL,                         INTENT(IN)        :: date0   ! initial day
    INTEGER,                      INTENT(IN)        :: itap    ! time step
    INTEGER,                      INTENT(IN)        :: jour    ! current day of the year
    LOGICAL,                      INTENT(IN)        :: debut   ! true if first run step
    LOGICAL,                      INTENT(IN)        :: lafin   ! true if last run step
    REAL, DIMENSION(klon),        INTENT(IN)        :: rlon    ! longitudes in degrees
    REAL, DIMENSION(klon),        INTENT(IN)        :: rlat    ! latitudes in degrees
    REAL, DIMENSION(klon),        INTENT(IN)        :: rugoro  ! rugosity length
    REAL, DIMENSION(klon),        INTENT(IN)        :: rmu0    ! cosine of solar zenith angle
    REAL, DIMENSION(klon),        INTENT(IN)        :: rain_f  ! rain fall
    REAL, DIMENSION(klon),        INTENT(IN)        :: snow_f  ! snow fall
    REAL, DIMENSION(klon),        INTENT(IN)        :: solsw_m ! net shortwave radiation at mean surface
    REAL, DIMENSION(klon),        INTENT(IN)        :: sollw_m ! net longwave radiation at mean surface
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: t       ! temperature (K)
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: q       ! water vapour (kg/kg)
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: u       ! u speed
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: v       ! v speed
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: pplay   ! mid-layer pression (Pa)
    REAL, DIMENSION(klon,klev+1), INTENT(IN)        :: paprs   ! pression between layers (Pa) 
    REAL, DIMENSION(klon, nbsrf), INTENT(IN)        :: pctsrf  ! sub-surface fraction
! Martin
    REAL, DIMENSION(klon),        INTENT(IN)        :: zsig    ! slope
    REAL, DIMENSION(klon),        INTENT(IN)        :: lwdown_m ! downward longwave radiation at mean s    
    REAL, DIMENSION(klon),        INTENT(IN)        :: gustiness ! gustiness

    REAL, DIMENSION(klon),        INTENT(IN)        :: cldt    ! total cloud fraction
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: pphi    ! geopotential (m2/s2)
! Martin

!!! nrlmd+jyg le 02/05/2011 et le 20/02/2012
!!    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: t_x       ! Temp\'erature hors poche froide
!!    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: t_w       ! Temp\'erature dans la poches froide
!!    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: q_x       ! 
!!    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: q_w       ! Pareil pour l'humidit\'e
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: wake_dlt  !temperature difference between (w) and (x) (K)
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: wake_dlq  !humidity difference between (w) and (x) (K)
    REAL, DIMENSION(klon),        INTENT(IN)        :: wake_s    ! Fraction de poches froides
    REAL, DIMENSION(klon),        INTENT(IN)        :: wake_cstar! Vitesse d'expansion des poches froides
    REAL, DIMENSION(klon),        INTENT(IN)        :: wake_dens
!!!

! Input/Output variables
!****************************************************************************************
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)     :: ts      ! temperature at surface (K)
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)     :: delta_tsurf !surface temperature difference between
                                                                   !wake and off-wake regions
!albedo SB >>>
    REAL, DIMENSIOn(6),intent(in) :: SFRWL
    REAL, DIMENSION(klon, nsw, nbsrf), INTENT(INOUT)     :: alb_dir,alb_dif
!albedo SB <<<
!jyg Pourquoi ustar et wstar sont-elles INOUT ?
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)     :: ustar   ! u* (m/s)
    REAL, DIMENSION(klon, nbsrf+1), INTENT(INOUT)   :: wstar   ! w* (m/s)
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)     :: u10m    ! u speed at 10m
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)     :: v10m    ! v speed at 10m
!jyg<
!!    REAL, DIMENSION(klon, klev+1, nbsrf+1), INTENT(INOUT) :: tke
    REAL, DIMENSION(klon, klev+1, nbsrf+1), INTENT(INOUT) :: tke_x
!>jyg 

!!! nrlmd+jyg le 02/05/2011 et le 20/02/2012
    REAL, DIMENSION(klon, klev+1, nbsrf+1), INTENT(INOUT) :: wake_dltke ! TKE_w - TKE_x
!!!

! Output variables
!****************************************************************************************
    REAL, DIMENSION(klon),        INTENT(OUT)       :: cdragh     ! drag coefficient for T and Q
    REAL, DIMENSION(klon),        INTENT(OUT)       :: cdragm     ! drag coefficient for wind
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zu1        ! u wind speed in first layer
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zv1        ! v wind speed in first layer
!albedo SB >>>
    REAL, DIMENSION(klon, nsw),   INTENT(OUT)       :: alb_dir_m,alb_dif_m
!albedo SB <<<
    ! Martin
    REAL, DIMENSION(klon),        INTENT(OUT)       :: alb3_lic
    ! Martin
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxsens     ! sensible heat flux at surface with inversed sign 
                                                                  ! (=> positive sign upwards)
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxevap     ! water vapour flux at surface, positiv upwards
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxtsol     ! temperature at surface, mean for each grid point
!!! jyg le ???
    REAL, DIMENSION(klon,klev),   INTENT(OUT)       :: d_t_w      !   !
    REAL, DIMENSION(klon,klev),   INTENT(OUT)       :: d_q_w      !      !  Tendances dans les poches
    REAL, DIMENSION(klon,klev),   INTENT(OUT)       :: d_t_x      !   !
    REAL, DIMENSION(klon,klev),   INTENT(OUT)       :: d_q_x      !      !  Tendances hors des poches
!!! jyg
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxfluxlat  ! latent flux, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zt2m       ! temperature at 2m, mean for each grid point
    INTEGER, DIMENSION(klon, 6),  INTENT(OUT)       :: zn2mout ! nb of times t2m is out of the [tsol,temp]
    REAL, DIMENSION(klon),        INTENT(OUT)       :: qsat2m
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_t        ! change in temperature 
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_t_diss       ! change in temperature 
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_q        ! change in water vapour
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_u        ! change in u speed
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_v        ! change in v speed

    REAL, INTENT(OUT):: zcoefh(:, :, :) ! (klon, klev, nbsrf + 1)
    ! coef for turbulent diffusion of T and Q, mean for each grid point

    REAL, INTENT(OUT):: zcoefm(:, :, :) ! (klon, klev, nbsrf + 1)
    ! coef for turbulent diffusion of U and V (?), mean for each grid point

!!! nrlmd+jyg le 02/05/2011 et le 20/02/2012
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxsens_x   ! Flux sensible hors poche
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxsens_w   ! Flux sensible dans la poche
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxfluxlat_x! Flux latent hors poche
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxfluxlat_w! Flux latent dans la poche
!!    REAL, DIMENSION(klon,klev),   INTENT(OUT)       :: d_wake_dlt
!!    REAL, DIMENSION(klon,klev),   INTENT(OUT)       :: d_wake_dlq

! Output only for diagnostics
    REAL, DIMENSION(klon),        INTENT(OUT)       :: cdragh_x
    REAL, DIMENSION(klon),        INTENT(OUT)       :: cdragh_w
    REAL, DIMENSION(klon),        INTENT(OUT)       :: cdragm_x
    REAL, DIMENSION(klon),        INTENT(OUT)       :: cdragm_w
    REAL, DIMENSION(klon),        INTENT(OUT)       :: kh
    REAL, DIMENSION(klon),        INTENT(OUT)       :: kh_x
    REAL, DIMENSION(klon),        INTENT(OUT)       :: kh_w
!!! 
    REAL, DIMENSION(klon),        INTENT(OUT)       :: slab_wfbils! heat balance at surface only for slab at ocean points
    REAL, DIMENSION(klon),        INTENT(OUT)       :: qsol     ! water height in the soil (mm)
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zq2m       ! water vapour at 2m, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_pblh     ! height of the planetary boundary layer(HPBL)
!!! jyg le 08/02/2012
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_pblh_x   ! height of the PBL in the off-wake region
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_pblh_w   ! height of the PBL in the wake region
!!!
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_plcl     ! condensation level
!!! jyg le 08/02/2012
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_plcl_x   ! condensation level in the off-wake region
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_plcl_w   ! condensation level in the wake region
!!!
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_capCL    ! CAPE of PBL
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_oliqCL   ! liquid water intergral of PBL
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_cteiCL   ! cloud top instab. crit. of PBL
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_pblT     ! temperature at PBLH
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_therm    ! thermal virtual temperature excess
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_trmb1    ! deep cape, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_trmb2    ! inhibition, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_trmb3    ! point Omega, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zustar     ! u*
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zu10m      ! u speed at 10m, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zv10m      ! v speed at 10m, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: fder_print ! fder for printing (=fder(i) + dflux_t(i) + dflux_q(i))
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxqsurf    ! humidity at surface, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: rh2m       ! relative humidity at 2m
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: zxfluxu    ! u wind tension, mean for each grid point
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: zxfluxv    ! v wind tension, mean for each grid point
    REAL, DIMENSION(klon, nbsrf+1), INTENT(INOUT)     :: z0m,z0h      ! rugosity length (m)
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)       :: agesno   ! age of snow at surface
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: solsw      ! net shortwave radiation at surface 
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: sollw      ! net longwave radiation at surface
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: d_ts       ! change in temperature at surface
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)       :: evap     ! evaporation at surface
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: fluxlat    ! latent flux
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: t2m        ! temperature at 2 meter height
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: wfbils     ! heat balance at surface
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: wfbilo     ! water balance at surface weighted by srf
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: wfevap     ! water balance (evap) at surface weighted by srf
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: wfrain     ! water balance (rain) at surface weighted by srf
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: wfsnow     ! water balance (snow) at surface weighted by srf
    REAL, DIMENSION(klon, klev, nbsrf), INTENT(OUT) :: flux_t     ! sensible heat flux (CpT) J/m**2/s (W/m**2)
                                                                  ! positve orientation downwards
    REAL, DIMENSION(klon, klev, nbsrf), INTENT(OUT) :: flux_u     ! u wind tension (kg m/s)/(m**2 s) or Pascal
    REAL, DIMENSION(klon, klev, nbsrf), INTENT(OUT) :: flux_v     ! v wind tension (kg m/s)/(m**2 s) or Pascal
!FC
    REAL, DIMENSION(klon, klev, nbsrf), INTENT(INOUT)  :: treedrg      ! tree drag (m)               


! Output not needed
    REAL, DIMENSION(klon),       INTENT(OUT)        :: dflux_t    ! change of sensible heat flux 
    REAL, DIMENSION(klon),       INTENT(OUT)        :: dflux_q    ! change of water vapour flux
    REAL, DIMENSION(klon),       INTENT(OUT)        :: zxsnow     ! snow at surface, mean for each grid point
    REAL, DIMENSION(klon, klev), INTENT(OUT)        :: zxfluxt    ! sensible heat flux, mean for each grid point
    REAL, DIMENSION(klon, klev), INTENT(OUT)        :: zxfluxq    ! water vapour flux, mean for each grid point
    REAL, DIMENSION(klon, nbsrf),INTENT(OUT)        :: q2m        ! water vapour at 2 meter height
    REAL, DIMENSION(klon, klev, nbsrf), INTENT(OUT) :: flux_q     ! water vapour flux(latent flux) (kg/m**2/s)

! Martin
! sisvat
    REAL, DIMENSION(klon),       INTENT(OUT)        :: qsnow      ! snow water content
    REAL, DIMENSION(klon),       INTENT(OUT)        :: snowhgt    ! snow height
    REAL, DIMENSION(klon),       INTENT(OUT)        :: to_ice     ! snow passed to ice
    REAL, DIMENSION(klon),       INTENT(OUT)        :: sissnow    ! snow in snow model
    REAL, DIMENSION(klon),       INTENT(OUT)        :: runoff     ! runoff on land ice
! Martin

! Local variables with attribute SAVE
!****************************************************************************************
    INTEGER, SAVE                            :: nhoridbg, nidbg   ! variables for IOIPSL
!$OMP THREADPRIVATE(nhoridbg, nidbg)
    LOGICAL, SAVE                            :: debugindex=.FALSE.
!$OMP THREADPRIVATE(debugindex)
    LOGICAL, SAVE                            :: first_call=.TRUE.
!$OMP THREADPRIVATE(first_call)
    CHARACTER(len=8), DIMENSION(nbsrf), SAVE :: cl_surf
!$OMP THREADPRIVATE(cl_surf)
    REAL, SAVE                               :: beta_land         ! beta for wx_dts
!$OMP THREADPRIVATE(beta_land)

! Other local variables
!****************************************************************************************
! >> PC
    INTEGER                            :: ierr
    INTEGER                            :: n 
! << PC
    INTEGER                            :: iflag_split 
    INTEGER                            :: i, k, nsrf 
    INTEGER                            :: knon, j
    INTEGER                            :: idayref
    INTEGER , DIMENSION(klon)          :: ni
    REAL                               :: yt1_new
    REAL                               :: zx_alf1, zx_alf2 !valeur ambiante par extrapola
    REAL                               :: amn, amx
    REAL                               :: f1 ! fraction de longeurs visibles parmi tout SW intervalle
    REAL, DIMENSION(klon)              :: r_co2_ppm     ! taux CO2 atmosphere
    REAL, DIMENSION(klon)              :: yts, yz0m, yz0h, ypct
!albedo SB >>>
    REAL, DIMENSION(klon)              :: yalb,yalb_vis
!albedo SB <<<
    REAL, DIMENSION(klon)              :: yt1, yq1, yu1, yv1
    REAL, DIMENSION(klon)              :: ysnow, yqsurf, yagesno, yqsol
    REAL, DIMENSION(klon)              :: yrain_f, ysnow_f
    REAL, DIMENSION(klon)              :: ysolsw, ysollw
    REAL, DIMENSION(klon)              :: yfder
    REAL, DIMENSION(klon)              :: yrugoro
    REAL, DIMENSION(klon)              :: yfluxlat
    REAL, DIMENSION(klon)              :: y_d_ts
    REAL, DIMENSION(klon)              :: y_flux_t1, y_flux_q1
    REAL, DIMENSION(klon)              :: y_dflux_t, y_dflux_q
    REAL, DIMENSION(klon)              :: y_flux_u1, y_flux_v1
    REAL, DIMENSION(klon)              :: yt2m, yq2m, yu10m
    INTEGER, DIMENSION(klon, nbsrf, 6) :: yn2mout, yn2mout_x, yn2mout_w
    INTEGER, DIMENSION(klon, nbsrf, 6) :: n2mout, n2mout_x, n2mout_w
    REAL, DIMENSION(klon)              :: yustar
    REAL, DIMENSION(klon)              :: ywstar
    REAL, DIMENSION(klon)              :: ywindsp
    REAL, DIMENSION(klon)              :: yt10m, yq10m
    REAL, DIMENSION(klon)              :: ypblh
    REAL, DIMENSION(klon)              :: ylcl
    REAL, DIMENSION(klon)              :: ycapCL
    REAL, DIMENSION(klon)              :: yoliqCL
    REAL, DIMENSION(klon)              :: ycteiCL
    REAL, DIMENSION(klon)              :: ypblT
    REAL, DIMENSION(klon)              :: ytherm
    REAL, DIMENSION(klon)              :: ytrmb1
    REAL, DIMENSION(klon)              :: ytrmb2
    REAL, DIMENSION(klon)              :: ytrmb3
    REAL, DIMENSION(klon)              :: uzon, vmer
    REAL, DIMENSION(klon)              :: tair1, qair1, tairsol
    REAL, DIMENSION(klon)              :: psfce, patm
    REAL, DIMENSION(klon)              :: qairsol, zgeo1, speed, zri1, pref !speed, zri1, pref, added by Fuxing WANG, 04/03/2015
    REAL, DIMENSION(klon)              :: yz0h_oupas
    REAL, DIMENSION(klon)              :: yfluxsens
    REAL, DIMENSION(klon)              :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon)              :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon)              :: ypsref
    REAL, DIMENSION(klon)              :: yevap, ytsurf_new, yalb3_new
!albedo SB >>>
    REAL, DIMENSION(klon,nsw)          :: yalb_dir_new, yalb_dif_new
!albedo SB <<<
    REAL, DIMENSION(klon)              :: ztsol
    REAL, DIMENSION(klon)              :: meansqT ! mean square deviation of subsurface temperatures
    REAL, DIMENSION(klon)              :: alb_m  ! mean albedo for whole SW interval
    REAL, DIMENSION(klon,klev)         :: y_d_t, y_d_q, y_d_t_diss
    REAL, DIMENSION(klon,klev)         :: y_d_u, y_d_v
    REAL, DIMENSION(klon,klev)         :: y_flux_t, y_flux_q
    REAL, DIMENSION(klon,klev)         :: y_flux_u, y_flux_v
    REAL, DIMENSION(klon,klev)         :: ycoefh, ycoefm,ycoefq
    REAL, DIMENSION(klon)              :: ycdragh, ycdragm
    REAL, DIMENSION(klon,klev)         :: yu, yv
    REAL, DIMENSION(klon,klev)         :: yt, yq
    REAL, DIMENSION(klon,klev)         :: ypplay, ydelp
    REAL, DIMENSION(klon,klev)         :: delp
    REAL, DIMENSION(klon,klev+1)       :: ypaprs
    REAL, DIMENSION(klon,klev+1)       :: ytke
    REAL, DIMENSION(klon,nsoilmx)      :: ytsoil
!FC 
    REAL, DIMENSION(klon,nvm_lmdz)          :: yveget
    REAL, DIMENSION(klon,nvm_lmdz)          :: ylai
    REAL, DIMENSION(klon,nvm_lmdz)          :: yheight
    REAL, DIMENSION(klon,klev)              :: y_d_u_frein
    REAL, DIMENSION(klon,klev)              :: y_d_v_frein
    REAL, DIMENSION(klon,klev)              :: y_treedrg
!FC

    CHARACTER(len=80)                  :: abort_message
    CHARACTER(len=20)                  :: modname = 'pbl_surface'
    LOGICAL, PARAMETER                 :: zxli=.FALSE. ! utiliser un jeu de fonctions simples
    LOGICAL, PARAMETER                 :: check=.FALSE.

!!! nrlmd le 02/05/2011
!!! jyg le 07/02/2012
    REAL, DIMENSION(klon)              :: ywake_s, ywake_cstar, ywake_dens
!!!
    REAL, DIMENSION(klon,klev+1)       :: ytke_x, ytke_w
    REAL, DIMENSION(klon,klev+1)       :: ywake_dltke
    REAL, DIMENSION(klon,klev)         :: yu_x, yv_x, yu_w, yv_w
    REAL, DIMENSION(klon,klev)         :: yt_x, yq_x, yt_w, yq_w
    REAL, DIMENSION(klon,klev)         :: ycoefh_x, ycoefm_x, ycoefh_w, ycoefm_w
    REAL, DIMENSION(klon,klev)         :: ycoefq_x, ycoefq_w
    REAL, DIMENSION(klon)              :: ycdragh_x, ycdragm_x, ycdragh_w, ycdragm_w
    REAL, DIMENSION(klon)              :: AcoefH_x, AcoefQ_x, BcoefH_x, BcoefQ_x
    REAL, DIMENSION(klon)              :: AcoefH_w, AcoefQ_w, BcoefH_w, BcoefQ_w
    REAL, DIMENSION(klon)              :: AcoefU_x, AcoefV_x, BcoefU_x, BcoefV_x
    REAL, DIMENSION(klon)              :: AcoefU_w, AcoefV_w, BcoefU_w, BcoefV_w
    REAL, DIMENSION(klon)              :: y_flux_t1_x, y_flux_q1_x, y_flux_t1_w, y_flux_q1_w
    REAL, DIMENSION(klon)              :: y_flux_u1_x, y_flux_v1_x, y_flux_u1_w, y_flux_v1_w
    REAL, DIMENSION(klon,klev)         :: y_flux_t_x, y_flux_q_x, y_flux_t_w, y_flux_q_w
    REAL, DIMENSION(klon,klev)         :: y_flux_u_x, y_flux_v_x, y_flux_u_w, y_flux_v_w
    REAL, DIMENSION(klon)              :: yfluxlat_x, yfluxlat_w
    REAL, DIMENSION(klon,klev)         :: y_d_t_x, y_d_q_x, y_d_t_w, y_d_q_w
    REAL, DIMENSION(klon,klev)         :: y_d_t_diss_x, y_d_t_diss_w
    REAL, DIMENSION(klon,klev)         :: d_t_diss_x, d_t_diss_w
    REAL, DIMENSION(klon,klev)         :: y_d_u_x, y_d_v_x, y_d_u_w, y_d_v_w
    REAL, DIMENSION(klon, klev, nbsrf) :: flux_t_x, flux_q_x, flux_t_w, flux_q_w
    REAL, DIMENSION(klon, klev, nbsrf) :: flux_u_x, flux_v_x, flux_u_w, flux_v_w
    REAL, DIMENSION(klon, nbsrf)       :: fluxlat_x, fluxlat_w
    REAL, DIMENSION(klon, klev)        :: zxfluxt_x, zxfluxq_x, zxfluxt_w, zxfluxq_w
    REAL, DIMENSION(klon, klev)        :: zxfluxu_x, zxfluxv_x, zxfluxu_w, zxfluxv_w
    REAL                               :: zx_qs_surf, zcor_surf, zdelta_surf
    REAL, DIMENSION(klon)              :: ytsurf_th, yqsatsurf
    REAL, DIMENSION(klon)              :: ybeta
    REAL, DIMENSION(klon, klev)        :: d_u_x
    REAL, DIMENSION(klon, klev)        :: d_u_w
    REAL, DIMENSION(klon, klev)        :: d_v_x
    REAL, DIMENSION(klon, klev)        :: d_v_w 

    REAL, DIMENSION(klon,klev)         :: CcoefH, CcoefQ, DcoefH, DcoefQ
    REAL, DIMENSION(klon,klev)         :: CcoefU, CcoefV, DcoefU, DcoefV
    REAL, DIMENSION(klon,klev)         :: CcoefH_x, CcoefQ_x, DcoefH_x, DcoefQ_x
    REAL, DIMENSION(klon,klev)         :: CcoefH_w, CcoefQ_w, DcoefH_w, DcoefQ_w
    REAL, DIMENSION(klon,klev)         :: CcoefU_x, CcoefV_x, DcoefU_x, DcoefV_x
    REAL, DIMENSION(klon,klev)         :: CcoefU_w, CcoefV_w, DcoefU_w, DcoefV_w
    REAL, DIMENSION(klon,klev)         :: Kcoef_hq, Kcoef_m, gama_h, gama_q
    REAL, DIMENSION(klon,klev)         :: Kcoef_hq_x, Kcoef_m_x, gama_h_x, gama_q_x
    REAL, DIMENSION(klon,klev)         :: Kcoef_hq_w, Kcoef_m_w, gama_h_w, gama_q_w
    REAL, DIMENSION(klon)              :: alf_1, alf_2, alf_1_x, alf_2_x, alf_1_w, alf_2_w
!!!
!!!jyg le 08/02/2012
    REAL, DIMENSION(klon, nbsrf)       :: windsp
!
    REAL, DIMENSION(klon, nbsrf)       :: t2m_x
    REAL, DIMENSION(klon, nbsrf)       :: q2m_x
    REAL, DIMENSION(klon)              :: rh2m_x
    REAL, DIMENSION(klon)              :: qsat2m_x
    REAL, DIMENSION(klon, nbsrf)       :: u10m_x
    REAL, DIMENSION(klon, nbsrf)       :: v10m_x
    REAL, DIMENSION(klon, nbsrf)       :: ustar_x
    REAL, DIMENSION(klon, nbsrf)       :: wstar_x
!              
    REAL, DIMENSION(klon, nbsrf)       :: pblh_x
    REAL, DIMENSION(klon, nbsrf)       :: plcl_x
    REAL, DIMENSION(klon, nbsrf)       :: capCL_x
    REAL, DIMENSION(klon, nbsrf)       :: oliqCL_x
    REAL, DIMENSION(klon, nbsrf)       :: cteiCL_x
    REAL, DIMENSION(klon, nbsrf)       :: pblt_x
    REAL, DIMENSION(klon, nbsrf)       :: therm_x
    REAL, DIMENSION(klon, nbsrf)       :: trmb1_x
    REAL, DIMENSION(klon, nbsrf)       :: trmb2_x
    REAL, DIMENSION(klon, nbsrf)       :: trmb3_x
!
    REAL, DIMENSION(klon, nbsrf)       :: t2m_w
    REAL, DIMENSION(klon, nbsrf)       :: q2m_w
    REAL, DIMENSION(klon)              :: rh2m_w
    REAL, DIMENSION(klon)              :: qsat2m_w
    REAL, DIMENSION(klon, nbsrf)       :: u10m_w
    REAL, DIMENSION(klon, nbsrf)       :: v10m_w
    REAL, DIMENSION(klon, nbsrf)       :: ustar_w
    REAL, DIMENSION(klon, nbsrf)       :: wstar_w
!                           
    REAL, DIMENSION(klon, nbsrf)       :: pblh_w
    REAL, DIMENSION(klon, nbsrf)       :: plcl_w
    REAL, DIMENSION(klon, nbsrf)       :: capCL_w
    REAL, DIMENSION(klon, nbsrf)       :: oliqCL_w
    REAL, DIMENSION(klon, nbsrf)       :: cteiCL_w
    REAL, DIMENSION(klon, nbsrf)       :: pblt_w
    REAL, DIMENSION(klon, nbsrf)       :: therm_w
    REAL, DIMENSION(klon, nbsrf)       :: trmb1_w
    REAL, DIMENSION(klon, nbsrf)       :: trmb2_w
    REAL, DIMENSION(klon, nbsrf)       :: trmb3_w
!
    REAL, DIMENSION(klon)       :: yt2m_x
    REAL, DIMENSION(klon)       :: yq2m_x
    REAL, DIMENSION(klon)       :: yt10m_x
    REAL, DIMENSION(klon)       :: yq10m_x
    REAL, DIMENSION(klon)       :: yu10m_x
    REAL, DIMENSION(klon)       :: yv10m_x
    REAL, DIMENSION(klon)       :: yustar_x
    REAL, DIMENSION(klon)       :: ywstar_x
!              
    REAL, DIMENSION(klon)       :: ypblh_x
    REAL, DIMENSION(klon)       :: ylcl_x
    REAL, DIMENSION(klon)       :: ycapCL_x
    REAL, DIMENSION(klon)       :: yoliqCL_x
    REAL, DIMENSION(klon)       :: ycteiCL_x
    REAL, DIMENSION(klon)       :: ypblt_x
    REAL, DIMENSION(klon)       :: ytherm_x
    REAL, DIMENSION(klon)       :: ytrmb1_x
    REAL, DIMENSION(klon)       :: ytrmb2_x
    REAL, DIMENSION(klon)       :: ytrmb3_x
!
    REAL, DIMENSION(klon)       :: yt2m_w
    REAL, DIMENSION(klon)       :: yq2m_w
    REAL, DIMENSION(klon)       :: yt10m_w
    REAL, DIMENSION(klon)       :: yq10m_w
    REAL, DIMENSION(klon)       :: yu10m_w
    REAL, DIMENSION(klon)       :: yv10m_w
    REAL, DIMENSION(klon)       :: yustar_w
    REAL, DIMENSION(klon)       :: ywstar_w
!                       
    REAL, DIMENSION(klon)       :: ypblh_w
    REAL, DIMENSION(klon)       :: ylcl_w
    REAL, DIMENSION(klon)       :: ycapCL_w
    REAL, DIMENSION(klon)       :: yoliqCL_w
    REAL, DIMENSION(klon)       :: ycteiCL_w
    REAL, DIMENSION(klon)       :: ypblt_w
    REAL, DIMENSION(klon)       :: ytherm_w
    REAL, DIMENSION(klon)       :: ytrmb1_w
    REAL, DIMENSION(klon)       :: ytrmb2_w
    REAL, DIMENSION(klon)       :: ytrmb3_w
!
    REAL, DIMENSION(klon)       :: uzon_x, vmer_x, speed_x, zri1_x, pref_x !speed_x, zri1_x, pref_x, added by Fuxing WANG, 04/03/2015
    REAL, DIMENSION(klon)       :: zgeo1_x, tair1_x, qair1_x, tairsol_x
!
    REAL, DIMENSION(klon)       :: uzon_w, vmer_w, speed_w, zri1_w, pref_w !speed_w, zri1_w, pref_w, added by Fuxing WANG, 04/03/2015
    REAL, DIMENSION(klon)       :: zgeo1_w, tair1_w, qair1_w, tairsol_w

!!! jyg le 25/03/2013
!!    Variables intermediaires pour le raccord des deux colonnes \`a la surface
!jyg<
!!    REAL   ::   dd_Ch
!!    REAL   ::   dd_Cm
!!    REAL   ::   dd_Kh
!!    REAL   ::   dd_Km
!!    REAL   ::   dd_u 
!!    REAL   ::   dd_v 
!!    REAL   ::   dd_t 
!!    REAL   ::   dd_q 
!!    REAL   ::   dd_AH
!!    REAL   ::   dd_AQ
!!    REAL   ::   dd_AU
!!    REAL   ::   dd_AV
!!    REAL   ::   dd_BH
!!    REAL   ::   dd_BQ
!!    REAL   ::   dd_BU
!!    REAL   ::   dd_BV
!!
!!    REAL   ::   dd_KHp
!!    REAL   ::   dd_KQp
!!    REAL   ::   dd_KUp
!!    REAL   ::   dd_KVp
!>jyg

!!!
!!! nrlmd le 13/06/2011
    REAL, DIMENSION(klon)              :: y_delta_flux_t1, y_delta_flux_q1, y_delta_flux_u1, y_delta_flux_v1
    REAL, DIMENSION(klon)              :: y_delta_tsurf,delta_coef,tau_eq
    REAL, PARAMETER                    :: facteur=2./sqrt(3.14)
    REAL, PARAMETER                    :: inertia=2000.
    REAL, DIMENSION(klon)              :: ytsurf_th_x,ytsurf_th_w,yqsatsurf_x,yqsatsurf_w
    REAL, DIMENSION(klon)              :: ydtsurf_th
    REAL                               :: zdelta_surf_x,zdelta_surf_w,zx_qs_surf_x,zx_qs_surf_w
    REAL                               :: zcor_surf_x,zcor_surf_w
    REAL                               :: mod_wind_x, mod_wind_w
    REAL                               :: rho1
    REAL, DIMENSION(klon)              :: Kech_h           ! Coefficient d'echange pour l'energie
    REAL, DIMENSION(klon)              :: Kech_h_x, Kech_h_w 
    REAL, DIMENSION(klon)              :: Kech_m
    REAL, DIMENSION(klon)              :: Kech_m_x, Kech_m_w 
    REAL, DIMENSION(klon)              :: yts_x,yts_w
!jyg<
!!    REAL, DIMENSION(klon)              :: Kech_Hp, Kech_H_xp, Kech_H_wp
!!    REAL, DIMENSION(klon)              :: Kech_Qp, Kech_Q_xp, Kech_Q_wp
!!    REAL, DIMENSION(klon)              :: Kech_Up, Kech_U_xp, Kech_U_wp
!!    REAL, DIMENSION(klon)              :: Kech_Vp, Kech_V_xp, Kech_V_wp
!>jyg
!jyg<
    REAL, DIMENSION(klon)              :: ah, bh     ! coefficients of the delta_Tsurf equation
!>jyg

    REAL                               :: vent
!
! For debugging with IOIPSL
    INTEGER, DIMENSION(nbp_lon*nbp_lat)    :: ndexbg
    REAL                               :: zjulian
    REAL, DIMENSION(klon)              :: tabindx
    REAL, DIMENSION(nbp_lon,nbp_lat)         :: zx_lon, zx_lat
    REAL, DIMENSION(nbp_lon,nbp_lat)         :: debugtab


    REAL, DIMENSION(klon,nbsrf)        :: pblh         ! height of the planetary boundary layer
    REAL, DIMENSION(klon,nbsrf)        :: plcl         ! condensation level
    REAL, DIMENSION(klon,nbsrf)        :: capCL
    REAL, DIMENSION(klon,nbsrf)        :: oliqCL
    REAL, DIMENSION(klon,nbsrf)        :: cteiCL
    REAL, DIMENSION(klon,nbsrf)        :: pblT
    REAL, DIMENSION(klon,nbsrf)        :: therm
    REAL, DIMENSION(klon,nbsrf)        :: trmb1        ! deep cape
    REAL, DIMENSION(klon,nbsrf)        :: trmb2        ! inhibition
    REAL, DIMENSION(klon,nbsrf)        :: trmb3        ! point Omega
    REAL, DIMENSION(klon,nbsrf)        :: zx_rh2m, zx_qsat2m
    REAL, DIMENSION(klon,nbsrf)        :: zx_t1
    REAL, DIMENSION(klon, nbsrf)       :: alb          ! mean albedo for whole SW interval
    REAL, DIMENSION(klon)              :: ylwdown      ! jg : temporary (ysollwdown)
    REAL, DIMENSION(klon)              :: ygustiness      ! jg : temporary (ysollwdown)

    REAL                               :: zx_qs1, zcor1, zdelta1 

    ! Martin
    REAL, DIMENSION(klon, nbsrf)       :: sollwd ! net longwave radiation at surface
    REAL, DIMENSION(klon)              :: ytoice
    REAL, DIMENSION(klon)              :: ysnowhgt, yqsnow, ysissnow, yrunoff
    REAL, DIMENSION(klon)              :: yzsig
    REAL, DIMENSION(klon,klev)         :: ypphi
    REAL, DIMENSION(klon)              :: ycldt
    REAL, DIMENSION(klon)              :: yrmu0
    ! Martin

!****************************************************************************************
! End of declarations
!****************************************************************************************

      IF (prt_level >=10) print *,' -> pbl_surface, itap ',itap
!
!!jyg      iflag_split = mod(iflag_pbl_split,2)
      iflag_split = mod(iflag_pbl_split,10)

!****************************************************************************************
! 1) Initialisation and validation tests 
!    Only done first time entering this subroutine
!
!****************************************************************************************

    IF (first_call) THEN

       iflag_new_t2mq2m=1
       CALL getin_p('iflag_new_t2mq2m',iflag_new_t2mq2m)
       print*,'pbl_iflag_new_t2mq2m=',iflag_new_t2mq2m

       print*,'PBL SURFACE AVEC GUSTINESS'
       first_call=.FALSE.
      
       ! Initialize ok_flux_surf (for 1D model)
       if (klon_glo>1) ok_flux_surf=.FALSE.

       ! intialize beta_land
       beta_land = 0.5
       call getin_p('beta_land', beta_land)
       
       ! Initilize debug IO
       IF (debugindex .AND. mpi_size==1) THEN 
          ! initialize IOIPSL output
          idayref = day_ini
          CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
          CALL grid1dTo2d_glo(rlon,zx_lon)
          DO i = 1, nbp_lon
             zx_lon(i,1) = rlon(i+1)
             zx_lon(i,nbp_lat) = rlon(i+1)
          ENDDO
          CALL grid1dTo2d_glo(rlat,zx_lat)
          CALL histbeg("sous_index",nbp_lon,zx_lon(:,1),nbp_lat,zx_lat(1,:), &
               1,nbp_lon,1,nbp_lat, &
               itau_phy,zjulian,dtime,nhoridbg,nidbg) 
          ! no vertical axis
          cl_surf(1)='ter'
          cl_surf(2)='lic'
          cl_surf(3)='oce'
          cl_surf(4)='sic'
          DO nsrf=1,nbsrf
             CALL histdef(nidbg, cl_surf(nsrf),cl_surf(nsrf), "-",nbp_lon, &
                  nbp_lat,nhoridbg, 1, 1, 1, -99, 32, "inst", dtime,dtime) 
          ENDDO

          CALL histend(nidbg)
          CALL histsync(nidbg)

       ENDIF
       
    ENDIF
          
!****************************************************************************************
! Force soil water content to qsol0 if qsol0>0 and VEGET=F (use bucket
! instead of ORCHIDEE)
    IF (qsol0>=0.) THEN
      PRINT*,'WARNING : On impose qsol=',qsol0
      qsol(:)=qsol0
    ENDIF
!****************************************************************************************

!****************************************************************************************
! 2) Initialization to zero 
!****************************************************************************************
!
! 2a) Initialization of all argument variables with INTENT(OUT)
!****************************************************************************************
 cdragh(:)=0. ; cdragm(:)=0.
 zu1(:)=0. ; zv1(:)=0.
!albedo SB >>>
  alb_dir_m=0. ; alb_dif_m=0. ; alb3_lic(:)=0.
!albedo SB <<<
 zxsens(:)=0. ; zxevap(:)=0. ; zxtsol(:)=0.
 d_t_w(:,:)=0. ; d_q_w(:,:)=0. ; d_t_x(:,:)=0. ; d_q_x(:,:)=0.
 zxfluxlat(:)=0.
 zt2m(:)=0. ; zq2m(:)=0. ; qsat2m(:)=0. ; rh2m(:)=0.
 zn2mout(:,:)=0 ;
 d_t(:,:)=0. ; d_t_diss(:,:)=0. ; d_q(:,:)=0. ; d_u(:,:)=0. ; d_v(:,:)=0.
 zcoefh(:,:,:)=0. ; zcoefm(:,:,:)=0.
 zxsens_x(:)=0. ; zxsens_w(:)=0. ; zxfluxlat_x(:)=0. ; zxfluxlat_w(:)=0.
 cdragh_x(:)=0. ; cdragh_w(:)=0. ; cdragm_x(:)=0. ; cdragm_w(:)=0.
 kh(:)=0. ; kh_x(:)=0. ; kh_w(:)=0.
 slab_wfbils(:)=0.
 s_pblh(:)=0. ; s_pblh_x(:)=0. ; s_pblh_w(:)=0.
 s_plcl(:)=0. ; s_plcl_x(:)=0. ; s_plcl_w(:)=0.
 s_capCL(:)=0. ; s_oliqCL(:)=0. ; s_cteiCL(:)=0. ; s_pblT(:)=0.
 s_therm(:)=0.
 s_trmb1(:)=0. ; s_trmb2(:)=0. ; s_trmb3(:)=0.
 zustar(:)=0.
 zu10m(:)=0. ; zv10m(:)=0.
 fder_print(:)=0.
 zxqsurf(:)=0.
 zxfluxu(:,:)=0. ; zxfluxv(:,:)=0.
 solsw(:,:)=0. ; sollw(:,:)=0.
 d_ts(:,:)=0.
 evap(:,:)=0.
 fluxlat(:,:)=0.
 wfbils(:,:)=0. ; wfbilo(:,:)=0.
 wfevap(:,:)=0. ; wfrain(:,:)=0. ; wfsnow(:,:)=0.
 flux_t(:,:,:)=0. ; flux_q(:,:,:)=0. ; flux_u(:,:,:)=0. ; flux_v(:,:,:)=0.
 dflux_t(:)=0. ; dflux_q(:)=0.
 zxsnow(:)=0.
 zxfluxt(:,:)=0. ; zxfluxq(:,:)=0.
 qsnow(:)=0. ; snowhgt(:)=0. ; to_ice(:)=0. ; sissnow(:)=0.
 runoff(:)=0.
    IF (iflag_pbl<20.or.iflag_pbl>=30) THEN
       zcoefh(:,:,:) = 0.0
       zcoefh(:,1,:) = 999999. ! zcoefh(:,k=1) should never be used
       zcoefm(:,:,:) = 0.0
       zcoefm(:,1,:) = 999999. !
    ELSE
      zcoefm(:,:,is_ave)=0.
      zcoefh(:,:,is_ave)=0.
    ENDIF
!!
!  The components "is_ave" of tke_x and wake_deltke are "OUT" variables
!jyg<
!!    tke(:,:,is_ave)=0.
    tke_x(:,:,is_ave)=0.
    wake_dltke(:,:,is_ave)=0.
!>jyg
!!! jyg le 23/02/2013
    t2m(:,:)       = 999999.     ! t2m and q2m are meaningfull only over sub-surfaces
    q2m(:,:)       = 999999.     ! actually present in the grid cell.
!!!
    rh2m(:) = 0. ; qsat2m(:) = 0.
!!!
!!! jyg le 10/02/2012
    rh2m_x(:) = 0. ; qsat2m_x(:) = 0. ; rh2m_w(:) = 0. ; qsat2m_w(:) = 0.
!!!

! 2b) Initialization of all local variables that will be compressed later
!****************************************************************************************
!!    cdragh = 0.0  ; cdragm = 0.0     ; dflux_t = 0.0   ; dflux_q = 0.0
    ypct = 0.0    ; yts = 0.0        ; ysnow = 0.0
!!    zv1 = 0.0     ; yqsurf = 0.0
!albedo SB >>>
    yqsurf = 0.0  ; yalb = 0.0 ; yalb_vis = 0.0
!albedo SB <<<
    yrain_f = 0.0 ; ysnow_f = 0.0    ; yfder = 0.0     ; ysolsw = 0.0    
    ysollw = 0.0  ; yz0m = 0.0 ; yz0h = 0.0    ; yu1 = 0.0    
    yv1 = 0.0     ; ypaprs = 0.0     ; ypplay = 0.0
    ydelp = 0.0   ; yu = 0.0         ; yv = 0.0        ; yt = 0.0         
    yq = 0.0      ; y_dflux_t = 0.0  ; y_dflux_q = 0.0 
    yrugoro = 0.0 ; ywindsp = 0.0   
!!    d_ts = 0.0    ; yfluxlat=0.0     ; flux_t = 0.0    ; flux_q = 0.0     
    yfluxlat=0.0 
!!    flux_u = 0.0  ; flux_v = 0.0     ; d_t = 0.0       ; d_q = 0.0      
!!    d_t_diss= 0.0 ;d_u = 0.0     ; d_v = 0.0 
    yqsol = 0.0    
    ytherm = 0.0  ; ytke=0.
!FC 
    y_treedrg=0.

    ! Martin
    ysnowhgt = 0.0; yqsnow = 0.0     ; yrunoff = 0.0   ; ytoice =0.0
    yalb3_new = 0.0  ; ysissnow = 0.0 
    ypphi = 0.0   ; ycldt = 0.0      ; yrmu0 = 0.0
    ! Martin

!!! nrlmd+jyg le 02/05/2011 et le 20/02/2012
    ytke_x=0.     ; ytke_w=0.        ; ywake_dltke=0.
    y_d_t_x=0.    ; y_d_t_w=0.       ; y_d_q_x=0.      ; y_d_q_w=0.
!!    d_t_w=0.      ; d_q_w=0.         
!!    d_t_x=0.      ; d_q_x=0.
!!    d_wake_dlt=0.    ; d_wake_dlq=0.
    yfluxlat_x=0. ; yfluxlat_w=0.
    ywake_s=0.    ; ywake_cstar=0.   ;ywake_dens=0.
!!!
!!! nrlmd le 13/06/2011
    tau_eq=0.     ; delta_coef=0.
    y_delta_flux_t1=0.
    ydtsurf_th=0.
    yts_x=0.      ; yts_w=0.
    y_delta_tsurf=0.
!!!
    ytsoil = 999999. 
!FC
    y_d_u_frein(:,:)=0.
    y_d_v_frein(:,:)=0.
!FC

! >> PC
!the yfields_out variable is defined in (klon,nbcf_out) even if it is used on
!the ORCHIDEE grid and as such should be defined in yfields_out(knon,nbcf_out) but
!the knon variable is not known at that level of pbl_surface_mod

!the yfields_in variable is defined in (klon,nbcf_in) even if it is used on the
!ORCHIDEE grid and as such should be defined in yfields_in(knon,nbcf_in) but the
!knon variable is not known at that level of pbl_surface_mod

   yfields_out(:,:) = 0.
! << PC


! 2c) Initialization of all local variables computed within the subsurface loop and used later on
!****************************************************************************************
    d_t_diss_x(:,:) = 0. ;        d_t_diss_w(:,:) = 0.
    d_u_x(:,:)=0. ;               d_u_w(:,:)=0. 
    d_v_x(:,:)=0. ;               d_v_w(:,:)=0.
    flux_t_x(:,:,:)=0. ;          flux_t_w(:,:,:)=0. 
    flux_q_x(:,:,:)=0. ;          flux_q_w(:,:,:)=0.
!
!jyg<
    flux_u_x(:,:,:)=0. ;          flux_u_w(:,:,:)=0.
    flux_v_x(:,:,:)=0. ;          flux_v_w(:,:,:)=0.
    fluxlat_x(:,:)=0. ;           fluxlat_w(:,:)=0. 
!>jyg
!
!jyg<
! pblh,plcl,capCL,cteiCL ... are meaningfull only over sub-surfaces
! actually present in the grid cell  ==> value set to 999999.
!                           
!jyg<
       ustar(:,:)   = 999999.
       wstar(:,:)   = 999999.
       windsp(:,:)  = SQRT(u10m(:,:)**2 + v10m(:,:)**2 )
       u10m(:,:)    = 999999. 
       v10m(:,:)    = 999999. 
!>jyg
!
       pblh(:,:)   = 999999.        ! Hauteur de couche limite
       plcl(:,:)   = 999999.        ! Niveau de condensation de la CLA
       capCL(:,:)  = 999999.        ! CAPE de couche limite
       oliqCL(:,:) = 999999.        ! eau_liqu integree de couche limite
       cteiCL(:,:) = 999999.        ! cloud top instab. crit. couche limite
       pblt(:,:)   = 999999.        ! T a la Hauteur de couche limite
       therm(:,:)  = 999999.
       trmb1(:,:)  = 999999.        ! deep_cape
       trmb2(:,:)  = 999999.        ! inhibition 
       trmb3(:,:)  = 999999.        ! Point Omega
!
       t2m_x(:,:)    = 999999. 
       q2m_x(:,:)    = 999999. 
       ustar_x(:,:)   = 999999.
       wstar_x(:,:)   = 999999.
       u10m_x(:,:)   = 999999. 
       v10m_x(:,:)   = 999999. 
!                           
       pblh_x(:,:)   = 999999.      ! Hauteur de couche limite
       plcl_x(:,:)   = 999999.      ! Niveau de condensation de la CLA
       capCL_x(:,:)  = 999999.      ! CAPE de couche limite
       oliqCL_x(:,:) = 999999.      ! eau_liqu integree de couche limite
       cteiCL_x(:,:) = 999999.      ! cloud top instab. crit. couche limite
       pblt_x(:,:)   = 999999.      ! T a la Hauteur de couche limite
       therm_x(:,:)  = 999999.      
       trmb1_x(:,:)  = 999999.      ! deep_cape
       trmb2_x(:,:)  = 999999.      ! inhibition 
       trmb3_x(:,:)  = 999999.      ! Point Omega
!
       t2m_w(:,:)    = 999999. 
       q2m_w(:,:)    = 999999. 
       ustar_w(:,:)   = 999999.
       wstar_w(:,:)   = 999999.
       u10m_w(:,:)   = 999999. 
       v10m_w(:,:)   = 999999. 
                           
       pblh_w(:,:)   = 999999.      ! Hauteur de couche limite
       plcl_w(:,:)   = 999999.      ! Niveau de condensation de la CLA
       capCL_w(:,:)  = 999999.      ! CAPE de couche limite
       oliqCL_w(:,:) = 999999.      ! eau_liqu integree de couche limite
       cteiCL_w(:,:) = 999999.      ! cloud top instab. crit. couche limite
       pblt_w(:,:)   = 999999.      ! T a la Hauteur de couche limite
       therm_w(:,:)  = 999999.      
       trmb1_w(:,:)  = 999999.      ! deep_cape
       trmb2_w(:,:)  = 999999.      ! inhibition 
       trmb3_w(:,:)  = 999999.      ! Point Omega
!!!      
!
!!!
!****************************************************************************************
! 3) - Calculate pressure thickness of each layer
!    - Calculate the wind at first layer
!    - Mean calculations of albedo
!    - Calculate net radiance at sub-surface
!****************************************************************************************
    DO k = 1, klev
       DO i = 1, klon
          delp(i,k) = paprs(i,k)-paprs(i,k+1)
       ENDDO
    ENDDO

!****************************************************************************************
! Test for rugos........ from physiq.. A la fin plutot???
!
!****************************************************************************************

    DO nsrf = 1, nbsrf
       DO i = 1, klon
          z0m(i,nsrf) = MAX(z0m(i,nsrf),z0min)
          z0h(i,nsrf) = MAX(z0h(i,nsrf),z0min)
       ENDDO
    ENDDO

! Mean calculations of albedo
!
! * alb  : mean albedo for whole SW interval
!
! Mean albedo for grid point
! * alb_m  : mean albedo at whole SW interval

    alb_dir_m(:,:) = 0.0
    alb_dif_m(:,:) = 0.0
    DO k = 1, nsw
     DO nsrf = 1, nbsrf
       DO i = 1, klon
          alb_dir_m(i,k) = alb_dir_m(i,k) + alb_dir(i,k,nsrf) * pctsrf(i,nsrf)
          alb_dif_m(i,k) = alb_dif_m(i,k) + alb_dif(i,k,nsrf) * pctsrf(i,nsrf)
       ENDDO
     ENDDO
    ENDDO

! We here suppose the fraction f1 of incoming radiance of visible radiance 
! as a fraction of all shortwave radiance 
    f1 = 0.5
!    f1 = 1    ! put f1=1 to recreate old calculations

!f1 is already included with SFRWL values in each surf files
    alb=0.0
    DO k=1,nsw
      DO nsrf = 1, nbsrf
        DO i = 1, klon
            alb(i,nsrf) = alb(i,nsrf) + alb_dir(i,k,nsrf)*SFRWL(k)
        ENDDO
      ENDDO
    ENDDO

    alb_m=0.0
    DO k = 1,nsw
      DO i = 1, klon
        alb_m(i) = alb_m(i) + alb_dir_m(i,k)*SFRWL(k)
      ENDDO
    ENDDO
!albedo SB <<<



! Calculation of mean temperature at surface grid points
    ztsol(:) = 0.0
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          ztsol(i) = ztsol(i) + ts(i,nsrf)*pctsrf(i,nsrf)
       ENDDO
    ENDDO

! Linear distrubution on sub-surface of long- and shortwave net radiance
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          sollw(i,nsrf) = sollw_m(i) + 4.0*RSIGMA*ztsol(i)**3 * (ztsol(i)-ts(i,nsrf))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! Martin
! Apparently introduced for sisvat but not used
!         sollwd(i,nsrf)= sollwd_m(i)
!         ! Martin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          solsw(i,nsrf) = solsw_m(i) * (1.-alb(i,nsrf)) / (1.-alb_m(i))
       ENDDO
    ENDDO
!
!<al1: second order corrections
!- net = dwn -up; up=sig( T4 + 4sum%T3T' + 6sum%T2T'2 +...)
   IF (iflag_order2_sollw == 1) THEN
    meansqT(:) = 0. ! as working buffer
    DO nsrf = 1, nbsrf
     DO i = 1, klon
      meansqT(i) = meansqT(i)+(ts(i,nsrf)-ztsol(i))**2 *pctsrf(i,nsrf)
     ENDDO
    ENDDO
    DO nsrf = 1, nbsrf
     DO i = 1, klon
      sollw(i,nsrf) = sollw(i,nsrf) &
                + 6.0*RSIGMA*ztsol(i)**2 *(meansqT(i)-(ztsol(i)-ts(i,nsrf))**2)
     ENDDO
    ENDDO
   ENDIF   ! iflag_order2_sollw == 1
!>al1

! >> PC
   IF (carbon_cycle_cpl .AND. carbon_cycle_tr .AND. nbcf_out.GT.0 ) THEN
       r_co2_ppm(:) = co2_send(:)
       DO n=1,nbcf_out
           IF (cfname_out(n) == "atmco2" ) fields_out(:,n) = co2_send(:)
       ENDDO
   ENDIF
   IF ( .NOT. carbon_cycle_tr .AND. nbcf_out.GT.0 ) THEN
       r_co2_ppm(:) = co2_ppm     ! Constant field
       DO n=1,nbcf_out
           IF (cfname_out(n) == "atmco2" ) fields_out(:,n) = co2_ppm
       ENDDO
   ENDIF
! << PC

!****************************************************************************************
! 4) Loop over different surfaces
!
! Only points containing a fraction of the sub surface will be treated.
! 
!****************************************************************************************
   
    loop_nbsrf: DO nsrf = 1, nbsrf
       IF (prt_level >=10) print *,' Loop nsrf ',nsrf

! Search for index(ni) and size(knon) of domaine to treat
       ni(:) = 0
       knon  = 0
       DO i = 1, klon
          IF (pctsrf(i,nsrf) > 0.) THEN
             knon = knon + 1
             ni(knon) = i
          ENDIF
       ENDDO

!!! jyg le 19/08/2012
!       IF (knon <= 0) THEN
!         IF (prt_level >= 10) print *,' no grid point for nsrf= ',nsrf
!         cycle loop_nbsrf
!       ENDIF
!!!

       ! write index, with IOIPSL
       IF (debugindex .AND. mpi_size==1) THEN 
          tabindx(:)=0.
          DO i=1,knon
             tabindx(i)=REAL(i)
          ENDDO
          debugtab(:,:) = 0.
          ndexbg(:) = 0
          CALL gath2cpl(tabindx,debugtab,knon,ni)
          CALL histwrite(nidbg,cl_surf(nsrf),itap,debugtab,nbp_lon*nbp_lat, ndexbg)
       ENDIF
       
!****************************************************************************************
! 5) Compress variables 
!
!****************************************************************************************

       DO j = 1, knon
          i = ni(j)
          ypct(j)    = pctsrf(i,nsrf)
          yts(j)     = ts(i,nsrf)
          ysnow(j)   = snow(i,nsrf)
          yqsurf(j)  = qsurf(i,nsrf)
          yalb(j)    = alb(i,nsrf)
!albedo SB >>>
          yalb_vis(j) = alb_dir(i,1,nsrf)
          IF (nsw==6) THEN
            yalb_vis(j)=(alb_dir(i,1,nsrf)*SFRWL(1)+alb_dir(i,2,nsrf)*SFRWL(2) &
              +alb_dir(i,3,nsrf)*SFRWL(3))/(SFRWL(1)+SFRWL(2)+SFRWL(3))
          ENDIF
!albedo SB <<<
          yrain_f(j) = rain_f(i)
          ysnow_f(j) = snow_f(i)
          yagesno(j) = agesno(i,nsrf)
          yfder(j)   = fder(i)
          ylwdown(j) = lwdown_m(i)
          ygustiness(j) = gustiness(i)
          ysolsw(j)  = solsw(i,nsrf)
          ysollw(j)  = sollw(i,nsrf)
          yz0m(j)  = z0m(i,nsrf)
          yz0h(j)  = z0h(i,nsrf)
          yrugoro(j) = rugoro(i)
          yu1(j)     = u(i,1)
          yv1(j)     = v(i,1)
          ypaprs(j,klev+1) = paprs(i,klev+1)
!jyg<
!!          ywindsp(j) = SQRT(u10m(i,nsrf)**2 + v10m(i,nsrf)**2 )
          ywindsp(j) = windsp(i,nsrf)
!>jyg
          ! Martin
          yzsig(j)   = zsig(i)
          ycldt(j)   = cldt(i)
          yrmu0(j)   = rmu0(i)
          ! Martin
!!! nrlmd le 13/06/2011
          y_delta_tsurf(j)=delta_tsurf(i,nsrf)
!!!
       ENDDO
! >> PC
!--compressing fields_out onto ORCHIDEE grid
!--these fields are shared and used directly surf_land_orchidee_mod
       DO n = 1, nbcf_out
         DO j = 1, knon
           i = ni(j)
           yfields_out(j,n) = fields_out(i,n)
         ENDDO
       ENDDO
! << PC
       DO k = 1, klev
          DO j = 1, knon
             i = ni(j)
             ypaprs(j,k) = paprs(i,k)
             ypplay(j,k) = pplay(i,k)
             ydelp(j,k)  = delp(i,k)
          ENDDO
       ENDDO
!
!!! jyg le 07/02/2012 et le 10/04/2013
        DO k = 1, klev+1
          DO j = 1, knon
             i = ni(j)
!jyg<
!!             ytke(j,k)   = tke(i,k,nsrf)
             ytke(j,k)   = tke_x(i,k,nsrf)
          ENDDO
        ENDDO
!>jyg
        DO k = 1, klev
          DO j = 1, knon
             i = ni(j)
             y_treedrg(j,k) =  treedrg(i,k,nsrf)
             yu(j,k) = u(i,k)
             yv(j,k) = v(i,k)
             yt(j,k) = t(i,k)
             yq(j,k) = q(i,k)
          ENDDO
        ENDDO
!
       IF (iflag_split.GE.1) THEN
!!! nrlmd le 02/05/2011
        DO k = 1, klev
          DO j = 1, knon
             i = ni(j)
             yu_x(j,k) = u(i,k)
             yv_x(j,k) = v(i,k)
             yt_x(j,k) = t(i,k)-wake_s(i)*wake_dlt(i,k)
             yq_x(j,k) = q(i,k)-wake_s(i)*wake_dlq(i,k)
             yu_w(j,k) = u(i,k)
             yv_w(j,k) = v(i,k)
             yt_w(j,k) = t(i,k)+(1.-wake_s(i))*wake_dlt(i,k)
             yq_w(j,k) = q(i,k)+(1.-wake_s(i))*wake_dlq(i,k)
!!!
          ENDDO
        ENDDO

        IF (prt_level .ge. 10) THEN
          print *,'pbl_surface, wake_s(1), wake_dlt(1,:) ', wake_s(1), wake_dlt(1,:)
          print *,'pbl_surface, wake_s(1), wake_dlq(1,:) ', wake_s(1), wake_dlq(1,:)
        ENDIF

!!! nrlmd le 02/05/2011
        DO k = 1, klev+1
          DO j = 1, knon
             i = ni(j)
!jyg<
!!             ytke_x(j,k) = tke(i,k,nsrf)-wake_s(i)*wake_dltke(i,k,nsrf)
!!             ytke_w(j,k) = tke(i,k,nsrf)+(1.-wake_s(i))*wake_dltke(i,k,nsrf)
!!             ywake_dltke(j,k) = wake_dltke(i,k,nsrf)
!!             ytke(j,k)     = tke(i,k,nsrf)
!
             ytke_x(j,k)      = tke_x(i,k,nsrf)
             ytke(j,k)        = tke_x(i,k,nsrf)+wake_s(i)*wake_dltke(i,k,nsrf)
             ytke_w(j,k)      = tke_x(i,k,nsrf)+wake_dltke(i,k,nsrf)
             ywake_dltke(j,k) = wake_dltke(i,k,nsrf)
!>jyg
          ENDDO
        ENDDO
!!!
!!! jyg le 07/02/2012
        DO j = 1, knon
          i = ni(j)
          ywake_s(j)=wake_s(i)
          ywake_cstar(j)=wake_cstar(i)
          ywake_dens(j)=wake_dens(i)
        ENDDO
!!!
!!! nrlmd le 13/06/2011
        DO j=1,knon
         yts_x(j)=yts(j)-ywake_s(j)*y_delta_tsurf(j)
         yts_w(j)=yts(j)+(1.-ywake_s(j))*y_delta_tsurf(j)
        ENDDO
!!!
       ENDIF  ! (iflag_split .ge.1)
!!!
       DO k = 1, nsoilmx
          DO j = 1, knon
             i = ni(j)
             ytsoil(j,k) = ftsoil(i,k,nsrf)
          ENDDO
       ENDDO
       
       ! qsol(water height in soil) only for bucket continental model
       IF ( nsrf .EQ. is_ter .AND. .NOT. ok_veget ) THEN 
          DO j = 1, knon
             i = ni(j)
             yqsol(j) = qsol(i)
          ENDDO
       ENDIF
       
!****************************************************************************************
! 6a) Calculate coefficients for turbulent diffusion at surface, cdragh et cdragm.
!
!****************************************************************************************

!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
!!!
!!! nrlmd & jyg les 02/05/2011, 13/06/2011, 05/02/2012
! Faire disparaitre les lignes commentees fin 2015 (le temps des tests)
!       CALL clcdrag( knon, nsrf, ypaprs, ypplay, &
!           yu(:,1), yv(:,1), yt(:,1), yq(:,1), &
!           yts, yqsurf, yrugos, &
!           ycdragm, ycdragh )
! Fuxing WANG, 04/03/2015, replace the clcdrag by the merged version: cdrag
        DO i = 1, knon
!          print*,'PBL ',i,RD
!          print*,'PBL ',yt(i,1),ypaprs(i,1),ypplay(i,1)
           zgeo1(i) = RD * yt(i,1) / (0.5*(ypaprs(i,1)+ypplay(i,1))) &
                * (ypaprs(i,1)-ypplay(i,1))
           speed(i) = SQRT(yu(i,1)**2+yv(i,1)**2)
        ENDDO
        CALL cdrag(knon, nsrf, &
            speed, yt(:,1), yq(:,1), zgeo1, ypaprs(:,1),&
            yts, yqsurf, yz0m, yz0h, &
            ycdragm, ycdragh, zri1, pref )

! --- special Dice: on force cdragm ( a defaut de forcer ustar) MPL 05082013
     IF (ok_prescr_ust) THEN
      DO i = 1, knon
       print *,'ycdragm avant=',ycdragm(i)
       vent= sqrt(yu(i,1)*yu(i,1)+yv(i,1)*yv(i,1))
!      ycdragm(i) = ust*ust/(1.+(yu(i,1)*yu(i,1)+yv(i,1)*yv(i,1)))
!      ycdragm(i) = ust*ust/((1.+sqrt(yu(i,1)*yu(i,1)+yv(i,1)*yv(i,1))) &
!     *sqrt(yu(i,1)*yu(i,1)+yv(i,1)*yv(i,1)))
       ycdragm(i) = ust*ust/(1.+vent)/vent
!      print *,'ycdragm ust yu yv apres=',ycdragm(i),ust,yu(i,1),yv(i,1)
      ENDDO
     ENDIF
        IF (prt_level >=10) print *,'clcdrag -> ycdragh ', ycdragh
       ELSE  !(iflag_split .eq.0)

! Faire disparaitre les lignes commentees fin 2015 (le temps des tests)
!       CALL clcdrag( knon, nsrf, ypaprs, ypplay, &
!           yu_x(:,1), yv_x(:,1), yt_x(:,1), yq_x(:,1), &
!           yts_x, yqsurf, yrugos, &
!           ycdragm_x, ycdragh_x )
! Fuxing WANG, 04/03/2015, replace the clcdrag by the merged version: cdrag
        DO i = 1, knon
           zgeo1_x(i) = RD * yt_x(i,1) / (0.5*(ypaprs(i,1)+ypplay(i,1))) &
                * (ypaprs(i,1)-ypplay(i,1))
           speed_x(i) = SQRT(yu_x(i,1)**2+yv_x(i,1)**2)
        ENDDO
        CALL cdrag(knon, nsrf, &
            speed_x, yt_x(:,1), yq_x(:,1), zgeo1_x, ypaprs(:,1),&
            yts_x, yqsurf, yz0m, yz0h, &
            ycdragm_x, ycdragh_x, zri1_x, pref_x )

! --- special Dice. JYG+MPL 25112013
        IF (ok_prescr_ust) THEN
         DO i = 1, knon
!         print *,'ycdragm_x avant=',ycdragm_x(i)
          vent= sqrt(yu_x(i,1)*yu_x(i,1)+yv_x(i,1)*yv_x(i,1))
          ycdragm_x(i) = ust*ust/(1.+vent)/vent
!         print *,'ycdragm_x ust yu yv apres=',ycdragm_x(i),ust,yu_x(i,1),yv_x(i,1)
         ENDDO
        ENDIF
        IF (prt_level >=10) print *,'clcdrag -> ycdragh_x ', ycdragh_x
!
! Faire disparaitre les lignes commentees fin 2015 (le temps des tests)
!        CALL clcdrag( knon, nsrf, ypaprs, ypplay, &
!            yu_w(:,1), yv_w(:,1), yt_w(:,1), yq_w(:,1), &
!            yts_w, yqsurf, yz0m, &
!            ycdragm_w, ycdragh_w )
! Fuxing WANG, 04/03/2015, replace the clcdrag by the merged version: cdrag
        DO i = 1, knon
           zgeo1_w(i) = RD * yt_w(i,1) / (0.5*(ypaprs(i,1)+ypplay(i,1))) &
                * (ypaprs(i,1)-ypplay(i,1))
           speed_w(i) = SQRT(yu_w(i,1)**2+yv_w(i,1)**2)
        ENDDO
        CALL cdrag(knon, nsrf, &
            speed_w, yt_w(:,1), yq_w(:,1), zgeo1_w, ypaprs(:,1),&
            yts_w, yqsurf, yz0m, yz0h, &
            ycdragm_w, ycdragh_w, zri1_w, pref_w )
!
        zgeo1(:) = wake_s(:)*zgeo1_w(:) + (1.-wake_s(:))*zgeo1_x(:)

! --- special Dice. JYG+MPL 25112013 puis BOMEX
        IF (ok_prescr_ust) THEN
         DO i = 1, knon
!         print *,'ycdragm_w avant=',ycdragm_w(i)
          vent= sqrt(yu_w(i,1)*yu_w(i,1)+yv_w(i,1)*yv_w(i,1))
          ycdragm_w(i) = ust*ust/(1.+vent)/vent
!         print *,'ycdragm_w ust yu yv apres=',ycdragm_w(i),ust,yu_w(i,1),yv_w(i,1)
         ENDDO
        ENDIF
        IF (prt_level >=10) print *,'clcdrag -> ycdragh_w ', ycdragh_w
!!!
       ENDIF  ! (iflag_split .eq.0)
!!!
       

!****************************************************************************************
! 6b) Calculate coefficients for turbulent diffusion in the atmosphere, ycoefh et ycoefm.
!
!****************************************************************************************

!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
!!!
!!! nrlmd & jyg les 02/05/2011, 13/06/2011, 05/02/2012
      IF (prt_level >=10) THEN
      print *,' args coef_diff_turb: yu ',  yu  
      print *,' args coef_diff_turb: yv ',  yv  
      print *,' args coef_diff_turb: yq ',  yq  
      print *,' args coef_diff_turb: yt ',  yt  
      print *,' args coef_diff_turb: yts ', yts  
      print *,' args coef_diff_turb: yz0m ', yz0m  
      print *,' args coef_diff_turb: yqsurf ', yqsurf  
      print *,' args coef_diff_turb: ycdragm ', ycdragm
      print *,' args coef_diff_turb: ycdragh ', ycdragh
      print *,' args coef_diff_turb: ytke ', ytke
       ENDIF
        CALL coef_diff_turb(dtime, nsrf, knon, ni,  &
            ypaprs, ypplay, yu, yv, yq, yt, yts, yqsurf, ycdragm, &
            ycoefm, ycoefh, ytke, y_treedrg)
!            ycoefm, ycoefh, ytke)
!FC y_treedrg ajoute
       IF (iflag_pbl>=20.AND.iflag_pbl<30) THEN
! In this case, coef_diff_turb is called for the Cd only
       DO k = 2, klev
          DO j = 1, knon
             i = ni(j)
             ycoefh(j,k)   = zcoefh(i,k,nsrf)
             ycoefm(j,k)   = zcoefm(i,k,nsrf)
          ENDDO
       ENDDO
       ENDIF
        IF (prt_level >=10) print *,'coef_diff_turb -> ycoefh ',ycoefh
!
       ELSE  !(iflag_split .eq.0)
      IF (prt_level >=10) THEN
      print *,' args coef_diff_turb: yu_x ',  yu_x  
      print *,' args coef_diff_turb: yv_x ',  yv_x  
      print *,' args coef_diff_turb: yq_x ',  yq_x  
      print *,' args coef_diff_turb: yt_x ',  yt_x  
      print *,' args coef_diff_turb: yts_x ', yts_x  
      print *,' args coef_diff_turb: yqsurf ', yqsurf  
      print *,' args coef_diff_turb: ycdragm_x ', ycdragm_x
      print *,' args coef_diff_turb: ycdragh_x ', ycdragh_x
      print *,' args coef_diff_turb: ytke_x ', ytke_x
       ENDIF
        CALL coef_diff_turb(dtime, nsrf, knon, ni,  &
            ypaprs, ypplay, yu_x, yv_x, yq_x, yt_x, yts_x, yqsurf, ycdragm_x, &
            ycoefm_x, ycoefh_x, ytke_x,y_treedrg)
!            ycoefm_x, ycoefh_x, ytke_x)
!FC doit on le mettre ( on ne l utilise pas si il y a du spliting)
       IF (iflag_pbl>=20.AND.iflag_pbl<30) THEN
! In this case, coef_diff_turb is called for the Cd only
       DO k = 2, klev
          DO j = 1, knon
             i = ni(j)
             ycoefh_x(j,k)   = zcoefh(i,k,nsrf)
             ycoefm_x(j,k)   = zcoefm(i,k,nsrf)
          ENDDO
       ENDDO
       ENDIF
        IF (prt_level >=10) print *,'coef_diff_turb -> ycoefh_x ',ycoefh_x
!
      IF (prt_level >=10) THEN
      print *,' args coef_diff_turb: yu_w ',  yu_w  
      print *,' args coef_diff_turb: yv_w ',  yv_w  
      print *,' args coef_diff_turb: yq_w ',  yq_w  
      print *,' args coef_diff_turb: yt_w ',  yt_w  
      print *,' args coef_diff_turb: yts_w ', yts_w  
      print *,' args coef_diff_turb: yqsurf ', yqsurf  
      print *,' args coef_diff_turb: ycdragm_w ', ycdragm_w
      print *,' args coef_diff_turb: ycdragh_w ', ycdragh_w
      print *,' args coef_diff_turb: ytke_w ', ytke_w
       ENDIF
        CALL coef_diff_turb(dtime, nsrf, knon, ni,  &
            ypaprs, ypplay, yu_w, yv_w, yq_w, yt_w, yts_w, yqsurf, ycdragm_w, &
            ycoefm_w, ycoefh_w, ytke_w,y_treedrg)
!            ycoefm_w, ycoefh_w, ytke_w)
       IF (iflag_pbl>=20.AND.iflag_pbl<30) THEN
! In this case, coef_diff_turb is called for the Cd only
       DO k = 2, klev
          DO j = 1, knon
             i = ni(j)
             ycoefh_w(j,k)   = zcoefh(i,k,nsrf)
             ycoefm_w(j,k)   = zcoefm(i,k,nsrf)
          ENDDO
       ENDDO
       ENDIF
        IF (prt_level >=10) print *,'coef_diff_turb -> ycoefh_w ',ycoefh_w
!
!!!jyg le 10/04/2013
!!   En attendant de traiter le transport des traceurs dans les poches froides, formule
!!   arbitraire pour ycoefh et ycoefm
      DO k = 2,klev
        DO j = 1,knon
         ycoefh(j,k) = ycoefh_x(j,k) + ywake_s(j)*(ycoefh_w(j,k) - ycoefh_x(j,k))
         ycoefm(j,k) = ycoefm_x(j,k) + ywake_s(j)*(ycoefm_w(j,k) - ycoefm_x(j,k))
        ENDDO
      ENDDO
!!!
       ENDIF  ! (iflag_split .eq.0)
!!!
       
!****************************************************************************************
! 
! 8) "La descente" - "The downhill"
!  
!  climb_hq_down and climb_wind_down calculate the coefficients
!  Ccoef_X et Dcoef_X for X=[H, Q, U, V].
!  Only the coefficients at surface for H and Q are returned.
!
!****************************************************************************************

! - Calculate the coefficients Ccoef_H, Ccoef_Q, Dcoef_H and Dcoef_Q 
!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
!!!
!!! nrlmd & jyg les 02/05/2011, 13/06/2011, 05/02/2012
        CALL climb_hq_down(knon, ycoefh, ypaprs, ypplay, &
            ydelp, yt, yq, dtime, &
!!! jyg le 09/05/2011
            CcoefH, CcoefQ, DcoefH, DcoefQ, &
            Kcoef_hq, gama_q, gama_h, &
!!!
            AcoefH, AcoefQ, BcoefH, BcoefQ)
       ELSE  !(iflag_split .eq.0)
        CALL climb_hq_down(knon, ycoefh_x, ypaprs, ypplay, &
            ydelp, yt_x, yq_x, dtime, &
!!! nrlmd le 02/05/2011
            CcoefH_x, CcoefQ_x, DcoefH_x, DcoefQ_x, &
            Kcoef_hq_x, gama_q_x, gama_h_x, &
!!!
            AcoefH_x, AcoefQ_x, BcoefH_x, BcoefQ_x)
!!!
       IF (prt_level >=10) THEN
         PRINT *,'pbl_surface (climb_hq_down.x->) AcoefH_x ',AcoefH_x
         PRINT *,'pbl_surface (climb_hq_down.x->) AcoefQ_x ',AcoefQ_x
         PRINT *,'pbl_surface (climb_hq_down.x->) BcoefH_x ',BcoefH_x
         PRINT *,'pbl_surface (climb_hq_down.x->) BcoefQ_x ',BcoefQ_x
       ENDIF
!
        CALL climb_hq_down(knon, ycoefh_w, ypaprs, ypplay, &
            ydelp, yt_w, yq_w, dtime, &
!!! nrlmd le 02/05/2011
            CcoefH_w, CcoefQ_w, DcoefH_w, DcoefQ_w, &
            Kcoef_hq_w, gama_q_w, gama_h_w, &
!!!
            AcoefH_w, AcoefQ_w, BcoefH_w, BcoefQ_w)
!!!
       IF (prt_level >=10) THEN
         PRINT *,'pbl_surface (climb_hq_down.w->) AcoefH_w ',AcoefH_w
         PRINT *,'pbl_surface (climb_hq_down.w->) AcoefQ_w ',AcoefQ_w
         PRINT *,'pbl_surface (climb_hq_down.w->) BcoefH_w ',BcoefH_w
         PRINT *,'pbl_surface (climb_hq_down.w->) BcoefQ_w ',BcoefQ_w
       ENDIF
!!!
       ENDIF  ! (iflag_split .eq.0)
!!!

! - Calculate the coefficients Ccoef_U, Ccoef_V, Dcoef_U and Dcoef_V
!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
!!! nrlmd & jyg les 02/05/2011, 13/06/2011, 05/02/2012
        CALL climb_wind_down(knon, dtime, ycoefm, ypplay, ypaprs, yt, ydelp, yu, yv, &
!!! jyg le 09/05/2011
            CcoefU, CcoefV, DcoefU, DcoefV, &
            Kcoef_m, alf_1, alf_2, &
!!!
            AcoefU, AcoefV, BcoefU, BcoefV)
       ELSE  ! (iflag_split .eq.0)
        CALL climb_wind_down(knon, dtime, ycoefm_x, ypplay, ypaprs, yt_x, ydelp, yu_x, yv_x, &
!!! nrlmd le 02/05/2011
            CcoefU_x, CcoefV_x, DcoefU_x, DcoefV_x, &
            Kcoef_m_x, alf_1_x, alf_2_x, &
!!!
            AcoefU_x, AcoefV_x, BcoefU_x, BcoefV_x)
!
        CALL climb_wind_down(knon, dtime, ycoefm_w, ypplay, ypaprs, yt_w, ydelp, yu_w, yv_w, &
!!! nrlmd le 02/05/2011
            CcoefU_w, CcoefV_w, DcoefU_w, DcoefV_w, &
            Kcoef_m_w, alf_1_w, alf_2_w, &
!!!
            AcoefU_w, AcoefV_w, BcoefU_w, BcoefV_w)
!!!      
       ENDIF  ! (iflag_split .eq.0)
!!!

!****************************************************************************************
! 9) Small calculations
!
!****************************************************************************************

! - Reference pressure is given the values at surface level          
       ypsref(:) = ypaprs(:,1)  

! - CO2 field on 2D grid to be sent to ORCHIDEE
!   Transform to compressed field
       IF (carbon_cycle_cpl) THEN
          DO i=1,knon
             r_co2_ppm(i) = co2_send(ni(i))
          ENDDO
       ELSE
          r_co2_ppm(:) = co2_ppm     ! Constant field
       ENDIF

!!! nrlmd le 02/05/2011  -----------------------On raccorde les 2 colonnes dans la couche 1 
!----------------------------------------------------------------------------------------
!!! jyg le 07/02/2012
!!! jyg le 01/02/2017
       IF (iflag_split .eq. 0) THEN
         yt1(:) = yt(:,1)
         yq1(:) = yq(:,1)
!!       ELSE IF (iflag_split .eq. 1) THEN
!!!
!jyg<
!!         CALL wx_pbl_fuse_no_dts(knon, dtime, ypplay, ywake_s, &
!!                                 yt_x, yt_w, yq_x, yq_w, &
!!                                 yu_x, yu_w, yv_x, yv_w, &
!!                                 ycdragh_x, ycdragh_w, ycdragm_x, ycdragm_w, &
!!                                 AcoefH_x, AcoefH_w, AcoefQ_x, AcoefQ_w, &
!!                                 AcoefU_x, AcoefU_w, AcoefV_x, AcoefV_w, &
!!                                 BcoefH_x, BcoefH_w, BcoefQ_x, BcoefQ_w, &
!!                                 BcoefU_x, BcoefU_w, BcoefV_x, BcoefV_w, &
!!                                 AcoefH, AcoefQ, AcoefU, AcoefV, &
!!                                 BcoefH, BcoefQ, BcoefU, BcoefV, &
!!                                 ycdragh, ycdragm, &
!!                                 yt1, yq1, yu1, yv1 &
!!                                 )
       ELSE IF (iflag_split .ge. 1) THEN
         CALL wx_pbl0_fuse(knon, dtime, ypplay, ywake_s, &
                         yt_x, yt_w, yq_x, yq_w, &
                         yu_x, yu_w, yv_x, yv_w, &
                         ycdragh_x, ycdragh_w, ycdragm_x, ycdragm_w, &
                         AcoefH_x, AcoefH_w, AcoefQ_x, AcoefQ_w, &
                         AcoefU_x, AcoefU_w, AcoefV_x, AcoefV_w, &
                         BcoefH_x, BcoefH_w, BcoefQ_x, BcoefQ_w, &
                         BcoefU_x, BcoefU_w, BcoefV_x, BcoefV_w, &
                         AcoefH, AcoefQ, AcoefU, AcoefV, &
                         BcoefH, BcoefQ, BcoefU, BcoefV, &
                         ycdragh, ycdragm, &
                         yt1, yq1, yu1, yv1 &
                         )
!!       ELSE IF (iflag_split .ge.2) THEN
!!!    Provisoire
!!         ah(:) = 0.
!!         bh(:) = 0.
!!         IF (nsrf == is_oce) THEN
!!           ybeta(:) = 1.
!!         ELSE
!!           ybeta(:) = beta_land
!!         ENDIF
!!         ycdragh(:) = ywake_s(:)*ycdragh_w(:) + (1.-ywake_s(:))*ycdragh_x(:)
!!         CALL wx_dts(knon, nsrf, ywake_cstar, ywake_s, ywake_dens, &
!!                     yts, ypplay(:,1), ybeta, ycdragh , ypaprs(:,1), &
!!                     yq(:,1), yt(:,1), yu(:,1), yv(:,1), ygustiness, &
!!                     ah, bh &
!!                     )
!!!
!!         CALL wx_pbl_fuse(knon, dtime, ypplay, ywake_s, &
!!                         yt_x, yt_w, yq_x, yq_w, &
!!                         yu_x, yu_w, yv_x, yv_w, &
!!                         ycdragh_x, ycdragh_w, ycdragm_x, ycdragm_w, &
!!                         AcoefH_x, AcoefH_w, AcoefQ_x, AcoefQ_w, &
!!                         AcoefU_x, AcoefU_w, AcoefV_x, AcoefV_w, &
!!                         BcoefH_x, BcoefH_w, BcoefQ_x, BcoefQ_w, &
!!                         BcoefU_x, BcoefU_w, BcoefV_x, BcoefV_w, &
!!                         ah, bh, &
!!                         AcoefH, AcoefQ, AcoefU, AcoefV, &
!!                         BcoefH, BcoefQ, BcoefU, BcoefV, &
!!                         ycdragh, ycdragm, &
!!                         yt1, yq1, yu1, yv1 &
!!                         )
!>jyg
!!!
         ENDIF  ! (iflag_split .eq.0)
!!!
       IF (prt_level >=10) THEN
         PRINT *,'pbl_surface (fuse->): yt(1,:) ',yt(1,:)
         PRINT *,'pbl_surface (fuse->): yq(1,:) ',yq(1,:)
         PRINT *,'pbl_surface (fuse->): yu(1,:) ',yu(1,:)
         PRINT *,'pbl_surface (fuse->): yv(1,:) ',yv(1,:)
         PRINT *,'pbl_surface (fuse->): AcoefH(1) ',AcoefH(1)
         PRINT *,'pbl_surface (fuse->): BcoefH(1) ',BcoefH(1)
       ENDIF

!****************************************************************************************
!
! Calulate t2m and q2m for the case of calculation at land grid points 
! t2m and q2m are needed as input to ORCHIDEE
!
!****************************************************************************************
       IF (nsrf == is_ter) THEN

          DO i = 1, knon
             zgeo1(i) = RD * yt(i,1) / (0.5*(ypaprs(i,1)+ypplay(i,1))) &
                  * (ypaprs(i,1)-ypplay(i,1))
          ENDDO

          ! Calculate the temperature et relative humidity at 2m and the wind at 10m 
          IF (iflag_new_t2mq2m==1) THEN
           CALL stdlevvarn(klon, knon, is_ter, zxli, &
               yu(:,1), yv(:,1), yt(:,1), yq(:,1), zgeo1, &
               yts, yqsurf, yz0m, yz0h, ypaprs(:,1), ypplay(:,1), &
               yt2m, yq2m, yt10m, yq10m, yu10m, yustar, &
               yn2mout(:, nsrf, :))
          ELSE
          CALL stdlevvar(klon, knon, is_ter, zxli, &
               yu(:,1), yv(:,1), yt(:,1), yq(:,1), zgeo1, &
               yts, yqsurf, yz0m, yz0h, ypaprs(:,1), ypplay(:,1), &
               yt2m, yq2m, yt10m, yq10m, yu10m, yustar)
          ENDIF
       ENDIF

!****************************************************************************************
!
! 10) Switch according to current surface
!     It is necessary to start with the continental surfaces because the ocean
!     needs their run-off.
!
!****************************************************************************************
       SELECT CASE(nsrf)
     
       CASE(is_ter)
!          print*,"DEBUGTS",yts(knon/2),ylwdown(knon/2)
          CALL surf_land(itap, dtime, date0, jour, knon, ni,&
               rlon, rlat, yrmu0, &
               debut, lafin, ydelp(:,1), r_co2_ppm, ysolsw, ysollw, yalb, &
!!jyg               yts, ypplay(:,1), ycdragh, ycdragm, yrain_f, ysnow_f, yt(:,1), yq(:,1),&
               yts, ypplay(:,1), ycdragh, ycdragm, yrain_f, ysnow_f, yt1, yq1,&
               AcoefH, AcoefQ, BcoefH, BcoefQ, & 
               AcoefU, AcoefV, BcoefU, BcoefV, & 
               ypsref, yu1, yv1, ygustiness, yrugoro, pctsrf, &
               ylwdown, yq2m, yt2m, &
               ysnow, yqsol, yagesno, ytsoil, &
               yz0m, yz0h, SFRWL, yalb_dir_new, yalb_dif_new, yevap, yfluxsens,yfluxlat,&
               yqsurf, ytsurf_new, y_dflux_t, y_dflux_q, &
               y_flux_u1, y_flux_v1, &
               yveget,ylai,yheight )
 
!FC quid qd yveget ylai yheight ne sont pas definit
!FC  yveget,ylai,yheight, &
            IF (ifl_pbltree .ge. 1) THEN
              CALL   freinage(knon, yu, yv, yt, &
!                yveget,ylai, yheight,ypaprs,ypplay,y_d_u_frein,y_d_v_frein)
                yveget,ylai, yheight,ypaprs,ypplay,y_treedrg, y_d_u_frein,y_d_v_frein)
            ENDIF

               
! Special DICE MPL 05082013 puis BOMEX
       IF (ok_prescr_ust) THEN
          DO j=1,knon
!         ysnow(:)=0.
!         yqsol(:)=0.
!         yagesno(:)=50.
!         ytsoil(:,:)=300.
!         yz0_new(:)=0.001
!         yevap(:)=flat/RLVTT
!         yfluxlat(:)=-flat
!         yfluxsens(:)=-fsens
!         yqsurf(:)=0.
!         ytsurf_new(:)=tg
!         y_dflux_t(:)=0.
!         y_dflux_q(:)=0.
          y_flux_u1(j)=ycdragm(j)*(1.+sqrt(yu(j,1)*yu(j,1)+yv(j,1)*yv(j,1)))*yu(j,1)*ypplay(j,1)/RD/yt(j,1)
          y_flux_v1(j)=ycdragm(j)*(1.+sqrt(yu(j,1)*yu(j,1)+yv(j,1)*yv(j,1)))*yv(j,1)*ypplay(j,1)/RD/yt(j,1)
          ENDDO
      ENDIF
     
       CASE(is_lic)
          ! Martin
          CALL surf_landice(itap, dtime, knon, ni, &
               rlon, rlat, debut, lafin, &
               yrmu0, ylwdown, yalb, ypphi(:,1), &
               ysolsw, ysollw, yts, ypplay(:,1), &
!!jyg               ycdragh, ycdragm, yrain_f, ysnow_f, yt(:,1), yq(:,1),&
               ycdragh, ycdragm, yrain_f, ysnow_f, yt1, yq1,&
               AcoefH, AcoefQ, BcoefH, BcoefQ, &
               AcoefU, AcoefV, BcoefU, BcoefV, &
               ypsref, yu1, yv1, ygustiness, yrugoro, pctsrf, &
               ysnow, yqsurf, yqsol, yagesno, &
               ytsoil, yz0m, yz0h, SFRWL, yalb_dir_new, yalb_dif_new, yevap,yfluxsens,yfluxlat, &
               ytsurf_new, y_dflux_t, y_dflux_q, &
               yzsig, ycldt, &
               ysnowhgt, yqsnow, ytoice, ysissnow, &
               yalb3_new, yrunoff, &
               y_flux_u1, y_flux_v1)

!jyg<
!!          alb3_lic(:)=0.
!>jyg
          DO j = 1, knon
             i = ni(j)
             alb3_lic(i) = yalb3_new(j)
             snowhgt(i)   = ysnowhgt(j)
             qsnow(i)     = yqsnow(j)
             to_ice(i)    = ytoice(j)
             sissnow(i)   = ysissnow(j)
             runoff(i)    = yrunoff(j)
          ENDDO
          ! Martin
! Special DICE MPL 05082013 puis BOMEX MPL 20150410
       IF (ok_prescr_ust) THEN
          DO j=1,knon
          y_flux_u1(j)=ycdragm(j)*(1.+sqrt(yu(j,1)*yu(j,1)+yv(j,1)*yv(j,1)))*yu(j,1)*ypplay(j,1)/RD/yt(j,1)
          y_flux_v1(j)=ycdragm(j)*(1.+sqrt(yu(j,1)*yu(j,1)+yv(j,1)*yv(j,1)))*yv(j,1)*ypplay(j,1)/RD/yt(j,1)
          ENDDO
      ENDIF
          
       CASE(is_oce)
           CALL surf_ocean(rlon, rlat, ysolsw, ysollw, yalb_vis, &
               ywindsp, rmu0, yfder, yts, &
               itap, dtime, jour, knon, ni, &
!!jyg               ypplay(:,1), zgeo1/RG, ycdragh, ycdragm, yrain_f, ysnow_f, yt(:,1), yq(:,1),&
               ypplay(:,1), zgeo1(1:knon)/RG, ycdragh, ycdragm, yrain_f, ysnow_f, yt(:,1), yq(:,1),&    ! ym missing init
               AcoefH, AcoefQ, BcoefH, BcoefQ, &
               AcoefU, AcoefV, BcoefU, BcoefV, &
               ypsref, yu1, yv1, ygustiness, yrugoro, pctsrf, &
               ysnow, yqsurf, yagesno, &
               yz0m, yz0h, SFRWL,yalb_dir_new, yalb_dif_new, yevap, yfluxsens,yfluxlat,&
               ytsurf_new, y_dflux_t, y_dflux_q, slab_wfbils, &
               y_flux_u1, y_flux_v1)
      IF (prt_level >=10) THEN
          print *,'arg de surf_ocean: ycdragh ',ycdragh
          print *,'arg de surf_ocean: ycdragm ',ycdragm
          print *,'arg de surf_ocean: yt ', yt
          print *,'arg de surf_ocean: yq ', yq
          print *,'arg de surf_ocean: yts ', yts
          print *,'arg de surf_ocean: AcoefH ',AcoefH
          print *,'arg de surf_ocean: AcoefQ ',AcoefQ
          print *,'arg de surf_ocean: BcoefH ',BcoefH
          print *,'arg de surf_ocean: BcoefQ ',BcoefQ
          print *,'arg de surf_ocean: yevap ',yevap
          print *,'arg de surf_ocean: yfluxsens ',yfluxsens
          print *,'arg de surf_ocean: yfluxlat ',yfluxlat
          print *,'arg de surf_ocean: ytsurf_new ',ytsurf_new
       ENDIF
! Special DICE MPL 05082013 puis BOMEX MPL 20150410
       IF (ok_prescr_ust) THEN
          DO j=1,knon
          y_flux_u1(j)=ycdragm(j)*(1.+sqrt(yu(j,1)*yu(j,1)+yv(j,1)*yv(j,1)))*yu(j,1)*ypplay(j,1)/RD/yt(j,1)
          y_flux_v1(j)=ycdragm(j)*(1.+sqrt(yu(j,1)*yu(j,1)+yv(j,1)*yv(j,1)))*yv(j,1)*ypplay(j,1)/RD/yt(j,1)
          ENDDO
      ENDIF
          
       CASE(is_sic)
          CALL surf_seaice( &
!albedo SB >>>
               rlon, rlat, ysolsw, ysollw, yalb_vis, yfder, &
!albedo SB <<<
               itap, dtime, jour, knon, ni, &
               lafin, &
!!jyg               yts, ypplay(:,1), ycdragh, ycdragm, yrain_f, ysnow_f, yt(:,1), yq(:,1),&
               yts, ypplay(:,1), ycdragh, ycdragm, yrain_f, ysnow_f, yt1, yq1,&
               AcoefH, AcoefQ, BcoefH, BcoefQ, &
               AcoefU, AcoefV, BcoefU, BcoefV, &
               ypsref, yu1, yv1, ygustiness, pctsrf, &
               ysnow, yqsurf, yqsol, yagesno, ytsoil, &
!albedo SB >>>
               yz0m, yz0h, SFRWL, yalb_dir_new, yalb_dif_new, yevap, yfluxsens,yfluxlat,&
!albedo SB <<<
               ytsurf_new, y_dflux_t, y_dflux_q, &
               y_flux_u1, y_flux_v1)
          
! Special DICE MPL 05082013 puis BOMEX MPL 20150410
       IF (ok_prescr_ust) THEN
          DO j=1,knon
          y_flux_u1(j)=ycdragm(j)*(1.+sqrt(yu(j,1)*yu(j,1)+yv(j,1)*yv(j,1)))*yu(j,1)*ypplay(j,1)/RD/yt(j,1)
          y_flux_v1(j)=ycdragm(j)*(1.+sqrt(yu(j,1)*yu(j,1)+yv(j,1)*yv(j,1)))*yv(j,1)*ypplay(j,1)/RD/yt(j,1)
          ENDDO
      ENDIF

       CASE DEFAULT
          WRITE(lunout,*) 'Surface index = ', nsrf
          abort_message = 'Surface index not valid'
          CALL abort_physic(modname,abort_message,1)
       END SELECT


!****************************************************************************************
! 11) - Calcul the increment of surface temperature
!
!****************************************************************************************

       IF (evap0>=0.) THEN
          yevap(:)=evap0
          yevap(:)=RLVTT*evap0
       ENDIF

       y_d_ts(1:knon)   = ytsurf_new(1:knon) - yts(1:knon)
 
!****************************************************************************************
!
! 12) "La remontee" - "The uphill"
!
!  The fluxes (y_flux_X) and tendancy (y_d_X) are calculated 
!  for X=H, Q, U and V, for all vertical levels.
!
!****************************************************************************************

!!!
!!! jyg le 10/04/2013
!!!
        IF (ok_flux_surf) THEN
          IF (prt_level >=10) THEN
           PRINT *,'pbl_surface: fsens flat RLVTT=',fsens,flat,RLVTT
          ENDIF
          y_flux_t1(:) =  fsens
          y_flux_q1(:) =  flat/RLVTT
          yfluxlat(:) =  flat
!
!!  Test sur iflag_split retire le 2/02/2018, sans vraiment comprendre la raison de ce test. (jyg)
!!          IF (iflag_split .eq.0) THEN
             DO j=1,knon
             Kech_h(j) = ycdragh(j) * (1.0+SQRT(yu(j,1)**2+yv(j,1)**2)) * &
                  ypplay(j,1)/(RD*yt(j,1))
             ENDDO
!!          ENDIF ! (iflag_split .eq.0)

          DO j = 1, knon
            yt1_new=(1./RCPD)*(AcoefH(j)+BcoefH(j)*y_flux_t1(j)*dtime)
            ytsurf_new(j)=yt1_new-y_flux_t1(j)/(Kech_h(j)*RCPD)
          ENDDO

          DO j=1,knon
          y_d_ts(j) = ytsurf_new(j) - yts(j) 
          ENDDO

        ELSE ! (ok_flux_surf)
          DO j=1,knon
          y_flux_t1(j) =  yfluxsens(j)
          y_flux_q1(j) = -yevap(j)
          ENDDO
        ENDIF

       IF (prt_level >=10) THEN
        DO j=1,knon
         print*,'y_flux_t1,yfluxlat,wakes' &
 &             ,  y_flux_t1(j), yfluxlat(j), ywake_s(j)
         print*,'beta,ytsurf_new', ybeta(j), ytsurf_new(j)
         print*,'inertia,facteur,cstar', inertia, facteur,wake_cstar(j)
        ENDDO
       ENDIF

!!! jyg le 07/02/2012 puis le 10/04/2013
!!       IF (iflag_split .eq.1) THEN
!!!!!
!!!jyg<
!!         CALL wx_pbl_split_no_dts(knon, ywake_s, &
!!                                AcoefH_x, AcoefH_w, &
!!                                AcoefQ_x, AcoefQ_w, &
!!                                AcoefU_x, AcoefU_w, &
!!                                AcoefV_x, AcoefV_w, &
!!                                y_flux_t1, y_flux_q1, y_flux_u1, y_flux_v1, &
!!                                y_flux_t1_x, y_flux_t1_w, &
!!                                y_flux_q1_x, y_flux_q1_w, &
!!                                y_flux_u1_x, y_flux_u1_w, &
!!                                y_flux_v1_x, y_flux_v1_w, &
!!                                yfluxlat_x, yfluxlat_w &
!!                                )
!!       ELSE IF (iflag_split .ge. 2) THEN
       IF (iflag_split .GE. 1) THEN
         CALL wx_pbl0_split(knon, dtime, ywake_s, &
                       y_flux_t1, y_flux_q1, y_flux_u1, y_flux_v1, &
                       y_flux_t1_x, y_flux_t1_w, &
                       y_flux_q1_x, y_flux_q1_w, &
                       y_flux_u1_x, y_flux_u1_w, &
                       y_flux_v1_x, y_flux_v1_w, &
                       yfluxlat_x, yfluxlat_w, &
                       y_delta_tsurf &
                       )
       ENDIF  ! (iflag_split .ge. 1)
!>jyg
!
 
!!jyg!!   A reprendre apres reflexion   ===============================================
!!jyg!!
!!jyg!!        DO j=1,knon
!!jyg!!!!! nrlmd le 13/06/2011
!!jyg!!
!!jyg!!!----Diffusion dans le sol dans le cas continental seulement
!!jyg!!       IF (nsrf.eq.is_ter) THEN
!!jyg!!!----Calcul du coefficient delta_coeff
!!jyg!!          tau_eq(j)=(ywake_s(j)/2.)*(1./max(wake_cstar(j),0.01))*sqrt(0.4/(3.14*max(wake_dens(j),8e-12)))
!!jyg!!
!!jyg!!!          delta_coef(j)=dtime/(inertia*sqrt(tau_eq(j)))
!!jyg!!          delta_coef(j)=facteur*sqrt(tau_eq(j))/inertia
!!jyg!!!          delta_coef(j)=0.
!!jyg!!       ELSE 
!!jyg!!         delta_coef(j)=0.
!!jyg!!       ENDIF
!!jyg!!
!!jyg!!!----Calcul de delta_tsurf
!!jyg!!         y_delta_tsurf(j)=delta_coef(j)*y_delta_flux_t1(j)
!!jyg!!
!!jyg!!!----Si il n'y a pas des poches...
!!jyg!!         IF (wake_cstar(j).le.0.01) THEN
!!jyg!!           y_delta_tsurf(j)=0.
!!jyg!!           y_delta_flux_t1(j)=0.
!!jyg!!         ENDIF
!!jyg!!
!!jyg!!!-----Calcul de ybeta (evap_r\'eelle/evap_potentielle)
!!jyg!!!!!!! jyg le 23/02/2012
!!jyg!!!!!!!
!!jyg!!!!        ybeta(j)=y_flux_q1(j)   /    &
!!jyg!!!! &        (Kech_h(j)*(yq(j,1)-yqsatsurf(j)))
!!jyg!!!!!!        ybeta(j)=-1.*yevap(j)   /    &
!!jyg!!!!!! &        (ywake_s(j)*Kech_h_w(j)*(yq_w(j,1)-yqsatsurf_w(j))+(1.-ywake_s(j))*Kech_h_x(j)*(yq_x(j,1)-yqsatsurf_x(j)))
!!jyg!!!!!!! fin jyg
!!jyg!!!!!
!!jyg!!
!!jyg!!       ENDDO
!!jyg!!
!!jyg!!!!! fin nrlmd le 13/06/2011
!!jyg!!
       IF (iflag_split .ge. 1) THEN
       IF (prt_level >=10) THEN
        DO j = 1, knon
         print*,'Chx,Chw,Ch', ycdragh_x(j), ycdragh_w(j), ycdragh(j)
         print*,'Khx,Khw,Kh', Kech_h_x(j), Kech_h_w(j), Kech_h(j)
!         print*,'tsurf_x,tsurf_w,tsurf,t1', ytsurf_th_x(j), ytsurf_th_w(j), ytsurf_th(j), yt(j,1)
         print*,'tsurf_x,t1x,tsurf_w,t1w,tsurf,t1,t1_ancien', &
 &               ytsurf_th_x(j), yt_x(j,1), ytsurf_th_w(j), yt_w(j,1), ytsurf_th(j), yt(j,1),t(j,1)
         print*,'qsatsurf,qsatsurf_x,qsatsurf_w', yqsatsurf(j), yqsatsurf_x(j), yqsatsurf_w(j)
         print*,'delta_coef,delta_flux,delta_tsurf,tau', delta_coef(j), y_delta_flux_t1(j), y_delta_tsurf(j), tau_eq(j) 
        ENDDO

        DO j=1,knon
         print*,'fluxT_x, fluxT_w, y_flux_t1, fluxQ_x, fluxQ_w, yfluxlat, wakes' &
 &             , y_flux_t1_x(j), y_flux_t1_w(j), y_flux_t1(j), y_flux_q1_x(j)*RLVTT, y_flux_q1_w(j)*RLVTT, yfluxlat(j), ywake_s(j)
         print*,'beta,ytsurf_new,yqsatsurf', ybeta(j), ytsurf_new(j), yqsatsurf(j)
         print*,'inertia,facteur,cstar', inertia, facteur,wake_cstar(j)
        ENDDO
       ENDIF  ! (prt_level >=10)

!!! jyg le 07/02/2012
       ENDIF  ! (iflag_split .ge.1)
!!!

!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
!!!
!!! nrlmd & jyg les 02/05/2011, 13/06/2011, 05/02/2012
        CALL climb_hq_up(knon, dtime, yt, yq, &
            y_flux_q1, y_flux_t1, ypaprs, ypplay, &
!!! jyg le 07/02/2012
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            CcoefH, CcoefQ, DcoefH, DcoefQ, &
            Kcoef_hq, gama_q, gama_h, &
!!!
            y_flux_q(:,:), y_flux_t(:,:), y_d_q(:,:), y_d_t(:,:))    
       ELSE  !(iflag_split .eq.0)
        CALL climb_hq_up(knon, dtime, yt_x, yq_x, &
            y_flux_q1_x, y_flux_t1_x, ypaprs, ypplay, &
!!! nrlmd le 02/05/2011
            AcoefH_x, AcoefQ_x, BcoefH_x, BcoefQ_x, &
            CcoefH_x, CcoefQ_x, DcoefH_x, DcoefQ_x, &
            Kcoef_hq_x, gama_q_x, gama_h_x, &
!!!
            y_flux_q_x(:,:), y_flux_t_x(:,:), y_d_q_x(:,:), y_d_t_x(:,:))    
!
       CALL climb_hq_up(knon, dtime, yt_w, yq_w, &
            y_flux_q1_w, y_flux_t1_w, ypaprs, ypplay, &
!!! nrlmd le 02/05/2011
            AcoefH_w, AcoefQ_w, BcoefH_w, BcoefQ_w, &
            CcoefH_w, CcoefQ_w, DcoefH_w, DcoefQ_w, &
            Kcoef_hq_w, gama_q_w, gama_h_w, &
!!!
            y_flux_q_w(:,:), y_flux_t_w(:,:), y_d_q_w(:,:), y_d_t_w(:,:))    
!!!
       ENDIF  ! (iflag_split .eq.0)
!!!

!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
!!!
!!! nrlmd & jyg les 02/05/2011, 13/06/2011, 05/02/2012
        CALL climb_wind_up(knon, dtime, yu, yv, y_flux_u1, y_flux_v1, &
!!! jyg le 07/02/2012
            AcoefU, AcoefV, BcoefU, BcoefV, &
            CcoefU, CcoefV, DcoefU, DcoefV, &
            Kcoef_m, &
!!!
            y_flux_u, y_flux_v, y_d_u, y_d_v)
     y_d_t_diss(:,:)=0.
     IF (iflag_pbl>=20 .and. iflag_pbl<30) THEN
        CALL yamada_c(knon,dtime,ypaprs,ypplay &
    &   ,yu,yv,yt,y_d_u,y_d_v,y_d_t,ycdragm,ytke,ycoefm,ycoefh,ycoefq,y_d_t_diss,yustar &
    &   ,iflag_pbl)
     ENDIF
!     print*,'yamada_c OK'

       ELSE  !(iflag_split .eq.0)
        CALL climb_wind_up(knon, dtime, yu_x, yv_x, y_flux_u1_x, y_flux_v1_x, &
!!! nrlmd le 02/05/2011
            AcoefU_x, AcoefV_x, BcoefU_x, BcoefV_x, &
            CcoefU_x, CcoefV_x, DcoefU_x, DcoefV_x, &
            Kcoef_m_x, &
!!!
            y_flux_u_x, y_flux_v_x, y_d_u_x, y_d_v_x)
!
     y_d_t_diss_x(:,:)=0.
     IF (iflag_pbl>=20 .and. iflag_pbl<30) THEN
        CALL yamada_c(knon,dtime,ypaprs,ypplay &
    &   ,yu_x,yv_x,yt_x,y_d_u_x,y_d_v_x,y_d_t_x,ycdragm_x,ytke_x,ycoefm_x,ycoefh_x &
        ,ycoefq_x,y_d_t_diss_x,yustar_x &
    &   ,iflag_pbl)
     ENDIF
!     print*,'yamada_c OK'

        CALL climb_wind_up(knon, dtime, yu_w, yv_w, y_flux_u1_w, y_flux_v1_w, &
!!! nrlmd le 02/05/2011
            AcoefU_w, AcoefV_w, BcoefU_w, BcoefV_w, &
            CcoefU_w, CcoefV_w, DcoefU_w, DcoefV_w, &
            Kcoef_m_w, &
!!!
            y_flux_u_w, y_flux_v_w, y_d_u_w, y_d_v_w)
!!!
     y_d_t_diss_w(:,:)=0.
     IF (iflag_pbl>=20 .and. iflag_pbl<30) THEN
        CALL yamada_c(knon,dtime,ypaprs,ypplay &
    &   ,yu_w,yv_w,yt_w,y_d_u_w,y_d_v_w,y_d_t_w,ycdragm_w,ytke_w,ycoefm_w,ycoefh_w &
        ,ycoefq_w,y_d_t_diss_w,yustar_w &
    &   ,iflag_pbl)
     ENDIF
!     print*,'yamada_c OK'
!
        IF (prt_level >=10) THEN
         print *, 'After climbing up, lfuxlat_x, fluxlat_w ', &
               yfluxlat_x, yfluxlat_w
        ENDIF
!
       ENDIF  ! (iflag_split .eq.0)
!!!

        DO j = 1, knon
          y_dflux_t(j) = y_dflux_t(j) * ypct(j)
          y_dflux_q(j) = y_dflux_q(j) * ypct(j)
        ENDDO

!****************************************************************************************
! 13) Transform variables for output format : 
!     - Decompress
!     - Multiply with pourcentage of current surface
!     - Cumulate in global variable
!
!****************************************************************************************


!!! jyg le 07/02/2012
       IF (iflag_split.EQ.0) THEN
!!!
        DO k = 1, klev
           DO j = 1, knon
             i = ni(j)
             y_d_t_diss(j,k)  = y_d_t_diss(j,k) * ypct(j)
             y_d_t(j,k)  = y_d_t(j,k) * ypct(j)
             y_d_q(j,k)  = y_d_q(j,k) * ypct(j)
             y_d_u(j,k)  = y_d_u(j,k) * ypct(j)
             y_d_v(j,k)  = y_d_v(j,k) * ypct(j)
!FC
             IF  (nsrf .EQ. is_ter .and. ifl_pbltree .GE. 1) THEN
!            if (y_d_u_frein(j,k).ne.0. ) then
!        print*, nsrf,'IS_TER ++', y_d_u_frein(j,k)*ypct(j),y_d_u(j,k),j,k
!            ENDIF
               y_d_u(j,k) =y_d_u(j,k) + y_d_u_frein(j,k)*ypct(j)
               y_d_v(j,k) =y_d_v(j,k) + y_d_v_frein(j,k)*ypct(j)
               treedrg(i,k,nsrf)=y_treedrg(j,k)
             ELSE 
               treedrg(i,k,nsrf)=0.
             ENDIF
!FC
             flux_t(i,k,nsrf) = y_flux_t(j,k)
             flux_q(i,k,nsrf) = y_flux_q(j,k)
             flux_u(i,k,nsrf) = y_flux_u(j,k)
             flux_v(i,k,nsrf) = y_flux_v(j,k)
           ENDDO
        ENDDO

       ELSE  !(iflag_split .eq.0)

! Tendances hors poches
        DO k = 1, klev
          DO j = 1, knon
            i = ni(j)
            y_d_t_diss_x(j,k)  = y_d_t_diss_x(j,k) * ypct(j)
            y_d_t_x(j,k)  = y_d_t_x(j,k) * ypct(j)
            y_d_q_x(j,k)  = y_d_q_x(j,k) * ypct(j)
            y_d_u_x(j,k)  = y_d_u_x(j,k) * ypct(j)
            y_d_v_x(j,k)  = y_d_v_x(j,k) * ypct(j)

            flux_t_x(i,k,nsrf) = y_flux_t_x(j,k)
            flux_q_x(i,k,nsrf) = y_flux_q_x(j,k)
            flux_u_x(i,k,nsrf) = y_flux_u_x(j,k)
            flux_v_x(i,k,nsrf) = y_flux_v_x(j,k)
          ENDDO
        ENDDO

! Tendances dans les poches
        DO k = 1, klev
          DO j = 1, knon
            i = ni(j)
            y_d_t_diss_w(j,k)  = y_d_t_diss_w(j,k) * ypct(j)
            y_d_t_w(j,k)  = y_d_t_w(j,k) * ypct(j)
            y_d_q_w(j,k)  = y_d_q_w(j,k) * ypct(j)
            y_d_u_w(j,k)  = y_d_u_w(j,k) * ypct(j)
            y_d_v_w(j,k)  = y_d_v_w(j,k) * ypct(j)

            flux_t_w(i,k,nsrf) = y_flux_t_w(j,k)
            flux_q_w(i,k,nsrf) = y_flux_q_w(j,k)
            flux_u_w(i,k,nsrf) = y_flux_u_w(j,k)
            flux_v_w(i,k,nsrf) = y_flux_v_w(j,k)
          ENDDO
        ENDDO

! Flux, tendances et Tke moyenne dans la maille
        DO k = 1, klev
          DO j = 1, knon
            i = ni(j)
            flux_t(i,k,nsrf) = flux_t_x(i,k,nsrf)+ywake_s(j)*(flux_t_w(i,k,nsrf)-flux_t_x(i,k,nsrf))
            flux_q(i,k,nsrf) = flux_q_x(i,k,nsrf)+ywake_s(j)*(flux_q_w(i,k,nsrf)-flux_q_x(i,k,nsrf))
            flux_u(i,k,nsrf) = flux_u_x(i,k,nsrf)+ywake_s(j)*(flux_u_w(i,k,nsrf)-flux_u_x(i,k,nsrf))
            flux_v(i,k,nsrf) = flux_v_x(i,k,nsrf)+ywake_s(j)*(flux_v_w(i,k,nsrf)-flux_v_x(i,k,nsrf))
          ENDDO
        ENDDO
        DO j=1,knon
          yfluxlat(j)=yfluxlat_x(j)+ywake_s(j)*(yfluxlat_w(j)-yfluxlat_x(j))
        ENDDO
        IF (prt_level >=10) THEN
          print *,' nsrf, flux_t(:,1,nsrf), flux_t_x(:,1,nsrf), flux_t_w(:,1,nsrf) ', &
                    nsrf, flux_t(:,1,nsrf), flux_t_x(:,1,nsrf), flux_t_w(:,1,nsrf)
        ENDIF

        DO k = 1, klev
          DO j = 1, knon
            y_d_t_diss(j,k) = y_d_t_diss_x(j,k)+ywake_s(j)*(y_d_t_diss_w(j,k) -y_d_t_diss_x(j,k))
            y_d_t(j,k) = y_d_t_x(j,k)+ywake_s(j)*(y_d_t_w(j,k) -y_d_t_x(j,k))
            y_d_q(j,k) = y_d_q_x(j,k)+ywake_s(j)*(y_d_q_w(j,k) -y_d_q_x(j,k))
            y_d_u(j,k) = y_d_u_x(j,k)+ywake_s(j)*(y_d_u_w(j,k) -y_d_u_x(j,k))
            y_d_v(j,k) = y_d_v_x(j,k)+ywake_s(j)*(y_d_v_w(j,k) -y_d_v_x(j,k))
          ENDDO
        ENDDO

       ENDIF  ! (iflag_split .eq.0)
!!!

!      print*,'Dans pbl OK1'

!jyg<
!!       evap(:,nsrf) = - flux_q(:,1,nsrf)
!>jyg
       DO j = 1, knon
          i = ni(j)
          evap(i,nsrf) = - flux_q(i,1,nsrf)                  !jyg
          d_ts(i,nsrf) = y_d_ts(j)
!albedo SB >>>
          DO k=1,nsw
            alb_dir(i,k,nsrf) = yalb_dir_new(j,k)
            alb_dif(i,k,nsrf) = yalb_dif_new(j,k)
          ENDDO
!albedo SB <<<
          snow(i,nsrf) = ysnow(j)  
          qsurf(i,nsrf) = yqsurf(j)
          z0m(i,nsrf) = yz0m(j)
          z0h(i,nsrf) = yz0h(j)
          fluxlat(i,nsrf) = yfluxlat(j)
          agesno(i,nsrf) = yagesno(j)  
          cdragh(i) = cdragh(i) + ycdragh(j)*ypct(j)
          cdragm(i) = cdragm(i) + ycdragm(j)*ypct(j)
          dflux_t(i) = dflux_t(i) + y_dflux_t(j)
          dflux_q(i) = dflux_q(i) + y_dflux_q(j)
       ENDDO

!      print*,'Dans pbl OK2'

!!! jyg le 07/02/2012
       IF (iflag_split .ge.1) THEN
!!!
!!! nrlmd le 02/05/2011
        DO j = 1, knon
          i = ni(j)
          fluxlat_x(i,nsrf) = yfluxlat_x(j)
          fluxlat_w(i,nsrf) = yfluxlat_w(j)
!!!
!!! nrlmd le 13/06/2011
!!jyg20170131          delta_tsurf(i,nsrf)=y_delta_tsurf(j)*ypct(j)
          delta_tsurf(i,nsrf)=y_delta_tsurf(j)
!
          cdragh_x(i) = cdragh_x(i) + ycdragh_x(j)*ypct(j)
          cdragh_w(i) = cdragh_w(i) + ycdragh_w(j)*ypct(j)
          cdragm_x(i) = cdragm_x(i) + ycdragm_x(j)*ypct(j)
          cdragm_w(i) = cdragm_w(i) + ycdragm_w(j)*ypct(j)
          kh(i) = kh(i) + Kech_h(j)*ypct(j)
          kh_x(i) = kh_x(i) + Kech_h_x(j)*ypct(j)
          kh_w(i) = kh_w(i) + Kech_h_w(j)*ypct(j)
!!!
        ENDDO
!!!      
       ENDIF  ! (iflag_split .ge.1)
!!!
!!! nrlmd le 02/05/2011
!!jyg le 20/02/2011
!!        tke_x(:,:,nsrf)=0.
!!        tke_w(:,:,nsrf)=0.
!!jyg le 20/02/2011
!!        DO k = 1, klev+1
!!          DO j = 1, knon
!!            i = ni(j)
!!            wake_dltke(i,k,nsrf) = ytke_w(j,k) - ytke_x(j,k)
!!            tke(i,k,nsrf)   = ytke_x(j,k) + ywake_s(j)*wake_dltke(i,k,nsrf)
!!          ENDDO
!!        ENDDO
!!jyg le 20/02/2011
!!        DO k = 1, klev+1
!!          DO j = 1, knon
!!            i = ni(j)
!!            tke(i,k,nsrf)=(1.-ywake_s(j))*tke_x(i,k,nsrf)+ywake_s(j)*tke_w(i,k,nsrf)
!!          ENDDO
!!        ENDDO
!!!
       IF (iflag_split .eq.0) THEN
        wake_dltke(:,:,nsrf) = 0.
        DO k = 1, klev+1
           DO j = 1, knon
              i = ni(j)
!jyg<
!!              tke(i,k,nsrf)    = ytke(j,k)
!!              tke(i,k,is_ave) = tke(i,k,is_ave) + ytke(j,k)*ypct(j)
              tke_x(i,k,nsrf)    = ytke(j,k)
              tke_x(i,k,is_ave) = tke_x(i,k,is_ave) + ytke(j,k)*ypct(j)
!>jyg
           ENDDO
        ENDDO

       ELSE  ! (iflag_split .eq.0)
        DO k = 1, klev+1
          DO j = 1, knon
            i = ni(j)
            wake_dltke(i,k,nsrf) = ytke_w(j,k) - ytke_x(j,k)
!jyg<
!!            tke(i,k,nsrf)   = ytke_x(j,k) + ywake_s(j)*wake_dltke(i,k,nsrf)
!!            tke(i,k,is_ave) = tke(i,k,is_ave) + tke(i,k,nsrf)*ypct(j)
            tke_x(i,k,nsrf)   = ytke_x(j,k)
            tke_x(i,k,is_ave)   = tke_x(i,k,is_ave) + tke_x(i,k,nsrf)*ypct(j)
            wake_dltke(i,k,is_ave)   = wake_dltke(i,k,is_ave) + wake_dltke(i,k,nsrf)*ypct(j)

!>jyg
          ENDDO
        ENDDO
       ENDIF  ! (iflag_split .eq.0)
!!!
       DO k = 2, klev
          DO j = 1, knon
             i = ni(j)
             zcoefh(i,k,nsrf) = ycoefh(j,k)
             zcoefm(i,k,nsrf) = ycoefm(j,k)
             zcoefh(i,k,is_ave) = zcoefh(i,k,is_ave) + ycoefh(j,k)*ypct(j)
             zcoefm(i,k,is_ave) = zcoefm(i,k,is_ave) + ycoefm(j,k)*ypct(j)
          ENDDO
       ENDDO

!      print*,'Dans pbl OK3'

       IF ( nsrf .EQ. is_ter ) THEN 
          DO j = 1, knon
             i = ni(j)
             qsol(i) = yqsol(j)
          ENDDO
       ENDIF
       
!jyg<
!!       ftsoil(:,:,nsrf) = 0.
!>jyg
       DO k = 1, nsoilmx
          DO j = 1, knon
             i = ni(j)
             ftsoil(i, k, nsrf) = ytsoil(j,k)
          ENDDO
       ENDDO
       
!!! jyg le 07/02/2012
       IF (iflag_split .ge.1) THEN
!!!
!!! nrlmd+jyg le 02/05/2011 et le 20/02/2012
        DO k = 1, klev
          DO j = 1, knon
           i = ni(j)
           d_t_diss_x(i,k) = d_t_diss_x(i,k) + y_d_t_diss_x(j,k)
           d_t_x(i,k) = d_t_x(i,k) + y_d_t_x(j,k)
           d_q_x(i,k) = d_q_x(i,k) + y_d_q_x(j,k)
           d_u_x(i,k) = d_u_x(i,k) + y_d_u_x(j,k)
           d_v_x(i,k) = d_v_x(i,k) + y_d_v_x(j,k)
!
           d_t_diss_w(i,k) = d_t_diss_w(i,k) + y_d_t_diss_w(j,k)
           d_t_w(i,k) = d_t_w(i,k) + y_d_t_w(j,k)
           d_q_w(i,k) = d_q_w(i,k) + y_d_q_w(j,k)
           d_u_w(i,k) = d_u_w(i,k) + y_d_u_w(j,k)
           d_v_w(i,k) = d_v_w(i,k) + y_d_v_w(j,k)
!
!!           d_wake_dlt(i,k) = d_wake_dlt(i,k) + y_d_t_w(i,k)-y_d_t_x(i,k)
!!           d_wake_dlq(i,k) = d_wake_dlq(i,k) + y_d_q_w(i,k)-y_d_q_x(i,k)
          ENDDO
        ENDDO
!!!
       ENDIF  ! (iflag_split .ge.1)
!!!
       
       DO k = 1, klev
          DO j = 1, knon
             i = ni(j)
             d_t_diss(i,k) = d_t_diss(i,k) + y_d_t_diss(j,k)
             d_t(i,k) = d_t(i,k) + y_d_t(j,k)
             d_q(i,k) = d_q(i,k) + y_d_q(j,k)
             d_u(i,k) = d_u(i,k) + y_d_u(j,k)
             d_v(i,k) = d_v(i,k) + y_d_v(j,k)
          ENDDO
       ENDDO

!      print*,'Dans pbl OK4'

       IF (prt_level >=10) THEN
         print *, 'pbl_surface tendencies for w: d_t_w, d_t_x, d_t ', &
          d_t_w(:,1), d_t_x(:,1), d_t(:,1)
       ENDIF

!****************************************************************************************
! 14) Calculate the temperature and relative humidity at 2m and the wind at 10m 
!     Call HBTM
!
!****************************************************************************************
!!!
!
#undef T2m     
#define T2m     
#ifdef T2m
! Calculations of diagnostic t,q at 2m and u, v at 10m

!      print*,'Dans pbl OK41'
!      print*,'tair1,yt(:,1),y_d_t(:,1)'
!      print*, tair1,yt(:,1),y_d_t(:,1)
!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
        DO j=1, knon
          uzon(j) = yu(j,1) + y_d_u(j,1)
          vmer(j) = yv(j,1) + y_d_v(j,1)
          tair1(j) = yt(j,1) + y_d_t(j,1) + y_d_t_diss(j,1)
          qair1(j) = yq(j,1) + y_d_q(j,1)
          zgeo1(j) = RD * tair1(j) / (0.5*(ypaprs(j,1)+ypplay(j,1))) &
               * (ypaprs(j,1)-ypplay(j,1))
          tairsol(j) = yts(j) + y_d_ts(j)
          qairsol(j) = yqsurf(j)
        ENDDO
       ELSE  ! (iflag_split .eq.0)
        DO j=1, knon
          uzon_x(j) = yu_x(j,1) + y_d_u_x(j,1)
          vmer_x(j) = yv_x(j,1) + y_d_v_x(j,1)
          tair1_x(j) = yt_x(j,1) + y_d_t_x(j,1) + y_d_t_diss_x(j,1)
          qair1_x(j) = yq_x(j,1) + y_d_q_x(j,1)
          zgeo1_x(j) = RD * tair1_x(j) / (0.5*(ypaprs(j,1)+ypplay(j,1))) &
               * (ypaprs(j,1)-ypplay(j,1))
          tairsol(j) = yts(j) + y_d_ts(j)
          tairsol_x(j) = tairsol(j) - ywake_s(j)*y_delta_tsurf(j)
          qairsol(j) = yqsurf(j)
        ENDDO
        DO j=1, knon
          uzon_w(j) = yu_w(j,1) + y_d_u_w(j,1)
          vmer_w(j) = yv_w(j,1) + y_d_v_w(j,1)
          tair1_w(j) = yt_w(j,1) + y_d_t_w(j,1) + y_d_t_diss_w(j,1)
          qair1_w(j) = yq_w(j,1) + y_d_q_w(j,1)
          zgeo1_w(j) = RD * tair1_w(j) / (0.5*(ypaprs(j,1)+ypplay(j,1))) &
               * (ypaprs(j,1)-ypplay(j,1))
          tairsol_w(j) = tairsol(j) + (1.- ywake_s(j))*y_delta_tsurf(j)
          qairsol(j) = yqsurf(j)
        ENDDO
!!!      
       ENDIF  ! (iflag_split .eq.0)
!!!
       DO j=1, knon
!         i = ni(j)
!         yz0h_oupas(j) = yz0m(j)
!         IF(nsrf.EQ.is_oce) THEN
!            yz0h_oupas(j) = z0m(i,nsrf)
!         ENDIF
          psfce(j)=ypaprs(j,1)
          patm(j)=ypplay(j,1)
       ENDDO

       IF (iflag_pbl_surface_t2m_bug==1) THEN
          yz0h_oupas(1:knon)=yz0m(1:knon)
       ELSE
          yz0h_oupas(1:knon)=yz0h(1:knon)
       ENDIF
       
!      print*,'Dans pbl OK42A'
!      print*,'tair1,yt(:,1),y_d_t(:,1)'
!      print*, tair1,yt(:,1),y_d_t(:,1)

! Calculate the temperature and relative humidity at 2m and the wind at 10m 
!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
        IF (iflag_new_t2mq2m==1) THEN
         CALL stdlevvarn(klon, knon, nsrf, zxli, &
            uzon, vmer, tair1, qair1, zgeo1, &
            tairsol, qairsol, yz0m, yz0h_oupas, psfce, patm, &
            yt2m, yq2m, yt10m, yq10m, yu10m, yustar, &
            yn2mout(:, nsrf, :))
        ELSE
        CALL stdlevvar(klon, knon, nsrf, zxli, &
            uzon, vmer, tair1, qair1, zgeo1, &
            tairsol, qairsol, yz0m, yz0h_oupas, psfce, patm, &
            yt2m, yq2m, yt10m, yq10m, yu10m, yustar)
        ENDIF
       ELSE  !(iflag_split .eq.0)
        IF (iflag_new_t2mq2m==1) THEN
         CALL stdlevvarn(klon, knon, nsrf, zxli, &
            uzon_x, vmer_x, tair1_x, qair1_x, zgeo1_x, &
            tairsol_x, qairsol, yz0m, yz0h_oupas, psfce, patm, &
            yt2m_x, yq2m_x, yt10m_x, yq10m_x, yu10m_x, yustar_x, &
            yn2mout_x(:, nsrf, :))
         CALL stdlevvarn(klon, knon, nsrf, zxli, &
            uzon_w, vmer_w, tair1_w, qair1_w, zgeo1_w, &
            tairsol_w, qairsol, yz0m, yz0h_oupas, psfce, patm, &
            yt2m_w, yq2m_w, yt10m_w, yq10m_w, yu10m_w, yustar_w, &
            yn2mout_w(:, nsrf, :))
        ELSE
        CALL stdlevvar(klon, knon, nsrf, zxli, &
            uzon_x, vmer_x, tair1_x, qair1_x, zgeo1_x, &
            tairsol_x, qairsol, yz0m, yz0h_oupas, psfce, patm, &
            yt2m_x, yq2m_x, yt10m_x, yq10m_x, yu10m_x, yustar_x)
        CALL stdlevvar(klon, knon, nsrf, zxli, &
            uzon_w, vmer_w, tair1_w, qair1_w, zgeo1_w, &
            tairsol_w, qairsol, yz0m, yz0h_oupas, psfce, patm, &
            yt2m_w, yq2m_w, yt10m_w, yq10m_w, yu10m_w, yustar_w)
        ENDIF
!!!
       ENDIF  ! (iflag_split .eq.0)
!!!
!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
        DO j=1, knon
          i = ni(j)
          t2m(i,nsrf)=yt2m(j)
          q2m(i,nsrf)=yq2m(j)
     ! u10m, v10m : composantes du vent a 10m sans spirale de Ekman
          ustar(i,nsrf)=yustar(j)
          u10m(i,nsrf)=(yu10m(j) * uzon(j))/SQRT(uzon(j)**2+vmer(j)**2)
          v10m(i,nsrf)=(yu10m(j) * vmer(j))/SQRT(uzon(j)**2+vmer(j)**2)
!
          DO k = 1, 6
           n2mout(i,nsrf,k) = yn2mout(j,nsrf,k)
          END DO
!
        ENDDO
       ELSE  !(iflag_split .eq.0)
        DO j=1, knon
          i = ni(j)
          t2m_x(i,nsrf)=yt2m_x(j)
          q2m_x(i,nsrf)=yq2m_x(j)
     ! u10m, v10m : composantes du vent a 10m sans spirale de Ekman
          ustar_x(i,nsrf)=yustar_x(j)
          u10m_x(i,nsrf)=(yu10m_x(j) * uzon_x(j))/SQRT(uzon_x(j)**2+vmer_x(j)**2)
          v10m_x(i,nsrf)=(yu10m_x(j) * vmer_x(j))/SQRT(uzon_x(j)**2+vmer_x(j)**2)
!
          DO k = 1, 6
           n2mout_x(i,nsrf,k) = yn2mout_x(j,nsrf,k)
          END DO
!
        ENDDO
        DO j=1, knon
          i = ni(j)
          t2m_w(i,nsrf)=yt2m_w(j)
          q2m_w(i,nsrf)=yq2m_w(j)
     ! u10m, v10m : composantes du vent a 10m sans spirale de Ekman
          ustar_w(i,nsrf)=yustar_w(j)
          u10m_w(i,nsrf)=(yu10m_w(j) * uzon_w(j))/SQRT(uzon_w(j)**2+vmer_w(j)**2)
          v10m_w(i,nsrf)=(yu10m_w(j) * vmer_w(j))/SQRT(uzon_w(j)**2+vmer_w(j)**2)
!
          ustar(i,nsrf) = ustar_x(i,nsrf) + wake_s(i)*(ustar_w(i,nsrf)-ustar_x(i,nsrf))
          u10m(i,nsrf) = u10m_x(i,nsrf) + wake_s(i)*(u10m_w(i,nsrf)-u10m_x(i,nsrf))
          v10m(i,nsrf) = v10m_x(i,nsrf) + wake_s(i)*(v10m_w(i,nsrf)-v10m_x(i,nsrf))
!
          DO k = 1, 6
           n2mout_w(i,nsrf,k) = yn2mout_w(j,nsrf,k)
          END DO
!
        ENDDO
!!!
       ENDIF  ! (iflag_split .eq.0)
!!!

!      print*,'Dans pbl OK43'
!IM Calcule de l'humidite relative a 2m (rh2m) pour diagnostique
!IM Ajoute dependance type surface
       IF (thermcep) THEN
!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
          DO j = 1, knon
             i=ni(j)
             zdelta1 = MAX(0.,SIGN(1., rtt-yt2m(j) ))
             zx_qs1  = r2es * FOEEW(yt2m(j),zdelta1)/paprs(i,1)
             zx_qs1  = MIN(0.5,zx_qs1)
             zcor1   = 1./(1.-RETV*zx_qs1)
             zx_qs1  = zx_qs1*zcor1
             
             rh2m(i)   = rh2m(i)   + yq2m(j)/zx_qs1 * pctsrf(i,nsrf)
             qsat2m(i) = qsat2m(i) + zx_qs1  * pctsrf(i,nsrf)
          ENDDO
       ELSE  ! (iflag_split .eq.0)
          DO j = 1, knon
             i=ni(j)
             zdelta1 = MAX(0.,SIGN(1., rtt-yt2m_x(j) ))
             zx_qs1  = r2es * FOEEW(yt2m_x(j),zdelta1)/paprs(i,1)
             zx_qs1  = MIN(0.5,zx_qs1)
             zcor1   = 1./(1.-RETV*zx_qs1)
             zx_qs1  = zx_qs1*zcor1
             
             rh2m_x(i)   = rh2m_x(i)   + yq2m_x(j)/zx_qs1 * pctsrf(i,nsrf)
             qsat2m_x(i) = qsat2m_x(i) + zx_qs1  * pctsrf(i,nsrf)
          ENDDO
          DO j = 1, knon
             i=ni(j)
             zdelta1 = MAX(0.,SIGN(1., rtt-yt2m_w(j) ))
             zx_qs1  = r2es * FOEEW(yt2m_w(j),zdelta1)/paprs(i,1)
             zx_qs1  = MIN(0.5,zx_qs1)
             zcor1   = 1./(1.-RETV*zx_qs1)
             zx_qs1  = zx_qs1*zcor1
             
             rh2m_w(i)   = rh2m_w(i)   + yq2m_w(j)/zx_qs1 * pctsrf(i,nsrf)
             qsat2m_w(i) = qsat2m_w(i) + zx_qs1  * pctsrf(i,nsrf)
          ENDDO
!!!      
       ENDIF  ! (iflag_split .eq.0)
!!!
       ENDIF
!
       IF (prt_level >=10) THEN
         print *, 'T2m, q2m, RH2m ', &
          t2m, q2m, rh2m
       ENDIF

!   print*,'OK pbl 5'
!
!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
        CALL hbtm(knon, ypaprs, ypplay, &
            yt2m,yt10m,yq2m,yq10m,yustar,ywstar, &
            y_flux_t,y_flux_q,yu,yv,yt,yq, &
            ypblh,ycapCL,yoliqCL,ycteiCL,ypblT, &
            ytherm,ytrmb1,ytrmb2,ytrmb3,ylcl)
          IF (prt_level >=10) THEN
       print *,' Arg. de HBTM: yt2m ',yt2m
       print *,' Arg. de HBTM: yt10m ',yt10m
       print *,' Arg. de HBTM: yq2m ',yq2m
       print *,' Arg. de HBTM: yq10m ',yq10m
       print *,' Arg. de HBTM: yustar ',yustar
       print *,' Arg. de HBTM: y_flux_t ',y_flux_t
       print *,' Arg. de HBTM: y_flux_q ',y_flux_q
       print *,' Arg. de HBTM: yu ',yu
       print *,' Arg. de HBTM: yv ',yv
       print *,' Arg. de HBTM: yt ',yt
       print *,' Arg. de HBTM: yq ',yq
          ENDIF
       ELSE  ! (iflag_split .eq.0)
        CALL HBTM(knon, ypaprs, ypplay, &
            yt2m_x,yt10m_x,yq2m_x,yq10m_x,yustar_x,ywstar_x, &
            y_flux_t_x,y_flux_q_x,yu_x,yv_x,yt_x,yq_x, &
            ypblh_x,ycapCL_x,yoliqCL_x,ycteiCL_x,ypblT_x, &
            ytherm_x,ytrmb1_x,ytrmb2_x,ytrmb3_x,ylcl_x)
          IF (prt_level >=10) THEN
       print *,' Arg. de HBTM: yt2m_x ',yt2m_x
       print *,' Arg. de HBTM: yt10m_x ',yt10m_x
       print *,' Arg. de HBTM: yq2m_x ',yq2m_x
       print *,' Arg. de HBTM: yq10m_x ',yq10m_x
       print *,' Arg. de HBTM: yustar_x ',yustar_x
       print *,' Arg. de HBTM: y_flux_t_x ',y_flux_t_x
       print *,' Arg. de HBTM: y_flux_q_x ',y_flux_q_x
       print *,' Arg. de HBTM: yu_x ',yu_x
       print *,' Arg. de HBTM: yv_x ',yv_x
       print *,' Arg. de HBTM: yt_x ',yt_x
       print *,' Arg. de HBTM: yq_x ',yq_x
          ENDIF
        CALL HBTM(knon, ypaprs, ypplay, &
            yt2m_w,yt10m_w,yq2m_w,yq10m_w,yustar_w,ywstar_w, &
            y_flux_t_w,y_flux_q_w,yu_w,yv_w,yt_w,yq_w, &
            ypblh_w,ycapCL_w,yoliqCL_w,ycteiCL_w,ypblT_w, &
            ytherm_w,ytrmb1_w,ytrmb2_w,ytrmb3_w,ylcl_w)
!!!      
       ENDIF  ! (iflag_split .eq.0)
!!!
       
!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
!!!
        DO j=1, knon
          i = ni(j)
          pblh(i,nsrf)   = ypblh(j)
          wstar(i,nsrf)  = ywstar(j)
          plcl(i,nsrf)   = ylcl(j)
          capCL(i,nsrf)  = ycapCL(j)
          oliqCL(i,nsrf) = yoliqCL(j)
          cteiCL(i,nsrf) = ycteiCL(j)
          pblT(i,nsrf)   = ypblT(j)
          therm(i,nsrf)  = ytherm(j)
          trmb1(i,nsrf)  = ytrmb1(j)
          trmb2(i,nsrf)  = ytrmb2(j)
          trmb3(i,nsrf)  = ytrmb3(j)
        ENDDO
        IF (prt_level >=10) THEN
          print *, 'After HBTM: pblh ', pblh
          print *, 'After HBTM: plcl ', plcl
          print *, 'After HBTM: cteiCL ', cteiCL
        ENDIF
       ELSE  !(iflag_split .eq.0)
        DO j=1, knon
          i = ni(j)
          pblh_x(i,nsrf)   = ypblh_x(j)
          wstar_x(i,nsrf)  = ywstar_x(j)
          plcl_x(i,nsrf)   = ylcl_x(j)
          capCL_x(i,nsrf)  = ycapCL_x(j)
          oliqCL_x(i,nsrf) = yoliqCL_x(j)
          cteiCL_x(i,nsrf) = ycteiCL_x(j)
          pblT_x(i,nsrf)   = ypblT_x(j)
          therm_x(i,nsrf)  = ytherm_x(j)
          trmb1_x(i,nsrf)  = ytrmb1_x(j)
          trmb2_x(i,nsrf)  = ytrmb2_x(j)
          trmb3_x(i,nsrf)  = ytrmb3_x(j)
        ENDDO
        IF (prt_level >=10) THEN
          print *, 'After HBTM: pblh_x ', pblh_x
          print *, 'After HBTM: plcl_x ', plcl_x
          print *, 'After HBTM: cteiCL_x ', cteiCL_x
        ENDIF
        DO j=1, knon
          i = ni(j)
          pblh_w(i,nsrf)   = ypblh_w(j)
          wstar_w(i,nsrf)  = ywstar_w(j)
          plcl_w(i,nsrf)   = ylcl_w(j)
          capCL_w(i,nsrf)  = ycapCL_w(j)
          oliqCL_w(i,nsrf) = yoliqCL_w(j)
          cteiCL_w(i,nsrf) = ycteiCL_w(j)
          pblT_w(i,nsrf)   = ypblT_w(j)
          therm_w(i,nsrf)  = ytherm_w(j)
          trmb1_w(i,nsrf)  = ytrmb1_w(j)
          trmb2_w(i,nsrf)  = ytrmb2_w(j)
          trmb3_w(i,nsrf)  = ytrmb3_w(j)
        ENDDO
        IF (prt_level >=10) THEN
          print *, 'After HBTM: pblh_w ', pblh_w
          print *, 'After HBTM: plcl_w ', plcl_w
          print *, 'After HBTM: cteiCL_w ', cteiCL_w
        ENDIF
!!!
       ENDIF  ! (iflag_split .eq.0)
!!!

!   print*,'OK pbl 6'
#else 
! T2m not defined
! No calculation
       PRINT*,' Warning !!! No T2m calculation. Output is set to zero.'
#endif

!****************************************************************************************
! 15) End of loop over different surfaces
!
!****************************************************************************************
    ENDDO loop_nbsrf

!****************************************************************************************
! 16) Calculate the mean value over all sub-surfaces for some variables
!
!****************************************************************************************
    
    z0m(:,nbsrf+1) = 0.0
    z0h(:,nbsrf+1) = 0.0
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          z0m(i,nbsrf+1) = z0m(i,nbsrf+1) + z0m(i,nsrf)*pctsrf(i,nsrf)
          z0h(i,nbsrf+1) = z0h(i,nbsrf+1) + z0h(i,nsrf)*pctsrf(i,nsrf)
       ENDDO
    ENDDO

!   print*,'OK pbl 7'
    zxfluxt(:,:) = 0.0 ; zxfluxq(:,:) = 0.0
    zxfluxu(:,:) = 0.0 ; zxfluxv(:,:) = 0.0
    zxfluxt_x(:,:) = 0.0 ; zxfluxq_x(:,:) = 0.0
    zxfluxu_x(:,:) = 0.0 ; zxfluxv_x(:,:) = 0.0
    zxfluxt_w(:,:) = 0.0 ; zxfluxq_w(:,:) = 0.0
    zxfluxu_w(:,:) = 0.0 ; zxfluxv_w(:,:) = 0.0

!!! jyg le 07/02/2012
       IF (iflag_split .ge.1) THEN
!!!
!!! nrlmd & jyg les 02/05/2011, 05/02/2012

        DO nsrf = 1, nbsrf
          DO k = 1, klev
            DO i = 1, klon
              zxfluxt_x(i,k) = zxfluxt_x(i,k) + flux_t_x(i,k,nsrf) * pctsrf(i,nsrf)
              zxfluxq_x(i,k) = zxfluxq_x(i,k) + flux_q_x(i,k,nsrf) * pctsrf(i,nsrf)
              zxfluxu_x(i,k) = zxfluxu_x(i,k) + flux_u_x(i,k,nsrf) * pctsrf(i,nsrf)
              zxfluxv_x(i,k) = zxfluxv_x(i,k) + flux_v_x(i,k,nsrf) * pctsrf(i,nsrf)
!
              zxfluxt_w(i,k) = zxfluxt_w(i,k) + flux_t_w(i,k,nsrf) * pctsrf(i,nsrf)
              zxfluxq_w(i,k) = zxfluxq_w(i,k) + flux_q_w(i,k,nsrf) * pctsrf(i,nsrf)
              zxfluxu_w(i,k) = zxfluxu_w(i,k) + flux_u_w(i,k,nsrf) * pctsrf(i,nsrf)
              zxfluxv_w(i,k) = zxfluxv_w(i,k) + flux_v_w(i,k,nsrf) * pctsrf(i,nsrf)
            ENDDO
          ENDDO
        ENDDO

    DO i = 1, klon
      zxsens_x(i) = - zxfluxt_x(i,1)
      zxsens_w(i) = - zxfluxt_w(i,1)
    ENDDO
!!!
       ENDIF  ! (iflag_split .ge.1)
!!!

    DO nsrf = 1, nbsrf
       DO k = 1, klev
          DO i = 1, klon
             zxfluxt(i,k) = zxfluxt(i,k) + flux_t(i,k,nsrf) * pctsrf(i,nsrf)
             zxfluxq(i,k) = zxfluxq(i,k) + flux_q(i,k,nsrf) * pctsrf(i,nsrf)
             zxfluxu(i,k) = zxfluxu(i,k) + flux_u(i,k,nsrf) * pctsrf(i,nsrf)
             zxfluxv(i,k) = zxfluxv(i,k) + flux_v(i,k,nsrf) * pctsrf(i,nsrf)
          ENDDO
       ENDDO
    ENDDO

    DO i = 1, klon
       zxsens(i)     = - zxfluxt(i,1) ! flux de chaleur sensible au sol
       zxevap(i)     = - zxfluxq(i,1) ! flux d'evaporation au sol
       fder_print(i) = fder(i) + dflux_t(i) + dflux_q(i)
    ENDDO
!!!
   
!
! Incrementer la temperature du sol
!
    zxtsol(:) = 0.0  ; zxfluxlat(:) = 0.0
    zt2m(:) = 0.0    ; zq2m(:) = 0.0 ; zn2mout(:,:) = 0
    zustar(:)=0.0 ; zu10m(:) = 0.0   ; zv10m(:) = 0.0
    s_pblh(:) = 0.0  ; s_plcl(:) = 0.0 
!!! jyg le 07/02/2012
     s_pblh_x(:) = 0.0  ; s_plcl_x(:) = 0.0 
     s_pblh_w(:) = 0.0  ; s_plcl_w(:) = 0.0 
!!!
    s_capCL(:) = 0.0 ; s_oliqCL(:) = 0.0
    s_cteiCL(:) = 0.0; s_pblT(:) = 0.0
    s_therm(:) = 0.0 ; s_trmb1(:) = 0.0
    s_trmb2(:) = 0.0 ; s_trmb3(:) = 0.0
    wstar(:,is_ave)=0.
    
!   print*,'OK pbl 9'
    
!!! nrlmd le 02/05/2011
    zxfluxlat_x(:) = 0.0  ;  zxfluxlat_w(:) = 0.0
!!!
    
    DO nsrf = 1, nbsrf
       DO i = 1, klon          
          ts(i,nsrf) = ts(i,nsrf) + d_ts(i,nsrf)
          
          wfbils(i,nsrf) = ( solsw(i,nsrf) + sollw(i,nsrf) &
               + flux_t(i,1,nsrf) + fluxlat(i,nsrf) ) * pctsrf(i,nsrf)
          wfbilo(i,nsrf) = (evap(i,nsrf)-(rain_f(i)+snow_f(i)))*pctsrf(i,nsrf)
          wfevap(i,nsrf) = evap(i,nsrf)*pctsrf(i,nsrf)
          wfrain(i,nsrf) = rain_f(i)*pctsrf(i,nsrf)
          wfsnow(i,nsrf) = snow_f(i)*pctsrf(i,nsrf)

          zxtsol(i)    = zxtsol(i)    + ts(i,nsrf)      * pctsrf(i,nsrf)
          zxfluxlat(i) = zxfluxlat(i) + fluxlat(i,nsrf) * pctsrf(i,nsrf)
       ENDDO
    ENDDO
!
!<al1 order 2 correction to zxtsol, for radiation computations (main atm effect of Ts)
   IF (iflag_order2_sollw == 1) THEN
    meansqT(:) = 0. ! as working buffer
    DO nsrf = 1, nbsrf
     DO i = 1, klon
      meansqT(i) = meansqT(i)+(ts(i,nsrf)-zxtsol(i))**2 *pctsrf(i,nsrf)
     ENDDO
    ENDDO
    zxtsol(:) = zxtsol(:)+1.5*meansqT(:)/zxtsol(:)
   ENDIF   ! iflag_order2_sollw == 1
!>al1
          
!!! jyg le 07/02/2012
       IF (iflag_split .eq.0) THEN
        DO nsrf = 1, nbsrf
         DO i = 1, klon          
          zt2m(i)  = zt2m(i)  + t2m(i,nsrf)  * pctsrf(i,nsrf)
          zq2m(i)  = zq2m(i)  + q2m(i,nsrf)  * pctsrf(i,nsrf)
!
          DO k = 1, 6
           zn2mout(i,k)  = zn2mout(i,k)  + n2mout(i,nsrf,k)  * pctsrf(i,nsrf)
          ENDDO
!
          zustar(i) = zustar(i) + ustar(i,nsrf) * pctsrf(i,nsrf)
          wstar(i,is_ave)=wstar(i,is_ave)+wstar(i,nsrf)*pctsrf(i,nsrf)
          zu10m(i) = zu10m(i) + u10m(i,nsrf) * pctsrf(i,nsrf)
          zv10m(i) = zv10m(i) + v10m(i,nsrf) * pctsrf(i,nsrf)

          s_pblh(i)   = s_pblh(i)   + pblh(i,nsrf)  * pctsrf(i,nsrf)
          s_plcl(i)   = s_plcl(i)   + plcl(i,nsrf)  * pctsrf(i,nsrf)
          s_capCL(i)  = s_capCL(i)  + capCL(i,nsrf) * pctsrf(i,nsrf)
          s_oliqCL(i) = s_oliqCL(i) + oliqCL(i,nsrf)* pctsrf(i,nsrf)
          s_cteiCL(i) = s_cteiCL(i) + cteiCL(i,nsrf)* pctsrf(i,nsrf)
          s_pblT(i)   = s_pblT(i)   + pblT(i,nsrf)  * pctsrf(i,nsrf)
          s_therm(i)  = s_therm(i)  + therm(i,nsrf) * pctsrf(i,nsrf)
          s_trmb1(i)  = s_trmb1(i)  + trmb1(i,nsrf) * pctsrf(i,nsrf)
          s_trmb2(i)  = s_trmb2(i)  + trmb2(i,nsrf) * pctsrf(i,nsrf)
          s_trmb3(i)  = s_trmb3(i)  + trmb3(i,nsrf) * pctsrf(i,nsrf)
         ENDDO
        ENDDO
       ELSE  !(iflag_split .eq.0)
        DO nsrf = 1, nbsrf
         DO i = 1, klon          
!!! nrlmd le 02/05/2011
          zxfluxlat_x(i) = zxfluxlat_x(i) + fluxlat_x(i,nsrf) * pctsrf(i,nsrf)
          zxfluxlat_w(i) = zxfluxlat_w(i) + fluxlat_w(i,nsrf) * pctsrf(i,nsrf)
!!!
!!! jyg le 08/02/2012
!!  Pour le moment, on sort les valeurs dans (x) et (w) de pblh et de plcl ; 
!!  pour zt2m, on fait la moyenne surfacique sur les sous-surfaces ;
!!  pour qsat2m, on fait la moyenne surfacique sur (x) et (w) ;
!!  pour les autres variables, on sort les valeurs de la region (x).
          zt2m(i)  = zt2m(i)  + (t2m_x(i,nsrf)+wake_s(i)*(t2m_w(i,nsrf)-t2m_x(i,nsrf))) * pctsrf(i,nsrf)
          zq2m(i)  = zq2m(i)  + q2m_x(i,nsrf)  * pctsrf(i,nsrf)
!
          DO k = 1, 6
           zn2mout(i,k)  = zn2mout(i,k)  + n2mout_x(i,nsrf,k)  * pctsrf(i,nsrf)
          ENDDO
!
          zustar(i) = zustar(i) + ustar_x(i,nsrf) * pctsrf(i,nsrf)
          wstar(i,is_ave)=wstar(i,is_ave)+wstar_x(i,nsrf)*pctsrf(i,nsrf)
          zu10m(i) = zu10m(i) + u10m_x(i,nsrf) * pctsrf(i,nsrf)
          zv10m(i) = zv10m(i) + v10m_x(i,nsrf) * pctsrf(i,nsrf)
!
          s_pblh(i)     = s_pblh(i)     + pblh_x(i,nsrf)  * pctsrf(i,nsrf)
          s_pblh_x(i)   = s_pblh_x(i)   + pblh_x(i,nsrf)  * pctsrf(i,nsrf)
          s_pblh_w(i)   = s_pblh_w(i)   + pblh_w(i,nsrf)  * pctsrf(i,nsrf)
!
          s_plcl(i)     = s_plcl(i)     + plcl_x(i,nsrf)  * pctsrf(i,nsrf)
          s_plcl_x(i)   = s_plcl_x(i)   + plcl_x(i,nsrf)  * pctsrf(i,nsrf)
          s_plcl_w(i)   = s_plcl_w(i)   + plcl_w(i,nsrf)  * pctsrf(i,nsrf)
!
          s_capCL(i)  = s_capCL(i)  + capCL_x(i,nsrf) * pctsrf(i,nsrf)
          s_oliqCL(i) = s_oliqCL(i) + oliqCL_x(i,nsrf)* pctsrf(i,nsrf)
          s_cteiCL(i) = s_cteiCL(i) + cteiCL_x(i,nsrf)* pctsrf(i,nsrf)
          s_pblT(i)   = s_pblT(i)   + pblT_x(i,nsrf)  * pctsrf(i,nsrf)
          s_therm(i)  = s_therm(i)  + therm_x(i,nsrf) * pctsrf(i,nsrf)
          s_trmb1(i)  = s_trmb1(i)  + trmb1_x(i,nsrf) * pctsrf(i,nsrf)
          s_trmb2(i)  = s_trmb2(i)  + trmb2_x(i,nsrf) * pctsrf(i,nsrf)
          s_trmb3(i)  = s_trmb3(i)  + trmb3_x(i,nsrf) * pctsrf(i,nsrf)
         ENDDO
        ENDDO
        DO i = 1, klon          
          qsat2m(i)= qsat2m_x(i)+ wake_s(i)*(qsat2m_x(i)-qsat2m_w(i))
        ENDDO
!!!
       ENDIF  ! (iflag_split .eq.0)
!!!

    IF (check) THEN
       amn=MIN(ts(1,is_ter),1000.)
       amx=MAX(ts(1,is_ter),-1000.)
       DO i=2, klon
          amn=MIN(ts(i,is_ter),amn)
          amx=MAX(ts(i,is_ter),amx)
       ENDDO
       PRINT*,' debut apres d_ts min max ftsol(ts)',itap,amn,amx
    ENDIF

!jg ?
!!$!
!!$! If a sub-surface does not exsist for a grid point, the mean value for all 
!!$! sub-surfaces is distributed.
!!$!
!!$    DO nsrf = 1, nbsrf
!!$       DO i = 1, klon
!!$          IF ((pctsrf_new(i,nsrf) .LT. epsfra) .OR. (t2m(i,nsrf).EQ.0.)) THEN
!!$             ts(i,nsrf)     = zxtsol(i)
!!$             t2m(i,nsrf)    = zt2m(i)
!!$             q2m(i,nsrf)    = zq2m(i)
!!$             u10m(i,nsrf)   = zu10m(i)
!!$             v10m(i,nsrf)   = zv10m(i)
!!$
!!$! Les variables qui suivent sont plus utilise, donc peut-etre pas la peine a les mettre ajour
!!$             pblh(i,nsrf)   = s_pblh(i)
!!$             plcl(i,nsrf)   = s_plcl(i)
!!$             capCL(i,nsrf)  = s_capCL(i)
!!$             oliqCL(i,nsrf) = s_oliqCL(i) 
!!$             cteiCL(i,nsrf) = s_cteiCL(i)
!!$             pblT(i,nsrf)   = s_pblT(i)
!!$             therm(i,nsrf)  = s_therm(i)
!!$             trmb1(i,nsrf)  = s_trmb1(i)
!!$             trmb2(i,nsrf)  = s_trmb2(i)
!!$             trmb3(i,nsrf)  = s_trmb3(i)
!!$          ENDIF
!!$       ENDDO
!!$    ENDDO


    DO i = 1, klon
       fder(i) = - 4.0*RSIGMA*zxtsol(i)**3 
    ENDDO
    
    zxqsurf(:) = 0.0
    zxsnow(:)  = 0.0
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          zxqsurf(i) = zxqsurf(i) + MAX(qsurf(i,nsrf),0.0) * pctsrf(i,nsrf)
          zxsnow(i)  = zxsnow(i)  + snow(i,nsrf)  * pctsrf(i,nsrf)
       ENDDO
    ENDDO

! Premier niveau de vent sortie dans physiq.F
    zu1(:) = u(:,1)
    zv1(:) = v(:,1)

  END SUBROUTINE pbl_surface
!
!****************************************************************************************
!
  SUBROUTINE pbl_surface_final(fder_rst, snow_rst, qsurf_rst, ftsoil_rst)

    USE indice_sol_mod

    INCLUDE "dimsoil.h"

! Ouput variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)                 :: fder_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)          :: snow_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)          :: qsurf_rst
    REAL, DIMENSION(klon, nsoilmx, nbsrf), INTENT(OUT) :: ftsoil_rst

 
!****************************************************************************************
! Return module variables for writing to restart file
!
!****************************************************************************************    
    fder_rst(:)       = fder(:)
    snow_rst(:,:)     = snow(:,:)
    qsurf_rst(:,:)    = qsurf(:,:)
    ftsoil_rst(:,:,:) = ftsoil(:,:,:)

!****************************************************************************************
! Deallocate module variables
!
!****************************************************************************************
!   DEALLOCATE(qsol, fder, snow, qsurf, evap, rugos, agesno, ftsoil)
    IF (ALLOCATED(fder)) DEALLOCATE(fder)
    IF (ALLOCATED(snow)) DEALLOCATE(snow)
    IF (ALLOCATED(qsurf)) DEALLOCATE(qsurf)
    IF (ALLOCATED(ftsoil)) DEALLOCATE(ftsoil)

!jyg<
!****************************************************************************************
! Deallocate variables for pbl splitting
!
!****************************************************************************************

    CALL wx_pbl_final
!>jyg

  END SUBROUTINE pbl_surface_final
!  
!****************************************************************************************
! 

!albedo SB >>>
SUBROUTINE pbl_surface_newfrac(itime, pctsrf_new, pctsrf_old, &
     evap, z0m, z0h, agesno,                                  &
     tsurf,alb_dir,alb_dif, ustar, u10m, v10m, tke)  
!albedo SB <<<
    ! Give default values where new fraction has appread

    USE indice_sol_mod

    INCLUDE "dimsoil.h"
    INCLUDE "clesphys.h"
    INCLUDE "compbl.h"

! Input variables
!****************************************************************************************
    INTEGER, INTENT(IN)                     :: itime
    REAL, DIMENSION(klon,nbsrf), INTENT(IN) :: pctsrf_new, pctsrf_old

! InOutput variables
!****************************************************************************************
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT)        :: tsurf
!albedo SB >>>
    REAL, DIMENSION(klon,nsw,nbsrf), INTENT(INOUT)       :: alb_dir, alb_dif 
    INTEGER :: k
!albedo SB <<<
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT)        :: ustar,u10m, v10m
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT)        :: evap, agesno
    REAL, DIMENSION(klon,nbsrf+1), INTENT(INOUT)        :: z0m,z0h
    REAL, DIMENSION(klon,klev+1,nbsrf+1), INTENT(INOUT) :: tke

! Local variables
!****************************************************************************************
    INTEGER           :: nsrf, nsrf_comp1, nsrf_comp2, nsrf_comp3, i
    CHARACTER(len=80) :: abort_message
    CHARACTER(len=20) :: modname = 'pbl_surface_newfrac'
    INTEGER, DIMENSION(nbsrf) :: nfois=0, mfois=0, pfois=0
!
! All at once !! 
!****************************************************************************************
    
    DO nsrf = 1, nbsrf
       ! First decide complement sub-surfaces
       SELECT CASE (nsrf)
       CASE(is_oce)
          nsrf_comp1=is_sic
          nsrf_comp2=is_ter
          nsrf_comp3=is_lic
       CASE(is_sic)
          nsrf_comp1=is_oce
          nsrf_comp2=is_ter
          nsrf_comp3=is_lic
       CASE(is_ter)
          nsrf_comp1=is_lic
          nsrf_comp2=is_oce
          nsrf_comp3=is_sic
       CASE(is_lic)
          nsrf_comp1=is_ter
          nsrf_comp2=is_oce
          nsrf_comp3=is_sic
       END SELECT

       ! Initialize all new fractions
       DO i=1, klon
          IF (pctsrf_new(i,nsrf) > 0. .AND. pctsrf_old(i,nsrf) == 0.) THEN
             
             IF (pctsrf_old(i,nsrf_comp1) > 0.) THEN
                ! Use the complement sub-surface, keeping the continents unchanged
                qsurf(i,nsrf) = qsurf(i,nsrf_comp1)
                evap(i,nsrf)  = evap(i,nsrf_comp1)
                z0m(i,nsrf) = z0m(i,nsrf_comp1)
                z0h(i,nsrf) = z0h(i,nsrf_comp1)
                tsurf(i,nsrf) = tsurf(i,nsrf_comp1)
!albedo SB >>>
                DO k=1,nsw
                 alb_dir(i,k,nsrf)=alb_dir(i,k,nsrf_comp1)
                 alb_dif(i,k,nsrf)=alb_dif(i,k,nsrf_comp1)
                ENDDO
!albedo SB <<<
                ustar(i,nsrf)  = ustar(i,nsrf_comp1)
                u10m(i,nsrf)  = u10m(i,nsrf_comp1)
                v10m(i,nsrf)  = v10m(i,nsrf_comp1)
                IF (iflag_pbl > 1) THEN
                 tke(i,:,nsrf) = tke(i,:,nsrf_comp1)
                ENDIF
                mfois(nsrf) = mfois(nsrf) + 1
                ! F. Codron sensible default values for ocean and sea ice
                IF (nsrf.EQ.is_oce) THEN
                   tsurf(i,nsrf) = 271.35
                   ! (temperature of sea water under sea ice, so that
                   ! is also the temperature of appearing sea water)
                   DO k=1,nsw
                      alb_dir(i,k,nsrf) = 0.06 ! typical Ocean albedo
                      alb_dif(i,k,nsrf) = 0.06
                   ENDDO
                ELSE IF (nsrf.EQ.is_sic) THEN
                   tsurf(i,nsrf) = 271.35
                   ! (Temperature at base of sea ice. Surface
                   ! temperature could be higher, up to 0 Celsius
                   ! degrees. We set it to -1.8 Celsius degrees for
                   ! consistency with the ocean slab model.)
                   DO k=1,nsw
                      alb_dir(i,k,nsrf) = 0.3 ! thin ice
                      alb_dif(i,k,nsrf) = 0.3
                   ENDDO
                ENDIF
             ELSE
                ! The continents have changed. The new fraction receives the mean sum of the existent fractions
                qsurf(i,nsrf) = qsurf(i,nsrf_comp2)*pctsrf_old(i,nsrf_comp2) + qsurf(i,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
                evap(i,nsrf)  = evap(i,nsrf_comp2) *pctsrf_old(i,nsrf_comp2) + evap(i,nsrf_comp3) *pctsrf_old(i,nsrf_comp3)
                z0m(i,nsrf) = z0m(i,nsrf_comp2)*pctsrf_old(i,nsrf_comp2) + z0m(i,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
                z0h(i,nsrf) = z0h(i,nsrf_comp2)*pctsrf_old(i,nsrf_comp2) + z0h(i,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
                tsurf(i,nsrf) = tsurf(i,nsrf_comp2)*pctsrf_old(i,nsrf_comp2) + tsurf(i,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
!albedo SB >>>
                DO k=1,nsw
                 alb_dir(i,k,nsrf)=alb_dir(i,k,nsrf_comp2)*pctsrf_old(i,nsrf_comp2)+&
                                        alb_dir(i,k,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
                 alb_dif(i,k,nsrf)=alb_dif(i,k,nsrf_comp2)*pctsrf_old(i,nsrf_comp2)+&
                                        alb_dif(i,k,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
                ENDDO
!albedo SB <<<
                ustar(i,nsrf)  = ustar(i,nsrf_comp2) *pctsrf_old(i,nsrf_comp2) + ustar(i,nsrf_comp3) *pctsrf_old(i,nsrf_comp3)
                u10m(i,nsrf)  = u10m(i,nsrf_comp2) *pctsrf_old(i,nsrf_comp2) + u10m(i,nsrf_comp3) *pctsrf_old(i,nsrf_comp3)
                v10m(i,nsrf)  = v10m(i,nsrf_comp2) *pctsrf_old(i,nsrf_comp2) + v10m(i,nsrf_comp3) *pctsrf_old(i,nsrf_comp3)
                IF (iflag_pbl > 1) THEN
                 tke(i,:,nsrf) = tke(i,:,nsrf_comp2)*pctsrf_old(i,nsrf_comp2) + tke(i,:,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
                ENDIF
            
                ! Security abort. This option has never been tested. To test, comment the following line.
!                abort_message='The fraction of the continents have changed!'
!                CALL abort_physic(modname,abort_message,1)
                nfois(nsrf) = nfois(nsrf) + 1
             ENDIF
             snow(i,nsrf)     = 0.
             agesno(i,nsrf)   = 0.
             ftsoil(i,:,nsrf) = tsurf(i,nsrf)
          ELSE
             pfois(nsrf) = pfois(nsrf)+ 1
          ENDIF
       ENDDO
       
    ENDDO

  END SUBROUTINE pbl_surface_newfrac
!  
!****************************************************************************************
!  
END MODULE pbl_surface_mod
