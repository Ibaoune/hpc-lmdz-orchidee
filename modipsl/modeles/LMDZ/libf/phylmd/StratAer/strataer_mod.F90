! $Id: strataer_mod.F90 3527 2019-05-30 13:43:48Z oboucher $
MODULE strataer_mod
! This module contains information about strato microphysic model parameters
  
  IMPLICIT NONE

  ! flag to constraint nucleation rate in a lat/pres box
  LOGICAL,SAVE :: flag_nuc_rate_box      ! Nucleation rate limit or not to a lat/pres
  !$OMP THREADPRIVATE(flag_nuc_rate_box)
  REAL,SAVE    :: nuclat_min             ! min lat to activate nuc rate
  REAL,SAVE    :: nuclat_max             ! max lat to activate nuc rate
  REAL,SAVE    :: nucpres_min            ! min pres to activate nuc rate
  REAL,SAVE    :: nucpres_max            ! max pres to activate nuc rate
  !$OMP THREADPRIVATE(nuclat_min, nuclat_max, nucpres_min, nucpres_max)
  
  ! flag for sulfur emission scenario: (0) background aer ; (1) volcanic eruption ; (2) strato aer injections (SAI)
  INTEGER,SAVE :: flag_sulf_emit
  !$OMP THREADPRIVATE(flag_sulf_emit)

  ! flag for sulfur emission altitude distribution: (0) gaussian; (1) uniform
  INTEGER,SAVE :: flag_sulf_emit_distrib
  !$OMP THREADPRIVATE(flag_sulf_emit_distrib)
  
  !--flag_sulf_emit=1 -- Volcanic eruption(s)
  INTEGER,SAVE             :: nErupt                    ! number of eruptions specs
  REAL,SAVE                :: injdur                    ! volcanic injection duration
  !$OMP THREADPRIVATE(nErupt, injdur)
  INTEGER,ALLOCATABLE,SAVE :: year_emit_vol(:)          ! year of emission date
  INTEGER,ALLOCATABLE,SAVE :: mth_emit_vol(:)           ! month of emission date
  INTEGER,ALLOCATABLE,SAVE :: day_emit_vol(:)           ! day of emission date
  !$OMP THREADPRIVATE(year_emit_vol, mth_emit_vol, day_emit_vol)
  REAL,ALLOCATABLE,SAVE    :: m_aer_emiss_vol(:)        ! emitted sulfur mass in kgS, e.g. 7Tg(S)=14Tg(SO2)
  REAL,ALLOCATABLE,SAVE    :: altemiss_vol(:)           ! emission altitude in m
  REAL,ALLOCATABLE,SAVE    :: sigma_alt_vol(:)          ! standard deviation of emission altitude in m
  !$OMP THREADPRIVATE(m_aer_emiss_vol, altemiss_vol, sigma_alt_vol)
  INTEGER,ALLOCATABLE,SAVE :: ponde_lonlat_vol(:)       ! lon/lat ponderation factor
  REAL,ALLOCATABLE,SAVE    :: xlat_min_vol(:)           ! min latitude of volcano in degree
  REAL,ALLOCATABLE,SAVE    :: xlat_max_vol(:)           ! max latitude of volcano in degree
  REAL,ALLOCATABLE,SAVE    :: xlon_min_vol(:)           ! min longitude of volcano in degree
  REAL,ALLOCATABLE,SAVE    :: xlon_max_vol(:)           ! max longitude of volcano in degree
  !$OMP THREADPRIVATE(ponde_lonlat_vol, xlat_min_vol, xlat_max_vol, xlon_min_vol, xlon_max_vol)
  
  !--flag_sulf_emit=2 --SAI
  REAL,SAVE    :: m_aer_emiss_sai        ! emitted sulfur mass in kgS, eg 1e9=1TgS, 1e10=10TgS
  REAL,SAVE    :: altemiss_sai           ! emission altitude in m
  REAL,SAVE    :: sigma_alt_sai          ! standard deviation of emission altitude in m
  !$OMP THREADPRIVATE(m_aer_emiss_sai, altemiss_sai, sigma_alt_sai)
  REAL,SAVE    :: xlat_sai               ! latitude of SAI in degree
  REAL,SAVE    :: xlon_sai               ! longitude of SAI in degree
  REAL,SAVE    :: dlat, dlon             ! delta latitude and d longitude of grid in degree
  !$OMP THREADPRIVATE(xlat_sai, xlon_sai, dlat, dlon)
  
  !--flag_sulf_emit=3 -- SAI
  REAL,SAVE    :: xlat_max_sai           ! maximum latitude of SAI in degrees
  REAL,SAVE    :: xlat_min_sai           ! minimum latitude of SAI in degrees
  !$OMP THREADPRIVATE(xlat_min_sai,xlat_max_sai)
  
CONTAINS
    
  SUBROUTINE strataer_init()
    USE ioipsl_getin_p_mod, ONLY : getin_p
    USE print_control_mod, ONLY : lunout
    USE mod_phys_lmdz_para, ONLY : is_master

    ! Local var
    INTEGER       :: ieru
 
    INTEGER :: i
    
    WRITE(lunout,*) 'IN STRATAER INIT WELCOME!'
   
    !Config Key  = flag_sulf_emit
    !Config Desc = aerosol emission mode
    ! - 0 = background aerosol
    ! - 1 = volcanic eruption
    ! - 2 = geo-ingeneering design
    ! - 3 = geo-engineering between two latitudes
    !Config Def  = 0
    !Config Help = Used in physiq.F
    !
    flag_sulf_emit = 0
    nErupt = 0 ! eruption number
    injdur = 0 ! init injection duration
    CALL getin_p('flag_sulf_emit',flag_sulf_emit)

    IF (flag_sulf_emit==1) THEN ! Volcano
       CALL getin_p('nErupt',nErupt)
       CALL getin_p('injdur',injdur)
    ELSEIF (flag_sulf_emit == 2) THEN ! SAI
       CALL getin_p('m_aer_emiss_sai',m_aer_emiss_sai)
       CALL getin_p('altemiss_sai',altemiss_sai)
       CALL getin_p('sigma_alt_sai',sigma_alt_sai)
       CALL getin_p('xlat_sai',xlat_sai)
       CALL getin_p('xlon_sai',xlon_sai)
       CALL getin_p('flag_sulf_emit_distrib',flag_sulf_emit_distrib)
    ELSEIF (flag_sulf_emit == 3) THEN ! SAI between latitudes
       CALL getin_p('m_aer_emiss_sai',m_aer_emiss_sai)
       CALL getin_p('altemiss_sai',altemiss_sai)
       CALL getin_p('sigma_alt_sai',sigma_alt_sai)
       CALL getin_p('xlon_sai',xlon_sai)
       CALL getin_p('xlat_max_sai',xlat_max_sai)
       CALL getin_p('xlat_min_sai',xlat_min_sai)
       CALL getin_p('flag_sulf_emit_distrib',flag_sulf_emit_distrib)
    ENDIF

    ALLOCATE(year_emit_vol(nErupt),mth_emit_vol(nErupt),day_emit_vol(nErupt))
    ALLOCATE(m_aer_emiss_vol(nErupt),altemiss_vol(nErupt),sigma_alt_vol(nErupt))
    ALLOCATE(xlat_min_vol(nErupt),xlon_min_vol(nErupt))
    ALLOCATE(xlat_max_vol(nErupt),xlon_max_vol(nErupt))
    
    IF (nErupt.GT.0) THEN 
      year_emit_vol=0 ; mth_emit_vol=0 ; day_emit_vol=0
      m_aer_emiss_vol=0. ; altemiss_vol=0. ; sigma_alt_vol=0.
      xlon_min_vol=0. ; xlon_max_vol=0.
      xlat_min_vol=0. ; xlat_max_vol=0.
    ENDIF
    
    CALL getin_p('year_emit_vol',year_emit_vol)
    CALL getin_p('mth_emit_vol',mth_emit_vol)
    CALL getin_p('day_emit_vol',day_emit_vol)
    CALL getin_p('m_aer_emiss_vol',m_aer_emiss_vol)
    CALL getin_p('altemiss_vol',altemiss_vol)
    CALL getin_p('sigma_alt_vol',sigma_alt_vol)
    CALL getin_p('xlon_min_vol',xlon_min_vol)
    CALL getin_p('xlon_max_vol',xlon_max_vol)
    CALL getin_p('xlat_min_vol',xlat_min_vol)
    CALL getin_p('xlat_max_vol',xlat_max_vol)
    !Config Key  = flag_nuc_rate_box
    !Config Desc = define or not a box for nucleation rate
    ! - F = global nucleation
    ! - T = 2D-box for nucleation need nuclat_min, nuclat_max, nucpres_min and
    ! nucpres_max
    !       to define its bounds.
    !Config Def  = F
    !Config Help = Used in physiq.F
    !
    flag_nuc_rate_box = .FALSE.
    CALL getin_p('flag_nuc_rate_box',flag_nuc_rate_box)
    CALL getin_p('nuclat_min',nuclat_min)
    CALL getin_p('nuclat_max',nuclat_max)
    CALL getin_p('nucpres_min',nucpres_min)
    CALL getin_p('nucpres_max',nucpres_max)

    !IF (is_master) THEN
    WRITE(lunout,*) 'flag_sulf_emit = ',flag_sulf_emit
    IF (flag_sulf_emit == 1) THEN
       WRITE(lunout,*) 'IN STRATAER nErupt: ',nErupt
       WRITE(lunout,*) 'IN STRATAER injdur: ',injdur
       WRITE(lunout,*) 'IN STRATAER : year_emit_vol',year_emit_vol
       WRITE(lunout,*) 'IN STRATAER : mth_emit_vol',mth_emit_vol
       WRITE(lunout,*) 'IN STRATAER : day_emit_vol',day_emit_vol
       WRITE(lunout,*) 'IN STRATAER : m_aer_emiss_vol',m_aer_emiss_vol
       WRITE(lunout,*) 'IN STRATAER : altemiss_vol',altemiss_vol
       WRITE(lunout,*) 'IN STRATAER : sigma_alt_vol',sigma_alt_vol
       WRITE(lunout,*) 'IN STRATAER : ponde_lonlat_vol',ponde_lonlat_vol
       WRITE(lunout,*) 'IN STRATAER : xlat_min_vol',xlat_min_vol
       WRITE(lunout,*) 'IN STRATAER : xlat_max_vol',xlat_max_vol
       WRITE(lunout,*) 'IN STRATAER : xlon_min_vol',xlon_min_vol
       WRITE(lunout,*) 'IN STRATAER : xlon_max_vol',xlon_max_vol
    ELSEIF (flag_sulf_emit == 2) THEN
       WRITE(lunout,*) 'IN STRATAER : m_aer_emiss_sai',m_aer_emiss_sai
       WRITE(lunout,*) 'IN STRATAER : altemiss_sai',altemiss_sai
       WRITE(lunout,*) 'IN STRATAER : sigma_alt_sai',sigma_alt_sai
       WRITE(lunout,*) 'IN STRATAER : xlat_sai',xlat_sai
       WRITE(lunout,*) 'IN STRATAER : xlon_sai',xlon_sai
       WRITE(lunout,*) 'flag_sulf_emit_distrib = ',flag_sulf_emit_distrib
    ELSEIF (flag_sulf_emit == 3) THEN
       WRITE(lunout,*) 'IN STRATAER : m_aer_emiss_sai',m_aer_emiss_sai
       WRITE(lunout,*) 'IN STRATAER : altemiss_sai',altemiss_sai
       WRITE(lunout,*) 'IN STRATAER : sigma_alt_sai',sigma_alt_sai
       WRITE(lunout,*) 'IN STRATAER : xlat_min_sai',xlat_min_sai
       WRITE(lunout,*) 'IN STRATAER : xlat_max_sai',xlat_max_sai
       WRITE(lunout,*) 'IN STRATAER : xlon_sai',xlon_sai
       WRITE(lunout,*) 'flag_sulf_emit_distrib = ',flag_sulf_emit_distrib
    ENDIF
    WRITE(lunout,*) 'IN STRATAER : flag_nuc_rate_box = ',flag_nuc_rate_box
    IF (flag_nuc_rate_box) THEN
       WRITE(lunout,*) 'IN STRATAER : nuclat_min = ',nuclat_min,', nuclat_max = ',nuclat_max
       WRITE(lunout,*) 'IN STRATAER : nucpres_min = ',nucpres_min,', nucpres_max = ',nucpres_max
    ENDIF
    !ENDIF

    CALL strataer_ponde_init
    WRITE(lunout,*) 'IN STRATAER INT2 END'

  END SUBROUTINE strataer_init
  
  ! Compute the ponderation to applicate in each grid point for all eruptions and init
  ! dlat & dlon variables
  SUBROUTINE strataer_ponde_init()
    
    USE regular_lonlat_mod, ONLY: lon_reg, lat_reg
    USE dimphy, ONLY: klon
    USE mod_grid_phy_lmdz, ONLY: nbp_lat, nbp_lon
    USE print_control_mod, ONLY : lunout
    USE YOMCST, ONLY : RPI

    ! local var
    REAL                :: pi,lat_reg_deg,lon_reg_deg! latitude and d longitude of grid in degree
    INTEGER             :: ieru, i, j
    
    ALLOCATE(ponde_lonlat_vol(nErupt))
    
    !Compute lon/lat ponderation for injection
    dlat=180./2./FLOAT(nbp_lat)   ! d latitude in degree
    dlon=360./2./FLOAT(nbp_lon)   ! d longitude in degree
    WRITE(lunout,*) 'IN STRATAER_INIT dlat=',dlat,'dlon=',dlon
    WRITE(lunout,*) 'IN STRATAER_INIT nErupt=',nErupt
    WRITE(lunout,*) 'IN STRATAER_INIT xlat_min=',xlat_min_vol,'xlat_max=',xlat_max_vol
    WRITE(lunout,*) 'IN STRATAER_INIT xlon_min=',xlon_min_vol,'xlon_max=',xlon_max_vol
    DO ieru=1, nErupt
       ponde_lonlat_vol(ieru) = 0
       DO i=1,nbp_lon
          lon_reg_deg = lon_reg(i)*180./RPI
          DO j=1,nbp_lat
             lat_reg_deg = lat_reg(j)*180./RPI
             IF  ( lat_reg_deg.GE.xlat_min_vol(ieru)-dlat .AND. lat_reg_deg.LT.xlat_max_vol(ieru)+dlat .AND. &
                  lon_reg_deg.GE.xlon_min_vol(ieru)-dlon .AND. lon_reg_deg.LT.xlon_max_vol(ieru)+dlon ) THEN
                ponde_lonlat_vol(ieru) = ponde_lonlat_vol(ieru) + 1
             ENDIF
          ENDDO
       ENDDO
       IF(ponde_lonlat_vol(ieru) == 0) THEN
          WRITE(lunout,*) 'STRATAER_INIT ERROR: no grid point found for eruption ieru=',ieru
       ENDIF
    ENDDO !ieru
    WRITE(lunout,*) 'IN STRATAER_INIT ponde_lonlat: ',ponde_lonlat_vol
    
  END SUBROUTINE strataer_ponde_init
  
END MODULE strataer_mod
