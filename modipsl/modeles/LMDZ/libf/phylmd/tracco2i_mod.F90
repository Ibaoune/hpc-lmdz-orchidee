MODULE tracco2i_mod
!
! This module does the work for the interactive CO2 tracers
! Authors: Patricia Cadule and Olivier Boucher
!
! Purpose and description:
!  -----------------------
! Main routine for the interactive carbon cycle
! Gather all carbon fluxes and emissions from ORCHIDEE, PISCES and fossil fuel
! Compute the net flux in source field which is used in phytrac
! Compute global CO2 mixing ratio for radiation scheme if option is activated
! Redistribute CO2 evenly over the atmosphere if transport is desactivated
!
CONTAINS

  SUBROUTINE tracco2i_init()
    ! This subroutine calls carbon_cycle_init needed to be done before first call to phys_output_write in physiq.
    USE carbon_cycle_mod, ONLY: carbon_cycle_init, carbon_cycle_cpl

    ! Initialize carbon_cycle_mod
    IF (carbon_cycle_cpl) THEN
       CALL carbon_cycle_init()
    ENDIF

  END SUBROUTINE tracco2i_init

  SUBROUTINE tracco2i(pdtphys, debutphy, &
       xlat, xlon, pphis, pphi, &
       t_seri, pplay, paprs, tr_seri, source)

    USE dimphy
    USE infotrac_phy
    USE geometry_mod, ONLY: cell_area
    USE carbon_cycle_mod, ONLY: id_CO2, nbcf_in, fields_in, cfname_in
    USE carbon_cycle_mod, ONLY: fco2_ocn_day, fco2_ff, fco2_bb, fco2_land, fco2_ocean
    USE carbon_cycle_mod, ONLY: fco2_land_nbp, fco2_land_nep, fco2_land_fLuc
    USE carbon_cycle_mod, ONLY: fco2_land_fwoodharvest, fco2_land_fHarvest
    USE carbon_cycle_mod, ONLY: carbon_cycle_cpl, carbon_cycle_tr, carbon_cycle_rad, RCO2_glo, RCO2_tot
    USE mod_grid_phy_lmdz
    USE mod_phys_lmdz_mpi_data, ONLY: is_mpi_root
    USE mod_phys_lmdz_para, ONLY: gather, bcast, scatter
    USE phys_cal_mod
    USE phys_state_var_mod, ONLY: pctsrf
    USE indice_sol_mod, ONLY: nbsrf, is_ter, is_lic, is_oce, is_sic

    IMPLICIT NONE

    INCLUDE "clesphys.h"
    INCLUDE "YOMCST.h"

! Input argument
!---------------
    REAL,INTENT(IN)    :: pdtphys    ! Pas d'integration pour la physique (seconde)
    LOGICAL,INTENT(IN) :: debutphy   ! le flag de l'initialisation de la physique

    REAL,DIMENSION(klon),INTENT(IN)        :: xlat    ! latitudes pour chaque point 
    REAL,DIMENSION(klon),INTENT(IN)        :: xlon    ! longitudes pour chaque point 
    REAL,DIMENSION(klon),INTENT(IN)        :: pphis   ! geopotentiel du sol
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pphi    ! geopotentiel de chaque couche

    REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
    REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs   ! pression pour chaque inter-couche (en Pa)
    REAL,DIMENSION(klon,nbtr),INTENT(INOUT):: source  ! flux de traceur [U/m2/s]

! Output argument
!----------------
    REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT)  :: tr_seri ! Concentration Traceur [U/kgA]  

! Local variables
!----------------

    INTEGER                        :: it, k, i, nb
    REAL, DIMENSION(klon,klev)     :: m_air     ! mass of air in every grid box [kg]
    REAL, DIMENSION(klon_glo,klev) :: co2_glo   ! variable temporaire sur la grille global
    REAL, DIMENSION(klon_glo,klev) :: m_air_glo ! variable temporaire sur la grille global

    LOGICAL, SAVE :: check_fCO2_nbp_in_cfname
!$OMP THREADPRIVATE(check_fCO2_nbp_in_cfname)
    INTEGER, SAVE :: day_pre=-1
!$OMP THREADPRIVATE(day_pre)

    IF (is_mpi_root) THEN
      PRINT *,'in tracco2i: date from phys_cal_mod =',year_cur,'-',mth_cur,'-',day_cur,'-',hour
    ENDIF

!--initialisation of CO2 field if not in restart file
!--dirty way of doing, do it better later
!--convert 280 ppm into kg CO2 / kg air
    IF (debutphy) THEN

! Initialisation de tr_seri(id_CO2) si pas initialise
      IF (MAXVAL(tr_seri(:,:,id_CO2)).LT.1.e-15) THEN 
        tr_seri(:,:,id_CO2)=co2_ppm0*1.e-6/RMD*RMCO2 !--initialised from co2_ppm0 in rdem
      ENDIF 

!--check if fCO2_nbp is in 
      check_fCO2_nbp_in_cfname=.FALSE.
      DO nb=1, nbcf_in
        IF (cfname_in(nb)=="fCO2_nbp") check_fCO2_nbp_in_cfname=.TRUE.
      ENDDO

    ENDIF

!--calculate mass of air in every grid box in kg air
    DO i=1, klon
    DO k=1, klev
      m_air(i,k)=(paprs(i,k)-paprs(i,k+1))/RG*cell_area(i)
    ENDDO
    ENDDO

!--call CO2 emission routine 
!--co2bb is zero for now
!--unit kg CO2 m-2 s-1 
    CALL co2_emissions(debutphy)

!--retrieving land and ocean CO2 flux
    fco2_land(:)=0.0
    fco2_ocean(:)=0.0
    fco2_land_nbp(:)=0.
    fco2_land_nep(:)=0.
    fco2_land_fLuc(:)=0.
    fco2_land_fwoodharvest(:)=0.
    fco2_land_fHarvest(:)=0.

    DO nb=1, nbcf_in

      SELECT CASE(cfname_in(nb))
!--dealing with the different fluxes coming from ORCHIDEE
!--fluxes come in unit of kg C m-2 s-1 is converted into kg CO2 m-2 s-1
      CASE("fCO2_nep")
          fco2_land_nep(:)=fields_in(:,nb)*RMCO2/RMC*pctsrf(:,is_ter)
      CASE("fCO2_fLuc")
          fco2_land_fLuc(:)=fields_in(:,nb)*RMCO2/RMC*pctsrf(:,is_ter)
      CASE("fCO2_fwoodharvest")
          fco2_land_fwoodharvest(:)=fields_in(:,nb)*RMCO2/RMC*pctsrf(:,is_ter)
      CASE("fCO2_fHarvest")
          fco2_land_fHarvest(:)=fields_in(:,nb)*RMCO2/RMC*pctsrf(:,is_ter)
      CASE("fCO2_nbp")
          fco2_land_nbp(:)=fields_in(:,nb)*RMCO2/RMC*pctsrf(:,is_ter)
!--fCO2_fco2_ocn comes in unit of mol C02 m-2 s-1 is converted into kg CO2 m-2 s-1 + change sign
      CASE("fCO2_fgco2") 
          fco2_ocean(:)=-1.*fco2_ocn_day(:)*RMCO2/1.e3*(pctsrf(:,is_oce)+pctsrf(:,is_sic))
      END SELECT

    ENDDO

!--if fCO2_nbp is transferred we use it, otherwise we use the sum of what has been passed from ORCHIDEE
    IF (check_fCO2_nbp_in_cfname)  THEN
       fco2_land(:)=fco2_land_nbp(:)
    ELSE
       fco2_land(:)=fco2_land_nep(:)+fco2_land_fLuc(:)+fco2_land_fwoodharvest(:)+fco2_land_fHarvest(:)
    ENDIF

!!--preparing the net anthropogenic flux at the surface for mixing layer
!!--unit kg CO2 / m2 / s
!    PRINT *, 'tracco2i_mod.F90 --- MAXVAL(fco2_ff) ',MAXVAL(fco2_ff)
!    PRINT *, 'tracco2i_mod.F90 --- MINVAL(fco2_ff) ',MINVAL(fco2_ff)
!    PRINT *, 'tracco2i_mod.F90 --- MAXVAL(fco2_bb) ',MAXVAL(fco2_bb)
!    PRINT *, 'tracco2i_mod.F90 --- MINVAL(fco2_bb) ',MINVAL(fco2_bb)
!    PRINT *, 'tracco2i_mod.F90 --- MAXVAL(fco2_land) ',MAXVAL(fco2_land)
!    PRINT *, 'tracco2i_mod.F90 --- MINVAL(fco2_land) ',MINVAL(fco2_land)
!    PRINT *, 'tracco2i_mod.F90 --- MAXVAL(fco2_ocean) ',MAXVAL(fco2_ocean)
!    PRINT *, 'tracco2i_mod.F90 --- MINVAL(fco2_ocean) ',MINVAL(fco2_ocean)
!    PRINT *, 'tracco2i_mod.F90 --- MAXVAL(source(:,id_CO2)) ',MAXVAL(source(:,id_CO2))
!    PRINT *, 'tracco2i_mod.F90 --- MINVAL(source(:,id_CO2)) ',MINVAL(source(:,id_CO2))
!
!--build final source term for CO2
    source(:,id_CO2)=fco2_ff(:)+fco2_bb(:)+fco2_land(:)+fco2_ocean(:)

!--computing global mean CO2 for radiation
!--for every timestep comment out the IF ENDIF statements 
!--otherwise this is updated every day
    IF (debutphy.OR.day_cur.NE.day_pre) THEN

      CALL gather(tr_seri(:,:,id_CO2),co2_glo)
      CALL gather(m_air,m_air_glo)

!$OMP MASTER

!--compute a global mean CO2 value and print its value in ppm
       IF (is_mpi_root) THEN
         RCO2_tot=SUM(co2_glo*m_air_glo)  !--unit kg CO2
         RCO2_glo=RCO2_tot/SUM(m_air_glo) !--unit kg CO2 / kg air
         PRINT *,'tracco2i: global CO2 in ppm =', RCO2_glo*1.e6*RMD/RMCO2
         PRINT *,'tracco2i: total CO2 in kg =', RCO2_tot
       ENDIF
!$OMP END MASTER
       CALL bcast(RCO2_glo)
       day_pre=day_cur
!--if not carbon_cycle_tr, then we reinitialize the CO2 each day to its global mean value
       IF (.NOT.carbon_cycle_tr) THEN 
         tr_seri(:,:,id_CO2)=RCO2_glo
       ENDIF
    ENDIF

  END SUBROUTINE tracco2i

  SUBROUTINE co2_emissions(debutphy)

    USE dimphy
    USE infotrac_phy
    USE geometry_mod, ONLY : cell_area
    USE mod_grid_phy_lmdz
    USE mod_phys_lmdz_mpi_data, ONLY: is_mpi_root
    USE mod_phys_lmdz_para, ONLY: gather, scatter
    USE phys_cal_mod

    USE netcdf95, ONLY: nf95_close, nf95_gw_var, nf95_inq_varid, nf95_open
    USE netcdf, ONLY: nf90_get_var, nf90_noerr, nf90_nowrite

    USE carbon_cycle_mod, ONLY : fco2_ff, fco2_bb, fco2_land, fco2_ocean

    IMPLICIT NONE

    INCLUDE "YOMCST.h"
    LOGICAL,INTENT(IN) :: debutphy

! For NetCDF:
    INTEGER ncid_in  ! IDs for input files
    INTEGER varid, ncerr

    INTEGER :: n_glo, n_month
    REAL, POINTER:: vector(:), time(:)
    REAL,ALLOCATABLE       :: flx_co2ff_glo(:,:) !  fossil-fuel CO2 
    REAL,ALLOCATABLE       :: flx_co2bb_glo(:,:) !  biomass-burning CO2 
    REAL,ALLOCATABLE, SAVE :: flx_co2ff(:,:)     !  fossil-fuel CO2 
    REAL,ALLOCATABLE, SAVE :: flx_co2bb(:,:)     !  biomass-burning CO2 
!$OMP THREADPRIVATE(flx_co2ff,flx_co2bb)

!! may be controlled via the .def later on
!! also co2bb for now comes from ORCHIDEE
    LOGICAL, PARAMETER :: readco2ff=.TRUE.
!! this should be left to FALSE for now
    LOGICAL, PARAMETER :: readco2bb=.FALSE. 

    CHARACTER (len = 20) :: modname = 'tracco2i.co2_emissions'
    CHARACTER (len = 80) :: abort_message

    IF (debutphy) THEN 

    ALLOCATE(flx_co2ff(klon,12))
    ALLOCATE(flx_co2bb(klon,12))

!$OMP MASTER
    IF (is_mpi_root) THEN
    
      IF (.NOT.ALLOCATED(flx_co2ff_glo)) ALLOCATE(flx_co2ff_glo(klon_glo,12))
      IF (.NOT.ALLOCATED(flx_co2bb_glo)) ALLOCATE(flx_co2bb_glo(klon_glo,12))

!--reading CO2 fossil fuel emissions
      IF (readco2ff) THEN 

        ! ... Open the COZff file
        CALL nf95_open("sflx_lmdz_co2_ff.nc", nf90_nowrite, ncid_in)

        CALL nf95_inq_varid(ncid_in, "vector", varid)
        CALL nf95_gw_var(ncid_in, varid, vector)
        n_glo = size(vector)
        IF (n_glo.NE.klon_glo) THEN
           abort_message='sflx_lmdz_co2_ff: le nombre de points n est pas egal a klon_glo'
           CALL abort_physic(modname,abort_message,1)
        ENDIF

        CALL nf95_inq_varid(ncid_in, "time", varid)
        CALL nf95_gw_var(ncid_in, varid, time)
        n_month = size(time)
        IF (n_month.NE.12) THEN
           abort_message='sflx_lmdz_co2_ff: le nombre de month n est pas egal a 12'
           CALL abort_physic(modname,abort_message,1)
        ENDIF

!--reading flx_co2 for fossil fuel
        CALL nf95_inq_varid(ncid_in, "flx_co2", varid)
        ncerr = nf90_get_var(ncid_in, varid, flx_co2ff_glo)

        CALL nf95_close(ncid_in)
    
      ELSE  !--co2ff not to be read
        flx_co2ff_glo(:,:)=0.0
      ENDIF

!--reading CO2 biomass burning emissions
!--using it will be inconsistent with treatment in ORCHIDEE
      IF (readco2bb) THEN 

      ! ... Open the CO2bb file
      CALL nf95_open("sflx_lmdz_co2_bb.nc", nf90_nowrite, ncid_in)

      CALL nf95_inq_varid(ncid_in, "vector", varid)
      CALL nf95_gw_var(ncid_in, varid, vector)
      n_glo = size(vector)
      IF (n_glo.NE.klon_glo) THEN
         abort_message='sflx_lmdz_co2_bb: le nombre de points n est pas egal a klon_glo'
         CALL abort_physic(modname,abort_message,1)
      ENDIF

      CALL nf95_inq_varid(ncid_in, "time", varid)
      CALL nf95_gw_var(ncid_in, varid, time)
      n_month = size(time)
      IF (n_month.NE.12) THEN
         abort_message='sflx_lmdz_co2_bb: le nombre de month n est pas egal a 12'
         CALL abort_physic(modname,abort_message,1)
      ENDIF

!--reading flx_co2 for biomass burning
      CALL nf95_inq_varid(ncid_in, "flx_co2", varid)
      ncerr = nf90_get_var(ncid_in, varid, flx_co2bb_glo)

      CALL nf95_close(ncid_in)
    
      ELSE  !--co2bb not to be read
        flx_co2bb_glo(:,:)=0.0
      ENDIF

    ENDIF
!$OMP END MASTER

    ! Allocation needed for all proc otherwise scatter might complain
    IF (.NOT.ALLOCATED(flx_co2ff_glo)) ALLOCATE(flx_co2ff_glo(0,0))
    IF (.NOT.ALLOCATED(flx_co2bb_glo)) ALLOCATE(flx_co2bb_glo(0,0))

!--scatter on all proc
    CALL scatter(flx_co2ff_glo,flx_co2ff)
    CALL scatter(flx_co2bb_glo,flx_co2bb)

   IF (ALLOCATED(flx_co2ff_glo)) DEALLOCATE(flx_co2ff_glo)
   IF (ALLOCATED(flx_co2bb_glo)) DEALLOCATE(flx_co2bb_glo)

  ENDIF !--end debuthy

!---select the correct month
  IF (mth_cur.LT.1.OR.mth_cur.GT.12) THEN
    PRINT *,'probleme avec le mois dans co2_ini =', mth_cur
  ENDIF

  fco2_ff(:) = flx_co2ff(:,mth_cur)
  fco2_bb(:) = flx_co2bb(:,mth_cur)

  END SUBROUTINE co2_emissions

END MODULE tracco2i_mod
