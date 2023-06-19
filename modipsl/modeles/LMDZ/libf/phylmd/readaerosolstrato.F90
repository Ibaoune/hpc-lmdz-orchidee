subroutine readaerosolstrato(debut)

    use netcdf95, only: nf95_close, nf95_gw_var, nf95_inq_dimid, & 
                        nf95_inq_varid, nf95_open
    use netcdf, only: nf90_get_var, nf90_noerr, nf90_nowrite

    USE phys_cal_mod, ONLY : mth_cur
    USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, klon_glo, &
                                 grid2dto1d_glo, grid_type, unstructured
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    USE mod_phys_lmdz_omp_data, ONLY :  is_omp_root
    USE mod_phys_lmdz_para 
    USE phys_state_var_mod
    USE phys_local_var_mod
    USE aero_mod
    USE dimphy
    USE print_control_mod, ONLY: prt_level,lunout
#ifdef CPP_XIOS
    USE xios
#endif
    implicit none

    include "YOMCST.h"

! Variable input
    logical debut

! Variables locales
    integer n_lat   ! number of latitudes in the input data
    integer n_lon   ! number of longitudes in the input data
    integer n_lev   ! number of levels in the input data
    integer n_month ! number of months in the input data
    real, pointer:: latitude(:)
    real, pointer:: longitude(:)
    real, pointer:: time(:)
    real, pointer:: lev(:)
    integer i, k, band, wave
    integer, save :: mth_pre=1
!$OMP THREADPRIVATE(mth_pre)

    real, allocatable, dimension(:,:), save :: tau_aer_strat
!$OMP THREADPRIVATE(tau_aer_strat)

! Champs reconstitues
    real, allocatable:: tauaerstrat(:, :, :, :)
    real, allocatable:: tauaerstrat_mois(:, :, :)
    real, allocatable:: tauaerstrat_mois_glo(:, :)
    real, allocatable:: tau_aer_strat_mpi(:, :)

! For NetCDF:
    integer ncid_in  ! IDs for input files
    integer varid, ncerr

! Stratospheric aerosols optical properties
! alpha_strat over the 2 bands is normalised by the 550 nm extinction coefficient
! alpha_strat_wave is *not* normalised by the 550 nm extinction coefficient
    real, dimension(nbands) :: alpha_strat, piz_strat, cg_strat
    data alpha_strat/0.9922547, 0.7114912 /
    data piz_strat  /0.9999998, 0.99762493/
    data cg_strat   /0.73107845,0.73229635/
    real, dimension(nwave_sw) :: alpha_strat_wave
    data alpha_strat_wave/3.36780953,3.34667683,3.20444202,3.0293026,2.82108808/

    CHARACTER (len = 20)                      :: modname = 'readaerosolstrato'
    CHARACTER (len = 80)                      :: abort_message

!--------------------------------------------------------

    IF (.not.ALLOCATED(tau_aer_strat)) ALLOCATE(tau_aer_strat(klon,klev))

!--only read file if beginning of run or start of new month
    IF (debut.OR.mth_cur.NE.mth_pre) THEN

!--only root reads
    IF (is_mpi_root.AND.is_omp_root) THEN

    IF (nbands.NE.2) THEN 
        abort_message='nbands doit etre egal a 2 dans readaerosolstrat'
        CALL abort_physic(modname,abort_message,1)
    ENDIF

    CALL nf95_open("taustrat.nc", nf90_nowrite, ncid_in)

    CALL nf95_inq_varid(ncid_in, "LEV", varid)
    CALL nf95_gw_var(ncid_in, varid, lev)
    n_lev = size(lev)
    IF (n_lev.NE.klev) THEN 
       abort_message='Le nombre de niveaux n est pas egal a klev'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    CALL nf95_inq_varid(ncid_in, "LAT", varid)
    CALL nf95_gw_var(ncid_in, varid, latitude)
    n_lat = size(latitude)
    WRITE(lunout,*) 'LAT aerosol strato=', n_lat, latitude
    IF (grid_type/=unstructured) THEN
      IF (n_lat.NE.nbp_lat) THEN 
         abort_message='Le nombre de lat n est pas egal a nbp_lat'
         CALL abort_physic(modname,abort_message,1)
      ENDIF
    ENDIF
    
    CALL nf95_inq_varid(ncid_in, "LON", varid)
    CALL nf95_gw_var(ncid_in, varid, longitude)
    n_lon = size(longitude)
    IF (grid_type/=unstructured) THEN
      WRITE(lunout,*) 'LON aerosol strato=', n_lon, longitude
      IF (n_lon.NE.nbp_lon) THEN 
         abort_message='Le nombre de lon n est pas egal a nbp_lon'
         CALL abort_physic(modname,abort_message,1)
      ENDIF
    ENDIF
    
    CALL nf95_inq_varid(ncid_in, "TIME", varid)
    CALL nf95_gw_var(ncid_in, varid, time)
    n_month = size(time)
    WRITE(lunout,*) 'TIME aerosol strato=', n_month, time
    IF (n_month.NE.12) THEN 
       abort_message='Le nombre de month n est pas egal a 12'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    IF (.not.ALLOCATED(tauaerstrat))          ALLOCATE(tauaerstrat(n_lon, n_lat, n_lev, n_month))
    IF (.not.ALLOCATED(tauaerstrat_mois))     ALLOCATE(tauaerstrat_mois(n_lon, n_lat, n_lev))
    IF (.not.ALLOCATED(tauaerstrat_mois_glo)) ALLOCATE(tauaerstrat_mois_glo(klon_glo, n_lev))

!--reading stratospheric AOD at 550 nm
    CALL nf95_inq_varid(ncid_in, "TAUSTRAT", varid)
    ncerr = nf90_get_var(ncid_in, varid, tauaerstrat)
    WRITE(lunout,*) 'code erreur readaerosolstrato=', ncerr, varid

    CALL nf95_close(ncid_in)

!---select the correct month
    IF (mth_cur.LT.1.OR.mth_cur.GT.12) THEN
     WRITE(lunout,*) 'probleme avec le mois dans readaerosolstrat =', mth_cur
    ENDIF
    tauaerstrat_mois(:,:,:) = tauaerstrat(:,:,:,mth_cur)

!---reduce to a klon_glo grid 
    CALL grid2dTo1d_glo(tauaerstrat_mois,tauaerstrat_mois_glo)

    ELSE
      ALLOCATE(tauaerstrat_mois(0,0,0))
    ENDIF !--is_mpi_root and is_omp_root

!$OMP BARRIER

    IF (grid_type==unstructured) THEN
#ifdef CPP_XIOS
      IF (is_omp_master) THEN
        CALL xios_send_field("taustrat_in",tauaerstrat_mois)
        ALLOCATE(tau_aer_strat_mpi(klon_mpi, klev))
        CALL xios_recv_field("taustrat_out",tau_aer_strat_mpi)
      ELSE
        ALLOCATE(tau_aer_strat_mpi(0,0))
      ENDIF
      CALL scatter_omp(tau_aer_strat_mpi,tau_aer_strat)
#endif
    ELSE  
!--scatter on all proc
      CALL scatter(tauaerstrat_mois_glo,tau_aer_strat)
    ENDIF

!--keep memory of previous month
    mth_pre=mth_cur
!
    IF (is_mpi_root.AND.is_omp_root) THEN
!
    DEALLOCATE(tauaerstrat)
    DEALLOCATE(tauaerstrat_mois)
    DEALLOCATE(tauaerstrat_mois_glo)
!
    ENDIF !-is_mpi_root and is_omp_root

!$OMP BARRIER

    ENDIF !--debut ou nouveau mois

!--total vertical aod at the 6 wavelengths
    DO wave=1, nwave_sw
    DO k=1, klev
    tausum_aero(:,wave,id_STRAT_phy)=tausum_aero(:,wave,id_STRAT_phy)+tau_aer_strat(:,k)*alpha_strat_wave(wave)/alpha_strat_wave(2)
    ENDDO
    ENDDO

!--weighted average for cg, piz and tau, adding strat aerosols on top of tropospheric ones
    DO band=1, nbands
!--anthropogenic aerosols bands 1 and 2
    cg_aero(:,:,3,band)  = ( cg_aero(:,:,3,band)*piz_aero(:,:,3,band)*tau_aero(:,:,3,band) +          &
                             cg_strat(band)*piz_strat(band)*alpha_strat(band)*tau_aer_strat(:,:) ) /  &
                             MAX( piz_aero(:,:,3,band)*tau_aero(:,:,3,band) +                         &
                                  piz_strat(band)*alpha_strat(band)*tau_aer_strat(:,:), 1.e-15 )
    piz_aero(:,:,3,band)  = ( piz_aero(:,:,3,band)*tau_aero(:,:,3,band) +                             &
                              piz_strat(band)*alpha_strat(band)*tau_aer_strat(:,:) ) /                &
                              MAX( tau_aero(:,:,3,band) + alpha_strat(band)*tau_aer_strat(:,:), 1.e-15 )
    tau_aero(:,:,3,band)  = tau_aero(:,:,3,band) + alpha_strat(band)*tau_aer_strat(:,:)
!--natural aerosols bands 1 and 2
    cg_aero(:,:,2,band)  = ( cg_aero(:,:,2,band)*piz_aero(:,:,2,band)*tau_aero(:,:,2,band) +          &
                             cg_strat(band)*piz_strat(band)*alpha_strat(band)*tau_aer_strat(:,:) ) /  &
                             MAX( piz_aero(:,:,2,band)*tau_aero(:,:,2,band) +                         &
                                  piz_strat(band)*alpha_strat(band)*tau_aer_strat(:,:), 1.e-15 )
    piz_aero(:,:,2,band)  = ( piz_aero(:,:,2,band)*tau_aero(:,:,2,band) +                             & 
                              piz_strat(band)*alpha_strat(band)*tau_aer_strat(:,:) ) /                &
                              MAX( tau_aero(:,:,2,band) + alpha_strat(band)*tau_aer_strat(:,:),1.e-15 )
    tau_aero(:,:,2,band)  = tau_aero(:,:,2,band) + alpha_strat(band)*tau_aer_strat(:,:)
    ENDDO

end subroutine readaerosolstrato
