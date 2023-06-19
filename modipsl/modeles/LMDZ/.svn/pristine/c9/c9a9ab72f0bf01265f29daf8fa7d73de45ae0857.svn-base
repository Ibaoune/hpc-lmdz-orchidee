!
! $Id: readaerosolstrato1_rrtm.F90 2526 2016-05-26 22:13:40Z oboucher $
!

SUBROUTINE readaerosolstrato1_rrtm(debut)

    USE netcdf95, ONLY: nf95_close, nf95_gw_var, nf95_inq_dimid, & 
                        nf95_inq_varid, nf95_open
    USE netcdf, ONLY: nf90_get_var, nf90_noerr, nf90_nowrite

    USE phys_cal_mod, ONLY : mth_cur
    USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, klon_glo, grid2dTo1d_glo, grid_type, unstructured
    USE mod_phys_lmdz_para 
    USE phys_state_var_mod
    USE phys_local_var_mod
    USE aero_mod
    USE dimphy
    USE YOERAD, ONLY : NLW
    USE YOMCST
#ifdef CPP_XIOS
    USE xios
#endif

    IMPLICIT NONE

! Variable input
    LOGICAL debut

! Variables locales
    INTEGER n_lat   ! number of latitudes in the input data
    INTEGER n_lon   ! number of longitudes in the input data
    INTEGER n_lev   ! number of levels in the input data
    INTEGER n_month ! number of months in the input data
    REAL, POINTER:: latitude(:)
    REAL, POINTER:: longitude(:)
    REAL, POINTER:: time(:)
    REAL, POINTER:: lev(:)
    INTEGER k, band, wave, i
    INTEGER, SAVE :: mth_pre=1
!$OMP THREADPRIVATE(mth_pre)

    REAL, ALLOCATABLE, DIMENSION(:,:), SAVE :: tau_aer_strat
!$OMP THREADPRIVATE(tau_aer_strat)

! Champs reconstitues
    REAL, ALLOCATABLE:: tauaerstrat(:, :, :, :)
    REAL, ALLOCATABLE:: tauaerstrat_mois(:, :, :)
    REAL, ALLOCATABLE:: tauaerstrat_mois_glo(:, :)
    REAL, ALLOCATABLE:: tauaerstrat_mpi(:, :)

! For NetCDF:
    INTEGER ncid_in  ! IDs for input files
    INTEGER varid, ncerr

! Stratospheric aerosols optical properties
! alpha_sw_strat over the 6 bands is normalised by the 550 nm extinction coefficient
    REAL, DIMENSION(nbands_sw_rrtm) :: alpha_sw_strat, piz_sw_strat, cg_sw_strat
    DATA alpha_sw_strat/0.8545564, 0.8451642, 0.9821724, 0.8145110, 0.3073565, 7.7966176E-02/
    DATA cg_sw_strat   /0.6997170, 0.6810035, 0.7403592, 0.7562674, 0.6676504, 0.3478689/
    DATA piz_sw_strat  /0.9999998, 0.9999998, 1.000000000, 0.9999958, 0.9977155, 0.4510679/
!
!--diagnostics AOD in the SW
! alpha_sw_strat_wave is *not* normalised by the 550 nm extinction coefficient
    REAL, DIMENSION(nwave_sw) :: alpha_sw_strat_wave
    DATA alpha_sw_strat_wave/3.708007,4.125824,4.136584,3.887478,3.507738/
!
!--diagnostics AOD in the LW at 10 um (not normalised by the 550 nm ext coefficient
    REAL :: alpha_lw_strat_wave(nwave_lw)
    DATA alpha_lw_strat_wave/0.2746812/
!
    REAL, DIMENSION(nbands_lw_rrtm) :: alpha_lw_abs_rrtm
    DATA alpha_lw_abs_rrtm/   8.8340312E-02, 6.9856711E-02, 6.2652975E-02, 5.7188231E-02, &
                              6.3157059E-02, 5.5072524E-02, 5.0571125E-02, 0.1349073, &    
                              0.1381676, 9.6506312E-02, 5.1312990E-02, 2.4256418E-02, &
                              2.7191756E-02, 3.3862915E-02, 1.6132960E-02, 1.4275438E-02/ ! calculated with Mie_SW_LW_RRTM_V2.4 (bimodal, corrected)
                                                                                          ! for r_0=/0.13E-6, 0.41E-6/ m, sigma_g=/1.26, 1.30/
                                                                                          ! order: increasing wavelength!
!--------------------------------------------------------

    IF (.not.ALLOCATED(tau_aer_strat)) ALLOCATE(tau_aer_strat(klon,klev))

!--we only read monthly strat aerosol data
    IF (debut.OR.mth_cur.NE.mth_pre) THEN

!--only root reads the data
    IF (is_mpi_root.AND.is_omp_root) THEN

    IF (nbands_sw_rrtm.NE.6) THEN 
        print *,'nbands_sw_rrtm doit etre egal a 6 dans readaerosolstrat_rrtm'
        STOP
    ENDIF

    CALL nf95_open("taustrat.nc", nf90_nowrite, ncid_in)

    CALL nf95_inq_varid(ncid_in, "LEV", varid)
    CALL nf95_gw_var(ncid_in, varid, lev)
    n_lev = size(lev)
    IF (n_lev.NE.klev) THEN 
       print *,'Le nombre de niveaux n est pas egal a klev'
       STOP
    ENDIF

    CALL nf95_inq_varid(ncid_in, "LAT", varid)
    CALL nf95_gw_var(ncid_in, varid, latitude)
    n_lat = size(latitude)
    print *, 'LAT aerosol strato=', n_lat, latitude

    IF (grid_type/=unstructured) THEN
      IF (n_lat.NE.nbp_lat) THEN 
         print *,'Le nombre de lat n est pas egal a nbp_lat'
         STOP
      ENDIF
    ENDIF
    
    CALL nf95_inq_varid(ncid_in, "LON", varid)
    CALL nf95_gw_var(ncid_in, varid, longitude)
    n_lon = size(longitude)
    print *, 'LON aerosol strato=', n_lon, longitude

    IF (grid_type/=unstructured) THEN
      IF (n_lon.NE.nbp_lon) THEN 
         print *,'Le nombre de lon n est pas egal a nbp_lon'
         STOP
      ENDIF
    ENDIF
    
    
    CALL nf95_inq_varid(ncid_in, "TIME", varid)
    CALL nf95_gw_var(ncid_in, varid, time)
    n_month = size(time)
    print *, 'TIME aerosol strato=', n_month, time
    IF (n_month.NE.12) THEN 
       print *,'Le nombre de month n est pas egal a 12'
       STOP
    ENDIF

    ALLOCATE(tauaerstrat(n_lon, n_lat, n_lev, n_month))
    ALLOCATE(tauaerstrat_mois(n_lon, n_lat, n_lev))
    ALLOCATE(tauaerstrat_mois_glo(klon_glo, n_lev))

!--reading stratospheric AOD at 550 nm
    CALL nf95_inq_varid(ncid_in, "TAUSTRAT", varid)
    ncerr = nf90_get_var(ncid_in, varid, tauaerstrat)
    print *,'code erreur readaerosolstrato=', ncerr, varid

    CALL nf95_close(ncid_in)

!---select the correct month
    IF (mth_cur.LT.1.OR.mth_cur.GT.12) THEN
      print *,'probleme avec le mois dans readaerosolstrat =', mth_cur
    ENDIF
    tauaerstrat_mois(:,:,:) = tauaerstrat(:,:,:,mth_cur)

!---reduce to a klon_glo grid 
    CALL grid2dTo1d_glo(tauaerstrat_mois,tauaerstrat_mois_glo)
    
    ELSE
      ALLOCATE(tauaerstrat_mois(0,0,0))
    ENDIF !--is_mpi_root and is_omp_root

!$OMP BARRIER

!--keep memory of previous month
    mth_pre=mth_cur

!--scatter on all proc
    
    IF (grid_type==unstructured) THEN
#ifdef CPP_XIOS
      IF (is_omp_master) THEN
        ALLOCATE(tauaerstrat_mpi(klon_mpi,klev))
        CALL xios_send_field("taustrat_in",tauaerstrat_mois)
        CALL xios_recv_field("taustrat_out",tauaerstrat_mpi)
      ELSE
        ALLOCATE(tauaerstrat_mpi(0,0))
      ENDIF
      CALL scatter_omp(tauaerstrat_mpi,tau_aer_strat)
#endif
    ELSE  
      CALL scatter(tauaerstrat_mois_glo,tau_aer_strat)
    ENDIF
    
    IF (is_mpi_root.AND.is_omp_root) THEN
!
    DEALLOCATE(tauaerstrat)
    DEALLOCATE(tauaerstrat_mois)
    DEALLOCATE(tauaerstrat_mois_glo)
!
    ENDIF !--is_mpi_root and is_omp_root

!$OMP BARRIER

    ENDIF !--debut ou nouveau mois

!--total vertical aod at the 5 SW wavelengths
    DO wave=1, nwave_sw
    DO k=1, klev
      tausum_aero(:,wave,id_STRAT_phy)=tausum_aero(:,wave,id_STRAT_phy)+ &
          tau_aer_strat(:,k)*alpha_sw_strat_wave(wave)/alpha_sw_strat_wave(2)
    ENDDO
    ENDDO

!--weighted average for cg, piz and tau, adding strat aerosols on top of tropospheric ones
    DO band=1, nbands_sw_rrtm
!--anthropogenic aerosols bands 1 to nbands_sw_rrtm 
    cg_aero_sw_rrtm(:,:,2,band)  = ( cg_aero_sw_rrtm(:,:,2,band)*piz_aero_sw_rrtm(:,:,2,band)*tau_aero_sw_rrtm(:,:,2,band) + &
                                  cg_sw_strat(band)*piz_sw_strat(band)*alpha_sw_strat(band)*tau_aer_strat(:,:) ) /           &
                             MAX( piz_aero_sw_rrtm(:,:,2,band)*tau_aero_sw_rrtm(:,:,2,band) +                                &
                                  piz_sw_strat(band)*alpha_sw_strat(band)*tau_aer_strat(:,:), 1.e-15 )
    piz_aero_sw_rrtm(:,:,2,band)  = ( piz_aero_sw_rrtm(:,:,2,band)*tau_aero_sw_rrtm(:,:,2,band) +                            &
                              piz_sw_strat(band)*alpha_sw_strat(band)*tau_aer_strat(:,:) ) /                                 &
                              MAX( tau_aero_sw_rrtm(:,:,2,band) + alpha_sw_strat(band)*tau_aer_strat(:,:), 1.e-15 )
    tau_aero_sw_rrtm(:,:,2,band)  = tau_aero_sw_rrtm(:,:,2,band) + alpha_sw_strat(band)*tau_aer_strat(:,:)
!--natural aerosols bands 1 to nbands_sw_rrtm
    cg_aero_sw_rrtm(:,:,1,band)  = ( cg_aero_sw_rrtm(:,:,1,band)*piz_aero_sw_rrtm(:,:,1,band)*tau_aero_sw_rrtm(:,:,1,band) + &
                             cg_sw_strat(band)*piz_sw_strat(band)*alpha_sw_strat(band)*tau_aer_strat(:,:) ) /                &
                             MAX( piz_aero_sw_rrtm(:,:,1,band)*tau_aero_sw_rrtm(:,:,1,band) +                                &
                                  piz_sw_strat(band)*alpha_sw_strat(band)*tau_aer_strat(:,:), 1.e-15 )
    piz_aero_sw_rrtm(:,:,1,band)  = ( piz_aero_sw_rrtm(:,:,1,band)*tau_aero_sw_rrtm(:,:,1,band) +                            &
                              piz_sw_strat(band)*alpha_sw_strat(band)*tau_aer_strat(:,:) ) /                                 &
                              MAX( tau_aero_sw_rrtm(:,:,1,band) + alpha_sw_strat(band)*tau_aer_strat(:,:),1.e-15 )
    tau_aero_sw_rrtm(:,:,1,band)  = tau_aero_sw_rrtm(:,:,1,band) + alpha_sw_strat(band)*tau_aer_strat(:,:)
!--no stratospheric aerosol in index 1 for these tests
!    cg_aero_sw_rrtm(:,:,1,band)  =  cg_aero_sw_rrtm(:,:,1,band)
!    piz_aero_sw_rrtm(:,:,1,band)  = piz_aero_sw_rrtm(:,:,1,band)
!    tau_aero_sw_rrtm(:,:,1,band)  = tau_aero_sw_rrtm(:,:,1,band)
    ENDDO

!--stratospheric AOD in LW
    IF (nbands_lw_rrtm .NE. NLW) then
      print*, 'different values for NLW (=',NLW,') and nbands_lw_rrtm (=', nbands_lw_rrtm, ')'
      STOP
    ENDIF 

!--total vertical aod at the 1 LW wavelength
    DO wave=1, nwave_lw
    DO k=1, klev
      tausum_aero(:,nwave_sw+wave,id_STRAT_phy)=tausum_aero(:,nwave_sw+wave,id_STRAT_phy)+ &
         tau_aer_strat(:,k)*alpha_lw_strat_wave(wave)/alpha_sw_strat_wave(2)
    ENDDO
    ENDDO

    DO band=1, nbands_lw_rrtm
    tau_aero_lw_rrtm(:,:,2,band)  = tau_aero_lw_rrtm(:,:,2,band) + alpha_lw_abs_rrtm(band)*tau_aer_strat(:,:)
    tau_aero_lw_rrtm(:,:,1,band)  = tau_aero_lw_rrtm(:,:,1,band) + alpha_lw_abs_rrtm(band)*tau_aer_strat(:,:) 
!--no stratospheric aerosols in index 1 for these tests
!    tau_aero_lw_rrtm(:,:,1,band)  = tau_aero_lw_rrtm(:,:,1,band) 
    ENDDO

!--default SSA value if there is no aerosol
!--to avoid 0 values that seems to cause some problem to RRTM
    WHERE (tau_aero_sw_rrtm.LT.1.e-14)
      piz_aero_sw_rrtm = 1.0
    ENDWHERE

!--in principle this should not be necessary 
!--as these variables have min values already but just in case
!--put 1e-15 min value to both SW and LW AOD
    tau_aero_sw_rrtm = MAX(tau_aero_sw_rrtm,1.e-15)
    tau_aero_lw_rrtm = MAX(tau_aero_lw_rrtm,1.e-15)

END SUBROUTINE readaerosolstrato1_rrtm
