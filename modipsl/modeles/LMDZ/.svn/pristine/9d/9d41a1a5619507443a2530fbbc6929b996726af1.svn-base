!
! $Id: readaerosolstrato2_rrtm.F90 2526 2016-05-26 22:13:40Z oboucher $
!
SUBROUTINE readaerosolstrato2_rrtm(debut, ok_volcan)

    USE netcdf95, ONLY: nf95_close, nf95_gw_var, nf95_inq_dimid, & 
                        nf95_inq_varid, nf95_open
    USE netcdf, ONLY: nf90_get_var, nf90_noerr, nf90_nowrite

    USE phys_cal_mod, ONLY : mth_cur
    USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, klon_glo, grid2dTo1d_glo, grid_type, unstructured
    USE mod_phys_lmdz_mpi_data
    USE mod_phys_lmdz_omp_data
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

    INCLUDE "clesphys.h"

    CHARACTER (len = 80) :: abort_message
    CHARACTER (LEN=20) :: modname = 'readaerosolstrato2'

! Variable input
    LOGICAL, INTENT(IN) ::  debut
    LOGICAL, INTENT(IN) ::  ok_volcan !activate volcanic diags

! Variables locales
    INTEGER n_lat   ! number of latitudes in the input data
    INTEGER n_lon   ! number of longitudes
    INTEGER n_lev   ! number of levels in the input data
    INTEGER n_month ! number of months in the input data
    INTEGER n_wav   ! number of wavelengths in the input data
    REAL, POINTER:: latitude(:)
    REAL, POINTER:: time(:)
    REAL, POINTER:: lev(:)
    REAL, POINTER:: wav(:)
    INTEGER i,k,wave,band
    INTEGER, SAVE :: mth_pre=1
!$OMP THREADPRIVATE(mth_pre)

    REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: tau_aer_strat
    REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: piz_aer_strat
    REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cg_aer_strat
    REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: taulw_aer_strat
!$OMP THREADPRIVATE(tau_aer_strat,piz_aer_strat,cg_aer_strat,taulw_aer_strat)

! Champs reconstitues
    REAL, ALLOCATABLE:: tauaerstrat(:, :, :, :)
    REAL, ALLOCATABLE:: pizaerstrat(:, :, :, :)
    REAL, ALLOCATABLE:: cgaerstrat(:, :, :, :)
    REAL, ALLOCATABLE:: taulwaerstrat(:, :, :, :)

    REAL, ALLOCATABLE:: tauaerstrat_mois(:, :, :, :)
    REAL, ALLOCATABLE:: pizaerstrat_mois(:, :, :, :)
    REAL, ALLOCATABLE:: cgaerstrat_mois(:, :, :, :)
    REAL, ALLOCATABLE:: taulwaerstrat_mois(:, :, :, :)

    REAL, ALLOCATABLE:: tauaerstrat_mois_glo(:, :, :)
    REAL, ALLOCATABLE:: pizaerstrat_mois_glo(:, :, :)
    REAL, ALLOCATABLE:: cgaerstrat_mois_glo(:, :, :)
    REAL, ALLOCATABLE:: taulwaerstrat_mois_glo(:, :, :)
    REAL, ALLOCATABLE:: tauaerstrat_mpi(:, :, :)
    REAL, ALLOCATABLE:: pizaerstrat_mpi(:, :, :)
    REAL, ALLOCATABLE:: cgaerstrat_mpi(:, :, :)
    REAL, ALLOCATABLE:: taulwaerstrat_mpi(:, :, :)

! For NetCDF:
    INTEGER ncid_in  ! IDs for input files
    INTEGER varid, ncerr

!--------------------------------------------------------

    IF (.not.ALLOCATED(tau_aer_strat)) ALLOCATE(tau_aer_strat(klon,klev,NSW))
    IF (.not.ALLOCATED(piz_aer_strat)) ALLOCATE(piz_aer_strat(klon,klev,NSW))
    IF (.not.ALLOCATED(cg_aer_strat))  ALLOCATE(cg_aer_strat(klon,klev,NSW))

    IF (.not.ALLOCATED(taulw_aer_strat)) ALLOCATE(taulw_aer_strat(klon,klev,NLW))

!--we only read monthly strat aerosol data
    IF (debut.OR.mth_cur.NE.mth_pre) THEN

!--only root reads the data
      IF (is_mpi_root.AND.is_omp_root) THEN

!--check mth_cur
        IF (mth_cur.LT.1.OR.mth_cur.GT.12) THEN
          print *,'probleme avec le mois dans readaerosolstrat =', mth_cur
        ENDIF

!--initialize n_lon as input data is 2D (lat-alt) only
        n_lon = nbp_lon

!--Starts with SW optical properties

        CALL nf95_open("tauswstrat.2D.nc", nf90_nowrite, ncid_in)

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

        IF (grid_type/=unstructured) THEN
           IF (n_lat.NE.nbp_lat) THEN 
             print *, 'latitude=', n_lat, nbp_lat
             abort_message='Le nombre de lat n est pas egal a nbp_lat'
             CALL abort_physic(modname,abort_message,1)
           ENDIF
        ENDIF

        CALL nf95_inq_varid(ncid_in, "TIME", varid)
        CALL nf95_gw_var(ncid_in, varid, time)
        n_month = size(time)
        IF (n_month.NE.12) THEN 
           abort_message='Le nombre de month n est pas egal a 12'
           CALL abort_physic(modname,abort_message,1)
        ENDIF

        CALL nf95_inq_varid(ncid_in, "WAV", varid)
        CALL nf95_gw_var(ncid_in, varid, wav)
        n_wav = size(wav)
        print *, 'WAV aerosol strato=', n_wav, wav
        IF (n_wav.NE.NSW) THEN 
           abort_message='Le nombre de wav n est pas egal a NSW'
           CALL abort_physic(modname,abort_message,1)
        ENDIF

        ALLOCATE(tauaerstrat(n_lat, n_lev, n_wav, n_month))
        ALLOCATE(pizaerstrat(n_lat, n_lev, n_wav, n_month))
        ALLOCATE(cgaerstrat(n_lat, n_lev, n_wav, n_month))

!--reading stratospheric aerosol tau per layer
        CALL nf95_inq_varid(ncid_in, "TAU_SUN", varid)
        ncerr = nf90_get_var(ncid_in, varid, tauaerstrat)
        print *,'code erreur readaerosolstrato=', ncerr, varid

!--reading stratospheric aerosol omega per layer
        CALL nf95_inq_varid(ncid_in, "OME_SUN", varid)
        ncerr = nf90_get_var(ncid_in, varid, pizaerstrat)
        print *,'code erreur readaerosolstrato=', ncerr, varid

!--reading stratospheric aerosol g per layer
        CALL nf95_inq_varid(ncid_in, "GGG_SUN", varid)
        ncerr = nf90_get_var(ncid_in, varid, cgaerstrat)
        print *,'code erreur readaerosolstrato sw=', ncerr, varid

        CALL nf95_close(ncid_in)

       
        IF (grid_type/=unstructured) THEN
          ALLOCATE(tauaerstrat_mois(n_lon, n_lat, n_lev, n_wav))
          ALLOCATE(pizaerstrat_mois(n_lon, n_lat, n_lev, n_wav))
          ALLOCATE(cgaerstrat_mois(n_lon, n_lat, n_lev, n_wav))

          ALLOCATE(tauaerstrat_mois_glo(klon_glo, n_lev, n_wav))
          ALLOCATE(pizaerstrat_mois_glo(klon_glo, n_lev, n_wav))
          ALLOCATE(cgaerstrat_mois_glo(klon_glo, n_lev, n_wav))
!--select the correct month
!--and copy into 1st longitude
          tauaerstrat_mois(1,:,:,:) = tauaerstrat(:,:,:,mth_cur)
          pizaerstrat_mois(1,:,:,:) = pizaerstrat(:,:,:,mth_cur)
          cgaerstrat_mois(1,:,:,:)  = cgaerstrat(:,:,:,mth_cur)

!--copy longitudes
          DO i=2, n_lon
           tauaerstrat_mois(i,:,:,:) = tauaerstrat_mois(1,:,:,:)
           pizaerstrat_mois(i,:,:,:) = pizaerstrat_mois(1,:,:,:)
           cgaerstrat_mois(i,:,:,:)  = cgaerstrat_mois(1,:,:,:)
          ENDDO

!---reduce to a klon_glo grid 
          DO band=1, NSW
            CALL grid2dTo1d_glo(tauaerstrat_mois(:,:,:,band),tauaerstrat_mois_glo(:,:,band))
            CALL grid2dTo1d_glo(pizaerstrat_mois(:,:,:,band),pizaerstrat_mois_glo(:,:,band))
            CALL grid2dTo1d_glo(cgaerstrat_mois(:,:,:,band),cgaerstrat_mois_glo(:,:,band))
          ENDDO
        ENDIF
!--Now LW optical properties
!

        CALL nf95_open("taulwstrat.2D.nc", nf90_nowrite, ncid_in)

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

        IF (grid_type/=unstructured) THEN
          IF (n_lat.NE.nbp_lat) THEN 
             abort_message='Le nombre de lat n est pas egal a nbp_lat'
             CALL abort_physic(modname,abort_message,1)
          ENDIF
        ENDIF
        
        CALL nf95_inq_varid(ncid_in, "TIME", varid)
        CALL nf95_gw_var(ncid_in, varid, time)
        n_month = size(time)
        IF (n_month.NE.12) THEN 
           abort_message='Le nombre de month n est pas egal a 12'
           CALL abort_physic(modname,abort_message,1)
        ENDIF

        CALL nf95_inq_varid(ncid_in, "WAV", varid)
        CALL nf95_gw_var(ncid_in, varid, wav)
        n_wav = size(wav)
        print *, 'WAV aerosol strato=', n_wav, wav
        IF (n_wav.NE.NLW) THEN 
           abort_message='Le nombre de wav n est pas egal a NLW'
           CALL abort_physic(modname,abort_message,1)
        ENDIF

        ALLOCATE(taulwaerstrat(n_lat, n_lev, n_wav, n_month))

!--reading stratospheric aerosol lw tau per layer
        CALL nf95_inq_varid(ncid_in, "TAU_EAR", varid)
        ncerr = nf90_get_var(ncid_in, varid, taulwaerstrat)
        print *,'code erreur readaerosolstrato lw=', ncerr, varid

        CALL nf95_close(ncid_in)

        IF (grid_type/=unstructured) THEN

          ALLOCATE(taulwaerstrat_mois(n_lon, n_lat, n_lev, n_wav))
          ALLOCATE(taulwaerstrat_mois_glo(klon_glo, n_lev, n_wav))

!--select the correct month
!--and copy into 1st longitude
          taulwaerstrat_mois(1,:,:,:) = taulwaerstrat(:,:,:,mth_cur)
!--copy longitudes
          DO i=2, n_lon
            taulwaerstrat_mois(i,:,:,:) = taulwaerstrat_mois(1,:,:,:)
          ENDDO

!---reduce to a klon_glo grid 
          DO band=1, NLW
            CALL grid2dTo1d_glo(taulwaerstrat_mois(:,:,:,band),taulwaerstrat_mois_glo(:,:,band))
          ENDDO
        ENDIF
      
      ELSE !--proc other than mpi_root and omp_root
           !--dummy allocation needed for debug mode

        ALLOCATE(tauaerstrat_mois_glo(1,1,1))
        ALLOCATE(pizaerstrat_mois_glo(1,1,1))
        ALLOCATE(cgaerstrat_mois_glo(1,1,1))
        ALLOCATE(taulwaerstrat_mois_glo(1,1,1))

        ALLOCATE(tauaerstrat(0,0,0,12))
        ALLOCATE(pizaerstrat(0,0,0,12))
        ALLOCATE(cgaerstrat(0,0,0,12))
        ALLOCATE(taulwaerstrat(0,0,0,12))


      ENDIF !--is_mpi_root and is_omp_root

!$OMP BARRIER

!--keep memory of previous month
      mth_pre=mth_cur

      IF (grid_type==unstructured) THEN

#ifdef CPP_XIOS

        IF (is_omp_master) THEN
          ALLOCATE(tauaerstrat_mpi(klon_mpi, klev, NSW))
          ALLOCATE(pizaerstrat_mpi(klon_mpi, klev, NSW))
          ALLOCATE(cgaerstrat_mpi(klon_mpi, klev, NSW))        
          ALLOCATE(taulwaerstrat_mpi(klon_mpi, klev, NLW))
         
          CALL xios_send_field("tauaerstrat_in",SPREAD(tauaerstrat(:,:,:,mth_cur),1,8))
          CALL xios_recv_field("tauaerstrat_out",tauaerstrat_mpi)
          CALL xios_send_field("pizaerstrat_in",SPREAD(pizaerstrat(:,:,:,mth_cur),1,8))
          CALL xios_recv_field("pizaerstrat_out",pizaerstrat_mpi)
          CALL xios_send_field("cgaerstrat_in",SPREAD(cgaerstrat(:,:,:,mth_cur),1,8))
          CALL xios_recv_field("cgaerstrat_out",cgaerstrat_mpi)
          CALL xios_send_field("taulwaerstrat_in",SPREAD(taulwaerstrat(:,:,:,mth_cur),1,8))
          CALL xios_recv_field("taulwaerstrat_out",taulwaerstrat_mpi)
        ELSE
          ALLOCATE(tauaerstrat_mpi(0, 0, 0))
          ALLOCATE(pizaerstrat_mpi(0, 0, 0))
          ALLOCATE(cgaerstrat_mpi(0, 0, 0))        
          ALLOCATE(taulwaerstrat_mpi(0, 0, 0))
        ENDIF  
        
        CALL scatter_omp(tauaerstrat_mpi,tau_aer_strat)
        CALL scatter_omp(pizaerstrat_mpi,piz_aer_strat)
        CALL scatter_omp(cgaerstrat_mpi,cg_aer_strat)
        CALL scatter_omp(taulwaerstrat_mpi,taulw_aer_strat)
#endif
      ELSE  
        
!--scatter on all proc
        CALL scatter(tauaerstrat_mois_glo,tau_aer_strat)
        CALL scatter(pizaerstrat_mois_glo,piz_aer_strat)
        CALL scatter(cgaerstrat_mois_glo,cg_aer_strat)
        CALL scatter(taulwaerstrat_mois_glo,taulw_aer_strat)
        IF (is_mpi_root.AND.is_omp_root)  DEALLOCATE(tauaerstrat_mois, pizaerstrat_mois, cgaerstrat_mois, taulwaerstrat_mois)

      ENDIF 

      IF (is_mpi_root.AND.is_omp_root) THEN
        DEALLOCATE(tauaerstrat, pizaerstrat, cgaerstrat,taulwaerstrat)
      ENDIF
      

!$OMP BARRIER

    ENDIF !--debut ou nouveau mois

!--total vertical aod at the 5 SW wavelengths
!--for now use band 3 AOD into all 5 wavelengths
!--it is only a reasonable approximation for 550 nm (wave=2)
    band=3
    DO i=1, klon
    DO k=1, klev
      IF (stratomask(i,k).GT.0.999999) THEN
        DO wave=1, nwave_sw
          tausum_aero(i,wave,id_STRAT_phy)=tausum_aero(i,wave,id_STRAT_phy)+tau_aer_strat(i,k,band)
        ENDDO
      ENDIF
    ENDDO
    ENDDO

    IF (.NOT. ok_volcan) THEN
!
!--this is the default case 
!--stratospheric aerosols are added to both index 2 and 1 for double radiation calls
!--weighted average for cg, piz and tau, adding strat aerosols on top of tropospheric ones
    DO band=1, NSW
      WHERE (stratomask.GT.0.999999)
!--strat aerosols are added to index 2 : natural and anthropogenic aerosols for bands 1 to NSW 
        cg_aero_sw_rrtm(:,:,2,band)  = ( cg_aero_sw_rrtm(:,:,2,band)*piz_aero_sw_rrtm(:,:,2,band)*tau_aero_sw_rrtm(:,:,2,band) + &
                                         cg_aer_strat(:,:,band)*piz_aer_strat(:,:,band)*tau_aer_strat(:,:,band) ) /              &
                                    MAX( piz_aero_sw_rrtm(:,:,2,band)*tau_aero_sw_rrtm(:,:,2,band) +                             &
                                         piz_aer_strat(:,:,band)*tau_aer_strat(:,:,band), 1.e-15 )
        piz_aero_sw_rrtm(:,:,2,band) = ( piz_aero_sw_rrtm(:,:,2,band)*tau_aero_sw_rrtm(:,:,2,band) +                             &
                                         piz_aer_strat(:,:,band)*tau_aer_strat(:,:,band) ) /                                     &
                                    MAX( tau_aero_sw_rrtm(:,:,2,band) + tau_aer_strat(:,:,band), 1.e-15 )
        tau_aero_sw_rrtm(:,:,2,band)  = tau_aero_sw_rrtm(:,:,2,band) + tau_aer_strat(:,:,band)
!--strat aerosols are added to index 1 : natural aerosols only for bands 1 to NSW
        cg_aero_sw_rrtm(:,:,1,band)  = ( cg_aero_sw_rrtm(:,:,1,band)*piz_aero_sw_rrtm(:,:,1,band)*tau_aero_sw_rrtm(:,:,1,band) + &
                cg_aer_strat(:,:,band)*piz_aer_strat(:,:,band)*tau_aer_strat(:,:,band) ) /              &
                MAX( piz_aero_sw_rrtm(:,:,1,band)*tau_aero_sw_rrtm(:,:,1,band) +                             &
                piz_aer_strat(:,:,band)*tau_aer_strat(:,:,band), 1.e-15 )
        piz_aero_sw_rrtm(:,:,1,band) = ( piz_aero_sw_rrtm(:,:,1,band)*tau_aero_sw_rrtm(:,:,1,band) +                             &
                piz_aer_strat(:,:,band)*tau_aer_strat(:,:,band) ) /                                     &
                MAX( tau_aero_sw_rrtm(:,:,1,band) + tau_aer_strat(:,:,band), 1.e-15 )
        tau_aero_sw_rrtm(:,:,1,band)  = tau_aero_sw_rrtm(:,:,1,band) + tau_aer_strat(:,:,band)
    ENDWHERE
    ENDDO
!
    ELSE
!
!--this is the VOLMIP case
!--stratospheric aerosols are only added to index 2 in this case 
!--weighted average for cg, piz and tau, adding strat aerosols on top of tropospheric ones
    DO band=1, NSW
      WHERE (stratomask.GT.0.999999)
!--strat aerosols are added to index 2 : natural and anthropogenic aerosols for bands 1 to NSW 
        cg_aero_sw_rrtm(:,:,2,band)  = ( cg_aero_sw_rrtm(:,:,2,band)*piz_aero_sw_rrtm(:,:,2,band)*tau_aero_sw_rrtm(:,:,2,band) + &
                                         cg_aer_strat(:,:,band)*piz_aer_strat(:,:,band)*tau_aer_strat(:,:,band) ) /              &
                                    MAX( piz_aero_sw_rrtm(:,:,2,band)*tau_aero_sw_rrtm(:,:,2,band) +                             &
                                         piz_aer_strat(:,:,band)*tau_aer_strat(:,:,band), 1.e-15 )
        piz_aero_sw_rrtm(:,:,2,band) = ( piz_aero_sw_rrtm(:,:,2,band)*tau_aero_sw_rrtm(:,:,2,band) +                             &
                                         piz_aer_strat(:,:,band)*tau_aer_strat(:,:,band) ) /                                     &
                                    MAX( tau_aero_sw_rrtm(:,:,2,band) + tau_aer_strat(:,:,band), 1.e-15 )
        tau_aero_sw_rrtm(:,:,2,band)  = tau_aero_sw_rrtm(:,:,2,band) + tau_aer_strat(:,:,band)
     ENDWHERE
  ENDDO
  ENDIF

!--total vertical aod at 10 um
!--this is approximated from band 7 of RRTM
    band=7
    DO i=1, klon
    DO k=1, klev
      IF (stratomask(i,k).GT.0.999999) THEN
        DO wave=1, nwave_lw
          tausum_aero(i,nwave_sw+wave,id_STRAT_phy)=tausum_aero(i,nwave_sw+wave,id_STRAT_phy)+taulw_aer_strat(i,k,band)
        ENDDO
      ENDIF
    ENDDO
    ENDDO

    IF (.NOT. ok_volcan) THEN
!--this is the default case 
!--stratospheric aerosols are added to both index 2 and 1
    DO band=1, NLW
      WHERE (stratomask.GT.0.999999)
        tau_aero_lw_rrtm(:,:,2,band)  = tau_aero_lw_rrtm(:,:,2,band) + taulw_aer_strat(:,:,band)
        tau_aero_lw_rrtm(:,:,1,band)  = tau_aero_lw_rrtm(:,:,1,band) + taulw_aer_strat(:,:,band)
      ENDWHERE
    ENDDO
!
    ELSE
!
!--this is the VOLMIP case
    DO band=1, NLW
!--stratospheric aerosols are not added to index 1 
!--and we copy index 2 in index 1 because we want the same dust aerosol LW properties as above
      tau_aero_lw_rrtm(:,:,1,band)  = tau_aero_lw_rrtm(:,:,2,band)
!
      WHERE (stratomask.GT.0.999999)
!--stratospheric aerosols are only added to index 2
        tau_aero_lw_rrtm(:,:,2,band)  = tau_aero_lw_rrtm(:,:,2,band) + taulw_aer_strat(:,:,band)
      ENDWHERE
    ENDDO
    ENDIF

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

END SUBROUTINE readaerosolstrato2_rrtm
