!
! $Id$
!
!--This routine is to be tested with MPI / OMP parallelism 
!--OB 26/03/2018

SUBROUTINE readchlorophyll(debut)

    USE netcdf95, ONLY: nf95_close, nf95_gw_var, nf95_inq_dimid, nf95_inq_varid, nf95_open
    USE netcdf, ONLY: nf90_get_var, nf90_noerr, nf90_nowrite
    USE phys_cal_mod, ONLY: mth_cur
    USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, klon_glo, grid2dto1d_glo
    USE mod_phys_lmdz_mpi_data, ONLY: is_mpi_root
    USE mod_phys_lmdz_omp_data, ONLY: is_omp_root
    USE mod_phys_lmdz_para, ONLY: scatter 
    USE phys_state_var_mod, ONLY: chl_con
    USE print_control_mod, ONLY: prt_level,lunout

    IMPLICIT NONE

    INCLUDE "YOMCST.h"

! Variable input
    LOGICAL debut

! Variables locales
    INTEGER n_lat   ! number of latitudes in the input data
    INTEGER n_lon   ! number of longitudes in the input data
    INTEGER n_lev   ! number of levels in the input data
    INTEGER n_month ! number of months in the input data
    REAL, POINTER :: latitude(:)
    REAL, POINTER :: longitude(:)
    REAL, POINTER :: time(:)
    INTEGER i, k
    INTEGER, SAVE :: mth_pre
!$OMP THREADPRIVATE(mth_pre)

! Champs reconstitues
    REAL, ALLOCATABLE :: chlorocon(:, :, :)
    REAL, ALLOCATABLE :: chlorocon_mois(:, :)
    REAL, ALLOCATABLE :: chlorocon_mois_glo(:)

! For NetCDF:
    INTEGER ncid_in  ! IDs for input files
    INTEGER varid, ncerr

!--------------------------------------------------------
    CHARACTER (len = 20)  :: modname = 'readchlorophyll'
    CHARACTER (len = 80)  :: abort_message

!--only read file if beginning of run or start of new month
    IF (debut.OR.mth_cur.NE.mth_pre) THEN

    IF (is_mpi_root.AND.is_omp_root) THEN

    CALL nf95_open("chlorophyll.nc", nf90_nowrite, ncid_in)

    CALL nf95_inq_varid(ncid_in, "lon", varid)
    CALL nf95_gw_var(ncid_in, varid, longitude)
    n_lon = size(longitude)
    IF (n_lon.NE.nbp_lon) THEN
       abort_message='Le nombre de lon n est pas egal a nbp_lon'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    CALL nf95_inq_varid(ncid_in, "lat", varid)
    CALL nf95_gw_var(ncid_in, varid, latitude)
    n_lat = size(latitude)
    IF (n_lat.NE.nbp_lat) THEN 
       abort_message='Le nombre de lat n est pas egal a jnbp_lat'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    CALL nf95_inq_varid(ncid_in, "time", varid)
    CALL nf95_gw_var(ncid_in, varid, time)
    n_month = size(time)
    IF (n_month.NE.12) THEN 
       abort_message='Le nombre de month n est pas egal a 12'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    IF (.not.ALLOCATED(chlorocon))          ALLOCATE(chlorocon(n_lon, n_lat, n_month))
    IF (.not.ALLOCATED(chlorocon_mois))     ALLOCATE(chlorocon_mois(n_lon, n_lat))
    IF (.not.ALLOCATED(chlorocon_mois_glo)) ALLOCATE(chlorocon_mois_glo(klon_glo))

!--reading stratospheric AOD at 550 nm
    CALL nf95_inq_varid(ncid_in, "CHL", varid)
    ncerr = nf90_get_var(ncid_in, varid, chlorocon)
    WRITE(lunout,*)'code erreur readchlorophyll=', ncerr, varid

    CALL nf95_close(ncid_in)

!---select the correct month
    IF (mth_cur.LT.1.OR.mth_cur.GT.12) THEN
      WRITE(lunout,*)'probleme avec le mois dans readchlorophyll =', mth_cur
    ENDIF
    chlorocon_mois(:,:) = chlorocon(:,:,mth_cur)

!---reduce to a klon_glo grid 
    CALL grid2dTo1d_glo(chlorocon_mois,chlorocon_mois_glo)

    WRITE(lunout,*)"chrolophyll current month",mth_cur
    DO i=1,klon_glo
!      if(isnan(chlorocon_mois_glo(i)))then ! isnan() is not in the Fortran standard...
!      Another way to check for NaN:
       IF (chlorocon_mois_glo(i).NE.chlorocon_mois_glo(i)) chlorocon_mois_glo(i)=0.
    ENDDO

!    DEALLOCATE(chlorocon)
!    DEALLOCATE(chlorocon_mois)
!    DEALLOCATE(chlorocon_mois_glo)
  
    ENDIF !--is_mpi_root and is_omp_root

!--scatter on all proc
    CALL scatter(chlorocon_mois_glo,chl_con)

!--keep memory of previous month
    mth_pre=mth_cur

    ENDIF !--debut ou nouveau mois

END SUBROUTINE readchlorophyll
