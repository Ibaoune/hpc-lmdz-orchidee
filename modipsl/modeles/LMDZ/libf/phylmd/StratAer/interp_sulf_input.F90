!
! $Id: interp_sulf_input.F90 3677 2020-05-06 15:18:32Z oboucher $
!
SUBROUTINE interp_sulf_input(debutphy,pdtphys,paprs,tr_seri)

  USE netcdf95, ONLY: nf95_close, nf95_gw_var, nf95_inq_dimid, & 
                      nf95_inq_varid, nf95_inquire_dimension, nf95_open
  USE netcdf, ONLY: nf90_get_var, nf90_noerr, nf90_nowrite

  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
  USE mod_phys_lmdz_omp_data, ONLY :  is_omp_root
  USE phys_local_var_mod, ONLY : budg_3D_backgr_ocs, budg_3D_backgr_so2
  USE phys_local_var_mod, ONLY : OCS_lifetime, SO2_lifetime
  USE mod_phys_lmdz_para 
  USE dimphy
  USE phys_cal_mod
  USE infotrac_phy
  USE aerophys
  USE YOMCST

  IMPLICIT NONE

  include "dimensions.h"

! Variable input
  REAL paprs(klon,klev+1)
  REAL tr_seri(klon,klev,nbtr)
  REAL, INTENT(IN)    :: pdtphys    ! Pas d'integration pour la physique (seconde)
  LOGICAL, INTENT(IN) :: debutphy   ! le flag de l'initialisation de la physique

! Variables locales
  INTEGER n_lat   ! number of latitudes in the input data
  INTEGER n_lon   ! number of longitudes in the input data
  INTEGER, SAVE :: n_lev   ! number of levels in the input data
!$OMP THREADPRIVATE(n_lev)
  INTEGER n_mth   ! number of months in the input data
  INTEGER, SAVE :: mth_pre
!$OMP THREADPRIVATE(mth_pre)

! Champs reconstitues
  REAL paprs_glo(klon_glo,klev+1)

  REAL, POINTER:: latitude(:)
! (of input data sorted in strictly ascending order)

  REAL, POINTER:: longitude(:)
! (of input data sorted in strictly ascending order)

  REAL, POINTER:: time(:)
! (of input data sorted in strictly ascending order)

  REAL, POINTER:: lev(:)
! levels of input data

  REAL, ALLOCATABLE :: OCS_clim_in(:, :, :, :)
  REAL, ALLOCATABLE :: SO2_clim_in(:, :, :, :)
  REAL, ALLOCATABLE :: OCS_clim_mth(:, :, :)
  REAL, ALLOCATABLE :: SO2_clim_mth(:, :, :)
  REAL, ALLOCATABLE :: OCS_clim_tmp(:, :)
  REAL, ALLOCATABLE :: SO2_clim_tmp(:, :)
  REAL OCS_clim_glo(klon_glo,klev)
  REAL SO2_clim_glo(klon_glo,klev)
  REAL, ALLOCATABLE :: OCS_lifetime_in(:, :, :, :)
  REAL, ALLOCATABLE :: SO2_lifetime_in(:, :, :, :)
  REAL, ALLOCATABLE :: OCS_lifetime_mth(:, :, :)
  REAL, ALLOCATABLE :: SO2_lifetime_mth(:, :, :)
  REAL, ALLOCATABLE :: OCS_lifetime_tmp(:, :)
  REAL, ALLOCATABLE :: SO2_lifetime_tmp(:, :)
  REAL OCS_lifetime_glo(klon_glo,klev)
  REAL SO2_lifetime_glo(klon_glo,klev)
!
  REAL, ALLOCATABLE, SAVE :: OCS_clim(:,:)
  REAL, ALLOCATABLE, SAVE :: SO2_clim(:,:)
!$OMP THREADPRIVATE(OCS_clim,SO2_clim)
!
  INTEGER i, k, kk, j
  REAL p_bound

! For NetCDF:
  INTEGER ncid_in  ! IDs for input files
  INTEGER varid, ncerr
    
  INTEGER, PARAMETER :: lev_input=17
!--pressure at interfaces of input data (in Pa)
  REAL, DIMENSION(lev_input+1), PARAMETER ::          & 
                    paprs_input=(/                    &
  1.00000002e+05,   6.06530673e+04,   3.67879449e+04, &
  2.23130165e+04,   1.35335286e+04,   8.20850004e+03, &
  4.97870695e+03,   3.01973841e+03,   1.83156393e+03, &
  1.11089968e+03,   6.73794715e+02,   4.08677153e+02, &
  2.47875223e+02,   1.50343923e+02,   9.11881985e+01, &
  5.53084382e+01,   3.35462635e+01,   0.0           /)
!
 IF (.NOT.ALLOCATED(OCS_clim)) ALLOCATE(OCS_clim(klon,klev))
 IF (.NOT.ALLOCATED(SO2_clim)) ALLOCATE(SO2_clim(klon,klev))

  IF (debutphy.OR.mth_cur.NE.mth_pre) THEN

!--preparation of global fields
  CALL gather(paprs, paprs_glo)

  IF (is_mpi_root.AND.is_omp_root) THEN

!--reading emission files
    CALL nf95_open("ocs_so2_annual_lmdz.nc", nf90_nowrite, ncid_in)

    CALL nf95_inq_varid(ncid_in, "LEV", varid)
    CALL nf95_gw_var(ncid_in, varid, lev)
    n_lev = size(lev)

    CALL nf95_inq_varid(ncid_in, "lat", varid)
    CALL nf95_gw_var(ncid_in, varid, latitude)
    n_lat = size(latitude)

    CALL nf95_inq_varid(ncid_in, "lon", varid)
    CALL nf95_gw_var(ncid_in, varid, longitude)
    n_lon = size(longitude)

    CALL nf95_inq_varid(ncid_in, "TIME", varid)
    CALL nf95_gw_var(ncid_in, varid, time)
    n_mth = size(time)

    IF (.NOT.ALLOCATED(OCS_clim_in))     ALLOCATE(OCS_clim_in(n_lon, n_lat, n_lev, n_mth))
    IF (.NOT.ALLOCATED(SO2_clim_in))     ALLOCATE(SO2_clim_in(n_lon, n_lat, n_lev, n_mth))
    IF (.NOT.ALLOCATED(OCS_lifetime_in)) ALLOCATE(OCS_lifetime_in(n_lon, n_lat, n_lev, n_mth))
    IF (.NOT.ALLOCATED(SO2_lifetime_in)) ALLOCATE(SO2_lifetime_in(n_lon, n_lat, n_lev, n_mth))

    CALL nf95_inq_varid(ncid_in, "OCS", varid)
    ncerr = nf90_get_var(ncid_in, varid, OCS_clim_in)
    print *,'code erreur OCS=', ncerr, varid

    CALL nf95_inq_varid(ncid_in, "SO2", varid)
    ncerr = nf90_get_var(ncid_in, varid, SO2_clim_in)
    print *,'code erreur SO2=', ncerr, varid

    CALL nf95_inq_varid(ncid_in, "OCS_LIFET", varid)
    ncerr = nf90_get_var(ncid_in, varid, OCS_lifetime_in)
    print *,'code erreur OCS_lifetime_in=', ncerr, varid

    CALL nf95_inq_varid(ncid_in, "SO2_LIFET", varid)
    ncerr = nf90_get_var(ncid_in, varid, SO2_lifetime_in)
    print *,'code erreur SO2_lifetime_in=', ncerr, varid

    CALL nf95_close(ncid_in)

    IF (.NOT.ALLOCATED(OCS_clim_mth)) ALLOCATE(OCS_clim_mth(n_lon, n_lat, n_lev))
    IF (.NOT.ALLOCATED(SO2_clim_mth)) ALLOCATE(SO2_clim_mth(n_lon, n_lat, n_lev))
    IF (.NOT.ALLOCATED(OCS_clim_tmp)) ALLOCATE(OCS_clim_tmp(klon_glo, n_lev))
    IF (.NOT.ALLOCATED(SO2_clim_tmp)) ALLOCATE(SO2_clim_tmp(klon_glo, n_lev))
    IF (.NOT.ALLOCATED(OCS_lifetime_mth)) ALLOCATE(OCS_lifetime_mth(n_lon, n_lat, n_lev))
    IF (.NOT.ALLOCATED(SO2_lifetime_mth)) ALLOCATE(SO2_lifetime_mth(n_lon, n_lat, n_lev))
    IF (.NOT.ALLOCATED(OCS_lifetime_tmp)) ALLOCATE(OCS_lifetime_tmp(klon_glo, n_lev))
    IF (.NOT.ALLOCATED(SO2_lifetime_tmp)) ALLOCATE(SO2_lifetime_tmp(klon_glo, n_lev))

!---select the correct month, undo multiplication with 1.e12 (precision reasons)
!---correct latitudinal order and convert input from volume mixing ratio to mass mixing ratio
    DO j=1,n_lat
      SO2_clim_mth(:,j,:) = 1.e-12*SO2_clim_in(:,n_lat+1-j,:,mth_cur)*mSO2mol/mAIRmol
      OCS_clim_mth(:,j,:) = 1.e-12*OCS_clim_in(:,n_lat+1-j,:,mth_cur)*mOCSmol/mAIRmol
      SO2_lifetime_mth(:,j,:) = SO2_lifetime_in(:,n_lat+1-j,:,mth_cur)
      OCS_lifetime_mth(:,j,:) = OCS_lifetime_in(:,n_lat+1-j,:,mth_cur)
    ENDDO

!---reduce to a klon_glo grid but keep the levels
    CALL grid2dTo1d_glo(OCS_clim_mth,OCS_clim_tmp)
    CALL grid2dTo1d_glo(SO2_clim_mth,SO2_clim_tmp)
    CALL grid2dTo1d_glo(OCS_lifetime_mth,OCS_lifetime_tmp)
    CALL grid2dTo1d_glo(SO2_lifetime_mth,SO2_lifetime_tmp)

!--set lifetime to very high value in uninsolated areas
    DO i=1, klon_glo
      DO kk=1, n_lev
        IF (OCS_lifetime_tmp(i,kk)==0.0) THEN
          OCS_lifetime_tmp(i,kk)=1.0e12
        ENDIF
        IF (SO2_lifetime_tmp(i,kk)==0.0) THEN
          SO2_lifetime_tmp(i,kk)=1.0e12
        ENDIF
      ENDDO
    ENDDO

  !---regrid weighted lifetime and climatologies
  DO i=1, klon_glo
    DO k=1, klev
     OCS_lifetime_glo(i,k)=0.0
     SO2_lifetime_glo(i,k)=0.0
     OCS_clim_glo(i,k)=0.0
     SO2_clim_glo(i,k)=0.0
     DO kk=1, n_lev
      OCS_lifetime_glo(i,k)=OCS_lifetime_glo(i,k)+ &
           MAX(0.0,MIN(paprs_glo(i,k),paprs_input(kk))-MAX(paprs_glo(i,k+1),paprs_input(kk+1))) &
           *OCS_lifetime_tmp(i,kk)/(paprs_glo(i,k)-paprs_glo(i,k+1))
      SO2_lifetime_glo(i,k)=SO2_lifetime_glo(i,k)+ &
           MAX(0.0,MIN(paprs_glo(i,k),paprs_input(kk))-MAX(paprs_glo(i,k+1),paprs_input(kk+1))) &
           *SO2_lifetime_tmp(i,kk)/(paprs_glo(i,k)-paprs_glo(i,k+1))
      OCS_clim_glo(i,k)=OCS_clim_glo(i,k)+ &
           MAX(0.0,MIN(paprs_glo(i,k),paprs_input(kk))-MAX(paprs_glo(i,k+1),paprs_input(kk+1))) &
           *OCS_clim_tmp(i,kk)/(paprs_glo(i,k)-paprs_glo(i,k+1))
      SO2_clim_glo(i,k)=SO2_clim_glo(i,k)+ &
           MAX(0.0,MIN(paprs_glo(i,k),paprs_input(kk))-MAX(paprs_glo(i,k+1),paprs_input(kk+1))) &
           *SO2_clim_tmp(i,kk)/(paprs_glo(i,k)-paprs_glo(i,k+1))
      ENDDO
    ENDDO
  ENDDO

  ENDIF !--is_mpi_root and is_omp_root

!--keep memory of previous month
  mth_pre=mth_cur

!--scatter global fields around
  CALL scatter(OCS_clim_glo, OCS_clim)
  CALL scatter(SO2_clim_glo, SO2_clim)
  CALL scatter(OCS_lifetime_glo, OCS_lifetime)
  CALL scatter(SO2_lifetime_glo, SO2_lifetime)

  IF (is_mpi_root.AND.is_omp_root) THEN
!
    DEALLOCATE(OCS_clim_in,SO2_clim_in)
    DEALLOCATE(OCS_clim_mth,SO2_clim_mth)
    DEALLOCATE(OCS_clim_tmp,SO2_clim_tmp)
    DEALLOCATE(OCS_lifetime_in,SO2_lifetime_in)
    DEALLOCATE(OCS_lifetime_mth,SO2_lifetime_mth)
    DEALLOCATE(OCS_lifetime_tmp,SO2_lifetime_tmp)
!
  ENDIF !--is_mpi_root and is_omp_root

  ENDIF ! debutphy.OR.new month

!--set to background value everywhere in the very beginning, later only in the troposphere
!--a little dangerous as the MAXVAL is not computed on the global field
  IF (debutphy.AND.MAXVAL(tr_seri).LT.1.e-30) THEN
    p_bound=0.0
  ELSE
    p_bound=50000.
  ENDIF

!--regridding tracer concentration on the vertical
  DO i=1, klon
    DO k=1, klev
      !
      !--OCS and SO2 prescribed back to their clim values below p_bound
      IF (paprs(i,k).GT.p_bound) THEN
        budg_3D_backgr_ocs(i,k)=OCS_clim(i,k)-tr_seri(i,k,id_OCS_strat)
        budg_3D_backgr_so2(i,k)=SO2_clim(i,k)-tr_seri(i,k,id_SO2_strat)
        tr_seri(i,k,id_OCS_strat)=OCS_clim(i,k)
        tr_seri(i,k,id_SO2_strat)=SO2_clim(i,k)
      ENDIF
    ENDDO
  ENDDO

  !convert SO2_backgr_tend from kg(SO2)/kgA to kg(S)/m2/layer/s for saving as diagnostic
  DO i=1, klon
    DO k=1, klev
      budg_3D_backgr_ocs(i,k)=budg_3D_backgr_ocs(i,k)*mSatom/mOCSmol*(paprs(i,k)-paprs(i,k+1))/RG/pdtphys
      budg_3D_backgr_so2(i,k)=budg_3D_backgr_so2(i,k)*mSatom/mSO2mol*(paprs(i,k)-paprs(i,k+1))/RG/pdtphys
    ENDDO
  ENDDO
 
  RETURN

END SUBROUTINE interp_sulf_input
