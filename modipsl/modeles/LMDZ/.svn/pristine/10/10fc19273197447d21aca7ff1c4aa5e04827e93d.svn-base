MODULE readaerosolstrato_m


CONTAINS

  SUBROUTINE init_readaerosolstrato(flag_aerosol_strat)
  USE mod_grid_phy_lmdz, ONLY: grid_type, unstructured
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: flag_aerosol_strat
    
    IF (grid_type==unstructured)  THEN 
   
      IF (flag_aerosol_strat == 1) THEN
        CALL init_readaerosolstrato1
      ELSE IF (flag_aerosol_strat == 2) THEN
        CALL init_readaerosolstrato2
      ENDIF
  
    ENDIF
  
  END SUBROUTINE init_readaerosolstrato
  


  SUBROUTINE init_readaerosolstrato1
#ifdef CPP_XIOS
  USE netcdf
  USE netcdf95, ONLY: nf95_close, nf95_gw_var, nf95_inq_dimid, & 
                      nf95_inq_varid, nf95_open
  USE mod_phys_lmdz_para
  USE xios
!  USE YOERAD, ONLY : NLW
  IMPLICIT NONE
  REAL, POINTER:: latitude(:)
  REAL, POINTER:: longitude(:)
  INTEGER :: nlat, nlon
  REAL    :: null_array(0)
  INTEGER :: ncid_in, varid
  
    IF (is_omp_master) THEN  
      IF (is_mpi_root) THEN
        CALL nf95_open("taustrat.nc", nf90_nowrite, ncid_in)
        CALL nf95_inq_varid(ncid_in, "LAT", varid)
        CALL nf95_gw_var(ncid_in, varid, latitude)
        CALL nf95_inq_varid(ncid_in, "LON", varid)
        CALL nf95_gw_var(ncid_in, varid, longitude)
        CALL nf95_close(ncid_in)
        nlat=size(latitude)
        nlon=size(longitude)
      ENDIF
      CALL bcast_mpi(nlat)
      CALL bcast_mpi(nlon)

      IF (is_mpi_root) THEN
        CALL xios_set_domain_attr("domain_taustrat",nj_glo=nlat, nj=nlat, jbegin=0, latvalue_1d=latitude)
        CALL xios_set_domain_attr("domain_taustrat",ni_glo=nlon, ni=nlon, ibegin=0, lonvalue_1d = longitude)
       ELSE
        CALL xios_set_domain_attr("domain_taustrat",nj_glo=nlat, nj=0, jbegin=0, latvalue_1d=null_array )
        CALL xios_set_domain_attr("domain_taustrat",ni_glo=nlon, ni=0, ibegin=0, lonvalue_1d=null_array)
      ENDIF    
  
      CALL xios_set_fieldgroup_attr("aerosol_strato1",enabled=.TRUE.)
    ENDIF
    
#endif
  END SUBROUTINE init_readaerosolstrato1
  
  SUBROUTINE init_readaerosolstrato2
#ifdef CPP_XIOS
  USE netcdf
  USE netcdf95, ONLY: nf95_close, nf95_gw_var, nf95_inq_dimid, & 
                      nf95_inq_varid, nf95_open
  USE mod_phys_lmdz_para
  USE xios
!  USE YOERAD, ONLY : NLW
  IMPLICIT NONE
  REAL, POINTER:: latitude(:)
  REAL, POINTER:: wav(:)
  INTEGER :: nlat,n_wav
  REAL    :: null_array(0)
  INTEGER :: ncid_in, varid

    IF (is_omp_master) THEN    
      IF (is_mpi_root) THEN
        CALL nf95_open("tauswstrat.2D.nc", nf90_nowrite, ncid_in)
        CALL nf95_inq_varid(ncid_in, "LAT", varid)
        CALL nf95_gw_var(ncid_in, varid, latitude)
        CALL nf95_inq_varid(ncid_in, "WAV", varid)
        CALL nf95_gw_var(ncid_in, varid, wav)
        CALL nf95_close(ncid_in)
        nlat=size(latitude)
        n_wav = size(wav)
      ENDIF
      CALL bcast_mpi(nlat)
      CALL bcast_mpi(n_wav)

      IF (is_mpi_root) THEN
        CALL xios_set_domain_attr("domain_tauswstrat",nj_glo=nlat, nj=nlat, jbegin=0, latvalue_1d=latitude)
        CALL xios_set_domain_attr("domain_tauswstrat",ni_glo=8, ni=8, ibegin=0, lonvalue_1d = (/ 0.,45.,90.,135.,180.,225.,270., 315. /))
        DEALLOCATE(latitude)
      ELSE
        CALL xios_set_domain_attr("domain_tauswstrat",nj_glo=nlat, nj=0, jbegin=0, latvalue_1d=null_array )
        CALL xios_set_domain_attr("domain_tauswstrat",ni_glo=8, ni=0, ibegin=0, lonvalue_1d=null_array)
      ENDIF    
      CALL  xios_set_axis_attr("nsw", n_glo=n_wav)
    
      IF (is_mpi_root) THEN
        CALL nf95_open("taulwstrat.2D.nc", nf90_nowrite, ncid_in)
        CALL nf95_inq_varid(ncid_in, "LAT", varid)
        CALL nf95_gw_var(ncid_in, varid, latitude)
        CALL nf95_inq_varid(ncid_in, "WAV", varid)
        CALL nf95_gw_var(ncid_in, varid, wav)
        CALL nf95_close(ncid_in)
        nlat=size(latitude)
        n_wav = size(wav)
      ENDIF
      CALL bcast_mpi(nlat)
      CALL bcast_mpi(n_wav)

      IF (is_mpi_root) THEN
        CALL xios_set_domain_attr("domain_taulwstrat",nj_glo=nlat, nj=nlat, jbegin=0, latvalue_1d=latitude)
        CALL xios_set_domain_attr("domain_taulwstrat",ni_glo=8, ni=8, ibegin=0, lonvalue_1d = (/ 0.,45.,90.,135.,180.,225.,270., 315. /))
        DEALLOCATE(latitude)
      ELSE
        CALL xios_set_domain_attr("domain_taulwstrat",nj_glo=nlat, nj=0, jbegin=0, latvalue_1d=null_array )
        CALL xios_set_domain_attr("domain_taulwstrat",ni_glo=8, ni=0, ibegin=0, lonvalue_1d=null_array)
      ENDIF    
    
      CALL  xios_set_axis_attr("nlw", n_glo=n_wav)
      CALL xios_set_fieldgroup_attr("aerosol_strato2",enabled=.TRUE.)

    ENDIF
#endif    
  END SUBROUTINE init_readaerosolstrato2
  

END MODULE readaerosolstrato_m

