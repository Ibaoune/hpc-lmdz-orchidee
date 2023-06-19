MODULE create_etat0_unstruct_mod






CONTAINS
  
  SUBROUTINE init_create_etat0_unstruct
#ifdef CPP_XIOS
  USE xios
  USE netcdf
  USE mod_phys_lmdz_para
  IMPLICIT NONE
  INTEGER :: file_id, iret
  
   ! for coupling activate ocean fraction reading from file "ocean_fraction.nc"
    IF (is_omp_master) THEN

      IF (NF90_OPEN("ocean_fraction.nc", NF90_NOWRITE, file_id)==NF90_NOERR) THEN
        CALL xios_set_file_attr("frac_ocean",enabled=.TRUE.)
        CALL xios_set_field_attr("mask",field_ref="frac_ocean_read")
        iret=NF90_CLOSE(file_id)
      ELSE IF (NF90_OPEN("land_water_0.05.nc", NF90_NOWRITE, file_id)==NF90_NOERR) THEN
        CALL xios_set_file_attr("land_water",name="land_water_0.05",enabled=.TRUE.)
        CALL xios_set_field_attr("mask",field_ref="land_water")
        iret=NF90_CLOSE(file_id)
      ELSE IF (NF90_OPEN("land_water_0.25.nc", NF90_NOWRITE, file_id)==NF90_NOERR) THEN
        CALL xios_set_file_attr("land_water",name="land_water_0.25",enabled=.TRUE.)
        CALL xios_set_field_attr("mask",field_ref="land_water")
        iret=NF90_CLOSE(file_id)
      ELSE IF (NF90_OPEN("land_water_0.50.nc", NF90_NOWRITE, file_id)==NF90_NOERR) THEN
        CALL xios_set_file_attr("land_water",name="land_water_0.50",enabled=.TRUE.)
        CALL xios_set_field_attr("mask",field_ref="land_water")
        iret=NF90_CLOSE(file_id)
      ENDIF

    ENDIF

#endif
  
  END SUBROUTINE init_create_etat0_unstruct
  
 
  SUBROUTINE create_etat0_unstruct
  USE dimphy
#ifdef CPP_XIOS
  USE xios
  USE infotrac_phy
  USE fonte_neige_mod
  USE pbl_surface_mod
  USE phys_state_var_mod
  USE indice_sol_mod
  USE mod_phys_lmdz_para
  USE print_control_mod, ONLY: lunout
  USE geometry_mod
  USE ioipsl_getin_p_mod, ONLY: getin_p

  IMPLICIT NONE
  INCLUDE 'dimsoil.h'

    LOGICAL :: no_ter_antartique   ! If true, no land points are allowed at Antartic
    REAL,    DIMENSION(klon)                 :: tsol
    REAL,    DIMENSION(klon)                 :: sn
    REAL,    DIMENSION(klon)                 :: rugmer
    REAL,    DIMENSION(klon)                 :: run_off_lic_0
    REAL,    DIMENSION(klon)                 :: lic
    REAL,    DIMENSION(klon)                 :: fder

    REAL,    DIMENSION(klon,nbsrf)           :: qsolsrf, snsrf
    REAL,    DIMENSION(klon,nsoilmx,nbsrf)   :: tsoil
    
    REAL,    DIMENSION(klon_mpi)             :: tsol_mpi, qsol_mpi, zmasq_mpi, lic_mpi
    REAL,    DIMENSION(klon_mpi)             :: zmea_mpi, zstd_mpi, zsig_mpi, zgam_mpi, zthe_mpi
    REAL,    DIMENSION(klon_mpi)             :: cell_area_mpi
    REAL,    DIMENSION(klon_mpi,nbsrf)       :: pctsrf_mpi

    INTEGER :: ji,j,i
 
    IF (is_omp_master) THEN
      CALL xios_recv_field("ts",tsol_mpi)
      CALL xios_recv_field("qs",qsol_mpi)
      CALL xios_recv_field("mask",zmasq_mpi)
      CALL xios_recv_field("landice",lic_mpi)
      CALL xios_recv_field("zmea",zmea_mpi)
      CALL xios_recv_field("zstd",zstd_mpi)
      CALL xios_recv_field("zsig",zsig_mpi)
      CALL xios_recv_field("zgam",zgam_mpi)
      CALL xios_recv_field("zthe",zthe_mpi)
    ENDIF
    CALL scatter_omp(tsol_mpi,tsol)
    CALL scatter_omp(qsol_mpi,qsol)
    CALL scatter_omp(zmasq_mpi,zmasq)
    CALL scatter_omp(lic_mpi,lic)
    CALL scatter_omp(zmea_mpi,zmea)
    CALL scatter_omp(zstd_mpi,zstd)
    CALL scatter_omp(zsig_mpi,zsig)
    CALL scatter_omp(zgam_mpi,zgam)
    CALL scatter_omp(zthe_mpi,zthe)

    radsol(:)   = 0.0
    rugmer(:) = 0.001
    sn(:)     = 0

    WHERE(qsol(:)<0) qsol(:)=0 
        
    WHERE(   zmasq(:)<EPSFRA) zmasq(:)=0.
    WHERE(1.-zmasq(:)<EPSFRA) zmasq(:)=1.

    pctsrf(:,:) = 0
    pctsrf(:,is_lic)=lic
    WHERE(pctsrf(:,is_lic)<EPSFRA) pctsrf(:,is_lic)=0. 
    WHERE(zmasq(:)<EPSFRA)         pctsrf(:,is_lic)=0.

    pctsrf(:,is_ter)=zmasq(:)

!--- Adequation with soil/sea mask
    DO ji=1,klon
      IF(zmasq(ji)>EPSFRA) THEN 
        IF(pctsrf(ji,is_lic)>=zmasq(ji)) THEN
          pctsrf(ji,is_lic)=zmasq(ji)
          pctsrf(ji,is_ter)=0.
        ELSE
          pctsrf(ji,is_ter)=zmasq(ji)-pctsrf(ji,is_lic)
          IF(pctsrf(ji,is_ter)<EPSFRA) THEN
            pctsrf(ji,is_ter)=0.
            pctsrf(ji,is_lic)=zmasq(ji)
          END IF 
        END IF 
      END IF 
    END DO 


  !--- Option no_ter_antartique removes all land fractions souther than 60S.
  !--- Land ice is set instead of the land fractions on these latitudes.
  !--- The ocean and sea-ice fractions are not changed.
  no_ter_antartique=.FALSE.
  CALL getin_p('no_ter_antartique',no_ter_antartique)
  WRITE(lunout,*)"no_ter_antartique=",no_ter_antartique
  IF (no_ter_antartique) THEN
     ! Remove all land fractions souther than 60S and set land-ice instead
     WRITE(lunout,*) "Remove land fractions souther than 60deg south by increasing"
     WRITE(lunout,*) "the continental ice fractions. No land can now be found at Antartic."
     DO ji=1, klon
        IF (latitude_deg(ji)<-60.0) THEN
           pctsrf(ji,is_lic) = pctsrf(ji,is_lic) + pctsrf(ji,is_ter)
           pctsrf(ji,is_ter) = 0
        END IF
     END DO
  END IF
    
! sub-surface ocean and sea ice (sea ice set to zero for start)
!*******************************************************************************
    pctsrf(:,is_oce)=(1.-zmasq(:))
    WHERE(pctsrf(:,is_oce)<EPSFRA) pctsrf(:,is_oce)=0.
    
    zval(:)=max(0.,zmea-2*zstd(:))
    zpic(:)=zmea+2*zstd(:)
    
!! WARNING    DON'T FORGET FOR LATER
!!ym  IF(couple) pctsrf(:,is_oce)=ocemask_fi(:)
!! 
    
! Init: tsol, qsol, sn, evap, tsoil, rain_fall, snow_fall, solsw, sollw, frugs
!*******************************************************************************
    DO i=1,nbsrf
     ftsol(:,i) = tsol
    END DO
  
    DO i=1,nbsrf
     snsrf(:,i) = sn
    END DO
!albedo SB >>>
!ym error : the sub surface dimension is the third not second
!    falb_dir(:,is_ter,:)=0.08; falb_dir(:,is_lic,:)=0.6
!    falb_dir(:,is_oce,:)=0.5;  falb_dir(:,is_sic,:)=0.6
    falb_dir(:,:,is_ter)=0.08; falb_dir(:,:,is_lic)=0.6
    falb_dir(:,:,is_oce)=0.5;  falb_dir(:,:,is_sic)=0.6

!ym falb_dif has been forgotten, initialize with defaukt value found in phyetat0 or 0 ?
!ym probably the uninitialized value was 0 for standard (regular grid) case
    falb_dif(:,:,:)=0

!albedo SB <<<
    fevap(:,:) = 0.
    DO i=1,nbsrf
     qsolsrf(:,i)=150.
    END DO
 
    DO i=1,nbsrf 
      DO j=1,nsoilmx
        tsoil(:,j,i) = tsol
      END DO
    END DO
 
    rain_fall = 0.; snow_fall = 0.
    solsw = 165.;   sollw = -53.
!ym warning missing init for sollwdown => set to 0
  sollwdown  = 0.
    
    
    t_ancien = 273.15
    u_ancien=0
    v_ancien=0
    q_ancien = 0.
    agesno = 0.

    z0m(:,is_oce) = rugmer(:)

   z0m(:,is_ter) = MAX(1.0e-05,zstd(:)*zsig(:)/2.0)
   z0m(:,is_lic) = MAX(1.0e-05,zstd(:)*zsig(:)/2.0)

   z0m(:,is_sic) = 0.001
   z0h(:,:)=z0m(:,:)

    fder = 0.0
    clwcon = 0.0
    rnebcon = 0.0
    ratqs = 0.0
    run_off_lic_0 = 0.0 
    rugoro = 0.0

! Before phyredem calling, surface modules and values to be saved in startphy.nc
! are initialized
!*******************************************************************************
    pbl_tke(:,:,:) = 1.e-8 
    zmax0(:) = 40.
    f0(:) = 1.e-5
    sig1(:,:) = 0.
    w01(:,:) = 0.
    wake_deltat(:,:) = 0.
    wake_deltaq(:,:) = 0.
    wake_s(:) = 0.
    wake_cstar(:) = 0.
    wake_fip(:) = 0.
    wake_pe = 0.
    fm_therm = 0.
    entr_therm = 0.
    detr_therm = 0.
    ale_bl = 0.
    ale_bl_trig =0. 
    alp_bl =0.
    CALL fonte_neige_init(run_off_lic_0)
    CALL pbl_surface_init( fder, snsrf, qsolsrf, tsoil )

    CALL gather_omp(cell_area,cell_area_mpi)
    CALL gather_omp(pctsrf,pctsrf_mpi)
    IF (is_omp_master) THEN
      CALL xios_send_field("area_ce0l",cell_area_mpi)
      CALL xios_send_field("fract_oce_ce0l",pctsrf_mpi(:,is_oce))
      CALL xios_send_field("fract_sic_ce0l",pctsrf_mpi(:,is_sic))
    ENDIF
    
    CALL phyredem( "startphy.nc" )

#endif
  END SUBROUTINE create_etat0_unstruct


END MODULE create_etat0_unstruct_mod
