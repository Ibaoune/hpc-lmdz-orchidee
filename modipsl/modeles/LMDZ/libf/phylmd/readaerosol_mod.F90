! $Id: readaerosol.F90 3436 2019-01-22 16:26:21Z emillour $
!
MODULE readaerosol_mod

  REAL, SAVE :: not_valid=-333.
  
  INTEGER, SAVE :: nbp_lon_src
!$OMP THREADPRIVATE(nbp_lon_src)  
  INTEGER, SAVE :: nbp_lat_src
!$OMP THREADPRIVATE(nbp_lat_src)  
  REAL, ALLOCATABLE, SAVE    :: psurf_interp(:,:)
!$OMP THREADPRIVATE(psurf_interp)  

CONTAINS

SUBROUTINE readaerosol(name_aero, type, filename, iyr_in, klev_src, pt_ap, pt_b, pt_out, psurf, load)

!****************************************************************************************
! This routine will read the aersosol from file. 
!
! Read a year data with get_aero_fromfile depending on aer_type : 
! - actuel   : read year 1980
! - preind   : read natural data
! - scenario : read one or two years and do eventually linare time interpolation
!
! Return pointer, pt_out, to the year read or result from interpolation
!****************************************************************************************
  USE dimphy
  USE print_control_mod, ONLY: lunout

  IMPLICIT NONE

  ! Input arguments
  CHARACTER(len=7), INTENT(IN) :: name_aero
  CHARACTER(len=*), INTENT(IN) :: type  ! actuel, annuel, scenario or preind
  CHARACTER(len=8), INTENT(IN) :: filename
  INTEGER, INTENT(IN)          :: iyr_in

  ! Output
  INTEGER, INTENT(OUT)            :: klev_src
  REAL, POINTER, DIMENSION(:)     :: pt_ap        ! Pointer for describing the vertical levels      
  REAL, POINTER, DIMENSION(:)     :: pt_b         ! Pointer for describing the vertical levels      
  REAL, POINTER, DIMENSION(:,:,:) :: pt_out       ! The massvar distributions, DIMENSION(klon, klev_src, 12)
  REAL, DIMENSION(klon,12), INTENT(OUT) :: psurf  ! Surface pression for 12 months
  REAL, DIMENSION(klon,12), INTENT(OUT) :: load   ! Aerosol mass load in each column for 12 months

  ! Local variables
  CHARACTER(len=4)                :: cyear
  REAL, POINTER, DIMENSION(:,:,:) :: pt_2
  REAL, DIMENSION(klon,12)        :: psurf2, load2
  INTEGER                         :: iyr1, iyr2, klev_src2
  INTEGER                         :: it, k, i
  LOGICAL, PARAMETER              :: lonlyone=.FALSE.

!****************************************************************************************
! Read data depending on aer_type
!
!****************************************************************************************

  IF (type == 'actuel') THEN
! Read and return data for year 1980
!****************************************************************************************
     cyear='1980'
     ! get_aero_fromfile returns pt_out allocated and initialized with data for 12 month 
     ! pt_out has dimensions (klon, klev_src, 12)
     CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, pt_out, psurf, load)
     

  ELSE IF (type == 'preind') THEN
! Read and return data from file with suffix .nat
!****************************************************************************************     
     cyear='.nat'
     ! get_aero_fromfile returns pt_out allocated and initialized with data for 12 month 
     ! pt_out has dimensions (klon, klev_src, 12)
     CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, pt_out, psurf, load)
     
  ELSE IF (type == 'annuel') THEN
! Read and return data from scenario annual files
!****************************************************************************************     
     WRITE(cyear,'(I4)') iyr_in
     WRITE(lunout,*) 'get_aero 3 iyr_in=', iyr_in,'   ',cyear
     ! get_aero_fromfile returns pt_out allocated and initialized with data for nbr_tsteps month 
     ! pt_out has dimensions (klon, klev_src, 12)
     CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, pt_out, psurf, load)
     
  ELSE IF (type == 'scenario') THEN
! Read data depending on actual year and interpolate if necessary
!****************************************************************************************
     IF (iyr_in .LT. 1850) THEN
        cyear='.nat'
        WRITE(lunout,*) 'get_aero 1 iyr_in=', iyr_in,'   ',cyear
        ! get_aero_fromfile returns pt_out allocated and initialized with data for 12 month 
        ! pt_out has dimensions (klon, klev_src, 12)
        CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, pt_out, psurf, load)
        
     ELSE IF (iyr_in .GE. 2100) THEN
        cyear='2100'
        WRITE(lunout,*) 'get_aero 2 iyr_in=', iyr_in,'   ',cyear
        ! get_aero_fromfile returns pt_out allocated and initialized with data for 12 month 
        ! pt_out has dimensions (klon, klev_src, 12)
        CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, pt_out, psurf, load)
        
     ELSE
        ! Read data from 2 decades and interpolate to actual year
        ! a) from actual 10-yr-period
        IF (iyr_in.LT.1900) THEN
           iyr1 = 1850
           iyr2 = 1900
        ELSE IF (iyr_in.GE.1900.AND.iyr_in.LT.1920) THEN
           iyr1 = 1900
           iyr2 = 1920
        ELSE 
           iyr1 = INT(iyr_in/10)*10
           iyr2 = INT(1+iyr_in/10)*10
        ENDIF
        
        WRITE(cyear,'(I4)') iyr1
        WRITE(lunout,*) 'get_aero 3 iyr_in=', iyr_in,'   ',cyear
        ! get_aero_fromfile returns pt_out allocated and initialized with data for 12 month 
        ! pt_out has dimensions (klon, klev_src, 12)
        CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, pt_out, psurf, load)
        
        ! If to read two decades:
        IF (.NOT.lonlyone) THEN 
           
           ! b) from the next following one
           WRITE(cyear,'(I4)') iyr2
           WRITE(lunout,*) 'get_aero 4 iyr_in=', iyr_in,'   ',cyear
           
           NULLIFY(pt_2)
           ! get_aero_fromfile returns pt_2 allocated and initialized with data for 12 month 
           ! pt_2 has dimensions (klon, klev_src, 12)
           CALL get_aero_fromfile(name_aero, cyear, filename, klev_src2, pt_ap, pt_b, pt_2, psurf2, load2)
           ! Test for same number of vertical levels
           IF (klev_src /= klev_src2) THEN
              WRITE(lunout,*) 'Two aerosols files with different number of vertical levels is not allowded'
              CALL abort_physic('readaersosol','Error in number of vertical levels',1)
           END IF
           
           ! Linare interpolate to the actual year:
           DO it=1,12
              DO k=1,klev_src
                 DO i = 1, klon
                    pt_out(i,k,it) = &
                         pt_out(i,k,it) - REAL(iyr_in-iyr1)/REAL(iyr2-iyr1) * &
                         (pt_out(i,k,it) - pt_2(i,k,it))
                 END DO
              END DO

              DO i = 1, klon
                 psurf(i,it) = &
                      psurf(i,it) - REAL(iyr_in-iyr1)/REAL(iyr2-iyr1) * &
                      (psurf(i,it) - psurf2(i,it))

                 load(i,it) = &
                      load(i,it) - REAL(iyr_in-iyr1)/REAL(iyr2-iyr1) * &
                      (load(i,it) - load2(i,it))
              END DO
           END DO

           ! Deallocate pt_2 no more needed
           DEALLOCATE(pt_2)
           
        END IF ! lonlyone
     END IF ! iyr_in .LT. 1850

  ELSE
     WRITE(lunout,*)'This option is not implemented : aer_type = ', type,' name_aero=',name_aero
     CALL abort_physic('readaerosol','Error : aer_type parameter not accepted',1)
  END IF ! type


END SUBROUTINE readaerosol


  SUBROUTINE init_aero_fromfile(flag_aerosol)
  USE netcdf
  USE mod_phys_lmdz_para
  USE mod_grid_phy_lmdz, ONLY: grid_type, unstructured
#ifdef CPP_XIOS
  USE xios
#endif
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: flag_aerosol
#ifdef CPP_XIOS
  REAL,ALLOCATABLE :: lat_src(:)
  REAL,ALLOCATABLE :: lon_src(:)
  CHARACTER(LEN=*),PARAMETER :: file_aerosol='aerosols.nat.nc'
  CHARACTER(LEN=*),PARAMETER :: file_so4='so4.nat.nc'
  INTEGER :: klev_src
  INTEGER :: ierr ,ncid, dimID, varid
  REAL :: null_array(0)

  IF (flag_aerosol>0 .AND. grid_type==unstructured) THEN
  
    IF (is_omp_root) THEN
  
      IF (is_mpi_root) THEN
    
        IF (nf90_open(TRIM(file_aerosol), NF90_NOWRITE, ncid) /= NF90_NOERR) THEN
          CALL check_err( nf90_open(TRIM(file_so4), NF90_NOWRITE, ncid), "pb open "//trim(file_so4) )
        ENDIF

        ! Read and test longitudes
        CALL check_err( nf90_inq_dimid(ncid, "lon", dimID),"pb inq dim lon") 
        CALL check_err( nf90_inquire_dimension(ncid, dimID, len = nbp_lon_src),"pb inq dim lon")
        CALL check_err( nf90_inq_varid(ncid, 'lon', varid),"pb inq lon" )
        ALLOCATE(lon_src(nbp_lon_src))
        CALL check_err( nf90_get_var(ncid, varid, lon_src(:)),"pb get lon" )

        ! Read and test latitudes
        CALL check_err( nf90_inq_dimid(ncid, "lat", dimID),"pb inq dim lat") 
        CALL check_err( nf90_inquire_dimension(ncid, dimID, len = nbp_lat_src),"pb inq dim lat")
        CALL check_err( nf90_inq_varid(ncid, 'lat', varid),"pb inq lat" )
        ALLOCATE(lat_src(nbp_lat_src))
        CALL check_err( nf90_get_var(ncid, varid, lat_src(:)),"pb get lat" )
        IF (nf90_inq_dimid(ncid, 'lev', dimid) /= NF90_NOERR) THEN
          IF (nf90_inq_dimid(ncid, 'presnivs', dimid)/= NF90_NOERR) THEN
             CALL check_err(nf90_inq_dimid(ncid, 'PRESNIVS', dimid),'dimension lev,PRESNIVS or presnivs not in file')
          ENDIF
        ENDIF
        CALL check_err( nf90_inquire_dimension(ncid, dimid, len = klev_src),"pb inq dim for PRESNIVS or lev" )
        CALL check_err( nf90_close(ncid),"pb in close" )   
      ENDIF

      CALL bcast_mpi(nbp_lat_src)
      CALL bcast_mpi(nbp_lon_src)
      CALL bcast_mpi(klev_src)

      IF (is_mpi_root ) THEN
        CALL xios_set_domain_attr("domain_aerosol",nj_glo=nbp_lat_src, nj=nbp_lat_src, jbegin=0, latvalue_1d=lat_src)
        CALL xios_set_domain_attr("domain_aerosol",ni_glo=nbp_lon_src, ni=nbp_lon_src, ibegin=0, lonvalue_1d=lon_src)
      ELSE
        CALL xios_set_domain_attr("domain_aerosol",nj_glo=nbp_lat_src, nj=0, jbegin=0, latvalue_1d=null_array )
        CALL xios_set_domain_attr("domain_aerosol",ni_glo=nbp_lon_src, ni=0, ibegin=0, lonvalue_1d=null_array)
      ENDIF
      CALL xios_set_axis_attr("axis_aerosol",n_glo=klev_src)
      CALL xios_set_fieldgroup_attr("aerosols", enabled=.TRUE.)
  
    ENDIF
    
  ENDIF    
#endif 
  END SUBROUTINE init_aero_fromfile


    

  SUBROUTINE get_aero_fromfile(varname, cyr, filename, klev_src, pt_ap, pt_b, pt_year, psurf_out, load_out)
!****************************************************************************************
! Read 12 month aerosol from file and distribute to local process on physical grid. 
! Vertical levels, klev_src, may differ from model levels if new file format.
!
! For mpi_root and master thread :
! 1) Open file 
! 2) Find vertical dimension klev_src
! 3) Read field month by month
! 4) Close file  
! 5) Transform the global field from 2D(nbp_lon,nbp_lat) to 1D(klon_glo)
!     - Also the levels and the latitudes have to be inversed
!
! For all processes and threads :
! 6) Scatter global field(klon_glo) to local process domain(klon)
! 7) Test for negative values
!****************************************************************************************

    USE netcdf
    USE dimphy
    USE mod_grid_phy_lmdz, ONLY: nbp_lon_=>nbp_lon, nbp_lat_=>nbp_lat, klon_glo, &
                                 grid2Dto1D_glo, grid_type, unstructured
    USE mod_phys_lmdz_para
    USE iophy, ONLY : io_lon, io_lat
    USE print_control_mod, ONLY: lunout
#ifdef CPP_XIOS
    USE xios
#endif
    IMPLICIT NONE
      
! Input argumets
    CHARACTER(len=7), INTENT(IN)          :: varname
    CHARACTER(len=4), INTENT(IN)          :: cyr
    CHARACTER(len=8), INTENT(IN)          :: filename

! Output arguments
    INTEGER, INTENT(OUT)                  :: klev_src     ! Number of vertical levels in file
    REAL, POINTER, DIMENSION(:)           :: pt_ap        ! Pointer for describing the vertical levels      
    REAL, POINTER, DIMENSION(:)           :: pt_b         ! Pointer for describing the vertical levels      
    REAL, POINTER, DIMENSION(:,:,:)       :: pt_year      ! Pointer-variabale from file, 12 month, grid : klon,klev_src
    REAL, POINTER, DIMENSION(:,:,:)       :: pt_year_mpi  ! Pointer-variabale from file, 12 month, grid : klon,klev_src
    REAL, DIMENSION(klon,12), INTENT(OUT) :: psurf_out    ! Surface pression for 12 months
    REAL, DIMENSION(klon_mpi,12)          :: psurf_out_mpi    ! Surface pression for 12 months
    REAL, DIMENSION(klon,12), INTENT(OUT) :: load_out     ! Aerosol mass load in each column
    REAL, DIMENSION(klon_mpi,12)          :: load_out_mpi     ! Aerosol mass load in each column
    INTEGER                               :: nbr_tsteps   ! number of month in file read

! Local variables
    CHARACTER(len=30)     :: fname
    CHARACTER(len=30)     :: cvar
    INTEGER               :: ncid, dimid, varid
    INTEGER               :: imth, i, j, k, ierr
    REAL                  :: npole, spole
    REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: varmth
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: varyear       ! Global variable read from file, 12 month
    REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: varyear_glo1D !(klon_glo, klev_src, 12)
    REAL, ALLOCATABLE, DIMENSION(:)       :: varktmp

    REAL,  ALLOCATABLE                    :: psurf_glo2D(:,:,:)   ! Surface pression for 12 months on dynamics global grid
    REAL, DIMENSION(klon_glo,12)          :: psurf_glo1D   ! -"- on physical global grid
    REAL,  ALLOCATABLE                    :: load_glo2D(:,:,:)    ! Load for 12 months on dynamics global grid
    REAL, DIMENSION(klon_glo,12)          :: load_glo1D    ! -"- on physical global grid
    REAL, ALLOCATABLE, DIMENSION(:,:)     :: vartmp
    REAL, ALLOCATABLE,DIMENSION(:)        :: lon_src              ! longitudes in file
    REAL, ALLOCATABLE,DIMENSION(:)        :: lat_src, lat_src_inv ! latitudes in file
    LOGICAL                               :: new_file             ! true if new file format detected
    LOGICAL                               :: invert_lat           ! true if the field has to be inverted for latitudes
    INTEGER                               :: nbp_lon, nbp_lat
    LOGICAL,SAVE                          :: first=.TRUE.
!$OMP THREADPRIVATE(first)
    
    IF (grid_type==unstructured) THEN
      nbp_lon=nbp_lon_src
      nbp_lat=nbp_lat_src
    ELSE
      nbp_lon=nbp_lon_
      nbp_lat=nbp_lat_
    ENDIF
    
    IF (is_mpi_root) THEN
    
      ALLOCATE(psurf_glo2D(nbp_lon,nbp_lat,12))
      ALLOCATE(load_glo2D(nbp_lon,nbp_lat,12))
      ALLOCATE(vartmp(nbp_lon,nbp_lat))
      ALLOCATE(lon_src(nbp_lon))
      ALLOCATE(lat_src(nbp_lat))
      ALLOCATE(lat_src_inv(nbp_lat))
    ELSE
      ALLOCATE(varyear(0,0,0,0))
      ALLOCATE(psurf_glo2D(0,0,0))
      ALLOCATE(load_glo2D(0,0,0))
    ENDIF
            
    ! Deallocate pointers
    IF (ASSOCIATED(pt_ap)) DEALLOCATE(pt_ap)
    IF (ASSOCIATED(pt_b))  DEALLOCATE(pt_b)

    IF (is_mpi_root .AND. is_omp_root) THEN

! 1) Open file 
!****************************************************************************************
! Add suffix to filename
       fname = trim(filename)//cyr//'.nc'
  
       WRITE(lunout,*) 'reading variable ',TRIM(varname),' in file ', TRIM(fname)
       CALL check_err( nf90_open(TRIM(fname), NF90_NOWRITE, ncid), "pb open "//trim(fname) )


       IF (grid_type/=unstructured) THEN

! Test for equal longitudes and latitudes in file and model
!****************************************************************************************
       ! Read and test longitudes
         CALL check_err( nf90_inq_varid(ncid, 'lon', varid),"pb inq lon" )
         CALL check_err( nf90_get_var(ncid, varid, lon_src(:)),"pb get lon" )
       
         IF (maxval(ABS(lon_src - io_lon)) > 0.001) THEN
            WRITE(lunout,*) 'Problem in longitudes read from file : ',TRIM(fname)
            WRITE(lunout,*) 'longitudes in file ', TRIM(fname),' : ', lon_src
            WRITE(lunout,*) 'longitudes in model :', io_lon
          
            CALL abort_physic('get_aero_fromfile', 'longitudes are not the same in file and model',1)
         END IF

         ! Read and test latitudes
         CALL check_err( nf90_inq_varid(ncid, 'lat', varid),"pb inq lat" )
         CALL check_err( nf90_get_var(ncid, varid, lat_src(:)),"pb get lat" )

         ! Invert source latitudes
         DO j = 1, nbp_lat
            lat_src_inv(j) = lat_src(nbp_lat +1 -j)
         END DO

         IF (maxval(ABS(lat_src - io_lat)) < 0.001) THEN
            ! Latitudes are the same
            invert_lat=.FALSE.
         ELSE IF (maxval(ABS(lat_src_inv - io_lat)) < 0.001) THEN
            ! Inverted source latitudes correspond to model latitudes
            WRITE(lunout,*) 'latitudes will be inverted for file : ',TRIM(fname)
            invert_lat=.TRUE.
         ELSE
            WRITE(lunout,*) 'Problem in latitudes read from file : ',TRIM(fname)
            WRITE(lunout,*) 'latitudes in file ', TRIM(fname),' : ', lat_src      
            WRITE(lunout,*) 'latitudes in model :', io_lat
            CALL abort_physic('get_aero_fromfile', 'latitudes do not correspond between file and model',1)
         END IF
       ENDIF

! 2) Check if old or new file is avalabale.
!    New type of file should contain the dimension 'lev'
!    Old type of file should contain the dimension 'PRESNIVS'
!****************************************************************************************
       ierr = nf90_inq_dimid(ncid, 'lev', dimid) 
       IF (ierr /= NF90_NOERR) THEN
          ! Coordinate axe lev not found. Check for presnivs.
          ierr = nf90_inq_dimid(ncid, 'presnivs', dimid) 
          IF (ierr /= NF90_NOERR) THEN
             ierr = nf90_inq_dimid(ncid, 'PRESNIVS', dimid)
             IF (ierr /= NF90_NOERR) THEN
                ! Dimension PRESNIVS not found either
                CALL abort_physic('get_aero_fromfile', 'dimension lev,PRESNIVS or presnivs not in file',1)
             ELSE 
                ! Old file found
                new_file=.FALSE.
                WRITE(lunout,*) 'Vertical interpolation for ',TRIM(varname),' will not be done'
             END IF
          ELSE
             ! New file found
             new_file=.TRUE.
             WRITE(lunout,*) 'Vertical interpolation for ',TRIM(varname),' will be done'
          ENDIF
       ELSE
          ! New file found
          new_file=.TRUE.
          WRITE(lunout,*) 'Vertical interpolation for ',TRIM(varname),' will be done'
       END IF
       
! 2) Find vertical dimension klev_src
!****************************************************************************************
       CALL check_err( nf90_inquire_dimension(ncid, dimid, len = klev_src),"pb inq dim for PRESNIVS or lev" )
       
     ! Allocate variables depending on the number of vertical levels
       ALLOCATE(varmth(nbp_lon,nbp_lat, klev_src), varyear(nbp_lon,nbp_lat, klev_src, 12), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('get_aero_fromfile', 'pb in allocation 1',1)

       ALLOCATE(pt_ap(klev_src), pt_b(klev_src), varktmp(klev_src), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('get_aero_fromfile', 'pb in allocation 2',1)

! 3) Read all variables from file
!    There is 2 options for the file structure :
!    new_file=TRUE  : read varyear, ps, pt_ap and pt_b
!    new_file=FALSE : read varyear month by month
!****************************************************************************************

       IF (new_file) THEN
! ++) Check number of month in file opened
!**************************************************************************************************
       ierr = nf90_inq_dimid(ncid, 'TIME',dimid)
       if (ierr /= NF90_NOERR) THEN
          ierr = nf90_inq_dimid(ncid, 'time_counter', dimid)
       ENDIF
       CALL check_err( nf90_inquire_dimension(ncid, dimid, len = nbr_tsteps),"pb inq dim TIME or time_counter" )
!       IF (nbr_tsteps /= 12 .AND. nbr_tsteps /= 14) THEN
       IF (nbr_tsteps /= 12 ) THEN
    CALL abort_physic('get_aero_fromfile', 'not the right number of months in aerosol file read (should be 12 for the moment)' &
                     ,1)
       ENDIF

! ++) Read the aerosol concentration month by month and concatenate to total variable varyear
!****************************************************************************************
          ! Get variable id
          !CALL check_err( nf90_inq_varid(ncid, TRIM(varname), varid),"pb inq var "//TRIM(varname) )
          print *,'readaerosol ', TRIM(varname)
          IF ( nf90_inq_varid(ncid, TRIM(varname), varid) /= NF90_NOERR ) THEN 
            ! Variable is not there
            WRITE(lunout,*) 'Attention '//TRIM(varname)//' is not in aerosol input file'
            varyear(:,:,:,:)=0.0
          ELSE 
            ! Get the variable
            CALL check_err( nf90_get_var(ncid, varid, varyear(:,:,:,:)),"pb get var "//TRIM(varname) )
          ENDIF
          
! ++) Read surface pression, 12 month in one variable
!****************************************************************************************
          ! Get variable id
          CALL check_err( nf90_inq_varid(ncid, "ps", varid),"pb inq var ps" )
          ! Get the variable
          CALL check_err( nf90_get_var(ncid, varid, psurf_glo2D),"pb get var ps" )
          
! ++) Read mass load, 12 month in one variable
!****************************************************************************************
          ! Get variable id
          !CALL check_err( nf90_inq_varid(ncid, "load_"//TRIM(varname), varid) ,"pb inq var load_"//TRIM(varname))
          IF ( nf90_inq_varid(ncid, "load_"//TRIM(varname), varid) /= NF90_NOERR ) THEN 
            WRITE(lunout,*) 'Attention load_'//TRIM(varname)//' is not in aerosol input file'
            load_glo2D(:,:,:)=0.0
          ELSE
            ! Get the variable
            CALL check_err( nf90_get_var(ncid, varid, load_glo2D),"pb get var load_"//TRIM(varname) )
          ENDIF
          
! ++) Read ap
!****************************************************************************************
          ! Get variable id
          CALL check_err( nf90_inq_varid(ncid, "ap", varid),"pb inq var ap" )
          ! Get the variable
          CALL check_err( nf90_get_var(ncid, varid, pt_ap),"pb get var ap" )

! ++) Read b
!****************************************************************************************
          ! Get variable id
          CALL check_err( nf90_inq_varid(ncid, "b", varid),"pb inq var b" )
          ! Get the variable
          CALL check_err( nf90_get_var(ncid, varid, pt_b),"pb get var b" )

          

       ELSE  ! old file

! ++) Read the aerosol concentration month by month and concatenate to total variable varyear
!****************************************************************************************
          DO imth=1, 12
             IF (imth.EQ.1) THEN
                cvar=TRIM(varname)//'JAN'
             ELSE IF (imth.EQ.2) THEN
                cvar=TRIM(varname)//'FEB'
             ELSE IF (imth.EQ.3) THEN
                cvar=TRIM(varname)//'MAR'
             ELSE IF (imth.EQ.4) THEN
                cvar=TRIM(varname)//'APR'
             ELSE IF (imth.EQ.5) THEN
                cvar=TRIM(varname)//'MAY'
             ELSE IF (imth.EQ.6) THEN
                cvar=TRIM(varname)//'JUN'
             ELSE IF (imth.EQ.7) THEN
                cvar=TRIM(varname)//'JUL'
             ELSE IF (imth.EQ.8) THEN
                cvar=TRIM(varname)//'AUG'
             ELSE IF (imth.EQ.9) THEN
                cvar=TRIM(varname)//'SEP'
             ELSE IF (imth.EQ.10) THEN
                cvar=TRIM(varname)//'OCT'
             ELSE IF (imth.EQ.11) THEN
                cvar=TRIM(varname)//'NOV'
             ELSE IF (imth.EQ.12) THEN
                cvar=TRIM(varname)//'DEC'
             END IF
             
             ! Get variable id
             CALL check_err( nf90_inq_varid(ncid, TRIM(cvar), varid),"pb inq var "//TRIM(cvar) )
             
             ! Get the variable
             CALL check_err( nf90_get_var(ncid, varid, varmth),"pb get var "//TRIM(cvar) )
             
             ! Store in variable for the whole year
             varyear(:,:,:,imth)=varmth(:,:,:)
             
          END DO
          
          ! Putting dummy 
          psurf_glo2D(:,:,:) = not_valid
          load_glo2D(:,:,:)  = not_valid
          pt_ap(:) = not_valid
          pt_b(:)  = not_valid

       END IF

! 4) Close file  
!****************************************************************************************
       CALL check_err( nf90_close(ncid),"pb in close" )
     

! 5) Transform the global field from 2D(nbp_lon,nbp_lat) to 1D(klon_glo)
!****************************************************************************************
! Test if vertical levels have to be inversed

       IF ((pt_b(1) < pt_b(klev_src)) .OR. .NOT. new_file) THEN
!          WRITE(lunout,*) 'Vertical axis in file ',TRIM(fname), ' needs to be inverted'
!          WRITE(lunout,*) 'before pt_ap = ', pt_ap
!          WRITE(lunout,*) 'before pt_b = ', pt_b
          
          ! Inverse vertical levels for varyear 
          DO imth=1, 12
             varmth(:,:,:) = varyear(:,:,:,imth) ! use varmth temporarly
             DO k=1, klev_src
                DO j=1, nbp_lat
                   DO i=1, nbp_lon
                      varyear(i,j,k,imth) = varmth(i,j,klev_src+1-k)
                   END DO
                END DO
             END DO
          END DO
           
          ! Inverte vertical axes for pt_ap and pt_b
          varktmp(:) = pt_ap(:)
          DO k=1, klev_src
             pt_ap(k) = varktmp(klev_src+1-k)
          END DO

          varktmp(:) = pt_b(:)
          DO k=1, klev_src
             pt_b(k) = varktmp(klev_src+1-k)
          END DO
          WRITE(lunout,*) 'after pt_ap = ', pt_ap
          WRITE(lunout,*) 'after pt_b = ', pt_b

       ELSE 
          WRITE(lunout,*) 'Vertical axis in file ',TRIM(fname), ' is ok, no vertical inversion is done'
          WRITE(lunout,*) 'pt_ap = ', pt_ap
          WRITE(lunout,*) 'pt_b = ', pt_b
       END IF


       IF (grid_type/=unstructured) THEN
  !     - Invert latitudes if necessary
         DO imth=1, 12
            IF (invert_lat) THEN

               ! Invert latitudes for the variable
               varmth(:,:,:) = varyear(:,:,:,imth) ! use varmth temporarly
               DO k=1,klev_src
                  DO j=1,nbp_lat
                     DO i=1,nbp_lon
                        varyear(i,j,k,imth) = varmth(i,nbp_lat+1-j,k)
                     END DO
                  END DO
               END DO
             
               ! Invert latitudes for surface pressure
               vartmp(:,:) = psurf_glo2D(:,:,imth)
               DO j=1,nbp_lat
                  DO i=1,nbp_lon
                     psurf_glo2D(i,j,imth)= vartmp(i,nbp_lat+1-j)
                  END DO
               END DO
             
               ! Invert latitudes for the load
               vartmp(:,:) = load_glo2D(:,:,imth)
               DO j=1,nbp_lat
                  DO i=1,nbp_lon
                     load_glo2D(i,j,imth)= vartmp(i,nbp_lat+1-j)
                  END DO
               END DO
            END IF ! invert_lat
             
            ! Do zonal mead at poles and distribut at whole first and last latitude
            DO k=1, klev_src
               npole=0.  ! North pole, j=1
               spole=0.  ! South pole, j=nbp_lat        
               DO i=1,nbp_lon
                  npole = npole + varyear(i,1,k,imth)
                  spole = spole + varyear(i,nbp_lat,k,imth)
               END DO
               npole = npole/REAL(nbp_lon)
               spole = spole/REAL(nbp_lon)
               varyear(:,1,    k,imth) = npole
               varyear(:,nbp_lat,k,imth) = spole
            END DO
         END DO ! imth
       
         ALLOCATE(varyear_glo1D(klon_glo, klev_src, 12), stat=ierr)
         IF (ierr /= 0) CALL abort_physic('get_aero_fromfile', 'pb in allocation 3',1)
       
         ! Transform from 2D to 1D field
         CALL grid2Dto1D_glo(varyear,varyear_glo1D)
         CALL grid2Dto1D_glo(psurf_glo2D,psurf_glo1D)
         CALL grid2Dto1D_glo(load_glo2D,load_glo1D)
      
      ENDIF
      
    ELSE
        ALLOCATE(varyear_glo1D(0,0,0))        
    END IF ! is_mpi_root .AND. is_omp_root

!$OMP BARRIER
  
! 6) Distribute to all processes
!    Scatter global field(klon_glo) to local process domain(klon)
!    and distribute klev_src to all processes
!****************************************************************************************

    ! Distribute klev_src
    CALL bcast(klev_src)

    ! Allocate and distribute pt_ap and pt_b
    IF (.NOT. ASSOCIATED(pt_ap)) THEN  ! if pt_ap is allocated also pt_b is allocated
       ALLOCATE(pt_ap(klev_src), pt_b(klev_src), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('get_aero_fromfile', 'pb in allocation 4',1)
    END IF
    CALL bcast(pt_ap)
    CALL bcast(pt_b)

    ! Allocate space for output pointer variable at local process
    IF (ASSOCIATED(pt_year)) DEALLOCATE(pt_year)
    ALLOCATE(pt_year(klon, klev_src, 12), stat=ierr)
    ALLOCATE(pt_year_mpi(klon_mpi, klev_src, 12), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('get_aero_fromfile', 'pb in allocation 5',1)

    IF (grid_type==unstructured) THEN
#ifdef CPP_XIOS
      IF (is_omp_master) THEN
        CALL xios_send_field(TRIM(varname)//"_in",varyear)
        CALL xios_recv_field(TRIM(varname)//"_out",pt_year_mpi)
        CALL xios_send_field("load_"//TRIM(varname)//"_in",load_glo2D)
        CALL xios_recv_field("load_"//TRIM(varname)//"_out",load_out_mpi)
        IF (first) THEN
          ALLOCATE(psurf_interp(klon_mpi,12))
          CALL xios_send_field("psurf_aerosol_in",psurf_glo2D)
          CALL xios_recv_field("psurf_aerosol_out",psurf_interp)
        ENDIF
      ENDIF
      CALL scatter_omp(pt_year_mpi,pt_year)
      CALL scatter_omp(load_out_mpi,load_out)
      CALL scatter_omp(psurf_interp,psurf_out)
      first=.FALSE.
#endif
    ELSE
      ! Scatter global field to local domain at local process
      CALL scatter(varyear_glo1D, pt_year)
      CALL scatter(psurf_glo1D, psurf_out)
      CALL scatter(load_glo1D,  load_out)
    ENDIF
! 7) Test for negative values
!****************************************************************************************
    IF (MINVAL(pt_year) < 0.) THEN
       WRITE(lunout,*) 'Warning! Negative values read from file :', fname
    END IF

  END SUBROUTINE get_aero_fromfile


  SUBROUTINE check_err(status,text)
    USE netcdf
    USE print_control_mod, ONLY: lunout
    IMPLICIT NONE

    INTEGER, INTENT (IN) :: status
    CHARACTER(len=*), INTENT (IN), OPTIONAL :: text

    IF (status /= NF90_NOERR) THEN
       WRITE(lunout,*) 'Error in get_aero_fromfile, netcdf error code = ',status
       IF (PRESENT(text)) THEN
          WRITE(lunout,*) 'Error in get_aero_fromfile : ',text
       END IF
       CALL abort_physic('get_aero_fromfile',trim(nf90_strerror(status)),1)
    END IF

  END SUBROUTINE check_err


END MODULE readaerosol_mod
