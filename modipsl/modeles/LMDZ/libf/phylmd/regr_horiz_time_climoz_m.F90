MODULE regr_horiz_time_climoz_m

  USE interpolation,     ONLY: locate
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, grid_type, unstructured
  USE nrtype,            ONLY: pi
  USE netcdf,   ONLY: NF90_CLOBBER, NF90_FLOAT,     NF90_GET_VAR, NF90_OPEN,   &
                      NF90_NOWRITE, NF90_NOERR,     NF90_GET_ATT, NF90_GLOBAL
  USE netcdf95, ONLY: NF95_DEF_DIM, NF95_INQ_DIMID, NF95_INQUIRE_DIMENSION,    &
                      NF95_DEF_VAR, NF95_INQ_VARID, NF95_INQUIRE_VARIABLE,     &
          NF95_OPEN,  NF95_CREATE,  NF95_GET_ATT,   NF95_GW_VAR,  HANDLE_ERR,  &
          NF95_CLOSE, NF95_ENDDEF,  NF95_PUT_ATT,   NF95_PUT_VAR, NF95_COPY_ATT
  USE print_control_mod, ONLY: lunout
  USE dimphy
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: regr_horiz_time_climoz
  REAL, PARAMETER :: deg2rad=pi/180.
  CHARACTER(LEN=13), PARAMETER :: vars_in(2)=['tro3         ','tro3_daylight']

  INTEGER :: nlat_ou, nlon_ou
  REAL, ALLOCATABLE :: latitude_glo(:)
!$OMP THREADPRIVATE(latitude_glo)
  INTEGER, ALLOCATABLE :: ind_cell_glo_glo(:)
!$OMP THREADPRIVATE(ind_cell_glo_glo)

CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE regr_horiz_time_climoz(read_climoz,interpt)
!
!-------------------------------------------------------------------------------
! Purpose: Regrid horizontally and in time zonal or 3D ozone climatologies.
!   * Read ozone climatology from netcdf file
!   * Regrid it horizontaly to LMDZ grid (quasi-conservative method)
!   * If interpt=T, interpolate linearly in time (one record each day)
!     If interpt=F, keep original time sampling  (14 months).
!   * Save it to a new netcdf file.
!-------------------------------------------------------------------------------
! Remarks:
!   * Up to 2 variables treated: "tro3" and "tro3_daylight" (if read_climoz=2)
!   * Input fields coordinates: (longitudes, latitudes, pressure_levels, time)
!   * Output grid cells centers coordinates given by [rlonv,] rlatu.
!   * Output grid cells edges   coordinates given by [rlonu,] rlatv.
!   * Input file [longitudes and] latitudes given in degrees.
!   * Input file pressure levels are given in Pa or hPa.
!   * All coordinates variables are stricly monotonic.
!   * Monthly fields are interpolated linearly in time to get daily values.
!   * Fields are known at the middle of the months, so interpolation requires an
!     additional record both for 1st half of january and 2nd half of december:
!     - For a 14-records "climoz.nc": records 1 and 14.
!     - For 12-records files:
!       record 12 of "climoz_m.nc" if available, or record 1  of "climoz.nc".
!       record 1  of "climoz_p.nc" if available, or record 12 of "climoz.nc".
!   * Calendar is taken into account to get one record each day (not 360 always).
!   * Missing values are filled in from sky to ground by copying lowest valid one.
!     Attribute "missing_value" or "_FillValue" must be present in input file.
!-------------------------------------------------------------------------------
  USE assert_m,           ONLY: assert
  USE cal_tools_m,        ONLY: year_len, mid_month
!!  USE control_mod,        ONLY: anneeref
  USE time_phylmdz_mod,   ONLY: annee_ref
  USE ioipsl,             ONLY: ioget_year_len, ioget_calendar
  USE regr_conserv_m,     ONLY: regr_conserv
  USE regr_lint_m,        ONLY: regr_lint
  USE regular_lonlat_mod, ONLY: boundslon_reg, boundslat_reg, south, west, east
  USE slopes_m,           ONLY: slopes
#ifdef CPP_XIOS
  USE xios
#endif
  USE mod_phys_lmdz_para, ONLY: is_mpi_root, is_master, is_omp_master, gather, gather_mpi, bcast_mpi, klon_mpi
  USE geometry_mod, ONLY : latitude_deg, ind_cell_glo
  USE mod_grid_phy_lmdz, ONLY: klon_glo

!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN) :: read_climoz ! read ozone climatology, 1 or 2
!                         1: read a single ozone climatology used day and night
!                         2: same + read also a daylight climatology
  LOGICAL, INTENT(IN) :: interpt     ! TRUE  => daily interpolation
                                     ! FALSE => no interpolation (14 months)
!-------------------------------------------------------------------------------
! Local variables:

!--- Input files variables
  INTEGER :: nlon_in                       ! Number of longitudes
  INTEGER :: nlat_in                       ! Number of latitudes
  INTEGER :: nlev_in                       ! Number of pressure levels
  INTEGER :: nmth_in                       ! Number of months
  REAL, POINTER     :: lon_in(:)           ! Longitudes   (ascending order, rad)
  REAL, POINTER     :: lat_in(:)           ! Latitudes    (ascending order, rad)
  REAL, POINTER     :: lev_in(:)           ! Pressure levels (ascen. order, hPa)
  REAL, ALLOCATABLE :: lon_in_edge(:)      ! Longitude intervals edges
                                           !              (ascending order,  / )
  REAL, ALLOCATABLE :: sinlat_in_edge(:)   ! Sinus of latitude intervals edges
                                           !              (ascending order,  / )
  LOGICAL :: ldec_lon, ldec_lat, ldec_lev  ! Decreasing order in input file
  CHARACTER(LEN=20) :: cal_in              ! Calendar
  REAL, ALLOCATABLE :: o3_in3(:,:,:,:,:)   ! Ozone climatologies
  REAL, ALLOCATABLE :: o3_in3bis(:,:,:,:,:)   ! Ozone climatologies
  REAL, ALLOCATABLE :: o3_in2  (:,:,:,:)   ! Ozone climatologies
  REAL, ALLOCATABLE :: o3_in2bis(:,:,:,:,:)   ! Ozone climatologies
  ! last index: 1 for the day-night average, 2 for the daylight field.
  REAL :: NaN

!--- Partially or totally regridded variables      (:,:,nlev_in,:,read_climoz)
  REAL, ALLOCATABLE :: o3_regr_lon   (:,:,:,:,:) ! (nlon_ou,nlat_in,:,0:13   ,:)
  REAL, ALLOCATABLE :: o3_regr_lonlat(:,:,:,:,:) ! (nlon_ou,nlat_ou,:,0:13   ,:)
  REAL, ALLOCATABLE :: o3_out3       (:,:,:,:,:) ! (nlon_ou,nlat_ou,:,ntim_ou,:)
  REAL, ALLOCATABLE :: o3_out3_glo   (:,:,:,:) !   (nbp_lat,:,ntim_ou,:)
  REAL, ALLOCATABLE :: o3_regr_lat     (:,:,:,:) !         (nlat_in,:,0:13   ,:)
  REAL, ALLOCATABLE :: o3_out2         (:,:,:,:) !         (nlat_ou,:,ntim_ou,:)
  REAL, ALLOCATABLE :: o3_out2_glo     (:,:,:,:) !         (nbp_lat,:,ntim_ou,:)
  REAL, ALLOCATABLE :: o3_out          (:,:,:,:) !         (nbp_lat,:,ntim_ou,:)
! Dimension number  | Interval                | Contains  | For variables:
!   1 (longitude)   | [rlonu(i-1), rlonu(i)]  | rlonv(i)  | all
!   2 (latitude)    | [rlatv(j), rlatv(j-1)]  | rlatu(j)  | all but o3_regr_lon
!   3 (press level) |                         |   lev(k)  | all
! Note that rlatv(0)=pi/2 and rlatv(nlat_ou)=-pi/2.
! Dimension 4 is: month number                             (all vars but o3_out)
!                 days elapsed since Jan. 1st 0h at mid-day (o3_out only)
  REAL, ALLOCATABLE :: v1(:)

!--- For NetCDF:
  INTEGER :: fID_in_m, fID_in, levID_ou, dimid, vID_in(read_climoz), ntim_ou
  INTEGER :: fID_in_p, fID_ou, timID_ou, varid, vID_ou(read_climoz), ndims, ncerr
  INTEGER, POINTER :: dIDs(:)
  CHARACTER(LEN=20) :: cal_ou     !--- Calendar; no time inter => same as input
  CHARACTER(LEN=80) :: press_unit !--- Pressure unit
  REAL    :: tmidmonth(0:13)      !--- Elapsed days since Jan-1 0h at mid-months
                                  ! Additional records 0, 13 for interpolation
  REAL, ALLOCATABLE :: tmidday(:) !--- Output times (mid-days since Jan 1st 0h)
  LOGICAL :: lprev, lnext         !--- Flags: previous/next files are present
  LOGICAL :: l3D, l2D             !--- Flag:  input fields are 3D or zonal
  INTEGER :: ii, i, j, k, l, m, dln, ib, ie, iv, dx1, dx2
  INTEGER, ALLOCATABLE :: sta(:), cnt(:)
  CHARACTER(LEN=80) :: sub, dim_nam, msg
  REAL :: null_array(0)
  LOGICAL,SAVE :: first=.TRUE.
!$OMP THREADPRIVATE(first)  
  REAL, ALLOCATABLE :: test_o3_in(:,:)
  REAL, ALLOCATABLE :: test_o3_out(:)


  IF (grid_type==unstructured) THEN 
    IF (first) THEN
      IF (is_master) THEN
        ALLOCATE(latitude_glo(klon_glo))
        ALLOCATE(ind_cell_glo_glo(klon_glo))
      ELSE
        ALLOCATE(latitude_glo(0))
        ALLOCATE(ind_cell_glo_glo(0))
      ENDIF
      CALL gather(latitude_deg,  latitude_glo)
      CALL gather(ind_cell_glo,  ind_cell_glo_glo)
    ENDIF 
  ENDIF
    
  IF (is_omp_master) THEN
    nlat_ou=nbp_lat
    nlon_ou=nbp_lon
    
   !-------------------------------------------------------------------------------
    IF (is_mpi_root) THEN
      sub="regr_horiz_time_climoz"
      WRITE(lunout,*)"Call sequence information: "//TRIM(sub)
      CALL assert(read_climoz == 1 .OR. read_climoz == 2, "regr_lat_time_climoz")

      CALL  NF95_OPEN("climoz.nc"  , NF90_NOWRITE, fID_in)
      lprev=NF90_OPEN("climoz_m.nc", NF90_NOWRITE, fID_in_m)==NF90_NOERR
      lnext=NF90_OPEN("climoz_p.nc", NF90_NOWRITE, fID_in_p)==NF90_NOERR

      !--- Get coordinates from the input file. Converts lon/lat in radians.
      !    Few inversions because "regr_conserv" and gcm need ascending vectors.
      CALL NF95_INQ_VARID(fID_in, vars_in(1), varid)
      CALL NF95_INQUIRE_VARIABLE(fID_in, varid, dimids=dIDs, ndims=ndims)
      l3D=ndims==4; l2D=ndims==3
      IF(l3D) WRITE(lunout,*)"Input files contain full 3D ozone fields."
      IF(l2D) WRITE(lunout,*)"Input files contain zonal 2D ozone fields."
      DO i=1,ndims
        CALL NF95_INQUIRE_DIMENSION(fID_in, dIDs(i), name=dim_nam, nclen=dln)
        CALL NF95_INQ_VARID(fID_in, dim_nam, varid)
        ii=i; IF(l2D) ii=i+1                              !--- ndims==3:NO LONGITUDE 
        SELECT CASE(ii)
          CASE(1)                                         !--- LONGITUDE
            CALL NF95_GW_VAR(fID_in, varid, lon_in)
            ldec_lon=lon_in(1)>lon_in(dln); IF(ldec_lon) lon_in=lon_in(dln:1:-1)
            nlon_in=dln; lon_in=lon_in*deg2rad
          CASE(2)                                         !--- LATITUDE
            CALL NF95_GW_VAR(fID_in, varid, lat_in)
            ldec_lat=lat_in(1)>lat_in(dln); IF(ldec_lat) lat_in=lat_in(dln:1:-1)
            nlat_in=dln; lat_in=lat_in*deg2rad
          CASE(3)                                         !--- PRESSURE LEVELS
            CALL NF95_GW_VAR(fID_in, varid, lev_in)
            ldec_lev=lev_in(1)>lev_in(dln); IF(ldec_lev) lev_in=lev_in(dln:1:-1)
            nlev_in=dln
            CALL NF95_GET_ATT(fID_in, varid, "units", press_unit)
            k=LEN_TRIM(press_unit)
            DO WHILE(ICHAR(press_unit(k:k))==0)
              press_unit(k:k)=' '; k=LEN_TRIM(press_unit) !--- REMOVE NULL END CHAR
            END DO
            IF(press_unit ==  "Pa") THEN
              lev_in = lev_in/100.                        !--- CONVERT TO hPa
            ELSE IF(press_unit /= "hPa") THEN
              CALL abort_physic(sub, "the only recognized units are Pa and hPa.",1)
            END IF
          CASE(4)                                         !--- TIME
            CALL NF95_INQUIRE_DIMENSION(fID_in, dIDs(i), nclen=nmth_in)
            cal_in='gregorian'
            IF(NF90_GET_ATT(fID_in, varid, 'calendar', cal_in)/=NF90_NOERR)        &  
            WRITE(lunout,*)'WARNING: missing "calendar" attribute for "'//       &
            TRIM(dim_nam)//'" in "climoz.nc". Choosing default: "gregorian".'
            k=LEN_TRIM(cal_in)
            DO WHILE(ICHAR(cal_in(k:k))==0)
              cal_in(k:k)=' '; k=LEN_TRIM(cal_in)         !--- REMOVE NULL END CHAR
            END DO
        END SELECT
      END DO

      !--- Prepare quantities for time interpolation
      tmidmonth=mid_month(annee_ref, cal_in)
      IF(interpt) THEN
        ntim_ou=ioget_year_len(annee_ref)
        ALLOCATE(tmidday(ntim_ou))
        tmidday=[(REAL(k)-0.5,k=1,ntim_ou)]
        CALL ioget_calendar(cal_ou)
      ELSE
        ntim_ou=14
        cal_ou=cal_in
      END IF
    ENDIF

    IF (grid_type==unstructured) THEN
      CALL bcast_mpi(nlon_in)
      CALL bcast_mpi(nlat_in)
      CALL bcast_mpi(nlev_in)
      CALL bcast_mpi(l3d)
      CALL bcast_mpi(tmidmonth)
      CALL bcast_mpi(tmidday)
      CALL bcast_mpi(ntim_ou)

#ifdef CPP_XIOS    
      IF (is_mpi_root) THEN
        CALL xios_set_domain_attr("domain_climoz",nj_glo=nlat_in, nj=nlat_in, jbegin=0, latvalue_1d=lat_in/deg2rad)
        IF (l3D) THEN
          CALL xios_set_domain_attr("domain_climoz",ni_glo=nlon_in, ni=nlon_in, ibegin=0, lonvalue_1d=lon_in/deg2rad)
        ELSE
          CALL xios_set_domain_attr("domain_climoz",ni_glo=8, ni=8, ibegin=0, lonvalue_1d = (/ 0.,45.,90.,135.,180.,225.,270., 315. /))
        ENDIF
      ELSE
        CALL xios_set_domain_attr("domain_climoz",nj_glo=nlat_in, nj=0, jbegin=0, latvalue_1d=null_array )
        IF (l3D) THEN
          CALL xios_set_domain_attr("domain_climoz",ni_glo=nlon_in, ni=0, ibegin=0, lonvalue_1d=null_array)
        ELSE
          CALL xios_set_domain_attr("domain_climoz",ni_glo=8, ni=0, ibegin=0, lonvalue_1d=null_array)
        ENDIF
      ENDIF
      CALL  xios_set_axis_attr("axis_climoz", n_glo=nlev_in)
      CALL  xios_set_axis_attr("time_axis_climoz", n_glo=ntim_ou)
      CALL  xios_set_axis_attr("time_axis_climoz", n_glo=ntim_ou)
      CALL  xios_set_axis_attr("tr_climoz", n_glo=read_climoz)
      CALL  xios_set_field_attr("tro3_out", enabled=.TRUE.)
      CALL  xios_set_field_attr("tro3_out", enabled=.TRUE.)
#endif
      
      IF (first) THEN
        first=.FALSE.
        RETURN
      ENDIF
    ENDIF
    
    
    IF (is_mpi_root) THEN     
      !--- Longitudes management:
      !    * Need to shift data if the origin of input file longitudes /= -pi
      !    * Need to add some margin in longitude to ensure input interval contains
      !      all the output intervals => at least one longitudes slice has to be
      !      duplicated, possibly more for undersampling.
      IF(l3D) THEN
        IF (grid_type==unstructured) THEN
          dx2=0
        ELSE
          !--- Compute input edges longitudes vector (no end point yet)
          ALLOCATE(v1(nlon_in+1))
          v1(1)=(lon_in(nlon_in)+lon_in(1))/2.-pi
          FORALL(i=2:nlon_in) v1(i)=(lon_in(i-1)+lon_in(i))/2.
          v1(nlon_in+1)=v1(1)+2.*pi
          DEALLOCATE(lon_in)

          !--- Shift input longitudes vector until it contains first output point boundslon_reg(1,west)
          v1=v1+2*pi*REAL(FLOOR((boundslon_reg(1,west)-v1(1))/(2.*pi)))

          !--- Ensure first input longitudes interval contains first output point boundslon_reg(1,west)
          dx1=locate(v1,boundslon_reg(1,west))-1
          v1=CSHIFT(v1,SHIFT=dx1,DIM=1) 
          v1(nlon_in-dx1+2:)=v1(nlon_in-dx1+2:)+2.*pi
    
          !--- Extend input longitudes vector until last interval contains boundslon_reg(nlat_ou,east)
          dx2=0; DO WHILE(v1(1+dx2)+2.*pi<=boundslon_reg(nlon_ou,east)); dx2=dx2+1; END DO

          !--- Final edges longitudes vector (with margin and end point)
          ALLOCATE(lon_in_edge(nlon_in+dx2+1)); lon_in_edge=[v1,v1(2:1+dx2)+2.*pi]
          DEALLOCATE(v1)
        ENDIF
      END IF

      !--- Compute sinus of intervals edges latitudes:
      ALLOCATE(sinlat_in_edge(nlat_in+1))
      sinlat_in_edge(1) = -1. ; sinlat_in_edge(nlat_in+1) = 1.
      FORALL(j=2:nlat_in) sinlat_in_edge(j)=SIN((lat_in(j-1)+lat_in(j))/2.)
      DEALLOCATE(lat_in)



      !--- Check for contiguous years:
      ib=0; ie=13
      IF(nmth_in == 14) THEN; lprev=.FALSE.; lnext=.FALSE.
        WRITE(lunout,*)'Using 14 months ozone climatology "climoz.nc"...'
      ELSE  
        IF(     lprev) WRITE(lunout,*)'Using "climoz_m.nc" last record (previous year).'
        IF(.NOT.lprev) WRITE(lunout,*)"No previous year file ; assuming periodicity."
        IF(     lnext) WRITE(lunout,*)'Using "climoz_p.nc" first record (next year).'
        IF(.NOT.lnext) WRITE(lunout,*)"No next year file ; assuming periodicity."
        IF(.NOT.lprev) ib=1
        IF(.NOT.lnext) ie=12
      END IF
      ALLOCATE(sta(ndims),cnt(ndims)); sta(:)=1  
      IF(l3D) cnt=[nlon_in,nlat_in,nlev_in,1]
      IF(l2D) cnt=[        nlat_in,nlev_in,1]  
      IF(l3D) ALLOCATE(o3_in3(nlon_in+dx2,nlat_in,nlev_in,ib:ie,read_climoz))
      IF(l2D) ALLOCATE(o3_in2(            nlat_in,nlev_in,ib:ie,read_climoz))

      !--- Read full current file and one record each available contiguous file
      DO iv=1,read_climoz
        msg=TRIM(sub)//" NF90_GET_VAR "//TRIM(vars_in(iv))
        CALL NF95_INQ_VARID(fID_in, vars_in(1), vID_in(iv))
        IF(l3D) ncerr=NF90_GET_VAR(fID_in, vID_in(iv), o3_in3(1:nlon_in,:,:,1:12,iv))
        IF(l2D) ncerr=NF90_GET_VAR(fID_in, vID_in(iv), o3_in2(          :,:,1:12,iv))
        CALL handle_err(TRIM(msg), ncerr, fID_in)
        IF(lprev) THEN; sta(ndims)=12  
          CALL NF95_INQ_VARID(fID_in_m, vars_in(1), vID_in(iv))
          IF(l3D) ncerr=NF90_GET_VAR(fID_in_m,vID_in(iv),o3_in3(1:nlon_in,:,:, 0,iv),sta,cnt)
          IF(l2d) ncerr=NF90_GET_VAR(fID_in_m,vID_in(iv),o3_in2(          :,:, 0,iv),sta,cnt)
          CALL handle_err(TRIM(msg)//" previous", ncerr, fID_in_m)
        END IF
        IF(lnext) THEN; sta(ndims)=1  
          CALL NF95_INQ_VARID(fID_in_p, vars_in(1), vID_in(iv))
          IF(l3D) ncerr=NF90_GET_VAR(fID_in_p,vID_in(iv),o3_in3(1:nlon_in,:,:,13,iv),sta,cnt)
          IF(l2D) ncerr=NF90_GET_VAR(fID_in_p,vID_in(iv),o3_in2(          :,:,13,iv),sta,cnt)
          CALL handle_err(TRIM(msg)//" next", ncerr, fID_in_p)
        END IF
      END DO
      IF(lprev.OR.lnext) DEALLOCATE(sta,cnt)
      IF(lprev) CALL NF95_CLOSE(fID_in_m)
      IF(lnext) CALL NF95_CLOSE(fID_in_p)

      !--- Revert decreasing coordinates vector
      IF(l3D) THEN
        IF(ldec_lon) o3_in3(1:nlon_in,:,:,:,:) = o3_in3(nlon_in:1:-1,:,:,:,:)
        IF(ldec_lat) o3_in3 = o3_in3(:,nlat_in:1:-1,:,:,:)
        IF(ldec_lev) o3_in3 = o3_in3(:,:,nlev_in:1:-1,:,:)
        
        IF (grid_type /= unstructured) THEN
          !--- Shift values for longitude and duplicate some longitudes slices
          o3_in3(1:nlon_in,:,:,:,:)=CSHIFT(o3_in3(1:nlon_in,:,:,:,:),SHIFT=dx1,DIM=1)
          o3_in3(nlon_in+1:nlon_in+dx2,:,:,:,:)=o3_in3(1:dx2,:,:,:,:)
        ENDIF
      ELSE
        IF(ldec_lat) o3_in2 = o3_in2(  nlat_in:1:-1,:,:,:)
        IF(ldec_lev) o3_in2 = o3_in2(  :,nlev_in:1:-1,:,:)
      END IF

     !--- Deal with missing values 
      DO m=1, read_climoz
        WRITE(msg,'(a,i0)')"regr_lat_time_climoz: field Nr.",m
        IF(NF90_GET_ATT(fID_in,vID_in(m),"missing_value",NaN)/= NF90_NOERR) THEN
          IF(NF90_GET_ATT(fID_in, vID_in(m),"_FillValue",NaN)/= NF90_NOERR) THEN
            WRITE(lunout,*)TRIM(msg)//": no missing value attribute found."; CYCLE
          END IF
        END IF
        WRITE(lunout,*)TRIM(msg)//": missing value attribute found."
        WRITE(lunout,*)"Trying to fill in NaNs ; a full field would be better."

        !--- Check top layer contains no NaNs & search NaNs from top to ground
        msg=TRIM(sub)//": NaNs in top layer !"
        IF(l3D) THEN
          IF(ANY(o3_in3(:,:,1,:,m)==NaN)) CALL abort_physic(sub,msg,1)
          DO k = 2,nlev_in
            WHERE(o3_in3(:,:,k,:,m)==NaN) o3_in3(:,:,k,:,m)=o3_in3(:,:,k-1,:,m)
          END DO
        ELSE
          IF(ANY(o3_in2(  :,1,:,m)==NaN)) THEN
            WRITE(lunout,*)msg
            !--- Fill in latitudes where all values are missing
            DO l=1,nmth_in
              !--- Next to south pole
              j=1;       DO WHILE(o3_in2(j,1,l,m)==NaN); j=j+1; END DO
              IF(j>1) &
                o3_in2(:j-1,:,l,m)=SPREAD(o3_in2(j,:,l,m),DIM=1,ncopies=j-1)
              !--- Next to north pole
              j=nlat_in; DO WHILE(o3_in2(j,1,l,m)==NaN); j=j+1; END DO
              IF(j<nlat_in) &
                o3_in2(j+1:,:,l,m)=SPREAD(o3_in2(j,:,l,m),DIM=1,ncopies=nlat_in-j)
            END DO
          END IF

          !--- Fill in high latitudes missing values
          !--- Highest level been filled-in, so has always valid values.
          DO k = 2,nlev_in
            WHERE(o3_in2(:,k,:,m)==NaN) o3_in2(:,k,:,m)=o3_in2(:,k-1,:,m)
          END DO
        END IF
      END DO

    ENDIF
    
    !=============================================================================
    IF(l3D) THEN                                                   !=== 3D FIELDS
    !=============================================================================
     IF (grid_type==unstructured) THEN
#ifdef CPP_XIOS
       nlat_ou=klon_mpi
       
       IF (is_mpi_root) THEN 
         ALLOCATE(o3_in3bis(nlon_in,nlat_in,nlev_in,0:13,read_climoz)) 
         o3_in3bis(:,:,:,ib:ie,:)=o3_in3(1:nlon_in,:,:,ib:ie,:)
       ELSE
         ALLOCATE(o3_in3bis(0,0,0,0,read_climoz))
       ENDIF
       ALLOCATE(o3_regr_lonlat(1, nlat_ou, nlev_in, 0:13, read_climoz))
       
       CALL xios_send_field("tro3_in",o3_in3bis(:,:,:,:,:))
       CALL xios_recv_field("tro3_out",o3_regr_lonlat(1,:,:,:,:))
#endif
     ELSE
          
       !--- Regrid in longitude
        ALLOCATE(o3_regr_lon(nlon_ou, nlat_in, nlev_in, ie-ib+1, read_climoz))
        CALL regr_conserv(1, o3_in3, xs = lon_in_edge,                             &
                            xt = [boundslon_reg(1,west),boundslon_reg(:,east)],    &
                            vt = o3_regr_lon, slope = slopes(1,o3_in3, lon_in_edge))
        DEALLOCATE(o3_in3)

        !--- Regrid in latitude: averaging with respect to SIN(lat) is
        !                        equivalent to weighting by COS(lat)
        !--- (inverted indices in "o3_regr_lonlat" because "rlatu" is decreasing)
        ALLOCATE(o3_regr_lonlat(nlon_ou, nlat_ou, nlev_in, 0:13, read_climoz))
        CALL regr_conserv(2, o3_regr_lon, xs = sinlat_in_edge,                     &
                        xt = [- 1., SIN(boundslat_reg(nlat_ou-1:1:-1,south)), 1.], &
                        vt = o3_regr_lonlat(:,nlat_ou:1:- 1,:,ib:ie,:),            &
                   slope = slopes(2,o3_regr_lon, sinlat_in_edge))
        DEALLOCATE(o3_regr_lon)

     ENDIF

     !--- Duplicate previous/next record(s) if they are not available
     IF(.NOT.lprev) o3_regr_lonlat(:,:,:, 0,:) = o3_regr_lonlat(:,:,:,12,:)
     IF(.NOT.lnext) o3_regr_lonlat(:,:,:,13,:) = o3_regr_lonlat(:,:,:, 1,:)
      
     !--- Regrid in time by linear interpolation:
     ALLOCATE(o3_out3(nlon_ou, nlat_ou, nlev_in, ntim_ou, read_climoz))
     IF(     interpt) CALL regr_lint(4,o3_regr_lonlat,tmidmonth,tmidday,o3_out3)
     IF(.NOT.interpt) o3_out3=o3_regr_lonlat
     DEALLOCATE(o3_regr_lonlat)

     nlat_ou=nbp_lat
     IF (grid_type==unstructured) THEN
#ifdef CPP_XIOS
       CALL xios_send_field('o3_out',o3_out3)
       ndims=3 
       ALLOCATE(o3_out3_glo(nlat_ou, nlev_in, ntim_ou, read_climoz))
       CALL gather_mpi(o3_out3(1,:,:,:,:), o3_out3_glo)
#endif
     ENDIF

    !--- Create the output file and get the variable IDs:
    CALL prepare_out(fID_in,nlev_in,ntim_ou, fID_ou,levID_ou,timID_ou,vID_ou, &
                     ndims, cal_ou)

    IF (is_mpi_root) THEN
      !--- Write remaining coordinate variables:
      CALL NF95_PUT_VAR(fID_ou, levID_ou, lev_in); DEALLOCATE(lev_in)
      IF(     interpt) CALL NF95_PUT_VAR(fID_ou, timID_ou, tmidday)
      IF(.NOT.interpt) CALL NF95_PUT_VAR(fID_ou, timID_ou, tmidmonth)

      !--- Write to file (the order of "rlatu" is inverted in the output file):
        IF (grid_type==unstructured) THEN

          ALLOCATE(o3_out(nlat_ou, nlev_in, ntim_ou, read_climoz))
          DO i=1,klon_glo
            o3_out(ind_cell_glo_glo(i),:,:,:)=o3_out3_glo(i,:,:,:)
          ENDDO

          DO m = 1, read_climoz
            CALL NF95_PUT_VAR(fID_ou, vID_ou(m), o3_out(nlat_ou:1:-1,:,:,m))
          END DO
          
        ELSE
          DO m = 1, read_climoz
            CALL NF95_PUT_VAR(fID_ou, vID_ou(m), o3_out3(:,nlat_ou:1:-1,:,:,m))
          END DO
      ENDIF
      CALL NF95_CLOSE(fID_ou)


    ENDIF


    !=============================================================================
    ELSE                                                         !=== ZONAL FIELDS
    !=============================================================================
   
     IF (grid_type==unstructured) THEN
#ifdef CPP_XIOS
       nlat_ou=klon_mpi

       IF (is_mpi_root) THEN 
         ALLOCATE(o3_in2bis(8,nlat_in,nlev_in,0:13,read_climoz)) 
         o3_in2bis(:,:,:,ib:ie,:)=SPREAD(o3_in2,1,8)
       ELSE
         ALLOCATE(o3_in2bis(0,0,0,0,read_climoz))
       ENDIF
       ALLOCATE(o3_regr_lat(nlat_ou, nlev_in, 0:13, read_climoz))
       CALL xios_send_field("tro3_in",o3_in2bis(:,:,:,:,:))
       CALL xios_recv_field("tro3_out",o3_regr_lat(:,:,:,:))
       IF(.NOT.lprev) o3_regr_lat(:,:, 0, :) = o3_regr_lat(:,:,12,:)
       IF(.NOT.lnext) o3_regr_lat(:,:,13, :) = o3_regr_lat(:,:, 1,:)
#endif       
     
     ELSE 
        !--- Regrid in latitude: averaging with respect to SIN(lat) is
        !                        equivalent to weighting by COS(lat)
        !--- (inverted indices in "o3_regr_lat" because "rlatu" is decreasing)
        ALLOCATE(o3_regr_lat(nlat_ou, nlev_in, 0:13, read_climoz))
        CALL regr_conserv(1, o3_in2, xs = sinlat_in_edge,                          &
                        xt = [- 1., SIN(boundslat_reg(nlat_ou-1:1:-1,south)), 1.], &
                        vt = o3_regr_lat(nlat_ou:1:- 1,:,ib:ie,:),                 &
                     slope = slopes(1,o3_in2, sinlat_in_edge))
        DEALLOCATE(o3_in2)

        !--- Duplicate previous/next record(s) if they are not available
        IF(.NOT.lprev) o3_regr_lat(:,:, 0,:) = o3_regr_lat(:,:,12,:)
        IF(.NOT.lnext) o3_regr_lat(:,:,13,:) = o3_regr_lat(:,:, 1,:)

     ENDIF
     
      !--- Regrid in time by linear interpolation:
      ALLOCATE(o3_out2(nlat_ou, nlev_in, ntim_ou, read_climoz))
      IF(     interpt) CALL regr_lint(3,o3_regr_lat, tmidmonth, tmidday, o3_out2)
      IF(.NOT.interpt) o3_out2=o3_regr_lat
      DEALLOCATE(o3_regr_lat)

      nlat_ou=nbp_lat
    
      IF (grid_type==unstructured) THEN
        ndims=3 
        ALLOCATE(o3_out2_glo(nlat_ou, nlev_in, ntim_ou, read_climoz))
        CALL gather_mpi(o3_out2, o3_out2_glo)
      ENDIF
      
      !--- Create the output file and get the variable IDs:
      CALL prepare_out(fID_in,nlev_in,ntim_ou, fID_ou,levID_ou,timID_ou,vID_ou, &
                         ndims, cal_ou)

      IF (is_mpi_root) THEN
      
        !--- Write remaining coordinate variables:
        CALL NF95_PUT_VAR(fID_ou, levID_ou, lev_in); DEALLOCATE(lev_in)
        IF(     interpt) CALL NF95_PUT_VAR(fID_ou, timID_ou, tmidday)
        IF(.NOT.interpt) CALL NF95_PUT_VAR(fID_ou, timID_ou, tmidmonth)

        IF (grid_type==unstructured) THEN

          ALLOCATE(o3_out3_glo(nlat_ou, nlev_in, ntim_ou, read_climoz))
          DO i=1,klon_glo
            o3_out(ind_cell_glo_glo(i),:,:,:)=o3_out2_glo(i,:,:,:)
          ENDDO


          DO m = 1, read_climoz
            CALL NF95_PUT_VAR(fID_ou, vID_ou(m), o3_out(nlat_ou:1:-1,:,:,m))
          END DO
        ELSE
          !--- Write to file (the order of "rlatu" is inverted in the output file):
          DO m = 1, read_climoz
            CALL NF95_PUT_VAR(fID_ou, vID_ou(m), o3_out2(nlat_ou:1:-1,:,:,m))
          END DO
        ENDIF
        
        CALL NF95_CLOSE(fID_ou)
      
      ENDIF

    !=============================================================================
    END IF
    !=============================================================================

    IF (is_mpi_root) CALL NF95_CLOSE(fID_in)

  ENDIF ! is_omp_master

  first=.FALSE.
END SUBROUTINE regr_horiz_time_climoz
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE prepare_out(fID_in, nlev_in, ntim_ou, fID_ou, vlevID, vtimID, &
                       vID_ou, ndims, cal_ou)
!-------------------------------------------------------------------------------
! Purpose:  This subroutine creates the NetCDF output file, defines
!     dimensions and variables, and writes some of the coordinate variables.
!-------------------------------------------------------------------------------
  USE regular_lonlat_mod, ONLY: lon_reg, lat_reg
  USE regular_lonlat_mod, ONLY: lon_reg, lat_reg
  USE mod_phys_lmdz_para, ONLY: is_mpi_root
  USE mod_grid_phy_lmdz, ONLY: klon_glo
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER, INTENT(IN)  :: fID_in, nlev_in, ntim_ou
  INTEGER, INTENT(OUT) :: fID_ou, vlevID,  vtimID
  INTEGER, INTENT(OUT) :: vID_ou(:)      ! dim(1/2) 1: O3day&night 2: O3daylight
  INTEGER, INTENT(IN)  :: ndims          ! fields rank (3 or 4)
  CHARACTER(LEN=*), INTENT(IN) :: cal_ou ! calendar
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: dlonID, dlatID, dlevID, dtimID, dIDs(4)
  INTEGER :: vlonID, vlatID, ncerr,  is
  REAL,ALLOCATABLE    :: latitude_glo_(:)
  CHARACTER(LEN=80) :: sub
  INTEGER :: i


!-------------------------------------------------------------------------------
  
  IF (is_mpi_root) THEN  
    sub="prepare_out"
    WRITE(lunout,*)"CALL sequence information: "//TRIM(sub)
    CALL NF95_CREATE("climoz_LMDZ.nc", NF90_clobber, fID_ou)

  !--- Dimensions:
    IF(ndims==4) &
    CALL NF95_DEF_DIM(fID_ou, "rlonv", nlon_ou, dlonID)
    CALL NF95_DEF_DIM(fID_ou, "rlatu", nlat_ou, dlatID)
    CALL NF95_DEF_DIM(fID_ou, "plev",  nlev_in, dlevID)
    CALL NF95_DEF_DIM(fID_ou, "time",  ntim_ou, dtimID)

    !--- Define coordinate variables:
    IF(ndims==4) &
    CALL NF95_DEF_VAR(fID_ou, "rlonv", NF90_FLOAT, dlonID, vlonID)
    CALL NF95_DEF_VAR(fID_ou, "rlatu", NF90_FLOAT, dlatID, vlatID)
    CALL NF95_DEF_VAR(fID_ou, "plev",  NF90_FLOAT, dlevID, vlevID)
    CALL NF95_DEF_VAR(fID_ou, "time",  NF90_FLOAT, dtimID, vtimID)
    IF(ndims==4) &
    CALL NF95_PUT_ATT(fID_ou, vlonID, "units", "degrees_east")
    CALL NF95_PUT_ATT(fID_ou, vlatID, "units", "degrees_north")
    CALL NF95_PUT_ATT(fID_ou, vlevID, "units", "millibar")
    CALL NF95_PUT_ATT(fID_ou, vtimID, "units", "days since 2000-1-1")
    IF(ndims==4) &
    CALL NF95_PUT_ATT(fID_ou, vlonID, "standard_name", "longitude")
    CALL NF95_PUT_ATT(fID_ou, vlatID, "standard_name", "latitude")
    CALL NF95_PUT_ATT(fID_ou, vlevID, "standard_name", "air_pressure")
    CALL NF95_PUT_ATT(fID_ou, vtimID, "standard_name", "time")
    CALL NF95_PUT_ATT(fID_ou, vlevID, "long_name",     "air pressure")
    CALL NF95_PUT_ATT(fID_ou, vtimID, "calendar",      cal_ou)

  !--- Define the main variables:
    IF(ndims==3) dIDs(1:3) = [ dlatID, dlevID, dtimID]
    IF(ndims==4) dIDs=[dlonID, dlatID, dlevID, dtimID]
    CALL NF95_DEF_VAR(fID_ou, vars_in(1), NF90_FLOAT, dIDs(1:ndims), vID_ou(1))
    CALL NF95_PUT_ATT(fID_ou, vID_ou(1), "long_name", "ozone mole fraction")
    CALL NF95_PUT_ATT(fID_ou, vID_ou(1), "standard_name", "mole_fraction_of_ozone&
      &_in_air")
    IF(SIZE(vID_ou) == 2) THEN
      CALL NF95_DEF_VAR(fID_ou, vars_in(2), NF90_FLOAT, dIDs(1:ndims), vID_ou(2))
      CALL NF95_PUT_ATT(fID_ou, vID_ou(2), "long_name","ozone mole fraction in da&
        &ylight")
    END IF

  !--- Global attributes:
  ! The following commands, copying attributes, may fail. That is OK.
  ! It should just mean that the attribute is not defined in the input file.
    CALL NF95_COPY_ATT(fID_in,NF90_GLOBAL,"Conventions",fID_ou,NF90_GLOBAL, ncerr)
    CALL handle_err_copy_att("Conventions")
    CALL NF95_COPY_ATT(fID_in,NF90_GLOBAL,"title",      fID_ou,NF90_GLOBAL, ncerr)
    CALL handle_err_copy_att("title")
    CALL NF95_COPY_ATT(fID_in,NF90_GLOBAL,"institution",fID_ou,NF90_GLOBAL, ncerr)
    CALL handle_err_copy_att("institution")
    CALL NF95_COPY_ATT(fID_in,NF90_GLOBAL,"source",     fID_ou,NF90_GLOBAL, ncerr)
    CALL handle_err_copy_att("source")
    CALL NF95_PUT_ATT (fID_ou,NF90_GLOBAL,"comment", "Regridded for LMDZ")
    CALL NF95_ENDDEF(fID_ou)

    IF (grid_type==unstructured) THEN
      ALLOCATE(latitude_glo_(klon_glo))
      DO i=1,klon_glo
        latitude_glo_(ind_cell_glo_glo(i))=latitude_glo(i)
      ENDDO
      CALL NF95_PUT_VAR(fID_ou, vlatID, latitude_glo_)
    ELSE
      !--- Write one of the coordinate variables:
      IF(ndims==4) CALL NF95_PUT_VAR(fID_ou, vlonID, lon_reg/deg2rad)
      CALL NF95_PUT_VAR(fID_ou, vlatID, lat_reg(nlat_ou:1:-1)/deg2rad)
    !    (convert from rad to degrees and sort in ascending order)
    ENDIF
  ENDIF
  
CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE handle_err_copy_att(att_name)
!
!-------------------------------------------------------------------------------
  USE netcdf, ONLY: NF90_NOERR, NF90_strerror
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: att_name
!-------------------------------------------------------------------------------
  IF(ncerr /= NF90_NOERR) &
    WRITE(lunout,*)TRIM(sub)//" prepare_out NF95_COPY_ATT "//TRIM(att_name)//  &
                      " -- "//TRIM(NF90_strerror(ncerr))

END SUBROUTINE handle_err_copy_att
!
!-------------------------------------------------------------------------------

END SUBROUTINE prepare_out
!
!-------------------------------------------------------------------------------

END MODULE regr_horiz_time_climoz_m
!
!-------------------------------------------------------------------------------
