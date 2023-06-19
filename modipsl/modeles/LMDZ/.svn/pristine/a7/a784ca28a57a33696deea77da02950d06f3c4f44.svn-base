PROGRAM ce0l
!
!-------------------------------------------------------------------------------
! Purpose: Initial states and boundary conditions files creation:
!     * start.nc    for dynamics    (using etat0dyn     routine)
!     * startphy.nc for physics     (using etat0phys    routine)
!     * limit.nc    for forced runs (using limit_netcdf routine)
!-------------------------------------------------------------------------------
! Notes:
!     * extrap=.T. (default) for data extrapolation, like for the SSTs when file
!                   does contain ocean points only.
!     * "masque" can be:
!       - read from file "o2a.nc"          (for coupled runs).
!       - read from file "startphy0.nc"    (from a previous run).
!       - created in etat0phys or etat0dyn (for forced  runs).
!     It is then passed to limit_netcdf to ensure consistancy.
!-------------------------------------------------------------------------------
  USE ioipsl, ONLY: ioconf_calendar, getin, flininfo, flinopen, flinget, flinclo
  USE control_mod,    ONLY: day_step, dayref, nsplit_phys
  USE etat0dyn,       ONLY: etat0dyn_netcdf
  USE etat0phys,      ONLY: etat0phys_netcdf
  USE limit,          ONLY: limit_netcdf
  USE netcdf,         ONLY: NF90_OPEN, NF90_NOWRITE, NF90_CLOSE, NF90_NOERR,    &
         NF90_INQUIRE_DIMENSION, NF90_INQ_DIMID, NF90_INQ_VARID, NF90_GET_VAR
  USE infotrac,       ONLY: type_trac, infotrac_init
  USE dimphy,         ONLY: klon
  USE test_disvert_m, ONLY: test_disvert
  USE filtreg_mod,    ONLY: inifilr
  USE iniphysiq_mod,  ONLY: iniphysiq
  USE mod_const_mpi,  ONLY: comm_lmdz
#ifdef CPP_PARA
  USE mod_const_mpi,  ONLY: init_const_mpi
  USE parallel_lmdz,  ONLY: init_parallel, mpi_rank, omp_rank
  USE bands,          ONLY: read_distrib, distrib_phys
  USE mod_hallo,      ONLY: init_mod_hallo
  USE mod_interface_dyn_phys, ONLY: init_interface_dyn_phys
#endif
  USE comconst_mod, ONLY: cpp, daysec, dtphys, dtvr, g, kappa, omeg, r, rad, &
                          pi, jmp1
  USE logic_mod, ONLY: iflag_phys, ok_etat0, ok_limit
  USE comvert_mod, ONLY: pa, preff, pressure_exner
  USE temps_mod, ONLY: calend, day_ini, dt

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Local variables:
  include "dimensions.h"
  include "paramet.h"
  include "comgeom2.h"
  include "iniprint.h"
  REAL               :: masque(iip1,jjp1)             !--- CONTINENTAL MASK
  REAL               :: phis  (iip1,jjp1)             !--- GROUND GEOPOTENTIAL
  CHARACTER(LEN=256) :: modname, fmt, calnd           !--- CALENDAR TYPE
  LOGICAL            :: use_filtre_fft
  LOGICAL, PARAMETER :: extrap=.FALSE.

!--- Local variables for ocean mask reading:
  INTEGER            :: nid_o2a, iml_omask, jml_omask, j
  INTEGER            :: fid, iret, llm_tmp, ttm_tmp, itaul(1)
  REAL, ALLOCATABLE  :: lon_omask(:,:), dlon_omask(:), ocemask(:,:)
  REAL, ALLOCATABLE  :: lat_omask(:,:), dlat_omask(:), ocetmp (:,:)
  REAL               :: date, lev(1)

!--- Local variables for land mask from startphy0 file reading
  INTEGER            :: nid_sta, nid_nph, nid_msk, nphys
  REAL, ALLOCATABLE  :: masktmp(:)

#ifndef CPP_PARA
! for iniphysiq in serial mode
  INTEGER,PARAMETER :: mpi_rank=0
  INTEGER :: distrib_phys(mpi_rank:mpi_rank)=(jjm-1)*iim+2
#endif
!-------------------------------------------------------------------------------
  modname="ce0l"

!--- Constants
  pi     = 4. * ATAN(1.)
  rad    = 6371229.
  daysec = 86400.
  omeg   = 2.*pi/daysec
  g      = 9.8
  kappa  = 0.2857143
  cpp    = 1004.70885
  jmp1   = jjm + 1
  preff   = 101325.
  pa      = 50000.

  CALL conf_gcm( 99, .TRUE. )
  dtvr = daysec/REAL(day_step)
  WRITE(lunout,*)'dtvr',dtvr
  CALL iniconst()
  CALL inigeom()

!--- Calendar choice
#ifdef CPP_IOIPSL
  calnd='gregorian'
  SELECT CASE(calend)
    CASE('earth_360d');CALL ioconf_calendar('360d');   calnd='with 360 days/year'
    CASE('earth_365d');CALL ioconf_calendar('noleap'); calnd='with no leap year'
    CASE('earth_366d');CALL ioconf_calendar('366d');   calnd='with leap years only'
    CASE('gregorian'); CALL ioconf_calendar('gregorian')
    CASE('standard');  CALL ioconf_calendar('gregorian')
    CASE('julian');    CALL ioconf_calendar('julian'); calnd='julian'
    CASE('proleptic_gregorian'); CALL ioconf_calendar('gregorian')
  !--- DC Bof...  => IOIPSL a mettre a jour: proleptic_gregorian /= gregorian
    CASE DEFAULT
      CALL abort_gcm('ce0l','Bad choice for calendar',1)
  END SELECT
  WRITE(lunout,*)'CHOSEN CALENDAR: Earth '//TRIM(calnd)
#endif

#ifdef CPP_PARA
!--- Physical grid + parallel initializations
  CALL init_const_mpi()
  CALL init_parallel()
  CALL read_distrib()
  CALL init_mod_hallo()
#endif
  WRITE(lunout,*)'---> klon=',klon

!--- Tracers initializations
  CALL infotrac_init()

  CALL inifilr()
  CALL iniphysiq(iim,jjm,llm, &
                 distrib_phys(mpi_rank),comm_lmdz, &
                 daysec,day_ini,dtphys/nsplit_phys, &
                 rlatu,rlatv,rlonu,rlonv,aire,cu,cv,rad,g,r,cpp,iflag_phys)
  IF(pressure_exner) CALL test_disvert

#ifdef CPP_PARA
  IF (mpi_rank==0.AND.omp_rank==0) THEN
#endif
  use_filtre_fft=.FALSE.
  CALL getin('use_filtre_fft',use_filtre_fft)
  IF(use_filtre_fft) THEN
     WRITE(lunout,*)"FFT filter not available for sequential dynamics."
     WRITE(lunout,*)"Your setting of variable use_filtre_fft is not used."
  ENDIF

!--- LAND MASK. THREE CASES:
!   1) read from ocean model    file "o2a.nc"    (coupled runs)
!   2) read from previous run   file="startphy0.nc"
!   3) computed from topography file "Relief.nc" (masque(:,:)=-99999.)
! In the first case, the mask from the ocean model is used compute the
! weights to ensure ocean fractions are the same for atmosphere and ocean.
!*******************************************************************************
  IF(NF90_OPEN("o2a.nc", NF90_NOWRITE, nid_o2a)==NF90_NOERR) THEN
    iret=NF90_CLOSE(nid_o2a)
    WRITE(lunout,*)'BEWARE !! Ocean mask "o2a.nc" file found'
    WRITE(lunout,*)'Coupled run.'
    CALL flininfo("o2a.nc", iml_omask, jml_omask, llm_tmp, ttm_tmp, nid_o2a)
    IF(iml_omask/=iim .OR.jml_omask/=jjp1) THEN
      WRITE(lunout,*)'Mismatching dimensions for ocean mask'
      WRITE(lunout,*)'iim  = ',iim ,' iml_omask = ',iml_omask
      WRITE(lunout,*)'jjp1 = ',jjp1,' jml_omask = ',jml_omask
      CALL abort_gcm(modname,'',1)
    END IF
    ALLOCATE(ocemask(iim,jjp1),lon_omask(iim,jjp1),dlon_omask(iim ))
    ALLOCATE(ocetmp (iim,jjp1),lat_omask(iim,jjp1),dlat_omask(jjp1))
    CALL flinopen("o2a.nc", .FALSE.,iml_omask,jml_omask,llm_tmp,               &
                  lon_omask,lat_omask,lev,ttm_tmp,itaul,date,dt,fid)
    CALL flinget(fid, "OceMask",    iim,jjp1,llm_tmp,ttm_tmp,1,1,ocetmp)
    CALL flinclo(fid)
    dlon_omask(1:iim ) = lon_omask(1:iim,1)
    dlat_omask(1:jjp1) = lat_omask(1,1:jjp1)
    ocemask = ocetmp
    IF(dlat_omask(1)<dlat_omask(jml_omask)) THEN
      DO j=1,jjp1; ocemask(:,j) = ocetmp(:,jjp1-j+1); END DO
    END IF
    DEALLOCATE(ocetmp,lon_omask,lat_omask,dlon_omask,dlat_omask)
    IF(prt_level>=1) THEN
      WRITE(fmt,"(i4,'i1)')")iim ; fmt='('//ADJUSTL(fmt)
      WRITE(lunout,*)'OCEAN MASK :'
      WRITE(lunout,fmt) NINT(ocemask)
    END IF
    masque(1:iim,:)=1.-ocemask(:,:)
    masque(iip1 ,:)=masque(1,:)
    DEALLOCATE(ocemask)
  ELSE IF(NF90_OPEN("startphy0.nc", NF90_NOWRITE, nid_sta)==NF90_NOERR) THEN
    WRITE(lunout,*)'BEWARE !! File "startphy0.nc" found.'
    WRITE(lunout,*)'Getting the land mask from a previous run.'
    iret=NF90_INQ_DIMID(nid_sta,'points_physiques',nid_nph)
    iret=NF90_INQUIRE_DIMENSION(nid_sta,nid_nph,len=nphys)
    IF(nphys/=klon) THEN
      WRITE(lunout,*)'Mismatching dimensions for land mask'
      WRITE(lunout,*)'nphys  = ',nphys ,' klon = ',klon
      iret=NF90_CLOSE(nid_sta)
      CALL abort_gcm(modname,'',1)
    END IF
    ALLOCATE(masktmp(klon))
    iret=NF90_INQ_VARID(nid_sta,'masque',nid_msk)
    iret=NF90_GET_VAR(nid_sta,nid_msk,masktmp)
    iret=NF90_CLOSE(nid_sta)
    CALL gr_fi_dyn(1,klon,iip1,jjp1,masktmp,masque)
    IF(prt_level>=1) THEN
      WRITE(fmt,"(i4,'i1)')")iip1 ; fmt='('//ADJUSTL(fmt)
      WRITE(lunout,*)'LAND MASK :'
      WRITE(lunout,fmt) NINT(masque)
    END IF
    DEALLOCATE(masktmp)
  ELSE
    WRITE(lunout,*)'BEWARE !! No ocean mask "o2a.nc" file or "startphy0.nc" file found'
    WRITE(lunout,*)'Land mask will be built from the topography file.'
    masque(:,:)=-99999.
  END IF
  phis(:,:)=-99999.

  IF(ok_etat0) THEN
    WRITE(lunout,'(//)')
    WRITE(lunout,*) '  ************************  '
    WRITE(lunout,*) '  ***  etat0phy_netcdf ***  '
    WRITE(lunout,*) '  ************************  '
    CALL etat0phys_netcdf(masque,phis)
    WRITE(lunout,'(//)')
    WRITE(lunout,*) '  ************************  '
    WRITE(lunout,*) '  ***  etat0dyn_netcdf ***  '
    WRITE(lunout,*) '  ************************  '
    CALL etat0dyn_netcdf(masque,phis)
  END IF

  IF(ok_limit) THEN
    WRITE(lunout,'(//)')
    WRITE(lunout,*) '  *********************  '
    WRITE(lunout,*) '  ***  Limit_netcdf ***  '
    WRITE(lunout,*) '  *********************  '
    WRITE(lunout,'(//)')
    CALL limit_netcdf(masque,phis,extrap)
  END IF

  WRITE(lunout,'(//)')
  WRITE(lunout,*) '  ***************************  '
  WRITE(lunout,*) '  ***  grilles_gcm_netcdf ***  '
  WRITE(lunout,*) '  ***************************  '
  WRITE(lunout,'(//)')
  CALL grilles_gcm_netcdf_sub(masque,phis)

#ifdef CPP_PARA
  END IF
#endif

END PROGRAM ce0l
!
!-------------------------------------------------------------------------------
