MODULE limit
!
!*******************************************************************************
! Author : L. Fairhead, 27/01/94
!-------------------------------------------------------------------------------
! Purpose: Boundary conditions files building for new model using climatologies.
!          Both grids have to be regular.
!-------------------------------------------------------------------------------
! Note: This routine is designed to work for Earth
!-------------------------------------------------------------------------------
! Modification history:
!  * 23/03/1994: Z. X. Li
!  *    09/1999: L. Fairhead (netcdf reading in LMDZ.3.3)
!  *    07/2001: P. Le Van
!  *    11/2009: L. Guez     (ozone day & night climatos, see etat0_netcdf.F90)
!  *    12/2009: D. Cugnet   (f77->f90, calendars, files from coupled runs)
!-------------------------------------------------------------------------------

  USE ioipsl,             ONLY: flininfo, flinopen, flinget, flinclo
  USE assert_eq_m,        ONLY: assert_eq
  USE cal_tools_m,        ONLY: year_len, mid_month
  USE conf_dat_m,         ONLY: conf_dat2d, conf_dat3d
  USE dimphy,             ONLY: klon, zmasq
  USE geometry_mod,       ONLY: longitude_deg, latitude_deg
  USE phys_state_var_mod, ONLY: pctsrf
  USE control_mod,        ONLY: anneeref
  USE init_ssrf_m,        ONLY: start_init_subsurf

  INTEGER,           PARAMETER :: ns=256
  CHARACTER(LEN=ns), PARAMETER :: &
  fsst(5)=['amipbc_sst_1x1.nc   ','amip_sst_1x1.nc     ','cpl_atm_sst.nc      '&
          ,'histmth_sst.nc      ','sstk.nc             '],                     &
  fsic(5)=['amipbc_sic_1x1.nc   ','amip_sic_1x1.nc     ','cpl_atm_sic.nc      '&
          ,'histmth_sic.nc      ','ci.nc               '],                     &
  vsst(5)=['tosbcs    ','tos       ','SISUTESW  ','tsol_oce  ','sstk      '],  &
  vsic(5)=['sicbcs    ','sic       ','SIICECOV  ','pourc_sic ','ci        '],  &
  frugo='Rugos.nc  ', falbe='Albedo.nc ', frelf='Relief.nc ',                  &
   vrug='RUGOS     ',  valb='ALBEDO    ',  vrel='RELIEF    ',                  &
  DegK(11)=['degK          ','degree_K      ','degreeK       ','deg_K         '&
           ,'degsK         ','degrees_K     ','degreesK      ','degs_K        '&
           ,'degree_kelvin ','degrees_kelvin','K             '],               &
  DegC(10)=['degC          ','degree_C      ','degreeC       ','deg_C         '&
           ,'degsC         ','degrees_C     ','degreesC      ','degs_C        '&
           ,'degree_Celsius','celsius       '], &
  Perc(2) =['%             ','percent       '], &
  Frac(2) =['1.0           ','1             ']

CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE limit_netcdf(masque, phis, extrap)
!
!-------------------------------------------------------------------------------
! Author : L. Fairhead, 27/01/94
!-------------------------------------------------------------------------------
! Purpose: Boundary conditions files building for new model using climatologies.
!          Both grids have to be regular.
!-------------------------------------------------------------------------------
! Note: This routine is designed to work for Earth
!-------------------------------------------------------------------------------
! Modification history:
!  * 23/03/1994: Z. X. Li
!  *    09/1999: L. Fairhead (netcdf reading in LMDZ.3.3)
!  *    07/2001: P. Le Van
!  *    11/2009: L. Guez     (ozone day & night climatos, see etat0_netcdf.F90)
!  *    12/2009: D. Cugnet   (f77->f90, calendars, files from coupled runs)
!  *    04/2016: D. Cugnet   (12/14 recs SST/SIC files: cyclic/interannual runs)
!  *    05/2017: D. Cugnet   (linear time interpolation for BCS files)
!-------------------------------------------------------------------------------
#ifndef CPP_1D
  USE indice_sol_mod
  USE netcdf,             ONLY: NF90_OPEN,    NF90_CREATE,  NF90_CLOSE,        &
                  NF90_DEF_DIM, NF90_DEF_VAR, NF90_PUT_VAR, NF90_PUT_ATT,      &
                  NF90_NOERR,   NF90_NOWRITE, NF90_DOUBLE,  NF90_GLOBAL,       &
                  NF90_CLOBBER, NF90_ENDDEF,  NF90_UNLIMITED, NF90_FLOAT
  USE inter_barxy_m,      ONLY: inter_barxy
  USE netcdf95,           ONLY: nf95_def_var, nf95_put_att, nf95_put_var
  USE comconst_mod, ONLY: pi
  USE phys_cal_mod, ONLY: calend
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  include "iniprint.h"
  include "dimensions.h"
  include "paramet.h"
  REAL, DIMENSION(iip1,jjp1), INTENT(INOUT) :: masque ! land mask
  REAL, DIMENSION(iip1,jjp1), INTENT(INOUT) :: phis   ! ground geopotential
  LOGICAL,                    INTENT(IN)    :: extrap ! SST extrapolation flag
!-------------------------------------------------------------------------------
! Local variables:
  include "comgeom2.h"

!--- INPUT NETCDF FILES AND VARIABLES NAMES ------------------------------------
  CHARACTER(LEN=ns) :: icefile, sstfile, fnam, varname

!--- OUTPUT VARIABLES FOR NETCDF FILE ------------------------------------------
  REAL               :: fi_ice(klon)
  REAL, POINTER      :: phy_rug(:,:)=>NULL(), phy_ice(:,:)=>NULL()
  REAL, POINTER      :: phy_sst(:,:)=>NULL(), phy_alb(:,:)=>NULL()
  REAL, ALLOCATABLE  :: phy_bil(:,:), pctsrf_t(:,:,:)
  INTEGER            :: nbad

!--- VARIABLES FOR OUTPUT FILE WRITING -----------------------------------------
  INTEGER :: nid, ndim, ntim, k, dims(2), ix_sic, ix_sst
  INTEGER :: id_tim,  id_SST,  id_BILS, id_RUG, id_ALB
  INTEGER :: id_FOCE, id_FSIC, id_FTER, id_FLIC, varid_longitude, varid_latitude
  INTEGER :: NF90_FORMAT
  INTEGER :: ndays                   !--- Depending on the output calendar
  CHARACTER(LEN=ns) :: str

!--- INITIALIZATIONS -----------------------------------------------------------
#ifdef NC_DOUBLE
  NF90_FORMAT=NF90_DOUBLE
#else
  NF90_FORMAT=NF90_FLOAT
#endif
  CALL inigeom

!--- MASK, GROUND GEOPOT. & SUBSURFACES COMPUTATION (IN CASE ok_etat0==.FALSE.)
   IF(ALL(masque==-99999.)) THEN
    CALL start_init_orog0(rlonv,rlatu,phis,masque)
    CALL gr_dyn_fi(1,iip1,jjp1,klon,masque,zmasq)          !--- To physical grid
    ALLOCATE(pctsrf(klon,nbsrf))
    CALL start_init_subsurf(.FALSE.)
  !--- TO MATCH EXACTLY WHAT WOULD BE DONE IN etat0phys_netcdf
    WHERE(   masque(:,:)<EPSFRA) masque(:,:)=0.
    WHERE(1.-masque(:,:)<EPSFRA) masque(:,:)=1.
  END IF

!--- Beware: anneeref (from gcm.def) is used to determine output time sampling
  ndays=year_len(anneeref)

!--- RUGOSITY TREATMENT --------------------------------------------------------
  CALL msg(0,""); CALL msg(0," *** TRAITEMENT DE LA RUGOSITE ***")
  CALL get_2Dfield(frugo,vrug,'RUG',ndays,phy_rug,mask=masque(1:iim,:))

!--- OCEAN TREATMENT -----------------------------------------------------------
  CALL msg(0,""); CALL msg(0," *** TRAITEMENT DE LA GLACE OCEANIQUE ***")

! Input SIC file selection
! Open file only to test if available
  DO ix_sic=1,SIZE(fsic)
     IF ( NF90_OPEN(TRIM(fsic(ix_sic)),NF90_NOWRITE,nid)==NF90_NOERR ) THEN
        icefile=fsic(ix_sic); varname=vsic(ix_sic); EXIT
     END IF
  END DO
  IF(ix_sic==SIZE(fsic)+1) THEN
     WRITE(lunout,*) 'ERROR! No sea-ice input file was found.'
     WRITE(lunout,*) 'One of following files must be available : '
     DO k=1,SIZE(fsic); WRITE(lunout,*) TRIM(fsic(k)); END DO
     CALL abort_physic('limit_netcdf','No sea-ice file was found',1)
  END IF
  CALL ncerr(NF90_CLOSE(nid),icefile)
  CALL msg(0,'Fichier choisi pour la glace de mer:'//TRIM(icefile))

  CALL get_2Dfield(icefile,varname, 'SIC',ndays,phy_ice)

  ALLOCATE(pctsrf_t(klon,nbsrf,ndays))
  DO k=1,ndays
     fi_ice=phy_ice(:,k)
     WHERE(fi_ice>=1.0  ) fi_ice=1.0
     WHERE(fi_ice<EPSFRA) fi_ice=0.0
     pctsrf_t(:,is_ter,k)=pctsrf(:,is_ter)       ! land soil
     pctsrf_t(:,is_lic,k)=pctsrf(:,is_lic)       ! land ice
     SELECT CASE(ix_sic)
        CASE(3)                                   ! SIC=pICE*(1-LIC-TER) (CPL)
        pctsrf_t(:,is_sic,k)=fi_ice(:)*(1.-pctsrf(:,is_lic)-pctsrf(:,is_ter))
        CASE(4)                                   ! SIC=pICE            (HIST)
        pctsrf_t(:,is_sic,k)=fi_ice(:)
        CASE DEFAULT                              ! SIC=pICE-LIC   (AMIP,ERAI)
        pctsrf_t(:,is_sic,k)=fi_ice-pctsrf_t(:,is_lic,k)
     END SELECT
     WHERE(pctsrf_t(:,is_sic,k)<=0) pctsrf_t(:,is_sic,k)=0.
     WHERE(1.0-zmasq<EPSFRA)
        pctsrf_t(:,is_sic,k)=0.0
        pctsrf_t(:,is_oce,k)=0.0
     ELSEWHERE
        WHERE(pctsrf_t(:,is_sic,k)>=1.0-zmasq)
           pctsrf_t(:,is_sic,k)=1.0-zmasq
           pctsrf_t(:,is_oce,k)=0.0
        ELSEWHERE
           pctsrf_t(:,is_oce,k)=1.0-zmasq-pctsrf_t(:,is_sic,k)
           WHERE(pctsrf_t(:,is_oce,k)<EPSFRA)
              pctsrf_t(:,is_oce,k)=0.0
              pctsrf_t(:,is_sic,k)=1.0-zmasq
           END WHERE
        END WHERE
     END WHERE
     nbad=COUNT(pctsrf_t(:,is_oce,k)<0.0)
     IF(nbad>0) WRITE(lunout,*) 'pb sous maille pour nb points = ',nbad
     nbad=COUNT(ABS(SUM(pctsrf_t(:,:,k),DIM=2)-1.0)>EPSFRA)
     IF(nbad>0) WRITE(lunout,*) 'pb sous surface pour nb points = ',nbad
  END DO
  DEALLOCATE(phy_ice)

!--- SST TREATMENT -------------------------------------------------------------
  CALL msg(0,""); CALL msg(0," *** TRAITEMENT DE LA SST ***")

! Input SST file selection
! Open file only to test if available
  DO ix_sst=1,SIZE(fsst)
     IF ( NF90_OPEN(TRIM(fsst(ix_sst)),NF90_NOWRITE,nid)==NF90_NOERR ) THEN
       sstfile=fsst(ix_sst); varname=vsst(ix_sst); EXIT
     END IF
  END DO
  IF(ix_sst==SIZE(fsst)+1) THEN
     WRITE(lunout,*) 'ERROR! No sst input file was found.'
     WRITE(lunout,*) 'One of following files must be available : '
     DO k=1,SIZE(fsst); WRITE(lunout,*) TRIM(fsst(k)); END DO
     CALL abort_physic('limit_netcdf','No sst file was found',1)
  END IF
  CALL ncerr(NF90_CLOSE(nid),sstfile)
  CALL msg(0,'Fichier choisi pour la temperature de mer: '//TRIM(sstfile))

  CALL get_2Dfield(sstfile,varname,'SST',ndays,phy_sst,flag=extrap)

!--- ALBEDO TREATMENT ----------------------------------------------------------
  CALL msg(0,""); CALL msg(0," *** TRAITEMENT DE L'ALBEDO ***")
  CALL get_2Dfield(falbe,valb,'ALB',ndays,phy_alb)

!--- REFERENCE GROUND HEAT FLUX TREATMENT --------------------------------------
  ALLOCATE(phy_bil(klon,ndays)); phy_bil=0.0

!--- OUTPUT FILE WRITING -------------------------------------------------------
  CALL msg(0,""); CALL msg(0,' *** Ecriture du fichier limit : debut ***')
  fnam="limit.nc"

  !--- File creation
  CALL ncerr(NF90_CREATE(fnam,NF90_CLOBBER,nid),fnam)
  CALL ncerr(NF90_PUT_ATT(nid,NF90_GLOBAL,"title","Fichier conditions aux limites"),fnam)
  str='File produced using ce0l executable.'
  str=TRIM(str)//NEW_LINE(' ')//'Sea Ice Concentration built from'
  SELECT CASE(ix_sic)
    CASE(1); str=TRIM(str)//' Amip mid-month boundary condition (BCS).'
    CASE(2); str=TRIM(str)//' Amip monthly mean observations.'
    CASE(3); str=TRIM(str)//' IPSL coupled model outputs.'
    CASE(4); str=TRIM(str)//' LMDZ model outputs.'
    CASE(5); str=TRIM(str)//' ci.nc file.'
  END SELECT
  str=TRIM(str)//NEW_LINE(' ')//'Sea Surface Temperature built from'
  SELECT CASE(ix_sst)
    CASE(1); str=TRIM(str)//' Amip mid-month boundary condition (BCS).'
    CASE(2); str=TRIM(str)//' Amip monthly mean observations.'
    CASE(3); str=TRIM(str)//' IPSL coupled model outputs.'
    CASE(4); str=TRIM(str)//' LMDZ model outputs.'
    CASE(5); str=TRIM(str)//' sstk.nc file.'
  END SELECT
  CALL ncerr(NF90_PUT_ATT(nid,NF90_GLOBAL,"history",TRIM(str)),fnam)

  !--- Dimensions creation
  CALL ncerr(NF90_DEF_DIM(nid,"points_physiques",klon,ndim),fnam)
  CALL ncerr(NF90_DEF_DIM(nid,"time",NF90_UNLIMITED,ntim),fnam)

  dims=[ndim,ntim]

  !--- Variables creation
  CALL ncerr(NF90_DEF_VAR(nid,"TEMPS",NF90_FORMAT,[ntim],id_tim),fnam)
  CALL ncerr(NF90_DEF_VAR(nid,"FOCE", NF90_FORMAT,dims,id_FOCE),fnam)
  CALL ncerr(NF90_DEF_VAR(nid,"FSIC", NF90_FORMAT,dims,id_FSIC),fnam)
  CALL ncerr(NF90_DEF_VAR(nid,"FTER", NF90_FORMAT,dims,id_FTER),fnam)
  CALL ncerr(NF90_DEF_VAR(nid,"FLIC", NF90_FORMAT,dims,id_FLIC),fnam)
  CALL ncerr(NF90_DEF_VAR(nid,"SST",  NF90_FORMAT,dims,id_SST),fnam)
  CALL ncerr(NF90_DEF_VAR(nid,"BILS", NF90_FORMAT,dims,id_BILS),fnam)
  CALL ncerr(NF90_DEF_VAR(nid,"ALB",  NF90_FORMAT,dims,id_ALB),fnam)
  CALL ncerr(NF90_DEF_VAR(nid,"RUG",  NF90_FORMAT,dims,id_RUG),fnam)
  call nf95_def_var(nid, "longitude", NF90_FLOAT, ndim, varid_longitude)
  call nf95_def_var(nid, "latitude",  NF90_FLOAT, ndim, varid_latitude)

  !--- Attributes creation
  CALL ncerr(NF90_PUT_ATT(nid,id_tim, "title","Jour dans l annee"),fnam)
  CALL ncerr(NF90_PUT_ATT(nid,id_tim, "calendar",calend),fnam)
  CALL ncerr(NF90_PUT_ATT(nid,id_FOCE,"title","Fraction ocean"),fnam)
  CALL ncerr(NF90_PUT_ATT(nid,id_FSIC,"title","Fraction glace de mer"),fnam)
  CALL ncerr(NF90_PUT_ATT(nid,id_FTER,"title","Fraction terre"),fnam)
  CALL ncerr(NF90_PUT_ATT(nid,id_FLIC,"title","Fraction land ice"),fnam)
  CALL ncerr(NF90_PUT_ATT(nid,id_SST ,"title","Temperature superficielle de la mer"),fnam)
  CALL ncerr(NF90_PUT_ATT(nid,id_BILS,"title","Reference flux de chaleur au sol"),fnam)
  CALL ncerr(NF90_PUT_ATT(nid,id_ALB, "title","Albedo a la surface"),fnam)
  CALL ncerr(NF90_PUT_ATT(nid,id_RUG, "title","Rugosite"),fnam)

  call nf95_put_att(nid, varid_longitude, "standard_name", "longitude")
  call nf95_put_att(nid, varid_longitude, "units", "degrees_east")

  call nf95_put_att(nid, varid_latitude, "standard_name", "latitude")
  call nf95_put_att(nid, varid_latitude, "units", "degrees_north")

  CALL ncerr(NF90_ENDDEF(nid),fnam)

  !--- Variables saving
  CALL ncerr(NF90_PUT_VAR(nid,id_tim,[(REAL(k),k=1,ndays)]),fnam)
  CALL ncerr(NF90_PUT_VAR(nid,id_FOCE,pctsrf_t(:,is_oce,:),[1,1],[klon,ndays]),fnam)
  CALL ncerr(NF90_PUT_VAR(nid,id_FSIC,pctsrf_t(:,is_sic,:),[1,1],[klon,ndays]),fnam)
  CALL ncerr(NF90_PUT_VAR(nid,id_FTER,pctsrf_t(:,is_ter,:),[1,1],[klon,ndays]),fnam)
  CALL ncerr(NF90_PUT_VAR(nid,id_FLIC,pctsrf_t(:,is_lic,:),[1,1],[klon,ndays]),fnam)
  CALL ncerr(NF90_PUT_VAR(nid,id_SST ,phy_sst(:,:),[1,1],[klon,ndays]),fnam)
  CALL ncerr(NF90_PUT_VAR(nid,id_BILS,phy_bil(:,:),[1,1],[klon,ndays]),fnam)
  CALL ncerr(NF90_PUT_VAR(nid,id_ALB ,phy_alb(:,:),[1,1],[klon,ndays]),fnam)
  CALL ncerr(NF90_PUT_VAR(nid,id_RUG ,phy_rug(:,:),[1,1],[klon,ndays]),fnam)
  call nf95_put_var(nid, varid_longitude, longitude_deg)
  call nf95_put_var(nid, varid_latitude, latitude_deg)

  CALL ncerr(NF90_CLOSE(nid),fnam)

  CALL msg(0,""); CALL msg(0,' *** Ecriture du fichier limit : fin ***')

  DEALLOCATE(pctsrf_t,phy_sst,phy_bil,phy_alb,phy_rug)


!===============================================================================
!
  CONTAINS
!
!===============================================================================


!-------------------------------------------------------------------------------
!
SUBROUTINE get_2Dfield(fnam, varname, mode, ndays, champo, flag, mask)
!
!-----------------------------------------------------------------------------
! Comments:
!   There are two assumptions concerning the NetCDF files, that are satisfied
!   with files that are conforming NC convention:
!     1) The last dimension of the variables used is the time record.
!     2) Dimensional variables have the same names as corresponding dimensions.
!-----------------------------------------------------------------------------
  USE netcdf, ONLY: NF90_OPEN, NF90_INQ_VARID, NF90_INQUIRE_VARIABLE, &
       NF90_CLOSE, NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, NF90_GET_VAR, &
       NF90_GET_ATT
  USE pchsp_95_m, only: pchsp_95
  USE pchfe_95_m, only: pchfe_95
  USE arth_m, only: arth
  USE indice_sol_mod

  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom2.h"
!-----------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*),  INTENT(IN)     :: fnam     ! NetCDF file name
  CHARACTER(LEN=*),  INTENT(IN)     :: varname  ! NetCDF variable name
  CHARACTER(LEN=*),  INTENT(IN)     :: mode     ! RUG, SIC, SST or ALB
  INTEGER,           INTENT(IN)     :: ndays    ! current year number of days
  REAL,    POINTER,  DIMENSION(:, :) :: champo  ! output field = f(t)
  LOGICAL, OPTIONAL, INTENT(IN)     :: flag     ! extrapol. (SST) old ice (SIC)
  REAL,    OPTIONAL, DIMENSION(iim, jjp1), INTENT(IN) :: mask
!------------------------------------------------------------------------------
! Local variables:
!--- NetCDF
  INTEGER           :: ncid, varid        ! NetCDF identifiers
  CHARACTER(LEN=ns) :: dnam               ! dimension name
!--- dimensions
  INTEGER           :: dids(4)            ! NetCDF dimensions identifiers
  REAL, ALLOCATABLE :: dlon_ini(:)        ! initial longitudes vector
  REAL, ALLOCATABLE :: dlat_ini(:)        ! initial latitudes  vector
  REAL, POINTER     :: dlon(:), dlat(:)   ! reordered lon/lat  vectors
!--- fields
  INTEGER :: imdep, jmdep, lmdep          ! dimensions of 'champ'
  REAL, ALLOCATABLE :: champ(:,:)         ! wanted field on initial grid
  REAL, ALLOCATABLE :: yder(:), timeyear(:)
  REAL              :: champint(iim,jjp1) ! interpolated field
  REAL, ALLOCATABLE :: champtime(:,:,:)
  REAL, ALLOCATABLE :: champan(:,:,:)
!--- input files
  CHARACTER(LEN=ns) :: fnam_m, fnam_p     ! previous/next files names
  CHARACTER(LEN=ns) :: cal_in             ! calendar
  CHARACTER(LEN=ns) :: units              ! attribute "units" in sic/sst file
  INTEGER           :: ndays_in           ! number of days
  REAL              :: value              ! mean/max value near equator
!--- misc
  INTEGER           :: i, j, k, l         ! loop counters
  REAL, ALLOCATABLE :: work(:,:)          ! used for extrapolation
  CHARACTER(LEN=ns) :: title, mess        ! for messages
  LOGICAL           :: is_bcs             ! flag for BCS data
  LOGICAL           :: extrp              ! flag for extrapolation
  LOGICAL           :: ll
  REAL              :: chmin, chmax, timeday, al
  INTEGER ierr, idx
  integer n_extrap ! number of extrapolated points
  logical skip

!------------------------------------------------------------------------------
!---Variables depending on keyword 'mode' -------------------------------------
  NULLIFY(champo)

  SELECT CASE(mode)
  CASE('RUG'); title='Rugosite'
  CASE('SIC'); title='Sea-ice'
  CASE('SST'); title='SST'
  CASE('ALB'); title='Albedo'
  END SELECT
  extrp=.FALSE.; IF(PRESENT(flag).AND.mode=='SST') extrp=flag
  is_bcs=(mode=='SIC'.AND.ix_sic==1).OR.(mode=='SST'.AND.ix_sst==1)
  idx=INDEX(fnam,'.nc')-1

!--- GETTING SOME DIMENSIONAL VARIABLES FROM FILE -----------------------------
  CALL msg(5,' Now reading file : '//TRIM(fnam))
  CALL ncerr(NF90_OPEN(fnam, NF90_NOWRITE, ncid),fnam)
  CALL ncerr(NF90_INQ_VARID(ncid, trim(varname), varid),fnam)
  CALL ncerr(NF90_INQUIRE_VARIABLE(ncid, varid, dimids=dids),fnam)

!--- Longitude
  CALL ncerr(NF90_INQUIRE_DIMENSION(ncid, dids(1), name=dnam, len=imdep),fnam)
  ALLOCATE(dlon_ini(imdep), dlon(imdep))
  CALL ncerr(NF90_INQ_VARID(ncid, dnam, varid), fnam)
  CALL ncerr(NF90_GET_VAR(ncid, varid, dlon_ini), fnam)
  CALL msg(5,'variable '//TRIM(dnam)//' dimension ', imdep)

!--- Latitude
  CALL ncerr(NF90_INQUIRE_DIMENSION(ncid, dids(2), name=dnam, len=jmdep),fnam)
  ALLOCATE(dlat_ini(jmdep), dlat(jmdep))
  CALL ncerr(NF90_INQ_VARID(ncid, dnam, varid), fnam)
  CALL ncerr(NF90_GET_VAR(ncid, varid, dlat_ini), fnam)
  CALL msg(5,'variable '//TRIM(dnam)//' dimension ', jmdep)

!--- Time (variable is not needed - it is rebuilt - but calendar is)
  CALL ncerr(NF90_INQUIRE_DIMENSION(ncid, dids(3), name=dnam, len=lmdep), fnam)
  ALLOCATE(timeyear(lmdep+2))
  CALL ncerr(NF90_INQ_VARID(ncid, dnam, varid), fnam)
  cal_in=' '
  IF(NF90_GET_ATT(ncid, varid, 'calendar', cal_in)/=NF90_NOERR) THEN
    SELECT CASE(mode)
      CASE('RUG', 'ALB'); cal_in='360d'
      CASE('SIC', 'SST'); cal_in='gregorian'
    END SELECT
    CALL msg(0,'WARNING: missing "calendar" attribute for "time" in '&
     &//TRIM(fnam)//'. Choosing default value.')
  END IF
  CALL strclean(cal_in)                     !--- REMOVE (WEIRD) NULL CHARACTERS
  CALL msg(0,'var, calendar, dim: '//TRIM(dnam)//' '//TRIM(cal_in), lmdep)
  
!--- Determining input file number of days, depending on calendar
  ndays_in=year_len(anneeref, cal_in)

!--- Rebuilding input time vector (field from input file might be unreliable)
  IF(lmdep==12) THEN
    timeyear=mid_month(anneeref, cal_in)
    CALL msg(0,'Monthly input file(s) for '//TRIM(title)//'.')
  ELSE IF(lmdep==ndays_in) THEN
    timeyear=[(REAL(k)-0.5,k=0,ndays_in+1)]
    CALL msg(0,'Daily input file (no time interpolation).')
  ELSE
    WRITE(mess,'(a,i3,a,i3,a)')'Mismatching input file: found',lmdep,        &
      ' records, 12/',ndays_in,' (monthly/daily needed).'
    CALL abort_physic('mid_month',TRIM(mess),1)
  END IF

!--- GETTING THE FIELD AND INTERPOLATING IT ----------------------------------
  ALLOCATE(champ(imdep, jmdep), champtime(iim, jjp1, lmdep+2))
  IF(extrp) ALLOCATE(work(imdep, jmdep))
  CALL msg(5,'')
  CALL msg(5,'READ AND INTERPOLATE HORIZONTALLY ', lmdep, ' FIELDS.')
  CALL ncerr(NF90_INQ_VARID(ncid, varname, varid), fnam)
  DO l=1, lmdep
    CALL ncerr(NF90_GET_VAR(ncid,varid,champ,[1,1,l],[imdep,jmdep,1]),fnam)
    CALL conf_dat2d(title, dlon_ini, dlat_ini, dlon, dlat, champ, .TRUE.)

    !--- FOR SIC/SST FIELDS ONLY
    IF(l==1.AND.is_in(mode,['SIC','SST'])) THEN

      !--- DETERMINE THE UNIT: READ FROM FILE OR ASSUMED USING FIELD VALUES
      ierr=NF90_GET_ATT(ncid, varid, 'units', units)
      IF(ierr==NF90_NOERR) THEN !--- ATTRIBUTE "units" FOUND IN THE FILE
        CALL strclean(units)
        IF(mode=='SIC'.AND.is_in(units,Perc)) units="%"
        IF(mode=='SIC'.AND.is_in(units,Frac)) units="1"
        IF(mode=='SST'.AND.is_in(units,DegC)) units="C"
        IF(mode=='SST'.AND.is_in(units,DegK)) units="K"
      ELSE                      !--- CHECK THE FIELD VALUES
        IF(mode=='SIC') value=MAXVAL(champ(:,:))
        IF(mode=='SST') value=   SUM(champ(:,jmdep/2),DIM=1)/REAL(imdep)
        IF(mode=='SIC') THEN; units="1"; IF(value>= 10.) units="%"; END IF
        IF(mode=='SST') THEN; units="C"; IF(value>=100.) units="K"; END IF
      END IF
      CALL msg(0,'INPUT FILE '//TRIM(title)//' UNIT IS: "'//TRIM(units)//'".')
      IF(ierr/=NF90_NOERR) CALL msg(0,'WARNING ! UNIT TO BE CHECKED ! '      &
        //'No "units" attribute, so only based on the fields values.')

      !--- CHECK VALUES ARE IN THE EXPECTED RANGE
      SELECT CASE(units)
        CASE('%'); ll=ANY(champ>100.0+EPSFRA); str='percentages > 100.'
        CASE('1'); ll=ANY(champ>  1.0+EPSFRA); str='fractions > 1.'
        CASE('C'); ll=ANY(champ<-100.).OR.ANY(champ> 60.); str='<-100 or >60 DegC'
        CASE('K'); ll=ANY(champ< 180.).OR.ANY(champ>330.); str='<180 or >330 DegK'
        CASE DEFAULT; CALL abort_physic(mode, 'Unrecognized '//TRIM(title)   &
                                                  //' unit: '//TRIM(units),1)
      END SELECT

      !--- DROPPED FOR BCS DATA (FRACTIONS CAN BE HIGHER THAN 1)
      IF(ll.AND.ix_sic/=1.AND.mode=='SIC') &
        CALL abort_physic(mode,'unrealistic '//TRIM(mode)//' found: '//TRIM(str))

    END IF

    IF(extrp) CALL extrapol(champ,imdep,jmdep,999999.,.TRUE.,.TRUE.,2,work)
    IF(l==1) THEN
      CALL msg(5,"--------------------------------------------------------")
      CALL msg(5,"$$$ Barycentric interpolation for "//TRIM(title)//" $$$")
      CALL msg(5,"--------------------------------------------------------")
    END IF
    IF(mode=='RUG') champ=LOG(champ)
    CALL inter_barxy(dlon,dlat(:jmdep-1),champ,rlonu(:iim),rlatv,champint)
    IF(mode=='RUG') THEN
      champint=EXP(champint)
      WHERE(NINT(mask)/=1) champint=0.001
    END IF
    champtime(:, :, l+1)=champint
  END DO
  CALL ncerr(NF90_CLOSE(ncid), fnam)

!--- FIRST RECORD: LAST ONE OF PREVIOUS YEAR (CURRENT YEAR IF UNAVAILABLE)
  fnam_m=fnam(1:idx)//'_m.nc'
  IF(NF90_OPEN(fnam_m,NF90_NOWRITE,ncid)==NF90_NOERR) THEN
    CALL msg(0,'Reading previous year file ("'//TRIM(fnam_m)//'") last record for '//TRIM(title))
    CALL ncerr(NF90_INQ_VARID(ncid, varname, varid),fnam_m)
    CALL ncerr(NF90_INQUIRE_VARIABLE(ncid, varid, dimids=dids),fnam_m)
    CALL ncerr(NF90_INQUIRE_DIMENSION(ncid, dids(3), len=l), fnam_m)
    CALL ncerr(NF90_GET_VAR(ncid,varid,champ,[1,1,l],[imdep,jmdep,1]),fnam_m)
    CALL ncerr(NF90_CLOSE(ncid), fnam_m)
    CALL conf_dat2d(title, dlon_ini, dlat_ini, dlon, dlat, champ, .TRUE.)
    IF(extrp) CALL extrapol(champ,imdep,jmdep,999999.,.TRUE.,.TRUE.,2,work)
    IF(mode=='RUG') champ=LOG(champ)
    CALL inter_barxy(dlon,dlat(:jmdep-1),champ,rlonu(:iim),rlatv,champint)
    IF(mode=='RUG') THEN
      champint=EXP(champint)
      WHERE(NINT(mask)/=1) champint=0.001
    END IF
    champtime(:, :, 1)=champint
  ELSE
    CALL msg(0,'Using current year file ("'//TRIM(fnam)//'") last record for '//TRIM(title))
    champtime(:, :, 1)=champtime(:, :, lmdep+1)
  END IF

!--- LAST RECORD: FIRST ONE OF NEXT YEAR (CURRENT YEAR IF UNAVAILABLE)
  fnam_p=fnam(1:idx)//'_p.nc'
  IF(NF90_OPEN(fnam_p,NF90_NOWRITE,ncid)==NF90_NOERR) THEN
    CALL msg(0,'Reading next year file ("'//TRIM(fnam_p)//'") first record for '//TRIM(title))
    CALL ncerr(NF90_INQ_VARID(ncid, varname, varid),fnam_p)
    CALL ncerr(NF90_GET_VAR(ncid,varid,champ,[1,1,1],[imdep,jmdep,1]),fnam_p)
    CALL ncerr(NF90_CLOSE(ncid), fnam_p)
    CALL conf_dat2d(title, dlon_ini, dlat_ini, dlon, dlat, champ, .TRUE.)
    IF(extrp) CALL extrapol(champ,imdep,jmdep,999999.,.TRUE.,.TRUE.,2,work)
    IF(mode=='RUG') champ=LOG(champ)
    CALL inter_barxy(dlon,dlat(:jmdep-1),champ,rlonu(:iim),rlatv,champint)
    IF(mode=='RUG') THEN
      champint=EXP(champint)
      WHERE(NINT(mask)/=1) champint=0.001
    END IF
    champtime(:, :, lmdep+2)=champint
  ELSE
    CALL msg(0,'Using current year file ("'//TRIM(fnam)//'") first record for '//TRIM(title))
    champtime(:, :, lmdep+2)=champtime(:, :, 2)
  END IF
  DEALLOCATE(dlon_ini, dlat_ini, dlon, dlat, champ)
  IF(extrp) DEALLOCATE(work)

!--- TIME INTERPOLATION ------------------------------------------------------
  IF(prt_level>0) THEN
     IF(ndays/=ndays_in) THEN
        WRITE(lunout,*)'DIFFERENT YEAR LENGTHS:'
        WRITE(lunout,*)' In the  input file: ',ndays_in
        WRITE(lunout,*)' In the output file: ',ndays
     END IF
     IF(lmdep==ndays_in) THEN
        WRITE(lunout, *)'NO TIME INTERPOLATION.'
        WRITE(lunout, *)' Daily input file.'
     ELSE
        IF(     is_bcs) WRITE(lunout, *)'LINEAR TIME INTERPOLATION.'
        IF(.NOT.is_bcs) WRITE(lunout, *)'SPLINES TIME INTERPOLATION.'
        WRITE(lunout, *)' Input time vector: ', timeyear
        WRITE(lunout, *)' Output time vector: from 0.5 to ', ndays-0.5
     END IF
  END IF
  ALLOCATE(champan(iip1, jjp1, ndays))

  IF(lmdep==ndays_in) THEN  !--- DAILY DATA: NO     TIME INTERPOLATION
    DO l=1,lmdep
      champan(1:iim,:,l)=champtime(:,:,l+1)
    END DO
  ELSE IF(is_bcs) THEN      !--- BCS   DATA: LINEAR TIME INTERPOLATION
    l=1
    DO k=1, ndays
      timeday = (REAL(k)-0.5)*REAL(ndays_in)/ndays
      IF(timeyear(l+1)<timeday) l=l+1
      al=(timeday-timeyear(l))/(timeyear(l+1)-timeyear(l))
      champan(1:iim,:,k) = champtime(1:iim,:,l)+al*(champtime(1:iim,:,l+1)-champtime(1:iim,:,l))
    END DO
  ELSE                      !--- AVE   DATA: SPLINE TIME INTERPOLATION
     skip = .false.
     n_extrap = 0
     ALLOCATE(yder(lmdep+2))
     DO j=1, jjp1
       DO i=1, iim
         yder = pchsp_95(timeyear, champtime(i, j, :), ibeg=2, iend=2, &
              vc_beg=0., vc_end=0.)
         CALL pchfe_95(timeyear, champtime(i, j, :), yder, skip, &
              arth(0.5, real(ndays_in) / ndays, ndays), champan(i, j, :), ierr)
         if (ierr < 0) stop 1
         n_extrap = n_extrap + ierr
       END DO
     END DO
     IF(n_extrap /= 0) WRITE(lunout,*) "get_2Dfield pchfe_95: n_extrap = ", n_extrap
     DEALLOCATE(yder)
  END IF
  champan(iip1, :, :)=champan(1, :, :)
  DEALLOCATE(champtime, timeyear)

!--- Checking the result
  DO j=1, jjp1
    CALL minmax(iip1, champan(1, j, 10), chmin, chmax)
    IF (prt_level>5) WRITE(lunout, *)' ',TRIM(title),' at time 10 ', chmin, chmax, j
  END DO

!--- SPECIAL FILTER FOR SST: SST>271.38 --------------------------------------
  IF(mode=='SST') THEN
    SELECT CASE(units)
      CASE("K"); CALL msg(0,'SST field is already in kelvins.')
      CASE("C"); CALL msg(0,'SST field converted from celcius degrees to kelvins.')
      champan(:, :, :)=champan(:, :, :)+273.15
    END SELECT
    CALL msg(0,'Filtering SST: Sea Surface Temperature >= 271.38')
    WHERE(champan<271.38) champan=271.38
  END IF

!--- SPECIAL FILTER FOR SIC: 0.0<SIC<1.0 -------------------------------------
  IF(mode=='SIC') THEN
    SELECT CASE(units)
      CASE("1"); CALL msg(0,'SIC field already in fraction of 1')
      CASE("%"); CALL msg(0,'SIC field converted from percentage to fraction of 1.')
       champan(:, :, :)=champan(:, :, :)/100.
    END SELECT
    CALL msg(0,'Filtering SIC: 0.0 <= Sea-ice <=1.0')
    WHERE(champan>1.0) champan=1.0
    WHERE(champan<0.0) champan=0.0
 END IF

!--- DYNAMICAL TO PHYSICAL GRID ----------------------------------------------
  ALLOCATE(champo(klon, ndays))
  DO k=1, ndays
    CALL gr_dyn_fi(1, iip1, jjp1, klon, champan(1, 1, k), champo(1, k))
  END DO
  DEALLOCATE(champan)

END SUBROUTINE get_2Dfield
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE start_init_orog0(lon_in,lat_in,phis,masque)
!
!-------------------------------------------------------------------------------
  USE grid_noro_m, ONLY: grid_noro0
  IMPLICIT NONE
!===============================================================================
! Purpose:  Compute "phis" just like it would be in start_init_orog.
!===============================================================================
! Arguments:
  REAL,             INTENT(IN)    :: lon_in(:), lat_in(:)   ! dim (iml) (jml)
  REAL,             INTENT(INOUT) :: phis(:,:), masque(:,:) ! dim (iml,jml)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=ns)  :: modname="start_init_orog0"
  INTEGER            :: fid, llm_tmp,ttm_tmp, iml,jml, iml_rel,jml_rel, itau(1)
  REAL               :: lev(1), date, dt, deg2rad
  REAL, ALLOCATABLE  :: lon_rad(:), lon_ini(:), lon_rel(:,:), relief_hi(:,:)
  REAL, ALLOCATABLE  :: lat_rad(:), lat_ini(:), lat_rel(:,:)
!-------------------------------------------------------------------------------
  iml=assert_eq(SIZE(lon_in),SIZE(phis,1),SIZE(masque,1),TRIM(modname)//" iml")
  jml=assert_eq(SIZE(lat_in),SIZE(phis,2),SIZE(masque,2),TRIM(modname)//" jml")
  IF(iml/=iip1) CALL abort_gcm(TRIM(modname),'iml/=iip1',1)
  IF(jml/=jjp1) CALL abort_gcm(TRIM(modname),'jml/=jjp1',1)
  pi=2.0*ASIN(1.0); deg2rad=pi/180.0
  IF(ANY(phis/=-99999.)) RETURN                  !--- phis ALREADY KNOWN

!--- HIGH RESOLUTION OROGRAPHY
  CALL flininfo(frelf, iml_rel, jml_rel, llm_tmp, ttm_tmp, fid)

  ALLOCATE(lat_rel(iml_rel,jml_rel),lon_rel(iml_rel,jml_rel))
  CALL flinopen(frelf, .FALSE., iml_rel, jml_rel, llm_tmp, lon_rel, lat_rel,   &
                lev, ttm_tmp, itau, date, dt, fid)
  ALLOCATE(relief_hi(iml_rel,jml_rel))
  CALL flinget(fid, vrel, iml_rel, jml_rel, llm_tmp, ttm_tmp, 1, 1, relief_hi)
  CALL flinclo(fid)

!--- IF ANGLES ARE IN DEGREES, THEY ARE CONVERTED INTO RADIANS
  ALLOCATE(lon_ini(iml_rel),lat_ini(jml_rel))
  lon_ini(:)=lon_rel(:,1); IF(MAXVAL(lon_rel)>pi) lon_ini=lon_ini*deg2rad
  lat_ini(:)=lat_rel(1,:); IF(MAXVAL(lat_rel)>pi) lat_ini=lat_ini*deg2rad

!--- FIELDS ARE PROCESSED TO BE ON STANDARD ANGULAR DOMAINS
  ALLOCATE(lon_rad(iml_rel),lat_rad(jml_rel))
  CALL conf_dat2d(vrel, lon_ini, lat_ini, lon_rad, lat_rad, relief_hi, .FALSE.)
  DEALLOCATE(lon_ini,lat_ini)

!--- COMPUTING SURFACE GEOPOTENTIAL USING ROUTINE grid_noro0
  WRITE(lunout,*)
  WRITE(lunout,*)'*** Compute surface geopotential ***'

!--- CALL OROGRAPHY MODULE (REDUCED VERSION) TO COMPUTE FIELDS
  CALL grid_noro0(lon_rad, lat_rad, relief_hi, lon_in, lat_in, phis, masque)
  phis = phis * 9.81
  phis(iml,:) = phis(1,:)
  DEALLOCATE(relief_hi,lon_rad,lat_rad)

END SUBROUTINE start_init_orog0
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE msg(lev,str1,i,str2)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,                    INTENT(IN) :: lev
  CHARACTER(LEN=*),           INTENT(IN) :: str1
  INTEGER,          OPTIONAL, INTENT(IN) :: i
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: str2
!-------------------------------------------------------------------------------
  IF(prt_level>=lev) THEN
    IF(PRESENT(str2)) THEN
      WRITE(lunout,*) TRIM(str1), i, TRIM(str2)
    ELSE IF(PRESENT(i)) THEN
      WRITE(lunout,*) TRIM(str1), i
    ELSE
      WRITE(lunout,*) TRIM(str1)
    END IF
  END IF

END SUBROUTINE msg
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE ncerr(ncres,fnam)
!
!-------------------------------------------------------------------------------
! Purpose: NetCDF errors handling.
!-------------------------------------------------------------------------------
  USE netcdf, ONLY : NF90_NOERR, NF90_STRERROR
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,          INTENT(IN) :: ncres
  CHARACTER(LEN=*), INTENT(IN) :: fnam
!-------------------------------------------------------------------------------
  IF(ncres/=NF90_NOERR) THEN
    WRITE(lunout,*)'Problem with file '//TRIM(fnam)//' in routine limit_netcdf.'
    CALL abort_physic('limit_netcdf',NF90_STRERROR(ncres),1)
  END IF

END SUBROUTINE ncerr
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE strclean(s)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Purpose: Remove tail null characters from the input string.
!-------------------------------------------------------------------------------
! Parameters:
  CHARACTER(LEN=*), INTENT(INOUT) :: s
!-------------------------------------------------------------------------------
! Local variable:
  INTEGER :: k
!-------------------------------------------------------------------------------
  k=LEN_TRIM(s); DO WHILE(ICHAR(s(k:k))==0); s(k:k)=' '; k=LEN_TRIM(s); END DO

END SUBROUTINE strclean
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
FUNCTION is_in(s1,s2) RESULT(res)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Purpose: Check wether s1 is present in the s2(:) list (case insensitive).
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: s1, s2(:)
  LOGICAL                      :: res
!-------------------------------------------------------------------------------
  res=.FALSE.; DO k=1,SIZE(s2); res=res.OR.strLow(s1)==strLow(s2(k)); END DO

END FUNCTION is_in
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
ELEMENTAL FUNCTION strLow(s) RESULT(res)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Purpose: Lower case conversion.
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: s
  CHARACTER(LEN=ns)            :: res
!-------------------------------------------------------------------------------
! Local variable:
  INTEGER :: k, ix
!-------------------------------------------------------------------------------
  res=s
  DO k=1,LEN(s); ix=IACHAR(s(k:k))
    IF(64<ix.AND.ix<91) res(k:k)=ACHAR(ix+32)
  END DO

END FUNCTION strLow
!
!-------------------------------------------------------------------------------

#endif
! of #ifndef CPP_1D
END SUBROUTINE limit_netcdf

END MODULE limit
!
!*******************************************************************************

