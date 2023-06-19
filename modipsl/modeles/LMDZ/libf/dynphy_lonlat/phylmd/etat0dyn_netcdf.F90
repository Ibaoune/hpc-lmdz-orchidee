MODULE etat0dyn
!
!*******************************************************************************
! Purpose: Create dynamical initial state using atmospheric fields from a
!          database of atmospheric to initialize the model.
!-------------------------------------------------------------------------------
! Comments:
!
!    *  This module is designed to work for Earth (and with ioipsl)
!
!    *  etat0dyn_netcdf routine can access to NetCDF data through the following
!  routine (to be called after restget):
!    CALL startget_dyn3d(varname, lon_in,  lat_in, pls, workvar,&
!                          champ, lon_in2, lat_in2)
!
!    *  Variables should have the following names in the NetCDF files:
!            'U'      : East ward wind              (in "ECDYN.nc")
!            'V'      : Northward wind              (in "ECDYN.nc")
!            'TEMP'   : Temperature                 (in "ECDYN.nc")
!            'R'      : Relative humidity           (in "ECDYN.nc")
!            'RELIEF' : High resolution orography   (in "Relief.nc") 
!
!    * The land mask and corresponding weights can be:
!      1) already known (in particular if etat0dyn has been called before) ;
!         in this case, ANY(masque(:,:)/=-99999.) = .TRUE.
!      2) computed using the ocean mask from the ocean model (to ensure ocean
!         fractions are the same for atmosphere and ocean) for coupled runs.
!         File name: "o2a.nc"  ;  variable name: "OceMask"
!      3) computed from topography file "Relief.nc" for forced runs.
!
!   *   There is a big mess with the longitude size. Should it be iml or iml+1 ?
!  I have chosen to use the iml+1 as an argument to this routine and we declare
!  internaly smaller fields when needed. This needs to be cleared once and for
!  all in LMDZ. A convention is required.
!-------------------------------------------------------------------------------
  USE ioipsl,         ONLY: flininfo, flinopen, flinget, flinclo, histclo
  USE assert_eq_m,    ONLY: assert_eq
  USE comconst_mod, ONLY: pi, cpp, kappa
  USE comvert_mod, ONLY: ap, bp, preff, pressure_exner
  USE temps_mod, ONLY: annee_ref, day_ref, itau_dyn, itau_phy, start_time
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: etat0dyn_netcdf

  include "iniprint.h"
  include "dimensions.h"
  include "paramet.h"
  include "comgeom2.h"
  include "comdissnew.h"
  REAL, SAVE :: deg2rad
  INTEGER,            SAVE      :: iml_dyn, jml_dyn, llm_dyn, ttm_dyn, fid_dyn
  REAL, ALLOCATABLE,  SAVE      :: lon_dyn(:,:), lat_dyn(:,:), levdyn_ini(:)
  CHARACTER(LEN=120), PARAMETER :: dynfname='ECDYN.nc'

CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE etat0dyn_netcdf(masque, phis)
!
!-------------------------------------------------------------------------------
! Purpose: Create dynamical initial states.
!-------------------------------------------------------------------------------
! Notes:  1) This routine is designed to work for Earth
!         2) If masque(:,:)/=-99999., masque and phis are already known.
!         Otherwise: compute it.
!-------------------------------------------------------------------------------
  USE control_mod
  USE regr_lat_time_coefoz_m, ONLY: regr_lat_time_coefoz
  USE regr_pr_o3_m,   ONLY: regr_pr_o3
  USE press_coefoz_m, ONLY: press_coefoz
  USE exner_hyb_m,    ONLY: exner_hyb
  USE exner_milieu_m, ONLY: exner_milieu
  USE infotrac,       ONLY: nqtot, tname
  USE filtreg_mod
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  REAL,    INTENT(INOUT) :: masque(iip1,jjp1)   !--- Land-ocean mask
  REAL,    INTENT(INOUT) :: phis  (iip1,jjp1)   !--- Ground geopotential
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname, fmt
  INTEGER            :: i, j, l, ji, itau, iday
  REAL               :: xpn, xps, time, phystep
  REAL, DIMENSION(iip1,jjp1)       :: psol
  REAL, DIMENSION(iip1,jjp1,llm+1) :: p3d
  REAL, DIMENSION(iip1,jjp1,llm)   :: uvent, t3d, tpot, qsat, qd
  REAL, DIMENSION(iip1,jjp1,llm)   :: pk, pls, y, masse
  REAL, DIMENSION(iip1,jjm ,llm)   :: vvent
  REAL, DIMENSION(ip1jm    ,llm)   :: pbarv
  REAL, DIMENSION(ip1jmp1  ,llm)   :: pbaru, phi, w
  REAL, DIMENSION(ip1jmp1)         :: pks
  REAL, DIMENSION(iim)             :: xppn, xpps
  REAL, ALLOCATABLE                :: q3d(:,:,:,:)
!-------------------------------------------------------------------------------
  modname='etat0dyn_netcdf'

  deg2rad = pi/180.0
  y(:,:,:)=0  !ym warning unitialized variable
  
! Compute psol AND tsol, knowing phis.
!*******************************************************************************
  CALL start_init_dyn(rlonv, rlatu, rlonu, rlatv, phis, psol)

! Mid-levels pressure computation
!*******************************************************************************
  CALL pression(ip1jmp1, ap, bp, psol, p3d)             !--- Update p3d
  IF(pressure_exner) THEN                               !--- Update pk, pks
    CALL exner_hyb   (ip1jmp1,psol,p3d,pks,pk)
  ELSE
    CALL exner_milieu(ip1jmp1,psol,p3d,pks,pk)
  END IF
  pls(:,:,:)=preff*(pk(:,:,:)/cpp)**(1./kappa)          !--- Update pls

! Update uvent, vvent, t3d and tpot
!*******************************************************************************
  uvent(:,:,:) = 0.0 ; vvent(:,:,:) = 0.0 ; t3d (:,:,:) = 0.0
  CALL startget_dyn3d('u'   ,rlonu,rlatu,pls,y ,uvent,rlonv,rlatv)
  CALL startget_dyn3d('v'   ,rlonv,rlatv,pls(:,:jjm,:),y(:,:jjm,:),vvent,      &
 &                           rlonu,rlatu(:jjm))
  CALL startget_dyn3d('t'   ,rlonv,rlatu,pls,y ,t3d ,rlonu,rlatv)
  tpot(:,:,:)=t3d(:,:,:)
  CALL startget_dyn3d('tpot',rlonv,rlatu,pls,pk,tpot,rlonu,rlatv)

  WRITE(lunout,*) 'T3D min,max:',MINVAL(t3d(:,:,:)),MAXVAL(t3d(:,:,:))
  WRITE(lunout,*) 'PLS min,max:',MINVAL(pls(:,:,:)),MAXVAL(pls(:,:,:))

! Humidity at saturation computation
!*******************************************************************************
  WRITE(lunout,*) 'avant q_sat'
  CALL q_sat(llm*jjp1*iip1, t3d, pls, qsat)
  WRITE(lunout,*) 'apres q_sat'
  WRITE(lunout,*) 'QSAT min,max:',MINVAL(qsat(:,:,:)),MAXVAL(qsat(:,:,:))
!  WRITE(lunout,*) 'QSAT :',qsat(10,20,:)
  qd (:,:,:) = 0.0
  CALL startget_dyn3d('q',rlonv,rlatu,pls,qsat,qd,rlonu,rlatv)
  ALLOCATE(q3d(iip1,jjp1,llm,nqtot)); q3d(:,:,:,:)=0.0 ; q3d(:,:,:,1)=qd(:,:,:)
  CALL flinclo(fid_dyn)

! Parameterization of ozone chemistry:
!*******************************************************************************
! Look for ozone tracer:
#ifndef INCA
  DO i=1,nqtot; IF(ANY(["O3","o3"]==tname(i))) EXIT; END DO
  IF(i/=nqtot+1) THEN
    CALL regr_lat_time_coefoz
    CALL press_coefoz
    CALL regr_pr_o3(p3d, q3d(:,:,:,i))
    q3d(:,:,:,i)=q3d(:,:,:,i)*48./ 29.                  !--- Mole->mass fraction         
  END IF
#endif
  q3d(iip1,:,:,:)=q3d(1,:,:,:)

! Writing
!*******************************************************************************
  CALL inidissip(lstardis, nitergdiv, nitergrot, niterh, tetagdiv, tetagrot,   &
       tetatemp, vert_prof_dissip)
  WRITE(lunout,*)'sortie inidissip'
  itau=0
  itau_dyn=0
  itau_phy=0
  iday=dayref+itau/day_step
  time=FLOAT(itau-(iday-dayref)*day_step)/day_step
  IF(time>1.) THEN
   time=time-1
   iday=iday+1
  END IF
  day_ref=dayref
  annee_ref=anneeref
  CALL geopot( ip1jmp1, tpot, pk, pks, phis, phi )
  WRITE(lunout,*)'sortie geopot'
  CALL caldyn0( itau, uvent, vvent, tpot, psol, masse, pk, phis,               &
                phi,  w, pbaru, pbarv, time+iday-dayref)
  WRITE(lunout,*)'sortie caldyn0'
  start_time = 0.
#ifdef CPP_PARA
  CALL dynredem0_loc( "start.nc", dayref, phis)
#else
  CALL dynredem0( "start.nc", dayref, phis)
#endif
  WRITE(lunout,*)'sortie dynredem0'
#ifdef CPP_PARA
  CALL dynredem1_loc( "start.nc", 0.0, vvent, uvent, tpot, q3d, masse, psol)
#else
  CALL dynredem1( "start.nc", 0.0, vvent, uvent, tpot, q3d, masse, psol)
#endif
  WRITE(lunout,*)'sortie dynredem1' 
  CALL histclo()

END SUBROUTINE etat0dyn_netcdf
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE startget_dyn3d(var, lon_in,  lat_in,  pls,  workvar,&
                        champ, lon_in2, lat_in2)
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!===============================================================================
! Purpose: Compute some quantities (u,v,t,q,tpot) using variables U,V,TEMP and R
!     (3D fields) of file dynfname.
!-------------------------------------------------------------------------------
! Note: An input auxilliary field "workvar" has to be specified in two cases:
!     * for "q":    the saturated humidity.
!     * for "tpot": the Exner function.
!===============================================================================
! Arguments:
  CHARACTER(LEN=*), INTENT(IN)    :: var
  REAL,             INTENT(IN)    :: lon_in(:)        ! dim (iml)
  REAL,             INTENT(IN)    :: lat_in(:)        ! dim (jml)
  REAL,             INTENT(IN)    :: pls    (:, :, :) ! dim (iml, jml, lml)
  REAL,             INTENT(IN)    :: workvar(:, :, :) ! dim (iml, jml, lml)
  REAL,             INTENT(INOUT) :: champ  (:, :, :) ! dim (iml, jml, lml)
  REAL,             INTENT(IN)    :: lon_in2(:)       ! dim (iml)
  REAL,             INTENT(IN)    :: lat_in2(:)       ! dim (jml2)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=10)  :: vname
  CHARACTER(LEN=256) :: msg, modname="startget_dyn3d"
  INTEGER            :: iml, jml, jml2, lml, il
  REAL               :: xppn, xpps
!-------------------------------------------------------------------------------
  iml=assert_eq([SIZE(lon_in),SIZE(pls,1),SIZE(workvar,1),SIZE(champ,1),       &
     &                                    SIZE(lon_in2)], TRIM(modname)//" iml")
  jml=assert_eq( SIZE(lat_in),SIZE(pls,2),SIZE(workvar,2),SIZE(champ,2),       &
     &                                                    TRIM(modname)//" jml")
  lml=assert_eq(              SIZE(pls,3),SIZE(workvar,3),SIZE(champ,3),       &
     &                                                    TRIM(modname)//" lml")
  jml2=SIZE(lat_in2)

!--- CHECK IF THE FIELD IS KNOWN
   SELECT CASE(var)
    CASE('u');    vname='U'
    CASE('v');    vname='V'
    CASE('t');    vname='TEMP'
    CASE('q');    vname='R';    msg='humidity as the saturated humidity'
    CASE('tpot'); msg='potential temperature as the Exner function'
    CASE DEFAULT; msg='No rule to extract variable '//TRIM(var)
      CALL abort_gcm(modname,TRIM(msg)//' from any data set',1)
  END SELECT

!--- CHECK IF SOMETHING IS MISSING
  IF((var=='tpot'.OR.var=='q').AND.MINVAL(workvar)==MAXVAL(workvar)) THEN
    msg='Could not compute '//TRIM(msg)//' is missing or constant.'
    CALL abort_gcm(modname,TRIM(msg),1)
  END IF

!--- INTERPOLATE 3D FIELD IF NEEDED
  IF(var/='tpot') CALL start_inter_3d(TRIM(vname),lon_in,lat_in,lon_in2,      &
                                                  lat_in2,pls,champ)

!--- COMPUTE THE REQUIRED FILED
  SELECT CASE(var)
    CASE('u'); DO il=1,lml; champ(:,:,il)=champ(:,:,il)*cu(:,1:jml); END DO
      champ(iml,:,:)=champ(1,:,:)                   !--- Eastward wind

    CASE('v'); DO il=1,lml; champ(:,:,il)=champ(:,:,il)*cv(:,1:jml); END DO
      champ(iml,:,:)=champ(1,:,:)                   !--- Northward wind

    CASE('tpot','q')
      IF(var=='tpot') THEN; champ=champ*cpp/workvar !--- Potential temperature
      ELSE;                 champ=champ*.01*workvar !--- Relative humidity
        WHERE(champ<0.) champ=1.0E-10
      END IF
      DO il=1,lml
        xppn = SUM(aire(:,1  )*champ(:,1  ,il))/apoln
        xpps = SUM(aire(:,jml)*champ(:,jml,il))/apols
        champ(:,1  ,il) = xppn
        champ(:,jml,il) = xpps
      END DO
  END SELECT

END SUBROUTINE startget_dyn3d
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE start_init_dyn(lon_in,lat_in,lon_in2,lat_in2,zs,psol)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!===============================================================================
! Purpose:   Compute psol, knowing phis.
!===============================================================================
! Arguments:
  REAL,    INTENT(IN)  :: lon_in (:),  lat_in (:)    ! dim (iml) (jml)
  REAL,    INTENT(IN)  :: lon_in2(:),  lat_in2(:)    ! dim (iml) (jml2)
  REAL,    INTENT(IN)  :: zs  (:,:)                  ! dim (iml,jml)
  REAL,    INTENT(OUT) :: psol(:,:)                  ! dim (iml,jml)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname='start_init_dyn'
  REAL               :: date, dt
  INTEGER            :: iml, jml, jml2, itau(1)
  REAL, ALLOCATABLE  :: lon_rad(:), lon_ini(:), var_ana(:,:)
  REAL, ALLOCATABLE  :: lat_rad(:), lat_ini(:)
  REAL, ALLOCATABLE  :: z(:,:), ps(:,:), ts(:,:)
!-------------------------------------------------------------------------------
  iml=assert_eq(SIZE(lon_in),SIZE(zs,1),SIZE(psol,1),SIZE(lon_in2),            &
      &                                              TRIM(modname)//" iml")
  jml=assert_eq(SIZE(lat_in),SIZE(zs,2),SIZE(psol,2),TRIM(modname)//" jml")
  jml2=SIZE(lat_in2)

  WRITE(lunout,*) 'Opening the surface analysis'
  CALL flininfo(dynfname, iml_dyn, jml_dyn, llm_dyn, ttm_dyn, fid_dyn)
  WRITE(lunout,*) 'Values read: ', iml_dyn, jml_dyn, llm_dyn, ttm_dyn

  ALLOCATE(lon_dyn(iml_dyn,jml_dyn), lat_dyn(iml_dyn,jml_dyn))
  ALLOCATE(levdyn_ini(llm_dyn))
  CALL flinopen(dynfname, .FALSE., iml_dyn, jml_dyn, llm_dyn,                  &
                lon_dyn,lat_dyn,levdyn_ini,ttm_dyn,itau,date,dt,fid_dyn)

!--- IF ANGLES ARE IN DEGREES, THEY ARE CONVERTED INTO RADIANS
  ALLOCATE(lon_ini(iml_dyn),lat_ini(jml_dyn))
  lon_ini(:)=lon_dyn(:,1); IF(MAXVAL(lon_dyn)>pi) lon_ini=lon_ini*deg2rad
  lat_ini(:)=lat_dyn(1,:); IF(MAXVAL(lat_dyn)>pi) lat_ini=lat_ini*deg2rad

  ALLOCATE(var_ana(iml_dyn,jml_dyn),lon_rad(iml_dyn),lat_rad(jml_dyn))
  CALL get_var_dyn('Z',z)                        !--- SURFACE GEOPOTENTIAL
  CALL get_var_dyn('SP',ps)                      !--- SURFACE PRESSURE
  CALL get_var_dyn('ST',ts)                      !--- SURFACE TEMPERATURE
!  CALL flinclo(fid_dyn)
  DEALLOCATE(var_ana,lon_rad,lat_rad,lon_ini,lat_ini)

!--- PSOL IS COMPUTED IN PASCALS
  psol(:iml-1,:) = ps(:iml-1,:)*(1.0+(z(:iml-1,:)-zs(:iml-1,:))/287.0          &
     &            /ts(:iml-1,:))
  psol(iml,:)=psol(1,:)
  DEALLOCATE(z,ps,ts)
  psol(:,1  )=SUM(aire(1:iml-1,1  )*psol(1:iml-1,1  ))/apoln  !--- NORTH POLE
  psol(:,jml)=SUM(aire(1:iml-1,jml)*psol(1:iml-1,jml))/apols  !--- SOUTH POLE

CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE get_var_dyn(title,field)
!
!-------------------------------------------------------------------------------
  USE conf_dat_m, ONLY: conf_dat2d
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*),  INTENT(IN)    :: title
  REAL, ALLOCATABLE, INTENT(INOUT) :: field(:,:)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: msg
  INTEGER :: tllm
!-------------------------------------------------------------------------------
  SELECT CASE(title)
    CASE('Z');     tllm=0;       msg='geopotential'
    CASE('SP');    tllm=0;       msg='surface pressure'
    CASE('ST');    tllm=llm_dyn; msg='temperature'
  END SELECT
  IF(.NOT.ALLOCATED(field)) THEN
    ALLOCATE(field(iml,jml))
    CALL flinget(fid_dyn, title, iml_dyn,jml_dyn, tllm, ttm_dyn, 1, 1, var_ana)
    CALL conf_dat2d(title, lon_ini, lat_ini, lon_rad, lat_rad, var_ana, .TRUE.)
    CALL interp_startvar(title, .TRUE., lon_rad,lat_rad, var_ana,              &
                                        lon_in, lat_in, lon_in2, lat_in2, field)
  ELSE IF(SIZE(field)/=SIZE(z)) THEN
    msg='The '//TRIM(msg)//' field we have does not have the right size'
    CALL abort_gcm(TRIM(modname),msg,1)
  END IF

END SUBROUTINE get_var_dyn
!
!-------------------------------------------------------------------------------

END SUBROUTINE start_init_dyn
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE start_inter_3d(var,lon_in,lat_in,lon_in2,lat_in2,pls_in,var3d)
!
!-------------------------------------------------------------------------------
  USE conf_dat_m, ONLY: conf_dat3d
  USE pchsp_95_m, ONLY: pchsp_95
  USE pchfe_95_m, ONLY: pchfe_95
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: var
  REAL,    INTENT(IN)  :: lon_in(:),  lat_in(:)   ! dim (iml) (jml)
  REAL,    INTENT(IN)  :: lon_in2(:), lat_in2(:)  ! dim (iml) (jml2)
  REAL,    INTENT(IN)  :: pls_in(:,:,:)           ! dim (iml,jml,lml)
  REAL,    INTENT(OUT) :: var3d (:,:,:)           ! dim (iml,jml,lml)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname='start_inter_3d'
  LOGICAL :: skip
  REAL    :: chmin, chmax
  INTEGER :: iml, jml, lml, jml2, ii, ij, il, ierr
  INTEGER :: n_extrap                             ! Extrapolated points number
  REAL, ALLOCATABLE :: ax(:), lon_rad(:), lon_ini(:), lev_dyn(:), yder(:)
  REAL, ALLOCATABLE :: ay(:), lat_rad(:), lat_ini(:), var_tmp3d(:,:,:)
  REAL, ALLOCATABLE, SAVE :: var_ana3d(:,:,:)
!-------------------------------------------------------------------------------
  iml=assert_eq(SIZE(lon_in),SIZE(lon_in2),SIZE(pls_in,1),SIZE(var3d,1),TRIM(modname)//" iml")
  jml=assert_eq(SIZE(lat_in),              SIZE(pls_in,2),SIZE(var3d,2),TRIM(modname)//" jml")
  lml=assert_eq(SIZE(pls_in,3),SIZE(var3d,3),TRIM(modname)//" lml"); jml2=SIZE(lat_in2)

  WRITE(lunout, *)'Going into flinget to extract the 3D field.'
  IF(.NOT.ALLOCATED(var_ana3d)) ALLOCATE(var_ana3d(iml_dyn, jml_dyn, llm_dyn))
  CALL flinget(fid_dyn,var,iml_dyn,jml_dyn,llm_dyn,ttm_dyn,1,1,var_ana3d)

!--- ANGLES IN DEGREES ARE CONVERTED INTO RADIANS
  ALLOCATE(lon_ini(iml_dyn), lat_ini(jml_dyn))
  lon_ini(:)=lon_dyn(:,1); IF(MAXVAL(lon_dyn)>pi) lon_ini=lon_ini*deg2rad
  lat_ini(:)=lat_dyn(1,:); IF(MAXVAL(lat_dyn)>pi) lat_ini=lat_ini*deg2rad

!--- FIELDS ARE PROCESSED TO BE ON STANDARD ANGULAR DOMAINS
  ALLOCATE(lon_rad(iml_dyn), lat_rad(jml_dyn), lev_dyn(llm_dyn))
  CALL conf_dat3d(var, lon_ini, lat_ini, levdyn_ini,                           &
                       lon_rad, lat_rad, lev_dyn, var_ana3d, .TRUE.)
  DEALLOCATE(lon_ini, lat_ini)

!--- COMPUTE THE REQUIRED FIELDS USING ROUTINE grid_noro
  ALLOCATE(var_tmp3d(iml,jml,llm_dyn))
  DO il = 1,llm_dyn
    CALL interp_startvar(var,il==1,lon_rad,lat_rad,var_ana3d(:,:,il),          &
                          lon_in,lat_in,lon_in2,lat_in2,var_tmp3d(:,:,il))
  END DO
  DEALLOCATE(lon_rad, lat_rad)

!--- VERTICAL INTERPOLATION FROM TOP OF ATMOSPHERE TO GROUND
  ALLOCATE(ax(llm_dyn),ay(llm_dyn),yder(llm_dyn))
  ax = lev_dyn(llm_dyn:1:-1) 
  skip = .FALSE.
  n_extrap = 0
  DO ij=1, jml
    DO ii=1, iml-1
      ay = var_tmp3d(ii, ij, llm_dyn:1:-1)
      yder = pchsp_95(ax, ay, ibeg=2, iend=2, vc_beg=0., vc_end=0.)
      CALL pchfe_95(ax, ay, yder, skip, pls_in(ii, ij, lml:1:-1),              &
           var3d(ii, ij, lml:1:-1), ierr)
      IF(ierr<0) CALL abort_gcm(TRIM(modname),'error in pchfe_95',1)
      n_extrap = n_extrap + ierr
    END DO
  END DO
  IF(n_extrap/=0) WRITE(lunout,*)TRIM(modname)//" pchfe_95: n_extrap=", n_extrap
  var3d(iml, :, :) = var3d(1, :, :) 

  DO il=1, lml
    CALL minmax(iml*jml, var3d(1, 1, il), chmin, chmax)
    WRITE(lunout, *)' '//TRIM(var)//'  min max l ', il, chmin, chmax
  END DO

END SUBROUTINE start_inter_3d
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE interp_startvar(nam,ibeg,lon,lat,vari,lon1,lat1,lon2,lat2,varo)
!
!-------------------------------------------------------------------------------
  USE inter_barxy_m, ONLY: inter_barxy
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*), INTENT(IN)  :: nam
  LOGICAL,          INTENT(IN)  :: ibeg
  REAL,             INTENT(IN)  :: lon(:), lat(:)   ! dim (ii) (jj)
  REAL,             INTENT(IN)  :: vari(:,:)        ! dim (ii,jj)
  REAL,             INTENT(IN)  :: lon1(:), lat1(:) ! dim (i1) (j1)
  REAL,             INTENT(IN)  :: lon2(:), lat2(:) ! dim (i1) (j2)
  REAL,             INTENT(OUT) :: varo(:,:)        ! dim (i1) (j1)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname="interp_startvar"
  INTEGER            :: ii, jj, i1, j1, j2
  REAL, ALLOCATABLE  :: vtmp(:,:)
!-------------------------------------------------------------------------------
  ii=assert_eq(SIZE(lon),            SIZE(vari,1),TRIM(modname)//" ii")
  jj=assert_eq(SIZE(lat),            SIZE(vari,2),TRIM(modname)//" jj")
  i1=assert_eq(SIZE(lon1),SIZE(lon2),SIZE(varo,1),TRIM(modname)//" i1")
  j1=assert_eq(SIZE(lat1),           SIZE(varo,2),TRIM(modname)//" j1")
  j2=SIZE(lat2)
  ALLOCATE(vtmp(i1-1,j1))
  IF(ibeg.AND.prt_level>1) THEN
    WRITE(lunout,*)"---------------------------------------------------------"
    WRITE(lunout,*)"$$$ Interpolation barycentrique pour "//TRIM(nam)//" $$$"
    WRITE(lunout,*)"---------------------------------------------------------"
  END IF
  CALL inter_barxy(lon, lat(:jj-1), vari, lon2(:i1-1), lat2, vtmp)
  CALL gr_int_dyn(vtmp, varo, i1-1, j1)

END SUBROUTINE interp_startvar
!
!-------------------------------------------------------------------------------

END MODULE etat0dyn
!
!*******************************************************************************
