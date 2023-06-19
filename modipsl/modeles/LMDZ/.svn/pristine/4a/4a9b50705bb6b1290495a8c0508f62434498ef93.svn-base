!
! $Id$
!
MODULE etat0phys
!
!*******************************************************************************
! Purpose: Create physical initial state using atmospheric fields from a
!          database of atmospheric to initialize the model.
!-------------------------------------------------------------------------------
! Comments:
!
!    *  This module is designed to work for Earth (and with ioipsl)
!
!    *  etat0phys_netcdf routine can access to NetCDF data through subroutines:
!         "start_init_phys" for variables contained in file "ECPHY.nc":
!            'ST'     : Surface temperature
!            'CDSW'   : Soil moisture
!         "start_init_orog" for variables contained in file "Relief.nc":
!            'RELIEF' : High resolution orography 
!
!    * The land mask and corresponding weights can be:
!      1) computed using the ocean mask from the ocean model (to ensure ocean
!         fractions are the same for atmosphere and ocean) for coupled runs.
!         File name: "o2a.nc"  ;  variable name: "OceMask"
!      2) computed from topography file "Relief.nc" for forced runs.
!
!    * Allowed values for read_climoz flag are 0, 1 and 2:
!      0: do not read an ozone climatology
!      1: read a single ozone climatology that will be used day and night
!      2: read two ozone climatologies, the average day and night climatology
!         and the daylight climatology
!-------------------------------------------------------------------------------
!    * There is a big mess with the longitude size. Should it be iml or iml+1 ?
!  I have chosen to use the iml+1 as an argument to this routine and we declare
!  internaly smaller fields when needed. This needs to be cleared once and for
!  all in LMDZ. A convention is required.
!-------------------------------------------------------------------------------

  USE ioipsl,             ONLY: flininfo, flinopen, flinget, flinclo
  USE assert_eq_m,        ONLY: assert_eq
  USE dimphy,             ONLY: klon
  USE conf_dat_m,         ONLY: conf_dat2d
  USE phys_state_var_mod, ONLY: zmea, zstd, zsig, zgam, zthe, zpic, zval, z0m, &
          solsw, radsol, t_ancien, wake_deltat, wake_s,  rain_fall, qsol, z0h, &
          sollw,sollwdown, rugoro, q_ancien, wake_deltaq, wake_pe, snow_fall, ratqs,w01, &
    sig1, ftsol, clwcon, fm_therm, wake_Cstar,  pctsrf,  entr_therm,radpas, f0,&
    zmax0,fevap, rnebcon,falb_dir, falb_dif, wake_fip,    agesno,  detr_therm, pbl_tke,  &
    phys_state_var_init, ql_ancien, qs_ancien, prlw_ancien, prsw_ancien, &
    prw_ancien, u10m,v10m, treedrg, u_ancien, v_ancien, wake_delta_pbl_TKE, wake_dens, &
    ale_bl, ale_bl_trig, alp_bl, &
    ale_wake, ale_bl_stat

  USE comconst_mod, ONLY: pi, dtvr

  PRIVATE
  PUBLIC :: etat0phys_netcdf

  include "iniprint.h"
  include "dimensions.h"
  include "paramet.h"
  include "comgeom2.h"
  include "dimsoil.h"
  include "clesphys.h"
  REAL, SAVE :: deg2rad
  REAL, SAVE, ALLOCATABLE :: tsol(:)
  INTEGER,            SAVE      :: iml_phys, jml_phys, llm_phys, ttm_phys, fid_phys
  REAL, ALLOCATABLE,  SAVE      :: lon_phys(:,:), lat_phys(:,:), levphys_ini(:)
  CHARACTER(LEN=256), PARAMETER :: oroparam="oro_params.nc"
  CHARACTER(LEN=256), PARAMETER :: orofname="Relief.nc", orogvar="RELIEF"
  CHARACTER(LEN=256), PARAMETER :: phyfname="ECPHY.nc",  psrfvar="SP"
  CHARACTER(LEN=256), PARAMETER :: qsolvar="CDSW",       tsrfvar="ST"


CONTAINS


!-------------------------------------------------------------------------------
!
SUBROUTINE etat0phys_netcdf(masque, phis)
!
!-------------------------------------------------------------------------------
! Purpose: Creates initial states
!-------------------------------------------------------------------------------
! Notes:  1) This routine is designed to work for Earth
!         2) If masque(:,:)/=-99999., masque and phis are already known.
!         Otherwise: compute it.
!-------------------------------------------------------------------------------
  USE control_mod
  USE fonte_neige_mod
  USE pbl_surface_mod
  USE regr_horiz_time_climoz_m, ONLY: regr_horiz_time_climoz
  USE indice_sol_mod
  USE conf_phys_m, ONLY: conf_phys
  USE init_ssrf_m, ONLY: start_init_subsurf
  !use ioipsl_getincom
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  REAL,    INTENT(INOUT) :: masque(:,:) !--- Land mask           dim(iip1,jjp1)
  REAL,    INTENT(INOUT) :: phis  (:,:) !--- Ground geopotential dim(iip1,jjp1)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname="etat0phys_netcdf", fmt
  INTEGER            :: i, j, l, ji, iml, jml
  LOGICAL            :: read_mask
  REAL               :: phystep, dummy
  REAL, DIMENSION(SIZE(masque,1),SIZE(masque,2)) :: masque_tmp,phiso
  REAL, DIMENSION(klon)               :: sn, rugmer, run_off_lic_0, fder
  REAL, DIMENSION(klon,nbsrf)         :: qsolsrf, snsrf
  REAL, DIMENSION(klon,nsoilmx,nbsrf) :: tsoil

!--- Arguments for conf_phys
  LOGICAL :: ok_journe, ok_mensuel, ok_instan, ok_hf, ok_LES, callstats
  REAL    :: solarlong0, seuil_inversion, fact_cldcon, facttemps
  LOGICAL :: ok_newmicro
  INTEGER :: iflag_radia, iflag_cldcon, iflag_ratqs
  REAL    :: ratqsbas, ratqshaut, tau_ratqs
  LOGICAL :: ok_ade, ok_aie, ok_volcan, ok_alw, ok_cdnc, aerosol_couple, chemistry_couple
  INTEGER :: flag_aerosol
  INTEGER :: flag_aerosol_strat
  LOGICAL :: flag_aer_feedback
  LOGICAL :: flag_bc_internal_mixture
  REAL    :: bl95_b0, bl95_b1
  INTEGER :: read_climoz                        !--- Read ozone climatology
  REAL    :: alp_offset
  LOGICAL :: filtre_oro=.false.

  deg2rad= pi/180.0
  iml=assert_eq(SIZE(masque,1),SIZE(phis,1),TRIM(modname)//" iml")
  jml=assert_eq(SIZE(masque,2),SIZE(phis,2),TRIM(modname)//" jml")

! Physics configuration
!*******************************************************************************
  CALL conf_phys(  ok_journe, ok_mensuel, ok_instan, ok_hf, ok_LES,     &
                   callstats,                                           &
                   solarlong0,seuil_inversion,                          &
                   fact_cldcon, facttemps,ok_newmicro,iflag_radia,      &
                   iflag_cldcon,                                        &
                   iflag_ratqs,ratqsbas,ratqshaut,tau_ratqs,            &
                   ok_ade, ok_aie, ok_alw, ok_cdnc, ok_volcan,          &
                   aerosol_couple,                                      &
                   chemistry_couple, flag_aerosol, flag_aerosol_strat,  & 
                   flag_aer_feedback,                                   &
                   flag_bc_internal_mixture, bl95_b0, bl95_b1, &
                   read_climoz, alp_offset)
  CALL phys_state_var_init(read_climoz)

!--- Initial atmospheric CO2 conc. from .def file
  co2_ppm0 = co2_ppm

! Compute ground geopotential, sub-cells quantities and possibly the mask.
!*******************************************************************************
  read_mask=ANY(masque/=-99999.); masque_tmp=masque
  CALL start_init_orog(rlonv, rlatu, phis, masque_tmp)

  CALL getin('filtre_oro',filtre_oro)
  IF (filtre_oro) CALL filtreoro(size(phis,1),size(phis,2),phis,masque_tmp,rlatu)

  WRITE(fmt,"(i4,'i1)')")iml ; fmt='('//ADJUSTL(fmt)
  IF(.NOT.read_mask) THEN                       !--- Keep mask form orography
    masque=masque_tmp
    IF(prt_level>=1) THEN
      WRITE(lunout,*)'BUILT MASK :'
      WRITE(lunout,fmt) NINT(masque)
    END IF
    WHERE(   masque(:,:)<EPSFRA) masque(:,:)=0.
    WHERE(1.-masque(:,:)<EPSFRA) masque(:,:)=1.
  END IF
  CALL gr_dyn_fi(1,iml,jml,klon,masque,zmasq) !--- Land mask to physical grid

! Compute tsol and qsol on physical grid, knowing phis on 2D grid.
!*******************************************************************************
  CALL start_init_phys(rlonu, rlatv, phis)

! Some initializations.
!*******************************************************************************
  sn    (:) = 0.0                               !--- Snow
  radsol(:) = 0.0                               !--- Net radiation at ground
  rugmer(:) = 0.001                             !--- Ocean rugosity
  !--- Ozone (zonal or 3D) interpolation in space and time (if 2nd arg is TRUE)
  IF(read_climoz>=1) CALL regr_horiz_time_climoz(read_climoz,ok_daily_climoz)

! Sub-surfaces initialization.
!*******************************************************************************
  CALL start_init_subsurf(read_mask)

! Write physical initial state
!*******************************************************************************
  WRITE(lunout,*)'phystep ',dtvr,iphysiq,nbapp_rad
  phystep = dtvr * FLOAT(iphysiq)
  radpas  = NINT (86400./phystep/ FLOAT(nbapp_rad) )
  WRITE(lunout,*)'phystep =', phystep, radpas

! Init: ftsol, snsrf, qsolsrf, tsoil, rain_fall, snow_fall, solsw, sollw, z0
!*******************************************************************************
  DO i=1,nbsrf; ftsol(:,i) = tsol; END DO
  DO i=1,nbsrf; snsrf(:,i) = sn;   END DO
  falb_dir(:, :, is_ter) = 0.08
  falb_dir(:, :, is_lic) = 0.6
  falb_dir(:, :, is_oce) = 0.5
  falb_dir(:, :, is_sic) = 0.6

!ym warning missing init for falb_dif => set to 0
  falb_dif(:,:,:)=0

  u10m(:,:)=0  
  v10m(:,:)=0  
  treedrg(:,:,:)=0

  fevap(:,:) = 0.
  DO i=1,nbsrf; qsolsrf(:,i)=150.; END DO
  DO i=1,nbsrf; DO j=1,nsoilmx; tsoil(:,j,i) = tsol; END DO; END DO
  rain_fall  = 0.
  snow_fall  = 0.
  solsw      = 165.
  sollw      = -53.
!ym warning missing init for sollwdown => set to 0
  sollwdown  = 0.
  t_ancien   = 273.15
  q_ancien   = 0.
  ql_ancien = 0.
  qs_ancien = 0.
  prlw_ancien = 0.
  prsw_ancien = 0.
  prw_ancien = 0.
  agesno     = 0.
  
  u_ancien = 0.
  v_ancien = 0.
  wake_delta_pbl_TKE(:,:,:)=0
  wake_dens(:)=0
  ale_bl = 0.
  ale_bl_trig =0.
  alp_bl=0.
  ale_wake=0.
  ale_bl_stat=0.
  
  z0m(:,:)=0 ! ym missing 5th subsurface initialization
 
  z0m(:,is_oce) = rugmer(:)
  z0m(:,is_ter) = MAX(1.0e-05,zstd(:)*zsig(:)/2.0)
  z0m(:,is_lic) = MAX(1.0e-05,zstd(:)*zsig(:)/2.0)
  z0m(:,is_sic) = 0.001
  z0h(:,:)=z0m(:,:)

  fder    = 0.0
  clwcon  = 0.0
  rnebcon = 0.0
  ratqs   = 0.0
  run_off_lic_0 = 0.0 
  rugoro  = 0.0

! Before phyredem calling, surface modules and values to be saved in startphy.nc
! are initialized
!*******************************************************************************
  dummy            = 1.0
  pbl_tke(:,:,:)   = 1.e-8 
  zmax0(:)         = 40.
  f0(:)            = 1.e-5
  sig1(:,:)        = 0.
  w01(:,:)         = 0.
  wake_deltat(:,:) = 0.
  wake_deltaq(:,:) = 0.
  wake_s(:)        = 0.
  wake_cstar(:)    = 0.
  wake_fip(:)      = 0.
  wake_pe          = 0.
  fm_therm         = 0.
  entr_therm       = 0.
  detr_therm       = 0.

  CALL fonte_neige_init(run_off_lic_0)
  CALL pbl_surface_init( fder, snsrf, qsolsrf, tsoil )
  CALL phyredem( "startphy.nc" )

!  WRITE(lunout,*)'CCCCCCCCCCCCCCCCCC REACTIVER SORTIE VISU DANS ETAT0'
!  WRITE(lunout,*)'entree histclo'
  CALL histclo()

END SUBROUTINE etat0phys_netcdf
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE start_init_orog(lon_in,lat_in,phis,masque)
!
!===============================================================================
! Comment:
!   This routine launch grid_noro, which computes parameters for SSO scheme as
!   described in LOTT & MILLER (1997) and LOTT(1999).
!   In case the file oroparam is present and the key read_orop is activated,
!   grid_noro is bypassed and sub-cell parameters are read from the file.
!===============================================================================
  USE grid_noro_m, ONLY: grid_noro, read_noro
  USE logic_mod,   ONLY: read_orop
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  REAL,    INTENT(IN)    :: lon_in(:), lat_in(:)   ! dim (iml) (jml)
  REAL,    INTENT(INOUT) :: phis(:,:), masque(:,:) ! dim (iml,jml)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname
  INTEGER            :: fid, llm_tmp,ttm_tmp, iml,jml, iml_rel,jml_rel, itau(1)
  INTEGER            :: ierr
  REAL               :: lev(1), date, dt
  REAL, ALLOCATABLE  :: lon_rad(:), lon_ini(:), lon_rel(:,:), relief_hi(:,:)
  REAL, ALLOCATABLE  :: lat_rad(:), lat_ini(:), lat_rel(:,:), tmp_var  (:,:)
  REAL, ALLOCATABLE  :: zmea0(:,:), zstd0(:,:), zsig0(:,:)
  REAL, ALLOCATABLE  :: zgam0(:,:), zthe0(:,:), zpic0(:,:), zval0(:,:)
!-------------------------------------------------------------------------------
  modname="start_init_orog"
  iml=assert_eq(SIZE(lon_in),SIZE(phis,1),SIZE(masque,1),TRIM(modname)//" iml")
  jml=assert_eq(SIZE(lat_in),SIZE(phis,2),SIZE(masque,2),TRIM(modname)//" jml")

!--- HIGH RESOLUTION OROGRAPHY
  CALL flininfo(orofname, iml_rel, jml_rel, llm_tmp, ttm_tmp, fid)

  ALLOCATE(lat_rel(iml_rel,jml_rel),lon_rel(iml_rel,jml_rel))
  CALL flinopen(orofname, .FALSE., iml_rel, jml_rel, llm_tmp, lon_rel, lat_rel,&
                lev, ttm_tmp, itau, date, dt, fid)
  ALLOCATE(relief_hi(iml_rel,jml_rel))
  CALL flinget(fid, orogvar, iml_rel, jml_rel, llm_tmp, ttm_tmp, 1,1, relief_hi)
  CALL flinclo(fid)

!--- IF ANGLES ARE IN DEGREES, THEY ARE CONVERTED INTO RADIANS
  ALLOCATE(lon_ini(iml_rel),lat_ini(jml_rel))
  lon_ini(:)=lon_rel(:,1); IF(MAXVAL(lon_rel)>pi) lon_ini=lon_ini*deg2rad
  lat_ini(:)=lat_rel(1,:); IF(MAXVAL(lat_rel)>pi) lat_ini=lat_ini*deg2rad

!--- FIELDS ARE PROCESSED TO BE ON STANDARD ANGULAR DOMAINS
  ALLOCATE(lon_rad(iml_rel),lat_rad(jml_rel))
  CALL conf_dat2d(orogvar, lon_ini, lat_ini, lon_rad, lat_rad, relief_hi,.FALSE.)
  DEALLOCATE(lon_ini,lat_ini)

!--- COMPUTING THE REQUIRED FIELDS USING ROUTINE grid_noro
  WRITE(lunout,*)
  WRITE(lunout,*)'*** Compute parameters needed for gravity wave drag code ***'

!--- ALLOCATIONS OF SUB-CELL SCALES QUANTITIES
  ALLOCATE(zmea0(iml,jml),zstd0(iml,jml)) !--- Mean orography and std deviation
  ALLOCATE(zsig0(iml,jml),zgam0(iml,jml)) !--- Slope and nisotropy
  zsig0(:,:)=0   !ym uninitialized variable
  zgam0(:,:)=0   !ym uninitialized variable
  ALLOCATE(zthe0(iml,jml))                !--- Highest slope orientation
  zthe0(:,:)=0   !ym uninitialized variable
  ALLOCATE(zpic0(iml,jml),zval0(iml,jml)) !--- Peaks and valley heights

!--- READ SUB-CELL SCALES PARAMETERS FROM A FILE (AT RIGHT RESOLUTION)
  OPEN(UNIT=66,FILE=oroparam,STATUS='OLD',IOSTAT=ierr)
  IF(ierr==0.AND.read_orop) THEN
    CLOSE(UNIT=66)
    CALL read_noro(lon_in,lat_in,oroparam,                                     &
                   phis,zmea0,zstd0,zsig0,zgam0,zthe0,zpic0,zval0,masque)
  ELSE
!--- CALL OROGRAPHY MODULE TO COMPUTE FIELDS
    CALL grid_noro(lon_rad,lat_rad,relief_hi,lon_in,lat_in,                    &
                   phis,zmea0,zstd0,zsig0,zgam0,zthe0,zpic0,zval0,masque)
  END IF
  phis = phis * 9.81
  phis(iml,:) = phis(1,:)
  DEALLOCATE(relief_hi,lon_rad,lat_rad)

!--- PUT QUANTITIES TO PHYSICAL GRID
  CALL gr_dyn_fi(1,iml,jml,klon,zmea0,zmea); DEALLOCATE(zmea0)
  CALL gr_dyn_fi(1,iml,jml,klon,zstd0,zstd); DEALLOCATE(zstd0)
  CALL gr_dyn_fi(1,iml,jml,klon,zsig0,zsig); DEALLOCATE(zsig0)
  CALL gr_dyn_fi(1,iml,jml,klon,zgam0,zgam); DEALLOCATE(zgam0)
  CALL gr_dyn_fi(1,iml,jml,klon,zthe0,zthe); DEALLOCATE(zthe0)
  CALL gr_dyn_fi(1,iml,jml,klon,zpic0,zpic); DEALLOCATE(zpic0)
  CALL gr_dyn_fi(1,iml,jml,klon,zval0,zval); DEALLOCATE(zval0)


END SUBROUTINE start_init_orog
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE start_init_phys(lon_in,lat_in,phis)
!
!===============================================================================
! Purpose:   Compute tsol and qsol, knowing phis.
!===============================================================================
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  REAL,    INTENT(IN) :: lon_in(:),  lat_in(:)       ! dim (iml) (jml2)
  REAL,    INTENT(IN) :: phis(:,:)                   ! dim (iml,jml)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname
  REAL               :: date, dt
  INTEGER            :: iml, jml, jml2, itau(1)
  REAL, ALLOCATABLE  :: lon_rad(:), lon_ini(:), var_ana(:,:)
  REAL, ALLOCATABLE  :: lat_rad(:), lat_ini(:)
  REAL, ALLOCATABLE  :: ts(:,:), qs(:,:)
!-------------------------------------------------------------------------------
  modname="start_init_phys"
  iml=assert_eq(SIZE(lon_in),SIZE(phis,1),TRIM(modname)//" iml")
  jml=SIZE(phis,2); jml2=SIZE(lat_in)

  WRITE(lunout,*)'Opening the surface analysis'
  CALL flininfo(phyfname, iml_phys, jml_phys, llm_phys, ttm_phys, fid_phys)
  WRITE(lunout,*) 'Values read: ',  iml_phys, jml_phys, llm_phys, ttm_phys

  ALLOCATE(lat_phys(iml_phys,jml_phys),lon_phys(iml_phys,jml_phys))
  ALLOCATE(levphys_ini(llm_phys))
  CALL flinopen(phyfname, .FALSE., iml_phys, jml_phys, llm_phys,              &
                lon_phys,lat_phys,levphys_ini,ttm_phys,itau,date,dt,fid_phys)

!--- IF ANGLES ARE IN DEGREES, THEY ARE CONVERTED INTO RADIANS
  ALLOCATE(lon_ini(iml_phys),lat_ini(jml_phys))
  lon_ini(:)=lon_phys(:,1); IF(MAXVAL(lon_phys)>pi) lon_ini=lon_ini*deg2rad
  lat_ini(:)=lat_phys(1,:); IF(MAXVAL(lat_phys)>pi) lat_ini=lat_ini*deg2rad

  ALLOCATE(var_ana(iml_phys,jml_phys),lon_rad(iml_phys),lat_rad(jml_phys))
  CALL get_var_phys(tsrfvar,ts)                   !--- SURFACE TEMPERATURE
  CALL get_var_phys(qsolvar,qs)                   !--- SOIL MOISTURE
  CALL flinclo(fid_phys)
  DEALLOCATE(var_ana,lon_rad,lat_rad,lon_ini,lat_ini)

!--- TSOL AND QSOL ON PHYSICAL GRID
  ALLOCATE(tsol(klon))
  CALL gr_dyn_fi(1,iml,jml,klon,ts,tsol)
  CALL gr_dyn_fi(1,iml,jml,klon,qs,qsol)
  DEALLOCATE(ts,qs)

CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE get_var_phys(title,field)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*),  INTENT(IN)    :: title
  REAL, ALLOCATABLE, INTENT(INOUT) :: field(:,:)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: tllm
!-------------------------------------------------------------------------------
  SELECT CASE(title)
    CASE(psrfvar);         tllm=0
    CASE(tsrfvar,qsolvar); tllm=llm_phys
  END SELECT
  IF(ALLOCATED(field)) RETURN
  ALLOCATE(field(iml,jml)); field(:,:)=0.
  CALL flinget(fid_phys,title,iml_phys,jml_phys,tllm,ttm_phys,1,1,var_ana)
  CALL conf_dat2d(title, lon_ini, lat_ini, lon_rad, lat_rad, var_ana, .TRUE.)
  CALL interp_startvar(title, .TRUE., lon_rad, lat_rad, var_ana,               &
                                      lon_in,  lat_in,  field)

END SUBROUTINE get_var_phys
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE start_init_phys
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE interp_startvar(nam,ibeg,lon,lat,vari,lon2,lat2,varo)
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
  REAL,             INTENT(IN)  :: lon2(:), lat2(:) ! dim (i1) (j2)
  REAL,             INTENT(OUT) :: varo(:,:)        ! dim (i1) (j1)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname
  INTEGER            :: ii, jj, i1, j1, j2
  REAL, ALLOCATABLE  :: vtmp(:,:)
!-------------------------------------------------------------------------------
  modname="interp_startvar"
  ii=assert_eq(SIZE(lon), SIZE(vari,1),TRIM(modname)//" ii")
  jj=assert_eq(SIZE(lat), SIZE(vari,2),TRIM(modname)//" jj")
  i1=assert_eq(SIZE(lon2),SIZE(varo,1),TRIM(modname)//" i1")
  j1=SIZE(varo,2); j2=SIZE(lat2)
  ALLOCATE(vtmp(i1-1,j1))
  IF(ibeg.AND.prt_level>1) THEN
    WRITE(lunout,*)"--------------------------------------------------------"
    WRITE(lunout,*)"$$$ Interpolation barycentrique pour "//TRIM(nam)//" $$$"
    WRITE(lunout,*)"--------------------------------------------------------"
  END IF
  CALL inter_barxy(lon, lat(:jj-1), vari, lon2(:i1-1), lat2, vtmp)
  CALL gr_int_dyn(vtmp, varo, i1-1, j1)

END SUBROUTINE interp_startvar
!
!-------------------------------------------------------------------------------
!
!*******************************************************************************

SUBROUTINE filtreoro(imp1,jmp1,phis,masque,rlatu)

IMPLICIT NONE

  INTEGER imp1,jmp1
  REAL, DIMENSION(imp1,jmp1) :: phis,masque
  REAL, DIMENSION(jmp1) :: rlatu
  REAL, DIMENSION(imp1) :: wwf
  REAL, DIMENSION(imp1,jmp1) :: phiso
  INTEGER :: ifiltre,ifi,ii,i,j
  REAL :: coslat0,ssz

  coslat0=0.5
  phiso=phis
  do j=2,jmp1-1
     print*,'avant if ',cos(rlatu(j)),coslat0
     if (cos(rlatu(j))<coslat0) then
         ! nb de pts affectes par le filtrage de part et d'autre du pt
         ifiltre=(coslat0/cos(rlatu(j))-1.)/2.
         wwf=0.
         do i=1,ifiltre
            wwf(i)=1.
         enddo
         wwf(ifiltre+1)=(coslat0/cos(rlatu(j))-1.)/2.-ifiltre
         do i=1,imp1-1
            if (masque(i,j)>0.9) then
               ssz=phis(i,j)
               do ifi=1,ifiltre+1
                  ii=i+ifi
                  if (ii>imp1-1) ii=ii-imp1+1
                  ssz=ssz+wwf(ifi)*phis(ii,j)
                  ii=i-ifi
                  if (ii<1) ii=ii+imp1-1
                  ssz=ssz+wwf(ifi)*phis(ii,j)
               enddo
               phis(i,j)=ssz*cos(rlatu(j))/coslat0
            endif
         enddo
         print*,'j=',j,coslat0/cos(rlatu(j)), (1.+2.*sum(wwf))*cos(rlatu(j))/coslat0
     endif
  enddo
  call dump2d(imp1,jmp1,phis,'phis ')
  call dump2d(imp1,jmp1,masque,'masque ')
  call dump2d(imp1,jmp1,phis-phiso,'dphis ')

END SUBROUTINE filtreoro


END MODULE etat0phys
