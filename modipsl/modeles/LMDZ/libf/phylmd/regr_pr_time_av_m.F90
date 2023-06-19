! $Id$
MODULE regr_pr_time_av_m

  USE write_field_phy
  USE mod_phys_lmdz_transfert_para, ONLY: bcast
  USE mod_phys_lmdz_para, ONLY: mpi_rank, omp_rank
  USE print_control_mod,  ONLY: prt_level
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Purpose: Regrid pressure with an averaging method. Operations done:
!  * on the root process: read a NetCDF 3D or 4D field at a given day.
!  * pack the fields to the LMDZ "horizontal "physics" grid and scatter
!    to all threads of all processes;
!  * in all the threads of all the processes, regrid the fields in pressure
!    to the LMDZ vertical grid.
!  * the forcing fields are stretched if the following arguments are present:
!     - "lat_in":  input file latitudes.
!     - "Ptrp_ou": target grid (LMDZ) tropopause pressure.
!   so that the tropopause is shifted to the position of the one of LMDZ.
!  Note that the ozone quantity conservation is not ensured.
!-------------------------------------------------------------------------------
! Initial routine: regr_pr_av_m module (L. Guez).
! Present version: David Cugnet ; corresponding additions:
!    - time interpolation
!    - 3D compliant
!    - vertical stretching of the field to allow tropopause matching between
!    input field and current lmdz field.
!-------------------------------------------------------------------------------
! Remarks:
!  * 3D fields are zonal means, with dimensions (latitude, pressure, julian day)
!  * 4D fields have the dimensions:  (longitude, latitude, pressure, julian day)
!  * All the fields are already on the horizontal grid (rlatu) or (rlatu,rlonv),
!    except that the latitudes are in ascending order in the input file.
!  * We assume that all the input fields have the same coordinates.
!  * The target vertical LMDZ grid is the grid of layer centers.
!  * Regridding in pressure can be:
!    - Ploc=='I': pressures provided at Interfaces    => conservative 2nd order
!         OK for ozone coefficients regridding in Cariolle routines.
!    - Ploc=='C': pressures provides at cells Centers => classical linear       
!         OK for ozone vertical regridding, especially when torpopause
!      adjustment is activated, to avoid "strairs shape effect" on profiles.
!  * All the fields are regridded as a single multi-dimensional array, so it
!    saves CPU time to call this procedure once for several NetCDF variables
!    rather than several times, each time for a single NetCDF variable.
!  * The input file pressure at tropopause can be (in decreasing priority):
!    1) read from the file if "tropopause_air_pressure" field is available.
!    2) computed using "tro3" and "tro3_at_tropopause' (if available).
!    3) computed using "tro3" and a fixed threshold otherwise, constant or
!    determined using an empirical three parameters law:
!         o3t(ppbV)=co1+co2*SIN(PI*(month-2)/6)*TANH(lat_deg/co3)
!       => co1 and co2 are in ppbV, and co3 in degrees.
!       => co3 allow a smooth transition between north and south hemispheres.
!  * If available, the field "ps" (input file pressure at surface) is used.
!    If not, the current LMDZ ground pressure is taken instead.
!  * Fields with suffix "m"/"p" are at the closest records earlier/later than
!  the mid-julian day "julien", on the global "dynamics" horizontal grid.
!  * Fields(i,j,k,l) are at longitude-latitude-name "rlonv(i)-rlatu(j)-nam(l)",
!    pressure level/interval (Ploc=="C"/"I") "pcen_in(k)"/"pcen_in(k:k+1)]".
!  * In the 2D file case, the values are the same for all longitudes.
!  * The tropopause correction works like this: the input fields (file) are
!  interpolated on output (LMDZ) pressure field, which is streched using a power
!  law in a limited zone made of 2 layers:
!    1) between lower bound (lower than lowest tropopause) and LMDZ tropopause
!    2) between LMDZ tropopause and upper bound (higher thzn highest tropopause)
!  The stretching function has the following form:
!        Sigma_str = Sigma^(1+alpha*phi(Sigma)), where:
!   * alpha=LOG(SigT_in/SigT_ou)/LOG(SigT_ou)
!     This value shifts the file tropopause to the height of the one of LMDZ.
!   * phi is quasi-linear (sections of 1/log function) in the adjacent intervals:
!       - from 0 to 1 in [Sig_top,SigT_ou]
!       - from 1 to 0 in [SigT_ou,Sig_bot]
!  This quasi-triangular localization function ponderates alpha-law from one near
!  the tropopause to zero each side apart.
!
! * The following fields are on the global "dynamics" grid, as read from files:
  REAL,    SAVE, ALLOCATABLE :: v1 (:,:,:,:)       !--- Current  time ozone fields
! v1: Field read/interpol at time "julien" on the global "dynamics" horiz. grid.
  REAL,    SAVE, ALLOCATABLE :: v1m(:,:,:,:)       !--- Previous time ozone fields
  REAL,    SAVE, ALLOCATABLE :: v1p(:,:,:,:)       !--- Next     time ozone fields
  REAL,    SAVE, ALLOCATABLE :: pgm(:,:), pgp(:,:) !--- Ground     pressure
  REAL,    SAVE, ALLOCATABLE :: ptm(:,:), ptp(:,:) !--- Tropopause pressure
  REAL,    SAVE, ALLOCATABLE :: otm(:,:), otp(:,:) !--- Tropopause o3 mix. ratio
  INTEGER, SAVE :: ntim_in                         !--- Records nb in input file
  INTEGER, SAVE :: itrp0                           !--- idx above chem tropop.
  INTEGER, SAVE :: irec                            !--- Current time index
!      * for daily   input files: current julian day number
!      * for monthly input files: julien is in [time_in(irec),time_in(irec+1)]
  LOGICAL, SAVE :: linterp                         !--- Interpolation in time
  LOGICAL, SAVE :: lPrSfile                        !--- Surface pressure flag
  LOGICAL, SAVE :: lPrTfile                        !--- Tropopause pressure flag
  LOGICAL, SAVE :: lO3Tfile                        !--- Tropopause ozone flag
  LOGICAL, SAVE :: lfirst=.TRUE.                   !--- First call flag
!$OMP THREADPRIVATE(lfirst)
  REAL,    PARAMETER :: pTropUp=9.E+3 !--- Value  <  tropopause pressure (Pa)
  REAL,    PARAMETER :: gamm   =0.4   !--- Max. stretched layer sigma ratio
  REAL,    PARAMETER :: rho    =1.4   !--- Max tropopauses sigma ratio
  REAL,    PARAMETER :: o3t0   =1.E-7 !--- Nominal O3 vmr at tropopause
  LOGICAL, PARAMETER :: lO3Tpara=.FALSE. !--- Parametrized O3 vmr at tropopause
  LOGICAL, PARAMETER :: ldebug=.FALSE.!--- Force writefield_phy multiple outputs
  REAL, PARAMETER :: ChemPTrMin=9.E+3 !--- Thresholds for minimum and maximum
  REAL, PARAMETER :: ChemPTrMax=4.E+4 !    chemical  tropopause pressure (Pa).
  REAL, PARAMETER :: DynPTrMin =8.E+3 !--- Thresholds for minimum and maximum
  REAL, PARAMETER :: DynPTrMax =4.E+4 !    dynamical tropopause pressure (Pa).

CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE regr_pr_time_av(fID, nam, julien, Ploc, Pre_in, Pre_ou, v3, Pgnd_ou,&
                                             time_in, lon_in, lat_in, Ptrp_ou)
!
!-------------------------------------------------------------------------------
  USE dimphy,         ONLY: klon
  USE netcdf95,       ONLY: NF95_INQ_VARID, NF95_INQUIRE_VARIABLE, handle_err, &
                            NF95_INQ_DIMID, NF95_INQUIRE_DIMENSION
  USE netcdf,         ONLY: NF90_INQ_VARID, NF90_GET_VAR, NF90_NOERR
  USE assert_m,       ONLY: assert
  USE assert_eq_m,    ONLY: assert_eq
!!  USE comvert_mod,    ONLY: scaleheight
  USE interpolation,  ONLY: locate
  USE regr_conserv_m, ONLY: regr_conserv
  USE regr_lint_m,    ONLY: regr_lint
  USE slopes_m,       ONLY: slopes
  USE mod_phys_lmdz_para,           ONLY: is_mpi_root,is_master
  USE mod_grid_phy_lmdz,            ONLY: nlon=>nbp_lon, nlat=>nbp_lat, nlev_ou=>nbp_lev, klon_glo, grid_type, unstructured
  USE mod_phys_lmdz_transfert_para, ONLY: scatter2d, scatter, gather
  USE phys_cal_mod,                 ONLY: calend, year_len, days_elapsed, jH_cur
  USE geometry_mod,                 ONLY: ind_cell_glo
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,           INTENT(IN) :: fID        !--- NetCDF file ID
  CHARACTER(LEN=13), INTENT(IN) :: nam(:)     !--- NetCDF variables names
  REAL,              INTENT(IN) :: julien     !--- Days since Jan 1st
  CHARACTER(LEN=1),  INTENT(IN) :: Ploc       !--- Pressures locations
  !--- File/LMDZ (resp. decreasing & increasing order) pressure, Pa
  !    At cells centers or interfaces depending on "Ploc" keyword (C/I)
  REAL,    INTENT(IN)  :: Pre_in(:)           !--- in:  file      (nlev_in[+1])
  REAL,    INTENT(IN)  :: Pre_ou(:,:)         !--- out: LMDZ (klon,nlev_ou[+1])
  REAL,    INTENT(OUT) :: v3(:,:,:)           !--- Regr. fld (klon,nlev_ou,n_var)
  REAL,    INTENT(IN), OPTIONAL :: Pgnd_ou(:) !--- LMDZ ground pressure   (klon)
  REAL,    INTENT(IN), OPTIONAL :: time_in(:) !--- Records times, in days
                                              !    since Jan 1 of current year
  REAL,    INTENT(IN), OPTIONAL :: lon_in(:)  !--- File longitudes vector (klon)
  REAL,    INTENT(IN), OPTIONAL :: lat_in(:)  !--- File latitudes  vector (klon)
  REAL,    INTENT(IN), OPTIONAL :: Ptrp_ou(:) !--- LMDZ tropopause pres   (klon)
!-------------------------------------------------------------------------------
! Local variables:
  include "clesphys.h"
  include "YOMCST.h"
  CHARACTER(LEN=80)  :: sub
  CHARACTER(LEN=320) :: str
  INTEGER :: vID, ncerr, n_var, ibot, iout, nn
  INTEGER :: i, nlev_in, n_dim, itop, itrp, i0
  LOGICAL :: lAdjTro                          !--- Need to adjust tropopause
  REAL    :: y_frac                           !--- Elapsed year fraction
  REAL    :: alpha, beta, al                  !--- For stretching/interpolation
  REAL    :: SigT_in, SigT_ou                 !--- Input and output tropopauses
  REAL    :: Sig_bot, Sig_top                 !--- Bounds of quasi-hat function
  REAL    :: Sig_bo0, Sig_to0                 !--- Lower/upper tropopauses
  REAL    :: Sig_in (SIZE(Pre_in))            !--- Input field sigma levels
  REAL    :: Sig_ou (SIZE(Pre_ou,2))          !--- Output LMDZ sigma levels
  REAL    :: phi    (SIZE(Pre_ou,2))          !--- Stretching exponent anomaly
  REAL    :: Pstr_ou(SIZE(Pre_ou,2))          !--- Stretched pressure levels
  REAL    :: Pres_ou(SIZE(Pre_ou,2))          !--- Pre_ou(i,:), reversed order
  REAL, DIMENSION(nlon, nlat) :: pg1,      pt1,      ot1
  REAL, DIMENSION(klon)       :: Pgnd_in,  Ptrp_in,  Otrp_in
  REAL, DIMENSION(klon)       :: Ptrop_ou, Pgrnd_ou
! * The following fields are scattered to the partial "physics" horizontal grid.
  REAL, POINTER :: v2(:,:,:)                  !--- Current  time ozone fields
!     In the 2D file case, the values are the same for all longitudes.
!     "v2(i, k, l)" is at longitude-latitude "xlon(i)-xlat(i)" and name "nam(l)"
! Both are:          * if Ploc=='I' in pressure interval "press_in_edg(k:k+1)"
!                    * if Ploc=='C' at pressure          "press_in_cen(k)"
  REAL, TARGET :: &
    v2i(klon,SIZE(Pre_in)-1,SIZE(nam)), &     !--- v2 in Ploc=='I' case
    v2c(klon,SIZE(Pre_in)  ,SIZE(nam))        !--- v2 in Ploc=='C' case
  INTEGER,ALLOCATABLE :: ind_cell_glo_glo(:)
  LOGICAL :: ll
!--- For debug
  REAL, DIMENSION(klon)             :: Ptrop_in, Ptrop_ef
  REAL, DIMENSION(klon)             :: dzStrain, dzStrain0
  REAL, DIMENSION(klon,SIZE(Pre_ou,2)) :: Pstrn_ou, phii
!-------------------------------------------------------------------------------
  sub="regr_pr_time_av"
  nlev_in=SIZE(Pre_in); IF(Ploc=='I') nlev_in=nlev_in-1
  IF(Ploc=='I') THEN; v2 => v2i; ELSE; v2 => v2c; END IF
  CALL assert(SIZE(Pre_ou,1)==klon,TRIM(sub)//" Pre_ou klon")
  CALL assert(SIZE(v3,1)==klon,    TRIM(sub)//" v3 klon")
  CALL assert(SIZE(v3,2)==nlev_ou, TRIM(sub)//" v3 nlev_ou")
  IF(Ploc=='I') CALL assert(SIZE(Pre_ou,2)==nlev_ou+1,TRIM(sub)//" Pre_ou nlev_ou+1")
  IF(Ploc=='C') CALL assert(SIZE(Pre_ou,2)==nlev_ou  ,TRIM(sub)//" Pre_ou nlev_ou")
  n_var = assert_eq(SIZE(nam),SIZE(v3,3),TRIM(sub)//" v3 n_var")
  IF(PRESENT(Pgnd_ou)) CALL assert(SIZE(Pgnd_ou)==klon,TRIM(sub)//" Pgnd_ou klon")
  IF(PRESENT(lon_in))  CALL assert(SIZE(lon_in )==klon,TRIM(sub)//" lon_in klon")
  IF(PRESENT(lat_in))  CALL assert(SIZE(lat_in )==klon,TRIM(sub)//" lat_in klon")
  IF(PRESENT(Ptrp_ou)) CALL assert(SIZE(Ptrp_ou)==klon,TRIM(sub)//" Ptrp_ou klon")
  lAdjTro=PRESENT(Ptrp_ou)
  IF(lAdjTro) THEN
    IF(.NOT.PRESENT(lat_in)) &
      CALL abort_physic(sub, 'Missing lat_in (required if adjust_tropopause=T)', 1)
    IF(.NOT.PRESENT(Pgnd_ou).AND.Ploc=='C') &
      CALL abort_physic(sub, 'Missing ground Pr(required if adjust_tropopause=T)', 1)
    IF(PRESENT(Pgnd_ou)) THEN; Pgrnd_ou=Pgnd_ou; ELSE; Pgrnd_ou=Pre_ou(:,1); END IF
  END IF

  !$OMP MASTER
  IF(is_mpi_root) THEN

    !=== CHECK WHICH FIELDS ARE AVAILABLE IN THE INPUT FILE
    IF(lfirst) THEN
      lPrSfile=lAdjTro.AND.NF90_INQ_VARID(fID,"ps"                     ,vID)==NF90_NOERR
      lPrTfile=lAdjTro.AND.NF90_INQ_VARID(fID,"tropopause_air_pressure",vID)==NF90_NOERR
      lO3Tfile=lAdjTro.AND.NF90_INQ_VARID(fID,"tro3_at_tropopause"     ,vID)==NF90_NOERR
      CALL NF95_INQ_DIMID(fID,"time",vID)
      CALL NF95_INQUIRE_DIMENSION(fID,vID,nclen=ntim_in)
      linterp=PRESENT(time_in).AND.ntim_in==14
      ALLOCATE(v1(nlon,nlat,nlev_in,n_var))
      IF(linterp) THEN
        ALLOCATE(v1m(nlon,nlat,nlev_in,n_var),v1p(nlon,nlat,nlev_in,n_var))
        IF(lPrSfile) ALLOCATE(pgm(nlon,nlat),pgp(nlon,nlat))
        IF(lPrTfile) ALLOCATE(ptm(nlon,nlat),ptp(nlon,nlat))
        IF(lO3Tfile) ALLOCATE(otm(nlon,nlat),otp(nlon,nlat))
      END IF
      !--- INITIAL INDEX: LOCATE A LAYER WELL ABOVE TROPOPAUSE (50hPa)
      IF(lAdjTro) itrp0=locate(Pre_in,pTropUp)
      CALL msg(linterp,'Monthly O3 files => ONLINE TIME INTERPOLATION.'    ,sub)
      CALL msg(lPrSfile,'Using GROUND PRESSURE from input O3 forcing file.',sub)
      CALL msg(lAdjTro ,'o3 forcing file tropopause location uses:'        ,sub)
      IF(lPrTfile)      THEN; str='    INPUT FILE PRESSURE'
      ELSE IF(lO3Tfile) THEN; str='    INPUT FILE O3 CONCENTRATION'
      ELSE IF(lO3Tpara) THEN; str='    PARAMETRIZED O3 concentration'
      ELSE;                   str='    CONSTANT O3 concentration'; END IF
      CALL msg(lAdjTro,TRIM(str)//' at tropopause')
    END IF

    !=== UPDATE (ALWAYS FOR DAILY FILES, EACH MONTH FOR MONTHLY FILES)
    CALL update_fields()

    !=== TIME INTERPOLATION FOR MONTHLY INPUT FILES
    IF(linterp) THEN
      WRITE(str,'(a,f12.8,2(a,f5.1))')'Interpolating O3 at julian day ',julien,&
        ' from fields at times ',time_in(irec),' and ', time_in(irec+1)
      CALL msg(.TRUE.,str,sub)
      al=(time_in(irec+1)-julien)/(time_in(irec+1)-time_in(irec))
      v1=al*v1m+(1.-al)*v1p
      IF(lPrSfile) pg1=al*pgm+(1.-al)*pgp
      IF(lPrTfile) pt1=al*ptm+(1.-al)*ptp
      IF(lO3Tfile) ot1=al*otm+(1.-al)*otp
    END IF
  ELSE
    IF (lfirst) ALLOCATE(v1(0,0,0,0))
  END IF
  !$OMP END MASTER
  !$OMP BARRIER
  IF(lfirst) THEN
    lfirst=.FALSE.;       CALL bcast(lfirst)
    IF(lAdjTro)           CALL bcast(itrp0)
    CALL bcast(lPrSfile); CALL bcast(lPrTfile)
    CALL bcast(lO3Tfile); CALL bcast(linterp)
  END IF
  
  IF (is_master) THEN 
    ALLOCATE(ind_cell_glo_glo(klon_glo))
  ELSE
    ALLOCATE(ind_cell_glo_glo(0))
  ENDIF
  CALL gather(ind_cell_glo,ind_cell_glo_glo)
  IF (is_master .AND. grid_type==unstructured) v1(:,:,:,:)=v1(:,ind_cell_glo_glo(:),:,:)
  
  CALL scatter2d(v1,v2)

  !--- No "ps" in input file => assumed to be equal to current LMDZ ground press
  IF(lPrSfile) THEN
    IF (is_master .AND. grid_type==unstructured) pg1(:,:)=pg1(:,ind_cell_glo_glo(:))
    CALL scatter2d(pg1,Pgnd_in)
  ELSE
    Pgnd_in=Pre_ou(:,1)
  END IF

  IF(lPrTfile) THEN
    IF (is_master .AND. grid_type==unstructured) pt1(:,:)=pt1(:,ind_cell_glo_glo(:))
    CALL scatter2d(pt1,Ptrp_in)
  ENDIF

  IF(lO3Tfile) THEN
    IF (is_master .AND. grid_type==unstructured) ot1(:,:)=ot1(:,ind_cell_glo_glo(:))
    CALL scatter2d(ot1,Otrp_in)
  ENDIF
  !--- No ground pressure in input file => choose it to be the one of LMDZ
  IF(lAdjTro.AND..NOT.lPrSfile) Pgnd_in(:)=Pgrnd_ou(:)

  !--- REGRID IN PRESSURE ; 3rd index inverted because "paprs" is decreasing
  IF(.NOT.lAdjTro) THEN
    DO i=1,klon
      Pres_ou=Pre_ou(i,SIZE(Pre_ou,2):1:-1)   !--- pplay & paprs are decreasing
      IF(Ploc=='C') CALL regr_lint   (1,v2(i,:,:), LOG(Pre_in(:)),             &
        LOG(Pres_ou(:)), v3(i,nlev_ou:1:-1,:))
      IF(Ploc=='I') CALL regr_conserv(1,v2(i,:,:),     Pre_in(:) ,             &
            Pres_ou(:) , v3(i,nlev_ou:1:-1,:), slopes(1,v2(i,:,:), Pre_in(:)))
    END DO
!-------------------------------------------------------------------------------
  ELSE                        !--- REGRID IN PRESSURE ; TROPOPAUSE ADJUSTMENT
!-------------------------------------------------------------------------------
    y_frac=(REAL(days_elapsed)+jH_cur)/year_len

    !--- OUTPUT SIGMA LEVELS
    DO i=1,klon

      !--- INPUT/OUTPUT (FILE/LMDZ) SIGMA LEVELS IN CURRENT COLUMN
      Pres_ou   = Pre_ou(i,SIZE(Pre_ou,2):1:-1)!--- pplay & paprs are decreasing
      Sig_in(:) = Pre_in (:)/Pgnd_in(i)            !--- increasing values
      Sig_ou(:) = Pres_ou(:)/Pgnd_ou(i)            !--- increasing values

      !--- INPUT (FILE) SIGMA LEVEL AT TROPOPAUSE ; extreme values are filtered
      ! to keep tropopause pressure realistic ; high values are usually due to
      ! ozone hole fooling the crude chemical tropopause detection algorithm.
      SigT_in = get_SigTrop(i,itrp)
      SigT_in=MIN(SigT_in,ChemPTrMax/Pgnd_in(i))   !--- too low  value filtered
      SigT_in=MAX(SigT_in,ChemPTrMin/Pgnd_ou(i))   !--- too high value filtered

      !--- OUTPUT (LMDZ) SIGMA LEVEL AT TROPOPAUSE ; too high variations of the
      ! dynamical tropopause (especially in filaments) are conterbalanced with
      ! a filter ensuring it stays within a certain distance around input (file)
      ! tropopause, hence avoiding avoid a too thick stretched region ; a final
      ! extra-safety filter keeps the tropopause pressure value realistic.
      SigT_ou = Ptrp_ou(i)/Pgnd_ou(i)
      IF(SigT_ou<SigT_in/rho) SigT_ou=SigT_in/rho  !--- too low  value w/r input
      IF(SigT_ou>SigT_in*rho) SigT_ou=SigT_in*rho  !--- too high value w/r input
      SigT_ou=MIN(SigT_ou,DynPTrMax/Pgnd_ou(i))    !--- too low  value filtered
      SigT_ou=MAX(SigT_ou,DynPTrMin/Pgnd_ou(i))    !--- too high value filtered
      Ptrop_ou(i)=SigT_ou*Pgnd_ou(i)
      iout = locate(Sig_ou(:),SigT_ou)

      !--- POWER LAW COEFFICIENT FOR TROPOPAUSES MATCHING
      alpha = LOG(SigT_in/SigT_ou)/LOG(SigT_ou)

      !--- DETERMINE STRETCHING DOMAIN UPPER AND LOWER BOUNDS
      Sig_bo0 = MAX(SigT_in,SigT_ou)               !--- lowest  tropopause
      Sig_to0 = MIN(SigT_in,SigT_ou)               !--- highest tropopause
      beta    = (Sig_bo0/Sig_to0)**gamm            !--- stretching exponent
      Sig_bot = MIN(Sig_bo0*beta,0.1*(9.+Sig_bo0)) !--- must be <1
      ibot = locate(Sig_ou(:),Sig_bot)             !--- layer index
      IF(ibot-iout<2) THEN                         !--- at least one layer thick
        ibot=MIN(iout+2,nlev_ou); Sig_bot=Sig_ou(ibot)
      END IF
      Sig_top = Sig_to0/beta                       !--- upper bound
      itop = locate(Sig_ou(:),Sig_top)             !--- layer index
      IF(iout-itop<2) THEN                         !--- at least one layer thick
        itop=MAX(iout-2,1); Sig_top=Sig_ou(itop)
      END IF

      !--- STRETCHING POWER LAW LOCALIZATION FUNCTION:
      !    0 in [0,Sig_top]    0->1 in [Sig_top,SigT_ou]
      !    0 in [Sig_bot,1]    1->0 in [SigT_ou, Sig_bot]
      phi(:)=0.
      phi(itop+1:iout) = (1.-LOG(Sig_top)/LOG(Sig_ou(itop+1:iout)))&
                            *LOG(SigT_ou)/LOG(SigT_ou/Sig_top)
      phi(iout+1:ibot) = (1.-LOG(Sig_bot)/LOG(Sig_ou(iout+1:ibot)))&
                            *LOG(SigT_ou)/LOG(SigT_ou/Sig_bot)

      !--- LOCALY STRECHED OUTPUT (LMDZ) PRESSURE PROFILES (INCREASING ORDER)
      Pstr_ou(:) = Pres_ou(:) * Sig_ou(:)**(alpha*phi(:))

      !--- REGRID INPUT PROFILE ON STRAINED VERTICAL OUTPUT LEVELS
      IF(Ploc=='C') CALL regr_lint   (1, v2(i,:,:), LOG(Pre_in(:)),            &
        LOG(Pstr_ou(:)), v3(i,nlev_ou:1:-1,:))
      IF(Ploc=='I') CALL regr_conserv(1, v2(i,:,:),     Pre_in(:) ,            &
            Pstr_ou(:) , v3(i,nlev_ou:1:-1,:), slopes(1,v2(i,:,:), Pre_in(:)))

      !--- CHECK CONCENTRATIONS. strato: 50ppbV-15ppmV ; tropo: 5ppbV-300ppbV.
      i0=nlev_ou-locate(Pres_ou(:),Ptrop_ou(i))+1
      ll=check_ozone(v3(i, 1:i0-1   ,1),lon_in(i),lat_in(i),1 ,'troposphere',  &
                     5.E-9,3.0E-7)
!     IF(ll) CALL abort_physic(sub, 'Inconsistent O3 values in troposphere', 1)
      ll=check_ozone(v3(i,i0:nlev_ou,1),lon_in(i),lat_in(i),i0,'stratosphere', &
                     5.E-8,1.5E-5)
!     IF(ll) CALL abort_physic(sub, 'Inconsistent O3 values in stratosphere', 1)

      IF(ldebug) THEN
        dzStrain0(i) = SIGN(7.*LOG(Sig_bo0/Sig_to0),SigT_in-SigT_ou)
        dzStrain (i) = SIGN(7.*LOG(Sig_bot/Sig_top),SigT_in-SigT_ou)
        Ptrop_in (i) = SigT_in*Pgnd_in(i)
        Pstrn_ou(i,:)= Pstr_ou
        phii(i,:)    = phi(:)
        Ptrop_ef(i)  = PTrop_chem(i, itrp, locate(Pres_ou(:),PTropUp),    &
                             Pres_ou(:), v3(:,nlev_ou:1:-1,1),o3trop=o3t0)
      END IF
    END DO
  END IF
  IF(ldebug.AND.lAdjTro) THEN
    CALL writefield_phy('PreSt_ou' ,Pstrn_ou,SIZE(Pre_ou,2)) !--- Strained Pres
    CALL writefield_phy('dzStrain' ,dzStrain ,1)     !--- Strained thickness
    CALL writefield_phy('dzStrain0',dzStrain0,1)     !--- Tropopauses distance
    CALL writefield_phy('phi',phii,nlev_ou)          !--- Localization function
    !--- Tropopauses pressures:
    CALL writefield_phy('PreTr_in',Ptrop_in,1)       !--- Input and effective
    CALL writefield_phy('PreTr_ou',Ptrop_ou,1)       !--- LMDz dyn tropopause
    CALL writefield_phy('PreTr_ef',Ptrop_ef,1)       !--- Effective chem tropop
  END IF
  IF(ldebug) THEN
    CALL writefield_phy('Ozone_in',v2(:,:,1),nlev_in)!--- Raw input O3 field
    CALL writefield_phy('Ozone_ou',v3(:,:,1),nlev_ou)!--- Output ozone field
    CALL writefield_phy('Pres_ou' ,Pre_ou,SIZE(Pre_ou,2))!--- LMDZ Pressure
  END IF

CONTAINS


!-------------------------------------------------------------------------------
!
SUBROUTINE update_fields()
!
!-------------------------------------------------------------------------------
  IF(.NOT.linterp) THEN                 !=== DAILY FILES: NO TIME INTERPOLATION
    CALL msg(.TRUE.,sub,'Updating Ozone forcing field: read from file.')
    irec=MIN(INT(julien)+1,ntim_in)     !--- irec is just the julian day number
    !--- MIN -> Security in the unlikely case of roundup errors.
    CALL get_3Dfields(v1)               !--- Read ozone field(s)
    IF(lAdjTro) THEN                    !--- Additional files for fields strain
      IF(lPrSfile) CALL get_2Dfield(pg1,"ps")
      IF(lPrTfile) CALL get_2Dfield(pt1,"tropopause_air_pressure")
      IF(lO3Tfile) CALL get_2Dfield(ot1,"tro3_at_tropopause")
    END IF
  ELSE                                  !=== MONTHLY FILES: GET 2 NEAREST RECS
    IF(lfirst) irec=locate(time_in,julien) !--- Need to locate surrounding times
    IF(.NOT.lfirst.AND.julien<time_in(irec+1)) RETURN
    CALL msg(.TRUE.,'Refreshing adjacent Ozone forcing fields.',sub)
    IF(lfirst) THEN                     !=== READ EARLIEST TIME FIELDS
      WRITE(str,'(a,i3,a,f12.8,a)')'Previous available field update (step 1): '&
      //'reading record ',irec,' (time ',time_in(irec),')'
      CALL msg(.TRUE.,str,sub)
      CALL get_3Dfields(v1m)            !--- Read ozone field(s)
      IF(lAdjTro) THEN                  !--- Additional files for fields strain
        IF(lPrSfile) CALL get_2Dfield(pgm,"ps")
        IF(lPrTfile) CALL get_2Dfield(ptm,"tropopause_air_pressure")
        IF(lO3Tfile) CALL get_2Dfield(otm,"tro3_at_tropopause")
      END IF
    ELSE                                !=== SHIFT FIELDS
      irec=irec+1
      WRITE(str,'(a,i3,a,f12.8,a)')'Previous available field update: shifting'&
      //' current next one (',irec,', time ',time_in(irec),')'
      CALL msg(.TRUE.,str,sub)
      v1m=v1p                           !--- Ozone fields
      IF(lAdjTro) THEN                  !--- Additional files for fields strain
        IF(lPrSfile) pgm=pgp             !--- Surface pressure
        IF(lPrTfile) ptm=ptp             !--- Tropopause pressure
        IF(lO3Tfile) otm=otp             !--- Tropopause ozone
      END IF
    END IF
    irec=irec+1
    WRITE(str,'(a,i3,a,f12.8,a)')'Next available field update: reading record'&
    ,irec,' (time ',time_in(irec),')'
    CALL msg(.TRUE.,str,sub)
    CALL get_3Dfields(v1p)              !--- Read ozone field(s)
    IF(lAdjTro) THEN                    !--- Additional files for fields strain
      IF(lPrSfile) CALL get_2Dfield(pgp,"ps")
      IF(lPrTfile) CALL get_2Dfield(ptp,"tropopause_air_pressure")
      IF(lO3Tfile) CALL get_2Dfield(otp,"tro3_at_tropopause")
    END IF
    irec=irec-1
  END IF

END SUBROUTINE update_fields
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE get_2Dfield(v,var)
!
!-------------------------------------------------------------------------------
! Purpose: Shortcut to get the "irec"th record of the surface field named "var"
!          from the input file.
! Remark: In case the field is zonal, it is duplicated along first dimension.
!-------------------------------------------------------------------------------
! Arguments:
  REAL,             INTENT(INOUT) :: v(:,:)
  CHARACTER(LEN=*), INTENT(IN)    :: var
!-------------------------------------------------------------------------------
  CALL NF95_INQ_VARID(fID, TRIM(var), vID)
  CALL NF95_INQUIRE_VARIABLE(fID, vID, ndims=n_dim)
  IF(n_dim==2) ncerr=NF90_GET_VAR(fID,vID,v(1,:), start=[  1,irec])
  IF(n_dim==3) ncerr=NF90_GET_VAR(fID,vID,v(:,:), start=[1,1,irec])
  CALL handle_err(TRIM(sub)//" NF90_GET_VAR "//TRIM(var),ncerr,fID)

  !--- Flip latitudes: ascending in input file, descending in "rlatu".
  IF(n_dim==3) THEN
    v(1,:) = v(1,nlat:1:-1)
    v(2:,:)= SPREAD(v(1,:),DIM=1,ncopies=nlon-1)  !--- Duplication
  ELSE
    v(:,:) = v(:,nlat:1:-1)
  END IF

END SUBROUTINE get_2Dfield
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE get_3Dfields(v)
!
!-------------------------------------------------------------------------------
! Purpose: Shortcut to get the "irec"th record of the 3D fields named "nam"
!          from the input file. Fields are stacked on fourth dimension.
! Remark: In case the field is zonal, it is duplicated along first dimension.
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(INOUT) :: v(:,:,:,:)
!-------------------------------------------------------------------------------
  DO i=1,SIZE(nam)
    CALL NF95_INQ_VARID(fID, TRIM(nam(i)), vID)
    CALL NF95_INQUIRE_VARIABLE(fID, vID, ndims=n_dim)
    IF(n_dim==3) ncerr=NF90_GET_VAR(fID,vID,v(1,:,:,i), start=[  1,1,irec])
    IF(n_dim==4) ncerr=NF90_GET_VAR(fID,vID,v(:,:,:,i), start=[1,1,1,irec])
    CALL handle_err(TRIM(sub)//" NF90_GET_VAR "//TRIM(nam(i)),ncerr,fID)
  END DO

  !--- Flip latitudes: ascending in input file, descending in "rlatu".
  IF(n_dim==3) THEN
    v(1,:,:,:) = v(1,nlat:1:-1,:,:)
    v(2:,:,:,:)= SPREAD(v(1,:,:,:),DIM=1,ncopies=nlon-1)  !--- Duplication
  ELSE
    v(:,:,:,:) = v(:,nlat:1:-1,:,:)
  END IF

END SUBROUTINE get_3Dfields
!
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
!
FUNCTION get_SigTrop(ih,it) RESULT(out)
!
!-------------------------------------------------------------------------------
! Arguments:
  REAL                 :: out
  INTEGER, INTENT(IN)  :: ih
  INTEGER, INTENT(OUT) :: it
!-------------------------------------------------------------------------------
  !--- Pressure at tropopause read from the forcing file
       IF(lPrTfile) THEN; out=Ptrp_in(ih)/Pgnd_in(ih); RETURN; END IF

  !--- Chemical tropopause definition based on a particular threshold
       IF(lO3Tfile) THEN; out=PTrop_chem(ih,it,itrp0,Pre_in,v2(:,:,1),Otrp_in(ih))
  ELSE IF(lO3Tpara) THEN; out=PTrop_chem(ih,it,itrp0,Pre_in,v2(:,:,1))
  ELSE                  ; out=PTrop_chem(ih,it,itrp0,Pre_in,v2(:,:,1),o3t0); END IF
  out=out/Pgnd_in(ih)

END FUNCTION get_SigTrop
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
FUNCTION PTrop_chem(ih,it,it0,pres,o3,o3trop) RESULT(out)
!
!-------------------------------------------------------------------------------
! Purpose: Determine the tropopause using chemical criterium, ie the pressure
!          at which the ozone concentration reaches a certain level.
!-------------------------------------------------------------------------------
! Remarks:
! * Input field is upside down (increasing pressure // increasing vertical idx)
!   The sweep is done from top to ground, starting at the 50hPa layer (idx it0),
!   where O3 is about 1,5 ppmV, until the first layer where o3<o3t is reached:
!   the O3 profile in this zone is decreasing with pressure.
! * Threshold o3t can be given as an input argument ("o3trop", in ppmV) or
!   determined using an empirical law roughly derived from ... & al.
!-------------------------------------------------------------------------------
! Arguments:
  REAL                        :: out           !--- Pressure at tropopause
  INTEGER,        INTENT(IN)  :: ih            !--- Horizontal index
  INTEGER,        INTENT(OUT) :: it            !--- Index of tropopause layer
  INTEGER,        INTENT(IN)  :: it0           !--- Idx: higher than tropopause
  REAL,           INTENT(IN)  :: pres(:)       !--- Pressure profile, increasing
  REAL,           INTENT(IN)  :: o3(:,:)       !--- Ozone field (pptV)
  REAL, OPTIONAL, INTENT(IN)  :: o3trop        !--- Ozone at tropopause
!-------------------------------------------------------------------------------
! Local variables:
  REAL :: o3t                                  !--- Ozone concent. at tropopause
  REAL :: al                                   !--- Interpolation coefficient
  REAL :: coef                                 !--- Coeff of latitude modulation
  REAL, PARAMETER :: co3(3)=[91.,28.,20.]      !--- Coeff for o3 at tropopause
!-------------------------------------------------------------------------------
  !--- DETERMINE THE OZONE CONCENTRATION AT TROPOPAUSE IN THE CURRENT COLUMN
  IF(PRESENT(o3trop)) THEN                     !=== THRESHOLD FROM ARGUMENTS
    o3t=o3trop
  ELSE                                         !=== EMPIRICAL LAW
    coef = TANH(lat_in(ih)/co3(3))             !--- Latitude  modulation
    coef = SIN (2.*RPI*(y_frac-1./6.)) * coef  !--- Seasonal  modulation
    o3t  = 1.E-9 * (co3(1)+co3(2)*coef)        !--- Value from parametrization
  END IF

  !--- FROM 50hPa, GO DOWN UNTIL "it" IS SUCH THAT o3(ih,it)>o3t>o3(ih,it+1)
  it=it0; DO WHILE(o3(ih,it+1)>=o3t); it=it+1; END DO
  al=(o3(ih,it)-o3t)/(o3(ih,it)-o3(ih,it+1))
  IF(Ploc=='C') out =      pres(it)**(1.-al) * pres(it+1)**al
  IF(Ploc=='I') out = SQRT(pres(it)**(1.-al) * pres(it+2)**al *pres(it+1))
  it = locate(pres(:), out)                    !--- pres(it)<Ptrop<pres(it+1)

END FUNCTION PTrop_chem
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
FUNCTION check_ozone(o3col, lon, lat, ilev0, layer, vmin, vmax) RESULT(out)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  LOGICAL                      :: out          !--- .T. => some wrong values
  REAL,             INTENT(IN) :: o3col(:), lon, lat
  INTEGER,          INTENT(IN) :: ilev0
  CHARACTER(LEN=*), INTENT(IN) :: layer
  REAL, OPTIONAL,   INTENT(IN) :: vmin
  REAL, OPTIONAL,   INTENT(IN) :: vmax
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: k
  LOGICAL :: lmin, lmax
  REAL    :: fac
  CHARACTER(LEN=6) :: unt
!-------------------------------------------------------------------------------
  !--- NOTHING TO DO
  lmin=.FALSE.; IF(PRESENT(vmin)) lmin=COUNT(o3col<vmin)/=0
  lmax=.FALSE.; IF(PRESENT(vmax)) lmax=COUNT(o3col>vmax)/=0
  out=lmin.OR.lmax; IF(.NOT.out.OR.prt_level>100) RETURN

  !--- SOME TOO LOW VALUES FOUND
  IF(lmin) THEN
    CALL unitc(vmin,unt,fac)
    DO k=1,SIZE(o3col); IF(o3col(k)>vmin) CYCLE
      WRITE(*,'(a,2f7.1,i3,a,2(f8.3,a))')'WARNING: inconsistent value in '     &
        //TRIM(layer)//': O3(',lon,lat,k+ilev0-1,')=',fac*o3col(k),unt//' < ', &
        fac*vmin,unt//' in '//TRIM(layer)
    END DO
  END IF

  !--- SOME TOO HIGH VALUES FOUND
  IF(lmax) THEN
    CALL unitc(vmax,unt,fac)
    DO k=1,SIZE(o3col); IF(o3col(k)<vmax) CYCLE
      WRITE(*,'(a,2f7.1,i3,a,2(f8.3,a))')'WARNING: inconsistent value in '     &
        //TRIM(layer)//': O3(',lon,lat,k+ilev0-1,')=',fac*o3col(k),unt//' > ', &
        fac*vmax,unt//' in '//TRIM(layer)
    END DO
  END IF

END FUNCTION check_ozone
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE unitc(val,unt,fac)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  REAL,             INTENT(IN)  :: val
  CHARACTER(LEN=6), INTENT(OUT) :: unt
  REAL,             INTENT(OUT) :: fac
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: ndgt
!-------------------------------------------------------------------------------
  ndgt=3*FLOOR(LOG10(val)/3.)
  SELECT CASE(ndgt)
    CASE(-6);     unt=' ppmV '; fac=1.E6
    CASE(-9);     unt=' ppbV '; fac=1.E9
    CASE DEFAULT; unt=' vmr  '; fac=1.0
  END SELECT

END SUBROUTINE unitc
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE msg(ll,str,sub)
!
!-------------------------------------------------------------------------------
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  LOGICAL,                    INTENT(IN) :: ll
  CHARACTER(LEN=*),           INTENT(IN) :: str
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: sub
!-------------------------------------------------------------------------------
  IF(.NOT.ll) RETURN
  IF(PRESENT(sub)) THEN
    WRITE(lunout,*)TRIM(sub)//': '//TRIM(str)
  ELSE
    WRITE(lunout,*)TRIM(str)
  END IF

END SUBROUTINE msg
!
!-------------------------------------------------------------------------------

END SUBROUTINE regr_pr_time_av
!
!-------------------------------------------------------------------------------

END MODULE regr_pr_time_av_m
