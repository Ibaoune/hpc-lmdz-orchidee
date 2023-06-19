!
! phys_local_var_mod.F90 1327 2010-03-17 15:33:56Z idelkadi $

MODULE phys_output_var_mod

  USE dimphy
  ! Variables outputs pour les ecritures des sorties
  !======================================================================
  !
  !
  !======================================================================
  ! Declaration des variables

  REAL, SAVE, ALLOCATABLE :: snow_o(:), zfra_o(:)
  !$OMP THREADPRIVATE(snow_o, zfra_o)
  REAL, SAVE, ALLOCATABLE :: sza_o(:) ! solar zenithal angle
  !$OMP THREADPRIVATE(sza_o)
  INTEGER, SAVE, ALLOCATABLE ::  itau_con(:)       ! Nombre de pas ou rflag <= 1
  !$OMP THREADPRIVATE(itau_con)
  REAL, SAVE, ALLOCATABLE :: bils_ec(:) ! Contribution of energy conservation
  REAL, SAVE, ALLOCATABLE :: bils_ech(:) ! Contribution of energy conservation
  REAL, SAVE, ALLOCATABLE :: bils_tke(:) ! Contribution of energy conservation
  REAL, SAVE, ALLOCATABLE :: bils_diss(:) ! Contribution of energy conservation
  REAL, SAVE, ALLOCATABLE :: bils_kinetic(:) ! bilan de chaleur au sol, kinetic
  REAL, SAVE, ALLOCATABLE :: bils_enthalp(:) ! bilan de chaleur au sol
  REAL, SAVE, ALLOCATABLE :: bils_latent(:) ! bilan de chaleur au sol
  !$OMP THREADPRIVATE(bils_ec,bils_ech,bils_tke,bils_diss,bils_kinetic,bils_enthalp,bils_latent)
  ! output variables for energy conservation tests, computed in add_phys_tend
  REAL, SAVE, ALLOCATABLE :: d_qw_col(:)      ! watter vapour mass budget for each column (kg/m2/s)
  REAL, SAVE, ALLOCATABLE :: d_ql_col(:)      ! liquid watter mass budget for each column (kg/m2/s)
  REAL, SAVE, ALLOCATABLE :: d_qs_col(:)      ! solid watter mass budget for each column (kg/m2/s)
  REAL, SAVE, ALLOCATABLE :: d_qt_col(:)      ! total watter mass budget for each column (kg/m2/s)
  REAL, SAVE, ALLOCATABLE :: d_ek_col(:)      ! kinetic energy budget for each column (W/m2)
  REAL, SAVE, ALLOCATABLE :: d_h_dair_col(:)  ! enthalpy budget of dry air for each column (W/m2)
  REAL, SAVE, ALLOCATABLE :: d_h_qw_col(:)    ! enthalpy budget of watter vapour for each column (W/m2)
  REAL, SAVE, ALLOCATABLE :: d_h_ql_col(:)    ! enthalpy budget of liquid watter for each column (W/m2)
  REAL, SAVE, ALLOCATABLE :: d_h_qs_col(:)    ! enthalpy budget of solid watter  for each column (W/m2)
  REAL, SAVE, ALLOCATABLE :: d_h_col(:)       ! total enthalpy budget for each column (W/m2)
  !$OMP THREADPRIVATE(d_qw_col, d_ql_col, d_qs_col, d_qt_col, d_ek_col, d_h_dair_col)
  !$OMP THREADPRIVATE(d_h_qw_col, d_h_ql_col, d_h_qs_col, d_h_col)

  ! Outputs used in cloudth_vert to extract the moments of the horizontal and 
  ! vertical PDFs
  REAL, SAVE, ALLOCATABLE :: cloudth_sth(:,:),cloudth_senv(:,:)
  !$OMP THREADPRIVATE(cloudth_sth,cloudth_senv)
  REAL, SAVE, ALLOCATABLE :: cloudth_sigmath(:,:),cloudth_sigmaenv(:,:)
  !$OMP THREADPRIVATE(cloudth_sigmath,cloudth_sigmaenv)

! Marine
! Variables de sortie du simulateur AIRS

  REAL, SAVE, ALLOCATABLE :: map_prop_hc(:),map_prop_hist(:),alt_tropo(:)
  !$OMP THREADPRIVATE(map_prop_hc,map_prop_hist,alt_tropo)
  REAL, SAVE, ALLOCATABLE :: map_emis_hc(:),map_iwp_hc(:),map_deltaz_hc(:), &
                       map_pcld_hc(:),map_tcld_hc(:)
  !$OMP THREADPRIVATE(map_emis_hc,map_iwp_hc,map_deltaz_hc,map_pcld_hc,map_tcld_hc)
  REAL, SAVE, ALLOCATABLE :: map_emis_hist(:),map_iwp_hist(:),map_deltaz_hist(:),map_rad_hist(:)         
  !$OMP THREADPRIVATE(map_emis_hist,map_iwp_hist,map_deltaz_hist,map_rad_hist)
  REAL, SAVE, ALLOCATABLE :: map_ntot(:),map_hc(:),map_hist(:)
  REAL, SAVE, ALLOCATABLE :: map_Cb(:),map_ThCi(:),map_Anv(:)
  !$OMP THREADPRIVATE(map_ntot,map_hc,map_hist,map_Cb,map_ThCi,map_Anv)
  REAL, SAVE, ALLOCATABLE :: map_emis_Cb(:),map_pcld_Cb(:),map_tcld_Cb(:)
  REAL, SAVE, ALLOCATABLE :: map_emis_ThCi(:),map_pcld_ThCi(:),map_tcld_ThCi(:)
  !$OMP THREADPRIVATE(map_emis_Cb,map_pcld_Cb,map_tcld_Cb,map_emis_ThCi)
  REAL, SAVE, ALLOCATABLE :: map_emis_Anv(:),map_pcld_Anv(:),map_tcld_Anv(:)
  !$OMP THREADPRIVATE(map_pcld_ThCi,map_tcld_ThCi,map_emis_Anv,map_pcld_Anv,map_tcld_Anv)              
   

  ! ug Plein de variables venues de phys_output_mod
  INTEGER, PARAMETER                           :: nfiles = 10
  LOGICAL, DIMENSION(nfiles), SAVE             :: clef_files
  LOGICAL, DIMENSION(nfiles), SAVE             :: clef_stations
  INTEGER, DIMENSION(nfiles), SAVE             :: lev_files
  INTEGER, DIMENSION(nfiles), SAVE             :: nid_files
  INTEGER, DIMENSION(nfiles), SAVE  :: nnid_files
  !$OMP THREADPRIVATE(clef_files, clef_stations, lev_files,nid_files,nnid_files)
  INTEGER, DIMENSION(nfiles), SAVE :: nnhorim

  INTEGER, DIMENSION(nfiles), SAVE :: nhorim, nvertm
  INTEGER, DIMENSION(nfiles), SAVE :: nvertap, nvertbp, nvertAlt
  REAL, DIMENSION(nfiles), SAVE                :: zoutm
  CHARACTER(LEN=20), DIMENSION(nfiles), SAVE   :: type_ecri
  !$OMP THREADPRIVATE(nnhorim, nhorim, nvertm, zoutm,type_ecri)
  CHARACTER(LEN=20), DIMENSION(nfiles), SAVE  :: type_ecri_files, phys_out_filetypes
  !$OMP THREADPRIVATE(type_ecri_files, phys_out_filetypes)
  CHARACTER(LEN=20), DIMENSION(nfiles), SAVE  :: phys_out_filenames
  !$OMP THREADPRIVATE(phys_out_filenames)

  ! swaero_diag : flag indicates if it is necessary to do calculation for some aerosol diagnostics
  ! swaerofree_diag : flag indicates if it is necessary to do calculation for some aerosol diagnostics
  ! dryaod_diag : flag indicates if it is necessary to do calculation for some aerosol diagnostics
  !--OB: this needs to be set to TRUE by default and changed back to FALSE after first radiation call
  !--    and corrected back to TRUE based on output requests
  LOGICAL, SAVE                                :: swaerofree_diag=.TRUE.
  LOGICAL, SAVE                                :: swaero_diag=.TRUE.
  LOGICAL, SAVE                                :: dryaod_diag=.TRUE.
  !$OMP THREADPRIVATE(swaerofree_diag, swaero_diag, dryaod_diag)
  ! ok_4xCO2atm : flag indicates if it is necessary to do a second call of
  ! radiation code with a 4xCO2 or another different GES to assess SW/LW
  ! in this case
  !--IM: as for swaero_diag or dryaod_diag this needs to be set to TRUE by default and
  !--    changed back to FALSE after first radiation call and corrected back to TRUE 
  !--    based on output requests
  LOGICAL, SAVE                                :: ok_4xCO2atm=.TRUE.
  !$OMP THREADPRIVATE(ok_4xCO2atm)

  INTEGER, SAVE:: levmin(nfiles) = 1
  INTEGER, SAVE:: levmax(nfiles)
  !$OMP THREADPRIVATE(levmin, levmax)

  REAL, SAVE                :: zdtime_moy
  !$OMP THREADPRIVATE(zdtime_moy)

  LOGICAL, SAVE :: vars_defined = .FALSE. ! ug PAS THREADPRIVATE ET C'EST NORMAL

  REAL, allocatable:: zustr_gwd_hines(:), zvstr_gwd_hines(:) ! (klon)
  REAL, allocatable:: zustr_gwd_front(:), zvstr_gwd_front(:) ! (klon)
  REAL, allocatable:: zustr_gwd_rando(:), zvstr_gwd_rando(:) ! (klon)
  !$OMP THREADPRIVATE(zustr_gwd_hines, zvstr_gwd_hines)
  !$OMP THREADPRIVATE(zustr_gwd_front, zvstr_gwd_front)
  !$OMP THREADPRIVATE(zustr_gwd_rando, zvstr_gwd_rando)

  TYPE ctrl_out
     INTEGER,DIMENSION(nfiles)            :: flag
     CHARACTER(len=20)                    :: name
     CHARACTER(len=150)                   :: description
     CHARACTER(len=20)                    :: unit
     CHARACTER(len=20),DIMENSION(nfiles)  :: type_ecrit
  END TYPE ctrl_out

  REAL, SAVE, ALLOCATABLE :: sens_prec_liq_o(:,:), sens_prec_sol_o(:,:)
  REAL, SAVE, ALLOCATABLE :: lat_prec_liq_o(:,:), lat_prec_sol_o(:,:)
 !$OMP THREADPRIVATE(sens_prec_liq_o, sens_prec_sol_o,lat_prec_liq_o,lat_prec_sol_o)

CONTAINS

  !======================================================================
  SUBROUTINE phys_output_var_init
    use dimphy

    IMPLICIT NONE

    include "clesphys.h"

    !------------------------------------------------

    allocate(snow_o(klon), zfra_o(klon))
    allocate(sza_o(klon) )
    allocate(itau_con(klon))
    allocate(sens_prec_liq_o(klon,2))
    allocate(sens_prec_sol_o(klon,2))
    allocate(lat_prec_liq_o(klon,2))
    allocate(lat_prec_sol_o(klon,2))
    sens_prec_liq_o = 0.0 ; sens_prec_sol_o = 0.0
    lat_prec_liq_o = 0.0 ; lat_prec_sol_o = 0.0

    allocate (bils_ec(klon),bils_ech(klon),bils_tke(klon),bils_diss(klon),bils_kinetic(klon),bils_enthalp(klon),bils_latent(klon))
    allocate (d_qw_col(klon), d_ql_col(klon), d_qs_col(klon), d_qt_col(klon), d_ek_col(klon), d_h_dair_col(klon) &
  &         , d_h_qw_col(klon), d_h_ql_col(klon), d_h_qs_col(klon), d_h_col(klon))
    d_qw_col=0. ; d_ql_col=0. ; d_qs_col=0. ; d_qt_col=0. ; d_ek_col=0. ; d_h_dair_col =0.
    d_h_qw_col=0. ; d_h_ql_col=0. ; d_h_qs_col=0. ; d_h_col=0.

    ! Outputs used in cloudth_vert
    allocate(cloudth_sth(klon,klev))
    allocate(cloudth_senv(klon,klev))
    cloudth_sth = 0. ; cloudth_senv = 0. 
    allocate(cloudth_sigmath(klon,klev))
    allocate(cloudth_sigmaenv(klon,klev))
    cloudth_sigmath = 0. ; cloudth_sigmaenv = 0.

! Marine
! Variables de sortie simulateur AIRS

!     if (ok_airs) then
      allocate (map_prop_hc(klon),map_prop_hist(klon))
      allocate (alt_tropo(klon))
      allocate (map_emis_hc(klon),map_iwp_hc(klon),map_deltaz_hc(klon))
      allocate (map_pcld_hc(klon),map_tcld_hc(klon))
      allocate (map_emis_hist(klon),map_iwp_hist(klon),map_deltaz_hist(klon))
      allocate (map_rad_hist(klon))
      allocate (map_ntot(klon),map_hc(klon),map_hist(klon))
      allocate (map_Cb(klon),map_ThCi(klon),map_Anv(klon))
      allocate (map_emis_Cb(klon),map_pcld_Cb(klon),map_tcld_Cb(klon))
      allocate (map_emis_ThCi(klon),map_pcld_ThCi(klon),map_tcld_ThCi(klon))
      allocate (map_emis_Anv(klon),map_pcld_Anv(klon),map_tcld_Anv(klon))
!     endif

    IF (ok_hines) allocate(zustr_gwd_hines(klon), zvstr_gwd_hines(klon))
    IF (.not.ok_hines.and.ok_gwd_rando) &
                  allocate(zustr_gwd_front(klon), zvstr_gwd_front(klon))
    IF (ok_gwd_rando) allocate(zustr_gwd_rando(klon), zvstr_gwd_rando(klon))

  END SUBROUTINE phys_output_var_init

  !======================================================================
  SUBROUTINE phys_output_var_end
    USE dimphy
    IMPLICIT NONE

    include "clesphys.h"

    deallocate(snow_o,zfra_o,itau_con)
    deallocate(sza_o)
    deallocate (bils_ec,bils_ech,bils_tke,bils_diss,bils_kinetic,bils_enthalp,bils_latent)
    deallocate (d_qw_col, d_ql_col, d_qs_col, d_qt_col, d_ek_col, d_h_dair_col &
  &           , d_h_qw_col, d_h_ql_col, d_h_qs_col, d_h_col)

    ! Outputs used in cloudth_vert
    deallocate(cloudth_sth)
    deallocate(cloudth_senv)
    deallocate(cloudth_sigmath)
    deallocate(cloudth_sigmaenv)

! Marine
! Variables de sortie simulateur AIRS

 !    if (ok_airs) then
      deallocate (map_prop_hc,map_prop_hist)
      deallocate (alt_tropo)
      deallocate (map_emis_hc,map_iwp_hc,map_deltaz_hc)
      deallocate (map_pcld_hc,map_tcld_hc)
      deallocate (map_emis_hist,map_iwp_hist,map_deltaz_hist)
      deallocate (map_rad_hist)
      deallocate (map_ntot,map_hc,map_hist)
      deallocate (map_Cb,map_ThCi,map_Anv)
      deallocate (map_emis_Cb,map_pcld_Cb,map_tcld_Cb)
      deallocate (map_emis_ThCi,map_pcld_ThCi,map_tcld_ThCi)
      deallocate (map_emis_Anv,map_pcld_Anv,map_tcld_Anv)
  !   endif

  END SUBROUTINE phys_output_var_end

END MODULE phys_output_var_mod
