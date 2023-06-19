!
! $Id: conf_phys.F90 1668 2012-10-12 10:47:37Z idelkadi $
!
!
!
MODULE conf_phys_m

  IMPLICIT NONE

CONTAINS

  SUBROUTINE conf_phys(ok_journe, ok_mensuel, ok_instan, ok_hf, &
       ok_LES,&
       callstats,&
       solarlong0,seuil_inversion, &
       fact_cldcon, facttemps,ok_newmicro,iflag_radia,&
       iflag_cld_th, &
       iflag_ratqs,ratqsbas,ratqshaut,tau_ratqs, &
       ok_ade, ok_aie, ok_alw, ok_cdnc, ok_volcan, aerosol_couple, chemistry_couple, &
       flag_aerosol, flag_aerosol_strat, flag_aer_feedback, &
       flag_bc_internal_mixture, bl95_b0, bl95_b1,&
       read_climoz, &
       alp_offset)

    USE IOIPSL
    USE surface_data
    USE phys_cal_mod
    USE carbon_cycle_mod,  ONLY: carbon_cycle_tr, carbon_cycle_cpl, carbon_cycle_rad, level_coupling_esm
    USE mod_grid_phy_lmdz, ONLY: klon_glo
    USE print_control_mod, ONLY: lunout
    USE phys_state_var_mod, ONLY: phys_tstep

    INCLUDE "conema3.h"
    INCLUDE "fisrtilp.h"
    INCLUDE "nuage.h"
    INCLUDE "YOMCST.h"
    INCLUDE "YOMCST2.h"
    INCLUDE "thermcell.h"

    !IM : on inclut/initialise les taux de CH4, N2O, CFC11 et CFC12
    INCLUDE "clesphys.h"
    INCLUDE "compbl.h"
    INCLUDE "comsoil.h"
    INCLUDE "YOEGWD.h"
    !
    ! Configuration de la "physique" de LMDZ a l'aide de la fonction
    ! GETIN de IOIPSL
    !
    ! LF 05/2001
    !
    ! type_ocean:      type d'ocean (force, slab, couple)
    ! version_ocean:   version d'ocean (opa8/nemo pour type_ocean=couple ou 
    !                                   sicOBS,sicINT,sicNO pour type_ocean=slab)
    ! ok_veget:   type de modele de vegetation
    ! ok_journe:  sorties journalieres
    ! ok_hf:  sorties haute frequence
    ! ok_mensuel: sorties mensuelles
    ! ok_instan:  sorties instantanees
    ! ok_ade, ok_aie: apply or not aerosol direct and indirect effects
    ! ok_alw: activate aerosol LW effect
    ! ok_cdnc, ok cloud droplet number concentration
    ! flag_aerosol_strat : flag pour les aerosols stratos
    ! flag_bc_internal_mixture : use BC internal mixture if true
    ! bl95_b*: parameters in the formula to link CDNC to aerosol mass conc 
    ! ok_volcan: activate volcanic diags (SW heat & LW cool rate, SW & LW flux)
    !

    ! Sortie:
    LOGICAL              :: ok_newmicro
    INTEGER              :: iflag_radia
    LOGICAL              :: ok_journe, ok_mensuel, ok_instan, ok_hf
    LOGICAL              :: ok_LES
    LOGICAL              :: callstats
    LOGICAL              :: ok_ade, ok_aie, ok_alw, ok_cdnc, ok_volcan
    LOGICAL              :: aerosol_couple, chemistry_couple
    INTEGER              :: flag_aerosol
    INTEGER              :: flag_aerosol_strat
    LOGICAL              :: flag_aer_feedback
    LOGICAL              :: flag_bc_internal_mixture
    REAL                 :: bl95_b0, bl95_b1
    REAL                 :: fact_cldcon, facttemps,ratqsbas,ratqshaut,tau_ratqs
    INTEGER              :: iflag_cld_th
    INTEGER              :: iflag_ratqs

    CHARACTER (len = 6), SAVE  :: type_ocean_omp, version_ocean_omp, ocean_omp
    CHARACTER (len = 10),SAVE  :: type_veget_omp
    CHARACTER (len = 8), SAVE  :: aer_type_omp
    LOGICAL, SAVE       :: ok_snow_omp
    LOGICAL, SAVE       :: ok_newmicro_omp
    LOGICAL, SAVE       :: ok_all_xml_omp
    LOGICAL, SAVE       :: ok_lwoff_omp
    LOGICAL, SAVE       :: ok_journe_omp, ok_mensuel_omp, ok_instan_omp, ok_hf_omp        
    LOGICAL, SAVE       :: ok_LES_omp   
    LOGICAL, SAVE       :: callstats_omp
    LOGICAL, SAVE       :: ok_ade_omp, ok_aie_omp, ok_alw_omp, ok_cdnc_omp, ok_volcan_omp
    LOGICAL, SAVE       :: aerosol_couple_omp, chemistry_couple_omp
    INTEGER, SAVE       :: flag_aerosol_omp
    INTEGER, SAVE       :: flag_aerosol_strat_omp
    LOGICAL, SAVE       :: flag_aer_feedback_omp
    LOGICAL, SAVE       :: flag_bc_internal_mixture_omp
    REAL,SAVE           :: bl95_b0_omp, bl95_b1_omp
    REAL,SAVE           :: freq_ISCCP_omp, ecrit_ISCCP_omp
    REAL,SAVE           :: freq_COSP_omp, freq_AIRS_omp 
    REAL,SAVE           :: fact_cldcon_omp, facttemps_omp,ratqsbas_omp
    REAL,SAVE           :: tau_cld_cv_omp, coefw_cld_cv_omp
    INTEGER, SAVE       :: iflag_cld_cv_omp

    REAL, SAVE          :: ratqshaut_omp
    REAL, SAVE          :: tau_ratqs_omp
    REAL, SAVE          :: t_coupl_omp
    INTEGER, SAVE       :: iflag_radia_omp
    INTEGER, SAVE       :: iflag_rrtm_omp
    INTEGER, SAVE       :: iflag_albedo_omp !albedo SB
    LOGICAL, SAVE       :: ok_chlorophyll_omp ! albedo SB  
    INTEGER, SAVE       :: NSW_omp
    INTEGER, SAVE       :: iflag_cld_th_omp, ip_ebil_phy_omp
    INTEGER, SAVE       :: iflag_ratqs_omp

    REAL, SAVE          :: f_cdrag_ter_omp,f_cdrag_oce_omp
    REAL, SAVE          :: f_rugoro_omp   , z0min_omp
    REAL, SAVE          :: z0m_seaice_omp,z0h_seaice_omp
    REAL, SAVE          :: min_wind_speed_omp,f_gust_wk_omp,f_gust_bl_omp,f_qsat_oce_omp, f_z0qh_oce_omp
    INTEGER, SAVE       :: iflag_gusts_omp,iflag_z0_oce_omp

    ! Local
    REAL                 :: zzz

    REAL :: seuil_inversion
    REAL,SAVE :: seuil_inversion_omp

    INTEGER,SAVE :: iflag_thermals_ed_omp,iflag_thermals_optflux_omp,iflag_thermals_closure_omp
    REAL, SAVE :: fact_thermals_ed_dz_omp
    INTEGER,SAVE :: iflag_thermals_omp,nsplit_thermals_omp
    REAL,SAVE :: tau_thermals_omp,alp_bl_k_omp
    ! nrlmd le 10/04/2012
    INTEGER,SAVE :: iflag_trig_bl_omp,iflag_clos_bl_omp
    INTEGER,SAVE :: tau_trig_shallow_omp,tau_trig_deep_omp
    REAL,SAVE    :: s_trig_omp
    ! fin nrlmd le 10/04/2012
    REAL :: alp_offset
    REAL, SAVE :: alp_offset_omp
    INTEGER,SAVE :: iflag_coupl_omp,iflag_clos_omp,iflag_wake_omp
    INTEGER,SAVE :: iflag_cvl_sigd_omp
    REAL, SAVE :: coef_clos_ls_omp
    REAL, SAVE :: supcrit1_omp, supcrit2_omp
    INTEGER, SAVE :: iflag_mix_omp
    INTEGER, SAVE :: iflag_mix_adiab_omp
    REAL, SAVE :: scut_omp, qqa1_omp, qqa2_omp, gammas_omp, Fmax_omp, alphas_omp
    REAL, SAVE :: tmax_fonte_cv_omp

    REAL,SAVE :: R_ecc_omp,R_peri_omp,R_incl_omp,solaire_omp
    REAL,SAVE      :: solaire_omp_init
    LOGICAL,SAVE :: ok_suntime_rrtm_omp
    REAL,SAVE :: co2_ppm_omp, RCO2_omp, co2_ppm_per_omp, RCO2_per_omp
    REAL,SAVE :: CH4_ppb_omp, RCH4_omp, CH4_ppb_per_omp, RCH4_per_omp
    REAL,SAVE :: N2O_ppb_omp, RN2O_omp, N2O_ppb_per_omp, RN2O_per_omp
    REAL,SAVE :: CFC11_ppt_omp,RCFC11_omp,CFC11_ppt_per_omp,RCFC11_per_omp
    REAL,SAVE :: CFC12_ppt_omp,RCFC12_omp,CFC12_ppt_per_omp,RCFC12_per_omp
    REAL,SAVE :: epmax_omp
    REAL,SAVE :: coef_epmax_cape_omp
    LOGICAL,SAVE :: ok_adj_ema_omp
    INTEGER,SAVE :: iflag_clw_omp
    REAL,SAVE :: cld_lc_lsc_omp,cld_lc_con_omp,cld_tau_lsc_omp,cld_tau_con_omp
    REAL,SAVE :: ffallv_lsc_omp, ffallv_con_omp,coef_eva_omp
    LOGICAL,SAVE :: reevap_ice_omp
    INTEGER,SAVE :: iflag_pdf_omp
    INTEGER,SAVE :: iflag_ice_thermo_omp
    INTEGER,SAVE :: iflag_t_glace_omp
    INTEGER,SAVE :: iflag_cloudth_vert_omp
    INTEGER,SAVE :: iflag_rain_incloud_vol_omp
    REAL,SAVE :: rad_froid_omp, rad_chau1_omp, rad_chau2_omp
    REAL,SAVE :: t_glace_min_omp, t_glace_max_omp
    REAL,SAVE :: exposant_glace_omp
    REAL,SAVE :: rei_min_omp, rei_max_omp
    INTEGER,SAVE :: iflag_sic_omp
    REAL,SAVE :: inertie_sol_omp,inertie_sno_omp,inertie_sic_omp
    REAL,SAVE :: inertie_lic_omp
    REAL,SAVE :: qsol0_omp
    REAL,SAVE :: evap0_omp
    REAL,SAVE :: albsno0_omp
    REAL      :: solarlong0
    REAL,SAVE :: solarlong0_omp
    INTEGER,SAVE :: top_height_omp,overlap_omp
    REAL,SAVE :: cdmmax_omp,cdhmax_omp,ksta_omp,ksta_ter_omp,f_ri_cd_min_omp
    LOGICAL,SAVE :: ok_kzmin_omp
    REAL, SAVE   :: pbl_lmixmin_alpha_omp
    REAL, SAVE ::  fmagic_omp, pmagic_omp
    INTEGER,SAVE :: iflag_pbl_omp,lev_histhf_omp,lev_histday_omp,lev_histmth_omp
    INTEGER,SAVE :: iflag_pbl_split_omp
!FC
    INTEGER,SAVE :: ifl_pbltree_omp
    REAL,SAVE :: Cd_frein_omp
!FC
    INTEGER,SAVE :: iflag_order2_sollw_omp
    INTEGER, SAVE :: lev_histins_omp, lev_histLES_omp 
    INTEGER, SAVE :: lev_histdayNMC_omp
    INTEGER, SAVE :: levout_histNMC_omp(3)
    LOGICAL, SAVE :: ok_histNMC_omp(3)
    REAL, SAVE    :: freq_outNMC_omp(3), freq_calNMC_omp(3)
    CHARACTER*4, SAVE :: type_run_omp
    LOGICAL, SAVE :: ok_cosp_omp, ok_airs_omp
    LOGICAL, SAVE :: ok_mensuelCOSP_omp,ok_journeCOSP_omp,ok_hfCOSP_omp
    REAL, SAVE    :: lonmin_ins_omp, lonmax_ins_omp, latmin_ins_omp, latmax_ins_omp
    REAL, SAVE    :: ecrit_hf_omp, ecrit_day_omp, ecrit_mth_omp, ecrit_reg_omp
    REAL, SAVE    :: ecrit_ins_omp
    REAL, SAVE    :: ecrit_LES_omp
    REAL, SAVE    :: ecrit_tra_omp
    REAL, SAVE    :: cvl_comp_threshold_omp
    REAL, SAVE    :: cvl_sig2feed_omp
    REAL, SAVE    :: cvl_corr_omp
    LOGICAL, SAVE :: ok_lic_melt_omp
    LOGICAL, SAVE :: ok_lic_cond_omp
    !
    INTEGER, SAVE :: iflag_cycle_diurne_omp
    LOGICAL, SAVE :: soil_model_omp,new_oliq_omp
    LOGICAL, SAVE :: ok_orodr_omp, ok_orolf_omp, ok_limitvrai_omp
    INTEGER, SAVE :: nbapp_rad_omp, iflag_con_omp
    INTEGER, SAVE :: nbapp_cv_omp, nbapp_wk_omp
    INTEGER, SAVE :: iflag_ener_conserv_omp
    LOGICAL, SAVE :: ok_conserv_q_omp
    INTEGER, SAVE :: iflag_fisrtilp_qsat_omp
    INTEGER, SAVE :: iflag_bergeron_omp
    LOGICAL,SAVE  :: ok_strato_omp
    LOGICAL,SAVE  :: ok_hines_omp, ok_gwd_rando_omp
    REAL, SAVE    :: gwd_rando_ruwmax_omp, gwd_rando_sat_omp
    REAL, SAVE    :: gwd_front_ruwmax_omp, gwd_front_sat_omp
    REAL, SAVE    :: sso_gkdrag_omp,sso_grahil_omp,sso_grcrit_omp
    REAL, SAVE    :: sso_gfrcri_omp,sso_gkwake_omp,sso_gklift_omp
    LOGICAL, SAVE :: ok_qch4_omp
    LOGICAL, SAVE :: carbon_cycle_tr_omp
    LOGICAL, SAVE :: carbon_cycle_cpl_omp
    LOGICAL, SAVE :: carbon_cycle_rad_omp
    INTEGER, SAVE :: level_coupling_esm_omp
    LOGICAL, SAVE :: adjust_tropopause_omp
    LOGICAL, SAVE :: ok_daily_climoz_omp

    INTEGER, INTENT(OUT):: read_climoz ! read ozone climatology, OpenMP shared
    ! Allowed values are 0, 1 and 2
    ! 0: do not read an ozone climatology
    ! 1: read a single ozone climatology that will be used day and night
    ! 2: read two ozone climatologies, the average day and night
    ! climatology and the daylight climatology

    !-----------------------------------------------------------------

    print*,'CONFPHYS ENTREE'
    !$OMP MASTER 
    !Config Key  = type_ocean 
    !Config Desc = Type d'ocean
    !Config Def  = force
    !Config Help = Type d'ocean utilise: force, slab,couple
    !
    type_ocean_omp = 'force '
    CALL getin('type_ocean', type_ocean_omp)
    !
    !Config Key  = version_ocean 
    !Config Desc = Version d'ocean
    !Config Def  = xxxxxx
    !Config Help = Version d'ocean utilise: opa8/nemo/sicOBS/xxxxxx
    !
    version_ocean_omp = 'xxxxxx'
    CALL getin('version_ocean', version_ocean_omp)

    !Config Key  = OCEAN
    !Config Desc = Old parameter name for type_ocean
    !Config Def  = yyyyyy
    !Config Help = This is only for testing purpose
    !
    ocean_omp = 'yyyyyy'
    CALL getin('OCEAN', ocean_omp)
    IF (ocean_omp /= 'yyyyyy') THEN
       WRITE(lunout,*)'ERROR! Old variable name OCEAN used in parmeter file.'
       WRITE(lunout,*)'Variable OCEAN has been replaced by the variable type_ocean.'
       WRITE(lunout,*)'You have to update your parameter file physiq.def to succed running'
       CALL abort_physic('conf_phys','Variable OCEAN no longer existing, use variable name type_ocean',1)
    ENDIF

    !Config Key  = t_coupl
    !Config Desc = Pas de temps du couplage atm/oce en sec.
    !Config Def  = 86400
    !Config Help = This is only for testing purpose
    !
    t_coupl_omp = 86400.
    CALL getin('t_coupl', t_coupl_omp)
    IF (t_coupl_omp == 0) THEN
       WRITE(lunout,*)'ERROR! Timestep of coupling between atmosphere and ocean'
       WRITE(lunout,*)'cannot be zero.'
       CALL abort_physic('conf_phys','t_coupl = 0.',1)
    ENDIF

    !
    !Config Key  = ok_all_xml 
    !Config Desc = utiliser les xml pourles définitions des champs pour xios
    !Config Def  = .FALSE.
    !Config Help = 
    !
    ok_all_xml_omp = .FALSE.
    CALL getin('ok_all_xml', ok_all_xml_omp)

    !
    !Config Key  = ok_lwoff
    !Config Desc = inhiber l effet radiatif LW des nuages
    !Config Def  = .FALSE.
    !Config Help = 
    !
    ok_lwoff_omp = .FALSE.
    CALL getin('ok_lwoff', ok_lwoff_omp)
    !

    !
    !Config Key  = VEGET 
    !Config Desc = Type de modele de vegetation
    !Config Def  = .FALSE.
    !Config Help = Type de modele de vegetation utilise
    !
    type_veget_omp ='orchidee'
    CALL getin('VEGET', type_veget_omp)
    !

    ! Martin
    !Config Key  = ok_snow
    !Config Desc = Flag to activate snow model SISVAT
    !Config Def  = .FALSE.
    ok_snow_omp = .FALSE.
    CALL getin('ok_snow', ok_snow_omp)
    ! Martin

    !Config Key  = OK_journe
    !Config Desc = Pour des sorties journalieres 
    !Config Def  = .FALSE.
    !Config Help = Pour creer le fichier histday contenant les sorties
    !              journalieres 
    !
    ok_journe_omp = .FALSE.
    CALL getin('OK_journe', ok_journe_omp)
    !
    !Config Key  = ok_hf
    !Config Desc = Pour des sorties haute frequence
    !Config Def  = .FALSE.
    !Config Help = Pour creer le fichier histhf contenant les sorties
    !              haute frequence ( 3h ou 6h)
    !
    ok_hf_omp = .FALSE.
    CALL getin('ok_hf', ok_hf_omp)
    !
    !Config Key  = OK_mensuel
    !Config Desc = Pour des sorties mensuelles 
    !Config Def  = .TRUE.
    !Config Help = Pour creer le fichier histmth contenant les sorties
    !              mensuelles 
    !
    ok_mensuel_omp = .TRUE.
    CALL getin('OK_mensuel', ok_mensuel_omp)
    !
    !Config Key  = OK_instan
    !Config Desc = Pour des sorties instantanees 
    !Config Def  = .FALSE.
    !Config Help = Pour creer le fichier histins contenant les sorties
    !              instantanees 
    !
    ok_instan_omp = .FALSE.
    CALL getin('OK_instan', ok_instan_omp)
    !
    !Config Key  = ok_ade
    !Config Desc = Aerosol direct effect or not?
    !Config Def  = .FALSE.
    !Config Help = Used in radlwsw.F
    !
    ok_ade_omp = .FALSE.
    CALL getin('ok_ade', ok_ade_omp)

    !Config Key  = ok_alw
    !Config Desc = Aerosol longwave effect or not?
    !Config Def  = .FALSE.
    !Config Help = Used in radlwsw.F
    !
    ok_alw_omp = .FALSE.
    CALL getin('ok_alw', ok_alw_omp)

    !
    !Config Key  = ok_aie
    !Config Desc = Aerosol indirect effect or not?
    !Config Def  = .FALSE.
    !Config Help = Used in nuage.F and radlwsw.F
    !
    ok_aie_omp = .FALSE.
    CALL getin('ok_aie', ok_aie_omp)

    !
    !Config Key  = ok_cdnc
    !Config Desc = ok cloud droplet number concentration
    !Config Def  = .FALSE.
    !Config Help = Used in newmicro.F
    !
    ok_cdnc_omp = .FALSE.
    CALL getin('ok_cdnc', ok_cdnc_omp)

    !
    !Config Key  = ok_volcan
    !Config Desc = ok to generate volcanic diags
    !Config Def  = .FALSE.
    !Config Help = Used in radlwsw_m.F
    !
    ok_volcan_omp = .FALSE.
    CALL getin('ok_volcan', ok_volcan_omp)

    !
    !Config Key  = aerosol_couple
    !Config Desc = read aerosol in file or calcul by inca
    !Config Def  = .FALSE.
    !Config Help = Used in physiq.F
    !
    aerosol_couple_omp = .FALSE.
    CALL getin('aerosol_couple',aerosol_couple_omp)
    !
    !Config Key  = chemistry_couple
    !Config Desc = read O3 chemistry in file or calcul by inca
    !Config Def  = .FALSE.
    !Config Help = Used in physiq.F
    !
    chemistry_couple_omp = .FALSE.
    CALL getin('chemistry_couple',chemistry_couple_omp)
    !
    !Config Key  = flag_aerosol
    !Config Desc = which aerosol is use for coupled model
    !Config Def  = 1
    !Config Help = Used in physiq.F
    !
    ! - flag_aerosol=0 => no aerosol
    ! - flag_aerosol=1 => so4 only (defaut) 
    ! - flag_aerosol=2 => bc  only 
    ! - flag_aerosol=3 => pom only
    ! - flag_aerosol=4 => seasalt only 
    ! - flag_aerosol=5 => dust only
    ! - flag_aerosol=6 => all aerosol
    ! - flag_aerosol=7 => natural aerosol + MACv2SP
    ! - (in this case aerosols.1980.nc should point to aerosols.nat.nc)

    flag_aerosol_omp = 0
    CALL getin('flag_aerosol',flag_aerosol_omp)

    !
    !Config Key  = flag_bc_internal_mixture
    !Config Desc = state of mixture for BC aerosols
    ! - n = external mixture
    ! - y = internal mixture
    !Config Def  = n
    !Config Help = Used in physiq.F / aeropt
    !
    flag_bc_internal_mixture_omp = .FALSE.
    CALL getin('flag_bc_internal_mixture',flag_bc_internal_mixture_omp)

    ! 
    !Config Key  = aer_type 
    !Config Desc = Use a constant field for the aerosols 
    !Config Def  = scenario 
    !Config Help = Used in readaerosol.F90 
    ! 
    aer_type_omp = 'scenario' 
    CALL getin('aer_type', aer_type_omp) 

    !
    !Config Key  = bl95_b0
    !Config Desc = Parameter in CDNC-maer link (Boucher&Lohmann 1995)
    !Config Def  = .FALSE.
    !Config Help = Used in nuage.F
    !
    bl95_b0_omp = 2.
    CALL getin('bl95_b0', bl95_b0_omp)

    !Config Key  = bl95_b1
    !Config Desc = Parameter in CDNC-maer link (Boucher&Lohmann 1995)
    !Config Def  = .FALSE.
    !Config Help = Used in nuage.F
    !
    bl95_b1_omp = 0.2
    CALL getin('bl95_b1', bl95_b1_omp)

    !Config Key  = freq_ISCCP
    !Config Desc = Frequence d'appel du simulateur ISCCP en secondes;
    !              par defaut 10800, i.e. 3 heures 
    !Config Def  = 10800.
    !Config Help = Used in ini_histISCCP.h
    !
    freq_ISCCP_omp = 10800.
    CALL getin('freq_ISCCP', freq_ISCCP_omp)
    !
    !Config Key  = ecrit_ISCCP
    !Config Desc = Frequence d'ecriture des resultats du simulateur ISCCP en nombre de jours;
    !              par defaut 1., i.e. 1 jour
    !Config Def  = 1.
    !Config Help = Used in ini_histISCCP.h
    !
    !
    ecrit_ISCCP_omp = 1.
    CALL getin('ecrit_ISCCP', ecrit_ISCCP_omp)

    !Config Key  = freq_COSP
    !Config Desc = Frequence d'appel du simulateur COSP en secondes;
    !              par defaut 10800, i.e. 3 heures
    !Config Def  = 10800.
    !Config Help = Used in ini_histdayCOSP.h
    !
    freq_COSP_omp = 10800.
    CALL getin('freq_COSP', freq_COSP_omp)

    !Config Key  = freq_AIRS
    !Config Desc = Frequence d'appel du simulateur AIRS en secondes;
    !              par defaut 10800, i.e. 3 heures
    !Config Def  = 10800.
    !Config Help = Used in ini_histdayAIRS.h
    !
    freq_AIRS_omp = 10800.
    CALL getin('freq_AIRS', freq_AIRS_omp)

    !
    !Config Key  = ip_ebil_phy
    !Config Desc = Niveau de sortie pour les diags bilan d'energie 
    !Config Def  = 0
    !Config Help = 
    !               
    ip_ebil_phy_omp = 0
    CALL getin('ip_ebil_phy', ip_ebil_phy_omp)
    IF (ip_ebil_phy_omp/=0) THEN
       CALL abort_physic('conf_phys','ip_ebil_phy_omp doit etre 0 sur cette version',1)
    ENDIF

    !
    !Config Key  = seuil_inversion
    !Config Desc = Seuil ur dTh pour le choix entre les schemas de CL
    !Config Def  = -0.1
    !Config Help = 
    !               
    seuil_inversion_omp = -0.1
    CALL getin('seuil_inversion', seuil_inversion_omp)

    !
    ! Constante solaire & Parametres orbitaux & taux gaz effet de serre BEG
    !
    !Config Key  = R_ecc
    !Config Desc = Excentricite
    !Config Def  = 0.016715
    !Config Help = 
    !               
    !valeur AMIP II
    R_ecc_omp = 0.016715
    CALL getin('R_ecc', R_ecc_omp)
    !
    !Config Key  = R_peri
    !Config Desc = Equinoxe
    !Config Def  = 
    !Config Help = 
    !               
    !
    !valeur AMIP II
    R_peri_omp = 102.7
    CALL getin('R_peri', R_peri_omp)
    !
    !Config Key  = R_incl
    !Config Desc = Inclinaison
    !Config Def  = 
    !Config Help = 
    !               
    !
    !valeur AMIP II
    R_incl_omp = 23.441
    CALL getin('R_incl', R_incl_omp)
    !
    !Config Key  = solaire
    !Config Desc = Constante solaire en W/m2
    !Config Def  = 1365.
    !Config Help = 
    !               
    !
    !valeur AMIP II
    solaire_omp = 1365.
    solaire_omp_init = solaire_omp     !--we keep track of the default value
    CALL getin('solaire', solaire_omp)
    !
    !Config Key  = co2_ppm
    !Config Desc = concentration du gaz carbonique en ppmv
    !Config Def  = 348.
    !Config Help = 
    !               
    !
    !valeur AMIP II
    co2_ppm_omp = 348.
    CALL getin('co2_ppm', co2_ppm_omp)
    !
    !Config Key  = RCO2
    !Config Desc = Concentration du CO2
    !Config Def  = co2_ppm * 1.0e-06  * 44.011/28.97
    !Config Def  = 348. * 1.0e-06  * 44.011/28.97
    !Config Help = 
    !               
    ! RCO2 = 5.286789092164308E-04
    !ancienne valeur
    RCO2_omp = co2_ppm_omp * 1.0e-06 * RMCO2 / RMD ! pour co2_ppm=348.

    !  CALL getin('RCO2', RCO2)
    !
    !Config Key  = RCH4
    !Config Desc = Concentration du CH4
    !Config Def  = 1.65E-06* 16.043/28.97
    !Config Help = 
    !               
    !
    !valeur AMIP II
    !OK  RCH4 = 1.65E-06* 16.043/28.97
    ! RCH4 = 9.137366240938903E-07
    !
    !ancienne valeur
    ! RCH4 = 1.72E-06* 16.043/28.97
    !OK CALL getin('RCH4', RCH4)
    zzz = 1650.
    CALL getin('CH4_ppb', zzz)
    CH4_ppb_omp = zzz
    RCH4_omp = CH4_ppb_omp * 1.0E-09 * RMCH4 / RMD
    !
    !Config Key  = RN2O
    !Config Desc = Concentration du N2O
    !Config Def  = 306.E-09* 44.013/28.97
    !Config Help = 
    !               
    !
    !valeur AMIP II
    !OK  RN2O = 306.E-09* 44.013/28.97
    ! RN2O = 4.648939592682085E-07
    !
    !ancienne valeur
    ! RN2O = 310.E-09* 44.013/28.97
    !OK  CALL getin('RN2O', RN2O)
    zzz=306.
    CALL getin('N2O_ppb', zzz)
    N2O_ppb_omp = zzz
    RN2O_omp = N2O_ppb_omp * 1.0E-09 * RMN2O / RMD
    !
    !Config Key  = RCFC11
    !Config Desc = Concentration du CFC11
    !Config Def  = 280.E-12* 137.3686/28.97
    !Config Help = 
    !               
    !
    !OK RCFC11 = 280.E-12* 137.3686/28.97
    zzz = 280.
    CALL getin('CFC11_ppt',zzz)
    CFC11_ppt_omp = zzz
    RCFC11_omp=CFC11_ppt_omp* 1.0E-12 * RMCFC11 / RMD
    ! RCFC11 = 1.327690990680013E-09
    !OK CALL getin('RCFC11', RCFC11)
    !
    !Config Key  = RCFC12
    !Config Desc = Concentration du CFC12
    !Config Def  = 484.E-12* 120.9140/28.97
    !Config Help = 
    !               
    !
    !OK RCFC12 = 484.E-12* 120.9140/28.97
    zzz = 484.
    CALL getin('CFC12_ppt',zzz)
    CFC12_ppt_omp = zzz
    RCFC12_omp = CFC12_ppt_omp * 1.0E-12 * RMCFC12 / RMD
    ! RCFC12 = 2.020102726958923E-09
    !OK CALL getin('RCFC12', RCFC12)

    !ajout CFMIP begin
    !
    !Config Key  = co2_ppm_per
    !Config Desc = concentration du co2_ppm_per
    !Config Def  = 348.
    !Config Help = 
    !               
    co2_ppm_per_omp = co2_ppm_omp
    CALL getin('co2_ppm_per', co2_ppm_per_omp)
    !
    !Config Key  = RCO2_per
    !Config Desc = Concentration du CO2_per
    !Config Def  = co2_ppm_per * 1.0e-06  * 44.011/28.97
    !Config Def  = 348. * 1.0e-06  * 44.011/28.97
    !Config Help = 
    !               
    RCO2_per_omp = co2_ppm_per_omp * 1.0e-06 * RMCO2 / RMD

    !Config Key  = ok_4xCO2atm
    !Config Desc = Calcul ou non effet radiatif 4xco2
    !Config Def  = .FALSE.
    !Config Help = 

    !Config Key  = RCH4_per
    !Config Desc = Concentration du CH4_per
    !Config Def  = 1.65E-06* 16.043/28.97
    !Config Help = 
    !               
    zzz = CH4_ppb_omp
    CALL getin('CH4_ppb_per', zzz)
    CH4_ppb_per_omp = zzz
    RCH4_per_omp = CH4_ppb_per_omp * 1.0E-09 * RMCH4 / RMD
    !
    !Config Key  = RN2O_per
    !Config Desc = Concentration du N2O_per
    !Config Def  = 306.E-09* 44.013/28.97
    !Config Help = 
    !               
    zzz = N2O_ppb_omp
    CALL getin('N2O_ppb_per', zzz)
    N2O_ppb_per_omp = zzz
    RN2O_per_omp = N2O_ppb_per_omp * 1.0E-09 * RMN2O / RMD
    !
    !Config Key  = RCFC11_per
    !Config Desc = Concentration du CFC11_per
    !Config Def  = 280.E-12* 137.3686/28.97
    !Config Help = 
    !               
    zzz = CFC11_ppt_omp
    CALL getin('CFC11_ppt_per',zzz)
    CFC11_ppt_per_omp = zzz
    RCFC11_per_omp=CFC11_ppt_per_omp* 1.0E-12 * RMCFC11 / RMD
    !
    !Config Key  = RCFC12_per
    !Config Desc = Concentration du CFC12_per
    !Config Def  = 484.E-12* 120.9140/28.97
    !Config Help = 
    !               
    zzz = CFC12_ppt_omp
    CALL getin('CFC12_ppt_per',zzz)
    CFC12_ppt_per_omp = zzz
    RCFC12_per_omp = CFC12_ppt_per_omp * 1.0E-12 * RMCFC12 / RMD
    !ajout CFMIP end

    !
    ! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
    ! Constantes precedemment dans dyn3d/conf_gcm

    !Config  Key  = iflag_cycle_diurne
    !Config  Desc = Cycle diurne
    !Config  Def  = 1
    !Config  Help = Cette option permet d'eteidre le cycle diurne.
    !Config         Peut etre util pour accelerer le code !
    iflag_cycle_diurne_omp = 1
    CALL getin('iflag_cycle_diurne',iflag_cycle_diurne_omp)

    !Config  Key  = soil_model
    !Config  Desc = Modele de sol
    !Config  Def  = y
    !Config  Help = Choix du modele de sol (Thermique ?)
    !Config         Option qui pourait un string afin de pouvoir
    !Config         plus de choix ! Ou meme une liste d'options !
    soil_model_omp = .TRUE.
    CALL getin('soil_model',soil_model_omp)

    !Config  Key  = new_oliq
    !Config  Desc = Nouvelle eau liquide
    !Config  Def  = y
    !Config  Help = Permet de mettre en route la
    !Config         nouvelle parametrisation de l'eau liquide !
    new_oliq_omp = .TRUE.
    CALL getin('new_oliq',new_oliq_omp)

    !Config  Key  = ok_orodr
    !Config  Desc = Orodr ???
    !Config  Def  = y
    !Config  Help = Y en a pas comprendre !
    !Config         
    ok_orodr_omp = .TRUE.
    CALL getin('ok_orodr',ok_orodr_omp)

    !Config  Key  =  ok_orolf
    !Config  Desc = Orolf ??
    !Config  Def  = y
    !Config  Help = Connais pas !
    ok_orolf_omp = .TRUE.
    CALL getin('ok_orolf', ok_orolf_omp)

    !Config  Key  = ok_limitvrai
    !Config  Desc = Force la lecture de la bonne annee
    !Config  Def  = n
    !Config  Help = On peut forcer le modele a lire le
    !Config         fichier SST de la bonne annee. C'est une tres bonne
    !Config         idee, pourquoi ne pas mettre toujours a y ???
    ok_limitvrai_omp = .FALSE.
    CALL getin('ok_limitvrai',ok_limitvrai_omp)

    !Config  Key  = nbapp_rad
    !Config  Desc = Frequence d'appel au rayonnement
    !Config  Def  = 12
    !Config  Help = Nombre  d'appels des routines de rayonnements
    !Config         par jour.
    nbapp_rad_omp = 12
    CALL getin('nbapp_rad',nbapp_rad_omp)

    !Config  Key  = iflag_con
    !Config  Desc = Flag de convection
    !Config  Def  = 2
    !Config  Help = Flag  pour la convection les options suivantes existent :
    !Config         1 pour LMD,
    !Config         2 pour Tiedtke,
    !Config         3 pour CCM(NCAR)  
    iflag_con_omp = 2
    CALL getin('iflag_con',iflag_con_omp)

    !Config  Key  = nbapp_cv
    !Config  Desc = Frequence d'appel a la convection
    !Config  Def  = 0
    !Config  Help = Nombre  d'appels des routines de convection
    !Config         par jour. Si =0, appel a chaque pas de temps physique.
    nbapp_cv_omp = 0
    CALL getin('nbapp_cv',nbapp_cv_omp)

    !Config  Key  = nbapp_wk
    !Config  Desc = Frequence d'appel aux wakes
    !Config  Def  = 0
    !Config  Help = Nombre  d'appels des routines de wakes
    !Config         par jour. Si =0, appel a chaque pas de temps physique.
    nbapp_wk_omp = 0
    CALL getin('nbapp_wk',nbapp_wk_omp)

    !Config  Key  = iflag_ener_conserv
    !Config  Desc = Flag de convection
    !Config  Def  = 1
    !Config  Help = Flag  pour la convection les options suivantes existent :
    !Config         -1 pour Kinetic energy correction
    !Config         1  conservation kinetic and enthalpy
    iflag_ener_conserv_omp = -1
    CALL getin('iflag_ener_conserv',iflag_ener_conserv_omp)

    !Config  Key  = ok_conserv_q
    !Config  Desc = Switch des corrections de conservation de l'eau
    !Config  Def  = y
    !Config  Help = Switch des corrections de conservation de l'eau
    !Config         y -> corrections activees
    !Config         n -> conformite avec versions anterieures au 1/4/2014
    ok_conserv_q_omp = .FALSE.
    CALL getin('ok_conserv_q',ok_conserv_q_omp)

    !Config  Key  = iflag_fisrtilp_qsat
    !Config  Desc = Flag de fisrtilp
    !Config  Def  = 0
    !Config  Help = Flag  pour la pluie grande-échelle les options suivantes existent :
    !Config         >1 nb iterations pour converger dans le calcul de qsat
    iflag_fisrtilp_qsat_omp = 0
    CALL getin('iflag_fisrtilp_qsat',iflag_fisrtilp_qsat_omp)

    !Config  Key  = iflag_bergeron
    !Config  Desc = Flag de fisrtilp
    !Config  Def  = 0
    !Config  Help = Flag  pour la pluie grande-échelle les options suivantes existent :
    !Config         0 pas d effet Bergeron
    !Config         1 effet Bergeron pour T<0
    iflag_bergeron_omp = 0
    CALL getin('iflag_bergeron',iflag_bergeron_omp)

    !
    !
    !
    ! Constante solaire & Parametres orbitaux & taux gaz effet de serre END
    !
    ! KE
    !

    !Config key  = cvl_comp_threshold
    !Config Desc = maximum fraction of convective points enabling compression
    !Config Def  = 1.00
    !Config Help = fields are compressed when less than a fraction cvl_comp_threshold
    !Config Help = of the points is convective.
    cvl_comp_threshold_omp = 1.00
    CALL getin('cvl_comp_threshold', cvl_comp_threshold_omp)

    !Config key  = cvl_sig2feed
    !Config Desc = sigma coordinate at top of feeding layer
    !Config Def  = 0.97
    !Config Help = deep convection is fed by the layer extending from the surface (pressure ps)
    !Config Help = and cvl_sig2feed*ps.
    cvl_sig2feed_omp = 0.97
    CALL getin('cvl_sig2feed', cvl_sig2feed_omp)

    !Config key  = cvl_corr
    !Config Desc = Facteur multiplication des precip convectives dans KE
    !Config Def  = 1.00
    !Config Help = 1.02 pour un moderne ou un pre-ind. A ajuster pour un glaciaire
    cvl_corr_omp = 1.00
    CALL getin('cvl_corr', cvl_corr_omp)


    !Config Key  = epmax
    !Config Desc = Efficacite precip
    !Config Def  = 0.993
    !Config Help = 
    !
    epmax_omp = .993
    CALL getin('epmax', epmax_omp)

    coef_epmax_cape_omp = 0.0   
    CALL getin('coef_epmax_cape', coef_epmax_cape_omp)        
    !
    !Config Key  = ok_adj_ema
    !Config Desc =  
    !Config Def  = FALSE
    !Config Help = 
    !
    ok_adj_ema_omp = .FALSE.
    CALL getin('ok_adj_ema',ok_adj_ema_omp)
    !
    !Config Key  = iflag_clw
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_clw_omp = 0
    CALL getin('iflag_clw',iflag_clw_omp)
    !
    !Config Key  = cld_lc_lsc 
    !Config Desc =  
    !Config Def  = 2.6e-4
    !Config Help = 
    !
    cld_lc_lsc_omp = 2.6e-4
    CALL getin('cld_lc_lsc',cld_lc_lsc_omp)
    !
    !Config Key  = cld_lc_con
    !Config Desc =  
    !Config Def  = 2.6e-4
    !Config Help = 
    !
    cld_lc_con_omp = 2.6e-4
    CALL getin('cld_lc_con',cld_lc_con_omp)
    !
    !Config Key  = cld_tau_lsc
    !Config Desc =  
    !Config Def  = 3600.
    !Config Help = 
    !
    cld_tau_lsc_omp = 3600.
    CALL getin('cld_tau_lsc',cld_tau_lsc_omp)
    !
    !Config Key  = cld_tau_con
    !Config Desc =  
    !Config Def  = 3600.
    !Config Help = 
    !
    cld_tau_con_omp = 3600.
    CALL getin('cld_tau_con',cld_tau_con_omp)
    !
    !Config Key  = ffallv_lsc
    !Config Desc =  
    !Config Def  = 1.
    !Config Help = 
    !
    ffallv_lsc_omp = 1.
    CALL getin('ffallv_lsc',ffallv_lsc_omp)
    !
    !Config Key  = ffallv_con
    !Config Desc =  
    !Config Def  = 1.
    !Config Help = 
    !
    ffallv_con_omp = 1.
    CALL getin('ffallv_con',ffallv_con_omp)
    !
    !Config Key  = coef_eva
    !Config Desc =  
    !Config Def  = 2.e-5
    !Config Help = 
    !
    coef_eva_omp = 2.e-5
    CALL getin('coef_eva',coef_eva_omp)
    !
    !Config Key  = reevap_ice
    !Config Desc =  
    !Config Def  = .FALSE.
    !Config Help = 
    !
    reevap_ice_omp = .FALSE.
    CALL getin('reevap_ice',reevap_ice_omp)

    !Config Key  = iflag_ratqs
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    iflag_ratqs_omp = 1
    CALL getin('iflag_ratqs',iflag_ratqs_omp)

    !
    !Config Key  = iflag_radia 
    !Config Desc =  
    !Config Def  = 1
    !Config Help = 
    !
    iflag_radia_omp = 1
    CALL getin('iflag_radia',iflag_radia_omp)

    !
    !Config Key  = iflag_rrtm 
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_rrtm_omp = 0
    CALL getin('iflag_rrtm',iflag_rrtm_omp)

    !
    !Config Key  = NSW 
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    NSW_omp = 2
    CALL getin('NSW',NSW_omp)
    !albedo SB >>>
    iflag_albedo_omp = 0
    CALL getin('iflag_albedo',iflag_albedo_omp)

    ok_chlorophyll_omp=.FALSE.
    CALL getin('ok_chlorophyll',ok_chlorophyll_omp)
    !albedo SB <<<
    !
    !Config Key  = ok_sun_time
    !Config Desc = oui ou non variabilite solaire
    !Config Def  = .FALSE.
    !Config Help =
    !
    !
    !valeur AMIP II
    ok_suntime_rrtm_omp = .FALSE.
    IF (iflag_rrtm_omp==1) THEN
      CALL getin('ok_suntime_rrtm',ok_suntime_rrtm_omp)
    ENDIF
    
    !Config Key  = flag_aerosol_strat
    !Config Desc = use stratospheric aerosols 0, 1, 2
    ! - 0 = no stratospheric aerosols 
    ! - 1 = stratospheric aerosols scaled from 550 nm AOD 
    ! - 2 = stratospheric aerosol properties from CMIP6
    !Option 2 is only available with RRTM, this is tested later on 
    !Config Def  = 0
    !Config Help = Used in physiq.F
    !
    flag_aerosol_strat_omp = 0
    CALL getin('flag_aerosol_strat',flag_aerosol_strat_omp)

    !Config Key  = flag_aer_feedback
    !Config Desc = (des)activate aerosol radiative feedback
    ! - F = no aerosol radiative feedback
    ! - T = aerosol radiative feedback
    !Config Def  = T
    !Config Help = Used in physiq.F
    !
    flag_aer_feedback_omp = .TRUE. 
    IF (iflag_rrtm_omp==1) THEN
       CALL getin('flag_aer_feedback',flag_aer_feedback_omp)
    ENDIF

    !Config Key  = iflag_cld_th 
    !Config Desc =  
    !Config Def  = 1
    !Config Help = 
    !
    iflag_cld_th_omp = 1
    ! On lit deux fois avec l'ancien et le nouveau nom 
    ! pour assurer une retrocompatiblite.
    ! A abandonner un jour
    CALL getin('iflag_cldcon',iflag_cld_th_omp)
    CALL getin('iflag_cld_th',iflag_cld_th_omp)
    iflag_cld_cv_omp = 0 
    CALL getin('iflag_cld_cv',iflag_cld_cv_omp)

    !
    !Config Key  = tau_cld_cv
    !Config Desc =
    !Config Def  = 10.
    !Config Help =
    !
    tau_cld_cv_omp = 10.
    CALL getin('tau_cld_cv',tau_cld_cv_omp)

    !
    !Config Key  = coefw_cld_cv
    !Config Desc =
    !Config Def  = 0.1
    !Config Help =
    !
    coefw_cld_cv_omp = 0.1
    CALL getin('coefw_cld_cv',coefw_cld_cv_omp)




    !
    !Config Key  = iflag_pdf 
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_pdf_omp = 0
    CALL getin('iflag_pdf',iflag_pdf_omp)
    !
    !Config Key  = fact_cldcon
    !Config Desc =  
    !Config Def  = 0.375
    !Config Help = 
    !
    fact_cldcon_omp = 0.375
    CALL getin('fact_cldcon',fact_cldcon_omp)

    !
    !Config Key  = facttemps
    !Config Desc =  
    !Config Def  = 1.e-4
    !Config Help = 
    !
    facttemps_omp = 1.e-4
    CALL getin('facttemps',facttemps_omp)

    !
    !Config Key  = ok_newmicro
    !Config Desc =  
    !Config Def  = .TRUE.
    !Config Help = 
    !
    ok_newmicro_omp = .TRUE.
    CALL getin('ok_newmicro',ok_newmicro_omp)
    !
    !Config Key  = ratqsbas
    !Config Desc =  
    !Config Def  = 0.01
    !Config Help = 
    !
    ratqsbas_omp = 0.01
    CALL getin('ratqsbas',ratqsbas_omp)
    !
    !Config Key  = ratqshaut
    !Config Desc =  
    !Config Def  = 0.3
    !Config Help = 
    !
    ratqshaut_omp = 0.3
    CALL getin('ratqshaut',ratqshaut_omp)

    !Config Key  = tau_ratqs
    !Config Desc =  
    !Config Def  = 1800.
    !Config Help = 
    !
    tau_ratqs_omp = 1800.
    CALL getin('tau_ratqs',tau_ratqs_omp)

    !
    !-----------------------------------------------------------------------
    ! Longitude solaire pour le calcul de l'ensoleillement en degre
    ! si on veut imposer la saison. Sinon, solarlong0=-999.999
    !Config Key  = solarlong0
    !Config Desc =  
    !Config Def  = -999.999 
    !Config Help = 
    !
    solarlong0_omp = -999.999
    CALL getin('solarlong0',solarlong0_omp)
    !
    !-----------------------------------------------------------------------
    !  Valeur imposee pour configuration idealisees
    !Config Key  = qsol0 pour le bucket, evap0 pour aquaplanetes, albsno0
    ! Default value -1 to activate the full computation
    qsol0_omp = -1.
    CALL getin('qsol0',qsol0_omp)
    evap0_omp = -1.
    CALL getin('evap0',evap0_omp)
    albsno0_omp = -1.
    CALL getin('albsno0',albsno0_omp)
    !
    !-----------------------------------------------------------------------
    !
    !Config Key  = iflag_sic
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_sic_omp = 0
    CALL getin('iflag_sic',iflag_sic_omp)
    !
    !Config Key  = inertie_sic
    !Config Desc =  
    !Config Def  = 2000.
    !Config Help = 
    !
    inertie_sic_omp = 2000.
    CALL getin('inertie_sic',inertie_sic_omp)
    !
    !Config Key  = inertie_lic
    !Config Desc =  
    !Config Def  = 2000.
    !Config Help = 
    !
    inertie_lic_omp = 2000.
    CALL getin('inertie_lic',inertie_lic_omp)
    !
    !Config Key  = inertie_sno
    !Config Desc =  
    !Config Def  = 2000.
    !Config Help = 
    !
    inertie_sno_omp = 2000.
    CALL getin('inertie_sno',inertie_sno_omp)
    !
    !Config Key  = inertie_sol
    !Config Desc =  
    !Config Def  = 2000.
    !Config Help = 
    !
    inertie_sol_omp = 2000.
    CALL getin('inertie_sol',inertie_sol_omp)

    !
    !Config Key  = rad_froid
    !Config Desc =  
    !Config Def  = 35.0
    !Config Help = 
    !
    rad_froid_omp = 35.0
    CALL getin('rad_froid',rad_froid_omp)

    !
    !Config Key  = rad_chau1
    !Config Desc =  
    !Config Def  = 13.0
    !Config Help = 
    !
    rad_chau1_omp = 13.0
    CALL getin('rad_chau1',rad_chau1_omp)

    !
    !Config Key  = rad_chau2
    !Config Desc =  
    !Config Def  = 9.0
    !Config Help = 
    !
    rad_chau2_omp = 9.0
    CALL getin('rad_chau2',rad_chau2_omp)

    !
    !Config Key  = t_glace_min
    !Config Desc =  
    !Config Def  = 258.
    !Config Help = 
    !
    t_glace_min_omp = 258.
    CALL getin('t_glace_min',t_glace_min_omp)

    !
    !Config Key  = t_glace_max
    !Config Desc =  
    !Config Def  = 273.13
    !Config Help = 
    !
    t_glace_max_omp = 273.13
    CALL getin('t_glace_max',t_glace_max_omp)

    !
    !Config Key  = exposant_glace
    !Config Desc =  
    !Config Def  = 2.
    !Config Help = 
    !
    exposant_glace_omp = 1.
    CALL getin('exposant_glace',exposant_glace_omp)

    !
    !Config Key  = iflag_t_glace
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_t_glace_omp = 0
    CALL getin('iflag_t_glace',iflag_t_glace_omp)

    !
    !Config Key  = iflag_cloudth_vert
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_cloudth_vert_omp = 0
    CALL getin('iflag_cloudth_vert',iflag_cloudth_vert_omp)

    !
    !Config Key  = iflag_rain_incloud_vol
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_rain_incloud_vol_omp = 0
    CALL getin('iflag_rain_incloud_vol',iflag_rain_incloud_vol_omp)

    !
    !Config Key  = iflag_ice_thermo
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_ice_thermo_omp = 0
    CALL getin('iflag_ice_thermo',iflag_ice_thermo_omp)

    !Config Key  = rei_min
    !Config Desc =  
    !Config Def  = 3.5
    !Config Help = 
    !
    rei_min_omp = 3.5
    CALL getin('rei_min',rei_min_omp)

    !
    !Config Key  = rei_max
    !Config Desc =  
    !Config Def  = 61.29
    !Config Help = 
    !
    rei_max_omp = 61.29
    CALL getin('rei_max',rei_max_omp)

    !
    !Config Key  = top_height
    !Config Desc =
    !Config Def  = 3
    !Config Help =
    !
    top_height_omp = 3
    CALL getin('top_height',top_height_omp)

    !
    !Config Key  = overlap
    !Config Desc =
    !Config Def  = 3
    !Config Help =
    !
    overlap_omp = 3
    CALL getin('overlap',overlap_omp)

    !
    !Config Key  = cdmmax
    !Config Desc =
    !Config Def  = 1.3E-3
    !Config Help =
    !
    cdmmax_omp = 1.3E-3
    CALL getin('cdmmax',cdmmax_omp)

    !
    !Config Key  = cdhmax
    !Config Desc =
    !Config Def  = 1.1E-3
    !Config Help =
    !
    cdhmax_omp = 1.1E-3
    CALL getin('cdhmax',cdhmax_omp)

    !261103
    !
    !Config Key  = ksta
    !Config Desc =
    !Config Def  = 1.0e-10
    !Config Help =
    !
    ksta_omp = 1.0e-10
    CALL getin('ksta',ksta_omp)

    !
    !Config Key  = ksta_ter
    !Config Desc =
    !Config Def  = 1.0e-10
    !Config Help =
    !
    ksta_ter_omp = 1.0e-10
    CALL getin('ksta_ter',ksta_ter_omp)

    !Config Key  = f_ri_cd_min
    !Config Desc =
    !Config Def  = 0.1
    !Config Help =
    !
    f_ri_cd_min_omp = 0.1
    CALL getin('f_ri_cd_min',f_ri_cd_min_omp)

    !
    !Config Key  = ok_kzmin
    !Config Desc =
    !Config Def  = .TRUE.
    !Config Help =
    !
    ok_kzmin_omp = .TRUE.
    CALL getin('ok_kzmin',ok_kzmin_omp)

    pbl_lmixmin_alpha_omp=0.0
    CALL getin('pbl_lmixmin_alpha',pbl_lmixmin_alpha_omp)

    !
    !Config Key  = fmagic
    !Config Desc = additionnal multiplicator factor used for albedo
    !Config Def  = 1.
    !Config Help = additionnal multiplicator factor used in albedo.F
    !
    fmagic_omp = 1.
    CALL getin('fmagic',fmagic_omp)

    !
    !Config Key  = pmagic
    !Config Desc = additional factor used for albedo
    !Config Def  = 0.
    !Config Help = additional factor used in albedo.F
    !
    pmagic_omp = 0.
    CALL getin('pmagic',pmagic_omp)


    !Config Key = ok_lic_melt
    !Config Desc = Prise en compte de la fonte de la calotte dans le bilan d'eau
    !Config Def  = .FALSE.
    !Config Help = mettre a .FALSE. pour assurer la conservation en eau
    ok_lic_melt_omp = .FALSE.
    CALL getin('ok_lic_melt', ok_lic_melt_omp)


    !Config Key = ok_lic_cond
    !Config Desc = Prise en compte depot de vapeur d'eau sur la calotte dans le bilan d'eau
    !Config Def  = .FALSE.
    !Config Help = mettre a .TRUE. pour assurer la conservation en eau
    ok_lic_cond_omp = .FALSE.
    CALL getin('ok_lic_cond', ok_lic_cond_omp)

    !
    ! PARAMETER FOR THE PLANETARY BOUNDARY LAYER
    !

    !Config Key  = iflag_pbl
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    iflag_pbl_omp = 1
    CALL getin('iflag_pbl',iflag_pbl_omp)

!FC
    !Config Key  = ifl_pbltree
    !Config Desc = drag from trees 0 no activated
    !Config Def  = 0
    !Config Help =
    !
    ifl_pbltree_omp = 0
    CALL getin('ifl_pbltree',ifl_pbltree_omp)
!FC
    !Config Key  = Cd_frein
    !Config Desc = drag from trees 
    !Config Def  = 7.5E-02 (valeur Masson mais fait planter avec des LAI eleves)
    !Config Help =
    !
    Cd_frein_omp = 7.5E-02
    CALL getin('Cd_frein',Cd_frein_omp)

    !
    !Config Key  = iflag_pbl_split
    !Config Desc = decimal flag: least signif digit = split vdf; next digit = split thermals
    !Config Def  = 0
    !Config Help = 0-> no splitting; 1-> vdf splitting; 10-> thermals splitting; 11-> full splitting
    !
    iflag_pbl_split_omp = 0
    call getin('iflag_pbl_split',iflag_pbl_split_omp)
    !
    !Config Key  = iflag_order2_sollw
    !Config Desc =
    !Config Def  = 0
    !Config Help =
    !
    iflag_order2_sollw_omp = 0
    CALL getin('iflag_order2_sollw',iflag_order2_sollw_omp)
    !
    !Config Key  = iflag_thermals
    !Config Desc =
    !Config Def  = 0
    !Config Help =
    !
    iflag_thermals_omp = 0
    CALL getin('iflag_thermals',iflag_thermals_omp)
    !
    !Config Key  = iflag_thermals_ed
    !Config Desc =
    !Config Def  = 0
    !Config Help =
    !
    fact_thermals_ed_dz_omp = 0.1

    CALL getin('fact_thermals_ed_dz',fact_thermals_ed_dz_omp)
    !
    !
    !Config Key  = iflag_thermals_ed
    !Config Desc =
    !Config Def  = 0
    !Config Help =
    !
    iflag_thermals_ed_omp = 0
    CALL getin('iflag_thermals_ed',iflag_thermals_ed_omp)
    !
    !
    !Config Key  = iflag_thermals_optflux
    !Config Desc =
    !Config Def  = 0
    !Config Help =
    !
    iflag_thermals_optflux_omp = 0
    CALL getin('iflag_thermals_optflux',iflag_thermals_optflux_omp)
    !
    !Config Key  = iflag_thermals_closure
    !Config Desc =
    !Config Def  = 0
    !Config Help =
    !
    iflag_thermals_closure_omp = 1
    CALL getin('iflag_thermals_closure',iflag_thermals_closure_omp)
    !
    !Config Key  = nsplit_thermals
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    nsplit_thermals_omp = 1
    CALL getin('nsplit_thermals',nsplit_thermals_omp)

    !Config Key  = alp_bl_k
    !Config Desc =
    !Config Def  = 0.
    !Config Help =
    !
    alp_bl_k_omp = 1.
    CALL getin('alp_bl_k',alp_bl_k_omp)

    ! nrlmd le 10/04/2012

    !Config Key  = iflag_trig_bl
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_trig_bl_omp = 0
    CALL getin('iflag_trig_bl',iflag_trig_bl_omp)

    !Config Key  = s_trig_bl
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    s_trig_omp = 2e7
    CALL getin('s_trig',s_trig_omp)

    !Config Key  = tau_trig_shallow
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    tau_trig_shallow_omp = 600
    CALL getin('tau_trig_shallow',tau_trig_shallow_omp)

    !Config Key  = tau_trig_deep
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    tau_trig_deep_omp = 1800
    CALL getin('tau_trig_deep',tau_trig_deep_omp)

    !Config Key  = iflag_clos_bl
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_clos_bl_omp = 0
    CALL getin('iflag_clos_bl',iflag_clos_bl_omp)

    ! fin nrlmd le 10/04/2012

    !
    !Config Key  = tau_thermals
    !Config Desc =
    !Config Def  = 0.
    !Config Help =
    !
    tau_thermals_omp = 0.
    CALL getin('tau_thermals',tau_thermals_omp)

    !
    !Config Key  = iflag_coupl
    !Config Desc =
    !Config Def  = 0
    !Config Help =
    !
    iflag_coupl_omp = 0
    CALL getin('iflag_coupl',iflag_coupl_omp)

    !
    !Config Key  = iflag_clos
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_clos_omp = 1
    CALL getin('iflag_clos',iflag_clos_omp)
    !
    !Config Key  = coef_clos_ls
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    coef_clos_ls_omp = 0.
    CALL getin('coef_clos_ls',coef_clos_ls_omp)

    !
    !Config Key  = iflag_cvl_sigd
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_cvl_sigd_omp = 0
    CALL getin('iflag_cvl_sigd',iflag_cvl_sigd_omp)

    !Config Key  = iflag_wake
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_wake_omp = 0
    CALL getin('iflag_wake',iflag_wake_omp)

    !Config Key  = alp_offset
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    alp_offset_omp = 0.
    CALL getin('alp_offset',alp_offset_omp)

    !
    !Config Key  = lev_histhf
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    lev_histhf_omp = 1
    CALL getin('lev_histhf',lev_histhf_omp)

    !
    !Config Key  = lev_histday
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    lev_histday_omp = 1
    CALL getin('lev_histday',lev_histday_omp)

    !
    !Config Key  = lev_histmth
    !Config Desc =
    !Config Def  = 2
    !Config Help =
    !
    lev_histmth_omp = 2
    CALL getin('lev_histmth',lev_histmth_omp)
    !
    !Config Key  = lev_histins
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    lev_histins_omp = 1
    CALL getin('lev_histins',lev_histins_omp)
    !
    !Config Key  = lev_histLES
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    lev_histLES_omp = 1
    CALL getin('lev_histLES',lev_histLES_omp)
    ! 
    !Config Key  = lev_histdayNMC
    !Config Desc =
    !Config Def  = 8
    !Config Help =
    !
    lev_histdayNMC_omp = 8
    CALL getin('lev_histdayNMC',lev_histdayNMC_omp)
    !
    !Config Key  = levout_histNMC
    !Config Desc =
    !Config Def  = 5
    !Config Help =
    !
    levout_histNMC_omp(1) = 5
    levout_histNMC_omp(2) = 5
    levout_histNMC_omp(3) = 5
    CALL getin('levout_histNMC',levout_histNMC_omp)
    !
    !histNMC BEG
    !Config Key  = ok_histNMC
    !Config Desc = ok_histNMC(1) = frequence de sortie fichiers histmthNMC
    !Config Desc = ok_histNMC(2) = frequence de sortie fichiers histdayNMC
    !Config Desc = ok_histNMC(3) = frequence de sortie fichiers histhfNMC
    !Config Def  = n, n, n
    !Config Help =
    !
    ok_histNMC_omp(1) = .FALSE.
    ok_histNMC_omp(2) = .FALSE.
    ok_histNMC_omp(3) = .FALSE.
    CALL getin('ok_histNMC',ok_histNMC_omp)
    !
    !Config Key  = freq_outNMC
    !Config Desc = freq_outNMC(1) = frequence de sortie fichiers histmthNMC
    !Config Desc = freq_outNMC(2) = frequence de sortie fichiers histdayNMC
    !Config Desc = freq_outNMC(3) = frequence de sortie fichiers histhfNMC
    !Config Def  = 2592000., 86400., 21600. (1mois, 1jour, 6h) 
    !Config Help =
    !
    freq_outNMC_omp(1) = mth_len
    freq_outNMC_omp(2) = 1.
    freq_outNMC_omp(3) = 1./4.
    CALL getin('freq_outNMC',freq_outNMC_omp)
    !
    !Config Key  = freq_calNMC
    !Config Desc = freq_calNMC(1) = frequence de calcul fichiers histmthNMC
    !Config Desc = freq_calNMC(2) = frequence de calcul fichiers histdayNMC
    !Config Desc = freq_calNMC(3) = frequence de calcul fichiers histhfNMC
    !Config Def  = phys_tstep
    !Config Help =
    !
    freq_calNMC_omp(1) = phys_tstep
    freq_calNMC_omp(2) = phys_tstep
    freq_calNMC_omp(3) = phys_tstep
    CALL getin('freq_calNMC',freq_calNMC_omp)
    !
    !Config Key  = type_run
    !Config Desc =
    !Config Def  = 'AMIP'/'CFMIP'  ou 'CLIM'/'ENSP'
    !Config Help =
    !
    type_run_omp = 'AMIP'
    CALL getin('type_run',type_run_omp)

    !
    !Config Key  = ok_cosp
    !Config Desc =
    !Config Def  = .FALSE.
    !Config Help =
    !
    ok_cosp_omp = .FALSE.
    CALL getin('ok_cosp',ok_cosp_omp)

    !
    !Config Key  = ok_airs
    !Config Desc =
    !Config Def  = .FALSE.
    !Config Help =
    !
    ok_airs_omp = .FALSE.
    CALL getin('ok_airs',ok_airs_omp)

    !
    !Config Key  = ok_mensuelCOSP
    !Config Desc =
    !Config Def  = .TRUE.
    !Config Help =
    !
    ok_mensuelCOSP_omp = .TRUE.
    CALL getin('ok_mensuelCOSP',ok_mensuelCOSP_omp)

    !
    !Config Key  = ok_journeCOSP
    !Config Desc =
    !Config Def  = .TRUE.
    !Config Help = 
    !
    ok_journeCOSP_omp = .TRUE.
    CALL getin('ok_journeCOSP',ok_journeCOSP_omp)

    !
    !Config Key  = ok_hfCOSP
    !Config Desc =
    !Config Def  = .FALSE.
    !Config Help =
    !
    ok_hfCOSP_omp = .FALSE.
    CALL getin('ok_hfCOSP',ok_hfCOSP_omp)

    !
    ! coordonnees (lonmin_ins, lonmax_ins, latmin_ins, latmax_ins) pour la zone 
    ! avec sorties instantannees tous les pas de temps de la physique => "histbilKP_ins.nc"
    !
    !Config Key  = lonmin_ins
    !Config Desc = 100.  
    !Config Def  = longitude minimale sorties "bilKP_ins"
    !Config Help = 
    !
    lonmin_ins_omp = 100.
    CALL getin('lonmin_ins',lonmin_ins_omp)
    !
    !Config Key  = lonmax_ins
    !Config Desc = 130. 
    !Config Def  = longitude maximale sorties "bilKP_ins"
    !Config Help =
    !
    lonmax_ins_omp = 130.
    CALL getin('lonmax_ins',lonmax_ins_omp)
    !
    !Config Key  = latmin_ins
    !Config Desc = -20.  
    !Config Def  = latitude minimale sorties "bilKP_ins"
    !Config Help = 
    !
    latmin_ins_omp = -20.
    CALL getin('latmin_ins',latmin_ins_omp)
    !
    !Config Key  = latmax_ins
    !Config Desc = 20. 
    !Config Def  = latitude maximale sorties "bilKP_ins"
    !Config Help =
    !
    latmax_ins_omp = 20.
    CALL getin('latmax_ins',latmax_ins_omp)
    !
    !Config Key  = ecrit_hf
    !Config Desc =
    !Config Def  = 1./8. !toutes les 3h
    !Config Help =
    !
    ecrit_hf_omp = 1./8.
    CALL getin('ecrit_hf',ecrit_hf_omp)
    !
    !Config Key  = ecrit_ins
    !Config Desc =
    !Config Def  = 1./48. ! toutes les 1/2 h
    !Config Help =
    !
    ecrit_ins_omp = 1./48.
    CALL getin('ecrit_ins',ecrit_ins_omp)
    !
    !Config Key  = ecrit_day
    !Config Desc =
    !Config Def  = 1.0 !tous les jours
    !Config Help = nombre de jours pour ecriture fichier histday.nc
    !
    ecrit_day_omp = 1.0
    CALL getin('ecrit_day',ecrit_day_omp)
    !
    !Config Key  = ecrit_mth
    !Config Desc =
    !Config Def  = 30. !tous les 30jours (1 fois par mois)
    !Config Help =
    !
    ecrit_mth_omp = 30.
    CALL getin('ecrit_mth',ecrit_mth_omp)
    !
    !Config Key  = ecrit_tra
    !Config Desc =
    !Config Def  = 30. !tous les 30jours (1 fois par mois)
    !Config Help =
    !
    ecrit_tra_omp = 0.
    CALL getin('ecrit_tra',ecrit_tra_omp)
    !
    !Config Key  = ecrit_reg
    !Config Desc =
    !Config Def  = 0.25  !4 fois par jour
    !Config Help =
    !
    ecrit_reg_omp = 0.25   !4 fois par jour
    CALL getin('ecrit_reg',ecrit_reg_omp)
    !
    !
    print*,'CONFPHYS OOK avant drag_ter'
    !
    ! PARAMETRES CDRAG
    !
    f_cdrag_ter_omp = 0.8
    CALL getin('f_cdrag_ter',f_cdrag_ter_omp)
    !
    f_cdrag_oce_omp = 0.8
    CALL getin('f_cdrag_oce',f_cdrag_oce_omp)
    !

    ! Gustiness flags
    f_z0qh_oce_omp = 1.
    CALL getin('f_z0qh_oce',f_z0qh_oce_omp)
    !
    f_qsat_oce_omp = 1.
    CALL getin('f_qsat_oce',f_qsat_oce_omp)
    !
    f_gust_bl_omp = 0.
    CALL getin('f_gust_bl',f_gust_bl_omp)
    !
    f_gust_wk_omp = 0.
    CALL getin('f_gust_wk',f_gust_wk_omp)
    !
    !Config Key  = iflag_z0_oce
    !Config Desc = 0 (z0h=z0m), 1 (diff. equ. for z0h and z0m), -1 (z0m=z0h=z0min)
    !Config Def  = 0   ! z0h = z0m
    !Config Help =
    !
    iflag_z0_oce_omp=0
    CALL getin('iflag_z0_oce',iflag_z0_oce_omp)
    !
    iflag_gusts_omp=0
    CALL getin('iflag_gusts',iflag_gusts_omp)
    !
    min_wind_speed_omp = 1.
    CALL getin('min_wind_speed',min_wind_speed_omp)

    z0m_seaice_omp = 0.002 ; CALL getin('z0m_seaice',z0m_seaice_omp)
    z0h_seaice_omp = 0.002 ; CALL getin('z0h_seaice',z0h_seaice_omp)

    f_rugoro_omp = 0.
    CALL getin('f_rugoro',f_rugoro_omp)

    z0min_omp = 0.000015
    CALL getin('z0min',z0min_omp)


    ! PARAMETERS FOR CONVECTIVE INHIBITION BY TROPOS. DRYNESS
    !
    !Config Key  = supcrit1
    !Config Desc =
    !Config Def  = .540
    !Config Help =
    !
    supcrit1_omp = .540
    CALL getin('supcrit1',supcrit1_omp)

    !
    !Config Key  = supcrit2
    !Config Desc =
    !Config Def  = .600
    !Config Help =
    !
    supcrit2_omp = .600
    CALL getin('supcrit2',supcrit2_omp)

    !
    ! PARAMETERS FOR THE MIXING DISTRIBUTION
    ! iflag_mix: 0=OLD, 
    !            1=NEW (JYG),            
    !            2=NEW + conv. depth inhib. by tropos. dryness
    ! '2' is NOT operationnal and should not be used.
    !
    !Config Key  = iflag_mix
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    iflag_mix_omp = 1
    CALL getin('iflag_mix',iflag_mix_omp)

!
    ! PARAMETERS FOR THE EROSION OF THE ADIABATIC ASCENTS
    ! iflag_mix_adiab: 0=OLD, 
    !                  1=NEW (CR),            
    !            
    !
    !Config Key  = iflag_mix_adiab
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    iflag_mix_adiab_omp = 0
    CALL getin('iflag_mix_adiab',iflag_mix_adiab_omp)

    !
    !Config Key  = scut
    !Config Desc =
    !Config Def  = 0.95
    !Config Help =
    !
    scut_omp = 0.95
    CALL getin('scut',scut_omp)

    !
    !Config Key  = qqa1
    !Config Desc =
    !Config Def  = 1.0
    !Config Help =
    !
    qqa1_omp = 1.0
    CALL getin('qqa1',qqa1_omp)

    !
    !Config Key  = qqa2
    !Config Desc =
    !Config Def  = 0.0
    !Config Help =
    !
    qqa2_omp = 0.0
    CALL getin('qqa2',qqa2_omp)

    !
    !Config Key  = gammas
    !Config Desc =
    !Config Def  = 0.05
    !Config Help =
    !
    gammas_omp = 0.05
    CALL getin('gammas',gammas_omp)

    !
    !Config Key  = Fmax
    !Config Desc =
    !Config Def  = 0.65
    !Config Help =
    !
    Fmax_omp = 0.65
    CALL getin('Fmax',Fmax_omp)

    !
    !Config Key  = tmax_fonte_cv
    !Config Desc =
    !Config Def  = 275.15
    !Config Help =
    !
    tmax_fonte_cv_omp = 275.15
    CALL getin('tmax_fonte_cv',tmax_fonte_cv_omp)

    !
    !Config Key  = alphas  
    !Config Desc =
    !Config Def  = -5.
    !Config Help =
    !
    alphas_omp = -5.
    CALL getin('alphas',alphas_omp)

    !Config key = ok_strato
    !Config  Desc = activation de la version strato
    !Config  Def  = .FALSE.
    !Config  Help = active la version stratosph\'erique de LMDZ de F. Lott
    !               Et la sponge layer (Runs Stratospheriques)

    ok_strato_omp=.FALSE.
    CALL getin('ok_strato',ok_strato_omp)

    !Config  key = ok_hines
    !Config  Desc = activation de la parametrisation de hines
    !Config  Def  = .FALSE.
    !Config  Help = Clefs controlant la parametrization de Hines

    ok_hines_omp=.FALSE.
    CALL getin('ok_hines',ok_hines_omp)

    !  Parametres pour les ondes de gravite
    !  
    !  Subgrid Scale Orography (Lott Miller (1997), Lott (1999))

    sso_gkdrag_omp = merge(0.1875, 0.2, ok_strato_omp)
    CALL getin('sso_gkdrag', sso_gkdrag_omp)

    sso_grahil_omp=merge(0.1,1.,ok_strato_omp)
    CALL getin('sso_grahil', sso_grahil_omp)

    sso_grcrit_omp =merge(1.,0.01,ok_strato_omp)
    CALL getin('sso_grcrit', sso_grcrit_omp)

    sso_gfrcri_omp = 1.
    CALL getin('sso_gfrcri', sso_gfrcri_omp)

    sso_gkwake_omp = 0.50
    CALL getin('sso_gkwake', sso_gkwake_omp)

    sso_gklift_omp = merge(0.25,0.50,ok_strato_omp)
    CALL getin('sso_gklift', sso_gklift_omp)

    ! Random gravity waves:

    ok_gwd_rando_omp = .FALSE.
    IF ( klon_glo == 1 ) THEN
       print*,'La parametrisation des ondes de gravites non orographiques'
       print*,'ne fonctionne pas en 1D'
    ELSE
       CALL getin('ok_gwd_rando', ok_gwd_rando_omp)
    ENDIF

    gwd_rando_ruwmax_omp = 2.00
    CALL getin('gwd_rando_ruwmax', gwd_rando_ruwmax_omp)

    gwd_rando_sat_omp = 0.25
    CALL getin('gwd_rando_sat', gwd_rando_sat_omp)

    gwd_front_ruwmax_omp = 2.50
    CALL getin('gwd_front_ruwmax', gwd_front_ruwmax_omp)

    gwd_front_sat_omp = 0.60
    CALL getin('gwd_front_sat', gwd_front_sat_omp)


    !Config  key = ok_qch4
    !Config  Desc = activation de la parametrisation du methane
    !Config  Def  = .FALSE.
    !Config  Help = Clef controlant l'activation de la parametrisation
    !               de l'humidite due a oxydation+photolyse du methane strato

    ok_qch4_omp=.FALSE.
    CALL getin('ok_qch4',ok_qch4_omp)

    !Config Key  = OK_LES                                               
    !Config Desc = Pour des sorties LES                                 
    !Config Def  = .FALSE.                                              
    !Config Help = Pour creer le fichier histLES contenant les sorties  
    !              LES                                                  
    !                                                                   
    ok_LES_omp = .FALSE.                                              
    CALL getin('OK_LES', ok_LES_omp)                                  

    !Config Key  = callstats                                               
    !Config Desc = Pour des sorties callstats                                 
    !Config Def  = .FALSE.                                              
    !Config Help = Pour creer le fichier stats contenant les sorties  
    !              stats                                                  
    !                                                                   
    callstats_omp = .FALSE.                                              
    CALL getin('callstats', callstats_omp)                                  
    !
    !Config Key  = ecrit_LES
    !Config Desc = Frequence d'ecriture des resultats du LES en nombre de jours;
    !              par defaut 1., i.e. 1 jour
    !Config Def  = 1./8.
    !Config Help = ... 
    !
    !
    adjust_tropopause = .FALSE.
    CALL getin('adjust_tropopause', adjust_tropopause_omp)
    !
    !Config Key  = adjust_tropopause
    !Config Desc = Adjust the ozone field from the climoz file by stretching its
    !              tropopause so that it matches the one of LMDZ.
    !Config Def  = .FALSE.
    !Config Help = Ensure tropospheric ozone column conservation.
    !
    !
    ok_daily_climoz = .FALSE.
    CALL getin('ok_daily_climoz', ok_daily_climoz_omp)
    !
    !Config Key  = ok_daily_climoz
    !Config Desc = Interpolate in time the ozone forcings within ce0l.
    !              .TRUE. if backward compatibility is needed.
    !Config Def  = .TRUE.
    !Config Help = .FALSE. ensure much fewer (no calendar dependency)
    !  and lighter monthly climoz files, inetrpolated in time at gcm run time.
    !
    ecrit_LES_omp = 1./8.
    CALL getin('ecrit_LES', ecrit_LES_omp)
    !
    read_climoz = 0 ! default value
    CALL getin('read_climoz', read_climoz)

    carbon_cycle_tr_omp=.FALSE.
    CALL getin('carbon_cycle_tr',carbon_cycle_tr_omp)

    carbon_cycle_cpl_omp=.FALSE.
    CALL getin('carbon_cycle_cpl',carbon_cycle_cpl_omp)

    carbon_cycle_rad_omp=.FALSE.
    CALL getin('carbon_cycle_rad',carbon_cycle_rad_omp)

    ! >> PC
    ! level_coupling_esm : level of coupling of the biogeochemical fields between LMDZ, ORCHIDEE and NEMO
    ! Definitions of level_coupling_esm in physiq.def
    ! level_coupling_esm = 0  ! No field exchange between LMDZ and ORCHIDEE models
    !                         ! No field exchange between LMDZ and NEMO
    ! level_coupling_esm = 1  ! Field exchange between LMDZ and ORCHIDEE models
    !                         ! No field exchange between LMDZ and NEMO models
    ! level_coupling_esm = 2  ! No field exchange between LMDZ and ORCHIDEE models
    !                         ! Field exchange between LMDZ and NEMO models
    ! level_coupling_esm = 3  ! Field exchange between LMDZ and ORCHIDEE models
    !                         ! Field exchange between LMDZ and NEMO models
    level_coupling_esm_omp=0 ! default value
    CALL getin('level_coupling_esm',level_coupling_esm_omp)
    ! << PC

    !$OMP END MASTER
    !$OMP BARRIER

    R_ecc = R_ecc_omp
    R_peri = R_peri_omp
    R_incl = R_incl_omp
    solaire = solaire_omp
    ok_suntime_rrtm = ok_suntime_rrtm_omp
    co2_ppm = co2_ppm_omp
    RCO2 = RCO2_omp
    CH4_ppb = CH4_ppb_omp
    RCH4 = RCH4_omp
    N2O_ppb = N2O_ppb_omp
    RN2O = RN2O_omp
    CFC11_ppt = CFC11_ppt_omp
    RCFC11 = RCFC11_omp
    CFC12_ppt = CFC12_ppt_omp
    RCFC12 = RCFC12_omp
    RCO2_act = RCO2
    RCH4_act = RCH4
    RN2O_act = RN2O
    RCFC11_act = RCFC11
    RCFC12_act = RCFC12
    RCO2_per = RCO2_per_omp
    RCH4_per = RCH4_per_omp
    RN2O_per = RN2O_per_omp
    RCFC11_per = RCFC11_per_omp
    RCFC12_per = RCFC12_per_omp

    iflag_cycle_diurne = iflag_cycle_diurne_omp
    soil_model = soil_model_omp
    new_oliq = new_oliq_omp
    ok_orodr = ok_orodr_omp
    ok_orolf = ok_orolf_omp
    ok_limitvrai = ok_limitvrai_omp
    nbapp_rad = nbapp_rad_omp
    iflag_con = iflag_con_omp
    nbapp_cv = nbapp_cv_omp
    nbapp_wk = nbapp_wk_omp
    iflag_ener_conserv = iflag_ener_conserv_omp
    ok_conserv_q = ok_conserv_q_omp
    iflag_fisrtilp_qsat = iflag_fisrtilp_qsat_omp
    iflag_bergeron = iflag_bergeron_omp

    epmax = epmax_omp
    coef_epmax_cape = coef_epmax_cape_omp
    ok_adj_ema = ok_adj_ema_omp
    iflag_clw = iflag_clw_omp
    cld_lc_lsc = cld_lc_lsc_omp
    cld_lc_con = cld_lc_con_omp
    cld_tau_lsc = cld_tau_lsc_omp
    cld_tau_con = cld_tau_con_omp
    ffallv_lsc = ffallv_lsc_omp
    ffallv_con = ffallv_con_omp
    coef_eva = coef_eva_omp
    reevap_ice = reevap_ice_omp
    iflag_pdf = iflag_pdf_omp
    solarlong0 = solarlong0_omp
    qsol0 = qsol0_omp
    evap0 = evap0_omp
    albsno0 = albsno0_omp
    iflag_sic = iflag_sic_omp
    inertie_sol = inertie_sol_omp
    inertie_sic = inertie_sic_omp
    inertie_lic = inertie_lic_omp
    inertie_sno = inertie_sno_omp
    rad_froid = rad_froid_omp
    rad_chau1 = rad_chau1_omp
    rad_chau2 = rad_chau2_omp
    t_glace_min = t_glace_min_omp
    t_glace_max = t_glace_max_omp
    exposant_glace = exposant_glace_omp
    iflag_t_glace = iflag_t_glace_omp
    iflag_cloudth_vert=iflag_cloudth_vert_omp
    iflag_rain_incloud_vol=iflag_rain_incloud_vol_omp
    iflag_ice_thermo = iflag_ice_thermo_omp
    rei_min = rei_min_omp
    rei_max = rei_max_omp
    top_height = top_height_omp
    overlap = overlap_omp
    cdmmax = cdmmax_omp
    cdhmax = cdhmax_omp
    ksta = ksta_omp
    ksta_ter = ksta_ter_omp
    f_ri_cd_min = f_ri_cd_min_omp
    ok_kzmin = ok_kzmin_omp
    pbl_lmixmin_alpha=pbl_lmixmin_alpha_omp
    fmagic = fmagic_omp
    pmagic = pmagic_omp
    iflag_pbl = iflag_pbl_omp
    iflag_pbl_split = iflag_pbl_split_omp
!FC
    ifl_pbltree = ifl_pbltree_omp
    Cd_frein    =Cd_frein_omp
    iflag_order2_sollw = iflag_order2_sollw_omp
    lev_histhf = lev_histhf_omp
    lev_histday = lev_histday_omp
    lev_histmth = lev_histmth_omp
    lev_histins = lev_histins_omp
    lev_histLES = lev_histLES_omp
    lev_histdayNMC = lev_histdayNMC_omp
    levout_histNMC = levout_histNMC_omp
    ok_histNMC(:) = ok_histNMC_omp(:)
    freq_outNMC(:) = freq_outNMC_omp(:)
    freq_calNMC(:) = freq_calNMC_omp(:)

    type_ocean = type_ocean_omp
    version_ocean = version_ocean_omp
    t_coupl = t_coupl_omp

    ok_veget=.TRUE.
    type_veget=type_veget_omp
    IF (type_veget=='n' .or. type_veget=='bucket' .or. type_veget=='betaclim') THEN
       ok_veget=.FALSE.
    ENDIF
    ! Martin
    ok_snow = ok_snow_omp
    ! Martin

    ok_all_xml = ok_all_xml_omp
    ok_lwoff = ok_lwoff_omp
    ok_newmicro = ok_newmicro_omp
    ok_journe = ok_journe_omp
    ok_hf = ok_hf_omp
    ok_mensuel = ok_mensuel_omp
    ok_instan = ok_instan_omp
    freq_ISCCP = freq_ISCCP_omp
    ecrit_ISCCP = ecrit_ISCCP_omp
    freq_COSP = freq_COSP_omp
    freq_AIRS = freq_AIRS_omp
    ok_ade = ok_ade_omp
    ok_aie = ok_aie_omp
    ok_alw = ok_alw_omp
    ok_cdnc = ok_cdnc_omp
    ok_volcan = ok_volcan_omp
    aerosol_couple = aerosol_couple_omp
    chemistry_couple = chemistry_couple_omp
    flag_aerosol=flag_aerosol_omp
    flag_aerosol_strat=flag_aerosol_strat_omp
    flag_aer_feedback=flag_aer_feedback_omp
    flag_bc_internal_mixture=flag_bc_internal_mixture_omp
    aer_type = aer_type_omp
    bl95_b0 = bl95_b0_omp
    bl95_b1 = bl95_b1_omp
    fact_cldcon = fact_cldcon_omp
    facttemps = facttemps_omp
    ratqsbas = ratqsbas_omp
    ratqshaut = ratqshaut_omp
    tau_ratqs = tau_ratqs_omp

    iflag_radia = iflag_radia_omp
    iflag_rrtm = iflag_rrtm_omp
    iflag_albedo = iflag_albedo_omp
    ok_chlorophyll = ok_chlorophyll_omp
    NSW = NSW_omp
    iflag_cld_th = iflag_cld_th_omp
    iflag_cld_cv = iflag_cld_cv_omp
    tau_cld_cv = tau_cld_cv_omp
    coefw_cld_cv = coefw_cld_cv_omp
    iflag_ratqs = iflag_ratqs_omp
    ip_ebil_phy = ip_ebil_phy_omp
    iflag_thermals = iflag_thermals_omp
    iflag_thermals_ed = iflag_thermals_ed_omp
    fact_thermals_ed_dz = fact_thermals_ed_dz_omp
    iflag_thermals_optflux = iflag_thermals_optflux_omp
    iflag_thermals_closure = iflag_thermals_closure_omp
    nsplit_thermals = nsplit_thermals_omp
    tau_thermals = tau_thermals_omp
    alp_bl_k = alp_bl_k_omp
    ! nrlmd le 10/04/2012
    iflag_trig_bl = iflag_trig_bl_omp
    s_trig = s_trig_omp
    tau_trig_shallow = tau_trig_shallow_omp
    tau_trig_deep = tau_trig_deep_omp
    iflag_clos_bl = iflag_clos_bl_omp
    ! fin nrlmd le 10/04/2012
    iflag_coupl = iflag_coupl_omp
    iflag_clos = iflag_clos_omp
    iflag_wake = iflag_wake_omp
    coef_clos_ls = coef_clos_ls_omp
    alp_offset = alp_offset_omp
    iflag_cvl_sigd = iflag_cvl_sigd_omp
    type_run = type_run_omp
    ok_cosp = ok_cosp_omp
    ok_airs = ok_airs_omp

    ok_mensuelCOSP = ok_mensuelCOSP_omp
    ok_journeCOSP = ok_journeCOSP_omp
    ok_hfCOSP = ok_hfCOSP_omp
    seuil_inversion=seuil_inversion_omp
    lonmin_ins = lonmin_ins_omp
    lonmax_ins = lonmax_ins_omp
    latmin_ins = latmin_ins_omp
    latmax_ins = latmax_ins_omp
    ecrit_hf   = ecrit_hf_omp
    ecrit_ins   = ecrit_ins_omp
    ecrit_day = ecrit_day_omp
    ecrit_mth = ecrit_mth_omp
    ecrit_tra = ecrit_tra_omp
    ecrit_reg = ecrit_reg_omp
    cvl_comp_threshold = cvl_comp_threshold_omp
    cvl_sig2feed = cvl_sig2feed_omp
    cvl_corr = cvl_corr_omp
    ok_lic_melt = ok_lic_melt_omp
    ok_lic_cond = ok_lic_cond_omp
    f_cdrag_ter=f_cdrag_ter_omp
    f_cdrag_oce=f_cdrag_oce_omp

    f_gust_wk=f_gust_wk_omp
    f_gust_bl=f_gust_bl_omp
    f_qsat_oce=f_qsat_oce_omp
    f_z0qh_oce=f_z0qh_oce_omp
    min_wind_speed=min_wind_speed_omp
    iflag_gusts=iflag_gusts_omp
    iflag_z0_oce=iflag_z0_oce_omp

    z0m_seaice=z0m_seaice_omp
    z0h_seaice=z0h_seaice_omp

    f_rugoro=f_rugoro_omp

    z0min=z0min_omp
    supcrit1 = supcrit1_omp
    supcrit2 = supcrit2_omp
    iflag_mix = iflag_mix_omp
    iflag_mix_adiab = iflag_mix_adiab_omp
    scut = scut_omp
    qqa1 = qqa1_omp
    qqa2 = qqa2_omp
    gammas = gammas_omp
    Fmax = Fmax_omp
    tmax_fonte_cv = tmax_fonte_cv_omp
    alphas = alphas_omp

    gkdrag=sso_gkdrag_omp
    grahilo=sso_grahil_omp
    grcrit=sso_grcrit_omp
    gfrcrit=sso_gfrcri_omp
    gkwake=sso_gkwake_omp 
    gklift=sso_gklift_omp 

    ok_strato = ok_strato_omp
    ok_hines = ok_hines_omp
    ok_gwd_rando = ok_gwd_rando_omp
    gwd_rando_ruwmax = gwd_rando_ruwmax_omp
    gwd_rando_sat = gwd_rando_sat_omp
    gwd_front_ruwmax = gwd_front_ruwmax_omp
    gwd_front_sat = gwd_front_sat_omp
    ok_qch4 = ok_qch4_omp
    ok_LES = ok_LES_omp
    callstats = callstats_omp
    ecrit_LES = ecrit_LES_omp
    adjust_tropopause = adjust_tropopause_omp
    ok_daily_climoz = ok_daily_climoz_omp
    carbon_cycle_tr = carbon_cycle_tr_omp
    carbon_cycle_cpl = carbon_cycle_cpl_omp
    carbon_cycle_rad = carbon_cycle_rad_omp
    level_coupling_esm = level_coupling_esm_omp

    ! Test of coherence between type_ocean and version_ocean
    IF (type_ocean=='couple' .AND. (version_ocean/='opa8' .AND. version_ocean/='nemo') ) THEN
       WRITE(lunout,*)' ERROR version_ocean=',version_ocean,' not valid in coupled configuration'
       CALL abort_physic('conf_phys','version_ocean not valid',1)
    ENDIF

    IF (type_ocean=='slab' .AND. version_ocean=='xxxxxx') THEN
       version_ocean='sicOBS'
    ELSE IF (type_ocean=='slab' .AND. version_ocean/='sicOBS' &
         .AND. version_ocean/='sicINT' .AND. version_ocean/='sicNO') THEN
       WRITE(lunout,*)' ERROR version_ocean=',version_ocean,' not valid with slab ocean'
       CALL abort_physic('conf_phys','version_ocean not valid',1)
    ENDIF

    !--test on radiative scheme 
    IF (iflag_rrtm .EQ. 0) THEN 
      IF (NSW.NE.2) THEN 
        WRITE(lunout,*) ' ERROR iflag_rrtm=0 and NSW<>2 not possible'
        CALL abort_physic('conf_phys','choice NSW not valid',1)
      ENDIF
    ELSE IF (iflag_rrtm .EQ. 1) THEN
      IF (NSW.NE.2.AND.NSW.NE.4.AND.NSW.NE.6) THEN
        WRITE(lunout,*) ' ERROR iflag_rrtm=1 and NSW<>2,4,6 not possible'
        CALL abort_physic('conf_phys','choice NSW not valid',1)
      ENDIF
    ELSE 
       WRITE(lunout,*) ' ERROR iflag_rrtm<>0,1'
       CALL abort_physic('conf_phys','choice iflag_rrtm not valid',1)
    ENDIF
    !--here we test that solaire has not been changed if ok_suntime_rrtm is activated
!    IF (ok_suntime_rrtm.AND.ABS(solaire-solaire_omp_init).GT.1.E-7) THEN 
!       WRITE(lunout,*) ' ERROR ok_suntime_rrtm=y and solaire is provided in def file'
!       CALL abort_physic('conf_phys','ok_suntime_rrtm=y and solaire is provided',1)
!    ENDIF
#ifdef CPP_StratAer
    IF (iflag_rrtm .NE. 1) THEN 
       WRITE(lunout,*) ' ERROR iflag_rrtm<>1 but StratAer activated'
       CALL abort_physic('conf_phys','iflag_rrtm not valid for StratAer',1)
    ENDIF
    IF (NSW .NE. 6) THEN 
       WRITE(lunout,*) ' ERROR NSW<>6 but StratAer activated'
       CALL abort_physic('conf_phys','NSW not valid for StratAer',1)
    ENDIF
#endif

    !--test on ocean surface albedo
    IF (iflag_albedo.LT.0.OR.iflag_albedo.GT.2) THEN
       WRITE(lunout,*) ' ERROR iflag_albedo<>0,1'
       CALL abort_physic('conf_phys','choice iflag_albedo not valid',1)
    ENDIF

    ! Test sur new_aod. Ce flag permet de retrouver les resultats de l'AR4
    ! il n'est utilisable que lors du couplage avec le SO4 seul 
    IF (ok_ade .OR. ok_aie) THEN 
       IF ( flag_aerosol .EQ. 0 ) THEN
          CALL abort_physic('conf_phys','flag_aerosol=0 not compatible avec ok_ade ou ok_aie=.TRUE.',1)
       ENDIF
    ENDIF

    ! Flag_aerosol cannot be to zero if we are in coupled mode for aerosol
    IF (aerosol_couple .AND. flag_aerosol .EQ. 0 ) THEN
       CALL abort_physic('conf_phys', 'flag_aerosol cannot be to zero if aerosol_couple=y ', 1)
    ENDIF

    ! Read_climoz need to be zero if we are in couple mode for chemistry 
    IF (chemistry_couple .AND. read_climoz .ne. 0) THEN 
       CALL abort_physic('conf_phys', 'read_climoz need to be to zero if chemistry_couple=y ', 1) 
    ENDIF

    ! flag_aerosol need to be different to zero if ok_cdnc is activated
    IF (ok_cdnc .AND. flag_aerosol .EQ. 0) THEN
       CALL abort_physic('conf_phys', 'flag_aerosol cannot be to zero if ok_cdnc is activated ', 1)
    ENDIF

    ! ok_cdnc must be set to y if ok_aie is activated
    IF (ok_aie .AND. .NOT. ok_cdnc) THEN
       CALL abort_physic('conf_phys', 'ok_cdnc must be set to y if ok_aie is activated',1)
    ENDIF

    ! flag_aerosol=7 => MACv2SP climatology 
    IF (flag_aerosol.EQ.7.AND. iflag_rrtm.NE.1) THEN
       CALL abort_physic('conf_phys', 'flag_aerosol=7 (MACv2SP) can only be activated with RRTM',1)
    ENDIF
    IF (flag_aerosol.EQ.7.AND. NSW.NE.6) THEN
       CALL abort_physic('conf_phys', 'flag_aerosol=7 (MACv2SP) can only be activated with NSW=6',1)
    ENDIF

    ! BC internal mixture is only possible with RRTM & NSW=6 & flag_aerosol=6 or aerosol_couple
    IF (flag_bc_internal_mixture .AND. NSW.NE.6) THEN 
       CALL abort_physic('conf_phys', 'flag_bc_internal_mixture can only be activated with NSW=6',1)
    ENDIF
    IF (flag_bc_internal_mixture .AND. iflag_rrtm.NE.1) THEN 
       CALL abort_physic('conf_phys', 'flag_bc_internal_mixture can only be activated with RRTM',1)
    ENDIF
    IF (flag_bc_internal_mixture .AND. flag_aerosol.NE.6) THEN 
       CALL abort_physic('conf_phys', 'flag_bc_internal_mixture can only be activated with flag_aerosol=6',1)
    ENDIF

    ! Test on carbon cycle
    IF (carbon_cycle_tr .AND. .NOT. carbon_cycle_cpl) THEN
       CALL abort_physic('conf_phys', 'carbon_cycle_cpl has to be TRUE if carbon_cycle_tr is on',1)
    ENDIF
    IF (carbon_cycle_rad .AND. .NOT. carbon_cycle_cpl) THEN
       CALL abort_physic('conf_phys', 'carbon_cycle_cpl has to be TRUE if carbon_cycle_rad is on',1)
    ENDIF

    ! ORCHIDEE must be activated for ifl_pbltree=1
    IF (.NOT. ok_veget .AND. ifl_pbltree==1) THEN
       WRITE(lunout,*)'Warning: ORCHIDEE must be activated for ifl_pbltree=1'
       WRITE(lunout,*)'ifl_pbltree is now changed to zero'
       ifl_pbltree=0
    ENDIF

    !$OMP MASTER

    WRITE(lunout,*) ' ##############################################'
    WRITE(lunout,*) ' Configuration des parametres de la physique: '
    WRITE(lunout,*) ' Type ocean = ', type_ocean
    WRITE(lunout,*) ' Version ocean = ', version_ocean
    WRITE(lunout,*) ' Config veget = ', ok_veget,type_veget
    WRITE(lunout,*) ' Snow model SISVAT : ok_snow = ', ok_snow
    WRITE(lunout,*) ' Config xml pour XIOS : ok_all_xml = ', ok_all_xml
    WRITE(lunout,*) ' Sortie journaliere = ', ok_journe
    WRITE(lunout,*) ' Sortie haute frequence = ', ok_hf
    WRITE(lunout,*) ' Sortie mensuelle = ', ok_mensuel
    WRITE(lunout,*) ' Sortie instantanee = ', ok_instan
    WRITE(lunout,*) ' Frequence appel simulateur ISCCP, freq_ISCCP =', freq_ISCCP
    WRITE(lunout,*) ' Frequence appel simulateur ISCCP, ecrit_ISCCP =', ecrit_ISCCP
    WRITE(lunout,*) ' Frequence appel simulateur COSP, freq_COSP =', freq_COSP
    WRITE(lunout,*) ' Frequence appel simulateur AIRS, freq_AIRS =', freq_AIRS
    WRITE(lunout,*) ' Sortie bilan d''energie, ip_ebil_phy =', ip_ebil_phy
    WRITE(lunout,*) ' Excentricite = ',R_ecc
    WRITE(lunout,*) ' Equinoxe = ',R_peri
    WRITE(lunout,*) ' Inclinaison =',R_incl
    WRITE(lunout,*) ' Constante solaire =',solaire
    WRITE(lunout,*) ' ok_suntime_rrtm =',ok_suntime_rrtm
    WRITE(lunout,*) ' co2_ppm =',co2_ppm
    WRITE(lunout,*) ' RCO2_act = ',RCO2_act
    WRITE(lunout,*) ' CH4_ppb =',CH4_ppb,' RCH4_act = ',RCH4_act
    WRITE(lunout,*) ' N2O_ppb =',N2O_ppb,' RN2O_act=  ',RN2O_act
    WRITE(lunout,*) ' CFC11_ppt=',CFC11_ppt,' RCFC11_act=  ',RCFC11_act
    WRITE(lunout,*) ' CFC12_ppt=',CFC12_ppt,' RCFC12_act=  ',RCFC12_act
    WRITE(lunout,*) ' RCO2_per = ',RCO2_per,' RCH4_per = ', RCH4_per
    WRITE(lunout,*) ' RN2O_per = ',RN2O_per,' RCFC11_per = ', RCFC11_per
    WRITE(lunout,*) ' RCFC12_per = ',RCFC12_per
    WRITE(lunout,*) ' cvl_comp_threshold=', cvl_comp_threshold
    WRITE(lunout,*) ' cvl_sig2feed=', cvl_sig2feed
    WRITE(lunout,*) ' cvl_corr=', cvl_corr
    WRITE(lunout,*) ' ok_lic_melt=', ok_lic_melt
    WRITE(lunout,*) ' ok_lic_cond=', ok_lic_cond
    WRITE(lunout,*) ' iflag_cycle_diurne=',iflag_cycle_diurne
    WRITE(lunout,*) ' soil_model=',soil_model
    WRITE(lunout,*) ' new_oliq=',new_oliq
    WRITE(lunout,*) ' ok_orodr=',ok_orodr
    WRITE(lunout,*) ' ok_orolf=',ok_orolf
    WRITE(lunout,*) ' ok_limitvrai=',ok_limitvrai
    WRITE(lunout,*) ' nbapp_rad=',nbapp_rad
    WRITE(lunout,*) ' iflag_con=',iflag_con
    WRITE(lunout,*) ' nbapp_cv=',nbapp_cv
    WRITE(lunout,*) ' nbapp_wk=',nbapp_wk
    WRITE(lunout,*) ' iflag_ener_conserv=',iflag_ener_conserv
    WRITE(lunout,*) ' ok_conserv_q=',ok_conserv_q
    WRITE(lunout,*) ' iflag_fisrtilp_qsat=',iflag_fisrtilp_qsat
    WRITE(lunout,*) ' iflag_bergeron=',iflag_bergeron
    WRITE(lunout,*) ' epmax = ', epmax
    WRITE(lunout,*) ' coef_epmax_cape = ', coef_epmax_cape
    WRITE(lunout,*) ' ok_adj_ema = ', ok_adj_ema
    WRITE(lunout,*) ' iflag_clw = ', iflag_clw
    WRITE(lunout,*) ' cld_lc_lsc = ', cld_lc_lsc
    WRITE(lunout,*) ' cld_lc_con = ', cld_lc_con
    WRITE(lunout,*) ' cld_tau_lsc = ', cld_tau_lsc
    WRITE(lunout,*) ' cld_tau_con = ', cld_tau_con
    WRITE(lunout,*) ' ffallv_lsc = ', ffallv_lsc
    WRITE(lunout,*) ' ffallv_con = ', ffallv_con
    WRITE(lunout,*) ' coef_eva = ', coef_eva
    WRITE(lunout,*) ' reevap_ice = ', reevap_ice
    WRITE(lunout,*) ' iflag_pdf = ', iflag_pdf
    WRITE(lunout,*) ' iflag_cld_th = ', iflag_cld_th
    WRITE(lunout,*) ' iflag_cld_cv = ', iflag_cld_cv
    WRITE(lunout,*) ' tau_cld_cv = ', tau_cld_cv
    WRITE(lunout,*) ' coefw_cld_cv = ', coefw_cld_cv
    WRITE(lunout,*) ' iflag_radia = ', iflag_radia
    WRITE(lunout,*) ' iflag_rrtm = ', iflag_rrtm
    WRITE(lunout,*) ' NSW = ', NSW
    WRITE(lunout,*) ' iflag_albedo = ', iflag_albedo !albedo SB
    WRITE(lunout,*) ' ok_chlorophyll =',ok_chlorophyll ! albedo SB
    WRITE(lunout,*) ' iflag_ratqs = ', iflag_ratqs
    WRITE(lunout,*) ' seuil_inversion = ', seuil_inversion
    WRITE(lunout,*) ' fact_cldcon = ', fact_cldcon
    WRITE(lunout,*) ' facttemps = ', facttemps
    WRITE(lunout,*) ' ok_newmicro = ',ok_newmicro 
    WRITE(lunout,*) ' ratqsbas = ',ratqsbas 
    WRITE(lunout,*) ' ratqshaut = ',ratqshaut 
    WRITE(lunout,*) ' tau_ratqs = ',tau_ratqs 
    WRITE(lunout,*) ' top_height = ',top_height 
    WRITE(lunout,*) ' rad_froid = ',rad_froid
    WRITE(lunout,*) ' rad_chau1 = ',rad_chau1
    WRITE(lunout,*) ' rad_chau2 = ',rad_chau2
    WRITE(lunout,*) ' t_glace_min = ',t_glace_min
    WRITE(lunout,*) ' t_glace_max = ',t_glace_max
    WRITE(lunout,*) ' exposant_glace = ',exposant_glace
    WRITE(lunout,*) ' iflag_t_glace = ',iflag_t_glace
    WRITE(lunout,*) ' iflag_cloudth_vert = ',iflag_cloudth_vert
    WRITE(lunout,*) ' iflag_rain_incloud_vol = ',iflag_rain_incloud_vol
    WRITE(lunout,*) ' iflag_ice_thermo = ',iflag_ice_thermo
    WRITE(lunout,*) ' rei_min = ',rei_min
    WRITE(lunout,*) ' rei_max = ',rei_max
    WRITE(lunout,*) ' overlap = ',overlap 
    WRITE(lunout,*) ' cdmmax = ',cdmmax 
    WRITE(lunout,*) ' cdhmax = ',cdhmax 
    WRITE(lunout,*) ' ksta = ',ksta 
    WRITE(lunout,*) ' ksta_ter = ',ksta_ter 
    WRITE(lunout,*) ' f_ri_cd_min = ',f_ri_cd_min 
    WRITE(lunout,*) ' ok_kzmin = ',ok_kzmin 
    WRITE(lunout,*) ' pbl_lmixmin_alpha = ',pbl_lmixmin_alpha
    WRITE(lunout,*) ' fmagic = ',fmagic
    WRITE(lunout,*) ' pmagic = ',pmagic
    WRITE(lunout,*) ' ok_ade = ',ok_ade
    WRITE(lunout,*) ' ok_volcan = ',ok_volcan
    WRITE(lunout,*) ' ok_aie = ',ok_aie
    WRITE(lunout,*) ' ok_alw = ',ok_alw
    WRITE(lunout,*) ' aerosol_couple = ', aerosol_couple
    WRITE(lunout,*) ' chemistry_couple = ', chemistry_couple
    WRITE(lunout,*) ' flag_aerosol = ', flag_aerosol
    WRITE(lunout,*) ' flag_aerosol_strat= ', flag_aerosol_strat
    WRITE(lunout,*) ' flag_aer_feedback= ', flag_aer_feedback
    WRITE(lunout,*) ' aer_type = ',aer_type
    WRITE(lunout,*) ' bl95_b0 = ',bl95_b0
    WRITE(lunout,*) ' bl95_b1 = ',bl95_b1
    WRITE(lunout,*) ' lev_histhf = ',lev_histhf 
    WRITE(lunout,*) ' lev_histday = ',lev_histday 
    WRITE(lunout,*) ' lev_histmth = ',lev_histmth 
    WRITE(lunout,*) ' lev_histins = ',lev_histins
    WRITE(lunout,*) ' lev_histLES = ',lev_histLES
    WRITE(lunout,*) ' lev_histdayNMC = ',lev_histdayNMC
    WRITE(lunout,*) ' levout_histNMC = ',levout_histNMC
    WRITE(lunout,*) ' ok_histNMC = ',ok_histNMC
    WRITE(lunout,*) ' freq_outNMC = ',freq_outNMC
    WRITE(lunout,*) ' freq_calNMC = ',freq_calNMC
    WRITE(lunout,*) ' iflag_pbl = ', iflag_pbl
!FC
    WRITE(lunout,*) ' ifl_pbltree = ', ifl_pbltree
    WRITE(lunout,*) ' Cd_frein = ', Cd_frein
    WRITE(lunout,*) ' iflag_pbl_split = ', iflag_pbl_split
    WRITE(lunout,*) ' iflag_order2_sollw = ', iflag_order2_sollw
    WRITE(lunout,*) ' iflag_thermals = ', iflag_thermals
    WRITE(lunout,*) ' iflag_thermals_ed = ', iflag_thermals_ed
    WRITE(lunout,*) ' fact_thermals_ed_dz = ', fact_thermals_ed_dz
    WRITE(lunout,*) ' iflag_thermals_optflux = ', iflag_thermals_optflux
    WRITE(lunout,*) ' iflag_thermals_closure = ', iflag_thermals_closure
    WRITE(lunout,*) ' iflag_clos = ', iflag_clos
    WRITE(lunout,*) ' coef_clos_ls = ', coef_clos_ls
    WRITE(lunout,*) ' type_run = ',type_run 
    WRITE(lunout,*) ' ok_cosp = ',ok_cosp
    WRITE(lunout,*) ' ok_airs = ',ok_airs

    WRITE(lunout,*) ' ok_mensuelCOSP = ',ok_mensuelCOSP
    WRITE(lunout,*) ' ok_journeCOSP = ',ok_journeCOSP
    WRITE(lunout,*) ' ok_hfCOSP =',ok_hfCOSP
    WRITE(lunout,*) ' solarlong0 = ', solarlong0
    WRITE(lunout,*) ' qsol0 = ', qsol0
    WRITE(lunout,*) ' evap0 = ', evap0
    WRITE(lunout,*) ' albsno0 = ', albsno0
    WRITE(lunout,*) ' iflag_sic = ', iflag_sic
    WRITE(lunout,*) ' inertie_sol = ', inertie_sol
    WRITE(lunout,*) ' inertie_sic = ', inertie_sic
    WRITE(lunout,*) ' inertie_lic = ', inertie_lic
    WRITE(lunout,*) ' inertie_sno = ', inertie_sno
    WRITE(lunout,*) ' f_cdrag_ter = ',f_cdrag_ter
    WRITE(lunout,*) ' f_cdrag_oce = ',f_cdrag_oce
    WRITE(lunout,*) ' f_rugoro = ',f_rugoro
    WRITE(lunout,*) ' z0min = ',z0min
    WRITE(lunout,*) ' supcrit1 = ', supcrit1
    WRITE(lunout,*) ' supcrit2 = ', supcrit2
    WRITE(lunout,*) ' iflag_mix = ', iflag_mix
    WRITE(lunout,*) ' iflag_mix_adiab = ', iflag_mix_adiab
    WRITE(lunout,*) ' scut = ', scut
    WRITE(lunout,*) ' qqa1 = ', qqa1
    WRITE(lunout,*) ' qqa2 = ', qqa2
    WRITE(lunout,*) ' gammas = ', gammas
    WRITE(lunout,*) ' Fmax = ', Fmax
    WRITE(lunout,*) ' tmax_fonte_cv = ', tmax_fonte_cv
    WRITE(lunout,*) ' alphas = ', alphas
    WRITE(lunout,*) ' iflag_wake = ', iflag_wake
    WRITE(lunout,*) ' alp_offset = ', alp_offset
    ! nrlmd le 10/04/2012
    WRITE(lunout,*) ' iflag_trig_bl = ', iflag_trig_bl
    WRITE(lunout,*) ' s_trig = ', s_trig
    WRITE(lunout,*) ' tau_trig_shallow = ', tau_trig_shallow
    WRITE(lunout,*) ' tau_trig_deep = ', tau_trig_deep
    WRITE(lunout,*) ' iflag_clos_bl = ', iflag_clos_bl
    ! fin nrlmd le 10/04/2012

    WRITE(lunout,*) ' lonmin lonmax latmin latmax bilKP_ins =',&
         lonmin_ins, lonmax_ins, latmin_ins, latmax_ins
    WRITE(lunout,*) ' ecrit_ hf, ins, day, mth, reg, tra, ISCCP, LES',&
         ecrit_hf, ecrit_ins, ecrit_day, ecrit_mth, ecrit_reg, ecrit_tra, ecrit_ISCCP, ecrit_LES

    WRITE(lunout,*) ' ok_strato = ', ok_strato
    WRITE(lunout,*) ' ok_hines = ',  ok_hines
    WRITE(lunout,*) ' ok_gwd_rando = ',  ok_gwd_rando
    WRITE(lunout,*) ' ok_qch4 = ',  ok_qch4
    WRITE(lunout,*) ' gwd_rando_ruwmax = ', gwd_rando_ruwmax
    WRITE(lunout,*) ' gwd_rando_sat = ', gwd_rando_sat
    WRITE(lunout,*) ' gwd_front_ruwmax = ', gwd_front_ruwmax
    WRITE(lunout,*) ' gwd_front_sat = ', gwd_front_sat
    WRITE(lunout,*) ' SSO gkdrag =',gkdrag
    WRITE(lunout,*) ' SSO grahilo=',grahilo
    WRITE(lunout,*) ' SSO grcrit=',grcrit
    WRITE(lunout,*) ' SSO gfrcrit=',gfrcrit
    WRITE(lunout,*) ' SSO gkwake=',gkwake
    WRITE(lunout,*) ' SSO gklift=',gklift
    WRITE(lunout,*) ' adjust_tropopause = ', adjust_tropopause
    WRITE(lunout,*) ' ok_daily_climoz = ',ok_daily_climoz
    WRITE(lunout,*) ' read_climoz = ', read_climoz
    WRITE(lunout,*) ' carbon_cycle_tr = ', carbon_cycle_tr
    WRITE(lunout,*) ' carbon_cycle_cpl = ', carbon_cycle_cpl
    WRITE(lunout,*) ' carbon_cycle_rad = ', carbon_cycle_rad
    WRITE(lunout,*) ' level_coupling_esm = ', level_coupling_esm

    !$OMP END MASTER

    RETURN

  END SUBROUTINE conf_phys

END MODULE conf_phys_m
!
!#################################################################
!

SUBROUTINE conf_interface(tau_calv)

  USE IOIPSL
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE
  ! Configuration de l'interace atm/surf
  !
  ! tau_calv:    temps de relaxation pour la fonte des glaciers
  !
  REAL          :: tau_calv
  REAL, SAVE    :: tau_calv_omp
  !
  !Config Key  = tau_calv
  !Config Desc = temps de relaxation pour fonte des glaciers en jours
  !Config Def  = 1 an 
  !Config Help = 
  !
  tau_calv_omp = 360.*10.
  !$OMP MASTER
  CALL getin('tau_calv',tau_calv_omp)
  !$OMP END MASTER
  !$OMP BARRIER
  !
  tau_calv=tau_calv_omp
  !
  !$OMP MASTER
  WRITE(lunout,*)' ##############################################'
  WRITE(lunout,*)' Configuration de l''interface atm/surfaces  : '
  WRITE(lunout,*)' tau_calv = ',tau_calv
  !$OMP END MASTER
  !
  RETURN

END SUBROUTINE conf_interface
