!
! $Id: phys_output_write_mod.F90 2298 2015-06-14 19:13:32Z fairhead $
!
MODULE phys_output_write_spl_mod

!JE20150620<<
!JE20150620>>
!JE20150620<<

  USE time_phylmdz_mod, ONLY: day_step_phy, start_time, itau_phy

  USE phytracr_spl_mod, ONLY : ok_chimeredust, id_prec, id_fine, id_coss, &
       id_codu, id_scdu , &
       d_tr_cl, d_tr_th, d_tr_cv, d_tr_lessi_impa, &
       d_tr_lessi_nucl, d_tr_insc, d_tr_bcscav, d_tr_evapls, d_tr_ls,  &
       d_tr_trsp, d_tr_sscav, d_tr_sat, d_tr_uscav ,&
       diff_aod550_tot,&
       diag_aod670_tot, diag_aod865_tot, &
       diff_aod550_tr2, diag_aod670_tr2, diag_aod865_tr2, &
       diag_aod550_ss, diag_aod670_ss, diag_aod865_ss, &
       diag_aod550_dust, diag_aod670_dust, diag_aod865_dust , &
       diag_aod550_dustsco, diag_aod670_dustsco, diag_aod865_dustsco, &
!       aod550_aqua, aod670_aqua, aod865_aqua, &
!       aod550_terra, aod670_terra, aod865_terra, &
       aod550_aqua,aod550_tr2_aqua,aod550_ss_aqua,aod550_dust_aqua,aod550_dustsco_aqua,&
       aod670_aqua,aod670_tr2_aqua,aod670_ss_aqua,aod670_dust_aqua,aod670_dustsco_aqua,&
       aod865_aqua,aod865_tr2_aqua,aod865_ss_aqua,aod865_dust_aqua,aod865_dustsco_aqua,&
       aod550_terra,aod550_tr2_terra,aod550_ss_terra,aod550_dust_terra,aod550_dustsco_terra,&
       aod670_terra,aod670_tr2_terra,aod670_ss_terra,aod670_dust_terra,aod670_dustsco_terra,&
       aod865_terra,aod865_tr2_terra,aod865_ss_terra,aod865_dust_terra,aod865_dustsco_terra,&
       trm01,trm02,trm03,trm04,trm05, &
       sconc01,sconc02,sconc03,sconc04,sconc05, &
       flux01,flux02,flux03,flux04,flux05,&
       ds01,ds02,ds03,ds04,ds05, &
       dh01,dh02,dh03,dh04,dh05, &
       dtrconv01,dtrconv02,dtrconv03,dtrconv04,dtrconv05, &
       dtherm01,dtherm02,dtherm03,dtherm04,dtherm05, &
       dhkecv01,dhkecv02,dhkecv03,dhkecv04,dhkecv05, &
       d_tr_ds01,d_tr_ds02,d_tr_ds03,d_tr_ds04,d_tr_ds05, &
       dhkelsc01,dhkelsc02,dhkelsc03,dhkelsc04,dhkelsc05, &
       d_tr_cv01,d_tr_cv02,d_tr_cv03,d_tr_cv04,d_tr_cv05, &
       d_tr_trsp01,d_tr_trsp02,d_tr_trsp03,d_tr_trsp04,d_tr_trsp05, &
       d_tr_sscav01,d_tr_sscav02,d_tr_sscav03,d_tr_sscav04,d_tr_sscav05, &
       d_tr_sat01,d_tr_sat02,d_tr_sat03,d_tr_sat04,d_tr_sat05, &
       d_tr_uscav01,d_tr_uscav02,d_tr_uscav03,d_tr_uscav04,d_tr_uscav05, &
       d_tr_insc01,d_tr_insc02,d_tr_insc03,d_tr_insc04,d_tr_insc05, &
       d_tr_bcscav01,d_tr_bcscav02,d_tr_bcscav03,d_tr_bcscav04,d_tr_bcscav05, &
       d_tr_evapls01,d_tr_evapls02,d_tr_evapls03,d_tr_evapls04,d_tr_evapls05, &
       d_tr_ls01,d_tr_ls02,d_tr_ls03,d_tr_ls04,d_tr_ls05, &
       d_tr_dyn01,d_tr_dyn02,d_tr_dyn03,d_tr_dyn04,d_tr_dyn05, &
       d_tr_cl01,d_tr_cl02,d_tr_cl03,d_tr_cl04,d_tr_cl05, &
       d_tr_th01,d_tr_th02,d_tr_th03,d_tr_th04,d_tr_th05, &
       sed_ss,sed_dust,sed_dustsco,his_g2pgas,his_g2paer, &
       sed_ss3D,sed_dust3D,sed_dustsco3D, &
       fluxbb, &
       fluxff,fluxbcbb,fluxbcff,fluxbcnff, &
       fluxbcba,fluxbc,fluxombb,fluxomff,fluxomnff, &
       fluxomba,fluxomnat,fluxom,fluxh2sff,fluxh2snff, &
       fluxso2ff,fluxso2nff,fluxso2bb,fluxso2vol,fluxso2ba, &
       fluxso2,fluxso4ff,fluxso4nff,fluxso4ba,fluxso4bb, &
       fluxso4,fluxdms,fluxh2sbio,fluxdustec,&
       fluxddfine,  &
       fluxddcoa,fluxddsco,fluxdd, &
       fluxssfine,fluxsscoa, &
       fluxss,flux_sparam_ind,flux_sparam_bb,flux_sparam_ff, &
       flux_sparam_ddfine,flux_sparam_ddcoa, &
       flux_sparam_ddsco,flux_sparam_ssfine, &
       flux_sparam_sscoa,u10m_ss,v10m_ss

  USE dustemission_mod, ONLY : m1dflux, m2dflux, m3dflux

!  USE phytrac_mod, ONLY : d_tr_cl, d_tr_th, d_tr_cv, d_tr_lessi_impa, &
!       d_tr_lessi_nucl, d_tr_insc, d_tr_bcscav, d_tr_evapls, d_tr_ls,  &
!       d_tr_trsp, d_tr_sscav, d_tr_sat, d_tr_uscav

!JE20150620>>

  ! Author: Abderrahmane IDELKADI (original include file)
  ! Author: Laurent FAIRHEAD (transformation to module/subroutine)
  ! Author: Ulysse GERARD (effective implementation)

CONTAINS 

  ! ug Routine pour définir (los du premier passageà) ET sortir les variables
  SUBROUTINE phys_output_write_spl(itap, pdtphys, paprs, pphis, &
       pplay, lmax_th, aerosol_couple,         &
       ok_ade, ok_aie, ivap, ok_sync, &
       ptconv, read_climoz, clevSTD, ptconvth, &
       d_t, qx, d_qx, d_tr_dyn, zmasse, flag_aerosol, flag_aerosol_strat, ok_cdnc)

    ! This subroutine does the actual writing of diagnostics that were
    ! defined and initialised in phys_output_mod.F90

    USE dimphy, ONLY: klon, klev, klevp1
    USE ocean_slab_mod, ONLY: nslay
    USE control_mod, ONLY: day_step, iphysiq
    USE phys_output_ctrlout_mod, ONLY: o_phis, o_aire, is_ter, is_lic, is_oce, &
         is_ave, is_sic, o_contfracATM, o_contfracOR, &
         o_aireTER, o_flat, o_slp, o_tsol, &
         o_t2m, o_t2m_min, o_t2m_max, &
         o_t2m_min_mon, o_t2m_max_mon, &
         o_q2m, o_ustar, o_u10m, o_v10m, &
         o_wind10m, o_wind10max, o_gusts, o_sicf, &
         o_psol, o_mass, o_qsurf, o_qsol, &
         o_precip, o_ndayrain, o_plul, o_pluc, &
         o_snow, o_msnow, o_fsnow, o_evap, &
         o_tops, o_tops0, o_topl, o_topl0, &
         o_SWupTOA, o_SWupTOAclr, o_SWdnTOA, &
         o_SWdnTOAclr, o_nettop, o_SWup200, &
         o_SWup200clr, o_SWdn200, o_SWdn200clr, &
         o_LWup200, o_LWup200clr, o_LWdn200, &
         o_LWdn200clr, o_sols, o_sols0, &
         o_soll, o_radsol, o_soll0, o_SWupSFC, &
         o_SWupSFCclr, o_SWdnSFC, o_SWdnSFCclr, &
         o_LWupSFC, o_LWdnSFC, o_LWupSFCclr, &
         o_LWdnSFCclr, o_bils, o_bils_diss, &
         o_bils_ec,o_bils_ech, o_bils_tke, o_bils_kinetic, &
         o_bils_latent, o_bils_enthalp, o_sens, &
         o_fder, o_ffonte, o_fqcalving, o_fqfonte, &
         o_taux, o_tauy, o_snowsrf, o_qsnow, &
         o_snowhgt, o_toice, o_sissnow, o_runoff, &
         o_albslw3, o_pourc_srf, o_fract_srf, &
         o_taux_srf, o_tauy_srf, o_tsol_srf, &
         o_evappot_srf, o_ustar_srf, o_u10m_srf, &
         o_v10m_srf, o_t2m_srf, o_evap_srf, &
         o_sens_srf, o_lat_srf, o_flw_srf, &
         o_fsw_srf, o_wbils_srf, o_wbilo_srf, &
         o_tke_srf, o_tke_max_srf,o_dltpbltke_srf, o_wstar, &
         o_cdrm, o_cdrh, o_cldl, o_cldm, o_cldh, &
         o_cldt, o_JrNt, o_cldljn, o_cldmjn, &
         o_cldhjn, o_cldtjn, o_cldq, o_lwp, o_iwp, &
         o_ue, o_ve, o_uq, o_vq, o_cape, o_pbase, &
         o_ptop, o_fbase, o_plcl, o_plfc, &
         o_wbeff, o_cape_max, o_upwd, o_Ma, &
         o_dnwd, o_dnwd0, o_ftime_con, o_mc, &
         o_prw, o_s_pblh, o_s_pblt, o_s_lcl, &
         o_s_therm, o_uSTDlevs, o_vSTDlevs, &
         o_wSTDlevs, o_zSTDlevs, o_qSTDlevs, &
         o_tSTDlevs, epsfra, o_t_oce_sic, &
         o_ale_bl, o_alp_bl, o_ale_wk, o_alp_wk, &
         o_ale, o_alp, o_cin, o_WAPE, o_wake_h, &
         o_wake_s, o_wake_deltat, o_wake_deltaq, &
         o_wake_omg, o_dtwak, o_dqwak, o_Vprecip, &
         o_ftd, o_fqd, o_wdtrainA, o_wdtrainM, &
         o_n2, o_s2, o_proba_notrig, &
         o_random_notrig, o_ale_bl_stat, &
         o_ale_bl_trig, o_alp_bl_det, &
         o_alp_bl_fluct_m, o_alp_bl_fluct_tke, &
         o_alp_bl_conv, o_alp_bl_stat, &
         o_slab_qflux, o_tslab, o_slab_bils, &
         o_slab_bilg, o_slab_sic, o_slab_tice, &
         o_weakinv, o_dthmin, o_cldtau, &
         o_cldemi, o_pr_con_l, o_pr_con_i, &
         o_pr_lsc_l, o_pr_lsc_i, o_re, o_fl, &
         o_rh2m, o_rh2m_min, o_rh2m_max, &
         o_qsat2m, o_tpot, o_tpote, o_SWnetOR, &
         o_LWdownOR, o_snowl, &
         o_solldown, o_dtsvdfo, o_dtsvdft, &
         o_dtsvdfg, o_dtsvdfi, o_z0m, o_z0h, o_od550aer, &
         o_od865aer, o_absvisaer, o_od550lt1aer, &
         o_sconcso4, o_sconcno3, o_sconcoa, o_sconcbc, &
         o_sconcss, o_sconcdust, o_concso4, o_concno3, &
         o_concoa, o_concbc, o_concss, o_concdust, &
         o_loadso4, o_loadoa, o_loadbc, o_loadss, &
         o_loaddust, o_tausumaero, o_tausumaero_lw, &
         o_topswad, o_topswad0, o_solswad, o_solswad0, &
         o_toplwad, o_toplwad0, o_sollwad, o_sollwad0, &
         o_swtoaas_nat, o_swsrfas_nat, &
         o_swtoacs_nat, o_swtoaas_ant, &
         o_swsrfas_ant, o_swtoacs_ant, &
         o_swsrfcs_ant, o_swtoacf_nat, &
         o_swsrfcf_nat, o_swtoacf_ant, &
         o_swsrfcs_nat, o_swsrfcf_ant, &
         o_swtoacf_zero, o_swsrfcf_zero, &
         o_topswai, o_solswai, o_scdnc, &
         o_cldncl, o_reffclws, o_reffclwc, &
         o_cldnvi, o_lcc, o_lcc3d, o_lcc3dcon, &
         o_lcc3dstra, o_reffclwtop, o_ec550aer, &
         o_lwcon, o_iwcon, o_temp, o_theta, &
         o_ovapinit, o_ovap, o_oliq, o_geop, &
         o_vitu, o_vitv, o_vitw, o_pres, o_paprs, &
         o_zfull, o_zhalf, o_rneb, o_rnebjn, o_rnebcon, &
         o_rnebls, o_rhum, o_ozone, o_ozone_light, &
         o_dtphy, o_dqphy, o_albe_srf, o_z0m_srf, o_z0h_srf, &
         o_ages_srf, o_snow_srf, o_alb1, o_alb2, o_tke, &
         o_tke_max, o_kz, o_kz_max, o_clwcon, &
         o_dtdyn, o_dqdyn, o_dudyn, o_dvdyn, &
         o_dtcon, o_tntc, o_ducon, o_dvcon, &
         o_dqcon, o_tnhusc, o_tnhusc, o_dtlsc, &
         o_dtlschr, o_dqlsc, o_beta_prec, &
         o_dtlscth, o_dtlscst, o_dqlscth, &
         o_dqlscst, o_plulth, o_plulst, &
         o_ptconvth, o_lmaxth, o_dtvdf, &
         o_dtdis, o_dqvdf, o_dteva, o_dqeva, &
         o_ptconv, o_ratqs, o_dtthe, & 
         o_duthe, o_dvthe, o_ftime_th, &
         o_f_th, o_e_th, o_w_th, o_q_th, &
         o_a_th, o_d_th, o_f0_th, o_zmax_th, &
         o_dqthe, o_dtajs, o_dqajs, o_dtswr, &
         o_dtsw0, o_dtlwr, o_dtlw0, o_dtec, &
         o_duvdf, o_dvvdf, o_duoro, o_dvoro, &
         o_dtoro, o_dulif, o_dvlif, o_dtlif, &
 !       o_duhin, o_dvhin, o_dthin, &
         o_dqch4, o_rsu, &
         o_rsd, o_rlu, o_rld, o_rsucs, o_rsdcs, &
         o_rlucs, o_rldcs, o_tnt, o_tntr, &
         o_tntscpbl, o_tnhus, o_tnhusscpbl, &
         o_evu, o_h2o, o_mcd, o_dmc, o_ref_liq, &
         o_ref_ice, o_rsut4co2, o_rlut4co2, &
         o_rsutcs4co2, o_rlutcs4co2, o_rsu4co2, &
         o_rlu4co2, o_rsucs4co2, o_rlucs4co2, &
         o_rsd4co2, o_rld4co2, o_rsdcs4co2, &
         o_rldcs4co2, o_tnondef, o_ta, o_zg, &
         o_hus, o_hur, o_ua, o_va, o_wap, &
         o_psbg, o_tro3, o_tro3_daylight, &
         o_uxv, o_vxq, o_vxT, o_wxq, o_vxphi, &
         o_wxT, o_uxu, o_vxv, o_TxT, o_trac, &
         o_dtr_vdf, o_dtr_the, o_dtr_con, &
         o_dtr_lessi_impa, o_dtr_lessi_nucl, &
         o_dtr_insc, o_dtr_bcscav, o_dtr_evapls, &
!        o_dtr_ls, o_dtr_dyn, o_dtr_cl, o_dtr_trsp, o_dtr_sscav, &
         o_dtr_ls, o_dtr_trsp, o_dtr_sscav, &
         o_dtr_sat, o_dtr_uscav, o_trac_cum, o_du_gwd_rando, o_dv_gwd_rando, &
!JE20150620<<
!         o_vstr_gwd_rando
         o_vstr_gwd_rando, &
         o_m1dflux,o_m2dflux,o_m3dflux, &
         o_taue550, &
         o_taue670,o_taue865, &
         o_taue550_tr2, o_taue670_tr2, o_taue865_tr2, &
         o_taue550_ss,o_taue670_ss, o_taue865_ss, &
         o_taue550_dust, o_taue670_dust, o_taue865_dust, &
         o_taue550_dustsco, o_taue670_dustsco, o_taue865_dustsco, &
         o_taue550_aqua, o_taue670_aqua, o_taue865_aqua, &
         o_taue550_terra, o_taue670_terra, o_taue865_terra, &
         o_taue550_fine_aqua     ,         o_taue670_fine_aqua     ,  &
         o_taue865_fine_aqua     ,         o_taue550_coss_aqua      ,  &
         o_taue670_coss_aqua      ,         o_taue865_coss_aqua      ,  &
         o_taue550_codu_aqua    ,         o_taue670_codu_aqua    ,  &
         o_taue865_codu_aqua    ,         o_taue670_scdu_aqua ,  &
         o_taue550_scdu_aqua ,         o_taue865_scdu_aqua ,  &
         o_taue550_fine_terra     ,         o_taue670_fine_terra     ,&
         o_taue865_fine_terra     ,         o_taue550_coss_terra      ,&
         o_taue670_coss_terra      ,         o_taue865_coss_terra      ,&
         o_taue550_codu_terra    ,         o_taue670_codu_terra    ,&
         o_taue865_codu_terra    ,         o_taue670_scdu_terra ,&
         o_taue550_scdu_terra ,         o_taue865_scdu_terra ,&
         o_trm01,o_trm02,o_trm03,o_trm04,o_trm05,&
         o_sconc01,o_sconc02,o_sconc03,o_sconc04,o_sconc05, &
         o_flux01,o_flux02,o_flux03,o_flux04,o_flux05, &
         o_ds01,o_ds02,o_ds03,o_ds04,o_ds05, &
         o_dh01,o_dh02,o_dh03,o_dh04,o_dh05, &
         o_dtrconv01,o_dtrconv02,o_dtrconv03,o_dtrconv04,o_dtrconv05, &
         o_dtherm01,o_dtherm02,o_dtherm03,o_dtherm04,o_dtherm05, &
         o_dhkecv01,o_dhkecv02,o_dhkecv03,o_dhkecv04,o_dhkecv05, &
         o_d_tr_ds01,o_d_tr_ds02,o_d_tr_ds03,o_d_tr_ds04,o_d_tr_ds05, &
         o_dhkelsc01,o_dhkelsc02,o_dhkelsc03,o_dhkelsc04,o_dhkelsc05, &
         o_d_tr_sat01,o_d_tr_cv01,o_d_tr_cv02,o_d_tr_cv03,o_d_tr_cv04,o_d_tr_cv05,&
         o_d_tr_trsp01,o_d_tr_trsp02,o_d_tr_trsp03,o_d_tr_trsp04,o_d_tr_trsp05,&
         o_d_tr_sscav01,o_d_tr_sscav02,o_d_tr_sscav03,o_d_tr_sscav04,o_d_tr_sscav05,&
         o_d_tr_sat02,o_d_tr_sat03,o_d_tr_sat04,o_d_tr_sat05,  &
         o_d_tr_uscav01,o_d_tr_uscav02,o_d_tr_uscav03,o_d_tr_uscav04,o_d_tr_uscav05,&
         o_d_tr_insc01,o_d_tr_insc02,o_d_tr_insc03,o_d_tr_insc04,o_d_tr_insc05,&
         o_d_tr_bcscav01,o_d_tr_bcscav02,o_d_tr_bcscav03,o_d_tr_bcscav04,o_d_tr_bcscav05,&
         o_d_tr_evapls01,o_d_tr_evapls02,o_d_tr_evapls03,o_d_tr_evapls04,o_d_tr_evapls05,&
         o_d_tr_ls01,o_d_tr_ls02,o_d_tr_ls03,o_d_tr_ls04,o_d_tr_ls05,&
         o_d_tr_dyn01,o_d_tr_dyn02,o_d_tr_dyn03,o_d_tr_dyn04,o_d_tr_dyn05,&
         o_d_tr_cl01,o_d_tr_cl02,o_d_tr_cl03,o_d_tr_cl04,o_d_tr_cl05,&
         o_d_tr_th01,o_d_tr_th02,o_d_tr_th03,o_d_tr_th04,o_d_tr_th05,&
         o_sed_ss,o_sed_dust,o_sed_dustsco,o_g2p_gas,o_g2p_aer, &
         o_sed_ss3D,o_sed_dust3D,o_sed_dustsco3D, &
         o_fluxbb, &
         o_fluxff    ,o_fluxbcbb  ,o_fluxbcff  ,o_fluxbcnff , &
         o_fluxbcba  ,o_fluxbc    ,o_fluxombb  ,o_fluxomff  , &
         o_fluxomnff ,o_fluxomba  ,o_fluxomnat ,o_fluxom    , &
         o_fluxh2sff ,o_fluxh2snff,o_fluxso2ff ,o_fluxso2nff, &
         o_fluxso2bb ,o_fluxso2vol,o_fluxso2ba ,o_fluxso2   , &
         o_fluxso4ff ,o_fluxso4nff,o_fluxso4bb ,o_fluxso4ba , &
         o_fluxso4   ,o_fluxdms   ,o_fluxh2sbio,o_fluxdustec, &
         o_fluxddfine,o_fluxddcoa ,o_fluxddsco ,o_fluxdd    ,&
         o_fluxssfine,o_fluxsscoa, o_fluxss, &
         o_flux_sparam_ind,o_flux_sparam_bb, &
         o_flux_sparam_ff ,o_flux_sparam_ddfine  ,o_flux_sparam_ddcoa, &
         o_flux_sparam_ddsco,o_flux_sparam_ssfine,o_flux_sparam_sscoa, &
         o_u10m_ss,o_v10m_ss

!JE20150620>>

    USE phys_state_var_mod, ONLY: pctsrf, paire_ter, rain_fall, snow_fall, &
         qsol, z0m, z0h, fevap, agesno, &
         nday_rain, rain_con, snow_con, &
         topsw, toplw, toplw0, swup, swdn, &
         topsw0, swup0, swdn0, SWup200, SWup200clr, &
         SWdn200, SWdn200clr, LWup200, LWup200clr, &
         LWdn200, LWdn200clr, solsw, solsw0, sollw, &
         radsol, sollw0, sollwdown, sollw, gustiness, &
         sollwdownclr, lwdn0, ftsol, ustar, u10m, &
         v10m, pbl_tke, wake_delta_pbl_TKE, &
         wstar, cape, ema_pcb, ema_pct, &
         ema_cbmf, Ma, fm_therm, ale_bl, alp_bl, ale, &
         alp, cin, wake_pe, wake_s, wake_deltat, &
         wake_deltaq, ftd, fqd, ale_bl_trig, albsol1, &
         rnebcon, wo, falb1, albsol2, coefh, clwcon0, &
         ratqs, entr_therm, zqasc, detr_therm, f0, &
         lwup, lwdn, lwup0, coefm, &
         swupp, lwupp, swup0p, lwup0p, swdnp, lwdnp, &
         swdn0p, lwdn0p, tnondef, O3sumSTD, uvsumSTD, &
         vqsumSTD, vTsumSTD, O3daysumSTD, wqsumSTD, &
         vphisumSTD, wTsumSTD, u2sumSTD, v2sumSTD, &
         T2sumSTD, nlevSTD, &
!        du_gwd_rando, dv_gwd_rando, &
         ulevSTD, vlevSTD, wlevSTD, philevSTD, qlevSTD, tlevSTD, &
         rhlevSTD, O3STD, O3daySTD, uvSTD, vqSTD, vTSTD, wqSTD, &
         vphiSTD, wTSTD, u2STD, v2STD, T2STD, missing_val_nf90

    USE phys_local_var_mod, ONLY: zxfluxlat, slp, zxtsol, zt2m, &
         t2m_min_mon, t2m_max_mon, evap, &
         zu10m, zv10m, zq2m, zustar, zxqsurf, &
         rain_lsc, snow_lsc, bils, sens, fder, &
         zxffonte, zxfqcalving, zxfqfonte, fluxu, &
         fluxv, zxsnow, qsnow, snowhgt, to_ice, &
         sissnow, runoff, albsol3_lic, evap_pot, &
         t2m, fluxt, fluxlat, fsollw, fsolsw, &
         wfbils, wfbilo, cdragm, cdragh, cldl, cldm, &
         cldh, cldt, JrNt, cldljn, cldmjn, cldhjn, &
         cldtjn, cldq, flwp, fiwp, ue, ve, uq, vq, &
         plcl, plfc, wbeff, upwd, dnwd, dnwd0, prw, &
         s_pblh, s_pblt, s_lcl, s_therm, uwriteSTD, &
         vwriteSTD, wwriteSTD, phiwriteSTD, qwriteSTD, &
         twriteSTD, ale_wake, alp_wake, wake_h, &
         wake_omg, d_t_wake, d_q_wake, Vprecip, &
         wdtrainA, wdtrainM, n2, s2, proba_notrig, &
         random_notrig, ale_bl_stat, &
         alp_bl_det, alp_bl_fluct_m, alp_bl_conv, &
         alp_bl_stat, alp_bl_fluct_tke, slab_wfbils, &
         weak_inversion, dthmin, cldtau, cldemi, &
         pmflxr, pmflxs, prfl, psfl, re, fl, rh2m, &
         qsat2m, tpote, tpot, d_ts, od550aer, &
         od865aer, absvisaer, od550lt1aer, sconcso4, sconcno3, &
         sconcoa, sconcbc, sconcss, sconcdust, concso4, concno3, &
         concoa, concbc, concss, concdust, loadso4, &
         loadoa, loadbc, loadss, loaddust, tausum_aero, &
         topswad_aero, topswad0_aero, solswad_aero, &
         solswad0_aero, topsw_aero, solsw_aero, &
         topsw0_aero, solsw0_aero, topswcf_aero, &
         solswcf_aero, topswai_aero, solswai_aero, &
         toplwad_aero, toplwad0_aero, sollwad_aero, &
         sollwad0_aero, toplwai_aero, sollwai_aero, &
         scdnc, cldncl, reffclws, reffclwc, cldnvi, &
         lcc, lcc3d, lcc3dcon, lcc3dstra, reffclwtop, &
         ec550aer, flwc, fiwc, t_seri, theta, q_seri, &
!jyg<
!!         ql_seri, zphi, u_seri, v_seri, omega, cldfra, &
         ql_seri, tr_seri, &
         zphi, u_seri, v_seri, omega, cldfra, &
!>jyg
         rneb, rnebjn, zx_rh, d_t_dyn, d_q_dyn, &
         d_u_dyn, d_v_dyn, d_t_con, d_t_ajsb, d_t_ajs, &
         d_u_ajs, d_v_ajs, &
         d_u_con, d_v_con, d_q_con, d_q_ajs, d_t_lsc, &
         d_t_lwr,d_t_lw0,d_t_swr,d_t_sw0, &
         d_t_eva, d_q_lsc, beta_prec, d_t_lscth, &
         d_t_lscst, d_q_lscth, d_q_lscst, plul_th, &
         plul_st, d_t_vdf, d_t_diss, d_q_vdf, d_q_eva, &
         zw2, fraca, zmax_th, d_q_ajsb, d_t_ec, d_u_vdf, &
         d_v_vdf, d_u_oro, d_v_oro, d_t_oro, d_u_lif, &
         d_v_lif, d_t_lif, &
!        d_u_hin, d_v_hin, d_t_hin, &
         d_q_ch4, pmfd, pmfu, ref_liq, ref_ice, rhwriteSTD

    USE phys_output_var_mod, ONLY: vars_defined, snow_o, zfra_o, bils_diss, &
         bils_ec,bils_ech, bils_tke, bils_kinetic, bils_latent, bils_enthalp, &
         itau_con, nfiles, clef_files, nid_files, zvstr_gwd_rando
    USE ocean_slab_mod, ONLY: tslab, slab_bils, slab_bilg, tice, seaice
    USE pbl_surface_mod, ONLY: snow
    USE indice_sol_mod, ONLY: nbsrf
    USE infotrac, ONLY: nqtot, nqo, nbtr, type_trac
    USE geometry_mod, ONLY: cell_area
    USE surface_data, ONLY: type_ocean, version_ocean, ok_veget, ok_snow
!    USE aero_mod, ONLY: naero_spc
    USE aero_mod, ONLY: naero_tot, id_STRAT_phy
    USE ioipsl, ONLY: histend, histsync
    USE iophy, ONLY: set_itau_iophy, histwrite_phy
    USE netcdf, ONLY: nf90_fill_real

#ifdef CPP_XIOS
    ! ug Pour les sorties XIOS
    USE xios, ONLY: xios_update_calendar
    USE wxios, ONLY: wxios_closedef, missing_val
#endif
    USE phys_cal_mod, ONLY : mth_len

    IMPLICIT NONE

!   INCLUDE "temps.h"
    INCLUDE "clesphys.h"
    INCLUDE "thermcell.h"
    INCLUDE "compbl.h"
    INCLUDE "YOMCST.h"
    INCLUDE "dimensions.h"
    include "iniprint.h"

    ! Input
    INTEGER :: itap, ivap, read_climoz
    INTEGER, DIMENSION(klon) :: lmax_th
    LOGICAL :: aerosol_couple, ok_sync
    LOGICAL, DIMENSION(klon, klev) :: ptconv, ptconvth
    REAL :: pdtphys
    CHARACTER (LEN=4), DIMENSION(nlevSTD) :: clevSTD
    REAL, DIMENSION(klon,nlevSTD) :: zx_tmp_fi3d_STD
    REAL, DIMENSION(klon) :: pphis
    REAL, DIMENSION(klon, klev) :: pplay, d_t
    REAL, DIMENSION(klon, klev+1) :: paprs
    REAL, DIMENSION(klon,klev,nqtot) :: qx, d_qx
    REAL,DIMENSION(klon,klev,nbtr),INTENT(IN)    :: d_tr_dyn
    REAL, DIMENSION(klon, llm) :: zmasse
    INTEGER :: flag_aerosol_strat
    INTEGER :: flag_aerosol 
    LOGICAL :: ok_cdnc
    REAL, DIMENSION(3) :: freq_moyNMC

    ! Local
    INTEGER, PARAMETER :: jjmp1=jjm+1-1/jjm
    INTEGER :: itau_w
    INTEGER :: i, iinit, iinitend=1, iff, iq, nsrf, k, ll, naero
    REAL, DIMENSION (klon) :: zx_tmp_fi2d
    REAL, DIMENSION (klon,klev) :: zx_tmp_fi3d, zpt_conv
    REAL, DIMENSION (klon,klev+1) :: zx_tmp_fi3d1
    CHARACTER (LEN=4)              :: bb2
    INTEGER, DIMENSION(iim*jjmp1)  :: ndex2d
    INTEGER, DIMENSION(iim*jjmp1*klev) :: ndex3d
    REAL, PARAMETER :: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2
!   REAL, PARAMETER :: missing_val=nf90_fill_real
#ifndef CPP_XIOS
    REAL :: missing_val
#endif
    REAL, PARAMETER :: un_jour=86400.

    ! On calcul le nouveau tau:
    itau_w = itau_phy + itap
    ! On le donne à iophy pour que les histwrite y aient accès:
    CALL set_itau_iophy(itau_w)

    IF (.NOT.vars_defined) THEN
       iinitend = 2
    ELSE
       iinitend = 1
    ENDIF

    ! ug la boucle qui suit ne sert qu'une fois, pour l'initialisation, sinon il n'y a toujours qu'un seul passage:
    DO iinit=1, iinitend
#ifdef CPP_XIOS
       !$OMP MASTER
       IF (vars_defined) THEN
          IF (prt_level >= 10) THEN
             write(lunout,*)"phys_output_write: call xios_update_calendar, itau_w=",itau_w
          ENDIF
!          CALL xios_update_calendar(itau_w)
          CALL xios_update_calendar(itap)
       ENDIF
       !$OMP END MASTER
       !$OMP BARRIER
#endif
       ! On procède à l'écriture ou à la définition des nombreuses variables:
!!! Champs 1D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL histwrite_phy(o_phis, pphis)
       CALL histwrite_phy(o_aire, cell_area)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=pctsrf(i,is_ter)+pctsrf(i,is_lic)
          ENDDO
       ENDIF

       CALL histwrite_phy(o_contfracATM, zx_tmp_fi2d)
       CALL histwrite_phy(o_contfracOR, pctsrf(:,is_ter))
       CALL histwrite_phy(o_aireTER, paire_ter)

!!! Champs 2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! JE20141223 <<
#include "spla_output_write.h"
! JE20141223 >>

       CALL histwrite_phy(o_flat, zxfluxlat)
       CALL histwrite_phy(o_slp, slp)
       CALL histwrite_phy(o_tsol, zxtsol)
       CALL histwrite_phy(o_t2m, zt2m)
       CALL histwrite_phy(o_t2m_min, zt2m)
       CALL histwrite_phy(o_t2m_max, zt2m)
       CALL histwrite_phy(o_t2m_max_mon, t2m_max_mon)
       CALL histwrite_phy(o_t2m_min_mon, t2m_min_mon)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=SQRT(zu10m(i)*zu10m(i)+zv10m(i)*zv10m(i))
          ENDDO
       ENDIF
       CALL histwrite_phy(o_wind10m, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=SQRT(zu10m(i)*zu10m(i)+zv10m(i)*zv10m(i))
          ENDDO
       ENDIF
       CALL histwrite_phy(o_wind10max, zx_tmp_fi2d)

       CALL histwrite_phy(o_gusts, gustiness)

       IF (vars_defined) THEN
          DO i = 1, klon
             zx_tmp_fi2d(i) = pctsrf(i,is_sic)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_sicf, zx_tmp_fi2d)
       CALL histwrite_phy(o_q2m, zq2m)
       CALL histwrite_phy(o_ustar, zustar)
       CALL histwrite_phy(o_u10m, zu10m)
       CALL histwrite_phy(o_v10m, zv10m)

       IF (vars_defined) THEN
          DO i = 1, klon
             zx_tmp_fi2d(i) = paprs(i,1)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_psol, zx_tmp_fi2d)
       CALL histwrite_phy(o_mass, zmasse)
       CALL histwrite_phy(o_qsurf, zxqsurf)

       IF (.NOT. ok_veget) THEN
          CALL histwrite_phy(o_qsol, qsol)
       ENDIF

       IF (vars_defined) THEN
          DO i = 1, klon
             zx_tmp_fi2d(i) = rain_fall(i) + snow_fall(i)
          ENDDO
       ENDIF

       CALL histwrite_phy(o_precip, zx_tmp_fi2d)
       CALL histwrite_phy(o_ndayrain, nday_rain)

       IF (vars_defined) THEN
          DO i = 1, klon
             zx_tmp_fi2d(i) = rain_lsc(i) + snow_lsc(i)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_plul, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i = 1, klon
             zx_tmp_fi2d(i) = rain_con(i) + snow_con(i)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_pluc, zx_tmp_fi2d)
       CALL histwrite_phy(o_snow, snow_fall)
       CALL histwrite_phy(o_msnow, zxsnow)
       CALL histwrite_phy(o_fsnow, zfra_o)
       CALL histwrite_phy(o_evap, evap)
       CALL histwrite_phy(o_tops, topsw)
       CALL histwrite_phy(o_tops0, topsw0)
       CALL histwrite_phy(o_topl, toplw)
       CALL histwrite_phy(o_topl0, toplw0)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1 : klon) = swup ( 1 : klon, klevp1 )
       ENDIF
       CALL histwrite_phy(o_SWupTOA, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1 : klon) = swup0 ( 1 : klon, klevp1 )
       ENDIF
       CALL histwrite_phy(o_SWupTOAclr, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1 : klon) = swdn ( 1 : klon, klevp1 )
       ENDIF
       CALL histwrite_phy(o_SWdnTOA, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1 : klon) = swdn0 ( 1 : klon, klevp1 )
       ENDIF
       CALL histwrite_phy(o_SWdnTOAclr, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = topsw(:)-toplw(:)
       ENDIF
       CALL histwrite_phy(o_nettop, zx_tmp_fi2d)
       CALL histwrite_phy(o_SWup200, SWup200)
       CALL histwrite_phy(o_SWup200clr, SWup200clr)
       CALL histwrite_phy(o_SWdn200, SWdn200)
       CALL histwrite_phy(o_SWdn200clr, SWdn200clr)
       CALL histwrite_phy(o_LWup200, LWup200)
       CALL histwrite_phy(o_LWup200clr, LWup200clr)
       CALL histwrite_phy(o_LWdn200, LWdn200)
       CALL histwrite_phy(o_LWdn200clr, LWdn200clr)
       CALL histwrite_phy(o_sols, solsw)
       CALL histwrite_phy(o_sols0, solsw0)
       CALL histwrite_phy(o_soll, sollw)
       CALL histwrite_phy(o_radsol, radsol)
       CALL histwrite_phy(o_soll0, sollw0)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1 : klon) = swup ( 1 : klon, 1 )
       ENDIF
       CALL histwrite_phy(o_SWupSFC, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1 : klon) = swup0 ( 1 : klon, 1 )
       ENDIF
       CALL histwrite_phy(o_SWupSFCclr, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1 : klon) = swdn ( 1 : klon, 1 )
       ENDIF
       CALL histwrite_phy(o_SWdnSFC, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1 : klon) = swdn0 ( 1 : klon, 1 )
       ENDIF
       CALL histwrite_phy(o_SWdnSFCclr, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1:klon)=sollwdown(1:klon)-sollw(1:klon)
       ENDIF
       CALL histwrite_phy(o_LWupSFC, zx_tmp_fi2d)
       CALL histwrite_phy(o_LWdnSFC, sollwdown)

       IF (vars_defined) THEN
          sollwdownclr(1:klon) = -1.*lwdn0(1:klon,1)
          zx_tmp_fi2d(1:klon)=sollwdownclr(1:klon)-sollw0(1:klon)
       ENDIF
       CALL histwrite_phy(o_LWupSFCclr, zx_tmp_fi2d)
       CALL histwrite_phy(o_LWdnSFCclr, sollwdownclr)
       CALL histwrite_phy(o_bils, bils)
       CALL histwrite_phy(o_bils_diss, bils_diss)
       CALL histwrite_phy(o_bils_ec, bils_ec)
       IF (iflag_ener_conserv>=1) THEN
         CALL histwrite_phy(o_bils_ech, bils_ech)
       ENDIF
       CALL histwrite_phy(o_bils_tke, bils_tke)
       CALL histwrite_phy(o_bils_kinetic, bils_kinetic)
       CALL histwrite_phy(o_bils_latent, bils_latent)
       CALL histwrite_phy(o_bils_enthalp, bils_enthalp)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1:klon)=-1*sens(1:klon)
       ENDIF
       CALL histwrite_phy(o_sens, zx_tmp_fi2d)
       CALL histwrite_phy(o_fder, fder)
       CALL histwrite_phy(o_ffonte, zxffonte)
       CALL histwrite_phy(o_fqcalving, zxfqcalving)
       CALL histwrite_phy(o_fqfonte, zxfqfonte)
       IF (vars_defined) THEN
          zx_tmp_fi2d=0.
          DO nsrf=1,nbsrf
             zx_tmp_fi2d(:)=zx_tmp_fi2d(:)+pctsrf(:,nsrf)*fluxu(:,1,nsrf)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_taux, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d=0.
          DO nsrf=1,nbsrf
             zx_tmp_fi2d(:)=zx_tmp_fi2d(:)+pctsrf(:,nsrf)*fluxv(:,1,nsrf)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_tauy, zx_tmp_fi2d)

       IF (ok_snow) THEN
          CALL histwrite_phy(o_snowsrf, snow_o)
          CALL histwrite_phy(o_qsnow, qsnow)
          CALL histwrite_phy(o_snowhgt,snowhgt)
          CALL histwrite_phy(o_toice,to_ice)
          CALL histwrite_phy(o_sissnow,sissnow)
          CALL histwrite_phy(o_runoff,runoff)
          CALL histwrite_phy(o_albslw3,albsol3_lic)
       ENDIF

       DO nsrf = 1, nbsrf
          IF (vars_defined)             zx_tmp_fi2d(1 : klon) = pctsrf( 1 : klon, nsrf)*100.
          CALL histwrite_phy(o_pourc_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)           zx_tmp_fi2d(1 : klon) = pctsrf( 1 : klon, nsrf)
          CALL histwrite_phy(o_fract_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = fluxu( 1 : klon, 1, nsrf)
          CALL histwrite_phy(o_taux_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = fluxv( 1 : klon, 1, nsrf)
          CALL histwrite_phy(o_tauy_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = ftsol( 1 : klon, nsrf)
          CALL histwrite_phy(o_tsol_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = evap_pot( 1 : klon, nsrf)
          CALL histwrite_phy(o_evappot_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)       zx_tmp_fi2d(1 : klon) = ustar(1 : klon, nsrf)
          CALL histwrite_phy(o_ustar_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)       zx_tmp_fi2d(1 : klon) = u10m(1 : klon, nsrf)
          CALL histwrite_phy(o_u10m_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)       zx_tmp_fi2d(1 : klon) = v10m(1 : klon, nsrf)
          CALL histwrite_phy(o_v10m_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)       zx_tmp_fi2d(1 : klon) = t2m(1 : klon, nsrf)
          CALL histwrite_phy(o_t2m_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)       zx_tmp_fi2d(1 : klon) = fevap(1 : klon, nsrf)
          CALL histwrite_phy(o_evap_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)        zx_tmp_fi2d(1 : klon) = fluxt( 1 : klon, 1, nsrf)
          CALL histwrite_phy(o_sens_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = fluxlat( 1 : klon, nsrf)
          CALL histwrite_phy(o_lat_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = fsollw( 1 : klon, nsrf)
          CALL histwrite_phy(o_flw_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = fsolsw( 1 : klon, nsrf)
          CALL histwrite_phy(o_fsw_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = wfbils( 1 : klon, nsrf)
          CALL histwrite_phy(o_wbils_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = wfbilo( 1 : klon, nsrf)
          CALL histwrite_phy(o_wbilo_srf(nsrf), zx_tmp_fi2d)

          IF (iflag_pbl > 1) THEN
             CALL histwrite_phy(o_tke_srf(nsrf),  pbl_tke(:,1:klev,nsrf))
             CALL histwrite_phy(o_tke_max_srf(nsrf),  pbl_tke(:,1:klev,nsrf))
          ENDIF
!jyg<
          IF (iflag_pbl > 1 .AND. iflag_wake>=1  .AND. iflag_pbl_split >=1) THEN
             CALL histwrite_phy(o_dltpbltke_srf(nsrf), wake_delta_pbl_TKE(:,1:klev,nsrf))
          ENDIF
!>jyg

       ENDDO
       DO nsrf=1,nbsrf+1
          CALL histwrite_phy(o_wstar(nsrf), wstar(1 : klon, nsrf))
       ENDDO

       CALL histwrite_phy(o_cdrm, cdragm)
       CALL histwrite_phy(o_cdrh, cdragh)
       CALL histwrite_phy(o_cldl, cldl)
       CALL histwrite_phy(o_cldm, cldm)
       CALL histwrite_phy(o_cldh, cldh)
       CALL histwrite_phy(o_cldt, cldt)
       CALL histwrite_phy(o_JrNt, JrNt)
       CALL histwrite_phy(o_cldljn, cldl*JrNt)
       CALL histwrite_phy(o_cldmjn, cldm*JrNt)
       CALL histwrite_phy(o_cldhjn, cldh*JrNt)
       CALL histwrite_phy(o_cldtjn, cldt*JrNt)
       CALL histwrite_phy(o_cldq, cldq)
       IF (vars_defined)       zx_tmp_fi2d(1:klon) = flwp(1:klon)
       CALL histwrite_phy(o_lwp, zx_tmp_fi2d)
       IF (vars_defined)       zx_tmp_fi2d(1:klon) = fiwp(1:klon)
       CALL histwrite_phy(o_iwp, zx_tmp_fi2d)
       CALL histwrite_phy(o_ue, ue)
       CALL histwrite_phy(o_ve, ve)
       CALL histwrite_phy(o_uq, uq)
       CALL histwrite_phy(o_vq, vq)
       IF (iflag_con.GE.3) THEN ! sb
          CALL histwrite_phy(o_cape, cape)
          CALL histwrite_phy(o_pbase, ema_pcb)
          CALL histwrite_phy(o_ptop, ema_pct)
          CALL histwrite_phy(o_fbase, ema_cbmf)
          IF (iflag_con /= 30) THEN
             CALL histwrite_phy(o_plcl, plcl)
             CALL histwrite_phy(o_plfc, plfc)
             CALL histwrite_phy(o_wbeff, wbeff)
          ENDIF

          CALL histwrite_phy(o_cape_max, cape)

          CALL histwrite_phy(o_upwd, upwd)
          CALL histwrite_phy(o_Ma, Ma)
          CALL histwrite_phy(o_dnwd, dnwd)
          CALL histwrite_phy(o_dnwd0, dnwd0)
          IF (vars_defined)         zx_tmp_fi2d=float(itau_con)/float(itap)
          CALL histwrite_phy(o_ftime_con, zx_tmp_fi2d)
          IF (vars_defined) THEN
             IF (iflag_thermals>=1)THEN
                zx_tmp_fi3d=dnwd+dnwd0+upwd+fm_therm(:,1:klev)
             ELSE
                zx_tmp_fi3d=dnwd+dnwd0+upwd
             ENDIF
          ENDIF
          CALL histwrite_phy(o_mc, zx_tmp_fi3d)
       ENDIF !iflag_con .GE. 3
       CALL histwrite_phy(o_prw, prw)
       CALL histwrite_phy(o_s_pblh, s_pblh)
       CALL histwrite_phy(o_s_pblt, s_pblt)
       CALL histwrite_phy(o_s_lcl, s_lcl)
       CALL histwrite_phy(o_s_therm, s_therm)
       !IM : Les champs suivants (s_capCL, s_oliqCL, s_cteiCL, s_trmb1, s_trmb2, s_trmb3) ne sont pas definis dans HBTM.F
       !       IF (o_s_capCL%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_capCL%name,itau_w,s_capCL)
       !       ENDIF
       !       IF (o_s_oliqCL%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_oliqCL%name,itau_w,s_oliqCL)
       !       ENDIF
       !       IF (o_s_cteiCL%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_cteiCL%name,itau_w,s_cteiCL)
       !       ENDIF
       !       IF (o_s_trmb1%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_trmb1%name,itau_w,s_trmb1)
       !       ENDIF
       !       IF (o_s_trmb2%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_trmb2%name,itau_w,s_trmb2)
       !       ENDIF
       !       IF (o_s_trmb3%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_trmb3%name,itau_w,s_trmb3)
       !       ENDIF

#ifdef CPP_IOIPSL
#ifndef CPP_XIOS
  IF (.NOT.ok_all_xml) THEN
       ! ATTENTION, LES ANCIENS HISTWRITE ONT ETES CONSERVES EN ATTENDANT MIEUX:
       ! Champs interpolles sur des niveaux de pression
       missing_val=missing_val_nf90
       DO iff=1, nfiles
          ll=0
          DO k=1, nlevSTD
             bb2=clevSTD(k) 
             IF (bb2.EQ."850".OR.bb2.EQ."700".OR. &
                  bb2.EQ."500".OR.bb2.EQ."200".OR. &
                  bb2.EQ."100".OR. &
                  bb2.EQ."50".OR.bb2.EQ."10") THEN

                ! a refaire correctement !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ll=ll+1
                CALL histwrite_phy(o_uSTDlevs(ll),uwriteSTD(:,k,iff), iff)
                CALL histwrite_phy(o_vSTDlevs(ll),vwriteSTD(:,k,iff), iff)
                CALL histwrite_phy(o_wSTDlevs(ll),wwriteSTD(:,k,iff), iff)
                CALL histwrite_phy(o_zSTDlevs(ll),phiwriteSTD(:,k,iff), iff)
                CALL histwrite_phy(o_qSTDlevs(ll),qwriteSTD(:,k,iff), iff)
                CALL histwrite_phy(o_tSTDlevs(ll),twriteSTD(:,k,iff), iff)

             ENDIF !(bb2.EQ."850".OR.bb2.EQ."700".OR.
          ENDDO
       ENDDO
  ENDIF
#endif
#endif

#ifdef CPP_XIOS
  IF (ok_all_xml) THEN
!XIOS  CALL xios_get_field_attr("u850",default_value=missing_val)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ll=0
          DO k=1, nlevSTD
             bb2=clevSTD(k) 
             IF (bb2.EQ."850".OR.bb2.EQ."700".OR. &
                bb2.EQ."500".OR.bb2.EQ."200".OR. &
                bb2.EQ."100".OR. &
                bb2.EQ."50".OR.bb2.EQ."10") THEN
                ll=ll+1
                CALL histwrite_phy(o_uSTDlevs(ll),ulevSTD(:,k))
                CALL histwrite_phy(o_vSTDlevs(ll),vlevSTD(:,k))
                CALL histwrite_phy(o_wSTDlevs(ll),wlevSTD(:,k))
                CALL histwrite_phy(o_zSTDlevs(ll),philevSTD(:,k))
                CALL histwrite_phy(o_qSTDlevs(ll),qlevSTD(:,k))
                CALL histwrite_phy(o_tSTDlevs(ll),tlevSTD(:,k))
             ENDIF !(bb2.EQ."850".OR.bb2.EQ."700".OR.
          ENDDO
  ENDIF
#endif
       IF (vars_defined) THEN
          DO i=1, klon
             IF (pctsrf(i,is_oce).GT.epsfra.OR. &
                  pctsrf(i,is_sic).GT.epsfra) THEN
                zx_tmp_fi2d(i) = (ftsol(i, is_oce) * pctsrf(i,is_oce)+ &
                     ftsol(i, is_sic) * pctsrf(i,is_sic))/ &
                     (pctsrf(i,is_oce)+pctsrf(i,is_sic))
             ELSE
                zx_tmp_fi2d(i) = 273.15
             ENDIF
          ENDDO
       ENDIF
       CALL histwrite_phy(o_t_oce_sic, zx_tmp_fi2d)

       ! Couplage convection-couche limite
       IF (iflag_con.GE.3) THEN
          IF (iflag_coupl>=1) THEN
             CALL histwrite_phy(o_ale_bl, ale_bl)
             CALL histwrite_phy(o_alp_bl, alp_bl)
          ENDIF !iflag_coupl>=1
       ENDIF !(iflag_con.GE.3)
       ! Wakes
       IF (iflag_con.EQ.3) THEN
          IF (iflag_wake>=1) THEN
             CALL histwrite_phy(o_ale_wk, ale_wake)
             CALL histwrite_phy(o_alp_wk, alp_wake)
             CALL histwrite_phy(o_ale, ale)
             CALL histwrite_phy(o_alp, alp)
             CALL histwrite_phy(o_cin, cin)
             CALL histwrite_phy(o_WAPE, wake_pe)
             CALL histwrite_phy(o_wake_h, wake_h)
             CALL histwrite_phy(o_wake_s, wake_s)
             CALL histwrite_phy(o_wake_deltat, wake_deltat)
             CALL histwrite_phy(o_wake_deltaq, wake_deltaq)
             CALL histwrite_phy(o_wake_omg, wake_omg)
             IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_wake(1:klon,1:klev) &
                  /pdtphys
             CALL histwrite_phy(o_dtwak, zx_tmp_fi3d)
             IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_wake(1:klon,1:klev)/pdtphys
             CALL histwrite_phy(o_dqwak, zx_tmp_fi3d)
          ENDIF ! iflag_wake>=1
          CALL histwrite_phy(o_ftd, ftd)
          CALL histwrite_phy(o_fqd, fqd)
       ENDIF !(iflag_con.EQ.3)
       IF (iflag_con.EQ.3.OR.iflag_con.EQ.30) THEN
          ! sortie RomP convection descente insaturee iflag_con=30
          ! etendue a iflag_con=3 (jyg)
          CALL histwrite_phy(o_Vprecip, Vprecip)
          CALL histwrite_phy(o_wdtrainA, wdtrainA)
          CALL histwrite_phy(o_wdtrainM, wdtrainM)
       ENDIF !(iflag_con.EQ.3.or.iflag_con.EQ.30)
!!! nrlmd le 10/04/2012
       IF (iflag_trig_bl>=1) THEN
          CALL histwrite_phy(o_n2, n2)
          CALL histwrite_phy(o_s2, s2)
          CALL histwrite_phy(o_proba_notrig, proba_notrig)
          CALL histwrite_phy(o_random_notrig, random_notrig)
          CALL histwrite_phy(o_ale_bl_stat, ale_bl_stat)
          CALL histwrite_phy(o_ale_bl_trig, ale_bl_trig)
       ENDIF  !(iflag_trig_bl>=1)
       IF (iflag_clos_bl>=1) THEN
          CALL histwrite_phy(o_alp_bl_det, alp_bl_det)
          CALL histwrite_phy(o_alp_bl_fluct_m, alp_bl_fluct_m)
          CALL histwrite_phy(o_alp_bl_fluct_tke,  &
               alp_bl_fluct_tke)
          CALL histwrite_phy(o_alp_bl_conv, alp_bl_conv)
          CALL histwrite_phy(o_alp_bl_stat, alp_bl_stat)
       ENDIF  !(iflag_clos_bl>=1)
!!! fin nrlmd le 10/04/2012
       ! Output of slab ocean variables
       IF (type_ocean=='slab ') THEN
          CALL histwrite_phy(o_slab_qflux, slab_wfbils)
          CALL histwrite_phy(o_slab_bils, slab_bils)
          IF (nslay.EQ.1) THEN
              zx_tmp_fi2d(:)=tslab(:,1)
              CALL histwrite_phy(o_tslab, zx_tmp_fi2d)
          ELSE
              CALL histwrite_phy(o_tslab, tslab)
          ENDIF
          IF (version_ocean=='sicINT') THEN
              CALL histwrite_phy(o_slab_bilg, slab_bilg)
              CALL histwrite_phy(o_slab_tice, tice)
              CALL histwrite_phy(o_slab_sic, seaice)
          ENDIF
       ENDIF !type_ocean == force/slab
       CALL histwrite_phy(o_weakinv, weak_inversion)
       CALL histwrite_phy(o_dthmin, dthmin)
       CALL histwrite_phy(o_cldtau, cldtau)
       CALL histwrite_phy(o_cldemi, cldemi)
       CALL histwrite_phy(o_pr_con_l, pmflxr(:,1:klev))
       CALL histwrite_phy(o_pr_con_i, pmflxs(:,1:klev))
       CALL histwrite_phy(o_pr_lsc_l, prfl(:,1:klev))
       CALL histwrite_phy(o_pr_lsc_i, psfl(:,1:klev))
       CALL histwrite_phy(o_re, re)
       CALL histwrite_phy(o_fl, fl)
       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=MIN(100.,rh2m(i)*100.)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_rh2m, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=MIN(100.,rh2m(i)*100.)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_rh2m_min, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=MIN(100.,rh2m(i)*100.)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_rh2m_max, zx_tmp_fi2d)

       CALL histwrite_phy(o_qsat2m, qsat2m)
       CALL histwrite_phy(o_tpot, tpot)
       CALL histwrite_phy(o_tpote, tpote)
       IF (vars_defined) zx_tmp_fi2d(1 : klon) = fsolsw( 1 : klon, is_ter)
       CALL histwrite_phy(o_SWnetOR,  zx_tmp_fi2d)
       CALL histwrite_phy(o_LWdownOR, sollwdown)
       CALL histwrite_phy(o_snowl, snow_lsc)
       CALL histwrite_phy(o_solldown, sollwdown)
       CALL histwrite_phy(o_dtsvdfo, d_ts(:,is_oce))
       CALL histwrite_phy(o_dtsvdft, d_ts(:,is_ter))
       CALL histwrite_phy(o_dtsvdfg,  d_ts(:,is_lic))
       CALL histwrite_phy(o_dtsvdfi, d_ts(:,is_sic))
       CALL histwrite_phy(o_z0m, z0m(:,nbsrf+1))
       CALL histwrite_phy(o_z0h, z0h(:,nbsrf+1))
       ! OD550 per species
!--OLIVIER
!This is warranted by treating INCA aerosols as offline aerosols
       IF (flag_aerosol.GT.0) THEN
          CALL histwrite_phy(o_od550aer, od550aer)
          CALL histwrite_phy(o_od865aer, od865aer)
          CALL histwrite_phy(o_absvisaer, absvisaer)
          CALL histwrite_phy(o_od550lt1aer, od550lt1aer)
          CALL histwrite_phy(o_sconcso4, sconcso4)
          CALL histwrite_phy(o_sconcno3, sconcno3)
          CALL histwrite_phy(o_sconcoa, sconcoa)
          CALL histwrite_phy(o_sconcbc, sconcbc)
          CALL histwrite_phy(o_sconcss, sconcss)
          CALL histwrite_phy(o_sconcdust, sconcdust)
          CALL histwrite_phy(o_concso4, concso4)
          CALL histwrite_phy(o_concno3, concno3)
          CALL histwrite_phy(o_concoa, concoa)
          CALL histwrite_phy(o_concbc, concbc)
          CALL histwrite_phy(o_concss, concss)
          CALL histwrite_phy(o_concdust, concdust)
          CALL histwrite_phy(o_loadso4, loadso4)
          CALL histwrite_phy(o_loadoa, loadoa)
          CALL histwrite_phy(o_loadbc, loadbc)
          CALL histwrite_phy(o_loadss, loadss)
          CALL histwrite_phy(o_loaddust, loaddust)
          !--STRAT AER
       ENDIF
       IF (flag_aerosol.GT.0.OR.flag_aerosol_strat>=1) THEN
!          DO naero = 1, naero_spc
!--correction mini bug OB
          DO naero = 1, naero_tot
             CALL histwrite_phy(o_tausumaero(naero), &
                  tausum_aero(:,2,naero) )
          ENDDO
       ENDIF
       IF (flag_aerosol_strat>=1) THEN
          CALL histwrite_phy(o_tausumaero_lw, &
               tausum_aero(:,6,id_STRAT_phy) )
       ENDIF
       IF (ok_ade) THEN
          CALL histwrite_phy(o_topswad, topswad_aero)
          CALL histwrite_phy(o_topswad0, topswad0_aero)
          CALL histwrite_phy(o_solswad, solswad_aero)
          CALL histwrite_phy(o_solswad0, solswad0_aero)
          CALL histwrite_phy(o_toplwad, toplwad_aero)
          CALL histwrite_phy(o_toplwad0, toplwad0_aero)
          CALL histwrite_phy(o_sollwad, sollwad_aero)
          CALL histwrite_phy(o_sollwad0, sollwad0_aero)
          !====MS forcing diagnostics
          CALL histwrite_phy(o_swtoaas_nat, topsw_aero(:,1))
          CALL histwrite_phy(o_swsrfas_nat, solsw_aero(:,1))
          CALL histwrite_phy(o_swtoacs_nat, topsw0_aero(:,1))
          CALL histwrite_phy(o_swsrfcs_nat, solsw0_aero(:,1))
          !ant
          CALL histwrite_phy(o_swtoaas_ant, topsw_aero(:,2))
          CALL histwrite_phy(o_swsrfas_ant, solsw_aero(:,2))
          CALL histwrite_phy(o_swtoacs_ant, topsw0_aero(:,2))
          CALL histwrite_phy(o_swsrfcs_ant, solsw0_aero(:,2))
          !cf
          IF (.not. aerosol_couple) THEN
             CALL histwrite_phy(o_swtoacf_nat, topswcf_aero(:,1))
             CALL histwrite_phy(o_swsrfcf_nat, solswcf_aero(:,1))
             CALL histwrite_phy(o_swtoacf_ant, topswcf_aero(:,2))
             CALL histwrite_phy(o_swsrfcf_ant, solswcf_aero(:,2))
             CALL histwrite_phy(o_swtoacf_zero,topswcf_aero(:,3))
             CALL histwrite_phy(o_swsrfcf_zero,solswcf_aero(:,3))
          ENDIF
          !====MS forcing diagnostics
       ENDIF
       IF (ok_aie) THEN
          CALL histwrite_phy(o_topswai, topswai_aero)
          CALL histwrite_phy(o_solswai, solswai_aero)
       ENDIF
       IF (flag_aerosol.GT.0.AND.ok_cdnc) THEN
          CALL histwrite_phy(o_scdnc, scdnc)
          CALL histwrite_phy(o_cldncl, cldncl)
          CALL histwrite_phy(o_reffclws, reffclws)
          CALL histwrite_phy(o_reffclwc, reffclwc)
          CALL histwrite_phy(o_cldnvi, cldnvi)
          CALL histwrite_phy(o_lcc, lcc)
          CALL histwrite_phy(o_lcc3d, lcc3d)
          CALL histwrite_phy(o_lcc3dcon, lcc3dcon)
          CALL histwrite_phy(o_lcc3dstra, lcc3dstra)
          CALL histwrite_phy(o_reffclwtop, reffclwtop)
       ENDIF
       ! Champs 3D:
       IF (ok_ade .OR. ok_aie) THEN
          CALL histwrite_phy(o_ec550aer, ec550aer)
       ENDIF
       CALL histwrite_phy(o_lwcon, flwc)
       CALL histwrite_phy(o_iwcon, fiwc)
       CALL histwrite_phy(o_temp, t_seri)
       CALL histwrite_phy(o_theta, theta)
       CALL histwrite_phy(o_ovapinit, qx(:,:,ivap))
       CALL histwrite_phy(o_ovap, q_seri)
       CALL histwrite_phy(o_oliq, ql_seri)
       CALL histwrite_phy(o_geop, zphi)
       CALL histwrite_phy(o_vitu, u_seri)
       CALL histwrite_phy(o_vitv, v_seri)
       CALL histwrite_phy(o_vitw, omega)
       CALL histwrite_phy(o_pres, pplay)
       CALL histwrite_phy(o_paprs, paprs(:,1:klev))
       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi3d1(i,1)= pphis(i)/RG
             !020611   zx_tmp_fi3d(i,1)= pphis(i)/RG
          ENDDO
          DO k=1, klev
             !020611        DO k=1, klev-1
             DO i=1, klon
                !020611         zx_tmp_fi3d(i,k+1)= zx_tmp_fi3d(i,k) - (t_seri(i,k) *RD * 
                zx_tmp_fi3d1(i,k+1)= zx_tmp_fi3d1(i,k) - (t_seri(i,k) *RD *  &
                     (paprs(i,k+1) - paprs(i,k))) / ( pplay(i,k) * RG ) 
             ENDDO
          ENDDO
       ENDIF
       CALL histwrite_phy(o_zfull,zx_tmp_fi3d1(:,2:klevp1))
       !020611    $o_zfull%name,itau_w,zx_tmp_fi3d)

       IF (vars_defined)  THEN
          DO i=1, klon
             zx_tmp_fi3d(i,1)= pphis(i)/RG - ( &
                  (t_seri(i,1)+zxtsol(i))/2. *RD * &
                  (pplay(i,1) - paprs(i,1)))/( (paprs(i,1)+pplay(i,1))/2.* RG)
          ENDDO
          DO k=1, klev-1
             DO i=1, klon
                zx_tmp_fi3d(i,k+1)= zx_tmp_fi3d(i,k) - ( &
                     (t_seri(i,k)+t_seri(i,k+1))/2. *RD *  &
                     (pplay(i,k+1) - pplay(i,k))) / ( paprs(i,k) * RG ) 
             ENDDO
          ENDDO
       ENDIF
       CALL histwrite_phy(o_zhalf, zx_tmp_fi3d)
       CALL histwrite_phy(o_rneb, cldfra)
       CALL histwrite_phy(o_rnebcon, rnebcon)
       CALL histwrite_phy(o_rnebls, rneb)
       IF (vars_defined)  THEN
          DO k=1, klev
             DO i=1, klon
                zx_tmp_fi3d(i,k)=cldfra(i,k)*JrNt(i)
             ENDDO
          ENDDO
       ENDIF
       CALL histwrite_phy(o_rnebjn, zx_tmp_fi3d)
       CALL histwrite_phy(o_rhum, zx_rh)
       CALL histwrite_phy(o_ozone, &
            wo(:, :, 1) * dobson_u * 1e3 / zmasse / rmo3 * rmd)

       IF (read_climoz == 2) THEN
          CALL histwrite_phy(o_ozone_light, &
               wo(:, :, 2) * dobson_u * 1e3 / zmasse / rmo3 * rmd)
       ENDIF

       CALL histwrite_phy(o_dtphy, d_t)
       CALL histwrite_phy(o_dqphy,  d_qx(:,:,ivap))
       DO nsrf=1, nbsrf
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = falb1( 1 : klon, nsrf)
          CALL histwrite_phy(o_albe_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = z0m( 1 : klon, nsrf)
          CALL histwrite_phy(o_z0m_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = z0h( 1 : klon, nsrf)
          CALL histwrite_phy(o_z0h_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = agesno( 1 : klon, nsrf)
          CALL histwrite_phy(o_ages_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = snow( 1 : klon, nsrf)
          CALL histwrite_phy(o_snow_srf(nsrf), zx_tmp_fi2d)
       ENDDO !nsrf=1, nbsrf
       CALL histwrite_phy(o_alb1, albsol1)
       CALL histwrite_phy(o_alb2, albsol2)
       !FH Sorties pour la couche limite
       IF (iflag_pbl>1) THEN
          zx_tmp_fi3d=0.
          IF (vars_defined) THEN
             DO nsrf=1,nbsrf
                DO k=1,klev
                   zx_tmp_fi3d(:,k)=zx_tmp_fi3d(:,k) &
                        +pctsrf(:,nsrf)*pbl_tke(:,k,nsrf)
                enddo
             enddo
          ENDIF
          CALL histwrite_phy(o_tke, zx_tmp_fi3d)

          CALL histwrite_phy(o_tke_max, zx_tmp_fi3d)
       ENDIF

       CALL histwrite_phy(o_kz, coefh(:,:,is_ave))

       CALL histwrite_phy(o_kz_max, coefh(:,:,is_ave))

       CALL histwrite_phy(o_clwcon, clwcon0)
       CALL histwrite_phy(o_dtdyn, d_t_dyn)
       CALL histwrite_phy(o_dqdyn, d_q_dyn)
       CALL histwrite_phy(o_dudyn, d_u_dyn)
       CALL histwrite_phy(o_dvdyn, d_v_dyn)

       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_t_con(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_dtcon, zx_tmp_fi3d)
       IF (iflag_thermals.eq.0)THEN
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=d_t_con(1:klon,1:klev)/pdtphys + &
                  d_t_ajsb(1:klon,1:klev)/pdtphys
          ENDIF
          CALL histwrite_phy(o_tntc, zx_tmp_fi3d)
       ELSEIF (iflag_thermals.ge.1.and.iflag_wake.EQ.1)THEN
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=d_t_con(1:klon,1:klev)/pdtphys + &
                  d_t_ajs(1:klon,1:klev)/pdtphys + &
                  d_t_wake(1:klon,1:klev)/pdtphys
          ENDIF
          CALL histwrite_phy(o_tntc, zx_tmp_fi3d)
       ENDIF
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_u_con(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_ducon, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_v_con(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dvcon, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_con(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dqcon, zx_tmp_fi3d)

       IF (iflag_thermals.EQ.0) THEN
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_con(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_tnhusc, zx_tmp_fi3d)
       ELSEIF (iflag_thermals.GE.1.AND.iflag_wake.EQ.1) THEN
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=d_q_con(1:klon,1:klev)/pdtphys + &
                  d_q_ajs(1:klon,1:klev)/pdtphys + &
                  d_q_wake(1:klon,1:klev)/pdtphys
          ENDIF
          CALL histwrite_phy(o_tnhusc, zx_tmp_fi3d)
       ENDIF

       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lsc(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtlsc, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon, 1:klev)=(d_t_lsc(1:klon,1:klev)+ &
            d_t_eva(1:klon,1:klev))/pdtphys
       CALL histwrite_phy(o_dtlschr, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_lsc(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dqlsc, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=beta_prec(1:klon,1:klev)
       CALL histwrite_phy(o_beta_prec, zx_tmp_fi3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Sorties specifiques a la separation thermiques/non thermiques
       IF (iflag_thermals>=1) THEN
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lscth(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dtlscth, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lscst(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dtlscst, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_lscth(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dqlscth, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_lscst(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dqlscst, zx_tmp_fi3d)
          CALL histwrite_phy(o_plulth, plul_th)
          CALL histwrite_phy(o_plulst, plul_st)
          IF (vars_defined) THEN
             DO k=1,klev
                DO i=1,klon
                   IF (ptconvth(i,k)) THEN
                      zx_tmp_fi3d(i,k)=1.
                   ELSE
                      zx_tmp_fi3d(i,k)=0.
                   ENDIF
                enddo
             enddo
          ENDIF
          CALL histwrite_phy(o_ptconvth, zx_tmp_fi3d)
          IF (vars_defined) THEN
             DO i=1,klon
                zx_tmp_fi2d(1:klon)=lmax_th(:)
             enddo
          ENDIF
          CALL histwrite_phy(o_lmaxth, zx_tmp_fi2d)
       ENDIF ! iflag_thermals>=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_vdf(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtvdf, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_diss(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtdis, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_vdf(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dqvdf, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_eva(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dteva, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_eva(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dqeva, zx_tmp_fi3d)
       zpt_conv = 0.
       WHERE (ptconv) zpt_conv = 1.
       CALL histwrite_phy(o_ptconv, zpt_conv)
       CALL histwrite_phy(o_ratqs, ratqs)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_t_ajs(1:klon,1:klev)/pdtphys - &
               d_t_ajsb(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_dtthe, zx_tmp_fi3d)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_u_ajs(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_duthe, zx_tmp_fi3d) 
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_v_ajs(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_dvthe, zx_tmp_fi3d)

       IF (iflag_thermals>=1) THEN
          ! Pour l instant 0 a y reflichir pour les thermiques
          zx_tmp_fi2d=0. 
          CALL histwrite_phy(o_ftime_th, zx_tmp_fi2d)
          CALL histwrite_phy(o_f_th, fm_therm)
          CALL histwrite_phy(o_e_th, entr_therm)
          CALL histwrite_phy(o_w_th, zw2)
          CALL histwrite_phy(o_q_th, zqasc)
          CALL histwrite_phy(o_a_th, fraca)
          CALL histwrite_phy(o_d_th, detr_therm)
          CALL histwrite_phy(o_f0_th, f0)
          CALL histwrite_phy(o_zmax_th, zmax_th)
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=d_q_ajs(1:klon,1:klev)/pdtphys - &
                  d_q_ajsb(1:klon,1:klev)/pdtphys
          ENDIF
          CALL histwrite_phy(o_dqthe, zx_tmp_fi3d)
       ENDIF !iflag_thermals
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_ajsb(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtajs, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_ajsb(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dqajs, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_swr(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtswr, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_sw0(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtsw0, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lwr(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtlwr, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lw0(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtlw0, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_ec(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtec, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_u_vdf(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_duvdf, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_v_vdf(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dvvdf, zx_tmp_fi3d)
       IF (ok_orodr) THEN
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_u_oro(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_duoro, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_v_oro(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dvoro, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_oro(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dtoro, zx_tmp_fi3d)
       ENDIF
       IF (ok_orolf) THEN
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_u_lIF (1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dulif, zx_tmp_fi3d)

          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_v_lIF (1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dvlif, zx_tmp_fi3d)

          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lIF (1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dtlif, zx_tmp_fi3d)
       ENDIF

!      IF (ok_hines) THEN
!         IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_u_hin(1:klon,1:klev)/pdtphys
!         CALL histwrite_phy(o_duhin, zx_tmp_fi3d)
!         IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_v_hin(1:klon,1:klev)/pdtphys
!         CALL histwrite_phy(o_dvhin, zx_tmp_fi3d)
!         IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_hin(1:klon,1:klev)/pdtphys
!         CALL histwrite_phy(o_dthin, zx_tmp_fi3d)
!      ENDIF

!      IF (ok_gwd_rando) THEN
!         CALL histwrite_phy(o_du_gwd_rando, du_gwd_ranDO / pdtphys)
!         CALL histwrite_phy(o_dv_gwd_rando, dv_gwd_ranDO / pdtphys)
!         CALL histwrite_phy(o_vstr_gwd_rando, zvstr_gwd_rando)
!      ENDIF

       IF (ok_qch4) THEN
          CALL histwrite_phy(o_dqch4, d_q_ch4 / pdtphys)
       ENDIF

       CALL histwrite_phy(o_rsu, swup)
       CALL histwrite_phy(o_rsd, swdn)
       CALL histwrite_phy(o_rlu, lwup)
       CALL histwrite_phy(o_rld, lwdn)
       CALL histwrite_phy(o_rsucs, swup0)
       CALL histwrite_phy(o_rsdcs, swdn0)
       CALL histwrite_phy(o_rlucs, lwup0)
       CALL histwrite_phy(o_rldcs, lwdn0)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_t(1:klon,1:klev)+ &
               d_t_dyn(1:klon,1:klev)
       ENDIF
       CALL histwrite_phy(o_tnt, zx_tmp_fi3d)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_t_swr(1:klon,1:klev)/pdtphys + &
               d_t_lwr(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_tntr, zx_tmp_fi3d)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)= (d_t_lsc(1:klon,1:klev)+ &
               d_t_eva(1:klon,1:klev)+ &
               d_t_vdf(1:klon,1:klev))/pdtphys
       ENDIF
       CALL histwrite_phy(o_tntscpbl, zx_tmp_fi3d)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_qx(1:klon,1:klev,ivap)+ &
               d_q_dyn(1:klon,1:klev)
       ENDIF
       CALL histwrite_phy(o_tnhus, zx_tmp_fi3d)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_q_lsc(1:klon,1:klev)/pdtphys+ &
               d_q_eva(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_tnhusscpbl, zx_tmp_fi3d)
       CALL histwrite_phy(o_evu, coefm(:,:,is_ave))
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=q_seri(1:klon,1:klev)+ &
               ql_seri(1:klon,1:klev) 
       ENDIF
       CALL histwrite_phy(o_h2o, zx_tmp_fi3d)
       IF (iflag_con >= 3) THEN
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=-1 * (dnwd(1:klon,1:klev)+ &
                  dnwd0(1:klon,1:klev)) 
          ENDIF
          CALL histwrite_phy(o_mcd, zx_tmp_fi3d)
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=upwd(1:klon,1:klev) + &
                  dnwd(1:klon,1:klev)+ dnwd0(1:klon,1:klev) 
          ENDIF
          CALL histwrite_phy(o_dmc, zx_tmp_fi3d)
       ELSEIF (iflag_con == 2) THEN
          CALL histwrite_phy(o_mcd,  pmfd)
          CALL histwrite_phy(o_dmc,  pmfu + pmfd)
       ENDIF
       CALL histwrite_phy(o_ref_liq, ref_liq)
       CALL histwrite_phy(o_ref_ice, ref_ice)
       IF (RCO2_per.NE.RCO2_act.OR.RCH4_per.NE.RCH4_act.OR. &
            RN2O_per.NE.RN2O_act.OR.RCFC11_per.NE.RCFC11_act.OR. &
            RCFC12_per.NE.RCFC12_act) THEN
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = swupp ( 1 : klon, klevp1 )
          CALL histwrite_phy(o_rsut4co2, zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = lwupp ( 1 : klon, klevp1 )
          CALL histwrite_phy(o_rlut4co2, zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = swup0p ( 1 : klon, klevp1 )
          CALL histwrite_phy(o_rsutcs4co2, zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = lwup0p ( 1 : klon, klevp1 )
          CALL histwrite_phy(o_rlutcs4co2, zx_tmp_fi2d)
          CALL histwrite_phy(o_rsu4co2, swupp)
          CALL histwrite_phy(o_rlu4co2, lwupp)
          CALL histwrite_phy(o_rsucs4co2, swup0p)
          CALL histwrite_phy(o_rlucs4co2, lwup0p)
          CALL histwrite_phy(o_rsd4co2, swdnp)
          CALL histwrite_phy(o_rld4co2, lwdnp)
          CALL histwrite_phy(o_rsdcs4co2, swdn0p)
          CALL histwrite_phy(o_rldcs4co2, lwdn0p)
       ENDIF
!!!!!!!!!!!! Sorties niveaux de pression NMC !!!!!!!!!!!!!!!!!!!!
#ifdef CPP_IOIPSL
#ifndef CPP_XIOS
  IF (.NOT.ok_all_xml) THEN 
       ! ATTENTION, LES ANCIENS HISTWRITE ONT ETES CONSERVES EN ATTENDANT MIEUX:
       ! Champs interpolles sur des niveaux de pression
       missing_val=missing_val_nf90
       DO iff=7, nfiles-1 !--here we deal with files 7,8 and 9

          CALL histwrite_phy(o_tnondef,tnondef(:,:,iff-6),iff)
          CALL histwrite_phy(o_ta,twriteSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_zg,phiwriteSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_hus,qwriteSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_hur,rhwriteSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_ua,uwriteSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_va,vwriteSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_wap,wwriteSTD(:,:,iff-6),iff)
          IF (vars_defined) THEN
             DO k=1, nlevSTD
                DO i=1, klon
                   IF (tnondef(i,k,iff-6).NE.missing_val) THEN
                      IF (freq_outNMC(iff-6).LT.0) THEN
                         freq_moyNMC(iff-6)=(mth_len*un_jour)/freq_calNMC(iff-6)
                      ELSE
                         freq_moyNMC(iff-6)=freq_outNMC(iff-6)/freq_calNMC(iff-6)
                      ENDIF
                      zx_tmp_fi3d_STD(i,k) = (100.*tnondef(i,k,iff-6))/freq_moyNMC(iff-6)
                   ELSE
                      zx_tmp_fi3d_STD(i,k) = missing_val
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          CALL histwrite_phy(o_psbg,zx_tmp_fi3d_STD,iff)
          IF (vars_defined) THEN
             DO k=1, nlevSTD
                DO i=1, klon
                   IF (O3sumSTD(i,k,iff-6).NE.missing_val) THEN
                      zx_tmp_fi3d_STD(i,k) = O3sumSTD(i,k,iff-6) * 1.e+9
                   ELSE
                      zx_tmp_fi3d_STD(i,k) = missing_val
                   ENDIF
                ENDDO
             ENDDO !k=1, nlevSTD
          ENDIF
          CALL histwrite_phy(o_tro3,zx_tmp_fi3d_STD,iff)
          IF (read_climoz == 2) THEN
             IF (vars_defined) THEN
                DO k=1, nlevSTD
                   DO i=1, klon
                      IF (O3daysumSTD(i,k,iff-6).NE.missing_val) THEN
                         zx_tmp_fi3d_STD(i,k) = O3daysumSTD(i,k,iff-6) * 1.e+9
                      ELSE
                         zx_tmp_fi3d_STD(i,k) = missing_val
                      ENDIF
                   ENDDO
                ENDDO !k=1, nlevSTD
             ENDIF
             CALL histwrite_phy(o_tro3_daylight,zx_tmp_fi3d_STD,iff)
          ENDIF
          CALL histwrite_phy(o_uxv,uvsumSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_vxq,vqsumSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_vxT,vTsumSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_wxq,wqsumSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_vxphi,vphisumSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_wxT,wTsumSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_uxu,u2sumSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_vxv,v2sumSTD(:,:,iff-6),iff)
          CALL histwrite_phy(o_TxT,T2sumSTD(:,:,iff-6),iff)
       ENDDO !nfiles
  ENDIF
#endif
#endif
#ifdef CPP_XIOS
  IF (ok_all_xml) THEN 
!      DO iff=7, nfiles

!         CALL histwrite_phy(o_tnondef,tnondef(:,:,3))
          CALL histwrite_phy(o_ta,tlevSTD(:,:))
          CALL histwrite_phy(o_zg,philevSTD(:,:))
          CALL histwrite_phy(o_hus,qlevSTD(:,:))
          CALL histwrite_phy(o_hur,rhlevSTD(:,:))
          CALL histwrite_phy(o_ua,ulevSTD(:,:))
          CALL histwrite_phy(o_va,vlevSTD(:,:))
          CALL histwrite_phy(o_wap,wlevSTD(:,:))
!         IF (vars_defined) THEN
!            DO k=1, nlevSTD
!               DO i=1, klon
!                  IF (tnondef(i,k,3).NE.missing_val) THEN
!                     IF (freq_outNMC(iff-6).LT.0) THEN
!                        freq_moyNMC(iff-6)=(mth_len*un_jour)/freq_calNMC(iff-6)
!                     ELSE
!                        freq_moyNMC(iff-6)=freq_outNMC(iff-6)/freq_calNMC(iff-6)
!                     ENDIF
!                     zx_tmp_fi3d_STD(i,k) = (100.*tnondef(i,k,3))/freq_moyNMC(iff-6)
!                  ELSE
!                     zx_tmp_fi3d_STD(i,k) = missing_val
!                  ENDIF
!               ENDDO
!            ENDDO
!         ENDIF
!         CALL histwrite_phy(o_psbg,zx_tmp_fi3d_STD)
          IF (vars_defined) THEN
             DO k=1, nlevSTD
                DO i=1, klon
                   IF (O3STD(i,k).NE.missing_val) THEN
                      zx_tmp_fi3d_STD(i,k) = O3STD(i,k) * 1.e+9
                   ELSE
                      zx_tmp_fi3d_STD(i,k) = missing_val
                   ENDIF
                ENDDO
             ENDDO !k=1, nlevSTD
          ENDIF
          CALL histwrite_phy(o_tro3,zx_tmp_fi3d_STD)
          IF (read_climoz == 2) THEN
             IF (vars_defined) THEN
                DO k=1, nlevSTD
                   DO i=1, klon
                      IF (O3daySTD(i,k).NE.missing_val) THEN
                         zx_tmp_fi3d_STD(i,k) = O3daySTD(i,k) * 1.e+9
                      ELSE
                         zx_tmp_fi3d_STD(i,k) = missing_val
                      ENDIF
                   ENDDO
                ENDDO !k=1, nlevSTD
             ENDIF
             CALL histwrite_phy(o_tro3_daylight,zx_tmp_fi3d_STD)
          ENDIF
          CALL histwrite_phy(o_uxv,uvSTD(:,:))
          CALL histwrite_phy(o_vxq,vqSTD(:,:))
          CALL histwrite_phy(o_vxT,vTSTD(:,:))
          CALL histwrite_phy(o_wxq,wqSTD(:,:))
          CALL histwrite_phy(o_vxphi,vphiSTD(:,:))
          CALL histwrite_phy(o_wxT,wTSTD(:,:))
          CALL histwrite_phy(o_uxu,u2STD(:,:))
          CALL histwrite_phy(o_vxv,v2STD(:,:))
          CALL histwrite_phy(o_TxT,T2STD(:,:))
!      ENDDO !nfiles
  ENDIF
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (nqtot.GE.nqo+1) THEN
            DO iq=nqo+1,nqtot
              IF (type_trac == 'lmdz' .OR. type_trac == 'repr') THEN

             CALL histwrite_phy(o_trac(iq-nqo), tr_seri(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_vdf(iq-nqo),d_tr_cl(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_the(iq-nqo),d_tr_th(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_con(iq-nqo),d_tr_cv(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_lessi_impa(iq-nqo),d_tr_lessi_impa(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_lessi_nucl(iq-nqo),d_tr_lessi_nucl(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_insc(iq-nqo),d_tr_insc(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_bcscav(iq-nqo),d_tr_bcscav(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_evapls(iq-nqo),d_tr_evapls(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_ls(iq-nqo),d_tr_ls(:,:,iq-nqo))
!            CALL histwrite_phy(o_dtr_dyn(iq-nqo),d_tr_dyn(:,:,iq-nqo))
!            CALL histwrite_phy(o_dtr_cl(iq-nqo),d_tr_cl(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_trsp(iq-nqo),d_tr_trsp(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_sscav(iq-nqo),d_tr_sscav(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_sat(iq-nqo),d_tr_sat(:,:,iq-nqo))
             CALL histwrite_phy(o_dtr_uscav(iq-nqo),d_tr_uscav(:,:,iq-nqo))
             zx_tmp_fi2d=0.
             IF (vars_defined) THEN
                DO k=1,klev
                   zx_tmp_fi2d(:)=zx_tmp_fi2d(:)+zmasse(:,k)*tr_seri(:,k,iq-nqo)
                ENDDO
             ENDIF
             CALL histwrite_phy(o_trac_cum(iq-nqo), zx_tmp_fi2d)
             ENDIF
          ENDDO
       ENDIF

       IF (.NOT.vars_defined) THEN
          !$OMP MASTER
#ifndef CPP_IOIPSL_NO_OUTPUT 
          DO iff=1,nfiles
             IF (clef_files(iff)) THEN
                CALL histend(nid_files(iff))
                ndex2d = 0
                ndex3d = 0

             ENDIF ! clef_files
          ENDDO !  iff
#endif
#ifdef CPP_XIOS
          !On finalise l'initialisation:
          CALL wxios_closedef()
#endif

          !$OMP END MASTER
          !$OMP BARRIER
          vars_defined = .TRUE.

       ENDIF

    ENDDO

    IF (vars_defined) THEN
       ! On synchronise les fichiers pour IOIPSL
#ifndef CPP_IOIPSL_NO_OUTPUT 
       !$OMP MASTER
       DO iff=1,nfiles
          IF (ok_sync .AND. clef_files(iff)) THEN
             CALL histsync(nid_files(iff))
          ENDIF
       ENDDO
       !$OMP END MASTER
#endif
    ENDIF

  END SUBROUTINE phys_output_write_spl

END MODULE phys_output_write_spl_mod
