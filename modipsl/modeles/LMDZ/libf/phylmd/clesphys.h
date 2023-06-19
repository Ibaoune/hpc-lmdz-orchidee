! $Id: clesphys.h 3435 2019-01-22 15:21:59Z fairhead $
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez \`a n'utiliser que des ! pour les commentaires
!                 et \`a bien positionner les & des lignes de continuation 
!                 (les placer en colonne 6 et en colonne 73)
!
!..include cles_phys.h
!
       INTEGER iflag_cycle_diurne
       LOGICAL soil_model,new_oliq,ok_orodr,ok_orolf 
       LOGICAL ok_limitvrai
       LOGICAL ok_all_xml
       LOGICAL ok_lwoff
       INTEGER nbapp_rad, iflag_con, nbapp_cv, nbapp_wk, iflag_ener_conserv
       REAL co2_ppm, co2_ppm0, solaire
!FC
       REAL Cd_frein
       LOGICAL ok_suntime_rrtm
       REAL(kind=8) RCO2, RCH4, RN2O, RCFC11, RCFC12  
       REAL(kind=8) RCO2_act, RCH4_act, RN2O_act, RCFC11_act, RCFC12_act  
       REAL(kind=8) CH4_ppb, N2O_ppb, CFC11_ppt, CFC12_ppt
!IM ajout CFMIP2/CMIP5
       REAL(kind=8) RCO2_per,RCH4_per,RN2O_per,RCFC11_per,RCFC12_per
       REAL(kind=8) CH4_ppb_per,N2O_ppb_per,CFC11_ppt_per,CFC12_ppt_per

!OM ---> correction du bilan d'eau global
!OM Correction sur precip KE
       REAL cvl_corr
!OM Fonte calotte dans bilan eau
       LOGICAL ok_lic_melt
!OB Depot de vapeur d eau sur la calotte pour le bilan eau
       LOGICAL ok_lic_cond

!IM simulateur ISCCP 
       INTEGER top_height, overlap
!IM seuils cdrm, cdrh
       REAL cdmmax, cdhmax
!IM param. stabilite s/ terres et en dehors
       REAL ksta, ksta_ter, f_ri_cd_min
!IM ok_kzmin : clef calcul Kzmin dans la CL de surface cf FH
       LOGICAL ok_kzmin
!IM, MAFo fmagic, pmagic : parametres - additionnel et multiplicatif - 
!                          pour regler l albedo sur ocean
       REAL pbl_lmixmin_alpha
       REAL fmagic, pmagic
! Hauteur (imposee) du contenu en eau du sol
           REAL qsol0,albsno0,evap0
! Frottement au sol (Cdrag)
       Real f_cdrag_ter,f_cdrag_oce
       REAL min_wind_speed,f_gust_wk,f_gust_bl,f_qsat_oce,f_z0qh_oce
       REAL z0m_seaice,z0h_seaice
       INTEGER iflag_gusts,iflag_z0_oce

! Rugoro
       Real f_rugoro,z0min

! tau_gl : constante de rappel de la temperature a la surface de la glace
       REAL tau_gl

!IM lev_histhf  : niveau sorties 6h
!IM lev_histday : niveau sorties journalieres
!IM lev_histmth : niveau sorties mensuelles
!IM lev_histdayNMC : on peut sortir soit sur 8 (comme AR5) ou bien 
!                    sur 17 niveaux de pression
       INTEGER lev_histhf, lev_histday, lev_histmth
       INTEGER lev_histdayNMC
       Integer lev_histins, lev_histLES  
!IM ok_histNMC  : sortie fichiers niveaux de pression (histmthNMC, histdayNMC, histhfNMC)
!IM freq_outNMC : frequences de sortie fichiers niveaux de pression (histmthNMC, histdayNMC, histhfNMC)
!IM freq_calNMC : frequences de calcul fis. hist*NMC.nc
       LOGICAL ok_histNMC(3)
       INTEGER levout_histNMC(3)
       REAL freq_outNMC(3) , freq_calNMC(3)
       CHARACTER(len=4) type_run
! aer_type: pour utiliser un fichier constant dans readaerosol 
       CHARACTER(len=8) :: aer_type 
       LOGICAL ok_regdyn
       REAL lonmin_ins, lonmax_ins, latmin_ins, latmax_ins
       REAL ecrit_ins, ecrit_hf, ecrit_day
       REAL ecrit_mth, ecrit_tra, ecrit_reg 
       REAL ecrit_LES
       REAL freq_ISCCP, ecrit_ISCCP
       REAL freq_COSP, freq_AIRS
       LOGICAL :: ok_cosp,ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP
       LOGICAL :: ok_airs
       INTEGER :: ip_ebil_phy, iflag_rrtm, iflag_ice_thermo, NSW, iflag_albedo
       LOGICAL :: ok_chlorophyll
       LOGICAL :: ok_strato
       LOGICAL :: ok_hines, ok_gwd_rando
       LOGICAL :: ok_qch4
       LOGICAL :: ok_conserv_q
       LOGICAL :: adjust_tropopause
       LOGICAL :: ok_daily_climoz
! flag to bypass or not the phytrac module
       INTEGER :: iflag_phytrac

       COMMON/clesphys/                                                 &
! REAL FIRST
     &       co2_ppm, solaire                                           &
     &     , RCO2, RCH4, RN2O, RCFC11, RCFC12                           &
     &     , RCO2_act, RCH4_act, RN2O_act, RCFC11_act, RCFC12_act       &
     &     , RCO2_per, RCH4_per, RN2O_per, RCFC11_per, RCFC12_per       &
     &     , CH4_ppb, N2O_ppb, CFC11_ppt, CFC12_ppt                     &
     &     , CH4_ppb_per, N2O_ppb_per, CFC11_ppt_per, CFC12_ppt_per     &
     &     , cdmmax,cdhmax,ksta,ksta_ter,f_ri_cd_min,pbl_lmixmin_alpha  &
     &     , fmagic, pmagic                                             &
     &     , f_cdrag_ter,f_cdrag_oce,f_rugoro,z0min,tau_gl              &
     &     , min_wind_speed,f_gust_wk,f_gust_bl,f_qsat_oce,f_z0qh_oce   &
     &     , z0m_seaice,z0h_seaice                                      &
     &     , freq_outNMC, freq_calNMC                                   &
     &     , lonmin_ins, lonmax_ins, latmin_ins, latmax_ins             &
     &     , freq_ISCCP, ecrit_ISCCP, freq_COSP, freq_AIRS              &
     &     , cvl_corr                                                   &
     &     , qsol0,albsno0,evap0                                        &
     &     , co2_ppm0                                                   &
!FC
     &     , Cd_frein                                                   &
     &     , ecrit_LES                                                  &
     &     , ecrit_ins, ecrit_hf, ecrit_day                             &
     &     , ecrit_mth, ecrit_tra, ecrit_reg                            &
! THEN INTEGER AND LOGICALS
     &     , top_height                                                 &
     &     , iflag_cycle_diurne, soil_model, new_oliq                         &
     &     , ok_orodr, ok_orolf, ok_limitvrai, nbapp_rad                &
     &     , iflag_con, nbapp_cv, nbapp_wk                              &
     &     , iflag_ener_conserv                                         &
     &     , ok_suntime_rrtm                                            & 
     &     , overlap                                                    &
     &     , ok_kzmin                                                   &
     &     , lev_histhf, lev_histday, lev_histmth                       &
     &     , lev_histins, lev_histLES, lev_histdayNMC, levout_histNMC   &
     &     , ok_histNMC                                                 &
     &     , type_run, ok_regdyn, ok_cosp, ok_airs                      &
     &     , ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP                     &
     &     , ip_ebil_phy                                                &
     &     , iflag_gusts ,iflag_z0_oce                                  &
     &     , ok_lic_melt, ok_lic_cond, aer_type                         &
     &     , iflag_rrtm, ok_strato,ok_hines, ok_qch4                    &
     &     , iflag_ice_thermo, ok_gwd_rando, NSW, iflag_albedo          &
     &     , ok_chlorophyll,ok_conserv_q, adjust_tropopause             &
     &     , ok_daily_climoz, ok_all_xml, ok_lwoff                      &
     &     , iflag_phytrac
     
       save /clesphys/
!$OMP THREADPRIVATE(/clesphys/)
