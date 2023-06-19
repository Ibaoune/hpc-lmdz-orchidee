#include "netcdf.inc"

! Declarations specifiques au cas Toga
        character*80 :: fich_toga
!        integer nlev_prof
!        parameter (nlev_prof = 41)
        integer nlev_toga, nt_toga
        parameter (nlev_toga=41, nt_toga=480)
        integer year_ini_toga, day_ini_toga, mth_ini_toga
        real day_ju_ini_toga   ! Julian day of toga coare first day
        parameter (year_ini_toga=1992) 
        parameter (mth_ini_toga=11)
        parameter (day_ini_toga=1)  !  1erNov1992
        real dt_toga
        parameter (dt_toga=6.*3600.)
!!
        integer year_print, month_print, day_print
        real    sec_print
!!
        real ts_toga(nt_toga)
        real plev_toga(nlev_toga,nt_toga),w_toga(nlev_toga,nt_toga)
        real t_toga(nlev_toga,nt_toga),q_toga(nlev_toga,nt_toga)
        real u_toga(nlev_toga,nt_toga),v_toga(nlev_toga,nt_toga)
        real ht_toga(nlev_toga,nt_toga),vt_toga(nlev_toga,nt_toga)
        real hq_toga(nlev_toga,nt_toga),vq_toga(nlev_toga,nt_toga)

        real ts_prof
        real plev_prof(nlev_toga),w_prof(nlev_toga)
        real t_prof(nlev_toga),q_prof(nlev_toga)
        real u_prof(nlev_toga),v_prof(nlev_toga)
        real ht_prof(nlev_toga),vt_prof(nlev_toga)
        real hq_prof(nlev_toga),vq_prof(nlev_toga)

        real w_mod(llm), t_mod(llm),q_mod(llm)
        real u_mod(llm),v_mod(llm), ht_mod(llm),vt_mod(llm),ug_mod(llm),vg_mod(llm)
	real temp_nudg_mod(llm),qv_nudg_mod(llm),u_nudg_mod(llm),v_nudg_mod(llm)
        real hq_mod(llm),vq_mod(llm),qv_mod(llm),ql_mod(llm),qt_mod(llm)
        real th_mod(llm)

        real ts_cur
        common /sst_forcing/ts_cur ! also in read_tsurf1d.F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declarations specifiques au cas RICO
        character*80 :: fich_rico
        integer nlev_rico

        parameter (nlev_rico=81)
        real ts_rico,ps_rico
        real w_rico(llm)
        real t_rico(llm),q_rico(llm)
        real u_rico(llm),v_rico(llm)
        real dth_rico(llm)
        real dqh_rico(llm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declarations specifiques au cas TWPice
        character*80 :: fich_twpice
        integer nlev_twpi, nt_twpi
        parameter (nlev_twpi=40, nt_twpi=215)
        integer year_ini_twpi, day_ini_twpi, mth_ini_twpi
        real heure_ini_twpi
        real day_ju_ini_twpi   ! Julian day of twpice first day
        parameter (year_ini_twpi=2006) 
        parameter (mth_ini_twpi=1)
        parameter (day_ini_twpi=17)  ! 17 = 17Jan2006
        parameter (heure_ini_twpi=10800.) !3h en secondes
        real dt_twpi
        parameter (dt_twpi=3.*3600.)

        real ts_twpi(nt_twpi)
        real plev_twpi(nlev_twpi,nt_twpi),w_twpi(nlev_twpi,nt_twpi)
        real t_twpi(nlev_twpi,nt_twpi),q_twpi(nlev_twpi,nt_twpi)
        real u_twpi(nlev_twpi,nt_twpi),v_twpi(nlev_twpi,nt_twpi)
        real ht_twpi(nlev_twpi,nt_twpi),vt_twpi(nlev_twpi,nt_twpi)
        real hq_twpi(nlev_twpi,nt_twpi),vq_twpi(nlev_twpi,nt_twpi)

        real ts_proftwp
        real plev_proftwp(nlev_twpi),w_proftwp(nlev_twpi)
        real t_proftwp(nlev_twpi),q_proftwp(nlev_twpi)
        real u_proftwp(nlev_twpi),v_proftwp(nlev_twpi)
        real ht_proftwp(nlev_twpi),vt_proftwp(nlev_twpi)
        real hq_proftwp(nlev_twpi),vq_proftwp(nlev_twpi)



!Declarations specifiques au cas FIRE
        character*80 :: fich_fire
        integer nlev_fire, nt_fire
        parameter (nlev_fire=120, nt_fire=1)  
        integer year_ini_fire, day_ini_fire, mth_ini_fire
        real heure_ini_fire
        parameter (year_ini_fire=1987) 
        parameter (mth_ini_fire=7)
        parameter (day_ini_fire=14)  ! 14 = 14Juil1987
        parameter (heure_ini_fire=0.) !0h en secondes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Declarations specifiques au cas GABLS4   (MPL 20141023)
!FHADETRUIRE
!       character*80 :: fich_gabls4
!       integer nlev_gabls4, nt_gabls4, nsol_gabls4
!       parameter (nlev_gabls4=90, nt_gabls4=37, nsol_gabls4=19)  
!       integer year_ini_gabls4, day_ini_gabls4, mth_ini_gabls4
!       real heure_ini_gabls4
!       real day_ju_ini_gabls4   ! Julian day of gabls4 first day
!       parameter (year_ini_gabls4=2009) 
!       parameter (mth_ini_gabls4=12)
!       parameter (day_ini_gabls4=11)  ! 11 = 11 decembre 2009
!       parameter (heure_ini_gabls4=0.) !0UTC en secondes
!       real dt_gabls4
!       parameter (dt_gabls4=3600.) ! 1 forcage ttes les heures
!
!profils initiaux:
!       real plev_gabls4(nlev_gabls4)
!       real zz_gabls4(nlev_gabls4)
!       real th_gabls4(nlev_gabls4),t_gabls4(nlev_gabls4),qv_gabls4(nlev_gabls4)
!       real u_gabls4(nlev_gabls4), v_gabls4(nlev_gabls4)
!       real depth_sn_gabls4(nsol_gabls4),tsnow_gabls4(nsol_gabls4),snow_dens_gabls4(nsol_gabls4)
!       real t_gabi(nlev_gabls4),qv_gabi(nlev_gabls4)
!       real u_gabi(nlev_gabls4), v_gabi(nlev_gabls4),ug_gabi(nlev_gabls4), vg_gabi(nlev_gabls4)
!       real ht_gabi(nlev_gabls4),hq_gabi(nlev_gabls4),poub(nlev_gabls4)
!       
!forcings
! Lignes a detruire ...
!       real ht_gabls4(nlev_gabls4,nt_gabls4),hq_gabls4(nlev_gabls4,nt_gabls4)
!       real ug_gabls4(nlev_gabls4,nt_gabls4),vg_gabls4(nlev_gabls4,nt_gabls4)
!       real tg_gabls4(nt_gabls4)
!       real ht_profg(nlev_gabls4),hq_profg(nlev_gabls4)
!       real ug_profg(nlev_gabls4),vg_profg(nlev_gabls4)
!       real tg_profg
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Declarations specifiques au cas DICE     (MPL 02072013)
        character*80 :: fich_dice
        integer nlev_dice, nt_dice
        parameter (nlev_dice=70, nt_dice=145)  
        integer year_ini_dice, day_ini_dice, mth_ini_dice
        real heure_ini_dice
        real day_ju_ini_dice   ! Julian day of dice first day
        parameter (year_ini_dice=1999) 
        parameter (mth_ini_dice=10)
        parameter (day_ini_dice=23)  ! 23 = 23 october 1999
        parameter (heure_ini_dice=68400.) !19UTC en secondes
        real dt_dice
        parameter (dt_dice=0.5*3600.) ! 1 forcage ttes les demi-heures

!profils initiaux:
        real plev_dice(nlev_dice)
        
        real zz_dice(nlev_dice)
        real t_dice(nlev_dice),qv_dice(nlev_dice)
        real u_dice(nlev_dice), v_dice(nlev_dice),o3_dice(nlev_dice)
        real ht_dice(nlev_dice,nt_dice)
        real hq_dice(nlev_dice,nt_dice), hu_dice(nlev_dice,nt_dice)
        real hv_dice(nlev_dice,nt_dice)
        real w_dice(nlev_dice,nt_dice),omega_dice(nlev_dice,nt_dice)
        real o3_mod(llm),hu_mod(llm),hv_mod(llm)
        real t_dicei(nlev_dice),qv_dicei(nlev_dice)
        real u_dicei(nlev_dice), v_dicei(nlev_dice),o3_dicei(nlev_dice)
        real ht_dicei(nlev_dice)
        real hq_dicei(nlev_dice), hu_dicei(nlev_dice)
        real hv_dicei(nlev_dice)
        real w_dicei(nlev_dice),omega_dicei(nlev_dice)

        
!forcings
        real shf_dice(nt_dice),lhf_dice(nt_dice)
        real lwup_dice(nt_dice),swup_dice(nt_dice)
        real tg_dice(nt_dice),ustar_dice(nt_dice),psurf_dice(nt_dice)
        real ug_dice(nt_dice),vg_dice(nt_dice)

        real shf_prof,lhf_prof,lwup_prof,swup_prof,tg_prof
        real ustar_prof,psurf_prof,cdrag
        real ht_profd(nlev_dice),hq_profd(nlev_dice),hu_profd(nlev_dice)
        real hv_profd(nlev_dice),w_profd(nlev_dice)
        real omega_profd(nlev_dice),ug_profd,vg_profd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declarations specifiques au cas GCSSold
        character*80 :: fich_gcssold_ctl
        character*80 :: fich_gcssold_dat
        real  ht_gcssold(llm),hq_gcssold(llm),hw_gcssold(llm)
        real  hu_gcssold(llm)
        real  hv_gcssold(llm)
        real  hthturb_gcssold(llm)
        real  hqturb_gcssold(llm)
        real  Ts_gcssold
        real  dtime_frcg
        logical :: Turb_fcg_gcssold

        common /turb_forcing/                                                   &
     &  dtime_frcg,hthturb_gcssold, hqturb_gcssold,Turb_fcg_gcssold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declarations specifiques au cas Arm_cu
        character*80 :: fich_armcu


        integer nlev_armcu, nt_armcu
        parameter (nlev_armcu=40, nt_armcu=31)
        integer year_ini_armcu, day_ini_armcu, mth_ini_armcu
        real  heure_ini_armcu
        real day_ju_ini_armcu                                ! Julian day of armcu case first day
        parameter (year_ini_armcu=1997) 
        parameter (mth_ini_armcu=6)
        parameter (day_ini_armcu=21)  ! 172 = 21 juin 1997
        parameter (heure_ini_armcu=41400)   ! 11:30 en secondes
        real dt_armcu
        parameter (dt_armcu=1.*1800.)   ! forcages donnes ttes les demi-heures par ifa_armcu.txt
        real sens_armcu(nt_armcu),flat_armcu(nt_armcu)
        real adv_theta_armcu(nt_armcu),rad_theta_armcu(nt_armcu)
        real adv_qt_armcu(nt_armcu)
        real theta_mod(llm),rv_mod(llm),play_mod(llm)
! profc comme "profil armcu"
        
! forcages interpoles dans le temps
        real adv_theta_prof,rad_theta_prof,adv_qt_prof
        real sens_prof,flat_prof,fact
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! declarations specifiques au cas Sandu
        character*80 :: fich_sandu
!        integer nlev_prof
!        parameter (nlev_prof = 41)
        integer nlev_sandu, nt_sandu
        parameter (nlev_sandu=87, nt_sandu=13)
        integer year_ini_sandu, day_ini_sandu, mth_ini_sandu
        real day_ju_ini_sandu                                ! Julian day of sandu case first day
        parameter (year_ini_sandu=2006)
        parameter (mth_ini_sandu=7)
        parameter (day_ini_sandu=15)  ! 196 = 15 juillet 2006
        real dt_sandu, tau_sandu
        logical  :: trouve_700=.true.
        parameter (dt_sandu=6.*3600.)   ! forcages donnes ttes les 6 heures par ifa_sandu.txt
        parameter (tau_sandu=30000*3600.)  ! temps de relaxation u,v,thetal,qt vers profil init et au dessus 700hPa
!!
        real ts_sandu(nt_sandu)
! profs comme "profil sandu"
        real plev_profs(nlev_sandu)
        real t_profs(nlev_sandu),thl_profs(nlev_sandu)
        real q_profs(nlev_sandu)
        real u_profs(nlev_sandu),v_profs(nlev_sandu),w_profs(nlev_sandu)
        real omega_profs(nlev_sandu),o3mmr_profs(nlev_sandu)

        real, dimension(llm) :: relax_u,relax_v,relax_thl
        real, dimension(llm,2) :: relax_q

        real thl_mod(llm),omega_mod(llm),o3mmr_mod(llm),tke_mod(llm)
!vertical advection computation
        real d_t_z(llm),d_th_z(llm), d_q_z(llm)
        real d_t_dyn_z(llm),d_th_dyn_z(llm), d_q_dyn_z(llm)
        real d_u_z(llm),d_v_z(llm)
        real d_u_dyn(llm),d_v_dyn(llm)
        real d_u_dyn_z(llm),d_v_dyn_z(llm)
        real d_u_adv(llm),d_v_adv(llm)
        real zz(llm)
        real zfact
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declarations specifiques au cas Astex
        character*80 :: fich_astex
        integer nlev_astex, nt_astex
        parameter (nlev_astex=34, nt_astex=49)
        integer year_ini_astex, day_ini_astex, mth_ini_astex
        real day_ju_ini_astex                                ! Julian day of astex case first day
        parameter (year_ini_astex=1992)
        parameter (mth_ini_astex=6)
        parameter (day_ini_astex=13)  ! 165 = 13 juin 1992
        real dt_astex
        parameter (dt_astex=3600.)    ! forcages donnes ttes les heures par ifa_astex.txt
        real ts_astex(nt_astex),div_astex(nt_astex),ug_astex(nt_astex)
        real vg_astex(nt_astex),ufa_astex(nt_astex),vfa_astex(nt_astex)
        real div_prof,ug_prof,vg_prof,ufa_prof,vfa_prof
! profa comme "profil astex"
        real plev_profa(nlev_astex)
        real t_profa(nlev_astex),thl_profa(nlev_astex)
        real qv_profa(nlev_astex),ql_profa(nlev_astex)
        real qt_profa(nlev_astex),o3mmr_profa(nlev_astex)
        real u_profa(nlev_astex),v_profa(nlev_astex),w_profa(nlev_astex)
        real tke_profa(nlev_astex)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Declarations specifiques au cas standard

        real w_mod_cas(llm), t_mod_cas(llm),q_mod_cas(llm)
        real theta_mod_cas(llm),thl_mod_cas(llm),thv_mod_cas(llm)
        real qv_mod_cas(llm),ql_mod_cas(llm),qi_mod_cas(llm)
        real ug_mod_cas(llm),vg_mod_cas(llm)
        real temp_nudg_mod_cas(llm),qv_nudg_mod_cas(llm),v_nudg_mod_cas(llm),u_nudg_mod_cas(llm)
        real u_mod_cas(llm),v_mod_cas(llm)
        real omega_mod_cas(llm)
        real ht_mod_cas(llm),vt_mod_cas(llm),dt_mod_cas(llm),dtrad_mod_cas(llm)
        real hth_mod_cas(llm),vth_mod_cas(llm),dth_mod_cas(llm)
        real hq_mod_cas(llm),vq_mod_cas(llm),dq_mod_cas(llm)
        real hu_mod_cas(llm),vu_mod_cas(llm),du_mod_cas(llm)
        real hv_mod_cas(llm),vv_mod_cas(llm),dv_mod_cas(llm)
        integer day_ini_cas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


