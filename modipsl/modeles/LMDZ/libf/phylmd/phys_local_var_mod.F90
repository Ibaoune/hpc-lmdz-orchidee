!
! $Id: phys_local_var_mod.F90 3823 2021-02-05 09:00:13Z musat $
!
      MODULE phys_local_var_mod
! Variables locales pour effectuer les appels en serie
!======================================================================
!
!
!======================================================================
! Declaration des variables

      REAL, SAVE, ALLOCATABLE :: t_seri(:,:), q_seri(:,:)
      !$OMP THREADPRIVATE(t_seri, q_seri)
      REAL, SAVE, ALLOCATABLE :: ql_seri(:,:),qs_seri(:,:)
      !$OMP THREADPRIVATE(ql_seri,qs_seri)
      REAL, SAVE, ALLOCATABLE :: u_seri(:,:), v_seri(:,:)
      !$OMP THREADPRIVATE(u_seri, v_seri)
      REAL, SAVE, ALLOCATABLE :: l_mixmin(:,:,:), l_mix(:,:,:)
      !$OMP THREADPRIVATE(l_mixmin, l_mix)

      REAL, SAVE, ALLOCATABLE :: tr_seri(:,:,:)
      !$OMP THREADPRIVATE(tr_seri)
      REAL, SAVE, ALLOCATABLE :: d_t_dyn(:,:), d_q_dyn(:,:)
      !$OMP THREADPRIVATE(d_t_dyn, d_q_dyn)
      REAL, SAVE, ALLOCATABLE :: d_ql_dyn(:,:), d_qs_dyn(:,:)
      !$OMP THREADPRIVATE(d_ql_dyn, d_qs_dyn)
      REAL, SAVE, ALLOCATABLE :: d_q_dyn2d(:), d_ql_dyn2d(:), d_qs_dyn2d(:)
      !$OMP THREADPRIVATE(d_q_dyn2d, d_ql_dyn2d, d_qs_dyn2d)
      REAL, SAVE, ALLOCATABLE :: d_u_dyn(:,:), d_v_dyn(:,:)
      !$OMP THREADPRIVATE(d_u_dyn, d_v_dyn)
      REAL, SAVE, ALLOCATABLE :: d_tr_dyn(:,:,:)
      !$OMP THREADPRIVATE(d_tr_dyn)
      REAL, SAVE, ALLOCATABLE :: d_t_con(:,:),d_q_con(:,:)
      !$OMP THREADPRIVATE(d_t_con,d_q_con)
      REAL, SAVE, ALLOCATABLE :: d_u_con(:,:),d_v_con(:,:)
      !$OMP THREADPRIVATE(d_u_con,d_v_con)
      REAL, SAVE, ALLOCATABLE :: d_t_wake(:,:),d_q_wake(:,:)
      !$OMP THREADPRIVATE( d_t_wake,d_q_wake)
      REAL, SAVE, ALLOCATABLE :: d_t_lsc(:,:),d_q_lsc(:,:),d_ql_lsc(:,:),d_qi_lsc(:,:)
      !$OMP THREADPRIVATE(d_t_lsc,d_q_lsc,d_ql_lsc,d_qi_lsc)
      REAL, SAVE, ALLOCATABLE :: d_t_lwr(:,:),d_t_lw0(:,:),d_t_swr(:,:),d_t_sw0(:,:)
      !$OMP THREADPRIVATE(d_t_lwr,d_t_lw0,d_t_swr,d_t_sw0)
      REAL, SAVE, ALLOCATABLE :: d_t_ajsb(:,:), d_q_ajsb(:,:)
      !$OMP THREADPRIVATE(d_t_ajsb, d_q_ajsb)
      REAL, SAVE, ALLOCATABLE :: d_t_ajs(:,:), d_q_ajs(:,:)
      !$OMP THREADPRIVATE(d_t_ajs, d_q_ajs)
      REAL, SAVE, ALLOCATABLE :: d_u_ajs(:,:), d_v_ajs(:,:)
      !$OMP THREADPRIVATE(d_u_ajs, d_v_ajs)
!nrlmd<
      REAL, SAVE, ALLOCATABLE :: d_t_ajs_w(:,:), d_q_ajs_w(:,:)
      !$OMP THREADPRIVATE(d_t_ajs_w, d_q_ajs_w)
      REAL, SAVE, ALLOCATABLE :: d_t_ajs_x(:,:), d_q_ajs_x(:,:)
      !$OMP THREADPRIVATE(d_t_ajs_x, d_q_ajs_x)
!>nrlmd
      REAL, SAVE, ALLOCATABLE :: d_t_eva(:,:),d_q_eva(:,:),d_ql_eva(:,:),d_qi_eva(:,:)
      !$OMP THREADPRIVATE(d_t_eva,d_q_eva,d_ql_eva,d_qi_eva)
      REAL, SAVE, ALLOCATABLE :: d_t_lscst(:,:),d_q_lscst(:,:)
      !$OMP THREADPRIVATE(d_t_lscst,d_q_lscst)
      REAL, SAVE, ALLOCATABLE :: d_t_lscth(:,:),d_q_lscth(:,:)
      !$OMP THREADPRIVATE(d_t_lscth,d_q_lscth)
      REAL, SAVE, ALLOCATABLE :: plul_th(:),plul_st(:)
      !$OMP THREADPRIVATE(plul_th,plul_st)
!tendances dues a oro et lif
      REAL, SAVE, ALLOCATABLE :: d_t_oli(:,:)
      !$OMP THREADPRIVATE(d_t_oli)
      REAL, SAVE, ALLOCATABLE :: d_u_oli(:,:), d_v_oli(:,:)
      !$OMP THREADPRIVATE(d_u_oli, d_v_oli)
      REAL, SAVE, ALLOCATABLE :: d_t_vdf(:,:), d_q_vdf(:,:), d_t_diss(:,:)
      !$OMP THREADPRIVATE( d_t_vdf, d_q_vdf,d_t_diss)
      REAL, SAVE, ALLOCATABLE :: d_u_vdf(:,:), d_v_vdf(:,:)
      !$OMP THREADPRIVATE(d_u_vdf, d_v_vdf)
!nrlmd+jyg<
      REAL, SAVE, ALLOCATABLE :: d_t_vdf_w(:,:), d_q_vdf_w(:,:)
      !$OMP THREADPRIVATE( d_t_vdf_w, d_q_vdf_w)
      REAL, SAVE, ALLOCATABLE :: d_t_vdf_x(:,:), d_q_vdf_x(:,:)
      !$OMP THREADPRIVATE( d_t_vdf_x, d_q_vdf_x)
!>nrlmd+jyg
      REAL, SAVE, ALLOCATABLE :: d_t_oro(:,:)
      !$OMP THREADPRIVATE(d_t_oro)
      REAL, SAVE, ALLOCATABLE :: d_u_oro(:,:), d_v_oro(:,:)
      !$OMP THREADPRIVATE(d_u_oro, d_v_oro)
      REAL, SAVE, ALLOCATABLE :: d_t_oro_gw(:,:)
      !$OMP THREADPRIVATE(d_t_oro_gw)
      REAL, SAVE, ALLOCATABLE :: d_u_oro_gw(:,:), d_v_oro_gw(:,:)
      !$OMP THREADPRIVATE(d_u_oro_gw, d_v_oro_gw)
      REAL, SAVE, ALLOCATABLE :: d_t_lif(:,:)
      !$OMP THREADPRIVATE(d_t_lif)
      REAL, SAVE, ALLOCATABLE :: d_u_lif(:,:), d_v_lif(:,:)
      !$OMP THREADPRIVATE(d_u_lif, d_v_lif)
! Tendances Ondes de G non oro (runs strato).
      REAL, SAVE, ALLOCATABLE :: du_gwd_hines(:,:)
      !$OMP THREADPRIVATE(du_gwd_hines)
      REAL, SAVE, ALLOCATABLE :: dv_gwd_hines(:,:)
      !$OMP THREADPRIVATE(dv_gwd_hines)
      REAL, SAVE, ALLOCATABLE :: dv_gwd_rando(:,:)
      !$OMP THREADPRIVATE(dv_gwd_rando)
      REAL, SAVE, ALLOCATABLE :: dv_gwd_front(:,:)
      !$OMP THREADPRIVATE(dv_gwd_front)
      REAL, SAVE, ALLOCATABLE :: east_gwstress(:,:)
      !$OMP THREADPRIVATE(east_gwstress)
      REAL, SAVE, ALLOCATABLE :: west_gwstress(:,:)
      !$OMP THREADPRIVATE(west_gwstress)
      REAL, SAVE, ALLOCATABLE :: d_t_hin(:,:)
      !$OMP THREADPRIVATE(d_t_hin)
! tendance due a l'oxydation du methane
      REAL, SAVE, ALLOCATABLE :: d_q_ch4(:,:)
      !$OMP THREADPRIVATE(d_q_ch4)

! tendance du a la conersion Ec -> E thermique
      REAL, SAVE, ALLOCATABLE :: d_t_ec(:,:)
      !$OMP THREADPRIVATE(d_t_ec)
      REAL, SAVE, ALLOCATABLE :: d_ts(:,:), d_tr(:,:,:)
      !$OMP THREADPRIVATE(d_ts, d_tr)

! diagnostique pour le rayonnement
      REAL, SAVE, ALLOCATABLE :: topswad_aero(:),  solswad_aero(:)      ! diag
      !$OMP THREADPRIVATE(topswad_aero,solswad_aero)
      REAL, SAVE, ALLOCATABLE :: topswai_aero(:),  solswai_aero(:)      ! diag
      !$OMP THREADPRIVATE(topswai_aero,solswai_aero)
      REAL, SAVE, ALLOCATABLE :: topswad0_aero(:), solswad0_aero(:)     ! diag
      !$OMP THREADPRIVATE(topswad0_aero,solswad0_aero)
      REAL, SAVE, ALLOCATABLE :: topsw_aero(:,:),  solsw_aero(:,:)      ! diag
      !$OMP THREADPRIVATE(topsw_aero,solsw_aero)
      REAL, SAVE, ALLOCATABLE :: topsw0_aero(:,:), solsw0_aero(:,:)     ! diag
      !$OMP THREADPRIVATE(topsw0_aero,solsw0_aero)
      REAL, SAVE, ALLOCATABLE :: topswcf_aero(:,:),  solswcf_aero(:,:)  ! diag
      !$OMP THREADPRIVATE(topswcf_aero,solswcf_aero)
! LW radiation diagnostics CK
      REAL, SAVE, ALLOCATABLE :: toplwad_aero(:),  sollwad_aero(:)      ! diag
      !$OMP THREADPRIVATE(toplwad_aero,sollwad_aero)
      REAL, SAVE, ALLOCATABLE :: toplwai_aero(:),  sollwai_aero(:)      ! diag
      !$OMP THREADPRIVATE(toplwai_aero,sollwai_aero)
      REAL, SAVE, ALLOCATABLE :: toplwad0_aero(:), sollwad0_aero(:)     ! diag
      !$OMP THREADPRIVATE(toplwad0_aero,sollwad0_aero)
! Special RRTM
      REAL, SAVE, ALLOCATABLE :: ZLWFT0_i(:,:),  ZSWFT0_i(:,:)      ! diag
      !$OMP THREADPRIVATE(ZLWFT0_i,ZSWFT0_i)
      REAL, SAVE, ALLOCATABLE :: ZFLDN0(:,:),  ZFLUP0(:,:)      ! diag
      !$OMP THREADPRIVATE(ZFLDN0,ZFLUP0)
      REAL, SAVE, ALLOCATABLE :: ZFSDN0(:,:),  ZFSUP0(:,:)      ! diag
      !$OMP THREADPRIVATE(ZFSDN0,ZFSUP0)
!
      REAL, SAVE, ALLOCATABLE :: stratomask(:,:) 
      !$OMP THREADPRIVATE(stratomask) 
      REAL, SAVE, ALLOCATABLE :: tausum_aero(:,:,:) 
      !$OMP THREADPRIVATE(tausum_aero) 
      REAL, SAVE, ALLOCATABLE :: drytausum_aero(:,:) 
      !$OMP THREADPRIVATE(drytausum_aero) 
      REAL, SAVE, ALLOCATABLE :: tau3d_aero(:,:,:,:) 
      !$OMP THREADPRIVATE(tau3d_aero) 
      REAL, SAVE, ALLOCATABLE :: scdnc(:,:)
      !$OMP THREADPRIVATE(scdnc)
      REAL, SAVE, ALLOCATABLE :: dNovrN(:) 
      !$OMP THREADPRIVATE(dNovrN) 
      REAL, SAVE, ALLOCATABLE :: cldncl(:)
      !$OMP THREADPRIVATE(cldncl)
      REAL, SAVE, ALLOCATABLE :: reffclwtop(:)
      !$OMP THREADPRIVATE(reffclwtop)
      REAL, SAVE, ALLOCATABLE :: lcc(:)
      !$OMP THREADPRIVATE(lcc)
      REAL, SAVE, ALLOCATABLE :: reffclws(:,:)
      !$OMP THREADPRIVATE(reffclws)
      REAL, SAVE, ALLOCATABLE :: reffclwc(:,:)
      !$OMP THREADPRIVATE(reffclwc)
      REAL, SAVE, ALLOCATABLE :: cldnvi(:) 
      !$OMP THREADPRIVATE(cldnvi)
      REAL, SAVE, ALLOCATABLE :: lcc3d(:,:)
      !$OMP THREADPRIVATE(lcc3d)
      REAL, SAVE, ALLOCATABLE :: lcc3dcon(:,:)
      !$OMP THREADPRIVATE(lcc3dcon)
      REAL, SAVE, ALLOCATABLE :: lcc3dstra(:,:)
      !$OMP THREADPRIVATE(lcc3dstra)
      REAL, SAVE, ALLOCATABLE :: icc3dcon(:,:)
      !$OMP THREADPRIVATE(icc3dcon)
      REAL, SAVE, ALLOCATABLE :: icc3dstra(:,:)
      !$OMP THREADPRIVATE(icc3dstra)
      REAL, SAVE, ALLOCATABLE :: zfice(:,:)
      !$OMP THREADPRIVATE(zfice)
      REAL, SAVE, ALLOCATABLE :: od443aer(:) 
      !$OMP THREADPRIVATE(od443aer) 
      REAL, SAVE, ALLOCATABLE :: od550aer(:) 
      !$OMP THREADPRIVATE(od550aer) 
      REAL, SAVE, ALLOCATABLE :: dryod550aer(:) 
      !$OMP THREADPRIVATE(dryod550aer) 
      REAL, SAVE, ALLOCATABLE :: abs550aer(:) 
      !$OMP THREADPRIVATE(abs550aer) 
      REAL, SAVE, ALLOCATABLE :: od865aer(:) 
      !$OMP THREADPRIVATE(od865aer) 
      REAL, SAVE, ALLOCATABLE :: ec550aer(:,:) 
      !$OMP THREADPRIVATE(ec550aer) 
      REAL, SAVE, ALLOCATABLE :: od550lt1aer(:) 
      !$OMP THREADPRIVATE(od550lt1aer) 
      REAL, SAVE, ALLOCATABLE :: sconcso4(:) 
      !$OMP THREADPRIVATE(sconcso4) 
      REAL, SAVE, ALLOCATABLE :: sconcno3(:) 
      !$OMP THREADPRIVATE(sconcno3) 
      REAL, SAVE, ALLOCATABLE :: sconcoa(:) 
      !$OMP THREADPRIVATE(sconcoa) 
      REAL, SAVE, ALLOCATABLE :: sconcbc(:) 
      !$OMP THREADPRIVATE(sconcbc) 
      REAL, SAVE, ALLOCATABLE :: sconcss(:) 
      !$OMP THREADPRIVATE(sconcss) 
      REAL, SAVE, ALLOCATABLE :: sconcdust(:) 
      !$OMP THREADPRIVATE(sconcdust) 
      REAL, SAVE, ALLOCATABLE :: concso4(:,:) 
      !$OMP THREADPRIVATE(concso4) 
      REAL, SAVE, ALLOCATABLE :: concno3(:,:) 
      !$OMP THREADPRIVATE(concno3) 
      REAL, SAVE, ALLOCATABLE :: concoa(:,:) 
      !$OMP THREADPRIVATE(concoa) 
      REAL, SAVE, ALLOCATABLE :: concbc(:,:) 
      !$OMP THREADPRIVATE(concbc) 
      REAL, SAVE, ALLOCATABLE :: concss(:,:) 
      !$OMP THREADPRIVATE(concss) 
      REAL, SAVE, ALLOCATABLE :: concdust(:,:) 
      !$OMP THREADPRIVATE(concdust) 
      REAL, SAVE, ALLOCATABLE :: loadso4(:) 
      !$OMP THREADPRIVATE(loadso4) 
      REAL, SAVE, ALLOCATABLE :: loadoa(:) 
      !$OMP THREADPRIVATE(loadoa) 
      REAL, SAVE, ALLOCATABLE :: loadbc(:) 
      !$OMP THREADPRIVATE(loadbc) 
      REAL, SAVE, ALLOCATABLE :: loadss(:) 
      !$OMP THREADPRIVATE(loadss) 
      REAL, SAVE, ALLOCATABLE :: loaddust(:) 
      !$OMP THREADPRIVATE(loaddust) 
      REAL, SAVE, ALLOCATABLE :: loadno3(:) 
      !$OMP THREADPRIVATE(loadno3) 
      REAL, SAVE, ALLOCATABLE :: load_tmp1(:), load_tmp2(:)
      !$OMP THREADPRIVATE(load_tmp1, load_tmp2) 
      REAL, SAVE, ALLOCATABLE :: load_tmp3(:), load_tmp4(:)
      !$OMP THREADPRIVATE(load_tmp3, load_tmp4) 
      REAL, SAVE, ALLOCATABLE :: load_tmp5(:), load_tmp6(:)
      !$OMP THREADPRIVATE(load_tmp5, load_tmp6) 
      REAL, SAVE, ALLOCATABLE :: load_tmp7(:), load_tmp8(:)
      !$OMP THREADPRIVATE(load_tmp7, load_tmp8) 
      REAL, SAVE, ALLOCATABLE :: load_tmp9(:), load_tmp10(:)
      !$OMP THREADPRIVATE(load_tmp9, load_tmp10) 

!IM ajout variables CFMIP2/CMIP5
      REAL,ALLOCATABLE,SAVE :: topswad_aerop(:), solswad_aerop(:)
!$OMP THREADPRIVATE(topswad_aerop, solswad_aerop)
      REAL,ALLOCATABLE,SAVE :: topswai_aerop(:), solswai_aerop(:)
!$OMP THREADPRIVATE(topswai_aerop, solswai_aerop)
      REAL,ALLOCATABLE,SAVE :: topswad0_aerop(:), solswad0_aerop(:)
!$OMP THREADPRIVATE(topswad0_aerop, solswad0_aerop)
      REAL,ALLOCATABLE,SAVE :: topsw_aerop(:,:), topsw0_aerop(:,:)  
!$OMP THREADPRIVATE(topsw_aerop, topsw0_aerop)
      REAL,ALLOCATABLE,SAVE :: solsw_aerop(:,:), solsw0_aerop(:,:)
!$OMP THREADPRIVATE(solsw_aerop, solsw0_aerop)
      REAL,ALLOCATABLE,SAVE :: topswcf_aerop(:,:), solswcf_aerop(:,:)
!$OMP THREADPRIVATE(topswcf_aerop, solswcf_aerop)

! additional LW variables CK
      REAL,ALLOCATABLE,SAVE :: toplwad_aerop(:), sollwad_aerop(:)
!$OMP THREADPRIVATE(toplwad_aerop, sollwad_aerop)
      REAL,ALLOCATABLE,SAVE :: toplwai_aerop(:), sollwai_aerop(:)
!$OMP THREADPRIVATE(toplwai_aerop, sollwai_aerop)
      REAL,ALLOCATABLE,SAVE :: toplwad0_aerop(:), sollwad0_aerop(:)
!$OMP THREADPRIVATE(toplwad0_aerop, sollwad0_aerop)

!Ajout de celles n??cessaires au phys_output_write_mod
      REAL, SAVE, ALLOCATABLE :: tal1(:), pal1(:), pab1(:), pab2(:)
!$OMP THREADPRIVATE(tal1, pal1, pab1, pab2)
      REAL, SAVE, ALLOCATABLE :: ptstar(:), pt0(:), slp(:)
!$OMP THREADPRIVATE(ptstar, pt0, slp)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: sens, flwp, fiwp 
!$OMP THREADPRIVATE(sens, flwp, fiwp)
!!
!!         Wake variables
      REAL,ALLOCATABLE,SAVE,DIMENSION(:)            :: alp_wake
!$OMP THREADPRIVATE(alp_wake)           
!!jyg!      REAL,ALLOCATABLE,SAVE,DIMENSION(:)            :: wake_h,wake_k
      REAL,ALLOCATABLE,SAVE,DIMENSION(:)            :: wake_h
      INTEGER,ALLOCATABLE,SAVE,DIMENSION(:)         :: wake_k
!$OMP THREADPRIVATE(wake_h,wake_k)                       
    REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)            :: wake_omg
!$OMP THREADPRIVATE(wake_omg)                     
    REAL, SAVE, ALLOCATABLE,DIMENSION(:,:)          :: d_deltat_wk, d_deltaq_wk
!$OMP THREADPRIVATE(d_deltat_wk, d_deltaq_wk)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:)            :: d_s_wk, d_dens_a_wk, d_dens_wk
!$OMP THREADPRIVATE(d_s_wk, d_dens_a_wk, d_dens_wk)
    REAL, SAVE, ALLOCATABLE,DIMENSION(:,:)          :: d_deltat_wk_gw, d_deltaq_wk_gw
!$OMP THREADPRIVATE(d_deltat_wk_gw, d_deltaq_wk_gw)
    REAL, SAVE, ALLOCATABLE,DIMENSION(:,:)          :: d_deltat_vdf, d_deltaq_vdf
!$OMP THREADPRIVATE(d_deltat_vdf, d_deltaq_vdf)
!!!      REAL,ALLOCATABLE,SAVE,DIMENSION(:)          :: d_s_vdf, d_dens_vdf
!!!OMP THREADPRIVATE(d_s_vdf, d_dens_vdf)
    REAL, SAVE, ALLOCATABLE,DIMENSION(:,:)          :: d_deltat_the, d_deltaq_the
!$OMP THREADPRIVATE(d_deltat_the, d_deltaq_the)
!!!      REAL,ALLOCATABLE,SAVE,DIMENSION(:)          :: d_s_the, d_dens_the
!!!OMP THREADPRIVATE(d_s_the, d_dens_the)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)           :: d_deltat_ajs_cv, d_deltaq_ajs_cv
!$OMP THREADPRIVATE(d_deltat_ajs_cv, d_deltaq_ajs_cv)                       
!!         End of Wake variables
!!
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: bils
!$OMP THREADPRIVATE(bils)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: cdragm, cdragh
!$OMP THREADPRIVATE(cdragm, cdragh)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: cldh, cldl, cldm, cldq, cldt, qsat2m
!$OMP THREADPRIVATE(cldh, cldl, cldm, cldq, cldt, qsat2m )
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: cldhjn, cldljn, cldmjn,cldtjn
!$OMP THREADPRIVATE(cldhjn, cldljn, cldmjn, cldtjn)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: JrNt
!$OMP THREADPRIVATE(JrNt)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: dthmin, evap, fder, plcl, plfc, prw, prlw, prsw
!$OMP THREADPRIVATE(dthmin, evap, fder, plcl, plfc, prw, prlw, prsw)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: zustar, zu10m, zv10m, rh2m
!$OMP THREADPRIVATE(zustar, zu10m, zv10m, rh2m)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: s_lcl, s_pblh, s_pblt, s_therm
!$OMP THREADPRIVATE(s_lcl, s_pblh, s_pblt, s_therm)
!
!nrlmd+jyg<
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: s_pblh_x, s_pblh_w
!$OMP THREADPRIVATE(s_pblh_x, s_pblh_w)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: s_lcl_x, s_lcl_w
!$OMP THREADPRIVATE(s_lcl_x, s_lcl_w)
!>nrlmd+jyg
!
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: slab_wfbils
!$OMP THREADPRIVATE(slab_wfbils)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: tpot, tpote, ue, uq, uwat, ve, vq, vwat, zxffonte
!$OMP THREADPRIVATE(tpot, tpote, ue, uq, uwat, ve, vq, vwat, zxffonte)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: zxfqcalving
!$OMP THREADPRIVATE(zxfqcalving)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: zxfluxlat, zxtsol, snow_lsc, zxfqfonte
!$OMP THREADPRIVATE(zxfluxlat, zxtsol, snow_lsc, zxfqfonte)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: zxrunofflic
!$OMP THREADPRIVATE(zxrunofflic)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: zxqsurf, rain_lsc, rain_num
!$OMP THREADPRIVATE(zxqsurf, rain_lsc, rain_num)
!
!jyg+nrlmd<
!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                          c
!       Declarations liees a la couche limite differentiee w-x             c
!                                                                          c
!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: sens_x, sens_w
!$OMP THREADPRIVATE(sens_x, sens_w)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: zxfluxlat_x, zxfluxlat_w
!$OMP THREADPRIVATE(zxfluxlat_x, zxfluxlat_w)
!jyg<
!!! Entrees supplementaires couche-limite
!!      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: t_x, t_w
!!!$OMP THREADPRIVATE(t_x, t_w)
!!      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: q_x, q_w
!!!$OMP THREADPRIVATE(q_x, q_w)
!>jyg
!!! Sorties ferret
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: dtvdf_x, dtvdf_w
!$OMP THREADPRIVATE(dtvdf_x, dtvdf_w)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: dqvdf_x, dqvdf_w
!$OMP THREADPRIVATE(dqvdf_x, dqvdf_w)
! Variables supplementaires dans physiq.F relative au splitting de la surface
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:) :: pbl_tke_input
!$OMP THREADPRIVATE(pbl_tke_input)
! Entree supplementaire Thermiques :
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: t_therm, q_therm, u_therm, v_therm
!$OMP THREADPRIVATE(t_therm, q_therm, u_therm, v_therm)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: cdragh_x, cdragh_w
!$OMP THREADPRIVATE(cdragh_x, cdragh_w)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: cdragm_x, cdragm_w
!$OMP THREADPRIVATE(cdragm_x, cdragm_w)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: kh, kh_x, kh_w
!$OMP THREADPRIVATE(kh, kh_x, kh_w)
!!!
!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      LOGICAL, SAVE, ALLOCATABLE :: ptconv(:,:)
      !$OMP THREADPRIVATE(ptconv)
      REAL, SAVE, ALLOCATABLE :: ratqsc(:,:)
      !$OMP THREADPRIVATE(ratqsc)
!>jyg+nrlmd
  !
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: wbeff, convoccur, zmax_th, zq2m, zt2m
!$OMP THREADPRIVATE(wbeff, convoccur, zmax_th, zq2m, zt2m)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: zt2m_min_mon, zt2m_max_mon
!$OMP THREADPRIVATE(zt2m_min_mon, zt2m_max_mon)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: t2m_min_mon, t2m_max_mon
!$OMP THREADPRIVATE(t2m_min_mon, t2m_max_mon)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: weak_inversion
!$OMP THREADPRIVATE(weak_inversion)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: alp_bl_conv,alp_bl_det
!$OMP THREADPRIVATE(alp_bl_conv,alp_bl_det)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: alp_bl_fluct_m,alp_bl_fluct_tke
!$OMP THREADPRIVATE(alp_bl_fluct_m,alp_bl_fluct_tke)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: alp_bl_stat, n2, s2
!$OMP THREADPRIVATE(alp_bl_stat, n2, s2)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: proba_notrig, random_notrig
!$OMP THREADPRIVATE(proba_notrig, random_notrig)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: cv_gen
!$OMP THREADPRIVATE(cv_gen)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: fsolsw, wfbils, wfbilo
!$OMP THREADPRIVATE(fsolsw, wfbils, wfbilo)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: wfevap, wfrain, wfsnow
!$OMP THREADPRIVATE(wfevap,wfrain,wfsnow)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: t2m, fluxlat, fsollw,evap_pot
!$OMP THREADPRIVATE(t2m, fluxlat, fsollw,evap_pot)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: dnwd0, omega
!$OMP THREADPRIVATE(dnwd0, omega)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: epmax_diag ! epmax_cape
!$OMP THREADPRIVATE(epmax_diag)
!
!  Deep convective variables used in phytrac
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: ep ! epmax_cape
!$OMP THREADPRIVATE(ep)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: dnwd, upwd
!$OMP THREADPRIVATE(dnwd, upwd)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: pmflxr, pmflxs
!$OMP THREADPRIVATE(pmflxr, pmflxs)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: wdtrainA, wdtrainS, wdtrainM
!$OMP THREADPRIVATE(wdtrainA, wdtrainS, wdtrainM)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: da, mp
!$OMP THREADPRIVATE(da, mp)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: wght_cvfd
!$OMP THREADPRIVATE(wght_cvfd)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:):: phi, phi2, elij, epmlmMm
!$OMP THREADPRIVATE(phi, phi2, elij, epmlmMm)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: d1a, dam
!$OMP THREADPRIVATE(d1a, dam)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: ev
!$OMP THREADPRIVATE(ev)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: qtaa
!$OMP THREADPRIVATE(qtaa)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: clw
!$OMP THREADPRIVATE(clw)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: eplaMm
!$OMP THREADPRIVATE(eplaMm)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:):: sij
!$OMP THREADPRIVATE(sij)
!
!      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: coefh, coefm, lambda_th
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: lambda_th
!$OMP THREADPRIVATE(lambda_th)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: cldemi, cldfra, cldtau, fiwc, fl, re, flwc
!$OMP THREADPRIVATE(cldemi, cldfra, cldtau, fiwc, fl, re, flwc)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: ref_liq, ref_ice, theta, zphi
!$OMP THREADPRIVATE(ref_liq, ref_ice, theta, zphi)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: ref_liq_pi, ref_ice_pi
!$OMP THREADPRIVATE(ref_liq_pi, ref_ice_pi)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: zx_rh
!$OMP THREADPRIVATE(zx_rh)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: prfl, psfl, fraca
!$OMP THREADPRIVATE(prfl, psfl, fraca)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: Vprecip, zw2
!$OMP THREADPRIVATE(Vprecip, zw2)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: pmfd, pmfu
!$OMP THREADPRIVATE(pmfd, pmfu)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:) :: fluxt, fluxu, fluxv
!$OMP THREADPRIVATE(fluxt, fluxu, fluxv)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:) :: uwriteSTD, vwriteSTD, wwriteSTD
!$OMP THREADPRIVATE(uwriteSTD, vwriteSTD, wwriteSTD)
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:) :: phiwriteSTD, qwriteSTD, twriteSTD, rhwriteSTD
!$OMP THREADPRIVATE(phiwriteSTD, qwriteSTD, twriteSTD, rhwriteSTD)

! ug et d'autres encore:
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: beta_prec
!$OMP THREADPRIVATE(beta_prec)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: rneb,rnebjn,rneblsvol
!$OMP THREADPRIVATE(rneb,rnebjn,rneblsvol)

! variables de sorties MM
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: zxsnow,snowhgt,qsnow,to_ice
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: sissnow,runoff,albsol3_lic
!$OMP THREADPRIVATE(zxsnow,snowhgt,qsnow,to_ice)
!$OMP THREADPRIVATE(sissnow,runoff,albsol3_lic)

      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: p_tropopause, z_tropopause, t_tropopause
!$OMP THREADPRIVATE(p_tropopause, z_tropopause, t_tropopause)

      INTEGER,ALLOCATABLE,SAVE,DIMENSION(:,:) :: zn2mout
!$OMP THREADPRIVATE(zn2mout)

#ifdef CPP_StratAer
!
! variables for stratospheric aerosol
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: R2SO4
!$OMP THREADPRIVATE(R2SO4)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: DENSO4
!$OMP THREADPRIVATE(DENSO4)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: f_r_wet
!$OMP THREADPRIVATE(f_r_wet)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: decfluxaer
!$OMP THREADPRIVATE(decfluxaer)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: mdw
!$OMP THREADPRIVATE(mdw)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: OCS_lifetime
!$OMP THREADPRIVATE(OCS_lifetime)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: SO2_lifetime
!$OMP THREADPRIVATE(SO2_lifetime)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: alpha_bin
!$OMP THREADPRIVATE(alpha_bin)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: piz_bin
!$OMP THREADPRIVATE(piz_bin)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: cg_bin
!$OMP THREADPRIVATE(cg_bin)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tau_strat_550
!$OMP THREADPRIVATE(tau_strat_550)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tau_strat_1020
!$OMP THREADPRIVATE(tau_strat_1020)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tausum_strat
!$OMP THREADPRIVATE(tausum_strat)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: surf_PM25_sulf
!$OMP THREADPRIVATE(surf_PM25_sulf)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: vsed_aer
!$OMP THREADPRIVATE(vsed_aer)
!
!---3D budget variables
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: budg_3D_nucl
!$OMP THREADPRIVATE(budg_3D_nucl)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: budg_3D_cond_evap
!$OMP THREADPRIVATE(budg_3D_cond_evap)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: budg_3D_ocs_to_so2
!$OMP THREADPRIVATE(budg_3D_ocs_to_so2)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: budg_3D_so2_to_h2so4
!$OMP THREADPRIVATE(budg_3D_so2_to_h2so4)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: budg_3D_backgr_ocs
!$OMP THREADPRIVATE(budg_3D_backgr_ocs)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: budg_3D_backgr_so2
!$OMP THREADPRIVATE(budg_3D_backgr_so2)
!
!---2D budget variables
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_dep_dry_ocs
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_dep_dry_so2
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_dep_dry_h2so4
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_dep_dry_part
!$OMP THREADPRIVATE(budg_dep_dry_ocs,budg_dep_dry_so2,budg_dep_dry_h2so4,budg_dep_dry_part)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_dep_wet_ocs
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_dep_wet_so2
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_dep_wet_h2so4
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_dep_wet_part
!$OMP THREADPRIVATE(budg_dep_wet_ocs,budg_dep_wet_so2,budg_dep_wet_h2so4,budg_dep_wet_part)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_emi_ocs
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_emi_so2
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_emi_h2so4
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_emi_part
!$OMP THREADPRIVATE(budg_emi_ocs,budg_emi_so2,budg_emi_h2so4,budg_emi_part)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_ocs_to_so2
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_so2_to_h2so4
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_h2so4_to_part
!$OMP THREADPRIVATE(budg_ocs_to_so2,budg_so2_to_h2so4,budg_h2so4_to_part)
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: budg_sed_part
!$OMP THREADPRIVATE(budg_sed_part)
#endif

CONTAINS

!======================================================================
SUBROUTINE phys_local_var_init
USE dimphy
USE infotrac_phy, ONLY : nbtr
USE aero_mod
USE indice_sol_mod
USE phys_output_var_mod
USE phys_state_var_mod

IMPLICIT NONE
      ALLOCATE(t_seri(klon,klev),q_seri(klon,klev),ql_seri(klon,klev),qs_seri(klon,klev))
      ALLOCATE(u_seri(klon,klev),v_seri(klon,klev))
      ALLOCATE(l_mixmin(klon,klev+1,nbsrf), l_mix(klon,klev+1,nbsrf))
      l_mix(:,:,:)=0. ; l_mixmin(:,:,:)=0. ! doit etre initialse car pas toujours remplis

      ALLOCATE(tr_seri(klon,klev,nbtr))
      ALLOCATE(d_t_dyn(klon,klev),d_q_dyn(klon,klev))
      ALLOCATE(d_ql_dyn(klon,klev),d_qs_dyn(klon,klev))
      ALLOCATE(d_q_dyn2d(klon),d_ql_dyn2d(klon),d_qs_dyn2d(klon))
      ALLOCATE(d_u_dyn(klon,klev),d_v_dyn(klon,klev))
      ALLOCATE(d_tr_dyn(klon,klev,nbtr))                   !RomP
      ALLOCATE(d_t_con(klon,klev),d_q_con(klon,klev))
      ALLOCATE(d_u_con(klon,klev),d_v_con(klon,klev))
      ALLOCATE(d_t_wake(klon,klev),d_q_wake(klon,klev))
      ALLOCATE(d_t_lsc(klon,klev),d_q_lsc(klon,klev))
      ALLOCATE(d_t_lwr(klon,klev),d_t_lw0(klon,klev))
      ALLOCATE(d_t_swr(klon,klev),d_t_sw0(klon,klev))
      ALLOCATE(d_ql_lsc(klon,klev),d_qi_lsc(klon,klev))
      ALLOCATE(d_t_ajsb(klon,klev),d_q_ajsb(klon,klev))
      ALLOCATE(d_t_ajs(klon,klev),d_q_ajs(klon,klev))
!nrlmd<
      ALLOCATE(d_t_ajs_w(klon,klev),d_q_ajs_w(klon,klev))
      ALLOCATE(d_t_ajs_x(klon,klev),d_q_ajs_x(klon,klev))
!>nrlmd
      ALLOCATE(d_u_ajs(klon,klev),d_v_ajs(klon,klev))
      ALLOCATE(d_t_eva(klon,klev),d_q_eva(klon,klev))
      ALLOCATE(d_ql_eva(klon,klev),d_qi_eva(klon,klev))
      ALLOCATE(d_t_lscst(klon,klev),d_q_lscst(klon,klev))
      ALLOCATE(d_t_lscth(klon,klev),d_q_lscth(klon,klev))
      ALLOCATE(plul_st(klon),plul_th(klon))
      ALLOCATE(d_t_vdf(klon,klev),d_q_vdf(klon,klev),d_t_diss(klon,klev))

      ALLOCATE(d_t_vdf_w(klon,klev),d_q_vdf_w(klon,klev))
      ALLOCATE(d_t_vdf_x(klon,klev),d_q_vdf_x(klon,klev))

      ALLOCATE(d_u_vdf(klon,klev),d_v_vdf(klon,klev))
      ALLOCATE(d_t_oli(klon,klev),d_t_oro(klon,klev))
      ALLOCATE(d_u_oli(klon,klev),d_v_oli(klon,klev))
      ALLOCATE(d_u_oro(klon,klev),d_v_oro(klon,klev))
      ALLOCATE(d_u_oro_gw(klon,klev),d_v_oro_gw(klon,klev))
      ALLOCATE(d_t_oro_gw(klon,klev))
      ALLOCATE(d_t_lif(klon,klev),d_t_ec(klon,klev))
      ALLOCATE(d_u_lif(klon,klev),d_v_lif(klon,klev))
      ALLOCATE(d_ts(klon,nbsrf), d_tr(klon,klev,nbtr))
! Special RRTM
      ALLOCATE(ZLWFT0_i(klon,klev+1),ZSWFT0_i(klon,klev+1),ZFLDN0(klon,klev+1))
      ZFLDN0= 0.
      ALLOCATE(ZFLUP0(klon,klev+1),ZFSDN0(klon,klev+1),ZFSUP0(klon,klev+1))
!
      ALLOCATE(topswad_aero(klon), solswad_aero(klon))
      ALLOCATE(topswai_aero(klon), solswai_aero(klon))
      ALLOCATE(topswad0_aero(klon), solswad0_aero(klon))
      ALLOCATE(toplwad_aero(klon), sollwad_aero(klon))
      ALLOCATE(toplwai_aero(klon), sollwai_aero(klon))
      ALLOCATE(toplwad0_aero(klon), sollwad0_aero(klon))
      ALLOCATE(topsw_aero(klon,naero_grp), solsw_aero(klon,naero_grp))
      ALLOCATE(topsw0_aero(klon,naero_grp), solsw0_aero(klon,naero_grp))
      ALLOCATE(topswcf_aero(klon,3), solswcf_aero(klon,3))
      ALLOCATE(du_gwd_hines(klon,klev),dv_gwd_hines(klon,klev))
      ALLOCATE(dv_gwd_rando(klon,klev),dv_gwd_front(klon,klev))
      ALLOCATE(east_gwstress(klon,klev),west_gwstress(klon,klev))
      east_gwstress(:,:)=0 !ym missing init
      west_gwstress(:,:)=0 !ym missing init
      ALLOCATE(d_t_hin(klon,klev))
      ALLOCATE(d_q_ch4(klon,klev))
      ALLOCATE(stratomask(klon,klev))
      ALLOCATE(tausum_aero(klon,nwave,naero_tot))
      ALLOCATE(drytausum_aero(klon,naero_tot))
      ALLOCATE(tau3d_aero(klon,klev,nwave,naero_tot)) 
      ALLOCATE(scdnc(klon, klev))
      ALLOCATE(dNovrN(klon))
      ALLOCATE(cldncl(klon))
      ALLOCATE(reffclwtop(klon))
      ALLOCATE(lcc(klon))
      ALLOCATE(reffclws(klon, klev))
      ALLOCATE(reffclwc(klon, klev))
      ALLOCATE(cldnvi(klon))
      ALLOCATE(lcc3d(klon, klev))
      ALLOCATE(lcc3dcon(klon, klev))
      ALLOCATE(lcc3dstra(klon, klev))
      ALLOCATE(icc3dcon(klon, klev))
      ALLOCATE(icc3dstra(klon, klev))
      ALLOCATE(zfice(klon, klev))
      ALLOCATE(od443aer(klon))
      ALLOCATE(od550aer(klon))
      ALLOCATE(od865aer(klon))
      ALLOCATE(dryod550aer(klon))
      dryod550aer(:) = 0.
      ALLOCATE(abs550aer(klon))
      abs550aer(:) = 0.
      ALLOCATE(ec550aer(klon,klev))
      ALLOCATE(od550lt1aer(klon))
      ALLOCATE(sconcso4(klon))
      ALLOCATE(sconcno3(klon))
      ALLOCATE(sconcoa(klon))
      ALLOCATE(sconcbc(klon))
      ALLOCATE(sconcss(klon))
      ALLOCATE(sconcdust(klon))
      ALLOCATE(concso4(klon,klev))
      ALLOCATE(concno3(klon,klev))
      ALLOCATE(concoa(klon,klev))
      ALLOCATE(concbc(klon,klev))
      ALLOCATE(concss(klon,klev))
      ALLOCATE(concdust(klon,klev))
      ALLOCATE(loadso4(klon))
      ALLOCATE(loadoa(klon))
      ALLOCATE(loadbc(klon))
      ALLOCATE(loadss(klon))
      ALLOCATE(loaddust(klon))
      ALLOCATE(loadno3(klon))
      ALLOCATE(load_tmp1(klon))
      ALLOCATE(load_tmp2(klon))
      ALLOCATE(load_tmp3(klon))
      ALLOCATE(load_tmp4(klon))
      ALLOCATE(load_tmp5(klon))
      ALLOCATE(load_tmp6(klon))
      ALLOCATE(load_tmp7(klon))
      ALLOCATE(load_tmp8(klon))
      ALLOCATE(load_tmp9(klon))
      ALLOCATE(load_tmp10(klon))

!IM ajout variables CFMIP2/CMIP5
      ALLOCATE(topswad_aerop(klon), solswad_aerop(klon))
      ALLOCATE(topswai_aerop(klon), solswai_aerop(klon))
      ALLOCATE(topswad0_aerop(klon), solswad0_aerop(klon))
      ALLOCATE(topsw_aerop(klon,naero_grp), topsw0_aerop(klon,naero_grp))
      ALLOCATE(solsw_aerop(klon,naero_grp), solsw0_aerop(klon,naero_grp))
      ALLOCATE(topswcf_aerop(klon,naero_grp), solswcf_aerop(klon,naero_grp))

! additional LW variables CK
      ALLOCATE(toplwad_aerop(klon), sollwad_aerop(klon))
      ALLOCATE(toplwai_aerop(klon), sollwai_aerop(klon))
      ALLOCATE(toplwad0_aerop(klon), sollwad0_aerop(klon))

! FH Ajout de celles necessaires au phys_output_write_mod

      ALLOCATE(tal1(klon), pal1(klon), pab1(klon), pab2(klon))
      ALLOCATE(ptstar(klon),pt0(klon),slp(klon))
!!
!!          Wake variables
      ALLOCATE(alp_wake(klon))
      ALLOCATE(wake_h(klon),wake_k(klon))
      ALLOCATE(wake_omg(klon, klev))
      ALLOCATE(d_deltat_wk(klon, klev), d_deltaq_wk(klon, klev))
      ALLOCATE(d_s_wk(klon), d_dens_a_wk(klon), d_dens_wk(klon))
      ALLOCATE(d_deltat_wk_gw(klon, klev), d_deltaq_wk_gw(klon, klev))
      ALLOCATE(d_deltat_vdf(klon, klev), d_deltaq_vdf(klon, klev))
!!      ALLOCATE( d_s_vdf(klon), d_dens_vdf(klon))
      ALLOCATE(d_deltat_the(klon, klev), d_deltaq_the(klon, klev))
!!      ALLOCATE( d_s_the(klon), d_dens_the(klon))
      ALLOCATE(d_deltat_ajs_cv(klon, klev), d_deltaq_ajs_cv(klon, klev))
!!         End of wake variables
!!
      ALLOCATE(bils(klon))
      ALLOCATE(cdragm(klon), cdragh(klon), cldh(klon), cldl(klon))
      ALLOCATE(cldm(klon), cldq(klon), cldt(klon), qsat2m(klon))
      ALLOCATE(cldhjn(klon), cldljn(klon), cldmjn(klon), cldtjn(klon))
      ALLOCATE(JrNt(klon))
      ALLOCATE(dthmin(klon), evap(klon), fder(klon), plcl(klon), plfc(klon))
      ALLOCATE(prw(klon), prlw(klon), prsw(klon), zustar(klon), zu10m(klon), zv10m(klon), rh2m(klon))
      ALLOCATE(s_lcl(klon))
      ALLOCATE(s_pblh(klon), s_pblt(klon), s_therm(klon))
!
!nrlmd+jyg<
      ALLOCATE(s_pblh_x(klon), s_pblh_w(klon))
      ALLOCATE(s_lcl_x(klon), s_lcl_w(klon))
!>nrlmd+jyg
!
      ALLOCATE(slab_wfbils(klon), tpot(klon), tpote(klon), ue(klon))
      ALLOCATE(uq(klon), ve(klon), vq(klon), zxffonte(klon))
      ALLOCATE(uwat(klon), vwat(klon))
      ALLOCATE(zxfqcalving(klon), zxfluxlat(klon))
      ALLOCATE(zxtsol(klon), snow_lsc(klon), zxfqfonte(klon), zxqsurf(klon))
      ALLOCATE(zxrunofflic(klon))
      ALLOCATE(rain_lsc(klon))
      ALLOCATE(rain_num(klon))
!
      ALLOCATE(sens_x(klon), sens_w(klon))
      ALLOCATE(zxfluxlat_x(klon), zxfluxlat_w(klon))
!jyg<
!!      ALLOCATE(t_x(klon,klev), t_w(klon,klev))
!!      ALLOCATE(q_x(klon,klev), q_w(klon,klev))
!>jyg
      ALLOCATE(dtvdf_x(klon,klev), dtvdf_w(klon,klev))
      dtvdf_x = 0 ; dtvdf_w=0 ;   !ym missing init
      ALLOCATE(dqvdf_x(klon,klev), dqvdf_w(klon,klev))
      dqvdf_x = 0 ; dqvdf_w=0 ;   !ym missing init
      ALLOCATE(pbl_tke_input(klon,klev+1,nbsrf))
      ALLOCATE(t_therm(klon,klev), q_therm(klon,klev),u_therm(klon,klev), v_therm(klon,klev))
      ALLOCATE(cdragh_x(klon), cdragh_w(klon))
      ALLOCATE(cdragm_x(klon), cdragm_w(klon))
      ALLOCATE(kh(klon), kh_x(klon), kh_w(klon))
!
      ALLOCATE(ptconv(klon,klev))
      ALLOCATE(ratqsc(klon,klev))
!
      ALLOCATE(wbeff(klon), convoccur(klon), zmax_th(klon))
      ALLOCATE(zq2m(klon), zt2m(klon), weak_inversion(klon))
      ALLOCATE(zt2m_min_mon(klon), zt2m_max_mon(klon))
      ALLOCATE(t2m_min_mon(klon), t2m_max_mon(klon))
      ALLOCATE(sens(klon), flwp(klon), fiwp(klon))
      ALLOCATE(alp_bl_conv(klon), alp_bl_det(klon))
      alp_bl_conv(:)=0 ;  alp_bl_det(:)=0
      ALLOCATE(alp_bl_fluct_m(klon), alp_bl_fluct_tke(klon))
      alp_bl_fluct_m(:)=0 ; alp_bl_fluct_tke(:)= 0.
      ALLOCATE(alp_bl_stat(klon), n2(klon), s2(klon))
      alp_bl_stat(:)=0
      ALLOCATE(proba_notrig(klon), random_notrig(klon))
      ALLOCATE(cv_gen(klon))

      ALLOCATE(dnwd0(klon, klev))
!      ALLOCATE(upwd(klon, klev), omega(klon, klev), coefh(klon, klev))
      ALLOCATE(omega(klon, klev))
      ALLOCATE(epmax_diag(klon)) ! epmax_cape
!      ALLOCATE(coefm(klon, klev), lambda_th(klon, klev), cldemi(klon, klev))
      ALLOCATE(lambda_th(klon, klev), cldemi(klon, klev))
      ALLOCATE(cldfra(klon, klev), cldtau(klon, klev), fiwc(klon, klev))
      ALLOCATE(fl(klon, klev), re(klon, klev), flwc(klon, klev))
      ALLOCATE(ref_liq(klon, klev), ref_ice(klon, klev), theta(klon, klev))
      ALLOCATE(ref_liq_pi(klon, klev), ref_ice_pi(klon, klev))
      ALLOCATE(zphi(klon, klev), zx_rh(klon, klev))
      ALLOCATE(pmfd(klon, klev), pmfu(klon, klev))

      ALLOCATE(t2m(klon, nbsrf), fluxlat(klon, nbsrf))
      ALLOCATE(fsollw(klon, nbsrf))
      ALLOCATE(fsolsw(klon, nbsrf), wfbils(klon, nbsrf), wfbilo(klon, nbsrf))
      ALLOCATE(wfevap(klon, nbsrf), wfrain(klon,nbsrf), wfsnow(klon, nbsrf))
      ALLOCATE(evap_pot(klon, nbsrf))
!
!  Deep convective variables used in phytrac
      ALLOCATE(pmflxr(klon, klev+1), pmflxs(klon, klev+1))
      ALLOCATE(wdtrainA(klon,klev),wdtrainS(klon,klev),wdtrainM(klon,klev))
      ALLOCATE(dnwd(klon, klev), upwd(klon, klev) )
      ALLOCATE(ep(klon,klev))                          ! epmax_cape
      ALLOCATE(da(klon,klev), mp(klon,klev) )
      ALLOCATE(phi(klon,klev,klev) )
      ALLOCATE(wght_cvfd(klon,klev) )
      ALLOCATE(phi2(klon,klev,klev) )
      ALLOCATE(d1a(klon,klev), dam(klon,klev))
      ALLOCATE(ev(klon,klev) )
      ALLOCATE(elij(klon,klev,klev) )
      ALLOCATE(qtaa(klon,klev) )
      ALLOCATE(clw(klon,klev) )
      ALLOCATE(epmlmMm(klon,klev,klev), eplaMm(klon,klev) )
      ALLOCATE(sij(klon,klev,klev) )

      ALLOCATE(prfl(klon, klev+1))
      ALLOCATE(psfl(klon, klev+1), fraca(klon, klev+1), Vprecip(klon, klev+1))
      ALLOCATE(zw2(klon, klev+1))

      ALLOCATE(fluxu(klon, klev, nbsrf), fluxv(klon, klev, nbsrf))
      ALLOCATE(fluxt(klon, klev, nbsrf))

      ALLOCATE(uwriteSTD(klon,nlevSTD,nfiles), vwriteSTD(klon,nlevSTD,nfiles))
      ALLOCATE(wwriteSTD(klon,nlevSTD,nfiles), phiwriteSTD(klon,nlevSTD,nfiles))
      ALLOCATE(qwriteSTD(klon,nlevSTD,nfiles), twriteSTD(klon,nlevSTD,nfiles))
      ALLOCATE(rhwriteSTD(klon,nlevSTD,nfiles))

! ug et d'autres encore:
      ALLOCATE(beta_prec(klon,klev))
      ALLOCATE(rneb(klon,klev),rnebjn(klon,klev),rneblsvol(klon,klev))


      ALLOCATE (zxsnow(klon),snowhgt(klon),qsnow(klon),to_ice(klon))
      ALLOCATE (sissnow(klon),runoff(klon),albsol3_lic(klon))

      ALLOCATE (p_tropopause(klon))
      ALLOCATE (z_tropopause(klon))
      ALLOCATE (t_tropopause(klon))

      ALLOCATE(zn2mout(klon,6))

#ifdef CPP_StratAer
      ALLOCATE (R2SO4(klon,klev))
      ALLOCATE (DENSO4(klon,klev))
      ALLOCATE (f_r_wet(klon,klev))
      ALLOCATE (decfluxaer(klon,nbtr))
      ALLOCATE (mdw(nbtr))
      ALLOCATE (budg_3D_nucl(klon,klev))
      ALLOCATE (budg_3D_cond_evap(klon,klev))
      ALLOCATE (budg_3D_ocs_to_so2(klon,klev))
      ALLOCATE (budg_3D_so2_to_h2so4(klon,klev))
      ALLOCATE (budg_3D_backgr_ocs(klon,klev))
      ALLOCATE (budg_3D_backgr_so2(klon,klev))
      ALLOCATE (OCS_lifetime(klon,klev))
      ALLOCATE (SO2_lifetime(klon,klev))
      ALLOCATE (alpha_bin(nbands_sw_rrtm+nbands_lw_rrtm+nwave,nbtr))
      ALLOCATE (piz_bin(nbands_sw_rrtm+nbands_lw_rrtm+nwave,nbtr))
      ALLOCATE (cg_bin(nbands_sw_rrtm+nbands_lw_rrtm+nwave,nbtr))
      ALLOCATE (tau_strat_550(klon,klev))
      ALLOCATE (tau_strat_1020(klon,klev))
      ALLOCATE (tausum_strat(klon,3))
      ALLOCATE (budg_dep_dry_ocs(klon))
      ALLOCATE (budg_dep_wet_ocs(klon))
      ALLOCATE (budg_dep_dry_so2(klon))
      ALLOCATE (budg_dep_wet_so2(klon))
      ALLOCATE (budg_dep_dry_h2so4(klon))
      ALLOCATE (budg_dep_wet_h2so4(klon))
      ALLOCATE (budg_dep_dry_part(klon))
      ALLOCATE (budg_dep_wet_part(klon))
      ALLOCATE (budg_emi_ocs(klon))
      ALLOCATE (budg_emi_so2(klon))
      ALLOCATE (budg_emi_h2so4(klon))
      ALLOCATE (budg_emi_part(klon))
      ALLOCATE (budg_ocs_to_so2(klon))
      ALLOCATE (budg_so2_to_h2so4(klon))
      ALLOCATE (budg_h2so4_to_part(klon))
      ALLOCATE (budg_sed_part(klon))
      ALLOCATE (surf_PM25_sulf(klon))
      ALLOCATE (vsed_aer(klon,klev))
#endif

END SUBROUTINE phys_local_var_init

!======================================================================
SUBROUTINE phys_local_var_end
USE dimphy
USE indice_sol_mod
IMPLICIT NONE
      DEALLOCATE(t_seri,q_seri,ql_seri,qs_seri)
      DEALLOCATE(u_seri,v_seri)
      DEALLOCATE(l_mixmin,l_mix)

      DEALLOCATE(tr_seri)
      DEALLOCATE(d_t_dyn,d_q_dyn)
      DEALLOCATE(d_ql_dyn,d_qs_dyn)
      DEALLOCATE(d_q_dyn2d,d_ql_dyn2d,d_qs_dyn2d)
      DEALLOCATE(d_u_dyn,d_v_dyn)
      DEALLOCATE(d_tr_dyn)                      !RomP
      DEALLOCATE(d_t_con,d_q_con)
      DEALLOCATE(d_u_con,d_v_con)
      DEALLOCATE(d_t_wake,d_q_wake)
      DEALLOCATE(d_t_lsc,d_q_lsc)
      DEALLOCATE(d_t_lwr,d_t_lw0)
      DEALLOCATE(d_t_swr,d_t_sw0)
      DEALLOCATE(d_ql_lsc,d_qi_lsc)
      DEALLOCATE(d_t_ajsb,d_q_ajsb)
      DEALLOCATE(d_t_ajs,d_q_ajs)
!nrlmd<
      DEALLOCATE(d_t_ajs_w,d_q_ajs_w)
      DEALLOCATE(d_t_ajs_x,d_q_ajs_x)
!>nrlmd
      DEALLOCATE(d_u_ajs,d_v_ajs)
      DEALLOCATE(d_t_eva,d_q_eva)
      DEALLOCATE(d_ql_eva,d_qi_eva)
      DEALLOCATE(d_t_lscst,d_q_lscst)
      DEALLOCATE(d_t_lscth,d_q_lscth)
      DEALLOCATE(plul_st,plul_th)
      DEALLOCATE(d_t_vdf,d_q_vdf,d_t_diss)
      DEALLOCATE(d_u_vdf,d_v_vdf)
      DEALLOCATE(d_t_oli,d_t_oro)
      DEALLOCATE(d_u_oli,d_v_oli)
      DEALLOCATE(d_u_oro,d_v_oro)
      DEALLOCATE(d_t_oro_gw)
      DEALLOCATE(d_u_oro_gw,d_v_oro_gw)
      DEALLOCATE(d_t_lif,d_t_ec)
      DEALLOCATE(d_u_lif,d_v_lif)
      DEALLOCATE(d_ts, d_tr)
      DEALLOCATE(topswad_aero,solswad_aero)
      DEALLOCATE(topswai_aero,solswai_aero)
      DEALLOCATE(topswad0_aero,solswad0_aero)
      DEALLOCATE(toplwad_aero,sollwad_aero)
      DEALLOCATE(toplwai_aero,sollwai_aero)
      DEALLOCATE(toplwad0_aero,sollwad0_aero)
      DEALLOCATE(topsw_aero,solsw_aero)
      DEALLOCATE(topsw0_aero,solsw0_aero)
      DEALLOCATE(topswcf_aero,solswcf_aero)
      DEALLOCATE(stratomask) 
      DEALLOCATE(tausum_aero) 
      DEALLOCATE(drytausum_aero) 
      DEALLOCATE(tau3d_aero) 
      DEALLOCATE(scdnc)
      DEALLOCATE(dNovrN)
      DEALLOCATE(cldncl)
      DEALLOCATE(reffclwtop)
      DEALLOCATE(lcc)
      DEALLOCATE(reffclws)
      DEALLOCATE(reffclwc)
      DEALLOCATE(cldnvi)
      DEALLOCATE(lcc3d)
      DEALLOCATE(lcc3dcon)
      DEALLOCATE(lcc3dstra)
      DEALLOCATE(icc3dcon)
      DEALLOCATE(icc3dstra)
      DEALLOCATE(zfice)
      DEALLOCATE(od443aer)
      DEALLOCATE(od550aer)
      DEALLOCATE(od865aer)
      DEALLOCATE(dryod550aer)
      DEALLOCATE(abs550aer)
      DEALLOCATE(ec550aer)
      DEALLOCATE(od550lt1aer)
      DEALLOCATE(sconcso4) 
      DEALLOCATE(sconcno3) 
      DEALLOCATE(sconcoa) 
      DEALLOCATE(sconcbc) 
      DEALLOCATE(sconcss) 
      DEALLOCATE(sconcdust) 
      DEALLOCATE(concso4) 
      DEALLOCATE(concno3) 
      DEALLOCATE(concoa) 
      DEALLOCATE(concbc) 
      DEALLOCATE(concss) 
      DEALLOCATE(concdust) 
      DEALLOCATE(loadso4) 
      DEALLOCATE(loadoa) 
      DEALLOCATE(loadbc) 
      DEALLOCATE(loadss) 
      DEALLOCATE(loaddust) 
      DEALLOCATE(loadno3) 
      DEALLOCATE(load_tmp1)
      DEALLOCATE(load_tmp2)
      DEALLOCATE(load_tmp3)
      DEALLOCATE(load_tmp4)
      DEALLOCATE(load_tmp5)
      DEALLOCATE(load_tmp6)
      DEALLOCATE(load_tmp7)
      DEALLOCATE(load_tmp8)
      DEALLOCATE(load_tmp9)
      DEALLOCATE(load_tmp10)
      DEALLOCATE(du_gwd_hines,dv_gwd_hines,d_t_hin)
      DEALLOCATE(d_q_ch4)
      DEALLOCATE(dv_gwd_rando,dv_gwd_front)
      DEALLOCATE(east_gwstress,west_gwstress)

!IM ajout variables CFMIP2/CMIP5
      DEALLOCATE(topswad_aerop, solswad_aerop)
      DEALLOCATE(topswai_aerop, solswai_aerop)
      DEALLOCATE(topswad0_aerop, solswad0_aerop)
      DEALLOCATE(topsw_aerop, topsw0_aerop)
      DEALLOCATE(solsw_aerop, solsw0_aerop)
      DEALLOCATE(topswcf_aerop, solswcf_aerop)

!CK LW diagnostics
      DEALLOCATE(toplwad_aerop, sollwad_aerop)
      DEALLOCATE(toplwai_aerop, sollwai_aerop)
      DEALLOCATE(toplwad0_aerop, sollwad0_aerop)

! FH Ajout de celles necessaires au phys_output_write_mod
      DEALLOCATE(tal1, pal1, pab1, pab2)
      DEALLOCATE(ptstar, pt0, slp)
!
      DEALLOCATE(alp_wake)
      DEALLOCATE(wake_h,wake_k)
      DEALLOCATE(wake_omg)
      DEALLOCATE(d_deltat_wk, d_deltaq_wk)
      DEALLOCATE(d_s_wk, d_dens_a_wk, d_dens_wk)
      DEALLOCATE(d_deltat_wk_gw, d_deltaq_wk_gw)
      DEALLOCATE(d_deltat_vdf, d_deltaq_vdf)
!!      DEALLOCATE( d_s_vdf, d_dens_vdf)
      DEALLOCATE(d_deltat_the, d_deltaq_the)
!!      DEALLOCATE( d_s_the, d_dens_the)
      DEALLOCATE(d_deltat_ajs_cv, d_deltaq_ajs_cv)
!
      DEALLOCATE(bils)
      DEALLOCATE(cdragm, cdragh, cldh, cldl)
      DEALLOCATE(cldm, cldq, cldt, qsat2m)
      DEALLOCATE(cldljn, cldmjn, cldhjn, cldtjn, JrNt)
      DEALLOCATE(dthmin, evap, fder, plcl, plfc)
      DEALLOCATE(prw, prlw, prsw, zustar, zu10m, zv10m, rh2m, s_lcl)
      DEALLOCATE(s_pblh, s_pblt, s_therm)
!
!nrlmd+jyg<
      DEALLOCATE(s_pblh_x, s_pblh_w)
      DEALLOCATE(s_lcl_x, s_lcl_w)
!>nrlmd+jyg
!
      DEALLOCATE(slab_wfbils, tpot, tpote, ue)
      DEALLOCATE(uq, ve, vq, zxffonte)
      DEALLOCATE(uwat, vwat)
      DEALLOCATE(zxfqcalving, zxfluxlat)
      DEALLOCATE(zxrunofflic)
      DEALLOCATE(zxtsol, snow_lsc, zxfqfonte, zxqsurf)
      DEALLOCATE(rain_lsc)
      DEALLOCATE(rain_num)
!
      DEALLOCATE(sens_x, sens_w)
      DEALLOCATE(zxfluxlat_x, zxfluxlat_w)
!jyg<
!!      DEALLOCATE(t_x, t_w)
!!      DEALLOCATE(q_x, q_w)
!>jyg
      DEALLOCATE(d_t_vdf_x, d_t_vdf_w)
      DEALLOCATE(d_q_vdf_x, d_q_vdf_w)
      DEALLOCATE(pbl_tke_input)
      DEALLOCATE(t_therm, q_therm, u_therm, v_therm)
      DEALLOCATE(cdragh_x, cdragh_w)
      DEALLOCATE(cdragm_x, cdragm_w)
      DEALLOCATE(kh, kh_x, kh_w)
!
      DEALLOCATE(ptconv)
      DEALLOCATE(ratqsc)
!
      DEALLOCATE(wbeff, convoccur, zmax_th)
      DEALLOCATE(zq2m, zt2m, weak_inversion)
      DEALLOCATE(zt2m_min_mon, zt2m_max_mon)
      DEALLOCATE(t2m_min_mon, t2m_max_mon)
      DEALLOCATE(sens, flwp, fiwp)
      DEALLOCATE(alp_bl_conv,alp_bl_det)
      DEALLOCATE(alp_bl_fluct_m,alp_bl_fluct_tke)
      DEALLOCATE(alp_bl_stat, n2, s2)
      DEALLOCATE(proba_notrig, random_notrig)
      DEALLOCATE(cv_gen)

      DEALLOCATE(dnwd0)
!      DEALLOCATE(upwd, omega, coefh)
      DEALLOCATE(omega)
      DEALLOCATE(epmax_diag)
!      DEALLOCATE(coefm, lambda_th, cldemi)
      DEALLOCATE(lambda_th, cldemi)
      DEALLOCATE(cldfra, cldtau, fiwc)
      DEALLOCATE(fl, re, flwc)
      DEALLOCATE(ref_liq, ref_ice, theta)
      DEALLOCATE(ref_liq_pi, ref_ice_pi)
      DEALLOCATE(zphi, zx_rh)
      DEALLOCATE(pmfd, pmfu)

      DEALLOCATE(t2m, fluxlat)
      DEALLOCATE(fsollw, evap_pot)
      DEALLOCATE(fsolsw, wfbils, wfbilo)
      DEALLOCATE(wfevap,wfrain,wfsnow)

      DEALLOCATE(pmflxr, pmflxs)
      DEALLOCATE(wdtrainA, wdtrainS, wdtrainM)
      DEALLOCATE(upwd, dnwd)
      DEALLOCATE(ep)
      DEALLOCATE(da, mp )
      DEALLOCATE(phi )
      DEALLOCATE(wght_cvfd )
      DEALLOCATE(phi2 )
      DEALLOCATE(d1a, dam)
      DEALLOCATE(ev )
      DEALLOCATE(elij )
      DEALLOCATE(qtaa )
      DEALLOCATE(clw )
      DEALLOCATE(epmlmMm, eplaMm )
      DEALLOCATE(sij )


      DEALLOCATE(prfl, psfl, fraca, Vprecip)
      DEALLOCATE(zw2)

      DEALLOCATE(fluxu, fluxv)
      DEALLOCATE(fluxt)

      DEALLOCATE(uwriteSTD, vwriteSTD)
      DEALLOCATE(wwriteSTD, phiwriteSTD)
      DEALLOCATE(qwriteSTD, twriteSTD, rhwriteSTD)

! ug et d'autres encore:
      DEALLOCATE(beta_prec)
      DEALLOCATE(rneb)
      DEALLOCATE (zxsnow,snowhgt,qsnow,to_ice,sissnow,runoff,albsol3_lic)

      DEALLOCATE (p_tropopause)
      DEALLOCATE (z_tropopause)
      DEALLOCATE (t_tropopause)

      DEALLOCATE(zn2mout)

#ifdef CPP_StratAer
! variables for strat. aerosol CK
      DEALLOCATE (R2SO4)
      DEALLOCATE (DENSO4)
      DEALLOCATE (f_r_wet)
      DEALLOCATE (decfluxaer)
      DEALLOCATE (mdw)
      DEALLOCATE (SO2_lifetime)
      DEALLOCATE (OCS_lifetime)
      DEALLOCATE (alpha_bin)
      DEALLOCATE (piz_bin)
      DEALLOCATE (cg_bin)
      DEALLOCATE (tau_strat_550)
      DEALLOCATE (tau_strat_1020)
      DEALLOCATE (tausum_strat)
      DEALLOCATE (surf_PM25_sulf)
      DEALLOCATE (vsed_aer)
      DEALLOCATE (budg_3D_ocs_to_so2)
      DEALLOCATE (budg_3D_so2_to_h2so4)
      DEALLOCATE (budg_3D_backgr_ocs)
      DEALLOCATE (budg_3D_backgr_so2)
      DEALLOCATE (budg_3D_nucl)
      DEALLOCATE (budg_3D_cond_evap)
      DEALLOCATE (budg_dep_dry_ocs)
      DEALLOCATE (budg_dep_wet_ocs)
      DEALLOCATE (budg_dep_dry_so2)
      DEALLOCATE (budg_dep_wet_so2)
      DEALLOCATE (budg_dep_dry_h2so4)
      DEALLOCATE (budg_dep_wet_h2so4)
      DEALLOCATE (budg_dep_dry_part)
      DEALLOCATE (budg_dep_wet_part)
      DEALLOCATE (budg_emi_ocs)
      DEALLOCATE (budg_emi_so2)
      DEALLOCATE (budg_emi_h2so4)
      DEALLOCATE (budg_emi_part)
      DEALLOCATE (budg_ocs_to_so2)
      DEALLOCATE (budg_so2_to_h2so4)
      DEALLOCATE (budg_h2so4_to_part)
      DEALLOCATE (budg_sed_part)
#endif

END SUBROUTINE phys_local_var_end

END MODULE phys_local_var_mod 
