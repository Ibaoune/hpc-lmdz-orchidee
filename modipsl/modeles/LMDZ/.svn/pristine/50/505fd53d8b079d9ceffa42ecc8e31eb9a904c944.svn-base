! $Id: physiq.F90 2298 2015-06-14 19:13:32Z fairhead $
!#define IO_DEBUG

MODULE phytracr_spl_mod


! Recuperation des morceaux de la physique de Jeronimo specifiques
! du modele d'aerosols d'Olivier n'co.
!
INCLUDE "chem.h"
INCLUDE "chem_spla.h"

  REAL,SAVE  :: scale_param_ssacc  !Scaling parameter for Fine Sea Salt
  REAL,SAVE ::  scale_param_sscoa  !Scaling parameter for Coarse Sea Salt



  REAL, DIMENSION(:),ALLOCATABLE,SAVE :: scale_param_ind !Scaling parameter for industrial emissions of SO2
  REAL, DIMENSION(:),ALLOCATABLE,SAVE :: scale_param_bb  !Scaling parameter for biomas burning (SO2,BC & OM)
  REAL, DIMENSION(:),ALLOCATABLE,SAVE :: scale_param_ff  !Scaling parameter for industrial emissions (fossil fuel)
  REAL, DIMENSION(:),ALLOCATABLE,SAVE :: scale_param_dustacc  !Scaling parameter for Fine Dust
  REAL, DIMENSION(:),ALLOCATABLE,SAVE :: scale_param_dustcoa  !Scaling parameter for Coarse Dust
  REAL, DIMENSION(:),ALLOCATABLE,SAVE :: scale_param_dustsco  !Scaling parameter for SCoarse Dust
  REAL, DIMENSION(:),ALLOCATABLE,SAVE :: param_wstarBLperregion  !parameter for ..
  REAL, DIMENSION(:),ALLOCATABLE,SAVE :: param_wstarWAKEperregion  !parameter for ..
  !$OMP THREADPRIVATE(scale_param_ind,scale_param_bb,scale_param_ff)
  !$OMP THREADPRIVATE(scale_param_dustacc,scale_param_dustcoa,scale_param_dustsco)
  !$OMP THREADPRIVATE(scale_param_ssacc,scale_param_sscoa)
  !$OMP THREADPRIVATE(param_wstarBLperregion,param_wstarWAKEperregion)
  REAL, DIMENSION(:),ALLOCATABLE,SAVE ::dust_ec, u10m_ec, v10m_ec
!$OMP THREADPRIVATE(dust_ec, u10m_ec, v10m_ec)

  CHARACTER*800 fileregionsdimsind
  CHARACTER*800 fileregionsdimsdust
  CHARACTER*800 fileregionsdimsbb
  CHARACTER*800 fileregionsdimswstar
!  CHARACTER*800 filescaleparamsind
!  CHARACTER*800 filescaleparamsdust
!  CHARACTER*800 filescaleparamsbb
  CHARACTER*100 paramname_ind
  CHARACTER*100 paramname_bb
  CHARACTER*100 paramname_ff
  CHARACTER*100 paramname_dustacc
  CHARACTER*100 paramname_dustcoa
  CHARACTER*100 paramname_dustsco
  CHARACTER*100 paramname_ssacc
  CHARACTER*100 paramname_sscoa
  CHARACTER*100 paramname_wstarBL
  CHARACTER*100 paramname_wstarWAKE


  CHARACTER*800 filescaleparams
  CHARACTER*800 paramsname


  !!------------------------ SULFUR emissions ----------------------------
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_so2volc_cont  ! emissions so2 volcan continuous
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_altvolc_cont  ! altitude  so2 volcan continuous
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_so2volc_expl  ! emissions so2 volcan explosive
!$OMP THREADPRIVATE( lmt_so2volc_cont,lmt_altvolc_cont,lmt_so2volc_expl )
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_altvolc_expl  ! altitude  so2 volcan explosive
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_so2ff_l       ! emissions so2 fossil fuel (low)
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_so2ff_h       ! emissions so2 fossil fuel (high)
!$OMP THREADPRIVATE( lmt_altvolc_expl,lmt_so2ff_l,lmt_so2ff_h )
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_so2nff        ! emissions so2 non-fossil fuel
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_so2ba         ! emissions de so2 bateau
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_so2bb_l       ! emissions de so2 biomass burning (low)
!$OMP THREADPRIVATE( lmt_so2nff,lmt_so2ba,lmt_so2bb_l )
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_so2bb_h       ! emissions de so2 biomass burning (high)
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_dmsconc       ! concentration de dms oceanique
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_dmsbio        ! emissions de dms bio
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_h2sbio        ! emissions de h2s bio
!$OMP THREADPRIVATE(lmt_so2bb_h,lmt_dmsconc,lmt_dmsbio,lmt_h2sbio )
  !------------------------- BLACK CARBON emissions ----------------------
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_bcff       ! emissions de BC fossil fuels
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_bcnff      ! emissions de BC non-fossil fuels
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_bcbb_l     ! emissions de BC biomass basses
!$OMP THREADPRIVATE( lmt_bcff,lmt_bcnff,lmt_bcbb_l)
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_bcbb_h     ! emissions de BC biomass hautes
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_bcba       ! emissions de BC bateau
!$OMP THREADPRIVATE(lmt_bcbb_h,lmt_bcba)
  !------------------------ ORGANIC MATTER emissions ---------------------
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_omff     ! emissions de OM fossil fuels
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_omnff    ! emissions de OM non-fossil fuels
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_ombb_l   ! emissions de OM biomass basses
!$OMP THREADPRIVATE( lmt_omff,lmt_omnff,lmt_ombb_l)
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_ombb_h   ! emissions de OM biomass hautes
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_omnat    ! emissions de OM Natural
  REAL , DIMENSION(:),ALLOCATABLE,SAVE :: lmt_omba     ! emissions de OM bateau
  REAL , DIMENSION(:,:),ALLOCATABLE,SAVE :: lmt_sea_salt    ! emissions de OM Natural
!$OMP THREADPRIVATE(lmt_ombb_h,lmt_omnat,lmt_omba,lmt_sea_salt)

!JE20141224 >>
  ! others
  REAL, DIMENSION(:),ALLOCATABLE,SAVE ::  tsol
!$OMP THREADPRIVATE(tsol)
  INTEGER :: ijulday
  LOGICAL , parameter :: edgar = .true.
  INTEGER , parameter :: flag_dms=4
  INTEGER*4  nbjour

      !
! Tracer tendencies, for outputs
!-------------------------------
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_cl  ! Td couche
!. limite/traceur
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_dec
!RomP
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_cv  ! Td
!onvection/traceur
! RomP >>>
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_insc
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_bcscav
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_evapls
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_ls
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_trsp
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_sscav
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_sat
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_uscav
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: qPr,qDi ! concentration tra
!dans pluie,air descente insaturee
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: qPa,qMel
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: qTrdi,dtrcvMA ! conc traceur
!descente air insaturee et td convective MA
!! RomP <<<
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_th  ! Td thermique
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_lessi_impa ! Td du
!lessivage par impaction
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_lessi_nucl ! Td du
!lessivage par nucleation
      REAL,DIMENSION(:,:),ALLOCATABLE,SAVE      :: qPrls      !jyg:
!oncentration tra dans pluie LS a la surf.
      REAL,DIMENSION(:,:),ALLOCATABLE,SAVE      :: d_tr_dry ! Td depot
!sec/traceur (1st layer),ALLOCATABLE,SAVE  jyg
      REAL,DIMENSION(:,:),ALLOCATABLE,SAVE      :: flux_tr_dry ! depot
!sec/traceur (surface),ALLOCATABLE,SAVE    jyg

! Index of each traceur
      INTEGER,SAVE :: id_prec, id_fine, id_coss, id_codu, id_scdu 

!$OMP THREADPRIVATE(d_tr_cl,d_tr_dec,d_tr_cv,d_tr_insc,d_tr_bcscav,d_tr_evapls)
!$OMP THREADPRIVATE(d_tr_ls,d_tr_trsp,d_tr_sscav,d_tr_sat,d_tr_uscav)
!$OMP THREADPRIVATE(qPr,qDi,qPa,qMel,qTrdi,dtrcvMA,d_tr_th,d_tr_lessi_impa)
!$OMP THREADPRIVATE(d_tr_lessi_nucl,qPrls,d_tr_dry,flux_tr_dry)
!$OMP THREADPRIVATE(id_prec,id_fine,id_coss,id_codu,id_scdu)

! JE20141224 <<

      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diff_aod550_tot  ! epaisseur optique total aerosol 550  nm
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod670_tot  ! epaisseur optique total aerosol 670 nm
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod865_tot  ! epaisseur optique total aerosol 865 nm
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diff_aod550_tr2  ! epaisseur optique Traceur 2 aerosol 550 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod670_tr2  ! epaisseur optique Traceur 2 aerosol 670 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod865_tr2  ! epaisseur optique Traceur 2 aerosol 865 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod550_ss  ! epaisseur optique Sels marins aerosol 550 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod670_ss  ! epaisseur optique Sels marins aerosol 670 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod865_ss   ! epaisseur optique Sels marins aerosol 865 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod550_dust ! epaisseur optique Dust aerosol 550 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod670_dust ! epaisseur optique Dust aerosol 670 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod865_dust ! epaisseur optique Dust aerosol 865 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod550_dustsco ! epaisseur optique Dust SCOarse aerosol 550 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod670_dustsco ! epaisseur optique Dust SCOarse aerosol 670 nm, diagnostic
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: diag_aod865_dustsco ! epaisseur optique Dust SCOarse aerosol 865 nm, diagnostic

!$OMP THREADPRIVATE(diff_aod550_tot,diag_aod670_tot,diag_aod865_tot)
!$OMP THREADPRIVATE(diff_aod550_tr2,diag_aod670_tr2,diag_aod865_tr2)
!$OMP THREADPRIVATE(diag_aod550_ss,diag_aod670_ss,diag_aod865_ss,diag_aod550_dust)
!$OMP THREADPRIVATE(diag_aod670_dust,diag_aod865_dust,diag_aod550_dustsco)
!$OMP THREADPRIVATE(diag_aod670_dustsco,diag_aod865_dustsco)


      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod550_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod550_tr2_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod550_ss_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod550_dust_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod550_dustsco_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod670_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod670_tr2_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod670_ss_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod670_dust_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod670_dustsco_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod865_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod865_tr2_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod865_ss_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod865_dust_terra  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod865_dustsco_terra  ! AOD at terra overpass time ( 10.30 local hour)


      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod550_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod550_tr2_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod550_ss_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod550_dust_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod550_dustsco_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod670_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod670_tr2_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod670_ss_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod670_dust_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod670_dustsco_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod865_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod865_tr2_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod865_ss_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod865_dust_aqua  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: aod865_dustsco_aqua  ! AOD at aqua overpass time ( 13.30 local hour)

!$OMP THREADPRIVATE(aod550_aqua,aod550_tr2_aqua,aod550_ss_aqua,aod550_dust_aqua,aod550_dustsco_aqua)
!$OMP THREADPRIVATE(aod670_aqua,aod670_tr2_aqua,aod670_ss_aqua,aod670_dust_aqua,aod670_dustsco_aqua)
!$OMP THREADPRIVATE(aod865_aqua,aod865_tr2_aqua,aod865_ss_aqua,aod865_dust_aqua,aod865_dustsco_aqua)
!$OMP THREADPRIVATE(aod550_terra,aod550_tr2_terra,aod550_ss_terra,aod550_dust_terra,aod550_dustsco_terra)
!$OMP THREADPRIVATE(aod670_terra,aod670_tr2_terra,aod670_ss_terra,aod670_dust_terra,aod670_dustsco_terra)
!$OMP THREADPRIVATE(aod865_terra,aod865_tr2_terra,aod865_ss_terra,aod865_dust_terra,aod865_dustsco_terra)


      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: sconc01 ! surface concentration
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: trm01   ! burden 
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: sconc02 ! surface concentration
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: trm02   ! burden 
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: sconc03 ! surface concentration
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: trm03   ! burden 
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: sconc04 ! surface concentration
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: trm04   ! burden 
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: sconc05 ! surface concentration
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: trm05   ! burden 
!$OMP THREADPRIVATE(sconc01,sconc02,sconc03,sconc04,sconc05)
!$OMP THREADPRIVATE(trm01,trm02,trm03,trm04,trm05)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux01        
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux02       
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux03       
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux04       
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux05       
!$OMP THREADPRIVATE(flux01,flux02,flux03,flux04,flux05)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: ds01         
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: ds02         
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: ds03         
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: ds04         
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: ds05         
!$OMP THREADPRIVATE(ds01,ds02,ds03,ds04,ds05)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dh01         
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dh02         
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dh03         
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dh04         
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dh05         
!$OMP THREADPRIVATE(dh01,dh02,dh03,dh04,dh05)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dtrconv01    
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dtrconv02    
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dtrconv03    
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dtrconv04    
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dtrconv05    
!$OMP THREADPRIVATE(dtrconv01,dtrconv02,dtrconv03,dtrconv04,dtrconv05)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dtherm01     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dtherm02     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dtherm03     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dtherm04     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dtherm05     
!$OMP THREADPRIVATE(dtherm01,dtherm02,dtherm03,dtherm04,dtherm05)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dhkecv01     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dhkecv02     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dhkecv03     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dhkecv04     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dhkecv05     
!$OMP THREADPRIVATE(dhkecv01,dhkecv02,dhkecv03,dhkecv04,dhkecv05)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: d_tr_ds01     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: d_tr_ds02     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: d_tr_ds03     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: d_tr_ds04     
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: d_tr_ds05     
!$OMP THREADPRIVATE(d_tr_ds01,d_tr_ds02,d_tr_ds03,d_tr_ds04,d_tr_ds05)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dhkelsc01    
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dhkelsc02    
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dhkelsc03    
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dhkelsc04    
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: dhkelsc05    
!$OMP THREADPRIVATE(dhkelsc01,dhkelsc02,dhkelsc03,dhkelsc04,dhkelsc05)
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_cv01    
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_cv02    
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_cv03    
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_cv04    
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_cv05    
!$OMP THREADPRIVATE(d_tr_cv01,d_tr_cv02,d_tr_cv03,d_tr_cv04,d_tr_cv05)
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_trsp01  
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_trsp02  
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_trsp03  
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_trsp04  
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_trsp05  
!$OMP THREADPRIVATE(d_tr_trsp01,d_tr_trsp02,d_tr_trsp03,d_tr_trsp04,d_tr_trsp05)
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_sscav01 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_sscav02 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_sscav03 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_sscav04 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_sscav05 
!$OMP THREADPRIVATE(d_tr_sscav01,d_tr_sscav02,d_tr_sscav03,d_tr_sscav04,d_tr_sscav05)
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_sat01   
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_sat02   
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_sat03   
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_sat04   
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_sat05   
!$OMP THREADPRIVATE(d_tr_sat01,d_tr_sat02,d_tr_sat03,d_tr_sat04,d_tr_sat05)
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_uscav01 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_uscav02 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_uscav03 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_uscav04 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_uscav05 
!$OMP THREADPRIVATE(d_tr_uscav01,d_tr_uscav02,d_tr_uscav03,d_tr_uscav04,d_tr_uscav05)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_insc01  
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_insc02  
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_insc03  
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_insc04  
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_insc05  
!$OMP THREADPRIVATE(d_tr_insc01,d_tr_insc02,d_tr_insc03,d_tr_insc04,d_tr_insc05)
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_bcscav01 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_bcscav02 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_bcscav03 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_bcscav04 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_bcscav05 
!$OMP THREADPRIVATE(d_tr_bcscav01,d_tr_bcscav02,d_tr_bcscav03,d_tr_bcscav04,d_tr_bcscav05)
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_evapls01   
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_evapls02   
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_evapls03   
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_evapls04   
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_evapls05   
!$OMP THREADPRIVATE(d_tr_evapls01,d_tr_evapls02,d_tr_evapls03,d_tr_evapls04,d_tr_evapls05)
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_ls01 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_ls02 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_ls03 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_ls04 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_ls05 
!$OMP THREADPRIVATE(d_tr_ls01,d_tr_ls02,d_tr_ls03,d_tr_ls04,d_tr_ls05)

      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_dyn01 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_dyn02 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_dyn03 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_dyn04 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_dyn05 
!$OMP THREADPRIVATE(d_tr_dyn01,d_tr_dyn02,d_tr_dyn03,d_tr_dyn04,d_tr_dyn05)

      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_cl01 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_cl02 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_cl03 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_cl04 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_cl05 
!$OMP THREADPRIVATE(d_tr_cl01,d_tr_cl02,d_tr_cl03,d_tr_cl04,d_tr_cl05)

      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_th01 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_th02 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_th03 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_th04 
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: d_tr_th05 
!$OMP THREADPRIVATE(d_tr_th01,d_tr_th02,d_tr_th03,d_tr_th04,d_tr_th05)

      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: sed_ss3D    ! corresponds to tracer 3
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: sed_dust3D  ! corresponds to tracer 4
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: sed_dustsco3D  ! corresponds to tracer 4
!$OMP THREADPRIVATE(sed_ss3D,sed_dust3D,sed_dustsco3D)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: sed_ss    ! corresponds to tracer 3
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: sed_dust  ! corresponds to tracer 4
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: sed_dustsco  ! corresponds to tracer 4
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: his_g2pgas  ! corresponds to tracer 4
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: his_g2paer  ! corresponds to tracer 4
!$OMP THREADPRIVATE(sed_ss,sed_dust,sed_dustsco,his_g2pgas,his_g2paer)

      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxbb
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxff
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxbcbb
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxbcff
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxbcnff
!$OMP THREADPRIVATE(fluxbb,fluxff,fluxbcbb,fluxbcff,fluxbcnff)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxbcba
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxbc
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxombb
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxomff
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxomnff
!$OMP THREADPRIVATE(fluxbcba,fluxbc,fluxombb,fluxomff,fluxomnff)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxomba
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxomnat
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxom
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxh2sff
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxh2snff
!$OMP THREADPRIVATE(fluxomba,fluxomnat,fluxom,fluxh2sff,fluxh2snff)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso2ff
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso2nff
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso2bb
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso2vol
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso2ba
!$OMP THREADPRIVATE(fluxso2ff,fluxso2nff,fluxso2bb,fluxso2vol,fluxso2ba)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso2
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso4ff
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso4nff
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso4bb
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso4ba
!$OMP THREADPRIVATE(fluxso2,fluxso4ff,fluxso4nff,fluxso4ba,fluxso4bb)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxso4
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxdms
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxh2sbio
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxdustec
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxddfine
!$OMP THREADPRIVATE(fluxso4,fluxdms,fluxh2sbio,fluxdustec,fluxddfine)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxddcoa
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxddsco
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxdd
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxssfine
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxsscoa
!$OMP THREADPRIVATE(fluxddcoa,fluxddsco,fluxdd,fluxssfine,fluxsscoa)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: fluxss
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux_sparam_ind
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux_sparam_bb
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux_sparam_ff
!$OMP THREADPRIVATE(fluxss,flux_sparam_ind,flux_sparam_bb,flux_sparam_ff)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux_sparam_ddfine
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux_sparam_ddcoa
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux_sparam_ddsco
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux_sparam_ssfine
!$OMP THREADPRIVATE(flux_sparam_ddfine,flux_sparam_ddcoa)
!$OMP THREADPRIVATE(flux_sparam_ddsco,flux_sparam_ssfine)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: flux_sparam_sscoa
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: u10m_ss
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: v10m_ss
!$OMP THREADPRIVATE(flux_sparam_sscoa,u10m_ss,v10m_ss)

! Select dust emission scheme ver the Sahara:
!      LOGICAL,PARAMETER,SAVE ::  ok_chimeredust=.FALSE.
      LOGICAL,PARAMETER ::  ok_chimeredust=.TRUE.
!!!!!! !$OMP THREADPRIVATE(ok_chimeredust)

!OH   REAL,SAVE :: scale_param_ssacc  !Scaling parameter for Fine Sea Salt
!OH   REAL,SAVE :: scale_param_sscoa  !Scaling parameter for Coarse Sea Salt
!OH   REAL,ALLOCATABLE,SAVE :: scale_param_ind(nbreg_ind) !Scaling parameter for industrial emissionsi of SO2
!OH   REAL,ALLOCATABLE,SAVE :: scale_param_bb(nbreg_bb)  !Scaling parameter for biomas burning (SO2, BC & OM)
!OH   REAL,ALLOCATABLE,SAVE :: scale_param_ff(nbreg_ind)  !Scaling parameter for industrial emissions (fossil fuel)
!OH   REAL,ALLOCATABLE,SAVE :: scale_param_dustacc(nbreg_dust)  !Scaling parameter for Fine Dust
!OH   REAL,ALLOCATABLE,SAVE :: scale_param_dustcoa(nbreg_dust)  !Scaling parameter for Coarse Dust
!OH   REAL,ALLOCATABLE,SAVE :: scale_param_dustsco(nbreg_dust)  !Scaling parameter for SCoarse Dust
!OH   REAL,ALLOCATABLE,SAVE :: param_wstarBLperregion(nbreg_wstardust)
!OH   REAL,ALLOCATABLE,SAVE :: param_wstarWAKEperregion(nbreg_wstardust)
!!!! !$OMP THREADPRIVATE( scale_param_ssacc, scale_param_sscoa, scale_param_ind, scale_param_bb, scale_param_ff, scale_param_dustacc, scale_param_dustcoa, scale_param_dustsco, param_wstarBLperregion, param_wstarWAKEperregion)


CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE phytracr_spl_ini(klon,nbreg_ind,nbreg_bb,nbreg_dust,nbreg_wstardust)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  IMPLICIT NONE
  INTEGER klon,nbreg_ind,nbreg_bb,nbreg_dust,nbreg_wstardust

  ALLOCATE(  tsol(klon)              )
  fileregionsdimsind='regions_ind_meta'
  fileregionsdimsdust='regions_dustacc_meta'
!  fileregionsdimsdust='regions_dust_meta'
  fileregionsdimsbb='regions_bb_meta'
  fileregionsdimswstar='regions_pwstarwake_meta'
  call  readregionsdims2_spl(nbreg_ind,fileregionsdimsind)
  call  readregionsdims2_spl(nbreg_dust,fileregionsdimsdust)
  call  readregionsdims2_spl(nbreg_bb,fileregionsdimsbb)
  call  readregionsdims2_spl(nbreg_wstardust,fileregionsdimswstar)

!readregions_spl()

  ALLOCATE(scale_param_ind(nbreg_ind))
  ALLOCATE(scale_param_bb(nbreg_bb))
  ALLOCATE(scale_param_ff(nbreg_ind))
  ALLOCATE(scale_param_dustacc(nbreg_dust))
  ALLOCATE(scale_param_dustcoa(nbreg_dust))
  ALLOCATE(scale_param_dustsco(nbreg_dust))
  ALLOCATE(param_wstarBLperregion(nbreg_wstardust))
  ALLOCATE(param_wstarWAKEperregion(nbreg_wstardust))
  ALLOCATE(  dust_ec(klon)           )
  ALLOCATE(  u10m_ec(klon)           )
  ALLOCATE(  v10m_ec(klon)           )
  ALLOCATE(  lmt_so2volc_cont(klon)  )
  ALLOCATE(  lmt_altvolc_cont(klon)  )
  ALLOCATE(  lmt_so2volc_expl(klon)  )
  ALLOCATE(  lmt_altvolc_expl(klon)  )
  ALLOCATE(  lmt_so2ff_l(klon)       )   
  ALLOCATE(  lmt_so2ff_h(klon)       )  
  ALLOCATE(  lmt_so2nff(klon)        )  
  ALLOCATE(  lmt_so2ba(klon)         )  
  ALLOCATE(  lmt_so2bb_l(klon)       )
  ALLOCATE(  lmt_so2bb_h(klon)       )  
  ALLOCATE(  lmt_dmsconc(klon)       )  
  ALLOCATE(  lmt_dmsbio(klon)        )  
  ALLOCATE(  lmt_h2sbio(klon)        )  
  ALLOCATE(  lmt_bcff(klon)          )
  ALLOCATE(  lmt_bcnff(klon)         )
  ALLOCATE(  lmt_bcbb_l(klon)        )
  ALLOCATE(  lmt_bcbb_h(klon)        )
  ALLOCATE(  lmt_bcba(klon)          )
  ALLOCATE(  lmt_omff(klon)          )  
  ALLOCATE(  lmt_omnff(klon)         )  
  ALLOCATE(  lmt_ombb_l(klon)        )  
  ALLOCATE(  lmt_ombb_h(klon)        )  
  ALLOCATE(  lmt_omnat(klon)         )  
  ALLOCATE(  lmt_omba(klon)          )           
  ALLOCATE(lmt_sea_salt(klon,ss_bins)) 




  !temporal hardcoded null inicialization of assimilation emmision factors
  scale_param_ssacc=1.
  scale_param_sscoa=1.
  scale_param_ind(:)=1.
  scale_param_bb(:)=1.
  scale_param_ff(:)=1.
  scale_param_dustacc(:)=1.
  scale_param_dustcoa(:)=1.
  scale_param_dustsco(:)=1.
  param_wstarBLperregion(:)=0.
  param_wstarWAKEperregion(:)=0.



RETURN
END SUBROUTINE phytracr_spl_ini




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE phytracr_spl ( debutphy,lafin,jD_cur,jH_cur,iflag_conv, &  ! I
                      pdtphys,ftsol,                                   &  ! I
                      t_seri,q_seri,paprs,pplay,RHcl,                  &  ! I
                      pmfu, pmfd, pen_u, pde_u, pen_d, pde_d,          &  ! I
                      coefh, cdragh, cdragm, yu1, yv1,                 &  ! I
                      u_seri, v_seri, rlat,rlon,                       &  ! I
                      pphis,pctsrf,pmflxr,pmflxs,prfl,psfl,            &  ! I
                      da,phi,phi2,d1a,dam,mp,ep,sigd,sij,clw,elij,     &  ! I
                      epmlmMm,eplaMm,upwd,dnwd,itop_con,ibas_con,      &  ! I
                      evapls,wdtrainA,  wdtrainM,wght_cvfd,              &  ! I
                      fm_therm, entr_therm, rneb,                      &  ! I
                      beta_fisrt,beta_v1,                              &  ! I
                      zu10m,zv10m,wstar,ale_bl,ale_wake,               &  ! I
                      d_tr_dyn,tr_seri)                                            ! O
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      USE IOIPSL
      USE dimphy
      USE infotrac
      USE indice_sol_mod
      USE write_field_phy
     

      USE mod_phys_lmdz_transfert_para

  USE phys_cal_mod, only: jD_1jan,year_len, mth_len, days_elapsed, jh_1jan, year_cur, &
       mth_cur, phys_cal_update

!
      IMPLICIT none
!

!======================================================================
! Auteur(s) FH
! Objet: Moniteur general des tendances traceurs
!
! Remarques en vrac:
! ------------------
! 1/ le call phytrac se fait avec nqmax-2 donc nous avons bien 
! les vrais traceurs (nbtr) dans phytrac (pas la vapeur ni eau liquide)
!======================================================================
#include "dimensions.h"
#include "chem.h"
#include "chem_spla.h"
#include "YOMCST.h"
#include "YOETHF.h"
#include "paramet.h"
#include "thermcell.h"

!======================================================================

! Arguments:
!
!  EN ENTREE:
!  ==========
!
!  divers:
!  -------
!
      real,intent(in) :: pdtphys  ! pas d'integration pour la physique (seconde)
      REAL, intent(in):: jD_cur, jH_cur
      real, intent(in) ::  ftsol(klon,nbsrf)  ! temperature du sol par type
      real, intent(in) ::  t_seri(klon,klev)  ! temperature
      real, intent(in) ::  u_seri(klon,klev)  ! vent
      real , intent(in) :: v_seri(klon,klev)  ! vent
      real , intent(in) :: q_seri(klon,klev)  ! vapeur d eau kg/kg

LOGICAL,  INTENT(IN)                          :: lafin

      real tr_seri(klon,klev,nbtr) ! traceur  
      real tmp_var(klon,klev) ! auxiliary variable to replace traceur  
      real tmp_var2(klon,nbtr) ! auxiliary variable to replace source
      real tmp_var3(klon,klev,nbtr) ! auxiliary variable 3D  
      real dummy1d ! JE auxiliary variable
      real aux_var2(klon) ! auxiliary variable to replace traceur  
      real aux_var3(klon,klev) ! auxiliary variable to replace traceur  
      real d_tr(klon,klev,nbtr)    ! traceur  tendance
      real sconc_seri(klon,nbtr) ! surface concentration of traceur  
!
      integer nbjour
      save nbjour
!$OMP THREADPRIVATE(nbjour)
!
      INTEGER  masque_aqua_cur(klon)
      INTEGER  masque_terra_cur(klon)
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: masque_aqua  !mask for 1 day
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: masque_terra !
!$OMP THREADPRIVATE(masque_aqua,masque_terra)
!!$OMP THREADPRIVATE(aod550_aqua,aod550_terra,aod670_aqua,aod670_terra)
!!$OMP THREADPRIVATE(aod865_aqua,aod865_terra)

  INTEGER, SAVE :: nbreg_dust, nbreg_ind, nbreg_bb, nbreg_ss,nbreg_wstardust
  !$OMP THREADPRIVATE(nbreg_dust, nbreg_ind, nbreg_bb,nbreg_ss,nbreg_wstardust)



      REAL lmt_dms(klon)           ! emissions de dms

!JE20150518<<
      REAL, DIMENSION(klon_glo)  :: aod550_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod550_tr2_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod550_ss_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod550_dust_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod550_dustsco_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod670_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod670_tr2_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod670_ss_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod670_dust_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod670_dustsco_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod865_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod865_tr2_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod865_ss_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod865_dust_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod865_dustsco_terra_glo  ! AOD at terra overpass time ( 10.30 local hour)

      REAL, DIMENSION(klon_glo)  :: aod550_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod550_tr2_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod550_ss_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod550_dust_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod550_dustsco_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod670_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod670_tr2_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod670_ss_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod670_dust_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod670_dustsco_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod865_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod865_tr2_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod865_ss_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod865_dust_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
      REAL, DIMENSION(klon_glo)  :: aod865_dustsco_aqua_glo  ! AOD at aqua overpass time ( 13.30 local hour)
!!!!!!!!!!!!!
!JE20150518>>




      real , intent(in) :: paprs(klon,klev+1)  ! pression pour chaque inter-couche (en Pa)
      real , intent(in) :: pplay(klon,klev)  ! pression pour le mileu de chaque couche (en Pa)
      real , intent(in) :: RHcl(klon,klev)  ! humidite relativen ciel clair
      real znivsig(klev)  ! indice des couches
      real paire(klon)
      real, intent(in) ::  pphis(klon)
      real, intent(in) ::  pctsrf(klon,nbsrf)
      logical , intent(in) :: debutphy   ! le flag de l'initialisation de la physique
!
!  Scaling Parameters:
!  ----------------------
!
      CHARACTER*50 c_Directory
      CHARACTER*80 c_FileName1 
      CHARACTER*80 c_FileName2
      CHARACTER*130 c_FullName1
      CHARACTER*130 c_FullName2
      INTEGER :: xidx, yidx
      INTEGER,DIMENSION(klon) :: mask_bbreg
      INTEGER,DIMENSION(klon) :: mask_ffso2reg
      INTEGER :: aux_mask1
      INTEGER :: aux_mask2
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: iregion_so4 !Defines regions for SO4
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: iregion_ind  !Defines regions for SO2, BC & OM
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: iregion_bb   !Defines regions for SO2, BC & OM
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: iregion_dust !Defines  dust regions
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: iregion_wstardust !Defines  dust regions
!$OMP THREADPRIVATE(iregion_so4,iregion_ind,iregion_bb,iregion_dust,iregion_wstardust)

!  Emissions:

!
!---------------------------- SEA SALT & DUST emissions ------------------------
      REAL lmt_sea_salt(klon,ss_bins) !Sea salt 0.03-8.0 um
      REAL u10m_ec1(klon),v10m_ec1(klon)
      REAL u10m_ec2(klon),v10m_ec2(klon),dust_ec2(klon)
      REAL dust_ec(klon)
!     new dust emission chimere je20140522
      REAL,DIMENSION(klon),INTENT(IN)                     :: zu10m
      REAL,DIMENSION(klon),INTENT(IN)                     :: zv10m
      REAL,DIMENSION(klon),INTENT(IN)  :: wstar,ale_bl,ale_wake


!
!  Rem : nbtr : nombre de vrais traceurs est defini dans dimphy.h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Dynamique
     !--------
      REAL,DIMENSION(klon,klev,nbtr),INTENT(IN)    :: d_tr_dyn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  convection:
!  -----------
!
      REAL , intent(in) :: pmfu(klon,klev)  ! flux de masse dans le panache montant
      REAL , intent(in) :: pmfd(klon,klev)  ! flux de masse dans le panache descendant
      REAL, intent(in) ::  pen_u(klon,klev) ! flux entraine dans le panache montant
      REAL, intent(in) ::  pde_u(klon,klev) ! flux detraine dans le panache montant
      REAL, intent(in) ::  pen_d(klon,klev) ! flux entraine dans le panache descendant
      REAL, intent(in) ::  pde_d(klon,klev) ! flux detraine dans le panache descendant
!
!  Convection KE scheme:
!  ---------------------
!
!! Variables pour le lessivage convectif
       REAL,DIMENSION(klon,klev),INTENT(IN)     :: da
       REAL,DIMENSION(klon,klev,klev),INTENT(IN):: phi
       REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: phi2
       REAL,DIMENSION(klon,klev),INTENT(IN)      :: d1a,dam
       REAL,DIMENSION(klon,klev),INTENT(IN)     :: mp
       REAL,DIMENSION(klon,klev),INTENT(IN)     :: upwd      ! saturated
!            updraft mass flux
       REAL,DIMENSION(klon,klev),INTENT(IN)     :: dnwd      ! saturated
!            downdraft mass flux
       INTEGER,DIMENSION(klon),INTENT(IN)     :: itop_con
       INTEGER,DIMENSION(klon),INTENT(IN)     :: ibas_con
       REAL,DIMENSION(klon,klev)      :: evapls
       REAL,DIMENSION(klon,klev),INTENT(IN)      :: wdtrainA
       REAL,DIMENSION(klon,klev),INTENT(IN)      :: wdtrainM


       REAL,DIMENSION(klon,klev),INTENT(IN)      :: ep
       REAL,DIMENSION(klon),INTENT(IN)           :: sigd
       REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: sij
       REAL,DIMENSION(klon,klev),INTENT(IN)      :: clw
       REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: elij
       REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: epmlmMm
       REAL,DIMENSION(klon,klev),INTENT(IN)      :: eplaMm
       REAL,DIMENSION(klon,klev),INTENT(IN)      :: wght_cvfd          !RL


!     KE: Tendances de traceurs (Td) et flux de traceurs:
!     ------------------------
       REAL,DIMENSION(klon,klev)      :: Mint
       REAL,DIMENSION(klon,klev,nbtr) :: zmfd1a
       REAL,DIMENSION(klon,klev,nbtr) :: zmfdam
       REAL,DIMENSION(klon,klev,nbtr) :: zmfphi2

!                                                        !tra dans pluie LS a la surf.
!      outputs for cvltr_spl
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_cv_o  
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_trsp_o
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_sscav_o
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_sat_o
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_uscav_o
     !!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_insc_o
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_bcscav_o
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_evapls_o
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_ls_o
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_dyn_o
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_cl_o
       REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_th_o
     !!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!

!$OMP THREADPRIVATE(d_tr_cv_o,d_tr_trsp_o,d_tr_sscav_o,d_tr_sat_o,d_tr_uscav_o)
!$OMP THREADPRIVATE(d_tr_insc_o,d_tr_bcscav_o,d_tr_evapls_o,d_tr_ls_o)
!$OMP THREADPRIVATE(d_tr_dyn_o,d_tr_cl_o,d_tr_th_o)


       INTEGER ::  nsplit
!

     

!
!  Lessivage
!  ---------
!
      REAL, intent(in) ::  pmflxr(klon,klev+1), pmflxs(klon,klev+1)   !--convection
      REAL, intent(in) ::  prfl(klon,klev+1),   psfl(klon,klev+1)     !--large-scale
! JE      REAL pmflxr(klon,klev), pmflxs(klon,klev)   !--convection       ! Titane
! JE      REAL prfl(klon,klev),   psfl(klon,klev)     !--large-scale      ! Titane
      REAL :: ql_incl ! contenu en eau liquide nuageuse dans le nuage ! ql_incl=oliq/rneb
      REAL  :: ql_incloud_ref    ! ref value of in-cloud condensed water content

       REAL,DIMENSION(klon,klev),INTENT(IN)   :: rneb    ! fraction nuageuse (grande echelle)
!

      REAL,DIMENSION(klon,klev) :: beta_fisrt ! taux de conversion
!                                                          ! de l'eau cond (de fisrtilp)
      REAL,DIMENSION(klon,klev) :: beta_v1    ! -- (originale version)
      INTEGER,SAVE  :: iflag_lscav_omp,iflag_lscav
!$OMP THREADPRIVATE(iflag_lscav_omp,iflag_lscav)




!Thermiques:
!----------
      REAL,DIMENSION(klon,klev+1),INTENT(IN)   :: fm_therm
      REAL,DIMENSION(klon,klev),INTENT(IN)     :: entr_therm


!
!  Couche limite:
!  --------------
!
      REAL , intent(in) :: coefh(klon,klev) ! coeff melange CL
      REAL , intent(in) :: cdragh(klon), cdragm(klon)
      REAL, intent(in) ::  yu1(klon)        ! vent dans la 1iere couche
      REAL, intent(in) ::  yv1(klon)        ! vent dans la 1iere couche
!
!
!----------------------------------------------------------------------
      REAL his_ds(klon,nbtr)
      REAL his_dh(klon,nbtr)
      REAL his_dhlsc(klon,nbtr)        ! in-cloud scavenging lsc
      REAL his_dhcon(klon,nbtr)       ! in-cloud scavenging con
      REAL his_dhbclsc(klon,nbtr)      ! below-cloud scavenging lsc
      REAL his_dhbccon(klon,nbtr)      ! below-cloud scavenging con
      REAL trm(klon,nbtr)
!
      REAL u10m_ec(klon), v10m_ec(klon)
!
      REAL his_th(klon,nbtr)
      REAL his_dhkecv(klon,nbtr)
      REAL his_dhkelsc(klon,nbtr)


!
!  Coordonnees
!  -----------
!
      REAL, intent(in) ::  rlat(klon)       ! latitudes pour chaque point 
      REAL, intent(in) ::  rlon(klon)       ! longitudes pour chaque point 
!
      INTEGER i, k, it, j, ig 
!
! DEFINITION OF DIAGNOSTIC VARIABLES
!
      REAL diag_trm(nbtr), diag_drydep(nbtr) 
      REAL diag_wetdep(nbtr), diag_cvtdep(nbtr)
      REAL diag_emissn(nbtr), diag_g2part
      REAL diag_sedimt
      REAL trm_aux(nbtr), src_aux(nbtr)
!
! Variables locales pour effectuer les appels en serie
!----------------------------------------------------
      REAL source_tr(klon,nbtr)
      REAL flux_tr(klon,nbtr)
      REAL m_conc(klon,klev)
!      REAL sed_ss(klon)    ! corresponds to tracer 3
!      REAL sed_dust(klon)  ! corresponds to tracer 4
!      REAL sed_dustsco(klon)  ! corresponds to tracer 4
      REAL henry(nbtr)  !--cste de Henry  mol/l/atm
      REAL kk(nbtr)     !--coefficient de var avec T (K)
      REAL alpha_r(nbtr)!--coefficient d'impaction pour la pluie
      REAL alpha_s(nbtr)!--coefficient d'impaction pour la neige
      REAL vdep_oce(nbtr), vdep_sic(nbtr)
      REAL vdep_ter(nbtr), vdep_lic(nbtr)
      REAL ccntrAA_spla(nbtr)
      REAL ccntrENV_spla(nbtr)
      REAL coefcoli_spla(nbtr)
      REAL dtrconv(klon,nbtr)
      REAL zrho(klon,klev), zdz(klon,klev)
      REAL zalt(klon,klev)
      REAL,DIMENSION(klon,klev)      :: zmasse    ! densité atmosphérique
!     .                                              Kg/m2
      REAL,DIMENSION(klon,klev)      :: ztra_th
      REAL qmin, qmax, aux
!      PARAMETER (qmin=0.0, qmax=1.e33)
      PARAMETER (qmin=1.e33, qmax=-1.e33)

! Variables to save data into file
!----------------------------------
   
      CHARACTER*2 str2
      LOGICAL ok_histrac
!JE2014124      PARAMETER (ok_histrac=.true.)
      PARAMETER (ok_histrac=.false.)
!      PARAMETER (ok_chimeredust=.false.) 
!      PARAMETER (ok_chimeredust=.true.) 
      INTEGER ndex2d(iim*(jjm+1)), ndex3d(iim*(jjm+1)*klev)
      INTEGER nhori1, nhori2, nhori3, nhori4, nhori5, nvert
      INTEGER nid_tra1, nid_tra2, nid_tra3, nid_tra4, nid_tra5
      SAVE nid_tra1, nid_tra2, nid_tra3, nid_tra4, nid_tra5
!$OMP THREADPRIVATE(nid_tra1, nid_tra2, nid_tra3, nid_tra4, nid_tra5)
      INTEGER itra
      SAVE itra                    ! compteur pour la physique
!$OMP THREADPRIVATE(itra)
      INTEGER ecrit_tra, ecrit_tra_h, ecrit_tra_m
      SAVE ecrit_tra, ecrit_tra_h, ecrit_tra_m
!$OMP THREADPRIVATE(ecrit_tra, ecrit_tra_h, ecrit_tra_m)
      REAL presnivs(klev) ! pressions approximat. des milieux couches ( en PA)
      REAL zx_tmp_2d(iim,jjm+1), zx_tmp_3d(iim,jjm+1,klev)
      REAL zx_tmp_fi2d(klon), zx_tmp_fi3d(klon, klev)
!      REAL zx_lon(iim,jjm+1), zx_lat(iim,jjm+1)
      REAL zx_lon_glo(nbp_lon,nbp_lat), zx_lat_glo(nbp_lon,nbp_lat)
      REAL zsto, zout, zout_h, zout_m, zjulian

!------Molar Masses
      REAL masse(nbtr)
!
      REAL fracso2emis                              !--fraction so2 emis en so2
      PARAMETER (fracso2emis=0.95) 
      REAL frach2sofso2                             !--fraction h2s from so2
      PARAMETER (frach2sofso2=0.0426)
!
!  Controles
!-------------
      LOGICAL convection,lessivage,lminmax,lcheckmass
      DATA convection,lessivage,lminmax,lcheckmass &
          /.true.,.true.,.true.,.false./
!
      REAL xconv(nbtr)
!
      LOGICAL anthropo, bateau, edgar
      DATA anthropo,bateau,edgar/.true.,.true.,.true./
!
!c bc_source
      INTEGER kminbc, kmaxbc
!JE20150715      PARAMETER (kminbc=3, kmaxbc=5)
      PARAMETER (kminbc=4, kmaxbc=7)
!
      REAL tr1_cont, tr2_cont, tr3_cont, tr4_cont
!
! JE for updating in  cltrac
      REAL,DIMENSION(klon,klev)             :: delp     ! epaisseur de couche (Pa)
!JE20140507      REAL,DIMENSION(klon,nbtr)       :: d_tr_dry ! Td depot sec/traceur (1st layer),ALLOCATABLE,SAVE  jyg
!JE20140507      REAL,DIMENSION(klon,nbtr)        ::  flux_tr_dry
!      SAVE  d_tr_dry
!! JE for include gas to particle conversion in output
!      REAL his_g2pgas(klon)      ! gastoparticle in gas units (check!)
!      REAL his_g2paer(klon)      ! gastoparticle in aerosol units (check!)
!
      INTEGER ,intent(in) :: iflag_conv
      LOGICAL iscm3  ! debug variable. for checkmass ! JE

!------------------------------------------------------------------------
!  only to compute time consumption of each process
!----
      INTEGER clock_start,clock_end,clock_rate,clock_start_spla
      INTEGER clock_end_outphytracr,clock_start_outphytracr
      INTEGER ti_init,dife,ti_inittype,ti_inittwrite
      INTEGER ti_spla,ti_emis,ti_depo,ti_cltr,ti_ther
      INTEGER ti_sedi,ti_gasp,ti_wetap,ti_cvltr,ti_lscs,ti_brop,ti_outs
      INTEGER ti_nophytracr,clock_per_max
      REAL tia_init,tia_inittype,tia_inittwrite
      REAL tia_spla,tia_emis,tia_depo,tia_cltr,tia_ther
      REAL tia_sedi,tia_gasp,tia_wetap,tia_cvltr,tia_lscs
      REAL tia_brop,tia_outs
      REAL tia_nophytracr
 
      SAVE tia_init,tia_inittype,tia_inittwrite
      SAVE tia_spla,tia_emis,tia_depo,tia_cltr,tia_ther
      SAVE tia_sedi,tia_gasp,tia_wetap,tia_cvltr,tia_lscs
      SAVE tia_brop,tia_outs
      SAVE ti_nophytracr
      SAVE tia_nophytracr
      SAVE clock_end_outphytracr,clock_start_outphytracr
      SAVE clock_per_max
      LOGICAL logitime
!$OMP THREADPRIVATE(tia_init,tia_inittype,tia_inittwrite)
!$OMP THREADPRIVATE(tia_spla,tia_emis,tia_depo,tia_cltr,tia_ther)
!$OMP THREADPRIVATE(tia_sedi,tia_gasp,tia_wetap,tia_cvltr,tia_lscs)
!$OMP THREADPRIVATE(tia_brop,tia_outs)
!$OMP THREADPRIVATE(ti_nophytracr)
!$OMP THREADPRIVATE(tia_nophytracr)
!$OMP THREADPRIVATE(clock_end_outphytracr,clock_start_outphytracr)
!$OMP THREADPRIVATE(clock_per_max)

!     utils parallelization
      REAL :: auxklon_glo(klon_glo)
      INTEGER :: iauxklon_glo(klon_glo)
      REAL, DIMENSION(klon_glo,nbp_lev) :: auxklonnbp_lev
      REAL, DIMENSION(klon_glo,nbp_lev,nbtr)  :: auxklonklevnbtr_glo
      REAL,DIMENSION(nbp_lon,nbp_lat) ::  zx_tmp_2d_glo
      REAL,DIMENSION(nbp_lon,nbp_lat,nbp_lev) :: zx_tmp_3d_glo
      REAL,DIMENSION(klon_glo) :: zx_tmp_fi2d_glo
      REAL,DIMENSION(klon_glo , nbp_lev) :: zx_tmp_fi3d_glo
      REAL,DIMENSION(klon_glo,nbtr) :: auxklonnbtr_glo



      source_tr=0.



      if (debutphy) then
#ifdef IOPHYS_DUST
         CALL iophys_ini
#endif
         nbreg_ind=1
         nbreg_bb=1
         nbreg_dust=1
         nbreg_wstardust=1
         CALL phytracr_spl_ini(klon,nbreg_ind,nbreg_bb,nbreg_dust,nbreg_wstardust)
      endif


#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRA'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif

  


  ijulday=jD_cur-jD_1jan+1
  nbjour = 1

  paramname_ind='ind'
  paramname_bb='bb'
  paramname_ff='ind'
  paramname_dustacc='dustacc'
  paramname_dustcoa='dustcoasco'
  paramname_dustsco='dustcoasco'
!  paramname_dustacc='dust'
!  paramname_dustcoa='dust'
!  paramname_dustsco='dust'
  paramname_wstarBL='pwstarbl'
  paramname_wstarWAKE='pwstarwake'
  paramname_ssacc='ssacc'
  paramname_sscoa='sscoa'

  filescaleparams='modvalues.nc'
  CALL readscaleparamsnc_spl(scale_param_ind,                        &
        nbreg_ind, paramname_ind,                                    &
        scale_param_ff, nbreg_ind,paramname_ff,                      &
        scale_param_bb, nbreg_bb,paramname_bb,                       &
        scale_param_dustacc, nbreg_dust,paramname_dustacc,           &
        scale_param_dustcoa, nbreg_dust,paramname_dustcoa,           &
        scale_param_dustsco, nbreg_dust,paramname_dustsco,           &
        param_wstarBLperregion, nbreg_wstardust, paramname_wstarBL, &
        param_wstarWAKEperregion, nbreg_wstardust, paramname_wstarWAKE, &
        scale_param_ssacc  ,  paramname_ssacc,                    &
        scale_param_sscoa  ,  paramname_sscoa,                    &
           filescaleparams,ijulday,jH_cur, pdtphys,debutphy)
! add seasalt

  print *,'JE : check scale_params'

  print *, 'nbreg_ind', nbreg_ind   
  print *, 'nbreg_dust', nbreg_dust  
  print *, 'nbreg_bb', nbreg_bb   
  print *, 'ind', scale_param_ind   
  print *, 'dustacc', scale_param_dustacc  
  print *, 'dustcoa', scale_param_dustcoa  
  print *, 'dustsco', scale_param_dustsco
  print *, 'wstardustBL', param_wstarBLperregion
  print *, 'wstardustWAKE', param_wstarWAKEperregion
  print *, 'ff', scale_param_ff  
  print *, 'bb', scale_param_bb  
  print *, 'ssacc', scale_param_ssacc
  print *, 'sscoa', scale_param_sscoa

  print *,'JE: before read_newemissions '
  print *,'JE: jD_cur:',jD_cur,' ijulday:',ijulday,' jH_cur:',jH_cur,' pdtphys:',pdtphys
  print *,'JE: now read_newemissions:'
  print *,'lmt_so2ff_l AVANT' , MINVAL(lmt_so2ff_l), MAXVAL(lmt_so2ff_l)
  call read_newemissions(ijulday,jH_cur ,edgar, flag_dms,debutphy, & !I
                         pdtphys, lafin, nbjour, pctsrf,  &       !I
                         t_seri, rlat, rlon, &                         !I
                         pmflxr, pmflxs, prfl, psfl, &            !I
                                 u10m_ec, v10m_ec, dust_ec, &     !O
                                 lmt_sea_salt, lmt_so2ff_l, &     !O
                                 lmt_so2ff_h, lmt_so2nff, &       !O
                                 lmt_so2ba, lmt_so2bb_l, lmt_so2bb_h, &  !O
                                 lmt_so2volc_cont, lmt_altvolc_cont, &   !O
                                 lmt_so2volc_expl, lmt_altvolc_expl, &   !O
                                 lmt_dmsbio, lmt_h2sbio, lmt_dmsconc, &  !O
                                 lmt_bcff, lmt_bcnff, lmt_bcbb_l, &      !O
                                 lmt_bcbb_h, lmt_bcba, lmt_omff, &       !O
                                 lmt_omnff, lmt_ombb_l, lmt_ombb_h, &    !O
                                 lmt_omnat, lmt_omba)                    !O


  print *,'Check emissions'
  print *,'lmt_so2ff_l' , MINVAL(lmt_so2ff_l), MAXVAL(lmt_so2ff_l)
  print *,'lmt_so2ff_h' , MINVAL(lmt_so2ff_h), MAXVAL(lmt_so2ff_h)
  print *,'lmt_so2nff' , MINVAL(lmt_so2nff), MAXVAL(lmt_so2nff)
  print *,'lmt_so2ba' , MINVAL(lmt_so2ba), MAXVAL(lmt_so2ba)
  print *,'lmt_so2bb_l' , MINVAL(lmt_so2bb_l), MAXVAL(lmt_so2bb_l)
  print *,'lmt_so2bb_h' , MINVAL(lmt_so2bb_h), MAXVAL(lmt_so2bb_h)
  print *,'lmt_so2volc_cont' , MINVAL(lmt_so2volc_cont), MAXVAL(lmt_so2volc_cont)
  print *,'lmt_altvolc_cont' , MINVAL(lmt_altvolc_cont), MAXVAL(lmt_altvolc_cont)
  print *,'lmt_so2volc_expl' , MINVAL(lmt_so2volc_expl), MAXVAL(lmt_so2volc_expl)
  print *,'lmt_altvolc_expl' , MINVAL(lmt_altvolc_expl), MAXVAL(lmt_altvolc_expl)
  print *,'lmt_dmsbio' , MINVAL(lmt_dmsbio), MAXVAL(lmt_dmsbio)
  print *,'lmt_h2sbio' , MINVAL(lmt_h2sbio), MAXVAL(lmt_h2sbio)
  print *,'lmt_dmsconc' , MINVAL(lmt_dmsconc), MAXVAL(lmt_dmsconc)
  print *,'lmt_bcff' , MINVAL(lmt_bcff), MAXVAL(lmt_bcff)
  print *,'lmt_bcnff' , MINVAL(lmt_bcnff), MAXVAL(lmt_bcnff)
  print *,'lmt_bcbb_l' , MINVAL(lmt_bcbb_l), MAXVAL(lmt_bcbb_l)
  print *,'lmt_bcbb_h' , MINVAL(lmt_bcbb_h), MAXVAL(lmt_bcbb_h)
  print *,'lmt_bcba' , MINVAL(lmt_bcba), MAXVAL(lmt_bcba)
  print *,'lmt_omff' , MINVAL(lmt_omff), MAXVAL(lmt_omff)
  print *,'lmt_omnff' , MINVAL(lmt_omnff), MAXVAL(lmt_omnff)
  print *,'lmt_ombb_l' , MINVAL(lmt_ombb_l), MAXVAL(lmt_ombb_l)
  print *,'lmt_ombb_h' , MINVAL(lmt_ombb_h), MAXVAL(lmt_ombb_h)
  print *,'lmt_omnat' , MINVAL(lmt_omnat), MAXVAL(lmt_omnat)
  print *,'lmt_omba' , MINVAL(lmt_omba), MAXVAL(lmt_omba)
  print *,'JE iflag_con',iflag_conv


!JE_dbg
   do i=1,klon
      tsol(i)=0.0
      do j=1,nbsrf
          tsol(i)=tsol(i)+ftsol(i,j)*pctsrf(i,j)
      enddo
   enddo


!======================================================================
!  INITIALISATIONS
!======================================================================
!             CALL checknanqfi(da(:,:),1.,-1.,' da_ before
!     . phytracr_inphytracr')

!
! computing time
!        logitime=.true.
        logitime=.false.
        IF (logitime) THEN
        clock_start=0
        clock_end=0
        clock_rate=0
       CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate,COUNT_MAX=clock_per_max)
        CALL SYSTEM_CLOCK(COUNT=clock_start_spla)
        clock_start=clock_start_spla
        clock_end_outphytracr=clock_start_spla
        ENDIF


! Definition of tracers index.
      print*,'OK ON PASSSE BIEN LA'
      CALL minmaxsource(source_tr,qmin,qmax,'A1 maxsource init phytracr')


      IF (debutphy) THEN
        id_prec=-1
        id_fine=-1
        id_coss=-1
        id_codu=-1
        id_scdu=-1
       !print *,nbtr
       do it=1,nbtr
        print *, it, tname(it+nqo)
        if (tname(it+nqo) == 'PREC' ) then
            id_prec=it
        endif
        if (tname(it+nqo) == 'FINE' ) then
            id_fine=it
        endif
        if (tname(it+nqo) == 'COSS' ) then
            id_coss=it
        endif
        if (tname(it+nqo) == 'CODU' ) then
            id_codu=it
        endif
        if (tname(it+nqo) == 'SCDU' ) then
            id_scdu=it
        endif
       enddo
       ! check consistency with dust emission scheme:
       if (ok_chimeredust) then
          if (.not.( id_scdu>0 .and. id_codu>0 .and. id_fine>0)) then
             call abort_gcm('phytracr_mod', 'pb in ok_chimdust 0',1)
          endif
       else 
          if (id_scdu>0) then 
       call abort_gcm('phytracr_mod', 'pb in ok_chimdust 1 SCDU',1)
          endif
          if ( (id_codu .le. 0) .or. ( id_fine.le.0)  ) then  
          call abort_gcm('phytracr_mod', 'pb in ok_chimdust 1',1) 
          endif
       endif


       !print *,id_prec,id_fine,id_coss,id_codu,id_scdu
       ENDIF






!---fraction of tracer that is convected (Tiedke)
      xconv(:)=0.
      if(id_prec>0)  xconv(id_prec)=0.8
      if(id_fine>0)  xconv(id_fine)=0.5
      if(id_coss>0)  xconv(id_coss)=0.5
      if(id_codu>0)  xconv(id_codu)=0.6
      if(id_scdu>0)  xconv(id_scdu)=0.6  !!JE fix

      masse(:)=1.
      if(id_prec>0)  masse(id_prec)=32.
      if(id_fine>0)  masse(id_fine)=6.02e23
      if(id_coss>0)  masse(id_coss)=6.02e23
      if(id_codu>0)  masse(id_codu)=6.02e23 
      if(id_scdu>0)  masse(id_scdu)=6.02e23 

      henry(:)=0.
      if(id_prec>0)  henry(id_prec)=1.4
      if(id_fine>0)  henry(id_fine)=0.0
      if(id_coss>0)  henry(id_coss)=0.0
      if(id_codu>0)  henry(id_codu)=0.0
      if(id_scdu>0)  henry(id_scdu)=0.0
      !henry= (/1.4, 0.0, 0.0, 0.0/)
      kk(:)=0.
      if(id_prec>0)  kk(id_prec)=2900.
      if(id_fine>0)  kk(id_fine)=0.0
      if(id_coss>0)  kk(id_coss)=0.0
      if(id_codu>0)  kk(id_codu)=0.0
      if(id_scdu>0)  kk(id_scdu)=0.0
      !kk = (/2900., 0., 0., 0./)
      alpha_r(:)=0.
      if(id_prec>0)  alpha_r(id_prec)=0.0
      if(id_fine>0)  alpha_r(id_fine)=0.001
      if(id_coss>0)  alpha_r(id_coss)=0.001
      if(id_codu>0)  alpha_r(id_codu)=0.001
      if(id_scdu>0)  alpha_r(id_scdu)=0.001  !JE fix
      alpha_s(:)=0.
      if(id_prec>0)  alpha_s(id_prec)=0.0
      if(id_fine>0)  alpha_s(id_fine)=0.01
      if(id_coss>0)  alpha_s(id_coss)=0.01
      if(id_codu>0)  alpha_s(id_codu)=0.01
      if(id_scdu>0)  alpha_s(id_scdu)=0.01  !JE fix

!      alpha_r =  (/0., 0.001, 0.001, 0.001/)
!      alpha_s = (/0., 0.01, 0.01, 0.01/)

! nhl      DATA vdep_oce /0.7, 0.05, 1.2, 1.2/
! nhl vdep_oce for tr1 is a weighted average of dms and so2 dep velocities
      !vdep_oce = (/0.28, 0.28, 1.2, 1.2/)
      vdep_oce(:)=0.
      if(id_prec>0)  vdep_oce(id_prec) = 0.28
      if(id_fine>0)  vdep_oce(id_fine) = 0.28
      if(id_coss>0)  vdep_oce(id_coss) = 1.2
      if(id_codu>0)  vdep_oce(id_codu) = 1.2
      if(id_scdu>0)  vdep_oce(id_scdu) = 1.2
      vdep_sic(:)=0.
      if(id_prec>0)  vdep_sic(id_prec) = 0.2
      if(id_fine>0)  vdep_sic(id_fine) = 0.17
      if(id_coss>0)  vdep_sic(id_coss) = 1.2
      if(id_codu>0)  vdep_sic(id_codu) = 1.2
      if(id_scdu>0)  vdep_sic(id_scdu) = 1.2

      !vdep_sic = (/0.2, 0.17, 1.2, 1.2/)     
      !vdep_ter = (/0.3, 0.14, 1.2, 1.2/)
      vdep_ter(:)=0.
      if(id_prec>0)  vdep_ter(id_prec) = 0.3
      if(id_fine>0)  vdep_ter(id_fine) = 0.14
      if(id_coss>0)  vdep_ter(id_coss) = 1.2
      if(id_codu>0)  vdep_ter(id_codu) = 1.2
      if(id_scdu>0)  vdep_ter(id_scdu) = 1.2

      vdep_lic(:)=0.
      if(id_prec>0)  vdep_lic(id_prec) = 0.2
      if(id_fine>0)  vdep_lic(id_fine) = 0.17
      if(id_coss>0)  vdep_lic(id_coss) = 1.2
      if(id_codu>0)  vdep_lic(id_codu) = 1.2
      if(id_scdu>0)  vdep_lic(id_scdu) = 1.2


      ! convective KE lessivage aer params:
      ccntrAA_spla(:)=0.
      if(id_prec>0)  ccntrAA_spla(id_prec)=-9999.
      if(id_fine>0)  ccntrAA_spla(id_fine)=0.7
      if(id_coss>0)  ccntrAA_spla(id_coss)=1.0
      if(id_codu>0)  ccntrAA_spla(id_codu)=0.7
      if(id_scdu>0)  ccntrAA_spla(id_scdu)=0.7

      ccntrENV_spla(:)=0.
      if(id_prec>0)  ccntrENV_spla(id_prec)=-9999.
      if(id_fine>0)  ccntrENV_spla(id_fine)=0.7
      if(id_coss>0)  ccntrENV_spla(id_coss)=1.0
      if(id_codu>0)  ccntrENV_spla(id_codu)=0.7
      if(id_scdu>0)  ccntrENV_spla(id_scdu)=0.7

      coefcoli_spla(:)=0.
      if(id_prec>0)  coefcoli_spla(id_prec)=-9999.
      if(id_fine>0)  coefcoli_spla(id_fine)=0.001
      if(id_coss>0)  coefcoli_spla(id_coss)=0.001
      if(id_codu>0)  coefcoli_spla(id_codu)=0.001
      if(id_scdu>0)  coefcoli_spla(id_scdu)=0.001

      !vdep_lic = (/0.2, 0.17, 1.2, 1.2/)      
!

      iscm3=.false.
      if (debutphy) then
!$OMP MASTER
         CALL suphel
         print *, 'let s check nbtr=', nbtr
! JE before put in zero
      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan init phytracr')
        ENDDO        
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'minmax init phytracr')
        ENDDO
        CALL minmaxsource(source_tr,qmin,qmax,'maxsource init phytracr')
      ENDIF
! JE   initializon to cero the tracers     
!         DO it=1, nbtr
!            tr_seri(:,:,it)=0.0
!         ENDDO
! JE end     
! Initializing to zero tr_seri for comparison purposes
!        tr_seri(:,:,:)=0.0
!
!        DO it=1,nbtr
!           trm_aux(it)=0.0
!           src_aux(it)=0.0
!           diag_trm(it)=0.0
!           diag_drydep(it)=0.0
!           diag_wetdep(it)=0.0
!           diag_cvtdep(it)=0.0
!           diag_emissn(it)=0.0
!        ENDDO
!        diag_g2part=0.0
         print *,'PREPARE FILES TO SAVE VARIABLES'
!
         nbjour=30
         ecrit_tra =   NINT(86400./pdtphys)                    !--1-day  average
         ecrit_tra_h = NINT(86400./pdtphys*0.25)               !--6-hour average
         ecrit_tra_m = NINT(86400./pdtphys*FLOAT(nbjour))      !--1-mth  average
         print *,'ecrit_tra=', pdtphys, ecrit_tra

         IF (ok_histrac) THEN
           IF (is_mpi_root .AND. is_omp_root) THEN
  
           itra=0
!
           CALL ymds2ju(1900, 1, 1, 0.0, zjulian)
!
!           print *, 'klon,iim,jjm+1 = ',klon,iim,jjm+1
           print *, 'glo klon,iim,jjm+1 = ',klon_glo,nbp_lon,nbp_lat
           CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,rlon,zx_lon_glo)
!
!           DO i = 1, iim
           DO i = 1, nbp_lon
             zx_lon_glo(i,1) = rlon(i+1)
             zx_lon_glo(i,nbp_lat) = rlon(i+1)
           ENDDO
!
      CALL histbeg("histrac_spl", nbp_lon,zx_lon_glo,            &
                       nbp_lat,zx_lat_glo,                       &
                       1,nbp_lon,1,nbp_lat, 0, zjulian, pdtphys, &
                       nhori1, nid_tra1)
!
      CALL histbeg("lessivage_spl", nbp_lon,zx_lon_glo,            &
                       nbp_lat,zx_lat_glo,                         &
                       1,nbp_lon,1,nbp_lat, 0, zjulian, pdtphys,   &
                       nhori2, nid_tra2)
! 
      CALL histbeg("traceur_spl", nbp_lon,zx_lon_glo,               &
                       nbp_lat,zx_lat_glo,                         &
                      1,nbp_lon,1,nbp_lat, 0, zjulian, pdtphys,    &
                       nhori3, nid_tra3)
!
      CALL histvert(nid_tra1, "presnivs", "Vertical levels", "mb",  &
                      nbp_lev, presnivs, nvert)
!
      CALL histvert(nid_tra2, "presnivs", "Vertical levels", "mb",  &
                      nbp_lev, presnivs, nvert)
!
      CALL histvert(nid_tra3, "presnivs", "Vertical levels", "mb",  &
                      nbp_lev, presnivs, nvert)
!
           zsto = pdtphys
           zout = pdtphys * FLOAT(ecrit_tra)
           zout_h = pdtphys * FLOAT(ecrit_tra_h)
           zout_m = pdtphys * FLOAT(ecrit_tra_m)
           print *,'zsto zout=', zsto, zout

!
!----------------- HISTORY FILES OF TRACER EMISSIONS -------------------
!
! HISTRAC
!
       CALL histdef(nid_tra1, "fluxbb", "Flux BB", "mg/m2/s",       &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,     &
                       "ave(X)", zsto,zout)                         
!                                                                   
      CALL histdef(nid_tra1, "fluxff", "Flux FF", "mg/m2/s",        & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,      &
                       "ave(X)", zsto,zout)                          
!                                                                    
      CALL histdef(nid_tra1, "fluxbcbb", "Flux BC-BB", "mg/m2/s",    &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,       &
                       "ave(X)", zsto,zout)                           
!                                                                     
      CALL histdef(nid_tra1, "fluxbcff", "Flux BC-FF", "mg/m2/s",     &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,        &
                       "ave(X)", zsto,zout)                            
!                                                                      
      CALL histdef(nid_tra1, "fluxbcnff", "Flux BC-NFF", "mg/m2/s",    &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                             
!                                                                       
      CALL histdef(nid_tra1, "fluxbcba", "Flux BC-BA", "mg/m2/s",       &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,          &
                       "ave(X)", zsto,zout)                              
!                                                                        
      CALL histdef(nid_tra1, "fluxbc", "Flux BC", "mg/m2/s",    &         
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,  &
                       "ave(X)", zsto,zout)                      
!                                                                
      CALL histdef(nid_tra1, "fluxombb", "Flux OM-BB", "mg/m2/s" ,  &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,      &
                       "ave(X)", zsto,zout)                          
!                                                                    
      CALL histdef(nid_tra1, "fluxomff", "Flux OM-FF", "mg/m2/s",    &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,       &
                       "ave(X)", zsto,zout)                           
!                                                                     
      CALL histdef(nid_tra1, "fluxomnff", "Flux OM-NFF", "mg/m2/s",  & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,       &
                       "ave(X)", zsto,zout)                           
!                                                                     
      CALL histdef(nid_tra1, "fluxomba", "Flux OM-BA", "mg/m2/s",    & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,       &
                       "ave(X)", zsto,zout)                           
!                                                                     
      CALL histdef(nid_tra1, "fluxomnat", "Flux OM-NT", "mg/m2/s",   & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,       &
                       "ave(X)", zsto,zout)                           
!                                                                     
      CALL histdef(nid_tra1, "fluxom", "Flux OM", "mg/m2/s",         & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,       &
                       "ave(X)", zsto,zout)                           
!                                                                     
      CALL histdef(nid_tra1,"fluxh2sff","Flux H2S FF","mgS/m2/s",    & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,       &
                       "ave(X)", zsto,zout)                           
!                                                                     
      CALL histdef(nid_tra1,"fluxh2snff","Flux H2S non-FF",          & 
                       "mgS/m2/s",nbp_lon,nbp_lat,nhori1, 1,1,1,     &
                        -99, 32,                                     & 
                       "ave(X)", zsto,zout)                           
!                                                                     
      CALL histdef(nid_tra1,"fluxso2ff","Flux SO2 FF","mgS/m2/s",    &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,       &
                       "ave(X)", zsto,zout)                           
!                                                                     
      CALL histdef(nid_tra1,"fluxso2nff","Flux SO2 non-FF",          & 
                       "mgS/m2/s",nbp_lon,nbp_lat,nhori1, 1,1,1,     &
                        -99, 32,                                     &
                       "ave(X)", zsto,zout)                           
!                                                                      
      CALL histdef(nid_tra1, "fluxso2bb", "Flux SO2 BB","mgS/m2/s",   & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,        &  
                       "ave(X)", zsto,zout)                           
!                                                                      
      CALL histdef(nid_tra1,"fluxso2vol","Flux SO2 Vol","mgS/m2/s",    &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1, "fluxso2ba", "Flux SO2 Ba","mgS/m2/s",    & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1, "fluxso2", "Flux SO2","mgS/m2/s",         & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1,"fluxso4ff","Flux SO4 FF","mgS/m2/s",      & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1,"fluxso4nff","Flux SO4 non-FF",            & 
                   "mgS/m2/s", nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32, &
                   "ave(X)", zsto,zout)                                
!                                                                       
      CALL histdef(nid_tra1, "fluxso4bb", "Flux SO4 BB","mgS/m2/s",    & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1, "fluxso4ba", "Flux SO4 Ba","mgS/m2/s",    & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1, "fluxso4", "Flux SO4","mgS/m2/s",         & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1, "fluxdms", "Flux DMS", "mgS/m2/s",        & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1,"fluxh2sbio","Flux H2S Bio","mgS/m2/s",    & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1, "fluxdustec",                             & 
                                      "Flux Dust EC", "mg/m2/s",       &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                             
!                                                                       
      CALL histdef(nid_tra1,"fluxddfine","DD Fine Mode","mg/m2/s",     &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1,"fluxddcoa","DD Coarse Mode","mg/m2/s",    & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1,"fluxddsco","DD SCoarse Mode","mg/m2/s",   & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1,"fluxdd","Flux DD","mg/m2/s",              & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1,"fluxssfine","SS Fine Mode","mg/m2/s",     & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1,"fluxsscoa","SS Coarse Mode","mg/m2/s",    & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
      CALL histdef(nid_tra1,"fluxss","Flux SS","mg/m2/s",              & 
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                            
!                                                                       
!nhl          CALL histdef(nid_tra1,"fluxso4chem","SO4 chem prod",      
!nhl    .                  "gAer/kgAir",
!nhl    .                  nbp_lon,nbp_lat,nhori1, nbp_lev,1,nbp_lev,nvert, 32,
!nhl    .                  "ave(X)", zsto,zout)
!
          CALL histdef(nid_tra1,"flux_sparam_ind","Ind emiss",      &
                       "mg/m2/s",                                   &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,      & 
                       "ave(X)", zsto,zout)                           
!                                                                    
          CALL histdef(nid_tra1,"flux_sparam_bb","BB emiss",        &
                       "mg/m2/s",                                   &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,      &
                       "ave(X)", zsto,zout)                          
!                                                                    
          CALL histdef(nid_tra1,"flux_sparam_ff","FF emiss",        &
                       "mg/m2/s",                                   &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,      &
                       "ave(X)", zsto,zout)                          
!                                                                    
          CALL histdef(nid_tra1,"flux_sparam_ddfine","DD fine emiss",  &
                       "mg/m2/s",                                      &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                             
!                                                                       
          CALL histdef(nid_tra1,"flux_sparam_ddcoa","DD coarse emiss",  & 
                       "mg/m2/s",                                       &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,          &
                       "ave(X)", zsto,zout)                              
!                                                                        
          CALL histdef(nid_tra1,"flux_sparam_ddsco","DD Scoarse emiss", &
                       "mg/m2/s",                                       &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,          &
                       "ave(X)", zsto,zout)                              
!                                                                        
          CALL histdef(nid_tra1,"flux_sparam_ssfine","SS fine emiss",   &
                       "mg/m2/s",                                       &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,          &
                       "ave(X)", zsto,zout)                              
!                                                                        
          CALL histdef(nid_tra1,"flux_sparam_sscoa","SS coarse emiss",  &
                       "mg/m2/s",                                       &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,          &
                       "ave(X)", zsto,zout)                              
!                                                                        
          CALL histdef(nid_tra1,"u10m","Zonal wind at 10 m",            &
                       "m/s",                                           &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,          &
                       "ave(X)", zsto,zout)                              
!                                                                        
          CALL histdef(nid_tra1,"v10m","Meridional wind at 10 m",       &
                       "m/s",                                           &
                       nbp_lon,nbp_lat,nhori1, 1,1,1, -99, 32,          &
                       "ave(X)", zsto,zout)                              
!                                                                        
!nhl          CALL histdef(nid_tra1,"flux_sparam_sulf","SO4 chem prod", 
!nhl    .                  "gAer/kgAir",
!nhl    .                  nbp_lon,nbp_lat,nhori1, nbp_lev,1,nbp_lev,nvert, 32,
!nhl    .                  "ave(X)", zsto,zout)
!
! TRACEUR
!
          CALL histdef(nid_tra3, "taue550", "Tau ext 550", " ",           &
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,            &
                       "ave(X)", zsto,zout)                                  
!                                                                           
          CALL histdef(nid_tra3, "taue670", "Tau ext 670", " ",            &  
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,             &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra3, "taue865", "Tau ext 865", " ",            & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,             &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra3, "taue550_tr2", "Tau ext 550tr2", " ",     & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,             &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra3, "taue670_tr2", "Tau ext 670tr2", " ",     & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,             &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra3, "taue865_tr2", "Tau ext 865tr2", " ",     & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,             &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra3, "taue550_ss", "Tau ext 550ss", " ",       & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,             &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra3, "taue670_ss", "Tau ext 670ss", " ",       & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,             &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra3, "taue865_ss", "Tau ext 865ss", " ",       & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,             &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra3, "taue550_dust", "Tau ext 550dust", " "    & 
                       ,nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,            &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra3, "taue670_dust", "Tau ext 670dust", " "    & 
                       ,nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,            &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra3, "taue865_dust", "Tau ext 865dust", " "    & 
                       ,nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,            &
                       "ave(X)", zsto,zout)                                 
                                                                            
          CALL histdef(nid_tra3, "taue550_dustsco",                     &   
                       "Tau ext 550dustsco", " "                        &
                       ,nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                              
!                                                                        
           CALL histdef(nid_tra3, "taue670_dustsco",                    &
                       "Tau ext 670dustsco", " "                        &
                       ,nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                              
!                                                                        
           CALL histdef(nid_tra3, "taue865_dustsco",                    &
                       "Tau ext 865dustsco", " "                        &
                       ,nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,         &
                       "ave(X)", zsto,zout)                              
                                                                         
                                                                        
        CALL histdef(nid_tra3, "taue550_aqua", "Tau ext 550 aqua", " ",   & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,            &
                       "inst(X)", zout,zout)                               
      CALL histdef(nid_tra3, "taue550_terra", "Tau ext 550 terra", " ",   & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,            & 
                       "inst(X)", zout,zout)                               
        CALL histdef(nid_tra3, "taue670_aqua", "Tau ext 670 aqua", " ",   & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,            & 
                       "inst(X)", zout,zout)                               
      CALL histdef(nid_tra3, "taue670_terra", "Tau ext 670 terra", " ",   & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,            & 
                       "inst(X)", zout,zout)                               
        CALL histdef(nid_tra3, "taue865_aqua", "Tau ext 865 aqua", " ",   & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,            & 
                       "inst(X)", zout,zout)                               
      CALL histdef(nid_tra3, "taue865_terra", "Tau ext 865 terra", " ",   & 
                       nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,            & 
                       "inst(X)", zout,zout)                               
                                                                           
                                                                           
          DO it=1, nbtr
!
          WRITE(str2,'(i2.2)') it
!
          CALL histdef(nid_tra3, "trm"//str2, "Burden No."//str2,         & 
                     "mgS/m2", nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,    &
                       "ave(X)", zsto,zout)                                
!                                                                          
          CALL histdef(nid_tra3, "sconc"//str2, "Surf Conc. No."//str2,   & 
                       "mg/m3", nbp_lon,nbp_lat,nhori3, 1,1,1, -99, 32,   &
                       "ave(X)", zsto,zout)                                
!                                                                          
! LESSIVAGE                                                                 
!
          CALL histdef(nid_tra2, "flux"//str2, "emission"//str2,           &  
                   "mgS/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1, -99, 32,     &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra2, "ds"//str2, "Depot sec No."//str2,        & 
                   "mgS/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1, -99, 32,     &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra2,"dh"//str2,                                 &
                    "Depot hum No."//str2,                                 &
                   "mgS/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1, -99, 32,     &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra2,"dtrconv"//str2,                           &
                     "Tiedke convective"//str2,                            &
                  "mgS/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1, -99, 32,      &
                       "ave(X)", zsto,zout)                                 
                                                                            
          CALL histdef(nid_tra2,"dtherm"//str2,                            &
                       "Thermals dtracer"//str2,                           &
                   "mgS/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1, -99, 32,     &
                       "ave(X)", zsto,zout)                                 
                                                                            
          CALL histdef(nid_tra2,"dhkecv"//str2,                            &
                       "KE dep hum convective"//str2,                      &
                   "mgS/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1, -99, 32,     &
                       "ave(X)", zsto,zout)                                 
          CALL histdef(nid_tra2,"dhkelsc"//str2,                            &
                       "KE dep hum large scale"//str2,                      &
                   "mgS/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1, -99, 32,      &
                       "ave(X)", zsto,zout)                                  
                                                                             
                                
          CALL histdef(nid_tra2,"d_tr_ds"//str2,                            &
                       " Tendance dep sec"//str2,                      &
                   "mgS/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1, -99, 32,     &
                       "ave(X)", zsto,zout)                                 

                                            
          CALL histdef(nid_tra2,"d_tr_cv"//str2,                          & 
                       "cvltr d_tr_cv"//str2,                             &
                       "mgS/m2/s",                                        &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,   &
                       "ave(X)", zsto,zout)                                 
          CALL histdef(nid_tra2,"d_tr_trsp"//str2                         & 
                       ,"cvltr d_tr_trsp"//str2,                          &
                       "mgS/m2/s",                                        &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,   &
                       "ave(X)", zsto,zout)                                
          CALL histdef(nid_tra2,"d_tr_sscav"//str2                        & 
                       ,"cvltr d_tr_sscav"//str2,"mgS/m2/s",                 &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,      &
                       "ave(X)", zsto,zout)                                  
          CALL histdef(nid_tra2,"d_tr_sat"//str2                            &  
                       ,"cvltr d_tr_sat"//str2,                             &  
                       "mgS/m2/s",                                          & 
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,     &
                       "ave(X)", zsto,zout)                                  
        CALL histdef(nid_tra2,"d_tr_uscav"//str2,                           & 
                    "cvltr d_tr_uscav"//str2,                               &
                       "mgS/m2/s",                                          &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,     &
                       "ave(X)", zsto,zout)                                  
        CALL histdef(nid_tra2,"d_tr_insc"//str2,                           &   !!!
                    "cvltr d_tr_insc"//str2,                               &
                       "mgS/m2/s",                                          &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,     &
                       "ave(X)", zsto,zout)                                  
        CALL histdef(nid_tra2,"d_tr_bcscav"//str2,                           & 
                    "cvltr d_tr_bcscav"//str2,                               &
                       "mgS/m2/s",                                          &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,     &
                       "ave(X)", zsto,zout)                                  
        CALL histdef(nid_tra2,"d_tr_evapls"//str2,                           & 
                    "cvltr d_tr_evapls"//str2,                               &
                       "mgS/m2/s",                                          &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,     &
                       "ave(X)", zsto,zout)                                  
        CALL histdef(nid_tra2,"d_tr_ls"//str2,                           & 
                    "cvltr d_tr_ls"//str2,                               &
                       "mgS/m2/s",                                          &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,     &
                       "ave(X)", zsto,zout)                                     !!
        CALL histdef(nid_tra2,"d_tr_dyn"//str2,                           & 
                    "large-scale d_tr_dyn"//str2,                               &
                       "mgS/m2/s",                                          &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,     &
                       "ave(X)", zsto,zout)                                     !!
        CALL histdef(nid_tra2,"d_tr_cl"//str2,                           & 
                    "cvltr d_tr_cl"//str2,                               &
                       "mgS/m2/s",                                          &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,     &
                       "ave(X)", zsto,zout)                                 !!
        CALL histdef(nid_tra2,"d_tr_th"//str2,                           & 
                    "cvltr d_tr_th"//str2,                               &
                       "mgS/m2/s",                                          &
                   nbp_lon,nbp_lat,nhori2, nbp_lev,1,nbp_lev,nvert, 32,     &
                       "ave(X)", zsto,zout)                                 !!
                                                                             


!
          ENDDO
!
          CALL histdef(nid_tra2, "sed_ss", "Sedmet. Tr3",                   & 
                       "mg/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1, -99,       &
                         32,                                                &
                       "ave(X)", zsto,zout)                                  
!                                                                            
          CALL histdef(nid_tra2, "sed_dust", "Sedmet. Tr4",                 &
                       "mg/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1,            &
                        -99, 32,                                            &
                       "ave(X)", zsto,zout)                                  
!                                                                            
          CALL histdef(nid_tra2, "sed_dustsco", "Sedmet. Tr5",              &
                       "mg/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1,            &
                        -99, 32,                                            &
                       "ave(X)", zsto,zout)                                  
!                                                                            
          CALL histdef(nid_tra2, "g2p_gas", "Gas2particle gas sink",       & 
                   "mg-S/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1,-99, 32,     &
                       "ave(X)", zsto,zout)                                 
!                                                                           
          CALL histdef(nid_tra2, "g2p_aer", "Gas2particle tr2 src",        & 
                       "mg/m2/s", nbp_lon,nbp_lat,nhori2, 1,1,1,-99,32,    &
                       "ave(X)", zsto,zout)                                 
!                                                                           
!-------------------------------------------------------------------        
!
          CALL histend(nid_tra1)
!
          CALL histend(nid_tra2)
!
          CALL histend(nid_tra3)
!
!-------------------------------------------------------------------

!       nbjour=1
         ENDIF ! mpi root
         ENDIF !--ok_histrac

!
!        IF (.NOT.edgar.AND.bateau) THEN
!        PRINT *,'ATTENTION risque de compter double les bateaux'
!        STOP
!        ENDIF
!
!
!
!$OMP END MASTER
!$OMP BARRIER
      endif ! debutphy
!
!======================================================================
! Initialisations
!======================================================================
!
!
! je  KE init
      IF (debutphy) THEN
!$OMP MASTER

      ALLOCATE(d_tr_cl(klon,klev,nbtr),d_tr_dry(klon,nbtr))
      ALLOCATE(flux_tr_dry(klon,nbtr),d_tr_dec(klon,klev,nbtr))
      ALLOCATE(d_tr_cv(klon,klev,nbtr))
      ALLOCATE(d_tr_insc(klon,klev,nbtr),d_tr_bcscav(klon,klev,nbtr))
      ALLOCATE(d_tr_evapls(klon,klev,nbtr),d_tr_ls(klon,klev,nbtr))
      ALLOCATE(qPrls(klon,nbtr),d_tr_trsp(klon,klev,nbtr))
      ALLOCATE(d_tr_sscav(klon,klev,nbtr),d_tr_sat(klon,klev,nbtr))
      ALLOCATE(d_tr_uscav(klon,klev,nbtr),qPr(klon,klev,nbtr))
      ALLOCATE(qDi(klon,klev,nbtr))
      ALLOCATE(qPa(klon,klev,nbtr),qMel(klon,klev,nbtr))
      ALLOCATE(qTrdi(klon,klev,nbtr),dtrcvMA(klon,klev,nbtr))
      ALLOCATE(d_tr_th(klon,klev,nbtr))
      ALLOCATE(d_tr_lessi_impa(klon,klev,nbtr))
      ALLOCATE(d_tr_lessi_nucl(klon,klev,nbtr))

      ALLOCATE( diff_aod550_tot(klon)     )
      ALLOCATE( diag_aod670_tot(klon)     )
      ALLOCATE( diag_aod865_tot(klon)     )
      ALLOCATE( diff_aod550_tr2(klon)     )
      ALLOCATE( diag_aod670_tr2(klon)     )
      ALLOCATE( diag_aod865_tr2(klon)     )
      ALLOCATE( diag_aod550_ss(klon)      )
      ALLOCATE( diag_aod670_ss(klon)      )
      ALLOCATE( diag_aod865_ss(klon)      )
      ALLOCATE( diag_aod550_dust(klon)    )
      ALLOCATE( diag_aod670_dust(klon)    )
      ALLOCATE( diag_aod865_dust(klon)    )
      ALLOCATE( diag_aod550_dustsco(klon)  )
      ALLOCATE( diag_aod670_dustsco(klon)  )
      ALLOCATE( diag_aod865_dustsco(klon)  )


      ALLOCATE(  sconc01(klon)     )
      ALLOCATE(  trm01(klon)     )
      ALLOCATE(  sconc02(klon)     )
      ALLOCATE(  trm02(klon)     )
      ALLOCATE(  sconc03(klon)     )
      ALLOCATE(  trm03(klon)     )
      ALLOCATE(  sconc04(klon)     )
      ALLOCATE(  trm04(klon)     )
      ALLOCATE(  sconc05(klon)     )
      ALLOCATE(  trm05(klon)     )


      ALLOCATE(  flux01(klon)     )
      ALLOCATE(  flux02(klon)     )
      ALLOCATE(  flux03(klon)     )
      ALLOCATE(  flux04(klon)     )
      ALLOCATE(  flux05(klon)     )
      ALLOCATE(  ds01(klon)     )
      ALLOCATE(  ds02(klon)     )
      ALLOCATE(  ds03(klon)     )
      ALLOCATE(  ds04(klon)     )
      ALLOCATE(  ds05(klon)     )
      ALLOCATE(  dh01(klon)     )
      ALLOCATE(  dh02(klon)     )
      ALLOCATE(  dh03(klon)     )
      ALLOCATE(  dh04(klon)     )
      ALLOCATE(  dh05(klon)     )
      ALLOCATE(  dtrconv01(klon)     )
      ALLOCATE(  dtrconv02(klon)     )
      ALLOCATE(  dtrconv03(klon)     )
      ALLOCATE(  dtrconv04(klon)     )
      ALLOCATE(  dtrconv05(klon)     )
      ALLOCATE(  dtherm01(klon)     )
      ALLOCATE(  dtherm02(klon)     )
      ALLOCATE(  dtherm03(klon)     )
      ALLOCATE(  dtherm04(klon)     )
      ALLOCATE(  dtherm05(klon)     )
      ALLOCATE(  dhkecv01(klon)     )
      ALLOCATE(  dhkecv02(klon)     )
      ALLOCATE(  dhkecv03(klon)     )
      ALLOCATE(  dhkecv04(klon)     )
      ALLOCATE(  dhkecv05(klon)     )
      ALLOCATE(  d_tr_ds01(klon)     )
      ALLOCATE(  d_tr_ds02(klon)     )
      ALLOCATE(  d_tr_ds03(klon)     )
      ALLOCATE(  d_tr_ds04(klon)     )
      ALLOCATE(  d_tr_ds05(klon)     )
      ALLOCATE(  dhkelsc01(klon)     )
      ALLOCATE(  dhkelsc02(klon)     )
      ALLOCATE(  dhkelsc03(klon)     )
      ALLOCATE(  dhkelsc04(klon)     )
      ALLOCATE(  dhkelsc05(klon)     )
      ALLOCATE(  d_tr_cv01(klon,klev))
      ALLOCATE(  d_tr_cv02(klon,klev))
      ALLOCATE(  d_tr_cv03(klon,klev))
      ALLOCATE(  d_tr_cv04(klon,klev))
      ALLOCATE(  d_tr_cv05(klon,klev))
      ALLOCATE(  d_tr_trsp01(klon,klev))
      ALLOCATE(  d_tr_trsp02(klon,klev))
      ALLOCATE(  d_tr_trsp03(klon,klev))
      ALLOCATE(  d_tr_trsp04(klon,klev))
      ALLOCATE(  d_tr_trsp05(klon,klev))
      ALLOCATE(  d_tr_sscav01(klon,klev))
      ALLOCATE(  d_tr_sscav02(klon,klev))
      ALLOCATE(  d_tr_sscav03(klon,klev))
      ALLOCATE(  d_tr_sscav04(klon,klev))
      ALLOCATE(  d_tr_sscav05(klon,klev))
      ALLOCATE(  d_tr_sat01(klon,klev))
      ALLOCATE(  d_tr_sat02(klon,klev))
      ALLOCATE(  d_tr_sat03(klon,klev))
      ALLOCATE(  d_tr_sat04(klon,klev))
      ALLOCATE(  d_tr_sat05(klon,klev))
      ALLOCATE(  d_tr_uscav01(klon,klev))
      ALLOCATE(  d_tr_uscav02(klon,klev))
      ALLOCATE(  d_tr_uscav03(klon,klev))
      ALLOCATE(  d_tr_uscav04(klon,klev))
      ALLOCATE(  d_tr_uscav05(klon,klev))
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(  d_tr_insc01(klon,klev))
      ALLOCATE(  d_tr_insc02(klon,klev))
      ALLOCATE(  d_tr_insc03(klon,klev))
      ALLOCATE(  d_tr_insc04(klon,klev))
      ALLOCATE(  d_tr_insc05(klon,klev))
      ALLOCATE(  d_tr_bcscav01(klon,klev))
      ALLOCATE(  d_tr_bcscav02(klon,klev))
      ALLOCATE(  d_tr_bcscav03(klon,klev))
      ALLOCATE(  d_tr_bcscav04(klon,klev))
      ALLOCATE(  d_tr_bcscav05(klon,klev))
      ALLOCATE(  d_tr_evapls01(klon,klev))
      ALLOCATE(  d_tr_evapls02(klon,klev))
      ALLOCATE(  d_tr_evapls03(klon,klev))
      ALLOCATE(  d_tr_evapls04(klon,klev))
      ALLOCATE(  d_tr_evapls05(klon,klev))
      ALLOCATE(  d_tr_ls01(klon,klev))
      ALLOCATE(  d_tr_ls02(klon,klev))
      ALLOCATE(  d_tr_ls03(klon,klev))
      ALLOCATE(  d_tr_ls04(klon,klev))
      ALLOCATE(  d_tr_ls05(klon,klev))
      ALLOCATE(  d_tr_dyn01(klon,klev))
      ALLOCATE(  d_tr_dyn02(klon,klev))
      ALLOCATE(  d_tr_dyn03(klon,klev))
      ALLOCATE(  d_tr_dyn04(klon,klev))
      ALLOCATE(  d_tr_dyn05(klon,klev))
      ALLOCATE(  d_tr_cl01(klon,klev))
      ALLOCATE(  d_tr_cl02(klon,klev))
      ALLOCATE(  d_tr_cl03(klon,klev))
      ALLOCATE(  d_tr_cl04(klon,klev))
      ALLOCATE(  d_tr_cl05(klon,klev))
      ALLOCATE(  d_tr_th01(klon,klev))
      ALLOCATE(  d_tr_th02(klon,klev))
      ALLOCATE(  d_tr_th03(klon,klev))
      ALLOCATE(  d_tr_th04(klon,klev))
      ALLOCATE(  d_tr_th05(klon,klev))

      ALLOCATE( sed_ss3D(klon,klev))
      ALLOCATE( sed_dust3D(klon,klev))
      ALLOCATE( sed_dustsco3D(klon,klev))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE( sed_ss(klon))
      ALLOCATE( sed_dust(klon))
      ALLOCATE( sed_dustsco(klon))
      ALLOCATE( his_g2pgas(klon))
      ALLOCATE( his_g2paer(klon))

      ALLOCATE( fluxbb(klon))
      ALLOCATE( fluxff(klon))
      ALLOCATE( fluxbcbb(klon))
      ALLOCATE( fluxbcff(klon))
      ALLOCATE( fluxbcnff(klon))
      ALLOCATE( fluxbcba(klon))
      ALLOCATE( fluxbc(klon))
      ALLOCATE( fluxombb(klon))
      ALLOCATE( fluxomff(klon))
      ALLOCATE( fluxomnff(klon))
      ALLOCATE( fluxomba(klon))
      ALLOCATE( fluxomnat(klon))
      ALLOCATE( fluxom(klon))
      ALLOCATE( fluxh2sff(klon))
      ALLOCATE( fluxh2snff(klon))
      ALLOCATE( fluxso2ff(klon))
      ALLOCATE( fluxso2nff(klon))
      ALLOCATE( fluxso2bb(klon))
      ALLOCATE( fluxso2vol(klon))
      ALLOCATE( fluxso2ba(klon))
      ALLOCATE( fluxso2(klon))
      ALLOCATE( fluxso4ff(klon))
      ALLOCATE( fluxso4nff(klon))
      ALLOCATE( fluxso4bb(klon))
      ALLOCATE( fluxso4ba(klon))
      ALLOCATE( fluxso4(klon))
      ALLOCATE( fluxdms(klon))
      ALLOCATE( fluxh2sbio(klon))
      ALLOCATE( fluxdustec(klon))
      ALLOCATE( fluxddfine(klon))
      ALLOCATE( fluxddcoa(klon))
      ALLOCATE( fluxddsco(klon))
      ALLOCATE( fluxdd(klon))
      ALLOCATE( fluxssfine(klon))
      ALLOCATE( fluxsscoa(klon))
      ALLOCATE( fluxss(klon))
      ALLOCATE( flux_sparam_ind(klon))
      ALLOCATE( flux_sparam_bb(klon))
      ALLOCATE( flux_sparam_ff(klon))
      ALLOCATE( flux_sparam_ddfine(klon))
      ALLOCATE( flux_sparam_ddcoa(klon))
      ALLOCATE( flux_sparam_ddsco(klon))
      ALLOCATE( flux_sparam_ssfine(klon))
      ALLOCATE( flux_sparam_sscoa(klon))
      ALLOCATE( u10m_ss(klon))
      ALLOCATE( v10m_ss(klon))


       ALLOCATE(d_tr_cv_o(klon,klev,nbtr))
       ALLOCATE(d_tr_trsp_o(klon,klev,nbtr))
       ALLOCATE(d_tr_sscav_o(klon,klev,nbtr), &
                d_tr_sat_o(klon,klev,nbtr))
        ALLOCATE(d_tr_uscav_o(klon,klev,nbtr))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ALLOCATE(d_tr_insc_o(klon,klev,nbtr))
        ALLOCATE(d_tr_bcscav_o(klon,klev,nbtr))
        ALLOCATE(d_tr_evapls_o(klon,klev,nbtr))
        ALLOCATE(d_tr_ls_o(klon,klev,nbtr))
        ALLOCATE(d_tr_dyn_o(klon,klev,nbtr))
        ALLOCATE(d_tr_cl_o(klon,klev,nbtr))
        ALLOCATE(d_tr_th_o(klon,klev,nbtr))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ALLOCATE(iregion_so4(klon))
        ALLOCATE(iregion_bb(klon))
        ALLOCATE(iregion_ind(klon))
        ALLOCATE(iregion_dust(klon))
        ALLOCATE(iregion_wstardust(klon))

!JE20150518<<
       ALLOCATE(masque_aqua(klon))  
       ALLOCATE(masque_terra(klon))  
!       ALLOCATE(aod550_aqua(klon))  
!       ALLOCATE(aod550_terra(klon))  
!       ALLOCATE(aod670_aqua(klon))  
!       ALLOCATE(aod670_terra(klon))  
!       ALLOCATE(aod865_aqua(klon))  
!       ALLOCATE(aod865_terra(klon)) 

      ALLOCATE( aod550_terra(klon))  
      ALLOCATE( aod550_tr2_terra(klon))  
      ALLOCATE( aod550_ss_terra(klon))   
      ALLOCATE( aod550_dust_terra(klon))   
      ALLOCATE( aod550_dustsco_terra(klon))   
      ALLOCATE( aod670_terra(klon))   
      ALLOCATE( aod670_tr2_terra(klon))  
      ALLOCATE( aod670_ss_terra(klon))  
      ALLOCATE( aod670_dust_terra(klon))  
      ALLOCATE( aod670_dustsco_terra(klon))  
      ALLOCATE( aod865_terra(klon))   
      ALLOCATE( aod865_tr2_terra(klon))  
      ALLOCATE( aod865_ss_terra(klon))  
      ALLOCATE( aod865_dust_terra(klon))  
      ALLOCATE( aod865_dustsco_terra(klon))  

      ALLOCATE( aod550_aqua(klon))  
      ALLOCATE( aod550_tr2_aqua(klon))  
      ALLOCATE( aod550_ss_aqua(klon))   
      ALLOCATE( aod550_dust_aqua(klon))   
      ALLOCATE( aod550_dustsco_aqua(klon))   
      ALLOCATE( aod670_aqua(klon))   
      ALLOCATE( aod670_tr2_aqua(klon))  
      ALLOCATE( aod670_ss_aqua(klon))  
      ALLOCATE( aod670_dust_aqua(klon))  
      ALLOCATE( aod670_dustsco_aqua(klon))  
      ALLOCATE( aod865_aqua(klon))   
      ALLOCATE( aod865_tr2_aqua(klon))  
      ALLOCATE( aod865_ss_aqua(klon))  
      ALLOCATE( aod865_dust_aqua(klon))  
      ALLOCATE( aod865_dustsco_aqua(klon))  
 

       masque_aqua(:)=0
       masque_terra(:)=0
!       aod550_aqua(:)=0.
!       aod550_terra(:)=0.
!       aod670_aqua(:)=0.
!       aod670_terra(:)=0.
!       aod865_aqua(:)=0.
!       aod865_terra(:)=0.

      aod550_terra(:)=0.  
      aod550_tr2_terra(:)=0.  
      aod550_ss_terra(:)=0.   
      aod550_dust_terra(:)=0.   
      aod550_dustsco_terra(:)=0.   
      aod670_terra(:)=0.   
      aod670_tr2_terra(:)=0.  
      aod670_ss_terra(:)=0.  
      aod670_dust_terra(:)=0.  
      aod670_dustsco_terra(:)=0.  
      aod865_terra(:)=0.   
      aod865_tr2_terra(:)=0.  
      aod865_ss_terra(:)=0.  
      aod865_dust_terra(:)=0.  
      aod865_dustsco_terra(:)=0.  
      aod550_aqua(:)=0.  
      aod550_tr2_aqua(:)=0.  
      aod550_ss_aqua(:)=0.   
      aod550_dust_aqua(:)=0.   
      aod550_dustsco_aqua(:)=0.   
      aod670_aqua(:)=0.   
      aod670_tr2_aqua(:)=0.  
      aod670_ss_aqua(:)=0.  
      aod670_dust_aqua(:)=0.  
      aod670_dustsco_aqua(:)=0.  
      aod865_aqua(:)=0.   
      aod865_tr2_aqua(:)=0.  
      aod865_ss_aqua(:)=0.  
      aod865_dust_aqua(:)=0.  
      aod865_dustsco_aqua(:)=0.  
!JE20150518>>





!
!Config Key  = iflag_lscav
!Config Desc = Large scale scavenging parametrization: 0=none,
!1=old(Genthon92),
!              2=1+PHeinrich, 3=Reddy_Boucher2004, 4=3+RPilon.
!Config Def  = 4
!Config
        iflag_lscav_omp=4
        call getin('iflag_lscav', iflag_lscav_omp)
        iflag_lscav=iflag_lscav_omp
! initialiation for time computation

        tia_spla=0.
        tia_emis=0.
        tia_depo=0.
        tia_cltr=0.
        tia_ther=0.
        tia_sedi=0.
        tia_gasp=0.
        tia_wetap=0.
        tia_cvltr=0.
        tia_lscs=0.
        tia_brop=0.
        tia_outs=0.
        tia_nophytracr=0.
        clock_start_outphytracr=clock_end_outphytracr+1
!$OMP END MASTER
!$OMP BARRIER
       ENDIF ! debutphy
     
      lmt_dms(:)=0.0
      aux_var2(:)=0.0
      aux_var3(:,:)=0.0
      source_tr(:,:)=0.0
      flux_tr(:,:)=0.0
      flux_sparam_bb(:)=0.0
      flux_sparam_ff(:)=0.0
      flux_sparam_ind(:)=0.0
      flux_sparam_ddfine(:)=0.0
      flux_sparam_ddcoa(:)=0.0
      flux_sparam_ddsco(:)=0.0
      flux_sparam_ssfine(:)=0.0
      flux_sparam_sscoa(:)=0.0

! initialiation for time computation
        
        ti_spla=0
        ti_emis=0
        ti_depo=0
        ti_cltr=0
        ti_ther=0
        ti_sedi=0
        ti_gasp=0
        ti_wetap=0
        ti_cvltr=0
        ti_lscs=0
        ti_brop=0
        ti_outs=0


       DO k=1,klev
        DO i=1,klon
         Mint(i,k)=0.
        END DO
       END DO


!
      DO it=1, nbtr
       DO k=1,klev
        DO i=1,klon
         d_tr_cv(i,k,it)=0.
         d_tr_trsp(i,k,it)=0.
         d_tr_sscav(i,k,it)=0.
         d_tr_sat(i,k,it)=0.
         d_tr_uscav(i,k,it)=0.
         d_tr(i,k,it)=0.
         d_tr_insc(i,k,it)=0.
         d_tr_bcscav(i,k,it)=0.
         d_tr_evapls(i,k,it)=0.
         d_tr_ls(i,k,it)=0.
         d_tr_cl(i,k,it)=0.
         d_tr_th(i,k,it)=0.
  
         d_tr_cv_o(i,k,it)=0.
         d_tr_trsp_o(i,k,it)=0.
         d_tr_sscav_o(i,k,it)=0.
         d_tr_sat_o(i,k,it)=0.
         d_tr_uscav_o(i,k,it)=0.


         qDi(i,k,it)=0.
         qPr(i,k,it)=0.
         qPa(i,k,it)=0.
         qMel(i,k,it)=0.
         qTrdi(i,k,it)=0.
         dtrcvMA(i,k,it)=0.
         zmfd1a(i,k,it)=0.
         zmfdam(i,k,it)=0.
         zmfphi2(i,k,it)=0.
        END DO
       END DO
      END DO


      DO it=1, nbtr
       DO i=1,klon
          qPrls(i,it)=0.0
          dtrconv(i,it)=0.0
!JE20140507<<
          d_tr_dry(i,it)=0.0
          flux_tr_dry(i,it)=0.0
!JE20140507>>
       ENDDO
      ENDDO

      DO it=1, nbtr
      DO i=1, klon 
        his_dh(i,it)=0.0
        his_dhlsc(i,it)=0.0
        his_dhcon(i,it)=0.0
        his_dhbclsc(i,it)=0.0
        his_dhbccon(i,it)=0.0
        trm(i,it)=0.0
        his_th(i,it)=0.0
        his_dhkecv(i,it)=0.0
        his_ds(i,it)=0.0
        his_dhkelsc(i,it)=0.0

      ENDDO
      ENDDO
!JE:      
      DO i=1, klon 
         his_g2pgas(i) = 0.0
         his_g2paer(i) = 0.0
      ENDDO
! endJE
!

      DO k=1, klev
      DO i = 1, klon
        zrho(i,k)=pplay(i,k)/t_seri(i,k)/RD
        zdz(i,k)=(paprs(i,k)-paprs(i,k+1))/zrho(i,k)/RG
        zmasse(i,k)=(paprs(i,k)-paprs(i,k+1))/RG
      ENDDO
      ENDDO
!
      DO i = 1, klon
        zalt(i,1)=pphis(i)/RG
      ENDDO
      DO k=1, klev-1
      DO i = 1, klon
        zalt(i,k+1)=zalt(i,k)+zdz(i,k)
      ENDDO
      ENDDO



      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_init=dife*MAX(0,SIGN(1,dife)) &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_init=tia_init+REAL(ti_init)/REAL(clock_rate)
      ENDIF
      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF



       IF (debutphy) then

      c_FullName1='regions_dustacc'
      !c_FullName1='regions_dust'
      call readregions_spl(iregion_dust,c_FullName1)
      c_FullName1='regions_ind'
      call readregions_spl(iregion_ind,c_FullName1)
      c_FullName1='regions_bb'
      call readregions_spl(iregion_bb,c_FullName1)
      c_FullName1='regions_pwstarwake'
      call readregions_spl(iregion_wstardust,c_FullName1)

!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     
      OPEN(25,FILE='dustregions_pyvar_je.data')
      OPEN(55,FILE='indregions_pyvar_je.data')
      OPEN(75,FILE='bbregions_pyvar_je.data')
      OPEN(95,FILE='wstardustregions_pyvar_je.data')
      OPEN(76,FILE='xlat.data')
      OPEN(77,FILE='xlon.data')
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER

      CALL gather(iregion_dust,iauxklon_glo)
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      DO k=1,klon_glo
        WRITE(25,'(i10)') iauxklon_glo(k)
      ENDDO
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
      CALL gather(iregion_ind,iauxklon_glo)
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      DO k=1,klon_glo
        WRITE(55,'(i10)') iauxklon_glo(k)
      ENDDO
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
      CALL gather(iregion_bb,iauxklon_glo)
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      DO k=1,klon_glo
        WRITE(75,'(i10)') iauxklon_glo(k)
      ENDDO
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
      CALL gather(iregion_wstardust,iauxklon_glo)
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      DO k=1,klon_glo
        WRITE(95,'(i10)') iauxklon_glo(k)
      ENDDO
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER


      CALL gather(rlat,auxklon_glo)
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      DO k=1,klon_glo
        WRITE(76,*) auxklon_glo(k)
      ENDDO
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
      CALL gather(rlon,auxklon_glo)
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      DO k=1,klon_glo
        WRITE(77,*) auxklon_glo(k)
      ENDDO

      CLOSE(25)
      CLOSE(55)
      CLOSE(75)
      CLOSE(76)
      CLOSE(77)

      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER

      ENDIF  ! debutphy

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_inittype=dife*MAX(0,SIGN(1,dife)) &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_inittype=tia_inittype+REAL(ti_inittype)/REAL(clock_rate)
      ENDIF

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF

!
!=======================================================================
! SAVING SURFACE TYPE
!=======================================================================
      IF (debutphy) THEN
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN

      OPEN(35,FILE='surface_ocean.data')
      OPEN(45,FILE='surface_seaice.data')
      OPEN(65,FILE='surface_land.data')
      OPEN(85,FILE='surface_landice.data')
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
      do i = 1, klon
                aux_var2(i) = pctsrf(i,is_oce) 
      enddo
      call gather(aux_var2,auxklon_glo)
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      DO i = 1, klon_glo
         WRITE (35,103)  auxklon_glo(i)
      ENDDO
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER

      do i = 1, klon
                aux_var2(i) = pctsrf(i,is_sic) 
      enddo
      call gather(aux_var2,auxklon_glo)
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      DO i = 1, klon_glo
         WRITE (45,103)  auxklon_glo(i)
      ENDDO
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER

      do i = 1, klon
                aux_var2(i) = pctsrf(i,is_ter) 
      enddo
      call gather(aux_var2,auxklon_glo)
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      DO i = 1, klon_glo
         WRITE (65,103)  auxklon_glo(i)
      ENDDO
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER

      do i = 1, klon
                aux_var2(i) = pctsrf(i,is_lic) 
      enddo
      call gather(aux_var2,auxklon_glo)
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      DO i = 1, klon_glo
         WRITE (85,103)  auxklon_glo(i)
      ENDDO
!
!      DO i = 1, klon
!         WRITE (35,103) pctsrf(i,is_oce)
!         WRITE (45,103) pctsrf(i,is_sic)
!         WRITE (65,103) pctsrf(i,is_ter)
!         WRITE (85,103) pctsrf(i,is_lic)
!      ENDDO
      CLOSE(35)
      CLOSE(45)
      CLOSE(65)
      CLOSE(85)
103   FORMAT (f6.2)
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
      ENDIF ! debutphy

!      stop
!
!=======================================================================
!
      DO it=1, nbtr
        DO j=1,klev
        DO i=1,klon
           tmp_var(i,j)=tr_seri(i,j,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
        DO j=1,klev
        DO i=1,klon
           tr_seri(i,j,it)=tmp_var(i,j)
        ENDDO
        ENDDO
      ENDDO
      iscm3=.true.

!=======================================================================
!
      DO k=1, klev
      DO i=1, klon
        m_conc(i,k)=pplay(i,k)/t_seri(i,k)/RKBOL*1.e-6
      ENDDO
      ENDDO

!
!
      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_avt_coarem')
        ENDDO        
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'avt coarem')
        ENDDO
        CALL minmaxsource(source_tr,qmin,qmax,'src: avt coarem')
      ENDIF

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_inittwrite=dife*MAX(0,SIGN(1,dife))  &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_inittwrite=tia_inittwrite+REAL(ti_inittwrite)/REAL(clock_rate)
      ENDIF

!
!
!=======================================================================
!                     EMISSIONS OF COARSE AEROSOLS
!=======================================================================


      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF



!      
      print *,'Number of tracers = ',nbtr

      print *,'AT BEGINNING OF PHYTRACR_SPL'
!      print *,'tr_seri = ',SUM(tr_seri(:,:,3)),MINVAL(tr_seri(:,:,3)),
!     .                                         MAXVAL(tr_seri(:,:,3))
#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('sav'//str2,1,'SOURCE','',source_tr(:,it))
         call iophys_ecrit('fav'//str2,1,'SOURCE','',source_tr(:,it))
      enddo
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRB'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif


      CALL coarsemission(pctsrf,pdtphys,t_seri,                            &
                        pmflxr,pmflxs,prfl,psfl,                           &
                        rlat,rlon,debutphy,                                & 
                        zu10m,zv10m,wstar,ale_bl,ale_wake,                 &
                        scale_param_ssacc,scale_param_sscoa,               & 
                        scale_param_dustacc,scale_param_dustcoa,           &
                        scale_param_dustsco,                               & 
                        nbreg_dust,                                        &
                        iregion_dust,dust_ec,                              & 
                        param_wstarBLperregion,param_wstarWAKEperregion,   &
                        nbreg_wstardust,                                   &
                        iregion_wstardust,                                 &
                        lmt_sea_salt,qmin,qmax,                            &
                                  flux_sparam_ddfine,flux_sparam_ddcoa,    & 
                                  flux_sparam_ddsco,                       &
                                  flux_sparam_ssfine,flux_sparam_sscoa,    & 
                              id_prec,id_fine,id_coss,id_codu,id_scdu,     &
                              ok_chimeredust,                           &
                                                     source_tr,flux_tr)    

#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('sap'//str2,1,'SOURCE','',source_tr(:,it))
         call iophys_ecrit('fap'//str2,1,'SOURCE','',source_tr(:,it))
      enddo
#endif

      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_after_coarem')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after coarem')
        ENDDO
        CALL minmaxsource(source_tr,qmin,qmax,'src: after coarem')
      ENDIF

!
!
!
!======================================================================
!                   EMISSIONS OF AEROSOL PRECURSORS      
!======================================================================
!
#ifdef IOPHYS_DUST
      print *,'INPUT TO PRECUREMISSION'
         call iophys_ecrit('ftsol',4,'ftsol','',ftsol)
         call iophys_ecrit('u10m_ec',1,'u10m_ec','',u10m_ec)
         call iophys_ecrit('v10m_ec',1,'v10m_ec','',v10m_ec)
         call iophys_ecrit('pctsrf',4,'pctsrf','',pctsrf)
         call iophys_ecrit('u_seri',klev,'u_seri','',u_seri)
         call iophys_ecrit('v_seri',klev,'v_seri','',v_seri)
         call iophys_ecrit('paprs',klev,'paprs','',paprs)
         call iophys_ecrit('pplay',klev,'pplay','',pplay)
         call iophys_ecrit('cdragh',1,'cdragh','',cdragh)
         call iophys_ecrit('cdragm',1,'cdragm','',cdragm)
         call iophys_ecrit('t_seri',klev,'t_seri','',t_seri)
         call iophys_ecrit('q_seri',klev,'q_seri','',q_seri)
         call iophys_ecrit('tsol',1,'tsol','',tsol)
         print*,'fracso2emis,frach2sofso2,bateau',fracso2emis,frach2sofso2,bateau
         print*,'kminbc,kmaxbc,pdtphys',kminbc,kmaxbc,pdtphys
         print*,'scale_param_bb,scale_param_ind',scale_param_bb,scale_param_ind
         print*,'iregion_ind,iregion_bb,nbreg_ind, nbreg_bb',iregion_ind,iregion_bb,nbreg_ind, nbreg_bb
         print*,'id_prec,id_fine',id_prec,id_fine
         call iophys_ecrit('zdz',klev,'zdz','',zdz)
         call iophys_ecrit('zalt',klev,'zalt','',zalt)
         call iophys_ecrit('lmt_so2ff_l',1,'lmt_so2ff_l','',lmt_so2ff_l)
         call iophys_ecrit('lmt_so2ff_h',1,'lmt_so2ff_h','',lmt_so2ff_h)
         call iophys_ecrit('lmt_so2nff',1,'lmt_so2nff','',lmt_so2nff)
         call iophys_ecrit('lmt_so2ba',1,'lmt_so2ba','',lmt_so2ba)
         call iophys_ecrit('lmt_so2bb_l',1,'lmt_so2bb_l','',lmt_so2bb_l)
         call iophys_ecrit('lmt_so2bb_h',1,'lmt_so2bb_h','',lmt_so2bb_h)
         call iophys_ecrit('lmt_so2volc_cont',1,'lmt_so2volc_cont','',lmt_so2volc_cont)
         call iophys_ecrit('lmt_altvolc_cont',1,'lmt_altvolc_cont','',lmt_altvolc_cont)
         call iophys_ecrit('lmt_so2volc_expl',1,'lmt_so2volc_expl','',lmt_so2volc_expl)
         call iophys_ecrit('lmt_altvolc_expl',1,'lmt_altvolc_expl','',lmt_altvolc_expl)
         call iophys_ecrit('lmt_dmsbio',1,'lmt_dmsbio','',lmt_dmsbio)
         call iophys_ecrit('lmt_h2sbio',1,'lmt_h2sbio','',lmt_h2sbio)
         call iophys_ecrit('lmt_dmsconc',1,'lmt_dmsconc','',lmt_dmsconc)
         call iophys_ecrit('lmt_dms',1,'lmt_dms','',lmt_dms)
         call iophys_ecrit('flux_sparam_ind',1,'flux_sparam_ind','',flux_sparam_ind)
         call iophys_ecrit('flux_sparam_bb',1,'flux_sparam_bb','',flux_sparam_bb)
#endif



     print*,'ON PASSE DANS precuremission'
     CALL precuremission(ftsol,u10m_ec,v10m_ec,pctsrf,                  &
                         u_seri,v_seri,paprs,pplay,cdragh,cdragm,       &
                         t_seri,q_seri,tsol,fracso2emis,frach2sofso2,   &
                         bateau,zdz,zalt,kminbc,kmaxbc,pdtphys,         &
                         scale_param_bb,scale_param_ind,                &
                         iregion_ind, iregion_bb,                       &
                         nbreg_ind, nbreg_bb,                           &
                         lmt_so2ff_l,lmt_so2ff_h, lmt_so2nff,lmt_so2ba, &
                         lmt_so2bb_l,lmt_so2bb_h,                       &
                         lmt_so2volc_cont,lmt_altvolc_cont,             &
                         lmt_so2volc_expl,lmt_altvolc_expl,             &
                         lmt_dmsbio,lmt_h2sbio, lmt_dmsconc, lmt_dms,   &
                         id_prec,id_fine,                               &
                                       flux_sparam_ind, flux_sparam_bb, &
                                       source_tr,flux_tr,tr_seri)       
!
      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_after precur')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after precur')
        ENDDO
        CALL minmaxsource(source_tr,qmin,qmax,'src: after precur')
      ENDIF

!=======================================================================
!                      EMISSIONS OF FINE AEROSOLS
!=======================================================================
#ifdef IOPHYS_DUST
!
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('tpr'//str2,1,'SOURCE','',source_tr(:,it))
         call iophys_ecrit('fpr'//str2,1,'SOURCE','',flux_tr(:,it))
      enddo
#endif

      CALL finemission(zdz,pdtphys,zalt,kminbc,kmaxbc,                     & 
                      scale_param_bb,scale_param_ff,                       &
                      iregion_ind,iregion_bb,                              & 
                      nbreg_ind,nbreg_bb,                                  & 
                      lmt_bcff, lmt_bcnff, lmt_bcbb_l,lmt_bcbb_h,          & 
                      lmt_bcba, lmt_omff, lmt_omnff,                       & 
                      lmt_ombb_l, lmt_ombb_h, lmt_omnat, lmt_omba,         & 
                      id_fine,                                             & 
                                       flux_sparam_bb, flux_sparam_ff,     & 
                                             source_tr,flux_tr,tr_seri)     
!
!
      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_after_fineem')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after fineem')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,   &
           pplay,t_seri,iscm3,'after fineem')                  
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: after fineem')
      ENDIF

!

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_emis=dife*MAX(0,SIGN(1,dife))   &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_emis=tia_emis+REAL(ti_emis)/REAL(clock_rate)
      ENDIF


#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('t'//str2,1,'SOURCE','',source_tr(:,it))
         call iophys_ecrit('f'//str2,1,'SOURCE','',flux_tr(:,it))
      enddo
#endif
!
!



!
!=======================================================================
!                 DRY DEPOSITION AND BOUNDARY LAYER MIXING
!=======================================================================
!
!        DO it=1,nbtr
!         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,
!     .      pplay,t_seri,iscm3,'')
!        ENDDO

!======================================================================
!    -- Dry deposition --
!======================================================================
      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF

      DO it=1, nbtr
         DO j=1,klev
         DO i=1,klon
           tmp_var(i,j)=tr_seri(i,j,it)
         ENDDO 
         ENDDO
         CALL cm3_to_kg(pplay,t_seri,tmp_var)
         DO j=1,klev
         DO i=1,klon
           tr_seri(i,j,it)=tmp_var(i,j)
         ENDDO 
         ENDDO
      ENDDO
      iscm3=.false.
!----------------------------
      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_before_depo')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'before depo')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz, &
           pplay,t_seri,iscm3,'before depo')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: before depo')
      ENDIF

#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRC'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif

      CALL deposition(vdep_oce,vdep_sic,vdep_ter,vdep_lic,pctsrf,      &
                     zrho,zdz,pdtphys,RHcl,masse,t_seri,pplay,paprs,  &
                     lminmax,qmin,qmax,                               &
                              his_ds,source_tr,tr_seri)
!
      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_after_depo')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after depo')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,   &
           pplay,t_seri,iscm3,'after depo')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: after depo')
      ENDIF

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_depo=dife*MAX(0,SIGN(1,dife))                      &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_depo=tia_depo+REAL(ti_depo)/REAL(clock_rate)
      ENDIF


!
!======================================================================
!    -- Boundary layer mixing --
!======================================================================

#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRD'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif



      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF

!

       DO k = 1, klev
        DO i = 1, klon
         delp(i,k) = paprs(i,k)-paprs(i,k+1)
        END DO
      END DO
!
      DO it=1, nbtr
      DO j=1, klev
      DO i=1, klon
        tmp_var(i,j)=tr_seri(i,j,it)
        aux_var2(i)=source_tr(i,it)
      ENDDO
      ENDDO
      IF (iflag_conv.EQ.2) THEN
! Tiedke
      CALL cltrac_spl(pdtphys,coefh,yu1,yv1,t_seri,tmp_var,  &
                 aux_var2,paprs,pplay,aux_var3)

      ELSE IF (iflag_conv.GE.3) THEN
!KE
      CALL cltrac(pdtphys, coefh,t_seri,tmp_var,aux_var2,paprs,pplay,  &
                 delp,aux_var3,d_tr_dry,flux_tr_dry(:,it))
      ENDIF

      DO i=1, klon
      DO j=1, klev
        tr_seri(i,j,it)=tmp_var(i,j)
        d_tr(i,j,it)=aux_var3(i,j)
        d_tr_cl(i,j,it)=d_tr(i,j,it)
      ENDDO
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         tr_seri(i,k,it) = tr_seri(i,k,it) + d_tr(i,k,it)
      ENDDO
      ENDDO
      print *,' AFTER Cltrac'
      IF (lminmax) THEN
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after cltrac')
      ENDIF
      ENDDO !--end itr loop

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_cltr=dife*MAX(0,SIGN(1,dife))     &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_cltr=tia_cltr+REAL(ti_cltr)/REAL(clock_rate)
      ENDIF



!======================================================================
!    -- Calcul de l'effet des thermiques for KE--
!======================================================================

#ifdef IOPHYS_DUST
      print*,'iflag_conv=',iflag_conv
      call iophys_ecrit('coefh',klev,'coefh','',coefh)
      call iophys_ecrit('yu1',1,'yu1','',yu1)
      call iophys_ecrit('yv1',1,'yv1','',yv1)
      call iophys_ecrit('delp',klev,'delp','',delp)
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRE'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif



      IF (iflag_conv.GE.3) THEN

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF




      
       IF (lminmax) THEN
        DO it=1,nbtr
       CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_before therm')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'before therm')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,  &
           pplay,t_seri,iscm3,'before therm')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'before therm')
      ENDIF

      DO it=1,nbtr
         DO k=1,klev
            DO i=1,klon
               tmp_var3(i,k,it)=tr_seri(i,k,it)
               d_tr_th(i,k,it)=0.
               tr_seri(i,k,it)=MAX(tr_seri(i,k,it),0.)
!JE: precursor >>1e10         tr_seri(i,k,it)=MIN(tr_seri(i,k,it),1.e10)
            END DO
         END DO
      END DO

!JE  new implicit scheme 20140323
      DO it=1,nbtr
        CALL thermcell_dq(klon,klev,1,pdtphys,fm_therm,entr_therm,  &
                         zmasse,tr_seri(1:klon,1:klev,it),         &
                         d_tr(1:klon,1:klev,it),ztra_th,0 )

        DO k=1,klev
           DO i=1,klon
              d_tr(i,k,it)=pdtphys*d_tr(i,k,it)
              d_tr_th(i,k,it)=d_tr_th(i,k,it)+d_tr(i,k,it)
              tr_seri(i,k,it)=MAX(tr_seri(i,k,it)+d_tr(i,k,it),0.)
              END DO
        END DO

      ENDDO

! old scheme explicit
!       nsplit=10
!       DO it=1,nbtr
!          DO isplit=1,nsplit
!              CALL dqthermcell(klon,klev,pdtphys/nsplit, 
!     .            fm_therm,entr_therm,zmasse, 
!     .            tr_seri(1:klon,1:klev,it),
!     .            d_tr(1:klon,1:klev,it),ztra_th)
!            DO k=1,klev
!               DO i=1,klon
!                  d_tr(i,k,it)=pdtphys*d_tr(i,k,it)/nsplit
!                  d_tr_th(i,k,it)=d_tr_th(i,k,it)+d_tr(i,k,it)
!                  tr_seri(i,k,it)=MAX(tr_seri(i,k,it)+d_tr(i,k,it),0.)
!               END DO
!            END DO
!         END DO ! nsplit1
!      END DO ! it
!JE end modif 20140323

      DO it=1,nbtr
         DO k=1,klev
            DO i=1,klon
          tmp_var(i,k)=tr_seri(i,k,it)-tmp_var3(i,k,it)
            ENDDO
         ENDDO
       IF (lminmax) THEN
      IF (lcheckmass) THEN
         CALL checkmass(tmp_var(:,:),RNAVO,masse(it),zdz,  &
           pplay,t_seri,iscm3,'dtr therm ')
      ENDIF
       ENDIF
         CALL kg_to_cm3(pplay,t_seri,tmp_var)

         DO k=1,klev
            DO i=1,klon
               his_th(i,it)=his_th(i,it)+    &
                           (tmp_var(i,k))/RNAVO*   &
                     masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys
            END DO !klon
         END DO !klev

      END DO !it
       IF (lminmax) THEN
        DO it=1,nbtr
       CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_after therm')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after therm')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,   &
           pplay,t_seri,iscm3,'after therm')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'after therm')
       ENDIF

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_ther=dife*MAX(0,SIGN(1,dife))   &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_ther=tia_ther+REAL(ti_ther)/REAL(clock_rate)
      ENDIF


      ENDIF ! iflag_conv KE
!------------------------------------
!      Sedimentation
!-----------------------------------
      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF


      DO it=1,nbtr
      DO j=1,klev
      DO i=1,klon
         tmp_var(i,j)=tr_seri(i,j,it)
      ENDDO
      ENDDO
      CALL kg_to_cm3(pplay,t_seri,tmp_var)
      DO j=1,klev
      DO i=1,klon
         tr_seri(i,j,it)=tmp_var(i,j)
      ENDDO
      ENDDO
      ENDDO !--end itr loop
      iscm3=.true.
!--------------------------------------
      print *,' BEFORE Sediment'

      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_before_sedi')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'before sedi')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,   &
           pplay,t_seri,iscm3,'before sedi')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: before sedi')
      ENDIF

      print *,'SPLA VERSION OF SEDIMENTATION IS USED'
      CALL sediment_mod(t_seri,pplay,zrho,paprs,pdtphys,RHcl,   & 
                                     id_coss,id_codu,id_scdu,  &
                                     ok_chimeredust,           &
                         sed_ss,sed_dust,sed_dustsco,          &
                         sed_ss3D,sed_dust3D,sed_dustsco3D,tr_seri)
      CALL cm3_to_kg(pplay,t_seri,sed_ss3D)
      CALL cm3_to_kg(pplay,t_seri,sed_dust3D)
      CALL cm3_to_kg(pplay,t_seri,sed_dustsco3D)

      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_after_sedi')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after sedi')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,  &
           pplay,t_seri,iscm3,'after sedi')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: after sedi')
      ENDIF

!
!=======================================================================
#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRF'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif



!
      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_sedi=dife*MAX(0,SIGN(1,dife))   &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_sedi=tia_sedi+REAL(ti_sedi)/REAL(clock_rate)
      ENDIF

      DO it=1, nbtr
         DO j=1,klev
         DO i=1,klon
           tmp_var(i,j)=tr_seri(i,j,it)
         ENDDO 
         ENDDO
         CALL cm3_to_kg(pplay,t_seri,tmp_var)
         DO j=1,klev
         DO i=1,klon
           tr_seri(i,j,it)=tmp_var(i,j)
         ENDDO 
         ENDDO
      ENDDO
      iscm3=.false.
!
!
!======================================================================
!                      GAS TO PARTICLE CONVERSION      
!======================================================================
!

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF

      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_beforegastopar')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'before gastopar')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,  &
           pplay,t_seri,iscm3,'before gastopar')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: before gastopar')
      ENDIF

      CALL gastoparticle(pdtphys,zdz,zrho,rlat, &
                   pplay,t_seri,id_prec,id_fine, &
                   tr_seri,his_g2pgas ,his_g2paer) 
!
      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_after_gastopar')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after gastopar')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,  &
           pplay,t_seri,iscm3,'after gastopar')
        ENDDO
       ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: after gastopar')
      ENDIF

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_gasp=dife*MAX(0,SIGN(1,dife))   &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_gasp=tia_gasp+REAL(ti_gasp)/REAL(clock_rate)
      ENDIF


!
!======================================================================
!          EFFECT OF PRECIPITATION: iflag_conv=2
!======================================================================
!

#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRG'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif


      IF (iflag_conv.EQ.2) THEN

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF




       DO it=1, nbtr
        DO j=1,klev
        DO i=1,klon
           tmp_var(i,j)=tr_seri(i,j,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
        DO j=1,klev
        DO i=1,klon
           tr_seri(i,j,it)=tmp_var(i,j)
        ENDDO
        ENDDO
      ENDDO
       iscm3=.true.
!------------------------------

      print *,'iflag_conv bef lessiv',iflag_conv
      IF (lessivage) THEN
!
      print *,' BEFORE Incloud'

      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_before_incloud')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'before incloud')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,  &
           pplay,t_seri,iscm3,'before incloud')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: before incloud')
      ENDIF


!      CALL incloud_scav(lminmax,qmin,qmax,masse,henry,kk,prfl,
!     .                  psfl,pmflxr,pmflxs,zrho,zdz,t_seri,pdtphys,

!     .                                     his_dhlsc,his_dhcon,tr_seri)
      print *,'iflag_conv bef incloud',iflag_conv

        IF (iflag_conv.EQ.2) THEN
! Tiedke
      CALL incloud_scav(.false.,qmin,qmax,masse,henry,kk,prfl,          &
                       psfl,pmflxr,pmflxs,zrho,zdz,t_seri,pdtphys,     &
                                          his_dhlsc,his_dhcon,tr_seri)

!---------- to use this option please comment lsc_scav at the end
!        ELSE IF (iflag_conv.GE.3) THEN
!
!      CALL incloud_scav_lsc(.false.,qmin,qmax,masse,henry,kk,prfl,
!     .                  psfl,pmflxr,pmflxs,zrho,zdz,t_seri,pdtphys,
!     .                                     his_dhlsc,his_dhcon,tr_seri)
!--------------------------------------------------------------

        ENDIF
!
!
      print *,' BEFORE blcloud (after incloud)'
      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_before_blcloud')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'before blcloud')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,   &
           pplay,t_seri,iscm3,'before blcloud')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: before blcloud')
      ENDIF

!      CALL blcloud_scav(lminmax,qmin,qmax,pdtphys,prfl,psfl,
!     .                  pmflxr,pmflxs,zdz,alpha_r,alpha_s,masse,
!     .                                  his_dhbclsc,his_dhbccon,tr_seri)

        IF (iflag_conv.EQ.2) THEN
! Tiedke

      CALL blcloud_scav(.false.,qmin,qmax,pdtphys,prfl,psfl,     &
                       pmflxr,pmflxs,zdz,alpha_r,alpha_s,masse,  &
                                       his_dhbclsc,his_dhbccon,tr_seri)

!---------- to use this option please comment lsc_scav at the end
!           and comment IF iflag=2 after "EFFECT OF PRECIPITATION:"
!        
!
!        ELSE IF (iflag_conv.GE.3) THEN
!
!      CALL blcloud_scav_lsc(.false.,qmin,qmax,pdtphys,prfl,psfl,
!     .                  pmflxr,pmflxs,zdz,alpha_r,alpha_s,masse,
!     .                                  his_dhbclsc,his_dhbccon,tr_seri)
!
!----------------------------------------------------------------------
        ENDIF


      print *,' AFTER blcloud '

      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_after_blcloud')
        ENDDO                           
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after blcloud')
        ENDDO                                  
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,  &
           pplay,t_seri,iscm3,'after blcloud')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: after blcloud')
      ENDIF


      ENDIF !--lessivage

      DO it=1, nbtr
         DO j=1,klev
         DO i=1,klon
           tmp_var(i,j)=tr_seri(i,j,it)
         ENDDO 
         ENDDO
         CALL cm3_to_kg(pplay,t_seri,tmp_var)
         DO j=1,klev
         DO i=1,klon
           tr_seri(i,j,it)=tmp_var(i,j)
         ENDDO 
         ENDDO
      ENDDO
       iscm3=.false.
!
      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_wetap=dife*MAX(0,SIGN(1,dife))    &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_wetap=tia_wetap+REAL(ti_wetap)/REAL(clock_rate)
      ENDIF




      ENDIF ! iflag_conv=2

!
!
!======================================================================
!                         EFFECT OF CONVECTION 
!======================================================================
!
#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRH'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif


      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF


      IF (convection) THEN
!
      print *,' BEFORE trconvect'

      IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_before_trconve')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'before trconve')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,  &
           pplay,t_seri,iscm3,'before trconve')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: before trconve')
      ENDIF


! JE        CALL trconvect(pplay,t_seri,pdtphys,pmfu,pmfd,pen_u,pde_u,
!     .             pen_d,pde_d,paprs,zdz,xconv,qmin,qmax,lminmax,masse,
!     .                                                 dtrconv,tr_seri)
! -------------------------------------------------------------     
        IF (iflag_conv.EQ.2) THEN
! Tiedke
         CALL trconvect(pplay,t_seri,pdtphys,pmfu,pmfd,pen_u,pde_u,  &
                  pen_d,pde_d,paprs,zdz,xconv,qmin,qmax,.false.,masse, &
                                                      dtrconv,tr_seri)
         DO it=1, nbtr
           d_tr_cv(:,:,it)=0.
         ENDDO

        ELSE IF (iflag_conv.GE.3) THEN
! KE
         print *,'JE: KE in phytracr_spl'
         DO it=1, nbtr
             DO k = 1, klev
              DO i = 1, klon
               tmp_var3(i,k,it)=tr_seri(i,k,it)
              END DO
             END DO
         ENDDO

         DO it=1, nbtr
!          routine for aerosols . otherwise, check cvltrorig
         print *,'Check sum before cvltr it',it,SUM(tr_seri(:,:,it))
!           IF (.FALSE.) THEN 
           CALL cvltr_spl(pdtphys, da, phi,phi2,d1a,dam, mp,ep,    &
            sigd,sij,wght_cvfd,clw,elij,epmlmMm,eplaMm,           &
            pmflxr,pmflxs,evapls,t_seri,wdtrainA,wdtrainM,          &
!            paprs,it,tr_seri,upwd,dnwd,itop_con,ibas_con,        &
            paprs,it,tmp_var3,upwd,dnwd,itop_con,ibas_con,        &
            henry,kk,zrho,ccntrAA_spla,ccntrENV_spla,coefcoli_spla, &
            id_prec,id_fine,id_coss, id_codu, id_scdu,              &
            d_tr_cv,d_tr_trsp,d_tr_sscav,d_tr_sat,d_tr_uscav,qDi,qPr, &
            qPa,qMel,qTrdi,dtrcvMA,Mint,                            &
            zmfd1a,zmfphi2,zmfdam)
!           ENDIF
!
!           IF (.FALSE.) THEN 
!           CALL cvltr(pdtphys, da, phi,phi2,d1a,dam, mp,ep,
!     .       sigd,sij,wght_cvfd,clw,elij,epmlmMm,eplaMm,
!     .       pmflxr,pmflxs,evapls,t_seri,wdtrainA,wdtrainM,
!     .       paprs,it,tmp_var3,upwd,dnwd,itop_con,ibas_con,
!     .       d_tr_cv,d_tr_trsp,d_tr_sscav,d_tr_sat,d_tr_uscav,qDi,qPr,
!     .       qPa,qMel,qTrdi,dtrcvMA,Mint,
!     .       zmfd1a,zmfphi2,zmfdam)
!!  pas lessivage convectif pou n'est pas un aerosol (i/else with cvltr)
!           ENDIF



!!!!!!!         CALL cvltrorig(it,pdtphys, da, phi,mp,paprs,pplay,tr_seri,
!!!         CALL cvltrorig(it,pdtphys, da, phi,mp,paprs,pplay,tmp_var3,
!!!     .               upwd,dnwd,d_tr_cv)
!             print *,'justbefore cvltrnoscav it= ',it
!             CALL checknanqfi(da(:,:),1.,-1.,' da')
!             CALL checknanqfi(wght_cvfd(:,:),1.,-1.,'weigth ')
!             CALL checknanqfi(mp(:,:),1.,-1.,'mp ')
!             CALL checknanqfi(paprs(:,:),1.,-1.,'paprs ')
!             CALL checknanqfi(pplay(:,:),1.,-1.,'pplay ')
!             CALL checknanqfi(tmp_var3(:,:,it),1.,-1.,'tmp_var3 ')
!             CALL checknanqfi(upwd(:,:),1.,-1.,'upwd ')
!             CALL checknanqfi(dnwd(:,:),1.,-1.,'dnwd ')
!             CALL checknanqfi(d_tr_cv(:,:,it),1.,-1.,'d_tr_cv ')
!             IF (.TRUE.) THEN
!             CALL cvltr_noscav(it,pdtphys, da, phi,mp,wght_cvfd,paprs,
!     .            pplay,tmp_var3,upwd,dnwd,d_tr_cv)
!             ENDIF
             DO k = 1, klev
              DO i = 1, klon
!               tr_seri(i,k,it) = tr_seri(i,k,it) + d_tr_cv(i,k,it)
               tr_seri(i,k,it)=(tmp_var3(i,k,it)+d_tr_cv(i,k,it))
               tmp_var(i,k)=d_tr_cv(i,k,it)

              END DO
             END DO

        CALL kg_to_cm3(pplay,t_seri,tmp_var) !just for his_* computation

             DO k = 1, klev
              DO i = 1, klon
               dtrconv(i,it)=0.0
               his_dhkecv(i,it)=his_dhkecv(i,it)-tmp_var(i,k)  &
                     /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys
              END DO
             END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL kg_to_cm3(pplay,t_seri,tmp_var) !just for his_* computation

             DO k = 1, klev
              DO i = 1, klon
               dtrconv(i,it)=0.0
               his_ds(i,it)=his_ds(i,it)-tmp_var(i,k)  &
                     /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys
              END DO
             END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (lminmax) THEN

         print *,'Check sum after cvltr it',it,SUM(tr_seri(:,:,it))
        CALL minmaxqfi2(d_tr_cv(:,:,it),qmin,qmax,'d_tr_cv:')
        CALL minmaxqfi2(d_tr_trsp(:,:,it),qmin,qmax,'d_tr_trsp:')
        CALL minmaxqfi2(d_tr_sscav(:,:,it),qmin,qmax,'d_tr_sscav:')
        CALL minmaxqfi2(d_tr_sat(:,:,it),qmin,qmax,'d_tr_sat:')
        CALL minmaxqfi2(d_tr_uscav(:,:,it),qmin,qmax,'d_tr_uscav:')
      IF (lcheckmass) THEN
        CALL checkmass(d_tr_cv(:,:,it),RNAVO,masse(it),zdz,  &
           pplay,t_seri,.false.,'d_tr_cv:')
      ENDIF
       ENDIF
         ENDDO ! it=1,nbtr

        ENDIF ! iflag_conv
       IF (lminmax) THEN
        DO it=1,nbtr
        CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_after_trcon')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after trconv')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz, &
           pplay,t_seri,iscm3,'after trconv')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: after trconv')
      ENDIF
      ENDIF ! convection

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_cvltr=dife*MAX(0,SIGN(1,dife))   &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_cvltr=tia_cvltr+REAL(ti_cvltr)/REAL(clock_rate)
      ENDIF



!
!
!=======================================================================
!      LARGE SCALE SCAVENGING KE
!=======================================================================
!     
#ifdef IOPHYS_DUST
      call iophys_ecrit('da',klev,'da','',da)
      call iophys_ecrit('phi',klev,'phi','',phi)
      call iophys_ecrit('phi2',klev,'phi2','',phi2)
      call iophys_ecrit('d1a',klev,'d1a','',d1a)
      call iophys_ecrit('dam',klev,'dam','',dam)
      call iophys_ecrit('mp',klev,'mp','',mp)
      call iophys_ecrit('ep',klev,'ep','',ep)
      call iophys_ecrit('sigd',klev,'sigd','',sigd)
      call iophys_ecrit('sij',klev,'sij','',sij)
      call iophys_ecrit('wght_cvfd',klev,'wght_cvfd','',wght_cvfd)
      call iophys_ecrit('clw',klev,'clw','',clw)
      call iophys_ecrit('elij',klev,'elij','',elij)
      call iophys_ecrit('epmlmMm',klev,'epmlmMm','',epmlmMm)
      call iophys_ecrit('eplaMm',klev,'eplaMm','',eplaMm)
      call iophys_ecrit('pmflxr',klev,'pmflxr','',pmflxr)
      call iophys_ecrit('pmflxs',klev,'pmflxs','',pmflxs)
      call iophys_ecrit('evapls',klev,'evapls','',evapls)
      call iophys_ecrit('wdtrainA',klev,'wdtrainA','',wdtrainA)
      call iophys_ecrit('wdtrainM',klev,'wdtrainM','',wdtrainM)

      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRI'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif


       IF (iflag_conv.GE.3) THEN
       IF (logitime) THEN
       CALL SYSTEM_CLOCK(COUNT=clock_start)
       ENDIF


       IF (lessivage)  THEN
       print *,' BEFORE lsc_scav '
       IF (lminmax) THEN
        DO it=1,nbtr
       CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_before_lsc_scav')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'before lsc_scav')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz,  &
           pplay,t_seri,iscm3,'before lsc_scav')
        ENDDO
      ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: before lsc_scav')
      ENDIF



       ql_incloud_ref = 10.e-4
       ql_incloud_ref =  5.e-4
! calcul du contenu en eau liquide au sein du nuage
       ql_incl = ql_incloud_ref
! choix du lessivage
      IF (iflag_lscav .EQ. 3 .OR. iflag_lscav .EQ. 4) THEN
      print *,'JE iflag_lscav',iflag_lscav
       DO it = 1, nbtr

!       incloud scavenging and removal by large scale rain ! orig : ql_incl
!         was replaced by 0.5e-3 kg/kg
!          the value 0.5e-3 kg/kg is from Giorgi and Chameides (1986), JGR
!         Liu (2001) proposed to use 1.5e-3 kg/kg

!       CALL lsc_scav_orig(pdtphys,it,iflag_lscav,ql_incl,prfl,psfl,
!     .               rneb,beta_fisrt, beta_v1,pplay,paprs,
!     .               t_seri,tr_seri,d_tr_insc, 
!     .               d_tr_bcscav,d_tr_evapls,qPrls)
       CALL lsc_scav_spl(pdtphys,it,iflag_lscav,ql_incl,prfl,psfl,  &
                    rneb,beta_fisrt, beta_v1,pplay,paprs,      &
                    t_seri,tr_seri,d_tr_insc,                  &
                    alpha_r,alpha_s,kk, henry,                 &
                    id_prec,id_fine,id_coss, id_codu, id_scdu, &
                    d_tr_bcscav,d_tr_evapls,qPrls)

!large scale scavenging tendency
       DO k = 1, klev
        DO i = 1, klon
         d_tr_ls(i,k,it)=d_tr_insc(i,k,it)+d_tr_bcscav(i,k,it) &
                        +d_tr_evapls(i,k,it)
         tr_seri(i,k,it)=tr_seri(i,k,it)+d_tr_ls(i,k,it)
          tmp_var(i,k)=d_tr_ls(i,k,it)
        ENDDO
       ENDDO

       CALL kg_to_cm3(pplay,t_seri,tmp_var)
         DO k=1,klev
            DO i=1,klon
            his_dhkelsc(i,it)=his_dhkelsc(i,it)-tmp_var(i,k)    &
                     /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys
     
            END DO
         END DO

       END DO  !tr
      ELSE
        his_dhkelsc(i,it)=0.0
        print *,'WARNING: NO lsc_scav, Please choose iflag_lscav=3 or 4'
       ENDIF !iflag_lscav

       print *,' AFTER lsc_scav '
       IF (lminmax) THEN
        DO it=1,nbtr
       CALL checknanqfi(tr_seri(:,:,it),qmin,qmax,'nan_after_lsc_scav')
        ENDDO
        DO it=1,nbtr
        CALL minmaxqfi2(tr_seri(:,:,it),qmin,qmax,'after lsc_scav')
        ENDDO
      IF (lcheckmass) THEN
        DO it=1,nbtr
         CALL checkmass(tr_seri(:,:,it),RNAVO,masse(it),zdz, &
           pplay,t_seri,iscm3,'after lsc_scav')
        ENDDO
       ENDIF
        CALL minmaxsource(source_tr,qmin,qmax,'src: after lsc_scav')
      ENDIF

      ENDIF ! lessivage
 
      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_lscs=dife*MAX(0,SIGN(1,dife))   &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_lscs=tia_lscs+REAL(ti_lscs)/REAL(clock_rate)
      ENDIF



      ENDIF !iflag_conv

 
!=======================================================================
!                         COMPUTING THE BURDEN
!=======================================================================
#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRJ'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif

!   
      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF

  
      DO it=1, nbtr
        DO j=1,klev
        DO i=1,klon
           tmp_var(i,j)=tr_seri(i,j,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
        DO j=1,klev
        DO i=1,klon
           tr_seri(i,j,it)=tmp_var(i,j)
        ENDDO
        ENDDO
      ENDDO
       iscm3=.true.

!
! Computing burden in mg/m2
      DO it=1, nbtr
      DO k=1, klev
      DO i=1, klon
        trm(i,it)=trm(i,it)+tr_seri(i,k,it)*1.e6*zdz(i,k)*  &
                 masse(it)*1.e3/RNAVO     !--mg S/m2
      ENDDO
      ENDDO
      ENDDO
!
! Computing Surface concentration in ug/m3
!
      DO it=1, nbtr
      DO i=1, klon
        sconc_seri(i,it)=tr_seri(i,1,it)*1.e6* &
                 masse(it)*1.e3/RNAVO     !--mg/m3 (tr_seri ist in g/cm3)
      ENDDO
      ENDDO
!
!=======================================================================
!                  CALCULATION OF OPTICAL PROPERTIES
!=======================================================================
!      
      CALL aeropt_spl(zdz, tr_seri, RHcl,                                 &
                        id_prec, id_fine, id_coss, id_codu, id_scdu,     &
                        ok_chimeredust,                                 &
                    diff_aod550_tot, diag_aod670_tot, diag_aod865_tot,     &
                    diff_aod550_tr2, diag_aod670_tr2, diag_aod865_tr2,     &
                    diag_aod550_ss,  diag_aod670_ss,  diag_aod865_ss,        &
                    diag_aod550_dust,diag_aod670_dust,diag_aod865_dust,  &
           diag_aod550_dustsco,diag_aod670_dustsco,diag_aod865_dustsco)  



      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)
      dife=clock_end-clock_start
      ti_brop=dife*MAX(0,SIGN(1,dife))   &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_brop=tia_brop+REAL(ti_brop)/REAL(clock_rate)
      ENDIF


!=======================================================================
!   MODIS terra/aqua simulation output
!=======================================================================
      masque_aqua_cur(:)=0
      masque_terra_cur(:)=0

      CALL satellite_out_spla(jD_cur,jH_cur,pdtphys,rlat,rlon,   &
                              masque_aqua_cur, masque_terra_cur )
      IF (jH_cur-pdtphys/86400. .LT. 0.) THEN
       !new utc day: put in 0 everything
!JE20150518<<
!       aod550_aqua(:) =0.
!       aod550_terra(:) =0.
!       aod670_aqua(:) =0.
!       aod670_terra(:) =0.
!       aod865_aqua(:) =0.
!       aod865_terra(:) =0.
       masque_aqua(:) =0
       masque_terra(:) =0
       aod550_terra(:)=0.  
       aod550_tr2_terra(:)=0.  
       aod550_ss_terra(:)=0.   
       aod550_dust_terra(:)=0.   
       aod550_dustsco_terra(:)=0.   
       aod670_terra(:)=0.   
       aod670_tr2_terra(:)=0.  
       aod670_ss_terra(:)=0.  
       aod670_dust_terra(:)=0.  
       aod670_dustsco_terra(:)=0.  
       aod865_terra(:)=0.   
       aod865_tr2_terra(:)=0.  
       aod865_ss_terra(:)=0.  
       aod865_dust_terra(:)=0.  
       aod865_dustsco_terra(:)=0.  
       aod550_aqua(:)=0.  
       aod550_tr2_aqua(:)=0.  
       aod550_ss_aqua(:)=0.   
       aod550_dust_aqua(:)=0.   
       aod550_dustsco_aqua(:)=0.   
       aod670_aqua(:)=0.   
       aod670_tr2_aqua(:)=0.  
       aod670_ss_aqua(:)=0.  
       aod670_dust_aqua(:)=0.  
       aod670_dustsco_aqua(:)=0.  
       aod865_aqua(:)=0.   
       aod865_tr2_aqua(:)=0.  
       aod865_ss_aqua(:)=0.  
       aod865_dust_aqua(:)=0.  
       aod865_dustsco_aqua(:)=0.  
!JE20150518>>
      ENDIF

      DO i=1,klon
!         aod550_aqua(i)=aod550_aqua(i)+   &
!                       masque_aqua_cur(i)*diff_aod550_tot(i)
!         aod670_aqua(i)=aod670_aqua(i)+   &
!                        masque_aqua_cur(i)*diag_aod670_tot(i)
!         aod865_aqua(i)=aod865_aqua(i)+   &
!                       masque_aqua_cur(i)*diag_aod865_tot(i)

       aod550_terra(i)=aod550_terra(i)+   &
                       masque_terra_cur(i)*diff_aod550_tot(i)
       aod550_tr2_terra(i)= aod550_tr2_terra(i)+ &
                       masque_terra_cur(i)*diff_aod550_tr2(i)
       aod550_ss_terra(i)=aod550_ss_terra(i) + &
                       masque_terra_cur(i)*diag_aod550_ss(i)
       aod550_dust_terra(i)=  aod550_dust_terra(i) + &
                       masque_terra_cur(i)*diag_aod550_dust(i)
       aod550_dustsco_terra(i)= aod550_dustsco_terra(i) + &
                       masque_terra_cur(i)*diag_aod550_dustsco(i)
       aod670_terra(i)=aod670_terra(i)+   &
                       masque_terra_cur(i)*diag_aod670_tot(i)
       aod670_tr2_terra(i)= aod670_tr2_terra(i)+ &
                       masque_terra_cur(i)*diag_aod670_tr2(i)
       aod670_ss_terra(i)=aod670_ss_terra(i) + &
                       masque_terra_cur(i)*diag_aod670_ss(i)
       aod670_dust_terra(i)=  aod670_dust_terra(i) + &
                       masque_terra_cur(i)*diag_aod670_dust(i)
       aod670_dustsco_terra(i)= aod670_dustsco_terra(i) + &
                       masque_terra_cur(i)*diag_aod670_dustsco(i)
       aod865_terra(i)=aod865_terra(i)+   &
                       masque_terra_cur(i)*diag_aod865_tot(i)
       aod865_tr2_terra(i)= aod865_tr2_terra(i)+ &
                       masque_terra_cur(i)*diag_aod865_tr2(i)
       aod865_ss_terra(i)=aod865_ss_terra(i) + &
                       masque_terra_cur(i)*diag_aod865_ss(i)
       aod865_dust_terra(i)=  aod865_dust_terra(i) + &
                       masque_terra_cur(i)*diag_aod865_dust(i)
       aod865_dustsco_terra(i)= aod865_dustsco_terra(i) + &
                       masque_terra_cur(i)*diag_aod865_dustsco(i)



       aod550_aqua(i)=aod550_aqua(i)+   &
                       masque_aqua_cur(i)*diff_aod550_tot(i)
       aod550_tr2_aqua(i)= aod550_tr2_aqua(i)+ &
                       masque_aqua_cur(i)*diff_aod550_tr2(i)
       aod550_ss_aqua(i)=aod550_ss_aqua(i) + &
                       masque_aqua_cur(i)*diag_aod550_ss(i)
       aod550_dust_aqua(i)=  aod550_dust_aqua(i) + &
                       masque_aqua_cur(i)*diag_aod550_dust(i)
       aod550_dustsco_aqua(i)= aod550_dustsco_aqua(i) + &
                       masque_aqua_cur(i)*diag_aod550_dustsco(i)
       aod670_aqua(i)=aod670_aqua(i)+   &
                       masque_aqua_cur(i)*diag_aod670_tot(i)
       aod670_tr2_aqua(i)= aod670_tr2_aqua(i)+ &
                       masque_aqua_cur(i)*diag_aod670_tr2(i)
       aod670_ss_aqua(i)=aod670_ss_aqua(i) + &
                       masque_aqua_cur(i)*diag_aod670_ss(i)
       aod670_dust_aqua(i)=  aod670_dust_aqua(i) + &
                       masque_aqua_cur(i)*diag_aod670_dust(i)
       aod670_dustsco_aqua(i)= aod670_dustsco_aqua(i) + &
                       masque_aqua_cur(i)*diag_aod670_dustsco(i)
       aod865_aqua(i)=aod865_aqua(i)+   &
                       masque_aqua_cur(i)*diag_aod865_tot(i)
       aod865_tr2_aqua(i)= aod865_tr2_aqua(i)+ &
                       masque_aqua_cur(i)*diag_aod865_tr2(i)
       aod865_ss_aqua(i)=aod865_ss_aqua(i) + &
                       masque_aqua_cur(i)*diag_aod865_ss(i)
       aod865_dust_aqua(i)=  aod865_dust_aqua(i) + &
                       masque_aqua_cur(i)*diag_aod865_dust(i)
       aod865_dustsco_aqua(i)= aod865_dustsco_aqua(i) + &
                       masque_aqua_cur(i)*diag_aod865_dustsco(i)
!         aod550_terra(i)=aod550_terra(i)+  &
!                       masque_terra_cur(i)*diff_aod550_tot(i)
!         aod670_terra(i)=aod670_terra(i)+  &
!                       masque_terra_cur(i)*diag_aod670_tot(i)
!         aod865_terra(i)=aod865_terra(i)+   &
!                       masque_terra_cur(i)*diag_aod865_tot(i)
         masque_aqua(i)=masque_aqua(i)+masque_aqua_cur(i)
         masque_terra(i)=masque_terra(i)+masque_terra_cur(i)
      ENDDO

      IF (jH_cur+pdtphys/86400. .GE. 1.) THEN  
!          print *,'last step of the day'
          DO i=1,klon
               IF (masque_aqua(i).GT. 0) THEN
                   aod550_aqua(i)=aod550_aqua(i)/masque_aqua(i)
                   aod670_aqua(i)=aod670_aqua(i)/masque_aqua(i)
                   aod865_aqua(i)=aod865_aqua(i)/masque_aqua(i)
                   aod550_tr2_aqua(i)=aod550_tr2_aqua(i)/masque_aqua(i)
                   aod670_tr2_aqua(i)=aod670_tr2_aqua(i)/masque_aqua(i)
                   aod865_tr2_aqua(i)=aod865_tr2_aqua(i)/masque_aqua(i)
                   aod550_ss_aqua(i)=aod550_ss_aqua(i)/masque_aqua(i)
                   aod670_ss_aqua(i)=aod670_ss_aqua(i)/masque_aqua(i)
                   aod865_ss_aqua(i)=aod865_ss_aqua(i)/masque_aqua(i)
                   aod550_dust_aqua(i)=aod550_dust_aqua(i)/masque_aqua(i)
                   aod670_dust_aqua(i)=aod670_dust_aqua(i)/masque_aqua(i)
                   aod865_dust_aqua(i)=aod865_dust_aqua(i)/masque_aqua(i)
                   aod550_dustsco_aqua(i)=aod550_dustsco_aqua(i)/masque_aqua(i)
                   aod670_dustsco_aqua(i)=aod670_dustsco_aqua(i)/masque_aqua(i)
                   aod865_dustsco_aqua(i)=aod865_dustsco_aqua(i)/masque_aqua(i)
               ELSE 
                   aod550_aqua(i) = -999.
                   aod670_aqua(i) = -999.
                   aod865_aqua(i) = -999.
                   aod550_tr2_aqua(i)= -999.
                   aod670_tr2_aqua(i)= -999.
                   aod865_tr2_aqua(i)= -999.
                   aod550_ss_aqua(i)= -999.
                   aod670_ss_aqua(i)= -999.
                   aod865_ss_aqua(i)= -999.
                   aod550_dust_aqua(i)= -999.
                   aod670_dust_aqua(i)= -999.
                   aod865_dust_aqua(i)= -999.
                   aod550_dustsco_aqua(i)= -999.
                   aod670_dustsco_aqua(i)= -999.
                   aod865_dustsco_aqua(i)= -999.
               ENDIF
               IF (masque_terra(i).GT. 0) THEN
                   aod550_terra(i)=aod550_terra(i)/masque_terra(i)
                   aod670_terra(i)=aod670_terra(i)/masque_terra(i)
                   aod865_terra(i)=aod865_terra(i)/masque_terra(i)
                   aod550_tr2_terra(i)=aod550_tr2_terra(i)/masque_terra(i)
                   aod670_tr2_terra(i)=aod670_tr2_terra(i)/masque_terra(i)
                   aod865_tr2_terra(i)=aod865_tr2_terra(i)/masque_terra(i)
                   aod550_ss_terra(i)=aod550_ss_terra(i)/masque_terra(i)
                   aod670_ss_terra(i)=aod670_ss_terra(i)/masque_terra(i)
                   aod865_ss_terra(i)=aod865_ss_terra(i)/masque_terra(i)
                   aod550_dust_terra(i)=aod550_dust_terra(i)/masque_terra(i)
                   aod670_dust_terra(i)=aod670_dust_terra(i)/masque_terra(i)
                   aod865_dust_terra(i)=aod865_dust_terra(i)/masque_terra(i)
                   aod550_dustsco_terra(i)=aod550_dustsco_terra(i)/masque_terra(i)
                   aod670_dustsco_terra(i)=aod670_dustsco_terra(i)/masque_terra(i)
                   aod865_dustsco_terra(i)=aod865_dustsco_terra(i)/masque_terra(i)
               ELSE 
                   aod550_terra(i) = -999.
                   aod670_terra(i) = -999.
                   aod865_terra(i) = -999.
                   aod550_tr2_terra(i)= -999.
                   aod670_tr2_terra(i)= -999.
                   aod865_tr2_terra(i)= -999.
                   aod550_ss_terra(i)= -999.
                   aod670_ss_terra(i)= -999.
                   aod865_ss_terra(i)= -999.
                   aod550_dust_terra(i)= -999.
                   aod670_dust_terra(i)= -999.
                   aod865_dust_terra(i)= -999.
                   aod550_dustsco_terra(i)= -999.
                   aod670_dustsco_terra(i)= -999.
                   aod865_dustsco_terra(i)= -999.
               ENDIF
!              IF (masque_terra(i).GT. 0) THEN
!                   aod550_terra(i) = aod550_terra(i)/masque_terra(i)
!                   aod670_terra(i)=aod670_terra(i)/masque_terra(i)
!                   aod865_terra(i)=aod865_terra(i)/masque_terra(i)
!
!               ELSE
!                   aod550_terra(i) = -999.
!                   aod670_terra(i) = -999.
!                   aod865_terra(i) = -999.
!               ENDIF
          ENDDO          
!      !write  dbg
!       CALL writefield_phy("aod550_aqua",aod550_aqua,1) 
!       CALL writefield_phy("aod550_terra",aod550_terra,1) 
!       CALL writefield_phy("masque_aqua",float(masque_aqua),1)
!       CALL writefield_phy("masque_terra",float(masque_terra),1)


      IF (ok_histrac) THEN
!      write in output file
      call gather(aod550_aqua,aod550_aqua_glo)
      call gather(aod550_terra,aod550_terra_glo)
      call gather(aod670_aqua,aod670_aqua_glo)
      call gather(aod670_terra,aod670_terra_glo)
      call gather(aod865_aqua,aod865_aqua_glo)
      call gather(aod865_terra,aod865_terra_glo)

!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN

      CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, aod550_aqua_glo ,zx_tmp_2d)
      CALL histwrite(nid_tra3,"taue550_aqua",itra,zx_tmp_2d, &
                                      nbp_lon*(nbp_lat),ndex2d)

      CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, aod550_terra_glo ,zx_tmp_2d)
      CALL histwrite(nid_tra3,"taue550_terra",itra,zx_tmp_2d, &
                                      nbp_lon*(nbp_lat),ndex2d)
      CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, aod670_aqua_glo ,zx_tmp_2d)
      CALL histwrite(nid_tra3,"taue670_aqua",itra,zx_tmp_2d, &
                                      nbp_lon*(nbp_lat),ndex2d)
      CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, aod670_terra_glo ,zx_tmp_2d)
      CALL histwrite(nid_tra3,"taue670_terra",itra,zx_tmp_2d, &
                                      nbp_lon*(nbp_lat),ndex2d)

      CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, aod865_aqua_glo ,zx_tmp_2d)
      CALL histwrite(nid_tra3,"taue865_aqua",itra,zx_tmp_2d, &
                                      nbp_lon*(nbp_lat),ndex2d)
      CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, aod865_terra_glo ,zx_tmp_2d)
      CALL histwrite(nid_tra3,"taue865_terra",itra,zx_tmp_2d, &
                                      nbp_lon*(nbp_lat),ndex2d)
      ENDIF
!$OMP END MASTER
!$OMP BARRIER
      ENDIF
!       !put in 0 everything
!       aod550_aqua(:) =0.
!       aod550_terra(:) =0.
!       aod670_aqua(:) =0.
!       aod670_terra(:) =0.
!       aod865_aqua(:) =0.
!       aod865_terra(:) =0.
!       masque_aqua(:) =0 
!       masque_terra(:) =0 
      ENDIF 


!
!======================================================================
!  Stockage sur bande histoire
!======================================================================
#ifdef IOPHYS_DUST
      do it=1,nbtr
         write(str2,'(i2.2)') it
         call iophys_ecrit('TRK'//str2,klev,'SOURCE','',tr_seri(:,:,it))
      enddo
#endif


!
      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ENDIF

      DO it=1, nbtr
         DO j=1,klev
         DO i=1,klon
           tmp_var(i,j)=tr_seri(i,j,it)
         ENDDO 
         ENDDO
         CALL cm3_to_kg(pplay,t_seri,tmp_var)
         DO j=1,klev
         DO i=1,klon
           tr_seri(i,j,it)=tmp_var(i,j)
         ENDDO 
         ENDDO
      ENDDO
       iscm3=.false.

!
!
!======================================================================
!  SAVING AEROSOL RELATED VARIABLES INTO FILE
!======================================================================
!
!JE20141224      IF (ok_histrac) THEN
!
      ndex2d = 0
      ndex3d = 0 
!
      itra=itra+1

      print *,'SAVING VARIABLES FOR DAY ',itra
!
      fluxbb(:)=0.0
      fluxff(:)=0.0
      fluxbcbb(:)=0.0
      fluxbcff(:)=0.0
      fluxbcnff(:)=0.0
      fluxbcba(:)=0.0
      fluxbc(:)=0.0
      fluxombb(:)=0.0
      fluxomff(:)=0.0
      fluxomnat(:)=0.0
      fluxomba(:)=0.0
      fluxomnff(:)=0.0
      fluxom(:)=0.0
      fluxh2sff(:)=0.0
      fluxh2snff(:)=0.0
      fluxh2sbio(:)=0.0
      fluxso2ff(:)=0.0
      fluxso2nff(:)=0.0
      fluxso2bb(:)=0.0
      fluxso2vol(:)=0.0
      fluxso2ba(:)=0.0
      fluxso2(:)=0.0
      fluxso4ff(:)=0.0
      fluxso4nff(:)=0.0
      fluxso4bb(:)=0.0
      fluxso4ba(:)=0.0
      fluxso4(:)=0.0
      fluxdms(:)=0.0
      fluxdustec(:)=0.0
      fluxddfine(:)=0.0
      fluxddcoa(:)=0.0
      fluxddsco(:)=0.0
      fluxdd(:)=0.0
      fluxssfine(:)=0.0
      fluxsscoa(:)=0.0
      fluxss(:)=0.0
      DO i=1, klon
         IF (iregion_ind(i).GT.0) THEN           ! LAND
           ! SULFUR EMISSIONS
           fluxh2sff(i)= (lmt_so2ff_l(i)+lmt_so2ff_h(i))*frach2sofso2*  &       
                         scale_param_ind(iregion_ind(i))*               &
                                    1.e4/RNAVO*masse_s*1.e3         ! mgS/m2/s
           fluxso2ff(i)=scale_param_ind(iregion_ind(i)) * fracso2emis * &
                        (lmt_so2ff_l(i)+lmt_so2ff_h(i)) * 1.e4/RNAVO * &
                                                    masse_s * 1.e3  ! mgS/m2/s
           ! SULPHATE EMISSIONS
           fluxso4ff(i)=scale_param_ind(iregion_ind(i))*(1-fracso2emis)* &
                         (lmt_so2ff_l(i)+lmt_so2ff_h(i)) * 1.e4/RNAVO * &
                                                    masse_s * 1.e3  ! mgS/m2/s
           ! BLACK CARBON EMISSIONS
           fluxbcff(i)=scale_param_ff(iregion_ind(i))* &
                                             lmt_bcff(i)*1.e4*1.e3  !/g/m2/s
           ! ORGANIC MATTER EMISSIONS
           fluxomff(i)=scale_param_ff(iregion_ind(i))* &
                               (lmt_omff(i))*1.e4*1.e3  !/g/m2/s
           ! FOSSIL FUEL EMISSIONS
           fluxff(i)=fluxbcff(i)+fluxomff(i)
         ENDIF
         IF (iregion_bb(i).GT.0) THEN           ! LAND
           ! SULFUR EMISSIONS
           fluxso2bb(i) =scale_param_bb(iregion_bb(i)) * fracso2emis *  &
                      (lmt_so2bb_l(i)+lmt_so2bb_h(i))*                 &
                (1.-pctsrf(i,is_oce))*1.e4/RNAVO*masse_s*1.e3       ! mgS/m2/s
           ! SULPHATE EMISSIONS
           fluxso4bb(i) =scale_param_bb(iregion_bb(i))*(1-fracso2emis)* &
                      (lmt_so2bb_l(i)+lmt_so2bb_h(i))*                 &
                (1.-pctsrf(i,is_oce))*1.e4/RNAVO*masse_s*1.e3       ! mgS/m2/s
           ! BLACK CARBON EMISSIONS
           fluxbcbb(i)=scale_param_bb(iregion_bb(i))*                   &
                           (lmt_bcbb_l(i)+lmt_bcbb_h(i))*1.e4*1.e3  !mg/m2/s
           ! ORGANIC MATTER EMISSIONS
           fluxombb(i)=scale_param_bb(iregion_bb(i))*                   &
                           (lmt_ombb_l(i)+lmt_ombb_h(i))*1.e4*1.e3  !mg/m2/s
           ! BIOMASS BURNING EMISSIONS
           fluxbb(i)=fluxbcbb(i)+fluxombb(i)
         ENDIF
         ! H2S EMISSIONS
         fluxh2sbio(i)=lmt_h2sbio(i)*1.e4/RNAVO*masse_s*1.e3      ! mgS/m2/s
         fluxh2snff(i)= lmt_so2nff(i)*frach2sofso2*  &
                                    1.e4/RNAVO*masse_s*1.e3         ! mgS/m2/s
         ! SULFUR DIOXIDE EMISSIONS
         fluxso2nff(i)=fracso2emis * lmt_so2nff(i) * 1.e4/RNAVO *  &
                                                    masse_s * 1.e3  ! mgS/m2/s
         fluxso2vol(i)=(lmt_so2volc_cont(i)+lmt_so2volc_expl(i))  &
                      *1.e4/RNAVO*masse_s*1.e3        ! mgS/m2/s
         fluxso2ba(i) =lmt_so2ba(i)*1.e4/RNAVO*masse_s*1.e3*      &
                                                        fracso2emis ! mgS/m2/s
         fluxso2(i)=fluxso2ff(i)+fluxso2bb(i)+fluxso2nff(i)+   &
                   fluxso2vol(i)+fluxso2ba(i)
         ! DMS EMISSIONS
         fluxdms(i)=( lmt_dms(i)+lmt_dmsbio(i) )              &
                   *1.e4/RNAVO*masse_s*1.e3          ! mgS/m2/s
         ! SULPHATE EMISSIONS
         fluxso4ba(i) =lmt_so2ba(i)*1.e4/RNAVO*masse_s*1.e3        &
                      *(1-fracso2emis) ! mgS/m2/s
         fluxso4nff(i)=(1-fracso2emis)*lmt_so2nff(i) * 1.e4/RNAVO *  &
                                                    masse_s * 1.e3  ! mgS/m2/s
         fluxso4(i)=fluxso4ff(i)+fluxso4bb(i)+fluxso4ba(i)+fluxso4nff(i)
         ! BLACK CARBON EMISSIONS

         fluxbcnff(i)=lmt_bcnff(i)*1.e4*1.e3  !mg/m2/s
         fluxbcba(i)=lmt_bcba(i)*1.e4*1.e3    !mg/m2/s
         fluxbc(i)=fluxbcbb(i)+fluxbcff(i)+fluxbcnff(i)+fluxbcba(i)
         ! ORGANIC MATTER EMISSIONS
         fluxomnat(i)=lmt_omnat(i)*1.e4*1.e3  !mg/m2/s
         fluxomba(i)=lmt_omba(i)*1.e4*1.e3  !mg/m2/s
         fluxomnff(i)=lmt_omnff(i)*1.e4*1.e3  !mg/m2/s
         fluxom(i)=fluxombb(i)+fluxomff(i)+fluxomnat(i)+fluxomba(i)+  &
                  fluxomnff(i)
        ! DUST EMISSIONS 
         fluxdustec(i)=dust_ec(i)*1.e6 ! old dust emission scheme 
!JE20140605<<         old dust emission version
!         fluxddfine(i)=scale_param_dustacc(iregion_dust(i))
!     .                                  * dust_ec(i)*0.093*1.e6
!         fluxddcoa(i)=scale_param_dustcoa(iregion_dust(i))
!     .                                  * dust_ec(i)*0.905*1.e6
!         fluxdd(i)=fluxddfine(i)+fluxddcoa(i)
!JE20140605>>
         fluxddfine(i)=flux_sparam_ddfine(i)
         fluxddcoa(i)=flux_sparam_ddcoa(i)
         fluxddsco(i)=flux_sparam_ddsco(i)
         fluxdd(i)=fluxddfine(i)+fluxddcoa(i)+fluxddsco(i)
        ! SEA SALT EMISSIONS 
         fluxssfine(i)=scale_param_ssacc*lmt_sea_salt(i,1)*1.e4*1.e3
         fluxsscoa(i)=scale_param_sscoa*lmt_sea_salt(i,2)*1.e4*1.e3
         fluxss(i)=fluxssfine(i)+fluxsscoa(i)
      ENDDO
!      prepare outputs cvltr

      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_cv(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_cv_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO
      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_trsp(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_trsp_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO
      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_sscav(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_sscav_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO
      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_sat(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_sat_o(i,k,it)=tmp_var(i,k)   &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO
      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_uscav(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_uscav_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_insc(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_insc_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO
     

      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_bcscav(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_bcscav_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO


      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_evapls(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_evapls_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO


      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_ls(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_ls_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO


      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_dyn(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_dyn_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO


      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_cl(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_cl_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO


      DO it=1, nbtr
        DO k=1,klev
        DO i=1,klon
           tmp_var(i,k)=d_tr_th(i,k,it)
        ENDDO
        ENDDO
        CALL kg_to_cm3(pplay,t_seri,tmp_var)
       DO k=1,klev
        DO i=1,klon
          d_tr_th_o(i,k,it)=tmp_var(i,k)  &
                         /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  
        ENDDO
       ENDDO
      ENDDO
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     DO it=1,nbtr
      WRITE(str2,'(i2.2)') it
       DO i=1, klon                                                        
        his_dh(i,it)= his_dhlsc(i,it)+his_dhcon(i,it)+               &
                   his_dhbclsc(i,it)+his_dhbccon(i,it)

       ENDDO
      ENDDO

      IF (ok_histrac) THEN
!
! SAVING VARIABLES IN TRACEUR
!
     call gather(diff_aod550_tot  ,auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo ,zx_tmp_2d_glo)
     CALL histwrite(nid_tra3,"taue550",itra,zx_tmp_2d_glo,                 &
                                      nbp_lon*(nbp_lat),ndex2d)              
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod670_tot  , auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo,zx_tmp_2d_glo)     
     CALL histwrite(nid_tra3,"taue670",itra,zx_tmp_2d_glo,                 &   
                                      nbp_lon*(nbp_lat),ndex2d)              
!                                                                       
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod865_tot  , auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo,zx_tmp_2d_glo)     
     CALL histwrite(nid_tra3,"taue865",itra,zx_tmp_2d_glo,                 &  
                                      nbp_lon*(nbp_lat),ndex2d)              
!                                                                       
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather(  diff_aod550_tr2 , auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo,zx_tmp_2d_glo)     
     CALL histwrite(nid_tra3,"taue550_tr2",itra,zx_tmp_2d_glo,             &  
                                      nbp_lon*(nbp_lat),ndex2d)              
!                                                                       
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather(  diag_aod670_tr2 , auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo,zx_tmp_2d_glo)     
     CALL histwrite(nid_tra3,"taue670_tr2",itra,zx_tmp_2d_glo,             &  
                                      nbp_lon*(nbp_lat),ndex2d)              
!                                                                       
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod865_tr2  , auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo,zx_tmp_2d_glo)     
     CALL histwrite(nid_tra3,"taue865_tr2",itra,zx_tmp_2d_glo,             &  
                                      nbp_lon*(nbp_lat),ndex2d)              
!                                                                       
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather(  diag_aod550_ss, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo ,zx_tmp_2d_glo)      
     CALL histwrite(nid_tra3,"taue550_ss",itra,zx_tmp_2d_glo,              &  
                                      nbp_lon*(nbp_lat),ndex2d)              
!                                                                       
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod670_ss , auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo ,zx_tmp_2d_glo)      
     CALL histwrite(nid_tra3,"taue670_ss",itra,zx_tmp_2d_glo,              &  
                                      nbp_lon*(nbp_lat),ndex2d)              
!                                                                       
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod865_ss, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo ,zx_tmp_2d_glo)      
     CALL histwrite(nid_tra3,"taue865_ss",itra,zx_tmp_2d_glo,              &  
                                      nbp_lon*(nbp_lat),ndex2d)              
!                                                                       
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod550_dust, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo ,zx_tmp_2d_glo)    
     CALL histwrite(nid_tra3,"taue550_dust",itra,zx_tmp_2d_glo,             & 
                                      nbp_lon*(nbp_lat),ndex2d)               
!                                                                        
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod670_dust, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo ,zx_tmp_2d_glo)     
     CALL histwrite(nid_tra3,"taue670_dust",itra,zx_tmp_2d_glo,             &  
                                      nbp_lon*(nbp_lat),ndex2d)               
!                                                                        
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod865_dust, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo ,zx_tmp_2d_glo)     
     CALL histwrite(nid_tra3,"taue865_dust",itra,zx_tmp_2d_glo,             &  
                                      nbp_lon*(nbp_lat),ndex2d)               
!                                                                        
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod550_dustsco, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo ,zx_tmp_2d_glo)  
     CALL histwrite(nid_tra3,"taue550_dustsco",itra,zx_tmp_2d_glo,          &  
                                      nbp_lon*(nbp_lat),ndex2d)               
!                                                                        
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod670_dustsco, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo ,zx_tmp_2d_glo)  
     CALL histwrite(nid_tra3,"taue670_dustsco",itra,zx_tmp_2d_glo,          &  
                                      nbp_lon*(nbp_lat),ndex2d)               
!                                                                        
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( diag_aod865_dustsco, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1, klon_glo,nbp_lon,nbp_lat, auxklon_glo ,zx_tmp_2d_glo)  
     CALL histwrite(nid_tra3,"taue865_dustsco",itra,zx_tmp_2d_glo,          &  
                                      nbp_lon*(nbp_lat),ndex2d)               
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
                                                                         
!$OMP MASTER
     DO it=1,nbtr                                                        
!                                                                        
     WRITE(str2,'(i2.2)') it
! 
     call gather( trm, auxklonnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklonnbtr_glo(1,it) , zx_tmp_2d_glo)
     CALL histwrite(nid_tra3,"trm"//str2,itra,zx_tmp_2d_glo,              & 
                                         nbp_lon*(nbp_lat),ndex2d)          
!                                                                      
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
     call gather( sconc_seri, auxklonnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklonnbtr_glo(1,it) ,zx_tmp_2d_glo)     
     CALL histwrite(nid_tra3,"sconc"//str2,itra,zx_tmp_2d_glo,            &  
                                         nbp_lon*(nbp_lat),ndex2d)          
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
!                                                                      
! SAVING VARIABLES IN LESSIVAGE                                         
!                                                                       
     call gather( flux_tr, auxklonnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklonnbtr_glo(1,it) ,zx_tmp_2d_glo)
     CALL histwrite(nid_tra2,"flux"//str2,itra,zx_tmp_2d_glo,               & 
                    nbp_lon*(nbp_lat),ndex2d)                                 
!                                                                        
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
     call gather( his_ds, auxklonnbtr_glo )
!! $OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklonnbtr_glo(1,it) ,zx_tmp_2d_glo)           
     CALL histwrite(nid_tra2,"ds"//str2,itra,zx_tmp_2d_glo,                 &  
                    nbp_lon*(nbp_lat),ndex2d)                                 
!                                                                        
      ENDIF
! !$OMP END MASTER
! !$OMP BARRIER
      ENDDO

     DO it=1,nbtr
     WRITE(str2,'(i2.2)') it
      DO i=1, klon                                                        
       zx_tmp_fi2d(i) = his_dhlsc(i,it)+his_dhcon(i,it)+               &  
                        his_dhbclsc(i,it)+his_dhbccon(i,it)
       his_dh(i,it)= his_dhlsc(i,it)+his_dhcon(i,it)+               &  
                   his_dhbclsc(i,it)+his_dhbccon(i,it)

      ENDDO
!
     call gather( zx_tmp_fi2d, auxklon_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)
     CALL histwrite(nid_tra2,"dh"//str2,itra,zx_tmp_2d_glo,                  &
                    nbp_lon*(nbp_lat),ndex2d)                                  
!                                                                         
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
     call gather( his_dhkecv, auxklonnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklonnbtr_glo(1,it) ,zx_tmp_2d_glo)        
     CALL histwrite(nid_tra2,"dhkecv"//str2,itra,zx_tmp_2d_glo,              &  
                    nbp_lon*(nbp_lat),ndex2d)       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           
!                                                                         
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
     call gather( his_dhkelsc, auxklonnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklonnbtr_glo(1,it) ,zx_tmp_2d_glo)       
     CALL histwrite(nid_tra2,"dhkelsc"//str2,itra,zx_tmp_2d_glo,             &  
                    nbp_lon*(nbp_lat),ndex2d)                                  
!                                                                         
                                                                          
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
!    call gather( d_tr_cv_o,  auxklonklevnbtr_glo )
     call gather( d_tr_cv,  auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,             &  
                      zx_tmp_3d_glo)                                         
     CALL histwrite(nid_tra2,"d_tr_cv"//str2,itra,zx_tmp_3d_glo,             &  
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
!
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
     call gather( d_tr_trsp_o, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,           &    
                      zx_tmp_3d_glo)                                            
     CALL histwrite(nid_tra2,"d_tr_trsp"//str2,itra,zx_tmp_3d_glo,           &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
!
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
     call gather( d_tr_sscav_o, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,          &    
                      zx_tmp_3d_glo)                                            
     CALL histwrite(nid_tra2,"d_tr_sscav"//str2,itra,zx_tmp_3d_glo,          &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
!
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
     call gather( d_tr_sat_o, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,            &    
                      zx_tmp_3d_glo)                                            
     CALL histwrite(nid_tra2,"d_tr_sat"//str2,itra,zx_tmp_3d_glo,            &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
!
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
     call gather( d_tr_uscav_o, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,          &    
                       zx_tmp_3d_glo)                                           
     CALL histwrite(nid_tra2,"d_tr_uscav"//str2,itra,zx_tmp_3d_glo,          &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
!                                                                            
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call gather( d_tr_insc_o, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,          &    
                       zx_tmp_3d_glo)                                           
     CALL histwrite(nid_tra2,"d_tr_insc"//str2,itra,zx_tmp_3d_glo,          &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
!                                                                            
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call gather( d_tr_bcscav_o, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,          &    
                       zx_tmp_3d_glo)                                           
     CALL histwrite(nid_tra2,"d_tr_bcscav"//str2,itra,zx_tmp_3d_glo,          &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
!                                                                            
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call gather( d_tr_evapls_o, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,          &    
                       zx_tmp_3d_glo)                                           
     CALL histwrite(nid_tra2,"d_tr_evapls"//str2,itra,zx_tmp_3d_glo,          &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
!                                                                            
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    call gather( d_tr_ls_o, auxklonklevnbtr_glo )
     call gather( d_tr_ls, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,          &    
                       zx_tmp_3d_glo)                                           
     CALL histwrite(nid_tra2,"d_tr_ls"//str2,itra,zx_tmp_3d_glo,          &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
                                                                            
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    call gather( d_tr_dyn_o, auxklonklevnbtr_glo )
     call gather( d_tr_dyn, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,          &    
                       zx_tmp_3d_glo)                                           
     CALL histwrite(nid_tra2,"d_tr_dyn"//str2,itra,zx_tmp_3d_glo,          &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
                                                                            
      print*,'ECRTIURES TENDANCES MODIFIEES NON MAIS'
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    call gather( d_tr_cl_o, auxklonklevnbtr_glo )
     call gather( d_tr_cl, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,          &    
                       zx_tmp_3d_glo)                                           
     CALL histwrite(nid_tra2,"d_tr_cl"//str2,itra,zx_tmp_3d_glo,          &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
                                                                            
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    call gather( d_tr_th_o, auxklonklevnbtr_glo )
     call gather( d_tr_th, auxklonklevnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,auxklonklevnbtr_glo(1,1,it) ,          &    
                       zx_tmp_3d_glo)                                           
     CALL histwrite(nid_tra2,"d_tr_th"//str2,itra,zx_tmp_3d_glo,          &    
                                  nbp_lon*(nbp_lat)*nbp_lev,ndex3d)                  
                                                                            
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call gather( dtrconv,auxklonnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklonnbtr_glo(1,it) ,zx_tmp_2d_glo)           
     CALL histwrite(nid_tra2,"dtrconv"//str2,itra,zx_tmp_2d_glo,            & 
                    nbp_lon*(nbp_lat),ndex2d)                                 
!                                                                        
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
     call gather( his_th, auxklonnbtr_glo )
! !$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklonnbtr_glo(1,it) ,zx_tmp_2d_glo)           
     CALL histwrite(nid_tra2,"dtherm"//str2,itra,zx_tmp_2d_glo,             &  
                    nbp_lon*(nbp_lat),ndex2d)                                 
      ENDIF ! mpi root
! !$OMP END MASTER
! !$OMP BARRIER
!                                                                        
                                                                         
     ENDDO                                                               
!
!$OMP END MASTER
!$OMP BARRIER
     call gather( sed_ss, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)
     CALL histwrite(nid_tra2,"sed_ss",itra,zx_tmp_2d_glo,                & 
                    nbp_lon*(nbp_lat),ndex2d)                              
!                                                                     
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( sed_dust, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)            
     CALL histwrite(nid_tra2,"sed_dust",itra,zx_tmp_2d_glo,               & 
                    nbp_lon*(nbp_lat),ndex2d)                               
!                                                                      
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( sed_dustsco, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)          
     CALL histwrite(nid_tra2,"sed_dustsco",itra,zx_tmp_2d_glo,              &
                    nbp_lon*(nbp_lat),ndex2d)                                 
!                                                                        
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( his_g2pgas, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)             
     CALL histwrite(nid_tra2,"g2p_gas",itra,zx_tmp_2d_glo,                   & 
                    nbp_lon*(nbp_lat),ndex2d)                                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( his_g2paer, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)              
     CALL histwrite(nid_tra2,"g2p_aer",itra,zx_tmp_2d_glo,                   &  
                    nbp_lon*(nbp_lat),ndex2d)                                  
! SAVING VARIABLES IN HISTRAC                                              
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxbb, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)                  
      CALL histwrite(nid_tra1,"fluxbb",itra,zx_tmp_2d_glo,                   & 
                                    nbp_lon*(nbp_lat),ndex2d)                 
!                                                                        
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)                 
     CALL histwrite(nid_tra1,"fluxff",itra,zx_tmp_2d_glo,                   &  
                                    nbp_lon*(nbp_lat),ndex2d)                 
!                                                                        
! ======================== BC =============================              
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxbcbb, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxbcbb",itra,zx_tmp_2d_glo,                 & 
                                    nbp_lon*(nbp_lat),ndex2d)                 
!                                                                        
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxbcff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxbcff",itra,zx_tmp_2d_glo,                  & 
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxbcnff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxbcnff",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxbcba, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)                
     CALL histwrite(nid_tra1,"fluxbcba",itra,zx_tmp_2d_glo,                  &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxbc, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)                  
     CALL histwrite(nid_tra1,"fluxbc",itra,zx_tmp_2d_glo,                    &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
! ======================== OM =============================               
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxombb, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)                
     CALL histwrite(nid_tra1,"fluxombb",itra,zx_tmp_2d_glo,                  &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxomff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)                
     CALL histwrite(nid_tra1,"fluxomff",itra,zx_tmp_2d_glo,                  &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxomnff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxomnff",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxomba, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)                
     CALL histwrite(nid_tra1,"fluxomba",itra,zx_tmp_2d_glo,                  &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxomnat, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxomnat",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxom, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)                  
     CALL histwrite(nid_tra1,"fluxom",itra,zx_tmp_2d_glo,                    &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
! ======================== SO4 =============================              
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso4ff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxso4ff",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso4nff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)              
     CALL histwrite(nid_tra1,"fluxso4nff",itra,zx_tmp_2d_glo,                &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso4bb, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxso4bb",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso4ba, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxso4ba",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso4, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo ,zx_tmp_2d_glo)                 
     CALL histwrite(nid_tra1,"fluxso4",itra,zx_tmp_2d_glo,                   &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
! ======================== H2S =============================              
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxh2sff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxh2sff",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                  
!                                                                         
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxh2snff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)              
     CALL histwrite(nid_tra1,"fluxh2snff",itra,zx_tmp_2d_glo,                 & 
                                    nbp_lon*(nbp_lat),ndex2d)                   
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxh2sbio, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxh2sbio",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                   
! ======================== SO2 =============================               
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso2ff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                
     CALL histwrite(nid_tra1,"fluxso2ff",itra,zx_tmp_2d_glo,                  &  
                                    nbp_lon*(nbp_lat),ndex2d)                   
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso2nff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxso2nff",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                   
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso2bb, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                
     CALL histwrite(nid_tra1,"fluxso2bb",itra,zx_tmp_2d_glo,                  &  
                                    nbp_lon*(nbp_lat),ndex2d)                   
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso2vol, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxso2vol",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                   
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso2ba, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                
     CALL histwrite(nid_tra1,"fluxso2ba",itra,zx_tmp_2d_glo,                  &  
                                    nbp_lon*(nbp_lat),ndex2d)                   
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxso2, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                  
     CALL histwrite(nid_tra1,"fluxso2",itra,zx_tmp_2d_glo,                    &  
                                    nbp_lon*(nbp_lat),ndex2d)                   
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxdms, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                  
     CALL histwrite(nid_tra1,"fluxdms",itra,zx_tmp_2d_glo,                    &  
                                    nbp_lon*(nbp_lat),ndex2d)                   
! ======================== DD =============================                
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxdustec, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxdustec",itra,zx_tmp_2d_glo,                 &  
                                    nbp_lon*(nbp_lat),ndex2d)                   
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxddfine, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxddfine",itra,zx_tmp_2d_glo,                 &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxddcoa, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                
     CALL histwrite(nid_tra1,"fluxddcoa",itra,zx_tmp_2d_glo,                  &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxddsco, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                
     CALL histwrite(nid_tra1,"fluxddsco",itra,zx_tmp_2d_glo,                  &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxdd, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                   
     CALL histwrite(nid_tra1,"fluxdd",itra,zx_tmp_2d_glo,                     &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
! ======================== SS =============================                
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxssfine, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)               
     CALL histwrite(nid_tra1,"fluxssfine",itra,zx_tmp_2d_glo,                 &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxsscoa, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                
     CALL histwrite(nid_tra1,"fluxsscoa",itra,zx_tmp_2d_glo,                  &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( fluxss, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                   
     CALL histwrite(nid_tra1,"fluxss",itra,zx_tmp_2d_glo,                     &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER

!     call gather( , auxklon_glo )
!!!!      IF (is_mpi_root .AND. is_omp_root) THEN
!nhl     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,fluxso4chem,zx_tmp_3d_glo)       
!nhl     CALL histwrite(nid_tra1,"fluxso4chem",itra,zx_tmp_3d_glo,            &  
!nhl    .                             nbp_lon*(nbp_lat)*nbp_lev,ndex3d)            
!                                                                          
     call gather( flux_sparam_ind, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)          
     CALL histwrite(nid_tra1,"flux_sparam_ind",itra,zx_tmp_2d_glo,            &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( flux_sparam_bb, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)           
     CALL histwrite(nid_tra1,"flux_sparam_bb",itra,zx_tmp_2d_glo,             &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( flux_sparam_ff, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)           
     CALL histwrite(nid_tra1,"flux_sparam_ff",itra,zx_tmp_2d_glo,             &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( flux_sparam_ddfine, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)       
     CALL histwrite(nid_tra1,"flux_sparam_ddfine",itra,zx_tmp_2d_glo,         &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( flux_sparam_ddcoa, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)        
     CALL histwrite(nid_tra1,"flux_sparam_ddcoa",itra,zx_tmp_2d_glo,          &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( flux_sparam_ddsco, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)        
     CALL histwrite(nid_tra1,"flux_sparam_ddsco",itra,zx_tmp_2d_glo,          &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( flux_sparam_ssfine, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)       
     CALL histwrite(nid_tra1,"flux_sparam_ssfine",itra,zx_tmp_2d_glo,         &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( flux_sparam_sscoa, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)        
     CALL histwrite(nid_tra1,"flux_sparam_sscoa",itra,zx_tmp_2d_glo,          &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( u10m_ec, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                  
     CALL histwrite(nid_tra1,"u10m",itra,zx_tmp_2d_glo,                       &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER
     call gather( v10m_ec, auxklon_glo )
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
     CALL gr_fi_ecrit(1,klon_glo,nbp_lon,nbp_lat,auxklon_glo,zx_tmp_2d_glo)                  
     CALL histwrite(nid_tra1,"v10m",itra,zx_tmp_2d_glo,                       &  
                                  nbp_lon*(nbp_lat),ndex2d)                     
!                                                                          
!     call gather( , auxklon_glo )
!!!   !$OMP MASTER
!      IF (is_mpi_root .AND. is_omp_root) THEN
!nhl     CALL gr_fi_ecrit(nbp_lev,klon_glo,nbp_lon,nbp_lat,flux_sparam_sulf,zx_tmp_3d_glo)  
!nhl     CALL histwrite(nid_tra1,"flux_sparam_sulf",itra,zx_tmp_3d_glo,       &  
!nhl    .                             nbp_lon*(nbp_lat)*nbp_lev,ndex3d)            
!                                                                          
      ENDIF ! mpi root
!$OMP END MASTER
!$OMP BARRIER

      ENDIF ! ok_histrac                                                    
                                                                            



!JE20141224
! saving variables for output
! 2D outputs
      DO i=1, klon
       trm01(i)=0. 
       trm02(i)=0. 
       trm03(i)=0. 
       trm04(i)=0. 
       trm05(i)=0.
       sconc01(i)=0. 
       sconc02(i)=0. 
       sconc03(i)=0. 
       sconc04(i)=0. 
       sconc05(i)=0.
       flux01(i)=0. 
       flux02(i)=0. 
       flux03(i)=0. 
       flux04(i)=0. 
       flux05(i)=0.
       ds01(i)=0. 
       ds02(i)=0. 
       ds03(i)=0. 
       ds04(i)=0. 
       ds05(i)=0.
       dh01(i)=0. 
       dh02(i)=0. 
       dh03(i)=0. 
       dh04(i)=0. 
       dh05(i)=0.
       dtrconv01(i)=0. 
       dtrconv02(i)=0. 
       dtrconv03(i)=0. 
       dtrconv04(i)=0. 
       dtrconv05(i)=0.
       dtherm01(i)=0. 
       dtherm02(i)=0. 
       dtherm03(i)=0. 
       dtherm04(i)=0. 
       dtherm05(i)=0.
       dhkecv01(i)=0. 
       dhkecv02(i)=0. 
       dhkecv03(i)=0. 
       dhkecv04(i)=0. 
       dhkecv05(i)=0.
       d_tr_ds01(i)=0. 
       d_tr_ds02(i)=0. 
       d_tr_ds03(i)=0. 
       d_tr_ds04(i)=0. 
       d_tr_ds05(i)=0.
       dhkelsc01(i)=0. 
       dhkelsc02(i)=0. 
       dhkelsc03(i)=0. 
       dhkelsc04(i)=0. 
       dhkelsc05(i)=0.
!       u10m_ss(i)=u10m_ec(i)
!       v10m_ss(i)=v10m_ec(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(id_prec>0)  trm01(i)=trm(i,id_prec) 
      if(id_fine>0)  trm02(i)=trm(i,id_fine) 
      if(id_coss>0)  trm03(i)=trm(i,id_coss) 
      if(id_codu>0)  trm04(i)=trm(i,id_codu) 
      if(id_scdu>0)  trm05(i)=trm(i,id_scdu)
      if(id_prec>0)    sconc01(i)=sconc_seri(i,id_prec)
      if(id_fine>0)    sconc02(i)=sconc_seri(i,id_fine)
      if(id_coss>0)    sconc03(i)=sconc_seri(i,id_coss)
      if(id_codu>0)    sconc04(i)=sconc_seri(i,id_codu)
      if(id_scdu>0)    sconc05(i)=sconc_seri(i,id_scdu)
      if(id_prec>0)    flux01(i)=flux_tr(i,id_prec)
      if(id_fine>0)    flux02(i)=flux_tr(i,id_fine)
      if(id_coss>0)    flux03(i)=flux_tr(i,id_coss)
      if(id_codu>0)    flux04(i)=flux_tr(i,id_codu)
      if(id_scdu>0)    flux05(i)=flux_tr(i,id_scdu)
      if(id_prec>0)    ds01(i)=his_ds(i,id_prec)
      if(id_fine>0)    ds02(i)=his_ds(i,id_fine)
      if(id_coss>0)    ds03(i)=his_ds(i,id_coss)
      if(id_codu>0)    ds04(i)=his_ds(i,id_codu)
      if(id_scdu>0)    ds05(i)=his_ds(i,id_scdu)
      if(id_prec>0)    dh01(i)=his_dh(i,id_prec)
      if(id_fine>0)    dh02(i)=his_dh(i,id_fine)
      if(id_coss>0)    dh03(i)=his_dh(i,id_coss)
      if(id_codu>0)    dh04(i)=his_dh(i,id_codu)
      if(id_scdu>0)    dh05(i)=his_dh(i,id_scdu)
      if(id_prec>0)    dtrconv01(i)=dtrconv(i,id_prec)
      if(id_fine>0)    dtrconv02(i)=dtrconv(i,id_fine)
      if(id_coss>0)    dtrconv03(i)=dtrconv(i,id_coss)
      if(id_codu>0)    dtrconv04(i)=dtrconv(i,id_codu)
      if(id_scdu>0)    dtrconv05(i)=dtrconv(i,id_scdu)
      if(id_prec>0)    dtherm01(i)=his_th(i,id_prec)
      if(id_fine>0)    dtherm02(i)=his_th(i,id_fine)
      if(id_coss>0)    dtherm03(i)=his_th(i,id_coss)
      if(id_codu>0)    dtherm04(i)=his_th(i,id_codu)
      if(id_scdu>0)    dtherm05(i)=his_th(i,id_scdu)
      if(id_prec>0)    dhkecv01(i)=his_dhkecv(i,id_prec)
      if(id_fine>0)    dhkecv02(i)=his_dhkecv(i,id_fine)
      if(id_coss>0)    dhkecv03(i)=his_dhkecv(i,id_coss)
      if(id_codu>0)    dhkecv04(i)=his_dhkecv(i,id_codu)
      if(id_scdu>0)    dhkecv05(i)=his_dhkecv(i,id_scdu)
      if(id_prec>0)    d_tr_ds01(i)=his_ds(i,id_prec)
      if(id_fine>0)    d_tr_ds02(i)=his_ds(i,id_fine)
      if(id_coss>0)    d_tr_ds03(i)=his_ds(i,id_coss)
      if(id_codu>0)    d_tr_ds04(i)=his_ds(i,id_codu)
      if(id_scdu>0)    d_tr_ds05(i)=his_ds(i,id_scdu)
      if(id_prec>0)    dhkelsc01(i)=his_dhkelsc(i,id_prec)
      if(id_fine>0)    dhkelsc02(i)=his_dhkelsc(i,id_fine)
      if(id_coss>0)    dhkelsc03(i)=his_dhkelsc(i,id_coss)
      if(id_codu>0)    dhkelsc04(i)=his_dhkelsc(i,id_codu)
      if(id_scdu>0)    dhkelsc05(i)=his_dhkelsc(i,id_scdu)
       u10m_ss(i)=u10m_ec(i)
       v10m_ss(i)=v10m_ec(i)
      ENDDO
! 3D outs
      DO i=1, klon
        DO k=1,klev
      d_tr_cv01(i,k)   =0.
      d_tr_cv02(i,k)   =0.
      d_tr_cv03(i,k)   =0.
      d_tr_cv04(i,k)   =0.
      d_tr_cv05(i,k)   =0.
      d_tr_trsp01(i,k) =0.
      d_tr_trsp02(i,k) =0.
      d_tr_trsp03(i,k) =0.
      d_tr_trsp04(i,k) =0.
      d_tr_trsp05(i,k) =0.
      d_tr_sscav01(i,k)=0.
      d_tr_sscav02(i,k)=0.
      d_tr_sscav03(i,k)=0.
      d_tr_sscav04(i,k)=0.
      d_tr_sscav05(i,k)=0.
      d_tr_sat01(i,k)  =0.
      d_tr_sat02(i,k)  =0.
      d_tr_sat03(i,k)  =0.
      d_tr_sat04(i,k)  =0.
      d_tr_sat05(i,k)  =0.
      d_tr_uscav01(i,k)=0.
      d_tr_uscav02(i,k)=0.
      d_tr_uscav03(i,k)=0.
      d_tr_uscav04(i,k)=0.
      d_tr_uscav05(i,k)=0.
      d_tr_insc01(i,k)=0.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      d_tr_insc02(i,k)=0.
      d_tr_insc03(i,k)=0.
      d_tr_insc04(i,k)=0.
      d_tr_insc05(i,k)=0.
      d_tr_bcscav01(i,k)=0.
      d_tr_bcscav02(i,k)=0.
      d_tr_bcscav03(i,k)=0.
      d_tr_bcscav04(i,k)=0.
      d_tr_bcscav05(i,k)=0.
      d_tr_evapls01(i,k)=0.
      d_tr_evapls02(i,k)=0.
      d_tr_evapls03(i,k)=0.
      d_tr_evapls04(i,k)=0.
      d_tr_evapls05(i,k)=0.
      d_tr_ls01(i,k)=0.
      d_tr_ls02(i,k)=0.
      d_tr_ls03(i,k)=0.
      d_tr_ls04(i,k)=0.
      d_tr_ls05(i,k)=0.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      d_tr_dyn01(i,k)=0.
      d_tr_dyn02(i,k)=0.
      d_tr_dyn03(i,k)=0.
      d_tr_dyn04(i,k)=0.
      d_tr_dyn05(i,k)=0.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      d_tr_cl01(i,k)=0.
      d_tr_cl02(i,k)=0.
      d_tr_cl03(i,k)=0.
      d_tr_cl04(i,k)=0.
      d_tr_cl05(i,k)=0.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      d_tr_th01(i,k)=0.
      d_tr_th02(i,k)=0.
      d_tr_th03(i,k)=0.
      d_tr_th04(i,k)=0.
      d_tr_th05(i,k)=0.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(id_prec>0)        d_tr_cv01(i,k)   =d_tr_cv_o(i,k,id_prec)
      if(id_fine>0)        d_tr_cv02(i,k)   =d_tr_cv_o(i,k,id_fine)
      if(id_coss>0)        d_tr_cv03(i,k)   =d_tr_cv_o(i,k,id_coss)
      if(id_codu>0)        d_tr_cv04(i,k)   =d_tr_cv_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_cv05(i,k)   =d_tr_cv_o(i,k,id_scdu)
      if(id_prec>0)        d_tr_trsp01(i,k) =d_tr_trsp_o(i,k,id_prec)
      if(id_fine>0)        d_tr_trsp02(i,k) =d_tr_trsp_o(i,k,id_fine)
      if(id_coss>0)        d_tr_trsp03(i,k) =d_tr_trsp_o(i,k,id_coss)
      if(id_codu>0)        d_tr_trsp04(i,k) =d_tr_trsp_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_trsp05(i,k) =d_tr_trsp_o(i,k,id_scdu)
      if(id_prec>0)        d_tr_sscav01(i,k)=d_tr_sscav_o(i,k,id_prec)
      if(id_fine>0)        d_tr_sscav02(i,k)=d_tr_sscav_o(i,k,id_fine)
      if(id_coss>0)        d_tr_sscav03(i,k)=d_tr_sscav_o(i,k,id_coss)
      if(id_codu>0)        d_tr_sscav04(i,k)=d_tr_sscav_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_sscav05(i,k)=d_tr_sscav_o(i,k,id_scdu)
      if(id_prec>0)        d_tr_sat01(i,k)  =d_tr_sat_o(i,k,id_prec)
      if(id_fine>0)        d_tr_sat02(i,k)  =d_tr_sat_o(i,k,id_fine)
      if(id_coss>0)        d_tr_sat03(i,k)  =d_tr_sat_o(i,k,id_coss)
      if(id_codu>0)        d_tr_sat04(i,k)  =d_tr_sat_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_sat05(i,k)  =d_tr_sat_o(i,k,id_scdu)
      if(id_prec>0)        d_tr_uscav01(i,k)=d_tr_uscav_o(i,k,id_prec)
      if(id_fine>0)        d_tr_uscav02(i,k)=d_tr_uscav_o(i,k,id_fine)
      if(id_coss>0)        d_tr_uscav03(i,k)=d_tr_uscav_o(i,k,id_coss)
      if(id_codu>0)        d_tr_uscav04(i,k)=d_tr_uscav_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_uscav05(i,k)=d_tr_uscav_o(i,k,id_scdu)
      if(id_prec>0)        d_tr_insc01(i,k)=d_tr_insc_o(i,k,id_prec)
      if(id_fine>0)        d_tr_insc02(i,k)=d_tr_insc_o(i,k,id_fine)
      if(id_coss>0)        d_tr_insc03(i,k)=d_tr_insc_o(i,k,id_coss)
      if(id_codu>0)        d_tr_insc04(i,k)=d_tr_insc_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_insc05(i,k)=d_tr_insc_o(i,k,id_scdu)
      if(id_prec>0)        d_tr_bcscav01(i,k)=d_tr_bcscav_o(i,k,id_prec)
      if(id_fine>0)        d_tr_bcscav02(i,k)=d_tr_bcscav_o(i,k,id_fine)
      if(id_coss>0)        d_tr_bcscav03(i,k)=d_tr_bcscav_o(i,k,id_coss)
      if(id_codu>0)        d_tr_bcscav04(i,k)=d_tr_bcscav_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_bcscav05(i,k)=d_tr_bcscav_o(i,k,id_scdu)
      if(id_prec>0)        d_tr_evapls01(i,k)=d_tr_evapls_o(i,k,id_prec)
      if(id_fine>0)        d_tr_evapls02(i,k)=d_tr_evapls_o(i,k,id_fine)
      if(id_coss>0)        d_tr_evapls03(i,k)=d_tr_evapls_o(i,k,id_coss)
      if(id_codu>0)        d_tr_evapls04(i,k)=d_tr_evapls_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_evapls05(i,k)=d_tr_evapls_o(i,k,id_scdu)
        ENDDO
      ENDDO
      IF(1==0) THEN
      DO i=1, klon
        DO k=1,klev
      if(id_prec>0)        d_tr_ls01(i,k)=d_tr_ls_o(i,k,id_prec)
      if(id_fine>0)        d_tr_ls02(i,k)=d_tr_ls_o(i,k,id_fine)
      if(id_coss>0)        d_tr_ls03(i,k)=d_tr_ls_o(i,k,id_coss)
      if(id_codu>0)        d_tr_ls04(i,k)=d_tr_ls_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_ls05(i,k)=d_tr_ls_o(i,k,id_scdu)
      if(id_prec>0)        d_tr_dyn01(i,k)=d_tr_dyn_o(i,k,id_prec)
      if(id_fine>0)        d_tr_dyn02(i,k)=d_tr_dyn_o(i,k,id_fine)
      if(id_coss>0)        d_tr_dyn03(i,k)=d_tr_dyn_o(i,k,id_coss)
      if(id_codu>0)        d_tr_dyn04(i,k)=d_tr_dyn_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_dyn05(i,k)=d_tr_dyn_o(i,k,id_scdu)
      if(id_prec>0)        d_tr_cl01(i,k)=d_tr_cl_o(i,k,id_prec)
      if(id_fine>0)        d_tr_cl02(i,k)=d_tr_cl_o(i,k,id_fine)
      if(id_coss>0)        d_tr_cl03(i,k)=d_tr_cl_o(i,k,id_coss)
      if(id_codu>0)        d_tr_cl04(i,k)=d_tr_cl_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_cl05(i,k)=d_tr_cl_o(i,k,id_scdu)
      if(id_prec>0)        d_tr_th01(i,k)=d_tr_th_o(i,k,id_prec)
      if(id_fine>0)        d_tr_th02(i,k)=d_tr_th_o(i,k,id_fine)
      if(id_coss>0)        d_tr_th03(i,k)=d_tr_th_o(i,k,id_coss)
      if(id_codu>0)        d_tr_th04(i,k)=d_tr_th_o(i,k,id_codu)
      if(id_scdu>0)        d_tr_th05(i,k)=d_tr_th_o(i,k,id_scdu)
        ENDDO
      ENDDO
      ELSE
      DO i=1, klon
        DO k=1,klev
      if(id_prec>0)        d_tr_ls01(i,k)=d_tr_ls(i,k,id_prec)/pdtphys
      if(id_fine>0)        d_tr_ls02(i,k)=d_tr_ls(i,k,id_fine)/pdtphys
      if(id_coss>0)        d_tr_ls03(i,k)=d_tr_ls(i,k,id_coss)/pdtphys
      if(id_codu>0)        d_tr_ls04(i,k)=d_tr_ls(i,k,id_codu)/pdtphys
      if(id_scdu>0)        d_tr_ls05(i,k)=d_tr_ls(i,k,id_scdu)/pdtphys
      if(id_prec>0)        d_tr_dyn01(i,k)=d_tr_dyn(i,k,id_prec)/pdtphys
      if(id_fine>0)        d_tr_dyn02(i,k)=d_tr_dyn(i,k,id_fine)/pdtphys
      if(id_coss>0)        d_tr_dyn03(i,k)=d_tr_dyn(i,k,id_coss)/pdtphys
      if(id_codu>0)        d_tr_dyn04(i,k)=d_tr_dyn(i,k,id_codu)/pdtphys
      if(id_scdu>0)        d_tr_dyn05(i,k)=d_tr_dyn(i,k,id_scdu)/pdtphys
      if(id_prec>0)        d_tr_cl01(i,k)=d_tr_cl(i,k,id_prec)/pdtphys
      if(id_fine>0)        d_tr_cl02(i,k)=d_tr_cl(i,k,id_fine)/pdtphys
      if(id_coss>0)        d_tr_cl03(i,k)=d_tr_cl(i,k,id_coss)/pdtphys
      if(id_codu>0)        d_tr_cl04(i,k)=d_tr_cl(i,k,id_codu)/pdtphys
      if(id_scdu>0)        d_tr_cl05(i,k)=d_tr_cl(i,k,id_scdu)/pdtphys
      if(id_prec>0)        d_tr_th01(i,k)=d_tr_th(i,k,id_prec)/pdtphys
      if(id_fine>0)        d_tr_th02(i,k)=d_tr_th(i,k,id_fine)/pdtphys
      if(id_coss>0)        d_tr_th03(i,k)=d_tr_th(i,k,id_coss)/pdtphys
      if(id_codu>0)        d_tr_th04(i,k)=d_tr_th(i,k,id_codu)/pdtphys
      if(id_scdu>0)        d_tr_th05(i,k)=d_tr_th(i,k,id_scdu)/pdtphys
        ENDDO
      ENDDO
      ENDIF
     

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)

      dife=clock_end-clock_start
      ti_outs=dife*MAX(0,SIGN(1,dife))   &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_outs=tia_outs+REAL(ti_outs)/REAL(clock_rate)
      ENDIF

      IF (logitime) THEN
      CALL SYSTEM_CLOCK(COUNT=clock_end)

      dife=clock_end-clock_start_spla
      ti_spla=dife*MAX(0,SIGN(1,dife)) &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))


      tia_spla=tia_spla+REAL(ti_spla)/REAL(clock_rate)
  print *,'times for this timestep:timeproc,timeproc/time_pytracr_spl-'
     print *,'time spla',REAL(ti_spla)/REAL(clock_rate)                &
      ,REAL(ti_spla)/REAL(ti_spla)                                      
     print *,'time init',REAL(ti_init)/REAL(clock_rate)                &
      ,REAL(ti_init)/REAL(ti_spla)                                      
     print *,'time inittype',REAL(ti_inittype)/REAL(clock_rate)        &
      ,REAL(ti_inittype)/REAL(ti_spla)                                  
     print *,'time inittwrite',REAL(ti_inittwrite)/REAL(clock_rate)    &
      ,REAL(ti_inittwrite)/REAL(ti_spla)                                
     print *,'time emis',REAL(ti_emis)/REAL(clock_rate)                &
      ,REAL(ti_emis)/REAL(ti_spla)                                      
     print *,'time depo ',REAL(ti_depo)/REAL(clock_rate)               &
      ,REAL(ti_depo)/REAL(ti_spla)                                      
     print *,'time cltr',REAL(ti_cltr)/REAL(clock_rate)                &
      ,REAL(ti_cltr)/REAL(ti_spla)                                      
     print *,'time ther',REAL(ti_ther)/REAL(clock_rate)                &
      ,REAL(ti_ther)/REAL(ti_spla)                                      
     print *,'time sedi',REAL(ti_sedi)/REAL(clock_rate)                &
      ,REAL(ti_sedi)/REAL(ti_spla)                                      
     print *,'time gas to part',REAL(ti_gasp)/REAL(clock_rate)         &
      ,REAL(ti_gasp)/REAL(ti_spla)                                      
     print *,'time AP wet',REAL(ti_wetap)/REAL(clock_rate)             &
      ,REAL(ti_wetap)/REAL(ti_spla)                                     
     print *,'time convective',REAL(ti_cvltr)/REAL(clock_rate)         &
      ,REAL(ti_cvltr)/REAL(ti_spla)                                     
     print *,'time NP lsc scav',REAL(ti_lscs)/REAL(clock_rate)         &
      ,REAL(ti_lscs)/REAL(ti_spla)                                      
     print *,'time opt,brdn,etc',REAL(ti_brop)/REAL(clock_rate)        &
      ,REAL(ti_brop)/REAL(ti_spla)                                      
     print *,'time outputs',REAL(ti_outs)/REAL(clock_rate)             &
      ,REAL(ti_outs)/REAL(ti_spla)


  print *,'--time accumulated: time proc, time proc/time phytracr_spl--'
      print *,'time spla',tia_spla
      print *,'time init',tia_init,tia_init/tia_spla
      print *,'time inittype',tia_inittype,tia_inittype/tia_spla
      print *,'time inittwrite',tia_inittwrite,tia_inittwrite/tia_spla
      print *,'time emis',tia_emis,tia_emis/tia_spla
      print *,'time depo',tia_depo,tia_depo/tia_spla
      print *,'time cltr',tia_cltr,tia_cltr/tia_spla
      print *,'time ther',tia_ther,tia_ther/tia_spla
      print *,'time sedi',tia_sedi,tia_sedi/tia_spla
      print *,'time gas to part',tia_gasp,tia_gasp/tia_spla
      print *,'time AP wet',tia_wetap,tia_wetap/tia_spla
      print *,'time convective',tia_cvltr,tia_cvltr/tia_spla
      print *,'time NP lsc scav',tia_lscs,tia_lscs/tia_spla
      print *,'time opt,brdn,etc',tia_brop,tia_brop/tia_spla
      print *,'time outputs',tia_outs,tia_outs/tia_spla



      dife=clock_end_outphytracr-clock_start_outphytracr
      ti_nophytracr=dife*MAX(0,SIGN(1,dife))  &
      +(dife+clock_per_max)*MAX(0,SIGN(1,-dife))
      tia_nophytracr=tia_nophytracr+REAL(ti_nophytracr)/REAL(clock_rate)
      print *,'Time outside phytracr; Time accum outside phytracr'
      print*,REAL(ti_nophytracr)/REAL(clock_rate),tia_nophytracr

      clock_start_outphytracr=clock_end

      ENDIF      
      print *,'END PHYTRACR_SPL '
  print *,'lmt_so2ff_l FIN' , MINVAL(lmt_so2ff_l), MAXVAL(lmt_so2ff_l)

!      CALL abort_gcm('TEST1', 'OK1', 1)

      RETURN
      END SUBROUTINE phytracr_spl
 
      SUBROUTINE readregionsdims2_spl(nbreg,fileregions)

      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para

      IMPLICIT NONE
      CHARACTER*800 fileregions
      CHARACTER*800 auxstr
      INTEGER nbreg
 
      IF (is_mpi_root .AND. is_omp_root) THEN 

      OPEN (UNIT=1,FILE=trim(adjustl(fileregions)))
      READ(1,'(a)') auxstr
      READ(1,'(i10)') nbreg
      CLOSE(UNIT=1)
      ENDIF
      CALL bcast(nbreg)

      END SUBROUTINE readregionsdims2_spl

      SUBROUTINE readregionsdims_spl(nbreg_ind,fileregionsdimsind,   &
                                    nbreg_dust,fileregionsdimsdust,  &
                                    nbreg_bb,fileregionsdimsbb)     
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para

      IMPLICIT NONE
      CHARACTER*800 fileregionsdimsind
      CHARACTER*800 fileregionsdimsdust
      CHARACTER*800 fileregionsdimsbb
      CHARACTER*800 auxstr
      INTEGER nbreg_ind,nbreg_dust,nbreg_bb
 
      IF (is_mpi_root .AND. is_omp_root) THEN 

      OPEN (UNIT=1,FILE=trim(adjustl(fileregionsdimsind)))
      READ(1,'(a)') auxstr
      READ(1,'(i10)') nbreg_ind
      CLOSE(UNIT=1)

      OPEN (UNIT=1,FILE=trim(adjustl(fileregionsdimsdust)))
      READ(1,'(a)') auxstr
      READ(1,'(i10)') nbreg_dust
      CLOSE(UNIT=1)

      OPEN (UNIT=1,FILE=trim(adjustl(fileregionsdimsbb)))
      READ(1,'(a)') auxstr
      READ(1,'(i10)') nbreg_bb
      CLOSE(UNIT=1)
      

      ENDIF
      CALL bcast(nbreg_ind)
      CALL bcast(nbreg_dust)
      CALL bcast(nbreg_bb)

      END SUBROUTINE readregionsdims_spl

      SUBROUTINE readregions_spl(iregion,filenameregion)
      USE dimphy
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para

      IMPLICIT NONE
      CHARACTER*(*) filenameregion
      INTEGER iregion(klon)
      INTEGER iregion_glo(klon_glo)
      INTEGER k
     
      IF (is_mpi_root .AND. is_omp_root) THEN

      print *,trim(adjustl(filenameregion))
      OPEN(1,file=trim(adjustl(filenameregion)))
      DO k=1,klon_glo
      READ(1,'(i10)') iregion_glo(k)
      ENDDO
      CLOSE(UNIT=1)
      ENDIF
      CALL scatter(iregion_glo,iregion)

      END SUBROUTINE readregions_spl

      SUBROUTINE readscaleparams_spl(scale_param, nbreg, &
                                             filescaleparams)
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      IMPLICIT NONE

      CHARACTER*800 filescaleparams
      INTEGER nbreg
      REAL scale_param(nbreg)
      INTEGER k

      IF (is_mpi_root .AND. is_omp_root) THEN
      OPEN(1,file=trim(adjustl(filescaleparams)),form='unformatted')
      do k=1,nbreg
        read(1)  scale_param(k)
      enddo
      CLOSE(1)  
      ENDIF
      CALL bcast(scale_param)
!      print *,'holaaaaaaaaaaaa'
!      print *,scale_param

      END SUBROUTINE readscaleparams_spl

      SUBROUTINE readscaleparamsnc_spl(scale_param_ind,                 &
        nbreg_ind, paramname_ind,                                       &
        scale_param_ff, nbreg_ff,paramname_ff,                          &
        scale_param_bb, nbreg_bb,paramname_bb,                          &
        scale_param_dustacc, nbreg_dustacc,paramname_dustacc,           &
        scale_param_dustcoa, nbreg_dustcoa,paramname_dustcoa,           &
        scale_param_dustsco, nbreg_dustsco,paramname_dustsco,           &
        param_wstarBLperregion, nbreg_wstardustBL, paramname_wstarBL,     &
        param_wstarWAKEperregion, nbreg_wstardustWAKE, paramname_wstarWAKE, &
        scale_param_ssacc  ,  paramname_ssacc,             &
        scale_param_sscoa  ,  paramname_sscoa,             &
           filescaleparams,julien,jH_phys, pdtphys,debutphy)
!      SUBROUTINE readscaleparamsnc_spl(scale_param, nbreg, &
!                                        filescaleparams,paramname,&
!                                        julien,jH_phys, pdtphys,debutphy)
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      IMPLICIT NONE

      CHARACTER*800 filescaleparams
      CHARACTER*100 paramname_ind,paramname_ff,paramname_bb          
      CHARACTER*100 paramname_dustacc, paramname_dustcoa 
      CHARACTER*100 paramname_dustsco
      CHARACTER*100 paramname_ssacc
      CHARACTER*100 paramname_sscoa
      CHARACTER*100 paramname_wstarBL
      CHARACTER*100 paramname_wstarWAKE
     
      INTEGER nbreg,iday
      INTEGER nbreg_ind, nbreg_ff, nbreg_bb , nbreg_dustacc
      INTEGER nbreg_dustcoa , nbreg_dustsco, nbreg_wstardustBL
      INTEGER  nbreg_wstardustWAKE
      INTEGER,PARAMETER ::  nbreg_ssacc=1 
      INTEGER,PARAMETER :: nbreg_sscoa=1 
      REAL,PARAMETER :: sca_resol = 24. ! resolution of scalig params in hours
      REAL scale_param_ind(nbreg_ind)
      REAL scale_param_bb(nbreg_bb)
      REAL scale_param_ff(nbreg_ff)
      REAL scale_param_dustacc(nbreg_dustacc)
      REAL scale_param_dustcoa(nbreg_dustcoa)
      REAL scale_param_dustsco(nbreg_dustsco)
      REAL param_wstarBLperregion(nbreg_wstardustBL)
      REAL param_wstarWAKEperregion(nbreg_wstardustWAKE)
      REAL scale_param_ssacc
      REAL scale_param_ssacc_tmp(nbreg_ssacc)
      REAL scale_param_sscoa
      REAL scale_param_sscoa_tmp(nbreg_sscoa)

      INTEGER k,step_sca,test_sca
      REAL :: jH_phys,  pdtphys
      REAL,SAVE :: jH_sca, jH_ini
      INTEGER julien
      LOGICAL debutphy
      SAVE step_sca,test_sca,iday
!$OMP THREADPRIVATE(step_sca,test_sca,iday)
!$OMP THREADPRIVATE(jH_sca,jH_ini)

      IF (debutphy) THEN
        iday=julien
        step_sca=1
        test_sca=0   
        jH_ini=jH_phys
        jH_sca=jH_phys
      ENDIF

      IF (test_sca .EQ. 0 ) THEN
        ! READ file!!
        call read_scalenc(filescaleparams,paramname_ind,            &
                           nbreg_ind,step_sca,                      &
                           scale_param_ind)
        call read_scalenc(filescaleparams,paramname_bb,            &
                           nbreg_bb,step_sca,                      &
                           scale_param_bb)
        call read_scalenc(filescaleparams,paramname_ff,            &
                           nbreg_ff,step_sca,                      &
                           scale_param_ff)
        call read_scalenc(filescaleparams,paramname_dustacc,            &
                           nbreg_dustacc,step_sca,                      &
                           scale_param_dustacc)
        call read_scalenc(filescaleparams,paramname_dustcoa,            &
                           nbreg_dustcoa,step_sca,                      &
                           scale_param_dustcoa)
        call read_scalenc(filescaleparams,paramname_dustsco,            &
                           nbreg_dustsco,step_sca,                      &
                           scale_param_dustsco)
        call read_scalenc(filescaleparams,paramname_wstarBL,            &
                           nbreg_wstardustBL,step_sca,                    &
                           param_wstarBLperregion)
        call read_scalenc(filescaleparams,paramname_wstarWAKE,          &
                           nbreg_wstardustWAKE,step_sca,                    &
                           param_wstarWAKEperregion)
        call read_scalenc(filescaleparams,paramname_ssacc,              &
                           nbreg_ssacc,step_sca,                        &
                           scale_param_ssacc_tmp)
        call read_scalenc(filescaleparams,paramname_sscoa,              &
                           nbreg_sscoa,step_sca,                        &
                           scale_param_sscoa_tmp)
         scale_param_ssacc=scale_param_ssacc_tmp(1)
         scale_param_sscoa=scale_param_sscoa_tmp(1)

       !print *,'JEREADFILE',julien,jH_phys
        step_sca= step_sca + 1
        test_sca=1
      ENDIF

      jH_sca=jH_sca+pdtphys/(24.*3600.)
      IF (jH_sca.GT.(sca_resol)/24.) THEN
          test_sca=0
          jH_sca=jH_ini
      ENDIF

      END SUBROUTINE readscaleparamsnc_spl

      SUBROUTINE read_scalenc(filescaleparams,paramname,nbreg,step_sca, &
                          scale_param)

      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      IMPLICIT NONE

      include "netcdf.inc"

      CHARACTER*800 filescaleparams
      CHARACTER*100 paramname
      INTEGER nbreg, step_sca
      REAL scale_param(nbreg)
      !local vars
      integer nid,ierr,nvarid
      real rcode,auxreal
      integer start(4),count(4), status
!      local
      integer debutread,countread
      CHARACTER*104 varname
      CHARACTER*2 aux_2s
      integer i, j, ig
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
          !nci=NCOPN(trim(adjustl(filescaleparams)),NCNOWRIT,rcode)
         ierr = NF_OPEN (trim(adjustl(filescaleparams)),NF_NOWRITE, nid)
          if (ierr .EQ. NF_NOERR) THEN
          debutread=step_sca
          countread=1

           do i=1,nbreg
            WRITE(aux_2s,'(i2.2)') i
            varname= trim(adjustl(paramname))//aux_2s
            print *,varname
            ierr = NF_INQ_VARID (nid,trim(adjustl(varname)), nvarid)
            ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debutread,          &
                         countread, auxreal)
            IF (ierr .NE. NF_NOERR) THEN
             PRINT*, 'Pb de lecture pour modvalues'
       print *,'JE  scale_var, step_sca',trim(adjustl(varname)),step_sca
             CALL HANDLE_ERR(ierr)
             print *,'error ierr= ',ierr
             CALL exit(1) 
            call abort_gcm('read_scalenc','error reading variable',1)
      ENDIF

            print *,auxreal
            scale_param(i)=auxreal 
           enddo

            ierr = NF_CLOSE(nid)
          else
           print *,'File '//trim(adjustl(filescaleparams))//' not found'
            print *,'doing nothing...'
          endif

      ENDIF ! mpi_root
!$OMP END MASTER
!$OMP BARRIER
!      CALL scatter(var local _glo,var local) o algo asi
      call bcast(scale_param)
      END SUBROUTINE read_scalenc


      
      END MODULE
