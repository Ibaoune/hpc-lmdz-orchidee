!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Module pour construir et detruire les variables "cospIN", "cospstateIN" et "cospOUT"
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE LMDZ_COSP_CONSTRUCT_DESTROY_MOD

use cosp_kinds, only: wp
use mod_cosp,   only: cosp_optical_inputs,cosp_column_inputs,cosp_outputs
use lmdz_cosp_read_outputkeys, only: cosp_config

   USE mod_cosp_config, only : R_UNDEF, CLOUDSAT_DBZE_BINS, SR_BINS, PARASOL_NREFL, & 
                              LIDAR_NTEMP, LIDAR_NCAT, LIDAR_NTYPE, & 
                              numMODISReffIceBins,reffICE_binCenters, &
                              numMODISReffLiqBins, reffLIQ_binCenters, &
			      numISCCPTauBins,numISCCPPresBins, &
			      numMISRTauBins,numMISRHgtBins, &
			      numModisTauBins,numMODISPresBins

  implicit none

contains


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospIN(cfg,npoints,ncolumns,nlevels,y)
    ! Inputs
    integer,intent(in) :: &
         npoints,  & ! Number of horizontal gridpoints
         ncolumns, & ! Number of subcolumns
         nlevels     ! Number of vertical levels

     type(cosp_config),intent(in) :: &
          cfg         ! COSP output config (logical keys)

    ! Outputs 
    type(cosp_optical_inputs),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
    y%Npart    = 4
    y%Nrefl    = PARASOL_NREFL
    allocate(y%frac_out(npoints,       ncolumns,nlevels))

    if (cfg%Lmodis .or. cfg%Lmisr .or. cfg%Lisccp) then
       allocate(y%tau_067(npoints,        ncolumns,nlevels),&
                y%emiss_11(npoints,       ncolumns,nlevels))
    endif
    if (cfg%Lcalipso) then
       allocate(y%betatot_calipso(npoints,        ncolumns,nlevels),&
                y%betatot_ice_calipso(npoints,    ncolumns,nlevels),&
                y%betatot_liq_calipso(npoints,    ncolumns,nlevels),&
                y%tautot_calipso(npoints,         ncolumns,nlevels),&
                y%tautot_ice_calipso(npoints,     ncolumns,nlevels),&
                y%tautot_liq_calipso(npoints,     ncolumns,nlevels),&
                y%beta_mol_calipso(npoints,                nlevels),&
                y%tau_mol_calipso(npoints,                 nlevels),&
                y%tautot_S_ice(npoints,   ncolumns        ),&
                y%tautot_S_liq(npoints,   ncolumns        ))
    endif

    if (cfg%LgrLidar532) then
       allocate(y%beta_mol_grLidar532(npoints,          nlevels),& 
                y%betatot_grLidar532(npoints,  ncolumns,nlevels),& 
                y%tau_mol_grLidar532(npoints,           nlevels),& 
                y%tautot_grLidar532(npoints,   ncolumns,nlevels)) 
    endif

    if (cfg%Latlid) then
       allocate(y%beta_mol_atlid(npoints,             nlevels),& 
                y%betatot_atlid(npoints,     ncolumns,nlevels),& 
                y%tau_mol_atlid(npoints,              nlevels),& 
                y%tautot_atlid(npoints,      ncolumns,nlevels))
    endif 

    if (cfg%Lcloudsat) then
       allocate(y%z_vol_cloudsat(npoints,  ncolumns,nlevels),&
                y%kr_vol_cloudsat(npoints, ncolumns,nlevels),&
                y%g_vol_cloudsat(npoints,  ncolumns,nlevels),&
                y%fracPrecipIce(npoints,   ncolumns))
    endif
    if (cfg%Lmodis) then
       allocate(y%fracLiq(npoints,        ncolumns,nlevels),&
                y%asym(npoints,           ncolumns,nlevels),&
                y%ss_alb(npoints,         ncolumns,nlevels))
    endif
    

  end subroutine construct_cospIN  
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospstateIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  subroutine construct_cospstateIN(npoints,nlevels,nchan,y)
    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal gridpoints
         nlevels, & ! Number of vertical levels
         nchan      ! Number of channels
    ! Outputs
    type(cosp_column_inputs),intent(out) :: y         
    
    allocate(y%sunlit(npoints),y%skt(npoints),y%land(npoints),y%at(npoints,nlevels),     &
             y%pfull(npoints,nlevels),y%phalf(npoints,nlevels+1),y%qv(npoints,nlevels),  &
             y%o3(npoints,nlevels),y%hgt_matrix(npoints,nlevels),y%u_sfc(npoints),       &
             y%v_sfc(npoints),y%lat(npoints),y%lon(nPoints),y%emis_sfc(nchan),           &
             y%cloudIce(nPoints,nLevels),y%cloudLiq(nPoints,nLevels),y%surfelev(npoints),&
             y%fl_snow(nPoints,nLevels),y%fl_rain(nPoints,nLevels),y%seaice(npoints),    &
             y%tca(nPoints,nLevels),y%hgt_matrix_half(npoints,nlevels+1))

  end subroutine construct_cospstateIN

  ! ######################################################################################
  ! SUBROUTINE construct_cosp_outputs
  !
  ! This subroutine allocates output fields based on input logical flag switches.
  ! ######################################################################################  
  subroutine construct_cosp_outputs(cfg,Npoints,Ncolumns,Nlevels,Nlvgrid,Nchan,x)
     ! Inputs        
     integer,intent(in) :: &
          Npoints,         & ! Number of sampled points
          Ncolumns,        & ! Number of subgrid columns
          Nlevels,         & ! Number of model levels
          Nlvgrid,         & ! Number of levels in L3 stats computation
          Nchan              ! Number of RTTOV channels            
     type(cosp_config),intent(in) :: &
          cfg         ! COSP output config (logical keys)

     ! Outputs
     type(cosp_outputs),intent(out) :: &
          x           ! COSP output structure  
   
     ! ISCCP simulator outputs
    if (cfg%Lboxtauisccp)    allocate(x%isccp_boxtau(Npoints,Ncolumns)) 
    if (cfg%Lboxptopisccp)   allocate(x%isccp_boxptop(Npoints,Ncolumns))
    if (cfg%Lclisccp)        allocate(x%isccp_fq(Npoints,numISCCPTauBins,numISCCPPresBins))
    if (cfg%Lcltisccp)       allocate(x%isccp_totalcldarea(Npoints))
    if (cfg%Lpctisccp)       allocate(x%isccp_meanptop(Npoints))
    if (cfg%Ltauisccp)       allocate(x%isccp_meantaucld(Npoints))
    if (cfg%Lmeantbisccp)    allocate(x%isccp_meantb(Npoints))
    if (cfg%Lmeantbclrisccp) allocate(x%isccp_meantbclr(Npoints))
    if (cfg%Lalbisccp)       allocate(x%isccp_meanalbedocld(Npoints))
    
    ! MISR simulator
    if (cfg%LclMISR) then 
       allocate(x%misr_fq(Npoints,numMISRTauBins,numMISRHgtBins))
       ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
       !        they are still computed. Should probably have a logical to control these
       !        outputs.
       allocate(x%misr_dist_model_layertops(Npoints,numMISRHgtBins))
       allocate(x%misr_meanztop(Npoints))
       allocate(x%misr_cldarea(Npoints))    
    endif
    
    ! MODIS simulator
    if (cfg%Lcltmodis)     allocate(x%modis_Cloud_Fraction_Total_Mean(Npoints))
    if (cfg%Lclwmodis)     allocate(x%modis_Cloud_Fraction_Water_Mean(Npoints))
    if (cfg%Lclimodis)     allocate(x%modis_Cloud_Fraction_Ice_Mean(Npoints))
    if (cfg%Lclhmodis)     allocate(x%modis_Cloud_Fraction_High_Mean(Npoints))
    if (cfg%Lclmmodis)     allocate(x%modis_Cloud_Fraction_Mid_Mean(Npoints))
    if (cfg%Lcllmodis)     allocate(x%modis_Cloud_Fraction_Low_Mean(Npoints))
    if (cfg%Ltautmodis)    allocate(x%modis_Optical_Thickness_Total_Mean(Npoints))
    if (cfg%Ltauwmodis)    allocate(x%modis_Optical_Thickness_Water_Mean(Npoints))
    if (cfg%Ltauimodis)    allocate(x%modis_Optical_Thickness_Ice_Mean(Npoints))
    if (cfg%Ltautlogmodis) allocate(x%modis_Optical_Thickness_Total_LogMean(Npoints))
    if (cfg%Ltauwlogmodis) allocate(x%modis_Optical_Thickness_Water_LogMean(Npoints))
    if (cfg%Ltauilogmodis) allocate(x%modis_Optical_Thickness_Ice_LogMean(Npoints))
    if (cfg%Lreffclwmodis) allocate(x%modis_Cloud_Particle_Size_Water_Mean(Npoints))
    if (cfg%Lreffclimodis) allocate(x%modis_Cloud_Particle_Size_Ice_Mean(Npoints))
    if (cfg%Lpctmodis)     allocate(x%modis_Cloud_Top_Pressure_Total_Mean(Npoints))
    if (cfg%Llwpmodis)     allocate(x%modis_Liquid_Water_Path_Mean(Npoints))
    if (cfg%Liwpmodis)     allocate(x%modis_Ice_Water_Path_Mean(Npoints))
    if (cfg%Lclmodis) then
        allocate(x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nPoints,numModisTauBins,numMODISPresBins))
        allocate(x%modis_Optical_thickness_vs_ReffLIQ(nPoints,numMODISTauBins,numMODISReffLiqBins))   
        allocate(x%modis_Optical_Thickness_vs_ReffICE(nPoints,numMODISTauBins,numMODISReffIceBins))
    endif
    
    ! LIDAR simulator
    if (cfg%LlidarBetaMol532) allocate(x%calipso_beta_mol(Npoints,Nlevels))
    if (cfg%Latb532)          allocate(x%calipso_beta_tot(Npoints,Ncolumns,Nlevels))
    if (cfg%LcfadLidarsr532)  then
        allocate(x%calipso_srbval(SR_BINS+1))
        allocate(x%calipso_cfad_sr(Npoints,SR_BINS,Nlvgrid))
        allocate(x%calipso_betaperp_tot(Npoints,Ncolumns,Nlevels))  
    endif
    if (cfg%Lclcalipso)       allocate(x%calipso_lidarcld(Npoints,Nlvgrid))
    if (cfg%Lclhcalipso .or. cfg%Lclmcalipso .or. cfg%Lcllcalipso .or. cfg%Lcltcalipso) then
        allocate(x%calipso_cldlayer(Npoints,LIDAR_NCAT))
    endif   
    if (cfg%Lclcalipsoice .or. cfg%Lclcalipsoliq .or. cfg%Lclcalipsoun) then
        allocate(x%calipso_lidarcldphase(Npoints,Nlvgrid,6))
    endif
    if (cfg%Lclcalipsotmp .or. cfg%Lclcalipsotmpliq .or. cfg%Lclcalipsoice .or. &
                          cfg%Lclcalipsotmpun .or. cfg%Lclcalipsotmpice) then
        allocate(x%calipso_lidarcldtmp(Npoints,LIDAR_NTEMP,5))
    endif
    if (cfg%Lcllcalipsoice .or. cfg%Lclmcalipsoice .or. cfg%Lclhcalipsoice .or.             &
        cfg%Lcltcalipsoice .or. cfg%Lcllcalipsoliq .or. cfg%Lclmcalipsoliq .or.             &
        cfg%Lclhcalipsoliq .or. cfg%Lcltcalipsoliq .or. cfg%Lcllcalipsoun  .or.             &
        cfg%Lclmcalipsoun  .or. cfg%Lclhcalipsoun  .or. cfg%Lcltcalipsoun) then
        allocate(x%calipso_cldlayerphase(Npoints,LIDAR_NCAT,6))     
    endif
    if (cfg%Lclopaquecalipso .or. cfg%Lclthincalipso .or. cfg%Lclzopaquecalipso) then
        allocate(x%calipso_cldtype(Npoints,LIDAR_NTYPE))
    endif 
    if (cfg%Lclopaquetemp .or. cfg%Lclthintemp .or. cfg%Lclzopaquetemp) then 
        allocate(x%calipso_cldtypetemp(Npoints,LIDAR_NTYPE))  
    endif
    if (cfg%Lclopaquemeanz .or. cfg%Lclthinmeanz) then 
        allocate(x%calipso_cldtypemeanz(Npoints,2))
    endif 
    if (cfg%Lclopaquemeanzse .or. cfg%Lclthinmeanzse .or. cfg%Lclzopaquecalipsose) then 
        allocate(x%calipso_cldtypemeanzse(Npoints,3)) 
    endif 
    if (cfg%Lclthinemis) then 
        allocate(x%calipso_cldthinemis(Npoints))
    endif
    if (cfg%Lclcalipsoopaque .or. cfg%Lclcalipsothin .or. cfg%Lclcalipsozopaque .or. &
             cfg%Lclcalipsoopacity) then 
        allocate(x%calipso_lidarcldtype(Npoints,Nlvgrid,LIDAR_NTYPE+1))
    endif
    ! These 2 outputs are part of the calipso output type, but are not controlled by an 
    ! logical switch in the output namelist, so if all other fields are on, then allocate
    if (cfg%LlidarBetaMol532 .or. cfg%Latb532        .or. cfg%LcfadLidarsr532 .or. cfg%Lclcalipso  .or.  &
        cfg%Lclcalipsoice    .or. cfg%Lclcalipsoliq  .or. cfg%Lclcalipsoun    .or. cfg%Lclcalipso2 .or.  &
        cfg%Lclhcalipso      .or. cfg%Lclmcalipso    .or. cfg%Lcllcalipso     .or. cfg%Lcltcalipso .or.  &
        cfg%Lclcalipsotmp    .or. cfg%Lclcalipsoice  .or. cfg%Lclcalipsotmpun .or.                       &
        cfg%Lclcalipsotmpliq .or. cfg%Lcllcalipsoice .or. cfg%Lclmcalipsoice  .or.                       &
        cfg%Lclhcalipsoice   .or. cfg%Lcltcalipsoice .or. cfg%Lcllcalipsoliq  .or.                       &
        cfg%Lclmcalipsoliq   .or. cfg%Lclhcalipsoliq .or. cfg%Lcltcalipsoliq  .or.                       &
        cfg%Lcllcalipsoun    .or. cfg%Lclmcalipsoun  .or. cfg%Lclhcalipsoun   .or. cfg%Lcltcalipsoun) then
       allocate(x%calipso_tau_tot(Npoints,Ncolumns,Nlevels))       
       allocate(x%calipso_temp_tot(Npoints,Nlevels))               
    endif

    ! GROUND LIDAR @ 532NM simulator
    if (cfg%LlidarBetaMol532gr) allocate(x%grLidar532_beta_mol(Npoints,Nlevels))
    if (cfg%Latb532gr)          allocate(x%grLidar532_beta_tot(Npoints,Ncolumns,Nlevels))
    if (cfg%LcfadLidarsr532gr) then 
        allocate(x%grLidar532_srbval(SR_BINS+1)) 
        allocate(x%grLidar532_cfad_sr(Npoints,SR_BINS,Nlvgrid))
    endif
    if (cfg%LclgrLidar532)     allocate(x%grLidar532_lidarcld(Npoints,Nlvgrid)) 
    if (cfg%LclhgrLidar532 .or. cfg%LclmgrLidar532 .or. cfg%LcllgrLidar532 .or. cfg%LcltgrLidar532) then
        allocate(x%grLidar532_cldlayer(Npoints,LIDAR_NCAT))  
    endif
      
    ! ATLID simulator
    if (cfg%LlidarBetaMol355) allocate(x%atlid_beta_mol(Npoints,Nlevels))
    if (cfg%Latb355)          allocate(x%atlid_beta_tot(Npoints,Ncolumns,Nlevels))
    if (cfg%LcfadLidarsr355) then
        allocate(x%atlid_srbval(SR_BINS+1)) 
        allocate(x%atlid_cfad_sr(Npoints,SR_BINS,Nlvgrid))
    endif 
    if (cfg%Lclatlid)     allocate(x%atlid_lidarcld(Npoints,Nlvgrid)) 
    if (cfg%Lclhatlid .or. cfg%Lclmatlid .or. cfg%Lcllatlid .or. cfg%Lcltatlid) then
        allocate(x%atlid_cldlayer(Npoints,LIDAR_NCAT)) 
    endif
      
    ! PARASOL
    if (cfg%Lparasolrefl) then
        allocate(x%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL))
        allocate(x%parasolGrid_refl(Npoints,PARASOL_NREFL))
    endif 

    ! Cloudsat simulator
    if (cfg%Ldbze94)        allocate(x%cloudsat_Ze_tot(Npoints,Ncolumns,Nlevels))
    if (cfg%LcfadDbze94)    allocate(x%cloudsat_cfad_ze(Npoints,cloudsat_DBZE_BINS,Nlvgrid))
    if (cfg%Lptradarflag0 .or. cfg%Lptradarflag1 .or. cfg%Lptradarflag2 .or. cfg%Lptradarflag3 .or. &
        cfg%Lptradarflag4 .or. cfg%Lptradarflag5 .or. cfg%Lptradarflag6 .or. cfg%Lptradarflag7 .or. &
        cfg%Lptradarflag8 .or. cfg%Lptradarflag9) then
       allocate(x%cloudsat_precip_cover(Npoints,cloudsat_DBZE_BINS))
    endif
    if (cfg%Lradarpia) allocate(x%cloudsat_pia(Npoints))

    ! Combined CALIPSO/CLOUDSAT fields
    if (cfg%Lclcalipso2)    allocate(x%lidar_only_freq_cloud(Npoints,Nlvgrid))
    if (cfg%Lcltlidarradar) allocate(x%radar_lidar_tcc(Npoints))
    if (cfg%Lcloudsat_tcc)  allocate(x%cloudsat_tcc(Npoints))
    if (cfg%Lcloudsat_tcc2) allocate(x%cloudsat_tcc2(Npoints))
            
    ! RTTOV
    if (cfg%Ltbrttov) allocate(x%rttov_tbs(Npoints,Nchan))
 
  end subroutine construct_cosp_outputs

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cospIN(y)
    type(cosp_optical_inputs),intent(inout) :: y

    if (allocated(y%tau_067))       	  deallocate(y%tau_067)
    if (allocated(y%emiss_11))      	  deallocate(y%emiss_11)
    if (allocated(y%frac_out))    	  deallocate(y%frac_out)
    if (allocated(y%beta_mol_calipso))    deallocate(y%beta_mol_calipso)
    if (allocated(y%tau_mol_calipso))     deallocate(y%tau_mol_calipso)
    if (allocated(y%betatot_calipso))     deallocate(y%betatot_calipso)
    if (allocated(y%betatot_ice_calipso)) deallocate(y%betatot_ice_calipso)
    if (allocated(y%betatot_liq_calipso)) deallocate(y%betatot_liq_calipso)
    if (allocated(y%tautot_calipso))      deallocate(y%tautot_calipso)
    if (allocated(y%tautot_ice_calipso))  deallocate(y%tautot_ice_calipso)
    if (allocated(y%tautot_liq_calipso))  deallocate(y%tautot_liq_calipso)
    if (allocated(y%tautot_S_liq))	  deallocate(y%tautot_S_liq)
    if (allocated(y%tautot_S_ice))    	  deallocate(y%tautot_S_ice)
    if (allocated(y%z_vol_cloudsat))  	  deallocate(y%z_vol_cloudsat)
    if (allocated(y%kr_vol_cloudsat)) 	  deallocate(y%kr_vol_cloudsat)
    if (allocated(y%g_vol_cloudsat)) 	  deallocate(y%g_vol_cloudsat)
    if (allocated(y%asym))            	  deallocate(y%asym)
    if (allocated(y%ss_alb))         	  deallocate(y%ss_alb)
    if (allocated(y%fracLiq))       	  deallocate(y%fracLiq)
    if (allocated(y%beta_mol_grLidar532)) deallocate(y%beta_mol_grLidar532)
    if (allocated(y%betatot_grLidar532))  deallocate(y%betatot_grLidar532)
    if (allocated(y%tau_mol_grLidar532))  deallocate(y%tau_mol_grLidar532)
    if (allocated(y%tautot_grLidar532))   deallocate(y%tautot_grLidar532)
    if (allocated(y%beta_mol_atlid))      deallocate(y%beta_mol_atlid) 
    if (allocated(y%betatot_atlid))       deallocate(y%betatot_atlid) 
    if (allocated(y%tau_mol_atlid))       deallocate(y%tau_mol_atlid) 
    if (allocated(y%tautot_atlid))        deallocate(y%tautot_atlid)
    if (allocated(y%fracPrecipIce))       deallocate(y%fracPrecipIce)

  end subroutine destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospstateIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cospstateIN(y)
    type(cosp_column_inputs),intent(inout) :: y

    if (allocated(y%sunlit))          deallocate(y%sunlit)
    if (allocated(y%skt))             deallocate(y%skt)
    if (allocated(y%land))            deallocate(y%land)
    if (allocated(y%at))              deallocate(y%at)
    if (allocated(y%pfull))           deallocate(y%pfull)
    if (allocated(y%phalf))           deallocate(y%phalf)
    if (allocated(y%qv))              deallocate(y%qv)
    if (allocated(y%o3))              deallocate(y%o3)
    if (allocated(y%hgt_matrix))      deallocate(y%hgt_matrix)
    if (allocated(y%u_sfc))           deallocate(y%u_sfc)
    if (allocated(y%v_sfc))           deallocate(y%v_sfc)
    if (allocated(y%lat))             deallocate(y%lat)
    if (allocated(y%lon))             deallocate(y%lon)
    if (allocated(y%emis_sfc))        deallocate(y%emis_sfc)
    if (allocated(y%cloudIce))        deallocate(y%cloudIce)
    if (allocated(y%cloudLiq))        deallocate(y%cloudLiq)
    if (allocated(y%seaice))          deallocate(y%seaice)
    if (allocated(y%fl_rain))         deallocate(y%fl_rain)
    if (allocated(y%fl_snow))         deallocate(y%fl_snow)
    if (allocated(y%tca))             deallocate(y%tca)
    if (allocated(y%hgt_matrix_half)) deallocate(y%hgt_matrix_half)    
    if (allocated(y%surfelev))        deallocate(y%surfelev)
    
  end subroutine destroy_cospstateIN
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_outputs
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cosp_outputs(y)
     type(cosp_outputs),intent(inout) :: y

     ! Deallocate and nullify
     if (associated(y%calipso_beta_mol))          then
        deallocate(y%calipso_beta_mol)
        nullify(y%calipso_beta_mol)
     endif
     if (associated(y%calipso_temp_tot))          then
        deallocate(y%calipso_temp_tot)
        nullify(y%calipso_temp_tot)     
     endif
     if (associated(y%calipso_betaperp_tot))      then
        deallocate(y%calipso_betaperp_tot)
        nullify(y%calipso_betaperp_tot)     
     endif
     if (associated(y%calipso_beta_tot))          then
        deallocate(y%calipso_beta_tot)    
        nullify(y%calipso_beta_tot)     
     endif
     if (associated(y%calipso_tau_tot))           then
        deallocate(y%calipso_tau_tot) 
        nullify(y%calipso_tau_tot)     
     endif
     if (associated(y%calipso_lidarcldphase))     then
        deallocate(y%calipso_lidarcldphase)
        nullify(y%calipso_lidarcldphase)     
     endif
     if (associated(y%calipso_lidarcldtype))     then
        deallocate(y%calipso_lidarcldtype) 
        nullify(y%calipso_lidarcldtype) 
     endif
     if (associated(y%calipso_cldlayerphase))     then
        deallocate(y%calipso_cldlayerphase)
        nullify(y%calipso_cldlayerphase)     
     endif
     if (associated(y%calipso_lidarcldtmp))       then
        deallocate(y%calipso_lidarcldtmp)
        nullify(y%calipso_lidarcldtmp)     
     endif
     if (associated(y%calipso_cldlayer))          then
        deallocate(y%calipso_cldlayer)
        nullify(y%calipso_cldlayer)     
     endif
     if (associated(y%calipso_cldtype))          then
        deallocate(y%calipso_cldtype)
        nullify(y%calipso_cldtype)  
     endif  
     if (associated(y%calipso_cldtypetemp))      then
        deallocate(y%calipso_cldtypetemp) 
        nullify(y%calipso_cldtypetemp) 
     endif 
     if (associated(y%calipso_cldtypemeanz))     then
        deallocate(y%calipso_cldtypemeanz)
        nullify(y%calipso_cldtypemeanz)  
     endif  
     if (associated(y%calipso_cldtypemeanzse))   then 
        deallocate(y%calipso_cldtypemeanzse) 
        nullify(y%calipso_cldtypemeanzse)  
     endif 
     if (associated(y%calipso_cldthinemis))      then
        deallocate(y%calipso_cldthinemis)
        nullify(y%calipso_cldthinemis)
     endif 
     if (associated(y%calipso_lidarcld))         then
        deallocate(y%calipso_lidarcld)
        nullify(y%calipso_lidarcld)     
     endif
     if (associated(y%calipso_srbval))            then
        deallocate(y%calipso_srbval)
        nullify(y%calipso_srbval)     
     endif
     if (associated(y%calipso_cfad_sr))          then
        deallocate(y%calipso_cfad_sr)
        nullify(y%calipso_cfad_sr)     
     endif
     if (associated(y%grLidar532_beta_mol))     then
        deallocate(y%grLidar532_beta_mol)
        nullify(y%grLidar532_beta_mol)  
     endif
     if (associated(y%grLidar532_beta_tot))     then 
        deallocate(y%grLidar532_beta_tot) 
        nullify(y%grLidar532_beta_tot)
     endif 
     if (associated(y%grLidar532_cldlayer))     then 
        deallocate(y%grLidar532_cldlayer) 
        nullify(y%grLidar532_cldlayer) 
     endif 
     if (associated(y%grLidar532_lidarcld))     then 
        deallocate(y%grLidar532_lidarcld) 
        nullify(y%grLidar532_lidarcld) 
     endif 
     if (associated(y%grLidar532_cfad_sr))      then 
        deallocate(y%grLidar532_cfad_sr)  
        nullify(y%grLidar532_cfad_sr) 
     endif 
     if (associated(y%grLidar532_srbval))       then
        deallocate(y%grLidar532_srbval) 
        nullify(y%grLidar532_srbval) 
     endif  
     if (associated(y%atlid_beta_mol))           then
        deallocate(y%atlid_beta_mol) 
        nullify(y%atlid_beta_mol) 
     endif 
     if (associated(y%atlid_beta_tot))           then
        deallocate(y%atlid_beta_tot) 
        nullify(y%atlid_beta_tot) 
     endif
     if (associated(y%atlid_cldlayer))           then 
        deallocate(y%atlid_cldlayer)  
        nullify(y%atlid_cldlayer) 
     endif 
     if (associated(y%atlid_lidarcld))           then 
        deallocate(y%atlid_lidarcld)  
        nullify(y%atlid_lidarcld) 
     endif  
     if (associated(y%atlid_cfad_sr))            then
        deallocate(y%atlid_cfad_sr) 
        nullify(y%atlid_cfad_sr) 
     endif
     if (associated(y%atlid_srbval))             then 
        deallocate(y%atlid_srbval) 
        nullify(y%atlid_srbval) 
     endif
     if (associated(y%parasolPix_refl))           then
        deallocate(y%parasolPix_refl)
        nullify(y%parasolPix_refl)     
     endif
     if (associated(y%parasolGrid_refl))          then
        deallocate(y%parasolGrid_refl) 
        nullify(y%parasolGrid_refl)     
     endif
     if (associated(y%cloudsat_Ze_tot))           then
        deallocate(y%cloudsat_Ze_tot) 
        nullify(y%cloudsat_Ze_tot)  
     endif
     if (associated(y%cloudsat_cfad_ze))          then
        deallocate(y%cloudsat_cfad_ze)
        nullify(y%cloudsat_cfad_ze)     
     endif
     if (associated(y%cloudsat_precip_cover))     then
        deallocate(y%cloudsat_precip_cover)
        nullify(y%cloudsat_precip_cover)
     endif
     if (associated(y%cloudsat_pia))              then
        deallocate(y%cloudsat_pia)
        nullify(y%cloudsat_pia)
     endif
     if (associated(y%cloudsat_tcc))           then
        deallocate(y%cloudsat_tcc) 
        nullify(y%cloudsat_tcc)  
     endif
     if (associated(y%cloudsat_tcc2))           then
        deallocate(y%cloudsat_tcc2) 
        nullify(y%cloudsat_tcc2)  
     endif
     if (associated(y%radar_lidar_tcc))           then
        deallocate(y%radar_lidar_tcc) 
        nullify(y%radar_lidar_tcc)  
     endif
     if (associated(y%cloudsat_tcc))           then
        deallocate(y%cloudsat_tcc) 
        nullify(y%cloudsat_tcc)  
     endif
     if (associated(y%cloudsat_tcc2))           then
        deallocate(y%cloudsat_tcc2) 
        nullify(y%cloudsat_tcc2)  
     endif
     if (associated(y%lidar_only_freq_cloud))     then
        deallocate(y%lidar_only_freq_cloud)
        nullify(y%lidar_only_freq_cloud)     
     endif
     if (associated(y%isccp_totalcldarea))        then
        deallocate(y%isccp_totalcldarea) 
        nullify(y%isccp_totalcldarea)  
     endif
     if (associated(y%isccp_meantb))              then
        deallocate(y%isccp_meantb) 
        nullify(y%isccp_meantb)     
     endif
     if (associated(y%isccp_meantbclr))           then
        deallocate(y%isccp_meantbclr)
        nullify(y%isccp_meantbclr)  
     endif
     if (associated(y%isccp_meanptop))            then
        deallocate(y%isccp_meanptop)
        nullify(y%isccp_meanptop)     
     endif
     if (associated(y%isccp_meantaucld))          then
        deallocate(y%isccp_meantaucld) 
        nullify(y%isccp_meantaucld)       
     endif
     if (associated(y%isccp_meanalbedocld))       then
        deallocate(y%isccp_meanalbedocld)
        nullify(y%isccp_meanalbedocld)     
     endif
     if (associated(y%isccp_boxtau))              then
        deallocate(y%isccp_boxtau)
        nullify(y%isccp_boxtau)       
     endif
     if (associated(y%isccp_boxptop))             then
        deallocate(y%isccp_boxptop)
        nullify(y%isccp_boxptop)     
     endif
     if (associated(y%isccp_fq))                  then
        deallocate(y%isccp_fq)
        nullify(y%isccp_fq)       
     endif
     if (associated(y%misr_fq))                   then
        deallocate(y%misr_fq) 
        nullify(y%misr_fq)     
     endif
     if (associated(y%misr_dist_model_layertops)) then
        deallocate(y%misr_dist_model_layertops)
        nullify(y%misr_dist_model_layertops)       
     endif
     if (associated(y%misr_meanztop))             then
        deallocate(y%misr_meanztop)
        nullify(y%misr_meanztop)     
     endif
     if (associated(y%misr_cldarea))              then
        deallocate(y%misr_cldarea)
        nullify(y%misr_cldarea)      
     endif
     if (associated(y%rttov_tbs))                 then
        deallocate(y%rttov_tbs)
        nullify(y%rttov_tbs)     
     endif
     if (associated(y%modis_Cloud_Fraction_Total_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Total_Mean)       
        nullify(y%modis_Cloud_Fraction_Total_Mean)       
     endif
     if (associated(y%modis_Cloud_Fraction_Ice_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Ice_Mean)     
        nullify(y%modis_Cloud_Fraction_Ice_Mean)     
     endif
     if (associated(y%modis_Cloud_Fraction_Water_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Water_Mean)           
        nullify(y%modis_Cloud_Fraction_Water_Mean)           
     endif
     if (associated(y%modis_Cloud_Fraction_High_Mean))                       then
        deallocate(y%modis_Cloud_Fraction_High_Mean)     
        nullify(y%modis_Cloud_Fraction_High_Mean)     
     endif
     if (associated(y%modis_Cloud_Fraction_Mid_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Mid_Mean)       
        nullify(y%modis_Cloud_Fraction_Mid_Mean)       
     endif
     if (associated(y%modis_Cloud_Fraction_Low_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Low_Mean)     
        nullify(y%modis_Cloud_Fraction_Low_Mean)     
     endif
     if (associated(y%modis_Optical_Thickness_Total_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Total_Mean)  
        nullify(y%modis_Optical_Thickness_Total_Mean)  
     endif
     if (associated(y%modis_Optical_Thickness_Water_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Water_Mean)     
        nullify(y%modis_Optical_Thickness_Water_Mean)     
     endif
     if (associated(y%modis_Optical_Thickness_Ice_Mean))                     then
        deallocate(y%modis_Optical_Thickness_Ice_Mean)       
        nullify(y%modis_Optical_Thickness_Ice_Mean)       
     endif
     if (associated(y%modis_Optical_Thickness_Total_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Total_LogMean)    
        nullify(y%modis_Optical_Thickness_Total_LogMean)    
     endif
     if (associated(y%modis_Optical_Thickness_Water_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Water_LogMean)     
        nullify(y%modis_Optical_Thickness_Water_LogMean)     
     endif
     if (associated(y%modis_Optical_Thickness_Ice_LogMean))                  then
        deallocate(y%modis_Optical_Thickness_Ice_LogMean)     
        nullify(y%modis_Optical_Thickness_Ice_LogMean)     
     endif
     if (associated(y%modis_Cloud_Particle_Size_Water_Mean))                 then
        deallocate(y%modis_Cloud_Particle_Size_Water_Mean)       
        nullify(y%modis_Cloud_Particle_Size_Water_Mean)       
     endif
     if (associated(y%modis_Cloud_Particle_Size_Ice_Mean))                   then
        deallocate(y%modis_Cloud_Particle_Size_Ice_Mean)     
        nullify(y%modis_Cloud_Particle_Size_Ice_Mean)     
     endif
     if (associated(y%modis_Cloud_Top_Pressure_Total_Mean))                  then
        deallocate(y%modis_Cloud_Top_Pressure_Total_Mean)           
        nullify(y%modis_Cloud_Top_Pressure_Total_Mean)           
     endif
     if (associated(y%modis_Liquid_Water_Path_Mean))                         then
        deallocate(y%modis_Liquid_Water_Path_Mean)     
        nullify(y%modis_Liquid_Water_Path_Mean)     
     endif
     if (associated(y%modis_Ice_Water_Path_Mean))                            then
        deallocate(y%modis_Ice_Water_Path_Mean)       
        nullify(y%modis_Ice_Water_Path_Mean)       
     endif
     if (associated(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure))        then
        deallocate(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
        nullify(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffLIQ))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffLIQ)
        nullify(y%modis_Optical_thickness_vs_ReffLIQ)
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffICE))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffICE)
        nullify(y%modis_Optical_thickness_vs_ReffICE)
     endif
        
   end subroutine destroy_cosp_outputs

END MODULE LMDZ_COSP_CONSTRUCT_DESTROY_MOD

