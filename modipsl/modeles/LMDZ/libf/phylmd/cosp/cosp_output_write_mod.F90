!!!! Abderrahmane Idelkadi aout 2013 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module pour definir (au 1er appel) et ecrire les variables dans les fichiers de sortie cosp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   MODULE cosp_output_write_mod
  
   USE cosp_output_mod
  
   IMPLICIT NONE

   INTEGER, SAVE  :: itau_iocosp
!$OMP THREADPRIVATE(itau_iocosp)
   INTEGER, save        :: Nlevout, Ncolout
!$OMP THREADPRIVATE(Nlevout, Ncolout)

!  INTERFACE histwrite_cosp
!    MODULE PROCEDURE histwrite2d_cosp,histwrite3d_cosp
!  END INTERFACE

   CONTAINS

  SUBROUTINE cosp_output_write(Nlevlmdz, Npoints, Ncolumns, itap, dtime, freq_COSP, missing_cosp, &
                               cfg, gbx, vgrid, sglidar, sgradar, stlidar, stradar, &
                               isccp, misr, modis)

    USE ioipsl
    USE time_phylmdz_mod, ONLY: itau_phy, start_time, day_step_phy
    USE print_control_mod, ONLY: lunout,prt_level

#ifdef CPP_XIOS
    USE wxios, only: wxios_closedef
    USE xios, only: xios_update_calendar, xios_field_is_active
#endif
  IMPLICIT NONE  
!!! Variables d'entree
  integer               :: itap, Nlevlmdz, Ncolumns, Npoints
  real                  :: freq_COSP, dtime, missing_val, missing_cosp
  type(cosp_config)     :: cfg     ! Control outputs
  type(cosp_gridbox)    :: gbx     ! Gridbox information. Input for COSP
  type(cosp_sglidar)    :: sglidar ! Output from lidar simulator
  type(cosp_sgradar)    :: sgradar ! Output from radar simulator
  type(cosp_isccp)      :: isccp   ! Output from ISCCP simulator
  type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator
  type(cosp_radarstats) :: stradar
  type(cosp_misr)       :: misr    ! Output from MISR
  type(cosp_modis)      :: modis   ! Outputs from Modis
  type(cosp_vgrid)      :: vgrid   ! Information on vertical grid of stats

!!! Variables locales
  integer               :: icl,k,ip
  logical               :: ok_sync
  integer               :: itau_wcosp, iff
  real, dimension(Npoints,PARASOL_NREFL) :: parasolcrefl, Ncref

! Variables locals intermidiaires pour inverser les axes des champs 4D
! Compatibilite avec sorties CMIP
  real, dimension(Npoints,Nlevout,SR_BINS) :: tmp_fi4da_cfadL
  real, dimension(Npoints,Nlevout,DBZE_BINS) :: tmp_fi4da_cfadR
  real, dimension(Npoints,MISR_N_CTH,7) :: tmp_fi4da_misr

#ifdef CPP_XIOS
  missing_val=missing_cosp
#else
  missing_val=0.
#endif

  Nlevout = vgrid%Nlvgrid
  Ncolout = Ncolumns

! A refaire
       itau_wcosp = itau_phy + itap + start_time * day_step_phy
        if (prt_level >= 10) then
             WRITE(lunout,*)'itau_wcosp, itap, start_time, day_step_phy =', & 
                             itau_wcosp, itap, start_time, day_step_phy
        endif

! On le donne a  cosp_output_write_mod pour que les histwrite y aient acces:
       CALL set_itau_iocosp(itau_wcosp)
        if (prt_level >= 10) then
              WRITE(lunout,*)'itau_iocosp =',itau_iocosp
        endif

    ok_sync = .TRUE.
    
!DO iinit=1, iinitend
! AI sept 2014 cette boucle supprimee
! On n'ecrit pas quand itap=1 (cosp)

!   if (prt_level >= 10) then
!         WRITE(lunout,*)'DO iinit=1, iinitend ',iinitend
!   endif

!!#ifdef CPP_XIOS
! !$OMP MASTER
!IF (cosp_varsdefined) THEN
!   if (prt_level >= 10) then
!         WRITE(lunout,*)'Apell xios_update_calendar cosp_varsdefined iinitend ', &
!                         cosp_varsdefined,iinitend
!   endif 
!    CALL xios_update_calendar(itau_wcosp)
!ENDIF
!  !$OMP END MASTER
!  !$OMP BARRIER
!!#endif

!!!! Sorties Calipso
 if (cfg%Llidar_sim) then
!!! AI 02 2018 
! Traitement missing_val
   where(stlidar%lidarcld == R_UNDEF) stlidar%lidarcld = missing_val
   where(stlidar%proftemp == R_UNDEF) stlidar%proftemp = missing_val   !TIBO  
   where(stlidar%profSR == R_UNDEF) stlidar%profSR = missing_val       !TIBO2
   where(sglidar%beta_mol == R_UNDEF) sglidar%beta_mol = missing_val  
   where(sglidar%beta_tot == R_UNDEF) sglidar%beta_tot = missing_val 
   where(stlidar%cldlayer == R_UNDEF) stlidar%cldlayer = missing_val
   where(stlidar%cldtype == R_UNDEF) stlidar%cldtype = missing_val     !OPAQ
   where(stlidar%cfad_sr == R_UNDEF) stlidar%cfad_sr = missing_val
! AI 11 / 2015
   where(stlidar%parasolrefl == R_UNDEF) stlidar%parasolrefl = missing_val
   where(stlidar%lidarcldtmp == R_UNDEF) stlidar%lidarcldtmp = missing_val
   where(stlidar%cldlayerphase == R_UNDEF) stlidar%cldlayerphase = missing_val
   where(stlidar%lidarcldphase == R_UNDEF) stlidar%lidarcldphase = missing_val
   where(stlidar%lidarcldtype == R_UNDEF) stlidar%lidarcldtype = missing_val   !OPAQ
   where(stlidar%lidarcldtmp == R_UNDEF) stlidar%lidarcldtmp = missing_val
 
!   print*,'Appel histwrite2d_cosp'
   if (cfg%Lcllcalipso) CALL histwrite2d_cosp(o_cllcalipso,stlidar%cldlayer(:,1))
   if (cfg%Lclhcalipso) CALL histwrite2d_cosp(o_clhcalipso,stlidar%cldlayer(:,3))
   if (cfg%Lclmcalipso) CALL histwrite2d_cosp(o_clmcalipso,stlidar%cldlayer(:,2)) 
   if (cfg%Lcltcalipso) CALL histwrite2d_cosp(o_cltcalipso,stlidar%cldlayer(:,4))
   if (cfg%Lclcalipso) CALL histwrite3d_cosp(o_clcalipso,stlidar%lidarcld,nvert)
   if (cfg%Lclcalipsotmp) CALL histwrite3d_cosp(o_clcalipsotmp,stlidar%lidarcldtmp(:,:,1),nverttemp)

   if (cfg%Lcllcalipsoice) CALL histwrite2d_cosp(o_cllcalipsoice,stlidar%cldlayerphase(:,1,1))
   if (cfg%Lclhcalipsoice) CALL histwrite2d_cosp(o_clhcalipsoice,stlidar%cldlayerphase(:,3,1))
   if (cfg%Lclmcalipsoice) CALL histwrite2d_cosp(o_clmcalipsoice,stlidar%cldlayerphase(:,2,1))
   if (cfg%Lcltcalipsoice) CALL histwrite2d_cosp(o_cltcalipsoice,stlidar%cldlayerphase(:,4,1))
   if (cfg%Lclcalipsoice) CALL histwrite3d_cosp(o_clcalipsoice,stlidar%lidarcldphase(:,:,1),nvert)
   if (cfg%Lclcalipsotmpice) CALL histwrite3d_cosp(o_clcalipsotmpice,stlidar%lidarcldtmp(:,:,2),nverttemp)

   if (cfg%Lcllcalipsoliq) CALL histwrite2d_cosp(o_cllcalipsoliq,stlidar%cldlayerphase(:,1,2))
   if (cfg%Lclhcalipsoliq) CALL histwrite2d_cosp(o_clhcalipsoliq,stlidar%cldlayerphase(:,3,2))
   if (cfg%Lclmcalipsoliq) CALL histwrite2d_cosp(o_clmcalipsoliq,stlidar%cldlayerphase(:,2,2))
   if (cfg%Lcltcalipsoliq) CALL histwrite2d_cosp(o_cltcalipsoliq,stlidar%cldlayerphase(:,4,2))
   if (cfg%Lclcalipsoliq) CALL histwrite3d_cosp(o_clcalipsoliq,stlidar%lidarcldphase(:,:,2),nvert)
   if (cfg%Lclcalipsotmpliq) CALL histwrite3d_cosp(o_clcalipsotmpliq,stlidar%lidarcldtmp(:,:,3),nverttemp)

   if (cfg%Lcllcalipsoun) CALL histwrite2d_cosp(o_cllcalipsoun,stlidar%cldlayerphase(:,1,3))
   if (cfg%Lclhcalipsoun) CALL histwrite2d_cosp(o_clhcalipsoun,stlidar%cldlayerphase(:,3,3))
   if (cfg%Lclmcalipsoun) CALL histwrite2d_cosp(o_clmcalipsoun,stlidar%cldlayerphase(:,2,3))
   if (cfg%Lcltcalipsoun) CALL histwrite2d_cosp(o_cltcalipsoun,stlidar%cldlayerphase(:,4,3))
   if (cfg%Lclcalipsoun) CALL histwrite3d_cosp(o_clcalipsoun,stlidar%lidarcldphase(:,:,3),nvert)
   if (cfg%Lclcalipsotmpun) CALL histwrite3d_cosp(o_clcalipsotmpun,stlidar%lidarcldtmp(:,:,4),nverttemp)

   if (cfg%Lclopaquecalipso) CALL histwrite2d_cosp(o_clopaquecalipso,stlidar%cldtype(:,1))               !OPAQ
   if (cfg%Lclthincalipso) CALL histwrite2d_cosp(o_clthincalipso,stlidar%cldtype(:,2))                 !OPAQ
   if (cfg%Lclzopaquecalipso) CALL histwrite2d_cosp(o_clzopaquecalipso,stlidar%cldtype(:,3))              !OPAQ

   if (cfg%Lclcalipsoopaque) CALL histwrite3d_cosp(o_clcalipsoopaque,stlidar%lidarcldtype(:,:,1),nvert)  !OPAQ
   if (cfg%Lclcalipsothin) CALL histwrite3d_cosp(o_clcalipsothin,stlidar%lidarcldtype(:,:,2),nvert)    !OPAQ
   if (cfg%Lclcalipsozopaque) CALL histwrite3d_cosp(o_clcalipsozopaque,stlidar%lidarcldtype(:,:,3),nvert) !OPAQ
   if (cfg%Lclcalipsoopacity) CALL histwrite3d_cosp(o_clcalipsoopacity,stlidar%lidarcldtype(:,:,4),nvert) !OPAQ

   if (cfg%Lproftemp) CALL histwrite3d_cosp(o_proftemp,stlidar%proftemp,nvert)                    !TIBO

#ifdef CPP_XIOS
   do icl=1,SR_BINS
      tmp_fi4da_cfadL(:,:,icl)=stlidar%cfad_sr(:,icl,:)
   enddo
!   if (cfg%LcfadLidarsr532) CALL histwrite4d_cosp(o_cfad_lidarsr532,stlidar%cfad_sr)
   if (cfg%LcfadLidarsr532) CALL histwrite4d_cosp(o_cfad_lidarsr532,tmp_fi4da_cfadL)
   if (cfg%LprofSR) CALL histwrite4d_cosp(o_profSR,stlidar%profSR)                              !TIBO
#else
   if (cfg%LcfadLidarsr532) then
     do icl=1,SR_BINS
        CALL histwrite3d_cosp(o_cfad_lidarsr532,stlidar%cfad_sr(:,icl,:),nvert,icl)
     enddo
   endif
   if (cfg%LprofSR) then
     do icl=1,Ncolumns                                                              !TIBO
        CALL histwrite3d_cosp(o_profSR,stlidar%profSR(:,icl,:),nvert,icl)           !TIBO
     enddo                                                                          !TIBO
   endif
#endif
   if (cfg%LparasolRefl) CALL histwrite3d_cosp(o_parasol_refl,stlidar%parasolrefl,nvertp)

  if (cfg%LparasolRefl) then 
    do k=1,PARASOL_NREFL
     do ip=1, Npoints
      if (stlidar%cldlayer(ip,4).gt.0.01.and.stlidar%parasolrefl(ip,k).ne.missing_val) then
        parasolcrefl(ip,k)=(stlidar%parasolrefl(ip,k)-0.03*(1.-stlidar%cldlayer(ip,4)))/ &
                             stlidar%cldlayer(ip,4)
         Ncref(ip,k) = 1.
      else
         parasolcrefl(ip,k)=missing_val
         Ncref(ip,k) = 0.
      endif
     enddo
    enddo
    CALL histwrite3d_cosp(o_Ncrefl,Ncref,nvertp)
    CALL histwrite3d_cosp(o_parasol_crefl,parasolcrefl,nvertp)
  endif

#ifdef CPP_XIOS
   if (cfg%Latb532) CALL histwrite4d_cosp(o_atb532,sglidar%beta_tot)
#else
   if (cfg%Latb532) then  
     do icl=1,Ncolumns 
        CALL histwrite3d_cosp(o_atb532,sglidar%beta_tot(:,icl,:),nvertmcosp,icl)
     enddo
   endif 
#endif

   if (cfg%LlidarBetaMol532) CALL histwrite3d_cosp(o_beta_mol532,sglidar%beta_mol,nvertmcosp) 

 endif !Lidar

!!! Sorties Cloudsat
 if (cfg%Lradar_sim) then

   where(stradar%cfad_ze == R_UNDEF) stradar%cfad_ze = missing_val
#ifdef CPP_XIOS
   do icl=1,DBZE_BINS
     tmp_fi4da_cfadR(:,:,icl)=stradar%cfad_ze(:,icl,:)
   enddo
   if (cfg%Ldbze94) CALL histwrite4d_cosp(o_dbze94,sgradar%Ze_tot)
!   if (cfg%LcfadDbze94) CALL histwrite4d_cosp(o_cfadDbze94,stradar%cfad_ze)
   if (cfg%LcfadDbze94) CALL histwrite4d_cosp(o_cfadDbze94,tmp_fi4da_cfadR)
#else
   if (cfg%Ldbze94) then
    do icl=1,Ncolumns
       CALL histwrite3d_cosp(o_dbze94,sgradar%Ze_tot(:,icl,:),nvert,icl)
    enddo
   endif
   if (cfg%LcfadDbze94) then
    do icl=1,DBZE_BINS
    CALL histwrite3d_cosp(o_cfadDbze94,stradar%cfad_ze(:,icl,:),nvert,icl)
    enddo
   endif
#endif
 endif
! endif pour radar

!!! Sorties combinees Cloudsat et Calipso
 if (cfg%Llidar_sim .and. cfg%Lradar_sim) then
   where(stradar%lidar_only_freq_cloud == R_UNDEF) &
                           stradar%lidar_only_freq_cloud = missing_val
   if (cfg%Lclcalipso) CALL histwrite3d_cosp(o_clcalipso2,stradar%lidar_only_freq_cloud,nvert)
   where(stradar%radar_lidar_tcc == R_UNDEF) &
                           stradar%radar_lidar_tcc = missing_val
   if (cfg%Lcltlidarradar) CALL histwrite2d_cosp(o_cltlidarradar,stradar%radar_lidar_tcc)
 endif

!!! Sorties Isccp
 if (cfg%Lisccp_sim) then
  where(isccp%totalcldarea == R_UNDEF) isccp%totalcldarea = missing_val
  where(isccp%meanptop == R_UNDEF) isccp%meanptop = missing_val
  where(isccp%meantaucld == R_UNDEF) isccp%meantaucld = missing_val
  where(isccp%meanalbedocld == R_UNDEF) isccp%meanalbedocld = missing_val
  where(isccp%meantb == R_UNDEF) isccp%meantb = missing_val
  where(isccp%meantbclr == R_UNDEF) isccp%meantbclr = missing_val
  where(isccp%fq_isccp == R_UNDEF) isccp%fq_isccp = missing_val
  where(isccp%boxtau == R_UNDEF) isccp%boxtau = missing_val
  where(isccp%boxptop == R_UNDEF) isccp%boxptop = missing_val 

   CALL histwrite2d_cosp(o_sunlit,gbx%sunlit)
#ifdef CPP_XIOS
  if (cfg%Lclisccp) CALL histwrite4d_cosp(o_clisccp2,isccp%fq_isccp)
#else
   if (cfg%Lclisccp) then
     do icl=1,7
       CALL histwrite3d_cosp(o_clisccp2,isccp%fq_isccp(:,icl,:),nvertisccp,icl) 
     enddo
   endif
#endif

   if (cfg%Lboxtauisccp) CALL histwrite3d_cosp(o_boxtauisccp,isccp%boxtau,nvertcol)
   if (cfg%Lboxptopisccp) CALL histwrite3d_cosp(o_boxptopisccp,isccp%boxptop,nvertcol) 
   if (cfg%Lcltisccp) CALL histwrite2d_cosp(o_tclisccp,isccp%totalcldarea) 
   if (cfg%Lpctisccp) CALL histwrite2d_cosp(o_ctpisccp,isccp%meanptop) 
   if (cfg%Ltauisccp) CALL histwrite2d_cosp(o_tauisccp,isccp%meantaucld) 
   if (cfg%Lalbisccp) CALL histwrite2d_cosp(o_albisccp,isccp%meanalbedocld) 
   if (cfg%Lmeantbisccp) CALL histwrite2d_cosp(o_meantbisccp,isccp%meantb) 
   if (cfg%Lmeantbclrisccp) CALL histwrite2d_cosp(o_meantbclrisccp,isccp%meantbclr)
 endif ! Isccp

!!! MISR simulator
 if (cfg%Lmisr_sim) then
   where(misr%fq_MISR == R_UNDEF) misr%fq_MISR = missing_val

#ifdef CPP_XIOS
   do icl=1,MISR_N_CTH
      tmp_fi4da_misr(:,icl,:)=misr%fq_MISR(:,:,icl)
   enddo
!   if (cfg%LclMISR) CALL histwrite4d_cosp(o_clMISR,misr%fq_MISR)
   if (cfg%LclMISR) CALL histwrite4d_cosp(o_clMISR,tmp_fi4da_misr)
#else
   if (cfg%LclMISR) then
    do icl=1,7 
      CALL histwrite3d_cosp(o_clMISR,misr%fq_MISR(:,icl,:),nvertmisr,icl)
    enddo
   endif
#endif
 endif
! endif pour Misr

!!! Modis simulator
 if (cfg%Lmodis_sim) then
  where(modis%Cloud_Fraction_Low_Mean == R_UNDEF) &
        modis%Cloud_Fraction_Low_Mean = missing_val
  where(modis%Cloud_Fraction_High_Mean == R_UNDEF) &
        modis%Cloud_Fraction_High_Mean = missing_val
  where(modis%Cloud_Fraction_Mid_Mean == R_UNDEF) &
        modis%Cloud_Fraction_Mid_Mean = missing_val
  where(modis%Cloud_Fraction_Total_Mean == R_UNDEF) &
        modis%Cloud_Fraction_Total_Mean = missing_val
  where(modis%Cloud_Fraction_Water_Mean == R_UNDEF) &
        modis%Cloud_Fraction_Water_Mean = missing_val
  where(modis%Cloud_Fraction_Ice_Mean == R_UNDEF) &
        modis%Cloud_Fraction_Ice_Mean = missing_val
  where(modis%Optical_Thickness_Total_Mean == R_UNDEF) &
        modis%Optical_Thickness_Total_Mean = missing_val
  where(modis%Optical_Thickness_Water_Mean == R_UNDEF) &
        modis%Optical_Thickness_Water_Mean = missing_val
  where(modis%Optical_Thickness_Ice_Mean == R_UNDEF) &
        modis%Optical_Thickness_Ice_Mean = missing_val
  where(modis%Cloud_Particle_Size_Water_Mean == R_UNDEF) &
        modis%Cloud_Particle_Size_Water_Mean = missing_val
  where(modis%Cloud_Particle_Size_Ice_Mean == R_UNDEF) &
        modis%Cloud_Particle_Size_Ice_Mean = missing_val
  where(modis%Cloud_Top_Pressure_Total_Mean == R_UNDEF) &
        modis%Cloud_Top_Pressure_Total_Mean = missing_val
  where(modis%Liquid_Water_Path_Mean == R_UNDEF) &
        modis%Liquid_Water_Path_Mean = missing_val 
  where(modis%Ice_Water_Path_Mean == R_UNDEF) &
        modis%Ice_Water_Path_Mean = missing_val

  where(modis%Optical_Thickness_Total_LogMean == R_UNDEF) &
          modis%Optical_Thickness_Total_LogMean = missing_val
           
  where(modis%Optical_Thickness_Water_LogMean == R_UNDEF) &
          modis%Optical_Thickness_Water_LogMean = missing_val

  where(modis%Optical_Thickness_Ice_LogMean == R_UNDEF) &
          modis%Optical_Thickness_Ice_LogMean = missing_val
    
  if (cfg%Lcllmodis) CALL histwrite2d_cosp(o_cllmodis,modis%Cloud_Fraction_Low_Mean)
  if (cfg%Lclhmodis) CALL histwrite2d_cosp(o_clhmodis,modis%Cloud_Fraction_High_Mean)
  if (cfg%Lclmmodis) CALL histwrite2d_cosp(o_clmmodis,modis%Cloud_Fraction_Mid_Mean)
  if (cfg%Lcltmodis) CALL histwrite2d_cosp(o_cltmodis,modis%Cloud_Fraction_Total_Mean)
  if (cfg%Lclwmodis) CALL histwrite2d_cosp(o_clwmodis,modis%Cloud_Fraction_Water_Mean)
  if (cfg%Lclimodis) CALL histwrite2d_cosp(o_climodis,modis%Cloud_Fraction_Ice_Mean)
  if (cfg%Ltautmodis)  CALL histwrite2d_cosp(o_tautmodis,modis%Optical_Thickness_Total_Mean)
  if (cfg%Ltauwmodis) CALL histwrite2d_cosp(o_tauwmodis,modis%Optical_Thickness_Water_Mean)
  if (cfg%Ltauimodis) CALL histwrite2d_cosp(o_tauimodis,modis%Optical_Thickness_Ice_Mean)
  if (cfg%Ltautlogmodis) CALL histwrite2d_cosp(o_tautlogmodis,modis%Optical_Thickness_Total_LogMean)  
  if (cfg%Ltauwlogmodis) CALL histwrite2d_cosp(o_tauwlogmodis,modis%Optical_Thickness_Water_LogMean)
  if (cfg%Ltauilogmodis) CALL histwrite2d_cosp(o_tauilogmodis,modis%Optical_Thickness_Ice_LogMean)
  if (cfg%Lreffclwmodis) CALL histwrite2d_cosp(o_reffclwmodis,modis%Cloud_Particle_Size_Water_Mean)
  if (cfg%Lreffclimodis) CALL histwrite2d_cosp(o_reffclimodis,modis%Cloud_Particle_Size_Ice_Mean)
  if (cfg%Lpctmodis) CALL histwrite2d_cosp(o_pctmodis,modis%Cloud_Top_Pressure_Total_Mean)
  if (cfg%Llwpmodis) CALL histwrite2d_cosp(o_lwpmodis,modis%Liquid_Water_Path_Mean)
  if (cfg%Liwpmodis) CALL histwrite2d_cosp(o_iwpmodis,modis%Ice_Water_Path_Mean)

    where(modis%Optical_Thickness_vs_Cloud_Top_Pressure == R_UNDEF) &
          modis%Optical_Thickness_vs_Cloud_Top_Pressure = missing_val

#ifdef CPP_XIOS
   if (cfg%Lclmodis) CALL histwrite4d_cosp(o_clmodis,modis%Optical_Thickness_vs_Cloud_Top_Pressure)
#else
  if (cfg%Lclmodis) then
   do icl=1,7
   CALL histwrite3d_cosp(o_clmodis, &
     modis%Optical_Thickness_vs_Cloud_Top_Pressure(:,icl,:),nvertisccp,icl)           
   enddo
  endif 
#endif

    where(modis%Optical_Thickness_vs_ReffIce == R_UNDEF) &
          modis%Optical_Thickness_vs_ReffIce = missing_val

    where(modis%Optical_Thickness_vs_ReffLiq == R_UNDEF) &
          modis%Optical_Thickness_vs_ReffLiq = missing_val

#ifdef CPP_XIOS
  if (cfg%Lcrimodis) CALL histwrite4d_cosp(o_crimodis,modis%Optical_Thickness_vs_ReffIce)
  if (cfg%Lcrlmodis) CALL histwrite4d_cosp(o_crlmodis,modis%Optical_Thickness_vs_ReffLiq)
#else
  if (cfg%Lcrimodis) then
    do icl=1,7
     CALL histwrite3d_cosp(o_crimodis, &
          modis%Optical_Thickness_vs_ReffIce(:,icl,:),nvertReffIce,icl)
    enddo
  endif
  if (cfg%Lcrlmodis) then
    do icl=1,7 
     CALL histwrite3d_cosp(o_crlmodis, &
          modis%Optical_Thickness_vs_ReffLiq(:,icl,:),nvertReffLiq,icl)
    enddo
  endif 
#endif
 endif !modis

 IF(.NOT.cosp_varsdefined) THEN
!$OMP MASTER
#ifndef CPP_IOIPSL_NO_OUTPUT
            DO iff=1,3
                IF (cosp_outfilekeys(iff)) THEN
                  CALL histend(cosp_nidfiles(iff))
                ENDIF ! cosp_outfilekeys
            ENDDO !  iff
#endif
! Fermeture dans phys_output_write
!#ifdef CPP_XIOS
            !On finalise l'initialisation:
            !CALL wxios_closedef()
!#endif

!$OMP END MASTER
!$OMP BARRIER
            cosp_varsdefined = .TRUE.
 END IF

    IF(cosp_varsdefined) THEN
! On synchronise les fichiers pour IOIPSL
#ifndef CPP_IOIPSL_NO_OUTPUT 
!$OMP MASTER
     DO iff=1,3
         IF (ok_sync .AND. cosp_outfilekeys(iff)) THEN
             CALL histsync(cosp_nidfiles(iff))
         ENDIF
     END DO
!$OMP END MASTER
#endif
    ENDIF  !cosp_varsdefined

    END SUBROUTINE cosp_output_write

! ug Routine pour definir itau_iocosp depuis cosp_output_write_mod:
  SUBROUTINE set_itau_iocosp(ito)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ito
      itau_iocosp = ito
  END SUBROUTINE

  SUBROUTINE histdef2d_cosp (iff,var)

    USE ioipsl
    USE dimphy
    use iophy
    USE mod_phys_lmdz_para
    USE mod_grid_phy_lmdz, ONLY: nbp_lon
    USE print_control_mod, ONLY: lunout,prt_level
#ifdef CPP_XIOS
  USE wxios
#endif

    IMPLICIT NONE

    INCLUDE "clesphys.h"

    INTEGER                          :: iff
    TYPE(ctrl_outcosp)               :: var

    REAL zstophym
    CHARACTER(LEN=20) :: typeecrit

    ! ug On récupère le type écrit de la structure:
    !       Assez moche, �|  refaire si meilleure méthode...
    IF (INDEX(var%cosp_typeecrit(iff), "once") > 0) THEN
       typeecrit = 'once'
    ELSE IF(INDEX(var%cosp_typeecrit(iff), "t_min") > 0) THEN
       typeecrit = 't_min(X)'
    ELSE IF(INDEX(var%cosp_typeecrit(iff), "t_max") > 0) THEN
       typeecrit = 't_max(X)'
    ELSE IF(INDEX(var%cosp_typeecrit(iff), "inst") > 0) THEN
       typeecrit = 'inst(X)'
    ELSE
       typeecrit = cosp_outfiletypes(iff)
    ENDIF

    IF (typeecrit=='inst(X)'.OR.typeecrit=='once') THEN
       zstophym=zoutm_cosp(iff)
    ELSE
       zstophym=zdtimemoy_cosp
    ENDIF

#ifdef CPP_XIOS
     IF (.not. ok_all_xml) then
       IF ( var%cles(iff) ) THEN
         if (prt_level >= 10) then
              WRITE(lunout,*)'Appel wxios_add_field_to_file var%name =',var%name 
         endif
        CALL wxios_add_field_to_file(var%name, 2, cosp_nidfiles(iff), cosp_outfilenames(iff), &
                                     var%description, var%unit, 1, typeecrit)
       ENDIF
     ENDIF
#endif

#ifndef CPP_IOIPSL_NO_OUTPUT 
       IF ( var%cles(iff) ) THEN
          CALL histdef (cosp_nidfiles(iff), var%name, var%description, var%unit, &
               nbp_lon,jj_nb,nhoricosp(iff), 1,1,1, -99, 32, &
               typeecrit, zstophym,zoutm_cosp(iff))
       ENDIF
#endif

  END SUBROUTINE histdef2d_cosp

 SUBROUTINE histdef3d_cosp (iff,var,nvertsave,ncols)
    USE ioipsl
    USE dimphy
    use iophy
    USE mod_phys_lmdz_para
    USE mod_grid_phy_lmdz, ONLY: nbp_lon
    USE print_control_mod, ONLY: lunout,prt_level

#ifdef CPP_XIOS
  USE wxios
#endif


    IMPLICIT NONE

    INCLUDE "clesphys.h"

    INTEGER                        :: iff, klevs
    INTEGER, INTENT(IN), OPTIONAL  :: ncols ! ug RUSTINE POUR LES variables 4D
    INTEGER, INTENT(IN)           :: nvertsave
    TYPE(ctrl_outcosp)             :: var

    REAL zstophym
    CHARACTER(LEN=20) :: typeecrit, nomi
    CHARACTER(LEN=20) :: nom
    character(len=2) :: str2
    CHARACTER(len=20) :: nam_axvert

! Axe vertical
      IF (nvertsave.eq.nvertp(iff)) THEN
          klevs=PARASOL_NREFL
          nam_axvert="sza"
      ELSE IF (nvertsave.eq.nvertisccp(iff)) THEN
          klevs=7
          nam_axvert="pressure2"
      ELSE IF (nvertsave.eq.nvertcol(iff)) THEN
          klevs=Ncolout
          nam_axvert="column"
      ELSE IF (nvertsave.eq.nverttemp(iff)) THEN
          klevs=LIDAR_NTEMP
          nam_axvert="temp"
      ELSE IF (nvertsave.eq.nvertmisr(iff)) THEN
          klevs=MISR_N_CTH
          nam_axvert="cth16"
      ELSE IF (nvertsave.eq.nvertReffIce(iff)) THEN
          klevs= numMODISReffIceBins
          nam_axvert="ReffIce"
      ELSE IF (nvertsave.eq.nvertReffLiq(iff)) THEN
          klevs= numMODISReffLiqBins
          nam_axvert="ReffLiq"
      ELSE
           klevs=Nlevout
           nam_axvert="presnivs"
      ENDIF

! ug RUSTINE POUR LES Champs 4D
      IF (PRESENT(ncols)) THEN
               write(str2,'(i2.2)')ncols
               nomi=var%name
               nom="c"//str2//"_"//nomi
      ELSE
               nom=var%name
      END IF

    ! ug On récupère le type écrit de la structure:
    !       Assez moche, �|  refaire si meilleure méthode...
    IF (INDEX(var%cosp_typeecrit(iff), "once") > 0) THEN
       typeecrit = 'once'
    ELSE IF(INDEX(var%cosp_typeecrit(iff), "t_min") > 0) THEN
       typeecrit = 't_min(X)'
    ELSE IF(INDEX(var%cosp_typeecrit(iff), "t_max") > 0) THEN
       typeecrit = 't_max(X)'
    ELSE IF(INDEX(var%cosp_typeecrit(iff), "inst") > 0) THEN
       typeecrit = 'inst(X)'
    ELSE
       typeecrit = cosp_outfiletypes(iff)
    ENDIF

    IF (typeecrit=='inst(X)'.OR.typeecrit=='once') THEN
       zstophym=zoutm_cosp(iff)
    ELSE
       zstophym=zdtimemoy_cosp
    ENDIF

#ifdef CPP_XIOS
      IF (.not. ok_all_xml) then
        IF ( var%cles(iff) ) THEN
          if (prt_level >= 10) then
              WRITE(lunout,*)'Appel wxios_add_field_to_file 3d nom variable nam_axvert = ',nom, nam_axvert 
          endif
          CALL wxios_add_field_to_file(nom, 3, cosp_nidfiles(iff), cosp_outfilenames(iff), &
                                       var%description, var%unit, 1, typeecrit, nam_axvert)
        ENDIF
      ENDIF
#endif

#ifndef CPP_IOIPSL_NO_OUTPUT
       IF ( var%cles(iff) ) THEN
          CALL histdef (cosp_nidfiles(iff), nom, var%description, var%unit, &
               nbp_lon, jj_nb, nhoricosp(iff), klevs, 1, &
               klevs, nvertsave, 32, typeecrit, &
               zstophym, zoutm_cosp(iff))
       ENDIF
#endif

  END SUBROUTINE histdef3d_cosp

 SUBROUTINE histwrite2d_cosp(var,field)
  USE dimphy
  USE mod_phys_lmdz_para
  USE ioipsl
  use iophy
  USE mod_grid_phy_lmdz, ONLY: nbp_lon
  USE print_control_mod, ONLY: lunout,prt_level

#ifdef CPP_XIOS
  USE xios, only: xios_send_field
#endif

  IMPLICIT NONE
  INCLUDE 'clesphys.h'

    TYPE(ctrl_outcosp), INTENT(IN) :: var
    REAL, DIMENSION(:), INTENT(IN) :: field

    INTEGER :: iff

    REAL,DIMENSION(klon_mpi) :: buffer_omp
    INTEGER, allocatable, DIMENSION(:) :: index2d
    REAL :: Field2d(nbp_lon,jj_nb)
    CHARACTER(LEN=20) ::  nomi, nom
    character(len=2) :: str2
    LOGICAL, SAVE  :: firstx
!$OMP THREADPRIVATE(firstx)

    IF (prt_level >= 9) WRITE(lunout,*)'Begin histrwrite2d ',var%name

  ! On regarde si on est dans la phase de définition ou d'écriture:
  IF(.NOT.cosp_varsdefined) THEN
!$OMP MASTER
      !Si phase de définition.... on définit
      CALL conf_cospoutputs(var%name,var%cles)
      DO iff=1, 3
         IF (cosp_outfilekeys(iff)) THEN
            CALL histdef2d_cosp(iff, var)
         ENDIF
      ENDDO
!$OMP END MASTER
  ELSE
    !Et sinon on.... écrit
    IF (SIZE(field)/=klon) &
  CALL abort_physic('iophy::histwrite2d_cosp','Field first DIMENSION not equal to klon',1) 

    CALL Gather_omp(field,buffer_omp)
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,Field2d)

! La boucle sur les fichiers:
      firstx=.true.
      DO iff=1, 3
           IF (var%cles(iff) .AND. cosp_outfilekeys(iff)) THEN
                ALLOCATE(index2d(nbp_lon*jj_nb))
#ifndef CPP_IOIPSL_NO_OUTPUT
        CALL histwrite(cosp_nidfiles(iff),var%name,itau_iocosp,Field2d,nbp_lon*jj_nb,index2d) 
#endif
                deallocate(index2d)
#ifdef CPP_XIOS
              IF (.not. ok_all_xml) then
                 if (firstx) then
                  if (prt_level >= 10) then
                    WRITE(lunout,*)'xios_send_field variable ',var%name
                  endif
                  CALL xios_send_field(var%name, Field2d)
                   firstx=.false.
                 endif
              ENDIF
#endif
           ENDIF
      ENDDO 

#ifdef CPP_XIOS
      IF (ok_all_xml) THEN
        if (prt_level >= 1) then
              WRITE(lunout,*)'xios_send_field variable ',var%name
        endif
       CALL xios_send_field(var%name, Field2d)
      ENDIF
#endif

!$OMP END MASTER   
  ENDIF ! vars_defined
  IF (prt_level >= 9) WRITE(lunout,*)'End histrwrite2d_cosp ',var%name
  END SUBROUTINE histwrite2d_cosp

! ug NOUVELLE VERSION DES WRITE AVEC LA BOUCLE DO RENTREE
! AI sept 2013
  SUBROUTINE histwrite3d_cosp(var, field, nverts, ncols)
  USE dimphy
  USE mod_phys_lmdz_para
  USE ioipsl
  use iophy
  USE mod_grid_phy_lmdz, ONLY: nbp_lon
  USE print_control_mod, ONLY: lunout,prt_level

#ifdef CPP_XIOS
  USE xios, only: xios_send_field
#endif


  IMPLICIT NONE
  INCLUDE 'clesphys.h'

    TYPE(ctrl_outcosp), INTENT(IN)    :: var
    REAL, DIMENSION(:,:), INTENT(IN)  :: field ! --> field(klon,:)
    INTEGER, INTENT(IN), OPTIONAL     :: ncols ! ug RUSTINE POUR LES Champs 4D.....
    INTEGER, DIMENSION(3), INTENT(IN) :: nverts

    INTEGER :: iff, k

    REAL,DIMENSION(klon_mpi,SIZE(field,2)) :: buffer_omp
    REAL :: Field3d(nbp_lon,jj_nb,SIZE(field,2))
    INTEGER :: ip, n, nlev
    INTEGER, ALLOCATABLE, DIMENSION(:) :: index3d
    CHARACTER(LEN=20) ::  nomi, nom
    character(len=2) :: str2
    LOGICAL, SAVE  :: firstx
!$OMP THREADPRIVATE(firstx)

  IF (prt_level >= 9) write(lunout,*)'Begin histrwrite3d ',var%name

! ug RUSTINE POUR LES STD LEVS.....
      IF (PRESENT(ncols)) THEN
              write(str2,'(i2.2)')ncols
              nomi=var%name
              nom="c"//str2//"_"//nomi
      ELSE
               nom=var%name
      END IF
  ! On regarde si on est dans la phase de définition ou d'écriture:
  IF(.NOT.cosp_varsdefined) THEN
      !Si phase de définition.... on définit
!$OMP MASTER
      CALL conf_cospoutputs(var%name,var%cles)
      DO iff=1, 3
        IF (cosp_outfilekeys(iff)) THEN
          CALL histdef3d_cosp(iff, var, nverts(iff), ncols)
        ENDIF
      ENDDO
!$OMP END MASTER
  ELSE
    !Et sinon on.... écrit
    IF (SIZE(field,1)/=klon) &
   CALL abort_physic('iophy::histwrite3d','Field first DIMENSION not equal to klon',1)                                  
    nlev=SIZE(field,2)


    CALL Gather_omp(field,buffer_omp)
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,field3d)

! BOUCLE SUR LES FICHIERS
     firstx=.true.
     DO iff=1, 3
        IF (var%cles(iff) .AND. cosp_outfilekeys(iff)) THEN
           ALLOCATE(index3d(nbp_lon*jj_nb*nlev))
#ifndef CPP_IOIPSL_NO_OUTPUT
    CALL histwrite(cosp_nidfiles(iff),nom,itau_iocosp,Field3d,nbp_lon*jj_nb*nlev,index3d) 
#endif

#ifdef CPP_XIOS
          IF (.not. ok_all_xml) then
           IF (firstx) THEN
               CALL xios_send_field(nom, Field3d(:,:,1:nlev))
               IF (prt_level >= 9) WRITE(lunout,*)'xios_send_field ',var%name
               firstx=.FALSE.
           ENDIF
          ENDIF
#endif
         deallocate(index3d)
        ENDIF
      ENDDO
#ifdef CPP_XIOS
    IF (ok_all_xml) THEN
     CALL xios_send_field(nom, Field3d(:,:,1:nlev))
     IF (prt_level >= 1) WRITE(lunout,*)'xios_send_field ',var%name
    ENDIF
#endif

!$OMP END MASTER   
  ENDIF ! vars_defined
  IF (prt_level >= 9) write(lunout,*)'End histrwrite3d_cosp ',nom
  END SUBROUTINE histwrite3d_cosp

! ug NOUVELLE VERSION DES WRITE AVEC LA BOUCLE DO RENTREE
! AI sept 2013
  SUBROUTINE histwrite4d_cosp(var, field)
  USE dimphy
  USE mod_phys_lmdz_para
  USE ioipsl
  use iophy
  USE mod_grid_phy_lmdz, ONLY: nbp_lon
  USE print_control_mod, ONLY: lunout,prt_level

#ifdef CPP_XIOS
  USE xios, only: xios_send_field
#endif


  IMPLICIT NONE
  INCLUDE 'clesphys.h'

    TYPE(ctrl_outcosp), INTENT(IN)    :: var
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: field ! --> field(klon,:)

    INTEGER :: iff, k

    REAL,DIMENSION(klon_mpi,SIZE(field,2),SIZE(field,3)) :: buffer_omp
    REAL :: field4d(nbp_lon,jj_nb,SIZE(field,2),SIZE(field,3))
    INTEGER :: ip, n, nlev, nlev2
    INTEGER, ALLOCATABLE, DIMENSION(:) :: index4d
    CHARACTER(LEN=20) ::  nomi, nom

  IF (prt_level >= 9) write(lunout,*)'Begin histrwrite4d ',var%name

  IF(cosp_varsdefined) THEN
    !Et sinon on.... écrit
    IF (SIZE(field,1)/=klon) &
   CALL abort_physic('iophy::histwrite3d','Field first DIMENSION not equal to klon',1)            

    nlev=SIZE(field,2)
    nlev2=SIZE(field,3)
    CALL Gather_omp(field,buffer_omp)
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,field4d)

#ifdef CPP_XIOS
!    IF (ok_all_xml) THEN
     CALL xios_send_field(var%name, Field4d(:,:,1:nlev,1:nlev2))
     IF (prt_level >= 1) WRITE(lunout,*)'xios_send_field ',var%name
!    ENDIF
#endif

!$OMP END MASTER   
  ENDIF ! vars_defined
  IF (prt_level >= 9) write(lunout,*)'End histrwrite4d_cosp ',nom
  END SUBROUTINE histwrite4d_cosp

  SUBROUTINE conf_cospoutputs(nam_var,cles_var)
!!! Lecture des noms et cles de sortie des variables dans config.def
    !   en utilisant les routines getin de IOIPSL  
    use ioipsl
    USE print_control_mod, ONLY: lunout,prt_level

    IMPLICIT NONE

   CHARACTER(LEN=20)               :: nam_var, nnam_var
   LOGICAL, DIMENSION(3)           :: cles_var

! Lecture dans config.def ou output.def de cles_var et name_var
    CALL getin('cles_'//nam_var,cles_var)
    CALL getin('name_'//nam_var,nam_var)
    IF(prt_level>10) WRITE(lunout,*)'nam_var cles_var ',nam_var,cles_var(:)

  END SUBROUTINE conf_cospoutputs

 END MODULE cosp_output_write_mod
