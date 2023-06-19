!!!! Abderrahmane Idelkadi aout 2013 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module pour definir (au 1er appel) et ecrire les variables dans les fichiers de sortie cosp
!
! R.Guzman jan 2019 (mise a jour pour COSPv2)
! On change le nom du module a "lmdz_cosp_output_write_mod" et celui de la routine a "lmdz_cosp_output_write"
! pour qu on sache qu il s agit d un module specifique a l implementation de COSP dans LMDZ
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   MODULE lmdz_cosp_output_write_mod

   USE lmdz_cosp_output_mod
   USE mod_cosp_config, only : R_UNDEF, CLOUDSAT_DBZE_BINS, SR_BINS, PARASOL_NREFL, & 
                              isccp_histPresCenters,tau_binCenters, LIDAR_NTEMP, & 
                              LIDAR_PHASE_TEMP,misr_histHgtCenters,numMISRHgtBins, & 
                              numMODISReffIceBins,reffICE_binCenters, &
                              numMODISReffLiqBins, reffLIQ_binCenters

   IMPLICIT NONE

   INTEGER, SAVE  :: itau_iocosp
!$OMP THREADPRIVATE(itau_iocosp)
   INTEGER, save        :: Nlevout, Ncolout
!$OMP THREADPRIVATE(Nlevout, Ncolout)

!  INTERFACE histwrite_cosp
!    MODULE PROCEDURE histwrite2d_cosp,histwrite3d_cosp
!  END INTERFACE

   CONTAINS

  SUBROUTINE lmdz_cosp_output_write(Nlevlmdz, Npoints, Ncolumns, itap, dtime, freq_COSP, missing_cosp, &
  				    cfg, Nlvgrid, cospOUT)


    USE ioipsl
    USE time_phylmdz_mod, ONLY: itau_phy, start_time, day_step_phy
    USE print_control_mod, ONLY: lunout,prt_level
    USE lmdz_cosp_read_outputkeys, only: cosp_config
!COSPv2
  use cosp_kinds, only: wp
  use mod_cosp,   only: cosp_outputs

#ifdef CPP_XIOS
    USE wxios, only: wxios_closedef
    USE xios, only: xios_update_calendar, xios_field_is_active
#endif
  IMPLICIT NONE  
!!! Variables d'entree
  integer               :: itap, Nlevlmdz, Ncolumns, Npoints, Nlvgrid
  real                  :: freq_COSP, dtime, missing_val, missing_cosp
  type(cosp_config)     :: cfg     ! Control outputs
  type(cosp_outputs)    :: &
       cospOUT           ! COSP simulator outputs


!!! Variables locales
  integer               :: icl,k,ip
  logical               :: ok_sync
  integer               :: itau_wcosp, iff
  real, dimension(Npoints,PARASOL_NREFL) :: parasolcrefl, Ncref

! Variables locals intermidiaires pour inverser les axes des champs 4D
! Compatibilite avec sorties CMIP
  real, dimension(Npoints,Nlevout,SR_BINS) :: tmp_fi4da_cfadL, tmp_fi4da_cfadLgr, tmp_fi4da_cfadLatlid
  real, dimension(Npoints,Nlevout,CLOUDSAT_DBZE_BINS) :: tmp_fi4da_cfadR
  real, dimension(Npoints,numMISRHgtBins,7) :: tmp_fi4da_misr

#ifdef CPP_XIOS
  missing_val=missing_cosp
#else
  missing_val=0.
#endif

  Nlevout = Nlvgrid
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
 if (cfg%Lcalipso) then
!!! AI 02 2018 
! Traitement missing_val
!!!   where(stlidar%lidarcld == R_UNDEF) stlidar%lidarcld = missing_val
!!!   where(sglidar%beta_mol == R_UNDEF) sglidar%beta_mol = missing_val  
!!!   where(sglidar%beta_tot == R_UNDEF) sglidar%beta_tot = missing_val 
!!!   where(stlidar%cldlayer == R_UNDEF) stlidar%cldlayer = missing_val
!   where(stlidar%cldtype == R_UNDEF) stlidar%cldtype = missing_val     !OPAQ
!!!   where(stlidar%cfad_sr == R_UNDEF) stlidar%cfad_sr = missing_val
! AI 11 / 2015
!!!   where(stlidar%parasolrefl == R_UNDEF) stlidar%parasolrefl = missing_val
!!!   where(stlidar%lidarcldtmp == R_UNDEF) stlidar%lidarcldtmp = missing_val
!!!   where(stlidar%cldlayerphase == R_UNDEF) stlidar%cldlayerphase = missing_val
!!!   where(stlidar%lidarcldphase == R_UNDEF) stlidar%lidarcldphase = missing_val
!   where(stlidar%lidarcldtype == R_UNDEF) stlidar%lidarcldtype = missing_val   !OPAQ
!!!   where(stlidar%lidarcldtmp == R_UNDEF) stlidar%lidarcldtmp = missing_val

!!! missing values pour toutes les valeurs R_UNDEF des variables de CALIPSO
!   where(cospOUT%calipso_betaperp_tot == R_UNDEF) cospOUT%calipso_betaperp_tot = missing_val
   where(cospOUT%calipso_beta_tot == R_UNDEF) cospOUT%calipso_beta_tot = missing_val
   where(cospOUT%calipso_tau_tot == R_UNDEF) cospOUT%calipso_tau_tot = missing_val
   where(cospOUT%calipso_lidarcldphase == R_UNDEF) cospOUT%calipso_lidarcldphase = missing_val
   where(cospOUT%calipso_lidarcldtype == R_UNDEF) cospOUT%calipso_lidarcldtype = missing_val
   where(cospOUT%calipso_cldlayerphase == R_UNDEF) cospOUT%calipso_cldlayerphase = missing_val
   where(cospOUT%calipso_lidarcldtmp == R_UNDEF)    cospOUT%calipso_lidarcldtmp = missing_val
   where(cospOUT%calipso_cfad_sr == R_UNDEF) 	    cospOUT%calipso_cfad_sr = missing_val
   where(cospOUT%calipso_lidarcld == R_UNDEF) 	    cospOUT%calipso_lidarcld = missing_val
   where(cospOUT%calipso_cldlayer == R_UNDEF)	    cospOUT%calipso_cldlayer = missing_val
   where(cospOUT%calipso_cldtype == R_UNDEF)	    cospOUT%calipso_cldtype = missing_val
   where(cospOUT%calipso_cldtypetemp == R_UNDEF)    cospOUT%calipso_cldtypetemp = missing_val
   where(cospOUT%calipso_cldtypemeanz == R_UNDEF)   cospOUT%calipso_cldtypemeanz = missing_val
   where(cospOUT%calipso_cldtypemeanzse == R_UNDEF) cospOUT%calipso_cldtypemeanzse = missing_val
   where(cospOUT%calipso_beta_mol == R_UNDEF) 	    cospOUT%calipso_beta_mol = missing_val
   where(cospOUT%calipso_temp_tot == R_UNDEF) 	    cospOUT%calipso_temp_tot = missing_val
   where(cospOUT%calipso_cldthinemis == R_UNDEF)    cospOUT%calipso_cldthinemis = missing_val
   where(cospOUT%calipso_srbval == R_UNDEF) 	    cospOUT%calipso_srbval = missing_val


!   print*,'Appel histwrite2d_cosp'
   if (cfg%Lcllcalipso) CALL histwrite2d_cosp(o_cllcalipso,cospOUT%calipso_cldlayer(:,1))
   if (cfg%Lclhcalipso) CALL histwrite2d_cosp(o_clhcalipso,cospOUT%calipso_cldlayer(:,3))
   if (cfg%Lclmcalipso) CALL histwrite2d_cosp(o_clmcalipso,cospOUT%calipso_cldlayer(:,2)) 
   if (cfg%Lcltcalipso) CALL histwrite2d_cosp(o_cltcalipso,cospOUT%calipso_cldlayer(:,4))
   if (cfg%Lclcalipso) CALL histwrite3d_cosp(o_clcalipso,cospOUT%calipso_lidarcld,nvert)
   if (cfg%Lclcalipsotmp) CALL histwrite3d_cosp(o_clcalipsotmp,cospOUT%calipso_lidarcldtmp(:,:,1),nverttemp)

   if (cfg%Lcllcalipsoice) CALL histwrite2d_cosp(o_cllcalipsoice,cospOUT%calipso_cldlayerphase(:,1,1))
   if (cfg%Lclhcalipsoice) CALL histwrite2d_cosp(o_clhcalipsoice,cospOUT%calipso_cldlayerphase(:,3,1))
   if (cfg%Lclmcalipsoice) CALL histwrite2d_cosp(o_clmcalipsoice,cospOUT%calipso_cldlayerphase(:,2,1))
   if (cfg%Lcltcalipsoice) CALL histwrite2d_cosp(o_cltcalipsoice,cospOUT%calipso_cldlayerphase(:,4,1))
   if (cfg%Lclcalipsoice) CALL histwrite3d_cosp(o_clcalipsoice,cospOUT%calipso_lidarcldphase(:,:,1),nvert)
   if (cfg%Lclcalipsotmpice) CALL histwrite3d_cosp(o_clcalipsotmpice,cospOUT%calipso_lidarcldtmp(:,:,2),nverttemp)

   if (cfg%Lcllcalipsoliq) CALL histwrite2d_cosp(o_cllcalipsoliq,cospOUT%calipso_cldlayerphase(:,1,2))
   if (cfg%Lclhcalipsoliq) CALL histwrite2d_cosp(o_clhcalipsoliq,cospOUT%calipso_cldlayerphase(:,3,2))
   if (cfg%Lclmcalipsoliq) CALL histwrite2d_cosp(o_clmcalipsoliq,cospOUT%calipso_cldlayerphase(:,2,2))
   if (cfg%Lcltcalipsoliq) CALL histwrite2d_cosp(o_cltcalipsoliq,cospOUT%calipso_cldlayerphase(:,4,2))
   if (cfg%Lclcalipsoliq) CALL histwrite3d_cosp(o_clcalipsoliq,cospOUT%calipso_lidarcldphase(:,:,2),nvert)
   if (cfg%Lclcalipsotmpliq) CALL histwrite3d_cosp(o_clcalipsotmpliq,cospOUT%calipso_lidarcldtmp(:,:,3),nverttemp)

   if (cfg%Lcllcalipsoun) CALL histwrite2d_cosp(o_cllcalipsoun,cospOUT%calipso_cldlayerphase(:,1,3))
   if (cfg%Lclhcalipsoun) CALL histwrite2d_cosp(o_clhcalipsoun,cospOUT%calipso_cldlayerphase(:,3,3))
   if (cfg%Lclmcalipsoun) CALL histwrite2d_cosp(o_clmcalipsoun,cospOUT%calipso_cldlayerphase(:,2,3))
   if (cfg%Lcltcalipsoun) CALL histwrite2d_cosp(o_cltcalipsoun,cospOUT%calipso_cldlayerphase(:,4,3))
   if (cfg%Lclcalipsoun) CALL histwrite3d_cosp(o_clcalipsoun,cospOUT%calipso_lidarcldphase(:,:,3),nvert)
   if (cfg%Lclcalipsotmpun) CALL histwrite3d_cosp(o_clcalipsotmpun,cospOUT%calipso_lidarcldtmp(:,:,4),nverttemp)

   if (cfg%Lclopaquecalipso) CALL histwrite2d_cosp(o_clopaquecalipso,cospOUT%calipso_cldtype(:,1))
   if (cfg%Lclthincalipso) CALL histwrite2d_cosp(o_clthincalipso,cospOUT%calipso_cldtype(:,2))   
   if (cfg%Lclzopaquecalipso) CALL histwrite2d_cosp(o_clzopaquecalipso,cospOUT%calipso_cldtype(:,3))
   if (cfg%Lclcalipsoopaque) CALL histwrite3d_cosp(o_clcalipsoopaque,cospOUT%calipso_lidarcldtype(:,:,1),nvert)
   if (cfg%Lclcalipsothin) CALL histwrite3d_cosp(o_clcalipsothin,cospOUT%calipso_lidarcldtype(:,:,2),nvert)  
   if (cfg%Lclcalipsozopaque) CALL histwrite3d_cosp(o_clcalipsozopaque,cospOUT%calipso_lidarcldtype(:,:,3),nvert)
   if (cfg%Lclcalipsoopacity) CALL histwrite3d_cosp(o_clcalipsoopacity,cospOUT%calipso_lidarcldtype(:,:,4),nvert)

   if (cfg%Lclopaquetemp) CALL histwrite2d_cosp(o_clopaquetemp,cospOUT%calipso_cldtypetemp(:,1))
   if (cfg%Lclthintemp) CALL histwrite2d_cosp(o_clthintemp,cospOUT%calipso_cldtypetemp(:,2))
   if (cfg%Lclzopaquetemp) CALL histwrite2d_cosp(o_clzopaquetemp,cospOUT%calipso_cldtypetemp(:,3))
   if (cfg%Lclopaquemeanz) CALL histwrite2d_cosp(o_clopaquemeanz,cospOUT%calipso_cldtypemeanz(:,1))
   if (cfg%Lclthinmeanz) CALL histwrite2d_cosp(o_clthinmeanz,cospOUT%calipso_cldtypemeanz(:,2))
   if (cfg%Lclthinemis) CALL histwrite2d_cosp(o_clthinemis,cospOUT%calipso_cldthinemis)
   if (cfg%Lclopaquemeanzse) CALL histwrite2d_cosp(o_clopaquemeanzse,cospOUT%calipso_cldtypemeanzse(:,1))
   if (cfg%Lclthinmeanzse) CALL histwrite2d_cosp(o_clthinmeanzse,cospOUT%calipso_cldtypemeanzse(:,2))
   if (cfg%Lclzopaquecalipsose) CALL histwrite2d_cosp(o_clzopaquecalipsose,cospOUT%calipso_cldtypemeanzse(:,3))


#ifdef CPP_XIOS
   do icl=1,SR_BINS
      tmp_fi4da_cfadL(:,:,icl)=cospOUT%calipso_cfad_sr(:,icl,:)
   enddo
!   if (cfg%LcfadLidarsr532) CALL histwrite4d_cosp(o_cfad_lidarsr532,stlidar%cfad_sr)
   if (cfg%LcfadLidarsr532) CALL histwrite4d_cosp(o_cfadLidarsr532,tmp_fi4da_cfadL)   !!! "_" enleve
#else
   if (cfg%LcfadLidarsr532) then
     do icl=1,SR_BINS
        CALL histwrite3d_cosp(o_cfadLidarsr532,cospOUT%calipso_cfad_sr(:,icl,:),nvert,icl)   !!! "_" enleve
     enddo
   endif
#endif

#ifdef CPP_XIOS
   if (cfg%Latb532) CALL histwrite4d_cosp(o_atb532,cospOUT%calipso_beta_tot)
#else
   if (cfg%Latb532) then  
     do icl=1,Ncolumns 
        CALL histwrite3d_cosp(o_atb532,cospOUT%calipso_beta_tot(:,icl,:),nvertmcosp,icl)
     enddo
   endif 
#endif

   if (cfg%LlidarBetaMol532) CALL histwrite3d_cosp(o_lidarBetaMol532,cospOUT%calipso_beta_mol,nvertmcosp) 

 endif !Calipso


!!!! Sorties Ground Lidar
 if (cfg%LgrLidar532) then

   where(cospOUT%grLidar532_beta_tot == R_UNDEF) cospOUT%grLidar532_beta_tot = missing_val
   where(cospOUT%grLidar532_cfad_sr == R_UNDEF) cospOUT%grLidar532_cfad_sr = missing_val
   where(cospOUT%grLidar532_lidarcld == R_UNDEF) cospOUT%grLidar532_lidarcld = missing_val
   where(cospOUT%grLidar532_cldlayer == R_UNDEF) cospOUT%grLidar532_cldlayer = missing_val
   where(cospOUT%grLidar532_beta_mol == R_UNDEF) cospOUT%grLidar532_beta_mol = missing_val
   where(cospOUT%grLidar532_srbval == R_UNDEF) cospOUT%grLidar532_srbval = missing_val

   if (cfg%LcllgrLidar532) CALL histwrite2d_cosp(o_cllgrLidar532,cospOUT%grLidar532_cldlayer(:,1))
   if (cfg%LclmgrLidar532) CALL histwrite2d_cosp(o_clmgrLidar532,cospOUT%grLidar532_cldlayer(:,2))
   if (cfg%LclhgrLidar532) CALL histwrite2d_cosp(o_clhgrLidar532,cospOUT%grLidar532_cldlayer(:,3))
   if (cfg%LcltgrLidar532) CALL histwrite2d_cosp(o_cltgrLidar532,cospOUT%grLidar532_cldlayer(:,4))

   if (cfg%LclgrLidar532) CALL histwrite3d_cosp(o_clgrLidar532,cospOUT%grLidar532_lidarcld,nvert)
   if (cfg%LlidarBetaMol532gr) CALL histwrite3d_cosp(o_lidarBetaMol532gr,cospOUT%grLidar532_beta_mol,nvertmcosp) 

#ifdef CPP_XIOS
   do icl=1,SR_BINS
      tmp_fi4da_cfadLgr(:,:,icl)=cospOUT%grLidar532_cfad_sr(:,icl,:)
   enddo
   if (cfg%LcfadLidarsr532gr) CALL histwrite4d_cosp(o_cfadLidarsr532gr,tmp_fi4da_cfadLgr)
#else
   if (cfg%LcfadLidarsr532gr) then
     do icl=1,SR_BINS
        CALL histwrite3d_cosp(o_cfadLidarsr532gr,cospOUT%grLidar532_cfad_sr(:,icl,:),nvert,icl)
     enddo
   endif
#endif

#ifdef CPP_XIOS
   if (cfg%Latb532gr) CALL histwrite4d_cosp(o_atb532gr,cospOUT%grLidar532_beta_tot)
#else
   if (cfg%Latb532gr) then  
     do icl=1,Ncolumns 
        CALL histwrite3d_cosp(o_atb532gr,cospOUT%grLidar532_beta_tot(:,icl,:),nvertmcosp,icl)
     enddo
   endif 
#endif

endif ! Ground Lidar 532 nm


!!!! Sorties Atlid
 if (cfg%Latlid) then

   where(cospOUT%atlid_beta_tot == R_UNDEF) cospOUT%atlid_beta_tot = missing_val
   where(cospOUT%atlid_cfad_sr == R_UNDEF) cospOUT%atlid_cfad_sr = missing_val
   where(cospOUT%atlid_lidarcld == R_UNDEF) cospOUT%atlid_lidarcld = missing_val
   where(cospOUT%atlid_cldlayer == R_UNDEF) cospOUT%atlid_cldlayer = missing_val
   where(cospOUT%atlid_beta_mol == R_UNDEF) cospOUT%atlid_beta_mol = missing_val
   where(cospOUT%atlid_srbval == R_UNDEF) cospOUT%atlid_srbval = missing_val

   if (cfg%Lcllatlid) CALL histwrite2d_cosp(o_cllatlid,cospOUT%atlid_cldlayer(:,1))
   if (cfg%Lclmatlid) CALL histwrite2d_cosp(o_clmatlid,cospOUT%atlid_cldlayer(:,2))
   if (cfg%Lclhatlid) CALL histwrite2d_cosp(o_clhatlid,cospOUT%atlid_cldlayer(:,3))
   if (cfg%Lcltatlid) CALL histwrite2d_cosp(o_cltatlid,cospOUT%atlid_cldlayer(:,4))

   if (cfg%Lclatlid) CALL histwrite3d_cosp(o_clatlid,cospOUT%atlid_lidarcld,nvert)
   if (cfg%LlidarBetaMol355) CALL histwrite3d_cosp(o_lidarBetaMol355,cospOUT%atlid_beta_mol,nvertmcosp) 

#ifdef CPP_XIOS
   do icl=1,SR_BINS
      tmp_fi4da_cfadLatlid(:,:,icl)=cospOUT%atlid_cfad_sr(:,icl,:)
   enddo
   if (cfg%LcfadLidarsr355) CALL histwrite4d_cosp(o_cfadlidarsr355,tmp_fi4da_cfadLatlid)
#else
   if (cfg%LcfadLidarsr355) then
     do icl=1,SR_BINS
        CALL histwrite3d_cosp(o_cfadlidarsr355,cospOUT%atlid_cfad_sr(:,icl,:),nvert,icl)
     enddo
   endif
#endif

#ifdef CPP_XIOS
   if (cfg%Latb355) CALL histwrite4d_cosp(o_atb355,cospOUT%atlid_beta_tot)
#else
   if (cfg%Latb355) then  
     do icl=1,Ncolumns 
        CALL histwrite3d_cosp(o_atb355,cospOUT%atlid_beta_tot(:,icl,:),nvertmcosp,icl)
     enddo
   endif 
#endif

endif ! Atlid


 if (cfg%Lparasol) then
   if (cfg%LparasolRefl) then
! Ces 2 diagnostics sont controles par la clef logique "LparasolRefl"

!!!   if (cfg%LparasolRefl) CALL histwrite3d_cosp(o_parasolrefl,cospOUT%parasolrefl,nvertp)
     CALL histwrite3d_cosp(o_parasolGrid_refl,cospOUT%parasolGrid_refl,nvertp)

#ifdef CPP_XIOS
     CALL histwrite4d_cosp(o_parasolPix_refl,cospOUT%parasolPix_refl)
#else 
     do icl=1,Ncolumns 
        CALL histwrite3d_cosp(o_parasolPix_refl,cospOUT%parasolPix_refl(:,icl,:),nvertp,icl)
     enddo
#endif

    endif ! LparasolRefl
   endif ! Parasol

!  if (cfg%LparasolRefl) then 
!    do k=1,PARASOL_NREFL
!     do ip=1, Npoints
!      if (stlidar%cldlayer(ip,4).gt.1.and.stlidar%parasolrefl(ip,k).ne.missing_val) then
!        parasolcrefl(ip,k)=(stlidar%parasolrefl(ip,k)-0.03*(1.-stlidar%cldlayer(ip,4)/100.))/ &
!                             (stlidar%cldlayer(ip,4)/100.)
!         Ncref(ip,k) = 1.
!      else
!         parasolcrefl(ip,k)=missing_val
!         Ncref(ip,k) = 0.
!      endif
!     enddo
!    enddo
!    CALL histwrite3d_cosp(o_Ncrefl,Ncref,nvertp)
!    CALL histwrite3d_cosp(o_parasol_crefl,parasolcrefl,nvertp)
!  endif


!!! Sorties CloudSat
 if (cfg%Lcloudsat) then

   where(cospOUT%cloudsat_Ze_tot == R_UNDEF) cospOUT%cloudsat_Ze_tot = missing_val
   where(cospOUT%cloudsat_cfad_ze == R_UNDEF) cospOUT%cloudsat_cfad_ze = missing_val
   where(cospOUT%cloudsat_precip_cover == R_UNDEF) cospOUT%cloudsat_precip_cover = missing_val
   where(cospOUT%cloudsat_pia == R_UNDEF) cospOUT%cloudsat_pia = missing_val

   if (cfg%Lptradarflag0) CALL histwrite2d_cosp(o_ptradarflag0,cospOUT%cloudsat_precip_cover(:,1))
   if (cfg%Lptradarflag1) CALL histwrite2d_cosp(o_ptradarflag1,cospOUT%cloudsat_precip_cover(:,2))
   if (cfg%Lptradarflag2) CALL histwrite2d_cosp(o_ptradarflag2,cospOUT%cloudsat_precip_cover(:,3))
   if (cfg%Lptradarflag3) CALL histwrite2d_cosp(o_ptradarflag3,cospOUT%cloudsat_precip_cover(:,4))
   if (cfg%Lptradarflag4) CALL histwrite2d_cosp(o_ptradarflag4,cospOUT%cloudsat_precip_cover(:,5))
   if (cfg%Lptradarflag5) CALL histwrite2d_cosp(o_ptradarflag5,cospOUT%cloudsat_precip_cover(:,6))
   if (cfg%Lptradarflag6) CALL histwrite2d_cosp(o_ptradarflag6,cospOUT%cloudsat_precip_cover(:,7))
   if (cfg%Lptradarflag7) CALL histwrite2d_cosp(o_ptradarflag7,cospOUT%cloudsat_precip_cover(:,8))
   if (cfg%Lptradarflag8) CALL histwrite2d_cosp(o_ptradarflag8,cospOUT%cloudsat_precip_cover(:,9))
   if (cfg%Lptradarflag9) CALL histwrite2d_cosp(o_ptradarflag9,cospOUT%cloudsat_precip_cover(:,10))
   if (cfg%Lradarpia) CALL histwrite2d_cosp(o_radarpia,cospOUT%cloudsat_pia)

#ifdef CPP_XIOS
   do icl=1,CLOUDSAT_DBZE_BINS
     tmp_fi4da_cfadR(:,:,icl)=cospOUT%cloudsat_cfad_ze(:,icl,:)
   enddo
   if (cfg%Ldbze94) CALL histwrite4d_cosp(o_dbze94,cospOUT%cloudsat_Ze_tot)
!   if (cfg%LcfadDbze94) CALL histwrite4d_cosp(o_cfadDbze94,stradar%cfad_ze)
   if (cfg%LcfadDbze94) CALL histwrite4d_cosp(o_cfadDbze94,tmp_fi4da_cfadR)
#else
   if (cfg%Ldbze94) then
    do icl=1,Ncolumns
       CALL histwrite3d_cosp(o_dbze94,cospOUT%cloudsat_Ze_tot(:,icl,:),nvert,icl)
    enddo
   endif
   if (cfg%LcfadDbze94) then
    do icl=1,CLOUDSAT_DBZE_BINS
    CALL histwrite3d_cosp(o_cfadDbze94,cospOUT%cloudsat_cfad_ze(:,icl,:),nvert,icl)
    enddo
   endif
#endif
 endif
! endif pour CloudSat


!!! Sorties combinees Cloudsat et Calipso
 if (cfg%Lcalipso .and. cfg%Lcloudsat) then
   where(cospOUT%lidar_only_freq_cloud == R_UNDEF) &
                           cospOUT%lidar_only_freq_cloud = missing_val
   where(cospOUT%cloudsat_tcc == R_UNDEF) &
                           cospOUT%cloudsat_tcc = missing_val
   where(cospOUT%cloudsat_tcc2 == R_UNDEF) &
                           cospOUT%cloudsat_tcc2 = missing_val
   where(cospOUT%radar_lidar_tcc == R_UNDEF) &
                           cospOUT%radar_lidar_tcc = missing_val

   if (cfg%Lclcalipso2) CALL histwrite3d_cosp(o_clcalipso2,cospOUT%lidar_only_freq_cloud,nvert)
   if (cfg%Lcloudsat_tcc) CALL histwrite2d_cosp(o_cloudsat_tcc,cospOUT%cloudsat_tcc)
   if (cfg%Lcloudsat_tcc2) CALL histwrite2d_cosp(o_cloudsat_tcc2,cospOUT%cloudsat_tcc2)
   if (cfg%Lcltlidarradar) CALL histwrite2d_cosp(o_cltlidarradar,cospOUT%radar_lidar_tcc)
 endif


!!! Sorties Isccp
 if (cfg%Lisccp) then
  where(cospOUT%isccp_totalcldarea == R_UNDEF) cospOUT%isccp_totalcldarea = missing_val
  where(cospOUT%isccp_meanptop == R_UNDEF) cospOUT%isccp_meanptop = missing_val
  where(cospOUT%isccp_meantaucld == R_UNDEF) cospOUT%isccp_meantaucld = missing_val
  where(cospOUT%isccp_meanalbedocld == R_UNDEF) cospOUT%isccp_meanalbedocld = missing_val
  where(cospOUT%isccp_meantb == R_UNDEF) cospOUT%isccp_meantb = missing_val
  where(cospOUT%isccp_meantbclr == R_UNDEF) cospOUT%isccp_meantbclr = missing_val
  where(cospOUT%isccp_fq == R_UNDEF) cospOUT%isccp_fq = missing_val
  where(cospOUT%isccp_boxtau == R_UNDEF) cospOUT%isccp_boxtau = missing_val
  where(cospOUT%isccp_boxptop == R_UNDEF) cospOUT%isccp_boxptop = missing_val 

!   CALL histwrite2d_cosp(o_sunlit,gbx%sunlit)
#ifdef CPP_XIOS
  if (cfg%Lclisccp) CALL histwrite4d_cosp(o_clisccp,cospOUT%isccp_fq)
#else
   if (cfg%Lclisccp) then
     do icl=1,7
       CALL histwrite3d_cosp(o_clisccp,cospOUT%isccp_fq(:,icl,:),nvertisccp,icl) 
     enddo
   endif
#endif

   if (cfg%Lboxtauisccp) CALL histwrite3d_cosp(o_boxtauisccp,cospOUT%isccp_boxtau,nvertcol)
   if (cfg%Lboxptopisccp) CALL histwrite3d_cosp(o_boxptopisccp,cospOUT%isccp_boxptop,nvertcol) 
   if (cfg%Lcltisccp) CALL histwrite2d_cosp(o_cltisccp,cospOUT%isccp_totalcldarea) 
   if (cfg%Lpctisccp) CALL histwrite2d_cosp(o_pctisccp,cospOUT%isccp_meanptop) 
   if (cfg%Ltauisccp) CALL histwrite2d_cosp(o_tauisccp,cospOUT%isccp_meantaucld) 
   if (cfg%Lalbisccp) CALL histwrite2d_cosp(o_albisccp,cospOUT%isccp_meanalbedocld) 
   if (cfg%Lmeantbisccp) CALL histwrite2d_cosp(o_meantbisccp,cospOUT%isccp_meantb) 
   if (cfg%Lmeantbclrisccp) CALL histwrite2d_cosp(o_meantbclrisccp,cospOUT%isccp_meantbclr)
 endif ! Isccp


!!! MISR simulator
 if (cfg%Lmisr) then

   if (cfg%LclMISR) then
! Ces 3 diagnostics sont controles par la clef logique "LclMISR"
   where(cospOUT%misr_fq == R_UNDEF) cospOUT%misr_fq = missing_val
!   where(cospOUT%misr_dist_model_layertops == R_UNDEF) cospOUT%misr_dist_model_layertops = missing_val
   where(cospOUT%misr_meanztop == R_UNDEF) cospOUT%misr_meanztop = missing_val
   where(cospOUT%misr_cldarea == R_UNDEF) cospOUT%misr_cldarea = missing_val

#ifdef CPP_XIOS
   do icl=1,numMISRHgtBins
      tmp_fi4da_misr(:,icl,:)=cospOUT%misr_fq(:,:,icl)
   enddo
!   if (cfg%LclMISR) CALL histwrite4d_cosp(o_clMISR,misr%fq_MISR)
!   if (cfg%LclMISR) CALL histwrite4d_cosp(o_clMISR,tmp_fi4da_misr)
   CALL histwrite4d_cosp(o_misr_fq,tmp_fi4da_misr)
#else
    do icl=1,7 
      CALL histwrite3d_cosp(o_misr_fq,cospOUT%misr_fq(:,icl,:),nvertmisr,icl)
    enddo
#endif

   CALL histwrite2d_cosp(o_misr_meanztop,cospOUT%misr_meanztop) 
   CALL histwrite2d_cosp(o_misr_cldarea,cospOUT%misr_cldarea) 
  endif ! LclMISR

 endif ! Misr


!!! Modis simulator
 if (cfg%Lmodis) then
  where(cospOUT%modis_Cloud_Fraction_Low_Mean == R_UNDEF) &
        cospOUT%modis_Cloud_Fraction_Low_Mean = missing_val
  where(cospOUT%modis_Cloud_Fraction_High_Mean == R_UNDEF) &
        cospOUT%modis_Cloud_Fraction_High_Mean = missing_val
  where(cospOUT%modis_Cloud_Fraction_Mid_Mean == R_UNDEF) &
        cospOUT%modis_Cloud_Fraction_Mid_Mean = missing_val
  where(cospOUT%modis_Cloud_Fraction_Total_Mean == R_UNDEF) &
        cospOUT%modis_Cloud_Fraction_Total_Mean = missing_val
  where(cospOUT%modis_Cloud_Fraction_Water_Mean == R_UNDEF) &
        cospOUT%modis_Cloud_Fraction_Water_Mean = missing_val
  where(cospOUT%modis_Cloud_Fraction_Ice_Mean == R_UNDEF) &
        cospOUT%modis_Cloud_Fraction_Ice_Mean = missing_val
  where(cospOUT%modis_Optical_Thickness_Total_Mean == R_UNDEF) &
        cospOUT%modis_Optical_Thickness_Total_Mean = missing_val
  where(cospOUT%modis_Optical_Thickness_Water_Mean == R_UNDEF) &
        cospOUT%modis_Optical_Thickness_Water_Mean = missing_val
  where(cospOUT%modis_Optical_Thickness_Ice_Mean == R_UNDEF) &
        cospOUT%modis_Optical_Thickness_Ice_Mean = missing_val
  where(cospOUT%modis_Cloud_Particle_Size_Water_Mean == R_UNDEF) &
        cospOUT%modis_Cloud_Particle_Size_Water_Mean = missing_val
  where(cospOUT%modis_Cloud_Particle_Size_Ice_Mean == R_UNDEF) &
        cospOUT%modis_Cloud_Particle_Size_Ice_Mean = missing_val
  where(cospOUT%modis_Cloud_Top_Pressure_Total_Mean == R_UNDEF) &
        cospOUT%modis_Cloud_Top_Pressure_Total_Mean = missing_val
  where(cospOUT%modis_Liquid_Water_Path_Mean == R_UNDEF) &
        cospOUT%modis_Liquid_Water_Path_Mean = missing_val 
  where(cospOUT%modis_Ice_Water_Path_Mean == R_UNDEF) &
        cospOUT%modis_Ice_Water_Path_Mean = missing_val
  where(cospOUT%modis_Optical_Thickness_Total_LogMean == R_UNDEF) &
          cospOUT%modis_Optical_Thickness_Total_LogMean = missing_val
  where(cospOUT%modis_Optical_Thickness_Water_LogMean == R_UNDEF) &
          cospOUT%modis_Optical_Thickness_Water_LogMean = missing_val
  where(cospOUT%modis_Optical_Thickness_Ice_LogMean == R_UNDEF) &
          cospOUT%modis_Optical_Thickness_Ice_LogMean = missing_val
    
  if (cfg%Lcllmodis) CALL histwrite2d_cosp(o_cllmodis,cospOUT%modis_Cloud_Fraction_Low_Mean)
  if (cfg%Lclhmodis) CALL histwrite2d_cosp(o_clhmodis,cospOUT%modis_Cloud_Fraction_High_Mean)
  if (cfg%Lclmmodis) CALL histwrite2d_cosp(o_clmmodis,cospOUT%modis_Cloud_Fraction_Mid_Mean)
  if (cfg%Lcltmodis) CALL histwrite2d_cosp(o_cltmodis,cospOUT%modis_Cloud_Fraction_Total_Mean)
  if (cfg%Lclwmodis) CALL histwrite2d_cosp(o_clwmodis,cospOUT%modis_Cloud_Fraction_Water_Mean)
  if (cfg%Lclimodis) CALL histwrite2d_cosp(o_climodis,cospOUT%modis_Cloud_Fraction_Ice_Mean)
  if (cfg%Ltautmodis) CALL histwrite2d_cosp(o_tautmodis,cospOUT%modis_Optical_Thickness_Total_Mean)
  if (cfg%Ltauwmodis) CALL histwrite2d_cosp(o_tauwmodis,cospOUT%modis_Optical_Thickness_Water_Mean)
  if (cfg%Ltauimodis) CALL histwrite2d_cosp(o_tauimodis,cospOUT%modis_Optical_Thickness_Ice_Mean)
  if (cfg%Ltautlogmodis) CALL histwrite2d_cosp(o_tautlogmodis,cospOUT%modis_Optical_Thickness_Total_LogMean)  
  if (cfg%Ltauwlogmodis) CALL histwrite2d_cosp(o_tauwlogmodis,cospOUT%modis_Optical_Thickness_Water_LogMean)
  if (cfg%Ltauilogmodis) CALL histwrite2d_cosp(o_tauilogmodis,cospOUT%modis_Optical_Thickness_Ice_LogMean)
  if (cfg%Lreffclwmodis) CALL histwrite2d_cosp(o_reffclwmodis,cospOUT%modis_Cloud_Particle_Size_Water_Mean)
  if (cfg%Lreffclimodis) CALL histwrite2d_cosp(o_reffclimodis,cospOUT%modis_Cloud_Particle_Size_Ice_Mean)
  if (cfg%Lpctmodis) CALL histwrite2d_cosp(o_pctmodis,cospOUT%modis_Cloud_Top_Pressure_Total_Mean)
  if (cfg%Llwpmodis) CALL histwrite2d_cosp(o_lwpmodis,cospOUT%modis_Liquid_Water_Path_Mean)
  if (cfg%Liwpmodis) CALL histwrite2d_cosp(o_iwpmodis,cospOUT%modis_Ice_Water_Path_Mean)

  if (cfg%Lclmodis) then
! Ces 3 diagnostics sont controles par la clef logique "Lclmodis"
    where(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure == R_UNDEF) &
          cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure = missing_val
    where(cospOUT%modis_Optical_Thickness_vs_ReffICE == R_UNDEF) &
          cospOUT%modis_Optical_Thickness_vs_ReffICE = missing_val
    where(cospOUT%modis_Optical_thickness_vs_ReffLIQ == R_UNDEF) &
          cospOUT%modis_Optical_thickness_vs_ReffLIQ = missing_val

#ifdef CPP_XIOS
    CALL histwrite4d_cosp(o_modis_ot_vs_ctp,cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure)
    CALL histwrite4d_cosp(o_modis_ot_vs_reffice,cospOUT%modis_Optical_Thickness_vs_ReffICE)
    CALL histwrite4d_cosp(o_modis_ot_vs_reffliq,cospOUT%modis_Optical_thickness_vs_ReffLIQ)
#else
   do icl=1,7
   CALL histwrite3d_cosp(o_modis_ot_vs_ctp, &
     cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,icl,:),nvertisccp,icl)
   CALL histwrite3d_cosp(o_modis_ot_vs_reffice, &
     cospOUT%modis_Optical_Thickness_vs_ReffICE(:,icl,:),nvertReffIce,icl)
   CALL histwrite3d_cosp(o_modis_ot_vs_reffliq, &
     cospOUT%modis_Optical_thickness_vs_ReffLIQ(:,icl,:),nvertReffLiq,icl)
   enddo
#endif

!#ifdef CPP_XIOS
!  if (cfg%Lclmodis) CALL histwrite4d_cosp(o_crimodis,modis%Optical_Thickness_vs_ReffIce)
!  if (cfg%Lclmodis) CALL histwrite4d_cosp(o_crlmodis,modis%Optical_Thickness_vs_ReffLiq)
!#else
!  if (cfg%Lclmodis) then
!    do icl=1,7
!     CALL histwrite3d_cosp(o_crimodis, &
!          modis%Optical_Thickness_vs_ReffIce(:,icl,:),nvertReffIce,icl)
!    enddo
!  endif
!  if (cfg%Lclmodis) then
!    do icl=1,7 
!     CALL histwrite3d_cosp(o_crlmodis, &
!          modis%Optical_Thickness_vs_ReffLiq(:,icl,:),nvertReffLiq,icl)
!    enddo
!  endif 
!#endif
  endif ! Lclmodis

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

    END SUBROUTINE lmdz_cosp_output_write


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
          klevs=numMISRHgtBins
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

 END MODULE lmdz_cosp_output_write_mod
