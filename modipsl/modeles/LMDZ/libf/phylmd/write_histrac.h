!$Id $
!***************************************
!  ECRITURE DU FICHIER :  histrac.nc
!***************************************
  IF (ecrit_tra > 0.) THEN
     
     itau_w = itau_phy + nstep + start_time * day_step / iphysiq
     
     CALL histwrite_phy(nid_tra,.FALSE.,"phis",itau_w,pphis)
     CALL histwrite_phy(nid_tra,.FALSE.,"aire",itau_w,cell_area)
     CALL histwrite_phy(nid_tra,.FALSE.,"zmasse",itau_w,zmasse)
! RomP >>>
     CALL histwrite_phy(nid_tra,.FALSE.,"sourceBE",itau_w,sourceBE)
! RomP <<<

!TRACEURS
!----------------
     DO it=1,nbtr
!!        iiq=niadv(it+2)                                                           ! jyg
        iiq=niadv(it+nqo)                                                           ! jyg

! CONCENTRATIONS
        CALL histwrite_phy(nid_tra,.FALSE.,tname(iiq),itau_w,tr_seri(:,:,it))

! TD LESSIVAGE       
      IF (lessivage .AND. aerosol(it)) THEN
        CALL histwrite_phy(nid_tra,.FALSE.,"fl"//tname(iiq),itau_w,flestottr(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_ls_"//tname(iiq),itau_w,d_tr_ls(:,:,it))
        IF(iflag_lscav .EQ. 3 .OR. iflag_lscav .EQ. 4) then
          CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_insc_"//tname(iiq),itau_w,d_tr_insc(:,:,it))
          CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_bcscav_"//tname(iiq),itau_w,d_tr_bcscav(:,:,it))
          CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_evls_"//tname(iiq),itau_w,d_tr_evapls(:,:,it))
          CALL histwrite_phy(nid_tra,.FALSE.,"qpr_ls_"//tname(iiq),itau_w,qPrls(:,it))
        ENDIF
      ENDIF

! TD THERMIQUES
        IF (iflag_thermals.gt.0) THEN
           CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_th_"//tname(iiq),itau_w,d_tr_th(:,:,it))
        ENDIF

! TD CONVECTION
        IF (iflag_con.GE.2) THEN
           CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_cv_"//tname(iiq),itau_w,d_tr_cv(:,:,it))
        ENDIF

! TD COUCHE-LIMITE
      IF (iflag_vdf_trac>=0) THEN
        CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_cl_"//tname(iiq),itau_w,d_tr_cl(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_dry_"//tname(iiq),itau_w,d_tr_dry(:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"flux_tr_dry_"//tname(iiq),itau_w,flux_tr_dry(:,it))
      ENDIF

! TD radio-decroissance
        CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_dec_"//tname(iiq),itau_w,d_tr_dec(:,:,it))

! RomP >>>
        IF (iflag_con.EQ.30) THEN
        CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_cvMA_"//tname(iiq),itau_w,dtrcvMA(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_trsp_"//tname(iiq),itau_w,d_tr_trsp(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_sscav_"//tname(iiq),itau_w,d_tr_sscav(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_sat_"//tname(iiq),itau_w,d_tr_sat(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"d_tr_uscav_"//tname(iiq),itau_w,d_tr_uscav(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"tr_pr_"//tname(iiq),itau_w,qPr(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"tr_aa_"//tname(iiq),itau_w,qPa(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"tr_mel_"//tname(iiq),itau_w,qMel(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"tr_di_"//tname(iiq),itau_w,qDi(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"tr_trspdi_"//tname(iiq),itau_w,qTrdi(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"zmfd1a_"//tname(iiq),itau_w,zmfd1a(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"zmfphi2_"//tname(iiq),itau_w,zmfphi2(:,:,it))
        CALL histwrite_phy(nid_tra,.FALSE.,"zmfdam_"//tname(iiq),itau_w,zmfdam(:,:,it))
        ENDIF
        CALL histwrite_phy(nid_tra,.FALSE.,"dtrdyn_"//tname(iiq),itau_w,d_tr_dyn(:,:,it))
! RomP <<<
     ENDDO
!---------------
!
!
! VENT (niveau 1)   
     CALL histwrite_phy(nid_tra,.FALSE.,"pyu1",itau_w,yu1)
     CALL histwrite_phy(nid_tra,.FALSE.,"pyv1",itau_w,yv1)
!
! TEMPERATURE DU SOL
     zx_tmp_fi2d(:)=ftsol(:,1)         
     CALL histwrite_phy(nid_tra,.FALSE.,"ftsol1",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=ftsol(:,2)
     CALL histwrite_phy(nid_tra,.FALSE.,"ftsol2",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=ftsol(:,3)
     CALL histwrite_phy(nid_tra,.FALSE.,"ftsol3",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=ftsol(:,4)
     CALL histwrite_phy(nid_tra,.FALSE.,"ftsol4",itau_w,zx_tmp_fi2d)
!      
! NATURE DU SOL
     zx_tmp_fi2d(:)=pctsrf(:,1)
     CALL histwrite_phy(nid_tra,.FALSE.,"psrf1",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=pctsrf(:,2)
     CALL histwrite_phy(nid_tra,.FALSE.,"psrf2",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=pctsrf(:,3)
     CALL histwrite_phy(nid_tra,.FALSE.,"psrf3",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=pctsrf(:,4)
     CALL histwrite_phy(nid_tra,.FALSE.,"psrf4",itau_w,zx_tmp_fi2d)
 
! DIVERS    
     CALL histwrite_phy(nid_tra,.FALSE.,"Mint",itau_w,Mint(:,:))
     CALL histwrite_phy(nid_tra,.FALSE.,"frac_impa",itau_w,frac_impa(:,:))    
     CALL histwrite_phy(nid_tra,.FALSE.,"frac_nucl",itau_w,frac_nucl(:,:))


     CALL histwrite_phy(nid_tra,.FALSE.,"pplay",itau_w,pplay)     
     CALL histwrite_phy(nid_tra,.FALSE.,"T",itau_w,t_seri)     
     CALL histwrite_phy(nid_tra,.FALSE.,"mfu",itau_w,pmfu)
     CALL histwrite_phy(nid_tra,.FALSE.,"mfd",itau_w,pmfd)
     CALL histwrite_phy(nid_tra,.FALSE.,"en_u",itau_w,pen_u)
     CALL histwrite_phy(nid_tra,.FALSE.,"en_d",itau_w,pen_d)
     CALL histwrite_phy(nid_tra,.FALSE.,"de_d",itau_w,pde_d)
     CALL histwrite_phy(nid_tra,.FALSE.,"de_u",itau_w,pde_u)
     CALL histwrite_phy(nid_tra,.FALSE.,"coefh",itau_w,coefh)

     IF (ok_sync) THEN
!$OMP MASTER
        CALL histsync(nid_tra)
!$OMP END MASTER
     ENDIF

  ENDIF !ecrit_tra>0.

