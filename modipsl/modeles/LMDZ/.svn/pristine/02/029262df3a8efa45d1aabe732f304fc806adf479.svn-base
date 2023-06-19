!
! $Id $
!
  IF (ecrit_tra>0.) THEN
!$OMP MASTER 
!!!     CALL ymds2ju(annee_ref, 1, day_ref, 0.0, zjulian)
! correction pour l heure initiale                               !jyg 
!                                                               !jyg 
      CALL ymds2ju(annee_ref, 1, day_ref, hour, zjulian)         !jyg

     CALL histbeg_phy("histrac", itau_phy, zjulian, pdtphys,nhori, nid_tra)
     CALL histvert(nid_tra, "presnivs", "Vertical levels", "Pa",klev, presnivs, nvert,"down")

     zsto = pdtphys
     zout = ecrit_tra
     CALL histdef(nid_tra, "phis", "Surface geop. height", "-",   &
          iim,jj_nb,nhori, 1,1,1, -99, 32,"once",  zsto,zout)
     CALL histdef(nid_tra, "aire", "Grid area", "-",              &
          iim,jj_nb,nhori, 1,1,1, -99, 32,"once",  zsto,zout)
     CALL histdef(nid_tra, "zmasse", "column density of air in cell", &
          "kg m-2", iim, jj_nb, nhori, klev, 1, klev, nvert, 32, "ave(X)", &
          zsto,zout)
! RomP >>>
     CALL histdef(nid_tra, "sourceBE", "source 7Be", &
          "at/kgA/s", iim, jj_nb, nhori, klev, 1, klev, nvert, 32, "ave(X)", &
          zsto,zout)
! RomP <<<

!TRACEURS
!----------------
     DO it = 1,nbtr
!!        iiq = niadv(it+2)                                                         ! jyg
        iiq = niadv(it+nqo)                                                         ! jyg

! CONCENTRATIONS
        CALL histdef(nid_tra, tname(iiq), ttext(iiq), "U/kga",    &
             iim,jj_nb,nhori, klev,1,klev,nvert, 32,"ave(X)", zsto,zout)

! TD LESSIVAGE
        IF (lessivage .AND. aerosol(it)) THEN
           CALL histdef(nid_tra, "fl"//tname(iiq),"Flux "//ttext(iiq), &
                "at/m2/s",iim,jj_nb,nhori, klev,1,klev,nvert, 32,       &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "d_tr_ls_"//tname(iiq),      &
                "tendance lessivage large scale"// ttext(iiq), "?",&
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,       &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "d_tr_insc_"//tname(iiq),      &
                "tendance lessivage large scale"// ttext(iiq), "?",&
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,       &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "d_tr_bcscav_"//tname(iiq),      &
                "tendance lessivage large scale"// ttext(iiq), "?",&
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,       &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "d_tr_evls_"//tname(iiq),      &
                "tendance lessivage large scale"// ttext(iiq), "?",&
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,       &
                "ave(X)", zsto,zout)
!  Tracer concentration in LS precipitation at surface
           CALL histdef(nid_tra, "qpr_ls_"//tname(iiq),       &
                "concentration in LS precip"// ttext(iiq), "at/kgw", &
                iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
                "ave(X)", zsto,zout)     
                 END IF

! TD THERMIQUES
        IF (iflag_thermals.gt.0) THEN
           CALL histdef(nid_tra, "d_tr_th_"//tname(iiq),      &
                "tendance thermique"// ttext(iiq), "?",       &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,       &
                "ave(X)", zsto,zout)
        ENDIF

! TD CONVECTION
        IF (iflag_con.GE.2) THEN
           CALL histdef(nid_tra, "d_tr_cv_"//tname(iiq),   &
                "tendance convection"// ttext(iiq), "?",   &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
        ENDIF

! RomP >>>
        IF (iflag_con.EQ.30) THEN
           CALL histdef(nid_tra, "d_tr_cvMA_"//tname(iiq),   &
                "tendance convection"// ttext(iiq), "?",&
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "d_tr_trsp_"//tname(iiq),   &
                "tendance transport "// ttext(iiq), "at/kga",   &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "d_tr_sscav_"//tname(iiq),   &
                "tendance lessivage flux satures "// ttext(iiq), "at/kga",   &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "d_tr_sat_"//tname(iiq),   &
                "tendance flux satures "// ttext(iiq), "at/kga",  &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "d_tr_uscav_"//tname(iiq),  &
                "tendance flux insatures "// ttext(iiq), "at/kga", &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "tr_pr_"//tname(iiq),  &
                "concentration dans precip"// ttext(iiq), "at/kga", &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "tr_aa_"//tname(iiq),  &
                "concentration precip issu AA"// ttext(iiq), "at/kga", &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "tr_mel_"//tname(iiq),  &
                "concentration precip issu melange"// ttext(iiq), "at/kga", &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "tr_di_"//tname(iiq),  &
                "concentration dans descente insaturee"// ttext(iiq), "at/kga", &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "tr_trspdi_"//tname(iiq),  &
                "conc descente insaturee MA"// ttext(iiq), "at/kga", &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "zmfd1a_"//tname(iiq),  &
                "zmfd1a"// ttext(iiq), "_", &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "zmfphi2_"//tname(iiq),  &
                "zmfphi2"// ttext(iiq), "_", &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
           CALL histdef(nid_tra, "zmfdam_"//tname(iiq),  &
                "zmfdam"// ttext(iiq), "_", &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
          ENDIF
! RomP <<<
           CALL histdef(nid_tra, "dtrdyn_"//tname(iiq),  &
                "td dyn tra"// ttext(iiq), "at/kga", &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,  &
                "ave(X)", zsto,zout)
! TD decroissance radioactive
           CALL histdef(nid_tra, "d_tr_dec_"//tname(iiq),   &
                "tendance decroi radio "// ttext(iiq), "",  &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,  &
                "ave(X)", zsto,zout)

! TD COUCHE-LIMITE
      IF (iflag_vdf_trac>=0) THEN
        CALL histdef(nid_tra, "d_tr_cl_"//tname(iiq),      &
             "tendance couche limite"// ttext(iiq), "?",   &
             iim,jj_nb,nhori, klev,1,klev,nvert, 32,       &
             "ave(X)", zsto,zout)
!  Dry deposit (1st layer and surface)
        CALL histdef(nid_tra, "d_tr_dry_"//tname(iiq),       &
             "tendancy dry deposit"// ttext(iiq), "at/kga/step", &
             iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
             "ave(X)", zsto,zout)     
        CALL histdef(nid_tra, "flux_tr_dry_"//tname(iiq),       &
             "dry deposit at surf (downward)"// ttext(iiq), "at/m2/step", &
             iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
             "ave(X)", zsto,zout)     
      ENDIF
     ENDDO

     CALL histdef(nid_tra, "Mint", "Mint","",         &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "inst(X)", zout,zout)
     CALL histdef(nid_tra, "frac_impa", "frac_impa","",         &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "inst(X)", zout,zout)
     CALL histdef(nid_tra, "frac_nucl", "frac_nucl","",         &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "inst(X)", zout,zout)
!---------------   
!
! VENT (niveau 1)
     CALL histdef(nid_tra, "pyu1", "Vent niv 1", "-",      &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)     
     CALL histdef(nid_tra, "pyv1", "Vent niv 1", "-",      &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)

! TEMPERATURE DU SOL
     CALL histdef(nid_tra, "ftsol1", "temper sol", "-",    &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)
     CALL histdef(nid_tra, "ftsol2", "temper sol", "-",    &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)
     CALL histdef(nid_tra, "ftsol3", "temper sol", "-",    &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst",  zout,zout)
     CALL histdef(nid_tra, "ftsol4", "temper sol", "-",    &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)

! NATURE DU SOL
     CALL histdef(nid_tra, "psrf1", "nature sol", "-",     &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)
     CALL histdef(nid_tra, "psrf2", "nature sol", "-",     &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)
     CALL histdef(nid_tra, "psrf3", "nature sol", "-",     &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)
     CALL histdef(nid_tra, "psrf4", "nature sol", "-",     &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 & 
          "inst(X)",  zout,zout)
! DIVERS
     CALL histdef(nid_tra, "pplay", "pressure","-",        &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "inst(X)", zout,zout)
     CALL histdef(nid_tra, "T", "temperature","K",         &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "inst(X)", zout,zout)
     CALL histdef(nid_tra, "mfu", "flux u mont","-",       &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "mfd", "flux u decen","-",      &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "en_u", "flux u mont","-",      &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "en_d", "flux u mont","-",      &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "de_d", "flux u mont","-",      &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "de_u", "flux u decen","-",     &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "coefh", "turbulent coef","-",  &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)   
     
     CALL histend(nid_tra)
!$OMP END MASTER
  END IF ! ecrit_tra>0.
  
