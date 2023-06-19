!
! $Header$
!
!IM on initialise les variables 
!
!       missing_val=nf90_fill_real
!
       IF (.not. ok_all_xml) then
        CALL ini_undefSTD(itap,itapm1)
       ENDIF
!
!IM on interpole les champs sur les niveaux STD de pression 
!IM a chaque pas de temps de la physique
!
!-------------------------------------------------------c
! positionnement de l'argument logique a .false.        c
! pour ne pas recalculer deux fois la meme chose !      c
! a cet effet un appel a plevel_new a ete deplace       c
! a la fin de la serie d'appels                         c
! la boucle 'DO k=1, nlevSTD' a ete internalisee        c
! dans plevel_new, d'ou la creation de cette routine... c
!-------------------------------------------------------c
!
        CALL plevel_new(klon,klev,nlevSTD,.true.,pplay,rlevSTD, &
                    t_seri,tlevSTD)
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   u_seri,ulevSTD)
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   v_seri,vlevSTD)
!

!
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zphi/RG,philevSTD)
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   qx(:,:,ivap),qlevSTD)
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_rh*100.,rhlevSTD)
!
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=u_seri(i,l)*v_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,uvSTD)
!
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=v_seri(i,l)*q_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,vqSTD)
!
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=v_seri(i,l)*t_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,vTSTD)
!
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=omega(i,l)*qx(i,l,ivap)
         ENDDO !i 
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,wqSTD)
!
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=v_seri(i,l)*zphi(i,l)/RG
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,vphiSTD)
!
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=omega(i,l)*t_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,wTSTD)
!
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=u_seri(i,l)*u_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,u2STD)
!
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=v_seri(i,l)*v_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,v2STD)
!
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=t_seri(i,l)*t_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,T2STD)

!
      zx_tmp_fi3d(:,:)=wo(:,:,1) * dobson_u * 1e3 / zmasse / rmo3 * rmd
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,O3STD)
!
      if (read_climoz == 2) THEN
      zx_tmp_fi3d(:,:)=wo(:,:,2) * dobson_u * 1e3 / zmasse / rmo3 * rmd
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD, &
                   zx_tmp_fi3d,O3daySTD)
      endif
!
        DO l=1, klev
        DO i=1, klon
         zx_tmp_fi3d(i,l)=paprs(i,l)
        ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.true.,zx_tmp_fi3d,rlevSTD, &
                   omega,wlevSTD)
!
!IM on somme les valeurs toutes les freq_calNMC secondes
!IM on moyenne a la fin du mois, du jour ou toutes les 6h
!
       IF (.not. ok_all_xml) then
        CALL undefSTD(itap, read_climoz)
        CALL moy_undefSTD(itap,itapm1)
       ENDIF
!
       CALL plevel(klon,klev,.true.,pplay,50000., &
                    zphi/RG,geo500)

!IM on interpole a chaque pas de temps le SWup(clr) et SWdn(clr) a 200 hPa
!
      CALL plevel(klon,klevp1,.true.,paprs,20000., &
           swdn0,SWdn200clr)
      CALL plevel(klon,klevp1,.false.,paprs,20000., &
           swdn,SWdn200)
      CALL plevel(klon,klevp1,.false.,paprs,20000., &
           swup0,SWup200clr)
      CALL plevel(klon,klevp1,.false.,paprs,20000., &
           swup,SWup200)
!
      CALL plevel(klon,klevp1,.false.,paprs,20000., &
           lwdn0,LWdn200clr)
      CALL plevel(klon,klevp1,.false.,paprs,20000., &
           lwdn,LWdn200)
      CALL plevel(klon,klevp1,.false.,paprs,20000., &
           lwup0,LWup200clr)
      CALL plevel(klon,klevp1,.false.,paprs,20000., &
           lwup,LWup200)
!
      twriteSTD(:,:,1)=tsumSTD(:,:,1)
      qwriteSTD(:,:,1)=qsumSTD(:,:,1)
      rhwriteSTD(:,:,1)=rhsumSTD(:,:,1)
      phiwriteSTD(:,:,1)=phisumSTD(:,:,1)
      uwriteSTD(:,:,1)=usumSTD(:,:,1)
      vwriteSTD(:,:,1)=vsumSTD(:,:,1)
      wwriteSTD(:,:,1)=wsumSTD(:,:,1)

      twriteSTD(:,:,2)=tsumSTD(:,:,2)
      qwriteSTD(:,:,2)=qsumSTD(:,:,2)
      rhwriteSTD(:,:,2)=rhsumSTD(:,:,2)
      phiwriteSTD(:,:,2)=phisumSTD(:,:,2)
      uwriteSTD(:,:,2)=usumSTD(:,:,2)
      vwriteSTD(:,:,2)=vsumSTD(:,:,2)
      wwriteSTD(:,:,2)=wsumSTD(:,:,2)

      twriteSTD(:,:,3)=tlevSTD(:,:)
      qwriteSTD(:,:,3)=qlevSTD(:,:)
      rhwriteSTD(:,:,3)=rhlevSTD(:,:)
      phiwriteSTD(:,:,3)=philevSTD(:,:)
      uwriteSTD(:,:,3)=ulevSTD(:,:)
      vwriteSTD(:,:,3)=vlevSTD(:,:)
      wwriteSTD(:,:,3)=wlevSTD(:,:)

      twriteSTD(:,:,4)=tlevSTD(:,:)
      qwriteSTD(:,:,4)=qlevSTD(:,:)
      rhwriteSTD(:,:,4)=rhlevSTD(:,:)
      phiwriteSTD(:,:,4)=philevSTD(:,:)
      uwriteSTD(:,:,4)=ulevSTD(:,:)
      vwriteSTD(:,:,4)=vlevSTD(:,:)
      wwriteSTD(:,:,4)=wlevSTD(:,:)
!
!IM initialisation 5eme fichier de sortie 
      twriteSTD(:,:,5)=tlevSTD(:,:)
      qwriteSTD(:,:,5)=qlevSTD(:,:)
      rhwriteSTD(:,:,5)=rhlevSTD(:,:)
      phiwriteSTD(:,:,5)=philevSTD(:,:)
      uwriteSTD(:,:,5)=ulevSTD(:,:)
      vwriteSTD(:,:,5)=vlevSTD(:,:)
      wwriteSTD(:,:,5)=wlevSTD(:,:)
!
!IM initialisation 6eme fichier de sortie 
      twriteSTD(:,:,6)=tlevSTD(:,:)
      qwriteSTD(:,:,6)=qlevSTD(:,:)
      rhwriteSTD(:,:,6)=rhlevSTD(:,:)
      phiwriteSTD(:,:,6)=philevSTD(:,:)
      uwriteSTD(:,:,6)=ulevSTD(:,:)
      vwriteSTD(:,:,6)=vlevSTD(:,:)
      wwriteSTD(:,:,6)=wlevSTD(:,:)
!IM for NMC files
      DO n=1, nlevSTD3
       DO k=1, nlevSTD
        if(rlevSTD3(n).EQ.rlevSTD(k)) THEN
         twriteSTD3(:,n)=tlevSTD(:,k)
         qwriteSTD3(:,n)=qlevSTD(:,k)
         rhwriteSTD3(:,n)=rhlevSTD(:,k)
         phiwriteSTD3(:,n)=philevSTD(:,k)
         uwriteSTD3(:,n)=ulevSTD(:,k)
         vwriteSTD3(:,n)=vlevSTD(:,k)
         wwriteSTD3(:,n)=wlevSTD(:,k)
        endif !rlevSTD3(n).EQ.rlevSTD(k)
       ENDDO 
      ENDDO 
!
      DO n=1, nlevSTD8
       DO k=1, nlevSTD
        if(rlevSTD8(n).EQ.rlevSTD(k)) THEN
         tnondefSTD8(:,n)=tnondef(:,k,2)
         twriteSTD8(:,n)=tsumSTD(:,k,2)
         qwriteSTD8(:,n)=qsumSTD(:,k,2)
         rhwriteSTD8(:,n)=rhsumSTD(:,k,2)
         phiwriteSTD8(:,n)=phisumSTD(:,k,2)
         uwriteSTD8(:,n)=usumSTD(:,k,2)
         vwriteSTD8(:,n)=vsumSTD(:,k,2)
         wwriteSTD8(:,n)=wsumSTD(:,k,2)
        endif !rlevSTD8(n).EQ.rlevSTD(k)
       ENDDO 
      ENDDO 
