!
! $Header$
!
      IF (is_sequential) THEN
      
      IF (type_run.EQ."AMIP") THEN
!
      ndex2d = 0
      itau_w = itau_phy + itap + start_time * day_step_phy
!
! Champs 2D:
!
      pi = ACOS(-1.)
      pir = 4.0*ATAN(1.0) / 180.0
!
      DO i=1, klon
       zx_tmp_fi2d(i)=(topsw(i)-toplw(i))
      ENDDO
!
      ok_msk=.FALSE.
      msk(1:klon)=pctsrf(1:klon,is_ter)
      CALL moyglo_pondaire(klon, zx_tmp_fi2d, cell_area,  &
           ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
!
      CALL gr_fi_ecrit(1, klon,nbp_lon,nbp_lat, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"bilTOA",itau_w, &
                     zx_tmp_2d,nbp_lon*nbp_lat,ndex2d)
!
      ok_msk=.FALSE.
      CALL moyglo_pondaire(klon, bils, cell_area,  &
           ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
!
      CALL gr_fi_ecrit(1, klon,nbp_lon,nbp_lat, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"bils",itau_w, &
                     zx_tmp_2d,nbp_lon*nbp_lat,ndex2d)
!
      DO k=1, klev
      DO i=1, klon
!IM 080904    zx_tmp_fi3d(i,k)=u(i,k)**2+v(i,k)**2
       zx_tmp_fi3d(i,k)=(u(i,k)**2+v(i,k)**2)/2.
      ENDDO
      ENDDO
!
      CALL moyglo_pondaima(klon, klev, zx_tmp_fi3d,  &
           cell_area, paprs, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
!
      CALL gr_fi_ecrit(1, klon,nbp_lon,nbp_lat, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"ecin",itau_w, &
                     zx_tmp_2d,nbp_lon*nbp_lat,ndex2d) 
!
!
!
      CALL gr_fi_ecrit(1,klon,nbp_lon,nbp_lat,cell_area,zx_tmp_2d)
      airetot=0.
      DO i=1, klon
       airetot=airetot+cell_area(i)
      ENDDO
!     IF(itap.EQ.1) PRINT*,'airetotphy=',airetot
!
      airetot=0.
      DO j=1, nbp_lat
       DO i=1, nbp_lon
        airetot=airetot+zx_tmp_2d(i,j)
       ENDDO
      ENDDO
!
!     IF(itap.EQ.1) PRINT*,'airetotij=',airetot,
!    $ '4piR2',4.*pi*RA*RA
!
      zx_tmp_fi2d(1:klon)=aam/airetot
      CALL gr_fi_ecrit(1,klon,nbp_lon,nbp_lat,zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"momang",itau_w,zx_tmp_2d, &
                     nbp_lon*nbp_lat,ndex2d)
!
      zx_tmp_fi2d(1:klon)=torsfc/airetot
      CALL gr_fi_ecrit(1,klon,nbp_lon,nbp_lat,zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"torsfc",itau_w,zx_tmp_2d, &
                     nbp_lon*nbp_lat,ndex2d)
!
!IM 151004 END
!
      CALL moyglo_pondmass(klon, klev, t_seri, &
           cell_area, paprs, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
!
      CALL gr_fi_ecrit(1,klon,nbp_lon,nbp_lat,zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"tamv",itau_w, &
                     zx_tmp_2d,nbp_lon*nbp_lat,ndex2d)
!
      ok_msk=.FALSE.
      CALL moyglo_pondaire(klon, paprs(:,1), cell_area,  &
           ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
!
      CALL gr_fi_ecrit(1, klon,nbp_lon,nbp_lat, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"psol",itau_w, &
                     zx_tmp_2d,nbp_lon*nbp_lat,ndex2d)
!
      ok_msk=.FALSE.
      CALL moyglo_pondaire(klon, evap, cell_area,  &
           ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
!
      CALL gr_fi_ecrit(1, klon,nbp_lon,nbp_lat, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"evap",itau_w, &
                     zx_tmp_2d,nbp_lon*nbp_lat,ndex2d)
!
!     DO i=1, klon
!      zx_tmp_fi2d(i)=SnowFrac(i,is_ter)
!     ENDDO
!
!     ok_msk=.TRUE.
!     msk(1:klon)=pctsrf(1:klon,is_ter)
!     CALL moyglo_pondaire(klon, zx_tmp_fi2d, cell_area, 
!    .                     ok_msk, msk, moyglo)
!     zx_tmp_fi2d(1:klon)=moyglo
!
!     CALL gr_fi_ecrit(1, klon,nbp_lon,nbp_lat,zx_tmp_fi2d,zx_tmp_2d)
!     CALL histwrite(nid_day_seri,"SnowFrac",
!    .               itau_w,zx_tmp_2d,nbp_lon*nbp_lat,ndex2d) 
!
!     DO i=1, klon
!IM 080904    zx_tmp_fi2d(i)=zsnow_mass(i)/330.*rowl
!      zx_tmp_fi2d(i)=zsnow_mass(i)
!     ENDDO
!
!IM 140904   ok_msk=.FALSE.
!     ok_msk=.TRUE.
!     msk(1:klon)=pctsrf(1:klon,is_ter)
!     CALL moyglo_pondaire(klon, zx_tmp_fi2d, cell_area, 
!    .     ok_msk, msk, moyglo)
!     zx_tmp_fi2d(1:klon)=moyglo
!
!     CALL gr_fi_ecrit(1, klon,nbp_lon,nbp_lat,zx_tmp_fi2d,zx_tmp_2d)
!     CALL histwrite(nid_day_seri,"snow_depth",itau_w,
!    .               zx_tmp_2d,nbp_lon*nbp_lat,ndex2d)
!
      DO i=1, klon
       zx_tmp_fi2d(i)=ftsol(i,is_oce)
      ENDDO
!
      ok_msk=.TRUE.
      msk(1:klon)=pctsrf(1:klon,is_oce)
      CALL moyglo_pondaire(klon, zx_tmp_fi2d, cell_area,  &
           ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
!
      CALL gr_fi_ecrit(1, klon,nbp_lon,nbp_lat, zx_tmp_fi2d, zx_tmp_2d)
      CALL histwrite(nid_day_seri,"tsol_"//clnsurf(is_oce), &
                     itau_w,zx_tmp_2d,nbp_lon*nbp_lat,ndex2d) 
!
!=================================================================
!=================================================================
!=================================================================
!
      if (ok_sync) then
        call histsync(nid_day_seri)
      endif
!
      ENDIF !fin test sur type_run.EQ."AMIP"
      
      ENDIF  ! mono_cpu
