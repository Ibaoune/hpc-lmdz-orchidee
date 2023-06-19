!
! $Id: initphysto.F90 2343 2015-08-20 10:02:53Z emillour $
!
SUBROUTINE initphysto(infile,tstep,t_ops,t_wrt,fileid)
  
  USE dimphy
  USE mod_phys_lmdz_para
  USE IOIPSL
  USE iophy
  USE indice_sol_mod
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, nbp_lev
  USE time_phylmdz_mod, ONLY: day_ref, annee_ref
  
  IMPLICIT NONE

!
!   Routine d'initialisation des ecritures des fichiers histoires LMDZ
!   au format IOIPSL
!
!   Appels succesifs des routines: histbeg
!                                  histhori
!                                  histver
!                                  histdef
!                                  histend
!
!   Entree:
!
!      infile: nom du fichier histoire a creer
!      day0,anne0: date de reference
!      tstep: duree du pas de temps en seconde
!      t_ops: frequence de l'operation pour IOIPSL
!      t_wrt: frequence d'ecriture sur le fichier
!
!   Sortie:
!      fileid: ID du fichier netcdf cree
!
!   L. Fairhead, LMD, 03/99
!
! =====================================================================

!   Arguments
  CHARACTER(len=*), INTENT(IN) :: infile
  REAL, INTENT(IN)             :: tstep
  REAL, INTENT(IN)             :: t_ops
  REAL, INTENT(IN)             :: t_wrt
  INTEGER, INTENT(OUT)         :: fileid

! Variables locales
  INTEGER nhoriid, i
  INTEGER l,k
  REAL nivsigs(nbp_lev)
  INTEGER tau0
  REAL zjulian
  INTEGER iq
  INTEGER uhoriid, vhoriid, thoriid, zvertiid
  INTEGER ii,jj
  INTEGER zan, idayref
  LOGICAL ok_sync
  REAL zx_lon(nbp_lon,nbp_lat), zx_lat(nbp_lon,nbp_lat)
  CHARACTER(len=12) :: nvar

!  Initialisations
!
  ok_sync= .TRUE.
!
!  Appel a histbeg: creation du fichier netcdf et initialisations diverses
!         

  zan = annee_ref
  idayref = day_ref
  CALL ymds2ju(zan, 1, idayref, 0.0, zjulian)
  tau0 = 0
  
  CALL histbeg_phy(infile,tau0, zjulian, tstep, &
       nhoriid, fileid)

!$OMP MASTER	
!  Appel a histvert pour la grille verticale
!
  DO l=1,nbp_lev
     nivsigs(l)=REAL(l)
  ENDDO
  
  CALL histvert(fileid, 'sig_s', 'Niveaux sigma', &
       'sigma_level', &
       nbp_lev, nivsigs, zvertiid)
!
!  Appels a histdef pour la definition des variables a sauvegarder
!
  CALL histdef(fileid, "phis", "Surface geop. height", "-", &
       nbp_lon,jj_nb,nhoriid, 1,1,1, -99, 32, &
       "once", t_ops, t_wrt)
  
  CALL histdef(fileid, "aire", "Grid area", "-", &
       nbp_lon,jj_nb,nhoriid, 1,1,1, -99, 32, &
       "once", t_ops, t_wrt)

  CALL histdef(fileid, "longitudes", "longitudes", "-", &
       nbp_lon,jj_nb,nhoriid, 1,1,1, -99, 32, &
       "once", t_ops, t_wrt)

  CALL histdef(fileid, "latitudes", "latitudes", "-", &
       nbp_lon,jj_nb,nhoriid, 1,1,1, -99, 32, &
       "once", t_ops, t_wrt)
! T 
  CALL histdef(fileid, 't', 'Temperature', 'K', nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! mfu 
  CALL histdef(fileid, 'mfu', 'flx m. pan. mt', 'kg m/s',nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! mfd 
  CALL histdef(fileid, 'mfd', 'flx m. pan. des', 'kg m/s',nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! en_u 
  CALL histdef(fileid, 'en_u', 'flx ent pan mt', 'kg m/s', nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! de_u 
  CALL histdef(fileid, 'de_u', 'flx det pan mt', 'kg m/s',nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! en_d 
  CALL histdef(fileid, 'en_d', 'flx ent pan dt', 'kg m/s', nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! de_d 
  CALL histdef(fileid, 'de_d', 'flx det pan dt', 'kg m/s', nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! coefh
  CALL histdef(fileid, "coefh", " ", " ", nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid,32, "inst(X)", t_ops, t_wrt)
! fm_th
  CALL histdef(fileid, "fm_th", " ", " ",nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid,32, "inst(X)", t_ops, t_wrt)
! en_th
  CALL histdef(fileid, "en_th", " ", " ",nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid,32, "inst(X)", t_ops, t_wrt)
! frac_impa
  CALL histdef(fileid, 'frac_impa', ' ', ' ',nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! frac_nucl
  CALL histdef(fileid, 'frac_nucl', ' ', ' ',nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! pyu1
  CALL histdef(fileid, "pyu1", " ", " ", nbp_lon,jj_nb,nhoriid, &
       1,1,1, -99, 32, "inst(X)", t_ops, t_wrt)
! pyv1
  CALL histdef(fileid, "pyv1", " ", " ", nbp_lon,jj_nb,nhoriid, &
       1,1,1, -99, 32,"inst(X)", t_ops, t_wrt)    
! ftsol1
  CALL histdef(fileid, "ftsol1", " ", " ",nbp_lon, jj_nb, nhoriid, &
       1, 1,1, -99,32, "inst(X)", t_ops, t_wrt)
! ftsol2
  CALL histdef(fileid, "ftsol2", " ", " ",nbp_lon, jj_nb, nhoriid, &
       1, 1,1, -99,32, "inst(X)", t_ops, t_wrt)
! ftsol3
  CALL histdef(fileid, "ftsol3", " ", " ", nbp_lon, jj_nb, nhoriid, &
       1, 1,1, -99,32, "inst(X)", t_ops, t_wrt)
! ftsol4
  CALL histdef(fileid, "ftsol4", " ", " ",nbp_lon, jj_nb, nhoriid, &
       1, 1,1, -99, 32, "inst(X)", t_ops, t_wrt)
! psrf1
  CALL histdef(fileid, "psrf1", " ", " ",nbp_lon, jj_nb, nhoriid, &
       1, 1, 1, -99,32, "inst(X)", t_ops, t_wrt)
! psrf2
  CALL histdef(fileid, "psrf2", " ", " ",nbp_lon, jj_nb, nhoriid, &
       1, 1, 1, -99, 32, "inst(X)", t_ops, t_wrt)
! psrf3
  CALL histdef(fileid, "psrf3", " ", " ",nbp_lon, jj_nb, nhoriid, &
       1, 1, 1, -99, 32, "inst(X)", t_ops, t_wrt)
! psrf4
  CALL histdef(fileid, "psrf4", " ", " ", nbp_lon, jj_nb, nhoriid, &
       1, 1, 1, -99,32, "inst(X)", t_ops, t_wrt)
! sh
  CALL histdef(fileid, 'sh', '', '', nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! da
  CALL histdef(fileid, 'da', '', '', nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! mp
  CALL histdef(fileid, 'mp', '', '', nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! upwd
  CALL histdef(fileid, 'upwd', '', '', nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! dnwd
  CALL histdef(fileid, 'dnwd', '', '', nbp_lon, jj_nb, nhoriid, &
       nbp_lev, 1, nbp_lev, zvertiid, 32, 'inst(X)', t_ops, t_wrt)

! phi
  DO k=1,nbp_lev
     IF (k<10) THEN
        WRITE(nvar,'(i1)') k
     ELSE IF (k<100) THEN
        WRITE(nvar,'(i2)') k
     ELSE
        WRITE(nvar,'(i3)') k
     END IF
     nvar='phi_lev'//trim(nvar)
     
     CALL histdef(fileid, nvar, '', '', nbp_lon, jj_nb, nhoriid, &
          nbp_lev, 1, nbp_lev, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
  END DO

  CALL histend(fileid)
  IF (ok_sync) CALL histsync
!$OMP END MASTER
	
END SUBROUTINE initphysto
