
! $Id: read_pstoke.F90 2345 2015-08-21 09:57:36Z emillour $



SUBROUTINE read_pstoke(irec, zrec, zklono, zklevo, airefi, phisfi, t, mfu, &
    mfd, en_u, de_u, en_d, de_d, coefh, fm_therm, en_therm, frac_impa, &
    frac_nucl, pyu1, pyv1, ftsol, psrf)

  ! ******************************************************************************
  ! Frederic HOURDIN, Abderrahmane IDELKADI
  ! Lecture des parametres physique stockes online necessaires pour
  ! recalculer offline le transport de traceurs sur une grille 2x plus fine
  ! que
  ! celle online
  ! A FAIRE : une seule routine au lieu de 2 (lectflux, redecoupe)!
  ! ******************************************************************************

  USE netcdf
  USE dimphy
  USE indice_sol_mod
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, nbp_lev

  IMPLICIT NONE

  include "netcdf.inc"

  INTEGER klono, klevo, imo, jmo
!  PARAMETER (imo=iim/2, jmo=(jjm+1)/2)
!  PARAMETER (klono=(jmo-1)*imo+2, klevo=llm)
  REAL :: phisfi(((nbp_lat/2)-1)*(nbp_lon/2)+2) !phisfi(klono)
  REAL,ALLOCATABLE :: phisfi2(:,:) !phisfi2(imo,jmo+1)
  REAL,ALLOCATABLE :: airefi2(:,:) !airefi2(imo, jmo+1)

  REAL :: mfu(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) ! mfu(klono, klevo)
  REAL :: mfd(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) ! mfd(klono, klevo)
  REAL :: en_u(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !en_u(klono, klevo)
  REAL :: de_u(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !de_u(klono, klevo)
  REAL :: en_d(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !en_d(klono, klevo)
  REAL :: de_d(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !de_d(klono, klevo)
  REAL :: coefh(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !coefh(klono, klevo)
  REAL :: fm_therm(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !fm_therm(klono, klevo)
  REAL :: en_therm(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !en_therm(klono, klevo)

  REAL,ALLOCATABLE :: mfu2(:,:,:) !mfu2(imo, jmo+1, klevo)
  REAL,ALLOCATABLE :: mfd2(:,:,:) !mfd2(imo, jmo+1, klevo)
  REAL,ALLOCATABLE :: en_u2(:,:,:) !en_u2(imo, jmo+1, klevo)
  REAL,ALLOCATABLE :: de_u2(:,:,:) !de_u2(imo, jmo+1, klevo)
  REAL,ALLOCATABLE :: en_d2(:,:,:) !en_d2(imo, jmo+1, klevo)
  REAL,ALLOCATABLE :: de_d2(:,:,:) !de_d2(imo, jmo+1, klevo)
  REAL,ALLOCATABLE :: coefh2(:,:,:) !coefh2(imo, jmo+1, klevo)
  REAL,ALLOCATABLE :: fm_therm2(:,:,:) !fm_therm2(imo, jmo+1, klevo)
  REAL,ALLOCATABLE :: en_therm2(:,:,:) !en_therm2(imo, jmo+1, klevo)

  REAL,ALLOCATABLE :: pl(:) !pl(klevo)
  INTEGER irec
  INTEGER xid, yid, zid, tid
  REAL zrec, zklono, zklevo, zim, zjm
  INTEGER ncrec, ncklono, ncklevo, ncim, ncjm

  REAL :: airefi(((nbp_lat/2)-1)*(nbp_lon/2)+2) !airefi(klono)
  CHARACTER *20 namedim

  ! !! attention !!
  ! attention il y a aussi le pb de def klono
  ! dim de phis??


  REAL :: frac_impa(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !frac_impa(klono, klevo)
  REAL :: frac_nucl(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !frac_nucl(klono, klevo)
  REAL,ALLOCATABLE :: frac_impa2(:,:,:) !frac_impa2(imo, jmo+1, klevo)
  REAL,ALLOCATABLE :: frac_nucl2(:,:,:) !frac_nucl2(imo, jmo+1, klevo)
  REAL :: pyu1(((nbp_lat/2)-1)*(nbp_lon/2)+2) !pyu1(klono)
  REAL :: pyv1(((nbp_lat/2)-1)*(nbp_lon/2)+2) !pyv1(klono)
  REAL,ALLOCATABLE :: pyu12(:,:), pyv12(:,:) !pyu12(imo, jmo+1), pyv12(imo, jmo+1)
  REAL :: ftsol(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !ftsol(klono, nbsrf)
  REAL :: psrf(((nbp_lat/2)-1)*(nbp_lon/2)+2,nbp_lev) !psrf(klono, nbsrf)
  REAL,ALLOCATABLE :: ftsol1(:),ftsol2(:) !ftsol1(klono), ftsol2(klono)
  REAL,ALLOCATABLE :: ftsol3(:),ftsol4(:) !ftsol3(klono), ftsol4(klono)
  REAL,ALLOCATABLE :: psrf1(:), psrf2(:) !psrf1(klono), psrf2(klono)
  REAL,ALLOCATABLE :: psrf3(:), psrf4(:) !psrf3(klono), psrf4(klono)
  REAL,ALLOCATABLE :: ftsol12(:,:) !ftsol12(imo, jmo+1)
  REAL,ALLOCATABLE :: ftsol22(:,:) !ftsol22(imo, jmo+1)
  REAL,ALLOCATABLE :: ftsol32(:,:) !ftsol32(imo, jmo+1)
  REAL,ALLOCATABLE :: ftsol42(:,:) !ftsol42(imo, jmo+1)
  REAL,ALLOCATABLE :: psrf12(:,:) !psrf12(imo, jmo+1)
  REAL,ALLOCATABLE :: psrf22(:,:) !psrf22(imo, jmo+1)
  REAL,ALLOCATABLE :: psrf32(:,:) !psrf32(imo, jmo+1)
  REAL,ALLOCATABLE :: psrf42(:,:) !psrf42(imo, jmo+1)
  REAL :: t(((nbp_lon/2)-1)*(nbp_lat/2)+2,nbp_lev) !t(klono, klevo)
  REAL,ALLOCATABLE :: t2(:,:,:) !t2(imo, jmo+1, klevo)
  INTEGER,SAVE :: ncidp
  INTEGER,SAVE :: varidt
  INTEGER,SAVE :: varidmfu, varidmfd, varidps, varidenu, variddeu
  INTEGER,SAVE :: varidend, varidded, varidch, varidfi, varidfn
  INTEGER,SAVE :: varidfmth, varidenth
  INTEGER,SAVE :: varidyu1, varidyv1, varidpl, varidai, varididvt
  INTEGER,SAVE :: varidfts1, varidfts2, varidfts3, varidfts4
  INTEGER,SAVE :: varidpsr1, varidpsr2, varidpsr3, varidpsr4

  INTEGER l, i
  INTEGER start(4), count(4), status
  REAL rcode
  LOGICAL,SAVE :: first=.TRUE.

  ! Allocate arrays
  imo=nbp_lon/2
  jmo=nbp_lat/2
  klono=(jmo-1)*imo+2
  klevo=nbp_lev
  
  ALLOCATE(phisfi2(imo,jmo+1))
  ALLOCATE(airefi2(imo, jmo+1))
  ALLOCATE(mfu2(imo, jmo+1, klevo))
  ALLOCATE(mfd2(imo, jmo+1, klevo))
  ALLOCATE(en_u2(imo, jmo+1, klevo))
  ALLOCATE(de_u2(imo, jmo+1, klevo))
  ALLOCATE(en_d2(imo, jmo+1, klevo))
  ALLOCATE(de_d2(imo, jmo+1, klevo))
  ALLOCATE(coefh2(imo, jmo+1, klevo))
  ALLOCATE(fm_therm2(imo, jmo+1, klevo))
  ALLOCATE(en_therm2(imo, jmo+1, klevo))
  ALLOCATE(pl(klevo))
  ALLOCATE(frac_impa2(imo, jmo+1, klevo))
  ALLOCATE(frac_nucl2(imo, jmo+1, klevo))
  ALLOCATE(pyu12(imo, jmo+1), pyv12(imo, jmo+1))
  ALLOCATE(ftsol1(klono), ftsol2(klono))
  ALLOCATE(ftsol3(klono), ftsol4(klono))
  ALLOCATE(psrf1(klono), psrf2(klono))
  ALLOCATE(psrf3(klono), psrf4(klono))
  ALLOCATE(ftsol12(imo, jmo+1))
  ALLOCATE(ftsol22(imo, jmo+1))
  ALLOCATE(ftsol32(imo, jmo+1))
  ALLOCATE(ftsol42(imo, jmo+1))
  ALLOCATE(psrf12(imo, jmo+1))
  ALLOCATE(psrf22(imo, jmo+1))
  ALLOCATE(psrf32(imo, jmo+1))
  ALLOCATE(psrf42(imo, jmo+1))
  ALLOCATE(t2(imo, jmo+1, klevo))

  ! ---------------------------------------------
  ! Initialisation de la lecture des fichiers
  ! ---------------------------------------------

  IF (irec==0) THEN

    rcode = nf90_open('phystoke.nc', nf90_nowrite, ncidp)

    rcode = nf90_inq_varid(ncidp, 'phis', varidps)
    PRINT *, 'ncidp,varidps', ncidp, varidps

    rcode = nf90_inq_varid(ncidp, 'sig_s', varidpl)
    PRINT *, 'ncidp,varidpl', ncidp, varidpl

    rcode = nf90_inq_varid(ncidp, 'aire', varidai)
    PRINT *, 'ncidp,varidai', ncidp, varidai

    ! A FAIRE: Es-il necessaire de stocke t?
    rcode = nf90_inq_varid(ncidp, 't', varidt)
    PRINT *, 'ncidp,varidt', ncidp, varidt

    rcode = nf90_inq_varid(ncidp, 'mfu', varidmfu)
    PRINT *, 'ncidp,varidmfu', ncidp, varidmfu

    rcode = nf90_inq_varid(ncidp, 'mfd', varidmfd)
    PRINT *, 'ncidp,varidmfd', ncidp, varidmfd

    rcode = nf90_inq_varid(ncidp, 'en_u', varidenu)
    PRINT *, 'ncidp,varidenu', ncidp, varidenu

    rcode = nf90_inq_varid(ncidp, 'de_u', variddeu)
    PRINT *, 'ncidp,variddeu', ncidp, variddeu

    rcode = nf90_inq_varid(ncidp, 'en_d', varidend)
    PRINT *, 'ncidp,varidend', ncidp, varidend

    rcode = nf90_inq_varid(ncidp, 'de_d', varidded)
    PRINT *, 'ncidp,varidded', ncidp, varidded

    rcode = nf90_inq_varid(ncidp, 'coefh', varidch)
    PRINT *, 'ncidp,varidch', ncidp, varidch

    ! abder (pour thermiques)
    rcode = nf90_inq_varid(ncidp, 'fm_th', varidfmth)
    PRINT *, 'ncidp,varidfmth', ncidp, varidfmth

    rcode = nf90_inq_varid(ncidp, 'en_th', varidenth)
    PRINT *, 'ncidp,varidenth', ncidp, varidenth

    rcode = nf90_inq_varid(ncidp, 'frac_impa', varidfi)
    PRINT *, 'ncidp,varidfi', ncidp, varidfi

    rcode = nf90_inq_varid(ncidp, 'frac_nucl', varidfn)
    PRINT *, 'ncidp,varidfn', ncidp, varidfn

    rcode = nf90_inq_varid(ncidp, 'pyu1', varidyu1)
    PRINT *, 'ncidp,varidyu1', ncidp, varidyu1

    rcode = nf90_inq_varid(ncidp, 'pyv1', varidyv1)
    PRINT *, 'ncidp,varidyv1', ncidp, varidyv1

    rcode = nf90_inq_varid(ncidp, 'ftsol1', varidfts1)
    PRINT *, 'ncidp,varidfts1', ncidp, varidfts1

    rcode = nf90_inq_varid(ncidp, 'ftsol2', varidfts2)
    PRINT *, 'ncidp,varidfts2', ncidp, varidfts2

    rcode = nf90_inq_varid(ncidp, 'ftsol3', varidfts3)
    PRINT *, 'ncidp,varidfts3', ncidp, varidfts3

    rcode = nf90_inq_varid(ncidp, 'ftsol4', varidfts4)
    PRINT *, 'ncidp,varidfts4', ncidp, varidfts4

    rcode = nf90_inq_varid(ncidp, 'psrf1', varidpsr1)
    PRINT *, 'ncidp,varidpsr1', ncidp, varidpsr1

    rcode = nf90_inq_varid(ncidp, 'psrf2', varidpsr2)
    PRINT *, 'ncidp,varidpsr2', ncidp, varidpsr2

    rcode = nf90_inq_varid(ncidp, 'psrf3', varidpsr3)
    PRINT *, 'ncidp,varidpsr3', ncidp, varidpsr3

    rcode = nf90_inq_varid(ncidp, 'psrf4', varidpsr4)
    PRINT *, 'ncidp,varidpsr4', ncidp, varidpsr4

    ! ID pour les dimensions

    status = nf_inq_dimid(ncidp, 'y', yid)
    status = nf_inq_dimid(ncidp, 'x', xid)
    status = nf_inq_dimid(ncidp, 'sig_s', zid)
    status = nf_inq_dimid(ncidp, 'time_counter', tid)

    ! lecture des dimensions

    status = nf_inq_dim(ncidp, yid, namedim, ncjm)
    status = nf_inq_dim(ncidp, xid, namedim, ncim)
    status = nf_inq_dim(ncidp, zid, namedim, ncklevo)
    status = nf_inq_dim(ncidp, tid, namedim, ncrec)

    zrec = ncrec
    zklevo = ncklevo
    zim = ncim
    zjm = ncjm

    zklono = zim*(zjm-2) + 2

    WRITE (*, *) 'read_pstoke : zrec = ', zrec
    WRITE (*, *) 'read_pstoke : zklevo = ', zklevo
    WRITE (*, *) 'read_pstoke : zim = ', zim
    WRITE (*, *) 'read_pstoke : zjm = ', zjm
    WRITE (*, *) 'read_pstoke : zklono = ', zklono

    ! niveaux de pression
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidpl, 1, zklevo, pl)
#else
    status = nf_get_vara_real(ncidp, varidpl, 1, zklevo, pl)
#endif

    ! lecture de aire et phis

    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = 0

    count(1) = zim
    count(2) = zjm
    count(3) = 1
    count(4) = 0

    ! phis
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidps, start, count, phisfi2)
#else
    status = nf_get_vara_real(ncidp, varidps, start, count, phisfi2)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, phisfi2, phisfi)

    ! aire
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidai, start, count, airefi2)
#else
    status = nf_get_vara_real(ncidp, varidai, start, count, airefi2)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, airefi2, airefi)
  ELSE

    PRINT *, 'ok1'

    ! ---------------------
    ! lecture des champs
    ! ---------------------

    PRINT *, 'WARNING!!! Il n y a pas de test de coherence'
    PRINT *, 'sur le nombre de niveaux verticaux dans le fichier nc'

    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = irec

    count(1) = zim
    count(2) = zjm
    count(3) = zklevo
    count(4) = 1


    ! *** Lessivage******************************************************
    ! frac_impa
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfi, start, count, frac_impa2)
#else
    status = nf_get_vara_real(ncidp, varidfi, start, count, frac_impa2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, frac_impa2, frac_impa)

    ! frac_nucl
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfn, start, count, frac_nucl2)
#else
    status = nf_get_vara_real(ncidp, varidfn, start, count, frac_nucl2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, frac_nucl2, frac_nucl)

    ! *** Temperature ******************************************************
    ! abder t
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidt, start, count, t2)
#else
    status = nf_get_vara_real(ncidp, varidt, start, count, t2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, t2, t)

    ! *** Flux pour le calcul de la convection TIEDTK ***********************
    ! mfu
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidmfu, start, count, mfu2)
#else
    status = nf_get_vara_real(ncidp, varidmfu, start, count, mfu2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, mfu2, mfu)

    ! mfd
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidmfd, start, count, mfd2)
#else
    status = nf_get_vara_real(ncidp, varidmfd, start, count, mfd2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, mfd2, mfd)

    ! en_u
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidenu, start, count, en_u2)
#else
    status = nf_get_vara_real(ncidp, varidenu, start, count, en_u2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, en_u2, en_u)

    ! de_u
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, variddeu, start, count, de_u2)
#else
    status = nf_get_vara_real(ncidp, variddeu, start, count, de_u2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, de_u2, de_u)

    ! en_d
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidend, start, count, en_d2)
#else
    status = nf_get_vara_real(ncidp, varidend, start, count, en_d2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, en_d2, en_d)

    ! de_d
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidded, start, count, de_d2)
#else
    status = nf_get_vara_real(ncidp, varidded, start, count, de_d2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, de_d2, de_d)

    ! **** Coeffecient du mellange
    ! turbulent**********************************
    ! coefh
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidch, start, count, coefh2)
#else
    status = nf_get_vara_real(ncidp, varidch, start, count, coefh2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, coefh2, coefh)

    ! *** Flux ascendant et entrant pour les
    ! Thermiques************************
    ! abder thermiques
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfmth, start, count, fm_therm2)
#else
    status = nf_get_vara_real(ncidp, varidfmth, start, count, fm_therm2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, fm_therm2, fm_therm)

#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidenth, start, count, en_therm2)
#else
    status = nf_get_vara_real(ncidp, varidenth, start, count, en_therm2)
#endif
    CALL gr_ecrit_fi(klevo, klono, imo, jmo+1, en_therm2, en_therm)

    ! *** Vitesses aux sol
    ! ******************************************************
    start(3) = irec
    start(4) = 0
    count(3) = 1
    count(4) = 0
    ! pyu1
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidyu1, start, count, pyu12)
#else
    status = nf_get_vara_real(ncidp, varidyu1, start, count, pyu12)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, pyu12, pyu1)

    ! pyv1
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidyv1, start, count, pyv12)
#else
    status = nf_get_vara_real(ncidp, varidyv1, start, count, pyv12)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, pyv12, pyv1)

    ! *** Temperature au sol ********************************************
    ! ftsol1
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfts1, start, count, ftsol12)
#else
    status = nf_get_vara_real(ncidp, varidfts1, start, count, ftsol12)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, ftsol12, ftsol1)

    ! ftsol2
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfts2, start, count, ftsol22)
#else
    status = nf_get_vara_real(ncidp, varidfts2, start, count, ftsol22)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, ftsol22, ftsol2)

    ! ftsol3
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfts3, start, count, ftsol32)
#else
    status = nf_get_vara_real(ncidp, varidfts3, start, count, ftsol32)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, ftsol32, ftsol3)

    ! ftsol4
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfts4, start, count, ftsol42)
#else
    status = nf_get_vara_real(ncidp, varidfts4, start, count, ftsol42)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, ftsol42, ftsol4)

    ! *** Nature du sol **************************************************
    ! psrf1
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidpsr1, start, count, psrf12)
#else
    status = nf_get_vara_real(ncidp, varidpsr1, start, count, psrf12)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, psrf12, psrf1)

    ! psrf2
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidpsr2, start, count, psrf22)
#else
    status = nf_get_vara_real(ncidp, varidpsr2, start, count, psrf22)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, psrf22, psrf2)

    ! psrf3
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidpsr3, start, count, psrf32)
#else
    status = nf_get_vara_real(ncidp, varidpsr3, start, count, psrf32)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, psrf32, psrf3)

    ! psrf4
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidpsr4, start, count, psrf42)
#else
    status = nf_get_vara_real(ncidp, varidpsr4, start, count, psrf42)
#endif
    CALL gr_ecrit_fi(1, klono, imo, jmo+1, psrf42, psrf4)

    DO i = 1, klono

      psrf(i, 1) = psrf1(i)
      psrf(i, 2) = psrf2(i)
      psrf(i, 3) = psrf3(i)
      psrf(i, 4) = psrf4(i)

      ftsol(i, 1) = ftsol1(i)
      ftsol(i, 2) = ftsol2(i)
      ftsol(i, 3) = ftsol3(i)
      ftsol(i, 4) = ftsol4(i)

    END DO

  END IF

  RETURN

END SUBROUTINE read_pstoke

