
! $Id: read_pstoke0.F90 2345 2015-08-21 09:57:36Z emillour $



SUBROUTINE read_pstoke0(irec, zrec, zkon, zkev, airefi, phisfi, t, mfu, mfd, &
    en_u, de_u, en_d, de_d, coefh, fm_therm, en_therm, frac_impa, frac_nucl, &
    pyu1, pyv1, ftsol, psrf)

  ! ******************************************************************************
  ! Frederic HOURDIN, Abderrahmane IDELKADI
  ! Lecture des parametres physique stockes online necessaires pour
  ! recalculer offline le transport des traceurs sur la meme grille que
  ! online
  ! A FAIRE : une seule routine au lieu de 2 (lectflux, redecoupe)!
  ! ******************************************************************************

  USE netcdf
  USE dimphy
  USE indice_sol_mod
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, nbp_lev

  IMPLICIT NONE

  include "netcdf.inc"

  INTEGER kon, kev, zkon, zkev
!  PARAMETER (kon=iim*(jjm-1)+2, kev=llm)
  REAL :: phisfi(nbp_lon*(nbp_lat-2)+2) !phisfi(kon)
  REAL,ALLOCATABLE :: phisfi2(:,:) !phisfi2(nbp_lon, nbp_lat)
  REAL,ALLOCATABLE :: airefi2(:,:) !airefi2(nbp_lon, nbp_lat)

  REAL :: mfu(nbp_lon*(nbp_lat-2)+2,nbp_lev) !mfu(kon, kev)
  REAL :: mfd(nbp_lon*(nbp_lat-2)+2,nbp_lev) !mfd(kon, kev)
  REAL :: en_u(nbp_lon*(nbp_lat-2)+2,nbp_lev) !en_u(kon, kev)
  REAL :: de_u(nbp_lon*(nbp_lat-2)+2,nbp_lev) !de_u(kon, kev)
  REAL :: en_d(nbp_lon*(nbp_lat-2)+2,nbp_lev) !en_d(kon, kev)
  REAL :: de_d(nbp_lon*(nbp_lat-2)+2,nbp_lev) !de_d(kon, kev)
  REAL :: coefh(nbp_lon*(nbp_lat-2)+2,nbp_lev) !coefh(kon, kev)

  ! abd 25 11 02
  ! Thermiques
  REAL :: fm_therm(nbp_lon*(nbp_lat-2)+2,nbp_lev) !fm_therm(kon, kev)
  REAL :: en_therm(nbp_lon*(nbp_lat-2)+2,nbp_lev) !en_therm(kon, kev)
  REAL :: t(nbp_lon*(nbp_lat-2)+2,nbp_lev) !t(kon, kev)

  REAL,ALLOCATABLE :: mfu2(:,:,:) !mfu2(nbp_lon, nbp_lat, kev)
  REAL,ALLOCATABLE :: mfd2(:,:,:) !mfd2(nbp_lon, nbp_lat, kev)
  REAL,ALLOCATABLE :: en_u2(:,:,:) !en_u2(nbp_lon, nbp_lat, kev)
  REAL,ALLOCATABLE :: de_u2(:,:,:) !de_u2(nbp_lon, nbp_lat, kev)
  REAL,ALLOCATABLE :: en_d2(:,:,:) !en_d2(nbp_lon, nbp_lat, kev)
  REAL,ALLOCATABLE :: de_d2(:,:,:) !de_d2(nbp_lon, nbp_lat, kev)
  REAL,ALLOCATABLE :: coefh2(:,:,:) !coefh2(nbp_lon, nbp_lat, kev)
  REAL,ALLOCATABLE :: t2(:,:,:) !t2(nbp_lon, nbp_lat, kev)
  ! Thermiques
  REAL,ALLOCATABLE :: fm_therm2(:,:,:) !fm_therm2(nbp_lon, nbp_lat, kev)
  REAL,ALLOCATABLE :: en_therm2(:,:,:) !en_therm2(nbp_lon, nbp_lat, kev)

  REAL,ALLOCATABLE :: pl(:) !pl(kev)
  INTEGER irec
  INTEGER xid, yid, zid, tid
  INTEGER zrec, zim, zjm
  INTEGER ncrec, nckon, nckev, ncim, ncjm

  REAL :: airefi(nbp_lon*(nbp_lat-2)+2) !airefi(kon)
  CHARACTER *20 namedim

  ! !! attention !!
  ! attention il y a aussi le pb de def kon
  ! dim de phis??

  REAL :: frac_impa(nbp_lon*(nbp_lat-2)+2,nbp_lev) !frac_impa(kon, kev)
  REAL :: frac_nucl(nbp_lon*(nbp_lat-2)+2,nbp_lev) !frac_nucl(kon, kev)
  REAL,ALLOCATABLE :: frac_impa2(:,:,:) !frac_impa2(nbp_lon, nbp_lat, kev)
  REAL,ALLOCATABLE :: frac_nucl2(:,:,:) !frac_nucl2(nbp_lon, nbp_lat, kev)
  REAL :: pyu1(nbp_lon*(nbp_lat-2)+2) !pyu1(kon)
  REAL :: pyv1(nbp_lon*(nbp_lat-2)+2) !pyv1(kon)
  REAL,ALLOCATABLE :: pyu12(:,:), pyv12(:,:) !pyu12(nbp_lon, nbp_lat), pyv12(nbp_lon, nbp_lat)
  REAL :: ftsol(nbp_lon*(nbp_lat-2)+2,nbp_lev) !ftsol(kon, nbsrf)
  REAL :: psrf(nbp_lon*(nbp_lat-2)+2,nbp_lev) !psrf(kon, nbsrf)
  REAL,ALLOCATABLE :: ftsol1(:),ftsol2(:) !ftsol1(kon), ftsol2(kon)
  REAL,ALLOCATABLE :: ftsol3(:),ftsol4(:) !ftsol3(kon), ftsol4(kon)
  REAL,ALLOCATABLE :: psrf1(:), psrf2(:) !psrf1(kon), psrf2(kon)
  REAL,ALLOCATABLE :: psrf3(:), psrf4(:) !psrf3(kon), psrf4(kon)
  REAL,ALLOCATABLE :: ftsol12(:,:) !ftsol12(nbp_lon, nbp_lat)
  REAL,ALLOCATABLE :: ftsol22(:,:) !ftsol22(nbp_lon, nbp_lat)
  REAL,ALLOCATABLE :: ftsol32(:,:) !ftsol32(nbp_lon, nbp_lat)
  REAL,ALLOCATABLE :: ftsol42(:,:) !ftsol42(nbp_lon, nbp_lat)
  REAL,ALLOCATABLE :: psrf12(:,:) !psrf12(nbp_lon, nbp_lat)
  REAL,ALLOCATABLE :: psrf22(:,:) !psrf22(nbp_lon, nbp_lat)
  REAL,ALLOCATABLE :: psrf32(:,:) !psrf32(nbp_lon, nbp_lat)
  REAL,ALLOCATABLE :: psrf42(:,:) !psrf42(nbp_lon, nbp_lat)

  INTEGER,SAVE :: ncidp
  INTEGER,SAVE :: varidmfu, varidmfd, varidps, varidenu, variddeu
  INTEGER,SAVE :: varidt
  INTEGER,SAVE :: varidend, varidded, varidch, varidfi, varidfn
  ! therm
  INTEGER,SAVE :: varidfmth, varidenth
  INTEGER,SAVE :: varidyu1, varidyv1, varidpl, varidai, varididvt
  INTEGER,SAVE :: varidfts1, varidfts2, varidfts3, varidfts4
  INTEGER,SAVE :: varidpsr1, varidpsr2, varidpsr3, varidpsr4

  INTEGER l, i
  INTEGER start(4), count(4), status
  REAL rcode
  LOGICAL,SAVE :: first=.TRUE.

  ! Allocate arrays
  kon=nbp_lon*(nbp_lat-2)+2
  kev=nbp_lev

  ALLOCATE(phisfi2(nbp_lon, nbp_lat))
  ALLOCATE(airefi2(nbp_lon, nbp_lat))
  ALLOCATE(mfu2(nbp_lon, nbp_lat, kev))
  ALLOCATE(mfd2(nbp_lon, nbp_lat, kev))
  ALLOCATE(en_u2(nbp_lon, nbp_lat, kev))
  ALLOCATE(de_u2(nbp_lon, nbp_lat, kev))
  ALLOCATE(en_d2(nbp_lon, nbp_lat, kev))
  ALLOCATE(de_d2(nbp_lon, nbp_lat, kev))
  ALLOCATE(coefh2(nbp_lon, nbp_lat, kev))
  ALLOCATE(t2(nbp_lon, nbp_lat, kev))
  ALLOCATE(fm_therm2(nbp_lon, nbp_lat, kev))
  ALLOCATE(en_therm2(nbp_lon, nbp_lat, kev))
  ALLOCATE(pl(kev))
  ALLOCATE(frac_impa2(nbp_lon, nbp_lat, kev))
  ALLOCATE(frac_nucl2(nbp_lon, nbp_lat, kev))
  ALLOCATE(pyu12(nbp_lon, nbp_lat), pyv12(nbp_lon, nbp_lat))
  ALLOCATE(ftsol1(kon), ftsol2(kon))
  ALLOCATE(ftsol3(kon), ftsol4(kon))
  ALLOCATE(psrf1(kon), psrf2(kon))
  ALLOCATE(psrf3(kon), psrf4(kon))
  ALLOCATE(ftsol12(nbp_lon, nbp_lat))
  ALLOCATE(ftsol22(nbp_lon, nbp_lat))
  ALLOCATE(ftsol32(nbp_lon, nbp_lat))
  ALLOCATE(ftsol42(nbp_lon, nbp_lat))
  ALLOCATE(psrf12(nbp_lon, nbp_lat))
  ALLOCATE(psrf22(nbp_lon, nbp_lat))
  ALLOCATE(psrf32(nbp_lon, nbp_lat))
  ALLOCATE(psrf42(nbp_lon, nbp_lat))

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

    ! Thermiques
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
    status = nf_inq_dim(ncidp, zid, namedim, nckev)
    status = nf_inq_dim(ncidp, tid, namedim, ncrec)

    zrec = ncrec
    zkev = nckev
    zim = ncim
    zjm = ncjm

    zkon = zim*(zjm-2) + 2

    WRITE (*, *) 'read_pstoke : zrec = ', zrec
    WRITE (*, *) 'read_pstoke : kev = ', zkev
    WRITE (*, *) 'read_pstoke : zim = ', zim
    WRITE (*, *) 'read_pstoke : zjm = ', zjm
    WRITE (*, *) 'read_pstoke : kon = ', zkon

    ! niveaux de pression

    status = nf_get_vara_real(ncidp, varidpl, 1, kev, pl)

    ! lecture de aire et phis

    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = 0

    count(1) = zim
    count(2) = zjm
    count(3) = 1
    count(4) = 0


    ! **** Geopotentiel au sol ***************************************
    ! phis
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidps, start, count, phisfi2)
#else
    status = nf_get_vara_real(ncidp, varidps, start, count, phisfi2)
#endif
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, phisfi2, phisfi)

    ! **** Aires des mails aux sol ************************************
    ! aire
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidai, start, count, airefi2)
#else
    status = nf_get_vara_real(ncidp, varidai, start, count, airefi2)
#endif
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, airefi2, airefi)
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
    count(3) = kev
    count(4) = 1

    ! **** Temperature ********************************************
    ! A FAIRE : Es-ce necessaire ?

    ! abder t
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidt, start, count, t2)
#else
    status = nf_get_vara_real(ncidp, varidt, start, count, t2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, t2, t)

    ! **** Flux pour la convection (Tiedtk)
    ! ********************************************
    ! mfu
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidmfu, start, count, mfu2)
#else
    status = nf_get_vara_real(ncidp, varidmfu, start, count, mfu2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, mfu2, mfu)

    ! mfd
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidmfd, start, count, mfd2)
#else
    status = nf_get_vara_real(ncidp, varidmfd, start, count, mfd2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, mfd2, mfd)

    ! en_u
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidenu, start, count, en_u2)
#else
    status = nf_get_vara_real(ncidp, varidenu, start, count, en_u2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, en_u2, en_u)

    ! de_u
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, variddeu, start, count, de_u2)
#else
    status = nf_get_vara_real(ncidp, variddeu, start, count, de_u2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, de_u2, de_u)

    ! en_d
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidend, start, count, en_d2)
#else
    status = nf_get_vara_real(ncidp, varidend, start, count, en_d2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, en_d2, en_d)

    ! de_d
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidded, start, count, de_d2)
#else
    status = nf_get_vara_real(ncidp, varidded, start, count, de_d2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, de_d2, de_d)

    ! **** Coefficient de mellange turbulent
    ! *******************************************
    ! coefh
    PRINT *, 'LECTURE de coefh a irec =', irec
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidch, start, count, coefh2)
#else
    status = nf_get_vara_real(ncidp, varidch, start, count, coefh2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, coefh2, coefh)
    ! call dump2d(iip1,jjp1,coefh2(1,2),'COEFH2READ   ')
    ! call dump2d(iim ,jjm ,coefh (2,2),'COEFH2READ   ')

    ! **** Flux ascendants et entrant dans le thermique
    ! **********************************
    ! Thermiques
    PRINT *, 'LECTURE de fm_therm a irec =', irec
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfmth, start, count, fm_therm2)
#else
    status = nf_get_vara_real(ncidp, varidfmth, start, count, fm_therm2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, fm_therm2, fm_therm)
    PRINT *, 'LECTURE de en_therm a irec =', irec
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidenth, start, count, en_therm2)
#else
    status = nf_get_vara_real(ncidp, varidenth, start, count, en_therm2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, en_therm2, en_therm)

    ! **** Coefficients de lessivage
    ! *******************************************
    ! frac_impa
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfi, start, count, frac_impa2)
#else
    status = nf_get_vara_real(ncidp, varidfi, start, count, frac_impa2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, frac_impa2, frac_impa)

    ! frac_nucl

#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfn, start, count, frac_nucl2)
#else
    status = nf_get_vara_real(ncidp, varidfn, start, count, frac_nucl2)
#endif
    CALL gr_ecrit_fi(kev, kon, nbp_lon, nbp_lat, frac_nucl2, frac_nucl)

    ! **** Vents aux sol ********************************************

    start(3) = irec
    start(4) = 0
    count(3) = 1
    count(4) = 0

    ! pyu1
    PRINT *, 'LECTURE de yu1 a irec =', irec
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidyu1, start, count, pyu12)
#else
    status = nf_get_vara_real(ncidp, varidyu1, start, count, pyu12)
#endif
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, pyu12, pyu1)

    ! pyv1
    PRINT *, 'LECTURE de yv1 a irec =', irec
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidyv1, start, count, pyv12)
#else
    status = nf_get_vara_real(ncidp, varidyv1, start, count, pyv12)
#endif
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, pyv12, pyv1)

    ! **** Temerature au sol ********************************************
    ! ftsol1
    PRINT *, 'LECTURE de ftsol1 a irec =', irec
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfts1, start, count, ftsol12)
#else
    status = nf_get_vara_real(ncidp, varidfts1, start, count, ftsol12)
#endif
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, ftsol12, ftsol1)

    ! ftsol2
    PRINT *, 'LECTURE de ftsol2 a irec =', irec
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfts2, start, count, ftsol22)
#else
    status = nf_get_vara_real(ncidp, varidfts2, start, count, ftsol22)
#endif
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, ftsol22, ftsol2)

    ! ftsol3
    PRINT *, 'LECTURE de ftsol3 a irec =', irec
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfts3, start, count, ftsol32)
#else
    status = nf_get_vara_real(ncidp, varidfts3, start, count, ftsol32)
#endif
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, ftsol32, ftsol3)

    ! ftsol4
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidfts4, start, count, ftsol42)
#else
    status = nf_get_vara_real(ncidp, varidfts4, start, count, ftsol42)
#endif
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, ftsol42, ftsol4)

    ! **** Nature sol ********************************************
    ! psrf1
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidpsr1, start, count, psrf12)
#else
    status = nf_get_vara_real(ncidp, varidpsr1, start, count, psrf12)
#endif
    ! call dump2d(iip1-1,jjm+1,psrf12,'PSRF1NC')
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, psrf12, psrf1)

    ! psrf2
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidpsr2, start, count, psrf22)
#else
    status = nf_get_vara_real(ncidp, varidpsr2, start, count, psrf22)
#endif
    ! call dump2d(iip1-1,jjm+1,psrf22,'PSRF2NC')
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, psrf22, psrf2)

    ! psrf3
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidpsr3, start, count, psrf32)
#else
    status = nf_get_vara_real(ncidp, varidpsr3, start, count, psrf32)
#endif
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, psrf32, psrf3)

    ! psrf4
#ifdef NC_DOUBLE
    status = nf_get_vara_double(ncidp, varidpsr4, start, count, psrf42)
#else
    status = nf_get_vara_real(ncidp, varidpsr4, start, count, psrf42)
#endif
    CALL gr_ecrit_fi(1, kon, nbp_lon, nbp_lat, psrf42, psrf4)

    DO i = 1, kon

      psrf(i, 1) = psrf1(i)
      psrf(i, 2) = psrf2(i)
      psrf(i, 3) = psrf3(i)
      ! test abderr
      ! print*,'Dans read_pstoke psrf3 =',psrf3(i),i
      psrf(i, 4) = psrf4(i)

      ftsol(i, 1) = ftsol1(i)
      ftsol(i, 2) = ftsol2(i)
      ftsol(i, 3) = ftsol3(i)
      ftsol(i, 4) = ftsol4(i)

    END DO

  END IF

  RETURN

END SUBROUTINE read_pstoke0

