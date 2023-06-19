! $Header$

SUBROUTINE condsurf(jour, jourvrai, lmt_bils)
  USE dimphy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  USE indice_sol_mod
  USE time_phylmdz_mod, ONLY: annee_ref
  IMPLICIT NONE

  ! I. Musat 05.2005

  ! Lire chaque jour le bilan de chaleur au sol issu
  ! d'un run atmospherique afin de l'utiliser dans
  ! dans un run "slab" ocean
  ! -----------------------------------------
  ! jour     : input  , numero du jour a lire
  ! jourvrai : input  , vrai jour de la simulation

  ! lmt_bils: bilan chaleur au sol (a utiliser pour "slab-ocean")

  include "netcdf.inc"
  INTEGER nid, nvarid
  INTEGER debut(2)
  INTEGER epais(2)

  include "clesphys.h"

  INTEGER nannemax
  PARAMETER (nannemax=60)

  INTEGER jour, jourvrai
  REAL lmt_bils(klon) !bilan chaleur au sol

  ! Variables locales:
  INTEGER ig, i, kt, ierr
  LOGICAL ok
  INTEGER anneelim, anneemax
  CHARACTER *20 fich

  REAL :: lmt_bils_glo(klon_glo)

  ! c
  ! c   .....................................................................
  ! c
  ! c    Pour lire le fichier limit correspondant vraiment  a l'annee de la
  ! c     simulation en cours , il suffit de mettre  ok_limitvrai = .TRUE.
  ! c
  ! c
  ! ......................................................................



  IF (jour<0 .OR. jour>(360-1)) THEN
    PRINT *, 'Le jour demande n est pas correct: ', jour
    CALL abort_physic('condsurf', '', 1)
  END IF

  anneelim = annee_ref
  anneemax = annee_ref + nannemax


  IF (ok_limitvrai) THEN
    DO kt = 1, nannemax
      IF (jourvrai<=(kt-1)*360+359) THEN
        WRITE (fich, '("limit",i4,".nc")') anneelim
        ! PRINT *,' Fichier  Limite ',fich
        GO TO 100
      END IF
      anneelim = anneelim + 1
    END DO

    PRINT *, ' PBS ! Le jour a lire sur le fichier limit ne se '
    PRINT *, ' trouve pas sur les ', nannemax, ' annees a partir de '
    PRINT *, ' l annee de debut', annee_ref
    CALL abort_physic('condsurf', '', 1)

100 CONTINUE

  ELSE

    WRITE (fich, '("limitNEW.nc")')
    ! PRINT *,' Fichier  Limite ',fich
  END IF

  ! Ouvrir le fichier en format NetCDF:

  !$OMP MASTER
  IF (is_mpi_root) THEN
    ierr = nf_open(fich, nf_nowrite, nid)
    IF (ierr/=nf_noerr) THEN
      WRITE (6, *) ' Pb d''ouverture du fichier ', fich
      WRITE (6, *) ' Le fichier limit ', fich, ' (avec 4 chiffres , pour'
      WRITE (6, *) '       l an 2000 )  ,  n existe  pas !  '
      WRITE (6, *) ' ierr = ', ierr
      CALL abort_physic('condsurf', '', 1)
    END IF
    ! DO k = 1, jour
    ! La tranche de donnees a lire:

    debut(1) = 1
    debut(2) = jourvrai
    epais(1) = klon_glo
    epais(2) = 1
    ! Bilan flux de chaleur au sol:

    ierr = nf_inq_varid(nid, 'BILS', nvarid)
    IF (ierr/=nf_noerr) THEN
      CALL abort_physic('cond_surf', 'Le champ <BILS> est absent', 1)
    END IF
    PRINT *, 'debut,epais', debut, epais, 'jour,jourvrai', jour, jourvrai
#ifdef NC_DOUBLE
    ierr = nf_get_vara_double(nid, nvarid, debut, epais, lmt_bils_glo)
#else
    ierr = nf_get_vara_real(nid, nvarid, debut, epais, lmt_bils_glo)
#endif
    IF (ierr/=nf_noerr) THEN
      CALL abort_physic('condsurf', 'Lecture echouee pour <BILS>', 1)
    END IF
    ! ENDDO !k = 1, jour

    ! Fermer le fichier:

    ierr = nf_close(nid)

  END IF ! is_mpi_root==0

  !$OMP END MASTER
  CALL scatter(lmt_bils_glo, lmt_bils)



  ! PRINT*, 'lmt_bils est lu pour jour: ', jour

  RETURN
END SUBROUTINE condsurf
