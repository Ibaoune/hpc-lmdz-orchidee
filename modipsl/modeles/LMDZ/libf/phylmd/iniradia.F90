SUBROUTINE iniradia(klon, klev, pres)

  IMPLICIT NONE
  ! ======================================================================

  ! Auteur(s) MP Lefebvre        date: 20080827

  ! Objet: initialise le rayonnement RRTM
  ! ======================================================================
  ! Arguments:

  ! klon----input-I-nombre de points horizontaux
  ! klev----input-I-nombre de couches verticales
  ! pres----input-R-pression pour chaque inter-couche (en Pa)
  ! ======================================================================

  INTEGER klon
  INTEGER klev
  REAL pres(klev+1)

  include "clesphys.h"

  ! CALL suphel     ! initialiser constantes et parametres phys.
  ! print*,'Physiq: apres suphel '
#if CPP_RRTM
  if (iflag_rrtm .eq. 1) then 
     CALL suinit(klon, klev)
     PRINT *, 'iniradia: apres suinit '
     ! calcul des niveaux de pression de reference au bord des couches pour
     ! l'intialisation des aerosols. Momentannement, on passe un point de
     ! grille du profil de pression.
     CALL surayolmd          ! initialiser le rayonnement RRTM
     PRINT *, 'iniradia: apres surayolmd '
  endif
#endif

  RETURN
END SUBROUTINE iniradia
