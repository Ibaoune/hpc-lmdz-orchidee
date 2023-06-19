
! $Header$

SUBROUTINE moyglo_pondaire(nhori, champ, aire, ok_msk, msk, moyglo)

  USE dimphy
  IMPLICIT NONE

  ! ==================================================================
  ! I. Musat, 07.2004

  ! Calcul moyenne globale ponderee par l'aire totale, avec ou sans masque

  ! moyenne = Somme_(champ* aire)/Somme_aire

  ! ==================================================================

  INTEGER i, nhori
  REAL champ(klon), aire(klon), msk(klon)
  LOGICAL ok_msk
  REAL moyglo

  ! var locale
  REAL airetot

  ! PRINT*,'moyglo_pondaire nhori',nhori

  airetot = 0.
  moyglo = 0.

  IF (ok_msk) THEN
    DO i = 1, nhori
      ! IF(msk(i).EQ.1.) THEN
      IF (msk(i)>0.) THEN

        ! aire totale
        airetot = airetot + aire(i)*msk(i)

        ! ponderation par la masse
        moyglo = moyglo + champ(i)*aire(i)*msk(i)
      END IF
    END DO

  ELSE !ok_msk
    DO i = 1, nhori

      ! aire totale
      airetot = airetot + aire(i)

      ! ponderation par la masse
      moyglo = moyglo + champ(i)*aire(i)
    END DO

  END IF

  ! moyenne ponderee par l'aire
  moyglo = moyglo/airetot

  RETURN
END SUBROUTINE moyglo_pondaire

SUBROUTINE moyglo_pondaima(nhori, nvert, champ, aire, pbord, moyglo)
  USE dimphy
  IMPLICIT NONE
  ! ==================================================================
  ! I. Musat, 07.2004

  ! Calcul moyenne globale ponderee par la masse d'air, divisee par l'aire
  ! totale avec ou sans masque

  ! moyenne = Somme_(champ* masse_dair)/Somme_aire

  ! ==================================================================
  include "YOMCST.h"
  INTEGER i, k, nhori, nvert
  REAL champ(klon, klev), aire(klon)
  REAL pbord(klon, klev+1)
  REAL moyglo

  ! var locale
  REAL airetot

  ! PRINT*,'moyglo_pondaima RG, nhori, nvert',RG,nhori,nvert

  ! ponderation par la masse
  moyglo = 0.
  DO k = 1, nvert
    DO i = 1, nhori
      moyglo = moyglo + champ(i, k)*(pbord(i,k)-pbord(i,k+1))/rg*aire(i)
    END DO
  END DO

  ! aire totale
  airetot = 0.
  DO i = 1, nhori
    airetot = airetot + aire(i)
  END DO

  ! moyenne par mettre carre avec ponderation par la masse
  moyglo = moyglo/airetot

  RETURN
END SUBROUTINE moyglo_pondaima

SUBROUTINE moyglo_pondmass(nhori, nvert, champ, aire, pbord, moyglo)
  USE dimphy
  IMPLICIT NONE
  ! ==================================================================
  ! I. Musat, 07.2004

  ! Calcul moyenne globale ponderee par la masse d'air, divisee par la
  ! masse totale d'air, avec ou sans masque

  ! moyenne = Somme_(champ* masse_dair)/Somme_(masse_dair)

  ! ==================================================================
  include "YOMCST.h"
  INTEGER i, k, nhori, nvert
  REAL champ(klon, klev), aire(klon)
  REAL pbord(klon, klev+1)
  REAL moyglo

  ! var locale
  REAL massetot

  ! PRINT*,'moyglo_pondmass RG, nhori, nvert',RG,nhori,nvert

  ! ponderation par la masse
  moyglo = 0.
  DO k = 1, nvert
    DO i = 1, nhori
      moyglo = moyglo + champ(i, k)*(pbord(i,k)-pbord(i,k+1))/rg*aire(i)
    END DO
  END DO

  ! masse totale
  massetot = 0.
  DO k = 1, nvert
    DO i = 1, nhori
      massetot = massetot + (pbord(i,k)-pbord(i,k+1))/rg*aire(i)
    END DO
  END DO

  ! moyenne par mettre carre avec ponderation par la masse
  moyglo = moyglo/massetot

  RETURN
END SUBROUTINE moyglo_pondmass

