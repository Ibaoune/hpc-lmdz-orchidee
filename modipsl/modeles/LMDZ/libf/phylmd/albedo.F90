! $Id: albedo.F90 2413 2015-12-18 19:27:39Z oboucher $
module albedo

  IMPLICIT NONE

contains

  SUBROUTINE alboc(rjour, rlat, albedo)
    USE dimphy
    ! ======================================================================
    ! Auteur(s): Z.X. Li (LMD/CNRS) (adaptation du GCM du LMD)
    ! Date: le 16 mars 1995
    ! Objet: Calculer l'albedo sur l'ocean
    ! Methode: Integrer numeriquement l'albedo pendant une journee

    ! Arguments;
    ! rjour (in,R)  : jour dans l'annee (a compter du 1 janvier)
    ! rlat (in,R)   : latitude en degre
    ! albedo (out,R): albedo obtenu (de 0 a 1)
    ! ======================================================================
    include "YOMCST.h"
    include "clesphys.h"

    INTEGER npts ! il controle la precision de l'integration
    PARAMETER (npts=120) ! 120 correspond a l'interval 6 minutes

    REAL rlat(klon), rjour, albedo(klon)
    REAL zdist, zlonsun, zpi, zdeclin
    REAL rmu, alb, srmu, salb, fauxo, aa, bb
    INTEGER i, k
    ! ccIM
    LOGICAL ancien_albedo
    PARAMETER (ancien_albedo=.FALSE.)
    ! SAVE albedo

    IF (ancien_albedo) THEN

       zpi = 4.*atan(1.)

       ! Calculer la longitude vraie de l'orbite terrestre:
       CALL orbite(rjour, zlonsun, zdist)

       ! Calculer la declinaison du soleil (qui varie entre + et - R_incl):
       zdeclin = asin(sin(zlonsun*zpi/180.0)*sin(r_incl*zpi/180.0))

       DO i = 1, klon
          aa = sin(rlat(i)*zpi/180.0)*sin(zdeclin)
          bb = cos(rlat(i)*zpi/180.0)*cos(zdeclin)

          ! Midi local (angle du temps = 0.0):
          rmu = aa + bb*cos(0.0)
          rmu = max(0.0, rmu)
          fauxo = (1.47-acos(rmu))/.15
          alb = 0.03 + 0.630/(1.+fauxo*fauxo)
          srmu = rmu
          salb = alb*rmu

          ! Faire l'integration numerique de midi a minuit (le facteur 2
          ! prend en compte l'autre moitie de la journee):
          DO k = 1, npts
             rmu = aa + bb*cos(real(k)/real(npts)*zpi)
             rmu = max(0.0, rmu)
             fauxo = (1.47-acos(rmu))/.15
             alb = 0.03 + 0.630/(1.+fauxo*fauxo)
             srmu = srmu + rmu*2.0
             salb = salb + alb*rmu*2.0
          END DO
          IF (srmu/=0.0) THEN
             albedo(i) = salb/srmu
          ELSE ! nuit polaire (on peut prendre une valeur quelconque)
             albedo(i) = 1.0
          END IF
       END DO

       ! nouvel albedo

    ELSE

       zpi = 4.*atan(1.)

       ! Calculer la longitude vraie de l'orbite terrestre:
       CALL orbite(rjour, zlonsun, zdist)

       ! Calculer la declinaison du soleil (qui varie entre + et - R_incl):
       zdeclin = asin(sin(zlonsun*zpi/180.0)*sin(r_incl*zpi/180.0))

       DO i = 1, klon
          aa = sin(rlat(i)*zpi/180.0)*sin(zdeclin)
          bb = cos(rlat(i)*zpi/180.0)*cos(zdeclin)

          ! Midi local (angle du temps = 0.0):
          rmu = aa + bb*cos(0.0)
          rmu = max(0.0, rmu)
          ! IM cf. PB  alb = 0.058/(rmu + 0.30)
          ! alb = 0.058/(rmu + 0.30) * 1.5
          alb = 0.058/(rmu+0.30)*1.2
          ! alb = 0.058/(rmu + 0.30) * 1.3
          srmu = rmu
          salb = alb*rmu

          ! Faire l'integration numerique de midi a minuit (le facteur 2
          ! prend en compte l'autre moitie de la journee):
          DO k = 1, npts
             rmu = aa + bb*cos(real(k)/real(npts)*zpi)
             rmu = max(0.0, rmu)
             ! IM cf. PB      alb = 0.058/(rmu + 0.30)
             ! alb = 0.058/(rmu + 0.30) * 1.5
             alb = 0.058/(rmu+0.30)*1.2
             ! alb = 0.058/(rmu + 0.30) * 1.3
             srmu = srmu + rmu*2.0
             salb = salb + alb*rmu*2.0
          END DO
          IF (srmu/=0.0) THEN
             albedo(i) = salb/srmu
          ELSE ! nuit polaire (on peut prendre une valeur quelconque)
             albedo(i) = 1.0
          END IF
       END DO
    END IF
    RETURN
  END SUBROUTINE alboc
  ! =====================================================================
  SUBROUTINE alboc_cd(rmu0, albedo)
    USE dimphy

    ! Auteur(s): Z.X. Li (LMD/CNRS)
    ! date: 19940624
    ! Calculer l'albedo sur l'ocean en fonction de l'angle zenithal moyen
    ! Formule due a Larson and Barkstrom (1977) Proc. of the symposium
    ! on radiation in the atmosphere, 19-28 August 1976, science Press,
    ! 1977 pp 451-453, ou These de 3eme cycle de Sylvie Joussaume.

    ! Arguments
    ! rmu0    (in): cosinus de l'angle solaire zenithal
    ! albedo (out): albedo de surface de l'ocean
    ! ======================================================================
    include "clesphys.h"
    REAL, intent(in):: rmu0(klon)
    real, intent(out):: albedo(klon)

    REAL fauxo
    INTEGER i
    LOGICAL ancien_albedo
    PARAMETER (ancien_albedo=.FALSE.)

    IF (ancien_albedo) THEN
       DO i = 1, klon
          fauxo = (1.47-acos(max(rmu0(i), 0.0)))/0.15
          albedo(i) = 0.03+.630/(1.+fauxo*fauxo)
          albedo(i) = max(min(albedo(i),0.60), 0.04)
       END DO
    ELSE
       DO i = 1, klon
          albedo(i) = 0.058/(max(rmu0(i), 0.0)+0.30)
          albedo(i) = max(min(albedo(i),0.60), 0.04)
       END DO
    END IF

  END SUBROUTINE alboc_cd

end module albedo
