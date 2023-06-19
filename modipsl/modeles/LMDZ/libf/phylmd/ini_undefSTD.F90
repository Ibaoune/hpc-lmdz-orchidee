
! $Id: ini_undefSTD.F90 3435 2019-01-22 15:21:59Z fairhead $


SUBROUTINE ini_undefstd(itap, itapm1)
  USE dimphy
  USE phys_state_var_mod ! Variables sauvegardees de la physique
  USE phys_cal_mod, ONLY: mth_len
  IMPLICIT NONE
  include "clesphys.h"

  ! ====================================================================

  ! I. Musat : 09.2004

  ! Initialisation - a des frequences differentes :

  ! 1) des variables moyennees sur la journee "day" ou sur le mois "mth"
  ! calculees a partir des valeurs "instantannees" de la physique

  ! 2) des variables moyennes mensuelles "NMC" calculees a partir des val.
  ! toutes les 6 heures

  ! nout=1 !var. journaliere "day" moyenne sur tous les pas de temps
  ! ! de la physique
  ! nout=2 !var. mensuelle "mth" moyennee sur tous les pas de temps
  ! ! de la physique
  ! nout=3 !var. mensuelle "NMC" moyennee toutes les 6heures

  ! NB: mettre "inst(X)" dans le write_hist*NMC.h !
  ! ====================================================================

  ! ym #include "dimphy.h"
  ! variables Input/Output
  ! INTEGER nlevSTD, klevSTD, itap
  INTEGER itap, itapm1
  ! PARAMETER(klevSTD=17)
  ! REAL dtime

  ! variables locales
  ! INTEGER i, k, nout, n
  INTEGER i, k, n
  ! PARAMETER(nout=3) !nout=1 day/nout=2 mth/nout=3 NMC
  ! REAL freq_outNMC(nout)
  REAL un_jour
  PARAMETER (un_jour=86400.)

  ! variables Output
  ! REAL tnondef(klon,klevSTD,nout)
  ! REAL tsumSTD(klon,klevSTD,nout)

  DO n = 1, nout

    ! initialisation variables en debut du mois, de la journee ou des 6h

    IF (n==1 .AND. itap-itapm1==1 .OR. n>1 .AND. mod(itap,nint( &
        freq_outnmc(n)/phys_tstep))==1) THEN
      ! print*,'ini_undefSTD n itap',n,itap
      DO k = 1, nlevstd
        DO i = 1, klon
          tnondef(i, k, n) = 0.
          tsumstd(i, k, n) = 0.
          usumstd(i, k, n) = 0.
          vsumstd(i, k, n) = 0.
          wsumstd(i, k, n) = 0.
          phisumstd(i, k, n) = 0.
          qsumstd(i, k, n) = 0.
          rhsumstd(i, k, n) = 0.
          uvsumstd(i, k, n) = 0.
          vqsumstd(i, k, n) = 0.
          vtsumstd(i, k, n) = 0.
          wqsumstd(i, k, n) = 0.
          vphisumstd(i, k, n) = 0.
          wtsumstd(i, k, n) = 0.
          u2sumstd(i, k, n) = 0.
          v2sumstd(i, k, n) = 0.
          t2sumstd(i, k, n) = 0.
          o3sumstd(i, k, n) = 0.
          o3daysumstd(i, k, n) = 0.
        END DO !i
      END DO !k
    END IF !
  END DO !n
  RETURN
END SUBROUTINE ini_undefstd
