
! $Id: undefSTD.F90 3435 2019-01-22 15:21:59Z fairhead $

SUBROUTINE undefstd(itap, read_climoz)
  USE netcdf
  USE dimphy
#ifdef CPP_IOIPSL 
  USE phys_state_var_mod
#endif
#ifdef CPP_XIOS
  USE wxios, ONLY: missing_val
#endif

  IMPLICIT NONE
  include "clesphys.h"
#ifndef CPP_XIOS 
  REAL :: missing_val
#endif

  ! ====================================================================

  ! I. Musat : 09.2004

  ! Calcul * du nombre de pas de temps (FLOAT(ecrit_XXX)-tnondef))
  ! ou la variable tlevSTD est bien definie (.NE.missing_val),
  ! et
  ! * de la somme de tlevSTD => tsumSTD

  ! nout=1 !var. journaliere "day" moyenne sur tous les pas de temps
  ! ! de la physique
  ! nout=2 !var. mensuelle "mth" moyennee sur tous les pas de temps
  ! ! de la physique
  ! nout=3 !var. mensuelle "NMC" moyennee toutes les ecrit_hf


  ! NB: mettre "inst(X)" dans le write_hist*NMC.h !
  ! ====================================================================

  ! ym#include "dimphy.h"
  ! variables Input

  ! INTEGER nlevSTD, klevSTD, itap
  ! PARAMETER(klevSTD=17)
  INTEGER itap

  ! variables locales
  ! INTEGER i, k, nout, n
  ! PARAMETER(nout=3) !nout=1 : day; =2 : mth; =3 : NMC
  INTEGER i, k, n
  ! REAL freq_calNMC(nout)
  INTEGER read_climoz

  ! variables Output
  ! REAL tlevSTD(klon,klevSTD), tsumSTD(klon,klevSTD,nout)
  ! LOGICAL oknondef(klon,klevSTD,nout)
  ! REAL tnondef(klon,klevSTD,nout)

! REAL missing_val

! missing_val = nf90_fill_real
#ifndef CPP_XIOS
      missing_val=missing_val_nf90
#endif

  DO n = 1, nout


    ! calcul variables tous les freq_calNMC(n)/phys_tstep pas de temps
    ! de la physique

    IF (mod(itap,nint(freq_calnmc(n)/phys_tstep))==0) THEN
      DO k = 1, nlevstd
        DO i = 1, klon
          IF (tlevstd(i,k)==missing_val) THEN
            ! IF(oknondef(i,k,n)) THEN
            tnondef(i, k, n) = tnondef(i, k, n) + 1.
            ! ENDIF !oknondef(i,k)

          ELSE IF (tlevstd(i,k)/=missing_val) THEN
            tsumstd(i, k, n) = tsumstd(i, k, n) + tlevstd(i, k)
            usumstd(i, k, n) = usumstd(i, k, n) + ulevstd(i, k)
            vsumstd(i, k, n) = vsumstd(i, k, n) + vlevstd(i, k)
            wsumstd(i, k, n) = wsumstd(i, k, n) + wlevstd(i, k)
            phisumstd(i, k, n) = phisumstd(i, k, n) + philevstd(i, k)
            qsumstd(i, k, n) = qsumstd(i, k, n) + qlevstd(i, k)
            rhsumstd(i, k, n) = rhsumstd(i, k, n) + rhlevstd(i, k)
            uvsumstd(i, k, n) = uvsumstd(i, k, n) + uvstd(i, k)
            vqsumstd(i, k, n) = vqsumstd(i, k, n) + vqstd(i, k)
            vtsumstd(i, k, n) = vtsumstd(i, k, n) + vtstd(i, k)
            wqsumstd(i, k, n) = wqsumstd(i, k, n) + wqstd(i, k)
            vphisumstd(i, k, n) = vphisumstd(i, k, n) + vphistd(i, k)
            wtsumstd(i, k, n) = wtsumstd(i, k, n) + wtstd(i, k)
            u2sumstd(i, k, n) = u2sumstd(i, k, n) + u2std(i, k)
            v2sumstd(i, k, n) = v2sumstd(i, k, n) + v2std(i, k)
            t2sumstd(i, k, n) = t2sumstd(i, k, n) + t2std(i, k)
            o3sumstd(i, k, n) = o3sumstd(i, k, n) + o3std(i, k)
            IF (read_climoz==2) o3daysumstd(i, k, n) = o3daysumstd(i, k, n) + &
              o3daystd(i, k)

          END IF
        END DO !i
      END DO !k

    END IF !MOD(itap,NINT(freq_calNMC(n)/phys_tstep)).EQ.0

  END DO !n

  RETURN
END SUBROUTINE undefstd
