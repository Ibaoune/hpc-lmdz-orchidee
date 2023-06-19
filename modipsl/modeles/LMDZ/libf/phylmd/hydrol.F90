
! $Header$



SUBROUTINE hydrol(dtime, pctsrf, rain_fall, snow_fall, evap, agesno, tsol, &
    qsol, snow, runoff)
  USE dimphy
  USE indice_sol_mod

  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS)
  ! date: 19940414
  ! ======================================================================

  ! Traitement de l'hydrologie du sol
  ! ---------------------------------
  ! rain_fall: taux de pluie
  ! snow_fall: taux de neige
  ! agesno: age de la neige
  ! evap: taux d'evaporation
  ! tsol: temperature du sol
  ! qsol: humidite du sol
  ! snow: couverture neigeuse

  include "YOMCST.h"

  REAL chasno ! epaisseur du sol: 0.15 m
  PARAMETER (chasno=3.334E+05/(2.3867E+06*0.15))
  REAL mx_eau_sol
  PARAMETER (mx_eau_sol=150.0)

  REAL dtime
  REAL pctsrf(klon, nbsrf)
  REAL snow(klon, nbsrf), tsol(klon, nbsrf), qsol(klon, nbsrf)
  REAL snow_fall(klon), rain_fall(klon), evap(klon)
  REAL runoff(klon), agesno(klon)

  INTEGER i, is
  REAL subli, fsno
  ! -----------------------------------------------------------------------
  DO i = 1, klon

    runoff(i) = 0.0

    is = is_ter
    snow(i, is) = snow(i, is) + snow_fall(i)*dtime*pctsrf(i, is)
    IF (pctsrf(i,is)>epsfra) THEN
      subli = min(evap(i)*dtime, snow(i,is))
      snow(i, is) = snow(i, is) - subli
      fsno = min(max((tsol(i,is)-rtt)/chasno,0.0), snow(i,is))
      snow(i, is) = snow(i, is) - fsno
      tsol(i, is) = tsol(i, is) - fsno*chasno
      qsol(i, is) = qsol(i, is) + (rain_fall(i)-evap(i))*dtime + subli + fsno
      qsol(i, is) = max(qsol(i,is), 0.0)
      runoff(i) = runoff(i) + max(qsol(i,is)-mx_eau_sol, 0.0)*pctsrf(i, is)
      qsol(i, is) = min(qsol(i,is), mx_eau_sol)
      ! cc         ELSE
      ! cc            snow(i,is) = 0.0
      ! cc            qsol(i,is) = 0.0
      ! cc            tsol(i,is) = 0.0
    END IF

    is = is_lic
    snow(i, is) = snow(i, is) + snow_fall(i)*dtime*pctsrf(i, is)
    IF (pctsrf(i,is)>epsfra) THEN
      subli = min(evap(i)*dtime, snow(i,is))
      snow(i, is) = snow(i, is) - subli
      fsno = min(max((tsol(i,is)-rtt)/chasno,0.0), snow(i,is))
      snow(i, is) = snow(i, is) - fsno
      tsol(i, is) = tsol(i, is) - fsno*chasno
      qsol(i, is) = qsol(i, is) + (rain_fall(i)-evap(i))*dtime + subli + fsno
      qsol(i, is) = max(qsol(i,is), 0.0)
      runoff(i) = runoff(i) + max(qsol(i,is)-mx_eau_sol, 0.0)*pctsrf(i, is)
      qsol(i, is) = min(qsol(i,is), mx_eau_sol)
      ! je limite la temperature a RTT-1.8 (il faudrait aussi prendre l'eau
      ! de
      ! la fonte) (Laurent Li, le 14mars98):
      ! IM cf GK   tsol(i,is) = MIN(tsol(i,is),RTT-1.8)
      ! IM cf GK : la glace fond a 0C, non pas a -1.8
      tsol(i, is) = min(tsol(i,is), rtt)

      ! cc         ELSE
      ! cc            snow(i,is) = 0.0
      ! cc            qsol(i,is) = 0.0
      ! cc            tsol(i,is) = 0.0
    END IF

    is = is_sic
    qsol(i, is) = 0.0
    snow(i, is) = snow(i, is) + snow_fall(i)*dtime*pctsrf(i, is)
    IF (pctsrf(i,is)>epsfra) THEN
      subli = min(evap(i)*dtime, snow(i,is))
      snow(i, is) = snow(i, is) - subli
      fsno = min(max((tsol(i,is)-rtt)/chasno,0.0), snow(i,is))
      snow(i, is) = snow(i, is) - fsno
      tsol(i, is) = tsol(i, is) - fsno*chasno
      ! je limite la temperature a RTT-1.8 (il faudrait aussi prendre l'eau
      ! de
      ! la fonte) (Laurent Li, le 14mars98):
      ! IM cf GK   tsol(i,is) = MIN(tsol(i,is),RTT-1.8)
      ! IM cf GK : la glace fond a 0C, non pas a -1.8
      tsol(i, is) = min(tsol(i,is), rtt)

      ! cc         ELSE
      ! cc            snow(i,is) = 0.0
      ! cc            tsol(i,is) = 0.0
    END IF

    agesno(i) = (agesno(i)+(1.-agesno(i)/50.)*dtime/86400.)* &
      exp(-1.*max(0.0,snow_fall(i))*dtime/0.3)
    agesno(i) = max(agesno(i), 0.0)

  END DO

  RETURN
END SUBROUTINE hydrol
