
! $Header$

SUBROUTINE suphel

  IMPLICIT NONE

  include "YOMCST.h"
  include "YOETHF.h"
  ! IM cf. JLD
  LOGICAL firstcall
  SAVE firstcall
  !$OMP THREADPRIVATE(firstcall)
  DATA firstcall/.TRUE./

  IF (firstcall) THEN
    PRINT *, 'suphel initialise les constantes du GCM'
    firstcall = .FALSE.
  ELSE
    PRINT *, 'suphel DEJA APPELE '
    RETURN
  END IF
  ! -----------------------------------------------------------------

  ! *       1.    DEFINE FUNDAMENTAL CONSTANTS.
  ! -----------------------------

  WRITE (UNIT=6, FMT='(''0*** Constants of the ICM   ***'')')
  rpi = 2.*asin(1.)
  rclum = 299792458.
  rhpla = 6.6260755E-34
  rkbol = 1.380658E-23
  rnavo = 6.0221367E+23
  WRITE (UNIT=6, FMT='('' *** Fundamental constants ***'')')
  WRITE (UNIT=6, FMT='(''           PI = '',E13.7,'' -'')') rpi
  WRITE (UNIT=6, FMT='(''            c = '',E13.7,''m s-1'')') rclum
  WRITE (UNIT=6, FMT='(''            h = '',E13.7,''J s'')') rhpla
  WRITE (UNIT=6, FMT='(''            K = '',E13.7,''J K-1'')') rkbol
  WRITE (UNIT=6, FMT='(''            N = '',E13.7,''mol-1'')') rnavo

  ! ----------------------------------------------------------------

  ! *       2.    DEFINE ASTRONOMICAL CONSTANTS.
  ! ------------------------------

  rday = 86400.
  rea = 149597870000.
  repsm = 0.409093

  rsiyea = 365.25*rday*2.*rpi/6.283076
  rsiday = rday/(1.+rday/rsiyea)
  romega = 2.*rpi/rsiday

  ! exp1      R_ecc = 0.05
  ! exp1      R_peri = 102.04
  ! exp1      R_incl = 22.5
  ! exp1      print*, 'Parametres orbitaux modifies'
  ! ref      R_ecc = 0.016724
  ! ref      R_peri = 102.04
  ! ref      R_incl = 23.5

  ! IM 161002 : pour avoir les ctes AMIP II
  ! IM 161002   R_ecc = 0.016724
  ! IM 161002   R_peri = 102.04
  ! IM 161002   R_incl = 23.5
  ! IM on mets R_ecc, R_peri, R_incl dans conf_phys.F90
  ! R_ecc = 0.016715
  ! R_peri = 102.7
  ! R_incl = 23.441

  WRITE (UNIT=6, FMT='('' *** Astronomical constants ***'')')
  WRITE (UNIT=6, FMT='(''          day = '',E13.7,'' s'')') rday
  WRITE (UNIT=6, FMT='('' half g. axis = '',E13.7,'' m'')') rea
  WRITE (UNIT=6, FMT='('' mean anomaly = '',E13.7,'' -'')') repsm
  WRITE (UNIT=6, FMT='('' sideral year = '',E13.7,'' s'')') rsiyea
  WRITE (UNIT=6, FMT='(''  sideral day = '',E13.7,'' s'')') rsiday
  WRITE (UNIT=6, FMT='(''        omega = '',E13.7,'' s-1'')') romega
  ! write(unit=6,fmt='('' excentricite = '',e13.7,''-'')')R_ecc
  ! write(unit=6,fmt='(''     equinoxe = '',e13.7,''-'')')R_peri
  ! write(unit=6,fmt='(''  inclinaison = '',e13.7,''-'')')R_incl

  ! ------------------------------------------------------------------

  ! *       3.    DEFINE GEOIDE.
  ! --------------

  rg = 9.80665
  ra = 6371229.
  r1sa = sngl(1.D0/dble(ra))
  WRITE (UNIT=6, FMT='('' ***         Geoide         ***'')')
  WRITE (UNIT=6, FMT='(''      Gravity = '',E13.7,'' m s-2'')') rg
  WRITE (UNIT=6, FMT='('' Earth radius = '',E13.7,'' m'')') ra
  WRITE (UNIT=6, FMT='('' Inverse E.R. = '',E13.7,'' m'')') r1sa

  ! -----------------------------------------------------------------

  ! *       4.    DEFINE RADIATION CONSTANTS.
  ! ---------------------------

  ! z.x.li      RSIGMA=2. * RPI**5 * RKBOL**4 /(15.* RCLUM**2 * RHPLA**3)
  rsigma = 2.*rpi**5*(rkbol/rhpla)**3*rkbol/rclum/rclum/15.
  ! IM init. dans conf_phys.F90   RI0=1365.
  WRITE (UNIT=6, FMT='('' ***        Radiation       ***'')')
  WRITE (UNIT=6, FMT='('' Stefan-Bol.  = '',E13.7,'' W m-2 K-4'' &
    &                                                         &
    &         )') rsigma
  ! IM init. dans conf_phys.F90   WRITE(UNIT=6,FMT='('' Solar const. =
  ! '',E13.7,'' W m-2'')')
  ! IM init. dans conf_phys.F90  S      RI0

  ! -----------------------------------------------------------------

  ! *       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
  ! ------------------------------------------

  r = rnavo*rkbol
  rmd = 28.9644
  rmo3 = 47.9942
  rmco2 = 44.011
  rmch4 = 16.043
  rmn2o = 44.013
  rmcfc11 = 137.3686
  rmcfc12 = 120.9140
  rmc   = 12.0107
  rmv = 18.0153
  rd = 1000.*r/rmd
  rv = 1000.*r/rmv
  rcpd = 3.5*rd
  rcvd = rcpd - rd
  rcpv = 4.*rv
  rcvv = rcpv - rv
  rkappa = rd/rcpd
  retv = rv/rd - 1.
  WRITE (UNIT=6, FMT='('' *** Thermodynamic, gas     ***'')')
  WRITE (UNIT=6, FMT='('' Perfect gas  = '',e13.7)') r
  WRITE (UNIT=6, FMT='('' Dry air mass = '',e13.7)') rmd
  WRITE (UNIT=6, FMT='('' Ozone   mass = '',e13.7)') rmo3
  WRITE (UNIT=6, FMT='('' CO2     mass = '',e13.7)') rmco2
  WRITE (UNIT=6, FMT='('' C       mass = '',e13.7)') rmc
  WRITE (UNIT=6, FMT='('' CH4     mass = '',e13.7)') rmch4
  WRITE (UNIT=6, FMT='('' N2O     mass = '',e13.7)') rmn2o
  WRITE (UNIT=6, FMT='('' CFC11   mass = '',e13.7)') rmcfc11
  WRITE (UNIT=6, FMT='('' CFC12   mass = '',e13.7)') rmcfc12
  WRITE (UNIT=6, FMT='('' Vapour  mass = '',e13.7)') rmv
  WRITE (UNIT=6, FMT='('' Dry air cst. = '',e13.7)') rd
  WRITE (UNIT=6, FMT='('' Vapour  cst. = '',e13.7)') rv
  WRITE (UNIT=6, FMT='(''         Cpd  = '',e13.7)') rcpd
  WRITE (UNIT=6, FMT='(''         Cvd  = '',e13.7)') rcvd
  WRITE (UNIT=6, FMT='(''         Cpv  = '',e13.7)') rcpv
  WRITE (UNIT=6, FMT='(''         Cvv  = '',e13.7)') rcvv
  WRITE (UNIT=6, FMT='(''      Rd/Cpd  = '',e13.7)') rkappa
  WRITE (UNIT=6, FMT='(''     Rv/Rd-1  = '',e13.7)') retv

  ! ----------------------------------------------------------------

  ! *       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.
  ! ---------------------------------------------

  rcw = rcpv
  WRITE (UNIT=6, FMT='('' *** Thermodynamic, liquid  ***'')')
  WRITE (UNIT=6, FMT='(''         Cw   = '',E13.7)') rcw

  ! ----------------------------------------------------------------

  ! *       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.
  ! --------------------------------------------

  rcs = rcpv
  WRITE (UNIT=6, FMT='('' *** thermodynamic, solid   ***'')')
  WRITE (UNIT=6, FMT='(''         Cs   = '',E13.7)') rcs

  ! ----------------------------------------------------------------

  ! *       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.
  ! ----------------------------------------------------

  rtt = 273.16
  rlvtt = 2.5008E+6
  rlstt = 2.8345E+6
  rlmlt = rlstt - rlvtt
  ratm = 100000.
  WRITE (UNIT=6, FMT='('' *** Thermodynamic, trans.  ***'')')
  WRITE (UNIT=6, FMT='('' Fusion point  = '',E13.7)') rtt
  WRITE (UNIT=6, FMT='(''        RLvTt  = '',E13.7)') rlvtt
  WRITE (UNIT=6, FMT='(''        RLsTt  = '',E13.7)') rlstt
  WRITE (UNIT=6, FMT='(''        RLMlt  = '',E13.7)') rlmlt
  WRITE (UNIT=6, FMT='('' Normal press. = '',E13.7)') ratm
  WRITE (UNIT=6, FMT='('' Latent heat :  '')')

  ! ----------------------------------------------------------------

  ! *       9.    SATURATED VAPOUR PRESSURE.
  ! --------------------------

  restt = 611.14
  rgamw = (rcw-rcpv)/rv
  rbetw = rlvtt/rv + rgamw*rtt
  ralpw = log(restt) + rbetw/rtt + rgamw*log(rtt)
  rgams = (rcs-rcpv)/rv
  rbets = rlstt/rv + rgams*rtt
  ralps = log(restt) + rbets/rtt + rgams*log(rtt)
  rgamd = rgams - rgamw
  rbetd = rbets - rbetw
  ralpd = ralps - ralpw

  ! ------------------------------------------------------------------

  ! *       10.   CONSTANTS FOR THERMODYNAMICAL FUNCTIONS.
  ! ----------------------------------------

  rvtmp2 = rcpv/rcpd - 1.
  rhoh2o = ratm/100.
  r2es = restt*rd/rv
  r3les = 17.269
  r3ies = 21.875
  r4les = 35.86
  r4ies = 7.66
  r5les = r3les*(rtt-r4les)
  r5ies = r3ies*(rtt-r4ies)

  ! ------------------------------------------------------------------

  ! *       10.   CONSTANTS FOR METHANE OXIDATION AND PHOTOLYSIS.
  ! -----------------------------------------------

  CALL SUMETHOX()

  RETURN
END SUBROUTINE suphel
