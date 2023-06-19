SUBROUTINE iniorbit(paphelie, pperiheli, pyear_day, pperi_day, pobliq)
  IMPLICIT NONE

  ! =======================================================================

  ! Auteur:
  ! -------
  ! Frederic Hourdin      22 Fevrier 1991

  ! Objet:
  ! ------
  ! Initialisation du sous programme orbite qui calcule
  ! a une date donnee de l'annee de duree year_day commencant
  ! a l'equinoxe de printemps et dont le perihelie se situe
  ! a la date peri_day, la distance au soleil et la declinaison.

  ! Interface:
  ! ----------
  ! - Doit etre appele avant d'utiliser orbite.
  ! - initialise une partie du common planete.h

  ! Arguments:
  ! ----------

  ! Input:
  ! ------
  ! aphelie       \   aphelie et perihelie de l'orbite
  ! periheli      /   en millions de kilometres.

  ! =======================================================================

  ! -----------------------------------------------------------------------
  ! Declarations:
  ! -------------

  include "planete.h"
  include "YOMCST.h"

  ! Arguments:
  ! ----------

  REAL paphelie, pperiheli, pyear_day, pperi_day, pobliq

  ! Local:
  ! ------

  REAL zxref, zanom, zz, zx0, zdx, pi
  INTEGER iter

  ! -----------------------------------------------------------------------

  pi = 2.*asin(1.)

  aphelie = paphelie
  periheli = pperiheli
  year_day = pyear_day
  obliquit = pobliq
  peri_day = pperi_day

  PRINT *, 'Perihelie en Mkm  ', periheli
  PRINT *, 'Aphelie  en Mkm   ', aphelie
  PRINT *, 'obliquite en degres  :', obliquit
  PRINT *, 'Jours dans l annee : ', year_day
  PRINT *, 'Date perihelie : ', peri_day
  unitastr = 149.597870
  e_elips = (aphelie-periheli)/(periheli+aphelie)
  p_elips = 0.5*(periheli+aphelie)*(1-e_elips*e_elips)/unitastr

  PRINT *, 'e_elips', e_elips
  PRINT *, 'p_elips', p_elips

  ! -----------------------------------------------------------------------
  ! calcul de l'angle polaire et de la distance au soleil :
  ! -------------------------------------------------------

  ! calcul de l'zanomalie moyenne

  zz = (year_day-pperi_day)/year_day
  zanom = 2.*pi*(zz-nint(zz))
  zxref = abs(zanom)
  PRINT *, 'zanom  ', zanom

  ! resolution de l'equation horaire  zx0 - e * sin (zx0) = zxref
  ! methode de Newton

  zx0 = zxref + r_ecc*sin(zxref)
  DO iter = 1, 100
    zdx = -(zx0-r_ecc*sin(zx0)-zxref)/(1.-r_ecc*cos(zx0))
    IF (abs(zdx)<=(1.E-12)) GO TO 120
    zx0 = zx0 + zdx
  END DO
120 CONTINUE
  zx0 = zx0 + zdx
  IF (zanom<0.) zx0 = -zx0
  PRINT *, 'zx0   ', zx0

  ! zteta est la longitude solaire

  timeperi = 2.*atan(sqrt((1.+r_ecc)/(1.-r_ecc))*tan(zx0/2.))
  PRINT *, 'longitude solaire du perihelie timeperi = ', timeperi

  RETURN
END SUBROUTINE iniorbit
