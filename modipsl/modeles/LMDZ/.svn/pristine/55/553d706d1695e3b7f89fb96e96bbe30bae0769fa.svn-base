SUBROUTINE solarlong(pday, psollong, pdist_sol)

  USE ioipsl
  USE print_control_mod, ONLY: lunout

  IMPLICIT NONE

  ! =======================================================================

  ! Objet:
  ! ------

  ! Calcul de la distance soleil-planete et de la declinaison
  ! en fonction du jour de l'annee.


  ! Methode:
  ! --------

  ! Calcul complet de l'elipse

  ! Interface:
  ! ----------

  ! Uncommon comprenant les parametres orbitaux.

  ! Arguments:
  ! ----------

  ! Input:
  ! ------
  ! pday          jour de l'annee (le jour 0 correspondant a l'equinoxe)
  ! lwrite        clef logique pour sorties de controle

  ! Output:
  ! -------
  ! pdist_sol     distance entre le soleil et la planete
  ! ( en unite astronomique pour utiliser la constante
  ! solaire terrestre 1370 Wm-2 )
  ! pdecli        declinaison ( en radians )

  ! =======================================================================
  ! -----------------------------------------------------------------------
  ! Declarations:
  ! -------------

  include "planete.h"
  include "YOMCST.h"

  ! arguments:
  ! ----------

  REAL pday, pdist_sol, pdecli, psollong
  LOGICAL lwrite

  ! Local:
  ! ------

  REAL zanom, xref, zx0, zdx, zteta, zz, pi
  INTEGER iter
  REAL :: pyear_day, pperi_day
  REAL :: jd_eq, jd_peri
  LOGICAL, SAVE :: first = .TRUE.
  !$OMP THREADPRIVATE(first)

  ! -----------------------------------------------------------------------
  ! calcul de l'angle polaire et de la distance au soleil :
  ! -------------------------------------------------------

  ! Initialisation eventuelle:
  IF (first) THEN
    CALL ioget_calendar(pyear_day)
    CALL ymds2ju(2000, 3, 21, 0., jd_eq)
    CALL ymds2ju(2001, 1, 4, 0., jd_peri)
    pperi_day = jd_peri - jd_eq
    pperi_day = r_peri + 180.
    WRITE (lunout, *) ' Number of days in a year = ', pyear_day
    ! call iniorbit(249.22,206.66,669.,485.,25.2)
    CALL iniorbit(152.59, 146.61, pyear_day, pperi_day, r_incl)
    first = .FALSE.
  END IF

  ! calcul de l'zanomalie moyenne

  zz = (pday-peri_day)/year_day
  pi = 2.*asin(1.)
  zanom = 2.*pi*(zz-nint(zz))
  xref = abs(zanom)

  ! resolution de l'equation horaire  zx0 - e * sin (zx0) = xref
  ! methode de Newton

  ! zx0=xref+e_elips*sin(xref)
  zx0 = xref + r_ecc*sin(xref)
  DO iter = 1, 10
    ! zdx=-(zx0-e_elips*sin(zx0)-xref)/(1.-e_elips*cos(zx0))
    zdx = -(zx0-r_ecc*sin(zx0)-xref)/(1.-r_ecc*cos(zx0))
    IF (abs(zdx)<=(1.E-7)) GO TO 120
    zx0 = zx0 + zdx
  END DO
120 CONTINUE
  zx0 = zx0 + zdx
  IF (zanom<0.) zx0 = -zx0

  ! zteta est la longitude solaire

  ! zteta=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))
  zteta = 2.*atan(sqrt((1.+r_ecc)/(1.-r_ecc))*tan(zx0/2.))

  psollong = zteta - timeperi

  IF (psollong<0.) psollong = psollong + 2.*pi
  IF (psollong>2.*pi) psollong = psollong - 2.*pi

  psollong = psollong*180./pi

  ! distance soleil

  pdist_sol = (1-r_ecc*r_ecc)/(1+r_ecc*cos(pi/180.*(psollong- &
    (r_peri+180.0))))
  ! pdist_sol = (1-e_elips*e_elips)
  ! &      /(1+e_elips*COS(pi/180.*(psollong-(R_peri+180.0))))
  ! -----------------------------------------------------------------------
  ! sorties eventuelles:
  ! ---------------------

  ! IF (lwrite) THEN
  ! PRINT*,'jour de l"annee   :',pday
  ! PRINT*,'distance au soleil (en unite astronomique) :',pdist_sol
  ! PRINT*,'declinaison (en degres) :',pdecli*180./pi
  ! ENDIF

  RETURN
END SUBROUTINE solarlong
