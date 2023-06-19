! $Id$
MODULE inifis_mod

CONTAINS

  SUBROUTINE inifis(prad, pg, pr, pcpp)
  ! Initialize some physical constants and settings
  USE init_print_control_mod, ONLY: init_print_control
  USE comcstphy, ONLY: rradius, & ! planet radius (m)
                       rr, & ! recuced gas constant: R/molar mass of atm
                       rg, & ! gravity
                       rcpp  ! specific heat of the atmosphere
  IMPLICIT NONE

  REAL,INTENT(IN) :: prad, pg, pr, pcpp

  ! Initialize flags lunout, prt_level, debug
  CALL init_print_control

  ! copy some fundamental parameters to physics
  rradius=prad
  rg=pg
  rr=pr
  rcpp=pcpp

  END SUBROUTINE inifis
  
END MODULE inifis_mod
