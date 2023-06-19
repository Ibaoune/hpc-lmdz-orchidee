! $Id: inifis_mod.F90 3435 2019-01-22 15:21:59Z fairhead $
MODULE inifis_mod

CONTAINS

  SUBROUTINE inifis(punjours, prad, pg, pr, pcpp)
  ! Initialize some physical constants and settings
  USE init_print_control_mod, ONLY : init_print_control
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE

  include "YOMCST.h"
  REAL,INTENT(IN) :: prad, pg, pr, pcpp, punjours

  CHARACTER (LEN=20) :: modname = 'inifis'
  CHARACTER (LEN=80) :: abort_message

  ! Initialize flags lunout, prt_level, debug
  CALL init_print_control

  ! suphel => initialize some physical constants (orbital parameters,
  !           geoid, gravity, thermodynamical constants, etc.) in the
  !           physics
  CALL suphel

  ! check that physical constants set in 'suphel' are coherent
  ! with values set in the dynamics:
  IF (rday/=punjours) THEN
    WRITE (lunout, *) 'inifis: length of day discrepancy!!!'
    WRITE (lunout, *) '  in the dynamics punjours=', punjours
    WRITE (lunout, *) '   but in the physics RDAY=', rday
    IF (abs(rday-punjours)>0.01*punjours) THEN
        ! stop here if the relative difference is more than 1%
      abort_message = 'length of day discrepancy'
      CALL abort_physic(modname, abort_message, 1)
    END IF
  END IF
  IF (rg/=pg) THEN
    WRITE (lunout, *) 'inifis: gravity discrepancy !!!'
    WRITE (lunout, *) '     in the dynamics pg=', pg
    WRITE (lunout, *) '  but in the physics RG=', rg
    IF (abs(rg-pg)>0.01*pg) THEN
        ! stop here if the relative difference is more than 1%
      abort_message = 'gravity discrepancy'
      CALL abort_physic(modname, abort_message, 1)
    END IF
  END IF
  IF (ra/=prad) THEN
    WRITE (lunout, *) 'inifis: planet radius discrepancy !!!'
    WRITE (lunout, *) '   in the dynamics prad=', prad
    WRITE (lunout, *) '  but in the physics RA=', ra
    IF (abs(ra-prad)>0.01*prad) THEN
        ! stop here if the relative difference is more than 1%
      abort_message = 'planet radius discrepancy'
      CALL abort_physic(modname, abort_message, 1)
    END IF
  END IF
  IF (rd/=pr) THEN
    WRITE (lunout, *) 'inifis: reduced gas constant discrepancy !!!'
    WRITE (lunout, *) '     in the dynamics pr=', pr
    WRITE (lunout, *) '  but in the physics RD=', rd
    IF (abs(rd-pr)>0.01*pr) THEN
        ! stop here if the relative difference is more than 1%
      abort_message = 'reduced gas constant discrepancy'
      CALL abort_physic(modname, abort_message, 1)
    END IF
  END IF
  IF (rcpd/=pcpp) THEN
    WRITE (lunout, *) 'inifis: specific heat discrepancy !!!'
    WRITE (lunout, *) '     in the dynamics pcpp=', pcpp
    WRITE (lunout, *) '  but in the physics RCPD=', rcpd
    IF (abs(rcpd-pcpp)>0.01*pcpp) THEN
        ! stop here if the relative difference is more than 1%
      abort_message = 'specific heat discrepancy'
      CALL abort_physic(modname, abort_message, 1)
    END IF
  END IF

  END SUBROUTINE inifis
  
END MODULE inifis_mod
