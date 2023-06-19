      SUBROUTINE perturb_radlwsw(zxtsol,iflag_radia)
!
!Case-specific radiative setup
!
      use dimphy
      IMPLICIT none
      INCLUDE "flux_arp.h"     
!
! Arguments :
!------------
      REAL,DIMENSION(klon), INTENT(INOUT) :: zxtsol
      INTEGER, INTENT(IN)                 :: iflag_radia
!
!======================================================================
!
      IF (iflag_radia == 2) THEN
!
! Iflag_radia = 2 : DICE case : 
!               on force zxtsol=tg pour le rayonnement (MPL 20130806)
!Sonia : Cas dice lmdz1d force en flux : on impose Tg pour imposer le LWUP
!
        zxtsol(:) =  tg
      ENDIF
!
      RETURN
      END

