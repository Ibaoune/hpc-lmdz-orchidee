SUBROUTINE SUJFH

! Purpose :
! -------
!   *SUJFH* Read namelist NAMJFH & Initialize module YOMJFH

! Interface :
! ---------

! Externals :
! ---------
!   POSNAM

! Method :
! ------

! Reference :
! ---------

! Author :
! ------
!   15-Mar-2005 R. El Khatib  *METEO-FRANCE*

! Modifications :
! -------------

! End Modifications
!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN   , ONLY : NULOUT   ,NULNAM
USE YOMLUN   , ONLY : NULOUT   
USE YOMPLDSW , ONLY : LOPT_RS6K
USE YOMJFH   , ONLY : N_VMASS

IMPLICIT NONE

REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "posnam.intfb.h"

#include "namjfh.h"

!-----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUJFH',0,ZHOOK_HANDLE)

IF (LOPT_RS6K) THEN
  N_VMASS=8
! Ce qui concerne NAMJFH commente par MPL le 15.04.09
! CALL POSNAM(NULNAM,'NAMJFH')
! READ(NULNAM,NAMJFH)
ELSE
  N_VMASS=0
ENDIF

WRITE(UNIT=NULOUT,FMT='('' MODULE YOMJFH'')')
WRITE(UNIT=NULOUT,FMT='(''   N_VMASS ='',I2)') N_VMASS 


IF (LHOOK) CALL DR_HOOK('SUJFH',1,ZHOOK_HANDLE)

END SUBROUTINE SUJFH
