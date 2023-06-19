SUBROUTINE SU_AERW

!**** *SU_AERW*   - DEFINES INDICES AND PARAMETERS FOR VARIOUS AEROSOL VARIABLES

!     PURPOSE.
!     --------
!           INITIALIZE YOEAERATM, YOEAERSRC, YOEAERSNK, THE MODULES THAT CONTAINS INDICES
!           ALLOWING TO GET THE AEROSOL PARAMETERS RELEVANT FOR THE PROGNOSTIC AEROSOL
!           CONFIGURATION.

!**   INTERFACE.
!     ----------
!        *CALL* *SU_AERW

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        YOEAERW

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2005-07-08

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOEAERATM, ONLY : LAERGBUD, LAERSCAV ,LAERSEDIM, LAERSURF, LAER6SDIA, &
 & LAERCLIMG, LAERCLIMZ, LAERCLIST, LAERDRYDP, LAERNGAT, LAERPRNT, &
 & REPSCAER , INDBG

USE YOEAERSRC, ONLY : JKBIN, JKTYP, LEPAERO  , &
 & NBINAER, NINDAER , NMAXTAER, NTAER  , NTYPAER, RGELAV , RGEMUV , &
 & RDGLAV , RDGMUV  , RCLONV  , RSLONV , RDCLONV, RDSLONV, LAEREXTR, &
 & RSSFLX , RLATVOL , RLONVOL , NDDUST

USE YOEPHY   , ONLY : LE4ALB

USE YOEDBUG  , ONLY : KSTPDBG

USE YOMCST   , ONLY : RPI
USE YOMGC    , ONLY : GELAM
USE YOMGEM   , ONLY : NGPTOT
USE YOMLEG   , ONLY : RMU
! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN   , ONLY : NULNAM, NULOUT
USE YOMLUN   , ONLY : NULOUT

USE YOM_YGFL, ONLY : NAERO

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IAER, ICAER, ITAER
INTEGER(KIND=JPIM) :: J, JAER, JL

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!     ----------------------------------------------------------------

#include "posnam.intfb.h"

#include "su_aerp.intfb.h"
#include "su_aerop.intfb.h"
!      ----------------------------------------------------------------

#include "naeaer.h"

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SU_AERW',0,ZHOOK_HANDLE)

!*         1.       DEFAULT VALUES OF PARAMETERS
!                   ----------------------------

print *,'DANS SU_AERW'
NBINAER(:) = (/ 3, 3, 2, 2, 1, 1, 1, 1, 1 /)

NMAXTAER=9
NTYPAER(:) = 0
JKBIN(:) = 0
JKTYP(:) = 0

LEPAERO=.FALSE.
LAERCLIMG=.FALSE.
LAERCLIMZ=.FALSE.
LAERCLIST=.FALSE.
LAERDRYDP=.FALSE.
LAEREXTR=.FALSE.
LAERGBUD=.FALSE.
LAERNGAT=.FALSE.
LAERPRNT=.FALSE.
LAERSCAV=.FALSE.
LAERSEDIM=.FALSE.
LAERSURF=.FALSE.
LAER6SDIA=.FALSE.
INDBG=1
NTAER  =0
NDDUST =2

! the 9 types and assumed number of bins are:
!  NTYPAER    bins  type
!     1       1- 3  sea-salt  0.03 - 0.5 -  5  - 20 microns
!     2       4- 6  dust      0.03 - 0.5 - 0.9 - 20 microns
!     3       7- 8  POM	    hydrophilic, hydrophobic
!     4       9-10  BC	    hydrophilic, hydrophobic
!     5      11     sulfate
!     6      12     fly ash
!     7      13     pseudo-prognostic stratospheric aerosols
!     8      14     pseudo-prognostic volcanic aerosols
!     9      15     prognostic stratospheric aerosols
 
DO JAER=1,NMAXTAER
  NTYPAER(JAER)=0
ENDDO

RLATVOL=-999._JPRB
RLONVOL=-999._JPRB
RGELAV =-999._JPRB
RGEMUV =-999._JPRB
RDGLAV = 999._JPRB
RDGMUV = 999._JPRB
RCLONV =-999._JPRB
RSLONV =-999._JPRB
RDCLONV= 999._JPRB
RDSLONV= 999._JPRB
DO J=1,3
  KSTPDBG(J)=-999
ENDDO

REPSCAER=1.E-15_JPRB

!     ------------------------------------------------------------------

!*         2.       READ VALUES OF PROGNOSTIC AEROSOL CONFIGURATION
!                   -----------------------------------------------
! Ce qui concerne NAEAER commente par MPL le 15.04.09
!IF(NAERO > 0) THEN
! CALL POSNAM(NULNAM,'NAEAER')
! READ (NULNAM,NAEAER)
!ENDIF

IF (.NOT.LE4ALB) THEN
  NDDUST=2
ENDIF
!     ------------------------------------------------------------------

!*       3.    INITIALIZE PROGNOSTIC AEROSOL PHYSICAL AND OPTICAL PARAMETERS
!              -------------------------------------------------------------

  CALL SU_AERP
print *,'SU_AERW: apres SU_AERP'
  CALL SU_AEROP
print *,'SU_AERW: apres SU_AEROP'

IF (LEPAERO) THEN

! define a composite index for each bin of each different aerosol type to be used
! in source, sedimentation and deposition routines
 
  ICAER=0
  DO JAER=1,NMAXTAER
    IF (NTYPAER(JAER) /= 0) THEN
      NTAER=NTAER+1
      ITAER=NTYPAER(JAER)
      DO IAER=1,ITAER
        ICAER=ICAER+1
        NINDAER(ICAER)=JAER*10+IAER
        JKTYP(ICAER)=JAER
        JKBIN(ICAER)=IAER
      ENDDO
    ENDIF
  ENDDO

!-- if volcanic aerosols, define the model coordinates

  IF (NTYPAER(9) /= 0) THEN
    RGEMUV=(RLATVOL+90._JPRB)*RPI/180._JPRB
    RGELAV=RLONVOL*RPI/180._JPRB
    RCLONV=COS(RGELAV)
    RSLONV=SIN(RGELAV)
    DO J=1,NGPTOT-1
      IF (RGELAV > GELAM(J) .AND. RGELAV <= GELAM(J+1) .AND. &
        & RGEMUV < RMU(JL) .AND. RGEMUV >= RMU(JL+1) ) THEN
        RDGMUV=ABS( RMU(J+1) - RMU(J))
        RDGLAV=ABS( GELAM(J+1)-GELAM(J) )
        RDSLONV=ABS( SIN(GELAM(JL+1))-SIN(GELAM(JL)) )
        RDCLONV=ABS( COS(GELAM(JL+1))-COS(GELAM(JL)) )
      ENDIF
    ENDDO
  ENDIF  

!      ----------------------------------------------------------------

!*       4.    PRINT FINAL VALUES.
!              -------------------

  WRITE(UNIT=NULOUT,FMT='('' LEPAERO = '',L5 &
   & ,'' NTAER   = '',I2 ,'' NDDUST = '',I1,/&
   & ,'' NTYPAER = '',9I3,/ &
   & ,'' NBINAER = '',9I3,/ &
   & ,'' JKTYP   = '',15I3,/&
   & ,'' JKBIN   = '',15I3  &
   & )')&
   & LEPAERO,NTAER,(NTYPAER(JAER),JAER=1,9),(NBINAER(JAER),JAER=1,9), &
   & (JKTYP(JAER),JAER=1,15),(JKBIN(JAER),JAER=1,15)
  WRITE(UNIT=NULOUT,FMT='('' LAERGBUD = '',L3 &
   & ,'' LAERNGAT = '',L3 &
   & ,'' LAERDRYDP= '',L3 &
   & ,'' LAERSEDIM= '',L3 &
   & ,'' LAERSCAV = '',L3 &
   & ,'' LAER6SDIA= '',L3 &
   & ,'' LAERCLIMZ= '',L3 &
   & ,'' LAERCLIMG= '',L3 &
   & ,'' LAERCLIST= '',L3 &
   & )')&
   & LAERGBUD,LAERNGAT, LAERDRYDP,LAERSEDIM,LAERSCAV,LAER6SDIA,LAERCLIMZ,LAERCLIMG,LAERCLIST
  WRITE(UNIT=NULOUT,FMT='('' RSSFLX= '',10E10.3)') RSSFLX
  IF (NTYPAER(9) /= 0) THEN
    WRITE(UNIT=NULOUT,FMT='('' RLATVOL= '',F5.2 &
     & ,'' RLONVOL= '',F6.2,'' RGEMUV= '',F6.4,'' RGELAV= '',F6.4 &
     & ,'' RCLONV = '',F6.4,'' RSLONV= '',F6.4,'' RDGMUV= '',F6.4 &
     & ,'' RDGLAV = '',F6.4,'' RDCLONV= '',F6.4,'' RDSLONV= '',F6.4 &
     & )')&
     & RLATVOL,RLONVOL,RGEMUV,RGELAV,RCLONV,RSLONV,RDGMUV,RDGLAV,RDCLONV,RDSLONV
  ENDIF
ENDIF

!     ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SU_AERW',1,ZHOOK_HANDLE)
END SUBROUTINE SU_AERW




































