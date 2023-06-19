SUBROUTINE SUSAT

!**** *SUSAT*   - INITIALIZE COMMON YOESAT

!     PURPOSE.
!     --------
!           INITIALIZE YOESAT, THE COMMON THAT CONTROLS THE
!           SIMULATION OF SATELLITE RADIANCES

!**   INTERFACE.
!     ----------
!        *CALL* *SUSAT

!        EXPLICIT ARGUMENTS :
!        --------------------
!            NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOESAT

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "IN CORE MODEL"

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 88-12-15
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1        ,ONLY : JPIM     ,JPRB
USE YOMHOOK         ,ONLY : LHOOK,   DR_HOOK

USE YOMLUN_IFSAUX   , ONLY : NULOUT
USE YOMCST          , ONLY : RPI
USE YOESAT          , ONLY : NGEO     ,RGALT    ,RGNAD    ,RGNOR    ,&
 & RGSOU    ,RGWST    ,RGEAS    ,LGEOSE   ,LGEOSW   ,&
 & LGMS     ,LINDSA   ,LMTO  

IMPLICIT NONE

INTEGER(KIND=JPIM) :: ISATEL, JSATEL

REAL(KIND=JPRB) :: ZDEGRAD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!*CALL COMDOC
!----------------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

IF (LHOOK) CALL DR_HOOK('SUSAT',0,ZHOOK_HANDLE)
ISATEL=5
DO JSATEL = 1 , ISATEL
  RGALT(JSATEL) = 0.0_JPRB
  RGNAD(JSATEL) = 0.0_JPRB
  RGNOR(JSATEL) = 0.0_JPRB
  RGSOU(JSATEL) = 0.0_JPRB
  RGWST(JSATEL) = 0.0_JPRB
  RGEAS(JSATEL) = 0.0_JPRB
ENDDO

IF (NGEO /= 0) THEN

!      ----------------------------------------------------------------

!*       2.    MODIFY DEFAULT VALUES FOR THE VARIOUS GEO.SATELLITES
!              ----------------------------------------------------

  ISATEL = 0
  ZDEGRAD = RPI / 180._JPRB

  WRITE(UNIT=NULOUT,FMT='('' COMMON YOESAT '')')
  WRITE(UNIT=NULOUT,FMT='('' NGEO  = '',I1 )') NGEO
!      ----------------------------------------------------------------

!*       2.1   GOES EAST SATELLITE
!              -------------------

  IF (LGEOSE) THEN
    ISATEL = ISATEL + 1
    RGALT(ISATEL) = 0.0_JPRB
    RGALT(ISATEL) = 35793000._JPRB
    RGNAD(ISATEL) = 285._JPRB * ZDEGRAD
    RGNOR(ISATEL) = +70._JPRB * ZDEGRAD
    RGSOU(ISATEL) = -70._JPRB * ZDEGRAD
    RGWST(ISATEL) = RGNAD(ISATEL) -70._JPRB * ZDEGRAD
    RGEAS(ISATEL) = RGNAD(ISATEL) +70._JPRB * ZDEGRAD
    WRITE(UNIT=NULOUT,FMT='('' LGOESE = '',L5 &
     & ,'' ALTITUDE  ='',F10.0 &
     & ,'' LONG.NADIR='',F9.6 &
     & ,'' LIMFOV N. ='',F9.6 &
     & ,'' S. ='',F9.6 &
     & ,'' W. ='',F9.6 &
     & ,'' E. ='',F9.6 &
     & )')&
     & LGEOSE,RGALT(ISATEL),RGNAD(ISATEL)&
     & ,RGNOR(ISATEL),RGSOU(ISATEL),RGWST(ISATEL),RGEAS(ISATEL)  
  ENDIF

!      ----------------------------------------------------------------

!*       2.2   GOES WEST SATELLITE
!              -------------------

  IF (LGEOSW) THEN
    ISATEL = ISATEL + 1
    RGALT(ISATEL) = 0.0_JPRB
    RGALT(ISATEL) = 35793000._JPRB
    RGNAD(ISATEL) = 225._JPRB * ZDEGRAD
    RGNOR(ISATEL) = +70._JPRB * ZDEGRAD
    RGSOU(ISATEL) = -70._JPRB * ZDEGRAD
    RGWST(ISATEL) = RGNAD(ISATEL) -70._JPRB * ZDEGRAD
    RGEAS(ISATEL) = RGNAD(ISATEL) +70._JPRB * ZDEGRAD
    WRITE(UNIT=NULOUT,FMT='('' LGEOSW = '',L5 &
     & ,'' ALTITUDE  ='',F10.0 &
     & ,'' LONG.NADIR='',F9.6 &
     & ,'' LIMFOV N. ='',F9.6 &
     & ,'' S. ='',F9.6 &
     & ,'' W. ='',F9.6 &
     & ,'' E. ='',F9.6 &
     & )')&
     & LGEOSW,RGALT(ISATEL),RGNAD(ISATEL)&
     & ,RGNOR(ISATEL),RGSOU(ISATEL),RGWST(ISATEL),RGEAS(ISATEL)  
  ENDIF

!      ----------------------------------------------------------------

!*       2.3   G.M.S. SATELLITE
!              ----------------

  IF (LGMS) THEN
    ISATEL = ISATEL + 1
    RGALT(ISATEL) = 0.0_JPRB
    RGALT(ISATEL) = 35793000._JPRB
    RGNAD(ISATEL) = 140._JPRB * ZDEGRAD
    RGNOR(ISATEL) = +70._JPRB * ZDEGRAD
    RGSOU(ISATEL) = -70._JPRB * ZDEGRAD
    RGWST(ISATEL) = RGNAD(ISATEL) -70._JPRB * ZDEGRAD
    RGEAS(ISATEL) = RGNAD(ISATEL) +70._JPRB * ZDEGRAD
    WRITE(UNIT=NULOUT,FMT='('' LGMS   = '',L5 &
     & ,'' ALTITUDE  ='',F10.0 &
     & ,'' LONG.NADIR='',F9.6 &
     & ,'' LIMFOV N. ='',F9.6 &
     & ,'' S. ='',F9.6 &
     & ,'' W. ='',F9.6 &
     & ,'' E. ='',F9.6 &
     & )')&
     & LGMS,RGALT(ISATEL),RGNAD(ISATEL)&
     & ,RGNOR(ISATEL),RGSOU(ISATEL),RGWST(ISATEL),RGEAS(ISATEL)  
  ENDIF

!      ----------------------------------------------------------------

!*       2.4   INDSAT SATELLITE
!              ----------------

  IF (LINDSA) THEN
    ISATEL = ISATEL + 1
    RGALT(ISATEL) = 0.0_JPRB
    RGALT(ISATEL) = 35793000._JPRB
! ????      RGNAD(ISATEL) = 70. * ZDEGRAD
    RGNAD(ISATEL) = 0.0_JPRB
    RGNOR(ISATEL) = +70._JPRB * ZDEGRAD
    RGSOU(ISATEL) = -70._JPRB * ZDEGRAD
    RGWST(ISATEL) = 0.0_JPRB
    RGEAS(ISATEL) = 0.0_JPRB
    WRITE(UNIT=NULOUT,FMT='('' LINDSA = '',L5 &
     & ,'' ALTITUDE  ='',F10.0 &
     & ,'' LONG.NADIR='',F9.6 &
     & ,'' LIMFOV N. ='',F9.6 &
     & ,'' S. ='',F9.6 &
     & ,'' W. ='',F9.6 &
     & ,'' E. ='',F9.6 &
     & )')&
     & LINDSA,RGALT(ISATEL),RGNAD(ISATEL)&
     & ,RGNOR(ISATEL),RGSOU(ISATEL),RGWST(ISATEL),RGEAS(ISATEL)  
  ENDIF

!      ----------------------------------------------------------------

!*       2.5   METEOSAT SATELLITE
!              ------------------

  IF (LMTO) THEN
    ISATEL = ISATEL + 1
    RGALT(ISATEL) = 35793000._JPRB
    RGNAD(ISATEL) = 0.0_JPRB * ZDEGRAD
    RGNOR(ISATEL) = +70._JPRB * ZDEGRAD
    RGSOU(ISATEL) = -70._JPRB * ZDEGRAD
    RGWST(ISATEL) = 2.0_JPRB * RPI - 70._JPRB * ZDEGRAD
    RGEAS(ISATEL) = +70._JPRB * ZDEGRAD
    WRITE(UNIT=NULOUT,FMT='('' LMTO   = '',L5 &
     & ,'' ALTITUDE  ='',F10.0 &
     & ,'' LONG.NADIR='',F9.6 &
     & ,'' LIMFOV N. ='',F9.6 &
     & ,'' S. ='',F9.6 &
     & ,'' W. ='',F9.6 &
     & ,'' E. ='',F9.6 &
     & )')&
     & LMTO,RGALT(ISATEL),RGNAD(ISATEL)&
     & ,RGNOR(ISATEL),RGSOU(ISATEL),RGWST(ISATEL),RGEAS(ISATEL)  
  ENDIF

ENDIF

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSAT',1,ZHOOK_HANDLE)
END SUBROUTINE SUSAT
