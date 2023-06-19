SUBROUTINE GPPREF(KPROMA,KSTART,KPROF,KFLEV,PVAH,PVBH,PALPH,PRESH,PRESF)

!**** *GPPREF* - Computes full level pressure

!     Purpose.
!     --------
!           Computes pressures at half and full model levels.

!**   Interface.
!     ----------
!        *CALL* *GPPREF(...)

!        Explicit arguments :
!        --------------------
!                              KPROMA :  dimensioning
!                              KSTART :  start of work
!                              KPROF  :  depth of work
!                              KFLEV     : vert. dimensioning
!                              PVAH(KFLEV),PVBH(KFLEV)- vertical coordinate
!                              PALPH (KPROMA,KFLEV)  - COEFF OF THE HYDROST
!                              PRESH(KPROMA,0:KFLEV) - HALF LEVEL PRESSURE
!                              PRESF(KPROMA,KFLEV)   - FULL LEVEL PRESSURE
!
!        Implicit arguments :  NONE.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!                                PHk*ln(PHk) - PHk-1*ln(PHk-1)
!     Full level P: ln(PFk) = [ ------------------------------- - 1. ]
!                                        PHk - PHk-1

!     which simplifies to:  PFk = Pk+1/2 * exp(-ALPHA)

!     In case of NDLNPR=1 it becomes even simpler (no need of LAPRXP any
!     more in principle !) :
!                           PFk = Pk+1/2 * (1.-ALPHA) except at the top
!     level :
!                           PF1 = P1.5 / (2+Cv/R)

!     Author.
!     -------
!        Erik Andersson, Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 92-11-23
!        Modified : 95-01-31 by Radmila Bubnova: correction in the case
!                            of the other approximation of d (ln p).
!        Modified : 00-11-22 by Agathe Untch: modifications for vertical
!                            finite elements
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Modified : 04-11-15 by K. YESSAD: improve the hierarchy of tests
!        Modified : 15-Feb-2005 by K. YESSAD: ZTOPPRES becomes TOPPRES
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCT0   , ONLY : LAPRXPK
USE YOMDYN   , ONLY : NDLNPR
USE YOMCST   , ONLY : RD       ,RCVD
USE YOMCVER  , ONLY : LVERTFE
USE YOMGEM   , ONLY : VAF      ,VBF      ,TOPPRES

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)                  :: PVAH(0:KFLEV) ! Argument NOT used
REAL(KIND=JPRB)                  :: PVBH(0:KFLEV) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRESF(KPROMA,KFLEV) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IFIRST, JLEV, JLON
REAL(KIND=JPRB) :: ZMUL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPREF',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    Level to begin normal computations
!              ----------------------------------

! This is introduced to allow the use of GPPREF without the implicit
! assumption that the top level input for pressure is 0 hPa.
! This restriction is only necessary in the case of use of NDLNPR=1.
!
! LVERTFE : .T./.F. Finite element/conventional vertical discretisation.
! NDLNPR  : NDLNPR=0: conventional formulation of delta, i.e. ln(P(l)/P(l-1)).
!           NDLNPR=1: formulation of delta used in non hydrostatic model,
! LAPRXPK : way of computing full-levels pressures in primitive equation
!
LVERTFE=.TRUE.    !!!!! A REVOIR (MPL) comment faut-il vraiment calculer PRESF ?

IF ((.NOT.LVERTFE).AND.(NDLNPR == 1)) THEN
  IF(PRESH(KSTART,0) <= TOPPRES)THEN
    IFIRST=2
  ELSE
    IFIRST=1
    DO JLON=KSTART,KPROF
      IF(PRESH(JLON,0) <= TOPPRES)THEN
        IFIRST=2
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       2.    COMPUTES FULL LEVEL PRESSURES.
!              ------------------------------

IF (LVERTFE) THEN
  DO JLEV=1,KFLEV
!   print *,'GPPREF: LVERTFE KFLEV KSTART KPROF JLEV',LVERTFE,KFLEV,KSTART,KPROF,JLEV
    PRESF(KSTART:KPROF,JLEV)=VAF(JLEV)+VBF(JLEV)*PRESH(KSTART:KPROF,KFLEV)  
  ENDDO
ELSE
  IF (NDLNPR == 0) THEN
    IF (LAPRXPK) THEN
      DO JLEV=1,KFLEV
        DO JLON=KSTART,KPROF
          PRESF(JLON,JLEV)=(PRESH(JLON,JLEV-1)+PRESH(JLON,JLEV))*0.5_JPRB
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,KFLEV
        DO JLON=KSTART,KPROF
          PRESF(JLON,JLEV)=EXP(-PALPH(JLON,JLEV))*PRESH(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF
  ELSEIF (NDLNPR == 1) THEN
    DO JLEV=IFIRST,KFLEV
      DO JLON=KSTART,KPROF
        PRESF(JLON,JLEV)=(1.0_JPRB-PALPH(JLON,JLEV))*PRESH(JLON,JLEV)
      ENDDO
    ENDDO
    ZMUL=1.0_JPRB/(2.0_JPRB+RCVD/RD)
    DO JLEV=1,IFIRST-1
      DO JLON=KSTART,KPROF
        PRESF(JLON,JLEV)=PRESH(JLON,JLEV)*ZMUL
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPREF',1,ZHOOK_HANDLE)
END SUBROUTINE GPPREF
