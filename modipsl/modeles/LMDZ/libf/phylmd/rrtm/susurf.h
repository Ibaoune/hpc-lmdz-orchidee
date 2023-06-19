SUBROUTINE SUSURF(KSW,KCSS,KSIL,KTILES,KTSW,&
 & LD_LLCCNL,LD_LLCCNO,LD_LEOCWA,LD_LEOCCO,LD_LEOCSA,LD_LLE4ALB,&
 & LD_LSCMEC,LD_LROUGH,PEXTZ0M,PEXTZ0H,&
 & PTHRFRTI,PTSTAND,PXP,PRCCNSEA,PRCCNLND,&
 & PRSUN)

!**   *SUSURF* IS THE SET-UP ROUTINE FOR surface modules containing constants

!     PURPOSE
!     -------
!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOESOIL*

!     INTERFACE.
!     ----------
!     CALL *SUSURF* FROM *SUPHEC*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     Original    A.C.M. BELJAARS         E.C.M.W.F.      89/11/02
!     MODIFICATIONS
!     -------------
!     J.-J. MORCRETTE         E.C.M.W.F.      91/07/14
!     P. VITERBO              E.C.M.W.F.       8/10/93
!     P. Viterbo     99-03-26    Tiling of the land surface
!     C. Fischer 00-12-20 Meteo-France recode initialization of rdat to avoid
!                         memory overflow on SUN workstation
!     J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     P. Viterbo   ECMWF   03-12-2004  Include user-defined RTHRFRTI
!     P. Viterbo   ECMWF   May 2005    Externalise surf
!        JJMorcrette 20060511 MODIS albedo

!  INTERFACE: 

!    Integers (In):

!      KSW       : NUMBER OF SHORTWAVE SPECTRAL INTERVALS
!      KCSS      : Number of soil levels
!      KSIL      : NUMBER OF (infrared) SPECTRAL INTERVALS
!      KTILES    : Number of surface tiles
!      KTSW      : Maximum possible number of sw spectral intervals

!    Logicals (In):

!      LD_LLCCNL : .T. IF CCN CONCENTRATION OVER LAND IS DIAGNOSED
!      LD_LLCCNO : .T. IF CCN CONCENTRATION OVER OCEAN IS DIAGNOSED
!      LD_LLE4ALB: .T. IF MODIS ALBEDO IS USED

!    Reals (In):

!      PTHRFRTI  : ! MINIMUM THRESHOLD FOR TILE FRACTION
!      PTSTAND   : ! REFERENCE TEMPERATURE FOR TEMPERATURE DEPENDENCE
!      PXP       : ! POLYNOMIAL COEFFICIENTS OF PLANCK FUNCTION
!      PRCCNSEA  : ! NUMBER CONCENTRATION (CM-3) OF CCNs OVER SEA
!      PRCCNLND  : ! NUMBER CONCENTRATION (CM-3) OF CCNs OVER LAND
!      PRSUN     : ! SOLAR FRACTION IN SPECTRAL INTERVALS

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KSW 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTHRFRTI
INTEGER(KIND=JPIM),INTENT(IN)    :: KCSS  
INTEGER(KIND=JPIM),INTENT(IN)    :: KSIL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTSW 
LOGICAL           ,INTENT(IN)    :: LD_LLCCNL 
LOGICAL           ,INTENT(IN)    :: LD_LLCCNO  
LOGICAL           ,INTENT(IN)    :: LD_LEOCWA
LOGICAL           ,INTENT(IN)    :: LD_LEOCCO
LOGICAL           ,INTENT(IN)    :: LD_LEOCSA
LOGICAL           ,INTENT(IN)    :: LD_LLE4ALB
LOGICAL           ,INTENT(IN)    :: LD_LSCMEC
LOGICAL           ,INTENT(IN)    :: LD_LROUGH
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXTZ0M
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXTZ0H
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTAND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXP(6,6) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCCNSEA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCCNLND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSUN(:) 

!     ------------------------------------------------------------------

END SUBROUTINE SUSURF
