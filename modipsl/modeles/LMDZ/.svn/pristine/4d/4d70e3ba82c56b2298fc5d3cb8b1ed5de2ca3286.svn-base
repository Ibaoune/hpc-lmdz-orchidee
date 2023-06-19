!OPTIONS XOPT(NOEVAL)
SUBROUTINE SUECRAD15 (KULOUT, KLEV, PETAH )

!**** *SUECRAD15*   - INITIALIZE COMMONS YOMRxx15 CONTROLLING RADIATION
!****                 FROZEN VERSION (CYCLE 15) OF SUECRAD

!     PURPOSE.
!     --------
!           INITIALIZE YOMRAD15, THE COMMON THAT CONTROLS THE
!           RADIATION OF THE MODEL, AND YOMRDU15 THAT INCLUDES
!           ADJUSTABLE PARAMETERS FOR RADIATION COMPUTATIONS

!**   INTERFACE.
!     ----------
!        CALL *SUECRAD15* FROM *SUPHEC*
!              ---------        ------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMONS YOMRAD15, YOMRDU15

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        SUAER, SUAERH, SUAERV, SULW, SUSW, SUCLD, SUOCST, SUSAT

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        96-11: Ph. Dandin. Meteo-France
!        ORIGINAL : 88-12-15 BY JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        R. El Khatib 01-02-02 proper initialization of NFRRC moved in SUCFU
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Bouyssel 27-09-04 initialisation of NSW
!        A. Alias    29-09-05 Sulfate aerosols (Hu Rong Ming)
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCT0   , ONLY : NPRINTLEV
USE YOMDIM   , ONLY : NDLON    ,NSMAX   ,NGPBLKS  ,NFLEVG   ,NPROMA
USE YOMDYN   , ONLY : TSTEP
! Ce qui concerne NULNAM commente par MPL le 15.04.09
!USE YOMLUN   , ONLY : NULNAM
USE YOMCST   , ONLY : RDAY     ,RG       ,RCPD
USE YOMPHY   , ONLY : LRAYFM15
USE YOEPHY   , ONLY : LEPHYS   ,LERADI
USE YOMRAD15 , ONLY : NAER15   ,NFLUX15  ,NMODE15  ,NRAD15   ,&
 & NRADFR15 ,NRADPFR15,NRADPLA15,NRINT15  ,NRADNFR15,&
 & NRADSFR15,NOVLP15  ,NRPROMA15,NRADF2C15,NRADC2F15,&
 & LERAD6H15,LERADHS15,LRADAER15,LNEWAER15
USE YOERAD   , ONLY : NAER     ,NTSW
!USE YOERAD   , ONLY : NAER     ,NSW      ,NTSW
! NSW mis dans .def MPL 20140211
USE YOMRDU15 , ONLY : NUAER15  ,NTRAER15 ,RCDAY15  ,R10E15   ,&
 & REELOG15 ,REPSC15  ,REPSCO15 ,REPSCQ15 ,REPSCT15 ,&
 & REPSCW15 ,DIFF15  
USE YOMAERD15, ONLY : CVDAES15 ,CVDAEL15 ,CVDAEU15 ,CVDAED15 ,&
 & CVDAEF15 ,&
 & RCAEOPS15,RCAEOPL15,RCAEOPU15,RCAEOPD15,RCTRBGA15,&
 & RCAEOPF15,&
 & RCVOBGA15,RCSTBGA15,RCTRPT15 ,RCAEADM15,RCAEROS15,&
 & RCAEADK15  
USE YOMPRAD  , ONLY : LODBGRADI,LODBGRADL
USE YOMRADF  , ONLY : EMTD     ,EMTU      ,TRSW    ,RMOON

IMPLICIT NONE

include "clesphys.h"

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PETAH(KLEV+1)
LOGICAL :: LLP

#include "namrad15.h"
!      ----------------------------------------------------------------

LOGICAL :: LLMESS

INTEGER(KIND=JPIM) :: IRADFR, IST1HR, IST6HR


REAL(KIND=JPRB) :: ZSTPHR, ZTSTEP
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "posnam.intfb.h"
#include "suaer15.intfb.h"
#include "suaerv15.intfb.h"
#include "suecradi15.intfb.h"
#include "suecradl.intfb.h"
#include "sulw15.intfb.h"
#include "surdi15.intfb.h"
#include "susat.intfb.h"
#include "susw15.intfb.h"

!      ----------------------------------------------------------------

!*         1.       SET DEFAULT VALUES.
!                   -------------------

!*         1.1      PRESET INDICES IN *YOMRAD15*
!                   --------------------------

IF (LHOOK) CALL DR_HOOK('SUECRAD15',0,ZHOOK_HANDLE)
LLMESS=.FALSE.
LERAD6H15=.TRUE.
LERADHS15=.TRUE.
LRADAER15=.TRUE.
LNEWAER15=.FALSE.
NAER15   =1
NAER=0
NFLUX15  =6
NMODE15  =0
NRAD15   =1
NRADFR15 =-3
NRADPFR15=36
NRADPLA15=15
NRINT15  =4
NRADF2C15=1
NRADC2F15=1
NUAER15  = 24
NTRAER15 = 15
NOVLP15 = 1
NSW=2
NTSW=2
IF(NSMAX >= 106) THEN
  NRPROMA15 = 80
ELSEIF(NSMAX == 63) THEN
  NRPROMA15=48
ELSE
  NRPROMA15=20
ENDIF

!*         1.3      SET SECURITY PARAMETERS
!                   -----------------------

REPSC15  = 1.E-12_JPRB
REPSCO15 = 1.E-12_JPRB
REPSCQ15 = 1.E-12_JPRB
REPSCT15 = 1.E-12_JPRB
REPSCW15 = 1.E-12_JPRB
REELOG15 = 1.E-12_JPRB

!     ------------------------------------------------------------------

!*         2.       READ VALUES OF RADIATION CONFIGURATION
!                   --------------------------------------

! Ce qui concerne NAMRAD15 commente par MPL le 15.04.09
!CALL POSNAM(NULNAM,'NAMRAD15')
!READ (NULNAM,NAMRAD15)

!     INITIALISE DATA STRUCTURES REQUIRED FOR RADIATION INTERPOLATION

LODBGRADI=.FALSE.
CALL SUECRADI15

IF( LLMESS )THEN

!     INITIALISE DATA STRUCTURES REQUIRED FOR RADIATION COURSE GRID
!     LOAD BALANCING

  LODBGRADL=.FALSE.
! CALL SUECRADL   ! MPL 1.12.08
  CALL ABOR1('JUSTE APRES CALL SUECRADL COMMENTE')
ENDIF

!      ----------------------------------------------------------------

!*       3.    INITIALIZE RADIATION COEFFICIENTS.
!              ----------------------------------

RCDAY15  = RDAY * RG / RCPD
DIFF15   = 1.66_JPRB
R10E15   = 0.4342945_JPRB

CALL SURDI15

!      ----------------------------------------------------------------

!*       4.    INITIALIZE RADIATION ABSORPTION COEFFICIENTS
!              --------------------------------------------

CALL SULW15
CALL SUSW15

!      ----------------------------------------------------------------

!*       5.    INITIALIZE AEROSOL OPTICAL PARAMETERS AND DISTRIBUTION
!              ------------------------------------------------------

!     INITIALIZATION DONE IN BLOCK DATA SUAERHBD

!- optical properties
CALL SUAER15

!     CALL SUAERH

CALL SUAERV15 ( KLEV  , PETAH,&
 & CVDAES15, CVDAEL15, CVDAEU15, CVDAED15, CVDAEF15,&
 & RCTRBGA15, RCVOBGA15, RCSTBGA15, RCAEOPS15, RCAEOPL15, RCAEOPU15,&
 & RCAEOPF15,&
 & RCAEOPD15, RCTRPT15 , RCAEADK15, RCAEADM15, RCAEROS15        )  

!      ----------------------------------------------------------------

!*       6.    INITIALIZE SATELLITE GEOMETRICAL/RADIOMETRIC PARAMETERS
!              -------------------------------------------------------

IF (LEPHYS) THEN
  IF (NMODE15 > 1) THEN
    CALL SUSAT
  ENDIF
ENDIF

!      ----------------------------------------------------------------

!*       7.    INITIALIZE CLIMATOLOGICAL OZONE DISTRIBUTION
!              --------------------------------------------
!                  (not done here!!!  called from APLPAR as it depends
!                     on model pressure levels!)

!      ----------------------------------------------------------------

!*       8.    SET UP MODEL CONFIGURATION FOR TIME-SPACE INTERPOLATION
!              -------------------------------------------------------

ZTSTEP=MAX(TSTEP,1.0_JPRB)
ZSTPHR=3600._JPRB/ZTSTEP
IRADFR=NRADFR15
IF(NRADFR15 < 0) THEN
  NRADFR15=-NRADFR15*ZSTPHR+0.5_JPRB
ENDIF
NRADPFR15=NRADPFR15*NRADFR15
IF (MOD(NRADPLA15,2) == 0.AND. NRADPLA15 /= 0) THEN
  NRADPLA15=NRADPLA15+1
ENDIF

IST1HR=ZSTPHR+0.05_JPRB
IST6HR=6._JPRB*ZSTPHR+0.05_JPRB
IF (MOD(3600._JPRB,ZTSTEP) > 0.1_JPRB) THEN
  IST1HR=IST1HR+1
  DO WHILE (MOD(IST6HR,IST1HR) /= 0)
    IST1HR=IST1HR+1
  ENDDO
ENDIF
NRADSFR15=IST1HR
NRADNFR15=NRADFR15

IF(LRAYFM15) THEN
  NRPROMA15=NDLON+6+(1-MOD(NDLON,2))
ENDIF



!*       9.    ALLOCATE WORK ARRAYS
!               --------------------

LLP = NPRINTLEV >= 1

  ALLOCATE(EMTD(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(UNIT=KULOUT,FMT=9) 'EMTD     ',SIZE(EMTD     ),SHAPE(EMTD     )
  ALLOCATE(TRSW(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(UNIT=KULOUT,FMT=9) 'TRSW     ',SIZE(TRSW     ),SHAPE(TRSW     )
  ALLOCATE(EMTU(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(UNIT=KULOUT,FMT=9) 'EMTU     ',SIZE(EMTU     ),SHAPE(EMTU     )
  ALLOCATE(RMOON(NPROMA,NGPBLKS))
  IF(LLP)WRITE(UNIT=KULOUT,FMT=9) 'RMOON    ',SIZE(RMOON    ),SHAPE(RMOON    )

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!      ----------------------------------------------------------------

!*       9.    PRINT FINAL VALUES.
!              -------------------




IF (LLP) THEN
  WRITE(UNIT=KULOUT,FMT='('' COMMON YOMRAD15 '')')
  WRITE(UNIT=KULOUT,FMT='('' LERADI  = '',L5 &
   & ,'' LERAD6H15 = '',L5)')&
   & LERADI,LERAD6H15  
  WRITE(UNIT=KULOUT,FMT='('' NRADFR15  = '',I2 &
   & ,'' NRADPFR15 = '',I3 &
   & ,'' NRADPLA15 = '',I2 &
   & ,'' NRINT15   = '',I1 &
   & ,'' NRPROMA15 = '',I5 &
   & ,'' NRADF2C15 = '',I1 &
   & ,'' NRADC2F15 = '',I1 &
   & )')&
   & NRADFR15,NRADPFR15,NRADPLA15,NRINT15,&
   & NRPROMA15,NRADF2C15,NRADC2F15  
   
  WRITE(UNIT=KULOUT,FMT='('' LERADHS15= '',L5,'' LRADAER15= '',L5 &
   & ,'' LNEWAER15= '',L5 &
   & ,'' NMODE15 = '',I1 &
   & ,'' NAER15  = '',I1 &
   & ,'' NFLUX15 = '',I2 &
   & ,'' NRAD15  = '',I2 &
   & )')&
   & LERADHS15,LRADAER15,LNEWAER15,NMODE15,NAER15,NFLUX15,NRAD15  
  WRITE(KULOUT,FMT='('' WARNING! CLOUD OVERLAP ASSUMPTION IS''&
   & ,'' NOVLP15   = '',I2 &
   & )')&
   & NOVLP15  

  WRITE(UNIT=KULOUT,FMT='('' MODULE YOERAD '')')
  WRITE(UNIT=KULOUT,FMT='('' NSW = '',I2, '' NTSW = '',I2)') NSW,NTSW
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUECRAD15',1,ZHOOK_HANDLE)
END SUBROUTINE SUECRAD15
