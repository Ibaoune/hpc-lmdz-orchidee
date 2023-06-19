      subroutine Atm_RT_INI

!--------------------------------------------------------------------------+
!                                                     Mon 24-Jun-2013  MAR |
!     subroutine Atm_RT_INI initializes  MAR PHYsics interface             |
!           with ECMWF Solar/Infrared   Radiative  Transfer Scheme         |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  2-Apr-2013      |
!           Last Modification by H. Gallee,           Mon 24-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     Content:   CALL of - ECMWF Code initializing the Radiation Transfert |
!                        - ECMWF                       Radiation Transfert |
!                                                                          |
!     ECMWF Code Source:  J.-J. Morcrette, 28 nov 2002                     |
!                                                                          |
!--------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY_RT_dat
      use Mod_PHY____grd
      use Mod_PHY_S0_grd
      use Mod_PHY_RT_grd
      use Mod_PHY_DY_kkl
      use Mod_SISVAT_gpt



#include "tsmbkind.h"




! Global Variables (ECMWF)
! ========================

USE PARRTM1D , ONLY : JP_LON   ,JP_IDIA  ,JP_FDIA  ,JP_TDIA  ,&
 &            JP_LEV ,JP_LW    ,JP_SW    ,JP_NUA   ,JP_MODE  ,&
 &            JP_AER ,JP_LEVP1
USE YOMCST   , ONLY : RD       ,RG       ,RTT      ,RSIGMA   ,&
 &            RCPD   ,RPI      ,RDAY     ,REA      ,RI0      ,&
 &            REPSM  ,RMD      ,RKBOL    ,RNAVO    ,R        ,&
 &            RLVTT  ,RLSTT
USE YOERAD   , ONLY : NSW      ,NTSW     ,NRADFR             ,&
 &            LRRTM  ,LINHOM   ,LOIFUEC  ,LTEMPDS  ,LOWASYF  ,&
 &            LOWHSSS,LONEWSW  ,LNEWAER  ,LHVOLCA            ,&
 &            NRADIP ,NRADLP   ,NOZOCL                       ,&
 &            NICEOPT,NLIQOPT  ,NOVLP    ,NHOWINH  ,RMINICE
USE YOEAERD  , ONLY : CVDAES   ,CVDAEL   ,CVDAEU   ,CVDAED            ,&
 &            RCTRBGA,RCVOBGA  ,RCSTBGA  ,RCAEOPS  ,RCAEOPL  ,RCAEOPU ,&
 &            RCAEOPD,RCTRPT   ,RCAEADK  ,RCAEADM  ,RCAEROS
USE YOERDI   , ONLY : RCH4     ,RN2O     ,RO3      ,RCFC11   ,&
 &                                                  RCFC12
USE YOERDU   , ONLY : RCDAY    ,R10E     ,DIFF     ,REPLOG   ,&
 &            REPSC  ,REPSCO   ,REPSCQ   ,REPSCT   ,REPSCW   ,&
 &            NTRAER
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &            R4IES  ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 &            RALVDCP,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU




      IMPLICIT NONE




!  INTERNAL VARIABLES
!  ==================

      integer   :: i     ,j     ,ikl

!   For Use in radlsw
!   ^^^^^^^^^^^^^^^^^
      INTEGER_M :: KIDIA ,KFDIA ,KTDIA ,KLON  ,KLEV  
      INTEGER_M :: KMODE ,KAER  ,KSW

      INTEGER_M :: KBOX  ,NBOX
      INTEGER_M :: NDUMP ,ILWRAD
!
      REAL_B    ::    PRII05

      REAL_B    ::     PAER5(JP_LON,JP_AER,JP_LEV)   ! Aerosol Optical Depth
      REAL_B    ::    PALBD5(JP_LON,JP_SW)  
      REAL_B    ::    PALBP5(JP_LON,JP_SW)
!
      REAL_B    ::     PAPH5(JP_LON,JP_LEVP1) 
      REAL_B    ::      PAP5(JP_LON,JP_LEV)
!
      REAL_B    ::    PCCO25
      REAL_B    ::    PCLFR5(JP_LON,JP_LEV) 
      REAL_B    ::      PDP5(JP_LON,JP_LEV)
      REAL_B    ::    PEMIS5(JP_LON)
      REAL_B    ::    PEMIW5(JP_LON)
      REAL_B    ::     PLSM5(JP_LON)
      REAL_B    ::     PMU05(JP_LON)
      REAL_B    ::    POZON5(JP_LON,JP_LEV)
      REAL_B    ::       PQ5(JP_LON,JP_LEV)
!
      REAL_B    ::    PQIWP5(JP_LON,JP_LEV) 
      REAL_B    ::    PQLWP5(JP_LON,JP_LEV)          ! Dropplets    Concentration
      REAL_B    ::    PSQIW5(JP_LON,JP_LEV)          ! Ice Crystals Concentration
      REAL_B    ::    PSQLW5(JP_LON,JP_LEV)          ! 
      REAL_B    ::      PQS5(JP_LON,JP_LEV)
      REAL_B    ::   PQRAIN5(JP_LON,JP_LEV)
      REAL_B    ::   PRAINT5(JP_LON,JP_LEV)
      REAL_B    ::   PRLVRI5(JP_LON,JP_LEV)
      REAL_B    ::   PRLVRL5(JP_LON,JP_LEV)
      REAL_B    ::      PTH5(JP_LON,JP_LEVP1) 
      REAL_B    ::       PT5(JP_LON,JP_LEV) 
      REAL_B    ::      PTS5(JP_LON)
      REAL_B    ::    PNBAS5(JP_LON)        
      REAL_B    ::    PNTOP5(JP_LON)

!   For Use in SUCLD
!   ^^^^^^^^^^^^^^^^
      REAL_B    ::    ZETA(JP_LEV)      
      REAL_B    ::   ZETAH(JP_LEVP1)

!   For Use in SUOVLP
!   ^^^^^^^^^^^^^^^^^
      REAL_B    ::   ZTVIR              
      REAL_B    ::   ZFACT
      REAL_B    ::     ZAZ(JP_LEV)      
      REAL_B    ::    ZAZH(JP_LEVP1)



!  Local  Variables
!  ----------------

      REAL_B    ::    RTIMTR

      INTEGER_M ::    JL      , JK      , JAER    ,JNU
      INTEGER_M ::    KULOUT  , NINDAT  , NSSSSS  ,KPRTLEV
      INTEGER_M ::    IYR     , MONTH   , IDAY    ,IMINUT
      INTEGER_M ::    KPRINT  , SKSHIFT




! Load External Functions
! =======================

#include "fctast.h"
#include "fcttim.h"
#include "fcttre.h"


      KPRTLEV= 1
      KPRINT = 1
      SKSHIFT= 0




! Time Base
! =========

      NINDAT = min(yearTU,2004)*10000+mon_TU*100+Day_TU                 ! Date   in the form  yyyyMMdd
      NSSSSS =     HourTU      * 3600+minuTU* 60+sec_TU                 ! Nb of second since day Begin
      IYR    =   NAA(NINDAT)
      MONTH  =   NMM(NINDAT)
      IDAY   =   NDD(NINDAT)
      RTIMTR = RTIME(IYR,MONTH,IDAY,NSSSSS)
      IMINUT =   INT(FLOAT(NSSSSS)/60.)




! Basic Initialization
! ====================

!  Dimensions (auxiliary variables)
!  --------------------------------

        KIDIA  = 1                       ! DO'NT CHANGE
        KFDIA  = JP_LON                  ! Nb Columns
        KTDIA  = 1                       !
        KLON   = JP_LON                  ! Nb Columns
        KLEV   = JP_LEV                  ! Nb Levels
        KMODE  = JP_MODE                 ! Used in Planck Fcts Specification
        KAER   = JP_AER                  !

!   Nb of Solar Spectral Intervals
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        KSW    =  JP_SW       ! SW Nb of Spectral Intervals (max is JP_SW=6)
        NSW    =  JP_SW       ! SW Nb of Spectral Intervals (max is JP_SW=6)
        NTSW   =  JP_SW       ! SW Nb of Spectral Intervals (max is JP_SW=6)

        KBOX   = 0                       !                                      \VER
        NBOX   = 1                       !                                      \VER

        ILWRAD = 1       ! 0: Morcrette,     1991 operational before 20000627
                         ! 1: Mlawer et al., 1997 now ECMWF-operational
                         ! 2: Morcrette,     1991 original as in ERA-15'

        NDUMP  = 3       ! No Print
!       NDUMP  = 2       ! 1D Results
!       NDUMP  = 1       ! Debug
!       NDUMP  = 0       ! ALL
        KULOUT = 6       ! Output Device for SUCST


!  Verification of the Dimensions
!  ------------------------------

        write(6,*) 'INITIALISATION OF ECMWF RADIATIVE TRANSFERT: BEGIN'

        write(6,*) 'CONTROL'
        write(6,*) 'JP_LON=',JP_LON
        write(6,*) 'kcolp',kcolp
        IF    (kcolp .GT.JP_LON)                                    THEN
          write(6,6001) kcolp ,JP_LON
 6001     format(' @!#& BAD dimensions SET-UP  (kcolp > JP_LON);  = ('&
     &                                            ,i3,',',i3,')')
          STOP
        END IF
        IF    (mzp   .NE.JP_LEV)                                    THEN
          write(6,6002) mzp   ,JP_LEV
 6002     format(' @!#& BAD dimensions SET-UP  (mzp  /= JP_LEV);  = ('&
     &                                            ,i3,',',i3,')')
          STOP
        END IF
        IF    (naero .NE.JP_AER)                                    THEN
          write(6,6003) naero ,JP_AER
 6003     format(' @!#& BAD dimensions SET-UP  (naero/= JP_AER);  = ('&
     &                                            ,i3,',',i3,')')
          STOP
        END IF


!  Mathematical Constants
!  ----------------------

        REPLOG = 1.E-12                  ! Minimum Logarithm Argument
        REPSC  = 1.E-12
        REPSCO = 1.E-12
        REPSCQ = 1.E-12
        REPSCT = 1.E-12
        REPSCW = 1.E-12


!  Switches (general)
!  ------------------

        LONEWSW= .TRUE.  ! .TRUE. SWSN radiative routine is     used
       IF (ILWRAD.EQ.1)     THEN
        LRRTM  = .TRUE.  ! .TRUE. RRTM radiative routine is     used
       ELSE
        LRRTM  = .FALSE. ! .FALSE.RRTM radiative routine is NOT used
       END IF
        LTEMPDS= .FALSE. ! .TRUE. ALLOWS FOR SURF. T DISCONTIN. IN RAD.COMPUT.

        NTRAER = 19      ! NUMBER OF TRANSMISSION FUNCTIONS  W OR W/O AEROSOLS


!  Switches (Clouds Optical Properties)
!  ------------------------------------

        LINHOM = .FALSE. ! Tiedke (1995) correct. factor (0.7) of tau not used
        NHOWINH=  2      ! Tau correction factor:           exp(-(sig/tau)^2)
                         !  (used if LINHOM = .TRUE.) 
        LOIFUEC= .FALSE. ! .FALSE. IF ICE   CLOUDS AS EBERT-CURRY (LW & SW)
        LOWASYF= .FALSE. ! .FALSE. IF WATER CLOUDS AS FOUQUART         (SW)
        LOWHSSS= .FALSE. ! .FALSE. IF WATER CLOUDS AS SMITH-SHI   (LW)
        NRADIP =  3      ! Ice    effective Radius:
                         !  0   fixed        40 microns
                         !  1   f(T)   40 - 130 microns
                         !  2   f(T)   30 -  60 microns Jakob-Klein
                         !  3   f(T,IWC)                Sun-Rikus,     1999
        RMINICE=  15     ! Minimum Diameter for Ice Particles (micronm)
                         ! Needed only if   NRADIP = 3
                         ! (see von Walden et al., 2003 (Oct) Tab. 2 p.1393)
        NRADLP =  0      ! Liquid effective Radius: f(Pressure)
                         !  0 effective radius - liquid as f(Pressure)
                         !  1 fixed 10 microns over land, 13 over ocean
                         !  2 computed from LWC          Martin et al, 1994
        NLIQOPT=  1      ! Cloud Optical Properties (Water): 1=ECMWF Operat.
                         !  0  LW: Smith-Shi,   1992; SW: Fouquart,    1987
                         !  1  LW: Savijarvi,   1997; SW: Slingo  ,    1989
                         !  2  LW: Lindner,Li,  2000; SW: Slingo  ,    1989
        NICEOPT=  1      ! Cloud Optical Properties (Ice)  : 1=ECMWF Operat.
                         !  0  LW: Smith,Shi  , 1992; SW: Ebert-Curry, 1992
                         !  1  LW: Ebert,Curry, 1992; SW: Ebert-Curry, 1992
                         !  2  LW: Fu,Liou    , 1993; SW: Fu,Liou    , 1993
                         !  3  LW: Fu et al.  , 1998; SW: Fu         , 1996
        NOVLP  =  2      ! CLOUD OVERLAP CONFIGURATION:
                         !  1=MRN, 2=MAX, 3=RAN, 4=Hogan


!  Switches (Aerosols/O3)
!  ----------------------

        LNEWAER= .TRUE.  ! Climatology of Aerosols: TEGEN ET AL. 1997 / GADS
        LHVOLCA= .TRUE.  ! .TRUE. IF GISS HISTORY OF VOLCANIC AEROSOLS IS ON
!       NOZOCL = -1      ! TESTS the vertical quadrature       (NO absorber)
!       NOZOCL =  0      ! whatever is read for O3 as input profile
!       NOZOCL =  1      ! OLD         ECMWF O3 climatology and    aerosols
        NOZOCL =  2      ! Fortuin-Langematz O3 climatology and    aerosols
!       NOZOCL =  3      ! OLD         ECMWF O3 climatology and NO aerosols
!       NOZOCL =  4      ! Fortuin-Langematz O3 climatology and NO aerosols
 

!  BASIC CONSTANTS
!  ---------------

!            *****
        CALL SUCST (KULOUT, NINDAT, NSSSSS, KPRTLEV) ! Initialize common YOMCST
!            *****                                   !  (Basic Constants)
                                                     ! Initialize common YOMRIP
                                                     !  (only date and time)

!   YOENCST - THERMODYNAMIC TRANSITION OF PHASE
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        RTT=273.16                                   !

!   YOETHF  - DERIVED CONSTANTS SPECIFIC TO ECMWF THERMODYNAMICS
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        RTWAT=RTT                                    ! 
        RTICE=RTT-23.                                !

!   YOERDU  - CONTROL, PARAMETERS AND SECURITY IN RADIATION
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        RCDAY  = RDAY * RG / RCPD                    !
        R10E   = 0.4342945                           ! DECIMAL /  NATURAL
                                                     ! LOG.        FACTOR
        DIFF   = 1.66                                ! DIFFUSIVITY FACTOR


!  Space/Time Independant Coefficients
!  -----------------------------------

        CALL SURDI            ! ECMWF Surface Albedo, Emissivity
        CALL SULWN            ! Initialize common YOELW (new LW Coeff.)
        CALL SUOLW            ! Initialize common YOELW (old LW Coeff.)

!   Initialization routine for RRTM
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        CALL SURRTAB          ! AER'S RRTM LW RADIATION
        CALL SURRTPK          ! Initialize common YOERRTWN 
                              !  (k-coefficients in spectral intervals)
        CALL SURRTRF          ! Initialize common YOERRTRF
                              !  (RRTM Reference Atmosphere)
        CALL SURRTFTR         ! Initialize common YOERRTRF

!            RRTM routine     ! BAND   [cm-1] ! low           ! high
!       ----------------------+---------------+---------------+---------
        CALL RRTM_KGB1        !  1:   10- 250 ! H2O           ! H2O
        CALL RRTM_KGB2        !  2:  250- 500 ! H2O           ! H2O
        CALL RRTM_KGB3        !  3:  500- 630 ! H2O,CO2       ! H2O,CO2
        CALL RRTM_KGB4        !  4:  630- 700 ! H2O,CO2       ! O3,CO2
        CALL RRTM_KGB5        !  5:  700- 820 ! H2O,CO2       ! O3,CO2
        CALL RRTM_KGB6        !  6:  820- 980 ! H2O           ! nothing
        CALL RRTM_KGB7        !  7:  980-1080 ! H2O,O3        ! O3
        CALL RRTM_KGB8        !  8: 1080-1180 ! (i.e.>~300mb) ! O3
                              !               ! H2O           !
        CALL RRTM_KGB9        !  9: 1180-1390 ! H2O,CH4       ! CH4
        CALL RRTM_KGB10       ! 10: 1390-1480 ! H2O           ! H2O
        CALL RRTM_KGB11       ! 11: 1480-1800 ! H2O           ! H2O
        CALL RRTM_KGB12       ! 12: 1800-2080 ! H2O,CO2       ! nothing
        CALL RRTM_KGB13       ! 13: 2080-2250 ! H2O,N2O       ! nothing
        CALL RRTM_KGB14       ! 14: 2250-2380 ! CO2           ! CO2
        CALL RRTM_KGB15       ! 15: 2380-2600 ! N2O,CO2       ! nothing
        CALL RRTM_KGB16       ! 16: 2600-3000 ! H2O,CH4       ! nothing

!   Reduce absorption coefficient data from 256 to 140 g-points
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        CALL RRTM_INIT_140GP

!   Initialization routine for SW (6 spectral interval resolution)
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        CALL SUSWN ( NTSW  , KSW  )      ! Initialize common YOESW


        write(6,*) 'INITIALISATION OF ECMWF RADIATIVE TRANSFERT: END  '




! Radiation: Global (Time dependant) Parameters
! =============================================

        PRII05 = RI0                                 ! INSOLATION
!       PRII05 = RI0/(dST_UA*dST_UA)                 ! INSOLATION   (dST_UA: Distance Soleil-Terre [UA])
        PCCO25 = 360.E-06*44./29.                    ! CONCENTRATION IN CO2                     [Pa/Pa]




! Surface Properties
! ==================

!  Land/sea Mask
!  -------------

! Martin CONTROL
!PRINT*, 'CONTROL dans Atm_RT_INI'
!PRINT*,'SIZE(psa_DY)=)',SIZE(psa_DY)
!PRINT*, 'psa_DY=',psa_DY
! Martin CONTROL

  DO jl=1,KLON
        ikl               = min(jl,kcolp)
        i                 = ii__AP(ikl)
        j                 = jj__AP(ikl)
        PLSM5 (JL)    = 1 - MaskSV_gpt(ikl)   




! Atmospheric Thermodynamics (Time and Space dependant)
! =====================================================

!  Pressure      Distribution
!  --------------------------

       JK=1+KLEV
        PAPH5 (JL,JK)     = (psa_DY(ikl)             + pt__DY) * 1000.  ! Pressure (Layer Interface)[Pa]
    DO JK=1,KLEV
        PAPH5 (JL,JK)     = (psa_DY(ikl) * sigma(JK) + pt__DY) * 1000.  ! Pressure (Layer)          [Pa]
        PAP5  (JL,JK)     = (psa_DY(ikl) * sigmi(JK) + pt__DY) * 1000.  ! Pressure (Layer Interface)[Pa]
        PDP5  (JL,JK)     =  psa_DY(ikl) *dsigmi(JK)           * 1000.  ! Pressure (Layer Thickness)[Pa]


!  Temperature   Distribution
!  --------------------------

        PT5   (JL,JK)     =  Ta__DY(ikl,JK)


!  Water Species Distribution
!  --------------------------

        PQ5   (JL,JK)     =  qv__DY(ikl,JK)                             ! Water Vapor  Concentr. [kg/kg]

    ENDDO
  ENDDO




! Initialization (Climatologies, Time independant)
! ================================================

!  Aerosols Radiative Characteristics (YOEAER)
!  ----------------------------------

!            *******
        CALL SUAERL                  ! Aerosols LW Radiative Charact.
        CALL SUAERSN ( NTSW , KSW )  ! Aerosols SW Radiative Charact.
!            *******


!  Aerosols Optical Thickness Horizontal Distribution (model grid
!  --------------------------------------------------  independant)

!            ********
        CALL SUAERH                  ! 
        CALL SUECAEBC                ! BLACK CARBON (URBAN/FOREST FIRE ORIGIN)
        CALL SUECAEOR                ! ORGANIC-TYPE
        CALL SUECAESD                ! SOIL-DUST                       ORIGIN
        CALL SUECAESS                ! SEA -SALT                       ORIGIN
        CALL SUECAESU                ! SULFATE-TYPE
!            ********


!  Clouds (YOECLD)
!  ---------------

        DO jk=1,klev
          ZETA(jk) =PAP5(1,jk) /PAPH5(1,klev+1)
        ENDDO
        DO jk=1,klev+1
          ZETAH(jk)=PAPH5(1,jk)/PAPH5(1,klev+1)
        ENDDO
! Martin Control
!PRINT*, 'PAPH5(1,:)=',PAPH5(1,:)
!PRINT*, 'ZETAH=',ZETAH
! Martin Control

!            *****
        CALL SUCLD  ( KLEV  , ZETA ) 
!            *****


!  Cloud Optical Parameters SW/LW (all parameterizations)
!  ------------------------------------------------------

!            *******
        CALL SUCLOPN ( NTSW , KSW , KLEV )  ! Initialize YOECLOP
!            *******


!  Radar Reflectivity
!  ------------------

         ZAZH(KLEV+1)=     0.
          ZAZ(1)     =100000.
           JL=KIDIA
        DO jk=KLEV,2,-1
          ZTVIR      =       PT5(JL,jk)   /(1.-0.608*PQ5(jl,jk))
          ZFACT      = LOG(PAPH5(jl,jk+1))-    LOG(PAPH5(jl,jk))
          ZAZH(jk)   =      ZAZH(jk+1)    + R *    ZTVIR/(RMD*RG)*ZFACT
           ZAZ(jk)   = 0.5*(ZAZH(jk+1)+ZAZH(jk))*1000.
        END DO

!            ******
        CALL SUOVLP  ( KLEV , ZAZ )         ! Initialize ALPHA1 (%radar refl.)
!            ******                         !  (Hogan & Illingsworth, 1999)


!  NO Absorber
!  -----------

        IF (NOZOCL.EQ.-1) THEN
                RCH4             = 1.E-18
                RN2O             = 1.E-18
                RO3              = 1.E-18
                RCFC11           = 1.E-18
                RCFC12           = 1.E-18
                PCCO25           = 1.E-18
          DO jk=1,klev
              DO jl=KIDIA,KFDIA
                POZON5(JL,JK)    = 0.
              ENDDO
          ENDDO
          DO JK=1,KLEV
            DO JAER=1,KAER
              DO JL=KIDIA,KFDIA
                PAER5(JL,JAER,JK)= ZEPAER
              ENDDO
            ENDDO
          ENDDO
        END IF



      return
      end subroutine Atm_RT_INI
