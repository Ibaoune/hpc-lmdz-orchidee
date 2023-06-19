      subroutine Atm_RT_RUN(rklonr,iklOUT)

!--------------------------------------------------------------------------+
!                                                     Mon 24-Jun-2013  MAR |
!     subroutine Atm_RT_RUN interfaces MAR PHYsics                         |
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
      use Mod_PHY____dat
      use Mod_PHY_RT_dat
      use Mod_PHY____grd
      use Mod_PHY_S0_grd
      use Mod_PHY_RT_grd
      use Mod_PHY_S0_kkl
      use Mod_PHY_RT_kkl
      use Mod_PHY_DY_kkl
      use Mod_PHY_CM_kkl
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

      logical                      ::  RTi_qi_and_qs = .TRUE.           ! * as Very small Cloud Particle: T
!Gilles: T > F to have interactive cloud frac. & radiation
      logical                      ::  CFrCEP        = .FALSE.          ! ECMWF or CMiPhy Cloud Fraction: T OR F



!  INPUT (from MAR/ECMWF Interface)
!  -----

      integer                      ::  rklonr                           ! nb    des pts de grilles
      integer                      ::  iklOUT                           ! OUTPUT   (pt  de grille)

      integer                      ::  kj                               ! Index des niveaux
      integer                      ::  nAe                              ! Index d'aerosols

! #DB integer                      ::  kio
! #DB integer                      ::  kjo
      integer, dimension(rklonr)   ::  k2ii                             !
      integer, dimension(rklonr)   ::  k2jj                             !
! #DB integer                      ::  jkjllw      
! #DB integer                      ::  lijio



!  INTERNAL VARIABLES
!  ------------------

!   For Use in radlsw
!   ^^^^^^^^^^^^^^^^^
      REAL_B    ::   PGELAM5(JP_LON)
      REAL_B    ::    PGEMU5(JP_LON)
      REAL_B    ::    PSLON5(JP_LON)
      REAL_B    ::    PCLON5(JP_LON)
      REAL_B    ::    ZOZON5(JP_LON,JP_LEV)
      REAL_B    ::     ZAER5(JP_LON,JP_AER,JP_LEV)
!
!
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
!
      REAL_B    ::    PEMIT5(JP_LON)
      REAL_B    ::     PFCT5(JP_LON,JP_LEVP1)
      REAL_B    ::     PFLT5(JP_LON,JP_LEVP1)
      REAL_B    ::     PFCS5(JP_LON,JP_LEVP1)
      REAL_B    ::     PFLS5(JP_LON,JP_LEVP1)
      REAL_B    ::   PFRSOD5(JP_LON)
      REAL_B    ::    PSUDU5(JP_LON)
      REAL_B    ::    PUVDF5(JP_LON)
      REAL_B    ::    PPARF5(JP_LON)
      REAL_B    ::    PFDCT5(JP_LON,JP_LEVP1)
      REAL_B    ::    PFUCT5(JP_LON,JP_LEVP1)
      REAL_B    ::    PFDLT5(JP_LON,JP_LEVP1)
      REAL_B    ::    PFULT5(JP_LON,JP_LEVP1)
      REAL_B    ::    PFDCS5(JP_LON,JP_LEVP1)
      REAL_B    ::    PFUCS5(JP_LON,JP_LEVP1)
      REAL_B    ::    PFDLS5(JP_LON,JP_LEVP1)
      REAL_B    ::    PFULS5(JP_LON,JP_LEVP1)
!
      REAL_B    ::    ZTAU5 (JP_LON,JP_SW,JP_LEV)    ! Cloud Optical Depth
      REAL_B    ::    ZTAUI5(JP_LON)                 ! Cloud Optical Depth (vert.int.)
!
      REAL_B    ::    ASWBOX(JP_LON,100)
      REAL_B    ::    OLRBOX(JP_LON,100)
      REAL_B    ::    SLWBOX(JP_LON,100)
      REAL_B    ::    SSWBOX(JP_LON,100)
      REAL_B    ::    TAUBOX(JP_LON,100)
      REAL_B    ::    CLDBOX(JP_LON,100,JP_LEV)

!   For Use in SUCLD
!   ^^^^^^^^^^^^^^^^
      REAL_B    ::   ZETAH(JP_LEVP1)



!  Local  Variables
!  ----------------

      REAL_B    ::         RTIMTR  , ZTHETOZ , ZANGOZC
      REAL_B    ::         Zone_t                                       !  Time  Zone                               [hr]
      REAL_B    ::         qcloud                                       !  Cloud Particles Concentration          [g/kg]
      REAL_B    ::         fcloud                                       !  Cloud           Fraction                  [%]
      REAL_B    ::         LWUpwd                                       !  Surface Upward  Longwave Radiation     [W/m2]
      REAL_B    ::         Emi_Cl                                       !  Cloud Emissivity (diagnostic)             [-]

      INTEGER_M ::    i       , j       , ikl     ,k
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




! OUTPUT for DEBUGGING
! ====================

! #DB kio = 1
! #DB kjo = 1




! Date & Time
! ===========

      NINDAT = min(yearTU,2004)*10000+mon_TU*100+Day_TU                 ! Date   in the form  yyyyMMdd
      NSSSSS =     HourTU      * 3600+minuTU* 60+sec_TU                 ! Nb of second since day Begin

! VER   v
        write(6,*) ' NINDAT: ',NINDAT
        write(6,*) ' NSSSSS: ',NSSSSS
! VER   ^
      IYR    =   NAA(NINDAT)
      MONTH  =   NMM(NINDAT)
      IDAY   =   NDD(NINDAT)
      RTIMTR = RTIME(IYR,MONTH,IDAY,NSSSSS)
      IMINUT =   INT(FLOAT(NSSSSS)/60.)




! Basic Initialization
! ====================

!  Dimensions (auxiliary variables)
!  --------------------------------

        KIDIA  = 1                                   ! DO'NT CHANGE
        KFDIA  = JP_LON                              ! Nb Columns
        KTDIA  = 1                                   !
        KLON   = JP_LON                              ! Nb Columns
        KLEV   = JP_LEV                              ! Nb Levels
        KMODE  = JP_MODE                             ! Used in Planck Fcts Specification
        KAER   = JP_AER                              !

        DO     ikl  = 1 ,  kcolp
          k2ii(ikl) = ii__AP(ikl)
          k2jj(ikl) = jj__AP(ikl)
        ENDDO

!   Nb of Solar Spectral Intervals
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        KSW    =  JP_SW                              ! SW Nb of Spectral Intervals (max is JP_SW=6)
        NSW    =  JP_SW                              ! SW Nb of Spectral Intervals (max is JP_SW=6)
        NTSW   =  JP_SW                              ! SW Nb of Spectral Intervals (max is JP_SW=6)

        KBOX   = 0                                   !                                      \VER
        NBOX   = 1                                   !                                      \VER

        ILWRAD = 1                                   ! 0: Morcrette,     1991 operational before 20000627
                                                     ! 1: Mlawer et al., 1997 now ECMWF-operational
                                                     ! 2: Morcrette,     1991 original as in ERA-15'

        NDUMP  = 3                                   ! No Print
!       NDUMP  = 2                                   ! 1D Results
!       NDUMP  = 1                                   ! Debug
!       NDUMP  = 0                                   ! ALL
        KULOUT = 6                                   ! Output Device for SUCST




! Radiation: Global (Time dependant) Parameters
! =============================================

        PRII05 = RI0/(dST_UA*dST_UA)                 ! INSOLATION   (dST_UA: Distance Soleil-Terre [UA])
        PCCO25 = 360.E-06*44./29.                    ! CONCENTRATION IN CO2 (PA/PA)




! Surface Properties
! ==================

!  Surface Albedo
!  --------------

  DO jl=1,KLON
        ikl               = min(jl,kcolp)
        i                 = ii__AP(ikl)
        j                 = jj__AP(ikl)
    DO jnu=1,KSW
        PALBD5(JL,JNU)    = Alb_SV_gpt(jl)
        PALBP5(JL,JNU)    = Alb_SV_gpt(jl)
    ENDDO

!  Surface Emissivity
!  ------------------

        PEMIS5(JL)        = EmisSV_gpt(jl)
        PEMIW5(JL)        = EmisSV_gpt(jl)


!  Land/sea Mask (1.=land)
!  -------------

        PLSM5 (JL)    = 1 - MaskSV_gpt(jl)


!  Cosine (Solar zenithal Distance)
!  --------------------------------

        PMU05 (JL)    =     csz_S0(ikl)




! Atmospheric Thermodynamics (Time and Space dependant)
! =====================================================

!  Pressure
!  --------

! Martin control
!PRINT*,'sigma(:)=',sigma(:)
!PRINT*,'sigmi(:)=',sigmi(:)
! Martin control

       JK=1+KLEV
        PAPH5 (JL,JK)     = (psa_DY(ikl)             + pt__DY) * 1000.  ! Pressure (Layer Interface)[Pa]
    DO JK=1,KLEV
        PAPH5 (JL,JK)     = (psa_DY(ikl) * sigma(JK) + pt__DY) * 1000.  ! Pressure (Layer)          [Pa]
        PAP5  (JL,JK)     = (psa_DY(ikl) * sigmi(JK) + pt__DY) * 1000.  ! Pressure (Layer Interface)[Pa]
        PDP5  (JL,JK)     =  psa_DY(ikl) *dsigmi(JK)           * 1000.  ! Pressure (Layer Thickness)[Pa]


!  Water Species      Distributions
!  --------------------------------

          PQ5    (JL,JK)  = qv__DY(ikl,JK)                              ! Water Vapor  Concentr. [kg/kg]
        IF      (RTi_qi_and_qs)                                     THEN
          PQIWP5 (JL,JK)  = qi__CM(ikl,JK)                             &! Ice Crystals Concentr. [kg/kg]
     &    + (1.-min(1.,exp((Ta__DY(ikl,JK)-258.15)*0.1)))              &! AND
     &    * (               qs__CM(ikl,JK)        *0.33 )               ! Snow Particl.Concentr. [kg/kg]
        ELSE                                                            !
          PQIWP5 (JL,JK)  = qi__CM(ikl,JK)                              ! Ice Crystals Concentr. [kg/kg]
        END IF
          PQIWP5 (JL,JK)    = max(0.,PQIWP5(JL, JK)-1.E-9)              ! Ice Crystals Concentr. [kg/kg]

          PQLWP5 (JL,JK)    = max(0.,qw__CM(ikl,JK)-1.E-9)              ! Dropplet     Concentr. [kg/kg]

          PQS5   (JL,JK)    =        qvswCM(ikl,JK)                     ! Saturat. % Water       [kg/kg]
          PQRAIN5(JL,JK)    =        qr__CM(ikl,JK)                     ! Drops        Concentr. [kg/kg]
          PRAINT5(JL,JK)    = 0.                                        !                           \VER
          PRLVRI5(JL,JK)    = 0.                                        ! e-mail J.-J.M. 20031203
          PSQIW5 (JL,JK)    = 1.                                        !    exp(-PRLVRI5(JL,JK))
          PRLVRL5(JL,JK)    = 0.                                        ! e-mail J.-J.M. 20031203
          PSQLW5 (JL,JK)    = 1.                                        !    exp(-PRLVRL5(JL,JK))


!  Cloud Fraction
!  --------------

        IF (CFrCEP)                                                 THEN
          PCLFR5 (JL,JK)    =       (PQIWP5(JL,JK)+PQLWP5(JL,JK))      &! ECMWF Paramet.
     &                             /(FraQws       *  PQS5(JL,JK))       !  (VERY Crude)
          PCLFR5 (JL,JK)    =       min( _ONE_    ,PCLFR5(JL,JK))       !
          PCLFR5 (JL,JK)    =       max(1.0E-3    ,PCLFR5(JL,JK))      &! no small values
     &        *max(_ZERO_,sign(_ONE_,PQIWP5(JL,JK)+PQLWP5(JL,JK)       &!
     &                              -2.E-9_JPRB))                       !
          CFraCM (ikl,JK)   =        PCLFR5(JL,JK)                      !
        ELSE
          PCLFR5 (JL,JK)    =        CFraCM(ikl,JK)                     ! Cloud Fraction from CMiPhy [-]
        END IF
!         write(4,4) JL,JK,PQLWP5(JL,JK),PQIWP5(JL,JK)                 &
!    &                    ,  PQS5(JL,JK)                               &
!    &                    ,CFraCM(JL,JK),PCLFR5(JL,JK)
!4        format(2i6,' Cloud Liq.W.= ', f9.6,5x                        &
!    &              ,' Cloud Sol.W.= ', f9.6,5x                        &
!    &              ,' Satur.W.Vap.= ', f9.6,5x                        &
!    &              ,' Cloud Fract.= ',2f9.6)

!  Temperature        Distribution
!  -------------------------------

          PT5    (JL,JK)    =        Ta__DY(ikl,JK)
          PTH5   (JL,JK)    = 0.5 * (Ta__DY(ikl,JK)+Ta__DY(ikl,max(1,JK-1)))
    END DO
        PTH5   (JL,KLEV+1)  =        Ta__DY(ikl,mzpp)
        PTS5   (JL)         =        Ta__DY(ikl,mzpp)


!  Convective Layer
!  ----------------

        PNBAS5 (JL)       = 1.                !                                 \VER
        PNTOP5 (JL)       = 1.                !                                 \VER
  ENDDO




! Assignation    (Climatologies, Time   dependant)
! ================================================

!  Aerosols Optical Thickness Horizontal Distribution (model grid
!  --------------------------------------------------  independant)

!            *******
        CALL SUECAEC ( NINDAT, IMINUT )     ! TEGEN ET AL. (1997, JGR 102, 
!            *******                        !               pp23895-23915)
 

!  Aerosols Optical Thickness Vertical   Distribution
!  --------------------------------------------------
!            *******


! Martin modification for MAR-LMDZ:
DO jk=1,klev+1
   ZETAH(jk)=PAPH5(1,jk)/PAPH5(1,klev+1)
ENDDO
! Martin Control
!Print*, 'Dans Atm_RT_RUN:'
!PRINT*,'ZETAH=',ZETAH
!PRINT*,'=KLEV',KLEV
! Martin Control

! Martin Control

        CALL SUAERV                                                  &
                 & ( KLEV   ,ZETAH                                   &
                 & , CVDAES ,CVDAEL ,CVDAEU ,CVDAED                  &
                 & , RCTRBGA,RCVOBGA,RCSTBGA,RCAEOPS,RCAEOPL,RCAEOPU &
                 & , RCAEOPD,RCTRPT ,RCAEADK,RCAEADM,RCAEROS         &
                 & )
!            *******


!  O3
!  --   

!   Fortuin-Langematz O3 climatology
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!            *******
        CALL SUECOZC ( NINDAT , IMINUT )
!            *******

!   ECMWF   Geleyn    O3 Climatology
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ZTHETOZ=RTETA( RTIMTR  )
        ZANGOZC=  REL( ZTHETOZ ) - 1.7535

!            *******
        CALL SUECOZO ( ZANGOZC )
!            *******




! Interpolation on the MAR Grid
! =============================

        DO JL=1,KLON
          ikl=min(JL,kcolp)
          PGELAM5(JL)=    lon__r(ikl)
           PGEMU5(JL)=SIN(lat__r(ikl))
           PCLON5(JL)=COS(lon__r(ikl))
           PSLON5(JL)=SIN(lon__r(ikl))
        END DO

!              ******
         CALL  RADACA                                       &
          &( KIDIA , KLON   , KLON  , KTDIA , KLEV          &
          &, PAPH5 , PGELAM5, PGEMU5, PCLON5, PSLON5, PTH5  &
          &, ZAER5 , ZOZON5                                 &
          &  )
!              ******


!  OLD         ECMWF O3 Climatology
!  --------------------------------

        IF (NOZOCL.EQ.1 .OR. NOZOCL.EQ.3)                           THEN
          DO JK=1,KLEV
            DO JL=KIDIA,KFDIA
              POZON5(JL,JK)      = ZOZON5(JL,JK)
            ENDDO
          ENDDO
        END IF


!  FORTUIN LANGEMATZ O3 Climatology
!  --------------------------------

        IF (NOZOCL.EQ.2 .OR. NOZOCL.EQ.4)                           THEN

!              ******
          CALL RADOZC ( KIDIA ,KLON  , KLON   , KTDIA, KLEV            &
                     &, KPRINT,KLON  , SKSHIFT, PAPH5, PGEMU5, ZOZON5)
!              ******

          DO JK=1,KLEV
            DO JL=KIDIA,KFDIA
              POZON5(JL,JK)      = ZOZON5(JL,JK)
            ENDDO
          ENDDO
 
        END IF


!  AEROSOLS
!  --------

        IF (NOZOCL.EQ.1 .OR. NOZOCL.EQ.2)                           THEN
          DO jk=1,klev
            DO jaer=1,KAER
              DO jl=KIDIA,KFDIA
                PAER5(JL,JAER,JK)= ZAER5(JL,JAER,JK)
              ENDDO
            ENDDO
          ENDDO
        END IF


!  NO AEROSOLS
!  -----------

        IF (NOZOCL.GT.2)                                            THEN
          DO jk=1,klev
            DO jaer=1,KAER
              DO jl=KIDIA,KFDIA
                PAER5(JL,JAER,JK)= ZEPAER
              ENDDO
            ENDDO
          ENDDO
        END IF


!  SECURITY CHECK ON AEROSOL AMOUNTS
!  ---------------------------------

        DO JK=1,KLEV
          DO JAER=1,KAER
            DO JL=KIDIA,KFDIA
              PAER5(JL,JAER,JK)=MAX(ZEPAER,PAER5(JL,JAER,JK))
            ENDDO
          ENDDO
        ENDDO




! Transmission to MAR Variables
! =============================

      DO JL=1,KLON
              ikl=min(JL,kcolp)

        DO JK=1,KLEV
              O3__RT(ikl,jk     ) = ZOZON5(JL,JK)                       ! O3      Concentr.
        END DO

        DO JK=1,KLEV
          DO JAER=1,KAER
              AersRT(ikl,jk,jaer) =  PAER5(JL,JAER,JK)                  ! Aerosol Optical Thickness
          END DO
        END DO

      END DO




! Solar and IR Transfer through the Atmosphere
! ============================================


! VER   v
!       STOP 'Fin momentanee'
! VER   ^
!                ******
            CALL RADLSW                                                 &
     &    ( KIDIA , KFDIA , KLON  , KTDIA , KLEV   , KMODE , KAER,      &
     &      KBOX  , NBOX                                                &
     &    , NDUMP , ILWRAD                                              &
     &    , PRII05                                                      &
     &    , PAER5 , PALBD5, PALBP5, PAPH5 , PAP5                        &
     &    , PCCO25, PCLFR5, PDP5  , PEMIS5, PEMIW5 , PLSM5 ,  PMU05,    &
     &      POZON5                                                      &
     &    , PQ5   , PQIWP5, PQLWP5, PSQIW5, PSQLW5 , PQS5  ,  PQRAIN5,  &
     &      PRAINT5                                                     &
     &    , PRLVRI5,PRLVRL5,PTH5  , PT5   , PTS5   , PNBAS5,  PNTOP5    &
     &    , PEMIT5, PFCT5 , PFLT5 , PFCS5 , PFLS5  , PFRSOD5, PSUDU5,   &
     &      PUVDF5, PPARF5                                              &
     &    , PFDCT5, PFUCT5, PFDLT5, PFULT5, PFDCS5 , PFUCS5,  PFDLS5,   &
     &      PFULS5                                                      &
     &    , ZTAU5 , ZTAUI5                                              &
     &    , ASWBOX, OLRBOX, SLWBOX, SSWBOX, TAUBOX , CLDBOX             &
! #DB&                                                     ,  k2ii,k2jj &
     &    )
!                ******




! Radiative Fluxes   Distributions
! ================================

  DO JL=1,KLON
        ikl = min(JL,kcolp)
    DO JK=1,KLEV+1
        FIRn_c (ikl,jk)    = PFCT5 (JL,JK)                              ! CLEAR-SKY LW NET FLUXES [W/m2]
        FIRn_t (ikl,jk)    = PFLT5 (JL,JK)                              ! TOTAL-SKY LW NET FLUXES [W/m2]
        FSOn_c (ikl,jk)    = PFCS5 (JL,JK)                              ! CLEAR-SKY SW NET FLUXES [W/m2]
        FSOn_t (ikl,jk)    = PFLS5 (JL,JK)                              ! TOTAL-SKY SW NET FLUXES [W/m2]
    END DO




! Cloud Optical Depth
! ===================

    DO JK=1,KLEV
        kj              =  KLEV + 1    -JK                              ! 
        ODCzRT (ikl,kj) =  ZTAU5(JL,  1,JK)                             ! Cloud Optical Depth
       DO nAe=1,naero                                                   !
        ODAzRT (ikl,jk) = ODAzRT(ikl,jk   ) + PAER5(JL,nAe,JK)          ! Aeros.Optical Depth
       END DO                                                           !
    END DO                                                              !
                                                                        !
        ODC_RT (ikl   ) = ZTAUI5(JL)                                    ! Cloud Optical Depth (vert.integr., 
                                                                        !                      1st interval)
        ODA_RT (ikl   ) = 0.                                            !
    DO JK=1,KLEV                                                        !
        ODA_RT (ikl   ) = ODA_RT(ikl)       + ODAzRT(ikl,jk)            ! Aeros.Optical Depth
    END DO




! SURFACE RADIATIVE CHARACTERISTICS (SW)
! ======================================

        FSOs_t(ikl)     = PFRSOD5(JL)                                   ! TOTAL-SKY SRF SW DOWNWARD FLUX    [W/m2]
        FSOdir(ikl)     =  PSUDU5(JL)                                   ! SOLAR RADIANCE IN SUN'S DIRECT.   [W/m2]
        FSOsUV(ikl)     =  PUVDF5(JL)                                   ! SURFAC.DOWNWARD U.V. RADIATION    [W/m2]
        FSOeff(ikl)     =  PPARF5(JL)                                   ! PHOTOSYNTHET. ACTIVE RADIATION    [W/m2]




! SURFACE RADIATIVE CHARACTERISTICS (LW)
! ======================================
!       EmisSV_gpt(jl ) =  PEMIT5(JL)                                   ! TOTAL LW EMISSIVITY                \VER

  END DO




!  Grid  Point   Dependant Variables <-- Atm_RT "Vector"Variables
!  ==============================================================

        DO   ikl = 1,kcolp

              i =  ii__AP(ikl)
              j =  jj__AP(ikl)

!  OutgoingLongWave Radiation Fluxes
!  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              OLR_RT    (ikl)   =-FIRn_t(ikl,1   )                                            ! Atm Top LongWave Heat Flux (+)  (  Upward)

!  Surface Downward Radiative Fluxes
!  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              SWDsRT    (ikl)   = FSOs_t(ikl     )                                            ! Surface ShrtWave Heat Flux (+)  (Downward)
              SWAsRT    (ikl)   = FSOs_t(ikl     )*(1.-Alb_SV_gpt(ikl))                       ! Surface ShrtWave Heat Flux (+)  (Absorbed)
              LWDsRT    (ikl)   = FIRn_t(ikl,mzpp)                                           &!
     &                          + Ta__DY(ikl,mzpp)*    Ta__DY(ikl,mzpp)                      &!
     &                           *Ta__DY(ikl,mzpp)*    Ta__DY(ikl,mzpp)                      &!
     &                           *StefBo          *    EmisSV_gpt(ikl)                        !

!  Cloud   Fraction
!  ^^^^^^^^^^^^^^^^
              ClouRT    (ikl)   =  0.0
          DO  k  = 1,mzp
              ClouRT    (ikl)   =  max(CFraCM(ikl,k)   , ClouRT(ikl))

!  Radiative Heating
!  ^^^^^^^^^^^^^^^^^
              LWdTRT    (ikl,k) =  - ( FIRn_t(ikl,k+1) - FIRn_t(ikl,k) )                     &!
     &                  *  Grav_F  / ( CpdAir *1.e3 *psa_DY(ikl) *dsigmi(k))                  !
              SWdTRT    (ikl,k) =  - ( FSOn_t(ikl,k+1) - FSOn_t(ikl,k) )                     &!
     &                  *  Grav_F  / ( CpdAir *1.e3 *psa_DY(ikl) *dsigmi(k))                  !

              dpktRT    (ikl,k) =    ( LWdTRT(ikl,k) + SWdTRT(ikl,k) ) / ExnrDY(ikl,k)
              LWdTRT    (ikl,k) =      LWdTRT(ikl,k) * 86400.
              SWdTRT    (ikl,k) =      SWdTRT(ikl,k) * 86400.

          END DO

        END DO





! OUTPUT
! ======

! OUTPUT for Analysis
! -------------------

      IF          (iklOUT.GT.0)                                     THEN

        ikl =      iklOUT
        i = ii__AP(iklOUT)
        j = jj__AP(iklOUT)
        Zone_t = lon__h(ikl) - 12.

            write(4,401)Day_TU,Mon_TU,HourTU       ,minuTU,            &
     &                                Zone_t       , i ,j
  401       format(//,' +-- Radiative  Heat Fluxes --+',               &
     &       i4,'/',i2,i4,':',i2,' h.LT (',f4.0')',                    &
     &                '  (i,j) = (',i3,',',i3,')',                     &
     &        //,' | Altitud | Pressur.| Temper.|  Ozone  | W.Vapor|', &
     &      ' Clouds | Clouds |  Opt.D | Aer.OD. | So Warm.|',         &
     &       ' Emiss.| IR Cool.| IR NetF.|',                           &
     &         /,' |   [km]  |   [hPa] |   [K]  | [cmSTP] | [g/kg] |', &
     &      ' [g/kg] |   [%]  |   [-]  |   [-]   | [K/day] |',         &
     &       '  [-]  | [K/day] |  [W/m2] |',                           &
     &         /,' +---------+---------+--------+---------+--------+', &
     &              '--------+--------+--------+---------+---------+', &
     &       '-------+---------+---------+')

            write(4,404)                           pt__DY*1.e+1        &
     &           ,                                 OLR_RT(ikl)
  404       format(                                                    &
     &                     ' |  SOMMET |',f8.2,' |', 7x ,' |', 8x ,' |'&
     &        , 3( 7x ,' |'),    7x ,' |', 8x ,' |', 8x ,' |'          &
     &        ,                            6x ,' |', 8x ,' |',f8.1,' |'&
     &        ,/,' +---------+---------+--------+---------+--------+', &
     &              '--------+--------+--------+---------+---------+', &
     &       '-------+---------+---------+')


            DO k =1,mzp

               qcloud = 1.e+3 *(qw__CM(ikl,k)+qi__CM(ikl,k))
               fcloud = 1.e+2 * CFraCM(ikl,k)
               Emi_Cl = 1.-exp(-qcloud *0.5 *(Z___DY(ikl,max(1,k-1))-Z___DY(ikl,k+1))  &
     &                            * sigma(k)* psa_DY(ikl) /    (287.*Ta__DY(ikl,k  )))

               write(4,402)    1.e-3*Z___DY(ikl,k),PAP5  (ikl,k)*1.e-2 &
     &           ,Ta__DY(ikl,k)     ,O3__RT(ikl,k)                     &
     &           ,qv__DY(ikl,k)*1.e3,qcloud       ,fcloud              &
     &           ,ODCzRT(ikl,k)     ,ODAzRT(ikl,k),SWdTRT(ikl,k)       &
     &           ,Emi_Cl                          ,LWdTRT(ikl,k)       &
     &                                            ,FIRn_t(ikl,k)
  402          format(' | ',  f7.3,' |' ,f8.2,' |',f7.2,' |',e8.2,' |',&
     &         2(f7.3,' |'),2(f7.2,' |'),e8.2,' |',f8.4,' |'           &
     &          ,                         f6.3,' |',f8.2,' |',f8.1,' |')

            END DO


               LWUpwd =          Ta__DY(ikl,mzpp)*    Ta__DY(ikl,mzpp) &!
     &                          *Ta__DY(ikl,mzpp)*    Ta__DY(ikl,mzpp) &!
     &                          *StefBo          *    EmisSV_gpt( ikl)  !

            write(4,403)                     1.e+1*psa_DY(ikl)         &
     &           ,Ta__DY(ikl,mzpp)                                     &
     &           ,ODC_RT(ikl     )  ,ODA_RT(  ikl),FSOs_t(ikl)         &
     &                                            ,FIRn_t(ikl,mzpp)    &
     &                            ,Alb_SV_gpt(ikl),SWAsRT(ikl)         &
     &                            ,EmisSV_gpt(ikl),LWDsRT(ikl)         &
     &                             ,LWUpwd
  403       format(                                                    &
     &           ' +---------+---------+--------+---------+--------+', &
     &              '--------+--------+--------+---------+---------+', &
     &       '-------+---------+---------+',                           &
     &                   /,' | AIR-SOL |',f8.2,' |',f7.2,' |', 8x ,' |'&
     &        , 3( 7x ,' |'),   f7.2,' |',e8.2,' |',f8.1,' |'          &
     &        ,                            6x ,' |', 8x ,' |',f8.1,' |'&
     &        ,          /,' |     SOL |', 8x ,' |', 7x ,' |', 8x ,' |'&
     &        , 3( 7x ,' |'),    7x ,' |',f8.2,' |',f8.1,' |'          &
     &        ,                           f6.3,' |',f8.1,' |',f8.1,' |')

      END IF !    (iklOUT.GT.0)


! OUTPUT for Debugging
! --------------------

! #DB     jkjllw=0
! #DB DO JL=1,KLON
! #DB     lijio =0
! #DB     ikl   = min(JL,kcolp)
! #DB DO JK=1,KLEV
!!      write(6,*) k2ii(JL),k2jj(JL),jl,jk,' FIRn_t: ', PFLT5(jl,jk)
!!      write(6,*) k2ii(JL),k2jj(JL),jl,jk,' FSOn_t: ', PFLS5(jl,jk)
! #DB   IF ( PFLT5(jl,jk).GT. 500..OR. PFLS5(jl,jk).GT. 500. .OR.              &
! #DB&       PFLT5(jl,jk).LT.-500..OR. PFLS5(jl,jk).LT.-500. .OR.              &
! #DB&       (k2ii(jl).EQ.kio.AND.k2jj(jl).EQ.kjo)) lijio=1
! #DB END DO
! #DB IF   (lijio.EQ.1)                                             THEN
! #DB   DO JK=1,KLEV
! #DB     IF (mod(jkjllw,20).EQ.0)                                             &
! #DB&        write(6,600)
! #DB     600 format('IN   PHYrad2CEP: Radiative Fluxes ',/                    &
! #DB&              ,'    i    j   JL   JK',9x,'Ta',9x,'Qv',9x,'Qi',9x,'Qw'    &
! #DB&              ,9x,'O3',8x,'CLD',8x,'COD',8x,'AOD',8x,'SOn',8x,'IRn')
! #DB         jkjllw=jkjllw+1
! #DB         write(6,601) k2ii(JL),k2jj(JL),JL,JK                             &
! #DB&                ,Ta__DY(ikl,JK),PQ5   (jk,jl)                            &
! #DB&                ,PQIWP5(JL,JK) ,PQLWP5(JL,JK)                            &
! #DB&                ,O3__RT(ikl,jk),CFraCM(ikl,JK)                           &
! #DB&                ,ODCzRT(ikl,jk),ODAzRT(ikl,jk)                           &
! #DB&                ,FSOn_t(ikl,jk),FIRn_t(ikl,jk)
! #DB     601  format(4i5,10e11.3)
! #DB   END DO
! #DB ENDIF
!
! #DB END DO




      return
      end subroutine Atm_RT_RUN
