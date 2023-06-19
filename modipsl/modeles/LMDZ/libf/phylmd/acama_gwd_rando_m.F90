!
! $Id: acama_gwd_rando_m.F90 3531 2019-06-06 15:08:45Z fairhead $
!
module ACAMA_GWD_rando_m

  implicit none

contains

  SUBROUTINE ACAMA_GWD_rando(DTIME, pp, plat, tt, uu, vv, rot, &
       zustr, zvstr, d_u, d_v,east_gwstress,west_gwstress)

    ! Parametrization of the momentum flux deposition due to a discrete
    ! number of gravity waves. 
    ! Author: F. Lott, A. de la Camara
    ! July, 24th, 2014
    ! Gaussian distribution of the source, source is vorticity squared
    ! Reference: de la Camara and Lott (GRL, 2015, vol 42, 2071-2078 )
    ! Lott et al (JAS, 2010, vol 67, page 157-170)
    ! Lott et al (JAS, 2012, vol 69, page 2134-2151)

!  ONLINE:
    use dimphy, only: klon, klev
    use assert_m, only: assert
    USE ioipsl_getin_p_mod, ONLY : getin_p
    USE vertical_layers_mod, ONLY : presnivs

    include "YOMCST.h"
    include "clesphys.h"
!  OFFLINE:
!   include "dimensions.h"
!   include "dimphy.h"
!END DIFFERENCE
    include "YOEGWD.h"

    ! 0. DECLARATIONS:

    ! 0.1 INPUTS
    REAL, intent(in)::DTIME ! Time step of the Physics
    REAL, intent(in):: PP(:, :) ! (KLON, KLEV) Pressure at full levels
    REAL, intent(in):: ROT(:,:) ! Relative vorticity              
    REAL, intent(in):: TT(:, :) ! (KLON, KLEV) Temp at full levels 
    REAL, intent(in):: UU(:, :) ! (KLON, KLEV) Zonal wind at full levels
    REAL, intent(in):: VV(:, :) ! (KLON, KLEV) Merid wind at full levels
    REAL, intent(in):: PLAT(:) ! (KLON) LATITUDE                   

    ! 0.2 OUTPUTS
    REAL, intent(out):: zustr(:), zvstr(:) ! (KLON) Surface Stresses

    REAL, intent(inout):: d_u(:, :), d_v(:, :) 
    REAL, intent(inout):: east_gwstress(:, :) !  Profile of eastward stress
    REAL, intent(inout):: west_gwstress(:, :) !  Profile of westward stress 
    ! (KLON, KLEV) tendencies on winds

    ! O.3 INTERNAL ARRAYS
    REAL BVLOW(klon)  !  LOW LEVEL BV FREQUENCY
    REAL ROTBA(KLON),CORIO(KLON)  !  BAROTROPIC REL. VORTICITY AND PLANETARY
    REAL UZ(KLON, KLEV + 1)

    INTEGER II, JJ, LL

    ! 0.3.0 TIME SCALE OF THE LIFE CYCLE OF THE WAVES PARAMETERIZED

    REAL DELTAT

    ! 0.3.1 GRAVITY-WAVES SPECIFICATIONS

    INTEGER, PARAMETER:: NK = 2, NP = 2, NO = 2, NW = NK * NP * NO
    INTEGER JK, JP, JO, JW
    INTEGER, PARAMETER:: NA = 5  !number of realizations to get the phase speed
    REAL KMIN, KMAX ! Min and Max horizontal wavenumbers
    REAL CMIN, CMAX ! Min and Max absolute ph. vel.
    REAL CPHA ! absolute PHASE VELOCITY frequency
    REAL ZK(NW, KLON) ! Horizontal wavenumber amplitude
    REAL ZP(NW, KLON) ! Horizontal wavenumber angle 
    REAL ZO(NW, KLON) ! Absolute frequency !

    ! Waves Intr. freq. at the 1/2 lev surrounding the full level
    REAL ZOM(NW, KLON), ZOP(NW, KLON)

    ! Wave EP-fluxes at the 2 semi levels surrounding the full level
    REAL WWM(NW, KLON), WWP(NW, KLON)

    REAL RUW0(NW, KLON) ! Fluxes at launching level

    REAL RUWP(NW, KLON), RVWP(NW, KLON)
    ! Fluxes X and Y for each waves at 1/2 Levels

    INTEGER LAUNCH, LTROP ! Launching altitude and tropo altitude

    REAL XLAUNCH ! Controle the launching altitude
    REAL XTROP ! SORT of Tropopause altitude 
    REAL RUW(KLON, KLEV + 1) ! Flux x at semi levels
    REAL RVW(KLON, KLEV + 1) ! Flux y at semi levels

    REAL PRMAX ! Maximum value of PREC, and for which our linear formula
    ! for GWs parameterisation apply

    ! 0.3.2 PARAMETERS OF WAVES DISSIPATIONS

    REAL RDISS, ZOISEC ! COEFF DE DISSIPATION, SECURITY FOR INTRINSIC FREQ
    REAL CORSEC ! SECURITY FOR INTRINSIC CORIOLIS
    REAL RUWFRT,SATFRT

    ! 0.3.3 BACKGROUND FLOW AT 1/2 LEVELS AND VERTICAL COORDINATE

    REAL H0 ! Characteristic Height of the atmosphere
    REAL DZ ! Characteristic depth of the source!
    REAL PR, TR ! Reference Pressure and Temperature

    REAL ZH(KLON, KLEV + 1) ! Log-pressure altitude

    REAL UH(KLON, KLEV + 1), VH(KLON, KLEV + 1) ! Winds at 1/2 levels
    REAL PH(KLON, KLEV + 1) ! Pressure at 1/2 levels
    REAL PSEC ! Security to avoid division by 0 pressure
    REAL PHM1(KLON, KLEV + 1) ! 1/Press at 1/2 levels
    REAL BV(KLON, KLEV + 1) ! Brunt Vaisala freq. (BVF) at 1/2 levels
    REAL BVSEC ! Security to avoid negative BVF

    REAL, DIMENSION(klev+1) ::HREF
    LOGICAL, SAVE :: gwd_reproductibilite_mpiomp=.true.
    LOGICAL, SAVE :: firstcall = .TRUE.
  !$OMP THREADPRIVATE(firstcall,gwd_reproductibilite_mpiomp)

    CHARACTER (LEN=20) :: modname='acama_gwd_rando_m'
    CHARACTER (LEN=80) :: abort_message



  IF (firstcall) THEN
    ! Cle introduite pour resoudre un probleme de non reproductibilite
    ! Le but est de pouvoir tester de revenir a la version precedenete
    ! A eliminer rapidement
    CALL getin_p('gwd_reproductibilite_mpiomp',gwd_reproductibilite_mpiomp)
    IF (NW+4*(NA-1)+NA>=KLEV) THEN
       abort_message = 'NW+3*NA>=KLEV Probleme pour generation des ondes'
       CALL abort_physic (modname,abort_message,1)
    ENDIF
    firstcall=.false.
!    CALL iophys_ini
  ENDIF

    !-----------------------------------------------------------------

    ! 1. INITIALISATIONS

    ! 1.1 Basic parameter

    ! Are provided from elsewhere (latent heat of vaporization, dry
    ! gaz constant for air, gravity constant, heat capacity of dry air
    ! at constant pressure, earth rotation rate, pi).

    ! 1.2 Tuning parameters of V14

! Values for linear in rot (recommended):
!   RUWFRT=0.005 ! As RUWMAX but for frontal waves 
!   SATFRT=1.00  ! As SAT    but for frontal waves
! Values when rot^2 is used                          
!    RUWFRT=0.02  ! As RUWMAX but for frontal waves 
!    SATFRT=1.00  ! As SAT    but for frontal waves
!    CMAX = 30.   ! Characteristic phase speed
! Values when rot^2*EXP(-pi*sqrt(J)) is used                          
!   RUWFRT=2.5   ! As RUWMAX but for frontal waves ~ N0*F0/4*DZ
!   SATFRT=0.60   ! As SAT    but for frontal waves
    RUWFRT=gwd_front_ruwmax  
    SATFRT=gwd_front_sat
    CMAX = 50.    ! Characteristic phase speed
! Phase speed test
!   RUWFRT=0.01
!   CMAX = 50.   ! Characteristic phase speed (TEST)
! Values when rot^2 and exp(-m^2*dz^2) are used      
!   RUWFRT=0.03  ! As RUWMAX but for frontal waves 
!   SATFRT=1.00  ! As SAT    but for frontal waves
! CRUCIAL PARAMETERS FOR THE WIND FILTERING
    XLAUNCH=0.95 ! Parameter that control launching altitude
    RDISS = 0.5  ! Diffusion parameter 

    ! maximum of rain for which our theory applies (in kg/m^2/s)

    DZ = 1000. ! Characteristic depth of the source
    XTROP=0.2 ! Parameter that control tropopause altitude
    DELTAT=24.*3600. ! Time scale of the waves (first introduced in 9b)
!   DELTAT=DTIME     ! No AR-1 Accumulation, OR OFFLINE             

    KMIN = 2.E-5
    ! minimum horizontal wavenumber (inverse of the subgrid scale resolution)

    KMAX = 1.E-3  ! Max horizontal wavenumber
    CMIN = 1.     ! Min phase velocity

    TR = 240. ! Reference Temperature
    PR = 101300. ! Reference pressure
    H0 = RD * TR / RG ! Characteristic vertical scale height

    BVSEC = 5.E-3 ! Security to avoid negative BVF 
    PSEC = 1.E-6 ! Security to avoid division by 0 pressure
    ZOISEC = 1.E-6 ! Security FOR 0 INTRINSIC FREQ
    CORSEC = ROMEGA*2.*SIN(2.*RPI/180.)! Security for CORIO

!  ONLINE
    call assert(klon == (/size(pp, 1), size(tt, 1), size(uu, 1), &
         size(vv, 1), size(rot,1), size(zustr), size(zvstr), size(d_u, 1), &
         size(d_v, 1), &
        size(east_gwstress,1), size(west_gwstress,1) /), &
        "ACAMA_GWD_RANDO klon")
    call assert(klev == (/size(pp, 2), size(tt, 2), size(uu, 2), &
         size(vv, 2), size(d_u, 2), size(d_v, 2), &
         size(east_gwstress,2), size(west_gwstress,2) /), &
         "ACAMA_GWD_RANDO klev")
!  END ONLINE

    IF(DELTAT < DTIME)THEN
!       PRINT *, 'flott_gwd_rando: deltat < dtime!'
!       STOP 1
       abort_message=' deltat < dtime! '
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    IF (KLEV < NW) THEN
!       PRINT *, 'flott_gwd_rando: you will have problem with random numbers'
!       STOP 1
       abort_message=' you will have problem with random numbers'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    ! 2. EVALUATION OF THE BACKGROUND FLOW AT SEMI-LEVELS

    ! Pressure and Inv of pressure
    DO LL = 2, KLEV
       PH(:, LL) = EXP((LOG(PP(:, LL)) + LOG(PP(:, LL - 1))) / 2.)
       PHM1(:, LL) = 1. / PH(:, LL) 
    end DO

    PH(:, KLEV + 1) = 0. 
    PHM1(:, KLEV + 1) = 1. / PSEC
    PH(:, 1) = 2. * PP(:, 1) - PH(:, 2)

    ! Launching altitude

    IF (gwd_reproductibilite_mpiomp) THEN
       ! Reprend la formule qui calcule PH en fonction de PP=play
       DO LL = 2, KLEV
          HREF(LL) = EXP((LOG(presnivs(LL)) + LOG(presnivs(LL - 1))) / 2.)
       end DO
       HREF(KLEV + 1) = 0.
       HREF(1) = 2. * presnivs(1) - HREF(2)
    ELSE
       HREF(1:KLEV)=PH(KLON/2,1:KLEV)
    ENDIF

    LAUNCH=0
    LTROP =0
    DO LL = 1, KLEV
       IF (HREF(LL) / HREF(1) > XLAUNCH) LAUNCH = LL
    ENDDO
    DO LL = 1, KLEV
       IF (HREF(LL) / HREF(1) > XTROP) LTROP = LL
    ENDDO
    !LAUNCH=22 ; LTROP=33
!   print*,'LAUNCH=',LAUNCH,'LTROP=',LTROP


!   PRINT *,'LAUNCH IN ACAMARA:',LAUNCH

    ! Log pressure vert. coordinate
    DO LL = 1, KLEV + 1 
       ZH(:, LL) = H0 * LOG(PR / (PH(:, LL) + PSEC))
    end DO

    ! BV frequency
    DO LL = 2, KLEV
       ! BVSEC: BV Frequency (UH USED IS AS A TEMPORARY ARRAY DOWN TO WINDS)
       UH(:, LL) = 0.5 * (TT(:, LL) + TT(:, LL - 1)) &
            * RD**2 / RCPD / H0**2 + (TT(:, LL) &
            - TT(:, LL - 1)) / (ZH(:, LL) - ZH(:, LL - 1)) * RD / H0
    end DO
    BVLOW = 0.5 * (TT(:, LTROP )+ TT(:, LAUNCH)) &
         * RD**2 / RCPD / H0**2 + (TT(:, LTROP ) &
         - TT(:, LAUNCH))/(ZH(:, LTROP )- ZH(:, LAUNCH)) * RD / H0

    UH(:, 1) = UH(:, 2)
    UH(:, KLEV + 1) = UH(:, KLEV)
    BV(:, 1) = UH(:, 2)
    BV(:, KLEV + 1) = UH(:, KLEV)
    ! SMOOTHING THE BV HELPS
    DO LL = 2, KLEV
       BV(:, LL)=(UH(:, LL+1)+2.*UH(:, LL)+UH(:, LL-1))/4.
    end DO

    BV=MAX(SQRT(MAX(BV, 0.)), BVSEC)
    BVLOW=MAX(SQRT(MAX(BVLOW, 0.)), BVSEC)

    ! WINDS
    DO LL = 2, KLEV
       UH(:, LL) = 0.5 * (UU(:, LL) + UU(:, LL - 1)) ! Zonal wind
       VH(:, LL) = 0.5 * (VV(:, LL) + VV(:, LL - 1)) ! Meridional wind
       UZ(:, LL) = ABS((SQRT(UU(:, LL)**2+VV(:, LL)**2) & 
          - SQRT(UU(:,LL-1)**2+VV(:, LL-1)**2)) &
          /(ZH(:, LL)-ZH(:, LL-1)) )
    end DO
    UH(:, 1) = 0.
    VH(:, 1) = 0.
    UH(:, KLEV + 1) = UU(:, KLEV)
    VH(:, KLEV + 1) = VV(:, KLEV)

    UZ(:, 1) = UZ(:, 2)
    UZ(:, KLEV + 1) = UZ(:, KLEV)
    UZ(:, :) = MAX(UZ(:,:), PSEC)

   ! BAROTROPIC VORTICITY AND INTEGRATED CORIOLIS PARAMETER
    
    CORIO(:) = MAX(ROMEGA*2.*ABS(SIN(PLAT(:)*RPI/180.)),CORSEC)
    ROTBA(:)=0.
    DO LL = 1,KLEV-1
        !ROTBA(:) = ROTBA(:) + (ROT(:,LL)+ROT(:,LL+1))/2./RG*(PP(:,LL)-PP(:,LL+1))
        ! Introducing the complete formula (exp of Richardson number):
        ROTBA(:) = ROTBA(:) + &
                !((ROT(:,LL)+ROT(:,LL+1))/2.)**2 &
                (CORIO(:)*TANH(ABS(ROT(:,LL)+ROT(:,LL+1))/2./CORIO(:)))**2 &
                /RG*(PP(:,LL)-PP(:,LL+1)) &
                * EXP(-RPI*BV(:,LL+1)/UZ(:,LL+1)) &
!               * DZ*BV(:,LL+1)/4./ABS(CORIO(:)) 
                * DZ*BV(:,LL+1)/4./1.E-4           !  Changes after 1991
!ARRET
    ENDDO
    !   PRINT *,'MAX ROTBA:',MAXVAL(ROTBA)
    !   ROTBA(:)=(1.*ROTBA(:)  & ! Testing zone
    !           +0.15*CORIO(:)**2                &
    !           /(COS(PLAT(:)*RPI/180.)+0.02)  &
    !           )*DZ*0.01/0.0001/4. ! & ! Testing zone
    !   MODIF GWD4 AFTER 1985
    !            *(1.25+SIN(PLAT(:)*RPI/180.))/(1.05+SIN(PLAT(:)*RPI/180.))/1.25
    !          *1./(COS(PLAT(:)*RPI/180.)+0.02)
    !    CORIO(:) = MAX(ROMEGA*2.*ABS(SIN(PLAT(:)*RPI/180.)),ZOISEC)/RG*PP(:,1)

    ! 3 WAVES CHARACTERISTICS CHOSEN RANDOMLY AT THE LAUNCH ALTITUDE

    ! The mod functions of weird arguments are used to produce the
    ! waves characteristics in an almost stochastic way

    JW = 0
    DO JW = 1, NW
             ! Angle
             DO II = 1, KLON
                ! Angle (0 or PI so far)
                ! ZP(JW, II) = (SIGN(1., 0.5 - MOD(TT(II, JW) * 10., 1.)) + 1.) &
                !      * RPI / 2.
                ! Angle between 0 and pi
                  ZP(JW, II) = MOD(TT(II, JW) * 10., 1.) * RPI
! TEST WITH POSITIVE WAVES ONLY (Part I/II)
!               ZP(JW, II) = 0.
                ! Horizontal wavenumber amplitude
                ZK(JW, II) = KMIN + (KMAX - KMIN) * MOD(TT(II, JW) * 100., 1.)
                ! Horizontal phase speed
                CPHA = 0.
                DO JJ = 1, NA
                    CPHA = CPHA + &
         CMAX*2.*(MOD(TT(II, JW+4*(JJ-1)+JJ)**2, 1.)-0.5)*SQRT(3.)/SQRT(NA*1.)
                END DO
                IF (CPHA.LT.0.)  THEN
                   CPHA = -1.*CPHA
                   ZP(JW,II) = ZP(JW,II) + RPI
! TEST WITH POSITIVE WAVES ONLY (Part II/II)
!               ZP(JW, II) = 0.
                ENDIF
                CPHA = CPHA + CMIN !we dont allow |c|<1m/s
                ! Absolute frequency is imposed
                ZO(JW, II) = CPHA * ZK(JW, II) 
                ! Intrinsic frequency is imposed
                ZO(JW, II) = ZO(JW, II) &
                     + ZK(JW, II) * COS(ZP(JW, II)) * UH(II, LAUNCH) &
                     + ZK(JW, II) * SIN(ZP(JW, II)) * VH(II, LAUNCH)
                ! Momentum flux at launch lev 
                ! LAUNCHED RANDOM WAVES WITH LOG-NORMAL AMPLITUDE
                !  RIGHT IN THE SH (GWD4 after 1990)
                  RUW0(JW, II) = 0.
                 DO JJ = 1, NA
                    RUW0(JW, II) = RUW0(JW,II) + &
         2.*(MOD(TT(II, JW+4*(JJ-1)+JJ)**2, 1.)-0.5)*SQRT(3.)/SQRT(NA*1.)
                END DO
                RUW0(JW, II) = RUWFRT &
                          * EXP(RUW0(JW,II))/1250. &  ! 2 mpa at south pole
       *((1.05+SIN(PLAT(II)*RPI/180.))/(1.01+SIN(PLAT(II)*RPI/180.))-2.05/2.01)
                ! RUW0(JW, II) = RUWFRT
             ENDDO
    end DO

    ! 4. COMPUTE THE FLUXES

    ! 4.0 

    ! 4.1 Vertical velocity at launching altitude to ensure 
    ! the correct value to the imposed fluxes.

    DO JW = 1, NW

       ! Evaluate intrinsic frequency at launching altitude:
       ZOP(JW, :) = ZO(JW, :) &
            - ZK(JW, :) * COS(ZP(JW, :)) * UH(:, LAUNCH) &
            - ZK(JW, :) * SIN(ZP(JW, :)) * VH(:, LAUNCH) 

       ! VERSION WITH FRONTAL SOURCES

       ! Momentum flux at launch level imposed by vorticity sources

       ! tanh limitation for values above CORIO (inertial instability).
       ! WWP(JW, :) = RUW0(JW, :) &
       WWP(JW, :) = RUWFRT      &
       !     * (CORIO(:)*TANH(ROTBA(:)/CORIO(:)))**2 &
       !    * ABS((CORIO(:)*TANH(ROTBA(:)/CORIO(:)))*CORIO(:)) &
       !  CONSTANT FLUX
       !    * (CORIO(:)*CORIO(:)) &
       ! MODERATION BY THE DEPTH OF THE SOURCE (DZ HERE)
       !      *EXP(-BVLOW(:)**2/MAX(ABS(ZOP(JW, :)),ZOISEC)**2 &
       !      *ZK(JW, :)**2*DZ**2) &
       ! COMPLETE FORMULA:
            !* CORIO(:)**2*TANH(ROTBA(:)/CORIO(:)**2) &
            * ROTBA(:) &
       !  RESTORE DIMENSION OF A FLUX
       !     *RD*TR/PR
       !     *1. + RUW0(JW, :)
             *1.

       ! Factor related to the characteristics of the waves: NONE

       ! Moderation by the depth of the source (dz here): NONE

       ! Put the stress in the right direction:

        RUWP(JW, :) = SIGN(1., ZOP(JW, :))*COS(ZP(JW, :)) * WWP(JW, :)
        RVWP(JW, :) = SIGN(1., ZOP(JW, :))*SIN(ZP(JW, :)) * WWP(JW, :)

    end DO

    ! 4.2 Uniform values below the launching altitude

    DO LL = 1, LAUNCH
       RUW(:, LL) = 0
       RVW(:, LL) = 0
       DO JW = 1, NW
          RUW(:, LL) = RUW(:, LL) + RUWP(JW, :)
          RVW(:, LL) = RVW(:, LL) + RVWP(JW, :)
       end DO
    end DO

    ! 4.3 Loop over altitudes, with passage from one level to the next
    ! done by i) conserving the EP flux, ii) dissipating a little,
    ! iii) testing critical levels, and vi) testing the breaking.

    DO LL = LAUNCH, KLEV - 1
       ! Warning: all the physics is here (passage from one level
       ! to the next)
       DO JW = 1, NW
          ZOM(JW, :) = ZOP(JW, :)
          WWM(JW, :) = WWP(JW, :)
          ! Intrinsic Frequency
          ZOP(JW, :) = ZO(JW, :) - ZK(JW, :) * COS(ZP(JW, :)) * UH(:, LL + 1) &
               - ZK(JW, :) * SIN(ZP(JW, :)) * VH(:, LL + 1) 

          ! No breaking (Eq.6)
          ! Dissipation (Eq. 8)
          WWP(JW, :) = WWM(JW, :) * EXP(- 4. * RDISS * PR / (PH(:, LL + 1) &
               + PH(:, LL)) * ((BV(:, LL + 1) + BV(:, LL)) / 2.)**3 &
               / MAX(ABS(ZOP(JW, :) + ZOM(JW, :)) / 2., ZOISEC)**4 &
               * ZK(JW, :)**3 * (ZH(:, LL + 1) - ZH(:, LL)))

          ! Critical levels (forced to zero if intrinsic frequency changes sign)
          ! Saturation (Eq. 12)
          WWP(JW, :) = min(WWP(JW, :), MAX(0., &
               SIGN(1., ZOP(JW, :) * ZOM(JW, :))) * ABS(ZOP(JW, :))**3 &
          !    / BV(:, LL + 1) * EXP(- ZH(:, LL + 1) / H0) * SATFRT**2 * KMIN**2 &
               / BV(:, LL + 1) * EXP(- ZH(:, LL + 1) / H0) * KMIN**2 &
!              *(SATFRT*(2.5+1.5*TANH((ZH(:,LL+1)/H0-8.)/2.)))**2 &
               *SATFRT**2       &
               / ZK(JW, :)**4)
       end DO

       ! Evaluate EP-flux from Eq. 7 and give the right orientation to
       ! the stress

       DO JW = 1, NW
          RUWP(JW, :) = SIGN(1., ZOP(JW, :))*COS(ZP(JW, :)) * WWP(JW, :)
          RVWP(JW, :) = SIGN(1., ZOP(JW, :))*SIN(ZP(JW, :)) * WWP(JW, :)
       end DO

       RUW(:, LL + 1) = 0.
       RVW(:, LL + 1) = 0.

       DO JW = 1, NW
          RUW(:, LL + 1) = RUW(:, LL + 1) + RUWP(JW, :) 
          RVW(:, LL + 1) = RVW(:, LL + 1) + RVWP(JW, :) 
          EAST_GWSTRESS(:, LL)=EAST_GWSTRESS(:, LL)+MAX(0.,RUWP(JW,:))/FLOAT(NW)
          WEST_GWSTRESS(:, LL)=WEST_GWSTRESS(:, LL)+MIN(0.,RUWP(JW,:))/FLOAT(NW)
       end DO
    end DO

    ! 5 CALCUL DES TENDANCES:

    ! 5.1 Rectification des flux au sommet et dans les basses couches

    RUW(:, KLEV + 1) = 0.
    RVW(:, KLEV + 1) = 0.
    RUW(:, 1) = RUW(:, LAUNCH)
    RVW(:, 1) = RVW(:, LAUNCH)
    DO LL = 1, LAUNCH
       RUW(:, LL) = RUW(:, LAUNCH+1)
       RVW(:, LL) = RVW(:, LAUNCH+1)
       EAST_GWSTRESS(:, LL)=EAST_GWSTRESS(:, LAUNCH)
       WEST_GWSTRESS(:, LL)=WEST_GWSTRESS(:, LAUNCH)
    end DO

    ! AR-1 RECURSIVE FORMULA (13) IN VERSION 4
    DO LL = 1, KLEV
       D_U(:, LL) = (1.-DTIME/DELTAT) * D_U(:, LL) + DTIME/DELTAT/REAL(NW) * &
            RG * (RUW(:, LL + 1) - RUW(:, LL)) &
            / (PH(:, LL + 1) - PH(:, LL)) * DTIME
!  NO AR1 FOR MERIDIONAL TENDENCIES
!      D_V(:, LL) = (1.-DTIME/DELTAT) * D_V(:, LL) + DTIME/DELTAT/REAL(NW) * &
       D_V(:, LL) =                                            1./REAL(NW) * &
            RG * (RVW(:, LL + 1) - RVW(:, LL)) &
            / (PH(:, LL + 1) - PH(:, LL)) * DTIME
    ENDDO

    ! Cosmetic: evaluation of the cumulated stress
    ZUSTR = 0.
    ZVSTR = 0.
    DO LL = 1, KLEV
       ZUSTR = ZUSTR + D_U(:, LL) / RG * (PH(:, LL + 1) - PH(:, LL))/DTIME
!      ZVSTR = ZVSTR + D_V(:, LL) / RG * (PH(:, LL + 1) - PH(:, LL))/DTIME
    ENDDO
! COSMETICS TO VISUALIZE ROTBA
    ZVSTR = ROTBA

  END SUBROUTINE ACAMA_GWD_RANDO

end module ACAMA_GWD_rando_m
