!
! $Id: flott_gwd_rando_m.F90 3531 2019-06-06 15:08:45Z fairhead $
!
module FLOTT_GWD_rando_m

  implicit none

contains

  SUBROUTINE FLOTT_GWD_rando(DTIME, pp, tt, uu, vv, prec, zustr, zvstr, d_u, &
       d_v,east_gwstress,west_gwstress)

    ! Parametrization of the momentum flux deposition due to a discrete
    ! number of gravity waves. 
    ! Author: F. Lott
    ! July, 12th, 2012
    ! Gaussian distribution of the source, source is precipitation
    ! Reference: Lott (JGR, vol 118, page 8897, 2013)

    !ONLINE:
      use dimphy, only: klon, klev
      use assert_m, only: assert
      USE ioipsl_getin_p_mod, ONLY : getin_p
      USE vertical_layers_mod, ONLY : presnivs
      CHARACTER (LEN=20) :: modname='flott_gwd_rando'
      CHARACTER (LEN=80) :: abort_message

      include "YOMCST.h"
      include "clesphys.h"
    ! OFFLINE:
    ! include "dimensions.h"
    ! include "dimphy.h"
    ! END OF DIFFERENCE ONLINE-OFFLINE
    include "YOEGWD.h"

    ! 0. DECLARATIONS:

    ! 0.1 INPUTS
    REAL, intent(in)::DTIME ! Time step of the Physics
    REAL, intent(in):: pp(:, :) ! (KLON, KLEV) Pressure at full levels
    REAL, intent(in):: prec(:) ! (klon) Precipitation (kg/m^2/s) 
    REAL, intent(in):: TT(:, :) ! (KLON, KLEV) Temp at full levels 
    REAL, intent(in):: UU(:, :) ! (KLON, KLEV) Zonal wind at full levels
    REAL, intent(in):: VV(:, :) ! (KLON, KLEV) Merid wind at full levels

    ! 0.2 OUTPUTS
    REAL, intent(out):: zustr(:), zvstr(:) ! (KLON) Surface Stresses

    REAL, intent(inout):: d_u(:, :), d_v(:, :) 
    REAL, intent(inout):: east_gwstress(:, :) !  Profile of eastward stress
    REAL, intent(inout):: west_gwstress(:, :) !  Profile of westward stress 

    ! (KLON, KLEV) tendencies on winds

    ! O.3 INTERNAL ARRAYS
    REAL BVLOW(klon)
    REAL DZ   !  Characteristic depth of the Source

    INTEGER II, JJ, LL

    ! 0.3.0 TIME SCALE OF THE LIFE CYCLE OF THE WAVES PARAMETERIZED

    REAL DELTAT

    ! 0.3.1 GRAVITY-WAVES SPECIFICATIONS

    INTEGER, PARAMETER:: NK = 2, NP = 2, NO = 2, NW = NK * NP * NO
    INTEGER JK, JP, JO, JW
    INTEGER, PARAMETER:: NA = 5  !number of realizations to get the phase speed
    REAL KMIN, KMAX ! Min and Max horizontal wavenumbers
    REAL CMAX ! standard deviation of the phase speed distribution
    REAL RUWMAX,SAT  ! ONLINE SPECIFIED IN run.def
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

    ! 0.3.3 BACKGROUND FLOW AT 1/2 LEVELS AND VERTICAL COORDINATE

    REAL H0 ! Characteristic Height of the atmosphere
    REAL PR, TR ! Reference Pressure and Temperature

    REAL ZH(KLON, KLEV + 1) ! Log-pressure altitude

    REAL UH(KLON, KLEV + 1), VH(KLON, KLEV + 1) ! Winds at 1/2 levels
    REAL PH(KLON, KLEV + 1) ! Pressure at 1/2 levels
    REAL PSEC ! Security to avoid division by 0 pressure
    REAL BV(KLON, KLEV + 1) ! Brunt Vaisala freq. (BVF) at 1/2 levels
    REAL BVSEC ! Security to avoid negative BVF
    REAL RAN_NUM_1,RAN_NUM_2,RAN_NUM_3

    REAL, DIMENSION(klev+1) ::HREF

    LOGICAL, SAVE :: gwd_reproductibilite_mpiomp=.true.
    LOGICAL, SAVE :: firstcall = .TRUE.
  !$OMP THREADPRIVATE(firstcall,gwd_reproductibilite_mpiomp)


  IF (firstcall) THEN
    ! Cle introduite pour resoudre un probleme de non reproductibilite
    ! Le but est de pouvoir tester de revenir a la version precedenete
    ! A eliminer rapidement
    CALL getin_p('gwd_reproductibilite_mpiomp',gwd_reproductibilite_mpiomp)
    IF (NW+3*NA>=KLEV) THEN
       abort_message = 'NW+3*NA>=KLEV Probleme pour generation des ondes'
       CALL abort_physic (modname,abort_message,1)
    ENDIF
    firstcall=.false.
  ENDIF


    !-----------------------------------------------------------------

    ! 1. INITIALISATIONS

    ! 1.1 Basic parameter

    ! Are provided from elsewhere (latent heat of vaporization, dry
    ! gaz constant for air, gravity constant, heat capacity of dry air
    ! at constant pressure, earth rotation rate, pi).

    ! 1.2 Tuning parameters of V14

    
    RDISS = 0.5 ! Diffusion parameter
    ! ONLINE 
      RUWMAX=GWD_RANDO_RUWMAX
      SAT=gwd_rando_sat
    !END ONLINE
    ! OFFLINE
    ! RUWMAX= 1.75    ! Launched flux
    ! SAT=0.25     ! Saturation parameter
    ! END OFFLINE

    PRMAX = 20. / 24. /3600.
    ! maximum of rain for which our theory applies (in kg/m^2/s)

 ! Characteristic depth of the source
    DZ = 1000.
    XLAUNCH=0.5 ! Parameter that control launching altitude
    XTROP=0.2 ! Parameter that control tropopause altitude
    DELTAT=24.*3600. ! Time scale of the waves (first introduced in 9b)
    !  OFFLINE
    !  DELTAT=DTIME
    !  END OFFLINE

    KMIN = 2.E-5
    ! minimum horizontal wavenumber (inverse of the subgrid scale resolution)

    KMAX = 1.E-3 ! Max horizontal wavenumber
    CMAX = 30. ! Max phase speed velocity

    TR = 240. ! Reference Temperature
    PR = 101300. ! Reference pressure
    H0 = RD * TR / RG ! Characteristic vertical scale height

    BVSEC = 5.E-3 ! Security to avoid negative BVF 
    PSEC = 1.E-6 ! Security to avoid division by 0 pressure
    ZOISEC = 1.E-6 ! Security FOR 0 INTRINSIC FREQ

IF (1==0) THEN
    !ONLINE
        call assert(klon == (/size(pp, 1), size(tt, 1), size(uu, 1), &
         size(vv, 1), size(zustr), size(zvstr), size(d_u, 1), &
         size(d_v, 1), &
         size(east_gwstress, 1), size(west_gwstress, 1) /), &
         "FLOTT_GWD_RANDO klon")
     call assert(klev == (/size(pp, 2), size(tt, 2), size(uu, 2), &
          size(vv, 2), size(d_u, 2), size(d_v, 2), &
          size(east_gwstress,2), size(west_gwstress,2) /), &
          "FLOTT_GWD_RANDO klev")
    !END ONLINE
ENDIF

    IF(DELTAT < DTIME)THEN
       abort_message='flott_gwd_rando: deltat < dtime!'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    IF (KLEV < NW) THEN
       abort_message='flott_gwd_rando: you will have problem with random numbers'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    ! 2. EVALUATION OF THE BACKGROUND FLOW AT SEMI-LEVELS

    ! Pressure and Inv of pressure
    DO LL = 2, KLEV
       PH(:, LL) = EXP((LOG(PP(:, LL)) + LOG(PP(:, LL - 1))) / 2.)
    end DO
    PH(:, KLEV + 1) = 0. 
    PH(:, 1) = 2. * PP(:, 1) - PH(:, 2)

    ! Launching altitude

    !Pour revenir a la version non reproductible en changeant le nombre de process
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
    BVLOW(:) = 0.5 * (TT(:, LTROP )+ TT(:, LAUNCH)) &
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
    end DO
    UH(:, 1) = 0.
    VH(:, 1) = 0.
    UH(:, KLEV + 1) = UU(:, KLEV)
    VH(:, KLEV + 1) = VV(:, KLEV)

    ! 3 WAVES CHARACTERISTICS CHOSEN RANDOMLY AT THE LAUNCH ALTITUDE

    ! The mod functions of weird arguments are used to produce the
    ! waves characteristics in an almost stochastic way

    DO JW = 1, NW
             ! Angle
             DO II = 1, KLON
                ! Angle (0 or PI so far)
                RAN_NUM_1=MOD(TT(II, JW) * 10., 1.)
                RAN_NUM_2= MOD(TT(II, JW) * 100., 1.)
                ZP(JW, II) = (SIGN(1., 0.5 - RAN_NUM_1) + 1.) &
                     * RPI / 2.
                ! Horizontal wavenumber amplitude
                ZK(JW, II) = KMIN + (KMAX - KMIN) *RAN_NUM_2
                ! Horizontal phase speed
                CPHA = 0.
                DO JJ = 1, NA
                    RAN_NUM_3=MOD(TT(II, JW+3*JJ)**2, 1.)
                    CPHA = CPHA + &
                    CMAX*2.*(RAN_NUM_3 -0.5)*SQRT(3.)/SQRT(NA*1.)
                END DO
                IF (CPHA.LT.0.)  THEN
                   CPHA = -1.*CPHA
                   ZP(JW,II) = ZP(JW,II) + RPI
                ENDIF
                ! Absolute frequency is imposed
                ZO(JW, II) = CPHA * ZK(JW, II) 
                ! Intrinsic frequency is imposed
                ZO(JW, II) = ZO(JW, II) &
                     + ZK(JW, II) * COS(ZP(JW, II)) * UH(II, LAUNCH) &
                     + ZK(JW, II) * SIN(ZP(JW, II)) * VH(II, LAUNCH)
                ! Momentum flux at launch lev 
                RUW0(JW, II) = RUWMAX
             ENDDO
    ENDDO

    ! 4. COMPUTE THE FLUXES

    ! 4.1 Vertical velocity at launching altitude to ensure 
    ! the correct value to the imposed fluxes.

    DO JW = 1, NW

       ! Evaluate intrinsic frequency at launching altitude:
       ZOP(JW, :) = ZO(JW, :) &
            - ZK(JW, :) * COS(ZP(JW, :)) * UH(:, LAUNCH) &
            - ZK(JW, :) * SIN(ZP(JW, :)) * VH(:, LAUNCH) 

       ! VERSION WITH CONVECTIVE SOURCE

       ! Vertical velocity at launch level, value to ensure the
       ! imposed factor related to the convective forcing:
       ! precipitations.

       ! tanh limitation to values above prmax:
       WWP(JW, :) = RUW0(JW, :) &
            * (RD / RCPD / H0 * RLVTT * PRMAX * TANH(PREC(:) / PRMAX))**2

       ! Factor related to the characteristics of the waves:
       WWP(JW, :) = WWP(JW, :) * ZK(JW, :)**3 / KMIN / BVLOW(:)  &
            / MAX(ABS(ZOP(JW, :)), ZOISEC)**3 

       ! Moderation by the depth of the source (dz here):
       WWP(JW, :) = WWP(JW, :) &
            * EXP(- BVLOW(:)**2 / MAX(ABS(ZOP(JW, :)), ZOISEC)**2 * ZK(JW, :)**2 &
            * DZ**2)

       ! Put the stress in the right direction:
       RUWP(JW, :) = ZOP(JW, :) / MAX(ABS(ZOP(JW, :)), ZOISEC)**2 &
            * BV(:, LAUNCH) * COS(ZP(JW, :)) * WWP(JW, :)**2
       RVWP(JW, :) = ZOP(JW, :) / MAX(ABS(ZOP(JW, :)), ZOISEC)**2 &
            * BV(:, LAUNCH) * SIN(ZP(JW, :)) * WWP(JW, :)**2
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
               / BV(:, LL + 1) * EXP(- ZH(:, LL + 1) / H0) * KMIN**2  &
               * SAT**2 / ZK(JW, :)**4)
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
! OFFLINE ONLY
!   PRINT *,'SAT PROFILE:'
!   DO LL=1,KLEV
!   PRINT *,ZH(KLON/2,LL)/1000.,SAT*(2.+TANH(ZH(KLON/2,LL)/H0-8.))
!   ENDDO

    ! 5 CALCUL DES TENDANCES:

    ! 5.1 Rectification des flux au sommet et dans les basses couches

    RUW(:, KLEV + 1) = 0.
    RVW(:, KLEV + 1) = 0.
    RUW(:, 1) = RUW(:, LAUNCH)
    RVW(:, 1) = RVW(:, LAUNCH)
    DO LL = 1, LAUNCH
       RUW(:, LL) = RUW(:, LAUNCH+1)
       RVW(:, LL) = RVW(:, LAUNCH+1)
       EAST_GWSTRESS(:, LL)  = EAST_GWSTRESS(:, LAUNCH)
       WEST_GWSTRESS(:, LL)  = WEST_GWSTRESS(:, LAUNCH)
    end DO

    ! AR-1 RECURSIVE FORMULA (13) IN VERSION 4
    DO LL = 1, KLEV
       D_U(:, LL) = (1.-DTIME/DELTAT) * D_U(:, LL) + DTIME/DELTAT/REAL(NW) * &
            RG * (RUW(:, LL + 1) - RUW(:, LL)) &
            / (PH(:, LL + 1) - PH(:, LL)) * DTIME
       ! NO AR-1 FOR MERIDIONAL TENDENCIES
       D_V(:, LL) =                                            1./REAL(NW) * &
            RG * (RVW(:, LL + 1) - RVW(:, LL)) &
            / (PH(:, LL + 1) - PH(:, LL)) * DTIME
    ENDDO

    ! Cosmetic: evaluation of the cumulated stress
    ZUSTR = 0.
    ZVSTR = 0.
    DO LL = 1, KLEV
       ZUSTR = ZUSTR + D_U(:, LL) / RG * (PH(:, LL + 1) - PH(:, LL))/DTIME
       ZVSTR = ZVSTR + D_V(:, LL) / RG * (PH(:, LL + 1) - PH(:, LL))/DTIME
    ENDDO


  END SUBROUTINE FLOTT_GWD_RANDO

end module FLOTT_GWD_rando_m
