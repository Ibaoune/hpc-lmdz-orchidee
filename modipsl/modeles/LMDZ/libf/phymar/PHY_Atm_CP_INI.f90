      subroutine PHY_Atm_CP_INI(mzc,kcolc)

!------------------------------------------------------------------------------+
!                                                         Mon 17-Jun-2013  MAR |
!     subroutine PHY_Atm_CP_INI    interfaces                                  |
!                Bechtold et al. (2001) Convective Parameterization with MAR   |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Mon  8-Apr-2013      |
!           Last Modification by H. Gallee,               Mon 17-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+
!                                                                              |
!     INPUT                                                                    |
!     ^^^^^        it_EXP             : Experiment Iteration Counter           |
!                  dt__CP             : Mass Flux  Scheme:   Time Step         |
!                  dxHOST             : grid spacing  of  HOST MODEL       [m] |
!                  qs__CM(kcolp,mzpp) : air  snow  Particl. concentr.  [kg/kg] |
!                  qr__CM(kcolp,mzp ) : air  rain  drops    concentr.  [kg/kg] |
!                                                                              |
!     INPUT/OUTPUT                                                             |
!     ^^^^^^^^^^^^                                                             |
!                  qw__CM(kcolp,mzp ) : air  cloud droplets concentr.  [kg/kg] |
!                  qi__CM(kcolp,mzp ) : air  cloud crystals concentr.  [kg/kg] |
!                  snowCP(kcolp     ) : Snow (convective)                  [m] |
!                  snowCM(kcolp     ) : Snow (convective + stratiform)     [m] |
!                  rainCP(kcolp     ) : Rain (convective)                  [m] |
!                  rainCM(kcolp     ) : Rain (convective + stratiform)     [m] |
!                                                                              |
!     OUTPUT       dpktCP(kcolp,mzp ) : Reduc. Pot.Temperat.Tendency   [K/X/s] |
!     ^^^^^^       dqv_CP(kcolp,mzp ) : Specific   Humidity Tendency [kg/kg/s] |
!                  dqw_CP(kcolp,mzp ) : cloud dropl.Concent.Tendency [kg/kg/s] |
!                  dqi_CP(kcolp,mzp ) : cloud cryst.Concent.Tendency [kg/kg/s] |
!                                                                              |
!                  dss_CP(kcolp     ) : Snow (convective)   Tendency     [m/s] |
!                  drr_CP(kcolp     ) : Rain (convective)   Tendency     [m/s] |
!                                                                              |
!                  pkt_DY(kcolp,mzp ) : Reduced  Potential Temperature   [K/X] |
!                                                                              |
!                  CAPECP(kcolp     ) : Convective Avail.Potent.Energy         |
!                                                                              |
!     REFER. :  MesoNH MASS FLUX CONVECTION Routine                            |
!     ^^^^^^^^         (Bechtold et al., 2001, QJRMS 127, pp 869-886)          |
!                                                                              |
!   # OPTION : #EW  Energy and Water  Conservation                             |
!   # ^^^^^^^^                                                                 |
!                                                                              |
!     MODIF. HGallee: 18-11-2004: Adaptation to CVAmnh.f90.laurent             |
!     ^^^^^^                      (Argument kensbl of CONVECTION removed)      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY____kkl
      use Mod_PHY_CP_ctr
      use Mod_PHY_CP_dat
      use Mod_PHY_CP_grd
      use Mod_PHY_CP_kkl
      use Mod_PHY_DY_kkl


      IMPLICIT NONE



! Global Variables
! ================

      integer                                       ::  mzc                    !  Nb of levels
      integer                                       ::  kcolc                  !  Nb of columns




! Local  Variables
! ================

      logical                                       ::  ANAtun = .TRUE.        !  Parameterization of anabatic breeze
      real(kind=real8)                              ::  dANA                   !

      integer                                       ::  i                      !  x-Axis Index
      integer                                       ::  j                      !  y-Axis Index
      integer                                       ::  k                      !  Level  Index (from top    to bottom)
      integer                                       ::  ikl                    !  Column Index




! ALLOCATION
! ==========

!          ****************
      call PHY_Atm_CP_ALLOC
!          ****************




! SET UP CONVECTION Dimension: Verification
! =========================================

            IF   (mzc   .ne. mzp  )                               THEN
                write(6,601) mzc  ,mzp
 601            format('    ?!&~@|@[#@#]=!!!, mzc   =',i4,' .NE.  mzp   =',i4)
                STOP
            END IF

            IF   (kcolc .ne. kcolp)                               THEN
                write(6,602) kcolc,kcolp
 602            format('    ?!&~@|@[#@#]=!!!, kcolc =',i4,' .NE.  kcolp =',i4)
                STOP
            END IF


            kfdia0  =   kcolp                                           ! index of the last  column of the vector




! SET UP CONVECTION SWITCHES
! ==========================

        IF (Lod_CP)                                                 THEN
            Odeep  = Lod_CP
        ELSE
            Odeep  = Odeep0
            write(6,*) 'Deep    Convection Switch     set to ',Odeep
        END IF
        IF (Los_CP)                                                 THEN
            Oshal  = Los_CP
        ELSE
            Oshal  = Oshal0
            write(6,*) 'Shallow Convection Switch     set to ',Oshal
        END IF

        IF (dtd_CP.GT.0.)                                           THEN
            pdtCVx = dtd_CP
        ELSE
            pdtCVx = pdtCV0
            write(6,*) 'd(time) Convection CALL       set to ',pdtCVx
        END IF

        IF (t_d_CP.GT.0.)                                           THEN
            PTdcv  = t_d_CP
        ELSE
            PTdcv  = PTdcv0
            write(6,*) 'Deep    Convection Time Scale set to ',PTdcv
        END IF
        IF (t_s_CP.GT.0.)                                           THEN
            PTscv  = t_s_CP
        ELSE
            PTscv  = PTscv0
            write(6,*) 'Shallow Convection Time Scale set to ',PTscv
        END IF




! Mass Flux Scheme: Set Up Time Stepping
! ======================================

            pdtCV  =  pdtCVx
        if (pdtCV .lt.dt__CP)                                     then
            pdtCV  =  dt__CP
            jjtCV0 =   1
        else
            jjtCV0 =  pdtCV  / dt__CP
!  ..       jjtCV0 :  Number of  Convective Steps for 1    Update Step

            pdtCV  =  dt__CP * jjtCV0
!  ..       pdtCV  :  Calibrated Convection                  Time Step

        end if
            iitCV0 =  0




! Set UP Anabatic Breeze Parameterization
! =======================================

       IF     ( Lo_ANA )                                            THEN
          rANA      = 2.0d+0 * vANA / xANA
!  ..     rANA      : Subgrid Mountain Breeze: Horizontal Divergence
!                    (Factor 2 included  for 2 horizontal Directions)


! Simple Tuning
! -------------

        IF(ANAtun)                                                  THEN
         DO ikl = 1,kcolp
          hANA(ikl) = sh__AP(ikl) * 2.0d+0                              !  2.0 is a tuning parameter
         END DO


! More sophisticated
! ------------------
        ELSE
         DO ikl = 1,kcolp
          dANA      = sha_AP(ikl)
          hANA(ikl) = abs(dANA)*max(zer0,dxHOST/xANA-un_1)
!  ..     hANA: D("Subgrid Mountain" Height - "Resolved Mountain" Height)
         END DO
        END IF

       END IF ! Lo_ANA




! Set UP Tendencies
! =================

       IF (it_EXP.EQ.1)                                             THEN
        DO ikl = 1,kcolp
        DO k   = 1,mzp
          dpktCP(ikl,k) = 0.
          dqv_CP(ikl,k) = 0.
          dqw_CP(ikl,k) = 0.
          dqi_CP(ikl,k) = 0.
        END DO
        END DO
       END IF




      return
      end subroutine PHY_Atm_CP_INI
