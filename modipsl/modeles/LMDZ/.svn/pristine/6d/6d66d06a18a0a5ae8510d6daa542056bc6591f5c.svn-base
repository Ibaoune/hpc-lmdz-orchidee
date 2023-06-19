      subroutine PHY_Atm_CM_INI

!------------------------------------------------------------------------------+
!                                                         Sat 29-Jun-2013  MAR |
!   MAR          PHY_Atm_CM_INI                                                |
!     subroutine PHY_Atm_CM_INI initializes Cloud Microphysical Scheme CMiPhy  |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Thu 21-Mar-2013      |
!           Last Modification by H. Gallee,               Sat 29-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+
!                                                                              |
!                                                                              |
! # OPTIONS:                                                                   |
! # ^^^^^^^                                                                    |
! #          #kk  SCu  fraction  Limitation                                    |
! #          #hb  Snow      particles distrib. parameter  n0___s is BSnow value|
! #          #LA  Snow,Rain particles distrib. parameters n0___s,r        TUNED|
! #          #rc  MAX  allowed Relative Humidity is a function of Grid Size    |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____grd
      use Mod_PHY_CM_dat
      use Mod_PHY_CM_kkl


      IMPLICIT NONE

      integer        ::  i     ,j     ,k     ,ikl



!  DATA TUNING
!  ===========

! #kk SSImax = 900.0        ! Maximum        Sursaturation % ICE (900 ==> RH=1000%)       [%]
! #LA n0___r = 3.0e06       ! intercept parameter / rain gamam distribution ! Tuning
! #LA n0___s = 4.0e06       ! intercept parameter / snow gamma distribution ! Tuning
! #hb n0___s = 0.1d18       ! intercept parameter / snow gamma distribution (Blown Snow only)
                            ! DO NOT USE unless for specific sensivity experiments
! #hb IF (it_EXP.eq.1) write(6,6000)
! #hb 6000 format(/,' ****************************************************', &
! #hb&            /,' * cnos  = 0.1d18 for PURE BLOWING SNOW EXPERIMENTS *', &
! #hb&            /,' *             DO not USE  OTHERWISE                *', &
! #hb&            /,' ****************************************************', &
! #hb&            /)

! _hl qisMAX = 0.0008       ! compromise when graupels are not included
!                           ! Ref.: Emde & Kahlig 1989, Ann.Geoph.      7, p.408  (18) 
      qw_MAX = qw_MAXL      ! Cloud droplets MAX concentration before autoconversion
                            ! Ref.: Lin et al.    1983, JCAM           22, p.1076 (50)
!     qw_MAX = 0.0003       ! critical liquid water mixing ratio
!                           ! Ref.: Sundqvist     1988, Physically-Based Modelling and Simulation of Climate and Climatic Change,
!                           !                           M.E. Schlesinger, Ed., Reidel, 433-461.



! Maximum Relative Humidity
! =========================

! #rc RH_MAX = 0.90+0.08*sqrt(max(0.,100.-dxHOST*0.001)/95.)




! IOs
! ===

        DO  ipt_CM=1,npt_CM
            ikl0CM(ipt_CM) = ikl_AP(i0__CM(ipt_CM),j0__CM(ipt_CM))
        ENDDO




! Cloud Microphysics Initialization
! =================================

      IF (it_EXP.LE.1)                                              THEN


! Precipitation: Mod_PHY_CM_kkl:
! -----------------------------

          DO ikl=1,kcolp
            rai0CM    (ikl)   = 0.
            rainCM    (ikl)   = 0.
            sno0CM    (ikl)   = 0.
            snobCM    (ikl)   = 0.
            snowCM    (ikl)   = 0.
            Ice0CM    (ikl)   = 0.
            Ice_CM    (ikl)   = 0.



! Hyrometeores : Mod_PHY_CM_kkl:
! -----------------------------

          DO k=   1,mzp
            CCNiCM    (ikl,k) = 0.
            qw__CM    (ikl,k) = 0.
            qi__CM    (ikl,k) = 0.
            qs__CM    (ikl,k) = 0.
! #qg       qg__CM    (ikl,k) = 0.
            qr__CM    (ikl,k) = 0.
          END DO
          END DO

      END IF

      end subroutine PHY_Atm_CM_INI
