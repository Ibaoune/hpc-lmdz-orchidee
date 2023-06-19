      subroutine PHY_Atm_S0_INI

!------------------------------------------------------------------------------+
!                                                         Sat 15-Jun-2013  MAR |
!   MAR          PHY_Atm_S0_INI                                                |
!     subroutine PHY_Atm_S0_INI initializes Insolation Computation             |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Thu 25-Apr-2013      |
!           Last Modification by H. Gallee,               Sat 15-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+
!                                                                              |
!     REFER.:   Ch.  Tricot, personal communication                            |
!     ^^^^^^^   M.F. Loutre, personal communication and thesis (1993)          |
!                                                                              |
!     INPUT :   Mon_TU, Day_TU        : Month and Day of the Year              |
!     ^^^^^^^   HourTU, MinuTU, Sec_TU: Hour, Minute, and Second               |
!               lon__r(kcolp) : Latitude                             [radians] |
!               lat__h(kcolp) : Longitude                              [hours] |
!                                                                              |
!     OUTPUT:   rsunS0        : Insolation normal to Atmosphere Top (W/m2)     |
!     ^^^^^^^   csz0S0        : Cosinus of the Zenithal Distance               |
!                                                                              |
!------------------------------------------------------------------------------+


! Global Variables
! ================

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY____kkl
      use Mod_PHY_S0_ctr
      use Mod_PHY_S0_dat
      use Mod_PHY_S0_grd
      use Mod_PHY_S0_kkl


      IMPLICIT NONE




! LOCAL VARIABLES
! ===============

      integer                                              ::  ikl    !
      real(kind=real8)                                     ::  omenor !  Fall Line Azimuth                      [radian] INI




! ALLOCATION
! ==========

!          ****************
      CALL PHY_Atm_S0_ALLOC
!          ****************


! INITIALIZATION
! ==============


! Insolation Parameters (kBP is time in kyear BP, and is prescribed in Mod_PHY_S0_dat)
! ---------------------

        ecc    =    ecc_EO(kBP)                                         ! Earth Orbit Eccentricity                   [-]
        perh   =    perhEO(kBP)                                         ! Longitude of Perihelion               [degree]
        xob    =    xob_EO(kBP)                                         ! Obliquity                             [degree]

        pirr   =    Dg2Rad / 3600.

        xee    =      ecc  * ecc                                        ! Square    of Eccentricity                  [-]
        xse    = sqrt(un_1 - xee)                                       ! Square Root (1-Ecc**2)                     [-]
        xe3    =      xee  * ecc                                        !
        xl     =      perh +  180.0                                     ! Longitude of Aphelion                 [degree]
        xllp   =      xl   * Dg2Rad                                     ! Longitude of Aphelion                 [radian]
        so     =  sin(xob  * Dg2Rad)                                    ! sinus     of Obliquity

        xlam  =(ecc/2.0+ecc*xee/8.0)*(1.0    +xse)*sin(    xllp)       &! true long. sun for mean long. = 0     [radian]
     &        -(xee/4.0            )*(0.5    +xse)*sin(2.0*xllp)       &!
     &         +ecc*xee/8.0         *(1.0/3.0+xse)*sin(3.0*xllp)        !
        xlam  = 2.0*xlam/Dg2Rad                                         ! true long. sun for mean long. = 0     [degree]

        step   = 360.00 /Tyear                                          ! Advance on Earth Orbit in 1 day   [degree/day]




! Initialisation of Slope effect on Insolation: Slope Azimuth
! ===========================================================

        IF         (FaceS0)                                         THEN

          DO ikl=1,kcolp

                    slopS0(ikl) =      cos(atan(slopAP(ikl)))           !  Cosine of Fall Line Angle

            IF (abs(sloxAP(ikl)).gt.zer0)                           THEN
                    omenor = atan(sloyAP(ikl) / sloxAP(ikl))            !  Fall Line Azimuth   (Upslope Direction)
              IF   (sloxAP(ikl) .lt.zer0)                              &!
     &              omenor =  omenor + piNmbr                           !

              IF   (omenor.gt.         piNmbr)                         &!
     &              omenor =  -2.0d0 * piNmbr + omenor                  !
              IF   (omenor.lt.        -piNmbr)                         &!
     &              omenor =   2.0d0 * piNmbr + omenor                  !

            ELSE
              IF   (sloyAP(ikl).gt.zer0) then
                    omenor =   0.5d0 * piNmbr
              ELSE
                    omenor =   1.5d0 * piNmbr
              END IF
            END IF

                    omenS0(ikl) =      omenor - piNmbr                  !  Fall Line Azimuth (Downslope Direction)
                                                                        !                 (in MAR Reference Frame)
                                                                        !              (positive counterclockwise)
          END DO

        END IF


      return
      end subroutine PHY_Atm_S0_INI
