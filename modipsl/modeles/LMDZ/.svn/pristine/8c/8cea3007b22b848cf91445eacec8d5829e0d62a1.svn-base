      subroutine PHY_Atm_S0_RUN

!------------------------------------------------------------------------------+
!                                                         Mon 17-Jun-2013  MAR |
!   MAR          PHY_Atm_S0_RUN                                                |
!     subroutine PHY_Atm_S0_RUN performs    Insolation Computation             |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Thu 25-Apr-2013      |
!           Last Modification by H. Gallee,               Mon 17-Jun-2013      |
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



! LOCAL VARIABLES
! ===============

      use Mod_Atm_S0_RUN


      IMPLICIT NONE


      integer                                       ::   i, j  !
      integer                                       ::  ikl    !
      integer                                       ::     nj  !
      integer                                       ::  lhc    !

      real(kind=real8)                              ::  dlamm  ! mean long. sun for ma-ja
      real(kind=real8)                              ::    anm  !
      real(kind=real8)                              ::   ranm  !
      real(kind=real8)                              ::   ranv  ! anomalie  vraie                          [radian]
      real(kind=real8)                              ::    anv  ! anomalie  vraie                          [degree]
      real(kind=real8)                              ::    tls  ! longitude vraie                          [degree] 
      real(kind=real8)                              ::  rlam   ! longitude vraie                          [radian]
      real(kind=real8)                              ::    sd   ! sinus of solar declination angle (delta)      [-]
      real(kind=real8)                              ::    cd   ! cosin of solar declination angle (delta)      [-]
      real(kind=real8)                              ::  deltar ! Solar Declination (angle sun at equator) [radian]
      real(kind=real8)                              ::  delta  ! Solar Declination (angle sun at equator) [degree]

      real(kind=real8)                              ::  ddt    !
      real(kind=real8)                              ::  arg    !
      real(kind=real8)                              ::  et     !
      real(kind=real8)                              ::  c1     !
      real(kind=real8)                              ::  c2     !
      real(kind=real8)                              ::  c3     !
      real(kind=real8)                              ::  s1     !
      real(kind=real8)                              ::  s2     !
      real(kind=real8)                              ::  s3     !

      real(kind=real8)                              ::  argc   !
      real(kind=real8)                              ::  ahc    !

      real(kind=real8)                              ::  LenDay !  Day     Length                               [h]
      real(kind=real8)                              ::  SunRis !  Time of Sun Rise                             [h]
      real(kind=real8)                              ::  SunSet !  Time of Sun Set                              [h]
      real(kind=real8)                              ::  timh   !
      real(kind=real8)                              ::  ahor   !  Time Angle                                [hour]
      real(kind=real8)                              ::  ahorr  !  Time Angle                              [radian]
      real(kind=real8)                              ::  chor   !
      real(kind=real8)                              ::  zenitr !                                                   RUN

      real(kind=real8)                              ::  omesun !  Sun       Azimuth                       [radian] RUN

      real(kind=real8)                              ::  comes  !
      real(kind=real8)                              ::  somes  !
      real(kind=real8)                              ::  omeswr !

      real(kind=real8)                              ::  r_azim !  azimuth         , fine resolution                ALL

      integer                                       ::  k_azim !                                                   ALL
      integer                                       ::  knazim !                                                   RUN




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                  !
! ALLOCATION
! ==========

      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                        THEN !
          allocate ( cszMIN(kcolp) )                               !  MIN    cos(zenith.Dist), if Mountain MASK        RUN
          allocate ( sin_sz(kcolp) )                               !     sin(Sun Zenith.Dist)                          RUN
          allocate ( sin_ho(kcolp) )                               !     sin(Sun Hour  Angle)                          RUN
      ENDIF
!                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! Insolation at the Top of the Atmosphere (TIME       PARAMETERS)
! ===============================================================

! Solar declination : delta
! -------------------------

      nj      = Day_TU+njYear(Mon_TU)                                   ! Nb of Days since Begin of Year           [day]

      dlamm   = xlam  +      (nj-80) * step                             ! mean long. sun for ma-ja
        anm   = dlamm -       xl
       ranm   =   anm *       Dg2Rad
       ranv   =  ranm +(2.0*ecc-xe3/4.0)*sin(ranm)                     &! anomalie  vraie                       [radian]
     &             +5.0/4.0*ecc*ecc     *sin(2.0*ranm)                 &!
     &            +13.0/12.0   *xe3     *sin(3.0*ranm)                  !
      anv     =  ranv /       Dg2Rad                                    ! anomalie  vraie                       [degree]
      tls     =   anv +       xl                                        ! longitude vraie                       [degree] 
      rlam    =   tls *       Dg2Rad                                    ! longitude vraie                       [radian]

      sd      =    so *   sin(rlam)                                     ! sinus of solar declination angle (delta)   [-]
      cd      =          sqrt(un_1 - sd * sd)                           ! cosin of solar declination angle (delta)   [-]
                                                                        ! sinus delta = sin(obl)*sin(lambda)
                                                                        !                       with lambda = real longitude
                                                                        !(Phd.thesis Marie-France Loutre, ASTR-UCL, Belgium, 1993)
      deltar  =          atan(       sd / cd)                           !
      delta   = deltar/Dg2Rad                                           ! Solar Declination (degrees, angle sun at equator)


! Eccentricity Effect
! -------------------

      dST_UA  =(1.0-ecc*ecc)/(1.0+ecc*cos(ranv))
      ddt     = 1.0/dST_UA                                              ! 1 / normalized earth's sun distance


! Insolation normal to the atmosphere (W/m2)
! ------------------------------------------

      rsunS0  = ddt *ddt *1368.


! Time Equation (Should maybe be modified in case other than present 
! -------------  conditions are used, minor impact)

      arg = om*nj
      c1  = cos(     arg)
      c2  = cos(2.d0*arg)
      c3  = cos(3.d0*arg)
      s1  = sin(     arg)
      s2  = sin(2.d0*arg)
      s3  = sin(3.d0*arg)

      et=0.0072d0*c1 -0.0528d0*c2 -0.0012d0*c3                         &! difference between true solar and mean solar hour angles
     &  -0.1229d0*s1 -0.1565d0*s2 -0.0041d0*s3                          ! (connected to the earth orbital rotation speed)




! Insolation at the Top of the Troposphere (Auxiliary Variables)
! ==============================================================


! Day Length, Time Sunrise and Sunset at Sounding Grid Point (i_x0,j_y0)
! ----------------------------------------------------------------------

      IF (LDayS0)                                                   THEN
          i = i_x0
          j = j_y0

               argc   = -tan(lat__r(ikl0))*tan(deltar)
       if (abs(argc).gt. 1.d0)                                      then
           ahc        =  0.d0
        if    (argc .gt. 1.d0)                                      then
           lhc       =  -1
           LenDay    =   0.d0                                           ! Polar  Night
        else
           lhc       =   1
           LenDay    =  24.d0                                           ! Midnight Sun
        end if
           SunRis    =   0.d0
           SunSet    =   0.d0
       else
           ahc       =  acos(argc)
           lhc       =   0

        if(ahc.lt.0.d0) ahc  =  -ahc
           ahc       =           ahc * 12.d0 / piNmbr
           LenDay    =           ahc *  2.d0
           SunRis    =  12.d0  - ahc - et 
           SunSet    =  SunRis + LenDay
       end if

      END IF

! Time Angle
! ----------

         timh  = HourTU + MinuTU      / 60.d0                           ! time         [hour]

      DO ikl=1,kcolp
         ahor  = timh   + lon__h(ikl) - 12.d0 - et                      ! time angle   [hour]
         ahorr = ahor   * piNmbr      / 12.d0                           ! time angle [radian]
         chor  = cos(ahorr)




! Solar Zenithal Distance zenitr (radians) and
! Insolation (W/m2) at the Atmosphere Top  ===
! =======================================

         csz0S0(ikl) =      sinLat(ikl) *sd                             &
     &                    + cosLat(ikl) *cd *chor
! Martin modif to avoid cos(sza)=0 for LMDZ:
!          csz0SV(ikl) = max(1E-6,csz0S0(ikl))
         csz0S0(ikl) =  max(csz0S0(ikl) ,zer0)
         csz_S0(ikl) =      csz0S0(ikl)

! Martin control
!        PRINT*,'csz0S0(',ikl,')=',csz0S0(ikl)
! Martin control

      ENDDO



! Sun   Azimuth
! =============

      IF   (FaceS0)                                                 THEN
        DO ikl = 1,kcolp


! Slope Impact
! ------------

                  zenitr      = acos(csz0S0(ikl))                       !  Solar Zenithal Distance [radians]
                  sin_sz(ikl) =  sin(zenitr)
                  sin_ho(ikl) =  sin(ahorr)
                  sin_sz(ikl) =  max(eps6      ,sin_sz(ikl))

                  comes       =(sd-sinLat(ikl) *csz0S0(ikl))           &!  Cosine of Sun Azimuth
     &                        /(   cosLat(ikl) *sin_sz(ikl))            !
                  somes       =(cd*sin_ho(ikl)) /sin_sz(ikl)            !    Sine of Sun Azimuth

          IF (abs(comes).gt.zer0)                                   THEN
                  omesun = atan(somes/comes)
            IF   (comes .lt.zer0)                                      &
     &            omesun =  omesun + piNmbr

            IF   (omesun.gt.         piNmbr)                           &
     &            omesun =  -2.0d0 * piNmbr + omesun
            IF   (omesun.lt.        -piNmbr)                           &
     &            omesun =   2.0d0 * piNmbr + omesun

          ELSE
            IF   (somes .gt.zer0)                                   THEN
                  omesun =   0.5d0 * piNmbr
            ELSE
                  omesun =   1.5d0 * piNmbr
            END IF
          END IF

          IF (ikl.eq.ikl0)                                             &
      &           omeswr =  omesun / Dg2Rad
                  omesun =  -2.0d0 * piNmbr + omesun + DirAxX * Dg2Rad  !  Sun Azimuth (in MAR Reference Frame)
                                                                        !           (positive counterclockwise)

! Minimum Zenithal Distance
! ~~~~~~~~~~~~~~~~~~~~~~~~~
                  cszMIN(ikl) = 0.0d00
          IF       (MMskS0)                                         THEN
                                r_azim = omesun  / d_azim
                                k_azim =       int(r_azim)

! Mountain S0 MASKING
! ~~~~~~~~~~~~~~~~~~~
           IF(k_azim.le.     0)                                     THEN
                                r_azim = r_azim  + n_azim
                                k_azim = k_azim  + n_azim
           END IF
                                knazim = k_azim  + 1
           IF(knazim.gt.n_azim) knazim = knazim  - n_azim

                  cszMIN(ikl) = cszkS0(ikl,k_azim)+(r_azim  - k_azim)  &
     &                        *(cszkS0(ikl,knazim)-cszkS0(ikl,k_azim))
          END IF ! (MMskS0)

! Cosine of Solar Normal Angle
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  csz_S0(ikl) = slopS0(ikl)            *csz0S0(ikl)    &
     &          + sin_sz(ikl) * slopS0(ikl)*sqrt(un_1  -slopS0(ikl))   &
     &                        *              cos(omesun-omenS0(ikl))
                  csz_S0(ikl) =              max(zer0  ,csz_S0(ikl))
          IF     (csz0S0(ikl).le.cszMIN(ikl))                          &
     &            csz_S0(ikl) = 0.0d00
        END DO
      END IF     ! (FaceS0)




! Output Insolation Characteristics
! =================================

      IF (MinuTU .eq. 0 .and. Sec_TU .eq. 0)                           THEN

         ahor  = timh   + lon__h(ikl0)  - 12.d0 - et
         zenitr=     acos(csz0S0(ikl0)) / Dg2Rad
! #AZ    anormr=     acos(csz_S0(ikl0)) / Dg2Rad
! #AZ    omenwr= DirAxX - omenS0(ikl0)  / Dg2Rad
! #AZ    if (omenwr.lt.  0.) omenwr = omenwr + 360.d0
! #AZ    if (omenwr.gt.360.) omenwr = omenwr - 360.d0
! #AZ                        omeswr = 360.d0 - omeswr
! #AZ    if (omeswr.lt.  0.) omeswr = omeswr + 360.d0
! #AZ    if (omeswr.gt.360.) omeswr = omeswr - 360.d0

        write(4,1)           lat__r(ikl0)/Dg2Rad,lon__h(ikl0)           &
     &            ,Day_TU,Mon_TU,HourTU,MinuTU,Sec_TU
 1      format(/,' lat.=',f6.1,3x,'long.=',f7.1,4x,'date :',i3,'-',i2,  &
     &           ' / ',i2,' h.UT',i3,' min.',i3,' sec.')
        write(4,2) i_x0,j_y0,lat__r(ikl0)/Dg2Rad,lon__h(ikl0)
 2      format(' Sounding at (',i3,i3,') / (',f6.2,'dg,',f6.2,'ho)')
        write(4,3) rsunS0*csz_S0(ikl0)     ,ahor,zenitr                 &
! #AZ&            ,omeswr,omenwr           ,     anormr                 &
     &            ,delta
 3      format(' Insolation [W/m2]  = ',f7.2,'   Hor.Angle = ',f7.2,    &
     &         '   Zenith.Angle = ',f7.2                                &
! #AZ&      ,/,'                      ', 7x ,'   Sol.Azim. = ',f7.2     &
! #AZ&      ,/,'                      ', 7x ,'   Nrm.Azim. = ',f7.2     &
! #AZ&        ,'   Normal Angle = ',f7.2                                &
     &      ,/,' Solar Declination  = ',f7.2)

      IF       (LDayS0)                                             THEN
       if (lhc.eq.-1)                                                   &
     &  write(4,4) SunRis,LenDay,SunSet
 4      format(' Sun Rise Time [h]  = ',f7.2,'   Day Leng. = ',f7.2,    &
     &         '   Sun Set Time = ',f7.2,'  -- POLAR  NIGHT --')
       if (lhc.eq. 0)                                                   &
     &  write(4,5) SunRis,LenDay,SunSet
 5      format(' Sun Rise Time [h]  = ',f7.2,'   Day Leng. = ',f7.2,    &
     &         '   Sun Set Time = ',f7.2,'  -- SOLAR  TIME  --')
       if (lhc.eq. 1)                                                   &
     &  write(4,6) SunRis,LenDay,SunSet
 6      format(' Sun Rise Time [h]  = ',f7.2,'   Day Leng. = ',f7.2,    &
     &         '   Sun Set Time = ',f7.2,'  -- MIDNIGHT SUN --')
      END IF ! (LDayS0)

      END IF ! 




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                  !
! DE-ALLOCATION
! =============

      IF (FlagDALLOC)                                         THEN !
          deallocate ( cszMIN )                                    !  Sum of cos(zenith.Dist) on Azimut Intervals      RUN
          deallocate ( sin_sz )                                    !     sin(Sun zenith.Dist)                          RUN
          deallocate ( sin_ho )                                    !     sin(Sun Hour  Angle)                          RUN
      ENDIF
!                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      return
      end subroutine PHY_Atm_S0_RUN
