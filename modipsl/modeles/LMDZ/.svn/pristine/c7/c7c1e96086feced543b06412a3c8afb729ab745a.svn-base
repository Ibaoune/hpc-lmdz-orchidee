MODULE YOMCST_SISVAT

IMPLICIT NONE
                !
! $Header$
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez à n'utiliser que des ! pour les commentaires
!                 et à bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!
! A1.0 Fundamental constants
      REAL,PARAMETER :: RPI=3.141592653589793238462643e0,               &
     &                  RCLUM=299792458.,   RHPLA=6.6260755E-34,        &
     &                  RKBOL=1.380658E-23, RNAVO=6.0221367E+23

! A1.1 Astronomical constants
      REAL,PARAMETER :: RDAY=86400.,REA=149597870000.,REPSM=0.409093,   &
     &        RSIYEA=365.25*RDAY*2.*RPI/6.283076,                       &
     &        RSIDAY=RDAY/(1.+RDAY/RSIYEA), ROMEGA=2.*RPI/RSIDAY

! A1.1.bis Constantes concernant l'orbite de la Terre:
      REAL,PARAMETER :: R_ecc=0.016715,R_peri=102.7,R_incl=23.441

! A1.2 Geoide
      REAL,PARAMETER :: RG=9.80665, RA=6371229.,R1SA=1.D0/RA

! A1.3 Radiation
      REAL,PARAMETER :: RSIGMA = 2.*RPI**5 * (RKBOL/RHPLA)**3           &
     &                   * RKBOL/RCLUM/RCLUM/15.
!     REAL RSIGMA,RI0

! A1.4 Thermodynamic gas phase
      REAL,PARAMETER :: R=RNAVO*RKBOL, RMD=28.9644,    RMO3=47.9942,    &
     &                  RMV=18.0153,   RD=1000.*R/RMD, RV=1000.*R/RMV,  &
     &                  RCPD=3.5*RD,   RCVD=RCPD-RD,   RCPV=4. *RV,     &
     &                  RCVV=RCPV-RV,  RKAPPA=RD/RCPD, RETV=RV/RD-1.
! A1.5,6 Thermodynamic liquid,solid phases
      REAL,PARAMETER :: RCW=RCPV,      RCS=RCPV

! A1.7 Thermodynamic transition of phase
      REAL,PARAMETER :: RTT=273.16,    RLVTT=2.5008E+6,RLSTT=2.8345E+6, &
     &                  RLMLT=RLSTT-RLVTT,             RATM=100000.

! A1.8 Curve of saturation
      REAL,PARAMETER :: RESTT=611.14,  RGAMW=(RCW-RCPV)/RV,             &
     &                  RBETW=RLVTT/RV+RGAMW*RTT,                       & 
     &                  RGAMS=(RCS-RCPV)/RV,                            &    
     &                  RBETS=RLSTT/RV+RGAMS*RTT,                       &  
     &                  RGAMD=RGAMS-RGAMW, RBETD=RBETS-RBETW
      REAL           :: RALPW,RALPS,RALPD

!      RALPW = LOG(RESTT)+RBETW/RTT+RGAMW*LOG(RTT) 
!      RALPS = LOG(RESTT)+RBETS/RTT+RGAMS*LOG(RTT) 
!      RALPD = RALPS-RALPW 
   

END MODULE YOMCST_SISVAT
