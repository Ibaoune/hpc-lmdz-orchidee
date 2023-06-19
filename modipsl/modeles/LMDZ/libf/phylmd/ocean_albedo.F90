!
! $Id$
!

SUBROUTINE ocean_albedo(knon,zrmu0,knindex,pwind,SFRWL,alb_dir_new,alb_dif_new)
!!
!!****  *ALBEDO_RS14*  
!!
!!    PURPOSE
!!    -------
!!     computes the direct & diffuse albedo over open water
!!     
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!	R. Séférian           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    03/2014
!!                  05/2014 R. Séférian & B. Decharme :: Adaptation to spectral
!!                  computation for diffuse and direct albedo
!!                  08/2014 S. Baek :: for wider wavelength range 200-4000nm and
!!                  adaptation to LMDZ + whitecap effect by Koepke + chrolophyll
!!                  map from climatology file
!!                  10/2016 O. Boucher :: some optimisation following R.
!!                  Seferian's work in the CNRM Model
!!       
!-------------------------------------------------------------------------------
!
!*           DECLARATIONS
!            ------------
!
USE ocean_albedo_para
USE dimphy
USE phys_state_var_mod, ONLY : chl_con
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!              -------------------------
!
include "clesphys.h"
!
INTEGER, INTENT(IN) :: knon
INTEGER, DIMENSION(klon), INTENT(IN) :: knindex
REAL, DIMENSION(klon), INTENT(IN) :: zrmu0         !--cos(SZA) on full vector
REAL, DIMENSION(klon), INTENT(IN) :: pwind         !--wind speed on compressed vector
REAL, DIMENSION(6),INTENT(IN) :: SFRWL
REAL, DIMENSION(klon,nsw), INTENT(OUT) :: alb_dir_new, alb_dif_new
!
!*      0.2    declarations of local variables
!              -------------------------
!
REAL, DIMENSION(klon)           :: ZCHL        ! surface chlorophyll
REAL, DIMENSION(klon)           :: ZCOSZEN     ! Cosine of the zenith solar angle
!
INTEGER                         :: JWL, INU    ! indexes
INTEGER                         :: JI
REAL                            :: ZWL         ! input parameter: wavelength and diffuse/direct fraction of light
REAL:: ZCHLABS, ZAW, ZBW, ZREFM, ZYLMD, ZUE, ZUE2 ! scalar computation variables
!
REAL, DIMENSION(klon) :: ZAP, ZXX2, ZR00, ZRR0, ZRRR               ! computation variables
REAL, DIMENSION(klon) :: ZR22, ZR11DF                              ! computation variables
REAL, DIMENSION(klon) :: ZBBP, ZNU, ZHB                            ! computation variables
REAL, DIMENSION(klon) :: ZR11, ZRW, ZRWDF, ZRDF                    ! 4 components of the OSA
REAL, DIMENSION(klon) :: ZSIG, ZFWC, ZWORK1, ZWORK2, ZWORK3
! 
!--initialisations-------------
!

IF (knon==0) RETURN ! A verifier pourquoi on en a besoin...

alb_dir_new(:,:) = 0. 
alb_dif_new(:,:) = 0. 
!
! Initialisation of chlorophyll content
! ZCHL(:) = CHL_CON!0.05 ! averaged global values for surface chlorophyll
IF (ok_chlorophyll) THEN
  ZCHL(1:knon)=CHL_CON(knindex(1:knon))
ELSE 
  ZCHL(1:knon) = 0.05
ENDIF

! variables that do not depend on wavelengths
! loop over the grid points
! functions of chlorophyll content
ZWORK1(1:knon)= EXP(LOG(ZCHL(1:knon))*0.65)
ZWORK2(1:knon)= 0.416 * EXP(LOG(ZCHL(1:knon))*0.766)
ZWORK3(1:knon)= LOG10(ZCHL(1:knon))
! store the cosine of the solar zenith angle
ZCOSZEN(1:knon) = zrmu0(knindex(1:knon))
! Compute sigma derived from wind speed (Cox & Munk reflectance model)
ZSIG(1:knon)=SQRT(0.003+0.00512*PWIND(1:knon))
! original : correction for foam (Eq 16-17)
! has to be update once we have information from wave model (discussion with G. Madec)
ZFWC(1:knon)=3.97e-4*PWIND(1:knon)**1.59 ! Salisbury 2014 eq(2) at 37GHz, value in fraction
!
DO JWL=1,NNWL           ! loop over the wavelengths
  !
  !---------------------------------------------------------------------------------
  ! 0- Compute baseline values
  !---------------------------------------------------------------------------------
    
  ! Get refractive index for the correspoding wavelength
  ZWL=XAKWL(JWL)      !!!--------- wavelength value
  ZREFM= XAKREFM(JWL) !!!--------- refraction index value
  
  !---------------------------------------------------------------------------------
  ! 1- Compute direct surface albedo (ZR11)
  !---------------------------------------------------------------------------------
  !
  ZXX2(1:knon)=SQRT(1.0-(1.0-ZCOSZEN(1:knon)**2)/ZREFM**2)
  ZRR0(1:knon)=0.50*(((ZXX2(1:knon)-ZREFM*ZCOSZEN(1:knon))/(ZXX2(1:knon)+ZREFM*ZCOSZEN(1:knon)))**2 +  & 
               ((ZCOSZEN(1:knon)-ZREFM*ZXX2(1:knon))/(ZCOSZEN(1:knon)+ZREFM*ZXX2(1:knon)))**2)
  ZRRR(1:knon)=0.50*(((ZXX2(1:knon)-1.34*ZCOSZEN(1:knon))/(ZXX2(1:knon)+1.34*ZCOSZEN(1:knon)))**2 + & 
               ((ZCOSZEN(1:knon)-1.34*ZXX2(1:knon))/(ZCOSZEN(1:knon)+1.34*ZXX2(1:knon)))**2)
  ZR11(1:knon)=ZRR0(1:knon)-(0.0152-1.7873*ZCOSZEN(1:knon)+6.8972*ZCOSZEN(1:knon)**2-8.5778*ZCOSZEN(1:knon)**3+ & 
               4.071*ZSIG(1:knon)-7.6446*ZCOSZEN(1:knon)*ZSIG(1:knon)) *  & 
               EXP(0.1643-7.8409*ZCOSZEN(1:knon)-3.5639*ZCOSZEN(1:knon)**2-2.3588*ZSIG(1:knon)+ & 
               10.0538*ZCOSZEN(1:knon)*ZSIG(1:knon))*ZRR0(1:knon)/ZRRR(1:knon)
  ! 
  !---------------------------------------------------------------------------------
  ! 2- Compute surface diffuse albedo (ZRDF)
  !---------------------------------------------------------------------------------
  ! Diffuse albedo from Jin et al., 2006 + estimation from diffuse fraction of
  ! light (relying later on AOD). CNRM model has opted for Eq 5b
  ZRDF(1:knon)=-0.1482-0.012*ZSIG(1:knon)+0.1609*ZREFM-0.0244*ZSIG(1:knon)*ZREFM ! surface diffuse (Eq 5a)
  !!ZRDF(1:knon)=-0.1479+0.1502*ZREFM-0.0176*ZSIG(1:knon)*ZREFM   ! surface diffuse (Eq 5b) 
 
  !---------------------------------------------------------------------------------
  ! *- Determine absorption and backscattering
  ! coefficients to determine reflectance below the surface (Ro) once for all
  !
  ! *.1- Absorption by chlorophyll
  ZCHLABS= XAKACHL(JWL) 
  ! *.2- Absorption by seawater 
  ZAW= XAKAW3(JWL) 
  ! *.3- Backscattering by seawater
  ZBW= XAKBW(JWL) 
  ! *.4- Backscattering by chlorophyll
  ZYLMD = EXP(0.014*(440.0-ZWL))
  ZAP(1:knon) = 0.06*ZCHLABS*ZWORK1(1:knon) +0.2*(XAW440+0.06*ZWORK1(1:knon))*ZYLMD
   
!!  WHERE ( ZCHL(1:knon) > 0.02 )
!!    ZNU(:)=MIN(0.0,0.5*(ZWORK3(:)-0.3))
!!    ZBBP(:)=(0.002+0.01*(0.5-0.25*ZWORK3(:))*(ZWL/550.)**ZNU(:))*ZWORK2(:)
!!  ELSEWHERE
!!    ZBBP(:)=0.019*(550./ZWL)*ZWORK2(:)       !ZBBPf=0.0113 at chl<=0.02
!!  ENDWHERE

    do JI = 1, knon
      IF (ZCHL(JI) > 0.02) THEN
        ZNU(JI)=MIN(0.0,0.5*(ZWORK3(JI)-0.3))
        ZBBP(JI)=(0.002+0.01*(0.5-0.25*ZWORK3(JI))*(ZWL/550.)**ZNU(JI)) &
                  *ZWORK2(JI)
      ELSE
        ZBBP(JI)=0.019*(550./ZWL)*ZWORK2(JI)       !ZBBPf=0.0113 at chl<=0.02 
      ENDIF
    ENDDO

  ! Morel-Gentili(1991), Eq (12)
  ! ZHB=h/(h+2*ZBBPf*(1.-h))        
  ZHB(1:knon)=0.5*ZBW/(0.5*ZBW+ZBBP(1:knon))
   
  !---------------------------------------------------------------------------------
  ! 3- Compute direct water-leaving albedo (ZRW)
  !---------------------------------------------------------------------------------
  ! Based on Morel & Gentilli 1991 parametrization
  ZR22(1:knon)=0.48168549-0.014894708*ZSIG(1:knon)-0.20703885*ZSIG(1:knon)**2

  ! Use Morel 91 formula to compute the direct reflectance
  ! below the surface
  ZR00(1:knon)=(0.5*ZBW+ZBBP(1:knon))/(ZAW+ZAP(1:knon)) *  & 
               (0.6279-0.2227*ZHB(1:knon)-0.0513*ZHB(1:knon)**2 + & 
               (-0.3119+0.2465*ZHB(1:knon))*ZCOSZEN(1:knon))
  ZRW(1:knon)=ZR00(1:knon)*(1.-ZR22(1:knon))/(1.-ZR00(1:knon)*ZR22(1:knon))

  !---------------------------------------------------------------------------------
  ! 4- Compute diffuse water-leaving albedo (ZRWDF)
  !---------------------------------------------------------------------------------
  ! as previous water-leaving computation but assumes a uniform incidence of
  ! shortwave at surface (ue)
  ZUE=0.676               ! equivalent u_unif for diffuse incidence
  ZUE2=SQRT(1.0-(1.0-ZUE**2)/ZREFM**2)
  ZRR0(1:knon)=0.50*(((ZUE2-ZREFM*ZUE)/(ZUE2+ZREFM*ZUE))**2 +((ZUE-ZREFM*ZUE2)/(ZUE+ZREFM*ZUE2))**2)
  ZRRR(1:knon)=0.50*(((ZUE2-1.34*ZUE)/(ZUE2+1.34*ZUE))**2 +((ZUE-1.34*ZUE2)/(ZUE+1.34*ZUE2))**2)
  ZR11DF(1:knon)=ZRR0(1:knon)-(0.0152-1.7873*ZUE+6.8972*ZUE**2-8.5778*ZUE**3+4.071*ZSIG(1:knon)-7.6446*ZUE*ZSIG(1:knon)) * &
                 EXP(0.1643-7.8409*ZUE-3.5639*ZUE**2-2.3588*ZSIG(1:knon)+10.0538*ZUE*ZSIG(1:knon))*ZRR0(1:knon)/ZRRR(1:knon)

  ! Use Morel 91 formula to compute the diffuse
  ! reflectance below the surface
  ZR00(1:knon) = (0.5*ZBW+ZBBP(1:knon)) / (ZAW+ZAP(1:knon)) &
       * (0.6279-0.2227*ZHB(1:knon)-0.0513*ZHB(1:knon)**2 &
       + (-0.3119+0.2465*ZHB(1:knon))*ZUE)
  ZRWDF(1:knon)=ZR00(1:knon)*(1.-ZR22(1:knon))*(1.-ZR11DF(1:knon))/(1.-ZR00(1:knon)*ZR22(1:knon))
   
  ! get waveband index inu for each nsw band
  SELECT CASE(nsw)
  CASE(2)
    IF (JWL.LE.49) THEN       ! from 200  to 680 nm 
     inu=1
    ELSE                      ! from 690  to 4000 nm
     inu=2
    ENDIF
  CASE(4)
    IF (JWL.LE.49) THEN       ! from 200  to 680 nm 
     inu=1
    ELSE IF (JWL.LE.99) THEN  ! from 690  to 1180 nm
     inu=2
    ELSE IF (JWL.LE.218) THEN ! from 1190 to 2370 nm
     inu=3
    ELSE                      ! from 2380 to 4000 nm
     inu=4
    ENDIF
  CASE(6)
    IF (JWL.LE.5) THEN        ! from 200  to 240 nm 
     inu=1
    ELSE IF (JWL.LE.24) THEN  ! from 250  to 430 nm
     inu=2
    ELSE IF (JWL.LE.49) THEN  ! from 440  to 680 nm
     inu=3
    ELSE IF (JWL.LE.99) THEN  ! from 690  to 1180 nm
     inu=4
    ELSE IF (JWL.LE.218) THEN ! from 1190 to 2370 nm
     inu=5
    ELSE                      ! from 2380 to 4000 nm
     inu=6
    ENDIF
  END SELECT

  ! partitionning direct and diffuse albedo
  ! excluding diffuse albedo ZRW on ZDIR_ALB

  !--direct
  alb_dir_new(1:knon,inu)=alb_dir_new(1:knon,inu) + & 
                          ( XFRWL(JWL) * ((1.-ZFWC(1:knon)) * (ZR11(1:knon)+ZRW(1:knon))   + ZFWC(1:knon)*XRWC(JWL)) )/SFRWL(inu)
  !--diffuse
  alb_dif_new(1:knon,inu)=alb_dif_new(1:knon,inu) + & 
                          ( XFRWL(JWL) * ((1.-ZFWC(1:knon)) * (ZRDF(1:knon)+ZRWDF(1:knon)) + ZFWC(1:knon)*XRWC(JWL)) )/SFRWL(inu)

ENDDO ! ending loop over wavelengths

END SUBROUTINE ocean_albedo
