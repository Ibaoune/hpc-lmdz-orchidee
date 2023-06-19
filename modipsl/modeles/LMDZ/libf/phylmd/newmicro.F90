! $Id: newmicro.F90 3281 2018-03-16 18:26:14Z musat $

SUBROUTINE newmicro(flag_aerosol, ok_cdnc, bl95_b0, bl95_b1, paprs, pplay, t, pqlwp, pclc, &
    pcltau, pclemi, pch, pcl, pcm, pct, pctlwp, xflwp, xfiwp, xflwc, xfiwc, &
    mass_solu_aero, mass_solu_aero_pi, pcldtaupi, re, fl, reliq, reice, &
    reliq_pi, reice_pi)

  USE dimphy
  USE phys_local_var_mod, ONLY: scdnc, cldncl, reffclwtop, lcc, reffclws, &
      reffclwc, cldnvi, lcc3d, lcc3dcon, lcc3dstra, icc3dcon, icc3dstra,  & 
      zfice, dNovrN
  USE phys_state_var_mod, ONLY: rnebcon, clwcon
  USE icefrac_lsc_mod ! computes ice fraction (JBM 3/14)
  USE ioipsl_getin_p_mod, ONLY : getin_p
  USE print_control_mod, ONLY: lunout


  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930910
  ! O.   Boucher (LMD/CNRS) mise a jour en 201212
  ! I. Musat (LMD/CNRS) : prise en compte de la meme hypothese de recouvrement
  !                       pour les nuages que pour le rayonnement rrtm via
  !                       le parametre novlp de radopt.h : 20160721
  ! Objet: Calculer epaisseur optique et emmissivite des nuages
  ! ======================================================================
  ! Arguments:
  ! ok_cdnc-input-L-flag pour calculer les rayons a partir des aerosols

  ! t-------input-R-temperature
  ! pqlwp---input-R-eau liquide nuageuse dans l'atmosphere dans la partie
  ! nuageuse (kg/kg)
  ! pclc----input-R-couverture nuageuse pour le rayonnement (0 a 1)
  ! mass_solu_aero-----input-R-total mass concentration for all soluble
  ! aerosols[ug/m^3]
  ! mass_solu_aero_pi--input-R-ditto, pre-industrial value

  ! bl95_b0-input-R-a PARAMETER, may be varied for tests (s-sea, l-land)
  ! bl95_b1-input-R-a PARAMETER, may be varied for tests (    -"-      )

  ! re------output-R-Cloud droplet effective radius multiplied by fl [um]
  ! fl------output-R-Denominator to re, introduced to avoid problems in
  ! the averaging of the output. fl is the fraction of liquid
  ! water clouds within a grid cell

  ! pcltau--output-R-epaisseur optique des nuages
  ! pclemi--output-R-emissivite des nuages (0 a 1)
  ! pcldtaupi-output-R-pre-industrial value of cloud optical thickness,

  ! pcl-output-R-2D low-level cloud cover
  ! pcm-output-R-2D mid-level cloud cover
  ! pch-output-R-2D high-level cloud cover
  ! pct-output-R-2D total cloud cover
  ! ======================================================================

  include "YOMCST.h"
  include "nuage.h"
  include "radepsi.h"
  include "radopt.h"

  ! choix de l'hypothese de recouvrement nuageuse via radopt.h (IM, 19.07.2016)
  ! !novlp=1: max-random
  ! !novlp=2: maximum
  ! !novlp=3: random
! LOGICAL random, maximum_random, maximum
! PARAMETER (random=.FALSE., maximum_random=.TRUE., maximum=.FALSE.)

  LOGICAL, SAVE :: first = .TRUE.
  !$OMP THREADPRIVATE(FIRST)
  INTEGER flag_max

  ! threshold PARAMETERs
  REAL thres_tau, thres_neb
  PARAMETER (thres_tau=0.3, thres_neb=0.001)

  REAL phase3d(klon, klev)
  REAL tcc(klon), ftmp(klon), lcc_integrat(klon), height(klon)

  REAL paprs(klon, klev+1)
  REAL pplay(klon, klev)
  REAL t(klon, klev)
  REAL pclc(klon, klev)
  REAL pqlwp(klon, klev)
  REAL pcltau(klon, klev)
  REAL pclemi(klon, klev)
  REAL pcldtaupi(klon, klev)

  REAL pct(klon)
  REAL pcl(klon)
  REAL pcm(klon)
  REAL pch(klon)
  REAL pctlwp(klon)

  LOGICAL lo

  ! !Abderr modif JL mail du 19.01.2011 18:31
  ! REAL cetahb, cetamb
  ! PARAMETER (cetahb = 0.45, cetamb = 0.80)
  ! Remplacer
  ! cetahb*paprs(i,1) par  prmhc
  ! cetamb*paprs(i,1) par  prlmc
  REAL prmhc ! Pressure between medium and high level cloud in Pa
  REAL prlmc ! Pressure between low and medium level cloud in Pa
  PARAMETER (prmhc=440.*100., prlmc=680.*100.)

  INTEGER i, k
  REAL xflwp(klon), xfiwp(klon)
  REAL xflwc(klon, klev), xfiwc(klon, klev)

  REAL radius

  REAL coef_froi, coef_chau
  PARAMETER (coef_chau=0.13, coef_froi=0.09)

  REAL seuil_neb
  PARAMETER (seuil_neb=0.001)

! JBM (3/14) nexpo is replaced by exposant_glace
! INTEGER nexpo ! exponentiel pour glace/eau
! PARAMETER (nexpo=6)
! PARAMETER (nexpo=1)
! if iflag_t_glace=0, the old values are used:
  REAL, PARAMETER :: t_glace_min_old = 258.
  REAL, PARAMETER :: t_glace_max_old = 273.13

  REAL rel, tc, rei
  REAL k_ice0, k_ice, df
  PARAMETER (k_ice0=0.005) ! units=m2/g
  PARAMETER (df=1.66) ! diffusivity factor

  ! jq for the aerosol indirect effect
  ! jq introduced by Johannes Quaas (quaas@lmd.jussieu.fr), 27/11/2003
  ! jq
  REAL mass_solu_aero(klon, klev) ! total mass concentration for all soluble aerosols [ug m-3]
  REAL mass_solu_aero_pi(klon, klev) ! - " - (pre-industrial value)
  REAL cdnc(klon, klev) ! cloud droplet number concentration [m-3]
  REAL re(klon, klev) ! cloud droplet effective radius [um]
  REAL cdnc_pi(klon, klev) ! cloud droplet number concentration [m-3] (pi value)
  REAL re_pi(klon, klev) ! cloud droplet effective radius [um] (pi value)

  REAL fl(klon, klev) ! xliq * rneb (denominator to re; fraction of liquid water clouds
  ! within the grid cell)

  INTEGER flag_aerosol
  LOGICAL ok_cdnc
  REAL bl95_b0, bl95_b1 ! Parameter in B&L 95-Formula

  ! jq-end
  ! IM cf. CR:parametres supplementaires
  REAL zclear(klon)
  REAL zcloud(klon)
  REAL zcloudh(klon)
  REAL zcloudm(klon)
  REAL zcloudl(klon)
  REAL rhodz(klon, klev) !--rho*dz pour la couche
  REAL zrho(klon, klev) !--rho pour la couche
  REAL dh(klon, klev) !--dz pour la couche
  REAL rad_chaud(klon, klev) !--rayon pour les nuages chauds
  REAL rad_chaud_pi(klon, klev) !--rayon pour les nuages chauds pre-industriels
  REAL zflwp_var, zfiwp_var
  REAL d_rei_dt

  ! Abderrahmane oct 2009
  REAL reliq(klon, klev), reice(klon, klev)
  REAL reliq_pi(klon, klev), reice_pi(klon, klev)

  REAL,SAVE :: cdnc_min=-1.
  REAL,SAVE :: cdnc_min_m3
  !$OMP THREADPRIVATE(cdnc_min,cdnc_min_m3)
  REAL,SAVE :: cdnc_max=-1.
  REAL,SAVE :: cdnc_max_m3
  !$OMP THREADPRIVATE(cdnc_max,cdnc_max_m3)

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FH : 2011/05/24

  ! rei = ( rei_max - rei_min ) * T(°C) / 81.4 + rei_max
  ! to be used for a temperature in celcius T(°C) < 0
  ! rei=rei_min for T(°C) < -81.4

  ! Calcul de la pente de la relation entre rayon effective des cristaux
  ! et la température.
  ! Pour retrouver les résultats numériques de la version d'origine,
  ! on impose 0.71 quand on est proche de 0.71

  if (first) THEN
      call getin_p('cdnc_min',cdnc_min)
      cdnc_min_m3=cdnc_min*1.E6
      IF (cdnc_min_m3<0.) cdnc_min_m3=20.E6 ! astuce pour retrocompatibilite
      write(lunout,*)'cdnc_min=', cdnc_min_m3/1.E6
      call getin_p('cdnc_max',cdnc_max)
      cdnc_max_m3=cdnc_max*1.E6
      IF (cdnc_max_m3<0.) cdnc_max_m3=1000.E6 ! astuce pour retrocompatibilite
      write(lunout,*)'cdnc_max=', cdnc_max_m3/1.E6
  ENDIF

  d_rei_dt = (rei_max-rei_min)/81.4
  IF (abs(d_rei_dt-0.71)<1.E-4) d_rei_dt = 0.71
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculer l'epaisseur optique et l'emmissivite des nuages
  ! IM inversion des DO

  xflwp = 0.D0
  xfiwp = 0.D0
  xflwc = 0.D0
  xfiwc = 0.D0

  reliq = 0.
  reice = 0.
  reliq_pi = 0.
  reice_pi = 0.

  IF (iflag_t_glace.EQ.0) THEN
    DO k = 1, klev
      DO i = 1, klon
        ! -layer calculation
        rhodz(i, k) = (paprs(i,k)-paprs(i,k+1))/rg ! kg/m2
        zrho(i, k) = pplay(i, k)/t(i, k)/rd ! kg/m3
        dh(i, k) = rhodz(i, k)/zrho(i, k) ! m
        ! -Fraction of ice in cloud using a linear transition
        zfice(i, k) = 1.0 - (t(i,k)-t_glace_min_old)/(t_glace_max_old-t_glace_min_old)
        zfice(i, k) = min(max(zfice(i,k),0.0), 1.0)
        ! -IM Total Liquid/Ice water content
        xflwc(i, k) = (1.-zfice(i,k))*pqlwp(i, k)
        xfiwc(i, k) = zfice(i, k)*pqlwp(i, k)
      ENDDO
    ENDDO
  ELSE ! of IF (iflag_t_glace.EQ.0)
    DO k = 1, klev
        CALL icefrac_lsc(klon,t(:,k),pplay(:,k)/paprs(:,1),zfice(:,k))
 

        ! JBM: icefrac_lsc is now contained icefrac_lsc_mod
!       zfice(i, k) = icefrac_lsc(t(i,k), t_glace_min, &
!                                 t_glace_max, exposant_glace)
      DO i = 1, klon
        ! -layer calculation
        rhodz(i, k) = (paprs(i,k)-paprs(i,k+1))/rg ! kg/m2
        zrho(i, k) = pplay(i, k)/t(i, k)/rd ! kg/m3
        dh(i, k) = rhodz(i, k)/zrho(i, k) ! m
        ! -IM Total Liquid/Ice water content
        xflwc(i, k) = (1.-zfice(i,k))*pqlwp(i, k)
        xfiwc(i, k) = zfice(i, k)*pqlwp(i, k)
      ENDDO
    ENDDO
  ENDIF

  IF (ok_cdnc) THEN

    ! --we compute cloud properties as a function of the aerosol load

    DO k = 1, klev
      DO i = 1, klon
        ! Formula "D" of Boucher and Lohmann, Tellus, 1995
        ! Cloud droplet number concentration (CDNC) is restricted
        ! to be within [20, 1000 cm^3]

        ! --pre-industrial case
        cdnc_pi(i, k) = 10.**(bl95_b0+bl95_b1*log(max(mass_solu_aero_pi(i,k), &
          1.E-4))/log(10.))*1.E6 !-m-3
        cdnc_pi(i, k) = min(cdnc_max_m3, max(cdnc_min_m3,cdnc_pi(i,k)))

      ENDDO
    ENDDO

    !--flag_aerosol=7 => MACv2SP climatology 
    !--in this case there is an enhancement factor
    IF (flag_aerosol .EQ. 7) THEN 

      !--present-day
      DO k = 1, klev
        DO i = 1, klon
          cdnc(i, k) = cdnc_pi(i,k)*dNovrN(i)
        ENDDO
      ENDDO

    !--standard case
    ELSE

      DO k = 1, klev
        DO i = 1, klon

          ! Formula "D" of Boucher and Lohmann, Tellus, 1995
          ! Cloud droplet number concentration (CDNC) is restricted
          ! to be within [20, 1000 cm^3]

          ! --present-day case
          cdnc(i, k) = 10.**(bl95_b0+bl95_b1*log(max(mass_solu_aero(i,k), &
            1.E-4))/log(10.))*1.E6 !-m-3
          cdnc(i, k) = min(cdnc_max_m3, max(cdnc_min_m3,cdnc(i,k)))

        ENDDO
      ENDDO

    ENDIF !--flag_aerosol

    !--computing cloud droplet size
    DO k = 1, klev
      DO i = 1, klon

        ! --present-day case
        rad_chaud(i, k) = 1.1*((pqlwp(i,k)*pplay(i, &
          k)/(rd*t(i,k)))/(4./3*rpi*1000.*cdnc(i,k)))**(1./3.)
        rad_chaud(i, k) = max(rad_chaud(i,k)*1.E6, 5.)

        ! --pre-industrial case
        rad_chaud_pi(i, k) = 1.1*((pqlwp(i,k)*pplay(i, &
          k)/(rd*t(i,k)))/(4./3.*rpi*1000.*cdnc_pi(i,k)))**(1./3.)
        rad_chaud_pi(i, k) = max(rad_chaud_pi(i,k)*1.E6, 5.)

        ! --pre-industrial case
        ! --liquid/ice cloud water paths:
        IF (pclc(i,k)<=seuil_neb) THEN

          pcldtaupi(i, k) = 0.0

        ELSE

          zflwp_var = 1000.*(1.-zfice(i,k))*pqlwp(i, k)/pclc(i, k)* &
            rhodz(i, k)
          zfiwp_var = 1000.*zfice(i, k)*pqlwp(i, k)/pclc(i, k)*rhodz(i, k)
          tc = t(i, k) - 273.15
          rei = d_rei_dt*tc + rei_max
          IF (tc<=-81.4) rei = rei_min

          ! -- cloud optical thickness :
          ! [for liquid clouds, traditional formula,
          ! for ice clouds, Ebert & Curry (1992)]

          IF (zfiwp_var==0. .OR. rei<=0.) rei = 1.
          pcldtaupi(i, k) = 3.0/2.0*zflwp_var/rad_chaud_pi(i, k) + &
            zfiwp_var*(3.448E-03+2.431/rei)

        ENDIF

      ENDDO
    ENDDO

  ELSE !--not ok_cdnc

    ! -prescribed cloud droplet radius

    DO k = 1, min(3, klev)
      DO i = 1, klon
        rad_chaud(i, k) = rad_chau2
        rad_chaud_pi(i, k) = rad_chau2
      ENDDO
    ENDDO
    DO k = min(3, klev) + 1, klev
      DO i = 1, klon
        rad_chaud(i, k) = rad_chau1
        rad_chaud_pi(i, k) = rad_chau1
      ENDDO
    ENDDO

  ENDIF !--ok_cdnc

  ! --computation of cloud optical depth and emissivity
  ! --in the general case

  DO k = 1, klev
    DO i = 1, klon

      IF (pclc(i,k)<=seuil_neb) THEN

        ! effective cloud droplet radius (microns) for liquid water clouds:
        ! For output diagnostics cloud droplet effective radius [um]
        ! we multiply here with f * xl (fraction of liquid water
        ! clouds in the grid cell) to avoid problems in the averaging of the
        ! output.
        ! In the output of IOIPSL, derive the REAL cloud droplet
        ! effective radius as re/fl

        fl(i, k) = seuil_neb*(1.-zfice(i,k))
        re(i, k) = rad_chaud(i, k)*fl(i, k)
        rel = 0.
        rei = 0.
        pclc(i, k) = 0.0
        pcltau(i, k) = 0.0
        pclemi(i, k) = 0.0

      ELSE

        ! -- liquid/ice cloud water paths:

        zflwp_var = 1000.*(1.-zfice(i,k))*pqlwp(i, k)/pclc(i, k)*rhodz(i, k)
        zfiwp_var = 1000.*zfice(i, k)*pqlwp(i, k)/pclc(i, k)*rhodz(i, k)

        ! effective cloud droplet radius (microns) for liquid water clouds:
        ! For output diagnostics cloud droplet effective radius [um]
        ! we multiply here with f * xl (fraction of liquid water
        ! clouds in the grid cell) to avoid problems in the averaging of the
        ! output.
        ! In the output of IOIPSL, derive the REAL cloud droplet
        ! effective radius as re/fl

        fl(i, k) = pclc(i, k)*(1.-zfice(i,k))
        re(i, k) = rad_chaud(i, k)*fl(i, k)

        rel = rad_chaud(i, k)

        ! for ice clouds: as a function of the ambiant temperature
        ! [formula used by Iacobellis and Somerville (2000), with an
        ! asymptotical value of 3.5 microns at T<-81.4 C added to be
        ! consistent with observations of Heymsfield et al. 1986]:
        ! 2011/05/24 : rei_min = 3.5 becomes a free PARAMETER as well as
        ! rei_max=61.29

        tc = t(i, k) - 273.15
        rei = d_rei_dt*tc + rei_max
        IF (tc<=-81.4) rei = rei_min

        ! -- cloud optical thickness :
        ! [for liquid clouds, traditional formula,
        ! for ice clouds, Ebert & Curry (1992)]

        IF (zflwp_var==0.) rel = 1.
        IF (zfiwp_var==0. .OR. rei<=0.) rei = 1.
        pcltau(i, k) = 3.0/2.0*(zflwp_var/rel) + zfiwp_var*(3.448E-03+2.431/ &
          rei)

        ! -- cloud infrared emissivity:
        ! [the broadband infrared absorption coefficient is PARAMETERized
        ! as a function of the effective cld droplet radius]
        ! Ebert and Curry (1992) formula as used by Kiehl & Zender (1995):

        k_ice = k_ice0 + 1.0/rei

        pclemi(i, k) = 1.0 - exp(-coef_chau*zflwp_var-df*k_ice*zfiwp_var)

      ENDIF

      reice(i, k) = rei

      xflwp(i) = xflwp(i) + xflwc(i, k)*rhodz(i, k)
      xfiwp(i) = xfiwp(i) + xfiwc(i, k)*rhodz(i, k)

    ENDDO
  ENDDO

  ! --if cloud droplet radius is fixed, then pcldtaupi=pcltau

  IF (.NOT. ok_cdnc) THEN
    DO k = 1, klev
      DO i = 1, klon
        pcldtaupi(i, k) = pcltau(i, k)
        reice_pi(i, k) = reice(i, k)
      ENDDO
    ENDDO
  ENDIF

  DO k = 1, klev
    DO i = 1, klon
      reliq(i, k) = rad_chaud(i, k)
      reliq_pi(i, k) = rad_chaud_pi(i, k)
      reice_pi(i, k) = reice(i, k)
    ENDDO
  ENDDO

  ! COMPUTE CLOUD LIQUID PATH AND TOTAL CLOUDINESS
  ! IM cf. CR:test: calcul prenant ou non en compte le recouvrement
  ! initialisations

  DO i = 1, klon
    zclear(i) = 1.
    zcloud(i) = 0.
    zcloudh(i) = 0.
    zcloudm(i) = 0.
    zcloudl(i) = 0.
    pch(i) = 1.0
    pcm(i) = 1.0
    pcl(i) = 1.0
    pctlwp(i) = 0.0
  ENDDO

  ! --calculation of liquid water path

  DO k = klev, 1, -1
    DO i = 1, klon
      pctlwp(i) = pctlwp(i) + pqlwp(i, k)*rhodz(i, k)
    ENDDO
  ENDDO

  ! --calculation of cloud properties with cloud overlap

  IF (novlp==1) THEN
    DO k = klev, 1, -1
      DO i = 1, klon
        zclear(i) = zclear(i)*(1.-max(pclc(i,k),zcloud(i)))/(1.-min(real( &
          zcloud(i),kind=8),1.-zepsec))
        pct(i) = 1. - zclear(i)
        IF (paprs(i,k)<prmhc) THEN
          pch(i) = pch(i)*(1.-max(pclc(i,k),zcloudh(i)))/(1.-min(real(zcloudh &
            (i),kind=8),1.-zepsec))
          zcloudh(i) = pclc(i, k)
        ELSE IF (paprs(i,k)>=prmhc .AND. paprs(i,k)<prlmc) THEN
          pcm(i) = pcm(i)*(1.-max(pclc(i,k),zcloudm(i)))/(1.-min(real(zcloudm &
            (i),kind=8),1.-zepsec))
          zcloudm(i) = pclc(i, k)
        ELSE IF (paprs(i,k)>=prlmc) THEN
          pcl(i) = pcl(i)*(1.-max(pclc(i,k),zcloudl(i)))/(1.-min(real(zcloudl &
            (i),kind=8),1.-zepsec))
          zcloudl(i) = pclc(i, k)
        ENDIF
        zcloud(i) = pclc(i, k)
      ENDDO
    ENDDO
  ELSE IF (novlp==2) THEN
    DO k = klev, 1, -1
      DO i = 1, klon
        zcloud(i) = max(pclc(i,k), zcloud(i))
        pct(i) = zcloud(i)
        IF (paprs(i,k)<prmhc) THEN
          pch(i) = min(pclc(i,k), pch(i))
        ELSE IF (paprs(i,k)>=prmhc .AND. paprs(i,k)<prlmc) THEN
          pcm(i) = min(pclc(i,k), pcm(i))
        ELSE IF (paprs(i,k)>=prlmc) THEN
          pcl(i) = min(pclc(i,k), pcl(i))
        ENDIF
      ENDDO
    ENDDO
  ELSE IF (novlp==3) THEN
    DO k = klev, 1, -1
      DO i = 1, klon
        zclear(i) = zclear(i)*(1.-pclc(i,k))
        pct(i) = 1 - zclear(i)
        IF (paprs(i,k)<prmhc) THEN
          pch(i) = pch(i)*(1.0-pclc(i,k))
        ELSE IF (paprs(i,k)>=prmhc .AND. paprs(i,k)<prlmc) THEN
          pcm(i) = pcm(i)*(1.0-pclc(i,k))
        ELSE IF (paprs(i,k)>=prlmc) THEN
          pcl(i) = pcl(i)*(1.0-pclc(i,k))
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  DO i = 1, klon
    pch(i) = 1. - pch(i)
    pcm(i) = 1. - pcm(i)
    pcl(i) = 1. - pcl(i)
  ENDDO

  ! ========================================================
  ! DIAGNOSTICS CALCULATION FOR CMIP5 PROTOCOL
  ! ========================================================
  ! change by Nicolas Yan (LSCE)
  ! Cloud Droplet Number Concentration (CDNC) : 3D variable
  ! Fractionnal cover by liquid water cloud (LCC3D) : 3D variable
  ! Cloud Droplet Number Concentration at top of cloud (CLDNCL) : 2D variable
  ! Droplet effective radius at top of cloud (REFFCLWTOP) : 2D variable
  ! Fractionnal cover by liquid water at top of clouds (LCC) : 2D variable

  IF (ok_cdnc) THEN

    DO k = 1, klev
      DO i = 1, klon
        phase3d(i, k) = 1 - zfice(i, k)
        IF (pclc(i,k)<=seuil_neb) THEN
          lcc3d(i, k) = seuil_neb*phase3d(i, k)
        ELSE
          lcc3d(i, k) = pclc(i, k)*phase3d(i, k)
        ENDIF
        scdnc(i, k) = lcc3d(i, k)*cdnc(i, k) ! m-3
      ENDDO
    ENDDO

    DO i = 1, klon
      lcc(i) = 0.
      reffclwtop(i) = 0.
      cldncl(i) = 0.
      IF (novlp.EQ.3 .OR. novlp.EQ.1) tcc(i) = 1.
      IF (novlp.EQ.2) tcc(i) = 0.
    ENDDO

    DO i = 1, klon
      DO k = klev - 1, 1, -1 !From TOA down

          ! Test, if the cloud optical depth exceeds the necessary
          ! threshold:

        IF (pcltau(i,k)>thres_tau .AND. pclc(i,k)>thres_neb) THEN

          IF (novlp.EQ.2) THEN
            IF (first) THEN
              WRITE (*, *) 'Hypothese de recouvrement: MAXIMUM'
              first = .FALSE.
            ENDIF
            flag_max = -1.
            ftmp(i) = max(tcc(i), pclc(i,k))
          ENDIF

          IF (novlp.EQ.3) THEN
            IF (first) THEN
              WRITE (*, *) 'Hypothese de recouvrement: RANDOM'
              first = .FALSE.
            ENDIF
            flag_max = 1.
            ftmp(i) = tcc(i)*(1-pclc(i,k))
          ENDIF

          IF (novlp.EQ.1) THEN
            IF (first) THEN
              WRITE (*, *) 'Hypothese de recouvrement: MAXIMUM_ &
                &                                             &
                &                                          RANDOM'
              first = .FALSE.
            ENDIF
            flag_max = 1.
            ftmp(i) = tcc(i)*(1.-max(pclc(i,k),pclc(i,k+1)))/(1.-min(pclc(i, &
              k+1),1.-thres_neb))
          ENDIF
          ! Effective radius of cloud droplet at top of cloud (m)
          reffclwtop(i) = reffclwtop(i) + rad_chaud(i, k)*1.0E-06*phase3d(i, &
            k)*(tcc(i)-ftmp(i))*flag_max
          ! CDNC at top of cloud (m-3)
          cldncl(i) = cldncl(i) + cdnc(i, k)*phase3d(i, k)*(tcc(i)-ftmp(i))* &
            flag_max
          ! Liquid Cloud Content at top of cloud
          lcc(i) = lcc(i) + phase3d(i, k)*(tcc(i)-ftmp(i))*flag_max
          ! Total Cloud Content at top of cloud
          tcc(i) = ftmp(i)

        ENDIF ! is there a visible, not-too-small cloud?
      ENDDO ! loop over k

      IF (novlp.EQ.3 .OR. novlp.EQ.1) tcc(i) = 1. - tcc(i)

    ENDDO ! loop over i

    ! ! Convective and Stratiform Cloud Droplet Effective Radius (REFFCLWC
    ! REFFCLWS)
    DO i = 1, klon
      DO k = 1, klev
        ! Weight to be used for outputs: eau_liquide*couverture nuageuse
        lcc3dcon(i, k) = rnebcon(i, k)*phase3d(i, k)*clwcon(i, k) ! eau liquide convective
        lcc3dstra(i, k) = pclc(i, k)*pqlwp(i, k)*phase3d(i, k)
        lcc3dstra(i, k) = lcc3dstra(i, k) - lcc3dcon(i, k) ! eau liquide stratiforme
        lcc3dstra(i, k) = max(lcc3dstra(i,k), 0.0)
        !FC pour la glace (CAUSES)
        icc3dcon(i, k) = rnebcon(i, k)*(1-phase3d(i, k))*clwcon(i, k) !  glace convective
        icc3dstra(i, k)= pclc(i, k)*pqlwp(i, k)*(1-phase3d(i, k))
        icc3dstra(i, k) = icc3dstra(i, k) - icc3dcon(i, k) ! glace stratiforme
        icc3dstra(i, k) = max( icc3dstra(i, k), 0.0)
        !FC (CAUSES)

        ! Compute cloud droplet radius as above in meter
        radius = 1.1*((pqlwp(i,k)*pplay(i,k)/(rd*t(i,k)))/(4./3*rpi*1000.* &
          cdnc(i,k)))**(1./3.)
        radius = max(radius, 5.E-6)
        ! Convective Cloud Droplet Effective Radius (REFFCLWC) : variable 3D
        reffclwc(i, k) = radius
        reffclwc(i, k) = reffclwc(i, k)*lcc3dcon(i, k)
        ! Stratiform Cloud Droplet Effective Radius (REFFCLWS) : variable 3D
        reffclws(i, k) = radius
        reffclws(i, k) = reffclws(i, k)*lcc3dstra(i, k)
      ENDDO !klev
    ENDDO !klon

    ! Column Integrated Cloud Droplet Number (CLDNVI) : variable 2D

    DO i = 1, klon
      cldnvi(i) = 0.
      lcc_integrat(i) = 0.
      height(i) = 0.
      DO k = 1, klev
        cldnvi(i) = cldnvi(i) + cdnc(i, k)*lcc3d(i, k)*dh(i, k)
        lcc_integrat(i) = lcc_integrat(i) + lcc3d(i, k)*dh(i, k)
        height(i) = height(i) + dh(i, k)
      ENDDO ! klev
      lcc_integrat(i) = lcc_integrat(i)/height(i)
      IF (lcc_integrat(i)<=1.0E-03) THEN
        cldnvi(i) = cldnvi(i)*lcc(i)/seuil_neb
      ELSE
        cldnvi(i) = cldnvi(i)*lcc(i)/lcc_integrat(i)
      ENDIF
    ENDDO ! klon

    DO i = 1, klon
      DO k = 1, klev
        IF (scdnc(i,k)<=0.0) scdnc(i, k) = 0.0
        IF (reffclws(i,k)<=0.0) reffclws(i, k) = 0.0
        IF (reffclwc(i,k)<=0.0) reffclwc(i, k) = 0.0
        IF (lcc3d(i,k)<=0.0) lcc3d(i, k) = 0.0
        IF (lcc3dcon(i,k)<=0.0) lcc3dcon(i, k) = 0.0
        IF (lcc3dstra(i,k)<=0.0) lcc3dstra(i, k) = 0.0
!FC (CAUSES)
        IF (icc3dcon(i,k)<=0.0) icc3dcon(i, k) = 0.0
        IF (icc3dstra(i,k)<=0.0) icc3dstra(i, k) = 0.0
!FC (CAUSES)
      ENDDO
      IF (reffclwtop(i)<=0.0) reffclwtop(i) = 0.0
      IF (cldncl(i)<=0.0) cldncl(i) = 0.0
      IF (cldnvi(i)<=0.0) cldnvi(i) = 0.0
      IF (lcc(i)<=0.0) lcc(i) = 0.0
    ENDDO

  ENDIF !ok_cdnc

  first=.false. !to be sure

  RETURN

END SUBROUTINE newmicro
