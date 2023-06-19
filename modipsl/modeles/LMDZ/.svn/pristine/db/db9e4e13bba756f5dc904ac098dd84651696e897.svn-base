! IM ctes ds clesphys.h   SUBROUTINE SW(PSCT, RCO2, PRMU0, PFRAC,
SUBROUTINE sw_lmdar4(psct, prmu0, pfrac, ppmb, pdp, ppsol, palbd, palbp, &
    ptave, pwv, pqs, pozon, paer, pcldsw, ptau, pomega, pcg, pheat, pheat0, &
    palbpla, ptopsw, psolsw, ptopsw0, psolsw0, zfsup, zfsdn, zfsup0, zfsdn0, &
    tauae, pizae, cgae, ptaua, pomegaa, ptopswad, psolswad, ptopswai, &
    psolswai, ok_ade, ok_aie)
  USE dimphy
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE

  include "YOMCST.h"

  ! ------------------------------------------------------------------

  ! PURPOSE.
  ! --------

  ! THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
  ! SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

  ! METHOD.
  ! -------

  ! 1. COMPUTES ABSORBER AMOUNTS                 (SWU)
  ! 2. COMPUTES FLUXES IN 1ST SPECTRAL INTERVAL  (SW1S)
  ! 3. COMPUTES FLUXES IN 2ND SPECTRAL INTERVAL  (SW2S)

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
  ! DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! 95-01-01   J.-J. MORCRETTE  Direct/Diffuse Albedo
  ! 03-11-27   J. QUAAS Introduce aerosol forcings (based on BOUCHER)
  ! ------------------------------------------------------------------

  ! * ARGUMENTS:

  REAL (KIND=8) psct ! constante solaire (valeur conseillee: 1370)
  ! IM ctes ds clesphys.h   REAL(KIND=8) RCO2  ! concentration CO2 (IPCC:
  ! 353.E-06*44.011/28.97)
  include "clesphys.h"

  REAL (KIND=8) ppsol(kdlon) ! SURFACE PRESSURE (PA)
  REAL (KIND=8) pdp(kdlon, kflev) ! LAYER THICKNESS (PA)
  REAL (KIND=8) ppmb(kdlon, kflev+1) ! HALF-LEVEL PRESSURE (MB)

  REAL (KIND=8) prmu0(kdlon) ! COSINE OF ZENITHAL ANGLE
  REAL (KIND=8) pfrac(kdlon) ! fraction de la journee

  REAL (KIND=8) ptave(kdlon, kflev) ! LAYER TEMPERATURE (K)
  REAL (KIND=8) pwv(kdlon, kflev) ! SPECIFIC HUMIDITY (KG/KG)
  REAL (KIND=8) pqs(kdlon, kflev) ! SATURATED WATER VAPOUR (KG/KG)
  REAL (KIND=8) pozon(kdlon, kflev) ! OZONE CONCENTRATION (KG/KG)
  REAL (KIND=8) paer(kdlon, kflev, 5) ! AEROSOLS' OPTICAL THICKNESS

  REAL (KIND=8) palbd(kdlon, 2) ! albedo du sol (lumiere diffuse)
  REAL (KIND=8) palbp(kdlon, 2) ! albedo du sol (lumiere parallele)

  REAL (KIND=8) pcldsw(kdlon, kflev) ! CLOUD FRACTION
  REAL (KIND=8) ptau(kdlon, 2, kflev) ! CLOUD OPTICAL THICKNESS
  REAL (KIND=8) pcg(kdlon, 2, kflev) ! ASYMETRY FACTOR
  REAL (KIND=8) pomega(kdlon, 2, kflev) ! SINGLE SCATTERING ALBEDO

  REAL (KIND=8) pheat(kdlon, kflev) ! SHORTWAVE HEATING (K/DAY)
  REAL (KIND=8) pheat0(kdlon, kflev) ! SHORTWAVE HEATING (K/DAY) clear-sky
  REAL (KIND=8) palbpla(kdlon) ! PLANETARY ALBEDO
  REAL (KIND=8) ptopsw(kdlon) ! SHORTWAVE FLUX AT T.O.A.
  REAL (KIND=8) psolsw(kdlon) ! SHORTWAVE FLUX AT SURFACE
  REAL (KIND=8) ptopsw0(kdlon) ! SHORTWAVE FLUX AT T.O.A. (CLEAR-SKY)
  REAL (KIND=8) psolsw0(kdlon) ! SHORTWAVE FLUX AT SURFACE (CLEAR-SKY)

  ! * LOCAL VARIABLES:

  REAL, PARAMETER :: dobson_u = 2.1415E-05 ! Dobson unit, in kg m-2

  REAL (KIND=8) zoz(kdlon, kflev)
  ! column-density of ozone in layer, in kilo-Dobsons

  REAL (KIND=8) zaki(kdlon, 2)
  REAL (KIND=8) zcld(kdlon, kflev)
  REAL (KIND=8) zclear(kdlon)
  REAL (KIND=8) zdsig(kdlon, kflev)
  REAL (KIND=8) zfact(kdlon)
  REAL (KIND=8) zfd(kdlon, kflev+1)
  REAL (KIND=8) zfdown(kdlon, kflev+1)
  REAL (KIND=8) zfu(kdlon, kflev+1)
  REAL (KIND=8) zfup(kdlon, kflev+1)
  REAL (KIND=8) zrmu(kdlon)
  REAL (KIND=8) zsec(kdlon)
  REAL (KIND=8) zud(kdlon, 5, kflev+1)
  REAL (KIND=8) zcldsw0(kdlon, kflev)

  REAL (KIND=8) zfsup(kdlon, kflev+1)
  REAL (KIND=8) zfsdn(kdlon, kflev+1)
  REAL (KIND=8) zfsup0(kdlon, kflev+1)
  REAL (KIND=8) zfsdn0(kdlon, kflev+1)

  INTEGER inu, jl, jk, i, k, kpl1

  INTEGER swpas ! Every swpas steps, sw is calculated
  PARAMETER (swpas=1)

  INTEGER itapsw
  LOGICAL appel1er
  DATA itapsw/0/
  DATA appel1er/.TRUE./
  SAVE itapsw, appel1er
  !$OMP THREADPRIVATE(appel1er)
  !$OMP THREADPRIVATE(itapsw)
  ! jq-Introduced for aerosol forcings
  REAL (KIND=8) flag_aer
  LOGICAL ok_ade, ok_aie ! use aerosol forcings or not?
  REAL (KIND=8) tauae(kdlon, kflev, 2) ! aerosol optical properties
  REAL (KIND=8) pizae(kdlon, kflev, 2) ! (see aeropt.F)
  REAL (KIND=8) cgae(kdlon, kflev, 2) ! -"-
  REAL (KIND=8) ptaua(kdlon, 2, kflev) ! CLOUD OPTICAL THICKNESS (pre-industrial value)
  REAL (KIND=8) pomegaa(kdlon, 2, kflev) ! SINGLE SCATTERING ALBEDO
  REAL (KIND=8) ptopswad(kdlon) ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL DIR)
  REAL (KIND=8) psolswad(kdlon) ! SHORTWAVE FLUX AT SURFACE(+AEROSOL DIR)
  REAL (KIND=8) ptopswai(kdlon) ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL IND)
  REAL (KIND=8) psolswai(kdlon) ! SHORTWAVE FLUX AT SURFACE(+AEROSOL IND)
  ! jq - Fluxes including aerosol effects
  REAL (KIND=8), ALLOCATABLE, SAVE :: zfsupad(:, :)
  !$OMP THREADPRIVATE(ZFSUPAD)
  REAL (KIND=8), ALLOCATABLE, SAVE :: zfsdnad(:, :)
  !$OMP THREADPRIVATE(ZFSDNAD)
  REAL (KIND=8), ALLOCATABLE, SAVE :: zfsupai(:, :)
  !$OMP THREADPRIVATE(ZFSUPAI)
  REAL (KIND=8), ALLOCATABLE, SAVE :: zfsdnai(:, :)
  !$OMP THREADPRIVATE(ZFSDNAI)
  LOGICAL initialized
  ! ym      SAVE ZFSUPAD, ZFSDNAD, ZFSUPAI, ZFSDNAI ! aerosol fluxes
  ! rv
  SAVE flag_aer
  !$OMP THREADPRIVATE(flag_aer)
  DATA initialized/.FALSE./
  SAVE initialized
  !$OMP THREADPRIVATE(initialized)
  ! jq-end
  REAL tmp_

  IF (.NOT. initialized) THEN
    flag_aer = 0.
    initialized = .TRUE.
    ALLOCATE (zfsupad(kdlon,kflev+1))
    ALLOCATE (zfsdnad(kdlon,kflev+1))
    ALLOCATE (zfsupai(kdlon,kflev+1))
    ALLOCATE (zfsdnai(kdlon,kflev+1))

    zfsupad(:, :) = 0.
    zfsdnad(:, :) = 0.
    zfsupai(:, :) = 0.
    zfsdnai(:, :) = 0.
  END IF

  IF (appel1er) THEN
    WRITE (lunout, *) 'SW calling frequency : ', swpas
    WRITE (lunout, *) '   In general, it should be 1'
    appel1er = .FALSE.
  END IF
  ! ------------------------------------------------------------------
  IF (mod(itapsw,swpas)==0) THEN

    tmp_ = 1./(dobson_u*1E3*rg)
    ! cdir collapse
    DO jk = 1, kflev
      DO jl = 1, kdlon
        zcldsw0(jl, jk) = 0.0
        zoz(jl, jk) = pozon(jl, jk)*tmp_*pdp(jl, jk)
      END DO
    END DO


    ! clear-sky:
    ! IM ctes ds clesphys.h  CALL SWU(PSCT,RCO2,ZCLDSW0,PPMB,PPSOL,
    CALL swu_lmdar4(psct, zcldsw0, ppmb, ppsol, prmu0, pfrac, ptave, pwv, &
      zaki, zcld, zclear, zdsig, zfact, zrmu, zsec, zud)
    inu = 1
    CALL sw1s_lmdar4(inu, paer, flag_aer, tauae, pizae, cgae, palbd, palbp, &
      pcg, zcld, zclear, zcldsw0, zdsig, pomega, zoz, zrmu, zsec, ptau, zud, &
      zfd, zfu)
    inu = 2
    CALL sw2s_lmdar4(inu, paer, flag_aer, tauae, pizae, cgae, zaki, palbd, &
      palbp, pcg, zcld, zclear, zcldsw0, zdsig, pomega, zoz, zrmu, zsec, &
      ptau, zud, pwv, pqs, zfdown, zfup)
    DO jk = 1, kflev + 1
      DO jl = 1, kdlon
        zfsup0(jl, jk) = (zfup(jl,jk)+zfu(jl,jk))*zfact(jl)
        zfsdn0(jl, jk) = (zfdown(jl,jk)+zfd(jl,jk))*zfact(jl)
      END DO
    END DO

    flag_aer = 0.0
    CALL swu_lmdar4(psct, pcldsw, ppmb, ppsol, prmu0, pfrac, ptave, pwv, &
      zaki, zcld, zclear, zdsig, zfact, zrmu, zsec, zud)
    inu = 1
    CALL sw1s_lmdar4(inu, paer, flag_aer, tauae, pizae, cgae, palbd, palbp, &
      pcg, zcld, zclear, pcldsw, zdsig, pomega, zoz, zrmu, zsec, ptau, zud, &
      zfd, zfu)
    inu = 2
    CALL sw2s_lmdar4(inu, paer, flag_aer, tauae, pizae, cgae, zaki, palbd, &
      palbp, pcg, zcld, zclear, pcldsw, zdsig, pomega, zoz, zrmu, zsec, ptau, &
      zud, pwv, pqs, zfdown, zfup)

    ! cloudy-sky:

    DO jk = 1, kflev + 1
      DO jl = 1, kdlon
        zfsup(jl, jk) = (zfup(jl,jk)+zfu(jl,jk))*zfact(jl)
        zfsdn(jl, jk) = (zfdown(jl,jk)+zfd(jl,jk))*zfact(jl)
      END DO
    END DO


    IF (ok_ade) THEN

      ! cloudy-sky + aerosol dir OB
      flag_aer = 1.0
      CALL swu_lmdar4(psct, pcldsw, ppmb, ppsol, prmu0, pfrac, ptave, pwv, &
        zaki, zcld, zclear, zdsig, zfact, zrmu, zsec, zud)
      inu = 1
      CALL sw1s_lmdar4(inu, paer, flag_aer, tauae, pizae, cgae, palbd, palbp, &
        pcg, zcld, zclear, pcldsw, zdsig, pomega, zoz, zrmu, zsec, ptau, zud, &
        zfd, zfu)
      inu = 2
      CALL sw2s_lmdar4(inu, paer, flag_aer, tauae, pizae, cgae, zaki, palbd, &
        palbp, pcg, zcld, zclear, pcldsw, zdsig, pomega, zoz, zrmu, zsec, &
        ptau, zud, pwv, pqs, zfdown, zfup)
      DO jk = 1, kflev + 1
        DO jl = 1, kdlon
          zfsupad(jl, jk) = zfsup(jl, jk)
          zfsdnad(jl, jk) = zfsdn(jl, jk)
          zfsup(jl, jk) = (zfup(jl,jk)+zfu(jl,jk))*zfact(jl)
          zfsdn(jl, jk) = (zfdown(jl,jk)+zfd(jl,jk))*zfact(jl)
        END DO
      END DO

    END IF ! ok_ade

    IF (ok_aie) THEN

      ! jq   cloudy-sky + aerosol direct + aerosol indirect
      flag_aer = 1.0
      CALL swu_lmdar4(psct, pcldsw, ppmb, ppsol, prmu0, pfrac, ptave, pwv, &
        zaki, zcld, zclear, zdsig, zfact, zrmu, zsec, zud)
      inu = 1
      CALL sw1s_lmdar4(inu, paer, flag_aer, tauae, pizae, cgae, palbd, palbp, &
        pcg, zcld, zclear, pcldsw, zdsig, pomegaa, zoz, zrmu, zsec, ptaua, &
        zud, zfd, zfu)
      inu = 2
      CALL sw2s_lmdar4(inu, paer, flag_aer, tauae, pizae, cgae, zaki, palbd, &
        palbp, pcg, zcld, zclear, pcldsw, zdsig, pomegaa, zoz, zrmu, zsec, &
        ptaua, zud, pwv, pqs, zfdown, zfup)
      DO jk = 1, kflev + 1
        DO jl = 1, kdlon
          zfsupai(jl, jk) = zfsup(jl, jk)
          zfsdnai(jl, jk) = zfsdn(jl, jk)
          zfsup(jl, jk) = (zfup(jl,jk)+zfu(jl,jk))*zfact(jl)
          zfsdn(jl, jk) = (zfdown(jl,jk)+zfd(jl,jk))*zfact(jl)
        END DO
      END DO
    END IF ! ok_aie
    ! jq -end

    itapsw = 0
  END IF
  itapsw = itapsw + 1

  DO k = 1, kflev
    kpl1 = k + 1
    DO i = 1, kdlon
      pheat(i, k) = -(zfsup(i,kpl1)-zfsup(i,k)) - (zfsdn(i,k)-zfsdn(i,kpl1))
      pheat(i, k) = pheat(i, k)*rday*rg/rcpd/pdp(i, k)
      pheat0(i, k) = -(zfsup0(i,kpl1)-zfsup0(i,k)) - &
        (zfsdn0(i,k)-zfsdn0(i,kpl1))
      pheat0(i, k) = pheat0(i, k)*rday*rg/rcpd/pdp(i, k)
    END DO
  END DO
  DO i = 1, kdlon
    palbpla(i) = zfsup(i, kflev+1)/(zfsdn(i,kflev+1)+1.0E-20)

    psolsw(i) = zfsdn(i, 1) - zfsup(i, 1)
    ptopsw(i) = zfsdn(i, kflev+1) - zfsup(i, kflev+1)

    psolsw0(i) = zfsdn0(i, 1) - zfsup0(i, 1)
    ptopsw0(i) = zfsdn0(i, kflev+1) - zfsup0(i, kflev+1)
    ! -OB
    psolswad(i) = zfsdnad(i, 1) - zfsupad(i, 1)
    ptopswad(i) = zfsdnad(i, kflev+1) - zfsupad(i, kflev+1)

    psolswai(i) = zfsdnai(i, 1) - zfsupai(i, 1)
    ptopswai(i) = zfsdnai(i, kflev+1) - zfsupai(i, kflev+1)
    ! -fin
  END DO

  RETURN
END SUBROUTINE sw_lmdar4

! IM ctes ds clesphys.h   SUBROUTINE SWU
! (PSCT,RCO2,PCLDSW,PPMB,PPSOL,PRMU0,PFRAC,
SUBROUTINE swu_lmdar4(psct, pcldsw, ppmb, ppsol, prmu0, pfrac, ptave, pwv, &
    paki, pcld, pclear, pdsig, pfact, prmu, psec, pud)
  USE dimphy
  USE radiation_ar4_param, ONLY: zpdh2o, zpdumg, zprh2o, zprumg, rtdh2o, &
    rtdumg, rth2o, rtumg
  IMPLICIT NONE
  include "radepsi.h"
  include "radopt.h"
  include "YOMCST.h"

  ! * ARGUMENTS:

  REAL (KIND=8) psct
  ! IM ctes ds clesphys.h   REAL(KIND=8) RCO2
  include "clesphys.h"
  REAL (KIND=8) pcldsw(kdlon, kflev)
  REAL (KIND=8) ppmb(kdlon, kflev+1)
  REAL (KIND=8) ppsol(kdlon)
  REAL (KIND=8) prmu0(kdlon)
  REAL (KIND=8) pfrac(kdlon)
  REAL (KIND=8) ptave(kdlon, kflev)
  REAL (KIND=8) pwv(kdlon, kflev)

  REAL (KIND=8) paki(kdlon, 2)
  REAL (KIND=8) pcld(kdlon, kflev)
  REAL (KIND=8) pclear(kdlon)
  REAL (KIND=8) pdsig(kdlon, kflev)
  REAL (KIND=8) pfact(kdlon)
  REAL (KIND=8) prmu(kdlon)
  REAL (KIND=8) psec(kdlon)
  REAL (KIND=8) pud(kdlon, 5, kflev+1)

  ! * LOCAL VARIABLES:

  INTEGER iind(2)
  REAL (KIND=8) zc1j(kdlon, kflev+1)
  REAL (KIND=8) zclear(kdlon)
  REAL (KIND=8) zcloud(kdlon)
  REAL (KIND=8) zn175(kdlon)
  REAL (KIND=8) zn190(kdlon)
  REAL (KIND=8) zo175(kdlon)
  REAL (KIND=8) zo190(kdlon)
  REAL (KIND=8) zsign(kdlon)
  REAL (KIND=8) zr(kdlon, 2)
  REAL (KIND=8) zsigo(kdlon)
  REAL (KIND=8) zud(kdlon, 2)
  REAL (KIND=8) zrth, zrtu, zwh2o, zdsco2, zdsh2o, zfppw
  INTEGER jl, jk, jkp1, jkl, jklp1, ja

  ! ------------------------------------------------------------------

  ! *         1.     COMPUTES AMOUNTS OF ABSORBERS
  ! -----------------------------


  iind(1) = 1
  iind(2) = 2

  ! *         1.1    INITIALIZES QUANTITIES
  ! ----------------------


  DO jl = 1, kdlon
    pud(jl, 1, kflev+1) = 0.
    pud(jl, 2, kflev+1) = 0.
    pud(jl, 3, kflev+1) = 0.
    pud(jl, 4, kflev+1) = 0.
    pud(jl, 5, kflev+1) = 0.
    pfact(jl) = prmu0(jl)*pfrac(jl)*psct
    prmu(jl) = sqrt(1224.*prmu0(jl)*prmu0(jl)+1.)/35.
    psec(jl) = 1./prmu(jl)
    zc1j(jl, kflev+1) = 0.
  END DO

  ! *          1.3    AMOUNTS OF ABSORBERS
  ! --------------------


  DO jl = 1, kdlon
    zud(jl, 1) = 0.
    zud(jl, 2) = 0.
    zo175(jl) = ppsol(jl)**(zpdumg+1.)
    zo190(jl) = ppsol(jl)**(zpdh2o+1.)
    zsigo(jl) = ppsol(jl)
    zclear(jl) = 1.
    zcloud(jl) = 0.
  END DO

  DO jk = 1, kflev
    jkp1 = jk + 1
    jkl = kflev + 1 - jk
    jklp1 = jkl + 1
    DO jl = 1, kdlon
      zrth = (rth2o/ptave(jl,jk))**rtdh2o
      zrtu = (rtumg/ptave(jl,jk))**rtdumg
      zwh2o = max(pwv(jl,jk), zepscq)
      zsign(jl) = 100.*ppmb(jl, jkp1)
      pdsig(jl, jk) = (zsigo(jl)-zsign(jl))/ppsol(jl)
      zn175(jl) = zsign(jl)**(zpdumg+1.)
      zn190(jl) = zsign(jl)**(zpdh2o+1.)
      zdsco2 = zo175(jl) - zn175(jl)
      zdsh2o = zo190(jl) - zn190(jl)
      pud(jl, 1, jk) = 1./(10.*rg*(zpdh2o+1.))/(zprh2o**zpdh2o)*zdsh2o*zwh2o* &
        zrth
      pud(jl, 2, jk) = 1./(10.*rg*(zpdumg+1.))/(zprumg**zpdumg)*zdsco2*rco2* &
        zrtu
      zfppw = 1.6078*zwh2o/(1.+0.608*zwh2o)
      pud(jl, 4, jk) = pud(jl, 1, jk)*zfppw
      pud(jl, 5, jk) = pud(jl, 1, jk)*(1.-zfppw)
      zud(jl, 1) = zud(jl, 1) + pud(jl, 1, jk)
      zud(jl, 2) = zud(jl, 2) + pud(jl, 2, jk)
      zsigo(jl) = zsign(jl)
      zo175(jl) = zn175(jl)
      zo190(jl) = zn190(jl)

      IF (novlp==1) THEN
        zclear(jl) = zclear(jl)*(1.-max(pcldsw(jl,jkl),zcloud(jl)))/(1.-min( &
          zcloud(jl),1.-zepsec))
        zc1j(jl, jkl) = 1.0 - zclear(jl)
        zcloud(jl) = pcldsw(jl, jkl)
      ELSE IF (novlp==2) THEN
        zcloud(jl) = max(pcldsw(jl,jkl), zcloud(jl))
        zc1j(jl, jkl) = zcloud(jl)
      ELSE IF (novlp==3) THEN
        zclear(jl) = zclear(jl)*(1.-pcldsw(jl,jkl))
        zcloud(jl) = 1.0 - zclear(jl)
        zc1j(jl, jkl) = zcloud(jl)
      END IF
    END DO
  END DO
  DO jl = 1, kdlon
    pclear(jl) = 1. - zc1j(jl, 1)
  END DO
  DO jk = 1, kflev
    DO jl = 1, kdlon
      IF (pclear(jl)<1.) THEN
        pcld(jl, jk) = pcldsw(jl, jk)/(1.-pclear(jl))
      ELSE
        pcld(jl, jk) = 0.
      END IF
    END DO
  END DO

  ! *         1.4    COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS
  ! -----------------------------------------------


  DO ja = 1, 2
    DO jl = 1, kdlon
      zud(jl, ja) = zud(jl, ja)*psec(jl)
    END DO
  END DO

  CALL swtt1_lmdar4(2, 2, iind, zud, zr)

  DO ja = 1, 2
    DO jl = 1, kdlon
      paki(jl, ja) = -log(zr(jl,ja))/zud(jl, ja)
    END DO
  END DO


  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE swu_lmdar4
SUBROUTINE sw1s_lmdar4(knu, paer, flag_aer, tauae, pizae, cgae, palbd, palbp, &
    pcg, pcld, pclear, pcldsw, pdsig, pomega, poz, prmu, psec, ptau, pud, &
    pfd, pfu)
  USE dimphy
  USE radiation_ar4_param, ONLY: rsun, rray
  USE infotrac_phy, ONLY: type_trac
#ifdef REPROBUS
  USE chem_rep, ONLY: rsuntime, ok_suntime
  USE print_control_mod, ONLY: lunout
#endif

  IMPLICIT NONE

  ! ------------------------------------------------------------------
  ! PURPOSE.
  ! --------

  ! THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
  ! SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

  ! METHOD.
  ! -------

  ! 1. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
  ! CONTINUUM SCATTERING
  ! 2. MULTIPLY BY OZONE TRANSMISSION FUNCTION

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
  ! DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! 94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
  ! ------------------------------------------------------------------

  ! * ARGUMENTS:

  INTEGER knu
  ! -OB
  REAL (KIND=8) flag_aer
  REAL (KIND=8) tauae(kdlon, kflev, 2)
  REAL (KIND=8) pizae(kdlon, kflev, 2)
  REAL (KIND=8) cgae(kdlon, kflev, 2)
  REAL (KIND=8) paer(kdlon, kflev, 5)
  REAL (KIND=8) palbd(kdlon, 2)
  REAL (KIND=8) palbp(kdlon, 2)
  REAL (KIND=8) pcg(kdlon, 2, kflev)
  REAL (KIND=8) pcld(kdlon, kflev)
  REAL (KIND=8) pcldsw(kdlon, kflev)
  REAL (KIND=8) pclear(kdlon)
  REAL (KIND=8) pdsig(kdlon, kflev)
  REAL (KIND=8) pomega(kdlon, 2, kflev)
  REAL (KIND=8) poz(kdlon, kflev)
  REAL (KIND=8) prmu(kdlon)
  REAL (KIND=8) psec(kdlon)
  REAL (KIND=8) ptau(kdlon, 2, kflev)
  REAL (KIND=8) pud(kdlon, 5, kflev+1)

  REAL (KIND=8) pfd(kdlon, kflev+1)
  REAL (KIND=8) pfu(kdlon, kflev+1)

  ! * LOCAL VARIABLES:

  INTEGER iind(4)

  REAL (KIND=8) zcgaz(kdlon, kflev)
  REAL (KIND=8) zdiff(kdlon)
  REAL (KIND=8) zdirf(kdlon)
  REAL (KIND=8) zpizaz(kdlon, kflev)
  REAL (KIND=8) zrayl(kdlon)
  REAL (KIND=8) zray1(kdlon, kflev+1)
  REAL (KIND=8) zray2(kdlon, kflev+1)
  REAL (KIND=8) zrefz(kdlon, 2, kflev+1)
  REAL (KIND=8) zrj(kdlon, 6, kflev+1)
  REAL (KIND=8) zrj0(kdlon, 6, kflev+1)
  REAL (KIND=8) zrk(kdlon, 6, kflev+1)
  REAL (KIND=8) zrk0(kdlon, 6, kflev+1)
  REAL (KIND=8) zrmue(kdlon, kflev+1)
  REAL (KIND=8) zrmu0(kdlon, kflev+1)
  REAL (KIND=8) zr(kdlon, 4)
  REAL (KIND=8) ztauaz(kdlon, kflev)
  REAL (KIND=8) ztra1(kdlon, kflev+1)
  REAL (KIND=8) ztra2(kdlon, kflev+1)
  REAL (KIND=8) zw(kdlon, 4)

  INTEGER jl, jk, k, jaj, ikm1, ikl

  ! If running with Reporbus, overwrite default values of RSUN.
  ! Otherwise keep default values from radiation_AR4_param module.
  IF (type_trac=='repr') THEN
#ifdef REPROBUS
    IF (ok_suntime) THEN
      rsun(1) = rsuntime(1)
      rsun(2) = rsuntime(2)
    END IF
    WRITE (lunout, *) 'RSUN(1): ', rsun(1)
#endif
  END IF

  ! ------------------------------------------------------------------

  ! *         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
  ! ----------------------- ------------------



  ! *         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
  ! -----------------------------------------


  DO jl = 1, kdlon
    zrayl(jl) = rray(knu, 1) + prmu(jl)*(rray(knu,2)+prmu(jl)*(rray(knu, &
      3)+prmu(jl)*(rray(knu,4)+prmu(jl)*(rray(knu,5)+prmu(jl)*rray(knu,6)))))
  END DO


  ! ------------------------------------------------------------------

  ! *         2.    CONTINUUM SCATTERING CALCULATIONS
  ! ---------------------------------


  ! *         2.1   CLEAR-SKY FRACTION OF THE COLUMN
  ! --------------------------------


  CALL swclr_lmdar4(knu, paer, flag_aer, tauae, pizae, cgae, palbp, pdsig, &
    zrayl, psec, zcgaz, zpizaz, zray1, zray2, zrefz, zrj0, zrk0, zrmu0, &
    ztauaz, ztra1, ztra2)

  ! *         2.2   CLOUDY FRACTION OF THE COLUMN
  ! -----------------------------


  CALL swr_lmdar4(knu, palbd, pcg, pcld, pdsig, pomega, zrayl, psec, ptau, &
    zcgaz, zpizaz, zray1, zray2, zrefz, zrj, zrk, zrmue, ztauaz, ztra1, &
    ztra2)

  ! ------------------------------------------------------------------

  ! *         3.    OZONE ABSORPTION
  ! ----------------


  iind(1) = 1
  iind(2) = 3
  iind(3) = 1
  iind(4) = 3

  ! *         3.1   DOWNWARD FLUXES
  ! ---------------


  jaj = 2

  DO jl = 1, kdlon
    zw(jl, 1) = 0.
    zw(jl, 2) = 0.
    zw(jl, 3) = 0.
    zw(jl, 4) = 0.
    pfd(jl, kflev+1) = ((1.-pclear(jl))*zrj(jl,jaj,kflev+1)+pclear(jl)*zrj0( &
      jl,jaj,kflev+1))*rsun(knu)
  END DO
  DO jk = 1, kflev
    ikl = kflev + 1 - jk
    DO jl = 1, kdlon
      zw(jl, 1) = zw(jl, 1) + pud(jl, 1, ikl)/zrmue(jl, ikl)
      zw(jl, 2) = zw(jl, 2) + poz(jl, ikl)/zrmue(jl, ikl)
      zw(jl, 3) = zw(jl, 3) + pud(jl, 1, ikl)/zrmu0(jl, ikl)
      zw(jl, 4) = zw(jl, 4) + poz(jl, ikl)/zrmu0(jl, ikl)
    END DO

    CALL swtt1_lmdar4(knu, 4, iind, zw, zr)

    DO jl = 1, kdlon
      zdiff(jl) = zr(jl, 1)*zr(jl, 2)*zrj(jl, jaj, ikl)
      zdirf(jl) = zr(jl, 3)*zr(jl, 4)*zrj0(jl, jaj, ikl)
      pfd(jl, ikl) = ((1.-pclear(jl))*zdiff(jl)+pclear(jl)*zdirf(jl))* &
        rsun(knu)
    END DO
  END DO

  ! *         3.2   UPWARD FLUXES
  ! -------------


  DO jl = 1, kdlon
    pfu(jl, 1) = ((1.-pclear(jl))*zdiff(jl)*palbd(jl,knu)+pclear(jl)*zdirf(jl &
      )*palbp(jl,knu))*rsun(knu)
  END DO

  DO jk = 2, kflev + 1
    ikm1 = jk - 1
    DO jl = 1, kdlon
      zw(jl, 1) = zw(jl, 1) + pud(jl, 1, ikm1)*1.66
      zw(jl, 2) = zw(jl, 2) + poz(jl, ikm1)*1.66
      zw(jl, 3) = zw(jl, 3) + pud(jl, 1, ikm1)*1.66
      zw(jl, 4) = zw(jl, 4) + poz(jl, ikm1)*1.66
    END DO

    CALL swtt1_lmdar4(knu, 4, iind, zw, zr)

    DO jl = 1, kdlon
      zdiff(jl) = zr(jl, 1)*zr(jl, 2)*zrk(jl, jaj, jk)
      zdirf(jl) = zr(jl, 3)*zr(jl, 4)*zrk0(jl, jaj, jk)
      pfu(jl, jk) = ((1.-pclear(jl))*zdiff(jl)+pclear(jl)*zdirf(jl))* &
        rsun(knu)
    END DO
  END DO

  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE sw1s_lmdar4
SUBROUTINE sw2s_lmdar4(knu, paer, flag_aer, tauae, pizae, cgae, paki, palbd, &
    palbp, pcg, pcld, pclear, pcldsw, pdsig, pomega, poz, prmu, psec, ptau, &
    pud, pwv, pqs, pfdown, pfup)
  USE dimphy
  USE radiation_ar4_param, ONLY: rsun, rray
  USE infotrac_phy, ONLY: type_trac
#ifdef REPROBUS
  USE chem_rep, ONLY: rsuntime, ok_suntime
#endif

  IMPLICIT NONE
  include "radepsi.h"

  ! ------------------------------------------------------------------
  ! PURPOSE.
  ! --------

  ! THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE
  ! SECOND SPECTRAL INTERVAL FOLLOWING FOUQUART AND BONNEL (1980).

  ! METHOD.
  ! -------

  ! 1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
  ! CONTINUUM SCATTERING
  ! 2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
  ! A GREY MOLECULAR ABSORPTION
  ! 3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
  ! OF ABSORBERS
  ! 4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
  ! 5. MULTIPLY BY OZONE TRANSMISSION FUNCTION

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
  ! DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! 94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
  ! ------------------------------------------------------------------
  ! * ARGUMENTS:

  INTEGER knu
  ! -OB
  REAL (KIND=8) flag_aer
  REAL (KIND=8) tauae(kdlon, kflev, 2)
  REAL (KIND=8) pizae(kdlon, kflev, 2)
  REAL (KIND=8) cgae(kdlon, kflev, 2)
  REAL (KIND=8) paer(kdlon, kflev, 5)
  REAL (KIND=8) paki(kdlon, 2)
  REAL (KIND=8) palbd(kdlon, 2)
  REAL (KIND=8) palbp(kdlon, 2)
  REAL (KIND=8) pcg(kdlon, 2, kflev)
  REAL (KIND=8) pcld(kdlon, kflev)
  REAL (KIND=8) pcldsw(kdlon, kflev)
  REAL (KIND=8) pclear(kdlon)
  REAL (KIND=8) pdsig(kdlon, kflev)
  REAL (KIND=8) pomega(kdlon, 2, kflev)
  REAL (KIND=8) poz(kdlon, kflev)
  REAL (KIND=8) pqs(kdlon, kflev)
  REAL (KIND=8) prmu(kdlon)
  REAL (KIND=8) psec(kdlon)
  REAL (KIND=8) ptau(kdlon, 2, kflev)
  REAL (KIND=8) pud(kdlon, 5, kflev+1)
  REAL (KIND=8) pwv(kdlon, kflev)

  REAL (KIND=8) pfdown(kdlon, kflev+1)
  REAL (KIND=8) pfup(kdlon, kflev+1)

  ! * LOCAL VARIABLES:

  INTEGER iind2(2), iind3(3)
  REAL (KIND=8) zcgaz(kdlon, kflev)
  REAL (KIND=8) zfd(kdlon, kflev+1)
  REAL (KIND=8) zfu(kdlon, kflev+1)
  REAL (KIND=8) zg(kdlon)
  REAL (KIND=8) zgg(kdlon)
  REAL (KIND=8) zpizaz(kdlon, kflev)
  REAL (KIND=8) zrayl(kdlon)
  REAL (KIND=8) zray1(kdlon, kflev+1)
  REAL (KIND=8) zray2(kdlon, kflev+1)
  REAL (KIND=8) zref(kdlon)
  REAL (KIND=8) zrefz(kdlon, 2, kflev+1)
  REAL (KIND=8) zre1(kdlon)
  REAL (KIND=8) zre2(kdlon)
  REAL (KIND=8) zrj(kdlon, 6, kflev+1)
  REAL (KIND=8) zrj0(kdlon, 6, kflev+1)
  REAL (KIND=8) zrk(kdlon, 6, kflev+1)
  REAL (KIND=8) zrk0(kdlon, 6, kflev+1)
  REAL (KIND=8) zrl(kdlon, 8)
  REAL (KIND=8) zrmue(kdlon, kflev+1)
  REAL (KIND=8) zrmu0(kdlon, kflev+1)
  REAL (KIND=8) zrmuz(kdlon)
  REAL (KIND=8) zrneb(kdlon)
  REAL (KIND=8) zruef(kdlon, 8)
  REAL (KIND=8) zr1(kdlon)
  REAL (KIND=8) zr2(kdlon, 2)
  REAL (KIND=8) zr3(kdlon, 3)
  REAL (KIND=8) zr4(kdlon)
  REAL (KIND=8) zr21(kdlon)
  REAL (KIND=8) zr22(kdlon)
  REAL (KIND=8) zs(kdlon)
  REAL (KIND=8) ztauaz(kdlon, kflev)
  REAL (KIND=8) zto1(kdlon)
  REAL (KIND=8) ztr(kdlon, 2, kflev+1)
  REAL (KIND=8) ztra1(kdlon, kflev+1)
  REAL (KIND=8) ztra2(kdlon, kflev+1)
  REAL (KIND=8) ztr1(kdlon)
  REAL (KIND=8) ztr2(kdlon)
  REAL (KIND=8) zw(kdlon)
  REAL (KIND=8) zw1(kdlon)
  REAL (KIND=8) zw2(kdlon, 2)
  REAL (KIND=8) zw3(kdlon, 3)
  REAL (KIND=8) zw4(kdlon)
  REAL (KIND=8) zw5(kdlon)

  INTEGER jl, jk, k, jaj, ikm1, ikl, jn, jabs, jkm1
  INTEGER jref, jkl, jklp1, jajp, jkki, jkkp4, jn2j, iabs
  REAL (KIND=8) zrmum1, zwh2o, zcneb, zaa, zbb, zrki, zre11

  ! If running with Reporbus, overwrite default values of RSUN.
  ! Otherwise keep default values from radiation_AR4_param module.
  IF (type_trac=='repr') THEN
#ifdef REPROBUS
    IF (ok_suntime) THEN
      rsun(1) = rsuntime(1)
      rsun(2) = rsuntime(2)
    END IF
#endif
  END IF

  ! ------------------------------------------------------------------

  ! *         1.     SECOND SPECTRAL INTERVAL (0.68-4.00 MICRON)
  ! -------------------------------------------



  ! *         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
  ! -----------------------------------------


  DO jl = 1, kdlon
    zrmum1 = 1. - prmu(jl)
    zrayl(jl) = rray(knu, 1) + zrmum1*(rray(knu,2)+zrmum1*(rray(knu, &
      3)+zrmum1*(rray(knu,4)+zrmum1*(rray(knu,5)+zrmum1*rray(knu,6)))))
  END DO

  ! ------------------------------------------------------------------

  ! *         2.    CONTINUUM SCATTERING CALCULATIONS
  ! ---------------------------------


  ! *         2.1   CLEAR-SKY FRACTION OF THE COLUMN
  ! --------------------------------


  CALL swclr_lmdar4(knu, paer, flag_aer, tauae, pizae, cgae, palbp, pdsig, &
    zrayl, psec, zcgaz, zpizaz, zray1, zray2, zrefz, zrj0, zrk0, zrmu0, &
    ztauaz, ztra1, ztra2)

  ! *         2.2   CLOUDY FRACTION OF THE COLUMN
  ! -----------------------------


  CALL swr_lmdar4(knu, palbd, pcg, pcld, pdsig, pomega, zrayl, psec, ptau, &
    zcgaz, zpizaz, zray1, zray2, zrefz, zrj, zrk, zrmue, ztauaz, ztra1, &
    ztra2)

  ! ------------------------------------------------------------------

  ! *         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
  ! ------------------------------------------------------


  jn = 2

  DO jabs = 1, 2
    ! *         3.1  SURFACE CONDITIONS
    ! ------------------


    DO jl = 1, kdlon
      zrefz(jl, 2, 1) = palbd(jl, knu)
      zrefz(jl, 1, 1) = palbd(jl, knu)
    END DO

    ! *         3.2  INTRODUCING CLOUD EFFECTS
    ! -------------------------


    DO jk = 2, kflev + 1
      jkm1 = jk - 1
      ikl = kflev + 1 - jkm1
      DO jl = 1, kdlon
        zrneb(jl) = pcld(jl, jkm1)
        IF (jabs==1 .AND. zrneb(jl)>2.*zeelog) THEN
          zwh2o = max(pwv(jl,jkm1), zeelog)
          zcneb = max(zeelog, min(zrneb(jl),1.-zeelog))
          zbb = pud(jl, jabs, jkm1)*pqs(jl, jkm1)/zwh2o
          zaa = max((pud(jl,jabs,jkm1)-zcneb*zbb)/(1.-zcneb), zeelog)
        ELSE
          zaa = pud(jl, jabs, jkm1)
          zbb = zaa
        END IF
        zrki = paki(jl, jabs)
        zs(jl) = exp(-zrki*zaa*1.66)
        zg(jl) = exp(-zrki*zaa/zrmue(jl,jk))
        ztr1(jl) = 0.
        zre1(jl) = 0.
        ztr2(jl) = 0.
        zre2(jl) = 0.

        zw(jl) = pomega(jl, knu, jkm1)
        zto1(jl) = ptau(jl, knu, jkm1)/zw(jl) + ztauaz(jl, jkm1)/zpizaz(jl, &
          jkm1) + zbb*zrki

        zr21(jl) = ptau(jl, knu, jkm1) + ztauaz(jl, jkm1)
        zr22(jl) = ptau(jl, knu, jkm1)/zr21(jl)
        zgg(jl) = zr22(jl)*pcg(jl, knu, jkm1) + (1.-zr22(jl))*zcgaz(jl, jkm1)
        zw(jl) = zr21(jl)/zto1(jl)
        zref(jl) = zrefz(jl, 1, jkm1)
        zrmuz(jl) = zrmue(jl, jk)
      END DO

      CALL swde_lmdar4(zgg, zref, zrmuz, zto1, zw, zre1, zre2, ztr1, ztr2)

      DO jl = 1, kdlon

        zrefz(jl, 2, jk) = (1.-zrneb(jl))*(zray1(jl,jkm1)+zrefz(jl,2,jkm1)* &
          ztra1(jl,jkm1)*ztra2(jl,jkm1))*zg(jl)*zs(jl) + zrneb(jl)*zre1(jl)

        ztr(jl, 2, jkm1) = zrneb(jl)*ztr1(jl) + (ztra1(jl,jkm1))*zg(jl)*(1.- &
          zrneb(jl))

        zrefz(jl, 1, jk) = (1.-zrneb(jl))*(zray1(jl,jkm1)+zrefz(jl,1,jkm1)* &
          ztra1(jl,jkm1)*ztra2(jl,jkm1)/(1.-zray2(jl,jkm1)*zrefz(jl,1, &
          jkm1)))*zg(jl)*zs(jl) + zrneb(jl)*zre2(jl)

        ztr(jl, 1, jkm1) = zrneb(jl)*ztr2(jl) + (ztra1(jl,jkm1)/(1.-zray2(jl, &
          jkm1)*zrefz(jl,1,jkm1)))*zg(jl)*(1.-zrneb(jl))

      END DO
    END DO

    ! *         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
    ! -------------------------------------------------


    DO jref = 1, 2

      jn = jn + 1

      DO jl = 1, kdlon
        zrj(jl, jn, kflev+1) = 1.
        zrk(jl, jn, kflev+1) = zrefz(jl, jref, kflev+1)
      END DO

      DO jk = 1, kflev
        jkl = kflev + 1 - jk
        jklp1 = jkl + 1
        DO jl = 1, kdlon
          zre11 = zrj(jl, jn, jklp1)*ztr(jl, jref, jkl)
          zrj(jl, jn, jkl) = zre11
          zrk(jl, jn, jkl) = zre11*zrefz(jl, jref, jkl)
        END DO
      END DO
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         4.    INVERT GREY AND CONTINUUM FLUXES
  ! --------------------------------



  ! *         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
  ! ---------------------------------------------


  DO jk = 1, kflev + 1
    DO jaj = 1, 5, 2
      jajp = jaj + 1
      DO jl = 1, kdlon
        zrj(jl, jaj, jk) = zrj(jl, jaj, jk) - zrj(jl, jajp, jk)
        zrk(jl, jaj, jk) = zrk(jl, jaj, jk) - zrk(jl, jajp, jk)
        zrj(jl, jaj, jk) = max(zrj(jl,jaj,jk), zeelog)
        zrk(jl, jaj, jk) = max(zrk(jl,jaj,jk), zeelog)
      END DO
    END DO
  END DO

  DO jk = 1, kflev + 1
    DO jaj = 2, 6, 2
      DO jl = 1, kdlon
        zrj(jl, jaj, jk) = max(zrj(jl,jaj,jk), zeelog)
        zrk(jl, jaj, jk) = max(zrk(jl,jaj,jk), zeelog)
      END DO
    END DO
  END DO

  ! *         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
  ! ---------------------------------------------


  DO jk = 1, kflev + 1
    jkki = 1
    DO jaj = 1, 2
      iind2(1) = jaj
      iind2(2) = jaj
      DO jn = 1, 2
        jn2j = jn + 2*jaj
        jkkp4 = jkki + 4

        ! *         4.2.1  EFFECTIVE ABSORBER AMOUNTS
        ! --------------------------


        DO jl = 1, kdlon
          zw2(jl, 1) = log(zrj(jl,jn,jk)/zrj(jl,jn2j,jk))/paki(jl, jaj)
          zw2(jl, 2) = log(zrk(jl,jn,jk)/zrk(jl,jn2j,jk))/paki(jl, jaj)
        END DO

        ! *         4.2.2  TRANSMISSION FUNCTION
        ! ---------------------


        CALL swtt1_lmdar4(knu, 2, iind2, zw2, zr2)

        DO jl = 1, kdlon
          zrl(jl, jkki) = zr2(jl, 1)
          zruef(jl, jkki) = zw2(jl, 1)
          zrl(jl, jkkp4) = zr2(jl, 2)
          zruef(jl, jkkp4) = zw2(jl, 2)
        END DO

        jkki = jkki + 1
      END DO
    END DO

    ! *         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
    ! ------------------------------------------------------


    DO jl = 1, kdlon
      pfdown(jl, jk) = zrj(jl, 1, jk)*zrl(jl, 1)*zrl(jl, 3) + &
        zrj(jl, 2, jk)*zrl(jl, 2)*zrl(jl, 4)
      pfup(jl, jk) = zrk(jl, 1, jk)*zrl(jl, 5)*zrl(jl, 7) + &
        zrk(jl, 2, jk)*zrl(jl, 6)*zrl(jl, 8)
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
  ! ----------------------------------------



  ! *         5.1   DOWNWARD FLUXES
  ! ---------------


  jaj = 2
  iind3(1) = 1
  iind3(2) = 2
  iind3(3) = 3

  DO jl = 1, kdlon
    zw3(jl, 1) = 0.
    zw3(jl, 2) = 0.
    zw3(jl, 3) = 0.
    zw4(jl) = 0.
    zw5(jl) = 0.
    zr4(jl) = 1.
    zfd(jl, kflev+1) = zrj0(jl, jaj, kflev+1)
  END DO
  DO jk = 1, kflev
    ikl = kflev + 1 - jk
    DO jl = 1, kdlon
      zw3(jl, 1) = zw3(jl, 1) + pud(jl, 1, ikl)/zrmu0(jl, ikl)
      zw3(jl, 2) = zw3(jl, 2) + pud(jl, 2, ikl)/zrmu0(jl, ikl)
      zw3(jl, 3) = zw3(jl, 3) + poz(jl, ikl)/zrmu0(jl, ikl)
      zw4(jl) = zw4(jl) + pud(jl, 4, ikl)/zrmu0(jl, ikl)
      zw5(jl) = zw5(jl) + pud(jl, 5, ikl)/zrmu0(jl, ikl)
    END DO

    CALL swtt1_lmdar4(knu, 3, iind3, zw3, zr3)

    DO jl = 1, kdlon
      ! ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
      zfd(jl, ikl) = zr3(jl, 1)*zr3(jl, 2)*zr3(jl, 3)*zr4(jl)* &
        zrj0(jl, jaj, ikl)
    END DO
  END DO

  ! *         5.2   UPWARD FLUXES
  ! -------------


  DO jl = 1, kdlon
    zfu(jl, 1) = zfd(jl, 1)*palbp(jl, knu)
  END DO

  DO jk = 2, kflev + 1
    ikm1 = jk - 1
    DO jl = 1, kdlon
      zw3(jl, 1) = zw3(jl, 1) + pud(jl, 1, ikm1)*1.66
      zw3(jl, 2) = zw3(jl, 2) + pud(jl, 2, ikm1)*1.66
      zw3(jl, 3) = zw3(jl, 3) + poz(jl, ikm1)*1.66
      zw4(jl) = zw4(jl) + pud(jl, 4, ikm1)*1.66
      zw5(jl) = zw5(jl) + pud(jl, 5, ikm1)*1.66
    END DO

    CALL swtt1_lmdar4(knu, 3, iind3, zw3, zr3)

    DO jl = 1, kdlon
      ! ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
      zfu(jl, jk) = zr3(jl, 1)*zr3(jl, 2)*zr3(jl, 3)*zr4(jl)* &
        zrk0(jl, jaj, jk)
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         6.     INTRODUCTION OF OZONE AND H2O CONTINUUM ABSORPTION
  ! --------------------------------------------------

  iabs = 3

  ! *         6.1    DOWNWARD FLUXES
  ! ---------------

  DO jl = 1, kdlon
    zw1(jl) = 0.
    zw4(jl) = 0.
    zw5(jl) = 0.
    zr1(jl) = 0.
    pfdown(jl, kflev+1) = ((1.-pclear(jl))*pfdown(jl,kflev+1)+pclear(jl)*zfd( &
      jl,kflev+1))*rsun(knu)
  END DO

  DO jk = 1, kflev
    ikl = kflev + 1 - jk
    DO jl = 1, kdlon
      zw1(jl) = zw1(jl) + poz(jl, ikl)/zrmue(jl, ikl)
      zw4(jl) = zw4(jl) + pud(jl, 4, ikl)/zrmue(jl, ikl)
      zw5(jl) = zw5(jl) + pud(jl, 5, ikl)/zrmue(jl, ikl)
      ! ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
    END DO

    CALL swtt_lmdar4(knu, iabs, zw1, zr1)

    DO jl = 1, kdlon
      pfdown(jl, ikl) = ((1.-pclear(jl))*zr1(jl)*zr4(jl)*pfdown(jl,ikl)+ &
        pclear(jl)*zfd(jl,ikl))*rsun(knu)
    END DO
  END DO

  ! *         6.2    UPWARD FLUXES
  ! -------------

  DO jl = 1, kdlon
    pfup(jl, 1) = ((1.-pclear(jl))*zr1(jl)*zr4(jl)*pfup(jl,1)+pclear(jl)*zfu( &
      jl,1))*rsun(knu)
  END DO

  DO jk = 2, kflev + 1
    ikm1 = jk - 1
    DO jl = 1, kdlon
      zw1(jl) = zw1(jl) + poz(jl, ikm1)*1.66
      zw4(jl) = zw4(jl) + pud(jl, 4, ikm1)*1.66
      zw5(jl) = zw5(jl) + pud(jl, 5, ikm1)*1.66
      ! ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
    END DO

    CALL swtt_lmdar4(knu, iabs, zw1, zr1)

    DO jl = 1, kdlon
      pfup(jl, jk) = ((1.-pclear(jl))*zr1(jl)*zr4(jl)*pfup(jl,jk)+pclear(jl)* &
        zfu(jl,jk))*rsun(knu)
    END DO
  END DO

  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE sw2s_lmdar4
SUBROUTINE swclr_lmdar4(knu, paer, flag_aer, tauae, pizae, cgae, palbp, &
    pdsig, prayl, psec, pcgaz, ppizaz, pray1, pray2, prefz, prj, prk, prmu0, &
    ptauaz, ptra1, ptra2)
  USE dimphy
  USE radiation_ar4_param, ONLY: taua, rpiza, rcga
  IMPLICIT NONE
  include "radepsi.h"
  include "radopt.h"

  ! ------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
  ! CLEAR-SKY COLUMN

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
  ! DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 94-11-15
  ! ------------------------------------------------------------------
  ! * ARGUMENTS:

  INTEGER knu
  ! -OB
  REAL (KIND=8) flag_aer
  REAL (KIND=8) tauae(kdlon, kflev, 2)
  REAL (KIND=8) pizae(kdlon, kflev, 2)
  REAL (KIND=8) cgae(kdlon, kflev, 2)
  REAL (KIND=8) paer(kdlon, kflev, 5)
  REAL (KIND=8) palbp(kdlon, 2)
  REAL (KIND=8) pdsig(kdlon, kflev)
  REAL (KIND=8) prayl(kdlon)
  REAL (KIND=8) psec(kdlon)

  REAL (KIND=8) pcgaz(kdlon, kflev)
  REAL (KIND=8) ppizaz(kdlon, kflev)
  REAL (KIND=8) pray1(kdlon, kflev+1)
  REAL (KIND=8) pray2(kdlon, kflev+1)
  REAL (KIND=8) prefz(kdlon, 2, kflev+1)
  REAL (KIND=8) prj(kdlon, 6, kflev+1)
  REAL (KIND=8) prk(kdlon, 6, kflev+1)
  REAL (KIND=8) prmu0(kdlon, kflev+1)
  REAL (KIND=8) ptauaz(kdlon, kflev)
  REAL (KIND=8) ptra1(kdlon, kflev+1)
  REAL (KIND=8) ptra2(kdlon, kflev+1)

  ! * LOCAL VARIABLES:

  REAL (KIND=8) zc0i(kdlon, kflev+1)
  REAL (KIND=8) zcle0(kdlon, kflev)
  REAL (KIND=8) zclear(kdlon)
  REAL (KIND=8) zr21(kdlon)
  REAL (KIND=8) zr23(kdlon)
  REAL (KIND=8) zss0(kdlon)
  REAL (KIND=8) zscat(kdlon)
  REAL (KIND=8) ztr(kdlon, 2, kflev+1)

  INTEGER jl, jk, ja, jae, jkl, jklp1, jaj, jkm1, in
  REAL (KIND=8) ztray, zgar, zratio, zff, zfacoa, zcorae
  REAL (KIND=8) zmue, zgap, zww, zto, zden, zmu1, zden1
  REAL (KIND=8) zbmu0, zbmu1, zre11

  ! ------------------------------------------------------------------

  ! *         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
  ! --------------------------------------------


  ! cdir collapse
  DO jk = 1, kflev + 1
    DO ja = 1, 6
      DO jl = 1, kdlon
        prj(jl, ja, jk) = 0.
        prk(jl, ja, jk) = 0.
      END DO
    END DO
  END DO

  DO jk = 1, kflev
    ! -OB
    ! DO 104 JL = 1, KDLON
    ! PCGAZ(JL,JK) = 0.
    ! PPIZAZ(JL,JK) =  0.
    ! PTAUAZ(JL,JK) = 0.
    ! 104  CONTINUE
    ! -OB
    ! DO 106 JAE=1,5
    ! DO 105 JL = 1, KDLON
    ! PTAUAZ(JL,JK)=PTAUAZ(JL,JK)
    ! S        +PAER(JL,JK,JAE)*TAUA(KNU,JAE)
    ! PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL,JK,JAE)
    ! S        * TAUA(KNU,JAE)*RPIZA(KNU,JAE)
    ! PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL,JK,JAE)
    ! S        * TAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)
    ! 105  CONTINUE
    ! 106  CONTINUE
    ! -OB
    DO jl = 1, kdlon
      ptauaz(jl, jk) = flag_aer*tauae(jl, jk, knu)
      ppizaz(jl, jk) = flag_aer*pizae(jl, jk, knu)
      pcgaz(jl, jk) = flag_aer*cgae(jl, jk, knu)
    END DO

    IF (flag_aer>0) THEN
      ! -OB
      DO jl = 1, kdlon
        ! PCGAZ(JL,JK)=PCGAZ(JL,JK)/PPIZAZ(JL,JK)
        ! PPIZAZ(JL,JK)=PPIZAZ(JL,JK)/PTAUAZ(JL,JK)
        ztray = prayl(jl)*pdsig(jl, jk)
        zratio = ztray/(ztray+ptauaz(jl,jk))
        zgar = pcgaz(jl, jk)
        zff = zgar*zgar
        ptauaz(jl, jk) = ztray + ptauaz(jl, jk)*(1.-ppizaz(jl,jk)*zff)
        pcgaz(jl, jk) = zgar*(1.-zratio)/(1.+zgar)
        ppizaz(jl, jk) = zratio + (1.-zratio)*ppizaz(jl, jk)*(1.-zff)/(1.- &
          ppizaz(jl,jk)*zff)
      END DO
    ELSE
      DO jl = 1, kdlon
        ztray = prayl(jl)*pdsig(jl, jk)
        ptauaz(jl, jk) = ztray
        pcgaz(jl, jk) = 0.
        ppizaz(jl, jk) = 1. - repsct
      END DO
    END IF ! check flag_aer
    ! 107  CONTINUE
    ! PRINT 9107,JK,((PAER(JL,JK,JAE),JAE=1,5)
    ! $ ,PTAUAZ(JL,JK),PPIZAZ(JL,JK),PCGAZ(JL,JK),JL=1,KDLON)
    ! 9107 FORMAT(1X,'SWCLR_107',I3,8E12.5)

  END DO

  ! ------------------------------------------------------------------

  ! *         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
  ! ----------------------------------------------


  DO jl = 1, kdlon
    zr23(jl) = 0.
    zc0i(jl, kflev+1) = 0.
    zclear(jl) = 1.
    zscat(jl) = 0.
  END DO

  jk = 1
  jkl = kflev + 1 - jk
  jklp1 = jkl + 1
  DO jl = 1, kdlon
    zfacoa = 1. - ppizaz(jl, jkl)*pcgaz(jl, jkl)*pcgaz(jl, jkl)
    zcorae = zfacoa*ptauaz(jl, jkl)*psec(jl)
    zr21(jl) = exp(-zcorae)
    zss0(jl) = 1. - zr21(jl)
    zcle0(jl, jkl) = zss0(jl)

    IF (novlp==1) THEN
      ! * maximum-random
      zclear(jl) = zclear(jl)*(1.0-max(zss0(jl),zscat(jl)))/ &
        (1.0-min(zscat(jl),1.-zepsec))
      zc0i(jl, jkl) = 1.0 - zclear(jl)
      zscat(jl) = zss0(jl)
    ELSE IF (novlp==2) THEN
      ! * maximum
      zscat(jl) = max(zss0(jl), zscat(jl))
      zc0i(jl, jkl) = zscat(jl)
    ELSE IF (novlp==3) THEN
      ! * random
      zclear(jl) = zclear(jl)*(1.0-zss0(jl))
      zscat(jl) = 1.0 - zclear(jl)
      zc0i(jl, jkl) = zscat(jl)
    END IF
  END DO

  DO jk = 2, kflev
    jkl = kflev + 1 - jk
    jklp1 = jkl + 1
    DO jl = 1, kdlon
      zfacoa = 1. - ppizaz(jl, jkl)*pcgaz(jl, jkl)*pcgaz(jl, jkl)
      zcorae = zfacoa*ptauaz(jl, jkl)*psec(jl)
      zr21(jl) = exp(-zcorae)
      zss0(jl) = 1. - zr21(jl)
      zcle0(jl, jkl) = zss0(jl)

      IF (novlp==1) THEN
        ! * maximum-random
        zclear(jl) = zclear(jl)*(1.0-max(zss0(jl),zscat(jl)))/ &
          (1.0-min(zscat(jl),1.-zepsec))
        zc0i(jl, jkl) = 1.0 - zclear(jl)
        zscat(jl) = zss0(jl)
      ELSE IF (novlp==2) THEN
        ! * maximum
        zscat(jl) = max(zss0(jl), zscat(jl))
        zc0i(jl, jkl) = zscat(jl)
      ELSE IF (novlp==3) THEN
        ! * random
        zclear(jl) = zclear(jl)*(1.0-zss0(jl))
        zscat(jl) = 1.0 - zclear(jl)
        zc0i(jl, jkl) = zscat(jl)
      END IF
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
  ! -----------------------------------------------


  DO jl = 1, kdlon
    pray1(jl, kflev+1) = 0.
    pray2(jl, kflev+1) = 0.
    prefz(jl, 2, 1) = palbp(jl, knu)
    prefz(jl, 1, 1) = palbp(jl, knu)
    ptra1(jl, kflev+1) = 1.
    ptra2(jl, kflev+1) = 1.
  END DO

  DO jk = 2, kflev + 1
    jkm1 = jk - 1
    DO jl = 1, kdlon

      ! ------------------------------------------------------------------

      ! *         3.1  EQUIVALENT ZENITH ANGLE
      ! -----------------------


      zmue = (1.-zc0i(jl,jk))*psec(jl) + zc0i(jl, jk)*1.66
      prmu0(jl, jk) = 1./zmue

      ! ------------------------------------------------------------------

      ! *         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
      ! ----------------------------------------------------


      zgap = pcgaz(jl, jkm1)
      zbmu0 = 0.5 - 0.75*zgap/zmue
      zww = ppizaz(jl, jkm1)
      zto = ptauaz(jl, jkm1)
      zden = 1. + (1.-zww+zbmu0*zww)*zto*zmue + (1-zww)*(1.-zww+2.*zbmu0*zww) &
        *zto*zto*zmue*zmue
      pray1(jl, jkm1) = zbmu0*zww*zto*zmue/zden
      ptra1(jl, jkm1) = 1./zden

      zmu1 = 0.5
      zbmu1 = 0.5 - 0.75*zgap*zmu1
      zden1 = 1. + (1.-zww+zbmu1*zww)*zto/zmu1 + (1-zww)*(1.-zww+2.*zbmu1*zww &
        )*zto*zto/zmu1/zmu1
      pray2(jl, jkm1) = zbmu1*zww*zto/zmu1/zden1
      ptra2(jl, jkm1) = 1./zden1



      prefz(jl, 1, jk) = (pray1(jl,jkm1)+prefz(jl,1,jkm1)*ptra1(jl,jkm1)* &
        ptra2(jl,jkm1)/(1.-pray2(jl,jkm1)*prefz(jl,1,jkm1)))

      ztr(jl, 1, jkm1) = (ptra1(jl,jkm1)/(1.-pray2(jl,jkm1)*prefz(jl,1, &
        jkm1)))

      prefz(jl, 2, jk) = (pray1(jl,jkm1)+prefz(jl,2,jkm1)*ptra1(jl,jkm1)* &
        ptra2(jl,jkm1))

      ztr(jl, 2, jkm1) = ptra1(jl, jkm1)

    END DO
  END DO
  DO jl = 1, kdlon
    zmue = (1.-zc0i(jl,1))*psec(jl) + zc0i(jl, 1)*1.66
    prmu0(jl, 1) = 1./zmue
  END DO

  ! ------------------------------------------------------------------

  ! *         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
  ! -------------------------------------------------


  IF (knu==1) THEN
    jaj = 2
    DO jl = 1, kdlon
      prj(jl, jaj, kflev+1) = 1.
      prk(jl, jaj, kflev+1) = prefz(jl, 1, kflev+1)
    END DO

    DO jk = 1, kflev
      jkl = kflev + 1 - jk
      jklp1 = jkl + 1
      DO jl = 1, kdlon
        zre11 = prj(jl, jaj, jklp1)*ztr(jl, 1, jkl)
        prj(jl, jaj, jkl) = zre11
        prk(jl, jaj, jkl) = zre11*prefz(jl, 1, jkl)
      END DO
    END DO

  ELSE

    DO jaj = 1, 2
      DO jl = 1, kdlon
        prj(jl, jaj, kflev+1) = 1.
        prk(jl, jaj, kflev+1) = prefz(jl, jaj, kflev+1)
      END DO

      DO jk = 1, kflev
        jkl = kflev + 1 - jk
        jklp1 = jkl + 1
        DO jl = 1, kdlon
          zre11 = prj(jl, jaj, jklp1)*ztr(jl, jaj, jkl)
          prj(jl, jaj, jkl) = zre11
          prk(jl, jaj, jkl) = zre11*prefz(jl, jaj, jkl)
        END DO
      END DO
    END DO

  END IF

  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE swclr_lmdar4
SUBROUTINE swr_lmdar4(knu, palbd, pcg, pcld, pdsig, pomega, prayl, psec, &
    ptau, pcgaz, ppizaz, pray1, pray2, prefz, prj, prk, prmue, ptauaz, ptra1, &
    ptra2)
  USE dimphy
  IMPLICIT NONE
  include "radepsi.h"
  include "radopt.h"

  ! ------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
  ! CONTINUUM SCATTERING

  ! METHOD.
  ! -------

  ! 1. COMPUTES CONTINUUM FLUXES CORRESPONDING TO AEROSOL
  ! OR/AND RAYLEIGH SCATTERING (NO MOLECULAR GAS ABSORPTION)

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
  ! DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! ------------------------------------------------------------------
  ! * ARGUMENTS:

  INTEGER knu
  REAL (KIND=8) palbd(kdlon, 2)
  REAL (KIND=8) pcg(kdlon, 2, kflev)
  REAL (KIND=8) pcld(kdlon, kflev)
  REAL (KIND=8) pdsig(kdlon, kflev)
  REAL (KIND=8) pomega(kdlon, 2, kflev)
  REAL (KIND=8) prayl(kdlon)
  REAL (KIND=8) psec(kdlon)
  REAL (KIND=8) ptau(kdlon, 2, kflev)

  REAL (KIND=8) pray1(kdlon, kflev+1)
  REAL (KIND=8) pray2(kdlon, kflev+1)
  REAL (KIND=8) prefz(kdlon, 2, kflev+1)
  REAL (KIND=8) prj(kdlon, 6, kflev+1)
  REAL (KIND=8) prk(kdlon, 6, kflev+1)
  REAL (KIND=8) prmue(kdlon, kflev+1)
  REAL (KIND=8) pcgaz(kdlon, kflev)
  REAL (KIND=8) ppizaz(kdlon, kflev)
  REAL (KIND=8) ptauaz(kdlon, kflev)
  REAL (KIND=8) ptra1(kdlon, kflev+1)
  REAL (KIND=8) ptra2(kdlon, kflev+1)

  ! * LOCAL VARIABLES:

  REAL (KIND=8) zc1i(kdlon, kflev+1)
  REAL (KIND=8) zcleq(kdlon, kflev)
  REAL (KIND=8) zclear(kdlon)
  REAL (KIND=8) zcloud(kdlon)
  REAL (KIND=8) zgg(kdlon)
  REAL (KIND=8) zref(kdlon)
  REAL (KIND=8) zre1(kdlon)
  REAL (KIND=8) zre2(kdlon)
  REAL (KIND=8) zrmuz(kdlon)
  REAL (KIND=8) zrneb(kdlon)
  REAL (KIND=8) zr21(kdlon)
  REAL (KIND=8) zr22(kdlon)
  REAL (KIND=8) zr23(kdlon)
  REAL (KIND=8) zss1(kdlon)
  REAL (KIND=8) zto1(kdlon)
  REAL (KIND=8) ztr(kdlon, 2, kflev+1)
  REAL (KIND=8) ztr1(kdlon)
  REAL (KIND=8) ztr2(kdlon)
  REAL (KIND=8) zw(kdlon)

  INTEGER jk, jl, ja, jkl, jklp1, jkm1, jaj
  REAL (KIND=8) zfacoa, zfacoc, zcorae, zcorcd
  REAL (KIND=8) zmue, zgap, zww, zto, zden, zden1
  REAL (KIND=8) zmu1, zre11, zbmu0, zbmu1

  ! ------------------------------------------------------------------

  ! *         1.    INITIALIZATION
  ! --------------


  DO jk = 1, kflev + 1
    DO ja = 1, 6
      DO jl = 1, kdlon
        prj(jl, ja, jk) = 0.
        prk(jl, ja, jk) = 0.
      END DO
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
  ! ----------------------------------------------


  DO jl = 1, kdlon
    zr23(jl) = 0.
    zc1i(jl, kflev+1) = 0.
    zclear(jl) = 1.
    zcloud(jl) = 0.
  END DO

  jk = 1
  jkl = kflev + 1 - jk
  jklp1 = jkl + 1
  DO jl = 1, kdlon
    zfacoa = 1. - ppizaz(jl, jkl)*pcgaz(jl, jkl)*pcgaz(jl, jkl)
    zfacoc = 1. - pomega(jl, knu, jkl)*pcg(jl, knu, jkl)*pcg(jl, knu, jkl)
    zcorae = zfacoa*ptauaz(jl, jkl)*psec(jl)
    zcorcd = zfacoc*ptau(jl, knu, jkl)*psec(jl)
    zr21(jl) = exp(-zcorae)
    zr22(jl) = exp(-zcorcd)
    zss1(jl) = pcld(jl, jkl)*(1.0-zr21(jl)*zr22(jl)) + &
      (1.0-pcld(jl,jkl))*(1.0-zr21(jl))
    zcleq(jl, jkl) = zss1(jl)

    IF (novlp==1) THEN
      ! * maximum-random
      zclear(jl) = zclear(jl)*(1.0-max(zss1(jl),zcloud(jl)))/ &
        (1.0-min(zcloud(jl),1.-zepsec))
      zc1i(jl, jkl) = 1.0 - zclear(jl)
      zcloud(jl) = zss1(jl)
    ELSE IF (novlp==2) THEN
      ! * maximum
      zcloud(jl) = max(zss1(jl), zcloud(jl))
      zc1i(jl, jkl) = zcloud(jl)
    ELSE IF (novlp==3) THEN
      ! * random
      zclear(jl) = zclear(jl)*(1.0-zss1(jl))
      zcloud(jl) = 1.0 - zclear(jl)
      zc1i(jl, jkl) = zcloud(jl)
    END IF
  END DO

  DO jk = 2, kflev
    jkl = kflev + 1 - jk
    jklp1 = jkl + 1
    DO jl = 1, kdlon
      zfacoa = 1. - ppizaz(jl, jkl)*pcgaz(jl, jkl)*pcgaz(jl, jkl)
      zfacoc = 1. - pomega(jl, knu, jkl)*pcg(jl, knu, jkl)*pcg(jl, knu, jkl)
      zcorae = zfacoa*ptauaz(jl, jkl)*psec(jl)
      zcorcd = zfacoc*ptau(jl, knu, jkl)*psec(jl)
      zr21(jl) = exp(-zcorae)
      zr22(jl) = exp(-zcorcd)
      zss1(jl) = pcld(jl, jkl)*(1.0-zr21(jl)*zr22(jl)) + &
        (1.0-pcld(jl,jkl))*(1.0-zr21(jl))
      zcleq(jl, jkl) = zss1(jl)

      IF (novlp==1) THEN
        ! * maximum-random
        zclear(jl) = zclear(jl)*(1.0-max(zss1(jl),zcloud(jl)))/ &
          (1.0-min(zcloud(jl),1.-zepsec))
        zc1i(jl, jkl) = 1.0 - zclear(jl)
        zcloud(jl) = zss1(jl)
      ELSE IF (novlp==2) THEN
        ! * maximum
        zcloud(jl) = max(zss1(jl), zcloud(jl))
        zc1i(jl, jkl) = zcloud(jl)
      ELSE IF (novlp==3) THEN
        ! * random
        zclear(jl) = zclear(jl)*(1.0-zss1(jl))
        zcloud(jl) = 1.0 - zclear(jl)
        zc1i(jl, jkl) = zcloud(jl)
      END IF
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
  ! -----------------------------------------------


  DO jl = 1, kdlon
    pray1(jl, kflev+1) = 0.
    pray2(jl, kflev+1) = 0.
    prefz(jl, 2, 1) = palbd(jl, knu)
    prefz(jl, 1, 1) = palbd(jl, knu)
    ptra1(jl, kflev+1) = 1.
    ptra2(jl, kflev+1) = 1.
  END DO

  DO jk = 2, kflev + 1
    jkm1 = jk - 1
    DO jl = 1, kdlon
      zrneb(jl) = pcld(jl, jkm1)
      zre1(jl) = 0.
      ztr1(jl) = 0.
      zre2(jl) = 0.
      ztr2(jl) = 0.

      ! ------------------------------------------------------------------

      ! *         3.1  EQUIVALENT ZENITH ANGLE
      ! -----------------------


      zmue = (1.-zc1i(jl,jk))*psec(jl) + zc1i(jl, jk)*1.66
      prmue(jl, jk) = 1./zmue

      ! ------------------------------------------------------------------

      ! *         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
      ! ----------------------------------------------------


      zgap = pcgaz(jl, jkm1)
      zbmu0 = 0.5 - 0.75*zgap/zmue
      zww = ppizaz(jl, jkm1)
      zto = ptauaz(jl, jkm1)
      zden = 1. + (1.-zww+zbmu0*zww)*zto*zmue + (1-zww)*(1.-zww+2.*zbmu0*zww) &
        *zto*zto*zmue*zmue
      pray1(jl, jkm1) = zbmu0*zww*zto*zmue/zden
      ptra1(jl, jkm1) = 1./zden
      ! PRINT *,' LOOP 342 ** 3 ** JL=',JL,PRAY1(JL,JKM1),PTRA1(JL,JKM1)

      zmu1 = 0.5
      zbmu1 = 0.5 - 0.75*zgap*zmu1
      zden1 = 1. + (1.-zww+zbmu1*zww)*zto/zmu1 + (1-zww)*(1.-zww+2.*zbmu1*zww &
        )*zto*zto/zmu1/zmu1
      pray2(jl, jkm1) = zbmu1*zww*zto/zmu1/zden1
      ptra2(jl, jkm1) = 1./zden1

      ! ------------------------------------------------------------------

      ! *         3.3  EFFECT OF CLOUD LAYER
      ! ---------------------


      zw(jl) = pomega(jl, knu, jkm1)
      zto1(jl) = ptau(jl, knu, jkm1)/zw(jl) + ptauaz(jl, jkm1)/ppizaz(jl, &
        jkm1)
      zr21(jl) = ptau(jl, knu, jkm1) + ptauaz(jl, jkm1)
      zr22(jl) = ptau(jl, knu, jkm1)/zr21(jl)
      zgg(jl) = zr22(jl)*pcg(jl, knu, jkm1) + (1.-zr22(jl))*pcgaz(jl, jkm1)
      ! Modif PhD - JJM 19/03/96 pour erreurs arrondis
      ! machine
      ! PHD PROTECTION ZW(JL) = ZR21(JL) / ZTO1(JL)
      IF (zw(jl)==1. .AND. ppizaz(jl,jkm1)==1.) THEN
        zw(jl) = 1.
      ELSE
        zw(jl) = zr21(jl)/zto1(jl)
      END IF
      zref(jl) = prefz(jl, 1, jkm1)
      zrmuz(jl) = prmue(jl, jk)
    END DO

    CALL swde_lmdar4(zgg, zref, zrmuz, zto1, zw, zre1, zre2, ztr1, ztr2)

    DO jl = 1, kdlon

      prefz(jl, 1, jk) = (1.-zrneb(jl))*(pray1(jl,jkm1)+prefz(jl,1,jkm1)* &
        ptra1(jl,jkm1)*ptra2(jl,jkm1)/(1.-pray2(jl,jkm1)*prefz(jl,1, &
        jkm1))) + zrneb(jl)*zre2(jl)

      ztr(jl, 1, jkm1) = zrneb(jl)*ztr2(jl) + (ptra1(jl,jkm1)/(1.-pray2(jl, &
        jkm1)*prefz(jl,1,jkm1)))*(1.-zrneb(jl))

      prefz(jl, 2, jk) = (1.-zrneb(jl))*(pray1(jl,jkm1)+prefz(jl,2,jkm1)* &
        ptra1(jl,jkm1)*ptra2(jl,jkm1)) + zrneb(jl)*zre1(jl)

      ztr(jl, 2, jkm1) = zrneb(jl)*ztr1(jl) + ptra1(jl, jkm1)*(1.-zrneb(jl))

    END DO
  END DO
  DO jl = 1, kdlon
    zmue = (1.-zc1i(jl,1))*psec(jl) + zc1i(jl, 1)*1.66
    prmue(jl, 1) = 1./zmue
  END DO

  ! ------------------------------------------------------------------

  ! *         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
  ! -------------------------------------------------


  IF (knu==1) THEN
    jaj = 2
    DO jl = 1, kdlon
      prj(jl, jaj, kflev+1) = 1.
      prk(jl, jaj, kflev+1) = prefz(jl, 1, kflev+1)
    END DO

    DO jk = 1, kflev
      jkl = kflev + 1 - jk
      jklp1 = jkl + 1
      DO jl = 1, kdlon
        zre11 = prj(jl, jaj, jklp1)*ztr(jl, 1, jkl)
        prj(jl, jaj, jkl) = zre11
        prk(jl, jaj, jkl) = zre11*prefz(jl, 1, jkl)
      END DO
    END DO

  ELSE

    DO jaj = 1, 2
      DO jl = 1, kdlon
        prj(jl, jaj, kflev+1) = 1.
        prk(jl, jaj, kflev+1) = prefz(jl, jaj, kflev+1)
      END DO

      DO jk = 1, kflev
        jkl = kflev + 1 - jk
        jklp1 = jkl + 1
        DO jl = 1, kdlon
          zre11 = prj(jl, jaj, jklp1)*ztr(jl, jaj, jkl)
          prj(jl, jaj, jkl) = zre11
          prk(jl, jaj, jkl) = zre11*prefz(jl, jaj, jkl)
        END DO
      END DO
    END DO

  END IF

  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE swr_lmdar4
SUBROUTINE swde_lmdar4(pgg, pref, prmuz, pto1, pw, pre1, pre2, ptr1, ptr2)
  USE dimphy
  IMPLICIT NONE

  ! ------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY
  ! LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.

  ! METHOD.
  ! -------

  ! STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 88-12-15
  ! ------------------------------------------------------------------
  ! * ARGUMENTS:

  REAL (KIND=8) pgg(kdlon) ! ASSYMETRY FACTOR
  REAL (KIND=8) pref(kdlon) ! REFLECTIVITY OF THE UNDERLYING LAYER
  REAL (KIND=8) prmuz(kdlon) ! COSINE OF SOLAR ZENITH ANGLE
  REAL (KIND=8) pto1(kdlon) ! OPTICAL THICKNESS
  REAL (KIND=8) pw(kdlon) ! SINGLE SCATTERING ALBEDO
  REAL (KIND=8) pre1(kdlon) ! LAYER REFLECTIVITY (NO UNDERLYING-LAYER REFLECTION)
  REAL (KIND=8) pre2(kdlon) ! LAYER REFLECTIVITY
  REAL (KIND=8) ptr1(kdlon) ! LAYER TRANSMISSIVITY (NO UNDERLYING-LAYER REFLECTION)
  REAL (KIND=8) ptr2(kdlon) ! LAYER TRANSMISSIVITY

  ! * LOCAL VARIABLES:

  INTEGER jl
  REAL (KIND=8) zff, zgp, ztop, zwcp, zdt, zx1, zwm
  REAL (KIND=8) zrm2, zrk, zx2, zrp, zalpha, zbeta, zarg
  REAL (KIND=8) zexmu0, zarg2, zexkp, zexkm, zxp2p, zxm2p, zap2b, zam2b
  REAL (KIND=8) za11, za12, za13, za21, za22, za23
  REAL (KIND=8) zdena, zc1a, zc2a, zri0a, zri1a
  REAL (KIND=8) zri0b, zri1b
  REAL (KIND=8) zb21, zb22, zb23, zdenb, zc1b, zc2b
  REAL (KIND=8) zri0c, zri1c, zri0d, zri1d

  ! ------------------------------------------------------------------

  ! *         1.      DELTA-EDDINGTON CALCULATIONS


  DO jl = 1, kdlon
    ! *         1.1     SET UP THE DELTA-MODIFIED PARAMETERS


    zff = pgg(jl)*pgg(jl)
    zgp = pgg(jl)/(1.+pgg(jl))
    ztop = (1.-pw(jl)*zff)*pto1(jl)
    zwcp = (1-zff)*pw(jl)/(1.-pw(jl)*zff)
    zdt = 2./3.
    zx1 = 1. - zwcp*zgp
    zwm = 1. - zwcp
    zrm2 = prmuz(jl)*prmuz(jl)
    zrk = sqrt(3.*zwm*zx1)
    zx2 = 4.*(1.-zrk*zrk*zrm2)
    zrp = zrk/zx1
    zalpha = 3.*zwcp*zrm2*(1.+zgp*zwm)/zx2
    zbeta = 3.*zwcp*prmuz(jl)*(1.+3.*zgp*zrm2*zwm)/zx2
    zarg = min(ztop/prmuz(jl), 200._8)
    zexmu0 = exp(-zarg)
    zarg2 = min(zrk*ztop, 200._8)
    zexkp = exp(zarg2)
    zexkm = 1./zexkp
    zxp2p = 1. + zdt*zrp
    zxm2p = 1. - zdt*zrp
    zap2b = zalpha + zdt*zbeta
    zam2b = zalpha - zdt*zbeta

    ! *         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER


    za11 = zxp2p
    za12 = zxm2p
    za13 = zap2b
    za22 = zxp2p*zexkp
    za21 = zxm2p*zexkm
    za23 = zam2b*zexmu0
    zdena = za11*za22 - za21*za12
    zc1a = (za22*za13-za12*za23)/zdena
    zc2a = (za11*za23-za21*za13)/zdena
    zri0a = zc1a + zc2a - zalpha
    zri1a = zrp*(zc1a-zc2a) - zbeta
    pre1(jl) = (zri0a-zdt*zri1a)/prmuz(jl)
    zri0b = zc1a*zexkm + zc2a*zexkp - zalpha*zexmu0
    zri1b = zrp*(zc1a*zexkm-zc2a*zexkp) - zbeta*zexmu0
    ptr1(jl) = zexmu0 + (zri0b+zdt*zri1b)/prmuz(jl)

    ! *         1.3     WITH REFLECTION FROM THE UNDERLYING LAYER


    zb21 = za21 - pref(jl)*zxp2p*zexkm
    zb22 = za22 - pref(jl)*zxm2p*zexkp
    zb23 = za23 - pref(jl)*zexmu0*(zap2b-prmuz(jl))
    zdenb = za11*zb22 - zb21*za12
    zc1b = (zb22*za13-za12*zb23)/zdenb
    zc2b = (za11*zb23-zb21*za13)/zdenb
    zri0c = zc1b + zc2b - zalpha
    zri1c = zrp*(zc1b-zc2b) - zbeta
    pre2(jl) = (zri0c-zdt*zri1c)/prmuz(jl)
    zri0d = zc1b*zexkm + zc2b*zexkp - zalpha*zexmu0
    zri1d = zrp*(zc1b*zexkm-zc2b*zexkp) - zbeta*zexmu0
    ptr2(jl) = zexmu0 + (zri0d+zdt*zri1d)/prmuz(jl)

  END DO
  RETURN
END SUBROUTINE swde_lmdar4
SUBROUTINE swtt_lmdar4(knu, ka, pu, ptr)
  USE dimphy
  USE radiation_ar4_param, ONLY: apad, bpad, d
  IMPLICIT NONE

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
  ! ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
  ! INTERVALS.

  ! METHOD.
  ! -------

  ! TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
  ! AND HORNER'S ALGORITHM.

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 88-12-15
  ! -----------------------------------------------------------------------

  ! * ARGUMENTS

  INTEGER knu ! INDEX OF THE SPECTRAL INTERVAL
  INTEGER ka ! INDEX OF THE ABSORBER
  REAL (KIND=8) pu(kdlon) ! ABSORBER AMOUNT

  REAL (KIND=8) ptr(kdlon) ! TRANSMISSION FUNCTION

  ! * LOCAL VARIABLES:

  REAL (KIND=8) zr1(kdlon), zr2(kdlon)
  INTEGER jl, i, j

  ! -----------------------------------------------------------------------

  ! *         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION


  DO jl = 1, kdlon
    zr1(jl) = apad(knu, ka, 1) + pu(jl)*(apad(knu,ka,2)+pu(jl)*(apad(knu,ka, &
      3)+pu(jl)*(apad(knu,ka,4)+pu(jl)*(apad(knu,ka,5)+pu(jl)*(apad(knu,ka,6) &
      +pu(jl)*(apad(knu,ka,7)))))))

    zr2(jl) = bpad(knu, ka, 1) + pu(jl)*(bpad(knu,ka,2)+pu(jl)*(bpad(knu,ka, &
      3)+pu(jl)*(bpad(knu,ka,4)+pu(jl)*(bpad(knu,ka,5)+pu(jl)*(bpad(knu,ka,6) &
      +pu(jl)*(bpad(knu,ka,7)))))))

    ! *         2.      ADD THE BACKGROUND TRANSMISSION



    ptr(jl) = (zr1(jl)/zr2(jl))*(1.-d(knu,ka)) + d(knu, ka)
  END DO

  RETURN
END SUBROUTINE swtt_lmdar4
SUBROUTINE swtt1_lmdar4(knu, kabs, kind, pu, ptr)
  USE dimphy
  USE radiation_ar4_param, ONLY: apad, bpad, d
  IMPLICIT NONE

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
  ! ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
  ! INTERVALS.

  ! METHOD.
  ! -------

  ! TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
  ! AND HORNER'S ALGORITHM.

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 95-01-20
  ! -----------------------------------------------------------------------
  ! * ARGUMENTS:

  INTEGER knu ! INDEX OF THE SPECTRAL INTERVAL
  INTEGER kabs ! NUMBER OF ABSORBERS
  INTEGER kind(kabs) ! INDICES OF THE ABSORBERS
  REAL (KIND=8) pu(kdlon, kabs) ! ABSORBER AMOUNT

  REAL (KIND=8) ptr(kdlon, kabs) ! TRANSMISSION FUNCTION

  ! * LOCAL VARIABLES:

  REAL (KIND=8) zr1(kdlon)
  REAL (KIND=8) zr2(kdlon)
  REAL (KIND=8) zu(kdlon)
  INTEGER jl, ja, i, j, ia

  ! -----------------------------------------------------------------------

  ! *         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION


  DO ja = 1, kabs
    ia = kind(ja)
    DO jl = 1, kdlon
      zu(jl) = pu(jl, ja)
      zr1(jl) = apad(knu, ia, 1) + zu(jl)*(apad(knu,ia,2)+zu(jl)*(apad(knu, &
        ia,3)+zu(jl)*(apad(knu,ia,4)+zu(jl)*(apad(knu,ia,5)+zu(jl)*(apad(knu, &
        ia,6)+zu(jl)*(apad(knu,ia,7)))))))

      zr2(jl) = bpad(knu, ia, 1) + zu(jl)*(bpad(knu,ia,2)+zu(jl)*(bpad(knu, &
        ia,3)+zu(jl)*(bpad(knu,ia,4)+zu(jl)*(bpad(knu,ia,5)+zu(jl)*(bpad(knu, &
        ia,6)+zu(jl)*(bpad(knu,ia,7)))))))

      ! *         2.      ADD THE BACKGROUND TRANSMISSION


      ptr(jl, ja) = (zr1(jl)/zr2(jl))*(1.-d(knu,ia)) + d(knu, ia)
    END DO
  END DO

  RETURN
END SUBROUTINE swtt1_lmdar4
! IM ctes ds clesphys.h   SUBROUTINE LW(RCO2,RCH4,RN2O,RCFC11,RCFC12,
SUBROUTINE lw_lmdar4(ppmb, pdp, ppsol, pdt0, pemis, ptl, ptave, pwv, pozon, &
    paer, pcldld, pcldlu, pview, pcolr, pcolr0, ptoplw, psollw, ptoplw0, &
    psollw0, psollwdown, &         ! IM  .
                                   ! psollwdown,psollwdownclr,
  ! IM  .              ptoplwdown,ptoplwdownclr)
    plwup, plwdn, plwup0, plwdn0)
  USE dimphy
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE
  include "raddimlw.h"
  include "YOMCST.h"

  ! -----------------------------------------------------------------------
  ! METHOD.
  ! -------

  ! 1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
  ! ABSORBERS.
  ! 2. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
  ! GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
  ! 3. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
  ! TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
  ! BOUNDARIES.
  ! 4. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
  ! 5. INTRODUCES THE EFFECTS OF THE CLOUDS ON THE FLUXES.


  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! -----------------------------------------------------------------------
  ! IM ctes ds clesphys.h
  ! REAL(KIND=8) RCO2   ! CO2 CONCENTRATION (IPCC:353.E-06* 44.011/28.97)
  ! REAL(KIND=8) RCH4   ! CH4 CONCENTRATION (IPCC: 1.72E-06* 16.043/28.97)
  ! REAL(KIND=8) RN2O   ! N2O CONCENTRATION (IPCC: 310.E-09* 44.013/28.97)
  ! REAL(KIND=8) RCFC11 ! CFC11 CONCENTRATION (IPCC: 280.E-12*
  ! 137.3686/28.97)
  ! REAL(KIND=8) RCFC12 ! CFC12 CONCENTRATION (IPCC: 484.E-12*
  ! 120.9140/28.97)
  include "clesphys.h"
  REAL (KIND=8) pcldld(kdlon, kflev) ! DOWNWARD EFFECTIVE CLOUD COVER
  REAL (KIND=8) pcldlu(kdlon, kflev) ! UPWARD EFFECTIVE CLOUD COVER
  REAL (KIND=8) pdp(kdlon, kflev) ! LAYER PRESSURE THICKNESS (Pa)
  REAL (KIND=8) pdt0(kdlon) ! SURFACE TEMPERATURE DISCONTINUITY (K)
  REAL (KIND=8) pemis(kdlon) ! SURFACE EMISSIVITY
  REAL (KIND=8) ppmb(kdlon, kflev+1) ! HALF LEVEL PRESSURE (mb)
  REAL (KIND=8) ppsol(kdlon) ! SURFACE PRESSURE (Pa)
  REAL (KIND=8) pozon(kdlon, kflev) ! O3 mass fraction
  REAL (KIND=8) ptl(kdlon, kflev+1) ! HALF LEVEL TEMPERATURE (K)
  REAL (KIND=8) paer(kdlon, kflev, 5) ! OPTICAL THICKNESS OF THE AEROSOLS
  REAL (KIND=8) ptave(kdlon, kflev) ! LAYER TEMPERATURE (K)
  REAL (KIND=8) pview(kdlon) ! COSECANT OF VIEWING ANGLE
  REAL (KIND=8) pwv(kdlon, kflev) ! SPECIFIC HUMIDITY (kg/kg)

  REAL (KIND=8) pcolr(kdlon, kflev) ! LONG-WAVE TENDENCY (K/day)
  REAL (KIND=8) pcolr0(kdlon, kflev) ! LONG-WAVE TENDENCY (K/day) clear-sky
  REAL (KIND=8) ptoplw(kdlon) ! LONGWAVE FLUX AT T.O.A.
  REAL (KIND=8) psollw(kdlon) ! LONGWAVE FLUX AT SURFACE
  REAL (KIND=8) ptoplw0(kdlon) ! LONGWAVE FLUX AT T.O.A. (CLEAR-SKY)
  REAL (KIND=8) psollw0(kdlon) ! LONGWAVE FLUX AT SURFACE (CLEAR-SKY)
  ! Rajout LF
  REAL (KIND=8) psollwdown(kdlon) ! LONGWAVE downwards flux at surface
  ! Rajout IM
  ! IM   real(kind=8) psollwdownclr(kdlon) ! LONGWAVE CS downwards flux at
  ! surface
  ! IM   real(kind=8) ptoplwdown(kdlon)    ! LONGWAVE downwards flux at
  ! T.O.A.
  ! IM   real(kind=8) ptoplwdownclr(kdlon) ! LONGWAVE CS downwards flux at
  ! T.O.A.
  ! IM
  REAL (KIND=8) plwup(kdlon, kflev+1) ! LW up total sky
  REAL (KIND=8) plwup0(kdlon, kflev+1) ! LW up clear sky
  REAL (KIND=8) plwdn(kdlon, kflev+1) ! LW down total sky
  REAL (KIND=8) plwdn0(kdlon, kflev+1) ! LW down clear sky
  ! -------------------------------------------------------------------------
  REAL (KIND=8) zabcu(kdlon, nua, 3*kflev+1)

  REAL (KIND=8) zoz(kdlon, kflev)
  ! equivalent pressure of ozone in a layer, in Pa

  ! ym      REAL(KIND=8) ZFLUX(KDLON,2,KFLEV+1) ! RADIATIVE FLUXES (1:up;
  ! 2:down)
  ! ym      REAL(KIND=8) ZFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
  ! ym      REAL(KIND=8) ZBINT(KDLON,KFLEV+1)            ! Intermediate
  ! variable
  ! ym      REAL(KIND=8) ZBSUI(KDLON)                    ! Intermediate
  ! variable
  ! ym      REAL(KIND=8) ZCTS(KDLON,KFLEV)               ! Intermediate
  ! variable
  ! ym      REAL(KIND=8) ZCNTRB(KDLON,KFLEV+1,KFLEV+1)   ! Intermediate
  ! variable
  ! ym      SAVE ZFLUX, ZFLUC, ZBINT, ZBSUI, ZCTS, ZCNTRB
  REAL (KIND=8), ALLOCATABLE, SAVE :: zflux(:, :, :) ! RADIATIVE FLUXES (1:up; 2:down)
  REAL (KIND=8), ALLOCATABLE, SAVE :: zfluc(:, :, :) ! CLEAR-SKY RADIATIVE FLUXES
  REAL (KIND=8), ALLOCATABLE, SAVE :: zbint(:, :) ! Intermediate variable
  REAL (KIND=8), ALLOCATABLE, SAVE :: zbsui(:) ! Intermediate variable
  REAL (KIND=8), ALLOCATABLE, SAVE :: zcts(:, :) ! Intermediate variable
  REAL (KIND=8), ALLOCATABLE, SAVE :: zcntrb(:, :, :) ! Intermediate variable
  !$OMP THREADPRIVATE(ZFLUX, ZFLUC, ZBINT, ZBSUI, ZCTS, ZCNTRB)

  INTEGER ilim, i, k, kpl1

  INTEGER lw0pas ! Every lw0pas steps, clear-sky is done
  PARAMETER (lw0pas=1)
  INTEGER lwpas ! Every lwpas steps, cloudy-sky is done
  PARAMETER (lwpas=1)

  INTEGER itaplw0, itaplw
  LOGICAL appel1er
  SAVE appel1er, itaplw0, itaplw
  !$OMP THREADPRIVATE(appel1er, itaplw0, itaplw)
  DATA appel1er/.TRUE./
  DATA itaplw0, itaplw/0, 0/

  ! ------------------------------------------------------------------
  IF (appel1er) THEN
    WRITE (lunout, *) 'LW clear-sky calling frequency: ', lw0pas
    WRITE (lunout, *) 'LW cloudy-sky calling frequency: ', lwpas
    WRITE (lunout, *) '   In general, they should be 1'
    ! ym
    ALLOCATE (zflux(kdlon,2,kflev+1))
    ALLOCATE (zfluc(kdlon,2,kflev+1))
    ALLOCATE (zbint(kdlon,kflev+1))
    ALLOCATE (zbsui(kdlon))
    ALLOCATE (zcts(kdlon,kflev))
    ALLOCATE (zcntrb(kdlon,kflev+1,kflev+1))
    appel1er = .FALSE.
  END IF

  IF (mod(itaplw0,lw0pas)==0) THEN
    ! Compute equivalent pressure of ozone from mass fraction:
    DO k = 1, kflev
      DO i = 1, kdlon
        zoz(i, k) = pozon(i, k)*pdp(i, k)
      END DO
    END DO
    ! IM ctes ds clesphys.h   CALL LWU(RCO2,RCH4, RN2O, RCFC11, RCFC12,
    CALL lwu_lmdar4(paer, pdp, ppmb, ppsol, zoz, ptave, pview, pwv, zabcu)
    CALL lwbv_lmdar4(ilim, pdp, pdt0, pemis, ppmb, ptl, ptave, zabcu, zfluc, &
      zbint, zbsui, zcts, zcntrb)
    itaplw0 = 0
  END IF
  itaplw0 = itaplw0 + 1

  IF (mod(itaplw,lwpas)==0) THEN
    CALL lwc_lmdar4(ilim, pcldld, pcldlu, pemis, zfluc, zbint, zbsui, zcts, &
      zcntrb, zflux)
    itaplw = 0
  END IF
  itaplw = itaplw + 1

  DO k = 1, kflev
    kpl1 = k + 1
    DO i = 1, kdlon
      pcolr(i, k) = zflux(i, 1, kpl1) + zflux(i, 2, kpl1) - zflux(i, 1, k) - &
        zflux(i, 2, k)
      pcolr(i, k) = pcolr(i, k)*rday*rg/rcpd/pdp(i, k)
      pcolr0(i, k) = zfluc(i, 1, kpl1) + zfluc(i, 2, kpl1) - zfluc(i, 1, k) - &
        zfluc(i, 2, k)
      pcolr0(i, k) = pcolr0(i, k)*rday*rg/rcpd/pdp(i, k)
    END DO
  END DO
  DO i = 1, kdlon
    psollw(i) = -zflux(i, 1, 1) - zflux(i, 2, 1)
    ptoplw(i) = zflux(i, 1, kflev+1) + zflux(i, 2, kflev+1)

    psollw0(i) = -zfluc(i, 1, 1) - zfluc(i, 2, 1)
    ptoplw0(i) = zfluc(i, 1, kflev+1) + zfluc(i, 2, kflev+1)
    psollwdown(i) = -zflux(i, 2, 1)

    ! IM attention aux signes !; LWtop >0, LWdn < 0
    DO k = 1, kflev + 1
      plwup(i, k) = zflux(i, 1, k)
      plwup0(i, k) = zfluc(i, 1, k)
      plwdn(i, k) = zflux(i, 2, k)
      plwdn0(i, k) = zfluc(i, 2, k)
    END DO
  END DO
  ! ------------------------------------------------------------------
  RETURN
END SUBROUTINE lw_lmdar4
! IM ctes ds clesphys.h   SUBROUTINE LWU(RCO2, RCH4, RN2O, RCFC11, RCFC12,
SUBROUTINE lwu_lmdar4(paer, pdp, ppmb, ppsol, poz, ptave, pview, pwv, pabcu)
  USE dimphy
  USE radiation_ar4_param, ONLY: tref, rt1, raer, at, bt, oct
  USE infotrac_phy, ONLY: type_trac
#ifdef REPROBUS
  USE chem_rep, ONLY: rch42d, rn2o2d, rcfc112d, rcfc122d, ok_rtime2d
#endif

  IMPLICIT NONE
  include "raddimlw.h"
  include "YOMCST.h"
  include "radepsi.h"
  include "radopt.h"

  ! PURPOSE.
  ! --------
  ! COMPUTES ABSORBER AMOUNTS INCLUDING PRESSURE AND
  ! TEMPERATURE EFFECTS

  ! METHOD.
  ! -------

  ! 1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
  ! ABSORBERS.


  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! Voigt lines (loop 404 modified) - JJM & PhD - 01/96
  ! -----------------------------------------------------------------------
  ! * ARGUMENTS:
  ! IM ctes ds clesphys.h
  ! REAL(KIND=8) RCO2
  ! REAL(KIND=8) RCH4, RN2O, RCFC11, RCFC12
  include "clesphys.h"
  REAL (KIND=8) paer(kdlon, kflev, 5)
  REAL (KIND=8) pdp(kdlon, kflev)
  REAL (KIND=8) ppmb(kdlon, kflev+1)
  REAL (KIND=8) ppsol(kdlon)
  REAL (KIND=8) poz(kdlon, kflev)
  REAL (KIND=8) ptave(kdlon, kflev)
  REAL (KIND=8) pview(kdlon)
  REAL (KIND=8) pwv(kdlon, kflev)

  REAL (KIND=8) pabcu(kdlon, nua, 3*kflev+1) ! EFFECTIVE ABSORBER AMOUNTS

  ! -----------------------------------------------------------------------
  ! * LOCAL VARIABLES:
  REAL (KIND=8) zably(kdlon, nua, 3*kflev+1)
  REAL (KIND=8) zduc(kdlon, 3*kflev+1)
  REAL (KIND=8) zphio(kdlon)
  REAL (KIND=8) zpsc2(kdlon)
  REAL (KIND=8) zpsc3(kdlon)
  REAL (KIND=8) zpsh1(kdlon)
  REAL (KIND=8) zpsh2(kdlon)
  REAL (KIND=8) zpsh3(kdlon)
  REAL (KIND=8) zpsh4(kdlon)
  REAL (KIND=8) zpsh5(kdlon)
  REAL (KIND=8) zpsh6(kdlon)
  REAL (KIND=8) zpsio(kdlon)
  REAL (KIND=8) ztcon(kdlon)
  REAL (KIND=8) zphm6(kdlon)
  REAL (KIND=8) zpsm6(kdlon)
  REAL (KIND=8) zphn6(kdlon)
  REAL (KIND=8) zpsn6(kdlon)
  REAL (KIND=8) zssig(kdlon, 3*kflev+1)
  REAL (KIND=8) ztavi(kdlon)
  REAL (KIND=8) zuaer(kdlon, ninter)
  REAL (KIND=8) zxoz(kdlon)
  REAL (KIND=8) zxwv(kdlon)

  INTEGER jl, jk, jkj, jkjr, jkjp, ig1
  INTEGER jki, jkip1, ja, jj
  INTEGER jkl, jkp1, jkk, jkjpn
  INTEGER jae1, jae2, jae3, jae, jjpn
  INTEGER ir, jc, jcp1
  REAL (KIND=8) zdpm, zupm, zupmh2o, zupmco2, zupmo3, zu6, zup
  REAL (KIND=8) zfppw, ztx, ztx2, zzably
  REAL (KIND=8) zcah1, zcbh1, zcah2, zcbh2, zcah3, zcbh3
  REAL (KIND=8) zcah4, zcbh4, zcah5, zcbh5, zcah6, zcbh6
  REAL (KIND=8) zcac8, zcbc8
  REAL (KIND=8) zalup, zdiff

  REAL (KIND=8) pvgco2, pvgh2o, pvgo3

  REAL (KIND=8) r10e ! DECIMAL/NATURAL LOG.FACTOR
  PARAMETER (r10e=0.4342945)

  ! -----------------------------------------------------------------------

  IF (levoigt) THEN
    pvgco2 = 60.
    pvgh2o = 30.
    pvgo3 = 400.
  ELSE
    pvgco2 = 0.
    pvgh2o = 0.
    pvgo3 = 0.
  END IF

  ! *         2.    PRESSURE OVER GAUSS SUB-LEVELS
  ! ------------------------------


  DO jl = 1, kdlon
    zssig(jl, 1) = ppmb(jl, 1)*100.
  END DO

  DO jk = 1, kflev
    jkj = (jk-1)*ng1p1 + 1
    jkjr = jkj
    jkjp = jkj + ng1p1
    DO jl = 1, kdlon
      zssig(jl, jkjp) = ppmb(jl, jk+1)*100.
    END DO
    DO ig1 = 1, ng1
      jkj = jkj + 1
      DO jl = 1, kdlon
        zssig(jl, jkj) = (zssig(jl,jkjr)+zssig(jl,jkjp))*0.5 + &
          rt1(ig1)*(zssig(jl,jkjp)-zssig(jl,jkjr))*0.5
      END DO
    END DO
  END DO

  ! -----------------------------------------------------------------------


  ! *         4.    PRESSURE THICKNESS AND MEAN PRESSURE OF SUB-LAYERS
  ! --------------------------------------------------


  DO jki = 1, 3*kflev
    jkip1 = jki + 1
    DO jl = 1, kdlon
      zably(jl, 5, jki) = (zssig(jl,jki)+zssig(jl,jkip1))*0.5
      zably(jl, 3, jki) = (zssig(jl,jki)-zssig(jl,jkip1))/(10.*rg)
    END DO
  END DO

  DO jk = 1, kflev
    jkp1 = jk + 1
    jkl = kflev + 1 - jk
    DO jl = 1, kdlon
      zxwv(jl) = max(pwv(jl,jk), zepscq)
      zxoz(jl) = max(poz(jl,jk)/pdp(jl,jk), zepsco)
    END DO
    jkj = (jk-1)*ng1p1 + 1
    jkjpn = jkj + ng1
    DO jkk = jkj, jkjpn
      DO jl = 1, kdlon
        zdpm = zably(jl, 3, jkk)
        zupm = zably(jl, 5, jkk)*zdpm/101325.
        zupmco2 = (zably(jl,5,jkk)+pvgco2)*zdpm/101325.
        zupmh2o = (zably(jl,5,jkk)+pvgh2o)*zdpm/101325.
        zupmo3 = (zably(jl,5,jkk)+pvgo3)*zdpm/101325.
        zduc(jl, jkk) = zdpm
        zably(jl, 12, jkk) = zxoz(jl)*zdpm
        zably(jl, 13, jkk) = zxoz(jl)*zupmo3
        zu6 = zxwv(jl)*zupm
        zfppw = 1.6078*zxwv(jl)/(1.+0.608*zxwv(jl))
        zably(jl, 6, jkk) = zxwv(jl)*zupmh2o
        zably(jl, 11, jkk) = zu6*zfppw
        zably(jl, 10, jkk) = zu6*(1.-zfppw)
        zably(jl, 9, jkk) = rco2*zupmco2
        zably(jl, 8, jkk) = rco2*zdpm
      END DO
    END DO
  END DO

  ! -----------------------------------------------------------------------


  ! *         5.    CUMULATIVE ABSORBER AMOUNTS FROM TOP OF ATMOSPHERE
  ! --------------------------------------------------


  DO ja = 1, nua
    DO jl = 1, kdlon
      pabcu(jl, ja, 3*kflev+1) = 0.
    END DO
  END DO

  DO jk = 1, kflev
    jj = (jk-1)*ng1p1 + 1
    jjpn = jj + ng1
    jkl = kflev + 1 - jk

    ! *         5.1  CUMULATIVE AEROSOL AMOUNTS FROM TOP OF ATMOSPHERE
    ! --------------------------------------------------


    jae1 = 3*kflev + 1 - jj
    jae2 = 3*kflev + 1 - (jj+1)
    jae3 = 3*kflev + 1 - jjpn
    DO jae = 1, 5
      DO jl = 1, kdlon
        zuaer(jl, jae) = (raer(jae,1)*paer(jl,jkl,1)+raer(jae,2)*paer(jl,jkl, &
          2)+raer(jae,3)*paer(jl,jkl,3)+raer(jae,4)*paer(jl,jkl,4)+ &
          raer(jae,5)*paer(jl,jkl,5))/(zduc(jl,jae1)+zduc(jl,jae2)+zduc(jl, &
          jae3))
      END DO
    END DO

    ! *         5.2  INTRODUCES TEMPERATURE EFFECTS ON ABSORBER AMOUNTS
    ! --------------------------------------------------


    DO jl = 1, kdlon
      ztavi(jl) = ptave(jl, jkl)
      ztcon(jl) = exp(6.08*(296./ztavi(jl)-1.))
      ztx = ztavi(jl) - tref
      ztx2 = ztx*ztx
      zzably = zably(jl, 6, jae1) + zably(jl, 6, jae2) + zably(jl, 6, jae3)
      zup = min(max(0.5*r10e*log(zzably)+5.,0._8), 6._8)
      zcah1 = at(1, 1) + zup*(at(1,2)+zup*(at(1,3)))
      zcbh1 = bt(1, 1) + zup*(bt(1,2)+zup*(bt(1,3)))
      zpsh1(jl) = exp(zcah1*ztx+zcbh1*ztx2)
      zcah2 = at(2, 1) + zup*(at(2,2)+zup*(at(2,3)))
      zcbh2 = bt(2, 1) + zup*(bt(2,2)+zup*(bt(2,3)))
      zpsh2(jl) = exp(zcah2*ztx+zcbh2*ztx2)
      zcah3 = at(3, 1) + zup*(at(3,2)+zup*(at(3,3)))
      zcbh3 = bt(3, 1) + zup*(bt(3,2)+zup*(bt(3,3)))
      zpsh3(jl) = exp(zcah3*ztx+zcbh3*ztx2)
      zcah4 = at(4, 1) + zup*(at(4,2)+zup*(at(4,3)))
      zcbh4 = bt(4, 1) + zup*(bt(4,2)+zup*(bt(4,3)))
      zpsh4(jl) = exp(zcah4*ztx+zcbh4*ztx2)
      zcah5 = at(5, 1) + zup*(at(5,2)+zup*(at(5,3)))
      zcbh5 = bt(5, 1) + zup*(bt(5,2)+zup*(bt(5,3)))
      zpsh5(jl) = exp(zcah5*ztx+zcbh5*ztx2)
      zcah6 = at(6, 1) + zup*(at(6,2)+zup*(at(6,3)))
      zcbh6 = bt(6, 1) + zup*(bt(6,2)+zup*(bt(6,3)))
      zpsh6(jl) = exp(zcah6*ztx+zcbh6*ztx2)
      zphm6(jl) = exp(-5.81E-4*ztx-1.13E-6*ztx2)
      zpsm6(jl) = exp(-5.57E-4*ztx-3.30E-6*ztx2)
      zphn6(jl) = exp(-3.46E-5*ztx+2.05E-7*ztx2)
      zpsn6(jl) = exp(3.70E-3*ztx-2.30E-6*ztx2)
    END DO

    DO jl = 1, kdlon
      ztavi(jl) = ptave(jl, jkl)
      ztx = ztavi(jl) - tref
      ztx2 = ztx*ztx
      zzably = zably(jl, 9, jae1) + zably(jl, 9, jae2) + zably(jl, 9, jae3)
      zalup = r10e*log(zzably)
      zup = max(0._8, 5.0+0.5*zalup)
      zpsc2(jl) = (ztavi(jl)/tref)**zup
      zcac8 = at(8, 1) + zup*(at(8,2)+zup*(at(8,3)))
      zcbc8 = bt(8, 1) + zup*(bt(8,2)+zup*(bt(8,3)))
      zpsc3(jl) = exp(zcac8*ztx+zcbc8*ztx2)
      zphio(jl) = exp(oct(1)*ztx+oct(2)*ztx2)
      zpsio(jl) = exp(2.*(oct(3)*ztx+oct(4)*ztx2))
    END DO

    DO jkk = jj, jjpn
      jc = 3*kflev + 1 - jkk
      jcp1 = jc + 1
      DO jl = 1, kdlon
        zdiff = pview(jl)
        pabcu(jl, 10, jc) = pabcu(jl, 10, jcp1) + zably(jl, 10, jc)*zdiff
        pabcu(jl, 11, jc) = pabcu(jl, 11, jcp1) + zably(jl, 11, jc)*ztcon(jl) &
          *zdiff

        pabcu(jl, 12, jc) = pabcu(jl, 12, jcp1) + zably(jl, 12, jc)*zphio(jl) &
          *zdiff
        pabcu(jl, 13, jc) = pabcu(jl, 13, jcp1) + zably(jl, 13, jc)*zpsio(jl) &
          *zdiff

        pabcu(jl, 7, jc) = pabcu(jl, 7, jcp1) + zably(jl, 9, jc)*zpsc2(jl)* &
          zdiff
        pabcu(jl, 8, jc) = pabcu(jl, 8, jcp1) + zably(jl, 9, jc)*zpsc3(jl)* &
          zdiff
        pabcu(jl, 9, jc) = pabcu(jl, 9, jcp1) + zably(jl, 9, jc)*zpsc3(jl)* &
          zdiff

        pabcu(jl, 1, jc) = pabcu(jl, 1, jcp1) + zably(jl, 6, jc)*zpsh1(jl)* &
          zdiff
        pabcu(jl, 2, jc) = pabcu(jl, 2, jcp1) + zably(jl, 6, jc)*zpsh2(jl)* &
          zdiff
        pabcu(jl, 3, jc) = pabcu(jl, 3, jcp1) + zably(jl, 6, jc)*zpsh5(jl)* &
          zdiff
        pabcu(jl, 4, jc) = pabcu(jl, 4, jcp1) + zably(jl, 6, jc)*zpsh3(jl)* &
          zdiff
        pabcu(jl, 5, jc) = pabcu(jl, 5, jcp1) + zably(jl, 6, jc)*zpsh4(jl)* &
          zdiff
        pabcu(jl, 6, jc) = pabcu(jl, 6, jcp1) + zably(jl, 6, jc)*zpsh6(jl)* &
          zdiff

        pabcu(jl, 14, jc) = pabcu(jl, 14, jcp1) + zuaer(jl, 1)*zduc(jl, jc)* &
          zdiff
        pabcu(jl, 15, jc) = pabcu(jl, 15, jcp1) + zuaer(jl, 2)*zduc(jl, jc)* &
          zdiff
        pabcu(jl, 16, jc) = pabcu(jl, 16, jcp1) + zuaer(jl, 3)*zduc(jl, jc)* &
          zdiff
        pabcu(jl, 17, jc) = pabcu(jl, 17, jcp1) + zuaer(jl, 4)*zduc(jl, jc)* &
          zdiff
        pabcu(jl, 18, jc) = pabcu(jl, 18, jcp1) + zuaer(jl, 5)*zduc(jl, jc)* &
          zdiff



        IF (type_trac=='repr') THEN
#ifdef REPROBUS
          IF (ok_rtime2d) THEN
            pabcu(jl, 19, jc) = pabcu(jl, 19, jcp1) + &
              zably(jl, 8, jc)*rch42d(jl, jc)/rco2*zphm6(jl)*zdiff
            pabcu(jl, 20, jc) = pabcu(jl, 20, jcp1) + &
              zably(jl, 9, jc)*rch42d(jl, jc)/rco2*zpsm6(jl)*zdiff
            pabcu(jl, 21, jc) = pabcu(jl, 21, jcp1) + &
              zably(jl, 8, jc)*rn2o2d(jl, jc)/rco2*zphn6(jl)*zdiff
            pabcu(jl, 22, jc) = pabcu(jl, 22, jcp1) + &
              zably(jl, 9, jc)*rn2o2d(jl, jc)/rco2*zpsn6(jl)*zdiff

            pabcu(jl, 23, jc) = pabcu(jl, 23, jcp1) + &
              zably(jl, 8, jc)*rcfc112d(jl, jc)/rco2*zdiff
            pabcu(jl, 24, jc) = pabcu(jl, 24, jcp1) + &
              zably(jl, 8, jc)*rcfc122d(jl, jc)/rco2*zdiff
          ELSE
              ! Same calculation as for type_trac /= repr
            pabcu(jl, 19, jc) = pabcu(jl, 19, jcp1) + &
              zably(jl, 8, jc)*rch4/rco2*zphm6(jl)*zdiff
            pabcu(jl, 20, jc) = pabcu(jl, 20, jcp1) + &
              zably(jl, 9, jc)*rch4/rco2*zpsm6(jl)*zdiff
            pabcu(jl, 21, jc) = pabcu(jl, 21, jcp1) + &
              zably(jl, 8, jc)*rn2o/rco2*zphn6(jl)*zdiff
            pabcu(jl, 22, jc) = pabcu(jl, 22, jcp1) + &
              zably(jl, 9, jc)*rn2o/rco2*zpsn6(jl)*zdiff

            pabcu(jl, 23, jc) = pabcu(jl, 23, jcp1) + &
              zably(jl, 8, jc)*rcfc11/rco2*zdiff
            pabcu(jl, 24, jc) = pabcu(jl, 24, jcp1) + &
              zably(jl, 8, jc)*rcfc12/rco2*zdiff
          END IF
#endif
        ELSE
          pabcu(jl, 19, jc) = pabcu(jl, 19, jcp1) + &
            zably(jl, 8, jc)*rch4/rco2*zphm6(jl)*zdiff
          pabcu(jl, 20, jc) = pabcu(jl, 20, jcp1) + &
            zably(jl, 9, jc)*rch4/rco2*zpsm6(jl)*zdiff
          pabcu(jl, 21, jc) = pabcu(jl, 21, jcp1) + &
            zably(jl, 8, jc)*rn2o/rco2*zphn6(jl)*zdiff
          pabcu(jl, 22, jc) = pabcu(jl, 22, jcp1) + &
            zably(jl, 9, jc)*rn2o/rco2*zpsn6(jl)*zdiff

          pabcu(jl, 23, jc) = pabcu(jl, 23, jcp1) + &
            zably(jl, 8, jc)*rcfc11/rco2*zdiff
          pabcu(jl, 24, jc) = pabcu(jl, 24, jcp1) + &
            zably(jl, 8, jc)*rcfc12/rco2*zdiff
        END IF

      END DO
    END DO

  END DO


  RETURN
END SUBROUTINE lwu_lmdar4
SUBROUTINE lwbv_lmdar4(klim, pdp, pdt0, pemis, ppmb, ptl, ptave, pabcu, &
    pfluc, pbint, pbsui, pcts, pcntrb)
  USE dimphy
  IMPLICIT NONE
  include "raddimlw.h"
  include "YOMCST.h"

  ! PURPOSE.
  ! --------
  ! TO COMPUTE THE PLANCK FUNCTION AND PERFORM THE
  ! VERTICAL INTEGRATION. SPLIT OUT FROM LW FOR MEMORY
  ! SAVING

  ! METHOD.
  ! -------

  ! 1. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
  ! GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
  ! 2. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
  ! TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
  ! BOUNDARIES.
  ! 3. COMPUTES THE CLEAR-SKY COOLING RATES.

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! MODIFICATION : 93-10-15 M.HAMRUD (SPLIT OUT FROM LW TO SAVE
  ! MEMORY)
  ! -----------------------------------------------------------------------
  ! * ARGUMENTS:
  INTEGER klim

  REAL (KIND=8) pdp(kdlon, kflev)
  REAL (KIND=8) pdt0(kdlon)
  REAL (KIND=8) pemis(kdlon)
  REAL (KIND=8) ppmb(kdlon, kflev+1)
  REAL (KIND=8) ptl(kdlon, kflev+1)
  REAL (KIND=8) ptave(kdlon, kflev)

  REAL (KIND=8) pfluc(kdlon, 2, kflev+1)

  REAL (KIND=8) pabcu(kdlon, nua, 3*kflev+1)
  REAL (KIND=8) pbint(kdlon, kflev+1)
  REAL (KIND=8) pbsui(kdlon)
  REAL (KIND=8) pcts(kdlon, kflev)
  REAL (KIND=8) pcntrb(kdlon, kflev+1, kflev+1)

  ! -------------------------------------------------------------------------

  ! * LOCAL VARIABLES:
  REAL (KIND=8) zb(kdlon, ninter, kflev+1)
  REAL (KIND=8) zbsur(kdlon, ninter)
  REAL (KIND=8) zbtop(kdlon, ninter)
  REAL (KIND=8) zdbsl(kdlon, ninter, kflev*2)
  REAL (KIND=8) zga(kdlon, 8, 2, kflev)
  REAL (KIND=8) zgb(kdlon, 8, 2, kflev)
  REAL (KIND=8) zgasur(kdlon, 8, 2)
  REAL (KIND=8) zgbsur(kdlon, 8, 2)
  REAL (KIND=8) zgatop(kdlon, 8, 2)
  REAL (KIND=8) zgbtop(kdlon, 8, 2)

  INTEGER nuaer, ntraer
  ! ------------------------------------------------------------------
  ! * COMPUTES PLANCK FUNCTIONS:
  CALL lwb_lmdar4(pdt0, ptave, ptl, zb, pbint, pbsui, zbsur, zbtop, zdbsl, &
    zga, zgb, zgasur, zgbsur, zgatop, zgbtop)
  ! ------------------------------------------------------------------
  ! * PERFORMS THE VERTICAL INTEGRATION:
  nuaer = nua
  ntraer = ntra
  CALL lwv_lmdar4(nuaer, ntraer, klim, pabcu, zb, pbint, pbsui, zbsur, zbtop, &
    zdbsl, pemis, ppmb, ptave, zga, zgb, zgasur, zgbsur, zgatop, zgbtop, &
    pcntrb, pcts, pfluc)
  ! ------------------------------------------------------------------
  RETURN
END SUBROUTINE lwbv_lmdar4
SUBROUTINE lwc_lmdar4(klim, pcldld, pcldlu, pemis, pfluc, pbint, pbsuin, &
    pcts, pcntrb, pflux)
  USE dimphy
  IMPLICIT NONE
  include "radepsi.h"
  include "radopt.h"

  ! PURPOSE.
  ! --------
  ! INTRODUCES CLOUD EFFECTS ON LONGWAVE FLUXES OR
  ! RADIANCES

  ! EXPLICIT ARGUMENTS :
  ! --------------------
  ! ==== INPUTS ===
  ! PBINT  : (KDLON,0:KFLEV)     ; HALF LEVEL PLANCK FUNCTION
  ! PBSUIN : (KDLON)             ; SURFACE PLANCK FUNCTION
  ! PCLDLD : (KDLON,KFLEV)       ; DOWNWARD EFFECTIVE CLOUD FRACTION
  ! PCLDLU : (KDLON,KFLEV)       ; UPWARD EFFECTIVE CLOUD FRACTION
  ! PCNTRB : (KDLON,KFLEV+1,KFLEV+1); CLEAR-SKY ENERGY EXCHANGE
  ! PCTS   : (KDLON,KFLEV)       ; CLEAR-SKY LAYER COOLING-TO-SPACE
  ! PEMIS  : (KDLON)             ; SURFACE EMISSIVITY
  ! PFLUC
  ! ==== OUTPUTS ===
  ! PFLUX(KDLON,2,KFLEV)         ; RADIATIVE FLUXES :
  ! 1  ==>  UPWARD   FLUX TOTAL
  ! 2  ==>  DOWNWARD FLUX TOTAL

  ! METHOD.
  ! -------

  ! 1. INITIALIZES ALL FLUXES TO CLEAR-SKY VALUES
  ! 2. EFFECT OF ONE OVERCAST UNITY EMISSIVITY CLOUD LAYER
  ! 3. EFFECT OF SEMI-TRANSPARENT, PARTIAL OR MULTI-LAYERED
  ! CLOUDS

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! Voigt lines (loop 231 to 233)  - JJM & PhD - 01/96
  ! -----------------------------------------------------------------------
  ! * ARGUMENTS:
  INTEGER klim
  REAL (KIND=8) pfluc(kdlon, 2, kflev+1) ! CLEAR-SKY RADIATIVE FLUXES
  REAL (KIND=8) pbint(kdlon, kflev+1) ! HALF LEVEL PLANCK FUNCTION
  REAL (KIND=8) pbsuin(kdlon) ! SURFACE PLANCK FUNCTION
  REAL (KIND=8) pcntrb(kdlon, kflev+1, kflev+1) !CLEAR-SKY ENERGY EXCHANGE
  REAL (KIND=8) pcts(kdlon, kflev) ! CLEAR-SKY LAYER COOLING-TO-SPACE

  REAL (KIND=8) pcldld(kdlon, kflev)
  REAL (KIND=8) pcldlu(kdlon, kflev)
  REAL (KIND=8) pemis(kdlon)

  REAL (KIND=8) pflux(kdlon, 2, kflev+1)
  ! -----------------------------------------------------------------------
  ! * LOCAL VARIABLES:
  INTEGER imx(kdlon), imxp(kdlon)

  REAL (KIND=8) zclear(kdlon), zcloud(kdlon), zdnf(kdlon, kflev+1, kflev+1), &
    zfd(kdlon), zfn10(kdlon), zfu(kdlon), zupf(kdlon, kflev+1, kflev+1)
  REAL (KIND=8) zclm(kdlon, kflev+1, kflev+1)

  INTEGER jk, jl, imaxc, imx1, imx2, jkj, jkp1, jkm1
  INTEGER jk1, jk2, jkc, jkcp1, jcloud
  INTEGER imxm1, imxp1
  REAL (KIND=8) zcfrac

  ! ------------------------------------------------------------------

  ! *         1.     INITIALIZATION
  ! --------------


  imaxc = 0

  DO jl = 1, kdlon
    imx(jl) = 0
    imxp(jl) = 0
    zcloud(jl) = 0.
  END DO

  ! *         1.1    SEARCH THE LAYER INDEX OF THE HIGHEST CLOUD
  ! -------------------------------------------


  DO jk = 1, kflev
    DO jl = 1, kdlon
      imx1 = imx(jl)
      imx2 = jk
      IF (pcldlu(jl,jk)>zepsc) THEN
        imxp(jl) = imx2
      ELSE
        imxp(jl) = imx1
      END IF
      imaxc = max(imxp(jl), imaxc)
      imx(jl) = imxp(jl)
    END DO
  END DO
  ! GM*******
  imaxc = kflev
  ! GM*******

  DO jk = 1, kflev + 1
    DO jl = 1, kdlon
      pflux(jl, 1, jk) = pfluc(jl, 1, jk)
      pflux(jl, 2, jk) = pfluc(jl, 2, jk)
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.      EFFECT OF CLOUDINESS ON LONGWAVE FLUXES
  ! ---------------------------------------

  IF (imaxc>0) THEN

    imxp1 = imaxc + 1
    imxm1 = imaxc - 1

    ! *         2.0     INITIALIZE TO CLEAR-SKY FLUXES
    ! ------------------------------


    DO jk1 = 1, kflev + 1
      DO jk2 = 1, kflev + 1
        DO jl = 1, kdlon
          zupf(jl, jk2, jk1) = pfluc(jl, 1, jk1)
          zdnf(jl, jk2, jk1) = pfluc(jl, 2, jk1)
        END DO
      END DO
    END DO

    ! *         2.1     FLUXES FOR ONE OVERCAST UNITY EMISSIVITY CLOUD
    ! ----------------------------------------------


    DO jkc = 1, imaxc
      jcloud = jkc
      jkcp1 = jcloud + 1

      ! *         2.1.1   ABOVE THE CLOUD
      ! ---------------


      DO jk = jkcp1, kflev + 1
        jkm1 = jk - 1
        DO jl = 1, kdlon
          zfu(jl) = 0.
        END DO
        IF (jk>jkcp1) THEN
          DO jkj = jkcp1, jkm1
            DO jl = 1, kdlon
              zfu(jl) = zfu(jl) + pcntrb(jl, jk, jkj)
            END DO
          END DO
        END IF

        DO jl = 1, kdlon
          zupf(jl, jkcp1, jk) = pbint(jl, jk) - zfu(jl)
        END DO
      END DO

      ! *         2.1.2   BELOW THE CLOUD
      ! ---------------


      DO jk = 1, jcloud
        jkp1 = jk + 1
        DO jl = 1, kdlon
          zfd(jl) = 0.
        END DO

        IF (jk<jcloud) THEN
          DO jkj = jkp1, jcloud
            DO jl = 1, kdlon
              zfd(jl) = zfd(jl) + pcntrb(jl, jk, jkj)
            END DO
          END DO
        END IF
        DO jl = 1, kdlon
          zdnf(jl, jkcp1, jk) = -pbint(jl, jk) - zfd(jl)
        END DO
      END DO

    END DO

    ! *         2.2     CLOUD COVER MATRIX
    ! ------------------

    ! *    ZCLM(JK1,JK2) IS THE OBSCURATION FACTOR BY CLOUD LAYERS BETWEEN
    ! HALF-LEVELS JK1 AND JK2 AS SEEN FROM JK1


    DO jk1 = 1, kflev + 1
      DO jk2 = 1, kflev + 1
        DO jl = 1, kdlon
          zclm(jl, jk1, jk2) = 0.
        END DO
      END DO
    END DO

    ! *         2.4     CLOUD COVER BELOW THE LEVEL OF CALCULATION
    ! ------------------------------------------


    DO jk1 = 2, kflev + 1
      DO jl = 1, kdlon
        zclear(jl) = 1.
        zcloud(jl) = 0.
      END DO
      DO jk = jk1 - 1, 1, -1
        DO jl = 1, kdlon
          IF (novlp==1) THEN
            ! * maximum-random
            zclear(jl) = zclear(jl)*(1.0-max(pcldlu(jl, &
              jk),zcloud(jl)))/(1.0-min(zcloud(jl),1.-zepsec))
            zclm(jl, jk1, jk) = 1.0 - zclear(jl)
            zcloud(jl) = pcldlu(jl, jk)
          ELSE IF (novlp==2) THEN
            ! * maximum
            zcloud(jl) = max(zcloud(jl), pcldlu(jl,jk))
            zclm(jl, jk1, jk) = zcloud(jl)
          ELSE IF (novlp==3) THEN
            ! * random
            zclear(jl) = zclear(jl)*(1.0-pcldlu(jl,jk))
            zcloud(jl) = 1.0 - zclear(jl)
            zclm(jl, jk1, jk) = zcloud(jl)
          END IF
        END DO
      END DO
    END DO

    ! *         2.5     CLOUD COVER ABOVE THE LEVEL OF CALCULATION
    ! ------------------------------------------


    DO jk1 = 1, kflev
      DO jl = 1, kdlon
        zclear(jl) = 1.
        zcloud(jl) = 0.
      END DO
      DO jk = jk1, kflev
        DO jl = 1, kdlon
          IF (novlp==1) THEN
            ! * maximum-random
            zclear(jl) = zclear(jl)*(1.0-max(pcldld(jl, &
              jk),zcloud(jl)))/(1.0-min(zcloud(jl),1.-zepsec))
            zclm(jl, jk1, jk) = 1.0 - zclear(jl)
            zcloud(jl) = pcldld(jl, jk)
          ELSE IF (novlp==2) THEN
            ! * maximum
            zcloud(jl) = max(zcloud(jl), pcldld(jl,jk))
            zclm(jl, jk1, jk) = zcloud(jl)
          ELSE IF (novlp==3) THEN
            ! * random
            zclear(jl) = zclear(jl)*(1.0-pcldld(jl,jk))
            zcloud(jl) = 1.0 - zclear(jl)
            zclm(jl, jk1, jk) = zcloud(jl)
          END IF
        END DO
      END DO
    END DO

    ! *         3.      FLUXES FOR PARTIAL/MULTIPLE LAYERED CLOUDINESS
    ! ----------------------------------------------


    ! *         3.1     DOWNWARD FLUXES
    ! ---------------


    DO jl = 1, kdlon
      pflux(jl, 2, kflev+1) = 0.
    END DO

    DO jk1 = kflev, 1, -1

      ! *                 CONTRIBUTION FROM CLEAR-SKY FRACTION

      DO jl = 1, kdlon
        zfd(jl) = (1.-zclm(jl,jk1,kflev))*zdnf(jl, 1, jk1)
      END DO

      ! *                 CONTRIBUTION FROM ADJACENT CLOUD

      DO jl = 1, kdlon
        zfd(jl) = zfd(jl) + zclm(jl, jk1, jk1)*zdnf(jl, jk1+1, jk1)
      END DO

      ! *                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS

      DO jk = kflev - 1, jk1, -1
        DO jl = 1, kdlon
          zcfrac = zclm(jl, jk1, jk+1) - zclm(jl, jk1, jk)
          zfd(jl) = zfd(jl) + zcfrac*zdnf(jl, jk+2, jk1)
        END DO
      END DO

      DO jl = 1, kdlon
        pflux(jl, 2, jk1) = zfd(jl)
      END DO

    END DO

    ! *         3.2     UPWARD FLUX AT THE SURFACE
    ! --------------------------


    DO jl = 1, kdlon
      pflux(jl, 1, 1) = pemis(jl)*pbsuin(jl) - (1.-pemis(jl))*pflux(jl, 2, 1)
    END DO

    ! *         3.3     UPWARD FLUXES
    ! -------------


    DO jk1 = 2, kflev + 1

      ! *                 CONTRIBUTION FROM CLEAR-SKY FRACTION

      DO jl = 1, kdlon
        zfu(jl) = (1.-zclm(jl,jk1,1))*zupf(jl, 1, jk1)
      END DO

      ! *                 CONTRIBUTION FROM ADJACENT CLOUD

      DO jl = 1, kdlon
        zfu(jl) = zfu(jl) + zclm(jl, jk1, jk1-1)*zupf(jl, jk1, jk1)
      END DO

      ! *                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS

      DO jk = 2, jk1 - 1
        DO jl = 1, kdlon
          zcfrac = zclm(jl, jk1, jk-1) - zclm(jl, jk1, jk)
          zfu(jl) = zfu(jl) + zcfrac*zupf(jl, jk, jk1)
        END DO
      END DO

      DO jl = 1, kdlon
        pflux(jl, 1, jk1) = zfu(jl)
      END DO

    END DO


  END IF

  ! *         2.3     END OF CLOUD EFFECT COMPUTATIONS


  IF (.NOT. levoigt) THEN
    DO jl = 1, kdlon
      zfn10(jl) = pflux(jl, 1, klim) + pflux(jl, 2, klim)
    END DO
    DO jk = klim + 1, kflev + 1
      DO jl = 1, kdlon
        zfn10(jl) = zfn10(jl) + pcts(jl, jk-1)
        pflux(jl, 1, jk) = zfn10(jl)
        pflux(jl, 2, jk) = 0.0
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE lwc_lmdar4
SUBROUTINE lwb_lmdar4(pdt0, ptave, ptl, pb, pbint, pbsuin, pbsur, pbtop, &
    pdbsl, pga, pgb, pgasur, pgbsur, pgatop, pgbtop)
  USE dimphy
  USE radiation_ar4_param, ONLY: tintp, xp, ga, gb
  IMPLICIT NONE
  include "raddimlw.h"

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! COMPUTES PLANCK FUNCTIONS

  ! EXPLICIT ARGUMENTS :
  ! --------------------
  ! ==== INPUTS ===
  ! PDT0   : (KDLON)             ; SURFACE TEMPERATURE DISCONTINUITY
  ! PTAVE  : (KDLON,KFLEV)       ; TEMPERATURE
  ! PTL    : (KDLON,0:KFLEV)     ; HALF LEVEL TEMPERATURE
  ! ==== OUTPUTS ===
  ! PB     : (KDLON,Ninter,KFLEV+1); SPECTRAL HALF LEVEL PLANCK FUNCTION
  ! PBINT  : (KDLON,KFLEV+1)     ; HALF LEVEL PLANCK FUNCTION
  ! PBSUIN : (KDLON)             ; SURFACE PLANCK FUNCTION
  ! PBSUR  : (KDLON,Ninter)        ; SURFACE SPECTRAL PLANCK FUNCTION
  ! PBTOP  : (KDLON,Ninter)        ; TOP SPECTRAL PLANCK FUNCTION
  ! PDBSL  : (KDLON,Ninter,KFLEV*2); SUB-LAYER PLANCK FUNCTION GRADIENT
  ! PGA    : (KDLON,8,2,KFLEV); dB/dT-weighted LAYER PADE APPROXIMANTS
  ! PGB    : (KDLON,8,2,KFLEV); dB/dT-weighted LAYER PADE APPROXIMANTS
  ! PGASUR, PGBSUR (KDLON,8,2)   ; SURFACE PADE APPROXIMANTS
  ! PGATOP, PGBTOP (KDLON,8,2)   ; T.O.A. PADE APPROXIMANTS

  ! IMPLICIT ARGUMENTS :   NONE
  ! --------------------

  ! METHOD.
  ! -------

  ! 1. COMPUTES THE PLANCK FUNCTION ON ALL LEVELS AND HALF LEVELS
  ! FROM A POLYNOMIAL DEVELOPMENT OF PLANCK FUNCTION

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS           "

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14

  ! -----------------------------------------------------------------------

  ! ARGUMENTS:

  REAL (KIND=8) pdt0(kdlon)
  REAL (KIND=8) ptave(kdlon, kflev)
  REAL (KIND=8) ptl(kdlon, kflev+1)

  REAL (KIND=8) pb(kdlon, ninter, kflev+1) ! SPECTRAL HALF LEVEL PLANCK FUNCTION
  REAL (KIND=8) pbint(kdlon, kflev+1) ! HALF LEVEL PLANCK FUNCTION
  REAL (KIND=8) pbsuin(kdlon) ! SURFACE PLANCK FUNCTION
  REAL (KIND=8) pbsur(kdlon, ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
  REAL (KIND=8) pbtop(kdlon, ninter) ! TOP SPECTRAL PLANCK FUNCTION
  REAL (KIND=8) pdbsl(kdlon, ninter, kflev*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
  REAL (KIND=8) pga(kdlon, 8, 2, kflev) ! dB/dT-weighted LAYER PADE APPROXIMANTS
  REAL (KIND=8) pgb(kdlon, 8, 2, kflev) ! dB/dT-weighted LAYER PADE APPROXIMANTS
  REAL (KIND=8) pgasur(kdlon, 8, 2) ! SURFACE PADE APPROXIMANTS
  REAL (KIND=8) pgbsur(kdlon, 8, 2) ! SURFACE PADE APPROXIMANTS
  REAL (KIND=8) pgatop(kdlon, 8, 2) ! T.O.A. PADE APPROXIMANTS
  REAL (KIND=8) pgbtop(kdlon, 8, 2) ! T.O.A. PADE APPROXIMANTS

  ! -------------------------------------------------------------------------
  ! *  LOCAL VARIABLES:
  INTEGER indb(kdlon), inds(kdlon)
  REAL (KIND=8) zblay(kdlon, kflev), zblev(kdlon, kflev+1)
  REAL (KIND=8) zres(kdlon), zres2(kdlon), zti(kdlon), zti2(kdlon)

  INTEGER jk, jl, ic, jnu, jf, jg
  INTEGER jk1, jk2
  INTEGER k, j, ixtox, indto, ixtx, indt
  INTEGER indsu, indtp
  REAL (KIND=8) zdsto1, zdstox, zdst1, zdstx

  ! * Quelques parametres:
  REAL (KIND=8) tstand
  PARAMETER (tstand=250.0)
  REAL (KIND=8) tstp
  PARAMETER (tstp=12.5)
  INTEGER mxixt
  PARAMETER (mxixt=10)

  ! * Used Data Block:
  ! REAL*8 TINTP(11)
  ! SAVE TINTP
  ! c$OMP THREADPRIVATE(TINTP)
  ! REAL*8 GA(11,16,3), GB(11,16,3)
  ! SAVE GA, GB
  ! c$OMP THREADPRIVATE(GA, GB)
  ! REAL*8 XP(6,6)
  ! SAVE XP
  ! c$OMP THREADPRIVATE(XP)

  ! DATA TINTP / 187.5, 200., 212.5, 225., 237.5, 250.,
  ! S             262.5, 275., 287.5, 300., 312.5 /
  ! -----------------------------------------------------------------------
  ! -- WATER VAPOR -- INT.1 -- 0- 500 CM-1 -- FROM ABS225 ----------------




  ! -- R.D. -- G = - 0.2 SLA


  ! ----- INTERVAL = 1 ----- T =  187.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 1, 1,IC),IC=1,3) /
  ! S 0.63499072E-02,-0.99506586E-03, 0.00000000E+00/
  ! DATA (GB( 1, 1,IC),IC=1,3) /
  ! S 0.63499072E-02, 0.97222852E-01, 0.10000000E+01/
  ! DATA (GA( 1, 2,IC),IC=1,3) /
  ! S 0.77266491E-02,-0.11661515E-02, 0.00000000E+00/
  ! DATA (GB( 1, 2,IC),IC=1,3) /
  ! S 0.77266491E-02, 0.10681591E+00, 0.10000000E+01/

  ! ----- INTERVAL = 1 ----- T =  200.0

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 2, 1,IC),IC=1,3) /
  ! S 0.65566348E-02,-0.10184169E-02, 0.00000000E+00/
  ! DATA (GB( 2, 1,IC),IC=1,3) /
  ! S 0.65566348E-02, 0.98862238E-01, 0.10000000E+01/
  ! DATA (GA( 2, 2,IC),IC=1,3) /
  ! S 0.81323287E-02,-0.11886130E-02, 0.00000000E+00/
  ! DATA (GB( 2, 2,IC),IC=1,3) /
  ! S 0.81323287E-02, 0.10921298E+00, 0.10000000E+01/

  ! ----- INTERVAL = 1 ----- T =  212.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 3, 1,IC),IC=1,3) /
  ! S 0.67849730E-02,-0.10404730E-02, 0.00000000E+00/
  ! DATA (GB( 3, 1,IC),IC=1,3) /
  ! S 0.67849730E-02, 0.10061504E+00, 0.10000000E+01/
  ! DATA (GA( 3, 2,IC),IC=1,3) /
  ! S 0.86507620E-02,-0.12139929E-02, 0.00000000E+00/
  ! DATA (GB( 3, 2,IC),IC=1,3) /
  ! S 0.86507620E-02, 0.11198225E+00, 0.10000000E+01/

  ! ----- INTERVAL = 1 ----- T =  225.0

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 4, 1,IC),IC=1,3) /
  ! S 0.70481947E-02,-0.10621792E-02, 0.00000000E+00/
  ! DATA (GB( 4, 1,IC),IC=1,3) /
  ! S 0.70481947E-02, 0.10256222E+00, 0.10000000E+01/
  ! DATA (GA( 4, 2,IC),IC=1,3) /
  ! S 0.92776391E-02,-0.12445811E-02, 0.00000000E+00/
  ! DATA (GB( 4, 2,IC),IC=1,3) /
  ! S 0.92776391E-02, 0.11487826E+00, 0.10000000E+01/

  ! ----- INTERVAL = 1 ----- T =  237.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 5, 1,IC),IC=1,3) /
  ! S 0.73585943E-02,-0.10847662E-02, 0.00000000E+00/
  ! DATA (GB( 5, 1,IC),IC=1,3) /
  ! S 0.73585943E-02, 0.10475952E+00, 0.10000000E+01/
  ! DATA (GA( 5, 2,IC),IC=1,3) /
  ! S 0.99806312E-02,-0.12807672E-02, 0.00000000E+00/
  ! DATA (GB( 5, 2,IC),IC=1,3) /
  ! S 0.99806312E-02, 0.11751113E+00, 0.10000000E+01/

  ! ----- INTERVAL = 1 ----- T =  250.0

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 6, 1,IC),IC=1,3) /
  ! S 0.77242818E-02,-0.11094726E-02, 0.00000000E+00/
  ! DATA (GB( 6, 1,IC),IC=1,3) /
  ! S 0.77242818E-02, 0.10720986E+00, 0.10000000E+01/
  ! DATA (GA( 6, 2,IC),IC=1,3) /
  ! S 0.10709803E-01,-0.13208251E-02, 0.00000000E+00/
  ! DATA (GB( 6, 2,IC),IC=1,3) /
  ! S 0.10709803E-01, 0.11951535E+00, 0.10000000E+01/

  ! ----- INTERVAL = 1 ----- T =  262.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 7, 1,IC),IC=1,3) /
  ! S 0.81472693E-02,-0.11372949E-02, 0.00000000E+00/
  ! DATA (GB( 7, 1,IC),IC=1,3) /
  ! S 0.81472693E-02, 0.10985370E+00, 0.10000000E+01/
  ! DATA (GA( 7, 2,IC),IC=1,3) /
  ! S 0.11414739E-01,-0.13619034E-02, 0.00000000E+00/
  ! DATA (GB( 7, 2,IC),IC=1,3) /
  ! S 0.11414739E-01, 0.12069945E+00, 0.10000000E+01/

  ! ----- INTERVAL = 1 ----- T =  275.0

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 8, 1,IC),IC=1,3) /
  ! S 0.86227527E-02,-0.11687683E-02, 0.00000000E+00/
  ! DATA (GB( 8, 1,IC),IC=1,3) /
  ! S 0.86227527E-02, 0.11257633E+00, 0.10000000E+01/
  ! DATA (GA( 8, 2,IC),IC=1,3) /
  ! S 0.12058772E-01,-0.14014165E-02, 0.00000000E+00/
  ! DATA (GB( 8, 2,IC),IC=1,3) /
  ! S 0.12058772E-01, 0.12108524E+00, 0.10000000E+01/

  ! ----- INTERVAL = 1 ----- T =  287.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 9, 1,IC),IC=1,3) /
  ! S 0.91396814E-02,-0.12038314E-02, 0.00000000E+00/
  ! DATA (GB( 9, 1,IC),IC=1,3) /
  ! S 0.91396814E-02, 0.11522980E+00, 0.10000000E+01/
  ! DATA (GA( 9, 2,IC),IC=1,3) /
  ! S 0.12623992E-01,-0.14378639E-02, 0.00000000E+00/
  ! DATA (GB( 9, 2,IC),IC=1,3) /
  ! S 0.12623992E-01, 0.12084229E+00, 0.10000000E+01/

  ! ----- INTERVAL = 1 ----- T =  300.0

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA(10, 1,IC),IC=1,3) /
  ! S 0.96825438E-02,-0.12418367E-02, 0.00000000E+00/
  ! DATA (GB(10, 1,IC),IC=1,3) /
  ! S 0.96825438E-02, 0.11766343E+00, 0.10000000E+01/
  ! DATA (GA(10, 2,IC),IC=1,3) /
  ! S 0.13108146E-01,-0.14708488E-02, 0.00000000E+00/
  ! DATA (GB(10, 2,IC),IC=1,3) /
  ! S 0.13108146E-01, 0.12019005E+00, 0.10000000E+01/

  ! ----- INTERVAL = 1 ----- T =  312.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA(11, 1,IC),IC=1,3) /
  ! S 0.10233955E-01,-0.12817135E-02, 0.00000000E+00/
  ! DATA (GB(11, 1,IC),IC=1,3) /
  ! S 0.10233955E-01, 0.11975320E+00, 0.10000000E+01/
  ! DATA (GA(11, 2,IC),IC=1,3) /
  ! S 0.13518390E-01,-0.15006791E-02, 0.00000000E+00/
  ! DATA (GB(11, 2,IC),IC=1,3) /
  ! S 0.13518390E-01, 0.11932684E+00, 0.10000000E+01/



  ! --- WATER VAPOR --- INTERVAL 2 -- 500-800 CM-1--- FROM ABS225 ---------




  ! --- R.D.  ---  G = 0.02 + 0.50 / ( 1 + 4.5 U )


  ! ----- INTERVAL = 2 ----- T =  187.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 1, 3,IC),IC=1,3) /
  ! S 0.11644593E+01, 0.41243390E+00, 0.00000000E+00/
  ! DATA (GB( 1, 3,IC),IC=1,3) /
  ! S 0.11644593E+01, 0.10346097E+01, 0.10000000E+01/
  ! DATA (GA( 1, 4,IC),IC=1,3) /
  ! S 0.12006968E+01, 0.48318936E+00, 0.00000000E+00/
  ! DATA (GB( 1, 4,IC),IC=1,3) /
  ! S 0.12006968E+01, 0.10626130E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  200.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 2, 3,IC),IC=1,3) /
  ! S 0.11747203E+01, 0.43407282E+00, 0.00000000E+00/
  ! DATA (GB( 2, 3,IC),IC=1,3) /
  ! S 0.11747203E+01, 0.10433655E+01, 0.10000000E+01/
  ! DATA (GA( 2, 4,IC),IC=1,3) /
  ! S 0.12108196E+01, 0.50501827E+00, 0.00000000E+00/
  ! DATA (GB( 2, 4,IC),IC=1,3) /
  ! S 0.12108196E+01, 0.10716026E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  212.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 3, 3,IC),IC=1,3) /
  ! S 0.11837872E+01, 0.45331413E+00, 0.00000000E+00/
  ! DATA (GB( 3, 3,IC),IC=1,3) /
  ! S 0.11837872E+01, 0.10511933E+01, 0.10000000E+01/
  ! DATA (GA( 3, 4,IC),IC=1,3) /
  ! S 0.12196717E+01, 0.52409502E+00, 0.00000000E+00/
  ! DATA (GB( 3, 4,IC),IC=1,3) /
  ! S 0.12196717E+01, 0.10795108E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  225.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 4, 3,IC),IC=1,3) /
  ! S 0.11918561E+01, 0.47048604E+00, 0.00000000E+00/
  ! DATA (GB( 4, 3,IC),IC=1,3) /
  ! S 0.11918561E+01, 0.10582150E+01, 0.10000000E+01/
  ! DATA (GA( 4, 4,IC),IC=1,3) /
  ! S 0.12274493E+01, 0.54085277E+00, 0.00000000E+00/
  ! DATA (GB( 4, 4,IC),IC=1,3) /
  ! S 0.12274493E+01, 0.10865006E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  237.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 5, 3,IC),IC=1,3) /
  ! S 0.11990757E+01, 0.48586286E+00, 0.00000000E+00/
  ! DATA (GB( 5, 3,IC),IC=1,3) /
  ! S 0.11990757E+01, 0.10645317E+01, 0.10000000E+01/
  ! DATA (GA( 5, 4,IC),IC=1,3) /
  ! S 0.12343189E+01, 0.55565422E+00, 0.00000000E+00/
  ! DATA (GB( 5, 4,IC),IC=1,3) /
  ! S 0.12343189E+01, 0.10927103E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  250.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 6, 3,IC),IC=1,3) /
  ! S 0.12055643E+01, 0.49968044E+00, 0.00000000E+00/
  ! DATA (GB( 6, 3,IC),IC=1,3) /
  ! S 0.12055643E+01, 0.10702313E+01, 0.10000000E+01/
  ! DATA (GA( 6, 4,IC),IC=1,3) /
  ! S 0.12404147E+01, 0.56878618E+00, 0.00000000E+00/
  ! DATA (GB( 6, 4,IC),IC=1,3) /
  ! S 0.12404147E+01, 0.10982489E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  262.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 7, 3,IC),IC=1,3) /
  ! S 0.12114186E+01, 0.51214132E+00, 0.00000000E+00/
  ! DATA (GB( 7, 3,IC),IC=1,3) /
  ! S 0.12114186E+01, 0.10753907E+01, 0.10000000E+01/
  ! DATA (GA( 7, 4,IC),IC=1,3) /
  ! S 0.12458431E+01, 0.58047395E+00, 0.00000000E+00/
  ! DATA (GB( 7, 4,IC),IC=1,3) /
  ! S 0.12458431E+01, 0.11032019E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  275.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 8, 3,IC),IC=1,3) /
  ! S 0.12167192E+01, 0.52341830E+00, 0.00000000E+00/
  ! DATA (GB( 8, 3,IC),IC=1,3) /
  ! S 0.12167192E+01, 0.10800762E+01, 0.10000000E+01/
  ! DATA (GA( 8, 4,IC),IC=1,3) /
  ! S 0.12506907E+01, 0.59089894E+00, 0.00000000E+00/
  ! DATA (GB( 8, 4,IC),IC=1,3) /
  ! S 0.12506907E+01, 0.11076379E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  287.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 9, 3,IC),IC=1,3) /
  ! S 0.12215344E+01, 0.53365803E+00, 0.00000000E+00/
  ! DATA (GB( 9, 3,IC),IC=1,3) /
  ! S 0.12215344E+01, 0.10843446E+01, 0.10000000E+01/
  ! DATA (GA( 9, 4,IC),IC=1,3) /
  ! S 0.12550299E+01, 0.60021475E+00, 0.00000000E+00/
  ! DATA (GB( 9, 4,IC),IC=1,3) /
  ! S 0.12550299E+01, 0.11116160E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  300.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA(10, 3,IC),IC=1,3) /
  ! S 0.12259226E+01, 0.54298448E+00, 0.00000000E+00/
  ! DATA (GB(10, 3,IC),IC=1,3) /
  ! S 0.12259226E+01, 0.10882439E+01, 0.10000000E+01/
  ! DATA (GA(10, 4,IC),IC=1,3) /
  ! S 0.12589256E+01, 0.60856112E+00, 0.00000000E+00/
  ! DATA (GB(10, 4,IC),IC=1,3) /
  ! S 0.12589256E+01, 0.11151910E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  312.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA(11, 3,IC),IC=1,3) /
  ! S 0.12299344E+01, 0.55150227E+00, 0.00000000E+00/
  ! DATA (GB(11, 3,IC),IC=1,3) /
  ! S 0.12299344E+01, 0.10918144E+01, 0.10000000E+01/
  ! DATA (GA(11, 4,IC),IC=1,3) /
  ! S 0.12624402E+01, 0.61607594E+00, 0.00000000E+00/
  ! DATA (GB(11, 4,IC),IC=1,3) /
  ! S 0.12624402E+01, 0.11184188E+01, 0.10000000E+01/






  ! - WATER VAPOR - INT. 3 -- 800-970 + 1110-1250 CM-1 -- FIT FROM 215 IS -


  ! -- WATER VAPOR LINES IN THE WINDOW REGION (800-1250 CM-1)



  ! --- G = 3.875E-03 ---------------

  ! ----- INTERVAL = 3 ----- T =  187.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 1, 7,IC),IC=1,3) /
  ! S 0.10192131E+02, 0.80737799E+01, 0.00000000E+00/
  ! DATA (GB( 1, 7,IC),IC=1,3) /
  ! S 0.10192131E+02, 0.82623280E+01, 0.10000000E+01/
  ! DATA (GA( 1, 8,IC),IC=1,3) /
  ! S 0.92439050E+01, 0.77425778E+01, 0.00000000E+00/
  ! DATA (GB( 1, 8,IC),IC=1,3) /
  ! S 0.92439050E+01, 0.79342219E+01, 0.10000000E+01/

  ! ----- INTERVAL = 3 ----- T =  200.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 2, 7,IC),IC=1,3) /
  ! S 0.97258602E+01, 0.79171158E+01, 0.00000000E+00/
  ! DATA (GB( 2, 7,IC),IC=1,3) /
  ! S 0.97258602E+01, 0.81072291E+01, 0.10000000E+01/
  ! DATA (GA( 2, 8,IC),IC=1,3) /
  ! S 0.87567422E+01, 0.75443460E+01, 0.00000000E+00/
  ! DATA (GB( 2, 8,IC),IC=1,3) /
  ! S 0.87567422E+01, 0.77373458E+01, 0.10000000E+01/

  ! ----- INTERVAL = 3 ----- T =  212.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 3, 7,IC),IC=1,3) /
  ! S 0.92992890E+01, 0.77609605E+01, 0.00000000E+00/
  ! DATA (GB( 3, 7,IC),IC=1,3) /
  ! S 0.92992890E+01, 0.79523834E+01, 0.10000000E+01/
  ! DATA (GA( 3, 8,IC),IC=1,3) /
  ! S 0.83270144E+01, 0.73526151E+01, 0.00000000E+00/
  ! DATA (GB( 3, 8,IC),IC=1,3) /
  ! S 0.83270144E+01, 0.75467334E+01, 0.10000000E+01/

  ! ----- INTERVAL = 3 ----- T =  225.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 4, 7,IC),IC=1,3) /
  ! S 0.89154021E+01, 0.76087371E+01, 0.00000000E+00/
  ! DATA (GB( 4, 7,IC),IC=1,3) /
  ! S 0.89154021E+01, 0.78012527E+01, 0.10000000E+01/
  ! DATA (GA( 4, 8,IC),IC=1,3) /
  ! S 0.79528337E+01, 0.71711188E+01, 0.00000000E+00/
  ! DATA (GB( 4, 8,IC),IC=1,3) /
  ! S 0.79528337E+01, 0.73661786E+01, 0.10000000E+01/

  ! ----- INTERVAL = 3 ----- T =  237.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 5, 7,IC),IC=1,3) /
  ! S 0.85730084E+01, 0.74627112E+01, 0.00000000E+00/
  ! DATA (GB( 5, 7,IC),IC=1,3) /
  ! S 0.85730084E+01, 0.76561458E+01, 0.10000000E+01/
  ! DATA (GA( 5, 8,IC),IC=1,3) /
  ! S 0.76286839E+01, 0.70015571E+01, 0.00000000E+00/
  ! DATA (GB( 5, 8,IC),IC=1,3) /
  ! S 0.76286839E+01, 0.71974319E+01, 0.10000000E+01/

  ! ----- INTERVAL = 3 ----- T =  250.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 6, 7,IC),IC=1,3) /
  ! S 0.82685838E+01, 0.73239981E+01, 0.00000000E+00/
  ! DATA (GB( 6, 7,IC),IC=1,3) /
  ! S 0.82685838E+01, 0.75182174E+01, 0.10000000E+01/
  ! DATA (GA( 6, 8,IC),IC=1,3) /
  ! S 0.73477879E+01, 0.68442532E+01, 0.00000000E+00/
  ! DATA (GB( 6, 8,IC),IC=1,3) /
  ! S 0.73477879E+01, 0.70408543E+01, 0.10000000E+01/

  ! ----- INTERVAL = 3 ----- T =  262.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 7, 7,IC),IC=1,3) /
  ! S 0.79978921E+01, 0.71929934E+01, 0.00000000E+00/
  ! DATA (GB( 7, 7,IC),IC=1,3) /
  ! S 0.79978921E+01, 0.73878952E+01, 0.10000000E+01/
  ! DATA (GA( 7, 8,IC),IC=1,3) /
  ! S 0.71035818E+01, 0.66987996E+01, 0.00000000E+00/
  ! DATA (GB( 7, 8,IC),IC=1,3) /
  ! S 0.71035818E+01, 0.68960649E+01, 0.10000000E+01/

  ! ----- INTERVAL = 3 ----- T =  275.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 8, 7,IC),IC=1,3) /
  ! S 0.77568055E+01, 0.70697065E+01, 0.00000000E+00/
  ! DATA (GB( 8, 7,IC),IC=1,3) /
  ! S 0.77568055E+01, 0.72652133E+01, 0.10000000E+01/
  ! DATA (GA( 8, 8,IC),IC=1,3) /
  ! S 0.68903312E+01, 0.65644820E+01, 0.00000000E+00/
  ! DATA (GB( 8, 8,IC),IC=1,3) /
  ! S 0.68903312E+01, 0.67623672E+01, 0.10000000E+01/

  ! ----- INTERVAL = 3 ----- T =  287.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 9, 7,IC),IC=1,3) /
  ! S 0.75416266E+01, 0.69539626E+01, 0.00000000E+00/
  ! DATA (GB( 9, 7,IC),IC=1,3) /
  ! S 0.75416266E+01, 0.71500151E+01, 0.10000000E+01/
  ! DATA (GA( 9, 8,IC),IC=1,3) /
  ! S 0.67032875E+01, 0.64405267E+01, 0.00000000E+00/
  ! DATA (GB( 9, 8,IC),IC=1,3) /
  ! S 0.67032875E+01, 0.66389989E+01, 0.10000000E+01/

  ! ----- INTERVAL = 3 ----- T =  300.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA(10, 7,IC),IC=1,3) /
  ! S 0.73491694E+01, 0.68455144E+01, 0.00000000E+00/
  ! DATA (GB(10, 7,IC),IC=1,3) /
  ! S 0.73491694E+01, 0.70420667E+01, 0.10000000E+01/
  ! DATA (GA(10, 8,IC),IC=1,3) /
  ! S 0.65386461E+01, 0.63262376E+01, 0.00000000E+00/
  ! DATA (GB(10, 8,IC),IC=1,3) /
  ! S 0.65386461E+01, 0.65252707E+01, 0.10000000E+01/

  ! ----- INTERVAL = 3 ----- T =  312.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA(11, 7,IC),IC=1,3) /
  ! S 0.71767400E+01, 0.67441020E+01, 0.00000000E+00/
  ! DATA (GB(11, 7,IC),IC=1,3) /
  ! S 0.71767400E+01, 0.69411177E+01, 0.10000000E+01/
  ! DATA (GA(11, 8,IC),IC=1,3) /
  ! S 0.63934377E+01, 0.62210701E+01, 0.00000000E+00/
  ! DATA (GB(11, 8,IC),IC=1,3) /
  ! S 0.63934377E+01, 0.64206412E+01, 0.10000000E+01/


  ! -- WATER VAPOR -- 970-1110 CM-1 ----------------------------------------

  ! -- G = 3.6E-03

  ! ----- INTERVAL = 4 ----- T =  187.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 1, 9,IC),IC=1,3) /
  ! S 0.24870635E+02, 0.10542131E+02, 0.00000000E+00/
  ! DATA (GB( 1, 9,IC),IC=1,3) /
  ! S 0.24870635E+02, 0.10656640E+02, 0.10000000E+01/
  ! DATA (GA( 1,10,IC),IC=1,3) /
  ! S 0.24586283E+02, 0.10490353E+02, 0.00000000E+00/
  ! DATA (GB( 1,10,IC),IC=1,3) /
  ! S 0.24586283E+02, 0.10605856E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  200.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 2, 9,IC),IC=1,3) /
  ! S 0.24725591E+02, 0.10515895E+02, 0.00000000E+00/
  ! DATA (GB( 2, 9,IC),IC=1,3) /
  ! S 0.24725591E+02, 0.10630910E+02, 0.10000000E+01/
  ! DATA (GA( 2,10,IC),IC=1,3) /
  ! S 0.24441465E+02, 0.10463512E+02, 0.00000000E+00/
  ! DATA (GB( 2,10,IC),IC=1,3) /
  ! S 0.24441465E+02, 0.10579514E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  212.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 3, 9,IC),IC=1,3) /
  ! S 0.24600320E+02, 0.10492949E+02, 0.00000000E+00/
  ! DATA (GB( 3, 9,IC),IC=1,3) /
  ! S 0.24600320E+02, 0.10608399E+02, 0.10000000E+01/
  ! DATA (GA( 3,10,IC),IC=1,3) /
  ! S 0.24311657E+02, 0.10439183E+02, 0.00000000E+00/
  ! DATA (GB( 3,10,IC),IC=1,3) /
  ! S 0.24311657E+02, 0.10555632E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  225.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 4, 9,IC),IC=1,3) /
  ! S 0.24487300E+02, 0.10472049E+02, 0.00000000E+00/
  ! DATA (GB( 4, 9,IC),IC=1,3) /
  ! S 0.24487300E+02, 0.10587891E+02, 0.10000000E+01/
  ! DATA (GA( 4,10,IC),IC=1,3) /
  ! S 0.24196167E+02, 0.10417324E+02, 0.00000000E+00/
  ! DATA (GB( 4,10,IC),IC=1,3) /
  ! S 0.24196167E+02, 0.10534169E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  237.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 5, 9,IC),IC=1,3) /
  ! S 0.24384935E+02, 0.10452961E+02, 0.00000000E+00/
  ! DATA (GB( 5, 9,IC),IC=1,3) /
  ! S 0.24384935E+02, 0.10569156E+02, 0.10000000E+01/
  ! DATA (GA( 5,10,IC),IC=1,3) /
  ! S 0.24093406E+02, 0.10397704E+02, 0.00000000E+00/
  ! DATA (GB( 5,10,IC),IC=1,3) /
  ! S 0.24093406E+02, 0.10514900E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  250.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 6, 9,IC),IC=1,3) /
  ! S 0.24292341E+02, 0.10435562E+02, 0.00000000E+00/
  ! DATA (GB( 6, 9,IC),IC=1,3) /
  ! S 0.24292341E+02, 0.10552075E+02, 0.10000000E+01/
  ! DATA (GA( 6,10,IC),IC=1,3) /
  ! S 0.24001597E+02, 0.10380038E+02, 0.00000000E+00/
  ! DATA (GB( 6,10,IC),IC=1,3) /
  ! S 0.24001597E+02, 0.10497547E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  262.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 7, 9,IC),IC=1,3) /
  ! S 0.24208572E+02, 0.10419710E+02, 0.00000000E+00/
  ! DATA (GB( 7, 9,IC),IC=1,3) /
  ! S 0.24208572E+02, 0.10536510E+02, 0.10000000E+01/
  ! DATA (GA( 7,10,IC),IC=1,3) /
  ! S 0.23919098E+02, 0.10364052E+02, 0.00000000E+00/
  ! DATA (GB( 7,10,IC),IC=1,3) /
  ! S 0.23919098E+02, 0.10481842E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  275.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 8, 9,IC),IC=1,3) /
  ! S 0.24132642E+02, 0.10405247E+02, 0.00000000E+00/
  ! DATA (GB( 8, 9,IC),IC=1,3) /
  ! S 0.24132642E+02, 0.10522307E+02, 0.10000000E+01/
  ! DATA (GA( 8,10,IC),IC=1,3) /
  ! S 0.23844511E+02, 0.10349509E+02, 0.00000000E+00/
  ! DATA (GB( 8,10,IC),IC=1,3) /
  ! S 0.23844511E+02, 0.10467553E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  287.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA( 9, 9,IC),IC=1,3) /
  ! S 0.24063614E+02, 0.10392022E+02, 0.00000000E+00/
  ! DATA (GB( 9, 9,IC),IC=1,3) /
  ! S 0.24063614E+02, 0.10509317E+02, 0.10000000E+01/
  ! DATA (GA( 9,10,IC),IC=1,3) /
  ! S 0.23776708E+02, 0.10336215E+02, 0.00000000E+00/
  ! DATA (GB( 9,10,IC),IC=1,3) /
  ! S 0.23776708E+02, 0.10454488E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  300.0

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA(10, 9,IC),IC=1,3) /
  ! S 0.24000649E+02, 0.10379892E+02, 0.00000000E+00/
  ! DATA (GB(10, 9,IC),IC=1,3) /
  ! S 0.24000649E+02, 0.10497402E+02, 0.10000000E+01/
  ! DATA (GA(10,10,IC),IC=1,3) /
  ! S 0.23714816E+02, 0.10324018E+02, 0.00000000E+00/
  ! DATA (GB(10,10,IC),IC=1,3) /
  ! S 0.23714816E+02, 0.10442501E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  312.5

  ! -- INDICES FOR PADE APPROXIMATION     1   28   37   45
  ! DATA (GA(11, 9,IC),IC=1,3) /
  ! S 0.23943021E+02, 0.10368736E+02, 0.00000000E+00/
  ! DATA (GB(11, 9,IC),IC=1,3) /
  ! S 0.23943021E+02, 0.10486443E+02, 0.10000000E+01/
  ! DATA (GA(11,10,IC),IC=1,3) /
  ! S 0.23658197E+02, 0.10312808E+02, 0.00000000E+00/
  ! DATA (GB(11,10,IC),IC=1,3) /
  ! S 0.23658197E+02, 0.10431483E+02, 0.10000000E+01/



  ! -- H2O -- WEAKER PARTS OF THE STRONG BANDS  -- FROM ABS225 ----

  ! -- WATER VAPOR --- 350 - 500 CM-1

  ! -- G = - 0.2*SLA, 0.0 +0.5/(1+0.5U)

  ! ----- INTERVAL = 5 ----- T =  187.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 1, 5,IC),IC=1,3) /
  ! S 0.15750172E+00,-0.22159303E-01, 0.00000000E+00/
  ! DATA (GB( 1, 5,IC),IC=1,3) /
  ! S 0.15750172E+00, 0.38103212E+00, 0.10000000E+01/
  ! DATA (GA( 1, 6,IC),IC=1,3) /
  ! S 0.17770551E+00,-0.24972399E-01, 0.00000000E+00/
  ! DATA (GB( 1, 6,IC),IC=1,3) /
  ! S 0.17770551E+00, 0.41646579E+00, 0.10000000E+01/

  ! ----- INTERVAL = 5 ----- T =  200.0

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 2, 5,IC),IC=1,3) /
  ! S 0.16174076E+00,-0.22748917E-01, 0.00000000E+00/
  ! DATA (GB( 2, 5,IC),IC=1,3) /
  ! S 0.16174076E+00, 0.38913800E+00, 0.10000000E+01/
  ! DATA (GA( 2, 6,IC),IC=1,3) /
  ! S 0.18176757E+00,-0.25537247E-01, 0.00000000E+00/
  ! DATA (GB( 2, 6,IC),IC=1,3) /
  ! S 0.18176757E+00, 0.42345095E+00, 0.10000000E+01/

  ! ----- INTERVAL = 5 ----- T =  212.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 3, 5,IC),IC=1,3) /
  ! S 0.16548628E+00,-0.23269898E-01, 0.00000000E+00/
  ! DATA (GB( 3, 5,IC),IC=1,3) /
  ! S 0.16548628E+00, 0.39613651E+00, 0.10000000E+01/
  ! DATA (GA( 3, 6,IC),IC=1,3) /
  ! S 0.18527967E+00,-0.26025624E-01, 0.00000000E+00/
  ! DATA (GB( 3, 6,IC),IC=1,3) /
  ! S 0.18527967E+00, 0.42937476E+00, 0.10000000E+01/

  ! ----- INTERVAL = 5 ----- T =  225.0

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 4, 5,IC),IC=1,3) /
  ! S 0.16881124E+00,-0.23732392E-01, 0.00000000E+00/
  ! DATA (GB( 4, 5,IC),IC=1,3) /
  ! S 0.16881124E+00, 0.40222421E+00, 0.10000000E+01/
  ! DATA (GA( 4, 6,IC),IC=1,3) /
  ! S 0.18833348E+00,-0.26450280E-01, 0.00000000E+00/
  ! DATA (GB( 4, 6,IC),IC=1,3) /
  ! S 0.18833348E+00, 0.43444062E+00, 0.10000000E+01/

  ! ----- INTERVAL = 5 ----- T =  237.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 5, 5,IC),IC=1,3) /
  ! S 0.17177839E+00,-0.24145123E-01, 0.00000000E+00/
  ! DATA (GB( 5, 5,IC),IC=1,3) /
  ! S 0.17177839E+00, 0.40756010E+00, 0.10000000E+01/
  ! DATA (GA( 5, 6,IC),IC=1,3) /
  ! S 0.19100108E+00,-0.26821236E-01, 0.00000000E+00/
  ! DATA (GB( 5, 6,IC),IC=1,3) /
  ! S 0.19100108E+00, 0.43880316E+00, 0.10000000E+01/

  ! ----- INTERVAL = 5 ----- T =  250.0

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 6, 5,IC),IC=1,3) /
  ! S 0.17443933E+00,-0.24515269E-01, 0.00000000E+00/
  ! DATA (GB( 6, 5,IC),IC=1,3) /
  ! S 0.17443933E+00, 0.41226954E+00, 0.10000000E+01/
  ! DATA (GA( 6, 6,IC),IC=1,3) /
  ! S 0.19334122E+00,-0.27146657E-01, 0.00000000E+00/
  ! DATA (GB( 6, 6,IC),IC=1,3) /
  ! S 0.19334122E+00, 0.44258354E+00, 0.10000000E+01/

  ! ----- INTERVAL = 5 ----- T =  262.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 7, 5,IC),IC=1,3) /
  ! S 0.17683622E+00,-0.24848690E-01, 0.00000000E+00/
  ! DATA (GB( 7, 5,IC),IC=1,3) /
  ! S 0.17683622E+00, 0.41645142E+00, 0.10000000E+01/
  ! DATA (GA( 7, 6,IC),IC=1,3) /
  ! S 0.19540288E+00,-0.27433354E-01, 0.00000000E+00/
  ! DATA (GB( 7, 6,IC),IC=1,3) /
  ! S 0.19540288E+00, 0.44587882E+00, 0.10000000E+01/

  ! ----- INTERVAL = 5 ----- T =  275.0

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 8, 5,IC),IC=1,3) /
  ! S 0.17900375E+00,-0.25150210E-01, 0.00000000E+00/
  ! DATA (GB( 8, 5,IC),IC=1,3) /
  ! S 0.17900375E+00, 0.42018474E+00, 0.10000000E+01/
  ! DATA (GA( 8, 6,IC),IC=1,3) /
  ! S 0.19722732E+00,-0.27687065E-01, 0.00000000E+00/
  ! DATA (GB( 8, 6,IC),IC=1,3) /
  ! S 0.19722732E+00, 0.44876776E+00, 0.10000000E+01/

  ! ----- INTERVAL = 5 ----- T =  287.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 9, 5,IC),IC=1,3) /
  ! S 0.18097099E+00,-0.25423873E-01, 0.00000000E+00/
  ! DATA (GB( 9, 5,IC),IC=1,3) /
  ! S 0.18097099E+00, 0.42353379E+00, 0.10000000E+01/
  ! DATA (GA( 9, 6,IC),IC=1,3) /
  ! S 0.19884918E+00,-0.27912608E-01, 0.00000000E+00/
  ! DATA (GB( 9, 6,IC),IC=1,3) /
  ! S 0.19884918E+00, 0.45131451E+00, 0.10000000E+01/

  ! ----- INTERVAL = 5 ----- T =  300.0

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA(10, 5,IC),IC=1,3) /
  ! S 0.18276283E+00,-0.25673139E-01, 0.00000000E+00/
  ! DATA (GB(10, 5,IC),IC=1,3) /
  ! S 0.18276283E+00, 0.42655211E+00, 0.10000000E+01/
  ! DATA (GA(10, 6,IC),IC=1,3) /
  ! S 0.20029696E+00,-0.28113944E-01, 0.00000000E+00/
  ! DATA (GB(10, 6,IC),IC=1,3) /
  ! S 0.20029696E+00, 0.45357095E+00, 0.10000000E+01/

  ! ----- INTERVAL = 5 ----- T =  312.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA(11, 5,IC),IC=1,3) /
  ! S 0.18440117E+00,-0.25901055E-01, 0.00000000E+00/
  ! DATA (GB(11, 5,IC),IC=1,3) /
  ! S 0.18440117E+00, 0.42928533E+00, 0.10000000E+01/
  ! DATA (GA(11, 6,IC),IC=1,3) /
  ! S 0.20159300E+00,-0.28294180E-01, 0.00000000E+00/
  ! DATA (GB(11, 6,IC),IC=1,3) /
  ! S 0.20159300E+00, 0.45557797E+00, 0.10000000E+01/




  ! - WATER VAPOR - WINGS OF VIBRATION-ROTATION BAND - 1250-1450+1880-2820 -
  ! --- G = 0.0


  ! ----- INTERVAL = 6 ----- T =  187.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 1,11,IC),IC=1,3) /
  ! S 0.11990218E+02,-0.12823142E+01, 0.00000000E+00/
  ! DATA (GB( 1,11,IC),IC=1,3) /
  ! S 0.11990218E+02, 0.26681588E+02, 0.10000000E+01/
  ! DATA (GA( 1,12,IC),IC=1,3) /
  ! S 0.79709806E+01,-0.74805226E+00, 0.00000000E+00/
  ! DATA (GB( 1,12,IC),IC=1,3) /
  ! S 0.79709806E+01, 0.18377807E+02, 0.10000000E+01/

  ! ----- INTERVAL = 6 ----- T =  200.0

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 2,11,IC),IC=1,3) /
  ! S 0.10904073E+02,-0.10571588E+01, 0.00000000E+00/
  ! DATA (GB( 2,11,IC),IC=1,3) /
  ! S 0.10904073E+02, 0.24728346E+02, 0.10000000E+01/
  ! DATA (GA( 2,12,IC),IC=1,3) /
  ! S 0.75400737E+01,-0.56252739E+00, 0.00000000E+00/
  ! DATA (GB( 2,12,IC),IC=1,3) /
  ! S 0.75400737E+01, 0.17643148E+02, 0.10000000E+01/

  ! ----- INTERVAL = 6 ----- T =  212.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 3,11,IC),IC=1,3) /
  ! S 0.89126838E+01,-0.74864953E+00, 0.00000000E+00/
  ! DATA (GB( 3,11,IC),IC=1,3) /
  ! S 0.89126838E+01, 0.20551342E+02, 0.10000000E+01/
  ! DATA (GA( 3,12,IC),IC=1,3) /
  ! S 0.81804377E+01,-0.46188072E+00, 0.00000000E+00/
  ! DATA (GB( 3,12,IC),IC=1,3) /
  ! S 0.81804377E+01, 0.19296161E+02, 0.10000000E+01/

  ! ----- INTERVAL = 6 ----- T =  225.0

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 4,11,IC),IC=1,3) /
  ! S 0.85622405E+01,-0.58705980E+00, 0.00000000E+00/
  ! DATA (GB( 4,11,IC),IC=1,3) /
  ! S 0.85622405E+01, 0.19955244E+02, 0.10000000E+01/
  ! DATA (GA( 4,12,IC),IC=1,3) /
  ! S 0.10564339E+02,-0.40712065E+00, 0.00000000E+00/
  ! DATA (GB( 4,12,IC),IC=1,3) /
  ! S 0.10564339E+02, 0.24951120E+02, 0.10000000E+01/

  ! ----- INTERVAL = 6 ----- T =  237.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 5,11,IC),IC=1,3) /
  ! S 0.94892164E+01,-0.49305772E+00, 0.00000000E+00/
  ! DATA (GB( 5,11,IC),IC=1,3) /
  ! S 0.94892164E+01, 0.22227100E+02, 0.10000000E+01/
  ! DATA (GA( 5,12,IC),IC=1,3) /
  ! S 0.46896789E+02,-0.15295996E+01, 0.00000000E+00/
  ! DATA (GB( 5,12,IC),IC=1,3) /
  ! S 0.46896789E+02, 0.10957372E+03, 0.10000000E+01/

  ! ----- INTERVAL = 6 ----- T =  250.0

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 6,11,IC),IC=1,3) /
  ! S 0.13580937E+02,-0.51461431E+00, 0.00000000E+00/
  ! DATA (GB( 6,11,IC),IC=1,3) /
  ! S 0.13580937E+02, 0.31770288E+02, 0.10000000E+01/
  ! DATA (GA( 6,12,IC),IC=1,3) /
  ! S-0.30926524E+01, 0.43555255E+00, 0.00000000E+00/
  ! DATA (GB( 6,12,IC),IC=1,3) /
  ! S-0.30926524E+01,-0.67432659E+01, 0.10000000E+01/

  ! ----- INTERVAL = 6 ----- T =  262.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 7,11,IC),IC=1,3) /
  ! S-0.32050918E+03, 0.12373350E+02, 0.00000000E+00/
  ! DATA (GB( 7,11,IC),IC=1,3) /
  ! S-0.32050918E+03,-0.74061287E+03, 0.10000000E+01/
  ! DATA (GA( 7,12,IC),IC=1,3) /
  ! S 0.85742941E+00, 0.50380874E+00, 0.00000000E+00/
  ! DATA (GB( 7,12,IC),IC=1,3) /
  ! S 0.85742941E+00, 0.24550746E+01, 0.10000000E+01/

  ! ----- INTERVAL = 6 ----- T =  275.0

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 8,11,IC),IC=1,3) /
  ! S-0.37133165E+01, 0.44809588E+00, 0.00000000E+00/
  ! DATA (GB( 8,11,IC),IC=1,3) /
  ! S-0.37133165E+01,-0.81329826E+01, 0.10000000E+01/
  ! DATA (GA( 8,12,IC),IC=1,3) /
  ! S 0.19164038E+01, 0.68537352E+00, 0.00000000E+00/
  ! DATA (GB( 8,12,IC),IC=1,3) /
  ! S 0.19164038E+01, 0.49089917E+01, 0.10000000E+01/

  ! ----- INTERVAL = 6 ----- T =  287.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA( 9,11,IC),IC=1,3) /
  ! S 0.18890836E+00, 0.46548918E+00, 0.00000000E+00/
  ! DATA (GB( 9,11,IC),IC=1,3) /
  ! S 0.18890836E+00, 0.90279822E+00, 0.10000000E+01/
  ! DATA (GA( 9,12,IC),IC=1,3) /
  ! S 0.23513199E+01, 0.89437630E+00, 0.00000000E+00/
  ! DATA (GB( 9,12,IC),IC=1,3) /
  ! S 0.23513199E+01, 0.59008712E+01, 0.10000000E+01/

  ! ----- INTERVAL = 6 ----- T =  300.0

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA(10,11,IC),IC=1,3) /
  ! S 0.14209226E+01, 0.59121475E+00, 0.00000000E+00/
  ! DATA (GB(10,11,IC),IC=1,3) /
  ! S 0.14209226E+01, 0.37532746E+01, 0.10000000E+01/
  ! DATA (GA(10,12,IC),IC=1,3) /
  ! S 0.25566644E+01, 0.11127003E+01, 0.00000000E+00/
  ! DATA (GB(10,12,IC),IC=1,3) /
  ! S 0.25566644E+01, 0.63532616E+01, 0.10000000E+01/

  ! ----- INTERVAL = 6 ----- T =  312.5

  ! -- INDICES FOR PADE APPROXIMATION   1 35 40 45
  ! DATA (GA(11,11,IC),IC=1,3) /
  ! S 0.19817679E+01, 0.74676119E+00, 0.00000000E+00/
  ! DATA (GB(11,11,IC),IC=1,3) /
  ! S 0.19817679E+01, 0.50437916E+01, 0.10000000E+01/
  ! DATA (GA(11,12,IC),IC=1,3) /
  ! S 0.26555181E+01, 0.13329782E+01, 0.00000000E+00/
  ! DATA (GB(11,12,IC),IC=1,3) /
  ! S 0.26555181E+01, 0.65558627E+01, 0.10000000E+01/





  ! -- END WATER VAPOR


  ! -- CO2 -- INT.2 -- 500-800 CM-1 --- FROM ABS225 ----------------------



  ! -- FIU = 0.8 + MAX(0.35,(7-IU)*0.9)  , X/T,  9

  ! ----- INTERVAL = 2 ----- T =  187.5

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA( 1,13,IC),IC=1,3) /
  ! S 0.87668459E-01, 0.13845511E+01, 0.00000000E+00/
  ! DATA (GB( 1,13,IC),IC=1,3) /
  ! S 0.87668459E-01, 0.23203798E+01, 0.10000000E+01/
  ! DATA (GA( 1,14,IC),IC=1,3) /
  ! S 0.74878820E-01, 0.11718758E+01, 0.00000000E+00/
  ! DATA (GB( 1,14,IC),IC=1,3) /
  ! S 0.74878820E-01, 0.20206726E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  200.0

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA( 2,13,IC),IC=1,3) /
  ! S 0.83754276E-01, 0.13187042E+01, 0.00000000E+00/
  ! DATA (GB( 2,13,IC),IC=1,3) /
  ! S 0.83754276E-01, 0.22288925E+01, 0.10000000E+01/
  ! DATA (GA( 2,14,IC),IC=1,3) /
  ! S 0.71650966E-01, 0.11216131E+01, 0.00000000E+00/
  ! DATA (GB( 2,14,IC),IC=1,3) /
  ! S 0.71650966E-01, 0.19441824E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  212.5

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA( 3,13,IC),IC=1,3) /
  ! S 0.80460283E-01, 0.12644396E+01, 0.00000000E+00/
  ! DATA (GB( 3,13,IC),IC=1,3) /
  ! S 0.80460283E-01, 0.21515593E+01, 0.10000000E+01/
  ! DATA (GA( 3,14,IC),IC=1,3) /
  ! S 0.68979615E-01, 0.10809473E+01, 0.00000000E+00/
  ! DATA (GB( 3,14,IC),IC=1,3) /
  ! S 0.68979615E-01, 0.18807257E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  225.0

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA( 4,13,IC),IC=1,3) /
  ! S 0.77659686E-01, 0.12191543E+01, 0.00000000E+00/
  ! DATA (GB( 4,13,IC),IC=1,3) /
  ! S 0.77659686E-01, 0.20855896E+01, 0.10000000E+01/
  ! DATA (GA( 4,14,IC),IC=1,3) /
  ! S 0.66745345E-01, 0.10476396E+01, 0.00000000E+00/
  ! DATA (GB( 4,14,IC),IC=1,3) /
  ! S 0.66745345E-01, 0.18275618E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  237.5

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA( 5,13,IC),IC=1,3) /
  ! S 0.75257056E-01, 0.11809511E+01, 0.00000000E+00/
  ! DATA (GB( 5,13,IC),IC=1,3) /
  ! S 0.75257056E-01, 0.20288489E+01, 0.10000000E+01/
  ! DATA (GA( 5,14,IC),IC=1,3) /
  ! S 0.64857571E-01, 0.10200373E+01, 0.00000000E+00/
  ! DATA (GB( 5,14,IC),IC=1,3) /
  ! S 0.64857571E-01, 0.17825910E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  250.0

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA( 6,13,IC),IC=1,3) /
  ! S 0.73179175E-01, 0.11484154E+01, 0.00000000E+00/
  ! DATA (GB( 6,13,IC),IC=1,3) /
  ! S 0.73179175E-01, 0.19796791E+01, 0.10000000E+01/
  ! DATA (GA( 6,14,IC),IC=1,3) /
  ! S 0.63248495E-01, 0.99692726E+00, 0.00000000E+00/
  ! DATA (GB( 6,14,IC),IC=1,3) /
  ! S 0.63248495E-01, 0.17442308E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  262.5

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA( 7,13,IC),IC=1,3) /
  ! S 0.71369063E-01, 0.11204723E+01, 0.00000000E+00/
  ! DATA (GB( 7,13,IC),IC=1,3) /
  ! S 0.71369063E-01, 0.19367778E+01, 0.10000000E+01/
  ! DATA (GA( 7,14,IC),IC=1,3) /
  ! S 0.61866970E-01, 0.97740923E+00, 0.00000000E+00/
  ! DATA (GB( 7,14,IC),IC=1,3) /
  ! S 0.61866970E-01, 0.17112809E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  275.0

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA( 8,13,IC),IC=1,3) /
  ! S 0.69781812E-01, 0.10962918E+01, 0.00000000E+00/
  ! DATA (GB( 8,13,IC),IC=1,3) /
  ! S 0.69781812E-01, 0.18991112E+01, 0.10000000E+01/
  ! DATA (GA( 8,14,IC),IC=1,3) /
  ! S 0.60673632E-01, 0.96080188E+00, 0.00000000E+00/
  ! DATA (GB( 8,14,IC),IC=1,3) /
  ! S 0.60673632E-01, 0.16828137E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  287.5

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA( 9,13,IC),IC=1,3) /
  ! S 0.68381606E-01, 0.10752229E+01, 0.00000000E+00/
  ! DATA (GB( 9,13,IC),IC=1,3) /
  ! S 0.68381606E-01, 0.18658501E+01, 0.10000000E+01/
  ! DATA (GA( 9,14,IC),IC=1,3) /
  ! S 0.59637277E-01, 0.94657562E+00, 0.00000000E+00/
  ! DATA (GB( 9,14,IC),IC=1,3) /
  ! S 0.59637277E-01, 0.16580908E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  300.0

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA(10,13,IC),IC=1,3) /
  ! S 0.67139539E-01, 0.10567474E+01, 0.00000000E+00/
  ! DATA (GB(10,13,IC),IC=1,3) /
  ! S 0.67139539E-01, 0.18363226E+01, 0.10000000E+01/
  ! DATA (GA(10,14,IC),IC=1,3) /
  ! S 0.58732178E-01, 0.93430511E+00, 0.00000000E+00/
  ! DATA (GB(10,14,IC),IC=1,3) /
  ! S 0.58732178E-01, 0.16365014E+01, 0.10000000E+01/

  ! ----- INTERVAL = 2 ----- T =  312.5

  ! -- INDICES FOR PADE APPROXIMATION   1 30 38 45
  ! DATA (GA(11,13,IC),IC=1,3) /
  ! S 0.66032012E-01, 0.10404465E+01, 0.00000000E+00/
  ! DATA (GB(11,13,IC),IC=1,3) /
  ! S 0.66032012E-01, 0.18099779E+01, 0.10000000E+01/
  ! DATA (GA(11,14,IC),IC=1,3) /
  ! S 0.57936092E-01, 0.92363528E+00, 0.00000000E+00/
  ! DATA (GB(11,14,IC),IC=1,3) /
  ! S 0.57936092E-01, 0.16175164E+01, 0.10000000E+01/










  ! -- CARBON DIOXIDE LINES IN THE WINDOW REGION (800-1250 CM-1)


  ! -- G = 0.0


  ! ----- INTERVAL = 4 ----- T =  187.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 1,15,IC),IC=1,3) /
  ! S 0.13230067E+02, 0.22042132E+02, 0.00000000E+00/
  ! DATA (GB( 1,15,IC),IC=1,3) /
  ! S 0.13230067E+02, 0.22051750E+02, 0.10000000E+01/
  ! DATA (GA( 1,16,IC),IC=1,3) /
  ! S 0.13183816E+02, 0.22169501E+02, 0.00000000E+00/
  ! DATA (GB( 1,16,IC),IC=1,3) /
  ! S 0.13183816E+02, 0.22178972E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  200.0

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 2,15,IC),IC=1,3) /
  ! S 0.13213564E+02, 0.22107298E+02, 0.00000000E+00/
  ! DATA (GB( 2,15,IC),IC=1,3) /
  ! S 0.13213564E+02, 0.22116850E+02, 0.10000000E+01/
  ! DATA (GA( 2,16,IC),IC=1,3) /
  ! S 0.13189991E+02, 0.22270075E+02, 0.00000000E+00/
  ! DATA (GB( 2,16,IC),IC=1,3) /
  ! S 0.13189991E+02, 0.22279484E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  212.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 3,15,IC),IC=1,3) /
  ! S 0.13209140E+02, 0.22180915E+02, 0.00000000E+00/
  ! DATA (GB( 3,15,IC),IC=1,3) /
  ! S 0.13209140E+02, 0.22190410E+02, 0.10000000E+01/
  ! DATA (GA( 3,16,IC),IC=1,3) /
  ! S 0.13209485E+02, 0.22379193E+02, 0.00000000E+00/
  ! DATA (GB( 3,16,IC),IC=1,3) /
  ! S 0.13209485E+02, 0.22388551E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  225.0

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 4,15,IC),IC=1,3) /
  ! S 0.13213894E+02, 0.22259478E+02, 0.00000000E+00/
  ! DATA (GB( 4,15,IC),IC=1,3) /
  ! S 0.13213894E+02, 0.22268925E+02, 0.10000000E+01/
  ! DATA (GA( 4,16,IC),IC=1,3) /
  ! S 0.13238789E+02, 0.22492992E+02, 0.00000000E+00/
  ! DATA (GB( 4,16,IC),IC=1,3) /
  ! S 0.13238789E+02, 0.22502309E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  237.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 5,15,IC),IC=1,3) /
  ! S 0.13225963E+02, 0.22341039E+02, 0.00000000E+00/
  ! DATA (GB( 5,15,IC),IC=1,3) /
  ! S 0.13225963E+02, 0.22350445E+02, 0.10000000E+01/
  ! DATA (GA( 5,16,IC),IC=1,3) /
  ! S 0.13275017E+02, 0.22608508E+02, 0.00000000E+00/
  ! DATA (GB( 5,16,IC),IC=1,3) /
  ! S 0.13275017E+02, 0.22617792E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  250.0

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 6,15,IC),IC=1,3) /
  ! S 0.13243806E+02, 0.22424247E+02, 0.00000000E+00/
  ! DATA (GB( 6,15,IC),IC=1,3) /
  ! S 0.13243806E+02, 0.22433617E+02, 0.10000000E+01/
  ! DATA (GA( 6,16,IC),IC=1,3) /
  ! S 0.13316096E+02, 0.22723843E+02, 0.00000000E+00/
  ! DATA (GB( 6,16,IC),IC=1,3) /
  ! S 0.13316096E+02, 0.22733099E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  262.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 7,15,IC),IC=1,3) /
  ! S 0.13266104E+02, 0.22508089E+02, 0.00000000E+00/
  ! DATA (GB( 7,15,IC),IC=1,3) /
  ! S 0.13266104E+02, 0.22517429E+02, 0.10000000E+01/
  ! DATA (GA( 7,16,IC),IC=1,3) /
  ! S 0.13360555E+02, 0.22837837E+02, 0.00000000E+00/
  ! DATA (GB( 7,16,IC),IC=1,3) /
  ! S 0.13360555E+02, 0.22847071E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  275.0

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 8,15,IC),IC=1,3) /
  ! S 0.13291782E+02, 0.22591771E+02, 0.00000000E+00/
  ! DATA (GB( 8,15,IC),IC=1,3) /
  ! S 0.13291782E+02, 0.22601086E+02, 0.10000000E+01/
  ! DATA (GA( 8,16,IC),IC=1,3) /
  ! S 0.13407324E+02, 0.22949751E+02, 0.00000000E+00/
  ! DATA (GB( 8,16,IC),IC=1,3) /
  ! S 0.13407324E+02, 0.22958967E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  287.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA( 9,15,IC),IC=1,3) /
  ! S 0.13319961E+02, 0.22674661E+02, 0.00000000E+00/
  ! DATA (GB( 9,15,IC),IC=1,3) /
  ! S 0.13319961E+02, 0.22683956E+02, 0.10000000E+01/
  ! DATA (GA( 9,16,IC),IC=1,3) /
  ! S 0.13455544E+02, 0.23059032E+02, 0.00000000E+00/
  ! DATA (GB( 9,16,IC),IC=1,3) /
  ! S 0.13455544E+02, 0.23068234E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  300.0

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA(10,15,IC),IC=1,3) /
  ! S 0.13349927E+02, 0.22756246E+02, 0.00000000E+00/
  ! DATA (GB(10,15,IC),IC=1,3) /
  ! S 0.13349927E+02, 0.22765522E+02, 0.10000000E+01/
  ! DATA (GA(10,16,IC),IC=1,3) /
  ! S 0.13504450E+02, 0.23165146E+02, 0.00000000E+00/
  ! DATA (GB(10,16,IC),IC=1,3) /
  ! S 0.13504450E+02, 0.23174336E+02, 0.10000000E+01/

  ! ----- INTERVAL = 4 ----- T =  312.5

  ! -- INDICES FOR PADE APPROXIMATION     1   15   29   45
  ! DATA (GA(11,15,IC),IC=1,3) /
  ! S 0.13381108E+02, 0.22836093E+02, 0.00000000E+00/
  ! DATA (GB(11,15,IC),IC=1,3) /
  ! S 0.13381108E+02, 0.22845354E+02, 0.10000000E+01/
  ! DATA (GA(11,16,IC),IC=1,3) /
  ! S 0.13553282E+02, 0.23267456E+02, 0.00000000E+00/
  ! DATA (GB(11,16,IC),IC=1,3) /
  ! S 0.13553282E+02, 0.23276638E+02, 0.10000000E+01/

  ! ------------------------------------------------------------------
  ! DATA (( XP(  J,K),J=1,6),       K=1,6) /
  ! S 0.46430621E+02, 0.12928299E+03, 0.20732648E+03,
  ! S 0.31398411E+03, 0.18373177E+03,-0.11412303E+03,
  ! S 0.73604774E+02, 0.27887914E+03, 0.27076947E+03,
  ! S-0.57322111E+02,-0.64742459E+02, 0.87238280E+02,
  ! S 0.37050866E+02, 0.20498759E+03, 0.37558029E+03,
  ! S 0.17401171E+03,-0.13350302E+03,-0.37651795E+02,
  ! S 0.14930141E+02, 0.89161160E+02, 0.17793062E+03,
  ! S 0.93433860E+02,-0.70646020E+02,-0.26373150E+02,
  ! S 0.40386780E+02, 0.10855270E+03, 0.50755010E+02,
  ! S-0.31496190E+02, 0.12791300E+00, 0.18017770E+01,
  ! S 0.90811926E+01, 0.75073923E+02, 0.24654438E+03,
  ! S 0.39332612E+03, 0.29385281E+03, 0.89107921E+02 /



  ! *         1.0     PLANCK FUNCTIONS AND GRADIENTS
  ! ------------------------------


  ! cdir collapse
  DO jk = 1, kflev + 1
    DO jl = 1, kdlon
      pbint(jl, jk) = 0.
    END DO
  END DO
  DO jl = 1, kdlon
    pbsuin(jl) = 0.
  END DO

  DO jnu = 1, ninter

    ! *         1.1   LEVELS FROM SURFACE TO KFLEV
    ! ----------------------------


    DO jk = 1, kflev
      DO jl = 1, kdlon
        zti(jl) = (ptl(jl,jk)-tstand)/tstand
        zres(jl) = xp(1, jnu) + zti(jl)*(xp(2,jnu)+zti(jl)*(xp(3, &
          jnu)+zti(jl)*(xp(4,jnu)+zti(jl)*(xp(5,jnu)+zti(jl)*(xp(6,jnu))))))
        pbint(jl, jk) = pbint(jl, jk) + zres(jl)
        pb(jl, jnu, jk) = zres(jl)
        zblev(jl, jk) = zres(jl)
        zti2(jl) = (ptave(jl,jk)-tstand)/tstand
        zres2(jl) = xp(1, jnu) + zti2(jl)*(xp(2,jnu)+zti2(jl)*(xp(3, &
          jnu)+zti2(jl)*(xp(4,jnu)+zti2(jl)*(xp(5,jnu)+zti2(jl)*(xp(6,jnu)))) &
          ))
        zblay(jl, jk) = zres2(jl)
      END DO
    END DO

    ! *         1.2   TOP OF THE ATMOSPHERE AND SURFACE
    ! ---------------------------------


    DO jl = 1, kdlon
      zti(jl) = (ptl(jl,kflev+1)-tstand)/tstand
      zti2(jl) = (ptl(jl,1)+pdt0(jl)-tstand)/tstand
      zres(jl) = xp(1, jnu) + zti(jl)*(xp(2,jnu)+zti(jl)*(xp(3, &
        jnu)+zti(jl)*(xp(4,jnu)+zti(jl)*(xp(5,jnu)+zti(jl)*(xp(6,jnu))))))
      zres2(jl) = xp(1, jnu) + zti2(jl)*(xp(2,jnu)+zti2(jl)*(xp(3, &
        jnu)+zti2(jl)*(xp(4,jnu)+zti2(jl)*(xp(5,jnu)+zti2(jl)*(xp(6,jnu))))))
      pbint(jl, kflev+1) = pbint(jl, kflev+1) + zres(jl)
      pb(jl, jnu, kflev+1) = zres(jl)
      zblev(jl, kflev+1) = zres(jl)
      pbtop(jl, jnu) = zres(jl)
      pbsur(jl, jnu) = zres2(jl)
      pbsuin(jl) = pbsuin(jl) + zres2(jl)
    END DO

    ! *         1.3   GRADIENTS IN SUB-LAYERS
    ! -----------------------


    DO jk = 1, kflev
      jk2 = 2*jk
      jk1 = jk2 - 1
      DO jl = 1, kdlon
        pdbsl(jl, jnu, jk1) = zblay(jl, jk) - zblev(jl, jk)
        pdbsl(jl, jnu, jk2) = zblev(jl, jk+1) - zblay(jl, jk)
      END DO
    END DO

  END DO

  ! *         2.0   CHOOSE THE RELEVANT SETS OF PADE APPROXIMANTS
  ! ---------------------------------------------




  DO jl = 1, kdlon
    zdsto1 = (ptl(jl,kflev+1)-tintp(1))/tstp
    ixtox = max(1, min(mxixt,int(zdsto1+1.)))
    zdstox = (ptl(jl,kflev+1)-tintp(ixtox))/tstp
    IF (zdstox<0.5) THEN
      indto = ixtox
    ELSE
      indto = ixtox + 1
    END IF
    indb(jl) = indto
    zdst1 = (ptl(jl,1)-tintp(1))/tstp
    ixtx = max(1, min(mxixt,int(zdst1+1.)))
    zdstx = (ptl(jl,1)-tintp(ixtx))/tstp
    IF (zdstx<0.5) THEN
      indt = ixtx
    ELSE
      indt = ixtx + 1
    END IF
    inds(jl) = indt
  END DO

  DO jf = 1, 2
    DO jg = 1, 8
      DO jl = 1, kdlon
        indsu = inds(jl)
        pgasur(jl, jg, jf) = ga(indsu, 2*jg-1, jf)
        pgbsur(jl, jg, jf) = gb(indsu, 2*jg-1, jf)
        indtp = indb(jl)
        pgatop(jl, jg, jf) = ga(indtp, 2*jg-1, jf)
        pgbtop(jl, jg, jf) = gb(indtp, 2*jg-1, jf)
      END DO
    END DO
  END DO

  DO jk = 1, kflev
    DO jl = 1, kdlon
      zdst1 = (ptave(jl,jk)-tintp(1))/tstp
      ixtx = max(1, min(mxixt,int(zdst1+1.)))
      zdstx = (ptave(jl,jk)-tintp(ixtx))/tstp
      IF (zdstx<0.5) THEN
        indt = ixtx
      ELSE
        indt = ixtx + 1
      END IF
      indb(jl) = indt
    END DO

    DO jf = 1, 2
      DO jg = 1, 8
        DO jl = 1, kdlon
          indt = indb(jl)
          pga(jl, jg, jf, jk) = ga(indt, 2*jg, jf)
          pgb(jl, jg, jf, jk) = gb(indt, 2*jg, jf)
        END DO
      END DO
    END DO
  END DO

  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE lwb_lmdar4
SUBROUTINE lwv_lmdar4(kuaer, ktraer, klim, pabcu, pb, pbint, pbsuin, pbsur, &
    pbtop, pdbsl, pemis, ppmb, ptave, pga, pgb, pgasur, pgbsur, pgatop, &
    pgbtop, pcntrb, pcts, pfluc)
  USE dimphy
  IMPLICIT NONE
  include "raddimlw.h"
  include "YOMCST.h"

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! CARRIES OUT THE VERTICAL INTEGRATION TO GIVE LONGWAVE
  ! FLUXES OR RADIANCES

  ! METHOD.
  ! -------

  ! 1. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING BETWEEN
  ! CONTRIBUTIONS BY -  THE NEARBY LAYERS
  ! -  THE DISTANT LAYERS
  ! -  THE BOUNDARY TERMS
  ! 2. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! -----------------------------------------------------------------------

  ! * ARGUMENTS:
  INTEGER kuaer, ktraer, klim

  REAL (KIND=8) pabcu(kdlon, nua, 3*kflev+1) ! EFFECTIVE ABSORBER AMOUNTS
  REAL (KIND=8) pb(kdlon, ninter, kflev+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
  REAL (KIND=8) pbint(kdlon, kflev+1) ! HALF-LEVEL PLANCK FUNCTIONS
  REAL (KIND=8) pbsur(kdlon, ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
  REAL (KIND=8) pbsuin(kdlon) ! SURFACE PLANCK FUNCTION
  REAL (KIND=8) pbtop(kdlon, ninter) ! T.O.A. SPECTRAL PLANCK FUNCTION
  REAL (KIND=8) pdbsl(kdlon, ninter, kflev*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
  REAL (KIND=8) pemis(kdlon) ! SURFACE EMISSIVITY
  REAL (KIND=8) ppmb(kdlon, kflev+1) ! HALF-LEVEL PRESSURE (MB)
  REAL (KIND=8) ptave(kdlon, kflev) ! TEMPERATURE
  REAL (KIND=8) pga(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS
  REAL (KIND=8) pgb(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS
  REAL (KIND=8) pgasur(kdlon, 8, 2) ! PADE APPROXIMANTS
  REAL (KIND=8) pgbsur(kdlon, 8, 2) ! PADE APPROXIMANTS
  REAL (KIND=8) pgatop(kdlon, 8, 2) ! PADE APPROXIMANTS
  REAL (KIND=8) pgbtop(kdlon, 8, 2) ! PADE APPROXIMANTS

  REAL (KIND=8) pcntrb(kdlon, kflev+1, kflev+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
  REAL (KIND=8) pcts(kdlon, kflev) ! COOLING-TO-SPACE TERM
  REAL (KIND=8) pfluc(kdlon, 2, kflev+1) ! CLEAR-SKY RADIATIVE FLUXES
  ! -----------------------------------------------------------------------
  ! LOCAL VARIABLES:
  REAL (KIND=8) zadjd(kdlon, kflev+1)
  REAL (KIND=8) zadju(kdlon, kflev+1)
  REAL (KIND=8) zdbdt(kdlon, ninter, kflev)
  REAL (KIND=8) zdisd(kdlon, kflev+1)
  REAL (KIND=8) zdisu(kdlon, kflev+1)

  INTEGER jk, jl
  ! -----------------------------------------------------------------------

  DO jk = 1, kflev + 1
    DO jl = 1, kdlon
      zadjd(jl, jk) = 0.
      zadju(jl, jk) = 0.
      zdisd(jl, jk) = 0.
      zdisu(jl, jk) = 0.
    END DO
  END DO

  DO jk = 1, kflev
    DO jl = 1, kdlon
      pcts(jl, jk) = 0.
    END DO
  END DO

  ! * CONTRIBUTION FROM ADJACENT LAYERS

  CALL lwvn_lmdar4(kuaer, ktraer, pabcu, pdbsl, pga, pgb, zadjd, zadju, &
    pcntrb, zdbdt)
  ! * CONTRIBUTION FROM DISTANT LAYERS

  CALL lwvd_lmdar4(kuaer, ktraer, pabcu, zdbdt, pga, pgb, pcntrb, zdisd, &
    zdisu)

  ! * EXCHANGE WITH THE BOUNDARIES

  CALL lwvb_lmdar4(kuaer, ktraer, klim, pabcu, zadjd, zadju, pb, pbint, &
    pbsuin, pbsur, pbtop, zdisd, zdisu, pemis, ppmb, pga, pgb, pgasur, &
    pgbsur, pgatop, pgbtop, pcts, pfluc)


  RETURN
END SUBROUTINE lwv_lmdar4
SUBROUTINE lwvb_lmdar4(kuaer, ktraer, klim, pabcu, padjd, padju, pb, pbint, &
    pbsui, pbsur, pbtop, pdisd, pdisu, pemis, ppmb, pga, pgb, pgasur, pgbsur, &
    pgatop, pgbtop, pcts, pfluc)
  USE dimphy
  IMPLICIT NONE
  include "raddimlw.h"
  include "radopt.h"

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! INTRODUCES THE EFFECTS OF THE BOUNDARIES IN THE VERTICAL
  ! INTEGRATION

  ! METHOD.
  ! -------

  ! 1. COMPUTES THE ENERGY EXCHANGE WITH TOP AND SURFACE OF THE
  ! ATMOSPHERE
  ! 2. COMPUTES THE COOLING-TO-SPACE AND HEATING-FROM-GROUND
  ! TERMS FOR THE APPROXIMATE COOLING RATE ABOVE 10 HPA
  ! 3. ADDS UP ALL CONTRIBUTIONS TO GET THE CLEAR-SKY FLUXES

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! Voigt lines (loop 2413 to 2427)  - JJM & PhD - 01/96
  ! -----------------------------------------------------------------------

  ! *       0.1   ARGUMENTS
  ! ---------

  INTEGER kuaer, ktraer, klim

  REAL (KIND=8) pabcu(kdlon, nua, 3*kflev+1) ! ABSORBER AMOUNTS
  REAL (KIND=8) padjd(kdlon, kflev+1) ! CONTRIBUTION BY ADJACENT LAYERS
  REAL (KIND=8) padju(kdlon, kflev+1) ! CONTRIBUTION BY ADJACENT LAYERS
  REAL (KIND=8) pb(kdlon, ninter, kflev+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
  REAL (KIND=8) pbint(kdlon, kflev+1) ! HALF-LEVEL PLANCK FUNCTIONS
  REAL (KIND=8) pbsur(kdlon, ninter) ! SPECTRAL SURFACE PLANCK FUNCTION
  REAL (KIND=8) pbsui(kdlon) ! SURFACE PLANCK FUNCTION
  REAL (KIND=8) pbtop(kdlon, ninter) ! SPECTRAL T.O.A. PLANCK FUNCTION
  REAL (KIND=8) pdisd(kdlon, kflev+1) ! CONTRIBUTION BY DISTANT LAYERS
  REAL (KIND=8) pdisu(kdlon, kflev+1) ! CONTRIBUTION BY DISTANT LAYERS
  REAL (KIND=8) pemis(kdlon) ! SURFACE EMISSIVITY
  REAL (KIND=8) ppmb(kdlon, kflev+1) ! PRESSURE MB
  REAL (KIND=8) pga(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS
  REAL (KIND=8) pgb(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS
  REAL (KIND=8) pgasur(kdlon, 8, 2) ! SURFACE PADE APPROXIMANTS
  REAL (KIND=8) pgbsur(kdlon, 8, 2) ! SURFACE PADE APPROXIMANTS
  REAL (KIND=8) pgatop(kdlon, 8, 2) ! T.O.A. PADE APPROXIMANTS
  REAL (KIND=8) pgbtop(kdlon, 8, 2) ! T.O.A. PADE APPROXIMANTS

  REAL (KIND=8) pfluc(kdlon, 2, kflev+1) ! CLEAR-SKY RADIATIVE FLUXES
  REAL (KIND=8) pcts(kdlon, kflev) ! COOLING-TO-SPACE TERM

  ! * LOCAL VARIABLES:

  REAL (KIND=8) zbgnd(kdlon)
  REAL (KIND=8) zfd(kdlon)
  REAL (KIND=8) zfn10(kdlon)
  REAL (KIND=8) zfu(kdlon)
  REAL (KIND=8) ztt(kdlon, ntra)
  REAL (KIND=8) ztt1(kdlon, ntra)
  REAL (KIND=8) ztt2(kdlon, ntra)
  REAL (KIND=8) zuu(kdlon, nua)
  REAL (KIND=8) zcnsol(kdlon)
  REAL (KIND=8) zcntop(kdlon)

  INTEGER jk, jl, ja
  INTEGER jstra, jstru
  INTEGER ind1, ind2, ind3, ind4, in, jlim
  REAL (KIND=8) zctstr

  ! -----------------------------------------------------------------------

  ! *         1.    INITIALIZATION
  ! --------------



  ! *         1.2     INITIALIZE TRANSMISSION FUNCTIONS
  ! ---------------------------------


  DO ja = 1, ntra
    DO jl = 1, kdlon
      ztt(jl, ja) = 1.0
      ztt1(jl, ja) = 1.0
      ztt2(jl, ja) = 1.0
    END DO
  END DO

  DO ja = 1, nua
    DO jl = 1, kdlon
      zuu(jl, ja) = 1.0
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.      VERTICAL INTEGRATION
  ! --------------------


  ind1 = 0
  ind3 = 0
  ind4 = 1
  ind2 = 1

  ! *         2.3     EXCHANGE WITH TOP OF THE ATMOSPHERE
  ! -----------------------------------


  DO jk = 1, kflev
    in = (jk-1)*ng1p1 + 1

    DO ja = 1, kuaer
      DO jl = 1, kdlon
        zuu(jl, ja) = pabcu(jl, ja, in)
      END DO
    END DO


    CALL lwtt_lmdar4(pgatop(1,1,1), pgbtop(1,1,1), zuu, ztt)

    DO jl = 1, kdlon
      zcntop(jl) = pbtop(jl, 1)*ztt(jl, 1)*ztt(jl, 10) + &
        pbtop(jl, 2)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
        pbtop(jl, 3)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
        pbtop(jl, 4)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
        pbtop(jl, 5)*ztt(jl, 3)*ztt(jl, 14) + pbtop(jl, 6)*ztt(jl, 6)*ztt(jl, &
        15)
      zfd(jl) = zcntop(jl) - pbint(jl, jk) - pdisd(jl, jk) - padjd(jl, jk)
      pfluc(jl, 2, jk) = zfd(jl)
    END DO

  END DO

  jk = kflev + 1
  in = (jk-1)*ng1p1 + 1

  DO jl = 1, kdlon
    zcntop(jl) = pbtop(jl, 1) + pbtop(jl, 2) + pbtop(jl, 3) + pbtop(jl, 4) + &
      pbtop(jl, 5) + pbtop(jl, 6)
    zfd(jl) = zcntop(jl) - pbint(jl, jk) - pdisd(jl, jk) - padjd(jl, jk)
    pfluc(jl, 2, jk) = zfd(jl)
  END DO

  ! *         2.4     COOLING-TO-SPACE OF LAYERS ABOVE 10 HPA
  ! ---------------------------------------



  ! *         2.4.1   INITIALIZATION
  ! --------------


  jlim = kflev

  IF (.NOT. levoigt) THEN
    DO jk = kflev, 1, -1
      IF (ppmb(1,jk)<10.0) THEN
        jlim = jk
      END IF
    END DO
  END IF
  klim = jlim

  IF (.NOT. levoigt) THEN
    DO ja = 1, ktraer
      DO jl = 1, kdlon
        ztt1(jl, ja) = 1.0
      END DO
    END DO

    ! *         2.4.2   LOOP OVER LAYERS ABOVE 10 HPA
    ! -----------------------------


    DO jstra = kflev, jlim, -1
      jstru = (jstra-1)*ng1p1 + 1

      DO ja = 1, kuaer
        DO jl = 1, kdlon
          zuu(jl, ja) = pabcu(jl, ja, jstru)
        END DO
      END DO


      CALL lwtt_lmdar4(pga(1,1,1,jstra), pgb(1,1,1,jstra), zuu, ztt)

      DO jl = 1, kdlon
        zctstr = (pb(jl,1,jstra)+pb(jl,1,jstra+1))* &
          (ztt1(jl,1)*ztt1(jl,10)-ztt(jl,1)*ztt(jl,10)) + &
          (pb(jl,2,jstra)+pb(jl,2,jstra+1))*(ztt1(jl,2)*ztt1(jl,7)*ztt1(jl,11 &
          )-ztt(jl,2)*ztt(jl,7)*ztt(jl,11)) + (pb(jl,3,jstra)+pb(jl,3,jstra+1 &
          ))*(ztt1(jl,4)*ztt1(jl,8)*ztt1(jl,12)-ztt(jl,4)*ztt(jl,8)*ztt(jl,12 &
          )) + (pb(jl,4,jstra)+pb(jl,4,jstra+1))*(ztt1(jl,5)*ztt1(jl,9)*ztt1( &
          jl,13)-ztt(jl,5)*ztt(jl,9)*ztt(jl,13)) + (pb(jl,5,jstra)+pb(jl,5, &
          jstra+1))*(ztt1(jl,3)*ztt1(jl,14)-ztt(jl,3)*ztt(jl,14)) + &
          (pb(jl,6,jstra)+pb(jl,6,jstra+1))*(ztt1(jl,6)*ztt1(jl,15)-ztt(jl,6) &
          *ztt(jl,15))
        pcts(jl, jstra) = zctstr*0.5
      END DO
      DO ja = 1, ktraer
        DO jl = 1, kdlon
          ztt1(jl, ja) = ztt(jl, ja)
        END DO
      END DO
    END DO
  END IF
  ! Mise a zero de securite pour PCTS en cas de LEVOIGT
  IF (levoigt) THEN
    DO jstra = 1, kflev
      DO jl = 1, kdlon
        pcts(jl, jstra) = 0.
      END DO
    END DO
  END IF

  ! *         2.5     EXCHANGE WITH LOWER LIMIT
  ! -------------------------


  DO jl = 1, kdlon
    zbgnd(jl) = pbsui(jl)*pemis(jl) - (1.-pemis(jl))*pfluc(jl, 2, 1) - &
      pbint(jl, 1)
  END DO

  jk = 1
  in = (jk-1)*ng1p1 + 1

  DO jl = 1, kdlon
    zcnsol(jl) = pbsur(jl, 1) + pbsur(jl, 2) + pbsur(jl, 3) + pbsur(jl, 4) + &
      pbsur(jl, 5) + pbsur(jl, 6)
    zcnsol(jl) = zcnsol(jl)*zbgnd(jl)/pbsui(jl)
    zfu(jl) = zcnsol(jl) + pbint(jl, jk) - pdisu(jl, jk) - padju(jl, jk)
    pfluc(jl, 1, jk) = zfu(jl)
  END DO

  DO jk = 2, kflev + 1
    in = (jk-1)*ng1p1 + 1


    DO ja = 1, kuaer
      DO jl = 1, kdlon
        zuu(jl, ja) = pabcu(jl, ja, 1) - pabcu(jl, ja, in)
      END DO
    END DO


    CALL lwtt_lmdar4(pgasur(1,1,1), pgbsur(1,1,1), zuu, ztt)

    DO jl = 1, kdlon
      zcnsol(jl) = pbsur(jl, 1)*ztt(jl, 1)*ztt(jl, 10) + &
        pbsur(jl, 2)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
        pbsur(jl, 3)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
        pbsur(jl, 4)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
        pbsur(jl, 5)*ztt(jl, 3)*ztt(jl, 14) + pbsur(jl, 6)*ztt(jl, 6)*ztt(jl, &
        15)
      zcnsol(jl) = zcnsol(jl)*zbgnd(jl)/pbsui(jl)
      zfu(jl) = zcnsol(jl) + pbint(jl, jk) - pdisu(jl, jk) - padju(jl, jk)
      pfluc(jl, 1, jk) = zfu(jl)
    END DO


  END DO

  ! *         2.7     CLEAR-SKY FLUXES
  ! ----------------


  IF (.NOT. levoigt) THEN
    DO jl = 1, kdlon
      zfn10(jl) = pfluc(jl, 1, jlim) + pfluc(jl, 2, jlim)
    END DO
    DO jk = jlim + 1, kflev + 1
      DO jl = 1, kdlon
        zfn10(jl) = zfn10(jl) + pcts(jl, jk-1)
        pfluc(jl, 1, jk) = zfn10(jl)
        pfluc(jl, 2, jk) = 0.
      END DO
    END DO
  END IF

  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE lwvb_lmdar4
SUBROUTINE lwvd_lmdar4(kuaer, ktraer, pabcu, pdbdt, pga, pgb, pcntrb, pdisd, &
    pdisu)
  USE dimphy
  IMPLICIT NONE
  include "raddimlw.h"

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! CARRIES OUT THE VERTICAL INTEGRATION ON THE DISTANT LAYERS

  ! METHOD.
  ! -------

  ! 1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
  ! CONTRIBUTIONS OF THE DISTANT LAYERS USING TRAPEZOIDAL RULE

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! -----------------------------------------------------------------------
  ! * ARGUMENTS:

  INTEGER kuaer, ktraer

  REAL (KIND=8) pabcu(kdlon, nua, 3*kflev+1) ! ABSORBER AMOUNTS
  REAL (KIND=8) pdbdt(kdlon, ninter, kflev) ! LAYER PLANCK FUNCTION GRADIENT
  REAL (KIND=8) pga(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS
  REAL (KIND=8) pgb(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS

  REAL (KIND=8) pcntrb(kdlon, kflev+1, kflev+1) ! ENERGY EXCHANGE MATRIX
  REAL (KIND=8) pdisd(kdlon, kflev+1) !  CONTRIBUTION BY DISTANT LAYERS
  REAL (KIND=8) pdisu(kdlon, kflev+1) !  CONTRIBUTION BY DISTANT LAYERS

  ! * LOCAL VARIABLES:

  REAL (KIND=8) zglayd(kdlon)
  REAL (KIND=8) zglayu(kdlon)
  REAL (KIND=8) ztt(kdlon, ntra)
  REAL (KIND=8) ztt1(kdlon, ntra)
  REAL (KIND=8) ztt2(kdlon, ntra)

  INTEGER jl, jk, ja, ikp1, ikn, ikd1, jkj, ikd2
  INTEGER ikjp1, ikm1, ikj, jlk, iku1, ijkl, iku2
  INTEGER ind1, ind2, ind3, ind4, itt
  REAL (KIND=8) zww, zdzxdg, zdzxmg

  ! *         1.    INITIALIZATION
  ! --------------


  ! *         1.1     INITIALIZE LAYER CONTRIBUTIONS
  ! ------------------------------


  DO jk = 1, kflev + 1
    DO jl = 1, kdlon
      pdisd(jl, jk) = 0.
      pdisu(jl, jk) = 0.
    END DO
  END DO

  ! *         1.2     INITIALIZE TRANSMISSION FUNCTIONS
  ! ---------------------------------



  DO ja = 1, ntra
    DO jl = 1, kdlon
      ztt(jl, ja) = 1.0
      ztt1(jl, ja) = 1.0
      ztt2(jl, ja) = 1.0
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.      VERTICAL INTEGRATION
  ! --------------------


  ind1 = 0
  ind3 = 0
  ind4 = 1
  ind2 = 1

  ! *         2.2     CONTRIBUTION FROM DISTANT LAYERS
  ! ---------------------------------



  ! *         2.2.1   DISTANT AND ABOVE LAYERS
  ! ------------------------




  ! *         2.2.2   FIRST UPPER LEVEL
  ! -----------------


  DO jk = 1, kflev - 1
    ikp1 = jk + 1
    ikn = (jk-1)*ng1p1 + 1
    ikd1 = jk*ng1p1 + 1

    CALL lwttm_lmdar4(pga(1,1,1,jk), pgb(1,1,1,jk), pabcu(1,1,ikn), &
      pabcu(1,1,ikd1), ztt1)

    ! *         2.2.3   HIGHER UP
    ! ---------


    itt = 1
    DO jkj = ikp1, kflev
      IF (itt==1) THEN
        itt = 2
      ELSE
        itt = 1
      END IF
      ikjp1 = jkj + 1
      ikd2 = jkj*ng1p1 + 1

      IF (itt==1) THEN
        CALL lwttm_lmdar4(pga(1,1,1,jkj), pgb(1,1,1,jkj), pabcu(1,1,ikn), &
          pabcu(1,1,ikd2), ztt1)
      ELSE
        CALL lwttm_lmdar4(pga(1,1,1,jkj), pgb(1,1,1,jkj), pabcu(1,1,ikn), &
          pabcu(1,1,ikd2), ztt2)
      END IF

      DO ja = 1, ktraer
        DO jl = 1, kdlon
          ztt(jl, ja) = (ztt1(jl,ja)+ztt2(jl,ja))*0.5
        END DO
      END DO

      DO jl = 1, kdlon
        zww = pdbdt(jl, 1, jkj)*ztt(jl, 1)*ztt(jl, 10) + &
          pdbdt(jl, 2, jkj)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
          pdbdt(jl, 3, jkj)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
          pdbdt(jl, 4, jkj)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
          pdbdt(jl, 5, jkj)*ztt(jl, 3)*ztt(jl, 14) + &
          pdbdt(jl, 6, jkj)*ztt(jl, 6)*ztt(jl, 15)
        zglayd(jl) = zww
        zdzxdg = zglayd(jl)
        pdisd(jl, jk) = pdisd(jl, jk) + zdzxdg
        pcntrb(jl, jk, ikjp1) = zdzxdg
      END DO


    END DO
  END DO

  ! *         2.2.4   DISTANT AND BELOW LAYERS
  ! ------------------------




  ! *         2.2.5   FIRST LOWER LEVEL
  ! -----------------


  DO jk = 3, kflev + 1
    ikn = (jk-1)*ng1p1 + 1
    ikm1 = jk - 1
    ikj = jk - 2
    iku1 = ikj*ng1p1 + 1


    CALL lwttm_lmdar4(pga(1,1,1,ikj), pgb(1,1,1,ikj), pabcu(1,1,iku1), &
      pabcu(1,1,ikn), ztt1)

    ! *         2.2.6   DOWN BELOW
    ! ----------


    itt = 1
    DO jlk = 1, ikj
      IF (itt==1) THEN
        itt = 2
      ELSE
        itt = 1
      END IF
      ijkl = ikm1 - jlk
      iku2 = (ijkl-1)*ng1p1 + 1


      IF (itt==1) THEN
        CALL lwttm_lmdar4(pga(1,1,1,ijkl), pgb(1,1,1,ijkl), pabcu(1,1,iku2), &
          pabcu(1,1,ikn), ztt1)
      ELSE
        CALL lwttm_lmdar4(pga(1,1,1,ijkl), pgb(1,1,1,ijkl), pabcu(1,1,iku2), &
          pabcu(1,1,ikn), ztt2)
      END IF

      DO ja = 1, ktraer
        DO jl = 1, kdlon
          ztt(jl, ja) = (ztt1(jl,ja)+ztt2(jl,ja))*0.5
        END DO
      END DO

      DO jl = 1, kdlon
        zww = pdbdt(jl, 1, ijkl)*ztt(jl, 1)*ztt(jl, 10) + &
          pdbdt(jl, 2, ijkl)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
          pdbdt(jl, 3, ijkl)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
          pdbdt(jl, 4, ijkl)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
          pdbdt(jl, 5, ijkl)*ztt(jl, 3)*ztt(jl, 14) + &
          pdbdt(jl, 6, ijkl)*ztt(jl, 6)*ztt(jl, 15)
        zglayu(jl) = zww
        zdzxmg = zglayu(jl)
        pdisu(jl, jk) = pdisu(jl, jk) + zdzxmg
        pcntrb(jl, jk, ijkl) = zdzxmg
      END DO


    END DO
  END DO

  RETURN
END SUBROUTINE lwvd_lmdar4
SUBROUTINE lwvn_lmdar4(kuaer, ktraer, pabcu, pdbsl, pga, pgb, padjd, padju, &
    pcntrb, pdbdt)
  USE dimphy
  USE radiation_ar4_param, ONLY: wg1
  IMPLICIT NONE
  include "raddimlw.h"

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! CARRIES OUT THE VERTICAL INTEGRATION ON NEARBY LAYERS
  ! TO GIVE LONGWAVE FLUXES OR RADIANCES

  ! METHOD.
  ! -------

  ! 1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
  ! CONTRIBUTIONS OF THE ADJACENT LAYERS USING A GAUSSIAN QUADRATURE

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! -----------------------------------------------------------------------

  ! * ARGUMENTS:

  INTEGER kuaer, ktraer

  REAL (KIND=8) pabcu(kdlon, nua, 3*kflev+1) ! ABSORBER AMOUNTS
  REAL (KIND=8) pdbsl(kdlon, ninter, kflev*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
  REAL (KIND=8) pga(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS
  REAL (KIND=8) pgb(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS

  REAL (KIND=8) padjd(kdlon, kflev+1) ! CONTRIBUTION OF ADJACENT LAYERS
  REAL (KIND=8) padju(kdlon, kflev+1) ! CONTRIBUTION OF ADJACENT LAYERS
  REAL (KIND=8) pcntrb(kdlon, kflev+1, kflev+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
  REAL (KIND=8) pdbdt(kdlon, ninter, kflev) !  LAYER PLANCK FUNCTION GRADIENT

  ! * LOCAL ARRAYS:

  REAL (KIND=8) zglayd(kdlon)
  REAL (KIND=8) zglayu(kdlon)
  REAL (KIND=8) ztt(kdlon, ntra)
  REAL (KIND=8) ztt1(kdlon, ntra)
  REAL (KIND=8) ztt2(kdlon, ntra)
  REAL (KIND=8) zuu(kdlon, nua)

  INTEGER jk, jl, ja, im12, ind, inu, ixu, jg
  INTEGER ixd, ibs, idd, imu, jk1, jk2, jnu
  REAL (KIND=8) zwtr

  ! -----------------------------------------------------------------------

  ! *         1.    INITIALIZATION
  ! --------------


  ! *         1.1     INITIALIZE LAYER CONTRIBUTIONS
  ! ------------------------------


  DO jk = 1, kflev + 1
    DO jl = 1, kdlon
      padjd(jl, jk) = 0.
      padju(jl, jk) = 0.
    END DO
  END DO

  ! *         1.2     INITIALIZE TRANSMISSION FUNCTIONS
  ! ---------------------------------


  DO ja = 1, ntra
    DO jl = 1, kdlon
      ztt(jl, ja) = 1.0
      ztt1(jl, ja) = 1.0
      ztt2(jl, ja) = 1.0
    END DO
  END DO

  DO ja = 1, nua
    DO jl = 1, kdlon
      zuu(jl, ja) = 0.
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.      VERTICAL INTEGRATION
  ! --------------------



  ! *         2.1     CONTRIBUTION FROM ADJACENT LAYERS
  ! ---------------------------------


  DO jk = 1, kflev
    ! *         2.1.1   DOWNWARD LAYERS
    ! ---------------


    im12 = 2*(jk-1)
    ind = (jk-1)*ng1p1 + 1
    ixd = ind
    inu = jk*ng1p1 + 1
    ixu = ind

    DO jl = 1, kdlon
      zglayd(jl) = 0.
      zglayu(jl) = 0.
    END DO

    DO jg = 1, ng1
      ibs = im12 + jg
      idd = ixd + jg
      DO ja = 1, kuaer
        DO jl = 1, kdlon
          zuu(jl, ja) = pabcu(jl, ja, ind) - pabcu(jl, ja, idd)
        END DO
      END DO


      CALL lwtt_lmdar4(pga(1,1,1,jk), pgb(1,1,1,jk), zuu, ztt)

      DO jl = 1, kdlon
        zwtr = pdbsl(jl, 1, ibs)*ztt(jl, 1)*ztt(jl, 10) + &
          pdbsl(jl, 2, ibs)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
          pdbsl(jl, 3, ibs)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
          pdbsl(jl, 4, ibs)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
          pdbsl(jl, 5, ibs)*ztt(jl, 3)*ztt(jl, 14) + &
          pdbsl(jl, 6, ibs)*ztt(jl, 6)*ztt(jl, 15)
        zglayd(jl) = zglayd(jl) + zwtr*wg1(jg)
      END DO

      ! *         2.1.2   DOWNWARD LAYERS
      ! ---------------


      imu = ixu + jg
      DO ja = 1, kuaer
        DO jl = 1, kdlon
          zuu(jl, ja) = pabcu(jl, ja, imu) - pabcu(jl, ja, inu)
        END DO
      END DO


      CALL lwtt_lmdar4(pga(1,1,1,jk), pgb(1,1,1,jk), zuu, ztt)

      DO jl = 1, kdlon
        zwtr = pdbsl(jl, 1, ibs)*ztt(jl, 1)*ztt(jl, 10) + &
          pdbsl(jl, 2, ibs)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
          pdbsl(jl, 3, ibs)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
          pdbsl(jl, 4, ibs)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
          pdbsl(jl, 5, ibs)*ztt(jl, 3)*ztt(jl, 14) + &
          pdbsl(jl, 6, ibs)*ztt(jl, 6)*ztt(jl, 15)
        zglayu(jl) = zglayu(jl) + zwtr*wg1(jg)
      END DO

    END DO

    DO jl = 1, kdlon
      padjd(jl, jk) = zglayd(jl)
      pcntrb(jl, jk, jk+1) = zglayd(jl)
      padju(jl, jk+1) = zglayu(jl)
      pcntrb(jl, jk+1, jk) = zglayu(jl)
      pcntrb(jl, jk, jk) = 0.0
    END DO

  END DO

  DO jk = 1, kflev
    jk2 = 2*jk
    jk1 = jk2 - 1
    DO jnu = 1, ninter
      DO jl = 1, kdlon
        pdbdt(jl, jnu, jk) = pdbsl(jl, jnu, jk1) + pdbsl(jl, jnu, jk2)
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE lwvn_lmdar4
SUBROUTINE lwtt_lmdar4(pga, pgb, puu, ptt)
  USE dimphy
  IMPLICIT NONE
  include "raddimlw.h"

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
  ! ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
  ! INTERVALS.

  ! METHOD.
  ! -------

  ! 1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
  ! COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
  ! 2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
  ! 3. TRANSMISSION BY H2O CONTINUUM AND AEROSOLS FOLLOW AN
  ! A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 88-12-15

  ! -----------------------------------------------------------------------
  REAL (KIND=8) o1h, o2h
  PARAMETER (o1h=2230.)
  PARAMETER (o2h=100.)
  REAL (KIND=8) rpialf0
  PARAMETER (rpialf0=2.0)

  ! * ARGUMENTS:

  REAL (KIND=8) puu(kdlon, nua)
  REAL (KIND=8) ptt(kdlon, ntra)
  REAL (KIND=8) pga(kdlon, 8, 2)
  REAL (KIND=8) pgb(kdlon, 8, 2)

  ! * LOCAL VARIABLES:

  REAL (KIND=8) zz, zxd, zxn
  REAL (KIND=8) zpu, zpu10, zpu11, zpu12, zpu13
  REAL (KIND=8) zeu, zeu10, zeu11, zeu12, zeu13
  REAL (KIND=8) zx, zy, zsq1, zsq2, zvxy, zuxy
  REAL (KIND=8) zaercn, zto1, zto2, zxch4, zych4, zxn2o, zyn2o
  REAL (KIND=8) zsqn21, zodn21, zsqh42, zodh42
  REAL (KIND=8) zsqh41, zodh41, zsqn22, zodn22, zttf11, zttf12
  REAL (KIND=8) zuu11, zuu12, za11, za12
  INTEGER jl, ja

  ! ------------------------------------------------------------------

  ! *         1.     HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
  ! -----------------------------------------------



  ! cdir collapse
  DO ja = 1, 8
    DO jl = 1, kdlon
      zz = sqrt(puu(jl,ja))
      ! ZXD(JL,1)=PGB( JL, 1,1) + ZZ(JL, 1)*(PGB( JL, 1,2) + ZZ(JL, 1))
      ! ZXN(JL,1)=PGA( JL, 1,1) + ZZ(JL, 1)*(PGA( JL, 1,2) )
      ! PTT(JL,1)=ZXN(JL,1)/ZXD(JL,1)
      zxd = pgb(jl, ja, 1) + zz*(pgb(jl,ja,2)+zz)
      zxn = pga(jl, ja, 1) + zz*(pga(jl,ja,2))
      ptt(jl, ja) = zxn/zxd
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.     CONTINUUM, OZONE AND AEROSOL TRANSMISSION FUNCTIONS
  ! ---------------------------------------------------


  DO jl = 1, kdlon
    ptt(jl, 9) = ptt(jl, 8)

    ! -  CONTINUUM ABSORPTION: E- AND P-TYPE

    zpu = 0.002*puu(jl, 10)
    zpu10 = 112.*zpu
    zpu11 = 6.25*zpu
    zpu12 = 5.00*zpu
    zpu13 = 80.0*zpu
    zeu = puu(jl, 11)
    zeu10 = 12.*zeu
    zeu11 = 6.25*zeu
    zeu12 = 5.00*zeu
    zeu13 = 80.0*zeu

    ! -  OZONE ABSORPTION

    zx = puu(jl, 12)
    zy = puu(jl, 13)
    zuxy = 4.*zx*zx/(rpialf0*zy)
    zsq1 = sqrt(1.+o1h*zuxy) - 1.
    zsq2 = sqrt(1.+o2h*zuxy) - 1.
    zvxy = rpialf0*zy/(2.*zx)
    zaercn = puu(jl, 17) + zeu12 + zpu12
    zto1 = exp(-zvxy*zsq1-zaercn)
    zto2 = exp(-zvxy*zsq2-zaercn)

    ! -- TRACE GASES (CH4, N2O, CFC-11, CFC-12)

    ! * CH4 IN INTERVAL 800-970 + 1110-1250 CM-1

    ! NEXOTIC=1
    ! IF (NEXOTIC.EQ.1) THEN
    zxch4 = puu(jl, 19)
    zych4 = puu(jl, 20)
    zuxy = 4.*zxch4*zxch4/(0.103*zych4)
    zsqh41 = sqrt(1.+33.7*zuxy) - 1.
    zvxy = 0.103*zych4/(2.*zxch4)
    zodh41 = zvxy*zsqh41

    ! * N2O IN INTERVAL 800-970 + 1110-1250 CM-1

    zxn2o = puu(jl, 21)
    zyn2o = puu(jl, 22)
    zuxy = 4.*zxn2o*zxn2o/(0.416*zyn2o)
    zsqn21 = sqrt(1.+21.3*zuxy) - 1.
    zvxy = 0.416*zyn2o/(2.*zxn2o)
    zodn21 = zvxy*zsqn21

    ! * CH4 IN INTERVAL 1250-1450 + 1880-2820 CM-1

    zuxy = 4.*zxch4*zxch4/(0.113*zych4)
    zsqh42 = sqrt(1.+400.*zuxy) - 1.
    zvxy = 0.113*zych4/(2.*zxch4)
    zodh42 = zvxy*zsqh42

    ! * N2O IN INTERVAL 1250-1450 + 1880-2820 CM-1

    zuxy = 4.*zxn2o*zxn2o/(0.197*zyn2o)
    zsqn22 = sqrt(1.+2000.*zuxy) - 1.
    zvxy = 0.197*zyn2o/(2.*zxn2o)
    zodn22 = zvxy*zsqn22

    ! * CFC-11 IN INTERVAL 800-970 + 1110-1250 CM-1

    za11 = 2.*puu(jl, 23)*4.404E+05
    zttf11 = 1. - za11*0.003225

    ! * CFC-12 IN INTERVAL 800-970 + 1110-1250 CM-1

    za12 = 2.*puu(jl, 24)*6.7435E+05
    zttf12 = 1. - za12*0.003225

    zuu11 = -puu(jl, 15) - zeu10 - zpu10
    zuu12 = -puu(jl, 16) - zeu11 - zpu11 - zodh41 - zodn21
    ptt(jl, 10) = exp(-puu(jl,14))
    ptt(jl, 11) = exp(zuu11)
    ptt(jl, 12) = exp(zuu12)*zttf11*zttf12
    ptt(jl, 13) = 0.7554*zto1 + 0.2446*zto2
    ptt(jl, 14) = ptt(jl, 10)*exp(-zeu13-zpu13)
    ptt(jl, 15) = exp(-puu(jl,14)-zodh42-zodn22)
  END DO

  RETURN
END SUBROUTINE lwtt_lmdar4
SUBROUTINE lwttm_lmdar4(pga, pgb, puu1, puu2, ptt)
  USE dimphy
  IMPLICIT NONE
  include "raddimlw.h"

  ! ------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
  ! ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
  ! INTERVALS.

  ! METHOD.
  ! -------

  ! 1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
  ! COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
  ! 2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
  ! 3. TRANSMISSION BY H2O CONTINUUM AND AEROSOLS FOLLOW AN
  ! A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 88-12-15

  ! -----------------------------------------------------------------------
  REAL (KIND=8) o1h, o2h
  PARAMETER (o1h=2230.)
  PARAMETER (o2h=100.)
  REAL (KIND=8) rpialf0
  PARAMETER (rpialf0=2.0)

  ! * ARGUMENTS:

  REAL (KIND=8) pga(kdlon, 8, 2) ! PADE APPROXIMANTS
  REAL (KIND=8) pgb(kdlon, 8, 2) ! PADE APPROXIMANTS
  REAL (KIND=8) puu1(kdlon, nua) ! ABSORBER AMOUNTS FROM TOP TO LEVEL 1
  REAL (KIND=8) puu2(kdlon, nua) ! ABSORBER AMOUNTS FROM TOP TO LEVEL 2
  REAL (KIND=8) ptt(kdlon, ntra) ! TRANSMISSION FUNCTIONS

  ! * LOCAL VARIABLES:

  INTEGER ja, jl
  REAL (KIND=8) zz, zxd, zxn
  REAL (KIND=8) zpu, zpu10, zpu11, zpu12, zpu13
  REAL (KIND=8) zeu, zeu10, zeu11, zeu12, zeu13
  REAL (KIND=8) zx, zy, zuxy, zsq1, zsq2, zvxy, zaercn, zto1, zto2
  REAL (KIND=8) zxch4, zych4, zsqh41, zodh41
  REAL (KIND=8) zxn2o, zyn2o, zsqn21, zodn21, zsqh42, zodh42
  REAL (KIND=8) zsqn22, zodn22, za11, zttf11, za12, zttf12
  REAL (KIND=8) zuu11, zuu12

  ! ------------------------------------------------------------------

  ! *         1.     HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
  ! -----------------------------------------------




  ! CDIR ON_ADB(PUU1)
  ! CDIR ON_ADB(PUU2)
  ! CDIR COLLAPSE
  DO ja = 1, 8
    DO jl = 1, kdlon
      zz = sqrt(puu1(jl,ja)-puu2(jl,ja))
      zxd = pgb(jl, ja, 1) + zz*(pgb(jl,ja,2)+zz)
      zxn = pga(jl, ja, 1) + zz*(pga(jl,ja,2))
      ptt(jl, ja) = zxn/zxd
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.     CONTINUUM, OZONE AND AEROSOL TRANSMISSION FUNCTIONS
  ! ---------------------------------------------------


  DO jl = 1, kdlon
    ptt(jl, 9) = ptt(jl, 8)

    ! -  CONTINUUM ABSORPTION: E- AND P-TYPE

    zpu = 0.002*(puu1(jl,10)-puu2(jl,10))
    zpu10 = 112.*zpu
    zpu11 = 6.25*zpu
    zpu12 = 5.00*zpu
    zpu13 = 80.0*zpu
    zeu = (puu1(jl,11)-puu2(jl,11))
    zeu10 = 12.*zeu
    zeu11 = 6.25*zeu
    zeu12 = 5.00*zeu
    zeu13 = 80.0*zeu

    ! -  OZONE ABSORPTION

    zx = (puu1(jl,12)-puu2(jl,12))
    zy = (puu1(jl,13)-puu2(jl,13))
    zuxy = 4.*zx*zx/(rpialf0*zy)
    zsq1 = sqrt(1.+o1h*zuxy) - 1.
    zsq2 = sqrt(1.+o2h*zuxy) - 1.
    zvxy = rpialf0*zy/(2.*zx)
    zaercn = (puu1(jl,17)-puu2(jl,17)) + zeu12 + zpu12
    zto1 = exp(-zvxy*zsq1-zaercn)
    zto2 = exp(-zvxy*zsq2-zaercn)

    ! -- TRACE GASES (CH4, N2O, CFC-11, CFC-12)

    ! * CH4 IN INTERVAL 800-970 + 1110-1250 CM-1

    zxch4 = (puu1(jl,19)-puu2(jl,19))
    zych4 = (puu1(jl,20)-puu2(jl,20))
    zuxy = 4.*zxch4*zxch4/(0.103*zych4)
    zsqh41 = sqrt(1.+33.7*zuxy) - 1.
    zvxy = 0.103*zych4/(2.*zxch4)
    zodh41 = zvxy*zsqh41

    ! * N2O IN INTERVAL 800-970 + 1110-1250 CM-1

    zxn2o = (puu1(jl,21)-puu2(jl,21))
    zyn2o = (puu1(jl,22)-puu2(jl,22))
    zuxy = 4.*zxn2o*zxn2o/(0.416*zyn2o)
    zsqn21 = sqrt(1.+21.3*zuxy) - 1.
    zvxy = 0.416*zyn2o/(2.*zxn2o)
    zodn21 = zvxy*zsqn21

    ! * CH4 IN INTERVAL 1250-1450 + 1880-2820 CM-1

    zuxy = 4.*zxch4*zxch4/(0.113*zych4)
    zsqh42 = sqrt(1.+400.*zuxy) - 1.
    zvxy = 0.113*zych4/(2.*zxch4)
    zodh42 = zvxy*zsqh42

    ! * N2O IN INTERVAL 1250-1450 + 1880-2820 CM-1

    zuxy = 4.*zxn2o*zxn2o/(0.197*zyn2o)
    zsqn22 = sqrt(1.+2000.*zuxy) - 1.
    zvxy = 0.197*zyn2o/(2.*zxn2o)
    zodn22 = zvxy*zsqn22

    ! * CFC-11 IN INTERVAL 800-970 + 1110-1250 CM-1

    za11 = (puu1(jl,23)-puu2(jl,23))*4.404E+05
    zttf11 = 1. - za11*0.003225

    ! * CFC-12 IN INTERVAL 800-970 + 1110-1250 CM-1

    za12 = (puu1(jl,24)-puu2(jl,24))*6.7435E+05
    zttf12 = 1. - za12*0.003225

    zuu11 = -(puu1(jl,15)-puu2(jl,15)) - zeu10 - zpu10
    zuu12 = -(puu1(jl,16)-puu2(jl,16)) - zeu11 - zpu11 - zodh41 - zodn21
    ptt(jl, 10) = exp(-(puu1(jl,14)-puu2(jl,14)))
    ptt(jl, 11) = exp(zuu11)
    ptt(jl, 12) = exp(zuu12)*zttf11*zttf12
    ptt(jl, 13) = 0.7554*zto1 + 0.2446*zto2
    ptt(jl, 14) = ptt(jl, 10)*exp(-zeu13-zpu13)
    ptt(jl, 15) = exp(-(puu1(jl,14)-puu2(jl,14))-zodh42-zodn22)
  END DO

  RETURN
END SUBROUTINE lwttm_lmdar4
