!
! $Id$
!
SUBROUTINE SW_AEROAR4(PSCT, PRMU0, PFRAC, &
     PPMB, PDP, &
     PPSOL, PALBD, PALBP,&
     PTAVE, PWV, PQS, POZON, PAER,&
     PCLDSW, PTAU, POMEGA, PCG,&
     PHEAT, PHEAT0,&
     PALBPLA,PTOPSW,PSOLSW,PTOPSW0,PSOLSW0,&
     ZFSUP,ZFSDN,ZFSUP0,ZFSDN0,&
     tauaero, pizaero, cgaero,&
     PTAUA, POMEGAA,&
     PTOPSWADAERO,PSOLSWADAERO,&
     PTOPSWAD0AERO,PSOLSWAD0AERO,&
     PTOPSWAIAERO,PSOLSWAIAERO,&
     PTOPSWAERO,PTOPSW0AERO,&
     PSOLSWAERO,PSOLSW0AERO,&
     PTOPSWCFAERO,PSOLSWCFAERO,&
     ok_ade, ok_aie, flag_aerosol, flag_aerosol_strat )

  USE dimphy
  USE phys_output_mod, ONLY : swaero_diag
  USE print_control_mod, ONLY: lunout
  USE aero_mod, ONLY : naero_grp
  IMPLICIT NONE

#include "YOMCST.h"
#include "clesphys.h"
  !
  !     ------------------------------------------------------------------
  !
  !     PURPOSE.
  !     --------
  !
  !          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
  !     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).
  !
  !     METHOD.
  !     -------
  !
  !          1. COMPUTES ABSORBER AMOUNTS                 (SWU)
  !          2. COMPUTES FLUXES IN 1ST SPECTRAL INTERVAL  (SW1S)
  !          3. COMPUTES FLUXES IN 2ND SPECTRAL INTERVAL  (SW2S)
  !
  !     REFERENCE.
  !     ----------
  !
  !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
  !        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
  !
  !     AUTHOR.
  !     -------
  !        JEAN-JACQUES MORCRETTE  *ECMWF*
  !
  !     MODIFICATIONS.
  !     --------------
  !        ORIGINAL : 89-07-14
  !        1995-01-01  J.-J. MORCRETTE  Direct/Diffuse Albedo
  !        2003-11-27  J. QUAAS Introduce aerosol forcings (based on BOUCHER)
  !        2009-04     A. COZIC - C.DEANDREIS Indroduce NAT/BC/POM/DUST/SS aerosol forcing
  !        2012-09     O. BOUCHER - reorganise aerosol cases with ok_ade, ok_aie, flag_aerosol
  !     ------------------------------------------------------------------
  !
  !* ARGUMENTS:
  !
  REAL(KIND=8) PSCT  ! constante solaire (valeur conseillee: 1370)

  REAL(KIND=8) PPSOL(KDLON)        ! SURFACE PRESSURE (PA)
  REAL(KIND=8) PDP(KDLON,KFLEV)    ! LAYER THICKNESS (PA)
  REAL(KIND=8) PPMB(KDLON,KFLEV+1) ! HALF-LEVEL PRESSURE (MB)

  REAL(KIND=8) PRMU0(KDLON)  ! COSINE OF ZENITHAL ANGLE
  REAL(KIND=8) PFRAC(KDLON)  ! fraction de la journee

  REAL(KIND=8) PTAVE(KDLON,KFLEV)  ! LAYER TEMPERATURE (K)
  REAL(KIND=8) PWV(KDLON,KFLEV)    ! SPECIFI! HUMIDITY (KG/KG)
  REAL(KIND=8) PQS(KDLON,KFLEV)    ! SATURATED WATER VAPOUR (KG/KG)
  REAL(KIND=8) POZON(KDLON,KFLEV)  ! OZONE CONCENTRATION (KG/KG)
  REAL(KIND=8) PAER(KDLON,KFLEV,5) ! AEROSOLS' OPTICAL THICKNESS

  REAL(KIND=8) PALBD(KDLON,2)  ! albedo du sol (lumiere diffuse)
  REAL(KIND=8) PALBP(KDLON,2)  ! albedo du sol (lumiere parallele)

  REAL(KIND=8) PCLDSW(KDLON,KFLEV)    ! CLOUD FRACTION
  REAL(KIND=8) PTAU(KDLON,2,KFLEV)    ! CLOUD OPTICAL THICKNESS (pre-industrial value)
  REAL(KIND=8) PCG(KDLON,2,KFLEV)     ! ASYMETRY FACTOR
  REAL(KIND=8) POMEGA(KDLON,2,KFLEV)  ! SINGLE SCATTERING ALBEDO

  REAL(KIND=8) PHEAT(KDLON,KFLEV) ! SHORTWAVE HEATING (K/DAY)
  REAL(KIND=8) PHEAT0(KDLON,KFLEV)! SHORTWAVE HEATING (K/DAY) clear-sky
  REAL(KIND=8) PALBPLA(KDLON)     ! PLANETARY ALBEDO
  REAL(KIND=8) PTOPSW(KDLON)      ! SHORTWAVE FLUX AT T.O.A.
  REAL(KIND=8) PSOLSW(KDLON)      ! SHORTWAVE FLUX AT SURFACE
  REAL(KIND=8) PTOPSW0(KDLON)     ! SHORTWAVE FLUX AT T.O.A. (CLEAR-SKY)
  REAL(KIND=8) PSOLSW0(KDLON)     ! SHORTWAVE FLUX AT SURFACE (CLEAR-SKY)
  !
  !* LOCAL VARIABLES:
  !
  real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

  REAL(KIND=8) ZOZ(KDLON,KFLEV)
  ! column-density of ozone in layer, in kilo-Dobsons

  REAL(KIND=8) ZAKI(KDLON,2)     
  REAL(KIND=8) ZCLD(KDLON,KFLEV)
  REAL(KIND=8) ZCLEAR(KDLON) 
  REAL(KIND=8) ZDSIG(KDLON,KFLEV)
  REAL(KIND=8) ZFACT(KDLON)
  REAL(KIND=8) ZFD(KDLON,KFLEV+1)
  REAL(KIND=8) ZFDOWN(KDLON,KFLEV+1)
  REAL(KIND=8) ZFU(KDLON,KFLEV+1)
  REAL(KIND=8) ZFUP(KDLON,KFLEV+1)
  REAL(KIND=8) ZRMU(KDLON)
  REAL(KIND=8) ZSEC(KDLON)
  REAL(KIND=8) ZUD(KDLON,5,KFLEV+1)
  REAL(KIND=8) ZCLDSW0(KDLON,KFLEV)

  REAL(KIND=8) ZFSUP(KDLON,KFLEV+1)
  REAL(KIND=8) ZFSDN(KDLON,KFLEV+1)
  REAL(KIND=8) ZFSUP0(KDLON,KFLEV+1)
  REAL(KIND=8) ZFSDN0(KDLON,KFLEV+1)

  INTEGER inu, jl, jk, i, k, kpl1

  INTEGER swpas  ! Every swpas steps, sw is calculated
  PARAMETER(swpas=1)

  INTEGER, SAVE :: itapsw = 0
  !$OMP THREADPRIVATE(itapsw)
  LOGICAL, SAVE :: appel1er = .TRUE.
  !$OMP THREADPRIVATE(appel1er)
  LOGICAL, SAVE :: initialized = .FALSE.
  !$OMP THREADPRIVATE(initialized)

  !jq-local flag introduced for aerosol forcings
  REAL(KIND=8), SAVE :: flag_aer
  !$OMP THREADPRIVATE(flag_aer)

  LOGICAL ok_ade, ok_aie    ! use aerosol forcings or not?
  INTEGER flag_aerosol_strat ! use stratospehric aerosols
  INTEGER flag_aerosol      ! global flag for aerosol 0 (no aerosol) or 1-5 (aerosols)
  REAL(KIND=8) tauaero(kdlon,kflev,naero_grp,2)  ! aerosol optical properties
  REAL(KIND=8) pizaero(kdlon,kflev,naero_grp,2)  ! (see aeropt.F)
  REAL(KIND=8) cgaero(kdlon,kflev,naero_grp,2)   ! -"-
  REAL(KIND=8) PTAUA(KDLON,2,KFLEV)    ! CLOUD OPTICAL THICKNESS (present-day value)
  REAL(KIND=8) POMEGAA(KDLON,2,KFLEV)  ! SINGLE SCATTERING ALBEDO
  REAL(KIND=8) PTOPSWADAERO(KDLON)     ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL DIR)
  REAL(KIND=8) PSOLSWADAERO(KDLON)     ! SHORTWAVE FLUX AT SURFACE(+AEROSOL DIR)
  REAL(KIND=8) PTOPSWAD0AERO(KDLON)    ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL DIR)
  REAL(KIND=8) PSOLSWAD0AERO(KDLON)    ! SHORTWAVE FLUX AT SURFACE(+AEROSOL DIR)
  REAL(KIND=8) PTOPSWAIAERO(KDLON)     ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL IND)
  REAL(KIND=8) PSOLSWAIAERO(KDLON)     ! SHORTWAVE FLUX AT SURFACE(+AEROSOL IND)
  REAL(KIND=8) PTOPSWAERO(KDLON,9)     ! SW TOA AS DRF nat & ant 
  REAL(KIND=8) PTOPSW0AERO(KDLON,9)    ! SW SRF AS DRF nat & ant 
  REAL(KIND=8) PSOLSWAERO(KDLON,9)     ! SW TOA CS DRF nat & ant
  REAL(KIND=8) PSOLSW0AERO(KDLON,9)    ! SW SRF CS DRF nat & ant
  REAL(KIND=8) PTOPSWCFAERO(KDLON,3)   !  SW TOA AS cloudRF nat & ant 
  REAL(KIND=8) PSOLSWCFAERO(KDLON,3)   !  SW SRF AS cloudRF nat & ant 

  !jq - Fluxes including aerosol effects
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSUPAD_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSUPAD_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSDNAD_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSDNAD_AERO)
  !jq - Fluxes including aerosol effects
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSUPAD0_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSUPAD0_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSDNAD0_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSDNAD0_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSUPAI_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSUPAI_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSDNAI_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSDNAI_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE ::  ZFSUP_AERO(:,:,:)
  !$OMP THREADPRIVATE(ZFSUP_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE ::  ZFSDN_AERO(:,:,:)
  !$OMP THREADPRIVATE(ZFSDN_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE ::  ZFSUP0_AERO(:,:,:)
  !$OMP THREADPRIVATE(ZFSUP0_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE ::  ZFSDN0_AERO(:,:,:)
  !$OMP THREADPRIVATE(ZFSDN0_AERO)

! Key to define the aerosol effect acting on climate
! OB: AEROSOLFEEDBACK_ACTIVE is now a LOGICAL
! TRUE: fluxes use natural and/or anthropogenic aerosols according to ok_ade and ok_aie, DEFAULT
! FALSE: fluxes use no aerosols (case 1)

  LOGICAL,SAVE :: AEROSOLFEEDBACK_ACTIVE = .TRUE.
!$OMP THREADPRIVATE(AEROSOLFEEDBACK_ACTIVE)  

      CHARACTER (LEN=20) :: modname='sw_aeroAR4'
      CHARACTER (LEN=80) :: abort_message

  IF(.NOT.initialized) THEN
     flag_aer=0.
     initialized=.TRUE.
     ALLOCATE(ZFSUPAD_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSDNAD_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSUPAD0_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSDNAD0_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSUPAI_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSDNAI_AERO(KDLON,KFLEV+1))
!-OB decrease size of these arrays to what is needed
!                | direct effect 
!ind effect      | no aerosol   natural  total
!natural (PTAU)  |   1            3       2     --ZFSUP/ZFSDN
!total (PTAUA)   |                5       4     --ZFSUP/ZFSDN
!no cloud        |   1            3       2     --ZFSUP0/ZFSDN0
! so we need which case when ? 
! ok_ade and ok_aie = 4-5, 4-2 and 2 
! ok_ade and not ok_aie = 2-3 and 2 
! not ok_ade and ok_aie = 5-3 and 5 
! not ok_ade and not ok_aie = 3 
! therefore the cases have the folliwng switches 
! 3 = not ok_ade or not ok_aie 
! 4 = ok_ade and ok_aie 
! 2 = ok_ade 
! 5 = ok_aie 
     ALLOCATE(ZFSUP_AERO (KDLON,KFLEV+1,5))
     ALLOCATE(ZFSDN_AERO (KDLON,KFLEV+1,5))
     ALLOCATE(ZFSUP0_AERO(KDLON,KFLEV+1,3))
     ALLOCATE(ZFSDN0_AERO(KDLON,KFLEV+1,3))
! end OB modif
     ZFSUPAD_AERO(:,:)=0.
     ZFSDNAD_AERO(:,:)=0.
     ZFSUPAD0_AERO(:,:)=0.
     ZFSDNAD0_AERO(:,:)=0.
     ZFSUPAI_AERO(:,:)=0.
     ZFSDNAI_AERO(:,:)=0.
     ZFSUP_AERO (:,:,:)=0.
     ZFSDN_AERO (:,:,:)=0.
     ZFSUP0_AERO(:,:,:)=0.
     ZFSDN0_AERO(:,:,:)=0.
  ENDIF

  IF (appel1er) THEN
     WRITE(lunout,*)'SW calling frequency : ', swpas
     WRITE(lunout,*) "   In general, it should be 1"
     appel1er = .FALSE.
  ENDIF
  !     ------------------------------------------------------------------
  IF (MOD(itapsw,swpas).EQ.0) THEN

     DO JK = 1 , KFLEV
        DO JL = 1, KDLON
           ZCLDSW0(JL,JK) = 0.0
           ZOZ(JL,JK) = POZON(JL,JK)*46.6968/RG &
                *PDP(JL,JK)*(101325.0/PPSOL(JL))
        ENDDO
     ENDDO

! clear sky with no aerosols at all is computed IF ACTIVEFEEDBACK_ACTIVE is false or for extended diag
     IF ( swaero_diag .or. .not. AEROSOLFEEDBACK_ACTIVE .OR. flag_aerosol .EQ. 0 ) THEN    

     ! clear-sky: zero aerosol effect
     flag_aer=0.0
     CALL SWU_LMDAR4(PSCT,ZCLDSW0,PPMB,PPSOL,&
          PRMU0,PFRAC,PTAVE,PWV,&
          ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
     INU = 1
     CALL SW1S_LMDAR4(INU,PAER, flag_aer, &
          tauaero(:,:,1,:), pizaero(:,:,1,:), cgaero(:,:,1,:),&
          PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0,&
          ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
          ZFD, ZFU)
     INU = 2
     CALL SW2S_LMDAR4(INU, PAER, flag_aer, &
          tauaero(:,:,1,:), pizaero(:,:,1,:), cgaero(:,:,1,:),&
          ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0,&
          ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
          PWV, PQS,&
          ZFDOWN, ZFUP)
     DO JK = 1 , KFLEV+1
        DO JL = 1, KDLON
           ZFSUP0_AERO(JL,JK,1) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
           ZFSDN0_AERO(JL,JK,1) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
        ENDDO
     ENDDO
     ENDIF ! swaero_diag .or. .not. AEROSOLFEEDBACK_ACTIVE

! cloudy sky with no aerosols at all is either computed IF no indirect effect is asked for, or for extended diag
     IF ( swaero_diag .or. .not. AEROSOLFEEDBACK_ACTIVE .OR. flag_aerosol .EQ. 0 ) THEN    
     ! cloudy-sky: zero aerosol effect
     flag_aer=0.0
     CALL SWU_LMDAR4(PSCT,PCLDSW,PPMB,PPSOL,&
          PRMU0,PFRAC,PTAVE,PWV,&
          ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
     INU = 1
     CALL SW1S_LMDAR4(INU, PAER, flag_aer, &
          tauaero(:,:,1,:), pizaero(:,:,1,:), cgaero(:,:,1,:),&
          PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
          ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
          ZFD, ZFU)
     INU = 2
     CALL SW2S_LMDAR4(INU, PAER, flag_aer, &
          tauaero(:,:,1,:), pizaero(:,:,1,:), cgaero(:,:,1,:),&
          ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
          ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
          PWV, PQS,&
          ZFDOWN, ZFUP)

     DO JK = 1 , KFLEV+1
        DO JL = 1, KDLON
           ZFSUP_AERO(JL,JK,1) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
           ZFSDN_AERO(JL,JK,1) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
        ENDDO
     ENDDO
     ENDIF ! swaero_diag .or. .not. AEROSOLFEEDBACK_ACTIVE

     IF (flag_aerosol.GT.0 .OR. flag_aerosol_strat.GT.0) THEN

     IF (ok_ade.and.swaero_diag .or. .not. ok_ade) THEN

        ! clear sky direct effect natural aerosol
        ! CAS AER (3)
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,ZCLDSW0,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,3,:), pizaero(:,:,3,:), cgaero(:,:,3,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,3,:), pizaero(:,:,3,:), cgaero(:,:,3,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)

        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP0_AERO(JL,JK,3) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
              ZFSDN0_AERO(JL,JK,3) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
           ENDDO
        ENDDO
     ENDIF !--end not swaero_diag or not ok_ade

     IF (ok_ade) THEN

        ! clear sky direct effect of total aerosol
        ! CAS AER (2)
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,ZCLDSW0,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)

        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP0_AERO(JL,JK,2) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL) 
              ZFSDN0_AERO(JL,JK,2) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL) 
           ENDDO
        ENDDO

        ! cloudy-sky with natural aerosols for indirect effect 
        ! but total aerosols for direct effect
        ! PTAU
        ! CAS AER (2)
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,PCLDSW,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)

        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP_AERO(JL,JK,2) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL) 
              ZFSDN_AERO(JL,JK,2) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL) 
           ENDDO
        ENDDO

     ENDIF !-end ok_ade

     IF ( .not. ok_ade .or. .not. ok_aie ) THEN

        ! cloudy-sky with natural aerosols for indirect effect 
        ! and natural aerosols for direct effect
        ! PTAU
        ! CAS AER (3)
        ! cloudy-sky direct effect natural aerosol
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,PCLDSW,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,3,:), pizaero(:,:,3,:), cgaero(:,:,3,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,3,:), pizaero(:,:,3,:), cgaero(:,:,3,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)

        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP_AERO(JL,JK,3) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
              ZFSDN_AERO(JL,JK,3) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
           ENDDO
        ENDDO

     ENDIF  !--true/false or false/true

     IF (ok_ade .and. ok_aie) THEN

        ! cloudy-sky with total aerosols for indirect effect 
        ! and total aerosols for direct effect
        ! PTAUA
        ! CAS AER (2)
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,PCLDSW,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)

        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP_AERO(JL,JK,4) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
              ZFSDN_AERO(JL,JK,4) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL) 
           ENDDO
        ENDDO
  
      ENDIF ! ok_ade .and. ok_aie

     IF (ok_aie) THEN
        ! cloudy-sky with total aerosols for indirect effect 
        ! and natural aerosols for direct effect
        ! PTAUA
        ! CAS AER (3)
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,PCLDSW,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,3,:), pizaero(:,:,3,:), cgaero(:,:,3,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,3,:), pizaero(:,:,3,:), cgaero(:,:,3,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)
  
        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP_AERO(JL,JK,5) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
              ZFSDN_AERO(JL,JK,5) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL) 
           ENDDO
        ENDDO

     ENDIF ! ok_aie      

     ENDIF !--if flag_aerosol GT 0 OR flag_aerosol_strat GT 0

     itapsw = 0
  ENDIF
  itapsw = itapsw + 1

  IF  ( AEROSOLFEEDBACK_ACTIVE .AND. (flag_aerosol.GT.0 .OR. flag_aerosol_strat.GT.0) ) THEN
  IF ( ok_ade .and. ok_aie  ) THEN
    ZFSUP(:,:) =    ZFSUP_AERO(:,:,4)
    ZFSDN(:,:) =    ZFSDN_AERO(:,:,4)
    ZFSUP0(:,:) =   ZFSUP0_AERO(:,:,2)
    ZFSDN0(:,:) =   ZFSDN0_AERO(:,:,2)
  ENDIF

  IF ( ok_ade .and. (.not. ok_aie) )  THEN
    ZFSUP(:,:) =    ZFSUP_AERO(:,:,2)
    ZFSDN(:,:) =    ZFSDN_AERO(:,:,2)
    ZFSUP0(:,:) =   ZFSUP0_AERO(:,:,2)
    ZFSDN0(:,:) =   ZFSDN0_AERO(:,:,2)
  ENDIF

  IF ( (.not. ok_ade) .and. ok_aie  )  THEN
    ZFSUP(:,:) =    ZFSUP_AERO(:,:,5)
    ZFSDN(:,:) =    ZFSDN_AERO(:,:,5)
    ZFSUP0(:,:) =   ZFSUP0_AERO(:,:,3)
    ZFSDN0(:,:) =   ZFSDN0_AERO(:,:,3)
  ENDIF

  IF ((.not. ok_ade) .and. (.not. ok_aie)) THEN
    ZFSUP(:,:) =    ZFSUP_AERO(:,:,3)
    ZFSDN(:,:) =    ZFSDN_AERO(:,:,3)
    ZFSUP0(:,:) =   ZFSUP0_AERO(:,:,3)
    ZFSDN0(:,:) =   ZFSDN0_AERO(:,:,3)
  ENDIF

! MS the following allows to compute the forcing diagostics without
! letting the aerosol forcing act on the meteorology
! SEE logic above
  ELSE 
    ZFSUP(:,:) =    ZFSUP_AERO(:,:,1)
    ZFSDN(:,:) =    ZFSDN_AERO(:,:,1)
    ZFSUP0(:,:) =   ZFSUP0_AERO(:,:,1)
    ZFSDN0(:,:) =   ZFSDN0_AERO(:,:,1)
  ENDIF

! Now computes heating rates
  DO k = 1, KFLEV
     kpl1 = k+1
     DO i = 1, KDLON
        PHEAT(i,k) = -(ZFSUP(i,kpl1)-ZFSUP(i,k))-(ZFSDN(i,k)-ZFSDN(i,kpl1))
        PHEAT(i,k) = PHEAT(i,k) * RDAY*RG/RCPD / PDP(i,k)
        PHEAT0(i,k) = -(ZFSUP0(i,kpl1)-ZFSUP0(i,k))-(ZFSDN0(i,k)-ZFSDN0(i,kpl1))
        PHEAT0(i,k) = PHEAT0(i,k) * RDAY*RG/RCPD / PDP(i,k)
     ENDDO
  ENDDO

  DO i = 1, KDLON
! effective SW surface albedo calculation
     PALBPLA(i) = ZFSUP(i,KFLEV+1)/(ZFSDN(i,KFLEV+1)+1.0e-20)
     
! clear sky net fluxes at TOA and SRF
     PSOLSW0(i) = ZFSDN0(i,1) - ZFSUP0(i,1)
     PTOPSW0(i) = ZFSDN0(i,KFLEV+1) - ZFSUP0(i,KFLEV+1)

! cloudy sky net fluxes at TOA and SRF
     PSOLSW(i) = ZFSDN(i,1) - ZFSUP(i,1)
     PTOPSW(i) = ZFSDN(i,KFLEV+1) - ZFSUP(i,KFLEV+1)

! net anthropogenic forcing direct and 1st indirect effect diagnostics
! requires a natural aerosol field read and used 
! Difference of net fluxes from double call to radiation

IF (ok_ade) THEN

! indices 1: natural; 2 anthropogenic 

! TOA/SRF all sky natural forcing
     PSOLSWAERO(i,1) = (ZFSDN_AERO(i,1,3) - ZFSUP_AERO(i,1,3))-(ZFSDN_AERO(i,1,1) - ZFSUP_AERO(i,1,1))
     PTOPSWAERO(i,1) = (ZFSDN_AERO(i,KFLEV+1,3) - ZFSUP_AERO(i,KFLEV+1,3))- (ZFSDN_AERO(i,KFLEV+1,1) - ZFSUP_AERO(i,KFLEV+1,1))

! TOA/SRF clear sky natural forcing
     PSOLSW0AERO(i,1) = (ZFSDN0_AERO(i,1,3) - ZFSUP0_AERO(i,1,3))-(ZFSDN0_AERO(i,1,1) - ZFSUP0_AERO(i,1,1))
     PTOPSW0AERO(i,1) = (ZFSDN0_AERO(i,KFLEV+1,3) - ZFSUP0_AERO(i,KFLEV+1,3))-(ZFSDN0_AERO(i,KFLEV+1,1) - ZFSUP0_AERO(i,KFLEV+1,1))

   IF (ok_aie) THEN 

! TOA/SRF all sky anthropogenic forcing
     PSOLSWAERO(i,2) = (ZFSDN_AERO(i,1,4) - ZFSUP_AERO(i,1,4))-(ZFSDN_AERO(i,1,5) - ZFSUP_AERO(i,1,5))
     PTOPSWAERO(i,2) = (ZFSDN_AERO(i,KFLEV+1,4) - ZFSUP_AERO(i,KFLEV+1,4))- (ZFSDN_AERO(i,KFLEV+1,5) - ZFSUP_AERO(i,KFLEV+1,5))

   ELSE 

! TOA/SRF all sky anthropogenic forcing
     PSOLSWAERO(i,2) = (ZFSDN_AERO(i,1,2) - ZFSUP_AERO(i,1,2))-(ZFSDN_AERO(i,1,3) - ZFSUP_AERO(i,1,3))
     PTOPSWAERO(i,2) = (ZFSDN_AERO(i,KFLEV+1,2) - ZFSUP_AERO(i,KFLEV+1,2))- (ZFSDN_AERO(i,KFLEV+1,3) - ZFSUP_AERO(i,KFLEV+1,3))

   ENDIF

! TOA/SRF clear sky anthropogenic forcing
     PSOLSW0AERO(i,2) = (ZFSDN0_AERO(i,1,2) - ZFSUP0_AERO(i,1,2))-(ZFSDN0_AERO(i,1,3) - ZFSUP0_AERO(i,1,3))
     PTOPSW0AERO(i,2) = (ZFSDN0_AERO(i,KFLEV+1,2) - ZFSUP0_AERO(i,KFLEV+1,2))-(ZFSDN0_AERO(i,KFLEV+1,3) - ZFSUP0_AERO(i,KFLEV+1,3))

! direct anthropogenic forcing , as in old LMDzT, however differences of net fluxes
     PSOLSWADAERO(i) = PSOLSWAERO(i,2)
     PTOPSWADAERO(i) = PTOPSWAERO(i,2)
     PSOLSWAD0AERO(i) = PSOLSW0AERO(i,2)
     PTOPSWAD0AERO(i) = PTOPSW0AERO(i,2)

! OB: these diagnostics may not always work but who need them
! Cloud forcing indices 1: natural; 2 anthropogenic; 3: zero aerosol direct effect
! Instantaneously computed cloudy sky direct aerosol effect, cloud forcing due to aerosols above clouds
! natural
     PSOLSWCFAERO(i,1) = PSOLSWAERO(i,1) - PSOLSW0AERO(i,1)
     PTOPSWCFAERO(i,1) = PTOPSWAERO(i,1) - PTOPSW0AERO(i,1)

! Instantaneously computed cloudy SKY DIRECT aerosol effect, cloud forcing due to aerosols above clouds
! anthropogenic
     PSOLSWCFAERO(i,2) = PSOLSWAERO(i,2) - PSOLSW0AERO(i,2)
     PTOPSWCFAERO(i,2) = PTOPSWAERO(i,2) - PTOPSW0AERO(i,2)

! Cloudforcing without aerosol
! zero
     PSOLSWCFAERO(i,3) = (ZFSDN_AERO(i,1,1) - ZFSUP_AERO(i,1,1))-(ZFSDN0_AERO(i,1,1) - ZFSUP0_AERO(i,1,1))
     PTOPSWCFAERO(i,3) = (ZFSDN_AERO(i,KFLEV+1,1) - ZFSUP_AERO(i,KFLEV+1,1))- (ZFSDN0_AERO(i,KFLEV+1,1) - ZFSUP0_AERO(i,KFLEV+1,1))

ENDIF

IF (ok_aie) THEN
   IF (ok_ade) THEN 
     PSOLSWAIAERO(i) = (ZFSDN_AERO(i,1,4) - ZFSUP_AERO(i,1,4))-(ZFSDN_AERO(i,1,2) - ZFSUP_AERO(i,1,2))
     PTOPSWAIAERO(i) = (ZFSDN_AERO(i,KFLEV+1,4) - ZFSUP_AERO(i,KFLEV+1,4))-(ZFSDN_AERO(i,KFLEV+1,2) - ZFSUP_AERO(i,KFLEV+1,2))
   ELSE 
     PSOLSWAIAERO(i) = (ZFSDN_AERO(i,1,5) - ZFSUP_AERO(i,1,5))-(ZFSDN_AERO(i,1,3) - ZFSUP_AERO(i,1,3))
     PTOPSWAIAERO(i) = (ZFSDN_AERO(i,KFLEV+1,5) - ZFSUP_AERO(i,KFLEV+1,5))-(ZFSDN_AERO(i,KFLEV+1,3) - ZFSUP_AERO(i,KFLEV+1,3))
   ENDIF
ENDIF

ENDDO

END SUBROUTINE SW_AEROAR4
