! $Id: readaerosol_optic_rrtm.F90 3630 2020-02-10 10:04:40Z fairhead $
!
SUBROUTINE readaerosol_optic_rrtm(debut, aerosol_couple, ok_alw, ok_volcan, &
     flag_aerosol, flag_bc_internal_mixture, itap, rjourvrai, &
     pdtphys, pplay, paprs, t_seri, rhcl, presnivs, &
     tr_seri, mass_solu_aero, mass_solu_aero_pi, &
     tau_aero, piz_aero, cg_aero, &
     tausum_aero, drytausum_aero, tau3d_aero )

  ! This routine will :
  ! 1) recevie the aerosols(already read and interpolated) corresponding to flag_aerosol
  ! 2) calculate the optical properties for the aerosols
  !

  USE dimphy
  USE aero_mod
  USE phys_local_var_mod, only: sconcso4,sconcno3,sconcoa,sconcbc,sconcss,sconcdust, &
       concso4,concno3,concoa,concbc,concss,concdust,loadso4,loadoa,loadbc,loadss,loaddust, &
       loadno3,load_tmp1,load_tmp2,load_tmp3,load_tmp4,load_tmp5,load_tmp6,load_tmp7, & 
       load_tmp8,load_tmp9,load_tmp10

  USE infotrac_phy
  USE YOMCST

  IMPLICIT NONE

  include "clesphys.h"

  ! Input arguments
  !****************************************************************************************
  LOGICAL, INTENT(IN)                      :: debut
  LOGICAL, INTENT(IN)                      :: aerosol_couple
  LOGICAL, INTENT(IN)                      :: ok_alw
  LOGICAL, INTENT(IN)                      :: ok_volcan
  INTEGER, INTENT(IN)                      :: flag_aerosol
  LOGICAL, INTENT(IN)                      :: flag_bc_internal_mixture
  INTEGER, INTENT(IN)                      :: itap
  REAL, INTENT(IN)                         :: rjourvrai
  REAL, INTENT(IN)                         :: pdtphys
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: pplay
  REAL, DIMENSION(klon,klev+1), INTENT(IN) :: paprs
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: t_seri
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: rhcl   ! humidite relative ciel clair
  REAL, DIMENSION(klev), INTENT(IN)        :: presnivs
  REAL, DIMENSION(klon,klev,nbtr), INTENT(IN) :: tr_seri ! concentration tracer

  ! Output arguments
  !****************************************************************************************
  REAL, DIMENSION(klon,klev), INTENT(OUT)     :: mass_solu_aero    ! Total mass for all soluble aerosols
  REAL, DIMENSION(klon,klev), INTENT(OUT)     :: mass_solu_aero_pi !     -"-     preindustrial values
  REAL, DIMENSION(klon,klev,2,NSW), INTENT(OUT) :: tau_aero    ! Aerosol optical thickness
  REAL, DIMENSION(klon,klev,2,NSW), INTENT(OUT) :: piz_aero    ! Single scattering albedo aerosol
  REAL, DIMENSION(klon,klev,2,NSW), INTENT(OUT) :: cg_aero     ! asymmetry parameter aerosol
  REAL, DIMENSION(klon,nwave,naero_tot), INTENT(OUT)       :: tausum_aero
  REAL, DIMENSION(klon,naero_tot), INTENT(OUT)             :: drytausum_aero
  REAL, DIMENSION(klon,klev,nwave,naero_tot), INTENT(OUT)  :: tau3d_aero

  ! Local variables
  !****************************************************************************************
  REAL, DIMENSION(klon)        :: aerindex      ! POLDER aerosol index 
  REAL, DIMENSION(klon,klev)   :: sulfacc       ! SO4 accumulation concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: sulfcoarse    ! SO4 coarse concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: bcsol         ! BC soluble concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: bcins         ! BC insoluble concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: pomsol        ! POM soluble concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: pomins        ! POM insoluble concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: cidust        ! DUST aerosol concentration  [ug/m3]
  REAL, DIMENSION(klon,klev)   :: sscoarse      ! SS Coarse concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: sssupco       ! SS Super Coarse concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: ssacu         ! SS Acumulation concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: nitracc       ! nitrate accumulation concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: nitrcoarse    ! nitrate coarse concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: nitrinscoarse ! nitrate insoluble coarse concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: sulfacc_pi
  REAL, DIMENSION(klon,klev)   :: sulfcoarse_pi
  REAL, DIMENSION(klon,klev)   :: bcsol_pi
  REAL, DIMENSION(klon,klev)   :: bcins_pi
  REAL, DIMENSION(klon,klev)   :: pomsol_pi
  REAL, DIMENSION(klon,klev)   :: pomins_pi
  REAL, DIMENSION(klon,klev)   :: cidust_pi
  REAL, DIMENSION(klon,klev)   :: sscoarse_pi
  REAL, DIMENSION(klon,klev)   :: sssupco_pi
  REAL, DIMENSION(klon,klev)   :: ssacu_pi
  REAL, DIMENSION(klon,klev)   :: nitracc_pi
  REAL, DIMENSION(klon,klev)   :: nitrcoarse_pi
  REAL, DIMENSION(klon,klev)   :: nitrinscoarse_pi
  REAL, DIMENSION(klon,klev)   :: pdel, zrho
  REAL, DIMENSION(klon,klev,naero_tot) :: m_allaer
  REAL, DIMENSION(klon,klev,naero_tot) :: m_allaer_pi !RAF  

  integer :: id_ASBCM, id_ASPOMM, id_ASSO4M, id_ASMSAM, id_CSSO4M, id_CSMSAM, id_SSSSM
  integer :: id_CSSSM, id_ASSSM, id_CIDUSTM, id_AIBCM, id_AIPOMM, id_ASNO3M, id_CSNO3M, id_CINO3M
  INTEGER :: k, i

  !--air density
  zrho(:,:)=pplay(:,:)/t_seri(:,:)/RD                     !--kg/m3

  !****************************************************************************************
  ! 1) Get aerosol mass
  !    
  !****************************************************************************************
  !
  !
  IF (aerosol_couple) THEN   !--we get aerosols from tr_seri array from INCA
     !
     !--copy fields from INCA tr_seri 
     !--convert to ug m-3 unit for consistency with offline fields
     !
     DO i=1,nbtr
        SELECT CASE(trim(solsym(i)))
           CASE ("ASBCM")
              id_ASBCM = i
           CASE ("ASPOMM") 
              id_ASPOMM = i 
           CASE ("ASSO4M")
              id_ASSO4M = i 
           CASE ("ASMSAM")
              id_ASMSAM = i 
           CASE ("CSSO4M")
              id_CSSO4M = i 
           CASE ("CSMSAM")
              id_CSMSAM = i 
           CASE ("SSSSM")
              id_SSSSM = i 
           CASE ("CSSSM")
              id_CSSSM = i 
           CASE ("ASSSM")
              id_ASSSM = i 
           CASE ("CIDUSTM")
              id_CIDUSTM = i 
           CASE ("AIBCM")
              id_AIBCM = i 
           CASE ("AIPOMM")
              id_AIPOMM = i 
           CASE ("ASNO3M")
              id_ASNO3M = i 
           CASE ("CSNO3M")
              id_CSNO3M = i 
           CASE ("CINO3M")
              id_CINO3M = i 
           END SELECT
     ENDDO

     bcsol(:,:)        =   tr_seri(:,:,id_ASBCM)                         *zrho(:,:)*1.e9  ! ASBCM
     pomsol(:,:)       =   tr_seri(:,:,id_ASPOMM)                        *zrho(:,:)*1.e9  ! ASPOMM
     sulfacc(:,:)      =  (tr_seri(:,:,id_ASSO4M)+tr_seri(:,:,id_ASMSAM))*zrho(:,:)*1.e9  ! ASSO4M (=SO4) + ASMSAM (=MSA)
     sulfcoarse(:,:)   =  (tr_seri(:,:,id_CSSO4M)+tr_seri(:,:,id_CSMSAM))*zrho(:,:)*1.e9  ! CSSO4M (=SO4) + CSMSAM (=MSA)
     sssupco(:,:)      =   tr_seri(:,:,id_SSSSM)                         *zrho(:,:)*1.e9  ! SSSSM
     sscoarse(:,:)     =   tr_seri(:,:,id_CSSSM)                         *zrho(:,:)*1.e9  ! CSSSM
     ssacu(:,:)        =   tr_seri(:,:,id_ASSSM)                         *zrho(:,:)*1.e9  ! ASSSM
     cidust(:,:)       =   tr_seri(:,:,id_CIDUSTM)                       *zrho(:,:)*1.e9  ! CIDUSTM
     bcins(:,:)        =   tr_seri(:,:,id_AIBCM)                         *zrho(:,:)*1.e9  ! AIBCM
     pomins(:,:)       =   tr_seri(:,:,id_AIPOMM)                        *zrho(:,:)*1.e9  ! AIPOMM
     nitracc(:,:)      =   tr_seri(:,:,id_ASNO3M)                        *zrho(:,:)*1.e9  ! ASNO3M
     nitrcoarse(:,:)   =   tr_seri(:,:,id_CSNO3M)                        *zrho(:,:)*1.e9  ! CSNO3M
     nitrinscoarse(:,:)=   tr_seri(:,:,id_CINO3M)                        *zrho(:,:)*1.e9  ! CINO3M
     !
     bcsol_pi(:,:)        =   0.0 ! ASBCM pre-ind
     pomsol_pi(:,:)       =   0.0 ! ASPOMM pre-ind
     sulfacc_pi(:,:)      =   0.0 ! ASSO4M (=SO4) + ASMSAM (=MSA) pre-ind
     sulfcoarse_pi(:,:)   =   0.0 ! CSSO4M (=SO4) + CSMSAM (=MSA) pre-ind
     sssupco_pi(:,:)      =   0.0 ! SSSSM pre-ind
     sscoarse_pi(:,:)     =   0.0 ! CSSSM pre-ind
     ssacu_pi(:,:)        =   0.0 ! ASSSM pre-ind
     cidust_pi(:,:)       =   0.0 ! CIDUSTM pre-ind
     bcins_pi(:,:)        =   0.0 ! AIBCM pre-ind
     pomins_pi(:,:)       =   0.0 ! AIPOMM pre-ind
     nitracc_pi(:,:)      =   0.0 ! ASNO3M pre-ind
     nitrcoarse_pi(:,:)   =   0.0 ! CSNO3M pre-ind
     nitrinscoarse_pi(:,:)=   0.0 ! CINO3M
     !
  ELSE !--not aerosol_couple
     !
     ! Read and interpolate sulfate
     IF ( flag_aerosol .EQ. 1 .OR. flag_aerosol .EQ. 6 .OR. flag_aerosol .EQ. 7 ) THEN 

        CALL readaerosol_interp(id_ASSO4M_phy, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, sulfacc, sulfacc_pi,loadso4)
     ELSE
        sulfacc(:,:) = 0. ; sulfacc_pi(:,:) = 0.
        loadso4=0.
     ENDIF

     ! Read and interpolate bcsol and bcins
     IF ( flag_aerosol .EQ. 2 .OR. flag_aerosol .EQ. 6 .OR. flag_aerosol .EQ. 7 ) THEN 

        ! Get bc aerosol distribution 
        CALL readaerosol_interp(id_ASBCM_phy, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, bcsol, bcsol_pi, load_tmp1 )
        CALL readaerosol_interp(id_AIBCM_phy, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, bcins, bcins_pi, load_tmp2 )
        loadbc(:)=load_tmp1(:)+load_tmp2(:)
     ELSE
        bcsol(:,:) = 0. ; bcsol_pi(:,:) = 0.
        bcins(:,:) = 0. ; bcins_pi(:,:) = 0.
        loadbc=0.
     ENDIF

     ! Read and interpolate pomsol and pomins
     IF ( flag_aerosol .EQ. 3 .OR. flag_aerosol .EQ. 6 .OR. flag_aerosol .EQ. 7 ) THEN

        CALL readaerosol_interp(id_ASPOMM_phy, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, pomsol, pomsol_pi, load_tmp3)
        CALL readaerosol_interp(id_AIPOMM_phy, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, pomins, pomins_pi, load_tmp4)
        loadoa(:)=load_tmp3(:)+load_tmp4(:)
     ELSE
        pomsol(:,:) = 0. ; pomsol_pi(:,:) = 0.
        pomins(:,:) = 0. ; pomins_pi(:,:) = 0.
        loadoa=0.
     ENDIF

     ! Read and interpolate csssm, ssssm, assssm
     IF (flag_aerosol .EQ. 4 .OR. flag_aerosol .EQ. 6 .OR. flag_aerosol .EQ. 7 ) THEN 

        CALL readaerosol_interp(id_SSSSM_phy ,itap, pdtphys,rjourvrai, &
        debut, pplay, paprs, t_seri, sssupco, sssupco_pi, load_tmp5) 
        CALL readaerosol_interp(id_CSSSM_phy ,itap, pdtphys,rjourvrai, &
        debut, pplay, paprs, t_seri, sscoarse,sscoarse_pi, load_tmp6) 
        CALL readaerosol_interp(id_ASSSM_phy ,itap, pdtphys,rjourvrai, &
        debut, pplay, paprs, t_seri, ssacu, ssacu_pi, load_tmp7) 
        loadss(:)=load_tmp5(:)+load_tmp6(:)+load_tmp7(:)
     ELSE
        sscoarse(:,:) = 0. ; sscoarse_pi(:,:) = 0. 
        ssacu(:,:)    = 0. ; ssacu_pi(:,:) = 0. 
        sssupco(:,:)  = 0. ; sssupco_pi = 0. 
        loadss=0.
     ENDIF

     ! Read and interpolate cidustm
     IF (flag_aerosol .EQ. 5 .OR. flag_aerosol .EQ. 6 .OR. flag_aerosol .EQ. 7 ) THEN 

        CALL readaerosol_interp(id_CIDUSTM_phy, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, cidust, cidust_pi, loaddust) 

     ELSE
        cidust(:,:) = 0. ; cidust_pi(:,:) = 0. 
        loaddust=0.
     ENDIF
     !
     ! Read and interpolate asno3m, csno3m, cino3m
     IF (flag_aerosol .EQ. 6 .OR. flag_aerosol .EQ. 7 ) THEN 

        CALL readaerosol_interp(id_ASNO3M_phy, itap, pdtphys, rjourvrai, & 
        debut, pplay, paprs, t_seri, nitracc, nitracc_pi, load_tmp8) 
        CALL readaerosol_interp(id_CSNO3M_phy, itap, pdtphys, rjourvrai, & 
        debut, pplay, paprs, t_seri, nitrcoarse, nitrcoarse_pi, load_tmp9) 
        CALL readaerosol_interp(id_CINO3M_phy, itap, pdtphys, rjourvrai, & 
        debut, pplay, paprs, t_seri, nitrinscoarse, nitrinscoarse_pi, load_tmp10) 
        loadno3(:)=load_tmp8(:)+load_tmp9(:)+load_tmp10(:)

     ELSE
        nitracc(:,:)         =   0.0 ; nitracc_pi(:,:)      =   0.0 
        nitrcoarse(:,:)      =   0.0 ; nitrcoarse_pi(:,:)   =   0.0
        nitrinscoarse(:,:)   =   0.0 ; nitrinscoarse_pi(:,:)=   0.0
        loadno3(:)=0.0
     ENDIF
     !
     ! CSSO4M is set to 0 as not reliable
     sulfcoarse(:,:)      =   0.0 ! CSSO4M (=SO4) + CSMSAM (=MSA) 
     sulfcoarse_pi(:,:)   =   0.0 ! CSSO4M (=SO4) + CSMSAM (=MSA) pre-ind

  ENDIF !--not aerosol_couple

  !
  ! Store all aerosols in one variable
  !
  m_allaer(:,:,id_ASBCM_phy)  = bcsol(:,:)        ! ASBCM
  m_allaer(:,:,id_ASPOMM_phy) = pomsol(:,:)       ! ASPOMM
  m_allaer(:,:,id_ASSO4M_phy) = sulfacc(:,:)      ! ASSO4M (= SO4) 
  m_allaer(:,:,id_CSSO4M_phy) = sulfcoarse(:,:)   ! CSSO4M 
  m_allaer(:,:,id_SSSSM_phy)  = sssupco(:,:)      ! SSSSM
  m_allaer(:,:,id_CSSSM_phy)  = sscoarse(:,:)     ! CSSSM
  m_allaer(:,:,id_ASSSM_phy)  = ssacu(:,:)        ! ASSSM
  m_allaer(:,:,id_CIDUSTM_phy)= cidust(:,:)       ! CIDUSTM
  m_allaer(:,:,id_AIBCM_phy)  = bcins(:,:)        ! AIBCM
  m_allaer(:,:,id_ASNO3M_phy) = nitracc(:,:)      ! ASNO3M
  m_allaer(:,:,id_CSNO3M_phy) = nitrcoarse(:,:)   ! CSNO3M
  m_allaer(:,:,id_CINO3M_phy) = nitrinscoarse(:,:)! CINO3M
  m_allaer(:,:,id_AIPOMM_phy) = pomins(:,:)       ! AIPOMM
  m_allaer(:,:,id_STRAT_phy)  = 0.0

  !RAF
  m_allaer_pi(:,:,id_ASBCM_phy)  = bcsol_pi(:,:)        ! ASBCM pre-ind
  m_allaer_pi(:,:,id_ASPOMM_phy) = pomsol_pi(:,:)       ! ASPOMM pre-ind
  m_allaer_pi(:,:,id_ASSO4M_phy) = sulfacc_pi(:,:)      ! ASSO4M (= SO4) pre-ind
  m_allaer_pi(:,:,id_CSSO4M_phy) = sulfcoarse_pi(:,:)   ! CSSO4M pre-ind
  m_allaer_pi(:,:,id_SSSSM_phy)  = sssupco_pi(:,:)      ! SSSSM pre-ind
  m_allaer_pi(:,:,id_CSSSM_phy)  = sscoarse_pi(:,:)     ! CSSSM pre-ind
  m_allaer_pi(:,:,id_ASSSM_phy)  = ssacu_pi(:,:)        ! ASSSM pre-ind
  m_allaer_pi(:,:,id_CIDUSTM_phy)= cidust_pi(:,:)       ! CIDUSTM pre-ind
  m_allaer_pi(:,:,id_AIBCM_phy)  = bcins_pi(:,:)        ! AIBCM pre-ind
  m_allaer_pi(:,:,id_ASNO3M_phy) = nitracc_pi(:,:)      ! ASNO3M pre-ind
  m_allaer_pi(:,:,id_CSNO3M_phy) = nitrcoarse_pi(:,:)   ! CSNO3M pre-ind
  m_allaer_pi(:,:,id_CINO3M_phy) = nitrinscoarse_pi(:,:)! CINO3M pre-ind
  m_allaer_pi(:,:,id_AIPOMM_phy) = pomins_pi(:,:)       ! AIPOMM pre-ind
  m_allaer_pi(:,:,id_STRAT_phy)  = 0.0

  !
  ! Calculate the total mass of all soluble aersosols
  ! to be revisited for AR6
  mass_solu_aero(:,:)    = sulfacc(:,:)    + bcsol(:,:)    + pomsol(:,:)    + nitracc(:,:)    + ssacu(:,:)
  mass_solu_aero_pi(:,:) = sulfacc_pi(:,:) + bcsol_pi(:,:) + pomsol_pi(:,:) + nitracc_pi(:,:) + ssacu_pi(:,:)

  !****************************************************************************************
  ! 2) Calculate optical properties for the aerosols
  !
  !****************************************************************************************
  DO k = 1, klev
     DO i = 1, klon
        pdel(i,k) = paprs(i,k) - paprs (i,k+1)
     ENDDO
  ENDDO

!--new aerosol properties
  ! aeropt_6bands for rrtm
  CALL aeropt_6bands_rrtm(          &
       pdel, m_allaer, rhcl,        & 
       tau_aero, piz_aero, cg_aero, &
       m_allaer_pi, flag_aerosol,   &
       flag_bc_internal_mixture, zrho, ok_volcan ) 

  ! aeropt_5wv only for validation and diagnostics
  CALL aeropt_5wv_rrtm(              &
       pdel, m_allaer,               &
       rhcl, aerindex,               & 
       flag_aerosol,                 & 
       flag_bc_internal_mixture,     & 
       pplay, t_seri,                &
       tausum_aero, drytausum_aero, tau3d_aero )

  !--call LW optical properties for tropospheric aerosols 
  CALL aeropt_lw_rrtm(ok_alw, pdel, zrho, flag_aerosol, m_allaer, m_allaer_pi)

  ! Diagnostics calculation for CMIP5 protocol
  sconcso4(:)  =m_allaer(:,1,id_ASSO4M_phy)*1.e-9
  sconcno3(:)  =(m_allaer(:,1,id_ASNO3M_phy)+m_allaer(:,1,id_CSNO3M_phy)+m_allaer(:,1,id_CINO3M_phy))*1.e-9
  sconcoa(:)   =(m_allaer(:,1,id_ASPOMM_phy)+m_allaer(:,1,id_AIPOMM_phy))*1.e-9
  sconcbc(:)   =(m_allaer(:,1,id_ASBCM_phy)+m_allaer(:,1,id_AIBCM_phy))*1.e-9
  sconcss(:)   =(m_allaer(:,1,id_ASSSM_phy)+m_allaer(:,1,id_CSSSM_phy)+m_allaer(:,1,id_SSSSM_phy))*1.e-9
  sconcdust(:) =m_allaer(:,1,id_CIDUSTM_phy)*1.e-9
  concso4(:,:) =m_allaer(:,:,id_ASSO4M_phy)*1.e-9
  concno3(:,:) =(m_allaer(:,:,id_ASNO3M_phy)+m_allaer(:,:,id_CSNO3M_phy)+m_allaer(:,:,id_CINO3M_phy))*1.e-9
  concoa(:,:)  =(m_allaer(:,:,id_ASPOMM_phy)+m_allaer(:,:,id_AIPOMM_phy))*1.e-9
  concbc(:,:)  =(m_allaer(:,:,id_ASBCM_phy)+m_allaer(:,:,id_AIBCM_phy))*1.e-9
  concss(:,:)  =(m_allaer(:,:,id_ASSSM_phy)+m_allaer(:,:,id_CSSSM_phy)+m_allaer(:,:,id_SSSSM_phy))*1.e-9
  concdust(:,:)=m_allaer(:,:,id_CIDUSTM_phy)*1.e-9

END SUBROUTINE readaerosol_optic_rrtm
