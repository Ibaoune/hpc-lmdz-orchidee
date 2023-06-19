!$Id $
!
MODULE tracinca_mod
!
! This module prepares and calls the INCA main subroutines. 
!
  IMPLICIT NONE  

  CHARACTER(len=4),SAVE :: config_inca
!$OMP THREADPRIVATE(config_inca)
                     ! config_inca='none' => without INCA
                     ! config_inca='chem' => INCA with chemistry
                     ! config_inca='aero' => INCA with aerosols
CONTAINS

  SUBROUTINE tracinca_init(aerosol,lessivage)
    ! This subroutine initialize some control varaibles. 

    USE infotrac_phy, ONLY: nbtr
    USE ioipsl_getin_p_mod, ONLY: getin_p
    IMPLICIT NONE
    
    ! Output variables
    LOGICAL,DIMENSION(nbtr), INTENT(OUT) :: aerosol
    LOGICAL,INTENT(OUT) :: lessivage
    
    
    ! Initialization
    lessivage  =.FALSE.
    aerosol(:) = .FALSE.

  END SUBROUTINE tracinca_init

  SUBROUTINE tracinca(                                &
       nstep,    julien,   gmtime,         lafin,     &
       pdtphys,  t_seri,   paprs,          pplay,     &
       pmfu,     upwd,     ftsol,  pctsrf, pphis,     &
       pphi,     albsol,   sh,             ch, rh,    &
       cldfra,   rneb,     diafra,         cldliq,    &
       itop_con, ibas_con, pmflxr,         pmflxs,    &
       prfl,     psfl,     aerosol_couple, flxmass_w, &
       tau_aero, piz_aero, cg_aero,        ccm,       &
       rfname,                                        &
       tr_seri,  source)      

!========================================================
!    -- CHIMIE INCA --
!========================================================

    USE dimphy
    USE infotrac_phy, ONLY: nbtr
    USE vampir
    USE indice_sol_mod
    USE geometry_mod, ONLY: cell_area
    USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
    USE aero_mod, ONLY : naero_grp
    IMPLICIT NONE
    
!==========================================================================
!                   -- DESCRIPTION DES ARGUMENTS --
!==========================================================================


! EN ENTREE ...
!
!Configuration grille,temps:
    INTEGER,INTENT(IN) :: nstep      ! Appel physique
    INTEGER,INTENT(IN) :: julien     ! Jour julien
    REAL,INTENT(IN)    :: gmtime
    REAL,INTENT(IN)    :: pdtphys    ! Pas d'integration pour la physique (seconde)
    LOGICAL,INTENT(IN) :: lafin      ! le flag de la fin de la physique
    

!Physique: 
!--------
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: sh      ! humidite specifique
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: ch      ! eau liquide
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: rh      ! humidite relative
    REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs   ! pression pour chaque inter-couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pphi    ! geopotentiel
    REAL,DIMENSION(klon),INTENT(IN)        :: pphis
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: cldliq  ! eau liquide nuageuse
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: cldfra  ! fraction nuageuse (tous les nuages)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: diafra  ! fraction nuageuse (convection ou stratus artificiels)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: rneb    ! fraction nuageuse (grande echelle)
    INTEGER,DIMENSION(klon),INTENT(IN)     :: itop_con
    INTEGER,DIMENSION(klon),INTENT(IN)     :: ibas_con
    REAL,DIMENSION(klon),INTENT(IN)        :: albsol  ! albedo surface
!
!Convection:
!----------
    REAL,DIMENSION(klon,klev),INTENT(IN) :: pmfu  ! flux de masse dans le panache montant - Tiedtke
    REAL,DIMENSION(klon,klev),INTENT(IN) :: upwd  ! flux de masse dans le panache montant - Emanuel

!...Tiedke     
    REAL,DIMENSION(klon,klev+1),INTENT(IN)   :: pmflxr, pmflxs ! Flux precipitant de pluie, neige aux interfaces [convection]
    REAL,DIMENSION(klon,klev+1),INTENT(IN)   :: prfl, psfl ! Flux precipitant de pluie, neige aux interfaces [large-scale]

    LOGICAL,INTENT(IN)                       :: aerosol_couple
    REAL,DIMENSION(klon,klev),INTENT(IN)     :: flxmass_w
    REAL,DIMENSION(klon,klev,naero_grp,2),INTENT(IN) :: tau_aero
    REAL,DIMENSION(klon,klev,naero_grp,2),INTENT(IN) :: piz_aero
    REAL,DIMENSION(klon,klev,naero_grp,2),INTENT(IN) :: cg_aero
    CHARACTER(len=4),DIMENSION(naero_grp),INTENT(IN) :: rfname 
    REAL,DIMENSION(klon,klev,2),INTENT(IN)   :: ccm 

! Arguments necessaires pour les sources et puits de traceur:
    REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: ftsol  ! Temperature du sol (surf)(Kelvin)
    REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: pctsrf ! Pourcentage de sol f(nature du sol)


  ! InOutput argument
    REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT) :: tr_seri ! Concentration Traceur [U/KgA]  

  ! Output arguments
    REAL,DIMENSION(klon,nbtr), INTENT(OUT)        :: source  ! a voir lorsque le flux de surface est prescrit 

!=======================================================================================
!                        -- VARIABLES LOCALES TRACEURS --
!=======================================================================================

    INTEGER :: k
    REAL,DIMENSION(klon,klev) :: pdel
    REAL,DIMENSION(klon,klev) :: zpmfu  ! flux de masse dans le panache montant
    REAL    :: calday
    INTEGER :: ncsec

    CALL VTe(VTphysiq)
    CALL VTb(VTinca)
    
    calday = REAL(julien) + gmtime
    ncsec  = NINT (86400.*gmtime)
     
    DO k = 1, klev
       pdel(:,k) = paprs(:,k) - paprs (:,k+1)
    END DO
  
#ifdef INCA
    IF (config_inca == 'aero' .OR. config_inca == 'chem') THEN 
       zpmfu(:,:)=pmfu(:,:)       
    ELSE IF (config_inca == 'aeNP') THEN
       zpmfu(:,:)=upwd(:,:)
    ENDIF

    CALL aerosolmain(                    &
         aerosol_couple,tr_seri,pdtphys, &
         pplay,pdel,prfl,pmflxr,psfl,    &
         pmflxs,zpmfu,itop_con,ibas_con,  &
         pphi,cell_area,nstep,rneb,t_seri, &      
         rh,tau_aero,piz_aero,cg_aero,   &
         rfname,ccm,lafin)
#endif


#ifdef INCA
    CALL chemmain (tr_seri, &   !mmr
         nstep,      & !nstep
         calday,     & !calday
         julien,     & !ncdate
         ncsec,      & !ncsec
         1,          & !lat
         pdtphys,    & !delt
         paprs(1,1), & !ps
         pplay,      & !pmid
         pdel,       & !pdel
         cell_area,  &
         pctsrf(1,1),& !oro
         ftsol,      & !tsurf
         albsol,     & !albs
         pphi,       & !zma
         pphis,      & !phis
         cldfra,     & !cldfr
         rneb,       & !cldfr_st
         diafra,     & !cldfr_cv
         itop_con,   & !cldtop
         ibas_con,   & !cldbot
         cldliq,     & !cwat
         prfl,       & !flxrst
         pmflxr,     & !flxrcv
         psfl,       & !flxsst
         pmflxs,     & !flxscv
         zpmfu,      & !flxupd   !--now depends on whether AP or NP
         flxmass_w,  & !flxmass_w
         t_seri,     & !tfld
         sh,         & !sh
         ch,         & !ql
         rh,         & !rh
         nbp_lon,    & !nx
         nbp_lat,    & !ny
         source )
#endif
    
    CALL VTe(VTinca)
    CALL VTb(VTphysiq)
    
    
  END SUBROUTINE tracinca


END MODULE tracinca_mod
