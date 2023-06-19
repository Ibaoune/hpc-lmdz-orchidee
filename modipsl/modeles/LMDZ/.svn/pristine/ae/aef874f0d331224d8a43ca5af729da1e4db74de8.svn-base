MODULE YOMARPHY

IMPLICIT NONE

SAVE

!*
!    -----------------------------------------------------------------

!    VARIABLES DE CONTROLE DE LA PHYSIQUE AROME :
!    CONTROL VARIABLES FOR THE AROME PHYSICS:

!    LMPA    : CLE GENERALE POUR LA PHYSIQUE AROME
!            : GLOBAL SWITCH FOR AROME PHYSICS
!    LMICRO  : CLE D'APPEL DU SCHEMA MICROPHYSIQUE ICE3
!              SWITCH FOR THE "ICE3" MICROPHYSICS SCHEME
!    LTURB   : CLE D'APPEL DU SCHEMA DE TURBULENCE
!              SWITCH FOR THE TURBULENCE SCHEME
!    LMSE    : CLE D'APPEL DU SHEMA DE SURFACE EXTERNALISE
!              SWITCH FOR THE EXTERNALIZED SURFACE SCHEME
!    LKFBCONV: Control key to KFB convection scheme (deep and/or shallow)
!    LKFBD   : Control key to KFB deep convection 
!    LKFBS   : Control key to KFB shallow convection 
!    LUSECHEM: Contol key for calling the Gas Chemistry scheme
!    LORILAM : Contol key for calling the Aerosol Chemistry scheme
!    LRDUST   : Contol key for calling the desertic aerosols
!    LBUFLUX  : If TRUE fluxes are calculated in AROEND_BUDGET,
!               if FALSE, tendencies remain 
!   CCOUPLING : Type of SURFEX coupling. E - explicit, I - implicit 

LOGICAL :: LMPA
LOGICAL :: LMICRO
LOGICAL :: LTURB
LOGICAL :: LMSE
LOGICAL :: LKFBCONV
LOGICAL :: LKFBD,LKFBS
LOGICAL :: LUSECHEM
LOGICAL :: LORILAM
LOGICAL :: LRDUST
LOGICAL :: LBUFLUX
CHARACTER(LEN=1) :: CCOUPLING
!    -------------------------------------------------------------------
!$OMP THREADPRIVATE(ccoupling,lbuflux,lkfbconv,lkfbd,lkfbs,lmicro,lmpa,lmse,lorilam,lrdust,lturb,lusechem)
END MODULE YOMARPHY

