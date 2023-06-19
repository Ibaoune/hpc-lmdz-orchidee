MODULE YOMARAR

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!*
!     ------------------------------------------------------------------

!     VARIABLES pour utiliser la PHYSIQUE de meso_NH :
!     VARIABLES to use the MESO-NH physics:

!     KSPLITR :  Number of small time step
!                integration for  rain
!                sedimendation

INTEGER(KIND=JPIM) :: NSPLITR     ! ??? (missing comment)
INTEGER(KIND=JPIM) :: NRR, NRRL, NRRI   !number of microphysical species
INTEGER(KIND=JPIM) :: NSV         !number of passiv variables in MesoNH,
                                  ! always 0 in AROME
INTEGER(KIND=JPIM) :: NSWB_MNH    !number of SW bands for surface
                                  ! (must be equal to NSW !!)
INTEGER(KIND=JPIM) :: NGPAR       !number of fields in the buffer containing
                                  ! the 2D pseudo-historical variables.
INTEGER(KIND=JPIM) :: MINPRR      !pointer on INPRR 
INTEGER(KIND=JPIM) :: MACPRR      !pointer on ACPRR 
INTEGER(KIND=JPIM) :: MINPRS      !pointer on INPRS
INTEGER(KIND=JPIM) :: MACPRS      !pointer on ACPRS
INTEGER(KIND=JPIM) :: MINPRG      !pointer on INPRG
INTEGER(KIND=JPIM) :: MACPRG      !pointer on ACPRG
INTEGER(KIND=JPIM) :: MALBDIR     !pointer on ALBDIR
INTEGER(KIND=JPIM) :: MALBSCA     !pointer on ALBSCA
INTEGER(KIND=JPIM) :: MRAIN       !pointer on surface rain
INTEGER(KIND=JPIM) :: MSNOW       !pointer on surface snow
INTEGER(KIND=JPIM) :: MGZ0        !pointer on GZ0 
INTEGER(KIND=JPIM) :: MGZ0H       !pointer on GZ0H 
INTEGER(KIND=JPIM) :: MVQS        !pointer on surface moisture
INTEGER(KIND=JPIM) :: MVTS        !pointer on surface temperature 
INTEGER(KIND=JPIM) :: MVEMIS      !pointer on surface emissivity
INTEGER(KIND=JPIM) :: MSWDIR      !pointer on SW direct surface flux
INTEGER(KIND=JPIM) :: MSWDIF      !pointer on SW surface diffuse flux
INTEGER(KIND=JPIM) :: MSFU        !pointer on U surface flux 
INTEGER(KIND=JPIM) :: MSFV        !pointer on V surface flux
INTEGER(KIND=JPIM) :: MSFRV       !pointer on Qv  surface flux
INTEGER(KIND=JPIM) :: MSFTH       !pointer on Theta surface flux
INTEGER(KIND=JPIM) :: MFRTHDS     !pointer on IR downward surface flux
INTEGER(KIND=JPIM) :: MUM         !pointer on U surface wind
INTEGER(KIND=JPIM) :: MVM         !pointer on V surface wind
INTEGER(KIND=JPIM) :: MTM         !pointer on T at first model layer
INTEGER(KIND=JPIM) :: MQVM        !pointer on QV at first model layer
INTEGER(KIND=JPIM) :: MPABSM      !pointer on Pressure at first model layer
INTEGER(KIND=JPIM) :: MPSURF      !pointer on surface pressure
INTEGER(KIND=JPIM) :: MZZ         !pointer on height of first model layer
INTEGER(KIND=JPIM) :: MRHODREF    !pointer on RHODREF
INTEGER(KIND=JPIM) :: MSFSV       !pointer on surf. fluxes of scalars

REAL(KIND=JPRB), DIMENSION(:), ALLOCATABLE  :: XSW_BANDS  !SW spectral bands
					      ! for ext. surface scheme
LOGICAL :: LOSUBG_COND ! see OSUBG_COND in mesoNH
LOGICAL :: LOSUBG_AUCV ! see OSUBG_AUCV in mesoNH
LOGICAL :: LOWARM      ! see OWARM in mesoNH
LOGICAL :: LOSIGMAS    ! see OSIGMAS in mesoNH

! * for the squall line case:
LOGICAL :: LSQUALL ! use for the squall line case
INTEGER(KIND=JPIM) :: NREFROI1 !starting point for cooling 
INTEGER(KIND=JPIM) :: NREFROI2 !end point for cooling
REAL(KIND=JPRB) :: VSQUALL ! mean velocity displacement of the squall line.

! * for the MESO-NH physics printings:
INTEGER(KIND=JPIM) :: NPTP ! index in NPROMA paquet where the print will be done
INTEGER(KIND=JPIM) :: NPRINTFR !frequency of physical prints in apl_arome
                               ! (default 1/h) 

! * for the surface diagnostics.
INTEGER(KIND=JPIM) :: NDIAGFR  ! frequency of surface scheme diagnostics

!* for other diagnostics
!  wmax per vertical level 
LOGICAL :: LDIAGWMAX !activate print of WMAX in apl_arome
INTEGER(KIND=JPIM) :: NDIAGWMAX ! frequency of preceding prints (in time step)

!* for chemical scheme
! time step factor
INTEGER(KIND=JPIM) :: NDTCHEM
!* for MNH budget anlysis
LOGICAL :: LAROBU_ENABLE

! for budgets and DDH, number of processes in budget arrays 
INTEGER(KIND=JPIM),ALLOCATABLE ::  NBUPROC(:)
INTEGER(KIND=JPIM),ALLOCATABLE ::  NJBUDG1(:)
INTEGER(KIND=JPIM),ALLOCATABLE ::  NJBUDG2(:)

!* number of budgets used in AROME, without pasive scalars
INTEGER(KIND=JPIM), PARAMETER ::  JPAROBUD = 11

!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(larobu_enable,ldiagwmax,losigmas,losubg_aucv,losubg_cond,lowarm,lsquall)
!$OMP THREADPRIVATE(macprg,macprr,macprs,malbdir,malbsca,mfrthds,mgz0,mgz0h,minprg,minprr,minprs)
!$OMP THREADPRIVATE(mpabsm,mpsurf,mqvm,mrain,mrhodref,msfrv,msfsv,msfth,msfu,msfv,msnow,mswdif)
!$OMP THREADPRIVATE(mswdir,mtm,mum,mvemis,mvm,mvqs,mvts,mzz,ndiagfr,ndiagwmax,ndtchem,ngpar)
!$OMP THREADPRIVATE(nprintfr,nptp,nrefroi1,nrefroi2,nrr,nrri,nrrl,nsplitr,nsv,nswb_mnh,vsquall)
!$OMP THREADPRIVATE(nbuproc,njbudg1,njbudg2,xsw_bands)
END MODULE YOMARAR
