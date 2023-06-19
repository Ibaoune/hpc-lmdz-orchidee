MODULE YOMSIMPHL

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!   -----------------------------------------------------------------

!*    with logical switches for simplified physical 
!           parametrization

!    LSIMPH   : switch for simplified physical parametrization 
!    LTRAJPS  : switch for write out and read the trajectory at 
!               t-dt for simplified physical computations in file
!    LTRAJPST : switch for write out and read the trajectory at 
!               t-dt for fluxes and tendencies in file
!    LSMOOTHD : smoothing some functions in direct computation
!    LSMOOTHA : smoothing some functions in TL and AD computations
!    LSMOOTHB : modifications in simplified physical parametrizations
!               to stabilize the TL and AD codes

!    LCVRASP  : key for calling deep convection of simplified
!               phys. parametrization (ACCONV)
!    LGWDSP   : key for calling the gravity wave drag of simpl.
!               phys. parametrization (ACDRAGL)
!    LRAYSP   : key for calling radiation scheme of simplified
!               phys. parametrization (ACRADS)
!    LSTRASP  : key for calling stratiform precipitation of simpl.
!               phys. parametrization (ACQWLSR)
!    LVDIFSP  : key for calling the vertical turbulent diffusion
!               of simpl. phys. param. (ACTSEC,ACDIFSP)
!    LRRMESSP : key for calling the mesospheric drag
!               of simpl. phys. param. (ACDRME)
!    LCLOUDS  : key for cloud parametrization

LOGICAL :: LSIMPH
LOGICAL :: LTRAJPS
LOGICAL :: LTRAJPST
LOGICAL :: LSMOOTHD
LOGICAL :: LSMOOTHA
LOGICAL :: LSMOOTHB
LOGICAL :: LCVRASP
LOGICAL :: LGWDSP
LOGICAL :: LRAYSP
LOGICAL :: LSTRASP
LOGICAL :: LVDIFSP
LOGICAL :: LRRMESSP
LOGICAL :: LCLOUDS

!   ----------------------------------------------------------------
!$OMP THREADPRIVATE(lclouds,lcvrasp,lgwdsp,lraysp,lrrmessp,lsimph,lsmootha,lsmoothb,lsmoothd,lstrasp)
!$OMP THREADPRIVATE(ltrajps,ltrajpst,lvdifsp)
END MODULE YOMSIMPHL
