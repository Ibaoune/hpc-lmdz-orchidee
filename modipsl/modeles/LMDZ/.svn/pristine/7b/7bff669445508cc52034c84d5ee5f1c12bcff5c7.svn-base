MODULE YOEAERATM

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAERATM* - CONTROL PARAMETERS FOR AEROSOLS IN THE ATMOSPHERE
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: INDBG
INTEGER(KIND=JPIM) :: NDD1, NSS1

REAL(KIND=JPRB) :: RMFMIN
REAL(KIND=JPRB) :: RMASSE(15)

REAL(KIND=JPRB) :: REPSCAER

LOGICAL :: LAERCLIMG, LAERCLIMZ, LAERCLIST, LAERDRYDP, LAERGBUD 
LOGICAL :: LAERNGAT , LAERPRNT , LAERSCAV , LAERSEDIM, LAERSURF, LAER6SDIA

!     ------------------------------------------------------------------
! NDD1       : location of first bin for desert dust
! NSS1       : location of first bin for sea-salt
! RMFMIN     : minimum mass flux for convective aerosol transport
! RMASSE     : Molar mass: N.B.: either g/mol or Avogadro number
!
! REPSCAER   : security on aerosol concentration: always >= 1.E-15
!
! LAERCLIMG  : .T. to start prognostic aerosols with geographical monthly 
!                  mean climatology
! LAERCLIMZ  : .T. to start prognostic aerosols with zonal annual mean 
!                  climatology
! LAERCLIST  : .T. to start prognostic aerosols with geographical monthly 
!                  mean climatology for background stratospheric only
! LAERDRYDP  : .T. dry deposition is active
! LAERGBUD   : .T. "clean" global budget for aerosols
! LAERNGAT   : .T. prevents negative aerosol concentrations
! LAERPRNT   : .T. for on-line prints of aerosol physical processes (NSTEP<=1)
! LAERSCAV   : .T. in-cloud and below cloud scavenging is active
! LAERSEDIM  : .T. sedimentation is active
! LAERSURF   : .T. if surface emissions
! LAER6SDIA  : .T. if radiance diagnostics with 6S
!     ------------------------------------------------------------------


!$OMP THREADPRIVATE(indbg,laer6sdia,laerclimg,laerclimz,laerclist,laerdrydp)
!$OMP THREADPRIVATE(laergbud,laerngat,laerprnt,laerscav,laersedim,laersurf)
!$OMP THREADPRIVATE(ndd1,nss1,repscaer,rmasse,rmfmin)

END MODULE YOEAERATM

