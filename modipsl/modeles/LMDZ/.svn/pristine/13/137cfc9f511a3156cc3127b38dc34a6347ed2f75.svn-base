MODULE YOEAERSNK

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAERSNK* - CONTROL OPTIONS FOR AEROSOLS' SINKS
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: R_R, R_S
REAL(KIND=JPRB) :: RALPHAR(15), RALPHAS(15)
REAL(KIND=JPRB) :: RFRBC(2)
REAL(KIND=JPRB) :: RFRIF 
REAL(KIND=JPRB) :: RFRDD(3), RFRSS(3) 
REAL(KIND=JPRB) :: RFROM(2) 
REAL(KIND=JPRB) :: RFRSO4
REAL(KIND=JPRB) :: RFRAER, RFRGAS 

INTEGER(KIND=JPIM) :: NBRH
REAL(KIND=JPRB) :: RRHMAX, RRHTAB(12), RRHO_SS(12), RSSGROW(12)
REAL(KIND=JPRB) :: RMMD_SS(3)
REAL(KIND=JPRB) :: RMMD_DD(3), RRHO_DD(3)
REAL(KIND=JPRB) :: RHO_WAT, RHO_ICE

REAL(KIND=JPRB) :: RVDPOCE(15), RVDPSIC(15), RVDPLND(15), RVDPLIC(15)
REAL(KIND=JPRB) :: RVSEDOCE(15),RVSEDSIC(15),RVSEDLND(15),RVSEDLIC(15)
!     ------------------------------------------------------------------
! R_R        : mean radius for rain drops (m)
! R_S        : mean radius for snow crystals (m)
! RALPHAR    : impaction coefficients for rain
! RALPHAS    : IMPACTION coefficients for snow
! dissolution constants
! RFRBCn     : for black carbon (BC1: hydrophobic, BC2: hydrophylic)
! RFRDDn     : for dust     (3 bins: 0.03 - 0.55 - 0.9 - 20 um)
! RFRIF      : for inorganic fraction (i.e., fly-ash)
! RFROMn     : for organic matter (OM1: hydrophobic, OM2: hydrophylic)
! RFRSSn     : for sea-salt (3 bins: 0.03 - 0.5 - 5 - 20 um)

! re-evaporation constants
! RFRAER     : for aerosols
! RFRGAS     : for gases

! RRHMAX     : 
! RRHTAB     : 

! dry deposition velocity
! RVDPOCE    :
! RVDPSIC    :
! RVDPLND    :
! RVDPLIC    :

! sedimentation (gravitational settling) velocity
! RVSEDOCE   : for ocean
! RVSEDSIC   :
! RVSEDLND   :
! RVSEDLIC   :
!     -----------------------------------------------------------------

!$OMP THREADPRIVATE(nbrh,r_r,r_s,ralphar,ralphas,rfraer,rfrbc,rfrdd,rfrgas)
!$OMP THREADPRIVATE(rfrif,rfrom,rfrso4,rfrss,rho_ice,rho_wat,rmmd_dd,rmmd_ss)
!$OMP THREADPRIVATE(rrhmax,rrho_dd,rrho_ss,rrhtab,rssgrow,rvdplic,rvdplnd)
!$OMP THREADPRIVATE(rvdpoce,rvdpsic,rvsedlic,rvsedlnd,rvsedoce,rvsedsic)

END MODULE YOEAERSNK

