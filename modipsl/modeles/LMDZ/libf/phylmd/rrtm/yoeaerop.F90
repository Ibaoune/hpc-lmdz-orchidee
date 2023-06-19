MODULE YOEAEROP

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAEROP* - OPTICAL PROPERTIES FOR PROGNOSTIC AEROSOLS
!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ALF_BC(16)     , ASY_BC(16)     , OMG_BC(16)
REAL(KIND=JPRB) :: ALF_DD(3,16)   , ASY_DD(3,16)   , OMG_DD(3,16)
REAL(KIND=JPRB) :: ALF_FA(16)     , ASY_FA(16)     , OMG_FA(16)
REAL(KIND=JPRB) :: ALF_OM(12,16)  , ASY_OM(12,16)  , OMG_OM(12,16)
REAL(KIND=JPRB) :: ALF_SS(12,16,3), ASY_SS(12,16,3), OMG_SS(12,16,3)
REAL(KIND=JPRB) :: ALF_SU(12,16)  , ASY_SU(12,16)  , OMG_SU(12,16)

!     ------------------------------------------------------------------
! 3 refers to 3 bins, 16 to 16 SW channel radiances, 12 to 12 reference RH
! BC is for black carbon
! DD is for desert dust
! FA is fort flying ash
! OM is for organic matter
! SS is for sea-salt
! SU is for sulfate

! ALF is alpha , the mass extinction coefficient   in m2/g
! ASY is g     , the assymetry factor		     ND
! OMG is pizero, the single scattering albedo	     ND
!     ------------------------------------------------------------------

!$OMP THREADPRIVATE(alf_bc,alf_dd,alf_fa,alf_om,alf_ss,alf_su,asy_bc,asy_dd)
!$OMP THREADPRIVATE(asy_fa,asy_om,asy_ss,asy_su,omg_bc,omg_dd,omg_fa,omg_om,omg_ss,omg_su)

END MODULE YOEAEROP

