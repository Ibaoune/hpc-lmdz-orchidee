! $Id: aerophys.F90 3526 2019-05-28 13:00:44Z fairhead $
!
MODULE aerophys
!
  IMPLICIT NONE
!
  REAL,PARAMETER                         :: ropx=1500.0              ! default aerosol particle mass density [kg/m3]
  REAL,PARAMETER                         :: dens_aer_dry=1848.682308 ! dry aerosol particle mass density at T_0=293K[kg/m3]
  REAL,PARAMETER                         :: dens_aer_ref=1483.905336 ! aerosol particle mass density with 75% H2SO4 at T_0=293K[kg/m3]
  REAL,PARAMETER                         :: mdwmin=0.002e-6          ! dry diameter of smallest aerosol particles [m]
  REAL,PARAMETER                         :: V_rat=2.0                ! volume ratio of neighboring size bins
  REAL,PARAMETER                         :: mfrac_H2SO4=0.75         ! default mass fraction of H2SO4 in the aerosol
  REAL, PARAMETER                        :: mAIRmol=28.949*1.66E-27  ! Average mass of an air molecule [kg]
  REAL, PARAMETER                        :: mH2Omol=18.016*1.66E-27  ! Mass of an H2O molecule [kg]
  REAL, PARAMETER                        :: mH2SO4mol=98.082*1.66E-27! Mass of an H2SO4 molecule [kg]
  REAL, PARAMETER                        :: mSO2mol=64.06*1.66E-27   ! Mass of an SO2 molecule [kg]
  REAL, PARAMETER                        :: mSatom=32.06*1.66E-27    ! Mass of a S atom [kg]
  REAL, PARAMETER                        :: mOCSmol=60.07*1.66E-27   ! Mass of an OCS molecule [kg]
!
END MODULE aerophys
