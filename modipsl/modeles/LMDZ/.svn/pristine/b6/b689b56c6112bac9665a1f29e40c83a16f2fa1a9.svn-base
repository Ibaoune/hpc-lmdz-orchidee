! $Id$
!
MODULE aero_mod

  IMPLICIT NONE

! Declaration des indices pour les aerosols 

! 1/ Total number of aerosols for which an aerosol optical depth is provided
!--strat aerosols are only prescribed naero_tot = 10 ==> 11 
!--adding nitrate naero_tot = 14 OB

  INTEGER, PARAMETER :: naero_tot = 14

! Identification number used in aeropt_2bands and aeropt_5wv
! corresponding to naero_tot
  INTEGER, PARAMETER :: id_ASBCM_phy    = 1
  INTEGER, PARAMETER :: id_ASPOMM_phy   = 2
  INTEGER, PARAMETER :: id_ASSO4M_phy   = 3
  INTEGER, PARAMETER :: id_CSSO4M_phy   = 4
  INTEGER, PARAMETER :: id_SSSSM_phy    = 5
  INTEGER, PARAMETER :: id_CSSSM_phy    = 6
  INTEGER, PARAMETER :: id_ASSSM_phy    = 7
  INTEGER, PARAMETER :: id_CIDUSTM_phy  = 8
  INTEGER, PARAMETER :: id_AIBCM_phy    = 9
  INTEGER, PARAMETER :: id_AIPOMM_phy   = 10
  INTEGER, PARAMETER :: id_ASNO3M_phy   = 11
  INTEGER, PARAMETER :: id_CSNO3M_phy   = 12
  INTEGER, PARAMETER :: id_CINO3M_phy   = 13
  INTEGER, PARAMETER :: id_STRAT_phy    = 14

! Corresponding names for the aerosols
  CHARACTER(len=7),DIMENSION(naero_tot), PARAMETER :: name_aero_tau=(/&
       "ASBCM  ", &
       "ASPOMM ", &
       "ASSO4M ", &
       "CSSO4M ", &
       "SSSSM  ", &
       "CSSSM  ", &
       "ASSSM  ", &
       "CIDUSTM", &
       "AIBCM  ", &
       "AIPOMM ", &
       "ASNO3M ", & 
       "CSNO3M ", &
       "CINO3M ", &
       "STRAT  " /)

! 2/ Total number of aerosols for which an aerosol mass is provided

  INTEGER, PARAMETER :: naero_spc = 13

! Corresponding names for the aerosols
  CHARACTER(len=7),DIMENSION(naero_spc), PARAMETER :: name_aero=(/&
       "ASBCM  ", &
       "ASPOMM ", &
       "SO4    ", &
       "CSSO4M ", &
       "SSSSM  ", &
       "CSSSM  ", &
       "ASSSM  ", &
       "CIDUSTM", &
       "AIBCM  ", &
       "AIPOMM " ,&
       "ASNO3M ", & 
       "CSNO3M ", &
       "CINO3M " /)

! 3/ Number of aerosol groups
  INTEGER, PARAMETER :: naero_grp = 13
  ! if info_trac = inca
  ! 1 = ZERO    
  ! 2 = AER total    
  ! 3 = NAT    
  ! 4 = BC    
  ! 5 = SO4    
  ! 6 = POM    
  ! 7 = DUST    
  ! 8 = SS    
  ! 9 = FNO3    
  ! 10 = DNO3
  ! 11 = SNO3
  ! 12 = SOAA
  ! 13 = SOAB 
  ! else 
  ! 1 = ZERO    
  ! 2 = AER total    
  ! 3 = NAT    
  ! 4 = BC    
  ! 5 = SO4    
  ! 6 = POM    
  ! 7 = DUST    
  ! 8 = SS    
  ! 9 = NO3    

! Number of diagnostics wavelengths (5 SW + 1 LW @ 10 um)
  INTEGER, PARAMETER :: nwave_sw = 5
  INTEGER, PARAMETER :: nwave_lw = 1
  INTEGER, PARAMETER :: nwave = nwave_sw + nwave_lw

! Number of modes spectral bands
  INTEGER, parameter :: nbands = 2
  INTEGER, parameter :: nbands_sw_rrtm = 6
  INTEGER, parameter :: nbands_lw_rrtm = 16

END MODULE aero_mod
