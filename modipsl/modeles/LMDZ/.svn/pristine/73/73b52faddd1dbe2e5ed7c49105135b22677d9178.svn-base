! A.Idelkadi sept 2013 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Module pour declarer et initialiser les parametres de controle des fichiers de sorties et des champs a sortir
!! La routine cosp_output_open (appelee 1 seule fois dans phy_cosp.F90) permet :
!! de creer les fichiers avec leurs grilles horizontales et verticales
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE cosp_output_mod

  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  use MOD_COSP_Modis_Simulator, only : cosp_modis
  use mod_modis_sim, only : numMODISReffIceBins, reffICE_binCenters, &
                            numMODISReffLiqBins, reffLIQ_binCenters

     IMPLICIT NONE
! cosp_output_mod
      INTEGER, PRIVATE             :: i
!!!!!!! Controle des fichier de sorties Cosp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      LOGICAL, DIMENSION(3), SAVE  :: cosp_outfilekeys
      INTEGER, DIMENSION(3), SAVE  :: cosp_nidfiles
!$OMP THREADPRIVATE(cosp_outfilekeys, cosp_nidfiles)
      INTEGER, DIMENSION(3), SAVE  :: nhoricosp,nvert,nvertmcosp,nvertcol,nvertbze, &
                                      nvertsratio,nvertisccp,nvertp,nverttemp,nvertmisr, &
                                      nvertReffIce,nvertReffLiq,nverttau
      REAL, DIMENSION(3), SAVE                :: zoutm_cosp
!$OMP THREADPRIVATE(nhoricosp, nvert,nvertmcosp,nvertcol,nvertsratio,nvertbze,nvertisccp,nvertp,zoutm_cosp,nverttemp,nvertmisr)
!$OMP THREADPRIVATE(nvertReffIce,nvertReffLiq,nverttau)
      REAL, SAVE                   :: zdtimemoy_cosp
!$OMP THREADPRIVATE(zdtimemoy_cosp) 
      CHARACTER(LEN=20), DIMENSION(3), SAVE  :: cosp_outfiletypes
      CHARACTER(LEN=20), DIMENSION(3), SAVE  :: cosp_outfilenames
      REAL, DIMENSION(3), SAVE               :: cosp_ecritfiles 
!$OMP THREADPRIVATE(cosp_outfiletypes, cosp_outfilenames, cosp_ecritfiles)

!!!!  Controle des variables a sortir dans les fichiers !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE ctrl_outcosp
     LOGICAL,DIMENSION(3)                 :: cles             !!! Sortir ou non le champs
     CHARACTER(len=20)                    :: name       
     CHARACTER(len=150)                   :: description      !!! Nom
     CHARACTER(len=20)                    :: unit             !!! Unite 
     CHARACTER(len=20),DIMENSION(3)  :: cosp_typeecrit        !!! Operation (ave, inst, ...)
  END TYPE ctrl_outcosp

! CALIPSO vars
  TYPE(ctrl_outcosp), SAVE :: o_cllcalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cllcalipso", "Lidar Low-level Cloud Fraction", "1", (/ ('', i=1, 3) /))                                   
  TYPE(ctrl_outcosp), SAVE :: o_clmcalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clmcalipso", "Lidar Mid-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clhcalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clhcalipso", "Lidar High-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltcalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cltcalipso", "Lidar Total Cloud Fraction", "1", (/ ('', i=1, 3) /)) 
  TYPE(ctrl_outcosp), SAVE :: o_clcalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipso", "Lidar Cloud Fraction (532 nm)", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cfad_lidarsr532 = ctrl_outcosp((/ .FALSE., .FALSE., .FALSE. /), &
         "cfad_lidarsr532", "Lidar Scattering Ratio CFAD (532 nm)", "1", (/ ('', i=1, 3) /))   
  TYPE(ctrl_outcosp), SAVE :: o_parasol_refl = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "parasol_refl", "PARASOL-like mono-directional reflectance","1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_parasol_crefl = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &              
         "parasol_crefl", "PARASOL-like mono-directional reflectance (integral)","1", (/ ('', i=1, 3) /))                  
  TYPE(ctrl_outcosp), SAVE :: o_Ncrefl = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "Ncrefl", "Nb PARASOL-like mono-directional reflectance (integral)","1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_atb532 = ctrl_outcosp((/ .FALSE., .FALSE., .FALSE. /), &
         "atb532", "Lidar Attenuated Total Backscatter (532 nm)","1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_beta_mol532 = ctrl_outcosp((/ .FALSE., .FALSE., .FALSE. /), &
         "beta_mol532", "Lidar Molecular Backscatter (532 nm)","m-1 sr-1", (/ ('', i=1, 3) /))
!! AI  11 2015
  TYPE(ctrl_outcosp), SAVE :: o_cllcalipsoice = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), & 
         "cllcalipsoice", "CALIPSO Ice-Phase Low Level Cloud Fraction", "%", (/ ('', i=1, 3) /))  
  TYPE(ctrl_outcosp), SAVE :: o_cllcalipsoliq = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cllcalipsoliq", "CALIPSO Liq-Phase Low Level Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clmcalipsoice = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clmcalipsoice", "CALIPSO Ice-Phase Mid Level Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clmcalipsoliq = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clmcalipsoliq", "CALIPSO Liq-Phase Mid Level Cloud Fraction", "%", (/ ('', i=1, 3) /))	 	 
  TYPE(ctrl_outcosp), SAVE :: o_clhcalipsoice = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), & 
         "clhcalipsoice", "CALIPSO Ice-Phase High Level Cloud Fraction", "%", (/ ('', i=1, 3) /))  
  TYPE(ctrl_outcosp), SAVE :: o_clhcalipsoliq = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clhcalipsoliq", "CALIPSO Liq-Phase High Level Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltcalipsoice = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cltcalipsoice", "CALIPSO Ice-Phase Tot Level Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltcalipsoliq = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cltcalipsoliq", "CALIPSO Liq-Phase Tot Level Cloud Fraction", "%", (/ ('', i=1, 3) /))	 	 
  TYPE(ctrl_outcosp), SAVE :: o_cllcalipsoun = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), & 
         "cllcalipsoun", "CALIPSO Undefined-Phase Low Level Cloud Fraction", "%", (/ ('', i=1, 3) /))  
  TYPE(ctrl_outcosp), SAVE :: o_clmcalipsoun = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clmcalipsoun", "CALIPSO Undefined-Phase Mid Level Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clhcalipsoun = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clhcalipsoun", "CALIPSO Undefined-Phase High Level Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltcalipsoun = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cltcalipsoun", "CALIPSO Undefined-Phase Tot Level Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsoice = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsoice", "Lidar Ice-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsoliq = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsoliq", "Lidar Liq-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsoun = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsoun", "Lidar Undef-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))	 
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsotmpice = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsotmpice", "Lidar Ice-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsotmpliq = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsotmpliq", "Lidar Liq-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsotmpun = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsotmpun", "Lidar Undef-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsotmp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsotmp", "Lidar Cloud Fraction", "%", (/ ('', i=1, 3) /))

  TYPE(ctrl_outcosp), SAVE :: o_clopaquecalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &  !OPAQ
         "clopaquecalipso", "Lidar Opaque Cloud Fraction", "%", (/ ('', i=1, 3) /))             !OPAQ
  TYPE(ctrl_outcosp), SAVE :: o_clthincalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &    !OPAQ
         "clthincalipso", "Lidar Thin Cloud Fraction", "%", (/ ('', i=1, 3) /))                 !OPAQ
  TYPE(ctrl_outcosp), SAVE :: o_clzopaquecalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), & !OPAQ
         "clzopaquecalipso", "Lidar mean opacity altitude", "m", (/ ('', i=1, 3) /))            !OPAQ
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsoopaque = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &  !OPAQ
         "clcalipsoopaque", "Lidar Opaque profile Cloud Fraction", "%", (/ ('', i=1, 3) /))     !OPAQ
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsothin = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &    !OPAQ
         "clcalipsothin", "Lidar Thin profile Cloud Fraction", "%", (/ ('', i=1, 3) /))         !OPAQ
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsozopaque = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), & !OPAQ
         "clcalipsozopaque", "Lidar z_opaque Fraction", "%", (/ ('', i=1, 3) /))	        !OPAQ
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsoopacity = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), & !OPAQ
         "clcalipsoopacity", "Lidar opacity Fraction", "%", (/ ('', i=1, 3) /))	                !OPAQ

  TYPE(ctrl_outcosp), SAVE :: o_proftemp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &         !TIBO
         "proftemp", "Temperature profiles (40 lev)", "K", (/ ('', i=1, 3) /))                  !TIBO
  TYPE(ctrl_outcosp), SAVE :: o_profSR = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &           !TIBO
         "profSR", "Lidar Scattering Ratio profiles (532 nm)", "1", (/ ('', i=1, 3) /))         !TIBO

! Radar Cloudsat
  TYPE(ctrl_outcosp), SAVE :: o_cfadDbze94 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cfadDbze94", "CloudSat Radar Reflectivity CFAD", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_dbze94 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "dbze94", "CloudSat Radar Reflectivity", "%", (/ ('', i=1, 3) /))

! Calipso + Cloudsat
  TYPE(ctrl_outcosp), SAVE :: o_clcalipso2 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipso2", "CALIPSO Cloud Fraction Undetected by CloudSat", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltlidarradar = ctrl_outcosp((/ .TRUE., .TRUE.,.TRUE. /), &          
         "cltlidarradar", "Lidar and Radar Total Cloud Fraction", "%", (/ ('', i=1, 3) /))
     
! ISCCP vars
  TYPE(ctrl_outcosp), SAVE :: o_sunlit = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "sunlit", "1 for day points, 0 for nightime","1",(/ ('', i=1, 3) /))                   
  TYPE(ctrl_outcosp), SAVE :: o_clisccp2 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clisccp2", "Cloud Fraction as Calculated by the ISCCP Simulator","%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_boxtauisccp = ctrl_outcosp((/ .FALSE., .FALSE., .FALSE. /), &
         "boxtauisccp", "Optical Depth in Each Column as Calculated by the ISCCP Simulator","1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_boxptopisccp = ctrl_outcosp((/ .FALSE., .FALSE., .FALSE. /), &
         "boxptopisccp", "Cloud Top Pressure in Each Column as Calculated by the ISCCP Simulator","Pa", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_tclisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
          "tclisccp", "Total Cloud Fraction as Calculated by the ISCCP Simulator", "%", (/ ('', i=1, 3) /)) 
  TYPE(ctrl_outcosp), SAVE :: o_ctpisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
          "ctpisccp", "Mean Cloud Top Pressure as Calculated by the ISCCP Simulator", "Pa", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_tauisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
          "tauisccp", "Optical Depth as Calculated by the ISCCP Simulator", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_albisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
          "albisccp", "Mean Cloud Albedo as Calculated by the ISCCP Simulator", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_meantbisccp = ctrl_outcosp((/ .FALSE., .FALSE., .FALSE. /), &
          "meantbisccp", " Mean all-sky 10.5 micron brightness temperature as calculated &
           by the ISCCP Simulator","K", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_meantbclrisccp = ctrl_outcosp((/ .FALSE., .FALSE., .FALSE. /), &
          "meantbclrisccp", "Mean clear-sky 10.5 micron brightness temperature as calculated &
           by the ISCCP Simulator","K", (/ ('', i=1, 3) /))

! MISR simulator
  TYPE(ctrl_outcosp), SAVE :: o_clMISR = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clMISR", "Cloud Fraction as Calculated by the MISR Simulator","%", (/ ('', i=1, 3) /))

! MODIS simulator
  TYPE(ctrl_outcosp), SAVE :: o_cllmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cllmodis", "MODIS Low-level Cloud Fraction", "%", (/ ('', i=1, 3) /))                                   
  TYPE(ctrl_outcosp), SAVE :: o_clmmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clmmodis", "MODIS Mid-level Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clhmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clhmodis", "MODIS High-level Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cltmodis", "MODIS Total Cloud Fraction", "%", (/ ('', i=1, 3) /)) 
  TYPE(ctrl_outcosp), SAVE :: o_clwmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clwmodis", "MODIS Cloud Fraction water mean", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_climodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "climodis", "MODIS Cloud Fraction ice mean", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_tautmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "tautmodis", "MODIS Optical_Thickness_Total_Mean", "1", (/ ('', i=1, 3) /))                                   
  TYPE(ctrl_outcosp), SAVE :: o_tauwmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "tauwmodis", "MODIS Optical_Thickness_Water_Mean", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_tauimodis= ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "tauimodis", "MODIS Optical_Thickness_Ice_Mean", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_tautlogmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "tautlogmodis", "MODIS Optical_Thickness_Total_logMean", "1", (/ ('', i=1, 3) /))                                   
  TYPE(ctrl_outcosp), SAVE :: o_tauwlogmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "tauwlogmodis", "MODIS Optical_Thickness_Water_logMean", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_tauilogmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "tauilogmodis", "MODIS Optical_Thickness_Ice_logMean", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_reffclwmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "reffclwmodis", "Modis Cloud_Particle_Size_Water_Mean", "m", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_reffclimodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "reffclimodis", "Modis Cloud_Particle_Size_Ice_Mean", "m", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_pctmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "pctmodis", "Modis Cloud_Top_Pressure_Total_Mean", "Pa", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_lwpmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "lwpmodis", "Modis Liquid_Water_Path_Mean", "kg m-2", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_iwpmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "iwpmodis", "Modis Ice_Water_Path_Mean", "kg m-2", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clmodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clmodis", "MODIS Cloud Area Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_crimodis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "crimodis", "Optical_Thickness_vs_ReffIce from Modis", "%", (/ ('',i=1, 3) /))          
  TYPE(ctrl_outcosp), SAVE :: o_crlmodis = ctrl_outcosp((/ .TRUE., .TRUE.,.TRUE. /), &
         "crlmodis", "Optical_Thickness_vs_ReffLiq from Modis", "%", (/ ('',i=1, 3) /))         

! Rttovs simulator
  TYPE(ctrl_outcosp), SAVE :: o_tbrttov = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "tbrttov", "Rttovs Cloud Area Fraction", "%", (/ ('', i=1, 3) /))

! Scops and others
  TYPE(ctrl_outcosp), SAVE :: o_fracout = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "fracout", "Subcolumn output from SCOPS", "%", (/ ('', i=1, 3) /))

  LOGICAL, SAVE :: cosp_varsdefined = .FALSE. ! ug PAS THREADPRIVATE ET C'EST NORMAL
  REAL, SAVE  :: Cosp_fill_value
!$OMP THREADPRIVATE(Cosp_fill_value)
 

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Ouverture des fichier et definition des  axes!!!!!!!!
  !! histbeg, histvert
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  SUBROUTINE cosp_output_open(Nlevlmdz, Ncolumns, presnivs, dtime, freq_cosp, &
                              ok_mensuelCOSP, ok_journeCOSP, ok_hfCOSP, ok_all_xml,  &
                              ecrit_mth, ecrit_day, ecrit_hf, use_vgrid, vgrid, stlidar)

  USE iophy
  USE ioipsl
  USE phys_cal_mod
  USE time_phylmdz_mod, ONLY: day_ref, annee_ref, day_ini, start_time, itau_phy
  USE print_control_mod, ONLY: lunout

#ifdef CPP_XIOS
    ! ug Pour les sorties XIOS
    USE wxios
#endif

  IMPLICIT NONE

!!! Variables d'entree
  integer                  :: Nlevlmdz, Ncolumns      ! Number of levels
  real,dimension(Nlevlmdz) :: presnivs
  real                     :: dtime, freq_cosp, ecrit_day, ecrit_hf, ecrit_mth 
  logical                  :: ok_mensuelCOSP, ok_journeCOSP, ok_hfCOSP, use_vgrid, ok_all_xml                    
  type(cosp_vgrid)   :: vgrid   ! Information on vertical grid of stats
  type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator

!!! Variables locales
  integer                   :: idayref, iff, ii
  real                      :: zjulian,zjulian_start
  real,dimension(Ncolumns)  :: column_ax
  real,dimension(DBZE_BINS) ::  dbze_ax
  CHARACTER(LEN=20), DIMENSION(3)  :: chfreq = (/ '1day', '1d  ', '3h  ' /)            
  real,parameter,dimension(SR_BINS) :: sratio_ax = (/0.005, &
                                                  0.605,2.09,4.,6., & 
                                          8.5,12.5,17.5,22.5,27.5,35.,45.,55.,70.,50040./)

!!! Variables d'entree

#ifdef CPP_XIOS
    ! ug Variables utilisées pour récupérer le calendrier pour xios
    INTEGER :: x_an, x_mois, x_jour
    REAL :: x_heure
    INTEGER :: ini_an, ini_mois, ini_jour
    REAL :: ini_heure
#endif

    WRITE(lunout,*) 'Debut cosp_output_mod.F90'
    print*,'cosp_varsdefined',cosp_varsdefined
    ! Initialisations (Valeurs par defaut)

!! Definition valeurs axes
    do ii=1,Ncolumns
      column_ax(ii) = real(ii)
    enddo

    do i=1,DBZE_BINS
     dbze_ax(i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*(i - 0.5)
    enddo
 
    cosp_outfilenames(1) = 'histmthCOSP'
    cosp_outfilenames(2) = 'histdayCOSP'
    cosp_outfilenames(3) = 'histhfCOSP'

    cosp_outfiletypes(1) = 'ave(X)'
    cosp_outfiletypes(2) = 'ave(X)'
    cosp_outfiletypes(3) = 'ave(X)'

    cosp_outfilekeys(1) = ok_mensuelCOSP
    cosp_outfilekeys(2) = ok_journeCOSP
    cosp_outfilekeys(3) = ok_hfCOSP

    cosp_ecritfiles(1) = mth_len*86400.
    cosp_ecritfiles(2) = 1.*86400.
    cosp_ecritfiles(3) = 0.125*86400.

! Lecture des parametres dans output.def ou config.def

    CALL getin('cosp_outfilenames',cosp_outfilenames)
    CALL getin('cosp_outfilekeys',cosp_outfilekeys)
    CALL getin('cosp_ecritfiles',cosp_ecritfiles)
    CALL getin('cosp_outfiletypes',cosp_outfiletypes)

    WRITE(lunout,*)'cosp_outfilenames=',cosp_outfilenames
    WRITE(lunout,*)'cosp_outfilekeys=',cosp_outfilekeys
    WRITE(lunout,*)'cosp_ecritfiles=',cosp_ecritfiles
    WRITE(lunout,*)'cosp_outfiletypes=',cosp_outfiletypes
    
    idayref = day_ref
    CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
    CALL ymds2ju(annee_ref, 1, day_ini, start_time, zjulian_start)

#ifdef CPP_XIOS
    
! recuperer la valeur indefine Xios
!    CALL xios_get_field_attr("clcalipso",default_value=Cosp_fill_value)
!         Cosp_fill_value=missing_val
          Cosp_fill_value=0.
         print*,'Cosp_fill_value=',Cosp_fill_value
!    if (use_vgrid) then
!      print*,'vgrid%Nlvgrid, vgrid%z = ',vgrid%Nlvgrid, vgrid%z
        CALL wxios_add_vaxis("height", vgrid%Nlvgrid, vgrid%z)
     print*,'wxios_add_vaxis '
!    else
!         WRITE(lunout,*) 'wxios_add_vaxis "presnivs", vgrid%Nlvgrid ',vgrid%Nlvgrid
!        CALL wxios_add_vaxis("presnivs", vgrid%Nlvgrid, presnivs)
!    endif
    WRITE(lunout,*) 'wxios_add_vaxis height_mlev, Nlevlmdz ',Nlevlmdz
    CALL wxios_add_vaxis("height_mlev", Nlevlmdz, vgrid%mz)
    WRITE(lunout,*) 'wxios_add_vaxis sza, PARASOL_NREFL ',PARASOL_NREFL
    CALL wxios_add_vaxis("sza", PARASOL_NREFL, PARASOL_SZA)
    WRITE(lunout,*) 'wxios_add_vaxis pressure2 ',7
    CALL wxios_add_vaxis("pressure2", 7, ISCCP_PC)
    WRITE(lunout,*) 'wxios_add_vaxis column ',Ncolumns
    CALL wxios_add_vaxis("column", Ncolumns, column_ax)

! AI nov 2015
   CALL wxios_add_vaxis("temp", LIDAR_NTEMP, LIDAR_PHASE_TEMP)
   CALL wxios_add_vaxis("cth", MISR_N_CTH, MISR_CTH)
   CALL wxios_add_vaxis("dbze", DBZE_BINS, dbze_ax) 
   CALL wxios_add_vaxis("scatratio", SR_BINS, sratio_ax)
   CALL wxios_add_vaxis("ReffIce", numMODISReffIceBins, reffICE_binCenters)
   CALL wxios_add_vaxis("ReffLiq", numMODISReffLiqBins, reffLIQ_binCenters)
   print*,'reffICE_binCenters=',reffICE_binCenters
   CALL wxios_add_vaxis("tau", 7, ISCCP_TAU)

#endif
   
    zdtimemoy_cosp = freq_COSP         ! Frequence ou l on moyenne

    DO iff=1,3
       zoutm_cosp(iff) = cosp_ecritfiles(iff) ! Frequence ou l on ecrit en seconde

       IF (cosp_outfilekeys(iff)) THEN
           CALL histbeg_phy_all(cosp_outfilenames(iff),itau_phy,zjulian,&
             dtime,nhoricosp(iff),cosp_nidfiles(iff))
!           print*,'histbeg_phy nhoricosp(iff),cosp_nidfiles(iff)', &
!                    nhoricosp(iff),cosp_nidfiles(iff)

#ifdef CPP_XIOS
        IF (.not. ok_all_xml) then
         WRITE(lunout,*) 'wxios_add_file ',cosp_outfilenames(iff)
         CALL wxios_add_file(cosp_outfilenames(iff),chfreq(iff),10)
        ENDIF
#endif

#ifndef CPP_IOIPSL_NO_OUTPUT 
! Definition de l'axe vertical
       if (use_vgrid) then
! Axe vertical Cosp 40 niveaux (en m)
      CALL histvert(cosp_nidfiles(iff),"height","height","m",vgrid%Nlvgrid,vgrid%z,nvert(iff))
       else
! Axe vertical modele LMDZ presnivs
      CALL histvert(cosp_nidfiles(iff),"presnivs","Vertical levels","Pa",vgrid%Nlvgrid,presnivs,nvert(iff),"down")
       endif
! Axe vertical niveaux modele (en m)
      CALL histvert(cosp_nidfiles(iff),"height_mlev","height_mlev","m",Nlevlmdz,vgrid%mz,nvertmcosp(iff))

      CALL histvert(cosp_nidfiles(iff),"sza","solar_zenith_angle","degrees",PARASOL_NREFL,PARASOL_SZA,nvertp(iff))

      CALL histvert(cosp_nidfiles(iff),"pressure2","pressure","mb",7,ISCCP_PC,nvertisccp(iff),"down")

      CALL histvert(cosp_nidfiles(iff),"column","column","count",Ncolumns,column_ax,nvertcol(iff)) !DBUG

      CALL histvert(cosp_nidfiles(iff),"temp","temperature","C",LIDAR_NTEMP,LIDAR_PHASE_TEMP,nverttemp(iff))

      CALL histvert(cosp_nidfiles(iff),"cth","altitude","m",MISR_N_CTH,MISR_CTH,nvertmisr(iff))
  
      CALL histvert(cosp_nidfiles(iff),"ReffIce","Effective_particle_size_Ice","microns",numMODISReffIceBins, reffICE_binCenters, &
                    nvertReffIce(iff))                                         
     
      CALL histvert(cosp_nidfiles(iff),"ReffLiq","Effective_particle_size_Liq","microns",numMODISReffLiqBins, reffLIQ_binCenters, &                                  
                    nvertReffLiq(iff))

      CALL histvert(cosp_nidfiles(iff),"dbze","equivalent_reflectivity_factor","dBZ",DBZE_BINS,dbze_ax,nvertbze(iff))
     
      CALL histvert(cosp_nidfiles(iff),"scatratio","backscattering_ratio","1",SR_BINS,sratio_ax,nvertsratio(iff))

      CALL histvert(cosp_nidfiles(iff),"tau","cloud optical depth","1",7,ISCCP_TAU,nverttau(iff)) 
     
!!! Valeur indefinie en cas IOIPSL
     Cosp_fill_value=0.

#endif

      ENDIF
  ENDDO

    end SUBROUTINE cosp_output_open

 END MODULE cosp_output_mod
