! A.Idelkadi sept 2013 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Module pour declarer et initialiser les parametres de controle des fichiers de sorties et des champs a sortir
!! La routine cosp_output_open (appelee 1 seule fois dans phy_cosp.F90) permet :
!! de creer les fichiers avec leurs grilles horizontales et verticales
!
! R.Guzman jan 2019 (mise a jour pour COSPv2)
! On change le nom du module a "lmdz_cosp_output_mod" et celui de la routine a "lmdz_cosp_output_open"
! pour qu'on sache qu'il s'agit d'un module specifique a l'implementation de COSP dans LMDZ
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE lmdz_cosp_output_mod

     USE COSP_KINDS, ONLY: wp,dp
     IMPLICIT NONE
! cosp_output_mod
      INTEGER          :: i
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

! CALIPSO vars (43)
  TYPE(ctrl_outcosp), SAVE :: o_cllcalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cllcalipso", "CALIPSO Low-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clmcalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clmcalipso", "CALIPSO Mid-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clhcalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clhcalipso", "CALIPSO High-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltcalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cltcalipso", "CALIPSO Total Cloud Fraction", "1", (/ ('', i=1, 3) /)) 
  TYPE(ctrl_outcosp), SAVE :: o_clcalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipso", "CALIPSO Cloud Fraction (532 nm)", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cfadLidarsr532 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cfadLidarsr532", "CALIPSO Scattering Ratio CFAD (532 nm)", "1", (/ ('', i=1, 3) /))   
  TYPE(ctrl_outcosp), SAVE :: o_atb532 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "atb532", "CALIPSO Attenuated Total Backscatter (532 nm)","1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_lidarBetaMol532 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "lidarBetaMol532", "CALIPSO Molecular Backscatter (532 nm)","m-1 sr-1", (/ ('', i=1, 3) /))
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
         "clcalipsoice", "CALIPSO Ice-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsoliq = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsoliq", "CALIPSO Liq-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsoun = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsoun", "CALIPSO Undef-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))	 
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsotmpice = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsotmpice", "CALIPSO Ice-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsotmpliq = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsotmpliq", "CALIPSO Liq-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsotmpun = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsotmpun", "CALIPSO Undef-Phase Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsotmp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsotmp", "CALIPSO Cloud Fraction", "%", (/ ('', i=1, 3) /))

  TYPE(ctrl_outcosp), SAVE :: o_clopaquecalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clopaquecalipso", "CALIPSO Opaque Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clthincalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clthincalipso", "CALIPSO Thin Cloud Fraction", "%", (/ ('', i=1, 3) /))      
  TYPE(ctrl_outcosp), SAVE :: o_clzopaquecalipso = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clzopaquecalipso", "CALIPSO mean opacity altitude", "m", (/ ('', i=1, 3) /))    
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsoopaque = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsoopaque", "CALIPSO Opaque profile Cloud Fraction", "%", (/ ('', i=1, 3) /))  
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsothin = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), & 
         "clcalipsothin", "CALIPSO Thin profile Cloud Fraction", "%", (/ ('', i=1, 3) /))      
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsozopaque = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsozopaque", "CALIPSO z_opaque Fraction", "%", (/ ('', i=1, 3) /))	   
  TYPE(ctrl_outcosp), SAVE :: o_clcalipsoopacity = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipsoopacity", "CALIPSO opacity Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clopaquetemp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clopaquetemp", "CALIPSO Opaque Cloud Temperature", "K", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clthintemp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clthintemp", "CALIPSO Thin Cloud Temperature", "K", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clzopaquetemp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clzopaquetemp", "CALIPSO z_opaque Temperature", "K", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clopaquemeanz = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clopaquemeanz", "CALIPSO Opaque Cloud Altitude", "m", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clthinmeanz = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clthinmeanz", "CALIPSO Thin Cloud Altitude", "m", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clthinemis = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clthinemis", "CALIPSO Thin Cloud Emissivity", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clopaquemeanzse = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clopaquemeanzse", "CALIPSO Opaque Cloud Altitude with respect to Surface Elevation", "m", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clthinmeanzse = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clthinmeanzse", "CALIPSO Thin Cloud Altitude with respect to Surface Elevation", "m", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clzopaquecalipsose = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clzopaquecalipsose", "CALIPSO mean opacity Altitude with respect to Surface Elevation", "m", (/ ('', i=1, 3) /))

! GROUND LIDAR 532 nm vars (8)
  TYPE(ctrl_outcosp), SAVE :: o_cllgrLidar532 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cllgrLidar532", "GROUND LIDAR 532 nm Low-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clmgrLidar532 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clmgrLidar532", "GROUND LIDAR 532 nm Mid-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clhgrLidar532 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clhgrLidar532", "GROUND LIDAR 532 nm High-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltgrLidar532 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cltgrLidar532", "GROUND LIDAR 532 nm Total Cloud Fraction", "1", (/ ('', i=1, 3) /)) 
  TYPE(ctrl_outcosp), SAVE :: o_clgrLidar532 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clgrLidar532", "GROUND LIDAR 532 nm Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cfadLidarsr532gr = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cfadLidarsr532gr", "GROUND LIDAR 532 nm Scattering Ratio CFAD", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_atb532gr = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "atb532gr", "GROUND LIDAR 532 nm Attenuated Total Backscatter (532 nm)","m-1 sr-1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_lidarBetaMol532gr = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "lidarBetaMol532gr", "GROUND LIDAR 532 nm Molecular Backscatter (532 nm)","m-1 sr-1", (/ ('', i=1, 3) /))

! ATLID vars (8)
  TYPE(ctrl_outcosp), SAVE :: o_cllatlid = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cllatlid", "ATLID Low-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clmatlid = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clmatlid", "ATLID Mid-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_clhatlid = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clhatlid", "ATLID High-level Cloud Fraction", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltatlid = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cltatlid", "ATLID Total Cloud Fraction", "1", (/ ('', i=1, 3) /)) 
  TYPE(ctrl_outcosp), SAVE :: o_clatlid = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clatlid", "ATLID Cloud Fraction (532 nm)", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cfadLidarsr355 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cfadLidarsr355", "ATLID Scattering Ratio CFAD (532 nm)", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_atb355 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "atb355", "ATLID Attenuated Total Backscatter (532 nm)","m-1 sr-1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_lidarBetaMol355 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "lidarBetaMol355", "ATLID Molecular Backscatter (532 nm)","m-1 sr-1", (/ ('', i=1, 3) /))

! PARASOL vars (2)
  TYPE(ctrl_outcosp), SAVE :: o_parasolGrid_refl = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "parasolGrid_refl", "PARASOL Reflectance","1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_parasolPix_refl = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), & 
         "parasolPix_refl", "PARASOL Subcolumn Reflectance","1", (/ ('', i=1, 3) /))                  
!  TYPE(ctrl_outcosp), SAVE :: o_Ncrefl = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
!         "Ncrefl", "Nb PARASOL-like mono-directional reflectance (integral)","1", (/ ('', i=1, 3) /))


! Radar CloudSat vars (13)
  TYPE(ctrl_outcosp), SAVE :: o_cfadDbze94 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "cfadDbze94", "CloudSat Radar Reflectivity CFAD", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_dbze94 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "dbze94", "CloudSat Radar Reflectivity", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_ptradarflag0 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "ptradarflag0", "CloudSat Radar Precipitation frequency flag 0", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_ptradarflag1 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "ptradarflag1", "CloudSat Radar Precipitation frequency flag 1", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_ptradarflag2 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "ptradarflag2", "CloudSat Radar Precipitation frequency flag 2", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_ptradarflag3 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "ptradarflag3", "CloudSat Radar Precipitation frequency flag 3", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_ptradarflag4 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "ptradarflag4", "CloudSat Radar Precipitation frequency flag 4", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_ptradarflag5 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "ptradarflag5", "CloudSat Radar Precipitation frequency flag 5", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_ptradarflag6 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "ptradarflag6", "CloudSat Radar Precipitation frequency flag 6", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_ptradarflag7 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "ptradarflag7", "CloudSat Radar Precipitation frequency flag 7", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_ptradarflag8 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "ptradarflag8", "CloudSat Radar Precipitation frequency flag 8", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_ptradarflag9 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "ptradarflag9", "CloudSat Radar Precipitation frequency flag 9", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_radarpia = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "radarpia", "CloudSat Radar PIA", "1", (/ ('', i=1, 3) /))

! CALIPSO + CloudSat vars (4)
  TYPE(ctrl_outcosp), SAVE :: o_clcalipso2 = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clcalipso2", "CALIPSO Cloud Fraction Undetected by CloudSat", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltlidarradar = ctrl_outcosp((/ .TRUE., .TRUE.,.TRUE. /), &          
         "cltlidarradar", "CALIPSO and CloudSat Radar Total Cloud Fraction", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cloudsat_tcc = ctrl_outcosp((/ .TRUE., .TRUE.,.TRUE. /), &          
         "cloudsat_tcc", "CALIPSO and CloudSat Radar Total Cloud Fraction?", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cloudsat_tcc2 = ctrl_outcosp((/ .TRUE., .TRUE.,.TRUE. /), &          
         "cloudsat_tcc2", "CALIPSO and CloudSat Radar Total Cloud Fraction?", "%", (/ ('', i=1, 3) /))

! ISCCP vars (9) + sunlit (1)
!  TYPE(ctrl_outcosp), SAVE :: o_sunlit = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
!         "sunlit", "1 for day points, 0 for nightime", "1", (/ ('', i=1, 3) /))                   
  TYPE(ctrl_outcosp), SAVE :: o_clisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "clisccp", "Cloud Fraction as Calculated by the ISCCP Simulator","%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_boxtauisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "boxtauisccp", "Optical Depth in Each Column as Calculated by the ISCCP Simulator","1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_boxptopisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "boxptopisccp", "Cloud Top Pressure in Each Column as Calculated by the ISCCP Simulator","Pa", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_cltisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
          "cltisccp", "Total Cloud Fraction as Calculated by the ISCCP Simulator", "%", (/ ('', i=1, 3) /)) 
  TYPE(ctrl_outcosp), SAVE :: o_pctisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
          "pctisccp", "Mean Cloud Top Pressure as Calculated by the ISCCP Simulator", "Pa", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_tauisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
          "tauisccp", "Optical Depth as Calculated by the ISCCP Simulator", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_albisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
          "albisccp", "Mean Cloud Albedo as Calculated by the ISCCP Simulator", "1", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_meantbisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
          "meantbisccp", " Mean all-sky 10.5 micron brightness temperature as calculated &
           by the ISCCP Simulator","K", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_meantbclrisccp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
          "meantbclrisccp", "Mean clear-sky 10.5 micron brightness temperature as calculated &
           by the ISCCP Simulator","K", (/ ('', i=1, 3) /))

! MISR simulator var (3)
  TYPE(ctrl_outcosp), SAVE :: o_misr_fq = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "misr_fq", "MISR joint-PDF of cloud top pressure and optical depth","%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_misr_meanztop = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "misr_meanztop", "MISR Mean Cloud Top Height","%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_misr_cldarea = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "misr_cldarea", "MISR cloud cover","%", (/ ('', i=1, 3) /))

! MODIS simulator vars (20)
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
  TYPE(ctrl_outcosp), SAVE :: o_modis_ot_vs_ctp = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "modis_ot_vs_ctp", "MODIS joint-PDF of cloud top pressure and optical depth", "%", (/ ('', i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_modis_ot_vs_reffice = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "modis_ot_vs_reffice", "MODIS Joint-PDF of optical-depth and ice particle size", "%", (/ ('',i=1, 3) /))
  TYPE(ctrl_outcosp), SAVE :: o_modis_ot_vs_reffliq = ctrl_outcosp((/ .TRUE., .TRUE.,.TRUE. /), &
         "modis_ot_vs_reffliq", "MODIS Joint-PDF of optical-depth and liquid particle size", "%", (/ ('',i=1, 3) /))

! Rttovs simulator var (1)
  TYPE(ctrl_outcosp), SAVE :: o_tbrttov = ctrl_outcosp((/ .TRUE., .TRUE., .TRUE. /), &
         "tbrttov", "Rttovs Cloud Area Fraction", "%", (/ ('', i=1, 3) /))

! Scops and others var (1)
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

  SUBROUTINE lmdz_cosp_output_open(Nlevlmdz, Ncolumns, presnivs, dtime, freq_cosp, &
                              ok_mensuelCOSP, ok_journeCOSP, ok_hfCOSP, ok_all_xml,  &
                              ecrit_mth, ecrit_day, ecrit_hf, use_vgrid, Nlvgrid, vgrid_z_loc, &
			      vgrid_mz)

  use mod_cosp_config, only : CLOUDSAT_DBZE_BINS, SR_BINS, CLOUDSAT_CFAD_ZE_MIN, PARASOL_NREFL, &
                              CLOUDSAT_CFAD_ZE_WIDTH,PARASOL_SZA, &
                              isccp_histPresCenters,tau_binCenters, LIDAR_NTEMP, &
                              LIDAR_PHASE_TEMP,misr_histHgtCenters,numMISRHgtBins, &
                              numMODISReffIceBins,reffICE_binCenters, &
                              numMODISReffLiqBins, reffLIQ_binCenters, pres_binCenters, &
			      cloudsat_binCenters, calipso_binCenters

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
  integer                  :: Nlevlmdz, Ncolumns, Nlvgrid      ! Number of levels
  real,dimension(Nlevlmdz) :: presnivs, vgrid_mz
  real,dimension(Nlvgrid)  :: vgrid_z_loc
  real                     :: dtime, freq_cosp, ecrit_day, ecrit_hf, ecrit_mth 
  logical                  :: use_vgrid
  logical                  :: ok_mensuelCOSP, ok_journeCOSP, ok_hfCOSP, ok_all_xml

!!! Variables locales
  integer                   :: idayref, iff, ii
  real                      :: zjulian,zjulian_start
  real(wp),dimension(Ncolumns)  :: column_ax
  CHARACTER(LEN=20), DIMENSION(3)  :: chfreq = (/ '1day', '1d  ', '3h  ' /)            

!!! Variables d'entree

#ifdef CPP_XIOS
    ! ug Variables utilisées pour récupérer le calendrier pour xios
    INTEGER :: x_an, x_mois, x_jour
    REAL :: x_heure
    INTEGER :: ini_an, ini_mois, ini_jour
    REAL :: ini_heure
#endif

    WRITE(lunout,*) 'Debut lmdz_cosp_output_mod.F90'
    print*,'cosp_varsdefined',cosp_varsdefined
    ! Initialisations (Valeurs par defaut)


!! Definition valeurs axes
    do ii=1,Ncolumns
      column_ax(ii) = real(ii)
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

    CALL wxios_add_vaxis("height", Nlvgrid, vgrid_z_loc)
    print*,'wxios_add_vaxis Nlvgrid, vgrid_z_loc',Nlvgrid,vgrid_z_loc

    WRITE(lunout,*) 'wxios_add_vaxis height_mlev, Nlevlmdz vgrid_mz ', &
                     Nlevlmdz,vgrid_mz
    CALL wxios_add_vaxis("height_mlev", Nlevlmdz, vgrid_mz)

    WRITE(lunout,*) 'wxios_add_vaxis sza, PARASOL_NREFL ', &
                     PARASOL_NREFL, PARASOL_SZA
    CALL wxios_add_vaxis("sza", PARASOL_NREFL, PARASOL_SZA)

    WRITE(lunout,*) 'wxios_add_vaxis pressure2 ',7,pres_binCenters
    CALL wxios_add_vaxis("pressure2", 7, pres_binCenters)

    WRITE(lunout,*) 'wxios_add_vaxis column ',Ncolumns,column_ax
    CALL wxios_add_vaxis("column", Ncolumns, column_ax)

   WRITE(lunout,*) 'wxios_add_vaxis temp LIDAR_NTEMP, LIDAR_PHASE_TEMP ', &
                    LIDAR_NTEMP, LIDAR_PHASE_TEMP
   CALL wxios_add_vaxis("temp", LIDAR_NTEMP, LIDAR_PHASE_TEMP)

   WRITE(lunout,*) 'wxios_add_vaxis cth16 numMISRHgtBins, misr_histHgtCenters ', & 
                    numMISRHgtBins, misr_histHgtCenters 
   CALL wxios_add_vaxis("cth", numMISRHgtBins, misr_histHgtCenters)

   WRITE(lunout,*) 'wxios_add_vaxis dbze CLOUDSAT_DBZE_BINS, cloudsat_binCenters ', &
                    CLOUDSAT_DBZE_BINS, cloudsat_binCenters
   CALL wxios_add_vaxis("dbze", CLOUDSAT_DBZE_BINS, cloudsat_binCenters) 

   WRITE(lunout,*) 'wxios_add_vaxis scatratio SR_BINS, calipso_binCenters', &
                   SR_BINS, calipso_binCenters
   CALL wxios_add_vaxis("scatratio", SR_BINS, calipso_binCenters)

   WRITE(lunout,*) 'wxios_add_vaxis ReffIce numMODISReffIceBins, &
                   reffICE_binCenters',numMODISReffIceBins, reffICE_binCenters
   CALL wxios_add_vaxis("ReffIce", numMODISReffIceBins, reffICE_binCenters)

   WRITE(lunout,*) 'wxios_add_vaxis ReffLiq numMODISReffLiqBins, &
                    reffLIQ_binCenters', numMODISReffLiqBins, reffLIQ_binCenters
   CALL wxios_add_vaxis("ReffLiq", numMODISReffLiqBins, reffLIQ_binCenters)

   WRITE(lunout,*) 'wxios_add_vaxis 7, tau_binCenters', &
                    7, tau_binCenters
   CALL wxios_add_vaxis("tau", 7, tau_binCenters)

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
      CALL histvert(cosp_nidfiles(iff),"height","height","m",Nlvgrid,vgrid_z_loc,nvert(iff))
       else
! Axe vertical modele LMDZ presnivs
      CALL histvert(cosp_nidfiles(iff),"presnivs","Vertical levels","Pa",Nlevlmdz,presnivs,nvert(iff),"down")
       endif
! Axe vertical niveaux modele (en m)
      CALL histvert(cosp_nidfiles(iff),"height_mlev","height_mlev","m",Nlevlmdz,vgrid_mz,nvertmcosp(iff))

      CALL histvert(cosp_nidfiles(iff),"sza","solar_zenith_angle","degrees",PARASOL_NREFL,PARASOL_SZA,nvertp(iff))

      CALL histvert(cosp_nidfiles(iff),"pressure2","pressure","mb",7,pres_binCenters,nvertisccp(iff),"down")

      CALL histvert(cosp_nidfiles(iff),"column","column","count",Ncolumns,column_ax,nvertcol(iff)) !DBUG

      CALL histvert(cosp_nidfiles(iff),"temp","temperature","C",LIDAR_NTEMP,LIDAR_PHASE_TEMP,nverttemp(iff))

      CALL histvert(cosp_nidfiles(iff),"cth","altitude","m",numMISRHgtBins,misr_histHgtCenters,nvertmisr(iff))
  
      CALL histvert(cosp_nidfiles(iff),"ReffIce","Effective_particle_size_Ice","microns",numMODISReffIceBins, reffICE_binCenters, &
                    nvertReffIce(iff))                                         
     
      CALL histvert(cosp_nidfiles(iff),"ReffLiq","Effective_particle_size_Liq","microns",numMODISReffLiqBins, reffLIQ_binCenters, &                                  
                    nvertReffLiq(iff))

      CALL histvert(cosp_nidfiles(iff),"dbze","equivalent_reflectivity_factor","dBZ",CLOUDSAT_DBZE_BINS,cloudsat_binCenters,nvertbze(iff))
     
      CALL histvert(cosp_nidfiles(iff),"scatratio","backscattering_ratio","1",SR_BINS,calipso_binCenters,nvertsratio(iff))

      CALL histvert(cosp_nidfiles(iff),"tau","cloud optical depth","1",7,tau_binCenters,nverttau(iff)) 
     
!!! Valeur indefinie en cas IOIPSL
     Cosp_fill_value=0.

#endif

      ENDIF
  ENDDO

    end SUBROUTINE lmdz_cosp_output_open

 END MODULE lmdz_cosp_output_mod
