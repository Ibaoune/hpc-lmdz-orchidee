!!!=============================================================================
!!! AI mars 2018
!!  Module permettant de controler les cles de sorties cosp
!!  pour LMDZ
!! 1. on initialise les cles au 1er passage a cosp itap de la physique = 1
!! 2. on garde la routine de lecture du fichier namelist cosp_out...txt pour le
!!    cas non XIOS (ioipsl)
!! 3. on rajoutte une subroutine qui interoge XIOS si les champs sont demandes
!!    dans les xml alors on les active et on active les simulateurs
!!    correspondant 
!!
!!! RG janvier 2019
!!  Mise a jour pour avoir toutes les sorties de COSPv2
!!  Changement de nom a "lmdz_cosp_read_outputkeys"
!!!=============================================================================

module lmdz_cosp_read_outputkeys

  USE mod_phys_lmdz_para

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE COSP_CONFIG
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE COSP_CONFIG
     logical :: &
         Lwrite_output,    & ! Control for output file
         Lfracout,         & ! Control SCOPS cloud fraction out ?
! simulators
         Lcloudsat,        & ! CLOUDSAT simulator on/off switch 
         Lcalipso,         & ! CALIPSO simulator on/off switch 
	 LgrLidar532,      & ! GROUND LIDAR simulator on/off switch
         Latlid,           & ! ATLID simulator on/off switch
         Lisccp,           & ! ISCCP simulator on/off switch
         Lmodis,           & ! MODIS simulatoe on/off switch
         Lmisr,            & ! MISR simulator on/off switch 
         Lrttov,           & ! RTTOV simulator on/off switch 
         Lparasol,         & ! PARASOL simulator on/off switch

! diagnostics 
         Lpctisccp,        & ! ISCCP mean cloud top pressure
         Lclisccp,         & ! ISCCP cloud area fraction
         Lboxptopisccp,    & ! ISCCP CTP in each column
         Lboxtauisccp,     & ! ISCCP optical epth in each column
         Ltauisccp,        & ! ISCCP mean optical depth
         Lcltisccp,        & ! ISCCP total cloud fraction
         Lmeantbisccp,     & ! ISCCP mean all-sky 10.5micron brightness temperature
         Lmeantbclrisccp,  & ! ISCCP mean clear-sky 10.5micron brightness temperature
         Lalbisccp,        & ! ISCCP mean cloud albedo         
         LclMISR,          & ! MISR cloud fraction
         Lcltmodis,        & ! MODIS total cloud fraction
         Lclwmodis,        & ! MODIS liquid cloud fraction
         Lclimodis,        & ! MODIS ice cloud fraction
         Lclhmodis,        & ! MODIS high-level cloud fraction
         Lclmmodis,        & ! MODIS mid-level cloud fraction
         Lcllmodis,        & ! MODIS low-level cloud fraction
         Ltautmodis,       & ! MODIS total cloud optical thicknes
         Ltauwmodis,       & ! MODIS liquid optical thickness
         Ltauimodis,       & ! MODIS ice optical thickness
         Ltautlogmodis,    & ! MODIS total cloud optical thickness (log10 mean)
         Ltauwlogmodis,    & ! MODIS liquid optical thickness (log10 mean)
         Ltauilogmodis,    & ! MODIS ice optical thickness (log10 mean)
         Lreffclwmodis,    & ! MODIS liquid cloud particle size
         Lreffclimodis,    & ! MODIS ice particle size
         Lpctmodis,        & ! MODIS cloud top pressure
         Llwpmodis,        & ! MODIS cloud liquid water path
         Liwpmodis,        & ! MODIS cloud ice water path
         Lclmodis,         & ! MODIS cloud area fraction
         Latb532,          & ! CALIPSO attenuated total backscatter (532nm)
         Latb532gr,        & ! GROUND LIDAR @ 532NM attenuated total backscatter (532nm)
         Latb355,          & ! ATLID attenuated total backscatter (355nm)
         LlidarBetaMol532, & ! CALIPSO molecular backscatter (532nm)         
         LlidarBetaMol532gr,&! GROUND LIDAR @ 532NM molecular backscatter (532nm)
         LlidarBetaMol355, & ! ATLID molecular backscatter (355nm) 
         LcfadLidarsr532,  & ! CALIPSO scattering ratio CFAD
         LcfadLidarsr532gr,& ! GROUND LIDAR @ 532NM scattering ratio CFAD  
         LcfadLidarsr355,  & ! ATLID scattering ratio CFAD 
         Lclcalipso2,      & ! CALIPSO cloud fraction undetected by cloudsat
         Lclcalipso,       & ! CALIPSO cloud area fraction
         LclgrLidar532,   & ! GROUND LIDAR @ 532NM cloud area fraction 
         Lclatlid,         & ! ATLID cloud area fraction 
         Lclhcalipso,      & ! CALIPSO high-level cloud fraction
         Lcllcalipso,      & ! CALIPSO low-level cloud fraction
         Lclmcalipso,      & ! CALIPSO mid-level cloud fraction
         Lcltcalipso,      & ! CALIPSO total cloud fraction
         LclhgrLidar532,  & ! GROUND LIDAR @ 532NM high-level cloud fraction 
         LcllgrLidar532,  & ! GROUND LIDAR @ 532NM low-level cloud fraction 
         LclmgrLidar532,  & ! GROUND LIDAR @ 532NM mid-level cloud fraction
         LcltgrLidar532,  & ! GROUND LIDAR @ 532NM total cloud fraction
         Lclhatlid,        & ! ATLID high-level cloud fraction
         Lcllatlid,        & ! ATLID low-level cloud fraction  
         Lclmatlid,        & ! ATLID mid-level cloud fraction 
         Lcltatlid,        & ! ATLID total cloud fraction
         Lcltlidarradar,   & ! CALIPSO-CLOUDSAT total cloud fraction
         Lcloudsat_tcc,    & !
         Lcloudsat_tcc2,   & !
         Lclcalipsoliq,    & ! CALIPSO liquid cloud area fraction
         Lclcalipsoice,    & ! CALIPSO ice cloud area fraction 
         Lclcalipsoun,     & ! CALIPSO undetected cloud area fraction
         Lclcalipsotmp,    & ! CALIPSO undetected cloud area fraction
         Lclcalipsotmpliq, & ! CALIPSO liquid cloud area fraction
         Lclcalipsotmpice, & ! CALIPSO ice cloud area fraction
         Lclcalipsotmpun,  & ! CALIPSO undetected cloud area fraction
         Lcltcalipsoliq,   & ! CALIPSO liquid total cloud fraction
         Lcltcalipsoice,   & ! CALIPSO ice total cloud fraction
         Lcltcalipsoun,    & ! CALIPSO undetected total cloud fraction
         Lclhcalipsoliq,   & ! CALIPSO high-level liquid cloud fraction
         Lclhcalipsoice,   & ! CALIPSO high-level ice cloud fraction
         Lclhcalipsoun,    & ! CALIPSO high-level undetected cloud fraction
         Lclmcalipsoliq,   & ! CALIPSO mid-level liquid cloud fraction
         Lclmcalipsoice,   & ! CALIPSO mid-level ice cloud fraction
         Lclmcalipsoun,    & ! CALIPSO mid-level undetected cloud fraction
         Lcllcalipsoliq,   & ! CALIPSO low-level liquid cloud fraction
         Lcllcalipsoice,   & ! CALIPSO low-level ice cloud fraction
         Lcllcalipsoun,    & ! CALIPSO low-level undetected cloud fraction
         Lclopaquecalipso, & ! CALIPSO opaque cloud cover (2D Map)
         Lclthincalipso,   & ! CALIPSO thin cloud cover (2D Map)
         Lclzopaquecalipso,& ! CALIPSO z_opaque altitude (opaque clouds only, 2D Map)
         Lclcalipsoopaque, & ! CALIPSO opaque cloud profiles 3D fraction 
         Lclcalipsothin,   & ! CALIPSO thin cloud profiles 3D fraction 
         Lclcalipsozopaque,& ! CALIPSO z_opaque 3D fraction 
         Lclcalipsoopacity,& ! CALIPSO opacity 3D fraction 
         Lclopaquetemp,    & ! CALIPSO opaque cloud temperature 
         Lclthintemp,      & ! CALIPSO thin cloud temperature
         Lclzopaquetemp,   & ! CALIPSO z_opaque temperature  
         Lclopaquemeanz,   & ! CALIPSO opaque cloud altitude  
         Lclthinmeanz,     & ! CALIPSO thin cloud altitude 
         Lclthinemis,      & ! CALIPSO thin cloud emissivity
         Lclopaquemeanzse,   & ! CALIPSO opaque cloud altitude with respect to SE 
         Lclthinmeanzse,     & ! CALIPSO thin cloud altitude with respect to SE
         Lclzopaquecalipsose,& ! CALIPSO z_opaque altitude with respect to SE
         LcfadDbze94,      & ! CLOUDSAT radar reflectivity CFAD
         Ldbze94,          & ! CLOUDSAT radar reflectivity
         LparasolRefl,     & ! PARASOL reflectance
         Ltbrttov,         & ! RTTOV mean clear-sky brightness temperature
         Lptradarflag0,    & ! CLOUDSAT no precipitation flag
         Lptradarflag1,    & ! CLOUDSAT rain possible flag
         Lptradarflag2,    & ! CLOUDSAT rain probable flag
         Lptradarflag3,    & ! CLOUDSAT rain certain flag
         Lptradarflag4,    & ! CLOUDSAT snow possible flag
         Lptradarflag5,    & ! CLOUDSAT snow certain flag
         Lptradarflag6,    & ! CLOUDSAT mixed-precipitation possible flag
         Lptradarflag7,    & ! CLOUDSAT mixed-precipitation certain flag
         Lptradarflag8,    & ! CLOUDSAT heavy rain flag
         Lptradarflag9,    & ! CLOUDSAT default flag
         Lradarpia           ! CLOUDSAT PIA


!     character(len=32),dimension(:),allocatable :: out_list
      character(len=32) :: out_list(107)
  END TYPE COSP_CONFIG       


CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE READ_COSP_OUTPUT_NL -------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SUBROUTINE cosp_outputkeys_init(cfg)
  implicit none
  type(cosp_config),intent(inout) :: cfg 
  character(len=32) :: out_list(107)
  integer :: i
                
   do i=1,107
      cfg%out_list(i)=''
   enddo


   cfg%Lcloudsat=.false.
   cfg%Lcalipso=.false.
   cfg%LgrLidar532=.false.
   cfg%Latlid=.false.
   cfg%Lisccp=.false.
   cfg%Lmodis=.false.
   cfg%Lmisr=.false.
   cfg%Lrttov=.false.
   cfg%Lparasol=.false.
   cfg%Lwrite_output=.false.
   cfg%Lfracout=.false.


   cfg%Lpctisccp=.false.
   cfg%Lclisccp=.false.
   cfg%Lboxptopisccp=.false.
   cfg%Lboxtauisccp=.false.
   cfg%Ltauisccp=.false.
   cfg%Lcltisccp=.false.
   cfg%Lmeantbisccp=.false.
   cfg%Lmeantbclrisccp=.false.
   cfg%Lalbisccp=.false.
   cfg%LclMISR=.false.
   cfg%Lcltmodis=.false.
   cfg%Lclwmodis=.false.
   cfg%Lclimodis=.false.
   cfg%Lclhmodis=.false.
   cfg%Lclmmodis=.false.
   cfg%Lcllmodis=.false.
   cfg%Ltautmodis=.false.
   cfg%Ltauwmodis=.false.
   cfg%Ltauimodis=.false.
   cfg%Ltautlogmodis=.false.
   cfg%Ltauwlogmodis=.false.
   cfg%Ltauilogmodis=.false.
   cfg%Lreffclwmodis=.false.
   cfg%Lreffclimodis=.false.
   cfg%Lpctmodis=.false.
   cfg%Llwpmodis=.false.
   cfg%Liwpmodis=.false.
   cfg%Lclmodis=.false.
   cfg%Latb532=.false.
   cfg%Latb532gr=.false.
   cfg%Latb355=.false.
   cfg%LlidarBetaMol532=.false.
   cfg%LlidarBetaMol532gr=.false.
   cfg%LlidarBetaMol355=.false.
   cfg%LcfadLidarsr532=.false.
   cfg%LcfadLidarsr532gr=.false.
   cfg%LcfadLidarsr355 =.false.
   cfg%Lclcalipso2=.false.
   cfg%Lclcalipso=.false.
   cfg%LclgrLidar532=.false.
   cfg%Lclatlid=.false.
   cfg%Lclhcalipso=.false.
   cfg%Lcllcalipso=.false.
   cfg%Lclmcalipso=.false.
   cfg%Lcltcalipso=.false.
   cfg%LclhgrLidar532=.false.
   cfg%LcllgrLidar532=.false.
   cfg%LclmgrLidar532 =.false.
   cfg%LcltgrLidar532=.false.
   cfg%Lclhatlid=.false.
   cfg%Lcllatlid=.false.
   cfg%Lclmatlid=.false.
   cfg%Lcltatlid=.false.
   cfg%Lcltlidarradar=.false.
   cfg%Lcloudsat_tcc=.false.
   cfg%Lcloudsat_tcc2=.false.
   cfg%Lclcalipsoliq=.false.
   cfg%Lclcalipsoice=.false.
   cfg%Lclcalipsoun=.false.
   cfg%Lclcalipsotmp=.false.
   cfg%Lclcalipsotmpliq=.false.
   cfg%Lclcalipsotmpice=.false.
   cfg%Lclcalipsotmpun=.false.
   cfg%Lcltcalipsoliq=.false.
   cfg%Lcltcalipsoice=.false.
   cfg%Lcltcalipsoun=.false.
   cfg%Lclhcalipsoliq=.false.
   cfg%Lclhcalipsoice=.false.
   cfg%Lclhcalipsoun=.false.
   cfg%Lclmcalipsoliq=.false.
   cfg%Lclmcalipsoice=.false.
   cfg%Lclmcalipsoun=.false.
   cfg%Lcllcalipsoliq=.false.
   cfg%Lcllcalipsoice=.false.
   cfg%Lcllcalipsoun =.false.
   cfg%Lclopaquecalipso=.false.
   cfg%Lclthincalipso=.false.
   cfg%Lclzopaquecalipso=.false.
   cfg%Lclcalipsoopaque=.false.
   cfg%Lclcalipsothin=.false.
   cfg%Lclcalipsozopaque=.false.
   cfg%Lclcalipsoopacity=.false.
   cfg%Lclopaquetemp=.false.
   cfg%Lclthintemp=.false.
   cfg%Lclzopaquetemp=.false.
   cfg%Lclopaquemeanz=.false.
   cfg%Lclthinmeanz=.false.
   cfg%Lclthinemis=.false.
   cfg%Lclopaquemeanzse=.false.
   cfg%Lclthinmeanzse=.false.
   cfg%Lclzopaquecalipsose=.false.
   cfg%LcfadDbze94=.false.
   cfg%Ldbze94=.false.
   cfg%Lparasolrefl=.false.
   cfg%Ltbrttov=.false.
   cfg%Lptradarflag0=.false.
   cfg%Lptradarflag1=.false.
   cfg%Lptradarflag2=.false.
   cfg%Lptradarflag3=.false.
   cfg%Lptradarflag4=.false.
   cfg%Lptradarflag5=.false.
   cfg%Lptradarflag6=.false.
   cfg%Lptradarflag7=.false.
   cfg%Lptradarflag8=.false.
   cfg%Lptradarflag9=.false.
   cfg%Lradarpia=.false.

 end subroutine cosp_outputkeys_init


 SUBROUTINE cosp_outputkeys_test(cfg)
  implicit none
  type(cosp_config),intent(out) :: cfg
  character(len=32) :: out_list(107)
  integer :: i

   do i=1,107
      cfg%out_list(i)=''
   enddo

   cfg%Lcloudsat=.false.
   cfg%Lcalipso=.true.
   cfg%LgrLidar532=.false.
   cfg%Latlid=.false.
   cfg%Lisccp=.false.
   cfg%Lmodis=.false.
   cfg%Lmisr=.false.
   cfg%Lrttov=.false.
   cfg%Lparasol=.false.
   cfg%Lwrite_output=.false.
!   cfg%Lfracout=.false.

   cfg%Lpctisccp=.false.
   cfg%Lclisccp=.false.
   cfg%Lboxptopisccp=.false.
   cfg%Lboxtauisccp=.false.
   cfg%Ltauisccp=.false.
   cfg%Lcltisccp=.false.
   cfg%Lmeantbisccp=.false.
   cfg%Lmeantbclrisccp=.false.
   cfg%Lalbisccp=.false.
   cfg%LclMISR=.false.
   cfg%Lcltmodis=.false.
   cfg%Lclwmodis=.false.
   cfg%Lclimodis=.false.
   cfg%Lclhmodis=.false.
   cfg%Lclmmodis=.false.
   cfg%Lcllmodis=.false.
   cfg%Ltautmodis=.false.
   cfg%Ltauwmodis=.false.
   cfg%Ltauimodis=.false.
   cfg%Ltautlogmodis=.false.
   cfg%Ltauwlogmodis=.false.
   cfg%Ltauilogmodis=.false.
   cfg%Lreffclwmodis=.false.
   cfg%Lreffclimodis=.false.
   cfg%Lpctmodis=.false.
   cfg%Llwpmodis=.false.
   cfg%Liwpmodis=.false.
   cfg%Lclmodis=.false.
   cfg%Latb532=.false.
   cfg%Latb532gr=.false.
   cfg%Latb355=.false.
   cfg%LlidarBetaMol532=.false.
   cfg%LlidarBetaMol532gr=.false.
   cfg%LlidarBetaMol355=.false.
   cfg%LcfadLidarsr532=.false.
   cfg%LcfadLidarsr532gr=.false.
   cfg%LcfadLidarsr355 =.false.
   cfg%Lclcalipso2=.false.
   cfg%Lclcalipso=.false.
   cfg%LclgrLidar532=.false.
   cfg%Lclatlid=.false.
   cfg%Lclhcalipso=.true.
   cfg%Lcllcalipso=.true.
   cfg%Lclmcalipso=.true.
   cfg%Lcltcalipso=.true.
   cfg%LclhgrLidar532=.false.
   cfg%LcllgrLidar532=.false.
   cfg%LclmgrLidar532 =.false.
   cfg%LcltgrLidar532=.false.
   cfg%Lclhatlid=.false.
   cfg%Lcllatlid=.false.
   cfg%Lclmatlid=.false.
   cfg%Lcltatlid=.false.
   cfg%Lcltlidarradar=.false.
   cfg%Lcloudsat_tcc=.false.
   cfg%Lcloudsat_tcc2=.false.
   cfg%Lclcalipsoliq=.false.
   cfg%Lclcalipsoice=.false.
   cfg%Lclcalipsoun=.false.
   cfg%Lclcalipsotmp=.false.
   cfg%Lclcalipsotmpliq=.false.
   cfg%Lclcalipsotmpice=.false.
   cfg%Lclcalipsotmpun=.false.
   cfg%Lcltcalipsoliq=.false.
   cfg%Lcltcalipsoice=.false.
   cfg%Lcltcalipsoun=.false.
   cfg%Lclhcalipsoliq=.false.
   cfg%Lclhcalipsoice=.false.
   cfg%Lclhcalipsoun=.false.
   cfg%Lclmcalipsoliq=.false.
   cfg%Lclmcalipsoice=.false.
   cfg%Lclmcalipsoun=.false.
   cfg%Lcllcalipsoliq=.false.
   cfg%Lcllcalipsoice=.false.
   cfg%Lcllcalipsoun =.false.
   cfg%Lclopaquecalipso=.false.
   cfg%Lclthincalipso=.false.
   cfg%Lclzopaquecalipso=.false.
   cfg%Lclcalipsoopaque=.false.
   cfg%Lclcalipsothin=.false.
   cfg%Lclcalipsozopaque=.false.
   cfg%Lclcalipsoopacity=.false.
   cfg%Lclopaquetemp=.false.
   cfg%Lclthintemp=.false.
   cfg%Lclzopaquetemp=.false.
   cfg%Lclopaquemeanz=.false.
   cfg%Lclthinmeanz=.false.
   cfg%Lclthinemis=.false.
   cfg%Lclopaquemeanzse=.false.
   cfg%Lclthinmeanzse=.false.
   cfg%Lclzopaquecalipsose=.false.
   cfg%LcfadDbze94=.false.
   cfg%Ldbze94=.false.
   cfg%Lparasolrefl=.false.
   cfg%Ltbrttov=.false.
   cfg%Lptradarflag0=.false.
   cfg%Lptradarflag1=.false.
   cfg%Lptradarflag2=.false.
   cfg%Lptradarflag3=.false.
   cfg%Lptradarflag4=.false.
   cfg%Lptradarflag5=.false.
   cfg%Lptradarflag6=.false.
   cfg%Lptradarflag7=.false.
   cfg%Lptradarflag8=.false.
   cfg%Lptradarflag9=.false.
   cfg%Lradarpia=.false.

 end subroutine cosp_outputkeys_test


 SUBROUTINE READ_COSP_OUTPUT_NL(itap,cosp_nl,cfg)

#ifdef CPP_XIOS
    USE xios, ONLY: xios_field_is_active
#endif
  implicit none
  character(len=*),intent(in) :: cosp_nl
  type(cosp_config),intent(inout) :: cfg
  ! Local variables
  integer :: i, itap

!COSPv2
  logical, save :: &
       Lisccp      = .false. ,& ! Local on/off switch for simulators
       Lmodis      = .false., & !
       Lmisr       = .false., & !
       Lcalipso    = .false., & !
       LgrLidar532 = .false., & !
       Latlid      = .false., & !
       Lcloudsat   = .false., & !
       Lrttov      = .false., & !
       Lparasol    = .false.    !


 logical, save :: Lcfaddbze94,Ldbze94,Latb532,LcfadLidarsr532,Lclcalipso,Lclhcalipso,    &
             Lcllcalipso,Lclmcalipso,Lcltcalipso,LparasolRefl,Lclcalipsoliq,             &
             Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice, &
             Lclcalipsotmpun,Lclhcalipsoliq,Lcllcalipsoliq,Lclmcalipsoliq,Lcltcalipsoliq,& 
             Lclhcalipsoice,Lcllcalipsoice,Lclmcalipsoice,Lcltcalipsoice,Lclhcalipsoun,  & 
             Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lclopaquecalipso,Lclthincalipso,  & 
             Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,        & 
             Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,Lclzopaquetemp,Lclopaquemeanz,  &
             Lclthinmeanz,Lclthinemis,Lclopaquemeanzse,Lclthinmeanzse,Lclzopaquecalipsose,&
             LlidarBetaMol532gr,LcfadLidarsr532gr,Latb532gr,LclgrLidar532,LclhgrLidar532,&
             LcllgrLidar532,LclmgrLidar532,LcltgrLidar532,LlidarBetaMol355,              &
             LcfadLidarsr355,Latb355,Lclatlid,Lclhatlid,Lcllatlid,Lclmatlid,Lcltatlid,   &
             Lalbisccp,Lboxptopisccp,Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,&
             Lmeantbisccp,Lmeantbclrisccp,LclMISR,Lclcalipso2,Lcltlidarradar,            &
             Lcloudsat_tcc, Lcloudsat_tcc2,Lfracout,   &
             LlidarBetaMol532,Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,         &
             Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis,Ltauwlogmodis,     &
             Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,    &
             Lclmodis,Ltbrttov,Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,  &
             Lptradarflag4,Lptradarflag5,Lptradarflag6,Lptradarflag7,Lptradarflag8,      &
             Lptradarflag9,Lradarpia

  namelist/COSP_OUTPUT/Lcfaddbze94,Ldbze94,Latb532,LcfadLidarsr532,Lclcalipso,           &
                       Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,LparasolRefl,     &
                       Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,           &
                       Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,Lclhcalipsoliq, &
                       Lcllcalipsoliq,Lclmcalipsoliq,Lcltcalipsoliq,Lclhcalipsoice,      &
                       Lcllcalipsoice,Lclmcalipsoice,Lcltcalipsoice,Lclhcalipsoun,       &
                       Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lclopaquecalipso,       &
                       Lclthincalipso,Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin, &
                       Lclcalipsozopaque,Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,    &
                       Lclzopaquetemp,Lclopaquemeanz,Lclthinmeanz,Lclthinemis,           &
                       Lclopaquemeanzse,Lclthinmeanzse,Lclzopaquecalipsose,              &
                       LlidarBetaMol532gr,LcfadLidarsr532gr,Latb532gr,LclgrLidar532,     &
                       LclhgrLidar532,LcllgrLidar532,LclmgrLidar532,LcltgrLidar532,      &
                       LlidarBetaMol355,LcfadLidarsr355,Latb355,Lclatlid,                &
                       Lclhatlid,Lcllatlid,Lclmatlid,Lcltatlid,Lalbisccp,Lboxptopisccp,  &
                       Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,Lmeantbisccp, &
                       Lmeantbclrisccp,LclMISR,Lclcalipso2,Lcltlidarradar,               &
                       Lcloudsat_tcc, Lcloudsat_tcc2, Lfracout,                          &
                       LlidarBetaMol532,Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,         &
                       Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,             &
                       Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,          &
                       Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,Ltbrttov,    &
                       Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,          &
                       Lptradarflag4,Lptradarflag5,Lptradarflag6,Lptradarflag7,          &
                       Lptradarflag8,Lptradarflag9,Lradarpia

  do i=1,107
    cfg%out_list(i)=''
  enddo
  
! Lecture du fichier namelist
  IF (is_master) THEN
    open(10,file=cosp_nl,status='old')
    read(10,nml=cosp_output)
    close(10)
  ENDIF

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Which simulators need to be run? Look at which outputs are requested.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (Lpctisccp .or. Lclisccp .or. Lboxptopisccp .or.  Lboxtauisccp .or. Ltauisccp .or. &
       Lcltisccp .or. Lmeantbisccp .or. Lmeantbclrisccp .or. Lalbisccp) Lisccp = .true.
  if (LclMISR) Lmisr = .true.
  if (Lcltmodis .or. Lclwmodis .or. Lclimodis .or. Lclhmodis .or. Lclmmodis .or.         &
       Lcllmodis .or. Ltautmodis .or. Ltauwmodis .or. Ltauimodis .or. Ltautlogmodis .or. &
       Ltauwlogmodis .or. Ltauilogmodis .or. Lreffclwmodis .or. Lreffclimodis .or.       &
       Lpctmodis .or. Llwpmodis .or. Liwpmodis .or. Lclmodis) Lmodis = .true.
  if (Lclcalipso2 .or. Lclcalipso .or.  Lclhcalipso .or. Lcllcalipso .or. Lclmcalipso    &
       .or. Lcltcalipso .or. Lcltlidarradar .or. Lclcalipsoliq .or. Lclcalipsoice .or.   &
       Lclcalipsoun .or. Lclcalipsotmp .or. Lclcalipsotmpliq .or. Lclcalipsotmpice .or.  &
       Lclcalipsotmpun .or. Lcltcalipsoliq .or. Lcltcalipsoice .or. Lcltcalipsoun .or.   &
       Lclhcalipsoliq .or. Lclhcalipsoice .or. Lclhcalipsoun .or. Lclmcalipsoliq .or.    &
       Lclmcalipsoice .or. Lclmcalipsoun .or. Lcllcalipsoliq .or. Lcllcalipsoice .or.    &
       Lcllcalipsoun .or. LlidarBetaMol532 .or. LcfadLidarsr532 .or. Lcltlidarradar .or. &
       Lcltlidarradar .or. Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso   &
       .or. Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or.             &
       Lclcalipsoopacity .or. Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp .or.    &
       Lclopaquemeanz .or. Lclthinmeanz .or. Lclthinemis .or. Lclopaquemeanzse .or.      & 
       Lclthinmeanzse .or. Lclzopaquecalipsose) Lcalipso = .true. 

  if (LlidarBetaMol532gr .or. LcfadLidarsr532gr .or. Latb532gr .or. LclgrLidar532 .or.  & 
       LclhgrLidar532 .or. LcllgrLidar532 .or. LclmgrLidar532 .or. LcltgrLidar532)   & 
       LgrLidar532 = .true.

  if (LlidarBetaMol355 .or. LcfadLidarsr355 .or. Latb355 .or. Lclatlid .or.              & 
       Lclhatlid .or. Lcllatlid .or. Lclmatlid .or. Lcltatlid)                           & 
       Latlid = .true. 

  if (LcfadDbze94 .or. Ldbze94 .or. Lcltlidarradar .or. Lptradarflag0 .or. Lptradarflag1 &
       .or. Lptradarflag2 .or. Lptradarflag3 .or. Lptradarflag4 .or. Lptradarflag5 .or.  &
       Lptradarflag6 .or. Lptradarflag7 .or. Lptradarflag8 .or. Lptradarflag9 .or.       &
       Lradarpia) Lcloudsat = .true.
  if (Lparasolrefl) Lparasol = .true.
  if (Ltbrttov) Lrttov = .true.

  CALL bcast(Lcloudsat)
  CALL bcast(Lcalipso)
  CALL bcast(LgrLidar532)
  CALL bcast(Latlid)
  CALL bcast(Lparasol)
  CALL bcast(Lisccp)
  CALL bcast(Lmodis)
  CALL bcast(Lmisr)
  CALL bcast(Lrttov)

  CALL bcast(Lcfaddbze94)
  CALL bcast(Ldbze94)
  CALL bcast(Latb532)
  CALL bcast(LcfadLidarsr532)
  CALL bcast(Lclcalipso)
  CALL bcast(Lclhcalipso)
  CALL bcast(Lcllcalipso)
  CALL bcast(Lclmcalipso)
  CALL bcast(Lcltcalipso)
  CALL bcast(LparasolRefl)
  CALL bcast(Lclcalipsoliq)
  CALL bcast(Lclcalipsoice)
  CALL bcast(Lclcalipsoun)
  CALL bcast(Lclcalipsotmp)
  CALL bcast(Lclcalipsotmpliq)
  CALL bcast(Lclcalipsotmpice)
  CALL bcast(Lclcalipsotmpun)
  CALL bcast(Lclhcalipsoliq)
  CALL bcast(Lcllcalipsoliq)
  CALL bcast(Lclmcalipsoliq)
  CALL bcast(Lcltcalipsoliq)
  CALL bcast(Lclhcalipsoice)
  CALL bcast(Lcllcalipsoice)
  CALL bcast(Lclmcalipsoice)
  CALL bcast(Lcltcalipsoice)
  CALL bcast(Lclhcalipsoun)
  CALL bcast(Lcllcalipsoun)
  CALL bcast(Lclmcalipsoun)
  CALL bcast(Lcltcalipsoun)
  CALL bcast(Lclopaquecalipso)
  CALL bcast(Lclthincalipso)
  CALL bcast(Lclzopaquecalipso)
  CALL bcast(Lclcalipsoopaque)
  CALL bcast(Lclcalipsothin)
  CALL bcast(Lclcalipsozopaque)
  CALL bcast(Lclcalipsoopacity)
  CALL bcast(Lclopaquetemp)
  CALL bcast(Lclthintemp)
  CALL bcast(Lclzopaquetemp)
  CALL bcast(Lclopaquemeanz)
  CALL bcast(Lclthinmeanz)
  CALL bcast(Lclthinemis)
  CALL bcast(Lclopaquemeanzse)
  CALL bcast(Lclthinmeanzse)
  CALL bcast(Lclzopaquecalipsose)
  CALL bcast(LlidarBetaMol532gr)
  CALL bcast(LcfadLidarsr532gr)
  CALL bcast(Latb532gr)
  CALL bcast(LclgrLidar532)
  CALL bcast(LclhgrLidar532)
  CALL bcast(LcllgrLidar532)
  CALL bcast(LclmgrLidar532)
  CALL bcast(LcltgrLidar532)
  CALL bcast(LlidarBetaMol355)
  CALL bcast(LcfadLidarsr355)
  CALL bcast(Latb355)
  CALL bcast(Lclatlid)
  CALL bcast(Lclhatlid)
  CALL bcast(Lcllatlid)
  CALL bcast(Lclmatlid)
  CALL bcast(Lcltatlid)
  CALL bcast(Lalbisccp)
  CALL bcast(Lboxptopisccp)
  CALL bcast(Lboxtauisccp)
  CALL bcast(Lpctisccp)
  CALL bcast(Lclisccp)
  CALL bcast(Ltauisccp)
  CALL bcast(Lcltisccp)
  CALL bcast(Lmeantbisccp)
  CALL bcast(Lmeantbclrisccp)
  CALL bcast(LclMISR)
  CALL bcast(Lclcalipso2)
  CALL bcast(Lcltlidarradar)
  CALL bcast(Lcloudsat_tcc)
  CALL bcast(Lcloudsat_tcc2)
  CALL bcast(Lfracout)
  CALL bcast(LlidarBetaMol532)
  CALL bcast(Lcltmodis)
  CALL bcast(Lclwmodis)
  CALL bcast(Lclimodis)
  CALL bcast(Lclhmodis)
  CALL bcast(Lclmmodis)
  CALL bcast(Lcllmodis)
  CALL bcast(Ltautmodis)
  CALL bcast(Ltauwmodis)
  CALL bcast(Ltauimodis)
  CALL bcast(Ltautlogmodis)
  CALL bcast(Ltauwlogmodis)
  CALL bcast(Ltauilogmodis)
  CALL bcast(Lreffclwmodis)
  CALL bcast(Lreffclimodis)
  CALL bcast(Lpctmodis)
  CALL bcast(Llwpmodis)
  CALL bcast(Liwpmodis)
  CALL bcast(Lclmodis)
  CALL bcast(Ltbrttov)
  CALL bcast(Lptradarflag0)
  CALL bcast(Lptradarflag1)
  CALL bcast(Lptradarflag2)
  CALL bcast(Lptradarflag3)
  CALL bcast(Lptradarflag4)
  CALL bcast(Lptradarflag5)
  CALL bcast(Lptradarflag6)
  CALL bcast(Lptradarflag7)
  CALL bcast(Lptradarflag8)
  CALL bcast(Lptradarflag9)
  CALL bcast(Lradarpia)

!$OMP BARRIER

  ! Copy instrument flags to cfg structure
  cfg%Lcloudsat = Lcloudsat
  cfg%Lcalipso = Lcalipso
  cfg%LgrLidar532 = LgrLidar532
  cfg%Latlid = Latlid
  cfg%Lisccp = Lisccp
  cfg%Lmodis = Lmodis
  cfg%Lmisr = Lmisr
  cfg%Lrttov = Lrttov
  cfg%Lparasol = Lparasol

    print*,' cfg%Lcloudsat = ',cfg%Lcloudsat

  ! Flag to control output to file
  cfg%Lwrite_output = .false.
  if (cfg%Lcloudsat .or. cfg%Lcalipso .or. cfg%LgrLidar532 .or. cfg%Latlid .or. cfg%Lisccp .or. cfg%Lmodis .or. cfg%Lmisr .or. cfg%Lrttov .or. cfg%Lparasol .or. Lfracout) then
    cfg%Lwrite_output = .true.
  endif

  ! Output diagnostics
  i = 1

  if (Lcfaddbze94)        cfg%out_list(i) = 'cfaddbze94'
  i = i+1
  if (Ldbze94)        cfg%out_list(i) = 'dbze94'
  i = i+1
  if (Latb532)        cfg%out_list(i) = 'atb532'
  i = i+1
  if (LcfadLidarsr532)        cfg%out_list(i) = 'cfadLidarsr532'
  i = i+1
  if (Lclcalipso)        cfg%out_list(i) = 'clcalipso'
  i = i+1
  if (Lclhcalipso)        cfg%out_list(i) = 'clhcalipso'
  i = i+1
  if (Lcllcalipso)        cfg%out_list(i) = 'cllcalipso'
  i = i+1
  if (Lclmcalipso)        cfg%out_list(i) = 'clmcalipso'
  i = i+1
  if (Lcltcalipso)        cfg%out_list(i) = 'cltcalipso'
  i = i+1
  if (LparasolRefl)        cfg%out_list(i) = 'parasolRefl'
  i = i+1
  if (Lclcalipsoliq)        cfg%out_list(i) = 'clcalipsoliq'
  i = i+1
  if (Lclcalipsoice)        cfg%out_list(i) = 'clcalipsoice'
  i = i+1
  if (Lclcalipsoun)        cfg%out_list(i) = 'clcalipsoun'
  i = i+1
  if (Lclcalipsotmp)        cfg%out_list(i) = 'clcalipsotmp'
  i = i+1
  if (Lclcalipsotmpliq)        cfg%out_list(i) = 'clcalipsotmpliq'
  i = i+1
  if (Lclcalipsotmpice)        cfg%out_list(i) = 'clcalipsotmpice'
  i = i+1
  if (Lclcalipsotmpun)        cfg%out_list(i) = 'clcalipsotmpun'
  i = i+1
  if (Lclhcalipsoliq)        cfg%out_list(i) = 'clhcalipsoliq'
  i = i+1
  if (Lcllcalipsoliq)        cfg%out_list(i) = 'cllcalipsoliq'
  i = i+1
  if (Lclmcalipsoliq)        cfg%out_list(i) = 'clmcalipsoliq'
  i = i+1
  if (Lcltcalipsoliq)        cfg%out_list(i) = 'cltcalipsoliq'
  i = i+1
  if (Lclhcalipsoice)        cfg%out_list(i) = 'clhcalipsoice'
  i = i+1
  if (Lcllcalipsoice)        cfg%out_list(i) = 'cllcalipsoice'
  i = i+1
  if (Lclmcalipsoice)        cfg%out_list(i) = 'clmcalipsoice'
  i = i+1
  if (Lcltcalipsoice)        cfg%out_list(i) = 'cltcalipsoice'
  i = i+1
  if (Lclhcalipsoun)        cfg%out_list(i) = 'clhcalipsoun'
  i = i+1
  if (Lcllcalipsoun)        cfg%out_list(i) = 'cllcalipsoun'
  i = i+1
  if (Lclmcalipsoun)        cfg%out_list(i) = 'clmcalipsoun'
  i = i+1
  if (Lcltcalipsoun)        cfg%out_list(i) = 'cltcalipsoun'
  i = i+1
  if (Lclopaquecalipso)        cfg%out_list(i) = 'clopaquecalipso'
  i = i+1
  if (Lclthincalipso)        cfg%out_list(i) = 'clthincalipso'
  i = i+1
  if (Lclzopaquecalipso)        cfg%out_list(i) = 'clzopaquecalipso'
  i = i+1
  if (Lclcalipsoopaque)        cfg%out_list(i) = 'clcalipsoopaque'
  i = i+1
  if (Lclcalipsothin)        cfg%out_list(i) = 'clcalipsothin'
  i = i+1
  if (Lclcalipsozopaque)        cfg%out_list(i) = 'clcalipsozopaque'
  i = i+1
  if (Lclcalipsoopacity)        cfg%out_list(i) = 'clcalipsoopacity'
  i = i+1
  if (Lclopaquetemp)        cfg%out_list(i) = 'clopaquetemp'
  i = i+1
  if (Lclthintemp)        cfg%out_list(i) = 'clthintemp'
  i = i+1
  if (Lclzopaquetemp)        cfg%out_list(i) = 'clzopaquetemp'
  i = i+1
  if (Lclopaquemeanz)        cfg%out_list(i) = 'clopaquemeanz'
  i = i+1
  if (Lclthinmeanz)        cfg%out_list(i) = 'clthinmeanz'
  i = i+1
  if (Lclthinemis)        cfg%out_list(i) = 'clthinemis'
  i = i+1
  if (Lclopaquemeanzse)        cfg%out_list(i) = 'clopaquemeanzse'
  i = i+1
  if (Lclthinmeanzse)        cfg%out_list(i) = 'clthinmeanzse'
  i = i+1
  if (Lclzopaquecalipsose)        cfg%out_list(i) = 'clzopaquecalipsose'
  i = i+1
  if (LlidarBetaMol532gr)        cfg%out_list(i) = 'lidarBetaMol532gr'
  i = i+1
  if (LcfadLidarsr532gr)        cfg%out_list(i) = 'cfadLidarsr532gr'
  i = i+1
  if (Latb532gr)        cfg%out_list(i) = 'atb532gr'
  i = i+1
  if (LclgrLidar532)        cfg%out_list(i) = 'clgrLidar532'
  i = i+1
  if (LclhgrLidar532)        cfg%out_list(i) = 'clhgrLidar532'
  i = i+1
  if (LcllgrLidar532)        cfg%out_list(i) = 'cllgrLidar532'
  i = i+1
  if (LclmgrLidar532)        cfg%out_list(i) = 'clmgrLidar532'
  i = i+1
  if (LcltgrLidar532)        cfg%out_list(i) = 'cltgrLidar532'
  i = i+1
  if (LlidarBetaMol355)        cfg%out_list(i) = 'lidarBetaMol355'
  i = i+1
  if (LcfadLidarsr355)        cfg%out_list(i) = 'cfadLidarsr355'
  i = i+1
  if (Latb355)        cfg%out_list(i) = 'atb355'
  i = i+1
  if (Lclatlid)        cfg%out_list(i) = 'clatlid'
  i = i+1
  if (Lclhatlid)        cfg%out_list(i) = 'clhatlid'
  i = i+1
  if (Lcllatlid)        cfg%out_list(i) = 'cllatlid'
  i = i+1
  if (Lclmatlid)        cfg%out_list(i) = 'clmatlid'
  i = i+1
  if (Lcltatlid)        cfg%out_list(i) = 'cltatlid'
  i = i+1
  if (Lalbisccp)        cfg%out_list(i) = 'albisccp'
  i = i+1
  if (Lboxptopisccp)        cfg%out_list(i) = 'boxptopisccp'
  i = i+1
  if (Lboxtauisccp)        cfg%out_list(i) = 'boxtauisccp'
  i = i+1
  if (Lpctisccp)        cfg%out_list(i) = 'pctisccp'
  i = i+1
  if (Lclisccp)        cfg%out_list(i) = 'clisccp'
  i = i+1
  if (Ltauisccp)        cfg%out_list(i) = 'tauisccp'
  i = i+1
  if (Lcltisccp)        cfg%out_list(i) = 'cltisccp'
  i = i+1
  if (Lmeantbisccp)        cfg%out_list(i) = 'meantbisccp'
  i = i+1
  if (Lmeantbclrisccp)        cfg%out_list(i) = 'meantbclrisccp'
  i = i+1
  if (LclMISR)        cfg%out_list(i) = 'clMISR'
  i = i+1
  if (Lclcalipso2)        cfg%out_list(i) = 'clcalipso2'
  i = i+1
  if (Lcltlidarradar)        cfg%out_list(i) = 'cltlidarradar'
  i = i+1
  if (Lcloudsat_tcc)        cfg%out_list(i) = 'cloudsat_tcc'
  i = i+1
  if (Lcloudsat_tcc2)        cfg%out_list(i) = 'cloudsat_tcc2'
  i = i+1
  if (Lfracout)        cfg%out_list(i) = 'fracout'
  i = i+1
  if (LlidarBetaMol532)        cfg%out_list(i) = 'lidarBetaMol532'
  i = i+1
  if (Lcltmodis)        cfg%out_list(i) = 'cltmodis'
  i = i+1
  if (Lclwmodis)        cfg%out_list(i) = 'clwmodis'
  i = i+1
  if (Lclimodis)        cfg%out_list(i) = 'climodis'
  i = i+1
  if (Lclhmodis)        cfg%out_list(i) = 'clhmodis'
  i = i+1
  if (Lclmmodis)        cfg%out_list(i) = 'clmmodis'
  i = i+1
  if (Lcllmodis)        cfg%out_list(i) = 'cllmodis'
  i = i+1
  if (Ltautmodis)        cfg%out_list(i) = 'tautmodis'
  i = i+1
  if (Ltauwmodis)        cfg%out_list(i) = 'tauwmodis'
  i = i+1
  if (Ltauimodis)        cfg%out_list(i) = 'tauimodis'
  i = i+1
  if (Ltautlogmodis)        cfg%out_list(i) = 'tautlogmodis'
  i = i+1
  if (Ltauwlogmodis)        cfg%out_list(i) = 'tauwlogmodis'
  i = i+1
  if (Ltauilogmodis)        cfg%out_list(i) = 'tauilogmodis'
  i = i+1
  if (Lreffclwmodis)        cfg%out_list(i) = 'reffclwmodis'
  i = i+1
  if (Lreffclimodis)        cfg%out_list(i) = 'reffclimodis'
  i = i+1
  if (Lpctmodis)        cfg%out_list(i) = 'pctmodis'
  i = i+1
  if (Llwpmodis)        cfg%out_list(i) = 'lwpmodis'
  i = i+1
  if (Liwpmodis)        cfg%out_list(i) = 'iwpmodis'
  i = i+1
  if (Lclmodis)        cfg%out_list(i) = 'clmodis'
  i = i+1
  if (Ltbrttov)        cfg%out_list(i) = 'tbrttov'
  i = i+1
  if (Lptradarflag0)        cfg%out_list(i) = 'ptradarflag0'
  i = i+1
  if (Lptradarflag1)        cfg%out_list(i) = 'ptradarflag1'
  i = i+1
  if (Lptradarflag2)        cfg%out_list(i) = 'ptradarflag2'
  i = i+1
  if (Lptradarflag3)        cfg%out_list(i) = 'ptradarflag3'
  i = i+1
  if (Lptradarflag4)        cfg%out_list(i) = 'ptradarflag4'
  i = i+1
  if (Lptradarflag5)        cfg%out_list(i) = 'ptradarflag5'
  i = i+1
  if (Lptradarflag6)        cfg%out_list(i) = 'ptradarflag6'
  i = i+1
  if (Lptradarflag7)        cfg%out_list(i) = 'ptradarflag7'
  i = i+1
  if (Lptradarflag8)        cfg%out_list(i) = 'ptradarflag8'
  i = i+1
  if (Lptradarflag9)        cfg%out_list(i) = 'ptradarflag9'
  i = i+1
  if (Lradarpia)        cfg%out_list(i) = 'radarpia'
!  i = i+1 !si on laisse, 108 au lieu de 107

  if (i.gt.107) then
     print *, 'COSP_IO: wrong number of output diagnostics'
     print *, i,107
     stop
  endif


  ! Copy diagnostic flags to cfg structure

   cfg%Lpctisccp = Lpctisccp
   cfg%Lclisccp = Lclisccp
   cfg%Lboxptopisccp = Lboxptopisccp
   cfg%Lboxtauisccp = Lboxtauisccp
   cfg%Ltauisccp = Ltauisccp
   cfg%Lcltisccp = Lcltisccp
   cfg%Lmeantbisccp = Lmeantbisccp
   cfg%Lmeantbclrisccp = Lmeantbclrisccp
   cfg%Lalbisccp = Lalbisccp
   cfg%LclMISR = LclMISR
   cfg%Lcltmodis = Lcltmodis
   cfg%Lclwmodis = Lclwmodis
   cfg%Lclimodis = Lclimodis
   cfg%Lclhmodis = Lclhmodis
   cfg%Lclmmodis = Lclmmodis
   cfg%Lcllmodis = Lcllmodis
   cfg%Ltautmodis = Ltautmodis
   cfg%Ltauwmodis = Ltauwmodis
   cfg%Ltauimodis = Ltauimodis
   cfg%Ltautlogmodis = Ltautlogmodis
   cfg%Ltauwlogmodis = Ltauwlogmodis
   cfg%Ltauilogmodis = Ltauilogmodis
   cfg%Lreffclwmodis = Lreffclwmodis
   cfg%Lreffclimodis = Lreffclimodis
   cfg%Lpctmodis = Lpctmodis
   cfg%Llwpmodis = Llwpmodis
   cfg%Liwpmodis = Liwpmodis
   cfg%Lclmodis = Lclmodis
   cfg%Latb532 = Latb532
   cfg%Latb532gr = Latb532gr
   cfg%Latb355 = Latb355
   cfg%LlidarBetaMol532 = LlidarBetaMol532
   cfg%LlidarBetaMol532gr = LlidarBetaMol532gr
   cfg%LlidarBetaMol355 = LlidarBetaMol355
   cfg%LcfadLidarsr532 = LcfadLidarsr532
   cfg%LcfadLidarsr532gr = LcfadLidarsr532gr
   cfg%LcfadLidarsr355 = LcfadLidarsr355
   cfg%Lclcalipso2 = Lclcalipso2
   cfg%Lclcalipso = Lclcalipso
   cfg%LclgrLidar532 = LclgrLidar532
   cfg%Lclatlid = Lclatlid
   cfg%Lclhcalipso = Lclhcalipso
   cfg%Lcllcalipso = Lcllcalipso
   cfg%Lclmcalipso = Lclmcalipso
   cfg%Lcltcalipso = Lcltcalipso
   cfg%LclhgrLidar532 = LclhgrLidar532
   cfg%LcllgrLidar532 = LcllgrLidar532
   cfg%LclmgrLidar532 = LclmgrLidar532
   cfg%LcltgrLidar532 = LcltgrLidar532
   cfg%Lclhatlid = Lclhatlid
   cfg%Lcllatlid = Lcllatlid
   cfg%Lclmatlid = Lclmatlid
   cfg%Lcltatlid = Lcltatlid
   cfg%Lcltlidarradar = Lcltlidarradar
   cfg%Lcloudsat_tcc = Lcloudsat_tcc
   cfg%Lcloudsat_tcc2 = Lcloudsat_tcc2
   cfg%Lclcalipsoliq = Lclcalipsoliq
   cfg%Lclcalipsoice = Lclcalipsoice
   cfg%Lclcalipsoun = Lclcalipsoun
   cfg%Lclcalipsotmp = Lclcalipsotmp
   cfg%Lclcalipsotmpliq = Lclcalipsotmpliq
   cfg%Lclcalipsotmpice = Lclcalipsotmpice
   cfg%Lclcalipsotmpun = Lclcalipsotmpun
   cfg%Lcltcalipsoliq = Lcltcalipsoliq
   cfg%Lcltcalipsoice = Lcltcalipsoice
   cfg%Lcltcalipsoun = Lcltcalipsoun
   cfg%Lclhcalipsoliq = Lclhcalipsoliq
   cfg%Lclhcalipsoice = Lclhcalipsoice
   cfg%Lclhcalipsoun = Lclhcalipsoun
   cfg%Lclmcalipsoliq = Lclmcalipsoliq
   cfg%Lclmcalipsoice = Lclmcalipsoice
   cfg%Lclmcalipsoun = Lclmcalipsoun
   cfg%Lcllcalipsoliq = Lcllcalipsoliq
   cfg%Lcllcalipsoice = Lcllcalipsoice
   cfg%Lcllcalipsoun = Lcllcalipsoun
   cfg%Lclopaquecalipso = Lclopaquecalipso
   cfg%Lclthincalipso = Lclthincalipso
   cfg%Lclzopaquecalipso = Lclzopaquecalipso
   cfg%Lclcalipsoopaque = Lclcalipsoopaque
   cfg%Lclcalipsothin = Lclcalipsothin
   cfg%Lclcalipsozopaque = Lclcalipsozopaque
   cfg%Lclcalipsoopacity = Lclcalipsoopacity
   cfg%Lclopaquetemp = Lclopaquetemp
   cfg%Lclthintemp = Lclthintemp
   cfg%Lclzopaquetemp = Lclzopaquetemp
   cfg%Lclopaquemeanz = Lclopaquemeanz
   cfg%Lclthinmeanz = Lclthinmeanz
   cfg%Lclthinemis = Lclthinemis
   cfg%Lclopaquemeanzse = Lclopaquemeanzse
   cfg%Lclthinmeanzse = Lclthinmeanzse
   cfg%Lclzopaquecalipsose = Lclzopaquecalipsose
   cfg%LcfadDbze94 = LcfadDbze94
   cfg%Ldbze94 = Ldbze94
   cfg%Lparasolrefl = Lparasolrefl
   cfg%Ltbrttov = Ltbrttov
   cfg%Lptradarflag0 = Lptradarflag0
   cfg%Lptradarflag1 = Lptradarflag1
   cfg%Lptradarflag2 = Lptradarflag2
   cfg%Lptradarflag3 = Lptradarflag3
   cfg%Lptradarflag4 = Lptradarflag4
   cfg%Lptradarflag5 = Lptradarflag5
   cfg%Lptradarflag6 = Lptradarflag6
   cfg%Lptradarflag7 = Lptradarflag7
   cfg%Lptradarflag8 = Lptradarflag8
   cfg%Lptradarflag9 = Lptradarflag9
   cfg%Lradarpia = Lradarpia
   cfg%Lfracout = Lfracout
  
 END SUBROUTINE READ_COSP_OUTPUT_NL


 SUBROUTINE read_xiosfieldactive(cfg)

#ifdef CPP_XIOS
    USE xios, ONLY: xios_field_is_active
#endif
  implicit none
  type(cosp_config),intent(inout) :: cfg
  integer :: i

#ifdef CPP_XIOS

!COSPv2 local variables
  logical :: &
       Lisccp      = .false., & ! Local on/off switch for simulators
       Lmodis      = .false., & !
       Lmisr       = .false., & !
       Lcalipso    = .false., & !
       LgrLidar532 = .false., & !
       Latlid      = .false., & !
       Lcloudsat   = .false., & !
       Lrttov      = .false., & !
       Lparasol    = .false., & !
       Lwrite_output=.false.    ! ecrire un fichier de sortie COSP ! RG voir si on garde

  logical :: Lcfaddbze94,Ldbze94,Latb532,LcfadLidarsr532,Lclcalipso,Lclhcalipso,         &
             Lcllcalipso,Lclmcalipso,Lcltcalipso,LparasolRefl,Lclcalipsoliq,             &
             Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice, &
             Lclcalipsotmpun,Lclhcalipsoliq,Lcllcalipsoliq,Lclmcalipsoliq,Lcltcalipsoliq,& 
             Lclhcalipsoice,Lcllcalipsoice,Lclmcalipsoice,Lcltcalipsoice,Lclhcalipsoun,  & 
             Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lclopaquecalipso,Lclthincalipso,  & 
             Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,        & 
             Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,Lclzopaquetemp,Lclopaquemeanz,  &
             Lclthinmeanz,Lclthinemis,Lclopaquemeanzse,Lclthinmeanzse,Lclzopaquecalipsose,&
             LlidarBetaMol532gr,LcfadLidarsr532gr,Latb532gr,LclgrLidar532,LclhgrLidar532,&
             LcllgrLidar532,LclmgrLidar532,LcltgrLidar532,LlidarBetaMol355,              &
             LcfadLidarsr355,Latb355,Lclatlid,Lclhatlid,Lcllatlid,Lclmatlid,Lcltatlid,   &
             Lalbisccp,Lboxptopisccp,Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,&
             Lmeantbisccp,Lmeantbclrisccp,LclMISR,Lclcalipso2,Lcltlidarradar,            &
             Lcloudsat_tcc, Lcloudsat_tcc2,Lfracout,   &
             LlidarBetaMol532,Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,         &
             Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis,Ltauwlogmodis,     &
             Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,    &
             Lclmodis,Ltbrttov,Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,  &
             Lptradarflag4,Lptradarflag5,Lptradarflag6,Lptradarflag7,Lptradarflag8,      &
             Lptradarflag9,Lradarpia

        
  character(len=32) :: out_list(107)

  do i=1,107
    cfg%out_list(i)=''
  enddo

! diagnostics
   Lpctisccp=.false.
   Lclisccp=.false.
   Lboxptopisccp=.false.
   Lboxtauisccp=.false.
   Ltauisccp=.false.
   Lcltisccp=.false.
   Lmeantbisccp=.false.
   Lmeantbclrisccp=.false.
   Lalbisccp=.false.
   LclMISR=.false.
   Lcltmodis=.false.
   Lclwmodis=.false.
   Lclimodis=.false.
   Lclhmodis=.false.
   Lclmmodis=.false.
   Lcllmodis=.false.
   Ltautmodis=.false.
   Ltauwmodis=.false.
   Ltauimodis=.false.
   Ltautlogmodis=.false.
   Ltauwlogmodis=.false.
   Ltauilogmodis=.false.
   Lreffclwmodis=.false.
   Lreffclimodis=.false.
   Lpctmodis=.false.
   Llwpmodis=.false.
   Liwpmodis=.false.
   Lclmodis=.false.
   Latb532=.false.
   Latb532gr=.false.
   Latb355=.false.
   LlidarBetaMol532=.false.
   LlidarBetaMol532gr=.false.
   LlidarBetaMol355=.false.
   LcfadLidarsr532=.false.
   LcfadLidarsr532gr=.false.
   LcfadLidarsr355 =.false.
   Lclcalipso2=.false.
   Lclcalipso=.false.
   LclgrLidar532=.false.
   Lclatlid=.false.
   Lclhcalipso=.false.
   Lcllcalipso=.false.
   Lclmcalipso=.false.
   Lcltcalipso=.false.
   LclhgrLidar532=.false.
   LcllgrLidar532=.false.
   LclmgrLidar532 =.false.
   LcltgrLidar532=.false.
   Lclhatlid=.false.
   Lcllatlid=.false.
   Lclmatlid=.false.
   Lcltatlid=.false.
   Lcltlidarradar=.false.
   Lcloudsat_tcc=.false.
   Lcloudsat_tcc2=.false.
   Lclcalipsoliq=.false.
   Lclcalipsoice=.false.
   Lclcalipsoun=.false.
   Lclcalipsotmp=.false.
   Lclcalipsotmpliq=.false.
   Lclcalipsotmpice=.false.
   Lclcalipsotmpun=.false.
   Lcltcalipsoliq=.false.
   Lcltcalipsoice=.false.
   Lcltcalipsoun=.false.
   Lclhcalipsoliq=.false.
   Lclhcalipsoice=.false.
   Lclhcalipsoun=.false.
   Lclmcalipsoliq=.false.
   Lclmcalipsoice=.false.
   Lclmcalipsoun=.false.
   Lcllcalipsoliq=.false.
   Lcllcalipsoice=.false.
   Lcllcalipsoun =.false.
   Lclopaquecalipso=.false.
   Lclthincalipso=.false.
   Lclzopaquecalipso=.false.
   Lclcalipsoopaque=.false.
   Lclcalipsothin=.false.
   Lclcalipsozopaque=.false.
   Lclcalipsoopacity=.false.
   Lclopaquetemp=.false.
   Lclthintemp=.false.
   Lclzopaquetemp=.false.
   Lclopaquemeanz=.false.
   Lclthinmeanz=.false.
   Lclthinemis=.false.
   Lclopaquemeanzse=.false.
   Lclthinmeanzse=.false.
   Lclzopaquecalipsose=.false.
   LcfadDbze94=.false.
   Ldbze94=.false.
   Lparasolrefl=.false.
   Ltbrttov=.false.
   Lptradarflag0=.false.
   Lptradarflag1=.false.
   Lptradarflag2=.false.
   Lptradarflag3=.false.
   Lptradarflag4=.false.
   Lptradarflag5=.false.
   Lptradarflag6=.false.
   Lptradarflag7=.false.
   Lptradarflag8=.false.
   Lptradarflag9=.false.
   Lradarpia=.false.
   Lfracout=.false.


  IF (is_master) THEN
! VEREFIER LES CHAMPS DEMANDES DANS .XML
! 2. Si champs active dans .xml alors mettre la cles de sortie en true
 IF (xios_field_is_active("cllcalipso")) Lcllcalipso=.TRUE.
 IF (xios_field_is_active("clmcalipso")) Lclmcalipso=.TRUE.
 IF (xios_field_is_active("clhcalipso")) Lclhcalipso=.TRUE.
 IF (xios_field_is_active("cltcalipso")) Lcltcalipso=.TRUE.
 IF (xios_field_is_active("cllcalipsoice")) Lcllcalipsoice=.TRUE.
 IF (xios_field_is_active("clmcalipsoice")) Lclmcalipsoice=.TRUE.
 IF (xios_field_is_active("clhcalipsoice")) Lclhcalipsoice=.TRUE.
 IF (xios_field_is_active("cltcalipsoice")) Lcltcalipsoice=.TRUE.
 IF (xios_field_is_active("cllcalipsoliq")) Lcllcalipsoliq=.TRUE.
 IF (xios_field_is_active("clmcalipsoliq")) Lclmcalipsoliq=.TRUE.
 IF (xios_field_is_active("clhcalipsoliq")) Lclhcalipsoliq=.TRUE.
 IF (xios_field_is_active("cltcalipsoliq")) Lcltcalipsoliq=.TRUE.
 IF (xios_field_is_active("cllcalipsoun")) Lcllcalipsoun=.TRUE.
 IF (xios_field_is_active("clmcalipsoun")) Lclmcalipsoun=.TRUE.
 IF (xios_field_is_active("clhcalipsoun")) Lclhcalipsoun=.TRUE.
 IF (xios_field_is_active("cltcalipsoun")) Lcltcalipsoun=.TRUE.
 IF (xios_field_is_active("clcalipso")) Lclcalipso=.TRUE.
 IF (xios_field_is_active("clcalipsoice")) Lclcalipsoice=.TRUE.
 IF (xios_field_is_active("clcalipsoliq")) Lclcalipsoliq=.TRUE.
 IF (xios_field_is_active("clcalipsoun")) Lclcalipsoun=.TRUE.
 IF (xios_field_is_active("clcalipsotmp")) Lclcalipsotmp=.TRUE.
 IF (xios_field_is_active("clcalipsotmpice")) Lclcalipsotmpice=.TRUE.
 IF (xios_field_is_active("clcalipsotmpliq")) Lclcalipsotmpliq=.TRUE.
 IF (xios_field_is_active("clcalipsotmpun")) Lclcalipsotmpun=.TRUE.
! IF (xios_field_is_active("parasol_crefl")) cfg%LparasolRefl=.TRUE.
! IF (xios_field_is_active("Ncrefl")) cfg%LparasolRefl=.TRUE.
 IF (xios_field_is_active("cfadLidarsr532")) LcfadLidarsr532=.TRUE.
 IF (xios_field_is_active("atb532")) Latb532=.TRUE.
 IF (xios_field_is_active("lidarBetaMol532")) LlidarBetaMol532=.TRUE.
 IF (xios_field_is_active("clopaquecalipso")) Lclopaquecalipso=.TRUE.
 IF (xios_field_is_active("clthincalipso")) Lclthincalipso=.TRUE.
 IF (xios_field_is_active("clzopaquecalipso")) Lclzopaquecalipso=.TRUE.
 IF (xios_field_is_active("clcalipsoopaque")) Lclcalipsoopaque=.TRUE.
 IF (xios_field_is_active("clcalipsothin")) Lclcalipsothin=.TRUE.
 IF (xios_field_is_active("clcalipsozopaque")) Lclcalipsozopaque=.TRUE.
 IF (xios_field_is_active("clcalipsoopacity")) Lclcalipsoopacity=.TRUE.
 IF (xios_field_is_active("clopaquetemp")) Lclopaquetemp=.TRUE.
 IF (xios_field_is_active("clthintemp")) Lclthintemp=.TRUE.
 IF (xios_field_is_active("clzopaquetemp")) Lclzopaquetemp=.TRUE.
 IF (xios_field_is_active("clopaquemeanz")) Lclopaquemeanz=.TRUE.
 IF (xios_field_is_active("clthinmeanz")) Lclthinmeanz=.TRUE.
 IF (xios_field_is_active("clthinemis")) Lclthinemis=.TRUE.
 IF (xios_field_is_active("clopaquemeanzse")) Lclopaquemeanzse=.TRUE.
 IF (xios_field_is_active("clthinmeanzse")) Lclthinmeanzse=.TRUE.
 IF (xios_field_is_active("clzopaquecalipsose")) Lclzopaquecalipsose=.TRUE.
!!!! 43 champs CALIPSO

 IF (xios_field_is_active("parasolGrid_refl")) LparasolRefl=.TRUE.
 IF (xios_field_is_active("parasolPix_refl")) LparasolRefl=.TRUE.
!!!! 2 champs PARASOL controles par une seule clef logique

 IF (xios_field_is_active("lidarBetaMol532gr")) LlidarBetaMol532gr=.TRUE.
 IF (xios_field_is_active("cfadLidarsr532gr")) LcfadLidarsr532gr=.TRUE.
 IF (xios_field_is_active("atb532gr")) Latb532gr=.TRUE.
 IF (xios_field_is_active("clgrLidar532")) LclgrLidar532=.TRUE.
 IF (xios_field_is_active("clhgrLidar532")) LclhgrLidar532=.TRUE.
 IF (xios_field_is_active("cllgrLidar532")) LcllgrLidar532=.TRUE.
 IF (xios_field_is_active("clmgrLidar532")) LclmgrLidar532=.TRUE.
 IF (xios_field_is_active("cltgrLidar532")) LcltgrLidar532=.TRUE.
!!!! 8 champs GROUND LIDAR 532 nm

 IF (xios_field_is_active("lidarBetaMol355")) LlidarBetaMol355=.TRUE.
 IF (xios_field_is_active("cfadLidarsr355")) LcfadLidarsr355=.TRUE.
 IF (xios_field_is_active("atb355")) Latb355=.TRUE.
 IF (xios_field_is_active("clatlid")) Lclatlid=.TRUE.
 IF (xios_field_is_active("clhatlid")) Lclhatlid=.TRUE.
 IF (xios_field_is_active("cllatlid")) Lcllatlid=.TRUE.
 IF (xios_field_is_active("clmatlid")) Lclmatlid=.TRUE.
 IF (xios_field_is_active("cltatlid")) Lcltatlid=.TRUE.
!!!! 8 champs ATLID

 IF (xios_field_is_active("cfadDbze94")) LcfadDbze94=.TRUE.
 IF (xios_field_is_active("dbze94")) Ldbze94=.TRUE.
 IF (xios_field_is_active("ptradarflag0")) Lptradarflag0=.TRUE.
 IF (xios_field_is_active("ptradarflag1")) Lptradarflag1=.TRUE.
 IF (xios_field_is_active("ptradarflag2")) Lptradarflag2=.TRUE.
 IF (xios_field_is_active("ptradarflag3")) Lptradarflag3=.TRUE.
 IF (xios_field_is_active("ptradarflag4")) Lptradarflag4=.TRUE.
 IF (xios_field_is_active("ptradarflag5")) Lptradarflag5=.TRUE.
 IF (xios_field_is_active("ptradarflag6")) Lptradarflag6=.TRUE.
 IF (xios_field_is_active("ptradarflag7")) Lptradarflag7=.TRUE.
 IF (xios_field_is_active("ptradarflag8")) Lptradarflag8=.TRUE.
 IF (xios_field_is_active("ptradarflag9")) Lptradarflag9=.TRUE.
 IF (xios_field_is_active("radarpia")) Lradarpia=.TRUE.
!!! 13 champs CLOUDSAT

 IF (xios_field_is_active("cltlidarradar")) Lcltlidarradar=.TRUE.
 IF (xios_field_is_active("clcalipso2")) Lclcalipso2=.TRUE.
 IF (xios_field_is_active("cloudsat_tcc")) Lcloudsat_tcc=.TRUE. !COSPv2
 IF (xios_field_is_active("cloudsat_tcc2")) Lcloudsat_tcc2=.TRUE. !COSPv2
!!! 4 champs CLOUDSAT et CALIPSO

 IF (xios_field_is_active("clisccp")) Lclisccp=.TRUE.
 IF (xios_field_is_active("boxtauisccp")) Lboxtauisccp=.TRUE.
 IF (xios_field_is_active("boxptopisccp")) Lboxptopisccp=.TRUE.
 IF (xios_field_is_active("cltisccp")) Lcltisccp=.TRUE.
 IF (xios_field_is_active("pctisccp")) Lpctisccp=.TRUE.
 IF (xios_field_is_active("tauisccp")) Ltauisccp=.TRUE.
 IF (xios_field_is_active("albisccp")) Lalbisccp=.TRUE.
 IF (xios_field_is_active("meantbisccp")) Lmeantbisccp=.TRUE.
 IF (xios_field_is_active("meantbclrisccp")) Lmeantbclrisccp=.TRUE.
!!! 9 champs ISCCP

 IF (xios_field_is_active("misr_fq")) LclMISR=.TRUE.
 IF (xios_field_is_active("misr_meanztop")) LclMISR=.TRUE.
 IF (xios_field_is_active("misr_cldarea")) LclMISR=.TRUE.
!!! 3 champs MISR controles par une seule clef logique

 IF (xios_field_is_active("cllmodis")) Lcllmodis=.TRUE.
 IF (xios_field_is_active("clmmodis")) Lclmmodis=.TRUE.
 IF (xios_field_is_active("clhmodis")) Lclhmodis=.TRUE.
 IF (xios_field_is_active("cltmodis")) Lcltmodis=.TRUE.
 IF (xios_field_is_active("clwmodis")) Lclwmodis=.TRUE.
 IF (xios_field_is_active("climodis")) Lclimodis=.TRUE.
 IF (xios_field_is_active("tautmodis")) Ltautmodis=.TRUE.
 IF (xios_field_is_active("tauwmodis")) Ltauwmodis=.TRUE.
 IF (xios_field_is_active("tauimodis")) Ltauimodis=.TRUE.
 IF (xios_field_is_active("tautlogmodis")) Ltautlogmodis=.TRUE.
 IF (xios_field_is_active("tauilogmodis")) Ltauilogmodis=.TRUE.
 IF (xios_field_is_active("tauwlogmodis")) Ltauwlogmodis=.TRUE.
 IF (xios_field_is_active("reffclwmodis")) Lreffclwmodis=.TRUE.
 IF (xios_field_is_active("reffclimodis")) Lreffclimodis=.TRUE.
 IF (xios_field_is_active("pctmodis")) Lpctmodis=.TRUE.
 IF (xios_field_is_active("lwpmodis")) Llwpmodis=.TRUE.
 IF (xios_field_is_active("iwpmodis")) Liwpmodis=.TRUE.
 IF (xios_field_is_active("modis_ot_vs_ctp")) Lclmodis=.TRUE.
 IF (xios_field_is_active("modis_ot_vs_reffice")) Lclmodis=.TRUE.
 IF (xios_field_is_active("modis_ot_vs_reffliq")) Lclmodis=.TRUE.
! IF (xios_field_is_active("crimodis")) Lcrimodis=.TRUE. !COSPv2
! IF (xios_field_is_active("crlmodis")) Lcrlmodis=.TRUE. !COSPv2
!!! 20 champs MODIS dont 3 champs controles par la clef logique "clmodis"

! IF (xios_field_is_active("tbrttov")) cfg%Ltbrttov=.TRUE.
! IF (xios_field_is_active("fracout")) Lfracout=.TRUE.

!!! 106 diagnostics en tout dans COSPv2 + fracout (107)


! 2.  si champs demande alors activer le simulateur correspondant
   IF (xios_field_is_active("cllcalipso").OR. &
       xios_field_is_active("clmcalipso").OR. &
       xios_field_is_active("clhcalipso").OR. &
       xios_field_is_active("cltcalipso").OR. &
       xios_field_is_active("cllcalipsoice").OR. &
       xios_field_is_active("clmcalipsoice").OR. &
       xios_field_is_active("clhcalipsoice").OR. &
       xios_field_is_active("cltcalipsoice").OR. &
       xios_field_is_active("cllcalipsoliq").OR. &
       xios_field_is_active("clmcalipsoliq").OR. &
       xios_field_is_active("clhcalipsoliq").OR. &
       xios_field_is_active("cltcalipsoliq").OR. &
       xios_field_is_active("cllcalipsoun").OR. &
       xios_field_is_active("clmcalipsoun").OR. &
       xios_field_is_active("clhcalipsoun").OR. &
       xios_field_is_active("cltcalipsoun").OR. &
       xios_field_is_active("clcalipso").OR. &
       xios_field_is_active("clcalipsoice").OR. &
       xios_field_is_active("clcalipsoliq").OR. &
       xios_field_is_active("clcalipsoun").OR. &
       xios_field_is_active("clcalipsotmp").OR. &
       xios_field_is_active("clcalipsotmpice").OR. &
       xios_field_is_active("clcalipsotmpliq").OR. &
       xios_field_is_active("clcalipsotmpun").OR. &
       xios_field_is_active("cfadLidarsr532").OR. &
       xios_field_is_active("atb532").OR. &
       xios_field_is_active("lidarBetaMol532").OR. &
       xios_field_is_active("clopaquecalipso").OR. &
       xios_field_is_active("clthincalipso").OR. &
       xios_field_is_active("clzopaquecalipso").OR. &
       xios_field_is_active("clcalipsoopaque").OR. &
       xios_field_is_active("clcalipsothin").OR. &
       xios_field_is_active("clcalipsozopaque").OR. &
       xios_field_is_active("clcalipsoopacity").OR. &
       xios_field_is_active("clopaquetemp").OR. &
       xios_field_is_active("clthintemp").OR. &
       xios_field_is_active("clzopaquetemp").OR. &
       xios_field_is_active("clopaquemeanz").OR. &
       xios_field_is_active("clthinmeanz").OR. &
       xios_field_is_active("clthinemis").OR. &
       xios_field_is_active("clopaquemeanzse").OR. &
       xios_field_is_active("clthinmeanzse").OR. &
       xios_field_is_active("clzopaquecalipsose"))  Lcalipso=.TRUE.

    IF (xios_field_is_active("parasolGrid_refl").OR. &
       xios_field_is_active("parasolPix_refl")) Lparasol=.TRUE.

   IF (xios_field_is_active("lidarBetaMol532gr").OR. &
       xios_field_is_active("cfadLidarsr532gr").OR. &
       xios_field_is_active("atb532gr").OR. &
       xios_field_is_active("clgrLidar532").OR. &
       xios_field_is_active("clhgrLidar532").OR. &
       xios_field_is_active("cllgrLidar532").OR. &
       xios_field_is_active("clmgrLidar532").OR. &
       xios_field_is_active("cltgrLidar532"))  LgrLidar532=.TRUE.

   IF (xios_field_is_active("lidarBetaMol355").OR. &
       xios_field_is_active("cfadLidarsr355").OR. &
       xios_field_is_active("atb355").OR. &
       xios_field_is_active("clatlid").OR. &
       xios_field_is_active("clhatlid").OR. &
       xios_field_is_active("cllatlid").OR. &
       xios_field_is_active("clmatlid").OR. &
       xios_field_is_active("cltatlid"))  Latlid=.TRUE.

    IF (xios_field_is_active("cfadDbze94").OR. &
      xios_field_is_active("dbze94").OR. &
      xios_field_is_active("ptradarflag0").OR. &
      xios_field_is_active("ptradarflag1").OR. &
      xios_field_is_active("ptradarflag2").OR. &
      xios_field_is_active("ptradarflag3").OR. &
      xios_field_is_active("ptradarflag4").OR. &
      xios_field_is_active("ptradarflag5").OR. &
      xios_field_is_active("ptradarflag6").OR. &
      xios_field_is_active("ptradarflag7").OR. &
      xios_field_is_active("ptradarflag8").OR. &
      xios_field_is_active("ptradarflag9").OR. &
      xios_field_is_active("radarpia")) Lcloudsat=.TRUE.

    IF (xios_field_is_active("cltlidarradar").OR. &
      xios_field_is_active("clcalipso2").OR. &
      xios_field_is_active("cloudsat_tcc").OR. &
      xios_field_is_active("cloudsat_tcc2")) THEN
               Lcloudsat=.TRUE.
               Lcalipso=.TRUE.
    ENDIF

    IF (xios_field_is_active("clisccp").OR. &     !COSPv2 nouveau nom
       xios_field_is_active("boxtauisccp").OR. &
       xios_field_is_active("boxptopisccp").OR. &
       xios_field_is_active("cltisccp").OR. &     !COSPv2 nouveau nom
       xios_field_is_active("pctisccp").OR. &     !COSPv2 nouveau nom
       xios_field_is_active("tauisccp").OR. &
       xios_field_is_active("albisccp").OR. &
       xios_field_is_active("meantbisccp").OR. &
       xios_field_is_active("meantbclrisccp")) Lisccp=.TRUE.

    IF (xios_field_is_active("misr_fq").OR. &
       xios_field_is_active("misr_meanztop").OR. &
       xios_field_is_active("misr_cldarea")) Lmisr=.TRUE.

    IF (xios_field_is_active("cllmodis").OR. &
       xios_field_is_active("clmmodis").OR. &
       xios_field_is_active("clhmodis").OR. &
       xios_field_is_active("cltmodis").OR. &
       xios_field_is_active("clwmodis").OR. &
       xios_field_is_active("climodis").OR. &
       xios_field_is_active("tautmodis").OR. &
       xios_field_is_active("tauwmodis").OR. &
       xios_field_is_active("tauimodis").OR. &
       xios_field_is_active("tautlogmodis").OR. &
       xios_field_is_active("tauilogmodis").OR. &
       xios_field_is_active("tauwlogmodis").OR. &
       xios_field_is_active("reffclwmodis").OR. &
       xios_field_is_active("reffclimodis").OR. &
       xios_field_is_active("pctmodis").OR. &
       xios_field_is_active("lwpmodis").OR. &
       xios_field_is_active("iwpmodis").OR. &
       xios_field_is_active("modis_ot_vs_ctp").OR. &
       xios_field_is_active("modis_ot_vs_reffice").OR. &
       xios_field_is_active("modis_ot_vs_reffliq")) Lmodis=.TRUE.
!       xios_field_is_active("crimodis").OR. &              !COSPv2 n'existe plus
!       xios_field_is_active("crlmodis")) Lmodis_sim=.TRUE. !COSPv2 n'existe plus

  ENDIF !   (is_master) 

!$OMP BARRIER

  CALL bcast(Lcloudsat)
  CALL bcast(Lcalipso)
  CALL bcast(LgrLidar532)
  CALL bcast(Latlid)
  CALL bcast(Lisccp)
  CALL bcast(Lmodis)
  CALL bcast(Lmisr)
  CALL bcast(Lparasol)
  CALL bcast(Lrttov)

  CALL bcast(Lcfaddbze94)
  CALL bcast(Ldbze94)
  CALL bcast(Latb532)
  CALL bcast(LcfadLidarsr532)
  CALL bcast(Lclcalipso)
  CALL bcast(Lclhcalipso)
  CALL bcast(Lcllcalipso)
  CALL bcast(Lclmcalipso)
  CALL bcast(Lcltcalipso)
  CALL bcast(LparasolRefl)
  CALL bcast(Lclcalipsoliq)
  CALL bcast(Lclcalipsoice)
  CALL bcast(Lclcalipsoun)
  CALL bcast(Lclcalipsotmp)
  CALL bcast(Lclcalipsotmpliq)
  CALL bcast(Lclcalipsotmpice)
  CALL bcast(Lclcalipsotmpun)
  CALL bcast(Lclhcalipsoliq)
  CALL bcast(Lcllcalipsoliq)
  CALL bcast(Lclmcalipsoliq)
  CALL bcast(Lcltcalipsoliq)
  CALL bcast(Lclhcalipsoice)
  CALL bcast(Lcllcalipsoice)
  CALL bcast(Lclmcalipsoice)
  CALL bcast(Lcltcalipsoice)
  CALL bcast(Lclhcalipsoun)
  CALL bcast(Lcllcalipsoun)
  CALL bcast(Lclmcalipsoun)
  CALL bcast(Lcltcalipsoun)
  CALL bcast(Lclopaquecalipso)
  CALL bcast(Lclthincalipso)
  CALL bcast(Lclzopaquecalipso)
  CALL bcast(Lclcalipsoopaque)
  CALL bcast(Lclcalipsothin)
  CALL bcast(Lclcalipsozopaque)
  CALL bcast(Lclcalipsoopacity)
  CALL bcast(Lclopaquetemp)
  CALL bcast(Lclthintemp)
  CALL bcast(Lclzopaquetemp)
  CALL bcast(Lclopaquemeanz)
  CALL bcast(Lclthinmeanz)
  CALL bcast(Lclthinemis)
  CALL bcast(Lclopaquemeanzse)
  CALL bcast(Lclthinmeanzse)
  CALL bcast(Lclzopaquecalipsose)
  CALL bcast(LlidarBetaMol532gr)
  CALL bcast(LcfadLidarsr532gr)
  CALL bcast(Latb532gr)
  CALL bcast(LclgrLidar532)
  CALL bcast(LclhgrLidar532)
  CALL bcast(LcllgrLidar532)
  CALL bcast(LclmgrLidar532)
  CALL bcast(LcltgrLidar532)
  CALL bcast(LlidarBetaMol355)
  CALL bcast(LcfadLidarsr355)
  CALL bcast(Latb355)
  CALL bcast(Lclatlid)
  CALL bcast(Lclhatlid)
  CALL bcast(Lcllatlid)
  CALL bcast(Lclmatlid)
  CALL bcast(Lcltatlid)
  CALL bcast(Lalbisccp)
  CALL bcast(Lboxptopisccp)
  CALL bcast(Lboxtauisccp)
  CALL bcast(Lpctisccp)
  CALL bcast(Lclisccp)
  CALL bcast(Ltauisccp)
  CALL bcast(Lcltisccp)
  CALL bcast(Lmeantbisccp)
  CALL bcast(Lmeantbclrisccp)
  CALL bcast(LclMISR)
  CALL bcast(Lclcalipso2)
  CALL bcast(Lcltlidarradar)
  CALL bcast(Lcloudsat_tcc)
  CALL bcast(Lcloudsat_tcc2)
  CALL bcast(Lfracout)
  CALL bcast(LlidarBetaMol532)
  CALL bcast(Lcltmodis)
  CALL bcast(Lclwmodis)
  CALL bcast(Lclimodis)
  CALL bcast(Lclhmodis)
  CALL bcast(Lclmmodis)
  CALL bcast(Lcllmodis)
  CALL bcast(Ltautmodis)
  CALL bcast(Ltauwmodis)
  CALL bcast(Ltauimodis)
  CALL bcast(Ltautlogmodis)
  CALL bcast(Ltauwlogmodis)
  CALL bcast(Ltauilogmodis)
  CALL bcast(Lreffclwmodis)
  CALL bcast(Lreffclimodis)
  CALL bcast(Lpctmodis)
  CALL bcast(Llwpmodis)
  CALL bcast(Liwpmodis)
  CALL bcast(Lclmodis)
  CALL bcast(Ltbrttov)
  CALL bcast(Lptradarflag0)
  CALL bcast(Lptradarflag1)
  CALL bcast(Lptradarflag2)
  CALL bcast(Lptradarflag3)
  CALL bcast(Lptradarflag4)
  CALL bcast(Lptradarflag5)
  CALL bcast(Lptradarflag6)
  CALL bcast(Lptradarflag7)
  CALL bcast(Lptradarflag8)
  CALL bcast(Lptradarflag9)
  CALL bcast(Lradarpia)

!    IF (xios_field_is_active("tbrttov")) cfg%Lrttov_sim=.TRUE.


  ! Copy diagnostic flags to cfg structure
! ISCCP simulator  (9 diagnostics)
   cfg%Lpctisccp = Lpctisccp
   cfg%Lclisccp = Lclisccp
   cfg%Lboxptopisccp = Lboxptopisccp
   cfg%Lboxtauisccp = Lboxtauisccp
   cfg%Ltauisccp = Ltauisccp
   cfg%Lcltisccp = Lcltisccp
   cfg%Lmeantbisccp = Lmeantbisccp
   cfg%Lmeantbclrisccp = Lmeantbclrisccp
   cfg%Lalbisccp = Lalbisccp

! CloudSat simulator  (13 diagnostics) 
   cfg%Ldbze94 = Ldbze94
   cfg%LcfadDbze94 = LcfadDbze94
   cfg%Lptradarflag0 = Lptradarflag0
   cfg%Lptradarflag1 = Lptradarflag1
   cfg%Lptradarflag2 = Lptradarflag2
   cfg%Lptradarflag3 = Lptradarflag3
   cfg%Lptradarflag4 = Lptradarflag4
   cfg%Lptradarflag5 = Lptradarflag5
   cfg%Lptradarflag6 = Lptradarflag6
   cfg%Lptradarflag7 = Lptradarflag7
   cfg%Lptradarflag8 = Lptradarflag8
   cfg%Lptradarflag9 = Lptradarflag9
   cfg%Lradarpia = Lradarpia

! Cloudsat et Calipso  (4 diagnostics)
   cfg%Lclcalipso2 = Lclcalipso2
   cfg%Lcltlidarradar = Lcltlidarradar
   cfg%Lcloudsat_tcc = Lcloudsat_tcc
   cfg%Lcloudsat_tcc2 = Lcloudsat_tcc2

! CALIPSO simulator  (43 diagnostics)
   cfg%Latb532 = Latb532
   cfg%LlidarBetaMol532 = LlidarBetaMol532
   cfg%LcfadLidarsr532 = LcfadLidarsr532
   cfg%Lclcalipso = Lclcalipso
   cfg%Lclhcalipso = Lclhcalipso
   cfg%Lcllcalipso = Lcllcalipso
   cfg%Lclmcalipso = Lclmcalipso
   cfg%Lcltcalipso = Lcltcalipso
   cfg%Lclcalipsoliq = Lclcalipsoliq
   cfg%Lclcalipsoice = Lclcalipsoice
   cfg%Lclcalipsoun = Lclcalipsoun
   cfg%Lclcalipsotmp = Lclcalipsotmp
   cfg%Lclcalipsotmpliq = Lclcalipsotmpliq
   cfg%Lclcalipsotmpice = Lclcalipsotmpice
   cfg%Lclcalipsotmpun = Lclcalipsotmpun
   cfg%Lcltcalipsoliq = Lcltcalipsoliq
   cfg%Lcltcalipsoice = Lcltcalipsoice
   cfg%Lcltcalipsoun = Lcltcalipsoun
   cfg%Lclhcalipsoliq = Lclhcalipsoliq
   cfg%Lclhcalipsoice = Lclhcalipsoice
   cfg%Lclhcalipsoun = Lclhcalipsoun
   cfg%Lclmcalipsoliq = Lclmcalipsoliq
   cfg%Lclmcalipsoice = Lclmcalipsoice
   cfg%Lclmcalipsoun = Lclmcalipsoun
   cfg%Lcllcalipsoliq = Lcllcalipsoliq
   cfg%Lcllcalipsoice = Lcllcalipsoice
   cfg%Lcllcalipsoun = Lcllcalipsoun
   cfg%Lclopaquecalipso = Lclopaquecalipso
   cfg%Lclthincalipso = Lclthincalipso
   cfg%Lclzopaquecalipso = Lclzopaquecalipso
   cfg%Lclcalipsoopaque = Lclcalipsoopaque
   cfg%Lclcalipsothin = Lclcalipsothin
   cfg%Lclcalipsozopaque = Lclcalipsozopaque
   cfg%Lclcalipsoopacity = Lclcalipsoopacity
   cfg%Lclopaquetemp = Lclopaquetemp
   cfg%Lclthintemp = Lclthintemp
   cfg%Lclzopaquetemp = Lclzopaquetemp
   cfg%Lclopaquemeanz = Lclopaquemeanz
   cfg%Lclthinmeanz = Lclthinmeanz
   cfg%Lclthinemis = Lclthinemis
   cfg%Lclopaquemeanzse = Lclopaquemeanzse
   cfg%Lclthinmeanzse = Lclthinmeanzse
   cfg%Lclzopaquecalipsose = Lclzopaquecalipsose

! PARASOL simulator  (2 diagnostics, 1 clef)
   cfg%Lparasolrefl = Lparasolrefl

! GROUND LIDAR 532 nm simulator  (8 diagnostics)
   cfg%Latb532gr = Latb532gr
   cfg%LlidarBetaMol532gr = LlidarBetaMol532gr
   cfg%LcfadLidarsr532gr = LcfadLidarsr532gr
   cfg%LclgrLidar532 = LclgrLidar532
   cfg%LclhgrLidar532 = LclhgrLidar532
   cfg%LcllgrLidar532 = LcllgrLidar532
   cfg%LclmgrLidar532 = LclmgrLidar532
   cfg%LcltgrLidar532 = LcltgrLidar532

! ATLID simulator  (8 diagnostics)
   cfg%Latb355 = Latb355
   cfg%LlidarBetaMol355 = LlidarBetaMol355
   cfg%LcfadLidarsr355 = LcfadLidarsr355
   cfg%Lclatlid = Lclatlid
   cfg%Lclhatlid = Lclhatlid
   cfg%Lcllatlid = Lcllatlid
   cfg%Lclmatlid = Lclmatlid
   cfg%Lcltatlid = Lcltatlid

! MISR simulator   (3 diagnostics, 1 clef)
  cfg%LclMISR = LclMISR

! RTTOV  (1 diagnostic)
  cfg%Ltbrttov = Ltbrttov

! MODIS simulator   (20 diagnostics, 18 clefs)
  cfg%Lcltmodis=Lcltmodis
  cfg%Lclwmodis=Lclwmodis
  cfg%Lclimodis=Lclimodis
  cfg%Lclhmodis=Lclhmodis
  cfg%Lclmmodis=Lclmmodis
  cfg%Lcllmodis=Lcllmodis
  cfg%Ltautmodis=Ltautmodis
  cfg%Ltauwmodis=Ltauwmodis
  cfg%Ltauimodis=Ltauimodis
  cfg%Ltautlogmodis=Ltautlogmodis
  cfg%Ltauwlogmodis=Ltauwlogmodis
  cfg%Ltauilogmodis=Ltauilogmodis
  cfg%Lreffclwmodis=Lreffclwmodis
  cfg%Lreffclimodis=Lreffclimodis
  cfg%Lpctmodis=Lpctmodis
  cfg%Llwpmodis=Llwpmodis
  cfg%Liwpmodis=Liwpmodis
  cfg%Lclmodis=Lclmodis

! Others
  cfg%Lfracout=Lfracout

! Copy instrument flags to cfg structure
  cfg%Lcloudsat = Lcloudsat
  cfg%Lcalipso = Lcalipso
  cfg%LgrLidar532 = LgrLidar532
  cfg%Latlid = Latlid
  cfg%Lparasol = Lparasol
  cfg%Lisccp = Lisccp
  cfg%Lmodis = Lmodis
  cfg%Lmisr = Lmisr
  cfg%Lrttov = Lrttov
  if (cfg%Lcloudsat .or. cfg%Lcalipso .or. cfg%LgrLidar532 .or. cfg%Latlid .or. cfg%Lisccp .or. cfg%Lmodis .or. cfg%Lmisr .or. cfg%Lrttov .or. cfg%Lparasol .or. cfg%Lfracout) then
   cfg%Lwrite_output=.TRUE.
 endif 

#endif

  END SUBROUTINE read_xiosfieldactive

END MODULE lmdz_cosp_read_outputkeys
