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
!!!=============================================================================

module cosp_read_otputkeys

!  USE MOD_COSP_CONSTANTS
!  USE MOD_COSP_TYPES
  use MOD_COSP_INTERFACE_v1p4
  USE mod_phys_lmdz_para

CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE READ_COSP_OUTPUT_NL -------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SUBROUTINE cosp_outputkeys_init(cfg)
  implicit none
  type(cosp_config),intent(out) :: cfg 
  character(len=32) :: out_list(78)
  integer :: i

                
   do i=1,78
      cfg%out_list(i)=''
   enddo

   cfg%Llidar_sim=.false.
   cfg%Lradar_sim=.false.
   cfg%Lparasol_sim=.false.
   cfg%Lisccp_sim=.false.
   cfg%Lmodis_sim=.false.
   cfg%Lmisr_sim=.false.
   cfg%Lrttov_sim=.false.
   cfg%Lstats=.false.
   cfg%Lwrite_output=.false.
   cfg%Ltoffset=.false.
   cfg%Lfracout=.false.

  cfg%Lcllcalipso=.FALSE.
  cfg%Lclmcalipso=.FALSE.
  cfg%Lclhcalipso=.FALSE.
  cfg%Lcltcalipso=.FALSE.
  cfg%Lcllcalipsoice=.FALSE.
  cfg%Lclmcalipsoice=.FALSE.
  cfg%Lclhcalipsoice=.FALSE.
  cfg%Lcltcalipsoice=.FALSE.
  cfg%Lcllcalipsoliq=.FALSE.
  cfg%Lclmcalipsoliq=.FALSE.
  cfg%Lclhcalipsoliq=.FALSE.
  cfg%Lcltcalipsoliq=.FALSE.
  cfg%Lcllcalipsoun=.FALSE.
  cfg%Lclmcalipsoun=.FALSE.
  cfg%Lclhcalipsoun=.FALSE.
  cfg%Lcltcalipsoun=.FALSE.
  cfg%Lclcalipso=.FALSE.
  cfg%Lclcalipsoice=.FALSE.
  cfg%Lclcalipsoliq=.FALSE.
  cfg%Lclcalipsoun=.FALSE.
  cfg%Lclcalipsotmp=.FALSE.
  cfg%Lclcalipsotmpice=.FALSE.
  cfg%Lclcalipsotmpliq=.FALSE.
  cfg%Lclcalipsotmpun=.FALSE.
  cfg%LparasolRefl=.FALSE.
  cfg%LcfadLidarsr532=.FALSE.
  cfg%Latb532=.FALSE.
  cfg%LlidarBetaMol532=.FALSE.
!  cfg%Lclopaquecalipso=.FALSE.
!  cfg%Lclthincalipso=.FALSE.
!  cfg%Lclzopaquecalipso=.FALSE.
!  cfg%Lclcalipsoopaque=.FALSE.
!  cfg%Lclcalipsothin=.FALSE.
!  cfg%Lclcalipsozopaque=.FALSE.
!  cfg%Lclcalipsoopacity=.FALSE.
!  cfg%Lproftemp=.FALSE.
!  cfg%LprofSR=.FALSE.

  cfg%LcfadDbze94=.FALSE.
  cfg%Ldbze94=.FALSE.
  cfg%Lcltlidarradar=.FALSE.
  cfg%Lclcalipso2=.FALSE.

  cfg%Lclisccp=.FALSE.
  cfg%Lboxtauisccp=.FALSE.
  cfg%Lboxptopisccp=.FALSE.
  cfg%Lcltisccp=.FALSE.
  cfg%Lpctisccp=.FALSE.
  cfg%Ltauisccp=.FALSE.
  cfg%Lalbisccp=.FALSE.
  cfg%Lmeantbisccp=.FALSE.
  cfg%Lmeantbclrisccp=.FALSE.

  cfg%LclMISR=.FALSE.

  cfg%Lcllmodis=.FALSE.
  cfg%Lclmmodis=.FALSE.
  cfg%Lclhmodis=.FALSE.
  cfg%Lcltmodis=.FALSE.
  cfg%Lclwmodis=.FALSE.
  cfg%Lclimodis=.FALSE.
  cfg%Ltautmodis=.FALSE.
  cfg%Ltauwmodis=.FALSE.
  cfg%Ltauimodis=.FALSE.
  cfg%Ltautlogmodis=.FALSE.
  cfg%Ltauilogmodis=.FALSE.
  cfg%Ltauwlogmodis=.FALSE.
  cfg%Lreffclwmodis=.FALSE.
  cfg%Lreffclimodis=.FALSE.
  cfg%Lpctmodis=.FALSE.
  cfg%Llwpmodis=.FALSE.
  cfg%Liwpmodis=.FALSE.
  cfg%Lclmodis=.FALSE.
!  cfg%Lcrimodis=.FALSE.
!  cfg%Lcrlmodis=.FALSE.

  cfg%Ltbrttov=.FALSE.

 end subroutine cosp_outputkeys_init

 SUBROUTINE cosp_outputkeys_test(cfg)
  implicit none
  type(cosp_config),intent(out) :: cfg
  character(len=32) :: out_list(78)
  integer :: i


   do i=1,78
      cfg%out_list(i)=''
   enddo

   cfg%Llidar_sim=.true.
   cfg%Lradar_sim=.false.
   cfg%Lisccp_sim=.false.
   cfg%Lmodis_sim=.false.
   cfg%Lmisr_sim=.false.
   cfg%Lrttov_sim=.false.
   cfg%Lstats=.false.
   cfg%Lwrite_output=.false.
   cfg%Ltoffset=.false.
   cfg%Lfracout=.false.

  cfg%Lcllcalipso=.TRUE.
  cfg%Lclmcalipso=.TRUE.
  cfg%Lclhcalipso=.TRUE.
  cfg%Lcltcalipso=.TRUE.
  cfg%Lcllcalipsoice=.FALSE.
  cfg%Lclmcalipsoice=.FALSE.
  cfg%Lclhcalipsoice=.FALSE.
  cfg%Lcltcalipsoice=.FALSE.
  cfg%Lcllcalipsoliq=.FALSE.
  cfg%Lclmcalipsoliq=.FALSE.
  cfg%Lclhcalipsoliq=.FALSE.
  cfg%Lcltcalipsoliq=.FALSE.
  cfg%Lcllcalipsoun=.FALSE.
  cfg%Lclmcalipsoun=.FALSE.
  cfg%Lclhcalipsoun=.FALSE.
  cfg%Lcltcalipsoun=.FALSE.
  cfg%Lclcalipso=.FALSE.
  cfg%Lclcalipsoice=.FALSE.
  cfg%Lclcalipsoliq=.FALSE.
  cfg%Lclcalipsoun=.FALSE.
  cfg%Lclcalipsotmp=.FALSE.
  cfg%Lclcalipsotmpice=.FALSE.
  cfg%Lclcalipsotmpliq=.FALSE.
  cfg%Lclcalipsotmpun=.FALSE.
  cfg%LparasolRefl=.FALSE.
  cfg%LcfadLidarsr532=.FALSE.
  cfg%Latb532=.FALSE.
  cfg%LlidarBetaMol532=.FALSE.
!  cfg%Lclopaquecalipso=.FALSE.
!  cfg%Lclthincalipso=.FALSE.
!  cfg%Lclzopaquecalipso=.FALSE.
!  cfg%Lclcalipsoopaque=.FALSE.
!  cfg%Lclcalipsothin=.FALSE.
!  cfg%Lclcalipsozopaque=.FALSE.
!  cfg%Lclcalipsoopacity=.FALSE.
!  cfg%Lproftemp=.FALSE.
!  cfg%LprofSR=.FALSE.

  cfg%LcfadDbze94=.FALSE.
  cfg%Ldbze94=.FALSE.
  cfg%Lcltlidarradar=.FALSE.
  cfg%Lclcalipso2=.FALSE.

  cfg%Lclisccp=.FALSE.
  cfg%Lboxtauisccp=.FALSE.
  cfg%Lboxptopisccp=.FALSE.
  cfg%Lcltisccp=.FALSE.
  cfg%Lpctisccp=.FALSE.
  cfg%Ltauisccp=.FALSE.
  cfg%Lalbisccp=.FALSE.
  cfg%Lmeantbisccp=.FALSE.
  cfg%Lmeantbclrisccp=.FALSE.

  cfg%LclMISR=.FALSE.

  cfg%Lcllmodis=.FALSE.
  cfg%Lclmmodis=.FALSE.
  cfg%Lclhmodis=.FALSE.
  cfg%Lcltmodis=.FALSE.
  cfg%Lclwmodis=.FALSE.
  cfg%Lclimodis=.FALSE.
  cfg%Ltautmodis=.FALSE.
  cfg%Ltauwmodis=.FALSE.
  cfg%Ltauimodis=.FALSE.
  cfg%Ltautlogmodis=.FALSE.
  cfg%Ltauilogmodis=.FALSE.
  cfg%Ltauwlogmodis=.FALSE.
  cfg%Lreffclwmodis=.FALSE.
  cfg%Lreffclimodis=.FALSE.
  cfg%Lpctmodis=.FALSE.
  cfg%Llwpmodis=.FALSE.
  cfg%Liwpmodis=.FALSE.
  cfg%Lclmodis=.FALSE.
!  cfg%Lcrimodis=.FALSE.
!  cfg%Lcrlmodis=.FALSE.

  cfg%Ltbrttov=.FALSE.

 end subroutine cosp_outputkeys_test

 SUBROUTINE READ_COSP_OUTPUT_NL(itap,cosp_nl,cfg)

#ifdef CPP_XIOS
    USE xios, ONLY: xios_field_is_active
#endif
  implicit none
  character(len=*),intent(in) :: cosp_nl
  type(cosp_config),intent(out) :: cfg
  ! Local variables
  integer :: i, itap

 logical, save :: Lradar_sim,Llidar_sim,Lparasol_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim, Lstats, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
             LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp,Lcllcalipso, &
             Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp,Lcltisccp, &
             Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun, &
             Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun, &
             Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun, &
             Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun, &
             Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun, &
             Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun, &
             Ltoffset,LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfracout,LlidarBetaMol532,Ltbrttov, &
             Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
             Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
             Liwpmodis,Lclmodis,Lcrimodis,Lcrlmodis,Lclopaquecalipso,Lclthincalipso,      &           !OPAQ (2)
             Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,Lclcalipsoopacity, & !OPAQ (5)
             LprofSR,Lproftemp                                                                        !TIBO (2)

  namelist/COSP_OUTPUT/Lradar_sim,Llidar_sim,Lparasol_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
             LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp, &
             Lcllcalipso,Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp, &
             Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun, &
             Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun, &
             Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun, &
             Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun, &
             Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun, &
             Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun, &
             Lcltisccp,Ltoffset,LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfracout,LlidarBetaMol532,Ltbrttov, &
             Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
             Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
             Liwpmodis,Lclmodis,Lcrimodis,Lcrlmodis,Lclopaquecalipso,Lclthincalipso,      &           !OPAQ (2)
             Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,Lclcalipsoopacity, & !OPAQ (5)
             LprofSR,Lproftemp                                                                        !TIBO (2)
   
  do i=1,78
    cfg%out_list(i)=''
  enddo
  
! Lecture du fichier namelist
  IF (is_master) THEN
    open(10,file=cosp_nl,status='old')
    read(10,nml=cosp_output)
    close(10)
  ENDIF
  
  CALL bcast(Lradar_sim)
  CALL bcast(Lparasol_sim)
  CALL bcast(Llidar_sim)
  CALL bcast(Lisccp_sim)
  CALL bcast(Lmodis_sim)
  CALL bcast(Lmisr_sim)
  CALL bcast(Lrttov_sim)

  CALL bcast(Lstats)

  CALL bcast(Lalbisccp)
  CALL bcast(Latb532)
  CALL bcast(Lboxptopisccp)
  CALL bcast(Lboxtauisccp)
  CALL bcast(LcfadDbze94)
  CALL bcast(LcfadLidarsr532)
  CALL bcast(Lclcalipso2)
  CALL bcast(Lclcalipso)
  CALL bcast(Lclhcalipso)
  CALL bcast(Lclcalipsoliq)
  CALL bcast(Lclcalipsoice)
  CALL bcast(Lclcalipsoun)
  CALL bcast(Lclcalipsotmp)
  CALL bcast(Lclcalipsotmpliq)
  CALL bcast(Lclcalipsotmpice)
  CALL bcast(Lclcalipsotmpun)
  CALL bcast(Lcltcalipsoliq)
  CALL bcast(Lcltcalipsoice)
  CALL bcast(Lcltcalipsoun)
  CALL bcast(Lclhcalipsoliq)
  CALL bcast(Lclhcalipsoice)
  CALL bcast(Lclhcalipsoun)
  CALL bcast(Lclmcalipsoliq)
  CALL bcast(Lclmcalipsoice)
  CALL bcast(Lclmcalipsoun)
  CALL bcast(Lcllcalipsoliq)
  CALL bcast(Lcllcalipsoice) 
  CALL bcast(Lcllcalipsoun)
  CALL bcast(Lclisccp)
  CALL bcast(Lcllcalipso)
  CALL bcast(Lclmcalipso)
  CALL bcast(Lcltcalipso)
  CALL bcast(Lcltlidarradar)
  CALL bcast(Lpctisccp)
  CALL bcast(Ldbze94)
  CALL bcast(Ltauisccp)
  CALL bcast(Lcltisccp)
  CALL bcast(LparasolRefl)
  CALL bcast(LclMISR)
  CALL bcast(Lmeantbisccp)
  CALL bcast(Lmeantbclrisccp)
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
  CALL bcast(Lcrimodis)
  CALL bcast(Lcrlmodis)
  CALL bcast(Lclopaquecalipso)  !OPAQ
  CALL bcast(Lclthincalipso)    !OPAQ
  CALL bcast(Lclzopaquecalipso) !OPAQ
  CALL bcast(Lclcalipsoopaque)  !OPAQ
  CALL bcast(Lclcalipsothin)    !OPAQ
  CALL bcast(Lclcalipsozopaque) !OPAQ
  CALL bcast(Lclcalipsoopacity) !OPAQ
  CALL bcast(LprofSR)           !TIBO
  CALL bcast(Lproftemp)         !TIBO

!$OMP BARRIER

!  print*,' Cles sorties cosp :'
!  print*,' Lradar_sim,Llidar_sim,Lisccp_sim,Lmisr_sim,Lrttov_sim', &
!           Lradar_sim,Llidar_sim,Lisccp_sim,Lmisr_sim,Lrttov_sim


  ! Deal with dependencies
  if (.not.Lradar_sim) then
    LcfadDbze94   = .false.
    Lclcalipso2    = .false.
    Lcltlidarradar = .false. ! Needs radar & lidar
    Ldbze94        = .false.
    Lclcalipso2    = .false. ! Needs radar & lidar
  endif

  if (.not.Lparasol_sim) then
   LparasolRefl   = .false.
  endif 

  if (.not.Llidar_sim) then
    Latb532          = .false.
    LcfadLidarsr532  = .false.
    Lclcalipso2      = .false.
    Lclcalipso       = .false.
    Lclhcalipso      = .false.
    Lcllcalipso      = .false.
    Lclmcalipso      = .false.
    Lcltcalipso      = .false.
    Lcltlidarradar   = .false. ! Needs radar & lidar
    LparasolRefl     = .false.
    LlidarBetaMol532 = .false.
!! AI
    Lclcalipsoliq       = .false.
    Lclcalipsoice       = .false.
    Lclcalipsoun        = .false.
    Lclcalipsotmp       = .false.
    Lclcalipsotmpun     = .false.
    Lclcalipsotmpliq    = .false.
    Lclcalipsotmpice    = .false.
    Lclhcalipsoliq      = .false.
    Lcllcalipsoliq      = .false.
    Lclmcalipsoliq      = .false.
    Lcltcalipsoliq      = .false.
    Lclhcalipsoice      = .false.
    Lcllcalipsoice      = .false.
    Lclmcalipsoice      = .false.
    Lcltcalipsoice      = .false.
    Lclhcalipsoun       = .false.
    Lcllcalipsoun       = .false.
    Lclmcalipsoun       = .false.
    Lcltcalipsoun       = .false.
    Lclopaquecalipso    = .false. !OPAQ
    Lclthincalipso      = .false. !OPAQ
    Lclzopaquecalipso   = .false. !OPAQ
    Lclcalipsoopaque    = .false. !OPAQ
    Lclcalipsothin      = .false. !OPAQ
    Lclcalipsozopaque   = .false. !OPAQ
    Lclcalipsoopacity   = .false. !OPAQ
    LprofSR             = .false. !TIBO
    Lproftemp           = .false. !TIBO
  endif

  if (.not.Lisccp_sim) then
    Lalbisccp       = .false.
    Lboxptopisccp   = .false.
    Lboxtauisccp    = .false.
    Lclisccp        = .false.
    Lpctisccp       = .false.
    Ltauisccp       = .false.
    Lcltisccp       = .false.
    Lmeantbisccp    = .false.
    Lmeantbclrisccp = .false.
  endif

  if (.not.Lmisr_sim) then
    LclMISR = .false.
  endif
  if (.not.Lrttov_sim) then
    Ltbrttov = .false.
  endif
  if ((.not.Lradar_sim).and.(.not.Llidar_sim).and. &
      (.not.Lisccp_sim).and.(.not.Lmisr_sim)) then
    Lfracout = .false.
    Lstats = .false.
  endif
 if (.not.Lmodis_sim) then
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
    Lcrimodis=.false.
    Lcrlmodis=.false.
  endif
  if (Lmodis_sim) Lisccp_sim = .true.

  ! Diagnostics that use Radar and Lidar
  if (((Lclcalipso2).or.(Lcltlidarradar)).and.((Lradar_sim).or.(Llidar_sim))) then
    Lclcalipso2    = .true.
    Lcltlidarradar = .true.
    Llidar_sim     = .true.
    Lradar_sim     = .true.
  endif

  if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) Lstats = .true.

  ! Copy instrument flags to cfg structure
  cfg%Lradar_sim = Lradar_sim
  cfg%Llidar_sim = Llidar_sim
  cfg%Lisccp_sim = Lisccp_sim
  cfg%Lmodis_sim = Lmodis_sim
  cfg%Lmisr_sim  = Lmisr_sim
  cfg%Lrttov_sim = Lrttov_sim
  cfg%Lparasol_sim = Lparasol_sim

  cfg%Lstats = Lstats

  ! Flag to control output to file
  cfg%Lwrite_output = .false.
  if (cfg%Lstats.or.cfg%Lmisr_sim.or.cfg%Lrttov_sim) then
    cfg%Lwrite_output = .true.
  endif

  ! Output diagnostics
  i = 1
  if (Lalbisccp)        cfg%out_list(i) = 'albisccp'
  i = i+1
  if (Latb532)          cfg%out_list(i) = 'atb532'
  i = i+1
  if (Lboxptopisccp)    cfg%out_list(i) = 'boxptopisccp'
  i = i+1
  if (Lboxtauisccp)     cfg%out_list(i) = 'boxtauisccp'
  i = i+1
  if (LcfadDbze94)      cfg%out_list(i) = 'cfadDbze94'
  i = i+1
  if (LcfadLidarsr532)  cfg%out_list(i) = 'cfadLidarsr532'
  i = i+1
  if (Lclcalipso2)      cfg%out_list(i) = 'clcalipso2'
  i = i+1
  if (Lclcalipso)       cfg%out_list(i) = 'clcalipso'
  i = i+1
  if (Lclhcalipso)      cfg%out_list(i) = 'clhcalipso'
  i = i+1
  if (Lclisccp)         cfg%out_list(i) = 'clisccp'
  i = i+1
  if (Lcllcalipso)      cfg%out_list(i) = 'cllcalipso'
  i = i+1
  if (Lclmcalipso)      cfg%out_list(i) = 'clmcalipso'
  i = i+1
  if (Lcltcalipso)      cfg%out_list(i) = 'cltcalipso'
  i = i+1

  if (Lcllcalipsoice)      cfg%out_list(i) = 'cllcalipsoice'
  i = i+1
  if (Lclmcalipsoice)      cfg%out_list(i) = 'clmcalipsoice'
  i = i+1
  if (Lclhcalipsoice)      cfg%out_list(i) = 'clhcalipsoice'
  i = i+1
  if (Lcltcalipsoice)      cfg%out_list(i) = 'cltcalipsoice'
  i = i+1
  if (Lcllcalipsoliq)      cfg%out_list(i) = 'cllcalipsoliq'
  i = i+1
  if (Lclmcalipsoliq)      cfg%out_list(i) = 'clmcalipsoliq'
  i = i+1
  if (Lclhcalipsoliq)      cfg%out_list(i) = 'clhcalipsoliq'
  i = i+1
  if (Lcltcalipsoliq)      cfg%out_list(i) = 'cltcalipsoliq'
  i = i+1
  if (Lcllcalipsoun)      cfg%out_list(i) = 'cllcalipsoun'
  i = i+1
  if (Lclmcalipsoun)      cfg%out_list(i) = 'clmcalipsoun'
  i = i+1
  if (Lclhcalipsoun)      cfg%out_list(i) = 'clhcalipsoun'
  i = i+1
  if (Lcltcalipsoun)      cfg%out_list(i) = 'cltcalipsoun'
  i = i+1

  if (Lclcalipsoice)       cfg%out_list(i) = 'clcalipsoice'
  i = i+1
  if (Lclcalipsoliq)       cfg%out_list(i) = 'clcalipsoliq'
  i = i+1
  if (Lclcalipsoun)       cfg%out_list(i) = 'clcalipsoun'
  i = i+1

  if (Lclcalipsotmp)       cfg%out_list(i) = 'clcalipsotmp'
  i = i+1
  if (Lclcalipsotmpice)       cfg%out_list(i) = 'clcalipsotmpice'
  i = i+1
  if (Lclcalipsotmpliq)       cfg%out_list(i) = 'clcalipsotmpliq'
  i = i+1
  if (Lclcalipsotmpun)       cfg%out_list(i) = 'clcalipsotmpun'
  i = i+1
  if (Lcltlidarradar)   cfg%out_list(i) = 'cltlidarradar'
  i = i+1
  if (Lpctisccp)        cfg%out_list(i) = 'pctisccp'
  i = i+1
  if (Ldbze94)          cfg%out_list(i) = 'dbze94'
  i = i+1
  if (Ltauisccp)        cfg%out_list(i) = 'tauisccp'
  i = i+1
  if (Lcltisccp)        cfg%out_list(i) = 'cltisccp'
  i = i+1
  if (Ltoffset)         cfg%out_list(i) = 'toffset'
  i = i+1
  if (LparasolRefl)     cfg%out_list(i) = 'parasolRefl'
  i = i+1
  if (LclMISR)          cfg%out_list(i) = 'clMISR'
  i = i+1
  if (Lmeantbisccp)     cfg%out_list(i) = 'meantbisccp'
  i = i+1
  if (Lmeantbclrisccp)  cfg%out_list(i) = 'meantbclrisccp'
  i = i+1
  if (Lfracout)         cfg%out_list(i) = 'fracout'
  i = i+1
  if (LlidarBetaMol532) cfg%out_list(i) = 'lidarBetaMol532'
  i = i+1
  if (Ltbrttov)         cfg%out_list(i) = 'tbrttov'
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
  if (Ltautmodis)       cfg%out_list(i) = 'tautmodis'
  i = i+1
  if (Ltauwmodis)       cfg%out_list(i) = 'tauwmodis'
  i = i+1
  if (Ltauimodis)       cfg%out_list(i) = 'tauimodis'
  i = i+1
  if (Ltautlogmodis)    cfg%out_list(i) = 'tautlogmodis'
  i = i+1
  if (Ltauwlogmodis)    cfg%out_list(i) = 'tauwlogmodis'
  i = i+1
  if (Ltauilogmodis)    cfg%out_list(i) = 'tauilogmodis'
  i = i+1
  if (Lreffclwmodis)    cfg%out_list(i) = 'reffclwmodis'
  i = i+1
  if (Lreffclimodis)    cfg%out_list(i) = 'reffclimodis'
  i = i+1
  if (Lpctmodis)        cfg%out_list(i) = 'pctmodis'
  i = i+1
  if (Llwpmodis)        cfg%out_list(i) = 'lwpmodis'
  i = i+1
  if (Liwpmodis)        cfg%out_list(i) = 'iwpmodis'
  i = i+1
  if (Lclmodis)         cfg%out_list(i) = 'clmodis'
  i = i+1
  if (Lcrimodis)         cfg%out_list(i) = 'crimodis'
  i = i+1
  if (Lcrlmodis)         cfg%out_list(i) = 'crlmodis'

  i = i+1                                                            !OPAQ
  if (Lclopaquecalipso)         cfg%out_list(i) = 'clopaquecalipso'  !OPAQ
  i = i+1                                                            !OPAQ
  if (Lclthincalipso)           cfg%out_list(i) = 'clthincalipso'    !OPAQ
  i = i+1                                                            !OPAQ
  if (Lclzopaquecalipso)        cfg%out_list(i) = 'clzopaquecalipso' !OPAQ
  i = i+1                                                            !OPAQ
  if (Lclcalipsoopaque)         cfg%out_list(i) = 'clcalipsoopaque'  !OPAQ
  i = i+1                                                            !OPAQ
  if (Lclcalipsothin)           cfg%out_list(i) = 'clcalipsothin'    !OPAQ
  i = i+1                                                            !OPAQ
  if (Lclcalipsozopaque)        cfg%out_list(i) = 'clcalipsozopaque' !OPAQ
  i = i+1                                                            !OPAQ
  if (Lclcalipsoopacity)        cfg%out_list(i) = 'clcalipsoopacity' !OPAQ
  i = i+1                                                            !TIBO
  if (LprofSR)                  cfg%out_list(i) = 'profSR'           !TIBO
  i = i+1                                                            !TIBO
  if (Lproftemp)                cfg%out_list(i) = 'proftemp'         !TIBO
    
  if (i.gt.78) then
     print *, 'COSP_IO: wrong number of output diagnostics'
     print *, i,78
     stop
  endif

  ! Copy diagnostic flags to cfg structure
  ! ISCCP simulator  
  cfg%Lalbisccp = Lalbisccp
  cfg%Latb532 = Latb532
  cfg%Lboxptopisccp = Lboxptopisccp
  cfg%Lboxtauisccp = Lboxtauisccp
  cfg%Lmeantbisccp = Lmeantbisccp
  cfg%Lmeantbclrisccp = Lmeantbclrisccp
  cfg%Lclisccp = Lclisccp
  cfg%Lpctisccp = Lpctisccp
  cfg%Ltauisccp = Ltauisccp
  cfg%Lcltisccp = Lcltisccp
  ! CloudSat simulator  
  cfg%Ldbze94 = Ldbze94
  cfg%LcfadDbze94 = LcfadDbze94
  ! CALIPSO/PARASOL simulator  
  cfg%LcfadLidarsr532 = LcfadLidarsr532
  cfg%Lclcalipso2 = Lclcalipso2
  cfg%Lclcalipso = Lclcalipso
  cfg%Lclhcalipso = Lclhcalipso
  cfg%Lcllcalipso = Lcllcalipso
  cfg%Lclmcalipso = Lclmcalipso
  cfg%Lcltcalipso = Lcltcalipso
  cfg%Lclhcalipsoice = Lclhcalipsoice
  cfg%Lcllcalipsoice = Lcllcalipsoice
  cfg%Lclmcalipsoice = Lclmcalipsoice
  cfg%Lcltcalipsoice = Lcltcalipsoice
  cfg%Lclhcalipsoliq = Lclhcalipsoliq
  cfg%Lcllcalipsoliq = Lcllcalipsoliq
  cfg%Lclmcalipsoliq = Lclmcalipsoliq
  cfg%Lcltcalipsoliq = Lcltcalipsoliq
  cfg%Lclhcalipsoun = Lclhcalipsoun
  cfg%Lcllcalipsoun = Lcllcalipsoun
  cfg%Lclmcalipsoun = Lclmcalipsoun
  cfg%Lcltcalipsoun = Lcltcalipsoun
  cfg%Lclcalipsoice = Lclcalipsoice
  cfg%Lclcalipsoliq = Lclcalipsoliq
  cfg%Lclcalipsoun = Lclcalipsoun
  cfg%Lclcalipsotmp = Lclcalipsotmp
  cfg%Lclcalipsotmpice = Lclcalipsotmpice
  cfg%Lclcalipsotmpliq = Lclcalipsotmpliq
  cfg%Lclcalipsotmpun = Lclcalipsotmpun
  cfg%Lcltlidarradar = Lcltlidarradar
  cfg%LparasolRefl = LparasolRefl
!  cfg%Lclopaquecalipso  = Lclopaquecalipso  !OPAQ
!  cfg%Lclthincalipso    = Lclthincalipso    !OPAQ
!  cfg%Lclzopaquecalipso = Lclzopaquecalipso !OPAQ
!  cfg%Lclcalipsoopaque  = Lclcalipsoopaque  !OPAQ
!  cfg%Lclcalipsothin    = Lclcalipsothin    !OPAQ
!  cfg%Lclcalipsozopaque = Lclcalipsozopaque !OPAQ
!  cfg%Lclcalipsoopacity = Lclcalipsoopacity !OPAQ
!  cfg%LprofSR           = LprofSR           !TIBO
!  cfg%Lproftemp         = Lproftemp         !TIBO
  ! MISR simulator  
  cfg%LclMISR = LclMISR
  ! Other
  cfg%Ltoffset = Ltoffset
  cfg%Lfracout = Lfracout
  cfg%LlidarBetaMol532 = LlidarBetaMol532
  ! RTTOV
  cfg%Ltbrttov = Ltbrttov
  ! MODIS simulator  
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
!  cfg%Lcrimodis=Lcrimodis
!  cfg%Lcrlmodis=Lcrlmodis
  
 END SUBROUTINE READ_COSP_OUTPUT_NL

 SUBROUTINE read_xiosfieldactive(cfg)

!    USE MOD_COSP_CONSTANTS
!    USE MOD_COSP_TYPES
#ifdef CPP_XIOS
    USE xios, ONLY: xios_field_is_active
#endif
  implicit none
  type(cosp_config),intent(out) :: cfg
  integer :: i

#ifdef CPP_XIOS

 logical :: Lradar_sim,Llidar_sim,Lparasol_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim, Lstats, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
             LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp,Lcllcalipso, &
             Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp,Lcltisccp, & 
             Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun, &
             Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun, &
             Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun, &
             Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun, &
             Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun, &
             Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun, &
             Ltoffset,LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfracout,LlidarBetaMol532,Ltbrttov, &
             Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis, &
             Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
             Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, & 
             Liwpmodis,Lclmodis,Lcrimodis,Lcrlmodis,Lclopaquecalipso,Lclthincalipso, &
             Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,Lclcalipsoopacity, &
             LprofSR,Lproftemp
        
  character(len=32) :: out_list(78)

  do i=1,78
    cfg%out_list(i)=''
  enddo

    LcfadDbze94   = .false.
    Lclcalipso2    = .false.
    Lcltlidarradar = .false. ! Needs radar & lidar
    Ldbze94        = .false.
    Lclcalipso2    = .false. ! Needs radar & lidar

    Latb532          = .false.
    LcfadLidarsr532  = .false.
    Lclcalipso       = .false.
    Lclhcalipso      = .false.
    Lcllcalipso      = .false.
    Lclmcalipso      = .false.
    Lcltcalipso      = .false.
    LparasolRefl     = .false.
    LlidarBetaMol532 = .false.
    Lclcalipsoliq       = .false.
    Lclcalipsoice       = .false.
    Lclcalipsoun        = .false.
    Lclcalipsotmp       = .false.
    Lclcalipsotmpun     = .false.
    Lclcalipsotmpliq    = .false.
    Lclcalipsotmpice    = .false.
    Lclhcalipsoliq      = .false.
    Lcllcalipsoliq      = .false.
    Lclmcalipsoliq      = .false.
    Lcltcalipsoliq      = .false.
    Lclhcalipsoice      = .false.
    Lcllcalipsoice      = .false.
    Lclmcalipsoice      = .false.
    Lcltcalipsoice      = .false.
    Lclhcalipsoun       = .false.
    Lcllcalipsoun       = .false.
    Lclmcalipsoun       = .false.
    Lcltcalipsoun       = .false.
    Lclopaquecalipso    = .false. !OPAQ
    Lclthincalipso      = .false. !OPAQ
    Lclzopaquecalipso   = .false. !OPAQ
    Lclcalipsoopaque    = .false. !OPAQ
    Lclcalipsothin      = .false. !OPAQ
    Lclcalipsozopaque   = .false. !OPAQ
    Lclcalipsoopacity   = .false. !OPAQ
    LprofSR             = .false. !TIBO
    Lproftemp           = .false. !TIBO

    Lalbisccp       = .false.
    Lboxptopisccp   = .false.
    Lboxtauisccp    = .false.
    Lclisccp        = .false.
    Lpctisccp       = .false.
    Ltauisccp       = .false.
    Lcltisccp       = .false.
    Lmeantbisccp    = .false.
    Lmeantbclrisccp = .false.

    LclMISR = .false.

    Ltbrttov = .false.

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
    Lcrimodis=.false.
    Lcrlmodis=.false.

    Lradar_sim=.false.
    Llidar_sim=.false.
    Lparasol_sim=.false.
    Lisccp_sim=.false.
    Lmodis_sim=.false.
    Lmisr_sim=.false.
    Lrttov_sim=.false.

    Lstats=.false.
!    Ltoffset=.false.
!    Lfracout=.false.
!    Lwrite_output=.false.

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
 IF (xios_field_is_active("parasol_refl")) LparasolRefl=.TRUE.
! IF (xios_field_is_active("parasol_crefl")) cfg%LparasolRefl=.TRUE.
! IF (xios_field_is_active("Ncrefl")) cfg%LparasolRefl=.TRUE.
 IF (xios_field_is_active("cfad_lidarsr532")) LcfadLidarsr532=.TRUE.
 IF (xios_field_is_active("atb532")) Latb532=.TRUE.
 IF (xios_field_is_active("beta_mol532")) LlidarBetaMol532=.TRUE.
 IF (xios_field_is_active("clopaquecalipso")) Lclopaquecalipso=.TRUE.
 IF (xios_field_is_active("clthincalipso")) Lclthincalipso=.TRUE.
 IF (xios_field_is_active("clzopaquecalipso")) Lclzopaquecalipso=.TRUE.
 IF (xios_field_is_active("clcalipsoopaque")) Lclcalipsoopaque=.TRUE.
 IF (xios_field_is_active("clcalipsothin")) Lclcalipsothin=.TRUE.
 IF (xios_field_is_active("clcalipsozopaque")) Lclcalipsozopaque=.TRUE.
 IF (xios_field_is_active("clcalipsoopacity")) Lclcalipsoopacity=.TRUE.
 IF (xios_field_is_active("proftemp")) Lproftemp=.TRUE.
 IF (xios_field_is_active("profSR")) LprofSR=.TRUE.
!!!! 38 champ Calipso

 IF (xios_field_is_active("cfadDbze94")) LcfadDbze94=.TRUE.
 IF (xios_field_is_active("dbze94")) Ldbze94=.TRUE.
!!! 2 champs CLOUDSAT

 IF (xios_field_is_active("cltlidarradar")) Lcltlidarradar=.TRUE.
 IF (xios_field_is_active("clcalipso2")) Lclcalipso2=.TRUE.
!!! 2 champs CLOUDSAT et CALIPSO

 IF (xios_field_is_active("clisccp2")) Lclisccp=.TRUE.
 IF (xios_field_is_active("boxtauisccp")) Lboxtauisccp=.TRUE.
 IF (xios_field_is_active("boxptopisccp")) Lboxptopisccp=.TRUE.
 IF (xios_field_is_active("tclisccp")) Lcltisccp=.TRUE.
 IF (xios_field_is_active("ctpisccp")) Lpctisccp=.TRUE.
 IF (xios_field_is_active("tauisccp")) Ltauisccp=.TRUE.
 IF (xios_field_is_active("albisccp")) Lalbisccp=.TRUE.
 IF (xios_field_is_active("meantbisccp")) Lmeantbisccp=.TRUE.
 IF (xios_field_is_active("meantbclrisccp")) Lmeantbclrisccp=.TRUE.
!!! 9 champs ISCCP

 IF (xios_field_is_active("clMISR")) LclMISR=.TRUE.
!!! 1 champs MISR

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
 IF (xios_field_is_active("clmodis")) Lclmodis=.TRUE.
 IF (xios_field_is_active("crimodis")) Lcrimodis=.TRUE.
 IF (xios_field_is_active("crlmodis")) Lcrlmodis=.TRUE.
!!! 20 champs MODIS
! IF (xios_field_is_active("tbrttov")) cfg%Ltbrttov=.TRUE.

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
       xios_field_is_active("cfad_lidarsr532").OR. &
       xios_field_is_active("atb532").OR. &
       xios_field_is_active("beta_mol532").OR. &
       xios_field_is_active("clopaquecalipso").OR. &
       xios_field_is_active("clthincalipso").OR. &
       xios_field_is_active("clzopaquecalipso").OR. &
       xios_field_is_active("clcalipsoopaque").OR. &
       xios_field_is_active("clcalipsothin").OR. &
       xios_field_is_active("clcalipsozopaque").OR. &
       xios_field_is_active("clcalipsoopacity").OR. &
       xios_field_is_active("proftemp").OR. &
       xios_field_is_active("profSR")) Llidar_sim=.TRUE.

    IF (xios_field_is_active("parasol_refl")) Lparasol_sim=.TRUE.

    IF (xios_field_is_active("cfadDbze94").OR. &
      xios_field_is_active("dbze94")) & Lradar_sim=.TRUE.

    IF (xios_field_is_active("cltlidarradar").OR. &
      xios_field_is_active("clcalipso2")) THEN
               Lradar_sim=.TRUE.
               Llidar_sim=.TRUE.
    ENDIF

    IF (xios_field_is_active("clisccp2").OR. &
       xios_field_is_active("boxtauisccp").OR. &
       xios_field_is_active("boxptopisccp").OR. &
       xios_field_is_active("tclisccp").OR. &
       xios_field_is_active("ctpisccp").OR. &
       xios_field_is_active("tauisccp").OR. &
       xios_field_is_active("albisccp").OR. &
       xios_field_is_active("meantbisccp").OR. &
       xios_field_is_active("meantbclrisccp")) Lisccp_sim=.TRUE.

    IF (xios_field_is_active("clMISR")) Lmisr_sim=.TRUE.

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
       xios_field_is_active("clmodis").OR. &
       xios_field_is_active("crimodis").OR. &
       xios_field_is_active("crlmodis")) Lmodis_sim=.TRUE.

  ENDIF !   (is_master) 

!$OMP BARRIER

  CALL bcast(Lradar_sim)
  CALL bcast(Llidar_sim)
  CALL bcast(Lisccp_sim)
  CALL bcast(Lmodis_sim)
  CALL bcast(Lmisr_sim)
  CALL bcast(Lparasol_sim)
  CALL bcast(Lrttov_sim)

  CALL bcast(Lstats)

  CALL bcast(Lalbisccp)
  CALL bcast(Latb532)
  CALL bcast(Lboxptopisccp)
  CALL bcast(Lboxtauisccp)
  CALL bcast(LcfadDbze94)
  CALL bcast(LcfadLidarsr532)
  CALL bcast(Lclcalipso2)
  CALL bcast(Lclcalipso)
  CALL bcast(Lclhcalipso)
  CALL bcast(Lclcalipsoliq)
  CALL bcast(Lclcalipsoice)
  CALL bcast(Lclcalipsoun)
  CALL bcast(Lclcalipsotmp)
  CALL bcast(Lclcalipsotmpliq)
  CALL bcast(Lclcalipsotmpice)
  CALL bcast(Lclcalipsotmpun)
  CALL bcast(Lcltcalipsoliq)
  CALL bcast(Lcltcalipsoice)
  CALL bcast(Lcltcalipsoun)
  CALL bcast(Lclhcalipsoliq)
  CALL bcast(Lclhcalipsoice)
  CALL bcast(Lclhcalipsoun)
  CALL bcast(Lclmcalipsoliq)
  CALL bcast(Lclmcalipsoice)
  CALL bcast(Lclmcalipsoun)
  CALL bcast(Lcllcalipsoliq)
  CALL bcast(Lcllcalipsoice)
  CALL bcast(Lcllcalipsoun)
  CALL bcast(Lclisccp)
  CALL bcast(Lcllcalipso)
  CALL bcast(Lclmcalipso)
  CALL bcast(Lcltcalipso)
  CALL bcast(Lcltlidarradar)
  CALL bcast(Lpctisccp)
  CALL bcast(Ldbze94)
  CALL bcast(Ltauisccp)
  CALL bcast(Lcltisccp)
  CALL bcast(LparasolRefl)
  CALL bcast(LclMISR)
  CALL bcast(Lmeantbisccp)
  CALL bcast(Lmeantbclrisccp)
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
  CALL bcast(Lcrimodis)
  CALL bcast(Lcrlmodis)
  CALL bcast(Lclopaquecalipso)  !OPAQ
  CALL bcast(Lclthincalipso)    !OPAQ
  CALL bcast(Lclzopaquecalipso) !OPAQ
  CALL bcast(Lclcalipsoopaque)  !OPAQ
  CALL bcast(Lclcalipsothin)    !OPAQ
  CALL bcast(Lclcalipsozopaque) !OPAQ
  CALL bcast(Lclcalipsoopacity) !OPAQ
  CALL bcast(LprofSR)           !TIBO
  CALL bcast(Lproftemp)         !TIBO

    if (Lmodis_sim) Lisccp_sim = .true.
    if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) Lstats = .true.
!    IF (xios_field_is_active("tbrttov")) cfg%Lrttov_sim=.TRUE.

  ! Copy diagnostic flags to cfg structure
  ! ISCCP simulator  
  cfg%Lalbisccp = Lalbisccp
  cfg%Latb532 = Latb532
  cfg%Lboxptopisccp = Lboxptopisccp
  cfg%Lboxtauisccp = Lboxtauisccp
  cfg%Lmeantbisccp = Lmeantbisccp
  cfg%Lmeantbclrisccp = Lmeantbclrisccp
  cfg%Lclisccp = Lclisccp
  cfg%Lpctisccp = Lpctisccp
  cfg%Ltauisccp = Ltauisccp
  cfg%Lcltisccp = Lcltisccp

! CloudSat simulator  
  cfg%Ldbze94 = Ldbze94
  cfg%LcfadDbze94 = LcfadDbze94

! Cloudsat et Calipso
  cfg%Lclcalipso2 = Lclcalipso2
  cfg%Lcltlidarradar = Lcltlidarradar

! CALIPSO/PARASOL simulator  
  cfg%LcfadLidarsr532 = LcfadLidarsr532
  cfg%Lclcalipso = Lclcalipso
  cfg%Lclhcalipso = Lclhcalipso
  cfg%Lcllcalipso = Lcllcalipso
  cfg%Lclmcalipso = Lclmcalipso
  cfg%Lcltcalipso = Lcltcalipso
  cfg%Lclhcalipsoice = Lclhcalipsoice
  cfg%Lcllcalipsoice = Lcllcalipsoice
  cfg%Lclmcalipsoice = Lclmcalipsoice
  cfg%Lcltcalipsoice = Lcltcalipsoice
  cfg%Lclhcalipsoliq = Lclhcalipsoliq
  cfg%Lcllcalipsoliq = Lcllcalipsoliq
  cfg%Lclmcalipsoliq = Lclmcalipsoliq
  cfg%Lcltcalipsoliq = Lcltcalipsoliq
  cfg%Lclhcalipsoun = Lclhcalipsoun
  cfg%Lcllcalipsoun = Lcllcalipsoun
  cfg%Lclmcalipsoun = Lclmcalipsoun
  cfg%Lcltcalipsoun = Lcltcalipsoun
  cfg%Lclcalipsoice = Lclcalipsoice
  cfg%Lclcalipsoliq = Lclcalipsoliq
  cfg%Lclcalipsoun = Lclcalipsoun
  cfg%Lclcalipsotmp = Lclcalipsotmp
  cfg%Lclcalipsotmpice = Lclcalipsotmpice
  cfg%Lclcalipsotmpliq = Lclcalipsotmpliq
  cfg%Lclcalipsotmpun = Lclcalipsotmpun
  cfg%LparasolRefl = LparasolRefl
!  cfg%Lclopaquecalipso  = Lclopaquecalipso  !OPAQ
!  cfg%Lclthincalipso    = Lclthincalipso    !OPAQ
!  cfg%Lclzopaquecalipso = Lclzopaquecalipso !OPAQ
!  cfg%Lclcalipsoopaque  = Lclcalipsoopaque  !OPAQ
!  cfg%Lclcalipsothin    = Lclcalipsothin    !OPAQ
!  cfg%Lclcalipsozopaque = Lclcalipsozopaque !OPAQ
!  cfg%Lclcalipsoopacity = Lclcalipsoopacity !OPAQ
!  cfg%LprofSR           = LprofSR           !TIBO
!  cfg%Lproftemp         = Lproftemp         !TIBO
  cfg%LlidarBetaMol532 = LlidarBetaMol532

! MISR simulator  
  cfg%LclMISR = LclMISR

! RTTOV
  cfg%Ltbrttov = Ltbrttov

! MODIS simulator  
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
!  cfg%Lcrimodis=Lcrimodis
!  cfg%Lcrlmodis=Lcrlmodis

! Others
!  cfg%Lwrite_output=Lwrite_output
!  cfg%Ltoffset=Ltoffset
!  cfg%Lfracout=Lfracout
  cfg%Lstats = Lstats

! Copy instrument flags to cfg structure
  cfg%Lradar_sim = Lradar_sim
  cfg%Llidar_sim = Llidar_sim
  cfg%Lparasol_sim = Lparasol_sim
  cfg%Lisccp_sim = Lisccp_sim
  cfg%Lmodis_sim = Lmodis_sim
  cfg%Lmisr_sim  = Lmisr_sim
  cfg%Lrttov_sim = Lrttov_sim
 if (cfg%Lradar_sim.or.cfg%Llidar_sim.or.cfg%Lisccp_sim.or.cfg%Lmodis_sim.or.cfg%Lmisr_sim) then
   cfg%Lwrite_output=.TRUE.
 endif 

#endif

  END SUBROUTINE read_xiosfieldactive

END MODULE cosp_read_otputkeys
