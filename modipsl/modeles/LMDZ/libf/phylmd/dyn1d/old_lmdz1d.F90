!
! $Id: lmdz1d.F90 3540 2019-06-25 14:50:13Z fairhead $
!

      SUBROUTINE old_lmdz1d

   USE ioipsl, only: ju2ymds, ymds2ju, ioconf_calendar,getin
   USE phys_state_var_mod, ONLY : phys_state_var_init, phys_state_var_end, &
       clwcon, detr_therm, &
       qsol, fevap, z0m, z0h, agesno, &
       du_gwd_rando, du_gwd_front, entr_therm, f0, fm_therm, &
       falb_dir, falb_dif, &
       ftsol, pbl_tke, pctsrf, radsol, rain_fall, snow_fall, ratqs, &
       rnebcon, rugoro, sig1, w01, solaire_etat0, sollw, sollwdown, &
       solsw, t_ancien, q_ancien, u_ancien, v_ancien, wake_cstar, &
       wake_delta_pbl_TKE, delta_tsurf, wake_fip, wake_pe, &
       wake_deltaq, wake_deltat, wake_s, wake_dens, &
       zgam, zmax0, zmea, zpic, zsig, &
       zstd, zthe, zval, ale_bl, ale_bl_trig, alp_bl, ql_ancien, qs_ancien, &
       prlw_ancien, prsw_ancien, prw_ancien, &
       u10m,v10m,ale_wake,ale_bl_stat

  
   USE dimphy
   USE surface_data, only : type_ocean,ok_veget
   USE pbl_surface_mod, only : ftsoil, pbl_surface_init, &
                                 pbl_surface_final
   USE fonte_neige_mod, only : fonte_neige_init, fonte_neige_final 

   USE infotrac ! new
   USE control_mod
   USE indice_sol_mod
   USE phyaqua_mod
!  USE mod_1D_cases_read
   USE mod_1D_cases_read2
   USE mod_1D_amma_read
   USE print_control_mod, ONLY: lunout, prt_level
   USE iniphysiq_mod, ONLY: iniphysiq
   USE mod_const_mpi, ONLY: comm_lmdz
   USE physiq_mod, ONLY: physiq
   USE comvert_mod, ONLY: presnivs, ap, bp, dpres,nivsig, nivsigs, pa, &
                          preff, aps, bps, pseudoalt, scaleheight
   USE temps_mod, ONLY: annee_ref, calend, day_end, day_ini, day_ref, &
                        itau_dyn, itau_phy, start_time, year_len
   USE phys_cal_mod, ONLY : year_len_phys_cal_mod => year_len 

      implicit none
#include "dimensions.h"
#include "YOMCST.h"
!!#include "control.h"
#include "clesphys.h"
#include "dimsoil.h"
!#include "indicesol.h"

#include "compar1d.h"
#include "flux_arp.h"
#include "date_cas.h"
#include "tsoilnudge.h"
#include "fcg_gcssold.h"
!!!#include "fbforcing.h"
#include "compbl.h"

!=====================================================================
! DECLARATIONS
!=====================================================================

!---------------------------------------------------------------------
!  Externals
!---------------------------------------------------------------------
      external fq_sat
      real fq_sat

!---------------------------------------------------------------------
!  Arguments d' initialisations de la physique (USER DEFINE)
!---------------------------------------------------------------------

      integer, parameter :: ngrid=1
      real :: zcufi    = 1.
      real :: zcvfi    = 1.

!-      real :: nat_surf
!-      logical :: ok_flux_surf
!-      real :: fsens
!-      real :: flat
!-      real :: tsurf
!-      real :: rugos
!-      real :: qsol(1:2)
!-      real :: qsurf
!-      real :: psurf
!-      real :: zsurf
!-      real :: albedo
!-
!-      real :: time     = 0.
!-      real :: time_ini
!-      real :: xlat 
!-      real :: xlon 
!-      real :: wtsurf 
!-      real :: wqsurf 
!-      real :: restart_runoff 
!-      real :: xagesno 
!-      real :: qsolinp 
!-      real :: zpicinp 
!-
      real :: fnday 
      real :: day, daytime 
      real :: day1
      real :: heure
      integer :: jour
      integer :: mois
      integer :: an
 
!---------------------------------------------------------------------
!  Declarations related to forcing and initial profiles 
!---------------------------------------------------------------------

        integer :: kmax = llm
        integer llm700,nq1,nq2
        INTEGER, PARAMETER :: nlev_max=1000, nqmx=1000
        real timestep, frac
        real height(nlev_max),tttprof(nlev_max),qtprof(nlev_max)
        real  uprof(nlev_max),vprof(nlev_max),e12prof(nlev_max)
        real  ugprof(nlev_max),vgprof(nlev_max),wfls(nlev_max)
        real  dqtdxls(nlev_max),dqtdyls(nlev_max)
        real  dqtdtls(nlev_max),thlpcar(nlev_max)
        real  qprof(nlev_max,nqmx)

!        integer :: forcing_type
        logical :: forcing_les     = .false.
        logical :: forcing_armcu   = .false.
        logical :: forcing_rico    = .false.
        logical :: forcing_radconv = .false.
        logical :: forcing_toga    = .false.
        logical :: forcing_twpice  = .false.
        logical :: forcing_amma    = .false.
        logical :: forcing_dice    = .false.
        logical :: forcing_gabls4  = .false.

        logical :: forcing_GCM2SCM = .false.
        logical :: forcing_GCSSold = .false.
        logical :: forcing_sandu   = .false.
        logical :: forcing_astex   = .false.
        logical :: forcing_fire    = .false.
        logical :: forcing_case    = .false.
        logical :: forcing_case2   = .false.
        logical :: forcing_SCM   = .false.
        integer :: type_ts_forcing ! 0 = SST constant; 1 = SST read from a file
!                                                            (cf read_tsurf1d.F)

real wwww
!vertical advection computation
!       real d_t_z(llm), d_q_z(llm)
!       real d_t_dyn_z(llm), dq_dyn_z(llm)
!       real zz(llm)
!       real zfact

!flag forcings
        logical :: nudge_wind=.true.
        logical :: nudge_thermo=.false.
        logical :: cptadvw=.true.
!=====================================================================
! DECLARATIONS FOR EACH CASE
!=====================================================================
!
#include "old_1D_decl_cases.h"
!
!---------------------------------------------------------------------
!  Declarations related to nudging
!---------------------------------------------------------------------
     integer :: nudge_max
     parameter (nudge_max=9)
     integer :: inudge_RHT=1
     integer :: inudge_UV=2
     logical :: nudge(nudge_max)
     real :: t_targ(llm)
     real :: rh_targ(llm)
     real :: u_targ(llm)
     real :: v_targ(llm)
!
!---------------------------------------------------------------------
!  Declarations related to vertical discretization:
!---------------------------------------------------------------------
      real :: pzero=1.e5
      real :: play (llm),zlay (llm),sig_s(llm),plev(llm+1)
      real :: playd(llm),zlayd(llm),ap_amma(llm+1),bp_amma(llm+1)

!---------------------------------------------------------------------
!  Declarations related to variables
!---------------------------------------------------------------------

      real :: phi(llm)
      real :: teta(llm),tetal(llm),temp(llm),u(llm),v(llm),w(llm)
      REAL rot(1, llm) ! relative vorticity, in s-1
      real :: rlat_rad(1),rlon_rad(1)
      real :: omega(llm+1),omega2(llm),rho(llm+1)
      real :: ug(llm),vg(llm),fcoriolis
      real :: sfdt, cfdt
      real :: du_phys(llm),dv_phys(llm),dt_phys(llm)
      real :: dt_dyn(llm)
      real :: dt_cooling(llm),d_t_adv(llm),d_t_nudge(llm)
      real :: d_u_nudge(llm),d_v_nudge(llm)
      real :: du_adv(llm),dv_adv(llm)
      real :: du_age(llm),dv_age(llm)
      real :: alpha
      real :: ttt

      REAL, ALLOCATABLE, DIMENSION(:,:):: q
      REAL, ALLOCATABLE, DIMENSION(:,:):: dq
      REAL, ALLOCATABLE, DIMENSION(:,:):: dq_dyn
      REAL, ALLOCATABLE, DIMENSION(:,:):: d_q_adv
      REAL, ALLOCATABLE, DIMENSION(:,:):: d_q_nudge
!      REAL, ALLOCATABLE, DIMENSION(:):: d_th_adv

!---------------------------------------------------------------------
!  Initialization of surface variables
!---------------------------------------------------------------------
      real :: run_off_lic_0(1)
      real :: fder(1),snsrf(1,nbsrf),qsurfsrf(1,nbsrf)
      real :: tsoil(1,nsoilmx,nbsrf)
!     real :: agesno(1,nbsrf)

!---------------------------------------------------------------------
!  Call to phyredem
!---------------------------------------------------------------------
      logical :: ok_writedem =.true.
      real :: sollw_in = 0.
      real :: solsw_in = 0.
      
!---------------------------------------------------------------------
!  Call to physiq
!---------------------------------------------------------------------
      logical :: firstcall=.true.
      logical :: lastcall=.false.
      real :: phis(1)    = 0.0
      real :: dpsrf(1)

!---------------------------------------------------------------------
!  Initializations of boundary conditions 
!---------------------------------------------------------------------
      real, allocatable :: phy_nat (:)  ! 0=ocean libre,1=land,2=glacier,3=banquise
      real, allocatable :: phy_alb (:)  ! Albedo land only (old value condsurf_jyg=0.3)
      real, allocatable :: phy_sst (:)  ! SST (will not be used; cf read_tsurf1d.F)
      real, allocatable :: phy_bil (:)  ! Ne sert que pour les slab_ocean
      real, allocatable :: phy_rug (:) ! Longueur rugosite utilisee sur land only
      real, allocatable :: phy_ice (:) ! Fraction de glace
      real, allocatable :: phy_fter(:) ! Fraction de terre
      real, allocatable :: phy_foce(:) ! Fraction de ocean
      real, allocatable :: phy_fsic(:) ! Fraction de glace
      real, allocatable :: phy_flic(:) ! Fraction de glace

!---------------------------------------------------------------------
!  Fichiers et d'autres variables
!---------------------------------------------------------------------
      integer :: k,l,i,it=1,mxcalc
      integer :: nsrf
      integer jcode
      INTEGER read_climoz 
!
      integer :: it_end ! iteration number of the last call
!Al1
      integer ecrit_slab_oc !1=ecrit,-1=lit,0=no file
      data ecrit_slab_oc/-1/
!
!     if flag_inhib_forcing = 0, tendencies of forcing are added
!                           <> 0, tendencies of forcing are not added
      INTEGER :: flag_inhib_forcing = 0

!=====================================================================
! INITIALIZATIONS 
!=====================================================================
      du_phys(:)=0.
      dv_phys(:)=0.
      dt_phys(:)=0.
      dt_dyn(:)=0.
      dt_cooling(:)=0.
      d_t_adv(:)=0.
      d_t_nudge(:)=0.
      d_u_nudge(:)=0.
      d_v_nudge(:)=0.
      du_adv(:)=0.
      dv_adv(:)=0.
      du_age(:)=0.
      dv_age(:)=0.
      
! Initialization of Common turb_forcing
       dtime_frcg = 0.
       Turb_fcg_gcssold=.false.
       hthturb_gcssold = 0.
       hqturb_gcssold = 0.




!---------------------------------------------------------------------
! OPTIONS OF THE 1D SIMULATION (lmdz1d.def => unicol.def)
!---------------------------------------------------------------------
!Al1
        call conf_unicol
!Al1 moves this gcssold var from common fcg_gcssold to 
        Turb_fcg_gcssold = xTurb_fcg_gcssold
! --------------------------------------------------------------------
        close(1)
!Al1
        write(*,*) 'lmdz1d.def lu => unicol.def'

! forcing_type defines the way the SCM is forced:
!forcing_type = 0 ==> forcing_les = .true.
!             initial profiles from file prof.inp.001
!             no forcing by LS convergence ; 
!             surface temperature imposed ;
!             radiative cooling may be imposed (iflag_radia=0 in physiq.def)
!forcing_type = 1 ==> forcing_radconv = .true.
!             idem forcing_type = 0, but the imposed radiative cooling 
!             is set to 0 (hence, if iflag_radia=0 in physiq.def, 
!             then there is no radiative cooling at all)
!forcing_type = 2 ==> forcing_toga = .true.
!             initial profiles from TOGA-COARE IFA files 
!             LS convergence and SST imposed from TOGA-COARE IFA files 
!forcing_type = 3 ==> forcing_GCM2SCM = .true.
!             initial profiles from the GCM output
!             LS convergence imposed from the GCM output
!forcing_type = 4 ==> forcing_twpice = .true.
!             initial profiles from TWP-ICE cdf file 
!             LS convergence, omega and SST imposed from TWP-ICE files 
!forcing_type = 5 ==> forcing_rico = .true.
!             initial profiles from RICO files 
!             LS convergence imposed from RICO files 
!forcing_type = 6 ==> forcing_amma = .true.
!             initial profiles from AMMA nc file 
!             LS convergence, omega and surface fluxes imposed from AMMA file  
!forcing_type = 7 ==> forcing_dice = .true.
!             initial profiles and large scale forcings in dice_driver.nc
!             Different stages: soil model alone, atm. model alone
!             then both models coupled
!forcing_type = 8 ==> forcing_gabls4 = .true.
!             initial profiles and large scale forcings in gabls4_driver.nc
!forcing_type >= 100 ==> forcing_case = .true.
!             initial profiles and large scale forcings in cas.nc
!             LS convergence, omega and SST imposed from CINDY-DYNAMO files 
!             101=cindynamo
!             102=bomex
!forcing_type >= 100 ==> forcing_case2 = .true.
!             temporary flag while all the 1D cases are not whith the same cas.nc forcing file
!             103=arm_cu2 ie arm_cu with new forcing format
!             104=rico2 ie rico with new forcing format
!forcing_type = 40 ==> forcing_GCSSold = .true.
!             initial profile from GCSS file
!             LS convergence imposed from GCSS file
!forcing_type = 50 ==> forcing_fire = .true.
!             forcing from fire.nc
!forcing_type = 59 ==> forcing_sandu = .true.
!             initial profiles from sanduref file: see prof.inp.001
!             SST varying with time and divergence constante: see ifa_sanduref.txt file
!             Radiation has to be computed interactively
!forcing_type = 60 ==> forcing_astex = .true.
!             initial profiles from file: see prof.inp.001 
!             SST,divergence,ug,vg,ufa,vfa varying with time : see ifa_astex.txt file
!             Radiation has to be computed interactively
!forcing_type = 61 ==> forcing_armcu = .true.
!             initial profiles from file: see prof.inp.001 
!             sensible and latent heat flux imposed: see ifa_arm_cu_1.txt
!             large scale advective forcing & radiative tendencies applied below 1000m: see ifa_arm_cu_2.txt
!             use geostrophic wind ug=10m/s vg=0m/s. Duration of the case 53100s 
!             Radiation to be switched off
!
      if (forcing_type <=0) THEN
       forcing_les = .true.
      elseif (forcing_type .eq.1) THEN
       forcing_radconv = .true.
      elseif (forcing_type .eq.2) THEN
       forcing_toga    = .true.
      elseif (forcing_type .eq.3) THEN
       forcing_GCM2SCM = .true.
      elseif (forcing_type .eq.4) THEN
       forcing_twpice = .true.
      elseif (forcing_type .eq.5) THEN
       forcing_rico = .true.
      elseif (forcing_type .eq.6) THEN
       forcing_amma = .true.
      elseif (forcing_type .eq.7) THEN
       forcing_dice = .true.
      elseif (forcing_type .eq.8) THEN
       forcing_gabls4 = .true.
      elseif (forcing_type .eq.101) THEN ! Cindynamo starts 1-10-2011 0h
       forcing_case = .true.
       year_ini_cas=2011
       mth_ini_cas=10
       day_deb=1
       heure_ini_cas=0.
       pdt_cas=3*3600.         ! forcing frequency
      elseif (forcing_type .eq.102) THEN ! Bomex starts 24-6-1969 0h
       forcing_case = .true.
       year_ini_cas=1969
       mth_ini_cas=6
       day_deb=24
       heure_ini_cas=0.
       pdt_cas=1800.         ! forcing frequency
      elseif (forcing_type .eq.103) THEN ! Arm_cu starts 21-6-1997 11h30
       forcing_case2 = .true.
       year_ini_cas=1997
       mth_ini_cas=6
       day_deb=21
       heure_ini_cas=11.5
       pdt_cas=1800.         ! forcing frequency
      elseif (forcing_type .eq.104) THEN ! rico starts 16-12-2004 0h
       forcing_case2 = .true.
       year_ini_cas=2004
       mth_ini_cas=12
       day_deb=16
       heure_ini_cas=0.
       pdt_cas=1800.         ! forcing frequency
      elseif (forcing_type .eq.105) THEN ! bomex starts 16-12-2004 0h
       forcing_case2 = .true.
       year_ini_cas=1969
       mth_ini_cas=6
       day_deb=24
       heure_ini_cas=0.
       pdt_cas=1800.         ! forcing frequency
      elseif (forcing_type .eq.106) THEN ! ayotte_24SC starts 6-11-1992 0h
       forcing_case2 = .true.
       year_ini_cas=1992
       mth_ini_cas=11
       day_deb=6
       heure_ini_cas=10.
       pdt_cas=86400.        ! forcing frequency
      elseif (forcing_type .eq.113) THEN ! Arm_cu starts 21-6-1997 11h30
       forcing_SCM = .true.
       year_ini_cas=1997
       ! It is possible that those parameters are run twice.
       call getin('anneeref',year_ini_cas)
       call getin('dayref',day_deb)
       mth_ini_cas=1 ! pour le moment on compte depuis le debut de l'annee
       call getin('time_ini',heure_ini_cas)
      elseif (forcing_type .eq.40) THEN
       forcing_GCSSold = .true.
      elseif (forcing_type .eq.50) THEN
       forcing_fire = .true.
      elseif (forcing_type .eq.59) THEN
       forcing_sandu   = .true.
      elseif (forcing_type .eq.60) THEN
       forcing_astex   = .true.
      elseif (forcing_type .eq.61) THEN
       forcing_armcu = .true.
       IF(llm.NE.19.AND.llm.NE.40) stop 'Erreur nombre de niveaux !!'
      else
       write (*,*) 'ERROR : unknown forcing_type ', forcing_type
       stop 'Forcing_type should be 0,1,2,3,4,5,6 or 40,59,60,61'
      ENDIF
      print*,"forcing type=",forcing_type

! if type_ts_forcing=0, the surface temp of 1D simulation is constant in time
! (specified by tsurf in lmdz1d.def); if type_ts_forcing=1, the surface temperature
! varies in time according to a forcing (e.g. forcing_toga) and is passed to read_tsurf1d.F
! through the common sst_forcing.

        type_ts_forcing = 0
        if (forcing_toga.or.forcing_sandu.or.forcing_astex .or. forcing_dice)                 &
     &    type_ts_forcing = 1
!
! Initialization of the logical switch for nudging
     jcode = iflag_nudge
     do i = 1,nudge_max
       nudge(i) = mod(jcode,10) .ge. 1
       jcode = jcode/10
     enddo
!---------------------------------------------------------------------
!  Definition of the run
!---------------------------------------------------------------------

      call conf_gcm( 99, .TRUE. )
      
!-----------------------------------------------------------------------
      allocate( phy_nat (year_len))  ! 0=ocean libre,1=land,2=glacier,3=banquise
      phy_nat(:)=0.0
      allocate( phy_alb (year_len))  ! Albedo land only (old value condsurf_jyg=0.3)
      allocate( phy_sst (year_len))  ! SST (will not be used; cf read_tsurf1d.F)
      allocate( phy_bil (year_len))  ! Ne sert que pour les slab_ocean
      phy_bil(:)=1.0
      allocate( phy_rug (year_len)) ! Longueur rugosite utilisee sur land only
      allocate( phy_ice (year_len)) ! Fraction de glace
      phy_ice(:)=0.0
      allocate( phy_fter(year_len)) ! Fraction de terre
      phy_fter(:)=0.0
      allocate( phy_foce(year_len)) ! Fraction de ocean
      phy_foce(:)=0.0
      allocate( phy_fsic(year_len)) ! Fraction de glace
      phy_fsic(:)=0.0
      allocate( phy_flic(year_len)) ! Fraction de glace
      phy_flic(:)=0.0
!-----------------------------------------------------------------------
!   Choix du calendrier
!   -------------------

!      calend = 'earth_365d'
      if (calend == 'earth_360d') then
        call ioconf_calendar('360d')
        write(*,*)'CALENDRIER CHOISI: Terrestre a 360 jours/an'
      else if (calend == 'earth_365d') then
        call ioconf_calendar('noleap')
        write(*,*)'CALENDRIER CHOISI: Terrestre a 365 jours/an'
      else if (calend == 'earth_366d') then
        call ioconf_calendar('all_leap')
        write(*,*)'CALENDRIER CHOISI: Terrestre bissextile'
      else if (calend == 'gregorian') then
        stop 'gregorian calend should not be used by normal user'
        call ioconf_calendar('gregorian') ! not to be used by normal users
        write(*,*)'CALENDRIER CHOISI: Gregorien'
      else
        write (*,*) 'ERROR : unknown calendar ', calend
        stop 'calend should be 360d,earth_365d,earth_366d,gregorian'
      endif
!-----------------------------------------------------------------------
!
!c Date :
!      La date est supposee donnee sous la forme [annee, numero du jour dans 
!      l annee] ; l heure est donnee dans time_ini, lu dans lmdz1d.def.
!      On appelle ymds2ju pour convertir [annee, jour] en [jour Julien].
!      Le numero du jour est dans "day". L heure est traitee separement.
!      La date complete est dans "daytime" (l'unite est le jour).
      if (nday>0) then
         fnday=nday
      else
         fnday=-nday/float(day_step)
      endif
      print *,'fnday=',fnday
!     start_time doit etre en FRACTION DE JOUR
      start_time=time_ini/24.

! Special case for arm_cu which lasts less than one day : 53100s !! (MPL 20111026)
      IF(forcing_type .EQ. 61) fnday=53100./86400.
      IF(forcing_type .EQ. 103) fnday=53100./86400.
! Special case for amma which lasts less than one day : 64800s !! (MPL 20120216)
      IF(forcing_type .EQ. 6) fnday=64800./86400.
!     IF(forcing_type .EQ. 6) fnday=50400./86400.
 IF(forcing_type .EQ. 8 ) fnday=129600./86400. 
      annee_ref = anneeref
      mois = 1
      day_ref = dayref
      heure = 0.
      itau_dyn = 0
      itau_phy = 0
      call ymds2ju(annee_ref,mois,day_ref,heure,day)
      day_ini = int(day)
      day_end = day_ini + int(fnday)

      IF (forcing_type .eq.2) THEN
! Convert the initial date of Toga-Coare to Julian day
      call ymds2ju                                                          &
     & (year_ini_toga,mth_ini_toga,day_ini_toga,heure,day_ju_ini_toga)

      ELSEIF (forcing_type .eq.4) THEN
! Convert the initial date of TWPICE to Julian day
      call ymds2ju                                                          &
     & (year_ini_twpi,mth_ini_twpi,day_ini_twpi,heure_ini_twpi              &
     & ,day_ju_ini_twpi)
      ELSEIF (forcing_type .eq.6) THEN
! Convert the initial date of AMMA to Julian day
      call ymds2ju                                                          &
     & (year_ini_amma,mth_ini_amma,day_ini_amma,heure_ini_amma              &
     & ,day_ju_ini_amma)
      ELSEIF (forcing_type .eq.7) THEN
! Convert the initial date of DICE to Julian day
      call ymds2ju                                                         &
     & (year_ini_dice,mth_ini_dice,day_ini_dice,heure_ini_dice             & 
     & ,day_ju_ini_dice)
 ELSEIF (forcing_type .eq.8 ) THEN
! Convert the initial date of GABLS4 to Julian day
      call ymds2ju                                                         &
     & (year_ini_gabls4,mth_ini_gabls4,day_ini_gabls4,heure_ini_gabls4     & 
     & ,day_ju_ini_gabls4)
      ELSEIF (forcing_type .gt.100) THEN
! Convert the initial date to Julian day
      day_ini_cas=day_deb
      print*,'time case',year_ini_cas,mth_ini_cas,day_ini_cas
      call ymds2ju                                                         &
     & (year_ini_cas,mth_ini_cas,day_ini_cas,heure_ini_cas*3600            &
     & ,day_ju_ini_cas)
      print*,'time case 2',day_ini_cas,day_ju_ini_cas
      ELSEIF (forcing_type .eq.59) THEN
! Convert the initial date of Sandu case to Julian day
      call ymds2ju                                                          &
     &   (year_ini_sandu,mth_ini_sandu,day_ini_sandu,                       &
     &    time_ini*3600.,day_ju_ini_sandu)

      ELSEIF (forcing_type .eq.60) THEN
! Convert the initial date of Astex case to Julian day
      call ymds2ju                                                          &
     &   (year_ini_astex,mth_ini_astex,day_ini_astex,                        &
     &    time_ini*3600.,day_ju_ini_astex)

      ELSEIF (forcing_type .eq.61) THEN
! Convert the initial date of Arm_cu case to Julian day
      call ymds2ju                                                          &
     & (year_ini_armcu,mth_ini_armcu,day_ini_armcu,heure_ini_armcu          &
     & ,day_ju_ini_armcu)
      ENDIF

      IF (forcing_type .gt.100) THEN
      daytime = day + heure_ini_cas/24. ! 1st day and initial time of the simulation 
      ELSE
      daytime = day + time_ini/24. ! 1st day and initial time of the simulation 
      ENDIF
! Print out the actual date of the beginning of the simulation :
      call ju2ymds(daytime,year_print, month_print,day_print,sec_print)
      print *,' Time of beginning : ',                                      &
     &        year_print, month_print, day_print, sec_print

!---------------------------------------------------------------------
! Initialization of dimensions, geometry and initial state
!---------------------------------------------------------------------
!      call init_phys_lmdz(1,1,llm,1,(/1/)) ! job now done via iniphysiq
!     but we still need to initialize dimphy module (klon,klev,etc.)  here.
      call init_dimphy1D(1,llm)
      call suphel
      call infotrac_init 

      if (nqtot>nqmx) STOP 'Augmenter nqmx dans lmdz1d.F'
      allocate(q(llm,nqtot)) ; q(:,:)=0.
      allocate(dq(llm,nqtot)) 
      allocate(dq_dyn(llm,nqtot)) 
      allocate(d_q_adv(llm,nqtot)) 
      allocate(d_q_nudge(llm,nqtot)) 
!      allocate(d_th_adv(llm)) 

      q(:,:) = 0.
      dq(:,:) = 0.
      dq_dyn(:,:) = 0.
      d_q_adv(:,:) = 0.
      d_q_nudge(:,:) = 0.

!
!   No ozone climatology need be read in this pre-initialization
!          (phys_state_var_init is called again in physiq)
      read_climoz = 0
!
      call phys_state_var_init(read_climoz)

      if (ngrid.ne.klon) then
         print*,'stop in inifis'
         print*,'Probleme de dimensions :'
         print*,'ngrid = ',ngrid
         print*,'klon  = ',klon
         stop
      endif
!!!=====================================================================
!!! Feedback forcing values for Gateaux differentiation (al1)
!!!=====================================================================
!!! Surface Planck forcing bracketing call radiation
!!      surf_Planck = 0.
!!      surf_Conv   = 0.
!!      write(*,*) 'Gateaux-dif Planck,Conv:',surf_Planck,surf_Conv
!!! a mettre dans le lmdz1d.def ou autre
!!
!!
      qsol = qsolinp
      qsurf = fq_sat(tsurf,psurf/100.)
      day1= day_ini
      time=daytime-day
      ts_toga(1)=tsurf ! needed by read_tsurf1d.F
      rho(1)=psurf/(rd*tsurf*(1.+(rv/rd-1.)*qsurf)) 

!
!! mpl et jyg le 22/08/2012 : 
!!  pour que les cas a flux de surface imposes marchent
      IF(.NOT.ok_flux_surf.or.max(abs(wtsurf),abs(wqsurf))>0.) THEN
       fsens=-wtsurf*rcpd*rho(1)
       flat=-wqsurf*rlvtt*rho(1)
       print *,'Flux: ok_flux wtsurf wqsurf',ok_flux_surf,wtsurf,wqsurf
      ENDIF
      print*,'Flux sol ',fsens,flat
!!      ok_flux_surf=.false.
!!      fsens=-wtsurf*rcpd*rho(1)
!!      flat=-wqsurf*rlvtt*rho(1)
!!!!

! Vertical discretization and pressure levels at half and mid levels:

      pa   = 5e4
!!      preff= 1.01325e5
      preff = psurf
      IF (ok_old_disvert) THEN
        call disvert0(pa,preff,ap,bp,dpres,presnivs,nivsigs,nivsig)
        print *,'On utilise disvert0'
        aps(1:llm)=0.5*(ap(1:llm)+ap(2:llm+1))
        bps(1:llm)=0.5*(bp(1:llm)+bp(2:llm+1))
        scaleheight=8.
        pseudoalt(1:llm)=-scaleheight*log(presnivs(1:llm)/preff)
      ELSE
        call disvert()
        print *,'On utilise disvert'
!       Nouvelle version disvert permettant d imposer ap,bp (modif L.Guez) MPL 18092012
!       Dans ce cas, on lit ap,bp dans le fichier hybrid.txt
      ENDIF

      sig_s=presnivs/preff
      plev =ap+bp*psurf
      play = 0.5*(plev(1:llm)+plev(2:llm+1))
      zlay=-rd*300.*log(play/psurf)/rg ! moved after reading profiles

      IF (forcing_type .eq. 59) THEN
! pour forcing_sandu, on cherche l'indice le plus proche de 700hpa#3000m
      write(*,*) '***********************'
      do l = 1, llm
       write(*,*) 'l,play(l),presnivs(l): ',l,play(l),presnivs(l)
       if (trouve_700 .and. play(l).le.70000) then
         llm700=l
         print *,'llm700,play=',llm700,play(l)/100.
         trouve_700= .false.
       endif
      enddo
      write(*,*) '***********************'
      ENDIF

!
!=====================================================================
! EVENTUALLY, READ FORCING DATA :
!=====================================================================

#include "old_1D_read_forc_cases.h"

      if (forcing_GCM2SCM) then
        write (*,*) 'forcing_GCM2SCM not yet implemented'
        stop 'in initialization'
      endif ! forcing_GCM2SCM

      print*,'mxcalc=',mxcalc
!     print*,'zlay=',zlay(mxcalc)
      print*,'play=',play(mxcalc)

!Al1 pour SST forced, appell?? depuis ocean_forced_noice
      ts_cur = tsurf ! SST used in read_tsurf1d
!=====================================================================
! Initialisation de la physique : 
!=====================================================================

!  Rq: conf_phys.F90 lit tous les flags de physiq.def; conf_phys appele depuis physiq.F
!
! day_step, iphysiq lus dans gcm.def ci-dessus
! timestep: calcule ci-dessous from rday et day_step
! ngrid=1
! llm: defini dans .../modipsl/modeles/LMDZ4/libf/grid/dimension
! rday: defini dans suphel.F (86400.)
! day_ini: lu dans run.def (dayref)
! rlat_rad,rlon-rad: transformes en radian de rlat,rlon lus dans lmdz1d.def (en degres)
! airefi,zcufi,zcvfi initialises au debut de ce programme
! rday,ra,rg,rd,rcpd declares dans YOMCST.h et calcules dans suphel.F
      day_step = float(nsplit_phys)*day_step/float(iphysiq)
      write (*,*) 'Time step divided by nsplit_phys (=',nsplit_phys,')'
      timestep =rday/day_step
      dtime_frcg = timestep
!
      zcufi=airefi
      zcvfi=airefi
!
      rlat_rad(1)=xlat*rpi/180.
      rlon_rad(1)=xlon*rpi/180.

     ! iniphysiq will call iniaqua who needs year_len from phys_cal_mod
     year_len_phys_cal_mod=year_len
           
     ! Ehouarn: iniphysiq requires arrays related to (3D) dynamics grid,
     ! e.g. for cell boundaries, which are meaningless in 1D; so pad these 
     ! with '0.' when necessary
      call iniphysiq(iim,jjm,llm, &
           1,comm_lmdz, &
           rday,day_ini,timestep,  &
           (/rlat_rad(1),0./),(/0./), &
           (/0.,0./),(/rlon_rad(1),0./),  &
           (/ (/airefi,0./),(/0.,0./) /), &
           (/zcufi,0.,0.,0./), &
           (/zcvfi,0./), &
           ra,rg,rd,rcpd,1)
      print*,'apres iniphysiq'

! 2 PARAMETRES QUI DEVRAIENT ETRE LUS DANS run.def MAIS NE LE SONT PAS ICI:
      co2_ppm= 330.0
      solaire=1370.0

! Ecriture du startphy avant le premier appel a la physique.
! On le met juste avant pour avoir acces a tous les champs

      if (ok_writedem) then

!--------------------------------------------------------------------------
! pbl_surface_init (called here) and pbl_surface_final (called by phyredem)
! need : qsol fder snow qsurf evap rugos agesno ftsoil
!--------------------------------------------------------------------------

        type_ocean = "force"
        run_off_lic_0(1) = restart_runoff 
        call fonte_neige_init(run_off_lic_0)

        fder=0.
        snsrf(1,:)=snowmass ! masse de neige des sous surface
        qsurfsrf(1,:)=qsurf ! humidite de l'air des sous surface
        fevap=0.
        z0m(1,:)=rugos     ! couverture de neige des sous surface
        z0h(1,:)=rugosh    ! couverture de neige des sous surface
        agesno  = xagesno
        tsoil(:,:,:)=tsurf
!------ AMMA 2e run avec modele sol et rayonnement actif (MPL 23052012)
!       tsoil(1,1,1)=299.18
!       tsoil(1,2,1)=300.08
!       tsoil(1,3,1)=301.88
!       tsoil(1,4,1)=305.48
!       tsoil(1,5,1)=308.00
!       tsoil(1,6,1)=308.00
!       tsoil(1,7,1)=308.00
!       tsoil(1,8,1)=308.00
!       tsoil(1,9,1)=308.00
!       tsoil(1,10,1)=308.00
!       tsoil(1,11,1)=308.00
!-----------------------------------------------------------------------
        call pbl_surface_init(fder, snsrf, qsurfsrf, tsoil)

!------------------ prepare limit conditions for limit.nc -----------------
!--   Ocean force

        print*,'avant phyredem'
        pctsrf(1,:)=0.
          if (nat_surf.eq.0.) then
          pctsrf(1,is_oce)=1.
          pctsrf(1,is_ter)=0.
          pctsrf(1,is_lic)=0.
          pctsrf(1,is_sic)=0.
        else if (nat_surf .eq. 1) then 
          pctsrf(1,is_oce)=0.
          pctsrf(1,is_ter)=1.
          pctsrf(1,is_lic)=0.
          pctsrf(1,is_sic)=0.
        else if (nat_surf .eq. 2) then 
          pctsrf(1,is_oce)=0.
          pctsrf(1,is_ter)=0.
          pctsrf(1,is_lic)=1.
          pctsrf(1,is_sic)=0.
        else if (nat_surf .eq. 3) then 
          pctsrf(1,is_oce)=0.
          pctsrf(1,is_ter)=0.
          pctsrf(1,is_lic)=0.
          pctsrf(1,is_sic)=1.

     end if


        print*,'nat_surf,pctsrf(1,is_oce),pctsrf(1,is_ter)',nat_surf         &
     &        ,pctsrf(1,is_oce),pctsrf(1,is_ter)

        zmasq=pctsrf(1,is_ter)+pctsrf(1,is_lic)
        zpic = zpicinp
        ftsol=tsurf
        nsw=6 ! on met le nb de bandes SW=6, pour initialiser
              ! 6 albedo, mais on peut quand meme tourner avec
              ! moins. Seules les 2 ou 4 premiers seront lus
        falb_dir=albedo
        falb_dif=albedo
        rugoro=rugos 
        t_ancien(1,:)=temp(:)
        q_ancien(1,:)=q(:,1)
        ql_ancien = 0.
        qs_ancien = 0.
        prlw_ancien = 0.
        prsw_ancien = 0.
        prw_ancien = 0.
!jyg<
!!        pbl_tke(:,:,:)=1.e-8
        pbl_tke(:,:,:)=0.
        pbl_tke(:,2,:)=1.e-2
        PRINT *, ' pbl_tke dans lmdz1d '
        if (prt_level .ge. 5) then 
         DO nsrf = 1,4
           PRINT *,'pbl_tke(1,:,',nsrf,') ',pbl_tke(1,:,nsrf)
         ENDDO
        end if

!>jyg

        rain_fall=0.
        snow_fall=0.
        solsw=0.
        sollw=0.
        sollwdown=rsigma*tsurf**4
        radsol=0.
        rnebcon=0.
        ratqs=0.
        clwcon=0.
        zmax0 = 0.
        zmea=0.
        zstd=0.
        zsig=0.
        zgam=0.
        zval=0. 
        zthe=0. 
        sig1=0.
        w01=0.
        wake_cstar = 0.
        wake_deltaq = 0.
        wake_deltat = 0.
        wake_delta_pbl_TKE(:,:,:) = 0.
        delta_tsurf = 0.
        wake_fip = 0.
        wake_pe = 0.
        wake_s = 0.
        wake_dens = 0.
        ale_bl = 0.
        ale_bl_trig = 0.
        alp_bl = 0.
        IF (ALLOCATED(du_gwd_rando)) du_gwd_rando = 0.
        IF (ALLOCATED(du_gwd_front)) du_gwd_front = 0.
        entr_therm = 0.
        detr_therm = 0.
        f0 = 0.
        fm_therm = 0.
        u_ancien(1,:)=u(:)
        v_ancien(1,:)=v(:)

u10m=0.
v10m=0.
ale_wake=0.
ale_bl_stat=0.
 
!------------------------------------------------------------------------
! Make file containing restart for the physics (startphy.nc)
!
! NB: List of the variables to be written by phyredem (via put_field):
! rlon,rlat,zmasq,pctsrf(:,is_ter),pctsrf(:,is_lic),pctsrf(:,is_oce)
! pctsrf(:,is_sic),ftsol(:,nsrf),tsoil(:,isoil,nsrf),qsurf(:,nsrf)
! qsol,falb_dir(:,nsrf),falb_dif(:,nsrf),evap(:,nsrf),snow(:,nsrf)
! radsol,solsw,sollw, sollwdown,fder,rain_fall,snow_fall,frugs(:,nsrf)
! agesno(:,nsrf),zmea,zstd,zsig,zgam,zthe,zpic,zval,rugoro
! t_ancien,q_ancien,,frugs(:,is_oce),clwcon(:,1),rnebcon(:,1),ratqs(:,1)
! run_off_lic_0,pbl_tke(:,1:klev,nsrf), zmax0,f0,sig1,w01
! wake_deltat,wake_deltaq,wake_s,wake_dens,wake_cstar,
! wake_fip,wake_delta_pbl_tke(:,1:klev,nsrf)
!
! NB2: The content of the startphy.nc file depends on some flags defined in
! the ".def" files. However, since conf_phys is not called in lmdz1d.F90, these flags have 
! to be set at some arbitratry convenient values.
!------------------------------------------------------------------------
!Al1 =============== restart option ==========================
        if (.not.restart) then
          iflag_pbl = 5
          call phyredem ("startphy.nc")
        else
! (desallocations)
        print*,'callin surf final'
          call pbl_surface_final( fder, snsrf, qsurfsrf, tsoil)
        print*,'after surf final'
          CALL fonte_neige_final(run_off_lic_0)
        endif

        ok_writedem=.false.
        print*,'apres phyredem'

      endif ! ok_writedem
      
!------------------------------------------------------------------------
! Make file containing boundary conditions (limit.nc) **Al1->restartdyn***
! --------------------------------------------------
! NB: List of the variables to be written in limit.nc 
!     (by writelim.F, subroutine of 1DUTILS.h):
!        phy_nat,phy_alb,phy_sst,phy_bil,phy_rug,phy_ice,
!        phy_fter,phy_foce,phy_flic,phy_fsic)
!------------------------------------------------------------------------
      do i=1,year_len
        phy_nat(i)  = nat_surf
        phy_alb(i)  = albedo
        phy_sst(i)  = tsurf ! read_tsurf1d will be used instead
        phy_rug(i)  = rugos
        phy_fter(i) = pctsrf(1,is_ter)
        phy_foce(i) = pctsrf(1,is_oce)
        phy_fsic(i) = pctsrf(1,is_sic)
        phy_flic(i) = pctsrf(1,is_lic)
      enddo

! fabrication de limit.nc
      call writelim (1,phy_nat,phy_alb,phy_sst,phy_bil,phy_rug,             &
     &               phy_ice,phy_fter,phy_foce,phy_flic,phy_fsic)


      call phys_state_var_end
!Al1
      if (restart) then
        print*,'call to restart dyn 1d'
        Call dyn1deta0("start1dyn.nc",plev,play,phi,phis,presnivs,          &
     &              u,v,temp,q,omega2)

       print*,'fnday,annee_ref,day_ref,day_ini',                            &
     &     fnday,annee_ref,day_ref,day_ini
!**      call ymds2ju(annee_ref,mois,day_ini,heure,day)
       day = day_ini
       day_end = day_ini + nday
       daytime = day + time_ini/24. ! 1st day and initial time of the simulation 

! Print out the actual date of the beginning of the simulation :
       call ju2ymds(daytime, an, mois, jour, heure)
       print *,' Time of beginning : y m d h',an, mois,jour,heure/3600.

       day = int(daytime)
       time=daytime-day
  
       print*,'****** intialised fields from restart1dyn *******'
       print*,'plev,play,phi,phis,presnivs,u,v,temp,q,omega2'
       print*,'temp(1),q(1,1),u(1),v(1),plev(1),phis :'
       print*,temp(1),q(1,1),u(1),v(1),plev(1),phis
! raz for safety
       do l=1,llm
         dq_dyn(l,1) = 0.
       enddo
      endif
!Al1 ================  end restart =================================
      IF (ecrit_slab_oc.eq.1) then
         open(97,file='div_slab.dat',STATUS='UNKNOWN')
       elseif (ecrit_slab_oc.eq.0) then
         open(97,file='div_slab.dat',STATUS='OLD')
       endif
!
!---------------------------------------------------------------------
!    Initialize target profile for RHT nudging if needed
!---------------------------------------------------------------------
      if (nudge(inudge_RHT)) then
        call nudge_RHT_init(plev,play,temp,q(:,1),t_targ,rh_targ)
      endif
      if (nudge(inudge_UV)) then
        call nudge_UV_init(plev,play,u,v,u_targ,v_targ)
      endif
!
!=====================================================================
       CALL iophys_ini
! START OF THE TEMPORAL LOOP :
!=====================================================================
           
      it_end = nint(fnday*day_step)
!test JLD     it_end = 10
      do while(it.le.it_end)

       if (prt_level.ge.1) then
         print*,'XXXXXXXXXXXXXXXXXXX ITAP,day,time=',                       &
     &             it,day,time,it_end,day_step
         print*,'PAS DE TEMPS ',timestep
       endif
!Al1 demande de restartphy.nc
       if (it.eq.it_end) lastcall=.True.

!---------------------------------------------------------------------
! Interpolation of forcings in time and onto model levels
!---------------------------------------------------------------------

#include "old_1D_interp_cases.h"

      if (forcing_GCM2SCM) then
        write (*,*) 'forcing_GCM2SCM not yet implemented'
        stop 'in time loop'
      endif ! forcing_GCM2SCM

!---------------------------------------------------------------------
!  Geopotential :
!---------------------------------------------------------------------

        phi(1)=RD*temp(1)*(plev(1)-play(1))/(.5*(plev(1)+play(1)))
        do l = 1, llm-1
          phi(l+1)=phi(l)+RD*(temp(l)+temp(l+1))*                           &
     &    (play(l)-play(l+1))/(play(l)+play(l+1))
        enddo

!---------------------------------------------------------------------
! Listing output for debug prt_level>=1
!---------------------------------------------------------------------
       if (prt_level>=1) then
         print *,' avant physiq : -------- day time ',day,time
         write(*,*) 'firstcall,lastcall,phis',                               &
     &               firstcall,lastcall,phis
       end if
       if (prt_level>=5) then
         write(*,'(a10,2a4,4a13)') 'BEFOR1 IT=','it','l',                   &
     &        'presniv','plev','play','phi'
         write(*,'(a10,2i4,4f13.2)') ('BEFOR1 IT= ',it,l,                   &
     &         presnivs(l),plev(l),play(l),phi(l),l=1,llm)
         write(*,'(a11,2a4,a11,6a8)') 'BEFOR2','it','l',                    &
     &         'presniv','u','v','temp','q1','q2','omega2'
         write(*,'(a11,2i4,f11.2,5f8.2,e10.2)') ('BEFOR2 IT= ',it,l,         &
     &   presnivs(l),u(l),v(l),temp(l),q(l,1),q(l,2),omega2(l),l=1,llm)
       endif

!---------------------------------------------------------------------
!   Call physiq :
!---------------------------------------------------------------------
       call physiq(ngrid,llm, &
                    firstcall,lastcall,timestep, &
                    plev,play,phi,phis,presnivs, &
                    u,v, rot, temp,q,omega2, &
                    du_phys,dv_phys,dt_phys,dq,dpsrf)
                firstcall=.false.

!---------------------------------------------------------------------
! Listing output for debug
!---------------------------------------------------------------------
        if (prt_level>=5) then
          write(*,'(a11,2a4,4a13)') 'AFTER1 IT=','it','l',                  &
     &        'presniv','plev','play','phi'
          write(*,'(a11,2i4,4f13.2)') ('AFTER1 it= ',it,l,                  &
     &    presnivs(l),plev(l),play(l),phi(l),l=1,llm)
          write(*,'(a11,2a4,a11,6a8)') 'AFTER2','it','l',                   &
     &         'presniv','u','v','temp','q1','q2','omega2'
          write(*,'(a11,2i4,f11.2,5f8.2,e10.2)') ('AFTER2 it= ',it,l,       &
     &    presnivs(l),u(l),v(l),temp(l),q(l,1),q(l,2),omega2(l),l=1,llm)
          write(*,'(a11,2a4,a11,5a8)') 'AFTER3','it','l',                   &
     &         'presniv','du_phys','dv_phys','dt_phys','dq1','dq2'    
           write(*,'(a11,2i4,f11.2,5f8.2)') ('AFTER3 it= ',it,l,            &
     &      presnivs(l),86400*du_phys(l),86400*dv_phys(l),                   &
     &       86400*dt_phys(l),86400*dq(l,1),dq(l,2),l=1,llm)
          write(*,*) 'dpsrf',dpsrf
        endif
!---------------------------------------------------------------------
!   Add physical tendencies :
!---------------------------------------------------------------------

       fcoriolis=2.*sin(rpi*xlat/180.)*romega
       if (forcing_radconv .or. forcing_fire) then
         fcoriolis=0.0
         dt_cooling=0.0
         d_t_adv=0.0
         d_q_adv=0.0
       endif
!      print*, 'calcul de fcoriolis ', fcoriolis

       if (forcing_toga .or. forcing_GCSSold .or. forcing_twpice            &
     &    .or.forcing_amma .or. forcing_type.eq.101) then
         fcoriolis=0.0 ; ug=0. ; vg=0.
       endif

       if(forcing_rico) then
          dt_cooling=0.
       endif

!CRio:Attention modif sp??cifique cas de Caroline
      if (forcing_type==-1) then
         fcoriolis=0.
!Nudging
        
!on calcule dt_cooling
        do l=1,llm
        if (play(l).ge.20000.) then
            dt_cooling(l)=-1.5/86400.
        elseif ((play(l).ge.10000.).and.((play(l).lt.20000.))) then
            dt_cooling(l)=-1.5/86400.*(play(l)-10000.)/(10000.)-1./86400.*(20000.-play(l))/10000.*(temp(l)-200.)
        else
            dt_cooling(l)=-1.*(temp(l)-200.)/86400.
        endif
        enddo

      endif     
!RC
      if (forcing_sandu) then
         ug(1:llm)=u_mod(1:llm)
         vg(1:llm)=v_mod(1:llm)
      endif

      IF (prt_level >= 5) print*, 'fcoriolis, xlat,mxcalc ', &
                                   fcoriolis, xlat,mxcalc

!       print *,'u-ug=',u-ug

!!!!!!!!!!!!!!!!!!!!!!!!
! Geostrophic wind
! Le calcul ci dessous est insuffisamment precis
!      du_age(1:mxcalc)=fcoriolis*(v(1:mxcalc)-vg(1:mxcalc))
!      dv_age(1:mxcalc)=-fcoriolis*(u(1:mxcalc)-ug(1:mxcalc))
!!!!!!!!!!!!!!!!!!!!!!!!
       sfdt = sin(0.5*fcoriolis*timestep)
       cfdt = cos(0.5*fcoriolis*timestep)
!       print *,'fcoriolis,sfdt,cfdt,timestep',fcoriolis,sfdt,cfdt,timestep
!
        du_age(1:mxcalc)= -2.*sfdt/timestep*                                &
     &          (sfdt*(u(1:mxcalc)-ug(1:mxcalc)) -                          &
     &           cfdt*(v(1:mxcalc)-vg(1:mxcalc))  )
!!     : fcoriolis*(v(1:mxcalc)-vg(1:mxcalc))
!
       dv_age(1:mxcalc)= -2.*sfdt/timestep*                                 &
     &          (cfdt*(u(1:mxcalc)-ug(1:mxcalc)) +                           &
     &           sfdt*(v(1:mxcalc)-vg(1:mxcalc))  )
!!     : -fcoriolis*(u(1:mxcalc)-ug(1:mxcalc))
!
!!!!!!!!!!!!!!!!!!!!!!!!
!  Nudging
!!!!!!!!!!!!!!!!!!!!!!!!
      d_t_nudge(:) = 0.
      d_q_nudge(:,:) = 0.
      d_u_nudge(:) = 0.
      d_v_nudge(:) = 0.
      if (nudge(inudge_RHT)) then
        call nudge_RHT(timestep,plev,play,t_targ,rh_targ,temp,q(:,1),     &
    &                  d_t_nudge,d_q_nudge(:,1))
      endif
      if (nudge(inudge_UV)) then
        call nudge_UV(timestep,plev,play,u_targ,v_targ,u,v,     &
    &                  d_u_nudge,d_v_nudge)
      endif
!
       if (forcing_fire) THEN

!let ww=if ( alt le 1100 ) then alt*-0.00001 else 0
!let wt=if ( alt le 1100 ) then min( -3.75e-5 , -7.5e-8*alt)  else 0
!let wq=if ( alt le 1100 ) then max( 1.5e-8 , 3e-11*alt)  else 0
           d_t_adv=0.
           d_q_adv=0.
           teta=temp*(pzero/play)**rkappa
           d_t_adv=0.
           d_q_adv=0.
           do l=2,llm-1
              if (zlay(l)<=1100) then
                  wwww=-0.00001*zlay(l)
                  d_t_adv(l)=-wwww*(teta(l)-teta(l+1))/(zlay(l)-zlay(l+1)) /(pzero/play(l))**rkappa
                  d_q_adv(l,1:2)=-wwww*(q(l,1:2)-q(l+1,1:2))/(zlay(l)-zlay(l+1))
                  d_t_adv(l)=d_t_adv(l)+min(-3.75e-5 , -7.5e-8*zlay(l))
                  d_q_adv(l,1)=d_q_adv(l,1)+max( 1.5e-8 , 3e-11*zlay(l))
              endif
           enddo

        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         call  writefield_phy('dv_age' ,dv_age,llm)
!         call  writefield_phy('du_age' ,du_age,llm)
!         call  writefield_phy('du_phys' ,du_phys,llm)
!         call  writefield_phy('u_tend' ,u,llm)
!         call  writefield_phy('u_g' ,ug,llm)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Increment state variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*,'OLDLMDZ1D AANT'
    IF (flag_inhib_forcing == 0) then ! if tendency of forcings should be added
print*,'OLDLMDZ1D ARES' , forcing_sandu

! pour les cas sandu et astex, on reclacule u,v,q,temp et teta dans 1D_nudge_sandu_astex.h
! au dessus de 700hpa, on relaxe vers les profils initiaux
      if (forcing_sandu .OR. forcing_astex) then
#include "1D_nudge_sandu_astex.h"
      else
        u(1:mxcalc)=u(1:mxcalc) + timestep*(                                &
     &              du_phys(1:mxcalc)                                       &
     &             +du_age(1:mxcalc)+du_adv(1:mxcalc)                       &
     &             +d_u_nudge(1:mxcalc) )            
        v(1:mxcalc)=v(1:mxcalc) + timestep*(                                 &
     &              dv_phys(1:mxcalc)                                       &
     &             +dv_age(1:mxcalc)+dv_adv(1:mxcalc)                       &
     &             +d_v_nudge(1:mxcalc) )
        q(1:mxcalc,:)=q(1:mxcalc,:)+timestep*(                              &
     &                dq(1:mxcalc,:)                                        &
     &               +d_q_adv(1:mxcalc,:)                                   &
     &               +d_q_nudge(1:mxcalc,:) )

        if (prt_level.ge.3) then
          print *,                                                          &
     &    'physiq-> temp(1),dt_phys(1),d_t_adv(1),dt_cooling(1) ',         &
     &              temp(1),dt_phys(1),d_t_adv(1),dt_cooling(1)
           print* ,'dv_phys=',dv_phys
           print* ,'dv_age=',dv_age
           print* ,'dv_adv=',dv_adv
           print* ,'d_v_nudge=',d_v_nudge
           print*, v
           print*, vg
        endif

        temp(1:mxcalc)=temp(1:mxcalc)+timestep*(                            &
     &              dt_phys(1:mxcalc)                                       &
     &             +d_t_adv(1:mxcalc)                                      &
     &             +d_t_nudge(1:mxcalc)                                      &
     &             +dt_cooling(1:mxcalc))  ! Taux de chauffage ou refroid.

print*,'OLDLMDZ1D IOPH'
      CALL iophys_ecrit('d_t_adv',klev,'d_t_adv','m/s',d_t_adv)
      CALL iophys_ecrit('d_t_nudge',klev,'d_t_nudge','m/s',d_t_nudge)

      endif  ! forcing_sandu or forcing_astex

        teta=temp*(pzero/play)**rkappa
!
!---------------------------------------------------------------------
!   Nudge soil temperature if requested
!---------------------------------------------------------------------

      IF (nudge_tsoil .AND. .NOT. lastcall) THEN
       ftsoil(1,isoil_nudge,:) = ftsoil(1,isoil_nudge,:)                     &
     &  -timestep/tau_soil_nudge*(ftsoil(1,isoil_nudge,:)-Tsoil_nudge)
      ENDIF

!---------------------------------------------------------------------
!   Add large-scale tendencies (advection, etc) :
!---------------------------------------------------------------------

!cc nrlmd
!cc        tmpvar=teta
!cc        call advect_vert(llm,omega,timestep,tmpvar,plev)
!cc
!cc        teta(1:mxcalc)=tmpvar(1:mxcalc)
!cc        tmpvar(:)=q(:,1)
!cc        call advect_vert(llm,omega,timestep,tmpvar,plev)
!cc        q(1:mxcalc,1)=tmpvar(1:mxcalc)
!cc        tmpvar(:)=q(:,2)
!cc        call advect_vert(llm,omega,timestep,tmpvar,plev)
!cc        q(1:mxcalc,2)=tmpvar(1:mxcalc)

   END IF ! end if tendency of tendency should be added

!---------------------------------------------------------------------
!   Air temperature :
!---------------------------------------------------------------------        
        if (lastcall) then
          print*,'Pas de temps final ',it
          call ju2ymds(daytime, an, mois, jour, heure)
          print*,'a la date : a m j h',an, mois, jour ,heure/3600.
        endif

!  incremente day time
!        print*,'daytime bef',daytime,1./day_step
        daytime = daytime+1./day_step
!Al1dbg
        day = int(daytime+0.1/day_step)
!        time = max(daytime-day,0.0)
!Al1&jyg: correction de bug
!cc        time = real(mod(it,day_step))/day_step
        time = time_ini/24.+real(mod(it,day_step))/day_step
!        print*,'daytime nxt time',daytime,time
        it=it+1

      enddo

!Al1
      if (ecrit_slab_oc.ne.-1) close(97)

!Al1 Call to 1D equivalent of dynredem (an,mois,jour,heure ?)
! -------------------------------------
       call dyn1dredem("restart1dyn.nc",                                    &
     &              plev,play,phi,phis,presnivs,                            &
     &              u,v,temp,q,omega2)

        CALL abort_gcm ('lmdz1d   ','The End  ',0)

        END SUBROUTINE old_lmdz1d

#include "old_1DUTILS_read_interp.h"
