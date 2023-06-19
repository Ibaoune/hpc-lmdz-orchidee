! Simulateur COSP : Cfmip Observation Simulator Package

! ISCCP, Radar (QuickBeam), Lidar et Parasol (ACTSIM), MISR, RTTOVS
!Idelkadi Abderrahmane Aout-Septembre 2009 First Version
!Idelkadi Abderrahmane Nov 2015 version v1.4.0

  subroutine phys_cosp2( itap,dtime,freq_cosp, &
                        ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP, &
                        ecrit_mth,ecrit_day,ecrit_hf, ok_all_xml, missing_val, &
                        Nptslmdz,Nlevlmdz,lon,lat, presnivs,overlaplmdz,sunlit, &
                        ref_liq,ref_ice,fracTerLic,u_wind,v_wind,phis,phi,ph,p,skt,t, &
                        sh,rh,tca,cca,mr_lsliq,mr_lsice,fl_lsrainI,fl_lssnowI, &
                        fl_ccrainI,fl_ccsnowI,mr_ozone,dtau_s,dem_s)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Inputs :
! itap,                                 !Increment de la physiq
! dtime,                                !Pas de temps physiq
! overlap,                              !Overlap type in SCOPS
! Npoints,                              !Nb de points de la grille physiq
! Nlevels,                              !Nb de niveaux verticaux
! Ncolumns,                             !Number of subcolumns
! lon,lat,                              !Longitudes et latitudes de la grille LMDZ
! ref_liq,ref_ice,                      !Rayons effectifs des particules liq et ice (en microm)
! fracTerLic,                           !Fraction terre a convertir en masque
! u_wind,v_wind,                        !Vents a 10m ???
! phi,                                  !Geopotentiel
! phis,                                 !Geopotentiel sol
! ph,                                   !pression pour chaque inter-couche
! p,                                    !Pression aux milieux des couches
! skt,t,                                !Temp au sol et temp 3D
! sh,                                   !Humidite specifique
! rh,                                   !Humidite relatif
! tca,                                  !Fraction nuageuse
! cca                                   !Fraction nuageuse convective
! mr_lsliq,                             !Liq Cloud water content
! mr_lsice,                             !Ice Cloud water content
! mr_ccliq,                             !Convective Cloud Liquid water content  
! mr_ccice,                             !Cloud ice water content
! fl_lsrain,                            !Large scale precipitation lic
! fl_lssnow,                            !Large scale precipitation ice
! fl_ccrain,                            !Convective precipitation lic
! fl_ccsnow,                            !Convective precipitation ice
! mr_ozone,                             !Concentration ozone (Kg/Kg)
! dem_s                                 !Cloud optical emissivity
! dtau_s               			!Cloud optical thickness
! emsfc_lw = 1.        			!Surface emissivity dans radlwsw.F90

!!! Outputs :
! calipso2D,                            !Lidar Low/heigh/Mean/Total-level Cloud Fraction
! calipso3D,                            !Lidar Cloud Fraction (532 nm)
! cfadlidar,                            !Lidar Scattering Ratio CFAD (532 nm)
! parasolrefl,                          !PARASOL-like mono-directional reflectance
! atb,                                  !Lidar Attenuated Total Backscatter (532 nm)
! betamol,                              !Lidar Molecular Backscatter (532 nm)
! cfaddbze,                             !Radar Reflectivity Factor CFAD (94 GHz)
! clcalipso2,                           !Cloud frequency of occurrence as seen by CALIPSO but not CloudSat
! dbze,                                 !Efective_reflectivity_factor
! cltlidarradar,                        !Lidar and Radar Total Cloud Fraction
! clMISR,                               !Cloud Fraction as Calculated by the MISR Simulator
! clisccp2,                             !Cloud Fraction as Calculated by the ISCCP Simulator
! boxtauisccp,                          !Optical Depth in Each Column as Calculated by the ISCCP Simulator
! boxptopisccp,                         !Cloud Top Pressure in Each Column as Calculated by the ISCCP Simulator
! tclisccp,                             !Total Cloud Fraction as Calculated by the ISCCP Simulator
! ctpisccp,                             !Mean Cloud Top Pressure as Calculated by the ISCCP Simulator
! tauisccp,                             !Mean Optical Depth as Calculated by the ISCCP Simulator
! albisccp,                             !Mean Cloud Albedo as Calculated by the ISCCP Simulator
! meantbisccp,                          !Mean all-sky 10.5 micron brightness temperature as calculated by the ISCCP Simulator
! meantbclrisccp                        !Mean clear-sky 10.5 micron brightness temperature as calculated by the ISCCP Simulator

!!! AI rajouter les nouvelles sorties
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! AI rajouter
!#include "cosp_defs.h" 
!  USE MOD_COSP_CONSTANTS
!  USE MOD_COSP_TYPES
!  USE MOD_COSP
  USE mod_phys_lmdz_para
  USE mod_grid_phy_lmdz
  use ioipsl
  use iophy
  use cosp_output_mod
  use cosp_output_write_mod
!  use MOD_COSP_Modis_Simulator, only : cosp_modis 
!  use mod_cosp_config, only : vgrid_zl,vgrid_zu,vgrid_z
  USE MOD_COSP_INTERFACE_v1p4, ONLY: cosp                   => cosp_interface_v1p4,      &
                                     cosp_gridbox,construct_cosp_vgrid,                  &
                                     construct_cosp_gridbox,                             &
                                     free_cosp_gridbox      => destroy_cosp_gridbox,     &
                                     free_cosp_sgradar      => destroy_cosp_sgradar,     &
                                     free_cosp_radarstats   => destroy_cosp_radarstats,  &
                                     free_cosp_sglidar      => destroy_cosp_sglidar,     &
                                     free_cosp_lidarstats   => destroy_cosp_lidarstats,  &
                                     free_cosp_isccp        => destroy_cosp_isccp,       &
                                     free_cosp_misr         => destroy_cosp_misr,        &
                                     free_cosp_rttov        => destroy_cosp_rttov,       &
                                     free_cosp_modis        => destroy_cosp_modis,       &
                                     free_cosp_vgrid        => destroy_cosp_vgrid,       &
                                     free_cosp_subgrid      => destroy_cosp_subgrid,     &
                                     construct_cosp_subgrid,cosp_config,cosp_subgrid,    &
                                     cosp_sglidar,cosp_lidarstats,                       &
                                     construct_cosp_lidarstats,construct_cosp_sglidar,   &
                                     cosp_isccp,construct_cosp_isccp,cosp_misr,          &
                                     construct_cosp_misr,cosp_rttov,construct_cosp_rttov,&
                                     cosp_sgradar,cosp_radarstats,                       &
                                     construct_cosp_radarstats,construct_cosp_sgradar,   &
                                     cosp_modis,construct_cosp_modis,                    &
                                     cosp_vgrid,I_CVCLIQ,I_LSCLIQ,I_CVCICE,I_LSCICE,     &
                                     I_LSRAIN,I_LSSNOW,I_LSGRPL,I_CVRAIN,I_CVSNOW

  use cosp_read_otputkeys

  IMPLICIT NONE

  ! Local variables
  character(len=64),PARAMETER  :: cosp_input_nl='cosp_input_nl.txt'
  character(len=64),PARAMETER  :: cosp_output_nl='cosp_output_nl.txt'
  integer, save :: isccp_topheight,isccp_topheight_direction,overlap
  integer,save  :: Ncolumns     ! Number of subcolumns in SCOPS
  integer, save :: Npoints      ! Number of gridpoints
!$OMP THREADPRIVATE(Npoints)
  integer, save :: Nlevels      ! Number of levels
  Integer :: Nptslmdz,Nlevlmdz ! Nb de points issus de physiq.F
  integer, save :: Nlr          ! Number of levels in statistical outputs
  integer, save :: Npoints_it   ! Max number of gridpoints to be processed in one iteration
!  integer :: i
  type(cosp_config),save :: cfg   ! Configuration options
!$OMP THREADPRIVATE(cfg)
  type(cosp_gridbox) :: gbx ! Gridbox information. Input for COSP
  type(cosp_subgrid) :: sgx     ! Subgrid outputs
  type(cosp_sgradar) :: sgradar ! Output from radar simulator
  type(cosp_sglidar) :: sglidar ! Output from lidar simulator
  type(cosp_isccp)   :: isccp   ! Output from ISCCP simulator
!! AI rajout modis
  type(cosp_modis)   :: modis   ! Output from MODIS simulator
!!
  type(cosp_misr)    :: misr    ! Output from MISR simulator
!! AI rajout rttovs
  type(cosp_rttov)   :: rttov   ! Output from RTTOV
!!
  type(cosp_vgrid)   :: vgrid   ! Information on vertical grid of stats
  type(cosp_radarstats) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator

  integer :: t0,t1,count_rate,count_max
  integer :: Nlon,Nlat
  real(wp),save :: radar_freq,k2,ZenAng,co2,ch4,n2o,co,emsfc_lw
  
!$OMP THREADPRIVATE(emsfc_lw)
  integer,dimension(RTTOV_MAX_CHANNELS),save :: Channels
  real(wp),dimension(RTTOV_MAX_CHANNELS),save :: Surfem
  integer, save :: surface_radar,use_mie_tables,use_gas_abs,do_ray,melt_lay
  integer, save :: Nprmts_max_hydro,Naero,Nprmts_max_aero,lidar_ice_type
  integer, save :: platform,satellite,Instrument,Nchannels
  logical, save :: use_vgrid,csat_vgrid,use_precipitation_fluxes,use_reff

! Declaration necessaires pour les sorties IOIPSL
  integer :: ii
  real    :: ecrit_day,ecrit_hf,ecrit_mth, missing_val
  logical :: ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP, ok_all_xml

  logical, save :: debut_cosp=.true.
!$OMP THREADPRIVATE(debut_cosp)

  logical, save :: first_write=.true.
!$OMP THREADPRIVATE(first_write)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Input variables from LMDZ-GCM
  integer                         :: overlaplmdz   !  overlap type: 1=max, 2=rand, 3=max/rand ! cosp input (output lmdz)
  real,dimension(Nptslmdz,Nlevlmdz) :: height,phi,p,ph,T,sh,rh,tca,cca,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice, & 
                                     fl_lsrain,fl_lssnow,fl_ccrain,fl_ccsnow,fl_lsgrpl, &
                                     zlev,zlev_half,mr_ozone,radliq,radice,dtau_s,dem_s,ref_liq,ref_ice
  real,dimension(Nptslmdz,Nlevlmdz) ::  fl_lsrainI,fl_lssnowI,fl_ccrainI,fl_ccsnowI
  real,dimension(Nptslmdz)        :: lon,lat,skt,fracTerLic,u_wind,v_wind,phis,sunlit         
  real,dimension(Nlevlmdz)        :: presnivs
  integer                         :: itap,k,ip
  real                            :: dtime,freq_cosp
  real,dimension(2)               :: time_bnds

  double precision                            :: d_dtime
  double precision,dimension(2)               :: d_time_bnds

   namelist/COSP_INPUT/overlap,isccp_topheight,isccp_topheight_direction, &
              npoints_it,ncolumns,use_vgrid,Nlr,csat_vgrid, &
              radar_freq,surface_radar,use_mie_tables, &
              use_gas_abs,do_ray,melt_lay,k2,Nprmts_max_hydro,Naero,Nprmts_max_aero, &
              lidar_ice_type,use_precipitation_fluxes,use_reff, &
              platform,satellite,Instrument,Nchannels, &
              Channels,Surfem,ZenAng,co2,ch4,n2o,co

!---------------- End of declaration of variables --------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Read namelist with COSP inputs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 print*,'Entree phys_cosp2'
 if (debut_cosp) then
  NPoints=Nptslmdz
  Nlevels=Nlevlmdz
  
! Lecture du namelist input 
!  CALL read_cosp_input
   IF (is_master) THEN
      OPEN(10,file=cosp_input_nl,status='old')
      READ(10,nml=cosp_input)
      CLOSE(10)
   ENDIF
!$OMP BARRIER 
    CALL bcast(overlap)
    CALL bcast(isccp_topheight)
    CALL bcast(isccp_topheight_direction)
    CALL bcast(npoints_it)
    CALL bcast(ncolumns)
    CALL bcast(use_vgrid)
    CALL bcast(Nlr)
    CALL bcast(csat_vgrid)
    CALL bcast(radar_freq)
    CALL bcast(surface_radar)
    CALL bcast(use_mie_tables)
    CALL bcast(use_gas_abs)
    CALL bcast(do_ray)
    CALL bcast(melt_lay)
    CALL bcast(k2)
    CALL bcast(Nprmts_max_hydro)
    CALL bcast(Naero)
    CALL bcast(Nprmts_max_aero)
    CALL bcast(lidar_ice_type)
    CALL bcast(use_precipitation_fluxes)
    CALL bcast(use_reff)
    CALL bcast(platform)
    CALL bcast(satellite)
    CALL bcast(Instrument)
    CALL bcast(Nchannels)
    CALL bcast(Channels)
    CALL bcast(Surfem)
    CALL bcast(ZenAng)
    CALL bcast(co2)
    CALL bcast(ch4)
    CALL bcast(n2o)
    CALL bcast(co)

    print*,'ok read  cosp_input_nl'

! Clefs Outputs initialisation
#ifdef CPP_XIOS
  call cosp_outputkeys_init(cfg)
#else
   call read_cosp_output_nl(itap,cosp_output_nl,cfg)
#endif

  print*,' Cles des differents simulateurs cosp a itap :',itap
  print*,'Lradar_sim,Llidar_sim,Lisccp_sim,Lmisr_sim,Lmodis_sim,Lrttov_sim,Lstats', &
          cfg%Lradar_sim,cfg%Llidar_sim,cfg%Lisccp_sim,cfg%Lmisr_sim,cfg%Lmodis_sim, &
          cfg%Lrttov_sim,cfg%Lstats

    if (overlaplmdz.ne.overlap) then
       print*,'Attention overlaplmdz different de overlap lu dans namelist '
    endif
   print*,'Fin lecture Namelists, debut_cosp =',debut_cosp

  endif ! debut_cosp

!!! Ici on modifie les cles logiques pour les outputs selon les champs actives dans les .xml
  if ((itap.gt.1).and.(first_write))then
#ifdef CPP_XIOS
    call read_xiosfieldactive(cfg)
#endif
    first_write=.false.

    print*,' Cles des differents simulateurs cosp a itap :',itap
    print*,'Lradar_sim,Llidar_sim,Lisccp_sim,Lmisr_sim,Lmodis_sim,Lrttov_sim,Lstats', &
          cfg%Lradar_sim,cfg%Llidar_sim,cfg%Lisccp_sim,cfg%Lmisr_sim,cfg%Lmodis_sim, &
          cfg%Lrttov_sim,cfg%Lstats
  endif

  time_bnds(1) = dtime-dtime/2.
  time_bnds(2) = dtime+dtime/2.

  d_time_bnds=time_bnds
  d_dtime=dtime

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Allocate memory for gridbox type
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AI mars 2017
!        print *, 'Allocating memory for gridbox type...'

! Surafce emissivity
        emsfc_lw = 1.

        call construct_cosp_gridbox(d_dtime,d_time_bnds,radar_freq,surface_radar,use_mie_tables,use_gas_abs, &
                                    do_ray,melt_lay,k2, &
                                    Npoints,Nlevels,Ncolumns,N_HYDRO,Nprmts_max_hydro,Naero,Nprmts_max_aero,Npoints_it, &
                                    lidar_ice_type,isccp_topheight,isccp_topheight_direction,overlap,emsfc_lw, &
                                    use_precipitation_fluxes,use_reff, &
                                    Platform,Satellite,Instrument,Nchannels,ZenAng, &
                                    channels(1:Nchannels),surfem(1:Nchannels),co2,ch4,n2o,co,gbx)
        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Here code to populate input structure
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!        print *, 'Populating input structure...'
        gbx%longitude = lon
        gbx%latitude = lat

        gbx%p = p !
        gbx%ph = ph
        gbx%zlev = phi/9.81

        zlev_half(:,1) = phis(:)/9.81
        do k = 2, Nlevels
          do ip = 1, Npoints
           zlev_half(ip,k) = phi(ip,k)/9.81 + &
               (phi(ip,k)-phi(ip,k-1))/9.81 * (ph(ip,k)-p(ip,k)) / (p(ip,k)-p(ip,k-1))
          enddo
        enddo
        gbx%zlev_half = zlev_half

        gbx%T = T
        gbx%q = rh*100.
        gbx%sh = sh
! On ne veut pas que cosp distingue les nuages stratiformes et convectifs
! on passe les contenus totaux (conv+strat)
        gbx%cca = 0. !convective_cloud_amount (1)
        gbx%tca = tca ! total_cloud_amount (1)
        gbx%psfc = ph(:,1) !pression de surface
        gbx%skt  = skt !Skin temperature (K)

        do ip = 1, Npoints
          if (fracTerLic(ip).ge.0.5) then
             gbx%land(ip) = 1.
          else
             gbx%land(ip) = 0.
          endif
        enddo
        gbx%mr_ozone  = mr_ozone !mass_fraction_of_ozone_in_air (kg/kg)
! A voir l equivalent LMDZ (u10m et v10m)
        gbx%u_wind  = u_wind !eastward_wind (m s-1)
        gbx%v_wind  = v_wind !northward_wind

! sunlit calcule a partir de la fraction d ensoleillement par jour
!      do ip = 1, Npoints
!        if (sunlit(ip).le.0.) then
!           gbx%sunlit(ip)=0.
!        else
!           gbx%sunlit(ip)=1.
!        endif
!      enddo
       gbx%sunlit=sunlit

! A voir l equivalent LMDZ
  mr_ccliq = 0.0
  mr_ccice = 0.0
        gbx%mr_hydro(:,:,I_LSCLIQ) = mr_lsliq !mixing_ratio_large_scale_cloud_liquid (kg/kg)
        gbx%mr_hydro(:,:,I_LSCICE) = mr_lsice !mixing_ratio_large_scale_cloud_ic
        gbx%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq !mixing_ratio_convective_cloud_liquid
        gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice !mixing_ratio_convective_cloud_ice
! A revoir
        fl_lsrain = fl_lsrainI + fl_ccrainI
        fl_lssnow = fl_lssnowI + fl_ccsnowI
        gbx%rain_ls = fl_lsrain !flux_large_scale_cloud_rain (kg m^-2 s^-1)
        gbx%snow_ls = fl_lssnow !flux_large_scale_cloud_snow
!  A voir l equivalent LMDZ
        fl_lsgrpl=0.
        fl_ccsnow = 0.
        fl_ccrain = 0.
        gbx%grpl_ls = fl_lsgrpl  !flux_large_scale_cloud_graupel
        gbx%rain_cv = fl_ccrain  !flux_convective_cloud_rain
        gbx%snow_cv = fl_ccsnow  !flux_convective_cloud_snow

     gbx%Reff(:,:,I_LSCLIQ) = ref_liq*1e-6
     gbx%Reff(:,:,I_LSCICE) = ref_ice*1e-6
!! AI A revoir
     gbx%Reff(:,:,I_CVCLIQ) = ref_liq*1e-6
     gbx%Reff(:,:,I_CVCICE) = ref_ice*1e-6

        ! ISCCP simulator
        gbx%dtau_s   = dtau_s
        gbx%dtau_c   = 0.
        gbx%dem_s    = dem_s
        gbx%dem_c   = 0.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Define new vertical grid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        print *, 'Defining new vertical grid...'
        call construct_cosp_vgrid(gbx,Nlr,use_vgrid,csat_vgrid,vgrid)
        print*,'Nlvgrid,use_vgrid,csat_vgrid ',Nlr,use_vgrid,csat_vgrid
        print*,'vgrid%z, vgrid%mz ',vgrid%z, vgrid%mz
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! Allocate memory for other types
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        print *, 'Allocating memory for other types...'
        call construct_cosp_subgrid(Npoints, Ncolumns, Nlevels, sgx)
        if (cfg%Lradar_sim) call construct_cosp_sgradar(Npoints,Ncolumns,Nlevels,N_HYDRO,sgradar) 
        if (cfg%Lradar_sim) call construct_cosp_radarstats(Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,stradar)
        if (cfg%Llidar_sim) call construct_cosp_sglidar(Npoints,Ncolumns,Nlevels,N_HYDRO,PARASOL_NREFL,sglidar)
        if (cfg%Llidar_sim) call construct_cosp_lidarstats(Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,stlidar)
        if (cfg%Lisccp_sim) call construct_cosp_isccp(Npoints,Ncolumns,Nlevels,isccp)
!! AI rajout
        if (cfg%Lmodis_sim) call construct_cosp_modis(Npoints,modis)
!!
        if (cfg%Lmisr_sim) call construct_cosp_misr(Npoints,misr)
!        call construct_cosp_rttov(cfg,Npoints,Nchannels,rttov)

    if (debut_cosp) then
      !$OMP MASTER
        print *, ' Open outpts files and define axis'
        call cosp_output_open(Nlevlmdz, Ncolumns, presnivs, dtime, freq_cosp, &
                              ok_mensuelCOSP, ok_journeCOSP, ok_hfCOSP, ok_all_xml, &
                              ecrit_mth, ecrit_day, ecrit_hf, use_vgrid, vgrid)
      !$OMP END MASTER
      !$OMP BARRIER
       debut_cosp=.false.
      endif ! debut_cosp

!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Call simulator
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      print*,'call simulateur'
       call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,   &
                                 isccp,misr,modis,rttov,stradar,stlidar)


!!!!!!!!!!!!!!!!!! Ecreture des sorties Cosp !!!!!!!!!!!!!!r!!!!!!:!!!!!

       print *, 'Calling write output'
        call cosp_output_write(Nlevlmdz, Npoints, Ncolumns, itap, dtime, freq_COSP, missing_val, &
                               cfg, gbx, vgrid, sglidar, sgradar, stlidar, stradar, & 
                               isccp, misr, modis)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Deallocate memory in derived types
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      print *, 'Deallocating memory...'
     call free_cosp_gridbox(gbx)
      print *, 'free_cosp_gridbox'
     call free_cosp_subgrid(sgx)
      print *, 'ffree_cosp_subgrid'
     call free_cosp_vgrid(vgrid)
      print *, 'free_cosp_vgrid'
     if (cfg%Lradar_sim) call free_cosp_sgradar(sgradar)
     if (cfg%Lradar_sim) call free_cosp_radarstats(stradar)
     if (cfg%Llidar_sim) call free_cosp_sglidar(sglidar)
     if (cfg%Llidar_sim) call free_cosp_lidarstats(stlidar)
     if (cfg%Lisccp_sim) call free_cosp_isccp(isccp)
     if (cfg%Lmisr_sim)  call free_cosp_misr(misr)
     if (cfg%Lmodis_sim) call free_cosp_modis(modis)
     if (cfg%Lrttov_sim) call free_cosp_rttov(rttov)
     print *, 'End phys_cosp2'
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Time in s. Only for testing purposes
!  call system_clock(t1,count_rate,count_max)
!  print *,(t1-t0)*1.0/count_rate
 
end subroutine phys_cosp2
