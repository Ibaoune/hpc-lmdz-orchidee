!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Nouveau code d'interface entre LMDZ et COSPv2 (version 2)
! L'ancienne interface s'appelait "phys_cosp2" et avait ete concue pour COSPv1.4.
! Dans cette nouvelle version de COSP, le code a ete restructure pour optimiser les calculs
! des differents simulateurs et pour proposer de nouvelles fonctionnalites (par exemple,
! intervenir sur les profils sous-maille, ou subcolumns, donnes en entre a COSP afin que
! leur definition soit coherente avec les parametrisations du modele hote).
! Cette version de COSP propose aussi de nombreux nouveaux diagnostics, notamment pour
! le simulateur lidar (diagnostics CALIPSO-OPAQ, lidar sol 532nm et lidar ATLID 355nm).
!
! Interface reecrite par R.Guzman (01/2019), a partir de l'interface initiale concue
! et ecrite par A.Idelkadi 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  subroutine phys_cosp2( itap,dtime,freq_cosp, &
  subroutine lmdz_cosp_interface(itap, dtime, freq_cosp, ok_mensuelCOSP, ok_journeCOSP,   &
                         ok_hfCOSP, ecrit_mth, ecrit_day, ecrit_hf, ok_all_xml,   &
                         missing_val, Nptslmdz, Nlevlmdz, lon, lat, presnivs,     &
                         overlaplmdz, sunlit, ref_liq, ref_ice, fracTerLic,       &
                         u_wind, v_wind, phis, phi, ph, p, skt, t, sh, rh,        &
                         tca, cca, mr_lsliq, mr_lsice, fl_lsrainI, fl_lssnowI,    &
                         fl_ccrainI, fl_ccsnowI, mr_ozone, dtau_s, dem_s)

!--------------  Inputs  ---------------
! itap,                                 !Increment de la physiq
! dtime,                                !Pas de temps physiq
! overlaplmdz,                          !Type Overlap venant de LMDZ
! Npoints,                              !Nb de points de la grille physiq
! Nlevels,                              !Nb de niveaux verticaux
! Ncolumns,                             !Number of subcolumns
! lon,lat,                              !Longitudes et latitudes de la grille LMDZ
! ref_liq, ref_ice,                     !Rayons effectifs des particules liq et ice (en micron)
! fracTerLic,                           !Fraction terre a convertir en masque
! u_wind, v_wind,                       !Vents a 10m ???
! phi,                                  !Geopotentiel
! phis,                                 !Geopotentiel sol
! ph,                                   !pression pour chaque inter-couche
! p,                                    !Pression aux milieux des couches
! skt, t,                               !Temp au sol et temp 3D
! sh,                                   !Humidite specifique
! rh,                                   !Humidite relative
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


!--------------  Outputs  --------------
! La liste complete des diagnostics de sortie (observables simulees) que l'on peut
! avoir avec COSPv2 se trouve au debut du fichier : cosp_read_otputkeys.F90


!!! Modules specifiques a l'interface LMDZ-COSP
  use mod_phys_lmdz_para
  use mod_grid_phy_lmdz
  use ioipsl
  use iophy
  use lmdz_cosp_output_mod
  use lmdz_cosp_output_write_mod
  use lmdz_cosp_read_outputkeys
  use lmdz_cosp_subsample_and_optics_mod, only : subsample_and_optics
  use lmdz_cosp_construct_destroy_mod

!!! Modules faisant partie du code source de COSPv2
  use cosp_kinds,		  only: wp                         
  use MOD_COSP_CONFIG,     		  only: N_HYDRO,RTTOV_MAX_CHANNELS, &
      			   	 		niv_sorties, vgrid_z_in
  use mod_quickbeam_optics,		  only: size_distribution,hydro_class_init, &
                                 		quickbeam_optics_init
  use quickbeam,           		  only: radar_cfg
  use mod_cosp,            		  only: cosp_init,cosp_optical_inputs, &
      						cosp_column_inputs,cosp_outputs, &
						cosp_simulator


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Declaration des variables 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  IMPLICIT NONE

  ! Local variables
  character(len=64),PARAMETER	:: cosp_input_nl  = 'cospv2_input_nl.txt'
  character(len=64),PARAMETER  	:: cosp_output_nl = 'cospv2_output_nl.txt'

  integer, save 	       	:: isccp_topheight, isccp_topheight_direction, overlap
  integer, save  	       	:: Ncolumns    	    	! Number of subcolumns in SCOPS
  integer, save 	       	:: Npoints            	! Number of gridpoints
!$OMP THREADPRIVATE(Npoints)
  integer, save			:: Nlevels		! Number of model vertical levels
  integer  			:: Nptslmdz, Nlevlmdz  	! Nb de points issus de physiq.F
  integer, save 		:: Npoints_it   	! Max number of gridpoints to be 
  	   			   			! processed in one iteration
  type(cosp_config), save 	:: cfg          	! Variable qui contient les cles
  		     		   			! logiques des simulateurs et des
							! diagnostics, definie dans:
							! lmdz_cosp_construct_destroy_mod
!$OMP THREADPRIVATE(cfg)

  integer			:: t0, t1, count_rate, count_max
  real(wp), save 		:: cloudsat_radar_freq, cloudsat_k2, rttov_ZenAng, co2, &
  	    			   ch4, n2o, co, emsfc_lw
!$OMP THREADPRIVATE(emsfc_lw)

  integer, dimension(RTTOV_MAX_CHANNELS), save	 :: rttov_Channels
  real(wp), dimension(RTTOV_MAX_CHANNELS), save  :: rttov_Surfem
  integer, save 			   	 :: surface_radar, use_mie_tables, &
  	   					    cloudsat_use_gas_abs, cloudsat_do_ray, &
						    melt_lay
  integer, save 				 :: lidar_ice_type
  integer, save 				 :: rttov_platform, rttov_satellite, &
  	   					    rttov_Instrument, rttov_Nchannels
  logical, save 				 :: use_vgrid_in, csat_vgrid_in, &
  	   					    use_precipitation_fluxes

! Declaration necessaires pour les sorties IOIPSL
  real        	:: ecrit_day, ecrit_hf, ecrit_mth, missing_val
  logical 	:: ok_mensuelCOSP, ok_journeCOSP, ok_hfCOSP, ok_all_xml
  logical, save :: debut_cosp=.true.
!$OMP THREADPRIVATE(debut_cosp)

  logical, save :: first_write=.true.
!$OMP THREADPRIVATE(first_write)

  integer, save :: cosp_init_flag = 0
!$OMP THREADPRIVATE(cosp_init_flag)


!-----------------------------  Input variables from LMDZ-GCM  -------------------------------
  integer                             	     :: overlaplmdz   ! overlap type: 1=max, 
  					       		      ! 2=rand, 3=max/rand
  real, dimension(Nptslmdz,Nlevlmdz) 	     :: phi, p, ph, T, sh, rh, tca, cca, mr_lsliq,   &
  					        mr_lsice, mr_ccliq, mr_ccice, fl_lsrain,     & 
                                       	        fl_lssnow, fl_ccrain, fl_ccsnow, fl_lsgrpl,  &
                                       	        zlev, zlev_half, mr_ozone, radliq, radice,   &
					        dtau_s, dem_s, dtau_c, dem_c, ref_liq, ref_ice
  real, dimension(Nptslmdz,Nlevlmdz) 	     :: fl_lsrainI, fl_lssnowI, fl_ccrainI, fl_ccsnowI
  real, dimension(Nptslmdz)        	     :: lon, lat, skt, fracTerLic, u_wind, v_wind, &
  					        phis, sunlit
  real, dimension(Nptslmdz)        	     :: land ! variables intermediaire pour masque TerLic
  real, dimension(Nlevlmdz)        	     :: presnivs
  integer                         	     :: itap, k, ip
  real                            	     :: dtime, freq_cosp
  real, dimension(2)               	     :: time_bnds

  double precision                           :: d_dtime
  double precision, dimension(2)             :: d_time_bnds


  ! ######################################################################################
  ! Declarations specific to COSP2
  ! ######################################################################################

  ! Local variables
  logical :: &
       Lsingle     = .true.,  & ! True if using MMF_v3_single_moment CLOUDSAT
       		     	      	! microphysical scheme (default)
       Ldouble     = .false.    ! True if using MMF_v3.5_two_moment CLOUDSAT 
       		     		! microphysical scheme
  type(size_distribution), save   	     :: sd            ! Hydrometeor description
!$OMP THREADPRIVATE(sd)
  type(radar_cfg), save			     :: rcfg_cloudsat ! Radar configuration
!$OMP THREADPRIVATE(rcfg_cloudsat)
  real, dimension(Nptslmdz,Nlevlmdz,N_HYDRO) :: Reff	      ! Liquid and Ice particles
  					     		      ! effective radius
  type(cosp_outputs)			     :: cospOUT       ! COSP simulator outputs
  type(cosp_optical_inputs) 		     :: cospIN        ! COSP optical (or derived?)
  					     		      ! fields needed by simulators
  type(cosp_column_inputs) 		     :: cospstateIN   ! COSP model fields needed
  					     		      ! by simulators
  character(len=256), dimension(100)	     :: cosp_status
  character(len=64), save      		     :: cloudsat_micro_scheme

  ! Indices to address arrays of LS and CONV hydrometeors
  integer,parameter :: &
       I_LSCLIQ = 1, & ! Large-scale (stratiform) liquid
       I_LSCICE = 2, & ! Large-scale (stratiform) ice
       I_LSRAIN = 3, & ! Large-scale (stratiform) rain
       I_LSSNOW = 4, & ! Large-scale (stratiform) snow
       I_CVCLIQ = 5, & ! Convective liquid
       I_CVCICE = 6, & ! Convective ice
       I_CVRAIN = 7, & ! Convective rain
       I_CVSNOW = 8, & ! Convective snow
       I_LSGRPL = 9    ! Large-scale (stratiform) groupel

! Parametres qui sont lus a partir du fichier "cosp_input_nl.txt"
   namelist/COSP_INPUT/overlap, isccp_topheight, isccp_topheight_direction,	  &
              npoints_it, ncolumns, use_vgrid_in, csat_vgrid_in,		  &
              cloudsat_radar_freq, surface_radar, use_mie_tables,		  &
              cloudsat_use_gas_abs, cloudsat_do_ray, melt_lay, cloudsat_k2, 	  &
              cloudsat_micro_scheme, lidar_ice_type, use_precipitation_fluxes, 	  &
              rttov_platform, rttov_satellite, rttov_Instrument, rttov_Nchannels, &
              rttov_Channels, rttov_Surfem, rttov_ZenAng, co2, ch4, n2o, co

!------------------------   Fin declaration des variables   ------------------------


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 1) Lecture du fichier "cosp_input_nl.txt", parametres d'entree pour COSP
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 print*,'Entree lmdz_cosp_interface' !phys_cosp2'
 if (debut_cosp) then
  NPoints=Nptslmdz
  Nlevels=Nlevlmdz
! Surface emissivity
        emsfc_lw = 1.
  
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
    CALL bcast(use_vgrid_in)
    CALL bcast(csat_vgrid_in)
    CALL bcast(cloudsat_radar_freq)
    CALL bcast(surface_radar)
    CALL bcast(cloudsat_use_gas_abs)
    CALL bcast(cloudsat_do_ray)
    CALL bcast(cloudsat_k2)
    CALL bcast(lidar_ice_type)
    CALL bcast(use_precipitation_fluxes)
    CALL bcast(rttov_platform)
    CALL bcast(rttov_satellite)
    CALL bcast(rttov_Instrument)
    CALL bcast(rttov_Nchannels)
    CALL bcast(rttov_Channels)
    CALL bcast(rttov_Surfem)
    CALL bcast(rttov_ZenAng)
    CALL bcast(co2)
    CALL bcast(ch4)
    CALL bcast(n2o)
    CALL bcast(co)
    CALL bcast(cloudsat_micro_scheme)

    print*,'ok read  cosp_input_nl'


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 2) Initialisation de COSPv2
!
! Il y a 2 options possibles pour ecrire les sorties : XIOS ou sorties standard
!
! Si le modele a ete compile AVEC l'option xios, le programme passe par les
! bouts de code delimites par "#ifdef CPP_XIOS" et "#else". Dans ce cas,
! l'initialisation de COSP se fait au deuxieme appel de cette interface.
!
! Si le modele a ete compile SANS l'option xios, le programme passe par les
! bouts de code delimites par "#else" et "#endif". Dans ce cas,
! l'initialisation de COSP se fait au premier appel de cette interface.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Clefs Outputs initialisation
#ifdef CPP_XIOS
  call cosp_outputkeys_init(cfg)
#else
   call read_cosp_output_nl(itap,cosp_output_nl,cfg)
#endif

  print*,' Cles des differents simulateurs cosp a itap :',itap
print*,'cfg%Lcloudsat, cfg%Lcalipso, cfg%LgrLidar532, cfg%Latlid, cfg%Lparasol, &
	cfg%Lisccp, cfg%Lmisr, cfg%Lmodis, cfg%Lrttov', &
        cfg%Lcloudsat, cfg%Lcalipso, cfg%LgrLidar532, cfg%Latlid, cfg%Lparasol, &
	cfg%Lisccp, cfg%Lmisr, cfg%Lmodis, cfg%Lrttov

    if (overlaplmdz.ne.overlap) then
       print*,'Attention overlaplmdz different de overlap lu dans namelist '
    endif

#ifdef CPP_XIOS
   print*,'On passe par ifdef CPP_XIOS'
#else
if (cosp_init_flag .eq. 0) then

    ! Initialize the distributional parameters for hydrometeors in radar simulator.
    ! In COSPv1.4, this was declared in cosp_defs.f.
    if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
       ldouble = .true. 
       lsingle = .false.
    endif
    call hydro_class_init(lsingle,ldouble,sd)
    call quickbeam_optics_init()

 print*,' just before call COSP_INIT, cosp_init_flag =', cosp_init_flag
  call COSP_INIT(cfg%Lisccp, cfg%Lmodis, cfg%Lmisr, cfg%Lcloudsat, cfg%Lcalipso, &
       		 cfg%LgrLidar532, cfg%Latlid, cfg%Lparasol, cfg%Lrttov,          &
       		 cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs,         &
       		 cloudsat_do_ray, isccp_topheight, isccp_topheight_direction,    &
		 surface_radar, rcfg_cloudsat, use_vgrid_in, csat_vgrid_in,      &
       		 niv_sorties, Nlevels, cloudsat_micro_scheme)
cosp_init_flag = 1
 print*,' just after call COSP_INIT, cosp_init_flag =', cosp_init_flag
endif
#endif

   print*,'Fin lecture Namelists, debut_cosp =',debut_cosp

  endif ! debut_cosp


!!! Ici on modifie les cles logiques pour les outputs selon les champs actives dans les .xml
  if ((itap.gt.1).and.(first_write))then
#ifdef CPP_XIOS
    call read_xiosfieldactive(cfg)
#endif
    first_write=.false.

if (cosp_init_flag .eq. 0) then

    ! Initialize the distributional parameters for hydrometeors in radar simulator.
    ! In COSPv1.4, this was declared in cosp_defs.f.
    if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
       ldouble = .true. 
       lsingle = .false.
    endif
    call hydro_class_init(lsingle,ldouble,sd)
    call quickbeam_optics_init()

 print*,' just before call COSP_INIT, cosp_init_flag =', cosp_init_flag
  call COSP_INIT(cfg%Lisccp, cfg%Lmodis, cfg%Lmisr, cfg%Lcloudsat, cfg%Lcalipso, &
       		 cfg%LgrLidar532, cfg%Latlid, cfg%Lparasol, cfg%Lrttov,          &
       		 cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs,         &
       		 cloudsat_do_ray, isccp_topheight, isccp_topheight_direction,    &
		 surface_radar, rcfg_cloudsat, use_vgrid_in, csat_vgrid_in,      &
       		 niv_sorties, Nlevels, cloudsat_micro_scheme)
cosp_init_flag = 1
 print*,' just after call COSP_INIT, cosp_init_flag =', cosp_init_flag
endif ! cosp_init_flag


print*,' Cles des differents simulateurs cosp a itap :',itap
print*,'cfg%Lcloudsat, cfg%Lcalipso, cfg%LgrLidar532, cfg%Latlid, cfg%Lparasol, &
	cfg%Lisccp, cfg%Lmisr, cfg%Lmodis, cfg%Lrttov', &
        cfg%Lcloudsat, cfg%Lcalipso, cfg%LgrLidar532, cfg%Latlid, cfg%Lparasol, &
	cfg%Lisccp, cfg%Lmisr, cfg%Lmodis, cfg%Lrttov

  endif !(itap.gt.1).and.(first_write)

  time_bnds(1) = dtime-dtime/2.
  time_bnds(2) = dtime+dtime/2.

  d_time_bnds=time_bnds
  d_dtime=dtime

!-------------------------   Fin initialisation de COSP   --------------------------


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 3) Calculs des champs d'entree COSP a partir des variables LMDZ
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! 0) Create ptop/ztop for gbx%pf and gbx%zlev are for the the interface, 
    !    also reverse CAM height/pressure values for input into CSOP
    !    CAM state%pint from top to surface, COSP wants surface to top.

! 0) Altitudes du modele calculees a partir de la variable geopotentiel phi et phis    
        zlev = phi/9.81

        zlev_half(:,1) = phis(:)/9.81
        do k = 2, Nlevels
          do ip = 1, Npoints
           zlev_half(ip,k) = phi(ip,k)/9.81 + &
               (phi(ip,k)-phi(ip,k-1))/9.81 * (ph(ip,k)-p(ip,k)) / (p(ip,k)-p(ip,k-1))
          enddo
        enddo

! 1) Quantite de nuages (couverture?), convectif (=0) et total
        cca = 0._wp ! convective_cloud_amount (1)
        tca = tca   ! total_cloud_amount (1)

! 2) Humidite relative est donnee tel quel (variable rh)

! 3) Masque terre/mer a partir de la variable fracTerLic
        do ip = 1, Npoints
          if (fracTerLic(ip).ge.0.5) then
             land(ip) = 1.
          else
             land(ip) = 0.
          endif
        enddo


! A voir l equivalent LMDZ
  mr_ccliq = 0.0
  mr_ccice = 0.0
!!!        gbx%mr_hydro(:,:,I_LSCLIQ) = mr_lsliq !mixing_ratio_large_scale_cloud_liquid (kg/kg)
!!!        gbx%mr_hydro(:,:,I_LSCICE) = mr_lsice !mixing_ratio_large_scale_cloud_ic
!!!        gbx%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq !mixing_ratio_convective_cloud_liquid
!!!        gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice !mixing_ratio_convective_cloud_ice
! A revoir
        fl_lsrain = fl_lsrainI + fl_ccrainI
        fl_lssnow = fl_lssnowI + fl_ccsnowI
!!!        gbx%rain_ls = fl_lsrain !flux_large_scale_cloud_rain (kg m^-2 s^-1)
!!!        gbx%snow_ls = fl_lssnow !flux_large_scale_cloud_snow
!  A voir l equivalent LMDZ
        fl_lsgrpl = 0.
        fl_ccsnow = 0.
        fl_ccrain = 0.
!!!        gbx%grpl_ls = fl_lsgrpl  !flux_large_scale_cloud_graupel
!!!        gbx%rain_cv = fl_ccrain  !flux_convective_cloud_rain
!!!        gbx%snow_cv = fl_ccsnow  !flux_convective_cloud_snow

        ! ISCCP simulator
        dtau_c  = 0.
        dem_c   = 0.

! note: reff_cosp dimensions should be same as cosp (reff_cosp has 9 hydrometeor dimension)
     Reff(1:Npoints,1:Nlevels,1:N_HYDRO) = 0.
     Reff(:,:,I_LSCLIQ) = ref_liq*1e-6
     Reff(:,:,I_LSCICE) = ref_ice*1e-6
     Reff(:,:,I_CVCLIQ) = ref_liq*1e-6
     Reff(:,:,I_CVCICE) = ref_ice*1e-6


if (cosp_init_flag .eq. 1) then      ! cosp_init_flag = 1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 4a) On construit la variable cospOUT qui contient tous les diagnostics de sortie.
! Elle sera remplie lors de l'appel du simulateur COSP
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call construct_cosp_outputs(cfg,Npoints,Ncolumns,Nlevels,niv_sorties,0,cospOUT)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 4b) On construit la variable cospstateIN que l'on va remplir avec les champs LMDZ
! Les champ verticaux doivent etre donnes a l'envers, c-a-d : (Nlevels:1) = (TOA:SFC)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call construct_cospstateIN(Npoints,Nlevels,0,cospstateIN)

    cospstateIN%lat                 		        = lat(1:Npoints)
    cospstateIN%lon                     		= lon(1:Npoints) 
    cospstateIN%at              	      		= t(1:Npoints,Nlevels:1:-1) 
    cospstateIN%qv              	     		= sh(1:Npoints,Nlevels:1:-1)
    cospstateIN%o3              	      		= mr_ozone(1:Npoints,Nlevels:1:-1)  
    cospstateIN%sunlit          	      		= sunlit(1:Npoints)
    cospstateIN%skt             	       		= skt(1:Npoints)
    cospstateIN%land            	     		= land(1:Npoints)
    cospstateIN%surfelev            	     		= zlev_half(1:Npoints,1)
    cospstateIN%pfull           	      		= p(1:Npoints,Nlevels:1:-1)
    cospstateIN%phalf(1:Npoints,1)	                = 0._wp
    cospstateIN%phalf(1:Npoints,2:Nlevels+1)	        = ph(1:Npoints,Nlevels:1:-1)  
    cospstateIN%hgt_matrix                  		= zlev(1:Npoints,Nlevels:1:-1) 
    cospstateIN%hgt_matrix_half(1:Npoints,Nlevels+1) 	= 0._wp
    cospstateIN%hgt_matrix_half(1:Npoints,1:Nlevels) 	= zlev_half(1:Npoints,Nlevels:1:-1) 


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 4c) On construit la variable cospIN qui contient les proprietes optiques subcolumn
! pour COSP. Elle sera essentiellement remplie dans la subroutine subsample_and_optics
! ou sont appeles SCOPS, PREC_SCOPS et les subroutines qui calculent les signaux
! simules pour chaque simulateur actif.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call construct_cospIN(cfg,Npoints,Ncolumns,Nlevels,cospIN)
    cospIN%emsfc_lw  = emsfc_lw
    if (cfg%Lcloudsat) cospIN%rcfg_cloudsat = rcfg_cloudsat


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 5) Appel de subsample_and_optics : Les champs verticaux doivent etre donnes a
! l'envers comme pour le remplissage de cospstateIN, c-a-d : (Nlevels:1) = (TOA:SFC)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call subsample_and_optics(cfg, Npoints, Nlevels, Ncolumns, N_HYDRO,overlap,	&
        		      use_precipitation_fluxes, lidar_ice_type, sd, &
			      tca(1:Npoints,Nlevels:1:-1), cca(1:Npoints,Nlevels:1:-1), &
           		      fl_lsrain(1:Npoints,Nlevels:1:-1), &
			      fl_lssnow(1:Npoints,Nlevels:1:-1), &
			      fl_lsgrpl(1:Npoints,Nlevels:1:-1), &
			      fl_ccrain(1:Npoints,Nlevels:1:-1), &
          		      fl_ccsnow(1:Npoints,Nlevels:1:-1), &
			      mr_lsliq(1:Npoints,Nlevels:1:-1), &
			      mr_lsice(1:Npoints,Nlevels:1:-1), &
			      mr_ccliq(1:Npoints,Nlevels:1:-1), &
			      mr_ccice(1:Npoints,Nlevels:1:-1), &
             		      Reff(1:Npoints,Nlevels:1:-1,:), &
			      dtau_c(1:Npoints,Nlevels:1:-1), &
			      dtau_s(1:Npoints,Nlevels:1:-1), &
			      dem_c(1:Npoints,Nlevels:1:-1), &
			      dem_s(1:Npoints,Nlevels:1:-1), cospstateIN, cospIN)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 6) On appelle le simulateur COSPv2 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      print*,'call simulateur'

cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT, 1, Npoints, debug=.false.)

endif ! cosp_init_flag = 1


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 7a) Ecriture des sorties 1: on cree d'abord les fichiers NCDF pour ecrire les sorties
! en appelant lmdz_cosp_output_open (lors du premier appel de cette interface pour les
! 2 options d'ecriture), ou sont definis les axes et les caracteristiques
! des fichiers de sortie avec les diagnostics.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (debut_cosp) then
      !$OMP MASTER

        print *, ' Open outpts files and define axis'
        call lmdz_cosp_output_open(Nlevlmdz, Ncolumns, presnivs, dtime, freq_cosp,       &
                                   ok_mensuelCOSP, ok_journeCOSP, ok_hfCOSP, ok_all_xml, &
                                   ecrit_mth, ecrit_day, ecrit_hf, use_vgrid_in,         &
			           niv_sorties, vgrid_z_in, zlev(1,:))

      !$OMP END MASTER
      !$OMP BARRIER
       debut_cosp=.false.
      endif ! debut_cosp

if (cosp_init_flag .eq. 1) then

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 7b) Ecriture des sorties 2: le remplissage des fichiers de sortie se fait a chaque
! appel de cette interface avec une difference entre les 2 options d'ecriture:
!
! AVEC xios, on commence a remplir les fichiers de sortie a partir du DEUXIEME
! appel de cette interface (lorsque cosp_init_flag = 1).
!
! SANS xios, on commence a remplir les fichiers de sortie a partir du PREMIER
! appel de cette interface (lorsque cosp_init_flag = 1).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       print *, 'Calling write output'
        call lmdz_cosp_output_write(Nlevlmdz, Npoints, Ncolumns, itap, dtime, freq_COSP, &
			            missing_val, cfg, niv_sorties, cospOUT)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 8) On libere la memoire allouee lors de cet appel a l'interface
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call destroy_cospIN(cospIN)
    call destroy_cospstateIN(cospstateIN)
    call destroy_cosp_outputs(cospOUT) 

endif ! cosp_init_flag = 1
 
end subroutine lmdz_cosp_interface
