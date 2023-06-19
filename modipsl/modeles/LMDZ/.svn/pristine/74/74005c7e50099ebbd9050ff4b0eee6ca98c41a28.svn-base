MODULE surf_sisvat_mod
  USE dimphy
  IMPLICIT NONE

  INTEGER, PARAMETER    :: nsnowmx=35
  INTEGER, PARAMETER    :: nsismx=46         ! = nsnowmx + nsoilmx 

CONTAINS

                                       

  SUBROUTINE surf_sisvat(knon,rlon,rlat, ikl2i, itime, dtime, debut, lafin, &
             rmu0, swdown, lwdown, pexner, ps, p1lay, &
             precip_rain, precip_snow, precip_snow_adv, snow_adv, &
             bl_height, wind_velo, temp_air, dens_air, spechum, tsurf, &
             rugos, snow_cont_air, alb_soil, slope, cloudf, &
             radsol, qsol, tsoil, snow, snowhgt, qsnow, to_ice, sissnow, agesno, & 
             AcoefH, AcoefQ, BcoefH, BcoefQ, cdragh, &
             runoff_lic, evap, fluxsens, fluxlat, dflux_s, dflux_l, &      
             tsurf_new, alb1, alb2, alb3, &
             emis_new, z0_new, qsurf)       
                                                                             
! +------------------------------------------------------------------------+   
! |                                                                        |   
! |   SubRoutine surf_sisvat: Interface between LMDZ and landice scheme    |
! |     of the SISVAT (Soil/Ice Snow Vegetation Atmosphere Transfer Scheme)|   
! |                                                                        |
! |   Author: Heinz Juergen Punge, LSCE                June 2009           |
! |     based on the MAR-SISVAT interface by Hubert Gallee                 |
! |                                                                        |   
! +------------------------------------------------------------------------+   
! |   
! |   In the current setup, SISVAT is used only to model the land ice      |
! |   part of the surface; hence it is called with the compressed variables| 
! |   from pbl_surface, and only by the surf_landice routine.              |
! |                                                                        |   
! |   In this interface it is assumed that the partitioning of the soil,   |
! |   and hence the number of grid points is constant during a simulation, |
! |   hence eg. snow properties remain stored in the global SISVAT         |
! |   variables between the calls and don't need to be handed over as      |
! |   arguments. When the partitioning is supposed to change, make sure to |
! |   update the variables.                                                | 
! |                                                                        |  
! |   INPUT                                                                |  
! |   ^^^^^     VegMod: SISVAT    is set up when .T.                       |  
! |             SnoMod: Snow Pack is set up when .T.                       |  
! |             reaLBC: Update Bound.Condit.when .T.                       |   
! |                                                                        |  
! |   INPUT    (via MODULES VARxSV, VARySV, VARtSV)                        |  
! |   ^^^^^     xxxxSV: SISVAT/LMDZ interfacing variables                  |   
! |                                                                        |  
! |   Preprocessing  Option: SISVAT PHYSICS                                |   
! |   ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^                                |   
! | #                       #HY                                            |   
! | #                       #SN: Snow         Model                        |   
! | #                       #BS: Blowing Snow Parameterization             |   
! |                                                                        |   
! | #                       #DS: diffuse radiation differing from direct   |   
! |                             (variable RADsod must still be included)   |
! | #                       #CP: SBL,                       Col de Porte   |   
! | #                       #cp  Solar Radiation,           Col de Porte   |   
! | #                       #AG: Snow Ageing,               Col de Porte   |   
! |                                                                        |   
! |   Preprocessing  Option: SISVAT IO                                     |   
! |   ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |   
! |   FILE                 |      CONTENT                                  |   
! |   ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |   
! | # ANI.yyyymmdd.LAB.nc  | #NC: OUTPUT on NetCDF File (Stand Alone EXP.) |   
! | # SISVAT_iii_jjj_n     | #WV: OUTPUT on ASCII  File (SISVAT Variables) |   
! | # SISVATtroubles       | #VF: OUTPUT on ASCII  File (SISVAT Troubles)  |   
! |                        |                                               |   
! | #                      | #ES: OUTPUT/Verification: Energy Conservation |   
! | #                      | #E2: OUTPUT/Verification: Energy Consrv.2e pt.|   
! |                        |                           (no premature stop) |   
! | #                      | #MW: OUTPUT/Verification: H2O    Conservation |   
! | #                      | #MS: OUTPUT/Verification: * Mass Conservation |   
! | #                      | #MI: OUTPUT/Verification: SeaIce Conservation |   
! |                        |                                               |   
! | # SISVAT__zSn.OUT      | #as: OUTPUT/Verification: Snow Layers Agrega. |   
! | # SISVAT__SnO.OUT      | #aw: OUTPUT/Verification: Albedo Parameteriz. |   
! | # SISVATe_qSn.OUT      | #em: OUTPUT/Verification: Energy/Water Budget |   
! | # SISVATu_qSn.OUT      | #su: OUTPUT/Verification: Slush  Parameteriz. |   
! | # SISVATw_qSo.OUT      | #mw: OUTPUT/Verif+Detail: H2O    Conservation |   
! |                        |                                               |   
! | # SISVAT__GSn.OUT      | #VP: OUTPUT/Verification: Snow   Properties   |   
! | # SISVAT__wEq.OUT      | #EQ: OUTPUT/Verification: Snow/Ice Water Eqv. |   
! |                        |                                               |   
! | #                      | #VR: VERIFICATION OUTPUT                      |   
! | #                      | #WR: Additional   OUTPUT                      |   
! |                                                                        |   
! +------------------------------------------------------------------------+ 
                                                       
    USE VAR_SV 
    USE VARdSV 
    USE VARxSV 
    USE VARySV
    USE VARtSV
 
    USE VARdCP
    USE VARphy, ra_earth=>ra
    USE YOMCST_SISVAT
                                        
    IMPLICIT NONE    
                                                                
! +--INTERFACE Variables                                                     
! +  =================== 
                                            
    include  "dimsoil.h"                                     
                                                                    

! +--Global Variables                                                          
! +  ================  
! Input Variables for SISVAT
    INTEGER,               INTENT(IN)      :: knon
    INTEGER,               INTENT(IN)      :: itime   
    REAL,                  INTENT(IN)      :: dtime 
    LOGICAL,               INTENT(IN)      :: debut     ! true if first step
    LOGICAL,               INTENT(IN)      :: lafin     ! true if last step

    INTEGER, DIMENSION(klon), INTENT(IN)   :: ikl2i     ! Index Decompression
    REAL, DIMENSION(klon), INTENT(IN)      :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)      :: rmu0      ! cos sol. zenith angle
    REAL, DIMENSION(klon), INTENT(IN)      :: swdown    !
    REAL, DIMENSION(klon), INTENT(IN)      :: lwdown    ! 
    REAL, DIMENSION(klon), INTENT(IN)      :: pexner    ! Exner potential
    REAL, DIMENSION(klon), INTENT(IN)      :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)      :: precip_snow_adv, snow_adv 
                                                        !Snow Drift
    REAL, DIMENSION(klon), INTENT(IN)      :: bl_height, wind_velo
    REAL, DIMENSION(klon), INTENT(IN)      :: temp_air, spechum, ps,p1lay
    REAL, DIMENSION(klon), INTENT(IN)      :: dens_air, tsurf            
    REAL, DIMENSION(klon), INTENT(IN)      :: rugos,snow_cont_air
    REAL, DIMENSION(klon), INTENT(IN)      :: alb_soil, slope 
    REAL, DIMENSION(klon), INTENT(IN)      :: cloudf   
    REAL, DIMENSION(klon), INTENT(IN)      :: AcoefH, AcoefQ
    REAL, DIMENSION(klon), INTENT(IN)      :: BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)      :: cdragh

! Variables exchanged between LMDZ and SISVAT
    REAL, DIMENSION(klon,nsoilmx), INTENT(OUT) :: tsoil ! Soil Temperature 
    REAL, DIMENSION(klon), INTENT(OUT)     :: qsol      ! Soil Water Content   
    REAL, DIMENSION(klon), INTENT(INOUT)   :: snow      ! Tot snow mass [kg/m2]
    REAL, DIMENSION(klon), INTENT(IN)      :: radsol    ! Surface absorbed rad.

! Output Variables from SISVAT
    REAL, DIMENSION(klon), INTENT(OUT)     :: alb1      ! Albedo SW 
    REAL, DIMENSION(klon), INTENT(OUT)     :: alb2,alb3      ! Albedo LW
    REAL, DIMENSION(klon), INTENT(OUT)     :: emis_new  ! Surface Emissivity
    REAL, DIMENSION(klon), INTENT(OUT)     :: z0_new    ! Momentum Roughn Lgt
    REAL, DIMENSION(klon), INTENT(OUT)     :: runoff_lic ! Runoff            
    REAL, DIMENSION(klon), INTENT(OUT)     :: dflux_s   ! d/dT sens. ht flux  
    REAL, DIMENSION(klon), INTENT(OUT)     :: dflux_l   ! d/dT latent ht flux 
    REAL, DIMENSION(klon), INTENT(OUT)     :: fluxsens  ! Sensible ht flux    
    REAL, DIMENSION(klon), INTENT(OUT)     :: fluxlat   ! Latent heat flux
    REAL, DIMENSION(klon), INTENT(OUT)     :: evap      ! Evaporation
    REAL, DIMENSION(klon), INTENT(OUT)     :: agesno    ! Snow age (top layer)
    REAL, DIMENSION(klon), INTENT(OUT)     :: tsurf_new ! Surface Temperature 
    REAL, DIMENSION(klon), INTENT(OUT)     :: qsurf     ! Surface Humidity 
    REAL, DIMENSION(klon), INTENT(OUT)     :: qsnow     ! Total H2O snow[kg/m2]
    REAL, DIMENSION(klon), INTENT(OUT)     :: snowhgt   ! Snow height (m)
    REAL, DIMENSION(klon), INTENT(OUT)     :: to_ice    ! Snow passed to ice 
    REAL, DIMENSION(klon), INTENT(OUT)     :: sissnow   ! Snow in model (kg/m2)
                                                                           
! +--OUTPUT for NetCDF File                                         
! +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                             
!   REAL, DIMENSION(klonv)                 :: SOsoKL    ! Abs Solar Radiation 
!   REAL, DIMENSION(klonv)                 :: IRsoKL    ! Abs IR    Radiation 
!   REAL, DIMENSION(klonv)                 :: HSsoKL    ! Abs Sensible Ht Flux 
!   REAL, DIMENSION(klonv)                 :: HLsoKL    ! Abs Latent Heat Flux 
!   REAL, DIMENSION(klonv)                 :: HLs_KL    ! Evaporation          
!   REAL, DIMENSION(klonv)                 :: HLv_KL    ! Transpiration        
!                              
!   REAL, DIMENSION(klonv)                 :: SOsoNC    ! Abs Solar Radiation  
!   REAL, DIMENSION(klonv)                 :: IRsoNC    ! Abs IR    Radiation  
!   REAL, DIMENSION(klonv)                 :: HSsoNC    ! Abs Sensible Ht Flux
!   REAL, DIMENSION(klonv)                 :: HLsoNC    ! Abs Latent Heat Flux
!   REAL, DIMENSION(klonv)                 :: HLs_NC    ! Evaporation      
!   REAL, DIMENSION(klonv)                 :: HLv_NC    ! Transpiration     
!                                      
!   REAL, DIMENSION(klonv,nsoilmx)           :: eta_NC    ! nsoilmx=nsol+1
!   REAL, DIMENSION(klonv,nsoilmx)           :: tsolNC    !   
!   REAL, DIMENSION(klonv)                   :: snowNC    !
!   INTEGER, DIMENSION(klonv)                :: isnoNC    !
!   INTEGER, DIMENSION(klonv)                :: ispiNC    !     
!   INTEGER, DIMENSION(klonv)                :: iiceNC    !           
!   REAL, DIMENSION(klonv)                   :: swaNC     !
!   REAL, DIMENSION(klonv)                   :: swsNC     !   
!   INTEGER, DIMENSION(klonv,nsno)           :: istoNC    ! 
!   REAL, DIMENSION(klonv,nsno)              :: dzsnNC    ! 
!   REAL, DIMENSION(klonv,nsno)              :: rhosnNC   !      
!   REAL, DIMENSION(klonv,nsno)              :: etasnNC   !  
!   REAL, DIMENSION(klonv,nsno)              :: tsnNC     ! 
!   REAL, DIMENSION(klonv,nsno)              :: g1snNC    !      
!   REAL, DIMENSION(klonv,nsno)              :: g2snNC    !                    
!   REAL, DIMENSION(klonv,nsno)              :: agsnNC    ! 
!   REAL, DIMENSION(klonv)                   :: meltNC    !      
!   REAL, DIMENSION(klonv)                   :: refrNC    !                
!   REAL, DIMENSION(klonv)                   :: alb1NC    !      
!   REAL, DIMENSION(klonv)                   :: alb2NC    !         
!   REAL, DIMENSION(klonv)                   :: rnofNC    ! 
                                                                          
! + Optional Variables:                                                 
! +--V,  dT(a-s)    Time Moving Averages                                       
! +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~           

! #AW real             V__mem(klon,ntaver)   ! ntaver defined in LMDZ_SL.inc  
! #AW real             VVmmem(klon)          !                                 
! #AW common/SVeSBLmem/V__mem,VVmmem         !                                 
! #AH real             T__mem(klon,ntaver)   !                                 
! #AH real             dTmmem(klon)          !                                 
! #AH common/STeSBLmem/T__mem,dTmmem         !                      
! +--u*, u*T*, u*s* Time Moving Averages                                    
! +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                       
! #AM real             u__mem(klon,ntaver)   ! ntaver defined in LMDZ_SL.inc  
! #AT real             uT_mem(klon,ntaver)   !                                 
! #AS real             us_mem(klon,ntaver)   !                                 
! #AM common/S_eSBLmem/u__mem                !                                 
! #AT.                ,uT_mem                !                                 
! #AS.                ,us_mem                !
!! #AH    REAL, DIMENSION(klonv,ntaver)            :: T__mem    !
!! #AH    REAL, DIMENSION(klonv)                   :: dTmmem    ! 
!! #AW    REAL, DIMENSION(klonv,ntaver)            :: V__mem    !      
!! #AM    REAL, DIMENSION(klonv)                   :: u__mem    !             
!! #AT    REAL, DIMENSION(klonv)                   :: uT_mem    ! 
!! #AS    REAL, DIMENSION(klonv)                   :: us_mem    !     
!! #AW    REAL, DIMENSION(klonv)                   :: VVmmem    !    

                                                                    
! +--Internal  Variables                                                       
! +  ===================                                          

    CHARACTER(len=20)               :: fichnom, fn_outfor ! Name for output file
    INTEGER                         :: i, ig, ikl, isl, isn, nt
    INTEGER                         :: gp_outfor, un_outfor

    REAL, PARAMETER                 :: f1=0.5
    REAL, PARAMETER                 :: sn_upp=5000.,sn_low=500.
    REAL, PARAMETER                 :: sn_add=400.,sn_div=2.
                                             ! snow mass upper,lower limit, 
                                             ! added mass/division lowest layer
    REAL, PARAMETER                 :: c1_zuo=12.960e+4, c2_zuo=2.160e+6
    REAL, PARAMETER                 :: c3_zuo=1.400e+2,  czemin=1.e-3  
                                             ! Parameters for drainage
! c1_zuo/ 2.796e+4/,c2_zuo/2.160e+6/,c3_zuo/1.400e+2/ !     Tuning
! +...        Run Off Parameters                                              
! +           86400*1.5 day     ...*25 days (Modif. ETH Camp: 86400*0.3day)    
! +           (Zuo and Oerlemans 1996, J.Glacio. 42, 305--317)             

    REAL, DIMENSION(klon)           :: eps0SL          ! surface Emissivity 
    REAL                            :: zsigma, Ua_min, Us_min
    REAL                            :: lambda          ! Par. soil discret.
    REAL, DIMENSION(nsoilmx), SAVE  :: dz1,dz2         ! Soil layer thicknesses
!$OMP THREADPRIVATE(dz1)
    LOGICAL, SAVE                   :: firstcall=.TRUE.,SnoMod, ok_outfor=.FALSE.!$OMP THREADPRIVATE(firstcall)                



! + Optional:  
!c #BW INTEGER                      ::  noUNIT                  
!c #BW REAL                         ::  BlowST,SnowSB   
                                          
                                        
! +--Internal Variables
! +  ==================

    INTEGER                         ::  ivg,iso

!========================================================================
 
      SnoMod=.true.
      zsigma=1000.
      dt__SV=dtime
      IF (ok_outfor) THEN
        un_outfor=51                 ! unit number for point output
        gp_outfor=79! 633 !79        ! grid point number for point output
        fn_outfor='outfor_SV.dat' 
      END IF

!     write(*,*)'Start of simulation? ',debut        !hj
      IF (debut) THEN 
        firstcall=.TRUE. 
        INI_SV=.false.
      ELSE
        firstcall=.false.
        INI_SV=.true.
      END IF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + INITIALISATION: BEGIN +++ 
! + -------------------------
! +
! + Compute soil discretization (as for LMDZ)
! + -----------------------------------------  
      IF (firstcall) THEN

! +--Array size
        klonv=klon
        knonv=knon
        write(*,*)'klon',klon,'klonv',klonv,'knon',knon 


        CALL INIT_VARtSV
        CALL INIT_VARxSV
        CALL INIT_VARySV

        eps0SL(:)=0.
! +--Soil layer thickness                                                   
! +  -----------------------  
!        write(*,'(/a)') 'Start SISVAT init: soil discretization ', nsoilmx
        CALL get_soil_levels(dz1,dz2,lambda)
        lambSV=lambda
        dz1_SV(1:knon,1:) = 0.     
        dz2_SV(1:knon,1:) = 0.
                      
        DO isl =   -nsol,0   
          dz_dSV(isl) = 0.5e-3*dz2(1-isl)           ! Soil layer thickness 
          DO ikl=1,knon
            dz1_SV(ikl,isl) = dz1(1-isl)    !1.e-3* 
            dz2_SV(ikl,isl) = dz2(1-isl)    !1.e-3*
          END DO
!          IF (knon > 0) THEN
!            write(*,*)'level:',dz_dSV(isl),dz1_SV(1,isl),dz2_SV(1,isl)
!          END IF
        END DO

        DO ikl=1,knon                                
          eps0SL(ikl )= 1.
          alb0SV(ikl) = alb_soil(ikl)              ! Soil Albedo        

! + Soil Upward IR Flux, Water Fluxes, roughness length  
          IRs_SV(ikl) =                               &
              -eps0SL(ikl)* rsigma*temp_air(ikl)      &  ! Upward IR Flux
              *temp_air(ikl)*temp_air(ikl)*temp_air(ikl)         
          TvegSV(ikl) = temp_air(ikl) 
 ! + Soil
        DO isl =   -nsol,0    
          TsisSV(ikl,isl) = temp_air(ikl)  !tsoil(ikl,1-isl)  Soil Temperature
          eta_SV(ikl,isl) = 0.0001         !etasoil(ikl,1-isl)Soil Water[m3/m3]
          ro__SV(ikl,isl) = 50.      !rosoil(ikl,1-isl) soil water volumic mass
        END DO            
        END DO                                            

! +--Surface Fall Line Slope                                                   
! +  -----------------------  
        IF (SnoMod)  THEN                
          DO ikl=1,knon  
            slopSV(ikl) = slope(ikl)
            SWf_SV(ikl) =             &   ! Normalized Decay of the  
              exp(-dt__SV             &   ! Surficial Water Content  
              /(c1_zuo                &   !(Zuo and Oerlemans 1996,  
            +c2_zuo*exp(-c3_zuo*abs(slopSV(ikl)))))  ! J.Glacio. 42, 305--317)
          END DO                                     
        END IF           
                             
! +--SBL  Characteristics history, not stored at restart.                     
        DO ikl=1,knon  
         DO nt=1,ntaver                                       ! #AW
           V__mem(ikl,nt)=wind_velo(ikl)                      ! #AH
           T__mem(ikl,nt)=temp_air(ikl)-tsurf(ikl) 
         END DO                                         
        END DO                                     


! + SISVAT_ini (as for use with MAR, but not computing soil layers)
! + -------------------------------------------------------------
!        write(*,'(/a)') 'Start SISVAT initialization: SISVAT_ini' 
        CALL SISVAT_ini(knon)

! open output file 
        IF (ok_outfor) THEN
          open(unit=un_outfor,status='new',file=fn_outfor)         
          rewind un_outfor      
          ikl=gp_outfor     ! index sur la grille land ice
          ig=611            ! index sur la grille globale
          write(un_outfor,501) fn_outfor, ikl, rlon(ig),rlat(ig)  
501     format(/,a18,/,'Grid point ',i4,' Long',f9.4,' Lat ',f9.4       &
     &            ,/,'++++++++++++++++++++++++++++++++++++++++++++++',  &
     & '++++++++++++++++++++++++++++++++++++++++++++++',                &
     & '++++++++++++++++++++++++++++++++++++++++++++++',                &
     &             /,' SWdown + IRdown + Wind   + Temp.  + Humid.   ',  &
     & '+ Press  +Precip_l+Precip_s+ Tsrf   + Clouds +'                 &
     & '+ Zenith + BLhgt  + Densair+ Exner +'                           &
     &            ,/,' sol_SV + IRd_SV + VV__SV + TaT_SV + QaT_SV   ',  &
     & '+ ps__SV + drr_SV + dsn_SV + Tsf_SV + cld_SV +'                 &
     & '+ coszSV + za__SV + rhT_SV + ExnrSV+'                           &
     &            ,/,' W/m2   + W/m2   + m/s    + K      + kg/kg    ',  &
     & '+ Pa     + kg/m2/s+ kg/m2/s+ K      + /1     +'        &
     & '+ -      + m      + kg/m3  +       +'                           &
     &            ,/,'++++++++++++++++++++++++++++++++++++++++++++++',  &
     & '++++++++++++++++++++++++++++++++++++++++++++++',                &
     & '++++++++++++++++++++++++++++++++++++++++++++++')
        END IF

! +--Read restart file
! +  =================================================   
        
! Martin
       PRINT*, 'On debranche sisvatetat0'
! Martin

       ! CALL sisvatetat0("startsis.nc",ikl2i)

 
      END IF  ! firstcall                        
! +                                
! +  +++  INITIALISATION:  END  +++                               
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + READ FORCINGS  
! + ------------------------ 
! + Update Forcings for SISVAT given by the LMDZ model.
! +
      DO ikl=1,knon

! +--Atmospheric Forcing                                    (INPUT)            
! +  ^^^^^^^^^^^^^^^^^^^                                     ^^^^^ 
        zSBLSV      = 1000.                         ! [m]                
        za__SV(ikl) = bl_height(ikl)                ! Height boundary layer [m]
        Ua_min      = epsi                          !
! #VM   Ua_min      = 0.2 * sqrt(za__SV(ikl)   )    !                    
        Ua_min      = 0.2 * sqrt(za__SV(ikl)   )    !                    
        VV__SV(ikl) = max(Ua_min, wind_velo(ikl))   ! Wind velocity       [m/s]
        Us_min      = 0.01
        us__SV(ikl) = max(Us_min, us__SV(ikl) ) 
        TaT_SV(ikl) = temp_air(ikl)                 ! BL top Temperature    [K]
        ExnrSV(ikl) = pexner(ikl)                   ! Exner potential          
        rhT_SV(ikl) = dens_air(ikl)                 ! Air density           
        QaT_SV(ikl) = spechum(ikl)                  ! Specific humidity 
! #VX   dQa_SV(ikl) = 0. !hj   dtDiff/zsigma(mz)    ! Water Vapor Flux
        ps__SV(ikl) = ps(ikl)                       ! surface pressure     [Pa]
        p1l_SV(ikl) = p1lay(ikl)                  ! lowest atm. layer press[Pa]

! +--Energy Fluxes                                          (INPUT)           
! +  ^^^^^^^^^^^^^                                           ^^^^^             
        coszSV(ikl) = max(czemin,rmu0(ikl))         ! cos(zenith.Dist.)  
        sol_SV(ikl) = swdown(ikl)                   ! downward Solar  
        IRd_SV(ikl) = lwdown(ikl)                   ! downward IR    
        rsolSV(ikl) = radsol(ikl)                   ! surface absorbed rad.   
!hj 110511
!hj 121011        IRs_SV(ikl) = min(IRs_SV(ikl),-epsi)        ! check upward IR    
!hj
! +--Water  Fluxes                                          (INPUT)           
! +  ^^^^^^^^^^^^^                                           ^^^^^             
        drr_SV(ikl) = precip_rain(ikl)              ! Rain fall rate  [kg/m2/s]
        dsn_SV(ikl) = precip_snow(ikl)              ! Snow fall rate  [kg/m2/s]
!c #BS  dbsnow      = -SLussl(i,j,n)                ! Erosion    
!c #BS.               *dtPhys     *rhT_SV(ikl) /ro_Wat                    
!c #BS  dsnbSV(ikl) = snow_adv(ikl)  ! min(max(zero,dbsnow)              
!c #BS.                    /    max(epsi,d_snow),unun)                    
!c #BS  dbs_SV(ikl) = snow_cont_air(ikl)
!c #BS                  blowSN(i,j,n)               !          [kg/m2]  
                                                                              
! +--Soil/Canopy                                            (INPUT)            
! +  ^^^^^^^^^^^                                             ^^^^^            
        alb0SV(ikl) = alb_soil(ikl)                 ! Soil background Albedo 
        AcoHSV(ikl) = AcoefH(ikl)  
        BcoHSV(ikl) = BcoefH(ikl)                     
        AcoQSV(ikl) = AcoefQ(ikl)  
        BcoQSV(ikl) = BcoefQ(ikl)
        cdH_SV(ikl) = cdragh(ikl)                         
        Tsf_SV(ikl) = tsurf(ikl)        !hj 12 03 2010

! +--Energy Fluxes                                          (INPUT/OUTPUT)    
! +  ^^^^^^^^^^^^^                                           ^^^^^^^^^^^^    
        IF (.not.firstcall) THEN  
! active hj 110411
        cld_SV(ikl) = cloudf(ikl)                    ! Cloudiness         
        END IF
 
! +--Time Averages of wind and surface-atmosphere temperature difference 
! +--for turbulence calculations                                       ! #AA
! Update stored arrays
      DO nt=1,ntaver-1                                                 ! #AW   
        V__mem(ikl,nt    ) = V__mem(ikl,nt+1)                          ! #AH   
        T__mem(ikl,nt    ) = T__mem(ikl,nt+1)                          ! #AA  
      ENDDO                                                            ! #AW
        V__mem(ikl,ntaver)=wind_velo(ikl)                              ! #AH
        T__mem(ikl,ntaver)=temp_air(ikl)-tsurf(ikl)                    ! #AW
! Calculate averages
        VVmmem(ikl)        = 0.0                                       ! #AH
        dTmmem(ikl)        = 0.0                                       ! #AA
      DO nt=1,ntaver                                                   ! #AW
        VVmmem(ikl)        = VVmmem(ikl)+V__mem(ikl,nt)                ! #AH
        dTmmem(ikl)        = dTmmem(ikl)+T__mem(ikl,nt)                ! #AA 
      ENDDO                                                            ! #AW 
        VVmmem(ikl)        = VVmmem(ikl)/ntaver                        ! #AH 
        dTmmem(ikl)        = dTmmem(ikl)/ntaver                   
      END DO

!                            
! +  +++  READ FORCINGS:  END  +++    
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +  OUTPUT FORCING
!
      IF (ok_outfor) THEN
        ikl=gp_outfor
        write(un_outfor,5000) sol_SV(ikl), IRd_SV(ikl), VV__SV(ikl),    &    
     &                        TaT_SV(ikl), QaT_SV(ikl), ps__SV(ikl),    &
     &                        drr_SV(ikl), dsn_SV(ikl), Tsf_SV(ikl),    &
     &                        cld_SV(ikl), coszSV(ikl), za__SV(ikl),    &
     &                        rhT_SV(ikl), ExnrSV(ikl) 
      
 5000     format(f8.3,' ',f8.3,' ',f8.4,' ',f8.4,' ',f10.8,' ',f8.1,' ', &
     &           f8.6,' ',f8.6,' ',f8.4,' ',f8.6,' ',f8.6,' ',f8.3,' ',  &
     &           f8.4,' ',f8.6)

      ENDIF
!                            
! +  OUTPUT FORCINGS:  END  +++    
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +--SISVAT EXECUTION                                                          
! +  ----------------                                                          
! +             
!      write(*,*) '  Start SISVAT execution!'

      call  SISVAT(SnoMod,.false.,1) 
!                          BloMod,jjtime)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + RETURN RESULTS  
! + -------------- 
! + Return (compressed) SISVAT variables to LMDZ             
! + 
      DO  ikl=1,knon                  ! use only 1:knon (actual ice sheet..) 
        runoff_lic(ikl)    = RnofSV(ikl)*dtime   ! RunOFF: intensity* time step
        dflux_s(ikl)       = dSdTSV(ikl)         ! Sens.H.Flux T-Der. 
        dflux_l(ikl)       = dLdTSV(ikl)         ! Latn.H.Flux T-Der. 
        fluxsens(ikl)      = HSs_sv(ikl)         ! HS                 
        fluxlat(ikl)       = HLs_sv(ikl)         ! HL                 
        evap(ikl)          = HLs_sv(ikl)/RLVTT   ! Evaporation RLVTT=LhvH2O   
        z0_new(ikl)        = Z0h_SV(ikl)         ! Moment.Roughn.L.


        snow(ikl)          = 0.
        snowhgt(ikl)       = 0.
        qsnow(ikl)         = 0.
        qsol(ikl)          = 0. 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +   Check snow thickness, add if too thin, substract if too thick     

        sissnow(ikl)       = 0.  !()
      DO  isn = 1,isnoSV(ikl)                                               
        sissnow(ikl)       = sissnow(ikl)+dzsnSV(ikl,isn)* ro__SV(ikl,isn)   
      END DO

      IF (sissnow(ikl) .LE. sn_low) THEN  !add snow 
      IF (isnoSV(ikl).GE.1) THEN
        dzsnSV(ikl,1)      = dzsnSV(ikl,1) + sn_add/max(ro__SV(ikl,1),epsi)  
        toicSV(ikl)        = toicSV(ikl)   - sn_add
      ELSE
        write(*,*) 'Attention, bare ice... point ',ikl
        isnoSV(ikl)        = 1
        istoSV(ikl,1)      = 100
!ym        istoSV(ikl)     = 100
        ro__SV(ikl,1)      = 400.      
        dzsnSV(ikl,1)      = sn_add/max(ro__SV(ikl,1),epsi)  ! 1.
        eta_SV(ikl,1)      = epsi
        TsisSV(ikl,1)      = min(TsisSV(ikl,0),TfSnow-0.2)   
        G1snSV(ikl,1)      = 99.
        G2snSV(ikl,1)      = 0.3 
        agsnSV(ikl,1)      = 10.    
        toicSV(ikl)        = toicSV(ikl)   - sn_add
      END IF
      END IF

      IF (sissnow(ikl) .ge. sn_upp) THEN  !thinnen snow layer below
        dzsnSV(ikl,1)      = dzsnSV(ikl,1)/sn_div
        toicSV(ikl) = toicSV(ikl)+dzsnSV(ikl,1)*ro__SV(ikl,1)/sn_div
      END IF

        sissnow(ikl)       = 0.  !()
      DO  isn = 1,isnoSV(ikl)                                               
        sissnow(ikl) = sissnow(ikl)+dzsnSV(ikl,isn)* ro__SV(ikl,isn)           
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        snowhgt(ikl) = snowhgt(ikl)+dzsnSV(ikl,isn)              
        qsnow(ikl)   = qsnow(ikl)+1e03*eta_SV(ikl,isn)*dzsnSV(ikl,isn)   
      END DO
        snow(ikl)    = sissnow(ikl)+toicSV(ikl)
        to_ice(ikl)  = toicSV(ikl)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      DO  isl =   -nsol,0    
        tsoil(ikl,1-isl)   = TsisSV(ikl,isl)       ! Soil Temperature  
        qsol(ikl)          = qsol(ikl)                      &    
                       +eta_SV(ikl,isl) * dz_dSV(isl)  
      END DO                                               
        agesno(ikl)        = agsnSV(ikl,isnoSV(ikl))        !          [day] 

        alb1(ikl)          = alb1sv(ikl)             ! Albedo VIS  
        alb2(ikl)          = ((So1dSV-f1)*alb1sv(ikl)                   &
     &                       +So2dSV*alb2sv(ikl)+So3dSV*alb3sv(ikl))/f1    
                                                     ! Albedo NIR 
        alb3(ikl)          = alb3sv(ikl)             ! Albedo FIR  
        tsurf_new(ikl)     = TsfnSV(ikl)
!hj220711          Tsrfsv(ikl)              ! Surf.Temperature   
!                TsisSV(ikl,0) *(1-min(1,isnoSV(ikl))) &                   
!                +TsisSV(ikl,max(1,isnoSV(ikl))) * min(1,isnoSV(ikl)) 
!        tsurf_new(ikl)   = max(Ts_Min,tsurf_new(ikl) )
!        tsurf_new(ikl)   = min(Ts_Max,tsurf_new(ikl) )
        qsurf(ikl)         = QaT_SV(ikl) 
        emis_new(ikl)      = eps0SL(ikl)  

!!!hjp 230611 sorties
!!        qsnow(ikl)=TsisSV(ikl,isnoSV(ikl))
!!        sissnow(ikl)=TsisSV(ikl,isnoSV(ikl)-1)
!!        snowhgt(ikl)=TsisSV(ikl,isnoSV(ikl)-2)
!!        qsol(ikl)=dzsnSV(ikl,isnoSV(ikl)-1)
!!        agesno(ikl)=dzsnSV(ikl,isnoSV(ikl))


      END DO 


! +  -----------------------------                              
! +  END --- RETURN RESULTS    
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + Uncompress SISVAT output variables & store
! +  -----------------------------            
! +
!      DO  ikl = 1,knon                                                     
!        i   = ikl2i(ikl)                             ! Compression index  
!      IF (i.GE.0) THEN
!        SOsoNC(i)       = SOsoKL(ikl)                ! Absorb.Sol.Rad.    
!        IRsoNC(i)       = IRsoKL(ikl)                ! Absorb.IR  Rad.    
!        HSsoNC(i)       = HSsoKL(ikl)                ! HS                 
!        HLsoNC(i)       = HLsoKL(ikl)                ! HL                 
!        HLs_NC(i)       = HLs_KL(ikl)                ! Evaporation          
!        HLv_NC(i)       = HLv_KL(ikl)                ! Transpiration      
!        
!        DO isl =   -nsol,0      
!          eta_NC(i,1-isl) = eta_SV(ikl,isl)          ! Soil Humidity      
!          tsolNC(i,1-isl) = TsisSV(ikl,isl)          ! Soil Temperature   
!        END DO                                                  
!        snowNC(i)       = snow(ikl)                  ! Snow mass    
!        isnoNC(i)       = isnoSV(ikl)                ! Nb Snow/Ice Lay.   
!        ispiNC(i)       = ispiSV(ikl)                ! Nb Supr.Ice Lay.   
!        iiceNC(i)       = iiceSV(ikl)                ! Nb      Ice Lay.   
!        swaNC(i)        = rusnSV(ikl)                ! Surficial Water    
!        swsNC(i)        = SWS_SV(ikl)                ! Surficial Wat.St.    
!        DO  isn = 1,nsno                                                   
!          istoNC(i,isn)   = istoSV(ikl,isn)          !            [-]     
!          dzsnNC(i,isn)   = dzsnSV(ikl,isn)          !            [m]     
!          rhosnNC(i,isn)  = ro__SV(ikl,isn)          !        [kg/m3]     
!          etasnNC(i,isn)  = eta_SV(ikl,isn)          !        [m3/m3]     
!          tsnNC(i,isn)    = TsisSV(ikl,isn)          !            [K]     
!          g1snNC(i,isn)   = G1snSV(ikl,isn)          ! [-]        [-]     
!          g2snNC(i,isn)   = G2snSV(ikl,isn)          ! [-] [0.0001 m]     
!          agsnNC(i,isn)   = agsnSV(ikl,isn)          !          [day]   
!        END DO                                                            
!!?c #IB  depsubNC(i)    = wes_SV(ikl)                ! Depo. / Subli.     
!        meltNC(i)       = wem_SV(ikl)                ! Melting            
!        refrNC(i)       = wer_SV(ikl)                ! Refreezing         
!        alb1NC(i)       = alb1sv(ikl)                ! Albedo SW  
!        alb2NC(i)       = (alb1sv(ikl)+3*alb2sv(ikl)+alb3sv(ikl))/5  
                                                      ! Albedo LW
!        rnofNC(i)       = RnofSV(ikl)              ! Run OFF Intensity
!   
!!?        SL_z0(i)       = Z0m_SV(ikl)                 ! Moment.Roughn.L.   
!!?        SL_r0(i,j,n)   = Z0h_SV(ikl)                 ! Heat   Roughn.L.  
!!!        zWE_NC(i)       = zWE_SV(ikl)                ! Current  *Thick.   
!!!        zWEcNC(i)       = zWEcSV(ikl)                ! Non-Erod.*Thick.   
!!!        hSalNC(i)       = hSalSV(ikl)                ! Salt.Layer Height  
!!!        hsenSL(i) = -SLuts(i,j)  *  rhAir    *cp     ! Sensible Heat Flux
!!!        hlatSL(i) = -SLuqs(i,j)  *  rhAir    *Lv_H2O ! Latent   Heat Flux 
!      END IF
!      END DO   

!hj + Put storage to Output file here!                                     
      IF (lafin) THEN

        fichnom = "restartsis.nc"
        CALL sisvatredem("restartsis.nc",ikl2i,rlon,rlat)            
        IF (ok_outfor) THEN
          close(unit=un_outfor)                    
        END IF
      END IF                                                          
! +  -----------------------------                              
! +  END --- RETURN RESULTS    
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END SUBROUTINE surf_sisvat







  SUBROUTINE get_soil_levels(dz1, dz2, lambda)
! ======================================================================
! Routine to compute the vertical discretization of the soil in analogy
! to LMDZ. In LMDZ it is done in soil.F, which is not used in the case 
! of SISVAT, therefore it's needed here.
!
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    USE mod_phys_lmdz_para

    INCLUDE "dimsoil.h"

    REAL, DIMENSION(nsoilmx), INTENT(OUT) :: dz2, dz1
    REAL, INTENT(OUT)                     :: lambda


!-----------------------------------------------------------------------
!   Depthts:
!   --------
    REAL fz,rk,fz1,rk1,rk2
    REAL min_period, dalph_soil
    INTEGER ierr,jk

    fz(rk)=fz1*(dalph_soil**rk-1.)/(dalph_soil-1.)

!    write(*,*)'Start soil level computation' 
!-----------------------------------------------------------------------
! Calculation of some constants
! NB! These constants do not depend on the sub-surfaces
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!   ground levels 
!   grnd=z/l where l is the skin depth of the diurnal cycle:
!-----------------------------------------------------------------------

     min_period=1800. ! en secondes
     dalph_soil=2.    ! rapport entre les epaisseurs de 2 couches succ.
! !$OMP MASTER
!     IF (is_mpi_root) THEN
!        OPEN(99,file='soil.def',status='old',form='formatted',iostat=ierr)
!        IF (ierr == 0) THEN ! Read file only if it exists
!           READ(99,*) min_period
!           READ(99,*) dalph_soil
!           PRINT*,'Discretization for the soil model'
!           PRINT*,'First level e-folding depth',min_period, &
!                '   dalph',dalph_soil
!           CLOSE(99)
!        END IF
!     ENDIF
! !$OMP END MASTER
!     CALL bcast(min_period)
!     CALL bcast(dalph_soil)

!   la premiere couche represente un dixieme de cycle diurne
     fz1=SQRT(min_period/3.14)
     
     DO jk=1,nsoilmx
        rk1=jk
        rk2=jk-1
        dz2(jk)=fz(rk1)-fz(rk2) 
     ENDDO
     DO jk=1,nsoilmx-1
        rk1=jk+.5
        rk2=jk-.5
        dz1(jk)=1./(fz(rk1)-fz(rk2))
     ENDDO
     lambda=fz(.5)*dz1(1)
     PRINT*,'full layers, intermediate layers (seconds)'
     DO jk=1,nsoilmx
        rk=jk
        rk1=jk+.5
        rk2=jk-.5
        PRINT *,'fz=', &
             fz(rk1)*fz(rk2)*3.14,fz(rk)*fz(rk)*3.14
     ENDDO

  END SUBROUTINE get_soil_levels

  SUBROUTINE SISVAT_ini(knon)                                                     
                                                                               
!C +------------------------------------------------------------------------+  
!C | MAR          SISVAT_ini                             Jd 11-10-2007  MAR | 
!C |   SubRoutine SISVAT_ini generates non time dependant SISVAT parameters |  
!C +------------------------------------------------------------------------+  
!C |   PARAMETERS:  klonv: Total Number of columns =                        |  
!C |   ^^^^^^^^^^        = Total Number of continental     grid boxes       |  
!C |                     X       Number of Mosaic Cell per grid box         |  
!C |                                                                        |  
!C |   INPUT:   dt__SV   : Time  Step                                   [s] | 
!C |   ^^^^^    dz_dSV   : Layer Thickness                              [m] |  
!C |                                                                        |  
!C |   OUTPUT:  RF__SV   : Root Fraction in Layer isl                   [-] |  
!C |   ^^^^^^   rocsSV   : Soil Contrib. to (ro c)_s exclud.Water  [J/kg/K] |  
!C |            etamSV   : Soil Minimum Humidity                    [m3/m3] |  
!C |                      (based on a prescribed Soil Relative Humidity)    |  
!C |            s1__SV   : Factor of eta**( b+2) in Hydraul.Diffusiv.       |  
!C |            s2__SV   : Factor of eta**( b+2) in Hydraul.Conduct.        |  
!C |            aKdtSV   : KHyd: Piecewise Linear Profile:  a * dt    [m]   |  
!C |            bKdtSV   : KHyd: Piecewise Linear Profile:  b * dt    [m/s] |  
!C |            dzsnSV(0): Soil first Layer Thickness                   [m] |  
!C |            dzmiSV   : Distance between two contiguous levels       [m] |  
!C |            dz78SV   : 7/8 (Layer Thickness)                        [m] | 
!C |            dz34SV   : 3/4 (Layer Thickness)                        [m] | 
!C |            dz_8SV   : 1/8 (Layer Thickness)                        [m] | 
!C |            dzAvSV   : 1/8  dz_(i-1) + 3/4 dz_(i) + 1/8 dz_(i+1)    [m] | 
!C |            dtz_SV   : dt/dz                                      [s/m] | 
!C |            OcndSV   : Swab Ocean / Soil Ratio                      [-] |
!C |            Implic   : Implicit Parameter  (0.5:  Crank-Nicholson)      |  
!C |            Explic   : Explicit Parameter = 1.0 - Implic                |  
!C |                                                                        | 
!C | # OPTIONS: #ER: Richards Equation is not smoothed                      |
!C | # ^^^^^^^  #kd: De Ridder   Discretization                             |
!C | #          #SH: Hapex-Sahel Values                                     !  
!C |                                                                        | 
!C +------------------------------------------------------------------------+  
!                                                                              
!                                                                             
                                                                              
!C +--Global Variables                                                         
!C +  ================         

      USE VARphy, ra_earth=>ra                                           
      USE VAR_SV                                                     
      USE VARdSV                                                         
      USE VAR0SV                                                           
      USE VARxSV
      USE VARtSV
      USE VARxSV
      USE VARySV
      IMPLICIT NONE                                                           
                                                                              
                                                                              
                                                                              
!C +--Arguments                                                     
!C +  ==================                                                      
       INTEGER,INTENT(IN) ::  knon                                       

!C +--Internal Variables                                                     
!C +  ==================                                                      
                                                                               
      INTEGER ::  ivt   ,ist   ,ivg   ,ikl   ,isl   ,isn   ,ikh               
      INTEGER ::  misl_2,nisl_2                                              
      REAL    ::  zDepth                                                       
      REAL    ::  d__eta,eta__1,eta__2,Khyd_1,Khyd_2                           
      REAL,PARAMETER  ::  RHsMin=  0.001        ! Min.Soil Relative Humidity   
      REAL    ::  PsiMax                        ! Max.Soil Water    Potential 
      REAL    ::  a_Khyd,b_Khyd                 ! Piecewis.Water Conductivity 


!c #WR REAL    ::  Khyd_x,Khyd_y                                               
                                                                              
                                                                              
                                                                 
!C +--Non Time Dependant SISVAT parameters                                    
!C +  ====================================                                
                                                                              
!C +--Soil Discretization                                                     
!C +  -------------------                                                     
                                                                              
!C +--Numerical Scheme Parameters                                             
!C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^                                             
        Implic = 0.75                           ! 0.5  <==> Crank-Nicholson  
        Explic = 1.00 - Implic                  !                            
                                                                             
!C +--Soil/Snow Layers Indices                                                
!C +  ^^^^^^^^^^^^^^^^^^^^^^^^                                                
      DO  isl=-nsol,0                                                          
        islpSV(isl) =           isl+1                                         
        islpSV(isl) = min(      islpSV(isl),0)                                
        islmSV(isl) =           isl-1                                         
        islmSV(isl) = max(-nsol,islmSV(isl))                                 
      END DO                                                                  
                                                                               
      DO  isn=1,nsno                                                           
        isnpSV(isn) =           isn+1                                          
        isnpSV(isn) = min(      isnpSV(isn),nsno)                            
      END DO                                                                 
                                                                              
!C +--Soil      Layers Thicknesses                                             
!C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
! Not used here as LMDZ method is applied, see SUBROUTINE get_soil_levels!    
!c #kd IF (nsol.gt.4)                                              THEN        
!c #kd   DO isl=-5,-nsol,-1                                                    
!c #kd     dz_dSV(isl)=   1.                                                  
!c #kd   END DO                                                                
!c #kd END IF                                                                 
!                                                                             
!      IF (nsol.ne.4)                                              THEN        
!        DO isl= 0,-nsol,-1                                                   
!          misl_2 =     -mod(isl,2)                                          
!          nisl_2 =         -isl/2                                            
!          dz_dSV(isl)=(((1-misl_2) * 0.001                                   
!     .                  +  misl_2  * 0.003) * 10**(nisl_2)) * 4.             
!C +...    dz_dSV(0)  =         Hapex-Sahel Calibration:       4 mm            
!                                                                             
!c +SH     dz_dSV(isl)=(((1-misl_2) * 0.001                                   
!c +SH.                  +  misl_2  * 0.003) * 10**(nisl_2)) * 1.              
!                                                                             
!c #05     dz_dSV(isl)=(((1-misl_2) * 0.001                                    
!c #05.                  +  misl_2  * 0.008) * 10**(nisl_2)) * 0.5             
!        END DO                                                                
!          dz_dSV(0)  =               0.001                                    
!          dz_dSV(-1) = dz_dSV(-1)  - dz_dSV(0)              + 0.004         
!      END IF                                                                 
                                                                               
        zz_dSV      = 0.                                                      
      DO  isl=-nsol,0                                                        
        dzmiSV(isl) = 0.500*(dz_dSV(isl)        +dz_dSV(islmSV(isl)))        
        dziiSV(isl) = 0.500* dz_dSV(isl)        /dzmiSV(isl)                  
        dzi_SV(isl) = 0.500* dz_dSV(islmSV(isl))/dzmiSV(isl)                   
        dtz_SV(isl) =        dt__SV             /dz_dSV(isl)                  
        dz78SV(isl) = 0.875* dz_dSV(isl)                                      
        dz34SV(isl) = 0.750* dz_dSV(isl)                                      
        dz_8SV(isl) = 0.125* dz_dSV(isl)                                       
        dzAvSV(isl) = 0.125* dz_dSV(islmSV(isl))                        &
     &              + 0.750* dz_dSV(isl)                                &
     &              + 0.125* dz_dSV(islpSV(isl))                          
!c #ER   dz78SV(isl) =        dz_dSV(isl)                                     
!c #ER   dz34SV(isl) =        dz_dSV(isl)                                      
!c #ER   dz_8SV(isl) = 0.                                                      
!c #ER   dzAvSV(isl) =        dz_dSV(isl)                                      
        zz_dSV      = zz_dSV+dz_dSV(isl)                                      
      END DO                                                                  
      DO ikl=1,knon !v                                                          
        dzsnSV(ikl,0) =      dz_dSV(0)                                        
      END DO                                                                  
                                                                              
!C +--Conversion to a 50 m Swab Ocean Discretization                           
!C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                           
        OcndSV = 0.                                                           
      DO isl=-nsol,0                                                          
        OcndSV = OcndSV +dz_dSV(isl)                                           
      END DO                                                                   
        OcndSV = 50.    /OcndSV                                                
                                                                               
                                                                               
!C +--Secondary Vegetation Parameters                                          
!C +  -------------------------------                                          
!                                                                              
!C +--Minimum Stomatal Resistance (Hapex Sahel Data)                           
!C +  (Taylor et al. 1997, J.Hydrol 188-189, p.1047)                           
!C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                           
!c #SH DO ivg=1,3                       !                                     
!c #SH   StodSV(ivg) = 210.             ! Millet                               
!c #SH END DO                           !                                      
!c #SH   StodSV(  4) = 120.             ! Sparse Tiger Bush                    
!c #SH DO ivg=5,6                       !                                     
!c #SH   StodSV(ivg) =  80.             ! Dense  Tiger Bush                   
!c #SH END DO                           !                                     
!c #SH   StodSV(  7) =  80.             ! Low    Trees (Fallow)               
!c #SH   StodSV( 10) =  80.             !                                     
!                                                                             
!C +--Minimum Stomatal Resistance (Tropical Forest)                          
!C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                           
!c #SH   StodSV(  8) =  60.             ! Medium Trees                        
!c #SH   StodSV( 11) =  60.             !                                     
!c #SH   StodSV(  9) =  40.             ! High   Trees                        
!c #SH   StodSV( 12) =  40.             !                                     
                                                                              
!C +--Root Fraction                                                           
!C +  ^^^^^^^^^^^^^                                                          
!C +  * GENERAL REFERENCE                                                     
!C +    Jackson et al., 1996: A global analysis of root distributions for     
!C +    terrestrial biomes. In Oecologia, 108, 389-411.                       
                                                                             
!C +  * ROOT PROFILE                                                           
!C +    The cumulative root fraction Y is given by                             
!C +        Y = 1 - beta**d   with d    the depth (in cm),                    
!C +                               beta a coefficient (vegetation dependent).
                                                                             
!C +  * BETA VALUES (for 11 world biomes)                                      
!C +  1  boreal forest                0.943                                    
!C +  2  crops                        0.961                                    
!C +  3  desert                       0.975                                   
!C +  4  sclerophyllous shrubs        0.964                                   
!C +  5  temperate coniferous forest  0.976                                   
!C +  6  temperate deciduous forest   0.966                                   
!C +  7  temperate grassland          0.943                                    
!C +  8  tropical deciduous forest    0.961                                    
!C +  9  tropical evergreen forest    0.962                                   
!C +  10 tropical grassland savanna   0.972                                    
!C +  11 tundra                       0.914                                   
!                                                                              
!C +  * ADVISED BETA VALUES FOR MAR                                           
!C +    (see 'block data SISVAT_dat', variable rbtdSV)                         
!C +                                                                           
!C +    SVAT veg. type         default      West Africa                        
!C +    0  barren soil         0.000        0.000                              
!C +    1  crops low           0.961 (2)    0.961 (2)                          
!C +    2  crops medium        0.961 (2)    0.961 (2)                          
!C +    3  crops high          0.961 (2)    0.961 (2)                          
!C +    4  grass low           0.943 (7)    0.943 (7)                          
!C +    5  grass medium        0.943 (7)    0.964 (4)                         
!C +    6  grass high          0.943 (7)    0.972 (10)                        
!C +    7  broadleaf low       0.966 (6)    0.968 (4,10)                      
!C +    8  broadleaf medium    0.966 (6)    0.962 (8,9)                       
!C +    9  broadleaf high      0.966 (6)    0.962 (8,9)                        
!C +    10 needleleaf low      0.976 (5)    0.971 (5,6)                       
!C +    11 needleleaf medium   0.976 (5)    0.976 (5)                         
!C +    12 needleleaf high     0.976 (5)    0.976 (5)                          
                                                                             
!C +    Numbers between brackets refer to Jackson's biomes. For more details  
!C +    about some choices, see the correspondance between the IGBP and SVAT   
!C +    vegetation classes (i.e. in NESTOR).                                 
                                                                              
!C +  * WARNING                                                                
!C +    Most of the roots are located in the first 2 m of soil. The root       
!C +    fraction per layer depends on the definition of the soil layer         
!C +    thickness. It will get wrong if a thick layer is defined around 2 m    
!C +    deep.                                                                  
                                                                               
      write(*,'(/a)') 'ROOT PROFILES (Jackson, 1996) :'                        
                                                                              
      DO ivt = 0, nvgt                                                        
        zDepth = 0.                                                           
        DO isl = 0, -nsol, -1                                                 
          IF (ivt .ne. 0) THEN                                                
            RF__SV(ivt,isl) =  rbtdSV(ivt)**zDepth *                     &
     &                         (1. - rbtdSV(ivt)**(dz_dSV(isl)*100) )        
            zDepth = zDepth + dz_dSV(isl)*100  !in cm                        
          ELSE                                                                 
            RF__SV(ivt,isl) = 0.                                              
          END IF                                                              
        END DO                                                                
        write(*,'(a,i2,a,i3,a,99f10.5:)')                                &
     &       '  RF__SV(', ivt, ',', -nsol, ':0) =', RF__SV(ivt,:)            
      END DO                                                                  
!      write(6, format(                                                  &
!     &  '  NOTE: If root fraction is not close to 0  around 2 m deep,', &
!     &/,'        Then you should redefine the soil layer thicknesses.', &
!     &/,'        See the code for more details.'))                         
                                                                              
                                                                               
!C +--Secondary Soil       Parameters                                          
!C +  -------------------------------                                         
                                                                               
      DO  ist=0,nsot                                                           
         rocsSV(ist)=(1.0-etadSV(ist))*1.2E+6  ! Soil Contrib. to (ro c)_s    
         s1__SV(ist)=     bCHdSV(ist)          & ! Factor of (eta)**(b+2)     
     &  *psidSV(ist)     *Ks_dSV(ist)          & !    in DR97, Eqn.(3.36)     
     & /(etadSV(ist)**(   bCHdSV(ist)+3.))     !                              
         s2__SV(ist)=     Ks_dSV(ist)          & ! Factor of (eta)**(2b+3)     
     & /(etadSV(ist)**(2.*bCHdSV(ist)+3.))     !    in DR97, Eqn.(3.35)       
                                                                             
!C +--Soil Minimum Humidity (from a prescribed minimum relative Humidity)     
!C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
         Psimax = -(log(RHsMin))/7.2E-5        ! DR97, Eqn 3.15 Inversion   
         etamSV(ist) =  etadSV(ist)                                      &
     &         *(PsiMax/psidSV(ist))**(-min(10.,1./bCHdSV(ist)))             
      END DO                                                                  
         etamSV(12)  =  0.                                                   
                                                                             
!C +--Piecewise Hydraulic Conductivity Profiles                               
!C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                               
      DO   ist=0,nsot                                                          
!c #WR      write(6, format(' Type |    etaSat | No |    eta__1 |    eta__2 |', &
!c #WR&           '    Khyd_1 |    Khyd_x |    Khyd_2 |    Khyd_y |'        &
!c #WR&         /,' -----+-----------+----+-----------+-----------+',       &
!c #WR&           '-----------+-----------+-----------+-----------+'))        
                                                                             
          d__eta          =  etadSV(ist)/nkhy                                 
          eta__1          =  0.                                              
          eta__2          =  d__eta                                           
        DO ikh=0,nkhy                                                         
          Khyd_1          =  s2__SV(ist)             & ! DR97, Eqn.(3.35)      
     &  *(eta__1      **(2. *bCHdSV(ist)+3.))        !                       
          Khyd_2          =  s2__SV(ist)             &!                       
     &  *(eta__2      **(2. *bCHdSV(ist)+3.))        !                       
                                                                             
          a_Khyd          = (Khyd_2-Khyd_1)/d__eta   !                       
          b_Khyd          =  Khyd_1-a_Khyd *eta__1   !                       
!c #WR     Khyd_x          =  a_Khyd*eta__1 +b_Khyd   !                       
!c #WR     Khyd_y          =  a_Khyd*eta__2 +b_Khyd   !                       
          aKdtSV(ist,ikh) =  a_Khyd       * dt__SV   !                       
          bKdtSV(ist,ikh) =  b_Khyd       * dt__SV   !                         
!c #WR     write(6,format(i5,' |',e10.2,' |',i3,' |', 6(e10.2,' |'))      & 
!c #WR&           ist,etadSV(ist),ikh,eta__1,                             &
!c #WR&          eta__2,Khyd_1,Khyd_x,Khyd_2,Khyd_y                       
                                                                          
          eta__1          = eta__1  + d__eta                              
          eta__2          = eta__2  + d__eta                              
        END DO                                                             
      END DO                                                               


!c
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + INITIALISATION: ARRAYS on grid  +++ 
! + -----------------------------------
! +                                                      
        DO ikl=1,knon                                

! + Water Fluxes, roughness length  
          us__SV(ikl) = 0.25                             ! Frict. Velocity  
                                      ! (approx ETH-Camp Lefebre etal CD 2005
          uts_SV(ikl) = 1.                               ! u*T*  arbitrary  
          uqs_SV(ikl) = 0.25 ! turb_vel(ikl,3)           ! u*q*    " 
          uss_SV(ikl) = 0.25 ! turb_vel(ikl,4)           ! u*s*    "
!!c #AE   usthSV(ikl) = 0. 
          Z0h_SV(ikl) = 0. ! rlength(ikl,9)
          Z0m_SV(ikl) = 0.                               ! Moment.Roughn.L.  
          Z0h_SV(ikl) = 0.                               ! Heat   Roughn.L.  
!!c #OR   Z0roSV(ikl) = 0.                               ! Orogr. Roughn.L.
          LMO_SV(ikl) = 0.
          VVs_SV(ikl) = 0.
          RRs_SV(ikl) = 0.
          DDs_SV(ikl) = 0.

! + Snow
          isnoSV(ikl) = 0                        ! # snow layers 
          BufsSV(ikl) = 0.    
          zzsnsv(ikl,:) = 0.                         ! Snow pack thickness
          qsnoSV(ikl) = 0.                        ! BL snow content  
          zWEcSV(ikl) = 0.
          dbs_SV(ikl) = 0.
          dsnbSV(ikl) = 0.
          esnbSV(ikl) = 0.
          BrosSV(ikl) = 0.
          BG1sSV(ikl) = 0.          
          BG2sSV(ikl) = 0.
          SWS_SV(ikl) = 0.
          RnofSV(ikl) = 0.                             ! RunOFF Intensity 
          toicSV(ikl) = 0.         !hj preli
! + Clouds !hj1602
          cld_SV(ikl) = 0.            

! + Vegetation
          LSmask(ikl) = 1                          ! Land/Sea   Mask    
          isotSV(ikl) = 12                         ! Soil       Type    
          iWaFSV(ikl) = 1                          ! Soil Drainage          
          Rootsv(ikl,:) = 0.
          ivgtSV(ikl) = 0                          ! Vegetation Type    
          LAI0SV(ikl) = 0                          ! LAI                
          glf0SV(ikl) = 0                          ! Green Leaf Frac.
!          TVegSV(ikl) = 0. 
          SnCaSV(ikl) = 0. 
          rrCaSV(ikl) = 0. 
          pSivSV(ikl) = 0. 

! + z0(Sastrugi)                                        
! #SZ     Z0SaSL(ikl) = 0.                                           
! #ZM   DO nt=1,ntavSL                                                      
! #ZM     SLn_z0(ikl,nt) = 0.5e-6                                       
! #ZM     SLn_b0(ikl,nt) = 0.5e-6                                       
! #ZM     SLn_r0(ikl,nt) = 0.5e-6                                       
! #ZM   END DO 


! +--z0(Orography Roughness)                                 
! #OR     SL_z0 (ikl) = min(SL_z0 (ikl),zsigma/3.)                    
! #OR     SLzoro(ikl) = min(SLzoro(ikl),zsigma/3.)                    

! +--SBL  Characteristics                                                          
! #AA   DO nt=  1,ntaver                                                     
! #AW     V_0aSL(ikl,nt)  = ssvSL(ikl)    
! #AH     dT0aSL(ikl,nt)  = temp_air(ikl)-tsurf(ikl)                     
! #AA   END DO                                                               


        END DO


                                                                           
      return                                                               

  END SUBROUTINE SISVAT_ini







!***************************************************************************

    SUBROUTINE sisvatetat0 (fichnom,ikl2i)

    USE dimphy
    USE mod_grid_phy_lmdz
    USE mod_phys_lmdz_para

    USE iostart
    USE VAR_SV
    USE VARdSV
    USE VARxSV          
    USE VARtSV
    USE indice_sol_mod

      IMPLICIT none
!======================================================================
! Auteur(s) HJ PUNGE (LSCE) date: 07/2009
! Objet: Ecriture de l'etat de redemarrage pour SISVAT
!======================================================================
    include "netcdf.inc"
!    include "indicesol.h"

    include "dimsoil.h"
    include "clesphys.h"
    include "thermcell.h"
    include "compbl.h"

!======================================================================
    CHARACTER(LEN=*) :: fichnom


    INTEGER, DIMENSION(klon), INTENT(IN) :: ikl2i
    REAL, DIMENSION(klon) :: rlon
    REAL, DIMENSION(klon) :: rlat

! les variables globales ecrites dans le fichier restart
    REAL, DIMENSION(klon) :: isno
    REAL, DIMENSION(klon) :: ispi
    REAL, DIMENSION(klon) :: iice
    REAL, DIMENSION(klon) :: rusn 
    REAL, DIMENSION(klon, nsno) :: isto

    REAL, DIMENSION(klon, nsismx) :: Tsis
    REAL, DIMENSION(klon, nsismx) :: eta
    REAL, DIMENSION(klon, nsismx) :: ro 

    REAL, DIMENSION(klon, nsno) :: dzsn      
    REAL, DIMENSION(klon, nsno) :: G1sn
    REAL, DIMENSION(klon, nsno) :: G2sn
    REAL, DIMENSION(klon, nsno) :: agsn

    REAL, DIMENSION(klon) :: toic

!    REAL, DIMENSION(klon) :: IRs
!    REAL, DIMENSION(klon) :: LMO
!    REAL, DIMENSION(klon) :: Bufs
!    REAL, DIMENSION(klon, 9) :: rlength
!    REAL, DIMENSION(klon, 5) :: turb_vel

    INTEGER  :: isl, ikl, i, isn , errT, erreta, errro, errdz, snopts
    CHARACTER (len=2) :: str2
    LOGICAL :: found
  
    errT=0
    errro=0
    erreta=0
    errdz=0
    snopts=0
! Ouvrir le fichier contenant l'etat initial:

    CALL open_startphy(fichnom)

! Lecture des latitudes, longitudes (coordonnees):

      CALL get_field("latitude",rlat,found)
      CALL get_field("longitude",rlon,found)

      CALL get_field("n_snows", isno,found)
      IF (.NOT. found) THEN
        PRINT*, 'phyetat0: Le champ <n_snows> est absent'
        PRINT *, 'fichier startsisvat non compatible avec sisvatetat0'
      ENDIF

      CALL get_field("n_ice_top",ispi,found)
      CALL get_field("n_ice",iice,found)
      CALL get_field("surf_water",rusn,found)
!      IF (.NOT. found) THEN
!        PRINT*, 'phyetat0: Le champ <surf_water> est absent'
!        rusn(:)=0.  
!      ENDIF


      CALL get_field("to_ice",toic,found)
      IF (.NOT. found) THEN
        PRINT*, 'phyetat0: Le champ <to_ice> est absent'
        toic(:)=0.  
      ENDIF


!      CALL get_field("IR_soil",IRs,found)
!      CALL get_field("LMO",LMO,found)
!      CALL get_field("snow_buffer",Bufs,found)
!      DO i = 1, 5
!            WRITE(str2,'(i2.2)') i
!            CALL get_field("turb_veloc"//str2, &
!                          turb_vel(:,i),found)
!      ENDDO
!      DO i = 1, 9
!            WRITE(str2,'(i2.2)') i
!            CALL get_field("rough_length"//str2, &
!                          rlength(:,i),found)
!      ENDDO
      DO isn = 1,nsno
        IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL get_field("AGESNOW"//str2, &
                          agsn(:,isn),found)
        ELSE
            PRINT*, "Trop de couches"
            CALL abort
        ENDIF
      ENDDO
      DO isn = 1,nsno
        IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL get_field("DZSNOW"//str2, &
                          dzsn(:,isn),found)
        ELSE
            PRINT*, "Trop de couches"
            CALL abort
        ENDIF
      ENDDO
      DO isn = 1,nsno
        IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL get_field("G2SNOW"//str2, &
                          G2sn(:,isn),found)
        ELSE
            PRINT*, "Trop de couches"
            CALL abort
        ENDIF
      ENDDO
      DO isn = 1,nsno
        IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL get_field("G1SNOW"//str2, &
                          G1sn(:,isn),found)
        ELSE
            PRINT*, "Trop de couches"
            CALL abort
        ENDIF
      ENDDO
      DO isn = 1,nsismx
        IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL get_field("ETA"//str2, &
                          eta(:,isn),found)
        ELSE
            PRINT*, "Trop de couches"
            CALL abort
        ENDIF
      ENDDO
      DO isn = 1,nsismx
        IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL get_field("RO"//str2, &
                          ro(:,isn),found)
        ELSE
            PRINT*, "Trop de couches"
            CALL abort
        ENDIF
      ENDDO
      DO isn = 1,nsismx
        IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL get_field("TSS"//str2, &
                          Tsis(:,isn),found)
        ELSE
            PRINT*, "Trop de couches"
            CALL abort
        ENDIF
      ENDDO
      DO isn = 1,nsno
        IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL get_field("HISTORY"//str2, &
                          isto(:,isn),found)
        ELSE
            PRINT*, "Trop de couches"
            CALL abort
        ENDIF
      ENDDO
      write(*,*)'Read ',fichnom,' finished!!'

!*********************************************************************************
! Compress restart file variables for SISVAT


      isnoSV(:)        = 0       
      ispiSV(:)        = 0 
      iiceSV(:)        = 0                                
       
      eta_SV(:,1:nsno) = 0.     
      TsisSV(:,1:nsno) = 0.         
      istoSV(:,:)      = 0 
      ro__SV(:,1:nsno) = 0.        
      G1snSV(:,:)      = 0.    
      G2snSV(:,:)      = 0. 
      dzsnSV(:,1:nsno) = 0.        
      agsnSV(:,:)      = 0.
      rusnSV(:)        = 0.   
      toicSV(:)        = 0.

      DO  ikl = 1,klon                                                   
          i   = ikl2i(ikl)   
          IF (i > 0) THEN
              isnoSV(ikl)     = INT(isno(i))          ! Nb Snow/Ice Lay.   
              ispiSV(ikl)     = INT(ispi(i))          ! Nb Supr.Ice Lay.  
              iiceSV(ikl)     = INT(iice(i))          ! Nb      Ice Lay.   
                                                                              
!!              IRs_SV(ikl)     = IRs(i)
!              LMO_SV(ikl)     = LMO(i)               !?              
!              us__SV(ikl)     = turb_vel(i,1)        
!              uts_SV(ikl)     = turb_vel(i,2)
!              uqs_SV(ikl)     = turb_vel(i,3) 
!              uss_SV(ikl)     = turb_vel(i,4) 
!              usthSV(ikl)     = turb_vel(i,5)
!              Z0m_SV(ikl)     = rlength(i,1)         ! Moment.Roughn.L. 
!              Z0mmSV(ikl)     = rlength(i,2)         ! Moment.Roughn.L.   
!              Z0mnSV(ikl)     = rlength(i,3)         ! Moment.Roughn.L.   
!              Z0SaSV(ikl)     = rlength(i,4)         ! Moment.Roughn.L.   
!              Z0e_SV(ikl)     = rlength(i,5)         ! Moment.Roughn.L.   
!              Z0emSV(ikl)     = rlength(i,6)         ! Moment.Roughn.L.   
!              Z0enSV(ikl)     = rlength(i,7)         ! Moment.Roughn.L.   
!              Z0roSV(ikl)     = rlength(i,8)         ! Moment.Roughn.L.   
!!              Z0h_SV(ikl)     = rlength(i,9)         ! Moment.Roughn.L.   

            DO isl =   -nsol,0    
              ro__SV(ikl,isl) = ro(i,nsno+1-isl)       !                    
              eta_SV(ikl,isl) = eta(i,nsno+1-isl)         ! Soil Humidity      
!hjp 15/10/2010
              IF (eta_SV(ikl,isl) <= 1.e-6) THEN          !hj check
                eta_SV(ikl,isl) = 1.e-6
              ENDIF
              TsisSV(ikl,isl) = Tsis(i,nsno+1-isl)        ! Soil Temperature  
              IF (TsisSV(ikl,isl) <= 1.) THEN          !hj check
!                errT=errT+1
                TsisSV(ikl,isl) = 273.15
              ENDIF

            END DO        
            write(*,*)'Copy histo', ikl

 
              istoSV(ikl,0) = 0                     ! Snow     History 

              G1snSV(ikl,0) = 0.5                   !      [-]     
              G2snSV(ikl,0) = 3.                    ! [-] [0.0001 m] 
              dzsnSV(ikl,0) = dz_dSV(0)             !            [m] 
              agsnSV(ikl,0) = 0.                    !          [day]     
   
            DO  isn = 1,isnoSV(ikl) !nsno      
              snopts=snopts+1 
              IF (isto(i,isn) > 10.) THEN          !hj check
                write(*,*)'Irregular isto',ikl,i,isn,isto(i,isn)
                isto(i,isn) = 1.
              ENDIF

              istoSV(ikl,isn) = INT(isto(i,isn))     ! Snow     History 
              ro__SV(ikl,isn) = ro(i,isn)            !        [kg/m3]     
              eta_SV(ikl,isn) = eta(i,isn)           !        [m3/m3]  
              TsisSV(ikl,isn) = Tsis(i,isn)          !            [K]  

             IF (TsisSV(ikl,isn) <= 1.) THEN          !hj check
              errT=errT+1
              TsisSV(ikl,isn) = TsisSV(ikl,0)
             ENDIF  
             IF (TsisSV(ikl,isn) <= 1.) THEN          !hj check
              TsisSV(ikl,isn) = 263.15
             ENDIF
             IF (eta_SV(ikl,isn) < 1.e-9) THEN          !hj check
              eta_SV(ikl,isn) = 1.e-6  
              erreta=erreta+1
             ENDIF  
             IF (ro__SV(ikl,isn) <= 10.) THEN          !hj check
              ro__SV(ikl,isn) = 11.
              errro=errro+1
             ENDIF 
              write(*,*)ikl,i,isn,Tsis(i,isn),G1sn(i,isn)
              G1snSV(ikl,isn) = G1sn(i,isn)          ! [-]        [-]     
              G2snSV(ikl,isn) = G2sn(i,isn)          ! [-] [0.0001 m] 
              dzsnSV(ikl,isn) = dzsn(i,isn)          !            [m] 
!             IF (dzsnSV(ikl,isn) < 5.e-6) THEN          !hj check
!              dzsnSV(ikl,isn) = 0.000005
!              errdz=errdz+1
!             ENDIF         
              agsnSV(ikl,isn) = agsn(i,isn)          !          [day]     
            END DO  
              rusnSV(ikl)     = rusn(i)              ! Surficial Water    
              toicSV(ikl)     = toic(i)              ! bilan snow to ice   
!!              BufsSV(ikl)     = Bufs(i)    
          END IF                       
        END DO    
!        write(*,*)snopts,'snow pts',errT,' T errors,',errro,'ro,'
!        write(*,*)'   ',erreta,'eta,',errdz,'dz' 
    END SUBROUTINE sisvatetat0




    SUBROUTINE sisvatredem (fichnom,ikl2i,rlon,rlat)
!======================================================================
! Auteur(s) HJ PUNGE (LSCE) date: 07/2009
! Objet: Ecriture de l'etat de redemarrage pour SISVAT
!======================================================================
    USE mod_grid_phy_lmdz
    USE mod_phys_lmdz_para
    USE iostart
    USE VAR_SV
    USE VARxSV          
    USE VARySV !hj tmp 12 03 2010
    USE VARtSV
    USE indice_sol_mod

    IMPLICIT none

    include "netcdf.inc"
!    include "indicesol.h" 
    include "dimsoil.h"
    include "clesphys.h" 
    include "thermcell.h"
    include "compbl.h"

!======================================================================

    CHARACTER(LEN=*) :: fichnom
    INTEGER, DIMENSION(klon), INTENT(IN) :: ikl2i
    REAL, DIMENSION(klon), INTENT(IN) :: rlon
    REAL, DIMENSION(klon), INTENT(IN) :: rlat

! les variables globales ecrites dans le fichier restart
    REAL, DIMENSION(klon) :: isno
    REAL, DIMENSION(klon) :: ispi
    REAL, DIMENSION(klon) :: iice
    REAL, DIMENSION(klon, nsnowmx) :: isto

    REAL, DIMENSION(klon, nsismx) :: Tsis
    REAL, DIMENSION(klon, nsismx) :: eta
    REAL, DIMENSION(klon, nsnowmx) :: dzsn
    REAL, DIMENSION(klon, nsismx) :: ro        
    REAL, DIMENSION(klon, nsnowmx) :: G1sn
    REAL, DIMENSION(klon, nsnowmx) :: G2sn
    REAL, DIMENSION(klon, nsnowmx) :: agsn
    REAL, DIMENSION(klon) :: IRs
    REAL, DIMENSION(klon) :: LMO
    REAL, DIMENSION(klon) :: rusn 
    REAL, DIMENSION(klon) :: toic
    REAL, DIMENSION(klon) :: Bufs
    REAL, DIMENSION(klon) :: alb1,alb2,alb3
    REAL, DIMENSION(klon, 9) :: rlength
    REAL, DIMENSION(klon, 5) :: turb_vel

    INTEGER isl, ikl, i, isn
    CHARACTER (len=2) :: str2
    INTEGER           :: pass

      isno(:)       = 0       
      ispi(:)       = 0 
      iice(:)       = 0                                
      IRs(:)        = 0.
      LMO(:)        = 0.                             
      turb_vel(:,:) = 0.
      rlength(:,:)  = 0.                                  
      eta(:,:)      = 0.     
      Tsis(:,:)     = 0.         
      isto(:,:)     = 0 
      ro(:,:)       = 0.        
      G1sn(:,:)     = 0.    
      G2sn(:,:)     = 0. 
      dzsn(:,:)     = 0.        
      agsn(:,:)     = 0. 
      rusn(:)       = 0.   
      toic(:)       = 0.   
      Bufs(:)       = 0.   
      alb1(:)       = 0.
      alb2(:)       = 0.
      alb3(:)       = 0.

!***************************************************************************
! Uncompress SISVAT output variables for storage
           
      DO  ikl = 1,klon                                                   
        i   = ikl2i(ikl)
      IF (i > 0) THEN
        isno(i)       = 1.*isnoSV(ikl)               ! Nb Snow/Ice Lay.   
        ispi(i)       = 1.*ispiSV(ikl)               ! Nb Supr.Ice Lay.   
        iice(i)       = 1.*iiceSV(ikl)               ! Nb      Ice Lay.       
   
!        IRs(i)        = IRs_SV(ikl)
!        LMO(i)        = LMO_SV(ikl)                             
!        turb_vel(i,1) = us__SV(ikl) 
!        turb_vel(i,2) = uts_SV(ikl) 
!        turb_vel(i,3) = uqs_SV(ikl) 
!        turb_vel(i,4) = uss_SV(ikl) 
!        turb_vel(i,5) = usthSV(ikl) 
!        rlength(i,1)  = Z0m_SV(ikl)                  ! Moment.Roughn.L. 
!        rlength(i,2)  = Z0mmSV(ikl)                  !  
!        rlength(i,3)  = Z0mnSV(ikl)                  !  
!        rlength(i,4)  = Z0SaSV(ikl)                  !   
!        rlength(i,5)  = Z0e_SV(ikl)                  !    
!        rlength(i,6)  = Z0emSV(ikl)                  !   
!        rlength(i,7)  = Z0enSV(ikl)                  !    
!        rlength(i,8)  = Z0roSV(ikl)                  !   
!        rlength(i,9)  = Z0h_SV(ikl)                  !    

        DO isl =   -nsol,0                           !                    
          eta(i,nsno+1-isl)  = eta_SV(ikl,isl)            ! Soil Humidity      
          Tsis(i,nsno+1-isl) = TsisSV(ikl,isl)            ! Soil Temperature   
          ro(i,nsno+1-isl)   = ro__SV(ikl,isl)            !        [kg/m3]    
        END DO        
 
  
        DO  isn = 1,nsno              
          isto(i,isn)   = 1.*istoSV(ikl,isn)         ! Snow     History 
          ro(i,isn)     = ro__SV(ikl,isn)            !        [kg/m3]     
          eta(i,isn)    = eta_SV(ikl,isn)            !        [m3/m3]  
          Tsis(i,isn)   = TsisSV(ikl,isn)            !            [K]    
          G1sn(i,isn)   = G1snSV(ikl,isn)            ! [-]        [-]     
          G2sn(i,isn)   = G2snSV(ikl,isn)            ! [-] [0.0001 m] 
          dzsn(i,isn)   = dzsnSV(ikl,isn)            !            [m]         
          agsn(i,isn)   = agsnSV(ikl,isn)            !          [day]     
        END DO  
        rusn(i)       = rusnSV(ikl)                  ! Surficial Water  
        toic(i)       = toicSV(ikl)                  ! Surficial Water  
        alb1(i)       = alb1sv(ikl) 
        alb2(i)       = alb2sv(ikl)
        alb3(i)       = alb3sv(ikl)
!        Bufs(i)       = BufsSV(ikl)      
      END IF                     
      END DO                                                

      CALL open_restartphy(fichnom)
      DO pass = 1, 2
        CALL put_field(pass,"longitude", &
                    "Longitudes de la grille physique",rlon)     
        CALL put_field(pass,"latitude","Latitudes de la grille physique",rlat)
  
        CALL put_field(pass,"n_snows", "number of snow/ice layers",isno)
        CALL put_field(pass,"n_ice_top", "number of top ice layers",ispi)
        CALL put_field(pass,"n_ice", "number of ice layers",iice)
        CALL put_field(pass,"IR_soil", "Soil IR flux",IRs)
        CALL put_field(pass,"LMO", "Monin-Obukhov Scale",LMO)
        CALL put_field(pass,"surf_water", "Surficial water",rusn)
        CALL put_field(pass,"snow_buffer", "Snow buffer layer",Bufs)
        CALL put_field(pass,"alb_1", "albedo sw",alb1)
        CALL put_field(pass,"alb_2", "albedo nIR",alb2)
        CALL put_field(pass,"alb_3", "albedo fIR",alb3)
        CALL put_field(pass,"to_ice", "Snow passed to ice",toic)

 !     DO i = 1, 5
 !       WRITE(str2,'(i2.2)') i
 !       CALL put_field(pass,"turb_veloc"//str2, &
 !                      "various turbulent velocities"//str2, &
 !                      turb_vel(:,i))
 !     ENDDO
 !     DO i = 1, 9
 !       WRITE(str2,'(i2.2)') i
 !       CALL put_field(pass,"rough_length"//str2, &
 !                      "various roughness lengths"//str2, &
 !                      rlength(:,i))
 !     ENDDO
        DO isn = 1,nsno
          IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL put_field(pass,"AGESNOW"//str2, &
                         "Age de la neige layer No."//str2, &
                         agsn(:,isn))
          ELSE
            PRINT*, "Trop de couches"
            CALL abort
          ENDIF
        ENDDO
        DO isn = 1,nsno
          IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL put_field(pass,"DZSNOW"//str2, &
                         "Snow/ice thickness layer No."//str2, &
                         dzsn(:,isn))
          ELSE
            PRINT*, "Trop de couches"
            CALL abort
          ENDIF
        ENDDO
        DO isn = 1,nsno
          IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL put_field(pass,"G2SNOW"//str2, &
                         "Snow Property 2, layer No."//str2, &
                         G2sn(:,isn))
          ELSE
            PRINT*, "Trop de couches"
            CALL abort
          ENDIF
        ENDDO
        DO isn = 1,nsno
          IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL put_field(pass,"G1SNOW"//str2, &
                         "Snow Property 1, layer No."//str2, &
                         G1sn(:,isn))
          ELSE
            PRINT*, "Trop de couches"
            CALL abort
          ENDIF
        ENDDO
        DO isn = 1,nsismx
          IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL put_field(pass,"ETA"//str2, &
                         "Soil/snow water content layer No."//str2, &
                         eta(:,isn))
          ELSE
            PRINT*, "Trop de couches"
            CALL abort
          ENDIF
        ENDDO
        DO isn = 1,nsismx   !nsno
          IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL put_field(pass,"RO"//str2, &
                           "Snow density layer No."//str2, &
                           ro(:,isn))
          ELSE
            PRINT*, "Trop de couches"
            CALL abort
          ENDIF
        ENDDO
        DO isn = 1,nsismx
          IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL put_field(pass,"TSS"//str2, &
                           "Soil/snow temperature layer No."//str2, &
                           Tsis(:,isn))
          ELSE
            PRINT*, "Trop de couches"
            CALL abort
          ENDIF
        ENDDO
        DO isn = 1,nsno
          IF (isn.LE.99) THEN
            WRITE(str2,'(i2.2)') isn
            CALL put_field(pass,"HISTORY"//str2, &
                           "Snow history layer No."//str2, &
                           isto(:,isn))
          ELSE
            PRINT*, "Trop de couches"
            CALL abort
          ENDIF
        ENDDO
      ENDDO

  END SUBROUTINE sisvatredem

END MODULE surf_sisvat_mod
