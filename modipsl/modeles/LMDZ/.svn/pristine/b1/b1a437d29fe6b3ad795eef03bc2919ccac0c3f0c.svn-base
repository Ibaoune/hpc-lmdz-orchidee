      subroutine PHY_SISVAT_ALLOC

!------------------------------------------------------------------------------+
!                                                         Thu 27-Jun-2013  MAR |
!                                                                              |
!     subroutine PHY_SISVAT_ALLOC  allocates prognostic variables of           |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme           |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 26-Feb-2013      |
!           Last Modification by H. Gallee,               Thu 27-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


      use Mod_Real
      use Mod_PHY____grd
      use Mod_SISVAT_dim
      use Mod_SISVAT_grd
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
      use Mod_SISVAT_dzS
      use Mod_SISVAT_cdf
      use Mod_SISVAT_gpt

      use Mod_SISVAT_flx
! #AW use Mod_SISVAT_xAW
! #AH use Mod_SISVAT_xAH
! #ZM use Mod_SISVAT_xZM
! #WL use Mod_SISVAT_xWL


      IMPLICIT NONE


      integer  ::  ikl   ,i     ,j     ,n        !




! ================================
! Initialization of Mod_SISVAT_grd
! ================================

      n2    =                      min(2,mwp)    !
      kcolv =(mxpp-ixp1+1)*(mypp-jyp1+1)*mwp     !

      nvege = nvgt                               !      Nb of Vegetation Types
      nsnow = nsno                               !  MAX Nb of Snow Layers
      nsoil = nsol                               !      Nb of Soil Layers
      nbPts = nb_Pts                             !
      nbwri = nb_wri                             !
      ntave = ntaver                             !      Nb of Time Steps used in V & dT(a-s) Time Moving Average
      ntavz = ntavez                             !      Nb of Time Steps used in z0, r0, ... Time Moving Average
      nLimi = nLimit                             !      Nb of Time Steps used in Water Vapor Flux Limit. Average



! VERIFICATION: a BAD SET-UP WILL LEAD to an EMERGENCY STOP
! ---------------------------------------------------------

        IF (nbwri.gt.mzp)                                           THEN
          write(6,600) nbwri,mzp
 600      format(/,'### MAR_SISVAT CRASH, nbwri =',i6,                  &
     &                             ' .GT. mzp       =',i3,' ',2x,' ###',&
     &           /,'    ?!&~@|@[#@#]=!!!',23x,'EMERGENCY STOP')
          stop
        END IF

        IF (nbwri.gt.mwp*NbPts)                                     THEN
          write(6,601) nb_wri,mwp,Nb_Pts
 601      format(/,'### MAR_SISVAT CRASH, nbwri =',i6,                  &
     &                             ' .GT. mwp*NbPts =',i3,'*',i2,' ###',&
     &           /,'    ?!&~@|@[#@#]=!!!',23x,'EMERGENCY STOP')
          stop
        END IF





! =================================
! ALLOCATION Mod_SISVAT_kkl - BEGIN
! =================================

! SISVAT INPUT        Variables
! -----------------------------

      allocate            ( LSmask(kcolp,mwp) )               ! Land-Sea   Mask
      allocate            ( isotSV(kcolp,mwp) )               ! Soil       Type
      allocate            ( iWaFSV(kcolp,mwp) )               ! Soil       Drainage:(1,0)=(y,n)
      allocate            ( ivgtSV(kcolp,mwp) )               ! Vegetation Type

      allocate            ( FracSV(kcolp,mwp) )               !  Grid Cell Fraction (Mosaic)                             [-]

      allocate            ( coszSV(kcolp,mwp) )               ! Cosine of Sun zenithal Angle
      allocate            ( sol_SV(kcolp,mwp) )               ! Downward  Solar    Radiation
      allocate            ( IRd_SV(kcolp,mwp) )               ! Downward  Longwave Radiation

      allocate            ( drr_SV(kcolp,mwp) )               ! Rain  Intensity                                    [kg/m2/s]
      allocate            ( dsn_SV(kcolp,mwp) )               ! Snow  Intensity                                    [kg/m2/s]
      allocate            ( dsnbSV(kcolp,mwp) )               ! Idem, fraction, from Drift                               [-]
      allocate            ( esnbSV(kcolp,mwp) )               ! Idem, fraction, from Drift                               [-]
      allocate            ( dbs_SV(kcolp,mwp) )               ! Drift Amount                                         [kg/m2]
      allocate            ( BrosSV(kcolp,mwp) )               ! Buffer Snow Layer Density
      allocate            ( BG1sSV(kcolp,mwp) )               ! Buffer Snow Layer Dendr/Sphe                             [-]
      allocate            ( BG2sSV(kcolp,mwp) )               ! Buffer Snow Layer Spher/Size                             [-] [0.0001 m]
      allocate            ( dz0_SV(kcolp,mwp) )               ! dz0(Sastrugi dh)                                         [m]

      allocate            ( cld_SV(kcolp,mwp) )               ! Cloudiness (seen from SBL)
      allocate            ( za__SV(kcolp,mwp) )               ! SBL Height
      allocate            ( VV__SV(kcolp,mwp) )               !(SBL Top)  Wind Velocity
      allocate            ( Ua__SV(kcolp,mwp) )               !(SBL Top)  Wind Velocity, x-Direction, t                [m/s]
      allocate            ( Ua0_SV(kcolp,mwp) )               !(SBL Top)  Wind Velocity, x-Direction, t -dt            [m/s]
      allocate            ( Va__SV(kcolp,mwp) )               !(SBL Top)  Wind Velocity, y-Direction, t                [m/s]
      allocate            ( Va0_SV(kcolp,mwp) )               !(SBL Top)  Wind Velocity, y-Direction, t -dt            [m/s]
      allocate            ( VV10SV(kcolp,mwp) )               ! 10-m      Wind Velocity
      allocate            ( VVs_SV(kcolp,mwp) )               !(Sastr,V)  Relevance
      allocate            ( RRsxSV(kcolp,mwp) )               !(Sastr,V)  Counter
      allocate            ( DDsxSV(kcolp,mwp) )               !(Sastr,V)  Angle
      allocate            ( DDs_SV(kcolp,mwp) )               !(Sastr,V)  Angle
      allocate            ( rhT_SV(kcolp,mwp) )               ! SBL Top   Air  Density
      allocate            ( TaT_SV(kcolp,mwp) )               ! SBL Top   Temperature                                    [K]
      allocate            ( Ts__SV(kcolp,mwp) )               ! Surface   Air Temperature (Mosaic)                       [K]
      allocate            ( pkPaSV(kcolp,mwp) )               ! Surface   Pressure                                     [kPa]
      allocate            ( WindSV(kcolp,mwp,mzp) )           ! Wind         Speed                                     [m/s]
      allocate            ( zza_SV(kcolp,mwp,mzp) )           ! Atmospheric  Levels HEIGHTS                              [m]
      allocate            ( roa_SV(kcolp,mwp,mzp) )           ! Air       Volumic   Mass                              [T/m3]
      allocate            ( Kz__SV(kcolp,mwp,mzp) )           ! Turbulent Diffusion Coefficients                      [m2/s]
      allocate            ( pktaSV(kcolp,mwp,mzp) )           ! Temperature / Exner Potential (Current Value in (PHY_)SISVAT)
      allocate            ( pkt0SV(kcolp,mwp,mzp) )           ! Temperature / Exner Potential (INPUT   Value of (PHY_)SISVAT)
      allocate            ( ExnrSV(kcolp,mwp) )               ! Surface       Exner Potential
      allocate            ( qv__SV(kcolp,mwp,mzp) )           ! Atmosph.  Specific Humidity                          [kg/kg]
      allocate            ( QaT_SV(kcolp,mwp) )               ! SBL Top   Specific Humidity                          [kg/kg]
      allocate            ( dQa_SV(kcolp,mwp) )               ! SBL Flux  Limitation of Qa
      allocate            ( SHumSV(kcolp,mwp) )               ! Surface   Specific Humidity                          [kg/kg]
      allocate            ( dSdTSV(kcolp,mwp) )               ! Sensible Heat Flux T Derivat.
      allocate            ( dLdTSV(kcolp,mwp) )               ! Latent   Heat Flux T Derivat.
      allocate            ( qsnoSV(kcolp,mwp) )               ! SBL Mean  Snow       Content                         [kg/kg]

      allocate            ( LAI0SV(kcolp,mwp) )               ! Nominal Leaf Area Index
      allocate            ( glf0SV(kcolp,mwp) )               ! Green   Leaf Fraction

      allocate            ( alb0SV(kcolp,mwp) )               ! Soil    Albedo
      allocate            ( slopSV(kcolp,mwp) )               ! Snow/Ice/Soil-Water Surf. Slope                          [-]
      allocate            ( slorSV(kcolp,mwp) )               ! Snow/Ice/Soil-Water Surf. Slope                     [radian]
      allocate            ( z0__SV(kcolp,mwp) )               ! Roughness Length Momentum   (Mosaic)                     [m]
      allocate            ( ROF_SV(kcolp,mwp) )               ! Cumulative Run-Off          (Mosaic)               [mm w.e.]


! SISVAT INPUT/OUTPUT Variables
! -----------------------------

      allocate            ( isnoSV(kcolp,mwp) )               ! Nb of Ice/Snow Layers
      allocate            ( ispiSV(kcolp,mwp) )               ! Uppermost superimposed ice
      allocate            ( iiceSV(kcolp,mwp) )               ! Nb of Ice      Layers
      allocate            ( istoSV(kcolp,mwp,    0:nsnow) )   ! Snow Layer     History

      allocate            ( albcSV(kcolp,mwp) )               ! Coupl. Surface Albedo (Surface-Canopy / Ocean)
      allocate            ( alb_SV(kcolp,mwp) )               ! Surface-Canopy Albedo
      allocate            ( emi_SV(kcolp,mwp) )               ! Surface-Canopy Emissivity
      allocate            ( IRs_SV(kcolp,mwp) )               ! Soil           IR Flux
      allocate            ( LMO_SV(kcolp,mwp) )               ! Monin-Obukhov  Scale
      allocate            ( us__SV(kcolp,mwp) )               ! Friction       Velocity
      allocate            ( uts_SV(kcolp,mwp) )               ! Temperature  Turbulent Scale
      allocate            ( cutsSV(kcolp,mwp) )               ! Temperature  Turbulent Scale C.
      allocate            ( uqs_SV(kcolp,mwp) )               ! Spec.Humid.  Turbulent Scale
      allocate            ( ussbSV(kcolp,mwp) )               ! Blowing Snow Eroded    Buffer                      [mm w.e.]
      allocate            ( uss_SV(kcolp,mwp) )               ! Blowing Snow Turbulent Scale
! #BS allocate            ( ussxSV(kcolp,mwp) )               ! Blowing Snow Turbulent Scale    (modified)
! #BD allocate            ( uds_SV(kcolp,mwp) )               ! Blowing Dust Flux Turbulent Scale                [kg/kg m/s]
      allocate            ( usthSV(kcolp,mwp) )               ! Blowing Snow Erosion Thresh.
      allocate            ( rCDmSV(kcolp,mwp) )               ! Square  Root Contribut. Drag_m
      allocate            ( rCDhSV(kcolp,mwp) )               ! Square  Root Contribut. Drag_h
      allocate            ( Z0m_SV(kcolp,mwp) )               ! Momentum     Roughness Length
      allocate            ( Z0mmSV(kcolp,mwp) )               !  z0(Momentum,    Time Mean)                              [m]
      allocate            ( Z0mnSV(kcolp,mwp) )               !  z0(Momentum,    instanta.)                              [m]
      allocate            ( Z0roSV(kcolp,mwp) )               ! Subgrid Topo Roughness Length
      allocate            ( Z0SaSV(kcolp,mwp) )               !  z0(Sastrugi  h)                                         [m]
      allocate            ( Z0e_SV(kcolp,mwp) )               !  z0(Snow eroded)                                         [m]
      allocate            ( Z0emSV(kcolp,mwp) )               !  z0(Snow eroded, Time Mean)                              [m]
      allocate            ( Z0enSV(kcolp,mwp) )               !  z0(Snow eroded, instanta.)                              [m]
      allocate            ( Z0h_SV(kcolp,mwp) )               ! Heat         Roughness Length
      allocate            ( Z0hmSV(kcolp,mwp) )               !  z0(Heat,        Time Mean)                              [m]
      allocate            ( Z0hnSV(kcolp,mwp) )               !  z0(Heat,        instanta.)                              [m]

      allocate            ( snCaSV(kcolp,mwp) )               ! Canopy  Snow   Thickness
      allocate            ( rrCaSV(kcolp,mwp) )               ! Canopy  Water  Content
      allocate            ( psivSV(kcolp,mwp) )               ! Leaf    Water  Potential 
      allocate            ( TvegSV(kcolp,mwp) )               ! Vegetation     Temperature

      allocate            ( TsisSV(kcolp,mwp,-nsoil:nsnow) )  ! Snow/Ice/Soil-Water Temperature
      allocate            ( ro__SV(kcolp,mwp,-nsoil:nsnow) )  ! Snow/Ice/Soil-Water VolumicMass
      allocate            ( eta_SV(kcolp,mwp,-nsoil:nsnow) )  ! Snow/Ice/Soil     Water Content
      allocate            ( G1snSV(kcolp,mwp,     0:nsnow) )  ! Snow Dendricity/Sphericity
      allocate            ( G2snSV(kcolp,mwp,     0:nsnow) )  ! Snow Sphericity/Size
      allocate            ( dzsnSV(kcolp,mwp,     0:nsnow) )  ! Snow Layer  Thickness
      allocate            ( agsnSV(kcolp,mwp,     0:nsnow) )  ! Snow Age
      allocate            ( BufsSV(kcolp,mwp) )               ! Snow Buffer Layer
      allocate            ( rusnSV(kcolp,mwp) )               ! Surficial   Water
      allocate            ( SWf_SV(kcolp,mwp) )               ! Normalized  Decay
      allocate            ( SWS_SV(kcolp,mwp) )               ! Surficial Water Status
      allocate            ( HFraSV(kcolp,mwp) )               ! Frazil      Thickness

      allocate            ( zWE_SV(kcolp,mwp) )               ! Current   Snow Thickness                 [mm w.e.]
      allocate            ( zWEcSV(kcolp,mwp) )               ! Compacted Snow Thickness                 [mm w.e.]
      allocate            ( dwemSV(kcolp,mwp) )               ! Only Melting  over dt__SV   (Mosaic)     [mm w.e.]
      allocate            ( dwerSV(kcolp,mwp) )               ! Refreezing    over dt__SV   (Mosaic)     [mm w.e.]
      allocate            ( dwesSV(kcolp,mwp) )               ! Sublimation   over dt__SV   (Mosaic)     [mm w.e.]
      allocate            ( wem0SV(kcolp,mwp) )               ! Only Melting       Budget   (Mosaic)     [mm w.e.]
      allocate            ( wem_SV(kcolp,mwp) )               ! Only Melting       Budget   (Mosaic)     [mm w.e.]
      allocate            ( wer0SV(kcolp,mwp) )               ! Refreezing         Budget   (Mosaic)     [mm w.e.]
      allocate            ( wer_SV(kcolp,mwp) )               ! Refreezing         Budget   (Mosaic)     [mm w.e.]
      allocate            ( wes0SV(kcolp,mwp) )               ! Sublimation        Budget   (Mosaic)     [mm w.e.]
      allocate            ( wes_SV(kcolp,mwp) )               ! Sublimation        Budget   (Mosaic)     [mm w.e.]
      allocate            ( wee_SV(kcolp,mwp) )               ! Evapotranspiration Budget   (Mosaic)     [mm w.e.]


! SISVAT OUTPUT       Variables
! -----------------------------

      allocate            ( no__SV(nbwri) )                   ! OUTPUT file Unit Number
      allocate            ( IOi_SV(NbPts) )                   ! OUTPUT point   i Coordinate (independant txt file)
      allocate            ( IOj_SV(NbPts) )                   ! OUTPUT point   j Coordinate (independant txt file)
      allocate            ( i___SV(nbwri) )                   ! OUTPUT point   i Coordinate
      allocate            ( j___SV(nbwri) )                   ! OUTPUT point   j Coordinate
      allocate            ( n___SV(nbwri) )                   ! OUTPUT point   n Coordinate
      allocate            ( lwriSV(kcolp,mwp) )               ! OUTPUT point vec Index

      allocate            ( IRu_SV(kcolp,mwp) )               ! UPward    IR Flux (effective)
      allocate            ( hSalSV(kcolp,mwp) )               ! Saltating Layer Height
      allocate            ( qSalSV(kcolp,mwp) )               ! Saltating Snow  Concentration
      allocate            ( RnofSV(kcolp,mwp) )               ! RunOFF    Intensity

! =================================
! ALLOCATION Mod_SISVAT_kkl -   END
! =================================



! =================================
! ALLOCATION Mod_SISVAT_grd - BEGIN
! =================================

      allocate            ( ii__SV(kcolv)       )             ! Mosaic point   i Coordinate
      allocate            ( jj__SV(kcolv)       )             ! Mosaic point   j Coordinate
      allocate            ( nn__SV(kcolv)       )             ! Mosaic point   n Coordinate
      allocate            ( ikp_SV(kcolv)       )             ! Grid Cell Column Index of a SISVAT Column
      allocate            ( ikl_SV(mxp,myp,mwp) )             ! SISVAT    Column Index 

! =================================
! ALLOCATION Mod_SISVAT_grd -   END
! =================================



! =================================
! ALLOCATION Mod_SISVAT_loc - BEGIN
! =================================


      allocate            ( NLaysv(kcolp,mwp) )               ! New   Snow     Layer   Switch
      allocate            ( i_thin(kcolp,mwp) )               ! Index of the thinest Layer
      allocate            ( LIndsv(kcolp,mwp) )               ! Contiguous Layer relative Index

      allocate            ( albisv(kcolp,mwp) )               ! Integrated Surface Albedo
      allocate            ( albssv(kcolp,mwp) )               ! Soil               Albedo [-]
      allocate            ( SoCasv(kcolp,mwp) )               ! Canopy  Absorbed Solar Radiat.
      allocate            ( SoSosv(kcolp,mwp) )               ! Surface Absorbed Solar Radiat.
      allocate            ( IRv_sv(kcolp,mwp) )               ! Vegetation IR Flux  [W/m2]
      allocate            ( Evg_sv(kcolp,mwp) )               ! Emissivity of Vegetation+Snow
      allocate            ( Eso_sv(kcolp,mwp) )               ! Emissivity of       Soil+Snow
      allocate            ( tau_sv(kcolp,mwp) )               ! Transmited Radiation Fraction
      allocate            ( rrMxsv(kcolp,mwp) )               ! Canopy Maximum Intercepted Rain
      allocate            ( LAIesv(kcolp,mwp) )               ! effective LAI for transpirati.
      allocate            ( LAI_sv(kcolp,mwp) )               ! corrected LAI in case of snow
      allocate            ( glf_sv(kcolp,mwp) )               ! Green  Leaf Fraction
      allocate            ( Sigmsv(kcolp,mwp) )               ! Canopy Ventilation  Factor
      allocate            ( HSv_sv(kcolp,mwp) )               ! Sensible Heat Flux  [W/m2]
      allocate            ( HLv_sv(kcolp,mwp) )               ! Latent   Heat Flux  [W/m2]
      allocate            ( HSs_sv(kcolp,mwp) )               ! Sensible Heat Flux (t)
      allocate            ( HLs_sv(kcolp,mwp) )               ! Latent   Heat Flux (t)
      allocate            ( sqrCm0(kcolp,mwp) )               ! in Neutral Drag Coef.Moment.
      allocate            ( sqrCh0(kcolp,mwp) )               ! in Neutral Drag Coef.Heat
      allocate            ( Lx_H2O(kcolp,mwp) )               ! Latent Heat of Vaporiz./Sublim.
      allocate            ( ram_sv(kcolp,mwp) )               ! Aerodyn.Resistance (Moment.)
      allocate            ( rah_sv(kcolp,mwp) )               ! Aerodyn.Resistance (Heat)
      allocate            ( Fh__sv(kcolp,mwp) )               ! Stability Function
      allocate            ( dFh_sv(kcolp,mwp) )               ! Stability Function (Deriv.)
      allocate            ( Evp_sv(kcolp,mwp) )               ! Evaporation        [kg/m2]
      allocate            ( EvT_sv(kcolp,mwp) )               ! Evapotranspiration [kg/m2]
      allocate            ( LSdzsv(kcolp,mwp) )               ! Land/Sea Vert. Discretiz. Fact.
      allocate            ( Tsrfsv(kcolp,mwp) )               ! Surface    Temperature
      allocate            ( sEX_sv(kcolp,mwp,-nsoil:nsnow+1) )! Verticaly Integr.Extinct.Coef.
      allocate            ( zzsnsv(kcolp,mwp,     0:nsnow)   )! Snow  Pack Thickness      [m]
      allocate            ( psi_sv(kcolp,mwp,-nsoil:0   )   ) ! Soil   Water        Potential
      allocate            ( Khydsv(kcolp,mwp,-nsoil:0   )   ) ! Soil   Hydraulic    Conductiv.
      allocate            ( Rootsv(kcolp,mwp,-nsoil:0)      ) ! Root Water Pump      [kg/m2/s]
      allocate            ( EExcsv(kcolp,mwp) )               ! Energy in Excess, current


! =================================
! ALLOCATION Mod_SISVAT_loc -   END
! =================================



! =================================
! ALLOCATION Mod_SISVAT_dzS - BEGIN
! =================================


      allocate            ( islpSV(         -nsoil:0) )   !
      allocate            ( isnpSV(   nsnow) )            !
      allocate            ( islmSV(         -nsoil:0) )   !

      allocate            ( dzmiSV(         -nsoil:0) )   ! dz_(i-1/2)
      allocate            ( dzi_SV(         -nsoil:0) )   ! dz_(i-1)/(dz_(i)+dz_(i-1))
      allocate            ( dziiSV(         -nsoil:0) )   ! dz_(i)  /(dz_(i)+dz_(i-1))
      allocate            ( dtz_SV(         -nsoil:0) )   ! dt / dz
      allocate            ( dz78SV(         -nsoil:0) )   ! 7/8 (dz)
      allocate            ( dz34SV(         -nsoil:0) )   ! 3/4 (dz)
      allocate            ( dz_8SV(         -nsoil:0) )   ! 1/8 (dz)
      allocate            ( dzAvSV(         -nsoil:0) )   ! 1/8dz_(-1)+3/4dz+1/8dz_(+1)
      allocate            ( RF__SV( 0:nvege,-nsoil:0) )   ! Root Fraction


! =================================
! ALLOCATION Mod_SISVAT_dzS -   END
! =================================



! =================================
! ALLOCATION Mod_SISVAT_gpt - BEGIN
! =================================

      allocate  ( sst_SB    (kcolp) )                         !  Ocean  FORCING (SST)                         [K]
      allocate  ( sif_SB    (kcolp) )                         !  Ocean  FORCING (Sea-Ice Fraction)            [-]
      allocate  ( MaskSV_gpt(kcolp) )                         !  Land(1)-Sea(0) Mask         (Cell value)     [-]
      allocate  ( Alb_SV_gpt(kcolp) )                         !  Surface  Albedo             (Cell average)   [-]
      allocate  ( EmisSV_gpt(kcolp) )                         !  Surface  LongWave Emissivity(Cell average)   [-]
      allocate  ( Tas_SV_gpt(kcolp) )                         !  Surface  Air    Temperature (Cell average)   [K]
      allocate  ( HSenSV_gpt(kcolp) )                         !  Sensible Heat   Flux   (+ => Upward)      [W/m2]
      allocate  ( HLatSV_gpt(kcolp) )                         !  Latent   Heat   Flux   (+ => Upward)      [W/m2]
      allocate  ( LMO_SV_gpt(kcolp) )                         !  Obukhov  Length             (Cell average)   [m]
      allocate  ( us__SV_gpt(kcolp) )                         !  Friction Velocity                          [m/s]
      allocate  ( uts_SV_gpt(kcolp) )                         !  Sensible Heat   Flux Turbulent Scale     [K m/s]
      allocate  ( uqs_SV_gpt(kcolp) )                         !  Latent   Heat   Flux Turbulent Scale [kg/kg m/s]
      allocate  ( WE2aSV_gpt(kcolp) )                         !  Cumulative H2O  Flux from the Surface  [mm w.e.]
      allocate  ( hFraSV_gpt(kcolp) )                         !  Frazil   Thickness                           [m]
      allocate  ( dpktSV_gpt(kcolp,mzp) )                     !  Reduced  Potential Temperature Tendency   [KX/s]


! =================================
! ALLOCATION Mod_SISVAT_gpt -   END
! =================================



! =================================
! ALLOCATION Mod_SISVAT_flx - BEGIN
! =================================

! OUTPUT for Stand Alone NetCDF File
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      allocate            ( SOsoKL(kcolp,mwp) )               ! Absorbed Solar    Radiation                    [W/m2]
      allocate            ( IRsoKL(kcolp,mwp) )               ! Absorbed IR       Radiation                    [W/m2]
      allocate            ( HSsoKL(kcolp,mwp) )               ! Absorbed Sensible Heat Flux                    [W/m2]
      allocate            ( HLsoKL(kcolp,mwp) )               ! Absorbed Latent   Heat Flux                    [W/m2]
      allocate            ( HLs_KL(kcolp,mwp) )               ! Evaporation                               [mm w.e./s]
      allocate            ( HLv_KL(kcolp,mwp) )               ! Transpiration                                  [W/m2]


! =================================
! ALLOCATION Mod_SISVAT_flx -   END
! =================================



! =================================
! ALLOCATION Mod_SISVAT_cdf - BEGIN
! =================================

! OUTPUT for Stand Alone NetCDF File
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      allocate  ( SOsoNC_xyn(ixp1:mxpp,jyp1:mypp,mwp) )   ! Absorbed Solar    Radiation                    [W/m2]
      allocate  ( IRsoNC_xyn(ixp1:mxpp,jyp1:mypp,mwp) )   ! Absorbed IR       Radiation                    [W/m2]
      allocate  ( HSsoNC_xyn(ixp1:mxpp,jyp1:mypp,mwp) )   ! Absorbed Sensible Heat Flux                    [W/m2]
      allocate  ( HLsoNC_xyn(ixp1:mxpp,jyp1:mypp,mwp) )   ! Absorbed Latent   Heat Flux                    [W/m2]
      allocate  ( HLs_NC_xyn(ixp1:mxpp,jyp1:mypp,mwp) )   ! Evaporation                               [mm w.e./s]
      allocate  ( HLv_NC_xyn(ixp1:mxpp,jyp1:mypp,mwp) )   ! Transpiration                                  [W/m2]
      allocate  ( eta_NC_xyn(ixp1:mxpp,jyp1:mypp,mwp) )   ! Soil Humidity                                 [kg/kg]


! =================================
! ALLOCATION Mod_SISVAT_cdf -   END
! =================================



! =================================
! ALLOCATION Mod_SISVAT_xAW - BEGIN
! =================================

! #AW allocate            ( V__mem(kcolp,mwp,ntave) )         ! V       Time           Steps
! #AW allocate            ( VVmmem(kcolp,mwp      ) )         ! V       Time Moving Averages


! =================================
! ALLOCATION Mod_SISVAT_xAW -   END
! =================================



! =================================
! ALLOCATION Mod_SISVAT_xAH - BEGIN
! =================================

! #AH allocate            ( T__mem(kcolp,mwp,ntave) )         ! dT(a-s) Time           Steps
! #AH allocate            ( dTmmem(kcolp,mwp      ) )         ! dT(a-s) Time Moving Averages


! =================================
! ALLOCATION Mod_SISVAT_xAH -   END
! =================================



! =================================
! ALLOCATION Mod_SISVAT_xZM - BEGIN
! =================================

! #ZM allocate            ( z0_mem(kcolp,mwp,ntavz) )         ! z0      Time           Steps
! #ZM allocate            ( r0_mem(kcolp,mwp,ntavz) )         ! r0      Time           Steps
! #ZM allocate            ( b0_mem(kcolp,mwp,ntavz) )         ! b0      Time           Steps


! =================================
! ALLOCATION Mod_SISVAT_xZM -   END
! =================================



! =================================
! ALLOCATION Mod_SISVAT_xWL - BEGIN
! =================================

! #WL allocate            ( WL_mem(kcolp,mwp,nLimi) )         ! WVLimit Time           Steps
! #WL allocate            ( WLmmem(kcolp,mwp)       )         ! WVLimit Time moving Averages


! =================================
! ALLOCATION Mod_SISVAT_xWL -   END
! =================================



! ============================================================================================
! Initialization of the Correspondance between 1-D Hor.SISVAT Grid and the 3-D Horizontal Grid
! ============================================================================================

      DO i=ixp1,mxpp
      DO j=jyp1,mypp
      DO n=1,mwp

                ikl    = ((n-1)*(mypp-jyp1+1)+j-jyp1) *(mxpp-ixp1+1) + i -ixp1+1
         ii__SV(ikl)   =                                               i
         jj__SV(ikl)   =                      j
         nn__SV(ikl)   =   n
         ikl_SV(i,j,n) =   ikl



! ============================================================================================
! Initialization of the Correspondance between 1-D Hor.SISVAT Grid and the 1-D Hor.AtmPHY Grid
! ============================================================================================

         ikp_SV(ikl)   =   ikl_AP(i,j) 

      ENDDO
      ENDDO
      ENDDO

      return
      end subroutine PHY_SISVAT_ALLOC
