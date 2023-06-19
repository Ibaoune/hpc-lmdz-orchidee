      subroutine SISVAT_ini
 
!--------------------------------------------------------------------------+
!                                                                          |
!   MAR          SISVAT_ini                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_ini generates non time dependant SISVAT parameters |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   dt__SV   : Time  Step                                   [s] |
!     ^^^^^    dz_dSV   : Layer Thickness                              [m] |
!                                                                          |
!     OUTPUT:  RF__SV   : Root Fraction in Layer isl                   [-] |
!     ^^^^^^   rocsSV   : Soil Contrib. to (ro c)_s exclud.Water  [J/kg/K] |
!              etamSV   : Soil Minimum Humidity                    [m3/m3] |
!                        (based on a prescribed Soil Relative Humidity)    |
!              s1__SV   : Factor of eta**( b+2) in Hydraul.Diffusiv.       |
!              s2__SV   : Factor of eta**( b+2) in Hydraul.Conduct.        |
!              aKdtSV   : KHyd: Piecewise Linear Profile:  a * dt    [m]   |
!              bKdtSV   : KHyd: Piecewise Linear Profile:  b * dt    [m/s] |
!              dzsnSV(0): Soil first Layer Thickness                   [m] |
!              dzmiSV   : Distance between two contiguous levels       [m] |
!              dz78SV   : 7/8 (Layer Thickness)                        [m] |
!              dz34SV   : 3/4 (Layer Thickness)                        [m] |
!              dz_8SV   : 1/8 (Layer Thickness)                        [m] |
!              dzAvSV   : 1/8  dz_(i-1) + 3/4 dz_(i) + 1/8 dz_(i+1)    [m] |
!              dtz_SV   : dt/dz                                      [s/m] |
!              OcndSV   : Swab Ocean / Soil Ratio                      [-] |
!              Implic   : Implicit Parameter  (0.5:  Crank-Nicholson)      |
!              Explic   : Explicit Parameter = 1.0 - Implic                |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Possibility                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                          |
!     #SH: Soil /Vegetation Model: Hapex-Sahel   Vegetation     DATA       |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #ER: Richards Equation is not smoothed                               |
!     #kd: Soil: De Ridder Discretization is forced                        |
!                                                                          |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
      use Mod_SISVAT_ctr
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
 
 
      IMPLICIT NONE
 
 
 
! Internal Variables
! ==================
 
      integer           ::  ivt   ,ist   ,ikl,ikv   ,isl   ,isn   ,ikh
      integer           ::  misl_2,nisl_2
      real(kind=real8)  ::  zDepth
      real(kind=real8)  ::  d__eta,eta__1,eta__2,Khyd_1,Khyd_2
      real(kind=real8)  ::  RHsMin=0.001                    ! Min.Soil Relative Humidity
      real(kind=real8)  ::  PsiMax                          ! Max.Soil Water    Potential
      real(kind=real8)  ::  a_Khyd,b_Khyd                   ! Piecewis.Water Conductivity
 
! OUTPUT/Verification: Soil Vertic.Discret.
! #kw real(kind=real8)  ::  Khyd_x,Khyd_y
 
 
 
! Non Time Dependant SISVAT parameters
! ====================================
 
! Decay of Angle(Wind,Sastrugi) Influence on z0 (Andreas, 1995, CCREL report 95-16)
! ---------------------------------------------
 
! #Za   Adz0dt = exp(-dt__SV/43200.)
 
 
! Soil Discretization
! -------------------
 
! Numerical Scheme Parameters
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^
        Implic = 0.75                           ! 0.5  <==> Crank-Nicholson
        Explic = 1.00 - Implic                  !
 
! Soil/Snow Layers Indices
! ^^^^^^^^^^^^^^^^^^^^^^^^
      DO  isl=-nsoil,0
        islpSV(isl) =            isl+1
        islpSV(isl) = min(       islpSV(isl),0)
        islmSV(isl) =            isl-1
        islmSV(isl) = max(-nsoil,islmSV(isl))
      END DO
 
      DO  isn=1,nsnow
        isnpSV(isn) =           isn+1
        isnpSV(isn) = min(      isnpSV(isn),nsnow)
      END DO
 
! Soil      Layers Thicknesses: De Ridder discretization
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^^^^^
! #kd IF (nsoil.gt.4)                                             THEN
! #kd   DO isl=-5,-nsoil,-1
! #kd     dz_dSV(isl)=   1.
! #kd   END DO
! #kd END IF
 
! Soil      Layers Thicknesses: standard  discretization
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^^^^^
      IF (nsoil.ne.4)                                             THEN
        DO isl= 0,-nsoil,-1
          misl_2 =     -mod(isl,2)
          nisl_2 =         -isl/2
          dz_dSV(isl)=(((1-misl_2) * 0.001                              &
     &                  +  misl_2  * 0.003) * 10**(nisl_2)) * 4.
!         dz_dSV(0)  =         Hapex-Sahel Calibration:       4 mm
 
        END DO
! tun     dz_dSV(0)  =               0.001
! tun     dz_dSV(-1) = dz_dSV(-1)  - dz_dSV(0)          + 0.004
      END IF
 
        zz_dSV      = 0.
      DO  isl=-nsoil,0
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
 
! Richards Equation is not smoothed
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #ER   dz78SV(isl) =        dz_dSV(isl)
! #ER   dz34SV(isl) =        dz_dSV(isl)
! #ER   dz_8SV(isl) = 0.
! #ER   dzAvSV(isl) =        dz_dSV(isl)
 
        zz_dSV      = zz_dSV+dz_dSV(isl)
      END DO
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
        dzsnSV(ikl,ikv,0) =      dz_dSV(0)
      END DO
      END DO
 
! Conversion to a 50 m Swab Ocean Discretization
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OcndSV = 0.
      DO isl=-nsoil,0
        OcndSV = OcndSV +dz_dSV(isl)
      END DO
        OcndSV = 50.    /OcndSV
 
 
! Secondary Vegetation Parameters
! -------------------------------
 
! Minimum Stomatal Resistance (Hapex Sahel Data)
! (Taylor et al. 1997, J.Hydrol 188-189, p.1047)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO ivg=1,3                       !
        StodSV(ivg) = 210.             ! Millet
      END DO                           !
        StodSV(  4) = 120.             ! Sparse Tiger Bush
      DO ivg=5,6                       !
        StodSV(ivg) =  80.             ! Dense  Tiger Bush
      END DO                           !
        StodSV(  7) =  80.             ! Low    Trees (Fallow)
        StodSV( 10) =  80.             !
 
! Minimum Stomatal Resistance (Tropical Forest)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        StodSV(  8) =  60.             ! Medium Trees
        StodSV( 11) =  60.             !
        StodSV(  9) =  40.             ! High   Trees
        StodSV( 12) =  40.             !
 
! Root Fraction
! ^^^^^^^^^^^^^
!    * GENERAL REFERENCE
!      Jackson et al., 1996: A global analysis of root distributions for
!      terrestrial biomes. In Oecologia, 108, 389-411.
 
!    * ROOT PROFILE
!      The cumulative root fraction Y is given by
!          Y = 1 - beta**d   with d    the depth (in cm),
!                                 beta a coefficient (vegetation dependent).
 
!    * BETA VALUES (for 11 world biomes)
!    1  boreal forest                0.943
!    2  crops                        0.961
!    3  desert                       0.975
!    4  sclerophyllous shrubs        0.964
!    5  temperate coniferous forest  0.976
!    6  temperate deciduous forest   0.966
!    7  temperate grassland          0.943
!    8  tropical deciduous forest    0.961
!    9  tropical evergreen forest    0.962
!    10 tropical grassland savanna   0.972
!    11 tundra                       0.914
 
!    * ADVISED BETA VALUES FOR MAR
!      (see 'block data SISVAT_dat', variable rbtdSV)
!
!      SVAT veg. type         default      West Africa
!      0  barren soil         0.000        0.000
!      1  crops low           0.961 (2)    0.961 (2)
!      2  crops medium        0.961 (2)    0.961 (2)
!      3  crops high          0.961 (2)    0.961 (2)
!      4  grass low           0.943 (7)    0.943 (7)
!      5  grass medium        0.943 (7)    0.964 (4)
!      6  grass high          0.943 (7)    0.972 (10)
!      7  broadleaf low       0.966 (6)    0.968 (4,10)
!      8  broadleaf medium    0.966 (6)    0.962 (8,9)
!      9  broadleaf high      0.966 (6)    0.962 (8,9)
!      10 needleleaf low      0.976 (5)    0.971 (5,6)
!      11 needleleaf medium   0.976 (5)    0.976 (5)
!      12 needleleaf high     0.976 (5)    0.976 (5)
 
!      Numbers between brackets refer to Jackson's biomes. For more details
!      about some choices, see the correspondance between the IGBP and SVAT
!      vegetation classes (i.e. in NESTOR).
 
!    * WARNING
!      Most of the roots are located in the first 2 m of soil. The root
!      fraction per layer depends on the definition of the soil layer
!      thickness. It will get wrong if a thick layer is defined around 2 m
!      deep.
 
      write(*,'(/a)') 'ROOT PROFILES (Jackson, 1996) :'
 
      DO ivt = 0, nvgt
        zDepth = 0.
        DO isl = 0, -nsoil, -1
          IF (ivt .ne. 0) THEN
            RF__SV(ivt,isl) =  rbtdSV(ivt)**zDepth *                    &
     &                         (1. - rbtdSV(ivt)**(dz_dSV(isl)*100) )
            zDepth = zDepth + dz_dSV(isl)*100  !in cm
          ELSE
            RF__SV(ivt,isl) = 0.
          END IF
        END DO
        write(*,'(a,i2,a,i3,a,99f10.5:)')                               &
     &       '  RF__SV(', ivt, ',',-nsoil, ':0) =', RF__SV(ivt,:)
      END DO
      write(6,6600)
 6600 format(                                                           &
     &  '  NOTE: If root fraction is not close to 0  around 2 m deep,', &
     &/,'        Then you should redefine the soil layer thicknesses.', &
     &/,'        See the code for more details.')
 
 
! Secondary Soil       Parameters
! -------------------------------
 
      DO  ist=0,nsot
         rocsSV(ist)=(1.0-etadSV(ist))*1.2E+6  ! Soil Contrib. to (ro c)_s
         s1__SV(ist)=     bCHdSV(ist)         &! Factor of (eta)**(b+2)
     &  *psidSV(ist)     *Ks_dSV(ist)         &!    in DR97, Eqn.(3.36)
     & /(etadSV(ist)**(   bCHdSV(ist)+3.))     !
         s2__SV(ist)=     Ks_dSV(ist)         &! Factor of (eta)**(2b+3)
     & /(etadSV(ist)**(2.*bCHdSV(ist)+3.))     !    in DR97, Eqn.(3.35)
 
! Soil Minimum Humidity (from a prescribed minimum relative Humidity)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         Psimax = -(log(RHsMin))/7.2E-5        ! DR97, Eqn 3.15 Inversion
         etamSV(ist) =  etadSV(ist)           &!
     &         *(PsiMax/psidSV(ist))**(-min(10.,1./bCHdSV(ist)))
      END DO
         etamSV(12)  =  0.
 
! Piecewise Hydraulic Conductivity Profiles
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO   ist=0,nsot
 
! OUTPUT/Verification: Soil Vertic.Discret.
! #kw     write(6,6000)
 6000     format(' Type |    etaSat | No |    eta__1 |    eta__2 |',    &
     &           '    Khyd_1 |    Khyd_x |    Khyd_2 |    Khyd_y |',    &
     &         /,' -----+-----------+----+-----------+-----------+',    &
     &           '-----------+-----------+-----------+-----------+')
 
          d__eta          =  etadSV(ist)/nkhy
          eta__1          =  0.
          eta__2          =  d__eta
        DO ikh=0,nkhy
          Khyd_1          =  s2__SV(ist)            &! DR97, Eqn.(3.35)
     &  *(eta__1      **(2. *bCHdSV(ist)+3.))        !
          Khyd_2          =  s2__SV(ist)            &!
     &  *(eta__2      **(2. *bCHdSV(ist)+3.))        !
 
          a_Khyd          = (Khyd_2-Khyd_1)/d__eta   !
          b_Khyd          =  Khyd_1-a_Khyd *eta__1   !
 
          aKdtSV(ist,ikh) =  a_Khyd       * dt__SV   !
          bKdtSV(ist,ikh) =  b_Khyd       * dt__SV   !
 
! OUTPUT/Verification: Soil Vertic.Discret.
! #kw     Khyd_x          =  a_Khyd*eta__1 +b_Khyd   !
! #kw     Khyd_y          =  a_Khyd*eta__2 +b_Khyd   !
! #kw     write(6,6001) ist,etadSV(ist),ikh,eta__1, &!
! #kw&          eta__2,Khyd_1,Khyd_x,Khyd_2,Khyd_y   !
 6001     format(i5,' |',e10.2,' |',i3,' |',        &!
     &                 6(e10.2,' |'))
 
          eta__1          = eta__1  + d__eta
          eta__2          = eta__2  + d__eta
        END DO
      END DO
 
 
      return
      end subroutine SISVAT_ini
 
 
 
      subroutine SISVAT(jjtime,kcolw)
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT                               Sat 29-Jun-2013  MAR |
!     SubRoutine SISVAT contains the fortran 77 code of the                |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Thu 14-Feb-2013      |
!           Last Modification by H. Gallee,           Sat 29-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   daHost   : Date Host Model                                  |
!     ^^^^^                                                                |
!                                                                          |
!     INPUT:   LSmask   : 1:          Land       MASK                      |
!     ^^^^^               0:          Sea        MASK                      |
!              ivgtSV   = 0,...,12:   Vegetation Type                      |
!              isotSV   = 0,...,12:   Soil       Type                      |
!                         0:          Water,          Liquid (Sea, Lake)   |
!                        12:          Water, Solid           (Ice)         |
!                                                                          |
!     INPUT:   coszSV   : Cosine of the Sun Zenithal Distance          [-] |
!     ^^^^^    sol_SV   : Surface Downward  Solar      Radiation    [W/m2] |
!              IRd_SV   : Surface Downward  Longwave   Radiation    [W/m2] |
!              drr_SV   : Rain  Intensity                        [kg/m2/s] |
!              dsn_SV   : Snow  Intensity                      [mm w.e./s] |
!              dsnbSV   : Snow  Intensity,  Drift Fraction             [-] |
!              dbs_SV   : Drift Amount                           [mm w.e.] |
!              za__SV   : Surface Boundary Layer (SBL) Height          [m] |
!              VV__SV   :(SBL Top)   Wind Velocity                   [m/s] |
!              TaT_SV   : SBL Top    Temperature                       [K] |
!              rhT_SV   : SBL Top    Air  Density                  [kg/m3] |
!              QaT_SV   : SBL Top    Specific  Humidity            [kg/kg] |
!              qsnoSV   : SBL Mean   Snow      Content             [kg/kg] |
!              LAI0SV   : Leaf Area  Index                             [-] |
!              glf0SV   : Green Leaf Fraction                          [-] |
!              alb0SV   : Soil Basic Albedo                            [-] |
!              dt__SV   : Time  Step                                   [s] |
!                                                                          |
!     INPUT /  isnoSV   = total Nb of Ice/Snow Layers                      |
!     OUTPUT:  ispiSV   = 0,...,nsno: Uppermost Superimposed Ice Layer     |
!     ^^^^^^   iiceSV   = total Nb of Ice      Layers                      |
!              istoSV   = 0,...,5 :   Snow     History (see istdSV data)   |
!                                                                          |
!     INPUT /  alb_SV   : Surface-Canopy Albedo                        [-] |
!     OUTPUT:  emi_SV   : Surface-Canopy Emissivity                    [-] |
!     ^^^^^^   IRs_SV   : Soil           IR Flux  (negative)        [W/m2] |
!              LMO_SV   : Monin-Obukhov               Scale            [m] |
!              us__SV   : Friction          Velocity                 [m/s] |
!              uts_SV   : Temperature       Turbulent Scale          [m/s] |
!              uqs_SV   : Specific Humidity Velocity                 [m/s] |
!              uss_SV   : Blowing Snow      Turbulent Scale          [m/s] |
!              usthSV   : Blowing Snow      Erosion   Threshold      [m/s] |
!              Z0m_SV   : Momentum     Roughness Length                [m] |
!              Z0mmSV   : Momentum     Roughness Length (time mean)    [m] |
!              Z0mnSV   : Momentum     Roughness Length (instantaneous)[m] |
!              Z0SaSV   : Sastrugi     Roughness Length                [m] |
!              Z0e_SV   : Erosion Snow Roughness Length                [m] |
!              Z0emSV   : Erosion Snow Roughness Length (time mean)    [m] |
!              Z0enSV   : Erosion Snow Roughness Length (instantaneous)[m] |
!              Z0roSV   : Subgrid Topo Roughness Length                [m] |
!              Z0h_SV   : Heat         Roughness Length                [m] |
!              snCaSV   : Canopy   Snow     Thickness            [mm w.e.] |
!              rrCaSV   : Canopy   Water    Content                [kg/m2] |
!              psivSV   : Leaf     Water    Potential                  [m] |
!              TvegSV   : Canopy   Temperature                         [K] |
!              TsisSV   : Soil/Ice Temperatures (layers -nsoil ,...,0)     |
!                       & Snow     Temperatures (layers  1,2,..,nsnow) [K] |
!              ro__SV   : Soil/Snow Volumic Mass                   [kg/m3] |
!              eta_SV   : Soil/Snow Water   Content                [m3/m3] |
!              G1snSV   : snow dendricity/sphericity                       |
!              G2snSV   : snow sphericity/grain size                       |
!              dzsnSV   : Snow Layer        Thickness                  [m] |
!              agsnSV   : Snow       Age                             [day] |
!              BufsSV   : Snow Buffer Layer              [kg/m2] .OR. [mm] |
!              BrosSV   : Snow Buffer Layer Density      [kg/m3]           |
!              BG1sSV   : Snow Buffer Layer Dendricity / Sphericity    [-] |
!              BG2sSV   : Snow Buffer Layer Sphericity / Size [-] [0.1 mm] |
!              rusnSV   : Surficial   Water              [kg/m2] .OR. [mm] |
!                                                                          |
!     OUTPUT:  no__SV   : OUTPUT file Unit Number                      [-] |
!     ^^^^^^   i___SV   : OUTPUT point   i Coordinate                  [-] |
!              j___SV   : OUTPUT point   j Coordinate                  [-] |
!              n___SV   : OUTPUT point   n Coordinate                  [-] |
!              lwriSV   : OUTPUT point vec Index                       [-] |
!                                                                          |
!     OUTPUT:  IRu_SV   : Upward     IR Flux (+, upw., effective)      [K] |
!     ^^^^^^   hSalSV   : Saltating Layer Height                       [m] |
!              qSalSV   : Saltating Snow  Concentration            [kg/kg] |
!              RnofSV   : RunOFF Intensity                       [kg/m2/s] |
!                                                                          |
!     Internal Variables:                                                  |
!     ^^^^^^^^^^^^^^^^^^                                                   |
!              NLaysv   = New            Snow Layer Switch             [-] |
!              albisv   : Snow/Ice/Water/Soil Integrated Albedo        [-] |
!              SoCasv   : Absorbed Solar Radiation by Canopy (Normaliz)[-] |
!              SoSosv   : Absorbed Solar Radiation by Surfac.(Normaliz)[-] |
!              tau_sv   : Fraction of Radiation transmitted by Canopy  [-] |
!              TBr_sv   : Brightness Temperature                       [K] |
!              IRupsv   : Upward     IR Flux (-, upw.)              [W/m2] |
!              IRv_sv   : Vegetation IR Flux                        [W/m2] |
!              rrMxsv   : Canopy Maximum Intercepted Rain          [kg/m2] |
!              Sigmsv   : Canopy Ventilation Factor                    [-] |
!              ram_sv   : Aerodynamic Resistance for Momentum        [s/m] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!              HSv_sv   : Vegetation Sensible Heat Flux             [W/m2] |
!              HLv_sv   : Vegetation Latent   Heat Flux             [W/m2] |
!              Rootsv   : Root Water Pump                        [kg/m2/s] |
!              Evp_sv   : Evaporation                              [kg/m2] |
!              EvT_sv   : Evapotranspiration                       [kg/m2] |
!              HSs_sv   : Surface    Sensible Heat Flux + => absorb.[W/m2] |
!              HLs_sv   : Surface    Latent   Heat Flux + => absorb.[W/m2] |
!              Lx_H2O   : Latent Heat of Vaporization/Sublimation   [J/kg] |
!              Tsrfsv   : Surface    Temperature                       [K] |
!              LAI_sv   : Leaf Area  Index (snow included)             [-] |
!              LAIesv   : Leaf Area  Index (effective / transpiration) [-] |
!              glf_sv   : Green Leaf Fraction of NOT fallen Leaves     [-] |
!              sEX_sv   : Verticaly Integrated Extinction Coefficient  [-] |
!              LSdzsv   : Vertical   Discretization Factor             [-] |
!                       =    1. Soil                                       |
!                       = 1000. Ocean                                      |
!              z_snsv   : Snow Pack  Thickness                         [m] |
!              zzsnsv   : Snow Pack  Thickness                         [m] |
!              albssv   : Soil       Albedo                            [-] |
!              Evg_sv   : Soil+Vegetation Emissivity                   [-] |
!              Eso_sv   : Soil+Snow       Emissivity                   [-] |
!              psi_sv   : Soil       Water    Potential                [m] |
!              Khydsv   : Soil   Hydraulic    Conductivity           [m/s] |
!                                                                          |
!              ETVg_d   : VegetationEnergy Power         Forcing    [W/m2] |
!              ETSo_0   : Snow/Soil Energy Power, before Forcing    [W/m2] |
!              ETSo_1   : Snow/Soil Energy Power, after  Forcing    [W/m2] |
!              ETSo_d   : Snow/Soil Energy Power         Forcing    [W/m2] |
!              EqSn_0   : Snow      Energy, before Phase Change     [J/m2] |
!              EqSn_1   : Snow      Energy, after  Phase Change     [J/m2] |
!              EqSn_d   : Snow      Energy,       net    Forcing    [J/m2] |
!              Enrsvd   : SVAT      Energy Power         Forcing    [W/m2] |
!              Enrbal   : SVAT      Energy Balance                  [W/m2] |
!              Wats_0   : Soil Water,  before Forcing                 [mm] |
!              Wats_1   : Soil Water,  after  Forcing                 [mm] |
!              Wats_d   : Soil Water          Forcing                 [mm] |
!              SIWm_0   : Snow        initial Mass               [mm w.e.] |
!              SIWm_1   : Snow        final   Mass               [mm w.e.] |
!              SIWa_i   : Snow Atmos. initial Forcing            [mm w.e.] |
!              SIWa_f   : Snow Atmos. final   Forcing(noConsumed)[mm w.e.] |
!              SIWe_i   : SnowErosion initial Forcing            [mm w.e.] |
!              SIWe_f   : SnowErosion final   Forcing(noConsumed)[mm w.e.] |
!              SIsubl   : Snow sublimed/deposed  Mass            [mm w.e.] |
!              SImelt   : Snow Melted            Mass            [mm w.e.] |
!              SIrnof   : Surficial Water + Run OFF Change       [mm w.e.] |
!              SIvAcr   : Sea-Ice    vertical Acretion           [mm w.e.] |
!              Watsvd   : SVAT Water          Forcing                 [mm] |
!              Watbal   : SVAT Water  Balance                       [W/m2] |
!                                                                          |
!              dsn_Ca,snCa_n :     Snow Contribution to the Canopy[m w.e.] |
!              drr_Ca,rrCa_n,drip: Rain Contribution to the Canopy [kg/m2] |
!              vk2      : Square of Von Karman Constant                [-] |
!              sqrCm0   : Factor of   Neutral Drag Coeffic.Momentum  [s/m] |
!              sqrCh0   : Factor of   Neutral Drag Coeffic.Heat      [s/m] |
!              EmiVeg   : Vegetation    Emissivity                     [-] |
!              EmiSol   : Soil          Emissivity                     [-] |
!              EmiSno   : Snow          Emissivity                     [-] |
!              EmiWat   : Water         Emissivity                     [-] |
!              Z0mSea   :          Sea  Roughness Length               [m] |
!              Z0mLnd   :          Land Roughness Length               [m] |
!              sqrrZ0   : u*t/u*                                           |
!              f_eff    : Marticorena & B. 1995 JGR (20)                   |
!              A_Fact   : Fundamental * Roughness                          |
!              Z0mBSn   :         BSnow Roughness Length               [m] |
!              Z0mBS0   : Mimimum BSnow Roughness Length (blown* )     [m] |
!              Z0m_Sn   :          Snow Roughness Length (surface)     [m] |
!              Z0m_S0   : Mimimum  Snow Roughness Length               [m] |
!              Z0m_S1   : Maximum  Snow Roughness Length               [m] |
!              Z0_GIM   : Minimum GIMEX Roughness Length               [m] |
!              Z0_ICE   : Sea Ice ISW   Roughness Length               [m] |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Possibility                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                          |
!     #AE: TURBULENCE: Aerosols Erosion / Turbulent Diffusion Coeff.       |
!     #BD: TraCer   Aeolian Erosion  Submodel           is turned ON       |
!     #BS: Explicit Cloud MICROPHYSICS: Blow. *(Snow)         Model)       |
!     #SN: SNOW Model                               may be turned ON       |
!     #NP: SNOW Model: Snow Properties may be those of Polar Snow          |
!     #ZG: SNOW Model: ETH-Camp & Greenland 3D simulations                 |
!     #MB: SNOW Model: Erosion Efficiency (Marticorena & Berga.1995)       |
!     #SI: SISVAT: Sea-Ice Fraction calculated from prescribed SST         |
!     #MT: SISVAT: Monin-Obukhov Theory is linearized (Garrat schem)       |
!     #SH: Soil /Vegetation Model: Hapex-Sahel   Vegetation     DATA       |
!     #ZO: SBL: Orography Roughness included from SL_z0 in MARdom          |
!     #ZW: SBL: Mom.: Roughn.Length= F(u*) Wang  MWR 129     ,  Sea        |
!     #ZT: SBL: Mom.: Roughn.Length= Typical value in polar models         |
!     #ZS: SBL: Mom.: Roughn.Length= F(u*) Andreas     (1995)  Snow        |
!     #Za: SBL: Mom.: Roughn.Length= F(u*) Andreas &al.(2004), Snow (modif.|
!     #ZA: SBL: Mom.: Roughn.Length= F(u*) Andreas &al.(2004), Snow (native|
!     #RS: SBL: Heat: Roughn.Length= F(u*,z0)  Andreas (1987)  Snow, Ice   |
!     #ZM: SBL: M/H   Roughn.Length: Box Moving Average (in Time)          |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Col de Porte                         |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^                                      |
!     #CP: Col de Porte  Turbulence        Parameterization                |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #Zw: SBL: Mom.: Roughn.Length= F(u*) Wang  MWR 129 bis ,  Sea        |
!     #ZN: SBL: Mom.: Roughn.Length= F(u*) Shao  & Lin (1999), Snow        |
!     #ZL: SBL: Z0mL  Roughn.Length= F(glf)                                |
!     #FL: SISVAT: LAI Assignation and Fallen Leaves Correction            |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_iii_jjj_n     | #e0: OUTPUT on ASCII  File (SISVAT Variables) |
!   #                      |      Energy Budg. Verif.: Soil+(Sea-Ice)+Snow |
!   #                      |(#e0  MUST BE PREPROCESSED BEFORE #e1 & #e2 !) |
!   # SISVAT_iii_jjj_n     | #m0: OUTPUT/Verification: H2O    Conservation |
!                          |                                               |
!   # stdout               | #s0: OUTPUT of Snow Buffer Layer              |
!                          |      unit  6, SubRoutine  SISVAT     **ONLY** |
!   # stdout               | #s2: OUTPUT of SnowFall, Snow Buffer          |
!                          |      unit  6, SubRoutine  SISVAT_BSn, _qSn    |
!   # stdout               | #b0: OUTPUT of Snow Erosion Statistics        |
!                          |      unit  6, SubRoutine  SISVAT_BSn **ONLY** |
!   # stdout               | #sf: OUTPUT of SnowFall, Z0 and Drag Coeff.   |
!                          |      unit  6, SubRoutines PHY_SISVAT, SISVAT  |
!   # stdout               | #sz: OUTPUT of Roughness Length & Drag Coeff. |
!                          |      unit  6, SubRoutine  SISVAT     **ONLY** |
!                                                                          |
!     SUGGESTIONS of MODIFICATIONS: see lines beginning with "C +!!!"      |
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                         |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_ctr
      use Mod_SISVAT_CdP
      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
      use Mod_SISVAT_log
      use Mod_SISVAT_loc
      use Mod_SISVAT_aux
! #BS use Mod_SISVAT_BSn
 
      use Mod_SISVATLmmm
 
 
 
      IMPLICIT NONE
 
 
 
      integer           ::  jjtime
      integer           ::  kcolw
 
 
 
! Internal Variables
! ==================
 
! Non Local
! ---------
 
      real(kind=real8)                                 d_Bufs,Bufs_N   ! Buffer Snow Layer Increment
      real(kind=real8)                                 Buf_ro,Bros_N   ! Buffer Snow Layer Density
! #NP real(kind=real8)                                 BufPro          ! Buffer Snow Layer Density
      real(kind=real8)                                 Buf_G1,BG1__N   ! Buffer Snow Layer Dendr/Sphe[-]
      real(kind=real8)                                 Buf_G2,BG2__N   ! Buffer Snow Layer Spher/Size[-]
 
! Energy         Budget
! ~~~~~~~~~~~~~~~~~~~~~
! #e1 real(kind=real8),              dimension(kcolw)::ETSo_0          ! Soil/Snow Power, before Forcing
! #e1 real(kind=real8),              dimension(kcolw)::ETSo_1          ! Soil/Snow Power, after  Forcing
! #e1 real(kind=real8),              dimension(kcolw)::ETSo_d          ! Soil/Snow Power, Forcing
 
 
! Local
! -----
 
! #e0 character(len=16)                            ::  FilNam          !
! #e0 integer                                      ::  noUNIT          ! OUTPUT File  Unit Number
 
      integer                                      ::  ikl   ,ikv      !
      integer                                      ::  isn   ,isl      !
      integer                                      ::  ist   ,n        !
      integer                                      ::  ist__s,ist__w   ! Soil/Water Body Identifier
      integer                                      ::  growth          ! Seasonal               Mask
      integer                                      ::  LISmsk          ! Land+Ice / Open    Sea Mask
      integer                                      ::  LSnMsk          ! Snow-Ice / No Snow-Ice Mask
      integer                                      ::  IceMsk          !      Ice               Mask
      integer                                      ::  SnoMsk          ! Snow     / No Snow     Mask
 
      real(kind=real8)                             ::  drip            ! Rain Contribution to the Canopy
      real(kind=real8)                             ::  drr_Ca,rrCa_n   ! Rain Contribution to the Canopy
      real(kind=real8)                             ::  dsn_Ca,snCa_n   ! Snow Contribution to the Canopy
      real(kind=real8)                             ::  roSMin =  30.   ! Minimum Snow Density
      real(kind=real8)                             ::  roSn_1 = 109.   ! Fallen  Snow Density, Indep. Param. (PAHAUT)
      real(kind=real8)                             ::  roSn_2 =   6.   ! Fallen  Snow Density, Temper.Param. (PAHAUT)
      real(kind=real8)                             ::  roSn_3 =  26.   ! Fallen  Snow Density, Wind   Param. (PAHAUT)
      real(kind=real8)                             ::  Dendr1 =  17.12 ! Fallen  Snow Dendric. Wind 1/Param. (GIRAUD)
      real(kind=real8)                             ::  Dendr2 = 128.   ! Fallen  Snow Dendric. Wind 2/Param. (GIRAUD)
      real(kind=real8)                             ::  Dendr3 = -20.   ! Fallen  Snow Dendric. Indep. Param. (GIRAUD)
 
      real(kind=real8)                             ::  Spher1 =   7.87 ! Fallen  Snow Spheric.,Wind 1/Param.
      real(kind=real8)                             ::  Spher2 =  38.   ! Fallen  Snow Spheric.,Wind 2/Param.
      real(kind=real8)                             ::  Spher3 =  50.   ! Fallen  Snow Spheric.,Wind 3/Param.
      real(kind=real8)                             ::  Spher4 =  90.   ! Fallen  Snow Spheric.,Indep. Param. (GIRAUD)
      real(kind=real8)                             ::  Polair          ! Polar   Snow Switch
! #BS real(kind=real8)                             ::  PorSno          !
      real(kind=real8)                             ::  Salt_f,PorRef   !
 
! For Diffusion of Surficial Water in the Snow Pack
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #DW real(kind=real8)                             ::  PorVol,rWater   ! Pore Volume, retained Water
! #DW real(kind=real8)                             ::  rusNEW,rdzNEW,etaNEW
 
! #BS real(kind=real8)                             ::  ro_new          !
! #NP real(kind=real8)                             ::  TaPole = 263.15 ! Maximum     Polar Temperature
      real(kind=real8)                             ::  T__Min = 200.00 ! Minimum realistic Temperature
 
! DATA      Emissivities                                               ! Pielke, 1984, pp. 383,409
! ^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real8)                             ::  EmiVeg =   0.98 ! Emissivity of Vegetation
      real(kind=real8)                             ::  EmiSol =   0.94 ! Emissivity of       Soil
      real(kind=real8)                             ::  EmiSno =   0.99 ! Emissivity of            Snow
      real(kind=real8)                             ::  EmiWat =   0.99 ! Emissivity of a Water Area
 
      real(kind=real8)                             ::  epsLMO = 1.e-18 ! minimum absolute value of LMo
      real(kind=real8)                             ::  vk2             ! Square of Von Karman Constant
      real(kind=real8)                             ::  u2star          !(u*)**2
! #ZL real(kind=real8)                             ::  fallen =   0.00 ! Fallen   Leaves         Switch
      real(kind=real8)                             ::  Z0mSea,Z0hSea   !          Sea  Roughness Length
      real(kind=real8)                             ::  Z0mLnd          !          Land Roughness Length
! #ZN real(kind=real8)                             ::  sqrrZ0          ! u*t/u*
! #MB real(kind=real8)                             ::  f_eff           ! Marticorena & B. 1995 JGR (20)
      real(kind=real8)                             ::  A_Fact          ! Fundamental * Roughness
 
      real(kind=real8)                             ::  Z0m_nu          ! Smooth R Snow Roughness Length
      real(kind=real8)                             ::  Z0mBSn          !         BSnow Roughness Length
      real(kind=real8)                             ::  Z0mBS0 = 0.5e-6 ! MINimum BSnow Roughness Length, Momentum
                                                                       ! Gallee et al. 2001 BLM 99 (19)
!     real(kind=real8)                             ::  Z0m_S0 = 0.00005! MINimum  Snow Roughness Length
!     real(kind=real8)                             ::  Z0m_S1 = 0.030  ! MAXimum  Snow Roughness Length, Sastrugis
! #ZS real(kind=real8)                             ::  Z0Sa_N          ! Regime   Snow Roughness Length
! #ZS real(kind=real8)                             ::  Z0SaSi          ! 1.IF Rgm Snow Roughness Length
! #ZG real(kind=real8)                             ::  Z0_GIM = 0.0013 ! Mimimum GIMEX Roughness Length  Ice Min Z0 = 0.0013 m (Broeke)
                                                                       !                                 Old Ice Z0 = 0.0500 m (Bruce)
                                                                       !                                              0.0500 m (Smeets)
                                                                       !                                              0.1200 m (Broeke)
      real(kind=real8)                             ::  Z0_ICE = 0.0010 ! Sea-Ice ISW   Roughness Length                        (Andreas)
      real(kind=real8)                             ::  Z0m_Sn          ! Snow  Surface Roughness Length
! #Za real(kind=real8)                             ::  Z0m_90          ! Snow  Surface Roughness Length
      real(kind=real8)                             ::  SnoWat          ! Snow Layer    Switch
! #RS real(kind=real8)                             ::  rstar,alors     !
! #RS real(kind=real8)                             ::  rstar0,rstar1   !
! #RS real(kind=real8)                             ::  rstar2          !
      real(kind=real8)                             ::  SameOK          ! 1. => Same Type of Grains
      real(kind=real8)                             ::  G1same          ! Averaged G1,  same Grains
      real(kind=real8)                             ::  G2same          ! Averaged G2,  same Grains
      real(kind=real8)                             ::  typ__1          ! 1. => Lay1 Type: Dendritic
      real(kind=real8)                             ::  zroNEW          ! dz X ro, if fresh Snow
      real(kind=real8)                             ::  G1_NEW          ! G1,      if fresh Snow
      real(kind=real8)                             ::  G2_NEW          ! G2,      if fresh Snow
      real(kind=real8)                             ::  zroOLD          ! dz X ro, if old   Snow
      real(kind=real8)                             ::  G1_OLD          ! G1,      if old   Snow
      real(kind=real8)                             ::  G2_OLD          ! G2,      if old   Snow
      real(kind=real8)                             ::  SizNEW          ! Size,    if fresh Snow
      real(kind=real8)                             ::  SphNEW          ! Spheric.,if fresh Snow
      real(kind=real8)                             ::  SizOLD          ! Size,    if old   Snow
      real(kind=real8)                             ::  SphOLD          ! Spheric.,if old   Snow
      real(kind=real8)                             ::  Siz_av          ! Averaged    Grain Size
      real(kind=real8)                             ::  Sph_av          ! Averaged    Grain Spher.
      real(kind=real8)                             ::  Den_av          ! Averaged    Grain Dendr.
      real(kind=real8)                             ::  DendOK          ! 1. => Average is  Dendr.
      real(kind=real8)                             ::  G1diff          ! Averaged G1, diff. Grains
      real(kind=real8)                             ::  G2diff          ! Averaged G2, diff. Grains
      real(kind=real8)                             ::  G1              ! Averaged G1
      real(kind=real8)                             ::  G2              ! Averaged G2
 
! Energy and Mass Budget
! ~~~~~~~~~~~~~~~~~~~~~~
! #e1 real(kind=real8)                             ::    EnsBal        ! Soil+Snow   , Power  Balance
! #e1 real(kind=real8)                             ::    EnvBal        !      Vegetat, Power  Balance
 
                                                                       ! H2O    Conservation
! #m0 real(kind=real8)                                   Watbal        ! Soil+Vegetat, Water  Balance
 
                                                                       ! * Mass Conservation
! #m1 real(kind=real8)                             ::    SnoBal        ! Snow Pack     Mass   Balance
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
! ALLOCATION                                              !
! ==========                                              !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)               THEN !
 
      allocate            ( TBr_sv(kcolp,mwp) )           ! Brightness Temperature
      allocate            ( IRdwsv(kcolp,mwp) )           ! DOWNward   IR Flux
      allocate            ( IRupsv(kcolp,mwp) )           ! UPward     IR Flux
      allocate            ( Bdzssv(kcolp,mwp) )           ! Buffer Snow Layer Thickness
      allocate            ( z_snsv(kcolp,mwp) )           ! Snow-Ice, current Thickness
 
! Energy         Budget
! ~~~~~~~~~~~~~~~~~~~~~
! #e1 allocate            ( ETVg_d(kcolp,mwp) )           ! VegetationPower, Forcing
! #e1 allocate            ( EqSn_0(kcolp,mwp) )           ! Snow Energy, befor Phase Change
! #e1 allocate            ( EqSn_1(kcolp,mwp) )           ! Snow Energy, after Phase Change
! #e1 allocate            ( EqSn_d(kcolp,mwp) )           ! Energy in Excess
 
! OUTPUT/Verification: H2O    Conservation
! #m0 allocate            ( Wats_0(kcolp,mwp) )           ! Soil Water,  before Forcing
! #m0 allocate            ( Wats_1(kcolp,mwp) )           ! Soil Water,  after  Forcing
! #m0 allocate            ( Wats_d(kcolp,mwp) )           ! Soil Water,         Forcing
 
! OUTPUT/Verification: * Mass Conservation
! #m1 allocate            ( SIsubl(kcolp,mwp) )           ! Snow Sublimed/Deposed Mass
! #m1 allocate            ( SImelt(kcolp,mwp) )           ! Snow Melted           Mass
! #m1 allocate            ( SIrnof(kcolp,mwp) )           ! Local Surficial Water + Run OFF
 
! OUTPUT/Verification: SeaIce Conservation
! #m2 allocate            ( SIvAcr(kcolp,mwp) )           ! Sea-Ice      Vertical Acretion
 
! Energy and Mass Budget
! ~~~~~~~~~~~~~~~~~~~~~~
! #e1 allocate            ( Enrsvd(kcolp,mwp) )           ! Soil+Vegetat  Power  Forcing
 
                                                          ! H2O    Conservation
! #m0 allocate            ( Watsv0(kcolp,mwp) )           ! Soil+Vegetat, before Forcing
! #m0 allocate            ( Watsvd(kcolp,mwp) )           ! Soil+Vegetat  Water  Forcing
 
                                                          ! * Mass Conservation
! #m1 allocate            ( SIWm_0(kcolp,mwp) )           ! Snow Initial              Mass
! #m1 allocate            ( SIWm_1(kcolp,mwp) )           ! Snow Final                Mass
! #m1 allocate            ( SIWa_i(kcolp,mwp) )           ! Snow Initial       ATM Forcing
! #m1 allocate            ( SIWa_f(kcolp,mwp) )           ! Snow Final         ATM Forcing
! #m1 allocate            ( SIWe_i(kcolp,mwp) )           ! Snow Initial       BLS Forcing
! #m1 allocate            ( SIWe_f(kcolp,mwp) )           ! Snow Final         BLS Forcing
 
 
      allocate            ( IcIndx(kcolp,mwp) )           ! No   Ice               Mask
 
      allocate            ( FallOK(kcolp,mwp) )           ! Snow Contribution to the Canopy
 
      END IF                                              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! Internal DATA
! =============
 
      vk2    =  vonKrm  *  vonKrm             ! Square of Von Karman Constant
! #FL fallen =             1.                 ! Fallen  Leaves         Switch
 
! #ZD Z0m_S0 =             0.00200            ! MINimum Snow Roughness Length
                                              ! MegaDunes    included
 
 
! BEGIN.main.
! ===========
 
      IF (.not.iniOUT)                                              THEN
               iniOUT = .true.
 
! Snow Pack Thickness
! -------------------
 
          DO ikl=1,kcolp
          DO ikv=1,mwp
            z_snsv(ikl,ikv)     = 0.0
          END DO
          END DO
        DO   isn=1,nsnow
          DO ikl=1,kcolp
          DO ikv=1,mwp
            z_snsv(ikl,ikv)     = z_snsv(ikl,ikv) + dzsnSV(ikl,ikv,isn)
            zzsnsv(ikl,ikv,isn) = z_snsv(ikl,ikv)
          END DO
          END DO
        END DO
 
 
! SISVAT Forcing VERIFICATION
! ---------------------------
 
        IF (IRs_SV(1,1).gt.-eps6)                                       &
     &  write(6,600)
 600    format(/,'### SISVAT ERROR, Soil IR Upward  not defined ###',   &
     &         /,'###               Initialize and Store IRs_SV ###')
 
 
! OUTPUT
! ------
 
                    FilLab              ='SISVAT'
                    SepLab              ='_'
                    nwUNIT              = 51
      END IF
 
! #e0 DO ikl=1,kcolp
! #e0 DO ikv=1,mwp
! #e0   IF   (lwriSV(ikl,ikv).ne.0.AND.no__SV(lwriSV(ikl,ikv)).eq.0)THEN
! #e0                nwUNIT              = nwUNIT+1
! #e0                no__SV(lwriSV(ikl,ikv)) = nwUNIT
! #e0      write(FilNam,'(a6,a1,2(i3.3,a1),i1)')                        &
! #e0&           FilLab,SepLab,i___SV(lwriSV(ikl,ikv)),                 &
! #e0&                  SepLab,j___SV(lwriSV(ikl,ikv)),                 &
! #e0&                  SepLab,n___SV(lwriSV(ikl,ikv))
! #e0      open(unit=nwUNIT,status='unknown',file=FilNam)
! #e0      rewind    nwUNIT
! #e0   END IF
! #e0 END DO
! #e0 END DO
 
! #e0 DO ikl=1,kcolp
! #e0 DO ikv=1,mwp
! #e0   IF (lwriSV(ikl,ikv).ne.0)                                   THEN
! #e0           noUNIT=no__SV(lwriSV(ikl,ikv))
! #e0     write(noUNIT,5000) daHost,i___SV(lwriSV(ikl,ikv)),            &
! #e0&                              j___SV(lwriSV(ikl,ikv)),            &
! #e0&                              n___SV(lwriSV(ikl,ikv)),            &
! #e0&                                     Z0m_SV(ikl,ikv) ,            &
! #e0&                                     albisv(ikl,ikv)
 5000     format(                                                       &
     &       /,              a18,'|           Grid Point ',2i4,         &
     &                                           ' (',i2,')',           &
     &         '    | Z0m =',f12.6,' | Albedo = ',f6.3,' |',            &
     &       /,' -------+',7('---------+'),2('--------+'))
! #e0   END IF
! #e0 END DO
! #e0 END DO
 
 
! "Soil" Humidity of Water Bodies
! ===============================
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
          ist    =      isotSV(ikl,ikv)                   ! Soil Type
          ist__s =  min(ist, 1)                           ! 1 => Soil
          ist__w =  1 - ist__s                            ! 1 => Water Body
        DO isl=-nsoil,0
          eta_SV(ikl,ikv,isl) = eta_SV(ikl,ikv,isl) * ist__s &! Soil
     &                    + etadSV(ist)     * ist__w      ! Water Body
        END DO
 
 
! Vertical Discretization Factor
! ==============================
 
          LSdzsv(ikl,ikv)     =                   ist__s &! Soil
     &                    + OcndSV          * ist__w      ! Water Body
      END DO
      END DO
 
 
! Vegetation Temperature Limits
! =============================
 
          DO ikl=1,kcolp
          DO ikv=1,mwp
            TvegSV(ikl,ikv) = max(TvegSV(ikl,ikv),T__Min) ! T__Min = 200.K
 
 
! LAI Assignation and Fallen Leaves Correction (#FL)! Dead Leaves fall
! ==================================================! => LAI = Green Leaves only
                                                    ! => GLF = 1
 
            LAI0SV(ikl,ikv) =     LAI0SV(ikl,ikv)*min(1,ivgtSV(ikl,ikv)) ! NO LAI if
!                                                            ! no vegetation
            glf_sv(ikl,ikv) =     glf0SV(ikl,ikv)
! #FL       glf_sv(ikl,ikv) =     1.                         ! #FL
            LAI_sv(ikl,ikv) =     LAI0SV(ikl,ikv)        &   !
! #FL&               *        glf0SV(ikl,ikv)            &   ! #FL
     &               +        0.
          END DO
          END DO
 
 
! LAI in Presence of Snow
! =======================
 
!         ASSUMPTION: LAI decreases   when Snow Thickness increases,
!         ^^^^^^^^^^      becoming  0 when Snow Thickn. = Displac.Height
          DO ikl=1,kcolp
          DO ikv=1,mwp
            LAI_sv(ikl,ikv) =     LAI_sv(ikl,ikv)                       &
     &               * (1.0 - zzsnsv(       ikl,ikv, isnoSV(ikl,ikv))   &
     &                      /(DH_dSV(ivgtSV(ikl,ikv))+eps6)      )
            LAI_sv(ikl,ikv) = max(LAI_sv(ikl,ikv),zer0)
            LAI_sv(ikl,ikv) = min(LAI_sv(ikl,ikv),ea_Max)
          END DO
          END DO
 
 
! Interception of Rain by the Canopy
! ==================================
 
! OUTPUT/Verification: H2O    Conservation: Vegetation Forcing
! #m0     DO ikl=1,kcolp
! #m0     DO ikv=1,mwp
! #m0       Watsv0(ikl,ikv) =      rrCaSV(ikl,ikv)   ! Canopy Water Cont.
! #m0       Watsvd(ikl,ikv) =      drr_SV(ikl,ikv)   ! Precipitation
! #m0     END DO
! #m0     END DO
 
 
! New Canopy Water Content
! ------------------------
 
          DO ikl=1,kcolp
          DO ikv=1,mwp
            rrMxsv(ikl,ikv) = 0.2*max(     eps6,LAI_sv(ikl,ikv)        )  ! Precip. Max. Intercept.
            Sigmsv(ikl,ikv) = 1.0-exp(-min(half*LAI_sv(ikl,ikv),ea_MAX))  ! Canopy Ventilation Coe.
!                                                                 ! (DR97, eqn 3.6)
            drr_Ca      = drr_SV(ikl,ikv)  *Sigmsv(ikl,ikv)      &! Intercepted Rain
     &                                 *dt__SV                    !
            rrCa_n      = rrCaSV(ikl,ikv)  +drr_Ca                ! New Canopy Water Contnt
                                                                  ! (DR97, eqn 3.28)
            drip        = rrCa_n       -rrMxsv(ikl,ikv)           ! Water  Drip
            drip        =      max(zer0,drip)                     !
            rrCa_n      = rrCa_n       -drip                      !
            IF (rrCa_n.LT.1.e-30)       rrCa_n = 0.               !
            drr_SV(ikl,ikv) = drr_SV(ikl,ikv) +(rrCaSV(ikl,ikv)  &! Update Rain  Contribut.
     &                                 -rrCa_n     )             &!
     &                                 /dt__SV                    !
            rrCaSV(ikl,ikv) = rrCa_n                              ! Upd.Canopy Water Contnt
 
 
! Interception of Snow by the Canopy
! ==================================
 
            dsn_Ca      = dsn_SV(ikl,ikv)  *Sigmsv(ikl,ikv) &! Intercepted Snow
     &                                 *dt__SV       !
            snCa_n      = snCaSV(ikl,ikv)  +dsn_Ca   ! New Canopy Snow Thickn.
            drip        = snCa_n       -rrMxsv(ikl,ikv)  !
            drip        =      max(zer0,drip)        !
            snCa_n      = snCa_n       -drip         !
            dsn_SV(ikl,ikv) = dsn_SV(ikl,ikv) +(snCaSV(ikl,ikv) &! Update Snow  Contribut.
     &                                 -snCa_n     )&!
     &                                 /dt__SV       !
            snCaSV(ikl,ikv) = snCa_n                 ! Upd.Canopy Snow Thickn.
          END DO
          END DO
 
 
! Snow Fall from the Canopy
! =========================
 
!         ASSUMPTION: snow fall from the canopy,
!         ^^^^^^^^^^  when the temperature of the vegetation is positive
!               (.OR. when snow over the canopy is saturated  with water)
 
          DO ikl=1,kcolp
          DO ikv=1,mwp
            FallOK(ikl,ikv) =  max(zer0,sign(un_1,TvegSV(ikl,ikv)-Tf_Sno+eps6)) &
     &                  *  max(zer0,sign(un_1,snCaSV(ikl,ikv)       -eps6))
            dsn_SV(ikl,ikv) =      dsn_SV(ikl,ikv)   +snCaSV(ikl,ikv)*FallOK(ikl,ikv)   &
     &                                       /dt__SV
            snCaSV(ikl,ikv) =      snCaSV(ikl,ikv) * (1.         -FallOK(ikl,ikv))
 
 
! Blowing Particles Threshold Friction velocity
! =============================================
 
! #AE       usthSV(ikl,ikv) =                     1.0e+2
          END DO
          END DO
 
 
! Contribution of Snow to the Surface Snow Pack
! =============================================
 
      IF (SnoMod)                                                 THEN
 
 
! OUTPUT/Verification: * Mass Conservation
! #m1   DO ikl=1,kcolp
! #m1   DO ikv=1,mwp
! #m1     SIWa_i(ikl,ikv) =(drr_SV(ikl,ikv) + dsn_SV(ikl,ikv))   *dt__SV ![mm w.e.]
! #m1     SIWe_i(ikl,ikv) = dbs_SV(ikl,ikv)                          !
! #m1     SIWm_0(ikl,ikv) = BufsSV(ikl,ikv) + HFraSV(ikl,ikv)    *rhoIce !
! #m1   DO isn=1,nsnow                                               !
! #m1     SIWm_0(ikl,ikv) = SIWm_0(ikl,ikv) + dzsnSV(ikl,ikv,isn)*ro__SV(ikl,ikv,isn)!
! #m1   END DO                                                       !
! #m1   END DO
! #m1   END DO                                                       !
 
 
! Blowing Snow
! ------------
 
!                         **********
        IF (BloMod)  call SISVAT_BSn
!                         **********
 
!                         **********
! #ve                call SISVAT_wEq('_BSn  ',1)
!                         **********
 
 
! Sea Ice
! -------
 
!            **********
! #SI   call SISVAT_SIc(                                                &
! #m2&                  SIvAcr                                          &
! #SI&                 )
!            **********
 
!            **********
! #ve   call SISVAT_wEq('_SIc  ',0)
!            **********
 
 
! Buffer Layer
! ------------
 
          DO ikl=1,kcolp
          DO ikv=1,mwp
            BufsSV(ikl,ikv) =      BufsSV(ikl,ikv)      !     [mm w.e.]
            d_Bufs      =  max(dsn_SV(ikl,ikv) *dt__SV,0.)  ! i.e., [kg/m2]
            dsn_SV(ikl,ikv) =      0.                   !
            Bufs_N      =      BufsSV(ikl,ikv) +d_Bufs  !
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND.&
! #s0&          ikv        .EQ.nwr_SV)                                 &
! #s0&      write(6,6601)      BufsSV(ikl,ikv) ,d_Bufs,Bufs_N
 6601       format(/,'Buffer *: ',3e15.6)
 
! Snow Density
! ^^^^^^^^^^^^
            Polair      =      0.00                     !
! #NP       Polair      =  max(zer0,                   &!
! #NP&                         sign(un_1,TaPole        &!
! #NP&                                  -TaT_SV(ikl,ikv)))  !
            Buf_ro      =  max( rosMin,                &! Fallen Snow Density
     &      roSn_1+roSn_2*     (TaT_SV(ikl,ikv)-Tf_Sno)&! [kg/m3]
     &            +roSn_3*sqrt( VV10SV(ikl,ikv)))       ! Pahaut    (CEN)
! #NP       BufPro      =  max( rosMin,                &! Fallen Snow Density
! #NP&         104. *sqrt( max( VV10SV(ikl,ikv)-6.0,0.0)))  ! Kotlyakov (1961)
            Bros_N      = (1. - Polair) *   Buf_ro     &! Temperate Snow
! #NP&                        + Polair  *   BufPro     &! Polar     Snow
     &                  +  0.                           !
 
! Instantaneous Density of deposited blown Snow (de Montmollin, 1978)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #BS       PorSno =      1.0d00     -  BSnoRo                          &
! #BS&                               /  rhoIce
! #BS       Salt_f =      usthSV(ikl,ikv)/  max(eps6,   us__SV(ikl,ikv))
! #BS       Salt_f =  min(Salt_f     ,  un_1)
! #BS       PorRef =      PorSno     /  max(eps6,1.-PorSno)             &
! #BS&               +log(Salt_f)
! #BS       Por_BS =      PorRef     /          (1.+PorRef)
! #BS       ro_new =      rhoIce     *          (1.-Por_BS)
! #BS       ro_new =  max(ro_new     ,  BSnoRo)
! #BS       Bros_N      = Bros_N     * (1.0-dsnbSV(ikl,ikv))            &
! #BS&                  + ro_new     *      dsnbSV(ikl,ikv)
 
! Instantaneous Density IF deposited blown Snow (Melted* from Canopy)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Bros_N      = Bros_N     * (1.0-FallOK(ikl,ikv))            &!
     &                  + 300.       *      FallOK(ikl,ikv)              !
 
! Time averaged Density of deposited blown Snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            BrosSV(ikl,ikv) =(Bros_N     *      d_Bufs                  &!
     &                   +BrosSV(ikl,ikv)*      BufsSV(ikl,ikv))        &!
     &                   /         max(eps6,Bufs_N)                      !
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND. &
! #s0&          ikv        .EQ.nwr_SV)                                  &
! #s0&      write(6,6602) Buf_ro,Bros_N,BrosSV(ikl,ikv),dsnbSV(ikl,ikv)
 6602       format('rho    *: ',3e15.6,'    dsnbSV: ',e15.6)
 
!  S.Falling Snow Properties (computed as in SISVAT_zAg)
!    ^^^^^^^^^^^^^^^^^^^^^^^
            Buf_G1      =  max(-G1_dSV,                                 &! Temperate Snow
     &               min(Dendr1*VV__SV(ikl,ikv)-Dendr2,                 &!     Dendricity
     &                   Dendr3                   ))                     !
            Buf_G2      =  min( Spher4,                                 &! Temperate Snow
     &               max(Spher1*VV__SV(ikl,ikv)+Spher2,                 &!     Sphericity
     &                   Spher3                   ))                     !
            Buf_G1      = (1. - Polair) *   Buf_G1                      &! Temperate Snow
     &                        + Polair  *   G1_dSV                       ! Polar     Snow
            Buf_G2      = (1. - Polair) *   Buf_G2                      &! Temperate Snow
     &                        + Polair  *   ADSdSV                       ! Polar     Snow
                G1      =                   Buf_G1                       ! NO  Blown Snow
                G2      =                   Buf_G2                       ! NO  Blown Snow
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND. &
! #s0&          ikv        .EQ.nwr_SV)                                  &
! #s0&      write(6,6603)       BG1sSV(ikl,ikv),BG2sSV(ikl,ikv)
 6603       format('G1,G2  *: ',3e15.6)
 
! S.1. Meme  Type  de Neige  / same Grain Type
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #BS       SameOK  =  max(zer0,                                        &
! #BS&                     sign(un_1,    Buf_G1             *G1_dSV     &
! #BS&                                 - eps_21                    ))
! #BS       G1same  = ((1.0-dsnbSV(ikl,ikv))*Buf_G1+dsnbSV(ikl,ikv) *G1_dSV)
! #BS       G2same  = ((1.0-dsnbSV(ikl,ikv))*Buf_G2+dsnbSV(ikl,ikv) *ADSdSV)
!           Blowing Snow Properties:                         G1_dSV, ADSdSV
 
! S.2. Types differents / differents Types
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #BS       typ__1  =  max(zer0,sign(un_1,eps6-Buf_G1))   ! =1.=> Dendritic
! #BS       zroNEW  =     typ__1  *(1.0-dsnbSV(ikl,ikv)) &! fract.Dendr.Lay.
! #BS&              + (1.-typ__1) *     dsnbSV(ikl,ikv)   !
! #BS       G1_NEW  =     typ__1  *Buf_G1                &! G1 of Dendr.Lay.
! #BS&              + (1.-typ__1) *G1_dSV                 !
! #BS       G2_NEW  =     typ__1  *Buf_G2                &! G2 of Dendr.Lay.
! #BS&              + (1.-typ__1) *ADSdSV                 !
! #BS       zroOLD  = (1.-typ__1) *(1.0-dsnbSV(ikl,ikv)) &! fract.Spher.Lay.
! #BS&              +     typ__1  *     dsnbSV(ikl,ikv)   !
! #BS       G1_OLD  = (1.-typ__1) *Buf_G1                &! G1 of Spher.Lay.
! #BS&              +     typ__1  *G1_dSV                 !
! #BS       G2_OLD  = (1.-typ__1) *Buf_G2                &! G2 of Spher.Lay.
! #BS&              +     typ__1  *ADSdSV                 !
! #BS       SizNEW  =    -G1_NEW  *DDcdSV/G1_dSV         &! Size  Dendr.Lay.
! #BS&               +(1.+G1_NEW         /G1_dSV)        &!
! #BS&                  *(G2_NEW  *DScdSV/G1_dSV         &!
! #BS&               +(1.-G2_NEW         /G1_dSV)*DFcdSV) !
! #BS       SphNEW  =     G2_NEW         /G1_dSV          ! Spher.Dendr.Lay.
! #BS       SizOLD  =     G2_OLD                          ! Size  Spher.Lay.
! #BS       SphOLD  =     G1_OLD         /G1_dSV          ! Spher.Spher.Lay.
! #BS       Siz_av =     (zroNEW*SizNEW+zroOLD*SizOLD)    ! Averaged Size
! #BS       Sph_av = min( zroNEW*SphNEW+zroOLD*SphOLD    &!
! #BS&                   ,  un_1)                         ! Averaged Sphericity
! #BS       Den_av = min((Siz_av -(    Sph_av *DScdSV    &!
! #BS&                            +(1.-Sph_av)*DFcdSV))  &!
! #BS&                 / (DDcdSV -(    Sph_av *DScdSV    &!
! #BS&                            +(1.-Sph_av)*DFcdSV))  &!
! #BS&                   ,  un_1)                         !
! #BS       DendOK  = max(zer0,                           !
! #BS&                    sign(un_1,     Sph_av *DScdSV  &! Small   Grains
! #BS&                              +(1.-Sph_av)*DFcdSV  &! Faceted Grains
! #BS&                              -    Siz_av        )) !
!           REMARQUE: le  type moyen (dendritique ou non) depend
!           ^^^^^^^^  de la  comparaison avec le diametre optique
!                     d'une neige recente de   dendricite nulle
!           REMARK:   the mean type  (dendritic   or not) depends
!           ^^^^^^    on the comparaison with the optical diameter
!                     of a recent snow    having zero dendricity
 
! #BS       G1diff  =(   -DendOK *Den_av                 &!
! #BS&               +(1.-DendOK)*Sph_av) *G1_dSV         !
! #BS       G2diff  =     DendOK *Sph_av  *G1_dSV        &!
! #BS&               +(1.-DendOK)*Siz_av                  !
! #BS       G1      =     SameOK *G1same                 &!
! #BS&               +(1.-SameOK)*G1diff                  !
! #BS       G2      =     SameOK *G2same                 &!
! #BS&               +(1.-SameOK)*G2diff                  !
 
            BG1__N      =((1. - FallOK(ikl,ikv))*   G1   &!
     &                        + FallOK(ikl,ikv) *   99.) &! Melted *  from Canopy
     &                  *       d_Bufs/max(eps6,d_Bufs)   !
            BG2__N      =((1. - FallOK(ikl,ikv))*   G2   &!
     &                        + FallOK(ikl,ikv) *   30.) &! Melted *  from Canopy
     &                  *       d_Bufs/max(eps6,d_Bufs)   !
 
!  S.Buffer  Snow Properties (computed as in SISVAT_zAg)
!    ^^^^^^^^^^^^^^^^^^^^^^^
            Buf_G1      =       BG1__N                    ! Falling   Snow
            Buf_G2      =       BG2__N                    ! Falling   Snow
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND. &
! #s0&          ikv        .EQ.nwr_SV)                                  &
! #s0&      write(6,6604)      Buf_G1      ,Buf_G2         ,FallOK(ikl,ikv)  &
! #s0&                                                     ,TvegSV(ikl,ikv)
 6604       format('G1,G2 F*: ',3e15.6,'    T__Veg: ',e15.6)
 
! S.1. Meme  Type  de Neige  / same Grain Type
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            SameOK  =  max(zer0,                                        &
     &                     sign(un_1,    Buf_G1 *BG1sSV(ikl,ikv)        &
     &                                 - eps_21                    ))
            G1same  = (d_Bufs*Buf_G1+BufsSV(ikl,ikv)*BG1sSV(ikl,ikv))   &
     &                     /max(eps6,Bufs_N)
            G2same  = (d_Bufs*Buf_G2+BufsSV(ikl,ikv)*BG2sSV(ikl,ikv))   &
     &                     /max(eps6,Bufs_N)
 
! S.2. Types differents / differents Types
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            typ__1  =  max(zer0,sign(un_1,eps6-Buf_G1))   ! =1.=> Dendritic
            zroNEW  =(    typ__1  *d_Bufs                &! fract.Dendr.Lay.
     &              + (1.-typ__1) *BufsSV(ikl,ikv))      &!
     &                   /max(eps6,Bufs_N)                !
            G1_NEW  =     typ__1  *Buf_G1                &! G1 of Dendr.Lay.
     &              + (1.-typ__1) *BG1sSV(ikl,ikv)        !
            G2_NEW  =     typ__1  *Buf_G2                &! G2 of Dendr.Lay.
     &              + (1.-typ__1) *BG2sSV(ikl,ikv)        !
            zroOLD  =((1.-typ__1) *d_Bufs                &! fract.Spher.Lay.
     &              +     typ__1  *BufsSV(ikl,ikv))      &!
     &                   /max(eps6,Bufs_N)                !
            G1_OLD  = (1.-typ__1) *Buf_G1                &! G1 of Spher.Lay.
     &              +     typ__1  *BG1sSV(ikl,ikv)        !
            G2_OLD  = (1.-typ__1) *Buf_G2                &! G2 of Spher.Lay.
     &              +     typ__1  *BG2sSV(ikl,ikv)        !
            SizNEW  =    -G1_NEW  *DDcdSV/G1_dSV         &! Size  Dendr.Lay.
     &               +(1.+G1_NEW         /G1_dSV)        &!
     &                  *(G2_NEW  *DScdSV/G1_dSV         &!
     &               +(1.-G2_NEW         /G1_dSV)*DFcdSV) !
            SphNEW  =     G2_NEW         /G1_dSV          ! Spher.Dendr.Lay.
            SizOLD  =     G2_OLD                          ! Size  Spher.Lay.
            SphOLD  =     G1_OLD         /G1_dSV          ! Spher.Spher.Lay.
            Siz_av  =   ( zroNEW  *SizNEW+zroOLD*SizOLD)  ! Averaged Size
            Sph_av = min( zroNEW  *SphNEW+zroOLD*SphOLD  &!
     &                  ,   un_1                       )  ! Averaged Sphericity
            Den_av = min((Siz_av  - (    Sph_av *DScdSV  &!
     &                              +(1.-Sph_av)*DFcdSV))&!
     &                 / (DDcdSV  - (    Sph_av *DScdSV  &!
     &                              +(1.-Sph_av)*DFcdSV))&!
     &                  ,   un_1                         )!
            DendOK  = max(zer0,                          &!
     &                    sign(un_1,     Sph_av *DScdSV  &! Small   Grains
     &                              +(1.-Sph_av)*DFcdSV  &! Faceted Grains
     &                              -    Siz_av        )) !
!           REMARQUE: le  type moyen (dendritique ou non) depend
!           ^^^^^^^^  de la  comparaison avec le diametre optique
!                     d'une neige recente de   dendricite nulle
!           REMARK:   the mean type  (dendritic   or not) depends
!           ^^^^^^    on the comparaison with the optical diameter
!                     of a recent snow    having zero dendricity
 
            G1diff  =(   -DendOK *Den_av                                &
     &               +(1.-DendOK)*Sph_av) *G1_dSV
            G2diff  =     DendOK *Sph_av  *G1_dSV                       &
     &               +(1.-DendOK)*Siz_av
            G1      =     SameOK *G1same                                &
     &               +(1.-SameOK)*G1diff
            G2      =     SameOK *G2same                                &
     &               +(1.-SameOK)*G2diff
 
            BG1sSV(ikl,ikv) =                       G1                  &!
     &                  *       Bufs_N/max(eps6,Bufs_N)                  !
            BG2sSV(ikl,ikv) =                       G2                  &!
     &                  *       Bufs_N/max(eps6,Bufs_N)                  !
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND. &
! #s0&          ikv        .EQ.nwr_SV)                                  &
! #s0&      write(6,6605) Buf_G1     ,typ__1                            &
! #s0&                   ,DendOK     ,Den_av     ,Sph_av     ,Siz_av    &
! #s0&                   ,G1same     ,G1diff     ,G1
 6605       format('B1,Typ  : ',2e15.6,11x,'OK,Den,Sph,Siz: ',4e15.6    &
     &          ,/,'          ',30x   ,11x,'sam,dif,G1    : ',3e15.6)
 
! Update of Buffer Layer Content & Decision about creating a new snow layer
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            BufsSV(ikl,ikv) =       Bufs_N              !     [mm w.e.]
            NLaysv(ikl,ikv) = min(un_1,                &!
     &                    max(zer0,                    &! Allows to create
     &                        sign(un_1,BufsSV(ikl,ikv)&! a new snow Layer
     &                                 -SMndSV     ))  &! if Buffer > SMndSV
     &                   *max(zer0,                    &! Except if * Erosion
     &                        sign(un_1,half           &! dominates
     &                                 -dsnbSV(ikl,ikv)))  &!
     &                   +max(zer0,                    &! Allows to create
     &                        sign(un_1,BufsSV(ikl,ikv)&! a new snow Layer
     &                                 -SMndSV*3.00)))  ! is Buffer > SMndSV*3
 
            Bdzssv(ikl,ikv) = 1.e-3*BufsSV(ikl,ikv)*rhoWat &! [mm w.e.] -> [m w.e.]
     &                            /max(eps6,BrosSV(ikl,ikv))!& [m w.e.] -> [m]
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Buffer G1, G2 variables
! #s0       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND. &
! #s0&          ikv        .EQ.nwr_SV)                                  &
! #s0&      write(6,6606) BG1sSV(ikl,ikv),BG2sSV(ikl,ikv)               &
! #s0&                   ,NLaysv(ikl,ikv),BdzsSV(ikl,ikv)
 6606       format('G1,G2 N*: ',2e15.6,i15,e27.6)
 
          END DO
          END DO
 
 
! Snow Pack Discretization
! ========================
 
!            **********
        call SISVAT_zSn
!            **********
 
!            **********
! #ve   call SISVAT_wEq('_zSn  ',0)
!            **********
 
! OUTPUT in SISVAT for ikl = 1 (preferably for Stand Alone Version)
! OUTPUT           for SnowFall and Snow Buffer
! #s2   IF          (isnoSV(1,1) .GT. 0)                                &
! #s2&  write(6,6004)isnoSV(1,1),  dsn_SV(1) *dt__SV + BufsSV(1),       &
! #s2&              (dzsnSV(1,isn)*ro__SV(1,isn),isn=1,isnoSV(1,1))
 6004   format(i3,'  dsn+Buf=',f6.2,6x,'z dz *ro =',10f6.2,             &
     &                                       (/,35x,10f6.2))
 
 
! Add a new Snow Layer
! ====================
 
          DO ikl=1,kcolp
          DO ikv=1,mwp
            isnoSV(ikl,ikv)     = isnoSV(ikl,ikv)              +NLaysv(ikl,ikv)
            isn                 = isnoSV(ikl,ikv)
            dzsnSV(ikl,ikv,isn) = dzsnSV(ikl,ikv,isn) * float(1-NLaysv(ikl,ikv))   &
     &                          + Bdzssv(ikl,ikv)     * float(  NLaysv(ikl,ikv))
            TsisSV(ikl,ikv,isn) = TsisSV(ikl,ikv,isn) * float(1-NLaysv(ikl,ikv))   &
     &                      + min(TaT_SV(ikl,ikv),Tf_Sno)*float(NLaysv(ikl,ikv))
            ro__SV(ikl,ikv,isn) = ro__SV(ikl,ikv,isn) * float(1-NLaysv(ikl,ikv))   &
     &                          + Brossv(ikl,ikv)     * float(  NLaysv(ikl,ikv))
            eta_SV(ikl,ikv,isn) = eta_SV(ikl,ikv,isn) * float(1-NLaysv(ikl,ikv))! + 0.
            agsnSV(ikl,ikv,isn) = agsnSV(ikl,ikv,isn) * float(1-NLaysv(ikl,ikv))! + 0.
            G1snSV(ikl,ikv,isn) = G1snSV(ikl,ikv,isn) * float(1-NLaysv(ikl,ikv))   &
     &                          + BG1ssv(ikl,ikv)     * float(  NLaysv(ikl,ikv))
            G2snSV(ikl,ikv,isn) = G2snSV(ikl,ikv,isn) * float(1-NLaysv(ikl,ikv))   &
     &                          + BG2ssv(ikl,ikv)     *    NLaysv(ikl,ikv)
            istoSV(ikl,ikv,isn) = istoSV(ikl,ikv,isn) * (1-NLaysv(ikl,ikv))        &
     &       + max(zer0,sign(un_1,TaT_SV(ikl,ikv)                                  &
     &                           -Tf_Sno-eps_21)) *         istdSV(2)              &
     &                                            *         NLaysv(ikl,ikv)
            BufsSV(ikl,ikv)     = BufsSV(ikl,ikv) * float(1-NLaysv(ikl,ikv))
            NLaysv(ikl,ikv)     = 0
          END DO
          END DO
 
 
! Snow Pack Thickness
! -------------------
 
          DO ikl=1,kcolp
          DO ikv=1,mwp
            z_snsv(ikl,ikv)     = 0.0
          END DO
          END DO
        DO   isn=1,nsnow
          DO ikl=1,kcolp
          DO ikv=1,mwp
            z_snsv(ikl,ikv)     = z_snsv(ikl,ikv) + dzsnSV(ikl,ikv,isn)
            zzsnsv(ikl,ikv,isn) = z_snsv(ikl,ikv)
          END DO
          END DO
        END DO
 
 
! Diffusion of Surficial Water in the Snow Pack
! ---------------------------------------------
 
! #DW     DO isn=1,nsnow
! #DW     DO ikl=1,kcolp
! #DW     DO ikv=1,mwp
! #DW       PorVol      = 1.     - ro__SV(ikl,ikv,isn) /  rhoIce         !
! #DW       PorVol      =      max(PorVol            ,zer0  )            !
! #DW       rWater      = ws0dSV * PorVol     *rhoWat*dzsnSV(ikl,ikv,isn)   &
! #DW&                  * max(zer0,                                     &
! #DW&                   sign(un_1,rusnSV(ikl,ikv)/rhoWat-zzsnsv(ikl,ikv,isn)   &
! #DW&                                               +dzsnSV(ikl,ikv,isn)))
! #DW       rusNEW      =      max(rusnSV(ikl,ikv)-rWater,zer0  )
! #DW       rWater      =          rusnSV(ikl,ikv)-rusNEW
! #DW       rdzNEW          =      rWater                               &
! #DW&                           + ro__SV(ikl,ikv,isn) * dzsnSV(ikl,ikv,isn)
! #DW       etaNEW          =      rWater / max(eps6,rdzNEW)
! #DW       rusnSV(ikl,ikv) =          rusNEW
! #DW       ro__SV(ikl,ikv,isn) =      rdzNEW / max(eps6,dzsnSV(ikl,ikv,isn))
! #DW       eta_SV(ikl,ikv,isn) =      eta_SV(ikl,ikv,isn)  +etaNEW
! #DW     END DO
! #DW     END DO
! #DW     END DO
 
      END IF
 
! OUTPUT in SISVAT for ikl = 1 (preferably for Stand Alone Version)
! OUTPUT           for SnowFall and Snow Buffer
! #s2   IF          (isnoSV(1,1) .GT. 0)                                  &
! #s2&  write(6,6006)isnoSV(1,1),  dsn_SV(1) *dt__SV + BufsSV(1),       &
! #s2&              (dzsnSV(1,isn)*ro__SV(1,isn),isn=1,isnoSV(1,1))
 6006   format(i3,'  dsn+Buf=',f6.2,6x,'* dz *ro =',10f6.2,             &
     &                                       (/,35x,10f6.2))
 
 
! Blowing Dust
! ============
 
! #BD   IF (BloMod)                                                 THEN
 
!         ***************
! #BD     call SISVAT_BDu
!         ***************
 
! #BD   END IF
 
 
 
! Soil      Albedo: Soil Humidity Correction
! ==========================================
 
!         REFERENCE: McCumber and Pielke (1981), Pielke (1984)
!         ^^^^^^^^^
          DO ikl=1,kcolp
          DO ikv=1,mwp
            albssv(ikl,ikv) =                                           &
     &      alb0SV(ikl,ikv) *(1.0-min(half,eta_SV(       ikl,ikv,0)     &
     &                                /etadSV(isotSV(ikl,ikv))))
!         REMARK:    Albedo of Water Surfaces (isotSV=0):
!         ^^^^^^     alb0SV := 2  X  effective value, while
!                    eta_SV :=          etadSV
          END DO
          END DO
 
 
! Snow Pack Optical Properties
! ============================
 
      IF (SnoMod)                                                   THEN
 
!            ******
        call SnOptP(                                                   &
! #AG&              jjtime                                             &
     &             )
!            ******
 
      ELSE
        DO ikl=1,kcolp
        DO ikv=1,mwp
          sEX_sv(ikl,ikv,1) = 1.0
          sEX_sv(ikl,ikv,0) = 0.0
          albisv(ikl,ikv)   = albssv(ikl,ikv)
        END DO
        END DO
      END IF
 
!            **********
! #ve   call SISVAT_wEq('SnOptP',0)
!            **********
 
 
! Solar Radiation Absorption and Effective Leaf Area Index
! ========================================================
 
!            ******
        call VgOptP
!            ******
 
 
! Surface-Canopy Emissivity
! =========================
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
            LSnMsk         =     min( 1,isnoSV(ikl,ikv))
            tau_sv(ikl,ikv)=     exp(  -LAI_sv(ikl,ikv))                ! Veg Transmit.Frac.
            Evg_sv(ikl,ikv)=  EmiVeg*(1-LSnMsk)+EmiSno*LSnMsk           ! Veg+Sno Emissivity
            Eso_sv(ikl,ikv)=  EmiSol*(1-LSnMsk)+EmiSno*LSnMsk           ! Sol+Sno Emissivity
            emi_SV(ikl,ikv)=                                           &
     &   (((EmiSol*     tau_sv(ikl,ikv)                                &
     &     +EmiVeg*(1.0-tau_sv(ikl,ikv))) *LSmask(ikl,ikv))            &
     &    + EmiWat                     *(1-LSmask(ikl,ikv)))*(1-LSnMsk)&
     &   +  EmiSno                                             *LSnMsk
        END DO
        END DO
 
 
! Soil/Vegetation Forcing/ Upward IR (INPUT, from previous time step)
! ===================================================================
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
! #e1     Enrsvd(ikl,ikv) =    - IRs_SV(ikl,ikv)
          IRupsv(ikl,ikv) =      IRs_SV(ikl,ikv) *     tau_sv(ikl,ikv) ! Upward   IR
        END DO
        END DO
 
 
! Turbulence
! ==========
 
! Latent Heat of Vaporization/Sublimation
! ---------------------------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          SnoWat      =                     min(isnoSV(ikl,ikv),0)
          Lx_H2O(ikl,ikv) =                                             &
     &    (1.-SnoWat) * LhvH2O                                          &
     &  +     SnoWat  *(LhsH2O * (1.-eta_SV(ikl,ikv,isnoSV(ikl,ikv)))   &
     &                 +LhvH2O *     eta_SV(ikl,ikv,isnoSV(ikl,ikv)) )
        END DO
        END DO
 
 
! Roughness Length for Momentum
! -----------------------------
 
! Land+Sea-Ice / Ice-free Sea Mask
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,kcolp
        DO ikv=1,mwp
          IcIndx(ikl,ikv) = 0
        END DO
        END DO
        DO isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
          IcIndx(ikl,ikv) = max(IcIndx(ikl,ikv),                        &
     &                      isn*max(0,                                  &
     &                              sign(1,                             &
     &                                   int(ro__SV(ikl,ikv,isn)-900.))))
        END DO
        END DO
        END DO
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          LISmsk    =     min(iiceSV(ikl,ikv)  ,1     )
          LISmsk    =     max(LSmask(ikl,ikv),LISmsk)
          IceMsk    =     max(0,sign(1       ,IcIndx(ikl,ikv)-1)  )
          SnoMsk    = max(min(isnoSV(ikl,ikv)-iiceSV(ikl,ikv),1),0)
 
! Sea  Roughness Length
! ^^^^^^^^^^^^^^^^^^^^^
          Z0mSea =       0.0002
          Z0hSea =       0.000049
 
! #Zw     Z0mSea =       0.0185*us__SV(ikl,ikv)*us__SV(ikl,ikv)   ! Doyle MWR 130
! #Zw&                         *Grav_I                   &! p.3088 2e col
 
! #ZW     Z0mSea =       0.016 *us__SV(ikl,ikv)*us__SV(ikl,ikv)  &! Wang  MWR 129
! #ZW&                         *Grav_I                   &! p.1377 (21)
! #ZW&           +       0.11  *A_MolV                   &!
! #ZW&                         /  max(eps6 ,us__SV(ikl,ikv))  !
 
! #Zw     Z0mSea =       0.0185*us__SV(ikl,ikv)*us__SV(ikl,ikv)  &! Wang  MWR 129
! #Zw&                         *Grav_I                   &! p.1377 (21)
! #Zw&           +       0.135 *A_MolV                   &!   (adapted)
! #Zw&                         /  max(eps6 ,us__SV(ikl,ikv))  !
 
! #ZW     Z0hSea =   max(0.000049,                       &! Wang  MWR 129
! #ZW&                   0.20  *A_MolV                   &! p.1377 (22)
! #ZW&                         /  max(eps6 ,us__SV(ikl,ikv)))
 
! #ZW     Z0mSea =   max(Z0mSea,eps6)                     !
 
! Land Roughness Length, Snow Contribution excluded
! ^^^^^^^^^^^^^^^^^^^^^^ Ice  Contribution included
!                        ^^^^^^^^^^^^^^^^^^^^^^^^^^
! If vegetation Seasonal Cycle described by  LAI     :
          growth      =min(max(0,7-ivgtSV(ikl,ikv)),1)
          Z0mLnd      =     Z0mdSV(ivgtSV(ikl,ikv))*LAI_sv(ikl,ikv)*growth  &
     &                                         /LAIdSV                  &
     &                +     Z0mdSV(ivgtSV(ikl,ikv))*         (1-growth)
 
! If vegetation Seasonal Cycle described by  GLF only:
! #ZL     Z0mLnd      =                                                 &
! #ZL&             fallen * Z0mLnd                                      &
! #ZL&        +(1.-fallen)* Z0mdSV(ivgtSV(ikl,ikv))*glf_sv(ikl,ikv)*growth  &
! #ZL&                 +    Z0mdSV(ivgtSV(ikl,ikv))*         (1-growth)
 
! Land Roughness Length, Influence of the Masking by Snow
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0mLnd      =max( Z0mLnd   ,                                  &
     &                      Z0mdSV(0)*(1-IceMsk)                        &
     &                     +Z0_ICE   *   IceMsk )
          Z0mLnd      =     Z0mLnd                                      &
     &                    -(zzsnsv(ikl,ikv,    isnoSV(ikl,ikv))         &
     &                     -zzsnsv(ikl,ikv,max(IcIndx(ikl,ikv),0)))/7.
          Z0mLnd      =max( Z0mLnd    ,    5.e-5  )  ! Min set := Z0 on *
!         Roughness  disappears under Snow
!         Assumption Height/Roughness Length =  7 is used
 
! Z0 Smooth Regime over Snow (Andreas 1995, CRREL Report 95-16, p. 8)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0m_nu =       5.e-5 ! z0s~(10-d)*exp(-vonKrm/sqrt(1.1e-03))
 
! Z0 Saltat.Regime over Snow (Gallee  et al., 2001, BLM 99 (19) p.11)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          u2star =       us__SV(ikl,ikv) *us__SV(ikl,ikv)
          Z0mBSn =       u2star      *0.536e-3   -  61.8e-6
          Z0mBSn =   max(Z0mBS0      ,Z0mBSn)
 
! Z0 Smooth + Saltat. Regime
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0enSV(ikl,ikv) =  Z0m_nu                                     &
     &                +  Z0mBSn
 
! Rough   Snow Surface Roughness Length (Typical Value)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0m_Sn =    0.250e-3 ! Andreas 1995, CRREL Report 95-16, fig.1&p.2
                               ! z0r~(10-d)*exp(-vonKrm/sqrt(1.5e-03))-5.e-5
          Z0m_Sn =    2.000e-3 ! Calibration    of MAR
! #ZT     Z0m_Sn =    1.000e-3 ! Exemple Tuning in RACMO
! #ZT     Z0m_Sn =    0.500e-3 ! Exemple Tuning in MAR
 
! Rough   Snow Surface Roughness Length (Variable Sastrugi Height)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          A_Fact      =  1.0000        ! Andreas et al., 2004, p.4
                                       ! ams.confex.com/ams/pdfpapers/68601.pdf
 
!                                                                        ! 0050=.003/.6
! #ZS     Z0Sa_N =                   (us__SV(ikl,ikv) -0.2)*0.0050      &! 0053=TUNING
! #ZS&           * max(zer0,sign(un_1,Tf_Sno-eps9                       &!
! #ZS&                               -TsisSV(ikl,ikv , isnoSV(ikl,ikv))))
!!#ZS     Z0SaSi = max(zer0,sign(un_1,Z0Sa_N                  ))         ! 1 if erosion
! #ZS     Z0SaSi = max(zer0,sign(un_1,zer0  -eps9 -uss_SV(ikl,ikv)))     !
! #ZS     Z0Sa_N = max(zer0,          Z0Sa_N)
! #ZS     Z0SaSV(ikl,ikv) =                                             &!
! #ZS&             max(Z0SaSV(ikl,ikv)   ,Z0SaSV(ikl,ikv)               &!
! #ZS&               + Z0SaSi*(Z0Sa_N-Z0SaSV(ikl,ikv))*exp(-dt__SV/43200.)) &!
! #ZS&               -            min(dz0_SV(ikl,ikv) ,     Z0SaSV(ikl,ikv)) !
 
! #ZS     A_Fact      =               Z0SaSV(ikl,ikv) *  5.0/0.15        ! A=5 if h~10cm
!         CAUTION: The influence of the sastrugi direction is not yet included
 
! #ZS     Z0m_Sn =                    Z0SaSV(ikl,ikv)                   &!
! #ZS&                              - Z0m_nu                             !
 
! Z0                         (Shao & Lin, 1999, BLM 91 (46)  p.222)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
! Z0 Saltat.Regime over Snow (Shao & Lin, 1999, BLM 91 (46)  p.222)
! #ZN     sqrrZ0 =       usthSV(ikl,ikv)/max( us__SV(ikl,ikv),0.001)
! #ZN     sqrrZ0 =                   min( sqrrZ0     ,0.999)
! #ZN     Z0mBSn =       0.55 *0.55 *exp(-sqrrZ0     *sqrrZ0)           &!
! #ZN&                  *us__SV(ikl,ikv)*     us__SV(ikl,ikv)*Grav_I*0.5 !
 
! Z0 Smooth + Saltat. Regime (Shao & Lin, 1999, BLM 91 (46)  p.222)
! #ZN     Z0enSV(ikl,ikv) = (Z0m_nu     **    sqrrZ0 )                  &!
! #ZN&                * (Z0mBSn     **(1.-sqrrZ0))
! #ZN     Z0enSV(ikl,ikv) =  max(Z0enSV(ikl,ikv), Z0m_nu)
 
! Z0                         (Andreas etAl., 2004
! ^^^^^^^^^^^^^^^^^^^^^^^^^^  ams.confex.com/ams/pdfpapers/68601.pdf)
! Z0 Smooth Regime over Snow (Andreas etAl., 2004
! #Za     Z0m_nu = 0.135*A_MolV / max(us__SV(ikl,ikv) , eps6)
 
! Z0 Saltat.Regime over Snow (Andreas etAl., 2004
! #Za     Z0mBSn = 0.035*u2star      *Grav_I
 
! Z0 Smooth + Saltat. Regime (Andreas etAl., 2004
! #Za     Z0enSV(ikl,ikv) =  Z0m_nu                                     &!
! #Za&                +  Z0mBSn                                          !
 
! Z0 Rough  Regime over Snow (Andreas etAl., 2004
! (.NOT. used by Erosion)     ams.confex.com/ams/pdfpapers/68601.pdf)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
! #Za     Z0m_90 =(10.-0.025*VVs_SV(ikl,ikv)/5.)                        &!
! #Za&            *exp(-0.4/sqrt(.00275+.00001*max(0.,VVs_SV(ikl,ikv)-5.)))  !
! #Za     Z0m_Sn =           DDs_SV(ikl,ikv)* Z0m_90 / 45.              &!
! #Za&         - DDs_SV(ikl,ikv)*DDs_SV(ikl,ikv)* Z0m_90 /(90.*90.)      !
 
! #ZA     u2star =      (us__SV(ikl,ikv) -0.1800)     / 0.1
! #ZA     Z0m_Sn =A_Fact*Z0mBSn *exp(-u2star*u2star)
 
! Z0 Rough  Regime over Snow (Andreas etAl., 2004
! #Za     u2star =      (us__SV(ikl,ikv) -0.1800)     / 0.1
! #Za     Z0m_Sn =A_Fact*Z0mBSn *exp(-u2star*u2star)
 
! Z0 Smooth + Saltat. Regime + Rough  Regime over Snow (Andreas etAl., 2004)
! #Za     Z0enSV(ikl,ikv) =  Z0enSV(ikl,ikv)                            &!
! #Za&                +  Z0m_Sn                                         !
 
! Z0               over Snow (instantaneous or time average)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0e_SV(ikl,ikv) =  Z0enSV(ikl,ikv)
! #ZM     Z0e_SV(ikl,ikv) =  Z0emSV(ikl,ikv)
 
! Momentum  Roughness Length
! ^^^^^^^^^^^^^^^^^^^^^^^^^^                                 ! Contribution of
          Z0mnSV(ikl,ikv) =  Z0mLnd                         &! Vegetation Form
     &                + (Z0m_Sn                             &! Sastrugi   Form
     &                +  Z0enSV(ikl,ikv))   *SnoMsk          ! Snow    Erosion
 
! Mom. Roughness Length, Discrimination among Ice/Land  and Ice-Free Ocean
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0mnSV(ikl,ikv) =  Z0mnSV(ikl,ikv)    *LISmsk     &! Ice and  Land
     &                  +Z0mSea      *(1-LISmsk)            &! Ice-Free Ocean
! #ZO&                  +Z0roSV(ikl,ikv)                    &! Subgrid  Topogr.
     &                  +0.
 
! GIS  Roughness Length
! ^^^^^^^^^^^^^^^^^^^^^
! #ZG     Z0mnSV(ikl,ikv) =                                            &!
! #ZG&      (1-LSmask(ikl,ikv)) *     Z0mnSV(ikl,ikv)                  &!
! #ZG&    +    LSmask(ikl,ikv)  * max(Z0mnSV(ikl,ikv),max(Z0_GIM,      &!
! #ZG&                                                    Z0_GIM+      &!
! #ZG&      (0.0032-Z0_GIM)*(ro__SV(ikl,ikv,isnoSV(ikl,ikv))-600.)     &!
! #ZG&                     /(920.00                 -600.)))            !
 
! Mom. Roughness Length, Instantaneous OR Box Moving Average in Time
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0m_SV(ikl,ikv) =  Z0mnSV(ikl,ikv)                   ! Z0mnSV  instant.
! #ZM     Z0m_SV(ikl,ikv) =  Z0mmSV(ikl,ikv)                   ! Z0mnSV  Average
 
! Corrected Threshold Friction Velocity before Erosion         ! Marticorena and
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^         ! Bergametti 1995
! #BS     Z0e_SV(ikl,ikv) =   min(Z0m_SV(ikl,ikv),Z0e_SV(ikl,ikv)) !
! #MB     f_eff=    log(0.35*(0.1        /Z0e_SV(ikl,ikv))**0.8)   ! JGR 100
! #MB     f_eff=1.-(log(      Z0m_SV(ikl,ikv)/Z0e_SV(ikl,ikv)      )) &! (20) p. 16420
! #MB&            /(max(      f_eff      ,eps6             ))  ! p.16426 2nd ?
! #MB     f_eff=    max(      f_eff      ,eps6              )  ! CONTROL
! #Mb     f_eff=2.0*max(      f_eff      ,eps6              )  ! TUNING
! #MB     f_eff=    min(      f_eff      ,un_1              )  !
! #MB     usthSV(ikl,ikv) =       usthSV(ikl,ikv)/f_eff        !
 
 
! Roughness Length for Scalars
! ----------------------------
 
          Z0hnSV(ikl,ikv) =     Z0mnSV(ikl,ikv)/  7.4
 
! Roughness Length for Scalars: Modification from Hapex-Sahel data
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Z0hnSV(ikl,ikv) =     Z0mnSV(ikl,ikv)/100.0
!                           Z0h = Z0m  /100.0   over the Sahel
!                                              (Taylor & Clark, QJRMS 127,p864)
 
! Roughness Length for Scalars: Modification for  Snow & Ice (Andrea, 1987)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #RS     rstar       =     Z0mnSV(ikl,ikv) * us__SV(ikl,ikv) / A_MolV
! #RS     rstar       = max(eps6,min(rstar,R_1000))
! #RS     alors       =          log(rstar)
! #RS     rstar0      = 1.250e0 * max(zer0,sign(un_1,0.135e0 - rstar))   &
! #RS&                +(1.      - max(zer0,sign(un_1,0.135e0 - rstar)))  &
! #RS&                *(0.149e0 * max(zer0,sign(un_1,2.500e0 - rstar))   &
! #RS&                + 0.317e0                                          &
! #RS&                *(1.      - max(zer0,sign(un_1,2.500e0 - rstar))))
! #RS     rstar1      = 0.      * max(zer0,sign(un_1,0.135e0 - rstar))   &
! #RS&                +(1.      - max(zer0,sign(un_1,0.135e0 - rstar)))  &
! #RS&                *(-0.55e0 * max(zer0,sign(un_1,2.500e0 - rstar))   &
! #RS&                - 0.565                                            &
! #RS&                *(1.      - max(zer0,sign(un_1,2.500e0 - rstar))))
! #RS     rstar2      = 0.      * max(zer0,sign(un_1,0.135e0 - rstar))   &
! #RS&                +(1.      - max(zer0,sign(un_1,0.135e0 - rstar)))  &
! #RS&                *(0.      * max(zer0,sign(un_1,2.500e0 - rstar))   &
! #RS&                - 0.183                                            &
! #RS&                *(1.00    - max(zer0,sign(un_1,2.500e0 - rstar))))
! #RS     Z0hnSV(ikl,ikv) = max(zer0                                     &
! #RS&                , sign(un_1,zzsnsv(ikl,ikv,isnoSV(ikl,ikv))-eps6)) &
! #RS&                * exp(rstar0+rstar1*alors+rstar2*alors*alors)      &
! #RS&                * 0.001e0 + Z0hnSV(ikl,ikv) * ( 1. - max(zer0      &
! #RS&                , sign(un_1,zzsnsv(ikl,ikv,isnoSV(ikl,ikv))-eps6)))
 
          Z0hnSV(ikl,ikv) =     Z0hSea             *(1-LISmsk)&! Ice-free Ocean
     &                +     Z0hnSV(ikl,ikv)        *   LISmsk  ! Ice and  Land
 
          Z0h_SV(ikl,ikv) =     Z0hnSV(ikl,ikv)
! #ZM     Z0h_SV(ikl,ikv) =     Z0hmSV(ikl,ikv)
 
 
! Contributions of the Roughness Lenghths to the neutral Drag Coefficient
! -----------------------------------------------------------------------
 
          IF (Garrat)                                                   &
     &    Z0m_SV(ikl,ikv) = max(2.0e-6     ,Z0m_SV(ikl,ikv)) ! Min Z0_m (Garrat Scheme)
          Z0m_SV(ikl,ikv) = min(Z0m_SV(ikl,ikv),za__SV(ikl,ikv)*0.3333)
          sqrCm0(ikl,ikv) = log(za__SV(ikl,ikv)/Z0m_SV(ikl,ikv))
! Martin control
!          PRINT*,'za__SV(:,:)=',za__SV(:,:)
!          PRINT*,'Z0h_SV=',Z0h_SV
! Martin control
          sqrCh0(ikl,ikv) = log(za__SV(ikl,ikv)/Z0h_SV(ikl,ikv))
 
! OUTPUT of SnowFall, Roughness Length and Drag Coefficients
! #sf     IF (ikl,ikv.EQ.1) write(6,6661) dsn_SV(ikl,ikv),us__SV(ikl,ikv),Z0SaSi&
! #sf&                        ,Z0Sa_N,Z0SaSV(ikl,ikv),Z0m_Sn,Z0m_SV(ikl,ikv)
 6661     format(20x,7f9.6)
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           of Roughness Length and Drag Coefficients
! #sz     IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ. jwr_SV .AND.  &
! #sz&        ikv        .EQ.nwr_SV)                                    &
! #sz&    write(6,6600)  za__SV(ikl,ikv) , Z0m_SV(ikl,ikv)              &
! #sz&                  ,sqrCm0(ikl,ikv) , za__SV(ikl,ikv)/Z0m_SV(ikl,ikv)  &
! #sz&                  ,Z0SaSV(ikl,ikv) , Z0h_SV(ikl,ikv)              &
! #sz&                  ,sqrCh0(ikl,ikv) , za__SV(ikl,ikv)/Z0h_SV(ikl,ikv)
 6600     format(/,' ** SISVAT     *0  '                                &
     &            ,'  za__SV  = ',e12.4,'  Z0m_SV  = ',e12.4            &
     &            ,'  sqrCm0  = ',e12.4,'  Za/Z0m  = ',e12.4            &
     &          ,/,'                   '                                &
     &            ,'  Z0SaSV  = ',e12.4,'  Z0h_SV  = ',e12.4            &
     &            ,'  sqrCh0  = ',e12.4,'  Za/Z0h  = ',e12.4)
 
 
! Vertical Stability Correction
! -----------------------------
 
! Surface/Canopy Temperature
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          Tsrfsv(ikl,ikv) = Sigmsv(ikl,ikv) * TvegSV(ikl,ikv)           &
     &          + (1. - Sigmsv(ikl,ikv))* TsisSV(ikl,ikv,isnoSV(ikl,ikv))
        END DO
        END DO
 
! Aerodynamic Resistance
! ^^^^^^^^^^^^^^^^^^^^^^
        IF            (SnoMod.AND.ColPrt)                           THEN
 
!                **********
            call ColPrt_SBL
!                **********
 
        ELSE
         IF           (Garrat)                                      THEN
 
!                **********
            call SISVAT_SBL
!                **********
 
         ELSE
 
!                **********
            call SISVATeSBL
!                **********
 
         END IF
 
         DO ikl=1,kcolp
         DO ikv=1,mwp
           IF (LMO_SV(ikl,ikv) .GT. 0.) LMO_SV(ikl,ikv) = max( epsLMO,LMO_SV(ikl,ikv))
           IF (LMO_SV(ikl,ikv) .LT. 0.) LMO_SV(ikl,ikv) = min(-epsLMO,LMO_SV(ikl,ikv))
         END DO
         END DO
        END IF
 
 
! Canopy Energy Balance
! =====================
 
!            **********
        call SISVAT_TVg(                                                &
! #e1&                  ETVg_d                                          &
     &                 )
!            **********
 
 
! Surface/Canopy Temperature
! ==========================
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          Tsrfsv(ikl,ikv) = Sigmsv(ikl,ikv) * TvegSV(ikl,ikv)           &
     &          + (1. - Sigmsv(ikl,ikv))* TsisSV(ikl,ikv,isnoSV(ikl,ikv))
        END DO
        END DO
 
 
! Soil   Energy Balance
! =====================
 
 
!            **********
        call SISVAT_TSo(                                                &
! #e1&                  ETSo_0     ,ETSo_1     ,ETSo_d     ,kcolw       &
     &                 )
!            **********
 
!            **********
! #ve   call SISVAT_wEq('_TSo  ',0)
!            **********
 
 
 
 
! Canopy Water  Balance
! =====================
 
! Soil Water     Potential
! ------------------------
 
      DO   isl=-nsoil,0
        DO ikl=1,kcolp
        DO ikv=1,mwp
          ist             =     isotSV(ikl,ikv)    ! Soil Type
          psi_sv(ikl,ikv,isl) =     psidSV(ist)   &! DR97, Eqn.(3.34)
     &  *(etadSV(ist) /max(eps6,eta_SV(ikl,ikv,isl))) &!
     &  **bCHdSV(ist)                              !
 
 
! Soil Hydraulic Conductivity
! ---------------------------
 
          Khydsv(ikl,ikv,isl) =    s2__SV(ist)    &! DR97, Eqn.(3.35)
     &  *(eta_SV(ikl,ikv,isl)**(2.*bCHdSV(ist)+3.))!
        END DO
        END DO
      END DO
 
!            **********
        call SISVAT_qVg
!            **********
 
 
! OUTPUT/Verification: H2O    Conservation: Vegetation Forcing
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Watsvd(ikl,ikv) =     (Watsvd(ikl,ikv)  &! Canopy Precip. IN
! #m0&                      -drr_SV(ikl,ikv)      &! Canopy Precip. OUT
! #m0&                      -Evp_sv(ikl,ikv))* dt__SV  ! Canopy Water Evap.
! #m0   END DO
! #m0   END DO
 
 
! Melting / Refreezing in the Snow Pack
! =====================================
 
      IF (SnoMod)                                                   THEN
 
!            **********
        call SISVAT_qSn(                                                &
! #e1&                  EqSn_0,EqSn_1,EqSn_d                            &
! #m1&                 ,SIsubl,SImelt,SIrnof                            &
     &                 )
!            **********
 
!            **********
! #ve   call SISVAT_wEq('_qSn  ',0)
!            **********
 
! OUTPUT in SISVAT for ikl = 1 (preferably for Stand Alone Version)
! OUTPUT           for SnowFall and Snow Buffer
! #s2   IF          (isnoSV(1,1) .GT. 0)                                &
! #s2&  write(6,6007)isnoSV(1,1),  dsn_SV(1) *dt__SV + BufsSV(1,1),     &
! #s2&              (dzsnSV(1,isn)*ro__SV(1,isn),isn=1,isnoSV(1,1))
 6007   format(i3,'  dsn+Buf=',f6.2,6x,'q dz *ro =',10f6.2,             &
     &                                       (/,35x,10f6.2))
 
 
! Snow Pack Thickness
! -------------------
 
          DO ikl=1,kcolp
          DO ikv=1,mwp
            z_snsv(ikl,ikv)     = 0.0
          END DO
          END DO
        DO   isn=1,nsnow
          DO ikl=1,kcolp
          DO ikv=1,mwp
            z_snsv(ikl,ikv)     = z_snsv(ikl,ikv) + dzsnSV(ikl,ikv,isn)
            zzsnsv(ikl,ikv,isn) = z_snsv(ikl,ikv)
          END DO
          END DO
        END DO
 
 
! Energy in Excess is added to the first Soil Layer
! -------------------------------------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
            z_snsv(ikl,ikv)   = max(zer0,                               &
     &                          sign(un_1,eps6-z_snsv(ikl,ikv)))
            TsisSV(ikl,ikv,0) = TsisSV(ikl,ikv,0)    + EExcsv(ikl,ikv)  &
     &                                       /(rocsSV(isotSV(ikl,ikv))  &
     &                                        +rcwdSV*eta_SV(ikl,ikv,0))
            EExcsv(ikl,ikv)   = 0.
        END DO
        END DO
 
 
! OUTPUT/Verification: * Mass Conservation: Mass (below the Canopy) and Forcing
! #m1   DO ikl=1,kcolp
! #m1   DO ikv=1,mwp
! #m1     SIWa_f(ikl,ikv) =(drr_SV(ikl,ikv) + dsn_SV(ikl,ikv))   *dt__SV ![mm w.e.]
! #m1     SIWe_f(ikl,ikv) = dbs_SV(ikl,ikv)                          !
! #m1     SIWm_1(ikl,ikv) = BufsSV(ikl,ikv) + HFraSV(ikl,ikv)    *rhoIce !
! #m1   DO isn=1,nsnow                                               !
! #m1     SIWm_1(ikl,ikv) = SIWm_1(ikl,ikv) + dzsnSV(ikl,ikv,isn)*ro__SV(ikl,ikv,isn)!
! #m1   END DO                                                       !
! #m1   END DO
! #m1   END DO                                                       !
 
      END IF
 
 
! Soil   Water  Balance
! =====================
 
!            **********
        call SISVAT_qSo(                                                &
! #m0&                 (Wats_0,Wats_1,Wats_d                            &
     &                 )
!            **********
 
 
! Surface/Canopy Fluxes
! =====================
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          IRdwsv(ikl,ikv)=tau_sv(ikl,ikv) *IRd_SV(ikl,ikv)*Eso_sv(ikl,ikv)  &! Downward IR
     &          +(1.0-tau_sv(ikl,ikv))*IRd_SV(ikl,ikv)*Evg_sv(ikl,ikv)   !
          IRupsv(ikl,ikv) =      IRupsv(ikl,ikv)                &! Upward   IR
     &                + 0.5 *IRv_sv(ikl,ikv) * (1.-tau_sv(ikl,ikv))  !
          IRu_SV(ikl,ikv) =     -IRupsv(ikl,ikv)                &! Upward   IR
     &                      +IRd_SV(ikl,ikv)                    &! (effective)
     &                      -IRdwsv(ikl,ikv)                     ! (positive)
          TBr_sv(ikl,ikv) =sqrt(sqrt(IRu_SV(ikl,ikv)/StefBo))    ! Brightness
!                                                                ! Temperature
          uts_SV(ikl,ikv) =     (HSv_sv(ikl,ikv) +HSs_sv(ikl,ikv))  &! u*T*
     &                     /(rhT_SV(ikl,ikv) *CpdAir)            !
          uqs_SV(ikl,ikv) =     (HLv_sv(ikl,ikv) +HLs_sv(ikl,ikv))  &! u*q*
     &                     /(rhT_SV(ikl,ikv) *LhvH2O)            !
 
! Surface/Canopy Temperature
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          Tsrfsv(ikl,ikv) = Sigmsv(ikl,ikv) * TvegSV(ikl,ikv)           &
     &          + (1. - Sigmsv(ikl,ikv))* TsisSV(ikl,ikv,isnoSV(ikl,ikv))
        END DO
        END DO
 
 
! Snow Pack Properties (sphericity, dendricity, size)
! ===================================================
 
      IF (SnoMod)                                                   THEN
 
!            **********
        call SISVAT_GSn
!            **********
 
!            **********
! #ve   call SISVAT_wEq('_GSn  ',0)
!            **********
 
 
! Surficial Water Freezing, including that of a Water Surface (isotSV=0)
! ======================================================================
 
 
      END IF
 
 
! OUTPUT
! ======
 
      IF (kcolv.LE.mwp)                                             THEN
          ikl=1
        write(4,4) Day_TU,LabMon(Mon_TU),YearTU,HourTU,MinuTU          &
     &            ,(ivgtSV(ikl,ikv),DH_dSV(ivgtSV(ikl,ikv))            &
     &                           ,1.e3*Z0m_SV(ikl,ikv),ikv,ikv=1,mwp)
 4      format(i3,a3,i5,i3,'h',i3,'   Vegetation: ',9(i6,2f8.3,i3))
      END IF
 
! #e0   DO ikl=1,kcolp
! #e0   DO ikv=1,mwp
! #e0   IF (lwriSV(ikl,ikv).ne.0)                                   THEN
! #e0           noUNIT =  no__SV(lwriSV(ikl,ikv))
! #e0     write(noUNIT,5001)                                            &
! #e0&       (SoSosv(ikl,ikv)+SoCasv(ikl,ikv))*sol_SV(ikl,ikv),         &
! #e0&        IRdwsv(ikl,ikv),IRu_SV(ikl,ikv),                          &
! #e0&        HSv_sv(ikl,ikv)+HSs_sv(ikl,ikv),                          &
! #e0&        HLv_sv(ikl,ikv)+HLs_sv(ikl,ikv), TaT_SV(ikl,ikv),         &
! #e0&        dsn_SV(ikl,ikv)*3.6e3,       drr_SV(ikl,ikv)*3.6e3,       &
! #e0&        SoSosv(ikl,ikv)             *sol_SV(ikl,ikv),             &
! #e0&                    IRv_sv(ikl,ikv) *0.5,                         &
! #e0&        HSv_sv(ikl,ikv),HLv_sv(ikl,ikv), TvegSV(ikl,ikv),         &
! #e0&                    SoCasv(ikl,ikv) *sol_SV(ikl,ikv),             &
! #e0&        HSs_sv(ikl,ikv),HLs_sv(ikl,ikv), TsisSV(ikl,ikv,isnoSV(ikl,ikv))
 5001     format(                                                       &
     &         '        |Net Solar| IR Down | IR Up   | HS/Dwn=+|',     &
     &          ' HL/Dwn=+| Temper. |         |  Snow  |  Rain  |',     &
     &       /,'        | [W/m2]  | [W/m2]  | [W/m2]  | [W/m2]  |',     &
     &          ' [W/m2]  | [K]     |         | [mm/h] | [mm/h] |',     &
     &       /,' -------+',7('---------+'),2('--------+'),              &
     &       /,' SISVAT |',f8.1,' |',f8.1,' |',f8.1,' |',f8.1,' |',     &
     &                     f8.1,' |A',f7.2,' |', 8x ,' |',2(f7.2,' |'), &
     &       /,' Canopy |',f8.1,' |', 8x ,' |',f8.1,' |',f8.1,' |',     &
     &                     f8.1,' |',f8.2,' |', 8x ,' |',2( 7x ,' |')   &
     &       /,' Soil   |',f8.1,' |', 8x ,' |', 8x ,' |',f8.1,' |',     &
     &                     f8.1,' |',f8.2,' |', 8x ,' |',2( 7x ,' |'))
 
 
! OUTPUT/Verification: Energy/Water Budget
! #e1     Enrsvd(ikl,ikv) = Enrsvd(ikl,ikv)                              &! Up Surf. IR
! #e1&                + IRs_SV(ikl,ikv)                                  &! Offset
! #e1&        + (      (SoSosv(ikl,ikv)                                  &! Net   Solar
! #e1&                 +SoCasv(ikl,ikv)) *sol_SV(ikl,ikv)                &!
! #e1&          +                     IRdwsv(ikl,ikv)                    &! Downward IR
! #e1&          +                     IRupsv(ikl,ikv)                    &! Upward   IR
! #e1&          +                     HSv_sv(ikl,ikv)+HSs_sv(ikl,ikv)    &! Sensible
! #e1&          +                     HLv_sv(ikl,ikv)+HLs_sv(ikl,ikv))    ! Latent
 
! #e1     write(noUNIT,5002)           Enrsvd(ikl,ikv),                  &
! #e1&                   ETSo_0(ikl,ikv),  ETSo_d(ikl,ikv),              &
! #e1&                   ETSo_0(ikl,ikv)+  ETSo_d(ikl,ikv), ETSo_1(ikl,ikv), &
! #e1&                   EqSn_0(ikl,ikv)                            /dt__SV, &
! #e1&                   EqSn_d(ikl,ikv)                            /dt__SV, &
! #e1&                  (EqSn_1(ikl,ikv)-  EqSn_0(ikl,ikv)- EqSn_d(ikl,ikv))/dt__SV, &
! #e1&                   EqSn_1(ikl,ikv)                            /dt__SV
 5002     format(                                                        &!
     &           ' -----------------+-------------------+',              &!
     &            '-----------------+-+-----------------+',              &!
     &          '-------------------+',                                  &!
     &         /,' SOIL/SNOW/VEGET. |                   |',              &!
     &            ' Power,  Forcing |                   |',              &! Enrsvd
     &          '                   |',                                  &!
! #el&         /,' -----------------+-------------------+',              &!
! #el&            '-----------------+-------------------+',              &!
! #el&          '-------------------+',                                  &!
     &         /,'                  |',   11x ,'        |',              &!
     &                f9.2,' [W/m2] |',   11x ,'        |',              &! Enrsvd
     &                11x ,'        |',                                  &!
     &         /,' -----------------+-------------------+',              &!
     &            '-----------------+-------------------+',              &!
     &          '-------------------+',                                  &!
     &         /,' SOIL/SNOW  (TSo) | Energy/dt, Time 0 |',              &!        ETSo_0
     &            ' Power,  Forcing |   Sum Tim.0+Forc. |',              &! ETSo_d/ETSo_0+d
     &          ' Energy/dt, Time 1 |',                                  &! ETSo_1
! #el&         /,' -----------------+-------------------+',              &!
! #el&            '-----------------+-------------------+',              &!
! #el&          '-------------------+',                                  &!
     &         /,'                  |',  f11.2,' [W/m2] |',              &!        ETSo_0
     &                f9.2,' [W/m2] |',  f11.2,' [W/m2] |',              &! ETSo_d/ETSo_0+d
     &               f11.2,' [W/m2] |',                                  &! ETSo_1
     &         /,' -----------------+-------------------+',              &!
     &            '-----------------+-------------------+',              &!
     &          '-------------------+',                                  &!
     &         /,'      SNOW  (qSn) | Energy/dt, Time 0 |',              &! EqSn_0/dt
     &            ' Power,  Excess  |   D(Tim.1-0-Forc.)|',              &! EqSn_d/dt, 1-0-d
     &          ' Energy/dt, Time 1 |',                                  &! EqSn_1/dt
! #el&         /,' -----------------+-------------------+',              &!
! #el&            '-----------------+-------------------+',              &!
! #el&          '-------------------+',                                  &!
     &         /,'                  |',  f12.2, '[W/m2] |',              &! EqSn_0/dt
     &                f9.2,' [W/m2] |',  f11.2,' [W/m2] |',              &! EqSn_d/dt, 1-0-d
     &               f12.2, '[W/m2] | ',                                 &! EqSn_1/dt
     &         /,' -----------------+-------------------+',              &!
     &            '-----------------+-------------------+',              &!
     &          '-------------------+')                                   !
 
! #e1             EnsBal = ETSo_1(ikl,ikv)-(ETSo_0(ikl,ikv)+Enrsvd(ikl,ikv))
! #e1             EnvBal = Enrsvd(ikl,ikv)- ETVg_d(ikl,ikv)
! #e1     IF (abs(EnsBal).gt.5.e-1                                       &
! #e2&    .OR.lwriSV(ikl,ikv).eq.    2                                   &
! #e1&                            )                                    THEN
! #e1       write(6,6001) daHost,i___SV(lwriSV(ikl,ikv)),                &
! #e1&                           j___SV(lwriSV(ikl,ikv)),                &
! #e1&                           n___SV(lwriSV(ikl,ikv)),                &
! #e1&                           ETSo_1(ikl,ikv),ETSo_0(ikl,ikv),ETSo_d(ikl,ikv),&
! #e1&                           ETSo_1(ikl,ikv)-ETSo_0(ikl,ikv)-ETSo_d(ikl,ikv),&
! #e1&                           Enrsvd(ikl,ikv),ETVg_d(ikl,ikv),ETSo_d(ikl,ikv),&
! #e1&                           Enrsvd(ikl,ikv)-ETVg_d(ikl,ikv)-ETSo_d(ikl,ikv)
 6001       format(a18,3i4,' (EB1'           ,f15.6,                     &
     &                ')  - [(EB0           ',f15.6,')',                 &
     &               /,55x,'+(ATM->Snow/Soil',f15.6,')] ',               &
     &                     '= EBAL'          ,f15.6,' [W/m2]',           &
     &               /,55x,' (ATM->SISVAT'   ,f18.6,                     &
     &               /,55x,'- Veg. ImBal.',   f18.6,')  ',               &
     &               /,55x,'- ATM->SnoSol',   f18.6,')  ',               &
     &                     '= ????'          ,f15.6,' [W/m2]')
! #e1           noEBal = noEBal + 1
! #e2           noEBal = noEBal - 1
! #e1       IF (noEBal.GE.       10) stop 'TOO MUCH ENERGY IMBALANCES'
! #e1     END IF
 
 
! OUTPUT/Verification: * Mass Conservation: Budget [mm w.e.]
! #m1     write(noUNIT,5010)                                             &
! #m1&             SIWm_0(ikl,ikv),  SIWa_i(ikl,ikv)-SIWa_f(ikl,ikv)     &
! #m1&            ,SIWm_0(ikl,ikv)+  SIWa_i(ikl,ikv)-SIWa_f(ikl,ikv)     &
! #m1&                          +SIWe_i(ikl,ikv)-SIWe_f(ikl,ikv)         &
! #m1&                          +SIsubl(ikl,ikv)                         &
! #m1&                          -SImelt(ikl,ikv)                         &
! #m1&                          -SIrnof(ikl,ikv)                         &
! #m2&                          +SIvAcr(ikl,ikv)                         &
! #m1&            ,SIWm_1(ikl,ikv),  SIWe_i(ikl,ikv)-SIWe_f(ikl,ikv)     &
! #m1&            ,              SIsubl(ikl,ikv)                         &
! #m1&            ,             -SImelt(ikl,ikv)                         &
! #m1&            ,             -SIrnof(ikl,ikv)                         &
! #m2&            ,              SIvAcr(ikl,ikv)
 5010     format(' SNOW             |   Snow,   Time 0  |',              &
     &            ' Snow,   Forcing |           Sum     |',              &
     &          '   Snow,   Time 1  |',                                  &
! #el&         /,' -----------------+-------------------+',              &
! #el&            '-----------------+-------------------+',              &
! #el&          '-------------------+',                                  &
     &         /,'                  |',    f13.3,' [mm] |',              &
     &           ' A',  f9.3,' [mm] |',    f13.3,' [mm] |',              &
     &                 f13.3,' [mm] |',                                  &
     &         /,'                  |',     13x ,'      |',              &
     &           ' E',  f9.3,' [mm] |',     13x ,'      |',              &
     &                  13x ,'      |',                                  &
     &         /,'                  |',     13x ,'      |',              &
     &           ' S',  f9.3,' [mm] |',     13x ,'      |',              &
     &                  13x ,'      |',                                  &
     &         /,'                  |',     13x ,'      |',              &
     &           '(M',  f9.3,' [mm])|  (included in A)  |',              &
     &                  13x ,'      |',                                  &
     &         /,'                  |',     13x ,'      |',              &
     &           ' R',  f9.3,' [mm] |',     13x ,'      |',              &
     &                  13x ,'      |',                                  &
! #m2&         /,'                  |',     13x ,'      |',              &
! #m2&           ' O',  f9.3,' [mm] |',     13x ,'      |',              &
! #m2&                  13x ,'      |',                                  &
     &         /,' -----------------+-------------------+',              &
     &            '-----------------+-------------------+',              &
     &          '-------------------+')
! #m1             SnoBal = SIWm_1(ikl,ikv)-(SIWm_0(ikl,ikv)              &
! #m1&                                 +SIWa_i(ikl,ikv)-SIWa_f(ikl,ikv)  &
! #m1&                                 +SIWe_i(ikl,ikv)-SIWe_f(ikl,ikv)) &
! #m1&                                 -SIsubl(ikl,ikv)                  &
! #m1&                                 +SIrnof(ikl,ikv)                  &
! #m2&                                 -SIvAcr(ikl,ikv)
! #m1     IF (abs(SnoBal).gt.eps6)                                  THEN
! #m1       write(6,6010) daHost,i___SV(lwriSV(ikl,ikv)),                &
! #m1&                           j___SV(lwriSV(ikl,ikv)),                &
! #m1&                           n___SV(lwriSV(ikl,ikv)),                &
! #m1&                           SIWm_1(ikl,ikv),SIWm_0(ikl,ikv),        &
! #m1&                           SIWa_i(ikl,ikv),SIWa_f(ikl,ikv),        &
! #m1&                           SIWe_i(ikl,ikv),SIWe_f(ikl,ikv),        &
! #m1&                           SIsubl(ikl,ikv),SImelt(ikl,ikv),        &
! #m2&                           SIrnof(ikl,ikv),SIvAcr(ikl,ikv),        &
! #m1&                           SnoBal
 6010       format(a18,3i4,' (MB1'        ,f12.6,                        &
     &                 ') - [(MB0        ',f12.6,        15x,')',        &
     &               /,51x,'+(ATM Forcing',f12.6,' - ',f12.6,')',        &
     &               /,51x,'+(BLS Forcing',f12.6,' - ',f12.6,')',        &
     &               /,51x,'-(Depo/Sublim',f12.6,        15x,')',        &
     &               /,51x,' !Melting    ',f12.6,'  included in A!',     &
     &               /,51x,'+(Run  OFF   ',f12.6,        15x,')',        &
! #m2&               /,51x,'-(Sea-Ice Acr',f12.6,        15x,')',        &
     &               /,29x,'= *BAL'       ,f12.6,      ' [mm w.e.]')
! #m1           noSBal = noSBal + 1
! #m1       IF (noSBal.GE.       10) stop 'TOO MUCH SNOW MASS IMBALANCE'
! #m1     END IF
 
 
! OUTPUT/Verification: H2O    Conservation: Water  Budget
! #m0     Watsv0(ikl,ikv) =  Watsv0(ikl,ikv)           &! Canopy Water Cont.
! #m0&                 + Wats_0(ikl,ikv)                ! Soil   Water Cont.
! #m0     Watsvd(ikl,ikv) =  Watsvd(ikl,ikv)           &! Canopy Forcing
! #m0&                 + Wats_d(ikl,ikv)                ! Soil   Forcing
 
! #m0     write(noUNIT,5003)                                             &
! #m0&                   Wats_0(ikl,ikv),  Wats_d(ikl,ikv),              &
! #m0&                   Wats_0(ikl,ikv)+  Wats_d(ikl,ikv),    Wats_1(ikl,ikv),  &
! #m0&                   Watsv0(ikl,ikv),  Watsvd(ikl,ikv),              &
! #m0&                   Watsv0(ikl,ikv)+  Watsvd(ikl,ikv),    Wats_1(ikl,ikv)   &
! #m0&                                                +rrCaSV(ikl,ikv)
 5003     format(' SOIL/SNOW  (qSo) |   Water,  Time 0  |',              &
     &            ' Water,  Forcing |           Sum     |',              &
     &          '   Water,  Time 1  |',                                  &
! #el&         /,' -----------------+-------------------+',              &
! #el&            '-----------------+-------------------+',              &
! #el&          '-------------------+',                                  &
     &         /,'                  |',    f13.3,' [mm] |',              &
     &                 f11.3,' [mm] |',    f13.3,' [mm] |',              &
     &                 f13.3,' [mm] |',                                  &
     &         /,' -----------------+-------------------+',              &
     &            '-----------------+-------------------+',              &
     &          '-------------------+',                                  &
     &         /,' SOIL/SNOW/VEGET. |   Water,  Time 0  |',              &
     &            ' Water,  Forcing |           Sum     |',              &
     &          '   Water,  Time 1  |',                                  &
! #el&         /,' -----------------+-------------------+',              &
! #el&            '-----------------+-------------------+',              &
! #el&          '-------------------+',                                  &
     &         /,'                  |',    f13.3,' [mm] |',              &
     &                 f11.3,' [mm] |',    f13.3,' [mm] |',              &
     &                 f13.3,' [mm] |',                                  &
     &         /,' -----------------+-------------------+',              &
     &            '-----------------+-------------------+',              &
     &          '-------------------+')
 
! #m0             WatBal = Wats_1(ikl,ikv)+rrCaSV(ikl,ikv)               &
! #m0&                   -(Watsv0(ikl,ikv)+Watsvd(ikl,ikv))
! #m0     IF (abs(WatBal).gt.eps6)                                  THEN
! #m0       write(6,6002) daHost,i___SV(lwriSV(ikl,ikv)),                &
! #m0&                           j___SV(lwriSV(ikl,ikv)),                &
! #m0&                           n___SV(lwriSV(ikl,ikv)),                &
! #m0&                           Wats_1(ikl,ikv),rrCaSV(ikl,ikv),        &
! #m0&                           Watsv0(ikl,ikv),Watsvd(ikl,ikv),WatBal, &
! #m0&                           Wats_1(ikl,ikv),                        &
! #m0&                           Wats_0(ikl,ikv),Wats_d(ikl,ikv),        &
! #m0&               Wats_1(ikl,ikv)-Wats_0(ikl,ikv)-Wats_d(ikl,ikv)
 6002       format(30x,' NEW Soil Water',3x,' Canopy   Water',3x,        &
     &                 ' OLD SVAT Water',4x,' FRC SVAT Water',           &
     &           /,a18,3i4,f15.6,' + ' ,f15.6,' - ' ,f15.6,              &
     &                           ' -  ',f15.6,'    ', 15x ,'    ',       &
     &      /,31x,'= ',f12.6,' [mm] (Water Balance)',                    &
     &           /,30x,' NEW Soil Water',3x,'               ',3x,        &
     &                 ' OLD Soil Water',4x,' FRC Soil Water',           &
     &           /,30x,f15.6,'   ' , 15x ,' - ' ,f15.6,                  &
     &                       ' -  ',f15.6,'    ', 15x ,'    ',           &
     &      /,31x,'= ',f12.6,' [mm] (3 terms SUM)')
! #m0           noWBal = noWBal + 1
! #m0       IF (noWBal.GE.       10) stop 'TOO MUCH WATER  IMBALANCES'
! #m0     END IF
 
 
! Water/Temperature Profiles
! --------------------------
 
! #e0       write(noUNIT,5004)
 5004       format(' -----+--------+--+-----+--------+----+---+',        &
     &  '--------+----+---+--------+------+-+--------+--------+',        &
     &           /,'    n |     z  |     dz |     ro |    eta |',        &
     &  '     T  |     G1 |     G2 | Extinc |        | HISTORY|',        &
     &           /,'      |    [m] |    [m] | [kg/m3]| [m3/m3]|',        &
     &  '    [K] |    [-] |    [-] |    [-] |        |   [-]  |',        &
     &           /,' -----+--------+--------+--------+--------+',        &
     &  '--------+--------+--------+--------+--------+--------+')
! #e0       write(noUNIT,5005) rusnSV(ikl,ikv),albisv(ikl,ikv)
 5005       format('      |        |        |        |W',f6.3,' |',      &
     &  '        |        |        |A',f6.3,' |        |        |')
! #e0       write(noUNIT,5015)                                           &
! #e0&                    (isn,zzsnsv(ikl,ikv,isn),dzsnSV(ikl,ikv,isn),  &
! #e0&                         ro__SV(ikl,ikv,isn),eta_SV(ikl,ikv,isn),  &
! #e0&                         TsisSV(ikl,ikv,isn),                      &
! #e0&                         G1snSV(ikl,ikv,isn),G2snSV(ikl,ikv,isn),  &
! #e0&                         sEX_sv(ikl,ikv,isn),istoSV(ikl,ikv,isn),  &
! #e0&                     isn=isnoSV(ikl,ikv),1,-1)
 5015       format((i5,' |',2(f7.3,' |'),               f7.1,' |',       &
     &           f7.3,' |' ,  f7.2,' |', 2(f7.1,' |'),  f7.3,' |',       &
     &            7x ,' |' ,  i5,'   |'                          ))
! #e0       write(noUNIT,5006)
 5006       format(' -----+--------+--------+--------+--------+',        &
     &  '--------+--------+--------+--------+--------+--------+')
! #e0       write(noUNIT,5007) TBr_sv(ikl,ikv),
! #e0&                         TvegSV(ikl,ikv),rrCaSV(ikl,ikv)*1.e3,
! #e0&                         EvT_sv(ikl,ikv)*86.4e3
 5007       format(' Brgh |',4(8x,'|'),  f7.2,' | [micm] |',4(8x,'|'),   &
     &           /,' VEGE |',4(8x,'|'),2(f7.2,' |'),        2(8x,'|'),   &
     &                                   f7.3,' |',           8x,'|' )
! #e0       write(noUNIT,5014)
 5014       format(' -----+--------+--------+--------+--------+',        &
     &  '--------+--------+--------+--------+--------+--------+',        &
     &           /,'    n |        |     dz |        |    eta |',        &
     &  '     T  |        |        |        | Root W.| W.Flow |',        &
     &           /,'      |        |    [m] |        | [m3/m3]|',        &
     &  '    [K] |        |        |        | [mm/d] | [mm/h] |',        &
     &           /,' -----+--------+--------+--------+--------+',        &
     &  '--------+--------+--------+--------+--------+--------+')
 
! #e0       write(noUNIT,5008)                                           &
! #e0&                    (isl,    LSdzsv(ikl,ikv)*dz_dSV(    isl),      &
! #e0&                                         eta_SV(ikl,ikv,isl),      &
! #e0&                                         TsisSV(ikl,ikv,isl),      &
! #e0&                                  86.4e3*Rootsv(ikl,ikv,isl),      &
! #e0&                                   3.6e3*Khydsv(ikl,ikv,isl),      &
! #e0&                     isl=0,-nsoil,-1)
 5008       format((i5,' |',   7x ,' |' ,  f7.3,' |' ,   7x ,' |',       &
     &           f7.3,' |' ,  f7.2,' |', 2( 7x ,' |'),   7x ,' |',       &
     &                                     f7.3,' |' ,  f7.2,' |'))
! #e0       write(noUNIT,5006)
! #e0       write(noUNIT,5009) RnofSV(ikl,ikv)* 3.6e3
 5009       format('      |',9(8x,'|'),f7.3,' |')
! #e0       write(noUNIT,5006)
! #e0   END IF
! #e0   END DO
! #e0   END DO
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
! DE-ALLOCATION                                           !
! =============                                           !
 
      IF (FlagDALLOC)                                THEN !
 
      deallocate          ( TBr_sv )                      ! Brightness Temperature
      deallocate          ( IRdwsv )                      ! DOWNward   IR Flux
      deallocate          ( IRupsv )                      ! UPward     IR Flux
      deallocate          ( Bdzssv )                      ! Buffer Snow Layer Thickness
      deallocate          ( z_snsv )                      ! Snow-Ice, current Thickness
 
! Energy         Budget
! ~~~~~~~~~~~~~~~~~~~~~
! #e1 deallocate          ( ETVg_d )                      ! VegetationPower, Forcing
! #e1 deallocate          ( ETSo_1 )                      ! Soil/Snow Power, after  Forcing
! #e1 deallocate          ( EqSn_0 )                      ! Snow Energy, befor Phase Change
! #e1 deallocate          ( EqSn_1 )                      ! Snow Energy, after Phase Change
! #e1 deallocate          ( EqSn_d )                      ! Energy in Excess
 
! OUTPUT/Verification: H2O    Conservation
! #m0 deallocate          ( Wats_0 )                      ! Soil Water,  before Forcing
! #m0 deallocate          ( Wats_1 )                      ! Soil Water,  after  Forcing
! #m0 deallocate          ( Wats_d )                      ! Soil Water,         Forcing
 
! OUTPUT/Verification: * Mass Conservation
! #m1 deallocate          ( SIsubl )                      ! Snow Sublimed/Deposed Mass
! #m1 deallocate          ( SImelt )                      ! Snow Melted           Mass
! #m1 deallocate          ( SIrnof )                      ! Local Surficial Water + Run OFF
 
! OUTPUT/Verification: SeaIce Conservation
! #m2 deallocate          ( SIvAcr )                      ! Sea-Ice      Vertical Acretion
 
! Energy and Mass Budget
! ~~~~~~~~~~~~~~~~~~~~~~
! #e1 deallocate          ( Enrsvd )                      ! Soil+Vegetat  Power  Forcing
 
                                                          ! H2O    Conservation
! #m0 deallocate          ( Watsv0 )                      ! Soil+Vegetat, before Forcing
! #m0 deallocate          ( Watsvd )                      ! Soil+Vegetat  Water  Forcing
 
                                                          ! * Mass Conservation
! #m1 deallocate          ( SIWm_0       ,SIWm_1 )        ! Snow Initial/Final        Mass
! #m1 deallocate          ( SIWa_i       ,SIWa_f )        ! Snow Initial/Final ATM Forcing
! #m1 deallocate          ( SIWe_i       ,SIWe_f )        ! Snow Initial/Final BLS Forcing
 
 
      deallocate          ( IcIndx )                      ! No   Ice               Mask
 
      deallocate          ( FallOK )                      ! Snow Contribution to the Canopy
 
      END IF                                              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
! END  .main. (SISVAT)
 
 
      return
      end subroutine SISVAT
 
 
 
      subroutine SISVAT_BSn
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_BSn                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_BSn treats Snow Erosion and Deposition             |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Thu 14-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     Preprocessing  Option: STANDARD Possibility                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                          |
!     #BS: Explicit Cloud MICROPHYSICS: Blow. *(Snow)         Model        |
!     #BM: Explicit Cloud MICROPHYSICS: de Montmollin Parameterizat.       |
!     #MA: SNOW Model: Increased polar B* Mobility (Mann et al.2000)       |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #BA: Budd et al.            1966, Ant.Res.Ser.9    u* BS Threshold   |
!     #BY: Budd et al. 1966, 2~m Averag Blow. *(Snow)    Properties        |
!     #AG: Snow Aging Col de Porte     (Brun et al.1991) discard BS at CdP |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # stdout               | #s2: OUTPUT of SnowFall, Snow Buffer          |
!                          |      unit  6, SubRoutine  SISVAT_BSn, _qSn    |
!   # stdout               | #b0: OUTPUT of Snow Erosion Statistics        |
!                          |      unit  6, SubRoutine  SISVAT_BSn **ONLY** |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_ctr
      use Mod_SISVAT_dat
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
      use Mod_SISVAT_BSn
      use Mod_SISVATLBSn
 
 
 
      IMPLICIT NONE
 
 
 
! Local Variables
! ===============
 
      integer                                        ::  ikl,ikv   ,isn   ,isnMAX      !
      integer                                        ::  Mobilm,Mobiln                 !
 
      real(kind=real8)                               ::  DendOK                        ! Dendricity Switch
      real(kind=real8)                               ::  SaltOK                        ! Saltation  Switch
      real(kind=real8)                               ::  MeltOK                        ! Saltation  Switch (Melting Snow)
      real(kind=real8)                               ::  SnowOK                        ! Pack Top   Switch
      real(kind=real8)                               ::  SaltM1,SaltM2,SaltMo          ! Saltation  Parameters
      real(kind=real8)                               ::  SaltMx = -5.83e-2             !
      real(kind=real8)                               ::  ShearX                        ! Arg. Max Shear Stress
      real(kind=real8)                               ::  SaltSU,Salt_U                 !
      real(kind=real8)                               ::  ArgFac,Fac_Mo                 !
      real(kind=real8)                               ::  FacRBS =  2.868               !
      real(kind=real8)                               ::  FacTBS =  0.085               !
      real(kind=real8)                               ::  ArguSi                        !
      real(kind=real8)                               ::  hdrift =  1.00e+1             ! Inverse erodibl.Snow Lay.Thickn.
      real(kind=real8)                               ::  h_mmWE =  0.01e00             ! Eroded Snow Layer Min Thickness
!     real(kind=real8)                               ::  tfv_vk =  5.10e-1             ! * Fall Veloc. / Von Karman Cst
                                                                                       ! tfv (Terminal Fall Veloc. =.216)
                                                                                       ! /vk (Von Karman Constant  =.4  )
                                                                                       ! (Wamser & Lykosov,   1995
                                                                                       !  Contr.Atm.Phys. 68, p.90)
      real(kind=real8)                               ::  dzweqo,dzweqn,bsno_x          !
      real(kind=real8)                               ::                hsno_x          !
      real(kind=real8)                               ::  ro_new                        !
! #BM real(kind=real8)                               ::  PorSno,PorRef                 !
! #BS real(kind=real8)                               ::  Salt_f                        !
      real(kind=real8)                               ::  MIN_Mo                        ! Minimum Mobility Fresh Fallen *
! #MA real(kind=real8)                               ::  AgBlow = 1.00                 ! Snow Mobility    Time  Scale
                                                                                       ! 1 Day (F.Domine, pers.communic.)
! #BS real(kind=real8)                               ::  snofOK                        ! Threshd Snow Fall
 
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0 real(kind=real8)                               ::  Sno0WE,Sno1WE                 ! Snow Mass before/after Erosion
! #b0 real(kind=real8)                               ::  SnodWE                        ! Snow Mass              Erosion
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
! ALLOCATION                                              !
! ==========                                              !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)               THEN !
 
      allocate            ( Mobile(kcolp,mwp) )           !
      allocate            ( SaltSI(kcolp,mwp,nsnow) )     ! Snow Drift Index
      allocate            ( sdrift(kcolp,mwp,nsnow) )     !
      allocate            ( xdrift(kcolp,mwp) )           !
      allocate            ( zdrift(kcolp,mwp) )           !
      allocate            ( tdepos(kcolp,mwp) )           !
      allocate            ( zdepos(kcolp,mwp,nsnow) )     !
      allocate            ( dbsaux(kcolp,mwp) )           ! Drift Amount   (Dummy Variable)
      allocate            ( isagr1(kcolp,mwp) )           ! 1st     Layer History
      allocate            ( isagr2(kcolp,mwp) )           ! 2nd     Layer History
 
      allocate            ( WEagre(kcolp,mwp) )           ! Snow Water Equivalent Thickness
      allocate            ( Agrege(kcolp,mwp) )           ! 1. when Agregation constrained
      allocate            ( dzagr1(kcolp,mwp) )           ! 1st     Layer Thickness
      allocate            ( dzagr2(kcolp,mwp) )           ! 2nd     Layer Thickness
      allocate            ( T_agr1(kcolp,mwp) )           ! 1st     Layer Temperature
      allocate            ( T_agr2(kcolp,mwp) )           ! 2nd     Layer Temperature
      allocate            ( roagr1(kcolp,mwp) )           ! 1st     Layer Density
      allocate            ( roagr2(kcolp,mwp) )           ! 2nd     Layer Density
      allocate            ( etagr1(kcolp,mwp) )           ! 1st     Layer Water Content
      allocate            ( etagr2(kcolp,mwp) )           ! 2nd     Layer Water Content
      allocate            ( G1agr1(kcolp,mwp) )           ! 1st     Layer Dendricity/Spher.
      allocate            ( G1agr2(kcolp,mwp) )           ! 2nd     Layer Dendricity/Spher.
      allocate            ( G2agr1(kcolp,mwp) )           ! 1st     Layer Sphericity/Size
      allocate            ( G2agr2(kcolp,mwp) )           ! 2nd     Layer Sphericity/Size
      allocate            ( agagr1(kcolp,mwp) )           ! 1st     Layer Age
      allocate            ( agagr2(kcolp,mwp) )           ! 2nd     Layer Age
 
      END IF                                              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! DATA
! ====
 
! Initialization
! ==============
 
      IF (.NOT.BlowIn)                                             THEN
               BlowIn = .true.
               FacSBS =  1.             / FacRBS
               FacUBS =  1.             / FacTBS
               Por_BS =  1.             - BSnoRo/      rhoIce
               SheaBS =                   Por_BS/(1.00-Por_BS)
!              SheaBS =  Arg(sqrt(shear = max shear stress in snow)):
!              shear  =  3.420d00 * exp(-(Por_BS      +Por_BS)          &
!    &                                  /(1.00        -Por_BS))
!              SheaBS :  see de Montmollin         (1978),
!                        These Univ. Sci. Medic. Grenoble, Fig. 1 p. 124
 
             DO ikl=1,kcolp             ! Parameterization of u*th
             DO ikv=1,mwp
               rCd10n      =  1./  26.5 ! was developed from observations made
             END DO
             END DO                     ! during assumed neutral conditions
 
               write(6,5000)  1./  rCd10n
 5000          format(/,' Blowing Snow Model  Initialization     ',     &
     &                /,' Vt / u*t =',f8.2,' (Neutral Assumption)',     &
     &                /,'           ', 8x ,' (Budd assumes  26.5)',/)
      END IF
 
 
! Snow Age (Influence on Snow Erosion Threshold)
! ==============================================
 
! #BS DO isn=1,nsnow
! #BS DO ikl=1,kcolp
! #BS DO ikv=1,mwp
! #BS   agsnSV(ikl,ikv,isn) = agsnSV(ikl,ikv,isn) + dt__SV/86400.
! #BS END DO
! #BS END DO
! #BS END DO
! #BS DO ikl=1,kcolp
! #BS DO ikv=1,mwp
! #BS   isn    = max(1 ,        isnoSV(ikl,ikv))
! #BS   snofOK = max(0.,sign(1.,dsn_SV(ikl,ikv)-eps6))  !  Threshold=1.e-6
! #BS   agsnSV(ikl,ikv,isn) =   (1.-snofOK) *agsnSV(ikl,ikv,isn)! ~0.1 mm w.e./day
! #BS END DO
! #BS END DO
      IF (.NOT.BloMod)                                     GO TO 1000
! #AG STOP '?!&~@|@[#@#] --- INCONSISTANT SNOW AGE --- EMERGENCY STOP'
 1000 CONTINUE
 
 
! EROSION
! =======
 
      DO isn = 1,nsnow
      DO ikl = 1,kcolp
      DO ikv=1,mwp
 
! Below the high Snow Density Threshold  (ro__SV < BSnoRo)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DendOK   =  max(zer0,sign(un_1,eps6-G1snSV(ikl,ikv,isn)  ))  !
        SaltOK   =  min(1   , max(istdSV(2)-istoSV(ikl,ikv,isn),0))  !
        MeltOK   =     (un_1                                    &!
     &             -max(zer0,sign(un_1,Tf_Sno-eps6              &!
     &                                     -TsisSV(ikl,ikv,isn)  )))&! Melting Snow
     &           *  min(un_1,DendOK                             &!
     &                  +(1.-DendOK)                            &!
     &                      *sign(un_1,     G2snSV(ikl,ikv,isn)-1.0))! 1.0 for 1mm
        SnowOK   =  min(1   , max(isnoSV(ikl,ikv)      +1 -isn ,0))  ! Snow Switch
 
        G1snSV(ikl,ikv,isn) =      SnowOK *    G1snSV(ikl,ikv,isn)  &
     &                  + (1.- SnowOK)*min(G1snSV(ikl,ikv,isn),G1_dSV)
        G2snSV(ikl,ikv,isn) =      SnowOK *    G2snSV(ikl,ikv,isn)  &
     &                  + (1.- SnowOK)*min(G2snSV(ikl,ikv,isn),G1_dSV)
 
        SaltOK   =  min(un_1 , SaltOK +    MeltOK)       * SnowOK
        SaltM1   = -0.750e-2 * G1snSV(ikl,ikv,isn)              &
     &             -0.500e-2 * G2snSV(ikl,ikv,isn)+ 0.500e00
!       SaltM1   :  Guyomarc'h & Merindol, 1997, Ann. Glac.
!         CAUTION:  Guyomarc'h & Merindol Dendricity Sign is +
!         ^^^^^^^^                    MAR Dendricity Sign is -
        SaltM2   = -0.833d-2 * G1snSV(ikl,ikv,isn)              &
     &             -0.583d-2 * G2snSV(ikl,ikv,isn)+ 0.833d00
        SaltMo   = (DendOK   * SaltM1 + (1.-DendOK) *     SaltM2       )
 
! Increased Mobility of Deposed (blown) Snow (Mann et al., 2000, JGR 105,
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Fig.2 p.24496 & text below)
        MIN_Mo   =  0.
! #MA   MIN_Mo   =  0.6 * exp(-agsnSV(ikl,ikv,isn)                  /AgBlow)
        SaltMo   =                                    max(SaltMo,MIN_Mo)
 
        SaltMo   =  SaltOK   * SaltMo + (1.-SaltOK) * min(SaltMo,SaltMx)
!       SaltMo   =  SaltOK   * SaltMo - (1.-SaltOK) *     0.9500 ! Tuning
        SaltMo   =         max(SaltMo ,  eps6-un_1)
 
        SaltSU   =     (1.00d0+SaltMo)     *FacSBS
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0   Salt_U   =        -log(SaltSU)     *FacUBS
! #b0   IF (ii__AP(ikl).EQ.iwr_SV.AND.jj__AP(ikl).EQ.jwr_SV        .AND. &
! #b0&      ikv        .EQ.nwr_SV.AND.isn        .EQ.isnoSV(ikl,ikv))    &
! #b0&    write(6,6010) isnoSV(ikl,ikv),G1snSV(ikl,ikv,isn)/G1_dSV       &
! #b0&                             ,G2snSV(ikl,ikv,isn)/G1_dSV           &
! #b0&                             ,ro__SV(ikl,ikv,isn),agsnSV(ikl,ikv,isn)  &
! #b0&                             ,SaltM1, SaltM2, SaltMo, Salt_U       &
! #b0&                             ,us__SV(ikl,ikv)   / rCd10n
 6010     format(/,'SISVAT_BSn',6x                                       &
     &           ,6x,i3,2x,'G1         =',f6.3,'   G2         =',f7.3    &
     &           ,      '   ro [kg/m3] =',f9.3,'   Age* [Day] =',f9.3    &
     &           ,   /,27x,'SaltM1     =',f6.3,'   SaltM2     =',f7.3    &
     &           ,      '   Mobility I.=',f9.3,'   Vt   [m/s] =',f9.3    &
     &           ,   /,27x,'            ', 6x ,'               ', 7x     &
     &           ,      '               ', 9x ,'   Vn10 [m/s] =',f9.3)
 
! Above the high Snow Density Threshold  (ro__SV > BSnoRo)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Por_BS      =  1.000       - ro__SV(ikl,ikv,isn)     /rhoIce
        ShearX      =                Por_BS/max(eps6,un_1-Por_BS)
!       ShearX ==> Arg(sqrt(shear)) with shear = max shear stress in snow:
!       shear       =  3.420d00 * exp(-(Por_BS      +Por_BS)             &
!    &                                /max(eps6,un_1-Por_BS))
!                      see de Montmollin         (1978),
!                      These Univ. Sci. Medic. Grenoble, Fig. 1 p. 124
 
! Influence of Density on Shear Stress if ro__SV > BSnoRo
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ArgFac      =  max(zer0  ,SheaBS-ShearX)     !
!       Fac_Mo      =  exp(       ArgFac       )     ! ** NOT ** tuned
        Fac_Mo   =     exp(       ArgFac       )     ! = 1 if ro__SV < BSnoRo
                                                     ! < 1 if ro__SV > BSnoRo
! Snow Drift Index
! ~~~~~~~~~~~~~~~~
        SaltSU      =  max(eps6  ,    SaltSU)
        SaltSU      =  exp(Fac_Mo*log(SaltSU))
        ArguSi      =     -FacTBS              *us__SV(ikl,ikv)/rCd10n
        SaltSI(ikl,ikv,isn) = (SaltSU-exp(ArguSi)) *FacRBS
!       SaltSI          :  Generalization of the Snow Drift Index of
!                          Guyomarc'h & Merindol (1997, Ann.Glaciol.)
 
! Threshold Friction Velocity
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SnowOK   =  1 -min(1,iabs(isn-isnoSV(ikl,ikv)))
        Salt_U   =               -log(SaltSU)  *FacUBS
!       Salt_U   :  Guyomarc'h & Merindol, 1997, Ann. Glac.
 
        usthSV(ikl,ikv) =     SnowOK *   (Salt_U   *rCd10n)              &
     &              + (1.-SnowOK)*    usthSV(ikl,ikv)
 
! Threshold Friction Velocity (Budd et al., 1966)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #BA   usthSV(ikl,ikv) =     SnowOK *   (Salt_U   /26.5)                &
! #BA&              + (1.-SnowOK)*    usthSV(ikl,ikv)
!       Us(U10)     :  Budd et al.            1966, Ant.Res.Ser.9
!                 (see Pomeroy & Gray 1995 NHRI Sci.Rep.7(30)p.62)
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! OUTPUT           for Snow Erosion Variables
! #b0   IF (ii__AP(ikl).EQ.iwr_SV.AND.jj__AP(ikl).EQ.jwr_SV        .AND. &
! #b0&      ikv        .EQ.nwr_SV.AND.isn        .EQ.isnoSV(ikl,ikv))    &
! #b0&    write(6,6011)     Fac_Mo,Por_BS,SaltSI(ikl,ikv,isn),usthSV(ikl,ikv)
 6011     format(      27x,'Fac_Mo     =',f6.3,'   Por_BS     =',f7.3    &
     &           ,      '   Drift    I.=',f9.3,'   ut*_0[m/s] =',f9.3)
      END DO
      END DO
      END DO
 
 
! Deepest Mobile Snow Layer
! -------------------------
 
      DO ikl = 1,kcolp
      DO ikv=1,mwp
        Mobile(ikl,ikv) = nsnow+1
      END DO
      END DO
      DO isn =   nsnow,1,-1
      DO ikl = 1,kcolp
      DO ikv=1,mwp
        isnMAX      =       max(   1,          isnoSV(ikl,ikv)             )
        isnMAX      =       min( isn,          isnMAX                  )
        Mobiln      = isn * max(zer0,sign(un_1,SaltSI(ikl,ikv,isnMAX)))
        Mobilm      =   1 - min(1   ,          Mobile(ikl,ikv) -1 -Mobiln)
!       Mobilm      =   1   ONLY IF   Mobiln = Mobile(ikl) -1 (0 otherwise)
 
        Mobile(ikl,ikv) =                 Mobilm * Mobiln               &
     &              +              (1-Mobilm)* Mobile(ikl,ikv)
      END DO
      END DO
      END DO
 
 
! Weighting the Amount of Snow to erode
! -------------------------------------
 
      DO ikl = 1,kcolp
      DO ikv=1,mwp
        zdrift(ikl,ikv) = 0.0
        xdrift(ikl,ikv) = 0.0
        dbsaux(ikl,ikv) = dbs_SV(ikl,ikv)
      END DO
      END DO
 
      DO isn = 1,nsnow
      DO ikl = 1,kcolp
      DO ikv=1,mwp
        zdrift(ikl,ikv)     =  zdrift(ikl,ikv)                          &
     &            + 0.50 * dzsnSV(ikl,ikv,isn) * (3.25  -SaltSI(ikl,ikv,isn))
        sdrift(ikl,ikv,isn) =  SaltSI(ikl,ikv,isn)                      &
     &          *exp(  max(Ea_Min, -zdrift(ikl,ikv)     *hdrift     ))  &
     &          *min(1,max(0     ,  isn +1          -Mobile(ikl,ikv)))  &
     &          *min(1,max(0     ,  isnoSV(ikl,ikv)     -isn +1     ))  &
!                Last 2 Lines force sdrift = 0 outside mobile Snow Layers
     &          *      max(zer0, sign(un_1,         -dbs_SV(ikl,ikv)))
!                Erosion is allowed only if available Blowing Snow
        xdrift(ikl,ikv)     =           sdrift(ikl,ikv,isn) +xdrift(ikl,ikv)
        zdrift(ikl,ikv)     =  zdrift(ikl,ikv)                          &
     &            + 0.50 * dzsnSV(ikl,ikv,isn) * (3.25  -SaltSI(ikl,ikv,isn))
      END DO
      END DO
      END DO
 
! Normalization
! ~~~~~~~~~~~~~
      DO isn = 1,nsnow
      DO ikl = 1,kcolp
      DO ikv=1,mwp
        sdrift(ikl,ikv,isn) =  sdrift(ikl,ikv,isn) /max(eps6,xdrift(ikl,ikv))
      END DO
      END DO
      END DO
 
 
! Weighting the Amount of Snow to depose
! --------------------------------------
 
      DO ikl = 1,kcolp
      DO ikv=1,mwp
        zdrift(ikl,ikv) = 0.0
        tdepos(ikl,ikv) = 0.0
      END DO
      END DO
 
      DO isn = 1,nsnow
      DO ikl = 1,kcolp
      DO ikv=1,mwp
        zdepos(ikl,ikv,isn) =      exp(-zdrift(ikl,ikv)   )             &
     &          *min(1,max(0     ,  isn +1          -Mobile(ikl,ikv)))  &
     &          *min(1,max(0     ,  isnoSV(ikl,ikv    ) -isn +1     ))
!                Last 2 Lines force zdepos = 0 outside mobile Snow Layers
        tdepos(ikl,ikv) = tdepos(ikl,ikv) + zdepos(ikl,ikv,isn)
        zdrift(ikl,ikv) = zdrift(ikl,ikv) + dzsnSV(ikl,ikv,isn) *ro__SV(ikl,ikv,isn)&
     &                                              /rhoWat
      END DO
      END DO
      END DO
 
! Normalization
! ~~~~~~~~~~~~~
      DO isn = 1,nsnow
      DO ikl = 1,kcolp
      DO ikv=1,mwp
        zdepos(ikl,ikv,isn) = zdepos(ikl,ikv,isn) / max(eps6,tdepos(ikl,ikv))
      END DO
      END DO
      END DO
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0 DO ikl = 1,kcolp
! #b0 DO ikv=1,mwp
! #b0   IF (ii__AP(ikl).EQ.iwr_SV.AND.jj__AP(ikl).EQ.jwr_SV .AND.      &
! #b0&      ikv        .EQ.nwr_SV                          )        THEN
! #b0     Sno0WE =   0.
! #b0   DO isn=1,nsnow
! #b0     Sno0WE =   Sno0WE                                            &
! #b0&           +   dzsnSV(ikl,ikv,isn) *ro__SV(ikl,ikv,isn)
! #b0   END DO
! #b0     write(6,6005)   Sno0WE                    ,dbs_SV(ikl,ikv)
 6005     format(                                                      &
     &      18x,'MB0',6x,'Sno1WE [mm]=',f9.3,19x,'0  dbs_SV [mm]=',f9.6)
! #b0     SnodWE =   dbs_SV(ikl,ikv)
! #b0   END IF
! #b0 END DO
! #b0 END DO
 
 
! Weighted  Erosion (Erosion amount is distributed       ! dbs_SV decreases
! -----------------  over the upper Snow Pack)           ! dzsnSV decreases
 
      DO isn = 1,nsnow
      DO ikl = 1,kcolp
      DO ikv=1,mwp
        SnowOK      = min(1,max(isnoSV(ikl,ikv)+1-isn ,0))  ! Snow Switch
        dzweqo      = dzsnSV(ikl,ikv,isn) *ro__SV(ikl,ikv,isn)  ! [kg/m2, mm w.e.]
        bsno_x      = dbsaux(ikl,ikv)     *sdrift(ikl,ikv,isn)
        dzweqn      = dzweqo          +bsno_x
        dzweqn  = max(dzweqn,          h_mmWE *SnowOK)
        dbs_SV(ikl,ikv) = dbs_SV(ikl,ikv)    +(dzweqo -dzweqn)
        dzsnSV(ikl,ikv,isn) =              dzweqn                       &
     &                       /max(eps6,ro__SV(ikl,ikv,isn))
      END DO
      END DO
      END DO
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0 DO ikl = 1,kcolp
! #b0 DO ikv=1,mwp
! #b0   IF (ii__AP(ikl)    .EQ. 1 .AND.  jj__AP(ikl)     .EQ. 1)    THEN
! #b0     SnodWE =   SnodWE         -dbs_SV(ikl,ikv)
! #b0     Sno1WE =   0.
! #b0   DO isn=1,nsnow
! #b0     Sno1WE =   Sno1WE                                             &
! #b0&           +   dzsnSV(ikl,ikv,isn)*ro__SV(ikl,ikv,isn)
! #b0   END DO
! #b0     write(6,6006)Sno1WE      , dbs_SV(ikl,ikv)
 6006     format(                                                       &
     &      18x,'MB1',6x,'Sno1WE [mm]=',f9.3,19x,'1  dbs_SV [mm]=',f9.6)
! #b0     write(6,6007)Sno1WE    ,SnodWE   ,Sno0WE,                     &
! #b0&                (Sno1WE    -SnodWE   -Sno0WE)
 6007     format(                                                       &
     &      18x,'MB ',5x,'(After  [mm]=',f6.0, ')-(Erosion[mm]=', f7.3, &
     &                                         ')-(Before [mm]=', f9.3, &
     &                                         ')= Budget [mm]=', f9.6)
! #b0   END IF
! #b0 END DO
! #b0 END DO
 
 
! ACCUMULATION of BLOWN SNOW                             ! dsn_SV decreases
! --------------------------                             ! dzsnSV increases
 
        DO ikl = 1,kcolp
        DO ikv=1,mwp
          tdepos(ikl,ikv) = dsn_SV(ikl,ikv) * dsnbSV(ikl,ikv) * dt__SV
          WEagre(ikl,ikv) = 0.
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0     IF (ii__AP(ikl).EQ.iwr_SV.AND.jj__AP(ikl).EQ.jwr_SV      .AND. &
! #b0&        ikv        .EQ.nwr_SV.AND.0          .LT.isnoSV(ikl,ikv))  &
! #b0&        write(6,6003) tdepos(ikl,ikv)  ,Mobile(ikl,ikv)
 6003         format(/,41x,'tdepos [-] =',f6.3,40x,'Mobil',i3            &
     &              ,/,27x,'Salt.Index    sdrift'                        &
     &              ,      '    zdepos  ro__snow  ro_bsnow  roN_snow'    &
     &              ,                '  dz__snow  dz_bsnow  dzN_snow'    &
     &              ,                '  d___snow'                        &
     &              ,/,27x,'             [kg/m3]   [kg/m3]   [kg/m3]'    &
     &              ,                '       [m]       [m]       [m]'    &
     &              ,                '   [kg/m2]')
        END DO
        END DO
 
      DO isn =     nsnow,1,-1
        DO ikl = 1,kcolp
        DO ikv=1,mwp
          WEagre(ikl,ikv) = WEagre(ikl,ikv) + ro__SV(ikl,ikv,isn)*dzsnSV(ikl,ikv,isn)
          isagr1(ikl,ikv) = istoSV(ikl,ikv,isn)
          isagr2(ikl,ikv) = 0.
 
! Density of deposited blown Snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ro_new =                    BSnoRo
 
! Density of deposited blown Snow (de Montmollin, 1978)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #BM     PorSno =      1.0d00     -  ro__SV(ikl,ikv,isn)               &
! #BM&                             /  rhoIce
! #BM     Salt_f =      usthSV(ikl,ikv)/  max(eps6,   us__SV(ikl,ikv))
! #BM     Salt_f =  min(Salt_f     ,  un_1)
! #BM     PorRef =      PorSno     /  max(eps6,1.-PorSno)               &
! #BM&             +log(Salt_f)
! #BM     Por_BS =      PorRef     / (1.0d00 + PorRef)
! #BM     ro_new =      rhoIce     * (1.0d00 - Por_BS)
! #BM     ro_new =  max(ro_new     ,  BSnoRo)
 
          roagr1(ikl,ikv) = ro__SV(ikl,ikv,isn)
          roagr2(ikl,ikv) = ro_new
          hsno_x      = tdepos(ikl,ikv)*  zdepos(ikl,ikv,isn)
 
          dzagr1(ikl,ikv) = dzsnSV(ikl,ikv,isn)
          dzagr2(ikl,ikv) = hsno_x     /  ro_new
!         Conversion    [kg/m2, i.e., mm w.e.] -----> [mSnow]
 
          dsn_SV(ikl,ikv) = dsn_SV(ikl,ikv)-  hsno_x / dt__SV
 
! Other Snow Properties
! ~~~~~~~~~~~~~~~~~~~~~
          T_agr1(ikl,ikv) =    TsisSV(ikl,ikv,isn)
          T_agr2(ikl,ikv) =min(Tf_Sno,TaT_SV(ikl,ikv))
          etagr1(ikl,ikv) =    eta_SV(ikl,ikv,isn)
          etagr2(ikl,ikv) =    0.0
          G1agr1(ikl,ikv) =    G1snSV(ikl,ikv,isn)
          G1agr2(ikl,ikv) =    G1_dSV
          G2agr1(ikl,ikv) =    G2snSV(ikl,ikv,isn)
          G2agr2(ikl,ikv) =    ADSdSV
! #BY     G2agr2(ikl,ikv) =    0.87d0
!         Budd et al. 1966, 2~m Average /Table 5 p. 97
 
          agagr1(ikl,ikv) =    agsnSV(ikl,ikv,isn)
          agagr2(ikl,ikv) =    0.
          Agrege(ikl,ikv) =    1.
        END DO
        END DO
 
! Agregation
! ~~~~~~~~~~
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
!            **********
        call SISVAT_zAg                                                &
     &             (ikl,ikv,isagr1(ikl,ikv),isagr2(ikl,ikv),WEagre(ikl,ikv)&
     &                 ,dzagr1(ikl,ikv),dzagr2(ikl,ikv),T_agr1(ikl,ikv),T_agr2(ikl,ikv)&
     &                 ,roagr1(ikl,ikv),roagr2(ikl,ikv),etagr1(ikl,ikv),etagr2(ikl,ikv)&
     &                 ,G1agr1(ikl,ikv),G1agr2(ikl,ikv),G2agr1(ikl,ikv),G2agr2(ikl,ikv)&
     &                 ,agagr1(ikl,ikv),agagr2(ikl,ikv),Agrege(ikl,ikv)&
     &                 )
!            **********
 
        END DO
        END DO
 
 
        DO ikl = 1,kcolp
        DO ikv=1,mwp
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! OUTPUT           for Snow Erosion Variables
! #b0     IF (ii__AP(ikl).EQ.iwr_SV.AND.jj__AP(ikl).EQ.jwr_SV      .AND. &
! #b0&        ikv        .EQ.nwr_SV.AND.isn        .LE.isnoSV(ikl,ikv))  &
! #b0&        write(6,6004)   isn          ,SaltSI(ikl,ikv,isn)          &
! #b0&                     ,sdrift(ikl,ikv,isn),zdepos(ikl,ikv,isn)      &
! #b0&                     ,ro__SV(ikl,ikv,isn),roagr2(ikl,ikv),roagr1(ikl,ikv)  &
! #b0&                     ,dzsnSV(ikl,ikv,isn),dzagr2(ikl,ikv),dzagr1(ikl,ikv)  &
! #b0&                     ,dsn_SV(ikl,ikv)
 6004         format((27x,i3,f7.2,2f10.6,3f10.3,4f10.6))
 
          istoSV(ikl,ikv,isn) = isagr1(ikl,ikv)
          dzsnSV(ikl,ikv,isn) = dzagr1(ikl,ikv)
          TsisSV(ikl,ikv,isn) = T_agr1(ikl,ikv)
          ro__SV(ikl,ikv,isn) = roagr1(ikl,ikv)
          eta_SV(ikl,ikv,isn) = etagr1(ikl,ikv)
          G1snSV(ikl,ikv,isn) = G1agr1(ikl,ikv)
          G2snSV(ikl,ikv,isn) = G2agr1(ikl,ikv)
          agsnSV(ikl,ikv,isn) = agagr1(ikl,ikv)
 
        END DO
        END DO
 
      END DO
 
! OUTPUT in SISVAT for ikl = 1 (preferably for Stand Alone Version)
! OUTPUT           for SnowFall and Snow Buffer
! #s2   IF          (isnoSV(1,1) .GT. 0)                               &
! #s2&  write(6,6008)isnoSV(1,1),  dsn_SV(1) *dt__SV + BufsSV(1,1),    &
! #s2&              (dzsnSV(1,isn)*ro__SV(1,isn),isn=1,isnoSV(1,1))
 6008   format(i3,'  dsn+Buf=',f6.2,6x,'A dz *ro =',10f6.2,            &
     &                                       (/,35x,10f6.2))
 
        DO ikl = 1,kcolp
        DO ikv=1,mwp
          hdrift      =  tdepos(ikl,ikv)/dt__SV
          esnbSV(ikl,ikv) = (dsnbSV(ikl,ikv)-1.00)*hdrift/max(dsn_SV(ikl,ikv),eps6) &
     &                  +dsnbSV(ikl,ikv)
          dsnbSV(ikl,ikv) =          min(un_1,   max(zer0,esnbSV(ikl,ikv) )   )
!         dsnbSV is now the Blown Snow fraction of precipitating snow
!                will be used for characterizing the Buffer Layer
!               (see update of  Bros_N, G1same, G2same, zroOLD, zroNEW)
        END DO
        END DO
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
! DE-ALLOCATION                                           !
! =============                                           !
 
      IF (FlagDALLOC)                                THEN !
 
      deallocate          ( SaltSI )                      ! Snow Drift Index
      deallocate          ( sdrift )                      !
      deallocate          ( xdrift )                      !
      deallocate          ( zdrift )                      !
      deallocate          ( tdepos )                      !
      deallocate          ( zdepos )                      !
      deallocate          ( dbsaux )                      ! Drift Amount   (Dummy Variable)
      deallocate          ( isagr1 )                      ! 1st     Layer History
      deallocate          ( isagr2 )                      ! 2nd     Layer History
 
      deallocate          ( WEagre )                      ! Snow Water Equivalent Thickness
      deallocate          ( Agrege )                      ! 1. when Agregation constrained
      deallocate          ( dzagr1 )                      ! 1st     Layer Thickness
      deallocate          ( dzagr2 )                      ! 2nd     Layer Thickness
      deallocate          ( T_agr1 )                      ! 1st     Layer Temperature
      deallocate          ( T_agr2 )                      ! 2nd     Layer Temperature
      deallocate          ( roagr1 )                      ! 1st     Layer Density
      deallocate          ( roagr2 )                      ! 2nd     Layer Density
      deallocate          ( etagr1 )                      ! 1st     Layer Water Content
      deallocate          ( etagr2 )                      ! 2nd     Layer Water Content
      deallocate          ( G1agr1 )                      ! 1st     Layer Dendricity/Spher.
      deallocate          ( G1agr2 )                      ! 2nd     Layer Dendricity/Spher.
      deallocate          ( G2agr1 )                      ! 1st     Layer Sphericity/Size
      deallocate          ( G2agr2 )                      ! 2nd     Layer Sphericity/Size
      deallocate          ( agagr1 )                      ! 1st     Layer Age
      deallocate          ( agagr2 )                      ! 2nd     Layer Age
 
      END IF                                              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SISVAT_BSn
 
 
 
      subroutine SISVAT_BDu
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_BDu                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_BDu treats Dust Erosion                            |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Thu 14-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     OUTPUT:  usthSV   : Blowing Snow Erosion   Threshold           [m/s] |
!     ^^^^^^                                                               |
!                                                                          |
!     REFER. : Fecan, F., B. Marticorena and G. Bergametti, 1999   (Fal99) |
!     ^^^^^^^^ Ann. Geophysicae 17, 149--157                               |
!              u* threshold: adapted from Fig. 4 p. 153                    |
!              Clay Content:         from Tab. 2 p. 155                    |
!                                                                          |
!--------------------------------------------------------------------------+
 
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_dat
      use Mod_SISVAT_kkl
      use Mod_SISVAT_BDu
 
 
 
      IMPLICIT NONE
 
 
 
! Local   Variables
! =================
 
      integer           ::             ikl,ikv   ,  isot
 
      real(kind=real8)  ::             eta_Du,usthDu
 
 
! Initialisation
! ==============
 
      IF (.NOT.logust)                                              THEN
        DO isot=1,nsot
               etaust(isot) = 0.0014 * claypc(isot) * claypc(isot)      &! Fal99
     &                      + 0.17   * claypc(isot)                      ! Eqn.(14)
        END DO                                                           !  p. 154
               logust = .true.
      END IF
 
 
! Soil Erodibility
! ----------------
 
      DO ikl = 1,kcolp
      DO ikv = 1,mwp
        eta_Du      =  max(     eta_SV(ikl,ikv,0),etaust(isotSV(ikl,ikv)))   ! Fal99
        eta_Du      =  max(eps6,eta_SV(ikl,ikv,0)-eta_Du             )   ! Eqn.(15)
        usthDu      = sqrt(un_1+1.21*exp(0.68*   log(eta_Du)    ))      &! p.  155
     &              *                         ustdmn(isotSV(ikl,ikv))   &!
     &              *                         f__ust(ivgtSV(ikl,ikv))    !
        usthSV(ikl,ikv) =                                               &
     &                 usthSV(ikl,ikv)*(1-max(0,1-isnoSV(ikl,ikv))) +   &
     &                 usthDu     *   max(0,1-isnoSV(ikl,ikv))
      END DO
      END DO
 
 
      return
      end subroutine SISVAT_BDu
 
 
 
      subroutine SISVAT_SIc(                                            &
! #m2&                      SIvAcr                                      &
     &                     )
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_SIc                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_SIc treats Sea-Ice and Ocean Latent Heat Exchanges |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Thu 14-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     INPUT:   TaT_SV   : SBL Top    Temperature                       [K] |
!     ^^^^^    isnoSV   : total Nb of Ice/Snow Layers                  [-] |
!              LSmask   : Land-Sea   Mask                              [-] |
!              dsn_SV   : Snow  Intensity                      [mm w.e./s] |
!                                                                          |
!     INPUT /  TsisSV   : Snow/Ice/Soil-Water Temperature              [K] |
!     OUTPUT:  eta_SV   : Soil/Snow Water   Content                [m3/m3] |
!     ^^^^^^   dzsnSV   : Snow Layer        Thickness                  [m] |
!                                                                          |
!     OUTPUT:  HFraSV   : Frazil            Thickness                  [m] |
!     ^^^^^^                                                               |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #SJ: Sea-Ice Bottom   accretion  and  ocean cooling due to SnowFall  |
!                                                                          |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_ctr
      use Mod_SISVAT_dat
      use Mod_SISVAT_kkl
      use Mod_SISVAT_SIc
! #m2 use Mod_SISVATLSIc
 
 
 
      IMPLICIT NONE
 
 
 
! Local Variables
! ===============
 
      integer           ::            ikl   ,ikv   ,n
      real(kind=real8)  ::            OCN_OK
! #SJ real(kind=real8)  ::            SIceOK
! #SJ real(kind=real8)  ::            SIcFrz
! #SJ real(kind=real8)  ::            Twat_n
 
      real(kind=real8)  ::            SalIce = 10.                     ! Sea-Ice   Salinity
      real(kind=real8)  ::            SalWat = 35.                     ! Sea-Water Salinity
                                                                       !  Typical Salinities in Terra Nova Bay
                                                                       ! (Bromwich and Kurtz,   1984, JGR, p.3568;
                                                                       !  Cavalieri and Martin, 1985,      p. 248)
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
! #m2 IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
! #m2 allocate                                       ( SIvAcr(kcolp,mwp) )
 
! #m2 END IF
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! Initialisation
! ==============
 
      IF (.NOT.SIcINI)                                              THEN
               SIcINI =  .true.
               Crodzw =  hC_Wat*rhoWat *          dz_dSV(0)     ! [J/m2/K]
               Lro__I =  LhfH2O*rhoIce *(1.-1.e-3*SalIce       &! [J/m3]
     &                 -(SalIce/SalWat)*(1.-1.e-3*SalWat) )     !
 
! OUTPUT/Verification: Energy/Water Budget
! #e1          Lro__I =  LhfH2O*rhoIce
 
      END IF
 
 
! Snow Fall cools Sea Water
! =========================
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
        OCN_OK        =  (1   -LSmask(ikl,ikv) )                       &! Free Ocean
     &             *max(0,1   -isnoSV(ikl,ikv)   )                      !
! #SJ   TsisSV(ikl,ikv,0) =        TsisSV(ikl,ikv,0)                   &! [K]
! #SJ& -OCN_OK*(Cn_dSV*(Tf_Sno-TaT_SV(ikl,ikv)   )                     &! [J/kg]
! #SJ&         +LhfH2O*(1.    -eta_SV(ikl,ikv,0)))                     &! [J/kg]
! #SJ&        * dsn_SV(ikl,ikv)   *dt__SV          / Crodzw             ! [kg/m2]
 
 
! Sea-Ice Formation
! =================
 
! #SJ   Twat_n      =      max(TsisSV(ikl,ikv,0  )  ,Tf_Sea)    ! [K]
! #SJ   SIcFrz      =  (Twat_n-TsisSV(ikl,ikv,0  ) )*Crodzw/Lro__I &! [m]
! #SJ&                              * 0.75
!  ***  Hibler (1984), Ocean Heat Flux: 25% of cooling (ANTARCTIC Ocean)
!      (Hansen and Takahashi Eds)
!       Geophys. Monogr. 29, M. Ewing Vol. 5, AGU, p. 241
 
 
! Frazil  Formation
! -----------------
 
! #SJ   HFraSV(ikl,ikv) =          SIcFrz           *OCN_OK
 
 
! Growth of the Sea-Ice First Ice Floe
! ------------------------------------
 
! #SJ   SIceOK        =  (1   -LSmask(ikl,ikvp,n)  )         &! Ice Cover.Ocean
! #SJ&             *min(  1   ,isnoSV(ikl,ikv)     )          !
! #SJ   dzsnSV(ikl,ikv,1) =        dzsnSV(ikl,ikv,1)         &! Vertical Acret.
! #SJ&                +        SIcFrz           *SIceOK       !
 
 
! OUTPUT/Verification: SeaIce Conservation: Diagnostic of Surface Mass Balance
! #m2   SIvAcr(ikl,ikv) = rhoIce*SIcFrz     *(OCN_OK+SIceOK)           &
! #m2&              - dt__SV*dsn_SV(ikl,ikv)* OCN_OK
 
 
! Water Fluxes Update
! -------------------
 
        RnofSV(ikl,ikv) =          RnofSV(ikl,ikv)                     &
     &              +          dsn_SV(ikl,ikv) *     OCN_OK
        dsn_SV(ikl,ikv) =          dsn_SV(ikl,ikv) * (1.-OCN_OK)
 
      END DO
      END DO
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
! #m2 IF (FlagDALLOC)                                             THEN !
 
! #m2 deallocate                                     ( SIvAcr )
 
! #m2 END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SISVAT_SIc
 
 
 
      subroutine SISVAT_zSn
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_zSn                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_zSn manages the Snow Pack vertical Discretization  |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Mon  4-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT /  NLaysv   = New             Snow Layer  Switch               |
!     OUTPUT:  isnoSV   = total Nb of Ice/Snow Layers                      |
!     ^^^^^^   ispiSV   = 0,...,nsno: Uppermost Superimposed Ice Layer     |
!              iiceSV   = total Nb of Ice      Layers                      |
!              istoSV   = 0,...,5 :   Snow     History (see istdSV data)   |
!                                                                          |
!     INPUT /  TsisSV   : Soil/Ice Temperatures (layers -nsoil,-nsoil+1, 0)|
!     OUTPUT:           & Snow     Temperatures (layers  1,2,...,nsno) [K] |
!     ^^^^^^   ro__SV   : Soil/Snow Volumic Mass                   [kg/m3] |
!              eta_SV   : Soil/Snow Water   Content                [m3/m3] |
!              dzsnSV   : Snow Layer        Thickness                  [m] |
!              G1snSV   : Dendricity (<0) or Sphericity (>0) of Snow Layer |
!              G2snSV   : Sphericity (>0) or Size            of Snow Layer |
!              agsnSV   : Snow       Age                             [day] |
!                                                                          |
!     METHOD:  1) Agregate the thinest Snow Layer                          |
!     ^^^^^^      if a new Snow Layer has been precipitated   (NLaysv = 1) |
!              2) Divide   a too thick Snow Layer except                   |
!                 if the maximum Number of Layer is reached                |
!                 in this case forces                          NLay_s = 1  |
!              3) Agregate the thinest Snow Layer                          |
!                 in order to divide a too thick Snow Layer                |
!                 at next Time Step when                       NLay_s = 1  |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: #SX: Search Ice/Snow Interface in Snow Model  |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_zSn.vz        | #vz: OUTPUT/Verification: Snow Layers Agrega. |
!                          |      unit 41, SubRoutine  SISVAT_zSn **ONLY** |
!   # SISVAT_GSn.vp        | #vp: OUTPUT/Verification: Snow   Properties   |
!                          |      unit 47, SubRoutines SISVAT_zSn, _GSn    |
!   # stdout               | #s1: OUTPUT of Snow Layers Agregation         |
!                          |      unit  6, SubRoutine  SISVAT_zSn, _zAg    |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
      use Mod_SISVAT_zSn
      use Mod_SISVATLzSn
 
 
 
      IMPLICIT NONE
 
 
 
! Internal Variables
! ==================
 
      integer                                      ::  ikl,ikv   ,isn   ,i !
 
      integer                                      ::  LstLay          ! 0 ====> isnoSV = 1
      integer                                      ::  isno_n          ! Snow Normal.Profile
      integer                                      ::  iice_n          ! Ice  Normal.Profile
      integer                                      ::  iiceOK          ! Ice         Switch
      integer                                      ::  icemix = 0      ! 0 ====> Agregated Snow+Ice=Snow
                                                                       ! 1                          Ice
      real(kind=real8)                             ::  staggr          !              stagger  Switch
 
      real(kind=real8)                             ::  OKthin          ! Swich ON  a  new thinest layer
      real(kind=real8)                             ::  dz_dif          ! difference from ideal discret.
      real(kind=real8)                             ::  thickL          ! Thick Layer          Indicator
! #SX real(kind=real8)                             ::  OK_ICE          ! Swich ON   uppermost Ice Layer
 
      real(kind=real8)                             ::  dzepsi = 0.0015 ! Min Single Snw Layer Thickness
      real(kind=real8)                             ::  dzxmin = 0.0020 ! Min Acceptable Layer Thickness
      real(kind=real8)                             ::  dz_min = 0.0050 ! Min Local      Layer Thickness
      real(kind=real8)                             ::  dz_max = 0.0300 ! Min Gener.     Layer Thickness
!     CAUTION:                           dz_max > dz_min*2 is required ! Otherwise re-agregation is
                                                                       ! activated  after splitting
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
      allocate            ( NLay_s(kcolp,mwp) )                        ! Split Snow Layer         Switch
      allocate            ( isagr1(kcolp,mwp) )                        ! 1st     Layer History
      allocate            ( isagr2(kcolp,mwp) )                        ! 2nd     Layer History
      allocate            ( isn1  (kcolp,mwp) )                        ! 1st layer to stagger
      allocate            ( WEagre(kcolp,mwp) )                        ! Snow Water Equivalent Thickness
      allocate            ( dzthin(kcolp,mwp) )                        ! Thickness of the thinest layer
      allocate            ( Agrege(kcolp,mwp) )                        ! 1. when Agregation constrained
      allocate            ( dzagr1(kcolp,mwp) )                        ! 1st     Layer Thickness
      allocate            ( dzagr2(kcolp,mwp) )                        ! 2nd     Layer Thickness
      allocate            ( T_agr1(kcolp,mwp) )                        ! 1st     Layer Temperature
      allocate            ( T_agr2(kcolp,mwp) )                        ! 2nd     Layer Temperature
      allocate            ( roagr1(kcolp,mwp) )                        ! 1st     Layer Density
      allocate            ( roagr2(kcolp,mwp) )                        ! 2nd     Layer Density
      allocate            ( etagr1(kcolp,mwp) )                        ! 1st     Layer Water Content
      allocate            ( etagr2(kcolp,mwp) )                        ! 2nd     Layer Water Content
      allocate            ( G1agr1(kcolp,mwp) )                        ! 1st     Layer Dendricity/Spher.
      allocate            ( G1agr2(kcolp,mwp) )                        ! 2nd     Layer Dendricity/Spher.
      allocate            ( G2agr1(kcolp,mwp) )                        ! 1st     Layer Sphericity/Size
      allocate            ( G2agr2(kcolp,mwp) )                        ! 2nd     Layer Sphericity/Size
      allocate            ( agagr1(kcolp,mwp) )                        ! 1st     Layer Age
      allocate            ( agagr2(kcolp,mwp) )                        ! 2nd     Layer Age
 
! #vz allocate            ( dz_ref(nsnow) )                            ! Snow Reference Discretization
! #vz allocate            ( dzwdif(nsnow) )                            !
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz IF (.NOT.as_opn)                                              THEN
! #vz          as_opn=.true.
! #vz     open(unit=41,status='unknown',file='SISVAT_zSn.vz')
! #vz     rewind    41
! #vz END IF
 
! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp IF (.NOT.VP_opn)                                              THEN
! #vp          VP_opn=.true.
! #vp     open(unit=47,status='unknown',file='SISVAT_GSn.vp')
! #vp     rewind    47
! #vp END IF
 
 
! Constrains Agregation         of too thin  Layers
! =================================================
 
! Search the thinest  non-zero Layer
! ----------------------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          dzthin(ikl,ikv) = 0.                          ! Arbitrary unrealistic
        END DO
        END DO                                          !       Layer Thickness
      DO   isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isno_n    =             isnoSV(ikl,ikv)-isn+1 ! Snow Normal.Profile
          iice_n    =             iiceSV(ikl,ikv)-isn   ! Ice  Normal.Profile
          iiceOK    = min(1,max(0,iice_n         +1))   ! Ice         Switch
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz     dz_ref(isn) =                                &!
! #vz&          dz_min *((1-iiceOK)*isno_n*isno_n      &! Theoretical Profile
! #vz&                 +    iiceOK *    2**iice_n)     &!
! #vz&               /max(1,isnoSV(ikl,ikv))            !
 
          dz_dif      = max(zer0,                      &! Actual      Profile
     &          dz_min                                 &!
     &                 *((1-iiceOK)*isno_n*isno_n      &! Theoretical Profile
     &                 +    iiceOK *2.   **iice_n)     &!
     &        - dzsnSV(ikl,ikv, isn)                    )   ! Actual      Profile
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz     dzwdif(isn) =     dz_dif                      !
 
          OKthin      = max(zer0,                      &!
     &                      sign(un_1,                 &!
     &                           dz_dif-dzthin(ikl,ikv)))  &! 1.=> New thinest Lay.
     &                * max(0,                         &! 1 => .le. isnoSV
     &                  min(1,                         &! 1 => isn is in the
     &                      isnoSV(ikl,ikv)-isn +1 ))  &!          Snow Pack
     &                * min(un_1,                      &         !
!
!                       1st additional Condition to accept OKthin
     &                  max(zer0,                               &! combination
     &                      sign(un_1,G1snSV(ikl,ikv,      isn  )   &! G1 with same
     &                               *G1snSV(ikl,ikv,max(1,isn-1))))&!  sign => OK
!
!                       2nd additional Condition to accept OKthin
     &                + max(zer0,                               &! G1>0
     &                      sign(un_1,G1snSV(ikl,ikv,      isn   )))&!  =>OK
!
!                       3rd additional Condition to accept OKthin
     &                + max(zer0,                               &! dz too small
     &                      sign(un_1,dzxmin                    &!  =>OK
     &                               -dzsnSV(ikl,ikv,      isn   ))))!
 
          i_thin(ikl,ikv) =    (1. - OKthin)  * i_thin(ikl,ikv)&! Update   thinest Lay.
     &                         + OKthin   * isn         !                Index
          dzthin(ikl,ikv) =    (1. - OKthin)  * dzthin(ikl,ikv)&!
     &                         + OKthin   * dz_dif      !
        END DO
        END DO
      END DO
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz     write(41,4150) daHost     ,n___SV(  lwriSV(1,1))             &
! #vz&                  ,i_thin(1,1),dzsnSV(1,i_thin(1,1))
 4150     format(/,'-',a18,i5,' ',70('-'),                             &
     &           /,' Thinest ',i3,':',f9.3)
 
      DO   isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
          OKthin =      max(zer0,                       &!
     &                      sign(un_1,                  &!
     &                           dz_min                 &!
     &                          -dzsnSV(ikl,ikv,isn)))  &!
     &                * max(zer0,                       &! ON if dz > 0
     &                      sign(un_1,                  &!
     &                           dzsnSV(ikl,ikv,isn)-eps6)) &!
     &           *min(1,max(0,                          &! Multiple Snow    Lay.
     &                      min (1,                     &! Switch = 1
     &                           isnoSV(ikl,ikv)        &!   if isno > iice + 1
     &                          -iiceSV(ikl,ikv)-1))    &!
                                                         !
     &             +int(max(zer0,                       &!
     &                      sign(un_1,                  &!
     &                           dzepsi                 &! Minimum accepted for
     &                          -dzsnSV(ikl,ikv,isn)))) &! 1 Snow Layer over Ice
     &             *int(max(zer0,                       &! ON if dz > 0
     &                      sign(un_1,                  &!
     &                           dzsnSV(ikl,ikv,isn)-eps6)))&!
     &                 *(1 -min (abs(isnoSV(ikl,ikv)    &! Switch = 1
     &                              -iiceSV(ikl,ikv)-1),1)) &!   if isno = iice + 1
                                                         !
     &                 +max(0,                          &! Ice
     &                      min (1,                     &! Switch
     &                           iiceSV(ikl,ikv)+1-isn)))   &!
     &             *min(un_1,                                    &!
     &                  max(zer0,                                &! combination
     &                      sign(un_1,G1snSV(ikl,ikv,      isn  )&! G1>0 + G1<0
     &                               *G1snSV(ikl,ikv,max(1,isn-1)))) &! NO
     &                + max(zer0,                                &!
     &                      sign(un_1,G1snSV(ikl,ikv,      isn   ))) &!
     &                + max(zer0,                                &!
     &                      sign(un_1,dzxmin                     &!
     &                               -dzsnSV(ikl,ikv,      isn   )))) !
          i_thin(ikl,ikv) =    (1. - OKthin)  * i_thin(ikl,ikv) &! Update   thinest Lay.
     &                         + OKthin   * isn          !                Index
        END DO
        END DO
      END DO
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz     write(41,4151) i_thin(1,1),dzsnSV(1,i_thin(1,1))             &
! #vz&                  ,isnoSV(1,1),dzsnSV(1,isnoSV(1,1))
 4151     format(' Thinest ',i3,':',f9.3,'   Max   =',i3,f12.3)
 
! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp   write(47,470)(G1snSV(1,isn),isn=1,isnoSV(1,1))
 470    format('Before _zCr1: G1 = ',10f8.1,(/,19x,10f8.1))
! #vp   write(47,472)(G2snSV(1,isn),isn=1,isnoSV(1,1))
 472    format('              G2 = ',10f8.1,(/,19x,10f8.1))
 
 
! Index of the contiguous Layer to agregate
! -----------------------------------------
 
!          **********
      call SISVAT_zCr
!          **********
 
 
! Assign the 2 Layers to agregate
! -------------------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isn         =    i_thin(ikl,ikv)
          isagr1(ikl,ikv) =    istoSV(ikl,ikv,isn)
          isagr2(ikl,ikv) =    istoSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          dzagr1(ikl,ikv) =    dzsnSV(ikl,ikv,isn)
          dzagr2(ikl,ikv) =    dzsnSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          T_agr1(ikl,ikv) =    TsisSV(ikl,ikv,isn)
          T_agr2(ikl,ikv) =    TsisSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          roagr1(ikl,ikv) =    ro__SV(ikl,ikv,isn)
          roagr2(ikl,ikv) =    ro__SV(ikl,ikv,isn+LIndsv(ikl,ikv))
          etagr1(ikl,ikv) =    eta_SV(ikl,ikv,isn)
          etagr2(ikl,ikv) =    eta_SV(ikl,ikv,isn+LIndsv(ikl,ikv))
          G1agr1(ikl,ikv) =    G1snSV(ikl,ikv,isn)
          G1agr2(ikl,ikv) =    G1snSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          G2agr1(ikl,ikv) =    G2snSV(ikl,ikv,isn)
          G2agr2(ikl,ikv) =    G2snSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          agagr1(ikl,ikv) =    agsnSV(ikl,ikv,isn)
          agagr2(ikl,ikv) =    agsnSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          LstLay      = min(1,max(  0,isnoSV(ikl,ikv) -1))  ! 0  if single Layer
          isnoSV(ikl,ikv) =               isnoSV(ikl,ikv)  &! decrement   isnoSV
     &     -(1-LstLay)* max(zer0,                      &! if downmost  Layer
     &                      sign(un_1,eps_21           &! <  1.e-21 m
     &                               -dzsnSV(ikl,ikv,1)))   !
          isnoSV(ikl,ikv) = max(   0,     isnoSV(ikl,ikv)   )   !
          Agrege(ikl,ikv) = max(zer0,                  &!
     &                      sign(un_1,dz_min           &! No Agregation
     &                               -dzagr1(ikl,ikv)  ))  &! if too thick Layer
     &                               *LstLay           &! if  a single Layer
     &                * min( max(0   ,isnoSV(ikl,ikv)+1&! if Agregation
     &                               -i_thin(ikl,ikv)  &!    with    a Layer
     &                               -LIndsv(ikl,ikv)  ),1) !    above the Pack
 
          WEagre(ikl,ikv) = 0.
        END DO
        END DO
 
        DO isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
          WEagre(ikl,ikv) = WEagre(ikl,ikv) + ro__SV(ikl,ikv,isn)*dzsnSV(ikl,ikv,isn)  &
     &                                *min(1,max(0,i_thin(ikl,ikv)+1-isn))
        END DO
        END DO
        END DO
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz write(41,410)
 410  format(/,' Agregation of too THIN Layers')
! #vz write(41,411) (100.*dz_ref(  isn),isn=1,nsnow)
! #vz write(41,412) (100.*dzwdif(  isn),isn=1,nsnow)
! #vz write(41,413) (100.*dzsnSV(1,isn),isn=1,nsnow)
! #vz write(41,414) (              isn ,isn=1,nsnow)
 411  format(' dz_ref [cm]:',10f8.2   ,/,('             ',10f8.2) )
 412  format(' dz_dif [cm]:',10f8.2   ,/,('             ',10f8.2) )
 413  format(' dzsnSV [cm]:',10f8.2   ,/,('             ',10f8.2) )
 414  format('             ',10(i5,3x),/,('             ',10(i5,3x)))
! #vz write(41,4111)      isnoSV(1    )
! #vz write(41,4112)      i_thin(1    )
! #vz write(41,4113)      LIndsv(1    )
! #vz write(41,4114)      Agrege(1    )
! #vz write(41,4115) 1.e2*dzagr1(1    )
! #vz write(41,4116) 1.e2*dzagr2(1    )
 4111 format(' isnoSV     :',  i8  )
 4112 format(' i_thin     :',  i8  )
 4113 format(' LIndsv     :',  i8  )
 4114 format(' Agrege     :',  f8.2)
 4115 format(' dzagr1     :',  f8.2)
 4116 format(' dzagr2     :',  f8.2)
 
! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp   write(47,471)(G1snSV(1,isn),isn=1,isnoSV(1,1))
 471    format('Before _zAg1: G1 = ',10f8.1,(/,19x,10f8.1))
! #vp   write(47,472)(G2snSV(1,isn),isn=1,isnoSV(1,1))
 
 
! Agregates
! ---------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
!            **********
        call SISVAT_zAg                                                &
     &             (ikl,ikv,isagr1(ikl,ikv),isagr2(ikl,ikv),WEagre(ikl,ikv)&
     &                 ,dzagr1(ikl,ikv),dzagr2(ikl,ikv),T_agr1(ikl,ikv),T_agr2(ikl,ikv)&
     &                 ,roagr1(ikl,ikv),roagr2(ikl,ikv),etagr1(ikl,ikv),etagr2(ikl,ikv)&
     &                 ,G1agr1(ikl,ikv),G1agr2(ikl,ikv),G2agr1(ikl,ikv),G2agr2(ikl,ikv)&
     &                 ,agagr1(ikl,ikv),agagr2(ikl,ikv),Agrege(ikl,ikv)&
     &                 )
!            **********
 
        END DO
        END DO
 
 
! Rearranges the Layers
! ---------------------
 
! New (agregated) Snow layer
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isn     =             i_thin(ikl,ikv)
          isn     = min(isn,isn+LIndsv(ikl,ikv))
          isnoSV(ikl,ikv) =         isnoSV(ikl,ikv) -Agrege(ikl,ikv)
          iiceSV(ikl,ikv) =         iiceSV(ikl,ikv)                    &
     &            -max(0,sign(1,iiceSV(ikl,ikv) -isn +icemix))         &
     &                                      *Agrege(ikl,ikv)           &
     &            *max(0,sign(1,iiceSV(ikl,ikv) -1          ))
          istoSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*istoSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *isagr1(ikl,ikv)
          dzsnSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*dzsnSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *dzagr1(ikl,ikv)
          TsisSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*TsisSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *T_agr1(ikl,ikv)
          ro__SV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*ro__SV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *roagr1(ikl,ikv)
          eta_SV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*eta_SV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *etagr1(ikl,ikv)
          G1snSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*G1snSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *G1agr1(ikl,ikv)
          G2snSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*G2snSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *G2agr1(ikl,ikv)
          agsnSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*agsnSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *agagr1(ikl,ikv)
        END DO
        END DO
 
! Above
! ^^^^^
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isn1(ikl,ikv)=max(i_thin(ikl,ikv),i_thin(ikl,ikv)+LIndsv(ikl,ikv))
        END DO
        END DO
        DO i=  1,nsnow-1
        DO ikl=1,kcolp
        DO ikv=1,mwp
            staggr        =  min(1,max(0,i +1 -isn1(ikl,ikv)   ))
            istoSV(ikl,ikv,i) = (1.-staggr     )*istoSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*istoSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *istoSV(ikl,ikv,i+1))
            dzsnSV(ikl,ikv,i) = (1.-staggr     )*dzsnSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*dzsnSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *dzsnSV(ikl,ikv,i+1))
            TsisSV(ikl,ikv,i) = (1.-staggr     )*TsisSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*TsisSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *TsisSV(ikl,ikv,i+1))
            ro__SV(ikl,ikv,i) = (1.-staggr     )*ro__SV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*ro__SV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *ro__SV(ikl,ikv,i+1))
            eta_SV(ikl,ikv,i) = (1.-staggr     )*eta_SV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*eta_SV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *eta_SV(ikl,ikv,i+1))
            G1snSV(ikl,ikv,i) = (1.-staggr     )*G1snSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*G1snSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *G1snSV(ikl,ikv,i+1))
            G2snSV(ikl,ikv,i) = (1.-staggr     )*G2snSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*G2snSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *G2snSV(ikl,ikv,i+1))
            agsnSV(ikl,ikv,i) = (1.-staggr     )*agsnSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*agsnSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *agsnSV(ikl,ikv,i+1))
        END DO
        END DO
        END DO
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isn             = min(isnoSV(ikl,ikv) +1,nsnow)
          istoSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*istoSV(ikl,ikv,isn)
          dzsnSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*dzsnSV(ikl,ikv,isn)
          TsisSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*TsisSV(ikl,ikv,isn)
          ro__SV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*ro__SV(ikl,ikv,isn)
          eta_SV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*eta_SV(ikl,ikv,isn)
          G1snSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*G1snSV(ikl,ikv,isn)
          G2snSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*G2snSV(ikl,ikv,isn)
          agsnSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*agsnSV(ikl,ikv,isn)
        END DO
        END DO
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #s1   IF (ii__AP(ikl).EQ.iwr_SV.AND.jj__AP(ikl).EQ.jwr_SV      .AND. &
! #s1&      ikv        .EQ.nwr_SV                          )     THEN
! #s1     write(6,5991) i_thin(ikl,ikv)
 5991     format(/,'First Agregation / Layer',i3,                      &
     &           /,'  i',11x,'T',9x,'rho',10x,'dz',11x,'H')
! #s1     write(6,5995) (isn,TsisSV(ikl,ikv,isn),ro__SV(ikl,ikv,isn)   &
! #s1&                      ,dzsnSV(ikl,ikv,isn),istoSV(ikl,ikv,isn),  &
! #s1&                   isn=isnoSV(ikl,ikv),1,-1)
 5995     format(i3,3f12.3,i12)
! #s1   END IF
 
 
! Constrains Splitting          of too thick Layers
! =================================================
 
 
! Search the thickest non-zero Layer
! ----------------------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          dzthin(ikl,ikv) =   0.                        ! Arbitrary unrealistic
        END DO
        END DO                                          !       Layer Thickness
      DO   isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isno_n    =             isnoSV(ikl,ikv)-isn+1 ! Snow Normal.Profile
          iice_n    =             iiceSV(ikl,ikv)-isn   ! Ice  Normal.Profile
          iiceOK    = min(1,max(0,iice_n         +1))   ! Ice         Switch
          dz_dif    =(      dzsnSV(ikl,ikv,isn)        &! Actual      Profile
     &        - dz_max *((1-iiceOK)*isno_n*isno_n      &! Theoretical Profile
     &                 +    iiceOK *2.   **iice_n)  )  &!
     &                 /max(dzsnSV(ikl,ikv,isn),eps6)   !
          OKthin      = max(zer0,                      &!
     &                      sign(un_1,                 &!
     &                           dz_dif-dzthin(ikl,ikv)))  &! 1.=>New thickest Lay.
     &                * max(0,                         &! 1 =>.le. isnoSV
     &                  min(1,                         &!
     &                      isnoSV(ikl,ikv)-isn +1 ))   !
          i_thin(ikl,ikv) =    (1. - OKthin)  * i_thin(ikl,ikv)&!  Update thickest Lay.
     &                         + OKthin   * isn         !                Index
          dzthin(ikl,ikv) =    (1. - OKthin)  * dzthin(ikl,ikv)&!
     &                         + OKthin   * dz_dif      !
        END DO
        END DO
      END DO
 
      DO   ikl=1,kcolp
      DO ikv=1,mwp
          ThickL      = max(zer0,                      &! 1. => a too   thick
     &                      sign(un_1,dzthin(ikl,ikv)  &!         Layer exists
     &                               -eps6       ))    &!
     &        * max(0,1-max(0   ,     isnoSV(ikl,ikv)  &! No spliting allowed
     &                               -nsnow+3    ))     ! if isno > nsnow - 3
          Agrege(ikl,ikv) =               ThickL       &! 1. => effective split
     &        * max(0,1-max(0   ,     NLaysv(ikl,ikv)  &!
     &                               +isnoSV(ikl,ikv)  &!
     &                               -nsnow+1    ))     !
          NLay_s(ikl,ikv) =               ThickL       &! Agregation
     &        * max(0,1-max(0   ,     NLaysv(ikl,ikv)  &! to allow  Splitting
     &                               +isnoSV(ikl,ikv)  &!   at next Time Step
     &                               -nsnow      ))    &!
     &                               -Agrege(ikl,ikv)   !
          NLay_s(ikl,ikv) = max(0   ,     NLay_s(ikl,ikv))  ! Agregation effective
      END DO
      END DO
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz     write(41,4152) i_thin(1,1),dzthin(1,1),ThickL
 4152     format(/,' Thickest',i3,':',f9.3,'   Split =',f4.0)
 
 
! Rearranges the Layers
! ---------------------
 
      DO isn=nsnow,2,-1
      DO ikl=1,kcolp
      DO ikv=1,mwp
        IF (Agrege(ikl,ikv).gt.0..AND.i_thin(ikl,ikv).lt.isnoSV(ikl,ikv))   THEN
          staggr          =  min(1,max(0,isn-i_thin(ikl,ikv)    -1))   &
     &                    *  min(1,max(0,    isnoSV(ikl,ikv)-isn+2))
          istoSV(ikl,ikv,isn) =        staggr  * istoSV(ikl,ikv ,isn-1)&
     &                    + (1. -  staggr) * istoSV(ikl,ikv ,isn  )
          dzsnSV(ikl,ikv,isn) =        staggr  * dzsnSV(ikl,ikv ,isn-1)&
     &                    + (1. -  staggr) * dzsnSV(ikl,ikv ,isn  )
          TsisSV(ikl,ikv,isn) =        staggr  * TsisSV(ikl,ikv ,isn-1)&
     &                    + (1. -  staggr) * TsisSV(ikl,ikv ,isn  )
          ro__SV(ikl,ikv,isn) =        staggr  * ro__SV(ikl,ikv ,isn-1)&
     &                    + (1. -  staggr) * ro__SV(ikl,ikv ,isn  )
          eta_SV(ikl,ikv,isn) =        staggr  * eta_SV(ikl,ikv ,isn-1)&
     &                    + (1. -  staggr) * eta_SV(ikl,ikv ,isn  )
          G1snSV(ikl,ikv,isn) =        staggr  * G1snSV(ikl,ikv ,isn-1)&
     &                    + (1. -  staggr) * G1snSV(ikl,ikv ,isn  )
          G2snSV(ikl,ikv,isn) =        staggr  * G2snSV(ikl,ikv ,isn-1)&
     &                    + (1. -  staggr) * G2snSV(ikl,ikv ,isn  )
          agsnSV(ikl,ikv,isn) =        staggr  * agsnSV(ikl,ikv ,isn-1)&
     &                    + (1. -  staggr) * agsnSV(ikl,ikv ,isn  )
        END IF
      END DO
      END DO
      END DO
 
      DO  ikl=1,kcolp
      DO ikv=1,mwp
          isn             =     i_thin(ikl,ikv)
          dzsnSV(ikl,ikv,isn) = 0.5*Agrege(ikl,ikv) *dzsnSV(ikl,ikv,isn)   &
     &                    + (1.-Agrege(ikl,ikv))*dzsnSV(ikl,ikv,isn)
 
          isn             = min(i_thin(ikl,ikv) +1,nsnow)
          istoSV(ikl,ikv,isn) =     Agrege(ikl,ikv) *istoSV(ikl,ikv,isn-1) &
     &                    + (1.-Agrege(ikl,ikv))*istoSV(ikl,ikv,isn)
          dzsnSV(ikl,ikv,isn) =     Agrege(ikl,ikv) *dzsnSV(ikl,ikv,isn-1) &
     &                    + (1.-Agrege(ikl,ikv))*dzsnSV(ikl,ikv,isn)
          TsisSV(ikl,ikv,isn) =     Agrege(ikl,ikv) *TsisSV(ikl,ikv,isn-1) &
     &                    + (1.-Agrege(ikl,ikv))*TsisSV(ikl,ikv,isn)
          ro__SV(ikl,ikv,isn) =     Agrege(ikl,ikv) *ro__SV(ikl,ikv,isn-1) &
     &                    + (1.-Agrege(ikl,ikv))*ro__SV(ikl,ikv,isn)
          eta_SV(ikl,ikv,isn) =     Agrege(ikl,ikv) *eta_SV(ikl,ikv,isn-1) &
     &                    + (1.-Agrege(ikl,ikv))*eta_SV(ikl,ikv,isn)
          G1snSV(ikl,ikv,isn) =     Agrege(ikl,ikv) *G1snSV(ikl,ikv,isn-1) &
     &                    + (1.-Agrege(ikl,ikv))*G1snSV(ikl,ikv,isn)
          G2snSV(ikl,ikv,isn) =     Agrege(ikl,ikv) *G2snSV(ikl,ikv,isn-1) &
     &                    + (1.-Agrege(ikl,ikv))*G2snSV(ikl,ikv,isn)
          agsnSV(ikl,ikv,isn) =     Agrege(ikl,ikv) *agsnSV(ikl,ikv,isn-1) &
     &                    + (1.-Agrege(ikl,ikv))*agsnSV(ikl,ikv,isn)
          isnoSV(ikl,ikv)     =     Agrege(ikl,ikv) +isnoSV(ikl,ikv)
          iiceSV(ikl,ikv)     =                  iiceSV(ikl,ikv)       &
     &                    +     Agrege(ikl,ikv) *max(0,sign(1,iiceSV(ikl,ikv)  &
     &                                                   -isn +icemix))&
     &                                      *max(0,sign(1,iiceSV(ikl,ikv)  &
     &                                                   -1          ))
      END DO
      END DO
 
 
! Constrains Agregation in case of too much  Layers
! =================================================
 
! Search the thinest   non-zero Layer
! -----------------------------------
 
! OUTPUT/Verification: Snow Thinest Layer
! #sd     write( 6,*)   ' '
! #sd     write( 6,*)   'Agregation 2'
! #sd     write( 6,6000) NLaysv(1)
 6000     format(i3,6x,                                                &
     &          'dzsnSV      dz_min      dz_dif      ',                &
     &          'OKthin      dzthin   i_thin')
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          dzthin(ikl,ikv) =   0.                        ! Arbitrary unrealistic
        END DO
        END DO                                          !       Layer Thickness
      DO   isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isno_n    =             isnoSV(ikl,ikv)-isn+1 ! Snow Normal.Profile
          iice_n    =             iiceSV(ikl,ikv)-isn   ! Ice  Normal.Profile
          iiceOK    = min(1,max(0,iice_n         +1))   ! Ice         Switch
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz     dz_ref(isn) =                                &!
! #vz&          dz_min *((1-iiceOK)*isno_n*isno_n      &! Theoretical Profile
! #vz&                 +    iiceOK *    2**iice_n)     &!
! #vz&               /max(1,isnoSV(ikl,ikv))            !
 
          dz_dif      =     dz_min                     &! Actual      Profile
     &                    - dzsnSV(ikl,ikv    ,isn)    &!
     &        /max(eps6,((1-iiceOK)*isno_n*isno_n      &! Theoretical Profile
     &                 +    iiceOK *2.   **iice_n))     !
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz     dzwdif(isn) =     dz_dif                      !
 
          OKthin      = max(zer0,                       &!
     &                      sign(un_1,                  &!
     &                           dz_dif  - dzthin(ikl,ikv)))&! 1.=> New thinest Lay.
     &                * max(0,                          &! 1 => .le. isnoSV
     &                  min(1,                          &!
     &                      isnoSV(ikl,ikv)-isn +1 ))    !
          i_thin(ikl,ikv) =    (1. - OKthin) * i_thin(ikl,ikv)  &! Update   thinest Lay.
     &                         + OKthin  * isn           !                Index
          dzthin(ikl,ikv) =    (1. - OKthin) * dzthin(ikl,ikv)  &!
     &                         + OKthin  * dz_dif        !
 
! OUTPUT/Verification: Snow Thinest Layer
! #sd   IF(isn.LE.isnoSV(1,1).AND.ikl.EQ.1.AND.ikv.EQ.1)                     &
! #sd&    write( 6,6001) isn,dzsnSV(ikl,ikv,isn),dz_min*isno_n*isno_n,dz_dif &
! #sd&               ,OKthin,dzthin(ikl,ikv),    i_thin(ikl,ikv)
 6001     format(i3,5f12.6,i9)
 
        END DO
        END DO
      END DO
 
! OUTPUT/Verification: Snow Thinest Layer
! #sd     write( 6,*)   ' '
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz     write(41,4153) i_thin(1,1),dzsnSV(1,i_thin(1,1))
 4153     format(/,' Thinest ',i3,':',f9.3)
! #vz     write(41,4151) i_thin(1,1),dzsnSV(1,i_thin(1,1))             &
! #vz&                  ,isnoSV(1,1),dzsnSV(1,isnoSV(1,1))
 
! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp   write(47,473)(G1snSV(1,isn),isn=1,isnoSV(1,1))
 473    format('Before _zCr2: G1 = ',10f8.1,(/,19x,10f8.1))
! #vp   write(47,472)(G2snSV(1,isn),isn=1,isnoSV(1,1))
 
 
! Index of the contiguous Layer to agregate
! -----------------------------------------
 
!          **********
      call SISVAT_zCr
!          **********
 
 
! Assign the 2 Layers to agregate
! -------------------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isn         =    i_thin(ikl,ikv)
          isagr1(ikl,ikv) =    istoSV(ikl,ikv,isn)
          isagr2(ikl,ikv) =    istoSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          dzagr1(ikl,ikv) =    dzsnSV(ikl,ikv,isn)
          dzagr2(ikl,ikv) =    dzsnSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          T_agr1(ikl,ikv) =    TsisSV(ikl,ikv,isn)
          T_agr2(ikl,ikv) =    TsisSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          roagr1(ikl,ikv) =    ro__SV(ikl,ikv,isn)
          roagr2(ikl,ikv) =    ro__SV(ikl,ikv,isn+LIndsv(ikl,ikv))
          etagr1(ikl,ikv) =    eta_SV(ikl,ikv,isn)
          etagr2(ikl,ikv) =    eta_SV(ikl,ikv,isn+LIndsv(ikl,ikv))
          G1agr1(ikl,ikv) =    G1snSV(ikl,ikv,isn)
          G1agr2(ikl,ikv) =    G1snSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          G2agr1(ikl,ikv) =    G2snSV(ikl,ikv,isn)
          G2agr2(ikl,ikv) =    G2snSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          agagr1(ikl,ikv) =    agsnSV(ikl,ikv,isn)
          agagr2(ikl,ikv) =    agsnSV(ikl,ikv,isn+LIndsv(ikl,ikv))
          LstLay      = min(1,max(  0,    isnoSV(ikl,ikv)-1   ))
          Agrege(ikl,ikv) = min(1,                                     &
     &                  max(0,                                         &
     &                      NLaysv(ikl,ikv)   +isnoSV(ikl,ikv)-nsnow   &
     &                     +NLay_s(ikl,ikv)                    )       &
     &                                    *LstLay           )
          isnoSV(ikl,ikv) =                    isnoSV(ikl,ikv)         &
     &     -(1-LstLay)*max(zer0,                                       &
     &                     sign(un_1,      eps_21                      &
     &                                    -dzsnSV(ikl,ikv,1)   ))
          isnoSV(ikl,ikv) =max(   0,           isnoSV(ikl,ikv)      )
 
          WEagre(ikl,ikv) = 0.
        END DO
        END DO
 
        DO isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
          WEagre(ikl,ikv) = WEagre(ikl,ikv) + ro__SV(ikl,ikv,isn)*dzsnSV(ikl,ikv,isn)  &
     &                                *min(1,max(0,i_thin(ikl,ikv)+1-isn))
        END DO
        END DO
        END DO
 
! OUTPUT/Verification: Snow Layers Agregation
! #vz write(41,4120)
 4120 format(' Agregation of too MUCH Layers')
! #vz write(41,411) (100.*dz_ref(  isn),isn=1,nsnow)
! #vz write(41,412) (100.*dzwdif(  isn),isn=1,nsnow)
! #vz write(41,413) (100.*dzsnSV(1,isn),isn=1,nsnow)
! #vz write(41,414) (              isn ,isn=1,nsnow)
! #vz write(41,4111)      isnoSV(1    )
! #vz write(41,4112)      i_thin(1    )
! #vz write(41,4113)      LIndsv(1    )
! #vz write(41,4114)      Agrege(1    )
 
! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp   write(47,474)(G1snSV(1,isn),isn=1,isnoSV(1,1))
 474    format('Before _zAg2: G1 = ',10f8.1,(/,19x,10f8.1))
! #vp   write(47,472)(G2snSV(1,isn),isn=1,isnoSV(1,1))
 
 
! Agregates
! ---------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
!            **********
        call SISVAT_zAg                                                &
     &             (ikl,ikv,isagr1(ikl,ikv),isagr2(ikl,ikv),WEagre(ikl,ikv)&
     &                 ,dzagr1(ikl,ikv),dzagr2(ikl,ikv),T_agr1(ikl,ikv),T_agr2(ikl,ikv)&
     &                 ,roagr1(ikl,ikv),roagr2(ikl,ikv),etagr1(ikl,ikv),etagr2(ikl,ikv)&
     &                 ,G1agr1(ikl,ikv),G1agr2(ikl,ikv),G2agr1(ikl,ikv),G2agr2(ikl,ikv)&
     &                 ,agagr1(ikl,ikv),agagr2(ikl,ikv),Agrege(ikl,ikv)&
     &                 )
!            **********
 
        END DO
        END DO
 
 
! Rearranges the Layers
! ---------------------
 
! New (agregated) Snow layer
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isn     =             i_thin(ikl,ikv)
          isn     = min(isn,isn+LIndsv(ikl,ikv))
          isnoSV(ikl,ikv) =         isnoSV(ikl,ikv) -Agrege(ikl,ikv)
          iiceSV(ikl,ikv) =         iiceSV(ikl,ikv)                    &
     &            -max(0,sign(1,iiceSV(ikl,ikv) -isn +icemix))         &
     &                                      *Agrege(ikl,ikv)           &
     &            *max(0,sign(1,iiceSV(ikl,ikv) -1          ))
          istoSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*istoSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *isagr1(ikl,ikv)
          dzsnSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*dzsnSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *dzagr1(ikl,ikv)
          TsisSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*TsisSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *T_agr1(ikl,ikv)
          ro__SV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*ro__SV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *roagr1(ikl,ikv)
          eta_SV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*eta_SV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *etagr1(ikl,ikv)
          G1snSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*G1snSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *G1agr1(ikl,ikv)
          G2snSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*G2snSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *G2agr1(ikl,ikv)
          agsnSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*agsnSV(ikl,ikv,isn)   &
     &                      +   Agrege(ikl,ikv) *agagr1(ikl,ikv)
        END DO
        END DO
 
! Above
! ^^^^^
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isn1(ikl,ikv)=max(i_thin(ikl,ikv),i_thin(ikl,ikv)+LIndsv(ikl,ikv))
        END DO
        END DO
        DO i=  1,nsnow-1
        DO ikl=1,kcolp
        DO ikv=1,mwp
            staggr        =  min(1,max(0,i +1 -isn1(ikl,ikv)   ))
            istoSV(ikl,ikv,i) = (1.-staggr     )*istoSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*istoSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *istoSV(ikl,ikv,i+1))
            dzsnSV(ikl,ikv,i) = (1.-staggr     )*dzsnSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*dzsnSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *dzsnSV(ikl,ikv,i+1))
            TsisSV(ikl,ikv,i) = (1.-staggr     )*TsisSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*TsisSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *TsisSV(ikl,ikv,i+1))
            ro__SV(ikl,ikv,i) = (1.-staggr     )*ro__SV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*ro__SV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *ro__SV(ikl,ikv,i+1))
            eta_SV(ikl,ikv,i) = (1.-staggr     )*eta_SV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*eta_SV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *eta_SV(ikl,ikv,i+1))
            G1snSV(ikl,ikv,i) = (1.-staggr     )*G1snSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*G1snSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *G1snSV(ikl,ikv,i+1))
            G2snSV(ikl,ikv,i) = (1.-staggr     )*G2snSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*G2snSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *G2snSV(ikl,ikv,i+1))
            agsnSV(ikl,ikv,i) = (1.-staggr     )*agsnSV(ikl,ikv,i  )   &
     &            + staggr*((1.-Agrege(ikl,ikv))*agsnSV(ikl,ikv,i  )   &
     &                      +   Agrege(ikl,ikv) *agsnSV(ikl,ikv,i+1))
        END DO
        END DO
        END DO
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          isn             = min(isnoSV(ikl,ikv) +1,nsnow)
          istoSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*istoSV(ikl,ikv,isn)
          dzsnSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*dzsnSV(ikl,ikv,isn)
          TsisSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*TsisSV(ikl,ikv,isn)
          ro__SV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*ro__SV(ikl,ikv,isn)
          eta_SV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*eta_SV(ikl,ikv,isn)
          G1snSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*G1snSV(ikl,ikv,isn)
          G2snSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*G2snSV(ikl,ikv,isn)
          agsnSV(ikl,ikv,isn) = (1.-Agrege(ikl,ikv))*agsnSV(ikl,ikv,isn)
        END DO
        END DO
 
! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp   write(47,475)(G1snSV(1,isn),isn=1,isnoSV(1,1))
 475    format('At End _zSn : G1 = ',10f8.1,(/,19x,10f8.1))
! #vp   write(47,472)(G2snSV(1,isn),isn=1,isnoSV(1,1))
 
 
! Search new Ice/Snow Interface
! =============================
 
! #SX   DO ikl=1,kcolp
! #SX   DO ikv=1,mwp
! #SX     iiceSV(ikl,ikv) =  0
! #SX   END DO
! #SX   END DO
 
! #SX   DO isn=1,nsnow
! #SX   DO ikl=1,kcolp
! #SX   DO ikv=1,mwp
! #SX     OK_ICE      = max(zer0,sign(un_1,ro__SV(ikl,ikv,isn)-850.))  &
! #SX&                * max(zer0,sign(un_1,dzsnSV(ikl,ikv,isn)-eps6))
! #SX     iiceSV(ikl,ikv) = (1.-OK_ICE)       *iiceSV(ikl,ikv)         &
! #SX&                +     OK_ICE        *isn
! #SX   END DO
! #SX   END DO
! #SX   END DO
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
      IF (FlagDALLOC)                                             THEN !
 
      deallocate          ( NLay_s )                                   ! Split Snow Layer         Switch
      deallocate          ( isagr1 )                                   ! 1st     Layer History
      deallocate          ( isagr2 )                                   ! 2nd     Layer History
      deallocate          ( isn1   )                                   ! 1st layer to stagger
      deallocate          ( WEagre )                                   ! Snow Water Equivalent Thickness
      deallocate          ( dzthin )                                   ! Thickness of the thinest layer
      deallocate          ( Agrege )                                   ! 1. when Agregation constrained
      deallocate          ( dzagr1 )                                   ! 1st     Layer Thickness
      deallocate          ( dzagr2 )                                   ! 2nd     Layer Thickness
      deallocate          ( T_agr1 )                                   ! 1st     Layer Temperature
      deallocate          ( T_agr2 )                                   ! 2nd     Layer Temperature
      deallocate          ( roagr1 )                                   ! 1st     Layer Density
      deallocate          ( roagr2 )                                   ! 2nd     Layer Density
      deallocate          ( etagr1 )                                   ! 1st     Layer Water Content
      deallocate          ( etagr2 )                                   ! 2nd     Layer Water Content
      deallocate          ( G1agr1 )                                   ! 1st     Layer Dendricity/Spher.
      deallocate          ( G1agr2 )                                   ! 2nd     Layer Dendricity/Spher.
      deallocate          ( G2agr1 )                                   ! 1st     Layer Sphericity/Size
      deallocate          ( G2agr2 )                                   ! 2nd     Layer Sphericity/Size
      deallocate          ( agagr1 )                                   ! 1st     Layer Age
      deallocate          ( agagr2 )                                   ! 2nd     Layer Age
 
! #vz deallocate          ( dz_ref )                                   ! Snow Reference Discretization
! #vz deallocate          ( dzwdif )                                   !
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SISVAT_zSn
 
 
 
      subroutine SISVAT_zCr
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_zCr                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_zCr determines criteria for Layers Agregation      |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT /  isnoSV   = total Nb of Ice/Snow Layers                      |
!     OUTPUT:  iiceSV   = total Nb of Ice      Layers                      |
!     ^^^^^^   ispiSV   = 0,..,nsnow: Uppermost Superimposed Ice Layer     |
!              istoSV   = 0,..,5    : Snow     History (see istdSV data)   |
!                                                                          |
!     INPUT /  ro__SV   : Soil/Snow Volumic Mass                   [kg/m3] |
!     OUTPUT:           & Snow     Temperatures (layers  1,2,..,nsnow) [K] |
!     ^^^^^^   G1snSV   : Dendricity (<0) or Sphericity (>0) of Snow Layer |
!              G2snSV   : Sphericity (>0) or Size            of Snow Layer |
!              agsnSV   : Snow       Age                             [day] |
!                                                                          |
!     OUTPUT:  LIndsv   : Relative Index of a contiguous Layer to agregate |
!     ^^^^^^                                                               |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
 
 
 
      IMPLICIT NONE
 
 
 
! Internal Variables
! ==================
 
      integer           ::  ikl,ikv   ,isn   ,is0   ,is1
      integer           ::  isno_1                        ! Switch:  ! Snow Layer over Ice
      real(kind=real8)  ::  Dtyp_0,Dtyp_1                 ! Snow Grains Difference Measure
      real(kind=real8)  ::  DenSph                        ! 1. when contiguous spheric
!                                                         !     and dendritic  Grains
      real(kind=real8)  ::  DendOK                        ! 1. when dendritic  Grains
      real(kind=real8)  ::  dTypMx = 200.0                ! Grain Type Differ.
      real(kind=real8)  ::  dTypSp =   0.5                ! Sphericity Weight
      real(kind=real8)  ::  dTypRo =   0.5                ! Density    Weight
      real(kind=real8)  ::  dTypDi =  10.0                ! Grain Diam.Weight
      real(kind=real8)  ::  dTypHi = 100.0                ! History    Weight
 
 
! Agregation Criteria
! ===================
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
          i_thin(ikl,ikv) = min(i_thin(ikl,ikv),isnoSV(ikl,ikv))
          isn         = max(1          ,i_thin(ikl,ikv))
 
 
! Comparison with the downward Layer
! ----------------------------------
 
          is0    = max(1,        i_thin(ikl,ikv)-1 )    ! Downward Layer Index
          DenSph = max(zer0,                           &! isn/is1
     &                 sign(un_1,                      &! Dendricity/Sphericity
     &                      eps6-G1snSV(ikl,ikv,isn)   &!            Switch
     &                          *G1snSV(ikl,ikv,is0)))  !
          DendOK = max(zer0,                           &! Dendricity Switch
     &                 sign(un_1,                      &!
     &                      eps6-G1snSV(ikl,ikv,isn)))  !
 
          Dtyp_0 =                                     &!
     &         DenSph *      dTypMx                    &!
     &    +(1.-DenSph)                                 &!
     &    *    DendOK *((abs(G1snSV(ikl,ikv,isn)       &! Dendricity
     &                      -G1snSV(ikl,ikv,is0))      &!     Contribution
     &                  +abs(G2snSV(ikl,ikv,isn)       &! Sphericity
     &                      -G2snSV(ikl,ikv,is0))) *dTypSp &!     Contribution
     &                  +abs(ro__SV(ikl,ikv,isn)       &! Density
     &                      -ro__SV(ikl,ikv,is0))  *dTypRo)&!     Contribution
     &    +(1.-DenSph)                                 &!
     &    *(1.-DendOK)*((abs(G1snSV(ikl,ikv,isn)       &! Sphericity
     &                      -G1snSV(ikl,ikv,is0))      &!     Contribution
     &                  +abs(G2snSV(ikl,ikv,isn)       &! Size
     &                      -G2snSV(ikl,ikv,is0))) *dTypDi &!     Contribution
     &                  +abs(ro__SV(ikl,ikv,isn)       &! Density
     &                      -ro__SV(ikl,ikv,is0))  *dTypRo) !     Contribution
          Dtyp_0 =                                     &!
     &                   min(dTypMx,                   &!
     &                       Dtyp_0                    &!
     &                  +abs(istoSV(ikl,ikv,isn)       &! History
     &                      -istoSV(ikl,ikv,is0))  *dTypHi)&!     Contribution
     &        +             (1 -abs(isn-is0))  * 1.e+6 &!"Same Layer"Score
     &        +  max(0,1-abs(iiceSV(ikl,ikv)           &!"Ice /Snow
     &                                 -is0))  * 1.e+6  ! Interface" Score
 
 
! Comparison with the   upward Layer
! ----------------------------------
 
          is1    = min(          i_thin(ikl,ikv)+1,    &! Upward   Layer Index
     &                 max(1,    isnoSV(ikl,ikv)  ))    !
          DenSph = max(zer0,                           &! isn/is1
     &                 sign(un_1,                      &! Dendricity/Sphericity
     &                      eps6-G1snSV(ikl,ikv,isn)   &!            Switch
     &                          *G1snSV(ikl,ikv,is1)))  !
          DendOK = max(zer0,                           &! Dendricity Switch
     &                 sign(un_1,                      &!
     &                      eps6-G1snSV(ikl,ikv,isn)))  !
 
          Dtyp_1 =                                     &!
     &         DenSph *      dTypMx                    &!
     &    +(1.-DenSph)                                 &!
     &    *    DendOK *((abs(G1snSV(ikl,ikv,isn)       &! Dendricity
     &                      -G1snSV(ikl,ikv,is1))      &!     Contribution
     &                  +abs(G2snSV(ikl,ikv,isn)       &! Sphericity
     &                      -G2snSV(ikl,ikv,is1))) *dTypSp &!     Contribution
     &                  +abs(ro__SV(ikl,ikv,isn)       &! Density
     &                      -ro__SV(ikl,ikv,is1))  *dTypRo)&!     Contribution
     &    +(1.-DenSph)                                 &!
     &    *(1.-DendOK)*((abs(G1snSV(ikl,ikv,isn)       &! Sphericity
     &                      -G1snSV(ikl,ikv,is1))      &!     Contribution
     &                  +abs(G2snSV(ikl,ikv,isn)       &! Size
     &                      -G2snSV(ikl,ikv,is1))) *dTypDi &!     Contribution
     &                  +abs(ro__SV(ikl,ikv,isn)       &! Density
     &                      -ro__SV(ikl,ikv,is1))  *dTypRo) !     Contribution
          Dtyp_1 =                                     &!
     &                   min(dTypMx,                   &!
     &                       Dtyp_1                    &!
     &                  +abs(istoSV(ikl,ikv,isn)       &! History
     &                      -istoSV(ikl,ikv,is1))  *dTypHi)&!     Contribution
     &        +             (1 -abs(isn-is1))  * 1.e+6 &!"Same Layer"Score
     &        +  max(0,1-abs(iiceSV(ikl,ikv)           &!"Ice /Snow
     &                                 -isn))  * 1.e+6  ! Interface" Score
 
 
! Index of the Layer to agregate
! ==============================
 
          LIndsv(ikl,ikv) = sign(un_1,Dtyp_0           &!
     &                           -Dtyp_1)               !
          isno_1      = (1 -min (abs(isnoSV(ikl,ikv)   &! Switch = 1
     &                              -iiceSV(ikl,ikv)-1),1))&!   if isno = iice +1
     &                * (1 -min (abs(isnoSV(ikl,ikv)   &! Switch = 1
     &                              -i_thin(ikl,ikv)  ),1)) !   if isno = i_ithin
          LIndsv(ikl,ikv) = (1 -isno_1) *LIndsv(ikl,ikv)   &! Contiguous Layer is
     &                     -isno_1                      ! downward for top L.
          i_thin(ikl,ikv) =  max(1,   i_thin(ikl,ikv)   )   !
      END DO
      END DO
 
 
      return
      end subroutine SISVAT_zCr
 
 
 
      subroutine SISVAT_zAg                                           &!
     &                 (ikl,ikv,isagra,isagrb,WEagra                  &!
     &                     ,dzagra,dzagrb,T_agra,T_agrb               &!
     &                     ,roagra,roagrb,etagra,etagrb               &!
     &                     ,G1agra,G1agrb,G2agra,G2agrb               &!
     &                     ,agagra,agagrb,Agreg1                      &!
     &                     )
 
!--------------------------------------------------------------------------+
!   MAR SURFACE                                       Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_zAg aggregates two contiguous snow layers          |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   isagrb   : 2nd     Layer History                            |
!     ^^^^^                                                                |
!                                                                          |
!     INPUT:   dzagrb   : 2nd     Layer Thickness                          |
!     ^^^^^    T_agrb   : 2nd     Layer Temperature                        |
!              roagrb   : 2nd     Layer Density                            |
!              etagrb   : 2nd     Layer Water Content                      |
!              G1agrb   : 2nd     Layer Dendricity/Spher.                  |
!              G2agrb   : 2nd     Layer Sphericity/Size                    |
!              agagrb   : 2nd     Age                                      |
!              Agreg1   : 1. when Agregation constrained                   |
!                                                                          |
!     INPUT /  isagra   : 1st     Layer History                            |
!     OUTPUT:                                                              |
!     ^^^^^^                                                               |
!                                                                          |
!     INPUT /  dzagra   : 1st     Layer Thickness                          |
!     OUTPUT:  T_agra   : 1st     Layer Temperature                        |
!     ^^^^^^   roagra   : 1st     Layer Density                            |
!              etagra   : 1st     Layer Water Content                      |
!              G1agra   : 1st     Layer Dendricity/Spher.                  |
!              G2agra   : 1st     Layer Sphericity/Size                    |
!              agagra   : 1st     Age                                      |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # stdout               | #s1: OUTPUT of Snow Layers Agregation         |
!                          |      unit  6, SubRoutine  SISVAT_zSn, _zAg    |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
 
 
      IMPLICIT NONE
 
 
 
! Transferred Variables
! =====================
 
 
! INPUT
! -----
 
      integer           ::  ikl,ikv                       !    Column Index
 
      integer           ::  isagrb                        ! 2nd Layer History
      real(kind=real8)  ::  dzagrb                        ! 2nd Layer Thickness
      real(kind=real8)  ::  T_agrb                        ! 2nd Layer Temperature
      real(kind=real8)  ::  roagrb                        ! 2nd Layer Density
      real(kind=real8)  ::  etagrb                        ! 2nd Layer Water Content
      real(kind=real8)  ::  G1agrb                        ! 2nd Layer Dendricity/Spher.
      real(kind=real8)  ::  G2agrb                        ! 2nd Layer Sphericity/Size
      real(kind=real8)  ::  agagrb                        ! 2nd Layer Age
 
 
! INPUT/OUTPUT
! ------------
 
      integer           ::  isagra                        ! 1st Layer History
      real(kind=real8)  ::  WEagra                        ! 1st Layer Height    [mm w.e.]
      real(kind=real8)  ::  Agreg1                        ! 1. ===>   Agregates
      real(kind=real8)  ::  dzagra                        ! 1st Layer Thickness
      real(kind=real8)  ::  T_agra                        ! 1st Layer Temperature
      real(kind=real8)  ::  roagra                        ! 1st Layer Density
      real(kind=real8)  ::  etagra                        ! 1st Layer Water Content
      real(kind=real8)  ::  G1agra                        ! 1st Layer Dendricity/Spher.
      real(kind=real8)  ::  G2agra                        ! 1st Layer Sphericity/Size
      real(kind=real8)  ::  agagra                        ! 1st Layer Age
 
 
 
! Internal Variables
! ==================
 
      integer           ::  nh                            ! Averaged    Snow History
      integer           ::  nh__OK                        ! 1=>Conserve Snow History
      real(kind=real8)  ::  rh                            !
      real(kind=real8)  ::  dz                            ! Thickness
      real(kind=real8)  ::  dzro_1                        ! Thickness X Density, Lay.1
      real(kind=real8)  ::  dzro_2                        ! Thickness X Density, Lay.2
      real(kind=real8)  ::  dzro                          ! Thickness X Density, Aver.
      real(kind=real8)  ::  ro                            ! Averaged    Density
      real(kind=real8)  ::  wn                            ! Averaged    Water Content
      real(kind=real8)  ::  tn                            ! Averaged    Temperature
      real(kind=real8)  ::  ag                            ! Averaged    Snow Age
      real(kind=real8)  ::  SameOK                        ! 1. => Same Type of Grains
      real(kind=real8)  ::  G1same                        ! Averaged G1,  same Grains
      real(kind=real8)  ::  G2same                        ! Averaged G2,  same Grains
      real(kind=real8)  ::  typ__1                        ! 1. => Lay1 Type: Dendritic
      real(kind=real8)  ::  zroNEW                        ! dz X ro, if fresh Snow
      real(kind=real8)  ::  G1_NEW                        ! G1,      if fresh Snow
      real(kind=real8)  ::  G2_NEW                        ! G2,      if fresh Snow
      real(kind=real8)  ::  zroOLD                        ! dz X ro, if old   Snow
      real(kind=real8)  ::  G1_OLD                        ! G1,      if old   Snow
      real(kind=real8)  ::  G2_OLD                        ! G2,      if old   Snow
      real(kind=real8)  ::  SizNEW                        ! Size,    if fresh Snow
      real(kind=real8)  ::  SphNEW                        ! Spheric.,if fresh Snow
      real(kind=real8)  ::  SizOLD                        ! Size,    if old   Snow
      real(kind=real8)  ::  SphOLD                        ! Spheric.,if old   Snow
      real(kind=real8)  ::  Siz_av                        ! Averaged    Grain Size
      real(kind=real8)  ::  Sph_av                        ! Averaged    Grain Spher.
      real(kind=real8)  ::  Den_av                        ! Averaged    Grain Dendr.
      real(kind=real8)  ::  DendOK                        ! 1. => Average is  Dendr.
      real(kind=real8)  ::  G1diff                        ! Averaged G1, diff. Grains
      real(kind=real8)  ::  G2diff                        ! Averaged G2, diff. Grains
      real(kind=real8)  ::  G1                            ! Averaged G1
      real(kind=real8)  ::  G2                            ! Averaged G2
 
 
 
! Mean   Properties
! =================
 
!  1 Densite, Contenu en Eau, Temperature /
!    Density, Water Content,  Temperature
!    ------------------------------------
 
          dz      =  dzagra      + dzagrb
          dzro_1  =  roagra      * dzagra
          dzro_2  =  roagrb      * dzagrb
          dzro    =  dzro_1      + dzro_2
          ro      =  dzro                                              &
     &     /max(eps6,dz)
          wn      = (dzro_1*etagra      + dzro_2*etagrb     )          &
     &     /max(eps6,dzro)
          tn      = (dzro_1*T_agra      + dzro_2*T_agrb     )          &
     &     /max(eps6,dzro)
          ag      = (dzro_1*agagra      + dzro_2*agagrb     )          &
     &     /max(eps6,dzro)
 
          rh      =  max(zer0,sign(un_1,zWEcSV(ikl,ikv)                &
     &                                         -0.5*WEagra     ))
 
          nh__OK  =  rh
          nh      =                 max(isagra     ,isagrb     )       &
! #HB&            *  nh__OK                                            &
! #HB&          + (1-nh__OK)*       min(isagra     ,isagrb     )       &
     &          +  0.
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #s1   IF (ii__AP(ikl).EQ.iwr_SV.AND.jj__AP(ikl).EQ.jwr_SV      .AND. &
! #s1&      ikv        .EQ.nwr_SV                          )       THEN
! #s1     write(6,5995) zWEcSV(ikl,ikv),WEagra                         &
! #s1&                 ,isagra     ,isagrb                             &
! #s1&                 ,nh__OK     ,nh
 5995     format(' WE2,WEa =',2f9.1,'  nha,b =',2i2,'  nh__OK,nh =',2i2)
! #s1   END IF
 
 
!  2 Nouveaux Types de Grains /  new Grain Types
!    -------------------------------------------
 
!  2.1. Meme  Type  de Neige  / same Grain Type
!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          SameOK  =  max(zer0,                                         &
     &                   sign(un_1, G1agra      *G1agrb       - eps_21))
          G1same  = (dzro_1*G1agra      + dzro_2*G1agrb     )          &
     &     /max(eps6,dzro)
          G2same  = (dzro_1*G2agra      + dzro_2*G2agrb     )          &
     &     /max(eps6,dzro)
 
!  2.2. Types differents / differents Types
!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          typ__1  =  max(zer0,sign(un_1,eps6-G1agra     )) ! =1.=> Dendritic
          zroNEW  =     typ__1  *dzro_1                   &! ro of Dendr.Lay.
     &            + (1.-typ__1) *dzro_2                    !
          G1_NEW  =     typ__1  *G1agra                   &! G1 of Dendr.Lay.
     &            + (1.-typ__1) *G1agrb                    !
          G2_NEW  =     typ__1  *G2agra                   &! G2 of Dendr.Lay.
     &            + (1.-typ__1) *G2agrb                    !
          zroOLD  = (1.-typ__1) *dzro_1                   &! ro of Spher.Lay.
     &            +     typ__1  *dzro_2                    !
          G1_OLD  = (1.-typ__1) *G1agra                   &! G1 of Spher.Lay.
     &            +     typ__1  *G1agrb                    !
          G2_OLD  = (1.-typ__1) *G2agra                   &! G2 of Spher.Lay.
     &            +     typ__1  *G2agrb                    !
          SizNEW  =    -G1_NEW  *DDcdSV/G1_dSV            &! Size  Dendr.Lay.
     &             +(1.+G1_NEW         /G1_dSV)           &!
     &                *(G2_NEW  *DScdSV/G1_dSV            &!
     &             +(1.-G2_NEW         /G1_dSV)*DFcdSV)    !
          SphNEW  =     G2_NEW         /G1_dSV             ! Spher.Dendr.Lay.
          SizOLD  =     G2_OLD                             ! Size  Spher.Lay.
          SphOLD  =     G1_OLD         /G1_dSV             ! Spher.Spher.Lay.
          Siz_av  = (zroNEW*SizNEW+zroOLD*SizOLD)     &! Averaged Size
     &     /max(eps6,dzro)                             !
          Sph_av  = (zroNEW*SphNEW+zroOLD*SphOLD)     &! Averaged Sphericity
     &     /max(eps6,dzro)                             !
          Den_av  = (Siz_av -(    Sph_av *DScdSV      &!
     &                       +(1.-Sph_av)*DFcdSV))    &!
     &            / (DDcdSV -(    Sph_av *DScdSV      &!
     &                       +(1.-Sph_av)*DFcdSV))     !
          DendOK  = max(zer0,                         &!
     &                  sign(un_1,     Sph_av *DScdSV &! Small   Grains Contr.
     &                            +(1.-Sph_av)*DFcdSV &! Faceted Grains Contr.
     &                            -    Siz_av        ))!
!         REMARQUE: le  type moyen (dendritique ou non) depend
!         ^^^^^^^^  de la  comparaison avec le diametre optique
!                   d'une neige recente de   dendricite nulle
!         REMARK:   the mean type  (dendritic   or not) depends
!         ^^^^^^    on the comparaison with the optical diameter
!                   of a recent snow    having zero dendricity
 
          G1diff  =(   -DendOK *Den_av                                 &
     &             +(1.-DendOK)*Sph_av) *G1_dSV
          G2diff  =     DendOK *Sph_av  *G1_dSV                        &
     &             +(1.-DendOK)*Siz_av
          G1      =     SameOK *G1same                                 &
     &             +(1.-SameOK)*G1diff
          G2      =     SameOK *G2same                                 &
     &             +(1.-SameOK)*G2diff
 
 
! Assignation to new Properties
! =============================
 
          isagra        = Agreg1      *nh +(1.-Agreg1     ) *isagra
          dzagra        = Agreg1      *dz +(1.-Agreg1     ) *dzagra
          T_agra        = Agreg1      *tn +(1.-Agreg1     ) *T_agra
          roagra        = Agreg1      *ro +(1.-Agreg1     ) *roagra
          etagra        = Agreg1      *wn +(1.-Agreg1     ) *etagra
          G1agra        = Agreg1      *G1 +(1.-Agreg1     ) *G1agra
          G2agra        = Agreg1      *G2 +(1.-Agreg1     ) *G2agra
          agagra        = Agreg1      *ag +(1.-Agreg1     ) *agagra
 
 
      return
      end subroutine SISVAT_zAg
 
 
 
      subroutine SnOptP(                                               &
! #AG&                  jjtime                                         &
     &                 )
 
!--------------------------------------------------------------------------+
!   MAR/SISVAT   SnOptP                               Wed 26-Jun-2013  MAR |
!     SubRoutine SnOptP computes the Snow Pack optical Properties          |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     Grid Boxes       |
!                       X       Number of Mosaic Cell per Grid Box         |
!                                                                          |
!     INPUT:   isnoSV   = total Nb of Ice/Snow Layers                      |
!     ^^^^^    ispiSV   = 0,...,nsno: Uppermost Superimposed Ice Layer     |
!                                                                          |
!              ivgtSV   = 0,...,12:   Vegetation Type                      |
!                         0:          Water, Solid or Liquid               |
!                                                                          |
!     INPUT:   G1snSV   : Dendricity (<0) or Sphericity (>0) of Snow Layer |
!     ^^^^^    G2snSV   : Sphericity (>0) or Size            of Snow Layer |
!              agsnSV   : Snow       Age                             [day] |
!              ro__SV   : Snow/Soil  Volumic Mass                  [kg/m3] |
!              eta_SV   : Water      Content                       [m3/m3] |
!              rusnSV   : Surficial  Water   Thickness   [kg/m2] .OR. [mm] |
!              SWS_SV   : Surficial  Water   Status                        |
!              dzsnSV   : Snow       Layer   Thickness                 [m] |
!                                                                          |
!              albssv   : Soil       Albedo                            [-] |
!              zzsnsv   : Snow       Pack    Thickness                 [m] |
!                                                                          |
!     OUTPUT:  albisv   : Snow/Ice/Water/Soil Integrated Albedo        [-] |
!     ^^^^^^   sEX_sv   : Verticaly Integrated Extinction Coefficient      |
!                                                                          |
!     Internal Variables:                                                  |
!     ^^^^^^^^^^^^^^^^^^                                                   |
!              SnOpSV   : Snow Grain optical Size                      [m] |
!              EX1_sv   : Integrated Snow Extinction (0.3--0.8micr.m)      |
!              EX2_sv   : Integrated Snow Extinction (0.8--1.5micr.m)      |
!              EX3_sv   : Integrated Snow Extinction (1.5--2.8micr.m)      |
!                                                                          |
!     METHODE:    Calcul de la taille optique des grains ? partir de       |
!     ^^^^^^^    -leur type decrit par les deux variables descriptives     |
!                      continues sur la plage -99/+99 passees en appel.    |
!                -la taille optique (1/10mm) des etoiles,                  |
!                                            des grains fins et            |
!                                            des jeunes faces planes       |
!                                                                          |
!     METHOD:     Computation of the optical diameter of the grains        |
!     ^^^^^^      described with the CROCUS formalism G1snSV / G2snSV      |
!                                                                          |
!     REFERENCE: Brun et al.      1989, J. Glaciol 35 pp. 333--342         |
!     ^^^^^^^^^  Brun et al.      1992, J. Glaciol 38 pp.  13-- 22         |
!                Eric Martin Sept.1996                                     |
!                                                                          |
!     CAUTION: Vegetation is not taken into account in albedo computations |
!     ^^^^^^^  Suggestion: 1) Reduce the displacement height  and/or LAI   |
!              (when snow)    for radiative transfert through vegetation   |
!                          2) Adapt leaf optical parameters                |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Possibility                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                          |
!     #CZ: Albedo Correction (Zenith Angle) (Warren,      1982)            |
!     #cz: Albedo Correction (Zenith Angle) (Segal etAl., 1991) (obsolete) |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Col de Porte                         |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^                                      |
!     #cp: Col de Porte Integrated Snow/Ice Albedo                         |
!     #AG: Snow Aging Col de Porte     (Brun et al.1991)                   |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SnOptP____.va        | #va: OUTPUT/Verification: Albedo Parameteriz. |
!                          |      unit 46, SubRoutine  SnOptP     **ONLY** |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_ctr
      use Mod_SISVAT_CdP
      use Mod_SISVAT_dat
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
! #va use Mod_SISVAT_SnOptP
      use Mod_SISVATLSnO
 
 
 
      IMPLICIT NONE
 
 
 
! Internal Variables
! ==================
 
      real(kind=real8)                               ::  coalbm                ! weighted Coalbedo, mean
! #AG real(kind=real8)                               ::  agesno
 
! #AG integer                                        ::  jjtime                !
      integer                                        ::  isn   ,ikl,ikv        !
! #va integer                                        ::  isn1                  !
 
! For the computation of the solar irradiance extinction in snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(kind=real8)                               ::  sbeta1 = 0.0192
      real(kind=real8)                               ::  sbeta2 = 0.4000
      real(kind=real8)                               ::  sbeta3 = 0.1098
      real(kind=real8)                               ::  sbeta4 = 1.0000
      real(kind=real8)                               ::  sbeta5 = 2.00e1
 
! Snow Age Maximum (Taiga, e.g. Col de Porte)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #AG real(kind=real8)                               ::  AgeMax = 60.0         ! [day]
 
      real(kind=real8)                               ::  AlbMin =  0.94        !  Albedo Minimum / visible (0.3--0.8 micrometers)
      real(kind=real8)                               ::  AlbMax =  0.99        !  Albedo Maximum
      real(kind=real8)                               ::  HSnoSV =  0.01        !  Snow   Thickness over witch interpolate Albedo to Ice  Albedo
      real(kind=real8)                               ::  HIceSV =  0.10        !  Snow   Thickness over witch interpolate Albedo to Soil Albedo
 
      real(kind=real8)                               ::  doptmx =  2.3e-3      ! [m]     Maximum   optical Diameter    (pi * R**2)
 
      real(kind=real8)                               ::  SignG1,Sph_OK
      real(kind=real8)                               ::  dalbed                !
! #cz real(kind=real8)                               ::  dalbeS                !
! #CZ real(kind=real8)                               ::  dalbeW                !
 
! #CZ real(kind=real8)                               ::  bsegal =  4.00
! #CZ real(kind=real8)                               ::  czeMAX =  0.173648178 ! 80.deg (Segal et al., 1991 JAS)
! #CZ real(kind=real8)                               ::  CZ_eff
 
      real(kind=real8)                               ::  RoFrez,SignRo,SnowOK,OpSqrt
      real(kind=real8)                               ::  albSn1,a_SII1
      real(kind=real8)                               ::  albSn2,a_SII2
      real(kind=real8)                               ::  albSn3,a_SII3
      real(kind=real8)                               ::  albSno
! #va real(kind=real8)                               ::  albIce,albIc1,albIc2,albIc3
      real(kind=real8)                               ::  albSII,albWIc
      real(kind=real8)                               ::  doptic,Snow_H,SIce_H,SnownH,SIcenH
      real(kind=real8)                               ::  exarg1,exarg2,exarg3,sign_0,sExt_0
      real(kind=real8)                               ::  albedo_old
      real(kind=real8)                               ::  ro_ave,dz_ave
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
      allocate            ( coalb1(kcolp,mwp) )                        ! weighted Coalbedo, Vis.
      allocate            ( coalb2(kcolp,mwp) )                        ! weighted Coalbedo, nIR 1
      allocate            ( coalb3(kcolp,mwp) )                        ! weighted Coalbedo, nIR 2
      allocate            ( sExt_1(kcolp,mwp) )                        ! Extinction Coeff., Vis.
      allocate            ( sExt_2(kcolp,mwp) )                        ! Extinction Coeff., nIR 1
      allocate            ( sExt_3(kcolp,mwp) )                        ! Extinction Coeff., nIR 2
      allocate            ( SnOpSV(kcolp,mwp,nsnow) )                  ! Snow Grain optical Size
 
      allocate            ( alb1sv(kcolp,mwp) )                        !
      allocate            ( alb2sv(kcolp,mwp) )                        !
      allocate            ( alb3sv(kcolp,mwp) )                        !
 
      END IF
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! Snow Grain optical Size
! =======================
 
      DO   isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
          G2snSV(ikl,ikv,isn) =  max(eps6,G2snSV(ikl,ikv,isn))
!         Avoid non physical Values
 
          SignG1          = sign(un_1,G1snSV(ikl,ikv,isn))
          Sph_OK          =  max(zer0,SignG1)
 
          SnOpSV(ikl,ikv,isn) =   1.e-4 *                              &
!         SI:           (from 1/10 mm to m)
 
 
! Contribution of Non Dendritic Snow
! ----------------------------------
 
     &    (    Sph_OK *(      G2snSV(ikl,ikv,isn)*G1snSV(ikl,ikv,isn)/G1_dSV   &
     &              +max(half*G2snSV(ikl,ikv,isn),DFcdSV)              &
     &                 *(1.00-G1snSV(ikl,ikv,isn)                /G1_dSV)) &
 
 
! Contribution of     Dendritic Snow
! ----------------------------------
 
     &    +(1.-Sph_OK)*(     -G1snSV(ikl,ikv,isn)*DDcdSV         /G1_dSV   &
     &                 +(1.00+G1snSV(ikl,ikv,isn)                /G1_dSV)  &
     &                  *    (G2snSV(ikl,ikv,isn)*DScdSV         /G1_dSV   &
     &                 +(1.00-G2snSV(ikl,ikv,isn)                /G1_dSV)  &
     &                                       *DFcdSV                 )))
          SnOpSV(ikl,ikv,isn) =  max(zer0,SnOpSV(ikl,ikv,isn))
        END DO
        END DO
      END DO
 
 
! Snow/Ice Albedo
! ===============
 
! Snow Age (Influence on Albedo)
! ------------------------------
 
! #AG IF (iabs(mod(jjtime,86400)).lt.dt__SV)                        THEN
! #AG   DO isn=1,nsnow
! #AG   DO ikl=1,kcolp
! #AG   DO ikv=1,mwp
! #AG     agsnSV(ikl,ikv,isn) = agsnSV(ikl,ikv,isn) + 1.               &
! #AG&           + max(zer0,DH_dSV(ivgtSV(ikl,ikv))-DH_dSV(4)) ! High Vegetation
!                                                          ! Impurities
!                           CAUTION: crude parameterization
!                           ^^^^^^^
! #AG   END DO
! #AG   END DO
! #AG   END DO
! #AG END IF
 
 
! Uppermost effective Snow Layer
! ------------------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
           isn   =  max(1,isnoSV(ikl,ikv))
 
          SignRo = sign(un_1, rocdSV - ro__SV(ikl,ikv,isn))
          SnowOK =  max(zer0,SignRo) ! Ice Density Threshold
 
          OpSqrt = sqrt(SnOpSV(ikl,ikv,isn))
 
          albSn1 =  0.96-1.580*OpSqrt
          albSn1 =  max(albSn1,AlbMin)
 
          albSn1 =  max(albSn1,zer0)
          albSn1 =  min(albSn1,un_1)
 
          albSn2 =  0.95-15.40*OpSqrt
          albSn2 =  max(albSn2,zer0)
          albSn2 =  min(albSn2,un_1)
 
          doptic =  min(SnOpSV(ikl,ikv,isn),doptmx)
          albSn3 =  346.3*doptic -32.31*OpSqrt +0.88
          albSn3 =  max(albSn3,zer0)
          albSn3 =  min(albSn3,un_1)
 
          albSno =  So1dSV*albSn1                                      &
     &           +  So2dSV*albSn2                                      &
     &           +  So3dSV*albSn3
 
          SnowOK =  SnowOK*max(zer0,sign(un_1,albSno-aI3dSV))
                 !  Minimum snow albedo is aI3dSV
 
          albSn1 =  SnowOK*albSn1+(1.0-SnowOK)*max(albSno,aI3dSV)
          albSn2 =  SnowOK*albSn2+(1.0-SnowOK)*max(albSno,aI3dSV)
          albSn3 =  SnowOK*albSn3+(1.0-SnowOK)*max(albSno,aI3dSV)
 
 
! Snow/Ice Pack Thickness
! -----------------------
 
          isn    =     max(min(isnoSV(ikl,ikv) ,ispiSV(ikl,ikv)),0)
          Snow_H =  zzsnsv(ikl,ikv,isnoSV(ikl,ikv))-zzsnsv(ikl,ikv,isn)
          SIce_H =  zzsnsv(ikl,ikv,isnoSV(ikl,ikv))
          SnownH =  Snow_H  /  HSnoSV
          SnownH =  min(un_1,  SnownH)
          SIcenH =  SIce_H  / (HIceSV                                  &
     &           +  max(zer0,Z0mdSV(ivgtSV(ikl,ikv))                   &
     &           -           Z0mdSV(4)         ))
          SIcenH =  min(un_1,  SIcenH)
 
!      The value of SnownH is set to 1 in case of ice lenses above
!      1m of dry snow (ro<700kg/m3) for using CROCUS albedo
 
          ro_ave =  0.
          dz_ave =  0.
          SnowOK =  1.
       DO isn    =  isnoSV(ikl,ikv),1,-1
          ro_ave =  ro_ave + ro__SV(ikl,ikv,isn) * dzsnSV(ikl,ikv,isn) * SnowOK
          dz_ave =  dz_ave +                   dzsnSV(ikl,ikv,isn) * SnowOK
          SnowOK =  max(zer0,sign(un_1,1.-dz_ave))
       END DO
 
          ro_ave =  ro_ave / max(dz_ave,eps6)
          SnowOK =  max(zer0,sign(un_1,700.-ro_ave))
 
          SnownH =  SnowOK + SnownH * (1. - SnowOK)
 
 
! Integrated Snow/Ice Albedo: Case of Water on Bare Ice
! -----------------------------------------------------
 
          isn    =  max(min(isnoSV(ikl,ikv) ,ispiSV(ikl,ikv)),0)
 
          albWIc =  aI1dSV-(aI1dSV-aI2dSV)                &!
     &           *  exp(-rusnSV(ikl,ikv)                  &!
     &           *  (1. -SWS_SV(ikl,ikv)                  &! 0 <=> freezing
     &           *  (1  -min(1,iabs(isn-isnoSV(ikl,ikv)))))   &! 1 <=> isn=isnoSV
     &                  /ru_dSV)                           !
 
          SignRo = sign(un_1,rhoIce-1.-ro__SV(ikl,ikv,isn))! RoSN<920kg/m3
          SnowOK =  max(zer0,SignRo)
 
          albWIc = (1. - SnowOK) * albWIc + SnowOK        &!
     &           * (aI2dSV + (aI3dSV -aI2dSV)             &!
     &           * (ro__SV(ikl,ikv,isn)-rhoIce)/(rocdSV-rhoIce))
 
!    rocdSV < ro < rhoIce | aI2dSV< al >aI3dSV (fct of density))
!             ro > rhoIce | aI1dSV< al >aI2dSV (fct of superficial water content)s
 
 
! Integrated Snow/Ice      Albedo
! -------------------------------
 
          a_SII1      =     albWIc      +(albSn1-albWIc)     *SnownH
          a_SII1      = min(a_SII1       ,albSn1)
 
          a_SII2      =     albWIc      +(albSn2-albWIc)     *SnownH
          a_SII2      = min(a_SII2       ,albSn2)
 
          a_SII3      =     albWIc      +(albSn3-albWIc)     *SnownH
          a_SII3      = min(a_SII3       ,albSn3)
 
! #AG     agesno =      min(agsnSV(ikl,ikv,isn)          ,AgeMax)
! #AG     a_SII1      =     a_SII1      -0.175*agesno/AgeMax
!                                        Impurities: Col de Porte Parameter.
 
 
!    Zenith Angle Correction (Segal et al.,         1991, JAS 48, p.1025)
!    ----------------------- (Wiscombe & Warren, dec1980, JAS   , p.2723)
!                            (Warren,               1982,  RG   , p.  81)
!                            --------------------------------------------
 
 
          dalbed = 0.0
! #CZ     CZ_eff = max(czemax                   ,coszSV(ikl,ikv))
! #cz     dalbeS =   ((bsegal+1.00)/(1.00+2.0*bsegal*CZ_eff)           &
! #cz&                -       1.00                          )*0.32     &
! #cz&              /  bsegal
! #cz     dalbeS = max(dalbeS,zer0)
! #cz     dalbed =     dalbeS      *       min(1,isnoSV(ikl,ikv))
 
! #CZ     dalbeW =(0.64 - CZ_eff  )*0.0625  ! Warren 1982, RevGeo, fig.12b
                                            ! 0.0625 = 5% * 1/0.8,   p.81
                                            ! 0.64   = cos(50)
! #CZ     dalbed =     dalbeW      *       min(1,isnoSV(ikl,ikv))
 
! Col de Porte Integrated Snow/Ice Albedo
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (ColPrt.AND.TotSol.gt.0.)                              THEN
            albSII =  (((Dr_1SN*a_SII1+Dr_2SN*a_SII2+Dr_3SN*a_SII3)    &
     &                  +dalbed                                    )   &
     &                  *DirSol                                        &
     &                 +(Df_1SN*a_SII1+Df_2SN*a_SII2+Df_3SN*a_SII3)    &
     &                  *DifSol*(1.   -cld_SV(ikl,ikv))                &
     &                 +(Dfc1SN*a_SII1+Dfc2SN*a_SII2+Dfc3SN*a_SII3)    &
     &                  *DifSol*       cld_SV(ikl,ikv)                  )  &
     &                 / TotSol
 
! Elsewhere    Integrated Snow/Ice Albedo
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ELSE
            albSII =     So1dSV*a_SII1                                 &
     &                 + So2dSV*a_SII2                                 &
     &                 + So3dSV*a_SII3
          END IF
 
 
! Integrated Snow/Ice/Soil Albedo
! -------------------------------
 
            alb1sv(ikl,ikv) =     albssv(ikl,ikv) +(a_SII1-albssv(ikl,ikv))*SIcenH
            alb1sv(ikl,ikv) = min(alb1sv(ikl,ikv)  ,a_SII1)
 
            alb2sv(ikl,ikv) =     albssv(ikl,ikv) +(a_SII2-albssv(ikl,ikv))*SIcenH
            alb2sv(ikl,ikv) = min(alb2sv(ikl,ikv)  ,a_SII2)
 
            alb3sv(ikl,ikv) =     albssv(ikl,ikv) +(a_SII3-albssv(ikl,ikv))*SIcenH
            alb3sv(ikl,ikv) = min(alb3sv(ikl,ikv)  ,a_SII3)
 
            albisv(ikl,ikv) =     albssv(ikl,ikv) +(albSII-albssv(ikl,ikv))*SIcenH
            albisv(ikl,ikv) = min(albisv(ikl,ikv)  ,albSII)
 
 
! Integrated Snow/Ice/Soil Albedo: Clouds Correction! Greuell & all., 1994
! --------------------------------------------------! Glob.&t Planet.Change
                                                    ! (9):91-114
          IF (.NOT.ColPrt)                                          THEN
            alb1sv(ikl,ikv) = alb1sv(ikl,ikv) + 0.05 *(cld_SV(ikl,ikv)-0.5)*SIcenH &
! #CZ&                  + dalbed      *    (1.-cld_SV(ikl,ikv))        &
     &                  + 0.
            alb2sv(ikl,ikv) = alb2sv(ikl,ikv) + 0.05 *(cld_SV(ikl,ikv)-0.5)*SIcenH &
! #CZ&                  + dalbed      *    (1.-cld_SV(ikl,ikv))        &
     &                  + 0.
            alb3sv(ikl,ikv) = alb3sv(ikl,ikv) + 0.05 *(cld_SV(ikl,ikv)-0.5)*SIcenH &
! #CZ&                  + dalbed      *    (1.-cld_SV(ikl,ikv))        &
     &                  + 0.
            albisv(ikl,ikv) = albisv(ikl,ikv) + 0.05 *(cld_SV(ikl,ikv)-0.5)*SIcenH &
! #CZ&                  + dalbed      *    (1.-cld_SV(ikl,ikv))        &
     &                  + 0.
          END IF
 
 
! Integrated Snow/Ice/Soil Albedo: Minimum snow albedo = 40%
! ----------------------------------------------------------
 
            albedo_old  = albisv(ikl,ikv)
 
            albisv(ikl,ikv) = max(albisv(ikl,ikv),0.400    * SIcenH    &
     &                  + albssv(ikl,ikv) *(1.0        - SIcenH))
            alb1sv(ikl,ikv) = alb1sv(ikl,ikv) - 1.0/3.0                &! 33 %
     &                  * (albedo_old-albisv(ikl,ikv)) / So1dSV
            alb2sv(ikl,ikv) = alb2sv(ikl,ikv) - 1.0/3.0                &! 33 %
     &                  * (albedo_old-albisv(ikl,ikv)) / So2dSV
            alb3sv(ikl,ikv) = alb3sv(ikl,ikv) - 1.0/3.0                &! 33 %
     &                  * (albedo_old-albisv(ikl,ikv)) / So3dSV
 
 
! Integrated Snow/Ice/Soil Albedo: Maximum albedo = 99%
! -----------------------------------------------------
 
            albedo_old  = albisv(ikl,ikv)
            albisv(ikl,ikv) = min(albisv(ikl,ikv),0.99)
            alb1sv(ikl,ikv) = alb1sv(ikl,ikv) - 1.0/3.0                &! 33 %
     &                  * (albedo_old-albisv(ikl,ikv)) / So1dSV
            alb2sv(ikl,ikv) = alb2sv(ikl,ikv) - 1.0/3.0                &! 33 %
     &                  * (albedo_old-albisv(ikl,ikv)) / So2dSV
            alb3sv(ikl,ikv) = alb3sv(ikl,ikv) - 1.0/3.0                &! 33 %
     &                  * (albedo_old-albisv(ikl,ikv)) / So3dSV
 
            alb1sv(ikl,ikv) = min(max(zer0,alb1sv(ikl,ikv)),albmax)
            alb2sv(ikl,ikv) = min(max(zer0,alb2sv(ikl,ikv)),albmax)
            alb3sv(ikl,ikv) = min(max(zer0,alb3sv(ikl,ikv)),albmax)
 
        END DO
        END DO
 
 
! Extinction Coefficient: Exponential Factor
! ==========================================
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          sExt_1(ikl,ikv)         = 1.
          sExt_2(ikl,ikv)         = 1.
          sExt_3(ikl,ikv)         = 1.
          sEX_sv(ikl,ikv,nsnow+1) = 1.
 
          coalb1(ikl,ikv) = (1.          -alb1sv(ikl,ikv))*So1dSV
          coalb2(ikl,ikv) = (1.          -alb2sv(ikl,ikv))*So2dSV
          coalb3(ikl,ikv) = (1.          -alb3sv(ikl,ikv))*So3dSV
          coalbm      =  coalb1(ikl,ikv) +coalb2(ikl,ikv) +coalb3(ikl,ikv)
          coalb1(ikl,ikv) =  coalb1(ikl,ikv)              /coalbm
          coalb2(ikl,ikv) =  coalb2(ikl,ikv)              /coalbm
          coalb3(ikl,ikv) =  coalb3(ikl,ikv)              /coalbm
        END DO
        END DO
 
      DO   isn=  nsnow,1,-1
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
          SignRo = sign(un_1, rocdSV - ro__SV(ikl,ikv,isn))
          SnowOK =  max(zer0,SignRo) ! Ice Density Threshold
 
          RoFrez =  1.e-3      * ro__SV(ikl,ikv,isn) * (1.0-eta_SV(ikl,ikv,isn))
 
          OpSqrt = sqrt(max(eps6,SnOpSV(ikl,ikv,isn)))
          exarg1 =      SnowOK  *1.e2 *max(sbeta1*RoFrez/OpSqrt,sbeta2)&
     &            +(1.0-SnowOK)           *sbeta5
          exarg2 =      SnowOK  *1.e2 *max(sbeta3*RoFrez/OpSqrt,sbeta4)&
     &            +(1.0-SnowOK)           *sbeta5
          exarg3 =      SnowOK  *1.e2     *sbeta5                      &
     &            +(1.0-SnowOK)           *sbeta5
 
! Col de Porte Snow Extinction Coefficient
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (ColPrt.AND.TotSol.gt.0.)                              THEN
            exarg1 = exarg1*(Dr_1SN*DirSol                             &
     &                      +Df_1SN*DifSol*(1.-cld_SV(ikl,ikv))        &
     &                      +Dfc1SN*DifSol*    cld_SV(ikl,ikv) )       &
     &                     /(Dr_1SN*TotSol)
            exarg2 = exarg2*(Dr_2SN*DirSol                             &
     &                      +Df_2SN*DifSol*(1.-cld_SV(ikl,ikv))        &
     &                      +Dfc2SN*DifSol*    cld_SV(ikl,ikv) )       &
     &                     /(Dr_2SN*TotSol)
            exarg3 = exarg3*(Dr_3SN*DirSol                             &
     &                      +Df_3SN*DifSol*(1.-cld_SV(ikl,ikv))        &
     &                      +Dfc3SN*DifSol*    cld_SV(ikl,ikv) )       &
     &                     /(Dr_3SN*TotSol)
          END IF
 
 
! Integrated Extinction of Solar Irradiance (Normalized Value)
! ============================================================
 
          sExt_1(ikl,ikv) = sExt_1(ikl,ikv)                            &
     &                          * exp(min(0.0,-exarg1 *dzsnSV(ikl,ikv,isn)))
          sign_0      =              sign(un_1,epsn   -sExt_1(ikl,ikv))
          sExt_0      =               max(zer0,sign_0)*sExt_1(ikl,ikv)
          sExt_1(ikl,ikv) = sExt_1(ikl,ikv)                   -sExt_0
 
          sExt_2(ikl,ikv) = sExt_2(ikl,ikv)                            &
     &                          * exp(min(0.0,-exarg2 *dzsnSV(ikl,ikv,isn)))
          sign_0      =              sign(un_1,epsn   -sExt_2(ikl,ikv))
          sExt_0      =               max(zer0,sign_0)*sExt_2(ikl,ikv)
          sExt_2(ikl,ikv) = sExt_2(ikl,ikv)                   -sExt_0
 
          sExt_3(ikl,ikv) = sExt_3(ikl,ikv)                            &
     &                          * exp(min(0.0,-exarg3 *dzsnSV(ikl,ikv,isn)))
          sign_0      =              sign(un_1,epsn   -sExt_3(ikl,ikv))
          sExt_0      =               max(zer0,sign_0)*sExt_3(ikl,ikv)
          sExt_3(ikl,ikv) = sExt_3(ikl,ikv)                   -sExt_0
 
          sEX_sv(ikl,ikv,isn) = coalb1(ikl,ikv) * sExt_1(ikl,ikv)      &
     &                    + coalb2(ikl,ikv) * sExt_2(ikl,ikv)          &
     &                    + coalb3(ikl,ikv) * sExt_3(ikl,ikv)
        END DO
        END DO
      END DO
 
      DO   isn=0,-nsoil,-1
        DO ikl=1,kcolp
        DO ikv=1,mwp
          sEX_sv(ikl,ikv,isn) = 0.0
        END DO
        END DO
      END DO
 
 
! Albedo: IO
! ==========
 
! #va IF (.NOT.aw_opn)                                              THEN
! #va          aw_opn = .true.
! #va          open(unit=46,status='unknown',file='SnOptP____.va')
! #va          rewind(   46)
! #va END IF
 
! #va     ikl,ikv=1
! #va     write(46,460)daHost
! #va 460 format('---------------------------------+----+',            &
! #va&          '-------+-------+-------+-------+-------+-------+',    &
! #va&                                  '-------+-------+-------+',    &
! #va&         /,'Snow/Ice Pack ',a18,' |    |',                       &
! #va&          ' z [m] |0.3/0.8|0.8/1.5|1.5/2.8| Full  |Opt[mm]|',    &
! #va&                                  ' G1    | G2    | ro    |',    &
! #va&         /,'---------------------------------+----+',            &
! #va&          '-------+-------+-------+-------+-------+-------+',    &
! #va&                                  '-------+-------+-------+')
!         ______________________________________________________________
! #va     write(46,461)            SIce_H,                             &
! #va&                             alb1sv(ikl,ikv),alb2sv(ikl,ikv),alb3sv(ikl,ikv),&
! #va&                             albisv(ikl,ikv)
! #va 461 format('Integrated Snow/Ice/Soil  Albedo |',                 &
! #va&            3x,' |',  f6.3,' |' ,4(f6.3,' |'), 6x ,' |',         &
! #va&                                            3( 6x ,' |'))
!         ______________________________________________________________
! #va     write(46,462)ispiSV(ikl,ikv),a_SII1,a_SII2,a_SII3,albSII
! #va 462 format('Integrated Snow/Ice       Albedo |',                 &
! #va&            i3,' |',   6x ,' |' ,4(f6.3,' |'), 6x ,' |',         &
! #va&                                            3( 6x ,' |'))
!         ______________________________________________________________
! #va     write(46,463)            rusnSV(ikl,ikv),         albWIc,    &
! #va&                             SWS_SV(ikl,ikv)
! #va 463 format('Integrated Water/Bare Ice Albedo |',                 &
! #va&            3x,' |',  f6.3,'w|' ,3( 6x, ' |'),                   &
! #va&                                   f6.3,' |' ,f6.3,' |',         &
! #va&                                            3( 6x ,' |'))
!         ______________________________________________________________
! #va     write(46,465)isn1       ,zzsnsv(ikl,ikv,isn1),               &
! #va&                             albIc1,albIc2,albIc3,albIce,        &
! #va&                        1.e3*SnOpSV(ikl,ikv,max(1,isnoSV(ikl,ikv)-1)),   &
! #va&                             G1snSV(ikl,ikv,max(1,isnoSV(ikl,ikv)-1)),   &
! #va&                             G2snSV(ikl,ikv,max(1,isnoSV(ikl,ikv)-1)),   &
! #va&                             ro__SV(ikl,ikv,max(1,isnoSV(ikl,ikv)-1))&
! #va&                      *(1. - eta_SV(ikl,ikv,max(1,isnoSV(ikl,ikv)-1)))
! #va 465 format('Surficial       Ice Lense        |',                 &
! #va&            i3,' |', (f6.3,'i|'),4(f6.3,' |'),f6.3,' |',         &
! #va&                                            3(f6.1,' |'))
!         ______________________________________________________________
! #va     write(46,466)isnoSV(ikl,ikv),zzsnsv(ikl,ikv,isnoSV(ikl,ikv)),&
! #va&                             albSn1,albSn2,albSn3,albSno,        &
! #va&                        1.e3*SnOpSV(ikl,ikv,isnoSV(ikl,ikv)),    &
! #va&                             G1snSV(ikl,ikv,isnoSV(ikl,ikv)),    &
! #va&                             G2snSV(ikl,ikv,isnoSV(ikl,ikv)),    &
! #va&                             ro__SV(ikl,ikv,isnoSV(ikl,ikv))     &
! #va&                      *(1. - eta_SV(ikl,ikv,isnoSV(ikl,ikv)))
! #va 466 format('Uppermost  Effective Snow Layer  |',                 &
! #va&            i3,' |', (f6.3,'*|'),4(f6.3,' |'),f6.3,' |',         &
! #va&                                            3(f6.1,' |'))
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
      IF (FlagDALLOC)                                             THEN !
 
      deallocate          ( coalb1 )                                   ! weighted Coalbedo, Vis.
      deallocate          ( coalb2 )                                   ! weighted Coalbedo, nIR 1
      deallocate          ( coalb3 )                                   ! weighted Coalbedo, nIR 2
      deallocate          ( sExt_1 )                                   ! Extinction Coeff., Vis.
      deallocate          ( sExt_2 )                                   ! Extinction Coeff., nIR 1
      deallocate          ( sExt_3 )                                   ! Extinction Coeff., nIR 2
      deallocate          ( SnOpSV )                                   ! Snow Grain optical Size
 
      deallocate          ( alb1sv )                                   !
      deallocate          ( alb2sv )                                   !
      deallocate          ( alb3sv )                                   !
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SnOptP
 
 
 
      subroutine VgOptP
 
!--------------------------------------------------------------------------+
!   MAR/SISVAT   VgOptP                               Wed 26-Jun-2013  MAR |
!     SubRoutine VgOptP computes the Canopy    optical Properties          |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     Grid Boxes       |
!                       X       Number of Mosaic Cell per Grid Box         |
!                                                                          |
!     INPUT:   ivgtSV   = 0,...,12:   Vegetation Type                      |
!     ^^^^^               0:          Water, Solid or Liquid               |
!                                                                          |
!     INPUT:   coszSV   : Cosine of the Sun Zenithal Distance          [-] |
!     ^^^^^    sol_SV   : Surface Downward  Solar   Radiation       [W/m2] |
!              snCaSV   : Canopy     Snow      Thickness         [mm w.e.] |
!                                                                          |
!              LAI_sv   : Leaf Area  Index      (snow included)        [-] |
!              glf_sv   : Green Leaf Fraction of NOT fallen Leaves     [-] |
!              albisv   : Snow/Ice/Water/Soil Integrated Albedo        [-] |
!                                                                          |
!     OUTPUT:  alb_SV   : Surface-Canopy Albedo                        [-] |
!     ^^^^^^   SoCasv   : Absorbed Solar Radiation by Canopy (Normaliz)[-] |
!              SoSosv   : Absorbed Solar Radiation by Surfac (Normaliz)[-] |
!              LAIesv   : Effective Leaf Area  Index for Transpiration [-] |
!                                                                          |
!     Internal Variables: Normalized Values:                               |
!     ^^^^^^^^^^^^^^^^^^                                                   |
!              u0_Vis   : Upward   Visible Radiation at Top Canopy     [-] |
!              absg_V   : Absorbed Visible Radiation by the Ground     [-] |
!              absv_V   : Absorbed Visible Radiation by the Canopy     [-] |
!              u0_nIR   : Upward   Near IR Radiation at Top Canopy     [-] |
!              absgnI   : Absorbed Near IR Radiation by the Ground     [-] |
!              absv_V   : Absorbed Near IR Radiation by the Canopy     [-] |
!                                                                          |
!     REFERENCE:   De Ridder, 1997, unpublished thesis, chapter 2 (DR97,2) |
!     ^^^^^^^^^                                                            |
!                                                                          |
!     ASSUMPTIONS: Leaf Inclination Index chi_l (eqn2.49 DR97) set to zero |
!     ^^^^^^^^^^^                         for all vegetation types         |
!                  Radiation Fluxes are normalized                         |
!                        with respect to incoming solar radiation (=I0+D0) |
!                                                                          |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_TRV
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
 
 
 
! Internal Variables
! ==================
 
      use Mod_SISVATLVgO
 
 
      IMPLICIT NONE
 
 
      integer             ::  ikl,ikv   ,kri
 
      real(kind=real8)    ::  exdRad,k_drad
      real(kind=real8)    ::  e_prad,e1pRad
      real(kind=real8)    ::  zv_fac,zv1fac,deadLF
      real(kind=real8)    ::  T_Rad0,A_Rad0
      real(kind=real8)    ::                r0_Rad,t0_Rad,nu_Rad
      real(kind=real8)    ::  Tr_Rad,Re_Rad,r__Rad,t__Rad,t1_Rad
      real(kind=real8)    ::  arggam, gamma
      real(kind=real8)    ::  gammaL
      real(kind=real8)    ::  denSig,Sig__c
      real(kind=real8)    ::  DDifH1,DDifC1
      real(kind=real8)    ::  DDifH2,DDifC2
      real(kind=real8)    ::  denS_s,denS_a,den_c1,DDif_L
      real(kind=real8)    ::  u0_Vis,absg_V,absv_V
      real(kind=real8)    ::  u0_nIR,absgnI,absvnI
      real(kind=real8)    ::  argexg,argexk
      real(kind=real8)    ::  residu,d_DDif,dDDifs,dDDifa
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
      allocate              ( k___sv(kcolp,mwp) )
      allocate              ( A0__sv(kcolp,mwp) )
      allocate              ( gamasv(kcolp,mwp) )
      allocate              ( Sigcsv(kcolp,mwp) )
      allocate              ( C1__sv(kcolp,mwp) )
      allocate              ( C2__sv(kcolp,mwp) )
      allocate              ( criLAI(kcolp,mwp) )
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! General Parameters, Solar Radiation Absorption
! ==============================================
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
 
            k_dRad = 0.5 /max(coszSV(ikl,ikv),eps6)           ! absorbed irradiance fraction
            e_pRad = 2.5   *  coszSV(ikl,ikv)                 ! exponential argument,
                                                              ! V/nIR radiation partitioning,
                                                              ! DR97, 2, eqn (2.53) & (2.54)
            exdRad =    exp(-min(k_dRad*LAI_sv(ikl,ikv),ea_MAX))  ! exponential, Irradi. Absorpt.
            e1pRad = 1.-exp(-    e_pRad                    )  ! exponential, V/nIR Rad. Part.
 
            ivg    =                ivgtSV(ikl,ikv)           ! Vegetation Type
            zv_fac =    min( snCaSV(ikl,ikv)/snCaMx          &! Contribution of Snow to Leaf
     &                      ,  un_1)                          ! Reflectivity and Transmissiv.
            zv1fac = 1.     -       zv_fac                    !
            deadLF = 1.     -       glf_sv(ikl,ikv)           ! Dead Leaf Fraction
 
 
! Visible Part of the Solar Radiation Spectrum (V,   0.4--0.7mi.m)
! ================================================================
 
            A_Rad0 =      0.25 + 0.697 * e1pRad ! Absorbed    Vis. Radiation
            T_Rad0 = 1. - A_Rad0                ! Transmitted Vis  Radiation
 
! Reflectivity, Transmissivity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            Re_Rad = glf_sv(ikl,ikv) *ReVisL(ivg)                      &
     &             + deadLF      *ReVisD(ivg)
            Tr_Rad = glf_sv(ikl,ikv) *TrVisL(ivg)                      &
     &             + deadLF      *TrVisD(ivg)
 
! Adaptation to Snow
! ^^^^^^^^^^^^^^^^^^
            Re_Rad = zv1fac      *Re_Rad      + zv_fac *reVisS
            Tr_Rad = zv1fac      *Tr_Rad      + zv_fac *trVisS
 
! Scattering /DR97, 2, eqn (2.26) and (2.27)           ! Diffuse  Radiation:
! ^^^^^^^^^^                                           ! ^^^^^^^^^^^^^^^^^^
            r__Rad = (2. *Re_Rad +     Tr_Rad) / 3.    ! Upw.  Scatter.Fract.
            t__Rad = (    Re_Rad + 2. *Tr_Rad) / 3.    ! Downw.Scatter.Fract.
 
            t1_Rad =  1. -t__Rad                       !
            arggam =      t1_Rad*t1_Rad-r__Rad*r__Rad  !
            arggam =  max(arggam,zer0)                 !
            gamma  = sqrt(arggam)                      ! eqn (2.39)
            gammaL =  min( gamma*LAI_sv(ikl,ikv),40.0) !
            DDifH1 =  exp( gammaL           )          ! Downw.Diffus.Solut.1
            DDifH2 =  exp(-gammaL           )          ! Downw.Diffus.Solut.2
!           REMARK:  These 2 contributions are zero in case of 0 Reflectivity
!           ^^^^^^
 
! Scattering /DR97, 2, eqn (2.19) and (2.20)           ! Direct   Radiation:
! ^^^^^^^^^^                                           ! ^^^^^^^^^^^^^^^^^^
            r0_Rad = 0.5 *((Re_Rad+Tr_Rad) *k_dRad    &! Upw.  Scatter.Fract.
     &                    +(Re_Rad-Tr_Rad) /    3.)    !
            t0_Rad = 0.5 *((Re_Rad+Tr_Rad) *k_dRad    &! Downw.Scatter.Fract.
     &                    -(Re_Rad-Tr_Rad) /    3.)    !
 
            nu_Rad = t1_Rad-r__Rad*albisv(ikl,ikv)     ! nu coeff., eqn 2.43
            den_c1 =  gamma*(DDifH1+DDifH2)           &! eqn (2.43) Denomin.
     &              +nu_Rad*(DDifH1-DDifH2)            !(Constant for DDifH1)
 
            denSig =  gamma*gamma - k_dRad*k_dRad      ! eqn (2.40) Denomin.
            denS_s = sign(un_1,denSig)                 !
            denS_a =  abs(     denSig)                 !
            denSig =  max(eps6,denS_a) * denS_s        !
            Sig__c = (r__Rad* r0_Rad                  &! sigma_c, eqn (2.40)
     &               +t0_Rad*(k_dRad+t1_Rad)) / denSig !
 
            DDifC1 = ((gamma-nu_Rad)*(T_Rad0-Sig__c*A_Rad0)*DDifH2     &
     &             +((k_dRad-nu_Rad)* Sig__c                           &
     &               +t0_Rad+r__Rad * albisv(ikl,ikv)) *A_Rad0 *exdRad)&
     &           /max(den_c1,eps6)
            DDifC2 =  T_Rad0        - DDifC1-Sig__c*A_Rad0
 
! Visible Diffuse Fluxes
! ^^^^^^^^^^^^^^^^^^^^^^
            DDif_L =  DDifC1*DDifH1 + DDifC2*DDifH2   &! DOWNward,
     &             +  Sig__c*A_Rad0 *exdRad            ! Canopy Basis
            u0_Vis = ((gamma+t1_Rad)*DDifC1           &! UPward
     &               -(gamma-t1_Rad)*DDifC2           &! Canopy Top
     &             -((k_dRad-t1_Rad)*Sig__c           &!
     &               +t0_Rad               )*A_Rad0)  &!
     &          / max(r__Rad,eps6)                     !
            u0_Vis = min(0.99,max(eps6,u0_Vis))        ! ERROR
            absg_V = (1.-albisv(ikl,ikv))*(A_Rad0*exdRad  &! Ground Absorption
     &                                +DDif_L       )  !
            absv_V = (1.-u0_Vis     )- absg_V          ! Veget. Absorption
 
! Parameters for Computing Effective LAI for Transpiration
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            gamasv(ikl,ikv) = gamma
            C1__sv(ikl,ikv) = DDifC1
            C2__sv(ikl,ikv) = DDifC2
            Sigcsv(ikl,ikv) = Sig__c
            k___sv(ikl,ikv) = k_dRad
            A0__sv(ikl,ikv) = A_Rad0
 
 
! Near-IR Part of the Solar Radiation Spectrum (nIR, 0.7--2.8mi.m)
! ================================================================
 
            A_Rad0 =      0.80 + 0.185 * e1pRad ! Absorbed    nIR. Radiation
            T_Rad0 = 1. - A_Rad0                ! Transmitted nIR  Radiation
 
! Reflectivity, Transmissivity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            Re_Rad = glf_sv(ikl,ikv) *RenIRL(ivg)                      &
     &             + deadLF      *RenIRD(ivg)
            Tr_Rad = glf_sv(ikl,ikv) *TrnIRL(ivg)                      &
     &             + deadLF      *TrnIRD(ivg)
 
! Adaptation to Snow
! ^^^^^^^^^^^^^^^^^^
            Re_Rad = zv1fac      *Re_Rad      + zv_fac *renIRS
            Tr_Rad = zv1fac      *Tr_Rad      + zv_fac *trnIRS
 
! Scattering /DR97, 2, eqn (2.26) and (2.27)           ! Diffuse  Radiation:
! ^^^^^^^^^^                                           ! ^^^^^^^^^^^^^^^^^^
            r__Rad = (2. *Re_Rad +     Tr_Rad) / 3.    ! Upw.  Scatter.Fract.
            t__Rad = (    Re_Rad + 2. *Tr_Rad) / 3.    ! Downw.Scatter.Fract.
 
            t1_Rad =  1. -t__Rad                       !
            arggam =      t1_Rad*t1_Rad-r__Rad*r__Rad  !
            arggam =  max(arggam,zer0)                 !
            gamma  = sqrt(arggam)                      ! eqn (2.39)
            DDifH1 =  exp( gamma*LAI_sv(ikl,ikv))      ! Downw.Diffus.Solut.1
            DDifH2 =  exp(-gamma*LAI_sv(ikl,ikv))      ! Downw.Diffus.Solut.2
!           REMARK:  These 2 contributions are zero in case of 0 Reflectivity
!           ^^^^^^
 
! Scattering /DR97, 2, eqn (2.19) and (2.20)           ! Direct   Radiation:
! ^^^^^^^^^^                                           ! ^^^^^^^^^^^^^^^^^^
            r0_Rad = 0.5 *((Re_Rad+Tr_Rad) *k_dRad    &! Upw.  Scatter.Fract.
     &                    +(Re_Rad-Tr_Rad) /    3.)    !
            t0_Rad = 0.5 *((Re_Rad+Tr_Rad) *k_dRad    &! Downw.Scatter.Fract.
     &                    -(Re_Rad-Tr_Rad) /    3.)    !
 
            nu_Rad = t1_Rad-r__Rad*albisv(ikl,ikv)     ! nu coeff., eqn 2.43
            den_c1 =  gamma*(DDifH1+DDifH2)           &! eqn (2.43) Denomin.
     &              +nu_Rad*(DDifH1-DDifH2)            !(Constant for DDifH1)
 
            denSig =  gamma*gamma - k_dRad*k_dRad      ! eqn (2.40) Denomin.
            denS_s = sign(un_1,denSig)                 !
            denS_a =  abs(     denSig)                 !
            denSig =  max(eps6,denS_a) * denS_s        !
            Sig__c = (r__Rad* r0_Rad                  &! sigma_c, eqn (2.40)
     &               +t0_Rad*(k_dRad+t1_Rad)) / denSig !
 
            DDifC1 = ((gamma-nu_Rad)*(T_Rad0-Sig__c*A_Rad0)*DDifH2     &
     &             +((k_dRad-nu_Rad)* Sig__c                           &
     &               +t0_Rad+r__Rad * albisv(ikl,ikv)) *A_Rad0 *exdRad)&
     &           /max(den_c1,eps6)
            DDifC2 =  T_Rad0        - DDifC1-Sig__c*A_Rad0
 
! Near IR Diffuse Fluxes
! ^^^^^^^^^^^^^^^^^^^^^^
            DDif_L =  DDifC1*DDifH1 + DDifC2*DDifH2   &! DOWNward,
     &             +  Sig__c*A_Rad0 *exdRad            ! Canopy Basis
            u0_nIR = ((gamma+t1_Rad)*DDifC1           &! UPward
     &               -(gamma-t1_Rad)*DDifC2           &! Canopy Top
     &             -((k_dRad-t1_Rad)*Sig__c           &!
     &               +t0_Rad               )*A_Rad0)  &!
     &          / max(r__Rad,eps6)                     !
            u0_nIR = min(0.99,max(eps6,u0_nIR))        ! ERROR
            absgnI = (1.-albisv(ikl,ikv))*(A_Rad0*exdRad  &! Ground Absorption
     &                                +DDif_L       )  !
            absvnI = (1.-u0_nIR     )- absgnI          ! Veget. Absorption
 
 
! Surface-Canopy Albedo and Normalized Solar Radiation Absorption
! ===============================================================
 
            alb_SV(ikl,ikv) = (u0_Vis+u0_nIR)*0.5d0
            SoCasv(ikl,ikv) = (absv_V+absvnI)*0.5d0
            SoSosv(ikl,ikv) = (absg_V+absgnI)*0.5d0
 
      END DO
      END DO
 
 
! Effective LAI for Transpiration
! ===============================
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
              criLAI(ikl,ikv) = 2.              ! LAI for which D0_Vis > 20W/m2
                                                ! DR97, 2, eqn (2.57)
        END DO
        END DO
 
      DO   kri=1,10
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
              argexg      =  min(criLAI(ikl,ikv)*gamasv(ikl,ikv),     ea_Max)
              argexk      =  min(criLAI(ikl,ikv)*k___sv(ikl,ikv),     ea_Max)
              residu      =      C1__sv(ikl,ikv)            *exp( argexg)  &
     &                          +C2__sv(ikl,ikv)            *exp(-argexg)  &
     &                          +A0__sv(ikl,ikv)*gamasv(ikl,ikv)*exp(-argexk)  &
     &                          -CriStR /max(sol_SV(ikl,ikv),       eps6)
 
              d_DDif      =      C1__sv(ikl,ikv)*gamasv(ikl,ikv)*exp( argexg)  &
     &                          -C2__sv(ikl,ikv)*gamasv(ikl,ikv)*exp(-argexg)  &
     &                          -A0__sv(ikl,ikv)*k___sv(ikl,ikv)*exp(-argexk)
              dDDifs      = sign(un_1,d_DDif)
              dDDifa      =  abs(     d_DDif)
              d_DDif      =  max(eps6,dDDifa) * dDDifs
 
              criLAI(ikl,ikv) =      criLAI(ikl,ikv)-residu/d_DDif
              criLAI(ikl,ikv) =  max(criLAI(ikl,ikv),zer0       )
              criLAI(ikl,ikv) =  min(criLAI(ikl,ikv),LAI_sv(ikl,ikv))
 
        END DO
        END DO
      END DO
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
              LAIesv(ikl,ikv) = criLAI(ikl,ikv) +(exp(-min(k___sv(ikl,ikv)*criLAI(ikl,ikv),ea_MAX)) &
     &                                   -exp(-min(k___sv(ikl,ikv)*LAI_sv(ikl,ikv),ea_MAX)))&
     &                                  /          k___sv(ikl,ikv)
        END DO
        END DO
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
      IF (FlagDALLOC)                                             THEN !
 
      deallocate            ( k___sv )
      deallocate            ( A0__sv )
      deallocate            ( gamasv )
      deallocate            ( Sigcsv )
      deallocate            ( C1__sv )
      deallocate            ( C2__sv )
      deallocate            ( criLAI )
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine VgOptP
 
 
 
      subroutine ColPrt_SBL
 
!--------------------------------------------------------------------------+
!   MAR          ColPrt_SBL                           Wed 26-Jun-2013  MAR |
!     SubRoutine ColPrt_SBL generates Surface Boundary Layers Properties   |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns                          |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   za__SV   : Surface Boundary Layer (SBL) Height          [m] |
!     ^^^^^    VV__SV   :(SBL Top)   Wind Velocity                   [m/s] |
!              TaT_SV   : SBL Top    Temperature                       [K] |
!              rhT_SV   : SBL Top  Air  Density                    [kg/m3] |
!              uqs_SV   : Specific  Humidity  Turbulent Flux         [m/s] |
!              Tsrfsv   : Surface    Temperature                       [K] |
!                                                                          |
!     INPUT /  LMO_SV   : Monin-Obukhov       Scale                    [m] |
!     OUTPUT:  us__SV   : Friction  Velocity                         [m/s] |
!     ^^^^^^   uts_SV   : Temperature         Turbulent Flux       [K.m/s] |
!                                                                          |
!     OUTPUT:  ram_sv   : Aerodynamic Resistance for Momentum        [s/m] |
!     ^^^^^^   rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!                                                                          |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_dat
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
 
 
 
      IMPLICIT NONE
 
 
 
! Internal Variables
! ==================
 
      integer           ::  ikl,ikv    ,ist    ,ist__s ,ist__w
      real(kind=real8)  ::  d_TaTs ,CD_m
      real(kind=real8)  ::  uustar ,thstar ,qqstar
      real(kind=real8)  ::  thstarv,thstars,thstara
      real(kind=real8)  ::  zeta   ,zeta_S ,zeta_A
      real(kind=real8)  ::  fCdCdP =  3.09                  ! Drag Coefficient Factor, Col de Porte
      real(kind=real8)  ::  Cd_min =  1.05                  ! Drag Coefficient Minimum Col de Porte
      real(kind=real8)  ::  cCdUns = -5.00                  ! Drag Coefficient Correction for Unstability
      real(kind=real8)  ::  RapCm0
 
 
! Aerodynamic Resistances
! =======================
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
 
! Surface Type
! ~~~~~~~~~~~~
        ist    =      isotSV(ikl,ikv)                ! Soil Type
        ist__s =  min(ist, 1)                        ! 1 => Soil
        ist__w =  1 - ist__s                         ! 1 => Water Body
 
! Drag and Aerodynamic Resistance
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        d_TaTs =      TaT_SV(ikl,ikv)-Tsrfsv(ikl,ikv)
        RapCm0 =  log(za__SV(ikl,ikv)/Z0mdSV(4          )) &
     &         /  log(za__SV(ikl,ikv)/Z0mdSV(ivgtSV(ikl,ikv)))
        RapCm0 =      RapCm0     *RapCm0             ! Neutral Drag Coefficient
                                                     ! Vegetation   Correction
        CD_m   =  max(Cd_min*RapCm0,                &! Actual  Drag Coefficient
     &                fCdCdP*RapCm0*VV__SV(ikl,ikv)  )  &!         for  Momentum
     &          *(1.+max(min(d_TaTs,zer0),cCdUns)   &! Unstability  Correction
     &                                   /cCdUns )  &!
     &          * 1.5                                !
        ram_sv(ikl,ikv) = rhT_SV(ikl,ikv) *CpdAir/CD_m   !
        rah_sv(ikl,ikv) = ram_sv(ikl,ikv)            !
 
 
! Turbulent Scales
! ================
 
! Friction Velocity                   u*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        uustar      =      VV__SV(ikl,ikv) / ram_sv(ikl,ikv)
        us__SV(ikl,ikv) = sqrt(uustar)
 
! Real    Temperature Turbulent Scale theta*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        uts_SV(ikl,ikv) =    - d_TaTs      / rah_sv(ikl,ikv)
        thstar      =      uts_SV(ikl,ikv) / us__SV(ikl,ikv)
 
! Specific Humidity   Turbulent Scale qq*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        qqstar      =      uqs_SV(ikl,ikv) / us__SV(ikl,ikv)
 
! Virtual Temperature Turbulent Scale thetav*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        thstarv     = thstar      + TaT_SV(ikl,ikv) *(0.608*qqstar)
        thstars     =     sign(un_1,thstarv)
        thstara     =      abs(     thstarv)
        thstarv     =      max(eps6,thstara)    *thstars
 
! Monin Obukhov Scale Height
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        LMO_SV(ikl,ikv) = TaT_SV(ikl,ikv) * uustar                     &
     &              /(vonKrm      * Grav_F     * thstarv)
        zeta        = za__SV(ikl,ikv) / LMO_SV(ikl,ikv)
        zeta_S      =          sign(un_1  ,zeta)
        zeta_A      =           abs(       zeta)
        zeta        = zeta_S  * max(eps6  ,zeta_A)
        LMO_SV(ikl,ikv) = za__SV(ikl,ikv) / zeta
 
      END DO
      END DO
 
 
      return
      end subroutine ColPrt_SBL
 
 
 
      subroutine SISVATeSBL
 
!--------------------------------------------------------------------------+
!   MAR          SISVATeSBL                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVATeSBL generates Surface Boundary Layers Properties   |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns                          |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   za__SV   : Surface Boundary Layer (SBL) Height          [m] |
!     ^^^^^    VV__SV   :(SBL Top)   Wind Velocity                   [m/s] |
!              TaT_SV   : SBL Top    Temperature                       [K] |
!              qsnoSV   : SBL Mean   Snow      Content             [kg/kg] |
!              uqs_SV   : Specific   Humidity  Turbulent Flux        [m/s] |
!              usthSV   : Blowing Snow Erosion   Threshold           [m/s] |
!              Z0m_SV   : Momentum     Roughness Length                [m] |
!              Z0h_SV   : Heat         Roughness Length                [m] |
!              Tsrfsv   : Surface    Temperature                       [K] |
!              sqrCm0   : Contribution of Z0m to Neutral Drag Coefficient  |
!              sqrCh0   : Contribution of Z0h to Neutral Drag Coefficient  |
!                                                                          |
!     INPUT /  LMO_SV   : Monin-Obukhov       Scale                    [m] |
!     OUTPUT:  us__SV   : Friction  Velocity                         [m/s] |
!     ^^^^^^   uts_SV   : Temperature         Turbulent Flux       [K.m/s] |
!              uss_SV   : Blowing Snow        Turbulent Flux         [m/s] |
!                                                                          |
!     OUTPUT:  hSalSV   : Saltating Layer Height                       [m] |
!     ^^^^^^   qSalSV   : Saltating Snow  Concentration            [kg/kg] |
!              ram_sv   : Aerodynamic Resistance for Momentum        [s/m] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD                                      |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^                                      |
!     #AE: TURBULENCE: Aerosols Erosion / Turbulent Diffusion Coeff.       |
!                                                                          |
!     #AW  TURBULENCE: Wind Time Mean (BOX Moving Average)                 |
!     #AH  TURBULENCE: Ta-T Time Mean (BOX Moving Average)                 |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #ZX  TURBULENCE: Strong Stability Limit (King    et al. 1996)        |
!     #zx  TURBULENCE: Strong Stability Limit (Mahalov et al. 2004)        |
!     #AX  TURBULENCE: recurrence                                          |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # stdout               | #ss: OUTPUT of Blowing Snow Variables         |
!                          |      unit  6, SubRoutine  SISVATeSBL **ONLY** |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_dat
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
      use Mod_SISVAT_ctr
! #AW use Mod_SISVAT_xAW
! #AH use Mod_SISVAT_xAH
      use Mod_SISVATLSBL
 
 
      IMPLICIT NONE
 
 
! Internal Variables
! ==================
 
! V,  dT(a-s)    Time Moving Averages
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer                                      ::  ikl,ikv   ,icount
! #AE integer                                      ::  nit    =   5    ! us(is0,uth) recursivity: Nb Iterations
! #AE integer                                      ::  iit
 
      real(kind=real8)                             ::  VVa_OK          ! effective SBL wind speed
      real(kind=real8)                             ::  Theta0 = 288.0  ! Potential Reference Temperature
 
!     real(kind=real8)                             ::  LMOsgn          ! Monin-Obukhov Scale Sign
!     real(kind=real8)                             ::  LMOabs          ! Monin-Obukhov Scale Abs.Value
 
      real(kind=real8)                             ::  uustar,thstar   !
      real(kind=real8)                             ::  qqstar,ssstar   !
      real(kind=real8)                             ::  thstarv,thstars !
      real(kind=real8)                             ::  thstara         !
      real(kind=real8)                             ::  zetam ,zetah    !
      real(kind=real8)                             ::  zeta0m,zeta0h   !
      real(kind=real8)                             ::  psim_s,xpsimi   !
      real(kind=real8)                             ::  psim_i,psim_z   !
! #AE real(kind=real8)                             ::  psis_s,psis_z   !
! #AE real(kind=real8)                             ::  psis_0          !
      real(kind=real8)                             ::  psih_s,xpsihi   !
      real(kind=real8)                             ::  psih_i,psih_z   !
      real(kind=real8)                             ::  psim_0,psih_0   !
      real(kind=real8)                             ::  dustar,u0star   !
 
      real(kind=real8)                             ::  sss__F,sss__N   !
      real(kind=real8)                             ::  usuth0          !
! #AE real(kind=real8)                             ::  dusuth,signus   !
! #AE real(kind=real8)                             ::  sss__K,sss__G   !
! #AE real(kind=real8)                             ::  us_127,us_227   !
! #AE real(kind=real8)                             ::  us_327,us_427   !
! #AE real(kind=real8)                             ::  us_527          !
 
      real(kind=real8)                             ::  stab_s          !
      real(kind=real8)                             ::  zetMAX =  1.e6  ! Strong Stability Limit
      real(kind=real8)                             ::  coef_m = 20.    ! Stabil.Funct.for Moment.: unstab.coef.
      real(kind=real8)                             ::  coef_h = 15.    ! Stabil.Funct.for Heat:    unstab.coef.
! #AE real(kind=real8)                             ::  SblPom =  1.27  ! Lower Boundary Height Parameter for Suspension
                                                                       ! Pommeroy, Gray and Landine 1993, J. Hydrology, 144(8) p.169
!     real(kind=real8)                             ::  fac_Ri          !
!     real(kind=real8)                             ::  Kz_vun          !
 
! OUTPUT of Snow Erosion Turbulence
! #b1 real(kind=real8)                             ::  W_pLMO          ! Pseudo Obukhov Length  (WRITE)
! #b1 real(kind=real8)                             ::  W_psim          ! Pseudo psim(z)         (WRITE)
 
! OUTPUT of Snow Erosion Turbulence (2)
! #b2 real(kind=real8)                             ::  W_NUs1          ! Contrib to U* numerat.1(WRITE)
! #b2 real(kind=real8)                             ::  W_NUs2          ! Contrib to U* numerat.2(WRITE)
! #b2 real(kind=real8)                             ::  W_NUs3          ! Contrib to U* numerat.3(WRITE)
! #b2 real(kind=real8)                             ::  W_DUs1          ! Contrib to U* denomin.1(WRITE)
! #b2 real(kind=real8)                             ::  W_DUs2          ! Contrib to U* denomin.2(WRITE)
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
      allocate            ( VVaSBL(kcolp,mwp) )                        ! effective SBL wind speed
      allocate            ( dTa_Ts(kcolp,mwp) )                        ! effective SBL Temperature diff.
      allocate            ( LMOmom(kcolp,mwp) )                        ! Monin-Obukhov Scale Momentum
      allocate            ( CDm(kcolp,mwp)    )                        ! Drag Coefficient, Momentum
      allocate            ( CDs(kcolp,mwp)    )                        ! Drag Coefficient, Blown **
      allocate            (rCDs(kcolp,mwp)    )                        ! Drag Coefficient, Blown **
      allocate            ( CDh(kcolp,mwp)    )                        ! Drag Coefficient, Scalar
      allocate            ( Richar(kcolp,mwp) )                        ! Richardson Number
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! Internal DATA
! =============
 
! #zx                       zetMAX = 1.e0                 ! Strong Stability Limit
                                                          !(Mahalov et al. 2004, GRL  31 2004GL021055)
! #ZX                       zetMAX = 4.28                 ! Strong Stability Limit
                                                          !(King    et al. 1996, JGR 101(7) p.19121)
 
 
! Effective SBL variables
! =======================
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
        VVaSBL(ikl,ikv)   = VV__SV(ikl,ikv)
! #AW   VVaSBL(ikl,ikv)   = VVmmem(ikl,ikv)
        dTa_Ts(ikl,ikv)   = TaT_SV(ikl,ikv)-Tsrfsv(ikl,ikv)
! #AH   dTa_Ts(ikl,ikv)   = dTmmem(ikl,ikv)
      END DO
      END DO
 
 
! Convergence Criterion
! =====================
 
      icount = 0
 
! #AX 1        CONTINUE
      icount = icount + 1
      dustar = 0.
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
          u0star      = us__SV(ikl,ikv)
 
 
! Turbulent Scales from previous Time Step
! ----------------------------------------
 
          u0star      =      max(eps6,u0star)      ! Friction Velocity     u*
          uustar      = u0star      * u0star       ! Friction Velocity^2  uu*
          thstar      = uts_SV(ikl,ikv) / u0star   ! Temperature       theta*
          qqstar      = uqs_SV(ikl,ikv) / u0star   ! Specific Humidity    qq*
          ssstar      = uss_SV(ikl,ikv) / u0star   ! Blown    Snow        ss*
 
 
! Monin-Obukhov Stability Parameter for Momentum
! ----------------------------------------------
 
! Pseudo Virtual Temperature Turbulent Scale thetav*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          thstarv     = thstar      + Theta0      *(0.608*qqstar)      &
     &                             /(1.+0.608*QaT_SV(ikl,ikv)-qsnoSV(ikl,ikv))
          thstars     =     sign(un_1,thstarv)
          thstara     =      abs(     thstarv)
          thstarv     =      max(eps6,thstara)*thstars
 
! Pseudo Obukhov Length Scale        (Gall?e et al., 2001 BLM 99, (A2) p.17)
! Full   Obukhov Length Scale        (when Blowing * is ##NOT## switched ON)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          LMO_SV(ikl,ikv) = Theta0      * max(eps6,uustar)             &
     &                /(vonKrm      * Grav_F  *thstarv)
 
! OUTPUT of Snow Erosion Turbulence
! #b1     W_pLMO      = LMO_SV(ikl,ikv)
 
          zetah       = za__SV(ikl,ikv) / LMO_SV(ikl,ikv)
          zetam       =           min(zetMAX,zetah)! Strong Stability Limit
                                                   !(Mahalov et al. 2004
                                                   ! GRL 31 2004GL021055)
          LMOmom(ikl,ikv) = za__SV(ikl,ikv) /(max(eps6,abs(zetam))     &
     &                              *sign(un_1,    zetam ))
          zeta0m      = Z0m_SV(ikl,ikv) / LMOmom(ikl,ikv)
          zeta0h      = Z0h_SV(ikl,ikv) / LMO_SV(ikl,ikv)
 
! Momentum Pseudo Stability Function (Gall?e et al. 2001, BLM 99, (11) p. 7)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          stab_s      =  max(zer0,sign(un_1,zetam))
 
          psim_s      =  -A_Stab *zetam
          xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zetam)))
          psim_i      =   2. *log(half*(un_1+xpsimi))                  &
     &                       +log(half*(un_1+xpsimi*xpsimi))           &
     &                   -2.*atan(xpsimi)   +half*piNmbr
          psim_z      =    stab_s*psim_s+(1.-stab_s)*psim_i
 
! OUTPUT of Snow Erosion Turbulence
! #b1     W_psim      =           psim_z
 
          psim_s      =  -A_Stab *zeta0m
          xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zeta0m)))
          psim_i      =   2. *log(half*(un_1+xpsimi))                  &
     &                       +log(half*(un_1+xpsimi*xpsimi))           &
     &                   -2.*atan(xpsimi)   +half*piNmbr
          psim_0      =    stab_s*psim_s+(1.-stab_s)*psim_i
 
! Virtual Temperature Turbulent Scale thetav*    (ss* impact included   )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ needed for new ss*)
! #AE     thstarv     = thstar      + Theta0      *(0.608*qqstar       &
! #AE&                                                   -ssstar       &
! #AE&                                             )                   &
! #AE&                             /(1.+0.608*QaT_SV(ikl,ikv)-qsnoSV(ikl,ikv))
! #AE     thstars     =     sign(un_1,thstarv)
! #AE     thstara     =      abs(     thstarv)
! #AE     thstarv     =      max(eps6,thstara)    *thstars
 
! Full   Obukhov Length Scale        (Gall?e et al. 2001, BLM 99, (A1) p.16)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #AE     LMO_SV(ikl,ikv) = Theta0      * us__SV(ikl,ikv)* us__SV(ikl,ikv) &
! #AE&                /(vonKrm      * Grav_F     * thstarv)
 
! #AE     zetah       = za__SV(ikl,ikv) / LMO_SV(ikl,ikv)
! #AE     zetam       =           min(zetMAX,zetah)! Strong Stability Limit
                                                   !(Mahalov et al. 2004
                                                   ! GRL 31 2004GL021055)
! #AE     LMOmom(ikl,ikv) = za__SV(ikl,ikv) /(max(eps6,abs(zetam))     &
! #AE&                              *sign(un_1,    zetam ))
! #AE     zeta0m      = Z0m_SV(ikl,ikv) / LMOmom(ikl,ikv)
 
! Snow Erosion    Stability Function (Gall?e et al. 2001, BLM 99, (11) p. 7)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #AE     stab_s      =  max(zer0,sign(un_1,zetam))
 
! #AE     psis_s      =  -AsStab *zetam
! #AE     xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zetam)))
! #AE     psim_i      =   2. *log(half*(un_1+xpsimi))                  &
! #AE&                       +log(half*(un_1+xpsimi*xpsimi))           &
! #AE&                   -2.*atan(xpsimi)   +half*piNmbr
! #AE     psis_z      =    stab_s*psis_s+(1.-stab_s)*psim_i
 
! #AE     psis_s      =  -AsStab *zeta0m
! #AE     xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zeta0m)))
! #AE     psim_i      =   2. *log(half*(un_1+xpsimi))                  &
! #AE&                       +log(half*(un_1+xpsimi*xpsimi))           &
! #AE&                   -2.*atan(xpsimi)   +half*piNmbr
! #AE     psis_0      =    stab_s*psis_s+(1.-stab_s)*psim_i
 
! Square Roots of the Drag Coefficient for Snow Erosion Turbulent Flux
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #AE     rCDmSV(ikl,ikv) = vonKrm/(sqrCm0(ikl,ikv)-psim_z+psim_0)
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Martin control : on remplace les "! #ss" par rien au début de la ligne
          IF (ii__AP(ikl).EQ.iwr_SV.AND.jj__AP(ikl).EQ.jwr_SV   .AND.  &
     &        ikv        .EQ.nwr_SV                          )     &
     &    write(6,6600)  Z0m_SV(ikl,ikv) , psim_z                      &
     &                  ,LMO_SV(ikl,ikv) , uustar                      &
     &                  ,sqrCm0(ikl,ikv) , psim_0                      &
     &                  ,LMOmom(ikl,ikv) , thstarv
      6600      format(/,' ** SISVATeSBL *0  '                         &
     &            ,'  Z0m_SV  = ',e12.4,'  psim_z  = ',e12.4           &
     &            ,'  LMO_SV  = ',e12.4,'  uustar  = ',e12.4           &
     &          ,/,'                   '                               &
     &            ,'  sqrCm0  = ',e12.4,'  psim_0  = ',e12.4           &
     &            ,'  LMOmom  = ',e12.4,'  thstarv = ',e12.4)
! Martin control : on remplace les "! #ss" par rien au début de la ligne
 
 
! Momentum            Turbulent Scale  u*
! ---------------------------------------
 
! Momentum            Turbulent Scale  u*          in case of NO Blow. Snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VVa_OK      =  max(0.000001,       VVaSBL(ikl,ikv))
          sss__N      =  vonKrm      *       VVa_OK
          sss__F      = (sqrCm0(ikl,ikv) - psim_z + psim_0)
          usuth0      =  sss__N /sss__F                ! u* if NO Blow. Snow
 
! Momentum            Turbulent Scale  u*          in case of    Blow. Snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #AE     sss__G      =  0.27417     * Grav_F
 
! ______________               _____
! Newton-Raphson (! Iteration, BEGIN)
! ~~~~~~~~~~~~~~               ~~~~~
! #AE     DO iit=1,nit
! #AE     sss__K      =  Grav_F      * r_Stab * A_Stab *za__SV(ikl,ikv)&
! #AE&                                     *rCDmSV(ikl,ikv)*rCDmSV(ikl,ikv)&
! #AE&                           /(1.+0.608*QaT_SV(ikl,ikv)-qsnoSV(ikl,ikv))
! #AE     us_127      =  exp(    SblPom *log(us__SV(ikl,ikv)))
! #AE     us_227      =  us_127         *    us__SV(ikl,ikv)
! #AE     us_327      =  us_227         *    us__SV(ikl,ikv)
! #AE     us_427      =  us_327         *    us__SV(ikl,ikv)
! #AE     us_527      =  us_427         *    us__SV(ikl,ikv)
 
! #AE     us__SV(ikl,ikv) =  us__SV(ikl,ikv)                           &
! #AE&    - (  us_527     *sss__F     /sss__N                          &
! #AE&      -  us_427                                                  &
! #AE&      -  us_227     *qsnoSV(ikl,ikv)*sss__K                      &
! #AE&      + (us__SV(ikl,ikv)*us__SV(ikl,ikv)-usthSV(ikl,ikv)*usthSV(ikl,ikv))/sss__G)&
! #AE&     /(  us_427*5.27*sss__F     /sss__N                          &
! #AE&      -  us_327*4.27                                             &
! #AE&      -  us_127*2.27*qsnoSV(ikl,ikv)*sss__K                      &
! #AE&      +  us__SV(ikl,ikv)*2.0                                 /sss__G)
 
! #AE     us__SV(ikl,ikv)= min(us__SV(ikl,ikv),usuth0)
! #AE     us__SV(ikl,ikv)= max(us__SV(ikl,ikv),eps6  )
! #AE     rCDmSV(ikl,ikv)=     us__SV(ikl,ikv)/VVa_OK
! #Ae     sss__F     =     vonKrm     /rCDmSV(ikl,ikv)
! #AE     END DO
! ______________               ___
! Newton-Raphson (! Iteration, END  )
! ~~~~~~~~~~~~~~               ~~~
 
! #AE     us_127      =  exp(    SblPom *log(us__SV(ikl,ikv)))
! #AE     us_227      =  us_127         *    us__SV(ikl,ikv)
 
! Momentum            Turbulent Scale  u*: 0-Limit in case of no Blow. Snow
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #AE     dusuth      =  us__SV(ikl,ikv) - usthSV(ikl,ikv)   ! u* - uth*
! #AE     signus      =  max(sign(un_1,dusuth),zer0)     ! 1 <=> u* - uth* > 0
          us__SV(ikl,ikv) =                             &!
! #AE&                   us__SV(ikl,ikv)  *signus  +    &! u* (_BS)
     &                   usuth0                         &! u* (nBS)
! #AE&                            *(1.-signus)          &!
     &                +  0.
 
 
! Blowing Snow        Turbulent Scale ss*
! ---------------------------------------
 
! Blowing Snow Surface Boundary Condition
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #AE     hSalSV(ikl,ikv) = 8.436e-2  *exp(SblPom  *log(us__SV(ikl,ikv)))
! #AE     qSalSV(ikl,ikv) = (us__SV(ikl,ikv) * us__SV(ikl,ikv)         &
! #AE&                  -usthSV(ikl,ikv) * usthSV(ikl,ikv))*signus     &
! #AE&                / (sss__G      * us_227     )
 
! Blowing Snow Surface Boundary Condition (modification, .NOT. tested)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #ae     qSalSV(ikl,ikv) = (us__SV(ikl,ikv) * us__SV(ikl,ikv)         &
! #ae&                  -usthSV(ikl,ikv) * usthSV(ikl,ikv))            &
! #ae&                  *signus      * us__SV(ikl,ikv) *3.25           &
! #ae&                 /(hSalSV(ikl,ikv) * Grav_F           )
 
! #AE     ssstar      =  rCDmSV(ikl,ikv) *(qsnoSV(ikl,ikv) -qSalSV(ikl,ikv))   &
! #AE&                 * r_Stab
 
! #AE     uss_SV(ikl,ikv) =  min(zer0    , us__SV(ikl,ikv) *ssstar)
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #ss     IF (ii__AP(ikl).EQ.iwr_SV.AND.jj__AP(ikl).EQ.jwr_SV   .AND.  &
! #ss&        ikv        .EQ.nwr_SV                          )  THEN
! #ss     write(6,6000)  daHost     ,     icount     ,                 &
! #ss&                   us__SV(ikl,ikv),1.e3*hSalSV(ikl,ikv),         &
! #ss&              1.e3*Z0m_SV(ikl,ikv),                              &
! #ss&              1.e3*qsnoSV(ikl,ikv),1.e3*qSalSV(ikl,ikv)          &
! #ss&                  ,usthSV(ikl,ikv),     us__SV(ikl,ikv)-usthSV(ikl,ikv), &
! #ss&              1.e3*ssstar     ,1.e3*us__SV(ikl,ikv)*ssstar
! #ss 6000 format(a18,i3,6x,'u*   [m/s] =',f6.3,'   hSalt[mm]='  ,e9.3,&
! #ss&                  '   Z0m   [mm] =',f9.3,'   q   [g/kg] =',f9.3, &
! #ss&               /,91x,                    '   qSa [g/kg] =',f9.3, &
! #ss&               /,27x, 'ut*[m/s]='  ,e9.3,'   u*-ut*   ='  ,e9.3, &
! #ss&                  '   s*  [g/kg] =',f9.3,'   us* [mm/s] =',f9.3)
! #ss     END IF
 
 
! Virtual Temperature Turbulent Scale thetav*    (ss* impact included)
! --------------------------------------------------------------------
 
! #AE     thstarv     = thstar      + Theta0      *(0.608*qqstar       &
! #AE&                                                   -ssstar       &
! #AE&                                             )                   &
! #AE&                             /(1.+0.608*QaT_SV(ikl,ikv)-qsnoSV(ikl,ikv))
! #AE     thstars     =     sign(un_1,thstarv)
! #AE     thstara     =      abs(     thstarv)
! #AE     thstarv     =      max(eps6,thstara)    *thstars
 
 
! Full   Obukhov Length Scale (Gall?e et al., 2001, BLM 99, (A1) p.16)
! --------------------------------------------------------------------
 
! #AE     LMO_SV(ikl,ikv) = Theta0      * us__SV(ikl,ikv)* us__SV(ikl,ikv) &
! #AE&                /(vonKrm      * Grav_F     * thstarv)
 
! #AE     zetah       = za__SV(ikl,ikv) / LMO_SV(ikl,ikv)
! #AE     zetam       =           min(zetMAX,zetah)! Strong Stability Limit
                                                   !(Mahalov et al. 2004
                                                   ! GRL 31 2004GL021055)
! #AE     LMOmom(ikl,ikv) = za__SV(ikl,ikv) /(max(eps6,abs(zetam))     &
! #AE&                              *sign(un_1,    zetam ))
! #AE     zeta0m      = Z0m_SV(ikl,ikv) / LMOmom(ikl,ikv)
! #AE     zeta0h      = Z0h_SV(ikl,ikv) / LMO_SV(ikl,ikv)
 
! OUTPUT in SISVAT at specified i,j,k,n (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #ss     IF (ii__AP(ikl).EQ.iwr_SV.AND.jj__AP(ikl).EQ.jwr_SV   .AND.  &
! #ss&        ikv        .EQ.nwr_SV                          )  THEN
! #ss     write(6,6001)  LMO_SV(ikl,ikv)    ,    zetah
! #ss 6001      format(18x,9x,'LMO  [m]=',f9.1,'   zetah[-] =',f9.3)
! #ss     END IF
 
 
! Turbulent Scales
! ----------------
 
! Momentum Stability Function (Gall?e et al., 2001, BLM 99, (11) p. 7)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          stab_s      =  max(zer0,sign(un_1,zetam))
 
          psim_s      =  -A_Stab *zetam
          xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zetam)))
          psim_i      =   2. *log(half*(un_1+xpsimi))                  &
     &                       +log(half*(un_1+xpsimi*xpsimi))           &
     &                   -2.*atan(xpsimi)   +half*piNmbr
          psim_z      =    stab_s*psim_s+(1.-stab_s)*psim_i
 
          psim_s      =  -A_Stab *zeta0m
          xpsimi      = sqrt(sqrt(un_1-coef_m*min(zer0,zeta0m)))
          psim_i      =   2. *log(half*(un_1+xpsimi))                  &
     &                       +log(half*(un_1+xpsimi*xpsimi))           &
     &                   -2.*atan(xpsimi)   +half*piNmbr
          psim_0      =    stab_s*psim_s+(1.-stab_s)*psim_i
 
! Heat     Stability Function (Gallee et al., 2001, BLM 99, (11) p. 7)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          stab_s      =  max(zer0,sign(un_1,zetah))
 
          psih_s      =  -AhStab *zetah
          xpsihi      = sqrt(sqrt(un_1-coef_h*min(zer0,zetah)))
          psih_i      =   2. *log(half*(un_1+xpsihi))
          psih_z      =    stab_s*psih_s+(1.-stab_s)*psih_i
 
          psih_s      =  -AhStab *zeta0h
          xpsihi      = sqrt(sqrt(un_1-coef_h*min(zer0,zeta0h)))
          psih_i      =   2. *log(half*(un_1+xpsihi))
          psih_0      =    stab_s*psih_s+(1.-stab_s)*psih_i
 
! Square Roots of the Drag Coefficients
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          rCDhSV(ikl,ikv) = vonKrm/(sqrCh0(ikl,ikv)-psih_z+psih_0)
          rCDmSV(ikl,ikv) = vonKrm/(sqrCm0(ikl,ikv)-psim_z+psim_0)
 
! Drag Coefficients
! ~~~~~~~~~~~~~~~~~
          CDh(ikl,ikv)    = rCDmSV(ikl,ikv) * rCDhSV(ikl,ikv)
          CDm(ikl,ikv)    = rCDmSV(ikl,ikv) * rCDmSV(ikl,ikv)
 
! Real    Temperature Turbulent Scale theta*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          thstar      = rCDhSV(ikl,ikv) * dTa_Ts(ikl,ikv)
          uts_SV(ikl,ikv) = us__SV(ikl,ikv) * thstar
 
 
! Convergence Criterion
! =====================
 
          dustar      = max(dustar,abs(us__SV(ikl,ikv)-u0star))
 
! OUTPUT of Snow Erosion Turbulence
! #b1    IF (icount     .EQ.1  )                                    THEN
! #b1     write(6,6004)
! #b1 6004      format(122('-'))
! #b1    IF (mod(VVaSBL(ikl,ikv),4.).LT.0.1)                        THEN
! #b1     write(6,6003)
! #b1 6003      format('   V  Ta-Ts  Z0      It'                       &
! #b1&   ,' du*     u*    sss__F   CD       Qss       Qs*     '        &
! #b1&   ,' PseudOL Full-OL zetam   zetah   psim_z  psih_z')
! #b1     write(6,6004)
! #b1    END IF
! #b1    END IF
! #b1     write(6,6002) VVaSBL(ikl,ikv),dTa_Ts(ikl,ikv),Z0m_SV(ikl,ikv),icount &
! #b1&                 ,dustar     ,us__SV(ikl,ikv),sss__F             &
! #b1&                 ,   CDm(ikl,ikv),qSalSV(ikl,ikv),ssstar         &
! #b1&                 ,W_pLMO     ,LMO_SV(ikl,ikv)                    &
! #b1&                 ,zetam      ,zetah      ,W_psim     ,psih_z
! #b1 6002 format(2f6.1,f8.4,i3,f9.6,f6.3,f9.3,3f9.6,2f8.2,2f8.4,2f8.2)
 
! OUTPUT of Snow Erosion Turbulence (2): u*_AE
! #b2    IF (icount     .EQ.1  )                                    THEN
! #b2     write(6,6014)
! #b2 6014      format(100('-'))
! #b2    IF (mod(VVaSBL(ikl,ikv),4.).LT.0.1)                        THEN
! #b2     write(6,6013)
! #b2 6013     format('   V  Ta-Ts  Z0      It'                        &
! #b2&   ,' du*     u*    sss__F   W_NUs1   W_NUs2   W_NUs3      '     &
! #b2&   ,' W_DUs1     W_DUs2 ')
! #b2     write(6,6014)
! #b2    END IF
! #b2    END IF
! #b2     write(6,6012) VVaSBL(ikl,ikv),dTa_Ts(ikl,ikv),Z0m_SV(ikl,ikv),icount &
! #b2&                 ,dustar     ,us__SV(ikl,ikv),sss__F             &
! #b2&                 ,W_NUs1     ,W_NUs2     ,W_NUs3                 &
! #b2&                 ,W_DUs1     ,W_DUs2
! #b2 6012      format(2f6.1,f8.4,i3,f9.6,f6.3,f9.3,3f9.3,2f12.3)
 
        END DO
        END DO
 
! #AX IF (                     icount.lt. 3)                     GO TO 1
!     IF (dustar.gt.0.0001.AND.icount.lt. 6)                     GO TO 1
 
 
! Aerodynamic Resistances
! -----------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          ram_sv(ikl,ikv) = 1./(CDm(ikl,ikv)*max(VVaSBL(ikl,ikv),eps6))
          rah_sv(ikl,ikv) = 1./(CDh(ikl,ikv)*max(VVaSBL(ikl,ikv),eps6))
        END DO
        END DO
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
      IF (FlagDALLOC)                                             THEN !
 
      deallocate          ( VVaSBL )                                   ! effective SBL wind speed
      deallocate          ( dTa_Ts )                                   ! effective SBL Temperature diff.
      deallocate          ( LMOmom )                                   ! Monin-Obukhov Scale Momentum
      deallocate          ( CDm    )                                   ! Drag Coefficient, Momentum
      deallocate          ( CDs    )                                   ! Drag Coefficient, Blown **
      deallocate          (rCDs    )                                   ! Drag Coefficient, Blown **
      deallocate          ( CDh    )                                   ! Drag Coefficient, Scalar
      deallocate          ( Richar )                                   ! Richardson Number
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SISVATeSBL
 
 
 
      subroutine SISVAT_SBL
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_SBL                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_SBL generates Surface Boundary Layers Properties   |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns                          |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   za__SV   : Surface Boundary Layer (SBL) Height          [m] |
!     ^^^^^    VV__SV   :(SBL Top)   Wind Velocity                   [m/s] |
!              TaT_SV   : SBL Top    Temperature                       [K] |
!              uqs_SV   : Specific   Humidity  Turbulent Flux        [m/s] |
!              Z0m_SV   : Momentum   Roughness Length                  [m] |
!              Z0h_SV   : Heat       Roughness Length                  [m] |
!              Tsrfsv   : Surface    Temperature                       [K] |
!              sqrCm0   : Contribution of Z0m to Neutral Drag Coefficient  |
!              sqrCh0   : Contribution of Z0h to Neutral Drag Coefficient  |
!                                                                          |
!     INPUT /  LMO_SV   : Monin-Obukhov       Scale                    [m] |
!     OUTPUT:  us__SV   : Friction  Velocity                         [m/s] |
!     ^^^^^^   uts_SV   : Temperature         Turbulent Flux       [K.m/s] |
!                                                                          |
!     OUTPUT:  Fh__sv   : Stability Function                           [-] |
!     ^^^^^^   dFh_sv   : Stability Function (Derivative)              [-] |
!              ram_sv   : Aerodynamic Resistance for Momentum        [s/m] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!                                                                          |
!     WARNING: SISVAT_SBL blows up for too small z0m values & large z_SBL  |
!     ^^^^^^^                      (z0m = 1.8e-6 m for z_SBL = 20 m)       |
!                                                                          |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # stdout               | #sb: OUTPUT/Verification: SISVAT_SBL          |
!                          |      unit  6, SubRoutine  SISVAT_SBL **ONLY** |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_dat
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
 
 
      IMPLICIT NONE
 
 
 
! Internal Variables
! ==================
 
      integer           ::  ikl,ikv    ,ist    ,ist__s ,ist__w
      real(kind=real8)  ::  CD_m_0 ,CD_h_0 ,ram0   ,rah0  ,rahMIN
      real(kind=real8)  ::  d_TaTs ,RiB__D ,RiBulk
      real(kind=real8)  ::  bmstab ,Am1_FU ,Am2_FU ,Fm_Uns
      real(kind=real8)  ::  bhstab ,Ah1_FU ,Ah2_FU ,Fh_Uns,dFh_Un
      real(kind=real8)  ::  Aux_FS ,FStabl ,dFSdRi ,Stabil,Fm_loc
      real(kind=real8)  ::  uustar ,thstar ,qqstar
      real(kind=real8)  ::  thstarv,thstars,thstara
      real(kind=real8)  ::  zeta   ,zeta_S ,zeta_A
 
      real(kind=real8)  ::  zetMAX = 4.28 ! Strong Stability Limit
!                                         !(King et al. 1996, JGR 101(7) p.19121)
 
 
! Aerodynamic Resistances
! =======================
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
 
! Surface Type
! ~~~~~~~~~~~~
        ist    =      isotSV(ikl,ikv)                ! Soil Type
        ist__s =  min(ist, 1)                        ! 1 => Soil
        ist__w =  1 - ist__s                         ! 1 => Water Body
 
! Neutral Parameters
! ~~~~~~~~~~~~~~~~~~
        CD_m_0 =  0.16/    (sqrCm0(ikl,ikv)*sqrCm0(ikl,ikv)) ! Neutral Drag Coeff.Mom.
        CD_h_0 =  0.16/    (sqrCm0(ikl,ikv)*sqrCh0(ikl,ikv)) ! Neutral Drag Coeff.Heat
        ram0   =  1.0 /    (CD_m_0     *VV__SV(ikl,ikv)) ! Neutral Aero Resis.Mom.
        rah0   =  1.0 /    (CD_h_0     *VV__SV(ikl,ikv)) ! Neutral Aero Resis.Heat
 
! Bulk Richardson Number
! ~~~~~~~~~~~~~~~~~~~~~~
        RiB__D =          VV__SV(ikl,ikv) *VV__SV(ikl,ikv)             &
     &                   *TaT_SV(ikl,ikv)
        d_TaTs =         (TaT_SV(ikl,ikv)- Tsrfsv(ikl,ikv))
        RiBulk =  Grav_F *za__SV(ikl,ikv)* d_TaTs                      &
     &          / RiB__D
 
! OUTPUT/Verification: SISVAT_SBL
! #sb       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND.&
! #sb&          ikv        .GE.nwr_SV)                                 &
! #sb&      write(6,6600) Tsrfsv(ikl,ikv),TaT_SV(ikl,ikv),VV__SV(ikl,ikv)  &
! #sb&                  , d_TaTs     ,RiBulk
! #sb 6600  format(/,'Tem(s,a), Wind  , d_TaTs, RiBulk = ',5e15.6)
 
! Unstable Case
! ~~~~~~~~~~~~~
        bmstab =  ist__s * (13.7 -0.34 /sqrt(CD_m_0)) &! Momentum
     &          + ist__w *   4.9                       !
        bmstab =  10.    *  bmstab   * CD_m_0         &!
     &              *sqrt(za__SV(ikl,ikv)/ Z0m_SV(ikl,ikv))!
        Am1_FU =  bmstab *    sqrt(abs(RiBulk))        !
        Am2_FU =  Am1_FU +1.0 +10.*abs(RiBulk)         !
        Fm_Uns = (Am1_FU +1.0)/        Am2_FU          !
 
! OUTPUT/Verification: SISVAT_SBL
! #sb       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND.&
! #sb&          ikv        .GE.nwr_SV)                                 &
! #sb&      write(6,6601) CD_m_0     ,Z0m_SV(ikl,ikv),bmstab           &
! #sb&                  , ist__s     ,ist__w
! #sb 6601  format(/,'CD_m_0  , Z0m_SV, bmstab, ist/sw = ',3e15.6,2i15)
 
        bhstab =  ist__s * ( 6.3 -0.18 /sqrt(CD_h_0)) &! Heat
     &          + ist__w *   2.6                       !
        bhstab =  10.    *  bhstab   * CD_h_0         &!
     &              *sqrt(za__SV(ikl,ikv)/ Z0h_SV(ikl,ikv))!
        Ah1_FU =  bhstab *    sqrt(abs(RiBulk))        !
        Ah2_FU =  Ah1_FU +1.0 +10.*abs(RiBulk)         !
        Fh_Uns = (Ah1_FU +1.0)/        Ah2_FU          !
        dFh_Un =((Ah1_FU +2.0)/(Ah2_FU*Ah2_FU)) * 5.   !
 
! Stable   Case
! ~~~~~~~~~~~~~
        Aux_FS =          1.0 + 5.*    RiBulk
        FStabl =                Aux_FS*Aux_FS
        dFSdRi =                Aux_FS          *10.
 
! Effective Stability Functions and Derivatives
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Stabil      =          sign(un_1,d_TaTs)
        Fm_loc      =  FStabl * max(zer0,Stabil)                       &
     &               - Fm_Uns * min(zer0,Stabil)
        Fh__sv(ikl,ikv) =  FStabl * max(zer0,Stabil)                   &
     &               - Fh_Uns * min(zer0,Stabil)
        dFh_sv(ikl,ikv) =  dFSdRi * max(zer0,Stabil)                   &
     &               - dFh_Un * min(zer0,Stabil)
 
! OUTPUT/Verification: SISVAT_SBL
! #sb       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND.&
! #sb&          ikv        .GE.nwr_SV)                                 &
! #sb&      write(6,6602) FStabl     ,Stabil                           &
! #sb&                   ,Fm_Uns     ,Fm_loc
! #sb 6602  format(/,'FStabl  , Stabil, Fm_Uns, Fm_loc = ',4e15.6)
 
! Aerodynamic Resistances
! ~~~~~~~~~~~~~~~~~~~~~~~
        ram_sv(ikl,ikv) =  ram0   *          Fm_loc
        rah_sv(ikl,ikv) =  rah0   *          Fh__sv(ikl,ikv)
        rahMIN   = max(rah_sv(ikl,ikv),  abs(d_TaTs)*60./za__SV(ikl,ikv))
                                               ! 60 for 30dgC within 1/2 hour
        dFh_sv(ikl,ikv) =  rah0   *          dFh_sv(ikl,ikv)           &
     &              *  rahMIN /          rah_sv(ikl,ikv)
        rah_sv(ikl,ikv) =  rahMIN
 
 
! Square Root Contributions to the Drag Coefficients
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rCDmSV(ikl,ikv) =  sqrt(ram_sv(ikl,ikv) *VV__SV(ikl,ikv))
        rCDmSV(ikl,ikv) =  1.     / max(eps6,rCDmSV(ikl,ikv))
        rCDhSV(ikl,ikv) =       rah_sv(ikl,ikv) *VV__SV(ikl,ikv)       &
     &                                  *rCDmSV(ikl,ikv)
        rCDhSV(ikl,ikv) = (1.     / max(eps6,rCDhSV(ikl,ikv)))
 
! OUTPUT/Verification: SISVAT_SBL
! #sb       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND.&
! #sb&          ikv        .GE.nwr_SV)                                 &
! #sb&      write(6,6603) ram_sv(ikl,ikv),rah_sv(ikl,ikv)              &
! #sb&                   ,rCDmSV(ikl,ikv),rCDhSV(ikl,ikv)
! #sb 6603  format(/,'AeR(m,h), rCD(m,h)               = ',4e15.6)
 
 
! Turbulent Scales
! ================
 
! Friction Velocity                   u*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        uustar      =      VV__SV(ikl,ikv) / ram_sv(ikl,ikv)
        us__SV(ikl,ikv) = sqrt(uustar)
 
! Real    Temperature Turbulent Scale theta*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        uts_SV(ikl,ikv) =      d_TaTs      / rah_sv(ikl,ikv)
        thstar      =      uts_SV(ikl,ikv) / us__SV(ikl,ikv)
 
! Specific Humidity   Turbulent Scale qq*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        qqstar      =      uqs_SV(ikl,ikv) / us__SV(ikl,ikv)
 
! Virtual Temperature Turbulent Scale thetav*
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        thstarv     = thstar      + TaT_SV(ikl,ikv) *(0.608*qqstar     &
     &                                             )
        thstars     =     sign(un_1,thstarv)
        thstara     =      abs(     thstarv)
        thstarv     =      max(eps6,thstara)    *thstars
 
! Monin Obukhov Scale Height
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        LMO_SV(ikl,ikv) = TaT_SV(ikl,ikv) * uustar                     &
     &              /(vonKrm      * Grav_F     * thstarv)
        zeta        = za__SV(ikl,ikv) / LMO_SV(ikl,ikv)
        zeta        =           min(zetMAX,zeta)  ! Strong Stability Limit
!                                                 ! King et al.   1996
!                                                 ! JGR 101(7) p.19121
        zeta_S      =          sign(un_1  ,zeta)
        zeta_A      =           abs(       zeta)
        zeta        = zeta_S  * max(eps6  ,zeta_A)
        LMO_SV(ikl,ikv) = za__SV(ikl,ikv) / zeta
 
! OUTPUT/Verification: SISVAT_SBL
! #sb       IF (ii__AP(ikl).EQ.iwr_SV .AND. jj__AP(ikl).EQ.jwr_SV .AND.&
! #sb&          ikv        .GE.nwr_SV)                                 &
! #sb&      write(6,6604) us__SV(ikl,ikv),uts_SV(ikl,ikv)              &
! #sb&                   ,LMO_SV(ikl,ikv),zeta
! #sb 6604  format(/,'***(m,h), LMO   , zeta           = ',4e15.6)
 
      END DO
      END DO
 
 
      return
      end subroutine SISVAT_SBL
 
 
 
      subroutine SISVAT_TVg(                                           &
! #e1&                     (ETVg_d                                     &
     &                     )
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_TVg                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_TVg computes the Canopy Energy Balance             |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   ivgtSV   = 0,...,12:   Vegetation Type                      |
!     ^^^^^               0:          Water, Solid or Liquid               |
!              isnoSV   = total Nb of Ice/Snow Layers                      |
!                                                                          |
!     INPUT:   sol_SV   : Downward Solar Radiation                  [W/m2] |
!     ^^^^^    IRd_SV   : Surface  Downward Longwave Radiation      [W/m2] |
!              TaT_SV   : SBL Top  Temperature                         [K] |
!              rhT_SV   : SBL Top  Air  Density                    [kg/m3] |
!              QaT_SV   : SBL Top  Specific  Humidity              [kg/kg] |
!              psivSV   : Leaf     Water     Potential                 [m] |
!              IRs_SV   : Soil     IR Flux  (previous time step)    [W/m2] |
!              dt__SV   : Time     Step                                [s] |
!                                                                          |
!              SoCasv   : Absorbed Solar Radiation by Canopy (Normaliz)[-] |
!              tau_sv   : Fraction of Radiation transmitted by Canopy  [-] |
!              Evg_sv   : Soil+Vegetation Emissivity                   [-] |
!              Eso_sv   : Soil+Snow       Emissivity                   [-] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!              Sigmsv   : Canopy Ventilation Factor                    [-] |
!              LAI_sv   : Leaf Area  Index                             [-] |
!              LAIesv   : Leaf Area  Index (effective / transpiration) [-] |
!              glf_sv   : Green Leaf Fraction of NOT fallen Leaves     [-] |
!              rrMxsv   : Canopy Maximum Intercepted Rain          [kg/m2] |
!                                                                          |
!     INPUT /  TvegSV   : Canopy   Temperature                         [K] |
!     OUTPUT:  rrCaSV   : Canopy     Water     Content             [kg/m2] |
!     ^^^^^^                                                               |
!                                                                          |
!     OUTPUT:  IRv_sv   : Vegetation IR Flux                        [W/m2] |
!     ^^^^^^   HSv_sv   : Sensible Heat Flux                        [W/m2] |
!              HLv_sv   : Latent   Heat Flux                        [W/m2] |
!              Evp_sv   : Evaporation                              [kg/m2] |
!              EvT_sv   : Evapotranspiration                       [kg/m2] |
!              ETVg_d   : Vegetation  Energy Power Forcing          [W/m2] |
!                                                                          |
!     Internal Variables:                                                  |
!     ^^^^^^^^^^^^^^^^^^                                                   |
!                                                                          |
!     METHOD: The Newton-Raphson Scheme is preferable                      |
!     ^^^^^^  when computing over a long time step the heat content        |
!             of a medium having a very small or zero heat capacity.       |
!             This is to handle strong non linearities arising             |
!             in conjunction with rapid temperature variations.            |
!                                                                          |
!     REFERENCE: DR97: Koen de Ridder thesis, UCL, 1997                    |
!     ^^^^^^^^^                                                            |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #NN: Newton-Raphson Increment not added in last Iteration            |
!     #nc: OUTPUT Preparation for Stand Alone NetCDF File                  |
!                                                                          |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_ctr
      use Mod_SISVAT_dat
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
      use Mod_SISVAT_flx
 
 
 
! Internal Variables
! ==================
 
      use Mod_SISVATLTVg
 
 
      IMPLICIT NONE
 
 
! OUTPUT
! ------
 
 
      integer                                      ::  ikl,ikv         ! Grid Point Index
      integer                                      ::  nitmax =  5     ! Maximum  Iterations Number
      integer                                      ::  nit             !          Iterations Counter
      real(kind=real8)                             ::  d_Tveg          ! Canopy Temperat. Increment
      real(kind=real8)                             ::  dTvMAX =  5.    ! Canopy Temperat. Increment MAX
      real(kind=real8)                             ::  dHvdTv          ! Derivativ.of Canopy Energ.Budg.
      real(kind=real8)                             ::  Hv_Tv0          ! Imbalance of Canopy Energ.Budg.
      real(kind=real8)                             ::  Hv_MAX          ! MAX Imbal.of Canopy Energ.Budg.
      real(kind=real8)                             ::  Hv_MIN =  0.1   ! MIN Imbal.of Canopy Energ.Budg.
      real(kind=real8)                             ::  Hswich          ! Newton-Raphson         Switch
      real(kind=real8)                             ::  tau_Ca          ! Canopy IR Radiation Absorption
      real(kind=real8)                             ::  IR_net          ! InfraRed  NET(t)
      real(kind=real8)                             ::  EvFrac          ! Condensat./Transpirat. Switch
      real(kind=real8)                             ::  SnoMsk =  0.0   ! Canopy Snow            Switch
      real(kind=real8)                             ::  den_qs,arg_qs   !
!     real(kind=real8)                             ::  esat_i          ! Saturation Vapor Pressure       [hPa]
      real(kind=real8)                             ::  qsatvg          ! Canopy Saturat. Spec. Humidity
      real(kind=real8)                             ::  dqs_dT          ! d(qsatvg)/dTv
      real(kind=real8)                             ::  FacEvp,FacEvT   !
      real(kind=real8)                             ::  Fac_Ev          ! Evapo(transpi)ration Factor
      real(kind=real8)                             ::  F_Stom          ! Funct.  (Leaf Water Potential)
      real(kind=real8)                             ::  R0Stom          ! Minimum Stomatal Resistance
      real(kind=real8)                             ::  R_Stom          !         Stomatal Resistance
      real(kind=real8)                             ::  LAI_OK          ! 1. ==>  Leaves   exist
      real(kind=real8)                             ::  rrCaOK,snCaOK   !
      real(kind=real8)                             ::  dEvpOK          ! Positive Definiteness Correct.
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
! #e1 allocate            ( ETVg_d(kcolp,mwp) )                        ! VegetationPower, Forcing
      allocate            ( Tveg_0(kcolp,mwp) )                        ! Canopy Temperature, Previous t
      allocate            ( dIRdTv(kcolp,mwp) )                        ! InfraRed  NET(t), Derivative(t)
      allocate            ( dHSdTv(kcolp,mwp) )                        ! Sensible Heat FL. Derivative(t)
      allocate            ( dHLdTv(kcolp,mwp) )                        ! Latent   Heat FL. Derivative(t)
      allocate            ( dEvpdT(kcolp,mwp) )                        ! Evapo(transpi)ration Derivative
      allocate            ( dEvTdT(kcolp,mwp) )                        ! Evapo(transpi)ration Derivative
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! Newton-Raphson Scheme
! =====================
 
      nit    = 0
  101 CONTINUE
      nit    = nit + 1
      HV_MAX = 0.
 
 
! Temperature of the Previous Time Step
! -------------------------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          Tveg_0(ikl,ikv) = TvegSV(ikl,ikv)
 
 
! IR    Radiation Absorption
! --------------------------
 
          tau_Ca = 1.  - tau_sv(ikl,ikv)                    ! Canopy Absorption
          IRv_sv(ikl,ikv) = -2.0   *Evg_sv(ikl,ikv)  *StefBo   &!
     &                         *TvegSV(ikl,ikv)  *TvegSV(ikl,ikv)  &! Downward IR (OUT)
     &                         *TvegSV(ikl,ikv)  *TvegSV(ikl,ikv)   ! + Upward IR (OUT)
          dIRdTv(ikl,ikv) =                                &!
     &     -Evg_sv(ikl,ikv)*                               &!
     &                8.*StefBo*TvegSV(ikl,ikv)  *TvegSV(ikl,ikv)  &! Downward IR (OUT)
     &                         *TvegSV(ikl,ikv)             ! + Upward IR (OUT)
          IR_net =       tau_Ca                            &!
     &    *(Evg_sv(ikl,ikv)* IRd_SV(ikl,ikv)               &! Downward IR (IN)
     &     -             IRs_SV(ikl,ikv)                   &!   Upward IR (IN)
     &     +             IRv_sv(ikl,ikv))                   !          IR (OUT)
 
 
! Sensible Heat Flux
! ------------------
 
          dHSdTv(ikl,ikv) = rhT_SV(ikl,ikv)* Sigmsv(ikl,ikv) *CpdAir   &! Derivative, t(n)
     &                / rah_sv(ikl,ikv)                     !
          HSv_sv(ikl,ikv) = dHSdTv(ikl,ikv)                &! Value,      t(n)
     &                *(TaT_SV(ikl,ikv)-TvegSV(ikl,ikv))    !
 
 
! Latent   Heat Flux
! ------------------
 
! Canopy Saturation Specific Humidity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!       IF      (DeRidder)                                THEN !
          den_qs      =         TvegSV(ikl,ikv)    - 35.8      !
          arg_qs      = 17.27 *(TvegSV(ikl,ikv)    -273.16)   &!
     &                                   / den_qs              !
          qsatvg      = .0038 *        exp(arg_qs)     *0.875  !  0.875 = Tuning Hapex-Sahel
          dqs_dT      = qsatvg     * 4099.2   /(den_qs *den_qs)!
!       ELSE IF (Dudhia_MAR)                              THEN !
!         esat_i      = 6.107                                 &!
!    &     *exp(ExpIsv*(un_1/WatIsv -un_1/TvegSV(ikl)    ))    !
!         qsatvg      = 0.622    *     esat_i                 &!
!    &      / (10.*pkPaSV(ikl) - 0.378*esat_i)                 !
!         dqs_dT      = qsatvg                                &!
!    &     *(1.0+0.6077*qsatvg     )                          &!
!    &     *    ExpIsv/(TvegSV(ikl)    *TvegSV(ikl)    )       !
!       END IF
 
! Canopy Stomatal Resistance
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          R0Stom = min(         StodSV(ivgtSV(ikl,ikv))     &!
     &                /max(eps6,glf_sv(       ikl,ikv)),StxdSV)  ! Min Stomatal R.
          F_Stom = pscdSV / max(pscdSV-psivSV(ikl,ikv) ,eps6)! F(Leaf Wat.Pot.)
                                                             ! DR97, eqn. 3.22
          R_Stom =(R0Stom / max(LAIesv(ikl,ikv), R0Stom/StxdSV))&! Can.Stomatal R.
     &           * F_Stom                                    ! DR97, eqn. 3.21
 
! Evaporation / Evapotranspiration
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          SnoMsk = max(zer0, sign(un_1,snCaSV(ikl,ikv)-eps_21))  !
          EvFrac = max(zer0, sign(un_1,QaT_SV(ikl,ikv)-qsatvg))  ! Condensation/
          EvFrac = EvFrac                                   &! Transpiration
     &       + (1.-EvFrac)*((1-SnoMsk)*         rrCaSV(ikl,ikv) &!        Switch
     &                                         /rrMxsv(ikl,ikv) &!
     &                      +  SnoMsk *min(un_1,snCaSV(ikl,ikv) &!
     &                                         /rrMxsv(ikl,ikv)))!
          Fac_Ev = rhT_SV(ikl,ikv) *Sigmsv(ikl,ikv)          ! Idem,  Factor
          FacEvp = Fac_Ev      *EvFrac     / rah_sv(ikl,ikv) !
          Evp_sv(ikl,ikv) = FacEvp*(qsatvg     - QaT_SV(ikl,ikv))! Evaporation
          dEvpdT(ikl,ikv) = FacEvp* dqs_dT                   ! Evp Derivative
          FacEvt = Fac_Ev * (1.-EvFrac)    /(rah_sv(ikl,ikv)&!
     &                              +R_Stom *Sigmsv(ikl,ikv))!
          EvT_sv(ikl,ikv) = FacEvt*(qsatvg     - QaT_SV(ikl,ikv))! EvapoTranspir.
          dEvTdT(ikl,ikv) = FacEvt* dqs_dT                   ! EvT Derivative
          HLv_sv(ikl,ikv) =-LhvH2O*(Evp_sv(ikl,ikv)+ EvT_sv(ikl,ikv))   &! Latent   Heat
     &                 -LhfH2O* Evp_sv(ikl,ikv)* SnoMsk      !(Subli.Contrib.)
          dHLdTv(ikl,ikv) = LhvH2O*(dEvpdT(ikl,ikv)+ dEvTdT(ikl,ikv))   &!
     &                 +LhfH2O* dEvpdT(ikl,ikv)* SnoMsk      !
 
 
! Imbalance  of the Canopy  Energy Budget
! ---------------------------------------
 
          LAI_OK = max(zer0,                                &! NO Budget if
     &                 sign(un_1, LAI_sv(ikl,ikv)-eps_21))   ! no Leaves
          Hv_Tv0 =    (  SoCasv(ikl,ikv)         *sol_SV(ikl,ikv)   &! Absorbed Solar
     &                 + IR_net                             &! NET      IR
     &                 + HSv_sv(ikl,ikv)                    &! Sensible Heat
     &                 + HLv_sv(ikl,ikv)                    &! Latent   Heat
     &                ) *LAI_OK                              !
 
! OUTPUT/Verification: Energy/Water Budget
! #e1     ETVg_d(ikl,ikv) =  Hv_Tv0                          ! Veg.Energ.Bal.
 
          Hswich      =          1.00
! #NN     Hswich      = max(zer0,                           &! Newton-Raphson
! #NN&                      sign(un_1,    abs(Hv_Tv0     )  &!         Switch
! #NN&                                       -Hv_MIN      )) !
 
 
! Derivative of the Canopy  Energy Budget
! ---------------------------------------
 
          dHvdTv    =   dIRdTv(ikl,ikv) * max(eps_21,tau_Ca)           &
     &                - dHSdTv(ikl,ikv)                                &
     &                - dHLdTv(ikl,ikv)
 
 
! Update Canopy and Surface/Canopy Temperatures
! ---------------------------------------------
 
          d_Tveg      = Hv_Tv0       / dHvdTv               !
          d_Tveg      =      sign(un_1,d_Tveg)             &! Increment
     &                       *min( abs(d_Tveg)     ,dTvMAX) ! Limitor
          TvegSV(ikl,ikv) = TvegSV(ikl,ikv)  - Hswich      *d_Tveg  ! Newton-Raphson
          Hv_MAX      = max(Hv_MAX,abs(Hv_Tv0     ))        !
 
 
! Update Vegetation Fluxes
! ------------------------
 
! #NN     IRv_sv(ikl,ikv) = IRv_sv(ikl,ikv)-dIRdTv(ikl,ikv)    *d_Tveg  ! Emitted  IR
! #NN     HSv_sv(ikl,ikv) = HSv_sv(ikl,ikv)+dHSdTv(ikl,ikv)    *d_Tveg  ! Sensible Heat
! #NN     Evp_sv(ikl,ikv) = Evp_sv(ikl,ikv)-dEvpdT(ikl,ikv)    *d_Tveg  ! Evapotranspir.
! #NN     EvT_sv(ikl,ikv) = EvT_sv(ikl,ikv)-dEvTdT(ikl,ikv)    *d_Tveg  ! Evapotranspir.
! #NN     HLv_sv(ikl,ikv) = HLv_sv(ikl,ikv)+dHLdTv(ikl,ikv)    *d_Tveg  ! Latent   Heat
 
          IRv_sv(ikl,ikv) = IRv_sv(ikl,ikv)                *LAI_OK
          HSv_sv(ikl,ikv) = HSv_sv(ikl,ikv)                *LAI_OK
          Evp_sv(ikl,ikv) = Evp_sv(ikl,ikv)                *LAI_OK
          EvT_sv(ikl,ikv) = EvT_sv(ikl,ikv)                *LAI_OK
          HLv_sv(ikl,ikv) = HLv_sv(ikl,ikv)                *LAI_OK
        END DO
        END DO
 
! #AX IF (                     nit.lt.nitmax)    GO TO 101
      IF (Hv_MAX.gt.Hv_MIN.and.nit.lt.nitmax)    GO TO 101
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          IRv_sv(ikl,ikv) = IRv_sv(ikl,ikv)                &! Emitted  IR
     &   +dIRdTv(ikl,ikv) *(TvegSV(ikl,ikv)-Tveg_0(ikl,ikv))!
          HSv_sv(ikl,ikv) = HSv_sv(ikl,ikv)                &! Sensible Heat
     &   -dHSdTv(ikl,ikv) *(TvegSV(ikl,ikv)-Tveg_0(ikl,ikv))!
          Evp_sv(ikl,ikv) = Evp_sv(ikl,ikv)                &! Evaporation
     &   +dEvpdT(ikl,ikv) *(TvegSV(ikl,ikv)-Tveg_0(ikl,ikv))!
          EvT_sv(ikl,ikv) = EvT_sv(ikl,ikv)                &! Transpiration
     &   +dEvTdT(ikl,ikv) *(TvegSV(ikl,ikv)-Tveg_0(ikl,ikv))!
          HLv_sv(ikl,ikv) = HLv_sv(ikl,ikv)                &! Latent   Heat
     &   -dHLdTv(ikl,ikv) *(TvegSV(ikl,ikv)-Tveg_0(ikl,ikv))!
 
! OUTPUT for Stand Alone NetCDF File
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          HLv_KL(ikl,ikv) = HLv_sv(ikl,ikv)
 
 
! Update Canopy Water Content
! ---------------------------
 
          rrCaSV(ikl,ikv) = rrCaSV(ikl,ikv)-(1.-SnoMsk)*Evp_sv(ikl,ikv)*dt__SV
          snCaSV(ikl,ikv) = snCaSV(ikl,ikv)-    SnoMsk *Evp_sv(ikl,ikv)*dt__SV
 
! Correction for Positive Definiteness (see WKarea/EvpVeg/EvpVeg.f)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          rrCaOK =  max(rrCaSV(ikl,ikv), 0.)
          snCaOK =  max(snCaSV(ikl,ikv), 0.)
          dEvpOK =     (rrCaOK-rrCaSV(ikl,ikv)                         &
     &                 +snCaOK-snCaSV(ikl,ikv))/dt__SV
 
          Evp_sv(ikl,ikv) = Evp_sv(ikl,ikv)       - dEvpOK  ! Evaporation
          HLv_sv(ikl,ikv) = HLv_sv(ikl,ikv)                &! Latent   Heat
     &    +(1.-SnoMsk)* LhvH2O            * dEvpOK         &!
     &    +    SnoMsk *(LhvH2O+LhfH2O)    * dEvpOK          !
 
          rrCaSV(ikl,ikv) = rrCaOK
          snCaSV(ikl,ikv) = snCaOK
 
        END DO
        END DO
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
      IF (FlagDALLOC)                                             THEN !
 
! #e1 deallocate          ( ETVg_d )                                   ! VegetationPower, Forcing
      deallocate          ( Tveg_0 )                                   ! Canopy Temperature, Previous t
      deallocate          ( dIRdTv )                                   ! InfraRed  NET(t), Derivative(t)
      deallocate          ( dHSdTv )                                   ! Sensible Heat FL. Derivative(t)
      deallocate          ( dHLdTv )                                   ! Latent   Heat FL. Derivative(t)
      deallocate          ( dEvpdT )                                   ! Evapo(transpi)ration Derivative
      deallocate          ( dEvTdT )                                   ! Evapo(transpi)ration Derivative
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SISVAT_TVg
 
 
 
      subroutine SISVAT_TSo(                                           &
! #e1&                     ,ETSo_0,ETSo_1,ETSo_d,kcolw                 &
     &                     )
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_TSo                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_TSo computes the Soil/Snow Energy Balance          |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   isotSV   = 0,...,11:   Soil       Type                      |
!     ^^^^^               0:          Water, Solid or Liquid               |
!              isnoSV   = total Nb of Ice/Snow Layers                      |
!              dQa_SV   = Limitation of  Water Vapor  Turbulent Flux       |
!                                                                          |
!     INPUT:   sol_SV   : Downward Solar Radiation                  [W/m2] |
!     ^^^^^    IRd_SV   : Surface Downward  Longwave   Radiation    [W/m2] |
!              za__SV   : SBL Top    Height                            [m] |
!              VV__SV   : SBL Top    Wind Speed                      [m/s] |
!              TaT_SV   : SBL Top    Temperature                       [K] |
!              ExnrSV   : Exner      Potential                         [-] |
!              rhT_SV   : SBL Top    Air  Density                  [kg/m3] |
!              QaT_SV   : SBL Top    Specific  Humidity            [kg/kg] |
!              LSdzsv   : Vertical   Discretization Factor             [-] |
!                       =    1. Soil                                       |
!                       = 1000. Ocean                                      |
!              dzsnSV   : Snow Layer Thickness                         [m] |
!              ro__SV   : Snow/Soil  Volumic Mass                  [kg/m3] |
!              eta_SV   : Soil Water Content                       [m3/m3] |
!              dt__SV   : Time Step                                    [s] |
!                                                                          |
!              SoSosv   : Absorbed Solar Radiation by Surfac.(Normaliz)[-] |
!              IRv_sv   : Vegetation  IR Radiation                  [W/m2] |
!              tau_sv   : Fraction of Radiation transmitted by Canopy  [-] |
!              Evg_sv   : Soil+Vegetation Emissivity                   [-] |
!              Eso_sv   : Soil+Snow       Emissivity                   [-] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!              Lx_H2O   : Latent Heat of Vaporization/Sublimation   [J/kg] |
!              Sigmsv   : Canopy Ventilation Factor                    [-] |
!              sEX_sv   : Verticaly Integrated Extinction Coefficient  [-] |
!                                                                          |
!     INPUT /  TsisSV   : Soil/Ice Temperatures (layers -nsoil,-nsoil+1 ,0)|
!     OUTPUT:           & Snow     Temperatures (layers  1,2,...,nsno) [K] |
!     ^^^^^^                                                               |
!                                                                          |
!     OUTPUT:  IRs_SV   : Soil      IR Radiation                    [W/m2] |
!     ^^^^^^   HSs_sv   : Sensible  Heat Flux                       [W/m2] |
!              HLs_sv   : Latent    Heat Flux                       [W/m2] |
!              ETSo_0   : Snow/Soil Energy Power, before Forcing    [W/m2] |
!              ETSo_1   : Snow/Soil Energy Power, after  Forcing    [W/m2] |
!              ETSo_d   : Snow/Soil Energy Power         Forcing    [W/m2] |
!                                                                          |
!     METHOD: NO   Skin Surface Temperature                                |
!     ^^^^^^  Semi-Implicit Crank Nicholson Scheme                         |
!                                                                          |
!     REFERENCE: DR97: Koen de Ridder thesis, UCL, 1997                    |
!     ^^^^^^^^^                                                            |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #WL: TURBULENCE: u*q* limited to SBL Saturat.Specif.Humid.           |
!     #TR: TURBULENCE: Richardson Number:    T Derivative is used          |
!     #TL: TURBULENCE: Latent     Heat Flux: T Derivative is used          |
!     #nc: OUTPUT Preparation for Stand Alone NetCDF File                  |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_iii_jjj_n     | #e1: OUTPUT/Verification: Energy Conservation |
!                          |                                               |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_ctr
      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
      use Mod_SISVAT_flx
 
 
 
! Internal Variables
! ==================
 
      use Mod_SISVATLTSo
 
 
      IMPLICIT NONE
 
 
      real(kind=real8)                               ::  deltak
      real(kind=real8)                               ::  Exp_SA, Imp_SA
      real(kind=real8)                               ::  ExpTOP, ImpTOP
      real(kind=real8)                               ::  ExpHSL, ImpHSL
      real(kind=real8)                               ::  epsi15= 1.0e-15
      integer                                        ::  is1   , is2
 
 
! OUTPUT/Verification: Energy/Water Budget
! #e1 real(kind=real8),              dimension(kcolw)::  ETSo_0        ! Soil/Snow Power, before Forcing
! #e1 real(kind=real8),              dimension(kcolw)::  ETSo_1        ! Soil/Snow Power, after  Forcing
! #e1 real(kind=real8),              dimension(kcolw)::  ETSo_d        ! Soil/Snow Power, Forcing
 
      integer                                        ::  ikl,ikv   ,isl!
      integer                                        ::  jsl   ,isn    !
      integer                                        ::  isHigh        ! Order of the tridiagonal Matrix
      integer                                        ::  ist__s,ist__w,ist             ! Soil/Water  Body Identifier
      integer                                        ::  islsgn        ! Soil/Snow Surfac.Identifier
 
      real(kind=real8)                               ::  eps__3= 1.e-3 ! Arbitrary    Low Number
      real(kind=real8)                               ::  etaMid,psiMid ! Layer Interface's Humidity
      real(kind=real8)                               ::  mu_eta        !     Soil thermal Conductivity
      real(kind=real8)                               ::  mu_exp=-0.4343! arg Soil thermal Conductivity
      real(kind=real8)                               ::  mu_min= 0.172 ! Min Soil thermal Conductivity
      real(kind=real8)                               ::  mu_max= 2.000 ! Max Soil thermal Conductivity
      real(kind=real8)                               ::  mu_aux        !     Snow thermal Conductivity
      real(kind=real8)                               ::  dTSurf        ! Previous Surface Temperature
      real(kind=real8)                               ::  den_qs,arg_qs ! Soil   Saturat. Spec. Humidity
      real(kind=real8)                               ::  esat_i        ! Saturation Vapor Pressure       [hPa]
      real(kind=real8)                               ::  etaSol        ! Soil Surface          Humidity
      real(kind=real8)                               ::  d__eta        ! Soil Surface Humidity Increm.
      real(kind=real8)                               ::  Elem_A        !
      real(kind=real8)                               ::  ElemaA,ElemsA !   Diagonal Coefficients
      real(kind=real8)                               ::  Elem_C        !
      real(kind=real8)                               ::  ElemaC,ElemsC !   Diagonal Coefficients
      real(kind=real8)                               ::  Ts_Min = 175. ! Snow MIN Temperature
      real(kind=real8)                               ::  Ts_Max = 300. ! Snow MIN Temperature Acceptable
                                                                       ! both including Snow Melt Energy
 
! OUTPUT/Verification: Energy/Water Budget
! #e1 real(kind=real8)                               ::  Exist0        ! Existing Layer Switch
 
      integer,parameter                              ::  nt_srf=10     !
      integer                                        ::  it_srf,itEuBk ! HL: Surface Scheme
      real(kind=real8)                               ::  agpsrf,xgpsrf !
      real(kind=real8)                               ::  dt_srf,dt_ver !
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
      allocate            ( Tsisva(kcolp,mwp,-nsoil:nsnow+mzp) )       !
      allocate            ( Fsisva(kcolp,mwp) )                        !
      allocate            ( dza__1(kcolp,mwp) )                        !
      allocate            ( mu_sno(kcolp,mwp) )                        !     Snow thermal Conductivity
      allocate            ( mu__dz(kcolp,mwp,-nsoil:nsnow+1)   )       ! mu_(eta,sno)   / dz
      allocate            ( dtC_sv(kcolp,mwp,-nsoil:nsnow)     )       ! dt      / C
      allocate            ( IRs__D(kcolp,mwp) )                        ! UpwardIR Previous Iter.Contr.
      allocate            ( dIRsdT(kcolp,mwp) )                        ! UpwardIR           T Derivat.
      allocate            ( f_HSHL(kcolp,mwp) )                        ! Factor common to HS and HL
      allocate            ( dRidTs(kcolp,mwp) )                        ! d(Rib)/d(Ts)
      allocate            ( HS___D(kcolp,mwp) )                        ! Sensible Heat Flux Atm.Contr.
      allocate            ( f___HL(kcolp,mwp) )                        !
      allocate            ( HL___D(kcolp,mwp) )                        ! Latent   Heat Flux Atm.Contr.
      allocate            ( TSurf0(kcolp,mwp) )                        ! Previous Surface Temperature
      allocate            ( qsatsg(kcolp,mwp) )                        ! Soil   Saturat. Spec. Humidity
      allocate            ( dqs_dT(kcolp,mwp) )                        ! d(qsatsg)/dTv
      allocate            ( Psi(   kcolp,mwp) )                        ! 1st Soil Layer Water Potential
      allocate            ( RHuSol(kcolp,mwp) )                        ! Soil Surface Relative Humidity
      allocate            ( RHu_av(kcolp,mwp) )                        ! Soil Surface Relative Humidity
      allocate            ( Diag_A(kcolp,mwp,-nsoil:nsnow+mzp) )       ! A Diagonal
      allocate            ( Diag_B(kcolp,mwp,-nsoil:nsnow+mzp) )       ! B Diagonal
      allocate            ( Diag_C(kcolp,mwp,-nsoil:nsnow+mzp) )       ! C Diagonal
      allocate            ( Term_D(kcolp,mwp,-nsoil:nsnow+mzp) )       !   Independant Term
      allocate            ( Aux__P(kcolp,mwp,-nsoil:nsnow+mzp) )       ! P Auxiliary Variable
      allocate            ( Aux__Q(kcolp,mwp,-nsoil:nsnow+mzp) )       ! Q Auxiliary Variable
      allocate            ( etaBAK(kcolp,mwp) )                        !
      allocate            ( etaNEW(kcolp,mwp) )                        !
      allocate            ( etEuBk(kcolp,mwp) )                        !
      allocate            ( fac_dt(kcolp,mwp) )                        !
      allocate            ( faceta(kcolp,mwp) )                        !
      allocate            ( PsiArg(kcolp,mwp) )                        !
      allocate            ( SHuSol(kcolp,mwp) )                        !
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! Implicitness of the numerical Scheme at Surface - Atmosphere Interface
! ======================================================================
 
!        vvvvvvvv
      IF (SVaKzT)                                                  THEN
       IF(mzp.LT.2) SVaUBC = .TRUE.
          Exp_SA =  Explic
          Imp_SA =  Implic
      ELSE
          Exp_SA =  1.00
          Imp_SA =  0.00
      END IF
          ExpHSL =  0.00
          ImpHSL =  1.00
 
 
 
 
! Heat Conduction Coefficient (zero in the Layers over the highest one)
! ===========================
!                               ---------------- isl    eta_SV, rho C (isl)
!
! Soil                          ++++++++++++++++        etaMid,    mu (isl)
! ----
!                               ---------------- isl-1  eta_SV, rho C (isl-1)
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
!         __________
          isl=-nsoil
 
          mu__dz(ikl,ikv,isl) = 0.
 
          dtC_sv(ikl,ikv,isl) = dtz_SV(isl)              &! dt / (dz X rho C)
     &                   /((rocsSV(isotSV(ikl,ikv))      &! [s / (m.J/m3/K)]
     &                     +rcwdSV*eta_SV(ikl,ikv,isl))  &!
     &                     *LSdzsv(ikl,ikv)            )  !
! #kv   END DO
! #kv   END DO
!         ______________
       DO isl=-nsoil+1,0
! #kv   DO ikl=1,kcolp
! #kv   DO ikv=1,mwp
          ist    =      isotSV(ikl,ikv)                   ! Soil Type
          ist__s =  min(ist, 1)                           ! 1 => Soil
          ist__w =  1 - ist__s                            ! 1 => Water Body
 
          etaMid = 0.5*(dz_dSV(isl-1)*eta_SV(ikl,ikv,isl-1)  &! eta at layers
     &                 +dz_dSV(isl)  *eta_SV(ikl,ikv,isl)  ) &!     interface
     &                 /dzmiSV(isl)                       ! LSdzsv implicit !
          etaMid =  max(etaMid,eps6)
          psiMid =      psidSV(ist)                      &!
     &                *(etadSV(ist)/etaMid)**bCHdSV(ist)  !
          mu_eta =      3.82      *(psiMid)**mu_exp       ! Soil Thermal
          mu_eta =  min(max(mu_eta, mu_min), mu_max)      ! Conductivity
                                                          ! DR97 eq.3.31
          mu_eta =  ist__s *mu_eta +ist__w * vK_dSV       ! Water Bodies
                                                          ! Correction
          mu__dz(ikl,ikv,isl) = mu_eta/(dzmiSV(isl)      &!
     &                             *LSdzsv(ikl,ikv))      !
 
          dtC_sv(ikl,ikv,isl) = dtz_SV(isl)              &! dt / (dz X rho C)
     &                   /((rocsSV(isotSV(ikl,ikv))      &!
     &                     +rcwdSV*eta_SV(ikl,ikv,isl))  &!
     &                     *LSdzsv(ikl,ikv)            )  !
! #kv   END DO
! #kv   END DO
       END DO
 
 
! Soil/Snow Interface
! -------------------
 
! Soil Contribution
! ^^^^^^^^^^^^^^^^^
!         _____
          isl=1
! #kv   DO ikl=1,kcolp
! #kv   DO ikv=1,mwp
          ist    =      isotSV(ikl,ikv)                   ! Soil Type
          ist__s =  min(ist, 1)                           ! 1 => Soil
          ist__w =  1 - ist__s                            ! 1 => Water Body
          psiMid =      psidSV(ist)                       ! Snow => Saturation
          mu_eta =      3.82      *(psiMid)**mu_exp       ! Soil Thermal
          mu_eta =  min(max(mu_eta, mu_min), mu_max)      ! Conductivity
                                                          ! DR97 eq.3.31
          mu_eta =  ist__s *mu_eta +ist__w * vK_dSV       ! Water Bodies
 
! Snow Contribution
! ^^^^^^^^^^^^^^^^^
          mu_sno(ikl,ikv) =  CdidSV                      &!
     &                 *(ro__SV(ikl,ikv,isl) /rhoWat) ** 1.88 !
          mu_sno(ikl,ikv) =          max(eps6,mu_sno(ikl,ikv))!
!         mu_sno :  Snow Heat Conductivity Coefficient [Wm/K]
!                   (Yen 1981, CRREL Rep., 81-10)
 
! Combined Heat Conductivity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          mu__dz(ikl,ikv,isl) = 2./(dzsnSV(ikl,ikv,isl  )&! Combined Heat
     &                         /mu_sno(ikl,ikv)          &! Conductivity
     &                         +LSdzsv(ikl,ikv)          &!
     &                         *dz_dSV(    isl-1)/mu_eta) ! Coefficient
 
! Inverted Heat Capacity
! ^^^^^^^^^^^^^^^^^^^^^^
          dtC_sv(ikl,ikv,isl) = dt__SV/max(eps6,         &! dt / (dz X rho C)
     &    dzsnSV(ikl,ikv,isl) * ro__SV(ikl,ikv,isl) *Cn_dSV)  !
 
 
! Snow
! ----
 
!           _____________________
         DO isl=1,isnoSV(ikl,ikv)
          ro__SV(ikl,ikv,isl) =                          &!
     &                   ro__SV(ikl,ikv ,isl)            &!
     &       * max(0,min(isnoSV(ikl,ikv)-isl+1,1))        !
         END DO
 
!           _____________________
         DO isl=1,isnoSV(ikl,ikv)
 
! Combined Heat Conductivity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          mu_aux      =  CdidSV                          &!
     &                 *(ro__SV(ikl,ikv,isl) /rhoWat) ** 1.88 !
          mu__dz(ikl,ikv,isl) =                          &!
     &      2.                        *mu_aux*mu_sno(ikl,ikv)&! Combined Heat
     &     /max(eps6,dzsnSV(ikl,ikv,isl  )*mu_sno(ikl,ikv)   &! Conductivity
     &              +dzsnSV(ikl,ikv,isl-1)*mu_aux     )   ! For upper Layer
          mu_sno(ikl,ikv)     =            mu_aux         !
 
! Inverted Heat Capacity
! ^^^^^^^^^^^^^^^^^^^^^^
          dtC_sv(ikl,ikv,isl) = dt__SV/max(eps__3,       &! dt / (dz X rho C)
     &    dzsnSV(ikl,ikv,isl) * ro__SV(ikl,ikv,isl) *Cn_dSV)  !
         END DO
! #kv   END DO
! #kv   END DO
 
 
! Uppermost Effective Layer: NO conduction
! ----------------------------------------
 
! #kv   DO ikl=1,kcolp
! #kv   DO ikv=1,mwp
          mu__dz(ikl,ikv,isnoSV(ikl,ikv)+1) = 0.0
 
 
 
! Soil - Ice - Snow - Vegetation - Atmosphere Variable
! ====================================================
 
          dza__1(ikl,ikv) =(zza_SV(ikl,ikv,min(mzp,2))   + zza_SV(ikl,ikv,1)) * 0.5
 
!              vvvvvvvv
!           IF (SVaKzT)                                             THEN
                Fsisva(ikl,ikv) = p0_kap *   roa_SV(ikl,ikv,1) * dza__1(ikl,ikv)   &
     &                     *  dtC_sv(ikl,ikv,isnoSV(  ikl,ikv))* CpdAir /  dt__SV
!           END IF
 
! Soil/Snow Interior
! ^^^^^^^^^^^^^^^^^^
!           __________
         DO isl=-nsoil        ,isnoSV(ikl,ikv)
          Tsisva(ikl,ikv,isl) =        TsisSV(ikl,ikv,isl)              &
     &                * p0_kap /   ExnrSV(ikl,ikv)
         END DO
! #kv   END DO
! #kv   END DO
 
! Atmosphere
! ^^^^^^^^^^
!         vvvvvvvv
!      IF (SVaKzT)                                                  THEN
!         _________
       DO isl=1,mzp
! #kv   DO ikl=1,kcolp
! #kv   DO ikv=1,mwp
          isn=isnoSV(ikl,ikv)+isl
          Tsisva(ikl,ikv,isn) =   pktaSV(ikl,ikv,isl) * Fsisva(ikl,ikv)
! #kv   END DO
! #kv   END DO
       END DO
 
!      END IF
 
 
! OUTPUT/Verification: Energy/Water Budget: Energy Budget (IN)
! #e1  DO ikl=1,kcolp
! #e1  DO ikv=1,mwp
! #e1     ETSo_0(ikl,ikv) = 0.
! #e1  END DO
! #e1  END DO
! #e1  DO ikl=1,kcolp
! #e1  DO ikv=1,mwp
! #e1   DO isl= -nsoil,isnoSV(ikl,ikv)
! #e1     Exist0      = isl -           isnoSV(ikl,ikv)
! #e1     Exist0      = 1.  - max(zer0,min(un_1,Exist0))
! #e1     ETSo_0(ikl,ikv) = ETSo_0(ikl,ikv)                             &
! #e1&                +(TsisSV(ikl,ikv,isl)-Tf_Sno)*Exist0              &
! #e1&                                 /dtC_sv(ikl,ikv,isl)
! #e1   END DO
! #e1  END DO
! #e1  END DO
 
 
! Tridiagonal Elimination: Set Up
! ===============================
 
! #kv   DO ikl=1,kcolp
! #kv   DO ikv=1,mwp
!           _______________________________
         DO isl= -nsoil+1,isnoSV(ikl,ikv)-1
          Elem_A          =  dtC_sv(ikl,ikv,isl)         *mu__dz(ikl,ikv,isl)
          Elem_C          =  dtC_sv(ikl,ikv,isl)         *mu__dz(ikl,ikv,isl+1)
          Diag_A(ikl,ikv,isl) = -Elem_A  *Implic
          Diag_C(ikl,ikv,isl) = -Elem_C  *Implic
          Diag_B(ikl,ikv,isl) =  1.0d+0  -Diag_A(ikl,ikv,isl)-Diag_C(ikl,ikv,isl)
          Term_D(ikl,ikv,isl) =  Explic *(Elem_A         *Tsisva(ikl,ikv,isl-1)  &
     &                               +Elem_C         *Tsisva(ikl,ikv,isl+1)) &
     &             +(1.0d+0 -Explic *(Elem_A+Elem_C))*Tsisva(ikl,ikv,isl)&
     &  + dtC_sv(ikl,ikv,isl)           * sol_SV(ikl,ikv)    *SoSosv(ikl,ikv)&
     &                                              *(sEX_sv(ikl,ikv,isl+1)  &
     &                                               -sEX_sv(ikl,ikv,isl  )) &
     &   *p0_kap          /  ExnrSV(ikl,ikv)
         END DO
 
! #kv   END DO
! #kv   END DO
 
! Soil  lowest Layer
! ^^^^^^^^^^^^^^^^^^
!          ___________
           isl= -nsoil
! #kv   DO ikl=1,kcolp
! #kv   DO ikv=1,mwp
          Elem_A          =  0.
          Elem_C          =  dtC_sv(ikl,ikv,isl)         *mu__dz(ikl,ikv,isl+1)
          Diag_A(ikl,ikv,isl) =  0.
          Diag_C(ikl,ikv,isl) = -Elem_C  *Implic
          Diag_B(ikl,ikv,isl) =  1.0d+0  -Diag_A(ikl,ikv,isl)-Diag_C(ikl,ikv,isl)
          Term_D(ikl,ikv,isl) =  Explic * Elem_C         *Tsisva(ikl,ikv,isl+1)  &
     &             +(1.0d+0 -Explic * Elem_C)        *Tsisva(ikl,ikv,isl)  ! & !
!    &  + dtC_sv(ikl,isl)           * sol_SV(ikl)    *SoSosv(ikl)        & ! .NOT. needed
!    &                                              *(sEX_sv(ikl,isl+1)  & ! since sEX_sv = 0.
!    &                                               -sEX_sv(ikl,isl  )) & !
!    &   *p0_kap          /  ExnrSV(ikl)
! #kv   END DO
! #kv   END DO
 
! Surface: UPwardIR Heat Flux
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #kv   DO ikl=1,kcolp
! #kv   DO ikv=1,mwp
!         _________________________
          isl         = isnoSV(ikl,ikv)
          dIRsdT(ikl,ikv) = Eso_sv(ikl,ikv)* StefBo          * 4. &! - d(IR)/d(T)
     &                             * TsisSV(ikl,ikv,isl)      &!
     &                             * TsisSV(ikl,ikv,isl)      &!
     &                             * TsisSV(ikl,ikv,isl)       !
          IRs__D(ikl,ikv) = dIRsdT(ikl,ikv)* TsisSV(ikl,ikv,isl) * 0.75!
 
! Surface: Richardson Number:   T Derivative
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #TR     dRidTs(ikl,ikv) =-Grav_F      *    za__SV(ikl,ikv)  &!
! #TR&                              *(1.-Sigmsv(ikl,ikv))     &!
! #TR&                /(TaT_SV(ikl,ikv) *    VV__SV(ikl,ikv)  &!
! #TR&                              *    VV__SV(ikl,ikv))      !
 
! Surface: Turbulent Heat Flux: Factors
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          f_HSHL(ikl,ikv) = rhT_SV(ikl,ikv) *(1.-Sigmsv(ikl,ikv)) &!#common factor
     &                              /    rah_sv(ikl,ikv)       ! to  HS, HL
          f___HL(ikl,ikv) = f_HSHL(ikl,ikv) *    Lx_H2O(ikl,ikv)
 
! Surface: Sensible  Heat Flux: T Derivative
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          dSdTSV(ikl,ikv) = f_HSHL(ikl,ikv) *    CpdAir       &!#- d(HS)/d(T)
! #TR&         *(1.0  -(TsisSV(ikl,ikv,isl) -TaT_SV(ikl,ikv)) &!#Richardson
! #TR&         * dRidTs(ikl,ikv)*dFh_sv(ikl,ikv)/rah_sv(ikl,ikv)) &! Nb. Correct.
     &         + 0.
          HS___D(ikl,ikv) = dSdTSV(ikl,ikv) *    TaT_SV(ikl,ikv)   !
 
! Surface: Latent    Heat Flux: Saturation Specific Humidity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!        IF      (DeRidder)                               THEN !
          den_qs      =         TsisSV(ikl,ikv,isl)- 35.8      !
          arg_qs      = 17.27 *(TsisSV(ikl,ikv,isl)-273.16)   &!
     &                                   / den_qs              !
          qsatsg(ikl,ikv) = .0038 *        exp(arg_qs)     *0.875  !  0.875 = Tuning Hapex-Sahel
          dqs_dT(ikl,ikv) = qsatsg(ikl,ikv)* 4099.2   /(den_qs *den_qs)!
!        ELSE IF (Dudhia_MAR)                             THEN !
!         esat_i      = 6.107                                 &!
!    &     *exp(ExpIsv*(un_1/WatIsv -un_1/TsisSV(ikl,isl)))    !
!         qsatsg(ikl) = 0.622    *     esat_i                 &!
!    &      / (10.*pkPaSV(ikl) - 0.378*esat_i)                 !
!         dqs_dT(ikl) = qsatsg(ikl)                           &!
!    &     *(1.0+0.6077*qsatsg(ikl))                          &!
!    &     *    ExpIsv/(TsisSV(ikl,isl)*TsisSV(ikl,isl))       !
!        END IF
 
          fac_dt(ikl,ikv) = f_HSHL(ikl,ikv)/(rhoWat   * dz_dSV(0)) !
! #kv   END DO
! #kv   END DO
 
! Surface: Latent    Heat Flux: Surface    Relative Humidity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              xgpsrf       =   1.05                            !
              agpsrf       = dt__SV*(   1.0-xgpsrf        )   &!
     &                             /(   1.0-xgpsrf**nt_srf)    !
              dt_srf       = agpsrf                            !
              dt_ver       = 0.                                !
! #kv       DO ikl=1,kcolp
! #kv       DO ikv=1,mwp
              isl          =          isnoSV(ikl,ikv)          !
              etaBAK(ikl,ikv)  = max(eps6,eta_SV(ikl,ikv ,isl))!
              etaNEW(ikl,ikv)  =          etaBAK(ikl,ikv)      !
              etEuBk(ikl,ikv)  =          etaNEW(ikl,ikv)      !
              RHu_av(ikl,ikv)  =     0.00                      !
              SHumSV(ikl,ikv)  =     0.00                      !
! #kv       END DO
! #kv       END DO
!         _______________
       DO it_srf=1,nt_srf                                      !
              dt_ver       = dt_ver     +dt_srf                !
! #kv       DO ikl=1,kcolp
! #kv       DO ikv=1,mwp
              faceta(ikl,ikv)  = fac_dt(ikl,ikv)*dt_srf        !
! #WL         faceta(ikl,ikv)  = faceta(ikl,ikv)              &!
! #WL&                  /(1.+faceta(ikl,ikv)*dQa_SV(ikl,ikv))  !    Limitation
!             CAUTION:     Please VERIFY dQa_SV Set-Up         ! by Atm.Conten
!    &        *max(0,sign(1.,qsatsg(ikl)-QaT_SV(ikl))))        ! NO Limitation
! #kv       END DO
! #kv       END DO                                             ! of Downw.Flux
          DO itEuBk=1,2                                        !
! #kv       DO ikl=1,kcolp
! #kv       DO ikv=1,mwp
              ist    = max(0,isotSV(ikl,ikv)-100*isnoSV(ikl,ikv))  ! 0 if    H2O
                                                               !
              Psi(ikl,ikv) =                                  &!
     &                psidSV(ist)                             &! DR97, Eqn 3.34
     &              *(etadSV(ist)                             &!
     &           /max(etEuBk(ikl,ikv),eps6))                  &!
     &              **bCHdSV(ist)                              !
              PsiArg(ikl,ikv) = 7.2E-5*Psi(ikl,ikv)            !
              RHuSol(ikl,ikv) =   exp(-min(ea_Max,PsiArg(ikl,ikv)))!
              SHuSol(ikl,ikv) =     qsatsg(ikl,ikv)  *RHuSol(ikl,ikv)  ! DR97, Eqn 3.15
              etEuBk(ikl,ikv) =                               &!
     &       (etaNEW(ikl,ikv) + faceta(ikl,ikv)*(QaT_SV(ikl,ikv)  &!
     &                                  -SHuSol(ikl,ikv)      &!
     &                    *(1.          -bCHdSV(ist)          &!
     &                                  *PsiArg(ikl,ikv))       ))&!
     &      /(1.          + faceta(ikl,ikv)* SHuSol(ikl,ikv)  &!
     &                                  *bCHdSV(ist)          &!
     &                                  *PsiArg(ikl,ikv)      &!
     &                                  /etaNEW(ikl,ikv))      !
              etEuBk(ikl,ikv) = etEuBk(ikl,ikv) -Rootsv(ikl,ikv,0)&!
     &                    * dt_srf     /(rhoWat*dz_dSV(0))     !
! #kv       END DO
! #kv       END DO
          END DO                                               !
! #kv       DO ikl=1,kcolp
! #kv       DO ikv=1,mwp
              RHu_av(ikl,ikv) =      RHu_av(ikl,ikv)+RHuSol(ikl,ikv)*dt_srf!
              SHumSV(ikl,ikv) =      SHumSV(ikl,ikv)+SHuSol(ikl,ikv)*dt_srf!
              etaNEW(ikl,ikv) =  max(etEuBk(ikl,ikv),eps6)     !
! #kv       END DO
! #kv       END DO
              dt_srf      =      xgpsrf                 *dt_srf!
       END DO                                                  !
 
 
! #kv   DO ikl=1,kcolp
! #kv   DO ikv=1,mwp
              RHu_av(ikl,ikv) =      RHu_av(ikl,ikv)            /dt_ver!
! tune        RHu_av(ikl) = 0.80*RHu_av(ikl)            /dt_ver!
 
! Surface  Type :  liquid or solid Water (i.e., H2O)  .OR.  SOIL
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          isl        =  isnoSV(ikl,ikv)                        !
          ist   = max(0,isotSV(ikl,ikv)-100*isnoSV(ikl,ikv))   ! 0 if    H2O
          ist__s= min(1,ist)                                   ! 1 if no H2O
          ist__w=     1-ist__s                                 ! 1 if    H2O
 
! Surface  Specific  Humidity :
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          SHumSV(ikl,ikv)=  ist__s         *SHumSV(ikl,ikv)  /dt_ver  &!      no H2O
     &                 +ist__w         *qsatsg(ikl,ikv)        !         H2O
 
! Surface: Latent    Heat Flux: Soil/Water Surface Contributions
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          d__eta     =  eta_SV(ikl,ikv,isl)-etaNEW(ikl,ikv)    !
          HL___D(ikl,ikv)=( ist__s *rhoWat *dz_dSV(0)         &! Soil Contrib.
     &                *(etaNEW(ikl,ikv)    -etaBAK(ikl,ikv)) / dt__SV &!
     &                 +ist__w         *f_HSHL(ikl,ikv)       &! H2O  Contrib.
     &                *(QaT_SV(ikl,ikv)    -qsatsg(ikl,ikv))         )&!
     &                * Lx_H2O(ikl,ikv)                        ! common factor
 
! Surface: Latent    Heat Flux: T Derivative
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          dLdTSV(ikl,ikv) = 0.
! #TL     dLdTSV(ikl,ikv) = f___HL(ikl,ikv)   * dqs_dT(ikl,ikv)   &!
! #TL&                                * RHu_av(ikl,ikv)        ! - d(HL)/d(T)
 
! #TL     HL___D(ikl,ikv) = HL___D(ikl,ikv)                   &!
! #TL&                 +dLdTSV(ikl,ikv)   * TsisSV(ikl,ikv,isl)!
! #kv   END DO
! #kv   END DO
 
! Surface: Tridiagonal Matrix Set Up
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #kv   DO ikl=1,kcolp
! #kv   DO ikv=1,mwp
!         ______________________________
          isl             =  isnoSV(ikl,ikv)
 
          TSurf0(ikl,ikv)     =  TsisSV(ikl,ikv,isl)             ! Diagnostic
          Elem_A          =  dtC_sv(ikl,ikv,isl)*mu__dz(ikl,ikv,isl) !
          Elem_C          =  dSdTSV(ikl,ikv)                     ! HS/Surf.Contr.
          ElemsC          =  Elem_C        * dtC_sv(ikl,ikv,isl) !
          ElemaC          =  Elem_C        *(dt__SV/CpdAir)     &!
     &                     /(roa_SV(ikl,ikv,1  )*dza__1(ikl,ikv))!
          Diag_A(ikl,ikv,isl) = -Elem_A *Implic                  !
          Diag_C(ikl,ikv,isl) = -ElemaC *Imp_SA                  !
          Diag_B(ikl,ikv,isl) =  1.0d+0 -Diag_A(ikl,ikv,isl)    &!
     &                     + ElemsC *ImpHSL                     &!
     &  + dtC_sv(ikl,ikv,isl) * (dIRsdT(ikl,ikv)                &! Upw. Sol IR
     &                      +dLdTSV(ikl,ikv))                    ! HL/Surf.Contr.
          Term_D(ikl,ikv,isl) =  Explic *Elem_A *Tsisva(ikl,ikv,isl-1)  &!
     &             +(1.0d+0 -Explic *Elem_A                     &!
     &                      -ExpHSL *ElemsC)*Tsisva(ikl,ikv,isl)&!
     &                      +Exp_SA *ElemaC *Tsisva(ikl,ikv,isl+1)   !                !T
          Term_D(ikl,ikv,isl) =  Term_D(ikl,ikv,isl)            &!                !T
     &  + dtC_sv(ikl,ikv,isl) * (sol_SV(ikl,ikv)    *SoSosv(ikl,ikv)&! Absorbed
     &                                     *(sEX_sv(ikl,ikv,isl+1)  &! Solar
     &                                      -sEX_sv(ikl,ikv,isl  )) &!
     &        +      tau_sv(ikl,ikv)      *IRd_SV(ikl,ikv) *Eso_sv(ikl,ikv) &! Down Atm IR
     &         -(1.0-tau_sv(ikl,ikv)) *0.5*IRv_sv(ikl,ikv)      &! Down Veg IR
     &                                +IRs__D(ikl,ikv)          &! Upw. Sol IR
     &                                +HL___D(ikl,ikv)          &! HL/Atmo.Contr.
     &      ) *      p0_kap       /    ExnrSV(ikl,ikv)           !
 
! SBL:        Tridiagonal Matrix Set Up
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!            vvvvvvvv
          IF (SVaKzT)                                               THEN
             Kz__SV(ikl,ikv,  1) = 0.
           IF(SVaUBC)                                               THEN
             Kz__SV(ikl,ikv,mzp) = 0.
             ExpTOP          = Explic
             ImpTOP          = Implic
           ELSE
             ExpTOP          = 1.0000
             ImpTOP          = 0.0000
           END IF
 
             isl=1
             jsl=isnoSV(ikl,ikv)
             isn=isnoSV(ikl,ikv) +isl
             is1=   min(isl+1,mzp)
             is2=   min(isn+1,mzp+nsnow)
             Elem_A          =  dSdTSV(ikl,ikv)                    &!
     &                         +dLdTSV(ikl,ikv)                    &!
     &                         +0.0d+0
             ElemaA          =  Elem_A        *(dt__SV/CpdAir)     &!
     &                        /(roa_SV(ikl,ikv,1  )*dza__1(ikl,ikv))!
             ElemsA          =  Elem_A        * dtC_sv(ikl,ikv,jsl)
             Diag_A(ikl,ikv,isn) = -ElemsA *Imp_SA                  !
             deltak=                            dt__SV             &!
     &                        /(roa_SV(ikl,ikv,1  )*dza__1(ikl,ikv))!
             Elem_C=deltak                                         &!
     &             *0.5*(roa_SV(ikl,ikv,is1  )+roa_SV(ikl,ikv,isl))&!
     &                 * Kz__SV(ikl,ikv,is1  )                     &!
     &       /max(epsi15,zza_SV(ikl,ikv,is1  )-zza_SV(ikl,ikv,isl)) !
             Diag_C(ikl,ikv,isn) = -Elem_C *Implic                  !
             Diag_B(ikl,ikv,isn) =  1.0d+0 -Diag_C(ikl,ikv,isn)    &!
     &                         +ElemaA *Implic                      !
             Term_D(ikl,ikv,isn) =  ElemsA *Exp_SA *Tsisva(ikl,ikv,isn-1)  &!
     &                +(1.0d+0 -ElemaA *Explic                     &!
     &                         -Elem_C *Explic)*Tsisva(ikl,ikv,isn)&!
     &                         +Elem_C *Explic *Tsisva(ikl,ikv,is2) !
 
! Atmosphere: Tridiagonal Matrix Set Up
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          IF(mzp.GE.4)                                              THEN
             is1=2
             is2=mzp-2
!            ___________
          DO isl=is1,is2
             isn=isnoSV(ikl,ikv)+isl                                !
             Elem_A=Elem_C                                          !
             deltak=0.5* dt__SV          /(roa_SV(ikl,ikv,isl)     &!
     &                 *(zza_SV(ikl,ikv,isl+1)-zza_SV(ikl,ikv,isl-1)))  !
             Elem_C=deltak                                         &!
     &             *0.5*(roa_SV(ikl,ikv,isl+1)+roa_SV(ikl,ikv,isl))&!
     &                 * Kz__SV(ikl,ikv,isl+1)                     &!
     &                 /(zza_SV(ikl,ikv,isl+1)-zza_SV(ikl,ikv,isl)) !
             Diag_A(ikl,ikv,isn) = -Implic *Elem_A                  !
             Diag_C(ikl,ikv,isn) = -Implic *Elem_C                  !
             Diag_B(ikl,ikv,isn) =  1.0d+0 -Diag_A(ikl,ikv,isn)    &!
     &                                 -Diag_C(ikl,ikv,isn)         !
             Term_D(ikl,ikv,isn) =  Explic *Elem_A *Tsisva(ikl,ikv,isn-1)  &!
     &                +(1.0d+0 -Explic *Elem_A                     &!
     &                         -Explic *Elem_C)*Tsisva(ikl,ikv,isn)&!
     &                         +Explic *Elem_C *Tsisva(ikl,ikv,isn+1)   !
          END DO
 
          ENDIF
 
! Atmosphere: Tridiagonal Matrix Set Up
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          IF(mzp .GE.2)                                             THEN
             is1=max(1,mzp-1)
!            _______
             isl=is1
             isn=isnoSV(ikl,ikv)+isl                                !
             Elem_A=Elem_C                                          !
             deltak=0.5* dt__SV          /(roa_SV(ikl,ikv,isl)     &!
     &                 *(zza_SV(ikl,ikv,isl+1)-zza_SV(ikl,ikv,isl-1)))  !
             Elem_C=deltak                                         &!
     &             *0.5*(roa_SV(ikl,ikv,isl+1)+roa_SV(ikl,ikv,isl))&!
     &                 * Kz__SV(ikl,ikv,isl+1)                     &!
     &                 /(zza_SV(ikl,ikv,isl+1)-zza_SV(ikl,ikv,isl)) !
             Diag_A(ikl,ikv,isn) = -Implic *Elem_A                  !
             Diag_C(ikl,ikv,isn) = -ImpTOP *Elem_C                  !
             Diag_B(ikl,ikv,isn) =  1.0d+0 -Diag_A(ikl,ikv,isn)    &!
     &                         +Implic *Elem_C                      !
             Term_D(ikl,ikv,isn) =  Explic *Elem_A *Tsisva(ikl,ikv,isn-1)  &!
     &                +(1.0d+0 -Explic *Elem_A                     &!
     &                         -Explic *Elem_C)*Tsisva(ikl,ikv,isn)&!
     &                         +ExpTOP *Elem_C *Tsisva(ikl,ikv,isn+1)   !
          END IF
 
 
             isHigh          =  isnoSV(ikl,ikv)+max(1,mzp-1)
          ELSE
             isHigh          =  isnoSV(ikl,ikv)
          END IF
! #kv   END DO
! #kv   END DO
 
 
 
! Tridiagonal Elimination
! =======================
 
! Forward  Sweep
! ^^^^^^^^^^^^^^
! #kv   DO ikl=  1,kcolp
! #kv   DO ikv=1,mwp
          Aux__P(ikl,ikv,-nsoil) = Diag_B(ikl,ikv,-nsoil)
          Aux__Q(ikl,ikv,-nsoil) =-Diag_C(ikl,ikv,-nsoil)/Aux__P(ikl,ikv,-nsoil)
! #kv   END DO
! #kv   END DO
 
!         ___________________
       DO isl=-nsoil+1,isHigh
! #kv   DO ikl=      1,kcolp
! #kv   DO ikv=1,mwp
          Aux__P(ikl,ikv,isl)   = Diag_A(ikl,ikv,isl)  *Aux__Q(ikl,ikv,isl-1)  &
     &                       +Diag_B(ikl,ikv,isl)
          Aux__Q(ikl,ikv,isl)   =-Diag_C(ikl,ikv,isl)  /Aux__P(ikl,ikv,isl)
! #kv   END DO
! #kv   END DO
       END DO
 
! #kv   DO ikl=      1,kcolp
! #kv   DO ikv=1,mwp
          Tsisva(ikl,ikv,-nsoil) = Term_D(ikl,ikv,-nsoil)/Aux__P(ikl,ikv,-nsoil)
! #kv   END DO
! #kv   END DO
 
!         ___________________
       DO isl=-nsoil+1,isHigh
! #kv   DO ikl=      1,kcolp
! #kv   DO ikv=1,mwp
          Tsisva(ikl,ikv,isl)   =(Term_D(ikl,ikv,isl)                  &
     &                       -Diag_A(ikl,ikv,isl)  *Tsisva(ikl,ikv,isl-1)) &
     &                                         /Aux__P(ikl,ikv,isl)
! #kv   END DO
! #kv   END DO
       END DO
 
! Backward Sweep
! ^^^^^^^^^^^^^^
!         ______________________
       DO isl=isHigh-1,-nsoil,-1
! #kv   DO ikl=1,kcolp
! #kv   DO ikv=1,mwp
          Tsisva(ikl,ikv,isl)   = Aux__Q(ikl,ikv,isl)  *Tsisva(ikl,ikv,isl+1)  &
     &                                         +Tsisva(ikl,ikv,isl)
! #kv   END DO
! #kv   END DO
       END DO
 
! Go Back to Temperatures
! ^^^^^^^^^^^^^^^^^^^^^^^
! #kv  DO ikl=1,kcolp
! #kv  DO ikv=1,mwp
         DO isl=-nsoil   , isnoSV(ikl,ikv)
          TsisSV(ikl,ikv,isl) =    Tsisva(ikl,ikv,isl)                 &
     &                /(p0_kap /   ExnrSV(ikl,ikv))
         END DO
 
!         vvvvvvvv
       IF (SVaKzT)                                                  THEN
!         _________
       DO isl=1,mzp
          isn=isnoSV(ikl,ikv)+isl
          pktaSV(ikl,ikv,isl) =   Tsisva(ikl,ikv,isn) / Fsisva(ikl,ikv)
       END DO
 
       IF  (SVaUBC)                                                 THEN
          is1             =      max(  1,mzp-1)
          pktaSV(ikl,ikv,mzp) =   pktaSV(ikl,ikv,is1)
!      ELSE
!         pktaSV(ikl,mzp)  is implicitely .NOT. modified
       END IF
 
       END IF
 
! Temperature Limits (avoids problems in case of no Snow Layers)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!         _______________________________
          isl               = isnoSV(ikl,ikv)
 
          dTSurf            = TsisSV(ikl,ikv,isl) -        TSurf0(ikl,ikv)
          TsisSV(ikl,ikv,isl)   = TSurf0(ikl,ikv) + sign(un_1,dTSurf)  &! 180.0 dgC/hr
     &              * min(abs(dTSurf),5.e-2*dt__SV)                 ! =0.05 dgC/s
! #kv  END DO
! #kv  END DO
 
!         ___________________
       DO isl=nsnow,1     ,-1
! #kv   DO ikl=     1,kcolp
! #kv   DO ikv=1,mwp
          TsisSV(ikl,ikv,isl)   = max(Ts_Min,       TsisSV(ikl,ikv,isl))
          TsisSV(ikl,ikv,isl)   = min(Ts_Max,       TsisSV(ikl,ikv,isl))
! #kv   END DO
! #kv   END DO
       END DO
 
 
! Update Surface    Fluxes
! ========================
 
! #kv   DO ikl=      1,kcolp
! #kv   DO ikv=1,mwp
!         _________________________
          isl         = isnoSV(ikl,ikv)
 
          IRs_SV(ikl,ikv) = IRs__D(ikl,ikv)                   &!
     &                - dIRsdT(ikl,ikv) * TsisSV(ikl,ikv,isl)  !
          HSs_sv(ikl,ikv) = HS___D(ikl,ikv)                   &! Sensible Heat
     &                - dSdTSV(ikl,ikv) * TsisSV(ikl,ikv,isl)  ! Downward > 0
          HLs_sv(ikl,ikv) = HL___D(ikl,ikv)                   &! Latent   Heat
     &                - dLdTSV(ikl,ikv) * TsisSV(ikl,ikv,isl)  ! Downward > 0
 
! OUTPUT for Stand Alone NetCDF File
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          SOsoKL(ikl,ikv) = sol_SV(ikl,ikv) * SoSosv(ikl,ikv)  ! Absorbed Sol.
          IRsoKL(ikl,ikv) =               IRs_SV(ikl,ikv)     &! Up Surf. IR
     &        +     tau_sv(ikl,ikv)      *IRd_SV(ikl,ikv)*Eso_sv(ikl,ikv) &! Down Atm IR
     &        -(1.0-tau_sv(ikl,ikv)) *0.5*IRv_sv(ikl,ikv)      ! Down Veg IR
          HSsoKL(ikl,ikv) = HSs_sv(ikl,ikv)                    ! HS
          HLsoKL(ikl,ikv) = HLs_sv(ikl,ikv)                    ! HL
          HLs_KL(ikl,ikv) = HLs_sv(ikl,ikv) / LhvH2O           ! mm w.e./sec
        END DO
        END DO
 
 
! OUTPUT/Verification: Energy/Water Budget: Energy Budget (OUT)
! #e1   DO ikl=1,kcolp
! #e1   DO ikv=1,mwp
! #e1     ETSo_d(ikl,ikv) =                                   &!
! #e1&         (     SoSosv(ikl,ikv)      *sol_SV(ikl,ikv)    &! Net   Solar
! #e1&         +                       IRs_SV(ikl,ikv)        &! Up Surf. IR
! #e1&         +     tau_sv(ikl,ikv)      *IRd_SV(ikl,ikv)*Eso_sv(ikl,ikv)&! Down Atm IR
! #e1&         -(1.0-tau_sv(ikl,ikv)) *0.5*IRv_sv(ikl,ikv)    &! Down Veg IR
! #e1&                                +HSs_sv(ikl,ikv)        &! Sensible
! #e1&                                +HLs_sv(ikl,ikv)            )! Latent
! #e1     ETSo_1(ikl,ikv) = 0.
! #e1   END DO
! #e1   END DO
!         _________________
! #e1  DO isl= -nsoil,nsnow
! #e1   DO ikl=1,kcolp
! #e1   DO ikv=1,mwp
! #e1     Exist0      = isl -           isnoSV(ikl,ikv)
! #e1     Exist0      = 1.  - max(zer0,min(un_1,Exist0))
! #e1     ETSo_1(ikl,ikv) = ETSo_1(ikl,ikv)                   &!
! #e1&                +(TsisSV(ikl,ikv,isl)-Tf_Sno)*Exist0    &!
! #e1&                                 /dtC_sv(ikl,ikv,isl)
! #e1   END DO
! #e1   END DO
! #e1  END DO
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
      IF (FlagDALLOC)                                             THEN !
 
      deallocate          ( Tsisva )                                   !
      deallocate          ( Fsisva )                                   !
      deallocate          ( dza__1 )                                   !
      deallocate          ( mu_sno )                                   !     Snow thermal Conductivity
      deallocate          ( mu__dz )                                   ! mu_(eta,sno)   / dz
      deallocate          ( dtC_sv )                                   ! dt      / C
      deallocate          ( IRs__D )                                   ! UpwardIR Previous Iter.Contr.
      deallocate          ( dIRsdT )                                   ! UpwardIR           T Derivat.
      deallocate          ( f_HSHL )                                   ! Factor common to HS and HL
      deallocate          ( dRidTs )                                   ! d(Rib)/d(Ts)
      deallocate          ( HS___D )                                   ! Sensible Heat Flux Atm.Contr.
      deallocate          ( f___HL )                                   !
      deallocate          ( HL___D )                                   ! Latent   Heat Flux Atm.Contr.
      deallocate          ( TSurf0 )                                   ! Previous Surface Temperature
      deallocate          ( qsatsg )                                   ! Soil   Saturat. Spec. Humidity
      deallocate          ( dqs_dT )                                   ! d(qsatsg)/dTv
      deallocate          ( Psi    )                                   ! 1st Soil Layer Water Potential
      deallocate          ( RHuSol )                                   ! Soil Surface Relative Humidity
      deallocate          ( RHu_av )                                   ! Soil Surface Relative Humidity
      deallocate          ( Diag_A )                                   ! A Diagonal
      deallocate          ( Diag_B )                                   ! B Diagonal
      deallocate          ( Diag_C )                                   ! C Diagonal
      deallocate          ( Term_D )                                   !   Independant Term
      deallocate          ( Aux__P )                                   ! P Auxiliary Variable
      deallocate          ( Aux__Q )                                   ! Q Auxiliary Variable
      deallocate          ( etaBAK )                                   !
      deallocate          ( etaNEW )                                   !
      deallocate          ( etEuBk )                                   !
      deallocate          ( fac_dt )                                   !
      deallocate          ( faceta )                                   !
      deallocate          ( PsiArg )                                   !
      deallocate          ( SHuSol )                                   !
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SISVAT_TSo
 
 
 
      subroutine SISVAT_qVg
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_qVg                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_qVg computes the Canopy Water  Balance             |
!                                    including  Root   Extraction          |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   ivgtSV   = 0,...,12:   Vegetation Type                      |
!     ^^^^^               0:          Water, Solid or Liquid               |
!                                                                          |
!     INPUT:   rhT_SV   : SBL Top    Air  Density                  [kg/m3] |
!     ^^^^^    QaT_SV   : SBL Top    Specific  Humidity            [kg/kg] |
!                                                                          |
!              TvegSV   : Canopy   Temperature                         [K] |
!              rrCaSV   : Canopy     Water     Content             [kg/m2] |
!              rrMxsv   : Canopy Maximum Intercepted Rain          [kg/m2] |
!              rah_sv   : Aerodynamic Resistance for Heat            [s/m] |
!              EvT_sv   : EvapoTranspiration                       [kg/m2] |
!              Sigmsv   : Canopy Ventilation Factor                    [-] |
!              glf_sv   : Green Leaf Fraction of NOT fallen Leaves     [-] |
!              LAIesv   : Leaf Area  Index (effective / transpiration) [-] |
!              psi_sv   : Soil       Water    Potential                [m] |
!              Khydsv   : Soil   Hydraulic    Conductivity           [m/s] |
!                                                                          |
!     INPUT /  psivSV   : Leaf       Water    Potential                [m] |
!     OUTPUT:                                                              |
!     ^^^^^^                                                               |
!                                                                          |
!     OUTPUT:  Rootsv   : Root Water Pump                        [kg/m2/s] |
!     ^^^^^^                                                               |
!                                                                          |
!     REMARK: Water Extraction by roots calibrated by Evapotranspiration   |
!     ^^^^^^  (computed in the Canopy Energy Balance)                      |
!                                                                          |
!     REFERENCE: DR97: Koen de Ridder thesis, UCL, 1997                    |
!     ^^^^^^^^^                                                            |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #RW: Root Water Flow slowed down by Soil Hydraulic Conductivity      |
!                                                                          |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
 
 
 
! Internal Variables
! ==================
 
      use Mod_SISVATLqVg
 
 
      IMPLICIT NONE
 
      integer                                      ::  ikl,ikv   ,isl  ! Grid Point, Layer Indices
      integer                                      ::  nitmax =  5     ! Maximum  Iterations Number
      integer                                      ::  nit             !          Iterations Counter
      real(kind=real8)                             ::  psidif          ! Soil-Canopy  Water Pot. Differ.
      real(kind=real8)                             ::  Root_W          ! Root   Water     Flow
      real(kind=real8)                             ::  RootOK          ! Roots take Water in Soil Layer
      real(kind=real8)                             ::  d_psiv          ! Canopy Water     Increment
      real(kind=real8)                             ::  dpvMAX = 20.    ! Canopy Water     Increment MAX
      real(kind=real8)                             ::  BWater          ! Imbalance of Canopy Water Budg.
      real(kind=real8)                             ::  BW_MAX          ! MAX Imbal.of Canopy Water Budg.
      real(kind=real8)                             ::  BW_MIN =  4.e-8 ! MIN Imbal.of Canopy Water Budg.
      real(kind=real8)                             ::  dBwdpv          ! Derivativ.of Canopy Water Budg.
      real(kind=real8)                             ::  Bswich          ! Newton-Raphson         Switch
      real(kind=real8)                             ::  EvFrac          ! Condensat./Transpiration Switch
      real(kind=real8)                             ::  den_qs,arg_qs   !
!     real(kind=real8)                             ::  esat_i          ! Saturation Vapor Pressure [hPa]
      real(kind=real8)                             ::  qsatvg          ! Canopy Saturat. Spec. Humidity
      real(kind=real8)                             ::  EvTran          ! EvapoTranspiration
      real(kind=real8)                             ::  dEdpsi          ! Evapotranspiration  Derivative
      real(kind=real8)                             ::  Fac_Ev,FacEvT   ! Evapotranspiration  Factors
      real(kind=real8)                             ::  denomE          ! Evapotranspiration  Denominator
      real(kind=real8)                             ::  F_Stom          ! Funct.  (Leaf Water Potential)
      real(kind=real8)                             ::  dFdpsi          ! Derivative  of F_Stom
      real(kind=real8)                             ::  denomF          ! Denominator of F_Stom
      real(kind=real8)                             ::  F___OK          ! (psi>psi_c) => F_Stom swich  ON
      real(kind=real8)                             ::  R0Stom          ! Minimum Stomatal Resistance
      real(kind=real8)                             ::  R_Stom          !         Stomatal Resistance
      real(kind=real8)                             ::  dRdpsi          ! Derivat.Stomatal Resistance
      real(kind=real8)                             ::  numerR          ! Numerat.Stomatal Resistance
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
      allocate            ( PlantW(kcolp,mwp) )                        ! Plant  Water
      allocate            ( dPdPsi(kcolp,mwp) )                        ! Plant  Water psi Derivative
      allocate            ( psiv_0(kcolp,mwp) )                        ! Canopy Temperature,  Previous t
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! Newton-Raphson Scheme
! =====================
 
      nit    = 0
  101 CONTINUE
      nit    = nit + 1
      BW_MAX = 0.
 
 
! W.Potential of the Previous Time Step
! -------------------------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          psiv_0(ikl,ikv) = psivSV(ikl,ikv)
 
 
! Extraction of Soil Water through the Plant Roots
! ------------------------------------------------
 
          PlantW(ikl,ikv)     = 0.                          ! Plant Water
          dPdPsi(ikl,ikv)     = 0.                          ! Idem, Derivat.
        END DO
        END DO
      DO   isl= -nsoil,0
        DO ikl=1,kcolp
        DO ikv=1,mwp
          psidif = psivSV(ikl,ikv)-(DH_dSV(ivgtSV(ikl,ikv))&! Soil-Canopy Water
     &                         +psi_sv(       ikl,ikv ,isl))! Potential  Diff.
          Root_W = rhoWat     * RF__SV(ivgtSV(ikl,ikv),isl)&! If > 0, Contrib.
     &              /max(eps_21,PR_dSV(ivgtSV(ikl,ikv))    &!     to Root Water
! #RW&                         +Khydsv(ikl,ikv   ,isl )*1.e-4  &! (DR97, eqn.3.20)
     &                                                   )  !
!         Pas de prise en compte de la resistance sol/racine dans proto-svat
!         (DR97, eqn.3.20)
          RootOK =      max(zer0, sign(un_1,psidif))
          Rootsv(ikl,ikv,isl) = Root_W*max(zer0,psidif)     ! Root  Water
          PlantW(ikl,ikv)     = PlantW(ikl,ikv) + Rootsv(ikl,ikv ,isl)  ! Plant Water
          dPdPsi(ikl,ikv)     = dPdPsi(ikl,ikv) + RootOK*Root_W ! idem, Derivat.
        END DO
        END DO
      END DO
 
 
! Latent   Heat Flux
! ------------------
 
! Canopy Saturation Specific Humidity
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        DO ikl=1,kcolp
        DO ikv=1,mwp
!       IF      (DeRidder)                                THEN !
          den_qs      =         TvegSV(ikl,ikv)    - 35.8      !
          arg_qs      = 17.27 *(TvegSV(ikl,ikv)    -273.16)   &!
     &                                   / den_qs              !
          qsatvg      = .0038 *        exp(arg_qs)     *0.875  !  0.875 = Tuning Hapex-Sahel
!         dqs_dT      = qsatvg     * 4099.2   /(den_qs *den_qs)!
!       ELSE IF (Dudhia_MAR)                              THEN !
!         esat_i      = 6.107                                 &!
!    &     *exp(ExpIsv*(un_1/WatIsv -un_1/TvegSV(ikl)    ))    !
!         qsatvg      = 0.622    *     esat_i                 &!
!    &      / (10.*pkPaSV(ikl) - 0.378*esat_i)                 !
!         dqs_dT      = qsatvg                                &!
!    &     *(1.0+0.6077*qsatvg     )                          &!
!    &     *    ExpIsv/(TvegSV(ikl)    *TvegSV(ikl)    )       !
!       END IF
 
 
! Canopy Stomatal Resistance
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
          R0Stom = min(         StodSV(ivgtSV(ikl,ikv))     &!
     &                /max(eps6,glf_sv(       ikl,ikv)),StxdSV)  ! Min Stomatal R.
          denomF =              pscdSV-psivSV(ikl,ikv)
          F___OK =     max(zer0,sign(un_1,denomF))
          denomF =     max(eps6,          denomF)            !
          F_Stom =              pscdSV  / denomF             ! F(Leaf Wat.Pot.)
          dFdpsi =             -F_Stom  / denomF             !
                                                             ! DR97, eqn. 3.22
          numerR = R0Stom / max(LAIesv(ikl,ikv), R0Stom/StxdSV)  !
          R_Stom = numerR *                  F_Stom          ! Can.Stomatal R.
                                                             ! DR97, eqn. 3.21
          dRdpsi = R_Stom *                  dFdpsi          !
 
! Evaporation / Evapotranspiration
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          EvFrac = max(zer0, sign(un_1,QaT_SV(ikl,ikv)-qsatvg))  ! Condensation/
          EvFrac = EvFrac                                   &! Transpiration
     &       + (1.-EvFrac)     *rrCaSV(ikl,ikv)/ rrMxsv(ikl,ikv) !        Switch
          Fac_Ev = rhT_SV(ikl,ikv) *Sigmsv(ikl,ikv)          ! idem,  Factor
          denomE = rah_sv(ikl,ikv) +R_Stom     * Sigmsv(ikl,ikv)
          FacEvT = Fac_Ev * (1.-EvFrac)    / denomE          !
          EvTran = FacEvT     *(qsatvg     - QaT_SV(ikl,ikv))! EvapoTranspir.
          dEdpsi =(EvTran     / denomE)    * dRdpsi          ! EvT Derivative
 
 
! Imbalance  of the Canopy  Water  Budget
! ---------------------------------------
 
          BWater      =( PlantW(ikl,ikv)                    &! Available  Water
     &                 - EvTran     )*        F___OK         ! Transpired Water
 
          Bswich      = max(zer0,                           &! Newton-Raphson
     &                      sign(un_1,    abs(BWater)       &!         Switch
     &                                       -BW_MIN))       !
 
 
! Derivative of the Canopy  Water  Budget
! ---------------------------------------
 
          dBwdpv      = dPdpsi(ikl,ikv)                     &!
     &                - dEdpsi
          dBwdpv      = sign(  un_1,    dBwdpv)             &!
     &                *  max(eps_21,abs(dBwdpv))             !
 
 
! Update Canopy and Surface/Canopy Temperatures
! ---------------------------------------------
 
          d_psiv      = BWater       / dBwdpv                !
          d_psiv      =      sign(un_1,d_psiv)              &! Increment
     &                       *min( abs(d_psiv)     ,dpvMAX)  ! Limitor
          psivSV(ikl,ikv) = psivSV(ikl,ikv)  - Bswich      *d_psiv   ! Newton-Raphson
          BW_MAX      = max(BW_MAX,abs(BWater))
        END DO
        END DO
 
 
! Update Root Water Fluxes | := Evapotranspiration
! ------------------------------------------------
 
      DO   isl= -nsoil,0
        DO ikl=1,kcolp
        DO ikv=1,mwp
          Rootsv(ikl,ikv,isl) = Rootsv(ikl,ikv,isl)*EvT_sv(ikl,ikv) &! Root  Water
     &                          /max(eps_21,PlantW(ikl,ikv)) !
        END DO
        END DO
      END DO
 
      IF (BW_MAX.gt.BW_MIN.and.nit.lt.nitmax)    GO TO 101
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
      IF (FlagDALLOC)                                             THEN !
 
      deallocate          ( PlantW )                                   ! Plant  Water
      deallocate          ( dPdPsi )                                   ! Plant  Water psi Derivative
      deallocate          ( psiv_0 )                                   ! Canopy Temperature,  Previous t
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SISVAT_qVg
 
 
 
      subroutine SISVAT_qSn                                            &
     &                     (                                           &
! #e1&                      EqSn_0,EqSn_1,EqSn_d                       &
! #m1&                     ,SIsubl,SImelt,SIrnof                       &
     &                     )
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_qSn                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_qSn updates  the Snow Water Content                |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   isnoSV   = total Nb of Ice/Snow Layers                      |
!     ^^^^^                                                                |
!                                                                          |
!     INPUT:   TaT_SV   : SBL Top    Temperature                       [K] |
!     ^^^^^    dt__SV   : Time Step                                    [s] |
!                                                                          |
!     INPUT /  drr_SV   : Rain Intensity                         [kg/m2/s] |
!     OUTPUT:  dzsnSV   : Snow Layer Thickness                         [m] |
!     ^^^^^^   eta_SV   : Snow Water Content                       [m3/m3] |
!              ro__SV   : Snow/Soil  Volumic Mass                  [kg/m3] |
!              TsisSV   : Soil/Ice Temperatures (layers -nsoil,-nsoil+1, 0)|
!                       & Snow     Temperatures (layers  1,2,..,nsnow) [K] |
!                                                                          |
!     OUTPUT:  SWS_SV   : Surficial Water Status                           |
!     ^^^^^^                                                               |
!              EExcsv   : Snow Energy in Excess, initial Forcing    [J/m2] |
!              EqSn_d   : Snow Energy in Excess, remaining          [J/m2] |
!              EqSn_0   : Snow Energy, before Phase Change          [J/m2] |
!              EqSn_1   : Snow Energy, after  Phase Change          [J/m2] |
!              SIsubl   : Snow sublimed/deposed Mass             [mm w.e.] |
!              SImelt   : Snow Melted           Mass             [mm w.e.] |
!              SIrnof   : Surficial Water + Run OFF Change       [mm w.e.] |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: STANDARD Possibility                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                          |
!     #IB: OUTPUT: Ice-Sheet Surface Mass Balance  (on NetCDF File )       |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:           (PLEASE VERIFY before USE)          |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #SU: SLUSH : Alternative Parameterization                            |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_iii_jjj_n     | #e0: OUTPUT on ASCII  File (SISVAT Variables) |
!   #                      |(#e0  MUST BE PREPROCESSED BEFORE #e1 & #e2 !) |
!   # SISVAT_iii_jjj_n     | #e1: OUTPUT/Verification: Energy Conservation |
!   # SISVAT_iii_jjj_n     | #m1: OUTPUT/Verification: * Mass Conservation |
!                          |                                               |
!   # SISVAT_qSn.vm        | #e5: OUTPUT/Verification: Energy/Water Budget |
!                          |      unit 43, SubRoutine  SISVAT_qSn **ONLY** |
!   # SISVAT_qSn.vu        | #vu: OUTPUT/Verification: Slush  Parameteriz. |
!                          |      unit 44, SubRoutine  SISVAT_qSn **ONLY** |
!                          |                                               |
!   # stdout               | #s2: OUTPUT of SnowFall, Snow Buffer          |
!                          |      unit  6, SubRoutine  SISVAT_BSn, _qSn    |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_ctr
      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
      use Mod_SISVAT_qSn
 
 
 
! Internal Variables
! ==================
 
      use Mod_SISVATLqSn
 
 
      IMPLICIT NONE
 
 
      integer                                      ::  ikl,ikv   ,isn  !
      integer                                      ::  nh              ! Non erodible Snow: up.lay.Index
      integer                                      ::  LayrOK          ! 1 (0)  if In(Above) Snow Pack
      integer                                      ::  k_face          ! 1 (0)  if Crystal(no) faceted
      integer                                      ::  LastOK          ! 1 ==>  1! Snow Layer
      integer                                      ::  NOLayr          ! 1     Layer  Update
! #SU integer                                      ::  kSlush          ! Slush Switch
      real(kind=real8)                             ::  dTSnow          ! Temperature                  [C]
      real(kind=real8)                             ::  OKmelt          ! 1 (0)  if        (no) Melting
      real(kind=real8)                             ::  EnMelt          ! Energy in excess, for Melting
      real(kind=real8)                             ::  SnHLat          ! Energy consumed   in  Melting
! #e4 real(kind=real8)                             ::  AdEnrg          ! Additional Energy from  Vapor
! #e3 real(kind=real8)                             ::  B_Enrg          ! Additional Energy from  Vapor
      real(kind=real8)                             ::  dzVap0,dzVap1   ! Vaporized Thickness          [m]
      real(kind=real8)                             ::  rosDry          ! Snow volumic Mass if no Water in
      real(kind=real8)                             ::  PorVol          ! Pore volume
      real(kind=real8)                             ::  PClose          ! Pore Hole Close OFF Switch
!     real(kind=real8)                             ::  SGDiam          !      Snow Grain Diameter
!     real(kind=real8)                             ::  SGDmax = 0.003  ! Max. Snow Grain Diameter     [m]
                                                                       ! (Rowe et al. 1995, JGR p.16268)
      real(kind=real8)                             ::  rWater          ! Retained Water           [kg/m2]
      real(kind=real8)                             ::  drrNEW          ! New available Water      [kg/m2]
      real(kind=real8)                             ::  rdzNEW          ! Snow          Mass       [kg/m2]
      real(kind=real8)                             ::  rdzsno          ! Snow          Mass       [kg/m2]
      real(kind=real8)                             ::  EnFrez          ! Energy Release    in  Freezing
      real(kind=real8)                             ::  WaFrez          ! Water  consumed   in  Melting
      real(kind=real8)                             ::  RapdOK          ! 1. ==> Snow melts rapidly
      real(kind=real8)                             ::  ThinOK          ! 1. ==> Snow Layer is thin
      real(kind=real8)                             ::  dzepsi = 0.0001 ! Minim. Snow Layer Thickness (!)
      real(kind=real8)                             ::  dz_Min = 1.e-3  ! Minim. Snow Layer Thickness
!                                                      dz_Min = 0.005  !
! ... Warning: Too high for Col de Porte: precludes 1st snow (layer) apparition
 
      real(kind=real8)                             ::  z_Melt          ! Last (thin) Layer Melting
      real(kind=real8)                             ::  rusnew          ! Surficial Water Thickness   [mm]
! #SU real(kind=real8)                             ::  zWater          ! Max Slush Water Thickness   [mm]
! #SU real(kind=real8)                             ::  zSlush          !     Slush Water Thickness   [mm]
! #SU real(kind=real8)                             ::  ro_new          ! New Snow/ice  Density    [kg/m3]
      real(kind=real8)                             ::  zc,zt           ! Non erod.Snow Thickness[mm w.e.]
 
 
! OUTPUT of SISVAT Trace Statistics (see assignation in PHY_SISVAT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer                                      ::  isnnew,isinew,isnUpD,isnitr
 
 
! #e5 real(kind=real8)                             ::  hourer          !
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
! OUTPUT/Verification: Energy/Water Budget
! #e1 allocate            ( EqSn_d(kcolp,mwp) )                        ! Energy in Excess, initial
! #e1 allocate            ( EqSn_0(kcolp,mwp) )                        ! Snow Energy, befor Phase Change
! #e5 allocate            ( EqSn01(kcolp,mwp) )                        ! Snow Energy, after Phase Change
! #e5 allocate            ( EqSn02(kcolp,mwp) )                        ! Snow Energy, after Phase Change
                                                                       !              .AND. Last Melting
! #e1 allocate            ( EqSn_1(kcolp,mwp) )                        ! Snow Energy, after Phase Change
                                                                       !              .AND. Mass Redistr.
! OUTPUT/Verification: * Mass Conservation
! #m1 allocate            ( SIsubl(kcolp,mwp) )                        ! Snow Deposed Mass
! #m1 allocate            ( SImelt(kcolp,mwp) )                        ! Snow Melted  Mass
! #m1 allocate            ( SIrnof(kcolp,mwp) )                        ! Local Surficial Water + Run OFF
 
      allocate            ( noSnow(kcolp,mwp) )                        ! Nb of Layers Updater
      allocate            ( EExdum(kcolp,mwp) )                        ! Energy in Excess when no Snow
      allocate            ( dzMelt(kcolp,mwp) )                        ! Melted    Thickness          [m]
 
! OUTPUT/Verification: Energy/Water Budget
! #e5 allocate            ( WqSn_0(kcolp,mwp) )                        ! Snow Water+Forcing  Initial
! #e5 allocate            ( WqSn_1(kcolp,mwp) )                        ! Snow Water+Forcing, Final
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! OUTPUT/Verification: Energy/Water Budget: Energy Budget (IN)
! #e1   DO ikl=1,kcolp
! #e1   DO ikv=1,mwp
! #e1     EqSn_0(ikl,ikv) = 0.
! #e1   END DO
! #e1   END DO
! #e1 DO   isn=nsnow,1,-1
! #e1   DO ikl=1,kcolp
! #e1   DO ikv=1,mwp
! #e1     EqSn_0(ikl,ikv) = EqSn_0(ikl,ikv) + ro__SV(ikl,ikv,isn) *dzsnSV(ikl,ikv,isn) &
! #e1&                *(Cn_dSV      *(TsisSV(ikl,ikv,isn) -Tf_Sno         )&
! #e1&                 -LhfH2O      *(1.              -eta_SV(ikl,ikv,isn)))
! #e1   END DO
! #e1   END DO
! #e1 END DO
 
 
! OUTPUT/Verification: Energy/Water Budget: Water  Budget (IN)
! #e5   DO ikl=1,kcolp
! #e5   DO ikv=1,mwp
! #e5     WqSn_0(ikl,ikv) = drr_SV(ikl,ikv) * dt__SV                   &
! #e5&                 +rusnSV(ikl,ikv)
! #e5   END DO
! #e5   END DO
! #e5 DO   isn=nsnow,1,-1
! #e5   DO ikl=1,kcolp
! #e5   DO ikv=1,mwp
! #e5     WqSn_0(ikl,ikv) = WqSn_0(ikl,ikv) + ro__SV(ikl,ikv,isn) *dzsnSV(ikl,ikv,isn)
! #e5   END DO
! #e5   END DO
! #e5 END DO
 
 
! OUTPUT/Verification: * Mass Conservation
! #m1   DO ikl=1,kcolp
! #m1   DO ikv=1,mwp
! #m1     SImelt(ikl,ikv) = 0.
! #m1     SIrnof(ikl,ikv) = rusnSV(ikl,ikv) + RnofSV(ikl,ikv) * dt__SV
! #m1   END DO
! #m1   END DO
 
 
! Initialization
! ==============
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
        noSnow(ikl,ikv)   = 0               ! Nb of Layers Updater
        ispiSV(ikl,ikv)   = 0               ! Pore Hole Close OFF Index
                                            ! (assumed to be the Top of
                                            !  the surimposed Ice Layer)
! #IB   dwemSV(ikl,ikv)   = 0.
! #IB   dwerSV(ikl,ikv)   = 0.
      END DO
      END DO
 
 
! Melting/Freezing Energy
! =======================
 
!  REMARK: Snow liquid Water Temperature assumed = Tf_Sno
!  ^^^^^^
        DO ikl=1,kcolp
        DO ikv=1,mwp
          EExdum(ikl,ikv) = drr_SV(ikl,ikv)     * hC_Wat *(TaT_SV(ikl,ikv)-Tf_Sno) &
     &                                  * dt__SV
          EExcsv(ikl,ikv) = EExdum(ikl,ikv)     *    min(1,isnoSV(ikl,ikv)) ! Snow exists
          EExdum(ikl,ikv) = EExdum(ikl,ikv)     -          EExcsv(ikl,ikv)  !
 
! OUTPUT/Verification: Energy/Water Budget
! #e1     EqSn_d(ikl,ikv) = EExcsv(ikl,ikv)                     !
 
        END DO
        END DO
 
 
! Surficial Water Status
! ----------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          SWS_SV(ikl,ikv) = max(zer0,sign(un_1,Tf_Sno                  &
     &                                    -TsisSV(ikl,ikv,isnoSV(ikl,ikv))))
        END DO
        END DO
 
      DO   isn=nsnow,1,-1
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
! Energy, store Previous Content
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          dTSnow      = TsisSV(ikl,ikv,isn) -          Tf_Sno
          EExcsv(ikl,ikv) = EExcsv(ikl,ikv)                            &
     &                + ro__SV(ikl,ikv,isn) * Cn_dSV * dTSnow          &
     &                                           * dzsnSV(ikl,ikv,isn)
          TsisSV(ikl,ikv,isn) =                        Tf_Sno
 
! Water,  store Previous Content
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          drr_SV(ikl,ikv) = drr_SV(ikl,ikv)                            &
     &                + ro__SV(ikl,ikv,isn)          * eta_SV(ikl,ikv,isn) &
     &                                           * dzsnSV(ikl,ikv,isn) &
     &                / dt__SV
          ro__SV(ikl,ikv,isn) =                                        &
     &                  ro__SV(ikl,ikv,isn) *(1.     - eta_SV(ikl,ikv,isn))
          eta_SV(ikl,ikv,isn) =  0.
 
 
! Melting  if EExcsv > 0
! ======================
 
          EnMelt      =    max(zer0,          EExcsv(ikl,ikv) )
 
! Energy Consumption
! ^^^^^^^^^^^^^^^^^^
          SnHLat      = ro__SV(ikl,ikv,isn) * LhfH2O          !
          dzMelt(ikl,ikv) = EnMelt      / max(SnHLat,    eps6 )   !
          noSnow(ikl,ikv) = noSnow(ikl,ikv)                  &!
     &   +int(max(zer0  ,sign(un_1,dzMelt(ikl,ikv)           &!
     &                            -dzsnSV(ikl,ikv ,isn))))   &! 1 if full Melt
     &       *min(1     , max(0 ,1+isnoSV(ikl,ikv)-isn))      ! 1 in the  Pack
          dzMelt(ikl,ikv) =                                  &!
     &              min(dzsnSV(ikl,ikv, isn),dzMelt(ikl,ikv)) !
          dzsnSV(ikl,ikv,isn) =                              &!
     &                  dzsnSV(ikl,ikv,isn) -dzMelt(ikl,ikv)  !
          EExcsv(ikl,ikv) = EExcsv(ikl,ikv)     -dzMelt(ikl,ikv)*SnHLat   !
! #IB     dwemSV(ikl,ikv) = dwemSV(ikl,ikv)     -dzMelt(ikl,ikv)*ro__SV(ikl,ikv,isn)
 
! Water  Production
! ^^^^^^^^^^^^^^^^^
          drr_SV(ikl,ikv) = drr_SV(ikl,ikv)                            &
     &                + ro__SV(ikl,ikv,isn) * dzMelt(ikl,ikv)/dt__SV
 
! OUTPUT/Verification: * Mass Conservation
! #m1     SImelt(ikl,ikv) = SImelt(ikl,ikv)                            &
! #m1&                + ro__SV(ikl,ikv,isn) * dzMelt(ikl,ikv)
 
          OKmelt      =max(zer0,sign(un_1,drr_SV(ikl,ikv)-eps6))
 
! Snow History
! ^^^^^^^^^^^^
          k_face          =       min(    istoSV(ikl,ikv,isn),istdSV(1)) &! = 1  if
     &                           *max(0,2-istoSV(ikl,ikv,isn)          )  ! faceted
          istoSV(ikl,ikv,isn) =                                      &!
     &     int(1.-OKmelt) *               istoSV(ikl,ikv,isn)        &!
     &    +int   (OKmelt) *((1-k_face) *  istdSV(2)                  &!
     &                     +   k_face  *  istdSV(3)      )            !
 
 
! Freezing if EExcsv < 0
! ======================
 
          rdzsno      =          ro__SV(ikl,ikv,isn) * dzsnSV(ikl,ikv ,isn)
          LayrOK      = min(   1, max(0          , isnoSV(ikl,ikv)-isn+1))
          EnFrez      = min(zer0,                  EExcsv(ikl,ikv))
          WaFrez      =   -(     EnFrez          * LayrOK / LhfH2O)
          drrNEW      = max(zer0,drr_SV(ikl,ikv)     - WaFrez / dt__SV)
          WaFrez      =    (     drr_SV(ikl,ikv)     - drrNEW)* dt__SV
          drr_SV(ikl,ikv) =          drrNEW
          EExcsv(ikl,ikv) =          EExcsv(ikl,ikv)     + WaFrez * LhfH2O
          EnFrez      = min(zer0,EExcsv(ikl,ikv))    * LayrOK
          rdzNEW      = WaFrez + rdzsno
          ro__SV(ikl,ikv,isn) =      rdzNEW /max(eps6, dzsnSV(ikl,ikv,isn))
          TsisSV(ikl,ikv,isn) =      Tf_Sno                            &
     &                + EnFrez /(Cn_dSV *max(eps6, rdzNEW)        )
          EExcsv(ikl,ikv) =          EExcsv(ikl,ikv)     - EnFrez
! #IB     dwerSV(ikl,ikv) = WaFrez                                     &
! #IB&                + dwerSV(ikl,ikv)
 
 
! Snow Water Content
! ==================
 
! Pore   Volume [-]
! ^^^^^^^^^^^^^^^^^
          rosDry      =(1.     - eta_SV(ikl,ikv,isn))* ro__SV(ikl,ikv,isn) !
          PorVol      = 1.     - rosDry          / rhoIce          !
          PorVol      =      max(PorVol          , zer0  )         !
 
! Water  Retention
! ^^^^^^^^^^^^^^^^
          rWater      = ws0dSV * PorVol * rhoWat * dzsnSV(ikl,ikv,isn)
          drrNEW      = max(zer0,drr_SV(ikl,ikv)     - rWater /dt__SV)
          rWater      =    (     drr_SV(ikl,ikv)     - drrNEW)*dt__SV
          drr_SV(ikl,ikv)     =      drrNEW
          rdzNEW          =      rWater                                &
     &                         + rosDry          * dzsnSV(ikl,ikv,isn)
          eta_SV(ikl,ikv,isn) =      rWater / max(eps6,rdzNEW)
          ro__SV(ikl,ikv,isn) =      rdzNEW / max(eps6,dzsnSV(ikl,ikv,isn))
 
! Pore Hole Close OFF
! ^^^^^^^^^^^^^^^^^^^
          PClose = max(zer0,                                           &
     &                 sign(un_1,ro__SV(ikl,ikv,isn)                   &
     &                          -roCdSV         ))
          ispiSV(ikl,ikv) =          ispiSV(ikl,ikv)  * int(1.-PClose) &
     &                +      max(ispiSV(ikl,ikv),isn) *int(Pclose)
          PClose = max(0   ,                                           &! Water under SuPer.Ice
     &                 min (1   ,ispiSV(ikl,ikv)                       &! contributes to
     &                          -isn            ))                      ! Surficial   Water
          rusnSV(ikl,ikv) =          rusnSV(ikl,ikv)                   &
     &                +          drr_SV(ikl,ikv) *dt__SV * PClose
          drr_SV(ikl,ikv) =          drr_SV(ikl,ikv)      *(1.-PClose)
 
        END DO
        END DO
      END DO
 
 
! Remove Zero-Thickness Layers
! ============================
 
 1000 CONTINUE
           isnitr =          0
      DO   ikl=1,kcolp
      DO ikv=1,mwp
           isnUpD =          0
           isinew =          0
        DO isn=1,nsnow-1
           isnnew =                                                    &
     &          int(un_1-max(zer0  ,sign(un_1,dzsnSV(ikl,ikv,isn)-dzepsi)))&
     &             *     max(0     , min(1   ,isnoSV(ikl,ikv) +1 -isn ))
           isnUpD =      max(isnUpD,          isnnew)
           isnitr =      max(isnitr,          isnnew)
           isinew =      isn*isnUpD *max(0, 1-isinew)                  &! LowerMost  0-Layer
     &                                       +isinew                    ! Index
           dzsnSV(ikl,ikv,isn) =                  dzsnSV(ikl,ikv,isn+isnnew)
           ro__SV(ikl,ikv,isn) =                  ro__SV(ikl,ikv,isn+isnnew)
           TsisSV(ikl,ikv,isn) =                  TsisSV(ikl,ikv,isn+isnnew)
           eta_SV(ikl,ikv,isn) =                  eta_SV(ikl,ikv,isn+isnnew)
           G1snSV(ikl,ikv,isn) =                  G1snSV(ikl,ikv,isn+isnnew)
           G2snSV(ikl,ikv,isn) =                  G2snSV(ikl,ikv,isn+isnnew)
           dzsnSV(ikl,ikv,isn+isnnew) =(1-isnnew)*dzsnSV(ikl,ikv,isn+isnnew)
           ro__SV(ikl,ikv,isn+isnnew) =(1-isnnew)*ro__SV(ikl,ikv,isn+isnnew)
           eta_SV(ikl,ikv,isn+isnnew) =(1-isnnew)*eta_SV(ikl,ikv,isn+isnnew)
           G1snSV(ikl,ikv,isn+isnnew) =(1-isnnew)*G1snSV(ikl,ikv,isn+isnnew)
           G2snSV(ikl,ikv,isn+isnnew) =(1-isnnew)*G2snSV(ikl,ikv,isn+isnnew)
        END DO
           isnoSV(ikl,ikv)   =   isnoSV(ikl,ikv)-isnUpD                 ! Nb of Snow   Layer
           ispiSV(ikl,ikv)   =   ispiSV(ikl,ikv)                       &! Nb of SuperI Layer
     &    -isnUpD *max(0,min(ispiSV(ikl,ikv)-isinew,1))                 ! Update  if I=0
 
      END DO
      END DO
      IF  (isnitr.GT.0)                                       GO TO 1000
 
 
! New upper Limit of the non erodible Snow (istoSV .GT. 1)
! ========================================
 
      DO   ikl=1,kcolp
      DO ikv=1,mwp
           nh =     0
        DO isn=  nsnow,1,-1
           nh =    nh + isn* min(istoSV(ikl,ikv,isn)-1,1)*max(0,1-nh)
        END DO
           zc =     0.
           zt =     0.
        DO isn=1,nsnow
           zc =    zc +          dzsnSV(ikl,ikv,isn) *ro__SV(ikl,ikv,isn)  &
     &                     * max(0,min(1,nh+1-isn))
           zt =    zt +          dzsnSV(ikl,ikv,isn) *ro__SV(ikl,ikv,isn)
        END DO
           zWE_SV(ikl,ikv) =                 zt
           zWEcSV(ikl,ikv) = min(zWEcSV(ikl,ikv),zt)
           zWEcSV(ikl,ikv) = max(zWEcSV(ikl,ikv),zc)
      END DO
      END DO
 
 
! OUTPUT/Verification: Energy/Water Budget: Energy Budget (OUT)
! #e5   DO ikl=1,kcolp
! #e5   DO ikv=1,mwp
! #e5     EqSn01(ikl,ikv) =-EqSn_0(ikl,ikv)                            &
! #e5&                 -EExcsv(ikl,ikv)
! #e5   END DO
! #e5   END DO
! #e5 DO   isn=nsnow,1,-1
! #e5   DO ikl=1,kcolp
! #e5   DO ikv=1,mwp
! #e5     EqSn01(ikl,ikv) = EqSn01(ikl,ikv) + ro__SV(ikl,ikv,isn) *dzsnSV(ikl,ikv,isn) &
! #e5&                *(Cn_dSV      *(TsisSV(ikl,ikv,isn) -Tf_Sno         )&
! #e5&                 -LhfH2O      *(1.              -eta_SV(ikl,ikv,isn)))
! #e5   END DO
! #e5   END DO
! #e5 END DO
 
 
! "Negative Heat" from supercooled rain
!  ------------------------------------
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
          EExcsv(ikl,ikv) = EExcsv(ikl,ikv) + EExdum(ikl,ikv)
 
 
! Surficial Water Run OFF
! -----------------------
 
          rusnew      = rusnSV(ikl,ikv) * SWf_SV(ikl,ikv)
          RnofSV(ikl,ikv) = RnofSV(ikl,ikv)                            &
     &                +(rusnSV(ikl,ikv) - rusnew     ) / dt__SV
          rusnSV(ikl,ikv) = rusnew
      END DO
      END DO
 
 
! Percolation down the Continental Ice Pack
! -----------------------------------------
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          drr_SV(ikl,ikv) = drr_SV(ikl,ikv) + rusnSV(ikl,ikv)          &
     &                     * (1-min(1,ispiSV(ikl,ikv)))/ dt__SV
          rusnSV(ikl,ikv) = rusnSV(ikl,ikv)                            &
     &                     *    min(1,ispiSV(ikl,ikv))
        END DO
        END DO
 
 
! Slush Formation (CAUTION: ADD RunOff Possibility before Activation)
! ---------------  ^^^^^^^  ^^^
 
! OUTPUT/Verification: Slush  Parameterization
! #vu IF (.NOT.su_opn)                                              THEN
! #vu          su_opn=.true.
! #vu     open(unit=44,status='unknown',file='SISVAT_qSn.vu')
! #vu     rewind    44
! #vu END IF
! #vu     write(44,440) daHost
! #vu 440 format('iSupI    i       dz       ro      eta',              &
! #vu&            '   PorVol   zSlush     ro_n    eta_n',2x,a18)
 
! #SU DO   isn=1,nsnow
! #SU   DO ikl=1,kcolp
! #SU   DO ikv=1,mwp
! #SU     kSlush = min(1,max(0,isn+1-ispiSV(ikl,ikv)))    ! Slush Switch
 
! Available Additional Pore   Volume [-]
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #SU     PorVol = 1. - ro__SV(ikl,ikv,isn)                &! [--]
! #SU&           *(1. - eta_SV(ikl,ikv,isn))/ rhoIce       &!
! #SU&           -      eta_SV(ikl,ikv,isn)                &!
! #SU&                 *ro__SV(ikl,ikv,isn) / rhoWat        !
! #SU     PorVol =  max(PorVol          , zer0  )           !
! #SU     zWater =      dzsnSV(ikl,ikv,isn) * PorVol * 1000.   &! [mm] OR [kg/m2]
! #SU&           * (1. -SWS_SV(ikl,ikv)                    &! 0 <=> freezing
! #SU&                *(1 -min(1,iabs(isn-isnoSV(ikl,ikv)))))   ! 1 <=> isn=isnoSV
! #SU     zSlush =  min(rusnSV(ikl,ikv)     , zWater)       ! [mm] OR [kg/m2]
! #SU     rusnSV(ikl,ikv) = rusnSV(ikl,ikv)     - zSlush    ! [mm] OR [kg/m2]
! #SU     ro_new      =(dzsnSV(ikl,ikv,isn) * ro__SV(ikl,ikv,isn)  &!
! #SU&                 +zSlush                           ) &!
! #SU&            / max(dzsnSV(ikl,ikv,isn) , eps6           )  !
 
! OUTPUT/Verification: Slush  Parameterization
! #vu     rusnew          = eta_SV(ikl,ikv,isn)             !
 
! #SU     eta_SV(ikl,ikv,isn) =(ro_new - ro__SV(ikl,ikv,isn)   &!
! #SU&                    *(1.     - eta_SV(ikl,ikv,isn))) &!
! #SU&               / max (ro_new , eps6            )      !
 
! OUTPUT/Verification: Slush  Parameterization
! #vu     IF    (isn.le.isnoSV(ikl,ikv))                   &!
! #vu&    write(44,441) ispiSV(ikl,ikv),isn,dzsnSV(ikl,ikv,isn)&!
! #vu&                 ,ro__SV(ikl,ikv,isn),rusnew         &!
! #vu&                 ,PorVol         ,zSlush             &!
! #vu&                 ,ro_new         ,eta_SV(ikl,ikv,isn) !
! #vu 441 format(2i5,f9.3,f9.1,f9.6,f9.3,f9.6,f9.1,f9.6)    !
 
! #SU     ro__SV(ikl,ikv,isn) =      ro_new                 !
! #SU   END DO
! #SU   END DO
! #SU END DO
 
 
! Impact of the Sublimation/Deposition on the Surface Mass Balance
! ================================================================
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
          isn                     = isnoSV(ikl,ikv)
          dzVap0                  =                       dt__SV       &
     &  * HLs_sv(ikl,ikv)             * min(isn             , 1   )    &
     &  /(Lx_H2O(ikl,ikv)             * max(ro__SV(ikl,ikv,isn) , eps6))
          NOLayr=int(min(zer0,sign(un_1,dzsnSV(ikl,ikv,isn) + dzVap0)))
          dzVap1=    min(zer0,          dzsnSV(ikl,ikv,isn) + dzVap0)
 
 
! Additional Energy (CAUTION: Verification is not performed)
! -----------------
 
! OUTPUT/Verification: Energy Consrv. (HLS)
! #e4     AdEnrg = dzVap0 * ro__SV(ikl,ikv,isnoSV(ikl,ikv))            &! Water   Vapor
! #e4&            *hC_Wat *(TsisSV(ikl,ikv,isnoSV(ikl,ikv)) -Tf_Sno)    ! Sensible Heat
 
! OUTPUT/Verification: Energy Consrv. (HL)
! #e3     B_Enrg =(Cn_dSV      *(TsisSV(ikl,ikv,isn) -Tf_Sno         ) &
! #e3&            -LhfH2O      *(1.              -eta_SV(ikl,ikv,isn)))&
! #e3&           /(1.          + dzVap0 /max(eps6,dzsnSV(ikl,ikv,isn)))
! #e3     eta_SV(ikl,ikv,isn) =                                        &
! #e3&           max(zer0,un_1 +(B_Enrg                                &
! #e3&                         -(TsisSV(ikl,ikv,isn) -Tf_Sno)*Cn_dSV)  &
! #e3&                          /LhfH2O                          )
! #e3     TsisSV(ikl,ikv,isn) =    ( B_Enrg                            &
! #e3&                         +(1.              -eta_SV(ikl,ikv,isn)) &
! #e3&                          *LhfH2O                          )     &
! #e3&                         / Cn_dSV                                &
! #e3&                         + Tf_Sno
 
! OUTPUT/Verification: Energy Conservation
! #e1     STOP "PLEASE add Energy (#e3) from deposition/sublimation"
 
 
! Update of the upper Snow layer Thickness
! ----------------------------------------
 
          dzsnSV(ikl,ikv,isn) =                                        &
     &           max(zer0,  dzsnSV(ikl,ikv,isnoSV(ikl,ikv)) + dzVap0)
          isnoSV(ikl,ikv)     = isnoSV(ikl,ikv)             + NOLayr
          isn             = isnoSV(ikl,ikv)
          dzsnSV(ikl,ikv,isn) = dzsnSV(ikl,ikv,isn) + dzVap1
! #IB     dwesSV(ikl,ikv)     = ro__SV(ikl,ikv,isn) * dzVap0
      END DO
      END DO
 
 
! OUTPUT/Verification: Energy/Water Budget: Energy Budget (OUT)
! #e5   DO ikl=1,kcolp
! #e5   DO ikv=1,mwp
! #e5     EqSn02(ikl,ikv) =-EqSn_0(ikl,ikv)                            &
! #e5&                 -EExcsv(ikl,ikv)
! #e5   END DO
! #e5   END DO
! #e5 DO   isn=nsnow,1,-1
! #e5   DO ikl=1,kcolp
! #e5   DO ikv=1,mwp
! #e5     EqSn02(ikl,ikv) = EqSn01(ikl,ikv) + ro__SV(ikl,ikv,isn) *dzsnSV(ikl,ikv,isn) &
! #e5&                *(Cn_dSV      *(TsisSV(ikl,ikv,isn) -Tf_Sno         )&
! #e5&                 -LhfH2O      *(1.              -eta_SV(ikl,ikv,isn)))
! #e5   END DO
! #e5   END DO
! #e5 END DO
 
 
! OUTPUT/Verification: * Mass Conservation
! #m1   DO ikl=1,kcolp
! #m1   DO ikv=1,mwp
! #m1     SIsubl(ikl,ikv) = dt__SV*HLs_sv(ikl,ikv)*min(isnoSV(ikl,ikv),1)  &
! #m1&                        /Lx_H2O(ikl,ikv)
! #m1     SIrnof(ikl,ikv) = rusnSV(ikl,ikv) + RnofSV(ikl,ikv) * dt__SV &
! #m1&                - SIrnof(ikl,ikv)
! #m1   END DO
! #m1   END DO
 
 
! Anticipated Disappearance of a rapidly Melting too thin Last Snow Layer
! =======================================================================
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
        LastOK = min(1   , max(0   ,iiceSV(ikl,ikv)-isnoSV(ikl,ikv)+2) &
     &                    *min(1   ,isnoSV(ikl,ikv)-iiceSV(ikl,ikv))   &
     &                    +min(1   ,isnoSV(ikl,ikv))              )
        RapdOK = max(zer0,sign(un_1,dzMelt(ikl,ikv)-eps6         ))
        ThinOK = max(zer0,sign(un_1,dz_Min     -dzsnSV(ikl,ikv,1)))
        z_Melt = LastOK     *RapdOK*ThinOK
        noSnow(ikl,ikv)   = noSnow(ikl,ikv)+int(z_Melt)
        z_Melt        =                 z_Melt *dzsnSV(ikl,ikv,1)
        dzsnSV(ikl,ikv,1) = dzsnSV(ikl,ikv,1) - z_Melt
        EExcsv(ikl,ikv)   = EExcsv(ikl,ikv)   - z_Melt *ro__SV(ikl,ikv,1)  &
     &                                *(1.     -eta_SV(ikl,ikv,1))*LhfH2O
 
! Water  Production
! ^^^^^^^^^^^^^^^^^
        drr_SV(ikl,ikv)   = drr_SV(ikl,ikv)                            &
     &                + ro__SV(ikl,ikv,1) * z_Melt /dt__SV
      END DO
      END DO
 
 
! Update Nb of Layers
! ===================
 
! OUTPUT in SISVAT for ikl = 1 (preferably for Stand Alone Version)
! OUTPUT           for SnowFall and Snow Buffer
! #s2   IF          (isnoSV(1,1) .GT. 0)                               &
! #s2&  write(6,6005)noSnow(1,1)
! #s2   6005   format(i3,' (noSnow) ')
 
      DO ikl=1,kcolp
      DO ikv=1,mwp
        isnoSV(ikl,ikv)   = isnoSV(ikl,ikv)                            &
     &         * min(1,iabs(isnoSV(ikl,ikv)-noSnow(ikl,ikv)))
      END DO
      END DO
 
 
! OUTPUT/Verification: Energy Conservation: Energy Budget (OUT)
! #e1   DO ikl=1,kcolp
! #e1   DO ikv=1,mwp
! #e1     EqSn_1(ikl,ikv) = 0.
! #e1   END DO
! #e1   END DO
! #e1 DO   isn=nsnow,1,-1
! #e1   DO ikl=1,kcolp
! #e1   DO ikv=1,mwp
! #e1     EqSn_1(ikl,ikv) = EqSn_1(ikl,ikv) + ro__SV(ikl,ikv,isn) *dzsnSV(ikl,ikv,isn) &
! #e1&                *(Cn_dSV      *(TsisSV(ikl,ikv,isn) -Tf_Sno         )&
! #e1&                 -LhfH2O      *(1.              -eta_SV(ikl,ikv,isn)))
! #e1   END DO
! #e1   END DO
! #e1 END DO
 
 
! OUTPUT/Verification: Energy/Water Budget: Water  Budget (OUT)
! #e5   DO ikl=1,kcolp
! #e5   DO ikv=1,mwp
! #e5     WqSn_0(ikl,ikv) = WqSn_0(ikl,ikv)                            &
! #e5&                + HLs_sv(ikl,ikv)    * dt__SV                    &
! #e5&             *min(isnoSV(ikl,ikv),1) / Lx_H2O(ikl,ikv)
! #e5     WqSn_1(ikl,ikv) = drr_SV(ikl,ikv)    * dt__SV                &
! #e5&                + rusnSV(ikl,ikv)                                &
! #e5&                + RnofSV(ikl,ikv)    * dt__SV
! #e5   END DO
! #e5   END DO
! #e5 DO   isn=nsnow,1,-1
! #e5   DO ikl=1,kcolp
! #e5   DO ikv=1,mwp
! #e5     WqSn_1(ikl,ikv) = WqSn_1(ikl,ikv)                            &
! #e5&                + ro__SV(ikl,ikv,isn)* dzsnSV(ikl,ikv,isn)
! #e5   END DO
! #e5   END DO
! #e5 END DO
 
 
! OUTPUT/Verification: Energy/Water Budget
! #e5 IF (.NOT.emopen)                                              THEN
! #e5          emopen = .true.
! #e5          open(unit=43,status='unknown',file='SISVAT_qSn.vm')
! #e5          rewind 43
! #e5   write(43,43)
! #e5 43 format('SubRoutine SISVAT_qSn: Local Energy and Water Budgets',&
! #e5&       /,'=====================================================')
! #e5 END IF
! #e5 DO ikl=1,kcolp
! #e5 DO ikv=1,mwp
! #e5 IF (EqSn01(ikl,ikv).gt.1.e-3) write(43,431) dahost,EqSn01(ikl,ikv)
! #e5 431  format(' WARNING (1) in _qSn,',         a18,                &
! #e5&       ': Energy Unbalance in Phase Change = ',e15.6)
! #e5 END DO
! #e5 END DO
! #e5 DO ikl=1,kcolp
! #e5 DO ikv=1,mwp
! #e5 IF (EqSn02(ikl,ikv).gt.1.e-3) write(43,432) dahost,EqSn01(ikl,ikv)
! #e5 432  format(' WARNING (2) in _qSn,',         a18,                &
! #e5&       ': Energy Unbalance in Phase Change = ',e15.6)
! #e5 END DO
! #e5 END DO
! #e5         timeer=timeer + dt__SV
! #e5         hourer=3600.0
! #e5 IF (mod(no_err,11).eq.0)                                      THEN
! #e5         no_err=       1
! #e5   write(43,435)timeer/hourer
! #e5 435    format(11('-'),'----------+-',3('-'),'----------+-',      &
! #e5&          3('-'),'----------+-',3('-'),'----------+-',           &
! #e5&                 '----------------+----------------+',           &
! #e5&       /,f8.2,3x,'EqSn_0(1) | ',3x,'EqSn_d(1) | ',               &
! #e5&              3x,'EqSn_1(1) | ',3x,'EExcsv(1) | ',               &
! #e5&                 'E_0+E_d-E_1-EE  |   Water Budget |',           &
! #e5&       /,11('-'),'----------+-',3('-'),'----------+-',           &
! #e5&          3('-'),'----------+-',3('-'),'----------+-',           &
! #e5&                 '----------------+----------------+')
! #e5 END IF
! #e5 IF (abs(EqSn_0(1)+EqSn_d(1)-EqSn_1(1)-EExcsv(1)).gt.eps6.OR.     &
! #e5&    abs(WqSn_1(1)-WqSn_0(1))                    .gt.eps6    ) THEN
! #e5         no_err=no_err+1
! #e5   write(43,436)   EqSn_0(1),EqSn_d(1)                            &
! #e5&                 ,EqSn_1(1),EExcsv(1)                            &
! #e5&                 ,EqSn_0(1)+EqSn_d(1)-EqSn_1(1)-EExcsv(1)        &
! #e5&                 ,WqSn_1(1)-WqSn_0(1)
! #e5 436    format(8x,f12.0,' + ',f12.0,' - ',f12.0,' - ',f12.3,' = ' &
! #e5&             ,f12.3,'    | ',f15.9)
! #e5 END IF
 
! OUTPUT/Verification: Energy Conservation
! #e1   DO ikl=1,kcolp
! #e1   DO ikv=1,mwp
! #e1     EqSn_d(ikl,ikv) = EqSn_d(ikl,ikv) - EExcsv(ikl,ikv)
! #e1   END DO
! #e1   END DO
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
      IF (FlagDALLOC)                                             THEN !
 
! OUTPUT/Verification: Energy/Water Budget
! #e1 deallocate          ( EqSn_d )                                   ! Energy in Excess, initial
! #e1 deallocate          ( EqSn_0 )                                   ! Snow Energy, befor Phase Change
! #e5 deallocate          ( EqSn01 )                                   ! Snow Energy, after Phase Change
! #e5 deallocate          ( EqSn02 )                                   ! Snow Energy, after Phase Change
                                                                       !              .AND. Last Melting
! #e1 deallocate          ( EqSn_1 )                                   ! Snow Energy, after Phase Change
                                                                       !              .AND. Mass Redistr.
! OUTPUT/Verification: * Mass Conservation
! #m1 deallocate          ( SIsubl )                                   ! Snow Deposed Mass
! #m1 deallocate          ( SImelt )                                   ! Snow Melted  Mass
! #m1 deallocate          ( SIrnof )                                   ! Local Surficial Water + Run OFF
 
      deallocate          ( noSnow )                                   ! Nb of Layers Updater
      deallocate          ( EExdum )                                   ! Energy in Excess when no Snow
      deallocate          ( dzMelt )                                   ! Melted    Thickness          [m]
 
! OUTPUT/Verification: Energy/Water Budget
! #e5 deallocate          ( WqSn_0 )                                   ! Snow Water+Forcing  Initial
! #e5 deallocate          ( WqSn_1 )                                   ! Snow Water+Forcing, Final
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SISVAT_qSn
 
 
 
      subroutine SISVAT_GSn
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_GSn                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_GSn simulates SNOW Metamorphism                    |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT /  isnoSV   = total Nb of Ice/Snow Layers                      |
!     OUTPUT:  iiceSV   = total Nb of Ice      Layers                      |
!     ^^^^^^   istoSV   = 0,...,5 :   Snow     History (see istdSV data)   |
!                                                                          |
!     INPUT:   TsisSV   : Soil/Ice Temperatures (layers -nsoil,-nsoil+1, 0)|
!     ^^^^^             & Snow     Temperatures (layers  1,2,...,nsno) [K] |
!              ro__SV   : Soil/Snow Volumic Mass                   [kg/m3] |
!              eta_SV   : Soil/Snow Water   Content                [m3/m3] |
!              slorSV   : Surface Slope                           [radian] |
!              dzsnSV   : Snow Layer        Thickness                  [m] |
!              dt__SV   : Time  Step                                   [s] |
!                                                                          |
!     INPUT /  G1snSV   : Dendricity (<0) or Sphericity (>0) of Snow Layer |
!     OUTPUT:  G2snSV   : Sphericity (>0) or Size            of Snow Layer |
!     ^^^^^^                                                               |
!                                                                          |
!     Formalisme adopte pour la Representation des Grains:                 |
!     Formalism         for the Representation of  Grains:                 |
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                 |
!                                                                          |
!               1       - -1                 Neige Fraiche                 |
!              / \      |                    -------------                 |
!             /   \     |  Dendricite        decrite  par Dendricite       |
!            /     \    |  Dendricity                  et Sphericite       |
!           /       \   |                                                  |
!          2---------3  -  0                 described by Dendricity       |
!                                                     and Sphericity       |
!          |---------|                                                     |
!          0         1                                                     |
!          Sphericite                                                      |
!          Sphericity                                                      |
!                                                                          |
!          4---------5  -                                                  |
!          |         |  |                                                  |
!          |         |  |  Diametre (1/10eme de mm) (ou Taille)            |
!          |         |  |  Diameter (1/10th  of mm) (or Size  )            |
!          |         |  |                                                  |
!          |         |  |                    Neige non dendritique         |
!          6---------7  -                    ---------------------         |
!                                            decrite  par Sphericite       |
!                                                      et     Taille       |
!                                            described by Sphericity       |
!                                                     and       Size       |
!                                                                          |
!     Les Variables du Modele:                                             |
!     Model         Variables:                                             |
!     ^^^^^^^^^^^^^^^^^^^^^^^^                                             |
!       Cas Dendritique               Cas non Dendritique                  |
!                                                                          |
!       G1snSV        : Dendricite    G1snSV        : Sphericite           |
!       G2snSV        : Sphericite    G2snSV        : Taille (1/10e mm)    |
!                                                     Size                 |
!                                                                          |
!     Cas Dendritique/ Dendritic Case                                      |
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                      |
!     Dendricite(Dendricity) G1snSV                                        |
!              varie     de -G1_dSV (-99 par defaut / etoile)          a 0 |
!              division par -G1_dSV pour obtenir des valeurs entre 1  et 0 |
!              varies  from -G1_dSV (default -99    / fresh snow)     to 0 |
!              division  by -G1_dSV to obtain values       between 1 and 0 |
!                                                                          |
!     Sphericite(Sphericity) G2snSV                                        |
!              varie     de  0         (cas completement anguleux)         |
!                         a  G1_dSV (99 par defaut, cas spherique)         |
!              division par  G1_dSV pour obtenir des valeurs entre 0  et 1 |
!              varies  from  0      (full faceted)               to G1_dSV |
!                                                                          |
!     Cas non Dendritique / non Dendritic Case                             |
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                             |
!     Sphericite(Sphericity) G1snSV                                        |
!              varie     de  0         (cas completement anguleux)         |
!                         a  G1_dSV (99 par defaut, cas spherique)         |
!              division par  G1_dSV pour obtenir des valeurs entre 0  et 1 |
!              varies  from  0      (full faceted)               to G1_dSV |
!                                                                          |
!     Taille    (Size)       G2snSV                                        |
!              superieure a  ADSdSV (.4 mm) et ne fait que croitre         |
!              greater than  ADSdSV (.4 mm) always increases               |
!                                                                          |
!     Exemples: Points caracteristiques des Figures ci-dessus              |
!     ^^^^^^^^^                                                            |
!                                                                          |
!                 G1snSV    G2snSV     dendricite  sphericite  taille      |
!                                      dendricity  sphericity  size        |
!     ------------------------------------------------------------------   |
!                                                              [1/10 mm]   |
!       1        -G1_dSV    sph3SN          1           0.5                |
!       2           0         0             0           0                  |
!       3           0       G1_dSV          0           1                  |
!       4           0       ADSdSV                      0       4.         |
!       5         G1_dSV    ADSdSV-vsphe1               1       3.         |
!       6           0         --                        0       --         |
!       7         G1_dSV      --                        1       --         |
!                                                                          |
!     par defaut: G1_dSV=99.                                               |
!                           sph3SN=50.                                     |
!                           ADSdSV= 4.                                     |
!                                  vsphe1=1.                               |
!                                                                          |
!     Methode:                                                             |
!     ^^^^^^^^                                                             |
!     1. Evolution Types de Grains selon Lois de Brun et al. (1992):       |
!        Grain metamorphism according to         Brun et al. (1992):       |
!        Plusieurs Cas sont a distiguer  /  the different Cases are:       |
!         1.1 Metamorphose Neige humide  /  wet Snow                       |
!         1.2 Metamorphose Neige seche   /  dry Snow                       |
!           1.2.1 Gradient faible        /  low      Temperature Gradient  |
!           1.2.2 Gradient moyen         /  moderate Temperature Gradient  |
!           1.2.3 Gradient fort          /  high     Temperature Gradient  |
!        Dans chaque Cas on separe Neige Dendritique et non Dendritique    |
!                             le Passage Dendritique -> non Dendritique    |
!                             se fait lorsque  G1snSV devient > 0          |
!        the Case of Dentritic or non Dendritic Snow is treated separately |
!        the Limit   Dentritic -> non Dendritic is reached when G1snSV > 0 |
!                                                                          |
!     2. Tassement: Loi de Viscosite adaptee selon le Type de Grains       |
!        Snow Settling:    Viscosity depends on the   Grain Type           |
!                                                                          |
!     3. Update Variables historiques (cas non dendritique seulement)      |
!        nhSNow defaut                                                     |
!                  0    Cas normal                                         |
!        istdSV(1) 1    Grains anguleux / faceted cristal                  |
!        istdSV(2) 2    Grains ayant ete en presence d eau liquide         |
!                       mais n'ayant pas eu de caractere anguleux    /     |
!                       liquid water and no faceted cristals before        |
!        istdSV(3) 3    Grains ayant ete en presence d eau liquide         |
!                       ayant eu auparavant un caractere anguleux    /     |
!                       liquid water and    faceted cristals before        |
!                                                                          |
!     REFER. : Brun et al.      1989, J. Glaciol 35 pp. 333--342           |
!     ^^^^^^^^ Brun et al.      1992, J. Glaciol 38 pp.  13-- 22           |
!              (CROCUS Model, adapted to MAR at CEN by H.Gallee)           |
!                                                                          |
!     REFER. : Marbouty, D.     1980, J. Glaciol 26 pp. xxx--xxx           |
!     ^^^^^^^^ (CROCUS Model, adapted to MAR at CEN by H.Gallee)           |
!              (for angular shapes)                                        |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_GSn.vp        | #vp: OUTPUT/Verification: Snow   Properties   |
!                          |      unit 47, SubRoutines SISVAT_zSn, _GSn    |
!   # stdout               | #vs: OUTPUT/Verification: Snow   Properties   |
!                          |      unit  6, SubRoutine  SISVAT_GSn          |
!                                                                          |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_ctr
      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
      use Mod_SISVAT_GSn
 
 
 
! Local  Variables
! ================
 
      use Mod_SISVATLGSn
 
 
      IMPLICIT NONE
 
 
! OUTPUT
! ------
 
! #vp integer           ::  k
 
      logical           ::  vector = .true.                            ! Vectorization  Switch
      integer           ::  ikl,ikv                                    !
      integer           ::  isn   ,isnp                                !
      integer           ::  istoOK                                     !
      real(kind=real8)  ::  G1_bak,G2_bak                              ! Old Values of G1, G2
      real(kind=real8)  ::  dTsndz                                     ! Temperature Gradient
      real(kind=real8)  ::  sWater                                     !        Water Content       [%]
      real(kind=real8)  ::  exp1Wa                                     !
      real(kind=real8)  ::  dDENDR                                     ! Dendricity Increment
      real(kind=real8)  ::  DENDRn                                     ! Normalized Dendricity
      real(kind=real8)  ::  SPHERn                                     ! Normalized Sphericity
      real(kind=real8)  ::  Wet_OK                                     ! Wet Metamorphism Switch
      real(kind=real8)  ::  OK__DE                                     !
      real(kind=real8)  ::  OK__wd                                     ! New G*, from wet Dendritic
      real(kind=real8)  ::  G1__wd                                     ! New G1, from wet Dendritic
      real(kind=real8)  ::  G2__wd                                     ! New G2, from wet Dendritic
      real(kind=real8)  ::  OKlowT                                     !
      real(kind=real8)  ::  facVap                                     !
      real(kind=real8)  ::  OK_ldd                                     !
      real(kind=real8)  ::  G1_ldd                                     !
      real(kind=real8)  ::  G2_ldd                                     !
      real(kind=real8)  ::  DiamGx                                     !
      real(kind=real8)  ::  DiamOK                                     !
      real(kind=real8)  ::  No_Big                                     !
      real(kind=real8)  ::  dSPHER                                     !
      real(kind=real8)  ::  SPHER0                                     !
      real(kind=real8)  ::  SPHbig                                     !
      real(kind=real8)  ::  G1_lds                                     !
      real(kind=real8)  ::  OK_mdT                                     !
      real(kind=real8)  ::  OKmidT                                     !
      real(kind=real8)  ::  OKhigT                                     !
      real(kind=real8)  ::  OK_mdd                                     !
      real(kind=real8)  ::  G1_mdd                                     !
      real(kind=real8)  ::  G2_mdd                                     !
      real(kind=real8)  ::  G1_mds                                     !
      real(kind=real8)  ::  OK_hdd                                     !
      real(kind=real8)  ::  G1_hdd                                     !
      real(kind=real8)  ::  G2_hdd                                     !
      real(kind=real8)  ::  OK_hds                                     !
      real(kind=real8)  ::  G1_hds                                     !
      real(kind=real8)  ::  T1__OK,T2__OK                              !
      real(kind=real8)  ::  T3_xOK,T3__OK,T3_nOK                       !
      real(kind=real8)  ::  ro1_OK,ro2_OK                              !
      real(kind=real8)  ::  dT1_OK,dT2_OK,dT3xOK,dT3_OK                !
      real(kind=real8)  ::  dT4xOK,dT4_OK,dT4nOK,AngSno                !
      real(kind=real8)  ::  G2_hds,SphrOK,HISupd                       !
      real(kind=real8)  ::  H1a_OK,H1b_OK,H1__OK                       !
      real(kind=real8)  ::  H23aOK,H23bOK,H23_OK                       !
      real(kind=real8)  ::  H2__OK,H3__OK                              !
      real(kind=real8)  ::  H45_OK,H4__OK,H5__OK                       !
      real(kind=real8)  ::  ViscSn,OK_Liq,OK_Ang,OKxLiq                !
      real(kind=real8)  ::  dSnMas,dzsnew,rosnew,rosmax                !
 
      real(kind=real8)  ::  epsi5  = 1.0e-5                            ! Alpha ev67 single precision	
      real(kind=real8)  ::  epsi15 = 1.0e-15                           ! Minimal 'dry' Sphericity
!     real(kind=real8)  ::  vdiam1 = 4.0                               ! Small Grains Min.Diam.[.0001m]
      real(kind=real8)  ::  vdiam2 = 0.5                               ! Spher.Variat.Max Diam.    [mm]
      real(kind=real8)  ::  vdiam3 = 3.0                               ! Min.Diam.|Limit Spher.    [mm]
      real(kind=real8)  ::  vdiam4 = 2.0                               ! Min.Diam.|Viscosity Change
      real(kind=real8)  ::  vsphe1 = 1.0                               ! Max Sphericity
      real(kind=real8)  ::  vsphe2 = 1.0e9                             ! Low T Metamorphism  Coeff.
      real(kind=real8)  ::  vsphe3 = 0.5                               ! Max.Sphericity (history=1)
      real(kind=real8)  ::  vsphe4 = 0.1                               ! Min.Sphericity=>history=1
 
! DATA (Coefficient Fonction fort Gradient Marbouty)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(kind=real8)  ::  vtang1 = 40.0                              ! Temperature Contribution v
      real(kind=real8)  ::  vtang2 =  6.0                              !
      real(kind=real8)  ::  vtang3 = 22.0                              !
      real(kind=real8)  ::  vtang4 =  0.7                              !
      real(kind=real8)  ::  vtang5 =  0.3                              !
      real(kind=real8)  ::  vtang6 =  6.0                              !
      real(kind=real8)  ::  vtang7 =  1.0                              !
      real(kind=real8)  ::  vtang8 =  0.8                              !
      real(kind=real8)  ::  vtang9 = 16.0                              !
      real(kind=real8)  ::  vtanga =  0.2                              !
      real(kind=real8)  ::  vtangb =  0.2                              !
      real(kind=real8)  ::  vtangc = 18.0                              ! Temperature Contribution ^
 
      real(kind=real8)  ::  vrang1 =  0.40                             ! Density     Contribution v
      real(kind=real8)  ::  vrang2 =  0.15                             ! Density     Contribution ^
 
      real(kind=real8)  ::  vgang1 =  0.70                             ! Grad(T)     Contribution v
      real(kind=real8)  ::  vgang2 =  0.25                             !
      real(kind=real8)  ::  vgang3 =  0.40                             !
      real(kind=real8)  ::  vgang4 =  0.50                             !
      real(kind=real8)  ::  vgang5 =  0.10                             !
      real(kind=real8)  ::  vgang6 =  0.15                             !
      real(kind=real8)  ::  vgang7 =  0.10                             !
      real(kind=real8)  ::  vgang8 =  0.55                             !
      real(kind=real8)  ::  vgang9 =  0.65                             !
      real(kind=real8)  ::  vganga =  0.20                             !
      real(kind=real8)  ::  vgangb =  0.85                             !
      real(kind=real8)  ::  vgangc =  0.15                             ! Grad(T)     Contribution ^
 
      real(kind=real8)  ::  vgran6 = 51.                               ! Max.Sphericity for Settling
      real(kind=real8)  ::  vtelv1 =  5.e-1                            ! Threshold | history = 2, 3
      real(kind=real8)  ::   vvap1 = -6.e3                             ! Vapor Pressure Coefficient
      real(kind=real8)  ::   vvap2 =  0.4                              ! Vapor Pressure Exponent
      real(kind=real8)  ::  vgrat1 =  0.05                             ! Boundary weak/mid   grad(T)
      real(kind=real8)  ::  vgrat2 =  0.15                             ! Boundary mid/strong grad(T)
      real(kind=real8)  ::     vfi =  0.09                             ! PHI,         strong grad(T)
 
      real(kind=real8)  ::  vvisc1 =   0.70                            ! Viscosity Coefficients
      real(kind=real8)  ::  vvisc2 =   1.11e5                          !
      real(kind=real8)  ::  vvisc3 =  23.00                            !
      real(kind=real8)  ::  vvisc4 =   0.10                            ! Viscosity Coefficients
      real(kind=real8)  ::  vvisc5 =   1.00                            ! Viscosity Coefficients, wet Snow
      real(kind=real8)  ::  vvisc6 =   2.00                            ! Viscosity Coefficients, wet Snow
      real(kind=real8)  ::  vvisc7 =  10.00                            ! Viscosity Coefficients, wet Snow
      real(kind=real8)  ::  rovisc =   0.25                            ! Wet Snow Density  Influence
      real(kind=real8)  ::    vdz3 =   0.30                            ! Maximum Layer Densification
 
      real(kind=real8)  ::  OK__ws                                     ! New G2
      real(kind=real8)  ::  G1__ws                                     ! New G1, from wet Spheric
      real(kind=real8)  ::  G2__ws                                     ! New G2, from wet Spheric
      real(kind=real8)  ::  husi_0 =  20.                              ! Constant for New G2:   10  * 2
      real(kind=real8)  ::  husi_1 =   0.23873                         ! Constant for New G2: (3/4) /pi
      real(kind=real8)  ::  husi_2 =   4.18880                         ! Constant for New G2: (4/3) *pi
      real(kind=real8)  ::  husi_3 =   0.33333                         ! Constant for New G2:  1/3
      real(kind=real8)  ::  vtail1 =   1.28e-08                        ! Constant for New G2:  Wet Metamorphism
      real(kind=real8)  ::  vtail2 =   4.22e-10                        ! Constant for New G2: (NON Dendritic / Spheric)
      real(kind=real8)  ::  frac_j                                     ! Time Step            [Day]
 
      real(kind=real8)  ::   vdent1 =  2.e8                            ! Wet Snow Metamorphism
      integer           ::  nvdent1 =  3                               ! (Coefficients for
      integer           ::  nvdent2 = 16                               !           Dendricity)
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
      allocate            ( ro_dry(kcolp,mwp,nsnow) )                  ! Dry Density            [g/cm3]
      allocate            ( etaSno(kcolp,mwp,nsnow) )                  ! Liquid Water Content   [g/cm2]
      allocate            ( SnMass(kcolp,mwp)       )                  ! Snow   Mass            [kg/m2]
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
!  1. Metamorphoses dans les Strates
!     Metamorphism
!     ==============================
 
      frac_j = dt__SV / 86400.                        ! Time Step [Day]
 
 
!  1.1 Initialisation: teneur en eau liquide et gradient de temperature
!  ------------------  liquid water content and temperature gradient
 
      DO   isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
          ro_dry(ikl,ikv,isn) = 1.e-3 *ro__SV(ikl,ikv,isn)   &! Dry Density
     &                    *(1.    -eta_SV(ikl,ikv,isn))   !         [g/cm3]
          etaSno(ikl,ikv,isn) = 1.e-1 *dzsnSV(ikl,ikv,isn)   &! Liquid Water
     &                    *        ro__SV(ikl,ikv,isn)   &! Content [g/cm2]
     &                    *        eta_SV(ikl,ikv,isn)!
        END DO
        END DO
      END DO
 
      DO   isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
          isnp   = min(isn+1,isnoSV(ikl,ikv))
 
          dTsndz = abs( (TsisSV(ikl,ikv,isnp)-TsisSV(ikl,ikv,isn-1)) *2.e-2&
     &            /max(((dzsnSV(ikl,ikv,isnp)+dzsnSV(ikl,ikv,isn)  )   &
     &                 *(           isnp -           isn)              &
     &                 +(dzsnSV(ikl,ikv,isn )+dzsnSV(ikl,ikv,isn-1))),eps6))
!         Factor 1.d-2 for Conversion K/m --> K/cm
 
 
!  1.2 Metamorphose humide
!      Wet Snow Metamorphism
!      ---------------------
 
          Wet_OK = max(zer0,sign(un_1,eta_SV(ikl,ikv,isn)-eps6))
 
 
!      Vitesse de diminution de la dendricite
!      Rate of the dendricity decrease
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          sWater=1.d-1*ro__SV(ikl,ikv,isn)*eta_SV(ikl,ikv,isn)         &
     &       /max(eps6,ro_dry(ikl,ikv,isn))
! .    sWater:Water Content [%]
!             1.d-1= 1.d2(1->%) * 1.d-3(ro__SV*eta_SV:kg/m3->g/cm3)
 
          exp1Wa=   sWater**nvdent1
          dDENDR=max(exp1Wa/nvdent2,vdent1*exp(vvap1/Tf_Sno))
 
!  1.2.1 Cas dendritique/dendritic Case
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          OK__wd=max(zer0,                           &!
     &               sign(un_1,-G1snSV(ikl,ikv,isn)  &!
     &                         -eps6           ))     !
 
          DENDRn=-G1snSV(ikl,ikv,isn)/G1_dSV  ! Normalized Dendricity (+)
          SPHERn= G2snSV(ikl,ikv,isn)/G1_dSV  ! Normalized Sphericity
          DENDRn= DENDRn -dDENDR *frac_j  ! New        Dendricity (+)
          SPHERn= SPHERn +dDENDR *frac_j  ! New        Sphericity
 
          OK__DE=max(zer0,                           &! IF 1.,
     &               sign(un_1, DENDRn               &! NO change
     &                         -eps6           ))     ! Dendr. -> Spheric
 
          G1__wd=OK__DE *    (      -DENDRn*G1_dSV)  &! Dendritic
     &      +(1.-OK__DE)* min(G1_dSV,SPHERn*G1_dSV)   ! Dendr. -> Spheric
          G2__wd=OK__DE * min(G1_dSV,SPHERn*G1_dSV)  &! Spheric
     &      +(1.-OK__DE)*(ADSdSV-min(SPHERn,vsphe1))  ! Spher. -> Size
 
!  1.2.2 Cas non dendritique non completement spherique
!        Evolution de la Sphericite seulement.
!        Non dendritic and not completely spheric Case
!        Evolution of    Sphericity only (not size)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          OK__ws=max(zer0,                           &!
     &               sign(un_1, G1_dSV               &!
     &                         -epsi5                &!
     &                         -G1snSV(ikl,ikv,isn))) !
 
          SPHERn= G1snSV(ikl,ikv,isn)/G1_dSV
          SPHERn= SPHERn +dDENDR *frac_j
          G1__ws=         min(G1_dSV,SPHERn*G1_dSV)
 
!  1.2.3 Cas non dendritique et spherique / non dendritic and spheric
!        Evolution de la Taille seulement / Evolution of Size only
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          G2__ws =  husi_0                                             &
     &           *( husi_1                                             &
     &            *(husi_2 *( G2snSV(ikl,ikv,isn)/husi_0)**3           &
     &                      +(vtail1 +vtail2 *exp1Wa    )*dt__SV))     &
     &           ** husi_3
 
 
!  1.3 Metamorposes seches / Dry Metamorphism
!      --------------------------------------
 
 
!  1.3.1 Calcul Metamorphose faible/low Gradient (0.00-0.05 deg/cm)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          OKlowT=max(zer0,                           &!
     &               sign(un_1, vgrat1               &!
     &                         -dTsndz         ))     !
 
          facVap=exp(vvap1/TsisSV(ikl,ikv,isn))
 
!  1.3.1.1 Cas dendritique / dendritic Case
 
          OK_ldd=max(zer0,                           &!
     &               sign(un_1,-G1snSV(ikl,ikv,isn)  &!
     &                         -eps6           ))     !
 
          DENDRn=-G1snSV(ikl,ikv,isn)     /G1_dSV
          SPHERn= G2snSV(ikl,ikv,isn)     /G1_dSV
          DENDRn= DENDRn-vdent1*facVap*frac_j
          SPHERn= SPHERn+vsphe2*facVap*frac_j
 
          OK__DE=max(zer0,                           &! IF 1.,
     &               sign(un_1, DENDRn               &! NO change
     &                         -eps6           ))     ! Dendr. -> Spheric
 
          G1_ldd= OK__DE *    (      -DENDRn*G1_dSV) &! Dendritic
     &       +(1.-OK__DE)* min(G1_dSV,SPHERn*G1_dSV)  ! Dendr. -> Spheric
          G2_ldd= OK__DE * min(G1_dSV,SPHERn*G1_dSV) &! Spheric
     &       +(1.-OK__DE)*(ADSdSV-min(SPHERn,vsphe1)) ! Spher. -> Size
 
!  1.3.1.2 Cas non dendritique / non dendritic Case
 
          SPHERn=G1snSV(ikl,ikv,isn)/G1_dSV
          DiamGx=G2snSV(ikl,ikv,isn)*0.1
 
          istoOK=min( abs(istoSV(ikl,ikv,isn)-       &!
     &                    istdSV(1)      ),1)         ! zero if istoSV = 1
          DiamOK=max(zer0,  sign(un_1,vdiam2-DiamGx))
          No_Big=    istoOK+DiamOK
          No_Big=min(No_Big,un_1)
 
          dSPHER=           vsphe2*facVap*frac_j      !
          SPHER0=    SPHERn+dSPHER                    ! small grains
          SPHbig=    SPHERn+dSPHER                   &! big   grains
     &        *exp(min(zer0,vdiam3-G2snSV(ikl,ikv,isn)))  ! (history = 2 or 3)
          SPHbig=       min(vsphe3,SPHbig)            ! limited sphericity
          SPHERn= No_Big *  SPHER0                   &!
     &      + (1.-No_Big)*  SPHbig                    !
 
! HG v Precudes underflow
          SPHERn=       max(epsi15,SPHERn)            !
! HG ^ Precudes underflow
 
          G1_lds=       min(G1_dSV,SPHERn*G1_dSV)     !
 
!  1.3.2 Calcul Metamorphose Gradient Moyen/Moderate (0.05-0.15)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          OK_mdT=max(zer0,                           &!
     &               sign(un_1, vgrat2               &!
     &                         -dTsndz))              !
          OKmidT=               OK_mdT  *(1.-OKlowT)  !
          OKhigT=          (1. -OK_mdT) *(1.-OKlowT)  !
 
          facVap=vdent1*exp(vvap1/TsisSV(ikl,ikv,isn))   &!
     &                 *   (1.e2 *dTsndz)**vvap2      !
 
!  1.3.2.1 cas dendritique / dendritic case.
 
          OK_mdd=max(zer0,                           &!
     &               sign(un_1,-G1snSV(ikl,ikv,isn)  &!
     &                         -eps6           ))     !
 
          DENDRn=-G1snSV(ikl,ikv,isn)/G1_dSV
          SPHERn= G2snSV(ikl,ikv,isn)/G1_dSV
          DENDRn= DENDRn - facVap*frac_j
          SPHERn= SPHERn - facVap*frac_j
 
          OK__DE=max(zer0,                           &! IF 1.,
     &               sign(un_1, DENDRn               &! NO change
     &                         -eps6           ))     ! Dendr. -> Spheric
 
          G1_mdd= OK__DE *    (      -DENDRn*G1_dSV) &! Dendritic
     &       +(1.-OK__DE)* max(zer0  ,SPHERn*G1_dSV)  ! Dendr. -> Spheric
          G2_mdd= OK__DE * max(zer0  ,SPHERn*G1_dSV) &! Spheric
     &       +(1.-OK__DE)*(ADSdSV-max(SPHERn,zer0  )) ! Spher. -> Size
 
!  1.3.2.2 Cas non dendritique / non dendritic Case
 
          SPHERn=G1snSV(ikl,ikv,isn)/G1_dSV
          SPHERn=         SPHERn-facVap*frac_j
          G1_mds=max(zer0,SPHERn*G1_dSV)
 
!  1.3.3 Calcul Metamorphose fort / high Gradient
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          facVap=vdent1*exp(vvap1/TsisSV(ikl,ikv,isn))   &!
     &                 *   (1.e2 *dTsndz)**vvap2      !
 
!  1.3.3.1 Cas dendritique / dendritic Case
 
          OK_hdd=max(zer0,                           &!
     &               sign(un_1,-G1snSV(ikl,ikv,isn)  &!
     &                         -eps6           ))     !
 
          DENDRn=-G1snSV(ikl,ikv,isn)/G1_dSV          !
          SPHERn= G2snSV(ikl,ikv,isn)/G1_dSV          !
          DENDRn= DENDRn - facVap*frac_j              !
          SPHERn= SPHERn - facVap*frac_j              ! Non dendritic
                                                      ! and angular
          OK__DE=max(zer0,                           &! IF 1.,
     &               sign(un_1, DENDRn               &! NO change
     &                         -eps6  ))              ! Dendr. -> Spheric
 
          G1_hdd= OK__DE *    (      -DENDRn*G1_dSV) &! Dendritic
     &       +(1.-OK__DE)* max(zer0  ,SPHERn*G1_dSV)  ! Dendr. -> Spheric
          G2_hdd= OK__DE * max(zer0  ,SPHERn*G1_dSV) &! Spheric
     &       +(1.-OK__DE)*(ADSdSV-max(SPHERn,zer0  )) ! Spher. -> Size
 
!  1.3.3.2 Cas non dendritique non completement anguleux.
!          non dendritic and spericity gt. 0
 
          OK_hds=max(zer0,                           &!
     &               sign(un_1, G1snSV(ikl,ikv,isn)  &!
     &                         -eps6           ))     !
 
          SPHERn= G1snSV(ikl,ikv,isn)/G1_dSV
          SPHERn= SPHERn - facVap*frac_j
          G1_hds= max(zer0,SPHERn*G1_dSV)
 
!  1.3.3.3 Cas non dendritique et anguleux
!          dendritic and spericity = 0.
 
          T1__OK = max(zer0,sign(un_1,TsisSV(ikl,ikv,isn)-Tf_Sno+vtang1))
          T2__OK = max(zer0,sign(un_1,TsisSV(ikl,ikv,isn)-Tf_Sno+vtang2))
          T3_xOK = max(zer0,sign(un_1,TsisSV(ikl,ikv,isn)-Tf_Sno+vtang3))
          T3__OK =                    T3_xOK  * (1. - T2__OK)
          T3_nOK =              (1. - T3_xOK) * (1. - T2__OK)
          ro1_OK = max(zer0,sign(un_1,vrang1-ro_dry(ikl,ikv,isn)))
          ro2_OK = max(zer0,sign(un_1,ro_dry(ikl,ikv,isn)-vrang2))
          dT1_OK = max(zer0,sign(un_1,vgang1-dTsndz         ))
          dT2_OK = max(zer0,sign(un_1,vgang2-dTsndz         ))
          dT3xOK = max(zer0,sign(un_1,vgang3-dTsndz         ))
          dT3_OK =                    dT3xOK  * (1. - dT2_OK)
          dT4xOK = max(zer0,sign(un_1,vgang4-dTsndz         ))
          dT4_OK =                    dT4xOK  * (1. - dT3_OK)          &
     &                                        * (1. - dT2_OK)
          dT4nOK =              (1. - dT4xOK) * (1. - dT3_OK)          &
     &                                        * (1. - dT2_OK)
 
!  Influence de la Temperature /Temperature Influence
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          AngSno =                                                     &
     &      T1__OK                                                     & ! 11
     &    *(T2__OK*(vtang4+vtang5*(Tf_Sno       -TsisSV(ikl,ikv,isn))  & ! 12
     &                    /vtang6)                                     & !
     &     +T3__OK*(vtang7-vtang8*(Tf_Sno-vtang2-TsisSV(ikl,ikv,isn))  & ! 13
     &                    /vtang9)                                     & !
     &     +T3_nOK*(vtanga-vtangb*(Tf_Sno-vtang3-TsisSV(ikl,ikv,isn))  & ! 14
     &                    /vtangc))                                    & !
 
!  Influence de la Masse Volumique /Density Influence
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     &    * ro1_OK                                                     & !
     &        *(   ro2_OK*(1. - (ro_dry(ikl,ikv,isn)-vrang2)           & !
     &                                  /(vrang1-vrang2))              & !
     &         +1.-ro2_OK                                )             & !
 
!  Influence du Gradient de Temperature /Temperature Gradient Influence
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     &        *(   dT1_OK*(dT2_OK*vgang5*(dTsndz-vgang6)               & ! 15
     &                                  /(vgang2-vgang6)               & !
     &                    +dT3_OK*vgang7                               & ! 16
     &                    +dT4_OK*vgang9                               & ! 17
     &                    +dT4nOK*vgangb                )              & ! 18
     &         +1.-dT1_OK                                )             & !
     &    + ro1_OK                                                     & !
     &        *    dT1_OK*(dT3_OK*vgang8*(dTsndz-vgang2)               & !
     &                                  /(vgang3-vgang2)               & !
     &                    +dT4_OK*vganga*(dTsndz-vgang3)               & !
     &                                  /(vgang4-vgang3)               & !
     &                    +dT4nOK*vgangc*(dTsndz-vgang4)               & !
     &                                  /(vgang1-vgang4))                !
 
          G2_hds = G2snSV(ikl,ikv,isn) + 1.d2 *AngSno*vfi     *frac_j
 
 
! New Properties
! --------------
 
          G1_bak          = G1snSV(ikl,ikv,isn)
          G2_bak          = G2snSV(ikl,ikv,isn)
 
          G1snSV(ikl,ikv,isn) = Wet_OK * (    OK__wd             *G1__wd   & !  1
     &                               +(1.-OK__wd)*    OK__ws *G1__ws   & !  2
     &                               +(1.-OK__wd)*(1.-OK__ws)*G1_bak)  & !  3
     &               +(1. - Wet_OK)                                    & !
     &                *(    OKlowT  *(    OK_ldd             *G1_ldd   & !  4
     &                               +(1.-OK_ldd)            *G1_lds)  & !  5
     &                    + OKmidT  *(    OK_mdd             *G1_mdd   & !  6
     &                               +(1.-OK_mdd)            *G1_mds)  & !  7
     &                    + OKhigT  *(    OK_hdd             *G1_hdd   & !  8
     &                               +(1.-OK_hdd)*    OK_hds *G1_hds   & !  9
     &                               +(1.-OK_hdd)*(1.-OK_hds)*G1_bak))   ! 10
 
! XF
          IF (G1snSV(ikl,ikv,isn).GE.0.0.AND.G1snSV(ikl,ikv,isn)<0.1)  & !
! HG -------------------------^^^^^^
     &        G2_hds = G2snSV(ikl,ikv,isn) + 1.d1 *AngSno*vfi    *frac_j !
!             Previens chute exageree de l'albedo lorsque G1~0.
! XF
 
          G2snSV(ikl,ikv,isn) = Wet_OK * (    OK__wd             *G2__wd   & !  1
     &                               +(1.-OK__wd)*    OK__ws *G2_bak   & !  2
     &                               +(1.-OK__wd)*(1.-OK__ws)*G2__ws)  & !  3
     &               +(1. - Wet_OK)                                    & !
     &                *(    OKlowT  *(    OK_ldd             *G2_ldd   & !  4
     &                               +(1.-OK_ldd)            *G2_bak)  & !  5
     &                    + OKmidT  *(    OK_mdd             *G2_mdd   & !  6
     &                               +(1.-OK_mdd)            *G2_bak)  & !  7
     &                    + OKhigT  *(    OK_hdd             *G2_hdd   & !  8
     &                               +(1.-OK_hdd)*    OK_hds *G2_bak   & !  9
     &                               +(1.-OK_hdd)*(1.-OK_hds)*G2_hds))   ! 10
 
! OUTPUT/Verification: Snow Layers Agregation: Properties
! #vp     G_curr( 1) =     Wet_OK             *    OK__wd
! #vp     G_curr( 2) =     Wet_OK             *(1.-OK__wd)*    OK__ws
! #vp     G_curr( 3) =     Wet_OK             *(1.-OK__wd)*(1.-OK__ws)
! #vp     G_curr( 4) = (1.-Wet_OK)*    OKlowT *    OK_ldd
! #vp     G_curr( 5) = (1.-Wet_OK)*    OKlowT *(1.-OK_ldd)
! #vp     G_curr( 6) = (1.-Wet_OK)*    OKmidT *    OK_mdd
! #vp     G_curr( 7) = (1.-Wet_OK)*    OKmidT *(1.-OK_mdd)
! #vp     G_curr( 8) = (1.-Wet_OK)*    OKhigT *    OK_hdd
! #vp     G_curr( 9) = (1.-Wet_OK)*    OKhigT *(1.-OK_hdd)*    OK_hds
! #vp     G_curr(10) = (1.-Wet_OK)*    OKhigT *(1.-OK_hdd)*(1.-OK_hds)
! #vp     G_curr(11) =     T1__OK                         * G_curr(10)
! #vp     G_curr(12) =     T2__OK                         * G_curr(10)
! #vp     G_curr(13) =     T3__OK                         * G_curr(10)
! #vp     G_curr(14) =     T3_nOK                         * G_curr(10)
! #vp     G_curr(15) =     ro1_OK*     dT1_OK *    dT2_OK * G_curr(10)
! #vp     G_curr(16) =     ro1_OK*     dT1_OK *    dT3_OK * G_curr(10)
! #vp     G_curr(17) =     ro1_OK*     dT1_OK *    dT4_OK * G_curr(10)
! #vp     G_curr(18) =     ro1_OK*     dT1_OK *    dT4nOK * G_curr(10)
 
! #vp     Gcases( 1) = max(Gcases( 1),G_curr( 1))
! #vp     Gcases( 2) = max(Gcases( 2),G_curr( 2))
! #vp     Gcases( 3) = max(Gcases( 3),G_curr( 3))
! #vp     Gcases( 4) = max(Gcases( 4),G_curr( 4))
! #vp     Gcases( 5) = max(Gcases( 5),G_curr( 5))
! #vp     Gcases( 6) = max(Gcases( 6),G_curr( 6))
! #vp     Gcases( 7) = max(Gcases( 7),G_curr( 7))
! #vp     Gcases( 8) = max(Gcases( 8),G_curr( 8))
! #vp     Gcases( 9) = max(Gcases( 9),G_curr( 9))
! #vp     Gcases(10) = max(Gcases(10),G_curr(10))
! #vp     Gcases(11) = max(Gcases(11),G_curr(11))
! #vp     Gcases(12) = max(Gcases(12),G_curr(12))
! #vp     Gcases(13) = max(Gcases(13),G_curr(13))
! #vp     Gcases(14) = max(Gcases(14),G_curr(14))
! #vp     Gcases(15) = max(Gcases(15),G_curr(15))
! #vp     Gcases(16) = max(Gcases(16),G_curr(16))
! #vp     Gcases(17) = max(Gcases(17),G_curr(17))
! #vp     Gcases(18) = max(Gcases(18),G_curr(18))
 
! #vp     IF          (isn    .le.     isnoSV(ikl,ikv))                &
! #vp&    write(47,471)isn            ,isnoSV(ikl,ikv)                    ,&
! #vp&                 TsisSV(ikl,ikv,isn),ro__SV(ikl,ikv,isn),eta_SV(ikl,ikv,isn),&
! #vp&                 G1_bak         ,G2_bak         ,istoSV(ikl,ikv,isn),&
! #vp&                 dTsndz,                                         &
! #vp&                (       k ,k=1,18),                              &
! #vp&                (G_curr(k),k=1,18),                              &
! #vp&                (Gcases(k),k=1,18),                              &
! #vp&                 Wet_OK,OK__wd,G1__wd,G2__wd,                    &
! #vp&                     1.-OK__wd,OK__ws,G1__ws,1.-OK__ws,G2__ws,   &
! #vp&              1.-Wet_OK,OKlowT,OK_ldd,G1_ldd,          G2_ldd,   &
! #vp&                            1.-OK_ldd,G1_lds,                    &
! #vp&                        OKmidT,OK_mdd,G1_mdd,          G1_mdd,   &
! #vp&                            1.-OK_mdd,G1_mds,                    &
! #vp&                        OKhigT,OK_hdd,G1_hdd,          G2_hdd,   &
! #vp&                            1.-OK_hdd,OK_hds,          G1_hds,   &
! #vp&                                             1.-OK_hds,G2_hds,   &
! #vp&                 G1snSV(ikl,ikv,isn),                            &
! #vp&                 G2snSV(ikl,ikv,isn)
! #vp 471 format(                                                      &
! #vp&         /,' isn     =  ',i4,6x,'(MAX.:',i4,')',                 &
! #vp&         /,' T       =  ',f8.3,                                  &
! #vp&         /,' ro      =  ',f8.3,                                  &
! #vp&         /,' eta     =  ',f8.3,                                  &
! #vp&         /,' G1      =  ',f8.3,                                  &
! #vp&         /,' G2      =  ',f8.3,                                  &
! #vp&         /,' Histor. =  ',i4  ,                                  &
! #vp&         /,' Grad(T) =  ',f8.3,'                   ' ,18i3  ,    &
! #vp&/,         '                       Current    Case: ',18f3.0,    &
! #vp&/,         '                       Cases performed: ',18f3.0,    &
! #vp&/,' ------------------------------------------------------------',&
! #vp&             '-----------+------------------+------------------+',&
! #vp&/,' Status                                                      ',&
! #vp&             '           | G1               | G2               |',&
! #vp&/,' ------------------------------------------------------------',&
! #vp&             '-----------+------------------+------------------+',&
! #vp&/,'    Wet_OK: ',f8.3,'                  OK__wd: ',f8.3,'       ',&
! #vp&             '           | G1__wd: ',f8.3,' | G2__wd: ',f8.5,' |',&
! #vp&/,'                                   1.-OK__wd: ',f8.3,' OK__ws',&
! #vp&             ': ',f8.3,' | G1__ws: ',f8.3,' |                  |',&
! #vp&/,'                                                    1.-OK__ws',&
! #vp&             ': ',f8.3,' |                  | G2__ws: ',f8.5,' |',&
! #vp&/,' 1.-Wet_OK: ',f8.3,' OKlowT: ',f8.3,' OK_ldd: ',f8.3,'       ',&
! #vp&             '           | G1_ldd: ',f8.3,' | G2_ldd: ',f8.5,' |',&
! #vp&/,'                                   1.-OK_ldd: ',f8.3,'       ',&
! #vp&             '           | G1_lds: ',f8.3,' |                  |',&
! #vp&/,'                     OKmidT: ',f8.3,' OK_mdd: ',f8.3,'       ',&
! #vp&             '           | G1_mdd: ',f8.3,' | G2_mdd: ',f8.5,' |',&
! #vp&/,'                                   1.-OK_mdd: ',f8.3,'       ',&
! #vp&             '           | G1_mds: ',f8.3,' |                  |',&
! #vp&/,'                     OKhigT: ',f8.3,' OK_hdd: ',f8.3,'       ',&
! #vp&             '           | G1_hdd: ',f8.3,' | G2_hdd: ',f8.5,' |',&
! #vp&/,'                                   1.-OK_hdd: ',f8.3,' OK_hds',&
! #vp&             ': ',f8.3,' | G1_hds: ',f8.3,' |                  |',&
! #vp&/,'                                                    1.-OK_hds',&
! #vp&             ': ',f8.3,' |                  | G2_hds: ',f8.5,' |',&
! #vp&/,' ------------------------------------------------------------',&
! #vp&             '-----------+------------------+------------------+',&
! #vp&/,'                                                             ',&
! #vp&             '           |         ',f8.3,' |         ',f8.5,' |',&
! #vp&/,' ------------------------------------------------------------',&
! #vp&             '-----------+------------------+------------------+')
        END DO
        END DO
      END DO
 
 
!  2. Mise a Jour Variables Historiques (Cas non dendritique)
!     Update of the historical Variables
!     =======================================================
 
      IF (vector)                                                   THEN
        DO isn=1,nsnow
        DO ikl=1,kcolp
        DO ikv=1,mwp
 
          SphrOK = max(zer0,sign(un_1,       G1snSV(ikl,ikv,isn)))
          H1a_OK = max(zer0,sign(un_1,vsphe4-G1snSV(ikl,ikv,isn)))
          H1b_OK =     1   - min(1   ,       istoSV(ikl,ikv,isn))
          H1__OK =                    H1a_OK*H1b_OK
          H23aOK = max(zer0,sign(un_1,vsphe4-G1_dSV                    &
     &                                      +G1snSV(ikl,ikv,isn)))
          H23bOK = max(zer0,sign(un_1,etaSno(ikl,ikv,isn)              &
     &                      /max(eps6,dzsnSV(ikl,ikv,isn))             &
     &                                      -vtelv1         ))
          H23_OK =                    H23aOK*H23bOK
          H2__OK =     1   - min(1   ,       istoSV(ikl,ikv,isn))
          H3__OK =     1   - min(1   ,   abs(istoSV(ikl,ikv,isn)-istdSV(1)))
          H45_OK = max(zer0,sign(un_1,Tf_Sno-TsisSV(ikl,ikv,isn)+eps6))
          H4__OK =     1   - min(1   ,   abs(istoSV(ikl,ikv,isn)-istdSV(2)))
          H5__OK =     1   - min(1   ,   abs(istoSV(ikl,ikv,isn)-istdSV(3)))
 
          HISupd          =                                            &
     &    SphrOK*(H1__OK                             *istdSV(1)        &
     &       +(1.-H1__OK)*    H23_OK         *(H2__OK*istdSV(2)        &
     &                                        +H3__OK*istdSV(3))       &
     &       +(1.-H1__OK)*(1.-H23_OK) *H45_OK*(H4__OK*istdSV(4)        &
     &                                        +H5__OK*istdSV(5)))
          istoSV(ikl,ikv,isn)=int(HISupd) +                            &
     &        int(1.-min(un_1,HISupd))               *istoSV(ikl,ikv,isn)
        END DO
        END DO
        END DO
      ELSE
 
 
!  2. Mise a Jour Variables Historiques (Cas non dendritique)
!     Update of the historical Variables
!     =======================================================
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
        DO isn=iiceSV(ikl,ikv),isnoSV(ikl,ikv)
          IF  (G1snSV(ikl,ikv,isn).ge.0.)                           THEN
            IF(G1snSV(ikl,ikv,isn).lt.vsphe4.and.istoSV(ikl,ikv,isn).eq.0)  THEN
                   istoSV(ikl,ikv,isn)=istdSV(1)
            ELSEIF(G1_dSV-G1snSV(ikl,ikv,isn)         .lt.vsphe4.and.  &
     &             etaSno(ikl,ikv,isn)/dzsnSV(ikl,ikv,isn).gt.vtelv1)   THEN
              IF  (istoSV(ikl,ikv,isn).eq.0)                           &
     &             istoSV(ikl,ikv,isn)=   istdSV(2)
              IF  (istoSV(ikl,ikv,isn).eq.istdSV(1))                   &
     &             istoSV(ikl,ikv,isn)=   istdSV(3)
            ELSEIF(TsisSV(ikl,ikv,isn).lt.Tf_Sno)                   THEN
              IF  (istoSV(ikl,ikv,isn).eq.istdSV(2))                   &
     &             istoSV(ikl,ikv,isn)=   istdSV(4)
              IF  (istoSV(ikl,ikv,isn).eq.istdSV(3))                   &
     &             istoSV(ikl,ikv,isn)=   istdSV(5)
            END IF
          END IF
        END DO
        END DO
        END DO
      END IF
 
 
!  3. Tassement mecanique /mechanical Settlement
!     ==========================================
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          SnMass(ikl,ikv) = 0.
        END DO
        END DO
      DO isn=nsnow,1,-1
        DO ikl=1,kcolp
        DO ikv=1,mwp
          dSnMas     = 100.*dzsnSV(ikl,ikv,isn)*ro_dry(ikl,ikv,isn)
          SnMass(ikl,ikv)=      SnMass(ikl,ikv)+0.5*dSnMas
          ViscSn     =      vvisc1         *vvisc2                     &
     &               *exp(vvisc3           *ro_dry(ikl,ikv,isn)        &
     &                   +vvisc4*abs(Tf_Sno-TsisSV(ikl,ikv,isn)))      &
     &                                     *ro_dry(ikl,ikv,isn)/rovisc
 
!  Changement de Viscosite si Teneur en Eau liquide
!  Change of the Viscosity if liquid Water Content
!  ------------------------------------------------
 
          OK_Liq     =    max(zer0,sign(un_1,etaSno(ikl,ikv,isn)-eps6))
          OK_Ang     =    max(zer0,sign(un_1,vgran6-G1snSV(ikl,ikv,isn)))  &
     &                *(1-min(1   , abs(istoSV(ikl,ikv,isn)-istdSV(1))))
 
! OUTPUT/Verification: Snow   Properties
! #vs     IF (G1snSV(ikl,ikv,isn).gt.0..AND.G1snSV(ikl,ikv,isn).lt.vsphe4  &
! #vs&                             .AND.istoSV(ikl,ikv,isn).eq.     0) &
! #vs&    THEN
! #vs       write(6,*) ikl,ikv,isn,' G1,G2,hist,OK_Ang  ',             &
! #vs&          G1snSV(ikl,ikv,isn), G2snSV(ikl,ikv,isn),istoSV(ikl,ikv,isn),OK_Ang
! #vs       stop "Grains anguleux mal d?finis"
! #vs     END IF
 
          OKxLiq     =    max(zer0,sign(un_1,vtelv1-etaSno(ikl,ikv,isn)&
     &                                    /max(eps6,dzsnSV(ikl,ikv,isn)))) &
     &               *    max(0   ,sign(1   ,istoSV(ikl,ikv,isn)       &
     &                                      -istdSV(1)      ))
          ViscSn     =                                                 &
     &    ViscSn*(    OK_Liq/(vvisc5+vvisc6*etaSno(ikl,ikv,isn)        &
     &                            /max(eps6,dzsnSV(ikl,ikv,isn)))      &
     &           +(1.-OK_Liq)                               )          &
     &          *(    OK_Ang*exp(min(ADSdSV,G2snSV(ikl,ikv,isn)-vdiam4))   &
     &           +(1.-OK_Ang)                                       )  &
     &          *(    OKxLiq*        vvisc7                            &
     &           +(1.-OKxLiq)              )
 
 
!  Calcul nouvelle Epaisseur / new Thickness
!  -----------------------------------------
 
          dzsnew         =                                             &
     &    dzsnSV(ikl,ikv,isn)                                          &
     &     *max(vdz3,                                                  &
     &         (un_1-dt__SV*max(SnMass(ikl,ikv)*cos(slorSV(ikl,ikv)),un_1) &
     &                     /max(ViscSn                      ,eps6)))
          rosnew         = ro__SV(ikl,ikv,isn) *dzsnSV(ikl,ikv,isn)    &
     &                            /max(eps6,dzsnew)
          rosmax         = 1.d0   /( (1.d0 -eta_SV(ikl,ikv,isn)) /rhoIce   &
     &                               +      eta_SV(ikl,ikv,isn)  /rhoWat)
          rosnew         =                        min(rosnew ,rosmax)
          dzsnSV(ikl,ikv,isn)= dzsnSV(ikl,ikv,isn) *ro__SV(ikl,ikv,isn)&
     &                            /max(eps6,rosnew)
          ro__SV(ikl,ikv,isn)= rosnew
          ro_dry(ikl,ikv,isn)= ro__SV(ikl,ikv,isn)*(1.-eta_SV(ikl,ikv,isn))*1.e-3
!         ro_dry: Dry Density (g/cm3)
 
          SnMass(ikl,ikv)    = SnMass(ikl,ikv)+dSnMas*0.5
        END DO
        END DO
      END DO
 
 
! OUTPUT/Verification: Snow   Properties
! #vs DO ikl = 1,kcolp
! #vs DO ikv=1,mwp
! #vs DO isn = 1,isnoSV(ikl,ikv)
! #vs   IF (G1snSV(ikl,ikv,isn).gt.0. .AND. G2snSV(ikl,ikv,isn).gt.D__MAX) THEN
! #vs     write(6,6600) G1snSV(ikl,ikv,isn),G2snSV(ikl,ikv,isn),ikl,ikv,isn
! #vs6600 format(/,'WARNING in _GSn: G1,G2 =',2f9.3,'  (ikl,ikv,isn) =',2i4)
! #vs     D__MAX =                      G2snSV(ikl,ikv,isn)
! #vs   END IF
! #vs   IF (                            G2snSV(ikl,ikv,isn).lt.0.    ) THEN
! #vs     write(6,6601) G1snSV(ikl,ikv,isn),G2snSV(ikl,ikv,isn),ikl,ikv,isn
! #vs6601 format(/,'ERROR 1 in _GSn: G1,G2 =',2f9.3,'  (ikl,ikv,isn) =',2i4)
! #vs     STOP
! #vs   END IF
! #vs   IF (G1snSV(ikl,ikv,isn).gt.G1_dSV+eps6                       ) THEN
! #vs     write(6,6602) G1snSV(ikl,ikv,isn),G2snSV(ikl,ikv,isn),ikl,ikv,isn
! #vs6602 format(/,'ERROR 2 in _GSn: G1,G2 =',2f9.3,'  (ikl,ikv,isn) =',2i4)
! #vs     STOP
! #vs   END IF
! #vs   IF (G1snSV(ikl,ikv,isn).lt.0.                           .AND.  &
! #vs&      G2snSV(ikl,ikv,isn).gt.G1_dSV+eps6                       ) THEN
! #vs     write(6,6603) G1snSV(ikl,ikv,isn),G2snSV(ikl,ikv,isn),ikl,ikv,isn
! #vs6603 format(/,'ERROR 3 in _GSn: G1,G2 =',2f9.3,'  (ikl,ikv,isn) =',2i4)
! #vs     STOP
! #vs   END IF
! #vs END DO
! #vs END DO
! #vs END DO
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
      IF (FlagDALLOC)                                             THEN !
 
      deallocate          ( ro_dry )                                   ! Dry Density            [g/cm3]
      deallocate          ( etaSno )                                   ! Liquid Water Content   [g/cm2]
      deallocate          ( SnMass )                                   ! Snow   Mass            [kg/m2]
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SISVAT_GSn
 
 
 
      subroutine SISVAT_qSo(                                           &
! #m0&                      Wats_0,Wats_1,Wats_d                       &
     &                     )
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_qSo                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_qSo computes the Soil      Water  Balance          |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     PARAMETERS:  kcolv: Total Number of columns =                        |
!     ^^^^^^^^^^        = Total Number of continental     grid boxes       |
!                       X       Number of Mosaic Cell per grid box         |
!                                                                          |
!     INPUT:   isnoSV   = total Nb of Ice/Snow Layers                      |
!     ^^^^^    isotSV   = 0,...,11:   Soil       Type                      |
!                         0:          Water, Solid or Liquid               |
!                                                                          |
!     INPUT:   rhT_SV   : SBL Top    Air  Density                  [kg/m3] |
!     ^^^^^    drr_SV   : Rain   Intensity                       [kg/m2/s] |
!              LSdzsv   : Vertical   Discretization Factor             [-] |
!                       =    1. Soil                                       |
!                       = 1000. Ocean                                      |
!              dt__SV   : Time   Step                                  [s] |
!                                                                          |
!              Lx_H2O   : Latent Heat of Vaporization/Sublimation   [J/kg] |
!              HLs_sv   : Latent Heat  Flux                         [W/m2] |
!              Rootsv   : Root   Water Pump                      [kg/m2/s] |
!                                                                          |
!     INPUT /  eta_SV   : Water      Content                       [m3/m3] |
!     OUTPUT:  Khydsv   : Soil   Hydraulic    Conductivity           [m/s] |
!     ^^^^^^                                                               |
!                                                                          |
!     OUTPUT:  RnofSV   : RunOFF Intensity                       [kg/m2/s] |
!     ^^^^^^   Wats_0   : Soil Water,  before Forcing                 [mm] |
!              Wats_1   : Soil Water,  after  Forcing                 [mm] |
!              Wats_d   : Soil Water          Forcing                 [mm] |
!                                                                          |
!     Internal Variables:                                                  |
!     ^^^^^^^^^^^^^^^^^^                                                   |
!              z_Bump   : (Partly)Bumpy Layers Height                  [m] |
!              z0Bump   :         Bumpy Layers Height                  [m] |
!              dzBump   :  Lowest Bumpy Layer:                         [m] |
!              etBump   :         Bumps Layer Averaged Humidity    [m3/m3] |
!              etaMid   : Layer Interface's Humidity               [m3/m3] |
!              eta__f   : Layer             Humidity  (Water Front)[m3/m3] |
!              Dhyd_f   : Soil  Hydraulic Diffusivity (Water Front) [m2/s] |
!              Dhydif   : Soil  Hydraulic Diffusivity               [m2/s] |
!              WgFlow   : Water         gravitational     Flux   [kg/m2/s] |
!              Wg_MAX   : Water MAXIMUM gravitational     Flux   [kg/m2/s] |
!              SatRat   : Water         Saturation        Flux   [kg/m2/s] |
!              WExces   : Water         Saturation Excess Flux   [kg/m2/s] |
!              Dhydtz   : Dhydif * dt / dz                             [m] |
!              FreeDr   : Free Drainage Fraction                       [-] |
!              Elem_A   : A Diagonal Coefficient                           |
!              Elem_C   : C Diagonal Coefficient                           |
!              Diag_A   : A Diagonal                                       |
!              Diag_B   : B Diagonal                                       |
!              Diag_C   : C Diagonal                                       |
!              Term_D   :   Independant Term                               |
!              Aux__P   : P Auxiliary Variable                             |
!              Aux__Q   : Q Auxiliary Variable                             |
!                                                                          |
!     TUNING PARAMETER:                                                    |
!     ^^^^^^^^^^^^^^^^                                                     |
!              z0soil   : Soil Surface averaged Bumps Height           [m] |
!                                                                          |
!     METHOD: NO   Skin Surface Humidity                                   |
!     ^^^^^^  Semi-Implicit Crank Nicholson Scheme                         |
!             (Partial) free Drainage, Water Bodies excepted (Lakes, Sea)  |
!                                                                          |
!                                                                          |
!     Preprocessing  Option:                                               |
!     ^^^^^^^^^^^^^^^^^^^^^                                                |
!     #GF: Saturation Front                                                |
!     #GH: Saturation Front allows Horton Runoff                           |
!     #GA: Soil Humidity Geometric Average                                 |
!     #TB: Parameterization of Terrain Bumps                               |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_iii_jjj_n     | #m0: OUTPUT/Verification: H2O    Conservation |
!   # SISVAT_iii_jjj_n     | #m1: OUTPUT/Verification: * Mass Conservation |
!   # stdout               | #mw: OUTPUT/Verification: H2O    Conservation |
!                          |      unit  6, SubRoutine  SISVAT_qSo **ONLY** |
!   # SISVAT_qSo.vw        | #vw: OUTPUT/Verif+Detail: H2O    Conservation |
!                          |      unit 42, SubRoutine  SISVAT_qSo **ONLY** |
!   # stdout               | #vg: OUTPUT/Verification: Gravitational Front |
!                          |      unit  6, SubRoutine  SISVAT_qSo **ONLY** |
!                                                                          |
!     REMARQUE: Inclure possibilite de creer mare sur bedrock impermeable  |
!     ^^^^^^^^                                                             |
!--------------------------------------------------------------------------+
 
 
! Global  Variables
! =================
 
      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_grd
 
 
 
! General Variables
! =================
 
      use Mod_SISVAT_ctr
      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
 
 
 
! Internal Variables
! ==================
 
      use Mod_SISVATLqSo
 
 
      IMPLICIT NONE
 
 
      integer                                        ::  isl           !
      integer                                        ::  ist   ,ikl,ikv!
      integer                                        ::  ikm   ,ikp    !
      integer                                        ::  ik0           !
      integer                                        ::  ist__s,ist__w ! Soil/Water Body Identifier
 
 
! #TB real(kind=real8)                               ::  z0soil = 0.02 ! Soil Surface Bumps Height  [m]
! #TB real(kind=real8)                               ::  z_Bump        !(Partly)Bumpy Layers Height [m]
! #TB real(kind=real8)                               ::  z0Bump        !        Bumpy Layers Height [m]
! #TB real(kind=real8)                               ::  dzBump        ! Lowest Bumpy Layer:
 
      real(kind=real8)                               ::  etaMid        ! Layer Interface's Humidity
      real(kind=real8)                               ::  Dhydif        ! Hydraulic Diffusivity   [m2/s]
! #GH real(kind=real8)                               ::  eta__f        ! Water Front Soil Water Content
      real(kind=real8)                               ::  Khyd_f        ! Water Front Hydraulic Conduct.
      real(kind=real8)                               ::  Khydav        ! Hydraulic Conductivity   [m/s]
! #GF real(kind=real8)                               ::  WgFlow        ! Water gravitat. Flux [kg/m2/s]
      real(kind=real8)                               ::  Wg_MAX        ! Water MAX.grav. Flux [kg/m2/s]
      real(kind=real8)                               ::  SatRat        ! Saturation      Flux [kg/m2/s]
! #GF real(kind=real8)                               ::  WExces        ! Saturat. Excess Flux [kg/m2/s]
      real(kind=real8)                               ::  Elem_A        !   Diagonal Coefficients
      real(kind=real8)                               ::  Elem_B        !   Diagonal Coefficients
      real(kind=real8)                               ::  Elem_C        !   Diagonal Coefficients
      real(kind=real8)                               ::  FreeDr        ! Free Drainage Fraction (actual)
!     real(kind=real8)                               ::  FreeD0 = 1.00 ! Free Drainage Fraction (1=Full)
 
! OUTPUT
! ------
 
! OUTPUT/Verification: H2O    Conservation
! #mw real(kind=real8)                               ::  hourwr
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! ALLOCATION                                                           !
! ==========                                                           !
 
      IF (it_RUN.EQ.1 .OR. FlagDALLOC)                            THEN !
 
! OUTPUT/Verification: H2O    Conservation
! #m0 allocate            ( Wats_0(kcolp,mwp) )                        ! Soil Water,  before forcing
! #m0 allocate            ( Wats_1(kcolp,mwp) )                        ! Soil Water,  after  forcing
! #m0 allocate            ( Wats_d(kcolp,mwp) )                        ! Soil Water          forcing
 
! #TB allocate            ( etBump(kcolp,mwp) )                        ! Bumps Layer Averaged Humidity
 
      allocate            ( SoRnOF(kcolp,mwp) )                        ! Soil     Run    OFF
      allocate            ( Dhydtz(kcolp,mwp,-nsoil:0) )               ! Dhydif * dt / dz           [m]
      allocate            ( Diag_A(kcolp,mwp,-nsoil:0) )               ! A Diagonal
      allocate            ( Diag_B(kcolp,mwp,-nsoil:0) )               ! B Diagonal
      allocate            ( Diag_C(kcolp,mwp,-nsoil:0) )               ! C Diagonal
      allocate            ( Term_D(kcolp,mwp,-nsoil:0) )               !   Independant Term
      allocate            ( Aux__P(kcolp,mwp,-nsoil:0) )               ! P Auxiliary Variable
      allocate            ( Aux__Q(kcolp,mwp,-nsoil:0) )               ! Q Auxiliary Variable
      allocate            ( etaaux(kcolp,mwp,-nsoil:-nsoil+1) )        ! Soil Water Content     [m3/m3]
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
! OUTPUT/Verification: H2O    Conservation: Water  Budget (IN)
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Wats_0(ikl,ikv) = 0.                ! OLD RunOFF Contrib.
! #m0     Wats_d(ikl,ikv) = drr_SV(ikl,ikv)   ! Water Surface Forc.
! #m0   END DO
! #m0   END DO
 
! #m0      isl= -nsoil
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Wats_0(ikl,ikv) = Wats_0(ikl,ikv)                            &
! #m0&      + rhoWat *( eta_SV(ikl,ikv,isl)   *dz78SV(isl)             &
! #m0&                + eta_SV(ikl,ikv,isl+1) *dz_8SV(isl) ) * LSdzsv(ikl,ikv)
! #m0   END DO
! #m0   END DO
 
! #m0 DO   isl= -nsoil+1,-1
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Wats_0(ikl,ikv) = Wats_0(ikl,ikv)                            &
! #m0&      + rhoWat *( eta_SV(ikl,ikv,isl)   *dz34SV(isl)             &
! #m0&                +(eta_SV(ikl,ikv,isl-1)                          &
! #m0&                 +eta_SV(ikl,ikv,isl+1))*dz_8SV(isl) ) * LSdzsv(ikl,ikv)
! #m0   END DO
! #m0   END DO
! #m0 END DO
 
! #m0      isl=  0
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Wats_0(ikl,ikv) = Wats_0(ikl,ikv)                            &
! #m0&      + rhoWat *( eta_SV(ikl,ikv,isl)   *dz78SV(isl)             &
! #m0&                + eta_SV(ikl,ikv,isl-1) *dz_8SV(isl) ) * LSdzsv(ikl,ikv)
! #m0   END DO
! #m0   END DO
 
 
! Gravitational Flow
! ==================
 
! .    METHOD: Surface Water Flux saturates successively the soil layers
!      ^^^^^^  from up to below, but is limited by infiltration capacity.
!              Hydraulic Conductivity again contributes after this step,
!              not redundantly because of a constant (saturated) profile.
 
! Flux  Limitor
! ^^^^^^^^^^^^^
           isl=0
        DO ikl=1,kcolp
        DO ikv=1,mwp
          ist    = isotSV(ikl,ikv)                 ! Soil Type
          ist__s = min(ist, 1)                     ! 1 => Soil
          ist__w = 1 - ist__s                      ! 1 => Water Body
          Dhydif = s1__SV(ist)                    &!
     &               *max(eps6,eta_SV(ikl,ikv,isl))   &! Hydraulic Diffusivity
     &                      **(bCHdSV(ist)+2.)     ! DR97, Eqn.(3.36)
          Dhydif = ist__s    * Dhydif             &!
     &           + ist__w    * vK_dSV              ! Water Bodies
 
          Khydav = ist__s    * Ks_dSV(ist)        &! DR97  Assumption
     &           + ist__w    * vK_dSV              ! Water Bodies
 
          Wg_MAX = rhoWat     *Dhydif             &! MAXimum  Infiltration
     &           *(etadSV(ist)-eta_SV(ikl,ikv,isl))   &!          Rate
     &           /(dzAvSV(isl)*LSdzsv(ikl,ikv)    )   &!
     &          +  rhoWat     *Khydav              !
 
! Surface Horton RunOFF
! ^^^^^^^^^^^^^^^^^^^^^
          SoRnOF(ikl,ikv) =                       &!
     &                max(zer0,drr_SV(ikl,ikv)-Wg_MAX) !
          drr_SV(ikl,ikv) =        drr_SV(ikl,ikv)-SoRnOF(ikl,ikv)
        END DO
        END DO
 
! #GF DO   isl=0,-nsoil,-1
! #GF   DO ikl=1,kcolp
! #GF   DO ikv=1,mwp
! #GF     ist    = isotSV(ikl,ikv)                 ! Soil Type
! #GF     ist__s = min(ist, 1)                     ! 1 => Soil
! #GF     ist__w = 1 - ist__s                      ! 1 => Water Body
 
! Water Diffusion
! ^^^^^^^^^^^^^^^
! #GF     Dhydif = s1__SV(ist)                    &!
! #GF&               *max(eps6,eta_SV(ikl,ikv,isl))   &! Hydraulic Diffusivity
! #GF&                      **(bCHdSV(ist)+2.)     ! DR97, Eqn.(3.36)
! #GF     Dhydif = ist__s    * Dhydif             &!
! #GF&           + ist__w    * vK_dSV              ! Water Bodies
 
! Water Conduction (without Horton Runoff)
! ^^^^^^^^^^^^^^^^
! #GF     Khyd_f =             Ks_dSV(ist)
!         Uses saturated K ==> Horton Runoff ~0    !
 
! Water Conduction (with    Horton Runoff)
! ^^^^^^^^^^^^^^^^
! #GH     ik0    = nkhy   *int(eta_SV(ikl,ikv,isl)&!
! #GH&                        /etadSV(ist))        !
! #GH     eta__f         =            1.          &!
! #GH&   -aKdtSV(ist,ik0)/(2. *dzAvSV(isl)        &!
! #GH&                        *LSdzsv(ikl,ikv))    !
! #GH     eta__f         = max(eps_21,eta__f)
! #GH     eta__f         = min(etadSV(ist),       &!
! #GH&                         eta_SV(ikl,ikv,isl) +  &!
! #GH&   (aKdtSV(ist,ik0)     *eta_SV(ikl,ikv,isl)&!
! #GH&   +bKdtSV(ist,ik0))   /(dzAvSV(isl)        &!
! #GH&                        *LSdzsv(ikl,ikv))   &!
! #GH&                       / eta__f          )
 
! #GH     eta__f         = .5*(eta_SV(ikl,ikv,isl)&!
! #GH&                        +eta__f)             !
!         eta__f         =     eta_SV(ikl,isl)     ! Another Possibility
 
! #GH     ik0    = nkhy       *eta__f             &!
! #GH&                        /etadSV(ist)         !
! #GH     Khyd_f =                                &!
! #GH&   (aKdtSV(ist,ik0)     *eta__f             &!
! #GH&   +bKdtSV(ist,ik0))    /dt__SV              !
 
! #GF     Khydav = ist__s    * Khyd_f             &! DR97  Assumption
! #GF&           + ist__w    * vK_dSV              ! Water Bodies
 
! Gravitational Flow
! ^^^^^^^^^^^^^^^^^^
! #GF     Wg_MAX =                                &! MAXimum  Infiltration
! #GF&             rhoWat     *Dhydif             &!          Rate
! #GF&           *(etadSV(ist)-eta_SV(ikl,ikv,isl))   &!
! #GF&           /(dzAvSV(isl)*LSdzsv(ikl,ikv)    )   &!
! #GF&          +  rhoWat     *Khydav              !
 
! OUTPUT/Verification: Gravitational Front
! #vg     write(6,6001) isl,drr_SV(ikl,ikv)*3.6e3,Wg_MAX     *3.6e3
! #vg6001 format(i3,'  vRain ,Wg_MAX ',2e12.3,' mm/hr')
 
! #GF     WgFlow =  min(Wg_MAX,drr_SV(ikl,ikv))    ! Infiltration Rate
! #GF     WExces =  max(zer0  ,drr_SV(ikl,ikv)-WgFlow) ! Water Excess => RunOff
 
! OUTPUT/Verification: Gravitational Front
! #vg     write(6,6002)     WgFlow     *3.6e3,WExces     *3.6e3
! #vg6002 format(3x,'  WgFlow,WExces ',2e12.3,' mm/hr')
 
! #GF     SoRnOF(ikl,ikv) =        SoRnOF(ikl,ikv)+WExces  !
! #GF     drr_SV(ikl,ikv) =        WgFlow          !
 
! OUTPUT/Verification: Gravitational Front
! #vg     write(6,6003)     SoRnOF(ikl,ikv)*3.6e3,drr_SV(ikl,ikv)*3.6e3
! #vg6003 format(3x,'  SoRnOF,drr_SV ',2e12.3,' mm/hr')
 
! #GF     SatRat =(etadSV(ist)-eta_SV(ikl,ikv,isl))   &! Saturation   Rate
! #GF&            *rhoWat     *dzAvSV(isl)        &!
! #GF&                        *LSdzsv(ikl,ikv)/dt__SV  !
! #GF     SatRat =  min(SatRat,drr_SV(ikl,ikv))    !
! #GF     drr_SV(ikl,ikv) =        drr_SV(ikl,ikv)-SatRat  ! Water Flux for Below
 
! OUTPUT/Verification: Gravitational Front
! #vg     write(6,6004)     SatRat     *3.6e3,drr_SV(ikl,ikv)*3.6e3
! #vg6004 format(3x,'  SatRat,drr_SV ',2e12.3,' mm/hr')
! #vg     write(6,6005)     eta_SV(ikl,ikv,isl)*1.e3
 
! #GF     eta_SV(ikl,ikv,isl) =    eta_SV(ikl,ikv,isl)&! Saturation
! #GF&                        +SatRat*dt__SV      &!
! #GF&                       /(rhoWat*dzAvSV(isl) &!
! #GF&                               *LSdzsv(ikl,ikv)) !
 
! OUTPUT/Verification: Gravitational Front
! #vg     write(6,6005)     eta_SV(ikl,ikv,isl)*1.e3
! #vg6005 format(3x,'  eta_SV,       ',e12.3,' g/kg')
! #GF   END DO
! #GF   END DO
! #GF END DO
! #GF   DO ikl=1,kcolp
! #GF   DO ikv=1,mwp
! #GF     SoRnOF(ikl,ikv)     =    SoRnOF(ikl,ikv)&! RunOFF Intensity
! #GF&                    +    drr_SV(ikl,ikv)     ! [kg/m2/s]
!         Inclure la possibilite de creer une mare sur un bedrock impermeable
! #GF     drr_SV(ikl,ikv) = 0.
! #GF   END DO
! #GF   END DO
 
 
! Temperature Correction due to a changed Soil Energy Content
! ===========================================================
 
! REMARQUE: Mettre en oeuvre le couplage humidite-energie
! ^^^^^^^^
 
 
! Full Resolution of the Richard's Equation
! =========================================
 
!         METHOD: Water content evolution results from water fluxes
!         ^^^^^^  at the layer boundaries
!                 Conductivity is approximated by a piecewise linear profile.
!                 Semi-Implicit Crank-Nicholson scheme is used.
!                (Bruen, 1997, Sensitivity of hydrological processes
!                              at the land-atmosphere interface.
!                              Proc. Royal Irish Academy,  IGBP symposium
!                              on global change and the Irish Environment.
!                              Publ.: Maynooth)
 
!                        - - - - - - - -   isl+1/2   - -  ^
!                                                         |
!     eta_SV(isl)        ---------------   isl     -----  +--dz_dSV(isl)  ^
!                                                         |               |
!     Dhydtz(isl) etaMid - - - - - - - -   isl-1/2   - -  v  dzmiSV(isl)--+
!                                                                         |
!     eta_SV(isl-1)      ---------------   isl-1   -----                  v
 
! Transfert       Coefficients
! ----------------------------
 
      DO   isl=-nsoil+1,0
        DO ikl=1,kcolp
        DO ikv=1,mwp
          ist    =      isotSV(ikl,ikv)                   ! Soil Type
          ist__s =      min(ist, 1)                       ! 1 => Soil
          ist__w =      1 - ist__s                        ! 1 => Water Body
          etaMid =     (dz_dSV(isl)  *eta_SV(ikl,ikv,isl-1)  &! eta at layers
     &                 +dz_dSV(isl-1)*eta_SV(ikl,ikv,isl)  ) &!     interface
     &           /(2.0* dzmiSV(isl))                      ! LSdzsv implicit !
! #GA     etaMid = sqrt(dz_dSV(isl)  *eta_SV(ikl,ikv,isl-1)  &! Idem, geometric
! #GA&                 *dz_dSV(isl-1)*eta_SV(ikl,ikv,isl)  ) &!       average
! #GA&           /(2.0* dzmiSV(isl))                      ! (Vauclin&al.1979)
          Dhydif          =    s1__SV(ist)               &! Hydraul.Diffusi.
     &  *(etaMid         **(   bCHdSV(ist)+2.))           ! DR97, Eqn.(3.36)
          Dhydtz(ikl,ikv,isl) =    Dhydif*dt__SV         &!
     &                              /(dzmiSV(isl)        &!
     &                               *LSdzsv(ikl,ikv))    !
          Dhydtz(ikl,ikv,isl) =    Dhydtz(ikl,ikv,isl) * ist__s  &! Soil
     &        +0.5*dzmiSV(isl)*LSdzsv(ikl,ikv)     * ist__w   ! Water bodies
 
        END DO
        END DO
      END DO
           isl=-nsoil
        DO ikl=1,kcolp
        DO ikv=1,mwp
          Dhydtz(ikl,ikv,isl) =    0.0                    !
        END DO
        END DO
 
 
! Tridiagonal Elimination: Set Up
! -------------------------------
 
! Soil/Snow Interior
! ^^^^^^^^^^^^^^^^^^
      DO   isl=-nsoil,-nsoil+1
        DO ikl=1,kcolp
        DO ikv=1,mwp
          etaaux(ikl,ikv,isl) =  eta_SV(ikl,ikv,isl)
        END DO
        END DO
      END DO
 
      DO   isl=-nsoil+1,-1
        DO ikl=1,kcolp
        DO ikv=1,mwp
          ist      =         isotSV(ikl,ikv)
          ikm      =nkhy*int(eta_SV(ikl,ikv,isl-1) / etadSV(ist))
          ik0      =nkhy*int(eta_SV(ikl,ikv,isl)   / etadSV(ist))
          ikp      =nkhy*int(eta_SV(ikl,ikv,isl+1) / etadSV(ist))
          Elem_A   =         Dhydtz(ikl,ikv,isl)                       &
     &                    -  aKdtSV(ist,ikm)* dziiSV(isl)  *LSdzsv(ikl,ikv)
          Elem_B   =      - (Dhydtz(ikl,ikv,isl)                       &
     &                      +Dhydtz(ikl,ikv,isl+1)                     &
     &                      -aKdtSV(ist,ik0)*(dziiSV(isl+1)            &
     &                                       -dzi_SV(isl) )*LSdzsv(ikl,ikv))
          Elem_C   =         Dhydtz(ikl,ikv,isl+1)                     &
     &                    +  aKdtSV(ist,ikp)* dzi_SV(isl+1)*LSdzsv(ikl,ikv)
          Diag_A(ikl,ikv,isl) =  dz_8SV(isl)        *LSdzsv(ikl,ikv)   &
     &                      -Implic            * Elem_A
          Diag_B(ikl,ikv,isl) =  dz34SV(isl)        *LSdzsv(ikl,ikv)   &
     &                      -Implic            * Elem_B
          Diag_C(ikl,ikv,isl) =  dz_8SV(isl)        *LSdzsv(ikl,ikv)   &
     &                      -Implic            * Elem_C
 
          Term_D(ikl,ikv,isl) = (dz_8SV(isl) *LSdzsv(ikl,ikv)          &
     &                      +Explic      *Elem_A     )*eta_SV(ikl,ikv,isl-1)&
     &                    + (dz34SV(isl) *LSdzsv(ikl,ikv)              &
     &                      +Explic      *Elem_B     )*eta_SV(ikl,ikv,isl) &
     &                    + (dz_8SV(isl) *LSdzsv(ikl,ikv)              &
     &                      +Explic      *Elem_C     )*eta_SV(ikl,ikv,isl+1)&
     &                    + (bKdtSV(ist,ikp)* dzi_SV(isl+1)            &
     &                      +bKdtSV(ist,ik0)*(dziiSV(isl+1)            &
     &                                       -dzi_SV(isl)  )           &
     &                      -bKdtSV(ist,ikm)* dziiSV(isl)   )          &
     &                                      * LSdzsv(ikl,ikv)          &
     &                    -  dt__SV         * Rootsv(ikl,ikv,isl)/rhoWat
        END DO
        END DO
      END DO
 
           isl=-nsoil
        DO ikl=1,kcolp
        DO ikv=1,mwp
          ist      =         isotSV(ikl,ikv)
          FreeDr   =         iWaFSV(ikl,ikv)       *  min(ist,1)
!         FreeDr   =         FreeD0            *  min(ist,1) ! Free Drainage
          ik0      =nkhy*int(eta_SV(ikl,ikv,isl  ) / etadSV(ist))
          ikp      =nkhy*int(eta_SV(ikl,ikv,isl+1) / etadSV(ist))
          Elem_A   =         0.
          Elem_B   =      - (Dhydtz(ikl,ikv,isl+1)                     &
     &                      -aKdtSV(ist,ik0)*(dziiSV(isl+1)*LSdzsv(ikl,ikv)&
     &                                       -FreeDr                  ))
          Elem_C   =         Dhydtz(ikl,ikv,isl+1)                     &
     &                    +  aKdtSV(ist,ikp)* dzi_SV(isl+1)*LSdzsv(ikl,ikv)
          Diag_A(ikl,ikv,isl) =  0.
          Diag_B(ikl,ikv,isl) =  dz78SV(isl) *LSdzsv(ikl,ikv)          &
     &                      -Implic      *Elem_B
          Diag_C(ikl,ikv,isl) =  dz_8SV(isl) *LSdzsv(ikl,ikv)          &
     &                      -Implic      *Elem_C
 
          Term_D(ikl,ikv,isl) = (dz78SV(isl) *LSdzsv(ikl,ikv)          &
     &                      +Explic      *Elem_B     )*eta_SV(ikl,ikv,isl) &
     &                    + (dz_8SV(isl) *LSdzsv(ikl,ikv)              &
     &                      +Explic      *Elem_C     )*eta_SV(ikl,ikv,isl+1)&
     &                    + (bKdtSV(ist,ikp)* dzi_SV(isl+1)*LSdzsv(ikl,ikv)&
     &                      +bKdtSV(ist,ik0)*(dziiSV(isl+1)*LSdzsv(ikl,ikv)&
     &                                       -FreeDr                  ))&
     &                    -  dt__SV         * Rootsv(ikl,ikv,isl)/rhoWat
        END DO
        END DO
 
           isl=0
        DO ikl=1,kcolp
        DO ikv=1,mwp
          ist      =         isotSV(ikl,ikv)
          ikm      =nkhy*int(eta_SV(ikl,ikv,isl-1) / etadSV(ist))
          ik0      =nkhy*int(eta_SV(ikl,ikv,isl)   / etadSV(ist))
          Elem_A   =         Dhydtz(ikl,ikv,isl)                       &
     &                    -  aKdtSV(ist,ikm)* dziiSV(isl)*LSdzsv(ikl,ikv)
          Elem_B   =      - (Dhydtz(ikl,ikv,isl)                       &
     &                      +aKdtSV(ist,ik0)* dzi_SV(isl)*LSdzsv(ikl,ikv))
          Elem_C   =         0.
          Diag_A(ikl,ikv,isl) =  dz_8SV(isl) *LSdzsv(ikl,ikv)          &
     &                    -  Implic      *Elem_A
          Diag_B(ikl,ikv,isl) =  dz78SV(isl) *LSdzsv(ikl,ikv)          &
     &                    -  Implic      *Elem_B
          Diag_C(ikl,ikv,isl) =  0.
 
          Term_D(ikl,ikv,isl) = (dz_8SV(isl) *LSdzsv(ikl,ikv)          &
     &                      +Explic      *Elem_A     )*eta_SV(ikl,ikv,isl-1)&
     &                    + (dz78SV(isl) *LSdzsv(ikl,ikv)              &
     &                      +Explic      *Elem_B     )*eta_SV(ikl,ikv,isl) &
     &                    - (bKdtSV(ist,ik0)* dzi_SV(isl)              &
     &                      +bKdtSV(ist,ikm)* dziiSV(isl))*LSdzsv(ikl,ikv) &
     &            + dt__SV *(HLs_sv(ikl,ikv)    *     (1-min(1,isnoSV(ikl,ikv)))&
     &                                                   / Lx_H2O(ikl,ikv) &
     &                      +drr_SV(ikl,ikv)                           &
     &                      -Rootsv(ikl,ikv,isl)             )/rhoWat
        END DO
        END DO
 
 
! Tridiagonal Elimination
! =======================
 
! Forward  Sweep
! ^^^^^^^^^^^^^^
        DO ikl=  1,kcolp
        DO ikv=1,mwp
          Aux__P(ikl,ikv,-nsoil) = Diag_B(ikl,ikv,-nsoil)
          Aux__Q(ikl,ikv,-nsoil) =-Diag_C(ikl,ikv,-nsoil)/Aux__P(ikl,ikv,-nsoil)
        END DO
        END DO
 
      DO   isl=-nsoil+1,0
        DO ikl=      1,kcolp
        DO ikv=1,mwp
          Aux__P(ikl,ikv,isl)   = Diag_A(ikl,ikv,isl)  *Aux__Q(ikl,ikv,isl-1)  &
     &                       +Diag_B(ikl,ikv,isl)
          Aux__Q(ikl,ikv,isl)   =-Diag_C(ikl,ikv,isl)  /Aux__P(ikl,ikv,isl)
        END DO
        END DO
      END DO
 
        DO ikl=      1,kcolp
        DO ikv=1,mwp
          eta_SV(ikl,ikv,-nsoil) = Term_D(ikl,ikv,-nsoil)/Aux__P(ikl,ikv,-nsoil)
        END DO
        END DO
 
      DO   isl=-nsoil+1,0
        DO ikl=      1,kcolp
        DO ikv=1,mwp
          eta_SV(ikl,ikv,isl)   =(Term_D(ikl,ikv,isl)                  &
     &                       -Diag_A(ikl,ikv,isl)  *eta_SV(ikl,ikv,isl-1)) &
     &                                         /Aux__P(ikl,ikv,isl)
        END DO
        END DO
      END DO
 
! Backward Sweep
! ^^^^^^^^^^^^^^
      DO   isl=-1,-nsoil,-1
        DO ikl= 1,kcolp
        DO ikv=1,mwp
          eta_SV(ikl,ikv,isl)   = Aux__Q(ikl,ikv,isl)  *eta_SV(ikl,ikv,isl+1)  &
     &                                         +eta_SV(ikl,ikv,isl)
        END DO
        END DO
      END DO
 
 
! Horton RunOFF Intensity
! =======================
 
      DO   isl=0,-nsoil,-1
        DO ikl=1,kcolp
        DO ikv=1,mwp
          ist    =   isotSV(ikl,ikv)               ! Soil Type
          SatRat =  (eta_SV(ikl,ikv,isl)-etadSV(ist)) &! OverSaturation Rate
     &              *rhoWat         *dzAvSV(isl)  &!
     &                              *LSdzsv(ikl,ikv)  &!
     &                              /dt__SV        !
          SoRnOF(ikl,ikv)     =          SoRnOF(ikl,ikv)  &!
     &                    + max(zer0,SatRat)       !
          eta_SV(ikl,ikv,isl) = max(eps6          &!
     &                         ,eta_SV(ikl,ikv,isl))   !
          eta_SV(ikl,ikv,isl) = min(eta_SV(ikl,ikv,isl)   &!
     &                         ,etadSV(ist)    )   !
        END DO
        END DO
      END DO
 
! OUTPUT/Verification: Soil Vertic.Discret.
! #kw     write(6,6010)
! #kw6010 format(/,1x)
      DO   isl= 0,-nsoil,-1
        DO ikl= 1,kcolp
        DO ikv=1,mwp
          ist      =          isotSV(ikl,ikv)
          ikp      = nkhy*int(eta_SV(ikl,ikv,isl)  /etadSV(ist))
          Khydsv(ikl,ikv,isl)   =(aKdtSV(ist,ikp)  *eta_SV(ikl,ikv,isl)&
     &                           +bKdtSV(ist,ikp)) *2.0/dt__SV
! OUTPUT/Verification: Soil Vertic.Discret.
! #kw     write(6,6011) ikl,ikv,isl,eta_SV(ikl,ikv,isl)*1.e3,          &
! #kw&                  ikp,    aKdtSV(ist,ikp),bKdtSV(ist,ikp),       &
! #kw&                          Khydsv(ikl,ikv,isl)
! #kw6011 format(2i3,f8.1,i3,3e12.3)
 
        END DO
        END DO
      END DO
 
 
! Additional RunOFF Intensity
! ===========================
 
        DO ikl=1,kcolp
        DO ikv=1,mwp
          ist      =          isotSV(ikl,ikv)
          ik0      = nkhy*int(etaaux(ikl,ikv,-nsoil  ) /etadSV(ist))
          FreeDr   =          iWaFSV(ikl,ikv)       *  min(ist,1)
!         FreeDr   =          FreeD0            *  min(ist,1) ! Free Drainage
          SoRnOF(ikl,ikv) =  SoRnOF(ikl,ikv)                           &
     &                + (aKdtSV(ist,ik0)*(etaaux(ikl,ikv,-nsoil)*Explic &
     &                                   +eta_SV(ikl,ikv,-nsoil)*Implic)&
     &                 + bKdtSV(ist,ik0)                           )   &
     &                 * FreeDr          *rhoWat           /dt__SV
 
! Full Run OFF: Update
! ~~~~~~~~~~~~~~~~~~~~
          RnofSV(ikl,ikv) = RnofSV(ikl,ikv)     + SoRnOF(ikl,ikv)
        END DO
        END DO
 
 
! Temperature Correction due to a changed Soil Energy Content
! ===========================================================
 
! REMARQUE: Mettre en oeuvre le couplage humidite-energie
! ^^^^^^^^
 
 
! Bumps/Asperites Treatment
! =========================
 
! Average over Bump Depth (z0soil)
! --------------------------------
 
! #TB       z_Bump      = 0.
! #TB     DO ikl=1,kcolp
! #TB     DO ikv=1,mwp
! #TB       etBump(ikl,ikv) = 0.
! #TB     END DO
! #TB     END DO
 
! #TB DO     isl=0,-nsoil,-1
! #TB       z0Bump      = z_Bump
! #TB       z_Bump      = z_Bump      +  dzAvSV(isl)
! #TB   IF (z_Bump.lt.z0soil)                                       THEN
! #TB     DO ikl=1,kcolp
! #TB     DO ikv=1,mwp
! #TB       etBump(ikl,ikv) = etBump(ikl,ikv) +  dzAvSV(isl)   *eta_SV(ikl,ikv,isl)
! #TB     END DO
! #TB     END DO
! #TB   END IF
! #TB   IF (z_Bump.gt.z0soil.AND.z0Bump.lt.z0soil)                  THEN
! #TB     DO ikl=1,kcolp
! #TB     DO ikv=1,mwp
! #TB       etBump(ikl,ikv) = etBump(ikl,ikv) + (z0soil-z0Bump)*eta_SV(ikl,ikv,isl)
! #TB       etBump(ikl,ikv) = etBump(ikl,ikv) /  z0soil
! #TB     END DO
! #TB     END DO
! #TB   END IF
! #TB END DO
 
 
! Correction
! ----------
 
! #TB       z_Bump      = 0.
! #TB DO     isl=0,-nsoil,-1
! #TB       z0Bump =  z_Bump
! #TB       z_Bump =  z_Bump +dzAvSV(isl)
! #TB   IF (z_Bump.lt.z0soil)                                       THEN
! #TB     DO ikl=1,kcolp
! #TB     DO ikv=1,mwp
! #TB       eta_SV(ikl,ikv,isl) = etBump(ikl,ikv)
! #TB     END DO
! #TB     END DO
! #TB   END IF
! #TB   IF (z_Bump.gt.z0soil.AND.z0Bump.lt.z0soil)                  THEN
! #TB       dzBump          =    z_Bump -  z0soil
! #TB     DO ikl=1,kcolp
! #TB     DO ikv=1,mwp
! #TB       eta_SV(ikl,ikv,isl) =(etBump(ikl,ikv)    *(dzAvSV(isl)-dzBump) &
! #TB&                      + eta_SV(ikl,ikv,isl)*             dzBump) &
! #TB&                      /                  dzAvSV(isl)
! #TB     END DO
! #TB     END DO
! #TB   END IF
! #TB END DO
 
 
! Positive Definite
! =================
 
! #TB DO   isl= 0,-nsoil,-1
! #TB   DO ikl= 1,kcolp
! #TB   DO ikv=1,mwp
! #TB     eta_SV(ikl,ikv,isl)   =          max(eps6,eta_SV(ikl,ikv,isl))
! #TB   END DO
! #TB   END DO
! #TB END DO
 
 
! OUTPUT/Verification: H2O    Conservation: Water  Budget (OUT)
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Wats_d(ikl,ikv) = Wats_d(ikl,ikv)           &!
! #m0&                + drr_SV(ikl,ikv)     *0.00     &! Precipitation is
!                                        \______________ already included
! #m0&                + HLs_sv(ikl,ikv)               &!
! #m0&          *(1-min(isnoSV(ikl,ikv),1)) /Lx_H2O(ikl,ikv)  &! Evaporation
! #m0&                - SoRnOF(ikl,ikv)                ! Soil RunOFF Contrib.
! #m0     Wats_1(ikl,ikv) = 0.                         !
 
! OUTPUT/Verification: H2O    Conservation
! #mw     Evapor(ikl,ikv) = HLs_sv(ikl,ikv)     *dt__SV   &!
! #mw&          *(1-min(isnoSV(ikl,ikv),1)) /Lx_H2O(ikl,ikv)   !
 
! #m0   END DO
! #m0   END DO
 
! #m0 DO   isl= -nsoil,0
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Wats_d(ikl,ikv) = Wats_d(ikl,ikv)           &!
! #m0&                - Rootsv(ikl,ikv,isl)            ! Root Extract.
! #m0   END DO
! #m0   END DO
! #m0 END DO
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Wats_d(ikl,ikv) = Wats_d(ikl,ikv)     *dt__SV!
! #m0   END DO
! #m0   END DO
 
! #m0      isl= -nsoil
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Wats_1(ikl,ikv) = Wats_1(ikl,ikv)                            &
! #m0&      + rhoWat *( eta_SV(ikl,ikv,isl)   *dz78SV(isl)             &
! #m0&                + eta_SV(ikl,ikv,isl+1) *dz_8SV(isl) ) *LSdzsv(ikl,ikv)
! #m0   END DO
! #m0   END DO
 
! #m0 DO   isl= -nsoil+1,-1
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Wats_1(ikl,ikv) = Wats_1(ikl,ikv)                            &
! #m0&      + rhoWat *( eta_SV(ikl,ikv,isl)   *dz34SV(isl)             &
! #m0&                +(eta_SV(ikl,ikv,isl-1)                          &
! #m0&                 +eta_SV(ikl,ikv,isl+1))*dz_8SV(isl) ) *LSdzsv(ikl,ikv)
! #m0   END DO
! #m0   END DO
! #m0 END DO
 
! #m0      isl=  0
! #m0   DO ikl=1,kcolp
! #m0   DO ikv=1,mwp
! #m0     Wats_1(ikl,ikv) = Wats_1(ikl,ikv)                            &
! #m0&      + rhoWat *( eta_SV(ikl,ikv,isl)   *dz78SV(isl)             &
! #m0&                + eta_SV(ikl,ikv,isl-1) *dz_8SV(isl) ) *LSdzsv(ikl,ikv)
! #m0   END DO
! #m0   END DO
 
 
! OUTPUT/Verification: H2O    Conservation
! #mw IF (.NOT.mwopen)                                              THEN
! #mw          mwopen = .true.
! #mw          open(unit=42,status='unknown',file='SISVAT_qSo.vw')
! #mw          rewind 42
! #mw   write(42,42)
! #mw42 format('SubRoutine SISVAT_qSo: Local Water Budget',            &
! #mw&       /,'=========================================')
! #mw END IF
! #mw         timewr=timewr + dt__SV
! #mw         hourwr=3600.0
! #mw IF (mod(timewr,hourwr).lt.eps6)                                  &
! #mw&  write(42,420)timewr/hourwr
! #mw420 format(11('-'),'----------+--------------+-',                 &
! #mw&          3('-'),'----------+--------------+-',                  &
! #mw&                 '----------------+----------------+',           &
! #mw&       /,f8.2,3x,'Wats_0(1) |    Wats_d(1) | ',                  &
! #mw&              3x,'Wats_1(1) | W_0+W_d-W_1  | ',                  &
! #mw&                 '   Soil Run OFF |   Soil Evapor. |',           &
! #mw&       /,11('-'),'----------+--------------+-',                  &
! #mw&          3('-'),'----------+--------------+-',                  &
! #mw&                 '----------------+----------------+')
! #mw   write(42,421)   Wats_0(1),Wats_d(1)                            &
! #mw&                 ,Wats_1(1)                                      &
! #mw&                 ,Wats_0(1)+Wats_d(1)-Wats_1(1)                  &
! #mw&                 ,SoRnOF(1),Evapor(1)
! #mw421 format(8x,f12.6,' + ',f12.6,' - ',f12.6,' = ',f12.6,' | ',f12.6,&
! #mw&      '      ',f15.6)
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! DE-ALLOCATION                                                        !
! =============                                                        !
 
      IF (FlagDALLOC)                                             THEN !
 
! OUTPUT/Verification: H2O    Conservation
! #m0 deallocate          ( Wats_0 )                                   ! Soil Water,  before forcing
! #m0 deallocate          ( Wats_1 )                                   ! Soil Water,  after  forcing
! #m0 deallocate          ( Wats_d )                                   ! Soil Water          forcing
 
! #TB deallocate          ( etBump )                                   ! Bumps Layer Averaged Humidity
 
      deallocate          ( SoRnOF )                                   ! Soil     Run    OFF
      deallocate          ( Dhydtz )                                   ! Dhydif * dt / dz           [m]
      deallocate          ( Diag_A )                                   ! A Diagonal
      deallocate          ( Diag_B )                                   ! B Diagonal
      deallocate          ( Diag_C )                                   ! C Diagonal
      deallocate          ( Term_D )                                   !   Independant Term
      deallocate          ( Aux__P )                                   ! P Auxiliary Variable
      deallocate          ( Aux__Q )                                   ! Q Auxiliary Variable
      deallocate          ( etaaux )                                   ! Soil Water Content     [m3/m3]
 
      END IF                                                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
      return
      end subroutine SISVAT_qSo
 
 
 
      subroutine SISVAT_wEq( labWEq ,istart)
 
!--------------------------------------------------------------------------+
!   MAR          SISVAT_wEq                           Wed 26-Jun-2013  MAR |
!     SubRoutine SISVAT_wEq computes the Snow/Ice  Water  Equivalent       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Feb-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!                                                                          |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.) |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                     |
!     FILE                 |      CONTENT                                  |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
!   # SISVAT_wEq.ve        | #ve: OUTPUT/Verification: Snow/Ice Water Eqv. |
!                          |      unit 45, SubRoutine  SISVAT_wEq **ONLY** |
!--------------------------------------------------------------------------+
 
 
! General Variables
! =================
 
      use Mod_Real
      use Mod_SISVAT_ctr
      use Mod_SISVAT_grd
      use Mod_SISVAT_kkl
      use Mod_SISVAT_wEq
 
 
      IMPLICIT NONE
 
 
! Global Variables
! ================
 
      character(len=6)  ::  labWEq
      integer           ::  istart
 
 
! Local  Variables
! ================
 
      integer           ::  ikl,ikv   ,isn
      real(kind=real8)  ::  SnoWEQ,IceWEQ
 
 
! Switch Initialization
! =====================
 
      IF (.NOT.logWEq)                                              THEN
               logWEq = .true.
               open(unit=45,status='unknown',file='SISVAT_wEq.ve')
               rewind    45
      END IF
 
 
! Snow Water Equivalent
! =====================
 
           ikl   = 1
           ikv   = 1
      IF          (isnoSV(ikl,ikv).gt.iiceSV(ikl,ikv))              THEN
 
          SnoWEQ = 0.
        DO isn   = iiceSV(ikl,ikv)+1 ,isnoSV(ikl,ikv)
          SnoWEQ = SnoWEQ           + ro__SV(ikl,ikv,isn) * dzsnSV(ikl,ikv,isn)
        END DO
 
      END IF
 
 
! Ice  Water Equivalent
! =====================
 
      IF          (iiceSV(ikl,ikv).gt.0)                            THEN
 
          IceWEQ = 0.
        DO isn   = 1      , iiceSV(ikl,ikv)
          IceWEQ = IceWEQ + ro__SV(ikl,ikv,isn) * dzsnSV(ikl,ikv,isn)
        END DO
 
      END IF
 
 
! OUTPUT
! ======
 
      IF (istart.eq.1)                                              THEN
        write(45,45)dahost,i___SV(lwriSV(1,1)),j___SV(lwriSV(1,1))     &
     &                    ,n___SV(lwriSV(1,1))
 45     format(a18,10('-'),'Pt.',3i4,60('-'))
      END IF
 
      write(45,450) labWEq,IceWEQ,iiceSV(ikl,ikv),SnoWEQ               &
     &                    ,IceWEQ+SnoWEQ,isnoSV(ikl,ikv)               &
     &                                  ,drr_SV(ikl,ikv)*dt__SV        &
     &                                  ,dsn_SV(ikl,ikv)*dt__SV        &
     &                                  ,BufsSV(ikl,ikv)
 450  format(a6,3x,'  I+S =',f11.4,'(',i2,') +',f11.4,' =',            &
     &                       f11.4,'(',i2,')',                         &
     &             '  drr =', f7.4,                                    &
     &             '  dsn =', f7.4,                                    &
     &             '  Buf =', f7.4)
 
 
      return
      end subroutine SISVAT_wEq
