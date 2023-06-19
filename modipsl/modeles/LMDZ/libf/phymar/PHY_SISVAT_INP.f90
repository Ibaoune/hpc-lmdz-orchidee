      subroutine PHY_SISVAT_INP(FlagSV_SBC,TsHOST,kcolq)

!------------------------------------------------------------------------------+
!                                                         Mon  1-Jul-2013  MAR |
!     SubRoutine PHY_SISVAT_INP initializes                                    |
!                    SISVAT    (Soil                                           |
!                               Ice                                            |
!                               Snow                                           |
!                               Vegetation                                     |
!                               Atmosphere                                     |
!                               Transfer   Scheme)                             |
!                                                                              |
!                               USING DATA                                     |
!                                                                              |
!     Here DATA are academic                                                   |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Wed  1-May-2013      |
!           Last Modification by H. Gallee,               Mon  1-Jul-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


! Global  Variables
! =================

      use Mod_Real
      use Mod_PHY____grd
      use Mod_PHY____dat
      use Mod_SISVAT_dat
      use Mod_SISVAT_grd
      use Mod_SISVAT_kkl
      use Mod_SISVAT_gpt


      IMPLICIT NONE


      logical                                               ::  FlagSV_SBC         !  DATA     Surface   Flag
      integer                                               ::  kcolq              !
      real(kind=real8), dimension(kcolq)                    ::  TsHOST             ! (Initial) Surface   Temperature                          [K]
      



! Local   Variables
! =================

! FOR INPUT
! ---------

!     integer         ,ALLOCATABLE ,dimension(:,:,:)  ::  Mask___xyn  ! Land(1)-Sea(0) Mask         (Mosaic)         [-]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  Frac___xyn  ! Grid Cell Fraction          (Mosaic)         [-]
!     integer         ,ALLOCATABLE ,dimension(:,:,:)  ::  iVgT___xyn  ! Vegetat. Type               (Mosaic)         [-]
!     integer         ,ALLOCATABLE ,dimension(:,:,:)  ::  iSoT___xyn  !     Soil Type               (Mosaic)         [-]
!     integer         ,ALLOCATABLE ,dimension(:,:,:)  ::  iDra___xyn  !     Soil Drainage Switch    (Mosaic)         [-]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  Albs___xyn  ! Dry Soil Albedo             (Mosaic)         [-]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  LAI____xyn  ! Leaf Area Index             (Mosaic)     [m2/m2]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  GLF____xyn  ! Green Leaf Fraction         (Mosaic)         [-]



! FOR INPUT/OUTPUT (RESTART)
! ----------------

!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  FSbz___xyn  ! Fallen Snow Buffer          (Mosaic)   [mm w.e.]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  Tveg___xyn  ! Vegetation      Temperature (Mosaic)         [K]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  Sveg___xyn  ! Canopy intercepted Snow     (Mosaic)   [mm w.e.]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  Rveg___xyn  ! Canopy intercepted Rain     (Mosaic)   [mm w.e.]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  pveg___xyn  ! Vegetation Water Potential  (Mosaic)         [m]

!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  z0_____xyn  ! Roughness Length Momentum   (Mosaic)         [m]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  z0ro___xyn  ! Orographic Roughness Length (Mosaic)         [m]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  r0_____xyn  ! Roughness Length Scalars    (Mosaic)         [m]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  uss____xyn  ! Blown    Snow   Flux Turbulent Scale [kg/kg m/s]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  uds____xyn  ! Blown    Dust   Flux Turbulent Scale [kg/kg m/s]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  BSbz___xyn  ! Eroded Snow Buffer          (Mosaic)   [mm w.e.]
! #AE real(kind=real8),ALLOCATABLE ,dimension(:,:)    ::  usth___xyn  ! Friction Velocity Saltation Threshold      [m/s]

!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  Alb____xyn  ! Surface  Albedo             (Mosaic)         [-]
!     integer         ,ALLOCATABLE ,dimension(:,:,:)  ::  SNnb___xyn  ! Nb of  Snow Layers          (Mosaic)         [-]
!     integer         ,ALLOCATABLE ,dimension(:,:,:)  ::  SInb___xyn  ! Nb of S.Ice Layers          (Mosaic)         [-]
!     integer         ,ALLOCATABLE ,dimension(:,:,:)  ::  ICnb___xyn  ! Nb of   Ice Layers          (Mosaic)         [-]

!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  SnWE___xyn  ! Snow  Pack  Thickness       (Mosaic)   [mm w.e.]

!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  wes____xyn  ! Sublimation        Budget   (Mosaic)   [mm w.e.]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  wer____xyn  ! Melting            Budget   (Mosaic)   [mm w.e.]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  wem____xyn  ! Refreezing         Budget   (Mosaic)   [mm w.e.]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  wee____xyn  ! Evapotranspiration Budget   (Mosaic)   [mm w.e.]

!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:,:)::  Tsis__dxyn  ! Soil-Ice-Snow   Temperature (Mosaic)         [K]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:,:)::  Qsis__dxyn  ! Soil-Ice-Snow   Humidity    (Mosaic)     [kg/kg]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:,:)::  SNdz__dxyn  ! Snow  Layer Thickness       (Mosaic)         [m]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:,:)::  SNro__dxyn  ! Snow  Layer Density         (Mosaic)     [kg/m3]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:,:)::  SNag__dxyn  ! Snow  Layer Age             (Mosaic)       [day]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:,:)::  SNg1__dxyn  ! Snow  Layer Property (d./s.)(Mosaic)         [-]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:,:)::  SNg2__dxyn  ! Snow  Layer Property (s./s.)(Mosaic)         [-]
!     integer         ,ALLOCATABLE ,dimension(:,:,:,:)::  SNhi__dxyn  ! Snow  Layer History         (Mosaic)         [-]

!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  SWz____xyn  ! Surficial Water Mass        (Mosaic)   [mm w.e.]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  SWx____xyn  ! Surficial Water Status      (Mosaic)   [0,1=m,f]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  SWHy___xyn  ! Surficial Water: Normaliz.Decay (Mosaic)     [-]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  ROdt___xyn  ! Run-Off  Intensity          (Mosaic) [mm w.e./s]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  ROF____xyn  ! Cumulative Run-Off          (Mosaic)   [mm w.e.]


! A-O Coupling
! ------------

!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  Alb_AO_xyn  ! Ocean COUPLING (Surface Albedo) n=2 sea-ice  [-]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  s_T_AO_xyn  ! Ocean COUPLING (SST)            n=1 opn-watr [K]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  dSdTAO_xyn  ! Ocean COUPLING (d(SH Flux) / dT)        [W/m2/K]
!     real(kind=real8),ALLOCATABLE ,dimension(:,:,:)  ::  dLdTAO_xyn  ! Ocean COUPLING (d(LH Flux) / dT)        [W/m2/K]


! Dummy
! -----

      integer  ::  i     ,j     ,n     ,l     ,nm    ,ikp






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!     NON ACADEMIC INITIALIZATION                                      !
!                                                                      !
      IF (FlagSV_SBC)                                               THEN
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! +++++++++++++++++++++++++++++++++++++++++++++++++++++
! INITIALIZATION of a non-1D case, by using INPUT Files
!                                    (HERE FROM NESTOR)

! Martin CONTROL
Print*,'mxpp=',mxpp
Print*,'mypp=',mypp
Print*,'ixp1=',ixp1
Print*,'jyp1=',jyp1

        IF   (mxpp-ixp1.GT.0 .OR. mypp-jyp1.GT.0)                   THEN
! +++++++++++++++++++++++++++++++++++++++++++++++++++++




! ================================
! FIRST INITIALIZATION     (BEGIN)

          IF       (it_EXP .LE. 1)                                  THEN
! ================================


! ================================
          END IF ! (it_EXP .LE. 1)

! FIRST INITIALIZATION       (END)
! ================================




! Update Sea-Ice    Fraction
! ==========================

! #IP     IF (             inpSBC                    )              THEN

!                ******
! #IP       call INIsic(newsicSI )  ! ATTENTION 2 premiers Arguments supprimés
!                ******

! #IP     END IF




! Update LAI and/or Green Leaf Fraction
! =====================================

! #GP                                          glfFIX = .true.
! #CM     IF (vegmod .AND. inpSBC .AND.  .NOT. glfFIX)              THEN

!                ******
! #CM       call INIglf(newglfSIS)  ! ATTENTION 2 premiers Arguments supprimés
!                ******

! #CM     END IF






! +++++++++++++++++++++++++++++++++++++++++++++++++++++
! INITIALIZATION of a     1D case, by prescribing INPUT

        ELSE
! +++++++++++++++++++++++++++++++++++++++++++++++++++++




! ================================
! FIRST INITIALIZATION     (BEGIN)

          IF       (it_EXP .LE. 1)                                  THEN
! ================================



! ALLOCATION
! ==========

! FOR INPUT
! ---------

!         allocate  ( Mask___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Land(1)-Sea(0) Mask         (Mosaic)         [-]
!         allocate  ( Frac___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Grid Cell Fraction          (Mosaic)         [-]
!         allocate  ( iVgT___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Vegetat. Type               (Mosaic)         [-]
!         allocate  ( iSoT___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !      Soil Type               (Mosaic)         [-]
!         allocate  ( iDra___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !      Soil Drainage Switch    (Mosaic)         [-]
!         allocate  ( Albs___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Dry Soil Albedo             (Mosaic)         [-]
!         allocate  ( LAI____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Leaf Area Index             (Mosaic)     [m2/m2]
!         allocate  ( GLF____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Green Leaf Fraction         (Mosaic)         [-]


! FOR INPUT/OUTPUT (RESTART)
! ----------------

!         allocate  ( FSbz___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Fallen Snow Buffer          (Mosaic)   [mm w.e.]
!         allocate  ( Tveg___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Vegetation      Temperature (Mosaic)         [K]
!         allocate  ( Sveg___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Canopy intercepted Snow     (Mosaic)   [mm w.e.]
!         allocate  ( Rveg___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Canopy intercepted Water    (Mosaic)   [mm w.e.]
!         allocate  ( pveg___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Vegetation Water Potential  (Mosaic)         [m]

!         allocate  ( z0_____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Roughness Length Momentum   (Mosaic)         [m]
!         allocate  ( z0ro___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  z0(Orographic)                               [m]
!         allocate  ( r0_____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Roughness Length Scalars    (Mosaic)         [m]
!         allocate  ( uss____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Blown    Snow   Flux Turbulent Scale [kg/kg m/s]
!         allocate  ( BSbz___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Eroded Snow Buffer          (Mosaic)   [mm w.e.]
! #AE     allocate  ( usth___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Friction Velocity Saltation Threshold      [m/s]

!         allocate  ( Alb____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Surface  Albedo             (Mosaic)         [-]
!         allocate  ( SNnb___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Nb of  Snow Layers          (Mosaic)         [-]
!         allocate  ( SInb___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Nb of S.Ice Layers          (Mosaic)         [-]
!         allocate  ( ICnb___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Nb of   Ice Layers          (Mosaic)         [-]

!         allocate  ( SnWE___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Snow  Pack  Thickness       (Mosaic)   [mm w.e.]

!         allocate  ( wes____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Sublimation        Budget   (Mosaic)   [mm w.e.]
!         allocate  ( wer____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Melting            Budget   (Mosaic)   [mm w.e.]
!         allocate  ( wem____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Refreezing         Budget   (Mosaic)   [mm w.e.]
!         allocate  ( wee____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Evapotranspiration Budget   (Mosaic)   [mm w.e.]

!         allocate  ( Tsis__dxyn(ixp1:mxpp,jyp1:mypp,mwp,-nsoil:nsnow) )
!                                                             !  Soil-Ice-Snow   Temperature (Mosaic)         [K]
!         allocate  ( Qsis__dxyn(ixp1:mxpp,jyp1:mypp,mwp,-nsoil:nsnow) )
!                                                             !  Soil-Ice-Snow   Humidity    (Mosaic)     [kg/kg]
!         allocate  ( SNdz__dxyn(ixp1:mxpp,jyp1:mypp,mwp,nsnow) ) !  Snow  Layer Thickness       (Mosaic)         [m]
!         allocate  ( SNro__dxyn(ixp1:mxpp,jyp1:mypp,mwp,nsnow) ) !  Snow  Layer Density         (Mosaic)     [kg/m3]
!         allocate  ( SNag__dxyn(ixp1:mxpp,jyp1:mypp,mwp,nsnow) ) !  Snow  Layer Age             (Mosaic)       [day]
!         allocate  ( SNg1__dxyn(ixp1:mxpp,jyp1:mypp,mwp,nsnow) ) !  Snow  Layer Property (d./s.)(Mosaic)         [-]
!         allocate  ( SNg2__dxyn(ixp1:mxpp,jyp1:mypp,mwp,nsnow) ) !  Snow  Layer Property (s./s.)(Mosaic)         [-]
!         allocate  ( SNhi__dxyn(ixp1:mxpp,jyp1:mypp,mwp,nsnow) ) !  Snow  Layer History         (Mosaic)         [-]

!         allocate  ( SWz____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Surficial Water Mass        (Mosaic)   [mm w.e.]
!         allocate  ( SWx____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Surficial Water Status      (Mosaic)   [0,1=m,f]
!         allocate  ( SWHy___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Normalized Decay of Surficial Water Content  [-]
!         allocate  ( ROdt___xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Run-Off  Intensity          (Mosaic) [mm w.e./s]
!         allocate  ( ROF____xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Cumulative Run-Off          (Mosaic)   [mm w.e.]


! A-O Coupling
! ------------

!         allocate  ( Alb_AO_xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Ocean COUPLING (Surface Albedo) n=2 sea-ice  [-]
!         allocate  ( Alb_AO_xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Ocean COUPLING (Surface Albedo) n=2 sea-ice  [-]
!         allocate  ( dSdTAO_xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Ocean COUPLING (d(SH Flux) / dT)        [W/m2/K]
!         allocate  ( dLdTAO_xyn(ixp1:mxpp,jyp1:mypp,mwp) )       !  Ocean COUPLING (d(LH Flux) / dT)        [W/m2/K]


! Martin CONTROL
PRINT*,'On initialise des paramètres en dur dans PHY_SISVAT_INP'

! ------------------------------------------
! From 2-D Horizontal to 1-D Horizontal Grid
! ------------------------------------------

          DO ikp=1,kcolp
             i   = ii__AP(ikp)
             j   = jj__AP(ikp)


            !MaskSV_gpt(ikp)     = 1                     !  Land(1)-Sea(0) Mask         (Grid Cell)      [-]
            ! Martin : aquaplanète
             MaskSV_gpt(ikp)     = 0                     !  Land(1)-Sea(0) Mask         (Grid Cell)      [-]

          DO n=1,mwp
!           Mask___xyn(i,j,n)   = 1                     !  Land(1)-Sea(0) Mask         (Mosaic)         [-]
            !LSmask    (ikp,n)   = 1                     !  Land(1)-Sea(0) Mask         (Mosaic)         [-]
            LSmask    (ikp,n)   = 0                     !  Land(1)-Sea(0) Mask         (Mosaic)         [-]

            SlopSV    (ikp,n)   = 0.0                   !  Assignation of Slope to a Particular Value

!           Frac___xyn(i,j,n)   = 1.0/float(mwp)        !  Grid Cell Fraction          (Mosaic)         [-]
            FracSV    (ikp,n)   = 1.0/float(mwp)        !  Grid Cell Fraction          (Mosaic)         [-]

!           z0_____xyn(i,j,n)   = 1.e-6                 !  Roughness Length Momentum   (Mosaic)         [m]   
            !z0__SV    (ikp,n)   = 1.e-6                 !  Roughness Length Momentum   (Mosaic)         [m]   
            z0__SV    (ikp,n)   = 1.e-3                 !  Roughness Length Momentum   (Mosaic)         [m]   
                                                      ! (Sensitive if non-zero Orographic Roughness)
!           r0_____xyn(i,j,n)   = 1.e-6                 !  Roughness Length Scalars    (Mosaic)         [m]
            !Z0h_SV    (ikp,n)   = 1.e-6                 !  Roughness Length Scalars    (Mosaic)         [m]
            Z0h_SV    (ikp,n)   = 1.e-3                 !  Roughness Length Scalars    (Mosaic)         [m]


! Soil and Vegetation
! -------------------

!           iVgT___xyn(i,j,n)   =   5                   !  Vegetat. Type               (Mosaic)         [-]   (Sensitive)
            iVgTSV    (ikp,n)   =   0                   !  Vegetat. Type               (Mosaic)         [-]   (Sensitive)
!           iSoT___xyn(i,j,n)   =   5                   !      Soil Type               (Mosaic)         [-]
            iSoTSV    (ikp,n)   =   0                   !      Soil Type               (Mosaic)         [-]
!           iDra___xyn(i,j,n)   =   1                   !      Soil Drainage Switch    (Mosaic)         [-]  
            iWaFSV    (ikp,n)   =   1                   !      Soil Drainage Switch    (Mosaic)         [-]  
!           Albs___xyn(i,j,n)   =   0.15                !  Dry Soil Albedo             (Mosaic)         [-]
            alb0SV    (ikp,n)   =   0.15                !  Dry Soil Albedo             (Mosaic)         [-]
!           LAI0SV    (ikp,n)   = LAI____xyn(i,j,n)     ! LAI
!           glf0SV    (ikp,n)   = GLF____xyn(i,j,n)     ! Green Leaf Frac.


! Snow Pack Initialization: NO Snow Pack
! --------------------------------------

!           SNnb___xyn(i,j,n)   = 0                     !  Nb of  Snow Layers          (Mosaic)         [-]   (Sensitive, maybe)
            isnoSV    (ikp,n)   = 0                     !  Nb of  Snow Layers          (Mosaic)         [-]   (Sensitive, maybe)
!           SInb___xyn(i,j,n)   = 0                     !  Nb of S.Ice Layers          (Mosaic)         [-]
            ispiSV    (ikp,n)   = 0                     !  Nb of S.Ice Layers          (Mosaic)         [-]
!           ICnb___xyn(i,j,n)   = 0                     !  Nb of   Ice Layers          (Mosaic)         [-]   (Sensitive, maybe)
            iiceSV    (ikp,n)   = 0                     !  Nb of   Ice Layers          (Mosaic)         [-]   (Sensitive, maybe)

          DO l=-nsol,0
!           Tsis__dxyn(i,j,n,l) =     TsHOST(ikp)       !  Soil-Ice-Snow   Temperature (Mosaic)         [K]
!           TsisSV    (ikp,n,l) = Tsis__dxyn(i,j,n,l)   !
            TsisSV    (ikp,n,l) =     TsHOST(ikp)       !
!           Qsis__dxyn(i,j,n,l) =                      &!                  Humidity    (Mosaic)     [kg/kg]
!    &          etadSV(iVgTSV(ikp,n))                   ! (Saturation is   assumed)
!           eta_SV    (ikp,n,l) = Qsis__dxyn(i,j,n,l)   !
            eta_SV    (ikp,n,l) =                      &!
     &          etadSV(iVgTSV(ikp,n))                   !
          ENDDO

          DO l=1,nsno
!           Tsis__dxyn(i,j,n,l) = 273.15                !  Soil-Ice-Snow   Temperature (Mosaic)         [K]
!           TsisSV    (ikp,n,l) = Tsis__dxyn(i,j,n,l)   !
            TsisSV    (ikp,n,l) = 273.15                !
!           Qsis__dxyn(i,j,n,l) =   0.00                !  Soil-Ice-Snow   Humidity    (Mosaic)     [kg/kg]
!           eta_SV    (ikp,n,l) = Qsis__dxyn(i,j,n,l)   !
            eta_SV    (ikp,n,l) =   0.00                !
!           SNdz__dxyn(i,j,n,l) =   0.00                !  Snow  Layer Thickness       (Mosaic)         [m]
!           dzsnSV    (ikp,n,l) = SNdz__dxyn(i,j,n,l)   !
            dzsnSV    (ikp,n,l) =   0.00                !
!           SNro__dxyn(i,j,n,l) =   0.00                !  Snow  Layer Density         (Mosaic)     [kg/m3]
!           ro__SV    (ikp,n,l) = SNro__dxyn(i,j,n,l)   !
            ro__SV    (ikp,n,l) =   0.00                !
!           SNag__dxyn(i,j,n,l) =   0.00                !  Snow  Layer Age             (Mosaic)       [day]
!           agsnSV    (ikp,n,l) = SNag__dxyn(i,j,n,l)   !
            agsnSV    (ikp,n,l) =   0.00                !
!           SNg1__dxyn(i,j,n,l) =   0.00                !  Snow  Layer Property (d./s.)(Mosaic)         [-]
!           G1snSV    (ikp,n,l) = SNg1__dxyn(i,j,n,l)   !
            G1snSV    (ikp,n,l) =   0.00                !
!           SNg2__dxyn(i,j,n,l) =   0.00                !  Snow  Layer Property (s./s.)(Mosaic)         [-]
!           G2snSV    (ikp,n,l) = SNg2__dxyn(i,j,n,l)   !
            G2snSV    (ikp,n,l) =   0.00                !
!           SNhi__dxyn(i,j,n,l) =   0                   !  Snow  Layer History         (Mosaic) 
!           istoSV    (ikp,n,l) = SNhi__dxyn(i,j,n,l)   !
            istoSV    (ikp,n,l) =   0                   !
          ENDDO

          ENDDO
          ENDDO


! ================================
          END IF ! (it_EXP .LE. 1)

! FIRST INITIALIZATION       (END)
! ================================




! OUTPUT
! ======

!             dSdTAO_xyn(i,j,n) = dSdTSV(ikp,n)                         ! Sens.H.Flux T-Der.  (A-O Coupling)  [W/m2/K]
!             dLdTAO_xyn(i,j,n) = dLdTSV(ikp,n)                         ! Latn.H.Flux T-Der.  (A-O Coupling)  [W/m2/K]
!!#SN         SnWE___xyn(i,j,n) = zWE_SV(ikp,n)                         ! Current  *Thick.
!             ROdt___xyn(i,j,n) = RnofSV(ikp,n)                         ! Run OFF Intensity




! RESTART
! =======

!!#SN         BufsSV(ikp,n) = FSbz___xyn(i,j,n)                         ! Snow Buffer Lay.
!!#SN         rusnSV(ikp,n) = SWz____xyn(i,j,n)                         ! Surficial Water
!!#SN         SWS_SV(ikp,n) = SWx____xyn(i,j,n)                         ! Surficial Water Status
!!#SN         SWf_SV(ikp,n) = SWHySV_xyn(i,j,n)                         ! Normalized Decay




! +++++++++++++++++++++++++++++++++++++++++++++++++++++
        END IF
! +++++++++++++++++++++++++++++++++++++++++++++++++++++






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!     ACADEMIC INITIALIZATION                                          !
!                                                                      !
      ELSE ! IF (.NOT.FlagSV_SBC)                                      !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! ==============================
! FIRST INITIALIZATION   (BEGIN)
        PRINT*, 'CONTROL PHY_SISVAT_INP'
        PRINT*,'kcolq=',kcolq
        PRINT*,'kcolp=',kcolp
        PRINT*,'kmwp=',mwp

        IF       (it_EXP .LE. 1)                                    THEN
! ==============================

            DO ikp=1,kcolp

              MaskSV_gpt(ikp)      = 0                                  !   Land-Sea Mask

            DO nm=1,mwp

              !LSMask(ikp,nm)       = 1                                  !   Land-Sea Mask
              ! Test aquaplanète (Martin):
              LSMask(ikp,nm)       = 0                                  !   Land-Sea Mask
              SlopSV(ikp,nm)       = 0.                                 !   Surface Slope
              FracSV(ikp,nm)       = 1.      / float(mwp)               !   Mesh Fraction
              z0__SV(ikp,nm)       = 1.                              !   z0   Momentum
              z0h_SV(ikp,nm)       = 1.                            !   z0    Scalars
              !z0__SV(ikp,nm)       = 0.1                              !   z0   Momentum
              !z0h_SV(ikp,nm)       = 0.1                            !   z0    Scalars
              !z0__SV(ikp,nm)       = 0.001                              !   z0   Momentum
              !z0h_SV(ikp,nm)       = 0.00001                            !   z0    Scalars
              isnoSV(ikp,nm)       = 0                                  !   NO Snow
              ispiSV(ikp,nm)       = 0                                  !   NO Ice (Superimposed)
              iiceSV(ikp,nm)       = 0                                  !   NO Ice

             !iVgTSV(ikp,nm)       = 5                                  !   GRASS MEDIUM Vegetat. Type  (Mosaic)        [-] (see Mod_SISVAT_dat)
              iVgTSV(ikp,nm)       = 0                                  !   GRASS MEDIUM Vegetat. Type  (Mosaic)        [-] (see Mod_SISVAT_dat)
             !iSoTSV(ikp,nm)       = 5                                  !   LOAM  Soil            Type  (Mosaic)        [-] (see Mod_SISVAT_dat)
              iSoTSV(ikp,nm)       = 0                                  !   LOAM  Soil            Type  (Mosaic)        [-] (see Mod_SISVAT_dat)
              iWaFSV(ikp,nm)       = 1                                  !         Soil Drainage Switch  (Mosaic)(1=free)[-]
              alb0SV(ikp,nm)       = 0.15                               !   Dry   Soil Albedo                           [-]

            DO l=-nsol,0
              TsisSV(ikp,nm,l)     = TsHOST(ikp)                        !   Soil-Ice-Snow   Temperature (Mosaic)        [K]
              eta_SV(ikp,nm,l)     =                                   &!                   Humidity    (Mosaic)    [kg/kg]
     &        etadSV(iVgTSV(ikp,nm))                                    !  (Saturation is   assumed)
            ENDDO

            DO l=1,nsno
              TsisSV(ikp,nm,l)     = 273.15                             !  Soil-Ice-Snow   Temperature (Mosaic)         [K]
              eta_SV(ikp,nm,l)     =   0.00                             !  Soil-Ice-Snow   Humidity    (Mosaic)     [kg/kg]
              dzsnSV(ikp,nm,l)     =   0.00                             !  Snow  Layer Thickness       (Mosaic)         [m]
              ro__SV(ikp,nm,l)     =   0.00                             !  Snow  Layer Density         (Mosaic)     [kg/m3]
              agsnSV(ikp,nm,l)     =   0.00                             !  Snow  Layer Age             (Mosaic)       [day]
              G1snSV(ikp,nm,l)     =   0.00                             !  Snow  Layer Property (d./s.)(Mosaic)         [-]
              G2snSV(ikp,nm,l)     =   0.00                             !  Snow  Layer Property (s./s.)(Mosaic)         [-]
              istoSV(ikp,nm,l)     =   0.00                             !  Snow  Layer History         (Mosaic)         [-]
            ENDDO

            ENDDO
            ENDDO


! ==============================
        END IF ! (it_EXP .LE. 1)

! FIRST INITIALIZATION     (END)
! ==============================


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      END IF !  (.NOT.FlagSV_SBC)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      end subroutine PHY_SISVAT_INP
