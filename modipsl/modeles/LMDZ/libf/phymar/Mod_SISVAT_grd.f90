      module Mod_SISVAT_grd


!--------------------------------------------------------------------------+
!                                                     Fri  7-Jun-2013  MAR |
!     module Mod_SISVAT_grd contains the dimensions of the domain of       |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Wed 27-Feb-2013      |
!                    modified by H. Gallee,           Fri  7-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



      integer, SAVE                                       ::  mwp         ! Nb of mosaic in one         Grid Cell
      integer, SAVE                                       ::  n2          ! Nb of mosaic in one Oceanic Grid Cell (min is 2)
      integer, SAVE                                       ::  kcolv       ! Nb of interior      Vertical Columns
      integer, SAVE                                       ::  nsoil       ! Nb of Soil          Levels beneath the Surface Level
      integer, SAVE                                       ::  nvege       ! Nb of Vegetation    Types
      integer, SAVE                                       ::  nsnow       ! Max Nb of Snow      Layers
      integer, SAVE                                       ::  nbPts       ! Nb of dumped        Grid     Points 
      integer, SAVE                                       ::  nbwri       ! Nb of dumped        Vertical Columns
      integer, SAVE                                       ::  ntave       ! Nb of Time Steps over which SBL relevant parameters  are averaged
      integer, SAVE                                       ::  ntavz       ! Nb of Time Steps over which z0, r0, ...  parameters  are averaged
      integer, SAVE                                       ::  nLimi       ! Nb of Time Steps over which Water Vapor Flux to limit is averaged

      integer, SAVE                                       ::  jt__SV      ! Number of DYnamical Time Steps for one SISVAT                       [-]
      real(kind=real8), SAVE                              ::  dt__SV      ! Time Step of Surface     Physics      (SISVAT)                      [s]

      integer, SAVE                                       ::  k_zb
      real(kind=real8), SAVE                              ::  z_zb = 25.  ! Level of negligible blowing particles concentration

      integer, SAVE                                       ::  k_SL        ! Parameter used in the Interpolation of V(10 m)
      real(kind=real8), SAVE                              ::  r_SL10      ! Parameter used in the Interpolation of V(10 m)

      integer, SAVE                                       ::  iwr_SV
      integer, SAVE                                       ::  jwr_SV
      integer, SAVE                                       ::  nwr_SV
      integer, SAVE         ,ALLOCATABLE ,dimension(:)    ::  ii__SV      ! Mosaic point   i Coordinate
      integer, SAVE         ,ALLOCATABLE ,dimension(:)    ::  jj__SV      ! Mosaic point   j Coordinate
      integer, SAVE         ,ALLOCATABLE ,dimension(:)    ::  nn__SV      ! Mosaic point   n Coordinate
      integer, SAVE         ,ALLOCATABLE ,dimension(:)    ::  ikp_SV      ! Grid Cell Column Index of a SISVAT Column
      integer, SAVE         ,ALLOCATABLE ,dimension(:,:,:)::  ikl_SV      ! SISVAT    Column Index     

      end module Mod_SISVAT_grd
