      module Mod_SISVAT_cdf

!--------------------------------------------------------------------------+
!                                                     Thu 16-May-2013  MAR |
!     module Mod_SISVAT_cdf contains Surface Energy Balance   Variables of |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                OUTPUT on a netcdf file by SISVAT stand alone version     |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Thu 16-May-2013      |
!           Last Modification by H. Gallee,           Thu 16-May-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE


! SISVAT INPUT        Variables
! -----------------------------



! SISVAT INPUT/OUTPUT Variables
! -----------------------------



! SISVAT OUTPUT       Variables
! -----------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  SOsoNC_xyn  ! Absorbed Solar    Radiation               [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  IRsoNC_xyn  ! Absorbed IR       Radiation               [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  HSsoNC_xyn  ! Absorbed Sensible Heat Flux               [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  HLsoNC_xyn  ! Absorbed Latent   Heat Flux               [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  HLs_NC_xyn  ! Evaporation                          [mm w.e./s]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  HLv_NC_xyn  ! Transpiration                             [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  eta_NC_xyn  ! Soil              Humidity               [m3/m2]


      end module Mod_SISVAT_cdf
