      module Mod_PHY_RT_xyz

!--------------------------------------------------------------------------+
!                                                     Thu  9-May-2013  MAR |
!     module Mod_PHY_RT contains the representation in cartesian frame of  |
!                       Variables of the   Radiative Transfer Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Tue  5-Mar-2013      |
!           Last Modification by H. Gallee,           Thu  9-May-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE


! Radiative Transfer INPUT        Variables
! -----------------------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  LWUsRT_xy    !  Surface LongWave Heat Flux (+)    (Upward)   [W/m2]


! Radiative Transfer INPUT/OUTPUT Variables
! -----------------------------------------



! Radiative Transfer OUTPUT       Variables
! -----------------------------------------

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  O3__RT_xyz   !  Ozone    Concentration                             [Pa/Pa]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:,:)::  AersRT_xyza  !  Aerosol  Optical Thickness                             [-]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  ClouRT_xy    !  Total Cloudiness above lowest Atmospheric Level        [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  ODC_RT_xy    !  Clouds   Optical Thickness (vertically integrated)     [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  ODC_RT_xyz   !  Clouds   Optical Thickness (Layer z)                   [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  ODA_RT_xy    !  Aerosols Optical Thickness (vertically integrated)     [-]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:,:)  ::  ODA_RT_xyz   !  Aerosols Optical Thickness (Layer z)                   [-]

      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  OLR_RT_xy    !  OutgoingLongWave Radiation (+)  (  Upward)          [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  SWDsRT_xy    !  Surface ShrtWave Heat Flux (+)  (Downward)          [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  SWAsRT_xy    !  Surface ShrtWave Heat Flux (+)  (Absorbed)          [W/m2]
      real(kind=real8), SAVE,ALLOCATABLE ,dimension(:,:)    ::  LWDsRT_xy    !  Surface LongWave Heat Flux (+)  (Downward)          [W/m2]



      end module Mod_PHY_RT_xyz
