      module Mod_genTKE_RUN

!--------------------------------------------------------------------------+
!                                                     Mon 17-Jun-2013  MAR |
!     module Mod_genTKE_RUN contains local variables of PHY_genTKE_RUN     |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Mon 17-Jun-2013      |
!           Last Modification by H. Gallee,           Mon 17-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


! Global  Variables
! =================

      use Mod_Real



!  Local  Variables
!  ================

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  dukkp1    !      Difference (u(k) - u(k+1))                      [m/s]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  dvkkp1    !      Difference (v(k) - v(k+1))                      [m/s]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  kkp1dz    !  1 / Difference (Z(k) - Z(k+1))                      [1/m]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  zShear    !      Wind Shear Contribution to TKE                [m2/s3]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  REq_PT    !  Reduced (Equivalent) Potential Temperature            [K]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  c_Buoy    !  Buoyancy Coefficient (g/theta) X (dtheta/dz)       [1/s2]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  Ri__Nb    !  Richardson Number                                     [-]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  Prandtl   !  Prandtl    Number (Kzm/Kzh)                           [-]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  Ls_inv    !  1 / Ls                      (Therry & Lacarr, 1983) [1/m]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  ML_inv    !  1 / ML  (Mixing      Length, Therry & Lacarr, 1983) [1/m]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  DL_inv    !  1 / DL  (Dissipation Length, Therry & Lacarr, 1983) [1/m]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  Dissip    !           Dissipation                              [m2/s3]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  TKEvav    !           TKE         Vertical moving Average      [m2/s2]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  epsvav    !           Dissipation Vertical moving Average      [m2/s3]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  pkt       !           Reduced     Potential       Temperature      [X]


      end module Mod_genTKE_RUN
