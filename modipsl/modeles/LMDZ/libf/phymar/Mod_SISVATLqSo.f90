      module Mod_SISVATLqSo

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVATLqSo contains local variables of SISVAT_qSo         |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Mon 17-Jun-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


! Global  Variables
! =================

      use Mod_Real



! Internal Variables
! ==================

      IMPLICIT NONE


! OUTPUT
! ------

! OUTPUT/Verification: H2O    Conservation
! #m0 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  Wats_0        ! Soil Water,  before forcing
! #m0 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  Wats_1        ! Soil Water,  after  forcing
! #m0 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  Wats_d        ! Soil Water          forcing


! #TB real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  etBump        ! Bumps Layer Averaged Humidity

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  SoRnOF        ! Soil     Run    OFF
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Dhydtz        ! Dhydif * dt / dz           [m]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Diag_A        ! A Diagonal
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Diag_B        ! B Diagonal
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Diag_C        ! C Diagonal
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Term_D        !   Independant Term
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Aux__P        ! P Auxiliary Variable
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Aux__Q        ! Q Auxiliary Variable
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  etaaux        ! Soil Water Content     [m3/m3]

! OUTPUT/Verification: H2O    Conservation
! #mw real(kind=real8), SAVE                               ::  hourwr
! #mw real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  Evapor        !


      end module Mod_SISVATLqSo
