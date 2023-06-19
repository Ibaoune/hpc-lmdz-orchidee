      module Mod_SISVATLVgO

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVATLVgO contains local variables of SISVAT VgOptP      |
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

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  k___sv
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  A0__sv
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  gamasv
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  Sigcsv
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  C1__sv
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  C2__sv
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  criLAI


      end module Mod_SISVATLVgO
