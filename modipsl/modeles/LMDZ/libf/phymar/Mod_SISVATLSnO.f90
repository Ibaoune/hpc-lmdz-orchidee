      module Mod_SISVATLSnO

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVATLSnO contains local variables of SISVAT SnOpt       |
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


      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  coalb1                ! weighted Coalbedo, Vis.
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  coalb2                ! weighted Coalbedo, nIR 1
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  coalb3                ! weighted Coalbedo, nIR 2
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  sExt_1                ! Extinction Coeff., Vis.
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  sExt_2                ! Extinction Coeff., nIR 1
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  sExt_3                ! Extinction Coeff., nIR 2
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  SnOpSV                ! Snow Grain optical Size
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  alb1sv
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  alb2sv
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  alb3sv


      end module Mod_SISVATLSnO
