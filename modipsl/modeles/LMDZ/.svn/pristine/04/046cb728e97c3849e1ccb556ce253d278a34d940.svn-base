      module Mod_SISVATLBSn

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVATLBSn contains local variables of SISVAT_BSn         |
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



      integer, SAVE         , ALLOCATABLE, dimension(:,:)    ::  Mobile                        !

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  SaltSI                        ! Snow Drift Index
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  sdrift                        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  xdrift                        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  zdrift                        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  tdepos                        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  zdepos                        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  dbsaux                        ! Drift Amount   (Dummy Variable)
      integer, SAVE         , ALLOCATABLE, dimension(:,:)    ::  isagr1                        ! 1st     Layer History
      integer, SAVE         , ALLOCATABLE, dimension(:,:)    ::  isagr2                        ! 2nd     Layer History

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  WEagre                        ! Snow Water Equivalent Thickness
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  Agrege                        ! 1. when Agregation constrained
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  dzagr1                        ! 1st     Layer Thickness
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  dzagr2                        ! 2nd     Layer Thickness
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  T_agr1                        ! 1st     Layer Temperature
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  T_agr2                        ! 2nd     Layer Temperature
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  roagr1                        ! 1st     Layer Density
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  roagr2                        ! 2nd     Layer Density
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  etagr1                        ! 1st     Layer Water Content
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  etagr2                        ! 2nd     Layer Water Content
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  G1agr1                        ! 1st     Layer Dendricity/Spher.
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  G1agr2                        ! 2nd     Layer Dendricity/Spher.
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  G2agr1                        ! 1st     Layer Sphericity/Size
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  G2agr2                        ! 2nd     Layer Sphericity/Size
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  agagr1                        ! 1st     Layer Age
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  agagr2                        ! 2nd     Layer Age



      end module Mod_SISVATLBSn
