      module Mod_SISVATLqSn

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVATLqSn contains local variables of SISVAT_qSn         |
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

! OUTPUT/Verification: Energy/Water Budget
! #e1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  EqSn_d          ! Energy in Excess, initial
! #e1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  EqSn_0          ! Snow Energy, befor Phase Change
! #e5 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  EqSn01          ! Snow Energy, after Phase Change
! #e5 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  EqSn02          ! Snow Energy, after Phase Change
                                                                       !              .AND. Last Melting 
! #e1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  EqSn_1          ! Snow Energy, after Phase Change
                                                                       !              .AND. Mass Redistr.
! OUTPUT/Verification: * Mass Conservation
! #m1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  SIsubl          ! Snow Deposed Mass
! #m1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  SImelt          ! Snow Melted  Mass
! #m1 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  SIrnof          ! Local Surficial Water + Run OFF


      integer, SAVE         , ALLOCATABLE, dimension(:,:)  ::  noSnow          ! Nb of Layers Updater
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  EExdum          ! Energy in Excess when no Snow
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  dzMelt          ! Melted    Thickness          [m]

! OUTPUT/Verification: Energy/Water Budget
! #e5 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  WqSn_0          ! Snow Water+Forcing  Initial
! #e5 real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)  ::  WqSn_1          ! Snow Water+Forcing, Final


      end module Mod_SISVATLqSn
