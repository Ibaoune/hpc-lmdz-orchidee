      module Mod_Atm_CM_QSa

!--------------------------------------------------------------------------+
!                                                     Mon 17-Jun-2013  MAR |
!     module Mod_Atm_CM_QSa contains local variables of PHY_Atm_CM_QSat    |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Mon 17-Jun-2013      |
!           Last Modification by H. Gallee,           Mon 17-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


! Global  Variables
! =================

      use Mod_Real



! Internal Variables
! ==================

      IMPLICIT NONE



      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  pa_hPa                  !  Pressure                               [hPa]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  pr_b75                  !  (1/0)          if pa_hPa (>/<)      75 [hPa]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  ei_sat                  !  Saturation Vapor Pressure over Ice     [hPa]
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:)  ::  ew_sat                  !  Saturation Vapor Pressure over Water   [hPa]


      end module Mod_Atm_CM_QSa
