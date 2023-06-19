
      module Mod_PHY_CP_ctr


!--------------------------------------------------------------------------+
!                                                     Mon 22-Apr-2013  MAR |
!     module Mod_PHY_CP_ctr contains the control variables used to set up  |
!            Convective Mass Flux Parameterization                         |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Wed 10-Apr-2013      |
!           Last Modification by H. Gallee,           Mon 22-Apr-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


      use Mod_Real


      IMPLICIT NONE



      logical           ::  Lod_CP            !   Deep    Convection (user)      Switch
      logical           ::  Odeep0 = .TRUE.   !   Deep    Convection (Control)   Switch
      logical           ::  Odeep             !   Deep    Convection (Effective) Switch
      logical           ::  Oshal0 = .TRUE.   !   Shallow Convection (Control)   Switch
      logical           ::  Oshal             !   Shallow Convection (Effective) Switch
      logical           ::  Los_CP            !   Shallow Convection (user)      Switch
      logical           ::  Odown0 = .TRUE.   !   Flag to take or not convective downdrafts into account
      logical           ::  Orset0 = .TRUE.   !   Flag to refresh or not all tendencies at every call
      logical           ::  OsetA0 = .false.  !   Flag to set convective adjustment time by user
      logical           ::  OCvTC0 = .false.  !   Flag to compute convective transport for  chemical tracer
      logical           ::  Lo_ANA = .false.  !   Flag to compute Anabatic Wind Speed Contribution to Vertical Velocity

      real(kind=real8), SAVE  ::  t_d_CP            !   Deep    Convection Time Scale                                        [s]
      real(kind=real8), SAVE  ::  t_s_CP            !   Shallow Convection Time Scale                                        [s]

      real(kind=real8), SAVE  ::  dtd_CP            !   time interval between 2 CALLs of deep convection (user INPUT)        [s]
      real(kind=real8), SAVE  ::  pdtCV0 = 1200.    !   time interval between 2 CALLs of deep convection (control)           [s]
      real(kind=real8), SAVE  ::  pdtCVx            !   time interval between 2 CALLs of deep convection (chosen)            [s]
      REAL              ::  pdtCV             !   time interval between 2 CALLs of deep convection (effective)         [s]

      integer, SAVE           ::  jjtCV0            !   Number of  Convective Steps for 1    Update Step                     [-]
      integer, SAVE           ::  iitCV0            !   No     of  Convective Step                                           [-]

      real(kind=real8), SAVE  ::  PTdcv0 = 1800.    !   deep    adjustment time                          (control)           [s]
      REAL              ::  PTdcv             !   deep    adjustment time                          (effective)         [s]

      real(kind=real8), SAVE  ::  PTscv0 =  600.    !   shallow adjustment time                          (control)           [s]
      REAL              ::  PTscv             !   shallow adjustment time                          (effective)         [s]

      integer, SAVE           ::  kidia0 =    1     !   index of the first column of the vector                              [-]
      integer, SAVE           ::  kfdia0            !   index of the last  column of the vector                              [-]

      integer, SAVE           ::  kbdia0 =    1     !   vertical computations: lowest                 level                  [-]
      integer, SAVE           ::  ktdia0 =    1     !   vertical computations: over KLEV + 1 - ktdia0 levels                 [-]

      integer, SAVE           ::  kIce_0 =    1     !   flag for ice ( 1 / 0 ) = (included / not icluded)                    [-]
      integer, SAVE           ::  kensbl =    3     !   Nb of ensembles for a "climate" run                                  [-]







      end module Mod_PHY_CP_ctr
