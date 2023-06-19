      subroutine PHY_Atm_RT_RUN(rklonr,iklOUT)

!------------------------------------------------------------------------------+
!                                                         Tue  7-May-2013  MAR |
!   MAR          PHY_Atm_RT_RUN                                                |
!     subroutine PHY_Atm_RT_RUN performs   Radiative Vertical Transfer  Scheme |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue  7-May-2013      |
!           Last Modification by H. Gallee,               Tue  7-May-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      IMPLICIT NONE



      integer                      ::  rklonr                           ! nb    des pts de grilles
      integer                      ::  iklOUT                           ! OUTPUT   (pt  de grille)



! Initialization of DATA
! ======================

!                **********
      CALL       Atm_RT_RUN(rklonr,iklOUT)
!                **********



      end subroutine PHY_Atm_RT_RUN
