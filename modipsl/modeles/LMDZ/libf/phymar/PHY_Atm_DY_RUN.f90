      subroutine PHY_Atm_DY_RUN

!------------------------------------------------------------------------------+
!                                                         Sat  8-Jun-2013  MAR |
!   MAR          PHY_Atm_DY_RUN                                                |
!     subroutine PHY_Atm_DY_RUN  updates Variables coming from MAR Dynamics    |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Fri 22-Mar-2013      |
!           Last Modification by H. Gallee,               Sat  8-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_DY_kkl
      use Mod_PHY_AT_kkl
! #LD use Mod_PHY_CM_kkl


      IMPLICIT NONE



      integer    ::    i,     j,     k,     ikl



! Assignation    of Mod_PHY_DY_kkl
! ================================


        DO ikl = 1,kcolp

        DO   k = 1,mzpp  

! #LD     ld_H2O(ikl,k) = qv__DY(ikl,k)                                 &
! #LD&           - 1.64 *(qw__CM(ikl,k)+qi__CM(ikl,k)                   &
! #LD&                   +qs__CM(ikl,k)+qr__CM(ikl,k))                  &
!    &                   +Aerosols ...
! #LD&           + 0.00

          TmidDY(ikl,k) =(Ta__DY(ikl,k)+Ta__DY(ikl,min(k+1,mzpp))) * 0.5
        ENDDO


        ENDDO


      end subroutine PHY_Atm_DY_RUN
