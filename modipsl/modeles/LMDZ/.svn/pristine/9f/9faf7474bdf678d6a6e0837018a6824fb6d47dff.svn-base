      subroutine PHY_MAR_DAT(sh___HOST                                 &
     &                      ,sh_a_HOST                                 &
     &                      ,dx___HOST                                 &
     &                      ,dy___HOST                                 &
     &                      ,slopxHOST                                 &
     &                      ,slopyHOST                                 &
     &                      ,slopeHOST                                 &
     &                      ,MMaskHOST                                 &
     &                      ,ixq1,mxqq                                 &
     &                      ,jyq1,myqq                                 &
     &                      ,m_azim)

!------------------------------------------------------------------------------+
!                                                         Sat 15-Jun-2013  MAR |
!     subroutine PHY_MAR_DAT INPUTs PHY_MAR data and                           |
!                            computes horizontal dependencies                  |
!                                                                              |
!     Applied to: ...                                                          |
!                                                                              |
!                                                                              |
! # OPTIONS: #xy  ...                                                          |
! # ^^^^^^^^ #..  ...                                                          |
!                                                                              |
!                                                                              |
! # CAUTION: Highly simplified version of PHY_MAR_DAT for academic applications|
! # ^^^^^^^^                                                                   |
!                                                                              |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Fri 14-Jun-2013      |
!           Last Modification by H. Gallee,               Sat 15-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real


      IMPLICIT NONE


      real            , dimension(ixq1:mxqq,jyq1:myqq)         ::  sh___HOST
      real(kind=real8), dimension(ixq1:mxqq,jyq1:myqq)         ::  sh_a_HOST
      real(kind=real8), dimension(ixq1:mxqq,jyq1:myqq)         ::  dx___HOST
      real(kind=real8), dimension(ixq1:mxqq,jyq1:myqq)         ::  dy___HOST
      real(kind=real8), dimension(ixq1:mxqq,jyq1:myqq)         ::  slopxHOST
      real(kind=real8), dimension(ixq1:mxqq,jyq1:myqq)         ::  slopyHOST
      real(kind=real8), dimension(ixq1:mxqq,jyq1:myqq)         ::  slopeHOST
      real(kind=real8), dimension(ixq1:mxqq,jyq1:myqq,m_azim)  ::  MMaskHOST

      integer                                                  ::  ixq1,mxqq
      integer                                                  ::  jyq1,myqq
      integer                                                  ::  m_azim




! Local Variables
! ===============

      integer                                                  ::  i   ,j   ,k




! Topographie
! ===========

        DO j= jyq1,myqq
        DO i= ixq1,mxqq
!         sh___HOST(i,j)   = 0.00
          sh_a_HOST(i,j)   = 0.00




! Discretisation Horizontale
! ==========================

          dx___HOST(i,j)   = 1.e3
          dy___HOST(i,j)   = 1.e3




! Pente
! =====

          slopxHOST(i,j)   = 0.00
          slopyHOST(i,j)   = 0.00
          slopeHOST(i,j)   = 0.00




! Masque de Montagnes
! ===================

        DO k=    1,m_azim
          MMaskHOST(i,j,k) = 0.00
        ENDDO


        ENDDO
        ENDDO




      end subroutine PHY_MAR_DAT
