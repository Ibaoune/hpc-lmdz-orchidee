SUBROUTINE enercin_loc ( vcov, ucov, vcont, ucont, ecin )
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van.
!-------------------------------------------------------------------------------
! Purpose: Compute kinetic energy at sigma levels.
  USE parallel_lmdz
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: vcov    (ijb_v:ije_v,llm)
  REAL, INTENT(IN)  :: ucov    (ijb_u:ije_u,llm)
  REAL, INTENT(IN)  :: vcont   (ijb_v:ije_v,llm)
  REAL, INTENT(IN)  :: ucont   (ijb_u:ije_u,llm)
  REAL, INTENT(OUT) :: ecin    (ijb_u:ije_u,llm)
!===============================================================================
! Notes:
!                 . V
!                i,j-1
!
!      alpha4 .       . alpha1
!
!
!        U .      . P     . U
!       i-1,j    i,j      i,j
!
!      alpha3 .       . alpha2
!
!
!                 . V
!                i,j
!
! Kinetic energy at scalar point P(i,j) (excluding poles) is:
!       Ecin = 0.5 * U(i-1,j)**2 *( alpha3 + alpha4 )  +
!              0.5 * U(i  ,j)**2 *( alpha1 + alpha2 )  +
!              0.5 * V(i,j-1)**2 *( alpha1 + alpha4 )  +
!              0.5 * V(i,  j)**2 *( alpha2 + alpha3 )
!===============================================================================
! Local variables:
  INTEGER :: l, ij, i, ijb, ije
  REAL    :: ecinni(iim), ecinsi(iim), ecinpn, ecinps
!===============================================================================
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm

    ijb=ij_begin
    ije=ij_end+iip1

    IF(pole_nord) ijb=ij_begin+iip1
    IF(pole_sud)  ije=ij_end-iip1

    DO ij = ijb,ije-1
      ecin(ij+1,l)=0.5*(ucov(ij    ,l)*ucont(ij    ,l)*alpha3p4(ij +1)          &
                      + ucov(ij+1  ,l)*ucont(ij+1  ,l)*alpha1p2(ij +1)          &
                      + vcov(ij-iim,l)*vcont(ij-iim,l)*alpha1p4(ij +1)          &
                      + vcov(ij+1  ,l)*vcont(ij+1  ,l)*alpha2p3(ij +1) )
    END DO

    !--- Correction: ecin(1,j,l)= ecin(iip1,j,l)
    DO ij=ijb,ije,iip1; ecin(ij,l) = ecin(ij+iim,l); END DO

    !--- North pole
    IF(pole_nord) THEN
      ecinni(:) = vcov(1:iim,l)*vcont(1:iim,l)*aire(1:iim)
      ecinpn = 0.5*SUM(ecinni)/apoln
      ecin(1:iip1,l)=ecinpn
    END IF

    !--- South pole
    IF(pole_sud) THEN
      DO i=1,iim
        ecinsi(i) = vcov(i+ip1jmi1,l)*vcont(i+ip1jmi1,l)*aire(i+ip1jm)
      END DO
      ecinps = 0.5*SUM(ecinsi)/apols
      ecin(1+ip1jm:ip1jmp1,l)=ecinps
    END IF
  END DO
!$OMP END DO NOWAIT

END SUBROUTINE enercin_loc

