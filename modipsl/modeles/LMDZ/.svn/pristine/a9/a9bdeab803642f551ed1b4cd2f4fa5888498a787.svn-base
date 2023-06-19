SUBROUTINE flumass_loc(massebx,masseby, vcont, ucont, pbaru, pbarv )
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr. Hourdin.
!-------------------------------------------------------------------------------
! Purpose: Compute mass flux at s levels.
  USE parallel_lmdz
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: massebx(ijb_u:ije_u,llm)
  REAL, INTENT(IN)  :: masseby(ijb_v:ije_v,llm)
  REAL, INTENT(IN)  :: vcont  (ijb_v:ije_v,llm)
  REAL, INTENT(IN)  :: ucont  (ijb_u:ije_u,llm)
  REAL, INTENT(OUT) :: pbaru  (ijb_u:ije_u,llm)
  REAL, INTENT(OUT) :: pbarv  (ijb_v:ije_v,llm)
!===============================================================================
! Method used:   A 2 equations system is solved.
!   * 1st one describes divergence computation at pole point nr. i (i=1 to im):
!     (0.5*(pbaru(i)-pbaru(i-1))-pbarv(i))/aire(i) = - SUM(pbarv(n))/aire pole
!   * 2nd one specifies that mean mass flux at pole is equal to 0:
!     SUM(pbaru(n)*local_area(n))=0
! This way, we determine additive constant common to pbary elements representing
!   pbaru(0,j,l) in divergence computation equation for point i=1. (i=1 to im)
!===============================================================================
! Local variables:
  REAL    :: sairen, saireun, ctn, ctn0, apbarun(iim)
  REAL    :: saires, saireus, cts, cts0, apbarus(iim)
  INTEGER :: l, i, ij, ijb, ije
!===============================================================================
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
  DO l=1,llm

    ijb=ij_begin
    ije=ij_end+iip1
    IF(pole_nord) ijb=ij_begin+iip1
    IF(pole_sud)  ije=ij_end-iip1
    pbaru(ijb:ije,l)=massebx(ijb:ije,l)*ucont(ijb:ije,l)

    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF(pole_nord) ijb=ij_begin
    IF(pole_sud)  ije=ij_end-iip1
    pbarv(ijb:ije,l)=masseby(ijb:ije,l)*vcont(ijb:ije,l)

  END DO
!$OMP END DO NOWAIT

  !--- North pole
  IF(pole_nord) THEN
    sairen =SUM(aire (1:iim))
    saireun=SUM(aireu(1:iim))
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
    DO l=1,llm
      ctn=SUM(pbarv(1:iim,l))/sairen
      pbaru(1,l)= pbarv(1,l)-ctn*aire(1)
      DO i=2,iim
        pbaru(i,l)=pbaru(i-1,l)+pbarv(i,l)-ctn*aire(i)
      END DO
      apbarun(:)=aireu(1:iim)*pbaru(1:iim,l)
      ctn0 = -SUM(apbarun)/saireun
      pbaru(1:iim,l)=2.*(pbaru(1:iim,l)+ctn0)
      pbaru(iip1,l)=pbaru(1,l)
    END DO
!$OMP END DO NOWAIT              
  END IF

  !--- South pole
  IF(pole_sud) THEN
    saires =SUM(aire (ip1jm+1:ip1jmp1-1))
    saireus=SUM(aireu(ip1jm+1:ip1jmp1-1))
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
    DO l=1,llm
      cts=SUM(pbarv(1+ip1jmi1:ip1jm-1,l))/saires
      pbaru(1+ip1jm,l)=-pbarv(1+ip1jmi1,l)+cts*aire(1+ip1jm)
      DO i=2,iim
        pbaru(i+ip1jm,l)=pbaru(i-1+ip1jm,l)-pbarv(i+ip1jmi1,l)+cts*aire(i+ip1jm)
      END DO
      apbarus(:)=aireu(1+ip1jm:ip1jmp1-1)*pbaru(1+ip1jm:ip1jmp1-1,l)
      cts0 = -SUM(apbarus)/saireus
      pbaru(1+ip1jm:ip1jmp1-1,l)=2.*(pbaru(1+ip1jm:ip1jmp1-1,l)+cts0)
      pbaru(ip1jmp1,l)=pbaru(1+ip1jm,l)
    END DO
!$OMP END DO NOWAIT         
  END IF

END SUBROUTINE flumass_loc

