!
! $Id $
!
      SUBROUTINE gr_fi_ecrit(nfield,nlon,iim,jjmp1,fi,ecrit)
      IMPLICIT NONE
c
c Transfer a field from the full physics grid
c to a full 2D (lonxlat) grid suitable for outputs
c
      INTEGER,INTENT(IN) :: nfield,nlon,iim,jjmp1
      REAL,INTENT(IN) :: fi(nlon,nfield)
      REAL,INTENT(out) :: ecrit(iim*jjmp1,nfield)
c
      INTEGER :: i, n, ig, jjm
c
      jjm = jjmp1 - 1
      DO n = 1, nfield
         DO i=1,iim
            ecrit(i,n) = fi(1,n)
            ecrit(i+jjm*iim,n) = fi(nlon,n)
         ENDDO
         DO ig = 1, nlon - 2
           ecrit(iim+ig,n) = fi(1+ig,n)
         ENDDO
      ENDDO

      END
