subroutine water_int(klon,klev,field3d,mass,field2d)

!=============================================================
! Compute the 2D burden from 3D mixing ratios
!  OB (obolmd@lmd.jussieu.fr)
!=============================================================

IMPLICIT none

! Arguments
INTEGER, INTENT(IN) :: klon,klev
REAL, DIMENSION(klon,klev),INTENT(IN)  :: field3d, mass 
REAL, DIMENSION(klon),     INTENT(OUT) :: field2d

INTEGER k

field2d(:)=0.0
DO k=1, klev
field2d(:)=field2d(:)+field3d(:,k)*mass(:,k)
ENDDO

RETURN
END
