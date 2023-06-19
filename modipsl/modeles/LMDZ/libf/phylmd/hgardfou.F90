
! $Id: hgardfou.F90 2399 2015-11-20 16:23:28Z emillour $
SUBROUTINE hgardfou(t, tsol, text,abortphy)
  USE dimphy, ONLY: klon, klev
  USE phys_state_var_mod, ONLY: pctsrf
  USE geometry_mod, ONLY: longitude_deg, latitude_deg
  USE indice_sol_mod, ONLY: nbsrf
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE
  ! ======================================================================
  ! Verifier la temperature
  ! ======================================================================
  include "YOMCST.h"
  REAL t(klon, klev), tsol(klon, nbsrf)
  CHARACTER(len=*), intent(in):: text
  CHARACTER (LEN=20) :: modname = 'hgardfou'
  INTEGER abortphy

  INTEGER i, k, nsrf
  REAL zt(klon)
  INTEGER jadrs(klon), jbad
  LOGICAL ok

  LOGICAL firstcall
  SAVE firstcall
  DATA firstcall/.TRUE./
  !$OMP THREADPRIVATE(firstcall)

  IF (firstcall) THEN
    WRITE (lunout, *) 'hgardfou garantit la temperature dans [100,370] K'
    firstcall = .FALSE.
    ! DO i = 1, klon
    ! WRITE(lunout,*)'i=',i,'rlon=',rlon(i),'rlat=',rlat(i)
    ! ENDDO

  END IF

  ok = .TRUE.
  DO k = 1, klev
    DO i = 1, klon
      zt(i) = t(i, k)
    END DO
#ifdef CRAY
    CALL whenfgt(klon, zt, 1, 370.0, jadrs, jbad)
#else
    jbad = 0
    DO i = 1, klon
      IF (zt(i)>370.) THEN
        jbad = jbad + 1
        jadrs(jbad) = i
      END IF
    END DO
#endif
    IF (jbad>0) THEN
      ok = .FALSE.
      DO i = 1, jbad
        WRITE (lunout, *) 'i,k,temperature,lon,lat,pourc ter,lic,oce,sic =', &
          jadrs(i), k, zt(jadrs(i)), longitude_deg(jadrs(i)), &
          latitude_deg(jadrs(i)),(pctsrf(jadrs(i),nsrf), nsrf=1, nbsrf)
      END DO
    END IF
#ifdef CRAY
    CALL whenflt(klon, zt, 1, 100.0, jadrs, jbad)
#else
    jbad = 0
    DO i = 1, klon
      ! IF (zt(i).LT.100.0) THEN
      IF (zt(i)<50.0) THEN
        jbad = jbad + 1
        jadrs(jbad) = i
      END IF
    END DO
#endif
    IF (jbad>0) THEN
      ok = .FALSE.
      DO i = 1, jbad
        WRITE (lunout, *) 'i,k,temperature,lon,lat,pourc ter,lic,oce,sic =', &
          jadrs(i), k, zt(jadrs(i)), longitude_deg(jadrs(i)), &
          latitude_deg(jadrs(i)), (pctsrf(jadrs(i),nsrf), nsrf=1, nbsrf)
      END DO
    END IF
  END DO

  DO nsrf = 1, nbsrf
    DO i = 1, klon
      zt(i) = tsol(i, nsrf)
    END DO
#ifdef CRAY
    CALL whenfgt(klon, zt, 1, 370.0, jadrs, jbad)
#else
    jbad = 0
    DO i = 1, klon
      IF (zt(i)>370.0) THEN
        jbad = jbad + 1
        jadrs(jbad) = i
      END IF
    END DO
#endif
    IF (jbad>0) THEN
      ok = .FALSE.
      DO i = 1, jbad
        WRITE (lunout, *) &
          'i,nsrf,temperature,lon,lat,pourc ter,lic,oce,sic =', jadrs(i), &
          nsrf, zt(jadrs(i)), longitude_deg(jadrs(i)), &
          latitude_deg(jadrs(i)), pctsrf(jadrs(i), nsrf)
      END DO
    END IF
#ifdef CRAY
    CALL whenflt(klon, zt, 1, 100.0, jadrs, jbad)
#else
    jbad = 0
    DO i = 1, klon
      ! IF (zt(i).LT.100.0) THEN
      IF (zt(i)<50.0) THEN
        jbad = jbad + 1
        jadrs(jbad) = i
      END IF
    END DO
#endif
    IF (jbad>0) THEN
      ok = .FALSE.
      DO i = 1, jbad
        WRITE (lunout, *) &
          'i,nsrf,temperature,lon,lat,pourc ter,lic,oce,sic =', jadrs(i), &
          nsrf, zt(jadrs(i)), longitude_deg(jadrs(i)), &
          latitude_deg(jadrs(i)), pctsrf(jadrs(i), nsrf)
      END DO
    END IF
  END DO

!  IF (.NOT. ok) CALL abort_physic(modname, text, 1)
  IF (.NOT. ok) abortphy=1

END SUBROUTINE hgardfou
