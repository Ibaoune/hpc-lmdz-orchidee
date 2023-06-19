SUBROUTINE borne_var_surf(klon,klev,nbsrf,                   &
         iflag_bug_t2m_stab_ipslcm61,                        &
         t1,q1,u1,v1,                                        &
         ftsol, qsurf, pctsrf, paprs,                        &
         t2m, q2m, u10m, v10m,                               &
         zt2m_cor, zq2m_cor, zu10m_cor, zv10m_cor,           &
         zrh2m_cor, zqsat2m_cor)

IMPLICIT NONE

!==================================================================
! Declarations
!==================================================================

! arguments
INTEGER klon,klev,nbsrf,iflag_bug_t2m_stab_ipslcm61
REAL,DIMENSION(klon),INTENT(IN) :: t1, q1, u1, v1
REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: t2m, q2m, u10m, v10m
REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: ftsol, pctsrf
REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs 
REAL,DIMENSION(klon),INTENT(IN) :: qsurf
REAL,DIMENSION (klon),INTENT(OUT) :: zt2m_cor, zq2m_cor, zu10m_cor, zv10m_cor
REAL,DIMENSION (klon),INTENT(OUT)  :: zrh2m_cor, zqsat2m_cor


! local
INTEGER i,nsrf
REAL,DIMENSION (klon,nbsrf) :: t2m_cor, q2m_cor, u10m_cor, v10m_cor
REAL                               :: zx_qs1, zcor1, zdelta1
include "YOMCST.h"
include "YOETHF.h"
include "FCTTRE.h"
!==================================================================
! Correction of sub surface variables
!==================================================================

zrh2m_cor=0.
zqsat2m_cor=0.

DO nsrf=1,nbsrf
   DO i=1,klon
    t2m_cor(i,nsrf)=t2m(i,nsrf)
    q2m_cor(i,nsrf)=q2m(i,nsrf)
    u10m_cor(i,nsrf)=u10m(i,nsrf)
    v10m_cor(i,nsrf)=v10m(i,nsrf)
     IF(iflag_bug_t2m_stab_ipslcm61.EQ.-2.AND.q2m(i,nsrf).LT.0.) THEN
      t2m_cor(i,nsrf)=MIN(t2m(i,nsrf),MAX(t1(i),ftsol(i,nsrf)))
      t2m_cor(i,nsrf)=MAX(t2m_cor(i,nsrf),MIN(t1(i),ftsol(i,nsrf)))
      q2m_cor(i,nsrf)=MIN(q2m(i,nsrf),MAX(q1(i),qsurf(i)))
      q2m_cor(i,nsrf)=MAX(q2m_cor(i,nsrf),MIN(q1(i),qsurf(i)))
      q2m_cor(i,nsrf)=MAX(q2m_cor(i,nsrf),0.)
      u10m_cor(i,nsrf)=SIGN(MIN(ABS(u1(i)),ABS(u10m(i,nsrf))),u1(i))
      v10m_cor(i,nsrf)=SIGN(MIN(ABS(v1(i)),ABS(v10m(i,nsrf))),v1(i))
     ELSEIF(iflag_bug_t2m_stab_ipslcm61.EQ.-1.AND.(ftsol(i,nsrf).LE.t1(i).OR.q2m(i,nsrf).LT.0.)) THEN
      t2m_cor(i,nsrf)=MIN(t2m(i,nsrf),MAX(t1(i),ftsol(i,nsrf)))
      t2m_cor(i,nsrf)=MAX(t2m_cor(i,nsrf),MIN(t1(i),ftsol(i,nsrf)))
      q2m_cor(i,nsrf)=MIN(q2m(i,nsrf),MAX(q1(i),qsurf(i)))
      q2m_cor(i,nsrf)=MAX(q2m_cor(i,nsrf),MIN(q1(i),qsurf(i)))
      q2m_cor(i,nsrf)=MAX(q2m_cor(i,nsrf),0.)
      u10m_cor(i,nsrf)=SIGN(MIN(ABS(u1(i)),ABS(u10m(i,nsrf))),u1(i))
      v10m_cor(i,nsrf)=SIGN(MIN(ABS(v1(i)),ABS(v10m(i,nsrf))),v1(i))
     ELSEIF(iflag_bug_t2m_stab_ipslcm61.EQ.1.AND.ftsol(i,nsrf).LE.t1(i)) THEN
      t2m_cor(i,nsrf)=MIN(t2m(i,nsrf),MAX(t1(i),ftsol(i,nsrf)))
      t2m_cor(i,nsrf)=MAX(t2m_cor(i,nsrf),MIN(t1(i),ftsol(i,nsrf)))
      q2m_cor(i,nsrf)=MIN(q2m(i,nsrf),MAX(q1(i),qsurf(i)))
      q2m_cor(i,nsrf)=MAX(q2m_cor(i,nsrf),MIN(q1(i),qsurf(i)))
      q2m_cor(i,nsrf)=MAX(q2m_cor(i,nsrf),0.)
      u10m_cor(i,nsrf)=SIGN(MIN(ABS(u1(i)),ABS(u10m(i,nsrf))),u1(i))
      v10m_cor(i,nsrf)=SIGN(MIN(ABS(v1(i)),ABS(v10m(i,nsrf))),v1(i))
     ELSEIF(iflag_bug_t2m_stab_ipslcm61.EQ.0) THEN
      t2m_cor(i,nsrf)=MIN(t2m(i,nsrf),MAX(t1(i),ftsol(i,nsrf)))
      t2m_cor(i,nsrf)=MAX(t2m_cor(i,nsrf),MIN(t1(i),ftsol(i,nsrf)))
      q2m_cor(i,nsrf)=MIN(q2m(i,nsrf),MAX(q1(i),qsurf(i)))
      q2m_cor(i,nsrf)=MAX(q2m_cor(i,nsrf),MIN(q1(i),qsurf(i)))
      q2m_cor(i,nsrf)=MAX(q2m_cor(i,nsrf),0.)
      u10m_cor(i,nsrf)=SIGN(MIN(ABS(u1(i)),ABS(u10m(i,nsrf))),u1(i))
      v10m_cor(i,nsrf)=SIGN(MIN(ABS(v1(i)),ABS(v10m(i,nsrf))),v1(i))
     ENDIF
!!!
     zdelta1 = MAX(0.,SIGN(1., rtt-t2m_cor(i,nsrf) ))
     zx_qs1  = r2es * FOEEW(t2m_cor(i,nsrf),zdelta1)/paprs(i,1)
     zx_qs1  = MIN(0.5,zx_qs1)
     zcor1   = 1./(1.-RETV*zx_qs1)
     zx_qs1  = zx_qs1*zcor1
     zrh2m_cor(i) = zrh2m_cor(i) + q2m_cor(i,nsrf)/zx_qs1 * pctsrf(i,nsrf)
     zqsat2m_cor(i) = zqsat2m_cor(i) + zx_qs1  * pctsrf(i,nsrf)
!!!
   ENDDO
ENDDO

!==================================================================
! Agregation of sub surfaces
!==================================================================

zt2m_cor=0.
zq2m_cor=0.
zu10m_cor=0.
zv10m_cor=0.
DO nsrf = 1, nbsrf
   DO i = 1, klon
      zt2m_cor(i)  = zt2m_cor(i)  + t2m_cor(i,nsrf)  * pctsrf(i,nsrf)
      zq2m_cor(i)  = zq2m_cor(i)  + q2m_cor(i,nsrf)  * pctsrf(i,nsrf)
      zu10m_cor(i) = zu10m_cor(i) + u10m_cor(i,nsrf) * pctsrf(i,nsrf)
      zv10m_cor(i) = zv10m_cor(i) + v10m_cor(i,nsrf) * pctsrf(i,nsrf)
   ENDDO
ENDDO

RETURN
END


