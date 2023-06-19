MODULE YOE_UVRAD

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOE_UVRAD* - COEFFICIENTS FOR ULTRAVIOLET RADIATION PROCESSOR
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: NUV, NRADUV, NUVTIM
INTEGER(KIND=JPIM) :: IPUV(3000), JCOP(3000), JUVLAM(3000)
LOGICAL :: LUVPROC, LUVTDEP, LUVDBG

REAL(KIND=JPRB) :: RK250(3000), RTUV1(3000), RTUV2(3000), RSUVB0(3000), RAYUVB(3000)
REAL(KIND=JPRB) :: RASA(4), RASB(4), RASC(4), RASD(4), RASE(4), RASF(4)
REAL(KIND=JPRB) :: RFA0(4), RFA1(4), RFB0(4), RFB1(4), RFB2(4), RFB3(4), &
 & RFC0(4), RFC1(4), RFC2(4), RFC3(4), RFD0(4), RFD1(4), RFD2(4), RFD3(4)
REAL(KIND=JPRB) :: RTAUVA(4,6), RPIUVA(4,6), RCGUVA(4,6)
REAL(KIND=JPRB) :: RXPO(3), RXPL(3), RCIEAS(3000), RSUVB(3000), RUVLAM(3000)
REAL(KIND=JPRB) :: RFCAER, RFCOZO, RMUZUV
!     -----------------------------------------------------------------
!$OMP THREADPRIVATE(ipuv,jcop,juvlam,luvdbg,luvproc,luvtdep,nraduv,nuv,nuvtim,rasa,rasb)
!$OMP THREADPRIVATE(rasc,rasd,rase,rasf,rayuvb,rcguva,rcieas,rfa0,rfa1,rfb0,rfb1,rfb2,rfb3)
!$OMP THREADPRIVATE(rfc0,rfc1,rfc2,rfc3,rfcaer,rfcozo,rfd0,rfd1,rfd2,rfd3,rk250,rmuzuv,rpiuva)
!$OMP THREADPRIVATE(rsuvb,rsuvb0,rtauva,rtuv1,rtuv2,ruvlam,rxpl,rxpo)
END MODULE YOE_UVRAD

