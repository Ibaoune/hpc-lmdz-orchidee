MODULE YOEAERSRC

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAERSRC* - CONTROL OPTIONS FOR AEROSOLS' SOURCE
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: NMAXTAER
INTEGER(KIND=JPIM) :: NTAER
INTEGER(KIND=JPIM) :: NBINAER(9)
INTEGER(KIND=JPIM) :: NINDAER(15), JKTYP(15), JKBIN(15)
INTEGER(KIND=JPIM) :: NTYPAER(9)
INTEGER(KIND=JPIM) :: NDDUST

LOGICAL :: LEPAERO, LAEREXTR

REAL(KIND=JPRB) :: RGELAV, RGEMUV, RDGLAV, RDGMUV
REAL(KIND=JPRB) :: RCLONV, RSLONV, RDCLONV, RDSLONV
REAL(KIND=JPRB) :: RLATVOL,RLONVOL

REAL(KIND=JPRB) :: RSSFLX(3)

!     ------------------------------------------------------------------
! MMAXTAER   : MAXIMUM TOTAL NUMBER OF AEROSOLS
! NTAER      : TOTAL NUMBER OF AEROSOLS
! NTYPAER( ) : NBINAER( )
!        (1) :         3 FOR SEA-SALT 
!        (2) :         3 FOR DESERT DUST
!        (3) :         2 FOR ORGANIC MATTERS
!        (4) :         2 FOR BLACK CARBON
!        (5) :         1 FOR SULFATE
!        (6) :         1 FOR FLY ASH
!        (7) :         1 FOR PSEUDO-PROGNOSTIC STRATOSPHERIC AEROSOLS
!        (8) :         1 FOR PROGNOSTIC STRATOSPHERIC AEROSOLS
!        (9) :         1 FOR VOLCANIC AEROSOLS

! RSSFLX     : sea salt flux for 3-size bins (in mg m-2 s-1) for a 1 m s-1
!              wind speed, at 10 m height and 80% RH 
! Rlat/lonVOL: LAT/LON of a possible volcanic eruption
! JDDUST     : 1 =LSCE, 2 =based on MODIS
!     ------------------------------------------------------------------



!$OMP THREADPRIVATE(jkbin,jktyp,laerextr,lepaero,nbinaer,nddust,nindaer)
!$OMP THREADPRIVATE(nmaxtaer,ntaer,ntypaer,rclonv,rdclonv,rdglav,rdgmuv)
!$OMP THREADPRIVATE(rdslonv,rgelav,rgemuv,rlatvol,rlonvol,rslonv,rssflx)

END MODULE YOEAERSRC

