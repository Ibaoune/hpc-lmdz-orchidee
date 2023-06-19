MODULE YOMGC

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    Grid point boundaries

! RCORI      : Coriolis parameter "f = 2 Omega sin(theta)".
! RCORIC     : 2 Omega cos(theta).
! GEMU       : sine of geographic latitude "sin(theta)".
! GSQM2      : cosine of geographic latitude "cos(theta)".
! GELAM      : geographic longitude "lambda".
! GELAT      : geographic latitude "theta".
! GECLO      : cosine of geographic longitude "cos(lambda)".
! GESLO      : sine of geographic longitude "sin(lambda)".
! GM         : mapping factor "M".
! GOMVRL     : zonal component of vector "2 vec(Omega) wedge a vec(k)".
! GOMVRM     : meridian component of vector "2 vec(Omega) wedge a vec(k)".
! GNORDL     : zonal component "gnordl" of unit vector directed towards
!              the geographic Northern pole.
! GNORDM     : meridian component "gnordm" of unit vector directed towards
!              the geographic Northern pole.
! GNORDLCL   : zonal component of vector "rssavnabla gnordl".
! GNORDMCL   : zonal component of vector "rssavnabla gnordm".
! GNORDMCM   : meridian component of vector "rssavnabla gnordm".
! OROG       : grid-point surface orography "Phi_s".
! OROGL      : zonal component of "vnabla' Phi_s".
! OROGM      : zonal component of "vnabla' Phi_s".
! OROGLL,OROGMM,OROGLM: second-order reduced derivatives of surface orography. 
! GAW        : Gaussian weight.
! NGPLAT     : DM-global number of the Gaussian grid latitude.
! NUNIQUEGP  : pointer array (see computation in sugem2.F90).

REAL(KIND=JPRB),ALLOCATABLE:: RCORI(:)
REAL(KIND=JPRB),ALLOCATABLE:: RCORIC(:)
REAL(KIND=JPRB),ALLOCATABLE:: GEMU(:)
REAL(KIND=JPRB),ALLOCATABLE:: GSQM2(:)
REAL(KIND=JPRB),ALLOCATABLE:: GELAM(:)
REAL(KIND=JPRB),ALLOCATABLE:: GELAT(:)
REAL(KIND=JPRB),ALLOCATABLE:: GECLO(:)
REAL(KIND=JPRB),ALLOCATABLE:: GESLO(:)
REAL(KIND=JPRB),ALLOCATABLE:: GM(:)
REAL(KIND=JPRB),ALLOCATABLE:: GOMVRL(:)
REAL(KIND=JPRB),ALLOCATABLE:: GOMVRM(:)
REAL(KIND=JPRB),ALLOCATABLE:: GNORDL(:)
REAL(KIND=JPRB),ALLOCATABLE:: GNORDM(:)
REAL(KIND=JPRB),ALLOCATABLE:: GNORDLCL(:)
REAL(KIND=JPRB),ALLOCATABLE:: GNORDMCL(:)
REAL(KIND=JPRB),ALLOCATABLE:: GNORDMCM(:)
REAL(KIND=JPRB),ALLOCATABLE:: OROG(:)
REAL(KIND=JPRB),ALLOCATABLE:: OROGL(:)
REAL(KIND=JPRB),ALLOCATABLE:: OROGM(:)
REAL(KIND=JPRB),ALLOCATABLE:: OROGLL(:)
REAL(KIND=JPRB),ALLOCATABLE:: OROGMM(:)
REAL(KIND=JPRB),ALLOCATABLE:: OROGLM(:)
REAL(KIND=JPRB),ALLOCATABLE:: GAW(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NGPLAT(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NUNIQUEGP(:)

!      ----------------------------------------------------------------

!$OMP THREADPRIVATE(gaw,geclo,gelam,gelat,gemu,geslo,gm,gnordl,gnordlcl,gnordm,gnordmcl,gnordmcm,gomvrl)
!$OMP THREADPRIVATE(gomvrm,gsqm2,ngplat,nuniquegp,orog,orogl,orogll,oroglm,orogm,orogmm,rcori,rcoric)
END MODULE YOMGC
