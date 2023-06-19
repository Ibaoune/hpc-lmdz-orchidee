MODULE YOMPRAD

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!----------------------------------------------------------------
!*    new radiation grid data structures
!----------------------------------------------------------------

! NRESOL_ID : transformation package resolution identifier
! NSMAX     : spectral truncation
! NASM0     : address in a spectral array of (m, n=m)
! MYMS      : actual wave numbers handled by this processor
! NDGLG     : number of latitude rows
! NDLON     : length of a row of latitude near equator
! NDGSAL    : Local version of NDGSA , always 1
! NDGENL    : Number of latitude rows for which this process has grid points
! NDGSAH    : 1-NROWIDE
! NDGENH    : NDGENL+NROWIDE
! NDGSAG    : Global version of NDGSA
! NDGENG    : Global version of NDGEN
! NDSUR1    : over dimensioning of NDLON for technical reasons (at least 2)
! NDLSUR    : NDLON+NDSUR1
! NDGSUR    : number of additional rows at each pole for semi-lagrangian
!           : interpolation
! MYFRSTACTLAT : first actual lat on this PE in grid-point space,
!              : it is nfrstlat(my_region_ns)
! MYLSTACTLAT  : last actual lat on this PE in grid-point space,
!              : it is nlstlat(my_region_ns)
! NGPTOT    : number of local grid points
! NGPTOTG   : total number of grid points
! NRGRI     : number of grid points on a latitude row
! NLOENG    : as above but extended at poles for SL interpolation
! GELAM     : radiation grid geographic longitude "lambda"
! GELAT     : radiation grid geographic latitude "theta"
! GESLO     : sine of geographic longitude "sin(lambda)"
! GECLO     : cosine of geographic longitude "cos(lambda)
! GEMU      : sine of geographic latitude "sin(theta)"
! RMU       : mu              sin(theta)
! RSQM2     : SQRT(R1MU2)     cos(theta)
! RLATIG    : arcsin(mu)      theta  GLOBAL VIEW
! RLATI     : arcsin(mu)      theta
! RIPI0     : bi-cubic interpolation coefficients
! RIPI1     : bi-cubic interpolation coefficients
! RIPI2     : bi-cubic interpolation coefficients

TYPE RADIATION_GRID_STRUCT
INTEGER(KIND=JPIM)                            :: NRESOL_ID, NGPTOT, NGPTOTG, &
 & NGPTOTMX, NSPEC2, NSMAX, &
 & NPTRFLOFF, NUMP, NDLON, &
 & NDGSAL, NDGENL, NDGSAH, NDGENH, &
 & NDGLG, NDGSAG, NDGENG, NDLSUR, &
 & NFRSTLOFF, NDSUR1,NDGSUR, &
 & MYFRSTACTLAT, MYLSTACTLAT  
INTEGER(KIND=JPIM), POINTER, DIMENSION(:)     :: NRGRI, NLOENG, NPTRFRSTLAT, NFRSTLAT, &
 & NLSTLAT, MYMS, NASM0  
INTEGER(KIND=JPIM), POINTER, DIMENSION(:,:)   :: NSTA, NONL
REAL(KIND=JPRB),    POINTER, DIMENSION(:)     :: GELAM, GELAT, GECLO, GESLO, GEMU, &
 & RMU, RSQM2, RLATIG, RLATI, &
 & RIPI0, RIPI1, RIPI2  
END TYPE RADIATION_GRID_STRUCT

TYPE(RADIATION_GRID_STRUCT) :: RADGRID

!----------------------------------------------------------------
!*    radiation on demand comms data structures
!----------------------------------------------------------------

LOGICAL :: LRADONDEM
LOGICAL :: LRADONDEM_ACTIVE
INTEGER(KIND=JPIM) :: NFIXRADFLD(2)
INTEGER(KIND=JPIM),ALLOCATABLE :: MASK_RI1(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MASK_RI2(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MASK_RO1(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MASK_RO2(:)

!----------------------------------------------------------------
!*    radiation interpolation/load balancing data structures
!----------------------------------------------------------------

! NRIRINT  : INTERPOLATION INTERVAL FOR LOCAL LATITUDES
! NRFRSTOFF: OFFSET TO FIRST COURSE POINT FROM FIRST FINE POINT
! NRLASTOFF: OFFSET TO FIRST COURSE POINT FROM LAST  FINE POINT
! NRIMAX   : NUMBER OF COURSE POINTS FOR LOCAL LATITUDES
! NRIMAXT  : TOTAL NUMBER OF COURSE POINTS IN THIS PROCESSOR
! NRIMAXLT : TOTAL NUMBER OF COURSE POINTS TO BE PROCESSED BY THIS PROC
! NRLPRCS  : NUMBER OF PROCESSORS REQUIRED FOR LOAD BALANCING
! NRIMAXLA : LOAD BALANCING DATA - SETA OF PROCESSOR TO COMMUNICATE WITH
! NRIMAXLB : LOAD BALANCING DATA - SETB OF PROCESSOR TO COMMUNICATE WITH
! NRIMAXLN : LOAD BALANCING DATA - NUMBER OF COURSE POINTS TO SND/RCV
! NRCNEEDW : NUMBER OF COURSE POINTS NEEDED ON WESTERN END OF LOCAL LATS
! NRCNEEDE : NUMBER OF COURSE POINTS NEEDED ON EASTERN END OF LOCAL LATS
! NRCSNDW  : NUMBER OF WESTERN COURSE POINTS TO SEND TO LOCAL PROCS
! NRCSNDE  : NUMBER OF EASTERN COURSE POINTS TO SEND TO LOCAL PROCS
! NRCRCVW  : NUMBER OF WESTERN COURSE POINTS TO RCV FROM LOCAL PROCS
! NRCRCVE  : NUMBER OF WESTERN COURSE POINTS TO RCV FROM LOCAL PROCS
! NRCSNDT  : TOTAL NUMBER OF COURSE POINTS TO SEND TO LOCAL PROCESSORS
! NRCRCVT  : TOTAL NUMBER OF COURSE POINTS TO RCV FROM LOCAL PROCESSORS
! NRCRCVWO : OFFSET OF WESTERN BOUNDARY COURSE POINTS IN BOUNDARY ARRAY
! NRCRCVEO : OFFSET OF EASTERN BOUNDARY COURSE POINTS IN BOUNDARY ARRAY

! LODBGRADI: DEBUGGING FLAG, EACH PROCESSOR WILL DUMP ABOVE INFO TO A
!            PER PROCESSOR FILE CALLED 'debug_aaa_bbb' 
!            where aaa=my_region_ns and bbb=my_region_ew
! LODBGRADL: DEBUGGING FLAG, EACH PROCESSOR WILL DUMP LOAD BALANCING
!            ARRAYS as above

! NRLBCHUNKS: NUMBER OF BLOCKING LOOPS FOR LOAD-BALANCING COMMS
! NRLRCHUNKS: NUMBER OF BLOCKING LOOPS FOR LOAD-BALANCING COMMS (RESTORING)
! NRLBPOINTS: MAX NUMBER OF POINTS FOR LOAD-BALANCING COMMS
! NRLRPOINTS: MAX NUMBER OF POINTS FOR LOAD-BALANCING COMMS (RESTORING)
! NRLBDATA  : NUMBER OF DATA VALUES PER LOAD-BALANCING POINT 
! NRLRDATA  : NUMBER OF DATA VALUES PER LOAD-BALANCING POINT (RESTORING)

LOGICAL :: LODBGRADI
LOGICAL :: LODBGRADL
INTEGER(KIND=JPIM) :: NRLBCHUNKS
INTEGER(KIND=JPIM) :: NRLRCHUNKS
INTEGER(KIND=JPIM) :: NRLBPOINTS
INTEGER(KIND=JPIM) :: NRLRPOINTS
INTEGER(KIND=JPIM) :: NRLBDATA
INTEGER(KIND=JPIM) :: NRLRDATA
INTEGER(KIND=JPIM), PARAMETER :: JPMAXLOADP=32
INTEGER(KIND=JPIM),ALLOCATABLE:: NRIRINT(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRFRSTOFF(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRLASTOFF(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRIMAX(:,:)
INTEGER(KIND=JPIM) :: NRIMAXT
INTEGER(KIND=JPIM) :: NRIMAXLT
INTEGER(KIND=JPIM) :: NRLPRCS
INTEGER(KIND=JPIM) :: NRIMAXLA(JPMAXLOADP)
INTEGER(KIND=JPIM) :: NRIMAXLB(JPMAXLOADP)
INTEGER(KIND=JPIM) :: NRIMAXLN(JPMAXLOADP)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRCNEEDW(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRCNEEDE(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRCSNDW(:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRCSNDE(:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRCRCVW(:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRCRCVE(:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRCSNDT(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRCRCVT(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRCRCVWO(:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRCRCVEO(:,:,:)
!$OMP THREADPRIVATE(lodbgradi,lodbgradl,lradondem,lradondem_active,nfixradfld,nrimaxla,nrimaxlb)
!$OMP THREADPRIVATE(nrimaxln,nrimaxlt,nrimaxt,nrlbchunks,nrlbdata,nrlbpoints,nrlprcs,nrlrchunks)
!$OMP THREADPRIVATE(nrlrdata,nrlrpoints,radgrid)
!$OMP THREADPRIVATE(mask_ri1,mask_ri2,mask_ro1,mask_ro2,nrcneede,nrcneedw,nrcrcve,nrcrcveo,nrcrcvt)
!$OMP THREADPRIVATE(nrcrcvw,nrcrcvwo,nrcsnde,nrcsndt,nrcsndw,nrfrstoff,nrimax,nrirint,nrlastoff)
END MODULE YOMPRAD
