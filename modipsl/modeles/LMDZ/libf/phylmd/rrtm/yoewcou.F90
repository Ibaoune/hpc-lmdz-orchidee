MODULE YOEWCOU

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!*    ** *YOEWCOU* - VARIABLES FOR COUPLING WITH THE WAVE MODEL

!     P. VITERBO     E.C.M.W.F.       07/10/88
!     P. VITERBO     E.C.M.W.F.       03/02/92
!     J. DOYLE       E.C.M.W.F.       21/11/96 
!     J. BIDLOT      E.C.M.W.F.       13/06/97 
!     J. BIDLOT      E.C.M.W.F.       11/08/98 
!     G. MOZDZYNSKI  E.C.M.W.F.       14/01/05

!      NAME      TYPE       PURPOSE
!      ----      ----       -------

!     *NLONW*    INTEGER    NUMBER OF POINTS IN A LATITUDE LINE IN
!                           THE WAVE MODEL
!     *NLATW*    INTEGER    NUMBER OF LATITUDES IN THE WAVE MODEL.
!     *NLON1W*   INTEGER    *NLONW*
!     *NLAT1W*   INTEGER    TOTAL NUMBER OF LATITUDES WITH THE WAVE
!                           MODEL RESOLUTION.
!     *NNORXW*   INTEGER    NUMBER OF EXTRA POINTS NORTHWARDS OF THE
!                           NORTHERN BOUNDARY OF THE WAVE MODEL.
!     *CBEGDAT*  CHARACTER  INITIAL DATE OF FORECAST (YYYYMMDDHHmm)
!     *NDURAT*   INTEGER    DURATION (IN MINUTES) OF TOTAL WAVE INTEGRATION
!                              (OR FORECAST TIME IN MINUTES).
!     *NSTPW*    INTEGER    FREQUENCY OF CALL TO THE WAVE MODEL.
!     *NRESUM*   INTEGER    TIME STEP OF A RESTART EVENT
!     *LWCOU*    LOGICAL    TRUE IF THE WAVE MODEL IS TO BE RUN.
!     *LWCOU2W*  LOGICAL    TRUE IF TWO-WAY INTERACTION WITH THE WAVE MODEL.
!                           FALSE IF ONE-WAY INTERACTION WITH THE WAVE MODEL
!     *LWCOUNORMS* LOGICAL  TRUE IF NORMS OF COUPLED FIELDS ARE REQUIRED
!     *RSOUTW*   REAL       SOUTH BOUNDARY OF THE WAVE MODEL.
!     *RNORTW*   REAL       NORTH BOUNDARY OF THE WAVE MODEL.
!     *RDEGREW*  REAL       RESOLUTION OF THE WAVE MODEL (DEGREES).
!
!     *MASK_WAVE_IN*  INTEGER  COMMS MASK FOR INPUT TO WAVE MODEL
!     *MASK_WAVE_OUT* INTEGER  COMMS MASK FOR OUTPUT FROM WAVE MODEL
!
!     *LWVIN_MASK_NOT_SET* LOGICAL indicates whether mask_wave_in
!                           has been updated on the first call to the 
!                           wave model
!     *LWVIN_UNINITIALISED* LOGICAL indicates whether the mwvin_* data
!                           structures are initialised
!     *MWVIN_SENDCNT* INTEGER nproc sized array describing how many grid 
!                           points need to be sent (or copied) by the 
!                           local task to remote tasks (or this task)
!     *MWVIN_RECVCNT* INTEGER nproc sized array describing how many grid 
!                           points need to be received by the local task 
!                           from a remote task
!     *MWVIN_SENDTOT* INTEGER total number of grid points to be sent to 
!                           remote tasks
!     *MWVIN_RECVTOT* INTEGER total number of grid points to be received 
!                           from remote tasks
!     *MWVIN_SENDOFF* INTEGER nproc sized array containing offsets into 
!                           the MWVIN_SENDBUF and MWVIN_SENDIND arrays
!     *MWVIN_SENDBUF* INTEGER local indexes of data on remote tasks that
!                           the local task needs
!     *MWVIN_SENDIND* INTEGER global indexes of data on remote tasks that
!                           the local task needs
!
!     *MWVIN_RECVOFF* INTEGER nproc sized array containing offsets into 
!                           the MWVIN_RECVBUF array
!     *MWVIN_RECVBUF* INTEGER local indexes of data on the local task that
!                           remote tasks need
!


INTEGER(KIND=JPIM) :: NLONW
INTEGER(KIND=JPIM) :: NLATW
INTEGER(KIND=JPIM) :: NLON1W
INTEGER(KIND=JPIM) :: NLAT1W
INTEGER(KIND=JPIM) :: NNORXW
INTEGER(KIND=JPIM) :: NDURAT
INTEGER(KIND=JPIM) :: NSTPW
INTEGER(KIND=JPIM) :: NRESUM

REAL(KIND=JPRB) :: RSOUTW
REAL(KIND=JPRB) :: RNORTW
REAL(KIND=JPRB) :: RDEGREW

LOGICAL :: LWCOU
LOGICAL :: LWCOU2W
LOGICAL :: LWCOUNORMS
CHARACTER :: CBEGDAT*12

INTEGER(KIND=JPIM), ALLOCATABLE :: MASK_WAVE_IN(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MASK_WAVE_OUT(:,:)

LOGICAL :: LWVIN_MASK_NOT_SET
LOGICAL :: LWVIN_UNINITIALISED

INTEGER(KIND=JPIM) :: MWVIN_SENDTOT
INTEGER(KIND=JPIM) :: MWVIN_RECVTOT

INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_SENDCNT(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_RECVCNT(:)

INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_SENDOFF(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_RECVOFF(:)

INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_SENDBUF(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_RECVBUF(:)

INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_SENDIND(:)


!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(cbegdat,lwcou,lwcou2w,lwcounorms,lwvin_mask_not_set,lwvin_uninitialised,mwvin_recvtot)
!$OMP THREADPRIVATE(mwvin_sendtot,ndurat,nlat1w,nlatw,nlon1w,nlonw,nnorxw,nresum,nstpw,rdegrew,rnortw,rsoutw)
END MODULE YOEWCOU
