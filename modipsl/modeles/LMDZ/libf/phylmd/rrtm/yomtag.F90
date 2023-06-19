MODULE YOMTAG

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    Tag identifiers used in message passing communication

! MTAGLM: tag for transpositions done in TRLTOM.
! MTAGMS: tag for transpositions done in TRMTOS.
! MTAGSM: tag for transpositions done in TRSTOM.
! MTAGMV: tag for transpositions done in TRMTOV.
! MTAGVH: tag for transpositions done in TRVTOH.
! MTAGMN: tag for transpositions done in TRMTON.
! MTAGNM: tag for transpositions done in TRNTOM.
! MTAGSPNO: tag for communications done in COMMSPNORM and COMMSPNORM1.
! MTAGSLAG: tag for halo constitution (horizontal interpolations
!  in the semi-Lagrangian scheme, the observation interpolator or FULLPOS).
! MTAGTIDE: tag for communications done in COMMTIDE.
! MTAGRAD : tag for communications done in SUECRAD (ECMWF physics)
! MTAGRCBDY: tag for communications done in RADCBDY (ECMWF physics).
! MTAGRCLB: tag for communications done in RADCLB (ECMWF physics).
! MTAGRCLBI: tag for communications done in SUECRADL (ECMWF physics).
! MTAGPART: tag for communications done in DICOMOUT and GATHFLNM.
! MTAGDISTSP: tag for communications done in DISSPEC, DISSPEC0 and DIWRSPE.
! MTAGDISTGP: tag for communications done in
!  DISGRID, DISGRID_C, DISGRIDFP, DIWRGRFP, DIWRGRID, IRCVGPF, IRCVGPFFP,
!  ISNDGPF, ISNDGPFFP, ORCVGPF, ORCVGPFFP, OSNDGPF, OSNDGPFFP.
! MTAGCAIN: tag for communications done in GATHERSPA.
! MTAGCOST: tag for communications done in
!  GATHERCOST1, GATHERCOST2, GATHERCOSTO and GATHERJCVERT.
! MTAGGSUM: tag for communications done in CASND1, CASNDR1 and GATHERSUM.
! MTAGGLOBSI: tag for communications done in CAEXCO and CAUPDO.
! MTAGGLOBSR: tag for communications done in CAEXCO and CAUPDO.
! MTAGOBSEQ: tag for communications done in MPOBSEQ.
! MTAGOBSEQAD: tag for communications done in MPOBSEQAD.
! MTAGFCE: tag for communications done in
!  COMMFCE1, COMMFCE2, COMMJBBAL and COMMJBDAT.
! MTAGBDY: tag for communications done in GATHERBDY.
! MTAGSIG: tag for communications done in SIGCHECK.
! MTAGBRPR: tag for communications done in BRPTOB and GATHERT.
! MTAGGPNORM: tag for communications done in GPNORM1.
! MTAGDDHRES: tag for communications done in DDHRCV and DDHSND.
! MTAGDDH1: tag for communications done in DISTDDH.
! MTAGDDH2: tag for communications done in DLADDH.
! MTAGDDH3: tag for communications done in DMADDH.
! MTAGDDH4: tag for communications done in DRESDDH.
! MTAGGETV: tag for communications done in SUHESS.
! MTAGOZON: tag for communications done in UPDO3CH.
! MTAGREADVEC: tag for communications done in READVEC.
! MT_DISTRIBUTED_VECTOR: tag for communications done in SUMPINI.
! MTAGLCZ: tag for communications done in COMMNSEC1.
! MTAGGOM: tag for communications done in GATHERGOM.
! MTAGFREQ: tag for communications done in GATHERFREQ.
! MTAGEIGMD: tag for communications done in GATHEREIGMD.
! MTAGKE: tag for communications done in VMODEENERGY.
! MTAGDISTFO: tag for communications done in DISFOU and DIWRFOU.

!      YOMTAG

INTEGER(KIND=JPIM) :: MTAGLM
INTEGER(KIND=JPIM) :: MTAGMS
INTEGER(KIND=JPIM) :: MTAGSM
INTEGER(KIND=JPIM) :: MTAGMV
INTEGER(KIND=JPIM) :: MTAGVH
INTEGER(KIND=JPIM) :: MTAGSPNO
INTEGER(KIND=JPIM) :: MTAGSLAG
INTEGER(KIND=JPIM) :: MTAGTIDE
INTEGER(KIND=JPIM) :: MTAGRAD
INTEGER(KIND=JPIM) :: MTAGRCBDY
INTEGER(KIND=JPIM) :: MTAGRCLB
INTEGER(KIND=JPIM) :: MTAGRCLBI
INTEGER(KIND=JPIM) :: MTAGPART
INTEGER(KIND=JPIM) :: MTAGDISTSP
INTEGER(KIND=JPIM) :: MTAGDISTGP
INTEGER(KIND=JPIM) :: MTAGMN
INTEGER(KIND=JPIM) :: MTAGNM
INTEGER(KIND=JPIM) :: MTAGCAIN
INTEGER(KIND=JPIM) :: MTAGCOST
INTEGER(KIND=JPIM) :: MTAGGSUM
INTEGER(KIND=JPIM) :: MTAGGLOBSI
INTEGER(KIND=JPIM) :: MTAGGLOBSR
INTEGER(KIND=JPIM) :: MTAGOBSEQ
INTEGER(KIND=JPIM) :: MTAGOBSEQAD
INTEGER(KIND=JPIM) :: MTAGFCE
INTEGER(KIND=JPIM) :: MTAGBDY
INTEGER(KIND=JPIM) :: MTAGDDHRES
INTEGER(KIND=JPIM) :: MTAGSIG
INTEGER(KIND=JPIM) :: MTAGBRPR
INTEGER(KIND=JPIM) :: MTAGGPNORM
INTEGER(KIND=JPIM) :: MTAGDDH1
INTEGER(KIND=JPIM) :: MTAGDDH2
INTEGER(KIND=JPIM) :: MTAGDDH3
INTEGER(KIND=JPIM) :: MTAGDDH4
INTEGER(KIND=JPIM) :: MTAGGETV
INTEGER(KIND=JPIM) :: MTAGOZON
INTEGER(KIND=JPIM) :: MTAGREADVEC
INTEGER(KIND=JPIM) :: MT_DISTRIBUTED_VECTOR
INTEGER(KIND=JPIM) :: MTAGLCZ
INTEGER(KIND=JPIM) :: MTAGGOM
INTEGER(KIND=JPIM) :: MTAGFREQ
INTEGER(KIND=JPIM) :: MTAGEIGMD
INTEGER(KIND=JPIM) :: MTAGKE
INTEGER(KIND=JPIM) :: MTAGDISTFO

!$OMP THREADPRIVATE(mt_distributed_vector,mtagbdy,mtagbrpr,mtagcain,mtagcost,mtagddh1,mtagddh2,mtagddh3,mtagddh4)
!$OMP THREADPRIVATE(mtagddhres,mtagdistfo,mtagdistgp,mtagdistsp,mtageigmd,mtagfce,mtagfreq,mtaggetv,mtagglobsi)
!$OMP THREADPRIVATE(mtagglobsr,mtaggom,mtaggpnorm,mtaggsum,mtagke,mtaglcz,mtaglm,mtagmn,mtagms,mtagmv,mtagnm)
!$OMP THREADPRIVATE(mtagobseq,mtagobseqad,mtagozon,mtagpart,mtagrad,mtagrcbdy,mtagrclb,mtagrclbi,mtagreadvec)
!$OMP THREADPRIVATE(mtagsig,mtagslag,mtagsm,mtagspno,mtagtide,mtagvh)
END MODULE YOMTAG
