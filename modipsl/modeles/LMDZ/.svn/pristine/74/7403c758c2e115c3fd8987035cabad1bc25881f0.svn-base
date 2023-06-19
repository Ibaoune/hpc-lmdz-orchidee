MODULE YOMGEM

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    * Number of grid points

!     NGPTOT   : Total number of grid columns on a PE
!     NGPTOT_CAP  : Size of grid points arrays for ALADIN
!     NGPTOTMX : Maximum number of grid columns on any of the PEs
!     NGPTOTG  : Total number of grid columns on the Globe
!     NGPTOTL(NPRGPNS,NPRGPEW)  : Total number of grid columns on on eacch PE

INTEGER(KIND=JPIM) :: NGPTOT
INTEGER(KIND=JPIM) :: NGPTOT_CAP
INTEGER(KIND=JPIM) :: NGPTOTMX
INTEGER(KIND=JPIM) :: NGPTOTG
INTEGER(KIND=JPIM),ALLOCATABLE, TARGET :: NGPTOTL(:,:)

!     ------------------------------------------------------------------

!*    * Defining the transformed sphere

!     RMUCEN : MU OF THE POLE OF STRETCHING
!     RLOCEN : LONGITUDE OF THE POLE OF STRETCHING
!     RSTRET : STRETCHING FACTOR
!     NSTTYP : 1 = POLE OF STRETCHING, POLE OF THE COLLOCATION GRID
!                 AT THE NORTHERN POLE OF THE REAL EARTH.
!              2 = THE POLE OF STRETCHING IS ANYWHERE ON THE REAL EARTH
!             AND ON THE EQUATOR OF THE COLLOCATION GRID ON THE MERIDIAN PI.
!                  THE EQUATOR OF THE COLLOCATION GRID IS TANGENT
!             TO A PARALLEL OF THE EARTH.

!     NHTYP  : 0 = regular grid
!            : 1 = number of points proportional to sqrt(1-mu**2)
!            : 2 = number of points read on namelist namrgri

!     RNLGINC: increment to get non-linear grid

!     R4JP    inverse de delta(teta) approche a l'ordre 1
!     RC2P1   RSTRET*RSTRET+1.
!     RC2M1   RSTRET*RSTRET-1.
!     RCOR0   COMPONENT (0,0) OF CORIOLIS
!     RCOR1   COMPONENT (0,1) OF CORIOLIS
!     RCOR2   COMPONENT (1,1) OF CORIOLIS

!     RCOLON(NGPTOT) cosine of longitude on transformed sphere
!     RSILON(NGPTOT)   sine        "             "         "
!     RINDX (NGPTOT) Longitude index
!     RINDY (NGPTOT) Latitude index
!     RATATH(NGPTOT) RA*TAN(THETA) on real sphere
!     RATATX(NGPTOT) Curvature term for LAM (for u eq.)

!     NLOEN(NDGSAG:NDGENG)  : number of active points on a parallel
!     NLOENG(NDGSAG:NDGENG) : global version of NLOEN
!     NMEN(NDGSAG:NDGENG)   : associated cut-off wave number
!     NMENTC(NDGSAG:NDGENG) : same as NMEN but for truncation NTCMAX.
!     NMENG(NDGSAG:NDGENG)   : global version of NMEN
!     NDGLU(0:MAX(NSMAX,NMSMAX)) : number of active points in an hemisphere
!                          for a given wave number m
!     NSTAGP(NGPTOT)     : start position of latitude data for boundary fields
!     NESTAGP(NGPTOT)    : start position of latitude data for boundary fields
!                          in extension zone (ALADIN).
!     NTSTAGP(NGPTOT)    : start position of latitude data for boundary fields
!                          in C+I+E zone (ALADIN).


!  REFLRHC :  reference length for critical relative humidity (cloud scheme)
!  REFLKUO, REFLCAPE, REFLRHC are ONLY to be used in the SETUP 
!  of TEQK and TEQC and TEQH

!  TEQK : ratio between REFLKUO and the model equivalent mesh size
!  TEQC : ratio between REFLCAPE and the model equivalent mesh size
!  TEQH : ratio between REFLRHC and the model equivalent mesh size

!  TEQK, TEQC, TEQH are to be used in the PHYSICS.

REAL(KIND=JPRB) :: RMUCEN
REAL(KIND=JPRB) :: RLOCEN
REAL(KIND=JPRB) :: RSTRET
INTEGER(KIND=JPIM) :: NSTTYP
INTEGER(KIND=JPIM) :: NHTYP
REAL(KIND=JPRB) :: RNLGINC
REAL(KIND=JPRB) :: R4JP
REAL(KIND=JPRB) :: RC2P1
REAL(KIND=JPRB) :: RC2M1
REAL(KIND=JPRB) :: RCOR0
REAL(KIND=JPRB) :: RCOR1
REAL(KIND=JPRB) :: RCOR2
REAL(KIND=JPRB),ALLOCATABLE:: RCOLON(:)
REAL(KIND=JPRB),ALLOCATABLE:: RSILON(:)
REAL(KIND=JPRB),ALLOCATABLE:: RINDX(:)
REAL(KIND=JPRB),ALLOCATABLE:: RINDY(:)
REAL(KIND=JPRB),ALLOCATABLE:: RATATH(:)
REAL(KIND=JPRB),ALLOCATABLE:: RATATX(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NLOEN(:)
INTEGER(KIND=JPIM),ALLOCATABLE,TARGET :: NLOENG(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NMEN(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NMENTC(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NMENG(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NDGLU(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NSTAGP(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NESTAGP(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NTSTAGP(:)

!     ------------------------------------------------------------------

!*    * Defining the transformed sphere: physics input

REAL(KIND=JPRB) :: REFLRHC
REAL(KIND=JPRB) :: TEQH
REAL(KIND=JPRB) :: REFLKUO
REAL(KIND=JPRB) :: REFLCAPE
REAL(KIND=JPRB) :: TEQK
REAL(KIND=JPRB) :: TEQC

!     ------------------------------------------------------------------

!*    * DEFINING THE VERTICAL COORDINATE

!     VP00  : REFERENCE PRESSURE FOR DEFINING VERTICAL COORDINATE
!     VALH  : (0:NFLEVG)
!     VBH   : (0:NFLEVG) : B of the vertical coordinate
!     VETAH : (0:NFLEVG) ; VERTICAL COORDINATE = VALH+VBH
!     VETAF : (0:NFLEVG+1) ; VERTICAL COORDINATE ON LAYERS.
!     VCUICO: is used to compute denominators of weights
!             for semi-Lagrangian vertical interpolations
!             applied to full-level variables.
!     VCUICOH:is used to compute denominators of weights 
!             for semi-Lagrangian vertical interpolations
!             applied to half-level variables.
!     VRLEVX: REAL(NRLEVX)
!     NVAUTF: NVAUTF(VRLEVX*eta) is the number of the layer (full level)
!             immediately above "eta", and is bounded by 1 and nflevg-1.
!     NVAUTH: NVAUTH(VRLEVX*eta) is the number of the interlayer (half level)
!             immediately above "eta", and is bounded by 0 and nflevg-1.
!     VAH   : (0:NFLEVG) ;  =VALH*VP00
!     VC    : (NFLEVG)   ;  =VAH(J)*VBH(J-1)-VAH(J-1)*VBH(J)
!     VDELB : (NFLEVG)   ;  =VBH(J)-VBH(J-1)
!     VDELA : (NFLEVG)   ;  =VAH(J)-VAH(J-1)
!     VAF   : like VAH but at full levels.
!     VBF   : like VBH but at full levels.
!     VRDETAH: 1/[Delta eta]
!     TOPPRES: REFERENCE "EVANESCENT" PRESSURE
!              TOPPRES allows to solve some calculations of singularities
!              when the top pressure of the model is zero (for ex. in
!              GPPREF, GPXYB, SUNHBMAT).

!     WE HAVE THEN FOR THE HALF LEVEL PRESSURE : VAH + VBH*(SURFACE PRESSURE)

!     NOTE THAT THE HALF LEVEL VALUE AT K+.5 IS VXXX(K)
!     (THE FULL LEVEL VALUES ARE FROM 1 TO NFLEVG)

REAL(KIND=JPRB) :: VP00
REAL(KIND=JPRB),ALLOCATABLE:: VALH(:)
REAL(KIND=JPRB),ALLOCATABLE:: VBH(:)
REAL(KIND=JPRB),ALLOCATABLE:: VETAH(:)
REAL(KIND=JPRB),ALLOCATABLE:: VETAF(:)
REAL(KIND=JPRB),ALLOCATABLE:: VCUICO(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: VCUICOH(:,:)
REAL(KIND=JPRB) :: VRLEVX
INTEGER(KIND=JPIM),ALLOCATABLE:: NVAUTF(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NVAUTH(:)
REAL(KIND=JPRB),ALLOCATABLE:: VAH(:)
REAL(KIND=JPRB),ALLOCATABLE:: VC(:)
REAL(KIND=JPRB),ALLOCATABLE:: VDELB(:)
REAL(KIND=JPRB),ALLOCATABLE:: VDELA(:)
REAL(KIND=JPRB),ALLOCATABLE:: VAF(:)
REAL(KIND=JPRB),ALLOCATABLE:: VBF(:)
REAL(KIND=JPRB),ALLOCATABLE:: VRDETAH(:)
REAL(KIND=JPRB) :: TOPPRES

!     ------------------------------------------------------------------

!*    * Miscellaneous

!     NBEEGP : ???
!     NBNEGP : ???

INTEGER(KIND=JPIM) :: NBEEGP
INTEGER(KIND=JPIM) :: NBNEGP

!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(nbeegp,nbnegp,ngptot,ngptot_cap,ngptotg,ngptotmx,nhtyp,nsttyp,r4jp,rc2m1,rc2p1,rcor0)
!$OMP THREADPRIVATE(rcor1,rcor2,reflcape,reflkuo,reflrhc,rlocen,rmucen,rnlginc,rstret,teqc,teqh,teqk,toppres,vp00,vrlevx)
!$OMP THREADPRIVATE(ndglu,nestagp,ngptotl,nloen,nloeng,nmen,nmeng,nmentc,nstagp,ntstagp,nvautf,nvauth,ratath,ratatx)
!$OMP THREADPRIVATE(rcolon,rindx,rindy,rsilon,vaf,vah,valh,vbf,vbh,vc,vcuico,vcuicoh,vdela,vdelb,vetaf,vetah,vrdetah)
END MODULE YOMGEM
