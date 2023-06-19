MODULE YOMDIM

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Dimensions of model working arrays

! === COLLOCATION GRID OF THE DYNAMICS ========================================

! NDGLG  : number of rows of latitudes
! NDGLL  : number of rows of latitudes for which this process is
!          performing Fourier Space calculations
! NDGNH  : number of rows in the northern hemisphere
! NDGSUR : number of additional rows at each pole for horizontal
!          interpolations.
! NDGSAG = -NDGSUR+1
! NDGSAL = Local version of NDGSAG.
! NDGSAH = 1-NSLWIDE in DM version.
! NDGSAFPH=1-NFPWIDE in DM version.
! NDGENG = NDGLG+NDGSUR
! NDGENL = Number of latitude rows for which this process has grid
!          point calculations to perform.
! NDGENH = NDGENL+NSLWIDE in DM version.
! NDGENFPH=NDGENL+NFPWIDE in DM version.
! NDGUNG : first row of the area of interest in Aladin
!        = NDGSAG in the global model
! NDGUXG : last  row of the area of interest in Aladin
!        = NDGENG in the global model
! NDGUNL : local first row in C+I zone in distributed memory Aladin
! NDGUXL : local last row in C+I zone in distributed memory Aladin
! NDLON  : length of a row of latitude near equator
! NDSUR1 : over dimensioning of NDLON for technical reasons (at least 2)
! NDLSUR = NDLON+NDSUR1
! NDLSM  = NDLSUR-1
! NDLUNG : first meridian of the area of interest in Aladin
!        = 1 in the global model
! NDLUXG : last  meridian of the area of interest in Aladin
!        = NDLON in the global model
! NDLUNL : local first meridian in C+I zone in distributed memory Aladin
! NDLUXL : local last meridian in C+I zone in distributed memory Aladin
! NPROMA : working dimension for grid-point computations
! NPROMB : working dimension for a 2nd call to DMN physics
! NPROMC : working dimension for a 2nd call to DMN physics
! NPROME : working dimension for ECMWF physics computations
! NPROMM : working dimension for DMN   physics computations
! NPROMNH: working dimension for non hydrostatic
! NPROMNH_GWADV: working dimension for non hydrostatic (arrays used only
!          when LGWADV=T)
! NPROMP : working dimension for physics computations
! NPROMV : working dimension for variational gridpoint fields
! NPROMVC: working dimension for some additional grid-point arrays
!          used only when "lvercor=.T.".
! NGPBLKS: number of grid point NPROMA-blocks.
! LOPTPROMA : .TRUE. NPROMA will be optimised
!           : .FALSE. NPROMA will not be optimised (forced by
!           : negative NPROMA in namelist)

INTEGER(KIND=JPIM) :: NDGLG
INTEGER(KIND=JPIM) :: NDGLL
INTEGER(KIND=JPIM) :: NDGNH
INTEGER(KIND=JPIM) :: NDGSUR
INTEGER(KIND=JPIM) :: NDGSAG
INTEGER(KIND=JPIM) :: NDGSAL
INTEGER(KIND=JPIM) :: NDGSAH
INTEGER(KIND=JPIM) :: NDGSAFPH
INTEGER(KIND=JPIM) :: NDGENG
INTEGER(KIND=JPIM) :: NDGENL
INTEGER(KIND=JPIM) :: NDGENH
INTEGER(KIND=JPIM) :: NDGENFPH
INTEGER(KIND=JPIM) :: NDGUNG
INTEGER(KIND=JPIM) :: NDGUXG
INTEGER(KIND=JPIM) :: NDGUNL
INTEGER(KIND=JPIM) :: NDGUXL
INTEGER(KIND=JPIM) :: NDLON
INTEGER(KIND=JPIM) :: NDSUR1
INTEGER(KIND=JPIM) :: NDLSUR
INTEGER(KIND=JPIM) :: NDLSM
INTEGER(KIND=JPIM) :: NDLUNG
INTEGER(KIND=JPIM) :: NDLUXG
INTEGER(KIND=JPIM),ALLOCATABLE:: NDLUNL(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NDLUXL(:,:)
INTEGER(KIND=JPIM) :: NPROMA
INTEGER(KIND=JPIM) :: NPROMB
INTEGER(KIND=JPIM) :: NPROMC
INTEGER(KIND=JPIM) :: NPROME
INTEGER(KIND=JPIM) :: NPROMM
INTEGER(KIND=JPIM) :: NPROMNH
INTEGER(KIND=JPIM) :: NPROMNH_GWADV 
INTEGER(KIND=JPIM) :: NPROMP
INTEGER(KIND=JPIM) :: NPROMV
INTEGER(KIND=JPIM) :: NPROMVC
INTEGER(KIND=JPIM) :: NGPBLKS
LOGICAL :: LOPTPROMA

! === VERTICAL RESOLUTION =====================================================

! NFLEVG : number of levels in grid point space
! NFLEVL : number of levels in Fourier and Legendre space
! NFLEVLMX : maximum NFLEVL among all PEs
! NFLSUR : over dimensioning of NFLEVL for technical reasons, always odd
! NFLSUL : number of additional levels for semi-lagrangian
! NFLSA  = 1    -NFLSUL
! NFLEN  = NFLEVG+NFLSUL

INTEGER(KIND=JPIM) :: NFLEVL
INTEGER(KIND=JPIM) :: NFLEVLMX
INTEGER(KIND=JPIM) :: NFLEVG
INTEGER(KIND=JPIM) :: NFLSUR
INTEGER(KIND=JPIM) :: NFLSUL
INTEGER(KIND=JPIM) :: NFLSA
INTEGER(KIND=JPIM) :: NFLEN

! === NUMBER OF FIELDS ========================================================

! NFTHER : number of spectral thermodynamic variables
! NFAUX  : number of auxillary variables in t+dt array
! NF3D   = number of 3D fields in the state of the model
! NFD2D  : number of 2D fields in the dynamics
! NFC2D  : number of 2D fields in the boundaries
! NPPM   : Number of interpolation methods in post-processing
! NFPPYX : Maximum number of modern dyn.met. post-processed fields
! NFPPYE : as long as SUPP is called after SUDIM/SUALLO,
!            NFPPYE = 1, if LMDYPP = .FALSE. (after SUALLO)
!            NFPPYE = NFPPYX if LMDYPP = .TRUE., but it could become more flex
! NFGPNH : number of (3D) fields in non hydrostatic gridpoint
! NS3D   : number of 3D fields in spectral space
! NS2D   : number of 2D fields in spectral space
! NS1D   : number of 1D fields in spectral space (for Aladin consistency -CF)
! NSAUX  : dimension of auxillary array == NFAUX

INTEGER(KIND=JPIM) :: NFTHER
INTEGER(KIND=JPIM) :: NFAUX
INTEGER(KIND=JPIM) :: NF3D
INTEGER(KIND=JPIM) :: NFD2D
INTEGER(KIND=JPIM) :: NFC2D
INTEGER(KIND=JPIM) :: NPPM
INTEGER(KIND=JPIM) :: NFPPYX
INTEGER(KIND=JPIM) :: NFPPYE
INTEGER(KIND=JPIM) :: NFGPNH
INTEGER(KIND=JPIM) :: NS3D
INTEGER(KIND=JPIM) :: NS2D
INTEGER(KIND=JPIM) :: NS1D
INTEGER(KIND=JPIM) :: NSAUX

! === GRID POINT ARRAYS =======================================================

! LVOR  : controls the allocation of vorticity
! LADER : controls the allocation of vor div and derivatives
! LUVDER: controls the allocation of derivatives for u and v
! LSPT  : .TRUE. if temperature variable as spectral field

LOGICAL :: LVOR
LOGICAL :: LADER
LOGICAL :: LUVDER
LOGICAL :: LSPT

! === SPECTRAL SPACE ==========================================================

! NSMAX  : truncation order
! NMSMAX  : truncation order in longitude
! NVARMAX: truncation order in 3d-var distributed direction
!          this is a priori longitude, so that nvarmax = nsmax in Arp/IFS
!          and nvarmax = nmsmax in Aladin
! NSEFRE : number of degrees of freedom in the spectral space
! NSPECG : number of complex spectral coefficients (global)
! NSPEC2G = 2*NSPECG
! NSPEC  : number of complex spectral coefficients (local, i.e. on this PE)
! NSPEC2 = 2*NSPEC
! NSPEC2MX : maximun NSPEC2 among all PEs
! NSMIN  : lower troncature for configurations 911 and 912.
! NTCMAX : truncation order for transmission coefficients.
! NMTCMAX: truncation order for transmission coefficients in longitude.
! NCMAX  : upper trunc. order for dilatation matrices (used in TRAGEO, fullpos)
! NCPEC  : number of complex spectral coefficients for truncation NCMAX (local)
! NCPEC2 : 2*NCPEC, where NCPEC (local)
! NXMAX  : truncation order for NMI
! NXPECG : number of complex spectral coefficients for NMI (global)
! NXPEC  : number of complex spectral coefficients for NMI (local)
! NTMAX  : truncation order for tendencies (on n, m<= NSMAX)
! NTPEC2 : 2*'number of complex spectral coefficients' for tendencies

INTEGER(KIND=JPIM) :: NSMAX
INTEGER(KIND=JPIM) :: NMSMAX
INTEGER(KIND=JPIM) :: NVARMAX
INTEGER(KIND=JPIM) :: NSEFRE
INTEGER(KIND=JPIM) :: NSPECG
INTEGER(KIND=JPIM) :: NSPEC2G
INTEGER(KIND=JPIM) :: NSPEC
INTEGER(KIND=JPIM) :: NSPEC2
INTEGER(KIND=JPIM) :: NSPEC2MX
INTEGER(KIND=JPIM) :: NSMIN
INTEGER(KIND=JPIM) :: NTCMAX
INTEGER(KIND=JPIM) :: NMTCMAX
INTEGER(KIND=JPIM) :: NCMAX
INTEGER(KIND=JPIM) :: NCPEC
INTEGER(KIND=JPIM) :: NCPEC2
INTEGER(KIND=JPIM) :: NXMAX
INTEGER(KIND=JPIM) :: NXPECG
INTEGER(KIND=JPIM) :: NXPEC
INTEGER(KIND=JPIM) :: NTMAX
INTEGER(KIND=JPIM) :: NTPEC2

! === DISTRIBUTED MEMORY DIMENSIONS ===========================================

! NUMP  :  Number of spectral waves handled by this processor
! NUMXP :  Same as NUMP, but related to NXMAX
! NUMCP :  Same as NUMP, but related to NCMAX
! NUMTP :  Same as NUMP, but related to NTMAX

INTEGER(KIND=JPIM) :: NUMP
INTEGER(KIND=JPIM) :: NUMXP
INTEGER(KIND=JPIM) :: NUMCP
INTEGER(KIND=JPIM) :: NUMTP

! === OTHER QUANTITIES ========================================================

! NRLEVX: dimension of NVAUTF in YOMGEM.
! NUNDEFLD: index value for unused/undefined fields (default=-9999)
!         : should be set to 1 when using compiler subscript checking

INTEGER(KIND=JPIM) :: NRLEVX
INTEGER(KIND=JPIM) :: NUNDEFLD

!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(lader,loptproma,lspt,luvder,lvor,ncmax,ncpec,ncpec2,ndgenfph)
!$OMP THREADPRIVATE(ndgeng,ndgenh,ndgenl,ndglg,ndgll,ndgnh,ndgsafph,ndgsag,ndgsah)
!$OMP THREADPRIVATE(ndgsal,ndgsur,ndgung,ndgunl,ndguxg,ndguxl,ndlon,ndlsm,ndlsur)
!$OMP THREADPRIVATE(ndlung,ndluxg,ndsur1,nf3d,nfaux,nfc2d,nfd2d,nfgpnh,nflen,nflevg)
!$OMP THREADPRIVATE(nflevl,nflevlmx,nflsa,nflsul,nflsur,nfppye,nfppyx,nfther,ngpblks)
!$OMP THREADPRIVATE(nmsmax,nmtcmax,nppm,nproma,npromb,npromc,nprome,npromm,npromnh)
!$OMP THREADPRIVATE(npromnh_gwadv,npromp,npromv,npromvc,nrlevx,ns1d,ns2d,ns3d,nsaux)
!$OMP THREADPRIVATE(nsefre,nsmax,nsmin,nspec,nspec2,nspec2g,nspec2mx,nspecg,ntcmax,ntmax)
!$OMP THREADPRIVATE(ntpec2,numcp,nump,numtp,numxp,nundefld,nvarmax,nxmax,nxpec,nxpecg)
!$OMP THREADPRIVATE(ndlunl,ndluxl)
END MODULE YOMDIM
