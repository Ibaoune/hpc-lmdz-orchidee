MODULE TYPE_GFLS

USE PARKIND1  ,ONLY : JPIM ,JPRB   

IMPLICIT NONE
SAVE

!-------------------------------------------------------------------------
! Derived types for describing the GFL structure. The descriptors themselves
! (YGFL and YGFLC) can be found in module yom_ygfl.F90.
!-------------------------------------------------------------------------
! Modifications:
! 03/07/09 C. Fischer - add Arome/Aladin attributes
! 03/10/01 C. Moussy  - add Arome/Aladin attributes coupling
! 03/10/31 M. Tudor   - add physics tendencies for predictor-corrector
! 05/10/10 J. Haseler - switch for I/O to trajectory structure
! 2004-Nov F. Vana    - update of CSLINT attribute
! 20-Feb-2005 Vivoda  - 3TL Eul PC scheme (GFLPC)
 
TYPE TYPE_GFLD
! Overall descriptor,dimensioning etc.
INTEGER(KIND=JPIM) :: NUMFLDS          ! Number of GFL fields
INTEGER(KIND=JPIM) :: NDERS            ! Number of horizontal derivatives fields
INTEGER(KIND=JPIM) :: NUMSPFLDS        ! Number of spectrally represented GFL fields
INTEGER(KIND=JPIM) :: NUMGPFLDS        ! Number of grid-point GFL fields
INTEGER(KIND=JPIM) :: NUMFLDS9         ! Number of GFL fields in (t-dt) part
INTEGER(KIND=JPIM) :: NUMFLDS1         ! Number of GFL fields in (t+dt) array
INTEGER(KIND=JPIM) :: NUMSPFLDS1       ! Number of spectrally represented GFL fields (t+dt)
INTEGER(KIND=JPIM) :: NUMFLDS5         ! Number of GFL fields (trajectory)
INTEGER(KIND=JPIM) :: NUMFLDSPHY       ! Number of GFL fields (phys.)
INTEGER(KIND=JPIM) :: NUMFLDS_SPL      ! Number of GFL fields (S.L. spline interpolation)
INTEGER(KIND=JPIM) :: NUMFLDS_SL1      ! Number of GFL fields in S.L. buffer 1
INTEGER(KIND=JPIM) :: NUMFLDSPC        ! Number of GFL fields (predictor/corrector)
INTEGER(KIND=JPIM) :: NDIM             ! Dimension of main array holding GFL fields(GFL)
INTEGER(KIND=JPIM) :: NUMFLDSPT        ! Number of GFL fields (phy. tend.)
INTEGER(KIND=JPIM) :: NDIM0            ! Dimension of t0 part of GFL
INTEGER(KIND=JPIM) :: NDIM9            ! Dimension of t-dt part of GFL
INTEGER(KIND=JPIM) :: NDIM1            ! Dimension of t+dt array (GFLT1)
INTEGER(KIND=JPIM) :: NDIM5            ! Dimension of traj. GFL array (GFL5)
INTEGER(KIND=JPIM) :: NDIMSLP          ! Diminsion of S.L. phys. GFL array (GFLSLP)
INTEGER(KIND=JPIM) :: NDIM_SPL         ! Dim. of arrays holding GFL fields (S.L.spline int.)
INTEGER(KIND=JPIM) :: NDIMPT           ! Dimension of phy. tend. GFL array (GFLPT)
INTEGER(KIND=JPIM) :: NDIMPC           ! Dimension of iterative scheme auxiliary array (GFLPC)

END TYPE TYPE_GFLD

TYPE TYPE_GFL_COMP ! Individual field descriptor

CHARACTER(LEN=16) :: CNAME      ! ARPEGE field name 
INTEGER(KIND=JPIM) :: IGRBCODE           ! GRIB code
LOGICAL :: LADV                 ! Field advected or not
INTEGER(KIND=JPIM) :: NREQIN    ! 1 if field requiered in input, 0 if not, -1 if initialised
                                ! with a reference value REFVALI
LOGICAL :: LREQOUT              ! T if field requiered in output
LOGICAL :: LGPINGP              ! GP field input as GP
LOGICAL :: LGP                  ! Field exists and of grid-point type
LOGICAL :: LSP                  ! Field exists and of spectral type
LOGICAL :: LCDERS               ! Derivatives required (spectral only)
LOGICAL :: LACTIVE              ! Field in use
LOGICAL :: LTHERMACT            ! Field thermodynamically active
REAL(KIND=JPRB) :: R
REAL(KIND=JPRB) :: RCP
LOGICAL :: LT9                  ! Field in t-dt GFL
LOGICAL :: LT1                  ! Field in t+dt GFL
LOGICAL :: LT5                  ! Field in trajectory GFL
LOGICAL :: LPHY                 ! Field in physics GFL
LOGICAL :: LPT                  ! Field in PC phy. tend. GFL (GFLPT)
LOGICAL :: LTRAJIO              ! Field written to and from trajectory structure
LOGICAL :: LPC                  ! Field in predictor/corrector time stepping (GFLPC)
REAL(KIND=JPRB) :: REFVALI      ! Reference value for init, used in case NREQIN==-1
! LAM specific attributes (Arome/Aladin)
LOGICAL :: LADJUST0             ! True if field is thermodynamically adjusted at t
                                  ! (immediatly after inverse spectral transforms)
LOGICAL :: LADJUST1             ! True if field is thermodynamically adjusted at t+dt
                                  ! (after SL interpolations and NL residuals)
INTEGER(KIND=JPIM) :: NCOUPLING ! 1 if field is coupled by Davies relaxation, 0 if not,
                                  ! -1 if coupled with reference value for coupling REFVALC
REAL(KIND=JPRB) :: REFVALC      ! Reference value for coupling, used in case NCOUPLING==-1
LOGICAL :: LBIPER               ! True if field must be biperiodised inside the transforms
! End LAM specific attributes (Arome/Aladin)
CHARACTER(LEN=12) :: CSLINT     ! S.L interpolaion "type"
INTEGER(KIND=JPIM) :: MP        ! Basic field "pointer"
INTEGER(KIND=JPIM) :: MPL       ! zonal derivative "pointer"
INTEGER(KIND=JPIM) :: MPM       ! Meridional derivative "pointer"
INTEGER(KIND=JPIM) :: MP9       ! Basic field "pointer" t-dt
INTEGER(KIND=JPIM) :: MP9_PH    ! Basic field "pointer" for Physics
INTEGER(KIND=JPIM) :: MP1       ! Basic field "pointer" t+dt
INTEGER(KIND=JPIM) :: MP5       ! Basic field "pointer" trajectory
INTEGER(KIND=JPIM) :: MP5L      ! zonal derivative "pointer" trajectory
INTEGER(KIND=JPIM) :: MP5M      ! Meridional derivative "pointer" trajectory
INTEGER(KIND=JPIM) :: MPSLP     ! Basic field "pointer" physics
INTEGER(KIND=JPIM) :: MPSP      ! Basic field "pointer" spectral space
INTEGER(KIND=JPIM) :: MP_SPL    ! Basic field "pointer" spline interpolation
INTEGER(KIND=JPIM) :: MP_SL1    ! Basic field "pointer" in SLBUF1
INTEGER(KIND=JPIM) :: MP_SLX    ! Basic field "pointer" in SLBUF1 for CPG_PT
INTEGER(KIND=JPIM) :: MPPT      ! Physics tendency "pointer"
INTEGER(KIND=JPIM) :: MPPC      ! Predictor/corrector auxiliary array "pointer"

TYPE(TYPE_GFL_COMP),POINTER :: PREVIOUS ! Pointer to previously def. field

END TYPE TYPE_GFL_COMP

TYPE TYPE_GFL_NAML ! Individual field descriptor for namelist input

CHARACTER(LEN=16) :: CNAME      ! ARPEGE field name 
INTEGER(KIND=JPIM) :: IGRBCODE  ! GRIB code
INTEGER(KIND=JPIM) :: NREQIN    ! 1 if field required in input, 0 if not, -1 if initialised
                                ! with a reference value REFVALI
REAL(KIND=JPRB) :: REFVALI      ! Reference value for initialisation, used in case NREQIN==-1
LOGICAL :: LREQOUT              ! T if field requiered in output
LOGICAL :: LGPINGP              ! GP field input as GP
LOGICAL :: LGP                  ! Field exists and of grid-point type
LOGICAL :: LSP                  ! Field exists and of spectral type
LOGICAL :: LCDERS               ! Derivatives required (spectral only)
LOGICAL :: LT9                  ! Field in t-dt GFL
LOGICAL :: LT1                  ! Field in t+dt GFL
LOGICAL :: LT5                  ! Field in trajectory GFL
LOGICAL :: LPHY                 ! Field with physics tendencies GFL
LOGICAL :: LPT                  ! Field in PC physics tendency GFLPT
LOGICAL :: LTRAJIO              ! Field written to and from trajectory structure
LOGICAL :: LPC                  ! Field in predictor/corrector time stepping GFLPC
LOGICAL :: LADV                 ! Field advected or not

LOGICAL :: LQM                  ! quasi-monotonous interpolation for field
LOGICAL :: LQMH                 ! quasi-monotonous interpolation in horizontal for field
LOGICAL :: LSLHD                ! Semi-lagrangian horizontal diffusion used for fiels
LOGICAL :: LRSPLINE             ! 12 points spline interpolation used for field
LOGICAL :: LHV                  ! Hermite vertical interpolation used for field (only ozone sofar)
LOGICAL :: LVSPLIP              ! vertical spline interpolation used for field (only ozone sofar)
INTEGER(KIND=JPIM) :: NCOUPLING ! 1 if field is coupled by Davies relaxation, 0 if not,
                                ! -1 if coupled with reference value for coupling REFVALC
REAL(KIND=JPRB) :: REFVALC      ! Reference value for coupling, used in case 
                                  ! NCOUPLING==-1
END TYPE TYPE_GFL_NAML

!-------------------------------------------------------------------------
END MODULE TYPE_GFLS
