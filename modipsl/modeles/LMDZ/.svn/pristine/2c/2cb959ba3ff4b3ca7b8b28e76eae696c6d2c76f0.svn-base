MODULE YOM_PHYS_GRID

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE

PRIVATE
PUBLIC TYPE_PHYS_POINT ,PHYS_GRID,DYN_GRID,DYN_SL,PHYS_SL, YPHYPOI, YDYNPOI, &
        & JPMXNEI, PHYS_GRID_STRUCT, SL_STRUCT, TYPE_DYN_POINT

!-------------------------------------------------------------------------
! Derived types for describing the coarse physics grid structure. 
!  The descriptors themselves
! (YGFL and YGFLC) can be found in module yom_ygfl.F90.
!-------------------------------------------------------------------------
! Modifications:

INTEGER(KIND=JPIM), PARAMETER :: JPMXNEI=36 ! maximum number of neighbouring
                                            ! points for averaging or
                                            ! interpolation
                                              
TYPE PHYS_GRID_STRUCT
  INTEGER(KIND=JPIM) :: NGPTOT          ! number of physics points in the task
  INTEGER(KIND=JPIM) :: NGPTOTG
  INTEGER(KIND=JPIM) :: NDGSAL, NDGENL
  INTEGER(KIND=JPIM) :: NDGSAG, NDGENG
  INTEGER(KIND=JPIM), POINTER :: NLOENG(:)
  REAL(KIND=JPRB), POINTER   :: RMU(:)
  REAL(KIND=JPRB), POINTER   :: RW(:)

  INTEGER(KIND=JPIM)   :: NRESOL_ID
  INTEGER(KIND=JPIM)   :: NGPTOTMX, NSPEC2, NSMAX
  INTEGER(KIND=JPIM)   :: NPTRFLOFF, NUMP, NDLON
  INTEGER(KIND=JPIM)   ::  NDGSAH, NDGENH
  INTEGER(KIND=JPIM)   ::  NDGLG, NDLSUR
  INTEGER(KIND=JPIM)   ::  NFRSTLOFF, NDSUR1,NDGSUR
  INTEGER(KIND=JPIM)   ::  MYFRSTACTLAT, MYLSTACTLAT
  INTEGER(KIND=JPIM)   :: NGPBLKS
  INTEGER(KIND=JPIM)   :: NPROMA

  INTEGER(KIND=JPIM), POINTER, DIMENSION(:)     :: NRGRI, NPTRFRSTLAT, NFRSTLAT
  INTEGER(KIND=JPIM), POINTER, DIMENSION(:)     :: NLSTLAT, MYMS, NASM0
  INTEGER(KIND=JPIM), POINTER, DIMENSION(:,:)   :: NSTA, NONL
  INTEGER(KIND=JPIM), POINTER, DIMENSION(:)     :: NSTAGP
  REAL(KIND=JPRB), POINTER   :: RSQM2(:), RLATIG(:), RLATI(:)
  REAL(KIND=JPRB), POINTER   :: RIPI0(:), RIPI1(:), RIPI2(:)
  INTEGER(KIND=JPIM), POINTER, DIMENSION(:)     :: NPTRLSTLAT

END TYPE PHYS_GRID_STRUCT

TYPE(PHYS_GRID_STRUCT) :: PHYS_GRID, DYN_GRID 

TYPE SL_STRUCT
  INTEGER(KIND=JPIM),POINTER :: NSLSTA(:)
  INTEGER(KIND=JPIM),POINTER :: NSLONL(:)
  INTEGER(KIND=JPIM),POINTER :: NSLOFF(:)
  INTEGER(KIND=JPIM),POINTER :: NSLEXT(:,:)
  INTEGER(KIND=JPIM),POINTER :: NSLSENDPOS(:)
  INTEGER(KIND=JPIM),POINTER :: NSLRECVPOS(:)
  INTEGER(KIND=JPIM),POINTER :: NSLSENDPTR(:)
  INTEGER(KIND=JPIM),POINTER :: NSLRECVPTR(:)
  INTEGER(KIND=JPIM),POINTER :: NSLCORE(:)
  INTEGER(KIND=JPIM),POINTER :: NSLCOMM(:)

  INTEGER(KIND=JPIM) :: NASLB1
  INTEGER(KIND=JPIM) :: NSLPROCS
  INTEGER(KIND=JPIM) :: NSLMPBUFSZ
  INTEGER(KIND=JPIM) :: NSLRPT
  INTEGER(KIND=JPIM) :: NSLSPT
  INTEGER(KIND=JPIM) :: NSLWIDEN
  INTEGER(KIND=JPIM) :: NSLWIDES
  INTEGER(KIND=JPIM) :: NSLWIDEE
  INTEGER(KIND=JPIM) :: NSLWIDEW
END TYPE SL_STRUCT

TYPE(SL_STRUCT) :: DYN_SL
TYPE(SL_STRUCT) :: PHYS_SL

 
TYPE TYPE_PHYS_POINT    ! Individual physics point characteristics

  REAL(KIND=JPRB) :: GELAM, GELAT, GEMU
  REAL(KIND=JPRB) :: GECLO, GESLO, GM, GAW
  REAL(KIND=JPRB) :: GNORDL, GNORDM, GSQM2
  REAL(KIND=JPRB) :: RCOLON, RSILON 
  REAL(KIND=JPRB) :: RINDX, RINDY 
  REAL(KIND=JPRB) :: OROG
  INTEGER(KIND=JPIM) :: NGPLAT       ! row number in the physics grid

  INTEGER(KIND=JPIM) :: NEIGH       ! number of neighbours in the dynamics grid 
                                    ! for going from the dynamics to the physics grid
  INTEGER(KIND=JPIM), POINTER :: NL0(:)    ! indexes in the interpolation buffer of 
                                           ! the dynamics neighbours
  REAL(KIND=JPRB), POINTER   :: WGT(:)     ! weights for every neighbouring 
                                           ! dynamics point 
END TYPE TYPE_PHYS_POINT

TYPE(TYPE_PHYS_POINT),ALLOCATABLE :: YPHYPOI(:)

TYPE TYPE_DYN_POINT     ! Individual dynamics point characteristics
  INTEGER(KIND=JPIM) :: NEIGH       ! number of neighbours in the physics grid 
                                    ! for going from the physics to the dynamics grid
  INTEGER(KIND=JPIM), POINTER :: NL0(:)    ! indexes in the interpolation buffer of 
                                           ! the physics neighbours
  REAL(KIND=JPRB), POINTER   :: WGT(:)     ! weights for every neighbouring 
                                           ! physics point 
END TYPE TYPE_DYN_POINT

TYPE(TYPE_DYN_POINT),ALLOCATABLE :: YDYNPOI(:)
 
!$OMP THREADPRIVATE(dyn_grid,dyn_sl,phys_grid,phys_sl)
!$OMP THREADPRIVATE(ydynpoi,yphypoi)
END MODULE YOM_PHYS_GRID
