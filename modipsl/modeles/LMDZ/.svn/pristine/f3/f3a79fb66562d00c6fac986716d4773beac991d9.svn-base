MODULE YOMMP

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

! ----------------------------------------------------------------------
!*    variables describing distributed memory parallelization

! ---------------------------------------

!  mp_type     :  1=blocked   (MPI_SEND/RECV)
!              :  2=buffered  (MPI_BSEND/MPI_BRECV)
!              :  3=immediate (MPI_ISEND/MPI_IRECV)
!  mbx_size    :  user-provided mailbox size

!  myproc      :  logical processor id (is in the range 1 to nproc)
!  myseta      :  own processor set a (is in the range 1 to nprgpns)
!  mysetb      :  own processor set b (is in the range 1 to nprgpew)
!  my_region_ns:  own processor set a (is in the range 1 to n_regions_ns)
!  my_region_ew:  own processor set b (is in the range 1 to n_regions_ew)
!  mysetw      :  own processor set a in wave space (1..nprtrw)   
!  mysetv      :  own processor set b in wave space (1..nprtrv)    
!  mysetm      :  own processor set a in spectral space (1..nprtrm)    
!  mysetn      :  own processor set b in spectral space (1..nprtrn)    
!  mysetaf     :  own processor set a in Fourier space (is in the range
!                   1 to nprocc)
!  ngpset2pe   :  grid point space processor mapping array (n_regions_ns,n_regions_ew)
!  nslpad      :  number of pad words initialised to a huge number at either
!                 of side of the sl halo, used to trap halo problems.
!                 The default is 0. 
!  nintype     :  type in input processing to be performed
!              :  1=pbio
!              :  2=mpi-io (future)
!  nouttype    :  type of output (post) processing to be performed
!              :  1=pbio
!              :  2=output to FDB
!              :  3=shared blocking MPI-I/O
!              :  4=shared blocking collective MPI-I/O
!              :  5=shared non-blocking MPI_I/O
!              :  6=shared non-blocking collective MPI_I/O
!  nstrin      :  number of processors required to perform input processing
!  nstrout     :  number of processors required to perform output processing
!  ngathout    :  to be described
!  nwrtout     :  to be described
!  nblkout     :  to be described
!  nfldin      :  number of input  fields to be buffered during distribution
!  nfldout     :  number of output fields to be buffered during gathering
!  nprcids(nproc) : array containing the process ids. It is the mapping
!                 between the process numbering in the application
!                 (from 1 to NPROC) and the numbering used by the
!                 underlying communication library.

!  lockio      :  io to be done in locked regions (.true.)

!  lsplit      :  true - latitudes are shared between a-sets
!                 false - a latitude belongs to only one a-set
!  leq_regions :  true - use new eq_regions partitioning
!                 false - use old NPRGPNS x NPRGPEW partitioning
!  lsplitout   :  output data provided in sequential files (.true.) or
!                 in directories (.false.)
!  limp        :  true: immediate message passing in transposition routines
!  limp_noolap :  true: isend/irecv with no overlap of message passing and 
!                       packing of buffers

INTEGER(KIND=JPIM),ALLOCATABLE:: NPRCIDS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NGPSET2PE(:,:)
LOGICAL :: LSPLIT
LOGICAL :: LEQ_REGIONS
LOGICAL :: LSPLITOUT
LOGICAL :: LOCKIO
LOGICAL :: LIMP
LOGICAL :: LIMP_NOOLAP

INTEGER(KIND=JPIM) :: MP_TYPE
INTEGER(KIND=JPIM) :: MBX_SIZE
INTEGER(KIND=JPIM) :: MYPROC
INTEGER(KIND=JPIM) :: MYSETA
INTEGER(KIND=JPIM) :: MYSETB
INTEGER(KIND=JPIM) :: MYSETW
INTEGER(KIND=JPIM) :: MYSETV
INTEGER(KIND=JPIM) :: MYSETM
INTEGER(KIND=JPIM) :: MYSETN
INTEGER(KIND=JPIM) :: MY_REGION_NS
INTEGER(KIND=JPIM) :: MY_REGION_EW
INTEGER(KIND=JPIM) :: NSTRIN
INTEGER(KIND=JPIM) :: NSTROUT
INTEGER(KIND=JPIM) :: NFLDIN
INTEGER(KIND=JPIM) :: NFLDOUT
INTEGER(KIND=JPIM) :: NSLPAD
INTEGER(KIND=JPIM) :: NINTYPE
INTEGER(KIND=JPIM) :: NOUTTYPE
INTEGER(KIND=JPIM) :: NGATHOUT
INTEGER(KIND=JPIM) :: NWRTOUT
INTEGER(KIND=JPIM) :: NBLKOUT

! ----------------------------------------------------------------------

!*    common block describing the partitioning of data

! ----------------------------------------------------

!  nprocm(0:ncmax) :  gives process which is responsible for Legendre
!             transforms, nmi, and spectral space calculations for a
!             certain wave number m
!  numprocfp(nfprgpg) : gives process which is responsible for FULL-POS
!             horizontal interpolation point. This is only used in
!             FULL-POS.
!  numpp(n_regions_ns) : the number of wave numbers each a-set is responsible
!             for. As aspecial case NUMP = NUMPP(MYSETA).
!  numxpp(n_regions_ns) : Similar to NUMPP() but for NXMAX.
!  nallms(0:max(nsmax,nmsmax)) :  wave numbers for all a-set concate-
!             nated together to give all wave numbers in a-set order.
!             Used when global spectral norms have to be gathered.
!  nptrms(n_regions_ns)  :  pointer to the first wave number of a given a-set
!             in nallms array.
!  mylats(1:ndgenl) if LMESSP else mylats(ndgsag:ndgeng) : mapping
!             between physical latitude number and local latitude number
!             in grid point space on this process. This is identical
!             for all processes within an a-set
!  nptrls(n_regions_ns) : pointer to first global latitude of each a-set
!             for which it performs the Fourier calculations
!  nptrlsf(n_regions_ns) : pointer to first global latitude of each a-set
!             for which it performs the Fourier calculations
!  nfrstlat(n_regions_ns) : first lat of each a-set in grid-point space
!  nfrstloff: offset for first lat of own a-set in grid-point space,
!             i.e. nfrstloff=nfrstlat(my_region_ns)-1
!  nlstlat(n_regions_ns) : last lat of each a-set in grid-point space
!  nptrfrstlat(n_regions_ns) : pointer to the first latitude of each a-set in
!             NSTA and NONL arrays
!  nptrlstlat(n_regions_ns) : pointer to the last latitude of each a-set in
!             NSTA and NONL arrays
!  nptrfloff    : offset for pointer to the first latitude of own a-set
!               NSTA and NONL arrays, i.e. nptrfrstlatf(my_region_ns)-1
!  nptrlat      : pointer to start of latitude in grid-point space
!  lsplitlat(ndglg) : true if latitude is split in grid point space
!              over two a-sets
!  myfrstactlat : first actual lat on this PE in grid-point space,
!                 it is nfrstlat(my_region_ns)
!  mylstactlat  : last actual lat on this PE in grid-point space,
!                 it is nlstlat(my_region_ns)
! ------------------------------------------------------------------
!  nptrsv(nprtrw+1) :  pointer to first spectral wave column to be
!             handled by each b-set. Used for semi-implicit calculations
!             and Jb vertical transforms, and only really if nprtrv>1.
!  nptrcv(nprtrv+1) :  As nptrsv but for ncmax arrays
!  nptrtv(nprtrv+1) :  As nptrsv but for ntmax arrays
!  nptrsvf(nprtrv+1) :  As nptrsv but for the case where full m-columns
!             have to be treated by one processor for the vertical
!             spectral calculations. This is the case if implicit
!             treatment of Coriolis terms is used and in other cases.
!  nptrmf(nprtrv+1)  :  Distribution of m-columns among b-sets used for
!             the full m-column cases where nptrsvf() is used.
!  nspstaf(0:nsmax) : pointer to where each m-column starts (used for
!             the full m-column cases where nptrsvf() is used.
!  numll(nprtrv+1) :  distribution of levels among b-sets for Legendre
!             transforms, FFT and horizontal diffusion.
!             To simplify coding numll(nprtrv+1) is defined to zero.
!  numvmo(nprtrv) :  number of vertical normal modes on each b-set
!  numvmojb(nprtrv) : number of vertical normal modes on each b-set for
!             Jb computations
!  nptrll(nprtrv+1) :  defines the first level treated on each b-set
!             To simplify coding nptrll(nprtrv+1)=nptrll(nprtrv)
!  npsp    :  =1 if surface pressure is handled by this processor for
!             the Legendre Trasforms and FFT calculations. npsp is
!             the same for all processors within a b-set.
!  npsurf(nprtrv)  :  contains the npsp-values for each b-set
!  nbsetlev(nflevg) :  the b-set on which a level belongs. Please use
!              global indexing.
!  nbsetsp :  the b-set on which the surface pressure belongs.
!  mylevs(nflevl) :  mapping between local and global numbering for the
!             levels handled by this process.
!  nvmodist(nvmodmxpp,nprtrv) : normal modes mapped to the different
!             b-sets. The same distribution strategy is used for NMI and
!             Jb calculations. The number of modes is usually larger
!             for Jb caluclations.
!  nspec2v :  number of spectral columns treated by this process for
!             semi-implicit calculations and other vertical transforms
!  ncpec2v :  like nspec2v for NCMAX arrays
!  ntpec2v :  like nspec2v for NTMAX arrays
!  nspec2vf:  number of spectral columns treated by this process for
!             semi-implicit calculations for the full m-column cases.
!             See nptrsvf().
!  nsta(ndgsag:ndgeng+n_regions_ns-1,n_regions_ew) :  Position of first grid column
!             for the latitudes on a processor. The information is
!             available for all processors. The b-sets are distinguished
!             by the last dimension of nsta(). The latitude band for
!             each a-set is addressed by nptrfrstlat(jaset),
!             nptrlstlat(jaset), and nptrfloff=nptrfrstlat(my_region_ns) on
!             this processors a-set. Each split latitude has two entries
!             in nsta(,:) which necessitates the rather complex
!             addressing of nsta(,:) and the overdimensioning of nsta by
!             n_regions_ns.
!  nonl(ndgsag:ndgeng+n_regions_ns-1,n_regions_ew)  :  number of grid columns for
!             the latitudes on a processor. Similar to nsta() in data
!             structure.
!             belong to it in fourier space. Available for all n_regions_ew
!             processors within this processors a-set.
!  napsets :  number of apple sets at the poles. Default is zero.
!  nglobalindex : mapping of local grid points to global grid points
!               : used for debugging
!  nglobalproc  : global data structure containing proc distribution
!                 an ngptotg array that maps owning proc
!  nlocalindex  : global data structure containing local index
!                 an ngptotg array that maps the local index into a
!                 ngptot array for the owning proc

!  -- SLCSET and SLRSET variables (based on NSLWIDE).
!  naslb1  :  local inner dimension of semi-Lagrangian buffer. It is
!             the number of columns in the core+halo region on this
!             processor.
!  nslprocs   : semi-Lagrangian communication :  number of processors
!             this processor needs to communicate with.
!  nslrpt     : the number of columns received from other PE's when
!             computing the halo for interpolations.
!  nslspt     : the number of columns sent to other PE's when
!             computing the halo for interpolations.
!  nslmpbufsz : size of semi-Lagrangian communication buffer in
!             slcomm.F. It is sized so the total requirement is kept
!             below ncombflen.
!  nslsta(ndgsal-nslwide:ndgenl+nslwide)  :  Start position in semi-
!             Lagrangian buffer ZSLBUF1 of grid columns for each local
!             and halo latitude.
!  nslonl(ndgsal-nslwide:ndgenl+nslwide)  :  number of grid columns on
!             each local and halo latitude in the semi-Lagrangian
!             buffer ZSLBUF1. Only used in dm version.
!  nsloff(ndgsal-nslwide:ndgenl+nslwide)  :  offset to the start of each
!             local and halo latitude in the semi-Lagrangian buffer
!             ZSLBUF1. Only used in dm version.
!  nslext(1-ndlon:ndlon+ndlon,1-nslwide:ndgenl+nslwide) in dm version
!  and nslext(nslext(0:ndlon+2,ndgsag:ndgeng) in sm version : pointer
!             that makes sure addressing of points in the east-west
!             extension zone is correct. It also handles the half
!             latitude shift of extension latitudes at the poles.
!             In the sm version this array is just the identity, but
!             used in order to keep sm and dm code in common.
!  nslsendpos: the addresses within the semi-Lagrangian buffer of point sent 
!            from this PE.
!  nslrecvpos: the addresses within the semi-Lagrangian buffer of point 
!            received on this PE.
!  nsendptr  : pointer to the first point for each of the PE's that has to 
!            receive semi-Lagrangian halo-data from this. 
!            Used for addressing nslsendpos().
!  nrecvptr  : pointer to the first point for each of the PE's that are sending
!            semi-Lagrangian halo-data to this PE. 
!            Used for addressing nslrecvpos().
!  nsendnum(nproc+1) : Pointing at the first semi-Lagrangian
!            halo data entry this processor is sending to each of the
!            other processors. The number of columns sent is equal to
!            nsendnum(irecver+1)-nsendnum(irecver), and might be zero.
!  nrecvnum(nproc+1) : Pointing at the first semi-Lagrangian
!            halo data entry this processor is receiving from each of
!            the other processors. The number of columns received is
!            equal to nrecvnum(isender+1)-nrecvnum(isender), it might
!            be zero.
!  nslcore(ngptot) :  Pointer to this processors core region points
!            within the semi-Lagrangian buffer
!  nslcomm(nslprocs)  : semi-Lagrangian communication : list of the
!             processors this proceesor has to communicate with.

!  -- SUFPCSET and SUFPRSET variables (based on NFPWIDE).
!  nafpb1      : FULL-POS version of naslb1
!  nfpprocs    : FULL-POS version of nslprocs
!  nfpmpbufsz  : FULL-POS version of nslmpbufsz
!  nfprpt      : FULL-POS version of nslrpt
!  nfpspt      : FULL-POS version of nslspt
!  nfpsta      : FULL-POS version of nslsta
!  nfponl      : FULL-POS version of nslonl
!  nfpoff      : FULL-POS version of nsloff
!  nfpext      : FULL-POS version of nslext
!  nfpsendpos  : FULL-POS version of nslsendpos
!  nfprecvpos  : FULL-POS version of nslrecvpos
!  nfpsendptr  : FULL-POS version of nsendptr
!  nfprecvptr  : FULL-POS version of nrecvptr
!  nfpcore     : FULL-POS version of nslcore
!  nfpcomm     : FULL-POS version of nslcomm

!   -- SLCSET variables (based on NOBWIDE)
!  nobsta      : observation version of nslsta
!  nobonl      : observation version of nslonl
!  noboff      : observation version of nsloff
 
!  -- SLCSET variables (based on NRIWIDE - model grid).
!  narib1      : Radiation input version of naslb1
!  nriprocs    : Radiation input version of nslprocs
!  nrimpbufsz  : Radiation input version of nslmpbufsz
!  nrirpt      : Radiation input version of nslrpt
!  nrispt      : Radiation input version of nslspt
!  nrista      : Radiation input version of nslsta
!  nrionl      : Radiation input version of nslonl
!  nrioff      : Radiation input version of nsloff
!  nriext      : Radiation input version of nslext
!  nrisendpos  : Radiation input version of nslsendpos
!  nrirecvpos  : Radiation input version of nslrecvpos
!  nrisendptr  : Radiation input version of nsendptr
!  nrirecvptr  : Radiation input version of nrecvptr
!  nricore     : Radiation input version of nslcore
!  nricomm     : Radiation input version of nslcomm

!  -- SLCSET variables (based on NROWIDE - radiation grid).
!  narob1      : Radiation input version of naslb1
!  nroprocs    : Radiation input version of nslprocs
!  nrompbufsz  : Radiation input version of nslmpbufsz
!  nrorpt      : Radiation input version of nslrpt
!  nrospt      : Radiation input version of nslspt
!  nrosta      : Radiation input version of nslsta
!  nroonl      : Radiation input version of nslonl
!  nrooff      : Radiation input version of nsloff
!  nroext      : Radiation input version of nslext
!  nrosendpos  : Radiation input version of nslsendpos
!  nrorecvpos  : Radiation input version of nslrecvpos
!  nrosendptr  : Radiation input version of nsendptr
!  nrorecvptr  : Radiation input version of nrecvptr
!  nrocore     : Radiation input version of nslcore
!  nrocomm     : Radiation input version of nslcomm

! ------------------------------------------------------------------

!  ncombflen : Size of communication buffer. This is the maximum per
!              processor buffer space (in words) that the IFS should use
!              for one or more sends before receives are issued from
!              destination processors.

INTEGER(KIND=JPIM),ALLOCATABLE:: NUMPP(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NUMXPP(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPROCM(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NUMPROCFP(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPTRMS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NALLMS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPTRLS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPTRSV(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPTRCV(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPTRTV(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPTRSVF(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPTRMF(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NSPSTAF(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NUMLL(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPTRLL(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NUMVMO(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NUMVMOJB(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: MYLEVS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPSURF(:)
INTEGER(KIND=JPIM),ALLOCATABLE,TARGET :: NSTA(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE,TARGET :: NONL(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE,TARGET :: NPTRFRSTLAT(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPTRLSTLAT(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NPTRLAT(:)
INTEGER(KIND=JPIM),ALLOCATABLE,TARGET :: NFRSTLAT(:)
INTEGER(KIND=JPIM),ALLOCATABLE,TARGET :: NLSTLAT(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NBSETLEV(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NGLOBALINDEX(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NGLOBALPROC(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NLOCALINDEX(:)

LOGICAL,ALLOCATABLE:: LSPLITLAT(:)

INTEGER(KIND=JPIM),ALLOCATABLE:: MYLATS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NVMODIST(:,:)

!     -- SLCSET and SLRSET variables (based on NSLWIDE).

INTEGER(KIND=JPIM),ALLOCATABLE:: NSLSTA(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NSLONL(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NSLOFF(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NSLEXT(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NSLSENDPOS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NSLRECVPOS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NSENDPTR(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRECVPTR(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NSLCORE(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NSLCOMM(:)

!     -- SUFPCSET and SUFPRSET variables (based on NFPWIDE).

INTEGER(KIND=JPIM),ALLOCATABLE:: NFPSTA(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFPONL(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFPOFF(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFPEXT(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFPSENDPOS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFPRECVPOS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFPSENDPTR(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFPRECVPTR(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFPCORE(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFPCOMM(:)

!     -- SLCSET variables (based on NOBWIDE)

INTEGER(KIND=JPIM),ALLOCATABLE:: NOBSTA(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NOBONL(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NOBOFF(:)

!     -- SLCSET variables (based on NRIWIDE).

INTEGER(KIND=JPIM),ALLOCATABLE:: NRISTA(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRIONL(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRIOFF(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRIEXT(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRISENDPOS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRIRECVPOS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRISENDPTR(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRIRECVPTR(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRICORE(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRICOMM(:)

!     -- SLCSET variables (based on NROWIDE).

INTEGER(KIND=JPIM),ALLOCATABLE:: NROSTA(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NROONL(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NROOFF(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NROEXT(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NROSENDPOS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRORECVPOS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NROSENDPTR(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRORECVPTR(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NROCORE(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NROCOMM(:)

INTEGER(KIND=JPIM) :: NUMXP
INTEGER(KIND=JPIM) :: NPSP
INTEGER(KIND=JPIM) :: NSPEC2V
INTEGER(KIND=JPIM) :: NCPEC2V
INTEGER(KIND=JPIM) :: NTPEC2V
INTEGER(KIND=JPIM) :: NSPEC2VF
INTEGER(KIND=JPIM) :: NBSETSP
INTEGER(KIND=JPIM) :: NFRSTLOFF
INTEGER(KIND=JPIM) :: MYFRSTACTLAT
INTEGER(KIND=JPIM) :: MYLSTACTLAT
INTEGER(KIND=JPIM) :: NAPSETS
INTEGER(KIND=JPIM) :: NPTRFLOFF
INTEGER(KIND=JPIM) :: NCOMBFLEN

!     -- scalar integers depending on NSLWIDE.

INTEGER(KIND=JPIM) :: NASLB1
INTEGER(KIND=JPIM) :: NSLPROCS
INTEGER(KIND=JPIM) :: NSLMPBUFSZ
INTEGER(KIND=JPIM) :: NSLRPT
INTEGER(KIND=JPIM) :: NSLSPT

!     -- scalar integers depending on NFPWIDE.

INTEGER(KIND=JPIM) :: NAFPB1
INTEGER(KIND=JPIM) :: NFPPROCS
INTEGER(KIND=JPIM) :: NFPMPBUFSZ
INTEGER(KIND=JPIM) :: NFPRPT
INTEGER(KIND=JPIM) :: NFPSPT

!     -- scalar integers depending on NRIWIDE.

INTEGER(KIND=JPIM) :: NARIB1
INTEGER(KIND=JPIM) :: NRIPROCS
INTEGER(KIND=JPIM) :: NRIMPBUFSZ
INTEGER(KIND=JPIM) :: NRIRPT
INTEGER(KIND=JPIM) :: NRISPT

!     -- scalar integers depending on NROWIDE.

INTEGER(KIND=JPIM) :: NAROB1
INTEGER(KIND=JPIM) :: NROPROCS
INTEGER(KIND=JPIM) :: NROMPBUFSZ
INTEGER(KIND=JPIM) :: NRORPT
INTEGER(KIND=JPIM) :: NROSPT

! ----------------------------------------------------------------------

!$OMP THREADPRIVATE(leq_regions,limp,limp_noolap,lockio,lsplit,lsplitout,mbx_size,mp_type,my_region_ew,my_region_ns)
!$OMP THREADPRIVATE(myfrstactlat,mylstactlat,myproc,myseta,mysetb,mysetm,mysetn,mysetv,mysetw,nafpb1,napsets,narib1)
!$OMP THREADPRIVATE(narob1,naslb1,nblkout,nbsetsp,ncombflen,ncpec2v,nfldin,nfldout,nfpmpbufsz,nfpprocs,nfprpt,nfpspt)
!$OMP THREADPRIVATE(nfrstloff,ngathout,nintype,nouttype,npsp,nptrfloff,nrimpbufsz,nriprocs,nrirpt,nrispt,nrompbufsz)
!$OMP THREADPRIVATE(nroprocs,nrorpt,nrospt,nslmpbufsz,nslpad,nslprocs,nslrpt,nslspt,nspec2v,nspec2vf,nstrin,nstrout)
!$OMP THREADPRIVATE(ntpec2v,numxp,nwrtout)
!$OMP THREADPRIVATE(lsplitlat,mylats,mylevs,nallms,nbsetlev,nfpcomm,nfpcore,nfpext,nfpoff,nfponl,nfprecvpos,nfprecvptr)
!$OMP THREADPRIVATE(nfpsendpos,nfpsendptr,nfpsta,nfrstlat,nglobalindex,nglobalproc,ngpset2pe,nlocalindex,nlstlat,noboff)
!$OMP THREADPRIVATE(nobonl,nobsta,nonl,nprcids,nprocm,npsurf,nptrcv,nptrfrstlat,nptrlat,nptrll,nptrls,nptrlstlat,nptrmf)
!$OMP THREADPRIVATE(nptrms,nptrsv,nptrsvf,nptrtv,nrecvptr,nricomm,nricore,nriext,nrioff,nrionl,nrirecvpos,nrirecvptr)
!$OMP THREADPRIVATE(nrisendpos,nrisendptr,nrista,nrocomm,nrocore,nroext,nrooff,nroonl,nrorecvpos,nrorecvptr,nrosendpos)
!$OMP THREADPRIVATE(nrosendptr,nrosta,nsendptr,nslcomm,nslcore,nslext,nsloff,nslonl,nslrecvpos,nslsendpos,nslsta)
!$OMP THREADPRIVATE(nspstaf,nsta,numll,numpp,numprocfp,numvmo,numvmojb,numxpp,nvmodist)
END MODULE YOMMP
