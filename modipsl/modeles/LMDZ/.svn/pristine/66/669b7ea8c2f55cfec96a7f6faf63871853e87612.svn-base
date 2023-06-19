MODULE YOMCT0

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Control variables for the job - constant within job

!========== TECHNICAL SWITCHES ================================================

! ----- configuration:
! NCONF      : configuration of the job

!                0- 99 : 3-D integration job
!              100-199 : variational job
!              200-299 : 2-D integration job
!              300-349 : KALMAN filter                    (currently obsolete)
!              350-399 : predictability model             (currently obsolete)
!              400-499 : test of the adjoint
!              500-599 : test of the tangent linear model
!              600-699 : eigenvalue/vector solvers
!              700-799 : optimal interpolation
!              800-899 : sensitivity
!              900-999 : miscellaneous other configurations.

!                    1 : 3-D primitive equation model

!                  131 : incremental 4-D VAR/3-D VAR

!                  201 : shallow-water model
!                  202 : vorticity equation model

!                  401 : test of adjoint with 3-D P.E. model
!                  421 : test of adjoint with shallow-water model
!                  422 : test of adjoint with vorticity equation model

!                  501 : test of tangent linear with 3-D P.E. model
!                  521 : test of tangent linear with shallow-water model
!                  522 : test of tangent linear with vorticity equation model

!                  601 : eigenvalue/vector solver for 3-D P.E. model

!                  701 : optimal interpolation with CANARI

!                  801 : sensitivity with 3-D P.E. model

!                  901 : set up initial conditions (CPREP1)
!                  903 : set up initial conditions (CPREP5)
!                  911 : computes dilatation matrices for spectral space
!                  912 : computes rotation matrices for spectral space
!                  923 : initialisation of climatologic files
!                  931 : creation of an ARPEGE file containing the NESDIS SST (INCLITC).
!                  932 : interpolates the sea-ice concentration field from
!                        satellite data to the ARPEGE grid (CSEAICE).
!                  940 : corrects the dry surface pressure in an ARPEGE file
!                        to keep total mass constant compared to a reference
!                        ARPEGE file (CORMASS)
!                  951 : difference between 2 model states (CPREP2)

! ----- ???:
! NTASKS     : ??? (no comment provided, currently used in CANARI only).

! ----- variables linked with diagnostics and outputs:
! LALLOPR    : .T. = print information about all allocations/deallocations
! NPRINTLEV  : 0 = "basic prints'; 1 = "more prints"; 2 = "debug prints"

! ----- type of file used:
! LFDBOP     : .T. = fields data base utilized
! LGRBOP     : .T. = output in GRIB (not ARPEGE)
! LFBDAP     : .T. = diagnostics written on trajectory files
! LARPEGEF   : .T. = use ARPEGE files
! LARPEGEF_TRAJHR   : .T. = use ARPEGE files for high resolution trajectory
! LARPEGEF_TRAJBG   : .T. = use ARPEGE files for background
! LARPEGEF_RDGP_INIT     : .T. = use grid-point ARPEGE files
! LARPEGEF_RDGP_TRAJHR   : .T. = use grid-point ARPEGE files for HR trajectory
! LARPEGEF_RDGP_TRAJBG   : .T. = use grid-point ARPEGE files for background
! CNDISPP    : directory of display files

! ----- variables linked to post-processing:
! LFPOS      : .T. = use of Full-POS rather than POS
! CFPNCF     : name of Full-POS control file (pseudo-configuration 927)
! LFPART2    : .T. = second part of interpolations for changing geometry
!              (pseudo-configuration 927)
! CFDIRLST   : path of postprocessing listing file
! CNPPATH    : directory of postprocessing namelist files

! ----- use of transmission coefficients stored in Fourier or spectral space:
! LRETCFOU   : .T. = reading Fourier transmission coefficients on file
!              for simplified physics.
! LWRTCFOU   : .T. = writing Fourier transmission coefficients on file.
! LRFOUTCNORM: activates diagnostics on Fourier transmission coefficients.
!              - LRFOUTCNORM(1)=.T.: global statistics.
!              - LRFOUTCNORM(2)=.T.: statistics per latitude.
!              - LRFOUTCNORM(3)=.T.: statistics per layer.
!              - LRFOUTCNORM(4)=.T.: statistics per wavenumber.
! LRGPTCNORM : activates diagnostics on grid-point transmission coefficients.
!              - LRGPTCNORM(1)=.T.: global statistics.
!              - LRGPTCNORM(2)=.T.: statistics per latitude.
!              - LRGPTCNORM(3)=.T.: statistics per layer.

! ----- other variables:
! LNF        : .T. = start, .F. = restart
! LSMSSIG    : .T. = send signals to SMS (ECMWF supervisor)
! LOPDIS     : .T. = calls OPDIS
! NCYCLE     : cycle identifier
! CNMEXP     : name of the experiment
!              An experiment is identified by its name (16 characters)
!              and its cycle (typically same experiment but without a bug)
! CFCLASS    : class for use by FDB
! CTYPE      : type for use by FDB
! LBACKG     : ???
! LMINIM     : ???

!========== MODEL SWITCHES ====================================================
! Remark: this section containes some dynamical variables which cannot be
!         put in NAMDYN because they are used in routines called before SUDYN
!         and sometimes to allocate some arrays.

! ----- advection scheme, semi-Lagrangian scheme:
! LSLAG     : .TRUE. = semi-lagrangian on
! LTWOTL    : .TRUE. if two-time-level semi-Lagrangian scheme will be used
! LRFRIC    : .TRUE. = Rayleigh friction in zonal wind (for p < 9.9hPa)
! LVERCOR   : .T./.F.: thin layer hypothesis relaxed/applied.

! ----- vertical discretisation, vertical boundaries:
! LAPRXPK   : way of computing full-levels pressures in primitive equation
!             hydrostatic model.
!             .T.: full levels are computed by PK=(PK+1/2 + PK-1/2)*0.5
!             .F.: full levels are computed by a more complicated formula
!                  consistent with "alpha" in geopotential formula.
! LREGETA   : .T.: for the interlayer L, ETA(L)=L/NFLEVG
!             .F.: for the interlayer L, ETA(L)=A(L)/P0+B(L)
! LRUBC     : .T.: if radiational upper boundary condition

! ----- numbering of timesteps:
! NSTART    : first timestep of model
! NSTOP     : last timestep of model

! ----- quadrature:
! NQUAD     : 1 ====> GAUSS

! ----- type of equations:
! LNHDYN    : .T. if 3-D non-hydrostatic dynamics is active

! ----- way of computing initial state:
! N2DINI    : 1 initialization for 2D with Haurwitz wave
!           : 2 initialization for 2D with real fields
! N3DINI    : 1 initialization for 3D with standard atmosphere
!             (not available since cycle 12)
!           : 0 initialization for 3D with real fields

! ----- semi-implicit scheme:
! LSITRIC   : .T.: if "tridiagonal" solver is used for the vertically
!             coupled semi-implicit equations

! ----- control of variables which are transformed in spectral space:
! LSPRT     : .T.: if R*T/Rd "virtual temperature" as spectral variable

! ----- diagnostics and frequencies:
! NSPPR     : 0: no spectrum printed in spnorm; only global norms averaged
!             on the vertical are printed
!           : 1: no spectrum printed in spnorm; only global norms averaged
!             on the vertical and global norms on each layer are printed
!           : 2: both total wavenumber spectrum and zonal wavenumber spectrum
!             are printed
! NFRPOS    : frequency of post-processing events (time-steps)
! NFRISP    : frequency of isp (animation !) events (time-steps)
! NFRCO     : frequency of coupled fields (time-steps)
! NFRCORM   : mass correction frequency (>0 time-steps, <0 hours, 0 no
!             correction)
! NFRHIS    : frequency of history write_ups (time-steps)
! NFRMASSCON: frequency of mass conservation fixup (time-steps)
! NFRGDI    : frequency of grid-point space diagnostics
! NFRSDI    : frequency of spectral space diagnostics
! NFRDHFG   : write-up frequency of global DDH
! NFRDHFZ   : write-up frequency of zonal DDH
! NFRDHFD   : write-up frequency of "limited domain" DDH
! NFRDHP    : write-up frequency of DDH files
! NPOSTS    : array containing postprocessing steps
! NPISPS    : array containing isp (animation !) steps
! NHISTS    : array containing history write-up steps
! NMASSCONS : array containing mass conservation fixup steps
! NGDITS    : array containing grid point diagnostics steps
! NSDITS    : array containing spectral diagnostics steps
! NDHFGTS   : array containing write out steps for global DDH
! NDHFZTS   : array containing write out steps for zonal means DDH
! NDHFDTS   : array containing write out steps for limited areas DDH
! NDHPTS    : array containing write out steps for DDH
! Explanation for N[XXX]TS:
!             1) if N[XXX]TS(0)=0 action if MOD(JSTEP,NFR[XXX])=0
!             2) if N[XXX]TS(0)>0 N[XXX]TS(0) significant numbers in 
!                N[XXX]TS are then considered and:
!                action for JSTEP=N[XXX]TS(.)*NFR[XXX]
!             3) IF N[XXX]TS(0)<0
!                action for JSTEP=(N[XXX]TS(.)/DELTAT)*NFR[XXX]

!========== ASSIMILATION SWITCHES =============================================

! LNOBGON : .F. if eq. .T. no stop if unsuccessful reading of CMA files
! LCANARI : .T. = term to control French OI
! LCASIG  : .T. = interpolation of the model errors (French OI)
! LGUESS  : .T. = term of first guess included
! LOBS    : .T. = term of observations included
! LSIMOB  : .T. = if simulated observations
! LOBSC1  : .T. = term of observations included in configuration 1
! LSCREEN : .T. = observation screening for variational assimilation
! LSCREEN_OPENMP : .T. = 4DVAR screening runs in OpenMP-parallel mode over timeslots
! L_SPLIT_SCREEN .T. = to split screenng
! L_SCREEN_CALL: .T. = call to screening routine SCREEN
! LOBSREF : .T. = comparison to observation for the trajectory (NCONF=131)
! LIFSMIN : .T. = if running minimisation
! LIFSTRAJ: .T. = if running high resolution trajectory integration
! NCNTVAR : Definition of the control variable of a variational job.
!           = 1 ===> control variables are model variables
!           = 2 ===> control variables are normalized departures of
!                    model variables from the background field
!           = 3 ===> ..........
! LOLDPP  : .T. use "old" p.p. of T,Q,U and V
! NSTEPINI: Initial step in hours for the initial conditions
!           at the beginning of 4D-Var trajectory (usually 3 hours).
!           It is used to update the step while saving the FCs along
!           the first trajectory.
! NINTERPTRAJ : Interpolation method applied to increments
! NINTERPINCR : Interpolation method applied to increments
!           = 1 ===> Bi-linear interpolation (default)
!           = 2 ===> Bi-cubic Rinterpolation
!           = 3 ===> Conserving interpolation
! For Conserving interpolation the style is defined by:
! NINTERPTRAJORDER : order of interpolation for trajectories
! NINTERPINCRORDER : order of interpolation for increments
! NINTERPTRAJLIMIT : 0 if limiter is NOT used (default), 1 if used for trajectory interpolation
! NINTERPINCRLIMIT : 0 if limiter is NOT used (default), 1 if used for increment interpolation

!========== ALADIN SWITCHES ===================================================

! LELAM   : .T. = limited area model with coupling or fully biperiodic
!                 model (torus)
!           .F. = global model (sphere)
! LRPLANE : .T. = plane cartesian geometry
!           .F. = spherical geometry
! LTENC   : TRUE if tendency coupling of surface pressure
!           FALSE if no tend. coupling of surf. pressure
! LALLTC  : used together with LTENC when LTENC=.T.
!          - no meaning for quadratic tendency coupling, where just t1 coupling
!          is applied at every NEFRCL time step
!          - for lin. tendency coupling:
!          TRUE if tendency coupling of surf. pres. at every step
!          FALSE if tend. coupl., except at every NEFRCL time steps
!          when just t1 coupling
! RTENC : multiplier of EALFA in the tendency coupling scheme
!  for stability reasons (RTENC<=1. close to 1)

!========== AROME SWITCH ======================================================

! LAROME  : .T. = AROME limited area model

!========== ECMWF Single Column Model =========================================
! LSCMEC  : .T. = ECMWF Single Column Model
! LSFCFLX : .T. = forcing with surface fluxes (latent and sensible).
! REXTSHF : externally supplied sensible heat flux [W/m^2]
! REXTLHF : externally supplied latent   heat flux [W/m^2]
! LROUGH  : .T. = surface roughness length is externally specified 
! REXTZ0M : externally supplied roughness length for momentum [m]
! REXTZ0H : externally supplied roughness length for heat [m]

!========== DISTRIBUTED MEMORY SWITCHES =======================================

! LMPOFF  : .T. = switch off the message passing library initialisation
!                 requested for special cases LMESSP=.T. and NPROC=1
!         : .F. = (default) full message passing features
! NPROC   : Total number of processors requested for this run
! N_REGIONS_NS : Number of regions (LEQ_REGIONS=T) or NPRGPNS (LEQ_REGIONS=F)
! N_REGIONS_EW : Maximum number of partitions for all regions
! N_REGIONS : Number of partitions in each region
! NPRGPNS : Number of processors used during grid-point phase in North-South
!           direction (previously known as NPROCA)
! NPRGPEW : Number of processors used during grid-point phase in East-West
!           direction (previously known as NPROCB)
! NPRTRNS : Number of processors used during transform phase in North-South
!           direction (previously the same as NPROCA and now implemented
!           such that NPRTRNS=NPTRW)
! NPRTRM  : Number of processors used during spectral computations in
!           meridional wave direction (previously the same as NPROCA
!           and now implemented such that NPRTRM=NPTRW)
! NPRTRN  : Number of processors used during spectral computations in total
!           wave direction (previously the same as NPROCA and now implemented
!           such that NPRTRN=NPTRV)
! NPRTRW  : Number of processors used during transform phase in wave space
!           (previously the same as NPROCA)
! NPRTRV  : Number of processors used during transform phase in vertical
!           direction (previously known as NPROCB)
! LOUTPUT : .T. = diagnostic output requested on this PE
! NOUTPUT : 0 = No diagnostic output
!           1 = Only diagnostic output from PE1 ( default )
!           2 = Diagnostic output from all PEs into separate files
! LREFOUT : .T. compare to reference run
! LREFGEN : .T. to generate reference file
! LMPDIAG : .T. = extensive message passing diagnostic output requested

!========== PC SCHEMES ========================================================
!----------------------------------------------
! PC SCHEMES QUANTITIES CONSTANT DURING INTEGRATION
!----------------------------------------------

! LPC_FULL - full PC scheme switch (with reiterations of trajectories)

! LPC_OLD - iterative scheme without trajectories recalculations
!         - with iteration of D3 terms (Bubnova et all., 1995)

! LPC_NESC - non-extrapolating two-time level SL SI scheme
!          - X(t+dt/2) = X(t) during predictor

!========== FORCING SWITCH ====================================================
! LSFORC - switch to activate the large scale forcings in setup and cpg
!==============================================================================

! * Parameters:
INTEGER(KIND=JPIM), PARAMETER :: JPNPST=240

! * Technical switches:
INTEGER(KIND=JPIM) :: NCONF
INTEGER(KIND=JPIM) :: NTASKS ! ??????
LOGICAL :: LALLOPR
INTEGER(KIND=JPIM) :: NPRINTLEV
LOGICAL :: LFDBOP
LOGICAL :: LGRBOP
LOGICAL :: LFBDAP
LOGICAL :: LARPEGEF
LOGICAL :: LARPEGEF_TRAJHR
LOGICAL :: LARPEGEF_TRAJBG
LOGICAL :: LARPEGEF_RDGP_INIT
LOGICAL :: LARPEGEF_RDGP_TRAJHR
LOGICAL :: LARPEGEF_RDGP_TRAJBG
CHARACTER (LEN = 120) ::  CNDISPP
LOGICAL :: LFPOS
CHARACTER (LEN = 6) ::  CFPNCF
LOGICAL :: LFPART2
CHARACTER (LEN = 120) ::  CFDIRLST
CHARACTER (LEN = 120) ::  CNPPATH
LOGICAL :: LRETCFOU
LOGICAL :: LWRTCFOU
LOGICAL :: LRFOUTCNORM(4)
LOGICAL :: LRGPTCNORM(3)
LOGICAL :: LNF
LOGICAL :: LSMSSIG
LOGICAL :: LOPDIS
INTEGER(KIND=JPIM) :: NCYCLE=128
CHARACTER (LEN = 16) ::  CNMEXP
CHARACTER (LEN = 2) ::  CFCLASS
CHARACTER (LEN = 2) ::  CTYPE
LOGICAL :: LBACKG
LOGICAL :: LMINIM

! * Model switches:
LOGICAL :: LSLAG
LOGICAL :: LTWOTL
LOGICAL :: LRFRIC
LOGICAL :: LVERCOR
LOGICAL :: LAPRXPK
LOGICAL :: LREGETA
LOGICAL :: LRUBC
INTEGER(KIND=JPIM) :: NSTART
INTEGER(KIND=JPIM) :: NSTOP
INTEGER(KIND=JPIM) :: NQUAD
LOGICAL :: LNHDYN
INTEGER(KIND=JPIM) :: N2DINI
INTEGER(KIND=JPIM) :: N3DINI
LOGICAL :: LSITRIC
LOGICAL :: LSPRT
INTEGER(KIND=JPIM) :: NSPPR
INTEGER(KIND=JPIM) :: NFRPOS
INTEGER(KIND=JPIM) :: NFRISP
INTEGER(KIND=JPIM) :: NFRCORM
INTEGER(KIND=JPIM) :: NFRCO
INTEGER(KIND=JPIM) :: NFRHIS
INTEGER(KIND=JPIM) :: NFRMASSCON
INTEGER(KIND=JPIM) :: NFRGDI
INTEGER(KIND=JPIM) :: NFRSDI
INTEGER(KIND=JPIM) :: NFRDHFG
INTEGER(KIND=JPIM) :: NFRDHFZ
INTEGER(KIND=JPIM) :: NFRDHFD
INTEGER(KIND=JPIM) :: NFRDHP
INTEGER(KIND=JPIM) :: NPOSTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NPISPS(0:JPNPST)
INTEGER(KIND=JPIM) :: NHISTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NMASSCONS(0:JPNPST)
INTEGER(KIND=JPIM) :: NGDITS(0:JPNPST)
INTEGER(KIND=JPIM) :: NSDITS(0:JPNPST)
INTEGER(KIND=JPIM) :: NDHFGTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NDHFZTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NDHFDTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NDHPTS(0:JPNPST)

! * Assimilation:
LOGICAL :: LNOBGON
LOGICAL :: LCANARI
LOGICAL :: LCASIG
LOGICAL :: LGUESS
LOGICAL :: LOBS
LOGICAL :: LSIMOB
LOGICAL :: LOBSC1
LOGICAL :: LSCREEN
LOGICAL :: LSCREEN_OPENMP
LOGICAL :: L_SPLIT_SCREEN
LOGICAL :: L_SCREEN_CALL
LOGICAL :: LOBSREF
LOGICAL :: LIFSMIN
LOGICAL :: LIFSTRAJ
INTEGER(KIND=JPIM) :: NCNTVAR
LOGICAL :: LOLDPP
INTEGER(KIND=JPIM) :: NSTEPINI
INTEGER(KIND=JPIM) :: NINTERPTRAJ
INTEGER(KIND=JPIM) :: NINTERPINCR
INTEGER(KIND=JPIM) :: NINTERPTRAJLIMIT
INTEGER(KIND=JPIM) :: NINTERPINCRLIMIT
INTEGER(KIND=JPIM) :: NINTERPTRAJORDER
INTEGER(KIND=JPIM) :: NINTERPINCRORDER

! * ALADIN:
LOGICAL :: LELAM
LOGICAL :: LRPLANE
LOGICAL :: LTENC
LOGICAL :: LALLTC
REAL(KIND=JPRB) :: RTENC

! * AROME:
LOGICAL :: LAROME

! * ECMWF Single Column Model:
LOGICAL :: LSCMEC
LOGICAL :: LSFCFLX
REAL(KIND=JPRB) :: REXTSHF
REAL(KIND=JPRB) :: REXTLHF
LOGICAL :: LROUGH
REAL(KIND=JPRB) :: REXTZ0M
REAL(KIND=JPRB) :: REXTZ0H

! * Distributed memory:
LOGICAL :: LMPOFF
INTEGER(KIND=JPIM) :: NPROC
INTEGER(KIND=JPIM) :: N_REGIONS_NS
INTEGER(KIND=JPIM) :: N_REGIONS_EW
INTEGER(KIND=JPIM),ALLOCATABLE :: N_REGIONS(:)
INTEGER(KIND=JPIM) :: NPRGPNS
INTEGER(KIND=JPIM) :: NPRGPEW
INTEGER(KIND=JPIM) :: NPRTRNS
INTEGER(KIND=JPIM) :: NPRTRM
INTEGER(KIND=JPIM) :: NPRTRN
INTEGER(KIND=JPIM) :: NPRTRW
INTEGER(KIND=JPIM) :: NPRTRV
LOGICAL :: LOUTPUT
INTEGER(KIND=JPIM) :: NOUTPUT
LOGICAL :: LREFOUT
LOGICAL :: LREFGEN
LOGICAL :: LMPDIAG

! * PC schemes.
LOGICAL ::  LPC_FULL
LOGICAL ::  LPC_OLD
LOGICAL ::  LPC_NESC
!
! * FORCING
LOGICAL ::  LSFORC
!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(cfclass,cfdirlst,cfpncf,cndispp,cnmexp,cnppath,ctype,l_screen_call,l_split_screen,lallopr,lalltc)
!$OMP THREADPRIVATE(laprxpk,larome,larpegef,larpegef_rdgp_init,larpegef_rdgp_trajbg,larpegef_rdgp_trajhr)
!$OMP THREADPRIVATE(larpegef_trajbg,larpegef_trajhr,lbackg,lcanari,lcasig,lelam,lfbdap,lfdbop,lfpart2,lfpos)
!$OMP THREADPRIVATE(lgrbop,lguess,lifsmin,lifstraj,lminim,lmpdiag,lmpoff,lnf,lnhdyn,lnobgon,lobs,lobsc1,lobsref)
!$OMP THREADPRIVATE(loldpp,lopdis,loutput,lpc_full,lpc_nesc,lpc_old,lrefgen,lrefout,lregeta,lretcfou,lrfoutcnorm)
!$OMP THREADPRIVATE(lrfric,lrgptcnorm,lrough,lrplane,lrubc,lscmec,lscreen,lscreen_openmp,lsfcflx,lsforc,lsimob)
!$OMP THREADPRIVATE(lsitric,lslag,lsmssig,lsprt,ltenc,ltwotl,lvercor,lwrtcfou,n2dini,n3dini,n_regions_ew,n_regions_ns)
!$OMP THREADPRIVATE(ncntvar,nconf,ncycle,ndhfdts,ndhfgts,ndhfzts,ndhpts,nfrco,nfrcorm,nfrdhfd,nfrdhfg,nfrdhfz,nfrdhp)
!$OMP THREADPRIVATE(nfrgdi,nfrhis,nfrisp,nfrmasscon,nfrpos,nfrsdi,ngdits,nhists,ninterpincr,ninterpincrlimit)
!$OMP THREADPRIVATE(ninterpincrorder,ninterptraj,ninterptrajlimit,ninterptrajorder,nmasscons,noutput,npisps,nposts)
!$OMP THREADPRIVATE(nprgpew,nprgpns,nprintlev,nproc,nprtrm,nprtrn,nprtrns,nprtrv,nprtrw,nquad,nsdits,nsppr,nstart)
!$OMP THREADPRIVATE(nstepini,nstop,ntasks,rextlhf,rextshf,rextz0h,rextz0m,rtenc)
!$OMP THREADPRIVATE(n_regions)
END MODULE YOMCT0
