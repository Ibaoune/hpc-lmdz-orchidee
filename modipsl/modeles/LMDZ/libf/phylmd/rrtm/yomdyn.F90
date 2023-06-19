MODULE YOMDYN

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -------------------------------------------------------------------------

!*    Control variables for the DYNAMICS

!=========== TIME STEPPING ====================================================

! TSTEP   : length of the timestep in seconds
! TDT     : For leap-frog scheme:
!            2*TSTEP except at the first time step where it is TSTEP
!           For a two-time level scheme (semi-Lagrangian), TDT is always TSTEP.
! REPS1   : timefiltering constant applied to t-1
! REPS2   : timefiltering constant applied to t+1
! REPSM1  : timefiltering constant applied to t-1 (moisture vars.)
! REPSM2  : timefiltering constant applied to t+1 (moisture vars.)
! REPSP1  : timefiltering constant applied to t-1 for all surface fields

REAL(KIND=JPRB) :: TSTEP
REAL(KIND=JPRB) :: TDT
REAL(KIND=JPRB) :: REPS1
REAL(KIND=JPRB) :: REPS2
REAL(KIND=JPRB) :: REPSM1
REAL(KIND=JPRB) :: REPSM2
REAL(KIND=JPRB) :: REPSP1

!====== MAIN HORIZONTAL DIFFUSION SCHEME ======================================

! * CHARACTERISTIC TIMES:
! HDIRVOR  : for diffusion of vorticity.
! HDIRDIV  : for diffusion of divergence.
! HDIRT    : for diffusion of temperature.
! HDIRQ    : for diffusion of humidity.
! HDIRO3   : for diffusion of ozone.
! HDIRPD   : for diffusion of pressure departure (non hydrostatic).
! HDIRVD   : for diffusion of vertical divergence (non hydrostatic).
! HDIRSP   : for diffusion of surface pressure.

! * REVERSE OF CHARACTERISTIC TIMES:
! HRDIRVOR  : for diffusion of vorticity.
! HRDIRDIV  : for diffusion of divergence.
! HRDIRT    : for diffusion of temperature.
! HRDIRQ    : for diffusion of humidity.
! HRDIRO3   : for diffusion of ozone.
! HRDIRPD   : for diffusion of pressure departure (non hydrostatic).
! HRDIRVD   : for diffusion of vertical divergence (non hydrostatic).
! HRDIRSP   : for diffusion of surface pressure.

! RRDXTAU  : overall intensity of HD 
! RDAMPVOR : local enhancing coefficient for diffusion of vorticity.
! RDAMPDIV : local enhancing coefficient for diffusion of divergence.
! RDAMPT   : local enhancing coefficient for diffusion of temperature.
! RDAMPQ   : local enhancing coefficient for diffusion of humidity.
! RDAMPO3  : local enhancing coefficient for diffusion of ozone.
! RDAMPPD  : local enhancing coefficient for diffusion of pressure departure.
! RDAMPVD  : local enhancing coefficient for diffusion of vertical divergence.
! RDAMPSP  : local enhancing coefficient for diffusion of surface pressure.
! LREPHD   : key to reproduce HD consistency:
!            if TRUE the consistency of HDIR[x] is ensured, while namelist
!                    values of RRDAMP[x] can be slightly modified;
!            if FALSE the HD is driven exactly by RRDXTAU and RRDAMP[x]
!                     but consistency of HDIR[x] is not guarranted
! LNEWHD   : only for ECMWF: "new" or "historical" values of HD set-up

! REXPDH   : order of the diffusion 
!            (exponent for the wavenumber dependency).
! FRANDH   : threshold for the wavenumber dependency.
! SLEVDH   : first threshold for the pressure dependency scaled by VP00.
! SLEVDH2  : second threshold for the pressure dependency scaled by VP00.
! SLEVDH3  : third threshold for the pressure dependency scaled by VP00
!            (used to bound the vertical increase of diffusion in the
!            upper stratosphere).
! NSREFDH  : threshold for the truncation dependency.

! * LEVEL AND WAVENUMBER DEPENDENT INVERSE CHARACTERISTIC TIMES:
! RDIVOR   : for diffusion of vorticity.
! RDIDIV   : for diffusion of divergence.
! RDITG    : for diffusion of temperature.
! RDIGFL   : for diffusion of GFL vars.
! RDIPD    : for diffusion of pressure departure (NH).
! RDIVD    : for diffusion of vertical divergence (NH).
! RDISP    : for diffusion of surface pressure.

! GMR      : coefficients for spectral multiplication by GM.
! RDHI     : main horizontal diffusion operator used for stretched ARPEGE.

! LSTRHD   : .T.: main horizontal diffusion operator adapted to stretched ARP.
! HDTIME_STRHD: TDT (if not, the main horizontal diffusion operator
!            used for stretched ARPEGE is recomputed).

REAL(KIND=JPRB) :: HDIRVOR
REAL(KIND=JPRB) :: HDIRDIV
REAL(KIND=JPRB) :: HDIRT
REAL(KIND=JPRB) :: HDIRQ
REAL(KIND=JPRB) :: HDIRO3
REAL(KIND=JPRB) :: HDIRPD
REAL(KIND=JPRB) :: HDIRVD
REAL(KIND=JPRB) :: HDIRSP
REAL(KIND=JPRB) :: HRDIRVOR
REAL(KIND=JPRB) :: HRDIRDIV
REAL(KIND=JPRB) :: HRDIRT
REAL(KIND=JPRB) :: HRDIRQ
REAL(KIND=JPRB) :: HRDIRO3
REAL(KIND=JPRB) :: HRDIRPD
REAL(KIND=JPRB) :: HRDIRVD
REAL(KIND=JPRB) :: HRDIRSP
REAL(KIND=JPRB) :: RRDXTAU
REAL(KIND=JPRB) :: RDAMPVOR
REAL(KIND=JPRB) :: RDAMPDIV
REAL(KIND=JPRB) :: RDAMPT
REAL(KIND=JPRB) :: RDAMPQ
REAL(KIND=JPRB) :: RDAMPO3
REAL(KIND=JPRB) :: RDAMPPD
REAL(KIND=JPRB) :: RDAMPVD
REAL(KIND=JPRB) :: RDAMPSP
LOGICAL :: LREPHD
LOGICAL :: LNEWHD
REAL(KIND=JPRB) :: REXPDH
REAL(KIND=JPRB) :: FRANDH
REAL(KIND=JPRB) :: SLEVDH
REAL(KIND=JPRB) :: SLEVDH2
REAL(KIND=JPRB) :: SLEVDH3
INTEGER(KIND=JPIM) :: NSREFDH
REAL(KIND=JPRB),ALLOCATABLE:: RDIVOR(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDIDIV(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDITG(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDIGFL(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDIPD(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDIVD(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDISP(:)
REAL(KIND=JPRB),ALLOCATABLE:: GMR(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDHI(:,:,:)
LOGICAL :: LSTRHD
REAL(KIND=JPRB) :: HDTIME_STRHD

!====== SEMI-LAGRANGIAN HORIZONTAL DIFFUSION SCHEME (SLHD) ====================

! * FOR SLHD INTERPOLATIONS:
! SLHDA   :    Scaling factor of the deformation in f(d) function
!              (including the model resolution correction)
! SLHDA0  :    Namelist variable allowing to compute SLHDA
!              (scaling factor of the deformation in f(d) function
!              without the model resolution correction)
! SLHDB   :    Exponent of the deformation in f(d) function
! SLHDD0  :    Treshold for deformation tensor enhancement
! ALPHINT :    Limit for the interval of enhancing linear
!              S-L interpolation by smoother (should be
!              within the interval <0,0.5>)
! GAMMAX  :    Maximum value for the Gamma function (the weight
!              of the smoother for the diffusive S-L interpolator),
!              including the timestep correction.
! GAMMAX0 :    Namelist variable allowing to compute GAMMAX
!              (maximum value for the Gamma function,
!              without the timestep correction).
! SLHDKMAX:    Maximum value for the Kappa function

! * THE "HDS" CHARACTERISTIC TIMES (obsolete):
! HDSRVOR : for diffusion of vorticity.
! HDSRDIV : for diffusion of divergence.
! HDSRVD  : for diffusion of vertical divergence (NH).

! * REVERSE OF THE "HDS" CHARACTERISTIC TIMES:
! HRDSRVOR : for diffusion of vorticity.
! HRDSRDIV : for diffusion of divergence.
! HRDSRVD  : for diffusion of vertical divergence (NH).

! RDAMPVORS: local enhancing coefficient for HDS diffusion of vorticity
! RDAMPDIVS: local enhancing coefficient for HDS diffusion of divergence
! RDAMPVDS : local enhancing coefficient for HDS diffusion of vert. divergence
! RDAMPHDS : ratio HRDSRDIV/HRDIRDIV.

! REXPDHS  : order of the diffusion 
!            (exponent for the wavenumber dependency).
! SLEVDHS  : first threshold for the pressure dependency scaled by VP00.
! SLEVDHS2 : second threshold for the pressure dependency scaled by VP00.
! SDRED    : variable modifying the vertical profile based on SLEVDH
!            ( g(l) becomes g(l)-SDRED in the "main" diffusion).

! * "HDS" LEVEL AND WAVENUMBER DEPENDENT INVERSE CHARACTERISTIC TIMES:
! RDSVOR   : for diffusion of vorticity.
! RDSDIV   : for diffusion of divergence.
! RDSVD    : for diffusion of NH vertical divergence variable.
! RDHS     : SLHD additional horizontal diffusion operator used for stretched ARPEGE.

REAL(KIND=JPRB),ALLOCATABLE :: SLHDA(:)
REAL(KIND=JPRB) :: SLHDA0
REAL(KIND=JPRB) :: SLHDB
REAL(KIND=JPRB),ALLOCATABLE :: SLHDD0(:)
REAL(KIND=JPRB) :: ALPHINT
REAL(KIND=JPRB) :: GAMMAX
REAL(KIND=JPRB) :: GAMMAX0
REAL(KIND=JPRB) :: SLHDKMAX
REAL(KIND=JPRB) :: HDSRVOR
REAL(KIND=JPRB) :: HDSRDIV
REAL(KIND=JPRB) :: HDSRVD
REAL(KIND=JPRB) :: HRDSRVOR
REAL(KIND=JPRB) :: HRDSRDIV
REAL(KIND=JPRB) :: HRDSRVD
REAL(KIND=JPRB) :: RDAMPVORS
REAL(KIND=JPRB) :: RDAMPDIVS
REAL(KIND=JPRB) :: RDAMPVDS
REAL(KIND=JPRB) :: RDAMPHDS
REAL(KIND=JPRB) :: REXPDHS
REAL(KIND=JPRB) :: SLEVDHS
REAL(KIND=JPRB) :: SLEVDHS2
REAL(KIND=JPRB) :: SDRED
REAL(KIND=JPRB),ALLOCATABLE:: RDSVOR(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDSDIV(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDSVD(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDHS(:,:,:)

!================== SPECTRAL ENHANCED DIFFUSION ===============================

! LFREIN  : switch to use spectral "enhanced diffusion" (.TRUE. if active)
! LFREINF : same as LFREIN but computed only at STEPO 0 of non-linear run
! LCHDIF  : change diffusion coefficients if LFREINF
! FLCCRI  : critical value of CFL criterion
! RFREIN  : constant for spectral "enhanced diffusion".

LOGICAL :: LFREIN
LOGICAL :: LFREINF
LOGICAL :: LCHDIF
REAL(KIND=JPRB) :: FLCCRI
REAL(KIND=JPRB) :: RFREIN

!======  QUANTITIES TO CHANGE THE VARIABLE IN THE T-EQN =======================

! RCORDIT(NFLEVG)    : correction term at full-levels for diffusion of T.
! RCORDIH(0:NFLEVG)  : correction term at half-levels for SL T-eqn if RCMSMP0/=0
! RCORDIF(NFLEVG)    : correction term at full-levels for SL T-eqn if RCMSMP0/=0

REAL(KIND=JPRB),ALLOCATABLE:: RCORDIT(:)
REAL(KIND=JPRB),ALLOCATABLE:: RCORDIH(:)
REAL(KIND=JPRB),ALLOCATABLE:: RCORDIF(:)

!==== MAXIMUM V-WINDS ALLOWED IN THE SEMI-LAGRANGIAN MODEL ====================

! VMAX1   : if V>VMAX1 (SM) or SQRT(U**2+V**2)>VMAX1 (DM),
!           warning in the SL scheme.
! VMAX2   : if V>VMAX2 (SM) or SQRT(U**2+V**2)>VMAX2 (DM),
!           abort in the SL scheme.

REAL(KIND=JPRB) :: VMAX1
REAL(KIND=JPRB) :: VMAX2

!================== DELTA FORMULATION =========================================

! NDLNPR : NDLNPR=0: conventional formulation of delta, i.e. ln(P(l)/P(l-1)).
!          NDLNPR=1: formulation of delta used in non hydrostatic model,
!                    i.e. (P(l)-P(l-1))/SQRT(P(l)*P(l-1)).

INTEGER(KIND=JPIM) :: NDLNPR

!==== RAYLEIGH FRICTION =======================================================

! RKRF(NFLEVG) : coefficient of Rayleigh friction

REAL(KIND=JPRB),ALLOCATABLE:: RKRF(:)
 
!==== VERTICAL FILTER ========================================================

! LVERFLT : switch to use filter in the vertical 
! REPSVFVO: coefficient for 2-del-eta vertical filter on vorticity
! REPSVFDI: coefficient for 2-del-eta vertical filter on divergence
! NLEVVF  : vertical filter applied for levs 1 to NLEVVF
! LVERAVE_HLUV: switch to filter (vertically) the half-level wind
!               which is computed in routine GPHLUV.

LOGICAL :: LVERFLT
REAL(KIND=JPRB) :: REPSVFVO
REAL(KIND=JPRB) :: REPSVFDI
INTEGER(KIND=JPIM) :: NLEVVF
LOGICAL :: LVERAVE_HLUV

!==== UPPER RADIATIVE BOUNDARY CONDITION ======================================

! RHYDR0 - upper boundary contition for hydrostatic
! RTEMRB - tuning temperature for upper radiative b. c. (LRUBC)
! NRUBC   : control of radiative upper boundary condition :
!           =0 <=> non computation
!           =1 <=> computation on the forecast field
!           =2 <=> computation on the departure of the forecast from the coupling field

REAL(KIND=JPRB) :: RHYDR0
REAL(KIND=JPRB) :: RTEMRB
INTEGER(KIND=JPIM) :: NRUBC

!==== SEMI-IMPLICIT SCHEME, VERTICAL EIGENMODES, PC SCHEMES ===================

! LSIDG   : .F.: Semi-implicit-scheme with reduced divergence.
!           .T.: Semi-implicit scheme with not reduced divergence.

! BETADT  : coefficient for the semi-implicit treatment of divergence,
!           temperature, continuity (and NH if required) equations.
! REFGEO  : reference geopotentiel for shallow-water model.
! SIPR    : reference surface pressure.
! SITR    : reference temperature.
! SITRA   : acoustic reference temperature.
! SITRUB : ref. temper. for SI corr. of temper.(for LRUBC=.T.)
! SIPRUB : coef. for SI corr. of surf. press.  (for LRUBC=.T.)
! SITIME  : =TDT (if not, Helmholtz matrices are recomputed in CNT4).
! SIRPRG  : auxiliary variable for SIGAM,SIGAMA.
! SIRPRN  : auxiliary variable for SITNU,SITNUA
! NSITER  : number of iterations to treat the non linear semi-implicit terms
!           in the non-hydrostatic scheme.
! NCURRENT_ITER : for LNHDYN with PC scheme - current iteration: 
!                   0                 - predictor
!                   1, 2, ..., NSITER - correctors
! LRHDI_LASTITERPC: T (resp. F): when a PC scheme is activated (for example
!  LPC_FULL=.T.), the horizontal diffusion is done at the last iteration
!  of the corrector step (resp. all iterations of the predictor-corrector
!  scheme).

! * PRESSURES LINKED TO A REFERENCE PRESSURE = SIPR
! SIALPH(NFLEVG)  : coefficients "alpha" of hydrostatics.
! SILNPR(NFLEVG)  : Log of ratio of pressures between levels.
! SIDELP(NFLEVG)  : pressure differences across layers.
! SIRDEL(NFLEVG)  : their inverse.
! SITLAH(0:NFLEVG): half-level pressures.
! SITLAF(NFLEVG)  : full-level pressures.
! SIDPHI(NFLEVG)  : geopotential differences across layers.

! SCGMAP((NSMAX+1)*(NSMAX+2)/2,3): coefficients for multiplication by (GM**2)
!                                  in spectral space.
! SIB(NFLEVG,NFLEVG)   : operator "B" of the SI scheme (DIV ===> DP/DT=B.DIV).
! SIMO(NFLEVG,NFLEVG)  : eigenvectors of "B".
! SIMI(NFLEVG,NFLEVG)  : SIMO**-1
! SIVP(NFLEVG)         : eigenvalues of "B".
! SIHEG(NFLEVG,(NSMAX+1)*(NSMAX+2)/2,3), SIHEG2(NFLEVG,NSMAX+1,2:3):
!  Helmholtz operator in case of SI computations with not reduced divergence. 
! SIHEGB(NFLEVG,(NSMAX+1)*(NSMAX+2)/2,3), SIHEGB2(NFLEVG,NSMAX+1,2:3):
!  Additional operators in case of LSIDG=T SI computations in the NH model.
! SITRICA(NSMAX,NFLEVG):   )  coefficients used in tridiagonal solver
! SITRICB(NSMAX,NFLEVG):   )  for the vertically-coupled semi-implicit
! SITRICC(NSMAX,NFLEVG):   )  equations (case LSITRIC=T).

! SIRUB(0:NFLEVG)             :  Kernel of the operator
!                  SIGAM    SITNU
!           (T,ps) -----> P -----> (T,ps)
!  0 is for surface pressure (or its log)
!  1 to NFLEVG is for temperature
!                               t
! S2ETA(NFLEVG)              : S S SIRUB, where S is a Laplacian operator
!    used to eliminate the 2 delta eta wave in the vertical temperature field

! SIFAC : [ 1 - beta**2 (Delta t)**2 C**2 (SITR/SITRA) (LLstar/H**2) ]
!         for NH model.
! SIFACI: [ 1 - beta**2 (Delta t)**2 C**2 (SITR/SITRA) (LLstar/H**2) ]**(-1)
!         for NH model.

! VNORM : constant for new scaling.

LOGICAL :: LSIDG
REAL(KIND=JPRB) :: BETADT
REAL(KIND=JPRB) :: REFGEO
REAL(KIND=JPRB) :: SIPR
REAL(KIND=JPRB) :: SITR
REAL(KIND=JPRB) :: SITRA
REAL(KIND=JPRB) :: SITRUB
REAL(KIND=JPRB) :: SIPRUB
REAL(KIND=JPRB) :: SITIME
REAL(KIND=JPRB) :: SIRPRG
REAL(KIND=JPRB) :: SIRPRN
INTEGER(KIND=JPIM) :: NSITER
INTEGER(KIND=JPIM) :: NCURRENT_ITER
LOGICAL :: LRHDI_LASTITERPC

REAL(KIND=JPRB),ALLOCATABLE:: SIALPH(:)
REAL(KIND=JPRB),ALLOCATABLE:: SILNPR(:)
REAL(KIND=JPRB),ALLOCATABLE:: SIDELP(:)
REAL(KIND=JPRB),ALLOCATABLE:: SIRDEL(:)
REAL(KIND=JPRB),ALLOCATABLE:: SITLAH(:)
REAL(KIND=JPRB),ALLOCATABLE:: SITLAF(:)
REAL(KIND=JPRB),ALLOCATABLE:: SIDPHI(:)
REAL(KIND=JPRB),ALLOCATABLE:: SCGMAP(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIB(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIMO(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIMI(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIVP(:)
REAL(KIND=JPRB),ALLOCATABLE:: SIHEG(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIHEG2(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIHEGB(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIHEGB2(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SITRICA(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SITRICB(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SITRICC(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIRUB(:)
REAL(KIND=JPRB),ALLOCATABLE:: S2ETA(:)
REAL(KIND=JPRB),ALLOCATABLE:: SIFAC(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIFACI(:,:)
REAL(KIND=JPRB) :: VNORM

!=========== SEMI-LAGRANGIAN SWITCHES AND WEIGHTS =============================
!=========== + ADDITIONAL "ADVECTION" SWITCHES ALSO USED IN EULERIAN ========== 

! * Switches NxLAG:
! NVLAG   :  switch for formulation or discretisation of continuity equation.
! NWLAG   :  switch for formulation or discretisation of momentum equations.
! NTLAG   :  switch for formulation or discretisation of temperature equation.
! NSPDLAG :  switch for formulation or discretisation of P-hat equation.
! NSVDLAG :  switch for formulation or discretisation of d-hat equation.
! Remarks about NxLAG:
! a) possible value for NxLAG:
!    NxLAG=1 -> interpolation of R.H.S. of the corresponding eq.
!               to the middle of the trajectory
!    NxLAG=2 -> averaging of R.H.S. of the corresponding eq.
!               along the trajectory with the part corresponding
!               to the departure point added to the t-dt term
!    NxLAG=3 -> averaging of R.H.S. of the corresponding eq.
!               along the trajectory with the part corresponding
!               to the departure point interpolated linearly
! c) For NVLAG and 2D model: 
!    NVLAG>0 stands for the conventional formulation of continuity equation.
!    NVLAG<0 stands for the Lagrangian formulation of continuity equation:
!     in this case the remark a) is valid for ABS(NVLAG).

! * Research of semi-Lagrangian trajectory:
! NITMP   : Number of iterations for computing the medium point of the
!           semi-lagrangian trajectory.
! VETAON  : VETAON*eta(layer nr 1)+(1.-VETAON)*eta(top) is the lower
!           value allowed for ETA of the origin/anterior point in
!           the 3D model.
! VETAOX  : VETAOX*eta(bottom layer)+(1.-VETAOX)*eta(ground) is the
!           upper value allowed for ETA of the origin/anterior point
!           in the 3D model.
! LSETTLS : type of extrapolations needed in the algorithm of trajectory
!           research in the 2TL SL scheme.
!           .F.: linear extrapolations (conventional algorithm).
!           .T.: stable extrapolations combining spatio-temporal extrapolations.
! LELTRA  : if TRUE then use "elegant" algorithm to find departure point
!           (only applicable in 2TL scheme for the shallow-water equations)
! RW2TLFF : when computing the refined position of the origin point for
!           Coriolis term, the new wind used is:
!           0.5*RW2TLFF*(V(F)+V(O)) + (1-RW2TLFF)*V(M)

! * Uncentering factor in the semi-Lagrangian scheme:
! VESL    : first order uncentering factor applied to non linear and linear
!           terms.
! XIDT    : pseudo-second order uncentering factor applied to linear terms,
!           when an alternative second-order averaging is required in the 
!           2TL SL scheme.
! LPC_XIDT: pseudo second order decentering in LPC_FULL PC scheme
!           key used to allocate special buffer for needed quantities
!           to transfer informations from predictor to corrector.

! * Switches for use of quasi-monotone interpolations:
! LQMW    :  Use quasi-monotone three-dimensional interpolations for wind
! LQMHW   :  Use quasi-monotone interpolations in the horizontal for wind
! LQMT    :  Use quasi-monotone three-dimensional interpolations for temperature
! LQMHT   :  Use quasi-monotone interpolations in the horizontal for temperature
! LQMP    :  Use quasi-monotone three-dimensional interpolations for cont. eq
! LQMHP   :  Use quasi-monotone interpolations in the horizontal for cont. eq
! LQMPD   :  Use quasi-monotone three-dimensional interpolations for P-hat eqn.
! LQMHPD  :  Use quasi-monotone interpolations in the horizontal for P-hat eqn.
! LQMVD   :  Use quasi-monotone three-dimensional interpolations for d-hat eqn.
! LQMHVD  :  Use quasi-monotone interpolations in the horizontal for d-hat eqn.

! * Switches for use of spline interpolations:
! LRSPLINE_W   : Use of spline for wind
! LRSPLINE_T   : Use of spline for temperature
! LRSPLINE_P   : Use of spline for continuity equation
! LRSPLINE_SPD : Use of spline for pressure departure
! LRSPLINE_SVD : Use of spline for vertical divergence


! * Treatment of Coriolis term:
! LADVF   : if TRUE then use "advective" treatment of Coriolis terms (SL);
!           in this case 2*Omega*Vec*r is computed analytically.
! LIMPF   : if TRUE then use implicit treatment of Coriolis terms (EUL and SL)
! L2TLFF  : if TRUE then use refined treatment of Coriolis term in 2TLSL scheme
!           (can be currently used also with the 3TL SL vertical interpolating
!           scheme).

! * Change variable with an Eulerian treatment of orography:
! RCMSLP0 : Real for tuning of the Tanguay/Ritchie correction in SL continuity
!           and temperature equations for 3D model.

! * Treatment of MF simplified physics in the semi-Lagrangian TL and AD codes.

! LSL_UNLPHY_F : if TRUE diabatic terms are evaluated at the final point F.
!                if FALSE diabatic terms are evaluated at the orig point O.
!                Remark: this variable is involved only in MF physics.

! * Switch for computation of Moisture Convergence for French deep convection scheme

! NCOMP_CVGQ   :  0 ==> Compute the CVGQ in an Eulerian manner, using spectral
!                       moisture stored in the YQ GFL variable.
!                       In this case YQ must be spectral and
!                       horizontal derivatives are used.
!                 1 ==> Compute the CVGQ in an Eulerian manner, using spectral
!                       moisture stored in the YCVGQ GFL spectral variable and
!                       its horizontal derivatives.
!                       This case is well designed for the case where YQ is
!                       a purely grid-point GFL.
!                 2 ==> Compute the CVGQ in a semi-Lagrangian manner
!                       (Lagrangian tendency - Eulerian tendency), using data
!                       stored in the YCVGQ grid-point variable.
!                       This case is well designed for the case where YQ is
!                       a purely grid-point GFL, and where LSLAG=T.
! remark ky: better to move this variable in SUDYNA/NAMDYNA/YOMDYNA in the
!  future to make it available in SUDIM1 when reading NAMGFL.

INTEGER(KIND=JPIM) :: NVLAG
INTEGER(KIND=JPIM) :: NWLAG
INTEGER(KIND=JPIM) :: NTLAG
INTEGER(KIND=JPIM) :: NSPDLAG
INTEGER(KIND=JPIM) :: NSVDLAG
INTEGER(KIND=JPIM) :: NITMP
REAL(KIND=JPRB) :: VETAON
REAL(KIND=JPRB) :: VETAOX
LOGICAL :: LSETTLS
LOGICAL :: LELTRA
REAL(KIND=JPRB) :: RW2TLFF
REAL(KIND=JPRB) :: VESL
REAL(KIND=JPRB) :: XIDT
LOGICAL ::  LPC_XIDT
LOGICAL :: LQMW
LOGICAL :: LQMHW
LOGICAL :: LQMT
LOGICAL :: LQMHT
LOGICAL :: LQMP
LOGICAL :: LQMHP
LOGICAL :: LQMPD
LOGICAL :: LQMHPD
LOGICAL :: LQMVD
LOGICAL :: LQMHVD
LOGICAL :: LADVF
LOGICAL :: LRSPLINE_W
LOGICAL :: LRSPLINE_T
LOGICAL :: LRSPLINE_P
LOGICAL :: LRSPLINE_SPD
LOGICAL :: LRSPLINE_SVD
LOGICAL :: LIMPF
LOGICAL :: L2TLFF
REAL(KIND=JPRB) :: RCMSLP0
LOGICAL :: LSL_UNLPHY_F
INTEGER(KIND=JPIM) :: NCOMP_CVGQ

!=========== RELAXATION OF THIN LAYER HYPOTHESIS ==============================
! (for more details about "rs", "Ts" see routines gpvcrs.F90 and gpvcts.F90)

! VCPR    : reference pressure (the pressure layer where "rs=a")
! VCTR    : reference temperature (VCTR=Ts(pressure=VCPR))
! VCAK    : coefficient alpha_K used in tha analytic formula of "Ts".
! LADVFW  : as LADVF but for term "-2 Omega vec W k".

REAL(KIND=JPRB) :: VCPR
REAL(KIND=JPRB) :: VCTR
REAL(KIND=JPRB) :: VCAK
LOGICAL :: LADVFW

!     ------------------------------------------------------------------
! LDRY_ECMWF     : .TRUE.  = COMPUTE Cp, R AND R/Cp WITHOUT Q REALTED TERMS
! LDRY_ECMWF     : .FALSE. = COMPUTE Cp, R AND R/Cp WITH    Q REALTED TERMS

LOGICAL :: LDRY_ECMWF

!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(alphint,betadt,flccri,frandh,gammax,gammax0,hdirdiv,hdiro3,hdirpd,hdirq,hdirsp,hdirt)
!$OMP THREADPRIVATE(hdirvd,hdirvor,hdsrdiv,hdsrvd,hdsrvor,hdtime_strhd,hrdirdiv,hrdiro3,hrdirpd,hrdirq)
!$OMP THREADPRIVATE(hrdirsp,hrdirt,hrdirvd,hrdirvor,hrdsrdiv,hrdsrvd,hrdsrvor,l2tlff,ladvf,ladvfw,lchdif)
!$OMP THREADPRIVATE(ldry_ecmwf,leltra,lfrein,lfreinf,limpf,lnewhd,lpc_xidt,lqmhp,lqmhpd,lqmht,lqmhvd,lqmhw)
!$OMP THREADPRIVATE(lqmp,lqmpd,lqmt,lqmvd,lqmw,lrephd,lrhdi_lastiterpc,lrspline_p,lrspline_spd,lrspline_svd)
!$OMP THREADPRIVATE(lrspline_t,lrspline_w,lsettls,lsidg,lsl_unlphy_f,lstrhd,lverave_hluv,lverflt,ncomp_cvgq)
!$OMP THREADPRIVATE(ncurrent_iter,ndlnpr,nitmp,nlevvf,nrubc,nsiter,nspdlag,nsrefdh,nsvdlag,ntlag,nvlag,nwlag)
!$OMP THREADPRIVATE(rcmslp0,rdampdiv,rdampdivs,rdamphds,rdampo3,rdamppd,rdampq,rdampsp,rdampt,rdampvd,rdampvds)
!$OMP THREADPRIVATE(rdampvor,rdampvors,refgeo,reps1,reps2,repsm1,repsm2,repsp1,repsvfdi,repsvfvo,rexpdh,rexpdhs)
!$OMP THREADPRIVATE(rfrein,rhydr0,rrdxtau,rtemrb,rw2tlff,sdred,sipr,siprub,sirprg,sirprn,sitime,sitr,sitra,sitrub)
!$OMP THREADPRIVATE(slevdh,slevdh2,slevdh3,slevdhs,slevdhs2,slhda0,slhdb,slhdkmax,tdt,tstep,vcak,vcpr,vctr,vesl)
!$OMP THREADPRIVATE(vetaon,vetaox,vmax1,vmax2,vnorm,xidt)
!$OMP THREADPRIVATE(gmr,rcordif,rcordih,rcordit,rdhi,rdhs,rdidiv,rdigfl,rdipd,rdisp,rditg,rdivd,rdivor,rdsdiv)
!$OMP THREADPRIVATE(rdsvd,rdsvor,rkrf,s2eta,scgmap,sialph,sib,sidelp,sidphi,sifac,sifaci,siheg,siheg2,sihegb)
!$OMP THREADPRIVATE(sihegb2,silnpr,simi,simo,sirdel,sirub,sitlaf,sitlah,sitrica,sitricb,sitricc,sivp,slhda,slhdd0)
END MODULE YOMDYN
