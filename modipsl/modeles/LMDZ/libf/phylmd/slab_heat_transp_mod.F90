!
MODULE slab_heat_transp_mod
!
! Slab ocean : temperature tendencies due to horizontal diffusion 
! and / or Ekman transport

USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, klon_glo
IMPLICIT NONE

  ! Variables copied over from dyn3d dynamics:
  REAL,SAVE,ALLOCATABLE :: fext(:) ! Coriolis f times cell area
  !$OMP THREADPRIVATE(fext)
  REAL,SAVE,ALLOCATABLE :: beta(:) ! df/dy
  !$OMP THREADPRIVATE(beta)
  REAL,SAVE,ALLOCATABLE :: unsairez(:) ! 1/cell area
  !$OMP THREADPRIVATE(unsairez)
  REAL,SAVE,ALLOCATABLE :: unsaire(:)
  !$OMP THREADPRIVATE(unsaire)
  REAL,SAVE,ALLOCATABLE :: cu(:) ! cell longitude dim (m)
  !$OMP THREADPRIVATE(cu)
  REAL,SAVE,ALLOCATABLE :: cv(:) ! cell latitude dim (m)
  !$OMP THREADPRIVATE(cv)
  REAL,SAVE,ALLOCATABLE :: cuvsurcv(:) ! cu/cv (v points)
  !$OMP THREADPRIVATE(cuvsurcv)
  REAL,SAVE,ALLOCATABLE :: cvusurcu(:) ! cv/cu (u points)
  !$OMP THREADPRIVATE(cvusurcu)
  REAL,SAVE,ALLOCATABLE :: aire(:) ! cell area
  !$OMP THREADPRIVATE(aire)
  REAL,SAVE :: apoln ! area of north pole points
  !$OMP THREADPRIVATE(apoln)
  REAL,SAVE :: apols ! area of south pole points
  !$OMP THREADPRIVATE(apols)
  REAL,SAVE,ALLOCATABLE :: aireu(:) ! area of u cells 
  !$OMP THREADPRIVATE(aireu)
  REAL,SAVE,ALLOCATABLE :: airev(:) ! area of v cells
  !$OMP THREADPRIVATE(airev)

  ! Local parameters for slab transport
  LOGICAL,SAVE :: alpha_var ! variable coef for deep temp (1 layer)
  !$OMP THREADPRIVATE(alpha_var)
  LOGICAL,SAVE :: slab_upstream ! upstream scheme ? (1 layer)
  !$OMP THREADPRIVATE(slab_upstream)
  LOGICAL,SAVE :: slab_sverdrup ! use wind stress curl at equator
  !$OMP THREADPRIVATE(slab_sverdrup)
  LOGICAL,SAVE :: slab_tyeq ! use merid wind stress at equator
  !$OMP THREADPRIVATE(slab_tyeq)
  LOGICAL,SAVE :: ekman_zonadv ! use zonal advection by Ekman currents
  !$OMP THREADPRIVATE(ekman_zonadv)
  LOGICAL,SAVE :: ekman_zonavg ! zonally average wind stress
  !$OMP THREADPRIVATE(ekman_zonavg)

  REAL,SAVE :: alpham
  !$OMP THREADPRIVATE(alpham)
  REAL,SAVE :: gmkappa
  !$OMP THREADPRIVATE(gmkappa)
  REAL,SAVE :: gm_smax
  !$OMP THREADPRIVATE(gm_smax)

! geometry variables : f, beta, mask...
  REAL,SAVE,ALLOCATABLE :: zmasqu(:) ! continent mask for zonal mass flux
  !$OMP THREADPRIVATE(zmasqu)
  REAL,SAVE,ALLOCATABLE :: zmasqv(:) ! continent mask for merid mass flux
  !$OMP THREADPRIVATE(zmasqv)
  REAL,SAVE,ALLOCATABLE :: unsfv(:) ! 1/f, v points
  !$OMP THREADPRIVATE(unsfv)
  REAL,SAVE,ALLOCATABLE :: unsbv(:) ! 1/beta
  !$OMP THREADPRIVATE(unsbv)
  REAL,SAVE,ALLOCATABLE :: unsev(:) ! 1/epsilon (drag)
  !$OMP THREADPRIVATE(unsev)
  REAL,SAVE,ALLOCATABLE :: unsfu(:) ! 1/F, u points
  !$OMP THREADPRIVATE(unsfu)
  REAL,SAVE,ALLOCATABLE :: unseu(:)
  !$OMP THREADPRIVATE(unseu)

  ! Routines from dyn3d, valid on global dynamics grid only:
  PRIVATE :: gr_fi_dyn, gr_dyn_fi ! to go between 1D nd 2D horiz grid
  PRIVATE :: gr_scal_v,gr_v_scal,gr_scal_u ! change on 2D grid U,V, T points
  PRIVATE :: grad,diverg

CONTAINS
 
  SUBROUTINE ini_slab_transp_geom(ip1jm,ip1jmp1,unsairez_,fext_,unsaire_,&
                                  cu_,cuvsurcv_,cv_,cvusurcu_, &
                                  aire_,apoln_,apols_, &
                                  aireu_,airev_,rlatv, rad, omeg)
    ! number of points in lon, lat
    IMPLICIT NONE
    ! Routine copies some geometry variables from the dynamical core
    ! see global vars for meaning
    INTEGER,INTENT(IN) :: ip1jm
    INTEGER,INTENT(IN) :: ip1jmp1
    REAL,INTENT(IN) :: unsairez_(ip1jm)
    REAL,INTENT(IN) :: fext_(ip1jm)
    REAL,INTENT(IN) :: unsaire_(ip1jmp1)
    REAL,INTENT(IN) :: cu_(ip1jmp1)
    REAL,INTENT(IN) :: cuvsurcv_(ip1jm)
    REAL,INTENT(IN) :: cv_(ip1jm)
    REAL,INTENT(IN) :: cvusurcu_(ip1jmp1)
    REAL,INTENT(IN) :: aire_(ip1jmp1)
    REAL,INTENT(IN) :: apoln_
    REAL,INTENT(IN) :: apols_
    REAL,INTENT(IN) :: aireu_(ip1jmp1)
    REAL,INTENT(IN) :: airev_(ip1jm)
    REAL,INTENT(IN) :: rlatv(nbp_lat-1)
    REAL,INTENT(IN) :: rad
    REAL,INTENT(IN) :: omeg

    CHARACTER (len = 20) :: modname = 'slab_heat_transp'
    CHARACTER (len = 80) :: abort_message
    
    ! Sanity check on dimensions
    if ((ip1jm.ne.((nbp_lon+1)*(nbp_lat-1))).or. &
        (ip1jmp1.ne.((nbp_lon+1)*nbp_lat))) then
      abort_message="ini_slab_transp_geom Error: wrong array sizes"
      CALL abort_physic(modname,abort_message,1)
    endif
! Allocations could be done only on master process/thread...
    allocate(unsairez(ip1jm))
    unsairez(:)=unsairez_(:)
    allocate(fext(ip1jm))
    fext(:)=fext_(:)
    allocate(unsaire(ip1jmp1))
    unsaire(:)=unsaire_(:)
    allocate(cu(ip1jmp1))
    cu(:)=cu_(:)
    allocate(cuvsurcv(ip1jm))
    cuvsurcv(:)=cuvsurcv_(:)
    allocate(cv(ip1jm))
    cv(:)=cv_(:)
    allocate(cvusurcu(ip1jmp1))
    cvusurcu(:)=cvusurcu_(:)
    allocate(aire(ip1jmp1))
    aire(:)=aire_(:)
    apoln=apoln_
    apols=apols_
    allocate(aireu(ip1jmp1))
    aireu(:)=aireu_(:)
    allocate(airev(ip1jm))
    airev(:)=airev_(:) 
    allocate(beta(nbp_lat-1))
    beta(:)=2*omeg*cos(rlatv(:))/rad

  END SUBROUTINE ini_slab_transp_geom

  SUBROUTINE ini_slab_transp(zmasq)

!    USE ioipsl_getin_p_mod, only: getin_p
    USE IOIPSL, ONLY : getin 
    IMPLICIT NONE

    REAL zmasq(klon_glo) ! ocean / continent mask, 1=continent
    REAL zmasq_2d((nbp_lon+1)*nbp_lat)
    REAL ff((nbp_lon+1)*(nbp_lat-1)) ! Coriolis parameter
    REAL eps ! epsilon friction timescale (s-1)
    INTEGER :: slab_ekman 
    INTEGER i
    INTEGER :: iim,iip1,jjp1,ip1jm,ip1jmp1

! Some definition for grid size
    ip1jm=(nbp_lon+1)*(nbp_lat-1)
    ip1jmp1=(nbp_lon+1)*nbp_lat
    iim=nbp_lon
    iip1=nbp_lon+1
    jjp1=nbp_lat
    ip1jm=(nbp_lon+1)*(nbp_lat-1)
    ip1jmp1=(nbp_lon+1)*nbp_lat

! Options for Heat transport
    ! Alpha variable?
      alpha_var=.FALSE.
      CALL getin('slab_alphav',alpha_var)
      print *,'alpha variable',alpha_var
!  centered ou upstream scheme for meridional transport
      slab_upstream=.FALSE.
      CALL getin('slab_upstream',slab_upstream)
      print *,'upstream slab scheme',slab_upstream
! Sverdrup balance at equator ?
      slab_sverdrup=.TRUE.
      CALL getin('slab_sverdrup',slab_sverdrup)
      print *,'Sverdrup balance',slab_sverdrup
! Use tauy for meridional flux at equator ?
      slab_tyeq=.TRUE.
      CALL getin('slab_tyeq',slab_tyeq)
      print *,'Tauy forcing at equator',slab_tyeq
! Use tauy for meridional flux at equator ?
      ekman_zonadv=.TRUE.
      CALL getin('slab_ekman_zonadv',ekman_zonadv)
      print *,'Use Ekman flow in zonal direction',ekman_zonadv
! Use tauy for meridional flux at equator ?
      ekman_zonavg=.FALSE.
      CALL getin('slab_ekman_zonavg',ekman_zonavg)
      print *,'Use zonally-averaged wind stress ?',ekman_zonavg
! Value of alpha
      alpham=2./3.
      CALL getin('slab_alpha',alpham)
      print *,'slab_alpha',alpham
! GM k coefficient (m2/s) for 2-layers
      gmkappa=1000.
      CALL getin('slab_gmkappa',gmkappa)
      print *,'slab_gmkappa',gmkappa
! GM k coefficient (m2/s) for 2-layers
      gm_smax=2e-3
      CALL getin('slab_gm_smax',gm_smax)
      print *,'slab_gm_smax',gm_smax
! -----------------------------------------------------------
! Define ocean / continent mask (no flux into continent cell)
    allocate(zmasqu(ip1jmp1))
    allocate(zmasqv(ip1jm))
    zmasqu=1.
    zmasqv=1.

    ! convert mask to 2D grid
    CALL gr_fi_dyn(1,iip1,jjp1,zmasq,zmasq_2d)
    ! put flux mask to 0 at boundaries of continent cells
    DO i=1,ip1jmp1-1
      IF (zmasq_2d(i).GT.1e-5 .OR. zmasq_2d(i+1).GT.1e-5) THEN
              zmasqu(i)=0.
      ENDIF
    END DO
    DO i=iip1,ip1jmp1,iip1
      zmasqu(i)=zmasqu(i-iim)
    END DO
    DO i=1,ip1jm
      IF (zmasq_2d(i).GT.1e-5 .OR. zmasq_2d(i+iip1).GT.1e-5) THEN
              zmasqv(i)=0.
      END IF
    END DO

! -----------------------------------------------------------
! Coriolis and friction for Ekman transport
    slab_ekman=2
    CALL getin("slab_ekman",slab_ekman)
    IF (slab_ekman.GT.0) THEN
      allocate(unsfv(ip1jm))
      allocate(unsev(ip1jm))
      allocate(unsfu(ip1jmp1))
      allocate(unseu(ip1jmp1))
      allocate(unsbv(ip1jm))

      eps=1e-5 ! Drag
      CALL getin('slab_eps',eps)
      print *,'epsilon=',eps
      ff=fext*unsairez ! Coriolis
      ! coefs to convert tau_x, tau_y to Ekman mass fluxes
      ! on 2D grid v points
      ! Compute correction factor [0 1] near the equator (f<<eps)
      IF (slab_sverdrup) THEN
         ! New formulation, sharper near equator, when eps gives Rossby Radius
         DO i=1,ip1jm
           unsev(i)=exp(-ff(i)*ff(i)/eps**2)
         ENDDO
      ELSE
         DO i=1,ip1jm
           unsev(i)=eps**2/(ff(i)*ff(i)+eps**2)
         ENDDO
      END IF ! slab_sverdrup
      ! 1/beta
      DO i=1,jjp1-1
        unsbv((i-1)*iip1+1:i*iip1)=unsev((i-1)*iip1+1:i*iip1)/beta(i)
      END DO
      ! 1/f
      ff=SIGN(MAX(ABS(ff),eps/100.),ff) ! avoid value 0 at equator...
      DO i=1,ip1jm
        unsfv(i)=(1.-unsev(i))/ff(i)
      END DO
      ! compute values on 2D u grid
      ! 1/eps
      unsev(:)=unsev(:)/eps
      CALL gr_v_scal(1,unsfv,unsfu)
      CALL gr_v_scal(1,unsev,unseu)
    END IF
  
  END SUBROUTINE ini_slab_transp

  SUBROUTINE divgrad_phy(nlevs,temp,delta)
! Computes temperature tendency due to horizontal diffusion :
! T Laplacian, later multiplied by diffusion coef and time-step

    IMPLICIT NONE
 
    INTEGER, INTENT(IN) :: nlevs ! nlevs : slab layers
    REAL, INTENT(IN)   :: temp(klon_glo,nlevs) ! slab temperature
    REAL , INTENT(OUT) :: delta(klon_glo,nlevs) ! temp laplacian (heat flux div.)
    REAL :: delta_2d((nbp_lon+1)*nbp_lat,nlevs)
    REAL ghx((nbp_lon+1)*nbp_lat,nlevs), ghy((nbp_lon+1)*(nbp_lat-1),nlevs)
    INTEGER :: ll,iip1,jjp1
 
    iip1=nbp_lon+1
    jjp1=nbp_lat
    
    ! transpose temp to 2D horiz. grid
    CALL gr_fi_dyn(nlevs,iip1,jjp1,temp,delta_2d)
    ! computes gradient (proportional to heat flx)
    CALL grad(nlevs,delta_2d,ghx,ghy)
    ! put flux to 0 at ocean / continent boundary
    DO ll=1,nlevs
        ghx(:,ll)=ghx(:,ll)*zmasqu
        ghy(:,ll)=ghy(:,ll)*zmasqv
    END DO
    ! flux divergence
    CALL diverg(nlevs,ghx,ghy,delta_2d)
    ! laplacian back to 1D grid
    CALL gr_dyn_fi(nlevs,iip1,jjp1,delta_2d,delta)
 
    RETURN
  END SUBROUTINE divgrad_phy

  SUBROUTINE slab_ekman1(tx_phy,ty_phy,ts_phy,dt_phy)
! 1.5 Layer Ekman transport temperature tendency
! note : tendency dt later multiplied by (delta t)/(rho.H)
! to convert from divergence of heat fluxes to T

      IMPLICIT NONE
      
      ! tx, ty : wind stress (different grids)
      ! fluxm, fluz : mass *or heat* fluxes
      ! dt : temperature tendency
      INTEGER ij

      ! ts surface temp, td deep temp (diagnosed)
      REAL ts_phy(klon_glo)
      REAL ts((nbp_lon+1)*nbp_lat),td((nbp_lon+1)*nbp_lat)
      ! zonal and meridional wind stress 
      REAL tx_phy(klon_glo),ty_phy(klon_glo)
      REAL tyu((nbp_lon+1)*nbp_lat),txu((nbp_lon+1)*nbp_lat)
      REAL txv((nbp_lon+1)*(nbp_lat-1)),tyv((nbp_lon+1)*(nbp_lat-1))
      REAL tcurl((nbp_lon+1)*(nbp_lat-1))
      ! zonal and meridional Ekman mass fluxes at u, v points (2D grid)
      REAL fluxz((nbp_lon+1)*nbp_lat),fluxm((nbp_lon+1)*(nbp_lat-1))
      ! vertical  and absolute mass fluxes (to estimate alpha)
      REAL fluxv((nbp_lon+1)*nbp_lat),fluxt((nbp_lon+1)*(nbp_lat-1))
      ! temperature tendency
      REAL dt((nbp_lon+1)*nbp_lat),dt_phy(klon_glo)
      REAL alpha((nbp_lon+1)*nbp_lat) ! deep temperature coef

      INTEGER iim,iip1,iip2,jjp1,ip1jm,ip1jmi1,ip1jmp1

! Grid definitions
      iim=nbp_lon
      iip1=nbp_lon+1
      iip2=nbp_lon+2
      jjp1=nbp_lat
      ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm
      ip1jmi1=(nbp_lon+1)*(nbp_lat-1)-(nbp_lon+1) ! = ip1jm - iip1
      ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1

! Convert taux,y to 2D  scalar grid
! Note: 2D grid size = iim*jjm. iip1=iim+1
! First and last points in zonal direction are the same 
! we use 1 index ij from 1 to (iim+1)*(jjm+1)
      ! north and south poles
      tx_phy(1)=0.
      tx_phy(klon_glo)=0.
      ty_phy(1)=0.
      ty_phy(klon_glo)=0.
      CALL gr_fi_dyn(1,iip1,jjp1,tx_phy,txu)
      CALL gr_fi_dyn(1,iip1,jjp1,ty_phy,tyu)
! convert to u,v grid (Arakawa C)
! Multiply by f or eps to get mass flux
      ! Meridional mass flux
      CALL gr_scal_v(1,txu,txv) ! wind stress at v points
      IF (slab_sverdrup) THEN ! Sverdrup bal. near equator
        tcurl=(txu(1:ip1jm)-txu(iip2:ip1jmp1))/cv(:)
        fluxm=-tcurl*unsbv-txv*unsfv ! in kg.s-1.m-1 (zonal distance)
      ELSE 
        CALL gr_scal_v(1,tyu,tyv)
        fluxm=tyv*unsev-txv*unsfv ! in kg.s-1.m-1 (zonal distance)
      ENDIF
      ! Zonal mass flux
      CALL gr_scal_u(1,txu,txu) ! wind stress at u points
      CALL gr_scal_u(1,tyu,tyu)
      fluxz=tyu*unsfu+txu*unseu
            
! Correct flux: continent mask and horiz grid size
      ! multiply m-flux by mask and dx: flux in kg.s-1
      fluxm=fluxm*cv*cuvsurcv*zmasqv 
      ! multiply z-flux by mask and dy: flux in kg.s-1
      fluxz=fluxz*cu*cvusurcu*zmasqu

! Compute vertical  and absolute mass flux (for variable alpha)
      IF (alpha_var) THEN
        DO ij=iip2,ip1jm
        fluxv(ij)=fluxz(ij)-fluxz(ij-1)-fluxm(ij)+fluxm(ij-iip1)
        fluxt(ij)=ABS(fluxz(ij))+ABS(fluxz(ij-1)) &
               +ABS(fluxm(ij))+ABS(fluxm(ij-iip1))
        ENDDO
        DO ij=iip1,ip1jmi1,iip1
            fluxt(ij+1)=fluxt(ij+iip1)
            fluxv(ij+1)=fluxv(ij+iip1)
        END DO
        fluxt(1)=SUM(ABS(fluxm(1:iim)))
        fluxt(ip1jmp1)=SUM(ABS(fluxm(ip1jm-iim:ip1jm-1)))
        fluxv(1)=-SUM(fluxm(1:iim))
        fluxv(ip1jmp1)=SUM(fluxm(ip1jm-iim:ip1jm-1))
        fluxt=MAX(fluxt,1.e-10)
      ENDIF

! Compute alpha coefficient. 
! Tdeep = Tsurf * alpha + 271.35 * (1-alpha)
      IF (alpha_var) THEN
          ! increase alpha (and Tdeep) in downwelling regions
          ! and decrease in upwelling regions
          ! to avoid "hot spots" where there is surface convergence
          DO ij=iip2,ip1jm
              alpha(ij)=alpham-fluxv(ij)/fluxt(ij)*(1.-alpham)
          ENDDO
          alpha(1:iip1)=alpham-fluxv(1)/fluxt(1)*(1.-alpham)
          alpha(ip1jm+1:ip1jmp1)=alpham-fluxv(ip1jmp1)/fluxt(ip1jmp1)*(1.-alpham)
      ELSE
          alpha(:)=alpham
          ! Tsurf-Tdeep ~ 10° in the Tropics
      ENDIF

! Estimate deep temperature
      CALL gr_fi_dyn(1,iip1,jjp1,ts_phy,ts)
      DO ij=1,ip1jmp1
         td(ij)=271.35+(ts(ij)-271.35)*alpha(ij)
         td(ij)=MIN(td(ij),ts(ij))
      END DO
       
! Meridional heat flux: multiply mass flux by (ts-td)
! flux in K.kg.s-1
      IF (slab_upstream) THEN
        ! upstream scheme to avoid hot spots
        DO ij=1,ip1jm
        IF (fluxm(ij).GE.0.) THEN
           fluxm(ij)=fluxm(ij)*(ts(ij+iip1)-td(ij))
        ELSE
           fluxm(ij)=fluxm(ij)*(ts(ij)-td(ij+iip1))
        END IF
        END DO
      ELSE
        ! centered scheme better in mid-latitudes
        DO ij=1,ip1jm
          fluxm(ij)=fluxm(ij)*(ts(ij+iip1)+ts(ij)-td(ij)-td(ij+iip1))/2.
        END DO
      ENDIF

! Zonal heat flux
! upstream scheme
      DO ij=iip2,ip1jm
          fluxz(ij)=fluxz(ij)*(ts(ij)+ts(ij+1)-td(ij+1)-td(ij))/2.
      END DO
      DO ij=iip1*2,ip1jmp1,iip1
           fluxz(ij)=fluxz(ij-iim)
      END DO

! temperature tendency = divergence of heat fluxes
! dt in K.s-1.kg.m-2 (T trend times mass/horiz surface)
      DO ij=iip2,ip1jm
         dt(ij)=(fluxz(ij-1)-fluxz(ij)+fluxm(ij)-fluxm(ij-iip1)) &
               /aire(ij) ! aire : grid area 
      END DO
      DO ij=iip1,ip1jmi1,iip1
         dt(ij+1)=dt(ij+iip1) 
      END DO
! special treatment at the Poles
      dt(1)=SUM(fluxm(1:iim))/apoln
      dt(ip1jmp1)=-SUM(fluxm(ip1jm-iim:ip1jm-1))/apols
      dt(2:iip1)=dt(1)
      dt(ip1jm+1:ip1jmp1)=dt(ip1jmp1)

! tendencies back to 1D grid
      CALL gr_dyn_fi(1,iip1,jjp1,dt,dt_phy)

      RETURN
  END SUBROUTINE slab_ekman1

  SUBROUTINE slab_ekman2(tx_phy,ty_phy,ts_phy,dt_phy_ek,dt_phy_gm,slab_gm)
! Temperature tendency for 2-layers slab ocean
! note : tendency dt later multiplied by (delta time)/(rho.H)
! to convert from divergence of heat fluxes to T

! 11/16 : Inclusion of GM-like eddy advection

      IMPLICIT NONE
      
      LOGICAL,INTENT(in) :: slab_gm
      ! Here, temperature and flux variables are on 2 layers
      INTEGER ij

      ! wind stress variables
      REAL tx_phy(klon_glo),ty_phy(klon_glo)
      REAL txv((nbp_lon+1)*(nbp_lat-1)), tyv((nbp_lon+1)*(nbp_lat-1))
      REAL tyu((nbp_lon+1)*nbp_lat),txu((nbp_lon+1)*nbp_lat)
      REAL tcurl((nbp_lon+1)*(nbp_lat-1))
      ! slab temperature on  1D, 2D grid
      REAL ts_phy(klon_glo,2), ts((nbp_lon+1)*nbp_lat,2)
      ! Temperature gradient, v-points
      REAL dty((nbp_lon+1)*(nbp_lat-1)),dtx((nbp_lon+1)*nbp_lat)
      ! Vertical temperature difference, V-points
      REAL dtz((nbp_lon+1)*(nbp_lat-1))
      ! zonal and meridional mass fluxes at u, v points (2D grid)
      REAL fluxz((nbp_lon+1)*nbp_lat), fluxm((nbp_lon+1)*(nbp_lat-1))
      ! vertical mass flux between the 2 layers
      REAL fluxv_ek((nbp_lon+1)*nbp_lat)
      REAL fluxv_gm((nbp_lon+1)*nbp_lat)
      ! zonal and meridional heat fluxes
      REAL fluxtz((nbp_lon+1)*nbp_lat,2)
      REAL fluxtm((nbp_lon+1)*(nbp_lat-1),2)
      ! temperature tendency (in K.s-1.kg.m-2)
      REAL dt_ek((nbp_lon+1)*nbp_lat,2), dt_phy_ek(klon_glo,2)
      REAL dt_gm((nbp_lon+1)*nbp_lat,2), dt_phy_gm(klon_glo,2)
      ! helper vars
      REAL zonavg, fluxv
      REAL, PARAMETER :: sea_den=1025. ! sea water density

      INTEGER iim,iip1,iip2,jjp1,ip1jm,ip1jmi1,ip1jmp1

! Grid definitions
      iim=nbp_lon
      iip1=nbp_lon+1
      iip2=nbp_lon+2
      jjp1=nbp_lat
      ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm
      ip1jmi1=(nbp_lon+1)*(nbp_lat-1)-(nbp_lon+1) ! = ip1jm - iip1
      ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1
! Convert temperature to 2D grid 
      CALL gr_fi_dyn(2,iip1,jjp1,ts_phy,ts)

! ------------------------------------
! Ekman mass fluxes and Temp tendency
! ------------------------------------
! Convert taux,y to 2D  scalar grid
      ! north and south poles tx,ty no meaning
      tx_phy(1)=0.
      tx_phy(klon_glo)=0.
      ty_phy(1)=0.
      ty_phy(klon_glo)=0.
      CALL gr_fi_dyn(1,iip1,jjp1,tx_phy,txu)
      CALL gr_fi_dyn(1,iip1,jjp1,ty_phy,tyu)
      IF (ekman_zonavg) THEN ! use zonal average of wind stress
        DO ij=1,jjp1-2
          zonavg=SUM(txu(ij*iip1+1:ij*iip1+iim))/iim
          txu(ij*iip1+1:(ij+1)*iip1)=zonavg
          zonavg=SUM(tyu(ij*iip1+1:ij*iip1+iim))/iim
          tyu(ij*iip1+1:(ij+1)*iip1)=zonavg
        END DO
      END IF
          
! Divide taux,y by f or eps, and convert to 2D u,v grids
! (Arakawa C grid)
      ! Meridional flux
      CALL gr_scal_v(1,txu,txv) ! wind stress at v points
      fluxm=-txv*unsfv ! in kg.s-1.m-1 (zonal distance)
      IF (slab_sverdrup) THEN ! Sverdrup bal. near equator
        tcurl=(txu(1:ip1jm)-txu(iip2:ip1jmp1))/cv(:) ! dtx/dy
        !poles curl = 0
        tcurl(1:iip1)=0.
        tcurl(ip1jmi1+1:ip1jm)=0.
        fluxm=fluxm-tcurl*unsbv
      ENDIF
      IF (slab_tyeq) THEN ! meridional wind forcing at equator
        CALL gr_scal_v(1,tyu,tyv)
        fluxm=fluxm+tyv*unsev ! in kg.s-1.m-1 (zonal distance)
      ENDIF
!  apply continent mask, multiply by horiz grid dimension
      fluxm=fluxm*cv*cuvsurcv*zmasqv

      ! Zonal flux
      IF (ekman_zonadv) THEN
        CALL gr_scal_u(1,txu,txu) ! wind stress at u points
        CALL gr_scal_u(1,tyu,tyu)
        fluxz=tyu*unsfu+txu*unseu
        !  apply continent mask, multiply by horiz grid dimension
        fluxz=fluxz*cu*cvusurcu*zmasqu
      END IF
            
!  Vertical mass flux from mass budget (divergence of horiz fluxes)
      IF (ekman_zonadv) THEN
         DO ij=iip2,ip1jm
           fluxv_ek(ij)=fluxz(ij)-fluxz(ij-1)-fluxm(ij)+fluxm(ij-iip1)
         ENDDO
      ELSE 
         DO ij=iip2,ip1jm
           fluxv_ek(ij)=-fluxm(ij)+fluxm(ij-iip1)
         ENDDO
      END IF
      DO ij=iip1,ip1jmi1,iip1
         fluxv_ek(ij+1)=fluxv_ek(ij+iip1)
      END DO
!  vertical mass flux at Poles
      fluxv_ek(1)=-SUM(fluxm(1:iim))     
      fluxv_ek(ip1jmp1)=SUM(fluxm(ip1jm-iim:ip1jm-1))

! Meridional heat fluxes 
      DO ij=1,ip1jm
          ! centered scheme
          fluxtm(ij,1)=fluxm(ij)*(ts(ij+iip1,1)+ts(ij,1))/2.
          fluxtm(ij,2)=-fluxm(ij)*(ts(ij+iip1,2)+ts(ij,2))/2.
      END DO

! Zonal heat fluxes
! Schema upstream      
      IF (ekman_zonadv) THEN
        DO ij=iip2,ip1jm
          IF (fluxz(ij).GE.0.) THEN
                 fluxtz(ij,1)=fluxz(ij)*ts(ij,1)
                 fluxtz(ij,2)=-fluxz(ij)*ts(ij+1,2)
          ELSE
                 fluxtz(ij,1)=fluxz(ij)*ts(ij+1,1)
                 fluxtz(ij,2)=-fluxz(ij)*ts(ij,2)
          ENDIF
        ENDDO
        DO ij=iip1*2,ip1jmp1,iip1
               fluxtz(ij,:)=fluxtz(ij-iim,:)
        END DO
      ELSE
        fluxtz(:,:)=0.
      ENDIF        
                   
! Temperature tendency, horizontal advection:
      DO ij=iip2,ip1jm
         dt_ek(ij,:)=fluxtz(ij-1,:)-fluxtz(ij,:) &
                 +fluxtm(ij,:)-fluxtm(ij-iip1,:)
      END DO
! Poles
      dt_ek(1,:)=SUM(fluxtm(1:iim,:),dim=1)
      dt_ek(ip1jmp1,:)=-SUM(fluxtm(ip1jm-iim:ip1jm-1,:),dim=1)

! ------------------------------------
! GM mass fluxes and Temp tendency
! ------------------------------------
     IF (slab_gm) THEN
! Vertical Temperature difference T1-T2 on v-grid points
      CALL gr_scal_v(1,ts(:,1)-ts(:,2),dtz)
      dtz(:)=MAX(dtz(:),0.25) 
! Horizontal Temperature differences 
      CALL grad(1,(ts(:,1)+ts(:,2))/2.,dtx,dty)
! Meridional flux = -k.s (s=dyT/dzT)
! Continent mask, multiply by dz/dy 
      fluxm=dty/dtz*500.*cuvsurcv*zmasqv
! slope limitation, multiply by kappa
      fluxm=-gmkappa*SIGN(MIN(ABS(fluxm),gm_smax*cv*cuvsurcv),dty)
! convert to kg/s
      fluxm(:)=fluxm(:)*sea_den
! Zonal flux = 0. (temporary)
      fluxz(:)=0.
!  Vertical mass flux from mass budget (divergence of horiz fluxes)
      DO ij=iip2,ip1jm
        fluxv_gm(ij)=fluxz(ij)-fluxz(ij-1)-fluxm(ij)+fluxm(ij-iip1)
      ENDDO
      DO ij=iip1,ip1jmi1,iip1
         fluxv_gm(ij+1)=fluxv_gm(ij+iip1)
      END DO
!  vertical mass flux at Poles
      fluxv_gm(1)=-SUM(fluxm(1:iim))     
      fluxv_gm(ip1jmp1)=SUM(fluxm(ip1jm-iim:ip1jm-1))

! Meridional heat fluxes 
      DO ij=1,ip1jm
          ! centered scheme
          fluxtm(ij,1)=fluxm(ij)*(ts(ij+iip1,1)+ts(ij,1))/2.
          fluxtm(ij,2)=-fluxm(ij)*(ts(ij+iip1,2)+ts(ij,2))/2.
      END DO

! Zonal heat fluxes
! Schema upstream      
      DO ij=iip2,ip1jm
      IF (fluxz(ij).GE.0.) THEN
             fluxtz(ij,1)=fluxz(ij)*ts(ij,1)
             fluxtz(ij,2)=-fluxz(ij)*ts(ij+1,2)
      ELSE
             fluxtz(ij,1)=fluxz(ij)*ts(ij+1,1)
             fluxtz(ij,2)=-fluxz(ij)*ts(ij,2)
      ENDIF
      ENDDO
      DO ij=iip1*2,ip1jmp1,iip1
             fluxtz(ij,:)=fluxtz(ij-iim,:)
      END DO
                   
! Temperature tendency :
! divergence of horizontal heat fluxes
      DO ij=iip2,ip1jm
         dt_gm(ij,:)=fluxtz(ij-1,:)-fluxtz(ij,:) &
                 +fluxtm(ij,:)-fluxtm(ij-iip1,:)
      END DO
! Poles
      dt_gm(1,:)=SUM(fluxtm(1:iim,:),dim=1)
      dt_gm(ip1jmp1,:)=-SUM(fluxtm(ip1jm-iim:ip1jm-1,:),dim=1)
     ELSE
       dt_gm(:,:)=0.
       fluxv_gm(:)=0.
     ENDIF ! slab_gm

! ------------------------------------
! Temp tendency from vertical advection
! Divide by cell area
! ------------------------------------
! vertical heat flux = mass flux * T, upstream scheme
      DO ij=iip2,ip1jm
         fluxv=fluxv_ek(ij)+fluxv_gm(ij) ! net flux, needed for upstream scheme
         IF (fluxv.GT.0.) THEN
           dt_ek(ij,1)=dt_ek(ij,1)+fluxv_ek(ij)*ts(ij,2)
           dt_ek(ij,2)=dt_ek(ij,2)-fluxv_ek(ij)*ts(ij,2)
           dt_gm(ij,1)=dt_gm(ij,1)+fluxv_gm(ij)*ts(ij,2)
           dt_gm(ij,2)=dt_gm(ij,2)-fluxv_gm(ij)*ts(ij,2)
         ELSE
           dt_ek(ij,1)=dt_ek(ij,1)+fluxv_ek(ij)*ts(ij,1)
           dt_ek(ij,2)=dt_ek(ij,2)-fluxv_ek(ij)*ts(ij,1)
           dt_gm(ij,1)=dt_gm(ij,1)+fluxv_gm(ij)*ts(ij,1)
           dt_gm(ij,2)=dt_gm(ij,2)-fluxv_gm(ij)*ts(ij,1)
         ENDIF
         ! divide by cell area
         dt_ek(ij,:)=dt_ek(ij,:)/aire(ij)
         dt_gm(ij,:)=dt_gm(ij,:)/aire(ij)
      END DO
      ! North Pole
      fluxv=fluxv_ek(1)+fluxv_gm(1)
        IF (fluxv.GT.0.) THEN
          dt_ek(1,1)=dt_ek(1,1)+fluxv_ek(1)*ts(1,2)
          dt_ek(1,2)=dt_ek(1,2)-fluxv_ek(1)*ts(1,2)
          dt_gm(1,1)=dt_gm(1,1)+fluxv_gm(1)*ts(1,2)
          dt_gm(1,2)=dt_gm(1,2)-fluxv_gm(1)*ts(1,2)
        ELSE
          dt_ek(1,1)=dt_ek(1,1)+fluxv_ek(1)*ts(1,1)
          dt_ek(1,2)=dt_ek(1,2)-fluxv_ek(1)*ts(1,1)
          dt_gm(1,1)=dt_gm(1,1)+fluxv_gm(1)*ts(1,1)
          dt_gm(1,2)=dt_gm(1,2)-fluxv_gm(1)*ts(1,1)
        ENDIF
      dt_ek(1,:)=dt_ek(1,:)/apoln
      dt_gm(1,:)=dt_gm(1,:)/apoln
      ! South pole
      fluxv=fluxv_ek(ip1jmp1)+fluxv_gm(ip1jmp1)
        IF (fluxv.GT.0.) THEN
          dt_ek(ip1jmp1,1)=dt_ek(ip1jmp1,1)+fluxv_ek(ip1jmp1)*ts(ip1jmp1,2)
          dt_ek(ip1jmp1,2)=dt_ek(ip1jmp1,2)-fluxv_ek(ip1jmp1)*ts(ip1jmp1,2)
          dt_gm(ip1jmp1,1)=dt_gm(ip1jmp1,1)+fluxv_gm(ip1jmp1)*ts(ip1jmp1,2)
          dt_gm(ip1jmp1,2)=dt_gm(ip1jmp1,2)-fluxv_gm(ip1jmp1)*ts(ip1jmp1,2)
        ELSE
          dt_ek(ip1jmp1,1)=dt_ek(ip1jmp1,1)+fluxv_ek(ip1jmp1)*ts(ip1jmp1,1)
          dt_ek(ip1jmp1,2)=dt_ek(ip1jmp1,2)-fluxv_ek(ip1jmp1)*ts(ip1jmp1,1)
          dt_gm(ip1jmp1,1)=dt_gm(ip1jmp1,1)+fluxv_gm(ip1jmp1)*ts(ip1jmp1,1)
          dt_gm(ip1jmp1,2)=dt_gm(ip1jmp1,2)-fluxv_gm(ip1jmp1)*ts(ip1jmp1,1)
        ENDIF
      dt_ek(ip1jmp1,:)=dt_ek(ip1jmp1,:)/apols
      dt_gm(ip1jmp1,:)=dt_gm(ip1jmp1,:)/apols
      
      dt_ek(2:iip1,1)=dt_ek(1,1)
      dt_ek(2:iip1,2)=dt_ek(1,2)
      dt_gm(2:iip1,1)=dt_gm(1,1)
      dt_gm(2:iip1,2)=dt_gm(1,2)
      dt_ek(ip1jm+1:ip1jmp1,1)=dt_ek(ip1jmp1,1)
      dt_ek(ip1jm+1:ip1jmp1,2)=dt_ek(ip1jmp1,2)
      dt_gm(ip1jm+1:ip1jmp1,1)=dt_gm(ip1jmp1,1)
      dt_gm(ip1jm+1:ip1jmp1,2)=dt_gm(ip1jmp1,2)

      DO ij=iip1,ip1jmi1,iip1
         dt_gm(ij+1,:)=dt_gm(ij+iip1,:) 
         dt_ek(ij+1,:)=dt_ek(ij+iip1,:) 
      END DO

! T tendency back to 1D grid...
      CALL gr_dyn_fi(2,iip1,jjp1,dt_ek,dt_phy_ek)
      CALL gr_dyn_fi(2,iip1,jjp1,dt_gm,dt_phy_gm)

      RETURN
  END SUBROUTINE slab_ekman2

  SUBROUTINE slab_gmdiff(ts_phy,dt_phy)
! Temperature tendency for 2-layers slab ocean
! Due to Gent-McWilliams type eddy-induced advection
      
      IMPLICIT NONE
      
      ! Here, temperature and flux variables are on 2 layers
      INTEGER ij
      ! Temperature gradient, v-points
      REAL dty((nbp_lon+1)*(nbp_lat-1)),dtx((nbp_lon+1)*nbp_lat)
      ! Vertical temperature difference, V-points
      REAL dtz((nbp_lon+1)*(nbp_lat-1))
      ! slab temperature on  1D, 2D grid
      REAL ts_phy(klon_glo,2),ts((nbp_lon+1)*nbp_lat,2)
      ! zonal and meridional mass fluxes at u, v points (2D grid)
      REAL fluxz((nbp_lon+1)*nbp_lat), fluxm((nbp_lon+1)*(nbp_lat-1))
      ! vertical mass flux between the 2 layers
      REAL fluxv((nbp_lon+1)*nbp_lat)
      ! zonal and meridional heat fluxes
      REAL fluxtz((nbp_lon+1)*nbp_lat,2)
      REAL fluxtm((nbp_lon+1)*(nbp_lat-1),2)
      ! temperature tendency (in K.s-1.kg.m-2)
      REAL dt((nbp_lon+1)*nbp_lat,2), dt_phy(klon_glo,2)

      INTEGER iim,iip1,iip2,jjp1,ip1jm,ip1jmi1,ip1jmp1

! Grid definitions
      iim=nbp_lon
      iip1=nbp_lon+1
      iip2=nbp_lon+2
      jjp1=nbp_lat
      ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm
      ip1jmi1=(nbp_lon+1)*(nbp_lat-1)-(nbp_lon+1) ! = ip1jm - iip1
      ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1

! Convert temperature to 2D grid 
      CALL gr_fi_dyn(2,iip1,jjp1,ts_phy,ts)
! Vertical Temperature difference T1-T2 on v-grid points
      CALL gr_scal_v(1,ts(:,1)-ts(:,2),dtz)
      dtz(:)=MAX(dtz(:),0.25) 
! Horizontal Temperature differences 
      CALL grad(1,(ts(:,1)+ts(:,2))/2.,dtx,dty)
! Meridional flux = -k.s (s=dyT/dzT)
! Continent mask, multiply by dz/dy 
      fluxm=dty/dtz*500.*cuvsurcv*zmasqv
! slope limitation, multiply by kappa
      fluxm=-gmkappa*SIGN(MIN(ABS(fluxm),gm_smax*cv*cuvsurcv),dty)
! Zonal flux = 0. (temporary)
      fluxz(:)=0.
!  Vertical mass flux from mass budget (divergence of horiz fluxes)
      DO ij=iip2,ip1jm
        fluxv(ij)=fluxz(ij)-fluxz(ij-1)-fluxm(ij)+fluxm(ij-iip1)
      ENDDO
      DO ij=iip1,ip1jmi1,iip1
         fluxv(ij+1)=fluxv(ij+iip1)
      END DO
!  vertical mass flux at Poles
      fluxv(1)=-SUM(fluxm(1:iim))     
      fluxv(ip1jmp1)=SUM(fluxm(ip1jm-iim:ip1jm-1))
      fluxv=fluxv

! Meridional heat fluxes 
      DO ij=1,ip1jm
          ! centered scheme
          fluxtm(ij,1)=fluxm(ij)*(ts(ij+iip1,1)+ts(ij,1))/2.
          fluxtm(ij,2)=-fluxm(ij)*(ts(ij+iip1,2)+ts(ij,2))/2.
      END DO

! Zonal heat fluxes
! Schema upstream      
      DO ij=iip2,ip1jm
      IF (fluxz(ij).GE.0.) THEN
             fluxtz(ij,1)=fluxz(ij)*ts(ij,1)
             fluxtz(ij,2)=-fluxz(ij)*ts(ij+1,2)
      ELSE
             fluxtz(ij,1)=fluxz(ij)*ts(ij+1,1)
             fluxtz(ij,2)=-fluxz(ij)*ts(ij,2)
      ENDIF
      ENDDO
      DO ij=iip1*2,ip1jmp1,iip1
             fluxtz(ij,:)=fluxtz(ij-iim,:)
      END DO
                   
! Temperature tendency :
      DO ij=iip2,ip1jm
! divergence of horizontal heat fluxes
         dt(ij,:)=fluxtz(ij-1,:)-fluxtz(ij,:) &
                 +fluxtm(ij,:)-fluxtm(ij-iip1,:)
! + vertical heat flux (mass flux * T, upstream scheme)
         IF (fluxv(ij).GT.0.) THEN
           dt(ij,1)=dt(ij,1)+fluxv(ij)*ts(ij,2)
           dt(ij,2)=dt(ij,2)-fluxv(ij)*ts(ij,2)
         ELSE
           dt(ij,1)=dt(ij,1)+fluxv(ij)*ts(ij,1)
           dt(ij,2)=dt(ij,2)-fluxv(ij)*ts(ij,1)
         ENDIF
         ! divide by cell area
         dt(ij,:)=dt(ij,:)/aire(ij)
      END DO
      DO ij=iip1,ip1jmi1,iip1
         dt(ij+1,:)=dt(ij+iip1,:) 
      END DO
! Poles
      dt(1,:)=SUM(fluxtm(1:iim,:),dim=1)
        IF (fluxv(1).GT.0.) THEN
          dt(1,1)=dt(1,1)+fluxv(1)*ts(1,2)
          dt(1,2)=dt(1,2)-fluxv(1)*ts(1,2)
        ELSE
          dt(1,1)=dt(1,1)+fluxv(1)*ts(1,1)
          dt(1,2)=dt(1,2)-fluxv(1)*ts(1,1)
        ENDIF
      dt(1,:)=dt(1,:)/apoln
      dt(ip1jmp1,:)=-SUM(fluxtm(ip1jm-iim:ip1jm-1,:),dim=1)
       IF (fluxv(ip1jmp1).GT.0.) THEN
         dt(ip1jmp1,1)=dt(ip1jmp1,1)+fluxv(ip1jmp1)*ts(ip1jmp1,2)
         dt(ip1jmp1,2)=dt(ip1jmp1,2)-fluxv(ip1jmp1)*ts(ip1jmp1,2)
       ELSE
         dt(ip1jmp1,1)=dt(ip1jmp1,1)+fluxv(ip1jmp1)*ts(ip1jmp1,1)
         dt(ip1jmp1,2)=dt(ip1jmp1,2)-fluxv(ip1jmp1)*ts(ip1jmp1,1)
       ENDIF
      dt(ip1jmp1,:)=dt(ip1jmp1,:)/apols
      dt(2:iip1,1)=dt(1,1)
      dt(2:iip1,2)=dt(1,2)
      dt(ip1jm+1:ip1jmp1,1)=dt(ip1jmp1,1)
      dt(ip1jm+1:ip1jmp1,2)=dt(ip1jmp1,2)

! T tendency back to 1D grid...
      CALL gr_dyn_fi(2,iip1,jjp1,dt,dt_phy)

      RETURN
  END SUBROUTINE slab_gmdiff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE gr_fi_dyn(nfield,im,jm,pfi,pdyn)
  ! Transfer a variable from 1D "physics" grid to 2D "dynamics" grid
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: im,jm,nfield
  REAL,INTENT(IN) :: pfi(klon_glo,nfield) ! on 1D grid
  REAL,INTENT(OUT) :: pdyn(im,jm,nfield) ! on 2D grid

  INTEGER :: i,j,ifield,ig

  DO ifield=1,nfield
    ! Handle poles
    DO i=1,im
      pdyn(i,1,ifield)=pfi(1,ifield)
      pdyn(i,jm,ifield)=pfi(klon_glo,ifield)
    ENDDO
    ! Other points
    DO j=2,jm-1
      ig=2+(j-2)*(im-1)
      CALL SCOPY(im-1,pfi(ig,ifield),1,pdyn(1,j,ifield),1)
      pdyn(im,j,ifield)=pdyn(1,j,ifield)
    ENDDO
  ENDDO ! of DO ifield=1,nfield

  END SUBROUTINE gr_fi_dyn

  SUBROUTINE gr_dyn_fi(nfield,im,jm,pdyn,pfi)
  ! Transfer a variable from 2D "dynamics" grid to 1D "physics" grid
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: im,jm,nfield
  REAL,INTENT(IN) :: pdyn(im,jm,nfield) ! on 2D grid
  REAL,INTENT(OUT) :: pfi(klon_glo,nfield) ! on 1D grid

  INTEGER j,ifield,ig

  CHARACTER (len = 20)                      :: modname = 'slab_heat_transp'
  CHARACTER (len = 80)                      :: abort_message

  ! Sanity check:
  IF(klon_glo.NE.2+(jm-2)*(im-1)) THEN
    abort_message="gr_dyn_fi error, wrong sizes"
    CALL abort_physic(modname,abort_message,1)
  ENDIF

  ! Handle poles
  CALL SCOPY(nfield,pdyn,im*jm,pfi,klon_glo)
  CALL SCOPY(nfield,pdyn(1,jm,1),im*jm,pfi(klon_glo,1),klon_glo)
  ! Other points
  DO ifield=1,nfield
    DO j=2,jm-1
      ig=2+(j-2)*(im-1)
      CALL SCOPY(im-1,pdyn(1,j,ifield),1,pfi(ig,ifield),1)
    ENDDO
  ENDDO

  END SUBROUTINE gr_dyn_fi

  SUBROUTINE  grad(klevel,pg,pgx,pgy)
  ! compute the covariant components pgx,pgy of the gradient of pg
  ! pgx = d(pg)/dx * delta(x) = delta(pg)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: klevel
  REAL,INTENT(IN) :: pg((nbp_lon+1)*nbp_lat,klevel)
  REAL,INTENT(OUT) :: pgx((nbp_lon+1)*nbp_lat,klevel)
  REAL,INTENT(OUT) :: pgy((nbp_lon+1)*(nbp_lat-1),klevel)

  INTEGER :: l,ij
  INTEGER :: iim,iip1,ip1jm,ip1jmp1

  iim=nbp_lon
  iip1=nbp_lon+1
  ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm
  ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1

  DO l=1,klevel
    DO ij=1,ip1jmp1-1
      pgx(ij,l)=pg(ij+1,l)-pg(ij,l)
    ENDDO
    ! correction for pgx(ip1,j,l) ...
    ! ... pgx(iip1,j,l)=pgx(1,j,l) ...
    DO ij=iip1,ip1jmp1,iip1
      pgx(ij,l)=pgx(ij-iim,l)
    ENDDO
    DO ij=1,ip1jm
      pgy(ij,l)=pg(ij,l)-pg(ij+iip1,l)
    ENDDO
  ENDDO

  END SUBROUTINE grad

  SUBROUTINE diverg(klevel,x,y,div)
  ! computes the divergence of a vector field of components
  ! x,y. x and y being covariant components
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: klevel
  REAL,INTENT(IN) :: x((nbp_lon+1)*nbp_lat,klevel)
  REAL,INTENT(IN) :: y((nbp_lon+1)*(nbp_lat-1),klevel)
  REAL,INTENT(OUT) :: div((nbp_lon+1)*nbp_lat,klevel)

  INTEGER :: l,ij
  INTEGER :: iim,iip1,iip2,ip1jm,ip1jmp1,ip1jmi1

  REAL :: aiy1(nbp_lon+1),aiy2(nbp_lon+1)
  REAL :: sumypn,sumyps
  REAL,EXTERNAL :: SSUM

  iim=nbp_lon
  iip1=nbp_lon+1
  iip2=nbp_lon+2
  ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm
  ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1
  ip1jmi1=(nbp_lon+1)*(nbp_lat-1)-(nbp_lon+1) ! = ip1jm - iip1

  DO l=1,klevel
    DO ij=iip2,ip1jm-1
      div(ij+1,l)= &
        cvusurcu(ij+1)*x(ij+1,l)-cvusurcu(ij)*x(ij,l)+ &
        cuvsurcv(ij-iim)*y(ij-iim,l)-cuvsurcv(ij+1)*y(ij+1,l)
    ENDDO
    ! correction for div(1,j,l) ...
    ! ... div(1,j,l)= div(iip1,j,l) ...
    DO ij=iip2,ip1jm,iip1
      div(ij,l)=div(ij+iim,l)
    ENDDO
    ! at the poles
    DO ij=1,iim
      aiy1(ij)=cuvsurcv(ij)*y(ij,l)
      aiy2(ij)=cuvsurcv(ij+ip1jmi1)*y(ij+ip1jmi1,l)
    ENDDO
    sumypn=SSUM(iim,aiy1,1)/apoln
    sumyps=SSUM(iim,aiy2,1)/apols
    DO ij=1,iip1
      div(ij,l)=-sumypn
      div(ij+ip1jm,l)=sumyps
    ENDDO
    ! End (poles)
  ENDDO ! of DO l=1,klevel

  !!! CALL filtreg( div, jjp1, klevel, 2, 2, .TRUE., 1 )
  DO l=1,klevel
    DO ij=iip2,ip1jm
      div(ij,l)=div(ij,l)*unsaire(ij)
    ENDDO
  ENDDO

  END SUBROUTINE diverg

  SUBROUTINE gr_v_scal(nx,x_v,x_scal)
  ! convert values from v points to scalar points on C-grid
  ! used to  compute unsfu, unseu (u points, but depends only on latitude)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nx ! number of levels or fields
  REAL,INTENT(IN) :: x_v((nbp_lon+1)*(nbp_lat-1),nx)
  REAL,INTENT(OUT) :: x_scal((nbp_lon+1)*nbp_lat,nx)

  INTEGER :: l,ij
  INTEGER :: iip1,iip2,ip1jm,ip1jmp1

  iip1=nbp_lon+1
  iip2=nbp_lon+2
  ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm
  ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1

  DO l=1,nx
    DO ij=iip2,ip1jm
      x_scal(ij,l)= &
                   (airev(ij-iip1)*x_v(ij-iip1,l)+airev(ij)*x_v(ij,l)) &
                  /(airev(ij-iip1)+airev(ij))
    ENDDO
    DO ij=1,iip1
      x_scal(ij,l)=0.
    ENDDO
    DO ij=ip1jm+1,ip1jmp1
      x_scal(ij,l)=0.
    ENDDO
  ENDDO

  END SUBROUTINE gr_v_scal
 
  SUBROUTINE gr_scal_v(nx,x_scal,x_v)
  ! convert values from scalar points to v points on C-grid
  ! used to compute wind stress at V points
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nx ! number of levels or fields
  REAL,INTENT(OUT) :: x_v((nbp_lon+1)*(nbp_lat-1),nx)
  REAL,INTENT(IN) :: x_scal((nbp_lon+1)*nbp_lat,nx)

  INTEGER :: l,ij
  INTEGER :: iip1,ip1jm

  iip1=nbp_lon+1
  ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm

      DO l=1,nx
        DO ij=1,ip1jm
          x_v(ij,l)= &
            (cu(ij)*cvusurcu(ij)*x_scal(ij,l)+ &
            cu(ij+iip1)*cvusurcu(ij+iip1)*x_scal(ij+iip1,l)) &
            /(cu(ij)*cvusurcu(ij)+cu(ij+iip1)*cvusurcu(ij+iip1))
        ENDDO
      ENDDO

  END SUBROUTINE gr_scal_v

  SUBROUTINE gr_scal_u(nx,x_scal,x_u)
  ! convert values from scalar points to U points on C-grid
  ! used to compute wind stress at U points
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nx
  REAL,INTENT(OUT) :: x_u((nbp_lon+1)*nbp_lat,nx)
  REAL,INTENT(IN) :: x_scal((nbp_lon+1)*nbp_lat,nx)

  INTEGER :: l,ij
  INTEGER :: iip1,jjp1,ip1jmp1

  iip1=nbp_lon+1
  jjp1=nbp_lat
  ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1

  DO l=1,nx
     DO ij=1,ip1jmp1-1
        x_u(ij,l)= &
         (aire(ij)*x_scal(ij,l)+aire(ij+1)*x_scal(ij+1,l)) &
         /(aire(ij)+aire(ij+1))
     ENDDO
  ENDDO

  CALL SCOPY(nx*jjp1,x_u(1,1),iip1,x_u(iip1,1),iip1)

  END SUBROUTINE gr_scal_u

END MODULE slab_heat_transp_mod
