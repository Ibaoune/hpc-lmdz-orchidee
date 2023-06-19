MODULE dustemission_mod

  IMPLICIT NONE
  !Parameters   
!  INTEGER, PARAMETER     :: nbins=12  ! number of aerosol bins: original 
!  INTEGER, PARAMETER     :: nbins=800  ! number of aerosol bins: for spla
!  INTEGER, PARAMETER     :: nbins=8000  ! number of aerosol bins: for spla
  
  INTEGER, PARAMETER     :: flag_feff=1 ! 0: deactivate feff (drag partition scheme)
  INTEGER, PARAMETER     :: nbins=800  ! number of aerosol bins: for spla
  INTEGER, PARAMETER     :: nmode=3   ! number of soil-dust modes 
  INTEGER, PARAMETER     :: ntyp=5   ! number of soil types 
  INTEGER, PARAMETER     :: nwb=12   ! number of points for the 10m wind
! speed weibull distribution (>=2)
  real   ,parameter     :: z10m=1000. !10m in cm
  REAL,PARAMETER         :: kref=3. !weibull parameter
  INTEGER, PARAMETER    :: nats=14 !number of mineral types (14 here for sand,
                                   ! silt, clay etc.)
  integer, parameter :: nclass=200000


  real   , parameter :: dmin=0.0001
  real   , parameter :: dmax=0.2
  integer, parameter :: nspe=nmode*3+1
  real   ,parameter     :: vkarm=0.41
!JE20150202 : updating scheme to chimere13b <<<
! original values
!  integer, parameter :: div1=3.
!  integer, parameter :: div2=3.
!  integer, parameter :: div3=3.
!  real   , parameter :: e1=3.61/div1
!  real   , parameter :: e2=3.52/div2
!  real   , parameter :: e3=3.46/div3
!  real   , parameter :: rop=2.65 ! particle density g/m3
!  real   , parameter :: roa=0.001227  ! air density g/m3
!  real   , parameter :: pi=3.14159  !!
!  real   , parameter :: gravity=981. !! cm!!
!  real   , parameter :: cd=1.*roa/gravity
! new values
!  logical, parameter :: ok_splatuning=.true.
! Div=3 from S. Alfaro (Sow et al ACPD 2011)
!JE 20150206
!  integer, parameter :: div1=3.
!  integer, parameter :: div2=3.
!  integer, parameter :: div3=3.
  integer, parameter :: div1=6.
  integer, parameter :: div2=6.
  integer, parameter :: div3=6.
  real   , parameter :: e1=3.61/div1
  real   , parameter :: e2=3.52/div2
  real   , parameter :: e3=3.46/div3
  real   , parameter :: rop=2.65 ! particle density g/m3
  real   , parameter :: roa=0.001227  ! air density g/m3
  real   , parameter :: pi=3.14159  !!
  real   , parameter :: gravity=981. !! cm!!
! C=2.61 from Marticorena and Bergametti 1995 instead of Gillete and Chen 2001
! (recommended C=1.1  in supply-limited dust source area.. )
   real   , parameter :: cd=2.61*roa/gravity
!  real   , parameter :: cd=1.0*roa/gravity
!JE20150202>>>>
  real,parameter     :: beta=16300.
  real, parameter, dimension(3) :: diam=(/1.5,6.7,14.2/)
  INTEGER, PARAMETER     :: ndistb=3
  real, parameter, dimension(3) :: sig=(/1.7,1.6,1.5/)

!   INTEGER, PARAMETER     :: nbinsHR=3000 !original
   INTEGER, PARAMETER     :: nbinsHR=30000
!min and max dust size in um
      REAL, PARAMETER     :: sizedustmin=0.0599 ! for spla
      REAL, PARAMETER     :: sizedustmax=63.
!      REAL, PARAMETER     :: sizedustmin=0.09  ! original
!      REAL, PARAMETER     :: sizedustmax=63.


  ! Calc variables
  REAL,DIMENSION(:,:), ALLOCATABLE,SAVE  :: massfrac
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: binsHR
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: binsHRcm
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: itv
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: sizedust !size dust bin (in um)
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: szdcm !the same but (in cm)

  !soil inputs from file donnees_lisa.nc . Should be in the same grid as the
  !model (regridded with nearest neighborhood from surfnew.nc)
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: sol
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: P
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: zos
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: z01
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: z02
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: D
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: A
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: solspe
  INTEGER,DIMENSION(:,:),ALLOCATABLE,SAVE :: masklisa
!!!  INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: maskdust
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: feff
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: feffdbg

  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: sizeclass 
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: sizeclass2 
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: uth
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: uth2
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: srel
  REAL,DIMENSION(:,:),  ALLOCATABLE,SAVE :: srel2

  INTEGER :: nat ! SOL data inside the loop (use as soil type index?)
  REAL :: ustarsalt 
  REAL :: var3a,var3b
  INTEGER  :: ns,nd,nsi,npi,ni ! counters
  INTEGER :: ncl

! outputs
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: m1dflux !fluxes for each soil mode
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: m2dflux
  REAL,DIMENSION(:),  ALLOCATABLE,SAVE :: m3dflux



!$OMP THREADPRIVATE(m1dflux)
!$OMP THREADPRIVATE(m2dflux)
!$OMP THREADPRIVATE(m3dflux)
!$OMP THREADPRIVATE(massfrac) 
!$OMP THREADPRIVATE(binsHR) 
!$OMP THREADPRIVATE(binsHRcm) 
!$OMP THREADPRIVATE(itv) 
!$OMP THREADPRIVATE(sizedust) 
!$OMP THREADPRIVATE(szdcm) 
!$OMP THREADPRIVATE(sol) 
!$OMP THREADPRIVATE(P) 
!$OMP THREADPRIVATE(zos) 
!$OMP THREADPRIVATE(z01) 
!$OMP THREADPRIVATE(z02) 
!$OMP THREADPRIVATE(D) 
!$OMP THREADPRIVATE(A) 
!$OMP THREADPRIVATE(solspe) 
!$OMP THREADPRIVATE(masklisa) 
!!!!$OMP THREADPRIVATE(maskdust) 
!$OMP THREADPRIVATE(feff) 
!$OMP THREADPRIVATE(feffdbg) 
!$OMP THREADPRIVATE(sizeclass) 
!$OMP THREADPRIVATE(sizeclass2) 
!$OMP THREADPRIVATE(uth) 
!$OMP THREADPRIVATE(uth2) 
!$OMP THREADPRIVATE(srel) 
!$OMP THREADPRIVATE(srel2) 



  CONTAINS

!--------------------------------------------------------------------------------------
!======================================================================================
!**************************************************************************************
!======================================================================================
!--------------------------------------------------------------------------------------

  SUBROUTINE dustemission( debutphy, xlat, xlon, &    !Input
                          pctsrf,zu10m,zv10m,wstar, & !Input
                          ale_bl,ale_wake, &          !Input
                          param_wstarBL, param_wstarWAKE, &  !Input
                          emdustacc,emdustcoa,emdustsco,maskdust)    !Output
  USE dimphy
  USE infotrac
  USE write_field_phy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  USE indice_sol_mod

  IMPLICIT NONE
  ! input :
  ! output: flux_sparam_ddfine,flux_sparam_ddcoa,
  ! first: 
  ! Model grid parameters
  REAL,DIMENSION(klon),     INTENT(IN)     :: xlat
  REAL,DIMENSION(klon),     INTENT(IN)     :: xlon
  REAL,DIMENSION(klon,nbsrf), INTENT(IN)     :: pctsrf
  REAL,DIMENSION(klon),INTENT(IN)          :: zu10m   ! 10m zonal wind
  REAL,DIMENSION(klon),INTENT(IN)          :: zv10m   ! meridional 10m wind
  REAL,DIMENSION(klon),INTENT(IN)          :: wstar
  REAL,DIMENSION(klon),INTENT(IN)          :: ale_bl
  REAL,DIMENSION(klon),INTENT(IN)          :: ale_wake
  REAL,DIMENSION(klon), INTENT(IN) :: param_wstarWAKE
  REAL,DIMENSION(klon), INTENT(IN) :: param_wstarBL
 
 
  LOGICAL  :: debutphy ! First physiqs run or not
  ! Intermediate variable: 12 bins emissions
  REAL,DIMENSION(:,:), ALLOCATABLE,SAVE  :: emisbinloc ! vertical emission fluxes

  !OUT variables
  REAL,DIMENSION(klon)  :: emdustacc,emdustcoa,emdustsco ! emission in spla dust bins:
! old radio : acc=0.03-0.5 micrometers, ccoa:0.5-10 micrometers
! new acc=0.03-0.5 micrometers, coa:0.5-3 micrometers ,sco:3-15 um
  INTEGER,DIMENSION(klon) :: maskdust ! where the emissions were calculated
!  INTEGER,DIMENSION(klon_glo) :: maskdust_glo ! auxiliar 
!  REAL,DIMENSION(klon_glo) :: raux_klon_glo ! auxiliar 

!$OMP THREADPRIVATE(emisbinloc)
!!!!!!$OMP THREADPRIVATE(maskdust)
  IF (debutphy) THEN
     ALLOCATE( emisbinloc(klon,nbins) )
  ENDIF

  IF( debutphy ) THEN  
     CALL initdust(xlat,xlon,pctsrf)
  ENDIF

!JE20141124  CALL  calcdustemission(debutphy,zu10m,zv10m,wstar,ale_bl,ale_wake,emisbinloc)
  CALL  calcdustemission(debutphy,zu10m,zv10m,wstar,ale_bl,ale_wake,param_wstarBL,param_wstarWAKE, & !I
                         emisbinloc)   !O

  CALL makemask(maskdust)

  IF( debutphy ) THEN  
!       call gather(maskdust,maskdust_glo)
!     !$OMP MASTER
!     IF (is_mpi_root .AND. is_omp_root) THEN
!       CALL writefield_phy("maskdust",float(maskdust_glo),1)
       CALL writefield_phy("maskdust",float(maskdust),1)
!     ENDIF
!     !$OMP END MASTER
!     !$OMP BARRIER
  ENDIF

  !CALL adaptdustemission(debutphy,emisbinloc,emdustacc,emdustcoa,emdustsco)
  CALL adaptdustemission(debutphy,emisbinloc,emdustacc,emdustcoa,emdustsco,maskdust,pctsrf)
  ! output in kg/m2/s


  END SUBROUTINE dustemission

!--------------------------------------------------------------------------------------
!======================================================================================
!**************************************************************************************
!======================================================================================
!--------------------------------------------------------------------------------------

 SUBROUTINE makemask(maskdustloc) 
  USE dimphy
  USE infotrac
  IMPLICIT NONE
  !Input
  INTEGER,DIMENSION(klon) :: maskdustloc
  INTEGER :: i,j,k
  integer :: iaux
 

  do k=1,klon
      maskdustloc(k)=0
      do i=1,ntyp
         if (masklisa(k,i)>0) then
             maskdustloc(k)=1
         endif
      enddo
  enddo

  END SUBROUTINE makemask

!--------------------------------------------------------------------------------------
!======================================================================================
!**************************************************************************************
!======================================================================================
!--------------------------------------------------------------------------------------

 SUBROUTINE adaptdustemission(debutphy,emisbinlocal, &
                emdustacc,emdustcoa,emdustsco,maskdust,pctsrf)
!                               emdustacc,emdustcoa,emdustsco)

  USE dimphy
  USE infotrac
  USE write_field_phy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  USE indice_sol_mod

  IMPLICIT NONE
  !Input
  REAL,DIMENSION(klon) :: emdustacc,emdustcoa,emdustsco
  REAL, DIMENSION(klon,nbins) :: emisbinlocal
  ! Local
  INTEGER  :: i,j,k
!!!  INTEGER,DIMENSION(2) :: iminacclow,iminacchigh,imincoalow,imincoahigh ! in
!case of small nbins.... not ready
  INTEGER,SAVE ::iminacclow,iminacchigh,imincoalow
  INTEGER,SAVE ::imincoahigh,iminscohigh,iminscolow
  INTEGER,DIMENSION(klon) :: maskdust ! where the emissions were calculated
  REAL,DIMENSION(klon,nbsrf),     INTENT(IN)     :: pctsrf
!  real,parameter :: sizeacclow=0.03
!  real,parameter :: sizeacchigh=0.5
!  real,parameter :: sizecoalow=0.5
!  real,parameter :: sizecoahigh=10.  ! in micrometers
  real,parameter :: sizeacclow=0.06
  real,parameter :: sizeacchigh=1.0
  real,parameter :: sizecoalow=1.0
  real,parameter :: sizecoahigh=6.  !20 ! diameter in micrometers
  real,parameter :: sizescolow=6.
  real,parameter :: sizescohigh=30.  ! in micrometers
!--------------------------------
!  real,parameter :: tuningfactorfine=0.9  ! factor for fine bins!!! important!!
  real,parameter :: tuningfactorfine=0.8  ! factor for fine bins!!! important!!
!  real,parameter :: tuningfactorfine=4.5  ! factor for fine bins!!! important!!
!  real,parameter :: tuningfactorcoa=3.6 ! factor for coarse bins!!! important!!
  real,parameter :: tuningfactorcoa=3.25 ! factor for coarse bins!!! important!!
!  real,parameter :: tuningfactorcoa=4.5  ! factor for coarse bins!!! important!!
!  real,parameter :: tuningfactorsco=3.6  ! factor for supercoarse bins!!! important!!
  real,parameter :: tuningfactorsco=3.25  ! factor for supercoarse bins!!! important!!
!  real,parameter :: tuningfactorsco=4.5  ! factor for supercoarse bins!!! important!!
  real,parameter :: basesumemission= 0.0  !1.e-6  ! emissions to SUM to each land pixel FOR ASSIMILATION ONLY important!!  in mg/m2/s, per bin
 !basesumemission = 1.e-6 increase the AOD in about 12%  (0.03 of AOD) , 
 !while 1e-8 increase in about 0.12%  (0.003 of AOD)

  real,dimension(klon) :: basesumacc,basesumcoa,basesumsco
!--------------------------------
!JE20140915  real,parameter :: sizeacclow=0.06
!JE20140915  real,parameter :: sizeacchigh=1.0
!JE20140915  real,parameter :: sizecoalow=1.0
!JE20140915  real,parameter :: sizecoahigh=10.  !20 ! diameter in micrometers
!JE20140915  real,parameter :: sizescolow=10.
!JE20140915  real,parameter :: sizescohigh=30.  ! in micrometers



  logical ::  debutphy
  real :: diff, auxr1,auxr2,auxr3,auxr4
  real,dimension(klon,nbins) :: itvmean
  real,dimension(klon,nbins+1) :: itv2
!  real,dimension(klon_glo,nbins) :: itvmean_glo
!  real,dimension(:,:) , allocatable  :: itvmean_glo
!  real,dimension(:,:), allocatable :: itv2_glo
  
  integer, save :: counter,counter1 !dbg
  REAL, DIMENSION(:,:),ALLOCATABLE,SAVE :: emisbinlocalmean,emisbinlocalmean2 !dbg
  REAL, DIMENSION(:,:),ALLOCATABLE :: emisbinlocalmean2_glo 
  logical :: writeaerosoldistrib
!$OMP THREADPRIVATE(iminacclow,iminacchigh,imincoalow,imincoahigh)

writeaerosoldistrib=.false.
if (debutphy) then

  if (sizedustmin>sizeacclow .or. sizedustmax<sizescohigh) then 
   call abort_gcm('adaptdustemission', 'Dust range problem',1)
  endif
  print *,'FINE DUST BIN: tuning EMISSION factor= ',tuningfactorfine
  print *,'COA DUST BIN: tuning EMISSION factor= ',tuningfactorcoa
  print *,'SCO DUST BIN: tuning EMISSION factor= ',tuningfactorsco
  print *,'ALL DUST BIN: SUM to the emissions (mg/m2/s) = ',basesumemission
  auxr1=9999.
  auxr2=9999.
  auxr3=9999.
  auxr4=9999.
  do i=1,nbins+1
   if (abs(sizeacclow-itv(i))<auxr1) then
          auxr1=abs( sizeacclow-itv(i))
          iminacclow=i
   endif
   if (abs(sizeacchigh-itv(i))<auxr2) then
          auxr2=abs( sizeacchigh-itv(i))
          iminacchigh=i
          imincoalow=i
   endif
   if (abs(sizecoahigh-itv(i))<auxr3) then
          auxr3=abs( sizecoahigh-itv(i))
          imincoahigh=i
          iminscolow=i
   endif
   if (abs(sizescohigh-itv(i))<auxr4) then
          auxr4=abs( sizescohigh-itv(i))
          iminscohigh=i
   endif
  enddo
if (writeaerosoldistrib) then
!JEdbg<<
  do j=1,klon
  do i=1,nbins
    itvmean(j,i)=(itv(i)+itv(i+1))/2.
    itv2(j,i)=itv(i)
    !print*, itv(i),itvmean(i),itv(i+1)
    !print*, sizedust(i)
  enddo
  itv2(j,nbins+1)=itv(nbins+1)
  enddo
!  allocate(itv2_glo(klon_glo,nbins+1))
!  allocate(itvmean_glo(klon_glo,nbins))
!  ALLOCATE(emisbinlocalmean2_glo(klon_glo,nbins))
!  call gather(itv2,itv2_glo)
!  call gather(itvmean,itvmean_glo)
!!$OMP MASTER
!  IF (is_mpi_root .AND. is_omp_root) THEN
!  CALL writefield_phy("itv2",itv2_glo,nbins+1)
!  CALL writefield_phy("itvmean",itvmean_glo,nbins)
  CALL writefield_phy("itv2",itv2,nbins+1)
  CALL writefield_phy("itvmean",itvmean,nbins)
!  ENDIF
!!$OMP END MASTER
!!$OMP BARRIER
  ALLOCATE(emisbinlocalmean(klon,nbins))
  ALLOCATE(emisbinlocalmean2(klon,nbins))
  do i=1,nbins
    do j=1,klon
       emisbinlocalmean(j,i)=0.0
       emisbinlocalmean2(j,i)=0.0
    enddo
  enddo
  counter=1
  counter1=0
!JEdbg>>
endif
endif


!print *,'JE'
!print *,iminacclow,iminacchigh,imincoalow,imincoahigh

! estimate and integrate bins into only accumulation and coarse
do k=1,klon
  basesumacc(k)=basesumemission*(pctsrf(k,is_ter))*1.e-6 ! from mg/m2/s
  basesumcoa(k)=basesumemission*(pctsrf(k,is_ter))*1.e-6
  basesumsco(k)=basesumemission*(pctsrf(k,is_ter))*1.e-6
enddo


do k=1,klon
auxr1=0.0
auxr2=0.0
auxr3=0.0
  do i=iminacclow,iminacchigh-1
   auxr1=auxr1+emisbinlocal(k,i)
  enddo
  emdustacc(k)=(auxr1 + basesumacc(k))*tuningfactorfine
  do i=imincoalow,imincoahigh-1
    auxr2=auxr2+emisbinlocal(k,i)
  enddo
  emdustcoa(k)=(auxr2 + basesumcoa(k))*tuningfactorcoa
  do i=iminscolow,iminscohigh-1
    auxr3=auxr3+emisbinlocal(k,i)
  enddo
  emdustsco(k)=(auxr3 + basesumsco(k))*tuningfactorsco
enddo


!do k=1,klon
!auxr1=0.0
!auxr2=0.0
!auxr3=0.0
!  do i=iminacclow,iminacchigh-1
!   auxr1=auxr1+emisbinlocal(k,i)
!  enddo
!  emdustacc(k)= auxr1*tuningfactor
!  do i=imincoalow,imincoahigh-1
!    auxr2=auxr2+emisbinlocal(k,i)
!  enddo
!  emdustcoa(k)=auxr2*tuningfactor
!  do i=iminscolow,iminscohigh-1
!    auxr3=auxr3+emisbinlocal(k,i)
!  enddo
!  emdustsco(k)=auxr3*tuningfactor
!enddo
!




!JEdbg<<
if (writeaerosoldistrib) then
  do i=1,nbins
    do j=1,klon
        emisbinlocalmean(j,i)=emisbinlocalmean(j,i)+emisbinlocal(j,i)
    enddo 
  enddo
counter1=counter1+1
! 4 timesteps of physics=4*15min=1hour..
! 1440 = 15 days
! 480 = 5 days
if (MOD(counter,1440).eq. 0) THEN
   !if (MOD(counter,480).eq. 0) THEN
   do k = 1,klon
   do i=1,nbins
     emisbinlocalmean2(k,i)=emisbinlocalmean(k,i)/float(counter1)
   enddo
   enddo
   counter1=0
!   call gather(emisbinlocalmean2,emisbinlocalmean2_glo)
!!$OMP MASTER
!   IF (is_mpi_root .AND. is_omp_root) THEN
!   CALL writefield_phy("emisbinlocalmean2",emisbinlocalmean2_glo,nbins)
   CALL writefield_phy("emisbinlocalmean2",emisbinlocalmean2,nbins)
!   ENDIF
!!$OMP END MASTER
!!$OMP BARRIER
   do i=1,nbins
     do j=1,klon
        emisbinlocalmean(j,i)=0.0
     enddo
   enddo
endif
counter=counter+1
endif
!JEdbg>>

!CALL abort_gcm('calcdustemission', 'OK1',1)
 
END SUBROUTINE adaptdustemission


!--------------------------------------------------------------------------------------
!======================================================================================
!**************************************************************************************
!======================================================================================
!--------------------------------------------------------------------------------------


!--------------------------------------------------------
  SUBROUTINE initdust(xlat,xlon,pctsrf)
  USE dimphy
  USE infotrac
  USE write_field_phy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  USE indice_sol_mod

  IMPLICIT NONE
  !Inputs
  REAL,DIMENSION(klon),     INTENT(IN)     :: xlat
  REAL,DIMENSION(klon),     INTENT(IN)     :: xlon
  ! JE20150605<< better to read
  ! REAL,DIMENSION(klon),     INTENT(IN)     :: pctsrf
  REAL,DIMENSION(klon,nbsrf),     INTENT(IN)     :: pctsrf
  ! JE20150605>>

  !Local
  REAL,DIMENSION(klon,ntyp)                :: solini
  REAL,DIMENSION(klon,ntyp)                :: Pini
  REAL,DIMENSION(klon,ntyp)                :: zosini
  REAL,DIMENSION(klon,ntyp)                :: z01ini
  REAL,DIMENSION(klon,ntyp)                :: z02ini
  REAL,DIMENSION(klon,ntyp)                :: Dini
  REAL,DIMENSION(klon,ntyp)                :: Aini

  REAL,DIMENSION(nclass)                   :: newstep

 !TEMPORAL/ MAKE MASK!!!
  REAL, PARAMETER        :: longmin=-20.
  REAL, PARAMETER        :: longmax=77.          !35.  
  REAL, PARAMETER        :: latmin=10.
  REAL, PARAMETER        :: latmax=40.
 !TEMPORAL/ MAKE MASK!!!
  !Local
  REAL, PARAMETER        :: eps=1.e-5
  REAL, PARAMETER        :: aeff=0.35
  REAL, PARAMETER        :: xeff=10.
  REAL, PARAMETER        :: trctunt=0.310723
 
  INTEGER     :: i,j,k,nts
  REAL :: dp,dstep
! Parametres pour le calcul de la repartition de l energie feff(klon,ntyp)
  REAL :: aa,bb,cc,dd,ee,ff
  REAL,DIMENSION(klon,ntyp) :: tmp1
  REAL :: p3t,var1,var2,et

! Parametres pour le calcul de la surface couverte par les particule de diametre
! dp srel(nats,nclass)
  REAL :: stotale,su,su_loc,xk,xl,xm,xn
  REAL,DIMENSION(nclass)                   :: subsoildist

! isolog and massfrac calcs
  INTEGER :: nbinsout,nb,miniso,nbins1,nbins2,istart,no
  REAL :: b1,b2,stepbin,minisograd,exden,d2min,d2max,numin,numax
  REAL                   :: delta1,delta2,ldmd
!  REAL, PARAMETER     :: sizedustmin=0.09
  REAL,DIMENSION(nbinsHR):: binsISOGRAD
  REAL,DIMENSION(nbinsHR):: vdISOGRAD
  REAL,DIMENSION(nbinsHR):: logvdISOGRAD
  REAL                   :: curvd
  REAL                   :: avdhr
  REAL                   :: bvdhr
  REAL,DIMENSION(nbins)  :: dmn1
  REAL,DIMENSION(nbins)  :: dmn2
  REAL,DIMENSION(nbins)  :: dmn3
  REAL,DIMENSION(nbinsHR) :: vdHR
  REAL,DIMENSION(nbinsHR)  :: vdout



! Parametres pour le calcul de la vitesse de friction seuil uth(nclass)
! calcul de la vitesse seuil selon la parametrisation de Shao and Lu
! 2000(icuth=1).
!  INTEGER                :: ich1
  REAL                   :: an
  REAL                   :: gam
  REAL                   :: sigshao
  REAL                   :: x1
  REAL                   :: x2
! Cas de Iversen and White 1982 (icuth=0) si necessaire.
  REAL, PARAMETER        :: adust=1331.647
  REAL, PARAMETER        :: bdust=0.38194
  REAL, PARAMETER        :: xdust=1.561228
  REAL, PARAMETER        :: ddust=0.006/(rop*gravity)

  CHARACTER(LEN=10) :: varname
  CHARACTER(LEN=80) :: fnsolspe
  CHARACTER(LEN=200) :: line
  

! nats: number of mineral types (14 here for sand, silt, clay etc.)
    ALLOCATE( binsHR(nbinsHR) )
    ALLOCATE( binsHRcm(nbinsHR) )
    ALLOCATE( itv(nbins+1) )
    ALLOCATE( sizedust(nbins) )
    ALLOCATE( szdcm(nbins) )
    ALLOCATE( massfrac(ndistb,nbins) )
    ALLOCATE( sol(klon,ntyp) )
    ALLOCATE( P(klon,ntyp) )
    ALLOCATE( zos(klon,ntyp) )
    ALLOCATE( z01(klon,ntyp) )
    ALLOCATE( z02(klon,ntyp) )
    ALLOCATE( D(klon,ntyp) )
    ALLOCATE( A(klon,ntyp) )
    ALLOCATE( masklisa(klon,ntyp) )
    ALLOCATE( solspe(nats,nspe) )
    ALLOCATE( feff(klon,ntyp) ) 
    ALLOCATE( feffdbg(klon,ntyp) ) 
    ALLOCATE( sizeclass(nclass) )
    ALLOCATE( sizeclass2(nclass) )
    ALLOCATE( uth(nclass) )
    ALLOCATE( uth2(nclass) )
    ALLOCATE( srel(nats,nclass) )
    ALLOCATE( srel2(nats,nclass) )
    ALLOCATE( m1dflux(klon) )
    ALLOCATE( m2dflux(klon) )
    ALLOCATE( m3dflux(klon) )



  ! read input data from donnees_lisa.nc
  varname='SOL'
  CALL read_surface(varname,solini)
  varname='P'
  CALL read_surface(varname,Pini)
  varname='ZOS_'
  CALL read_surface(varname,zosini)
  varname='ZO1_'
  CALL read_surface(varname,z01ini)
  varname='ZO2_'
  CALL read_surface(varname,z02ini)
  varname='D'
  CALL read_surface(varname,Dini)
  varname='A'
  CALL read_surface(varname,Aini)
print *,'beforewritephy',mpi_rank,omp_rank
  CALL writefield_phy("SOLinit",solini,5)
  CALL writefield_phy("Pinit",Pini,5)
  CALL writefield_phy("ZOSinit",zosini,5)
  CALL writefield_phy("ZO1init",z01ini,5)
  CALL writefield_phy("ZO2init",z02ini,5)
  CALL writefield_phy("Dinit",Dini,5)
  CALL writefield_phy("Ainit",Aini,5)
print *,'afterwritephy',mpi_rank,omp_rank

  DO i=1,klon
    DO nts=1,ntyp
        masklisa(i,nts)=0
        IF(Pini(i,nts)>=0.) THEN
               masklisa(i,nts)=1
        ENDIF 
    ENDDO
  ENDDO
!print *,'JEOK1',mpi_rank,omp_rank
  DO i=1,klon
    !print *,Pini(i,1),Pini(i,2),Pini(i,3),Pini(i,4),Pini(i,5)
    DO nts=1,ntyp
      !IF(xlon(i).ge.longmin.and.xlon(i).le.longmax.and. &
      !&      xlat(i).ge.latmin.and.xlat(i).le.latmax    &
      !&      .and.pctsrf(i)>0.5.and.Pini(i,nts)>0.)THEN
  ! JE20150605<< easier to read
  !    IF(pctsrf(i)>0.5.and.Pini(i,nts)>0.)THEN
      IF(pctsrf(i,is_ter)>0.5.and.Pini(i,nts)>0.)THEN
  ! JE20150605>>
           sol(i,nts) = solini(i,nts)
             P(i,nts) = Pini(i,nts)
           zos(i,nts) = zosini(i,nts)
           z01(i,nts) = z01ini(i,nts)
           z02(i,nts) = z02ini(i,nts)
             D(i,nts) = Dini(i,nts)
             A(i,nts) = Aini(i,nts)
!             masklisa(i,nts)=1
      ELSE
              sol(i,nts) =0.0
               P(i,nts)  =0.0
              zos(i,nts) =0.0
              z01(i,nts) =0.0
              z02(i,nts) =0.0
               D(i,nts)  =0.0
               A(i,nts)  =0.0
!            masklisa(i,nts)=0
      ENDIF
    ENDDO
  ENDDO

! print *,'JEOK2',mpi_rank,omp_rank
if ( 1==1 ) then

! print *,'JEOK4',mpi_rank,omp_rank
   CALL writefield_phy("SOL",sol,5)
   CALL writefield_phy("P"  ,P  ,5)
   CALL writefield_phy("ZOS",zos,5)
   CALL writefield_phy("ZO1",z01,5)
   CALL writefield_phy("ZO2",z02,5)
   CALL writefield_phy("D"  ,D  ,5)
   CALL writefield_phy("A"  ,A  ,5)
   CALL writefield_phy("masklisa",float(masklisa),5)
   CALL writefield_phy("pctsrf",pctsrf,1)
   CALL writefield_phy("xlon",xlon,1)
   CALL writefield_phy("xlat",xlat,1)
!print *,'JEOK5',mpi_rank,omp_rank
!print *,'JEOK6',mpi_rank,omp_rank

endif

  !CALL abort_gcm('initdustemission', 'OK1',1)

! Lecture des donnees du LISA indiquant le type de sol utilise.
       ! in lines: landuse
! nats: number of mineral types (14 here for sand, silt, clay etc.)
! data format in columns
! mmd1  sigma1  p1  mmd2  sigma2  p2 mmd3 ... alpha
!print *,'JEOK7',mpi_rank,omp_rank
!$OMP MASTER
IF (is_mpi_root .AND. is_omp_root) THEN
!print *,'JEOK9',mpi_rank,omp_rank
 fnsolspe='SOILSPEC.data'
  PRINT*,'  o Reading ',fnsolspe(1:40)
  OPEN(10,file=fnsolspe)
  READ(10,*)line
  DO i=1,nats
     READ(10,*)line
     READ(10,*)(solspe(i,j),j=1,nspe)
!JE     alfa(i)=solspe(i,10)
!     PRINT*,'i,alfa(i)',i,alfa(i)
!     WRITE(18,*)i,alfa(i)
  END DO
!     print*,'solspe(14,10)= ',solspe(14,10)
  CLOSE(10)
ENDIF
!$OMP END MASTER
!$OMP BARRIER
!print *,'JEOK10',mpi_rank,omp_rank
  call bcast(solspe)
! Calcul de la distribution en taille des particules de Dust
! Elle depend du nombre de classe des particules nclass.
!c making full resolution soil diameter table
      dstep=0.0005
      dp=dmin
      do i=1,nclass
         dp=dp*exp(dstep)
         sizeclass(i)=dp
         if(dp.ge.dmax+eps)goto 30
         newstep(i)=dstep
        ! WRITE(18,*)i,sizeclass(i)
      enddo
30   continue
      print*,'IK5'
      ncl=i-1
      print*,'   soil size classes used   ',ncl,' / ',nclass
      print*,'   soil size min: ',sizeclass(1),' soil size max: ',sizeclass(ncl)
      if(ncl.gt.nclass)stop

! Threshold velocity:
if (.false.) then
!if (.true.) then
!c 0: Iversen and White 1982
       print *,'Using  Iversen and White 1982 Uth'
         do i=1,ncl
            bb=adust*(sizeclass(i)**xdust)+bdust
            cc=sqrt(1+ddust*(sizeclass(i)**(-2.5)))
            xk=sqrt(abs(rop*gravity*sizeclass(i)/roa))
            if (bb.lt.10.) then
               dd=sqrt(1.928*(bb**0.092)-1.)
               uth(i)=0.129*xk*cc/dd
            else
               ee=-0.0617*(bb-10.)
               ff=1.-0.0858*exp(ee)
               uth(i)=0.12*xk*cc*ff
            endif
         enddo
endif
if(.true.) then
! 1: Shao and Lu 2000
       print *,'Using  Shao and Lu 2000 Uth'
            an=0.0123
            gam=0.3
            do i=1,ncl
               sigshao=rop/roa
               x1=sigshao*gravity*sizeclass(i)
               x2=gam/(roa*sizeclass(i))
               uth(i)=sqrt(an*(x1+x2))
            enddo
endif



!Calcul de la surface couverte par les particules de diametre dp
!srel(nats,nclass)

!  OPEN(18,file='srel.dat',form='formatted',access='sequential')
  DO ns=1,nats
     stotale=0.
     DO i=1,ncl
        su=0.
        DO j=1,nmode
           nd=((j-1)*3)+1
           nsi=((j-1)*3)+2
           npi=((j-1)*3)+3
           IF (solspe(ns,nd).EQ.0.)THEN
              su_loc=0.
           ELSE
              xk=solspe(ns,npi)/(sqrt(2.*pi)*log(solspe(ns,nsi)))
              xl=((log(sizeclass(i))-log(solspe(ns,nd)))**2) &
     &              /(2.*(log(solspe(ns,nsi)))**2)
              xm=xk*exp(-xl)
              xn=rop*(2./3.)*(sizeclass(i)/2.)
              su_loc=xm*newstep(i)/xn
           END IF
           su=su+su_loc
        END DO
        subsoildist(i)=su
        stotale=stotale+su
     END DO
     DO i=1,ncl
        IF (subsoildist(i).gt.0..and.stotale.gt.0.)THEN
            srel(ns,i)=subsoildist(i)/stotale

        ELSE
            srel(ns,i)=0.
        END IF
!        WRITE(18,*)i,srel(1,i)
     END DO
!     PRINT*,'ns , srel(10000) ',ns,srel(ns,10000)
  END DO


! Calcul de la repartition de l energie feff(klon,ntyp)

  !efficient fraction for the friction velocities as a function
  !of z0s, Zo1 et Zo2 to retrieve the apparent threshold
  !wind friction velocity.
      !    feff(:,:)=0.
       do i=1,klon
          do k=1,ntyp
     !     print*,'IKKK ',i,klon,k,ntyp
             if (zos(i,k).eq.0..or.z01(i,k).eq.0.) then
     !       if (zos(i,k)<=0..or.z01(i,k)<=0.) then
!              if (zos(i,k)<0..or.z01(i,k)<0.) then
     !            print*,'INI DUST WARNING zos ou z01<0',zos(i,k),z01(i,k)
!              endif
              feff(i,k)=0.
              feffdbg(i,k)=0.
 !         print*,'IKKK A ',i,klon,k,ntyp
            else
! drag partition betzeen the erodable surface and zo1
     !     print*,'IKKK B0 ',i,klon,k,ntyp,z01(i,k),zos(i,k),xeff,aeff
              aa=log(z01(i,k)/zos(i,k))
              tmp1(i,k)=aa
              bb=log(aeff*(xeff/zos(i,k))**0.8)
              cc=1.-aa/bb
              feffdbg(i,k)=cc
       !   print*,'IKKK B1 ',i,klon,k,ntyp
! drag partition between zo1 and zo2
! feff: total efficient fraction
              if(D(i,k).eq.0.)then
                 feff(i,k)=cc
   !       print*,'IKKK C ',i,klon,k,ntyp
              else
                 dd=log(z02(i,k)/z01(i,k))
                 ee=log(aeff*(D(i,k)/z01(i,k))**0.8)
                 feff(i,k)=(1.-dd/ee)*cc
   !       print*,'IKKK D ',i,klon,k,ntyp
              endif
              if (feff(i,k).lt.0.)feff(i,k)=0.
              if (feffdbg(i,k).lt.0.)feffdbg(i,k)=0.
              if (feff(i,k).gt.1.)feff(i,k)=1.
              if (feffdbg(i,k).gt.1.)feffdbg(i,k)=1.
    !      print*,'IKKK E ',i,klon,k,ntyp
            endif
          enddo
        enddo
! JE20150120<<
  if (flag_feff .eq. 0) then
    print *,'JE_dbg FORCED deactivated feff'
    do i=1,klon
      do k=1,ntyp
        feff(i,k)=1.
      enddo
    enddo
  endif
! JE20150120>>


if (1==1) then
!  !  CALL writefield_phy("AA",tmp1(1:klon,1:5),5)
!
    CALL writefield_phy("REPART5",feff(1:klon,1:5),5)
    CALL writefield_phy("REPART5dbg",feffdbg(1:klon,1:5),5)
endif


!c---------------FOR def_prepag01modif(var3a,var3b)-----------------------
        p3t=0.0001
        var1=e3**(1./3.)
        var2=(rop*pi)**(1./3.)
        var3a=trctunt*var1/var2
        et=e1+(sqrt(p3t*(e3*e3+e1*e2-e2*e3-e1*e3))/p3t)
        var1=et**(1./3.)
        var2=(rop*pi)**(1./3.)
        var3b=trctunt*var1/var2

! JE20140902: comment all isograd distributions and make own massfrac:


!  if (.false.) then
!!**************L718
!
!!c------------------------------------------------------------------------
!! isolog distrib and massfrac calculations.
!
!      nbinsout=nbins+1
!      b1=log(sizedustmin)
!      b2=log(sizedustmax)
!! restricted ISOLOG bins distributions
!
!!      step=(b2-b1)/(nbinsout-1)
!!      DO ni=1,nbinsout
!!         itv(ni)=exp(b1+(ni-1)*step)
!!      ENDDO
!!OPEN(18,file='vdepothrsizbin.dat',form='formatted',access='sequential')
!! Restricted ISOGRADIENT bins distributions
!!!!!!!Making HIGH RESOLUTION dust size distribution (in um)!!!!!!
!  stepbin=(b2-b1)/(nbinsHR-1)
!  DO nb=1,nbinsHR
!     binsHR(nb)=exp(b1+(nb-1)*stepbin)
!  END DO
!
!  DO nb=1,nbinsHR
!     binsHRcm(nb)=binsHR(nb)*1.e-4
!  END DO
!! Making HIGH RESOLUTION dry deposition velocity
!     CALL calcvd(vdout)
!
!
! DO nb=1,nbinsHR
!     vdHR(nb)=vdout(nb)
!!  WRITE(18,*),binsHR(nb),vdHR(nb)
!  END DO
!
!   !searching for minimum value of dry deposition velocity
!  minisograd=1.e20
!  DO nb=1,nbinsHR
!     IF(vdHR(nb).le.minisograd)THEN
!       minisograd=vdHR(nb)
!       miniso=nb
!     END IF
!  END DO
!
!! searching for optimal number of bins in positive slope Vd part
!
!  nbins1=1
!  nbins2=nbinsout-1
!     DO k=1,1000
!        delta1=abs(log(vdHR(1))-log(vdHR(miniso)) )/nbins1
!        delta2=abs(log(vdHR(nbinsHR))-log(vdHR(miniso)))/(nbins2-1)
!        IF(delta2.GE.delta1)THEN
!        GOTO 50
!
!        ELSE
!           nbins2=nbins2-1
!           nbins1=nbins1+1
!        END IF
!        IF(nbins2.EQ.1)THEN
!           PRINT*,'!! warning: lower limit was found: STOP'
!           STOP
!        END IF
!     END DO
! 50  CONTINUE
!print*,'IK10'
!! building the optimized distribution
!  logvdISOGRAD(1)=log(vdHR(1))
!  binsISOGRAD(1)=binsHR(1)
!  DO k=2,nbins1
!     logvdISOGRAD(k)=logvdISOGRAD(1)-(k-1)*delta1
!  END DO
!
!  logvdISOGRAD(nbins1+1)=log(minisograd)
!
!  DO k=1,nbins2
!     logvdISOGRAD(nbins1+1+k)=logvdISOGRAD(nbins1+1)+k*delta2
!  END DO
!
!  DO k=1,nbinsout
!     vdISOGRAD(k)=exp(logvdISOGRAD(k))
!  END DO
!! sequential search of the correspondig particle diameter
!  istart=1
!  DO k=2,nbinsout-1
!     curvd=vdISOGRAD(k)
!     DO nb=istart,nbinsHR
!        avdhr=vdHR(nb)
!        bvdhr=vdHR(nb+1)
!        IF(nb.LT.(miniso-1).AND.curvd.LT.avdhr.AND. &
!           curvd.GT.bvdhr)THEN
!           binsISOGRAD(k)=binsHR(nb)
!           istart=nb
!           GOTO 60
!        END IF
!        IF(nb.eq.miniso)THEN
!           binsISOGRAD(k)=binsHR(nb)
!           istart=nb+1
!           GOTO 60
!        END IF
!        IF(nb.GT.miniso.AND.curvd.GT.avdhr.AND. &
!                              curvd.LT.bvdhr)THEN
!           binsISOGRAD(k)=binsHR(nb)
!           istart=nb
!           GOTO 60
!        END IF
!     END DO
! 60  CONTINUE
!  END DO
!print*,'IK11'
!  binsISOGRAD(nbinsout)=binsHR(nbinsHR)
!  vdISOGRAD(nbinsout)=vdHR(nbinsHR)
!  DO nb=1,nbinsout
!     itv(nb)=binsISOGRAD(nb)
!  END DO
!  end ! JE20140902 if false

! Making dust size distribution (in um)
!
      nbinsout=nbins+1
      b1=log(sizedustmin)
      b2=log(sizedustmax)
      stepbin=(b2-b1)/(nbinsout-1)
     DO ni=1,nbinsout
        itv(ni)=exp(b1+(ni-1)*stepbin)
     ENDDO

!  stepbin=(b2-b1)/(nbinsHR-1)
!  DO nb=1,nbinsHR
!     binsHR(nb)=exp(b1+(nb-1)*stepbin)
!  END DO
!
!  DO nb=1,nbinsHR
!     binsHRcm(nb)=binsHR(nb)*1.e-4
!  END DO



      DO nb=1,nbins
         ldmd=(log(itv(nb+1))-log(itv(nb)))/2.
         sizedust(nb)=exp(log(itv(nb))+ldmd)
        PRINT*,nb, itv(nb),' < ',sizedust(nb),' < ',itv(nb+1)
      ENDDO
! Making dust size distribution (in cm) ???
      DO nb=1,nbins
         szdcm(nb)=sizedust(nb)*1.e-4
      ENDDO
!c preparation of emissions reaffectation.
      DO k=1,ndistb
         exden=sqrt(2.)*log(sig(k))
         DO nb=1,nbins
            d2min=itv(nb)
            d2max=itv(nb+1)
            numin=log(d2min/diam(k))
            numax=log(d2max/diam(k))
            massfrac(k,nb)=0.5*(erf(numax/exden)-erf(numin/exden))
            !print *,k,nb,massfrac(k,nb)
         ENDDO
      ENDDO

!$OMP MASTER
IF (is_mpi_root .AND. is_omp_root) THEN
      OPEN (unit=15001, file='massfrac')
      DO k=1,ndistb
        DO nb=1,nbins
            write(15001,*),k,nb,massfrac(k,nb)
        ENDDO
      ENDDO
ENDIF
!$OMP END MASTER
!$OMP BARRIER

 !CALL abort_gcm('calcdustemission', 'OK1',1)

  !------------


  END SUBROUTINE initdust

!--------------------------------------------------------------------------------------
!======================================================================================
!**************************************************************************************
!======================================================================================
!--------------------------------------------------------------------------------------

  SUBROUTINE calcdustemission(debutphy,zu10m,zv10m,wstar, &
                              ale_bl,ale_wake,param_wstarBL,param_wstarWAKE, &
                              emisbin)
  ! emisions over 12 dust bin
  USE dimphy
  USE infotrac

  IMPLICIT NONE
  ! Input
  LOGICAL, INTENT(IN)                   :: debutphy ! First physiqs run or not
  REAL,DIMENSION(klon),INTENT(IN)          :: zu10m   ! 10m zonal wind
  REAL,DIMENSION(klon),INTENT(IN)          :: zv10m   ! meridional 10m wind
  REAL,DIMENSION(klon),INTENT(IN)          :: wstar   
  REAL,DIMENSION(klon),INTENT(IN)          :: ale_bl
  REAL,DIMENSION(klon),INTENT(IN)          :: ale_wake
  
  ! Local variables
!  INTEGER, PARAMETER :: flag_wstar=150
!  INTEGER, PARAMETER :: flag_wstar=120
!  INTEGER, PARAMETER :: flag_wstar=125
  REAL,DIMENSION(klon), INTENT(IN) :: param_wstarWAKE
  REAL,DIMENSION(klon), INTENT(IN) :: param_wstarBL
  REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: fluxdust ! horizonal emission fluxes in UNITS for the nmod soil aerosol modes
  REAL,DIMENSION(:), ALLOCATABLE,SAVE   :: wind10ms   ! 10m wind distribution in m/s
  REAL,DIMENSION(:), ALLOCATABLE,SAVE   :: wind10cm   ! 10m wind distribution in cm/s
  REAL,DIMENSION(klon)                  :: zwstar   
  REAL,DIMENSION(nwb)                :: probu
!  REAL, DIMENSION(nmode) :: fluxN,ftN,adN,fdpN,pN,eN ! in the original code N=1,2,3
  REAL :: flux1,flux2,flux3,ft1,ft2,ft3


  REAL, PARAMETER        :: umin=21.
  REAL, PARAMETER        :: woff=0.5  ! min value of 10m wind speed accepted for emissions
  REAL :: pdfcum,U10mMOD,pdfu,weilambda
  REAL :: z0salt,ceff,cerod,cpcent
  REAL :: cdnms,ustarns,modwm,utmin
!JE20150202   REAL :: cdnms,ustarns,modwm
  REAL :: fdp1,fdp2,ad1,ad2,ad3,flux_diam
  REAL :: dfec1,dfec2,dfec3,t1,t2,t3,p1,p2,p3,dec,ec
  ! auxiliar counters
  INTEGER                               :: kwb
  INTEGER                               :: i,j,k,l,n
  INTEGER  :: kfin,ideb,ifin,kfin2,istep
  REAL :: auxreal

  ! Output
  !REAL,DIMENSION(:,:), ALLOCATABLE,SAVE  :: emisbin ! vertical emission fluxes in UNITS for the 12 bins
  REAL,DIMENSION(klon,nbins)  :: emisbin ! vertical emission fluxes in UNITS for the 12 bins
!$OMP THREADPRIVATE(fluxdust) 
!$OMP THREADPRIVATE(wind10ms) 
!$OMP THREADPRIVATE(wind10cm) 

  !----------------------------------------------------
  ! initialization
  !----------------------------------------------------
  IF (debutphy) THEN
!   ALLOCATE( emisbin(klon,nbins) )
   ALLOCATE( fluxdust(klon,nmode) )
   ALLOCATE( wind10ms(nwb) )
   ALLOCATE( wind10cm(nwb) )
  ENDIF !debutphy

  
  DO i=1,klon
   DO j=1,nbins
    emisbin(i,j)  = 0.0 
   ENDDO !j,nbind
   DO j=1,nmode
    fluxdust(i,j) = 0.0
   ENDDO !j,nmode
  zwstar(i) = 0.0
  ENDDO !i,klon
  !----------------------------------------------------
  ! calculations
  !----------------------------------------------------

  ! including BL processes..
!JE201041120  zwstar=sqrt(2.*(ale_bl+0.01*(flag_wstar-100)*ale_wake))
!JE20141124  zwstar=sqrt(2.*(flag_wstarBL*ale_bl+0.01*(flag_wstar-100)*ale_wake))
!  print
!*,'zwstar=sqrt(2.*(',flag_wstarBL,'ale_bl+0.01*(',flag_wstar,'-100)*ale_wake))'
  !
    DO i=1,klon  ! main loop
     zwstar(i)=sqrt(2.*(param_wstarBL(i)*ale_bl(i)+param_wstarWAKE(i)*ale_wake(i)))
     U10mMOD=MAX(woff,sqrt(zu10m(i)*zu10m(i)+zv10m(i)*zv10m(i)))
     pdfcum=0.
     ! Wind weibull distribution:

           DO kwb=1,nwb
                flux1=0.
                flux2=0.
                flux3=0.
!JE20141205<< differences in weibull distribution: expectance of the distribution is
!0.89*U10mMODD instead of U10mMOD
! now: lambda parameter of weibull ditribution is estimated as
! lambda=U10mMOD/gamma(1+1/kref)
! gamma function estimated with stirling formula
                auxreal=1.+1./kref
                weilambda = U10mMOD/exp(auxreal*log(auxreal)-auxreal &
                         - 0.5*log(auxreal/(2.*pi))+1./(12.*auxreal) &
                         -1./(360.*(auxreal**3.))+1./(1260.*(auxreal**5.)))
                IF(nwb.gt.1)THEN
                   wind10ms(kwb)=kwb*2.*U10mMOD/nwb
!original
!                   pdfu=(kref/U10mMOD)*(wind10ms(kwb)/U10mMOD)**(kref-1) &
!                      *exp(-(wind10ms(kwb)/U10mMOD)**kref)
                   pdfu=(kref/weilambda)*(wind10ms(kwb)/weilambda)**(kref-1) &
                      *exp(-(wind10ms(kwb)/weilambda)**kref)
!                   !print *,'JEdbg  U10mMOD weilambda  ',U10mMOD,weilambda
!JE20141205>>

                   probu(kwb)=pdfu*2.*U10mMOD/nwb
                   pdfcum=pdfcum+probu(kwb)
                      IF(probu(kwb).le.1.e-2)GOTO 70
                ELSE
                   wind10ms(kwb)=U10mMOD
                   probu(kwb)=1.
                ENDIF
             wind10cm(kwb)=wind10ms(kwb)*100.
             DO n=1,ntyp
                   ft1=0.
                   ft2=0.
                   ft3=0.
!JE20150129<<<<

             IF(.FALSE.) THEN
!                  nat=int(sol(i,n))
!                    print *,i,n
                    IF(sol(i,n).gt.1..and.sol(i,n).lt.15.) nat=int(sol(i,n))
!JE20140526<<
!                    print *,'JE: WARNING: nat=0 forced to nat=99!! and doing nothing'
                   IF(sol(i,n).lt.0.5) THEN
                      nat=99
                      GOTO 80
                    ENDIF
!JE20140526>>
 

                 !IF(n.eq.1.and.nat.eq.99)GOTO 80
             !      if(n.eq.1) print*,'nat1=',nat,'sol1=',sol(i,n)
                   IF(n.eq.1.and.nat.eq.99)GOTO 80

             ENDIF
             IF(.TRUE.) THEN
                nat=int(sol(i,n))
                if(n == 1 .and. nat >= 14 .or. nat < 1 .or. nat > 19) GOTO 80
             ENDIF
!JE20150129>>>>

                      z0salt=z02(i,n)
                      ceff=feff(i,n)
                      cerod=A(i,n)
                      cpcent=P(i,n)
                      ustarsalt=0.
                   IF(ceff.le.0..or.z0salt.eq.0.)GOTO 80
                   IF(cerod.eq.0.or.cpcent.eq.0.)GOTO 80
! in cm: utmin, umin, z10m, z0salt, ustarns
! in meters: modwm
! no dimension for: cdnms
! Cas ou wsta=0.
                      cdnms=vkarm/(log(z10m/z0salt))
                      modwm=sqrt((wind10ms(kwb)**2)+(1.2*zwstar(i))**2)
                      ustarns=cdnms*modwm*100.
!JE20150202 <<
! Do not have too much sense.. and is not anymore in the chimere14b version.
!
!                      utmin=umin/(cdnms*ceff)
!                   IF(wind10cm(kwb).ge.utmin)THEN
!                      ustarsalt=ustarns+  &
!                    (0.3*(wind10cm(kwb)/100.-utmin/100.)**2.)
!                   ELSE
!                      ustarsalt=ustarns
!                   ENDIF
! ustarsalt should be :
                    ustarsalt=ustarns
!JE20150202 >>


                   IF(ustarsalt.lt.umin/ceff)GOTO 80
!                      print*,'ustarsalt = ',ustarsalt
!----------------------------------------
                    CALL def_copyncl(kfin)

!                  CALL def_ag01(kfin,ft1,ft2,ft3)
      do ni=1,kfin
         fdp1=1.-(uth2(ni)/(ceff*ustarsalt))
         if (fdp1.le.0..or.srel2(nat,ni).eq.0.) then
            ad1=0.
            ad2=0.
            ad3=0.
         else
       
            fdp2=(1.+(uth2(ni)/(ceff*ustarsalt)))**2.
            flux_diam=cd*srel2(nat,ni)*fdp1*fdp2*(ustarsalt**3)
            dec=flux_diam*beta
            ec=(pi/3.)*1.e+2*rop* (sizeclass2(ni)**3) *(ustarsalt**2)
            dfec1=(ec-e1)
            dfec2=(ec-e2)
            dfec3=(ec-e3)
            t1=0.
            t2=0.
            t3=0.
            if(ec.ge.e1)t1=1.
            if(ec.ge.e2)t2=1.
            if(ec.ge.e3)t3=1.
            if(dfec3.ne.0.)then
               p1=t1*dfec1/dfec3
               p2=t2*(1.-p1)*dfec2/dfec3
            else
               p1=0.
               p2=0.
            endif
            p3=t3*(1.-p1-p2)
            ad1=(pi/6.)*dec*rop*p1*((diam(1)*1.e-4)**3.)/e1
            ad2=(pi/6.)*dec*rop*p2*((diam(2)*1.e-4)**3.)/e2
            ad3=(pi/6.)*dec*rop*p3*((diam(3)*1.e-4)**3.)/e3

         endif
         ft1=ft1+ad1
         ft2=ft2+ad2
         ft3=ft3+ad3
      enddo ! of loop over nc

!....

        flux1=flux1+ft1*cpcent*cerod
        flux2=flux2+ft2*cpcent*cerod
        flux3=flux3+ft3*cpcent*cerod
!       print *,'JEflux :',kwb,n,flux1,flux2,flux3
80 CONTINUE
             ENDDO !n=1,ntyp
70 CONTINUE
        fluxdust(i,1)=fluxdust(i,1)+flux1*probu(kwb)
        fluxdust(i,2)=fluxdust(i,2)+flux2*probu(kwb)
        fluxdust(i,3)=fluxdust(i,3)+flux3*probu(kwb)
   ENDDO !kwb=1,nwb
      m1dflux(i)=10.*fluxdust(i,1)
      m2dflux(i)=10.*fluxdust(i,2)          ! tous en Kg/m2/s
      m3dflux(i)=10.*fluxdust(i,3)



    ENDDO  ! i, klon



   
  !from horizontal to  vertical fluxes for each bin
  DO k=1,klon
   DO i=1,ndistb
    DO j=1,nbins
!JE20150202 <<
    emisbin(k,j) = emisbin(k,j)+10*fluxdust(k,i)*massfrac(i,j)  
!JE20150202 >>
    ENDDO !j, nbind
   ENDDO  !i, nmode
  ENDDO !k,klon


!print *,' JE fluxdust in calcdust'
! DO k=1,klon
!   DO i=1,ndistb
!!print *,k,i,fluxdust(k,i)
!enddo
!enddo
!print *,' JE emisbin in calcdust'
!do k=1,klon
!do j=1,nbins
!!print *,k,j,emisbin(k,j)
!enddo
!enddo 
!  CALL abort_gcm('calcdustemission', 'OK1',1)

  END SUBROUTINE calcdustemission
!--------------------------------------------------------------------------------------
!======================================================================================
!**************************************************************************************
!======================================================================================
!--------------------------------------------------------------------------------------



SUBROUTINE def_copyncl(kfin)
      implicit none

 integer i,n,kfin,ideb,ifin,istep,kfin2
    real dsmin,dsmax

! estimation of the reduced soil size distribution
         dsmin=var3a*(ustarsalt**(-2./3.))
         dsmax=var3b*(ustarsalt**(-2./3.))
!      print*,'ustarsalt = ',ustarsalt,'dsmin=',dsmin,'dsmax=',dsmax
! dichotomy
         call def_dichotomy(sizeclass,nclass,1,ncl,dsmin,ideb)
   !      print*,'ideb = ',ideb
         call def_dichotomy(sizeclass,nclass,ideb,ncl,dsmax,ifin)
   !      print*,'ifin = ',ifin
! readaptation of large sizes particles
         kfin=0
         do i=ideb,ifin
            kfin=kfin+1
            sizeclass2(kfin)=sizeclass(i)
            uth2(kfin)=uth(i)
            srel2(nat,kfin)=srel(nat,i)
         enddo
  !          print*,'je suis la'
         kfin2=kfin
         istep=50
         do i=ifin,ncl,istep
            kfin=kfin+1
            sizeclass2(kfin)=sizeclass(i)
            uth2(kfin)=uth(i)
            srel2(nat,kfin)=srel(nat,i)*istep
         enddo
         if(kfin.ge.nclass)then
            print*,'$$$$ Tables dimension problem:',kfin,'>',nclass
         endif
!---------------        

      return
END SUBROUTINE def_copyncl

!--------------------------------------------------------------------------------------
!======================================================================================
!**************************************************************************************
!======================================================================================
!--------------------------------------------------------------------------------------

subroutine def_dichotomy(siz,nclass,i1,i2,ds,iout)
!c---------------------------------------------------------------
!c 'size' is the table to scan
!c      of dimension 'nclass', and reduced size of interest [i1:i2]
!c 'ds' is the value to find in 'size'
!c 'iout' is the number in the table 'size' correspondig to 'ds'
!c---------------------------------------------------------------

      implicit none
      integer i1,i2,nclass,iout,ismin,ismax,k2,ihalf,idiff
      real siz(nclass),ds
!c-----------------------------
      iout=0
      ismin=i1
      ismax=i2
      ihalf=int((ismax+ismin)/2.)
      do k2=1,1000000
          if(ds.gt.siz(ihalf))then
             ismin=ihalf
          else
             ismax=ihalf
          endif
          ihalf=int((ismax+ismin)/2.)
          idiff=ismax-ismin
          if(idiff.le.1)then
             iout=ismin
             goto 52
          endif
      enddo
 52   continue
      if(iout.eq.0)then
        print*,'$$$$ Tables dimension problem: ',iout
      endif

      end subroutine def_dichotomy

!--------------------------------------------------------------------------------------
!======================================================================================
!**************************************************************************************
!======================================================================================
!--------------------------------------------------------------------------------------


  SUBROUTINE calcvd(vdout)
  IMPLICIT NONE
  INTEGER                     :: nb
!JE already def in module  INTEGER,PARAMETER           :: nbinsHR=3000
  INTEGER,PARAMETER           :: idiffusi=1
  INTEGER,PARAMETER           :: idrydep=2
  REAL                        :: dmn1,dmn2,dmn3,ra,rb
  REAL                        :: rhod,muair,nuair,R,airmolw
  REAL                        :: bolz,rhoair,gravity,karman,zed,z0m
  REAL                        :: ustarbin,temp,pres,rac,lambda,ccexp
  REAL                        :: Cc,rexp,pi,St
  REAL,DIMENSION(nbinsHR)     :: setvel
  REAL,DIMENSION(nbinsHR)     :: diffmol1
  REAL,DIMENSION(nbinsHR)     :: diffmol2
  REAL,DIMENSION(nbinsHR)     :: schmidtnumb
  REAL,DIMENSION(nbinsHR)     :: diffmole

  REAL,DIMENSION(nbinsHR),   INTENT(OUT)              :: vdout
! users physical parameters
       temp=273.15+25. ! in Kelvin degrees
       pres=1013. ! in hPa
! calculation of Ra
       zed=10.
       z0m=0.05
       karman=0.41
       ustarbin=30. ! cm/s
       ra=log(zed/z0m)/(karman*ustarbin)
!c physical constants (all must in units g / cm / s)
       rhod=2.65                 ! g.cm-3
       muair=1.72e-4             ! g.cm-1.s-1
       nuair=0.1461              ! air kinematic viscosity [cm2.s-1]
       R=8.314                   ! molar gas constant J.mol-1.K-1
       airmolw=28.8              ! molar weight of air
       bolz=1.38e-16                    ! g.cm2.s-2
       rhoair=1.225e-3                  ! g.cm-3
       gravity=981.                ! gravity [cm.s-2]
       pi=3.14159
       rac=sqrt(8.0*airmolw/(pi*R*temp))
       lambda=(2.*muair)/(pres*rac)
! c Cc and Vs
!       binsHRcm=binsize*1.e-4
       DO nb=1,nbinsHR
       ccexp=exp((-0.55*binsHRcm(nb))/lambda)
       Cc=1.+(2.*lambda)/binsHRcm(nb)*(1.257+0.4*ccexp)
       setvel(nb)=(1./18.)*((binsHRcm(nb))**2*rhod*gravity*Cc)/muair
!c Molecular diffusivity (s/cm2) and Schmidt number
!c Note that the molecular diffusivity uses diameter in micro-m
!c to give a result directly in s/cm2
       dmn1=2.38e-7/binsHR(nb)
       dmn2=0.163/binsHR(nb)
       dmn3=0.0548*exp(-6.66*binsHR(nb))/binsHR(nb)
       diffmol1(nb)=dmn1*(1.+dmn2+dmn3)
       diffmol2(nb)=bolz*temp*Cc/(3.*pi*muair*binsHRcm(nb))
       IF(idiffusi.EQ.1)diffmole(nb)=diffmol1(nb)
       IF(idiffusi.EQ.2)diffmole(nb)=diffmol2(nb)
       schmidtnumb(nb)=nuair/diffmole(nb)
       St=setvel(nb)*ustarbin*ustarbin/(gravity*nuair)
       rb=1./(ustarbin*((schmidtnumb(nb))**(-2./3.)+10.**(-3./St)))
!c wesely (primarily designed for gases)
       IF(idrydep.EQ.1)THEN
          vdout(nb)=1./(ra+rb+ra*rb*setvel(nb))+setvel(nb)
       END IF
!c venkatram and pleim (more adaptated to particles but numerically unstable)
       IF(idrydep.EQ.2)THEN
        rexp=exp(-(ra+rb)*setvel(nb))
        vdout(nb)=setvel(nb)/(1.-rexp)
       END IF
       END DO
!c-----------------
       RETURN
      END SUBROUTINE calcvd



END MODULE dustemission_mod

