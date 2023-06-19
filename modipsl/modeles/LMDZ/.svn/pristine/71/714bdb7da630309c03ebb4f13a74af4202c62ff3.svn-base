!
! $Id: iniacademic.F90 1625 2012-05-09 13:14:48Z lguez $
!
SUBROUTINE iniacademic_loc(vcov,ucov,teta,q,masse,ps,phis,time_0)

  USE filtreg_mod, ONLY: inifilr
  use exner_hyb_m, only: exner_hyb
  use exner_milieu_m, only: exner_milieu
  USE infotrac, ONLY: nqtot,niso_possibles,ok_isotopes,iqpere,ok_iso_verif,tnat,alpha_ideal, &
        & iqiso,phase_num,iso_indnum,iso_num,zone_num
  USE control_mod, ONLY: day_step,planet_type
  USE parallel_lmdz, ONLY: ijb_u, ije_u, ijb_v, ije_v
#ifdef CPP_IOIPSL
  USE IOIPSL, ONLY: getin
#else
  ! if not using IOIPSL, we still need to use (a local version of) getin
  USE ioipsl_getincom, ONLY: getin
#endif
  USE Write_Field
  USE comconst_mod, ONLY: cpp, kappa, g, daysec, dtvr, pi, im, jm
  USE logic_mod, ONLY: iflag_phys, read_start
  USE comvert_mod, ONLY: ap, bp, preff, presnivs, pressure_exner
  USE temps_mod, ONLY: annee_ref, day_ini, day_ref
  USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0


  !   Author:    Frederic Hourdin      original: 15/01/93
  ! The forcing defined here is from Held and Suarez, 1994, Bulletin
  ! of the American Meteorological Society, 75, 1825.

  IMPLICIT NONE

  !   Declararations:
  !   ---------------

  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
  include "academic.h"
  include "iniprint.h"

  !   Arguments:
  !   ----------

  REAL,INTENT(OUT) :: time_0

  !   fields
  REAL,INTENT(OUT) :: vcov(ijb_v:ije_v,llm) ! meridional covariant wind
  REAL,INTENT(OUT) :: ucov(ijb_u:ije_u,llm) ! zonal covariant wind
  REAL,INTENT(OUT) :: teta(ijb_u:ije_u,llm) ! potential temperature (K)
  REAL,INTENT(OUT) :: q(ijb_u:ije_u,llm,nqtot) ! advected tracers (.../kg_of_air)
  REAL,INTENT(OUT) :: ps(ijb_u:ije_u) ! surface pressure (Pa)
  REAL,INTENT(OUT) :: masse(ijb_u:ije_u,llm) ! air mass in grid cell (kg)
  REAL,INTENT(OUT) :: phis(ijb_u:ije_u) ! surface geopotential

  !   Local:
  !   ------

  REAL,ALLOCATABLE :: vcov_glo(:,:),ucov_glo(:,:),teta_glo(:,:)
  REAL,ALLOCATABLE :: q_glo(:,:),masse_glo(:,:),ps_glo(:)
  REAL,ALLOCATABLE :: phis_glo(:)
  REAL p (ip1jmp1,llmp1  )               ! pression aux interfac.des couches
  REAL pks(ip1jmp1)                      ! exner au  sol
  REAL pk(ip1jmp1,llm)                   ! exner au milieu des couches
  REAL phi(ip1jmp1,llm)                  ! geopotentiel
  REAL ddsin,zsig,tetapv,w_pv  ! variables auxiliaires
  real tetastrat ! potential temperature in the stratosphere, in K
  real tetajl(jjp1,llm)
  INTEGER i,j,l,lsup,ij

  REAL teta0,ttp,delt_y,delt_z,eps ! Constantes pour profil de T
  REAL k_f,k_c_a,k_c_s         ! Constantes de rappel
  LOGICAL ok_geost             ! Initialisation vent geost. ou nul
  LOGICAL ok_pv                ! Polar Vortex
  REAL phi_pv,dphi_pv,gam_pv   ! Constantes pour polar vortex 

  real zz,ran1
  integer idum

  REAL zdtvr
  
  character(len=*),parameter :: modname="iniacademic"
  character(len=80) :: abort_message

  ! Sanity check: verify that options selected by user are not incompatible
  if ((iflag_phys==1).and. .not. read_start) then
    write(lunout,*) trim(modname)," error: if read_start is set to ", &
    " false then iflag_phys should not be 1"
    write(lunout,*) "You most likely want an aquaplanet initialisation", &
    " (iflag_phys >= 100)"
    call abort_gcm(modname,"incompatible iflag_phys==1 and read_start==.false.",1)
  endif
  
  !-----------------------------------------------------------------------
  ! 1. Initializations for Earth-like case
  ! --------------------------------------
  !
  ! initialize planet radius, rotation rate,...
  call conf_planete

  time_0=0.
  day_ref=1
  ! annee_ref=0

  im         = iim
  jm         = jjm
  day_ini    = 1
  dtvr    = daysec/REAL(day_step)
  zdtvr=dtvr
  etot0      = 0.
  ptot0      = 0.
  ztot0      = 0.
  stot0      = 0.
  ang0       = 0.      

  if (llm == 1) then
     ! specific initializations for the shallow water case
     kappa=1
  endif

  CALL iniconst
  CALL inigeom
  CALL inifilr

  if (llm == 1) then
     ! initialize fields for the shallow water case, if required
     if (.not.read_start) then
        phis(ijb_u:ije_u)=0.
        q(ijb_u:ije_u,1:llm,1:nqtot)=0
        CALL sw_case_williamson91_6_loc(vcov,ucov,teta,masse,ps)
     endif
  endif

  academic_case: if (iflag_phys >= 2) then
     ! initializations

     ! 1. local parameters
     ! by convention, winter is in the southern hemisphere
     ! Geostrophic wind or no wind?
     ok_geost=.TRUE.
     CALL getin('ok_geost',ok_geost)
     ! Constants for Newtonian relaxation and friction
     k_f=1.                !friction 
     CALL getin('k_j',k_f)
     k_f=1./(daysec*k_f)
     k_c_s=4.  !cooling surface
     CALL getin('k_c_s',k_c_s)
     k_c_s=1./(daysec*k_c_s)
     k_c_a=40. !cooling free atm
     CALL getin('k_c_a',k_c_a)
     k_c_a=1./(daysec*k_c_a)
     ! Constants for Teta equilibrium profile
     teta0=315.     ! mean Teta (S.H. 315K)
     CALL getin('teta0',teta0)
     ttp=200.       ! Tropopause temperature (S.H. 200K)
     CALL getin('ttp',ttp)
     eps=0.         ! Deviation to N-S symmetry(~0-20K)
     CALL getin('eps',eps)
     delt_y=60.     ! Merid Temp. Gradient (S.H. 60K)
     CALL getin('delt_y',delt_y)
     delt_z=10.     ! Vertical Gradient (S.H. 10K)
     CALL getin('delt_z',delt_z)
     ! Polar vortex
     ok_pv=.false.
     CALL getin('ok_pv',ok_pv)
     phi_pv=-50.            ! Latitude of edge of vortex
     CALL getin('phi_pv',phi_pv)
     phi_pv=phi_pv*pi/180.
     dphi_pv=5.             ! Width of the edge
     CALL getin('dphi_pv',dphi_pv)
     dphi_pv=dphi_pv*pi/180.
     gam_pv=4.              ! -dT/dz vortex (in K/km)
     CALL getin('gam_pv',gam_pv)

     ! 2. Initialize fields towards which to relax
     ! Friction
     knewt_g=k_c_a
     DO l=1,llm
        zsig=presnivs(l)/preff
        knewt_t(l)=(k_c_s-k_c_a)*MAX(0.,(zsig-0.7)/0.3)
        kfrict(l)=k_f*MAX(0.,(zsig-0.7)/0.3)
     ENDDO
     DO j=1,jjp1
        clat4((j-1)*iip1+1:j*iip1)=cos(rlatu(j))**4
     ENDDO

     ! Potential temperature 
     DO l=1,llm
        zsig=presnivs(l)/preff
        tetastrat=ttp*zsig**(-kappa)
        tetapv=tetastrat
        IF ((ok_pv).AND.(zsig.LT.0.1)) THEN
           tetapv=tetastrat*(zsig*10.)**(kappa*cpp*gam_pv/1000./g)
        ENDIF
        DO j=1,jjp1
           ! Troposphere
           ddsin=sin(rlatu(j))
           tetajl(j,l)=teta0-delt_y*ddsin*ddsin+eps*ddsin &
                -delt_z*(1.-ddsin*ddsin)*log(zsig)
           if (planet_type=="giant") then
             tetajl(j,l)=teta0+(delt_y*                   &
                ((sin(rlatu(j)*3.14159*eps+0.0001))**2)   &
                / ((rlatu(j)*3.14159*eps+0.0001)**2))     &
                -delt_z*log(zsig)
           endif
           ! Profil stratospherique isotherme (+vortex)
           w_pv=(1.-tanh((rlatu(j)-phi_pv)/dphi_pv))/2.
           tetastrat=tetastrat*(1.-w_pv)+tetapv*w_pv             
           tetajl(j,l)=MAX(tetajl(j,l),tetastrat)  
        ENDDO
     ENDDO

     !          CALL writefield('theta_eq',tetajl)

     do l=1,llm
        do j=1,jjp1
           do i=1,iip1
              ij=(j-1)*iip1+i
              tetarappel(ij,l)=tetajl(j,l)
           enddo
        enddo
     enddo

     ! 3. Initialize fields (if necessary)
     IF (.NOT. read_start) THEN
       ! allocate global fields:
!       allocate(vcov_glo(ip1jm,llm))
       allocate(ucov_glo(ip1jmp1,llm))
       allocate(teta_glo(ip1jmp1,llm))
       allocate(ps_glo(ip1jmp1))
       allocate(masse_glo(ip1jmp1,llm))
       allocate(phis_glo(ip1jmp1))

        ! surface pressure
        if (iflag_phys>2) then
           ! specific value for CMIP5 aqua/terra planets
           ! "Specify the initial dry mass to be equivalent to
           !  a global mean surface pressure (101325 minus 245) Pa."
           ps_glo(:)=101080.  
        else
           ! use reference surface pressure
           ps_glo(:)=preff
        endif
        
        ! ground geopotential
        phis_glo(:)=0.

        CALL pression ( ip1jmp1, ap, bp, ps_glo, p       )
        if (pressure_exner) then
          CALL exner_hyb( ip1jmp1, ps_glo, p, pks, pk )
        else
          call exner_milieu(ip1jmp1,ps_glo,p,pks,pk)
        endif
        CALL massdair(p,masse_glo)

        ! bulk initialization of temperature
        teta_glo(:,:)=tetarappel(:,:)

        ! geopotential
        CALL geopot(ip1jmp1,teta_glo,pk,pks,phis_glo,phi)

        ! winds
        if (ok_geost) then
           call ugeostr(phi,ucov_glo)
        else
           ucov_glo(:,:)=0.
        endif
        vcov(ijb_v:ije_v,1:llm)=0.

        ! bulk initialization of tracers
        if (planet_type=="earth") then
           ! Earth: first two tracers will be water

           do i=1,nqtot
              if (i == 1) q(ijb_u:ije_u,:,i)=1.e-10
              if (i == 2) q(ijb_u:ije_u,:,i)=1.e-15
              if (i.gt.2) q(ijb_u:ije_u,:,i)=0.

              ! CRisi: init des isotopes
              ! distill de Rayleigh très simplifiée
              if (ok_isotopes) then
                if ((iso_num(i).gt.0).and.(zone_num(i).eq.0)) then          
                   q(ijb_u:ije_u,:,i)=q(ijb_u:ije_u,:,iqpere(i))       &
      &                  *tnat(iso_num(i))                             &
      &                  *(q(ijb_u:ije_u,:,iqpere(i))/30.e-3)                              &
     &                   **(alpha_ideal(iso_num(i))-1)
                endif                
                if ((iso_num(i).gt.0).and.(zone_num(i).eq.1)) then
                  q(ijb_u:ije_u,:,i)=q(ijb_u:ije_u,:,iqiso(iso_indnum(i),phase_num(i)))
                endif
              endif !if (ok_isotopes) then

           enddo
        else
           q(ijb_u:ije_u,:,:)=0
        endif ! of if (planet_type=="earth")

        if (ok_iso_verif) then
           call check_isotopes(q,ijb_u,ije_u,'iniacademic_loc')
        endif !if (ok_iso_verif) then

        ! add random perturbation to temperature
        idum  = -1
        zz = ran1(idum)
        idum  = 0
        do l=1,llm
           do ij=iip2,ip1jm
              teta_glo(ij,l)=teta_glo(ij,l)*(1.+0.005*ran1(idum))
           enddo
        enddo

        ! maintain periodicity in longitude
        do l=1,llm
           do ij=1,ip1jmp1,iip1
              teta_glo(ij+iim,l)=teta_glo(ij,l)
           enddo
        enddo

        ! copy data from global array to local array:
        teta(ijb_u:ije_u,:)=teta_glo(ijb_u:ije_u,:)
        ucov(ijb_u:ije_u,:)=ucov_glo(ijb_u:ije_u,:)
!        vcov(ijb_v:ije_v,:)=vcov_glo(ijb_v:ije_v,:)
        masse(ijb_u:ije_u,:)=masse_glo(ijb_u:ije_u,:)
        ps(ijb_u:ije_u)=ps_glo(ijb_u:ije_u)
        phis(ijb_u:ije_u)=phis_glo(ijb_u:ije_u)

        deallocate(teta_glo)
        deallocate(ucov_glo)
!        deallocate(vcov_glo)
        deallocate(masse_glo)
        deallocate(ps_glo)
        deallocate(phis_glo)
     ENDIF ! of IF (.NOT. read_start)
  endif academic_case

END SUBROUTINE iniacademic_loc
