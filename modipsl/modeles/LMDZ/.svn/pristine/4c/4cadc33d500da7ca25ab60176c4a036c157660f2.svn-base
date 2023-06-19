! $Id: physiq.F 1565 2011-08-31 12:53:29Z jghattas $
MODULE physiq_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE physiq (nlon,nlev, &
     &            debut,lafin,pdtphys, &
     &            paprs,pplay,pphi,pphis,presnivs, &
     &            u,v,t,qx, &
     &            flxmass_w, &
     &            d_u, d_v, d_t, d_qx, d_ps)

      USE dimphy, only : klon,klev
      USE infotrac_phy, only : nqtot
      USE geometry_mod, only : latitude
      USE comcstphy, only : rg
      USE iophy, only : histbeg_phy,histwrite_phy
      USE ioipsl, only : getin,histvert,histdef,histend,ymds2ju
      USE mod_phys_lmdz_para, only : jj_nb
      USE phys_state_var_mod, only : phys_state_var_init
      USE mod_grid_phy_lmdz, ONLY: nbp_lon,nbp_lat

#ifdef CPP_XIOS
      USE xios, ONLY: xios_update_calendar
      USE wxios, only: wxios_add_vaxis, wxios_set_cal, wxios_closedef
      USE iophy, ONLY: histwrite_phy
#endif

      IMPLICIT none
!
! Routine argument:
!
      integer,intent(in) :: nlon ! number of atmospheric colums
      integer,intent(in) :: nlev ! number of vertical levels (should be =klev)
      logical,intent(in) :: debut ! signals first call to physics
      logical,intent(in) :: lafin ! signals last call to physics
      real,intent(in) :: pdtphys ! physics time step (s)
      real,intent(in) :: paprs(klon,klev+1) ! interlayer pressure (Pa)
      real,intent(in) :: pplay(klon,klev) ! mid-layer pressure (Pa)
      real,intent(in) :: pphi(klon,klev) ! geopotential at mid-layer
      real,intent(in) :: pphis(klon) ! surface geopotential
      real,intent(in) :: presnivs(klev) ! pseudo-pressure (Pa) of mid-layers
      real,intent(in) :: u(klon,klev) ! eastward zonal wind (m/s)
      real,intent(in) :: v(klon,klev) ! northward meridional wind (m/s)
      real,intent(in) :: t(klon,klev) ! temperature (K)
      real,intent(in) :: qx(klon,klev,nqtot) ! tracers (.../kg_air)
      real,intent(in) :: flxmass_w(klon,klev) ! vertical mass flux
      real,intent(out) :: d_u(klon,klev) ! physics tendency on u (m/s/s)
      real,intent(out) :: d_v(klon,klev) ! physics tendency on v (m/s/s)
      real,intent(out) :: d_t(klon,klev) ! physics tendency on t (K/s)
      real,intent(out) :: d_qx(klon,klev,nqtot) ! physics tendency on tracers
      real,intent(out) :: d_ps(klon) ! physics tendency on surface pressure

integer,save :: itau=0 ! counter to count number of calls to physics
!$OMP THREADPRIVATE(itau)
real :: temp_newton(klon,klev)
integer :: k
logical, save :: first=.true.
!$OMP THREADPRIVATE(first)

! For I/Os
integer :: itau0
real :: zjulian
real :: dtime
integer :: nhori ! horizontal coordinate ID
integer,save :: nid_hist ! output file ID
!$OMP THREADPRIVATE(nid_hist)
integer :: zvertid ! vertical coordinate ID
integer,save :: iwrite_phys ! output every iwrite_phys physics step
!$OMP THREADPRIVATE(iwrite_phys)
integer,save :: iwrite_phys_omp ! intermediate variable to read iwrite_phys
                                ! (must be shared by all threads)
real :: t_ops ! frequency of the IOIPSL operations (eg average over...)
real :: t_wrt ! frequency of the IOIPSL outputs

! initializations
if (debut) then ! Things to do only for the first call to physics 
! load initial conditions for physics (including the grid)
  call phys_state_var_init() ! some initializations, required before calling phyetat0
  call phyetat0("startphy.nc")

! Initialize outputs:
  itau0=0
!$OMP MASTER
  iwrite_phys_omp=1 !default: output every physics timestep
  ! NB: getin() is not threadsafe; only one thread should call it.
  call getin("iwrite_phys",iwrite_phys_omp)
!$OMP END MASTER
!$OMP BARRIER
  iwrite_phys=iwrite_phys_omp
  t_ops=pdtphys*iwrite_phys ! frequency of the IOIPSL operation
  t_wrt=pdtphys*iwrite_phys ! frequency of the outputs in the file
  ! compute zjulian for annee0=1979 and month=1 dayref=1 and hour=0.0
  !CALL ymds2ju(annee0, month, dayref, hour, zjulian)
  call ymds2ju(1979, 1, 1, 0.0, zjulian)
  dtime=pdtphys
#ifndef CPP_IOIPSL_NO_OUTPUT
  ! Initialize IOIPSL output file
  call histbeg_phy("histins.nc",itau0,zjulian,dtime,nhori,nid_hist)
#endif

!$OMP MASTER

#ifndef CPP_IOIPSL_NO_OUTPUT 
! IOIPSL
  ! define vertical coordinate
  call histvert(nid_hist,"presnivs","Vertical levels","Pa",klev, &
                presnivs,zvertid,'down')
  ! define variables which will be written in "histins.nc" file
  call histdef(nid_hist,'temperature','Atmospheric temperature','K', &
               nbp_lon,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'u','Eastward Zonal Wind','m/s', &
               nbp_lon,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'v','Northward Meridional Wind','m/s', &
               nbp_lon,jj_nb,nhori,klev,1,klev,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  call histdef(nid_hist,'ps','Surface Pressure','Pa', &
               nbp_lon,jj_nb,nhori,1,1,1,zvertid,32, &
               'inst(X)',t_ops,t_wrt)
  ! end definition sequence
  call histend(nid_hist)
#endif

#ifdef CPP_XIOS
!XIOS
    ! Declare available vertical axes to be used in output files:    
    CALL wxios_add_vaxis("presnivs", klev, presnivs)

    ! Declare calendar and time step
    CALL wxios_set_cal(dtime,"earth_360d",1,1,1,0.0,1,1,1,0.0)
    
    !Finalize the context:
    CALL wxios_closedef()
#endif
!$OMP END MASTER
!$OMP BARRIER
endif ! of if (debut)

! increment local time counter itau
itau=itau+1

! set all tendencies to zero
d_u(1:klon,1:klev)=0.
d_v(1:klon,1:klev)=0.
d_t(1:klon,1:klev)=0.
d_qx(1:klon,1:klev,1:nqtot)=0.
d_ps(1:klon)=0.

! compute tendencies to return to the dynamics:
! "friction" on the first layer
d_u(1:klon,1)=-u(1:klon,1)/86400.
d_v(1:klon,1)=-v(1:klon,1)/86400.
! newtonian relaxation towards temp_newton()
do k=1,klev
  temp_newton(1:klon,k)=280.+cos(latitude(1:klon))*40.-pphi(1:klon,k)/rg*6.e-3
  d_t(1:klon,k)=(temp_newton(1:klon,k)-t(1:klon,k))/1.e5
enddo


print*,'PHYDEV: itau=',itau

! write some outputs:
! IOIPSL
#ifndef CPP_IOIPSL_NO_OUTPUT 
if (modulo(itau,iwrite_phys)==0) then
  call histwrite_phy(nid_hist,.false.,"temperature",itau,t)
  call histwrite_phy(nid_hist,.false.,"u",itau,u)
  call histwrite_phy(nid_hist,.false.,"v",itau,v)
  call histwrite_phy(nid_hist,.false.,"ps",itau,paprs(:,1))
endif
#endif

!XIOS
#ifdef CPP_XIOS
!$OMP MASTER
    !Increment XIOS time
    CALL xios_update_calendar(itau)
!$OMP END MASTER
!$OMP BARRIER

    !Send fields to XIOS: (NB these fields must also be defined as
    ! <field id="..." /> in iodef.xml to be correctly used
    CALL histwrite_phy("temperature",t)
    CALL histwrite_phy("temp_newton",temp_newton)
    CALL histwrite_phy("u",u)
    CALL histwrite_phy("v",v)
    CALL histwrite_phy("ps",paprs(:,1))
#endif

! if lastcall, then it is time to write "restartphy.nc" file
if (lafin) then
  call phyredem("restartphy.nc")
endif

end subroutine physiq

END MODULE physiq_mod
