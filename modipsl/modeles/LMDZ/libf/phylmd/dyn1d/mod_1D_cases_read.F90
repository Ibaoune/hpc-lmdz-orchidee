!
! $Id: mod_1D_cases_read.F90 2373 2015-10-13 17:28:01Z jyg $
!
MODULE mod_1D_cases_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Declarations specifiques au cas standard
        character*80 :: fich_cas
! Discr?tisation 
        integer nlev_cas, nt_cas


!       integer year_ini_cas, day_ini_cas, mth_ini_cas
!       real heure_ini_cas
!       real day_ju_ini_cas   ! Julian day of case first day
!       parameter (year_ini_cas=2011)
!       parameter (year_ini_cas=1969)
!       parameter (mth_ini_cas=10)
!       parameter (mth_ini_cas=6)
!       parameter (day_ini_cas=1)  ! 10 = 10Juil2006
!       parameter (day_ini_cas=24)  ! 24 = 24 juin 1969
!       parameter (heure_ini_cas=0.) !0h en secondes
!       real pdt_cas
!       parameter (pdt_cas=3.*3600)

!CR ATTENTION TEST AMMA
!        parameter (year_ini_cas=2006)
!        parameter (mth_ini_cas=7)
!        parameter (day_ini_cas=10)  ! 10 = 10Juil2006
!        parameter (heure_ini_cas=0.) !0h en secondes
!        parameter (pdt_cas=1800.)

!profils environnementaux
        real, allocatable::  plev_cas(:,:)

        real, allocatable::  z_cas(:,:)
        real, allocatable::  t_cas(:,:),q_cas(:,:),rh_cas(:,:)
        real, allocatable::  th_cas(:,:),rv_cas(:,:)
        real, allocatable::  u_cas(:,:)
        real, allocatable::  v_cas(:,:)

!forcing
        real, allocatable::  ht_cas(:,:),vt_cas(:,:),dt_cas(:,:),dtrad_cas(:,:)
        real, allocatable::  hth_cas(:,:),vth_cas(:,:),dth_cas(:,:)
        real, allocatable::  hq_cas(:,:),vq_cas(:,:),dq_cas(:,:)
        real, allocatable::  hr_cas(:,:),vr_cas(:,:),dr_cas(:,:)
        real, allocatable::  hu_cas(:,:),vu_cas(:,:),du_cas(:,:)
        real, allocatable::  hv_cas(:,:),vv_cas(:,:),dv_cas(:,:)
        real, allocatable::  vitw_cas(:,:)
        real, allocatable::  ug_cas(:,:),vg_cas(:,:)
        real, allocatable::  lat_cas(:),sens_cas(:),ts_cas(:),ustar_cas(:)
        real, allocatable::  uw_cas(:,:),vw_cas(:,:),q1_cas(:,:),q2_cas(:,:)

!champs interpoles
        real, allocatable::  plev_prof_cas(:)
        real, allocatable::  t_prof_cas(:)
        real, allocatable::  q_prof_cas(:)
        real, allocatable::  u_prof_cas(:)
        real, allocatable::  v_prof_cas(:)        

        real, allocatable::  vitw_prof_cas(:)
        real, allocatable::  ug_prof_cas(:)
        real, allocatable::  vg_prof_cas(:)
        real, allocatable::  ht_prof_cas(:)
        real, allocatable::  hq_prof_cas(:)
        real, allocatable::  vt_prof_cas(:)
        real, allocatable::  vq_prof_cas(:)
        real, allocatable::  dt_prof_cas(:)
        real, allocatable::  dtrad_prof_cas(:)
        real, allocatable::  dq_prof_cas(:)
        real, allocatable::  hu_prof_cas(:)
        real, allocatable::  hv_prof_cas(:)
        real, allocatable::  vu_prof_cas(:)
        real, allocatable::  vv_prof_cas(:)
        real, allocatable::  du_prof_cas(:)
        real, allocatable::  dv_prof_cas(:)
        real, allocatable::  uw_prof_cas(:)
        real, allocatable::  vw_prof_cas(:)
        real, allocatable::  q1_prof_cas(:)
        real, allocatable::  q2_prof_cas(:)


        real lat_prof_cas,sens_prof_cas,ts_prof_cas,ustar_prof_cas
      


CONTAINS

SUBROUTINE read_1D_cas
      implicit none

#include "netcdf.inc"

      INTEGER nid,rid,ierr
      INTEGER ii,jj

      fich_cas='setup/cas.nc'
      print*,'fich_cas ',fich_cas
      ierr = NF_OPEN(fich_cas,NF_NOWRITE,nid)
      print*,'fich_cas,NF_NOWRITE,nid ',fich_cas,NF_NOWRITE,nid
      if (ierr.NE.NF_NOERR) then
         write(*,*) 'ERROR: GROS Pb opening forcings nc file '
         write(*,*) NF_STRERROR(ierr)
         stop ""
      endif
!.......................................................................
      ierr=NF_INQ_DIMID(nid,'lat',rid)
      IF (ierr.NE.NF_NOERR) THEN
         print*, 'Oh probleme lecture dimension lat'
      ENDIF
      ierr=NF_INQ_DIMLEN(nid,rid,ii)
      print*,'OK1 nid,rid,lat',nid,rid,ii
!.......................................................................
      ierr=NF_INQ_DIMID(nid,'lon',rid)
      IF (ierr.NE.NF_NOERR) THEN
         print*, 'Oh probleme lecture dimension lon'
      ENDIF
      ierr=NF_INQ_DIMLEN(nid,rid,jj)
      print*,'OK2 nid,rid,lat',nid,rid,jj
!.......................................................................
      ierr=NF_INQ_DIMID(nid,'lev',rid)
      IF (ierr.NE.NF_NOERR) THEN
         print*, 'Oh probleme lecture dimension zz'
      ENDIF
      ierr=NF_INQ_DIMLEN(nid,rid,nlev_cas)
      print*,'OK3 nid,rid,nlev_cas',nid,rid,nlev_cas
!.......................................................................
      ierr=NF_INQ_DIMID(nid,'time',rid)
      print*,'nid,rid',nid,rid
      nt_cas=0
      IF (ierr.NE.NF_NOERR) THEN
        stop 'probleme lecture dimension sens'
      ENDIF
      ierr=NF_INQ_DIMLEN(nid,rid,nt_cas)
      print*,'OK4 nid,rid,nt_cas',nid,rid,nt_cas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!profils moyens:
        allocate(plev_cas(nlev_cas,nt_cas))        
        allocate(z_cas(nlev_cas,nt_cas))
        allocate(t_cas(nlev_cas,nt_cas),q_cas(nlev_cas,nt_cas),rh_cas(nlev_cas,nt_cas))
        allocate(th_cas(nlev_cas,nt_cas),rv_cas(nlev_cas,nt_cas))
        allocate(u_cas(nlev_cas,nt_cas))
        allocate(v_cas(nlev_cas,nt_cas))

!forcing
        allocate(ht_cas(nlev_cas,nt_cas),vt_cas(nlev_cas,nt_cas),dt_cas(nlev_cas,nt_cas),dtrad_cas(nlev_cas,nt_cas))
        allocate(hq_cas(nlev_cas,nt_cas),vq_cas(nlev_cas,nt_cas),dq_cas(nlev_cas,nt_cas))
        allocate(hth_cas(nlev_cas,nt_cas),vth_cas(nlev_cas,nt_cas),dth_cas(nlev_cas,nt_cas))
        allocate(hr_cas(nlev_cas,nt_cas),vr_cas(nlev_cas,nt_cas),dr_cas(nlev_cas,nt_cas))
        allocate(hu_cas(nlev_cas,nt_cas),vu_cas(nlev_cas,nt_cas),du_cas(nlev_cas,nt_cas))
        allocate(hv_cas(nlev_cas,nt_cas),vv_cas(nlev_cas,nt_cas),dv_cas(nlev_cas,nt_cas))
        allocate(vitw_cas(nlev_cas,nt_cas))
        allocate(ug_cas(nlev_cas,nt_cas))
        allocate(vg_cas(nlev_cas,nt_cas))
        allocate(lat_cas(nt_cas),sens_cas(nt_cas),ts_cas(nt_cas),ustar_cas(nt_cas))
        allocate(uw_cas(nlev_cas,nt_cas),vw_cas(nlev_cas,nt_cas),q1_cas(nlev_cas,nt_cas),q2_cas(nlev_cas,nt_cas))


!champs interpoles
        allocate(plev_prof_cas(nlev_cas))
        allocate(t_prof_cas(nlev_cas))
        allocate(q_prof_cas(nlev_cas))
        allocate(u_prof_cas(nlev_cas))
        allocate(v_prof_cas(nlev_cas))

        allocate(vitw_prof_cas(nlev_cas))
        allocate(ug_prof_cas(nlev_cas))
        allocate(vg_prof_cas(nlev_cas))
        allocate(ht_prof_cas(nlev_cas))
        allocate(hq_prof_cas(nlev_cas))
        allocate(hu_prof_cas(nlev_cas))
        allocate(hv_prof_cas(nlev_cas))
        allocate(vt_prof_cas(nlev_cas))
        allocate(vq_prof_cas(nlev_cas))
        allocate(vu_prof_cas(nlev_cas))
        allocate(vv_prof_cas(nlev_cas))
        allocate(dt_prof_cas(nlev_cas))
        allocate(dtrad_prof_cas(nlev_cas))
        allocate(dq_prof_cas(nlev_cas))
        allocate(du_prof_cas(nlev_cas))
        allocate(dv_prof_cas(nlev_cas))
        allocate(uw_prof_cas(nlev_cas))
        allocate(vw_prof_cas(nlev_cas))
        allocate(q1_prof_cas(nlev_cas))
        allocate(q2_prof_cas(nlev_cas))

        print*,'Allocations OK'
        call read_cas(nid,nlev_cas,nt_cas                                       &
     &     ,z_cas,plev_cas,t_cas,q_cas,rh_cas,th_cas,rv_cas,u_cas,v_cas         &
     &     ,ug_cas,vg_cas,vitw_cas,du_cas,hu_cas,vu_cas,dv_cas,hv_cas,vv_cas    &
     &     ,dt_cas,dtrad_cas,ht_cas,vt_cas,dq_cas,hq_cas,vq_cas                 &
     &     ,dth_cas,hth_cas,vth_cas,dr_cas,hr_cas,vr_cas,sens_cas,lat_cas,ts_cas&
     &     ,ustar_cas,uw_cas,vw_cas,q1_cas,q2_cas)
        print*,'Read cas OK'


END SUBROUTINE read_1D_cas



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE deallocate_1D_cases
!profils environnementaux:
        deallocate(plev_cas)
        
        deallocate(z_cas)
        deallocate(t_cas,q_cas,rh_cas)
        deallocate(th_cas,rv_cas)
        deallocate(u_cas)
        deallocate(v_cas)
        
!forcing
        deallocate(ht_cas,vt_cas,dt_cas,dtrad_cas)
        deallocate(hq_cas,vq_cas,dq_cas)
        deallocate(hth_cas,vth_cas,dth_cas)
        deallocate(hr_cas,vr_cas,dr_cas)
        deallocate(hu_cas,vu_cas,du_cas)
        deallocate(hv_cas,vv_cas,dv_cas)
        deallocate(vitw_cas)
        deallocate(ug_cas)
        deallocate(vg_cas)
        deallocate(lat_cas,sens_cas,ts_cas,ustar_cas,uw_cas,vw_cas,q1_cas,q2_cas)

!champs interpoles
        deallocate(plev_prof_cas)
        deallocate(t_prof_cas)
        deallocate(q_prof_cas)
        deallocate(u_prof_cas)
        deallocate(v_prof_cas)

        deallocate(vitw_prof_cas)
        deallocate(ug_prof_cas)
        deallocate(vg_prof_cas)
        deallocate(ht_prof_cas)
        deallocate(hq_prof_cas)
        deallocate(hu_prof_cas)
        deallocate(hv_prof_cas)
        deallocate(vt_prof_cas)
        deallocate(vq_prof_cas)
        deallocate(vu_prof_cas)
        deallocate(vv_prof_cas)
        deallocate(dt_prof_cas)
        deallocate(dtrad_prof_cas)
        deallocate(dq_prof_cas)
        deallocate(du_prof_cas)
        deallocate(dv_prof_cas)
        deallocate(t_prof_cas)
        deallocate(q_prof_cas)
        deallocate(u_prof_cas)
        deallocate(v_prof_cas)
        deallocate(uw_prof_cas)
        deallocate(vw_prof_cas)
        deallocate(q1_prof_cas)
        deallocate(q2_prof_cas)

END SUBROUTINE deallocate_1D_cases


END MODULE mod_1D_cases_read
!=====================================================================
      subroutine read_cas(nid,nlevel,ntime                          &
     &     ,zz,pp,temp,qv,rh,theta,rv,u,v,ug,vg,w,                   &
     &     du,hu,vu,dv,hv,vv,dt,dtrad,ht,vt,dq,hq,vq,                     &
     &     dth,hth,vth,dr,hr,vr,sens,flat,ts,ustar,uw,vw,q1,q2)

!program reading forcing of the case study
      implicit none
#include "netcdf.inc"

      integer ntime,nlevel

      real zz(nlevel,ntime)
      real pp(nlevel,ntime)
      real temp(nlevel,ntime),qv(nlevel,ntime),rh(nlevel,ntime)
      real theta(nlevel,ntime),rv(nlevel,ntime)
      real u(nlevel,ntime)
      real v(nlevel,ntime)
      real ug(nlevel,ntime)
      real vg(nlevel,ntime)
      real w(nlevel,ntime)
      real du(nlevel,ntime),hu(nlevel,ntime),vu(nlevel,ntime)
      real dv(nlevel,ntime),hv(nlevel,ntime),vv(nlevel,ntime)
      real dt(nlevel,ntime),ht(nlevel,ntime),vt(nlevel,ntime)
      real dtrad(nlevel,ntime)
      real dq(nlevel,ntime),hq(nlevel,ntime),vq(nlevel,ntime)
      real dth(nlevel,ntime),hth(nlevel,ntime),vth(nlevel,ntime)
      real dr(nlevel,ntime),hr(nlevel,ntime),vr(nlevel,ntime)
      real flat(ntime),sens(ntime),ts(ntime),ustar(ntime)
      real uw(nlevel,ntime),vw(nlevel,ntime),q1(nlevel,ntime),q2(nlevel,ntime)


      integer nid, ierr,rid
      integer nbvar3d
      parameter(nbvar3d=39)
      integer var3didin(nbvar3d)

       ierr=NF_INQ_VARID(nid,"zz",var3didin(1)) 
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'lev'
         endif
      
      ierr=NF_INQ_VARID(nid,"pp",var3didin(2)) 
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'plev'
         endif


      ierr=NF_INQ_VARID(nid,"temp",var3didin(3))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'temp'
         endif

      ierr=NF_INQ_VARID(nid,"qv",var3didin(4))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'qv'
         endif

      ierr=NF_INQ_VARID(nid,"rh",var3didin(5))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'rh'
         endif

      ierr=NF_INQ_VARID(nid,"theta",var3didin(6))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'theta'
         endif

      ierr=NF_INQ_VARID(nid,"rv",var3didin(7))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'rv'
         endif


      ierr=NF_INQ_VARID(nid,"u",var3didin(8))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'u'
         endif

      ierr=NF_INQ_VARID(nid,"v",var3didin(9))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'v'
         endif

       ierr=NF_INQ_VARID(nid,"ug",var3didin(10))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'ug'
         endif

      ierr=NF_INQ_VARID(nid,"vg",var3didin(11))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'vg'
         endif

      ierr=NF_INQ_VARID(nid,"w",var3didin(12))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'w'
         endif

      ierr=NF_INQ_VARID(nid,"advu",var3didin(13))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'advu'
         endif

      ierr=NF_INQ_VARID(nid,"hu",var3didin(14))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hu'
         endif

       ierr=NF_INQ_VARID(nid,"vu",var3didin(15))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'vu'
         endif

       ierr=NF_INQ_VARID(nid,"advv",var3didin(16))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'advv'
         endif

      ierr=NF_INQ_VARID(nid,"hv",var3didin(17))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hv'
         endif

       ierr=NF_INQ_VARID(nid,"vv",var3didin(18))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'vv'
         endif

      ierr=NF_INQ_VARID(nid,"advT",var3didin(19))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'advT'
         endif

      ierr=NF_INQ_VARID(nid,"hT",var3didin(20))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hT'
         endif

      ierr=NF_INQ_VARID(nid,"vT",var3didin(21))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'vT'
         endif

      ierr=NF_INQ_VARID(nid,"advq",var3didin(22))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'advq'
         endif
      
      ierr=NF_INQ_VARID(nid,"hq",var3didin(23))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hq'
         endif

      ierr=NF_INQ_VARID(nid,"vq",var3didin(24))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'vq'
         endif

      ierr=NF_INQ_VARID(nid,"advth",var3didin(25))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'advth'
         endif

      ierr=NF_INQ_VARID(nid,"hth",var3didin(26))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hth'
         endif

      ierr=NF_INQ_VARID(nid,"vth",var3didin(27))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'vth'
         endif

      ierr=NF_INQ_VARID(nid,"advr",var3didin(28))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'advr'
         endif
      
      ierr=NF_INQ_VARID(nid,"hr",var3didin(29))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hr'
         endif

      ierr=NF_INQ_VARID(nid,"vr",var3didin(30))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'vr'
         endif

      ierr=NF_INQ_VARID(nid,"radT",var3didin(31))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'radT'
         endif

      ierr=NF_INQ_VARID(nid,"sens",var3didin(32))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'sens'
         endif

      ierr=NF_INQ_VARID(nid,"flat",var3didin(33))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'flat'
         endif

      ierr=NF_INQ_VARID(nid,"ts",var3didin(34))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'ts'
         endif

      ierr=NF_INQ_VARID(nid,"ustar",var3didin(35))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'ustar'
         endif

      ierr=NF_INQ_VARID(nid,"uw",var3didin(36))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'uw'
         endif

      ierr=NF_INQ_VARID(nid,"vw",var3didin(37))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'vw'
         endif

      ierr=NF_INQ_VARID(nid,"q1",var3didin(38))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'q1'
         endif

      ierr=NF_INQ_VARID(nid,"q2",var3didin(39))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'q2'
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(1),zz)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(1),zz)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture z ok',zz

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(2),pp)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(2),pp)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture pp ok',pp


#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(3),temp)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(3),temp)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture T ok',temp

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(4),qv)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(4),qv)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture qv ok',qv
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(5),rh)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(5),rh)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture rh ok',rh

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(6),theta)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(6),theta)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture theta ok',theta

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(7),rv)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(7),rv)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture rv ok',rv

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(8),u)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(8),u)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture u ok',u

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(9),v)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(9),v)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture v ok',v

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(10),ug)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(10),ug)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture ug ok',ug

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(11),vg)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(11),vg)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture vg ok',vg

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(12),w)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(12),w)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture w ok',w

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(13),du)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(13),du)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture du ok',du

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(14),hu)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(14),hu)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture hu ok',hu

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(15),vu)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(15),vu)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture vu ok',vu

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(16),dv)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(16),dv)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dv ok',dv

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(17),hv)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(17),hv)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture hv ok',hv

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(18),vv)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(18),vv)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture vv ok',vv

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(19),dt)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(19),dt)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dt ok',dt

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(20),ht)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(20),ht)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture ht ok',ht

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(21),vt)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(21),vt)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture vt ok',vt

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(22),dq)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(22),dq)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dq ok',dq

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(23),hq)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(23),hq)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture hq ok',hq

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(24),vq)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(24),vq)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture vq ok',vq

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(25),dth)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(25),dth)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dth ok',dth

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(26),hth)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(26),hth)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture hth ok',hth

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(27),vth)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(27),vth)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture vth ok',vth

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(28),dr)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(28),dr)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dr ok',dr

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(29),hr)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(29),hr)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture hr ok',hr

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(30),vr)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(30),vr)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture vr ok',vr

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(31),dtrad)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(31),dtrad)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dtrad ok',dtrad

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(32),sens)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(32),sens)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture sens ok',sens

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(33),flat)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(33),flat)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture flat ok',flat

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(34),ts)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(34),ts)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture ts ok',ts

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(35),ustar)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(35),ustar)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture ustar ok',ustar

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(36),uw)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(36),uw)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture uw ok',uw

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(37),vw)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(37),vw)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture vw ok',vw

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(38),q1)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(38),q1)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture q1 ok',q1

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(39),q2)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(39),q2)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture q2 ok',q2


         return 
         end subroutine read_cas
!======================================================================
        SUBROUTINE interp_case_time(day,day1,annee_ref                &
!    &         ,year_cas,day_cas,nt_cas,pdt_forc,nlev_cas      &
     &         ,nt_cas,nlev_cas                                       &
     &         ,ts_cas,plev_cas,t_cas,q_cas,u_cas,v_cas               &
     &         ,ug_cas,vg_cas,vitw_cas,du_cas,hu_cas,vu_cas           &
     &         ,dv_cas,hv_cas,vv_cas,dt_cas,ht_cas,vt_cas,dtrad_cas   &
     &         ,dq_cas,hq_cas,vq_cas,lat_cas,sens_cas,ustar_cas       &
     &         ,uw_cas,vw_cas,q1_cas,q2_cas                           &
     &         ,ts_prof_cas,plev_prof_cas,t_prof_cas,q_prof_cas       &
     &         ,u_prof_cas,v_prof_cas,ug_prof_cas,vg_prof_cas         &
     &         ,vitw_prof_cas,du_prof_cas,hu_prof_cas,vu_prof_cas     &
     &         ,dv_prof_cas,hv_prof_cas,vv_prof_cas,dt_prof_cas       &
     &         ,ht_prof_cas,vt_prof_cas,dtrad_prof_cas,dq_prof_cas    &
     &         ,hq_prof_cas,vq_prof_cas,lat_prof_cas,sens_prof_cas    &
     &         ,ustar_prof_cas,uw_prof_cas,vw_prof_cas,q1_prof_cas,q2_prof_cas)
          

        implicit none

!---------------------------------------------------------------------------------------
! Time interpolation of a 2D field to the timestep corresponding to day
!
! day: current julian day (e.g. 717538.2)
! day1: first day of the simulation
! nt_cas: total nb of data in the forcing 
! pdt_cas: total time interval (in sec) between 2 forcing data
!---------------------------------------------------------------------------------------

#include "compar1d.h"
#include "date_cas.h"

! inputs:
        integer annee_ref
        integer nt_cas,nlev_cas
        real day, day1,day_cas
        real ts_cas(nt_cas)
        real plev_cas(nlev_cas,nt_cas)
        real t_cas(nlev_cas,nt_cas),q_cas(nlev_cas,nt_cas)
        real u_cas(nlev_cas,nt_cas),v_cas(nlev_cas,nt_cas)
        real ug_cas(nlev_cas,nt_cas),vg_cas(nlev_cas,nt_cas)
        real vitw_cas(nlev_cas,nt_cas)
        real du_cas(nlev_cas,nt_cas),hu_cas(nlev_cas,nt_cas),vu_cas(nlev_cas,nt_cas)
        real dv_cas(nlev_cas,nt_cas),hv_cas(nlev_cas,nt_cas),vv_cas(nlev_cas,nt_cas)
        real dt_cas(nlev_cas,nt_cas),ht_cas(nlev_cas,nt_cas),vt_cas(nlev_cas,nt_cas)
        real dtrad_cas(nlev_cas,nt_cas)
        real dq_cas(nlev_cas,nt_cas),hq_cas(nlev_cas,nt_cas),vq_cas(nlev_cas,nt_cas)
        real lat_cas(nt_cas)
        real sens_cas(nt_cas)
        real ustar_cas(nt_cas),uw_cas(nlev_cas,nt_cas),vw_cas(nlev_cas,nt_cas)
        real q1_cas(nlev_cas,nt_cas),q2_cas(nlev_cas,nt_cas)

! outputs:
        real plev_prof_cas(nlev_cas)
        real t_prof_cas(nlev_cas),q_prof_cas(nlev_cas)
        real u_prof_cas(nlev_cas),v_prof_cas(nlev_cas)
        real ug_prof_cas(nlev_cas),vg_prof_cas(nlev_cas)
        real vitw_prof_cas(nlev_cas)
        real du_prof_cas(nlev_cas),hu_prof_cas(nlev_cas),vu_prof_cas(nlev_cas)
        real dv_prof_cas(nlev_cas),hv_prof_cas(nlev_cas),vv_prof_cas(nlev_cas)
        real dt_prof_cas(nlev_cas),ht_prof_cas(nlev_cas),vt_prof_cas(nlev_cas)
        real dtrad_prof_cas(nlev_cas)
        real dq_prof_cas(nlev_cas),hq_prof_cas(nlev_cas),vq_prof_cas(nlev_cas)
        real lat_prof_cas,sens_prof_cas,ts_prof_cas,ustar_prof_cas
        real uw_prof_cas(nlev_cas),vw_prof_cas(nlev_cas),q1_prof_cas(nlev_cas),q2_prof_cas(nlev_cas)
! local:
        integer it_cas1, it_cas2,k
        real timeit,time_cas1,time_cas2,frac


        print*,'Check time',day1,day_ju_ini_cas,day_deb+1,pdt_cas

! On teste si la date du cas AMMA est correcte.
! C est pour memoire car en fait les fichiers .def
! sont censes etre corrects.
! A supprimer a terme (MPL 20150623)
!     if ((forcing_type.eq.10).and.(1.eq.0)) then
! Check that initial day of the simulation consistent with AMMA case:
!      if (annee_ref.ne.2006) then
!       print*,'Pour AMMA, annee_ref doit etre 2006'
!       print*,'Changer annee_ref dans run.def'
!       stop
!      endif
!      if (annee_ref.eq.2006 .and. day1.lt.day_cas) then
!       print*,'AMMA a debute le 10 juillet 2006',day1,day_cas
!       print*,'Changer dayref dans run.def'
!       stop
!      endif
!      if (annee_ref.eq.2006 .and. day1.gt.day_cas+1) then
!       print*,'AMMA a fini le 11 juillet'
!       print*,'Changer dayref ou nday dans run.def'
!       stop
!      endif
!      endif

! Determine timestep relative to the 1st day:
!       timeit=(day-day1)*86400.
!       if (annee_ref.eq.1992) then
!        timeit=(day-day_cas)*86400.
!       else
!        timeit=(day+61.-1.)*86400. ! 61 days between Nov01 and Dec31 1992
!       endif
      timeit=(day-day_ju_ini_cas)*86400
      print *,'day=',day
      print *,'day_ju_ini_cas=',day_ju_ini_cas
      print *,'pdt_cas=',pdt_cas
      print *,'timeit=',timeit
      print *,'nt_cas=',nt_cas

! Determine the closest observation times:
!       it_cas1=INT(timeit/pdt_cas)+1
!       it_cas2=it_cas1 + 1
!       time_cas1=(it_cas1-1)*pdt_cas
!       time_cas2=(it_cas2-1)*pdt_cas

       it_cas1=INT(timeit/pdt_cas)+1
       IF (it_cas1 .EQ. nt_cas) THEN
       it_cas2=it_cas1 
       ELSE
       it_cas2=it_cas1 + 1
       ENDIF
       time_cas1=(it_cas1-1)*pdt_cas
       time_cas2=(it_cas2-1)*pdt_cas
      print *,'it_cas1=',it_cas1
      print *,'it_cas2=',it_cas2
      print *,'time_cas1=',time_cas1
      print *,'time_cas2=',time_cas2

       if (it_cas1 .gt. nt_cas) then
        write(*,*) 'PB-stop: day, day_ju_ini_cas,it_cas1, it_cas2, timeit: '            &
     &        ,day,day_ju_ini_cas,it_cas1,it_cas2,timeit
        stop
       endif

! time interpolation:
       IF (it_cas1 .EQ. it_cas2) THEN
          frac=0.
       ELSE
          frac=(time_cas2-timeit)/(time_cas2-time_cas1)
          frac=max(frac,0.0)
       ENDIF

       lat_prof_cas = lat_cas(it_cas2)                                       &
     &          -frac*(lat_cas(it_cas2)-lat_cas(it_cas1)) 
       sens_prof_cas = sens_cas(it_cas2)                                     &
     &          -frac*(sens_cas(it_cas2)-sens_cas(it_cas1))
       ts_prof_cas = ts_cas(it_cas2)                                         &
     &          -frac*(ts_cas(it_cas2)-ts_cas(it_cas1))
       ustar_prof_cas = ustar_cas(it_cas2)                                   &
     &          -frac*(ustar_cas(it_cas2)-ustar_cas(it_cas1))

       do k=1,nlev_cas
        plev_prof_cas(k) = plev_cas(k,it_cas2)                               &
     &          -frac*(plev_cas(k,it_cas2)-plev_cas(k,it_cas1))
        t_prof_cas(k) = t_cas(k,it_cas2)                               &
     &          -frac*(t_cas(k,it_cas2)-t_cas(k,it_cas1))
        q_prof_cas(k) = q_cas(k,it_cas2)                               &
     &          -frac*(q_cas(k,it_cas2)-q_cas(k,it_cas1))
        u_prof_cas(k) = u_cas(k,it_cas2)                               &
     &          -frac*(u_cas(k,it_cas2)-u_cas(k,it_cas1))
        v_prof_cas(k) = v_cas(k,it_cas2)                               &
     &          -frac*(v_cas(k,it_cas2)-v_cas(k,it_cas1))
        ug_prof_cas(k) = ug_cas(k,it_cas2)                               &
     &          -frac*(ug_cas(k,it_cas2)-ug_cas(k,it_cas1))
        vg_prof_cas(k) = vg_cas(k,it_cas2)                               &
     &          -frac*(vg_cas(k,it_cas2)-vg_cas(k,it_cas1))
        vitw_prof_cas(k) = vitw_cas(k,it_cas2)                               &
     &          -frac*(vitw_cas(k,it_cas2)-vitw_cas(k,it_cas1))
        du_prof_cas(k) = du_cas(k,it_cas2)                                   &
     &          -frac*(du_cas(k,it_cas2)-du_cas(k,it_cas1))
        hu_prof_cas(k) = hu_cas(k,it_cas2)                                   &
     &          -frac*(hu_cas(k,it_cas2)-hu_cas(k,it_cas1))
        vu_prof_cas(k) = vu_cas(k,it_cas2)                                   &
     &          -frac*(vu_cas(k,it_cas2)-vu_cas(k,it_cas1))
        dv_prof_cas(k) = dv_cas(k,it_cas2)                                   &
     &          -frac*(dv_cas(k,it_cas2)-dv_cas(k,it_cas1))
        hv_prof_cas(k) = hv_cas(k,it_cas2)                                   &
     &          -frac*(hv_cas(k,it_cas2)-hv_cas(k,it_cas1))
        vv_prof_cas(k) = vv_cas(k,it_cas2)                                   &
     &          -frac*(vv_cas(k,it_cas2)-vv_cas(k,it_cas1))
        dt_prof_cas(k) = dt_cas(k,it_cas2)                                   &
     &          -frac*(dt_cas(k,it_cas2)-dt_cas(k,it_cas1))
        ht_prof_cas(k) = ht_cas(k,it_cas2)                                   &
     &          -frac*(ht_cas(k,it_cas2)-ht_cas(k,it_cas1))
        vt_prof_cas(k) = vt_cas(k,it_cas2)                                   &
     &          -frac*(vt_cas(k,it_cas2)-vt_cas(k,it_cas1))
        dtrad_prof_cas(k) = dtrad_cas(k,it_cas2)                                   &
     &          -frac*(dtrad_cas(k,it_cas2)-dtrad_cas(k,it_cas1))
        dq_prof_cas(k) = dq_cas(k,it_cas2)                                   &
     &          -frac*(dq_cas(k,it_cas2)-dq_cas(k,it_cas1))
        hq_prof_cas(k) = hq_cas(k,it_cas2)                                   &
     &          -frac*(hq_cas(k,it_cas2)-hq_cas(k,it_cas1))
        vq_prof_cas(k) = vq_cas(k,it_cas2)                                   &
     &          -frac*(vq_cas(k,it_cas2)-vq_cas(k,it_cas1))
       uw_prof_cas(k) = uw_cas(k,it_cas2)                                   &
     &          -frac*(uw_cas(k,it_cas2)-uw_cas(k,it_cas1))
       vw_prof_cas(k) = vw_cas(k,it_cas2)                                   &
     &          -frac*(vw_cas(k,it_cas2)-vw_cas(k,it_cas1))
       q1_prof_cas(k) = q1_cas(k,it_cas2)                                   &
     &          -frac*(q1_cas(k,it_cas2)-q1_cas(k,it_cas1))
       q2_prof_cas(k) = q2_cas(k,it_cas2)                                   &
     &          -frac*(q2_cas(k,it_cas2)-q2_cas(k,it_cas1))
        enddo

        return
        END

!**********************************************************************************************
