!
! $Id: mod_1D_cases_read.F90 2373 2015-10-13 17:28:01Z jyg $
!
MODULE mod_1D_cases_read_std

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Declarations specifiques au cas standard
        character*80 :: fich_cas
! Discr?tisation 
        integer nlev_cas, nt_cas


!profils environnementaux
        real, allocatable::  plev_cas(:,:),plevh_cas(:)
        real, allocatable::  ap_cas(:),bp_cas(:)

        real, allocatable::  z_cas(:,:),zh_cas(:)
        real, allocatable::  t_cas(:,:),q_cas(:,:),qv_cas(:,:),ql_cas(:,:),qi_cas(:,:),rh_cas(:,:)
        real, allocatable::  th_cas(:,:),thv_cas(:,:),thl_cas(:,:),rv_cas(:,:)
        real, allocatable::  u_cas(:,:),v_cas(:,:),vitw_cas(:,:),omega_cas(:,:)

!forcing
        real, allocatable::  ht_cas(:,:),vt_cas(:,:),dt_cas(:,:),dtrad_cas(:,:)
        real, allocatable::  hth_cas(:,:),vth_cas(:,:),dth_cas(:,:)
        real, allocatable::  hq_cas(:,:),vq_cas(:,:),dq_cas(:,:)
        real, allocatable::  hr_cas(:,:),vr_cas(:,:),dr_cas(:,:)
        real, allocatable::  hu_cas(:,:),vu_cas(:,:),du_cas(:,:)
        real, allocatable::  hv_cas(:,:),vv_cas(:,:),dv_cas(:,:)
        real, allocatable::  ug_cas(:,:),vg_cas(:,:)
        real, allocatable::  temp_nudg_cas(:,:),qv_nudg_cas(:,:),u_nudg_cas(:,:),v_nudg_cas(:,:)
        real, allocatable::  lat_cas(:),sens_cas(:),ts_cas(:),ps_cas(:),ustar_cas(:)
        real, allocatable::  uw_cas(:,:),vw_cas(:,:),q1_cas(:,:),q2_cas(:,:),tke_cas(:)

!champs interpoles
        real, allocatable::  plev_prof_cas(:)
        real, allocatable::  t_prof_cas(:)
        real, allocatable::  theta_prof_cas(:)
        real, allocatable::  thl_prof_cas(:)
        real, allocatable::  thv_prof_cas(:)
        real, allocatable::  q_prof_cas(:)
        real, allocatable::  qv_prof_cas(:)
        real, allocatable::  ql_prof_cas(:)
        real, allocatable::  qi_prof_cas(:)
        real, allocatable::  rh_prof_cas(:)
        real, allocatable::  rv_prof_cas(:)
        real, allocatable::  u_prof_cas(:)
        real, allocatable::  v_prof_cas(:)        
        real, allocatable::  vitw_prof_cas(:)
        real, allocatable::  omega_prof_cas(:)
        real, allocatable::  ug_prof_cas(:)
        real, allocatable::  vg_prof_cas(:)
        real, allocatable::  temp_nudg_prof_cas(:),qv_nudg_prof_cas(:),u_nudg_prof_cas(:),v_nudg_prof_cas(:)
        real, allocatable::  ht_prof_cas(:)
        real, allocatable::  hth_prof_cas(:)
        real, allocatable::  hq_prof_cas(:)
        real, allocatable::  vt_prof_cas(:)
        real, allocatable::  vth_prof_cas(:)
        real, allocatable::  vq_prof_cas(:)
        real, allocatable::  dt_prof_cas(:)
        real, allocatable::  dth_prof_cas(:)
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


        real lat_prof_cas,sens_prof_cas,ts_prof_cas,ps_prof_cas,ustar_prof_cas,tke_prof_cas
        real o3_cas,orog_cas,albedo_cas,emiss_cas,t_skin_cas,q_skin_cas,mom_rough,heat_rough,rugos_cas,sand_cas,clay_cas
      


CONTAINS


!**********************************************************************************************
SUBROUTINE read_SCM_cas
      implicit none

#include "netcdf.inc"
#include "date_cas.h"

      INTEGER nid,rid,ierr
      INTEGER ii,jj,timeid
      REAL, ALLOCATABLE :: time_val(:)

      print*,'ON EST VRAIMENT LA'
      fich_cas='cas.nc'
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
      print*,'OK1 read2: nid,rid,lat',nid,rid,ii
!.......................................................................
      ierr=NF_INQ_DIMID(nid,'lon',rid)
      IF (ierr.NE.NF_NOERR) THEN
         print*, 'Oh probleme lecture dimension lon'
      ENDIF
      ierr=NF_INQ_DIMLEN(nid,rid,jj)
      print*,'OK2 read2: nid,rid,lat',nid,rid,jj
!.......................................................................
      ierr=NF_INQ_DIMID(nid,'lev',rid)
      IF (ierr.NE.NF_NOERR) THEN
         print*, 'Oh probleme lecture dimension nlev'
      ENDIF
      ierr=NF_INQ_DIMLEN(nid,rid,nlev_cas)
      print*,'OK3 read2: nid,rid,nlev_cas',nid,rid,nlev_cas
      IF ( .NOT. ( nlev_cas > 10 .AND. nlev_cas < 200000 )) THEN
              print*,'Valeur de nlev_cas peu probable'
              STOP
      ENDIF
!.......................................................................
      ierr=NF_INQ_DIMID(nid,'time',rid)
      nt_cas=0
      IF (ierr.NE.NF_NOERR) THEN
        stop 'Oh probleme lecture dimension time'
      ENDIF
      ierr=NF_INQ_DIMLEN(nid,rid,nt_cas)
      print*,'OK4 read2: nid,rid,nt_cas',nid,rid,nt_cas
! Lecture de l'axe des temps
      print*,'LECTURE DU TEMPS'
      ierr=NF_INQ_VARID(nid,'time',timeid)
         if(ierr/=NF_NOERR) then
           print *,'Variable time manquante dans cas.nc:'
           ierr=NF_NOERR
         else
                 allocate(time_val(nt_cas))
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,timeid,time_val)
#else
           ierr = NF_GET_VAR_REAL(nid,timeid,time_val)
#endif
           if(ierr/=NF_NOERR) then
              print *,'Pb a la lecture de time cas.nc: '
           endif
   endif
   IF (nt_cas>1) THEN
           pdt_cas=time_val(2)-time_val(1)
   ELSE
           pdt_cas=0.
   ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!profils moyens:
        allocate(plev_cas(nlev_cas,nt_cas),plevh_cas(nlev_cas+1))        
        allocate(z_cas(nlev_cas,nt_cas),zh_cas(nlev_cas+1))
        allocate(ap_cas(nlev_cas+1),bp_cas(nt_cas+1))
        allocate(t_cas(nlev_cas,nt_cas),q_cas(nlev_cas,nt_cas),qv_cas(nlev_cas,nt_cas),ql_cas(nlev_cas,nt_cas), &
             qi_cas(nlev_cas,nt_cas),rh_cas(nlev_cas,nt_cas))
        allocate(th_cas(nlev_cas,nt_cas),thl_cas(nlev_cas,nt_cas),thv_cas(nlev_cas,nt_cas),rv_cas(nlev_cas,nt_cas))
        allocate(u_cas(nlev_cas,nt_cas),v_cas(nlev_cas,nt_cas),vitw_cas(nlev_cas,nt_cas),omega_cas(nlev_cas,nt_cas))

!forcing
        allocate(ht_cas(nlev_cas,nt_cas),vt_cas(nlev_cas,nt_cas),dt_cas(nlev_cas,nt_cas),dtrad_cas(nlev_cas,nt_cas))
        allocate(hq_cas(nlev_cas,nt_cas),vq_cas(nlev_cas,nt_cas),dq_cas(nlev_cas,nt_cas))
        allocate(hth_cas(nlev_cas,nt_cas),vth_cas(nlev_cas,nt_cas),dth_cas(nlev_cas,nt_cas))
        allocate(hr_cas(nlev_cas,nt_cas),vr_cas(nlev_cas,nt_cas),dr_cas(nlev_cas,nt_cas))
        allocate(hu_cas(nlev_cas,nt_cas),vu_cas(nlev_cas,nt_cas),du_cas(nlev_cas,nt_cas))
        allocate(hv_cas(nlev_cas,nt_cas),vv_cas(nlev_cas,nt_cas),dv_cas(nlev_cas,nt_cas))
        allocate(ug_cas(nlev_cas,nt_cas),vg_cas(nlev_cas,nt_cas))
        allocate(temp_nudg_cas(nlev_cas,nt_cas),qv_nudg_cas(nlev_cas,nt_cas))
        allocate(u_nudg_cas(nlev_cas,nt_cas),v_nudg_cas(nlev_cas,nt_cas))
        allocate(lat_cas(nt_cas),sens_cas(nt_cas),ts_cas(nt_cas),ps_cas(nt_cas),ustar_cas(nt_cas),tke_cas(nt_cas))
        allocate(uw_cas(nlev_cas,nt_cas),vw_cas(nlev_cas,nt_cas),q1_cas(nlev_cas,nt_cas),q2_cas(nlev_cas,nt_cas))



!champs interpoles
        allocate(plev_prof_cas(nlev_cas))
        allocate(t_prof_cas(nlev_cas))
        allocate(theta_prof_cas(nlev_cas))
        allocate(thl_prof_cas(nlev_cas))
        allocate(thv_prof_cas(nlev_cas))
        allocate(q_prof_cas(nlev_cas))
        allocate(qv_prof_cas(nlev_cas))
        allocate(ql_prof_cas(nlev_cas))
        allocate(qi_prof_cas(nlev_cas))
        allocate(rh_prof_cas(nlev_cas))
        allocate(rv_prof_cas(nlev_cas))
        allocate(u_prof_cas(nlev_cas))
        allocate(v_prof_cas(nlev_cas))
        allocate(vitw_prof_cas(nlev_cas))
        allocate(omega_prof_cas(nlev_cas))
        allocate(ug_prof_cas(nlev_cas))
        allocate(vg_prof_cas(nlev_cas))
        allocate(temp_nudg_prof_cas(nlev_cas),qv_nudg_prof_cas(nlev_cas))
        allocate(u_nudg_prof_cas(nlev_cas),v_nudg_prof_cas(nlev_cas))
        allocate(ht_prof_cas(nlev_cas))
        allocate(hth_prof_cas(nlev_cas))
        allocate(hq_prof_cas(nlev_cas))
        allocate(hu_prof_cas(nlev_cas))
        allocate(hv_prof_cas(nlev_cas))
        allocate(vt_prof_cas(nlev_cas))
        allocate(vth_prof_cas(nlev_cas))
        allocate(vq_prof_cas(nlev_cas))
        allocate(vu_prof_cas(nlev_cas))
        allocate(vv_prof_cas(nlev_cas))
        allocate(dt_prof_cas(nlev_cas))
        allocate(dth_prof_cas(nlev_cas))
        allocate(dtrad_prof_cas(nlev_cas))
        allocate(dq_prof_cas(nlev_cas))
        allocate(du_prof_cas(nlev_cas))
        allocate(dv_prof_cas(nlev_cas))
        allocate(uw_prof_cas(nlev_cas))
        allocate(vw_prof_cas(nlev_cas))
        allocate(q1_prof_cas(nlev_cas))
        allocate(q2_prof_cas(nlev_cas))

        print*,'Allocations OK'
        CALL read_SCM (nid,nlev_cas,nt_cas,                                                                     &
     &     ap_cas,bp_cas,z_cas,plev_cas,zh_cas,plevh_cas,t_cas,th_cas,thv_cas,thl_cas,qv_cas,                   &
     &     ql_cas,qi_cas,rh_cas,rv_cas,u_cas,v_cas,vitw_cas,omega_cas,ug_cas,vg_cas,                            &
     &     temp_nudg_cas,qv_nudg_cas,u_nudg_cas,v_nudg_cas,                                                     &
     &     du_cas,hu_cas,vu_cas,                                                                                &
     &     dv_cas,hv_cas,vv_cas,dt_cas,ht_cas,vt_cas,dq_cas,hq_cas,vq_cas,dth_cas,hth_cas,vth_cas,               &
     &     dr_cas,hr_cas,vr_cas,dtrad_cas,sens_cas,lat_cas,ts_cas,ps_cas,ustar_cas,tke_cas,                      &
     &     uw_cas,vw_cas,q1_cas,q2_cas,orog_cas,albedo_cas,emiss_cas,t_skin_cas,q_skin_cas,mom_rough,heat_rough, &
     &     o3_cas,rugos_cas,clay_cas,sand_cas)
        print*,'read_SCM cas OK'
        do ii=1,nlev_cas
        print*,'apres read2_SCM, plev_cas=',ii,plev_cas(ii,1)
        !print*,'apres read_SCM, plev_cas=',ii,omega_cas(ii,nt_cas/2+1)
        enddo


END SUBROUTINE read_SCM_cas


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE deallocate2_1D_cases
!profils environnementaux:
        deallocate(plev_cas,plevh_cas)
        
        deallocate(z_cas,zh_cas)
        deallocate(ap_cas,bp_cas)
        deallocate(t_cas,q_cas,qv_cas,ql_cas,qi_cas,rh_cas)
        deallocate(th_cas,thl_cas,thv_cas,rv_cas)
        deallocate(u_cas,v_cas,vitw_cas,omega_cas)
        
!forcing
        deallocate(ht_cas,vt_cas,dt_cas,dtrad_cas)
        deallocate(hq_cas,vq_cas,dq_cas)
        deallocate(hth_cas,vth_cas,dth_cas)
        deallocate(hr_cas,vr_cas,dr_cas)
        deallocate(hu_cas,vu_cas,du_cas)
        deallocate(hv_cas,vv_cas,dv_cas)
        deallocate(ug_cas)
        deallocate(vg_cas)
        deallocate(lat_cas,sens_cas,ts_cas,ps_cas,ustar_cas,tke_cas,uw_cas,vw_cas,q1_cas,q2_cas)

!champs interpoles
        deallocate(plev_prof_cas)
        deallocate(t_prof_cas)
        deallocate(theta_prof_cas)
        deallocate(thl_prof_cas)
        deallocate(thv_prof_cas)
        deallocate(q_prof_cas)
        deallocate(qv_prof_cas)
        deallocate(ql_prof_cas)
        deallocate(qi_prof_cas)
        deallocate(rh_prof_cas)
        deallocate(rv_prof_cas)
        deallocate(u_prof_cas)
        deallocate(v_prof_cas)
        deallocate(vitw_prof_cas)
        deallocate(omega_prof_cas)
        deallocate(ug_prof_cas)
        deallocate(vg_prof_cas)
        deallocate(temp_nudg_prof_cas,qv_nudg_prof_cas,u_nudg_prof_cas,v_nudg_prof_cas)
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
        deallocate(u_prof_cas)
        deallocate(v_prof_cas)
        deallocate(uw_prof_cas)
        deallocate(vw_prof_cas)
        deallocate(q1_prof_cas)
        deallocate(q2_prof_cas)

END SUBROUTINE deallocate2_1D_cases


!=====================================================================
      SUBROUTINE read_SCM(nid,nlevel,ntime,                                       &
     &     ap,bp,zz,pp,zzh,pph,temp,theta,thv,thl,qv,ql,qi,rh,rv,u,v,vitw,omega,ug,vg,&
     &     temp_nudg,qv_nudg,u_nudg,v_nudg,                                        &
     &     du,hu,vu,dv,hv,vv,dt,ht,vt,dq,hq,vq,                                    &
     &     dth,hth,vth,dr,hr,vr,dtrad,sens,flat,ts,ps,ustar,tket,uw,vw,q1,q2,       &
     &     orog_cas,albedo_cas,emiss_cas,t_skin_cas,q_skin_cas,mom_rough,          &
     &     heat_rough,o3_cas,rugos_cas,clay_cas,sand_cas)

!program reading forcing of the case study
      implicit none
#include "netcdf.inc"
#include "compar1d.h"

      integer ntime,nlevel,k,t

      real ap(nlevel+1),bp(nlevel+1)
      real zz(nlevel,ntime),zzh(nlevel+1)
      real pp(nlevel,ntime),pph(nlevel+1)
!profils initiaux
      real temp0(nlevel),qv0(nlevel),ql0(nlevel),qi0(nlevel),u0(nlevel),v0(nlevel),tke0(nlevel)
      real pp0(nlevel)    
      real temp(nlevel,ntime),qv(nlevel,ntime),ql(nlevel,ntime),qi(nlevel,ntime),rh(nlevel,ntime)
      real theta(nlevel,ntime),thv(nlevel,ntime),thl(nlevel,ntime),rv(nlevel,ntime)
      real u(nlevel,ntime),v(nlevel,ntime),tket(ntime)
      real temp_nudg(nlevel,ntime),qv_nudg(nlevel,ntime),u_nudg(nlevel,ntime),v_nudg(nlevel,ntime)
      real ug(nlevel,ntime),vg(nlevel,ntime)
      real vitw(nlevel,ntime),omega(nlevel,ntime)
      real du(nlevel,ntime),hu(nlevel,ntime),vu(nlevel,ntime)
      real dv(nlevel,ntime),hv(nlevel,ntime),vv(nlevel,ntime)
      real dt(nlevel,ntime),ht(nlevel,ntime),vt(nlevel,ntime)
      real dtrad(nlevel,ntime)
      real dq(nlevel,ntime),hq(nlevel,ntime),vq(nlevel,ntime)
      real dth(nlevel,ntime),hth(nlevel,ntime),vth(nlevel,ntime),hthl(nlevel,ntime)
      real dr(nlevel,ntime),hr(nlevel,ntime),vr(nlevel,ntime)
      real flat(ntime),sens(ntime),ustar(ntime)
      real uw(nlevel,ntime),vw(nlevel,ntime),q1(nlevel,ntime),q2(nlevel,ntime)
      real ts(ntime),ps(ntime)
      real orog_cas,albedo_cas,emiss_cas,t_skin_cas,q_skin_cas,mom_rough,heat_rough,o3_cas,rugos_cas,clay_cas,sand_cas
      real apbp(nlevel+1),resul(nlevel,ntime),resul1(nlevel),resul2(ntime),resul3


      integer nid, ierr,ierr1,ierr2,rid,i
      integer nbvar3d
      parameter(nbvar3d=74)
      integer var3didin(nbvar3d),missing_var(nbvar3d)
      character*13 name_var(1:nbvar3d)


      data name_var/ &
     ! coordonnees pression (n+1 niveaux) #4
     & 'coor_par_a','coor_par_b','height_h','pressure_h',& ! #1-#4
     ! coordonnees pression (n niveaux) #8
     &'temp','qv','ql','qi','u','v','tke','pressure',& ! #5-#12
     ! coordonnees pression + temps #42
     &'w','omega','ug','vg','uadv','uadvh','uadvv','vadv','vadvh','vadvv','temp_adv','tadvh','tadvv',& !  #13 - #25
     &'qv_adv','qvadvh','qvadvv','thadv','thadvh','thadvv','thladvh',                             & ! #26 - #32
     & 'radv','radvh','radvv','radcool','q1','q2','ustress','vstress',                           & ! #33 - #40
     & 'rh','temp_nudging','qv_nudging','u_nudging','v_nudging',                                       & ! #41-45
     &'height_f','pressure_forc','tempt','theta','thv','thl','qvt','qlt','qit','rv','ut','vt',   & ! #46-58
     ! coordonnees temps #12
     &'tket','sfc_sens_flx','sfc_lat_flx','ts','ps','ustar',&
     &'orog','albedo','emiss','t_skin','q_skin','mom_rough','heat_rough',&
     ! scalaires #4
     &'o3','rugos','clay','sand'/

!-----------------------------------------------------------------------
! Checking availability of variable #i in the cas.nc file
!     missing_var=1 if the variable is missing
!-----------------------------------------------------------------------

       do i=1,nbvar3d
         missing_var(i)=0.
         ierr=NF_INQ_VARID(nid,name_var(i),var3didin(i)) 
         if(ierr/=NF_NOERR) then
           print *,'Variable manquante dans cas.nc:',i,name_var(i)
           ierr=NF_NOERR
           missing_var(i)=1
         else

!-----------------------------------------------------------------------
! Activating keys depending on the presence of specific variables in cas.nc
!-----------------------------------------------------------------------
if ( 1 == 1 ) THEN
! A MODIFIER: il faudrait dire nudging_temp mais faut le declarer dans compar1d.h etc...        
!           if ( name_var(i) == 'temp_nudging' .and. nint(nudging_t)==0) stop 'Nudging inconsistency temp'
            if ( name_var(i) == 'qv_nudging' .and. nint(nudging_qv)==0) stop 'Nudging inconsistency qv'
            if ( name_var(i) == 'u_nudging' .and. nint(nudging_u)==0) stop 'Nudging inconsistency u'
            if ( name_var(i) == 'v_nudging' .and. nint(nudging_u)==0) stop 'Nudging inconsistency v'
    ELSE
             print*,'GUIDAGE : CONSISTENCY CHECK DEACTIVATED FOR TESTS of SANDU/REF'
    ENDIF

!-----------------------------------------------------------------------
! Reading variables 1D (N+1) vertical variables (nlevelp1,lat,lon)
!-----------------------------------------------------------------------
           if(i.LE.4) then
#ifdef NC_DOUBLE
           ierr = NF_GET_VAR_DOUBLE(nid,var3didin(i),apbp)
#else
           ierr = NF_GET_VAR_REAL(nid,var3didin(i),apbp)
#endif
           print *,'read2_cas(apbp), on a lu ',i,name_var(i)
           if(ierr/=NF_NOERR) then
              print *,'Pb a la lecture de cas.nc: ',name_var(i)
              stop "getvarup"
           endif

!-----------------------------------------------------------------------
!  Reading 1D (N) vertical varialbes    (nlevel,lat,lon)   
!-----------------------------------------------------------------------
           else if(i.gt.4.and.i.LE.12) then  
#ifdef NC_DOUBLE
           ierr = NF_GET_VAR_DOUBLE(nid,var3didin(i),resul1)
#else
           ierr = NF_GET_VAR_REAL(nid,var3didin(i),resul1)
#endif
           print *,'read2_cas(resul1), on a lu ',i,name_var(i)
           if(ierr/=NF_NOERR) then
              print *,'Pb a la lecture de cas.nc: ',name_var(i)
              stop "getvarup"
           endif
         print*,'Lecture de la variable #i ',i,name_var(i),minval(resul1),maxval(resul1)

!-----------------------------------------------------------------------
!  Reading 2D tim-vertical variables  (time,nlevel,lat,lon)
!  TBD : seems to be the same as above.
!-----------------------------------------------------------------------
           else if(i.gt.12.and.i.LE.57) then
#ifdef NC_DOUBLE
           ierr = NF_GET_VAR_DOUBLE(nid,var3didin(i),resul)
#else
           ierr = NF_GET_VAR_REAL(nid,var3didin(i),resul)
#endif
           print *,'read2_cas(resul), on a lu ',i,name_var(i)
           if(ierr/=NF_NOERR) then
              print *,'Pb a la lecture de cas.nc: ',name_var(i)
              stop "getvarup"
           endif
         print*,'Lecture de la variable #i ',i,name_var(i),minval(resul),maxval(resul)

!-----------------------------------------------------------------------
!  Reading 1D time variables (time,lat,lon)
!-----------------------------------------------------------------------
           else if (i.gt.57.and.i.LE.63) then
#ifdef NC_DOUBLE
           ierr = NF_GET_VAR_DOUBLE(nid,var3didin(i),resul2)
#else
           ierr = NF_GET_VAR_REAL(nid,var3didin(i),resul2)
#endif
           print *,'read2_cas(resul2), on a lu ',i,name_var(i)
           if(ierr/=NF_NOERR) then
              print *,'Pb a la lecture de cas.nc: ',name_var(i)
              stop "getvarup"
           endif
         print*,'Lecture de la variable #i  ',i,name_var(i),minval(resul2),maxval(resul2)

!-----------------------------------------------------------------------
! Reading scalar variables (lat,lon)
!-----------------------------------------------------------------------
           else
#ifdef NC_DOUBLE
           ierr = NF_GET_VAR_DOUBLE(nid,var3didin(i),resul3)
#else
           ierr = NF_GET_VAR_REAL(nid,var3didin(i),resul3)
#endif
           print *,'read2_cas(resul3), on a lu ',i,name_var(i)
           if(ierr/=NF_NOERR) then
              print *,'Pb a la lecture de cas.nc: ',name_var(i)
              stop "getvarup"
           endif
         print*,'Lecture de la variable #i ',i,name_var(i),resul3
           endif
         endif

!-----------------------------------------------------------------------
! Attributing variables
!-----------------------------------------------------------------------
         select case(i)
         !case(1) ; ap=apbp       ! donnees indexees en nlevel+1
         ! case(2) ; bp=apbp
           case(3) ; zzh=apbp
           case(4) ; pph=apbp
           case(5) ; temp0=resul1    ! donnees initiales
           case(6) ; qv0=resul1
           case(7) ; ql0=resul1
           case(8) ; qi0=resul1
           case(9) ; u0=resul1
           case(10) ; v0=resul1
           case(11) ; tke0=resul1
           case(12) ; pp0=resul1
           case(13) ; vitw=resul    ! donnees indexees en nlevel,time
           case(14) ; omega=resul
           case(15) ; ug=resul
           case(16) ; vg=resul
           case(17) ; du=resul
           case(18) ; hu=resul
           case(19) ; vu=resul
           case(20) ; dv=resul
           case(21) ; hv=resul
           case(22) ; vv=resul
           case(23) ; dt=resul
           case(24) ; ht=resul
           case(25) ; vt=resul
           case(26) ; dq=resul
           case(27) ; hq=resul
           case(28) ; vq=resul
           case(29) ; dth=resul
           case(30) ; hth=resul
           case(31) ; vth=resul
           case(32) ; hthl=resul
           case(33) ; dr=resul
           case(34) ; hr=resul
           case(35) ; vr=resul
           case(36) ; dtrad=resul
           case(37) ; q1=resul
           case(38) ; q2=resul
           case(39) ; uw=resul
           case(40) ; vw=resul
           case(41) ; rh=resul
           case(42) ; temp_nudg=resul
           case(43) ; qv_nudg=resul
           case(44) ; u_nudg=resul
           case(45) ; v_nudg=resul
           case(46) ; zz=resul      ! donnees en time,nlevel pour profil initial
           case(47) ; pp=resul
           case(48) ; temp=resul
           case(49) ; theta=resul
           case(50) ; thv=resul
           case(51) ; thl=resul
           case(52) ; qv=resul
           case(53) ; ql=resul
           case(54) ; qi=resul
           case(55) ; rv=resul
           case(56) ; u=resul
           case(57) ; v=resul
           case(58) ; tket=resul2   ! donnees indexees en time
           case(59) ; sens=resul2
           case(60) ; flat=resul2
           case(61) ; ts=resul2
           case(62) ; ps=resul2
           case(63) ; ustar=resul2
           case(64) ; orog_cas=resul3      ! constantes
           case(65) ; albedo_cas=resul3
           case(66) ; emiss_cas=resul3
           case(67) ; t_skin_cas=resul3
           case(68) ; q_skin_cas=resul3
           case(69) ; mom_rough=resul3
           case(70) ; heat_rough=resul3
           case(71) ; o3_cas=resul3        
           case(72) ; rugos_cas=resul3
           case(73) ; clay_cas=resul3
           case(74) ; sand_cas=resul3
         end select
         resul=0.
         resul1=0.
         resul2=0.
         resul3=0.
       enddo
         print*,'Lecture de la variable APRES ,sens ',minval(sens),maxval(sens)
         print*,'Lecture de la variable APRES ,flat ',minval(flat),maxval(flat)

!CR:ATTENTION EN ATTENTE DE REGLER LA QUESTION DU PAS DE TEMPS INITIAL
       do t=1,ntime
          do k=1,nlevel
             temp(k,t)=temp0(k)
             qv(k,t)=qv0(k)
             ql(k,t)=ql0(k)
             qi(k,t)=qi0(k)
             u(k,t)=u0(k)
             v(k,t)=v0(k)
             !tke(k,t)=tke0(k)
          enddo
       enddo
       !!!! TRAVAIL : EN FONCTION DES DECISIONS SUR LA SPECIFICATION DE W
       !!!omega=-vitw*pres*rg/(rd*temp)
!-----------------------------------------------------------------------

         return 
         END SUBROUTINE read_SCM
!======================================================================

!======================================================================

!**********************************************************************************************
        SUBROUTINE interp_case_time_std(day,day1,annee_ref                           &
!    &         ,year_cas,day_cas,nt_cas,pdt_forc,nlev_cas                         &
     &         ,nt_cas,nlev_cas                                                   &
     &         ,ts_cas,ps_cas,plev_cas,t_cas,theta_cas,thv_cas,thl_cas            &
     &         ,qv_cas,ql_cas,qi_cas,u_cas,v_cas                                  &
     &         ,ug_cas,vg_cas,temp_nudg_cas,qv_nudg_cas,u_nudg_cas,v_nudg_cas     &
     &         ,vitw_cas,omega_cas,du_cas,hu_cas,vu_cas             &
     &         ,dv_cas,hv_cas,vv_cas,dt_cas,ht_cas,vt_cas,dtrad_cas               &
     &         ,dq_cas,hq_cas,vq_cas,dth_cas,hth_cas,vth_cas                      &
     &         ,lat_cas,sens_cas,ustar_cas                                        &
     &         ,uw_cas,vw_cas,q1_cas,q2_cas,tke_cas                               &
!
     &         ,ts_prof_cas,ps_prof_cas,plev_prof_cas,t_prof_cas,theta_prof_cas               &
     &         ,thv_prof_cas,thl_prof_cas,qv_prof_cas,ql_prof_cas,qi_prof_cas     &
     &         ,u_prof_cas,v_prof_cas,ug_prof_cas,vg_prof_cas                     &
     &         ,temp_nudg_prof_cas,qv_nudg_prof_cas,u_nudg_prof_cas,v_nudg_prof_cas     &
     &         ,vitw_prof_cas,omega_prof_cas,du_prof_cas,hu_prof_cas,vu_prof_cas  &
     &         ,dv_prof_cas,hv_prof_cas,vv_prof_cas,dt_prof_cas                   &
     &         ,ht_prof_cas,vt_prof_cas,dtrad_prof_cas,dq_prof_cas                &
     &         ,hq_prof_cas,vq_prof_cas,dth_prof_cas,hth_prof_cas,vth_prof_cas    &
     &         ,lat_prof_cas,sens_prof_cas                                        &
     &         ,ustar_prof_cas,uw_prof_cas,vw_prof_cas,q1_prof_cas,q2_prof_cas,tke_prof_cas)
          

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
        real ts_cas(nt_cas),ps_cas(nt_cas)
        real plev_cas(nlev_cas,nt_cas)
        real t_cas(nlev_cas,nt_cas),theta_cas(nlev_cas,nt_cas),thv_cas(nlev_cas,nt_cas),thl_cas(nlev_cas,nt_cas)
        real qv_cas(nlev_cas,nt_cas),ql_cas(nlev_cas,nt_cas),qi_cas(nlev_cas,nt_cas)
        real u_cas(nlev_cas,nt_cas),v_cas(nlev_cas,nt_cas)
        real ug_cas(nlev_cas,nt_cas),vg_cas(nlev_cas,nt_cas)
        real temp_nudg_cas(nlev_cas,nt_cas),qv_nudg_cas(nlev_cas,nt_cas)
        real u_nudg_cas(nlev_cas,nt_cas),v_nudg_cas(nlev_cas,nt_cas)

        real vitw_cas(nlev_cas,nt_cas),omega_cas(nlev_cas,nt_cas)
        real du_cas(nlev_cas,nt_cas),hu_cas(nlev_cas,nt_cas),vu_cas(nlev_cas,nt_cas)
        real dv_cas(nlev_cas,nt_cas),hv_cas(nlev_cas,nt_cas),vv_cas(nlev_cas,nt_cas)
        real dt_cas(nlev_cas,nt_cas),ht_cas(nlev_cas,nt_cas),vt_cas(nlev_cas,nt_cas)
        real dth_cas(nlev_cas,nt_cas),hth_cas(nlev_cas,nt_cas),vth_cas(nlev_cas,nt_cas)
        real dtrad_cas(nlev_cas,nt_cas)
        real dq_cas(nlev_cas,nt_cas),hq_cas(nlev_cas,nt_cas),vq_cas(nlev_cas,nt_cas)
        real lat_cas(nt_cas),sens_cas(nt_cas),tke_cas(nt_cas)
        real ustar_cas(nt_cas),uw_cas(nlev_cas,nt_cas),vw_cas(nlev_cas,nt_cas)
        real q1_cas(nlev_cas,nt_cas),q2_cas(nlev_cas,nt_cas)

! outputs:
        real plev_prof_cas(nlev_cas)
        real t_prof_cas(nlev_cas),theta_prof_cas(nlev_cas),thl_prof_cas(nlev_cas),thv_prof_cas(nlev_cas)
        real qv_prof_cas(nlev_cas),ql_prof_cas(nlev_cas),qi_prof_cas(nlev_cas)
        real u_prof_cas(nlev_cas),v_prof_cas(nlev_cas)
        real ug_prof_cas(nlev_cas),vg_prof_cas(nlev_cas)
        real temp_nudg_prof_cas(nlev_cas),qv_nudg_prof_cas(nlev_cas)
        real u_nudg_prof_cas(nlev_cas),v_nudg_prof_cas(nlev_cas)

        real vitw_prof_cas(nlev_cas),omega_prof_cas(nlev_cas)
        real du_prof_cas(nlev_cas),hu_prof_cas(nlev_cas),vu_prof_cas(nlev_cas)
        real dv_prof_cas(nlev_cas),hv_prof_cas(nlev_cas),vv_prof_cas(nlev_cas)
        real dt_prof_cas(nlev_cas),ht_prof_cas(nlev_cas),vt_prof_cas(nlev_cas)
        real dth_prof_cas(nlev_cas),hth_prof_cas(nlev_cas),vth_prof_cas(nlev_cas)
        real dtrad_prof_cas(nlev_cas)
        real dq_prof_cas(nlev_cas),hq_prof_cas(nlev_cas),vq_prof_cas(nlev_cas)
        real lat_prof_cas,sens_prof_cas,tke_prof_cas,ts_prof_cas,ps_prof_cas,ustar_prof_cas
        real uw_prof_cas(nlev_cas),vw_prof_cas(nlev_cas),q1_prof_cas(nlev_cas),q2_prof_cas(nlev_cas)
! local:
        integer it_cas1, it_cas2,k
        real timeit,time_cas1,time_cas2,frac


        print*,'Check time',day1,day_ju_ini_cas,day_deb+1,pdt_cas
!       do k=1,nlev_cas
!       print*,'debut de interp2_case_time, plev_cas=',k,plev_cas(k,1)
!       enddo

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
!     print *,'timeit,pdt_cas,nt_cas=',timeit,pdt_cas,nt_cas
!     print *,'it_cas1,it_cas2,time_cas1,time_cas2=',it_cas1,it_cas2,time_cas1,time_cas2

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

       lat_prof_cas = lat_cas(it_cas2)                                   &
     &          -frac*(lat_cas(it_cas2)-lat_cas(it_cas1)) 
       sens_prof_cas = sens_cas(it_cas2)                                 &
     &          -frac*(sens_cas(it_cas2)-sens_cas(it_cas1))
       tke_prof_cas = tke_cas(it_cas2)                                   &
     &          -frac*(tke_cas(it_cas2)-tke_cas(it_cas1))
       ts_prof_cas = ts_cas(it_cas2)                                     &
     &          -frac*(ts_cas(it_cas2)-ts_cas(it_cas1))
       ps_prof_cas = ps_cas(it_cas2)                                     &
     &          -frac*(ps_cas(it_cas2)-ps_cas(it_cas1))
       ustar_prof_cas = ustar_cas(it_cas2)                               &
     &          -frac*(ustar_cas(it_cas2)-ustar_cas(it_cas1))

       do k=1,nlev_cas
        plev_prof_cas(k) = plev_cas(k,it_cas2)                           &     
     &          -frac*(plev_cas(k,it_cas2)-plev_cas(k,it_cas1))
        t_prof_cas(k) = t_cas(k,it_cas2)                                 &        
     &          -frac*(t_cas(k,it_cas2)-t_cas(k,it_cas1))
        !print *,'k,frac,plev_cas1,plev_cas2=',k,frac,plev_cas(k,it_cas1),plev_cas(k,it_cas2)
        theta_prof_cas(k) = theta_cas(k,it_cas2)                         &                      
     &          -frac*(theta_cas(k,it_cas2)-theta_cas(k,it_cas1))
        thv_prof_cas(k) = thv_cas(k,it_cas2)                             &          
     &          -frac*(thv_cas(k,it_cas2)-thv_cas(k,it_cas1))
        thl_prof_cas(k) = thl_cas(k,it_cas2)                             &              
     &          -frac*(thl_cas(k,it_cas2)-thl_cas(k,it_cas1))
        qv_prof_cas(k) = qv_cas(k,it_cas2)                               &
     &          -frac*(qv_cas(k,it_cas2)-qv_cas(k,it_cas1))
        ql_prof_cas(k) = ql_cas(k,it_cas2)                               &
     &          -frac*(ql_cas(k,it_cas2)-ql_cas(k,it_cas1))
        qi_prof_cas(k) = qi_cas(k,it_cas2)                               &
     &          -frac*(qi_cas(k,it_cas2)-qi_cas(k,it_cas1))
        u_prof_cas(k) = u_cas(k,it_cas2)                                 &
     &          -frac*(u_cas(k,it_cas2)-u_cas(k,it_cas1))
        v_prof_cas(k) = v_cas(k,it_cas2)                                 &
     &          -frac*(v_cas(k,it_cas2)-v_cas(k,it_cas1))
        ug_prof_cas(k) = ug_cas(k,it_cas2)                               &
     &          -frac*(ug_cas(k,it_cas2)-ug_cas(k,it_cas1))
        vg_prof_cas(k) = vg_cas(k,it_cas2)                               &
     &          -frac*(vg_cas(k,it_cas2)-vg_cas(k,it_cas1))
        temp_nudg_prof_cas(k) = temp_nudg_cas(k,it_cas2)                    &
     &          -frac*(temp_nudg_cas(k,it_cas2)-temp_nudg_cas(k,it_cas1))
        qv_nudg_prof_cas(k) = qv_nudg_cas(k,it_cas2)                        &
     &          -frac*(qv_nudg_cas(k,it_cas2)-qv_nudg_cas(k,it_cas1))
        u_nudg_prof_cas(k) = u_nudg_cas(k,it_cas2)                          &
     &          -frac*(u_nudg_cas(k,it_cas2)-u_nudg_cas(k,it_cas1))
        v_nudg_prof_cas(k) = v_nudg_cas(k,it_cas2)                          &
     &          -frac*(v_nudg_cas(k,it_cas2)-v_nudg_cas(k,it_cas1))
        vitw_prof_cas(k) = vitw_cas(k,it_cas2)                           &
     &          -frac*(vitw_cas(k,it_cas2)-vitw_cas(k,it_cas1))
        omega_prof_cas(k) = omega_cas(k,it_cas2)                         &
     &          -frac*(omega_cas(k,it_cas2)-omega_cas(k,it_cas1))
        du_prof_cas(k) = du_cas(k,it_cas2)                               &
     &          -frac*(du_cas(k,it_cas2)-du_cas(k,it_cas1))
        hu_prof_cas(k) = hu_cas(k,it_cas2)                               &
     &          -frac*(hu_cas(k,it_cas2)-hu_cas(k,it_cas1))
        vu_prof_cas(k) = vu_cas(k,it_cas2)                               &
     &          -frac*(vu_cas(k,it_cas2)-vu_cas(k,it_cas1))
        dv_prof_cas(k) = dv_cas(k,it_cas2)                               &
     &          -frac*(dv_cas(k,it_cas2)-dv_cas(k,it_cas1))
        hv_prof_cas(k) = hv_cas(k,it_cas2)                               &
     &          -frac*(hv_cas(k,it_cas2)-hv_cas(k,it_cas1))
        vv_prof_cas(k) = vv_cas(k,it_cas2)                               &
     &          -frac*(vv_cas(k,it_cas2)-vv_cas(k,it_cas1))
        dt_prof_cas(k) = dt_cas(k,it_cas2)                               &
     &          -frac*(dt_cas(k,it_cas2)-dt_cas(k,it_cas1))
        ht_prof_cas(k) = ht_cas(k,it_cas2)                               &
     &          -frac*(ht_cas(k,it_cas2)-ht_cas(k,it_cas1))
        vt_prof_cas(k) = vt_cas(k,it_cas2)                               &
     &          -frac*(vt_cas(k,it_cas2)-vt_cas(k,it_cas1))
        dth_prof_cas(k) = dth_cas(k,it_cas2)                             &
     &          -frac*(dth_cas(k,it_cas2)-dth_cas(k,it_cas1))
        hth_prof_cas(k) = hth_cas(k,it_cas2)                             &
     &          -frac*(hth_cas(k,it_cas2)-hth_cas(k,it_cas1))
        vth_prof_cas(k) = vth_cas(k,it_cas2)                             &
     &          -frac*(vth_cas(k,it_cas2)-vth_cas(k,it_cas1))
        dtrad_prof_cas(k) = dtrad_cas(k,it_cas2)                         &
     &          -frac*(dtrad_cas(k,it_cas2)-dtrad_cas(k,it_cas1))
        dq_prof_cas(k) = dq_cas(k,it_cas2)                               &
     &          -frac*(dq_cas(k,it_cas2)-dq_cas(k,it_cas1))
        hq_prof_cas(k) = hq_cas(k,it_cas2)                               &
     &          -frac*(hq_cas(k,it_cas2)-hq_cas(k,it_cas1))
        vq_prof_cas(k) = vq_cas(k,it_cas2)                               &
     &          -frac*(vq_cas(k,it_cas2)-vq_cas(k,it_cas1))
       uw_prof_cas(k) = uw_cas(k,it_cas2)                                &
     &          -frac*(uw_cas(k,it_cas2)-uw_cas(k,it_cas1))
       vw_prof_cas(k) = vw_cas(k,it_cas2)                                &
     &          -frac*(vw_cas(k,it_cas2)-vw_cas(k,it_cas1))
       q1_prof_cas(k) = q1_cas(k,it_cas2)                                &
     &          -frac*(q1_cas(k,it_cas2)-q1_cas(k,it_cas1))
       q2_prof_cas(k) = q2_cas(k,it_cas2)                                &
     &          -frac*(q2_cas(k,it_cas2)-q2_cas(k,it_cas1))
        enddo

        return
        END SUBROUTINE interp_case_time_std

!**********************************************************************************************
!=====================================================================
       SUBROUTINE interp2_case_vertical_std(play,nlev_cas,plev_prof_cas                                &
     &         ,t_prof_cas,th_prof_cas,thv_prof_cas,thl_prof_cas                                       &
     &         ,qv_prof_cas,ql_prof_cas,qi_prof_cas,u_prof_cas,v_prof_cas                              &
     &         ,ug_prof_cas,vg_prof_cas                                                                &
     &         ,temp_nudg_prof_cas,qv_nudg_prof_cas,u_nudg_prof_cas,v_nudg_prof_cas                    &
     &         ,vitw_prof_cas,omega_prof_cas                                                           &
     &         ,du_prof_cas,hu_prof_cas,vu_prof_cas,dv_prof_cas,hv_prof_cas,vv_prof_cas                &
     &         ,dt_prof_cas,ht_prof_cas,vt_prof_cas,dtrad_prof_cas,dq_prof_cas,hq_prof_cas,vq_prof_cas &
     &         ,dth_prof_cas,hth_prof_cas,vth_prof_cas                                                 &
!
     &         ,t_mod_cas,theta_mod_cas,thv_mod_cas,thl_mod_cas                                        &
     &         ,qv_mod_cas,ql_mod_cas,qi_mod_cas,u_mod_cas,v_mod_cas                                   &
     &         ,ug_mod_cas,vg_mod_cas                                                                  &
     &         ,temp_nudg_mod_cas,qv_nudg_mod_cas,u_nudg_mod_cas,v_nudg_mod_cas                        &
     &         ,w_mod_cas,omega_mod_cas                                                                &
     &         ,du_mod_cas,hu_mod_cas,vu_mod_cas,dv_mod_cas,hv_mod_cas,vv_mod_cas                      &
     &         ,dt_mod_cas,ht_mod_cas,vt_mod_cas,dtrad_mod_cas,dq_mod_cas,hq_mod_cas,vq_mod_cas        &
     &         ,dth_mod_cas,hth_mod_cas,vth_mod_cas,mxcalc)
 
       implicit none
 
#include "YOMCST.h"
#include "dimensions.h"

!-------------------------------------------------------------------------
! Vertical interpolation of generic case forcing data onto mod_casel levels
!-------------------------------------------------------------------------
 
       integer nlevmax
       parameter (nlevmax=41)
       integer nlev_cas,mxcalc
!       real play(llm), plev_prof(nlevmax) 
!       real t_prof(nlevmax),q_prof(nlevmax)
!       real u_prof(nlevmax),v_prof(nlevmax), w_prof(nlevmax)
!       real ht_prof(nlevmax),vt_prof(nlevmax)
!       real hq_prof(nlevmax),vq_prof(nlevmax)
 
       real play(llm), plev_prof_cas(nlev_cas) 
       real t_prof_cas(nlev_cas),th_prof_cas(nlev_cas),thv_prof_cas(nlev_cas),thl_prof_cas(nlev_cas)
       real qv_prof_cas(nlev_cas),ql_prof_cas(nlev_cas),qi_prof_cas(nlev_cas)
       real u_prof_cas(nlev_cas),v_prof_cas(nlev_cas)
       real ug_prof_cas(nlev_cas),vg_prof_cas(nlev_cas), vitw_prof_cas(nlev_cas),omega_prof_cas(nlev_cas)
       real temp_nudg_prof_cas(nlev_cas),qv_nudg_prof_cas(nlev_cas)
       real u_nudg_prof_cas(nlev_cas),v_nudg_prof_cas(nlev_cas)

       real du_prof_cas(nlev_cas),hu_prof_cas(nlev_cas),vu_prof_cas(nlev_cas)
       real dv_prof_cas(nlev_cas),hv_prof_cas(nlev_cas),vv_prof_cas(nlev_cas)
       real dt_prof_cas(nlev_cas),ht_prof_cas(nlev_cas),vt_prof_cas(nlev_cas),dtrad_prof_cas(nlev_cas)
       real dth_prof_cas(nlev_cas),hth_prof_cas(nlev_cas),vth_prof_cas(nlev_cas)
       real dq_prof_cas(nlev_cas),hq_prof_cas(nlev_cas),vq_prof_cas(nlev_cas)
 
       real t_mod_cas(llm),theta_mod_cas(llm),thv_mod_cas(llm),thl_mod_cas(llm)
       real qv_mod_cas(llm),ql_mod_cas(llm),qi_mod_cas(llm)
       real u_mod_cas(llm),v_mod_cas(llm)
       real ug_mod_cas(llm),vg_mod_cas(llm), w_mod_cas(llm),omega_mod_cas(llm)
       real temp_nudg_mod_cas(llm),qv_nudg_mod_cas(llm)
       real u_nudg_mod_cas(llm),v_nudg_mod_cas(llm)
       real du_mod_cas(llm),hu_mod_cas(llm),vu_mod_cas(llm)
       real dv_mod_cas(llm),hv_mod_cas(llm),vv_mod_cas(llm)
       real dt_mod_cas(llm),ht_mod_cas(llm),vt_mod_cas(llm),dtrad_mod_cas(llm)
       real dth_mod_cas(llm),hth_mod_cas(llm),vth_mod_cas(llm)
       real dq_mod_cas(llm),hq_mod_cas(llm),vq_mod_cas(llm)
 
       integer l,k,k1,k2
       real frac,frac1,frac2,fact
 
!       do l = 1, llm
!       print *,'debut interp2, play=',l,play(l)
!       enddo
!      do l = 1, nlev_cas
!      print *,'debut interp2, plev_prof_cas=',l,play(l),plev_prof_cas(l)
!      enddo

       do l = 1, llm

        if (play(l).ge.plev_prof_cas(nlev_cas)) then
 
        mxcalc=l
!        print *,'debut interp2, mxcalc=',mxcalc
         k1=0
         k2=0

         if (play(l).le.plev_prof_cas(1)) then

         do k = 1, nlev_cas-1
          if (play(l).le.plev_prof_cas(k).and. play(l).gt.plev_prof_cas(k+1)) then
            k1=k
            k2=k+1
          endif
         enddo

         if (k1.eq.0 .or. k2.eq.0) then
          write(*,*) 'PB! k1, k2 = ',k1,k2
          write(*,*) 'l,play(l) = ',l,play(l)/100
         do k = 1, nlev_cas-1
          write(*,*) 'k,plev_prof_cas(k) = ',k,plev_prof_cas(k)/100
         enddo
         endif

         frac = (plev_prof_cas(k2)-play(l))/(plev_prof_cas(k2)-plev_prof_cas(k1))
         t_mod_cas(l)= t_prof_cas(k2) - frac*(t_prof_cas(k2)-t_prof_cas(k1))
         theta_mod_cas(l)= th_prof_cas(k2) - frac*(th_prof_cas(k2)-th_prof_cas(k1))
         if(theta_mod_cas(l).NE.0) t_mod_cas(l)= theta_mod_cas(l)*(play(l)/100000.)**(RD/RCPD)
         thv_mod_cas(l)= thv_prof_cas(k2) - frac*(thv_prof_cas(k2)-thv_prof_cas(k1))
         thl_mod_cas(l)= thl_prof_cas(k2) - frac*(thl_prof_cas(k2)-thl_prof_cas(k1))
         qv_mod_cas(l)= qv_prof_cas(k2) - frac*(qv_prof_cas(k2)-qv_prof_cas(k1))
         ql_mod_cas(l)= ql_prof_cas(k2) - frac*(ql_prof_cas(k2)-ql_prof_cas(k1))
         qi_mod_cas(l)= qi_prof_cas(k2) - frac*(qi_prof_cas(k2)-qi_prof_cas(k1))
         u_mod_cas(l)= u_prof_cas(k2) - frac*(u_prof_cas(k2)-u_prof_cas(k1))
         v_mod_cas(l)= v_prof_cas(k2) - frac*(v_prof_cas(k2)-v_prof_cas(k1))
         ug_mod_cas(l)= ug_prof_cas(k2) - frac*(ug_prof_cas(k2)-ug_prof_cas(k1))
         vg_mod_cas(l)= vg_prof_cas(k2) - frac*(vg_prof_cas(k2)-vg_prof_cas(k1))
         temp_nudg_mod_cas(l)= temp_nudg_prof_cas(k2) - frac*(temp_nudg_prof_cas(k2)-temp_nudg_prof_cas(k1))
         qv_nudg_mod_cas(l)= qv_nudg_prof_cas(k2) - frac*(qv_nudg_prof_cas(k2)-qv_nudg_prof_cas(k1))
         u_nudg_mod_cas(l)= u_nudg_prof_cas(k2) - frac*(u_nudg_prof_cas(k2)-u_nudg_prof_cas(k1))
         v_nudg_mod_cas(l)= v_nudg_prof_cas(k2) - frac*(v_nudg_prof_cas(k2)-v_nudg_prof_cas(k1))
         w_mod_cas(l)= vitw_prof_cas(k2) - frac*(vitw_prof_cas(k2)-vitw_prof_cas(k1))
         omega_mod_cas(l)= omega_prof_cas(k2) - frac*(omega_prof_cas(k2)-omega_prof_cas(k1))
         du_mod_cas(l)= du_prof_cas(k2) - frac*(du_prof_cas(k2)-du_prof_cas(k1))
         hu_mod_cas(l)= hu_prof_cas(k2) - frac*(hu_prof_cas(k2)-hu_prof_cas(k1))
         vu_mod_cas(l)= vu_prof_cas(k2) - frac*(vu_prof_cas(k2)-vu_prof_cas(k1))
         dv_mod_cas(l)= dv_prof_cas(k2) - frac*(dv_prof_cas(k2)-dv_prof_cas(k1))
         hv_mod_cas(l)= hv_prof_cas(k2) - frac*(hv_prof_cas(k2)-hv_prof_cas(k1))
         vv_mod_cas(l)= vv_prof_cas(k2) - frac*(vv_prof_cas(k2)-vv_prof_cas(k1))
         dt_mod_cas(l)= dt_prof_cas(k2) - frac*(dt_prof_cas(k2)-dt_prof_cas(k1))
         ht_mod_cas(l)= ht_prof_cas(k2) - frac*(ht_prof_cas(k2)-ht_prof_cas(k1))
         vt_mod_cas(l)= vt_prof_cas(k2) - frac*(vt_prof_cas(k2)-vt_prof_cas(k1))
         dth_mod_cas(l)= dth_prof_cas(k2) - frac*(dth_prof_cas(k2)-dth_prof_cas(k1))
         hth_mod_cas(l)= hth_prof_cas(k2) - frac*(hth_prof_cas(k2)-hth_prof_cas(k1))
         vth_mod_cas(l)= vth_prof_cas(k2) - frac*(vth_prof_cas(k2)-vth_prof_cas(k1))
         dq_mod_cas(l)= dq_prof_cas(k2) - frac*(dq_prof_cas(k2)-dq_prof_cas(k1))
         hq_mod_cas(l)= hq_prof_cas(k2) - frac*(hq_prof_cas(k2)-hq_prof_cas(k1))
         vq_mod_cas(l)= vq_prof_cas(k2) - frac*(vq_prof_cas(k2)-vq_prof_cas(k1))
         dtrad_mod_cas(l)= dtrad_prof_cas(k2) - frac*(dtrad_prof_cas(k2)-dtrad_prof_cas(k1))
     
         else !play>plev_prof_cas(1)

         k1=1
         k2=2
         print *,'interp2_vert, k1,k2=',plev_prof_cas(k1),plev_prof_cas(k2)
         frac1 = (play(l)-plev_prof_cas(k2))/(plev_prof_cas(k1)-plev_prof_cas(k2))
         frac2 = (play(l)-plev_prof_cas(k1))/(plev_prof_cas(k1)-plev_prof_cas(k2))
         t_mod_cas(l)= frac1*t_prof_cas(k1) - frac2*t_prof_cas(k2)
         theta_mod_cas(l)= frac1*th_prof_cas(k1) - frac2*th_prof_cas(k2)
         if(theta_mod_cas(l).NE.0) t_mod_cas(l)= theta_mod_cas(l)*(play(l)/100000.)**(RD/RCPD)
         thv_mod_cas(l)= frac1*thv_prof_cas(k1) - frac2*thv_prof_cas(k2)
         thl_mod_cas(l)= frac1*thl_prof_cas(k1) - frac2*thl_prof_cas(k2)
         qv_mod_cas(l)= frac1*qv_prof_cas(k1) - frac2*qv_prof_cas(k2)
         ql_mod_cas(l)= frac1*ql_prof_cas(k1) - frac2*ql_prof_cas(k2)
         qi_mod_cas(l)= frac1*qi_prof_cas(k1) - frac2*qi_prof_cas(k2)
         u_mod_cas(l)= frac1*u_prof_cas(k1) - frac2*u_prof_cas(k2)
         v_mod_cas(l)= frac1*v_prof_cas(k1) - frac2*v_prof_cas(k2)
         ug_mod_cas(l)= frac1*ug_prof_cas(k1) - frac2*ug_prof_cas(k2)
         vg_mod_cas(l)= frac1*vg_prof_cas(k1) - frac2*vg_prof_cas(k2)
         temp_nudg_mod_cas(l)= frac1*temp_nudg_prof_cas(k1) - frac2*temp_nudg_prof_cas(k2)
         qv_nudg_mod_cas(l)= frac1*qv_nudg_prof_cas(k1) - frac2*qv_nudg_prof_cas(k2)
         u_nudg_mod_cas(l)= frac1*u_nudg_prof_cas(k1) - frac2*u_nudg_prof_cas(k2)
         v_nudg_mod_cas(l)= frac1*v_nudg_prof_cas(k1) - frac2*v_nudg_prof_cas(k2)
         w_mod_cas(l)= frac1*vitw_prof_cas(k1) - frac2*vitw_prof_cas(k2)
         omega_mod_cas(l)= frac1*omega_prof_cas(k1) - frac2*omega_prof_cas(k2)
         du_mod_cas(l)= frac1*du_prof_cas(k1) - frac2*du_prof_cas(k2)
         hu_mod_cas(l)= frac1*hu_prof_cas(k1) - frac2*hu_prof_cas(k2)
         vu_mod_cas(l)= frac1*vu_prof_cas(k1) - frac2*vu_prof_cas(k2)
         dv_mod_cas(l)= frac1*dv_prof_cas(k1) - frac2*dv_prof_cas(k2)
         hv_mod_cas(l)= frac1*hv_prof_cas(k1) - frac2*hv_prof_cas(k2)
         vv_mod_cas(l)= frac1*vv_prof_cas(k1) - frac2*vv_prof_cas(k2)
         dt_mod_cas(l)= frac1*dt_prof_cas(k1) - frac2*dt_prof_cas(k2)
         ht_mod_cas(l)= frac1*ht_prof_cas(k1) - frac2*ht_prof_cas(k2)
         vt_mod_cas(l)= frac1*vt_prof_cas(k1) - frac2*vt_prof_cas(k2)
         dth_mod_cas(l)= frac1*dth_prof_cas(k1) - frac2*dth_prof_cas(k2)
         hth_mod_cas(l)= frac1*hth_prof_cas(k1) - frac2*hth_prof_cas(k2)
         vth_mod_cas(l)= frac1*vth_prof_cas(k1) - frac2*vth_prof_cas(k2)
         dq_mod_cas(l)= frac1*dq_prof_cas(k1) - frac2*dq_prof_cas(k2)
         hq_mod_cas(l)= frac1*hq_prof_cas(k1) - frac2*hq_prof_cas(k2)
         vq_mod_cas(l)= frac1*vq_prof_cas(k1) - frac2*vq_prof_cas(k2)
         dtrad_mod_cas(l)= frac1*dtrad_prof_cas(k1) - frac2*dtrad_prof_cas(k2)

         endif ! play.le.plev_prof_cas(1)

        else ! above max altitude of forcing file
 
!jyg
         fact=20.*(plev_prof_cas(nlev_cas)-play(l))/plev_prof_cas(nlev_cas) !jyg
         fact = max(fact,0.)                                           !jyg
         fact = exp(-fact)                                             !jyg
         t_mod_cas(l)= t_prof_cas(nlev_cas)                            !jyg
         theta_mod_cas(l)= th_prof_cas(nlev_cas)                       !jyg
         thv_mod_cas(l)= thv_prof_cas(nlev_cas)                        !jyg
         thl_mod_cas(l)= thl_prof_cas(nlev_cas)                        !jyg
         qv_mod_cas(l)= qv_prof_cas(nlev_cas)*fact                     !jyg
         ql_mod_cas(l)= ql_prof_cas(nlev_cas)*fact                     !jyg
         qi_mod_cas(l)= qi_prof_cas(nlev_cas)*fact                     !jyg
         u_mod_cas(l)= u_prof_cas(nlev_cas)*fact                       !jyg
         v_mod_cas(l)= v_prof_cas(nlev_cas)*fact                       !jyg
         ug_mod_cas(l)= ug_prof_cas(nlev_cas)                          !jyg
         vg_mod_cas(l)= vg_prof_cas(nlev_cas)                          !jyg
         temp_nudg_mod_cas(l)= temp_nudg_prof_cas(nlev_cas)            !jyg
         qv_nudg_mod_cas(l)= qv_nudg_prof_cas(nlev_cas)                !jyg
         u_nudg_mod_cas(l)= u_nudg_prof_cas(nlev_cas)                  !jyg
         v_nudg_mod_cas(l)= v_nudg_prof_cas(nlev_cas)                  !jyg
         thv_mod_cas(l)= thv_prof_cas(nlev_cas)                        !jyg
         w_mod_cas(l)= 0.0                                             !jyg
         omega_mod_cas(l)= 0.0                                         !jyg
         du_mod_cas(l)= du_prof_cas(nlev_cas)*fact 
         hu_mod_cas(l)= hu_prof_cas(nlev_cas)*fact                     !jyg
         vu_mod_cas(l)= vu_prof_cas(nlev_cas)*fact                     !jyg
         dv_mod_cas(l)= dv_prof_cas(nlev_cas)*fact
         hv_mod_cas(l)= hv_prof_cas(nlev_cas)*fact                     !jyg
         vv_mod_cas(l)= vv_prof_cas(nlev_cas)*fact                     !jyg
         dt_mod_cas(l)= dt_prof_cas(nlev_cas) 
         ht_mod_cas(l)= ht_prof_cas(nlev_cas)                          !jyg
         vt_mod_cas(l)= vt_prof_cas(nlev_cas)                          !jyg
         dth_mod_cas(l)= dth_prof_cas(nlev_cas) 
         hth_mod_cas(l)= hth_prof_cas(nlev_cas)                        !jyg
         vth_mod_cas(l)= vth_prof_cas(nlev_cas)                        !jyg
         dq_mod_cas(l)= dq_prof_cas(nlev_cas)*fact 
         hq_mod_cas(l)= hq_prof_cas(nlev_cas)*fact                     !jyg
         vq_mod_cas(l)= vq_prof_cas(nlev_cas)*fact                     !jyg
         dtrad_mod_cas(l)= dtrad_prof_cas(nlev_cas)*fact               !jyg
 
        endif ! play
 
       enddo ! l

          return
          end SUBROUTINE interp2_case_vertical_std
!***************************************************************************** 





END MODULE mod_1D_cases_read_std
