!#define DEBUG_IO
MODULE call_calfis_mod

    REAL,POINTER,SAVE :: ucov(:,:)
    REAL,POINTER,SAVE :: vcov(:,:) 
    REAL,POINTER,SAVE :: teta(:,:) 
    REAL,POINTER,SAVE :: masse(:,:) 
    REAL,POINTER,SAVE :: ps(:) 
    REAL,POINTER,SAVE :: phis(:) 
    REAL,POINTER,SAVE :: q(:,:,:) 
    REAL,POINTER,SAVE :: flxw(:,:) 

    REAL,POINTER,SAVE :: p(:,:) 
    REAL,POINTER,SAVE :: pks(:) 
    REAL,POINTER,SAVE :: pk(:,:) 
    REAL,POINTER,SAVE :: pkf(:,:) 
    REAL,POINTER,SAVE :: phi(:,:) 
    REAL,POINTER,SAVE :: du(:,:) 
    REAL,POINTER,SAVE :: dv(:,:) 
    REAL,POINTER,SAVE :: dteta(:,:) 
    REAL,POINTER,SAVE :: dq(:,:,:) 
    REAL,POINTER,SAVE :: dufi(:,:) 
    REAL,POINTER,SAVE :: dvfi(:,:) 
    REAL,POINTER,SAVE :: dtetafi(:,:) 
    REAL,POINTER,SAVE :: dqfi(:,:,:) 
    REAL,POINTER,SAVE :: dpfi(:) 
   
    
    
    
    
CONTAINS

  SUBROUTINE call_calfis_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  USE infotrac
  IMPLICIT NONE
    TYPE(distrib),POINTER :: d
    d=>distrib_physic

    CALL allocate_u(ucov,llm,d)
    CALL allocate_v(vcov,llm,d)
    CALL allocate_u(teta,llm,d)
    CALL allocate_u(masse,llm,d)
    CALL allocate_u(ps,d)
    CALL allocate_u(phis,d)
    CALL allocate_u(q,llm,nqtot,d)
    CALL allocate_u(flxw,llm,d)
    CALL allocate_u(p,llmp1,d)
    CALL allocate_u(pks,d)
    pks(:)=0
    CALL allocate_u(pk,llm,d)
    pk(:,:)=0
    CALL allocate_u(pkf,llm,d)
    CALL allocate_u(phi,llm,d)
    CALL allocate_u(du,llm,d)
    CALL allocate_v(dv,llm,d)
    CALL allocate_u(dteta,llm,d)
    CALL allocate_u(dq,llm,nqtot,d)
    CALL allocate_u(dufi,llm,d)
    CALL allocate_v(dvfi,llm,d)
    CALL allocate_u(dtetafi,llm,d)
    CALL allocate_u(dqfi,llm,nqtot,d)
    CALL allocate_u(dpfi,d)
  
  END SUBROUTINE call_calfis_allocate
  
  
  SUBROUTINE call_calfis(itau,lafin,ucov_dyn,vcov_dyn,teta_dyn,masse_dyn,ps_dyn, &
                         phis_dyn,q_dyn,flxw_dyn)
  USE dimensions_mod
  use exner_hyb_loc_m, only: exner_hyb_loc
  use exner_milieu_loc_m, only: exner_milieu_loc
  USE parallel_lmdz
  USE times
  USE mod_hallo
  USE Bands
  USE vampir
  USE infotrac
  USE control_mod
  USE write_field_loc
  USE write_field
  USE comconst_mod, ONLY: dtphys
  USE logic_mod, ONLY: leapf, forward, ok_strato
  USE comvert_mod, ONLY: ap, bp, pressure_exner
  USE temps_mod, ONLY: day_ini, day_ref, jd_ref, jh_ref, start_time
  
  IMPLICIT NONE
    INCLUDE "iniprint.h"

    INTEGER,INTENT(IN) :: itau ! (time) iteration step number
    LOGICAL,INTENT(IN) :: lafin ! .true. if final time step
    REAL,INTENT(INOUT) :: ucov_dyn(ijb_u:ije_u,llm) ! covariant zonal wind
    REAL,INTENT(INOUT) :: vcov_dyn(ijb_v:ije_v,llm) ! covariant meridional wind
    REAL,INTENT(INOUT) :: teta_dyn(ijb_u:ije_u,llm) ! potential temperature
    REAL,INTENT(INOUT) :: masse_dyn(ijb_u:ije_u,llm) ! air mass
    REAL,INTENT(INOUT) :: ps_dyn(ijb_u:ije_u) ! surface pressure
    REAL,INTENT(INOUT) :: phis_dyn(ijb_u:ije_u) ! surface geopotential
    REAL,INTENT(INOUT) :: q_dyn(ijb_u:ije_u,llm,nqtot) ! advected tracers
    REAL,INTENT(INOUT) :: flxw_dyn(ijb_u:ije_u,llm) ! vertical mass flux

    REAL :: dufi_tmp(iip1,llm)    
    REAL :: dvfi_tmp(iip1,llm)  
    REAL :: dtetafi_tmp(iip1,llm)
    REAL :: dpfi_tmp(iip1)
    REAL :: dqfi_tmp(iip1,llm,nqtot)

    REAL :: jD_cur, jH_cur
    CHARACTER(LEN=15) :: ztit
    TYPE(Request),SAVE :: Request_physic
!$OMP THREADPRIVATE(Request_physic )
    INTEGER :: ijb,ije,l,j
    
    
#ifdef DEBUG_IO    
    CALL WriteField_u('ucovfi',ucov)
    CALL WriteField_v('vcovfi',vcov)
    CALL WriteField_u('tetafi',teta)
    CALL WriteField_u('pfi',p)
    CALL WriteField_u('pkfi',pk)
    DO j=1,nqtot
      CALL WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
    ENDDO
#endif

!
!     .......   Ajout   P.Le Van ( 17/04/96 )   ...........
!


  !$OMP MASTER
    CALL suspend_timer(timer_caldyn)
    IF (prt_level >= 10) THEN
      WRITE(lunout,*) 'leapfrog_p: Entree dans la physique : Iteration No ',itau
    ENDIF
  !$OMP END MASTER
   
           jD_cur = jD_ref + day_ini - day_ref                           &
     &        + (itau+1)/day_step

           IF (planet_type .eq."generic") THEN
              ! AS: we make jD_cur to be pday
              jD_cur = int(day_ini + itau/day_step)
           ENDIF

           jH_cur = jH_ref + start_time +                                &
     &              mod(itau+1,day_step)/float(day_step) 
    if (jH_cur > 1.0 ) then
      jD_cur = jD_cur +1.
      jH_cur = jH_cur -1.
    endif

!   Inbterface avec les routines de phylmd (phymars ... )
!   -----------------------------------------------------

!+jld

!  Diagnostique de conservation de l'energie : initialisation
 
!-jld
  !$OMP BARRIER
  !$OMP MASTER
    CALL VTb(VThallo)
  !$OMP END MASTER

#ifdef DEBUG_IO    
    CALL WriteField_u('ucovfi',ucov)
    CALL WriteField_v('vcovfi',vcov)
    CALL WriteField_u('tetafi',teta)
    CALL WriteField_u('pfi',p)
    CALL WriteField_u('pkfi',pk)
#endif
    
    CALL SetTag(Request_physic,800)
    CALL Register_SwapField_u(ucov_dyn,ucov,distrib_physic,Request_physic,up=2,down=2)
    CALL Register_SwapField_v(vcov_dyn,vcov,distrib_physic,Request_physic,up=2,down=2)
    CALL Register_SwapField_u(teta_dyn,teta,distrib_physic,Request_physic,up=2,down=2)
    CALL Register_SwapField_u(masse_dyn,masse,distrib_physic,Request_physic,up=1,down=2)
    CALL Register_SwapField_u(ps_dyn,ps,distrib_physic,Request_physic,up=2,down=2)
    CALL Register_SwapField_u(phis_dyn,phis,distrib_physic,Request_physic,up=2,down=2)
    CALL Register_SwapField_u(q_dyn,q,distrib_physic,Request_physic,up=2,down=2)
    CALL Register_SwapField_u(flxw_dyn,flxw,distrib_physic,Request_physic,up=2,down=2)
 
    CALL SendRequest(Request_Physic)
  !$OMP BARRIER
    CALL WaitRequest(Request_Physic)       

  !$OMP BARRIER
  !$OMP MASTER
    CALL Set_Distrib(distrib_Physic)
    CALL VTe(VThallo)
        
    CALL VTb(VTphysiq)
  !$OMP END MASTER
  !$OMP BARRIER

    CALL pression_loc (  ip1jmp1, ap, bp, ps,  p      )

  !$OMP BARRIER
    CALL exner_hyb_loc(  ip1jmp1, ps, p, pks, pk, pkf )
  !$OMP BARRIER
    CALL geopot_loc  ( ip1jmp1, teta  , pk , pks,  phis  , phi   )


    CALL Register_Hallo_u(p,llmp1,2,2,2,2,Request_physic)
    CALL Register_Hallo_u(pk,llm,2,2,2,2,Request_physic)
    CALL Register_Hallo_u(phi,llm,2,2,2,2,Request_physic)
        
    CALL SendRequest(Request_Physic)
  !$OMP BARRIER
    CALL WaitRequest(Request_Physic)
             
  !$OMP BARRIER
  
  
#ifdef DEBUG_IO    
    CALL WriteField_u('ucovfi',ucov)
    CALL WriteField_v('vcovfi',vcov)
    CALL WriteField_u('tetafi',teta)
    CALL WriteField_u('pfi',p)
    CALL WriteField_u('pkfi',pk)
    DO j=1,nqtot
      CALL WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
    ENDDO
#endif

  !$OMP BARRIER

#ifdef CPP_PHYS
    CALL calfis_loc(lafin ,jD_cur, jH_cur,                       &
                     ucov,vcov,teta,q,masse,ps,p,pk,phis,phi ,   &
                     du,dv,dteta,dq,                             &
                     flxw, dufi,dvfi,dtetafi,dqfi,dpfi  )
#endif
    ijb=ij_begin
    ije=ij_end  
    IF ( .not. pole_nord) THEN
  
    !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
        dufi_tmp(1:iip1,l)   = dufi(ijb:ijb+iim,l) 
        dvfi_tmp(1:iip1,l)   = dvfi(ijb:ijb+iim,l)  
        dtetafi_tmp(1:iip1,l)= dtetafi(ijb:ijb+iim,l)  
        dqfi_tmp(1:iip1,l,:) = dqfi(ijb:ijb+iim,l,:)  
      ENDDO
    !$OMP END DO NOWAIT

    !$OMP MASTER
      dpfi_tmp(1:iip1)     = dpfi(ijb:ijb+iim)  
    !$OMP END MASTER
    
    ENDIF ! of if ( .not. pole_nord)

  !$OMP BARRIER
  !$OMP MASTER
    CALL Set_Distrib(distrib_Physic_bis)
    CALL VTb(VThallo)
  !$OMP END MASTER
  !$OMP BARRIER
 
    CALL Register_Hallo_u(dufi,llm,1,0,0,1,Request_physic)
    CALL Register_Hallo_v(dvfi,llm,1,0,0,1,Request_physic)
    CALL Register_Hallo_u(dtetafi,llm,1,0,0,1,Request_physic)
    CALL Register_Hallo_u(dpfi,1,1,0,0,1,Request_physic)

    DO j=1,nqtot
      CALL Register_Hallo_u(dqfi(:,:,j),llm,1,0,0,1,Request_physic)
    ENDDO
        
    CALL SendRequest(Request_Physic)
  !$OMP BARRIER
    CALL WaitRequest(Request_Physic)
             
  !$OMP BARRIER
  !$OMP MASTER
    CALL VTe(VThallo)
    CALL Set_Distrib(distrib_Physic)
  !$OMP END MASTER
  !$OMP BARRIER        
    ijb=ij_begin
    IF (.not. pole_nord) THEN
        
    !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
        dufi(ijb:ijb+iim,l) = dufi(ijb:ijb+iim,l)+dufi_tmp(1:iip1,l)
        dvfi(ijb:ijb+iim,l) = dvfi(ijb:ijb+iim,l)+dvfi_tmp(1:iip1,l) 
        dtetafi(ijb:ijb+iim,l) = dtetafi(ijb:ijb+iim,l)+dtetafi_tmp(1:iip1,l)
        dqfi(ijb:ijb+iim,l,:) = dqfi(ijb:ijb+iim,l,:) + dqfi_tmp(1:iip1,l,:)
      ENDDO
    !$OMP END DO NOWAIT

    !$OMP MASTER
      dpfi(ijb:ijb+iim)   = dpfi(ijb:ijb+iim)+ dpfi_tmp(1:iip1)
    !$OMP END MASTER
          
    endif ! of if (.not. pole_nord)
        
        
#ifdef DEBUG_IO           
    CALL WriteField_u('dufi',dufi)
    CALL WriteField_v('dvfi',dvfi) 
    CALL WriteField_u('dtetafi',dtetafi)
    CALL WriteField_u('dpfi',dpfi)
    DO j=1,nqtot
      CALL WriteField_u('dqfi'//trim(int2str(j)),dqfi(:,:,j))
    ENDDO
#endif

  !$OMP BARRIER

!      ajout des tendances physiques:
!      ------------------------------
#ifdef DEBUG_IO    
    CALL WriteField_u('ucovfi',ucov)
    CALL WriteField_v('vcovfi',vcov)
    CALL WriteField_u('tetafi',teta)
    CALL WriteField_u('psfi',ps)
    DO j=1,nqtot
      CALL WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
    ENDDO
#endif

#ifdef DEBUG_IO           
    CALL WriteField_u('ucovfi',ucov)
    CALL WriteField_v('vcovfi',vcov)
    CALL WriteField_u('tetafi',teta)
    CALL WriteField_u('psfi',ps)
    DO j=1,nqtot
      CALL WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
    ENDDO
#endif

    CALL addfi_loc( dtphys, leapf, forward   ,              &
                    ucov, vcov, teta , q   ,ps ,            &
                    dufi, dvfi, dtetafi , dqfi ,dpfi  )
    ! since addfi updates ps(), also update p(), masse() and pk()
    CALL pression_loc(ip1jmp1,ap,bp,ps,p)
!$OMP BARRIER
    CALL massdair_loc(p,masse)
!$OMP BARRIER
    if (pressure_exner) then
      CALL exner_hyb_loc(ijnb_u,ps,p,pks,pk,pkf)
    else 
      CALL exner_milieu_loc(ijnb_u,ps,p,pks,pk,pkf)
    endif
!$OMP BARRIER

#ifdef DEBUG_IO    
    CALL WriteField_u('ucovfi',ucov)
    CALL WriteField_v('vcovfi',vcov)
    CALL WriteField_u('tetafi',teta)
    CALL WriteField_u('psfi',ps)
    DO j=1,nqtot
      CALL WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
    ENDDO
#endif

    IF (ok_strato) THEN
!      CALL top_bound_loc( vcov,ucov,teta,masse,dufi,dvfi,dtetafi)
      CALL top_bound_loc(vcov,ucov,teta,masse,dtphys)
    ENDIF

  !$OMP BARRIER
  !$OMP MASTER
    CALL VTe(VTphysiq)
    CALL VTb(VThallo)
  !$OMP END MASTER

    CALL SetTag(Request_physic,800)
    CALL Register_SwapField_u(ucov,ucov_dyn,distrib_caldyn,Request_physic)
    CALL Register_SwapField_v(vcov,vcov_dyn,distrib_caldyn,Request_physic)
    CALL Register_SwapField_u(teta,teta_dyn,distrib_caldyn,Request_physic)
    CALL Register_SwapField_u(masse,masse_dyn,distrib_caldyn,Request_physic)
    CALL Register_SwapField_u(ps,ps_dyn,distrib_caldyn,Request_physic)
    CALL Register_SwapField_u(q,q_dyn,distrib_caldyn,Request_physic)
    CALL SendRequest(Request_Physic)
  !$OMP BARRIER
    CALL WaitRequest(Request_Physic)     

  !$OMP BARRIER
  !$OMP MASTER
    CALL VTe(VThallo)
    CALL set_distrib(distrib_caldyn)
  !$OMP END MASTER
  !$OMP BARRIER

!
!  Diagnostique de conservation de l'energie : difference
    IF (ip_ebil_dyn.ge.1 ) THEN 
      ztit='bil phys'
!      CALL diagedyn(ztit,2,1,1,dtphys,ucov, vcov , ps, p ,pk , teta , q(:,:,1), q(:,:,2))
      write(lunout,*)"call_calfis: diagedyn disabled in dyn3dmem !!"
    ENDIF 

#ifdef DEBUG_IO    
    CALL WriteField_u('ucovfi',ucov_dyn)
    CALL WriteField_v('vcovfi',vcov_dyn)
    CALL WriteField_u('tetafi',teta_dyn)
    CALL WriteField_u('psfi',ps_dyn)
    DO j=1,nqtot
      CALL WriteField_u('qfi'//trim(int2str(j)),q_dyn(:,:,j))
    ENDDO
#endif


!-jld
    !$OMP MASTER
      CALL resume_timer(timer_caldyn)
    !$OMP END MASTER

  END SUBROUTINE call_calfis
  
END MODULE call_calfis_mod
