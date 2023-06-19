MODULE call_dissip_mod

    REAL,POINTER,SAVE :: ucov(:,:)
    REAL,POINTER,SAVE :: vcov(:,:)
    REAL,POINTER,SAVE :: teta(:,:)
    REAL,POINTER,SAVE :: p(:,: )
    REAL,POINTER,SAVE :: pk(:,:)

    REAL,POINTER,SAVE :: ucont(:,:)
    REAL,POINTER,SAVE :: vcont(:,:)
    REAL,POINTER,SAVE :: ecin(:,:)
    REAL,POINTER,SAVE :: ecin0(:,:)
    REAL,POINTER,SAVE :: dudis(:,:)
    REAL,POINTER,SAVE :: dvdis(:,:)
    REAL,POINTER,SAVE :: dtetadis(:,:)
    REAL,POINTER,SAVE :: dtetaecdt(:,:)



CONTAINS
  
  SUBROUTINE call_dissip_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  USE dissip_mod, ONLY : dissip_allocate
  IMPLICIT NONE
    TYPE(distrib),POINTER :: d
    d=>distrib_dissip

    CALL allocate_u(ucov,llm,d)
    ucov(:,:)=0
    CALL allocate_v(vcov,llm,d)
    vcov(:,:)=0
    CALL allocate_u(teta,llm,d)
    CALL allocate_u(p,llmp1,d)
    CALL allocate_u(pk,llm,d)
    CALL allocate_u(ucont,llm,d)
    CALL allocate_v(vcont,llm,d)
    CALL allocate_u(ecin,llm,d)
    CALL allocate_u(ecin0,llm,d)
    CALL allocate_u(dudis,llm,d)
    CALL allocate_v(dvdis,llm,d)
    CALL allocate_u(dtetadis,llm,d)
    CALL allocate_u(dtetaecdt,llm,d)
    
    
    CALL dissip_allocate
    
  END SUBROUTINE call_dissip_allocate
  
  SUBROUTINE call_dissip_switch_dissip(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  USE dissip_mod, ONLY : dissip_switch_dissip
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL switch_u(ucov,distrib_dissip,dist)
    CALL switch_v(vcov,distrib_dissip,dist)
    CALL switch_u(teta,distrib_dissip,dist)
    CALL switch_u(p,distrib_dissip,dist)
    CALL switch_u(pk,distrib_dissip,dist)
    CALL switch_u(ucont,distrib_dissip,dist)
    CALL switch_v(vcont,distrib_dissip,dist)
    CALL switch_u(ecin,distrib_dissip,dist)
    CALL switch_u(ecin0,distrib_dissip,dist)
    CALL switch_u(dudis,distrib_dissip,dist)
    CALL switch_v(dvdis,distrib_dissip,dist)
    CALL switch_u(dtetadis,distrib_dissip,dist)
    CALL switch_u(dtetaecdt,distrib_dissip,dist)

    CALL dissip_switch_dissip(dist)
    
  END SUBROUTINE call_dissip_switch_dissip  
  

  
  SUBROUTINE call_dissip(ucov_dyn,vcov_dyn,teta_dyn,p_dyn,pk_dyn,ps_dyn)
  USE dimensions_mod
  USE parallel_lmdz
  USE times
  USE mod_hallo
  USE Bands
  USE vampir
  USE write_field_loc
  IMPLICIT NONE
    INCLUDE 'comgeom.h'
    REAL,INTENT(INOUT) :: ucov_dyn(ijb_u:ije_u,llm) ! covariant zonal wind
    REAL,INTENT(INOUT) :: vcov_dyn(ijb_v:ije_v,llm) ! covariant meridional wind
    REAL,INTENT(INOUT) :: teta_dyn(ijb_u:ije_u,llm) ! covariant meridional wind
    REAL,INTENT(INOUT) :: p_dyn(ijb_u:ije_u,llmp1 ) ! pressure at interlayer
    REAL,INTENT(INOUT) :: pk_dyn(ijb_u:ije_u,llm) ! Exner at midlayer
    REAL,INTENT(INOUT) :: ps_dyn(ijb_u:ije_u) ! surface pressure
    REAL :: tppn(iim),tpps(iim)
    REAL :: tpn,tps

    REAL  SSUM
    LOGICAL,PARAMETER :: dissip_conservative=.TRUE.
    TYPE(Request),SAVE :: Request_dissip 
!$OMP THREADPRIVATE(Request_dissip )    
    INTEGER :: ij,l,ijb,ije 
  
    
  !$OMP MASTER
    CALL suspend_timer(timer_caldyn)
        
!       print*,'Entree dans la dissipation : Iteration No ',true_itau
!   calcul de l'energie cinetique avant dissipation
!       print *,'Passage dans la dissipation'

    CALL VTb(VThallo)
  !$OMP END MASTER

  !$OMP BARRIER

    CALL Register_SwapField_u(ucov_dyn,ucov,distrib_dissip, Request_dissip,up=1,down=1)
    CALL Register_SwapField_v(vcov_dyn,vcov,distrib_dissip, Request_dissip,up=1,down=1)
    CALL Register_SwapField_u(teta_dyn,teta,distrib_dissip, Request_dissip)
    CALL Register_SwapField_u(p_dyn,p,distrib_dissip,Request_dissip)
    CALL Register_SwapField_u(pk_dyn,pk,distrib_dissip,Request_dissip)

    CALL SendRequest(Request_dissip)       
  !$OMP BARRIER
    CALL WaitRequest(Request_dissip)       

  !$OMP BARRIER
  !$OMP MASTER
    CALL set_distrib(distrib_dissip)
    CALL VTe(VThallo)
    CALL VTb(VTdissipation)
    CALL start_timer(timer_dissip)
  !$OMP END MASTER
  !$OMP BARRIER

    CALL covcont_loc(llm,ucov,vcov,ucont,vcont)
    CALL enercin_loc(vcov,ucov,vcont,ucont,ecin0)

!   dissipation

!        CALL FTRACE_REGION_BEGIN("dissip")
    CALL dissip_loc(vcov,ucov,teta,p,dvdis,dudis,dtetadis)

#ifdef DEBUG_IO    
    CALL WriteField_u('dudis',dudis)
    CALL WriteField_v('dvdis',dvdis)
    CALL WriteField_u('dtetadis',dtetadis)
#endif
 
!      CALL FTRACE_REGION_END("dissip")
         
    ijb=ij_begin
    ije=ij_end
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)        
    DO l=1,llm
      ucov(ijb:ije,l)=ucov(ijb:ije,l)+dudis(ijb:ije,l)
    ENDDO
  !$OMP END DO NOWAIT        

    IF (pole_sud) ije=ije-iip1
   
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)        
    DO l=1,llm
      vcov(ijb:ije,l)=vcov(ijb:ije,l)+dvdis(ijb:ije,l)
    ENDDO
  !$OMP END DO NOWAIT        

!       teta=teta+dtetadis


!------------------------------------------------------------------------
    IF (dissip_conservative) THEN
!       On rajoute la tendance due a la transform. Ec -> E therm. cree
!       lors de la dissipation
    !$OMP BARRIER
    !$OMP MASTER
      CALL suspend_timer(timer_dissip)
      CALL VTb(VThallo)
    !$OMP END MASTER
      CALL Register_Hallo_u(ucov,llm,1,1,1,1,Request_Dissip)
      CALL Register_Hallo_v(vcov,llm,1,1,1,1,Request_Dissip)
      CALL SendRequest(Request_Dissip)
    !$OMP BARRIER
      CALL WaitRequest(Request_Dissip)
    !$OMP MASTER
      CALL VTe(VThallo)
      CALL resume_timer(timer_dissip)
    !$OMP END MASTER
    !$OMP BARRIER            
      CALL covcont_loc(llm,ucov,vcov,ucont,vcont)
      CALL enercin_loc(vcov,ucov,vcont,ucont,ecin)
            
      ijb=ij_begin
      ije=ij_end
    !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
      DO l=1,llm
        DO ij=ijb,ije
           dtetaecdt(ij,l)= (ecin0(ij,l)-ecin(ij,l))/ pk(ij,l)
           dtetadis(ij,l)=dtetadis(ij,l)+dtetaecdt(ij,l)
        ENDDO
      ENDDO
    !$OMP END DO NOWAIT            

    ENDIF

    ijb=ij_begin
    ije=ij_end

  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
    DO l=1,llm
      DO ij=ijb,ije
         teta(ij,l)=teta(ij,l)+dtetadis(ij,l)
      ENDDO
    ENDDO
  !$OMP END DO NOWAIT         

!------------------------------------------------------------------------


!    .......        P. Le Van (  ajout  le 17/04/96  )   ...........
!   ...      Calcul de la valeur moyenne, unique de h aux poles  .....
!

    ijb=ij_begin
    ije=ij_end
         
    IF (pole_nord) THEN
  
   !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l  =  1, llm
        DO ij =  1,iim
          tppn(ij)  = aire(  ij    ) * teta(  ij    ,l)
        ENDDO
        tpn  = SSUM(iim,tppn,1)/apoln

        DO ij = 1, iip1
          teta(  ij    ,l) = tpn
        ENDDO
      ENDDO
    !$OMP END DO NOWAIT

         if (1 == 0) then
!!! Ehouarn: lines here 1) kill 1+1=2 in the dynamics
!!!                     2) should probably not be here anyway
!!! but are kept for those who would want to revert to previous behaviour
    !$OMP MASTER               
      DO ij =  1,iim
        tppn(ij)  = aire(  ij    ) * ps_dyn (  ij    )
      ENDDO
      tpn  = SSUM(iim,tppn,1)/apoln
  
      DO ij = 1, iip1
        ps_dyn(  ij    ) = tpn
      ENDDO
    !$OMP END MASTER
    
    ENDIF ! of if (1 == 0)
    endif ! of of (pole_nord)
        
    IF (pole_sud) THEN

    !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l  =  1, llm
        DO ij =  1,iim
          tpps(ij)  = aire(ij+ip1jm) * teta(ij+ip1jm,l)
        ENDDO
        
        tps  = SSUM(iim,tpps,1)/apols

        DO ij = 1, iip1
          teta(ij+ip1jm,l) = tps
        ENDDO
      ENDDO
    !$OMP END DO NOWAIT

    if (1 == 0) then
!!! Ehouarn: lines here 1) kill 1+1=2 in the dynamics
!!!                     2) should probably not be here anyway
!!! but are kept for those who would want to revert to previous behaviour
    !$OMP MASTER               
      DO ij =  1,iim
        tpps(ij)  = aire(ij+ip1jm) * ps_dyn (ij+ip1jm)
      ENDDO
      tps  = SSUM(iim,tpps,1)/apols
  
      DO ij = 1, iip1
        ps_dyn(ij+ip1jm) = tps
      ENDDO
    !$OMP END MASTER
    ENDIF ! of if (1 == 0)
    endif ! of if (pole_sud)


  !$OMP BARRIER
  !$OMP MASTER
    CALL VTe(VTdissipation)
    CALL stop_timer(timer_dissip)
    CALL VTb(VThallo)
  !$OMP END MASTER
 
    CALL Register_SwapField_u(ucov,ucov_dyn,distrib_caldyn,Request_dissip)
    CALL Register_SwapField_v(vcov,vcov_dyn,distrib_caldyn,Request_dissip)
    CALL Register_SwapField_u(teta,teta_dyn,distrib_caldyn,Request_dissip)
    CALL Register_SwapField_u(p,p_dyn,distrib_caldyn,Request_dissip)
    CALL Register_SwapField_u(pk,pk_dyn,distrib_caldyn,Request_dissip)

    CALL SendRequest(Request_dissip)       

  !$OMP BARRIER
    CALL WaitRequest(Request_dissip)       
  !$OMP BARRIER
  !$OMP MASTER
    CALL set_distrib(distrib_caldyn)
    CALL VTe(VThallo)
    CALL resume_timer(timer_caldyn)
!        print *,'fin dissipation'
  !$OMP END MASTER
  !$OMP BARRIER
  
  
  END SUBROUTINE call_dissip

END MODULE call_dissip_mod
