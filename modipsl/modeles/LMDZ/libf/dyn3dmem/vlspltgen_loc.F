!
! $Header$
!
       SUBROUTINE vlspltgen_loc( q,iadv,pente_max,masse,w,pbaru,pbarv,
     &                           pdt, p,pk,teta                 )
     
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget, F.Codron 
c
c    ********************************************************************
c          Shema  d'advection " pseudo amont " .
c      + test sur humidite specifique: Q advecte< Qsat aval
c                   (F. Codron, 10/99)
c    ********************************************************************
c     q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c
c     pente_max facteur de limitation des pentes: 2 en general
c                                                0 pour un schema amont
c     pbaru,pbarv,w flux de masse en u ,v ,w
c     pdt pas de temps
c
c     teta temperature potentielle, p pression aux interfaces,
c     pk exner au milieu des couches necessaire pour calculer Qsat
c   --------------------------------------------------------------------
      USE parallel_lmdz
      USE mod_hallo
      USE Write_Field_loc
      USE VAMPIR
      ! CRisi: on rajoute variables utiles d'infotrac  
      USE infotrac, ONLY : nqtot,nqperes,nqdesc,nqfils,iqfils, 
     &    ok_iso_verif
      USE vlspltgen_mod
      USE comconst_mod, ONLY: cpp
      IMPLICIT NONE

c
      include "dimensions.h"
      include "paramet.h"

c
c   Arguments:
c   ----------
      INTEGER iadv(nqtot)
      REAL masse(ijb_u:ije_u,llm),pente_max
      REAL pbaru( ijb_u:ije_u,llm ),pbarv( ijb_v:ije_v,llm)
      REAL q(ijb_u:ije_u,llm,nqtot)
      REAL w(ijb_u:ije_u,llm),pdt
      REAL p(ijb_u:ije_u,llmp1),teta(ijb_u:ije_u,llm)
      REAL pk(ijb_u:ije_u,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l
c
      REAL zzpbar, zzw

      REAL qmin,qmax
      DATA qmin,qmax/0.,1.e33/

c--pour rapport de melange saturant--

      REAL rtt,retv,r2es,r3les,r3ies,r4les,r4ies,play
      REAL ptarg,pdelarg,foeew,zdelta
      REAL tempe(ijb_u:ije_u)
      INTEGER ijb,ije,iq,iq2,ifils
      LOGICAL, SAVE :: firstcall=.TRUE.
!$OMP THREADPRIVATE(firstcall)
      type(request),SAVE :: MyRequest1
!$OMP THREADPRIVATE(MyRequest1)
      type(request),SAVE :: MyRequest2
!$OMP THREADPRIVATE(MyRequest2)
c    fonction psat(T)

       FOEEW ( PTARG,PDELARG ) = EXP (
     *          (R3LES*(1.-PDELARG)+R3IES*PDELARG) * (PTARG-RTT)
     * / (PTARG-(R4LES*(1.-PDELARG)+R4IES*PDELARG)) )

        r2es  = 380.11733 
        r3les = 17.269
        r3ies = 21.875
        r4les = 35.86
        r4ies = 7.66
        retv = 0.6077667
        rtt  = 273.16

c Allocate variables depending on dynamic variable nqtot

         IF (firstcall) THEN
            firstcall=.FALSE.
         END IF
c-- Calcul de Qsat en chaque point
c-- approximation: au milieu des couches play(l)=(p(l)+p(l+1))/2
c   pour eviter une exponentielle.

      call SetTag(MyRequest1,100)
      call SetTag(MyRequest2,101)

        
        ijb=ij_begin-iip1
        ije=ij_end+iip1
        if (pole_nord) ijb=ij_begin
        if (pole_sud) ije=ij_end
        
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
        DO l = 1, llm
         DO ij = ijb, ije
          tempe(ij) = teta(ij,l) * pk(ij,l) /cpp
         ENDDO
         DO ij = ijb, ije
          zdelta = MAX( 0., SIGN(1., rtt - tempe(ij)) )
          play   = 0.5*(p(ij,l)+p(ij,l+1))
          qsat(ij,l) = MIN(0.5, r2es* FOEEW(tempe(ij),zdelta) / play )
          qsat(ij,l) = qsat(ij,l) / ( 1. - retv * qsat(ij,l) )
         ENDDO
        ENDDO
c$OMP END DO NOWAIT
c      PRINT*,'Debut vlsplt version debug sans vlyqs'

        zzpbar = 0.5 * pdt
        zzw    = pdt

      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ijb+iip1
      if (pole_sud)  ije=ije-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l=1,llm
        DO ij = ijb,ije
            mu(ij,l)=pbaru(ij,l) * zzpbar
         ENDDO
      ENDDO
c$OMP END DO NOWAIT
      
      ijb=ij_begin-iip1
      ije=ij_end
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij=ijb,ije
            mv(ij,l)=pbarv(ij,l) * zzpbar
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

      ijb=ij_begin
      ije=ij_end

      DO iq=1,nqtot
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l=1,llm
         DO ij=ijb,ije
            mw(ij,l,iq)=w(ij,l) * zzw
         ENDDO
      ENDDO
c$OMP END DO NOWAIT
      ENDDO

      DO iq=1,nqtot  
c$OMP MASTER
      DO ij=ijb,ije
         mw(ij,llm+1,iq)=0.
      ENDDO
c$OMP END MASTER
      ENDDO

c      CALL SCOPY(ijp1llm,q,1,zq,1)
c      CALL SCOPY(ijp1llm,masse,1,zm,1)

       ijb=ij_begin
       ije=ij_end

      DO iq=1,nqtot       
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
        DO l=1,llm
          zq(ijb:ije,l,iq)=q(ijb:ije,l,iq)
          zm(ijb:ije,l,iq)=masse(ijb:ije,l)
        ENDDO
c$OMP END DO NOWAIT
      ENDDO

#ifdef DEBUG_IO     
       CALL WriteField_u('mu',mu)
       CALL WriteField_v('mv',mv)
       CALL WriteField_u('mw',mw)
       CALL WriteField_u('qsat',qsat)
#endif

      ! verif temporaire
      ijb=ij_begin
      ije=ij_end  
      if (ok_iso_verif) then
        call check_isotopes(zq,ijb,ije,'vlspltgen_loc 191')
      endif !if (ok_iso_verif) then   

c$OMP BARRIER           
!      DO iq=1,nqtot
      DO iq=1,nqperes ! CRisi: on ne boucle que sur les pères= ceux qui sont transportés directement par l'air
       !write(*,*) 'vlspltgen 192: iq,iadv=',iq,iadv(iq)
#ifdef DEBUG_IO    
       CALL WriteField_u('zq',zq(:,:,iq))
       CALL WriteField_u('zm',zm(:,:,iq)) 
#endif
        if(iadv(iq) == 0) then
        
          cycle 
        
        else if (iadv(iq)==10) then

#ifdef _ADV_HALO        
! CRisi: on ajoute les nombres de fils et tableaux des fils
! On suppose qu'on ne peut advecter les fils que par le schéma 10.  
          call vlx_loc(zq,pente_max,zm,mu,
     &                     ij_begin,ij_begin+2*iip1-1,iq)
          call vlx_loc(zq,pente_max,zm,mu,
     &               ij_end-2*iip1+1,ij_end,iq)
#else
          call vlx_loc(zq,pente_max,zm,mu,
     &                     ij_begin,ij_end,iq)
#endif

c$OMP MASTER
          call VTb(VTHallo)
c$OMP END MASTER
          call Register_Hallo_u(zq(:,:,iq),llm,2,2,2,2,MyRequest1)
          call Register_Hallo_u(zm(:,:,iq),llm,1,1,1,1,MyRequest1)
! CRisi
          do ifils=1,nqdesc(iq)
            iq2=iqfils(ifils,iq)
            call Register_Hallo_u(zq(:,:,iq2),llm,2,2,2,2,MyRequest1)
            call Register_Hallo_u(zm(:,:,iq2),llm,1,1,1,1,MyRequest1)
          enddo

c$OMP MASTER
          call VTe(VTHallo)
c$OMP END MASTER
        else if (iadv(iq)==14) then

#ifdef _ADV_HALO           
          call vlxqs_loc(zq,pente_max,zm,mu,
     &                   qsat,ij_begin,ij_begin+2*iip1-1,iq)
          call vlxqs_loc(zq,pente_max,zm,mu,
     &                   qsat,ij_end-2*iip1+1,ij_end,iq)
#else
          call vlxqs_loc(zq,pente_max,zm,mu,
     &                   qsat,ij_begin,ij_end,iq)
#endif

c$OMP MASTER
          call VTb(VTHallo)
c$OMP END MASTER

          call Register_Hallo_u(zq(:,:,iq),llm,2,2,2,2,MyRequest1)
          call Register_Hallo_u(zm(:,:,iq),llm,1,1,1,1,MyRequest1)
          do ifils=1,nqdesc(iq)
            iq2=iqfils(ifils,iq)
            call Register_Hallo_u(zq(:,:,iq2),llm,2,2,2,2,MyRequest1)
            call Register_Hallo_u(zm(:,:,iq2),llm,1,1,1,1,MyRequest1)
          enddo

c$OMP MASTER
          call VTe(VTHallo)
c$OMP END MASTER 
        else
        
          stop 'vlspltgen_p : schema non parallelise'
      
        endif
      
      enddo !DO iq=1,nqperes
      
      
c$OMP BARRIER      
c$OMP MASTER      
      call VTb(VTHallo)
c$OMP END MASTER

      call SendRequest(MyRequest1)

c$OMP MASTER
      call VTe(VTHallo)
c$OMP END MASTER       
c$OMP BARRIER

      ! verif temporaire
      ijb=ij_begin-2*iip1
      ije=ij_end+2*iip1  
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end  
      if (ok_iso_verif) then
           call check_isotopes(zq,ij_begin,ij_end,'vlspltgen_loc 280')
      endif !if (ok_iso_verif) then

      do iq=1,nqperes
        !write(*,*) 'vlspltgen 279: iq=',iq

        if(iadv(iq) == 0) then
        
          cycle 
        
        else if (iadv(iq)==10) then

#ifdef _ADV_HALLO
          call vlx_loc(zq,pente_max,zm,mu,
     &                 ij_begin+2*iip1,ij_end-2*iip1,iq)
#endif        
        else if (iadv(iq)==14) then
#ifdef _ADV_HALLO
          call vlxqs_loc(zq,pente_max,zm,mu,
     &                    qsat,ij_begin+2*iip1,ij_end-2*iip1,iq)
#endif    
        else
        
          stop 'vlspltgen_p : schema non parallelise'
      
        endif
      
      enddo
c$OMP BARRIER      
c$OMP MASTER
      call VTb(VTHallo)
c$OMP END MASTER

!      call WaitRecvRequest(MyRequest1)
!      call WaitSendRequest(MyRequest1)
c$OMP BARRIER
       call WaitRequest(MyRequest1)


c$OMP MASTER
      call VTe(VTHallo)
c$OMP END MASTER
c$OMP BARRIER

      
      if (ok_iso_verif) then
           call check_isotopes(zq,ij_begin,ij_end,'vlspltgen_loc 326')
      endif !if (ok_iso_verif) then        
      if (ok_iso_verif) then
           ijb=ij_begin-2*iip1
           ije=ij_end+2*iip1
           if (pole_nord) ijb=ij_begin
           if (pole_sud)  ije=ij_end
           call check_isotopes(zq,ijb,ije,'vlspltgen_loc 336')
      endif !if (ok_iso_verif) then  

      do iq=1,nqperes
       !write(*,*) 'vlspltgen 321: iq=',iq
#ifdef DEBUG_IO    
       CALL WriteField_u('zq',zq(:,:,iq))
       CALL WriteField_u('zm',zm(:,:,iq))
#endif

        if(iadv(iq) == 0) then
        
          cycle 
        
        else if (iadv(iq)==10) then
        
          call vly_loc(zq,pente_max,zm,mv,iq)
  
        else if (iadv(iq)==14) then
      
          call vlyqs_loc(zq,pente_max,zm,mv,
     &                   qsat,iq)
 
        else
        
          stop 'vlspltgen_p : schema non parallelise'
      
        endif
       
       enddo

      if (ok_iso_verif) then
           call check_isotopes(zq,ij_begin,ij_end,'vlspltgen_loc 357')
      endif !if (ok_iso_verif) then

      do iq=1,nqperes
      !write(*,*) 'vlspltgen 349: iq=',iq
#ifdef DEBUG_IO    
       CALL WriteField_u('zq',zq(:,:,iq))
       CALL WriteField_u('zm',zm(:,:,iq))
#endif
        if(iadv(iq) == 0) then 
          
          cycle 
        
        else if (iadv(iq)==10 .or. iadv(iq)==14 ) then

c$OMP BARRIER        
#ifdef _ADV_HALLO
          call vlz_loc(zq,pente_max,zm,mw,
     &               ij_begin,ij_begin+2*iip1-1,iq)
          call vlz_loc(zq,pente_max,zm,mw,
     &               ij_end-2*iip1+1,ij_end,iq)
#else
          call vlz_loc(zq,pente_max,zm,mw,
     &               ij_begin,ij_end,iq)
#endif
c$OMP BARRIER

c$OMP MASTER
          call VTb(VTHallo)
c$OMP END MASTER

          call Register_Hallo_u(zq(:,:,iq),llm,2,2,2,2,MyRequest2)
          call Register_Hallo_u(zm(:,:,iq),llm,1,1,1,1,MyRequest2)
          ! CRisi
          do ifils=1,nqdesc(iq)
            iq2=iqfils(ifils,iq)
            call Register_Hallo_u(zq(:,:,iq2),llm,2,2,2,2,MyRequest2)
            call Register_Hallo_u(zm(:,:,iq2),llm,1,1,1,1,MyRequest2)
          enddo      
c$OMP MASTER
          call VTe(VTHallo)
c$OMP END MASTER        
c$OMP BARRIER
        else
        
          stop 'vlspltgen_p : schema non parallelise'
      
        endif
      
      enddo
c$OMP BARRIER      

c$OMP MASTER        
      call VTb(VTHallo)
c$OMP END MASTER

      call SendRequest(MyRequest2)

c$OMP MASTER
      call VTe(VTHallo)
c$OMP END MASTER        


      if (ok_iso_verif) then
           call check_isotopes(zq,ij_begin,ij_end,'vlspltgen_loc 429')
      endif !if (ok_iso_verif) then 

c$OMP BARRIER
      do iq=1,nqperes
      !write(*,*) 'vlspltgen 409: iq=',iq

        if(iadv(iq) == 0) then
          
          cycle 
        
        else if (iadv(iq)==10 .or. iadv(iq)==14 ) then
c$OMP BARRIER        

#ifdef _ADV_HALLO
          call vlz_loc(zq,pente_max,zm,mw,
     &               ij_begin+2*iip1,ij_end-2*iip1,iq)
#endif

c$OMP BARRIER        
        else
        
          stop 'vlspltgen_p : schema non parallelise'
      
        endif
      
      enddo
      !write(*,*) 'vlspltgen_loc 476'

c$OMP BARRIER
      !write(*,*) 'vlspltgen_loc 477'
c$OMP MASTER
      call VTb(VTHallo)
c$OMP END MASTER

!      call WaitRecvRequest(MyRequest2)
!      call WaitSendRequest(MyRequest2)
c$OMP BARRIER
       CALL WaitRequest(MyRequest2)

c$OMP MASTER
      call VTe(VTHallo)
c$OMP END MASTER
c$OMP BARRIER


      !write(*,*) 'vlspltgen_loc 494'
      if (ok_iso_verif) then
           call check_isotopes(zq,ij_begin,ij_end,'vlspltgen_loc 461')
      endif !if (ok_iso_verif) then

      do iq=1,nqperes
      !write(*,*) 'vlspltgen 449: iq=',iq
#ifdef DEBUG_IO    
       CALL WriteField_u('zq',zq(:,:,iq))
       CALL WriteField_u('zm',zm(:,:,iq))
#endif
        if(iadv(iq) == 0) then
        
          cycle 
        
        else if (iadv(iq)==10) then
        
          call vly_loc(zq,pente_max,zm,mv,iq)
  
        else if (iadv(iq)==14) then
      
          call vlyqs_loc(zq,pente_max,zm,mv,
     &                   qsat,iq)
 
        else
        
          stop 'vlspltgen_p : schema non parallelise'
      
        endif
       
       enddo !do iq=1,nqperes

      if (ok_iso_verif) then
           call check_isotopes(zq,ij_begin,ij_end,'vlspltgen_loc 493')
      endif !if (ok_iso_verif) then

      do iq=1,nqperes
      !write(*,*) 'vlspltgen 477: iq=',iq
#ifdef DEBUG_IO    
       CALL WriteField_u('zq',zq(:,:,iq))
       CALL WriteField_u('zm',zm(:,:,iq))
#endif
        if(iadv(iq) == 0) then 
          
          cycle 
        
        else if (iadv(iq)==10) then
        
          call vlx_loc(zq,pente_max,zm,mu,
     &               ij_begin,ij_end,iq)
  
        else if (iadv(iq)==14) then
      
          call vlxqs_loc(zq,pente_max,zm,mu,
     &                 qsat, ij_begin,ij_end,iq)
 
        else
        
          stop 'vlspltgen_p : schema non parallelise'
      
        endif
       
       enddo !do iq=1,nqperes

      !write(*,*) 'vlspltgen 550: apres derniere serie de call vlx'
      if (ok_iso_verif) then
           call check_isotopes(zq,ij_begin,ij_end,'vlspltgen_loc 521')
      endif !if (ok_iso_verif) then
     
      ijb=ij_begin
      ije=ij_end
      !write(*,*) 'vlspltgen_loc 557'
c$OMP BARRIER      

      !write(*,*) 'vlspltgen_loc 559'  
      DO iq=1,nqtot
       !write(*,*) 'vlspltgen_loc 561, iq=',iq  
#ifdef DEBUG_IO    
       CALL WriteField_u('zq',zq(:,:,iq))
       CALL WriteField_u('zm',zm(:,:,iq))
#endif
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)          
        DO l=1,llm
           DO ij=ijb,ije
c             print *,'zq-->',ij,l,iq,zq(ij,l,iq)
c             print *,'q-->',ij,l,iq,q(ij,l,iq)
             q(ij,l,iq)=zq(ij,l,iq)
           ENDDO
        ENDDO
c$OMP END DO NOWAIT    
      !write(*,*) 'vlspltgen_loc 575'      

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
           DO ij=ijb,ije-iip1+1,iip1
              q(ij+iim,l,iq)=q(ij,l,iq)
           ENDDO
        ENDDO
c$OMP END DO NOWAIT  
      !write(*,*) 'vlspltgen_loc 583'  
      ENDDO !DO iq=1,nqtot
        
      if (ok_iso_verif) then
           call check_isotopes(q,ij_begin,ij_end,'vlspltgen_loc 557')
      endif !if (ok_iso_verif) then

c$OMP BARRIER

cc$OMP MASTER      
c      call WaitSendRequest(MyRequest1) 
c      call WaitSendRequest(MyRequest2)
cc$OMP END MASTER
cc$OMP BARRIER

      !write(*,*) 'vlspltgen 597: sortie'  
      RETURN
      END
