!
! $Id: caladvtrac_p.F 1299 2010-01-20 14:27:21Z fairhead $
!
c
c
            SUBROUTINE caladvtrac_loc(q,pbaru,pbarv ,
     *                   p ,masse, dq ,  teta,
     *                   flxw, pk, iapptrac)
      USE parallel_lmdz 
      USE infotrac, ONLY : nqtot
      USE control_mod, ONLY : iapp_tracvl,planet_type
      USE caladvtrac_mod
      USE mod_hallo
      USE bands
      USE times
      USE Vampir
      USE write_field_loc
c
      IMPLICIT NONE
c
c     Auteurs:   F.Hourdin , P.Le Van, F.Forget, F.Codron  
c
c     F.Codron (10/99) : ajout humidite specifique pour eau vapeur
c=======================================================================
c
c       Shema de  Van Leer
c
c=======================================================================


      include "dimensions.h"
      include "paramet.h"

c   Arguments:
c   ----------
      REAL :: pbaru( ijb_u:ije_u,llm ),pbarv( ijb_v:ije_v,llm)
      REAL :: masse(ijb_u:ije_u,llm)
      REAL :: p( ijb_u:ije_u,llmp1)
      REAL :: q( ijb_u:ije_u,llm,nqtot),dq( ijb_u:ije_u,llm, nqtot )
      REAL :: teta( ijb_u:ije_u,llm),pk( ijb_u:ije_u,llm)
      REAL :: flxw(ijb_u:ije_u,llm)
      INTEGER :: iapptrac
c   Local:
c   ------
!      REAL :: pbarug(ijb_u:ije_u,llm)
!      REAL :: pbarvg(ijb_v:ije_v,llm)
!      REAL :: wg(ijb_u:ije_u,llm)
      
      REAL :: flxw_adv(distrib_vanleer%ijb_u:distrib_vanleer%ije_u,llm)
      INTEGER,SAVE :: iadvtr=0
!$OMP THREADPRIVATE(iadvtr)
      INTEGER ::  ijb,ije,ijbu,ijbv,ijeu,ijev,j
      INTEGER :: ij,l
      TYPE(Request),SAVE :: Request_vanleer
!$OMP THREADPRIVATE(Request_vanleer)

      !write(*,*) 'caladvtrac 58: entree'     
      ijbu=ij_begin
      ijeu=ij_end
      
      ijbv=ij_begin-iip1
      ijev=ij_end
      if (pole_nord) ijbv=ij_begin
      if (pole_sud)  ijev=ij_end-iip1

      IF(iadvtr.EQ.0) THEN
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)        
	DO l=1,llm   
          pbaruc(ijbu:ijeu,l)=0.
          pbarvc(ijbv:ijev,l)=0.
        ENDDO
c$OMP END DO NOWAIT  
      ENDIF

c   accumulation des flux de masse horizontaux
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij = ijbu,ijeu
            pbaruc(ij,l) = pbaruc(ij,l) + pbaru(ij,l)
         ENDDO
         DO ij = ijbv,ijev
            pbarvc(ij,l) = pbarvc(ij,l) + pbarv(ij,l)
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

c   selection de la masse instantannee des mailles avant le transport.
      IF(iadvtr.EQ.0) THEN

          ijb=ij_begin
          ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
       DO l=1,llm
          massem(ijb:ije,l)=masse(ijb:ije,l)
       ENDDO
c$OMP END DO NOWAIT

      ENDIF

      iadvtr   = iadvtr+1

c$OMP MASTER
      iapptrac = iadvtr
c$OMP END MASTER

c   Test pour savoir si on advecte a ce pas de temps

      IF ( iadvtr.EQ.iapp_tracvl ) THEN
      !write(*,*) 'caladvtrac 133'
c$OMP MASTER
        call suspend_timer(timer_caldyn)
c$OMP END MASTER
      
      ijb=ij_begin
      ije=ij_end
      
cc   ..  Modif P.Le Van  ( 20/12/97 )  ....
cc

c   traitement des flux de masse avant advection.
c     1. calcul de w
c     2. groupement des mailles pres du pole.

        CALL groupe_loc( massem, pbaruc,pbarvc, pbarug,pbarvg,wg )

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l=1,llm
        flxw(ijb:ije,l)=wg(ijb:ije,l)/REAL(iapp_tracvl)
      ENDDO
c$OMP ENDDO NOWAIT

#ifdef DEBUG_IO    
         CALL WriteField_u('pbarug1',pbarug)
         CALL WriteField_v('pbarvg1',pbarvg)
         CALL WriteField_u('wg1',wg)
#endif

c$OMP BARRIER


c$OMP MASTER
      call VTb(VTHallo)
c$OMP END MASTER

      call Register_SwapField_u(pbarug,pbarug_adv, distrib_vanleer,
     &                          Request_vanleer)
      call Register_SwapField_v(pbarvg,pbarvg_adv, distrib_vanleer,
     &                          Request_vanleer,up=1)
      call Register_SwapField_u(massem,massem_adv, distrib_vanleer,
     &                          Request_vanleer)
      call Register_SwapField_u(wg,wg_adv,distrib_vanleer,
     &                          Request_vanleer)
      call Register_SwapField_u(teta,teta_adv, distrib_vanleer,
     &                          Request_vanleer,up=1,down=1)
      call Register_SwapField_u(p,p_adv, distrib_vanleer,
     &                          Request_vanleer,up=1,down=1)
      call Register_SwapField_u(pk,pk_adv, distrib_vanleer,
     &                          Request_vanleer,up=1,down=1)
      call Register_SwapField_u(q,q_adv, distrib_vanleer,
     &                          Request_vanleer)

      call SendRequest(Request_vanleer)
c$OMP BARRIER
      call WaitRequest(Request_vanleer)


c$OMP BARRIER
c$OMP MASTER      
      call Set_Distrib(distrib_vanleer)
      call VTe(VTHallo)
      call VTb(VTadvection)
      call start_timer(timer_vanleer)
c$OMP END MASTER
c$OMP BARRIER
!      CALL WriteField_u('pbarug_adv',pbarug_adv)
!      CALL WriteField_u('',)
      
      
#ifdef DEBUG_IO
         CALL WriteField_u('pbarug1',pbarug_adv)
         CALL WriteField_v('pbarvg1',pbarvg_adv)
         CALL WriteField_u('wg1',wg_adv)
#endif        
      !write(*,*) 'caladvtrac 185'  
      CALL advtrac_loc( pbarug_adv,pbarvg_adv,wg_adv, 
     *             p_adv,  massem_adv,q_adv, teta_adv,
     .             pk_adv)     
      !write(*,*) 'caladvtrac 189'


c$OMP MASTER
	call VTe(VTadvection)
        call stop_timer(timer_vanleer)
        call VTb(VThallo)
c$OMP END MASTER

        call Register_SwapField_u(q_adv,q,distrib_caldyn,
     *                             Request_vanleer)

        call SendRequest(Request_vanleer)
c$OMP BARRIER
        call WaitRequest(Request_vanleer)      

c$OMP BARRIER
c$OMP MASTER
        call Set_Distrib(distrib_caldyn)
	call VTe(VThallo)
	call resume_timer(timer_caldyn)
c$OMP END MASTER
c$OMP BARRIER
          iadvtr=0
       ENDIF ! if iadvtr.EQ.iapp_tracvl

      END


