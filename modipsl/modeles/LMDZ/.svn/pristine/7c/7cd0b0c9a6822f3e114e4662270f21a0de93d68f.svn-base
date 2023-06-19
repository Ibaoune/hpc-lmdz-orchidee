      SUBROUTINE gradiv2_loc(klevel, xcov, ycov, ld, gdx_out, gdy_out )
c
c     P. Le Van
c
c   **********************************************************
c                                ld
c       calcul  de  (grad (div) )   du vect. v ....
c
c     xcov et ycov etant les composant.covariantes de v
c   **********************************************************
c     xcont , ycont et ld  sont des arguments  d'entree pour le s-prog
c      gdx   et  gdy       sont des arguments de sortie pour le s-prog
c
c
      USE parallel_lmdz
      USE times
      USE Write_field_p
      USE mod_hallo
      USE mod_filtreg_p
      USE gradiv2_mod
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comdissipn.h"
c
c     ........    variables en arguments      ........

      INTEGER klevel
      REAL  xcov( ijb_u:ije_u,klevel ), ycov( ijb_v:ije_v,klevel )
      REAL gdx_out( ijb_u:ije_u,klevel ), gdy_out( ijb_v:ije_v,klevel)
c
c     ........       variables locales       .........
c
      REAL      :: tmp_div2(ijb_u:ije_u,llm)
      REAL signe, nugrads
      INTEGER l,ij,iter,ld
      INTEGER :: ijb,ije,jjb,jje
      Type(Request),SAVE  :: request_dissip
!$OMP THREADPRIVATE(request_dissip)      
c    ........................................................
c
c
c      CALL SCOPY( ip1jmp1 * klevel, xcov, 1, gdx, 1 )
c      CALL SCOPY(   ip1jm * klevel, ycov, 1, gdy, 1 )
      
      ijb=ij_begin
      ije=ij_end
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO   l = 1, klevel
        gdx(ijb:ije,l)=xcov(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT      
      
      ijb=ij_begin
      ije=ij_end
      if(pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO   l = 1, klevel
        gdy(ijb:ije,l)=ycov(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT

c$OMP BARRIER
       call Register_Hallo_v(gdy,llm,1,0,0,1,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)
c$OMP BARRIER
c
c
      signe   = (-1.)**ld
      nugrads = signe * cdivu
c


      CALL    divergf_loc( klevel, gdx,   gdy , div )
c      call write_field3d_p('div1',reshape(div,(/iip1,jjp1,llm/)))

      IF( ld.GT.1 )   THEN
c$OMP BARRIER
       call Register_Hallo_u(div,llm,1,1,1,1,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)
c$OMP BARRIER
	CALL laplacien_loc( klevel, div,  div     )

c    ......  Iteration de l'operateur laplacien_gam   .......
c         call write_field3d_p('div2',reshape(div,(/iip1,jjp1,llm/)))

        DO iter = 1, ld -2
c$OMP BARRIER
       call Register_Hallo_u(div,llm,1,1,1,1,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)

c$OMP BARRIER

         CALL laplacien_gam_loc(klevel,cuvscvgam1,cvuscugam1,
     &                          unsair_gam1,unsapolnga1, unsapolsga1,
     &                          div, div       )
        ENDDO
c        call write_field3d_p('div3',reshape(div,(/iip1,jjp1,llm/)))
      ENDIF

       jjb=jj_begin
       jje=jj_end
       
       CALL filtreg_p( div   ,jjb_u,jje_u,jjb,jje, jjp1, 
     &                 klevel, 2, 1, .TRUE., 1 )
c       call exchange_Hallo(div,ip1jmp1,llm,0,1)
c$OMP BARRIER
       call Register_Hallo_u(div,llm,1,1,1,1,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)

c$OMP BARRIER


       CALL  grad_loc( klevel,  div,   gdx,  gdy )

c
      ijb=ij_begin
      ije=ij_end
         
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       DO   l = 1, klevel
         
         if (pole_sud) ije=ij_end
         DO  ij = ijb, ije
          gdx_out( ij,l ) = gdx( ij,l ) * nugrads
         ENDDO
         
         if (pole_sud) ije=ij_end-iip1
         DO  ij = ijb, ije
          gdy_out( ij,l ) = gdy( ij,l ) * nugrads
         ENDDO
       
       ENDDO
c$OMP END DO NOWAIT
c
       RETURN
       END
