!
! $Id$
!
c
c
#define DEBUG_IO
#undef DEBUG_IO
      SUBROUTINE advtrac_loc(pbarug,pbarvg ,wg,
     *                   p,  massem,q,teta,
     *                   pk   )

c     Auteur :  F. Hourdin
c
c     Modif. P. Le Van     (20/12/97)
c            F. Codron     (10/99)
c            D. Le Croller (07/2001)
c            M.A Filiberti (04/2002)
c
      USE parallel_lmdz
      USE Write_Field_loc
      USE Write_Field
      USE Bands
      USE mod_hallo
      USE Vampir
      USE times
      USE infotrac, ONLY: nqtot, iadv, ok_iso_verif
      USE control_mod, ONLY: iapp_tracvl, day_step, planet_type
      USE advtrac_mod, ONLY: finmasse
      USE comconst_mod, ONLY: dtvr
      IMPLICIT NONE
c
      include "dimensions.h"
      include "paramet.h"
      include "comdissip.h"
      include "comgeom2.h"
      include "description.h"

c-------------------------------------------------------------------
c     Arguments
c-------------------------------------------------------------------
c     Ajout PPM
c--------------------------------------------------------
      REAL massebx(ijb_u:ije_u,llm),masseby(ijb_v:ije_v,llm)
c--------------------------------------------------------
      INTEGER iapptrac
      REAL pbarug(ijb_u:ije_u,llm),pbarvg(ijb_v:ije_v,llm)
      REAL wg(ijb_u:ije_u,llm)
      REAL q(ijb_u:ije_u,llm,nqtot),massem(ijb_u:ije_u,llm)
      REAL p( ijb_u:ije_u,llmp1 ),teta(ijb_u:ije_u,llm)
      REAL pk(ijb_u:ije_u,llm)

c-------------------------------------------------------------
c     Variables locales
c-------------------------------------------------------------

      REAL zdp(ijb_u:ije_u)
      REAL (kind=kind(1.d0)) :: t_initial, t_final, tps_cpu
      INTEGER,SAVE :: iadvtr=0
c$OMP THREADPRIVATE(iadvtr)
      INTEGER ij,l,iq,iiq
      REAL zdpmin, zdpmax
c----------------------------------------------------------
c     Rajouts pour PPM
c----------------------------------------------------------
      INTEGER indice,n
      REAL dtbon ! Pas de temps adaptatif pour que CFL<1
      REAL CFLmaxz,aaa,bbb ! CFL maximum
      REAL psppm(iim,jjb_u:jje_u) ! pression  au sol
      REAL unatppm(iim,jjb_u:jje_u,llm),vnatppm(iim,jjb_u:jje_u,llm)
      REAL qppm(iim*jjnb_u,llm,nqtot)
      REAL fluxwppm(iim,jjb_u:jje_u,llm)
      REAL apppm(llmp1), bpppm(llmp1)
      LOGICAL dum,fill
      DATA fill/.true./
      DATA dum/.true./
      integer ijb,ije,ijbu,ijbv,ijeu,ijev,j
      type(Request),SAVE :: testRequest
!$OMP THREADPRIVATE(testRequest)

c  test sur l''eventuelle creation de valeurs negatives de la masse
         ijb=ij_begin
         ije=ij_end
         if (pole_nord) ijb=ij_begin+iip1
         if (pole_sud) ije=ij_end-iip1
         
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)          
         DO l=1,llm-1
            DO ij = ijb+1,ije
              zdp(ij) =    pbarug(ij-1,l)   - pbarug(ij,l)
     s                  - pbarvg(ij-iip1,l) + pbarvg(ij,l)
     s                  +       wg(ij,l+1)  - wg(ij,l)
            ENDDO
            
c            CALL SCOPY( jjm -1 ,zdp(iip1+iip1),iip1,zdp(iip2),iip1 )
c ym  ---> pourquoi jjm-1 et non jjm ? a cause du pole ?
            
            do ij=ijb,ije-iip1+1,iip1
              zdp(ij)=zdp(ij+iip1-1)
            enddo
            
            DO ij = ijb,ije
               zdp(ij)= zdp(ij)*dtvr/ massem(ij,l) 
            ENDDO 


c            CALL minmax ( ip1jm-iip1, zdp(iip2), zdpmin,zdpmax )
c  ym ---> eventuellement a revoir
            CALL minmax ( ije-ijb+1, zdp(ijb), zdpmin,zdpmax )
            
            IF(MAX(ABS(zdpmin),ABS(zdpmax)).GT.0.5) THEN
            PRINT*,'WARNING DP/P l=',l,'  MIN:',zdpmin,
     s        '   MAX:', zdpmax
            ENDIF

         ENDDO
c$OMP END DO NOWAIT

c-------------------------------------------------------------------
c   Advection proprement dite (Modification Le Croller (07/2001)
c-------------------------------------------------------------------

c----------------------------------------------------
c        Calcul des moyennes basées sur la masse
c----------------------------------------------------

cym      ----> Normalement, inutile pour les schémas classiques
cym      ----> Revérifier lors de la parallélisation des autres schemas
   
cym          call massbar_p(massem,massebx,masseby)          

#ifdef DEBUG_IO    
          CALL WriteField_u('massem',massem)
          CALL WriteField_u('wg',wg)
          CALL WriteField_u('pbarug',pbarug)
          CALL WriteField_v('pbarvg',pbarvg)
          CALL WriteField_u('p_tmp',p)
          CALL WriteField_u('pk_tmp',pk)
          CALL WriteField_u('teta_tmp',teta)
          do j=1,nqtot
            call WriteField_u('q_adv'//trim(int2str(j)),
     .                q(:,:,j))
          enddo
#endif

!          
!       call Register_Hallo_v(pbarvg,llm,1,1,1,1,TestRequest)
!
!       call SendRequest(TestRequest)
!c$OMP BARRIER
!       call WaitRequest(TestRequest)
c$OMP BARRIER
                 
          !write(*,*) 'advtrac 157: appel de vlspltgen_loc'
          call vlspltgen_loc( q,iadv, 2., massem, wg ,
     *                        pbarug,pbarvg,dtvr,p,
     *                        pk,teta )

          !write(*,*) 'advtrac 162: apres appel vlspltgen_loc'
      if (ok_iso_verif) then
           call check_isotopes(q,ijb_u,ije_u,'advtrac 162')
      endif !if (ok_iso_verif) then

#ifdef DEBUG_IO     
          do j=1,nqtot
            call WriteField_u('q_adv'//trim(int2str(j)),
     .                q(:,:,j))
          enddo
#endif
         
          GOTO 1234     
c-----------------------------------------------------------
c     Appel des sous programmes d'advection
c-----------------------------------------------------------
      do iq=1,nqtot
c        call clock(t_initial)
        if(iadv(iq) == 0) cycle 
c   ----------------------------------------------------------------
c   Schema de Van Leer I MUSCL
c   ----------------------------------------------------------------
        if(iadv(iq).eq.10) THEN
     
!LF            call vlsplt_p(q(1,1,iq),2.,massem,wg,pbarug,pbarvg,dtvr)

c   ----------------------------------------------------------------
c   Schema "pseudo amont" + test sur humidite specifique
C    pour la vapeur d'eau. F. Codron
c   ----------------------------------------------------------------
        else if(iadv(iq).eq.14) then
c
cym           stop 'advtrac : appel à vlspltqs :schema non parallelise'
!LF           CALL vlspltqs_p( q(1,1,1), 2., massem, wg ,
!LF     *                 pbarug,pbarvg,dtvr,p,pk,teta )
c   ----------------------------------------------------------------
c   Schema de Frederic Hourdin
c   ----------------------------------------------------------------
        else if(iadv(iq).eq.12) then
          stop 'advtrac : schema non parallelise'
c            Pas de temps adaptatif
           call adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           if (n.GT.1) then
           write(*,*) 'WARNING horizontal dt=',dtbon,'dtvr=',
     s             dtvr,'n=',n
           endif
           do indice=1,n
            call advn(q(1,1,iq),massem,wg,pbarug,pbarvg,dtbon,1)
           end do
        else if(iadv(iq).eq.13) then
          stop 'advtrac : schema non parallelise'
c            Pas de temps adaptatif
           call adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           if (n.GT.1) then
           write(*,*) 'WARNING horizontal dt=',dtbon,'dtvr=',
     s             dtvr,'n=',n
           endif
          do indice=1,n
            call advn(q(1,1,iq),massem,wg,pbarug,pbarvg,dtbon,2)
          end do
c   ----------------------------------------------------------------
c   Schema de pente SLOPES
c   ----------------------------------------------------------------
        else if (iadv(iq).eq.20) then
          stop 'advtrac : schema non parallelise'

            call pentes_ini (q(1,1,iq),wg,massem,pbarug,pbarvg,0)

c   ----------------------------------------------------------------
c   Schema de Prather
c   ----------------------------------------------------------------
        else if (iadv(iq).eq.30) then
          stop 'advtrac : schema non parallelise'
c            Pas de temps adaptatif
           call adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           if (n.GT.1) then
           write(*,*) 'WARNING horizontal dt=',dtbon,'dtvr=',
     s             dtvr,'n=',n
           endif
           call  prather(q(1,1,iq),wg,massem,pbarug,pbarvg,
     s                     n,dtbon)
c   ----------------------------------------------------------------
c   Schemas PPM Lin et Rood
c   ----------------------------------------------------------------
       else if (iadv(iq).eq.11.OR.(iadv(iq).GE.16.AND.
     s                     iadv(iq).LE.18)) then

           stop 'advtrac : schema non parallelise'

c        Test sur le flux horizontal
c        Pas de temps adaptatif
         call adaptdt(iadv(iq),dtbon,n,pbarug,massem)
         if (n.GT.1) then
           write(*,*) 'WARNING horizontal dt=',dtbon,'dtvr=',
     s             dtvr,'n=',n
         endif
c        Test sur le flux vertical
         CFLmaxz=0.
         do l=2,llm
           do ij=iip2,ip1jm
            aaa=wg(ij,l)*dtvr/massem(ij,l)
            CFLmaxz=max(CFLmaxz,aaa)
            bbb=-wg(ij,l)*dtvr/massem(ij,l-1)
            CFLmaxz=max(CFLmaxz,bbb)
           enddo
         enddo
         if (CFLmaxz.GE.1) then
            write(*,*) 'WARNING vertical','CFLmaxz=', CFLmaxz
         endif

c-----------------------------------------------------------
c        Ss-prg interface LMDZ.4->PPM3d
c-----------------------------------------------------------

          call interpre(q(1,1,iq),qppm(1,1,iq),wg,fluxwppm,massem,
     s                 apppm,bpppm,massebx,masseby,pbarug,pbarvg,
     s                 unatppm,vnatppm,psppm)

          do indice=1,n
c---------------------------------------------------------------------
c                         VL (version PPM) horiz. et PPM vert.
c---------------------------------------------------------------------
                if (iadv(iq).eq.11) then
c                  Ss-prg PPM3d de Lin
                  call ppm3d(1,qppm(1,1,iq),
     s                       psppm,psppm,
     s                       unatppm,vnatppm,fluxwppm,dtbon,2,2,2,1,
     s                       iim,jjp1,2,llm,apppm,bpppm,0.01,6400000,
     s                       fill,dum,220.)

c----------------------------------------------------------------------
c                           Monotonic PPM
c----------------------------------------------------------------------
               else if (iadv(iq).eq.16) then
c                  Ss-prg PPM3d de Lin
                  call ppm3d(1,qppm(1,1,iq),
     s                       psppm,psppm,
     s                       unatppm,vnatppm,fluxwppm,dtbon,3,3,3,1,
     s                       iim,jjp1,2,llm,apppm,bpppm,0.01,6400000,
     s                       fill,dum,220.)
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c                           Semi Monotonic PPM
c---------------------------------------------------------------------
               else if (iadv(iq).eq.17) then
c                  Ss-prg PPM3d de Lin
                  call ppm3d(1,qppm(1,1,iq),
     s                       psppm,psppm,
     s                       unatppm,vnatppm,fluxwppm,dtbon,4,4,4,1,
     s                       iim,jjp1,2,llm,apppm,bpppm,0.01,6400000,
     s                       fill,dum,220.)
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c                         Positive Definite PPM
c---------------------------------------------------------------------
                else if (iadv(iq).eq.18) then
c                  Ss-prg PPM3d de Lin
                  call ppm3d(1,qppm(1,1,iq),
     s                       psppm,psppm,
     s                       unatppm,vnatppm,fluxwppm,dtbon,5,5,5,1,
     s                       iim,jjp1,2,llm,apppm,bpppm,0.01,6400000,
     s                       fill,dum,220.)
c---------------------------------------------------------------------
                endif
            enddo
c-----------------------------------------------------------------
c               Ss-prg interface PPM3d-LMDZ.4
c-----------------------------------------------------------------
                  call interpost(q(1,1,iq),qppm(1,1,iq))
            endif
c----------------------------------------------------------------------

c-----------------------------------------------------------------
c On impose une seule valeur du traceur au pôle Sud j=jjm+1=jjp1
c et Nord j=1
c-----------------------------------------------------------------

c                  call traceurpole(q(1,1,iq),massem)

c calcul du temps cpu pour un schema donne


       end DO

1234  CONTINUE
c$OMP BARRIER

      if (planet_type=="earth") then

      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       DO l = 1, llm
         DO ij = ijb, ije
           finmasse(ij,l) =  p(ij,l) - p(ij,l+1) 
         ENDDO
       ENDDO
c$OMP END DO

        ! CRisi: on passe nqtot et non nq
       CALL qminimum_loc( q, nqtot, finmasse )

      endif ! of if (planet_type=="earth")

       RETURN
       END

