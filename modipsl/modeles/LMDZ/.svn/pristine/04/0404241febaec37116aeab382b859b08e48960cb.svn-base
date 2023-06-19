!
! $Id$
!
      RECURSIVE SUBROUTINE vlx_loc(q,pente_max,masse,u_m,ijb_x,ije_x,iq)

c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c
c
c   --------------------------------------------------------------------
      USE parallel_lmdz
      USE infotrac, ONLY : nqtot,nqfils,nqdesc,iqfils ! CRisi
      IMPLICIT NONE
c
      include "dimensions.h"
      include "paramet.h"
      include "iniprint.h"
c
c
c   Arguments:
c   ----------
      REAL masse(ijb_u:ije_u,llm,nqtot),pente_max
      REAL u_m( ijb_u:ije_u,llm),pbarv( iip1,jjb_v:jje_v,llm)
      REAL q(ijb_u:ije_u,llm,nqtot) ! CRisi: ajout dimension nqtot
      REAL w(ijb_u:ije_u,llm) 
      INTEGER iq ! CRisi
c
c      Local 
c   ---------
c
      INTEGER ij,l,j,i,iju,ijq,indu(ijnb_u),niju
      INTEGER n0,iadvplus(ijb_u:ije_u,llm),nl(llm)
c
      REAL new_m,zu_m,zdum(ijb_u:ije_u,llm)
      REAL sigu(ijb_u:ije_u),dxq(ijb_u:ije_u,llm),dxqu(ijb_u:ije_u)
      REAL zz(ijb_u:ije_u)
      REAL adxqu(ijb_u:ije_u),dxqmax(ijb_u:ije_u,llm)
      REAL u_mq(ijb_u:ije_u,llm)

      REAL Ratio(ijb_u:ije_u,llm,nqtot) ! CRisi
      INTEGER ifils,iq2 ! CRisi

      Logical extremum

      REAL      SSUM
      EXTERNAL  SSUM

      REAL z1,z2,z3

      INTEGER ijb,ije,ijb_x,ije_x
      
      !write(*,*) 'vlsplt 58: entree dans vlx_loc, iq,ijb_x=',
!     &   iq,ijb_x
c   calcul de la pente a droite et a gauche de la maille

      ijb=ijb_x
      ije=ije_x
        
      if (pole_nord.and.ijb==1) ijb=ijb+iip1
      if (pole_sud.and.ije==ip1jmp1)  ije=ije-iip1
         
      IF (pente_max.gt.-1.e-5) THEN
c       IF (pente_max.gt.10) THEN

c   calcul des pentes avec limitation, Van Leer scheme I:
c   -----------------------------------------------------
      ! on a besoin de q entre ijb et ije
c   calcul de la pente aux points u
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)         
         DO l = 1, llm
            
            DO ij=ijb,ije-1
               dxqu(ij)=q(ij+1,l,iq)-q(ij,l,iq)
c              IF(u_m(ij,l).lt.0.) stop'limx n admet pas les U<0'
c              sigu(ij)=u_m(ij,l)/masse(ij,l,iq)
            ENDDO
            DO ij=ijb+iip1-1,ije,iip1
               dxqu(ij)=dxqu(ij-iim)
c              sigu(ij)=sigu(ij-iim)
            ENDDO

            DO ij=ijb,ije
               adxqu(ij)=abs(dxqu(ij))
            ENDDO

c   calcul de la pente maximum dans la maille en valeur absolue

            DO ij=ijb+1,ije
               dxqmax(ij,l)=pente_max*
     ,      min(adxqu(ij-1),adxqu(ij))
c limitation subtile
c    ,      min(adxqu(ij-1)/sigu(ij-1),adxqu(ij)/(1.-sigu(ij)))
          

            ENDDO

            DO ij=ijb+iip1-1,ije,iip1
               dxqmax(ij-iim,l)=dxqmax(ij,l)
            ENDDO

            DO ij=ijb+1,ije
#ifdef CRAY
               dxq(ij,l)=
     ,         cvmgp(dxqu(ij-1)+dxqu(ij),0.,dxqu(ij-1)*dxqu(ij))
#else
               IF(dxqu(ij-1)*dxqu(ij).gt.0) THEN
                  dxq(ij,l)=dxqu(ij-1)+dxqu(ij)
               ELSE
c   extremum local
                  dxq(ij,l)=0.
               ENDIF
#endif
               dxq(ij,l)=0.5*dxq(ij,l)
               dxq(ij,l)=
     ,         sign(min(abs(dxq(ij,l)),dxqmax(ij,l)),dxq(ij,l))
            ENDDO

         ENDDO ! l=1,llm
c$OMP END DO NOWAIT
c	print*,'Ok calcul des pentes'

      ELSE ! (pente_max.lt.-1.e-5)

c   Pentes produits:
c   ----------------
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l = 1, llm
            DO ij=ijb,ije-1
               dxqu(ij)=q(ij+1,l,iq)-q(ij,l,iq)
            ENDDO
            DO ij=ijb+iip1-1,ije,iip1
               dxqu(ij)=dxqu(ij-iim)
            ENDDO

            DO ij=ijb+1,ije
               zz(ij)=dxqu(ij-1)*dxqu(ij)
               zz(ij)=zz(ij)+zz(ij)
               IF(zz(ij).gt.0) THEN
                  dxq(ij,l)=zz(ij)/(dxqu(ij-1)+dxqu(ij))
               ELSE
c   extremum local
                  dxq(ij,l)=0.
               ENDIF
            ENDDO

         ENDDO
c$OMP END DO NOWAIT
      ENDIF ! (pente_max.lt.-1.e-5)

      !write(*,*) 'vlx 156: iq,ijb_x=',iq,ijb_x

c   bouclage de la pente en iip1:
c   -----------------------------
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij=ijb+iip1-1,ije,iip1
            dxq(ij-iim,l)=dxq(ij,l)
         ENDDO
         DO ij=ijb,ije
            iadvplus(ij,l)=0
         ENDDO

      ENDDO
c$OMP END DO NOWAIT
c	 print*,'Bouclage en iip1'

c   calcul des flux a gauche et a droite

#ifdef CRAY
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
       DO ij=ijb,ije-1
          zdum(ij,l)=cvmgp(1.-u_m(ij,l)/masse(ij,l,iq),
     ,                     1.+u_m(ij,l)/masse(ij+1,l,iq),
     ,                     u_m(ij,l,iq))
          zdum(ij,l)=0.5*zdum(ij,l)
          u_mq(ij,l)=cvmgp(
     ,                q(ij,l,iq)+zdum(ij,l)*dxq(ij,l),
     ,                q(ij+1,l,iq)-zdum(ij,l)*dxq(ij+1,l),
     ,                u_m(ij,l))
          u_mq(ij,l)=u_m(ij,l)*u_mq(ij,l)
       ENDDO
      ENDDO
c$OMP END DO NOWAIT
#else
c   on cumule le flux correspondant a toutes les mailles dont la masse
c   au travers de la paroi pENDant le pas de temps.
c	print*,'Cumule ....'
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        ! on a besoin de masse entre ijb et ije
      DO l=1,llm
       DO ij=ijb,ije-1
c	print*,'masse(',ij,')=',masse(ij,l,iq)
          IF (u_m(ij,l).gt.0.) THEN
             zdum(ij,l)=1.-u_m(ij,l)/masse(ij,l,iq)
             u_mq(ij,l)=u_m(ij,l)*(q(ij,l,iq)
     :           +0.5*zdum(ij,l)*dxq(ij,l))
          ELSE
             zdum(ij,l)=1.+u_m(ij,l)/masse(ij+1,l,iq)
             u_mq(ij,l)=u_m(ij,l)*(q(ij+1,l,iq)
     :           -0.5*zdum(ij,l)*dxq(ij+1,l))
          ENDIF
       ENDDO
      ENDDO
c$OMP END DO NOWAIT
#endif
c	stop

c	go to 9999
c   detection des points ou on advecte plus que la masse de la
c   maille
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij=ijb,ije-1
            IF(zdum(ij,l).lt.0) THEN
               iadvplus(ij,l)=1
               u_mq(ij,l)=0.
            ENDIF
         ENDDO
      ENDDO
c$OMP END DO NOWAIT
c	print*,'Ok test 1'

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
       DO ij=ijb+iip1-1,ije,iip1
          iadvplus(ij,l)=iadvplus(ij-iim,l)
       ENDDO
      ENDDO
c$OMP END DO NOWAIT
c	 print*,'Ok test 2'
        

c   traitement special pour le cas ou on advecte en longitude plus que le
c   contenu de la maille.
c   cette partie est mal vectorisee.

c  calcul du nombre de maille sur lequel on advecte plus que la maille.

      n0=0
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         nl(l)=0
         DO ij=ijb,ije
            nl(l)=nl(l)+iadvplus(ij,l)
         ENDDO
         n0=n0+nl(l)
      ENDDO
c$OMP END DO NOWAIT
cym      IF(n0.gt.1) THEN
cym      IF(n0.gt.0) THEN

c      PRINT*,'Nombre de points pour lesquels on advect plus que le'
c     &       ,'contenu de la maille : ',n0
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)


         DO l=1,llm
            IF(nl(l).gt.0) THEN
               iju=0
c   indicage des mailles concernees par le traitement special
               DO ij=ijb,ije
                  IF(iadvplus(ij,l).eq.1.and.mod(ij,iip1).ne.0) THEN
                     iju=iju+1
                     indu(iju)=ij
                  ENDIF
               ENDDO
               niju=iju
               !PRINT*,'vlx 278, niju,nl',niju,nl(l)

c  traitement des mailles
               DO iju=1,niju
                  ij=indu(iju)
                  j=(ij-1)/iip1+1
                  zu_m=u_m(ij,l)
                  u_mq(ij,l)=0.
                  IF(zu_m.gt.0.) THEN
                     ijq=ij
                     i=ijq-(j-1)*iip1
c   accumulation pour les mailles completements advectees
                     do while(zu_m.gt.masse(ijq,l,iq))
                        u_mq(ij,l)=u_mq(ij,l)
     &                          +q(ijq,l,iq)*masse(ijq,l,iq)
                        zu_m=zu_m-masse(ijq,l,iq)
                        i=mod(i-2+iim,iim)+1
                        ijq=(j-1)*iip1+i
                     ENDDO
c   ajout de la maille non completement advectee
                     u_mq(ij,l)=u_mq(ij,l)+zu_m*
     &               (q(ijq,l,iq)+0.5*
     &               (1.-zu_m/masse(ijq,l,iq))*dxq(ijq,l))
                  ELSE
                     ijq=ij+1
                     i=ijq-(j-1)*iip1
c   accumulation pour les mailles completements advectees
                     do while(-zu_m.gt.masse(ijq,l,iq))
                        u_mq(ij,l)=u_mq(ij,l)-q(ijq,l,iq)
     &                           *masse(ijq,l,iq)
                        zu_m=zu_m+masse(ijq,l,iq)
                        i=mod(i,iim)+1
                        ijq=(j-1)*iip1+i
                     ENDDO
c   ajout de la maille non completement advectee
                     u_mq(ij,l)=u_mq(ij,l)+zu_m*(q(ijq,l,iq)-
     &               0.5*(1.+zu_m/masse(ijq,l,iq))*dxq(ijq,l))
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
c$OMP END DO NOWAIT
cym      ENDIF  ! n0.gt.0 
9999    continue

c   bouclage en latitude
c	print*,'Avant bouclage en latitude'
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
        DO ij=ijb+iip1-1,ije,iip1
           u_mq(ij,l)=u_mq(ij-iim,l)
        ENDDO
      ENDDO
c$OMP END DO NOWAIT

! CRisi: appel récursif de l'advection sur les fils.
! Il faut faire ça avant d'avoir mis à jour q et masse

      !write(*,*) 'vlsplt 326: iq,ijb_x,nqfils(iq)=',iq,ijb_x,nqfils(iq) 

      if (nqfils(iq).gt.0) then  
       do ifils=1,nqdesc(iq)
         iq2=iqfils(ifils,iq)
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l=1,llm
          DO ij=ijb,ije
           ! On a besoin de q et masse seulement entre ijb et ije. On ne
           ! les calcule donc que de ijb à ije
           masse(ij,l,iq2)=masse(ij,l,iq)*q(ij,l,iq)
           Ratio(ij,l,iq2)=q(ij,l,iq2)/q(ij,l,iq)
          enddo   
         enddo
c$OMP END DO NOWAIT
        enddo !do ifils=1,nqdesc(iq)
        do ifils=1,nqfils(iq)
         iq2=iqfils(ifils,iq)
         call vlx_loc(Ratio,pente_max,masse,u_mq,ijb_x,ije_x,iq2)
        enddo !do ifils=1,nqfils(iq)
      endif !if (nqfils(iq).gt.0) then
! end CRisi

      !write(*,*) 'vlsplt 360: iq,ijb_x=',iq,ijb_x

c   calcul des tENDances
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij=ijb+1,ije
            new_m=masse(ij,l,iq)+u_m(ij-1,l)-u_m(ij,l)
            q(ij,l,iq)=(q(ij,l,iq)*masse(ij,l,iq)+
     &        u_mq(ij-1,l)-u_mq(ij,l))
     &        /new_m
            masse(ij,l,iq)=new_m
         ENDDO
c   ModIF Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
         DO ij=ijb+iip1-1,ije,iip1
            q(ij-iim,l,iq)=q(ij,l,iq)
            masse(ij-iim,l,iq)=masse(ij,l,iq)
         ENDDO
      ENDDO
c$OMP END DO NOWAIT
      !write(*,*) 'vlsplt 380: iq,ijb_x=',iq,ijb_x

! retablir les fils en rapport de melange par rapport a l'air:
      ! On calcule q entre ijb+1 et ije -> on fait pareil pour ratio 
      ! puis on boucle en longitude
      if (nqfils(iq).gt.0) then  
       do ifils=1,nqdesc(iq)
         iq2=iqfils(ifils,iq)  
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)   
         DO l=1,llm
          DO ij=ijb+1,ije
            q(ij,l,iq2)=q(ij,l,iq)*Ratio(ij,l,iq2)            
          enddo
          DO ij=ijb+iip1-1,ije,iip1
             q(ij-iim,l,iq2)=q(ij,l,iq2)
          enddo ! DO ij=ijb+iip1-1,ije,iip1
         enddo !DO l=1,llm
c$OMP END DO NOWAIT
        enddo !do ifils=1,nqdesc(iq)
      endif !if (nqfils(iq).gt.0) then

      !write(*,*) 'vlsplt 399: iq,ijb_x=',iq,ijb_x
c     CALL SCOPY((jjm-1)*llm,q(iip1+iip1,1),iip1,q(iip2,1),iip1)
c     CALL SCOPY((jjm-1)*llm,masse(iip1+iip1,1),iip1,masse(iip2,1),iip1)


      RETURN
      END


      RECURSIVE SUBROUTINE vly_loc(q,pente_max,masse,masse_adv_v,iq)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     q,masse_adv_v,w sont des arguments d'entree  pour le s-pg ....
c     dq 	       sont des arguments de sortie pour le s-pg ....
c
c
c   --------------------------------------------------------------------
      USE parallel_lmdz
      USE infotrac, ONLY : nqtot,nqfils,nqdesc,iqfils ! CRisi
      USE comconst_mod, ONLY: pi
      IMPLICIT NONE
c
      include "dimensions.h"
      include "paramet.h"
      include "comgeom.h"
c
c
c   Arguments:
c   ----------
      REAL masse(ijb_u:ije_u,llm,nqtot),pente_max
      REAL masse_adv_v( ijb_v:ije_v,llm)
      REAL q(ijb_u:ije_u,llm,nqtot), dq( ijb_u:ije_u,llm)
      INTEGER iq ! CRisi
c
c      Local 
c   ---------
c
      INTEGER i,ij,l
c
      REAL airej2,airejjm,airescb(iim),airesch(iim)
      REAL dyq(ijb_u:ije_u,llm),dyqv(ijb_v:ije_v),zdvm(ijb_u:ije_u,llm)
      REAL adyqv(ijb_v:ije_v),dyqmax(ijb_u:ije_u)
      REAL qbyv(ijb_v:ije_v,llm)

      REAL qpns,qpsn,appn,apps,dyn1,dys1,dyn2,dys2,newmasse,fn,fs
c     REAL newq,oldmasse
      Logical extremum,first,testcpu
      REAL temps0,temps1,temps2,temps3,temps4,temps5,second
      SAVE temps0,temps1,temps2,temps3,temps4,temps5
c$OMP THREADPRIVATE(temps0,temps1,temps2,temps3,temps4,temps5)
      SAVE first,testcpu
c$OMP THREADPRIVATE(first,testcpu)

      REAL convpn,convps,convmpn,convmps
      real massepn,masseps,qpn,qps
      REAL sinlon(iip1),sinlondlon(iip1)
      REAL coslon(iip1),coslondlon(iip1)
      SAVE sinlon,coslon,sinlondlon,coslondlon
c$OMP THREADPRIVATE(sinlon,coslon,sinlondlon,coslondlon)
      SAVE airej2,airejjm
c$OMP THREADPRIVATE(airej2,airejjm)

      REAL Ratio(ijb_u:ije_u,llm,nqtot) ! CRisi
      INTEGER ifils,iq2 ! CRisi
c
c
      REAL      SSUM
      EXTERNAL  SSUM

      DATA first,testcpu/.true.,.false./
      DATA temps0,temps1,temps2,temps3,temps4,temps5/0.,0.,0.,0.,0.,0./
      INTEGER ijb,ije

      ijb=ij_begin-2*iip1
      ije=ij_end+2*iip1  
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end

      IF(first) THEN
         PRINT*,'Shema  Amont nouveau  appele dans  Vanleer   '
         first=.false.
         do i=2,iip1
            coslon(i)=cos(rlonv(i))
            sinlon(i)=sin(rlonv(i))
            coslondlon(i)=coslon(i)*(rlonu(i)-rlonu(i-1))/pi
            sinlondlon(i)=sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
         ENDDO
         coslon(1)=coslon(iip1)
         coslondlon(1)=coslondlon(iip1)
         sinlon(1)=sinlon(iip1)
         sinlondlon(1)=sinlondlon(iip1)
         airej2 = SSUM( iim, aire(iip2), 1 )
         airejjm= SSUM( iim, aire(ip1jm -iim), 1 ) 
      ENDIF

c
c	PRINT*,'CALCUL EN LATITUDE'

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1, llm
c
c   --------------------------------
c      CALCUL EN LATITUDE
c   --------------------------------

c   On commence par calculer la valeur du traceur moyenne sur le premier cercle
c   de latitude autour du pole (qpns pour le pole nord et qpsn pour
c    le pole nord) qui sera utilisee pour evaluer les pentes au pole.
      
      if (pole_nord) then
        DO i = 1, iim
          airescb(i) = aire(i+ iip1) * q(i+ iip1,l,iq)
        ENDDO
        qpns   = SSUM( iim,  airescb ,1 ) / airej2
      endif
      
      if (pole_sud) then
        DO i = 1, iim
          airesch(i) = aire(i+ ip1jm- iip1) * q(i+ ip1jm- iip1,l,iq)
        ENDDO
        qpsn   = SSUM( iim,  airesch ,1 ) / airejjm
      endif
      
c   calcul des pentes aux points v

      ijb=ij_begin-2*iip1
      ije=ij_end+iip1
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end-iip1
      
      ! on a besoin de q entre ij_begin-2*iip1 et ij_end+2*iip1
      ! Si pole sud, entre ij_begin-2*iip1 et ij_end
      ! Si pole Nord, entre ij_begin et ij_end+2*iip1
      DO ij=ijb,ije
         dyqv(ij)=q(ij,l,iq)-q(ij+iip1,l,iq)
         adyqv(ij)=abs(dyqv(ij))
      ENDDO
 

c   calcul des pentes aux points scalaires
      ijb=ij_begin-iip1
      ije=ij_end+iip1
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
      DO ij=ijb,ije
         dyq(ij,l)=.5*(dyqv(ij-iip1)+dyqv(ij))
         dyqmax(ij)=min(adyqv(ij-iip1),adyqv(ij))
         dyqmax(ij)=pente_max*dyqmax(ij)
      ENDDO

c   calcul des pentes aux poles
      IF (pole_nord) THEN
        DO ij=1,iip1
           dyq(ij,l)=qpns-q(ij+iip1,l,iq)
        ENDDO
        
        dyn1=0.
        dyn2=0.
        DO ij=1,iim
          dyn1=dyn1+sinlondlon(ij)*dyq(ij,l)
          dyn2=dyn2+coslondlon(ij)*dyq(ij,l)
        ENDDO
        DO ij=1,iip1
          dyq(ij,l)=dyn1*sinlon(ij)+dyn2*coslon(ij)
        ENDDO
        
        DO ij=1,iip1
         dyq(ij,l)=0.
        ENDDO
c ym tout cela ne sert pas a grand chose
      ENDIF
      
      IF (pole_sud) THEN

        DO ij=1,iip1
           dyq(ip1jm+ij,l)=q(ip1jm+ij-iip1,l,iq)-qpsn
        ENDDO

        dys1=0.
        dys2=0.

        DO ij=1,iim
          dys1=dys1+sinlondlon(ij)*dyq(ip1jm+ij,l)
          dys2=dys2+coslondlon(ij)*dyq(ip1jm+ij,l)
        ENDDO

        DO ij=1,iip1
          dyq(ip1jm+ij,l)=dys1*sinlon(ij)+dys2*coslon(ij)
        ENDDO
        
        DO ij=1,iip1
         dyq(ip1jm+ij,l)=0.
        ENDDO
c ym tout cela ne sert pas a grand chose
      ENDIF

c   filtrage de la derivee

c   calcul des pentes limites aux poles
c ym partie inutile
c      goto 8888
c      fn=1.
c      fs=1.
c      DO ij=1,iim
c         IF(pente_max*adyqv(ij).lt.abs(dyq(ij,l))) THEN
c            fn=min(pente_max*adyqv(ij)/abs(dyq(ij,l)),fn)
c         ENDIF
c      IF(pente_max*adyqv(ij+ip1jm-iip1).lt.abs(dyq(ij+ip1jm,l))) THEN
c         fs=min(pente_max*adyqv(ij+ip1jm-iip1)/abs(dyq(ij+ip1jm,l)),fs)
c         ENDIF
c      ENDDO
c      DO ij=1,iip1
c         dyq(ij,l)=fn*dyq(ij,l)
c         dyq(ip1jm+ij,l)=fs*dyq(ip1jm+ij,l)
c      ENDDO
c 8888    continue


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  En memoire de dIFferents tests sur la 
C  limitation des pentes aux poles.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PRINT*,dyq(1)
C     PRINT*,dyqv(iip1+1)
C     appn=abs(dyq(1)/dyqv(iip1+1))
C     PRINT*,dyq(ip1jm+1)
C     PRINT*,dyqv(ip1jm-iip1+1)
C     apps=abs(dyq(ip1jm+1)/dyqv(ip1jm-iip1+1))
C     DO ij=2,iim
C        appn=amax1(abs(dyq(ij)/dyqv(ij)),appn)
C        apps=amax1(abs(dyq(ip1jm+ij)/dyqv(ip1jm-iip1+ij)),apps)
C     ENDDO
C     appn=min(pente_max/appn,1.)
C     apps=min(pente_max/apps,1.)
C
C
C   cas ou on a un extremum au pole
C
C     IF(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
C    &   appn=0.
C     IF(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
C    &   dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
C    &   apps=0.
C
C   limitation des pentes aux poles
C     DO ij=1,iip1
C        dyq(ij)=appn*dyq(ij)
C        dyq(ip1jm+ij)=apps*dyq(ip1jm+ij)
C     ENDDO
C
C   test
C      DO ij=1,iip1
C         dyq(iip1+ij)=0.
C         dyq(ip1jm+ij-iip1)=0.
C      ENDDO
C      DO ij=1,ip1jmp1
C         dyq(ij)=dyq(ij)*cos(rlatu((ij-1)/iip1+1))
C      ENDDO
C
C changement 10 07 96
C     IF(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
C    &   THEN
C        DO ij=1,iip1
C           dyqmax(ij)=0.
C        ENDDO
C     ELSE
C        DO ij=1,iip1
C           dyqmax(ij)=pente_max*abs(dyqv(ij))
C        ENDDO
C     ENDIF
C
C     IF(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
C    & dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
C    &THEN
C        DO ij=ip1jm+1,ip1jmp1
C           dyqmax(ij)=0.
C        ENDDO
C     ELSE
C        DO ij=ip1jm+1,ip1jmp1
C           dyqmax(ij)=pente_max*abs(dyqv(ij-iip1))
C        ENDDO
C     ENDIF
C   fin changement 10 07 96
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c   calcul des pentes limitees
      ijb=ij_begin-iip1
      ije=ij_end+iip1
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1

      DO ij=ijb,ije
         IF(dyqv(ij)*dyqv(ij-iip1).gt.0.) THEN
            dyq(ij,l)=sign(min(abs(dyq(ij,l)),dyqmax(ij)),dyq(ij,l))
         ELSE
            dyq(ij,l)=0.
         ENDIF
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
          IF(masse_adv_v(ij,l).gt.0) THEN
              qbyv(ij,l)=q(ij+iip1,l,iq)+dyq(ij+iip1,l)*
     ,                   0.5*(1.-masse_adv_v(ij,l)
     ,                   /masse(ij+iip1,l,iq))
          ELSE
              qbyv(ij,l)=q(ij,l,iq)-dyq(ij,l)*
     ,                   0.5*(1.+masse_adv_v(ij,l)/masse(ij,l,iq))
          ENDIF
          qbyv(ij,l)=masse_adv_v(ij,l)*qbyv(ij,l)
       ENDDO
      ENDDO
c$OMP END DO NOWAIT

! CRisi: appel récursif de l'advection sur les fils.
! Il faut faire ça avant d'avoir mis à jour q et masse
      !write(*,*) 'vly 689: iq,nqfils(iq)=',iq,nqfils(iq) 

      ijb=ij_begin-2*iip1
      ije=ij_end+2*iip1
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end
   
      if (nqfils(iq).gt.0) then  
       do ifils=1,nqdesc(iq)
         iq2=iqfils(ifils,iq)
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
         DO l=1,llm
         DO ij=ijb,ije
           masse(ij,l,iq2)=masse(ij,l,iq)*q(ij,l,iq)
           Ratio(ij,l,iq2)=q(ij,l,iq2)/q(ij,l,iq)     
          enddo   
         enddo
c$OMP END DO NOWAIT
        enddo !do ifils=1,nqdesc(iq)

        do ifils=1,nqfils(iq)
         iq2=iqfils(ifils,iq)
         call vly_loc(Ratio,pente_max,masse,qbyv,iq2)
        enddo !do ifils=1,nqfils(iq)
      endif !if (nqfils(iq).gt.0) then
! end CRisi
      
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l=1,llm
         DO ij=ijb,ije
            newmasse=masse(ij,l,iq)
     &         +masse_adv_v(ij,l)-masse_adv_v(ij-iip1,l)

            q(ij,l,iq)=(q(ij,l,iq)*masse(ij,l,iq)+qbyv(ij,l)
     &         -qbyv(ij-iip1,l))/newmasse

            masse(ij,l,iq)=newmasse

         ENDDO


c.-. ancienne version
c        convpn=SSUM(iim,qbyv(1,l),1)/apoln
c        convmpn=ssum(iim,masse_adv_v(1,l),1)/apoln
         if (pole_nord) then
           convpn=SSUM(iim,qbyv(1,l),1)
           convmpn=ssum(iim,masse_adv_v(1,l),1)
           massepn=ssum(iim,masse(1,l,iq),1)
           qpn=0.
           do ij=1,iim
              qpn=qpn+masse(ij,l,iq)*q(ij,l,iq)
           enddo
           qpn=(qpn+convpn)/(massepn+convmpn)
           do ij=1,iip1
              q(ij,l,iq)=qpn
           enddo
         endif
         
c        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)/apols
c        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)/apols
         
         if (pole_sud) then
         
           convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
           convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
           masseps=ssum(iim, masse(ip1jm+1,l,iq),1)
           qps=0.
           do ij = ip1jm+1,ip1jmp1-1
              qps=qps+masse(ij,l,iq)*q(ij,l,iq)
           enddo
           qps=(qps+convps)/(masseps+convmps)
           do ij=ip1jm+1,ip1jmp1
              q(ij,l,iq)=qps
           enddo
         endif
c.-. fin ancienne version

c._. nouvelle version
c        convpn=SSUM(iim,qbyv(1,l),1)
c        convmpn=ssum(iim,masse_adv_v(1,l),1)
c        oldmasse=ssum(iim,masse(1,l),1)
c        newmasse=oldmasse+convmpn
c        newq=(q(1,l)*oldmasse+convpn)/newmasse
c        newmasse=newmasse/apoln
c        DO ij = 1,iip1
c           q(ij,l)=newq
c           masse(ij,l,iq)=newmasse*aire(ij)
c        ENDDO
c        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
c        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
c        oldmasse=ssum(iim,masse(ip1jm-iim,l),1)
c        newmasse=oldmasse+convmps
c        newq=(q(ip1jmp1,l)*oldmasse+convps)/newmasse
c        newmasse=newmasse/apols
c        DO ij = ip1jm+1,ip1jmp1
c           q(ij,l)=newq
c           masse(ij,l,iq)=newmasse*aire(ij)
c        ENDDO
c._. fin nouvelle version
      ENDDO
c$OMP END DO NOWAIT

! retablir les fils en rapport de melange par rapport a l'air:
      ijb=ij_begin
      ije=ij_end
!      if (pole_nord) ijb=ij_begin
!      if (pole_sud)  ije=ij_end

      if (nqfils(iq).gt.0) then  
       do ifils=1,nqdesc(iq)
         iq2=iqfils(ifils,iq)  
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)   
         DO l=1,llm
          DO ij=ijb,ije
            q(ij,l,iq2)=q(ij,l,iq)*Ratio(ij,l,iq2)            
          enddo
         enddo
c$OMP END DO NOWAIT
        enddo !do ifils=1,nqdesc(iq)
      endif !if (nqfils(iq).gt.0) then


      RETURN
      END
      
      
      
      RECURSIVE SUBROUTINE vlz_loc(q,pente_max,masse,w,ijb_x,ije_x,iq)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c    q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c     dq 	       sont des arguments de sortie pour le s-pg ....
c
c
c   --------------------------------------------------------------------
      USE parallel_lmdz
      USE vlz_mod
      USE infotrac, ONLY : nqtot,nqfils,nqdesc,iqfils ! CRisi  
      IMPLICIT NONE
c
      include "dimensions.h"
      include "paramet.h"
      include "iniprint.h"
c
c
c   Arguments:
c   ----------
      REAL masse(ijb_u:ije_u,llm,nqtot),pente_max
      REAL q(ijb_u:ije_u,llm,nqtot)
      REAL w(ijb_u:ije_u,llm+1,nqtot)
      INTEGER iq
c
c      Local 
c   ---------
c
      INTEGER i,ij,l,j,ii

      REAL,DIMENSION(ijb_u:ije_u,llm+1) :: wresi,morig,qorig,dzqorig
      INTEGER,DIMENSION(ijb_u:ije_u,llm+1) :: lorig
      INTEGER,SAVE :: countcfl
!$OMP THREADPRIVATE(countcfl)
c
      REAL newmasse

      REAL dzqmax
      REAL sigw

      LOGICAL testcpu
      SAVE testcpu
c$OMP THREADPRIVATE(testcpu)
      REAL temps0,temps1,temps2,temps3,temps4,temps5,second
      SAVE temps0,temps1,temps2,temps3,temps4,temps5
c$OMP THREADPRIVATE(temps0,temps1,temps2,temps3,temps4,temps5)

      REAL      SSUM
      EXTERNAL  SSUM

      DATA testcpu/.false./
      DATA temps0,temps1,temps2,temps3,temps4,temps5/0.,0.,0.,0.,0.,0./
      INTEGER ijb,ije,ijb_x,ije_x
      LOGICAL,SAVE :: first=.TRUE.
!$OMP THREADPRIVATE(first)

      !REAL masseq(ijb_u:ije_u,llm,nqtot),Ratio(ijb_u:ije_u,llm,nqtot) ! CRisi
      ! Ces varibles doivent être déclarées en pointer et en save dans
      ! vlz_loc si on veut qu'elles soient vues par tous les threads.  
      INTEGER ifils,iq2 ! CRisi


      IF (first) THEN
       first=.FALSE.
      ENDIF              
c    On oriente tout dans le sens de la pression c'est a dire dans le
c    sens de W

      !write(*,*) 'vlsplt 926: entree dans vlz_loc, iq=',iq
#ifdef BIDON
      IF(testcpu) THEN
         temps0=second(0.)
      ENDIF
#endif

      ijb=ijb_x
      ije=ije_x

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l=2,llm
         DO ij=ijb,ije
            dzqw(ij,l)=q(ij,l-1,iq)-q(ij,l,iq)
            adzqw(ij,l)=abs(dzqw(ij,l))
         ENDDO
      ENDDO
c$OMP END DO

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=2,llm-1
         DO ij=ijb,ije
#ifdef CRAY
            dzq(ij,l)=0.5*
     ,      cvmgp(dzqw(ij,l)+dzqw(ij,l+1),0.,dzqw(ij,l)*dzqw(ij,l+1))
#else
            IF(dzqw(ij,l)*dzqw(ij,l+1).gt.0.) THEN
                dzq(ij,l)=0.5*(dzqw(ij,l)+dzqw(ij,l+1))
            ELSE
                dzq(ij,l)=0.
            ENDIF
#endif
            dzqmax=pente_max*min(adzqw(ij,l),adzqw(ij,l+1))
            dzq(ij,l)=sign(min(abs(dzq(ij,l)),dzqmax),dzq(ij,l))
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

c$OMP MASTER
      DO ij=ijb,ije
         dzq(ij,1)=0.
         dzq(ij,llm)=0.
      ENDDO
c$OMP END MASTER
c$OMP BARRIER
#ifdef BIDON
      IF(testcpu) THEN
         temps1=temps1+second(0.)-temps0
      ENDIF
#endif

!--------------------------------------------------------
! On repere les points qui violent le CFL (|w| > masse)
!--------------------------------------------------------

      countcfl=0
!     print*,'vlz nouveau'
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 2,llm
         DO ij = ijb,ije
          IF(  (w(ij,l,iq)>0.AND.w(ij,l,iq)>masse(ij,l,iq))
     s    .OR. (w(ij,l,iq)<=0.AND.ABS(w(ij,l,iq))>masse(ij,l-1,iq)) )
     s    countcfl=countcfl+1
         ENDDO
      ENDDO
c$OMP END DO NOWAIT    

c ---------------------------------------------------------------
c  Identification des mailles ou on viole le CFL : w > masse
c ---------------------------------------------------------------

      IF (countcfl==0) THEN

c ---------------------------------------------------------------
c   .... calcul des termes d'advection verticale  .......
c     Dans le cas où le  |w| < masse partout.
c     Version d'origine
c     Pourrait etre enleve si on voit que le code plus general
c     est aussi rapide
c ---------------------------------------------------------------

c calcul de  - d( q   * w )/ d(sigma)    qu'on ajoute a  dq pour calculer dq

       !write(*,*) 'vlz 982,ijb,ije=',ijb,ije
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
       DO l = 1,llm-1
         do  ij = ijb,ije
          IF(w(ij,l+1,iq).gt.0.) THEN
             sigw=w(ij,l+1,iq)/masse(ij,l+1,iq)
             wq(ij,l+1,iq)=w(ij,l+1,iq)*(q(ij,l+1,iq)
     :           +0.5*(1.-sigw)*dzq(ij,l+1))
          ELSE
             sigw=w(ij,l+1,iq)/masse(ij,l,iq)
             wq(ij,l+1,iq)=w(ij,l+1,iq)*(q(ij,l,iq)
     :           -0.5*(1.+sigw)*dzq(ij,l))
          ENDIF
         ENDDO
       ENDDO
c$OMP END DO NOWAIT    
       !write(*,*) 'vlz 1001'    

      ELSE ! countcfl>=1

      IF (prt_level>9) THEN
        WRITE(lunout,*)'vlz passage dans le non local'
      ENDIF
c ---------------------------------------------------------------
c  Debut du traitement du cas ou on viole le CFL : w > masse
c ---------------------------------------------------------------

c Initialisation

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
       DO l = 2,llm
         DO ij = ijb,ije
            wresi(ij,l)=w(ij,l,iq)
            wq(ij,l,iq)=0.
            IF(w(ij,l,iq).gt.0.) THEN
               lorig(ij,l)=l
               morig(ij,l)=masse(ij,l,iq)
               qorig(ij,l)=q(ij,l,iq)
               dzqorig(ij,l)=dzq(ij,l)
            ELSE
               lorig(ij,l)=l-1
               morig(ij,l)=masse(ij,l-1,iq)
               qorig(ij,l)=q(ij,l-1,iq)
               dzqorig(ij,l)=dzq(ij,l-1)
            ENDIF
         ENDDO
       ENDDO
c$OMP END DO NO WAIT

c Reindicage vertical en accumulant les flux sur
c  les mailles qui viollent le CFL
c  on itère jusqu'à ce que tous les poins satisfassent
c  le critère
      DO WHILE (countcfl>=1) 
        IF (prt_level>9) THEN
          WRITE(lunout,*)'On viole le CFL Vertical sur ',countcfl,' pts'
        ENDIF
      countcfl=0

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 2,llm
         DO ij = ijb,ije
          IF (ABS(wresi(ij,l))>morig(ij,l)) THEN
             countcfl=countcfl+1
! rm : les 8 lignes ci dessous pourraient sans doute s'ecrire
! avec la fonction sign
             IF(w(ij,l,iq)>0.) THEN
                wresi(ij,l)=wresi(ij,l)-morig(ij,l)
                wq(ij,l,iq)=wq(ij,l,iq)+morig(ij,l)*qorig(ij,l)
                lorig(ij,l)=lorig(ij,l)+1
             ELSE
                wresi(ij,l)=wresi(ij,l)+morig(ij,l)
                wq(ij,l,iq)=wq(ij,l,iq)-morig(ij,l)*qorig(ij,l)
                lorig(ij,l)=lorig(ij,l)-1
             ENDIF
             morig(ij,l)=masse(ij,lorig(ij,l),iq)
             qorig(ij,l)=q(ij,lorig(ij,l),iq)
             dzqorig(ij,l)=dzq(ij,lorig(ij,l))
          ENDIF
         ENDDO
      ENDDO
c$OMP END DO NO WAIT

      ENDDO ! WHILE (countcfl>=1)

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
       DO l = 2,llm
         do  ij = ijb,ije
          sigw=wresi(ij,l)/morig(ij,l)
          IF(w(ij,l,iq).gt.0.) THEN
             wq(ij,l,iq)=wq(ij,l,iq)+wresi(ij,l)*(qorig(ij,l)
     :           +0.5*(1.-sigw)*dzqorig(ij,l))
          ELSE
             wq(ij,l,iq)=wq(ij,l,iq)+wresi(ij,l)*(qorig(ij,l)
     :           -0.5*(1.+sigw)*dzqorig(ij,l))
          ENDIF
         ENDDO
       ENDDO
c$OMP END DO NOWAIT    


       ENDIF ! councfl=0



c$OMP MASTER
       DO ij=ijb,ije
          wq(ij,llm+1,iq)=0.
          wq(ij,1,iq)=0.
       ENDDO
c$OMP END MASTER
c$OMP BARRIER

! CRisi: appel récursif de l'advection sur les fils.
! Il faut faire ça avant d'avoir mis à jour q et masse
      !write(*,*) 'vlsplt 942: iq,nqfils(iq)=',iq,nqfils(iq)
      if (nqfils(iq).gt.0) then  
       do ifils=1,nqdesc(iq)
         iq2=iqfils(ifils,iq)
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l=1,llm
          DO ij=ijb,ije
           masse(ij,l,iq2)=masse(ij,l,iq)*q(ij,l,iq)
           Ratio(ij,l,iq2)=q(ij,l,iq2)/q(ij,l,iq)
           !wq(ij,l,iq2)=wq(ij,l,iq) ! correction bug le 15mai2015
           w(ij,l,iq2)=wq(ij,l,iq)
          enddo   
         enddo
c$OMP END DO NOWAIT
        enddo !do ifils=1,nqdesc(iq)
c$OMP BARRIER

        do ifils=1,nqfils(iq)
         iq2=iqfils(ifils,iq)
         call vlz_loc(Ratio,pente_max,masse,w,ijb_x,ije_x,iq2)
        enddo !do ifils=1,nqfils(iq)
      endif !if (nqfils(iq).gt.0) then
! end CRisi  

! CRisi: On rajoute ici une barrière car on veut être sur que tous les
! wq soient synchronisés

c$OMP BARRIER
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij=ijb,ije
            newmasse=masse(ij,l,iq)+w(ij,l+1,iq)-w(ij,l,iq)
            q(ij,l,iq)=(q(ij,l,iq)*masse(ij,l,iq)
     &         +wq(ij,l+1,iq)-wq(ij,l,iq))
     &         /newmasse
            masse(ij,l,iq)=newmasse
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

      
! retablir les fils en rapport de melange par rapport a l'air:
      if (nqfils(iq).gt.0) then  
       do ifils=1,nqdesc(iq)
         iq2=iqfils(ifils,iq)  
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)   
         DO l=1,llm
          DO ij=ijb,ije
            q(ij,l,iq2)=q(ij,l,iq)*Ratio(ij,l,iq2)            
          enddo
         enddo
c$OMP END DO NOWAIT
        enddo !do ifils=1,nqdesc(iq)
      endif !if (nqfils(iq).gt.0) then

      RETURN
      END
c      SUBROUTINE minmaxq(zq,qmin,qmax,comment)
c
c#include "dimensions.h"
c#include "paramet.h"

c      CHARACTER*(*) comment
c      real qmin,qmax
c      real zq(ip1jmp1,llm)

c      INTEGER jadrs(ip1jmp1), jbad, k, i


c      DO k = 1, llm
c         jbad = 0
c         DO i = 1, ip1jmp1
c         IF (zq(i,k).GT.qmax .OR. zq(i,k).LT.qmin) THEN
c            jbad = jbad + 1
c            jadrs(jbad) = i
c         ENDIF
c         ENDDO
c         IF (jbad.GT.0) THEN
c         PRINT*, comment
c         DO i = 1, jbad
cc            PRINT*, "i,k,zq=", jadrs(i),k,zq(jadrs(i),k)
c         ENDDO
c         ENDIF
c      ENDDO

c      return
c      end




