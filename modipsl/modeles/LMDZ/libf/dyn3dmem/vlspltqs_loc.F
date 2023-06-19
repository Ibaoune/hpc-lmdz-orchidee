      SUBROUTINE vlxqs_loc(q,pente_max,masse,u_m,qsat,ijb_x,ije_x,iq)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d''advection " pseudo amont " .
c    ********************************************************************
c
c   --------------------------------------------------------------------
      USE parallel_lmdz
      USE infotrac, ONLY : nqtot,nqfils,nqdesc,iqfils ! CRisi  
      IMPLICIT NONE
c
      include "dimensions.h"
      include "paramet.h"
c
c
c   Arguments:
c   ----------
      REAL masse(ijb_u:ije_u,llm,nqtot),pente_max
      REAL u_m( ijb_u:ije_u,llm )
      REAL q(ijb_u:ije_u,llm,nqtot)
      REAL qsat(ijb_u:ije_u,llm)
      INTEGER iq ! CRisi
c
c      Local 
c   ---------
c
      INTEGER ij,l,j,i,iju,ijq,indu(ijnb_u),niju
      INTEGER n0,iadvplus(ijb_u:ije_u,llm),nl(llm)
c
      REAL new_m,zu_m,zdum(ijb_u:ije_u,llm)
      REAL dxq(ijb_u:ije_u,llm),dxqu(ijb_u:ije_u)
      REAL zz(ijb_u:ije_u)
      REAL adxqu(ijb_u:ije_u),dxqmax(ijb_u:ije_u,llm)
      REAL u_mq(ijb_u:ije_u,llm)
      REAL Ratio(ijb_u:ije_u,llm,nqtot) ! CRisi
      INTEGER ifils,iq2 ! CRisi


      REAL      SSUM


      INTEGER ijb,ije,ijb_x,ije_x
      
      !write(*,*) 'vlspltqs 58: entree vlxqs_loc, iq,ijb_x=',
!     &   iq,ijb_x

c   calcul de la pente a droite et a gauche de la maille

c      ijb=ij_begin
c      ije=ij_end

      ijb=ijb_x
      ije=ije_x
        
      if (pole_nord.and.ijb==1) ijb=ijb+iip1
      if (pole_sud.and.ije==ip1jmp1)  ije=ije-iip1
      
      IF (pente_max.gt.-1.e-5) THEN
c     IF (pente_max.gt.10) THEN

c   calcul des pentes avec limitation, Van Leer scheme I:
c   -----------------------------------------------------

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
      
      if (pole_nord) THEN
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm      
          iadvplus(1:iip1,l)=0
        ENDDO
c$OMP END DO NOWAIT
      endif
      
      if (pole_sud)  THEN
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm 
          iadvplus(ip1jm+1:ip1jmp1,l)=0
        ENDDO
c$OMP END DO NOWAIT
      endif
              
c   calcul des flux a gauche et a droite

#ifdef CRAY
c--pas encore modification sur Qsat
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
       DO ij=ijb,ije-1
          zdum(ij,l)=cvmgp(1.-u_m(ij,l)/masse(ij,l,iq),
     ,                     1.+u_m(ij,l)/masse(ij+1,l,iq),
     ,                     u_m(ij,l))
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
c   le rapport de melange de l''air advecte est min(q_vanleer, Qsat_downwind)
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
       DO ij=ijb,ije-1
          IF (u_m(ij,l).gt.0.) THEN
             zdum(ij,l)=1.-u_m(ij,l)/masse(ij,l,iq)
             u_mq(ij,l)=u_m(ij,l)*
     $         min(q(ij,l,iq)+0.5*zdum(ij,l)*dxq(ij,l),qsat(ij+1,l))
          ELSE
             zdum(ij,l)=1.+u_m(ij,l)/masse(ij+1,l,iq)
             u_mq(ij,l)=u_m(ij,l)*
     $         min(q(ij+1,l,iq)-0.5*zdum(ij,l)*dxq(ij+1,l),qsat(ij,l))
          ENDIF
       ENDDO
      ENDDO
c$OMP END DO NOWAIT
#endif


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

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
       DO ij=ijb+iip1-1,ije,iip1
          iadvplus(ij,l)=iadvplus(ij-iim,l)
       ENDDO
      ENDDO
c$OMP END DO NOWAIT



c   traitement special pour le cas ou on advecte en longitude plus que le
c   contenu de la maille.
c   cette partie est mal vectorisee.

c   pas d'influence de la pression saturante (pour l'instant)

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

cym ATTENTION ICI en OpenMP reduction pas forcement nécessaire
cym      IF(n0.gt.1) THEN
cym        IF(n0.gt.0) THEN
ccc      PRINT*,'Nombre de points pour lesquels on advect plus que le'
ccc     &       ,'contenu de la maille : ',n0
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
               !PRINT*,'vlxqs 280: niju,nl',niju,nl(l)

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
                        u_mq(ij,l)=u_mq(ij,l)+q(ijq,l,iq)
     &                     *masse(ijq,l,iq) 
                        zu_m=zu_m-masse(ijq,l,iq)
                        i=mod(i-2+iim,iim)+1
                        ijq=(j-1)*iip1+i
                     ENDDO
c   ajout de la maille non completement advectee
                     u_mq(ij,l)=u_mq(ij,l)+zu_m*(q(ijq,l,iq)
     &                 +0.5*(1.-zu_m/masse(ijq,l,iq))*dxq(ijq,l))
                  ELSE
                     ijq=ij+1
                     i=ijq-(j-1)*iip1
c   accumulation pour les mailles completements advectees
                     do while(-zu_m.gt.masse(ijq,l,iq))
                        u_mq(ij,l)=u_mq(ij,l)-q(ijq,l,iq)
     &                   *masse(ijq,l,iq)
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



c   bouclage en latitude
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
        DO ij=ijb+iip1-1,ije,iip1
           u_mq(ij,l)=u_mq(ij-iim,l)
        ENDDO
      ENDDO
c$OMP END DO NOWAIT

! CRisi: appel récursif de l'advection sur les fils.
! Il faut faire ça avant d'avoir mis à jour q et masse
      !write(*,*) 'vlspltqs 336: iq,ijb_x,nqfils(iq)=',
!     &     iq,ijb_x,nqfils(iq)  

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
        enddo !do ifils=1,nqfils(iq)
        do ifils=1,nqfils(iq)
         iq2=iqfils(ifils,iq)
         !write(*,*) 'vlxqs 349: on appelle vlx pour iq2=',iq2
         call vlx_loc(Ratio,pente_max,masse,u_mq,ijb_x,ije_x,iq2)
        enddo !do ifils=1,nqfils(iq)
      endif !if (nqfils(iq).gt.0) then
! end CRisi

      !write(*,*) 'vlspltqs 360: iq,ijb_x=',iq,ijb_x   

c   calcul des tendances
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij=ijb+1,ije
            new_m=masse(ij,l,iq)+u_m(ij-1,l)-u_m(ij,l)
            q(ij,l,iq)=(q(ij,l,iq)*masse(ij,l,iq)+
     &      u_mq(ij-1,l)-u_mq(ij,l))
     &      /new_m
            masse(ij,l,iq)=new_m
         ENDDO
c   Modif Fred 22 03 96 correction d''un bug (les scopy ci-dessous)
         DO ij=ijb+iip1-1,ije,iip1
            q(ij-iim,l,iq)=q(ij,l,iq)
            masse(ij-iim,l,iq)=masse(ij,l,iq)
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

      !write(*,*) 'vlspltqs 380: iq,ijb_x=',iq,ijb_x

! retablir les fils en rapport de melange par rapport a l'air:
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
         enddo
c$OMP END DO NOWAIT
        enddo !do ifils=1,nqdesc(iq)
      endif !if (nqfils(iq).gt.0) then

      !write(*,*) 'vlspltqs 399: iq,ijb_x=',iq,ijb_x

c     CALL SCOPY((jjm-1)*llm,q(iip1+iip1,1),iip1,q(iip2,1),iip1)
c     CALL SCOPY((jjm-1)*llm,masse(iip1+iip1,1,iq),iip1,masse(iip2,1,iq),iip1)


      RETURN
      END
      SUBROUTINE vlyqs_loc(q,pente_max,masse,masse_adv_v,qsat,iq)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     q,masse_adv_v,w sont des arguments d'entree  pour le s-pg ....
c     qsat                est   un argument de sortie pour le s-pg ....
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
      REAL q(ijb_u:ije_u,llm,nqtot)
      REAL qsat(ijb_u:ije_u,llm)
      INTEGER iq ! CRisi
c
c      Local 
c   ---------
c
      INTEGER i,ij,l
c
      REAL airej2,airejjm,airescb(iim),airesch(iim)
      REAL dyq(ijb_u:ije_u,llm),dyqv(ijb_v:ije_v)
      REAL adyqv(ijb_v:ije_v),dyqmax(ijb_u:ije_u)
      REAL qbyv(ijb_v:ije_v,llm,nqtot)

      REAL qpns,qpsn,dyn1,dys1,dyn2,dys2,newmasse,fn,fs
c     REAL newq,oldmasse
      Logical first
      SAVE first
c$OMP THREADPRIVATE(first)
      REAL convpn,convps,convmpn,convmps
      REAL sinlon(iip1),sinlondlon(iip1)
      REAL coslon(iip1),coslondlon(iip1)
      SAVE sinlon,coslon,sinlondlon,coslondlon
      SAVE airej2,airejjm
c$OMP THREADPRIVATE(sinlon,coslon,sinlondlon,coslondlon)
c$OMP THREADPRIVATE(airej2,airejjm)
c
c
      REAL Ratio(ijb_u:ije_u,llm,nqtot) ! CRisi
      INTEGER ifils,iq2 ! CRisi

      REAL      SSUM

      DATA first/.true./
      INTEGER ijb,ije

      ijb=ij_begin-2*iip1
      ije=ij_end+2*iip1  
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end
      ij=3525
      l=3
      if ((ij.ge.ijb).and.(ij.le.ije)) then
        !write(*,*) 'vlyqs 480: ij,l,iq,ijb,q(ij,l,:)=',
!     &             ij,l,iq,ijb,q(ij,l,:)
      endif  

      IF(first) THEN
         PRINT*,'Shema  Amont nouveau  appele dans  Vanleer   '
         PRINT*,'vlyqs_loc, iq=',iq
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
      
      IF (pole_nord) THEN

c   calcul des pentes aux poles
        DO ij=1,iip1
           dyq(ij,l)=qpns-q(ij+iip1,l,iq)
        ENDDO

c   filtrage de la derivee        
        dyn1=0.
        dyn2=0.
        DO ij=1,iim
          dyn1=dyn1+sinlondlon(ij)*dyq(ij,l)
          dyn2=dyn2+coslondlon(ij)*dyq(ij,l)
        ENDDO
        DO ij=1,iip1
          dyq(ij,l)=dyn1*sinlon(ij)+dyn2*coslon(ij)
        ENDDO

c   calcul des pentes limites aux poles
        fn=1.
        DO ij=1,iim
          IF(pente_max*adyqv(ij).lt.abs(dyq(ij,l))) THEN
            fn=min(pente_max*adyqv(ij)/abs(dyq(ij,l)),fn)
          ENDIF
        ENDDO
      
        DO ij=1,iip1
         dyq(ij,l)=fn*dyq(ij,l)
        ENDDO
          
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
        
c   calcul des pentes limites aux poles        
        fs=1.
        DO ij=1,iim
        IF(pente_max*adyqv(ij+ip1jm-iip1).lt.abs(dyq(ij+ip1jm,l))) THEN
         fs=min(pente_max*adyqv(ij+ip1jm-iip1)/abs(dyq(ij+ip1jm,l)),fs)
        ENDIF
        ENDDO
    
        DO ij=1,iip1
         dyq(ip1jm+ij,l)=fs*dyq(ip1jm+ij,l)
        ENDDO
        
      ENDIF


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
         IF( masse_adv_v(ij,l).GT.0. ) THEN
           qbyv(ij,l,iq)= MIN( qsat(ij+iip1,l), q(ij+iip1,l,iq )  +
     ,      dyq(ij+iip1,l)*0.5*(1.-masse_adv_v(ij,l)
     ,      /masse(ij+iip1,l,iq)))
         ELSE
              qbyv(ij,l,iq)= MIN( qsat(ij,l), q(ij,l,iq) - dyq(ij,l) *
     ,                   0.5*(1.+masse_adv_v(ij,l)/masse(ij,l,iq)) )
         ENDIF
          qbyv(ij,l,iq) = masse_adv_v(ij,l)*qbyv(ij,l,iq)
       ENDDO
      ENDDO
c$OMP END DO NOWAIT

! CRisi: appel récursif de l'advection sur les fils.
! Il faut faire ça avant d'avoir mis à jour q et masse
      !write(*,*) 'vlyqs 689: iq,nqfils(iq)=',iq,nqfils(iq)
      
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
     &      +masse_adv_v(ij,l)-masse_adv_v(ij-iip1,l)
            q(ij,l,iq)=(q(ij,l,iq)*masse(ij,l,iq)+qbyv(ij,l,iq)
     &         -qbyv(ij-iip1,l,iq))/newmasse
            masse(ij,l,iq)=newmasse
         ENDDO
c.-. ancienne version

         IF (pole_nord) THEN

           convpn=SSUM(iim,qbyv(1,l,iq),1)/apoln
           convmpn=ssum(iim,masse_adv_v(1,l),1)/apoln
           DO ij = 1,iip1
              newmasse=masse(ij,l,iq)+convmpn*aire(ij)
              q(ij,l,iq)=(q(ij,l,iq)*masse(ij,l,iq)+convpn*aire(ij))/
     &                 newmasse
              masse(ij,l,iq)=newmasse
           ENDDO
         
         ENDIF
         
         IF (pole_sud) THEN
         
           convps  = -SSUM(iim,qbyv(ip1jm-iim,l,iq),iq,1)/apols
           convmps = -SSUM(iim,masse_adv_v(ip1jm-iim,l),1)/apols
           DO ij = ip1jm+1,ip1jmp1
              newmasse=masse(ij,l,iq)+convmps*aire(ij)
              q(ij,l,iq)=(q(ij,l,iq)*masse(ij,l,iq)+convps*aire(ij))/
     &                 newmasse
              masse(ij,l,iq)=newmasse
           ENDDO
         
         ENDIF
c.-. fin ancienne version

c._. nouvelle version
c        convpn=SSUM(iim,qbyv(1,l,iq),1)
c        convmpn=ssum(iim,masse_adv_v(1,l),1)
c        oldmasse=ssum(iim,masse(1,l,iq),1)
c        newmasse=oldmasse+convmpn
c        newq=(q(1,l,iq)*oldmasse+convpn)/newmasse
c        newmasse=newmasse/apoln
c        DO ij = 1,iip1
c           q(ij,l,iq)=newq
c           masse(ij,l,iq)=newmasse*aire(ij)
c        ENDDO
c        convps=-SSUM(iim,qbyv(ip1jm-iim,l,iq),1)
c        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
c        oldmasse=ssum(iim,masse(ip1jm-iim,l,iq),1)
c        newmasse=oldmasse+convmps
c        newq=(q(ip1jmp1,l,iq)*oldmasse+convps)/newmasse
c        newmasse=newmasse/apols
c        DO ij = ip1jm+1,ip1jmp1
c           q(ij,l,iq)=newq
c           masse(ij,l,iq)=newmasse*aire(ij)
c        ENDDO
c._. fin nouvelle version
      ENDDO
c$OMP END DO NOWAIT

! retablir les fils en rapport de melange par rapport a l'air:
      ijb=ij_begin
      ije=ij_end
!      if (pole_nord) ijb=ij_begin+iip1
!      if (pole_sud)  ije=ij_end-iip1
  
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
