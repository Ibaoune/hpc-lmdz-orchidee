!
! $Header$
!
      SUBROUTINE yamada_c(ngrid,timestep,plev,play &
     &   ,pu,pv,pt,d_u,d_v,d_t,cd,q2,km,kn,kq,d_t_diss,ustar &
     &   ,iflag_pbl)
      USE dimphy, ONLY: klon, klev
      USE print_control_mod, ONLY: prt_level
      USE ioipsl_getin_p_mod, ONLY : getin_p

      IMPLICIT NONE
#include "YOMCST.h"
!
! timestep : pas de temps
! g  : g
! zlev : altitude a chaque niveau (interface inferieure de la couche
!        de meme indice)
! zlay : altitude au centre de chaque couche
! u,v : vitesse au centre de chaque couche
!       (en entree : la valeur au debut du pas de temps)
! teta : temperature potentielle au centre de chaque couche
!        (en entree : la valeur au debut du pas de temps)
! cd : cdrag
!      (en entree : la valeur au debut du pas de temps)
! q2 : $q^2$ au bas de chaque couche
!      (en entree : la valeur au debut du pas de temps)
!      (en sortie : la valeur a la fin du pas de temps)
! km : diffusivite turbulente de quantite de mouvement (au bas de chaque
!      couche)
!      (en sortie : la valeur a la fin du pas de temps)
! kn : diffusivite turbulente des scalaires (au bas de chaque couche)
!      (en sortie : la valeur a la fin du pas de temps)
!
!  iflag_pbl doit valoir entre 6 et 9
!      l=6, on prend  systematiquement une longueur d'equilibre
!    iflag_pbl=6 : MY 2.0
!    iflag_pbl=7 : MY 2.0.Fournier
!    iflag_pbl=8/9 : MY 2.5
!       iflag_pbl=8 with special obsolete treatments for convergence
!       with Cmpi5 NPv3.1 simulations
!    iflag_pbl=10/11 :  New scheme M2 and N2 explicit and dissiptation exact
!    iflag_pbl=12 = 11 with vertical diffusion off q2
!
!  2013/04/01 (FH hourdin@lmd.jussieu.fr)
!     Correction for very stable PBLs (iflag_pbl=10 and 11)
!     iflag_pbl=8 converges numerically with NPv3.1
!     iflag_pbl=11 -> the model starts with NP from start files created by ce0l
!                  -> the model can run with longer time-steps.
!.......................................................................

      REAL, DIMENSION(klon,klev) :: d_u,d_v,d_t
      REAL, DIMENSION(klon,klev) :: pu,pv,pt
      REAL, DIMENSION(klon,klev) :: d_t_diss

      REAL timestep
      real plev(klon,klev+1)
      real play(klon,klev)
      real ustar(klon)
      real kmin,qmin,pblhmin(klon),coriol(klon)
      REAL zlev(klon,klev+1)
      REAL zlay(klon,klev)
      REAL zu(klon,klev)
      REAL zv(klon,klev)
      REAL zt(klon,klev)
      REAL teta(klon,klev)
      REAL cd(klon)
      REAL q2(klon,klev+1),qpre
      REAL unsdz(klon,klev)
      REAL unsdzdec(klon,klev+1)

      REAL km(klon,klev)
      REAL kmpre(klon,klev+1),tmp2
      REAL mpre(klon,klev+1)
      REAL kn(klon,klev)
      REAL kq(klon,klev)
      real ff(klon,klev+1),delta(klon,klev+1)
      real aa(klon,klev+1),aa0,aa1
      integer iflag_pbl,ngrid
      integer nlay,nlev

      logical first
      integer ipas
      save first,ipas
!FH/IM     data first,ipas/.true.,0/
      data first,ipas/.false.,0/
!$OMP THREADPRIVATE( first,ipas)
       INTEGER, SAVE :: iflag_tke_diff=0
!$OMP THREADPRIVATE(iflag_tke_diff)


      integer ig,k


      real ri,zrif,zalpha,zsm,zsn
      real rif(klon,klev+1),sm(klon,klev+1),alpha(klon,klev)

      real m2(klon,klev+1),dz(klon,klev+1),zq,n2(klon,klev+1)
      REAL, DIMENSION(klon,klev+1) :: km2,kn2,sqrtq
      real dtetadz(klon,klev+1)
      real m2cstat,mcstat,kmcstat
      real l(klon,klev+1)
      real leff(klon,klev+1)
      real,allocatable,save :: l0(:)
!$OMP THREADPRIVATE(l0)      
      real sq(klon),sqz(klon),zz(klon,klev+1)
      integer iter

      real ric,rifc,b1,kap
      save ric,rifc,b1,kap
      data ric,rifc,b1,kap/0.195,0.191,16.6,0.4/
!$OMP THREADPRIVATE(ric,rifc,b1,kap)
      real frif,falpha,fsm
      real fl,zzz,zl0,zq2,zn2

      real rino(klon,klev+1),smyam(klon,klev),styam(klon,klev)
      real lyam(klon,klev),knyam(klon,klev)
      real w2yam(klon,klev),t2yam(klon,klev)
      logical,save :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)       
      CHARACTER(len=20),PARAMETER :: modname="yamada_c"
REAL, DIMENSION(klon,klev+1) :: fluxu,fluxv,fluxt
REAL, DIMENSION(klon,klev+1) :: dddu,dddv,dddt
REAL, DIMENSION(klon,klev) :: exner,masse
REAL, DIMENSION(klon,klev+1) :: masseb,q2old,q2neg
      LOGICAL okiophys

      frif(ri)=0.6588*(ri+0.1776-sqrt(ri*ri-0.3221*ri+0.03156))
      falpha(ri)=1.318*(0.2231-ri)/(0.2341-ri)
      fsm(ri)=1.96*(0.1912-ri)*(0.2341-ri)/((1.-ri)*(0.2231-ri))
      fl(zzz,zl0,zq2,zn2)= &
     &     max(min(l0(ig)*kap*zlev(ig,k)/(kap*zlev(ig,k)+l0(ig)) &
     &     ,0.5*sqrt(q2(ig,k))/sqrt(max(n2(ig,k),1.e-10))) ,1.)


      okiophys=klon==1
      if (firstcall) then
        CALL getin_p('iflag_tke_diff',iflag_tke_diff)
        allocate(l0(klon))
#define IOPHYS
#ifdef IOPHYS
!        call iophys_ini
#endif
        firstcall=.false.
      endif

   IF (ngrid<=0) RETURN ! Bizarre : on n a pas ce probeleme pour coef_diff_turb

#ifdef IOPHYS
if (okiophys) then
call iophys_ecrit('q2i',klev,'q2 debut my','m2/s2',q2(:,1:klev))
call iophys_ecrit('kmi',klev,'Kz debut my','m/s2',km(:,1:klev))
endif
#endif

      nlay=klev
      nlev=klev+1


!-------------------------------------------------------------------------
! Computation of conservative source terms from the turbulent tendencies
!-------------------------------------------------------------------------


   zalpha=0.5 ! Anciennement 0.5. Essayer de voir pourquoi ? 
   zu(:,:)=pu(:,:)+zalpha*d_u(:,:)
   zv(:,:)=pv(:,:)+zalpha*d_v(:,:)
   zt(:,:)=pt(:,:)+zalpha*d_t(:,:)

   do k=1,klev
      exner(:,k)=(play(:,k)/plev(:,1))**RKAPPA
      masse(:,k)=(plev(:,k)-plev(:,k+1))/RG
      teta(:,k)=zt(:,k)/exner(:,k)
   enddo

! Atmospheric mass at layer interfaces, where the TKE is computed
   masseb(:,:)=0.
   do k=1,klev
      masseb(:,k)=masseb(:,k)+masse(:,k)
      masseb(:,k+1)=masseb(:,k+1)+masse(:,k)
    enddo
    masseb(:,:)=0.5*masseb(:,:)

   zlev(:,1)=0.
   zlay(:,1)=RCPD*teta(:,1)*(1.-exner(:,1))
   do k=1,klev-1
      zlay(:,k+1)=zlay(:,k)+0.5*RCPD*(teta(:,k)+teta(:,k+1))*(exner(:,k)-exner(:,k+1))/RG
      zlev(:,k)=0.5*(zlay(:,k)+zlay(:,k+1)) ! PASBO
   enddo

   fluxu(:,klev+1)=0.
   fluxv(:,klev+1)=0.
   fluxt(:,klev+1)=0.

   do k=klev,1,-1
      fluxu(:,k)=fluxu(:,k+1)+masse(:,k)*d_u(:,k)
      fluxv(:,k)=fluxv(:,k+1)+masse(:,k)*d_v(:,k)
      fluxt(:,k)=fluxt(:,k+1)+masse(:,k)*d_t(:,k)/exner(:,k) ! Flux de theta
   enddo

   dddu(:,1)=2*zu(:,1)*fluxu(:,1)
   dddv(:,1)=2*zv(:,1)*fluxv(:,1)
   dddt(:,1)=(exner(:,1)-1.)*fluxt(:,1)

   do k=2,klev
      dddu(:,k)=(zu(:,k)-zu(:,k-1))*fluxu(:,k)
      dddv(:,k)=(zv(:,k)-zv(:,k-1))*fluxv(:,k)
      dddt(:,k)=(exner(:,k)-exner(:,k-1))*fluxt(:,k)
   enddo
   dddu(:,klev+1)=0.
   dddv(:,klev+1)=0.
   dddt(:,klev+1)=0.

#ifdef IOPHYS
if (okiophys) then
      call iophys_ecrit('zlay',klev,'Geop','m',zlay)
      call iophys_ecrit('teta',klev,'teta','K',teta)
      call iophys_ecrit('temp',klev,'temp','K',zt)
      call iophys_ecrit('pt',klev,'temp','K',pt)
      call iophys_ecrit('pu',klev,'u','m/s',pu)
      call iophys_ecrit('pv',klev,'v','m/s',pv)
      call iophys_ecrit('d_u',klev,'d_u','m/s2',d_u)
      call iophys_ecrit('d_v',klev,'d_v','m/s2',d_v)
      call iophys_ecrit('d_t',klev,'d_t','K/s',d_t)
      call iophys_ecrit('exner',klev,'exner','',exner)
      call iophys_ecrit('masse',klev,'masse','',masse)
      call iophys_ecrit('masseb',klev,'masseb','',masseb)
endif
#endif



      ipas=ipas+1


!.......................................................................
!  les increments verticaux
!.......................................................................
!
!!!!!! allerte !!!!!c
!!!!!! zlev n'est pas declare a nlev !!!!!c
!!!!!! ---->
                                                      DO ig=1,ngrid
            zlev(ig,nlev)=zlay(ig,nlay) &
     &             +( zlay(ig,nlay) - zlev(ig,nlev-1) )
                                                      ENDDO
!!!!!! <----
!!!!!! allerte !!!!!c
!
      DO k=1,nlay
                                                      DO ig=1,ngrid
        unsdz(ig,k)=1.E+0/(zlev(ig,k+1)-zlev(ig,k))
                                                      ENDDO
      ENDDO
                                                      DO ig=1,ngrid
      unsdzdec(ig,1)=1.E+0/(zlay(ig,1)-zlev(ig,1))
                                                      ENDDO
      DO k=2,nlay
                                                      DO ig=1,ngrid
        unsdzdec(ig,k)=1.E+0/(zlay(ig,k)-zlay(ig,k-1))
                                                     ENDDO
      ENDDO
                                                      DO ig=1,ngrid
      unsdzdec(ig,nlay+1)=1.E+0/(zlev(ig,nlay+1)-zlay(ig,nlay))
                                                     ENDDO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing M^2, N^2, Richardson numbers, stability functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do k=2,klev
                                                          do ig=1,ngrid
         dz(ig,k)=zlay(ig,k)-zlay(ig,k-1)
         m2(ig,k)=((zu(ig,k)-zu(ig,k-1))**2+(zv(ig,k)-zv(ig,k-1))**2)/(dz(ig,k)*dz(ig,k))
         dtetadz(ig,k)=(teta(ig,k)-teta(ig,k-1))/dz(ig,k)
         n2(ig,k)=RG*2.*dtetadz(ig,k)/(teta(ig,k-1)+teta(ig,k))
!        n2(ig,k)=0.
         ri=n2(ig,k)/max(m2(ig,k),1.e-10)
         if (ri.lt.ric) then
            rif(ig,k)=frif(ri)
         else
            rif(ig,k)=rifc
         endif
         if(rif(ig,k)<0.16) then
            alpha(ig,k)=falpha(rif(ig,k))
            sm(ig,k)=fsm(rif(ig,k))
         else
            alpha(ig,k)=1.12
            sm(ig,k)=0.085
         endif
         zz(ig,k)=b1*m2(ig,k)*(1.-rif(ig,k))*sm(ig,k)
                                                          enddo
      enddo



!====================================================================
!  Computing the mixing length
!====================================================================

!   Mise a jour de l0
      if (iflag_pbl==8.or.iflag_pbl==10) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Iterative computation of l0
! This version is kept for iflag_pbl only for convergence
! with NPv3.1 Cmip5 simulations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                                          do ig=1,ngrid
      sq(ig)=1.e-10
      sqz(ig)=1.e-10
                                                          enddo
      do k=2,klev-1
                                                          do ig=1,ngrid
        zq=sqrt(q2(ig,k))
        sqz(ig)=sqz(ig)+zq*zlev(ig,k)*(zlay(ig,k)-zlay(ig,k-1))
        sq(ig)=sq(ig)+zq*(zlay(ig,k)-zlay(ig,k-1))
                                                          enddo
      enddo
                                                          do ig=1,ngrid
      l0(ig)=0.2*sqz(ig)/sq(ig)
                                                          enddo
      do k=2,klev
                                                          do ig=1,ngrid
         l(ig,k)=fl(zlev(ig,k),l0(ig),q2(ig,k),n2(ig,k))
                                                          enddo
      enddo
!     print*,'L0 cas 8 ou 10 ',l0

      else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! In all other case, the assymptotic mixing length l0 is imposed (100m)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          l0(:)=150.
          do k=2,klev
                                                          do ig=1,ngrid
             l(ig,k)=fl(zlev(ig,k),l0(ig),q2(ig,k),n2(ig,k))
                                                          enddo
          enddo
!     print*,'L0 cas autres ',l0

      endif


#ifdef IOPHYS
if (okiophys) then
call iophys_ecrit('rif',klev,'Flux Richardson','m',rif(:,1:klev))
call iophys_ecrit('m2',klev,'m2 ','m/s',m2(:,1:klev))
call iophys_ecrit('Km2app',klev,'m2 conserv','m/s',km(:,1:klev)*m2(:,1:klev))
call iophys_ecrit('Km',klev,'Km','m2/s',km(:,1:klev))
endif
#endif


IF (iflag_pbl<20) then
      ! For diagnostics only
      RETURN

ELSE

!  print*,'OK1'

! Evolution of TKE under source terms K M2 and K N2
   leff(:,:)=max(l(:,:),1.)

!##################################################################
!#  IF (iflag_pbl==29) THEN
!#     STOP'Ne pas utiliser iflag_pbl=29'
!#     km2(:,:)=km(:,:)*m2(:,:)
!#     kn2(:,:)=kn2(:,:)*rif(:,:)
!#  ELSEIF (iflag_pbl==25) THEN
! VERSION AVEC LA TKE EN MILIEU DE COUCHE
!#     STOP'Ne pas utiliser iflag_pbl=25'
!#     DO k=1,klev
!#        km2(:,k)=-0.5*(dddu(:,k)+dddv(:,k)+dddu(:,k+1)+dddv(:,k+1)) &
!#        &        /(masse(:,k)*timestep)
!#        kn2(:,k)=rcpd*0.5*(dddt(:,k)+dddt(:,k+1))/(masse(:,k)*timestep)
!#        leff(:,k)=0.5*(leff(:,k)+leff(:,k+1))
!#     ENDDO
!#     km2(:,klev+1)=0. ; kn2(:,klev+1)=0.
!#  ELSE
!#################################################################

      km2(:,:)=-(dddu(:,:)+dddv(:,:))/(masseb(:,:)*timestep)
      kn2(:,:)=rcpd*dddt(:,:)/(masseb(:,:)*timestep)
!   ENDIF
   q2neg(:,:)=q2(:,:)+timestep*(km2(:,:)-kn2(:,:))
   q2(:,:)=min(max(q2neg(:,:),1.e-10),1.e4)

 
#ifdef IOPHYS
if (okiophys) then
      call iophys_ecrit('km2',klev,'m2 conserv','m/s',km2(:,1:klev))
      call iophys_ecrit('kn2',klev,'n2 conserv','m/s',kn2(:,1:klev))
endif
#endif

! Dissipation of TKE
   q2old(:,:)=q2(:,:)
   q2(:,:)=1./(1./sqrt(q2(:,:))+timestep/(2*leff(:,:)*b1))
   q2(:,:)=q2(:,:)*q2(:,:)
!  IF (iflag_pbl<=24) THEN
      DO k=1,klev
         d_t_diss(:,k)=(masseb(:,k)*(q2neg(:,k)-q2(:,k))+masseb(:,k+1)*(q2neg(:,k+1)-q2(:,k+1)))/(2.*rcpd*masse(:,k))
      ENDDO

!###################################################################
!  ELSE IF (iflag_pbl<=27) THEN
!     DO k=1,klev
!        d_t_diss(:,k)=(q2neg(:,k)-q2(:,k))/rcpd
!     ENDDO
!  ENDIF
!  print*,'iflag_pbl ',d_t_diss
!###################################################################


! Compuation of stability functions
!   IF (iflag_pbl/=29) THEN
      DO k=1,klev
      DO ig=1,ngrid
         IF (ABS(km2(ig,k))<=1.e-20) THEN
            rif(ig,k)=0.
         ELSE
            rif(ig,k)=min(kn2(ig,k)/km2(ig,k),rifc)
         ENDIF
         IF (rif(ig,k).lt.0.16) THEN
            alpha(ig,k)=falpha(rif(ig,k))
            sm(ig,k)=fsm(rif(ig,k))
         else
            alpha(ig,k)=1.12
            sm(ig,k)=0.085
         endif
      ENDDO
      ENDDO
!    ENDIF

! Computation of turbulent diffusivities
!  IF (25<=iflag_pbl.and.iflag_pbl<=28) THEN
!    DO k=2,klev
!       sqrtq(:,k)=sqrt(0.5*(q2(:,k)+q2(:,k-1)))
!    ENDDO
!  ELSE
   kq(:,:)=0.
   DO k=1,klev
      ! Coefficient au milieu des couches pour diffuser la TKE
      kq(:,k)=0.5*leff(:,k)*sqrt(q2(:,k))*0.2
   ENDDO

#ifdef IOPHYS
if (okiophys) then
call iophys_ecrit('q2b',klev,'KTE inter','m2/s',q2(:,1:klev))
endif
#endif

  IF (iflag_tke_diff==1) THEN
    CALL vdif_q2(timestep, RG, RD, ngrid, plev, pt, kq, q2)
  ENDIF

   km(:,:)=0.
   kn(:,:)=0.
   DO k=1,klev
      km(:,k)=leff(:,k)*sqrt(q2(:,k))*sm(:,k)
      kn(:,k)=km(:,k)*alpha(:,k)
   ENDDO


#ifdef IOPHYS
if (okiophys) then
call iophys_ecrit('mixingl',klev,'Mixing length','m',leff(:,1:klev))
call iophys_ecrit('rife',klev,'Flux Richardson','m',rif(:,1:klev))
call iophys_ecrit('q2f',klev,'KTE finale','m2/s',q2(:,1:klev))
call iophys_ecrit('q2neg',klev,'KTE non bornee','m2/s',q2neg(:,1:klev))
call iophys_ecrit('alpha',klev,'alpha','',alpha(:,1:klev))
call iophys_ecrit('sm',klev,'sm','',sm(:,1:klev))
call iophys_ecrit('q2f',klev,'KTE finale','m2/s',q2(:,1:klev))
call iophys_ecrit('kmf',klev,'Kz final','m2/s',km(:,1:klev))
call iophys_ecrit('knf',klev,'Kz final','m2/s',kn(:,1:klev))
call iophys_ecrit('kqf',klev,'Kz final','m2/s',kq(:,1:klev))
endif
#endif


ENDIF


!  print*,'OK2'
      RETURN
      END
