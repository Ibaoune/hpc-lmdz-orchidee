      subroutine thermcell_dq(ngrid,nlay,impl,ptimestep,fm,entr,  &
     &           masse,q,dq,qa,lev_out)
      USE print_control_mod, ONLY: prt_level
      implicit none

!=======================================================================
!
!   Calcul du transport verticale dans la couche limite en presence
!   de "thermiques" explicitement representes
!   calcul du dq/dt une fois qu'on connait les ascendances
!
! Modif 2013/01/04 (FH hourdin@lmd.jussieu.fr)
!  Introduction of an implicit computation of vertical advection in
!  the environment of thermal plumes in thermcell_dq
!  impl =     0 : explicit, 1 : implicit, -1 : old version
!
!=======================================================================

      integer ngrid,nlay,impl

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real entr(ngrid,nlay)
      real q(ngrid,nlay)
      real dq(ngrid,nlay)
      integer lev_out                           ! niveau pour les print

      real qa(ngrid,nlay),detr(ngrid,nlay),wqd(ngrid,nlay+1)

      real zzm

      integer ig,k
      real cfl

      real qold(ngrid,nlay),fqa(ngrid,nlay+1)
      integer niter,iter
      CHARACTER (LEN=20) :: modname='thermcell_dq'
      CHARACTER (LEN=80) :: abort_message


! Old explicite scheme
      if (impl<=-1) then
         call thermcell_dq_o(ngrid,nlay,impl,ptimestep,fm,entr,  &
     &           masse,q,dq,qa,lev_out)
         return
      endif

! Calcul du critere CFL pour l'advection dans la subsidence
      cfl = 0.
      do k=1,nlay
         do ig=1,ngrid
            zzm=masse(ig,k)/ptimestep
            cfl=max(cfl,fm(ig,k)/zzm)
            if (entr(ig,k).gt.zzm) then
               print*,'entr*dt>m,1',k,entr(ig,k)*ptimestep,masse(ig,k)
               abort_message = 'entr dt > m, 1st'
               CALL abort_physic (modname,abort_message,1)
            endif
         enddo
      enddo

      qold=q


      if (prt_level.ge.1) print*,'Q2 THERMCEL_DQ 0'

!   calcul du detrainement
      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
!           print*,'Q2 DQ ',detr(ig,k),fm(ig,k),entr(ig,k)
!test
            if (detr(ig,k).lt.0.) then
               entr(ig,k)=entr(ig,k)-detr(ig,k)
               detr(ig,k)=0.
!               print*,'detr2<0!!!','ig=',ig,'k=',k,'f=',fm(ig,k),
!     s         'f+1=',fm(ig,k+1),'e=',entr(ig,k),'d=',detr(ig,k)
            endif
            if (fm(ig,k+1).lt.0.) then
!               print*,'fm2<0!!!'
            endif
            if (entr(ig,k).lt.0.) then
!               print*,'entr2<0!!!'
            endif
         enddo
      enddo

! Computation of tracer concentrations in the ascending plume
      do ig=1,ngrid
         qa(ig,1)=q(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.  &
     &         1.e-5*masse(ig,k)) then
         qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))  &
     &         /(fm(ig,k+1)+detr(ig,k))
            else
               qa(ig,k)=q(ig,k)
            endif
            if (qa(ig,k).lt.0.) then
!               print*,'qa<0!!!'
            endif
            if (q(ig,k).lt.0.) then
!               print*,'q<0!!!'
            endif
         enddo
      enddo

! Plume vertical flux
      do k=2,nlay-1
         fqa(:,k)=fm(:,k)*qa(:,k-1)
      enddo
      fqa(:,1)=0. ; fqa(:,nlay)=0.


! Trace species evolution
   if (impl==0) then
      do k=1,nlay-1
         q(:,k)=q(:,k)+(fqa(:,k)-fqa(:,k+1)-fm(:,k)*q(:,k)+fm(:,k+1)*q(:,k+1)) &
     &               *ptimestep/masse(:,k)
      enddo
   else
      do k=nlay-1,1,-1
! FH debut de modif : le calcul ci dessous modifiait numériquement
! la concentration quand le flux de masse etait nul car on divisait
! puis multipliait par masse/ptimestep.
!        q(:,k)=(masse(:,k)*q(:,k)/ptimestep+fqa(:,k)-fqa(:,k+1)+fm(:,k+1)*q(:,k+1)) &
!    &               /(fm(:,k)+masse(:,k)/ptimestep)
         q(:,k)=(q(:,k)+ptimestep/masse(:,k)*(fqa(:,k)-fqa(:,k+1)+fm(:,k+1)*q(:,k+1))) &
      &               /(1.+fm(:,k)*ptimestep/masse(:,k))
! FH fin de modif.
      enddo
   endif

! Tendencies
      do k=1,nlay
         do ig=1,ngrid
            dq(ig,k)=(q(ig,k)-qold(ig,k))/ptimestep
            q(ig,k)=qold(ig,k)
         enddo
      enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Obsolete version kept for convergence with Cmip5 NPv3.1 simulations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine thermcell_dq_o(ngrid,nlay,impl,ptimestep,fm,entr,  &
     &           masse,q,dq,qa,lev_out)
      USE print_control_mod, ONLY: prt_level
      implicit none

!=======================================================================
!
!   Calcul du transport verticale dans la couche limite en presence
!   de "thermiques" explicitement representes
!   calcul du dq/dt une fois qu'on connait les ascendances
!
!=======================================================================

      integer ngrid,nlay,impl

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real entr(ngrid,nlay)
      real q(ngrid,nlay)
      real dq(ngrid,nlay)
      integer lev_out                           ! niveau pour les print

      real qa(ngrid,nlay),detr(ngrid,nlay),wqd(ngrid,nlay+1)

      real zzm

      integer ig,k
      real cfl

      real qold(ngrid,nlay)
      real ztimestep
      integer niter,iter
      CHARACTER (LEN=20) :: modname='thermcell_dq'
      CHARACTER (LEN=80) :: abort_message



! Calcul du critere CFL pour l'advection dans la subsidence
      cfl = 0.
      do k=1,nlay
         do ig=1,ngrid
            zzm=masse(ig,k)/ptimestep
            cfl=max(cfl,fm(ig,k)/zzm)
            if (entr(ig,k).gt.zzm) then
               print*,'entr*dt>m,2',k,entr(ig,k)*ptimestep,masse(ig,k)
               abort_message = 'entr dt > m, 2nd'
               CALL abort_physic (modname,abort_message,1)
            endif
         enddo
      enddo

!IM 090508     print*,'CFL CFL CFL CFL ',cfl

#undef CFL
#ifdef CFL
! On subdivise le calcul en niter pas de temps.
      niter=int(cfl)+1
#else
      niter=1
#endif

      ztimestep=ptimestep/niter
      qold=q


do iter=1,niter
      if (prt_level.ge.1) print*,'Q2 THERMCEL_DQ 0'

!   calcul du detrainement
      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
!           print*,'Q2 DQ ',detr(ig,k),fm(ig,k),entr(ig,k)
!test
            if (detr(ig,k).lt.0.) then
               entr(ig,k)=entr(ig,k)-detr(ig,k)
               detr(ig,k)=0.
!               print*,'detr2<0!!!','ig=',ig,'k=',k,'f=',fm(ig,k),
!     s         'f+1=',fm(ig,k+1),'e=',entr(ig,k),'d=',detr(ig,k)
            endif
            if (fm(ig,k+1).lt.0.) then
!               print*,'fm2<0!!!'
            endif
            if (entr(ig,k).lt.0.) then
!               print*,'entr2<0!!!'
            endif
         enddo
      enddo

!   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         qa(ig,1)=q(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ztimestep.gt.  &
     &         1.e-5*masse(ig,k)) then
         qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))  &
     &         /(fm(ig,k+1)+detr(ig,k))
            else
               qa(ig,k)=q(ig,k)
            endif
            if (qa(ig,k).lt.0.) then
!               print*,'qa<0!!!'
            endif
            if (q(ig,k).lt.0.) then
!               print*,'q<0!!!'
            endif
         enddo
      enddo

! Calcul du flux subsident

      do k=2,nlay
         do ig=1,ngrid
#undef centre
#ifdef centre
             wqd(ig,k)=fm(ig,k)*0.5*(q(ig,k-1)+q(ig,k))
#else

#define plusqueun
#ifdef plusqueun
! Schema avec advection sur plus qu'une maille.
            zzm=masse(ig,k)/ztimestep
            if (fm(ig,k)>zzm) then
               wqd(ig,k)=zzm*q(ig,k)+(fm(ig,k)-zzm)*q(ig,k+1)
            else
               wqd(ig,k)=fm(ig,k)*q(ig,k)
            endif
#else
            wqd(ig,k)=fm(ig,k)*q(ig,k)
#endif
#endif

            if (wqd(ig,k).lt.0.) then
!               print*,'wqd<0!!!'
            endif
         enddo
      enddo
      do ig=1,ngrid
         wqd(ig,1)=0.
         wqd(ig,nlay+1)=0.
      enddo
     

! Calcul des tendances
      do k=1,nlay
         do ig=1,ngrid
            q(ig,k)=q(ig,k)+(detr(ig,k)*qa(ig,k)-entr(ig,k)*q(ig,k)  &
     &               -wqd(ig,k)+wqd(ig,k+1))  &
     &               *ztimestep/masse(ig,k)
!            if (dq(ig,k).lt.0.) then
!               print*,'dq<0!!!'
!            endif
         enddo
      enddo


enddo


! Calcul des tendances
      do k=1,nlay
         do ig=1,ngrid
            dq(ig,k)=(q(ig,k)-qold(ig,k))/ptimestep
            q(ig,k)=qold(ig,k)
         enddo
      enddo

      return
      end
