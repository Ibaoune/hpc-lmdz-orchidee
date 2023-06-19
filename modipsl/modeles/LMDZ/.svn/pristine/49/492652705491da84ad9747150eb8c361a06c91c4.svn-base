SUBROUTINE calcratqs(klon,klev,prt_level,lunout,       &
           iflag_ratqs,iflag_con,iflag_cld_th,pdtphys, &
           ratqsbas,ratqshaut,ratqsp0,ratqsdp, &
           tau_ratqs,fact_cldcon,   &
           ptconv,ptconvth,clwcon0th, rnebcon0th,      &
           paprs,pplay,q_seri,zqsat,fm_therm,          &
           ratqs,ratqsc)

implicit none

!========================================================================
! Computation of ratqs, the width of the subrid scale water distribution
! (normalized by the mean value)
! Various options controled by flags iflag_con and iflag_ratqs
! F Hourdin 2012/12/06
!========================================================================

! Declarations

! Input
integer,intent(in) :: klon,klev,prt_level,lunout
integer,intent(in) :: iflag_con,iflag_cld_th,iflag_ratqs
real,intent(in) :: pdtphys,ratqsbas,ratqshaut,fact_cldcon,tau_ratqs
real,intent(in) :: ratqsp0, ratqsdp
real, dimension(klon,klev+1),intent(in) :: paprs
real, dimension(klon,klev),intent(in) :: pplay,q_seri,zqsat,fm_therm
logical, dimension(klon,klev),intent(in) :: ptconv
real, dimension(klon,klev),intent(in) :: rnebcon0th,clwcon0th

! Output
real, dimension(klon,klev),intent(inout) :: ratqs,ratqsc
logical, dimension(klon,klev),intent(inout) :: ptconvth

! local
integer i,k
real, dimension(klon,klev) :: ratqss
real facteur,zfratqs1,zfratqs2

!-------------------------------------------------------------------------
!  Caclul des ratqs
!-------------------------------------------------------------------------

!      print*,'calcul des ratqs'
!   ratqs convectifs a l'ancienne en fonction de q(z=0)-q / q
!   ----------------
!   on ecrase le tableau ratqsc calcule par clouds_gno
      if (iflag_cld_th.eq.1) then
         do k=1,klev
         do i=1,klon
            if(ptconv(i,k)) then
              ratqsc(i,k)=ratqsbas &
              +fact_cldcon*(q_seri(i,1)-q_seri(i,k))/q_seri(i,k)
            else
               ratqsc(i,k)=0.
            endif
         enddo
         enddo

!-----------------------------------------------------------------------
!  par nversion de la fonction log normale
!-----------------------------------------------------------------------
      else if (iflag_cld_th.eq.4) then
         ptconvth(:,:)=.false.
         ratqsc(:,:)=0.
         if(prt_level.ge.9) print*,'avant clouds_gno thermique'
         call clouds_gno &
         (klon,klev,q_seri,zqsat,clwcon0th,ptconvth,ratqsc,rnebcon0th)
         if(prt_level.ge.9) print*,' CLOUDS_GNO OK'
       
       endif

!   ratqs stables
!   -------------

      if (iflag_ratqs.eq.0) then

! Le cas iflag_ratqs=0 correspond a la version IPCC 2005 du modele.
         do k=1,klev
            do i=1, klon
               ratqss(i,k)=ratqsbas+(ratqshaut-ratqsbas)* &
               min((paprs(i,1)-pplay(i,k))/(paprs(i,1)-30000.),1.) 
            enddo 
         enddo

! Pour iflag_ratqs=1 ou 2, le ratqs est constant au dessus de 
! 300 hPa (ratqshaut), varie lineariement en fonction de la pression
! entre 600 et 300 hPa et est soit constant (ratqsbas) pour iflag_ratqs=1
! soit lineaire (entre 0 a la surface et ratqsbas) pour iflag_ratqs=2
! Il s'agit de differents tests dans la phase de reglage du modele
! avec thermiques.

      else if (iflag_ratqs.eq.1) then

         do k=1,klev
            do i=1, klon
               if (pplay(i,k).ge.60000.) then
                  ratqss(i,k)=ratqsbas
               else if ((pplay(i,k).ge.30000.).and.(pplay(i,k).lt.60000.)) then
                  ratqss(i,k)=ratqsbas+(ratqshaut-ratqsbas)*(60000.-pplay(i,k))/(60000.-30000.)
               else
                  ratqss(i,k)=ratqshaut
               endif
            enddo
         enddo

      else if (iflag_ratqs.eq.2) then

         do k=1,klev
            do i=1, klon
               if (pplay(i,k).ge.60000.) then
                  ratqss(i,k)=ratqsbas*(paprs(i,1)-pplay(i,k))/(paprs(i,1)-60000.)
               else if ((pplay(i,k).ge.30000.).and.(pplay(i,k).lt.60000.)) then
                    ratqss(i,k)=ratqsbas+(ratqshaut-ratqsbas)*(60000.-pplay(i,k))/(60000.-30000.)
               else
                    ratqss(i,k)=ratqshaut
               endif
            enddo
         enddo

      else if (iflag_ratqs==3) then
         do k=1,klev
           ratqss(:,k)=ratqsbas+(ratqshaut-ratqsbas) &
           *min( ((paprs(:,1)-pplay(:,k))/70000.)**2 , 1. )
         enddo

      else if (iflag_ratqs==4) then
         do k=1,klev
           ratqss(:,k)=ratqsbas+0.5*(ratqshaut-ratqsbas) &
!          *( tanh( (50000.-pplay(:,k))/20000.) + 1.)
           *( tanh( (ratqsp0-pplay(:,k))/ratqsdp) + 1.)
         enddo

      endif




!  ratqs final
!  -----------

      if (iflag_cld_th.eq.1 .or.iflag_cld_th.eq.2.or.iflag_cld_th.eq.4) then

! On ajoute une constante au ratqsc*2 pour tenir compte de 
! fluctuations turbulentes de petite echelle

         do k=1,klev
            do i=1,klon
               if ((fm_therm(i,k).gt.1.e-10)) then
                  ratqsc(i,k)=sqrt(ratqsc(i,k)**2+0.05**2)
               endif
            enddo
         enddo

!   les ratqs sont une combinaison de ratqss et ratqsc
       if(prt_level.ge.9) write(lunout,*)'PHYLMD NOUVEAU TAU_RATQS ',tau_ratqs

         if (tau_ratqs>1.e-10) then
            facteur=exp(-pdtphys/tau_ratqs)
         else
            facteur=0.
         endif
         ratqs(:,:)=ratqsc(:,:)*(1.-facteur)+ratqs(:,:)*facteur
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH 22/09/2009
! La ligne ci-dessous faisait osciller le modele et donnait une solution
! assymptotique bidon et dépendant fortement du pas de temps.
!        ratqs(:,:)=sqrt(ratqs(:,:)**2+ratqss(:,:)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ratqs(:,:)=max(ratqs(:,:),ratqss(:,:))
      else if (iflag_cld_th<=6) then
!   on ne prend que le ratqs stable pour fisrtilp
         ratqs(:,:)=ratqss(:,:)
      else
          zfratqs1=exp(-pdtphys/10800.)
          zfratqs2=exp(-pdtphys/10800.)
          do k=1,klev
             do i=1,klon
                if (ratqsc(i,k).gt.1.e-10) then
                   ratqs(i,k)=ratqs(i,k)*zfratqs2+(iflag_cld_th/100.)*ratqsc(i,k)*(1.-zfratqs2)
                endif
                ratqs(i,k)=min(ratqs(i,k)*zfratqs1+ratqss(i,k)*(1.-zfratqs1),0.5)
             enddo
          enddo
      endif


return
end
