! $Id: thermcell_main.F90 2351 2015-08-25 15:14:59Z emillour $
!
      SUBROUTINE thermcell_alp(ngrid,nlay,ptimestep  &
     &                  ,pplay,pplev  &
     &                  ,fm0,entr0,lmax  &
     &                  ,ale_bl,alp_bl,lalim_conv,wght_th &
     &                  ,zw2,fraca &
!!! ncessaire en plus
     &                  ,pcon,rhobarz,wth3,wmax_sec,lalim,fm,alim_star,zmax &
!!! nrlmd le 10/04/2012
     &                  ,pbl_tke,pctsrf,omega,airephy &
     &                  ,zlcl,fraca0,w0,w_conv,therm_tke_max0,env_tke_max0 &
     &                  ,n2,s2,ale_bl_stat &
     &                  ,therm_tke_max,env_tke_max &
     &                  ,alp_bl_det,alp_bl_fluct_m,alp_bl_fluct_tke &
     &                  ,alp_bl_conv,alp_bl_stat &
!!! fin nrlmd le 10/04/2012
     &)

      USE dimphy
      USE indice_sol_mod
      IMPLICIT NONE

!=======================================================================
!   Auteurs: Frederic Hourdin, Catherine Rio, Anne Mathieu
!   Version du 09.02.07
!   Calcul du transport vertical dans la couche limite en presence
!   de "thermiques" explicitement representes avec processus nuageux
!
!   Reecriture a partir d'un listing papier a Habas, le 14/02/00
!
!   le thermique est suppose homogene et dissipe par melange avec
!   son environnement. la longueur l_mix controle l'efficacite du
!   melange
!
!   Le calcul du transport des differentes especes se fait en prenant
!   en compte:
!     1. un flux de masse montant
!     2. un flux de masse descendant
!     3. un entrainement
!     4. un detrainement
!
! Modif 2013/01/04 (FH hourdin@lmd.jussieu.fr)
!    Introduction of an implicit computation of vertical advection in
!    the environment of thermal plumes in thermcell_dq
!    impl =     0 : explicit, 1 : implicit, -1 : old version
!    controled by iflag_thermals =
!       15, 16 run with impl=-1 : numerical convergence with NPv3
!       17, 18 run with impl=1  : more stable
!    15 and 17 correspond to the activation of the stratocumulus "bidouille"
!
!=======================================================================
!-----------------------------------------------------------------------
!   declarations:
!   -------------

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "thermcell.h"

!   arguments:
!   ----------

!IM 140508

      INTEGER ngrid,nlay
      real ptimestep
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)

!   local:
!   ------


      REAL susqr2pi, reuler

      INTEGER ig,k,l
      INTEGER lmax(klon),lalim(klon)
      real zmax(klon),zw2(klon,klev+1)

!on garde le zmax du pas de temps precedent


      real fraca(klon,klev+1)
      real wth3(klon,klev)
! FH probleme de dimensionnement avec l'allocation dynamique
!     common/comtherm/thetath2,wth2
      real rhobarz(klon,klev)

      real wmax_sec(klon)
      real fm0(klon,klev+1),entr0(klon,klev)
      real fm(klon,klev+1)

!niveau de condensation
      real pcon(klon)

      real alim_star(klon,klev)

!!! nrlmd le 10/04/2012

!------Entrées
      real pbl_tke(klon,klev+1,nbsrf)
      real pctsrf(klon,nbsrf)
      real omega(klon,klev)
      real airephy(klon)
!------Sorties
      real zlcl(klon),fraca0(klon),w0(klon),w_conv(klon)
      real therm_tke_max0(klon),env_tke_max0(klon)
      real n2(klon),s2(klon)
      real ale_bl_stat(klon)
      real therm_tke_max(klon,klev),env_tke_max(klon,klev)
      real alp_bl_det(klon),alp_bl_fluct_m(klon),alp_bl_fluct_tke(klon),alp_bl_conv(klon),alp_bl_stat(klon)
!------Local
      integer nsrf
      real rhobarz0(klon)                    ! Densité au LCL
      logical ok_lcl(klon)                   ! Existence du LCL des thermiques
      integer klcl(klon)                     ! Niveau du LCL
      real interp(klon)                      ! Coef d'interpolation pour le LCL
!--Triggering
      real Su                                ! Surface unité: celle d'un updraft élémentaire
      parameter(Su=4e4)
      real hcoef                             ! Coefficient directeur pour le calcul de s2
      parameter(hcoef=1)
      real hmincoef                          ! Coefficient directeur pour l'ordonnée à l'origine pour le calcul de s2 
      parameter(hmincoef=0.3)
      real eps1                              ! Fraction de surface occupée par la population 1 : eps1=n1*s1/(fraca0*Sd)
      parameter(eps1=0.3)
      real hmin(ngrid)                       ! Ordonnée à l'origine pour le calcul de s2
      real zmax_moy(ngrid)                   ! Hauteur moyenne des thermiques : zmax_moy = zlcl + 0.33 (zmax-zlcl)
      real zmax_moy_coef
      parameter(zmax_moy_coef=0.33)
      real depth(klon)                       ! Epaisseur moyenne du cumulus
      real w_max(klon)                       ! Vitesse max statistique 
      real s_max(klon)
!--Closure
      real pbl_tke_max(klon,klev)            ! Profil de TKE moyenne 
      real pbl_tke_max0(klon)                ! TKE moyenne au LCL
      real w_ls(klon,klev)                   ! Vitesse verticale grande échelle (m/s)
      real coef_m                            ! On considère un rendement pour alp_bl_fluct_m
      parameter(coef_m=1.)
      real coef_tke                          ! On considère un rendement pour alp_bl_fluct_tke
      parameter(coef_tke=1.)

!!! fin nrlmd le 10/04/2012

!
      !nouvelles variables pour la convection
      real ale_bl(klon)
      real alp_bl(klon)
      real alp_int(klon),dp_int(klon),zdp
      real fm_tot(klon)
      real wght_th(klon,klev)
      integer lalim_conv(klon)
!v1d     logical therm
!v1d     save therm


!------------------------------------------------------------
!  Initialize output arrays related to stochastic triggering
!------------------------------------------------------------
  DO ig = 1,klon
     zlcl(ig) = 0.
     fraca0(ig) = 0.
     w0(ig) = 0.
     w_conv(ig) = 0.
     therm_tke_max0(ig) = 0.
     env_tke_max0(ig) = 0.
     n2(ig) = 0.
     s2(ig) = 0.
     ale_bl_stat(ig) = 0.
     alp_bl_det(ig) = 0.
     alp_bl_fluct_m(ig) = 0.
     alp_bl_fluct_tke(ig) = 0.
     alp_bl_conv(ig) = 0.
     alp_bl_stat(ig) = 0.
  ENDDO
  DO l = 1,klev
    DO ig = 1,klon
     therm_tke_max(ig,l) = 0.
     env_tke_max(ig,l) = 0.
    ENDDO
  ENDDO
!------------------------------------------------------------


!------------Test sur le LCL des thermiques
    do ig=1,ngrid
      ok_lcl(ig)=.false.
      if ( (pcon(ig) .gt. pplay(ig,klev-1)) .and. (pcon(ig) .lt. pplay(ig,1)) ) ok_lcl(ig)=.true.
    enddo

!------------Localisation des niveaux entourant le LCL et du coef d'interpolation
    do l=1,nlay-1
      do ig=1,ngrid
        if (ok_lcl(ig)) then 
!ATTENTION,zw2 calcule en pplev
!          if ((pplay(ig,l) .ge. pcon(ig)) .and. (pplay(ig,l+1) .le. pcon(ig))) then
!          klcl(ig)=l
!          interp(ig)=(pcon(ig)-pplay(ig,klcl(ig)))/(pplay(ig,klcl(ig)+1)-pplay(ig,klcl(ig)))
!          endif
          if ((pplev(ig,l) .ge. pcon(ig)) .and. (pplev(ig,l+1) .le. pcon(ig))) then
          klcl(ig)=l
          interp(ig)=(pcon(ig)-pplev(ig,klcl(ig)))/(pplev(ig,klcl(ig)+1)-pplev(ig,klcl(ig)))
          endif
        endif
      enddo
    enddo

!------------Hauteur des thermiques
!!jyg le 27/04/2012
!!    do ig =1,ngrid
!!    rhobarz0(ig)=rhobarz(ig,klcl(ig))+(rhobarz(ig,klcl(ig)+1) &
!! &               -rhobarz(ig,klcl(ig)))*interp(ig)
!!    zlcl(ig)=(pplev(ig,1)-pcon(ig))/(rhobarz0(ig)*RG)
!!      if ( (.not.ok_lcl(ig)) .or. (zlcl(ig).gt.zmax(ig)) ) zlcl(ig)=zmax(ig) ! Si zclc > zmax alors on pose zlcl = zmax
!!    enddo
    do ig =1,ngrid
!CR:REHABILITATION ZMAX CONTINU
     if (ok_lcl(ig)) then 
      rhobarz0(ig)=rhobarz(ig,klcl(ig))+(rhobarz(ig,klcl(ig)+1) &
 &               -rhobarz(ig,klcl(ig)))*interp(ig)
      zlcl(ig)=(pplev(ig,1)-pcon(ig))/(rhobarz0(ig)*RG)
      zlcl(ig)=min(zlcl(ig),zmax(ig))   ! Si zlcl > zmax alors on pose zlcl = zmax
     else
      rhobarz0(ig)=0.
      zlcl(ig)=zmax(ig)
     endif
    enddo
!!jyg fin

!------------Calcul des propriétés du thermique au LCL 
  IF ( (iflag_trig_bl.ge.1) .or. (iflag_clos_bl.ge.1) ) THEN 

  !-----Initialisation de la TKE moyenne 
   do l=1,nlay
    do ig=1,ngrid
     pbl_tke_max(ig,l)=0.
    enddo
   enddo

!-----Calcul de la TKE moyenne 
   do nsrf=1,nbsrf
    do l=1,nlay
     do ig=1,ngrid
     pbl_tke_max(ig,l)=pctsrf(ig,nsrf)*pbl_tke(ig,l,nsrf)+pbl_tke_max(ig,l)
     enddo
    enddo
   enddo

!-----Initialisations des TKE dans et hors des thermiques 
   do l=1,nlay
    do ig=1,ngrid
    therm_tke_max(ig,l)=pbl_tke_max(ig,l)
    env_tke_max(ig,l)=pbl_tke_max(ig,l)
    enddo
   enddo

!-----Calcul de la TKE transportée par les thermiques : therm_tke_max
   call thermcell_tke_transport(ngrid,nlay,ptimestep,fm0,entr0,  &
  &           rg,pplev,therm_tke_max)
!   print *,' thermcell_tke_transport -> '   !!jyg

!-----Calcul des profils verticaux de TKE hors thermiques : env_tke_max, et de la vitesse verticale grande échelle : W_ls
   do l=1,nlay
    do ig=1,ngrid
     pbl_tke_max(ig,l)=fraca(ig,l)*therm_tke_max(ig,l)+(1.-fraca(ig,l))*env_tke_max(ig,l)         !  Recalcul de TKE moyenne aprés transport de TKE_TH
     env_tke_max(ig,l)=(pbl_tke_max(ig,l)-fraca(ig,l)*therm_tke_max(ig,l))/(1.-fraca(ig,l))       !  Recalcul de TKE dans  l'environnement aprés transport de TKE_TH
     w_ls(ig,l)=-1.*omega(ig,l)/(RG*rhobarz(ig,l))                                                !  Vitesse verticale de grande échelle
    enddo
   enddo
!    print *,' apres w_ls = '   !!jyg

  do ig=1,ngrid
   if (ok_lcl(ig)) then
     fraca0(ig)=fraca(ig,klcl(ig))+(fraca(ig,klcl(ig)+1) &
 &             -fraca(ig,klcl(ig)))*interp(ig)
     w0(ig)=zw2(ig,klcl(ig))+(zw2(ig,klcl(ig)+1) &
 &         -zw2(ig,klcl(ig)))*interp(ig)
     w_conv(ig)=w_ls(ig,klcl(ig))+(w_ls(ig,klcl(ig)+1) &
 &             -w_ls(ig,klcl(ig)))*interp(ig)
     therm_tke_max0(ig)=therm_tke_max(ig,klcl(ig)) &
 &                     +(therm_tke_max(ig,klcl(ig)+1)-therm_tke_max(ig,klcl(ig)))*interp(ig)
     env_tke_max0(ig)=env_tke_max(ig,klcl(ig))+(env_tke_max(ig,klcl(ig)+1) &
 &                   -env_tke_max(ig,klcl(ig)))*interp(ig)
     pbl_tke_max0(ig)=pbl_tke_max(ig,klcl(ig))+(pbl_tke_max(ig,klcl(ig)+1) &
 &                   -pbl_tke_max(ig,klcl(ig)))*interp(ig)
     if (therm_tke_max0(ig).ge.20.) therm_tke_max0(ig)=20.
     if (env_tke_max0(ig).ge.20.) env_tke_max0(ig)=20.
     if (pbl_tke_max0(ig).ge.20.) pbl_tke_max0(ig)=20.
   else 
     fraca0(ig)=0.
     w0(ig)=0.
!!jyg le 27/04/2012
!!     zlcl(ig)=0.
!!
   endif
  enddo

  ENDIF ! IF ( (iflag_trig_bl.ge.1) .or. (iflag_clos_bl.ge.1) )
!  print *,'ENDIF  ( (iflag_trig_bl.ge.1) .or. (iflag_clos_bl.ge.1) ) '    !!jyg

!------------Triggering------------------
  IF (iflag_trig_bl.ge.1) THEN 

!-----Initialisations
   depth(:)=0.
   n2(:)=0.
   s2(:)=100. ! some low value, arbitrary
   s_max(:)=0.

!-----Epaisseur du nuage (depth) et détermination de la queue du spectre de panaches (n2,s2) et du panache le plus gros (s_max)
   do ig=1,ngrid
     zmax_moy(ig)=zlcl(ig)+zmax_moy_coef*(zmax(ig)-zlcl(ig))
     depth(ig)=zmax_moy(ig)-zlcl(ig)
     hmin(ig)=hmincoef*zlcl(ig)
     if (depth(ig).ge.10.) then 
       s2(ig)=(hcoef*depth(ig)+hmin(ig))**2
       n2(ig)=(1.-eps1)*fraca0(ig)*airephy(ig)/s2(ig)
!!
!!jyg le 27/04/2012
!!       s_max(ig)=s2(ig)*log(n2(ig))
!!       if (n2(ig) .lt. 1) s_max(ig)=0.
       s_max(ig)=s2(ig)*log(max(n2(ig),1.))
!!fin jyg
     else
       n2(ig)=0.
       s_max(ig)=0.
     endif
   enddo
!   print *,'avant Calcul de Wmax '    !!jyg

!-----Calcul de Wmax et ALE_BL_STAT associée
!!jyg le 30/04/2012
!!   do ig=1,ngrid
!!     if ( (depth(ig).ge.10.) .and. (s_max(ig).gt.1.) ) then
!!     w_max(ig)=w0(ig)*(1.+sqrt(2.*log(s_max(ig)/su)-log(2.*3.14)-log(2.*log(s_max(ig)/su)-log(2.*3.14))))
!!     ale_bl_stat(ig)=0.5*w_max(ig)**2
!!     else
!!     w_max(ig)=0.
!!     ale_bl_stat(ig)=0.
!!     endif
!!   enddo
   susqr2pi=su*sqrt(2.*Rpi)
   reuler=exp(1.)
   do ig=1,ngrid
     if ( (depth(ig).ge.10.) .and. (s_max(ig).gt.susqr2pi*reuler) ) then
      w_max(ig)=w0(ig)*(1.+sqrt(2.*log(s_max(ig)/susqr2pi)-log(2.*log(s_max(ig)/susqr2pi))))
      ale_bl_stat(ig)=0.5*w_max(ig)**2
     else
      w_max(ig)=0.
      ale_bl_stat(ig)=0.
     endif
   enddo

  ENDIF ! iflag_trig_bl
!  print *,'ENDIF  iflag_trig_bl'    !!jyg

!------------Closure------------------

  IF (iflag_clos_bl.ge.2) THEN 

!-----Calcul de ALP_BL_STAT
  do ig=1,ngrid
  alp_bl_det(ig)=0.5*coef_m*rhobarz0(ig)*(w0(ig)**3)*fraca0(ig)*(1.-2.*fraca0(ig))/((1.-fraca0(ig))**2)
  alp_bl_fluct_m(ig)=1.5*rhobarz0(ig)*fraca0(ig)*(w_conv(ig)+coef_m*w0(ig))* &
 &                   (w0(ig)**2)
  alp_bl_fluct_tke(ig)=3.*coef_m*rhobarz0(ig)*w0(ig)*fraca0(ig)*(therm_tke_max0(ig)-env_tke_max0(ig)) &
 &                    +3.*rhobarz0(ig)*w_conv(ig)*pbl_tke_max0(ig)
    if (iflag_clos_bl.ge.2) then 
    alp_bl_conv(ig)=1.5*coef_m*rhobarz0(ig)*fraca0(ig)*(fraca0(ig)/(1.-fraca0(ig)))*w_conv(ig)* &
 &                   (w0(ig)**2)
    else
    alp_bl_conv(ig)=0.
    endif
  alp_bl_stat(ig)=alp_bl_det(ig)+alp_bl_fluct_m(ig)+alp_bl_fluct_tke(ig)+alp_bl_conv(ig)
  enddo

!-----Sécurité ALP infinie
  do ig=1,ngrid
   if (fraca0(ig).gt.0.98) alp_bl_stat(ig)=2.
  enddo

  ENDIF ! (iflag_clos_bl.ge.2)

!!! fin nrlmd le 10/04/2012

!      print*,'avant calcul ale et alp' 
!calcul de ALE et ALP pour la convection
      alp_bl(:)=0.
      ale_bl(:)=0.
!          print*,'ALE,ALP ,l,zw2(ig,l),ale_bl(ig),alp_bl(ig)'
      do l=1,nlay
      do ig=1,ngrid
           alp_bl(ig)=max(alp_bl(ig),0.5*rhobarz(ig,l)*wth3(ig,l) )
           ale_bl(ig)=max(ale_bl(ig),0.5*zw2(ig,l)**2)
!          print*,'ALE,ALP',l,zw2(ig,l),ale_bl(ig),alp_bl(ig)
      enddo
      enddo

! ale sec (max de wmax/2 sous la zone d'inhibition) dans
! le cas iflag_trig_bl=3
      IF (iflag_trig_bl==3) ale_bl(:)=0.5*wmax_sec(:)**2

!test:calcul de la ponderation des couches pour KE
!initialisations

      fm_tot(:)=0.
      wght_th(:,:)=1.
      lalim_conv(:)=lalim(:)

      do k=1,klev
         do ig=1,ngrid
            if (k<=lalim_conv(ig)) fm_tot(ig)=fm_tot(ig)+fm(ig,k)
         enddo
      enddo

! assez bizarre car, si on est dans la couche d'alim et que alim_star et
! plus petit que 1.e-10, on prend wght_th=1.
      do k=1,klev
         do ig=1,ngrid
            if (k<=lalim_conv(ig).and.alim_star(ig,k)>1.e-10) then
               wght_th(ig,k)=alim_star(ig,k)
            endif
         enddo
      enddo

!      print*,'apres wght_th'
!test pour prolonger la convection
      do ig=1,ngrid
!v1d  if ((alim_star(ig,1).lt.1.e-10).and.(therm)) then
      if ((alim_star(ig,1).lt.1.e-10)) then
      lalim_conv(ig)=1
      wght_th(ig,1)=1.
!      print*,'lalim_conv ok',lalim_conv(ig),wght_th(ig,1)
      endif
      enddo

!------------------------------------------------------------------------
! Modif CR/FH 20110310 : alp integree sur la verticale.
! Integrale verticale de ALP.
! wth3 etant aux niveaux inter-couches, on utilise d play comme masse des
! couches
!------------------------------------------------------------------------

      alp_int(:)=0.
      dp_int(:)=0.
      do l=2,nlay
        do ig=1,ngrid
           if(l.LE.lmax(ig)) THEN
           zdp=pplay(ig,l-1)-pplay(ig,l)
           alp_int(ig)=alp_int(ig)+0.5*rhobarz(ig,l)*wth3(ig,l)*zdp
           dp_int(ig)=dp_int(ig)+zdp
           endif
        enddo
      enddo

      if (iflag_coupl>=3 .and. iflag_coupl<=5) then
      do ig=1,ngrid
!valeur integree de alp_bl * 0.5:
        if (dp_int(ig)>0.) then
        alp_bl(ig)=alp_int(ig)/dp_int(ig)
        endif
      enddo!
      endif


! Facteur multiplicatif sur alp_bl
      alp_bl(:)=alp_bl_k*alp_bl(:)

!------------------------------------------------------------------------



      return
      end
