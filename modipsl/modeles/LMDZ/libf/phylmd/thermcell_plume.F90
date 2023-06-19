!
! $Id: thermcell_plume.F90 3074 2017-11-15 13:31:44Z fhourdin $
!
      SUBROUTINE thermcell_plume(itap,ngrid,klev,ptimestep,ztv,zthl,po,zl,rhobarz,  &
     &           zlev,pplev,pphi,zpspsk,alim_star,alim_star_tot,  &
     &           lalim,f0,detr_star,entr_star,f_star,csc,ztva,  &
     &           ztla,zqla,zqta,zha,zw2,w_est,ztva_est,zqsatth,lmix,lmix_bis,linter &
    &           ,lev_out,lunout1,igout)
!     &           ,lev_out,lunout1,igout,zbuoy,zbuoyjam)
!--------------------------------------------------------------------------
! Auhtors : Catherine Rio, Frédéric Hourdin, Arnaud Jam
!
!thermcell_plume: calcule les valeurs de qt, thetal et w dans l ascendance
!   This versions starts from a cleaning of thermcell_plume_6A (2019/01/20)
!   thermcell_plume_6A is activate for flag_thermas_ed < 10
!   thermcell_plume_5B for flag_thermas_ed < 20
!   thermcell_plume for flag_thermals_ed>= 20
!   Various options are controled by the flag_thermals_ed parameter
!   = 20 : equivalent to thermcell_plume_6A with flag_thermals_ed=8
!   = 21 : the Jam strato-cumulus modif is not activated in detrainment
!   = 29 : an other way to compute the modified buoyancy (to be tested)
!--------------------------------------------------------------------------
USE IOIPSL, ONLY : getin
USE ioipsl_getin_p_mod, ONLY : getin_p

       USE print_control_mod, ONLY: prt_level
       IMPLICIT NONE

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "thermcell.h"

      INTEGER itap
      INTEGER lunout1,igout
      INTEGER ngrid,klev
      REAL ptimestep
      REAL ztv(ngrid,klev)
      REAL zthl(ngrid,klev)
      REAL po(ngrid,klev)
      REAL zl(ngrid,klev)
      REAL rhobarz(ngrid,klev)
      REAL zlev(ngrid,klev+1)
      REAL pplev(ngrid,klev+1)
      REAL pphi(ngrid,klev)
      REAL zpspsk(ngrid,klev)
      REAL alim_star(ngrid,klev)
      REAL f0(ngrid)
      INTEGER lalim(ngrid)
      integer lev_out                           ! niveau pour les print
      integer nbpb
    
      real alim_star_tot(ngrid)

      REAL ztva(ngrid,klev)
      REAL ztla(ngrid,klev)
      REAL zqla(ngrid,klev)
      REAL zqta(ngrid,klev)
      REAL zha(ngrid,klev)

      REAL detr_star(ngrid,klev)
      REAL coefc
      REAL entr_star(ngrid,klev)
      REAL detr(ngrid,klev)
      REAL entr(ngrid,klev)

      REAL csc(ngrid,klev)

      REAL zw2(ngrid,klev+1)
      REAL w_est(ngrid,klev+1)
      REAL f_star(ngrid,klev+1)
      REAL wa_moy(ngrid,klev+1)

      REAL ztva_est(ngrid,klev)
      REAL ztv_est(ngrid,klev)
      REAL zqla_est(ngrid,klev)
      REAL zqsatth(ngrid,klev)
      REAL zta_est(ngrid,klev)
      REAL ztemp(ngrid),zqsat(ngrid)
      REAL zdw2,zdw2bis
      REAL zw2modif
      REAL zw2fact,zw2factbis
      REAL zeps(ngrid,klev)

      REAL linter(ngrid)
      INTEGER lmix(ngrid)
      INTEGER lmix_bis(ngrid)
      REAL    wmaxa(ngrid)

      INTEGER ig,l,k,lt,it,lm

      real zdz,zbuoy(ngrid,klev),zalpha,gamma(ngrid,klev),zdqt(ngrid,klev),zw2m
      real zbuoyjam(ngrid,klev),zdqtjam(ngrid,klev)
      real zdz2,zdz3,lmel,entrbis,zdzbis
      real d_temp(ngrid)
      real ztv1,ztv2,factinv,zinv,zlmel
      real zlmelup,zlmeldwn,zlt,zltdwn,zltup
      real atv1,atv2,btv1,btv2
      real ztv_est1,ztv_est2
      real zcor,zdelta,zcvm5,qlbef
      real zbetalpha, coefzlmel
      real eps
      REAL REPS,RLvCp,DDT0
      PARAMETER (DDT0=.01)
      logical Zsat
      LOGICAL active(ngrid),activetmp(ngrid)
      REAL fact_gamma,fact_gamma2,fact_epsilon2

      REAL, SAVE :: fact_epsilon=0.002
      REAL, SAVE :: betalpha=0.9
      REAL, SAVE :: afact=2./3.
      REAL, SAVE :: fact_shell=1.
      REAL,SAVE :: detr_min=1.e-5
      REAL,SAVE :: entr_min=1.e-5
      REAL,SAVE :: detr_q_coef=0.012
      REAL,SAVE :: detr_q_power=0.5
      REAL,SAVE :: mix0=0.
      INTEGER,SAVE :: thermals_flag_alim=0

!$OMP THREADPRIVATE(fact_epsilon, betalpha, afact, fact_shell)
!$OMP THREADPRIVATE(detr_min, entr_min, detr_q_coef, detr_q_power)
!$OMP THREADPRIVATE( mix0, thermals_flag_alim)

      LOGICAL, SAVE :: first=.true.
  !$OMP THREADPRIVATE(first)


      REAL c2(ngrid,klev)

      if (ngrid==1) print*,'THERMCELL PLUME MODIFIE 2014/07/11'
      Zsat=.false.
! Initialisation

      RLvCp = RLVTT/RCPD
      IF (first) THEN

     CALL getin_p('thermals_fact_epsilon',fact_epsilon)
     CALL getin_p('thermals_betalpha',betalpha)
     CALL getin_p('thermals_afact',afact)
     CALL getin_p('thermals_fact_shell',fact_shell)
     CALL getin_p('thermals_detr_min',detr_min)
     CALL getin_p('thermals_entr_min',entr_min)
     CALL getin_p('thermals_detr_q_coef',detr_q_coef)
     CALL getin_p('thermals_detr_q_power',detr_q_power)
     CALL getin_p('thermals_mix0',mix0)
     CALL getin_p('thermals_flag_alim',thermals_flag_alim)


      first=.false.
      ENDIF

      zbetalpha=betalpha/(1.+betalpha)


! Initialisations des variables r?elles
if (1==1) then
      ztva(:,:)=ztv(:,:)
      ztva_est(:,:)=ztva(:,:)
      ztv_est(:,:)=ztv(:,:)
      ztla(:,:)=zthl(:,:)
      zqta(:,:)=po(:,:)
      zqla(:,:)=0.
      zha(:,:) = ztva(:,:)
else
      ztva(:,:)=0.
      ztv_est(:,:)=0.
      ztva_est(:,:)=0.
      ztla(:,:)=0.
      zqta(:,:)=0.
      zha(:,:) =0.
endif

      zqla_est(:,:)=0.
      zqsatth(:,:)=0.
      zqla(:,:)=0.
      detr_star(:,:)=0.
      entr_star(:,:)=0.
      alim_star(:,:)=0.
      alim_star_tot(:)=0.
      csc(:,:)=0.
      detr(:,:)=0.
      entr(:,:)=0.
      zw2(:,:)=0.
      zbuoy(:,:)=0.
      zbuoyjam(:,:)=0.
      gamma(:,:)=0.
      zeps(:,:)=0.
      w_est(:,:)=0.
      f_star(:,:)=0.
      wa_moy(:,:)=0.
      linter(:)=1.
!     linter(:)=1.
! Initialisation des variables entieres
      lmix(:)=1
      lmix_bis(:)=2
      wmaxa(:)=0.


!-------------------------------------------------------------------------
! On ne considere comme actif que les colonnes dont les deux premieres
! couches sont instables.
!-------------------------------------------------------------------------

      active(:)=ztv(:,1)>ztv(:,2)
      d_temp(:)=0. ! Pour activer un contraste de temperature a la base
                   ! du panache
!  Cet appel pourrait être fait avant thermcell_plume dans thermcell_main
      CALL thermcell_alim(thermals_flag_alim,ngrid,klev,ztv,d_temp,zlev,alim_star,lalim)

!------------------------------------------------------------------------------
! Calcul dans la premiere couche
! On decide dans cette version que le thermique n'est actif que si la premiere
! couche est instable.
! Pourrait etre change si on veut que le thermiques puisse se d??clencher
! dans une couche l>1
!------------------------------------------------------------------------------
do ig=1,ngrid
! Le panache va prendre au debut les caracteristiques de l'air contenu
! dans cette couche.
    if (active(ig)) then
    ztla(ig,1)=zthl(ig,1) 
    zqta(ig,1)=po(ig,1)
    zqla(ig,1)=zl(ig,1)
!cr: attention, prise en compte de f*(1)=1
    f_star(ig,2)=alim_star(ig,1)
    zw2(ig,2)=2.*RG*(ztv(ig,1)-ztv(ig,2))/ztv(ig,2)  &
&                     *(zlev(ig,2)-zlev(ig,1))  &
&                     *0.4*pphi(ig,1)/(pphi(ig,2)-pphi(ig,1))
    w_est(ig,2)=zw2(ig,2)
    endif
enddo
!

!==============================================================================
!boucle de calcul de la vitesse verticale dans le thermique
!==============================================================================
do l=2,klev-1
!==============================================================================


! On decide si le thermique est encore actif ou non
! AFaire : Il faut sans doute ajouter entr_star a alim_star dans ce test
    do ig=1,ngrid
       active(ig)=active(ig) &
&                 .and. zw2(ig,l)>1.e-10 &
&                 .and. f_star(ig,l)+alim_star(ig,l)>1.e-10
    enddo



!---------------------------------------------------------------------------
! calcul des proprietes thermodynamiques et de la vitesse de la couche l
! sans tenir compte du detrainement et de l'entrainement dans cette
! couche
! C'est a dire qu'on suppose 
! ztla(l)=ztla(l-1) et zqta(l)=zqta(l-1)
! Ici encore, on doit pouvoir ajouter entr_star (qui peut etre calculer
! avant) a l'alimentation pour avoir un calcul plus propre
!---------------------------------------------------------------------------

   ztemp(:)=zpspsk(:,l)*ztla(:,l-1)
   call thermcell_qsat(ngrid,active,pplev(:,l),ztemp,zqta(:,l-1),zqsat(:))
    do ig=1,ngrid 
!       print*,'active',active(ig),ig,l
        if(active(ig)) then 
        zqla_est(ig,l)=max(0.,zqta(ig,l-1)-zqsat(ig))
        ztva_est(ig,l) = ztla(ig,l-1)*zpspsk(ig,l)+RLvCp*zqla_est(ig,l)
        zta_est(ig,l)=ztva_est(ig,l)
        ztva_est(ig,l) = ztva_est(ig,l)/zpspsk(ig,l)
        ztva_est(ig,l) = ztva_est(ig,l)*(1.+RETV*(zqta(ig,l-1)  &
     &      -zqla_est(ig,l))-zqla_est(ig,l))
 

!Modif AJAM

        zbuoy(ig,l)=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l) 
        zdz=zlev(ig,l+1)-zlev(ig,l)         
        lmel=fact_thermals_ed_dz*zlev(ig,l)
!        lmel=0.09*zlev(ig,l)
        zlmel=zlev(ig,l)+lmel
        zlmelup=zlmel+(zdz/2)
        zlmeldwn=zlmel-(zdz/2)

        lt=l+1
        zlt=zlev(ig,lt)
        zdz3=zlev(ig,lt+1)-zlt
        zltdwn=zlt-zdz3/2
        zltup=zlt+zdz3/2
         
!=========================================================================
! 3. Calcul de la flotabilite modifie par melange avec l'air au dessus
!=========================================================================

!--------------------------------------------------
        lt=l+1
        zlt=zlev(ig,lt)
        zdz2=zlev(ig,lt)-zlev(ig,l)

        do while (lmel.gt.zdz2)
           lt=lt+1
           zlt=zlev(ig,lt)
           zdz2=zlt-zlev(ig,l)
        enddo
        zdz3=zlev(ig,lt+1)-zlt
        zltdwn=zlev(ig,lt)-zdz3/2
        zlmelup=zlmel+(zdz/2)
        coefzlmel=Min(1.,(zlmelup-zltdwn)/zdz)
        zbuoyjam(ig,l)=1.*RG*(coefzlmel*(ztva_est(ig,l)- & 
    &   ztv(ig,lt))/ztv(ig,lt)+(1.-coefzlmel)*(ztva_est(ig,l)- &
    &   ztv(ig,lt-1))/ztv(ig,lt-1))+0.*zbuoy(ig,l)

!------------------------------------------------
!AJAM:nouveau calcul de w?  
!------------------------------------------------
        zdz=zlev(ig,l+1)-zlev(ig,l)
        zdzbis=zlev(ig,l)-zlev(ig,l-1)
        zbuoy(ig,l)=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)
        zw2fact=fact_epsilon*2.*zdz/(1.+betalpha)
        zw2factbis=fact_epsilon*2.*zdzbis/(1.+betalpha)
        zdw2=afact*zbuoy(ig,l)/fact_epsilon
        zdw2bis=afact*zbuoy(ig,l-1)/fact_epsilon
        lm=Max(1,l-2)
        w_est(ig,l+1)=Max(0.0001,exp(-zw2fact)*(w_est(ig,l)-zdw2)+zdw2)
       endif
    enddo


!-------------------------------------------------
!calcul des taux d'entrainement et de detrainement
!-------------------------------------------------

     do ig=1,ngrid
        if (active(ig)) then

!          zw2m=max(0.5*(w_est(ig,l)+w_est(ig,l+1)),0.1)
          zw2m=w_est(ig,l+1)
          zdz=zlev(ig,l+1)-zlev(ig,l)
          zbuoy(ig,l)=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)
          zalpha=f0(ig)*f_star(ig,l)/sqrt(w_est(ig,l+1))/rhobarz(ig,l)
          zdqt(ig,l)=max(zqta(ig,l-1)-po(ig,l),0.)/po(ig,l)

!=========================================================================
! 4. Calcul de l'entrainement et du detrainement
!=========================================================================

          detr_star(ig,l)=f_star(ig,l)*zdz             &
    &     *( mix0 * 0.1 / (zalpha+0.001)               &
    &     + MAX(detr_min, -afact*zbetalpha*zbuoyjam(ig,l)/zw2m   &
    &     + detr_q_coef*(zdqt(ig,l)/zw2m)**detr_q_power))

          if ( iflag_thermals_ed == 20 ) then
             entr_star(ig,l)=f_star(ig,l)*zdz* (         &
    &          mix0 * 0.1 / (zalpha+0.001)               &
    &        + zbetalpha*MAX(entr_min,                   &
    &        afact*zbuoyjam(ig,l)/zw2m - fact_epsilon))
          else
             entr_star(ig,l)=f_star(ig,l)*zdz* (         &
    &          mix0 * 0.1 / (zalpha+0.001)               &
    &        + zbetalpha*MAX(entr_min,                   &
    &        afact*zbuoy(ig,l)/zw2m - fact_epsilon))
          endif
          
! En dessous de lalim, on prend le max de alim_star et entr_star pour
! alim_star et 0 sinon
        if (l.lt.lalim(ig)) then
          alim_star(ig,l)=max(alim_star(ig,l),entr_star(ig,l))
          entr_star(ig,l)=0.
        endif
        f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)  &
     &              -detr_star(ig,l)

      endif
   enddo


!============================================================================
! 5. calcul de la vitesse verticale en melangeant Tl et qt du thermique
!===========================================================================

   activetmp(:)=active(:) .and. f_star(:,l+1)>1.e-10
   do ig=1,ngrid
       if (activetmp(ig)) then 
           Zsat=.false.
           ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+  &
     &            (alim_star(ig,l)+entr_star(ig,l))*zthl(ig,l))  &
     &            /(f_star(ig,l+1)+detr_star(ig,l))
           zqta(ig,l)=(f_star(ig,l)*zqta(ig,l-1)+  &
     &            (alim_star(ig,l)+entr_star(ig,l))*po(ig,l))  &
     &            /(f_star(ig,l+1)+detr_star(ig,l))

        endif
    enddo

   ztemp(:)=zpspsk(:,l)*ztla(:,l)
   call thermcell_qsat(ngrid,activetmp,pplev(:,l),ztemp,zqta(:,l),zqsatth(:,l))
   do ig=1,ngrid
      if (activetmp(ig)) then
! on ecrit de maniere conservative (sat ou non)
!          T = Tl +Lv/Cp ql
           zqla(ig,l)=max(0.,zqta(ig,l)-zqsatth(ig,l))
           ztva(ig,l) = ztla(ig,l)*zpspsk(ig,l)+RLvCp*zqla(ig,l)
           ztva(ig,l) = ztva(ig,l)/zpspsk(ig,l)
!on rajoute le calcul de zha pour diagnostiques (temp potentielle)
           zha(ig,l) = ztva(ig,l)
           ztva(ig,l) = ztva(ig,l)*(1.+RETV*(zqta(ig,l)  &
     &              -zqla(ig,l))-zqla(ig,l))
           zbuoy(ig,l)=RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
           zdz=zlev(ig,l+1)-zlev(ig,l)
           zdzbis=zlev(ig,l)-zlev(ig,l-1)
           zeps(ig,l)=(entr_star(ig,l)+alim_star(ig,l))/(f_star(ig,l)*zdz)
           zw2fact=fact_epsilon*2.*zdz/(1.+betalpha)
           zw2factbis=fact_epsilon*2.*zdzbis/(1.+betalpha)
           zdw2= afact*zbuoy(ig,l)/(fact_epsilon)
           zdw2bis= afact*zbuoy(ig,l-1)/(fact_epsilon)
           zw2(ig,l+1)=Max(0.0001,exp(-zw2fact)*(zw2(ig,l)-zdw2)+zdw2)
      endif
   enddo

   if (prt_level.ge.20) print*,'coucou calcul detr 460: ig, l',ig, l
!
!===========================================================================
! 6. initialisations pour le calcul de la hauteur du thermique, de l'inversion et de la vitesse verticale max 
!===========================================================================

   nbpb=0
   do ig=1,ngrid
            if (zw2(ig,l+1)>0. .and. zw2(ig,l+1).lt.1.e-10) then
!               stop'On tombe sur le cas particulier de thermcell_dry'
!               print*,'On tombe sur le cas particulier de thermcell_plume'
                nbpb=nbpb+1
                zw2(ig,l+1)=0.
                linter(ig)=l+1
            endif

        if (zw2(ig,l+1).lt.0.) then 
           linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))  &
     &               -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
           zw2(ig,l+1)=0.
!+CR:04/05/12:correction calcul linter pour calcul de zmax continu
        elseif (f_star(ig,l+1).lt.0.) then
           linter(ig)=(l*(f_star(ig,l+1)-f_star(ig,l))  &
     &               -f_star(ig,l))/(f_star(ig,l+1)-f_star(ig,l))
           zw2(ig,l+1)=0.
!fin CR:04/05/12
        endif

           wa_moy(ig,l+1)=sqrt(zw2(ig,l+1)) 

        if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
!   lmix est le niveau de la couche ou w (wa_moy) est maximum
!on rajoute le calcul de lmix_bis
            if (zqla(ig,l).lt.1.e-10) then
               lmix_bis(ig)=l+1
            endif
            lmix(ig)=l+1
            wmaxa(ig)=wa_moy(ig,l+1)
        endif
   enddo

   if (nbpb>0) then
   print*,'WARNING on tombe ',nbpb,' x sur un pb pour l=',l,' dans thermcell_plume'
   endif

!=========================================================================
! FIN DE LA BOUCLE VERTICALE
      enddo
!=========================================================================

!on recalcule alim_star_tot
       do ig=1,ngrid
          alim_star_tot(ig)=0.
       enddo
       do ig=1,ngrid
          do l=1,lalim(ig)-1
          alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
          enddo
       enddo
       

        if (prt_level.ge.20) print*,'coucou calcul detr 470: ig, l', ig, l

#undef wrgrads_thermcell
#ifdef wrgrads_thermcell
         call wrgradsfi(1,klev,entr_star(igout,1:klev),'esta      ','esta      ')
         call wrgradsfi(1,klev,detr_star(igout,1:klev),'dsta      ','dsta      ')
         call wrgradsfi(1,klev,zbuoy(igout,1:klev),'buoy      ','buoy      ')
         call wrgradsfi(1,klev,zdqt(igout,1:klev),'dqt      ','dqt      ')
         call wrgradsfi(1,klev,w_est(igout,1:klev),'w_est     ','w_est     ')
         call wrgradsfi(1,klev,w_est(igout,2:klev+1),'w_es2     ','w_es2     ')
         call wrgradsfi(1,klev,zw2(igout,1:klev),'zw2A      ','zw2A      ')
#endif


     return 
     end
