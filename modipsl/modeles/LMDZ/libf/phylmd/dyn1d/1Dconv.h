!
! $Id: 1Dconv.h 2310 2015-06-24 10:04:31Z fairhead $
!
        subroutine get_uvd(itap,dtime,file_forctl,file_fordat,                  &
     &       ht,hq,hw,hu,hv,hthturb,hqturb,                                     &
     &       Ts,imp_fcg,ts_fcg,Tp_fcg,Turb_fcg)                                 
!
        implicit none
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cette routine permet d obtenir u_convg,v_convg,ht,hq et ainsi de
! pouvoir calculer la convergence et le cisaillement dans la physiq
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#include "YOMCST.h"

      INTEGER klev
      REAL play(100)  !pression en Pa au milieu de chaque couche GCM
      INTEGER JM(100) !pression en Pa au milieu de chaque couche GCM
      REAL coef1(100) !coefficient d interpolation
      REAL coef2(100) !coefficient d interpolation

      INTEGER nblvlm !nombre de niveau de pression du mesoNH
      REAL playm(100)  !pression en Pa au milieu de chaque couche Meso-NH
      REAL hplaym(100) !pression en hPa milieux des couches Meso-NH

      integer i,j,k,ll,in

      CHARACTER*80 file_forctl,file_fordat

      COMMON/com1_phys_gcss/play,coef1,coef2,JM,klev
      COMMON/com2_phys_gcss/playm,hplaym,nblvlm

!======================================================================
! methode: on va chercher les donnees du mesoNH de meteo france, on y
!          a acces a tout pas detemps grace a la routine rdgrads qui
!          est une boucle lisant dans ces fichiers.
!          Puis on interpole ces donnes sur les 11 niveaux du gcm et
!          et sur les pas de temps de ce meme gcm
!----------------------------------------------------------------------
! input:
!       pasmax     :nombre de pas de temps maximum du mesoNH
!       dt         :pas de temps du meso_NH (en secondes)
!----------------------------------------------------------------------
      integer pasmax,dt
      save pasmax,dt
!----------------------------------------------------------------------
! arguments:
!           itap   :compteur de la physique(le nombre de ces pas est
!                   fixe dans la subroutine calcul_ini_gcm de interpo
!                   -lation
!           dtime  :pas detemps du gcm (en secondes)
!           ht     :convergence horizontale de temperature(K/s)
!           hq     :    "         "       d humidite (kg/kg/s)
!           hw     :vitesse verticale moyenne (m/s**2)
!           hu     :convergence horizontale d impulsion le long de x
!                  (kg/(m^2 s^2)
!           hv     : idem le long de y.
!           Ts     : Temperature de surface (K)
!           imp_fcg: var. logical .eq. T si forcage en impulsion
!           ts_fcg: var. logical .eq. T si forcage en Ts present dans fichier
!           Tp_fcg: var. logical .eq. T si forcage donne en Temp potentielle
!           Turb_fcg: var. logical .eq. T si forcage turbulent present dans fichier
!----------------------------------------------------------------------
        integer itap
        real dtime
        real ht(100)
        real hq(100)
        real hu(100)
        real hv(100)
        real hw(100)
        real hthturb(100)
        real hqturb(100)
        real Ts, Ts_subr
        logical imp_fcg
        logical ts_fcg
        logical Tp_fcg
        logical Turb_fcg
!----------------------------------------------------------------------
! Variables internes de get_uvd (note : l interpolation temporelle
! est faite entre les pas de temps before et after, sur les variables
! definies sur la grille du SCM; on atteint exactement les valeurs Meso
! aux milieux des pas de temps Meso)
!     time0     :date initiale en secondes
!     time      :temps associe a chaque pas du SCM
!     pas       :numero du pas du meso_NH (on lit en pas : le premier pas
!                 des donnees est duplique)
!     pasprev   :numero du pas de lecture precedent
!     htaft     :advection horizontale de temp. au pas de temps after
!     hqaft     :    "         "      d humidite        "
!     hwaft     :vitesse verticalle moyenne  au pas de temps after
!     huaft,hvaft :advection horizontale d impulsion au pas de temps after
!     tsaft     : surface temperature 'after time step'
!     htbef     :idem htaft, mais pour le pas de temps before
!     hqbef     :voir hqaft
!     hwbef     :voir hwaft
!     hubef,hvbef : idem huaft,hvaft, mais pour before
!     tsbef     : surface temperature 'before time step'
!----------------------------------------------------------------------
        integer time0,pas,pasprev
        save time0,pas,pasprev
        real time
        real htaft(100),hqaft(100),hwaft(100),huaft(100),hvaft(100)
        real hthturbaft(100),hqturbaft(100)
        real Tsaft
        save htaft,hqaft,hwaft,huaft,hvaft,hthturbaft,hqturbaft
        real htbef(100),hqbef(100),hwbef(100),hubef(100),hvbef(100)
        real hthturbbef(100),hqturbbef(100)
        real Tsbef
        save htbef,hqbef,hwbef,hubef,hvbef,hthturbbef,hqturbbef
!
        real timeaft,timebef
        save timeaft,timebef
        integer temps
        character*4 string
!----------------------------------------------------------------------
! variables arguments de la subroutine rdgrads
!---------------------------------------------------------------------
        integer icompt,icomp1 !compteurs de rdgrads
        real z(100)         ! altitude (grille Meso)
        real ht_mes(100)    !convergence horizontale de temperature
                            !-(grille Meso)
        real hq_mes(100)    !convergence horizontale d humidite
                            !(grille Meso)
        real hw_mes(100)    !vitesse verticale moyenne
                            !(grille Meso)
        real hu_mes(100),hv_mes(100)    !convergence horizontale d impulsion
                                        !(grille Meso)
        real hthturb_mes(100) !tendance horizontale de T_pot, due aux
                              !flux turbulents
        real hqturb_mes(100) !tendance horizontale d humidite, due aux
                              !flux turbulents
!
!---------------------------------------------------------------------
! variable argument de la subroutine copie
!---------------------------------------------------------------------
! SB        real pplay(100)    !pression en milieu de couche du gcm
! SB                            !argument de la physique
!---------------------------------------------------------------------
! variables destinees a la lecture du pas de temps du fichier de donnees
!---------------------------------------------------------------------
       character*80 aaa,atemps,spaces,apasmax
       integer nch,imn,ipa
!---------------------------------------------------------------------
!  procedures appelees
        external rdgrads    !lire en iterant dans forcing.dat
!---------------------------------------------------------------------
               print*,'le pas itap est:',itap
!*** on determine le pas du meso_NH correspondant au nouvel itap ***
!*** pour aller chercher les champs dans rdgrads                 ***
!
        time=time0+itap*dtime
!c        temps=int(time/dt+1)
!c        pas=min(temps,pasmax)
        temps = 1 + int((dt + 2*time)/(2*dt))
        pas=min(temps,pasmax-1)
             print*,'le pas Meso est:',pas
!
!
!===================================================================
!
!*** on remplit les champs before avec les champs after du pas   ***
!*** precedent en format gcm                                     ***
        if(pas.gt.pasprev)then
          do i=1,klev
             htbef(i)=htaft(i)
             hqbef(i)=hqaft(i)
             hwbef(i)=hwaft(i)
             hubef(i)=huaft(i)
             hvbef(i)=hvaft(i)
             hThTurbbef(i)=hThTurbaft(i)
             hqTurbbef(i)=hqTurbaft(i)
          enddo
          tsbef = tsaft
          timebef=pasprev*dt
          timeaft=timebef+dt
          icomp1 = nblvlm*4
          IF (ts_fcg) icomp1 = icomp1 + 1
          IF (imp_fcg) icomp1 = icomp1 + nblvlm*2
          IF (Turb_fcg) icomp1 = icomp1 + nblvlm*2
          icompt = icomp1*pas
         print *, 'imp_fcg,ts_fcg,Turb_fcg,pas,nblvlm,icompt'
         print *, imp_fcg,ts_fcg,Turb_fcg,pas,nblvlm,icompt
                       print*,'le pas pas est:',pas
!*** on va chercher les nouveaux champs after dans toga.dat     ***
!*** champs en format meso_NH                                   ***
          open(99,FILE=file_fordat,FORM='UNFORMATTED',                        &
     &             ACCESS='DIRECT',RECL=8)
          call rdgrads(99,icompt,nblvlm,z,ht_mes,hq_mes,hw_mes                &
     &                  ,hu_mes,hv_mes,hthturb_mes,hqturb_mes                 &
     &                  ,ts_fcg,ts_subr,imp_fcg,Turb_fcg)
!

               if(Tp_fcg) then
!     (le forcage est donne en temperature potentielle)
         do i = 1,nblvlm
           ht_mes(i) = ht_mes(i)*(hplaym(i)/1000.)**rkappa
         enddo
               endif ! Tp_fcg
        if(Turb_fcg) then
         do i = 1,nblvlm
           hThTurb_mes(i) = hThTurb_mes(i)*(hplaym(i)/1000.)**rkappa
         enddo
        endif  ! Turb_fcg
!
               print*,'ht_mes ',(ht_mes(i),i=1,nblvlm)
               print*,'hq_mes ',(hq_mes(i),i=1,nblvlm)
               print*,'hw_mes ',(hw_mes(i),i=1,nblvlm)
                  if(imp_fcg) then
               print*,'hu_mes ',(hu_mes(i),i=1,nblvlm)
               print*,'hv_mes ',(hv_mes(i),i=1,nblvlm)
                  endif
                  if(Turb_fcg) then
               print*,'hThTurb_mes ',(hThTurb_mes(i),i=1,nblvlm)
               print*,'hqTurb_mes ',(hqTurb_mes(i),i=1,nblvlm)
                  endif
          IF (ts_fcg) print*,'ts_subr', ts_subr
!*** on interpole les champs meso_NH sur les niveaux de pression***
!*** gcm . on obtient le nouveau champ after                    ***
            do k=1,klev
             if (JM(k) .eq. 0) then
         htaft(k)=              ht_mes(jm(k)+1)
         hqaft(k)=              hq_mes(jm(k)+1)
         hwaft(k)=              hw_mes(jm(k)+1)
               if(imp_fcg) then
           huaft(k)=              hu_mes(jm(k)+1)
           hvaft(k)=              hv_mes(jm(k)+1)
               endif ! imp_fcg
               if(Turb_fcg) then
           hThTurbaft(k)=         hThTurb_mes(jm(k)+1)
           hqTurbaft(k)=          hqTurb_mes(jm(k)+1)
               endif ! Turb_fcg
             else ! JM(k) .eq. 0
           htaft(k)=coef1(k)*ht_mes(jm(k))+coef2(k)*ht_mes(jm(k)+1)
           hqaft(k)=coef1(k)*hq_mes(jm(k))+coef2(k)*hq_mes(jm(k)+1)
           hwaft(k)=coef1(k)*hw_mes(jm(k))+coef2(k)*hw_mes(jm(k)+1)
               if(imp_fcg) then
           huaft(k)=coef1(k)*hu_mes(jm(k))+coef2(k)*hu_mes(jm(k)+1)
           hvaft(k)=coef1(k)*hv_mes(jm(k))+coef2(k)*hv_mes(jm(k)+1)
               endif ! imp_fcg
               if(Turb_fcg) then
           hThTurbaft(k)=coef1(k)*hThTurb_mes(jm(k))                            &
     &               +coef2(k)*hThTurb_mes(jm(k)+1)
           hqTurbaft(k) =coef1(k)*hqTurb_mes(jm(k))                             &
     &               +coef2(k)*hqTurb_mes(jm(k)+1)
               endif ! Turb_fcg
             endif ! JM(k) .eq. 0
            enddo
            tsaft = ts_subr
            pasprev=pas
         else ! pas.gt.pasprev
            print*,'timebef est:',timebef
         endif ! pas.gt.pasprev    fin du bloc relatif au passage au pas
                                  !de temps (meso) suivant
!*** si on atteint le pas max des donnees experimentales ,on     ***
!*** on conserve les derniers champs calcules                    ***
      if(temps.ge.pasmax)then
          do ll=1,klev
               ht(ll)=htaft(ll)
               hq(ll)=hqaft(ll)
               hw(ll)=hwaft(ll)
               hu(ll)=huaft(ll)
               hv(ll)=hvaft(ll)
               hThTurb(ll)=hThTurbaft(ll)
               hqTurb(ll)=hqTurbaft(ll)
          enddo
          ts_subr = tsaft
      else ! temps.ge.pasmax
!*** on interpole sur les pas de temps de 10mn du gcm a partir   ***
!** des pas de temps de 1h du meso_NH                            ***
         do j=1,klev
         ht(j)=((timeaft-time)*htbef(j)+(time-timebef)*htaft(j))/dt
         hq(j)=((timeaft-time)*hqbef(j)+(time-timebef)*hqaft(j))/dt
         hw(j)=((timeaft-time)*hwbef(j)+(time-timebef)*hwaft(j))/dt
             if(imp_fcg) then
         hu(j)=((timeaft-time)*hubef(j)+(time-timebef)*huaft(j))/dt
         hv(j)=((timeaft-time)*hvbef(j)+(time-timebef)*hvaft(j))/dt
             endif ! imp_fcg
             if(Turb_fcg) then
         hThTurb(j)=((timeaft-time)*hThTurbbef(j)                           &
     &           +(time-timebef)*hThTurbaft(j))/dt
         hqTurb(j)= ((timeaft-time)*hqTurbbef(j)                            &
     &           +(time-timebef)*hqTurbaft(j))/dt
             endif ! Turb_fcg
         enddo
         ts_subr = ((timeaft-time)*tsbef + (time-timebef)*tsaft)/dt
       endif ! temps.ge.pasmax
!
        print *,' time,timebef,timeaft',time,timebef,timeaft
        print *,' ht,htbef,htaft,hthturb,hthturbbef,hthturbaft'
        do j= 1,klev
           print *, j,ht(j),htbef(j),htaft(j),                              &
     &             hthturb(j),hthturbbef(j),hthturbaft(j)
        enddo
        print *,' hq,hqbef,hqaft,hqturb,hqturbbef,hqturbaft'
        do j= 1,klev
           print *, j,hq(j),hqbef(j),hqaft(j),                              &
     &             hqturb(j),hqturbbef(j),hqturbaft(j)
        enddo
!
!-------------------------------------------------------------------
!
         IF (Ts_fcg) Ts = Ts_subr
         return
!
!-----------------------------------------------------------------------
! on sort les champs de "convergence" pour l instant initial 'in'
! ceci se passe au pas temps itap=0 de la physique
!-----------------------------------------------------------------------
        entry get_uvd2(itap,dtime,file_forctl,file_fordat,                  &
     &           ht,hq,hw,hu,hv,hThTurb,hqTurb,ts,                          &
     &           imp_fcg,ts_fcg,Tp_fcg,Turb_fcg)
             print*,'le pas itap est:',itap
!
!===================================================================
!
       write(*,'(a)') 'OPEN '//file_forctl
       open(97,FILE=file_forctl,FORM='FORMATTED')
!
!------------------
      do i=1,1000
      read(97,1000,end=999) string
 1000 format (a4)
      if (string .eq. 'TDEF') go to 50
      enddo
 50   backspace(97)
!-------------------------------------------------------------------
!   *** on lit le pas de temps dans le fichier de donnees ***
!   *** "forcing.ctl" et pasmax                           ***
!-------------------------------------------------------------------
      read(97,2000) aaa
 2000  format (a80)
         print*,'aaa est',aaa
      aaa=spaces(aaa,1)
         print*,'aaa',aaa
      call getsch(aaa,' ',' ',5,atemps,nch)
         print*,'atemps est',atemps
        atemps=atemps(1:nch-2)
         print*,'atemps',atemps
        read(atemps,*) imn
        dt=imn*60
         print*,'le pas de temps dt',dt
      call getsch(aaa,' ',' ',2,apasmax,nch)
        apasmax=apasmax(1:nch)
        read(apasmax,*) ipa
        pasmax=ipa
         print*,'pasmax est',pasmax
      CLOSE(97)
!------------------------------------------------------------------
! *** on lit le pas de temps initial de la simulation ***
!------------------------------------------------------------------
                  in=itap
!c                  pasprev=in
!c                  time0=dt*(pasprev-1)
                  pasprev=in-1
                  time0=dt*pasprev
!
          close(98)
!
      write(*,'(a)') 'OPEN '//file_fordat
      open(99,FILE=file_fordat,FORM='UNFORMATTED',                          &
     &          ACCESS='DIRECT',RECL=8)
          icomp1 = nblvlm*4
          IF (ts_fcg) icomp1 = icomp1 + 1
          IF (imp_fcg) icomp1 = icomp1 + nblvlm*2
          IF (Turb_fcg) icomp1 = icomp1 + nblvlm*2
          icompt = icomp1*(in-1)
          call rdgrads(99,icompt,nblvlm,z,ht_mes,hq_mes,hw_mes              &
     &                   ,hu_mes,hv_mes,hthturb_mes,hqturb_mes              &
     &                   ,ts_fcg,ts_subr,imp_fcg,Turb_fcg)
          print *, 'get_uvd : rdgrads ->'
          print *, tp_fcg
!
! following commented out because we have temperature already in ARM case
!   (otherwise this is the potential temperature )
!------------------------------------------------------------------------
               if(Tp_fcg) then
          do i = 1,nblvlm
            ht_mes(i) = ht_mes(i)*(hplaym(i)/1000.)**rkappa
          enddo
               endif ! Tp_fcg
           print*,'ht_mes ',(ht_mes(i),i=1,nblvlm)
           print*,'hq_mes ',(hq_mes(i),i=1,nblvlm)
           print*,'hw_mes ',(hw_mes(i),i=1,nblvlm)
              if(imp_fcg) then
           print*,'hu_mes ',(hu_mes(i),i=1,nblvlm)
           print*,'hv_mes ',(hv_mes(i),i=1,nblvlm)
           print*,'t',ts_subr
              endif ! imp_fcg
              if(Turb_fcg) then
                 print*,'hThTurb_mes ',(hThTurb_mes(i),i=1,nblvlm)
                 print*,'hqTurb ',     (hqTurb_mes(i),i=1,nblvlm)
              endif ! Turb_fcg
!----------------------------------------------------------------------
! on a obtenu des champs initiaux sur les niveaux du meso_NH
! on interpole sur les niveaux du gcm(niveau pression bien sur!)
!-----------------------------------------------------------------------
            do k=1,klev
             if (JM(k) .eq. 0) then
!FKC bug? ne faut il pas convertir tsol en tendance ????
!RT bug taken care of by removing the stuff
           htaft(k)=              ht_mes(jm(k)+1)
           hqaft(k)=              hq_mes(jm(k)+1)
           hwaft(k)=              hw_mes(jm(k)+1)
               if(imp_fcg) then
           huaft(k)=              hu_mes(jm(k)+1)
           hvaft(k)=              hv_mes(jm(k)+1)
               endif ! imp_fcg
               if(Turb_fcg) then
           hThTurbaft(k)=         hThTurb_mes(jm(k)+1)
           hqTurbaft(k)=          hqTurb_mes(jm(k)+1)
               endif ! Turb_fcg
             else ! JM(k) .eq. 0
           htaft(k)=coef1(k)*ht_mes(jm(k))+coef2(k)*ht_mes(jm(k)+1)
           hqaft(k)=coef1(k)*hq_mes(jm(k))+coef2(k)*hq_mes(jm(k)+1)
           hwaft(k)=coef1(k)*hw_mes(jm(k))+coef2(k)*hw_mes(jm(k)+1)
               if(imp_fcg) then
           huaft(k)=coef1(k)*hu_mes(jm(k))+coef2(k)*hu_mes(jm(k)+1)
           hvaft(k)=coef1(k)*hv_mes(jm(k))+coef2(k)*hv_mes(jm(k)+1)
               endif ! imp_fcg
               if(Turb_fcg) then
           hThTurbaft(k)=coef1(k)*hThTurb_mes(jm(k))                        &
     &                  +coef2(k)*hThTurb_mes(jm(k)+1)
           hqTurbaft(k) =coef1(k)*hqTurb_mes(jm(k))                         &
     &                  +coef2(k)*hqTurb_mes(jm(k)+1)
               endif ! Turb_fcg
             endif ! JM(k) .eq. 0
            enddo
            tsaft = ts_subr
! valeurs initiales des champs de convergence
          do k=1,klev
             ht(k)=htaft(k)
             hq(k)=hqaft(k)
             hw(k)=hwaft(k)
                if(imp_fcg) then
             hu(k)=huaft(k)
             hv(k)=hvaft(k)
                endif ! imp_fcg
                if(Turb_fcg) then
             hThTurb(k)=hThTurbaft(k)
             hqTurb(k) =hqTurbaft(k)
                endif ! Turb_fcg
          enddo
          ts_subr = tsaft
          close(99)
          close(98)
!
!-------------------------------------------------------------------
!
!
 100      IF (Ts_fcg) Ts = Ts_subr
        return
!
999     continue
        stop 'erreur lecture, file forcing.ctl'
        end

      SUBROUTINE advect_tvl(dtime,zt,zq,vu_f,vv_f,t_f,q_f                   &
     &                     ,d_t_adv,d_q_adv)
      use dimphy
      implicit none

#include "dimensions.h"
!cccc#include "dimphy.h"

      integer k
      real dtime, fact, du, dv, cx, cy, alx, aly
      real zt(klev), zq(klev,3)
      real vu_f(klev), vv_f(klev), t_f(klev), q_f(klev,3)

      real d_t_adv(klev), d_q_adv(klev,3)

! Velocity of moving cell
      data cx,cy /12., -2./

! Dimensions of moving cell
      data alx,aly /100000.,150000./

      do k = 1, klev
            du = abs(vu_f(k)-cx)/alx
            dv = abs(vv_f(k)-cy)/aly
            fact = dtime *(du+dv-du*dv*dtime)
            d_t_adv(k) = fact * (t_f(k)-zt(k))
            d_q_adv(k,1) = fact * (q_f(k,1)-zq(k,1))
            d_q_adv(k,2) = fact * (q_f(k,2)-zq(k,2))
            d_q_adv(k,3) = fact * (q_f(k,3)-zq(k,3))
      enddo

      return
      end

      SUBROUTINE copie(klevgcm,playgcm,psolgcm,file_forctl)
      implicit none

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cette routine remplit les COMMON com1_phys_gcss et com2_phys_gcss.h
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER klev !nombre de niveau de pression du GCM
      REAL play(100)  !pression en Pa au milieu de chaque couche GCM
      INTEGER JM(100)
      REAL coef1(100)   !coefficient d interpolation
      REAL coef2(100)   !coefficient d interpolation

      INTEGER nblvlm !nombre de niveau de pression du mesoNH
      REAL playm(100) !pression en Pa au milieu de chaque couche Meso-NH
      REAL hplaym(100)!pression en hecto-Pa des milieux de couche Meso-NH

      COMMON/com1_phys_gcss/play,coef1,coef2,JM,klev
      COMMON/com2_phys_gcss/playm,hplaym,nblvlm

      integer k,klevgcm
      real playgcm(klevgcm) ! pression en milieu de couche du gcm
      real psolgcm
      character*80 file_forctl

      klev = klevgcm

!---------------------------------------------------------------------
! pression au milieu des couches du gcm dans la physiq
! (SB: remplace le call conv_lipress_gcm(playgcm) )
!---------------------------------------------------------------------

       do k = 1, klev
        play(k) = playgcm(k)
        print*,'la pression gcm est:',play(k)
       enddo

!----------------------------------------------------------------------
! lecture du descripteur des donnees Meso-NH (forcing.ctl):
!  -> nb niveaux du meso.NH (nblvlm) + pressions meso.NH
! (on remplit le COMMON com2_phys_gcss)
!----------------------------------------------------------------------

      call mesolupbis(file_forctl)

      print*,'la valeur de nblvlm est:',nblvlm

!----------------------------------------------------------------------
! etude de la correspondance entre les niveaux meso.NH et GCM;
! calcul des coefficients d interpolation coef1 et coef2
! (on remplit le COMMON com1_phys_gcss)
!----------------------------------------------------------------------

      call corresbis(psolgcm)

!---------------------------------------------------------
! TEST sur le remplissage de com1_phys_gcss et com2_phys_gcss:
!---------------------------------------------------------
 
      write(*,*) ' '
      write(*,*) 'TESTS com1_phys_gcss et com2_phys_gcss dans copie.F'
      write(*,*) '--------------------------------------'
      write(*,*) 'GCM: nb niveaux:',klev,' et pression, coeffs:'
      do k = 1, klev
      write(*,*) play(k), coef1(k), coef2(k)
      enddo
      write(*,*) 'MESO-NH: nb niveaux:',nblvlm,' et pression:'
      do k = 1, nblvlm
      write(*,*) playm(k), hplaym(k)
      enddo
      write(*,*) ' '

      end
      SUBROUTINE mesolupbis(file_forctl)
      implicit none 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Lecture descripteur des donnees MESO-NH (forcing.ctl):
! -------------------------------------------------------
!
!     Cette subroutine lit dans le fichier de controle "essai.ctl"
!     et affiche le nombre de niveaux du Meso-NH ainsi que les valeurs
!     des pressions en milieu de couche du Meso-NH (en Pa puis en hPa).
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      INTEGER nblvlm !nombre de niveau de pression du mesoNH
      REAL playm(100)  !pression en Pa milieu de chaque couche Meso-NH
      REAL hplaym(100) !pression en hPa des milieux de couche Meso-NH
      COMMON/com2_phys_gcss/playm,hplaym,nblvlm

      INTEGER i,lu,mlz,mlzh
 
      character*80 file_forctl

      character*4 a
      character*80 aaa,anblvl,spaces
      integer nch

      lu=9
      open(lu,file=file_forctl,form='formatted')
!
      do i=1,1000
      read(lu,1000,end=999) a
      if (a .eq. 'ZDEF') go to 100
      enddo
!
 100  backspace(lu)
      print*,'  DESCRIPTION DES 2 MODELES : '
      print*,' '
!
      read(lu,2000) aaa
 2000  format (a80)
       aaa=spaces(aaa,1)
       call getsch(aaa,' ',' ',2,anblvl,nch)
         read(anblvl,*) nblvlm

!
      print*,'nbre de niveaux de pression Meso-NH :',nblvlm
      print*,' '
      print*,'pression en Pa de chaque couche du meso-NH :'
!
      read(lu,*) (playm(mlz),mlz=1,nblvlm)
!      Si la pression est en HPa, la multiplier par 100
      if (playm(1) .lt. 10000.) then
        do mlz = 1,nblvlm
         playm(mlz) = playm(mlz)*100.
        enddo
      endif
      print*,(playm(mlz),mlz=1,nblvlm)
!
 1000 format (a4)
 1001 format(5x,i2)
!
      print*,' '
      do mlzh=1,nblvlm
      hplaym(mlzh)=playm(mlzh)/100.
      enddo
!
      print*,'pression en hPa de chaque couche du meso-NH: '
      print*,(hplaym(mlzh),mlzh=1,nblvlm)
!
      close (lu)
      return
!
 999  stop 'erreur lecture des niveaux pression des donnees'
      end

      SUBROUTINE rdgrads(itape,icount,nl,z,ht,hq,hw,hu,hv,hthtur,hqtur,     &
     &  ts_fcg,ts,imp_fcg,Turb_fcg)
      IMPLICIT none
      INTEGER itape,icount,icomp, nl
      real z(nl),ht(nl),hq(nl),hw(nl),hu(nl),hv(nl)
      real hthtur(nl),hqtur(nl)
      real ts
!
      INTEGER k
!
      LOGICAL imp_fcg,ts_fcg,Turb_fcg
!
      icomp = icount
!
!
         do k=1,nl
            icomp=icomp+1
            read(itape,rec=icomp)z(k)
            print *,'icomp,k,z(k) ',icomp,k,z(k)
         enddo
         do k=1,nl
            icomp=icomp+1
            read(itape,rec=icomp)hT(k)
             print*, hT(k), k
         enddo
         do k=1,nl
            icomp=icomp+1
            read(itape,rec=icomp)hQ(k)
         enddo
!
             if(turb_fcg) then
         do k=1,nl
            icomp=icomp+1
           read(itape,rec=icomp)hThTur(k)
         enddo
         do k=1,nl
            icomp=icomp+1
           read(itape,rec=icomp)hqTur(k)
         enddo
             endif
         print *,' apres lecture hthtur, hqtur'
!
          if(imp_fcg) then

         do k=1,nl
            icomp=icomp+1
           read(itape,rec=icomp)hu(k)
         enddo
         do k=1,nl
            icomp=icomp+1
            read(itape,rec=icomp)hv(k)
         enddo

          endif
!
         do k=1,nl
            icomp=icomp+1
            read(itape,rec=icomp)hw(k)
         enddo
!
              if(ts_fcg) then
         icomp=icomp+1
         read(itape,rec=icomp)ts
              endif
!
      print *,' rdgrads ->'

      RETURN
      END

      SUBROUTINE corresbis(psol)
      implicit none

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Cette subroutine calcule et affiche les valeurs des coefficients
! d interpolation qui serviront dans la formule d interpolation elle-
! meme.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER klev    !nombre de niveau de pression du GCM
      REAL play(100)  !pression en Pa au milieu de chaque couche GCM
      INTEGER JM(100)
      REAL coef1(100) !coefficient d interpolation
      REAL coef2(100) !coefficient d interpolation

      INTEGER nblvlm !nombre de niveau de pression du mesoNH
      REAL playm(100) !pression en Pa milieu de chaque couche Meso-NH
      REAL hplaym(100)!pression en hPa des milieux de couche Meso-NH

      COMMON/com1_phys_gcss/play,coef1,coef2,JM,klev
      COMMON/com2_phys_gcss/playm,hplaym,nblvlm

      REAL psol
      REAL val
      INTEGER k, mlz


      do k=1,klev
         val=play(k)
       if (val .gt. playm(1)) then
          mlz = 0
          JM(1) = mlz
          coef1(1)=(playm(mlz+1)-val)/(playm(mlz+1)-psol)
          coef2(1)=(val-psol)/(playm(mlz+1)-psol)
       else if (val .gt. playm(nblvlm)) then
         do mlz=1,nblvlm
          if (     val .le. playm(mlz).and. val .gt. playm(mlz+1))then
           JM(k)=mlz
           coef1(k)=(playm(mlz+1)-val)/(playm(mlz+1)-playm(mlz))
           coef2(k)=(val-playm(mlz))/(playm(mlz+1)-playm(mlz))
          endif
         enddo
       else
         JM(k) = nblvlm-1
         coef1(k) = 0.
         coef2(k) = 0.
       endif
      enddo
!
!c      if (play(klev) .le. playm(nblvlm)) then
!c         mlz=nblvlm-1
!c         JM(klev)=mlz
!c         coef1(klev)=(playm(mlz+1)-val)
!c     *            /(playm(mlz+1)-playm(mlz))
!c         coef2(klev)=(val-playm(mlz))
!c     *            /(playm(mlz+1)-playm(mlz))
!c      endif
!
      print*,' '
      print*,'         INTERPOLATION  : '
      print*,' '
      print*,'correspondance de 9 niveaux du GCM sur les 53 du meso-NH:'
      print*,(JM(k),k=1,klev)
      print*,'correspondance de 9 niveaux du GCM sur les 53 du meso-NH:'
      print*,(JM(k),k=1,klev)
      print*,' '
      print*,'vals du premier coef d"interpolation pour les 9 niveaux: '
      print*,(coef1(k),k=1,klev)
      print*,' '
      print*,'valeurs du deuxieme coef d"interpolation pour les 9 niveaux:'
      print*,(coef2(k),k=1,klev)
!
      return
      end
      SUBROUTINE GETSCH(STR,DEL,TRM,NTH,SST,NCH)
!***************************************************************
!*                                                             *
!*                                                             *
!* GETSCH                                                      *
!*                                                             *
!*                                                             *
!* modified by :                                               *
!***************************************************************
!*   Return in SST the character string found between the NTH-1 and NTH
!*   occurence of the delimiter 'DEL' but before the terminator 'TRM' in
!*   the input string 'STR'. If TRM=DEL then STR is considered unlimited.
!*   NCH=Length of the string returned in SST or =-1 if NTH is <1 or if
!*   NTH is greater than the number of delimiters in STR.
      IMPLICIT INTEGER (A-Z)
      CHARACTER STR*(*),DEL*1,TRM*1,SST*(*)
      NCH=-1
      SST=' '
      IF(NTH.GT.0) THEN
        IF(TRM.EQ.DEL) THEN
          LENGTH=LEN(STR)
        ELSE
          LENGTH=INDEX(STR,TRM)-1
          IF(LENGTH.LT.0) LENGTH=LEN(STR)
        ENDIF
!*     Find beginning and end of the NTH DEL-limited substring in STR
        END=-1
        DO 1,N=1,NTH
        IF(END.EQ.LENGTH) RETURN
        BEG=END+2
        END=BEG+INDEX(STR(BEG:LENGTH),DEL)-2
        IF(END.EQ.BEG-2) END=LENGTH
!*        PRINT *,'NTH,LENGTH,N,BEG,END=',NTH,LENGTH,N,BEG,END
    1   CONTINUE
        NCH=END-BEG+1
        IF(NCH.GT.0) SST=STR(BEG:END)
      ENDIF
      END
      CHARACTER*(*) FUNCTION SPACES(STR,NSPACE)
!
! CERN PROGLIB# M433    SPACES          .VERSION KERNFOR  4.14  860211
! ORIG.  6/05/86 M.GOOSSENS/DD
!
!-    The function value SPACES returns the character string STR with
!-    leading blanks removed and each occurence of one or more blanks
!-    replaced by NSPACE blanks inside the string STR
!
      CHARACTER*(*) STR
!
      LENSPA = LEN(SPACES)
      SPACES = ' '
      IF (NSPACE.LT.0) NSPACE = 0
      IBLANK = 1
      ISPACE = 1
  100 INONBL = INDEXC(STR(IBLANK:),' ')
      IF (INONBL.EQ.0) THEN
          SPACES(ISPACE:) = STR(IBLANK:)
                                                    GO TO 999
      ENDIF
      INONBL = INONBL + IBLANK - 1
      IBLANK = INDEX(STR(INONBL:),' ')
      IF (IBLANK.EQ.0) THEN
          SPACES(ISPACE:) = STR(INONBL:)
                                                    GO TO 999
      ENDIF
      IBLANK = IBLANK + INONBL - 1
      SPACES(ISPACE:) = STR(INONBL:IBLANK-1)
      ISPACE = ISPACE + IBLANK - INONBL + NSPACE
      IF (ISPACE.LE.LENSPA)                         GO TO 100
  999 END
      FUNCTION INDEXC(STR,SSTR)
!
! CERN PROGLIB# M433    INDEXC          .VERSION KERNFOR  4.14  860211
! ORIG. 26/03/86 M.GOOSSENS/DD
!
!-    Find the leftmost position where substring SSTR does not match
!-    string STR scanning forward
!
      CHARACTER*(*) STR,SSTR
!
      LENS   = LEN(STR)
      LENSS  = LEN(SSTR)
!
      DO 10 I=1,LENS-LENSS+1
          IF (STR(I:I+LENSS-1).NE.SSTR) THEN
              INDEXC = I
                                         GO TO 999
          ENDIF
   10 CONTINUE
      INDEXC = 0
!
  999 END
