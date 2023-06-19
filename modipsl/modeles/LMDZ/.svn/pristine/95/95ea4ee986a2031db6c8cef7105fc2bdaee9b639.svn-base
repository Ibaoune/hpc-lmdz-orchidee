!
! $Id: 1D_interp_cases.h 3537 2019-06-19 08:29:16Z fhourdin $
!
!---------------------------------------------------------------------
! Forcing_LES case: constant dq_dyn
!---------------------------------------------------------------------
      if (forcing_LES) then
        DO l = 1,llm
          d_q_adv(l,1) = dq_dyn(l,1)
        ENDDO
      endif ! forcing_LES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Interpolation forcing in time and onto model levels
!---------------------------------------------------------------------
      if (forcing_GCSSold) then

       call get_uvd(it,timestep,fich_gcssold_ctl,fich_gcssold_dat,              &
     &               ht_gcssold,hq_gcssold,hw_gcssold,                          &
     &               hu_gcssold,hv_gcssold,                                     &
     &               hthturb_gcssold,hqturb_gcssold,Ts_gcssold,                 &
     &               imp_fcg_gcssold,ts_fcg_gcssold,                            &
     &               Tp_fcg_gcssold,Turb_fcg_gcssold)
       if (prt_level.ge.1) then
         print *,' get_uvd -> hqturb_gcssold ',it,hqturb_gcssold
       endif
! large-scale forcing :
!!!      tsurf = ts_gcssold
      do l = 1, llm
!       u(l) = hu_gcssold(l) !  on prescrit le vent
!       v(l) = hv_gcssold(l)    !  on prescrit le vent
!       omega(l) = hw_gcssold(l)
!       rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
!       omega2(l)=-rho(l)*omega(l)
       omega(l) = hw_gcssold(l)
       omega2(l)= omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
       d_t_adv(l) = ht_gcssold(l)
       d_q_adv(l,1) = hq_gcssold(l)
       dt_cooling(l) = 0.0
      enddo

      endif ! forcing_GCSSold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Interpolation Toga forcing 
!---------------------------------------------------------------------
      if (forcing_toga) then

       if (prt_level.ge.1) then
        print*,                                                             &
     & '#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/dt_toga=',     &
     &    day,day1,(day-day1)*86400.,(day-day1)*86400/dt_toga
       endif

! time interpolation:
        CALL interp_toga_time(daytime,day1,annee_ref                        &
     &             ,year_ini_toga,day_ju_ini_toga,nt_toga,dt_toga           &
     &             ,nlev_toga,ts_toga,plev_toga,t_toga,q_toga,u_toga        &
     &             ,v_toga,w_toga,ht_toga,vt_toga,hq_toga,vq_toga           &
     &             ,ts_prof,plev_prof,t_prof,q_prof,u_prof,v_prof,w_prof    &
     &             ,ht_prof,vt_prof,hq_prof,vq_prof)

        if (type_ts_forcing.eq.1) ts_cur = ts_prof ! SST used in read_tsurf1d

! vertical interpolation:
      CALL interp_toga_vertical(play,nlev_toga,plev_prof                    &
     &         ,t_prof,q_prof,u_prof,v_prof,w_prof                          &
     &         ,ht_prof,vt_prof,hq_prof,vq_prof                             &
     &         ,t_mod,q_mod,u_mod,v_mod,w_mod                               &
     &         ,ht_mod,vt_mod,hq_mod,vq_mod,mxcalc)

! large-scale forcing :
      tsurf = ts_prof
      do l = 1, llm
       u(l) = u_mod(l) ! sb: on prescrit le vent
       v(l) = v_mod(l) ! sb: on prescrit le vent
!       omega(l) = w_prof(l)
!       rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
!       omega2(l)=-rho(l)*omega(l)
       omega(l) = w_mod(l)
       omega2(l)= omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
       d_t_adv(l) = alpha*omega(l)/rcpd-(ht_mod(l)+vt_mod(l))
       d_q_adv(l,1) = -(hq_mod(l)+vq_mod(l))
       dt_cooling(l) = 0.0
      enddo

      endif ! forcing_toga
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolation DICE forcing 
!---------------------------------------------------------------------
      if (forcing_dice) then

       if (prt_level.ge.1) then
        print*,'#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/dt_dice=',&
     &    day,day1,(day-day1)*86400.,(day-day1)*86400/dt_dice
       endif

! time interpolation:
      CALL interp_dice_time(daytime,day1,annee_ref                    &
     &             ,year_ini_dice,day_ju_ini_dice,nt_dice,dt_dice     & 
     &             ,nlev_dice,shf_dice,lhf_dice,lwup_dice,swup_dice   &
     &             ,tg_dice,ustar_dice,psurf_dice,ug_dice,vg_dice     &
     &             ,ht_dice,hq_dice,hu_dice,hv_dice,w_dice,omega_dice &
     &             ,shf_prof,lhf_prof,lwup_prof,swup_prof,tg_prof     &
     &             ,ustar_prof,psurf_prof,ug_profd,vg_profd           &
     &             ,ht_profd,hq_profd,hu_profd,hv_profd,w_profd       &
     &             ,omega_profd)
!     do l = 1, llm
!     print *,'llm l omega_profd',llm,l,omega_profd(l)
!     enddo

        if (type_ts_forcing.eq.1) ts_cur = tg_prof ! SST used in read_tsurf1d

! vertical interpolation:
      CALL interp_dice_vertical(play,nlev_dice,nt_dice,plev_dice        &
     &         ,t_dice,qv_dice,u_dice,v_dice,o3_dice                   &
     &         ,ht_profd,hq_profd,hu_profd,hv_profd,w_profd,omega_profd &
     &         ,t_mod,qv_mod,u_mod,v_mod,o3_mod                        &
     &         ,ht_mod,hq_mod,hu_mod,hv_mod,w_mod,omega_mod,mxcalc)
!     do l = 1, llm
!      print *,'llm l omega_mod',llm,l,omega_mod(l)
!     enddo

! Les forcages DICE sont donnes /jour et non /seconde !
      ht_mod(:)=ht_mod(:)/86400.
      hq_mod(:)=hq_mod(:)/86400.
      hu_mod(:)=hu_mod(:)/86400.
      hv_mod(:)=hv_mod(:)/86400.

!calcul de l'advection verticale a partir du omega (repris cas TWPICE, MPL 05082013)
!Calcul des gradients verticaux
!initialisation
      d_t_z(:)=0.
      d_q_z(:)=0.
      d_u_z(:)=0.
      d_v_z(:)=0.
      DO l=2,llm-1
       d_t_z(l)=(temp(l+1)-temp(l-1))/(play(l+1)-play(l-1))
       d_q_z(l)=(q(l+1,1)-q(l-1,1)) /(play(l+1)-play(l-1))
       d_u_z(l)=(u(l+1)-u(l-1))/(play(l+1)-play(l-1))
       d_v_z(l)=(v(l+1)-v(l-1))/(play(l+1)-play(l-1))
      ENDDO
      d_t_z(1)=d_t_z(2)
      d_q_z(1)=d_q_z(2)
!     d_u_z(1)=u(2)/(play(2)-psurf)/5.
!     d_v_z(1)=v(2)/(play(2)-psurf)/5.
      d_u_z(1)=0.
      d_v_z(1)=0.
      d_t_z(llm)=d_t_z(llm-1)
      d_q_z(llm)=d_q_z(llm-1)
      d_u_z(llm)=d_u_z(llm-1)
      d_v_z(llm)=d_v_z(llm-1)

!Calcul de l advection verticale: 
! utiliser omega (Pa/s) et non w (m/s) !! MP 20131108
      d_t_dyn_z(:)=omega_mod(:)*d_t_z(:)
      d_q_dyn_z(:)=omega_mod(:)*d_q_z(:)
      d_u_dyn_z(:)=omega_mod(:)*d_u_z(:)
      d_v_dyn_z(:)=omega_mod(:)*d_v_z(:)

! large-scale forcing :
!     tsurf = tg_prof    MPL 20130925 commente
      psurf = psurf_prof
! For this case, fluxes are imposed
      fsens=-1*shf_prof
      flat=-1*lhf_prof
      ust=ustar_prof
      tg=tg_prof
      print *,'ust= ',ust
      do l = 1, llm
       ug(l)= ug_profd
       vg(l)= vg_profd
!       omega(l) = w_prof(l)
!      rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
!       omega2(l)=-rho(l)*omega(l)
!      omega(l) = w_mod(l)*(-rg*rho(l))
       omega(l) = omega_mod(l)
       omega2(l)= omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
       d_t_adv(l) = alpha*omega(l)/rcpd+ht_mod(l)-d_t_dyn_z(l)
       d_q_adv(l,1) = hq_mod(l)-d_q_dyn_z(l)
       d_u_adv(l) = hu_mod(l)-d_u_dyn_z(l)
       d_v_adv(l) = hv_mod(l)-d_v_dyn_z(l)
       dt_cooling(l) = 0.0
      enddo

      endif ! forcing_dice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolation gabls4 forcing 
!---------------------------------------------------------------------
      if (forcing_gabls4 ) then

       if (prt_level.ge.1) then
        print*,'#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/dt_gabls4=',&
     &    day,day1,(day-day1)*86400.,(day-day1)*86400/dt_gabls4
       endif

! time interpolation:
      CALL interp_gabls4_time(daytime,day1,annee_ref                                     &
     &             ,year_ini_gabls4,day_ju_ini_gabls4,nt_gabls4,dt_gabls4,nlev_gabls4  & 
     &             ,ug_gabls4,vg_gabls4,ht_gabls4,hq_gabls4,tg_gabls4                            &
     &             ,ug_profg,vg_profg,ht_profg,hq_profg,tg_profg)

        if (type_ts_forcing.eq.1) ts_cur = tg_prof ! SST used in read_tsurf1d

! vertical interpolation:
! on re-utilise le programme interp_dice_vertical: les transformations sur 
! plev_gabls4,th_gabls4,qv_gabls4,u_gabls4,v_gabls4 ne sont pas prises en compte.
! seules celles sur ht_profg,hq_profg,ug_profg,vg_profg sont prises en compte.

      CALL interp_dice_vertical(play,nlev_gabls4,nt_gabls4,plev_gabls4         &
!    &         ,t_gabls4,qv_gabls4,u_gabls4,v_gabls4,poub            &
     &         ,poub,poub,poub,poub,poub                             &
     &         ,ht_profg,hq_profg,ug_profg,vg_profg,poub,poub        &
     &         ,t_mod,qv_mod,u_mod,v_mod,o3_mod                      &
     &         ,ht_mod,hq_mod,ug_mod,vg_mod,w_mod,omega_mod,mxcalc)

      do l = 1, llm
       ug(l)= ug_mod(l)
       vg(l)= vg_mod(l)
       d_t_adv(l)=ht_mod(l)
       d_q_adv(l,1)=hq_mod(l)
      enddo

      endif ! forcing_gabls4
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Interpolation forcing TWPice
!---------------------------------------------------------------------
      if (forcing_twpice) then

        print*,                                                             &
     & '#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/dt_twpi=',     &
     &    daytime,day1,(daytime-day1)*86400.,                               &
     &    (daytime-day1)*86400/dt_twpi

! time interpolation:
        CALL interp_toga_time(daytime,day1,annee_ref                        &
     &       ,year_ini_twpi,day_ju_ini_twpi,nt_twpi,dt_twpi,nlev_twpi       &
     &       ,ts_twpi,plev_twpi,t_twpi,q_twpi,u_twpi,v_twpi,w_twpi          &
     &       ,ht_twpi,vt_twpi,hq_twpi,vq_twpi                               &
     &       ,ts_proftwp,plev_proftwp,t_proftwp,q_proftwp,u_proftwp         &
     &       ,v_proftwp,w_proftwp                                           &
     &       ,ht_proftwp,vt_proftwp,hq_proftwp,vq_proftwp)

! vertical interpolation:
      CALL interp_toga_vertical(play,nlev_twpi,plev_proftwp                 &
     &         ,t_proftwp,q_proftwp,u_proftwp,v_proftwp,w_proftwp           &
     &         ,ht_proftwp,vt_proftwp,hq_proftwp,vq_proftwp                 &
     &         ,t_mod,q_mod,u_mod,v_mod,w_mod                               &
     &         ,ht_mod,vt_mod,hq_mod,vq_mod,mxcalc)


!calcul de l'advection verticale a partir du omega
!Calcul des gradients verticaux
!initialisation
      d_t_z(:)=0.
      d_q_z(:)=0.
      d_t_dyn_z(:)=0.
      d_q_dyn_z(:)=0.
      DO l=2,llm-1
       d_t_z(l)=(temp(l+1)-temp(l-1))/(play(l+1)-play(l-1))
       d_q_z(l)=(q(l+1,1)-q(l-1,1))/(play(l+1)-play(l-1))
      ENDDO
      d_t_z(1)=d_t_z(2)
      d_q_z(1)=d_q_z(2)
      d_t_z(llm)=d_t_z(llm-1)
      d_q_z(llm)=d_q_z(llm-1)

!Calcul de l advection verticale
      d_t_dyn_z(:)=w_mod(:)*d_t_z(:)
      d_q_dyn_z(:)=w_mod(:)*d_q_z(:)

!wind nudging above 500m with a 2h time scale
        do l=1,llm
        if (nudge_wind) then
!           if (phi(l).gt.5000.) then
        if (phi(l).gt.0.) then
        u(l)=u(l)+timestep*(u_mod(l)-u(l))/(2.*3600.)
        v(l)=v(l)+timestep*(v_mod(l)-v(l))/(2.*3600.)
           endif    
        else
        u(l) = u_mod(l) 
        v(l) = v_mod(l)
        endif
        enddo

!CR:nudging of q and theta with a 6h time scale above 15km
        if (nudge_thermo) then
        do l=1,llm
           zz(l)=phi(l)/9.8
           if ((zz(l).le.16000.).and.(zz(l).gt.15000.)) then
             zfact=(zz(l)-15000.)/1000.
        q(l,1)=q(l,1)+timestep*(q_mod(l)-q(l,1))/(6.*3600.)*zfact
        temp(l)=temp(l)+timestep*(t_mod(l)-temp(l))/(6.*3600.)*zfact
           else if (zz(l).gt.16000.) then 
        q(l,1)=q(l,1)+timestep*(q_mod(l)-q(l,1))/(6.*3600.)
        temp(l)=temp(l)+timestep*(t_mod(l)-temp(l))/(6.*3600.)
           endif
        enddo    
        endif

      do l = 1, llm
       omega(l) = w_mod(l)
       omega2(l)= omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq
       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!calcul de l'advection totale
        if (cptadvw) then
        d_t_adv(l) = alpha*omega(l)/rcpd+ht_mod(l)-d_t_dyn_z(l)
!        print*,'temp vert adv',l,ht_mod(l),vt_mod(l),-d_t_dyn_z(l)
        d_q_adv(l,1) = hq_mod(l)-d_q_dyn_z(l)
!        print*,'q vert adv',l,hq_mod(l),vq_mod(l),-d_q_dyn_z(l)
        else
        d_t_adv(l) = alpha*omega(l)/rcpd+(ht_mod(l)+vt_mod(l))
        d_q_adv(l,1) = (hq_mod(l)+vq_mod(l))
        endif
       dt_cooling(l) = 0.0
      enddo

      endif ! forcing_twpice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Interpolation forcing AMMA
!---------------------------------------------------------------------

       if (forcing_amma) then

        print*,                                                             &
     & '#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/dt_amma=',     &
     &    daytime,day1,(daytime-day1)*86400.,                               &
     &    (daytime-day1)*86400/dt_amma

! time interpolation using TOGA interpolation routine
        CALL interp_amma_time(daytime,day1,annee_ref                        &
     &       ,year_ini_amma,day_ju_ini_amma,nt_amma,dt_amma,nlev_amma       &
     &       ,vitw_amma,ht_amma,hq_amma,lat_amma,sens_amma                  &
     &       ,vitw_profamma,ht_profamma,hq_profamma,lat_profamma            &
     &       ,sens_profamma)

      print*,'apres interpolation temporelle AMMA'

      do k=1,nlev_amma
         th_profamma(k)=0.
         q_profamma(k)=0.
         u_profamma(k)=0.
         v_profamma(k)=0.
         vt_profamma(k)=0.
         vq_profamma(k)=0.
       enddo
! vertical interpolation using TOGA interpolation routine:
!      write(*,*)'avant interp vert', t_proftwp
      CALL interp_toga_vertical(play,nlev_amma,plev_amma                      &
     &         ,th_profamma,q_profamma,u_profamma,v_profamma                 &
     &         ,vitw_profamma                                               &
     &         ,ht_profamma,vt_profamma,hq_profamma,vq_profamma             &
     &         ,t_mod,q_mod,u_mod,v_mod,w_mod                               &
     &         ,ht_mod,vt_mod,hq_mod,vq_mod,mxcalc)
       write(*,*) 'Profil initial forcing AMMA interpole'


!calcul de l'advection verticale a partir du omega
!Calcul des gradients verticaux
!initialisation
      do l=1,llm
      d_t_z(l)=0.
      d_q_z(l)=0.
      enddo

      DO l=2,llm-1
       d_t_z(l)=(temp(l+1)-temp(l-1))/(play(l+1)-play(l-1))
       d_q_z(l)=(q(l+1,1)-q(l-1,1))/(play(l+1)-play(l-1))
      ENDDO
      d_t_z(1)=d_t_z(2)
      d_q_z(1)=d_q_z(2)
      d_t_z(llm)=d_t_z(llm-1)
      d_q_z(llm)=d_q_z(llm-1)


      do l = 1, llm
       rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
       omega(l) = w_mod(l)*(-rg*rho(l))
       omega2(l)= omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq
       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!calcul de l'advection totale
!        d_t_adv(l) = alpha*omega(l)/rcpd+ht_mod(l)-omega(l)*d_t_z(l)
!attention: on impose dth
        d_t_adv(l) = alpha*omega(l)/rcpd+                                  &
     &         ht_mod(l)*(play(l)/pzero)**rkappa-omega(l)*d_t_z(l)
!        d_t_adv(l) = 0.
!        print*,'temp vert adv',l,ht_mod(l),vt_mod(l),-d_t_dyn_z(l)
        d_q_adv(l,1) = hq_mod(l)-omega(l)*d_q_z(l)
!        d_q_adv(l,1) = 0.
!        print*,'q vert adv',l,hq_mod(l),vq_mod(l),-d_q_dyn_z(l)
   
       dt_cooling(l) = 0.0
      enddo


!     ok_flux_surf=.false.
      fsens=-1.*sens_profamma
      flat=-1.*lat_profamma

      endif ! forcing_amma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Interpolation forcing Rico
!---------------------------------------------------------------------
      if (forcing_rico) then
!      call lstendH(llm,omega,dt_dyn,dq_dyn,du_dyn, dv_dyn,q,temp,u,v,play)
       call lstendH(llm,nqtot,omega,dt_dyn,dq_dyn,q,temp,u,v,play)

        do l=1,llm
       d_t_adv(l) =  (dth_rico(l) +  dt_dyn(l))
       d_q_adv(l,1) = (dqh_rico(l) +  dq_dyn(l,1))
       d_q_adv(l,2) = 0.
        enddo
      endif  ! forcing_rico
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Interpolation forcing Arm_cu
!---------------------------------------------------------------------
      if (forcing_armcu) then

        print*,                                                             &
     & '#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/dt_armcu=',    &
     &    day,day1,(day-day1)*86400.,(day-day1)*86400/dt_armcu

! time interpolation:
! ATTENTION, cet appel ne convient pas pour TOGA !!
! revoir 1DUTILS.h et les arguments
      CALL interp_armcu_time(daytime,day1,annee_ref                         &
     &            ,year_ini_armcu,day_ju_ini_armcu,nt_armcu,dt_armcu        &
     &            ,nlev_armcu,sens_armcu,flat_armcu,adv_theta_armcu          &
     &            ,rad_theta_armcu,adv_qt_armcu,sens_prof,flat_prof         &
     &            ,adv_theta_prof,rad_theta_prof,adv_qt_prof)

! vertical interpolation:
! No vertical interpolation if nlev imposed to 19 or 40

! For this case, fluxes are imposed
       fsens=-1*sens_prof
       flat=-1*flat_prof

! Advective forcings are given in K or g/kg ... BY HOUR
      do l = 1, llm
       ug(l)= u_mod(l)
       vg(l)= v_mod(l)
       IF((phi(l)/RG).LT.1000) THEN
         d_t_adv(l) = (adv_theta_prof + rad_theta_prof)/3600.
         d_q_adv(l,1) = adv_qt_prof/1000./3600.
         d_q_adv(l,2) = 0.0
!        print *,'INF1000: phi dth dq1 dq2',
!    :  phi(l)/RG,d_t_adv(l),d_q_adv(l,1),d_q_adv(l,2)
       ELSEIF ((phi(l)/RG).GE.1000.AND.(phi(l)/RG).lt.3000) THEN
         fact=((phi(l)/RG)-1000.)/2000.
         fact=1-fact
         d_t_adv(l) = (adv_theta_prof + rad_theta_prof)*fact/3600.
         d_q_adv(l,1) = adv_qt_prof*fact/1000./3600.
         d_q_adv(l,2) = 0.0
!        print *,'SUP1000: phi fact dth dq1 dq2',
!    :  phi(l)/RG,fact,d_t_adv(l),d_q_adv(l,1),d_q_adv(l,2)
       ELSE
         d_t_adv(l) = 0.0
         d_q_adv(l,1) = 0.0
         d_q_adv(l,2) = 0.0
!        print *,'SUP3000: phi dth dq1 dq2',
!    :  phi(l)/RG,d_t_adv(l),d_q_adv(l,1),d_q_adv(l,2)
       ENDIF
      dt_cooling(l) = 0.0  
!     print *,'Interp armcu: phi dth dq1 dq2',
!    :  l,phi(l),d_t_adv(l),d_q_adv(l,1),d_q_adv(l,2)
      enddo
      endif ! forcing_armcu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Interpolation forcing in time and onto model levels
!---------------------------------------------------------------------
      if (forcing_sandu) then

        print*,                                                             &
     & '#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/dt_sandu=',    &
     &    day,day1,(day-day1)*86400.,(day-day1)*86400/dt_sandu 

! time interpolation:
! ATTENTION, cet appel ne convient pas pour TOGA !!
! revoir 1DUTILS.h et les arguments
      CALL interp_sandu_time(daytime,day1,annee_ref                         &
     &             ,year_ini_sandu,day_ju_ini_sandu,nt_sandu,dt_sandu       &
     &             ,nlev_sandu                                              &
     &             ,ts_sandu,ts_prof)

        if (type_ts_forcing.eq.1) ts_cur = ts_prof ! SST used in read_tsurf1d

! vertical interpolation:
      CALL interp_sandu_vertical(play,nlev_sandu,plev_profs                 &
     &         ,t_profs,thl_profs,q_profs,u_profs,v_profs,w_profs           &
     &         ,omega_profs,o3mmr_profs                                     &
     &         ,t_mod,thl_mod,q_mod,u_mod,v_mod,w_mod                       &
     &         ,omega_mod,o3mmr_mod,mxcalc)
!calcul de l'advection verticale
!Calcul des gradients verticaux
!initialisation
      d_t_z(:)=0.
      d_q_z(:)=0.
      d_t_dyn_z(:)=0.
      d_q_dyn_z(:)=0.
! schema centre
!     DO l=2,llm-1
!      d_t_z(l)=(temp(l+1)-temp(l-1))
!    &          /(play(l+1)-play(l-1))
!      d_q_z(l)=(q(l+1,1)-q(l-1,1))
!    &          /(play(l+1)-play(l-1))
! schema amont
      DO l=2,llm-1
       d_t_z(l)=(temp(l+1)-temp(l))/(play(l+1)-play(l))
       d_q_z(l)=(q(l+1,1)-q(l,1))/(play(l+1)-play(l))
!     print *,'l temp2 temp0 play2 play0 omega_mod',
!    & temp(l+1),temp(l-1),play(l+1),play(l-1),omega_mod(l)
      ENDDO
      d_t_z(1)=d_t_z(2)
      d_q_z(1)=d_q_z(2)
      d_t_z(llm)=d_t_z(llm-1)
      d_q_z(llm)=d_q_z(llm-1)

!  calcul de l advection verticale
! Confusion w (m/s) et omega (Pa/s) !!
      d_t_dyn_z(:)=omega_mod(:)*d_t_z(:)
      d_q_dyn_z(:)=omega_mod(:)*d_q_z(:)
!     do l=1,llm
!      print *,'d_t_dyn omega_mod d_t_z d_q_dyn d_q_z',
!    :l,d_t_dyn_z(l),omega_mod(l),d_t_z(l),d_q_dyn_z(l),d_q_z(l)
!     enddo


! large-scale forcing : pour le cas Sandu ces forcages sont la SST
! et une divergence constante -> profil de omega
      tsurf = ts_prof
      write(*,*) 'SST suivante: ',tsurf
      do l = 1, llm
       omega(l) = omega_mod(l)
       omega2(l)= omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!
!      d_t_adv(l) = 0.0
!      d_q_adv(l,1) = 0.0
!CR:test advection=0
!calcul de l'advection verticale
        d_t_adv(l) = alpha*omega(l)/rcpd-d_t_dyn_z(l)
!        print*,'temp adv',l,-d_t_dyn_z(l)
        d_q_adv(l,1) = -d_q_dyn_z(l)
!        print*,'q adv',l,-d_q_dyn_z(l)
       dt_cooling(l) = 0.0
      enddo
      endif ! forcing_sandu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Interpolation forcing in time and onto model levels
!---------------------------------------------------------------------
      if (forcing_astex) then

        print*,                                                             &
     & '#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/dt_astex=',    &
     &    day,day1,(day-day1)*86400.,(day-day1)*86400/dt_astex

! time interpolation:
! ATTENTION, cet appel ne convient pas pour TOGA !!
! revoir 1DUTILS.h et les arguments
      CALL interp_astex_time(daytime,day1,annee_ref                         &
     &             ,year_ini_astex,day_ju_ini_astex,nt_astex,dt_astex       &
     &             ,nlev_astex,div_astex,ts_astex,ug_astex,vg_astex         &
     &             ,ufa_astex,vfa_astex,div_prof,ts_prof,ug_prof,vg_prof    &
     &             ,ufa_prof,vfa_prof)

        if (type_ts_forcing.eq.1) ts_cur = ts_prof ! SST used in read_tsurf1d

! vertical interpolation:
      CALL interp_astex_vertical(play,nlev_astex,plev_profa                 &
     &         ,t_profa,thl_profa,qv_profa,ql_profa,qt_profa                &
     &         ,u_profa,v_profa,w_profa,tke_profa,o3mmr_profa               &
     &         ,t_mod,thl_mod,qv_mod,ql_mod,qt_mod,u_mod,v_mod,w_mod        &
     &         ,tke_mod,o3mmr_mod,mxcalc)
!calcul de l'advection verticale
!Calcul des gradients verticaux
!initialisation
      d_t_z(:)=0.
      d_q_z(:)=0.
      d_t_dyn_z(:)=0.
      d_q_dyn_z(:)=0.
! schema centre
!     DO l=2,llm-1
!      d_t_z(l)=(temp(l+1)-temp(l-1))
!    &          /(play(l+1)-play(l-1))
!      d_q_z(l)=(q(l+1,1)-q(l-1,1))
!    &          /(play(l+1)-play(l-1))
! schema amont
      DO l=2,llm-1
       d_t_z(l)=(temp(l+1)-temp(l))/(play(l+1)-play(l))
       d_q_z(l)=(q(l+1,1)-q(l,1))/(play(l+1)-play(l))
!     print *,'l temp2 temp0 play2 play0 omega_mod',
!    & temp(l+1),temp(l-1),play(l+1),play(l-1),omega_mod(l)
      ENDDO
      d_t_z(1)=d_t_z(2)
      d_q_z(1)=d_q_z(2)
      d_t_z(llm)=d_t_z(llm-1)
      d_q_z(llm)=d_q_z(llm-1)

!  calcul de l advection verticale
! Confusion w (m/s) et omega (Pa/s) !!
      d_t_dyn_z(:)=w_mod(:)*d_t_z(:)
      d_q_dyn_z(:)=w_mod(:)*d_q_z(:)
!     do l=1,llm
!      print *,'d_t_dyn omega_mod d_t_z d_q_dyn d_q_z',
!    :l,d_t_dyn_z(l),omega_mod(l),d_t_z(l),d_q_dyn_z(l),d_q_z(l)
!     enddo


! large-scale forcing : pour le cas Astex ces forcages sont la SST
! la divergence,ug,vg,ufa,vfa
      tsurf = ts_prof
      write(*,*) 'SST suivante: ',tsurf
      do l = 1, llm
       omega(l) = w_mod(l)
       omega2(l)= omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!
!      d_t_adv(l) = 0.0
!      d_q_adv(l,1) = 0.0
!CR:test advection=0
!calcul de l'advection verticale
        d_t_adv(l) = alpha*omega(l)/rcpd-d_t_dyn_z(l)
!        print*,'temp adv',l,-d_t_dyn_z(l)
        d_q_adv(l,1) = -d_q_dyn_z(l)
!        print*,'q adv',l,-d_q_dyn_z(l)
       dt_cooling(l) = 0.0
      enddo
      endif ! forcing_astex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Interpolation forcing standard case
!---------------------------------------------------------------------
      if (forcing_case) then

         print*,'FORCING CASE forcing_case'

        print*,                                                             &
     & '#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/pdt_cas=',     &
     &    daytime,day1,(daytime-day1)*86400.,                               &
     &    (daytime-day1)*86400/pdt_cas

! time interpolation:
        CALL interp_case_time(daytime,day1,annee_ref                                        &
!    &       ,year_ini_cas,day_ju_ini_cas,nt_cas,pdt_cas,nlev_cas                           &
     &       ,nt_cas,nlev_cas                                                               &
     &       ,ts_cas,plev_cas,t_cas,q_cas,u_cas,v_cas,ug_cas,vg_cas                         &
     &       ,vitw_cas,du_cas,hu_cas,vu_cas                                                 &
     &       ,dv_cas,hv_cas,vv_cas,dt_cas,ht_cas,vt_cas,dtrad_cas                           &
     &       ,dq_cas,hq_cas,vq_cas,lat_cas,sens_cas,ustar_cas                               &
     &       ,uw_cas,vw_cas,q1_cas,q2_cas                                                   &
     &       ,ts_prof_cas,plev_prof_cas,t_prof_cas,q_prof_cas,u_prof_cas,v_prof_cas         &
     &       ,ug_prof_cas,vg_prof_cas,vitw_prof_cas,du_prof_cas,hu_prof_cas,vu_prof_cas     &
     &       ,dv_prof_cas,hv_prof_cas,vv_prof_cas,dt_prof_cas,ht_prof_cas,vt_prof_cas       &
     &       ,dtrad_prof_cas,dq_prof_cas,hq_prof_cas,vq_prof_cas,lat_prof_cas               &
     &       ,sens_prof_cas,ustar_prof_cas,uw_prof_cas,vw_prof_cas,q1_prof_cas,q2_prof_cas)

             ts_cur = ts_prof_cas 
             psurf=plev_prof_cas(1)

! vertical interpolation:
      CALL interp_case_vertical(play,nlev_cas,plev_prof_cas            &
     &         ,t_prof_cas,q_prof_cas,u_prof_cas,v_prof_cas,ug_prof_cas,vg_prof_cas,vitw_prof_cas                         &
     &         ,du_prof_cas,hu_prof_cas,vu_prof_cas,dv_prof_cas,hv_prof_cas,vv_prof_cas           &
     &         ,dt_prof_cas,ht_prof_cas,vt_prof_cas,dtrad_prof_cas,dq_prof_cas,hq_prof_cas,vq_prof_cas           &
     &         ,t_mod_cas,q_mod_cas,u_mod_cas,v_mod_cas,ug_mod_cas,vg_mod_cas,w_mod_cas                              &
     &         ,du_mod_cas,hu_mod_cas,vu_mod_cas,dv_mod_cas,hv_mod_cas,vv_mod_cas               &
     &         ,dt_mod_cas,ht_mod_cas,vt_mod_cas,dtrad_mod_cas,dq_mod_cas,hq_mod_cas,vq_mod_cas,mxcalc)


!calcul de l'advection verticale a partir du omega
!Calcul des gradients verticaux
!initialisation
      d_t_z(:)=0.
      d_q_z(:)=0.
      d_u_z(:)=0.
      d_v_z(:)=0.
      d_t_dyn_z(:)=0.
      d_q_dyn_z(:)=0.
      d_u_dyn_z(:)=0.
      d_v_dyn_z(:)=0.
      DO l=2,llm-1
       d_t_z(l)=(temp(l+1)-temp(l-1))/(play(l+1)-play(l-1))
       d_q_z(l)=(q(l+1,1)-q(l-1,1))/(play(l+1)-play(l-1))
       d_u_z(l)=(u(l+1)-u(l-1))/(play(l+1)-play(l-1))
       d_v_z(l)=(v(l+1)-v(l-1))/(play(l+1)-play(l-1))
      ENDDO
      d_t_z(1)=d_t_z(2)
      d_q_z(1)=d_q_z(2)
      d_u_z(1)=d_u_z(2)
      d_v_z(1)=d_v_z(2)
      d_t_z(llm)=d_t_z(llm-1)
      d_q_z(llm)=d_q_z(llm-1)
      d_u_z(llm)=d_u_z(llm-1)
      d_v_z(llm)=d_v_z(llm-1)

!Calcul de l advection verticale

      d_t_dyn_z(:)=w_mod_cas(:)*d_t_z(:)

      d_q_dyn_z(:)=w_mod_cas(:)*d_q_z(:)
      d_u_dyn_z(:)=w_mod_cas(:)*d_u_z(:)
      d_v_dyn_z(:)=w_mod_cas(:)*d_v_z(:)

!wind nudging 
      if (nudge_u.gt.0.) then
        do l=1,llm
           u(l)=u(l)+timestep*(u_mod_cas(l)-u(l))/(nudge_u)
        enddo
      else
        do l=1,llm
        u(l) = u_mod_cas(l) 
        enddo
      endif

      if (nudge_v.gt.0.) then
        do l=1,llm
           v(l)=v(l)+timestep*(v_mod_cas(l)-v(l))/(nudge_v)
        enddo
      else
        do l=1,llm
        v(l) = v_mod_cas(l) 
        enddo
      endif

      if (nudge_w.gt.0.) then
        do l=1,llm
           w(l)=w(l)+timestep*(w_mod_cas(l)-w(l))/(nudge_w)
        enddo
      else
        do l=1,llm
        w(l) = w_mod_cas(l) 
        enddo
      endif

!nudging of q and temp 
      if (nudge_t.gt.0.) then
        do l=1,llm
           temp(l)=temp(l)+timestep*(t_mod_cas(l)-temp(l))/(nudge_t)
        enddo
      endif
      if (nudge_q.gt.0.) then
        do l=1,llm
           q(l,1)=q(l,1)+timestep*(q_mod_cas(l)-q(l,1))/(nudge_q)
        enddo
      endif

      do l = 1, llm
       omega(l) = w_mod_cas(l)  ! juste car w_mod_cas en Pa/s (MPL 20170310)
       omega2(l)= omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq
       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)

!calcul advection
        if ((tend_u.eq.1).and.(tend_w.eq.0)) then
           d_u_adv(l)=du_mod_cas(l)
        else if ((tend_u.eq.1).and.(tend_w.eq.1)) then
           d_u_adv(l)=hu_mod_cas(l)-d_u_dyn_z(l)
        endif

        if ((tend_v.eq.1).and.(tend_w.eq.0)) then
           d_v_adv(l)=dv_mod_cas(l)
        else if ((tend_v.eq.1).and.(tend_w.eq.1)) then
           d_v_adv(l)=hv_mod_cas(l)-d_v_dyn_z(l)
        endif

        if ((tend_t.eq.1).and.(tend_w.eq.0)) then
!           d_t_adv(l)=alpha*omega(l)/rcpd+dt_mod_cas(l)
           d_t_adv(l)=alpha*omega(l)/rcpd-dt_mod_cas(l)
        else if ((tend_t.eq.1).and.(tend_w.eq.1)) then
!           d_t_adv(l)=alpha*omega(l)/rcpd+ht_mod_cas(l)-d_t_dyn_z(l)
           d_t_adv(l)=alpha*omega(l)/rcpd-ht_mod_cas(l)-d_t_dyn_z(l)
        endif

        if ((tend_q.eq.1).and.(tend_w.eq.0)) then
!           d_q_adv(l,1)=dq_mod_cas(l)
           d_q_adv(l,1)=-1*dq_mod_cas(l)
        else if ((tend_q.eq.1).and.(tend_w.eq.1)) then
!           d_q_adv(l,1)=hq_mod_cas(l)-d_q_dyn_z(l)
           d_q_adv(l,1)=-1*hq_mod_cas(l)-d_q_dyn_z(l)
        endif
         
        if (tend_rayo.eq.1) then
           dt_cooling(l) = dtrad_mod_cas(l)
!          print *,'dt_cooling=',dt_cooling(l)
        else
           dt_cooling(l) = 0.0
        endif
      enddo

! Faut-il multiplier par -1 ? (MPL 20160713)
      IF(ok_flux_surf) THEN
       fsens=sens_prof_cas
       flat=lat_prof_cas
      ENDIF
!
      IF (ok_prescr_ust) THEN
       ust=ustar_prof_cas
       print *,'ust=',ust
      ENDIF
      endif ! forcing_case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Interpolation forcing standard case
!---------------------------------------------------------------------
      if (forcing_case2 .OR. forcing_SCM) then

         print*,'FORCING CASE forcing_case2'
        print*,                                                             &
     & '#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/pdt_cas=',     &
     &    daytime,day1,(daytime-day1)*86400.,                               &
     &    (daytime-day1)*86400/pdt_cas

! time interpolation:
        CALL interp2_case_time(daytime,day1,annee_ref                                       &
!    &       ,year_ini_cas,day_ju_ini_cas,nt_cas,pdt_cas,nlev_cas                           &
     &       ,nt_cas,nlev_cas                                                               &
     &       ,ts_cas,ps_cas,plev_cas,t_cas,th_cas,thv_cas,thl_cas,qv_cas,ql_cas,qi_cas      &
     &       ,u_cas,v_cas,ug_cas,vg_cas,vitw_cas,omega_cas,du_cas,hu_cas,vu_cas             &
     &       ,dv_cas,hv_cas,vv_cas,dt_cas,ht_cas,vt_cas,dtrad_cas                           &
     &       ,dq_cas,hq_cas,vq_cas,dth_cas,hth_cas,vth_cas,lat_cas,sens_cas,ustar_cas       &
     &       ,uw_cas,vw_cas,q1_cas,q2_cas,tke_cas                                           &
!
     &       ,ts_prof_cas,plev_prof_cas,t_prof_cas,theta_prof_cas,thv_prof_cas  &
     &       ,thl_prof_cas,qv_prof_cas,ql_prof_cas,qi_prof_cas                              &
     &       ,u_prof_cas,v_prof_cas,ug_prof_cas,vg_prof_cas,vitw_prof_cas,omega_prof_cas    &
     &       ,du_prof_cas,hu_prof_cas,vu_prof_cas                                           &
     &       ,dv_prof_cas,hv_prof_cas,vv_prof_cas,dt_prof_cas,ht_prof_cas,vt_prof_cas       &
     &       ,dtrad_prof_cas,dq_prof_cas,hq_prof_cas,vq_prof_cas                            &
     &       ,dth_prof_cas,hth_prof_cas,vth_prof_cas,lat_prof_cas                           &
     &       ,sens_prof_cas,ustar_prof_cas,uw_prof_cas,vw_prof_cas,q1_prof_cas,q2_prof_cas,tke_prof_cas)

             ts_cur = ts_prof_cas 
!            psurf=plev_prof_cas(1)
             psurf=ps_prof_cas

! vertical interpolation:
      CALL interp2_case_vertical(play,nlev_cas,plev_prof_cas                                              &
     &         ,t_prof_cas,theta_prof_cas,thv_prof_cas,thl_prof_cas                                          &
     &         ,qv_prof_cas,ql_prof_cas,qi_prof_cas,u_prof_cas,v_prof_cas                                 &
     &         ,ug_prof_cas,vg_prof_cas,vitw_prof_cas,omega_prof_cas                                      &
     &         ,du_prof_cas,hu_prof_cas,vu_prof_cas,dv_prof_cas,hv_prof_cas,vv_prof_cas                   &
     &         ,dt_prof_cas,ht_prof_cas,vt_prof_cas,dtrad_prof_cas,dq_prof_cas,hq_prof_cas,vq_prof_cas    &
     &         ,dth_prof_cas,hth_prof_cas,vth_prof_cas                                                    &
!
     &         ,t_mod_cas,theta_mod_cas,thv_mod_cas,thl_mod_cas,qv_mod_cas,ql_mod_cas,qi_mod_cas          &
     &         ,u_mod_cas,v_mod_cas,ug_mod_cas,vg_mod_cas,w_mod_cas,omega_mod_cas                         &
     &         ,du_mod_cas,hu_mod_cas,vu_mod_cas,dv_mod_cas,hv_mod_cas,vv_mod_cas                         &
     &         ,dt_mod_cas,ht_mod_cas,vt_mod_cas,dtrad_mod_cas,dq_mod_cas,hq_mod_cas,vq_mod_cas           &
     &         ,dth_mod_cas,hth_mod_cas,vth_mod_cas,mxcalc)


      DO l=1,llm
      teta(l)=temp(l)*(100000./play(l))**(rd/rcpd)
      ENDDO
!calcul de l'advection verticale a partir du omega
!Calcul des gradients verticaux
!initialisation
      d_t_z(:)=0.
      d_th_z(:)=0.
      d_q_z(:)=0.
      d_u_z(:)=0.
      d_v_z(:)=0.
      d_t_dyn_z(:)=0.
      d_th_dyn_z(:)=0.
      d_q_dyn_z(:)=0.
      d_u_dyn_z(:)=0.
      d_v_dyn_z(:)=0.
      DO l=2,llm-1
       d_t_z(l)=(temp(l+1)-temp(l-1))/(play(l+1)-play(l-1))
       d_th_z(l)=(teta(l+1)-teta(l-1))/(play(l+1)-play(l-1))
       d_q_z(l)=(q(l+1,1)-q(l-1,1))/(play(l+1)-play(l-1))
       d_u_z(l)=(u(l+1)-u(l-1))/(play(l+1)-play(l-1))
       d_v_z(l)=(v(l+1)-v(l-1))/(play(l+1)-play(l-1))
      ENDDO
      d_t_z(1)=d_t_z(2)
      d_th_z(1)=d_th_z(2)
      d_q_z(1)=d_q_z(2)
      d_u_z(1)=d_u_z(2)
      d_v_z(1)=d_v_z(2)
      d_t_z(llm)=d_t_z(llm-1)
      d_th_z(llm)=d_th_z(llm-1)
      d_q_z(llm)=d_q_z(llm-1)
      d_u_z(llm)=d_u_z(llm-1)
      d_v_z(llm)=d_v_z(llm-1)

!Calcul de l advection verticale
! Modif w_mod_cas -> omega_mod_cas (MM+MPL 20170310)
      d_t_dyn_z(:)=omega_mod_cas(:)*d_t_z(:)
      d_th_dyn_z(:)=omega_mod_cas(:)*d_th_z(:)
      d_q_dyn_z(:)=omega_mod_cas(:)*d_q_z(:)
      d_u_dyn_z(:)=omega_mod_cas(:)*d_u_z(:)
      d_v_dyn_z(:)=omega_mod_cas(:)*d_v_z(:)

!geostrophic wind 
      if (forc_geo.eq.1) then
        do l=1,llm
        ug(l) = ug_mod_cas(l) 
        vg(l) = vg_mod_cas(l) 
        enddo
      endif
!wind nudging 
      if (nudging_u.gt.0.) then
        do l=1,llm
           u(l)=u(l)+timestep*(u_mod_cas(l)-u(l))/(nudge_u)
        enddo
!     else
!       do l=1,llm
!          u(l) = u_mod_cas(l) 
!       enddo
      endif

      if (nudging_v.gt.0.) then
        do l=1,llm
           v(l)=v(l)+timestep*(v_mod_cas(l)-v(l))/(nudge_v)
        enddo
!     else
!       do l=1,llm
!          v(l) = v_mod_cas(l) 
!       enddo
      endif

      if (nudging_w.gt.0.) then
        do l=1,llm
           w(l)=w(l)+timestep*(w_mod_cas(l)-w(l))/(nudge_w)
        enddo
 !    else
 !      do l=1,llm
 !         w(l) = w_mod_cas(l) 
 !      enddo
      endif

!nudging of q and temp 
      if (nudging_t.gt.0.) then
        do l=1,llm
           temp(l)=temp(l)+timestep*(t_mod_cas(l)-temp(l))/(nudge_t)
        enddo
      endif
      if (nudging_qv.gt.0.) then
        do l=1,llm
           q(l,1)=q(l,1)+timestep*(q_mod_cas(l)-q(l,1))/(nudge_q)
        enddo
      endif

      do l = 1, llm
! Modif w_mod_cas -> omega_mod_cas (MM+MPL 20170309)
       omega(l) = omega_mod_cas(l)
       omega2(l)= omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq
       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)

!calcul advections
        if ((forc_u.eq.1).and.(forc_w.eq.0)) then
           d_u_adv(l)=du_mod_cas(l)
        else if ((forc_u.eq.1).and.(forc_w.eq.1)) then
           d_u_adv(l)=hu_mod_cas(l)-d_u_dyn_z(l)
        endif

        if ((forc_v.eq.1).and.(forc_w.eq.0)) then
           d_v_adv(l)=dv_mod_cas(l)
        else if ((forc_v.eq.1).and.(forc_w.eq.1)) then
           d_v_adv(l)=hv_mod_cas(l)-d_v_dyn_z(l)
        endif

! Puisque dth a ete converti en dt, on traite de la meme facon 
! les flags tadv et thadv
        if ((tadv.eq.1.or.thadv.eq.1) .and. (forc_w.eq.0)) then
!          d_t_adv(l)=alpha*omega(l)/rcpd-dt_mod_cas(l)
           d_t_adv(l)=alpha*omega(l)/rcpd+dt_mod_cas(l)
        else if ((tadv.eq.1.or.thadv.eq.1) .and. (forc_w.eq.1)) then
!          d_t_adv(l)=alpha*omega(l)/rcpd-ht_mod_cas(l)-d_t_dyn_z(l)
           d_t_adv(l)=alpha*omega(l)/rcpd+ht_mod_cas(l)-d_t_dyn_z(l)
        endif

!       if ((thadv.eq.1) .and. (forc_w.eq.0)) then
!          d_t_adv(l)=alpha*omega(l)/rcpd-dth_mod_cas(l)
!          d_t_adv(l)=alpha*omega(l)/rcpd+dth_mod_cas(l)
!       else if ((thadv.eq.1) .and. (forc_w.eq.1)) then
!          d_t_adv(l)=alpha*omega(l)/rcpd-hth_mod_cas(l)-d_t_dyn_z(l)
!          d_t_adv(l)=alpha*omega(l)/rcpd+hth_mod_cas(l)-d_t_dyn_z(l)
!       endif

        if ((qadv.eq.1) .and. (forc_w.eq.0)) then
           d_q_adv(l,1)=dq_mod_cas(l)
!          d_q_adv(l,1)=-1*dq_mod_cas(l)
        else if ((qadv.eq.1) .and. (forc_w.eq.1)) then
           d_q_adv(l,1)=hq_mod_cas(l)-d_q_dyn_z(l)
!          d_q_adv(l,1)=-1*hq_mod_cas(l)-d_q_dyn_z(l)
        endif
         
        if (trad.eq.1) then
           tend_rayo=1
           dt_cooling(l) = dtrad_mod_cas(l)
!          print *,'dt_cooling=',dt_cooling(l)
        else
           dt_cooling(l) = 0.0
        endif
      enddo

! Faut-il multiplier par -1 ? (MPL 20160713)
      IF(ok_flux_surf) THEN
       fsens=-1.*sens_prof_cas
       flat=-1.*lat_prof_cas
       print *,'1D_interp: sens,flat',fsens,flat
      ENDIF
!
      IF (ok_prescr_ust) THEN
       ust=ustar_prof_cas
       print *,'ust=',ust
      ENDIF
      endif ! forcing_case2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

