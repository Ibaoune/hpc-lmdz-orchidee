!
! $Id: 1D_read_forc_cases.h 3537 2019-06-19 08:29:16Z fhourdin $
!
!----------------------------------------------------------------------
! forcing_les = .T. : Impose a constant cooling 
! forcing_radconv = .T. : Pure radiative-convective equilibrium:
!----------------------------------------------------------------------


      nq1=0
      nq2=0

      if (forcing_les .or. forcing_radconv                                      &
     &    .or. forcing_GCSSold .or. forcing_fire) then

      if (forcing_fire) then
!----------------------------------------------------------------------
!read fire forcings from fire.nc
!----------------------------------------------------------------------
      fich_fire='fire.nc'
      call read_fire(fich_fire,nlev_fire,nt_fire                                &
     &     ,height,tttprof,qtprof,uprof,vprof,e12prof                           &
     &     ,ugprof,vgprof,wfls,dqtdxls                                          &
     &     ,dqtdyls,dqtdtls,thlpcar)
      write(*,*) 'Forcing FIRE lu'
      kmax=120            ! nombre de niveaux dans les profils et forcages
      else 
!----------------------------------------------------------------------
! Read profiles from files: prof.inp.001 and lscale.inp.001
! (repris de readlesfiles)
!----------------------------------------------------------------------

      call readprofiles(nlev_max,kmax,nqtot,height,                             &
     &           tttprof,qtprof,uprof,vprof,                                    &
     &           e12prof,ugprof,vgprof,                                         &
     &           wfls,dqtdxls,dqtdyls,dqtdtls,                                  &
     &           thlpcar,qprof,nq1,nq2)
      endif

! compute altitudes of play levels.
      zlay(1) =zsurf +  rd*tsurf*(psurf-play(1))/(rg*psurf)
      do l = 2,llm
        zlay(l) = zlay(l-1)+rd*tsurf*(psurf-play(1))/(rg*psurf)
      enddo

!----------------------------------------------------------------------
! Interpolation of the profiles given on the input file to
! model levels
!----------------------------------------------------------------------
      zlay(1) = zsurf +  rd*tsurf*(psurf-play(1))/(rg*psurf)
      do l=1,llm
        ! Above the max altutide of the input file

        if (zlay(l)<height(kmax)) mxcalc=l

        frac = (height(kmax)-zlay(l))/(height (kmax)-height(kmax-1))
        ttt =tttprof(kmax)-frac*(tttprof(kmax)-tttprof(kmax-1))
       if ((forcing_GCSSold .AND. tp_ini_GCSSold) .OR. forcing_fire)then ! pot. temp. in initial profile
          temp(l) = ttt*(play(l)/pzero)**rkappa
          teta(l) = ttt
       else
          temp(l) = ttt
          teta(l) = ttt*(pzero/play(l))**rkappa
       endif
          print *,' temp,teta ',l,temp(l),teta(l)
        q(l,1)  = qtprof(kmax)-frac*( qtprof(kmax)- qtprof(kmax-1))
        u(l)    =  uprof(kmax)-frac*(  uprof(kmax)-  uprof(kmax-1))
        v(l)    =  vprof(kmax)-frac*(  vprof(kmax)-  vprof(kmax-1))
        ug(l)   = ugprof(kmax)-frac*( ugprof(kmax)- ugprof(kmax-1))
        vg(l)   = vgprof(kmax)-frac*( vgprof(kmax)- vgprof(kmax-1))
        IF (nq2>0) q(l,nq1:nq2)=qprof(kmax,nq1:nq2)                         &
     &               -frac*(qprof(kmax,nq1:nq2)-qprof(kmax-1,nq1:nq2))
        omega(l)=   wfls(kmax)-frac*(   wfls(kmax)-   wfls(kmax-1))

        dq_dyn(l,1) = dqtdtls(kmax)-frac*(dqtdtls(kmax)-dqtdtls(kmax-1))
        dt_cooling(l)=thlpcar(kmax)-frac*(thlpcar(kmax)-thlpcar(kmax-1))
        do k=2,kmax
          print *,'k l height(k) height(k-1) zlay(l) frac=',k,l,height(k),height(k-1),zlay(l),frac
          frac = (height(k)-zlay(l))/(height(k)-height(k-1))
          if(l==1) print*,'k, height, tttprof',k,height(k),tttprof(k)
          if(zlay(l)>height(k-1).and.zlay(l)<height(k)) then
            ttt =tttprof(k)-frac*(tttprof(k)-tttprof(k-1))
       if ((forcing_GCSSold .AND. tp_ini_GCSSold) .OR. forcing_fire)then ! pot. temp. in initial profile
          temp(l) = ttt*(play(l)/pzero)**rkappa
          teta(l) = ttt
       else
          temp(l) = ttt
          teta(l) = ttt*(pzero/play(l))**rkappa
       endif
          print *,' temp,teta ',l,temp(l),teta(l)
            q(l,1)  = qtprof(k)-frac*( qtprof(k)- qtprof(k-1))
            u(l)    =  uprof(k)-frac*(  uprof(k)-  uprof(k-1))
            v(l)    =  vprof(k)-frac*(  vprof(k)-  vprof(k-1))
            ug(l)   = ugprof(k)-frac*( ugprof(k)- ugprof(k-1))
            vg(l)   = vgprof(k)-frac*( vgprof(k)- vgprof(k-1))
            IF (nq2>0) q(l,nq1:nq2)=qprof(k,nq1:nq2)                        &
     &                   -frac*(qprof(k,nq1:nq2)-qprof(k-1,nq1:nq2))
            omega(l)=   wfls(k)-frac*(   wfls(k)-   wfls(k-1))
            dq_dyn(l,1)=dqtdtls(k)-frac*(dqtdtls(k)-dqtdtls(k-1))
            dt_cooling(l)=thlpcar(k)-frac*(thlpcar(k)-thlpcar(k-1))
          elseif(zlay(l)<height(1)) then ! profils uniformes pour z<height(1)
            ttt =tttprof(1)
       if ((forcing_GCSSold .AND. tp_ini_GCSSold) .OR. forcing_fire)then ! pot. temp. in initial profile
          temp(l) = ttt*(play(l)/pzero)**rkappa
          teta(l) = ttt
       else
          temp(l) = ttt
          teta(l) = ttt*(pzero/play(l))**rkappa
       endif
            q(l,1)  = qtprof(1)
            u(l)    =  uprof(1)
            v(l)    =  vprof(1)
            ug(l)   = ugprof(1)
            vg(l)   = vgprof(1)
            omega(l)=   wfls(1)
            IF (nq2>0) q(l,nq1:nq2)=qprof(1,nq1:nq2)
            dq_dyn(l,1)  =dqtdtls(1)
            dt_cooling(l)=thlpcar(1)
          endif
        enddo 

        temp(l)=max(min(temp(l),350.),150.)
        rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
        if (l .lt. llm) then
          zlay(l+1) = zlay(l) + (play(l)-play(l+1))/(rg*rho(l))
        endif
        omega2(l)=-rho(l)*omega(l)
        omega(l)= omega(l)*(-rg*rho(l)) !en Pa/s
        if (l>1) then
        if(zlay(l-1)>height(kmax)) then
           omega(l)=0.0
           omega2(l)=0.0
        endif   
        endif
        if(q(l,1)<0.) q(l,1)=0.0
        q(l,2)  = 0.0
      enddo

      endif ! forcing_les .or. forcing_GCSSold .or. forcing_fire
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing for GCSSold:
!---------------------------------------------------------------------
      if (forcing_GCSSold) then
       fich_gcssold_ctl = './forcing.ctl'
       fich_gcssold_dat = './forcing8.dat'
       call copie(llm,play,psurf,fich_gcssold_ctl)
       call get_uvd2(it,timestep,fich_gcssold_ctl,fich_gcssold_dat,         &
     &               ht_gcssold,hq_gcssold,hw_gcssold,                      &
     &               hu_gcssold,hv_gcssold,                                 &
     &               hthturb_gcssold,hqturb_gcssold,Ts_gcssold,             &
     &               imp_fcg_gcssold,ts_fcg_gcssold,                        &
     &               Tp_fcg_gcssold,Turb_fcg_gcssold)
       print *,' get_uvd2 -> hqturb_gcssold ',hqturb_gcssold
      endif ! forcing_GCSSold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing for RICO:
!---------------------------------------------------------------------
      if (forcing_rico) then

!       call writefield_phy('omega', omega,llm+1)
      fich_rico = 'rico.txt'
       call read_rico(fich_rico,nlev_rico,ps_rico,play                      &
     &             ,ts_rico,t_rico,q_rico,u_rico,v_rico,w_rico              &
     &             ,dth_rico,dqh_rico)
        print*, ' on a lu et prepare RICO'

       mxcalc=llm
       print *, airefi, ' airefi '
       do l = 1, llm
       rho(l)  = play(l)/(rd*t_rico(l)*(1.+(rv/rd-1.)*q_rico(l)))
       temp(l) = t_rico(l)
       q(l,1) = q_rico(l)
       q(l,2) = 0.0
       u(l) = u_rico(l)
       v(l) = v_rico(l)
       ug(l)=u_rico(l)
       vg(l)=v_rico(l)
       omega(l) = -w_rico(l)*rg
       omega2(l) = omega(l)/rg*airefi
       enddo
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing from TOGA-COARE experiment (Ciesielski et al. 2002) :
!---------------------------------------------------------------------

      if (forcing_toga) then

! read TOGA-COARE forcing (native vertical grid, nt_toga timesteps):
      fich_toga = './d_toga/ifa_toga_coare_v21_dime.txt'
      CALL read_togacoare(fich_toga,nlev_toga,nt_toga                       &
     &         ,ts_toga,plev_toga,t_toga,q_toga,u_toga,v_toga,w_toga        &
     &         ,ht_toga,vt_toga,hq_toga,vq_toga)

       write(*,*) 'Forcing TOGA lu'

! time interpolation for initial conditions:
      write(*,*) 'AVT 1ere INTERPOLATION: day,day1 = ',day,day1
      CALL interp_toga_time(daytime,day1,annee_ref                          &
     &             ,year_ini_toga,day_ju_ini_toga,nt_toga,dt_toga           &
     &             ,nlev_toga,ts_toga,plev_toga,t_toga,q_toga,u_toga        &
     &             ,v_toga,w_toga,ht_toga,vt_toga,hq_toga,vq_toga           &
     &             ,ts_prof,plev_prof,t_prof,q_prof,u_prof,v_prof,w_prof    &
     &             ,ht_prof,vt_prof,hq_prof,vq_prof)

! vertical interpolation:
      CALL interp_toga_vertical(play,nlev_toga,plev_prof                    &
     &         ,t_prof,q_prof,u_prof,v_prof,w_prof                          &
     &         ,ht_prof,vt_prof,hq_prof,vq_prof                             &
     &         ,t_mod,q_mod,u_mod,v_mod,w_mod                               &
     &         ,ht_mod,vt_mod,hq_mod,vq_mod,mxcalc)
       write(*,*) 'Profil initial forcing TOGA interpole'

! initial and boundary conditions :
      tsurf = ts_prof
      write(*,*) 'SST initiale: ',tsurf
      do l = 1, llm
       temp(l) = t_mod(l)
       q(l,1) = q_mod(l)
       q(l,2) = 0.0
       u(l) = u_mod(l)
       v(l) = v_mod(l)
       omega(l) = w_mod(l)
       omega2(l)=omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq
!?       rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
!?       omega2(l)=-rho(l)*omega(l)
       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
       d_t_adv(l) = alpha*omega(l)/rcpd-(ht_mod(l)+vt_mod(l))
       d_q_adv(l,1) = -(hq_mod(l)+vq_mod(l))
       d_q_adv(l,2) = 0.0
      enddo

      endif ! forcing_toga
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing from TWPICE experiment (Shaocheng et al. 2010) :
!---------------------------------------------------------------------

      if (forcing_twpice) then
!read TWP-ICE forcings
     fich_twpice='d_twpi/twp180iopsndgvarana_v2.1_C3.c1.20060117.000000.cdf'
      call read_twpice(fich_twpice,nlev_twpi,nt_twpi                        &
     &     ,ts_twpi,plev_twpi,t_twpi,q_twpi,u_twpi,v_twpi,w_twpi            &
     &     ,ht_twpi,vt_twpi,hq_twpi,vq_twpi)

      write(*,*) 'Forcing TWP-ICE lu'
!Time interpolation for initial conditions using TOGA interpolation routine
         write(*,*) 'AVT 1ere INTERPOLATION: day,day1 = ',daytime,day1   
      CALL interp_toga_time(daytime,day1,annee_ref                          &
     &          ,year_ini_twpi,day_ju_ini_twpi,nt_twpi,dt_twpi,nlev_twpi    &
     &             ,ts_twpi,plev_twpi,t_twpi,q_twpi,u_twpi,v_twpi,w_twpi    &
     &             ,ht_twpi,vt_twpi,hq_twpi,vq_twpi                         &
     &             ,ts_proftwp,plev_proftwp,t_proftwp,q_proftwp             &
     &             ,u_proftwp,v_proftwp,w_proftwp                           &
     &             ,ht_proftwp,vt_proftwp,hq_proftwp,vq_proftwp)

! vertical interpolation using TOGA interpolation routine:
!      write(*,*)'avant interp vert', t_proftwp
      CALL interp_toga_vertical(play,nlev_twpi,plev_proftwp                 &
     &         ,t_proftwp,q_proftwp,u_proftwp,v_proftwp,w_proftwp           &
     &         ,ht_proftwp,vt_proftwp,hq_proftwp,vq_proftwp                 &
     &         ,t_mod,q_mod,u_mod,v_mod,w_mod                               &
     &         ,ht_mod,vt_mod,hq_mod,vq_mod,mxcalc)
!       write(*,*) 'Profil initial forcing TWP-ICE interpole',t_mod

! initial and boundary conditions :
!      tsurf = ts_proftwp
      write(*,*) 'SST initiale: ',tsurf
      do l = 1, llm
       temp(l) = t_mod(l)
       q(l,1) = q_mod(l)
       q(l,2) = 0.0
       u(l) = u_mod(l)
       v(l) = v_mod(l)
       omega(l) = w_mod(l)
       omega2(l)=omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!on applique le forcage total au premier pas de temps
!attention: signe different de toga
       d_t_adv(l) = alpha*omega(l)/rcpd+(ht_mod(l)+vt_mod(l))
       d_q_adv(l,1) = (hq_mod(l)+vq_mod(l))
       d_q_adv(l,2) = 0.0
      enddo     
       
      endif !forcing_twpice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing from AMMA experiment (Couvreux et al. 2010) :
!---------------------------------------------------------------------

      if (forcing_amma) then

      call read_1D_cases

      write(*,*) 'Forcing AMMA lu'

!champs initiaux:
      do k=1,nlev_amma
         th_ammai(k)=th_amma(k)
         q_ammai(k)=q_amma(k)
         u_ammai(k)=u_amma(k)
         v_ammai(k)=v_amma(k)
         vitw_ammai(k)=vitw_amma(k,12)
         ht_ammai(k)=ht_amma(k,12)
         hq_ammai(k)=hq_amma(k,12)
         vt_ammai(k)=0.
         vq_ammai(k)=0.
      enddo   
      omega(:)=0.      
      omega2(:)=0.
      rho(:)=0.
! vertical interpolation using TOGA interpolation routine:
!      write(*,*)'avant interp vert', t_proftwp
      CALL interp_toga_vertical(play,nlev_amma,plev_amma                    &
     &         ,th_ammai,q_ammai,u_ammai,v_ammai,vitw_ammai                 &
     &         ,ht_ammai,vt_ammai,hq_ammai,vq_ammai                         &
     &         ,t_mod,q_mod,u_mod,v_mod,w_mod                               &
     &         ,ht_mod,vt_mod,hq_mod,vq_mod,mxcalc)
!       write(*,*) 'Profil initial forcing TWP-ICE interpole',t_mod

! initial and boundary conditions :
!      tsurf = ts_proftwp
      write(*,*) 'SST initiale mxcalc: ',tsurf,mxcalc
      do l = 1, llm
! Ligne du dessous ?? decommenter si on lit theta au lieu de temp
!      temp(l) = t_mod(l)*(play(l)/pzero)**rkappa 
       temp(l) = t_mod(l) 
       q(l,1) = q_mod(l)
       q(l,2) = 0.0
!      print *,'read_forc: l,temp,q=',l,temp(l),q(l,1)
       u(l) = u_mod(l)
       v(l) = v_mod(l)
       rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
       omega(l) = w_mod(l)*(-rg*rho(l))
       omega2(l)=omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!on applique le forcage total au premier pas de temps
!attention: signe different de toga
       d_t_adv(l) = alpha*omega(l)/rcpd+ht_mod(l)
!forcage en th
!       d_t_adv(l) = ht_mod(l)
       d_q_adv(l,1) = hq_mod(l)
       d_q_adv(l,2) = 0.0
       dt_cooling(l)=0.
      enddo     
       write(*,*) 'Prof initeforcing AMMA interpole temp39',temp(39)
      

!     ok_flux_surf=.false.
      fsens=-1.*sens_amma(12)
      flat=-1.*lat_amma(12)
       
      endif !forcing_amma


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing from DICE experiment (see file DICE_protocol_vn2-3.pdf)
!---------------------------------------------------------------------

      if (forcing_dice) then
!read DICE forcings
      fich_dice='dice_driver.nc'
      call read_dice(fich_dice,nlev_dice,nt_dice                    &
     &     ,zz_dice,plev_dice,t_dice,qv_dice,u_dice,v_dice,o3_dice &
     &     ,shf_dice,lhf_dice,lwup_dice,swup_dice,tg_dice,ustar_dice& 
     &     ,psurf_dice,ug_dice,vg_dice,ht_dice,hq_dice              &
     &     ,hu_dice,hv_dice,w_dice,omega_dice)

      write(*,*) 'Forcing DICE lu'

!champs initiaux:
      do k=1,nlev_dice
         t_dicei(k)=t_dice(k)
         qv_dicei(k)=qv_dice(k)
         u_dicei(k)=u_dice(k)
         v_dicei(k)=v_dice(k)
         o3_dicei(k)=o3_dice(k)
         ht_dicei(k)=ht_dice(k,1)
         hq_dicei(k)=hq_dice(k,1)
         hu_dicei(k)=hu_dice(k,1)
         hv_dicei(k)=hv_dice(k,1)
         w_dicei(k)=w_dice(k,1)
         omega_dicei(k)=omega_dice(k,1)
      enddo   
      omega(:)=0.      
      omega2(:)=0.
      rho(:)=0.
! vertical interpolation using TOGA interpolation routine:
!      write(*,*)'avant interp vert', t_proftwp
!
!     CALL interp_dice_time(daytime,day1,annee_ref
!    i             ,year_ini_dice,day_ju_ini_dice,nt_dice,dt_dice
!    i             ,nlev_dice,shf_dice,lhf_dice,lwup_dice,swup_dice
!    i             ,tg_dice,ustar_dice,psurf_dice,ug_dice,vg_dice
!    i             ,ht_dice,hq_dice,hu_dice,hv_dice,w_dice,omega_dice
!    o             ,shf_prof,lhf_prof,lwup_prof,swup_prof,tg_prof
!    o             ,ustar_prof,psurf_prof,ug_profd,vg_profd
!    o             ,ht_profd,hq_profd,hu_profd,hv_profd,w_profd
!    o             ,omega_profd)

      CALL interp_dice_vertical(play,nlev_dice,nt_dice,plev_dice       &
     &         ,t_dicei,qv_dicei,u_dicei,v_dicei,o3_dicei             &
     &         ,ht_dicei,hq_dicei,hu_dicei,hv_dicei,w_dicei,omega_dicei&
     &         ,t_mod,qv_mod,u_mod,v_mod,o3_mod                       &
     &         ,ht_mod,hq_mod,hu_mod,hv_mod,w_mod,omega_mod,mxcalc)

! Pour tester les advections horizontales de T et Q, on met w_mod et omega_mod ?? zero (MPL 20131108)
!     w_mod(:,:)=0.
!     omega_mod(:,:)=0.

!       write(*,*) 'Profil initial forcing DICE interpole',t_mod
! Les forcages DICE sont donnes /jour et non /seconde !
      ht_mod(:)=ht_mod(:)/86400.
      hq_mod(:)=hq_mod(:)/86400.
      hu_mod(:)=hu_mod(:)/86400.
      hv_mod(:)=hv_mod(:)/86400.

! initial and boundary conditions :
      write(*,*) 'SST initiale mxcalc: ',tsurf,mxcalc
      do l = 1, llm
! Ligne du dessous ?? decommenter si on lit theta au lieu de temp
!      temp(l) = th_mod(l)*(play(l)/pzero)**rkappa 
       temp(l) = t_mod(l) 
       q(l,1) = qv_mod(l)
       q(l,2) = 0.0
!      print *,'read_forc: l,temp,q=',l,temp(l),q(l,1)
       u(l) = u_mod(l)
       v(l) = v_mod(l)
       ug(l)=ug_dice(1)
       vg(l)=vg_dice(1)
       rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
!      omega(l) = w_mod(l)*(-rg*rho(l))
       omega(l) = omega_mod(l)
       omega2(l)=omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!on applique le forcage total au premier pas de temps
!attention: signe different de toga
       d_t_adv(l) = alpha*omega(l)/rcpd+ht_mod(l)
!forcage en th
!       d_t_adv(l) = ht_mod(l)
       d_q_adv(l,1) = hq_mod(l)
       d_q_adv(l,2) = 0.0
       dt_cooling(l)=0.
      enddo     
       write(*,*) 'Profil initial forcing DICE interpole temp39',temp(39)
      

!     ok_flux_surf=.false.
      fsens=-1.*shf_dice(1)
      flat=-1.*lhf_dice(1)
! Le cas Dice doit etre force avec ustar mais on peut simplifier en forcant par 
! le coefficient de trainee en surface cd**2=ustar*vent(k=1)
! On commence ici a stocker ustar dans cdrag puis on terminera le calcul dans pbl_surface
! MPL 05082013
      ust=ustar_dice(1)
      tg=tg_dice(1)
      print *,'ust= ',ust
      IF (tsurf .LE. 0.) THEN
       tsurf= tg_dice(1)
      ENDIF
      psurf= psurf_dice(1)
      solsw_in = (1.-albedo)/albedo*swup_dice(1)
      sollw_in = (0.7*RSIGMA*temp(1)**4)-lwup_dice(1)
      PRINT *,'1D_READ_FORC : solsw, sollw',solsw_in,sollw_in
      endif !forcing_dice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing from GABLS4 experiment
!---------------------------------------------------------------------

!!!! Si la temperature de surface n'est pas impos??e:
  
      if (forcing_gabls4) then
!read GABLS4 forcings
      
      fich_gabls4='gabls4_driver.nc'
      
       	
      call read_gabls4(fich_gabls4,nlev_gabls4,nt_gabls4,nsol_gabls4,zz_gabls4,depth_sn_gabls4,ug_gabls4,vg_gabls4 &
     & ,plev_gabls4,th_gabls4,t_gabls4,qv_gabls4,u_gabls4,v_gabls4,ht_gabls4,hq_gabls4,tg_gabls4,tsnow_gabls4,snow_dens_gabls4)

      write(*,*) 'Forcing GABLS4 lu'

!champs initiaux:
      do k=1,nlev_gabls4
         t_gabi(k)=t_gabls4(k)
         qv_gabi(k)=qv_gabls4(k)
         u_gabi(k)=u_gabls4(k)
         v_gabi(k)=v_gabls4(k)
         poub(k)=0.
         ht_gabi(k)=ht_gabls4(k,1)
         hq_gabi(k)=hq_gabls4(k,1)
         ug_gabi(k)=ug_gabls4(k,1)
         vg_gabi(k)=vg_gabls4(k,1)
      enddo 
  
      omega(:)=0.      
      omega2(:)=0.
      rho(:)=0.
! vertical interpolation using TOGA interpolation routine:
!      write(*,*)'avant interp vert', t_proftwp
!
!     CALL interp_dice_time(daytime,day1,annee_ref
!    i             ,year_ini_dice,day_ju_ini_dice,nt_dice,dt_dice
!    i             ,nlev_dice,shf_dice,lhf_dice,lwup_dice,swup_dice
!    i             ,tg_dice,ustar_dice,psurf_dice,ug_dice,vg_dice
!    i             ,ht_dice,hq_dice,hu_dice,hv_dice,w_dice,omega_dice
!    o             ,shf_prof,lhf_prof,lwup_prof,swup_prof,tg_prof
!    o             ,ustar_prof,psurf_prof,ug_profd,vg_profd
!    o             ,ht_profd,hq_profd,hu_profd,hv_profd,w_profd
!    o             ,omega_profd)

      CALL interp_dice_vertical(play,nlev_gabls4,nt_gabls4,plev_gabls4       &
     &         ,t_gabi,qv_gabi,u_gabi,v_gabi,poub                  &
     &         ,ht_gabi,hq_gabi,ug_gabi,vg_gabi,poub,poub          &
     &         ,t_mod,qv_mod,u_mod,v_mod,o3_mod                    &
     &         ,ht_mod,hq_mod,ug_mod,vg_mod,w_mod,omega_mod,mxcalc)

! Les forcages GABLS4 ont l air d etre en K/S quoiqu en dise le fichier gabls4_driver.nc !? MPL 20141024
!     ht_mod(:)=ht_mod(:)/86400.
!     hq_mod(:)=hq_mod(:)/86400.

! initial and boundary conditions :
      write(*,*) 'SST initiale mxcalc: ',tsurf,mxcalc
      do l = 1, llm
! Ligne du dessous ?? decommenter si on lit theta au lieu de temp
!      temp(l) = th_mod(l)*(play(l)/pzero)**rkappa 
       temp(l) = t_mod(l) 
       q(l,1) = qv_mod(l)
       q(l,2) = 0.0
!      print *,'read_forc: l,temp,q=',l,temp(l),q(l,1)
       u(l) = u_mod(l)
       v(l) = v_mod(l)
       ug(l)=ug_mod(l)
       vg(l)=vg_mod(l)
       
!
!	tg=tsurf
!       

       print *,'***** tsurf=',tsurf
       rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
!      omega(l) = w_mod(l)*(-rg*rho(l))
       omega(l) = omega_mod(l)
       omega2(l)=omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq
       
    

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!on applique le forcage total au premier pas de temps
!attention: signe different de toga
!      d_t_adv(l) = alpha*omega(l)/rcpd+ht_mod(l)
!forcage en th
       d_t_adv(l) = ht_mod(l)
       d_q_adv(l,1) = hq_mod(l)
       d_q_adv(l,2) = 0.0
       dt_cooling(l)=0.
      enddo     

!--------------- Residus forcages du cas Dice (a supprimer) MPL 20141024---------------
! Le cas Dice doit etre force avec ustar mais on peut simplifier en forcant par 
! le coefficient de trainee en surface cd**2=ustar*vent(k=1)
! On commence ici a stocker ustar dans cdrag puis on terminera le calcul dans pbl_surface
! MPL 05082013
!     ust=ustar_dice(1)
!     tg=tg_dice(1)
!     print *,'ust= ',ust
!     IF (tsurf .LE. 0.) THEN
!      tsurf= tg_dice(1)
!     ENDIF
!     psurf= psurf_dice(1)
!     solsw_in = (1.-albedo)/albedo*swup_dice(1)
!     sollw_in = (0.7*RSIGMA*temp(1)**4)-lwup_dice(1)
!     PRINT *,'1D_READ_FORC : solsw, sollw',solsw_in,sollw_in
!--------------------------------------------------------------------------------------
      endif !forcing_gabls4



! Forcing from Arm_Cu case                   
! For this case, ifa_armcu.txt contains sensible, latent heat fluxes
! large scale advective forcing,radiative forcing 
! and advective tendency of theta and qt to be applied 
!---------------------------------------------------------------------

      if (forcing_armcu) then
! read armcu forcing :
       write(*,*) 'Avant lecture Forcing Arm_Cu'
      fich_armcu = './ifa_armcu.txt'
      CALL read_armcu(fich_armcu,nlev_armcu,nt_armcu,                       &
     & sens_armcu,flat_armcu,adv_theta_armcu,                               &
     & rad_theta_armcu,adv_qt_armcu)
       write(*,*) 'Forcing Arm_Cu lu'

!----------------------------------------------------------------------
! Read profiles from file: prof.inp.19 or prof.inp.40 
! For this case, profiles are given for two vertical resolution
! 19 or 40 levels
!
! Comment from: http://www.knmi.nl/samenw/eurocs/ARM/profiles.html
! Note that the initial profiles contain no liquid water! 
! (so potential temperature can be interpreted as liquid water 
! potential temperature and water vapor as total water)
! profiles are given at full levels
!----------------------------------------------------------------------

      call readprofile_armcu(nlev_max,kmax,height,play_mod,u_mod,           &
     &           v_mod,theta_mod,t_mod,qv_mod,rv_mod,ap,bp)

! time interpolation for initial conditions:
      write(*,*) 'AVT 1ere INTERPOLATION: day,day1 = ',day,day1

      print *,'Avant interp_armcu_time'
      print *,'daytime=',daytime
      print *,'day1=',day1
      print *,'annee_ref=',annee_ref
      print *,'year_ini_armcu=',year_ini_armcu
      print *,'day_ju_ini_armcu=',day_ju_ini_armcu
      print *,'nt_armcu=',nt_armcu
      print *,'dt_armcu=',dt_armcu
      print *,'nlev_armcu=',nlev_armcu
      CALL interp_armcu_time(daytime,day1,annee_ref                         &
     &            ,year_ini_armcu,day_ju_ini_armcu,nt_armcu,dt_armcu        &
     &            ,nlev_armcu,sens_armcu,flat_armcu,adv_theta_armcu         &
     &            ,rad_theta_armcu,adv_qt_armcu,sens_prof,flat_prof         &
     &            ,adv_theta_prof,rad_theta_prof,adv_qt_prof)
       write(*,*) 'Forcages interpoles dans temps'

! No vertical interpolation if nlev imposed to 19 or 40
! The vertical grid stops at 4000m # 600hPa
      mxcalc=llm

! initial and boundary conditions :
!     tsurf = ts_prof
! tsurf read in lmdz1d.def
      write(*,*) 'Tsurf initiale: ',tsurf
      do l = 1, llm
       play(l)=play_mod(l)*100.
       presnivs(l)=play(l)
       zlay(l)=height(l)
       temp(l) = t_mod(l)
       teta(l)=theta_mod(l)
       q(l,1) = qv_mod(l)/1000.
! No liquid water in the initial profil
       q(l,2) = 0.
       u(l) = u_mod(l)
       ug(l)= u_mod(l)
       v(l) = v_mod(l)
       vg(l)= v_mod(l)
! Advective forcings are given in K or g/kg ... per HOUR
!      IF(height(l).LT.1000) THEN
!        d_t_adv(l) = (adv_theta_prof + rad_theta_prof)/3600.
!        d_q_adv(l,1) = adv_qt_prof/1000./3600.
!        d_q_adv(l,2) = 0.0
!      ELSEIF (height(l).GE.1000.AND.height(l).LT.3000) THEN
!        d_t_adv(l) = (adv_theta_prof + rad_theta_prof)*
!    :               (1-(height(l)-1000.)/2000.)
!        d_t_adv(l) = d_t_adv(l)/3600.
!        d_q_adv(l,1) = adv_qt_prof*(1-(height(l)-1000.)/2000.)
!        d_q_adv(l,1) = d_q_adv(l,1)/1000./3600.
!        d_q_adv(l,2) = 0.0
!      ELSE
!        d_t_adv(l) = 0.0
!        d_q_adv(l,1) = 0.0
!        d_q_adv(l,2) = 0.0
!      ENDIF
      enddo
! plev at half levels is given in proh.inp.19 or proh.inp.40 files
      plev(1)= ap(llm+1)+bp(llm+1)*psurf
      do l = 1, llm
      plev(l+1) = ap(llm-l+1)+bp(llm-l+1)*psurf
      print *,'Read_forc: l height play plev zlay temp',                    &
     &   l,height(l),play(l),plev(l),zlay(l),temp(l)
      enddo
! For this case, fluxes are imposed
       fsens=-1*sens_prof
       flat=-1*flat_prof

      endif ! forcing_armcu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing from transition case of Irina Sandu                  
!---------------------------------------------------------------------

      if (forcing_sandu) then
       write(*,*) 'Avant lecture Forcing SANDU'

! read sanduref forcing :
      fich_sandu = './ifa_sanduref.txt'
      CALL read_sandu(fich_sandu,nlev_sandu,nt_sandu,ts_sandu)

       write(*,*) 'Forcing SANDU lu'

!----------------------------------------------------------------------
! Read profiles from file: prof.inp.001 
!----------------------------------------------------------------------

      call readprofile_sandu(nlev_max,kmax,height,plev_profs,t_profs,       &
     &           thl_profs,q_profs,u_profs,v_profs,                         &
     &           w_profs,omega_profs,o3mmr_profs)

! time interpolation for initial conditions:
      write(*,*) 'AVT 1ere INTERPOLATION: day,day1 = ',day,day1
! ATTENTION, cet appel ne convient pas pour le cas SANDU !!
! revoir 1DUTILS.h et les arguments

      print *,'Avant interp_sandu_time'
      print *,'daytime=',daytime
      print *,'day1=',day1
      print *,'annee_ref=',annee_ref
      print *,'year_ini_sandu=',year_ini_sandu
      print *,'day_ju_ini_sandu=',day_ju_ini_sandu
      print *,'nt_sandu=',nt_sandu
      print *,'dt_sandu=',dt_sandu
      print *,'nlev_sandu=',nlev_sandu
      CALL interp_sandu_time(daytime,day1,annee_ref                         &
     &             ,year_ini_sandu,day_ju_ini_sandu,nt_sandu,dt_sandu       &
     &             ,nlev_sandu                                              &
     &             ,ts_sandu,ts_prof)

! vertical interpolation:
      print *,'Avant interp_vertical: nlev_sandu=',nlev_sandu
      CALL interp_sandu_vertical(play,nlev_sandu,plev_profs                 &
     &         ,t_profs,thl_profs,q_profs,u_profs,v_profs,w_profs           &
     &         ,omega_profs,o3mmr_profs                                     &
     &         ,t_mod,thl_mod,q_mod,u_mod,v_mod,w_mod                       &
     &         ,omega_mod,o3mmr_mod,mxcalc)
       write(*,*) 'Profil initial forcing SANDU interpole'

! initial and boundary conditions :
      tsurf = ts_prof
      write(*,*) 'SST initiale: ',tsurf
      do l = 1, llm
       temp(l) = t_mod(l)
       tetal(l)=thl_mod(l)
       q(l,1) = q_mod(l)
       q(l,2) = 0.0
       u(l) = u_mod(l)
       v(l) = v_mod(l)
       w(l) = w_mod(l)
       omega(l) = omega_mod(l)
       omega2(l)=omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq
!?       rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
!?       omega2(l)=-rho(l)*omega(l)
       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!      d_t_adv(l) = alpha*omega(l)/rcpd+vt_mod(l)
!      d_q_adv(l,1) = vq_mod(l)
       d_t_adv(l) = alpha*omega(l)/rcpd
       d_q_adv(l,1) = 0.0
       d_q_adv(l,2) = 0.0
      enddo

      endif ! forcing_sandu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing from Astex case
!---------------------------------------------------------------------

      if (forcing_astex) then
       write(*,*) 'Avant lecture Forcing Astex'

! read astex forcing :
      fich_astex = './ifa_astex.txt'
      CALL read_astex(fich_astex,nlev_astex,nt_astex,div_astex,ts_astex,    &
     &  ug_astex,vg_astex,ufa_astex,vfa_astex)

       write(*,*) 'Forcing Astex lu'

!----------------------------------------------------------------------
! Read profiles from file: prof.inp.001
!----------------------------------------------------------------------

      call readprofile_astex(nlev_max,kmax,height,plev_profa,t_profa,       &
     &           thl_profa,qv_profa,ql_profa,qt_profa,u_profa,v_profa,      &
     &           w_profa,tke_profa,o3mmr_profa)

! time interpolation for initial conditions:
      write(*,*) 'AVT 1ere INTERPOLATION: day,day1 = ',day,day1
! ATTENTION, cet appel ne convient pas pour le cas SANDU !!
! revoir 1DUTILS.h et les arguments

      print *,'Avant interp_astex_time'
      print *,'daytime=',daytime
      print *,'day1=',day1
      print *,'annee_ref=',annee_ref
      print *,'year_ini_astex=',year_ini_astex
      print *,'day_ju_ini_astex=',day_ju_ini_astex
      print *,'nt_astex=',nt_astex
      print *,'dt_astex=',dt_astex
      print *,'nlev_astex=',nlev_astex
      CALL interp_astex_time(daytime,day1,annee_ref                         &
     &             ,year_ini_astex,day_ju_ini_astex,nt_astex,dt_astex       &
     &             ,nlev_astex,div_astex,ts_astex,ug_astex,vg_astex         &
     &             ,ufa_astex,vfa_astex,div_prof,ts_prof,ug_prof,vg_prof    &
     &             ,ufa_prof,vfa_prof)

! vertical interpolation:
      print *,'Avant interp_vertical: nlev_astex=',nlev_astex
      CALL interp_astex_vertical(play,nlev_astex,plev_profa                 &
     &         ,t_profa,thl_profa,qv_profa,ql_profa,qt_profa                &
     &         ,u_profa,v_profa,w_profa,tke_profa,o3mmr_profa               &
     &         ,t_mod,thl_mod,qv_mod,ql_mod,qt_mod,u_mod,v_mod,w_mod        &
     &         ,tke_mod,o3mmr_mod,mxcalc)
       write(*,*) 'Profil initial forcing Astex interpole'

! initial and boundary conditions :
      tsurf = ts_prof
      write(*,*) 'SST initiale: ',tsurf
      do l = 1, llm
       temp(l) = t_mod(l)
       tetal(l)=thl_mod(l)
       q(l,1) = qv_mod(l)
       q(l,2) = ql_mod(l)
       u(l) = u_mod(l)
       v(l) = v_mod(l)
       w(l) = w_mod(l)
       omega(l) = w_mod(l)
!      omega2(l)=omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq
!      rho(l)  = play(l)/(rd*temp(l)*(1.+(rv/rd-1.)*q(l,1)))
!      omega2(l)=-rho(l)*omega(l)
       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!      d_t_adv(l) = alpha*omega(l)/rcpd+vt_mod(l)
!      d_q_adv(l,1) = vq_mod(l)
       d_t_adv(l) = alpha*omega(l)/rcpd
       d_q_adv(l,1) = 0.0
       d_q_adv(l,2) = 0.0
      enddo

      endif ! forcing_astex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing from standard case :
!---------------------------------------------------------------------

      if (forcing_case) then

         write(*,*),'avant call read_1D_cas'
         call read_1D_cas
         write(*,*) 'Forcing read'

!Time interpolation for initial conditions using TOGA interpolation routine
         write(*,*) 'AVT 1ere INTERPOLATION: day,day1 = ',daytime,day1   
      CALL interp_case_time(day,day1,annee_ref                                                              &
!    &         ,year_ini_cas,day_ju_ini_cas,nt_cas,pdt_cas,nlev_cas                                         &
     &         ,nt_cas,nlev_cas                                                                             &
     &         ,ts_cas,plev_cas,t_cas,q_cas,u_cas,v_cas                                                     &
     &         ,ug_cas,vg_cas,vitw_cas,du_cas,hu_cas,vu_cas                                                 &
     &         ,dv_cas,hv_cas,vv_cas,dt_cas,ht_cas,vt_cas,dtrad_cas                                         &
     &         ,dq_cas,hq_cas,vq_cas,lat_cas,sens_cas,ustar_cas                                             &
     &         ,uw_cas,vw_cas,q1_cas,q2_cas                                                                 &
     &         ,ts_prof_cas,plev_prof_cas,t_prof_cas,q_prof_cas,u_prof_cas,v_prof_cas                       &
     &         ,ug_prof_cas,vg_prof_cas,vitw_prof_cas,du_prof_cas,hu_prof_cas,vu_prof_cas                   &
     &         ,dv_prof_cas,hv_prof_cas,vv_prof_cas,dt_prof_cas,ht_prof_cas,vt_prof_cas,dtrad_prof_cas      &
     &         ,dq_prof_cas,hq_prof_cas,vq_prof_cas,lat_prof_cas,sens_prof_cas,ustar_prof_cas               &
     &         ,uw_prof_cas,vw_prof_cas,q1_prof_cas,q2_prof_cas)

! vertical interpolation using TOGA interpolation routine:
!      write(*,*)'avant interp vert', t_prof
      CALL interp_case_vertical(play,nlev_cas,plev_prof_cas            &
     &         ,t_prof_cas,q_prof_cas,u_prof_cas,v_prof_cas,ug_prof_cas,vg_prof_cas,vitw_prof_cas    &
     &         ,du_prof_cas,hu_prof_cas,vu_prof_cas,dv_prof_cas,hv_prof_cas,vv_prof_cas           &
     &         ,dt_prof_cas,ht_prof_cas,vt_prof_cas,dtrad_prof_cas,dq_prof_cas,hq_prof_cas,vq_prof_cas           &
     &         ,t_mod_cas,q_mod_cas,u_mod_cas,v_mod_cas,ug_mod_cas,vg_mod_cas,w_mod_cas           &
     &         ,du_mod_cas,hu_mod_cas,vu_mod_cas,dv_mod_cas,hv_mod_cas,vv_mod_cas               &
     &         ,dt_mod_cas,ht_mod_cas,vt_mod_cas,dtrad_mod_cas,dq_mod_cas,hq_mod_cas,vq_mod_cas,mxcalc)
!       write(*,*) 'Profil initial forcing case interpole',t_mod

! initial and boundary conditions :
!      tsurf = ts_prof_cas
      ts_cur = ts_prof_cas 
      psurf=plev_prof_cas(1)
      write(*,*) 'SST initiale: ',tsurf
      do l = 1, llm
       temp(l) = t_mod_cas(l)
       q(l,1) = q_mod_cas(l)
       q(l,2) = 0.0
       u(l) = u_mod_cas(l)
       v(l) = v_mod_cas(l)
       omega(l) = w_mod_cas(l)
       omega2(l)=omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!on applique le forcage total au premier pas de temps
!attention: signe different de toga
       d_t_adv(l) = alpha*omega(l)/rcpd+(ht_mod_cas(l)+vt_mod_cas(l))
       d_q_adv(l,1) = (hq_mod_cas(l)+vq_mod_cas(l))
       d_q_adv(l,2) = 0.0
       d_u_adv(l) = (hu_mod_cas(l)+vu_mod_cas(l))
! correction bug d_u -> d_v (MM+MPL 20170310)
       d_v_adv(l) = (hv_mod_cas(l)+vv_mod_cas(l))
      enddo     

! In case fluxes are imposed
       IF (ok_flux_surf) THEN
       fsens=sens_prof_cas
       flat=lat_prof_cas
       ENDIF
       IF (ok_prescr_ust) THEN
       ust=ustar_prof_cas
       print *,'ust=',ust
       ENDIF

      endif !forcing_case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing from standard case :
!---------------------------------------------------------------------

      if (forcing_case2) then

         write(*,*),'avant call read2_1D_cas'
         call read2_1D_cas
         write(*,*) 'Forcing read'

!Time interpolation for initial conditions using interpolation routine
         write(*,*) 'AVT 1ere INTERPOLATION: day,day1 = ',daytime,day1   
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

      do l = 1, nlev_cas
      print *,'apres 1ere interp: plev_cas, plev_prof_cas=',l,plev_cas(l,1),plev_prof_cas(l)
      enddo

! vertical interpolation using interpolation routine:
!      write(*,*)'avant interp vert', t_prof
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

!       write(*,*) 'Profil initial forcing case interpole',t_mod

! initial and boundary conditions :
!      tsurf = ts_prof_cas
      ts_cur = ts_prof_cas 
      psurf=plev_prof_cas(1)
      write(*,*) 'SST initiale: ',tsurf
      do l = 1, llm
       temp(l) = t_mod_cas(l)
       q(l,1) = qv_mod_cas(l)
       q(l,2) = ql_mod_cas(l)
       u(l) = u_mod_cas(l)
       ug(l)= ug_mod_cas(l)
       v(l) = v_mod_cas(l)
       vg(l)= vg_mod_cas(l)
! Modif w_mod_cas -> omega_mod_cas (MM+MPL 20170309)
       omega(l) = omega_mod_cas(l)
       omega2(l)=omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!on applique le forcage total au premier pas de temps
!attention: signe different de toga
       d_t_adv(l) = alpha*omega(l)/rcpd+(ht_mod_cas(l)+vt_mod_cas(l))
       d_t_adv(l) = alpha*omega(l)/rcpd+(ht_mod_cas(l)+vt_mod_cas(l))
!      d_q_adv(l,1) = (hq_mod_cas(l)+vq_mod_cas(l))
       d_q_adv(l,1) = dq_mod_cas(l)
       d_q_adv(l,2) = 0.0
!      d_u_adv(l) = (hu_mod_cas(l)+vu_mod_cas(l))
       d_u_adv(l) = du_mod_cas(l)
!      d_v_adv(l) = (hv_mod_cas(l)+vv_mod_cas(l))
! correction bug d_u -> d_v (MM+MPL 20170310)
       d_v_adv(l) = dv_mod_cas(l)
      enddo     

! Faut-il multiplier par -1 ? (MPL 20160713)
       IF (ok_flux_surf) THEN
       fsens=-1.*sens_prof_cas
       flat=-1.*lat_prof_cas
       ENDIF
!
       IF (ok_prescr_ust) THEN
       ust=ustar_prof_cas
       print *,'ust=',ust
       ENDIF

      endif !forcing_case2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
! Forcing from standard case :
!---------------------------------------------------------------------

      if (forcing_SCM) then

         write(*,*),'avant call old_read_SCM_cas'
         call old_read_SCM_cas
         write(*,*) 'Forcing read' 

!Time interpolation for initial conditions using interpolation routine
         write(*,*) 'AVT 1ere INTERPOLATION: day,day1 = ',daytime,day1   
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

      do l = 1, nlev_cas
      print *,'apres 1ere interp: plev_cas, plev_prof_cas=',l,plev_cas(l,1),plev_prof_cas(l)
      enddo

! vertical interpolation using interpolation routine:
!      write(*,*)'avant interp vert', t_prof
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

!       write(*,*) 'Profil initial forcing case interpole',t_mod

! initial and boundary conditions :
!      tsurf = ts_prof_cas
      ts_cur = ts_prof_cas 
      psurf=plev_prof_cas(1)
      write(*,*) 'SST initiale: ',tsurf
      do l = 1, llm
       temp(l) = t_mod_cas(l)
       q(l,1) = qv_mod_cas(l)
       q(l,2) = ql_mod_cas(l)
       u(l) = u_mod_cas(l)
       ug(l)= ug_mod_cas(l)
       v(l) = v_mod_cas(l)
       vg(l)= vg_mod_cas(l)
! Modif w_mod_cas -> omega_mod_cas (MM+MPL 20170309)
       omega(l) = omega_mod_cas(l)
       omega2(l)=omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq

       alpha = rd*temp(l)*(1.+(rv/rd-1.)*q(l,1))/play(l)
!on applique le forcage total au premier pas de temps
!attention: signe different de toga
       d_t_adv(l) = alpha*omega(l)/rcpd+(ht_mod_cas(l)+vt_mod_cas(l))
       d_t_adv(l) = alpha*omega(l)/rcpd+(ht_mod_cas(l)+vt_mod_cas(l))
!      d_q_adv(l,1) = (hq_mod_cas(l)+vq_mod_cas(l))
       d_q_adv(l,1) = dq_mod_cas(l)
       d_q_adv(l,2) = 0.0
!      d_u_adv(l) = (hu_mod_cas(l)+vu_mod_cas(l))
       d_u_adv(l) = du_mod_cas(l)
!      d_v_adv(l) = (hv_mod_cas(l)+vv_mod_cas(l))
! correction bug d_u -> d_v (MM+MPL 20170310)
       d_v_adv(l) = dv_mod_cas(l)
      enddo     

! Faut-il multiplier par -1 ? (MPL 20160713)
       IF (ok_flux_surf) THEN
       fsens=-1.*sens_prof_cas
       flat=-1.*lat_prof_cas
       ENDIF
!
       IF (ok_prescr_ust) THEN
       ust=ustar_prof_cas
       print *,'ust=',ust
       ENDIF

      endif !forcing_SCM
