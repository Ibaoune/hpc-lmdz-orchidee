
         print*,'FORCING CASE forcing_case2'
!       print*,                                                             &
!    & '#### ITAP,day,day1,(day-day1)*86400,(day-day1)*86400/pdt_cas=',     &
!    &    daytime,day1,(daytime-day1)*86400.,                               &
!    &    (daytime-day1)*86400/pdt_cas

! time interpolation:
        CALL interp_case_time_std(daytime,day1,annee_ref                                       &
!    &       ,year_ini_cas,day_ju_ini_cas,nt_cas,pdt_cas,nlev_cas                           &
     &       ,nt_cas,nlev_cas                                                               &
     &       ,ts_cas,ps_cas,plev_cas,t_cas,th_cas,thv_cas,thl_cas,qv_cas,ql_cas,qi_cas      &
     &       ,u_cas,v_cas,ug_cas,vg_cas                                                     &
     &       ,temp_nudg_cas,qv_nudg_cas,u_nudg_cas,v_nudg_cas                               &
     &       ,vitw_cas,omega_cas,du_cas,hu_cas,vu_cas                                       &
     &       ,dv_cas,hv_cas,vv_cas,dt_cas,ht_cas,vt_cas,dtrad_cas                           &
     &       ,dq_cas,hq_cas,vq_cas,dth_cas,hth_cas,vth_cas,lat_cas,sens_cas,ustar_cas       &
     &       ,uw_cas,vw_cas,q1_cas,q2_cas,tke_cas                                           &
!
     &       ,ts_prof_cas,ps_prof_cas,plev_prof_cas,t_prof_cas,theta_prof_cas,thv_prof_cas  &
     &       ,thl_prof_cas,qv_prof_cas,ql_prof_cas,qi_prof_cas                              &
     &       ,u_prof_cas,v_prof_cas,ug_prof_cas,vg_prof_cas                                 &
     &       ,temp_nudg_prof_cas,qv_nudg_prof_cas,u_nudg_prof_cas,v_nudg_prof_cas           &
     &       ,vitw_prof_cas,omega_prof_cas                                                  &
     &       ,du_prof_cas,hu_prof_cas,vu_prof_cas                                           &
     &       ,dv_prof_cas,hv_prof_cas,vv_prof_cas,dt_prof_cas,ht_prof_cas,vt_prof_cas       &
     &       ,dtrad_prof_cas,dq_prof_cas,hq_prof_cas,vq_prof_cas                            &
     &       ,dth_prof_cas,hth_prof_cas,vth_prof_cas,lat_prof_cas                           &
     &       ,sens_prof_cas,ustar_prof_cas,uw_prof_cas,vw_prof_cas,q1_prof_cas,q2_prof_cas,tke_prof_cas)

             ts_cur = ts_prof_cas 
!            psurf=plev_prof_cas(1)
             psurf=ps_prof_cas

! vertical interpolation:
      CALL interp2_case_vertical_std(play,nlev_cas,plev_prof_cas                                              &
     &         ,t_prof_cas,theta_prof_cas,thv_prof_cas,thl_prof_cas                                       &
     &         ,qv_prof_cas,ql_prof_cas,qi_prof_cas,u_prof_cas,v_prof_cas                                 &
     &         ,ug_prof_cas,vg_prof_cas                                                                   &
     &         ,temp_nudg_prof_cas,qv_nudg_prof_cas,u_nudg_prof_cas,v_nudg_prof_cas                       &
     &         ,vitw_prof_cas,omega_prof_cas                                      &
     &         ,du_prof_cas,hu_prof_cas,vu_prof_cas,dv_prof_cas,hv_prof_cas,vv_prof_cas                   &
     &         ,dt_prof_cas,ht_prof_cas,vt_prof_cas,dtrad_prof_cas,dq_prof_cas,hq_prof_cas,vq_prof_cas    &
     &         ,dth_prof_cas,hth_prof_cas,vth_prof_cas                                                    &
!
     &         ,t_mod_cas,theta_mod_cas,thv_mod_cas,thl_mod_cas,qv_mod_cas,ql_mod_cas,qi_mod_cas          &
     &         ,u_mod_cas,v_mod_cas,ug_mod_cas,vg_mod_cas                                                 &
     &         ,temp_nudg_mod_cas,qv_nudg_mod_cas,u_nudg_mod_cas,v_nudg_mod_cas                           &
     &         ,w_mod_cas,omega_mod_cas                                                                   &
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
      if (1==0) then
         DO l=2,llm-1
          d_t_z(l)=(temp(l+1)-temp(l-1))/(play(l+1)-play(l-1))
          d_th_z(l)=(teta(l+1)-teta(l-1))/(play(l+1)-play(l-1))
          d_q_z(l)=(q(l+1,1)-q(l-1,1))/(play(l+1)-play(l-1))
          d_u_z(l)=(u(l+1)-u(l-1))/(play(l+1)-play(l-1))
          d_v_z(l)=(v(l+1)-v(l-1))/(play(l+1)-play(l-1))
         ENDDO
      else
         DO l=2,llm-1
	    IF (omega(l)>0.) THEN
             d_t_z(l)=(temp(l+1)-temp(l))/(play(l+1)-play(l))
             d_th_z(l)=(teta(l+1)-teta(l))/(play(l+1)-play(l))
             d_q_z(l)=(q(l+1,1)-q(l,1))/(play(l+1)-play(l))
             d_u_z(l)=(u(l+1)-u(l))/(play(l+1)-play(l))
             d_v_z(l)=(v(l+1)-v(l))/(play(l+1)-play(l))
	    ELSE
             d_t_z(l)=(temp(l-1)-temp(l))/(play(l-1)-play(l))
             d_th_z(l)=(teta(l-1)-teta(l))/(play(l-1)-play(l))
             d_q_z(l)=(q(l-1,1)-q(l,1))/(play(l-1)-play(l))
             d_u_z(l)=(u(l-1)-u(l))/(play(l-1)-play(l))
             d_v_z(l)=(v(l-1)-v(l))/(play(l-1)-play(l))
	    ENDIF
	 ENDDO
      endif
      d_t_z(1)=d_t_z(2)
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

! TRAVAIL : PRENDRE DES NOTATIONS COHERENTES POUR W
      do l = 1, llm
! Modif w_mod_cas -> omega_mod_cas (MM+MPL 20170309)
       omega(l) = -w_mod_cas(l)*play(l)*rg/(rd*temp(l))
      enddo

!Calcul de l advection verticale
! Modif w_mod_cas -> omega_mod_cas (MM+MPL 20170310)
      d_t_dyn_z(:)=omega(:)*d_t_z(:)
      d_th_dyn_z(:)=omega(:)*d_th_z(:)
      d_q_dyn_z(:)=omega(:)*d_q_z(:)
      d_u_dyn_z(:)=omega(:)*d_u_z(:)
      d_v_dyn_z(:)=omega(:)*d_v_z(:)

!geostrophic wind 
      if (forc_geo.eq.1) then
        do l=1,llm
        ug(l) = ug_mod_cas(l) 
        vg(l) = vg_mod_cas(l) 
        enddo
      endif

      do l = 1, llm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modif w_mod_cas -> omega_mod_cas (MM+MPL 20170309)
       !!! omega2(l)= omega(l)/rg*airefi ! flxmass_w calcule comme ds physiq
       omega(l) = omega_mod_cas(l)
       omega2(l)= omega_mod_cas(l)/rg*airefi ! flxmass_w calcule comme ds physiq

! On effectue la somme du forcage total et de la decomposition
! horizontal/vertical en supposant que soit l'un soit l'autre
! sont remplis mais jamais les deux

       d_t_adv(l) = dt_mod_cas(l)+ht_mod_cas(l)+vt_mod_cas(l)
       d_q_adv(l,1) = dq_mod_cas(l)+hq_mod_cas(l)+vq_mod_cas(l)
       d_q_adv(l,2) = 0.0
       d_u_adv(l) = du_mod_cas(l)+hu_mod_cas(l)+vu_mod_cas(l)
       d_v_adv(l) = dv_mod_cas(l)+hv_mod_cas(l)+vv_mod_cas(l)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CONSERVE EN ATTENDANT QUE LE CAS EN QUESTION FONCTIONNE EN STD !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (forc_w==1) then
!          d_q_adv(l,1)=d_q_adv(l,1)-d_q_dyn_z(l)
!          d_t_adv(l)=d_t_adv(l)-d_t_dyn_z(l)
!          d_v_adv(l)=d_v_adv(l)-d_v_dyn_z(l)
!          d_u_adv(l)=d_u_adv(l)-d_u_dyn_z(l)
!       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
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
