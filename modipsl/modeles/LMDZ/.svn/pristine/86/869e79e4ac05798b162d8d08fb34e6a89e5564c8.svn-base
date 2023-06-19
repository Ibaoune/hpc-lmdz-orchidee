!
! $Id: compar1d.h 2010-08-04 17:02:56Z lahellec $
!
      integer :: forcing_type
      integer :: tend_u,tend_v,tend_w,tend_t,tend_q,tend_rayo
      real :: nudge_u,nudge_v,nudge_w,nudge_t,nudge_q
      integer :: iflag_nudge
      real :: nat_surf
      real :: tsurf
      real :: rugos
      real :: rugosh
      real :: xqsol(1:2)
      real :: qsurf
      real :: psurf
      real :: zsurf
      real :: albedo
      real :: snowmass

      real :: time
      real :: time_ini
      real :: xlat 
      real :: xlon 
      real :: airefi
      real :: wtsurf 
      real :: wqsurf 
      real :: restart_runoff 
      real :: xagesno 
      real :: qsolinp 
      real :: zpicinp 

      logical :: restart
      logical :: ok_old_disvert

! Pour les forcages communs: ces entiers valent 0 ou 1
! tadv= advection tempe, tadvv= adv tempe verticale, tadvh= adv tempe horizontale
! idem pour l advection en theta
! qadv= advection q, qadvv= adv q verticale, qadvh= adv q horizontale
! trad= 0 (rayonnement actif) ou 1 (prescrit par tend_rad) ou adv (prescir et contenu dans les tadv)
! forcages en omega, w, vent geostrophique ou ustar
! Parametres de nudging en u,v,t,q valent 0 ou 1 ou le temps de nudging

      integer :: tadv, tadvv, tadvh, qadv, qadvv, qadvh, thadv, thadvv, thadvh, trad
      integer :: forc_omega, forc_u, forc_v, forc_w, forc_geo, forc_ustar
      real    :: nudging_u, nudging_v, nudging_w, nudging_t, nudging_qv
      real    :: p_nudging_u, p_nudging_v, p_nudging_w, p_nudging_t, p_nudging_qv
      common/com_par1d/                                                 &
     & nat_surf,tsurf,rugos,rugosh,                                     &
     & xqsol,qsurf,psurf,zsurf,albedo,time,time_ini,xlat,xlon,airefi,   &
     & wtsurf,wqsurf,restart_runoff,xagesno,qsolinp,zpicinp,            &
     & forcing_type,tend_u,tend_v,tend_w,tend_t,tend_q,tend_rayo,       &
     & nudge_u,nudge_v,nudge_w,nudge_t,nudge_q,                         &
     & iflag_nudge,snowmass,                                            &
     & restart,ok_old_disvert,                                          &
     & tadv, tadvv, tadvh, qadv, qadvv, qadvh, thadv, thadvv, thadvh,   &
     & trad, forc_omega, forc_u, forc_v, forc_w, forc_geo, forc_ustar,  &
     & nudging_u, nudging_v, nudging_t, nudging_qv, nudging_w,          &
     & p_nudging_u, p_nudging_v, p_nudging_t, p_nudging_qv, p_nudging_w

!$OMP THREADPRIVATE(/com_par1d/)











