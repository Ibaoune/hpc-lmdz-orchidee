!
! $Header$
!-- Modified by : Filiberti M-A 06/2005
!
      real epmax             ! 0.993
      real coef_epmax_cape             ! 0.993
!jyg<
      REAL  cvl_comp_threshold     ! 0.
!>jyg
      logical ok_adj_ema      ! F
      integer iflag_clw      ! 0
      integer iflag_cvl_sigd
      real cvl_sig2feed      ! 0.97

!jyg<
!!      common/comconema1/epmax,coef_epmax_cape,ok_adj_ema,iflag_clw,sig1feed,sig2feed
!!      common/comconema2/iflag_cvl_sigd
      common/comconema1/epmax,coef_epmax_cape, cvl_comp_threshold, cvl_sig2feed
      common/comconema2/iflag_cvl_sigd, iflag_clw, ok_adj_ema
!>jyg

!      common/comconema/epmax,coef_epmax_cape,ok_adj_ema,iflag_clw
!$OMP THREADPRIVATE(/comconema1/)
!$OMP THREADPRIVATE(/comconema2/)

