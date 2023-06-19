      integer            :: iflag_thermals,nsplit_thermals

!!! nrlmd le 10/04/2012
      integer            :: iflag_trig_bl,iflag_clos_bl
      integer            :: tau_trig_shallow,tau_trig_deep
      real               :: s_trig
!!! fin nrlmd le 10/04/2012

      real,parameter     :: r_aspect_thermals=2.,l_mix_thermals=30.
      real               :: alp_bl_k
      real               :: tau_thermals,fact_thermals_ed_dz
      integer,parameter  :: w2di_thermals=0
      integer            :: isplit

      integer            :: iflag_coupl,iflag_clos,iflag_wake
      integer            :: iflag_thermals_ed,iflag_thermals_optflux,iflag_thermals_closure

      common/ctherm1/iflag_thermals,nsplit_thermals,iflag_thermals_closure
      common/ctherm2/tau_thermals,alp_bl_k,fact_thermals_ed_dz
      common/ctherm4/iflag_coupl,iflag_clos,iflag_wake
      common/ctherm5/iflag_thermals_ed,iflag_thermals_optflux

!!! nrlmd le 10/04/2012
      common/ctherm6/iflag_trig_bl,iflag_clos_bl
      common/ctherm7/tau_trig_shallow,tau_trig_deep
      common/ctherm8/s_trig
!!! fin nrlmd le 10/04/2012

!$OMP THREADPRIVATE(/ctherm1/,/ctherm2/,/ctherm4/,/ctherm5/)
!$OMP THREADPRIVATE(/ctherm6/,/ctherm7/,/ctherm8/)
