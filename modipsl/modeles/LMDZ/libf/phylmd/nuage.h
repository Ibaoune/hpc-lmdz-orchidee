!
! $Id: nuage.h 2945 2017-07-12 14:20:24Z jbmadeleine $
!
      REAL rad_froid, rad_chau1, rad_chau2, t_glace_max, t_glace_min
      REAL exposant_glace
      REAL rei_min,rei_max
      REAL tau_cld_cv,coefw_cld_cv

      REAL tmax_fonte_cv

      INTEGER iflag_t_glace, iflag_cloudth_vert, iflag_cld_cv
      INTEGER iflag_rain_incloud_vol

      common /nuagecom/ rad_froid,rad_chau1, rad_chau2,t_glace_max,     &
     &                  t_glace_min,exposant_glace,rei_min,rei_max,     &
     &                  tau_cld_cv,coefw_cld_cv,                        &
     &                  tmax_fonte_cv,                                  &
     &                  iflag_t_glace,iflag_cloudth_vert,iflag_cld_cv,  &
     &                  iflag_rain_incloud_vol
!$OMP THREADPRIVATE(/nuagecom/)
