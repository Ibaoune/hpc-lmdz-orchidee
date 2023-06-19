!
! $Header$
!
      logical cvflag_grav
      logical cvflag_ice
      logical ok_optim_yield
      logical ok_entrain
      logical ok_homo_tend
      logical ok_convstop
      logical ok_intermittent
      logical cvflag_prec_eject
      logical qsat_depends_on_qt
      logical adiab_ascent_mass_flux_depends_on_ejectliq
      logical keepbug_ice_frac
      integer icvflag_Tpa

      COMMON /cvflag/ icvflag_Tpa, &
                      cvflag_grav, cvflag_ice, &
                      ok_optim_yield, &
                      ok_entrain, &
                      ok_homo_tend, &
                      ok_convstop, ok_intermittent, &
                      cvflag_prec_eject, &
                      qsat_depends_on_qt, &
                      adiab_ascent_mass_flux_depends_on_ejectliq, &
                      keepbug_ice_frac
!$OMP THREADPRIVATE(/cvflag/)
