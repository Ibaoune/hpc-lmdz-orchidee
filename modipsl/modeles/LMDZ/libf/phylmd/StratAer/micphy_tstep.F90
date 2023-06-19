!
! $Id: micphy_tstep.F90 3677 2020-05-06 15:18:32Z oboucher $
!
SUBROUTINE micphy_tstep(pdtphys,tr_seri,t_seri,pplay,paprs,rh,is_strato)

  USE geometry_mod, ONLY : latitude_deg !NL- latitude corr. to local domain
  USE dimphy, ONLY : klon,klev
  USE aerophys
  USE infotrac_phy
  USE phys_local_var_mod, ONLY: mdw, budg_3D_nucl, budg_3D_cond_evap, budg_h2so4_to_part, R2SO4, DENSO4, f_r_wet
  USE nucleation_tstep_mod
  USE cond_evap_tstep_mod
  USE sulfate_aer_mod, ONLY : STRAACT
  USE YOMCST, ONLY : RPI, RD, RG
  USE print_control_mod, ONLY: lunout
  USE strataer_mod
  
  IMPLICIT NONE

  !--------------------------------------------------------

  ! transfer variables when calling this routine
  REAL,INTENT(IN)                               :: pdtphys ! Pas d'integration pour la physique (seconde)
  REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT)  :: tr_seri ! Concentration Traceur [U/KgA]
  REAL,DIMENSION(klon,klev),INTENT(IN)          :: t_seri  ! Temperature
  REAL,DIMENSION(klon,klev),INTENT(IN)          :: pplay   ! pression pour le mileu de chaque couche (en Pa)
  REAL,DIMENSION(klon,klev+1),INTENT(IN)        :: paprs   ! pression pour chaque inter-couche (en Pa)
  REAL,DIMENSION(klon,klev),INTENT(IN)          :: rh      ! humidite relative
  LOGICAL,DIMENSION(klon,klev),INTENT(IN)       :: is_strato

  ! local variables in coagulation routine
  INTEGER, PARAMETER        :: nbtstep=4  ! Max number of time steps in microphysics per time step in physics
  INTEGER                   :: it,ilon,ilev,count_tstep
  REAL                      :: rhoa !H2SO4 number density [molecules/cm3]
  REAL                      :: ntot !total number of molecules in the critical cluster (ntot>4)
  REAL                      :: x    ! molefraction of H2SO4 in the critical cluster     
  REAL Vbin(nbtr_bin)
  REAL a_xm, b_xm, c_xm
  REAL PDT, dt
  REAL H2SO4_init
  REAL ACTSO4(klon,klev)
  REAL RRSI(nbtr_bin)
  REAL nucl_rate
  REAL cond_evap_rate
  REAL evap_rate
  REAL FL(nbtr_bin)
  REAL ASO4(nbtr_bin)
  REAL DNDR(nbtr_bin)
  REAL H2SO4_sat(nbtr_bin)

  DO it=1,nbtr_bin
    Vbin(it)=4.0*RPI*((mdw(it)/2.)**3)/3.0
  ENDDO

  !coefficients for H2SO4 density parametrization used for nucleation if ntot<4
  a_xm = 0.7681724 + 1.*(2.1847140 + 1.*(7.1630022 + 1.*(-44.31447 + &
       & 1.*(88.75606 + 1.*(-75.73729 + 1.*23.43228)))))
  b_xm = 1.808225e-3 + 1.*(-9.294656e-3 + 1.*(-0.03742148 + 1.*(0.2565321 + &
       & 1.*(-0.5362872 + 1.*(0.4857736 - 1.*0.1629592)))))
  c_xm = -3.478524e-6 + 1.*(1.335867e-5 + 1.*(5.195706e-5 + 1.*(-3.717636e-4 + &
       & 1.*(7.990811e-4 + 1.*(-7.458060e-4 + 1.*2.58139e-4 )))))

  ! STRAACT (R2SO4, t_seri -> H2SO4 activity coefficient (ACTSO4)) for cond/evap
  CALL STRAACT(ACTSO4)

  ! compute particle radius in cm RRSI from diameter in m 
  DO it=1,nbtr_bin
    RRSI(it)=mdw(it)/2.*100.
  ENDDO

  DO ilon=1, klon
!
!--initialisation of diagnostic
  budg_h2so4_to_part(ilon)=0.0
!
  DO ilev=1, klev
!
!--initialisation of diagnostic
  budg_3D_nucl(ilon,ilev)=0.0
  budg_3D_cond_evap(ilon,ilev)=0.0
!
  ! only in the stratosphere
  IF (is_strato(ilon,ilev)) THEN
    ! initialize sulfur fluxes
    H2SO4_init=tr_seri(ilon,ilev,id_H2SO4_strat)
    ! adaptive timestep for nucleation and condensation
    PDT=pdtphys
    count_tstep=0
    DO WHILE (PDT>0.0)
      count_tstep=count_tstep+1
      IF (count_tstep .GT. nbtstep)  EXIT
      ! convert tr_seri(GASH2SO4) (in kg/kgA) to H2SO4 number density (in molecules/cm3) 
      rhoa=tr_seri(ilon,ilev,id_H2SO4_strat) &
          & *pplay(ilon,ilev)/t_seri(ilon,ilev)/RD/1.E6/mH2SO4mol
      ! compute nucleation rate in kg(H2SO4)/kgA/s
      CALL nucleation_rate(rhoa,t_seri(ilon,ilev),pplay(ilon,ilev),rh(ilon,ilev), &
           & a_xm,b_xm,c_xm,nucl_rate,ntot,x)
      !NL - add nucleation box (if flag on)
      IF (flag_nuc_rate_box) THEN
         IF (latitude_deg(ilon).LE.nuclat_min .OR. latitude_deg(ilon).GE.nuclat_max &
              .OR. pplay(ilon,ilev).GE.nucpres_max .AND. pplay(ilon,ilev).LE.nucpres_min) THEN
            nucl_rate=0.0
         ENDIF
      ENDIF
      ! compute cond/evap rate in kg(H2SO4)/kgA/s
      CALL condens_evapor_rate(rhoa,t_seri(ilon,ilev),pplay(ilon,ilev), &
             & ACTSO4(ilon,ilev),R2SO4(ilon,ilev),DENSO4(ilon,ilev),f_r_wet(ilon,ilev), &
             & RRSI,Vbin,FL,ASO4,DNDR)
      ! consider only condensation (positive FL)
      DO it=1,nbtr_bin
        FL(it)=MAX(FL(it),0.)
      ENDDO
      ! compute total H2SO4 cond flux for all particles
      cond_evap_rate=0.0
      DO it=1, nbtr_bin
        cond_evap_rate=cond_evap_rate+tr_seri(ilon,ilev,it+nbtr_sulgas)*FL(it)*mH2SO4mol
      ENDDO
      ! determine appropriate time step
      dt=(H2SO4_init-H2SO4_sat(nbtr_bin))/float(nbtstep)/MAX(1.e-30, nucl_rate+cond_evap_rate) !cond_evap_rate pos. for cond. and neg. for evap.
      IF (dt.LT.0.0) THEN
        dt=PDT
      ENDIF
      dt=MIN(dt,PDT)
      ! update H2SO4 concentration
      tr_seri(ilon,ilev,id_H2SO4_strat)=MAX(0.,tr_seri(ilon,ilev,id_H2SO4_strat)-(nucl_rate+cond_evap_rate)*dt)
      ! apply cond to bins
      CALL cond_evap_part(dt,FL,ASO4,f_r_wet(ilon,ilev),RRSI,Vbin,tr_seri(ilon,ilev,:))
      ! apply nucl. to bins
      CALL nucleation_part(nucl_rate,ntot,x,dt,Vbin,tr_seri(ilon,ilev,:))
      ! compute fluxes as diagnostic in [kg(S)/m2/layer/s] (now - for evap and + for cond)
      budg_3D_cond_evap(ilon,ilev)=budg_3D_cond_evap(ilon,ilev)+mSatom/mH2SO4mol &
               & *cond_evap_rate*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG*dt/pdtphys
      budg_3D_nucl(ilon,ilev)=budg_3D_nucl(ilon,ilev)+mSatom/mH2SO4mol &
               & *nucl_rate*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG*dt/pdtphys
      ! update time step
      PDT=PDT-dt
    ENDDO
    ! convert tr_seri(GASH2SO4) (in kg/kgA) to H2SO4 number density (in molecules/cm3) 
    rhoa=tr_seri(ilon,ilev,id_H2SO4_strat) &
        & *pplay(ilon,ilev)/t_seri(ilon,ilev)/RD/1.E6/mH2SO4mol
    ! compute cond/evap rate in kg(H2SO4)/kgA/s (now only evap for pdtphys)
    CALL condens_evapor_rate(rhoa,t_seri(ilon,ilev),pplay(ilon,ilev), &
           & ACTSO4(ilon,ilev),R2SO4(ilon,ilev),DENSO4(ilon,ilev),f_r_wet(ilon,ilev), &
           & RRSI,Vbin,FL,ASO4,DNDR)
    ! limit evaporation (negative FL) over one physics time step to H2SO4 content of the droplet
    DO it=1,nbtr_bin
      FL(it)=MAX(FL(it)*pdtphys,0.-ASO4(it))/pdtphys
      ! consider only evap (negative FL)
      FL(it)=MIN(FL(it),0.)
    ENDDO
    ! compute total H2SO4 evap flux for all particles
    evap_rate=0.0
    DO it=1, nbtr_bin
      evap_rate=evap_rate+tr_seri(ilon,ilev,it+nbtr_sulgas)*FL(it)*mH2SO4mol
    ENDDO
    ! update H2SO4 concentration after evap
    tr_seri(ilon,ilev,id_H2SO4_strat)=MAX(0.,tr_seri(ilon,ilev,id_H2SO4_strat)-evap_rate*pdtphys)
    ! apply evap to bins
    CALL cond_evap_part(pdtphys,FL,ASO4,f_r_wet(ilon,ilev),RRSI,Vbin,tr_seri(ilon,ilev,:))
    ! compute fluxes as diagnostic in [kg(S)/m2/layer/s] (now - for evap and + for cond)
    budg_3D_cond_evap(ilon,ilev)=budg_3D_cond_evap(ilon,ilev)+mSatom/mH2SO4mol &
             & *evap_rate*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG
    ! compute vertically integrated flux due to the net effect of nucleation and condensation/evaporation
    budg_h2so4_to_part(ilon)=budg_h2so4_to_part(ilon)+(H2SO4_init-tr_seri(ilon,ilev,id_H2SO4_strat)) &
             & *mSatom/mH2SO4mol*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG/pdtphys
  ENDIF
  ENDDO
  ENDDO

  IF (MINVAL(tr_seri).LT.0.0) THEN
    DO ilon=1, klon
    DO ilev=1, klev    
    DO it=1, nbtr
      IF (tr_seri(ilon,ilev,it).LT.0.0) THEN
         WRITE(lunout,*) 'micphy_tstep: negative concentration', tr_seri(ilon,ilev,it), ilon, ilev, it
      ENDIF
    ENDDO
    ENDDO
    ENDDO
  ENDIF

END SUBROUTINE micphy_tstep
