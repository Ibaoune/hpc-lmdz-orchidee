!
! $Id: traccoag_mod.F90 3677 2020-05-06 15:18:32Z oboucher $
!
MODULE traccoag_mod
!
! This module calculates the concentration of aerosol particles in certain size bins
! considering coagulation and sedimentation.
!
CONTAINS

  SUBROUTINE traccoag(pdtphys, gmtime, debutphy, julien, &
       presnivs, xlat, xlon, pphis, pphi, &
       t_seri, pplay, paprs, sh, rh, tr_seri)

    USE phys_local_var_mod, ONLY: mdw, R2SO4, DENSO4, f_r_wet, surf_PM25_sulf, & 
        & budg_emi_ocs, budg_emi_so2, budg_emi_h2so4, budg_emi_part 

    USE dimphy
    USE infotrac_phy
    USE aerophys
    USE geometry_mod, ONLY : cell_area, boundslat
    USE mod_grid_phy_lmdz
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    USE mod_phys_lmdz_para, only: gather, scatter
    USE phys_cal_mod
    USE sulfate_aer_mod
    USE phys_local_var_mod, ONLY: stratomask
    USE YOMCST
    USE print_control_mod, ONLY: lunout
    USE strataer_mod
    USE phys_cal_mod, ONLY : year_len

    IMPLICIT NONE

! Input argument
!---------------
    REAL,INTENT(IN)    :: pdtphys    ! Pas d'integration pour la physique (seconde)
    REAL,INTENT(IN)    :: gmtime     ! Heure courante
    LOGICAL,INTENT(IN) :: debutphy   ! le flag de l'initialisation de la physique
    INTEGER,INTENT(IN) :: julien     ! Jour julien

    REAL,DIMENSION(klev),INTENT(IN)        :: presnivs! pressions approximat. des milieux couches (en PA)
    REAL,DIMENSION(klon),INTENT(IN)        :: xlat    ! latitudes pour chaque point 
    REAL,DIMENSION(klon),INTENT(IN)        :: xlon    ! longitudes pour chaque point 
    REAL,DIMENSION(klon),INTENT(IN)        :: pphis   ! geopotentiel du sol
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pphi    ! geopotentiel de chaque couche

    REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
    REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs   ! pression pour chaque inter-couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: sh      ! humidite specifique
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: rh      ! humidite relative   

! Output argument
!----------------
    REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT)  :: tr_seri ! Concentration Traceur [U/KgA]  

! Local variables
!----------------
    REAL                                   :: m_aer_emiss_vol_daily ! daily injection mass emission
    INTEGER                                :: it, k, i, ilon, ilev, itime, i_int, ieru
    LOGICAL,DIMENSION(klon,klev)           :: is_strato           ! true = above tropopause, false = below
    REAL,DIMENSION(klon,klev)              :: m_air_gridbox       ! mass of air in every grid box [kg]
    REAL,DIMENSION(klon_glo,klev,nbtr)     :: tr_seri_glo         ! Concentration Traceur [U/KgA]  
    REAL,DIMENSION(klev+1)                 :: altLMDz             ! altitude of layer interfaces in m
    REAL,DIMENSION(klev)                   :: f_lay_emiss         ! fraction of emission for every vertical layer
    REAL                                   :: f_lay_sum           ! sum of layer emission fractions
    REAL                                   :: alt                 ! altitude for integral calculation
    INTEGER,PARAMETER                      :: n_int_alt=10        ! number of subintervals for integration over Gaussian emission profile
    REAL,DIMENSION(nbtr_bin)               :: r_bin               ! particle radius in size bin [m]
    REAL,DIMENSION(nbtr_bin)               :: r_lower             ! particle radius at lower bin boundary [m]
    REAL,DIMENSION(nbtr_bin)               :: r_upper             ! particle radius at upper bin boundary [m]
    REAL,DIMENSION(nbtr_bin)               :: m_part_dry          ! mass of one dry particle in size bin [kg]
    REAL                                   :: zrho                ! Density of air [kg/m3]
    REAL                                   :: zdz                 ! thickness of atm. model layer in m
    REAL,DIMENSION(klev)                   :: zdm                 ! mass of atm. model layer in kg
    REAL,DIMENSION(klon,klev)              :: dens_aer            ! density of aerosol particles [kg/m3 aerosol] with default H2SO4 mass fraction
    REAL                                   :: emission            ! emission
    REAL                                   :: theta_min, theta_max ! for SAI computation between two latitudes
    REAL                                   :: dlat_loc

    IF (is_mpi_root) THEN
       WRITE(lunout,*) 'in traccoag: date from phys_cal_mod =',year_cur,'-',mth_cur,'-',day_cur,'-',hour
       WRITE(lunout,*) 'IN traccoag flag_sulf_emit: ',flag_sulf_emit
    ENDIF
    
    DO it=1, nbtr_bin
      r_bin(it)=mdw(it)/2.
    ENDDO

!--set boundaries of size bins
    DO it=1, nbtr_bin
    IF (it.EQ.1) THEN
      r_upper(it)=sqrt(r_bin(it+1)*r_bin(it))
      r_lower(it)=r_bin(it)**2./r_upper(it)
    ELSEIF (it.EQ.nbtr_bin) THEN
      r_lower(it)=sqrt(r_bin(it)*r_bin(it-1))
      r_upper(it)=r_bin(it)**2./r_lower(it)
    ELSE
      r_lower(it)=sqrt(r_bin(it)*r_bin(it-1))
      r_upper(it)=sqrt(r_bin(it+1)*r_bin(it))
    ENDIF
    ENDDO

    IF (debutphy .and. is_mpi_root) THEN
      DO it=1, nbtr_bin
        WRITE(lunout,*) 'radius bin', it, ':', r_bin(it), '(from',  r_lower(it), 'to', r_upper(it), ')'
      ENDDO
    ENDIF

!--initialising logical is_strato from stratomask
    is_strato(:,:)=.FALSE. 
    WHERE (stratomask.GT.0.5) is_strato=.TRUE. 

! STRACOMP (H2O, P, t_seri -> aerosol composition (R2SO4)) 
! H2SO4 mass fraction in aerosol (%)
    CALL stracomp(sh,t_seri,pplay)

! aerosol density (gr/cm3)
    CALL denh2sa(t_seri)

! compute factor for converting dry to wet radius (for every grid box)
    f_r_wet(:,:) = (dens_aer_dry/(DENSO4(:,:)*1000.)/(R2SO4(:,:)/100.))**(1./3.)

!--calculate mass of air in every grid box
    DO ilon=1, klon
    DO ilev=1, klev
      m_air_gridbox(ilon,ilev)=(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG*cell_area(ilon)
    ENDDO
    ENDDO

!    IF (debutphy) THEN
!      CALL gather(tr_seri, tr_seri_glo)
!      IF (MAXVAL(tr_seri_glo).LT.1.e-30) THEN
!--initialising tracer concentrations to zero
!        DO it=1, nbtr
!        tr_seri(:,:,it)=0.0
!        ENDDO
!      ENDIF
!    ENDIF

!--initialise emission diagnostics
    budg_emi_ocs(:)=0.0
    budg_emi_so2(:)=0.0
    budg_emi_h2so4(:)=0.0
    budg_emi_part(:)=0.0

!--sulfur emission, depending on chosen scenario (flag_sulf_emit)
    SELECT CASE(flag_sulf_emit)

    CASE(0) ! background aerosol
      ! do nothing (no emission)

    CASE(1) ! volcanic eruption
      !--only emit on day of eruption
      ! stretch emission over one day of Pinatubo eruption
       DO ieru=1, nErupt
          IF (year_cur==year_emit_vol(ieru).AND.mth_cur==mth_emit_vol(ieru).AND.&
               day_cur>=day_emit_vol(ieru).AND.day_cur<(day_emit_vol(ieru)+injdur)) THEN
             !
             ! daily injection mass emission - NL
             m_aer_emiss_vol_daily = m_aer_emiss_vol(ieru)/(REAL(injdur)*REAL(ponde_lonlat_vol(ieru)))
             WRITE(lunout,*) 'IN traccoag DD m_aer_emiss_vol(ieru)=',m_aer_emiss_vol(ieru), &
                  'ponde_lonlat_vol(ieru)=',ponde_lonlat_vol(ieru),'(injdur*ponde_lonlat_vol(ieru))', &
                  (injdur*ponde_lonlat_vol(ieru)),'m_aer_emiss_vol_daily=',m_aer_emiss_vol_daily,'ieru=',ieru
             WRITE(lunout,*) 'IN traccoag, dlon=',dlon
             DO i=1,klon
                !Pinatubo eruption at 15.14N, 120.35E
                dlat_loc=180./RPI/2.*(boundslat(i,1)-boundslat(i,3)) ! dlat = half difference of boundary latitudes
                WRITE(lunout,*) 'IN traccoag, dlat=',dlat_loc
                IF ( xlat(i).GE.xlat_min_vol(ieru)-dlat_loc .AND. xlat(i).LT.xlat_max_vol(ieru)+dlat_loc .AND. &
                     xlon(i).GE.xlon_min_vol(ieru)-dlon .AND. xlon(i).LT.xlon_max_vol(ieru)+dlon ) THEN
                   !
                   WRITE(lunout,*) 'coordinates of volcanic injection point=',xlat(i),xlon(i),day_cur,mth_cur,year_cur
                   WRITE(lunout,*) 'DD m_aer_emiss_vol_daily=',m_aer_emiss_vol_daily
                   !         compute altLMDz
                   altLMDz(:)=0.0
                   DO k=1, klev
                      zrho=pplay(i,k)/t_seri(i,k)/RD       !air density in kg/m3
                      zdm(k)=(paprs(i,k)-paprs(i,k+1))/RG  !mass of layer in kg
                      zdz=zdm(k)/zrho                      !thickness of layer in m
                      altLMDz(k+1)=altLMDz(k)+zdz          !altitude of interface
                   ENDDO

                   SELECT CASE(flag_sulf_emit_distrib)
                   
                   CASE(0) ! Gaussian distribution
                   !compute distribution of emission to vertical model layers (based on Gaussian peak in altitude)
                   f_lay_sum=0.0
                   DO k=1, klev
                      f_lay_emiss(k)=0.0
                      DO i_int=1, n_int_alt
                         alt=altLMDz(k)+float(i_int)*(altLMDz(k+1)-altLMDz(k))/float(n_int_alt)
                         f_lay_emiss(k)=f_lay_emiss(k)+1./(sqrt(2.*RPI)*sigma_alt_vol(ieru))* &
                              &              exp(-0.5*((alt-altemiss_vol(ieru))/sigma_alt_vol(ieru))**2.)*   &
                              &              (altLMDz(k+1)-altLMDz(k))/float(n_int_alt)
                      ENDDO
                      f_lay_sum=f_lay_sum+f_lay_emiss(k)
                   ENDDO
                   
                   CASE(1) ! Uniform distribution
                   ! In this case, parameter sigma_alt_vol(ieru) is considered to be half the
                   ! height of the injection, centered around altemiss_vol(ieru)
                   DO k=1, klev
                      f_lay_emiss(k)=max(min(altemiss_vol(ieru)+sigma_alt_vol(ieru),altLMDz(k+1))- &
                      & max(altemiss_vol(ieru)-sigma_alt_vol(ieru),altLMDz(k)),0.)/(2.*sigma_alt_vol(ieru))
                      f_lay_sum=f_lay_sum+f_lay_emiss(k)
                   ENDDO

                   END SELECT        ! End CASE over flag_sulf_emit_distrib)

                   WRITE(lunout,*) "IN traccoag m_aer_emiss_vol=",m_aer_emiss_vol(ieru)
                   WRITE(lunout,*) "IN traccoag f_lay_emiss=",f_lay_emiss
                   !correct for step integration error
                   f_lay_emiss(:)=f_lay_emiss(:)/f_lay_sum
                   !emission as SO2 gas (with m(SO2)=64/32*m_aer_emiss)
                   !vertically distributed emission
                   DO k=1, klev
                      ! stretch emission over one day of Pinatubo eruption
                      emission=m_aer_emiss_vol_daily*(mSO2mol/mSatom)/m_air_gridbox(i,k)*f_lay_emiss(k)/1./(86400.-pdtphys)
                      tr_seri(i,k,id_SO2_strat)=tr_seri(i,k,id_SO2_strat)+emission*pdtphys
                      budg_emi_so2(i)=budg_emi_so2(i)+emission*zdm(k)*mSatom/mSO2mol
                   ENDDO
                ENDIF ! emission grid cell
             ENDDO ! klon loop
             WRITE(lunout,*) "IN traccoag (ieru=",ieru,") m_aer_emiss_vol_daily=",m_aer_emiss_vol_daily
          ENDIF ! emission period
       ENDDO ! eruption number
       
    CASE(2) ! stratospheric aerosol injections (SAI)
!
      DO i=1,klon
!       SAI standard scenario with continuous emission from 1 grid point at the equator
!       SAI emission on single month
!       SAI continuous emission o
        dlat_loc=180./RPI/2.*(boundslat(i,1)-boundslat(i,3)) ! dlat = half difference of boundary latitudes
        IF  ( xlat(i).GE.xlat_sai-dlat_loc .AND. xlat(i).LT.xlat_sai+dlat_loc .AND. &
          &   xlon(i).GE.xlon_sai-dlon .AND. xlon(i).LT.xlon_sai+dlon ) THEN
!
!         compute altLMDz
          altLMDz(:)=0.0
          DO k=1, klev
            zrho=pplay(i,k)/t_seri(i,k)/RD       !air density in kg/m3
            zdm(k)=(paprs(i,k)-paprs(i,k+1))/RG  !mass of layer in kg
            zdz=zdm(k)/zrho                      !thickness of layer in m
            altLMDz(k+1)=altLMDz(k)+zdz          !altitude of interface
          ENDDO

          SELECT CASE(flag_sulf_emit_distrib)

          CASE(0) ! Gaussian distribution
          !compute distribution of emission to vertical model layers (based on Gaussian peak in altitude)
          f_lay_sum=0.0
               DO k=1, klev
                     f_lay_emiss(k)=0.0
                     DO i_int=1, n_int_alt
                         alt=altLMDz(k)+float(i_int)*(altLMDz(k+1)-altLMDz(k))/float(n_int_alt)
                         f_lay_emiss(k)=f_lay_emiss(k)+1./(sqrt(2.*RPI)*sigma_alt_sai)* &
                         &              exp(-0.5*((alt-altemiss_sai)/sigma_alt_sai)**2.)*   & 
                         &              (altLMDz(k+1)-altLMDz(k))/float(n_int_alt) 
                     ENDDO
                     f_lay_sum=f_lay_sum+f_lay_emiss(k)
               ENDDO

          CASE(1) ! Uniform distribution
          f_lay_sum=0.0
          ! In this case, parameter sigma_alt_vol(ieru) is considered to be half
          ! the height of the injection, centered around altemiss_sai
               DO k=1, klev
                    f_lay_emiss(k)=max(min(altemiss_sai+sigma_alt_sai,altLMDz(k+1))- &
                    & max(altemiss_sai-sigma_alt_sai,altLMDz(k)),0.)/(2.*sigma_alt_sai)
                    f_lay_sum=f_lay_sum+f_lay_emiss(k)
               ENDDO

          END SELECT ! Gaussian or uniform distribution

          !correct for step integration error
          f_lay_emiss(:)=f_lay_emiss(:)/f_lay_sum
          !emission as SO2 gas (with m(SO2)=64/32*m_aer_emiss)
          !vertically distributed emission
          DO k=1, klev
            ! stretch emission over whole year (360d)
            emission=m_aer_emiss_sai*(mSO2mol/mSatom)/m_air_gridbox(i,k)*f_lay_emiss(k)/FLOAT(year_len)/86400.  
            tr_seri(i,k,id_SO2_strat)=tr_seri(i,k,id_SO2_strat)+emission*pdtphys
            budg_emi_so2(i)=budg_emi_so2(i)+emission*zdm(k)*mSatom/mSO2mol
          ENDDO

!          !emission as monodisperse particles with 0.1um dry radius (BIN21)
!          !vertically distributed emission
!          DO k=1, klev
!            ! stretch emission over whole year (360d)
!            emission=m_aer_emiss*(mH2SO4mol/mSatom)/m_part_dry(21)/m_air_gridbox(i,k)*f_lay_emiss(k)/FLOAT(year_len)/86400.
!            tr_seri(i,k,id_BIN01_strat+20)=tr_seri(i,k,id_BIN01_strat+20)+emission*pdtphys 
!            budg_emi_part(i)=budg_emi_part(i)+emission*zdm(k)*mSatom/mH2SO4mol
!          ENDDO
        ENDIF ! emission grid cell
      ENDDO ! klon loop

    CASE(3) ! --- SAI injection over a single band of longitude and between
            !     lat_min and lat_max

    DO i=1,klon
!       SAI scenario with continuous emission
        dlat_loc=180./RPI/2.*(boundslat(i,1)-boundslat(i,3)) ! dlat = half difference of boundary latitudes
        theta_min = max(xlat(i)-dlat_loc,xlat_min_sai)
        theta_max = min(xlat(i)+dlat_loc,xlat_max_sai)
        IF  ( xlat(i).GE.xlat_min_sai-dlat_loc .AND. xlat(i).LT.xlat_max_sai+dlat_loc .AND. &
          &   xlon(i).GE.xlon_sai-dlon .AND. xlon(i).LT.xlon_sai+dlon ) THEN
!
!         compute altLMDz
          altLMDz(:)=0.0
          DO k=1, klev
            zrho=pplay(i,k)/t_seri(i,k)/RD       !air density in kg/m3
            zdm(k)=(paprs(i,k)-paprs(i,k+1))/RG  !mass of layer in kg
            zdz=zdm(k)/zrho                      !thickness of layer in m
            altLMDz(k+1)=altLMDz(k)+zdz          !altitude of interface
          ENDDO

          SELECT CASE(flag_sulf_emit_distrib)

          CASE(0) ! Gaussian distribution
          !compute distribution of emission to vertical model layers (based on
          !Gaussian peak in altitude)
          f_lay_sum=0.0
               DO k=1, klev
                     f_lay_emiss(k)=0.0
                     DO i_int=1, n_int_alt
                         alt=altLMDz(k)+float(i_int)*(altLMDz(k+1)-altLMDz(k))/float(n_int_alt)
                         f_lay_emiss(k)=f_lay_emiss(k)+1./(sqrt(2.*RPI)*sigma_alt_sai)* &
                         & exp(-0.5*((alt-altemiss_sai)/sigma_alt_sai)**2.)*   &
                         & (altLMDz(k+1)-altLMDz(k))/float(n_int_alt)
                     ENDDO
                     f_lay_sum=f_lay_sum+f_lay_emiss(k)
               ENDDO

          CASE(1) ! Uniform distribution
          f_lay_sum=0.0
          ! In this case, parameter sigma_alt_vol(ieru) is considered to be half
          ! the height of the injection, centered around altemiss_sai
               DO k=1, klev
                    f_lay_emiss(k)=max(min(altemiss_sai+sigma_alt_sai,altLMDz(k+1))- &
                    & max(altemiss_sai-sigma_alt_sai,altLMDz(k)),0.)/(2.*sigma_alt_sai)
                    f_lay_sum=f_lay_sum+f_lay_emiss(k)
               ENDDO

          END SELECT ! Gaussian or uniform distribution

          !correct for step integration error
          f_lay_emiss(:)=f_lay_emiss(:)/f_lay_sum
          !emission as SO2 gas (with m(SO2)=64/32*m_aer_emiss)
          !vertically distributed emission
          DO k=1, klev
            ! stretch emission over whole year (360d)
            emission=m_aer_emiss_sai*(mSO2mol/mSatom)/m_air_gridbox(i,k)*f_lay_emiss(k)/ &
                      & FLOAT(year_len)/86400.*(sin(theta_max/180.*RPI)-sin(theta_min/180.*RPI))/ & 
                      & (sin(xlat_max_sai/180.*RPI)-sin(xlat_min_sai/180.*RPI))
            tr_seri(i,k,id_SO2_strat)=tr_seri(i,k,id_SO2_strat)+emission*pdtphys
            budg_emi_so2(i)=budg_emi_so2(i)+emission*zdm(k)*mSatom/mSO2mol
          ENDDO

!          !emission as monodisperse particles with 0.1um dry radius (BIN21)
!          !vertically distributed emission
!          DO k=1, klev
!            ! stretch emission over whole year (360d)
!            emission=m_aer_emiss*(mH2SO4mol/mSatom)/m_part_dry(21)/m_air_gridbox(i,k)*f_lay_emiss(k)/year_len/86400 
!            tr_seri(i,k,id_BIN01_strat+20)=tr_seri(i,k,id_BIN01_strat+20)+emission*pdtphys 
!            budg_emi_part(i)=budg_emi_part(i)+emission*zdm(k)*mSatom/mH2SO4mol
!          ENDDO
        ENDIF ! emission grid cell
      ENDDO ! klon loop

    END SELECT ! emission scenario (flag_sulf_emit)

!--read background concentrations of OCS and SO2 and lifetimes from input file
!--update the variables defined in phys_local_var_mod
    CALL interp_sulf_input(debutphy,pdtphys,paprs,tr_seri)

!--convert OCS to SO2 in the stratosphere
    CALL ocs_to_so2(pdtphys,tr_seri,t_seri,pplay,paprs,is_strato)

!--convert SO2 to H2SO4
    CALL so2_to_h2so4(pdtphys,tr_seri,t_seri,pplay,paprs,is_strato)

!--common routine for nucleation and condensation/evaporation with adaptive timestep
    CALL micphy_tstep(pdtphys,tr_seri,t_seri,pplay,paprs,rh,is_strato)

!--call coagulation routine 
    CALL coagulate(pdtphys,mdw,tr_seri,t_seri,pplay,dens_aer,is_strato)

!--call sedimentation routine 
    CALL aer_sedimnt(pdtphys, t_seri, pplay, paprs, tr_seri, dens_aer)

!--compute mass concentration of PM2.5 sulfate particles (wet diameter and mass) at the surface for health studies
    surf_PM25_sulf(:)=0.0
    DO i=1,klon
      DO it=1, nbtr_bin
        IF (mdw(it) .LT. 2.5e-6) THEN
          !surf_PM25_sulf(i)=surf_PM25_sulf(i)+tr_seri(i,1,it+nbtr_sulgas)*m_part(i,1,it) &
          !assume that particles consist of ammonium sulfate at the surface (132g/mol) 
          !and are dry at T = 20 deg. C and 50 perc. humidity
          surf_PM25_sulf(i)=surf_PM25_sulf(i)+tr_seri(i,1,it+nbtr_sulgas) &
                           & *132./98.*dens_aer_dry*4./3.*RPI*(mdw(it)/2.)**3 &
                           & *pplay(i,1)/t_seri(i,1)/RD*1.e9
        ENDIF
      ENDDO
    ENDDO

!    CALL minmaxsimple(tr_seri(:,:,id_SO2_strat),0.0,0.0,'fin SO2')
!    DO it=1, nbtr_bin
!      CALL minmaxsimple(tr_seri(:,:,nbtr_sulgas+it),0.0,0.0,'fin bin ')
!    ENDDO

  END SUBROUTINE traccoag

END MODULE traccoag_mod
