SUBROUTINE MIECALC_AER(tau_strat, piz_strat, cg_strat, tau_strat_wave, tau_lw_abs_rrtm, paprs, debut)

!-------Mie computations for a size distribution 
!       of homogeneous spheres.
!
!==========================================================
!--Ref : Toon and Ackerman, Applied Optics, 1981
!        Stephens, CSIRO, 1979
! Attention : surdimensionement des tableaux
! to be compiled with double precision option (-r8 on Sun)
! AUTHOR: Olivier Boucher, Christoph Kleinschmitt
!-------SIZE distribution properties----------------
!--sigma_g : geometric standard deviation 
!--r_0     : geometric number mean radius (um)/modal radius
!--Ntot    : total concentration in m-3 

  USE phys_local_var_mod, ONLY: tr_seri, mdw, alpha_bin, piz_bin, cg_bin
  USE aerophys
  USE aero_mod
  USE infotrac_phy, ONLY : nbtr, nbtr_bin, nbtr_sulgas, id_SO2_strat
  USE dimphy
  USE YOMCST  , ONLY : RG, RPI
  USE mod_phys_lmdz_para, only: gather, scatter, bcast
  USE mod_grid_phy_lmdz, ONLY : klon_glo
  USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
  USE print_control_mod, ONLY: prt_level, lunout

  IMPLICIT NONE

! Variable input
  LOGICAL,INTENT(IN) :: debut   ! le flag de l'initialisation de la physique
  REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs   ! pression pour chaque inter-couche (en Pa)

! Stratospheric aerosols optical properties
  REAL, DIMENSION(klon,klev,nbands_sw_rrtm) :: tau_strat, piz_strat, cg_strat
  REAL, DIMENSION(klon,klev,nwave_sw+nwave_lw) :: tau_strat_wave
  REAL, DIMENSION(klon,klev,nbands_lw_rrtm) :: tau_lw_abs_rrtm

!!  REAL,DIMENSION(klon_glo,klev,nbtr)     :: tr_seri_glo         ! Concentration Traceur [U/KgA]  

! local variables
  REAL Ntot
  PARAMETER (Ntot=1.0)
  LOGICAL, PARAMETER :: refr_ind_interpol = .TRUE. ! set interpolation of refractive index
  REAL r_0    ! aerosol particle radius [m]
  INTEGER bin_number, ilon, ilev
  REAL masse,volume,surface
  REAL rmin, rmax    !----integral bounds in  m

!-------------------------------------

  COMPLEX m          !----refractive index m=n_r-i*n_i
  INTEGER Nmax,Nstart !--number of iterations
  COMPLEX k2, k3, z1, z2
  COMPLEX u1,u5,u6,u8
  COMPLEX a(1:21000), b(1:21000)
  COMPLEX I 
  INTEGER n  !--loop index
  REAL nnn
  COMPLEX nn
  REAL Q_ext, Q_abs, Q_sca, g, omega   !--parameters for radius r
  REAL x, x_old  !--size parameter
  REAL r, r_lower, r_upper  !--radius
  REAL sigma_sca, sigma_ext, sigma_abs
  REAL omegatot,  gtot !--averaged parameters
  COMPLEX ksiz2(-1:21000), psiz2(1:21000)
  COMPLEX nu1z1(1:21010), nu1z2(1:21010)
  COMPLEX nu3z2(0:21000)
  REAL number, deltar
  INTEGER bin, Nbin, it
  PARAMETER (Nbin=10)
  LOGICAL smallx

!---wavelengths STREAMER
  INTEGER Nwv, NwvmaxSW
  PARAMETER (NwvmaxSW=24)
  REAL lambda(1:NwvmaxSW+1)
  DATA lambda/0.28E-6, 0.30E-6, 0.33E-6, 0.36E-6, 0.40E-6, &
              0.44E-6, 0.48E-6, 0.52E-6, 0.57E-6, 0.64E-6, &
              0.69E-6, 0.75E-6, 0.78E-6, 0.87E-6, 1.00E-6, &
              1.10E-6, 1.19E-6, 1.28E-6, 1.53E-6, 1.64E-6, &
              2.13E-6, 2.38E-6, 2.91E-6, 3.42E-6, 4.00E-6/

!---wavelengths de references
!---be careful here the 5th wavelength is 1020 nm
  INTEGER nb
  REAL lambda_ref(nwave_sw+nwave_lw)
  DATA lambda_ref /0.443E-6,0.550E-6,0.670E-6,  &
                   0.765E-6,1.020E-6,10.E-6/

!--LW 
  INTEGER NwvmaxLW
  PARAMETER (NwvmaxLW=500)
  REAL Tb, hh, cc, kb
  PARAMETER (Tb=220.0, hh=6.62607e-34)
  PARAMETER (cc=2.99792e8, kb=1.38065e-23)

!---TOA fluxes - Streamer Cs
  REAL weight(1:NwvmaxSW), weightLW
!c        DATA weight/0.839920E1, 0.231208E2, 0.322393E2, 0.465058E2,
!c     .              0.678199E2, 0.798964E2, 0.771359E2, 0.888472E2,
!c     .              0.115281E3, 0.727565E2, 0.816992E2, 0.336172E2,
!c     .              0.914603E2, 0.112706E3, 0.658840E2, 0.524470E2,
!c     .              0.391067E2, 0.883864E2, 0.276672E2, 0.681812E2,
!c     .              0.190966E2, 0.250766E2, 0.128704E2, 0.698720E1/
!---TOA fluxes - Tad
  DATA weight/ 4.20, 11.56, 16.12, 23.25, 33.91, 39.95, 38.57, &
              44.42, 57.64, 29.36, 47.87, 16.81, 45.74, 56.35, &
              32.94, 26.22, 19.55, 44.19, 13.83, 34.09,  9.55, &
              12.54,  6.44,  3.49/
!C---BOA fluxes - Tad
!c        DATA weight/ 0.01,  4.05, 9.51,  15.99, 26.07, 33.10, 33.07,
!c     .              39.91, 52.67, 27.89, 43.60, 13.67, 42.22, 40.12,
!c     .              32.70, 14.44, 19.48, 14.23, 13.43, 16.42,  8.33,
!c     .               0.95,  0.65, 2.76/

  REAL lambda_int(1:NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW), ll
  REAL dlambda_int(1:NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW), dl

  REAL n_r(1:NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW)
  REAL n_i(1:NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW)

  REAL ilambda, ilambda_prev, ilambda_max, ilambda_min
  REAL n_r_h2so4, n_i_h2so4
  REAL n_r_h2so4_prev, n_i_h2so4_prev

  REAL final_a(1:NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW)
  REAL final_g(1:NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW)
  REAL final_w(1:NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW)

  INTEGER band, bandSW, bandLW, wavenumber

!---wavelengths SW RRTM
  REAL wv_rrtm_SW(nbands_sw_rrtm+1)
  DATA wv_rrtm_SW/  0.185E-6, 0.25E-6, 0.44E-6, 0.69E-6,  &
                     1.19E-6, 2.38E-6, 4.00E-6/

!---wavenumbers and wavelengths LW RRTM
  REAL wn_rrtm(nbands_lw_rrtm+1), wv_rrtm(nbands_lw_rrtm+1)
  DATA wn_rrtm/  10.,  250.,  500.,  630.,  700.,  820.,  &
                980., 1080., 1180., 1390., 1480., 1800.,  &
               2080., 2250., 2380., 2600., 3000./

!--GCM results
  REAL gcm_a(nbands_sw_rrtm+nbands_lw_rrtm)
  REAL gcm_g(nbands_sw_rrtm+nbands_lw_rrtm)
  REAL gcm_w(nbands_sw_rrtm+nbands_lw_rrtm)
  REAL gcm_weight_a(nbands_sw_rrtm+nbands_lw_rrtm)
  REAL gcm_weight_g(nbands_sw_rrtm+nbands_lw_rrtm)
  REAL gcm_weight_w(nbands_sw_rrtm+nbands_lw_rrtm)

  REAL ss_a(nbands_sw_rrtm+nbands_lw_rrtm+nwave_sw+nwave_lw)
  REAL ss_w(nbands_sw_rrtm+nbands_lw_rrtm+nwave_sw+nwave_lw)
  REAL ss_g(nbands_sw_rrtm+nbands_lw_rrtm+nwave_sw+nwave_lw)

  !-- fichier h2so4_0.75_300.00_hummel_1988_p_q.dat
  ! -- wavenumber (cm-1), wavelength (um), n_r, n_i
  INTEGER, PARAMETER :: nb_lambda_h2so4=62
  REAL, DIMENSION (nb_lambda_h2so4,4) :: ref_ind

!---------------------------------------------------------

  IF (debut) THEN   

     ref_ind = RESHAPE( (/ &
      200.000,   50.0000,   2.01000,   6.5000E-01, &
      250.000,   40.0000,   1.94000,   6.3000E-01, & 
      285.714,   35.0000,   1.72000,   5.2000E-01, & 
      333.333,   30.0000,   1.73000,   2.9000E-01, & 
      358.423,   27.9000,   1.78000,   2.5000E-01, & 
      400.000,   25.0000,   1.84000,   2.4000E-01, &
      444.444,   22.5000,   1.82000,   2.9000E-01, &
      469.484,   21.3000,   1.79000,   2.5000E-01, &
      500.000,   20.0000,   1.81000,   2.3000E-01, &
      540.541,   18.5000,   1.92700,   3.0200E-01, &
      555.556,   18.0000,   1.95000,   4.1000E-01, &
      581.395,   17.2000,   1.72400,   5.9000E-01, &
      609.756,   16.4000,   1.52000,   4.1400E-01, &
      666.667,   15.0000,   1.59000,   2.1100E-01, &
      675.676,   14.8000,   1.61000,   2.0500E-01, &
      714.286,   14.0000,   1.64000,   1.9500E-01, &
      769.231,   13.0000,   1.69000,   1.9500E-01, &
      800.000,   12.5000,   1.74000,   1.9800E-01, &
      869.565,   11.5000,   1.89000,   3.7400E-01, &
      909.091,   11.0000,   1.67000,   4.8500E-01, &
      944.198,   10.5910,   1.72000,   3.4000E-01, &
     1000.000,   10.0000,   1.89000,   4.5500E-01, &
     1020.408,    9.8000,   1.91000,   6.8000E-01, &
     1052.632,    9.5000,   1.67000,   7.5000E-01, &
     1086.957,    9.2000,   1.60000,   5.8600E-01, &
     1111.111,    9.0000,   1.65000,   6.3300E-01, &
     1149.425,    8.7000,   1.53000,   7.7200E-01, &
     1176.471,    8.5000,   1.37000,   7.5500E-01, &
     1219.512,    8.2000,   1.20000,   6.4500E-01, &
     1265.823,    7.9000,   1.14000,   4.8800E-01, &
     1388.889,    7.2000,   1.21000,   1.7600E-01, &
     1538.462,    6.5000,   1.37000,   1.2800E-01, &
     1612.903,    6.2000,   1.42400,   1.6500E-01, &
     1666.667,    6.0000,   1.42500,   1.9500E-01, &
     1818.182,    5.5000,   1.33700,   1.8300E-01, &
     2000.000,    5.0000,   1.36000,   1.2100E-01, &
     2222.222,    4.5000,   1.38500,   1.2000E-01, &
     2500.000,    4.0000,   1.39800,   1.2600E-01, &
     2666.667,    3.7500,   1.39600,   1.3100E-01, &
     2857.143,    3.5000,   1.37600,   1.5800E-01, &
     2948.113,    3.3920,   1.35200,   1.5900E-01, &
     3125.000,    3.2000,   1.31100,   1.3500E-01, &
     3333.333,    3.0000,   1.29300,   9.5500E-02, &
     3703.704,    2.7000,   1.30300,   5.7000E-03, &
     4000.000,    2.5000,   1.34400,   3.7600E-03, &
     4444.444,    2.2500,   1.37000,   1.8000E-03, &
     5000.000,    2.0000,   1.38400,   1.2600E-03, &
     5555.556,    1.8000,   1.39000,   5.5000E-04, &
     6510.417,    1.5360,   1.40300,   1.3700E-04, &
     7692.308,    1.3000,   1.41000,   1.0000E-05, &
     9433.962,    1.0600,   1.42000,   1.5000E-06, &
    11627.907,    0.8600,   1.42500,   1.7900E-07, &
    14409.222,    0.6940,   1.42800,   1.9900E-08, &
    15797.788,    0.6330,   1.42900,   1.4700E-08, &
    18181.818,    0.5500,   1.43000,   1.0000E-08, &
    19417.476,    0.5150,   1.43100,   1.0000E-08, &
    20491.803,    0.4880,   1.43200,   1.0000E-08, &
    25000.000,    0.4000,   1.44000,   1.0000E-08, &
    29673.591,    0.3370,   1.45900,   1.0000E-08, &
    33333.333,    0.3000,   1.46900,   1.0000E-08, &
    40000.000,    0.2500,   1.48400,   1.0000E-08, &
    50000.000,    0.2000,   1.49800,   1.0000E-08 /), (/nb_lambda_h2so4,4/), order=(/2,1/) )

  !--initialising dry diameters to geometrically spaced mass/volume (see Jacobson 1994)
    mdw(1)=mdwmin
    IF (V_rat.LT.1.62) THEN ! compensate for dip in second bin for lower volume ratio
      mdw(2)=mdw(1)*2.**(1./3.)
      DO it=3, nbtr_bin
        mdw(it)=mdw(it-1)*V_rat**(1./3.)
      ENDDO 
    ELSE
      DO it=2, nbtr_bin
        mdw(it)=mdw(it-1)*V_rat**(1./3.)
      ENDDO
    ENDIF
    WRITE(lunout,*) 'init mdw=', mdw

    !--compute particle radius for a composition of 75% H2SO4 / 25% H2O at T=293K
    DO bin_number=1, nbtr_bin
      r_0=(dens_aer_dry/dens_aer_ref/0.75)**(1./3.)*mdw(bin_number)/2.
    !--integral boundaries set to bin boundaries
      rmin=r_0/sqrt(V_rat**(1./3.))
      rmax=r_0*sqrt(V_rat**(1./3.))

    !--set up SW
      DO Nwv=1, NwvmaxSW
        lambda_int(Nwv)=( lambda(Nwv)+lambda(Nwv+1) ) /2.
      ENDDO

      DO nb=1, nwave_sw+nwave_lw
        lambda_int(NwvmaxSW+nb)=lambda_ref(nb)
      ENDDO

    !--set up LW
    !--conversion wavenumber in cm-1 to wavelength in m
      DO Nwv=1, nbands_lw_rrtm+1
        wv_rrtm(Nwv)=10000./wn_rrtm(Nwv)*1.e-6
      ENDDO

      DO Nwv=1, NwvmaxLW
        lambda_int(NwvmaxSW+nwave_sw+nwave_lw+Nwv)= &
          exp( log(wv_rrtm(1))+float(Nwv-1)/float(NwvmaxLW-1)* &
          (log(wv_rrtm(nbands_lw_rrtm+1))-log(wv_rrtm(1))) )
      ENDDO

!--computing the dlambdas
      Nwv=1
      dlambda_int(NwvmaxSW+nwave_sw+nwave_lw+Nwv)= &
      &  lambda_int(NwvmaxSW+nwave_sw+nwave_lw+Nwv)- &
      &  lambda_int(NwvmaxSW+nwave_sw+nwave_lw+Nwv+1)
      DO Nwv=2, NwvmaxLW-1
      dlambda_int(NwvmaxSW+nwave_sw+nwave_lw+Nwv)= &
      &  (lambda_int(NwvmaxSW+nwave_sw+nwave_lw+Nwv-1)- &
      &  lambda_int(NwvmaxSW+nwave_sw+nwave_lw+Nwv+1))/2.
      ENDDO
      Nwv=NwvmaxLW
      dlambda_int(NwvmaxSW+nwave_sw+nwave_lw+Nwv)= &
      &  lambda_int(NwvmaxSW+nwave_sw+nwave_lw+Nwv-1)- &
      &  lambda_int(NwvmaxSW+nwave_sw+nwave_lw+Nwv)

      IF (refr_ind_interpol) THEN
 
        ilambda_max=ref_ind(1,2)/1.e6 !--in m
        ilambda_min=ref_ind(nb_lambda_h2so4,2)/1.e6 !--in m
        DO Nwv=1, NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW
          IF (lambda_int(Nwv).GT.ilambda_max) THEN
            !for lambda out of data range, take boundary values
            n_r(Nwv)=ref_ind(1,3)
            n_i(Nwv)=ref_ind(1,4)
          ELSEIF (lambda_int(Nwv).LE.ilambda_min) THEN
            n_r(Nwv)=ref_ind(nb_lambda_h2so4,3)
            n_i(Nwv)=ref_ind(nb_lambda_h2so4,4)
          ELSE
            DO nb=2,nb_lambda_h2so4
              ilambda=ref_ind(nb,2)/1.e6
              ilambda_prev=ref_ind(nb-1,2)/1.e6
              n_r_h2so4=ref_ind(nb,3)
              n_r_h2so4_prev=ref_ind(nb-1,3)
              n_i_h2so4=ref_ind(nb,4)
              n_i_h2so4_prev=ref_ind(nb-1,4)
              IF (lambda_int(Nwv).GT.ilambda.AND. &
                lambda_int(Nwv).LE.ilambda_prev) THEN !--- linear interpolation
                n_r(Nwv)=n_r_h2so4+(lambda_int(Nwv)-ilambda)/ &
                     (ilambda_prev-ilambda)*(n_r_h2so4_prev-n_r_h2so4)
                n_i(Nwv)=n_i_h2so4+(lambda_int(Nwv)-ilambda)/ &
                     (ilambda_prev-ilambda)*(n_i_h2so4_prev-n_i_h2so4)
              ENDIF
            ENDDO
          ENDIF
        ENDDO

      ELSE  !-- no refr_ind_interpol, closest neighbour from upper wavelength

        DO Nwv=1, NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW
          n_r(Nwv)=ref_ind(1,3)
          n_i(Nwv)=ref_ind(1,4)
          DO nb=2,nb_lambda_h2so4
            IF (ref_ind(nb,2)/1.e6.GT.lambda_int(Nwv)) THEN !--- step function
              n_r(Nwv)=ref_ind(nb,3)
              n_i(Nwv)=ref_ind(nb,4)
            ENDIF  
          ENDDO
        ENDDO
      ENDIF

    !---Loop on wavelengths
      DO Nwv=1, NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW

      m=CMPLX(n_r(Nwv),-n_i(Nwv))

      I=CMPLX(0.,1.)

      sigma_sca=0.0
      sigma_ext=0.0
      sigma_abs=0.0
      gtot=0.0
      omegatot=0.0
      masse = 0.0
      volume=0.0
      surface=0.0

      DO bin=1, Nbin !---loop on size bins

      r_lower=exp(log(rmin)+FLOAT(bin-1)/FLOAT(Nbin)*(log(rmax)-log(rmin)))
      r_upper=exp(log(rmin)+FLOAT(bin)/FLOAT(Nbin)*(log(rmax)-log(rmin)))
      deltar=r_upper-r_lower

      r=sqrt(r_lower*r_upper)
      x=2.*RPI*r/lambda_int(Nwv)

!we impose a minimum value for x and extrapolate quantities for small x values
      smallx = .FALSE.
      IF (x.LT.0.001) THEN
        smallx = .TRUE.
        x_old = x
        x = 0.001
      ENDIF

      number=Ntot*deltar/(rmax-rmin) !dN/dr constant over tracer bin
!      masse=masse  +4./3.*RPI*(r**3)*number*deltar*ropx*1.E3  !--g/m3
      volume=volume+4./3.*RPI*(r**3)*number*deltar
      surface=surface+4.*RPI*r**2*number*deltar

      k2=m
      k3=CMPLX(1.0,0.0)

      z2=CMPLX(x,0.0)
      z1=m*z2

      IF (0.0.LE.x.AND.x.LE.8.) THEN
        Nmax=INT(x+4*x**(1./3.)+1.)+2
      ELSEIF (8..LT.x.AND.x.LT.4200.) THEN
        Nmax=INT(x+4.05*x**(1./3.)+2.)+1
      ELSEIF (4200..LE.x.AND.x.LE.20000.) THEN
        Nmax=INT(x+4*x**(1./3.)+2.)+1
      ELSE
        WRITE(lunout,*) 'x out of bound, x=', x
        STOP
      ENDIF

      Nstart=Nmax+100

    !-----------loop for nu1z1, nu1z2

      nu1z1(Nstart)=CMPLX(0.0,0.0)
      nu1z2(Nstart)=CMPLX(0.0,0.0)
      DO n=Nstart-1, 1 , -1
        nn=CMPLX(FLOAT(n),0.0)
        nu1z1(n)=(nn+1.)/z1 - 1./( (nn+1.)/z1 + nu1z1(n+1) ) 
        nu1z2(n)=(nn+1.)/z2 - 1./( (nn+1.)/z2 + nu1z2(n+1) )
      ENDDO

    !------------loop for nu3z2

      nu3z2(0)=-I
      DO n=1, Nmax
        nn=CMPLX(FLOAT(n),0.0)
        nu3z2(n)=-nn/z2 + 1./ (nn/z2 - nu3z2(n-1) ) 
      ENDDO

    !-----------loop for psiz2 and ksiz2 (z2)
      ksiz2(-1)=COS(REAL(z2))-I*SIN(REAL(z2))
      ksiz2(0)=SIN(REAL(z2))+I*COS(REAL(z2))
      DO n=1,Nmax
       nn=CMPLX(FLOAT(n),0.0)
       ksiz2(n)=(2.*nn-1.)/z2 * ksiz2(n-1) - ksiz2(n-2)
       psiz2(n)=CMPLX(REAL(ksiz2(n)),0.0)
      ENDDO

    !-----------loop for a(n) and b(n)

      DO n=1, Nmax
        u1=k3*nu1z1(n) - k2*nu1z2(n)
        u5=k3*nu1z1(n) - k2*nu3z2(n)
        u6=k2*nu1z1(n) - k3*nu1z2(n)
        u8=k2*nu1z1(n) - k3*nu3z2(n)
        a(n)=psiz2(n)/ksiz2(n) * u1/u5
        b(n)=psiz2(n)/ksiz2(n) * u6/u8
      ENDDO

    !-----------------final loop--------------
      Q_ext=0.0
      Q_sca=0.0
      g=0.0

      DO n=Nmax-1,1,-1
        nnn=FLOAT(n)
        Q_ext=Q_ext+ (2.*nnn+1.) * REAL( a(n)+b(n) ) 
        Q_sca=Q_sca+ (2.*nnn+1.) *  &
                   REAL( a(n)*CONJG(a(n)) + b(n)*CONJG(b(n)) )
        g=g + nnn*(nnn+2.)/(nnn+1.) *   &
           REAL( a(n)*CONJG(a(n+1))+b(n)*CONJG(b(n+1)) )  +   &
              (2.*nnn+1.)/nnn/(nnn+1.) * REAL(a(n)*CONJG(b(n)))
      ENDDO

      Q_ext=2./x**2 * Q_ext
      Q_sca=2./x**2 * Q_sca
    !--extrapolation in case of small x values
      IF (smallx) THEN
        Q_ext = x_old/x * Q_ext
        Q_sca = x_old/x * Q_sca
      ENDIF

      Q_abs=Q_ext-Q_sca

      IF (AIMAG(m).EQ.0.0) Q_abs=0.0
      omega=Q_sca/Q_ext

    ! g is wrong in the smallx case (but that does not matter as long as we ignore LW scattering)
      g=g*4./x**2/Q_sca

      sigma_sca=sigma_sca+r**2*Q_sca*number
      sigma_abs=sigma_abs+r**2*Q_abs*number
      sigma_ext=sigma_ext+r**2*Q_ext*number
      omegatot=omegatot+r**2*Q_ext*omega*number
      gtot    =gtot+r**2*Q_sca*g*number
 
      ENDDO   !---bin

    !------------------------------------------------------------------

      sigma_sca=RPI*sigma_sca
      sigma_abs=RPI*sigma_abs
      sigma_ext=RPI*sigma_ext
      gtot=RPI*gtot/sigma_sca
      omegatot=RPI*omegatot/sigma_ext

      final_g(Nwv)=gtot
      final_w(Nwv)=omegatot
!      final_a(Nwv)=sigma_ext/masse
      final_a(Nwv)=sigma_ext !extinction/absorption cross section per particle

      ENDDO  !--loop on wavelength

    !---averaging over LMDZ wavebands

      DO band=1, nbands_sw_rrtm+nbands_lw_rrtm
        gcm_a(band)=0.0
        gcm_g(band)=0.0
        gcm_w(band)=0.0
        gcm_weight_a(band)=0.0
        gcm_weight_g(band)=0.0
        gcm_weight_w(band)=0.0
      ENDDO

    !---band 1 is now in the UV, so we take the first wavelength as being representative 
      DO Nwv=1,1
        bandSW=1
        gcm_a(bandSW)=gcm_a(bandSW)+final_a(Nwv)*weight(Nwv)
        gcm_weight_a(bandSW)=gcm_weight_a(bandSW)+weight(Nwv)
        gcm_w(bandSW)=gcm_w(bandSW)+final_w(Nwv)*final_a(Nwv)*weight(Nwv)
        gcm_weight_w(bandSW)=gcm_weight_w(bandSW)+final_a(Nwv)*weight(Nwv)
        gcm_g(bandSW)=gcm_g(bandSW)+final_g(Nwv)*final_a(Nwv)*final_w(Nwv)*weight(Nwv)
        gcm_weight_g(bandSW)=gcm_weight_g(bandSW)+final_a(Nwv)*final_w(Nwv)*weight(Nwv)
      ENDDO

      DO Nwv=1,NwvmaxSW

        IF (lambda_int(Nwv).LE.wv_rrtm_SW(3)) THEN      !--RRTM spectral interval 2
          bandSW=2
        ELSEIF (lambda_int(Nwv).LE.wv_rrtm_SW(4)) THEN  !--RRTM spectral interval 3
          bandSW=3
        ELSEIF (lambda_int(Nwv).LE.wv_rrtm_SW(5)) THEN  !--RRTM spectral interval 4
          bandSW=4
        ELSEIF (lambda_int(Nwv).LE.wv_rrtm_SW(6)) THEN  !--RRTM spectral interval 5
          bandSW=5
        ELSE                                            !--RRTM spectral interval 6
          bandSW=6
        ENDIF

        gcm_a(bandSW)=gcm_a(bandSW)+final_a(Nwv)*weight(Nwv)
        gcm_weight_a(bandSW)=gcm_weight_a(bandSW)+weight(Nwv)
        gcm_w(bandSW)=gcm_w(bandSW)+final_w(Nwv)*final_a(Nwv)*weight(Nwv)
        gcm_weight_w(bandSW)=gcm_weight_w(bandSW)+final_a(Nwv)*weight(Nwv)
        gcm_g(bandSW)=gcm_g(bandSW)+final_g(Nwv)*final_a(Nwv)*final_w(Nwv)*weight(Nwv)
        gcm_weight_g(bandSW)=gcm_weight_g(bandSW)+final_a(Nwv)*final_w(Nwv)*weight(Nwv)

      ENDDO

      DO Nwv=NwvmaxSW+nwave_sw+nwave_lw+1,NwvmaxSW+nwave_sw+nwave_lw+NwvmaxLW
        ll=lambda_int(Nwv)
        dl=dlambda_int(Nwv)
        weightLW=1./ll**5./(exp(hh*cc/kb/Tb/ll)-1.)*dl
        bandLW=1  !--default value starting from the highest lambda
        DO band=2, nbands_lw_rrtm
          IF (ll.LT.wv_rrtm(band)) THEN   !--as long as
            bandLW=band
          ENDIF
        ENDDO
        gcm_a(nbands_sw_rrtm+bandLW)=gcm_a(nbands_sw_rrtm+bandLW)+final_a(Nwv)*   &
             (1.-final_w(Nwv))*weightLW
        gcm_weight_a(nbands_sw_rrtm+bandLW)=gcm_weight_a(nbands_sw_rrtm+bandLW)+weightLW

        gcm_w(nbands_sw_rrtm+bandLW)=gcm_w(nbands_sw_rrtm+bandLW)+final_w(Nwv)*   &
             final_a(Nwv)*weightLW
        gcm_weight_w(nbands_sw_rrtm+bandLW)=gcm_weight_w(nbands_sw_rrtm+bandLW)+  &
             final_a(Nwv)*weightLW

        gcm_g(nbands_sw_rrtm+bandLW)=gcm_g(nbands_sw_rrtm+bandLW)+final_g(Nwv)*   &
             final_a(Nwv)*final_w(Nwv)*weightLW
        gcm_weight_g(nbands_sw_rrtm+bandLW)=gcm_weight_g(nbands_sw_rrtm+bandLW)+  &
             final_a(Nwv)*final_w(Nwv)*weightLW
      ENDDO

      DO band=1, nbands_sw_rrtm+nbands_lw_rrtm
        gcm_a(band)=gcm_a(band)/gcm_weight_a(band)
        gcm_w(band)=gcm_w(band)/gcm_weight_w(band)
        gcm_g(band)=gcm_g(band)/gcm_weight_g(band)
        ss_a(band)=gcm_a(band)
        ss_w(band)=gcm_w(band)
        ss_g(band)=gcm_g(band)
      ENDDO

      DO nb=1, nwave_sw+nwave_lw
        ss_a(nbands_sw_rrtm+nbands_lw_rrtm+nb)=final_a(NwvmaxSW+nb)
        ss_w(nbands_sw_rrtm+nbands_lw_rrtm+nb)=final_w(NwvmaxSW+nb)
        ss_g(nbands_sw_rrtm+nbands_lw_rrtm+nb)=final_g(NwvmaxSW+nb)
      ENDDO

      DO nb=1,nbands_sw_rrtm+nbands_lw_rrtm+nwave_sw+nwave_lw
        alpha_bin(nb,bin_number)=ss_a(nb) !extinction/absorption cross section per particle
        piz_bin(nb,bin_number)=ss_w(nb)
        cg_bin(nb,bin_number)=ss_g(nb)
      ENDDO

    ENDDO !loop over tracer bins

!!$OMP END MASTER
!  CALL bcast(alpha_bin)
!  CALL bcast(piz_bin)
!  CALL bcast(cg_bin)
!!$OMP BARRIER

    !set to default values at first time step (tr_seri still zero)
    tau_strat(:,:,:)=1.e-15
    piz_strat(:,:,:)=1.0
    cg_strat(:,:,:)=0.0
    tau_lw_abs_rrtm(:,:,:)=1.e-15
    tau_strat_wave(:,:,:)=1.e-15

  ELSE  !-- not debut

  !--compute optical properties of actual size distribution (from tr_seri)
    DO ilon=1,klon
      DO ilev=1, klev
        DO nb=1,nbands_sw_rrtm
          tau_strat(ilon,ilev,nb)=0.0
          DO bin_number=1, nbtr_bin
            tau_strat(ilon,ilev,nb)=tau_strat(ilon,ilev,nb)+alpha_bin(nb,bin_number) &
                                *tr_seri(ilon,ilev,bin_number+nbtr_sulgas)*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG
          ENDDO

          piz_strat(ilon,ilev,nb)=0.0
          DO bin_number=1, nbtr_bin
            piz_strat(ilon,ilev,nb)=piz_strat(ilon,ilev,nb)+piz_bin(nb,bin_number)*alpha_bin(nb,bin_number) &
                                *tr_seri(ilon,ilev,bin_number+nbtr_sulgas)*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG
          ENDDO
          piz_strat(ilon,ilev,nb)=piz_strat(ilon,ilev,nb)/MAX(tau_strat(ilon,ilev,nb),1.e-15)

          cg_strat(ilon,ilev,nb)=0.0
          DO bin_number=1, nbtr_bin
            cg_strat(ilon,ilev,nb)=cg_strat(ilon,ilev,nb)+cg_bin(nb,bin_number)*piz_bin(nb,bin_number)*alpha_bin(nb,bin_number) &
                                *tr_seri(ilon,ilev,bin_number+nbtr_sulgas)*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG
          ENDDO
          cg_strat(ilon,ilev,nb)=cg_strat(ilon,ilev,nb)/MAX(tau_strat(ilon,ilev,nb)*piz_strat(ilon,ilev,nb),1.e-15)
        ENDDO
        DO nb=1,nbands_lw_rrtm
          tau_lw_abs_rrtm(ilon,ilev,nb)=0.0
          DO bin_number=1, nbtr_bin
            tau_lw_abs_rrtm(ilon,ilev,nb)=tau_lw_abs_rrtm(ilon,ilev,nb)+alpha_bin(nbands_sw_rrtm+nb,bin_number) &
                                *tr_seri(ilon,ilev,bin_number+nbtr_sulgas)*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG
          ENDDO
        ENDDO
        DO nb=1,nwave_sw+nwave_lw
          tau_strat_wave(ilon,ilev,nb)=0.0
          DO bin_number=1, nbtr_bin
            tau_strat_wave(ilon,ilev,nb)=tau_strat_wave(ilon,ilev,nb)+alpha_bin(nbands_sw_rrtm+nbands_lw_rrtm+nb,bin_number) &
                                *tr_seri(ilon,ilev,bin_number+nbtr_sulgas)*(paprs(ilon,ilev)-paprs(ilon,ilev+1))/RG
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  ENDIF !debut

END SUBROUTINE MIECALC_AER
