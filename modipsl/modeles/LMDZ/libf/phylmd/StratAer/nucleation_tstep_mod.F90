!
! $Id: nucleation_tstep_mod.F90 3690 2020-05-30 14:14:30Z oboucher $
!
MODULE nucleation_tstep_mod

CONTAINS

SUBROUTINE nucleation_rate(rhoa,t_seri,pplay,rh,a_xm,b_xm,c_xm,nucl_rate,ntot_n,x_n)

  USE aerophys
  USE infotrac_phy
  USE YOMCST, ONLY : RPI, RD, RMD, RKBOL, RNAVO

  IMPLICIT NONE

  ! input variables
  LOGICAL, PARAMETER :: flag_new_nucl=.TRUE.
  REAL rhoa    ! H2SO4 number density [molecules/cm3]
  REAL t_seri  ! temperature (K)
  REAL pplay   ! pressure (Pa)
  REAL rh      ! relative humidity (between 0 and 1)
  REAL a_xm, b_xm, c_xm

  ! output variables
  REAL nucl_rate
  REAL ntot_n ! total number of molecules in the critical cluster for neutral nucleation
  REAL x_n    ! mole fraction of H2SO4 in the critical cluster for neutral nucleation

  ! local variables
  REAL jnuc_n ! nucleation rate in 1/cm3s (10^-7-10^10 1/cm3s) for neutral nucleation
  REAL rc_n   ! radius of the critical cluster in nm for neutral nucleation
  REAL na_n   ! sulfuric acid molecules in the neutral critical cluster (NOT IN USE)
  LOGICAL kinetic_n ! true if kinetic neutral nucleation (NOT IN USE)
  LOGICAL kinetic_i ! true if kinetic ion-induced nucleation (NOT IN USE)
  REAL rhoatres ! threshold concentration of H2SO4 (1/cm^3) for neutral kinetic nucleation (NOT IN USE)
  REAL VH2SO4mol
  REAL ntot_i, x_i, jnuc_i, rc_i, na_i, n_i ! quantities for charged nucleation (NOT IN USE)
  REAL csi     ! Ion condensation sink (s-1) NOT IN USE
  REAL airn    ! Air molecule concentration in (cm-3) NOT IN USE
  REAL ipr     ! Ion pair production rate (cm-3 s-1) NOT IN USE

  ! call nucleation routine
  IF (.NOT.flag_new_nucl) THEN
    ! Use older routine from Hanna Vehkamäki (FMI)
    CALL binapara(t_seri,rh,rhoa,jnuc_n,x_n,ntot_n,rc_n)
    ! when total number of molecules is too small 
    ! then set jnuc_n to collision rate of two H2SO4 molecules (following personal communication of Ulrike Niemeier and Hanna Vehkamäki)
    IF (ntot_n < 4.0) THEN
      VH2SO4mol=mH2SO4mol/(1.E-3*(a_xm+t_seri*(b_xm+t_seri*c_xm))) !cm3
      jnuc_n = rhoa**2. *(3./4.*RPI)**(1./6.) *(12.*RKBOL*t_seri/mH2SO4mol)**0.5 &
           & *100.*(2.*VH2SO4mol**(1./3.))**2. !1/(cm3s)
      ntot_n=2.0
      x_n=1.0
    ENDIF
  ELSE
    ! Use new routine from Anni Maattanen (LATMOS)
    csi=0.0   ! no charged nucleation for now
    ipr=-1.0  ! dummy value to make sure charged nucleation does not occur
    airn=0.0  ! NOT IN USE
!   airn=pplay/t_seri/RD/1.E3*RNAVO/RMD ! molec cm-3 (for future use, to be confirmed)
    CALL newbinapara(t_seri,rh,rhoa,csi,airn,ipr,jnuc_n,ntot_n,jnuc_i,ntot_i, &
                   & x_n,x_i,na_n,na_i,rc_n,rc_i,n_i,kinetic_n,kinetic_i,rhoatres)
  ENDIF

  ! convert jnuc_n from particles/cm3/s to kg(H2SO4)/kgA/s
  nucl_rate=jnuc_n*ntot_n*x_n*mH2SO4mol/(pplay/t_seri/RD/1.E6)

END SUBROUTINE nucleation_rate

!--------------------------------------------------------------------------------------------------

SUBROUTINE nucleation_part(nucl_rate,ntot,x,dt,Vbin,tr_seri)

  USE aerophys
  USE infotrac_phy

  IMPLICIT NONE

  ! input variable
  REAL nucl_rate
  REAL ntot ! total number of molecules in the critical cluster
  REAL x    ! mole raction of H2SO4 in the critical cluster
  REAL dt
  REAL Vbin(nbtr_bin)

  ! output variable
  REAL tr_seri(nbtr)

  ! local variable
  INTEGER k
  REAL Vnew
  REAL ff(nbtr_bin)

  ! dry volume of nucleated particle (at reference temperature)
  Vnew=mH2SO4mol*ntot*x/dens_aer_dry

  ! compute distribution factor for particles of intermediate size (from Jacobson 1994, equation 13)
  ff(:)=0.0 
  DO k=1, nbtr_bin 
  ! CK 20160531: bug fix for first bin
    IF (k.LE.nbtr_bin-1) THEN 
      IF (Vbin(k).LE.Vnew.AND.Vnew.LT.Vbin(k+1)) THEN
        ff(k)= Vbin(k)/Vnew*(Vbin(k+1)-Vnew)/(Vbin(k+1)-Vbin(k))
      ENDIF
    ENDIF
    IF (k.EQ.1.AND.Vnew.LE.Vbin(k)) THEN
      ff(k)= 1.
    ENDIF
    IF (k.GT.1) THEN 
      IF (Vbin(k-1).LT.Vnew.AND.Vnew.LT.Vbin(k)) THEN
        ff(k)= 1.-ff(k-1)
      ENDIF
    ENDIF
    IF (k.EQ.nbtr_bin.AND.Vnew.GE.Vbin(k)) THEN
      ff(k)= 1.
    ENDIF
  ENDDO
  ! correction of ff for volume conservation
  DO k=1, nbtr_bin         
    ff(k)=ff(k)*Vnew/Vbin(k)
  ENDDO

  ! add nucleated H2SO4 to corresponding aerosol bin
  DO k=1, nbtr_bin
    tr_seri(k+nbtr_sulgas)=tr_seri(k+nbtr_sulgas)+nucl_rate*dt/(mH2SO4mol*ntot*x)*ff(k)
  ENDDO

END SUBROUTINE nucleation_part

!---------------------------------------------------------------------------------------------------

SUBROUTINE binapara(pt,prh,rhoa_in,jnuc,x,ntot,rc)


  !    Fortran 90 subroutine binapara
  !
  !    Calculates parametrized values of nucleation rate,
  !    mole fraction of sulphuric acid
  !    total number of particles, and the radius of the critical cluster
  !    in H2O-H2SO4 system IF temperature, saturatio ratio of water and 
  !    sulfuric acid concentration  are given. 
  !
  !    Copyright (C) 2002 Hanna Vehkamäki
  !
  !    Division of Atmospheric Sciences
  !    Department of Physical Sciences
  !    P.O. Box 64
  !    FIN-00014 University of Helsinki
  !    Finland
  !    
  !    hanna.vehkamaki@helsinki.fi

  IMPLICIT NONE 

  REAL :: pt,t     !temperature in K (190.15-300.15K)
  REAL :: prh,rh    !saturation ratio of water (0.0001-1)
  REAL,INTENT(IN) :: rhoa_in    !sulfuric acid concentration in 1/cm3 (10^4-10^11 1/cm3)
  REAL,INTENT(OUT) :: jnuc    !nucleation rate in 1/cm3s (10^-7-10^10 1/cm3s)
  REAL,INTENT(OUT) :: ntot !total number of molecules in the critical cluster (ntot>4)
  REAL,INTENT(OUT) :: x    ! mole fraction of H2SO4 in the critical cluster       
  REAL,INTENT(OUT) :: rc    !radius of the critical cluster in nm     
  REAL :: rhotres    ! threshold concentration of h2so4 (1/cm^3)
                     ! which produces nucleation rate   1/(cm^3 s) as a function of rh and t 
  REAL rhoa

! CK: use intermediate variable to avoid overwriting
  t=pt
  rh=prh
  rhoa=rhoa_in

  IF (t < 190.15) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): temperature < 190.15 K, using 190.15 K (T=',t,'K)'
     t=190.15
  ENDIF

  IF (t > 300.15) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): temperature > 300.15 K, using 300.15 K'
     t=300.15
  ENDIF

  IF (rh < 0.0001) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): saturation ratio of water < 0.0001, using 0.0001'
     rh=0.0001
  ENDIF

  IF (rh > 1.0) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): saturation ratio of water > 1 using 1'
!     print *, 'rh=',rh
     rh=1.0
  ENDIF

  IF (rhoa < 1.e4) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): sulfuric acid concentration < 1e4 1/cm3, using 1e4 1/cm3'
     rhoa=1.e4
  ENDIF

  IF (rhoa > 1.e11) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): sulfuric acid concentration > 1e11 1/cm3, using 1e11 1/cm3'
     rhoa=1.e11
  ENDIF

  x=  0.7409967177282139 - 0.002663785665140117*t + 0.002010478847383187*LOG(rh)  &
       & - 0.0001832894131464668*t*LOG(rh) + 0.001574072538464286*LOG(rh)**2      &
       & - 0.00001790589121766952*t*LOG(rh)**2 + 0.0001844027436573778*LOG(rh)**3 &
       & -  1.503452308794887E-6*t*LOG(rh)**3 - 0.003499978417957668*LOG(rhoa)    &
       & + 0.0000504021689382576*t*LOG(rhoa)

  jnuc= 0.1430901615568665 + 2.219563673425199*t - 0.02739106114964264*t**2 + &
       &  0.00007228107239317088*t**3 + 5.91822263375044/x +                  &
       &  0.1174886643003278*LOG(rh) + 0.4625315047693772*t*LOG(rh) -          &
       &  0.01180591129059253*t**2*LOG(rh) +                                  &
       &  0.0000404196487152575*t**3*LOG(rh) + (15.79628615047088*LOG(rh))/x - &
       &  0.215553951893509*LOG(rh)**2 - 0.0810269192332194*t*LOG(rh)**2 +    &
       &  0.001435808434184642*t**2*LOG(rh)**2 -                               &
       &  4.775796947178588E-6*t**3*LOG(rh)**2 -                               &
       &  (2.912974063702185*LOG(rh)**2)/x - 3.588557942822751*LOG(rh)**3 +   &
       &  0.04950795302831703*t*LOG(rh)**3 -                                   &
       &  0.0002138195118737068*t**2*LOG(rh)**3 +                             &
       &  3.108005107949533E-7*t**3*LOG(rh)**3 -                               &
       &  (0.02933332747098296*LOG(rh)**3)/x +                                &
       &  1.145983818561277*LOG(rhoa) -                                        &
       &  0.6007956227856778*t*LOG(rhoa) +                                    &
       &  0.00864244733283759*t**2*LOG(rhoa) -                                 &
       &  0.00002289467254710888*t**3*LOG(rhoa) -                              &
       &  (8.44984513869014*LOG(rhoa))/x +                                    &
       &  2.158548369286559*LOG(rh)*LOG(rhoa) +                               &
       &  0.0808121412840917*t*LOG(rh)*LOG(rhoa) -                             &
       &  0.0004073815255395214*t**2*LOG(rh)*LOG(rhoa) -                       &
       &  4.019572560156515E-7*t**3*LOG(rh)*LOG(rhoa) +                       &
       &  (0.7213255852557236*LOG(rh)*LOG(rhoa))/x +                          &
       &  1.62409850488771*LOG(rh)**2*LOG(rhoa) -                              &
       &  0.01601062035325362*t*LOG(rh)**2*LOG(rhoa) +                        &
       &  0.00003771238979714162*t**2*LOG(rh)**2*LOG(rhoa) +                  &
       &  3.217942606371182E-8*t**3*LOG(rh)**2*LOG(rhoa) -                     &
       &  (0.01132550810022116*LOG(rh)**2*LOG(rhoa))/x +                      &
       &  9.71681713056504*LOG(rhoa)**2 -                                      &
       &  0.1150478558347306*t*LOG(rhoa)**2 +                                 &
       &  0.0001570982486038294*t**2*LOG(rhoa)**2 +                           &
       &  4.009144680125015E-7*t**3*LOG(rhoa)**2 +                            &
       &  (0.7118597859976135*LOG(rhoa)**2)/x -                                &
       &  1.056105824379897*LOG(rh)*LOG(rhoa)**2 +                            &
       &  0.00903377584628419*t*LOG(rh)*LOG(rhoa)**2 -                         &
       &  0.00001984167387090606*t**2*LOG(rh)*LOG(rhoa)**2 +                  &
       &  2.460478196482179E-8*t**3*LOG(rh)*LOG(rhoa)**2 -                     &
       &  (0.05790872906645181*LOG(rh)*LOG(rhoa)**2)/x -                       &
       &  0.1487119673397459*LOG(rhoa)**3 +                                   &
       &  0.002835082097822667*t*LOG(rhoa)**3 -                                &
       &  9.24618825471694E-6*t**2*LOG(rhoa)**3 +                             &
       &  5.004267665960894E-9*t**3*LOG(rhoa)**3 -                             &
       &  (0.01270805101481648*LOG(rhoa)**3)/x
  jnuc=EXP(jnuc) !1/(cm3s)

  ntot =-0.002954125078716302 - 0.0976834264241286*t + 0.001024847927067835*t**2 - 2.186459697726116E-6*t**3 -    &
       &   0.1017165718716887/x - 0.002050640345231486*LOG(rh) - 0.007585041382707174*t*LOG(rh) +                &
       &   0.0001926539658089536*t**2*LOG(rh) - 6.70429719683894E-7*t**3*LOG(rh) -                                &
       &   (0.2557744774673163*LOG(rh))/x + 0.003223076552477191*LOG(rh)**2 + 0.000852636632240633*t*LOG(rh)**2 - &
       &   0.00001547571354871789*t**2*LOG(rh)**2 + 5.666608424980593E-8*t**3*LOG(rh)**2 +                       &
       &   (0.03384437400744206*LOG(rh)**2)/x + 0.04743226764572505*LOG(rh)**3 -                                  &
       &   0.0006251042204583412*t*LOG(rh)**3 + 2.650663328519478E-6*t**2*LOG(rh)**3 -                            &
       &   3.674710848763778E-9*t**3*LOG(rh)**3 - (0.0002672510825259393*LOG(rh)**3)/x -                          &
       &   0.01252108546759328*LOG(rhoa) + 0.005806550506277202*t*LOG(rhoa) -                                     &
       &   0.0001016735312443444*t**2*LOG(rhoa) + 2.881946187214505E-7*t**3*LOG(rhoa) +                          &
       &   (0.0942243379396279*LOG(rhoa))/x - 0.0385459592773097*LOG(rh)*LOG(rhoa) -                              &
       &   0.0006723156277391984*t*LOG(rh)*LOG(rhoa) + 2.602884877659698E-6*t**2*LOG(rh)*LOG(rhoa) +             &
       &   1.194163699688297E-8*t**3*LOG(rh)*LOG(rhoa) - (0.00851515345806281*LOG(rh)*LOG(rhoa))/x -              &
       &   0.01837488495738111*LOG(rh)**2*LOG(rhoa) + 0.0001720723574407498*t*LOG(rh)**2*LOG(rhoa) -              &
       &   3.717657974086814E-7*t**2*LOG(rh)**2*LOG(rhoa) -                                                       &
       &   5.148746022615196E-10*t**3*LOG(rh)**2*LOG(rhoa) +                                                     &
       &   (0.0002686602132926594*LOG(rh)**2*LOG(rhoa))/x - 0.06199739728812199*LOG(rhoa)**2 +                   &
       &   0.000906958053583576*t*LOG(rhoa)**2 - 9.11727926129757E-7*t**2*LOG(rhoa)**2 -                          &
       &   5.367963396508457E-9*t**3*LOG(rhoa)**2 - (0.007742343393937707*LOG(rhoa)**2)/x +                      &
       &   0.0121827103101659*LOG(rh)*LOG(rhoa)**2 - 0.0001066499571188091*t*LOG(rh)*LOG(rhoa)**2 +              &
       &   2.534598655067518E-7*t**2*LOG(rh)*LOG(rhoa)**2 -                                                       &
       &   3.635186504599571E-10*t**3*LOG(rh)*LOG(rhoa)**2 +                                                     &
       &   (0.0006100650851863252*LOG(rh)*LOG(rhoa)**2)/x + 0.0003201836700403512*LOG(rhoa)**3 -                  &
       &   0.0000174761713262546*t*LOG(rhoa)**3 + 6.065037668052182E-8*t**2*LOG(rhoa)**3 -                        &
       &   1.421771723004557E-11*t**3*LOG(rhoa)**3 + (0.0001357509859501723*LOG(rhoa)**3)/x
  ntot=EXP(ntot)

  rc=EXP(-1.6524245+0.42316402*x+0.33466487*LOG(ntot)) !nm

  IF (jnuc < 1.E-7) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev'): nucleation rate < 1E-7/cm3s, using 0.0/cm3s,'
     jnuc=0.0
  ENDIF

  IF (jnuc > 1.e10) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev, &
!        & '): nucleation rate > 1e10/cm3s, using 1e10/cm3s'
     jnuc=1.e10
  ENDIF

  rhotres=EXP( -279.2430007512709 + 11.73439886096903*rh + 22700.92970508331/t &
       & - (1088.644983466801*rh)/t + 1.144362942094912*t                      &
       & - 0.03023314602163684*rh*t - 0.001302541390154324*t**2                &
       & - 6.386965238433532*LOG(rh) + (854.980361026715*LOG(rh))/t            &
       & + 0.00879662256826497*t*LOG(rh)) !1/cm3

  RETURN

END SUBROUTINE binapara

!---------------------------------------------------------------------------------------------------

SUBROUTINE newbinapara(t,satrat,rhoa,csi,airn,ipr,jnuc_n_real,ntot_n_real,jnuc_i_real,ntot_i_real,        &
                   &   x_n_real,x_i_real,na_n_real,na_i_real,rc_n_real,rc_i_real,n_i_real,                & 
                   &   kinetic_n,kinetic_i,rhoatres_real)

  !    Fortran 90 subroutine newbinapara
  !
  !    Calculates parametrized values for neutral and ion-induced sulfuric acid-water particle formation rate
  !    of critical clusters, 
  !    number of particle in the critical clusters, the radii of the critical clusters
  !    in H2O-H2SO4-ion system if temperature, saturation ratio of water, sulfuric acid concentration, 
  !    and, optionally, either condensation sink due to pre-existing particle and ion pair production rate,
  !    or atmospheric concentration of negative ions are given. 
  !
  !    The code calculates also the kinetic limit and the particle formation rate
  !    above this limit (in which case we set ntot=1 and na=1)
  !
  !    Copyright (C)2018 Määttänen et al. 2018
  !    
  !    anni.maattanen@latmos.ipsl.fr
  !    joonas.merikanto@fmi.fi
  !    hanna.vehkamaki@helsinki.fi
  !
  !    References
  !    A. Määttänen, J. Merikanto, H. Henschel, J. Duplissy, R. Makkonen, 
  !    I. K. Ortega and H. Vehkamäki (2018), New parameterizations for 
  !    neutral and ion-induced sulfuric acid-water particle formation in 
  !    nucleation and kinetic regimes, J. Geophys. Res. Atmos., 122, doi:10.1002/2017JD027429.
  !
  !    Brasseur, G., and A.  Chatel (1983),  paper  presented  at  the  9th  Annual  Meeting  of  the  
  !    European Geophysical Society, Leeds, Great Britain, August 1982.
  !  
  !    Dunne, Eimear M., et al.(2016), Global atmospheric particle formation from CERN CLOUD measurements,
  !    Science 354.6316, 1119-1124.   
  !

  USE aerophys
  USE YOMCST, ONLY : RPI, RKBOL

  IMPLICIT NONE 

  !----------------------------------------------------
  
  !Global intent in
  REAL,INTENT(IN) :: t         ! temperature in K 
  REAL,INTENT(IN) :: satrat    ! saturatio ratio of water (between zero and 1)
  REAL,INTENT(IN) :: rhoa      ! sulfuric acid concentration in 1/cm3
  REAL,INTENT(IN) :: csi       ! Ion condensation sink (s-1)
  REAL,INTENT(IN) :: airn      ! Air molecule concentration in (cm-3)
  REAL,INTENT(IN) :: ipr       ! Ion pair production rate (cm-3 s-1)
  !Global intent out
  REAL,INTENT(OUT) :: jnuc_n_real   ! Neutral nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
  REAL,INTENT(OUT) :: ntot_n_real   ! total number of molecules in the neutral critical cluster
  REAL,INTENT(OUT) :: jnuc_i_real   ! Charged nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
  REAL,INTENT(OUT) :: ntot_i_real   ! total number of molecules in the charged critical cluster
  REAL,INTENT(OUT) :: x_n_real      ! mole fraction of H2SO4 in the neutral critical cluster 
  REAL,INTENT(OUT) :: x_i_real      ! mole fraction of H2SO4 in the charged critical cluster 
                                           ! (note that x_n=x_i in nucleation regime) 
  REAL,INTENT(OUT) :: na_n_real     ! sulfuric acid molecules in the neutral critical cluster
  REAL,INTENT(OUT) :: na_i_real     ! sulfuric molecules in the charged critical cluster
  REAL,INTENT(OUT) :: rc_n_real     ! radius of the charged critical cluster in nm 
  REAL,INTENT(OUT) :: rc_i_real     ! radius of the charged critical cluster in nm 
  REAL,INTENT(OUT) :: n_i_real      ! number of ion pairs in air (cm-3) 
  REAL,INTENT(OUT) :: rhoatres_real ! threshold concentration of H2SO4 (1/cm^3) for neutral kinetic nucleation
  LOGICAL,INTENT(OUT)  :: kinetic_n        ! true if kinetic neutral nucleation
  LOGICAL,INTENT(OUT)  :: kinetic_i        ! true if kinetic ion-induced nucleation

  ! Local
  DOUBLE PRECISION :: jnuc_n      ! Neutral nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
  DOUBLE PRECISION :: ntot_n      ! total number of molecules in the neutral critical cluster
  DOUBLE PRECISION :: jnuc_i      ! Charged nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
  DOUBLE PRECISION :: ntot_i      ! total number of molecules in the charged critical cluster
  DOUBLE PRECISION :: x_n         ! mole fraction of H2SO4 in the neutral critical cluster 
  DOUBLE PRECISION :: x_i         ! mole fraction of H2SO4 in the charged critical cluster 
                                              ! (note that x_n=x_i in nucleation regime) 
  DOUBLE PRECISION :: na_n        ! sulfuric acid molecules in the neutral critical cluster
  DOUBLE PRECISION :: na_i        ! sulfuric molecules in the charged critical cluster
  DOUBLE PRECISION :: rc_n        ! radius of the charged critical cluster in nm 
  DOUBLE PRECISION :: rc_i        ! radius of the charged critical cluster in nm 
  DOUBLE PRECISION :: n_i         ! number of ion pairs in air (cm-3) 
  DOUBLE PRECISION :: rhoatres    ! threshold concentration of H2SO4 (1/cm^3) for neutral kinetic nucleation
  DOUBLE PRECISION :: x           ! mole fraction of H2SO4 in the critical cluster 
  DOUBLE PRECISION :: satratln    ! bounded water saturation ratio for neutral case (between 5.E-6 - 1.0)
  DOUBLE PRECISION :: satratli    ! bounded water saturation ratio for ion-induced case (between 1.E-7 - 0.95)
  DOUBLE PRECISION :: rhoaln      ! bounded concentration of h2so4 for neutral case (between 10^10 - 10^19 m-3)
  DOUBLE PRECISION :: rhoali      ! bounded concentration of h2so4 for ion-induced case (between 10^10 - 10^22 m-3)
  DOUBLE PRECISION :: tln         ! bounded temperature for neutral case (between 165-400 K)
  DOUBLE PRECISION :: tli         ! bounded temperature for ion-induced case (195-400 K)
  DOUBLE PRECISION :: kinrhotresn ! threshold sulfuric acid for neutral kinetic nucleation   
  DOUBLE PRECISION :: kinrhotresi ! threshold sulfuric acid for ion-induced kinetic nucleation 
  DOUBLE PRECISION :: jnuc_i1     ! Ion-induced rate for n_i=1 cm-3 
  DOUBLE PRECISION :: xloss       ! Ion loss rate 
  DOUBLE PRECISION :: recomb      ! Ion-ion recombination rate 

  !--- 0) Initializations:

  kinetic_n=.FALSE.
  kinetic_i=.FALSE.
  jnuc_n=0.0
  jnuc_i=0.0
  ntot_n=0.0
  ntot_i=0.0
  na_n=0.0
  na_i=0.0
  rc_n=0.0
  rc_i=0.0
  x=0.0
  x_n=0.0
  x_i=0.0
  satratln=satrat
  satratli=satrat
  rhoaln=rhoa
  rhoali=rhoa
  tln=t
  tli=t
  n_i=0.0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Boundary values according to parameterization limits   

  !Temperature bounds
  IF (t.LE.165.) THEN
!     print *,'Warning: temperature < 165.0 K, using 165.0 K in neutral nucleation calculation'
     tln=165.0
  ENDIF
  IF (t.LE.195.) THEN
!     print *,'Warning: temperature < 195.0 K, using 195.0 K in ion-induced nucleation calculation'
     tli=195.0
  ENDIF
  IF (t.GE.400.) THEN
!     print *,'Warning: temperature > 400. K, using 400. K in nucleation calculations'
     tln=400.
     tli=400.
  ENDIF
  
  ! Saturation ratio bounds
  IF (satrat.LT.1.E-7) THEN
!     print *,'Warning: saturation ratio of water < 1.E-7, using 1.E-7 in ion-induced nucleation calculation'
     satratli=1.E-7
  ENDIF
  IF (satrat.LT.1.E-5) THEN
!     print *,'Warning: saturation ratio of water < 1.E-5, using 1.E-5 in neutral nucleation calculation'
     satratln=1.E-5
  ENDIF
  IF (satrat.GT.0.95) THEN
!     print *,'Warning: saturation ratio of water > 0.95, using 0.95 in ion-induced nucleation calculation'
     satratli=0.95
  ENDIF
  IF (satrat.GT.1.0) THEN
!     print *,'Warning: saturation ratio of water > 1 using 1 in neutral nucleation calculation'
     satratln=1.0
  ENDIF
  
  ! Sulfuric acid concentration bounds
  IF (rhoa.LE.1.E4) THEN
!     print *,'Warning: sulfuric acid < 1e4 1/cm3, using 1e4 1/cm3 in nucleation calculation'
     rhoaln=1.E4
     rhoali=1.E4
  ENDIF
  IF (rhoa.GT.1.E13) THEN
!     print *,'Warning: sulfuric acid > 1e13 1/cm3, using 1e13 1/cm3 in neutral nucleation calculation'
     rhoaln=1.E13
  ENDIF
  IF (rhoa.GT.1.E16) THEN
!     print *,'Warning: sulfuric acid concentration > 1e16 1/cm3, using 1e16 1/cm3 in ion-induced nucleation calculation'
     rhoali=1.E16
  ENDIF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Critical cluster composition (valid for both cases, bounds not used here) 
  x_n=  7.9036365428891719E-1 - 2.8414059650092153E-3*tln + 1.4976802556584141E-2*LOG(satratln)  &
       & - 2.4511581740839115E-4*tln*LOG(satratln) + 3.4319869471066424E-3 *LOG(satratln)**2     &
       & - 2.8799393617748428E-5*tln*LOG(satratln)**2 + 3.0174314126331765E-4*LOG(satratln)**3   &
       & - 2.2673492408841294E-6*tln*LOG(satratln)**3 - 4.3948464567032377E-3*LOG(rhoaln)        &
       & + 5.3305314722492146E-5*tln*LOG(rhoaln)
  x_i=  7.9036365428891719E-1 - 2.8414059650092153E-3*tli + 1.4976802556584141E-2*LOG(satratli)  &
       & - 2.4511581740839115E-4*tli*LOG(satratli) + 3.4319869471066424E-3 *LOG(satratli)**2     &
       & - 2.8799393617748428E-5*tli*LOG(satratli)**2 + 3.0174314126331765E-4*LOG(satratli)**3   &
       & - 2.2673492408841294E-6*tli*LOG(satratli)**3 - 4.3948464567032377E-3*LOG(rhoali)        &
       & + 5.3305314722492146E-5*tli*LOG(rhoali)
       
  x_n=MIN(MAX(x_n,1.E-30),1.) 
  x_i=MIN(MAX(x_i,1.E-30),1.) 
  
  !Neutral nucleation
  
  !Kinetic limit check
  IF (satratln .GE. 1.E-2 .AND. satratln .LE. 1.) THEN
     kinrhotresn=EXP(7.8920778706888086E+1 + 7.3665492897447082*satratln - 1.2420166571163805E+4/tln &
          & + (-6.1831234251470971E+2*satratln)/tln - 2.4501159970109945E-2*tln                      &
          & -1.3463066443605762E-2*satratln*tln + 8.3736373989909194E-06*tln**2                      &
          & -1.4673887785408892*LOG(satratln) + (-3.2141890006517094E+1*LOG(satratln))/tln           &
          & + 2.7137429081917556E-3*tln*LOG(satratln)) !1/cm3     
     IF (kinrhotresn.LT.rhoaln) kinetic_n=.TRUE.
  ENDIF

  IF (satratln .GE. 1.E-4  .AND. satratln .LT. 1.E-2) THEN     
     kinrhotresn=EXP(7.9074383049843647E+1 - 2.8746005462158347E+1*satratln - 1.2070272068458380E+4/tln &
          & + (-5.9205040320056632E+3*satratln)/tln - 2.4800372593452726E-2*tln                         &
          & -4.3983007681295948E-2*satratln*tln + 2.5943854791342071E-5*tln**2                          &
          & -2.3141363245211317*LOG(satratln) + (9.9186787997857735E+1*LOG(satratln))/tln               &
          & + 5.6819382556144681E-3*tln*LOG(satratln)) !1/cm3
     IF (kinrhotresn.LT.rhoaln) kinetic_n=.TRUE.
  ENDIF

  IF (satratln .GE. 5.E-6  .AND. satratln .LT. 1.E-4) THEN
     kinrhotresn=EXP(8.5599712000361677E+1 + 2.7335119660796581E+3*satratln - 1.1842350246291651E+4/tln &
          & + (-1.2439843468881438E+6*satratln)/tln - 5.4536964974944230E-2*tln                         &
          & + 5.0886987425326087*satratln*tln + 7.1964722655507067E-5*tln**2                            &
          & -2.4472627526306372*LOG(satratln) + (1.7561478001423779E+2*LOG(satratln))/tln               &
          & + 6.2640132818141811E-3*tln*LOG(satratln)) !1/cm3
     IF (kinrhotresn.LT.rhoaln) kinetic_n=.TRUE. 
  ENDIF
  
  IF (kinetic_n) THEN    
     ! Dimer formation rate
     jnuc_n=1.E6*(2.*0.3E-9)**2.*SQRT(8.*RPI*RKBOL*(1./mH2SO4mol+1./mH2SO4mol))/2.*SQRT(t)*rhoa**2.
!     jnuc_n=1.E6*(2.*0.3E-9)**2.*SQRT(8.*3.141593*1.38E-23*(1./(1.661E-27*98.07)+1./(1.661E-27*98.07)))/2.*SQRT(t)*rhoa**2.
     ntot_n=1. !set to 1 
     na_n=1.   ! The critical cluster contains one molecules but the produced cluster contains 2 molecules
     x_n=na_n/ntot_n  ! so also set this to 1
     rc_n=0.3E-9
  ELSE
     jnuc_n= 2.1361182605986115E-1 + 3.3827029855551838*tln -3.2423555796175563E-2*tln**2 +        &
          &  7.0120069477221989E-5*tln**3 +8.0286874752695141/x_n +                                &
          &  (-2.6939840579762231E-1)*LOG(satratln) +1.6079879299099518*tln*LOG(satratln) +        &
          &  (-1.9667486968141933E-2)*tln**2*LOG(satratln) +                                       &
          &  5.5244755979770844E-5*tln**3*LOG(satratln) + (7.8884704837892468*LOG(satratln))/x_n + &
          &  4.6374659198909596*LOG(satratln)**2 - 8.2002809894792153E-2*tln*LOG(satratln)**2 +    &
          &  8.5077424451172196E-4*tln**2*LOG(satratln)**2 +                                       &
          &  (-2.6518510168987462E-6)*tln**3*LOG(satratln)**2 +                                    &
          &  (-1.4625482500575278*LOG(satratln)**2)/x_n - 5.2413002989192037E-1*LOG(satratln)**3 + &
          &  5.2755117653715865E-3*tln*LOG(satratln)**3 +                                          &
          &  (-2.9491061332113830E-6)*tln**2*LOG(satratln)**3 +                                    &
          &  (-2.4815454194486752E-8)*tln**3*LOG(satratln)**3 +                                    &
          &  (-5.2663760117394626E-2*LOG(satratln)**3)/x_n +                                       &
          &  1.6496664658266762*LOG(rhoaln) +                                                      &
          &  (-8.0809397859218401E-1)*tln*LOG(rhoaln) +                                            &
          &  8.9302927091946642E-3*tln**2*LOG(rhoaln) +                                            &
          &  (-1.9583649496497497E-5)*tln**3*LOG(rhoaln) +                                         &
          &  (-8.9505572676891685*LOG(rhoaln))/x_n +                                               &
          &  (-3.0025283601622881E+1)*LOG(satratln)*LOG(rhoaln) +                                  &
          &  3.0783365644763633E-1*tln*LOG(satratln)*LOG(rhoaln) +                                 &
          &  (-7.4521756337984706E-4)*tln**2*LOG(satratln)*LOG(rhoaln) +                           &
          &  (-5.7651433870681853E-7)*tln**3*LOG(satratln)*LOG(rhoaln) +                           &
          &  (1.2872868529673207*LOG(satratln)*LOG(rhoaln))/x_n +                                  &
          &  (-6.1739867501526535E-1)*LOG(satratln)**2*LOG(rhoaln) +                               &
          &  7.2347385705333975E-3*tln*LOG(satratln)**2*LOG(rhoaln) +                              &
          &  (-3.0640494530822439E-5)*tln**2*LOG(satratln)**2*LOG(rhoaln) +                        &
          &  6.5944609194346214E-8*tln**3*LOG(satratln)**2*LOG(rhoaln) +                           &
          &  (-2.8681650332461055E-2*LOG(satratln)**2*LOG(rhoaln))/x_n +                           &
          &  6.5213802375160306*LOG(rhoaln)**2 +                                                   &
          &  (-4.7907162004793016E-2)*tln*LOG(rhoaln)**2 +                                         &
          &  (-1.0727890114215117E-4)*tln**2*LOG(rhoaln)**2 +                                      &
          &  5.6401818280534507E-7*tln**3*LOG(rhoaln)**2 +                                         &
          &  (5.4113070888923009E-1*LOG(rhoaln)**2)/x_n +                                          &
          &  5.2062808476476330E-1*LOG(satratln)*LOG(rhoaln)**2 +                                  &
          &  (-6.0696882500824584E-3)*tln*LOG(satratln)*LOG(rhoaln)**2 +                           &
          &  2.3851383302608477E-5*tln**2*LOG(satratln)*LOG(rhoaln)**2 +                           &
          &  (-1.5243837103067096E-8)*tln**3*LOG(satratln)*LOG(rhoaln)**2 +                        &
          &  (-5.6543192378015687E-2*LOG(satratln)*LOG(rhoaln)**2)/x_n +                           &
          &  (-1.1630806410696815E-1)*LOG(rhoaln)**3 +                                             &
          &  1.3806404273119610E-3*tln*LOG(rhoaln)**3 +                                            &
          &  (-2.0199865087650833E-6)*tln**2*LOG(rhoaln)**3 +                                      &
          &  (-3.0200284885763192E-9)*tln**3*LOG(rhoaln)**3 +                                      &
          &  (-6.9425267104126316E-3*LOG(rhoaln)**3)/x_n
     jnuc_n=EXP(jnuc_n) 
     
     ntot_n =-3.5863435141979573E-3 - 1.0098670235841110E-1*tln + 8.9741268319259721E-4*tln**2 - 1.4855098605195757E-6*tln**3  &
          &   - 1.2080330016937095E-1/x_n + 1.1902674923928015E-3*LOG(satratln) - 1.9211358507172177E-2*tln*LOG(satratln) +    &
          &   2.4648094311204255E-4*tln**2*LOG(satratln) - 7.5641448594711666E-7*tln**3*LOG(satratln) +                        &
          &   (-2.0668639384228818E-02*LOG(satratln))/x_n - 3.7593072011595188E-2*LOG(satratln)**2 +                           &
          &   9.0993182774415718E-4 *tln*LOG(satratln)**2 +                                                                    &
          &   (-9.5698412164297149E-6)*tln**2*LOG(satratln)**2 + 3.7163166416110421E-8*tln**3*LOG(satratln)**2 +               &
          &   (1.1026579525210847E-2*LOG(satratln)**2)/x_n + 1.1530844115561925E-2 *LOG(satratln)**3 +                         &
          &   (-1.8083253906466668E-4)*tln*LOG(satratln)**3 + 8.0213604053330654E-7*tln**2*LOG(satratln)**3 +                  &
          &   (-8.5797885383051337E-10)*tln**3*LOG(satratln)**3 + (1.0243693899717402E-3*LOG(satratln)**3)/x_n +               &
          &   (-1.7248695296299649E-2)*LOG(rhoaln) + 1.1294004162437157E-2*tln*LOG(rhoaln) +                                   &
          &   (-1.2283640163189278E-4)*tln**2*LOG(rhoaln) + 2.7391732258259009E-7*tln**3*LOG(rhoaln) +                         &
          &   (6.8505583974029602E-2*LOG(rhoaln))/x_n +2.9750968179523635E-1*LOG(satratln)*LOG(rhoaln) +                       &
          &   (-3.6681154503992296E-3)*tln*LOG(satratln)*LOG(rhoaln) + 1.0636473034653114E-5*tln**2*LOG(satratln)*LOG(rhoaln)+ &
          &   5.8687098466515866E-9*tln**3*LOG(satratln)*LOG(rhoaln) + (-5.2028866094191509E-3*LOG(satratln)*LOG(rhoaln))/x_n+ &
          &   7.6971988880587231E-4*LOG(satratln)**2*LOG(rhoaln) - 2.4605575820433763E-5*tln*LOG(satratln)**2*LOG(rhoaln) +    &
          &   2.3818484400893008E-7*tln**2*LOG(satratln)**2*LOG(rhoaln) +                                                      &
          &   (-8.8474102392445200E-10)*tln**3*LOG(satratln)**2*LOG(rhoaln) +                                                  &
          &   (-1.6640566678168968E-4*LOG(satratln)**2*LOG(rhoaln))/x_n - 7.7390093776705471E-2*LOG(rhoaln)**2 +               &
          &   5.8220163188828482E-4*tln*LOG(rhoaln)**2 + 1.2291679321523287E-6*tln**2*LOG(rhoaln)**2 +                         &
          &   (-7.4690997508075749E-9)*tln**3*LOG(rhoaln)**2 + (-5.6357941220497648E-3*LOG(rhoaln)**2)/x_n +                   &
          &   (-4.7170109625089768E-3)*LOG(satratln)*LOG(rhoaln)**2 + 6.9828868534370193E-5*tln*LOG(satratln)*LOG(rhoaln)**2 + &
          &   (-3.1738912157036403E-7)*tln**2*LOG(satratln)*LOG(rhoaln)**2 +                                                   &
          &   2.3975538706787416E-10*tln**3*LOG(satratln)*LOG(rhoaln)**2 +                                                     &
          &   (4.2304213386288567E-4*LOG(satratln)*LOG(rhoaln)**2)/x_n + 1.3696520973423231E-3*LOG(rhoaln)**3 +                &
          &   (-1.6863387574788199E-5)*tln*LOG(rhoaln)**3 + 2.7959499278844516E-8*tln**2*LOG(rhoaln)**3 +                      &
          &   3.9423927013227455E-11*tln**3*LOG(rhoaln)**3 + (8.6136359966337272E-5*LOG(rhoaln)**3)/x_n
     ntot_n=EXP(ntot_n)
     
     rc_n=EXP(-22.378268374023630+0.44462953606125100*x_n+0.33499495707849131*LOG(ntot_n)) !in meters
     
     na_n=x_n*ntot_n
     IF (na_n .LT. 1.) THEN
        print *, 'Warning: number of acid molecules < 1 in nucleation regime, setting na_n=1'
        na_n=1.0
     ENDIF
  ENDIF
  
  ! Set the neutral nucleation rate to 0.0 if less than 1.0E-7      
  IF (jnuc_n.LT.1.E-7) THEN
     jnuc_n=0.0
  ENDIF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Threshold neutral nucleation rate (j > 1/cm3s) parameterization (can be commented out if not needed)
  IF (tln .GE. 310.) THEN
     rhoatres=EXP(-2.8220714121794250 + 1.1492362322651116E+1*satratln -3.3034839106184218E+3/tln &
          & + (-7.1828571490168133E+2*satratln)/tln + 1.4649510835204091E-1*tln                   &
          & -3.0442736551916524E-2*satratln*tln -9.3258567137451497E-5*tln**2                     &
          & -1.1583992506895649E+1*LOG(satratln) + (1.5184848765906165E+3*LOG(satratln))/tln      &
          & + 1.8144983916747057E-2*tln*LOG(satratln)) !1/cm3
  ENDIF

  IF (tln .GT. 190. .AND. tln .LT. 310.) THEN
     rhoatres=EXP(-3.1820396091231999E+2 + 7.2451289153199676*satratln + 2.6729355170089486E+4/tln &
          & + (-7.1492506076423069E+2*satratln)/tln + 1.2617291148391978*tln                       &
          & - 1.6438112080468487E-2*satratln*tln -1.4185518234553220E-3*tln**2                     &
          & -9.2864597847386694*LOG(satratln) + (1.2607421852455602E+3*LOG(satratln))/tln          &
          & + 1.3324434472218746E-2*tln*LOG(satratln)) !1/cm3
  ENDIF

  IF (tln .LT. 185. .AND. tln .GT. 155.) THEN
     rhoatres=1.1788859232398459E+5 - 1.0244255702550814E+4*satratln +  &
          & 4.6815029684321962E+3*satratln**2 -1.6755952338499657E+2*tln
  ENDIF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Ion-induced nucleation parameterization 
 
  IF (ipr.GT.0.0) THEN ! if the ion production rate is above zero
     
     ! Calculate the ion induced nucleation rate wrt. concentration of 1 ion/cm3
     
     kinrhotresi = 5.3742280876674478E1  - 6.6837931590012266E-3 *LOG(satratli)**(-2)                                     &
          & - 1.0142598385422842E-01 * LOG(satratli)**(-1) - 6.4170597272606873E+00 * LOG(satratli)                       &
          & - 6.4315798914824518E-01 * LOG(satratli)**2 - 2.4428391714772721E-02 * LOG(satratli)**3                       &
          & - 3.5356658734539019E-04 * LOG(satratli)**4 + 2.5400015099140506E-05 * tli * LOG(satratli)**(-2)              &
          & - 2.7928900816637790E-04 * tli * LOG(satratli)**(-1) + 4.4108573484923690E-02 * tli * LOG(satratli)           &
          & + 6.3943789012475532E-03 * tli * LOG(satratli)**(2) + 2.3164296174966580E-04 * tli * LOG(satratli)**(3)       &
          & + 3.0372070669934950E-06 * tli * LOG(satratli)**4 + 3.8255873977423475E-06 * tli**2 * LOG(satratli)**(-1)     &
          & - 1.2344793083561629E-04 * tli**2 * LOG(satratli) - 1.7959048869810192E-05 * tli**2 * LOG(satratli)**(2)      &
          & - 3.2165622558722767E-07 * tli**2 * LOG(satratli)**3 - 4.7136923780988659E-09 * tli**3 * LOG(satratli)**(-1)  &
          & + 1.1873317184482216E-07 * tli**3 * LOG(satratli) + 1.5685860354866621E-08 * tli**3 * LOG(satratli)**2        &
          & - 1.4329645891059557E+04 * tli**(-1) + 1.3842599842575321E-01 * tli                                           &
          & - 4.1376265912842938E-04 * tli**(2) + 3.9147639775826004E-07 * tli**3
     
     kinrhotresi=EXP(kinrhotresi) !1/cm3
     
     IF (kinrhotresi.LT.rhoali) kinetic_i=.true.
     
     IF (kinetic_i) THEN    
!        jnuc_i1=1.0E6*(0.3E-9 + 0.487E-9)**2.*SQRT(8.*3.141593*1.38E-23*(1./(1.661E-27*98.07)+1./(1.661E-27*98.07)))*  &
        jnuc_i1=1.0E6*(0.3E-9 + 0.487E-9)**2.*SQRT(8.*RPI*RKBOL*(1./mH2SO4mol+1./mH2SO4mol))*  &
             &  SQRT(tli)*rhoali !1/cm3s  
        ntot_i=1. !set to 1 
        na_i=1.
        x_i=na_i/ntot_i  ! so also set this to 1
        rc_i=0.487E-9
     ELSE
        jnuc_i1 = 3.0108954259038608E+01+tli*6.1176722090512577E+01+(tli**2)*8.7240333618891663E-01+(tli**3)*                &
             & (-4.6191788649375719E-03)+(tli**(-1))*8.3537059107024481E-01 +                                                &
             & (1.5028549216690628E+01+tli*(-1.9310989753720623E-01)+(tli**2)*8.0155514634860480E-04+(tli**3)*               &
             & (-1.0832730707799128E-06)+(tli**(-1))*1.7577660457989019)*(LOG(satratli)**(-2)) +                             &
             & (-2.0487870170216488E-01 + tli * 1.3263949252910405E-03 + (tli**2) *(-8.4195688402450274E-06) +               &
             & (tli**3)*1.6154895940993287E-08 + (tli**(-1))*3.8734212545203874E+01) * (LOG(satratli)**(-2)*LOG(rhoali)) +   &
             & (1.4955918863858371 + tli * 9.2290004245522454E+01 + (tli**2) *(-8.9006965195392618E-01) +                    &
             & (tli**3) * 2.2319123411013099E-03 + (tli**(-1)) * 4.0180079996840852E-03) *                                   &
             & (LOG(satratli)**(-1) * LOG(rhoali)**(-1)) +                                                                   &
             & (7.9018031228561085 + tli *(-1.1649433968658949E+01) + (tli**2) * 1.1400827854910951E-01 +                    &
             & (tli**3) *(-3.1941526492127755E-04) + (tli**(-1)) *(-3.7662115740271446E-01)) * (LOG(satratli)**(-1)) +       &
             & (1.5725237111225979E+02 + tli *(-1.0051649979836277) + (tli**2) * 1.1866484014507624E-03 +                    &
             & (tli**3) * 7.3557614998540389E-06 + (tli**(-1)) * 2.6270197023115189) * (LOG(satratli)**(-1) * LOG(rhoali)) + &
             & (-1.6973840122470968E+01 + tli * 1.1258423691432135E-01 + (tli**2) *(-2.9850139351463793E-04) + (tli**3) *    &
             & 1.4301286324827064E-07 + (tli**(-1)) * 1.3163389235253725E+01) * (LOG(satratli)**(-1) * LOG(rhoali)**2) +     &
             & (-1.0399591631839757 + tli * 2.7022055588257691E-03 + (tli**2) *(-2.1507467231330936E-06) + (tli**3) *        &
             & 3.8059489037584171E-10 + (tli**(-1)) * 1.5000492788553410E+02) * (LOG(satratli)**(-1) * LOG(rhoali)**3) +     &
             & (1.2250990965305315 + tli * 3.0495946490079444E+01 + (tli**2) * 2.1051563135187106E+01 + (tli**3) *           &
             & (-8.2200682916580878E-02) + (tli**(-1)) * 2.9965871386685029E-02) * (LOG(rhoali)**(-2)) +                     &
             & (4.8281605955680433 + tli * 1.7346551710836445E+02 + (tli**2) *(-1.0113602140796010E+01) + (tli**3) *         &
             & 3.7482518458685089E-02 + (tli**(-1)) *(-1.4449998158558205E-01)) * (LOG(rhoali)**(-1)) +                      &
             & (2.3399230964451237E+02 + tli *(-2.3099267235261948E+01) + (tli**2) * 8.0122962140916354E-02 +                &
             & (tli**3) * 6.1542576994557088E-05 + (tli**(-1)) * 5.3718413254843007) * (LOG(rhoali)) +                       &
             & (1.0299715519499360E+02 + tli *(-6.4663357203364136E-02) + (tli**2) *(-2.0487150565050316E-03) +              &
             & (tli**3) * 8.7935289055530897E-07 + (tli**(-1)) * 3.6013204601215229E+01) * (LOG(rhoali)**2) +                &
             & (-3.5452115439584042 + tli * 1.7083445731159330E-02 + (tli**2) *(-1.2552625290862626E-05) + (tli**3) *        &
             & 1.2968447449182847E-09 + (tli**(-1)) * 1.5748687512056560E+02) * (LOG(rhoali)**3) +                           &
             & (2.2338490119517975 + tli * 1.0229410216045540E+02 + (tli**2) *(-3.2103611955174052) + (tli**3) *             &
             & 1.3397152304977591E-02 + (tli**(-1)) *(-2.4155187776460030E-02)) * (LOG(satratli)* LOG(rhoali)**(-2)) +       &
             & (3.7592282990713963 + tli *(-1.5257988769009816E+02) + (tli**2) * 2.6113805420558802 + (tli**3) *             &
             & (-9.0380721653694363E-03) + (tli**(-1)) *(-1.3974197138171082E-01)) * (LOG(satratli)* LOG(rhoali)**(-1)) +    &
             & (1.8293600730573988E+01 + tli * 1.8344728606002992E+01 + (tli**2) *(-4.0063363221106751E-01) + (tli**3)       &
             & * 1.4842749371258522E-03 + (tli**(-1)) * 1.1848846003282287) * (LOG(satratli)) +                              &
             & (-1.7634531623032314E+02 + tli * 4.9011762441271278 + (tli**2) *(-1.3195821562746339E-02) + (tli**3) *        &
             & (-2.8668619526430859E-05) + (tli**(-1)) *(-2.9823396976393551E-01)) * (LOG(satratli)* LOG(rhoali)) +          &
             & (-3.2944043694275727E+01 + tli * 1.2517571921051887E-01 + (tli**2) * 8.3239769771186714E-05 + (tli**3) *      &
             & 2.8191859341519507E-07 + (tli**(-1)) *(-2.7352880736682319E+01)) * (LOG(satratli)* LOG(rhoali)**2) +          &
             & (-1.1451811137553243 + tli * 2.0625997485732494E-03 + (tli**2) *(-3.4225389469233624E-06) + (tli**3) *        &
             & 4.4437613496984567E-10 + (tli**(-1)) * 1.8666644332606754E+02) * (LOG(satratli)* LOG(rhoali)**3) +            &
             & (3.2270897099493567E+01 + tli * 7.7898447327513687E-01 + (tli**2) *(-6.5662738484679626E-03) + (tli**3) *     &
             & 3.7899330796456790E-06 + (tli**(-1)) * 7.1106427501756542E-01) * (LOG(satratli)**2 * LOG(rhoali)**(-1)) +     &
             & (-2.8901906781697811E+01 + tli *(-1.5356398793054860) + (tli**2) * 1.9267271774384788E-02 + (tli**3) *        &
             & (-5.3886270475516162E-05) + (tli**(-1)) * 5.0490415975693426E-01) * (LOG(satratli)**2) +                      &
             & (3.3365683645733924E+01 + tli *(-3.6114561564894537E-01) + (tli**2) * 9.2977354471929262E-04 + (tli**3) *     &
             & 1.9549769069511355E-07 + (tli**(-1)) *(-8.8865930095112855)) * (LOG(satratli)**2 * LOG(rhoali)) +             &
             & (2.4592563042806375 + tli *(-8.3227071743101084E-03) + (tli**2) * 8.2563338043447783E-06 + (tli**3) *         &
             & (-8.4374976698593496E-09) + (tli**(-1)) *(-2.0938173949893473E+02)) * (LOG(satratli)**2 * LOG(rhoali)**2) +   &
             & (4.4099823444352317E+01 + tli * 2.5915665826835252 + (tli**2) *(-1.6449091819482634E-02) + (tli**3) *         &
             & 2.6797249816144721E-05 + (tli**(-1)) * 5.5045672663909995E-01)* satratli
        jnuc_i1=EXP(jnuc_i1)
        
        ntot_i = ABS((-4.8324296064013375E+04 + tli * 5.0469120697428906E+02 + (tli**2) *(-1.1528940488496042E+00) +         &
             & (tli**(-1)) *(-8.6892744676239192E+02) + (tli**(3)) * 4.0030302028120469E-04) +                               &
             & (-6.7259105232039847E+03 + tli * 1.9197488157452008E+02 + (tli**2) *(-1.3602976930126354E+00) +               &
             & (tli**(-1)) *(-1.1212637938360332E+02) + (tli**(3)) * 2.8515597265933207E-03) *                               &
             & LOG(satratli)**(-2) * LOG(rhoali)**(-2) +                                                                     &
             & (2.6216455217763342E+02 + tli *(-2.3687553252750821E+00) + (tli**2) * 7.4074554767517521E-03 +                &
             & (tli**(-1)) *(-1.9213956820114927E+03) + (tli**(3)) *(-9.3839114856129453E-06)) * LOG(satratli)**(-2) +       &
             & (3.9652478944137344E+00 + tli * 1.2469375098256536E-02 + (tli**2) *(-9.9837754694045633E-05) + (tli**(-1)) *  &
             & (-5.1919499210175138E+02) + (tli**(3)) * 1.6489001324583862E-07) * LOG(satratli)**(-2) * LOG(rhoali) +        &
             & (2.4975714429096206E+02 + tli * 1.7107594562445172E+02 + (tli**2) *(-7.8988711365135289E-01) + (tli**(-1)) *  &
             & (-2.2243599782483177E+01) + (tli**(3)) *(-1.6291523004095427E-04)) * LOG(satratli)**(-1) * LOG(rhoali)**(-2) +&
             & (-8.9270715592533611E+02 + tli * 1.2053538883338946E+02 + (tli**2) *(-1.5490408828541018E+00) + (tli**(-1)) * &
             & (-1.1243275579419826E+01) + (tli**(3)) * 4.8053105606904655E-03) * LOG(satratli)**(-1) * LOG(rhoali)**(-1) +  &
             & (7.6426441642091631E+03 + tli *(-7.1785462414656578E+01) + (tli**2) * 2.3851864923199523E-01 + (tli**(-1)) *  &
             & 8.5591775688708395E+01 + (tli**(3)) *(-3.7000473243342858E-04)) * LOG(satratli)**(-1) +                       &
             & (-5.1516826398607911E+01 + tli * 9.1385720811460558E-01 + (tli**2) *(-3.5477100262158974E-03) +               &
             & (tli**(-1)) * 2.7545544507625586E+03 + (tli**(3)) * 5.4708262093640928E-06) * LOG(satratli)**(-1) * LOG(rhoali) + &
             & (-3.0386767129196176E+02 + tli *(-1.1033438883583569E+04) + (tli**2) * 8.1296859732896067E+01 + (tli**(-1)) * &
             & 1.2625883141097162E+01 + (tli**(3)) *(-1.2728497822219101E-01)) * LOG(rhoali)**(-2) +                         &
             & (-3.3763494256461472E+03 + tli * 3.1916579136391006E+03 + (tli**2) *(-2.7234339474441143E+01) + (tli**(-1)) * &
             & (-2.1897653262707397E+01) + (tli**(3)) * 5.1788505812259071E-02) * LOG(rhoali)**(-1) +                        &
             & (-1.8817843873687068E+03 + tli * 4.3038072285882070E+00 + (tli**2) * 6.6244087689671860E-03 + (tli**(-1)) *   &
             & (-2.7133073605696295E+03) + (tli**(3)) *(-1.7951557394285043E-05)) * LOG(rhoali) +                            &
             & (-1.7668827539244447E+02 + tli * 4.8160932330629913E-01 + (tli**2) *(-6.3133007671100293E-04) + (tli**(-1)) * &
             & 2.5631774669873157E+04 + (tli**(3)) * 4.1534484127873519E-07) * LOG(rhoali)**(2) +                            &
             & (-1.6661835889222382E+03 + tli * 1.3708900504682877E+03 + (tli**2) *(-1.7919060052198969E+01) + (tli**(-1)) * &
             & (-3.5145029804436405E+01) + (tli**(3)) * 5.1047240947371224E-02) * LOG(satratli)* LOG(rhoali)**(-2) +         &
             & (1.0843549363030939E+04 + tli *(-7.3557073636139577E+01) + (tli**2) * 1.2054625131778862E+00 + (tli**(-1)) *  &
             & 1.9358737917864391E+02 + (tli**(3)) *(-4.2871620775911338E-03)) * LOG(satratli)* LOG(rhoali)**(-1) +          &
             & (-2.4269802549752835E+03 + tli * 1.1348265061941714E+01 + (tli**2) *(-5.0430423939495157E-02) + (tli**(-1)) * &
             & 2.3709874548950634E+03 + (tli**(3)) * 1.4091851828620244E-04) * LOG(satratli) +                               &
             & (5.2745372575251588E+02 + tli *(-2.6080675912627314E+00) + (tli**2) * 5.6902218056670145E-03 + (tli**(-1)) *  &
             & (-3.2149319482897838E+04) + (tli**(3)) *(-5.4121996056745853E-06)) * LOG(satratli)* LOG(rhoali) +             &
             & (-1.6401959518360403E+01 + tli * 2.4322962162439640E-01 + (tli**2) * 1.1744366627725344E-03 + (tli**(-1)) *   &
             & (-8.2694427518413195E+03) + (tli**(3)) *(-5.0028379203873102E-06))* LOG(satratli)**(2) +                      &
             & (-2.7556572017167782E+03 + tli * 4.9293344495058264E+01 + (tli**2) *(-2.6503456520676050E-01) + (tli**(-1)) * &
             & 1.2130698030982167E+03 + (tli**(3)) * 4.3530610668042957E-04)* LOG(satratli)**2 * LOG(rhoali)**(-1) +         &
             & (-6.3419182228959192E+00 + tli * 4.0636212834605827E-02 + (tli**2) *(-1.0450112687842742E-04) + (tli**(-1)) * &
             & 3.1035882189759656E+02 + (tli**(3)) * 9.4328418657873500E-08)* LOG(satratli)**(-3) +                          &
             & (3.0189213304689042E+03 + tli *(-2.3804654203861684E+01) + (tli**2) * 6.8113013411972942E-02 + (tli**(-1)) *  &
             & 6.3112071081188913E+02 + (tli**(3)) *(-9.4460854261685723E-05))* (satratli) * LOG(rhoali) +                   &
             & (1.1924791930673702E+04 + tli *(-1.1973824959206000E+02) + (tli**2) * 1.6888713097971020E-01 + (tli**(-1)) *  &
             & 1.8735938211539585E+02 + (tli**(3)) * 5.0974564680442852E-04)* (satratli) +                                   &
             & (3.6409071302482083E+01 + tli * 1.7919859306449623E-01 + (tli**2) *(-1.0020116255895206E-03) + (tli**(-1)) *  &
             & (-8.3521083354432303E+03) + (tli**(3)) * 1.5879900546795635E-06)* satratli * LOG(rhoali)**(2))
          
        rc_i = (-3.6318550637865524E-08 + tli * 2.1740704135789128E-09   + (tli**2) *                            &
             & (-8.5521429066506161E-12) + (tli**3) *(-9.3538647454573390E-15)) +                                &
             & (2.1366936839394922E-08 + tli *(-2.4087168827395623E-10) + (tli**2) * 8.7969869277074319E-13 +    &
             & (tli**3) *(-1.0294466881303291E-15))* LOG(satratli)**(-2) * LOG(rhoali)**(-1) +                   &
             & (-7.7804007761164303E-10 + tli * 1.0327058173517932E-11 + (tli**2) *(-4.2557697639692428E-14) +   &
             & (tli**3) * 5.4082507061618662E-17)* LOG(satratli)**(-2) +                                         &
             & (3.2628927397420860E-12 + tli *(-7.6475692919751066E-14) + (tli**2) * 4.1985816845259788E-16 +    &
             & (tli**3) *(-6.2281395889592719E-19))* LOG(satratli)**(-2) * LOG(rhoali) +                         &
             & (2.0442205540818555E-09 + tli * 4.0441858911249830E-08 + (tli**2) *(-3.3423487629482825E-10) +    &
             & (tli**3) * 6.8000404742985678E-13)* LOG(satratli)**(-1) * LOG(rhoali)**(-2) +                     &
             & (1.8381489183824627E-08 + tli *(-8.9853322951518919E-09) + (tli**2) * 7.5888799566036185E-11 +    &
             & (tli**3) *(-1.5823457864755549E-13))* LOG(satratli)**(-1) * LOG(rhoali)**(-1) +                   &
             & (1.1795760639695057E-07 + tli *(-8.1046722896375875E-10) + (tli**2) * 9.1868604369041857E-14 +    &
             & (tli**3) * 4.7882428237444610E-15)* LOG(satratli)**(-1) +                                         &
             & (-4.4028846582545952E-09 + tli * 4.6541269232626618E-11 + (tli**2) *(-1.1939929984285194E-13) +   &
             & (tli**3) * 2.3602037016614437E-17)* LOG(satratli)**(-1) * LOG(rhoali) +                           &
             & (2.7885056884209128E-11 + tli *(-4.5167129624119121E-13) + (tli**2) * 1.6558404997394422E-15 +    &
             & (tli**3) *(-1.2037336621218054E-18))* LOG(satratli)**(-1) * LOG(rhoali)**2 +                      &
             & (-2.3719627171699983E-09 + tli *(-1.5260127909292053E-07) + (tli**2) * 1.7177017944754134E-09 +   &
             & (tli**3) *(-4.7031737537526395E-12))* LOG(rhoali)**(-2) +                                         &
             & (-5.6946433724699646E-09 + tli * 8.4629788237081735E-09 + (tli**2) *(-1.7674135187061521E-10) +   &
             & (tli**3) * 6.6236547903091862E-13)* LOG(rhoali)**(-1) +                                           &
             & (-2.2808617930606012E-08 + tli * 1.4773376696847775E-10 + (tli**2) *(-1.3076953119957355E-13) +   &
             & (tli**3) * 2.3625301497914000E-16)* LOG(rhoali) +                                                 &
             & (1.4014269939947841E-10 + tli *(-2.3675117757377632E-12) + (tli**2) * 5.1514033966707879E-15 +    &
             & (tli**3) *(-4.8864233454747856E-18))* LOG(rhoali)**2 +                                            &
             & (6.5464943868885886E-11 + tli * 1.6494354816942769E-08 + (tli**2) *(-1.7480097393483653E-10) +    &
             & (tli**3) * 4.7460075628523984E-13)* LOG(satratli)* LOG(rhoali)**(-2) +                            &
             & (8.4737893183927871E-09 + tli *(-6.0243327445597118E-09) + (tli**2) * 5.8766070529814883E-11 +    &
             & (tli**3) *(-1.4926748560042018E-13))* LOG(satratli)* LOG(rhoali)**(-1) +                          &
             & (1.0761964135701397E-07 + tli *(-1.0142496009071148E-09) + (tli**2) * 2.1337312466519190E-12 +    &
             & (tli**3) * 1.6376014957685404E-15)* LOG(satratli) +                                               &
             & (-3.5621571395968670E-09 + tli * 4.1175339587760905E-11 + (tli**2) *(-1.3535372357998504E-13) +   &
             & (tli**3) * 8.9334219536920720E-17)* LOG(satratli)* LOG(rhoali) +                                  &
             & (2.0700482083136289E-11 + tli *(-3.9238944562717421E-13) + (tli**2) * 1.5850961422040196E-15 +    &
             & (tli**3) *(-1.5336775610911665E-18))* LOG(satratli)* LOG(rhoali)**2 +                             &
             & (1.8524255464416206E-09 + tli *(-2.1959816152743264E-11) + (tli**2) *(-6.4478119501677012E-14) +  &
             & (tli**3) * 5.5135243833766056E-16)* LOG(satratli)**2 * LOG(rhoali)**(-1) +                        &
             & (1.9349488650922679E-09 + tli *(-2.2647295919976428E-11) + (tli**2) * 9.2917479748268751E-14 +    &
             & (tli**3) *(-1.2741959892173170E-16))* LOG(satratli)**2 +                                          &
             & (2.1484978031650972E-11 + tli *(-9.3976642475838013E-14) + (tli**2) *(-4.8892738002751923E-16) +  &
             & (tli**3) * 1.4676120441783832E-18)* LOG(satratli)**2 * LOG(rhoali) +                              &
             & (6.7565715216420310E-13 + tli *(-3.5421162549480807E-15) + (tli**2) *(-3.4201196868693569E-18) +  &
             & (tli**3) * 2.2260187650412392E-20)* LOG(satratli)**3 * LOG(rhoali)
                    
        na_i=x_i*ntot_i
        IF (na_i .LT. 1.) THEN
!           print *, 'Warning: number of acid molecules < 1 in nucleation regime, setting na_n=1'
           na_n=1.0
        ENDIF
     ENDIF
    
     jnuc_i=jnuc_i1 
     ! Ion loss rate (1/s)
     xloss=csi+jnuc_i
     
     ! Recombination (here following Brasseur and Chatel, 1983)   
     recomb=6.0E-8*SQRT(300./tli)+6.0E-26*airn*(300./tli)**4
     
     ! Small ion concentration in air (1/cm3) (following Dunne et al., 2016)
     ! max function is to avoid n_i to go practically zero at very high J_ion 
     n_i=MAX(0.01,(SQRT(xloss**2.0+4.0*recomb*ipr)-xloss)/(2.0*recomb))
     
     ! Ion-induced nucleation rate
     ! Min function is to ensure that max function above does not cause J_ion to overshoot 
     jnuc_i=MIN(ipr,n_i*jnuc_i1)
     ! Set the ion-induced nucleation rate to 0.0 if less than 1.0E-7      
     IF (jnuc_i.LT.1.E-7) THEN
        jnuc_i=0.0
     ENDIF

  ENDIF

!--conversion from double precision to float in case the model is run in single precision
  jnuc_n_real   = jnuc_n
  ntot_n_real   = ntot_n
  jnuc_i_real   = jnuc_i
  ntot_i_real   = ntot_i
  x_n_real      = x_n
  x_i_real      = x_i
  na_n_real     = na_n
  na_i_real     = na_i
  rc_n_real     = rc_n
  rc_i_real     = rc_i
  n_i_real      = n_i
  rhoatres_real = rhoatres
  
END SUBROUTINE newbinapara

END MODULE nucleation_tstep_mod
