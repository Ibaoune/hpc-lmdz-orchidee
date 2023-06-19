
! $Id: hines_gwd.F90 3102 2017-12-03 20:27:42Z oboucher $

SUBROUTINE hines_gwd(nlon, nlev, dtime, paphm1x, papm1x, rlat, tx, ux, vx, &
    zustrhi, zvstrhi, d_t_hin, d_u_hin, d_v_hin)

  ! ########################################################################
  ! Parametrization of the momentum flux deposition due to a broad band
  ! spectrum of gravity waves, following Hines (1997a,b), as coded by
  ! McLANDRESS (1995). Modified by McFARLANE and MANZINI (1995-1997)
  ! MAECHAM model stand alone version
  ! ########################################################################


  USE dimphy
  IMPLICIT NONE

  include "YOEGWD.h"
  include "YOMCST.h"

  INTEGER nazmth
  PARAMETER (nazmth=8)

  ! INPUT ARGUMENTS.
  ! ----- ----------

  ! - 2D
  ! PAPHM1   : HALF LEVEL PRESSURE (T-DT)
  ! PAPM1    : FULL LEVEL PRESSURE (T-DT)
  ! PTM1     : TEMPERATURE (T-DT)
  ! PUM1     : ZONAL WIND (T-DT)
  ! PVM1     : MERIDIONAL WIND (T-DT)


  ! REFERENCE.
  ! ----------
  ! SEE MODEL DOCUMENTATION

  ! AUTHOR.
  ! -------

  ! N. MCFARLANE   DKRZ-HAMBURG   MAY 1995
  ! STAND ALONE E. MANZINI MPI-HAMBURG FEBRUARY 1997

  ! BASED ON A COMBINATION OF THE OROGRAPHIC SCHEME BY N.MCFARLANE 1987
  ! AND THE HINES SCHEME AS CODED BY C. MCLANDRESS 1995.



  ! ym      INTEGER KLEVM1

  REAL paphm1(klon, klev+1), papm1(klon, klev)
  REAL ptm1(klon, klev), pum1(klon, klev), pvm1(klon, klev)
  REAL prflux(klon)
  ! 1
  ! 1
  ! 1
  REAL rlat(klon), coslat(klon)

  REAL th(klon, klev), utendgw(klon, klev), vtendgw(klon, klev), &
    pressg(klon), uhs(klon, klev), vhs(klon, klev), zpr(klon)

  ! * VERTICAL POSITIONING ARRAYS.

  REAL sgj(klon, klev), shj(klon, klev), shxkj(klon, klev), dsgj(klon, klev)

  ! * LOGICAL SWITCHES TO CONTROL ROOF DRAG, ENVELOP GW DRAG AND
  ! * HINES' DOPPLER SPREADING EXTROWAVE GW DRAG.
  ! * LOZPR IS TRUE FOR ZPR ENHANCEMENT


  ! * WORK ARRAYS.

  REAL m_alpha(klon, klev, nazmth), v_alpha(klon, klev, nazmth), &
    sigma_alpha(klon, klev, nazmth), sigsqh_alpha(klon, klev, nazmth), &
    drag_u(klon, klev), drag_v(klon, klev), flux_u(klon, klev), &
    flux_v(klon, klev), heat(klon, klev), diffco(klon, klev), &
    bvfreq(klon, klev), density(klon, klev), sigma_t(klon, klev), &
    visc_mol(klon, klev), alt(klon, klev), sigsqmcw(klon, klev, nazmth), &
    sigmatm(klon, klev), ak_alpha(klon, nazmth), k_alpha(klon, nazmth), &
    mmin_alpha(klon, nazmth), i_alpha(klon, nazmth), rmswind(klon), &
    bvfbot(klon), densbot(klon)
  REAL smoothr1(klon, klev), smoothr2(klon, klev)
  REAL sigalpmc(klon, klev, nazmth)
  REAL f2mod(klon, klev)

  ! * THES ARE THE INPUT PARAMETERS FOR HINES ROUTINE AND
  ! * ARE SPECIFIED IN ROUTINE HINES_SETUP. SINCE THIS IS CALLED
  ! * ONLY AT FIRST CALL TO THIS ROUTINE THESE VARIABLES MUST BE SAVED
  ! * FOR USE AT SUBSEQUENT CALLS. THIS CAN BE AVOIDED BY CALLING
  ! * HINES_SETUP IN MAIN PROGRAM AND PASSING THE PARAMETERS AS
  ! * SUBROUTINE ARGUEMENTS.


  REAL rmscon
  INTEGER nmessg, iprint, ilrms
  INTEGER ifl

  INTEGER naz, icutoff, nsmax, iheatcal
  REAL slope, f1, f2, f3, f5, f6, kstar(klon), alt_cutoff, smco

  ! PROVIDED AS INPUT

  INTEGER nlon, nlev

  REAL dtime
  REAL paphm1x(nlon, nlev+1), papm1x(nlon, nlev)
  REAL ux(nlon, nlev), vx(nlon, nlev), tx(nlon, nlev)

  ! VARIABLES FOR OUTPUT


  REAL d_t_hin(nlon, nlev), d_u_hin(nlon, nlev), d_v_hin(nlon, nlev)
  REAL zustrhi(nlon), zvstrhi(nlon)


  ! * LOGICAL SWITCHES TO CONTROL PRECIP ENHANCEMENT AND
  ! * HINES' DOPPLER SPREADING EXTROWAVE GW DRAG.
  ! * LOZPR IS TRUE FOR ZPR ENHANCEMENT

  LOGICAL lozpr, lorms(klon)

  ! LOCAL PARAMETERS TO MAKE THINGS WORK (TEMPORARY VARIABLE)

  REAL rhoh2o, zpcons, rgocp, zlat, dttdsf, ratio, hscal
  INTEGER i, j, l, jl, jk, le, lref, lrefp, levbot

  ! DATA PARAMETERS NEEDED, EXPLAINED LATER

  REAL v0, vmin, dmpscal, taufac, hmin, apibt, cpart, fcrit
  REAL pcrit, pcons
  INTEGER iplev, ierror



  ! PRINT *,' IT IS STARTED HINES GOING ON...'




  ! *    COMPUTATIONAL CONSTANTS.
  ! ------------- ----------


  d_t_hin(:, :) = 0.

  rhoh2o = 1000.
  zpcons = (1000.*86400.)/rhoh2o
  ! ym      KLEVM1=KLEV-1


  DO jl = kidia, kfdia
    paphm1(jl, 1) = paphm1x(jl, klev+1)
    DO jk = 1, klev
      le = klev + 1 - jk
      paphm1(jl, jk+1) = paphm1x(jl, le)
      papm1(jl, jk) = papm1x(jl, le)
      ptm1(jl, jk) = tx(jl, le)
      pum1(jl, jk) = ux(jl, le)
      pvm1(jl, jk) = vx(jl, le)
    END DO
  END DO

  ! Define constants and arrays needed for the ccc/mam gwd scheme
  ! *Constants:

  rgocp = rd/rcpd
  lrefp = klev - 1
  lref = klev - 2
  ! 1
  ! 1    *Arrays
  ! 1
  DO jk = 1, klev
    DO jl = kidia, kfdia
      shj(jl, jk) = papm1(jl, jk)/paphm1(jl, klev+1)
      sgj(jl, jk) = papm1(jl, jk)/paphm1(jl, klev+1)
      dsgj(jl, jk) = (paphm1(jl,jk+1)-paphm1(jl,jk))/paphm1(jl, klev+1)
      shxkj(jl, jk) = (papm1(jl,jk)/paphm1(jl,klev+1))**rgocp
      th(jl, jk) = ptm1(jl, jk)
    END DO
  END DO

  ! C
  DO jl = kidia, kfdia
    pressg(jl) = paphm1(jl, klev+1)
  END DO


  DO jl = kidia, kfdia
    prflux(jl) = 0.0
    zpr(jl) = zpcons*prflux(jl)
    zlat = (rlat(jl)/180.)*rpi
    coslat(jl) = cos(zlat)
  END DO

  ! /#########################################################################
  ! /
  ! /

  ! * AUG. 14/95 - C. MCLANDRESS.
  ! * SEP.    95   N. MCFARLANE.

  ! * THIS ROUTINE CALCULATES THE HORIZONTAL WIND TENDENCIES
  ! * DUE TO MCFARLANE'S OROGRAPHIC GW DRAG SCHEME, HINES'
  ! * DOPPLER SPREAD SCHEME FOR "EXTROWAVES" AND ADDS ON
  ! * ROOF DRAG. IT IS BASED ON THE ROUTINE GWDFLX8.

  ! * LREFP IS THE INDEX OF THE MODEL LEVEL BELOW THE REFERENCE LEVEL
  ! * I/O ARRAYS PASSED FROM MAIN.
  ! * (PRESSG = SURFACE PRESSURE)




  ! * CONSTANTS VALUES DEFINED IN DATA STATEMENT ARE :
  ! * VMIN     = MIMINUM WIND IN THE DIRECTION OF REFERENCE LEVEL
  ! *            WIND BEFORE WE CONSIDER BREAKING TO HAVE OCCURED.
  ! * DMPSCAL  = DAMPING TIME FOR GW DRAG IN SECONDS.
  ! * TAUFAC   = 1/(LENGTH SCALE).
  ! * HMIN     = MIMINUM ENVELOPE HEIGHT REQUIRED TO PRODUCE GW DRAG.
  ! * V0       = VALUE OF WIND THAT APPROXIMATES ZERO.


  DATA vmin/5.0/, v0/1.E-10/, taufac/5.E-6/, hmin/40000./, dmpscal/6.5E+6/, &
    apibt/1.5708/, cpart/0.7/, fcrit/1./

  ! * HINES EXTROWAVE GWD CONSTANTS DEFINED IN DATA STATEMENT ARE:
  ! * RMSCON = ROOT MEAN SQUARE GRAVITY WAVE WIND AT LOWEST LEVEL (M/S).
  ! * NMESSG  = UNIT NUMBER FOR PRINTED MESSAGES.
  ! * IPRINT  = 1 TO DO PRINT OUT SOME HINES ARRAYS.
  ! * IFL     = FIRST CALL FLAG TO HINES_SETUP ("SAVE" IT)
  ! * PCRIT = CRITICAL VALUE OF ZPR (MM/D)
  ! * IPLEV = LEVEL OF APPLICATION OF PRCIT
  ! * PCONS = FACTOR OF ZPR ENHANCEMENT


  DATA pcrit/5./, pcons/4.75/

  iplev = lrefp - 1

  DATA rmscon/1.00/iprint/2/, nmessg/6/
  DATA ifl/0/

  lozpr = .FALSE.

  ! -----------------------------------------------------------------------



  ! * SET ERROR FLAG

  ierror = 0

  ! * SPECIFY VARIOUS PARAMETERS FOR HINES ROUTINE AT VERY FIRST CALL.
  ! * (NOTE THAT ARRAY K_ALPHA IS SPECIFIED SO MAKE SURE THAT
  ! * IT IS NOT OVERWRITTEN LATER ON).

  CALL hines_setup(naz, slope, f1, f2, f3, f5, f6, kstar, icutoff, &
    alt_cutoff, smco, nsmax, iheatcal, k_alpha, ierror, nmessg, klon, nazmth, &
    coslat)
  IF (ierror/=0) GO TO 999

  ! * START GWD CALCULATIONS.

  lref = lrefp - 1


  DO j = 1, nazmth
    DO l = 1, klev
      DO i = kidia, klon
        sigsqmcw(i, l, j) = 0.
      END DO
    END DO
  END DO



  ! * INITIALIZE NECESSARY ARRAYS.

  DO l = 1, klev
    DO i = kidia, kfdia
      utendgw(i, l) = 0.
      vtendgw(i, l) = 0.

      uhs(i, l) = 0.
      vhs(i, l) = 0.

    END DO
  END DO

  ! * IF USING HINES SCHEME THEN CALCULATE B V FREQUENCY AT ALL POINTS
  ! * AND SMOOTH BVFREQ.

  DO l = 2, klev
    DO i = kidia, kfdia
      dttdsf = (th(i,l)/shxkj(i,l)-th(i,l-1)/shxkj(i,l-1))/ &
        (shj(i,l)-shj(i,l-1))
      dttdsf = min(dttdsf, -5./sgj(i,l))
      bvfreq(i, l) = sqrt(-dttdsf*sgj(i,l)*(sgj(i,l)**rgocp)/rd)*rg/ptm1(i, l &
        )
    END DO
  END DO
  DO l = 1, klev
    DO i = kidia, kfdia
      IF (l==1) THEN
        bvfreq(i, l) = bvfreq(i, l+1)
      END IF
      IF (l>1) THEN
        ratio = 5.*log(sgj(i,l)/sgj(i,l-1))
        bvfreq(i, l) = (bvfreq(i,l-1)+ratio*bvfreq(i,l))/(1.+ratio)
      END IF
    END DO
  END DO

  ! * CALCULATE GW DRAG DUE TO HINES' EXTROWAVES
  ! * SET MOLECULAR VISCOSITY TO A VERY SMALL VALUE.
  ! * IF THE MODEL TOP IS GREATER THAN 100 KM THEN THE ACTUAL
  ! * VISCOSITY COEFFICIENT COULD BE SPECIFIED HERE.

  DO l = 1, klev
    DO i = kidia, kfdia
      visc_mol(i, l) = 1.5E-5
      drag_u(i, l) = 0.
      drag_v(i, l) = 0.
      flux_u(i, l) = 0.
      flux_v(i, l) = 0.
      heat(i, l) = 0.
      diffco(i, l) = 0.
    END DO
  END DO

  ! * ALTITUDE AND DENSITY AT BOTTOM.

  DO i = kidia, kfdia
    hscal = rd*ptm1(i, klev)/rg
    density(i, klev) = sgj(i, klev)*pressg(i)/(rg*hscal)
    alt(i, klev) = 0.
  END DO

  ! * ALTITUDE AND DENSITY AT REMAINING LEVELS.

  DO l = klev - 1, 1, -1
    DO i = kidia, kfdia
      hscal = rd*ptm1(i, l)/rg
      alt(i, l) = alt(i, l+1) + hscal*dsgj(i, l)/sgj(i, l)
      density(i, l) = sgj(i, l)*pressg(i)/(rg*hscal)
    END DO
  END DO


  ! * INITIALIZE SWITCHES FOR HINES GWD CALCULATION

  ilrms = 0

  DO i = kidia, kfdia
    lorms(i) = .FALSE.
  END DO


  ! * DEFILE BOTTOM LAUNCH LEVEL

  levbot = iplev

  ! * BACKGROUND WIND MINUS VALUE AT BOTTOM LAUNCH LEVEL.

  DO l = 1, levbot
    DO i = kidia, kfdia
      uhs(i, l) = pum1(i, l) - pum1(i, levbot)
      vhs(i, l) = pvm1(i, l) - pvm1(i, levbot)
    END DO
  END DO

  ! * SPECIFY ROOT MEAN SQUARE WIND AT BOTTOM LAUNCH LEVEL.

  DO i = kidia, kfdia
    rmswind(i) = rmscon
  END DO

  IF (lozpr) THEN
    DO i = kidia, kfdia
      IF (zpr(i)>pcrit) THEN
        rmswind(i) = rmscon + ((zpr(i)-pcrit)/zpr(i))*pcons
      END IF
    END DO
  END IF

  DO i = kidia, kfdia
    IF (rmswind(i)>0.0) THEN
      ilrms = ilrms + 1
      lorms(i) = .TRUE.
    END IF
  END DO

  ! * CALCULATE GWD (NOTE THAT DIFFUSION COEFFICIENT AND
  ! * HEATING RATE ONLY CALCULATED IF IHEATCAL = 1).

  IF (ilrms>0) THEN

    CALL hines_extro0(drag_u, drag_v, heat, diffco, flux_u, flux_v, uhs, vhs, &
      bvfreq, density, visc_mol, alt, rmswind, k_alpha, m_alpha, v_alpha, &
      sigma_alpha, sigsqh_alpha, ak_alpha, mmin_alpha, i_alpha, sigma_t, &
      densbot, bvfbot, 1, iheatcal, icutoff, iprint, nsmax, smco, alt_cutoff, &
      kstar, slope, f1, f2, f3, f5, f6, naz, sigsqmcw, sigmatm, kidia, klon, &
      1, levbot, klon, klev, nazmth, lorms, smoothr1, smoothr2, sigalpmc, &
      f2mod)

    ! * ADD ON HINES' GWD TENDENCIES TO OROGRAPHIC TENDENCIES AND
    ! * APPLY HINES' GW DRAG ON (UROW,VROW) WORK ARRAYS.

    DO l = 1, klev
      DO i = kidia, kfdia
        utendgw(i, l) = utendgw(i, l) + drag_u(i, l)
        vtendgw(i, l) = vtendgw(i, l) + drag_v(i, l)
      END DO
    END DO


    ! * END OF HINES CALCULATIONS.

  END IF

  ! -----------------------------------------------------------------------

  DO jl = kidia, kfdia
    zustrhi(jl) = flux_u(jl, 1)
    zvstrhi(jl) = flux_v(jl, 1)
    DO jk = 1, klev
      le = klev - jk + 1
      d_u_hin(jl, jk) = utendgw(jl, le)*dtime
      d_v_hin(jl, jk) = vtendgw(jl, le)*dtime
    END DO
  END DO

  ! PRINT *,'UTENDGW:',UTENDGW

  ! PRINT *,' HINES HAS BEEN COMPLETED (LONG ISNT IT...)'

  RETURN
999 CONTINUE

  ! * IF ERROR DETECTED THEN ABORT.

  WRITE (nmessg, 6000)
  WRITE (nmessg, 6010) ierror
6000 FORMAT (/' EXECUTION ABORTED IN GWDOREXV')
6010 FORMAT ('     ERROR FLAG =', I4)


  RETURN
END SUBROUTINE hines_gwd
! /
! /


SUBROUTINE hines_extro0(drag_u, drag_v, heat, diffco, flux_u, flux_v, vel_u, &
    vel_v, bvfreq, density, visc_mol, alt, rmswind, k_alpha, m_alpha, &
    v_alpha, sigma_alpha, sigsqh_alpha, ak_alpha, mmin_alpha, i_alpha, &
    sigma_t, densb, bvfb, iorder, iheatcal, icutoff, iprint, nsmax, smco, &
    alt_cutoff, kstar, slope, f1, f2, f3, f5, f6, naz, sigsqmcw, sigmatm, &
    il1, il2, lev1, lev2, nlons, nlevs, nazmth, lorms, smoothr1, smoothr2, &
    sigalpmc, f2mod)

  IMPLICIT NONE

  ! Main routine for Hines' "extrowave" gravity wave parameterization based
  ! on Hines' Doppler spread theory. This routine calculates zonal
  ! and meridional components of gravity wave drag, heating rates
  ! and diffusion coefficient on a longitude by altitude grid.
  ! No "mythical" lower boundary region calculation is made so it
  ! is assumed that lowest level winds are weak (i.e, approximately zero).

  ! Aug. 13/95 - C. McLandress
  ! SEPT. /95  - N.McFarlane

  ! Modifications:

  ! Output arguements:

  ! * DRAG_U = zonal component of gravity wave drag (m/s^2).
  ! * DRAG_V = meridional component of gravity wave drag (m/s^2).
  ! * HEAT   = gravity wave heating (K/sec).
  ! * DIFFCO = diffusion coefficient (m^2/sec)
  ! * FLUX_U = zonal component of vertical momentum flux (Pascals)
  ! * FLUX_V = meridional component of vertical momentum flux (Pascals)

  ! Input arguements:

  ! * VEL_U      = background zonal wind component (m/s).
  ! * VEL_V      = background meridional wind component (m/s).
  ! * BVFREQ     = background Brunt Vassala frequency (radians/sec).
  ! * DENSITY    = background density (kg/m^3)
  ! * VISC_MOL   = molecular viscosity (m^2/s)
  ! * ALT        = altitude of momentum, density, buoyancy levels (m)
  ! *              (NOTE: levels ordered so that ALT(I,1) > ALT(I,2), etc.)
  ! * RMSWIND   = root mean square gravity wave wind at lowest level (m/s).
  ! * K_ALPHA    = horizontal wavenumber of each azimuth (1/m).
  ! * IORDER	   = 1 means vertical levels are indexed from top down
  ! *              (i.e., highest level indexed 1 and lowest level NLEVS);
  ! *           .NE. 1 highest level is index NLEVS.
  ! * IHEATCAL   = 1 to calculate heating rates and diffusion coefficient.
  ! * IPRINT     = 1 to print out various arrays.
  ! * ICUTOFF    = 1 to exponentially damp GWD, heating and diffusion
  ! *              arrays above ALT_CUTOFF; otherwise arrays not modified.
  ! * ALT_CUTOFF = altitude in meters above which exponential decay applied.
  ! * SMCO       = smoothing factor used to smooth cutoff vertical
  ! *              wavenumbers and total rms winds in vertical direction
  ! *              before calculating drag or heating
  ! *              (SMCO >= 1 ==> 1:SMCO:1 stencil used).
  ! * NSMAX      = number of times smoother applied ( >= 1),
  ! *            = 0 means no smoothing performed.
  ! * KSTAR      = typical gravity wave horizontal wavenumber (1/m).
  ! * SLOPE      = slope of incident vertical wavenumber spectrum
  ! *              (SLOPE must equal 1., 1.5 or 2.).
  ! * F1 to F6   = Hines's fudge factors (F4 not needed since used for
  ! *              vertical flux of vertical momentum).
  ! * NAZ        = actual number of horizontal azimuths used.
  ! * IL1        = first longitudinal index to use (IL1 >= 1).
  ! * IL2        = last longitudinal index to use (IL1 <= IL2 <= NLONS).
  ! * LEV1       = index of first level for drag calculation.
  ! * LEV2       = index of last level for drag calculation
  ! *              (i.e., LEV1 < LEV2 <= NLEVS).
  ! * NLONS      = number of longitudes.
  ! * NLEVS      = number of vertical levels.
  ! * NAZMTH     = azimuthal array dimension (NAZMTH >= NAZ).

  ! Work arrays.

  ! * M_ALPHA      = cutoff vertical wavenumber (1/m).
  ! * V_ALPHA      = wind component at each azimuth (m/s) and if IHEATCAL=1
  ! *                holds vertical derivative of cutoff wavenumber.
  ! * SIGMA_ALPHA  = total rms wind in each azimuth (m/s).
  ! * SIGSQH_ALPHA = portion of wind variance from waves having wave
  ! *                normals in the alpha azimuth (m/s).
  ! * SIGMA_T      = total rms horizontal wind (m/s).
  ! * AK_ALPHA     = spectral amplitude factor at each azimuth
  ! *                (i.e.,{AjKj}) in m^4/s^2.
  ! * I_ALPHA      = Hines' integral.
  ! * MMIN_ALPHA   = minimum value of cutoff wavenumber.
  ! * DENSB        = background density at bottom level.
  ! * BVFB         = buoyancy frequency at bottom level and
  ! *                work array for ICUTOFF = 1.

  ! * LORMS       = .TRUE. for drag computation

  INTEGER naz, nlons, nlevs, nazmth, il1, il2, lev1, lev2
  INTEGER icutoff, nsmax, iorder, iheatcal, iprint
  REAL kstar(nlons), f1, f2, f3, f5, f6, slope
  REAL alt_cutoff, smco
  REAL drag_u(nlons, nlevs), drag_v(nlons, nlevs)
  REAL heat(nlons, nlevs), diffco(nlons, nlevs)
  REAL flux_u(nlons, nlevs), flux_v(nlons, nlevs)
  REAL vel_u(nlons, nlevs), vel_v(nlons, nlevs)
  REAL bvfreq(nlons, nlevs), density(nlons, nlevs)
  REAL visc_mol(nlons, nlevs), alt(nlons, nlevs)
  REAL rmswind(nlons), bvfb(nlons), densb(nlons)
  REAL sigma_t(nlons, nlevs), sigsqmcw(nlons, nlevs, nazmth)
  REAL sigma_alpha(nlons, nlevs, nazmth), sigmatm(nlons, nlevs)
  REAL sigsqh_alpha(nlons, nlevs, nazmth)
  REAL m_alpha(nlons, nlevs, nazmth), v_alpha(nlons, nlevs, nazmth)
  REAL ak_alpha(nlons, nazmth), k_alpha(nlons, nazmth)
  REAL mmin_alpha(nlons, nazmth), i_alpha(nlons, nazmth)
  REAL smoothr1(nlons, nlevs), smoothr2(nlons, nlevs)
  REAL sigalpmc(nlons, nlevs, nazmth)
  REAL f2mod(nlons, nlevs)

  LOGICAL lorms(nlons)

  ! Internal variables.

  INTEGER levbot, levtop, i, n, l, lev1p, lev2m
  INTEGER ilprt1, ilprt2
  ! -----------------------------------------------------------------------

  ! PRINT *,' IN HINES_EXTRO0'
  lev1p = lev1 + 1
  lev2m = lev2 - 1

  ! Index of lowest altitude level (bottom of drag calculation).

  levbot = lev2
  levtop = lev1
  IF (iorder/=1) THEN
    WRITE (6, 1)
1   FORMAT (2X, ' error: IORDER NOT ONE! ')
  END IF

  ! Buoyancy and density at bottom level.

  DO i = il1, il2
    bvfb(i) = bvfreq(i, levbot)
    densb(i) = density(i, levbot)
  END DO

  ! initialize some variables

  DO n = 1, naz
    DO l = lev1, lev2
      DO i = il1, il2
        m_alpha(i, l, n) = 0.0
      END DO
    END DO
  END DO
  DO l = lev1, lev2
    DO i = il1, il2
      sigma_t(i, l) = 0.0
    END DO
  END DO
  DO n = 1, naz
    DO i = il1, il2
      i_alpha(i, n) = 0.0
    END DO
  END DO

  ! Compute azimuthal wind components from zonal and meridional winds.

  CALL hines_wind(v_alpha, vel_u, vel_v, naz, il1, il2, lev1, lev2, nlons, &
    nlevs, nazmth)

  ! Calculate cutoff vertical wavenumber and velocity variances.

  CALL hines_wavnum(m_alpha, sigma_alpha, sigsqh_alpha, sigma_t, ak_alpha, &
    v_alpha, visc_mol, density, densb, bvfreq, bvfb, rmswind, i_alpha, &
    mmin_alpha, kstar, slope, f1, f2, f3, naz, levbot, levtop, il1, il2, &
    nlons, nlevs, nazmth, sigsqmcw, sigmatm, lorms, sigalpmc, f2mod)
  ! Smooth cutoff wavenumbers and total rms velocity in the vertical
  ! direction NSMAX times, using FLUX_U as temporary work array.

  IF (nsmax>0) THEN
    DO n = 1, naz
      DO l = lev1, lev2
        DO i = il1, il2
          smoothr1(i, l) = m_alpha(i, l, n)
        END DO
      END DO
      CALL vert_smooth(smoothr1, smoothr2, smco, nsmax, il1, il2, lev1, lev2, &
        nlons, nlevs)
      DO l = lev1, lev2
        DO i = il1, il2
          m_alpha(i, l, n) = smoothr1(i, l)
        END DO
      END DO
    END DO
    CALL vert_smooth(sigma_t, smoothr2, smco, nsmax, il1, il2, lev1, lev2, &
      nlons, nlevs)
  END IF

  ! Calculate zonal and meridional components of the
  ! momentum flux and drag.

  CALL hines_flux(flux_u, flux_v, drag_u, drag_v, alt, density, densb, &
    m_alpha, ak_alpha, k_alpha, slope, naz, il1, il2, lev1, lev2, nlons, &
    nlevs, nazmth, lorms)

  ! Cutoff drag above ALT_CUTOFF, using BVFB as temporary work array.

  IF (icutoff==1) THEN
    CALL hines_exp(drag_u, bvfb, alt, alt_cutoff, iorder, il1, il2, lev1, &
      lev2, nlons, nlevs)
    CALL hines_exp(drag_v, bvfb, alt, alt_cutoff, iorder, il1, il2, lev1, &
      lev2, nlons, nlevs)
  END IF

  ! Print out various arrays for diagnostic purposes.

  IF (iprint==1) THEN
    ilprt1 = 15
    ilprt2 = 16
    CALL hines_print(flux_u, flux_v, drag_u, drag_v, alt, sigma_t, &
      sigma_alpha, v_alpha, m_alpha, 1, 1, 6, ilprt1, ilprt2, lev1, lev2, &
      naz, nlons, nlevs, nazmth)
  END IF

  ! If not calculating heating rate and diffusion coefficient then finished.

  IF (iheatcal/=1) RETURN

  ! Calculate vertical derivative of cutoff wavenumber (store
  ! in array V_ALPHA) using centered differences at interior gridpoints
  ! and one-sided differences at first and last levels.

  DO n = 1, naz
    DO l = lev1p, lev2m
      DO i = il1, il2
        v_alpha(i, l, n) = (m_alpha(i,l+1,n)-m_alpha(i,l-1,n))/ &
          (alt(i,l+1)-alt(i,l-1))
      END DO
    END DO
    DO i = il1, il2
      v_alpha(i, lev1, n) = (m_alpha(i,lev1p,n)-m_alpha(i,lev1,n))/ &
        (alt(i,lev1p)-alt(i,lev1))
    END DO
    DO i = il1, il2
      v_alpha(i, lev2, n) = (m_alpha(i,lev2,n)-m_alpha(i,lev2m,n))/ &
        (alt(i,lev2)-alt(i,lev2m))
    END DO
  END DO

  ! Heating rate and diffusion coefficient.

  CALL hines_heat(heat, diffco, m_alpha, v_alpha, ak_alpha, k_alpha, bvfreq, &
    density, densb, sigma_t, visc_mol, kstar, slope, f2, f3, f5, f6, naz, &
    il1, il2, lev1, lev2, nlons, nlevs, nazmth)

  ! Finished.

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE hines_extro0

SUBROUTINE hines_wavnum(m_alpha, sigma_alpha, sigsqh_alpha, sigma_t, &
    ak_alpha, v_alpha, visc_mol, density, densb, bvfreq, bvfb, rms_wind, &
    i_alpha, mmin_alpha, kstar, slope, f1, f2, f3, naz, levbot, levtop, il1, &
    il2, nlons, nlevs, nazmth, sigsqmcw, sigmatm, lorms, sigalpmc, f2mod)
  IMPLICIT NONE
  ! This routine calculates the cutoff vertical wavenumber and velocity
  ! variances on a longitude by altitude grid for the Hines' Doppler
  ! spread gravity wave drag parameterization scheme.
  ! NOTE: (1) only values of four or eight can be used for # azimuths (NAZ).
  ! (2) only values of 1.0, 1.5 or 2.0 can be used for slope (SLOPE).

  ! Aug. 10/95 - C. McLandress

  ! Output arguements:

  ! * M_ALPHA      = cutoff wavenumber at each azimuth (1/m).
  ! * SIGMA_ALPHA  = total rms wind in each azimuth (m/s).
  ! * SIGSQH_ALPHA = portion of wind variance from waves having wave
  ! *                normals in the alpha azimuth (m/s).
  ! * SIGMA_T      = total rms horizontal wind (m/s).
  ! * AK_ALPHA     = spectral amplitude factor at each azimuth
  ! *                (i.e.,{AjKj}) in m^4/s^2.

  ! Input arguements:

  ! * V_ALPHA  = wind component at each azimuth (m/s).
  ! * VISC_MOL = molecular viscosity (m^2/s)
  ! * DENSITY  = background density (kg/m^3).
  ! * DENSB    = background density at model bottom (kg/m^3).
  ! * BVFREQ   = background Brunt Vassala frequency (radians/sec).
  ! * BVFB     = background Brunt Vassala frequency at model bottom.
  ! * RMS_WIND = root mean square gravity wave wind at lowest level (m/s).
  ! * KSTAR    = typical gravity wave horizontal wavenumber (1/m).
  ! * SLOPE    = slope of incident vertical wavenumber spectrum
  ! *            (SLOPE = 1., 1.5 or 2.).
  ! * F1,F2,F3 = Hines's fudge factors.
  ! * NAZ      = actual number of horizontal azimuths used (4 or 8).
  ! * LEVBOT   = index of lowest vertical level.
  ! * LEVTOP   = index of highest vertical level
  ! *            (NOTE: if LEVTOP < LEVBOT then level index
  ! *             increases from top down).
  ! * IL1      = first longitudinal index to use (IL1 >= 1).
  ! * IL2      = last longitudinal index to use (IL1 <= IL2 <= NLONS).
  ! * NLONS    = number of longitudes.
  ! * NLEVS    = number of vertical levels.
  ! * NAZMTH   = azimuthal array dimension (NAZMTH >= NAZ).

  ! * LORMS       = .TRUE. for drag computation

  ! Input work arrays:

  ! * I_ALPHA    = Hines' integral at a single level.
  ! * MMIN_ALPHA = minimum value of cutoff wavenumber.

  INTEGER naz, levbot, levtop, il1, il2, nlons, nlevs, nazmth
  REAL slope, kstar(nlons), f1, f2, f3, f2mfac
  REAL m_alpha(nlons, nlevs, nazmth)
  REAL sigma_alpha(nlons, nlevs, nazmth)
  REAL sigalpmc(nlons, nlevs, nazmth)
  REAL sigsqh_alpha(nlons, nlevs, nazmth)
  REAL sigsqmcw(nlons, nlevs, nazmth)
  REAL sigma_t(nlons, nlevs)
  REAL sigmatm(nlons, nlevs)
  REAL ak_alpha(nlons, nazmth)
  REAL v_alpha(nlons, nlevs, nazmth)
  REAL visc_mol(nlons, nlevs)
  REAL f2mod(nlons, nlevs)
  REAL density(nlons, nlevs), densb(nlons)
  REAL bvfreq(nlons, nlevs), bvfb(nlons), rms_wind(nlons)
  REAL i_alpha(nlons, nazmth), mmin_alpha(nlons, nazmth)

  LOGICAL lorms(nlons)

  ! Internal variables.

  INTEGER i, l, n, lstart, lend, lincr, lbelow
  REAL m_sub_m_turb, m_sub_m_mol, m_trial
  REAL visc, visc_min, azfac, sp1

  ! c      REAL  N_OVER_M(1000), SIGFAC(1000)

  REAL n_over_m(nlons), sigfac(nlons)
  DATA visc_min/1.E-10/
  ! -----------------------------------------------------------------------


  ! PRINT *,'IN HINES_WAVNUM'
  sp1 = slope + 1.

  ! Indices of levels to process.

  IF (levbot>levtop) THEN
    lstart = levbot - 1
    lend = levtop
    lincr = -1
  ELSE
    WRITE (6, 1)
1   FORMAT (2X, ' error: IORDER NOT ONE! ')
  END IF

  ! Use horizontal isotropy to calculate azimuthal variances at bottom level.

  azfac = 1./real(naz)
  DO n = 1, naz
    DO i = il1, il2
      sigsqh_alpha(i, levbot, n) = azfac*rms_wind(i)**2
    END DO
  END DO

  ! Velocity variances at bottom level.

  CALL hines_sigma(sigma_t, sigma_alpha, sigsqh_alpha, naz, levbot, il1, il2, &
    nlons, nlevs, nazmth)

  CALL hines_sigma(sigmatm, sigalpmc, sigsqmcw, naz, levbot, il1, il2, nlons, &
    nlevs, nazmth)

  ! Calculate cutoff wavenumber and spectral amplitude factor
  ! at bottom level where it is assumed that background winds vanish
  ! and also initialize minimum value of cutoff wavnumber.

  DO n = 1, naz
    DO i = il1, il2
      IF (lorms(i)) THEN
        m_alpha(i, levbot, n) = bvfb(i)/(f1*sigma_alpha(i,levbot,n)+f2* &
          sigma_t(i,levbot))
        ak_alpha(i, n) = sigsqh_alpha(i, levbot, n)/ &
          (m_alpha(i,levbot,n)**sp1/sp1)
        mmin_alpha(i, n) = m_alpha(i, levbot, n)
      END IF
    END DO
  END DO

  ! Calculate quantities from the bottom upwards,
  ! starting one level above bottom.

  DO l = lstart, lend, lincr

    ! Level beneath present level.

    lbelow = l - lincr

    ! Calculate N/m_M where m_M is maximum permissible value of the vertical
    ! wavenumber (i.e., m > m_M are obliterated) and N is buoyancy frequency.
    ! m_M is taken as the smaller of the instability-induced
    ! wavenumber (M_SUB_M_TURB) and that imposed by molecular viscosity
    ! (M_SUB_M_MOL). Since variance at this level is not yet known
    ! use value at level below.

    DO i = il1, il2
      IF (lorms(i)) THEN

        f2mfac = sigmatm(i, lbelow)**2
        f2mod(i, lbelow) = 1. + 2.*f2mfac/(f2mfac+sigma_t(i,lbelow)**2)

        visc = amax1(visc_mol(i,l), visc_min)
        m_sub_m_turb = bvfreq(i, l)/(f2*f2mod(i,lbelow)*sigma_t(i,lbelow))
        m_sub_m_mol = (bvfreq(i,l)*kstar(i)/visc)**0.33333333/f3
        IF (m_sub_m_turb<m_sub_m_mol) THEN
          n_over_m(i) = f2*f2mod(i, lbelow)*sigma_t(i, lbelow)
        ELSE
          n_over_m(i) = bvfreq(i, l)/m_sub_m_mol
        END IF
      END IF
    END DO

    ! Calculate cutoff wavenumber at this level.

    DO n = 1, naz
      DO i = il1, il2
        IF (lorms(i)) THEN

          ! Calculate trial value (since variance at this level is not yet
          ! known
          ! use value at level below). If trial value is negative or if it
          ! exceeds
          ! minimum value (not permitted) then set it to minimum value.

          m_trial = bvfb(i)/(f1*(sigma_alpha(i,lbelow,n)+sigalpmc(i,lbelow, &
            n))+n_over_m(i)+v_alpha(i,l,n))
          IF (m_trial<=0. .OR. m_trial>mmin_alpha(i,n)) THEN
            m_trial = mmin_alpha(i, n)
          END IF
          m_alpha(i, l, n) = m_trial

          ! Reset minimum value of cutoff wavenumber if necessary.

          IF (m_alpha(i,l,n)<mmin_alpha(i,n)) THEN
            mmin_alpha(i, n) = m_alpha(i, l, n)
          END IF

        END IF
      END DO
    END DO

    ! Calculate the Hines integral at this level.

    CALL hines_intgrl(i_alpha, v_alpha, m_alpha, bvfb, slope, naz, l, il1, &
      il2, nlons, nlevs, nazmth, lorms)


    ! Calculate the velocity variances at this level.

    DO i = il1, il2
      sigfac(i) = densb(i)/density(i, l)*bvfreq(i, l)/bvfb(i)
    END DO
    DO n = 1, naz
      DO i = il1, il2
        sigsqh_alpha(i, l, n) = sigfac(i)*ak_alpha(i, n)*i_alpha(i, n)
      END DO
    END DO
    CALL hines_sigma(sigma_t, sigma_alpha, sigsqh_alpha, naz, l, il1, il2, &
      nlons, nlevs, nazmth)

    CALL hines_sigma(sigmatm, sigalpmc, sigsqmcw, naz, l, il1, il2, nlons, &
      nlevs, nazmth)

    ! End of level loop.

  END DO

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE hines_wavnum

SUBROUTINE hines_wind(v_alpha, vel_u, vel_v, naz, il1, il2, lev1, lev2, &
    nlons, nlevs, nazmth)
    IMPLICIT NONE
  ! This routine calculates the azimuthal horizontal background wind
  ! components
  ! on a longitude by altitude grid for the case of 4 or 8 azimuths for
  ! the Hines' Doppler spread GWD parameterization scheme.

  ! Aug. 7/95 - C. McLandress

  ! Output arguement:

  ! * V_ALPHA   = background wind component at each azimuth (m/s).
  ! *             (note: first azimuth is in eastward direction
  ! *              and rotate in counterclockwise direction.)

  ! Input arguements:

  ! * VEL_U     = background zonal wind component (m/s).
  ! * VEL_V     = background meridional wind component (m/s).
  ! * NAZ       = actual number of horizontal azimuths used (must be 4 or 8).
  ! * IL1       = first longitudinal index to use (IL1 >= 1).
  ! * IL2       = last longitudinal index to use (IL1 <= IL2 <= NLONS).
  ! * LEV1      = first altitude level to use (LEV1 >=1).
  ! * LEV2      = last altitude level to use (LEV1 < LEV2 <= NLEVS).
  ! * NLONS     = number of longitudes.
  ! * NLEVS     = number of vertical levels.
  ! * NAZMTH    = azimuthal array dimension (NAZMTH >= NAZ).

  ! Constants in DATA statements.

  ! * COS45 = cosine of 45 degrees.
  ! * UMIN  = minimum allowable value for zonal or meridional
  ! *         wind component (m/s).

  ! Subroutine arguements.

  INTEGER naz, il1, il2, lev1, lev2
  INTEGER nlons, nlevs, nazmth
  REAL v_alpha(nlons, nlevs, nazmth)
  REAL vel_u(nlons, nlevs), vel_v(nlons, nlevs)

  ! Internal variables.

  INTEGER i, l
  REAL u, v, cos45, umin

  DATA cos45/0.7071068/
  DATA umin/0.001/
  ! -----------------------------------------------------------------------

  ! Case with 4 azimuths.


  ! PRINT *,'IN HINES_WIND'
  IF (naz==4) THEN
    DO l = lev1, lev2
      DO i = il1, il2
        u = vel_u(i, l)
        v = vel_v(i, l)
        IF (abs(u)<umin) u = umin
        IF (abs(v)<umin) v = umin
        v_alpha(i, l, 1) = u
        v_alpha(i, l, 2) = v
        v_alpha(i, l, 3) = -u
        v_alpha(i, l, 4) = -v
      END DO
    END DO
  END IF

  ! Case with 8 azimuths.

  IF (naz==8) THEN
    DO l = lev1, lev2
      DO i = il1, il2
        u = vel_u(i, l)
        v = vel_v(i, l)
        IF (abs(u)<umin) u = umin
        IF (abs(v)<umin) v = umin
        v_alpha(i, l, 1) = u
        v_alpha(i, l, 2) = cos45*(v+u)
        v_alpha(i, l, 3) = v
        v_alpha(i, l, 4) = cos45*(v-u)
        v_alpha(i, l, 5) = -u
        v_alpha(i, l, 6) = -v_alpha(i, l, 2)
        v_alpha(i, l, 7) = -v
        v_alpha(i, l, 8) = -v_alpha(i, l, 4)
      END DO
    END DO
  END IF

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE hines_wind

SUBROUTINE hines_flux(flux_u, flux_v, drag_u, drag_v, alt, density, densb, &
    m_alpha, ak_alpha, k_alpha, slope, naz, il1, il2, lev1, lev2, nlons, &
    nlevs, nazmth, lorms)
    IMPLICIT NONE
  ! Calculate zonal and meridional components of the vertical flux
  ! of horizontal momentum and corresponding wave drag (force per unit mass)
  ! on a longitude by altitude grid for the Hines' Doppler spread
  ! GWD parameterization scheme.
  ! NOTE: only 4 or 8 azimuths can be used.

  ! Aug. 6/95 - C. McLandress

  ! Output arguements:

  ! * FLUX_U = zonal component of vertical momentum flux (Pascals)
  ! * FLUX_V = meridional component of vertical momentum flux (Pascals)
  ! * DRAG_U = zonal component of drag (m/s^2).
  ! * DRAG_V = meridional component of drag (m/s^2).

  ! Input arguements:

  ! * ALT       = altitudes (m).
  ! * DENSITY   = background density (kg/m^3).
  ! * DENSB     = background density at bottom level (kg/m^3).
  ! * M_ALPHA   = cutoff vertical wavenumber (1/m).
  ! * AK_ALPHA  = spectral amplitude factor (i.e., {AjKj} in m^4/s^2).
  ! * K_ALPHA   = horizontal wavenumber (1/m).
  ! * SLOPE     = slope of incident vertical wavenumber spectrum.
  ! * NAZ       = actual number of horizontal azimuths used (must be 4 or 8).
  ! * IL1       = first longitudinal index to use (IL1 >= 1).
  ! * IL2       = last longitudinal index to use (IL1 <= IL2 <= NLONS).
  ! * LEV1      = first altitude level to use (LEV1 >=1).
  ! * LEV2      = last altitude level to use (LEV1 < LEV2 <= NLEVS).
  ! * NLONS     = number of longitudes.
  ! * NLEVS     = number of vertical levels.
  ! * NAZMTH    = azimuthal array dimension (NAZMTH >= NAZ).

  ! * LORMS       = .TRUE. for drag computation

  ! Constant in DATA statement.

  ! * COS45 = cosine of 45 degrees.

  ! Subroutine arguements.

  INTEGER naz, il1, il2, lev1, lev2
  INTEGER nlons, nlevs, nazmth
  REAL slope
  REAL flux_u(nlons, nlevs), flux_v(nlons, nlevs)
  REAL drag_u(nlons, nlevs), drag_v(nlons, nlevs)
  REAL alt(nlons, nlevs), density(nlons, nlevs), densb(nlons)
  REAL m_alpha(nlons, nlevs, nazmth)
  REAL ak_alpha(nlons, nazmth), k_alpha(nlons, nazmth)

  LOGICAL lorms(nlons)

  ! Internal variables.

  INTEGER i, l, lev1p, lev2m, lev2p
  REAL cos45, prod2, prod4, prod6, prod8, dendz, dendz2
  DATA cos45/0.7071068/
  ! -----------------------------------------------------------------------

  lev1p = lev1 + 1
  lev2m = lev2 - 1
  lev2p = lev2 + 1

  ! Sum over azimuths for case where SLOPE = 1.

  IF (slope==1.) THEN

    ! Case with 4 azimuths.

    IF (naz==4) THEN
      DO l = lev1, lev2
        DO i = il1, il2
          flux_u(i, l) = ak_alpha(i, 1)*k_alpha(i, 1)*m_alpha(i, l, 1) - &
            ak_alpha(i, 3)*k_alpha(i, 3)*m_alpha(i, l, 3)
          flux_v(i, l) = ak_alpha(i, 2)*k_alpha(i, 2)*m_alpha(i, l, 2) - &
            ak_alpha(i, 4)*k_alpha(i, 4)*m_alpha(i, l, 4)
        END DO
      END DO
    END IF

    ! Case with 8 azimuths.

    IF (naz==8) THEN
      DO l = lev1, lev2
        DO i = il1, il2
          prod2 = ak_alpha(i, 2)*k_alpha(i, 2)*m_alpha(i, l, 2)
          prod4 = ak_alpha(i, 4)*k_alpha(i, 4)*m_alpha(i, l, 4)
          prod6 = ak_alpha(i, 6)*k_alpha(i, 6)*m_alpha(i, l, 6)
          prod8 = ak_alpha(i, 8)*k_alpha(i, 8)*m_alpha(i, l, 8)
          flux_u(i, l) = ak_alpha(i, 1)*k_alpha(i, 1)*m_alpha(i, l, 1) - &
            ak_alpha(i, 5)*k_alpha(i, 5)*m_alpha(i, l, 5) + &
            cos45*(prod2-prod4-prod6+prod8)
          flux_v(i, l) = ak_alpha(i, 3)*k_alpha(i, 3)*m_alpha(i, l, 3) - &
            ak_alpha(i, 7)*k_alpha(i, 7)*m_alpha(i, l, 7) + &
            cos45*(prod2+prod4-prod6-prod8)
        END DO
      END DO
    END IF

  END IF

  ! Sum over azimuths for case where SLOPE not equal to 1.

  IF (slope/=1.) THEN

    ! Case with 4 azimuths.

    IF (naz==4) THEN
      DO l = lev1, lev2
        DO i = il1, il2
          flux_u(i, l) = ak_alpha(i, 1)*k_alpha(i, 1)* &
            m_alpha(i, l, 1)**slope - ak_alpha(i, 3)*k_alpha(i, 3)*m_alpha(i, &
            l, 3)**slope
          flux_v(i, l) = ak_alpha(i, 2)*k_alpha(i, 2)* &
            m_alpha(i, l, 2)**slope - ak_alpha(i, 4)*k_alpha(i, 4)*m_alpha(i, &
            l, 4)**slope
        END DO
      END DO
    END IF

    ! Case with 8 azimuths.

    IF (naz==8) THEN
      DO l = lev1, lev2
        DO i = il1, il2
          prod2 = ak_alpha(i, 2)*k_alpha(i, 2)*m_alpha(i, l, 2)**slope
          prod4 = ak_alpha(i, 4)*k_alpha(i, 4)*m_alpha(i, l, 4)**slope
          prod6 = ak_alpha(i, 6)*k_alpha(i, 6)*m_alpha(i, l, 6)**slope
          prod8 = ak_alpha(i, 8)*k_alpha(i, 8)*m_alpha(i, l, 8)**slope
          flux_u(i, l) = ak_alpha(i, 1)*k_alpha(i, 1)* &
            m_alpha(i, l, 1)**slope - ak_alpha(i, 5)*k_alpha(i, 5)*m_alpha(i, &
            l, 5)**slope + cos45*(prod2-prod4-prod6+prod8)
          flux_v(i, l) = ak_alpha(i, 3)*k_alpha(i, 3)* &
            m_alpha(i, l, 3)**slope - ak_alpha(i, 7)*k_alpha(i, 7)*m_alpha(i, &
            l, 7)**slope + cos45*(prod2+prod4-prod6-prod8)
        END DO
      END DO
    END IF

  END IF

  ! Calculate flux from sum.

  DO l = lev1, lev2
    DO i = il1, il2
      flux_u(i, l) = flux_u(i, l)*densb(i)/slope
      flux_v(i, l) = flux_v(i, l)*densb(i)/slope
    END DO
  END DO

  ! Calculate drag at intermediate levels using centered differences

  DO l = lev1p, lev2m
    DO i = il1, il2
      IF (lorms(i)) THEN
        ! cc       DENDZ2 = DENSITY(I,L) * ( ALT(I,L+1) - ALT(I,L-1) )
        dendz2 = density(i, l)*(alt(i,l-1)-alt(i,l))
        ! cc       DRAG_U(I,L) = - ( FLUX_U(I,L+1) - FLUX_U(I,L-1) ) / DENDZ2
        drag_u(i, l) = -(flux_u(i,l-1)-flux_u(i,l))/dendz2
        ! cc       DRAG_V(I,L) = - ( FLUX_V(I,L+1) - FLUX_V(I,L-1) ) / DENDZ2
        drag_v(i, l) = -(flux_v(i,l-1)-flux_v(i,l))/dendz2

      END IF
    END DO
  END DO

  ! Drag at first and last levels using one-side differences.

  DO i = il1, il2
    IF (lorms(i)) THEN
      dendz = density(i, lev1)*(alt(i,lev1)-alt(i,lev1p))
      drag_u(i, lev1) = flux_u(i, lev1)/dendz
      drag_v(i, lev1) = flux_v(i, lev1)/dendz
    END IF
  END DO
  DO i = il1, il2
    IF (lorms(i)) THEN
      dendz = density(i, lev2)*(alt(i,lev2m)-alt(i,lev2))
      drag_u(i, lev2) = -(flux_u(i,lev2m)-flux_u(i,lev2))/dendz
      drag_v(i, lev2) = -(flux_v(i,lev2m)-flux_v(i,lev2))/dendz
    END IF
  END DO
  IF (nlevs>lev2) THEN
    DO i = il1, il2
      IF (lorms(i)) THEN
        dendz = density(i, lev2p)*(alt(i,lev2)-alt(i,lev2p))
        drag_u(i, lev2p) = -flux_u(i, lev2)/dendz
        drag_v(i, lev2p) = -flux_v(i, lev2)/dendz
      END IF
    END DO
  END IF

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE hines_flux

SUBROUTINE hines_heat(heat, diffco, m_alpha, dmdz_alpha, ak_alpha, k_alpha, &
    bvfreq, density, densb, sigma_t, visc_mol, kstar, slope, f2, f3, f5, f6, &
    naz, il1, il2, lev1, lev2, nlons, nlevs, nazmth)
  IMPLICIT NONE
  ! This routine calculates the gravity wave induced heating and
  ! diffusion coefficient on a longitude by altitude grid for
  ! the Hines' Doppler spread gravity wave drag parameterization scheme.

  ! Aug. 6/95 - C. McLandress

  ! Output arguements:

  ! * HEAT   = gravity wave heating (K/sec).
  ! * DIFFCO = diffusion coefficient (m^2/sec)

  ! Input arguements:

  ! * M_ALPHA     = cutoff vertical wavenumber (1/m).
  ! * DMDZ_ALPHA  = vertical derivative of cutoff wavenumber.
  ! * AK_ALPHA    = spectral amplitude factor of each azimuth
  ! (i.e., {AjKj} in m^4/s^2).
  ! * K_ALPHA     = horizontal wavenumber of each azimuth (1/m).
  ! * BVFREQ      = background Brunt Vassala frequency (rad/sec).
  ! * DENSITY     = background density (kg/m^3).
  ! * DENSB       = background density at bottom level (kg/m^3).
  ! * SIGMA_T     = total rms horizontal wind (m/s).
  ! * VISC_MOL    = molecular viscosity (m^2/s).
  ! * KSTAR       = typical gravity wave horizontal wavenumber (1/m).
  ! * SLOPE       = slope of incident vertical wavenumber spectrum.
  ! * F2,F3,F5,F6 = Hines's fudge factors.
  ! * NAZ         = actual number of horizontal azimuths used.
  ! * IL1         = first longitudinal index to use (IL1 >= 1).
  ! * IL2         = last longitudinal index to use (IL1 <= IL2 <= NLONS).
  ! * LEV1        = first altitude level to use (LEV1 >=1).
  ! * LEV2        = last altitude level to use (LEV1 < LEV2 <= NLEVS).
  ! * NLONS       = number of longitudes.
  ! * NLEVS       = number of vertical levels.
  ! * NAZMTH      = azimuthal array dimension (NAZMTH >= NAZ).

  INTEGER naz, il1, il2, lev1, lev2, nlons, nlevs, nazmth
  REAL kstar(nlons), slope, f2, f3, f5, f6
  REAL heat(nlons, nlevs), diffco(nlons, nlevs)
  REAL m_alpha(nlons, nlevs, nazmth), dmdz_alpha(nlons, nlevs, nazmth)
  REAL ak_alpha(nlons, nazmth), k_alpha(nlons, nazmth)
  REAL bvfreq(nlons, nlevs), density(nlons, nlevs), densb(nlons)
  REAL sigma_t(nlons, nlevs), visc_mol(nlons, nlevs)

  ! Internal variables.

  INTEGER i, l, n
  REAL m_sub_m_turb, m_sub_m_mol, m_sub_m, heatng
  REAL visc, visc_min, cpgas, sm1

  ! specific heat at constant pressure

  DATA cpgas/1004./

  ! minimum permissible viscosity

  DATA visc_min/1.E-10/
  ! -----------------------------------------------------------------------

  ! Initialize heating array.

  DO l = 1, nlevs
    DO i = 1, nlons
      heat(i, l) = 0.
    END DO
  END DO

  ! Perform sum over azimuths for case where SLOPE = 1.

  IF (slope==1.) THEN
    DO n = 1, naz
      DO l = lev1, lev2
        DO i = il1, il2
          heat(i, l) = heat(i, l) + ak_alpha(i, n)*k_alpha(i, n)*dmdz_alpha(i &
            , l, n)
        END DO
      END DO
    END DO
  END IF

  ! Perform sum over azimuths for case where SLOPE not 1.

  IF (slope/=1.) THEN
    sm1 = slope - 1.
    DO n = 1, naz
      DO l = lev1, lev2
        DO i = il1, il2
          heat(i, l) = heat(i, l) + ak_alpha(i, n)*k_alpha(i, n)*m_alpha(i, l &
            , n)**sm1*dmdz_alpha(i, l, n)
        END DO
      END DO
    END DO
  END IF

  ! Heating and diffusion.

  DO l = lev1, lev2
    DO i = il1, il2

      ! Maximum permissible value of cutoff wavenumber is the smaller
      ! of the instability-induced wavenumber (M_SUB_M_TURB) and
      ! that imposed by molecular viscosity (M_SUB_M_MOL).

      visc = amax1(visc_mol(i,l), visc_min)
      m_sub_m_turb = bvfreq(i, l)/(f2*sigma_t(i,l))
      m_sub_m_mol = (bvfreq(i,l)*kstar(i)/visc)**0.33333333/f3
      m_sub_m = amin1(m_sub_m_turb, m_sub_m_mol)

      heatng = -heat(i, l)*f5*bvfreq(i, l)/m_sub_m*densb(i)/density(i, l)
      diffco(i, l) = f6*heatng**0.33333333/m_sub_m**1.33333333
      heat(i, l) = heatng/cpgas

    END DO
  END DO

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE hines_heat

SUBROUTINE hines_sigma(sigma_t, sigma_alpha, sigsqh_alpha, naz, lev, il1, &
    il2, nlons, nlevs, nazmth)
  IMPLICIT NONE
  ! This routine calculates the total rms and azimuthal rms horizontal
  ! velocities at a given level on a longitude by altitude grid for
  ! the Hines' Doppler spread GWD parameterization scheme.
  ! NOTE: only four or eight azimuths can be used.

  ! Aug. 7/95 - C. McLandress

  ! Output arguements:

  ! * SIGMA_T      = total rms horizontal wind (m/s).
  ! * SIGMA_ALPHA  = total rms wind in each azimuth (m/s).

  ! Input arguements:

  ! * SIGSQH_ALPHA = portion of wind variance from waves having wave
  ! *                normals in the alpha azimuth (m/s).
  ! * NAZ       = actual number of horizontal azimuths used (must be 4 or 8).
  ! * LEV       = altitude level to process.
  ! * IL1       = first longitudinal index to use (IL1 >= 1).
  ! * IL2       = last longitudinal index to use (IL1 <= IL2 <= NLONS).
  ! * NLONS     = number of longitudes.
  ! * NLEVS     = number of vertical levels.
  ! * NAZMTH    = azimuthal array dimension (NAZMTH >= NAZ).

  ! Subroutine arguements.

  INTEGER lev, naz, il1, il2
  INTEGER nlons, nlevs, nazmth
  REAL sigma_t(nlons, nlevs)
  REAL sigma_alpha(nlons, nlevs, nazmth)
  REAL sigsqh_alpha(nlons, nlevs, nazmth)

  ! Internal variables.

  INTEGER i, n
  REAL sum_even, sum_odd
  ! -----------------------------------------------------------------------

  ! Calculate azimuthal rms velocity for the 4 azimuth case.

  IF (naz==4) THEN
    DO i = il1, il2
      sigma_alpha(i, lev, 1) = sqrt(sigsqh_alpha(i,lev,1)+sigsqh_alpha(i,lev, &
        3))
      sigma_alpha(i, lev, 2) = sqrt(sigsqh_alpha(i,lev,2)+sigsqh_alpha(i,lev, &
        4))
      sigma_alpha(i, lev, 3) = sigma_alpha(i, lev, 1)
      sigma_alpha(i, lev, 4) = sigma_alpha(i, lev, 2)
    END DO
  END IF

  ! Calculate azimuthal rms velocity for the 8 azimuth case.

  IF (naz==8) THEN
    DO i = il1, il2
      sum_odd = (sigsqh_alpha(i,lev,1)+sigsqh_alpha(i,lev,3)+ &
        sigsqh_alpha(i,lev,5)+sigsqh_alpha(i,lev,7))/2.
      sum_even = (sigsqh_alpha(i,lev,2)+sigsqh_alpha(i,lev,4)+ &
        sigsqh_alpha(i,lev,6)+sigsqh_alpha(i,lev,8))/2.
      sigma_alpha(i, lev, 1) = sqrt(sigsqh_alpha(i,lev,1)+sigsqh_alpha(i,lev, &
        5)+sum_even)
      sigma_alpha(i, lev, 2) = sqrt(sigsqh_alpha(i,lev,2)+sigsqh_alpha(i,lev, &
        6)+sum_odd)
      sigma_alpha(i, lev, 3) = sqrt(sigsqh_alpha(i,lev,3)+sigsqh_alpha(i,lev, &
        7)+sum_even)
      sigma_alpha(i, lev, 4) = sqrt(sigsqh_alpha(i,lev,4)+sigsqh_alpha(i,lev, &
        8)+sum_odd)
      sigma_alpha(i, lev, 5) = sigma_alpha(i, lev, 1)
      sigma_alpha(i, lev, 6) = sigma_alpha(i, lev, 2)
      sigma_alpha(i, lev, 7) = sigma_alpha(i, lev, 3)
      sigma_alpha(i, lev, 8) = sigma_alpha(i, lev, 4)
    END DO
  END IF

  ! Calculate total rms velocity.

  DO i = il1, il2
    sigma_t(i, lev) = 0.
  END DO
  DO n = 1, naz
    DO i = il1, il2
      sigma_t(i, lev) = sigma_t(i, lev) + sigsqh_alpha(i, lev, n)
    END DO
  END DO
  DO i = il1, il2
    sigma_t(i, lev) = sqrt(sigma_t(i,lev))
  END DO

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE hines_sigma

SUBROUTINE hines_intgrl(i_alpha, v_alpha, m_alpha, bvfb, slope, naz, lev, &
    il1, il2, nlons, nlevs, nazmth, lorms)
  IMPLICIT NONE
  ! This routine calculates the vertical wavenumber integral
  ! for a single vertical level at each azimuth on a longitude grid
  ! for the Hines' Doppler spread GWD parameterization scheme.
  ! NOTE: (1) only spectral slopes of 1, 1.5 or 2 are permitted.
  ! (2) the integral is written in terms of the product QM
  ! which by construction is always less than 1. Series
  ! solutions are used for small |QM| and analytical solutions
  ! for remaining values.

  ! Aug. 8/95 - C. McLandress

  ! Output arguement:

  ! * I_ALPHA = Hines' integral.

  ! Input arguements:

  ! * V_ALPHA = azimuthal wind component (m/s).
  ! * M_ALPHA = azimuthal cutoff vertical wavenumber (1/m).
  ! * BVFB    = background Brunt Vassala frequency at model bottom.
  ! * SLOPE   = slope of initial vertical wavenumber spectrum
  ! *           (must use SLOPE = 1., 1.5 or 2.)
  ! * NAZ     = actual number of horizontal azimuths used.
  ! * LEV     = altitude level to process.
  ! * IL1     = first longitudinal index to use (IL1 >= 1).
  ! * IL2     = last longitudinal index to use (IL1 <= IL2 <= NLONS).
  ! * NLONS   = number of longitudes.
  ! * NLEVS   = number of vertical levels.
  ! * NAZMTH  = azimuthal array dimension (NAZMTH >= NAZ).

  ! * LORMS       = .TRUE. for drag computation

  ! Constants in DATA statements:

  ! * QMIN = minimum value of Q_ALPHA (avoids indeterminant form of integral)
  ! * QM_MIN = minimum value of Q_ALPHA * M_ALPHA (used to avoid numerical
  ! *          problems).

  INTEGER lev, naz, il1, il2, nlons, nlevs, nazmth
  REAL i_alpha(nlons, nazmth)
  REAL v_alpha(nlons, nlevs, nazmth)
  REAL m_alpha(nlons, nlevs, nazmth)
  REAL bvfb(nlons), slope

  LOGICAL lorms(nlons)

  ! Internal variables.

  INTEGER i, n
  REAL q_alpha, qm, sqrtqm, q_min, qm_min

  DATA q_min/1.0/, qm_min/0.01/
  ! -----------------------------------------------------------------------

  ! For integer value SLOPE = 1.

  IF (slope==1.) THEN

    DO n = 1, naz
      DO i = il1, il2
        IF (lorms(i)) THEN

          q_alpha = v_alpha(i, lev, n)/bvfb(i)
          qm = q_alpha*m_alpha(i, lev, n)

          ! If |QM| is small then use first 4 terms series of Taylor series
          ! expansion of integral in order to avoid indeterminate form of
          ! integral,
          ! otherwise use analytical form of integral.

          IF (abs(q_alpha)<q_min .OR. abs(qm)<qm_min) THEN
            IF (q_alpha==0.) THEN
              i_alpha(i, n) = m_alpha(i, lev, n)**2/2.
            ELSE
              i_alpha(i, n) = (qm**2/2.+qm**3/3.+qm**4/4.+qm**5/5.)/ &
                q_alpha**2
            END IF
          ELSE
            i_alpha(i, n) = -(alog(1.-qm)+qm)/q_alpha**2
          END IF

        END IF
      END DO
    END DO

  END IF

  ! For integer value SLOPE = 2.

  IF (slope==2.) THEN

    DO n = 1, naz
      DO i = il1, il2
        IF (lorms(i)) THEN

          q_alpha = v_alpha(i, lev, n)/bvfb(i)
          qm = q_alpha*m_alpha(i, lev, n)

          ! If |QM| is small then use first 4 terms series of Taylor series
          ! expansion of integral in order to avoid indeterminate form of
          ! integral,
          ! otherwise use analytical form of integral.

          IF (abs(q_alpha)<q_min .OR. abs(qm)<qm_min) THEN
            IF (q_alpha==0.) THEN
              i_alpha(i, n) = m_alpha(i, lev, n)**3/3.
            ELSE
              i_alpha(i, n) = (qm**3/3.+qm**4/4.+qm**5/5.+qm**6/6.)/ &
                q_alpha**3
            END IF
          ELSE
            i_alpha(i, n) = -(alog(1.-qm)+qm+qm**2/2.)/q_alpha**3
          END IF

        END IF
      END DO
    END DO

  END IF

  ! For real value SLOPE = 1.5

  IF (slope==1.5) THEN

    DO n = 1, naz
      DO i = il1, il2
        IF (lorms(i)) THEN

          q_alpha = v_alpha(i, lev, n)/bvfb(i)
          qm = q_alpha*m_alpha(i, lev, n)

          ! If |QM| is small then use first 4 terms series of Taylor series
          ! expansion of integral in order to avoid indeterminate form of
          ! integral,
          ! otherwise use analytical form of integral.

          IF (abs(q_alpha)<q_min .OR. abs(qm)<qm_min) THEN
            IF (q_alpha==0.) THEN
              i_alpha(i, n) = m_alpha(i, lev, n)**2.5/2.5
            ELSE
              i_alpha(i, n) = (qm/2.5+qm**2/3.5+qm**3/4.5+qm**4/5.5)* &
                m_alpha(i, lev, n)**1.5/q_alpha
            END IF
          ELSE
            qm = abs(qm)
            sqrtqm = sqrt(qm)
            IF (q_alpha>=0.) THEN
              i_alpha(i, n) = (alog((1.+sqrtqm)/(1.-sqrtqm))-2.*sqrtqm*(1.+qm &
                /3.))/q_alpha**2.5
            ELSE
              i_alpha(i, n) = 2.*(atan(sqrtqm)+sqrtqm*(qm/3.-1.))/ &
                abs(q_alpha)**2.5
            END IF
          END IF

        END IF
      END DO
    END DO

  END IF

  ! If integral is negative (which in principal should not happen) then
  ! print a message and some info since execution will abort when calculating
  ! the variances.

  ! DO 80 N = 1,NAZ
  ! DO 70 I = IL1,IL2
  ! IF (I_ALPHA(I,N).LT.0.)  THEN
  ! WRITE (6,*)
  ! WRITE (6,*) '******************************'
  ! WRITE (6,*) 'Hines integral I_ALPHA < 0 '
  ! WRITE (6,*) '  longitude I=',I
  ! WRITE (6,*) '  azimuth   N=',N
  ! WRITE (6,*) '  level   LEV=',LEV
  ! WRITE (6,*) '  I_ALPHA =',I_ALPHA(I,N)
  ! WRITE (6,*) '  V_ALPHA =',V_ALPHA(I,LEV,N)
  ! WRITE (6,*) '  M_ALPHA =',M_ALPHA(I,LEV,N)
  ! WRITE (6,*) '  Q_ALPHA =',V_ALPHA(I,LEV,N) / BVFB(I)
  ! WRITE (6,*) '  QM      =',V_ALPHA(I,LEV,N) / BVFB(I)
  ! ^                                * M_ALPHA(I,LEV,N)
  ! WRITE (6,*) '******************************'
  ! END IF
  ! 70     CONTINUE
  ! 80   CONTINUE

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE hines_intgrl

SUBROUTINE hines_setup(naz, slope, f1, f2, f3, f5, f6, kstar, icutoff, &
    alt_cutoff, smco, nsmax, iheatcal, k_alpha, ierror, nmessg, nlons, &
    nazmth, coslat)
  IMPLICIT NONE
  ! This routine specifies various parameters needed for the
  ! the Hines' Doppler spread gravity wave drag parameterization scheme.

  ! Aug. 8/95 - C. McLandress

  ! Output arguements:

  ! * NAZ        = actual number of horizontal azimuths used
  ! *              (code set up presently for only NAZ = 4 or 8).
  ! * SLOPE      = slope of incident vertical wavenumber spectrum
  ! *              (code set up presently for SLOPE 1., 1.5 or 2.).
  ! * F1         = "fudge factor" used in calculation of trial value of
  ! *              azimuthal cutoff wavenumber M_ALPHA (1.2 <= F1 <= 1.9).
  ! * F2         = "fudge factor" used in calculation of maximum
  ! *              permissible instabiliy-induced cutoff wavenumber
  ! *              M_SUB_M_TURB (0.1 <= F2 <= 1.4).
  ! * F3         = "fudge factor" used in calculation of maximum
  ! *              permissible molecular viscosity-induced cutoff wavenumber
  ! *              M_SUB_M_MOL (0.1 <= F2 <= 1.4).
  ! * F5         = "fudge factor" used in calculation of heating rate
  ! *              (1 <= F5 <= 3).
  ! * F6         = "fudge factor" used in calculation of turbulent
  ! *              diffusivity coefficient.
  ! * KSTAR      = typical gravity wave horizontal wavenumber (1/m)
  ! *              used in calculation of M_SUB_M_TURB.
  ! * ICUTOFF    = 1 to exponentially damp off GWD, heating and diffusion
  ! *              arrays above ALT_CUTOFF; otherwise arrays not modified.
  ! * ALT_CUTOFF = altitude in meters above which exponential decay applied.
  ! * SMCO       = smoother used to smooth cutoff vertical wavenumbers
  ! *              and total rms winds before calculating drag or heating.
  ! *              (==> a 1:SMCO:1 stencil used; SMCO >= 1.).
  ! * NSMAX      = number of times smoother applied ( >= 1),
  ! *            = 0 means no smoothing performed.
  ! * IHEATCAL   = 1 to calculate heating rates and diffusion coefficient.
  ! *            = 0 means only drag and flux calculated.
  ! * K_ALPHA    = horizontal wavenumber of each azimuth (1/m) which
  ! *              is set here to KSTAR.
  ! * IERROR     = error flag.
  ! *            = 0 no errors.
  ! *            = 10 ==> NAZ > NAZMTH
  ! *            = 20 ==> invalid number of azimuths (NAZ must be 4 or 8).
  ! *            = 30 ==> invalid slope (SLOPE must be 1., 1.5 or 2.).
  ! *            = 40 ==> invalid smoother (SMCO must be >= 1.)

  ! Input arguements:

  ! * NMESSG  = output unit number where messages to be printed.
  ! * NLONS   = number of longitudes.
  ! * NAZMTH  = azimuthal array dimension (NAZMTH >= NAZ).

  INTEGER naz, nlons, nazmth, iheatcal, icutoff
  INTEGER nmessg, nsmax, ierror
  REAL kstar(nlons), slope, f1, f2, f3, f5, f6, alt_cutoff, smco
  REAL k_alpha(nlons, nazmth), coslat(nlons)
  REAL ksmin, ksmax

  ! Internal variables.

  INTEGER i, n
  ! -----------------------------------------------------------------------

  ! Specify constants.

  naz = 8
  slope = 1.
  f1 = 1.5
  f2 = 0.3
  f3 = 1.0
  f5 = 3.0
  f6 = 1.0
  ksmin = 1.E-5
  ksmax = 1.E-4
  DO i = 1, nlons
    kstar(i) = ksmin/(coslat(i)+(ksmin/ksmax))
  END DO
  icutoff = 1
  alt_cutoff = 105.E3
  smco = 2.0
  ! SMCO       = 1.0
  nsmax = 5
  ! NSMAX      = 2
  iheatcal = 0

  ! Print information to output file.

  ! WRITE (NMESSG,6000)
  ! 6000 FORMAT (/' Subroutine HINES_SETUP:')
  ! WRITE (NMESSG,*)  '  SLOPE = ', SLOPE
  ! WRITE (NMESSG,*)  '  NAZ = ', NAZ
  ! WRITE (NMESSG,*)  '  F1,F2,F3  = ', F1, F2, F3
  ! WRITE (NMESSG,*)  '  F5,F6     = ', F5, F6
  ! WRITE (NMESSG,*)  '  KSTAR     = ', KSTAR
  ! >           ,'  COSLAT     = ', COSLAT
  ! IF (ICUTOFF .EQ. 1)  THEN
  ! WRITE (NMESSG,*) '  Drag exponentially damped above ',
  ! &                       ALT_CUTOFF/1.E3
  ! END IF
  ! IF (NSMAX.LT.1 )  THEN
  ! WRITE (NMESSG,*) '  No smoothing of cutoff wavenumbers, etc'
  ! ELSE
  ! WRITE (NMESSG,*) '  Cutoff wavenumbers and sig_t smoothed:'
  ! WRITE (NMESSG,*) '    SMCO  =', SMCO
  ! WRITE (NMESSG,*) '    NSMAX =', NSMAX
  ! END IF

  ! Check that things are setup correctly and log error if not

  ierror = 0
  IF (naz>nazmth) ierror = 10
  IF (naz/=4 .AND. naz/=8) ierror = 20
  IF (slope/=1. .AND. slope/=1.5 .AND. slope/=2.) ierror = 30
  IF (smco<1.) ierror = 40

  ! Use single value for azimuthal-dependent horizontal wavenumber.

  DO n = 1, naz
    DO i = 1, nlons
      k_alpha(i, n) = kstar(i)
    END DO
  END DO

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE hines_setup

SUBROUTINE hines_print(flux_u, flux_v, drag_u, drag_v, alt, sigma_t, &
    sigma_alpha, v_alpha, m_alpha, iu_print, iv_print, nmessg, ilprt1, &
    ilprt2, levprt1, levprt2, naz, nlons, nlevs, nazmth)
  IMPLICIT NONE
  ! Print out altitude profiles of various quantities from
  ! Hines' Doppler spread gravity wave drag parameterization scheme.
  ! (NOTE: only for NAZ = 4 or 8).

  ! Aug. 8/95 - C. McLandress

  ! Input arguements:

  ! * IU_PRINT = 1 to print out values in east-west direction.
  ! * IV_PRINT = 1 to print out values in north-south direction.
  ! * NMESSG   = unit number for printed output.
  ! * ILPRT1   = first longitudinal index to print.
  ! * ILPRT2   = last longitudinal index to print.
  ! * LEVPRT1  = first altitude level to print.
  ! * LEVPRT2  = last altitude level to print.

  INTEGER naz, ilprt1, ilprt2, levprt1, levprt2
  INTEGER nlons, nlevs, nazmth
  INTEGER iu_print, iv_print, nmessg
  REAL flux_u(nlons, nlevs), flux_v(nlons, nlevs)
  REAL drag_u(nlons, nlevs), drag_v(nlons, nlevs)
  REAL alt(nlons, nlevs), sigma_t(nlons, nlevs)
  REAL sigma_alpha(nlons, nlevs, nazmth)
  REAL v_alpha(nlons, nlevs, nazmth), m_alpha(nlons, nlevs, nazmth)

  ! Internal variables.

  INTEGER n_east, n_west, n_north, n_south
  INTEGER i, l
  ! -----------------------------------------------------------------------

  ! Azimuthal indices of cardinal directions.

  n_east = 1
  IF (naz==4) THEN
    n_west = 3
    n_north = 2
    n_south = 4
  ELSE IF (naz==8) THEN
    n_west = 5
    n_north = 3
    n_south = 7
  END IF

  ! Print out values for range of longitudes.

  DO i = ilprt1, ilprt2

    ! Print east-west wind, sigmas, cutoff wavenumbers, flux and drag.

    IF (iu_print==1) THEN
      WRITE (nmessg, *)
      WRITE (nmessg, 6001) i
      WRITE (nmessg, 6005)
6001  FORMAT ('Hines GW (east-west) at longitude I =', I3)
6005  FORMAT (15X, ' U ', 2X, 'sig_E', 2X, 'sig_T', 3X, 'm_E', 4X, 'm_W', 4X, &
        'fluxU', 5X, 'gwdU')
      DO l = levprt1, levprt2
        WRITE (nmessg, 6701) alt(i, l)/1.E3, v_alpha(i, l, n_east), &
          sigma_alpha(i, l, n_east), sigma_t(i, l), &
          m_alpha(i, l, n_east)*1.E3, m_alpha(i, l, n_west)*1.E3, &
          flux_u(i, l)*1.E5, drag_u(i, l)*24.*3600.
      END DO
6701  FORMAT (' z=', F7.2, 1X, 3F7.1, 2F7.3, F9.4, F9.3)
    END IF

    ! Print north-south winds, sigmas, cutoff wavenumbers, flux and drag.

    IF (iv_print==1) THEN
      WRITE (nmessg, *)
      WRITE (nmessg, 6002) i
6002  FORMAT ('Hines GW (north-south) at longitude I =', I3)
      WRITE (nmessg, 6006)
6006  FORMAT (15X, ' V ', 2X, 'sig_N', 2X, 'sig_T', 3X, 'm_N', 4X, 'm_S', 4X, &
        'fluxV', 5X, 'gwdV')
      DO l = levprt1, levprt2
        WRITE (nmessg, 6701) alt(i, l)/1.E3, v_alpha(i, l, n_north), &
          sigma_alpha(i, l, n_north), sigma_t(i, l), &
          m_alpha(i, l, n_north)*1.E3, m_alpha(i, l, n_south)*1.E3, &
          flux_v(i, l)*1.E5, drag_v(i, l)*24.*3600.
      END DO
    END IF

  END DO

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE hines_print

SUBROUTINE hines_exp(data, data_zmax, alt, alt_exp, iorder, il1, il2, lev1, &
    lev2, nlons, nlevs)
  IMPLICIT NONE
  ! This routine exponentially damps a longitude by altitude array
  ! of data above a specified altitude.

  ! Aug. 13/95 - C. McLandress

  ! Output arguements:

  ! * DATA = modified data array.

  ! Input arguements:

  ! * DATA    = original data array.
  ! * ALT     = altitudes.
  ! * ALT_EXP = altitude above which exponential decay applied.
  ! * IORDER	= 1 means vertical levels are indexed from top down
  ! *           (i.e., highest level indexed 1 and lowest level NLEVS);
  ! *           .NE. 1 highest level is index NLEVS.
  ! * IL1     = first longitudinal index to use (IL1 >= 1).
  ! * IL2     = last longitudinal index to use (IL1 <= IL2 <= NLONS).
  ! * LEV1    = first altitude level to use (LEV1 >=1).
  ! * LEV2    = last altitude level to use (LEV1 < LEV2 <= NLEVS).
  ! * NLONS   = number of longitudes.
  ! * NLEVS   = number of vertical

  ! Input work arrays:

  ! * DATA_ZMAX = data values just above altitude ALT_EXP.

  INTEGER iorder, il1, il2, lev1, lev2, nlons, nlevs
  REAL alt_exp
  REAL data(nlons, nlevs), data_zmax(nlons), alt(nlons, nlevs)

  ! Internal variables.

  INTEGER levbot, levtop, lincr, i, l
  REAL hscale
  DATA hscale/5.E3/
  ! -----------------------------------------------------------------------

  ! Index of lowest altitude level (bottom of drag calculation).

  levbot = lev2
  levtop = lev1
  lincr = 1
  IF (iorder/=1) THEN
    levbot = lev1
    levtop = lev2
    lincr = -1
  END IF

  ! Data values at first level above ALT_EXP.

  DO i = il1, il2
    DO l = levtop, levbot, lincr
      IF (alt(i,l)>=alt_exp) THEN
        data_zmax(i) = data(i, l)
      END IF
    END DO
  END DO

  ! Exponentially damp field above ALT_EXP to model top at L=1.

  DO l = 1, lev2
    DO i = il1, il2
      IF (alt(i,l)>=alt_exp) THEN
        data(i, l) = data_zmax(i)*exp((alt_exp-alt(i,l))/hscale)
      END IF
    END DO
  END DO

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE hines_exp

SUBROUTINE vert_smooth(data, work, coeff, nsmooth, il1, il2, lev1, lev2, &
    nlons, nlevs)
  IMPLICIT NONE
  ! Smooth a longitude by altitude array in the vertical over a
  ! specified number of levels using a three point smoother.

  ! NOTE: input array DATA is modified on output!

  ! Aug. 3/95 - C. McLandress

  ! Output arguement:

  ! * DATA    = smoothed array (on output).

  ! Input arguements:

  ! * DATA    = unsmoothed array of data (on input).
  ! * WORK    = work array of same dimension as DATA.
  ! * COEFF   = smoothing coefficient for a 1:COEFF:1 stencil.
  ! *           (e.g., COEFF = 2 will result in a smoother which
  ! *           weights the level L gridpoint by two and the two
  ! *           adjecent levels (L+1 and L-1) by one).
  ! * NSMOOTH = number of times to smooth in vertical.
  ! *           (e.g., NSMOOTH=1 means smoothed only once,
  ! *           NSMOOTH=2 means smoothing repeated twice, etc.)
  ! * IL1     = first longitudinal index to use (IL1 >= 1).
  ! * IL2     = last longitudinal index to use (IL1 <= IL2 <= NLONS).
  ! * LEV1    = first altitude level to use (LEV1 >=1).
  ! * LEV2    = last altitude level to use (LEV1 < LEV2 <= NLEVS).
  ! * NLONS   = number of longitudes.
  ! * NLEVS   = number of vertical levels.

  ! Subroutine arguements.

  INTEGER nsmooth, il1, il2, lev1, lev2, nlons, nlevs
  REAL coeff
  REAL data(nlons, nlevs), work(nlons, nlevs)

  ! Internal variables.

  INTEGER i, l, ns, lev1p, lev2m
  REAL sum_wts
  ! -----------------------------------------------------------------------

  ! Calculate sum of weights.

  sum_wts = coeff + 2.

  lev1p = lev1 + 1
  lev2m = lev2 - 1

  ! Smooth NSMOOTH times

  DO ns = 1, nsmooth

    ! Copy data into work array.

    DO l = lev1, lev2
      DO i = il1, il2
        work(i, l) = data(i, l)
      END DO
    END DO

    ! Smooth array WORK in vertical direction and put into DATA.

    DO l = lev1p, lev2m
      DO i = il1, il2
        data(i, l) = (work(i,l+1)+coeff*work(i,l)+work(i,l-1))/sum_wts
      END DO
    END DO

  END DO

  RETURN
  ! -----------------------------------------------------------------------
END SUBROUTINE vert_smooth






