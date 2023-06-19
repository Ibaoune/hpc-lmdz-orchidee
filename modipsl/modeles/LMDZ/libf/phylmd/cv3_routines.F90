
! $Id: cv3_routines.F90 3670 2020-04-27 08:49:09Z jyg $




SUBROUTINE cv3_param(nd, k_upper, delt)

  USE ioipsl_getin_p_mod, ONLY : getin_p
  use mod_phys_lmdz_para
  IMPLICIT NONE

!------------------------------------------------------------
!Set parameters for convectL for iflag_con = 3
!------------------------------------------------------------


!***  PBCRIT IS THE CRITICAL CLOUD DEPTH (MB) BENEATH WHICH THE ***
!***      PRECIPITATION EFFICIENCY IS ASSUMED TO BE ZERO     ***
!***  PTCRIT IS THE CLOUD DEPTH (MB) ABOVE WHICH THE PRECIP. ***
!***            EFFICIENCY IS ASSUMED TO BE UNITY            ***
!***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
!***  SPFAC IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE ***
!***                        OF CLOUD                         ***

![TAU: CHARACTERISTIC TIMESCALE USED TO COMPUTE ALPHA & BETA]
!***    ALPHA AND BETA ARE PARAMETERS THAT CONTROL THE RATE OF ***
!***                 APPROACH TO QUASI-EQUILIBRIUM           ***
!***    (THEIR STANDARD VALUES ARE 1.0 AND 0.96, RESPECTIVELY) ***
!***           (BETA MUST BE LESS THAN OR EQUAL TO 1)        ***

!***    DTCRIT IS THE CRITICAL BUOYANCY (K) USED TO ADJUST THE ***
!***                 APPROACH TO QUASI-EQUILIBRIUM           ***
!***                     IT MUST BE LESS THAN 0              ***

  include "cv3param.h"
  include "cvflag.h"
  include "conema3.h"

  INTEGER, INTENT(IN)              :: nd
  INTEGER, INTENT(IN)              :: k_upper
  REAL, INTENT(IN)                 :: delt ! timestep (seconds)

! Local variables
  CHARACTER (LEN=20) :: modname = 'cv3_param'
  CHARACTER (LEN=80) :: abort_message

  LOGICAL, SAVE :: first = .TRUE.
!$OMP THREADPRIVATE(first)

!glb  noff: integer limit for convection (nd-noff)
! minorig: First level of convection

! -- limit levels for convection:

!jyg<
!  noff is chosen such that nl = k_upper so that upmost loops end at about 22 km
!
  noff = min(max(nd-k_upper, 1), (nd+1)/2)
!!  noff = 1
!>jyg
  minorig = 1
  nl = nd - noff
  nlp = nl + 1
  nlm = nl - 1

  IF (first) THEN
! -- "microphysical" parameters:
! IM beg: ajout fis. reglage ep
! CR+JYG: shedding coefficient (used when iflag_mix_adiab=1)
! IM lu dans physiq.def via conf_phys.F90     epmax  = 0.993

    omtrain = 45.0 ! used also for snow (no disctinction rain/snow)
! -- misc:
    dtovsh = -0.2 ! dT for overshoot
! cc      dttrig = 5.   ! (loose) condition for triggering
    dttrig = 10. ! (loose) condition for triggering
    dtcrit = -2.0
! -- end of convection
! -- interface cloud parameterization:
    delta = 0.01 ! cld
! -- interface with boundary-layer (gust factor): (sb)
    betad = 10.0 ! original value (from convect 4.3)

! Var interm pour le getin
     cv_flag_feed=1
     CALL getin_p('cv_flag_feed',cv_flag_feed)
     T_top_max = 1000.
     CALL getin_p('t_top_max',T_top_max)
     dpbase=-40.
     CALL getin_p('dpbase',dpbase)
     pbcrit=150.0
     CALL getin_p('pbcrit',pbcrit)
     ptcrit=500.0
     CALL getin_p('ptcrit',ptcrit)
     sigdz=0.01
     CALL getin_p('sigdz',sigdz)
     spfac=0.15
     CALL getin_p('spfac',spfac)
     tau=8000.
     CALL getin_p('tau',tau)
     flag_wb=1
     CALL getin_p('flag_wb',flag_wb)
     wbmax=6.
     CALL getin_p('wbmax',wbmax)
     ok_convstop=.False.
     CALL getin_p('ok_convstop',ok_convstop)
     tau_stop=15000.
     CALL getin_p('tau_stop',tau_stop)
     ok_intermittent=.False.
     CALL getin_p('ok_intermittent',ok_intermittent)
     ok_optim_yield=.False.
     CALL getin_p('ok_optim_yield',ok_optim_yield)
     ok_homo_tend=.TRUE.
     CALL getin_p('ok_homo_tend',ok_homo_tend)
     ok_entrain=.TRUE.
     CALL getin_p('ok_entrain',ok_entrain)

     coef_peel=0.25
     CALL getin_p('coef_peel',coef_peel)

     flag_epKEorig=1
     CALL getin_p('flag_epKEorig',flag_epKEorig)
     elcrit=0.0003
     CALL getin_p('elcrit',elcrit)
     tlcrit=-55.0
     CALL getin_p('tlcrit',tlcrit)
     ejectliq=0.
     CALL getin_p('ejectliq',ejectliq)
     ejectice=0.
     CALL getin_p('ejectice',ejectice)
     cvflag_prec_eject = .FALSE.
     CALL getin_p('cvflag_prec_eject',cvflag_prec_eject)
     qsat_depends_on_qt = .FALSE.
     CALL getin_p('qsat_depends_on_qt',qsat_depends_on_qt)
     adiab_ascent_mass_flux_depends_on_ejectliq = .FALSE.
     CALL getin_p('adiab_ascent_mass_flux_depends_on_ejectliq',adiab_ascent_mass_flux_depends_on_ejectliq)
     keepbug_ice_frac = .TRUE.
     CALL getin_p('keepbug_ice_frac', keepbug_ice_frac)

    WRITE (*, *) 't_top_max=', t_top_max
    WRITE (*, *) 'dpbase=', dpbase
    WRITE (*, *) 'pbcrit=', pbcrit
    WRITE (*, *) 'ptcrit=', ptcrit
    WRITE (*, *) 'sigdz=', sigdz
    WRITE (*, *) 'spfac=', spfac
    WRITE (*, *) 'tau=', tau
    WRITE (*, *) 'flag_wb=', flag_wb
    WRITE (*, *) 'wbmax=', wbmax
    WRITE (*, *) 'ok_convstop=', ok_convstop
    WRITE (*, *) 'tau_stop=', tau_stop
    WRITE (*, *) 'ok_intermittent=', ok_intermittent
    WRITE (*, *) 'ok_optim_yield =', ok_optim_yield
    WRITE (*, *) 'coef_peel=', coef_peel

    WRITE (*, *) 'flag_epKEorig=', flag_epKEorig
    WRITE (*, *) 'elcrit=', elcrit
    WRITE (*, *) 'tlcrit=', tlcrit
    WRITE (*, *) 'ejectliq=', ejectliq
    WRITE (*, *) 'ejectice=', ejectice
    WRITE (*, *) 'cvflag_prec_eject =', cvflag_prec_eject 
    WRITE (*, *) 'qsat_depends_on_qt =', qsat_depends_on_qt 
    WRITE (*, *) 'adiab_ascent_mass_flux_depends_on_ejectliq =', adiab_ascent_mass_flux_depends_on_ejectliq
    WRITE (*, *) 'keepbug_ice_frac =', keepbug_ice_frac 

    first = .FALSE.
  END IF ! (first)

  beta = 1.0 - delt/tau
  alpha1 = 1.5E-3
!JYG    Correction bug alpha
  alpha1 = alpha1*1.5
  alpha = alpha1*delt/tau
!JYG    Bug
! cc increase alpha to compensate W decrease:
! c      alpha  = alpha*1.5

  noconv_stop = max(2.,tau_stop/delt)

  RETURN
END SUBROUTINE cv3_param

SUBROUTINE cv3_incrcount(len, nd, delt, sig)

IMPLICIT NONE

! =====================================================================
!  Increment the counter sig(nd)
! =====================================================================

  include "cv3param.h"
  include "cvflag.h"

!inputs:
  INTEGER, INTENT(IN)                     :: len
  INTEGER, INTENT(IN)                     :: nd
  REAL, INTENT(IN)                        :: delt ! timestep (seconds)

!input/output
  REAL, DIMENSION(len,nd), INTENT(INOUT)  :: sig

!local variables
  INTEGER il

!    print *,'cv3_incrcount : noconv_stop ',noconv_stop
!    print *,'cv3_incrcount in, sig(1,nd) ',sig(1,nd)
    IF(ok_convstop) THEN
      DO il = 1, len
        sig(il, nd) = sig(il, nd) + 1.
        sig(il, nd) = min(sig(il,nd), noconv_stop+0.1)
      END DO
    ELSE
      DO il = 1, len
        sig(il, nd) = sig(il, nd) + 1.
        sig(il, nd) = min(sig(il,nd), 12.1)
      END DO
    ENDIF  ! (ok_convstop)
!    print *,'cv3_incrcount out, sig(1,nd) ',sig(1,nd)

  RETURN
END SUBROUTINE cv3_incrcount

SUBROUTINE cv3_prelim(len, nd, ndp1, t, q, p, ph, &
                      lv, lf, cpn, tv, gz, h, hm, th)
  IMPLICIT NONE

! =====================================================================
! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
! "ori": from convect4.3 (vectorized)
! "convect3": to be exactly consistent with convect3
! =====================================================================

! inputs:
  INTEGER len, nd, ndp1
  REAL t(len, nd), q(len, nd), p(len, nd), ph(len, ndp1)

! outputs:
  REAL lv(len, nd), lf(len, nd), cpn(len, nd), tv(len, nd)
  REAL gz(len, nd), h(len, nd), hm(len, nd)
  REAL th(len, nd)

! local variables:
  INTEGER k, i
  REAL rdcp
  REAL tvx, tvy ! convect3
  REAL cpx(len, nd)

  include "cvthermo.h"
  include "cv3param.h"


! ori      do 110 k=1,nlp
! abderr     do 110 k=1,nl ! convect3
  DO k = 1, nlp

    DO i = 1, len
! debug          lv(i,k)= lv0-clmcpv*(t(i,k)-t0)
      lv(i, k) = lv0 - clmcpv*(t(i,k)-273.15)
!!      lf(i, k) = lf0 - clmci*(t(i,k)-273.15)   ! erreur de signe !!
      lf(i, k) = lf0 + clmci*(t(i,k)-273.15)
      cpn(i, k) = cpd*(1.0-q(i,k)) + cpv*q(i, k)
      cpx(i, k) = cpd*(1.0-q(i,k)) + cl*q(i, k)
! ori          tv(i,k)=t(i,k)*(1.0+q(i,k)*epsim1)
      tv(i, k) = t(i, k)*(1.0+q(i,k)/eps-q(i,k))
      rdcp = (rrd*(1.-q(i,k))+q(i,k)*rrv)/cpn(i, k)
      th(i, k) = t(i, k)*(1000.0/p(i,k))**rdcp
    END DO
  END DO

! gz = phi at the full levels (same as p).

!!  DO i = 1, len                    !jyg
!!    gz(i, 1) = 0.0                 !jyg
!!  END DO                           !jyg
    gz(:,:) = 0.                     !jyg: initialization of the whole array
! ori      do 140 k=2,nlp
  DO k = 2, nl ! convect3
    DO i = 1, len
      tvx = t(i, k)*(1.+q(i,k)/eps-q(i,k))         !convect3
      tvy = t(i, k-1)*(1.+q(i,k-1)/eps-q(i,k-1))   !convect3
      gz(i, k) = gz(i, k-1) + 0.5*rrd*(tvx+tvy)* & !convect3
                 (p(i,k-1)-p(i,k))/ph(i, k)        !convect3

! c        print *,' gz(',k,')',gz(i,k),' tvx',tvx,' tvy ',tvy

! ori         gz(i,k)=gz(i,k-1)+hrd*(tv(i,k-1)+tv(i,k))
! ori    &         *(p(i,k-1)-p(i,k))/ph(i,k)
    END DO
  END DO

! h  = phi + cpT (dry static energy).
! hm = phi + cp(T-Tbase)+Lq

! ori      do 170 k=1,nlp
  DO k = 1, nl ! convect3
    DO i = 1, len
      h(i, k) = gz(i, k) + cpn(i, k)*t(i, k)
      hm(i, k) = gz(i, k) + cpx(i, k)*(t(i,k)-t(i,1)) + lv(i, k)*q(i, k)
    END DO
  END DO

  RETURN
END SUBROUTINE cv3_prelim

SUBROUTINE cv3_feed(len, nd, ok_conserv_q, &
                    t, q, u, v, p, ph, h, gz, &
                    p1feed, p2feed, wght, &
                    wghti, tnk, thnk, qnk, qsnk, unk, vnk, &
                    cpnk, hnk, nk, icb, icbmax, iflag, gznk, plcl)

  USE mod_phys_lmdz_transfert_para, ONLY : bcast
  USE add_phys_tend_mod, ONLY: fl_cor_ebil
  USE print_control_mod, ONLY: prt_level
  IMPLICIT NONE

! ================================================================
! Purpose: CONVECTIVE FEED

! Main differences with cv_feed:
! - ph added in input
! - here, nk(i)=minorig
! - icb defined differently (plcl compared with ph instead of p)
! - dry static energy as argument instead of moist static energy

! Main differences with convect3:
! - we do not compute dplcldt and dplcldr of CLIFT anymore
! - values iflag different (but tests identical)
! - A,B explicitely defined (!...)
! ================================================================

  include "cv3param.h"
  include "cvthermo.h"

!inputs:
  INTEGER, INTENT (IN)                               :: len, nd
  LOGICAL, INTENT (IN)                               :: ok_conserv_q
  REAL, DIMENSION (len, nd), INTENT (IN)             :: t, q, p
  REAL, DIMENSION (len, nd), INTENT (IN)             :: u, v
  REAL, DIMENSION (len, nd), INTENT (IN)             :: h, gz
  REAL, DIMENSION (len, nd+1), INTENT (IN)           :: ph
  REAL, DIMENSION (len), INTENT (IN)                 :: p1feed
  REAL, DIMENSION (nd), INTENT (IN)                  :: wght
!input-output
  REAL, DIMENSION (len), INTENT (INOUT)              :: p2feed
!outputs:
  INTEGER, INTENT (OUT)                              :: icbmax
  INTEGER, DIMENSION (len), INTENT (OUT)             :: iflag, nk, icb
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: wghti
  REAL, DIMENSION (len), INTENT (OUT)                :: tnk, thnk, qnk, qsnk
  REAL, DIMENSION (len), INTENT (OUT)                :: unk, vnk
  REAL, DIMENSION (len), INTENT (OUT)                :: cpnk, hnk, gznk
  REAL, DIMENSION (len), INTENT (OUT)                :: plcl

!local variables:
  INTEGER i, k, iter, niter
  INTEGER ihmin(len)
  REAL work(len)
  REAL pup(len), plo(len), pfeed(len)
  REAL plclup(len), plcllo(len), plclfeed(len)
  REAL pfeedmin(len)
  REAL posit(len)
  LOGICAL nocond(len)

!jyg20140217<
  INTEGER iostat
  LOGICAL, SAVE :: first
  LOGICAL, SAVE :: ok_new_feed
  REAL, SAVE :: dp_lcl_feed
!$OMP THREADPRIVATE (first,ok_new_feed,dp_lcl_feed)
  DATA first/.TRUE./
  DATA dp_lcl_feed/2./

  IF (first) THEN
!$OMP MASTER
    ok_new_feed = ok_conserv_q
    OPEN (98, FILE='cv3feed_param.data', STATUS='old', FORM='formatted', IOSTAT=iostat)
    IF (iostat==0) THEN
      READ (98, *, END=998) ok_new_feed
998   CONTINUE
      CLOSE (98)
    END IF
    PRINT *, ' ok_new_feed: ', ok_new_feed
!$OMP END MASTER
    call bcast(ok_new_feed)
    first = .FALSE.   
  END IF
!jyg>
! -------------------------------------------------------------------
! --- Origin level of ascending parcels for convect3:
! -------------------------------------------------------------------

  DO i = 1, len
    nk(i) = minorig
    gznk(i) = gz(i, nk(i))
  END DO

! -------------------------------------------------------------------
! --- Adjust feeding layer thickness so that lifting up to the top of
! --- the feeding layer does not induce condensation (i.e. so that
! --- plcl < p2feed).
! --- Method : iterative secant method.
! -------------------------------------------------------------------

! 1- First bracketing of the solution : ph(nk+1), p2feed

! 1.a- LCL associated with p2feed
  DO i = 1, len
    pup(i) = p2feed(i)
  END DO
  IF (fl_cor_ebil >=2 ) THEN
    CALL cv3_estatmix(len, nd, iflag, p1feed, pup, p, ph, &
                     t, q, u, v, h, gz, wght, &
                     wghti, nk, tnk, thnk, qnk, qsnk, unk, vnk, plclup)
  ELSE
    CALL cv3_enthalpmix(len, nd, iflag, p1feed, pup, p, ph, &
                       t, q, u, v, wght, &
                       wghti, nk, tnk, thnk, qnk, qsnk, unk, vnk, plclup)
  ENDIF  ! (fl_cor_ebil >=2 ) 
! 1.b- LCL associated with ph(nk+1)
  DO i = 1, len
    plo(i) = ph(i, nk(i)+1)
  END DO
  IF (fl_cor_ebil >=2 ) THEN
    CALL cv3_estatmix(len, nd, iflag, p1feed, plo, p, ph, &
                     t, q, u, v, h, gz, wght, &
                     wghti, nk, tnk, thnk, qnk, qsnk, unk, vnk, plcllo)
  ELSE
    CALL cv3_enthalpmix(len, nd, iflag, p1feed, plo, p, ph, &
                       t, q, u, v, wght, &
                       wghti, nk, tnk, thnk, qnk, qsnk, unk, vnk, plcllo)
  ENDIF  ! (fl_cor_ebil >=2 ) 
! 2- Iterations
  niter = 5
  DO iter = 1, niter
    DO i = 1, len
      plcllo(i) = min(plo(i), plcllo(i))
      plclup(i) = max(pup(i), plclup(i))
      nocond(i) = plclup(i) <= pup(i)
    END DO
    DO i = 1, len
      IF (nocond(i)) THEN
        pfeed(i) = pup(i)
      ELSE
!JYG20140217<
        IF (ok_new_feed) THEN
          pfeed(i) = (pup(i)*(plo(i)-plcllo(i)-dp_lcl_feed)+  &
                      plo(i)*(plclup(i)-pup(i)+dp_lcl_feed))/ &
                     (plo(i)-plcllo(i)+plclup(i)-pup(i))
        ELSE
          pfeed(i) = (pup(i)*(plo(i)-plcllo(i))+  &
                      plo(i)*(plclup(i)-pup(i)))/ &
                     (plo(i)-plcllo(i)+plclup(i)-pup(i))
        END IF
!JYG>
      END IF
    END DO
!jyg20140217<
! For the last iteration, make sure that the top of the feeding layer
! and LCL are not in the same layer:
    IF (ok_new_feed) THEN
      IF (iter==niter) THEN
        DO i = 1,len                         !jyg
          pfeedmin(i) = ph(i,minorig+1)      !jyg
        ENDDO                                !jyg
        DO k = minorig+1, nl                 !jyg
!!        DO k = minorig, nl                 !jyg
          DO i = 1, len
            IF (ph(i,k)>=plclfeed(i)) pfeedmin(i) = ph(i, k)
          END DO
        END DO
        DO i = 1, len
          pfeed(i) = max(pfeedmin(i), pfeed(i))
        END DO
      END IF
    END IF
!jyg>

    IF (fl_cor_ebil >=2 ) THEN
      CALL cv3_estatmix(len, nd, iflag, p1feed, pfeed, p, ph, &
                       t, q, u, v, h, gz, wght, &
                       wghti, nk, tnk, thnk, qnk, qsnk, unk, vnk, plclfeed)
    ELSE
      CALL cv3_enthalpmix(len, nd, iflag, p1feed, pfeed, p, ph, &
                         t, q, u, v, wght, &
                         wghti, nk, tnk, thnk, qnk, qsnk, unk, vnk, plclfeed)
    ENDIF  ! (fl_cor_ebil >=2 ) 
!jyg20140217<
    IF (ok_new_feed) THEN
      DO i = 1, len
        posit(i) = (sign(1.,plclfeed(i)-pfeed(i)+dp_lcl_feed)+1.)*0.5
        IF (plclfeed(i)-pfeed(i)+dp_lcl_feed==0.) posit(i) = 1.
      END DO
    ELSE
      DO i = 1, len
        posit(i) = (sign(1.,plclfeed(i)-pfeed(i))+1.)*0.5
        IF (plclfeed(i)==pfeed(i)) posit(i) = 1.
      END DO
    END IF
!jyg>
    DO i = 1, len
! - posit = 1 when lcl is below top of feeding layer (plclfeed>pfeed)
! -               => pup=pfeed
! - posit = 0 when lcl is above top of feeding layer (plclfeed<pfeed)
! -               => plo=pfeed
      pup(i) = posit(i)*pfeed(i) + (1.-posit(i))*pup(i)
      plo(i) = (1.-posit(i))*pfeed(i) + posit(i)*plo(i)
      plclup(i) = posit(i)*plclfeed(i) + (1.-posit(i))*plclup(i)
      plcllo(i) = (1.-posit(i))*plclfeed(i) + posit(i)*plcllo(i)
    END DO
  END DO !  iter

  DO i = 1, len
    p2feed(i) = pfeed(i)
    plcl(i) = plclfeed(i)
  END DO

  DO i = 1, len
    cpnk(i) = cpd*(1.0-qnk(i)) + cpv*qnk(i)
    hnk(i) = gz(i, 1) + cpnk(i)*tnk(i)
  END DO

! -------------------------------------------------------------------
! --- Check whether parcel level temperature and specific humidity
! --- are reasonable
! -------------------------------------------------------------------
  IF (cv_flag_feed == 1) THEN
    DO i = 1, len
      IF (((tnk(i)<250.0)                       .OR.  &
           (qnk(i)<=0.0))                       .AND. &
          (iflag(i)==0)) iflag(i) = 7
    END DO
  ELSEIF (cv_flag_feed >= 2) THEN
! --- and demand that LCL be high enough
    DO i = 1, len
      IF (((tnk(i)<250.0)                       .OR.  &
           (qnk(i)<=0.0)                        .OR.  &
           (plcl(i)>min(0.99*ph(i,1),ph(i,3)))) .AND. &
          (iflag(i)==0)) iflag(i) = 7
    END DO
  ENDIF
  IF (prt_level .GE. 10) THEN
    print *,'cv3_feed : iflag(1), pfeed(1), plcl(1), wghti(1,k) ', &
                        iflag(1), pfeed(1), plcl(1), (wghti(1,k),k=1,10)
  ENDIF

! -------------------------------------------------------------------
! --- Calculate first level above lcl (=icb)
! -------------------------------------------------------------------

!@      do 270 i=1,len
!@       icb(i)=nlm
!@ 270  continue
!@c
!@      do 290 k=minorig,nl
!@        do 280 i=1,len
!@          if((k.ge.(nk(i)+1)).and.(p(i,k).lt.plcl(i)))
!@     &    icb(i)=min(icb(i),k)
!@ 280    continue
!@ 290  continue
!@c
!@      do 300 i=1,len
!@        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
!@ 300  continue

  DO i = 1, len
    icb(i) = nlm
  END DO

! la modification consiste a comparer plcl a ph et non a p:
! icb est defini par :  ph(icb)<plcl<ph(icb-1)
!@      do 290 k=minorig,nl
  DO k = 3, nl - 1 ! modif pour que icb soit sup/egal a 2
    DO i = 1, len
      IF (ph(i,k)<plcl(i)) icb(i) = min(icb(i), k)
    END DO
  END DO


! print*,'icb dans cv3_feed '
! write(*,'(64i2)') icb(2:len-1)
! call dump2d(64,43,'plcl dans cv3_feed ',plcl(2:len-1))

  DO i = 1, len
!@        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
    IF ((icb(i)==nlm) .AND. (iflag(i)==0)) iflag(i) = 9
  END DO

  DO i = 1, len
    icb(i) = icb(i) - 1 ! icb sup ou egal a 2
  END DO

! Compute icbmax.

  icbmax = 2
  DO i = 1, len
!!        icbmax=max(icbmax,icb(i))
    IF (iflag(i)<7) icbmax = max(icbmax, icb(i))     ! sb Jun7th02
  END DO

  RETURN
END SUBROUTINE cv3_feed

SUBROUTINE cv3_undilute1(len, nd, t, qs, gz, plcl, p, icb, tnk, qnk, gznk, &
                         tp, tvp, clw, icbs)
  IMPLICIT NONE

! ----------------------------------------------------------------
! Equivalent de TLIFT entre NK et ICB+1 inclus

! Differences with convect4:
!    - specify plcl in input
!    - icbs is the first level above LCL (may differ from icb)
!    - in the iterations, used x(icbs) instead x(icb)
!    - many minor differences in the iterations
!    - tvp is computed in only one time
!    - icbs: first level above Plcl (IMIN de TLIFT) in output
!    - if icbs=icb, compute also tp(icb+1),tvp(icb+1) & clw(icb+1)
! ----------------------------------------------------------------

  include "cvthermo.h"
  include "cv3param.h"

! inputs:
  INTEGER, INTENT (IN)                              :: len, nd
  INTEGER, DIMENSION (len), INTENT (IN)             :: icb
  REAL, DIMENSION (len, nd), INTENT (IN)            :: t, qs, gz
  REAL, DIMENSION (len), INTENT (IN)                :: tnk, qnk, gznk
  REAL, DIMENSION (len, nd), INTENT (IN)            :: p
  REAL, DIMENSION (len), INTENT (IN)                :: plcl              ! convect3

! outputs:
  INTEGER, DIMENSION (len), INTENT (OUT)            :: icbs
  REAL, DIMENSION (len, nd), INTENT (OUT)           :: tp, tvp, clw

! local variables:
  INTEGER i, k
  INTEGER icb1(len), icbsmax2                                            ! convect3
  REAL tg, qg, alv, s, ahg, tc, denom, es, rg
  REAL ah0(len), cpp(len)
  REAL ticb(len), gzicb(len)
  REAL qsicb(len)                                                        ! convect3
  REAL cpinv(len)                                                        ! convect3

! -------------------------------------------------------------------
! --- Calculates the lifted parcel virtual temperature at nk,
! --- the actual temperature, and the adiabatic
! --- liquid water content. The procedure is to solve the equation.
!     cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.
! -------------------------------------------------------------------


! ***  Calculate certain parcel quantities, including static energy   ***

  DO i = 1, len
    ah0(i) = (cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i) + qnk(i)*(lv0-clmcpv*(tnk(i)-273.15)) + gznk(i)
    cpp(i) = cpd*(1.-qnk(i)) + qnk(i)*cpv
    cpinv(i) = 1./cpp(i)
  END DO

! ***   Calculate lifted parcel quantities below cloud base   ***

  DO i = 1, len                                           !convect3
    icb1(i) = min(max(icb(i), 2), nl)
! if icb is below LCL, start loop at ICB+1:
! (icbs est le premier niveau au-dessus du LCL)
    icbs(i) = icb1(i)                                     !convect3
    IF (plcl(i)<p(i,icb1(i))) THEN
      icbs(i) = min(icbs(i)+1, nl)                        !convect3
    END IF
  END DO                                                  !convect3

  DO i = 1, len !convect3
    ticb(i) = t(i, icbs(i))                               !convect3
    gzicb(i) = gz(i, icbs(i))                             !convect3
    qsicb(i) = qs(i, icbs(i))                             !convect3
  END DO !convect3


! Re-compute icbsmax (icbsmax2):                          !convect3
!                                                         !convect3
  icbsmax2 = 2                                            !convect3
  DO i = 1, len                                           !convect3
    icbsmax2 = max(icbsmax2, icbs(i))                     !convect3
  END DO                                                  !convect3

! initialization outputs:

  DO k = 1, icbsmax2                                      ! convect3
    DO i = 1, len                                         ! convect3
      tp(i, k) = 0.0                                      ! convect3
      tvp(i, k) = 0.0                                     ! convect3
      clw(i, k) = 0.0                                     ! convect3
    END DO                                                ! convect3
  END DO                                                  ! convect3

! tp and tvp below cloud base:

  DO k = minorig, icbsmax2 - 1
    DO i = 1, len
      tp(i, k) = tnk(i) - (gz(i,k)-gznk(i))*cpinv(i)
      tvp(i, k) = tp(i, k)*(1.+qnk(i)/eps-qnk(i))        !whole thing (convect3)
    END DO
  END DO

! ***  Find lifted parcel quantities above cloud base    ***

  DO i = 1, len
    tg = ticb(i)
! ori         qg=qs(i,icb(i))
    qg = qsicb(i) ! convect3
! debug         alv=lv0-clmcpv*(ticb(i)-t0)
    alv = lv0 - clmcpv*(ticb(i)-273.15)

! First iteration.

! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
    s = cpd*(1.-qnk(i)) + cl*qnk(i) + &                   ! convect3
        alv*alv*qg/(rrv*ticb(i)*ticb(i))                  ! convect3
    s = 1./s
! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
    ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gzicb(i) ! convect3
    tg = tg + s*(ah0(i)-ahg)
! ori          tg=max(tg,35.0)
! debug          tc=tg-t0
    tc = tg - 273.15
    denom = 243.5 + tc
    denom = max(denom, 1.0) ! convect3
! ori          if(tc.ge.0.0)then
    es = 6.112*exp(17.67*tc/denom)
! ori          else
! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori          endif
! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
    qg = eps*es/(p(i,icbs(i))-es*(1.-eps))

! Second iteration.


! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
! ori          s=1./s
! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
    ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gzicb(i) ! convect3
    tg = tg + s*(ah0(i)-ahg)
! ori          tg=max(tg,35.0)
! debug          tc=tg-t0
    tc = tg - 273.15
    denom = 243.5 + tc
    denom = max(denom, 1.0)                               ! convect3
! ori          if(tc.ge.0.0)then
    es = 6.112*exp(17.67*tc/denom)
! ori          else
! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori          end if
! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
    qg = eps*es/(p(i,icbs(i))-es*(1.-eps))

    alv = lv0 - clmcpv*(ticb(i)-273.15)

! ori c approximation here:
! ori         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i)
! ori     &   -gz(i,icb(i))-alv*qg)/cpd

! convect3: no approximation:
    tp(i, icbs(i)) = (ah0(i)-gz(i,icbs(i))-alv*qg)/(cpd+(cl-cpd)*qnk(i))

! ori         clw(i,icb(i))=qnk(i)-qg
! ori         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
    clw(i, icbs(i)) = qnk(i) - qg
    clw(i, icbs(i)) = max(0.0, clw(i,icbs(i)))

    rg = qg/(1.-qnk(i))
! ori         tvp(i,icb(i))=tp(i,icb(i))*(1.+rg*epsi)
! convect3: (qg utilise au lieu du vrai mixing ratio rg)
    tvp(i, icbs(i)) = tp(i, icbs(i))*(1.+qg/eps-qnk(i))   !whole thing

  END DO

! ori      do 380 k=minorig,icbsmax2
! ori       do 370 i=1,len
! ori         tvp(i,k)=tvp(i,k)-tp(i,k)*qnk(i)
! ori 370   continue
! ori 380  continue


! -- The following is only for convect3:

! * icbs is the first level above the LCL:
! if plcl<p(icb), then icbs=icb+1
! if plcl>p(icb), then icbs=icb

! * the routine above computes tvp from minorig to icbs (included).

! * to compute buoybase (in cv3_trigger.F), both tvp(icb) and tvp(icb+1)
! must be known. This is the case if icbs=icb+1, but not if icbs=icb.

! * therefore, in the case icbs=icb, we compute tvp at level icb+1
! (tvp at other levels will be computed in cv3_undilute2.F)


  DO i = 1, len
    ticb(i) = t(i, icb(i)+1)
    gzicb(i) = gz(i, icb(i)+1)
    qsicb(i) = qs(i, icb(i)+1)
  END DO

  DO i = 1, len
    tg = ticb(i)
    qg = qsicb(i) ! convect3
! debug         alv=lv0-clmcpv*(ticb(i)-t0)
    alv = lv0 - clmcpv*(ticb(i)-273.15)

! First iteration.

! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
    s = cpd*(1.-qnk(i)) + cl*qnk(i) &                         ! convect3
      +alv*alv*qg/(rrv*ticb(i)*ticb(i))                       ! convect3
    s = 1./s
! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
    ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gzicb(i)     ! convect3
    tg = tg + s*(ah0(i)-ahg)
! ori          tg=max(tg,35.0)
! debug          tc=tg-t0
    tc = tg - 273.15
    denom = 243.5 + tc
    denom = max(denom, 1.0)                                   ! convect3
! ori          if(tc.ge.0.0)then
    es = 6.112*exp(17.67*tc/denom)
! ori          else
! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori          endif
! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
    qg = eps*es/(p(i,icb(i)+1)-es*(1.-eps))

! Second iteration.


! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
! ori          s=1./s
! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
    ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gzicb(i)     ! convect3
    tg = tg + s*(ah0(i)-ahg)
! ori          tg=max(tg,35.0)
! debug          tc=tg-t0
    tc = tg - 273.15
    denom = 243.5 + tc
    denom = max(denom, 1.0)                                   ! convect3
! ori          if(tc.ge.0.0)then
    es = 6.112*exp(17.67*tc/denom)
! ori          else
! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori          end if
! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
    qg = eps*es/(p(i,icb(i)+1)-es*(1.-eps))

    alv = lv0 - clmcpv*(ticb(i)-273.15)

! ori c approximation here:
! ori         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i)
! ori     &   -gz(i,icb(i))-alv*qg)/cpd

! convect3: no approximation:
    tp(i, icb(i)+1) = (ah0(i)-gz(i,icb(i)+1)-alv*qg)/(cpd+(cl-cpd)*qnk(i))

! ori         clw(i,icb(i))=qnk(i)-qg
! ori         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
    clw(i, icb(i)+1) = qnk(i) - qg
    clw(i, icb(i)+1) = max(0.0, clw(i,icb(i)+1))

    rg = qg/(1.-qnk(i))
! ori         tvp(i,icb(i))=tp(i,icb(i))*(1.+rg*epsi)
! convect3: (qg utilise au lieu du vrai mixing ratio rg)
    tvp(i, icb(i)+1) = tp(i, icb(i)+1)*(1.+qg/eps-qnk(i))     !whole thing

  END DO

  RETURN
END SUBROUTINE cv3_undilute1

SUBROUTINE cv3_trigger(len, nd, icb, plcl, p, th, tv, tvp, thnk, &
                       pbase, buoybase, iflag, sig, w0)
  IMPLICIT NONE

! -------------------------------------------------------------------
! --- TRIGGERING

! - computes the cloud base
! - triggering (crude in this version)
! - relaxation of sig and w0 when no convection

! Caution1: if no convection, we set iflag=14
! (it used to be 0 in convect3)

! Caution2: at this stage, tvp (and thus buoy) are know up
! through icb only!
! -> the buoyancy below cloud base not (yet) set to the cloud base buoyancy
! -------------------------------------------------------------------

  include "cv3param.h"

! input:
  INTEGER len, nd
  INTEGER icb(len)
  REAL plcl(len), p(len, nd)
  REAL th(len, nd), tv(len, nd), tvp(len, nd)
  REAL thnk(len)

! output:
  REAL pbase(len), buoybase(len)

! input AND output:
  INTEGER iflag(len)
  REAL sig(len, nd), w0(len, nd)

! local variables:
  INTEGER i, k
  REAL tvpbase, tvbase, tdif, ath, ath1


! ***   set cloud base buoyancy at (plcl+dpbase) level buoyancy

  DO i = 1, len
    pbase(i) = plcl(i) + dpbase
    tvpbase = tvp(i, icb(i))  *(pbase(i)-p(i,icb(i)+1))/(p(i,icb(i))-p(i,icb(i)+1)) + &
              tvp(i, icb(i)+1)*(p(i,icb(i))-pbase(i))  /(p(i,icb(i))-p(i,icb(i)+1))
    tvbase = tv(i, icb(i))  *(pbase(i)-p(i,icb(i)+1))/(p(i,icb(i))-p(i,icb(i)+1)) + &
             tv(i, icb(i)+1)*(p(i,icb(i))-pbase(i))  /(p(i,icb(i))-p(i,icb(i)+1))
    buoybase(i) = tvpbase - tvbase
  END DO


! ***   make sure that column is dry adiabatic between the surface  ***
! ***    and cloud base, and that lifted air is positively buoyant  ***
! ***                         at cloud base                         ***
! ***       if not, return to calling program after resetting       ***
! ***                        sig(i) and w0(i)                       ***


! oct3      do 200 i=1,len
! oct3
! oct3       tdif = buoybase(i)
! oct3       ath1 = th(i,1)
! oct3       ath  = th(i,icb(i)-1) - dttrig
! oct3
! oct3       if (tdif.lt.dtcrit .or. ath.gt.ath1) then
! oct3         do 60 k=1,nl
! oct3            sig(i,k) = beta*sig(i,k) - 2.*alpha*tdif*tdif
! oct3            sig(i,k) = AMAX1(sig(i,k),0.0)
! oct3            w0(i,k)  = beta*w0(i,k)
! oct3   60    continue
! oct3         iflag(i)=4 ! pour version vectorisee
! oct3c convect3         iflag(i)=0
! oct3cccc         return
! oct3       endif
! oct3
! oct3200   continue

! -- oct3: on reecrit la boucle 200 (pour la vectorisation)

  DO k = 1, nl
    DO i = 1, len

      tdif = buoybase(i)
      ath1 = thnk(i)
      ath = th(i, icb(i)-1) - dttrig

      IF (tdif<dtcrit .OR. ath>ath1) THEN
        sig(i, k) = beta*sig(i, k) - 2.*alpha*tdif*tdif
        sig(i, k) = amax1(sig(i,k), 0.0)
        w0(i, k) = beta*w0(i, k)
        iflag(i) = 14 ! pour version vectorisee
! convect3         iflag(i)=0
      END IF

    END DO
  END DO

! fin oct3 --

  RETURN
END SUBROUTINE cv3_trigger

SUBROUTINE cv3_compress(len, nloc, ncum, nd, ntra, &
                        iflag1, nk1, icb1, icbs1, &
                        plcl1, tnk1, qnk1, gznk1, pbase1, buoybase1, &
                        t1, q1, qs1, u1, v1, gz1, th1, &
                        tra1, &
                        h1, lv1, cpn1, p1, ph1, tv1, tp1, tvp1, clw1, &
                        sig1, w01, &
                        iflag, nk, icb, icbs, &
                        plcl, tnk, qnk, gznk, pbase, buoybase, &
                        t, q, qs, u, v, gz, th, &
                        tra, &
                        h, lv, cpn, p, ph, tv, tp, tvp, clw, &
                        sig, w0)
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE

  include "cv3param.h"

!inputs:
  INTEGER len, ncum, nd, ntra, nloc
  INTEGER iflag1(len), nk1(len), icb1(len), icbs1(len)
  REAL plcl1(len), tnk1(len), qnk1(len), gznk1(len)
  REAL pbase1(len), buoybase1(len)
  REAL t1(len, nd), q1(len, nd), qs1(len, nd), u1(len, nd), v1(len, nd)
  REAL gz1(len, nd), h1(len, nd), lv1(len, nd), cpn1(len, nd)
  REAL p1(len, nd), ph1(len, nd+1), tv1(len, nd), tp1(len, nd)
  REAL tvp1(len, nd), clw1(len, nd)
  REAL th1(len, nd)
  REAL sig1(len, nd), w01(len, nd)
  REAL tra1(len, nd, ntra)

!outputs:
! en fait, on a nloc=len pour l'instant (cf cv_driver)
  INTEGER iflag(nloc), nk(nloc), icb(nloc), icbs(nloc)
  REAL plcl(nloc), tnk(nloc), qnk(nloc), gznk(nloc)
  REAL pbase(nloc), buoybase(nloc)
  REAL t(nloc, nd), q(nloc, nd), qs(nloc, nd), u(nloc, nd), v(nloc, nd)
  REAL gz(nloc, nd), h(nloc, nd), lv(nloc, nd), cpn(nloc, nd)
  REAL p(nloc, nd), ph(nloc, nd+1), tv(nloc, nd), tp(nloc, nd)
  REAL tvp(nloc, nd), clw(nloc, nd)
  REAL th(nloc, nd)
  REAL sig(nloc, nd), w0(nloc, nd)
  REAL tra(nloc, nd, ntra)

!local variables:
  INTEGER i, k, nn, j

  CHARACTER (LEN=20) :: modname = 'cv3_compress'
  CHARACTER (LEN=80) :: abort_message

  DO k = 1, nl + 1
    nn = 0
    DO i = 1, len
      IF (iflag1(i)==0) THEN
        nn = nn + 1
        sig(nn, k) = sig1(i, k)
        w0(nn, k) = w01(i, k)
        t(nn, k) = t1(i, k)
        q(nn, k) = q1(i, k)
        qs(nn, k) = qs1(i, k)
        u(nn, k) = u1(i, k)
        v(nn, k) = v1(i, k)
        gz(nn, k) = gz1(i, k)
        h(nn, k) = h1(i, k)
        lv(nn, k) = lv1(i, k)
        cpn(nn, k) = cpn1(i, k)
        p(nn, k) = p1(i, k)
        ph(nn, k) = ph1(i, k)
        tv(nn, k) = tv1(i, k)
        tp(nn, k) = tp1(i, k)
        tvp(nn, k) = tvp1(i, k)
        clw(nn, k) = clw1(i, k)
        th(nn, k) = th1(i, k)
      END IF
    END DO
  END DO

!AC!      do 121 j=1,ntra
!AC!ccccc      do 111 k=1,nl+1
!AC!      do 111 k=1,nd
!AC!       nn=0
!AC!      do 101 i=1,len
!AC!      if(iflag1(i).eq.0)then
!AC!       nn=nn+1
!AC!       tra(nn,k,j)=tra1(i,k,j)
!AC!      endif
!AC! 101  continue
!AC! 111  continue
!AC! 121  continue

  IF (nn/=ncum) THEN
    WRITE (lunout, *) 'strange! nn not equal to ncum: ', nn, ncum
    abort_message = ''
    CALL abort_physic(modname, abort_message, 1)
  END IF

  nn = 0
  DO i = 1, len
    IF (iflag1(i)==0) THEN
      nn = nn + 1
      pbase(nn) = pbase1(i)
      buoybase(nn) = buoybase1(i)
      plcl(nn) = plcl1(i)
      tnk(nn) = tnk1(i)
      qnk(nn) = qnk1(i)
      gznk(nn) = gznk1(i)
      nk(nn) = nk1(i)
      icb(nn) = icb1(i)
      icbs(nn) = icbs1(i)
      iflag(nn) = iflag1(i)
    END IF
  END DO

  RETURN
END SUBROUTINE cv3_compress

SUBROUTINE icefrac(t, clw, qi, nl, len)
  IMPLICIT NONE


!JAM--------------------------------------------------------------------
! Calcul de la quantité d'eau sous forme de glace
! --------------------------------------------------------------------
  INTEGER nl, len
  REAL qi(len, nl)
  REAL t(len, nl), clw(len, nl)
  REAL fracg
  INTEGER k, i

  DO k = 3, nl
    DO i = 1, len
      IF (t(i,k)>263.15) THEN
        qi(i, k) = 0.
      ELSE
        IF (t(i,k)<243.15) THEN
          qi(i, k) = clw(i, k)
        ELSE
          fracg = (263.15-t(i,k))/20
          qi(i, k) = clw(i, k)*fracg
        END IF
      END IF
! print*,t(i,k),qi(i,k),'temp,testglace'
    END DO
  END DO

  RETURN

END SUBROUTINE icefrac

SUBROUTINE cv3_undilute2(nloc, ncum, nd, iflag, icb, icbs, nk, &
                         tnk, qnk, gznk, hnk, t, q, qs, gz, &
                         p, ph, h, tv, lv, lf, pbase, buoybase, plcl, &
                         inb, tp, tvp, clw, hp, ep, sigp, buoy, &
                         frac_a, frac_s, qpreca, qta)
  USE print_control_mod, ONLY: prt_level
  IMPLICIT NONE

! ---------------------------------------------------------------------
! Purpose:
! FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
! &
! COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
! FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
! &
! FIND THE LEVEL OF NEUTRAL BUOYANCY

! Main differences convect3/convect4:
!   - icbs (input) is the first level above LCL (may differ from icb)
!   - many minor differences in the iterations
!   - condensed water not removed from tvp in convect3
!   - vertical profile of buoyancy computed here (use of buoybase)
!   - the determination of inb is different
!   - no inb1, only inb in output
! ---------------------------------------------------------------------

  include "cvthermo.h"
  include "cv3param.h"
  include "conema3.h"
  include "cvflag.h"
  include "YOMCST2.h"

!inputs:
  INTEGER, INTENT (IN)                               :: ncum, nd, nloc
  INTEGER, DIMENSION (nloc), INTENT (IN)             :: icb, icbs, nk
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: t, q, qs, gz
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: p
  REAL, DIMENSION (nloc, nd+1), INTENT (IN)          :: ph
  REAL, DIMENSION (nloc), INTENT (IN)                :: tnk, qnk, gznk
  REAL, DIMENSION (nloc), INTENT (IN)                :: hnk
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: lv, lf, tv, h
  REAL, DIMENSION (nloc), INTENT (IN)                :: pbase, buoybase, plcl

!input/outputs:
  REAL, DIMENSION (nloc, nd), INTENT (INOUT)         :: tp, tvp, clw   ! Input for k = 1, icb+1 (computed in cv3_undilute1)
                                                                       ! Output above
  INTEGER, DIMENSION (nloc), INTENT (INOUT)          :: iflag

!outputs:
  INTEGER, DIMENSION (nloc), INTENT (OUT)            :: inb
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: ep, sigp, hp
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: buoy
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: frac_a, frac_s
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: qpreca
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: qta

!local variables:
  INTEGER i, j, k
  REAL smallestreal
  REAL tg, qg, dqgdT, ahg, alv, alf, s, tc, es, esi, denom, rg, tca, elacrit
  REAL                                               :: phinu2p
  REAL                                               :: qhthreshold
  REAL                                               :: als
  REAL                                               :: qsat_new, snew
  REAL, DIMENSION (nloc,nd)                          :: qi
  REAL, DIMENSION (nloc,nd)                          :: ha    ! moist static energy of adiabatic ascents 
                                                              ! taking into account precip ejection
  REAL, DIMENSION (nloc,nd)                          :: hla   ! liquid water static energy of adiabatic ascents 
                                                              ! taking into account precip ejection
  REAL, DIMENSION (nloc,nd)                          :: qcld  ! specific cloud water
  REAL, DIMENSION (nloc,nd)                          :: qhsat    ! specific humidity at saturation
  REAL, DIMENSION (nloc,nd)                          :: dqhsatdT ! dqhsat/dT
  REAL, DIMENSION (nloc,nd)                          :: frac  ! ice fraction function of envt temperature
  REAL, DIMENSION (nloc,nd)                          :: qps   ! specific solid precipitation
  REAL, DIMENSION (nloc,nd)                          :: qpl   ! specific liquid precipitation
  REAL, DIMENSION (nloc)                             :: ah0, cape, capem, byp
  LOGICAL, DIMENSION (nloc)                          :: lcape
  INTEGER, DIMENSION (nloc)                          :: iposit
  REAL                                               :: denomm1
  REAL                                               :: by, defrac, pden, tbis
  REAL                                               :: fracg
  REAL                                               :: deltap
  REAL, SAVE                                         :: Tx, Tm
  DATA Tx/263.15/, Tm/243.15/
!$OMP THREADPRIVATE(Tx, Tm)
  REAL                                               :: aa, bb, dd, ddelta, discr
  REAL                                               :: ff, fp
  REAL                                               :: coefx, coefm, Zx, Zm, Ux, U, Um

  IF (prt_level >= 10) THEN
    print *,'cv3_undilute2.0. icvflag_Tpa, t(1,k), q(1,k), qs(1,k) ', &
                        icvflag_Tpa, (k, t(1,k), q(1,k), qs(1,k), k = 1,nl)
  ENDIF
  smallestreal=tiny(smallestreal)

! =====================================================================
! --- SOME INITIALIZATIONS
! =====================================================================

  DO k = 1, nl
    DO i = 1, ncum
      qi(i, k) = 0.
    END DO
  END DO


! =====================================================================
! --- FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
! =====================================================================

! ---       The procedure is to solve the equation.
!                cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.

! ***  Calculate certain parcel quantities, including static energy   ***


  DO i = 1, ncum
    ah0(i) = (cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i)+ & 
! debug          qnk(i)*(lv0-clmcpv*(tnk(i)-t0))+gznk(i)
             qnk(i)*(lv0-clmcpv*(tnk(i)-273.15)) + gznk(i)
  END DO
!
!  Ice fraction
!
  IF (cvflag_ice) THEN
    DO k = minorig, nl
      DO i = 1, ncum
          frac(i, k) = (Tx - t(i,k))/(Tx - Tm)
          frac(i, k) = min(max(frac(i,k),0.0), 1.0)
      END DO
    END DO
! Below cloud base, set ice fraction to cloud base value
    DO k = 1, nl
      DO i = 1, ncum
        IF (k<icb(i)) THEN
          frac(i,k) = frac(i,icb(i))
        END IF
      END DO
    END DO
  ELSE
    DO k = 1, nl
      DO i = 1, ncum
          frac(i,k) = 0.
      END DO
    END DO
  ENDIF ! (cvflag_ice)


  DO k = minorig, nl
    DO i = 1,ncum
      ha(i,k) = ah0(i)
      hla(i,k) = hnk(i)
      qta(i,k) = qnk(i)
      qpreca(i,k) = 0.
      frac_a(i,k) = 0.
      frac_s(i,k) = frac(i,k)
      qpl(i,k) = 0.
      qps(i,k) = 0.
      qhsat(i,k) = qs(i,k)
      qcld(i,k) = max(qta(i,k)-qhsat(i,k),0.)
      IF (k <= icb(i)+1) THEN
        qhsat(i,k) = qnk(i)-clw(i,k)
        qcld(i,k) = clw(i,k)
      ENDIF 
    ENDDO
  ENDDO

!jyg<
! =====================================================================
! --- SET THE THE FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
! =====================================================================
  DO k = 1, nl
    DO i = 1, ncum
      ep(i, k) = 0.0
      sigp(i, k) = spfac
    END DO
  END DO
!>jyg
!

! ***  Find lifted parcel quantities above cloud base    ***

!----------------------------------------------------------------------------
!
  IF (icvflag_Tpa == 2) THEN
!
!----------------------------------------------------------------------------
!
    DO k = minorig + 1, nl
      DO i = 1,ncum
        tp(i,k) = t(i,k)
      ENDDO
!!      alv = lv0 - clmcpv*(t(i,k)-273.15)
!!      alf = lf0 + clmci*(t(i,k)-273.15)
!!      als = alf + alv
      DO j = 1,4
        DO i = 1, ncum
! ori	    if(k.ge.(icb(i)+1))then
          IF (k>=(icbs(i)+1)) THEN                                ! convect3
            tg = tp(i, k)
            IF (tg .gt. Tx) THEN
              es = 6.112*exp(17.67*(tg - 273.15)/(tg + 243.5 - 273.15))
              qg = eps*es/(p(i,k)-es*(1.-eps))
            ELSE
              esi = exp(23.33086-(6111.72784/tg)+0.15215*log(tg))
              qg = eps*esi/(p(i,k)-esi*(1.-eps))
            ENDIF
! Ice fraction
            ff = 0.
            fp = 1./(Tx - Tm)
            IF (tg < Tx) THEN
              IF (tg > Tm) THEN
                ff = (Tx - tg)*fp
              ELSE
                ff = 1.
              ENDIF ! (tg > Tm)
            ENDIF ! (tg < Tx)
! Intermediate variables
            aa = cpd + (cl-cpd)*qnk(i) + lv(i,k)*lv(i,k)*qg/(rrv*tg*tg)
            ahg = (cpd + (cl-cpd)*qnk(i))*tg + lv(i,k)*qg - &
                  lf(i,k)*ff*(qnk(i) - qg) + gz(i,k)
            dd = lf(i,k)*lv(i,k)*qg/(rrv*tg*tg)
            ddelta = lf(i,k)*(qnk(i) - qg)
            bb = aa + ddelta*fp + dd*fp*(Tx-tg)
! Compute Zx and Zm
            coefx = aa
            coefm = aa + dd
            IF (tg .gt. Tx) THEN
              Zx = ahg            + coefx*(Tx - tg)
              Zm = ahg - ddelta   + coefm*(Tm - tg)
            ELSE
              IF (tg .gt. Tm) THEN
                Zx = ahg          + (coefx +fp*ddelta)*(Tx - Tg)
                Zm = ahg          + (coefm +fp*ddelta)*(Tm - Tg)
              ELSE
                Zx = ahg + ddelta + coefx*(Tx - tg)
                Zm = ahg          + coefm*(Tm - tg)
              ENDIF ! (tg .gt. Tm)
            ENDIF ! (tg .gt. Tx)
! Compute the masks Um, U, Ux
            Um = (sign(1., Zm-ah0(i))+1.)/2.
            Ux = (sign(1., ah0(i)-Zx)+1.)/2.
            U = (1. - Um)*(1. - Ux)
! Compute the updated parcell temperature Tp : 3 cases depending on tg value
            IF (tg .gt. Tx) THEN
              discr = bb*bb - 4*dd*fp*(ah0(i) - ahg + ddelta*fp*(Tx-tg))
              Tp(i,k) = tg + &
                  Um*  (ah0(i) - ahg + ddelta)           /(aa + dd) + &
                  U *2*(ah0(i) - ahg + ddelta*fp*(Tx-tg))/(bb + sqrt(discr)) + &
                  Ux*  (ah0(i) - ahg)                    /aa
            ELSEIF (tg .gt. Tm) THEN
              discr = bb*bb - 4*dd*fp*(ah0(i) - ahg)
              Tp(i,k) = tg + &
                  Um*  (ah0(i) - ahg + ddelta*fp*(tg-Tm))/(aa + dd) + &
                  U *2*(ah0(i) - ahg)                    /(bb + sqrt(discr)) + &
                  Ux*  (ah0(i) - ahg + ddelta*fp*(tg-Tx))/aa
            ELSE
              discr = bb*bb - 4*dd*fp*(ah0(i) - ahg + ddelta*fp*(Tm-tg))
              Tp(i,k) = tg + &
                  Um*  (ah0(i) - ahg)                    /(aa + dd) + &
                  U *2*(ah0(i) - ahg + ddelta*fp*(Tm-tg))/(bb + sqrt(discr)) + &
                  Ux*  (ah0(i) - ahg - ddelta)           /aa
            ENDIF ! (tg .gt. Tx)
!
!!     print *,' j, k, Um, U, Ux, aa, bb, discr, dd, ddelta ', j, k, Um, U, Ux, aa, bb, discr, dd, ddelta
!!     print *,' j, k, ah0(i), ahg, tg, qg, tp(i,k), ff ', j, k, ah0(i), ahg, tg, qg, tp(i,k), ff
          END IF ! (k>=(icbs(i)+1))
        END DO ! i = 1, ncum
      END DO ! j = 1,4
      DO i = 1, ncum
        IF (k>=(icbs(i)+1)) THEN                                ! convect3
          tg = tp(i, k)
          IF (tg .gt. Tx) THEN
            es = 6.112*exp(17.67*(tg - 273.15)/(tg + 243.5 - 273.15))
            qg = eps*es/(p(i,k)-es*(1.-eps))
          ELSE
            esi = exp(23.33086-(6111.72784/tg)+0.15215*log(tg))
            qg = eps*esi/(p(i,k)-esi*(1.-eps))
          ENDIF
          clw(i, k) = qnk(i) - qg
          clw(i, k) = max(0.0, clw(i,k))
          tvp(i, k) = max(0., tp(i,k)*(1.+qg/eps-qnk(i)))
! print*,tvp(i,k),'tvp'
          IF (clw(i,k)<1.E-11) THEN
            tp(i, k) = tv(i, k)
            tvp(i, k) = tv(i, k)
          END IF ! (clw(i,k)<1.E-11)
        END IF ! (k>=(icbs(i)+1))
      END DO ! i = 1, ncum
    END DO ! k = minorig + 1, nl
!----------------------------------------------------------------------------
!
  ELSE IF (icvflag_Tpa == 1) THEN  ! (icvflag_Tpa == 2)
!
!----------------------------------------------------------------------------
!
    DO k = minorig + 1, nl
      DO i = 1,ncum
        tp(i,k) = t(i,k)
      ENDDO
!!      alv = lv0 - clmcpv*(t(i,k)-273.15)
!!      alf = lf0 + clmci*(t(i,k)-273.15)
!!      als = alf + alv
      DO j = 1,4
        DO i = 1, ncum
! ori	    if(k.ge.(icb(i)+1))then
          IF (k>=(icbs(i)+1)) THEN                                ! convect3
            tg = tp(i, k)
            IF (tg .gt. Tx .OR. .NOT.cvflag_ice) THEN
              es = 6.112*exp(17.67*(tg - 273.15)/(tg + 243.5 - 273.15))
              qg = eps*es/(p(i,k)-es*(1.-eps))
              dqgdT = lv(i,k)*qg/(rrv*tg*tg)
            ELSE
              esi = exp(23.33086-(6111.72784/tg)+0.15215*log(tg))
              qg = eps*esi/(p(i,k)-esi*(1.-eps))
              dqgdT = (lv(i,k)+lf(i,k))*qg/(rrv*tg*tg)
            ENDIF
            IF (qsat_depends_on_qt) THEN
              dqgdT = dqgdT*(1.-qta(i,k-1))/(1.-qg)**2
              qg = qg*(1.-qta(i,k-1))/(1.-qg)            
            ENDIF
            ahg = (cpd + (cl-cpd)*qta(i,k-1))*tg + lv(i,k)*qg - &
                  lf(i,k)*frac(i,k)*(qta(i,k-1) - qg) + gz(i,k)
            Tp(i,k) = tg + (ah0(i) - ahg)/ &
                    (cpd + (cl-cpd)*qta(i,k-1) + (lv(i,k)+frac(i,k)*lf(i,k))*dqgdT)
!!   print *,'undilute2 iterations k, Tp(i,k), ah0(i), ahg ', &
!!                                 k, Tp(i,k), ah0(i), ahg
          END IF ! (k>=(icbs(i)+1))
        END DO ! i = 1, ncum
      END DO ! j = 1,4
      DO i = 1, ncum
        IF (k>=(icbs(i)+1)) THEN                                ! convect3
          tg = tp(i, k)
          IF (tg .gt. Tx .OR. .NOT.cvflag_ice) THEN
            es = 6.112*exp(17.67*(tg - 273.15)/(tg + 243.5 - 273.15))
            qg = eps*es/(p(i,k)-es*(1.-eps))
          ELSE
            esi = exp(23.33086-(6111.72784/tg)+0.15215*log(tg))
            qg = eps*esi/(p(i,k)-esi*(1.-eps))
          ENDIF
          IF (qsat_depends_on_qt) THEN
            qg = qg*(1.-qta(i,k-1))/(1.-qg)            
          ENDIF
          qhsat(i,k) = qg
        END IF ! (k>=(icbs(i)+1))
      END DO ! i = 1, ncum
      DO i = 1, ncum
        IF (k>=(icbs(i)+1)) THEN                                ! convect3
          clw(i, k) = qta(i,k-1) - qhsat(i,k)
          clw(i, k) = max(0.0, clw(i,k))
          tvp(i, k) = max(0., tp(i,k)*(1.+qhsat(i,k)/eps-qta(i,k-1)))
! print*,tvp(i,k),'tvp'
          IF (clw(i,k)<1.E-11) THEN
            tp(i, k) = tv(i, k)
            tvp(i, k) = tv(i, k)
          END IF ! (clw(i,k)<1.E-11)
        END IF ! (k>=(icbs(i)+1))
      END DO ! i = 1, ncum
!
      IF (cvflag_prec_eject) THEN
        DO i = 1, ncum
          IF (k>=(icbs(i)+1)) THEN                                ! convect3
!  Specific precipitation (liquid and solid) and ice content 
!  before ejection of precipitation                                                     !!jygprl
            elacrit = elcrit*min(max(1.-(tp(i,k)-T0)/Tlcrit, 0.), 1.)                   !!jygprl
!!!!            qcld(i,k) = min(clw(i,k), elacrit)                                          !!jygprl
            qhthreshold = elacrit*(1.-qta(i,k-1))/(1.-elacrit)
            qcld(i,k) = min(clw(i,k), qhthreshold)             !!jygprl
!!!!            phinu2p = max(qhsat(i,k-1) + qcld(i,k-1) - (qhsat(i,k) + qcld(i,k)),0.)   !!jygprl
            phinu2p = max(clw(i,k) - max(qta(i,k-1) - qhsat(i,k-1), qhthreshold), 0.)
            qpl(i,k) = qpl(i,k-1) + (1.-frac(i,k))*phinu2p                            !!jygprl
            qps(i,k) = qps(i,k-1) + frac(i,k)     *phinu2p                            !!jygprl
            qi(i,k) = (1.-ejectliq)*clw(i,k)*frac(i,k) + &                            !!jygprl
                     ejectliq*(qps(i,k-1) + frac(i,k)*(phinu2p+qcld(i,k)))            !!jygprl
!!
!  =====================================================================================
!  Ejection of precipitation from adiabatic ascents if requested (cvflag_prec_eject=True):
!  Compute the steps of total water (qta), of moist static energy (ha), of specific 
!  precipitation (qpl and qps) and of specific cloud water (qcld) associated with precipitation
!   ejection.
!  =====================================================================================
!  
!   Verif
            qpreca(i,k) = ejectliq*qpl(i,k) + ejectice*qps(i,k)                                   !!jygprl
            frac_a(i,k) = ejectice*qps(i,k)/max(qpreca(i,k),smallestreal)                         !!jygprl
            frac_s(i,k) = (1.-ejectliq)*frac(i,k) + &                                             !!jygprl
               ejectliq*(1. - (qpl(i,k)+(1.-frac(i,k))*qcld(i,k))/max(clw(i,k),smallestreal))     !!jygprl
!          
            denomm1 = 1./(1. - qpreca(i,k))
!          
            qta(i,k) = qta(i,k-1) - &
                      qpreca(i,k)*(1.-qta(i,k-1))*denomm1
            ha(i,k)  = ha(i,k-1) + &
                      ( qpreca(i,k)*(-(1.-qta(i,k-1))*(cl-cpd)*tp(i,k) + &
                                  lv(i,k)*qhsat(i,k) - lf(i,k)*(frac_s(i,k)*qcld(i,k)+qps(i,k))) + &
                        lf(i,k)*ejectice*qps(i,k))*denomm1
            hla(i,k) = hla(i,k-1) + &
                      ( qpreca(i,k)*(-(1.-qta(i,k-1))*(cpv-cpd)*tp(i,k) - &
                                  lv(i,k)*((1.-frac_s(i,k))*qcld(i,k)+qpl(i,k)) - &
                                  (lv(i,k)+lf(i,k))*(frac_s(i,k)*qcld(i,k)+qps(i,k))) + &
                        lv(i,k)*ejectliq*qpl(i,k) + (lv(i,k)+lf(i,k))*ejectice*qps(i,k))*denomm1
            qpl(i,k) = qpl(i,k)*(1.-ejectliq)*denomm1
            qps(i,k) = qps(i,k)*(1.-ejectice)*denomm1
            qcld(i,k) = qcld(i,k)*denomm1
            qhsat(i,k) = qhsat(i,k)*(1.-qta(i,k))/(1.-qta(i,k-1))
         END IF ! (k>=(icbs(i)+1))
        END DO ! i = 1, ncum
      ENDIF  ! (cvflag_prec_eject)
! 
    END DO ! k = minorig + 1, nl
!
!----------------------------------------------------------------------------
!
  ELSE IF (icvflag_Tpa == 0) THEN! (icvflag_Tpa == 2) ELSE IF(icvflag_Tpa == 1)
!
!----------------------------------------------------------------------------
!
  DO k = minorig + 1, nl
    DO i = 1, ncum
! ori	    if(k.ge.(icb(i)+1))then
      IF (k>=(icbs(i)+1)) THEN                                ! convect3
        tg = t(i, k)
        qg = qs(i, k)
! debug	      alv=lv0-clmcpv*(t(i,k)-t0)
        alv = lv0 - clmcpv*(t(i,k)-273.15)

! First iteration.

! ori	       s=cpd+alv*alv*qg/(rrv*t(i,k)*t(i,k))
        s = cpd*(1.-qnk(i)) + cl*qnk(i) + &                   ! convect3
            alv*alv*qg/(rrv*t(i,k)*t(i,k))                    ! convect3
        s = 1./s
! ori	       ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
        ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gz(i, k) ! convect3
        tg = tg + s*(ah0(i)-ahg)
! ori	       tg=max(tg,35.0)
! debug	       tc=tg-t0
        tc = tg - 273.15
        denom = 243.5 + tc
        denom = max(denom, 1.0)                               ! convect3
! ori	       if(tc.ge.0.0)then
        es = 6.112*exp(17.67*tc/denom)
! ori	       else
! ori			es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori	       endif
        qg = eps*es/(p(i,k)-es*(1.-eps))

! Second iteration.

! ori	       s=cpd+alv*alv*qg/(rrv*t(i,k)*t(i,k))
! ori	       s=1./s
! ori	       ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
        ahg = cpd*tg + (cl-cpd)*qnk(i)*tg + alv*qg + gz(i, k) ! convect3
        tg = tg + s*(ah0(i)-ahg)
! ori	       tg=max(tg,35.0)
! debug	       tc=tg-t0
        tc = tg - 273.15
        denom = 243.5 + tc
        denom = max(denom, 1.0)                               ! convect3
! ori	       if(tc.ge.0.0)then
        es = 6.112*exp(17.67*tc/denom)
! ori	       else
! ori			es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori	       endif
        qg = eps*es/(p(i,k)-es*(1.-eps))

! debug	       alv=lv0-clmcpv*(t(i,k)-t0)
        alv = lv0 - clmcpv*(t(i,k)-273.15)
! print*,'cpd dans convect2 ',cpd
! print*,'tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd'
! print*,tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd

! ori c approximation here:
! ori        tp(i,k)=(ah0(i)-(cl-cpd)*qnk(i)*t(i,k)-gz(i,k)-alv*qg)/cpd

! convect3: no approximation:
        IF (cvflag_ice) THEN
          tp(i, k) = max(0., (ah0(i)-gz(i,k)-alv*qg)/(cpd+(cl-cpd)*qnk(i)))
        ELSE
          tp(i, k) = (ah0(i)-gz(i,k)-alv*qg)/(cpd+(cl-cpd)*qnk(i))
        END IF

        clw(i, k) = qnk(i) - qg
        clw(i, k) = max(0.0, clw(i,k))
        rg = qg/(1.-qnk(i))
! ori               tvp(i,k)=tp(i,k)*(1.+rg*epsi)
! convect3: (qg utilise au lieu du vrai mixing ratio rg):
        tvp(i, k) = tp(i, k)*(1.+qg/eps-qnk(i)) ! whole thing
        IF (cvflag_ice) THEN
          IF (clw(i,k)<1.E-11) THEN
            tp(i, k) = tv(i, k)
            tvp(i, k) = tv(i, k)
          END IF
        END IF
!jyg<
!!      END IF  ! Endif moved to the end of the loop
!>jyg

      IF (cvflag_ice) THEN
!CR:attention boucle en klon dans Icefrac
! Call Icefrac(t,clw,qi,nl,nloc)
        IF (t(i,k)>263.15) THEN
          qi(i, k) = 0.
        ELSE
          IF (t(i,k)<243.15) THEN
            qi(i, k) = clw(i, k)
          ELSE
            fracg = (263.15-t(i,k))/20
            qi(i, k) = clw(i, k)*fracg
          END IF
        END IF
!CR: fin test
        IF (t(i,k)<263.15) THEN
!CR: on commente les calculs d'Arnaud car division par zero
! nouveau calcul propose par JYG
!       alv=lv0-clmcpv*(t(i,k)-273.15)
!       alf=lf0-clmci*(t(i,k)-273.15)
!       tg=tp(i,k)
!       tc=tp(i,k)-273.15
!       denom=243.5+tc
!       do j=1,3
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! il faudra que esi vienne en argument de la convection
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        tbis=t(i,k)+(tp(i,k)-tg)
!        esi=exp(23.33086-(6111.72784/tbis) + &
!                       0.15215*log(tbis))
!        qsat_new=eps*esi/(p(i,k)-esi*(1.-eps))
!        snew=cpd*(1.-qnk(i))+cl*qnk(i)+alv*alv*qsat_new/ &
!                                       (rrv*tbis*tbis)
!        snew=1./snew
!        print*,esi,qsat_new,snew,'esi,qsat,snew'
!        tp(i,k)=tg+(alf*qi(i,k)+alv*qg*(1.-(esi/es)))*snew
!        print*,k,tp(i,k),qnk(i),'avec glace'
!        print*,'tpNAN',tg,alf,qi(i,k),alv,qg,esi,es,snew
!       enddo

          alv = lv0 - clmcpv*(t(i,k)-273.15)
          alf = lf0 + clmci*(t(i,k)-273.15)
          als = alf + alv
          tg = tp(i, k)
          tp(i, k) = t(i, k)
          DO j = 1, 3
            esi = exp(23.33086-(6111.72784/tp(i,k))+0.15215*log(tp(i,k)))
            qsat_new = eps*esi/(p(i,k)-esi*(1.-eps))
            snew = cpd*(1.-qnk(i)) + cl*qnk(i) + alv*als*qsat_new/ &
                                                 (rrv*tp(i,k)*tp(i,k))
            snew = 1./snew
! c             print*,esi,qsat_new,snew,'esi,qsat,snew'
            tp(i, k) = tp(i, k) + &
                       ((cpd*(1.-qnk(i))+cl*qnk(i))*(tg-tp(i,k)) + &
                        alv*(qg-qsat_new)+alf*qi(i,k))*snew
! print*,k,tp(i,k),qsat_new,qnk(i),qi(i,k), &
!              'k,tp,q,qt,qi avec glace'
          END DO

!CR:reprise du code AJ
          clw(i, k) = qnk(i) - qsat_new
          clw(i, k) = max(0.0, clw(i,k))
          tvp(i, k) = max(0., tp(i,k)*(1.+qsat_new/eps-qnk(i)))
! print*,tvp(i,k),'tvp'
        END IF
        IF (clw(i,k)<1.E-11) THEN
          tp(i, k) = tv(i, k)
          tvp(i, k) = tv(i, k)
        END IF
      END IF ! (cvflag_ice)
!jyg<
      END IF ! (k>=(icbs(i)+1))
!>jyg
    END DO
  END DO

!----------------------------------------------------------------------------
!
  ENDIF ! (icvflag_Tpa == 2) ELSEIF (icvflag_Tpa == 1) ELSE (icvflag_Tpa == 0)
!
!----------------------------------------------------------------------------
!
! =====================================================================
! --- SET THE PRECIPITATION EFFICIENCIES 
! --- THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)
! =====================================================================
!
  IF (flag_epkeorig/=1) THEN
    DO k = 1, nl ! convect3
      DO i = 1, ncum
!jyg<
       IF(k>=icb(i)) THEN
!>jyg
         pden = ptcrit - pbcrit
         ep(i, k) = (plcl(i)-p(i,k)-pbcrit)/pden*epmax
         ep(i, k) = max(ep(i,k), 0.0)
         ep(i, k) = min(ep(i,k), epmax)
!!         sigp(i, k) = spfac  ! jyg
        ENDIF   ! (k>=icb(i))
      END DO
    END DO
  ELSE
    DO k = 1, nl
      DO i = 1, ncum
        IF(k>=icb(i)) THEN
!!        IF (k>=(nk(i)+1)) THEN
!>jyg
          tca = tp(i, k) - t0
          IF (tca>=0.0) THEN
            elacrit = elcrit
          ELSE
            elacrit = elcrit*(1.0-tca/tlcrit)
          END IF
          elacrit = max(elacrit, 0.0)
          ep(i, k) = 1.0 - elacrit/max(clw(i,k), 1.0E-8)
          ep(i, k) = max(ep(i,k), 0.0)
          ep(i, k) = min(ep(i,k), epmax)
!!          sigp(i, k) = spfac  ! jyg
        END IF  ! (k>=icb(i))
      END DO
    END DO
  END IF
!
!   =========================================================================
  IF (prt_level >= 10) THEN
    print *,'cv3_undilute2.1. tp(1,k), tvp(1,k) ', &
                          (k, tp(1,k), tvp(1,k), k = 1,nl)
  ENDIF
!
! =====================================================================
! --- CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL
! --- VIRTUAL TEMPERATURE
! =====================================================================

! dans convect3, tvp est calcule en une seule fois, et sans retirer
! l'eau condensee (~> reversible CAPE)

! ori      do 340 k=minorig+1,nl
! ori        do 330 i=1,ncum
! ori        if(k.ge.(icb(i)+1))then
! ori          tvp(i,k)=tvp(i,k)*(1.0-qnk(i)+ep(i,k)*clw(i,k))
! oric         print*,'i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)'
! oric         print*, i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)
! ori        endif
! ori 330    continue
! ori 340  continue

! ori      do 350 i=1,ncum
! ori       tvp(i,nlp)=tvp(i,nl)-(gz(i,nlp)-gz(i,nl))/cpd
! ori 350  continue

  DO i = 1, ncum                                           ! convect3
    tp(i, nlp) = tp(i, nl)                                 ! convect3
  END DO                                                   ! convect3

! =====================================================================
! --- EFFECTIVE VERTICAL PROFILE OF BUOYANCY (convect3 only):
! =====================================================================

! -- this is for convect3 only:

! first estimate of buoyancy:

!jyg : k-loop outside i-loop (07042015)
  DO k = 1, nl
    DO i = 1, ncum
      buoy(i, k) = tvp(i, k) - tv(i, k)
    END DO
  END DO

! set buoyancy=buoybase for all levels below base
! for safety, set buoy(icb)=buoybase

!jyg : k-loop outside i-loop (07042015)
  DO k = 1, nl
    DO i = 1, ncum
      IF ((k>=icb(i)) .AND. (k<=nl) .AND. (p(i,k)>=pbase(i))) THEN
        buoy(i, k) = buoybase(i)
      END IF
    END DO
  END DO
  DO i = 1, ncum
!    buoy(icb(i),k)=buoybase(i)
    buoy(i, icb(i)) = buoybase(i)
  END DO

! -- end convect3

! =====================================================================
! --- FIND THE FIRST MODEL LEVEL (INB) ABOVE THE PARCEL'S
! --- LEVEL OF NEUTRAL BUOYANCY
! =====================================================================

! -- this is for convect3 only:

  DO i = 1, ncum
    inb(i) = nl - 1
    iposit(i) = nl
  END DO


! --    iposit(i) = first level, above icb, with positive buoyancy
  DO k = 1, nl - 1
    DO i = 1, ncum
      IF (k>=icb(i) .AND. buoy(i,k)>0.) THEN
        iposit(i) = min(iposit(i), k)
      END IF
    END DO
  END DO

  DO i = 1, ncum
    IF (iposit(i)==nl) THEN
      iposit(i) = icb(i)
    END IF
  END DO

  DO k = 1, nl - 1
    DO i = 1, ncum
      IF ((k>=iposit(i)) .AND. (buoy(i,k)<dtovsh)) THEN
        inb(i) = min(inb(i), k)
      END IF
    END DO
  END DO

!CR fix computation of inb
!keep flag or modify in all cases?
  IF (iflag_mix_adiab.eq.1) THEN
  DO i = 1, ncum
     cape(i)=0.
     inb(i)=icb(i)+1
  ENDDO
  
  DO k = 2, nl 
    DO i = 1, ncum
       IF ((k>=iposit(i))) THEN
       deltap = min(plcl(i), ph(i,k-1)) - min(plcl(i), ph(i,k))
       cape(i) = cape(i) + rrd*buoy(i, k-1)*deltap/p(i, k-1)
       IF (cape(i).gt.0.) THEN
        inb(i) = max(inb(i), k)
       END IF
       ENDIF
    ENDDO
  ENDDO

!  DO i = 1, ncum
!     print*,"inb",inb(i)
!  ENDDO

  endif

! -- end convect3

! ori      do 510 i=1,ncum
! ori        cape(i)=0.0
! ori        capem(i)=0.0
! ori        inb(i)=icb(i)+1
! ori        inb1(i)=inb(i)
! ori 510  continue

! Originial Code

!    do 530 k=minorig+1,nl-1
!     do 520 i=1,ncum
!      if(k.ge.(icb(i)+1))then
!       by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
!       byp=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
!       cape(i)=cape(i)+by
!       if(by.ge.0.0)inb1(i)=k+1
!       if(cape(i).gt.0.0)then
!        inb(i)=k+1
!        capem(i)=cape(i)
!       endif
!      endif
!520    continue
!530  continue
!    do 540 i=1,ncum
!     byp=(tvp(i,nl)-tv(i,nl))*dph(i,nl)/p(i,nl)
!     cape(i)=capem(i)+byp
!     defrac=capem(i)-cape(i)
!     defrac=max(defrac,0.001)
!     frac(i)=-cape(i)/defrac
!     frac(i)=min(frac(i),1.0)
!     frac(i)=max(frac(i),0.0)
!540   continue

!    K Emanuel fix

!    call zilch(byp,ncum)
!    do 530 k=minorig+1,nl-1
!     do 520 i=1,ncum
!      if(k.ge.(icb(i)+1))then
!       by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
!       cape(i)=cape(i)+by
!       if(by.ge.0.0)inb1(i)=k+1
!       if(cape(i).gt.0.0)then
!        inb(i)=k+1
!        capem(i)=cape(i)
!        byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
!       endif
!      endif
!520    continue
!530  continue
!    do 540 i=1,ncum
!     inb(i)=max(inb(i),inb1(i))
!     cape(i)=capem(i)+byp(i)
!     defrac=capem(i)-cape(i)
!     defrac=max(defrac,0.001)
!     frac(i)=-cape(i)/defrac
!     frac(i)=min(frac(i),1.0)
!     frac(i)=max(frac(i),0.0)
!540   continue

! J Teixeira fix

! ori      call zilch(byp,ncum)
! ori      do 515 i=1,ncum
! ori        lcape(i)=.true.
! ori 515  continue
! ori      do 530 k=minorig+1,nl-1
! ori        do 520 i=1,ncum
! ori          if(cape(i).lt.0.0)lcape(i)=.false.
! ori          if((k.ge.(icb(i)+1)).and.lcape(i))then
! ori            by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
! ori            byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
! ori            cape(i)=cape(i)+by
! ori            if(by.ge.0.0)inb1(i)=k+1
! ori            if(cape(i).gt.0.0)then
! ori              inb(i)=k+1
! ori              capem(i)=cape(i)
! ori            endif
! ori          endif
! ori 520    continue
! ori 530  continue
! ori      do 540 i=1,ncum
! ori          cape(i)=capem(i)+byp(i)
! ori          defrac=capem(i)-cape(i)
! ori          defrac=max(defrac,0.001)
! ori          frac(i)=-cape(i)/defrac
! ori          frac(i)=min(frac(i),1.0)
! ori          frac(i)=max(frac(i),0.0)
! ori 540  continue

! --------------------------------------------------------------------
!   Prevent convection when top is too hot
! --------------------------------------------------------------------
  DO i = 1,ncum
    IF (t(i,inb(i)) > T_top_max) iflag(i) = 10
  ENDDO

! =====================================================================
! ---   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL
! =====================================================================

  DO k = 1, nl
    DO i = 1, ncum
      hp(i, k) = h(i, k)
    END DO
  END DO

!jyg : cvflag_ice test outside the loops (07042015)
!
  IF (cvflag_ice) THEN
!
  IF (cvflag_prec_eject) THEN
!!    DO k = minorig + 1, nl
!!      DO i = 1, ncum
!!        IF ((k>=icb(i)) .AND. (k<=inb(i))) THEN
!!          frac_s(i,k) = qi(i,k)/max(clw(i,k),smallestreal)    
!!          frac_s(i,k) = 1. - (qpl(i,k)+(1.-frac_s(i,k))*qcld(i,k))/max(clw(i,k),smallestreal)    
!!        END IF
!!      END DO
!!    END DO
  ELSE    ! (cvflag_prec_eject)
    DO k = minorig + 1, nl
      DO i = 1, ncum
        IF ((k>=icb(i)) .AND. (k<=inb(i))) THEN
!jyg< frac computation moved to beginning of cv3_undilute2.
!     kept here for compatibility test with CMip6 version
          frac_s(i, k) = 1. - (t(i,k)-243.15)/(263.15-243.15)
          frac_s(i, k) = min(max(frac_s(i,k),0.0), 1.0)
        END IF
      END DO
    END DO
  ENDIF  ! (cvflag_prec_eject) ELSE
    DO k = minorig + 1, nl
      DO i = 1, ncum
        IF ((k>=icb(i)) .AND. (k<=inb(i))) THEN
!!          hp(i, k) = hnk(i) + (lv(i,k)+(cpd-cpv)*t(i,k)+frac_s(i,k)*lf(i,k))* &     !!jygprl
!!                              ep(i, k)*clw(i, k)                                    !!jygprl
          hp(i, k) = hla(i,k-1) + (lv(i,k)+(cpd-cpv)*t(i,k)+frac_s(i,k)*lf(i,k))* &   !!jygprl
                              ep(i, k)*clw(i, k)                                      !!jygprl
        END IF
      END DO
    END DO
!
  ELSE   ! (cvflag_ice)
!
    DO k = minorig + 1, nl
      DO i = 1, ncum
        IF ((k>=icb(i)) .AND. (k<=inb(i))) THEN
!jyg<   (energy conservation tests)
!!          hp(i, k) = hnk(i) + (lv(i,k)+(cpd-cpv)*tp(i,k))*ep(i, k)*clw(i, k)
!!          hp(i, k) = ( hnk(i) + (lv(i,k)+(cpd-cpv)*t(i,k))*ep(i, k)*clw(i, k) ) / &
!!                     (1. - ep(i,k)*clw(i,k))
!!          hp(i, k) = ( hnk(i) + (lv(i,k)+(cpd-cl)*t(i,k))*ep(i, k)*clw(i, k) ) / &
!!                     (1. - ep(i,k)*clw(i,k))
          hp(i, k) = hnk(i) + (lv(i,k)+(cpd-cpv)*t(i,k))*ep(i, k)*clw(i, k)
        END IF
      END DO
    END DO
!
  END IF  ! (cvflag_ice)

  RETURN
END SUBROUTINE cv3_undilute2

SUBROUTINE cv3_closure(nloc, ncum, nd, icb, inb, &
                       pbase, p, ph, tv, buoy, &
                       sig, w0, cape, m, iflag)
  IMPLICIT NONE

! ===================================================================
! ---  CLOSURE OF CONVECT3
!
! vectorization: S. Bony
! ===================================================================

  include "cvthermo.h"
  include "cv3param.h"

!input:
  INTEGER ncum, nd, nloc
  INTEGER icb(nloc), inb(nloc)
  REAL pbase(nloc)
  REAL p(nloc, nd), ph(nloc, nd+1)
  REAL tv(nloc, nd), buoy(nloc, nd)

!input/output:
  REAL sig(nloc, nd), w0(nloc, nd)
  INTEGER iflag(nloc)

!output:
  REAL cape(nloc)
  REAL m(nloc, nd)

!local variables:
  INTEGER i, j, k, icbmax
  REAL deltap, fac, w, amu
  REAL dtmin(nloc, nd), sigold(nloc, nd)
  REAL cbmflast(nloc)


! -------------------------------------------------------
! -- Initialization
! -------------------------------------------------------

  DO k = 1, nl
    DO i = 1, ncum
      m(i, k) = 0.0
    END DO
  END DO

! -------------------------------------------------------
! -- Reset sig(i) and w0(i) for i>inb and i<icb
! -------------------------------------------------------

! update sig and w0 above LNB:

  DO k = 1, nl - 1
    DO i = 1, ncum
      IF ((inb(i)<(nl-1)) .AND. (k>=(inb(i)+1))) THEN
        sig(i, k) = beta*sig(i, k) + &
                    2.*alpha*buoy(i, inb(i))*abs(buoy(i,inb(i)))
        sig(i, k) = amax1(sig(i,k), 0.0)
        w0(i, k) = beta*w0(i, k)
      END IF
    END DO
  END DO

! compute icbmax:

  icbmax = 2
  DO i = 1, ncum
    icbmax = max(icbmax, icb(i))
  END DO

! update sig and w0 below cloud base:

  DO k = 1, icbmax
    DO i = 1, ncum
      IF (k<=icb(i)) THEN
        sig(i, k) = beta*sig(i, k) - &
                    2.*alpha*buoy(i, icb(i))*buoy(i, icb(i))
        sig(i, k) = max(sig(i,k), 0.0)
        w0(i, k) = beta*w0(i, k)
      END IF
    END DO
  END DO

!!      if(inb.lt.(nl-1))then
!!         do 85 i=inb+1,nl-1
!!            sig(i)=beta*sig(i)+2.*alpha*buoy(inb)*
!!     1              abs(buoy(inb))
!!            sig(i)=max(sig(i),0.0)
!!            w0(i)=beta*w0(i)
!!   85    continue
!!      end if

!!      do 87 i=1,icb
!!         sig(i)=beta*sig(i)-2.*alpha*buoy(icb)*buoy(icb)
!!         sig(i)=max(sig(i),0.0)
!!         w0(i)=beta*w0(i)
!!   87 continue

! -------------------------------------------------------------
! -- Reset fractional areas of updrafts and w0 at initial time
! -- and after 10 time steps of no convection
! -------------------------------------------------------------

  DO k = 1, nl - 1
    DO i = 1, ncum
      IF (sig(i,nd)<1.5 .OR. sig(i,nd)>12.0) THEN
        sig(i, k) = 0.0
        w0(i, k) = 0.0
      END IF
    END DO
  END DO

! -------------------------------------------------------------
! -- Calculate convective available potential energy (cape),
! -- vertical velocity (w), fractional area covered by
! -- undilute updraft (sig), and updraft mass flux (m)
! -------------------------------------------------------------

  DO i = 1, ncum
    cape(i) = 0.0
  END DO

! compute dtmin (minimum buoyancy between ICB and given level k):

  DO i = 1, ncum
    DO k = 1, nl
      dtmin(i, k) = 100.0
    END DO
  END DO

  DO i = 1, ncum
    DO k = 1, nl
      DO j = minorig, nl
        IF ((k>=(icb(i)+1)) .AND. (k<=inb(i)) .AND. (j>=icb(i)) .AND. (j<=(k-1))) THEN
          dtmin(i, k) = amin1(dtmin(i,k), buoy(i,j))
        END IF
      END DO
    END DO
  END DO

! the interval on which cape is computed starts at pbase :

  DO k = 1, nl
    DO i = 1, ncum

      IF ((k>=(icb(i)+1)) .AND. (k<=inb(i))) THEN

        deltap = min(pbase(i), ph(i,k-1)) - min(pbase(i), ph(i,k))
        cape(i) = cape(i) + rrd*buoy(i, k-1)*deltap/p(i, k-1)
        cape(i) = amax1(0.0, cape(i))
        sigold(i, k) = sig(i, k)

! dtmin(i,k)=100.0
! do 97 j=icb(i),k-1 ! mauvaise vectorisation
! dtmin(i,k)=AMIN1(dtmin(i,k),buoy(i,j))
! 97     continue

        sig(i, k) = beta*sig(i, k) + alpha*dtmin(i, k)*abs(dtmin(i,k))
        sig(i, k) = max(sig(i,k), 0.0)
        sig(i, k) = amin1(sig(i,k), 0.01)
        fac = amin1(((dtcrit-dtmin(i,k))/dtcrit), 1.0)
        w = (1.-beta)*fac*sqrt(cape(i)) + beta*w0(i, k)
        amu = 0.5*(sig(i,k)+sigold(i,k))*w
        m(i, k) = amu*0.007*p(i, k)*(ph(i,k)-ph(i,k+1))/tv(i, k)
        w0(i, k) = w
      END IF

    END DO
  END DO

  DO i = 1, ncum
    w0(i, icb(i)) = 0.5*w0(i, icb(i)+1)
    m(i, icb(i)) = 0.5*m(i, icb(i)+1)*(ph(i,icb(i))-ph(i,icb(i)+1))/(ph(i,icb(i)+1)-ph(i,icb(i)+2))
    sig(i, icb(i)) = sig(i, icb(i)+1)
    sig(i, icb(i)-1) = sig(i, icb(i))
  END DO

! ccc 3. Compute final cloud base mass flux and set iflag to 3 if
! ccc    cloud base mass flux is exceedingly small and is decreasing (i.e. if
! ccc    the final mass flux (cbmflast) is greater than the target mass flux
! ccc    (cbmf) ??).
! cc
! c      do i = 1,ncum
! c       cbmflast(i) = 0.
! c      enddo
! cc
! c      do k= 1,nl
! c       do i = 1,ncum
! c        IF (k .ge. icb(i) .and. k .le. inb(i)) THEN
! c         cbmflast(i) = cbmflast(i)+M(i,k)
! c        ENDIF
! c       enddo
! c      enddo
! cc
! c      do i = 1,ncum
! c       IF (cbmflast(i) .lt. 1.e-6) THEN
! c         iflag(i) = 3
! c       ENDIF
! c      enddo
! cc
! c      do k= 1,nl
! c       do i = 1,ncum
! c        IF (iflag(i) .ge. 3) THEN
! c         M(i,k) = 0.
! c         sig(i,k) = 0.
! c         w0(i,k) = 0.
! c        ENDIF
! c       enddo
! c      enddo
! cc
!!      cape=0.0
!!      do 98 i=icb+1,inb
!!         deltap = min(pbase,ph(i-1))-min(pbase,ph(i))
!!         cape=cape+rrd*buoy(i-1)*deltap/p(i-1)
!!         dcape=rrd*buoy(i-1)*deltap/p(i-1)
!!         dlnp=deltap/p(i-1)
!!         cape=max(0.0,cape)
!!         sigold=sig(i)

!!         dtmin=100.0
!!         do 97 j=icb,i-1
!!            dtmin=amin1(dtmin,buoy(j))
!!   97    continue

!!         sig(i)=beta*sig(i)+alpha*dtmin*abs(dtmin)
!!         sig(i)=max(sig(i),0.0)
!!         sig(i)=amin1(sig(i),0.01)
!!         fac=amin1(((dtcrit-dtmin)/dtcrit),1.0)
!!         w=(1.-beta)*fac*sqrt(cape)+beta*w0(i)
!!         amu=0.5*(sig(i)+sigold)*w
!!         m(i)=amu*0.007*p(i)*(ph(i)-ph(i+1))/tv(i)
!!         w0(i)=w
!!   98 continue
!!      w0(icb)=0.5*w0(icb+1)
!!      m(icb)=0.5*m(icb+1)*(ph(icb)-ph(icb+1))/(ph(icb+1)-ph(icb+2))
!!      sig(icb)=sig(icb+1)
!!      sig(icb-1)=sig(icb)

  RETURN
END SUBROUTINE cv3_closure

SUBROUTINE cv3_mixing(nloc, ncum, nd, na, ntra, icb, nk, inb, &
                      ph, t, rr, rs, u, v, tra, h, lv, lf, frac, qnk, &
                      unk, vnk, hp, tv, tvp, ep, clw, m, sig, &
                      ment, qent, uent, vent, nent, sij, elij, ments, qents, traent)
  IMPLICIT NONE

! ---------------------------------------------------------------------
! a faire:
! - vectorisation de la partie normalisation des flux (do 789...)
! ---------------------------------------------------------------------

  include "cvthermo.h"
  include "cv3param.h"
  include "cvflag.h"

!inputs:
  INTEGER, INTENT (IN)                               :: ncum, nd, na, ntra, nloc
  INTEGER, DIMENSION (nloc), INTENT (IN)             :: icb, inb, nk
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: sig
  REAL, DIMENSION (nloc), INTENT (IN)                :: qnk, unk, vnk
  REAL, DIMENSION (nloc, nd+1), INTENT (IN)          :: ph
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: t, rr, rs
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: u, v
  REAL, DIMENSION (nloc, nd, ntra), INTENT (IN)      :: tra               ! input of convect3
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: lv, h, hp
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: lf, frac
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: tv, tvp, ep, clw
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: m                 ! input of convect3

!outputs:
  REAL, DIMENSION (nloc, na, na), INTENT (OUT)        :: ment, qent
  REAL, DIMENSION (nloc, na, na), INTENT (OUT)        :: uent, vent
  REAL, DIMENSION (nloc, na, na), INTENT (OUT)        :: sij, elij
  REAL, DIMENSION (nloc, nd, nd, ntra), INTENT (OUT)  :: traent
  REAL, DIMENSION (nloc, nd, nd), INTENT (OUT)        :: ments, qents
  INTEGER, DIMENSION (nloc, nd), INTENT (OUT)         :: nent

!local variables:
  INTEGER i, j, k, il, im, jm
  INTEGER num1, num2
  REAL rti, bf2, anum, denom, dei, altem, cwat, stemp, qp
  REAL alt, smid, sjmin, sjmax, delp, delm
  REAL asij(nloc), smax(nloc), scrit(nloc)
  REAL asum(nloc, nd), bsum(nloc, nd), csum(nloc, nd)
  REAL sigij(nloc, nd, nd)
  REAL wgh
  REAL zm(nloc, na)
  LOGICAL lwork(nloc)

! =====================================================================
! --- INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
! =====================================================================

! ori        do 360 i=1,ncum*nlp
  DO j = 1, nl
    DO i = 1, ncum
      nent(i, j) = 0
! in convect3, m is computed in cv3_closure
! ori          m(i,1)=0.0
    END DO
  END DO

! ori      do 400 k=1,nlp
! ori       do 390 j=1,nlp
  DO j = 1, nl
    DO k = 1, nl
      DO i = 1, ncum
        qent(i, k, j) = rr(i, j)
        uent(i, k, j) = u(i, j)
        vent(i, k, j) = v(i, j)
        elij(i, k, j) = 0.0
!ym            ment(i,k,j)=0.0
!ym            sij(i,k,j)=0.0
      END DO
    END DO
  END DO

!ym
  ment(1:ncum, 1:nd, 1:nd) = 0.0
  sij(1:ncum, 1:nd, 1:nd) = 0.0

!AC!      do k=1,ntra
!AC!       do j=1,nd  ! instead nlp
!AC!        do i=1,nd ! instead nlp
!AC!         do il=1,ncum
!AC!            traent(il,i,j,k)=tra(il,j,k)
!AC!         enddo
!AC!        enddo
!AC!       enddo
!AC!      enddo
  zm(:, :) = 0.

! =====================================================================
! --- CALCULATE ENTRAINED AIR MASS FLUX (ment), TOTAL WATER MIXING
! --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
! --- FRACTION (sij)
! =====================================================================

  DO i = minorig + 1, nl

    DO j = minorig, nl
      DO il = 1, ncum
        IF ((i>=icb(il)) .AND. (i<=inb(il)) .AND. (j>=(icb(il)-1)) .AND. (j<=inb(il))) THEN

          rti = qnk(il) - ep(il, i)*clw(il, i)
          bf2 = 1. + lv(il, j)*lv(il, j)*rs(il, j)/(rrv*t(il,j)*t(il,j)*cpd)


          IF (cvflag_ice) THEN
! print*,cvflag_ice,'cvflag_ice dans do 700'
            IF (t(il,j)<=263.15) THEN
              bf2 = 1. + (lf(il,j)+lv(il,j))*(lv(il,j)+frac(il,j)* &
                   lf(il,j))*rs(il, j)/(rrv*t(il,j)*t(il,j)*cpd)
            END IF
          END IF

          anum = h(il, j) - hp(il, i) + (cpv-cpd)*t(il, j)*(rti-rr(il,j))
          denom = h(il, i) - hp(il, i) + (cpd-cpv)*(rr(il,i)-rti)*t(il, j)
          dei = denom
          IF (abs(dei)<0.01) dei = 0.01
          sij(il, i, j) = anum/dei
          sij(il, i, i) = 1.0
          altem = sij(il, i, j)*rr(il, i) + (1.-sij(il,i,j))*rti - rs(il, j)
          altem = altem/bf2
          cwat = clw(il, j)*(1.-ep(il,j))
          stemp = sij(il, i, j)
          IF ((stemp<0.0 .OR. stemp>1.0 .OR. altem>cwat) .AND. j>i) THEN

            IF (cvflag_ice) THEN
              anum = anum - (lv(il,j)+frac(il,j)*lf(il,j))*(rti-rs(il,j)-cwat*bf2)
              denom = denom + (lv(il,j)+frac(il,j)*lf(il,j))*(rr(il,i)-rti)
            ELSE
              anum = anum - lv(il, j)*(rti-rs(il,j)-cwat*bf2)
              denom = denom + lv(il, j)*(rr(il,i)-rti)
            END IF

            IF (abs(denom)<0.01) denom = 0.01
            sij(il, i, j) = anum/denom
            altem = sij(il, i, j)*rr(il, i) + (1.-sij(il,i,j))*rti - rs(il, j)
            altem = altem - (bf2-1.)*cwat
          END IF
          IF (sij(il,i,j)>0.0 .AND. sij(il,i,j)<0.95) THEN
            qent(il, i, j) = sij(il, i, j)*rr(il, i) + (1.-sij(il,i,j))*rti
            uent(il, i, j) = sij(il, i, j)*u(il, i) + (1.-sij(il,i,j))*unk(il)
            vent(il, i, j) = sij(il, i, j)*v(il, i) + (1.-sij(il,i,j))*vnk(il)
!!!!      do k=1,ntra
!!!!      traent(il,i,j,k)=sij(il,i,j)*tra(il,i,k)
!!!!     :      +(1.-sij(il,i,j))*tra(il,nk(il),k)
!!!!      end do
            elij(il, i, j) = altem
            elij(il, i, j) = max(0.0, elij(il,i,j))
            ment(il, i, j) = m(il, i)/(1.-sij(il,i,j))
            nent(il, i) = nent(il, i) + 1
          END IF
          sij(il, i, j) = max(0.0, sij(il,i,j))
          sij(il, i, j) = amin1(1.0, sij(il,i,j))
        END IF ! new
      END DO
    END DO

!AC!       do k=1,ntra
!AC!        do j=minorig,nl
!AC!         do il=1,ncum
!AC!          if( (i.ge.icb(il)).and.(i.le.inb(il)).and.
!AC!     :       (j.ge.(icb(il)-1)).and.(j.le.inb(il)))then
!AC!            traent(il,i,j,k)=sij(il,i,j)*tra(il,i,k)
!AC!     :            +(1.-sij(il,i,j))*tra(il,nk(il),k)
!AC!          endif
!AC!         enddo
!AC!        enddo
!AC!       enddo


! ***   if no air can entrain at level i assume that updraft detrains  ***
! ***   at that level and calculate detrained air flux and properties  ***


! @      do 170 i=icb(il),inb(il)

    DO il = 1, ncum
      IF ((i>=icb(il)) .AND. (i<=inb(il)) .AND. (nent(il,i)==0)) THEN
! @      if(nent(il,i).eq.0)then
        ment(il, i, i) = m(il, i)
        qent(il, i, i) = qnk(il) - ep(il, i)*clw(il, i)
        uent(il, i, i) = unk(il)
        vent(il, i, i) = vnk(il)
        elij(il, i, i) = clw(il, i)
! MAF      sij(il,i,i)=1.0
        sij(il, i, i) = 0.0
      END IF
    END DO
  END DO

!AC!      do j=1,ntra
!AC!       do i=minorig+1,nl
!AC!        do il=1,ncum
!AC!         if (i.ge.icb(il) .and. i.le.inb(il) .and. nent(il,i).eq.0) then
!AC!          traent(il,i,i,j)=tra(il,nk(il),j)
!AC!         endif
!AC!        enddo
!AC!       enddo
!AC!      enddo

  DO j = minorig, nl
    DO i = minorig, nl
      DO il = 1, ncum
        IF ((j>=(icb(il)-1)) .AND. (j<=inb(il)) .AND. (i>=icb(il)) .AND. (i<=inb(il))) THEN
          sigij(il, i, j) = sij(il, i, j)
        END IF
      END DO
    END DO
  END DO
! @      enddo

! @170   continue

! =====================================================================
! ---  NORMALIZE ENTRAINED AIR MASS FLUXES
! ---  TO REPRESENT EQUAL PROBABILITIES OF MIXING
! =====================================================================

  CALL zilch(asum, nloc*nd)
  CALL zilch(csum, nloc*nd)
  CALL zilch(csum, nloc*nd)

  DO il = 1, ncum
    lwork(il) = .FALSE.
  END DO

  DO i = minorig + 1, nl

    num1 = 0
    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il)) num1 = num1 + 1
    END DO
    IF (num1<=0) GO TO 789


    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il)) THEN
        lwork(il) = (nent(il,i)/=0)
        qp = qnk(il) - ep(il, i)*clw(il, i)

        IF (cvflag_ice) THEN

          anum = h(il, i) - hp(il, i) - (lv(il,i)+frac(il,i)*lf(il,i))* &
                       (qp-rs(il,i)) + (cpv-cpd)*t(il, i)*(qp-rr(il,i))
          denom = h(il, i) - hp(il, i) + (lv(il,i)+frac(il,i)*lf(il,i))* &
                       (rr(il,i)-qp) + (cpd-cpv)*t(il, i)*(rr(il,i)-qp)
        ELSE

          anum = h(il, i) - hp(il, i) - lv(il, i)*(qp-rs(il,i)) + &
                       (cpv-cpd)*t(il, i)*(qp-rr(il,i))
          denom = h(il, i) - hp(il, i) + lv(il, i)*(rr(il,i)-qp) + &
                       (cpd-cpv)*t(il, i)*(rr(il,i)-qp)
        END IF

        IF (abs(denom)<0.01) denom = 0.01
        scrit(il) = anum/denom
        alt = qp - rs(il, i) + scrit(il)*(rr(il,i)-qp)
        IF (scrit(il)<=0.0 .OR. alt<=0.0) scrit(il) = 1.0
        smax(il) = 0.0
        asij(il) = 0.0
      END IF
    END DO

    DO j = nl, minorig, -1

      num2 = 0
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. &
            j>=(icb(il)-1) .AND. j<=inb(il) .AND. &
            lwork(il)) num2 = num2 + 1
      END DO
      IF (num2<=0) GO TO 175

      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. &
            j>=(icb(il)-1) .AND. j<=inb(il) .AND. &
            lwork(il)) THEN

          IF (sij(il,i,j)>1.0E-16 .AND. sij(il,i,j)<0.95) THEN
            wgh = 1.0
            IF (j>i) THEN
              sjmax = max(sij(il,i,j+1), smax(il))
              sjmax = amin1(sjmax, scrit(il))
              smax(il) = max(sij(il,i,j), smax(il))
              sjmin = max(sij(il,i,j-1), smax(il))
              sjmin = amin1(sjmin, scrit(il))
              IF (sij(il,i,j)<(smax(il)-1.0E-16)) wgh = 0.0
              smid = amin1(sij(il,i,j), scrit(il))
            ELSE
              sjmax = max(sij(il,i,j+1), scrit(il))
              smid = max(sij(il,i,j), scrit(il))
              sjmin = 0.0
              IF (j>1) sjmin = sij(il, i, j-1)
              sjmin = max(sjmin, scrit(il))
            END IF
            delp = abs(sjmax-smid)
            delm = abs(sjmin-smid)
            asij(il) = asij(il) + wgh*(delp+delm)
            ment(il, i, j) = ment(il, i, j)*(delp+delm)*wgh
          END IF
        END IF
      END DO

175 END DO

    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il)) THEN
        asij(il) = max(1.0E-16, asij(il))
        asij(il) = 1.0/asij(il)
        asum(il, i) = 0.0
        bsum(il, i) = 0.0
        csum(il, i) = 0.0
      END IF
    END DO

    DO j = minorig, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. &
            j>=(icb(il)-1) .AND. j<=inb(il)) THEN
          ment(il, i, j) = ment(il, i, j)*asij(il)
        END IF
      END DO
    END DO

    DO j = minorig, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. &
            j>=(icb(il)-1) .AND. j<=inb(il)) THEN
          asum(il, i) = asum(il, i) + ment(il, i, j)
          ment(il, i, j) = ment(il, i, j)*sig(il, j)
          bsum(il, i) = bsum(il, i) + ment(il, i, j)
        END IF
      END DO
    END DO

    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il)) THEN
        bsum(il, i) = max(bsum(il,i), 1.0E-16)
        bsum(il, i) = 1.0/bsum(il, i)
      END IF
    END DO

    DO j = minorig, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. &
            j>=(icb(il)-1) .AND. j<=inb(il)) THEN
          ment(il, i, j) = ment(il, i, j)*asum(il, i)*bsum(il, i)
        END IF
      END DO
    END DO

    DO j = minorig, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. &
            j>=(icb(il)-1) .AND. j<=inb(il)) THEN
          csum(il, i) = csum(il, i) + ment(il, i, j)
        END IF
      END DO
    END DO

    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. &
          csum(il,i)<m(il,i)) THEN
        nent(il, i) = 0
        ment(il, i, i) = m(il, i)
        qent(il, i, i) = qnk(il) - ep(il, i)*clw(il, i)
        uent(il, i, i) = unk(il)
        vent(il, i, i) = vnk(il)
        elij(il, i, i) = clw(il, i)
! MAF        sij(il,i,i)=1.0
        sij(il, i, i) = 0.0
      END IF
    END DO ! il

!AC!      do j=1,ntra
!AC!       do il=1,ncum
!AC!        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
!AC!     :     .and. csum(il,i).lt.m(il,i) ) then
!AC!         traent(il,i,i,j)=tra(il,nk(il),j)
!AC!        endif
!AC!       enddo
!AC!      enddo
789 END DO

! MAF: renormalisation de MENT
  CALL zilch(zm, nloc*na)
  DO jm = 1, nl
    DO im = 1, nl
      DO il = 1, ncum
        zm(il, im) = zm(il, im) + (1.-sij(il,im,jm))*ment(il, im, jm)
      END DO
    END DO
  END DO

  DO jm = 1, nl
    DO im = 1, nl
      DO il = 1, ncum
        IF (zm(il,im)/=0.) THEN
          ment(il, im, jm) = ment(il, im, jm)*m(il, im)/zm(il, im)
        END IF
      END DO
    END DO
  END DO

  DO jm = 1, nl
    DO im = 1, nl
      DO il = 1, ncum
        qents(il, im, jm) = qent(il, im, jm)
        ments(il, im, jm) = ment(il, im, jm)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cv3_mixing

SUBROUTINE cv3_unsat(nloc, ncum, nd, na, ntra, icb, inb, iflag, &
                     t, rr, rs, gz, u, v, tra, p, ph, &
                     th, tv, lv, lf, cpn, ep, sigp, clw, frac_s, qpreca, frac_a, qta , &                       !!jygprl
                     m, ment, elij, delt, plcl, coef_clos, &
                     mp, rp, up, vp, trap, wt, water, evap, fondue, ice, &
                     faci, b, sigd, &
                     wdtrainA, wdtrainS, wdtrainM)                                      ! RomP
  USE print_control_mod, ONLY: prt_level, lunout
  IMPLICIT NONE


  include "cvthermo.h"
  include "cv3param.h"
  include "cvflag.h"
  include "nuage.h"

!inputs:
  INTEGER, INTENT (IN)                               :: ncum, nd, na, ntra, nloc
  INTEGER, DIMENSION (nloc), INTENT (IN)             :: icb, inb
  REAL, INTENT(IN)                                   :: delt
  REAL, DIMENSION (nloc), INTENT (IN)                :: plcl
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: t, rr, rs
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: gz
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: u, v
  REAL, DIMENSION (nloc, nd, ntra), INTENT(IN)       :: tra
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: p
  REAL, DIMENSION (nloc, nd+1), INTENT (IN)          :: ph
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: ep, sigp, clw   !adiab ascent shedding
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: frac_s          !ice fraction in adiab ascent shedding !!jygprl
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: qpreca          !adiab ascent precip                   !!jygprl
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: frac_a          !ice fraction in adiab ascent precip   !!jygprl
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: qta             !adiab ascent specific total water     !!jygprl
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: th, tv, lv, cpn
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: lf
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: m
  REAL, DIMENSION (nloc, na, na), INTENT (IN)        :: ment, elij
  REAL, DIMENSION (nloc), INTENT (IN)                :: coef_clos

!input/output
  INTEGER, DIMENSION (nloc), INTENT (INOUT)          :: iflag(nloc)

!outputs:
  REAL, DIMENSION (nloc, na), INTENT (OUT)           :: mp, rp, up, vp
  REAL, DIMENSION (nloc, na), INTENT (OUT)           :: water, evap, wt
  REAL, DIMENSION (nloc, na), INTENT (OUT)           :: ice, fondue
  REAL, DIMENSION (nloc, na), INTENT (OUT)           :: faci            ! ice fraction in precipitation
  REAL, DIMENSION (nloc, na, ntra), INTENT (OUT)     :: trap
  REAL, DIMENSION (nloc, na), INTENT (OUT)           :: b
  REAL, DIMENSION (nloc), INTENT (OUT)               :: sigd
! 25/08/10 - RomP---- ajout des masses precipitantes ejectees
! de l ascendance adiabatique et des flux melanges Pa et Pm.
! Distinction des wdtrain
! Pa = wdtrainA     Pm = wdtrainM
  REAL, DIMENSION (nloc, na), INTENT (OUT)           :: wdtrainA, wdtrainS, wdtrainM

!local variables
  INTEGER i, j, k, il, num1, ndp1
  REAL smallestreal
  REAL tinv, delti, coef
  REAL awat, afac, afac1, afac2, bfac
  REAL pr1, pr2, sigt, b6, c6, d6, e6, f6, revap, delth
  REAL amfac, amp2, xf, tf, fac2, ur, sru, fac, d, af, bf
  REAL ampmax, thaw
  REAL tevap(nloc)
  REAL, DIMENSION (nloc, na)      :: lvcp, lfcp
  REAL, DIMENSION (nloc, na)      :: h, hm
  REAL, DIMENSION (nloc, na)      :: ma
  REAL, DIMENSION (nloc, na)      :: frac          ! ice fraction in precipitation source
  REAL, DIMENSION (nloc, na)      :: fraci         ! provisionnal ice fraction in precipitation
  REAL, DIMENSION (nloc, na)      :: prec
  REAL wdtrain(nloc)
  LOGICAL lwork(nloc), mplus(nloc)


! ------------------------------------------------------
IF (prt_level .GE. 10) print *,' ->cv3_unsat, iflag(1) ', iflag(1)

smallestreal=tiny(smallestreal)

! =============================
! --- INITIALIZE OUTPUT ARRAYS 
! =============================
!  (loops up to nl+1)
mp(:,:) = 0.
rp(:,:) = 0.
up(:,:) = 0.
vp(:,:) = 0.
water(:,:) = 0.
evap(:,:) = 0.
wt(:,:) = 0.
ice(:,:) = 0.
fondue(:,:) = 0.
faci(:,:) = 0.
b(:,:) = 0.
sigd(:) = 0.
!! RomP >>>
wdtrainA(:,:) = 0.
wdtrainS(:,:) = 0.
wdtrainM(:,:) = 0.
!! RomP <<<

  DO i = 1, nlp
    DO il = 1, ncum
      rp(il, i) = rr(il, i)
      up(il, i) = u(il, i)
      vp(il, i) = v(il, i)
      wt(il, i) = 0.001
    END DO
  END DO

! ***  Set the fractionnal area sigd of precipitating downdraughts
  DO il = 1, ncum
    sigd(il) = sigdz*coef_clos(il)
  END DO

! =====================================================================
! --- INITIALIZE VARIOUS ARRAYS AND PARAMETERS USED IN THE COMPUTATIONS
! =====================================================================
!  (loops up to nl+1)

  delti = 1./delt
  tinv = 1./3.

  DO i = 1, nlp
    DO il = 1, ncum
      frac(il, i) = 0.0
      fraci(il, i) = 0.0
      prec(il, i) = 0.0
      lvcp(il, i) = lv(il, i)/cpn(il, i)
      lfcp(il, i) = lf(il, i)/cpn(il, i)
    END DO
  END DO

!AC!        do k=1,ntra
!AC!         do i=1,nd
!AC!          do il=1,ncum
!AC!           trap(il,i,k)=tra(il,i,k)
!AC!          enddo
!AC!         enddo
!AC!        enddo

! ***  check whether ep(inb)=0, if so, skip precipitating    ***
! ***             downdraft calculation                      ***


  DO il = 1, ncum
!!          lwork(il)=.TRUE.
!!          if(ep(il,inb(il)).lt.0.0001)lwork(il)=.FALSE.
!jyg<
!!    lwork(il) = ep(il, inb(il)) >= 0.0001
    lwork(il) = ep(il, inb(il)) >= 0.0001 .AND. iflag(il) <= 2
  END DO

!
! Get adiabatic ascent mass flux 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (adiab_ascent_mass_flux_depends_on_ejectliq) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Warning : this option leads to water conservation violation
!!!           Expert only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO il = 1, ncum
      ma(il, nlp) = 0.
      ma(il, 1)   = 0.
    END DO

  DO i = nl, 2, -1
      DO il = 1, ncum
        ma(il, i) = ma(il, i+1)*(1.-qta(il,i))/(1.-qta(il,i-1)) + m(il, i)
      END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ELSE ! (adiab_ascent_mass_flux_depends_on_ejectliq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO il = 1, ncum
      ma(il, nlp) = 0.
      ma(il, 1)   = 0.
    END DO

  DO i = nl, 2, -1
      DO il = 1, ncum
        ma(il, i) = ma(il, i+1) + m(il, i)
      END DO
  END DO

  ENDIF ! (adiab_ascent_mass_flux_depends_on_ejectliq) ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! ***                    begin downdraft loop                    ***
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  DO i = nl + 1, 1, -1

    num1 = 0
    DO il = 1, ncum
      IF (i<=inb(il) .AND. lwork(il)) num1 = num1 + 1
    END DO
    IF (num1<=0) GO TO 400

    CALL zilch(wdtrain, ncum)


! ***  integrate liquid water equation to find condensed water   ***
! ***                and condensed water flux                    ***
!
!
! ***              calculate detrained precipitation             ***


    DO il = 1, ncum                                                   
      IF (i<=inb(il) .AND. lwork(il)) THEN                            
        wdtrain(il) = grav*ep(il, i)*m(il, i)*clw(il, i)            
        wdtrainS(il, i) = wdtrain(il)/grav                                            !   Ps   jyg
!!        wdtrainA(il, i) = wdtrain(il)/grav                                          !   Ps   RomP
      END IF                                                          
    END DO                                                            

    IF (i>1) THEN
      DO j = 1, i - 1
        DO il = 1, ncum
          IF (i<=inb(il) .AND. lwork(il)) THEN
            awat = elij(il, j, i) - (1.-ep(il,i))*clw(il, i)
            awat = max(awat, 0.0)
            wdtrain(il) = wdtrain(il) + grav*awat*ment(il, j, i)
            wdtrainM(il, i) = wdtrain(il)/grav - wdtrainS(il, i)    !   Pm  jyg
!!            wdtrainM(il, i) = wdtrain(il)/grav - wdtrainA(il, i)  !   Pm  RomP
          END IF
        END DO
      END DO
    END IF

    IF (cvflag_prec_eject) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (adiab_ascent_mass_flux_depends_on_ejectliq) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Warning : this option leads to water conservation violation
!!!           Expert only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF ( i > 1) THEN
            DO il = 1, ncum
              IF (i<=inb(il) .AND. lwork(il)) THEN
                wdtrainA(il,i) = ma(il, i+1)*(qta(il, i-1)-qta(il,i))/(1. - qta(il, i-1))    !   Pa   jygprl
                wdtrain(il) = wdtrain(il) + grav*wdtrainA(il,i)
              END IF
            END DO
          ENDIF  ! ( i > 1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSE ! (adiab_ascent_mass_flux_depends_on_ejectliq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF ( i > 1) THEN
            DO il = 1, ncum
              IF (i<=inb(il) .AND. lwork(il)) THEN
                wdtrainA(il,i) = ma(il, i+1)*(qta(il, i-1)-qta(il,i))                        !   Pa   jygprl
                wdtrain(il) = wdtrain(il) + grav*wdtrainA(il,i)
              END IF
            END DO
          ENDIF  ! ( i > 1)

      ENDIF ! (adiab_ascent_mass_flux_depends_on_ejectliq) ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ENDIF  ! (cvflag_prec_eject)


! ***    find rain water and evaporation using provisional   ***
! ***              estimates of rp(i)and rp(i-1)             ***


    IF (cvflag_ice) THEN                                                                                !!jygprl
      IF (cvflag_prec_eject) THEN
        DO il = 1, ncum                                                                                   !!jygprl
          IF (i<=inb(il) .AND. lwork(il)) THEN                                                            !!jygprl
            frac(il, i) = (frac_a(il,i)*wdtrainA(il,i)+frac_s(il,i)*(wdtrainS(il,i)+wdtrainM(il,i))) / &  !!jygprl
                          max(wdtrainA(il,i)+wdtrainS(il,i)+wdtrainM(il,i),smallestreal)                  !!jygprl
            fraci(il, i) = frac(il, i)                                                                    !!jygprl
          END IF                                                                                          !!jygprl
        END DO                                                                                            !!jygprl
      ELSE  ! (cvflag_prec_eject)
        DO il = 1, ncum                                                                                   !!jygprl
          IF (i<=inb(il) .AND. lwork(il)) THEN                                                            !!jygprl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF (keepbug_ice_frac) THEN
              frac(il, i) = frac_s(il, i)
!       Ice fraction computed again here as a function of the temperature seen by unsaturated downdraughts
!       (i.e. the cold pool temperature) for compatibility with earlier versions.
              fraci(il, i) = 1. - (t(il,i)-243.15)/(263.15-243.15)
              fraci(il, i) = min(max(fraci(il,i),0.0), 1.0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ELSE  ! (keepbug_ice_frac)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              frac(il, i) = frac_s(il, i)
              fraci(il, i) = frac(il, i)                                                                    !!jygprl
            ENDIF  ! (keepbug_ice_frac)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          END IF                                                                                          !!jygprl
        END DO                                                                                            !!jygprl
      ENDIF  ! (cvflag_prec_eject)
    END IF                                                                                              !!jygprl


    DO il = 1, ncum
      IF (i<=inb(il) .AND. lwork(il)) THEN

        wt(il, i) = 45.0

        IF (i<inb(il)) THEN
          rp(il, i) = rp(il, i+1) + &
                      (cpd*(t(il,i+1)-t(il,i))+gz(il,i+1)-gz(il,i))/lv(il, i)
          rp(il, i) = 0.5*(rp(il,i)+rr(il,i))
        END IF
        rp(il, i) = max(rp(il,i), 0.0)
        rp(il, i) = amin1(rp(il,i), rs(il,i))
        rp(il, inb(il)) = rr(il, inb(il))

        IF (i==1) THEN
          afac = p(il, 1)*(rs(il,1)-rp(il,1))/(1.0E4+2000.0*p(il,1)*rs(il,1))
          IF (cvflag_ice) THEN
            afac1 = p(il, i)*(rs(il,1)-rp(il,1))/(1.0E4+2000.0*p(il,1)*rs(il,1))
          END IF
        ELSE
          rp(il, i-1) = rp(il, i) + (cpd*(t(il,i)-t(il,i-1))+gz(il,i)-gz(il,i-1))/lv(il, i)
          rp(il, i-1) = 0.5*(rp(il,i-1)+rr(il,i-1))
          rp(il, i-1) = amin1(rp(il,i-1), rs(il,i-1))
          rp(il, i-1) = max(rp(il,i-1), 0.0)
          afac1 = p(il, i)*(rs(il,i)-rp(il,i))/(1.0E4+2000.0*p(il,i)*rs(il,i))
          afac2 = p(il, i-1)*(rs(il,i-1)-rp(il,i-1))/(1.0E4+2000.0*p(il,i-1)*rs(il,i-1))
          afac = 0.5*(afac1+afac2)
        END IF
        IF (i==inb(il)) afac = 0.0
        afac = max(afac, 0.0)
        bfac = 1./(sigd(il)*wt(il,i))

!
    IF (prt_level >= 20) THEN
      Print*, 'cv3_unsat after provisional rp estimate: rp, afac, bfac ', &
          i, rp(1, i), afac,bfac
    ENDIF
!
!JYG1
! cc        sigt=1.0
! cc        if(i.ge.icb)sigt=sigp(i)
! prise en compte de la variation progressive de sigt dans
! les couches icb et icb-1:
! pour plcl<ph(i+1), pr1=0 & pr2=1
! pour plcl>ph(i),   pr1=1 & pr2=0
! pour ph(i+1)<plcl<ph(i), pr1 est la proportion a cheval
! sur le nuage, et pr2 est la proportion sous la base du
! nuage.
        pr1 = (plcl(il)-ph(il,i+1))/(ph(il,i)-ph(il,i+1))
        pr1 = max(0., min(1.,pr1))
        pr2 = (ph(il,i)-plcl(il))/(ph(il,i)-ph(il,i+1))
        pr2 = max(0., min(1.,pr2))
        sigt = sigp(il, i)*pr1 + pr2
!JYG2

!JYG----
!    b6 = bfac*100.*sigd(il)*(ph(il,i)-ph(il,i+1))*sigt*afac
!    c6 = water(il,i+1) + wdtrain(il)*bfac
!    c6 = prec(il,i+1) + wdtrain(il)*bfac
!    revap=0.5*(-b6+sqrt(b6*b6+4.*c6))
!    evap(il,i)=sigt*afac*revap
!    water(il,i)=revap*revap
!    prec(il,i)=revap*revap
!!        print *,' i,b6,c6,revap,evap(il,i),water(il,i),wdtrain(il) ', &
!!                 i,b6,c6,revap,evap(il,i),water(il,i),wdtrain(il)
!!---end jyg---

! --------retour à la formulation originale d'Emanuel.
        IF (cvflag_ice) THEN

!   b6=bfac*50.*sigd(il)*(ph(il,i)-ph(il,i+1))*sigt*afac
!   c6=prec(il,i+1)+bfac*wdtrain(il) &
!       -50.*sigd(il)*bfac*(ph(il,i)-ph(il,i+1))*evap(il,i+1)
!   if(c6.gt.0.0)then
!   revap=0.5*(-b6+sqrt(b6*b6+4.*c6))

!JAM  Attention: evap=sigt*E
!    Modification: evap devient l'évaporation en milieu de couche
!    car nécessaire dans cv3_yield
!    Du coup, il faut modifier pas mal d'équations...
!    et l'expression de afac qui devient afac1
!    revap=sqrt((prec(i+1)+prec(i))/2)

          b6 = bfac*50.*sigd(il)*(ph(il,i)-ph(il,i+1))*sigt*afac1
          c6 = prec(il, i+1) + 0.5*bfac*wdtrain(il)
! print *,'bfac,sigd(il),sigt,afac1 ',bfac,sigd(il),sigt,afac1
! print *,'prec(il,i+1),wdtrain(il) ',prec(il,i+1),wdtrain(il)
! print *,'b6,c6,b6*b6+4.*c6 ',b6,c6,b6*b6+4.*c6
          IF (c6>b6*b6+1.E-20) THEN
            revap = 2.*c6/(b6+sqrt(b6*b6+4.*c6))
          ELSE
            revap = (-b6+sqrt(b6*b6+4.*c6))/2.
          END IF
          prec(il, i) = max(0., 2.*revap*revap-prec(il,i+1))
! print*,prec(il,i),'neige'

!JYG    Dans sa formulation originale, Emanuel calcule l'evaporation par:
! c             evap(il,i)=sigt*afac*revap
! ce qui n'est pas correct. Dans cv_routines, la formulation a été modifiee.
! Ici,l'evaporation evap est simplement calculee par l'equation de
! conservation.
! prec(il,i)=revap*revap
! else
!JYG----   Correction : si c6 <= 0, water(il,i)=0.
! prec(il,i)=0.
! endif

!JYG---   Dans tous les cas, evaporation = [tt ce qui entre dans la couche i]
! moins [tt ce qui sort de la couche i]
! print *, 'evap avec ice'
          evap(il, i) = (wdtrain(il)+sigd(il)*wt(il,i)*(prec(il,i+1)-prec(il,i))) / &
                        (sigd(il)*(ph(il,i)-ph(il,i+1))*100.)
!
    IF (prt_level >= 20) THEN
      Print*, 'cv3_unsat after evap computation: wdtrain, sigd, wt, prec(i+1),prec(i) ', &
          i, wdtrain(1), sigd(1), wt(1,i), prec(1,i+1),prec(1,i)
    ENDIF
!

!jyg<
          d6 = prec(il,i)-prec(il,i+1)

!!          d6 = bfac*wdtrain(il) - 100.*sigd(il)*bfac*(ph(il,i)-ph(il,i+1))*evap(il, i)
!!          e6 = bfac*wdtrain(il)
!!          f6 = -100.*sigd(il)*bfac*(ph(il,i)-ph(il,i+1))*evap(il, i)
!>jyg
!CR:tmax_fonte_cv: T for which ice is totally melted (used to be 275.15)
          thaw = (t(il,i)-273.15)/(tmax_fonte_cv-273.15)
          thaw = min(max(thaw,0.0), 1.0)
!jyg<
          water(il, i) = water(il, i+1) + (1-fraci(il,i))*d6
          ice(il, i)   = ice(il, i+1)   + fraci(il, i)*d6
          water(il, i) = min(prec(il,i), max(water(il,i), 0.))
          ice(il, i)   = min(prec(il,i), max(ice(il,i),   0.))

!!          water(il, i) = water(il, i+1) + (1-fraci(il,i))*d6
!!          water(il, i) = max(water(il,i), 0.)
!!          ice(il, i) = ice(il, i+1) + fraci(il, i)*d6
!!          ice(il, i) = max(ice(il,i), 0.)
!>jyg
          fondue(il, i) = ice(il, i)*thaw
          water(il, i) = water(il, i) + fondue(il, i)
          ice(il, i) = ice(il, i) - fondue(il, i)

          IF (water(il,i)+ice(il,i)<1.E-30) THEN
            faci(il, i) = 0.
          ELSE
            faci(il, i) = ice(il, i)/(water(il,i)+ice(il,i))
          END IF

!           water(il,i)=water(il,i+1)+(1.-fraci(il,i))*e6+(1.-faci(il,i))*f6
!           water(il,i)=max(water(il,i),0.)
!           ice(il,i)=ice(il,i+1)+fraci(il,i)*e6+faci(il,i)*f6
!           ice(il,i)=max(ice(il,i),0.)
!           fondue(il,i)=ice(il,i)*thaw
!           water(il,i)=water(il,i)+fondue(il,i)
!           ice(il,i)=ice(il,i)-fondue(il,i)
            
!           if((water(il,i)+ice(il,i)).lt.1.e-30)then
!             faci(il,i)=0.
!           else
!             faci(il,i)=ice(il,i)/(water(il,i)+ice(il,i))
!           endif

        ELSE
          b6 = bfac*50.*sigd(il)*(ph(il,i)-ph(il,i+1))*sigt*afac
          c6 = water(il, i+1) + bfac*wdtrain(il) - &
               50.*sigd(il)*bfac*(ph(il,i)-ph(il,i+1))*evap(il, i+1)
          IF (c6>0.0) THEN
            revap = 0.5*(-b6+sqrt(b6*b6+4.*c6))
            water(il, i) = revap*revap
          ELSE
            water(il, i) = 0.
          END IF
! print *, 'evap sans ice'
          evap(il, i) = (wdtrain(il)+sigd(il)*wt(il,i)*(water(il,i+1)-water(il,i)))/ &
                        (sigd(il)*(ph(il,i)-ph(il,i+1))*100.)

        END IF
      END IF !(i.le.inb(il) .and. lwork(il))
    END DO
! ----------------------------------------------------------------

! cc
! ***  calculate precipitating downdraft mass flux under     ***
! ***              hydrostatic approximation                 ***

    DO il = 1, ncum
      IF (i<=inb(il) .AND. lwork(il) .AND. i/=1) THEN

        tevap(il) = max(0.0, evap(il,i))
        delth = max(0.001, (th(il,i)-th(il,i-1)))
        IF (cvflag_ice) THEN
          IF (cvflag_grav) THEN
            mp(il, i) = 100.*ginv*(lvcp(il,i)*sigd(il)*tevap(il)* &
                                               (p(il,i-1)-p(il,i))/delth + &
                                   lfcp(il,i)*sigd(il)*faci(il,i)*tevap(il)* &
                                               (p(il,i-1)-p(il,i))/delth + &
                                   lfcp(il,i)*sigd(il)*wt(il,i)/100.*fondue(il,i)* &
                                               (p(il,i-1)-p(il,i))/delth/(ph(il,i)-ph(il,i+1)))
          ELSE
            mp(il, i) = 10.*(lvcp(il,i)*sigd(il)*tevap(il)* &
                                                (p(il,i-1)-p(il,i))/delth + &
                             lfcp(il,i)*sigd(il)*faci(il,i)*tevap(il)* &
                                                (p(il,i-1)-p(il,i))/delth + &
                             lfcp(il,i)*sigd(il)*wt(il,i)/100.*fondue(il,i)* &
                                                (p(il,i-1)-p(il,i))/delth/(ph(il,i)-ph(il,i+1)))

          END IF
        ELSE
          IF (cvflag_grav) THEN
            mp(il, i) = 100.*ginv*lvcp(il, i)*sigd(il)*tevap(il)* &
                                                (p(il,i-1)-p(il,i))/delth
          ELSE
            mp(il, i) = 10.*lvcp(il, i)*sigd(il)*tevap(il)* &
                                                (p(il,i-1)-p(il,i))/delth
          END IF

        END IF

      END IF !(i.le.inb(il) .and. lwork(il) .and. i.ne.1)
      IF (prt_level .GE. 20) THEN
        PRINT *,'cv3_unsat, mp hydrostatic ', i, mp(il,i)
      ENDIF
    END DO
! ----------------------------------------------------------------

! ***           if hydrostatic assumption fails,             ***
! ***   solve cubic difference equation for downdraft theta  ***
! ***  and mass flux from two simultaneous differential eqns ***

    DO il = 1, ncum
      IF (i<=inb(il) .AND. lwork(il) .AND. i/=1) THEN

        amfac = sigd(il)*sigd(il)*70.0*ph(il, i)*(p(il,i-1)-p(il,i))* &
                         (th(il,i)-th(il,i-1))/(tv(il,i)*th(il,i))
        amp2 = abs(mp(il,i+1)*mp(il,i+1)-mp(il,i)*mp(il,i))

        IF (amp2>(0.1*amfac)) THEN
          xf = 100.0*sigd(il)*sigd(il)*sigd(il)*(ph(il,i)-ph(il,i+1))
          tf = b(il, i) - 5.0*(th(il,i)-th(il,i-1))*t(il, i) / &
                              (lvcp(il,i)*sigd(il)*th(il,i))
          af = xf*tf + mp(il, i+1)*mp(il, i+1)*tinv

          IF (cvflag_ice) THEN
            bf = 2.*(tinv*mp(il,i+1))**3 + tinv*mp(il, i+1)*xf*tf + &
                 50.*(p(il,i-1)-p(il,i))*xf*(tevap(il)*(1.+(lf(il,i)/lv(il,i))*faci(il,i)) + &
                (lf(il,i)/lv(il,i))*wt(il,i)/100.*fondue(il,i)/(ph(il,i)-ph(il,i+1)))
          ELSE

            bf = 2.*(tinv*mp(il,i+1))**3 + tinv*mp(il, i+1)*xf*tf + &
                                           50.*(p(il,i-1)-p(il,i))*xf*tevap(il)
          END IF

          fac2 = 1.0
          IF (bf<0.0) fac2 = -1.0
          bf = abs(bf)
          ur = 0.25*bf*bf - af*af*af*tinv*tinv*tinv
          IF (ur>=0.0) THEN
            sru = sqrt(ur)
            fac = 1.0
            IF ((0.5*bf-sru)<0.0) fac = -1.0
            mp(il, i) = mp(il, i+1)*tinv + (0.5*bf+sru)**tinv + &
                                           fac*(abs(0.5*bf-sru))**tinv
          ELSE
            d = atan(2.*sqrt(-ur)/(bf+1.0E-28))
            IF (fac2<0.0) d = 3.14159 - d
            mp(il, i) = mp(il, i+1)*tinv + 2.*sqrt(af*tinv)*cos(d*tinv)
          END IF
          mp(il, i) = max(0.0, mp(il,i))
          IF (prt_level .GE. 20) THEN
            PRINT *,'cv3_unsat, mp cubic ', i, mp(il,i)
          ENDIF

          IF (cvflag_ice) THEN
            IF (cvflag_grav) THEN
!JYG : il y a vraisemblablement une erreur dans la ligne 2 suivante:
! il faut diviser par (mp(il,i)*sigd(il)*grav) et non par (mp(il,i)+sigd(il)*0.1).
! Et il faut bien revoir les facteurs 100.
              b(il, i-1) = b(il, i) + 100.0*(p(il,i-1)-p(il,i))* &
                           (tevap(il)*(1.+(lf(il,i)/lv(il,i))*faci(il,i)) + &
                           (lf(il,i)/lv(il,i))*wt(il,i)/100.*fondue(il,i) / &
                           (ph(il,i)-ph(il,i+1))) / &
                           (mp(il,i)+sigd(il)*0.1) - &
                           10.0*(th(il,i)-th(il,i-1))*t(il, i) / &
                           (lvcp(il,i)*sigd(il)*th(il,i))
            ELSE
              b(il, i-1) = b(il, i) + 100.0*(p(il,i-1)-p(il,i))*&
                           (tevap(il)*(1.+(lf(il,i)/lv(il,i))*faci(il,i)) + &
                           (lf(il,i)/lv(il,i))*wt(il,i)/100.*fondue(il,i) / &
                           (ph(il,i)-ph(il,i+1))) / &
                           (mp(il,i)+sigd(il)*0.1) - &
                           10.0*(th(il,i)-th(il,i-1))*t(il, i) / &
                           (lvcp(il,i)*sigd(il)*th(il,i))
            END IF
          ELSE
            IF (cvflag_grav) THEN
              b(il, i-1) = b(il, i) + 100.0*(p(il,i-1)-p(il,i))*tevap(il) / &
                           (mp(il,i)+sigd(il)*0.1) - &
                           10.0*(th(il,i)-th(il,i-1))*t(il, i) / &
                           (lvcp(il,i)*sigd(il)*th(il,i))
            ELSE
              b(il, i-1) = b(il, i) + 100.0*(p(il,i-1)-p(il,i))*tevap(il) / &
                           (mp(il,i)+sigd(il)*0.1) - &
                           10.0*(th(il,i)-th(il,i-1))*t(il, i) / &
                           (lvcp(il,i)*sigd(il)*th(il,i))
            END IF
          END IF
          b(il, i-1) = max(b(il,i-1), 0.0)

        END IF !(amp2.gt.(0.1*amfac))

!jyg<    This part shifted 10 lines farther
!!! ***         limit magnitude of mp(i) to meet cfl condition      ***
!!
!!        ampmax = 2.0*(ph(il,i)-ph(il,i+1))*delti
!!        amp2 = 2.0*(ph(il,i-1)-ph(il,i))*delti
!!        ampmax = min(ampmax, amp2)
!!        mp(il, i) = min(mp(il,i), ampmax)
!>jyg

! ***      force mp to decrease linearly to zero                 ***
! ***       between cloud base and the surface                   ***


! c      if(p(il,i).gt.p(il,icb(il)))then
! c       mp(il,i)=mp(il,icb(il))*(p(il,1)-p(il,i))/(p(il,1)-p(il,icb(il)))
! c      endif
        IF (ph(il,i)>0.9*plcl(il)) THEN
          mp(il, i) = mp(il, i)*(ph(il,1)-ph(il,i))/(ph(il,1)-0.9*plcl(il))
        END IF

!jyg<    Shifted part
! ***         limit magnitude of mp(i) to meet cfl condition      ***

        ampmax = 2.0*(ph(il,i)-ph(il,i+1))*delti
        amp2 = 2.0*(ph(il,i-1)-ph(il,i))*delti
        ampmax = min(ampmax, amp2)
        mp(il, i) = min(mp(il,i), ampmax)
!>jyg

      END IF ! (i.le.inb(il) .and. lwork(il) .and. i.ne.1)
    END DO
! ----------------------------------------------------------------
!
    IF (prt_level >= 20) THEN
      Print*, 'cv3_unsat after mp computation: mp, b(i), b(i-1) ', &
          i, mp(1, i), b(1,i), b(1,max(i-1,1))
    ENDIF
!

! ***       find mixing ratio of precipitating downdraft     ***

    DO il = 1, ncum
      IF (i<inb(il) .AND. lwork(il)) THEN
        mplus(il) = mp(il, i) > mp(il, i+1)
      END IF ! (i.lt.inb(il) .and. lwork(il))
    END DO

    DO il = 1, ncum
      IF (i<inb(il) .AND. lwork(il)) THEN

        rp(il, i) = rr(il, i)

        IF (mplus(il)) THEN

          IF (cvflag_grav) THEN
            rp(il, i) = rp(il, i+1)*mp(il, i+1) + rr(il, i)*(mp(il,i)-mp(il,i+1)) + &
              100.*ginv*0.5*sigd(il)*(ph(il,i)-ph(il,i+1))*(evap(il,i+1)+evap(il,i))
          ELSE
            rp(il, i) = rp(il, i+1)*mp(il, i+1) + rr(il, i)*(mp(il,i)-mp(il,i+1)) + &
              5.*sigd(il)*(ph(il,i)-ph(il,i+1))*(evap(il,i+1)+evap(il,i))
          END IF
          rp(il, i) = rp(il, i)/mp(il, i)
          up(il, i) = up(il, i+1)*mp(il, i+1) + u(il, i)*(mp(il,i)-mp(il,i+1))
          up(il, i) = up(il, i)/mp(il, i)
          vp(il, i) = vp(il, i+1)*mp(il, i+1) + v(il, i)*(mp(il,i)-mp(il,i+1))
          vp(il, i) = vp(il, i)/mp(il, i)

        ELSE ! if (mplus(il))

          IF (mp(il,i+1)>1.0E-16) THEN
            IF (cvflag_grav) THEN
              rp(il, i) = rp(il,i+1) + 100.*ginv*0.5*sigd(il)*(ph(il,i)-ph(il,i+1)) * &
                                       (evap(il,i+1)+evap(il,i))/mp(il,i+1)
            ELSE
              rp(il, i) = rp(il,i+1) + 5.*sigd(il)*(ph(il,i)-ph(il,i+1)) * &
                                       (evap(il,i+1)+evap(il,i))/mp(il, i+1)
            END IF
            up(il, i) = up(il, i+1)
            vp(il, i) = vp(il, i+1)
          END IF ! (mp(il,i+1).gt.1.0e-16)
        END IF ! (mplus(il)) else if (.not.mplus(il))

        rp(il, i) = amin1(rp(il,i), rs(il,i))
        rp(il, i) = max(rp(il,i), 0.0)

      END IF ! (i.lt.inb(il) .and. lwork(il))
    END DO
! ----------------------------------------------------------------

! ***       find tracer concentrations in precipitating downdraft     ***

!AC!      do j=1,ntra
!AC!       do il = 1,ncum
!AC!       if (i.lt.inb(il) .and. lwork(il)) then
!AC!c
!AC!         if(mplus(il))then
!AC!          trap(il,i,j)=trap(il,i+1,j)*mp(il,i+1)
!AC!     :              +trap(il,i,j)*(mp(il,i)-mp(il,i+1))
!AC!          trap(il,i,j)=trap(il,i,j)/mp(il,i)
!AC!         else ! if (mplus(il))
!AC!          if(mp(il,i+1).gt.1.0e-16)then
!AC!           trap(il,i,j)=trap(il,i+1,j)
!AC!          endif
!AC!         endif ! (mplus(il)) else if (.not.mplus(il))
!AC!c
!AC!        endif ! (i.lt.inb(il) .and. lwork(il))
!AC!       enddo
!AC!      end do

400 END DO
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ***                    end of downdraft loop                    ***

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  RETURN

END SUBROUTINE cv3_unsat

SUBROUTINE cv3_yield(nloc, ncum, nd, na, ntra, ok_conserv_q, &
                     icb, inb, delt, &
                     t, rr, t_wake, rr_wake, s_wake, u, v, tra, &
                     gz, p, ph, h, hp, lv, lf, cpn, th, th_wake, &
                     ep, clw, qpreca, m, tp, mp, rp, up, vp, trap, &
                     wt, water, ice, evap, fondue, faci, b, sigd, &
                     ment, qent, hent, iflag_mix, uent, vent, &
                     nent, elij, traent, sig, &
                     tv, tvp, wghti, &
                     iflag, precip, Vprecip, Vprecipi, &     ! jyg: Vprecipi
                     ft, fr, fu, fv, ftra, &                 ! jyg
                     cbmf, upwd, dnwd, dnwd0, ma, mip, &
!!                     tls, tps,                             ! useless . jyg
                     qcondc, wd, &
                     ftd, fqd, qta, qtc, sigt, tau_cld_cv, coefw_cld_cv)

    USE print_control_mod, ONLY: lunout, prt_level
    USE add_phys_tend_mod, only : fl_cor_ebil

  IMPLICIT NONE

  include "cvthermo.h"
  include "cv3param.h"
  include "cvflag.h"
  include "conema3.h"

!inputs:
      INTEGER, INTENT (IN)                               :: iflag_mix
      INTEGER, INTENT (IN)                               :: ncum, nd, na, ntra, nloc
      LOGICAL, INTENT (IN)                               :: ok_conserv_q
      INTEGER, DIMENSION (nloc), INTENT (IN)             :: icb, inb
      REAL, INTENT (IN)                                  :: delt
      REAL, DIMENSION (nloc, nd), INTENT (IN)            :: t, rr, u, v
      REAL, DIMENSION (nloc, nd), INTENT (IN)            :: t_wake, rr_wake
      REAL, DIMENSION (nloc), INTENT (IN)                :: s_wake
      REAL, DIMENSION (nloc, nd, ntra), INTENT (IN)      :: tra
      REAL, DIMENSION (nloc, nd), INTENT (IN)            :: p
      REAL, DIMENSION (nloc, nd+1), INTENT (IN)          :: ph
      REAL, DIMENSION (nloc, na), INTENT (IN)            :: gz, h, hp
      REAL, DIMENSION (nloc, na), INTENT (IN)            :: th, tp
      REAL, DIMENSION (nloc, na), INTENT (IN)            :: lv, cpn, ep, clw
      REAL, DIMENSION (nloc, na), INTENT (IN)            :: lf
      REAL, DIMENSION (nloc, na), INTENT (IN)            :: rp, up
      REAL, DIMENSION (nloc, na), INTENT (IN)            :: vp
      REAL, DIMENSION (nloc, nd), INTENT (IN)            :: wt
      REAL, DIMENSION (nloc, nd, ntra), INTENT (IN)      :: trap
      REAL, DIMENSION (nloc, na), INTENT (IN)            :: water, evap, b
      REAL, DIMENSION (nloc, na), INTENT (IN)            :: fondue, faci, ice
      REAL, DIMENSION (nloc, na, na), INTENT (IN)        :: qent, uent
      REAL, DIMENSION (nloc, na, na), INTENT (IN)        :: hent
      REAL, DIMENSION (nloc, na, na), INTENT (IN)        :: vent, elij
      INTEGER, DIMENSION (nloc, nd), INTENT (IN)         :: nent
      REAL, DIMENSION (nloc, na, na, ntra), INTENT (IN)  :: traent
      REAL, DIMENSION (nloc, nd), INTENT (IN)            :: tv, tvp, wghti
      REAL, DIMENSION (nloc, nd), INTENT (IN)            :: qta
      REAL, DIMENSION (nloc, na),INTENT(IN)              :: qpreca
      REAL, INTENT(IN)                                   :: tau_cld_cv, coefw_cld_cv
!
!input/output:
      REAL, DIMENSION (nloc, na), INTENT (INOUT)         :: m, mp
      REAL, DIMENSION (nloc, na, na), INTENT (INOUT)     :: ment
      INTEGER, DIMENSION (nloc), INTENT (INOUT)          :: iflag
      REAL, DIMENSION (nloc, nd), INTENT (INOUT)         :: sig
      REAL, DIMENSION (nloc), INTENT (INOUT)             :: sigd
!
!outputs:
      REAL, DIMENSION (nloc), INTENT (OUT)               :: precip
      REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: ft, fr, fu, fv
      REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: ftd, fqd
      REAL, DIMENSION (nloc, nd, ntra), INTENT (OUT)     :: ftra
      REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: upwd, dnwd, ma
      REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: dnwd0, mip
      REAL, DIMENSION (nloc, nd+1), INTENT (OUT)         :: Vprecip
      REAL, DIMENSION (nloc, nd+1), INTENT (OUT)         :: Vprecipi
!!      REAL tls(nloc, nd), tps(nloc, nd)                    ! useless . jyg
      REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: qcondc                      ! cld
      REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: qtc, sigt                   ! cld
      REAL, DIMENSION (nloc), INTENT (OUT)               :: wd                          ! gust
      REAL, DIMENSION (nloc), INTENT (OUT)               :: cbmf
!
!local variables:
      INTEGER                                            :: i, k, il, n, j, num1
      REAL                                               :: rat, delti
      REAL                                               :: ax, bx, cx, dx, ex
      REAL                                               :: cpinv, rdcp, dpinv
      REAL                                               :: sigaq
      REAL, DIMENSION (nloc)                             ::  awat
      REAL, DIMENSION (nloc, nd)                         :: lvcp, lfcp              ! , mke ! unused . jyg
      REAL, DIMENSION (nloc)                             :: am, work, ad, amp1
!!      real up1(nloc), dn1(nloc)
      REAL, DIMENSION (nloc, nd, nd)                     :: up1, dn1
!jyg<
      REAL, DIMENSION (nloc, nd)                         :: up_to, up_from
      REAL, DIMENSION (nloc, nd)                         :: dn_to, dn_from
!>jyg
      REAL, DIMENSION (nloc)                             :: asum, bsum, csum, dsum
      REAL, DIMENSION (nloc)                             :: esum, fsum, gsum, hsum
      REAL, DIMENSION (nloc, nd)                         :: th_wake
      REAL, DIMENSION (nloc)                             :: alpha_qpos, alpha_qpos1
      REAL, DIMENSION (nloc, nd)                         :: qcond, nqcond, wa           ! cld
      REAL, DIMENSION (nloc, nd)                         :: siga, sax, mac              ! cld
      REAL, DIMENSION (nloc)                             :: sument
      REAL, DIMENSION (nloc, nd)                         :: sigment, qtment             ! cld
      REAL sumdq !jyg
! 
! -------------------------------------------------------------

! initialization:

  delti = 1.0/delt
! print*,'cv3_yield initialisation delt', delt

  DO il = 1, ncum
    precip(il) = 0.0
    wd(il) = 0.0 ! gust
  END DO

!   Fluxes are on a staggered grid : loops extend up to nl+1
  DO i = 1, nlp
    DO il = 1, ncum
      Vprecip(il, i) = 0.0
      Vprecipi(il, i) = 0.0                               ! jyg
      upwd(il, i) = 0.0
      dnwd(il, i) = 0.0
      dnwd0(il, i) = 0.0
      mip(il, i) = 0.0
    END DO
  END DO
  DO i = 1, nl
    DO il = 1, ncum
      ft(il, i) = 0.0
      fr(il, i) = 0.0
      fu(il, i) = 0.0
      fv(il, i) = 0.0
      ftd(il, i) = 0.0
      fqd(il, i) = 0.0
      qcondc(il, i) = 0.0 ! cld
      qcond(il, i) = 0.0 ! cld
      qtc(il, i) = 0.0 ! cld
      qtment(il, i) = 0.0 ! cld
      sigment(il, i) = 0.0 ! cld
      sigt(il, i) = 0.0 ! cld
      nqcond(il, i) = 0.0 ! cld
    END DO
  END DO
! print*,'cv3_yield initialisation 2'
!AC!      do j=1,ntra
!AC!       do i=1,nd
!AC!        do il=1,ncum
!AC!          ftra(il,i,j)=0.0
!AC!        enddo
!AC!       enddo
!AC!      enddo
! print*,'cv3_yield initialisation 3'
  DO i = 1, nl
    DO il = 1, ncum
      lvcp(il, i) = lv(il, i)/cpn(il, i)
      lfcp(il, i) = lf(il, i)/cpn(il, i)
    END DO
  END DO



! ***  calculate surface precipitation in mm/day     ***

  DO il = 1, ncum
    IF (ep(il,inb(il))>=0.0001 .AND. iflag(il)<=1) THEN
      IF (cvflag_ice) THEN
        precip(il) = wt(il, 1)*sigd(il)*(water(il,1)+ice(il,1)) &
                              *86400.*1000./(rowl*grav)
      ELSE
        precip(il) = wt(il, 1)*sigd(il)*water(il, 1) &
                              *86400.*1000./(rowl*grav)
      END IF
    END IF
  END DO
! print*,'cv3_yield apres calcul precip'


! ===  calculate vertical profile of  precipitation in kg/m2/s  ===

  DO i = 1, nl
    DO il = 1, ncum
      IF (ep(il,inb(il))>=0.0001 .AND. i<=inb(il) .AND. iflag(il)<=1) THEN
        IF (cvflag_ice) THEN
          Vprecip(il, i) = wt(il, i)*sigd(il)*(water(il,i)+ice(il,i))/grav
          Vprecipi(il, i) = wt(il, i)*sigd(il)*ice(il,i)/grav                   ! jyg
        ELSE
          Vprecip(il, i) = wt(il, i)*sigd(il)*water(il, i)/grav
          Vprecipi(il, i) = 0.                                                  ! jyg
        END IF
      END IF
    END DO
  END DO


! ***  Calculate downdraft velocity scale    ***
! ***  NE PAS UTILISER POUR L'INSTANT ***

!!      do il=1,ncum
!!        wd(il)=betad*abs(mp(il,icb(il)))*0.01*rrd*t(il,icb(il)) &
!!                                       /(sigd(il)*p(il,icb(il)))
!!      enddo


! ***  calculate tendencies of lowest level potential temperature  ***
! ***                      and mixing ratio                        ***

  DO il = 1, ncum
    work(il) = 1.0/(ph(il,1)-ph(il,2))
    cbmf(il) = 0.0
  END DO

! - Adiabatic ascent mass flux "ma" and cloud base mass flux "cbmf"
!-----------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (adiab_ascent_mass_flux_depends_on_ejectliq) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Warning : this option leads to water conservation violation
!!!           Expert only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO il = 1, ncum
    ma(il, nlp) = 0.
    ma(il, 1)   = 0.
  END DO
  DO k = nl, 2, -1
    DO il = 1, ncum
      ma(il, k) = ma(il, k+1)*(1.-qta(il, k))/(1.-qta(il, k-1)) + m(il, k)
      cbmf(il) = max(cbmf(il), ma(il,k))
    END DO
  END DO
  DO k = 2,nl
    DO il = 1, ncum
      IF (k <icb(il)) THEN
        ma(il, k) = ma(il, k-1) + wghti(il,k-1)*cbmf(il)
      ENDIF
    END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ELSE ! (adiab_ascent_mass_flux_depends_on_ejectliq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Line kept for compatibility with earlier versions
  DO k = 2, nl
    DO il = 1, ncum
      IF (k>=icb(il)) THEN
        cbmf(il) = cbmf(il) + m(il, k)
      END IF
    END DO
  END DO

  DO il = 1, ncum
    ma(il, nlp) = 0.
    ma(il, 1)   = 0.
  END DO
  DO k = nl, 2, -1
    DO il = 1, ncum
      ma(il, k) = ma(il, k+1) + m(il, k)
    END DO
  END DO
  DO k = 2,nl
    DO il = 1, ncum
      IF (k <icb(il)) THEN
        ma(il, k) = ma(il, k-1) + wghti(il,k-1)*cbmf(il)
      ENDIF
    END DO
  END DO

  ENDIF ! (adiab_ascent_mass_flux_depends_on_ejectliq) ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    print*,'cv3_yield avant ft'
! am is the part of cbmf taken from the first level
  DO il = 1, ncum
    am(il) = cbmf(il)*wghti(il, 1)
  END DO

  DO il = 1, ncum
    IF (iflag(il)<=1) THEN
! convect3      if((0.1*dpinv*am).ge.delti)iflag(il)=4
!JYG  Correction pour conserver l'eau
! cc       ft(il,1)=-0.5*lvcp(il,1)*sigd(il)*(evap(il,1)+evap(il,2))          !precip
      IF (cvflag_ice) THEN
        ft(il, 1) = -lvcp(il, 1)*sigd(il)*evap(il, 1) - &
                     lfcp(il, 1)*sigd(il)*evap(il, 1)*faci(il, 1) - &
                     lfcp(il, 1)*sigd(il)*(fondue(il,1)*wt(il,1)) / &
                       (100.*(ph(il,1)-ph(il,2)))                             !precip
      ELSE
        ft(il, 1) = -lvcp(il, 1)*sigd(il)*evap(il, 1)
      END IF

      ft(il, 1) = ft(il, 1) - 0.009*grav*sigd(il)*mp(il, 2)*t_wake(il, 1)*b(il, 1)*work(il)

      IF (cvflag_ice) THEN
        ft(il, 1) = ft(il, 1) + 0.01*sigd(il)*wt(il, 1)*(cl-cpd)*water(il, 2) * &
                                     (t_wake(il,2)-t_wake(il,1))*work(il)/cpn(il, 1) + &
                                0.01*sigd(il)*wt(il, 1)*(ci-cpd)*ice(il, 2) * &
                                     (t_wake(il,2)-t_wake(il,1))*work(il)/cpn(il, 1)
      ELSE
        ft(il, 1) = ft(il, 1) + 0.01*sigd(il)*wt(il, 1)*(cl-cpd)*water(il, 2) * &
                                     (t_wake(il,2)-t_wake(il,1))*work(il)/cpn(il, 1)
      END IF

      ftd(il, 1) = ft(il, 1)                                                  ! fin precip

      IF ((0.01*grav*work(il)*am(il))>=delti) iflag(il) = 1 !consist vect
!jyg<
        IF (fl_cor_ebil >= 2) THEN
          ft(il, 1) = ft(il, 1) + 0.01*grav*work(il)*am(il) * &
                    ((t(il,2)-t(il,1))*cpn(il,2)+gz(il,2)-gz(il,1))/cpn(il,1)
        ELSE
          ft(il, 1) = ft(il, 1) + 0.01*grav*work(il)*am(il) * &
                    (t(il,2)-t(il,1)+(gz(il,2)-gz(il,1))/cpn(il,1))
        ENDIF
!>jyg
    END IF ! iflag
  END DO


  DO j = 2, nl
    IF (iflag_mix>0) THEN
      DO il = 1, ncum
! FH WARNING a modifier :
        cpinv = 0.
! cpinv=1.0/cpn(il,1)
        IF (j<=inb(il) .AND. iflag(il)<=1) THEN
          ft(il, 1) = ft(il, 1) + 0.01*grav*work(il)*ment(il, j, 1) * &
                     (hent(il,j,1)-h(il,1)+t(il,1)*(cpv-cpd)*(rr(il,1)-qent(il,j,1)))*cpinv
        END IF ! j
      END DO
    END IF
  END DO
! fin sature


  DO il = 1, ncum
    IF (iflag(il)<=1) THEN
!JYG1  Correction pour mieux conserver l'eau (conformite avec CONVECT4.3)
      fr(il, 1) = 0.01*grav*mp(il, 2)*(rp(il,2)-rr_wake(il,1))*work(il) + &
                  sigd(il)*evap(il, 1)
!!!                  sigd(il)*0.5*(evap(il,1)+evap(il,2))

      fqd(il, 1) = fr(il, 1) !precip

      fr(il, 1) = fr(il, 1) + 0.01*grav*am(il)*(rr(il,2)-rr(il,1))*work(il)        !sature

      fu(il, 1) = fu(il, 1) + 0.01*grav*work(il)*(mp(il,2)*(up(il,2)-u(il,1)) + &
                                                  am(il)*(u(il,2)-u(il,1)))
      fv(il, 1) = fv(il, 1) + 0.01*grav*work(il)*(mp(il,2)*(vp(il,2)-v(il,1)) + &
                                                  am(il)*(v(il,2)-v(il,1)))
    END IF ! iflag
  END DO ! il


!AC!     do j=1,ntra
!AC!      do il=1,ncum
!AC!       if (iflag(il) .le. 1) then
!AC!       if (cvflag_grav) then
!AC!        ftra(il,1,j)=ftra(il,1,j)+0.01*grav*work(il)
!AC!    :                     *(mp(il,2)*(trap(il,2,j)-tra(il,1,j))
!AC!    :             +am(il)*(tra(il,2,j)-tra(il,1,j)))
!AC!       else
!AC!        ftra(il,1,j)=ftra(il,1,j)+0.1*work(il)
!AC!    :                     *(mp(il,2)*(trap(il,2,j)-tra(il,1,j))
!AC!    :             +am(il)*(tra(il,2,j)-tra(il,1,j)))
!AC!       endif
!AC!       endif  ! iflag
!AC!      enddo
!AC!     enddo

  DO j = 2, nl
    DO il = 1, ncum
      IF (j<=inb(il) .AND. iflag(il)<=1) THEN
        fr(il, 1) = fr(il, 1) + 0.01*grav*work(il)*ment(il, j, 1)*(qent(il,j,1)-rr(il,1))
        fu(il, 1) = fu(il, 1) + 0.01*grav*work(il)*ment(il, j, 1)*(uent(il,j,1)-u(il,1))
        fv(il, 1) = fv(il, 1) + 0.01*grav*work(il)*ment(il, j, 1)*(vent(il,j,1)-v(il,1))
      END IF ! j
    END DO
  END DO

!AC!      do k=1,ntra
!AC!       do j=2,nl
!AC!        do il=1,ncum
!AC!         if (j.le.inb(il) .and. iflag(il) .le. 1) then
!AC!
!AC!          if (cvflag_grav) then
!AC!           ftra(il,1,k)=ftra(il,1,k)+0.01*grav*work(il)*ment(il,j,1)
!AC!     :                *(traent(il,j,1,k)-tra(il,1,k))
!AC!          else
!AC!           ftra(il,1,k)=ftra(il,1,k)+0.1*work(il)*ment(il,j,1)
!AC!     :                *(traent(il,j,1,k)-tra(il,1,k))
!AC!          endif
!AC!
!AC!         endif
!AC!        enddo
!AC!       enddo
!AC!      enddo
! print*,'cv3_yield apres ft'

!jyg<
!-----------------------------------------------------------
           IF (ok_optim_yield) THEN                       !|
!-----------------------------------------------------------
!
!***                                                      ***
!***    Compute convective mass fluxes upwd and dnwd      ***

!
! =================================================
!              upward fluxes                      |
! ------------------------------------------------
!
upwd(:,:) = 0.
up_to(:,:) = 0.
up_from(:,:) = 0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (adiab_ascent_mass_flux_depends_on_ejectliq) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The decrease of the adiabatic ascent mass flux due to ejection of precipitation 
!! is taken into account. 
!! WARNING : in the present version, taking into account the mass-flux decrease due to 
!! precipitation ejection leads to water conservation violation.
!
! - Upward mass flux of mixed draughts
!---------------------------------------
DO i = 2, nl
  DO j = 1, i-1
    DO il = 1, ncum
      IF (i<=inb(il)) THEN
        up_to(il,i) = up_to(il,i) + ment(il,j,i)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!
DO j = 3, nl
  DO i = 2, j-1
    DO il = 1, ncum
      IF (j<=inb(il)) THEN
        up_from(il,i) = up_from(il,i) + ment(il,i,j)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!
! The difference between upwd(il,i) and upwd(il,i-1) is due to updrafts ending in layer 
!(i-1) (theses drafts cross interface (i-1) but not interface(i)) and to updrafts starting 
!from layer (i-1) (theses drafts cross interface (i) but not interface(i-1)): 
! 
DO i = 2, nlp
  DO il = 1, ncum
    IF (i<=inb(il)+1) THEN
      upwd(il,i) = max(0., upwd(il,i-1) - up_to(il,i-1) + up_from(il,i-1))
    ENDIF
  ENDDO
ENDDO
!
! - Total upward mass flux
!---------------------------
DO i = 2, nlp
  DO il = 1, ncum
    IF (i<=inb(il)+1) THEN
      upwd(il,i) = upwd(il,i) + ma(il,i)
    ENDIF
  ENDDO
ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ELSE ! (adiab_ascent_mass_flux_depends_on_ejectliq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The decrease of the adiabatic ascent mass flux due to ejection of precipitation 
!! is not taken into account. 
!
! - Upward mass flux
!-------------------
DO i = 2, nl
  DO il = 1, ncum
    IF (i<=inb(il)) THEN
      up_to(il,i) = m(il,i)
    ENDIF
  ENDDO
  DO j = 1, i-1
    DO il = 1, ncum
      IF (i<=inb(il)) THEN
        up_to(il,i) = up_to(il,i) + ment(il,j,i)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!
DO i = 1, nl
  DO il = 1, ncum
    IF (i<=inb(il)) THEN
      up_from(il,i) = cbmf(il)*wghti(il,i)
    ENDIF
  ENDDO
ENDDO
!
DO j = 3, nl
  DO i = 2, j-1
    DO il = 1, ncum
      IF (j<=inb(il)) THEN
        up_from(il,i) = up_from(il,i) + ment(il,i,j)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!
! The difference between upwd(il,i) and upwd(il,i-1) is due to updrafts ending in layer 
!(i-1) (theses drafts cross interface (i-1) but not interface(i)) and to updrafts starting 
!from layer (i-1) (theses drafts cross interface (i) but not interface(i-1)): 
! 
DO i = 2, nlp
  DO il = 1, ncum
    IF (i<=inb(il)+1) THEN
      upwd(il,i) = max(0., upwd(il,i-1) - up_to(il,i-1) + up_from(il,i-1))
    ENDIF
  ENDDO
ENDDO


  ENDIF ! (adiab_ascent_mass_flux_depends_on_ejectliq) ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! =================================================
!              downward fluxes                    |
! ------------------------------------------------
dnwd(:,:) = 0.
dn_to(:,:) = 0.
dn_from(:,:) = 0.
DO i = 1, nl
  DO j = i+1, nl
    DO il = 1, ncum
      IF (j<=inb(il)) THEN
        dn_to(il,i) = dn_to(il,i) + ment(il,j,i)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!
DO j = 1, nl
  DO i = j+1, nl
    DO il = 1, ncum
      IF (i<=inb(il)) THEN
        dn_from(il,i) = dn_from(il,i) + ment(il,i,j)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!
! The difference between dnwd(il,i) and dnwd(il,i+1) is due to downdrafts ending in layer 
!(i) (theses drafts cross interface (i+1) but not interface(i)) and to downdrafts 
!starting from layer (i) (theses drafts cross interface (i) but not interface(i+1)): 
!
DO i = nl-1, 1, -1
  DO il = 1, ncum
    dnwd(il,i) = max(0., dnwd(il,i+1) - dn_to(il,i) + dn_from(il,i))
  ENDDO
ENDDO
! =================================================
!
!-----------------------------------------------------------
        ENDIF !(ok_optim_yield)                           !|
!-----------------------------------------------------------
!>jyg

! ***  calculate tendencies of potential temperature and mixing ratio  ***
! ***               at levels above the lowest level                   ***

! ***  first find the net saturated updraft and downdraft mass fluxes  ***
! ***                      through each level                          ***


!jyg<
!!  DO i = 2, nl + 1 ! newvecto: mettre nl au lieu nl+1?
  DO i = 2, nl
!>jyg 

    num1 = 0
    DO il = 1, ncum
      IF (i<=inb(il) .AND. iflag(il)<=1) num1 = num1 + 1
    END DO
    IF (num1<=0) GO TO 500

!
!jyg<
!-----------------------------------------------------------
           IF (ok_optim_yield) THEN                       !|
!-----------------------------------------------------------
DO il = 1, ncum
   amp1(il) = upwd(il,i+1)
   ad(il) = dnwd(il,i)
ENDDO
!-----------------------------------------------------------
        ELSE !(ok_optim_yield)                            !|
!-----------------------------------------------------------
!>jyg
    DO il = 1,ncum
      amp1(il) = 0.
      ad(il) = 0.
    ENDDO

    DO k = 1, nl + 1
      DO il = 1, ncum
        IF (i>=icb(il)) THEN
          IF (k>=i+1 .AND. k<=(inb(il)+1)) THEN
            amp1(il) = amp1(il) + m(il, k)
          END IF
        ELSE
! AMP1 is the part of cbmf taken from layers I and lower
          IF (k<=i) THEN
            amp1(il) = amp1(il) + cbmf(il)*wghti(il, k)
          END IF
        END IF
      END DO
    END DO

    DO j = i + 1, nl + 1         
       DO k = 1, i
          !yor! reverted j and k loops 
          DO il = 1, ncum
!yor!        IF (i<=inb(il) .AND. j<=(inb(il)+1)) THEN ! the second condition implies the first !
             IF (j<=(inb(il)+1)) THEN  
                amp1(il) = amp1(il) + ment(il, k, j)
             END IF
          END DO
       END DO
    END DO

    DO k = 1, i - 1
!jyg<
!!      DO j = i, nl + 1 ! newvecto: nl au lieu nl+1?
      DO j = i, nl
!>jyg
        DO il = 1, ncum
!yor!        IF (i<=inb(il) .AND. j<=inb(il)) THEN ! the second condition implies the 1st !
             IF (j<=inb(il)) THEN   
            ad(il) = ad(il) + ment(il, j, k)
          END IF
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------
        ENDIF !(ok_optim_yield)                           !|
!-----------------------------------------------------------
!
!!   print *,'yield, i, amp1, ad', i, amp1(1), ad(1)

    DO il = 1, ncum
      IF (i<=inb(il) .AND. iflag(il)<=1) THEN
        dpinv = 1.0/(ph(il,i)-ph(il,i+1))
        cpinv = 1.0/cpn(il, i)

! convect3      if((0.1*dpinv*amp1).ge.delti)iflag(il)=4
        IF ((0.01*grav*dpinv*amp1(il))>=delti) iflag(il) = 1 ! vecto

! precip
! cc       ft(il,i)= -0.5*sigd(il)*lvcp(il,i)*(evap(il,i)+evap(il,i+1))
        IF (cvflag_ice) THEN
          ft(il, i) = -sigd(il)*lvcp(il, i)*evap(il, i) - &
                       sigd(il)*lfcp(il, i)*evap(il, i)*faci(il, i) - &
                       sigd(il)*lfcp(il, i)*fondue(il, i)*wt(il, i)/(100.*(p(il,i-1)-p(il,i)))
        ELSE
          ft(il, i) = -sigd(il)*lvcp(il, i)*evap(il, i)
        END IF

        rat = cpn(il, i-1)*cpinv

        ft(il, i) = ft(il, i) - 0.009*grav*sigd(il) * &
                     (mp(il,i+1)*t_wake(il,i)*b(il,i)-mp(il,i)*t_wake(il,i-1)*rat*b(il,i-1))*dpinv
        IF (cvflag_ice) THEN
          ft(il, i) = ft(il, i) + 0.01*sigd(il)*wt(il, i)*(cl-cpd)*water(il, i+1) * &
                                       (t_wake(il,i+1)-t_wake(il,i))*dpinv*cpinv + &
                                  0.01*sigd(il)*wt(il, i)*(ci-cpd)*ice(il, i+1) * &
                                       (t_wake(il,i+1)-t_wake(il,i))*dpinv*cpinv
        ELSE
          ft(il, i) = ft(il, i) + 0.01*sigd(il)*wt(il, i)*(cl-cpd)*water(il, i+1) * &
                                       (t_wake(il,i+1)-t_wake(il,i))*dpinv* &
            cpinv
        END IF

        ftd(il, i) = ft(il, i)
! fin precip

! sature
!jyg<
        IF (fl_cor_ebil >= 2) THEN
          ft(il, i) = ft(il, i) + 0.01*grav*dpinv * &
              ( amp1(il)*( (t(il,i+1)-t(il,i))*cpn(il,i+1) + gz(il,i+1)-gz(il,i))*cpinv - &
                ad(il)*( (t(il,i)-t(il,i-1))*cpn(il,i-1) + gz(il,i)-gz(il,i-1))*cpinv)
        ELSE
          ft(il, i) = ft(il, i) + 0.01*grav*dpinv * &
                     (amp1(il)*(t(il,i+1)-t(il,i) + (gz(il,i+1)-gz(il,i))*cpinv) - &
                      ad(il)*(t(il,i)-t(il,i-1)+(gz(il,i)-gz(il,i-1))*cpinv))
        ENDIF
!>jyg


        IF (iflag_mix==0) THEN
          ft(il, i) = ft(il, i) + 0.01*grav*dpinv*ment(il, i, i)*(hp(il,i)-h(il,i) + &
                                    t(il,i)*(cpv-cpd)*(rr(il,i)-qent(il,i,i)))*cpinv
        END IF
!
! sb: on ne fait pas encore la correction permettant de mieux
! conserver l'eau:
!JYG: correction permettant de mieux conserver l'eau:
! cc         fr(il,i)=0.5*sigd(il)*(evap(il,i)+evap(il,i+1))
        fr(il, i) = sigd(il)*evap(il, i) + 0.01*grav*(mp(il,i+1)*(rp(il,i+1)-rr_wake(il,i)) - &
                                                      mp(il,i)*(rp(il,i)-rr_wake(il,i-1)))*dpinv
        fqd(il, i) = fr(il, i)                                                                     ! precip

        fu(il, i) = 0.01*grav*(mp(il,i+1)*(up(il,i+1)-u(il,i)) - &
                               mp(il,i)*(up(il,i)-u(il,i-1)))*dpinv
        fv(il, i) = 0.01*grav*(mp(il,i+1)*(vp(il,i+1)-v(il,i)) - &
                               mp(il,i)*(vp(il,i)-v(il,i-1)))*dpinv


        fr(il, i) = fr(il, i) + 0.01*grav*dpinv*(amp1(il)*(rr(il,i+1)-rr(il,i)) - &
                                                 ad(il)*(rr(il,i)-rr(il,i-1)))
        fu(il, i) = fu(il, i) + 0.01*grav*dpinv*(amp1(il)*(u(il,i+1)-u(il,i)) - &
                                                 ad(il)*(u(il,i)-u(il,i-1)))
        fv(il, i) = fv(il, i) + 0.01*grav*dpinv*(amp1(il)*(v(il,i+1)-v(il,i)) - &
                                                 ad(il)*(v(il,i)-v(il,i-1)))

      END IF ! i
    END DO

!AC!      do k=1,ntra
!AC!       do il=1,ncum
!AC!        if (i.le.inb(il) .and. iflag(il) .le. 1) then
!AC!         dpinv=1.0/(ph(il,i)-ph(il,i+1))
!AC!         cpinv=1.0/cpn(il,i)
!AC!         if (cvflag_grav) then
!AC!           ftra(il,i,k)=ftra(il,i,k)+0.01*grav*dpinv
!AC!     :         *(amp1(il)*(tra(il,i+1,k)-tra(il,i,k))
!AC!     :           -ad(il)*(tra(il,i,k)-tra(il,i-1,k)))
!AC!         else
!AC!           ftra(il,i,k)=ftra(il,i,k)+0.1*dpinv
!AC!     :         *(amp1(il)*(tra(il,i+1,k)-tra(il,i,k))
!AC!     :           -ad(il)*(tra(il,i,k)-tra(il,i-1,k)))
!AC!         endif
!AC!        endif
!AC!       enddo
!AC!      enddo

    DO k = 1, i - 1

      DO il = 1, ncum
        awat(il) = elij(il, k, i) - (1.-ep(il,i))*clw(il, i)
        awat(il) = max(awat(il), 0.0)
      END DO

      IF (iflag_mix/=0) THEN
        DO il = 1, ncum
          IF (i<=inb(il) .AND. iflag(il)<=1) THEN
            dpinv = 1.0/(ph(il,i)-ph(il,i+1))
            cpinv = 1.0/cpn(il, i)
            ft(il, i) = ft(il, i) + 0.01*grav*dpinv*ment(il, k, i) * &
                 (hent(il,k,i)-h(il,i)+t(il,i)*(cpv-cpd)*(rr(il,i)+awat(il)-qent(il,k,i)))*cpinv
!
!
          END IF ! i
        END DO
      END IF

      DO il = 1, ncum
        IF (i<=inb(il) .AND. iflag(il)<=1) THEN
          dpinv = 1.0/(ph(il,i)-ph(il,i+1))
          cpinv = 1.0/cpn(il, i)
          fr(il, i) = fr(il, i) + 0.01*grav*dpinv*ment(il, k, i) * &
                                                       (qent(il,k,i)-awat(il)-rr(il,i))
          fu(il, i) = fu(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(uent(il,k,i)-u(il,i))
          fv(il, i) = fv(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(vent(il,k,i)-v(il,i))

! (saturated updrafts resulting from mixing)                                   ! cld
          qcond(il, i) = qcond(il, i) + (elij(il,k,i)-awat(il))                ! cld
          qtment(il, i) = qtment(il, i) + qent(il,k,i)                         ! cld
          nqcond(il, i) = nqcond(il, i) + 1.                                   ! cld
        END IF ! i
      END DO
    END DO

!AC!      do j=1,ntra
!AC!       do k=1,i-1
!AC!        do il=1,ncum
!AC!         if (i.le.inb(il) .and. iflag(il) .le. 1) then
!AC!          dpinv=1.0/(ph(il,i)-ph(il,i+1))
!AC!          cpinv=1.0/cpn(il,i)
!AC!          if (cvflag_grav) then
!AC!           ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv*ment(il,k,i)
!AC!     :        *(traent(il,k,i,j)-tra(il,i,j))
!AC!          else
!AC!           ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv*ment(il,k,i)
!AC!     :        *(traent(il,k,i,j)-tra(il,i,j))
!AC!          endif
!AC!         endif
!AC!        enddo
!AC!       enddo
!AC!      enddo

!jyg<
!!    DO k = i, nl + 1
    DO k = i, nl
!>jyg

      IF (iflag_mix/=0) THEN
        DO il = 1, ncum
          IF (i<=inb(il) .AND. k<=inb(il) .AND. iflag(il)<=1) THEN
            dpinv = 1.0/(ph(il,i)-ph(il,i+1))
            cpinv = 1.0/cpn(il, i)
            ft(il, i) = ft(il, i) + 0.01*grav*dpinv*ment(il, k, i) * &
                  (hent(il,k,i)-h(il,i)+t(il,i)*(cpv-cpd)*(rr(il,i)-qent(il,k,i)))*cpinv


          END IF ! i
        END DO
      END IF

      DO il = 1, ncum
        IF (i<=inb(il) .AND. k<=inb(il) .AND. iflag(il)<=1) THEN
          dpinv = 1.0/(ph(il,i)-ph(il,i+1))
          cpinv = 1.0/cpn(il, i)

          fr(il, i) = fr(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(qent(il,k,i)-rr(il,i))
          fu(il, i) = fu(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(uent(il,k,i)-u(il,i))
          fv(il, i) = fv(il, i) + 0.01*grav*dpinv*ment(il, k, i)*(vent(il,k,i)-v(il,i))
        END IF ! i and k
      END DO
    END DO

!AC!      do j=1,ntra
!AC!       do k=i,nl+1
!AC!        do il=1,ncum
!AC!         if (i.le.inb(il) .and. k.le.inb(il)
!AC!     $                .and. iflag(il) .le. 1) then
!AC!          dpinv=1.0/(ph(il,i)-ph(il,i+1))
!AC!          cpinv=1.0/cpn(il,i)
!AC!          if (cvflag_grav) then
!AC!           ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv*ment(il,k,i)
!AC!     :         *(traent(il,k,i,j)-tra(il,i,j))
!AC!          else
!AC!           ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv*ment(il,k,i)
!AC!     :             *(traent(il,k,i,j)-tra(il,i,j))
!AC!          endif
!AC!         endif ! i and k
!AC!        enddo
!AC!       enddo
!AC!      enddo

! sb: interface with the cloud parameterization:                               ! cld

    DO k = i + 1, nl
      DO il = 1, ncum
        IF (k<=inb(il) .AND. i<=inb(il) .AND. iflag(il)<=1) THEN               ! cld
! (saturated downdrafts resulting from mixing)                                 ! cld
          qcond(il, i) = qcond(il, i) + elij(il, k, i)                         ! cld
          qtment(il, i) = qent(il,k,i) + qtment(il,i)                          ! cld
          nqcond(il, i) = nqcond(il, i) + 1.                                   ! cld
        END IF ! cld
      END DO ! cld
    END DO ! cld

!ym BIG Warning : it seems that the k loop is missing !!!
!ym Strong advice to check this
!ym add a k loop temporary 

! (particular case: no detraining level is found)                              ! cld
! Verif merge Dynamico<<<<<<< .working
    DO il = 1, ncum                                                            ! cld
      IF (i<=inb(il) .AND. nent(il,i)==0 .AND. iflag(il)<=1) THEN              ! cld
        qcond(il, i) = qcond(il, i) + (1.-ep(il,i))*clw(il, i)                 ! cld
!jyg<   Bug correction 20180620
!      PROBLEM: Should not qent(il,i,i) be taken into account even if nent(il,i)/=0?
!!        qtment(il, i) = qent(il,k,i) + qtment(il,i)                            ! cld
        qtment(il, i) = qent(il,i,i) + qtment(il,i)                            ! cld
!>jyg
        nqcond(il, i) = nqcond(il, i) + 1.                                     ! cld
      END IF                                                                   ! cld
    END DO                                                                     ! cld
! Verif merge Dynamico =======
! Verif merge Dynamico     DO k = i + 1, nl
! Verif merge Dynamico       DO il = 1, ncum        !ym k loop added                                    ! cld
! Verif merge Dynamico         IF (i<=inb(il) .AND. nent(il,i)==0 .AND. iflag(il)<=1) THEN              ! cld
! Verif merge Dynamico           qcond(il, i) = qcond(il, i) + (1.-ep(il,i))*clw(il, i)                 ! cld
! Verif merge Dynamico           qtment(il, i) = qent(il,k,i) + qtment(il,i)                          ! cld
! Verif merge Dynamico           nqcond(il, i) = nqcond(il, i) + 1.                                     ! cld
! Verif merge Dynamico         END IF                                                                   ! cld
! Verif merge Dynamico       END DO
! Verif merge Dynamico     ENDDO                                                                     ! cld
! Verif merge Dynamico >>>>>>> .merge-right.r3413

    DO il = 1, ncum                                                            ! cld
      IF (i<=inb(il) .AND. nqcond(il,i)/=0 .AND. iflag(il)<=1) THEN            ! cld
        qcond(il, i) = qcond(il, i)/nqcond(il, i)                              ! cld
        qtment(il, i) = qtment(il,i)/nqcond(il, i)                             ! cld
      END IF                                                                   ! cld
    END DO

!AC!      do j=1,ntra
!AC!       do il=1,ncum
!AC!        if (i.le.inb(il) .and. iflag(il) .le. 1) then
!AC!         dpinv=1.0/(ph(il,i)-ph(il,i+1))
!AC!         cpinv=1.0/cpn(il,i)
!AC!
!AC!         if (cvflag_grav) then
!AC!          ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv
!AC!     :     *(mp(il,i+1)*(trap(il,i+1,j)-tra(il,i,j))
!AC!     :     -mp(il,i)*(trap(il,i,j)-trap(il,i-1,j)))
!AC!         else
!AC!          ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv
!AC!     :     *(mp(il,i+1)*(trap(il,i+1,j)-tra(il,i,j))
!AC!     :     -mp(il,i)*(trap(il,i,j)-trap(il,i-1,j)))
!AC!         endif
!AC!        endif ! i
!AC!       enddo
!AC!      enddo


500 END DO

!JYG<
!Conservation de l'eau
!   sumdq = 0.
!   DO k = 1, nl
!     sumdq = sumdq + fr(1, k)*100.*(ph(1,k)-ph(1,k+1))/grav
!   END DO
!   PRINT *, 'cv3_yield, apres 500, sum(dq), precip, somme ', sumdq, Vprecip(1, 1), sumdq + vprecip(1, 1)
!JYG>
! ***   move the detrainment at level inb down to level inb-1   ***
! ***        in such a way as to preserve the vertically        ***
! ***          integrated enthalpy and water tendencies         ***

! Correction bug le 18-03-09
  DO il = 1, ncum
    IF (iflag(il)<=1) THEN
      ax = 0.01*grav*ment(il, inb(il), inb(il))* &
           (hp(il,inb(il))-h(il,inb(il))+t(il,inb(il))*(cpv-cpd)*(rr(il,inb(il))-qent(il,inb(il),inb(il))))/ &
                                (cpn(il,inb(il))*(ph(il,inb(il))-ph(il,inb(il)+1)))
      ft(il, inb(il)) = ft(il, inb(il)) - ax
      ft(il, inb(il)-1) = ft(il, inb(il)-1) + ax*cpn(il, inb(il))*(ph(il,inb(il))-ph(il,inb(il)+1))/ &
                              (cpn(il,inb(il)-1)*(ph(il,inb(il)-1)-ph(il,inb(il))))

      bx = 0.01*grav*ment(il, inb(il), inb(il))*(qent(il,inb(il),inb(il))-rr(il,inb(il)))/ &
                                                 (ph(il,inb(il))-ph(il,inb(il)+1))
      fr(il, inb(il)) = fr(il, inb(il)) - bx
      fr(il, inb(il)-1) = fr(il, inb(il)-1) + bx*(ph(il,inb(il))-ph(il,inb(il)+1))/ &
                                                 (ph(il,inb(il)-1)-ph(il,inb(il)))

      cx = 0.01*grav*ment(il, inb(il), inb(il))*(uent(il,inb(il),inb(il))-u(il,inb(il)))/ &
                                                 (ph(il,inb(il))-ph(il,inb(il)+1))
      fu(il, inb(il)) = fu(il, inb(il)) - cx
      fu(il, inb(il)-1) = fu(il, inb(il)-1) + cx*(ph(il,inb(il))-ph(il,inb(il)+1))/ &
                                                 (ph(il,inb(il)-1)-ph(il,inb(il)))

      dx = 0.01*grav*ment(il, inb(il), inb(il))*(vent(il,inb(il),inb(il))-v(il,inb(il)))/ &
                                                 (ph(il,inb(il))-ph(il,inb(il)+1))
      fv(il, inb(il)) = fv(il, inb(il)) - dx
      fv(il, inb(il)-1) = fv(il, inb(il)-1) + dx*(ph(il,inb(il))-ph(il,inb(il)+1))/ &
                                                 (ph(il,inb(il)-1)-ph(il,inb(il)))
    END IF !iflag
  END DO

!JYG<
!Conservation de l'eau
!   sumdq = 0.
!   DO k = 1, nl
!     sumdq = sumdq + fr(1, k)*100.*(ph(1,k)-ph(1,k+1))/grav
!   END DO
!   PRINT *, 'cv3_yield, apres 503, sum(dq), precip, somme ', sumdq, Vprecip(1, 1), sumdq + vprecip(1, 1)
!JYG>

!AC!      do j=1,ntra
!AC!       do il=1,ncum
!AC!        IF (iflag(il) .le. 1) THEN
!AC!	IF (cvflag_grav) then
!AC!        ex=0.01*grav*ment(il,inb(il),inb(il))
!AC!     :      *(traent(il,inb(il),inb(il),j)-tra(il,inb(il),j))
!AC!     :      /(ph(i l,inb(il))-ph(il,inb(il)+1))
!AC!        ftra(il,inb(il),j)=ftra(il,inb(il),j)-ex
!AC!        ftra(il,inb(il)-1,j)=ftra(il,inb(il)-1,j)
!AC!     :       +ex*(ph(il,inb(il))-ph(il,inb(il)+1))
!AC!     :          /(ph(il,inb(il)-1)-ph(il,inb(il)))
!AC!	else
!AC!        ex=0.1*ment(il,inb(il),inb(il))
!AC!     :      *(traent(il,inb(il),inb(il),j)-tra(il,inb(il),j))
!AC!     :      /(ph(i l,inb(il))-ph(il,inb(il)+1))
!AC!        ftra(il,inb(il),j)=ftra(il,inb(il),j)-ex
!AC!        ftra(il,inb(il)-1,j)=ftra(il,inb(il)-1,j)
!AC!     :       +ex*(ph(il,inb(il))-ph(il,inb(il)+1))
!AC!     :          /(ph(il,inb(il)-1)-ph(il,inb(il)))
!AC!        ENDIF   !cvflag grav
!AC!        ENDIF    !iflag
!AC!       enddo
!AC!      enddo


! ***    homogenize tendencies below cloud base    ***


  DO il = 1, ncum
    asum(il) = 0.0
    bsum(il) = 0.0
    csum(il) = 0.0
    dsum(il) = 0.0
    esum(il) = 0.0
    fsum(il) = 0.0
    gsum(il) = 0.0
    hsum(il) = 0.0
  END DO

!do i=1,nl
!do il=1,ncum
!th_wake(il,i)=t_wake(il,i)*(1000.0/p(il,i))**rdcp
!enddo
!enddo

  DO i = 1, nl
    DO il = 1, ncum
      IF (i<=(icb(il)-1) .AND. iflag(il)<=1) THEN
!jyg  Saturated part : use T profile
        asum(il) = asum(il) + (ft(il,i)-ftd(il,i))*(ph(il,i)-ph(il,i+1))
!jyg<20140311
!Correction pour conserver l eau
        IF (ok_conserv_q) THEN
          bsum(il) = bsum(il) + (fr(il,i)-fqd(il,i))*(ph(il,i)-ph(il,i+1))
          csum(il) = csum(il) + (ph(il,i)-ph(il,i+1))

        ELSE
          bsum(il)=bsum(il)+(fr(il,i)-fqd(il,i))*(lv(il,i)+(cl-cpd)*(t(il,i)-t(il,1)))* &
                            (ph(il,i)-ph(il,i+1))
          csum(il)=csum(il)+(lv(il,i)+(cl-cpd)*(t(il,i)-t(il,1)))* &
                            (ph(il,i)-ph(il,i+1))
        ENDIF ! (ok_conserv_q)
!jyg>
        dsum(il) = dsum(il) + t(il, i)*(ph(il,i)-ph(il,i+1))/th(il, i)
!jyg  Unsaturated part : use T_wake profile
        esum(il) = esum(il) + ftd(il, i)*(ph(il,i)-ph(il,i+1))
!jyg<20140311
!Correction pour conserver l eau
        IF (ok_conserv_q) THEN
          fsum(il) = fsum(il) + fqd(il, i)*(ph(il,i)-ph(il,i+1))
          gsum(il) = gsum(il) + (ph(il,i)-ph(il,i+1))
        ELSE
          fsum(il)=fsum(il)+fqd(il,i)*(lv(il,i)+(cl-cpd)*(t_wake(il,i)-t_wake(il,1)))* &
                            (ph(il,i)-ph(il,i+1))
          gsum(il)=gsum(il)+(lv(il,i)+(cl-cpd)*(t_wake(il,i)-t_wake(il,1)))* &
                            (ph(il,i)-ph(il,i+1))
        ENDIF ! (ok_conserv_q)
!jyg>
        hsum(il) = hsum(il) + t_wake(il, i)*(ph(il,i)-ph(il,i+1))/th_wake(il, i)
      END IF
    END DO
  END DO

!!!!      do 700 i=1,icb(il)-1
  IF (ok_homo_tend) THEN
    DO i = 1, nl
      DO il = 1, ncum
        IF (i<=(icb(il)-1) .AND. iflag(il)<=1) THEN
          ftd(il, i) = esum(il)*t_wake(il, i)/(th_wake(il,i)*hsum(il))
          fqd(il, i) = fsum(il)/gsum(il)
          ft(il, i) = ftd(il, i) + asum(il)*t(il, i)/(th(il,i)*dsum(il))
          fr(il, i) = fqd(il, i) + bsum(il)/csum(il)
        END IF
      END DO
    END DO
  ENDIF

!jyg<
!Conservation de l'eau
!!  sumdq = 0.
!!  DO k = 1, nl
!!    sumdq = sumdq + fr(1, k)*100.*(ph(1,k)-ph(1,k+1))/grav
!!  END DO
!!  PRINT *, 'cv3_yield, apres hom, sum(dq), precip, somme ', sumdq, Vprecip(1, 1), sumdq + vprecip(1, 1)
!jyg>


! ***   Check that moisture stays positive. If not, scale tendencies
! in order to ensure moisture positivity
  DO il = 1, ncum
    alpha_qpos(il) = 1.
    IF (iflag(il)<=1) THEN
      IF (fr(il,1)<=0.) THEN
        alpha_qpos(il) = max(alpha_qpos(il), (-delt*fr(il,1))/(s_wake(il)*rr_wake(il,1)+(1.-s_wake(il))*rr(il,1)))
      END IF
    END IF
  END DO
  DO i = 2, nl
    DO il = 1, ncum
      IF (iflag(il)<=1) THEN
        IF (fr(il,i)<=0.) THEN
          alpha_qpos1(il) = max(1., (-delt*fr(il,i))/(s_wake(il)*rr_wake(il,i)+(1.-s_wake(il))*rr(il,i)))
          IF (alpha_qpos1(il)>=alpha_qpos(il)) alpha_qpos(il) = alpha_qpos1(il)
        END IF
      END IF
    END DO
  END DO
  DO il = 1, ncum
    IF (iflag(il)<=1 .AND. alpha_qpos(il)>1.001) THEN
      alpha_qpos(il) = alpha_qpos(il)*1.1
    END IF
  END DO
!
    IF (prt_level .GE. 5) THEN
      print *,' CV3_YIELD : alpha_qpos ',alpha_qpos(1)
    ENDIF
!
  DO il = 1, ncum
    IF (iflag(il)<=1) THEN
      sigd(il) = sigd(il)/alpha_qpos(il)
      precip(il) = precip(il)/alpha_qpos(il)
      cbmf(il) = cbmf(il)/alpha_qpos(il)
    END IF
  END DO
  DO i = 1, nl
    DO il = 1, ncum
      IF (iflag(il)<=1) THEN
        fr(il, i) = fr(il, i)/alpha_qpos(il)
        ft(il, i) = ft(il, i)/alpha_qpos(il)
        fqd(il, i) = fqd(il, i)/alpha_qpos(il)
        ftd(il, i) = ftd(il, i)/alpha_qpos(il)
        fu(il, i) = fu(il, i)/alpha_qpos(il)
        fv(il, i) = fv(il, i)/alpha_qpos(il)
        m(il, i) = m(il, i)/alpha_qpos(il)
        mp(il, i) = mp(il, i)/alpha_qpos(il)
        Vprecip(il, i) = Vprecip(il, i)/alpha_qpos(il)
        Vprecipi(il, i) = Vprecipi(il, i)/alpha_qpos(il)                     ! jyg
      END IF
    END DO
  END DO
!jyg<
!-----------------------------------------------------------
           IF (ok_optim_yield) THEN                       !|
!-----------------------------------------------------------
  DO i = 1, nl
    DO il = 1, ncum
      IF (iflag(il)<=1) THEN
        upwd(il, i) = upwd(il, i)/alpha_qpos(il)
        dnwd(il, i) = dnwd(il, i)/alpha_qpos(il)
      END IF
    END DO
  END DO
!-----------------------------------------------------------
        ENDIF !(ok_optim_yield)                           !|
!-----------------------------------------------------------
!>jyg
  DO j = 1, nl !yor! inverted i and j loops
     DO i = 1, nl
      DO il = 1, ncum
        IF (iflag(il)<=1) THEN
          ment(il, i, j) = ment(il, i, j)/alpha_qpos(il)
        END IF
      END DO
    END DO
  END DO

!AC!      DO j = 1,ntra
!AC!      DO i = 1,nl
!AC!       DO il = 1,ncum
!AC!        IF (iflag(il) .le. 1) THEN
!AC!         ftra(il,i,j) = ftra(il,i,j)/alpha_qpos(il)
!AC!        ENDIF
!AC!       ENDDO
!AC!      ENDDO
!AC!      ENDDO


! ***           reset counter and return           ***

! Reset counter only for points actually convective (jyg)
! In order take into account the possibility of changing the compression,
! reset m, sig and w0 to zero for non-convecting points.
  DO il = 1, ncum
    IF (iflag(il) < 3) THEN
      sig(il, nd) = 2.0
    ENDIF
  END DO


  DO i = 1, nl
    DO il = 1, ncum
      dnwd0(il, i) = -mp(il, i)
    END DO
  END DO
!jyg<  (loops stop at nl)
!!  DO i = nl + 1, nd
!!    DO il = 1, ncum
!!      dnwd0(il, i) = 0.
!!    END DO
!!  END DO
!>jyg


!jyg<
!-----------------------------------------------------------
           IF (.NOT.ok_optim_yield) THEN                  !|
!-----------------------------------------------------------
  DO i = 1, nl
    DO il = 1, ncum
      upwd(il, i) = 0.0
      dnwd(il, i) = 0.0
    END DO
  END DO

!!  DO i = 1, nl                                           ! useless; jyg
!!    DO il = 1, ncum                                      ! useless; jyg
!!      IF (i>=icb(il) .AND. i<=inb(il)) THEN              ! useless; jyg
!!        upwd(il, i) = 0.0                                ! useless; jyg
!!        dnwd(il, i) = 0.0                                ! useless; jyg
!!      END IF                                             ! useless; jyg
!!    END DO                                               ! useless; jyg
!!  END DO                                                 ! useless; jyg

  DO i = 1, nl
    DO k = 1, nl
      DO il = 1, ncum
        up1(il, k, i) = 0.0
        dn1(il, k, i) = 0.0
      END DO
    END DO
  END DO

!yor! commented original
!  DO i = 1, nl
!    DO k = i, nl
!      DO n = 1, i - 1
!        DO il = 1, ncum
!          IF (i>=icb(il) .AND. i<=inb(il) .AND. k<=inb(il)) THEN
!            up1(il, k, i) = up1(il, k, i) + ment(il, n, k)
!            dn1(il, k, i) = dn1(il, k, i) - ment(il, k, n)
!          END IF
!        END DO
!      END DO
!    END DO
!  END DO
!yor! replaced with
  DO i = 1, nl
    DO k = i, nl
      DO n = 1, i - 1
        DO il = 1, ncum
          IF (i>=icb(il) .AND. k<=inb(il)) THEN ! yor ! as i always <= k
             up1(il, k, i) = up1(il, k, i) + ment(il, n, k)
          END IF
        END DO
      END DO
    END DO
  END DO
  DO i = 1, nl
    DO n = 1, i - 1
      DO k = i, nl
        DO il = 1, ncum
          IF (i>=icb(il) .AND. k<=inb(il)) THEN ! yor !  i always <= k
             dn1(il, k, i) = dn1(il, k, i) - ment(il, k, n)
          END IF
        END DO
      END DO
    END DO
  END DO
!yor! end replace 

  DO i = 1, nl
    DO k = 1, nl
      DO il = 1, ncum
        IF (i>=icb(il)) THEN
          IF (k>=i .AND. k<=(inb(il))) THEN
            upwd(il, i) = upwd(il, i) + m(il, k)
          END IF
        ELSE
          IF (k<i) THEN
            upwd(il, i) = upwd(il, i) + cbmf(il)*wghti(il, k)
          END IF
        END IF
! c        print *,'cbmf',il,i,k,cbmf(il),wghti(il,k)
      END DO
    END DO
  END DO

  DO i = 2, nl
    DO k = i, nl
      DO il = 1, ncum
! test         if (i.ge.icb(il).and.i.le.inb(il).and.k.le.inb(il)) then
        IF (i<=inb(il) .AND. k<=inb(il)) THEN
          upwd(il, i) = upwd(il, i) + up1(il, k, i)
          dnwd(il, i) = dnwd(il, i) + dn1(il, k, i)
        END IF
! c         print *,'upwd',il,i,k,inb(il),upwd(il,i),m(il,k),up1(il,k,i)
      END DO
    END DO
  END DO


!!!!      DO il=1,ncum
!!!!      do i=icb(il),inb(il)
!!!!
!!!!      upwd(il,i)=0.0
!!!!      dnwd(il,i)=0.0
!!!!      do k=i,inb(il)
!!!!      up1=0.0
!!!!      dn1=0.0
!!!!      do n=1,i-1
!!!!      up1=up1+ment(il,n,k)
!!!!      dn1=dn1-ment(il,k,n)
!!!!      enddo
!!!!      upwd(il,i)=upwd(il,i)+m(il,k)+up1
!!!!      dnwd(il,i)=dnwd(il,i)+dn1
!!!!      enddo
!!!!      enddo
!!!!
!!!!      ENDDO

!!  DO i = 1, nlp
!!    DO il = 1, ncum
!!      ma(il, i) = 0
!!    END DO
!!  END DO
!!
!!  DO i = 1, nl
!!    DO j = i, nl
!!      DO il = 1, ncum
!!        ma(il, i) = ma(il, i) + m(il, j)
!!      END DO
!!    END DO
!!  END DO

!jyg<  (loops stop at nl)
!!  DO i = nl + 1, nd
!!    DO il = 1, ncum
!!      ma(il, i) = 0.
!!    END DO
!!  END DO
!>jyg

!!  DO i = 1, nl
!!    DO il = 1, ncum
!!      IF (i<=(icb(il)-1)) THEN
!!        ma(il, i) = 0
!!      END IF
!!    END DO
!!  END DO

!-----------------------------------------------------------
        ENDIF !(.NOT.ok_optim_yield)                      !|
!-----------------------------------------------------------
!>jyg

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! determination de la variation de flux ascendant entre
! deux niveau non dilue mip
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  DO i = 1, nl
    DO il = 1, ncum
      mip(il, i) = m(il, i)
    END DO
  END DO

!jyg<  (loops stop at nl)
!!  DO i = nl + 1, nd
!!    DO il = 1, ncum
!!      mip(il, i) = 0.
!!    END DO
!!  END DO
!>jyg


! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! icb represente de niveau ou se trouve la
! base du nuage , et inb le top du nuage
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!  DO i = 1, nd                                  ! unused . jyg
!!    DO il = 1, ncum                             ! unused . jyg
!!      mke(il, i) = upwd(il, i) + dnwd(il, i)    ! unused . jyg
!!    END DO                                      ! unused . jyg
!!  END DO                                        ! unused . jyg

!!  DO i = 1, nd                                                                 ! unused . jyg
!!    DO il = 1, ncum                                                            ! unused . jyg
!!      rdcp = (rrd*(1.-rr(il,i))-rr(il,i)*rrv)/(cpd*(1.-rr(il,i))+rr(il,i)*cpv) ! unused . jyg
!!      tls(il, i) = t(il, i)*(1000.0/p(il,i))**rdcp                             ! unused . jyg
!!      tps(il, i) = tp(il, i)                                                   ! unused . jyg
!!    END DO                                                                     ! unused . jyg
!!  END DO                                                                       ! unused . jyg


! *** diagnose the in-cloud mixing ratio   ***                       ! cld
! ***           of condensed water         ***                       ! cld
!! cld                                                               
                                                                     
  DO i = 1, nl+1                                                     ! cld
    DO il = 1, ncum                                                  ! cld
      mac(il, i) = 0.0                                               ! cld
      wa(il, i) = 0.0                                                ! cld
      siga(il, i) = 0.0                                              ! cld
      sax(il, i) = 0.0                                               ! cld
    END DO                                                           ! cld
  END DO                                                             ! cld
                                                                     
  DO i = minorig, nl                                                 ! cld
    DO k = i + 1, nl + 1                                             ! cld
      DO il = 1, ncum                                                ! cld
        IF (i<=inb(il) .AND. k<=(inb(il)+1) .AND. iflag(il)<=1) THEN ! cld
          mac(il, i) = mac(il, i) + m(il, k)                         ! cld
        END IF                                                       ! cld
      END DO                                                         ! cld
    END DO                                                           ! cld
  END DO                                                             ! cld

  DO i = 1, nl                                                       ! cld
    DO j = 1, i                                                      ! cld
      DO il = 1, ncum                                                ! cld
        IF (i>=icb(il) .AND. i<=(inb(il)-1) &                        ! cld
            .AND. j>=icb(il) .AND. iflag(il)<=1) THEN                ! cld
          sax(il, i) = sax(il, i) + rrd*(tvp(il,j)-tv(il,j)) &       ! cld
            *(ph(il,j)-ph(il,j+1))/p(il, j)                          ! cld
        END IF                                                       ! cld
      END DO                                                         ! cld
    END DO                                                           ! cld
  END DO                                                             ! cld

  DO i = 1, nl                                                       ! cld
    DO il = 1, ncum                                                  ! cld
      IF (i>=icb(il) .AND. i<=(inb(il)-1) &                          ! cld
          .AND. sax(il,i)>0.0 .AND. iflag(il)<=1) THEN               ! cld
        wa(il, i) = sqrt(2.*sax(il,i))                               ! cld
      END IF                                                         ! cld
    END DO                                                           ! cld
  END DO  
                                                           ! cld
  DO i = 1, nl  

! 14/01/15 AJ je remets les parties manquantes cf JYG
! Initialize sument to 0

    DO il = 1,ncum
     sument(il) = 0.
    ENDDO

! Sum mixed mass fluxes in sument

    DO k = 1,nl
      DO il = 1,ncum
        IF  (k<=inb(il) .AND. i<=inb(il) .AND. iflag(il)<=1) THEN   ! cld
          sument(il) =sument(il) + abs(ment(il,k,i))
        ENDIF
      ENDDO     ! il
    ENDDO       ! k

! 14/01/15 AJ delta n'a rien à faire là...                                                 
    DO il = 1, ncum                                                  ! cld
!!      IF (wa(il,i)>0.0 .AND. iflag(il)<=1) &                         ! cld
!!        siga(il, i) = mac(il, i)/(coefw_cld_cv*wa(il, i)) &          ! cld
!!        *rrd*tvp(il, i)/p(il, i)/100.                                ! cld
!!
!!      siga(il, i) = min(siga(il,i), 1.0)                             ! cld
      sigaq = 0.
      IF (wa(il,i)>0.0 .AND. iflag(il)<=1)  THEN                     ! cld
        siga(il, i) = mac(il, i)/(coefw_cld_cv*wa(il, i)) &          ! cld
                     *rrd*tvp(il, i)/p(il, i)/100.                   ! cld
        siga(il, i) = min(siga(il,i), 1.0)                           ! cld
        sigaq = siga(il,i)*qta(il,i-1)                               ! cld
      ENDIF

! IM cf. FH 
! 14/01/15 AJ ne correspond pas à ce qui a été codé par JYG et SB           
                                                         
      IF (iflag_clw==0) THEN                                         ! cld
        qcondc(il, i) = siga(il, i)*clw(il, i)*(1.-ep(il,i)) &       ! cld
          +(1.-siga(il,i))*qcond(il, i)                              ! cld


        sigment(il,i)=sument(il)*tau_cld_cv/(ph(il,i)-ph(il,i+1))    ! cld
        sigment(il, i) = min(1.e-4+sigment(il,i), 1.0 - siga(il,i))  ! cld
!!        qtc(il, i) = (siga(il,i)*qta(il,i-1)+sigment(il,i)*qtment(il,i)) & ! cld
        qtc(il, i) = (sigaq+sigment(il,i)*qtment(il,i)) & ! cld
                     /(siga(il,i)+sigment(il,i))                     ! cld
        sigt(il,i) = sigment(il, i) + siga(il, i)

!        qtc(il, i) = siga(il,i)*qta(il,i-1)+(1.-siga(il,i))*qtment(il,i) ! cld
!     print*,'BIGAUSSIAN CONV',siga(il,i),sigment(il,i),qtc(il,i)  
      		
      ELSE IF (iflag_clw==1) THEN                                    ! cld
        qcondc(il, i) = qcond(il, i)                                 ! cld
        qtc(il,i) = qtment(il,i)                                     ! cld
      END IF                                                         ! cld

    END DO                                                           ! cld
  END DO
! print*,'cv3_yield fin'

  RETURN
END SUBROUTINE cv3_yield

!AC! et !RomP >>>
SUBROUTINE cv3_tracer(nloc, len, ncum, nd, na, &
                      ment, sigij, da, phi, phi2, d1a, dam, &
                      ep, Vprecip, elij, clw, epmlmMm, eplaMm, &
                      icb, inb)
  IMPLICIT NONE

  include "cv3param.h"

!inputs:
  INTEGER, INTENT (IN)                               :: ncum, nd, na, nloc, len
  INTEGER, DIMENSION (len), INTENT (IN)              :: icb, inb
  REAL, DIMENSION (len, na, na), INTENT (IN)         :: ment, sigij, elij
  REAL, DIMENSION (len, nd), INTENT (IN)             :: clw
  REAL, DIMENSION (len, na), INTENT (IN)             :: ep
  REAL, DIMENSION (len, nd+1), INTENT (IN)           :: Vprecip
!ouputs:
  REAL, DIMENSION (len, na, na), INTENT (OUT)        :: phi, phi2, epmlmMm
  REAL, DIMENSION (len, na), INTENT (OUT)            :: da, d1a, dam, eplaMm
!
! variables pour tracer dans precip de l'AA et des mel
!local variables:
  INTEGER i, j, k
  REAL epm(nloc, na, na)

! variables d'Emanuel : du second indice au troisieme
! --->    tab(i,k,j) -> de l origine k a l arrivee j
! ment, sigij, elij
! variables personnelles : du troisieme au second indice
! --->    tab(i,j,k) -> de k a j
! phi, phi2

! initialisations

  da(:, :) = 0.
  d1a(:, :) = 0.
  dam(:, :) = 0.
  epm(:, :, :) = 0.
  eplaMm(:, :) = 0.
  epmlmMm(:, :, :) = 0.
  phi(:, :, :) = 0.
  phi2(:, :, :) = 0.

! fraction deau condensee dans les melanges convertie en precip : epm
! et eau condensée précipitée dans masse d'air saturé : l_m*dM_m/dzdz.dzdz
  DO j = 1, nl
    DO k = 1, nl
      DO i = 1, ncum
        IF (k>=icb(i) .AND. k<=inb(i) .AND. & 
!!jyg              j.ge.k.and.j.le.inb(i)) then
!!jyg             epm(i,j,k)=1.-(1.-ep(i,j))*clw(i,j)/elij(i,k,j)
            j>k .AND. j<=inb(i)) THEN
          epm(i, j, k) = 1. - (1.-ep(i,j))*clw(i, j)/max(elij(i,k,j), 1.E-16)
!!
          epm(i, j, k) = max(epm(i,j,k), 0.0)
        END IF
      END DO
    END DO
  END DO


  DO j = 1, nl
    DO k = 1, nl
      DO i = 1, ncum
        IF (k>=icb(i) .AND. k<=inb(i)) THEN
          eplaMm(i, j) = eplamm(i, j) + &
                         ep(i, j)*clw(i, j)*ment(i, j, k)*(1.-sigij(i,j,k))
        END IF
      END DO
    END DO
  END DO

  DO j = 1, nl
    DO k = 1, j - 1
      DO i = 1, ncum
        IF (k>=icb(i) .AND. k<=inb(i) .AND. j<=inb(i)) THEN
          epmlmMm(i, j, k) = epm(i, j, k)*elij(i, k, j)*ment(i, k, j)
        END IF
      END DO
    END DO
  END DO

! matrices pour calculer la tendance des concentrations dans cvltr.F90
  DO j = 1, nl
    DO k = 1, nl
      DO i = 1, ncum
        da(i, j) = da(i, j) + (1.-sigij(i,k,j))*ment(i, k, j)
        phi(i, j, k) = sigij(i, k, j)*ment(i, k, j)
        d1a(i, j) = d1a(i, j) + ment(i, k, j)*ep(i, k)*(1.-sigij(i,k,j))
        IF (k<=j) THEN
          dam(i, j) = dam(i, j) + ment(i, k, j)*epm(i, k, j)*(1.-ep(i,k))*(1.-sigij(i,k,j))
          phi2(i, j, k) = phi(i, j, k)*epm(i, j, k)
        END IF
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cv3_tracer
!AC! et !RomP <<<

SUBROUTINE cv3_uncompress(nloc, len, ncum, nd, ntra, idcum, &
                          iflag, &
                          precip, sig, w0, &
                          ft, fq, fu, fv, ftra, &
                          Ma, upwd, dnwd, dnwd0, qcondc, wd, cape, &
                          epmax_diag, & ! epmax_cape
                          iflag1, &
                          precip1, sig1, w01, &
                          ft1, fq1, fu1, fv1, ftra1, &
                          Ma1, upwd1, dnwd1, dnwd01, qcondc1, wd1, cape1, &
                          epmax_diag1) ! epmax_cape
  IMPLICIT NONE

  include "cv3param.h"

!inputs:
  INTEGER len, ncum, nd, ntra, nloc
  INTEGER idcum(nloc)
  INTEGER iflag(nloc)
  REAL precip(nloc)
  REAL sig(nloc, nd), w0(nloc, nd)
  REAL ft(nloc, nd), fq(nloc, nd), fu(nloc, nd), fv(nloc, nd)
  REAL ftra(nloc, nd, ntra)
  REAL ma(nloc, nd)
  REAL upwd(nloc, nd), dnwd(nloc, nd), dnwd0(nloc, nd)
  REAL qcondc(nloc, nd)
  REAL wd(nloc), cape(nloc)
  REAL epmax_diag(nloc)

!outputs:
  INTEGER iflag1(len)
  REAL precip1(len)
  REAL sig1(len, nd), w01(len, nd)
  REAL ft1(len, nd), fq1(len, nd), fu1(len, nd), fv1(len, nd)
  REAL ftra1(len, nd, ntra)
  REAL ma1(len, nd)
  REAL upwd1(len, nd), dnwd1(len, nd), dnwd01(len, nd)
  REAL qcondc1(nloc, nd)
  REAL wd1(nloc), cape1(nloc)
  REAL epmax_diag1(len) ! epmax_cape

!local variables:
  INTEGER i, k, j

  DO i = 1, ncum
    precip1(idcum(i)) = precip(i)
    iflag1(idcum(i)) = iflag(i)
    wd1(idcum(i)) = wd(i)
    cape1(idcum(i)) = cape(i)
    epmax_diag1(idcum(i))=epmax_diag(i) ! epmax_cape
  END DO

  DO k = 1, nl
    DO i = 1, ncum
      sig1(idcum(i), k) = sig(i, k)
      w01(idcum(i), k) = w0(i, k)
      ft1(idcum(i), k) = ft(i, k)
      fq1(idcum(i), k) = fq(i, k)
      fu1(idcum(i), k) = fu(i, k)
      fv1(idcum(i), k) = fv(i, k)
      ma1(idcum(i), k) = ma(i, k)
      upwd1(idcum(i), k) = upwd(i, k)
      dnwd1(idcum(i), k) = dnwd(i, k)
      dnwd01(idcum(i), k) = dnwd0(i, k)
      qcondc1(idcum(i), k) = qcondc(i, k)
    END DO
  END DO

  DO i = 1, ncum
    sig1(idcum(i), nd) = sig(i, nd)
  END DO


!AC!        do 2100 j=1,ntra
!AC!c oct3         do 2110 k=1,nl
!AC!         do 2110 k=1,nd ! oct3
!AC!          do 2120 i=1,ncum
!AC!            ftra1(idcum(i),k,j)=ftra(i,k,j)
!AC! 2120     continue
!AC! 2110    continue
!AC! 2100   continue
!
  RETURN
END SUBROUTINE cv3_uncompress


        subroutine cv3_epmax_fn_cape(nloc,ncum,nd &
                 , ep,hp,icb,inb,clw,nk,t,h,hnk,lv,lf,frac &
                 , pbase, p, ph, tv, buoy, sig, w0,iflag &
                 , epmax_diag)
        implicit none

        ! On fait varier epmax en fn de la cape
        ! Il faut donc recalculer ep, et hp qui a déjà été calculé et
        ! qui en dépend
        ! Toutes les autres variables fn de ep sont calculées plus bas.

  include "cvthermo.h"
  include "cv3param.h"  
  include "conema3.h"
  include "cvflag.h"

! inputs:
      INTEGER, INTENT (IN)                               :: ncum, nd, nloc
      INTEGER, DIMENSION (nloc), INTENT (IN)             :: icb, inb, nk
      REAL, DIMENSION (nloc), INTENT (IN)                :: hnk,pbase
      REAL, DIMENSION (nloc, nd), INTENT (IN)            :: t, lv, lf, tv, h
      REAL, DIMENSION (nloc, nd), INTENT (IN)            :: clw, buoy,frac
      REAL, DIMENSION (nloc, nd), INTENT (IN)            :: sig,w0
      INTEGER, DIMENSION (nloc), INTENT (IN)             :: iflag(nloc)
      REAL, DIMENSION (nloc, nd), INTENT (IN)            :: p
      REAL, DIMENSION (nloc, nd+1), INTENT (IN)          :: ph
! inouts:
      REAL, DIMENSION (nloc, nd), INTENT (INOUT)         :: ep,hp  
! outputs
      REAL, DIMENSION (nloc), INTENT (OUT)           :: epmax_diag

! local
      integer i,k    
!      real hp_bak(nloc,nd)
!      real ep_bak(nloc,nd)
      real m_loc(nloc,nd)
      real sig_loc(nloc,nd)
      real w0_loc(nloc,nd)
      integer iflag_loc(nloc)
      real cape(nloc)
        
        if (coef_epmax_cape.gt.1e-12) then

        ! il faut calculer la cape: on fait un calcule simple car tant qu'on ne
        ! connait pas ep, on ne connait pas les mélanges, ddfts etc... qui sont
        ! necessaires au calcul de la cape dans la nouvelle physique
        
!        write(*,*) 'cv3_routines check 4303'
        do i=1,ncum
        do k=1,nd
          sig_loc(i,k)=sig(i,k)
          w0_loc(i,k)=w0(i,k)
          iflag_loc(i)=iflag(i)
!          ep_bak(i,k)=ep(i,k)
        enddo ! do k=1,nd
        enddo !do i=1,ncum

!        write(*,*) 'cv3_routines check 4311'
!        write(*,*) 'nl=',nl
        CALL cv3_closure(nloc, ncum, nd, icb, inb, & ! na->nd
          pbase, p, ph, tv, buoy, &
          sig_loc, w0_loc, cape, m_loc,iflag_loc)

!        write(*,*) 'cv3_routines check 4316'
!        write(*,*) 'ep(1,:)=',ep(1,:)
        do i=1,ncum
           epmax_diag(i)=epmax-coef_epmax_cape*sqrt(cape(i))
           epmax_diag(i)=amax1(epmax_diag(i),0.0)
!           write(*,*) 'i,icb,inb,cape,epmax_diag=', &
!                i,icb(i),inb(i),cape(i),epmax_diag(i)
           do k=1,nl
                ep(i,k)=ep(i,k)/epmax*epmax_diag(i)
                ep(i,k)=amax1(ep(i,k),0.0)
                ep(i,k)=amin1(ep(i,k),epmax_diag(i))
           enddo
        enddo
 !       write(*,*) 'ep(1,:)=',ep(1,:)

      !write(*,*) 'cv3_routines check 4326'
! On recalcule hp:
!      do k=1,nl
!        do i=1,ncum
!	  hp_bak(i,k)=hp(i,k)
!	enddo
!      enddo
      do k=1,nl
        do i=1,ncum
          hp(i,k)=h(i,k)
        enddo
      enddo

  IF (cvflag_ice) THEN

      do k=minorig+1,nl
       do i=1,ncum
        if((k.ge.icb(i)).and.(k.le.inb(i)))then
          hp(i, k) = hnk(i) + (lv(i,k)+(cpd-cpv)*t(i,k)+frac(i,k)*lf(i,k))* &
                              ep(i, k)*clw(i, k)
        endif
       enddo
      enddo !do k=minorig+1,n
  ELSE !IF (cvflag_ice) THEN

      DO k = minorig + 1, nl
       DO i = 1, ncum
        IF ((k>=icb(i)) .AND. (k<=inb(i))) THEN
          hp(i,k)=hnk(i)+(lv(i,k)+(cpd-cpv)*t(i,k))*ep(i,k)*clw(i,k)
        endif
       enddo
      enddo !do k=minorig+1,n

  ENDIF !IF (cvflag_ice) THEN     
      !write(*,*) 'cv3_routines check 4345'
!      do i=1,ncum  
!       do k=1,nl
!        if ((abs(hp_bak(i,k)-hp(i,k))/hp_bak(i,k).gt.1e-1).or. &
!            ((abs(hp_bak(i,k)-hp(i,k))/hp_bak(i,k).gt.1e-4).and. &
!            (ep(i,k)-ep_bak(i,k).lt.1e-4))) then
!           write(*,*) 'i,k=',i,k
!           write(*,*) 'coef_epmax_cape=',coef_epmax_cape
!           write(*,*) 'epmax_diag(i)=',epmax_diag(i)
!           write(*,*) 'ep(i,k)=',ep(i,k)
!           write(*,*) 'ep_bak(i,k)=',ep_bak(i,k)
!           write(*,*) 'hp(i,k)=',hp(i,k)
!           write(*,*) 'hp_bak(i,k)=',hp_bak(i,k)
!           write(*,*) 'h(i,k)=',h(i,k)
!           write(*,*) 'nk(i)=',nk(i)
!           write(*,*) 'h(i,nk(i))=',h(i,nk(i))
!           write(*,*) 'lv(i,k)=',lv(i,k)
!           write(*,*) 't(i,k)=',t(i,k)
!           write(*,*) 'clw(i,k)=',clw(i,k)
!           write(*,*) 'cpd,cpv=',cpd,cpv
!           stop
!        endif
!       enddo !do k=1,nl
!      enddo !do i=1,ncum  
      endif !if (coef_epmax_cape.gt.1e-12) then
      !write(*,*) 'cv3_routines check 4367'

      return
      end subroutine cv3_epmax_fn_cape



