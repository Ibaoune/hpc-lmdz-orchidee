
! $Header$

SUBROUTINE hbtm2l(knon, paprs, pplay, t2m, t10m, q2m, q10m, ustar, flux_t, flux_q, u, v, t, q, pblh, therm, plcl, cape, &
    cin, eauliq, ctei, d_qt, d_thv, dlt_2, xhis, posint, omega, diagok)
  USE dimphy
  IMPLICIT NONE

  ! ***************************************************************
  ! *                                                             *
  ! * HBTM2L   D'apres Holstag&Boville et Troen&Mahrt             *
  ! *                 JAS 47              BLM                     *
  ! * Algorithmes These Anne Mathieu                              *
  ! * Critere d'Entrainement Peter Duynkerke (JAS 50)             *
  ! * written by  : Anne MATHIEU & Alain LAHELLEC, 22/11/99       *
  ! * features : implem. exces Mathieu                            *
  ! ***************************************************************
  ! * mods : decembre 99 passage th a niveau plus bas. voir fixer *
  ! ***************************************************************
  ! * fin therm a la HBTM passage a forme Mathieu 12/09/2001      *
  ! ***************************************************************
  ! AM Fev 2003 Adaptation a LMDZ version couplee                 *
  ! *
  ! Pour le moment on fait passer en argument les grdeurs de surface :
  ! flux, t,q2m, t,q10m, on va utiliser systematiquement les grdeurs a 2m ms
  ! on garde la possibilite de changer si besoin est (jusqu'a present la
  ! forme de HB avec le 1er niveau modele etait conservee)       *
  ! ***************************************************************
  ! * re-ecriture complete Alain Mars 2012 dans LMDZ5V5           *
  ! ***************************************************************
  include "YOMCST.h"
  REAL rlvcp, reps
  ! Arguments:

  INTEGER knon ! nombre de points a calculer
  ! AM
  REAL t2m(klon), t10m(klon) ! temperature a 2 et 10m
  REAL q2m(klon), q10m(klon) ! q a 2 et 10m
  REAL ustar(klon)
  REAL paprs(klon, klev+1) ! pression a inter-couche (Pa)
  REAL pplay(klon, klev) ! pression au milieu de couche (Pa)
  REAL flux_t(klon, klev), flux_q(klon, klev) ! Flux
  REAL u(klon, klev) ! vitesse U (m/s)
  REAL v(klon, klev) ! vitesse V (m/s)
  REAL t(klon, klev) ! temperature (K)
  REAL q(klon, klev) ! vapeur d'eau (kg/kg)

  INTEGER isommet
  REAL vk
  PARAMETER (vk=0.35) ! Von Karman => passer a .41 ! cf U.Olgstrom
  REAL ricr
  PARAMETER (ricr=0.4)
  REAL fak
  PARAMETER (fak=8.5) ! b calcul du Prandtl et de dTetas
  REAL fakn
  PARAMETER (fakn=7.2) ! a
  REAL onet
  PARAMETER (onet=1.0/3.0)
  REAL betam
  PARAMETER (betam=15.0) ! pour Phim / h dans la S.L stable
  REAL betah
  PARAMETER (betah=15.0)
  REAL betas
  PARAMETER (betas=5.0) ! Phit dans la S.L. stable (mais 2 formes / z/OBL<>1
  REAL sffrac
  PARAMETER (sffrac=0.1) ! S.L. = z/h < .1
  REAL binm
  PARAMETER (binm=betam*sffrac)
  REAL binh
  PARAMETER (binh=betah*sffrac)

  REAL q_star, t_star
  REAL b1, b2, b212, b2sr ! Lambert correlations T' q' avec T* q*
  PARAMETER (b1=70., b2=20.) ! b1 entre 70 et 100

  REAL z(klon, klev)
  ! AM
  REAL zref, dt0
  PARAMETER (zref=2.) ! Niveau de ref a 2m
  PARAMETER (dt0=0.1) ! convergence do while

  INTEGER i, k, j
  REAL khfs(klon) ! surface kinematic heat flux [mK/s]
  REAL kqfs(klon) ! sfc kinematic constituent flux [m/s]
  REAL heatv(klon) ! surface virtual heat flux
  REAL rhino(klon, klev) ! bulk Richardon no. mais en Theta_v
  LOGICAL unstbl(klon) ! pts w/unstbl pbl (positive virtual ht flx)
  LOGICAL check(klon) ! True=>chk if Richardson no.>critcal
  LOGICAL omegafl(klon) ! flag de prolongement cape pour pt Omega
  REAL obklen(klon) ! Monin-Obukhov lengh

  REAL pblh(klon) ! PBL H               (m)
  REAL therm(klon) ! exces du thermique  (K)
  REAL plcl(klon) ! Lifted Cnd Level (Pa)
  REAL cape(klon) ! Cape
  REAL cin(klon) ! Inhibition
  REAL eauliq(klon) ! Eau Liqu integree
  REAL ctei(klon) ! Cld Top Entr. Instab.
  REAL d_qt(klon) ! Saut de qT a l'inversion
  REAL d_thv(klon) !         Theta_e
  REAL dlt_2(klon) ! Ordonnee a gauche de courbe de melange
  REAL xhis(klon) ! fraction de melange pour flottab nulle
  REAL posint(klon) ! partie positive de l'int. de Peter
  REAL omega(klon) ! point ultime de l'ascention du thermique
  REAL diagok(klon) ! pour traiter les sous-mailles sans info
  ! Algorithme thermique
  REAL s(klon, klev) ! [P/Po]^Kappa milieux couches
  REAL th_th(klon) ! potential temperature of thermal
  REAL the_th(klon) ! equivalent potential temperature of thermal
  REAL qt_th(klon) ! total water of thermal
  REAL tbef(klon) ! T thermique niveau ou calcul precedent
  LOGICAL zsat(klon) ! le thermique est sature
  LOGICAL zcin(klon) ! calcul d'inhibition
  REAL kape(klon) ! Cape locale
  REAL kin(klon) ! Cin locale
  ! calcul de CTEI etc
  REAL the1, the2, aa, bb, zthvd, zthvu, qsat, chi, rh, zxt, zdu2
  REAL rnum, denom, th1, th2, tv1, tv2, thv1, thv2, ql1, ql2, dt
  REAL dqsat_dt, qsat2, qt1, q1, q2, t1, t2, tl1, te2, xnull, delt_the
  REAL delt_qt, quadsat, spblh, reduc
  ! diag      REAL dTv21(klon,klev)

  REAL phiminv(klon) ! inverse phi function for momentum
  REAL phihinv(klon) ! inverse phi function for heat
  REAL wm(klon) ! turbulent velocity scale for momentum
  REAL zm(klon) ! current level height
  REAL zp(klon) ! current level height + one level up
  REAL zcor, zdelta, zcvm5
  REAL fac, pblmin
  REAL missing_val

  include "YOETHF.h"
  include "FCTTRE.h"

  ! c      missing_val=nf90_fill_real (avec include netcdf)
  missing_val = 0.

  ! initialisations (Anne)
  isommet = klev
  b212 = sqrt(b1*b2)
  b2sr = sqrt(b2)

  ! Initialisation thermo
  rlvcp = rlvtt/rcpd
  reps = rd/rv
  ! raz
  q_star = 0.
  t_star = 0.
  cape(:) = missing_val
  kape(:) = 0.
  cin(:) = missing_val
  eauliq(:) = missing_val
  ctei(:) = missing_val
  d_qt(:) = missing_val
  d_thv(:) = missing_val
  dlt_2(:) = missing_val
  xhis(:) = missing_val
  posint(:) = missing_val
  kin(:) = missing_val
  omega(:) = missing_val
  diagok(:) = 0.
  ! diag      dTv21(:,:)= missing_val

  ! Calculer les hauteurs de chaque couche
  DO i = 1, knon
    z(i, 1) = rd*t(i, 1)/(0.5*(paprs(i,1)+pplay(i,1)))*(paprs(i,1)-pplay(i,1))/rg
    s(i, 1) = (pplay(i,1)/paprs(i,1))**rkappa
  END DO
  ! s(k) = [pplay(k)/ps]^kappa
  ! + + + + + + + + + pplay  <-> s(k)   t  dp=pplay(k-1)-pplay(k)

  ! -----------------  paprs <-> sig(k)

  ! + + + + + + + + + pplay  <-> s(k-1)


  ! + + + + + + + + + pplay  <-> s(1)   t  dp=paprs-pplay   z(1)

  ! -----------------  paprs <-> sig(1)

  DO k = 2, klev
    DO i = 1, knon
      z(i, k) = z(i, k-1) + rd*0.5*(t(i,k-1)+t(i,k))/paprs(i, k)*(pplay(i,k-1)-pplay(i,k))/rg
      s(i, k) = (pplay(i,k)/paprs(i,1))**rkappa
    END DO
  END DO
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++  Determination des grandeurs de surface  +++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO i = 1, knon
    ! AM Niveau de ref choisi a 2m
    zxt = t2m(i)

    ! ***************************************************
    ! attention, il doit s'agir de <w'theta'>
    ! ;Calcul de tcls virtuel et de w'theta'virtuel
    ! ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ! tcls=tcls*(1+.608*qcls)

    ! ;Pour avoir w'theta',
    ! ; il faut diviser par ro.Cp
    ! Cp=Cpd*(1+0.84*qcls)
    ! fcs=fcs/(ro_surf*Cp)
    ! ;On transforme w'theta' en w'thetav'
    ! Lv=(2.501-0.00237*(tcls-273.15))*1.E6
    ! xle=xle/(ro_surf*Lv)
    ! fcsv=fcs+.608*xle*tcls
    ! ***************************************************
    ! dif khfs est deja w't'_v / heatv(i) = khfs(i) + RETV*zxt*kqfs(i)
    ! AM calcul de Ro = paprs(i,1)/Rd zxt
    ! AM convention >0 vers le bas ds lmdz
    khfs(i) = -flux_t(i, 1)*zxt*rd/(rcpd*paprs(i,1))
    kqfs(i) = -flux_q(i, 1)*zxt*rd/(paprs(i,1))
    ! AM   verifier que khfs et kqfs sont bien de la forme w'l'
    heatv(i) = khfs(i) + retv*zxt*kqfs(i)
    ! a comparer aussi aux sorties de clqh : flux_T/RoCp et flux_q/RoLv
    ! AM ustar est en entree (calcul dans stdlevvar avec t2m q2m)
    ! Theta et qT du thermique sans exces
    qt_th(i) = q2m(i)
    ! Al1 Th_th restera la Theta du thermique sans exces jusqu'au 3eme calcul
    th_th(i) = t2m(i)
  END DO

  DO i = 1, knon
    rhino(i, 1) = 0.0 ! Global Richardson
    check(i) = .TRUE.
    pblh(i) = z(i, 1) ! on initialise pblh a l'altitude du 1er niveau
    ! Attention Plcl est pression ou altitude ?
    ! plcl(i) = 6000. ! m
    plcl(i) = 200. ! hPa
    IF (heatv(i)>0.0001) THEN
      ! Lambda = -u*^3 / (alpha.g.kvon.<w'Theta'v>
      obklen(i) = -t(i, 1)*ustar(i)**3/(rg*vk*heatv(i))
    ELSE
      ! set pblh to the friction high (cf + bas)
      pblh(i) = 700.0*ustar(i)
      check(i) = .FALSE.
    END IF
  END DO


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! PBL height calculation:
  ! Search for level of pbl. Scan upward until the Richardson number between
  ! the first level and the current level exceeds the "critical" value.
  ! (bonne idee Nu de separer le Ric et l'exces de temp du thermique)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fac = 100.0
  DO k = 2, isommet
    DO i = 1, knon
      IF (check(i)) THEN
        zdu2 = u(i, k)**2 + v(i, k)**2
        zdu2 = max(zdu2, 1.0E-20)
        ! Theta_v environnement
        zthvd = t(i, k)/s(i, k)*(1.+retv*q(i,k))
        zthvu = th_th(i)*(1.+retv*qt_th(i))
        ! Le Ri bulk par Theta_v
        rhino(i, k) = (z(i,k)-zref)*rg*(zthvd-zthvu)/(zdu2*0.5*(zthvd+zthvu))

        IF (rhino(i,k)>=ricr) THEN
          pblh(i) = z(i, k-1) + (z(i,k-1)-z(i,k))*(ricr-rhino(i,k-1))/(rhino(i,k-1)-rhino(i,k))
          ! test04 (la pblh est encore ici sous-estime'e)
          pblh(i) = pblh(i) + 100.
          ! pblT(i) = t(i,k-1) + (t(i,k)-t(i,k-1)) *
          ! .              (pblh(i)-z(i,k-1))/(z(i,k)-z(i,k-1))
          check(i) = .FALSE.
        END IF
      END IF
    END DO
  END DO


  ! Set pbl height to maximum value where computation exceeds number of
  ! layers allowed

  DO i = 1, knon
    IF (check(i)) pblh(i) = z(i, isommet)
  END DO

  ! Improve estimate of pbl height for the unstable points.
  ! Find unstable points (sensible heat flux is upward):

  DO i = 1, knon
    IF (heatv(i)>0.) THEN
      unstbl(i) = .TRUE.
      check(i) = .TRUE.
    ELSE
      unstbl(i) = .FALSE.
      check(i) = .FALSE.
    END IF
  END DO

  ! For the unstable case, compute velocity scale and the
  ! convective temperature excess:

  DO i = 1, knon
    IF (check(i)) THEN
      phiminv(i) = (1.-binm*pblh(i)/obklen(i))**onet
      ! ***************************************************
      ! Wm ? et W* ? c'est la formule pour z/h < .1
      ! ;Calcul de w* ;;
      ! ;;;;;;;;;;;;;;;;
      ! w_star=((g/tcls)*fcsv*z(ind))^(1/3.) [ou prendre la premiere approx de h)
      ! ;; CALCUL DE wm ;;
      ! ;;;;;;;;;;;;;;;;;;
      ! ; Ici on considerera que l'on est dans la couche de surf jusqu'a 100m
      ! ; On prend svt couche de surface=0.1*h mais on ne connait pas h
      ! ;;;;;;;;;;;Dans la couche de surface
      ! if (z(ind) le 20) then begin
      ! Phim=(1.-15.*(z(ind)/L))^(-1/3.)
      ! wm=u_star/Phim
      ! ;;;;;;;;;;;En dehors de la couche de surface
      ! endif else if (z(ind) gt 20) then begin
      ! wm=(u_star^3+c1*w_star^3)^(1/3.)
      ! endif
      ! ***************************************************
      wm(i) = ustar(i)*phiminv(i)
      ! ======================================================================
      ! valeurs de Dominique Lambert de la campagne SEMAPHORE :
      ! <T'^2> = 100.T*^2; <q'^2> = 20.q*^2 a 10m
      ! <Tv'^2> = (1+1.2q).100.T* + 1.2Tv.sqrt(20*100).T*.q* + (.608*Tv)^2*20.q*^2;
      ! et dTetavS = sqrt(<Tv'^2>) ainsi calculee.
      ! avec : T*=<w'T'>_s/w* et q*=<w'q'>/w*
      ! !!! on peut donc utiliser w* pour les fluctuations <-> Lambert
      ! (leur corellation pourrait dependre de beta par ex)
      ! if fcsv(i,j) gt 0 then begin
      ! dTetavs=b1*(1.+2.*.608*q_10(i,j))*(fcs(i,j)/wm(i,j))^2+$
      ! (.608*Thetav_10(i,j))^2*b2*(xle(i,j)/wm(i,j))^2+$
      ! 2.*.608*thetav_10(i,j)*sqrt(b1*b2)*(xle(i,j)/wm(i,j))*(fcs(i,j)/wm(i,j))
      ! dqs=b2*(xle(i,j)/wm(i,j))^2
      ! theta_s(i,j)=thetav_10(i,j)+sqrt(dTetavs)
      ! q_s(i,j)=q_10(i,j)+sqrt(dqs)
      ! endif else begin
      ! Theta_s(i,j)=thetav_10(i,j)
      ! q_s(i,j)=q_10(i,j)
      ! endelse
      ! leur reference est le niveau a 10m, mais on prend 2m ici.
      ! ======================================================================
      ! Premier calcul de l'exces tu thermique
      ! ======================================================================
      ! HBTM        therm(i) = heatv(i)*fak/wm(i)
      ! forme Mathieu :
      q_star = max(0., kqfs(i)/wm(i))
      t_star = max(0., khfs(i)/wm(i))
      ! Al1 Houston, we have a problem : il arrive en effet que heatv soit
      ! positif (=thermique instable) mais pas t_star : avec evaporation
      ! importante, il se peut qu'on refroidisse la 2m Que faire alors ?
      ! Garder le seul terme en q_star^2 ? ou rendre negatif le t_star^2 ?
      therm(i) = sqrt(b1*(1.+2.*retv*qt_th(i))*t_star**2+(retv*th_th(i))**2*b2*q_star*q_star+2.*retv*th_th(i)*b212* &
        q_star*t_star)

      ! Theta et qT du thermique (forme H&B) avec exces
      ! (attention, on ajoute therm(i) qui est virtuelle ...)
      ! pourquoi pas sqrt(b1)*t_star ?
      ! dqs = b2sr*kqfs(i)/wm(i)
      qt_th(i) = qt_th(i) + b2sr*q_star
      rhino(i, 1) = 0.0
    END IF
  END DO

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++ Improve pblh estimate for unstable conditions using the +++++++
  ! ++          convective temperature excess :                +++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  DO k = 2, isommet
    DO i = 1, knon
      IF (check(i)) THEN
        ! test     zdu2 = (u(i,k)-u(i,1))**2+(v(i,k)-v(i,1))**2+fac*ustar(i)**2
        zdu2 = u(i, k)**2 + v(i, k)**2
        zdu2 = max(zdu2, 1.0E-20)
        ! Theta_v environnement
        zthvd = t(i, k)/s(i, k)*(1.+retv*q(i,k))

        ! et therm Theta_v (avec hypothese de constance de H&B,
        ! qui assimile qT a vapeur)
        zthvu = th_th(i)*(1.+retv*qt_th(i)) + therm(i)


        ! Le Ri par Theta_v
        ! AM Niveau de ref 2m
        rhino(i, k) = (z(i,k)-zref)*rg*(zthvd-zthvu)/(zdu2*0.5*(zthvd+zthvu))

        ! Niveau critique atteint
        IF (rhino(i,k)>=ricr) THEN
          pblh(i) = z(i, k-1) + (z(i,k-1)-z(i,k))*(ricr-rhino(i,k-1))/(rhino(i,k-1)-rhino(i,k))
          ! test04
          pblh(i) = pblh(i) + 100.
          ! pblT(i) = t(i,k-1) + (t(i,k)-t(i,k-1)) *
          ! .              (pblh(i)-z(i,k-1))/(z(i,k)-z(i,k-1))
          check(i) = .FALSE.
        END IF
      END IF
    END DO
  END DO

  ! Set pbl height to maximum value where computation exceeds number of
  ! layers allowed (H&B)

  DO i = 1, knon
    IF (check(i)) pblh(i) = z(i, isommet)
  END DO

  ! PBL height must be greater than some minimum mechanical mixing depth
  ! Several investigators have proposed minimum mechanical mixing depth
  ! relationships as a function of the local friction velocity, u*.  We
  ! make use of a linear relationship of the form h = c u* where c=700.
  ! The scaling arguments that give rise to this relationship most often
  ! represent the coefficient c as some constant over the local coriolis
  ! parameter.  Here we make use of the experimental results of Koracin
  ! and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
  ! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
  ! latitude value for f so that c = 0.07/f = 700.  (H&B)
  ! Al1 calcul de pblT dans ce cas
  DO i = 1, knon
    pblmin = 700.0*ustar(i)
    IF (pblh(i)<pblmin) check(i) = .TRUE.
  END DO
  DO i = 1, knon
    IF (check(i)) THEN
      pblh(i) = 700.0*ustar(i)
      ! et par exemple :
      ! pblT(i) = t(i,2) + (t(i,3)-t(i,2)) *
      ! .              (pblh(i)-z(i,2))/(z(i,3)-z(i,2))
    END IF
  END DO

  ! ********************************************************************
  ! pblh is now available; do preparation for final calculations :
  ! ********************************************************************
  DO i = 1, knon
    check(i) = .TRUE.
    zsat(i) = .FALSE.
    zcin(i) = .FALSE.
    ! omegafl utilise pour prolongement CAPE
    omegafl(i) = .FALSE.

    ! Do additional preparation for unstable cases only, set temperature
    ! and moisture perturbations depending on stability.
    ! Rq: les formules sont prises dans leur forme Couche de Surface
    IF (unstbl(i)) THEN
      ! Al pblh a change', on recalcule :
      zxt = (th_th(i)-zref*0.5*rg/rcpd/(1.+rvtmp2*qt_th(i)))*(1.+retv*qt_th(i))
      phiminv(i) = (1.-binm*pblh(i)/obklen(i))**onet
      phihinv(i) = sqrt(1.-binh*pblh(i)/obklen(i))
      wm(i) = ustar(i)*phiminv(i)
    END IF
  END DO


  ! =======================================================
  ! last upward integration
  ! For all unstable layers, compute integral info and CTEI
  ! =======================================================

  ! 1/Recompute surface characteristics with the improved pblh
  ! ----------------------------------------------------------
  DO i = 1, knon
    IF (unstbl(i)) THEN
      diagok(i) = 1.
      ! from missing_value to zero
      cape(i) = 0.
      cin(i) = 0.
      eauliq(i) = 0.
      ctei(i) = 0.
      d_qt(i) = 0.
      d_thv(i) = 0.
      dlt_2(i) = 0.
      xhis(i) = 0.
      posint(i) = 0.
      kin(i) = 0.
      omega(i) = 0.

      phiminv(i) = (1.-binm*pblh(i)/obklen(i))**onet
      wm(i) = ustar(i)*phiminv(i)
      q_star = max(0., kqfs(i)/wm(i))
      t_star = max(0., khfs(i)/wm(i))
      therm(i) = sqrt(b1*(1.+2.*retv*qt_th(i))*t_star**2+(retv*th_th(i))**2*b2*q_star*q_star+2.*retv*th_th(i)*b212* &
        q_star*t_star)
      ! Al1diag
      ! trmb1(i) = b1*(1.+2.*RETV*qT_th(i))*t_star**2
      ! trmb2(i) = (RETV*Th_th(i))**2*b2*q_star*q_star
      ! trmb3(i) = 2.*RETV*Th_th(i)*b212*q_star*t_star

      ! Th_th will now be the thermal-theta (including exces)
      ! c         Th_th(i) = Th_th(i)+sqrt(b1)*max(0.,khfs(i)/wm(i))
      th_th(i) = th_th(i) + therm(i)
      ! al1diag
      ! trmb2(i) = wm(i)
      ! trmb3(i) = phiminv(i)
      ! and computes Theta_e for thermal
      the_th(i) = th_th(i) + rlvcp*qt_th(i)
    END IF ! unstbl
    ! Al1 compute a first guess of Plcl with the Bolton/Emanuel formula
    t2 = th_th(i)
    ! thermodyn functions
    zdelta = max(0., sign(1.,rtt-t2))
    qsat = r2es*foeew(t2, zdelta)/paprs(i, 1)
    qsat = min(0.5, qsat)
    zcor = 1./(1.-retv*qsat)
    qsat = qsat*zcor
    ! relative humidity of thermal at 2m
    rh = qt_th(i)/qsat
    chi = t2/(1669.0-122.0*rh-t2)
    plcl(i) = paprs(i, 1)*(rh**chi)
    ! al1diag
    ! ctei(i) = Plcl(i)
    ! cape(i) = T2
    ! trmb1(i)= Chi
    ! select unstable columns (=thermals)
    check(i) = .FALSE.
    IF (heatv(i)>0.) check(i) = .TRUE.
    ! diag
    ! dTv21(i,1) = T2*(1+RETV*qT_th(i))-t(i,1)*(1+RETV*q(i,1))
  END DO
  ! ----------------------------------
  ! 2/ upward integration for thermals
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ k loop
  DO k = 2, isommet
    DO i = 1, knon
      IF (check(i) .OR. omegafl(i)) THEN
        ! CC         if (pplay(i,k) .le. plcl(i)) then
        zm(i) = z(i, k-1)
        zp(i) = z(i, k)
        ! Environnement : calcul de Tv1 a partir de t(:,:)== T liquide
        ! ==============
        tl1 = t(i, k)
        t1 = tl1
        zdelta = max(0., sign(1.,rtt-t1))
        qsat = r2es*foeew(t1, zdelta)/pplay(i, k)
        qsat = min(0.5, qsat)
        zcor = 1./(1.-retv*qsat)
        qsat = qsat*zcor
        q1 = min(q(i,k), qsat)
        ql1 = max(0., q(i,k)-q1)
        ! thermodyn function (Tl2Tql)
        dt = rlvcp*ql1
        DO WHILE (abs(dt)>=dt0)
          t1 = t1 + dt
          zdelta = max(0., sign(1.,rtt-t1))
          zcvm5 = r5les*(1.-zdelta) + r5ies*zdelta
          qsat = r2es*foeew(t1, zdelta)/pplay(i, k)
          qsat = min(0.5, qsat)
          zcor = 1./(1.-retv*qsat)
          qsat = qsat*zcor
          dqsat_dt = foede(t1, zdelta, zcvm5, qsat, zcor)
          ! correction lineaire pour conserver Tl env
          ! << Tl = T1 + DT - RLvCp*(ql1 - dqsat/dT*DT >>
          denom = 1. + rlvcp*dqsat_dt
          q1 = min(q(i,k), qsat)
          ql1 = q(i, k) - q1 ! can be negative
          rnum = tl1 - t1 + rlvcp*ql1
          dt = rnum/denom
        END DO
        ql1 = max(0., ql1)
        tv1 = t1*(1.+retv*q1-ql1)
        ! Thermique    : on atteint le seuil B/E de condensation
        ! ==============

        IF (.NOT. zsat(i)) THEN
          ! first guess from The_th(i) = Th_th(i) + RLvCp* [qv=qT_th(i)]
          t2 = s(i, k)*the_th(i) - rlvcp*qt_th(i)
          zdelta = max(0., sign(1.,rtt-t2))
          qsat = r2es*foeew(t2, zdelta)/pplay(i, k)
          qsat = min(0.5, qsat)
          zcor = 1./(1.-retv*qsat)
          qsat = qsat*zcor
          q2 = min(qt_th(i), qsat)
          ql2 = max(0., qt_th(i)-q2)
          IF (ql2>0.0001) zsat(i) = .TRUE.
          tbef(i) = t2
          ! a PBLH non sature
          IF (zm(i)<pblh(i) .AND. zp(i)>=pblh(i)) THEN
            reduc = (pblh(i)-zm(i))/(zp(i)-zm(i))
            spblh = s(i, k-1) + reduc*(s(i,k)-s(i,k-1))
            ! lmdz : qT1 et Thv1
            t1 = (t(i,k-1)+reduc*(t(i,k)-t(i,k-1)))
            thv1 = t1*(1.+retv*q(i,k))/spblh
            ! on calcule pour le cas sans nuage un ctei en Delta Thv
            thv2 = t2/spblh*(1.+retv*qt_th(i))
            ctei(i) = thv1 - thv2
            tv2 = t2*(1.+retv*q2-ql2)
            ! diag
            ! dTv21(i,k) = Tv2-Tv1
            check(i) = .FALSE.
            omegafl(i) = .TRUE.
          END IF
        END IF

        IF (zsat(i)) THEN
          ! thermodyn functions (Te2Tqsat)
          t2 = tbef(i)
          dt = 1.
          te2 = s(i, k)*the_th(i)
          DO WHILE (abs(dt)>=dt0)
            zdelta = max(0., sign(1.,rtt-t2))
            zcvm5 = r5les*(1.-zdelta) + r5ies*zdelta
            qsat = r2es*foeew(t2, zdelta)/pplay(i, k)
            qsat = min(0.5, qsat)
            zcor = 1./(1.-retv*qsat)
            qsat = qsat*zcor
            dqsat_dt = foede(t2, zdelta, zcvm5, qsat, zcor)
            ! correction lineaire pour conserver Te_th
            ! << Te = T2 + DT + RLvCp*(qsatbef + dq/dT*DT >>
            denom = 1. + rlvcp*dqsat_dt
            rnum = te2 - t2 - rlvcp*qsat
            dt = rnum/denom
            t2 = t2 + dt
          END DO
          q2 = min(qt_th(i), qsat)
          ql2 = max(0., qt_th(i)-q2)
          ! jusqu'a PBLH y compris
          IF (zm(i)<pblh(i)) THEN

            ! mais a PBLH, interpolation et complements
            IF (zp(i)>=pblh(i)) THEN
              reduc = (pblh(i)-zm(i))/(zp(i)-zm(i))
              spblh = s(i, k-1) + reduc*(s(i,k)-s(i,k-1))
              ! CAPE et EauLiq a pblH
              cape(i) = kape(i) + reduc*(zp(i)-zm(i))*rg*.5/(tv2+tv1)*max(0., (tv2-tv1))
              eauliq(i) = eauliq(i) + reduc*(paprs(i,k-1)-paprs(i,k))*ql2/rg
              ! CTEI
              the2 = (t2+rlvcp*q2)/spblh
              ! T1 est en realite la Tl env (on a donc strict The1)
              t1 = (t(i,k-1)+reduc*(t(i,k)-t(i,k-1)))
              the1 = (t1+rlvcp*q(i,k))/spblh
              ! Calcul de la Cloud Top Entrainement Instability
              ! cf Mathieu Lahellec QJRMS (2005) Comments to DYCOMS-II
              ! saut a l'inversion :
              delt_the = the1 - the2 ! negatif
              delt_qt = q(i, k) - qt_th(i) ! negatif
              d_qt(i) = -delt_qt
              dlt_2(i) = .63*delt_the - the2*delt_qt
              ! init ctei(i)
              ctei(i) = dlt_2(i)
              IF (dlt_2(i)<-0.1) THEN
                ! integrale de Peter :
                aa = delt_the - delt_qt*(rlvcp-retv*the2)
                bb = (rlvcp-(1.+retv)*the2)*ql2
                d_thv(i) = aa - bb
                ! approx de Xhi_s et de l'integrale Xint=ctei(i)
                xhis(i) = bb/(aa-dlt_2(i))
                ! trmb1(i) = xhis
                ! trmb3(i) = dlt_2
                xnull = bb/aa
                IF (xhis(i)>0.1) THEN
                  ctei(i) = dlt_2(i)*xhis(i) + aa*(1.-xhis(i)) + bb*alog(xhis(i))
                ELSE
                  ctei(i) = .5*(dlt_2(i)+aa-bb)
                END IF
                IF (xnull>0.) THEN
                  posint(i) = aa - bb + bb*alog(xnull)
                ELSE
                  posint(i) = 0.
                END IF
              ELSE
                ctei(i) = 1.
                posint(i) = 1.
              END IF
              check(i) = .FALSE.
              omegafl(i) = .TRUE.
            END IF ! end a pblh
            IF (check(i)) eauliq(i) = eauliq(i) + (paprs(i,k)-paprs(i,k+1))*ql2/rg
          END IF

        END IF ! Zsat

        ! KAPE : thermique / environnement
        tv2 = t2*(1.+retv*q2-ql2)
        ! diag
        ! dTv21(i,k) = Tv2-Tv1
        ! Kape courante
        kape(i) = kape(i) + (zp(i)-zm(i))*rg*.5/(tv2+tv1)*max(0., (tv2-tv1))
        ! Cin
        IF (zcin(i) .AND. tv2-tv1>0.) THEN
          zcin(i) = .FALSE.
          cin(i) = kin(i)
        END IF
        IF (.NOT. zcin(i) .AND. tv2-tv1<0.) THEN
          zcin(i) = .TRUE.
          kin(i) = kin(i) + (zp(i)-zm(i))*rg*.5/(tv2+tv1)*min(0., (tv2-tv1))
        END IF
        IF (kape(i)+kin(i)<0.) THEN
          omega(i) = zm(i)
          ! trmb3(i) = paprs(i,k)
          omegafl(i) = .FALSE.
          ! diag
          ! print*,'Tv2-Tv1 (k): ',i,(dTv21(i,j),j=1,k)
        END IF
        ! CC         EndIf !plcl
      END IF ! check(i)
    END DO
  END DO ! end of level loop
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ end k loop
  RETURN
END SUBROUTINE hbtm2l
