
! $Header$


SUBROUTINE hbtm(knon, paprs, pplay, t2m, t10m, q2m, q10m, ustar, wstar, &
    flux_t, flux_q, u, v, t, q, pblh, cape, eauliq, ctei, pblt, therm, trmb1, &
    trmb2, trmb3, plcl)
  USE dimphy
  IMPLICIT NONE

  ! ***************************************************************
  ! *                                                             *
  ! * HBTM2   D'apres Holstag&Boville et Troen&Mahrt              *
  ! *                 JAS 47              BLM                     *
  ! * Algorithme These Anne Mathieu                               *
  ! * Critere d'Entrainement Peter Duynkerke (JAS 50)             *
  ! * written by  : Anne MATHIEU & Alain LAHELLEC, 22/11/99       *
  ! * features : implem. exces Mathieu                            *
  ! ***************************************************************
  ! * mods : decembre 99 passage th a niveau plus bas. voir fixer *
  ! * la prise du th a z/Lambda = -.2 (max Ray)                   *
  ! * Autre algo : entrainement ~ Theta+v =cste mais comment=>The?*
  ! * on peut fixer q a .7qsat(cf non adiab)=>T2 et The2          *
  ! * voir aussi //KE pblh = niveau The_e ou l = env.             *
  ! ***************************************************************
  ! * fin therm a la HBTM passage a forme Mathieu 12/09/2001      *
  ! ***************************************************************
  ! *


  ! AM Fev 2003
  ! Adaptation a LMDZ version couplee

  ! Pour le moment on fait passer en argument les grdeurs de surface :
  ! flux, t,q2m, t,q10m, on va utiliser systematiquement les grdeurs a 2m ms
  ! on garde la possibilite de changer si besoin est (jusqu'a present la
  ! forme de HB avec le 1er niveau modele etait conservee)





  include "YOMCST.h"
  REAL rlvcp, reps
  ! Arguments:

  INTEGER knon ! nombre de points a calculer
  ! AM
  REAL t2m(klon), t10m(klon) ! temperature a 2 et 10m
  REAL q2m(klon), q10m(klon) ! q a 2 et 10m
  REAL ustar(klon)
  REAL wstar(klon) ! w*, convective velocity scale
  REAL paprs(klon, klev+1) ! pression a inter-couche (Pa)
  REAL pplay(klon, klev) ! pression au milieu de couche (Pa)
  REAL flux_t(klon, klev), flux_q(klon, klev) ! Flux
  REAL u(klon, klev) ! vitesse U (m/s)
  REAL v(klon, klev) ! vitesse V (m/s)
  REAL t(klon, klev) ! temperature (K)
  REAL q(klon, klev) ! vapeur d'eau (kg/kg)
  ! AM      REAL cd_h(klon) ! coefficient de friction au sol pour chaleur
  ! AM      REAL cd_m(klon) ! coefficient de friction au sol pour vitesse

  INTEGER isommet
  ! um      PARAMETER (isommet=klev) ! limite max sommet pbl
  REAL, PARAMETER :: vk = 0.35 ! Von Karman => passer a .41 ! cf U.Olgstrom
  REAL, PARAMETER :: ricr = 0.4
  REAL, PARAMETER :: fak = 8.5 ! b calcul du Prandtl et de dTetas
  REAL, PARAMETER :: fakn = 7.2 ! a
  REAL, PARAMETER :: onet = 1.0/3.0
  REAL, PARAMETER :: t_coup = 273.15
  REAL, PARAMETER :: zkmin = 0.01
  REAL, PARAMETER :: betam = 15.0 ! pour Phim / h dans la S.L stable
  REAL, PARAMETER :: betah = 15.0
  REAL, PARAMETER :: betas = 5.0 ! Phit dans la S.L. stable (mais 2 formes / z/OBL<>1
  REAL, PARAMETER :: sffrac = 0.1 ! S.L. = z/h < .1
  REAL, PARAMETER :: usmin = 1.E-12
  REAL, PARAMETER :: binm = betam*sffrac
  REAL, PARAMETER :: binh = betah*sffrac
  REAL, PARAMETER :: ccon = fak*sffrac*vk
  REAL, PARAMETER :: b1 = 70., b2 = 20.
  REAL, PARAMETER :: zref = 2. ! Niveau de ref a 2m peut eventuellement
  ! etre choisi a 10m
  REAL q_star, t_star
  REAL b212, b2sr ! Lambert correlations T' q' avec T* q*

  REAL z(klon, klev)
  ! AM      REAL pcfm(klon,klev), pcfh(klon,klev)
  INTEGER i, k, j
  REAL zxt
  ! AM      REAL zxt, zxq, zxu, zxv, zxmod, taux, tauy
  ! AM      REAL zx_alf1, zx_alf2 ! parametres pour extrapolation
  REAL khfs(klon) ! surface kinematic heat flux [mK/s]
  REAL kqfs(klon) ! sfc kinematic constituent flux [m/s]
  REAL heatv(klon) ! surface virtual heat flux
  REAL rhino(klon, klev) ! bulk Richardon no. mais en Theta_v
  LOGICAL unstbl(klon) ! pts w/unstbl pbl (positive virtual ht flx)
  LOGICAL stblev(klon) ! stable pbl with levels within pbl
  LOGICAL unslev(klon) ! unstbl pbl with levels within pbl
  LOGICAL unssrf(klon) ! unstb pbl w/lvls within srf pbl lyr
  LOGICAL unsout(klon) ! unstb pbl w/lvls in outer pbl lyr
  LOGICAL check(klon) ! True=>chk if Richardson no.>critcal
  LOGICAL omegafl(klon) ! flag de prolongerment cape pour pt Omega
  REAL pblh(klon)
  REAL pblt(klon)
  REAL plcl(klon)
  ! AM      REAL cgh(klon,2:klev) ! counter-gradient term for heat [K/m]
  ! AM      REAL cgq(klon,2:klev) ! counter-gradient term for constituents
  ! AM      REAL cgs(klon,2:klev) ! counter-gradient star (cg/flux)
  REAL unsobklen(klon) ! Monin-Obukhov lengh
  ! AM      REAL ztvd, ztvu,
  REAL zdu2
  REAL therm(klon) ! thermal virtual temperature excess
  REAL trmb1(klon), trmb2(klon), trmb3(klon)
  ! Algorithme thermique
  REAL s(klon, klev) ! [P/Po]^Kappa milieux couches
  REAL th_th(klon) ! potential temperature of thermal
  REAL the_th(klon) ! equivalent potential temperature of thermal
  REAL qt_th(klon) ! total water  of thermal
  REAL tbef(klon) ! T thermique niveau precedent
  REAL qsatbef(klon)
  LOGICAL zsat(klon) ! le thermique est sature
  REAL cape(klon) ! Cape du thermique
  REAL kape(klon) ! Cape locale
  REAL eauliq(klon) ! Eau liqu integr du thermique
  REAL ctei(klon) ! Critere d'instab d'entrainmt des nuages de CL
  REAL the1, the2, aa, bb, zthvd, zthvu, xintpos, qqsat
  ! IM 091204 BEG
  REAL a1, a2, a3
  ! IM 091204 END
  REAL xhis, rnum, denom, th1, th2, thv1, thv2, ql2
  REAL dqsat_dt, qsat2, qt1, q2, t1, t2, xnull, delt_the
  REAL delt_qt, delt_2, quadsat, spblh, reduc

  REAL phiminv(klon) ! inverse phi function for momentum
  REAL phihinv(klon) ! inverse phi function for heat
  REAL wm(klon) ! turbulent velocity scale for momentum
  REAL fak1(klon) ! k*ustar*pblh
  REAL fak2(klon) ! k*wm*pblh
  REAL fak3(klon) ! fakn*wstar/wm
  REAL pblk(klon) ! level eddy diffusivity for momentum
  REAL pr(klon) ! Prandtl number for eddy diffusivities
  REAL zl(klon) ! zmzp / Obukhov length
  REAL zh(klon) ! zmzp / pblh
  REAL zzh(klon) ! (1-(zmzp/pblh))**2
  REAL zm(klon) ! current level height
  REAL zp(klon) ! current level height + one level up
  REAL zcor, zdelta, zcvm5
  ! AM      REAL zxqs
  REAL fac, pblmin, zmzp, term

  include "YOETHF.h"
  include "FCTTRE.h"



  ! initialisations (Anne)
  isommet = klev
  th_th(:) = 0.
  q_star = 0
  t_star = 0


  b212 = sqrt(b1*b2)
  b2sr = sqrt(b2)

  ! ============================================================
  ! Fonctions thermo implicites
  ! ============================================================
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Tetens : pression partielle de vap d'eau e_sat(T)
  ! =================================================
  ! ++ e_sat(T) = r2*exp( r3*(T-Tf)/(T-r4) ) id a r2*FOEWE
  ! ++ avec :
  ! ++ Tf = 273.16 K  (Temp de fusion de la glace)
  ! ++ r2 = 611.14 Pa
  ! ++ r3 = 17.269 (liquide) 21.875 (solide) adim
  ! ++ r4 = 35.86             7.66           Kelvin
  ! ++  q_sat = eps*e_sat/(p-(1-eps)*e_sat)
  ! ++ deriv� :
  ! ++ =========
  ! ++                   r3*(Tf-r4)*q_sat(T,p)
  ! ++ d_qsat_dT = --------------------------------
  ! ++             (T-r4)^2*( 1-(1-eps)*e_sat(T)/p )
  ! ++ pour zcvm5=Lv, c'est FOEDE
  ! ++ Rq :(1.-REPS)*esarg/Parg id a RETV*Qsat
  ! ------------------------------------------------------------------

  ! Initialisation
  rlvcp = rlvtt/rcpd
  reps = rd/rv


  ! DO i = 1, klon
  ! pcfh(i,1) = cd_h(i)
  ! pcfm(i,1) = cd_m(i)
  ! ENDDO
  ! DO k = 2, klev
  ! DO i = 1, klon
  ! pcfh(i,k) = zkmin
  ! pcfm(i,k) = zkmin
  ! cgs(i,k) = 0.0
  ! cgh(i,k) = 0.0
  ! cgq(i,k) = 0.0
  ! ENDDO
  ! ENDDO

  ! Calculer les hauteurs de chaque couche
  ! (geopotentielle Int_dp/ro = Int_[Rd.T.dp/p] z = geop/g)
  ! pourquoi ne pas utiliser Phi/RG ?
  DO i = 1, knon
    z(i, 1) = rd*t(i, 1)/(0.5*(paprs(i,1)+pplay(i,1)))*(paprs(i,1)-pplay(i,1) &
      )/rg
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
      z(i, k) = z(i, k-1) + rd*0.5*(t(i,k-1)+t(i,k))/paprs(i, k)*(pplay(i,k-1 &
        )-pplay(i,k))/rg
      s(i, k) = (pplay(i,k)/paprs(i,1))**rkappa
    END DO
  END DO
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++  Determination des grandeurs de surface  +++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO i = 1, knon
    ! AM         IF (thermcep) THEN
    ! AM           zdelta=MAX(0.,SIGN(1.,RTT-tsol(i)))
    ! zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
    ! zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*q(i,1))
    ! AM           zxqs= r2es * FOEEW(tsol(i),zdelta)/paprs(i,1)
    ! AM           zxqs=MIN(0.5,zxqs)
    ! AM           zcor=1./(1.-retv*zxqs)
    ! AM           zxqs=zxqs*zcor
    ! AM         ELSE
    ! AM           IF (tsol(i).LT.t_coup) THEN
    ! AM              zxqs = qsats(tsol(i)) / paprs(i,1)
    ! AM           ELSE
    ! AM              zxqs = qsatl(tsol(i)) / paprs(i,1)
    ! AM           ENDIF
    ! AM         ENDIF
    ! niveau de reference bulk; mais ici, c,a pourrait etre le niveau de ref
    ! du thermique
    ! AM        zx_alf1 = 1.0
    ! AM        zx_alf2 = 1.0 - zx_alf1
    ! AM        zxt = (t(i,1)+z(i,1)*RG/RCPD/(1.+RVTMP2*q(i,1)))
    ! AM     .        *(1.+RETV*q(i,1))*zx_alf1
    ! AM     .      + (t(i,2)+z(i,2)*RG/RCPD/(1.+RVTMP2*q(i,2)))
    ! AM     .        *(1.+RETV*q(i,2))*zx_alf2
    ! AM        zxu = u(i,1)*zx_alf1+u(i,2)*zx_alf2
    ! AM        zxv = v(i,1)*zx_alf1+v(i,2)*zx_alf2
    ! AM        zxq = q(i,1)*zx_alf1+q(i,2)*zx_alf2
    ! AM
    ! AMAM           zxu = u10m(i)
    ! AMAM           zxv = v10m(i)
    ! AMAM           zxmod = 1.0+SQRT(zxu**2+zxv**2)
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
    ! AM        khfs(i) = (tsol(i)*(1.+RETV*q(i,1))-zxt) *zxmod*cd_h(i)
    ! AM        kqfs(i) = (zxqs-zxq) *zxmod*cd_h(i) * beta(i)
    ! AM
    ! dif khfs est deja w't'_v / heatv(i) = khfs(i) + RETV*zxt*kqfs(i)
    ! AM calcule de Ro = paprs(i,1)/Rd zxt
    ! AM convention >0 vers le bas ds lmdz
    khfs(i) = -flux_t(i, 1)*zxt*rd/(rcpd*paprs(i,1))
    kqfs(i) = -flux_q(i, 1)*zxt*rd/(paprs(i,1))
    ! AM   verifier que khfs et kqfs sont bien de la forme w'l'
    heatv(i) = khfs(i) + 0.608*zxt*kqfs(i)
    ! a comparer aussi aux sorties de clqh : flux_T/RoCp et flux_q/RoLv
    ! AM        heatv(i) = khfs(i)
    ! AM ustar est en entree
    ! AM        taux = zxu *zxmod*cd_m(i)
    ! AM        tauy = zxv *zxmod*cd_m(i)
    ! AM        ustar(i) = SQRT(taux**2+tauy**2)
    ! AM        ustar(i) = MAX(SQRT(ustar(i)),0.01)
    ! Theta et qT du thermique sans exces (interpolin vers surf)
    ! chgt de niveau du thermique (jeudi 30/12/1999)
    ! (interpolation lineaire avant integration phi_h)
    ! AM        qT_th(i) = zxqs*beta(i) + 4./z(i,1)*(q(i,1)-zxqs*beta(i))
    ! AM        qT_th(i) = max(qT_th(i),q(i,1))
    qt_th(i) = q2m(i)
    ! n The_th restera la Theta du thermique sans exces jusqu'a 2eme calcul
    ! n reste a regler convention P) pour Theta
    ! The_th(i) = tsol(i) + 4./z(i,1)*(t(i,1)-tsol(i))
    ! -                      + RLvCp*qT_th(i)
    ! AM        Th_th(i) = tsol(i) + 4./z(i,1)*(t(i,1)-tsol(i))
    th_th(i) = t2m(i)
  END DO

  DO i = 1, knon
    rhino(i, 1) = 0.0 ! Global Richardson
    check(i) = .TRUE.
    pblh(i) = z(i, 1) ! on initialise pblh a l'altitude du 1er niveau
    plcl(i) = 6000.
    ! Lambda = -u*^3 / (alpha.g.kvon.<w'Theta'v>
    unsobklen(i) = -rg*vk*heatv(i)/(t(i,1)*max(ustar(i),usmin)**3)
    trmb1(i) = 0.
    trmb2(i) = 0.
    trmb3(i) = 0.
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
        ! pourquoi / niveau 1 (au lieu du sol) et le terme en u*^2 ?
        ! test     zdu2 =
        ! (u(i,k)-u(i,1))**2+(v(i,k)-v(i,1))**2+fac*ustar(i)**2
        zdu2 = u(i, k)**2 + v(i, k)**2
        zdu2 = max(zdu2, 1.0E-20)
        ! Theta_v environnement
        zthvd = t(i, k)/s(i, k)*(1.+retv*q(i,k))

        ! therm Theta_v sans exces (avec hypothese fausse de H&B, sinon,
        ! passer par Theta_e et virpot)
        ! zthvu=t(i,1)/s(i,1)*(1.+RETV*q(i,1))
        ! AM         zthvu = Th_th(i)*(1.+RETV*q(i,1))
        zthvu = th_th(i)*(1.+retv*qt_th(i))
        ! Le Ri par Theta_v
        ! AM         rhino(i,k) = (z(i,k)-z(i,1))*RG*(zthvd-zthvu)
        ! AM     .               /(zdu2*0.5*(zthvd+zthvu))
        ! AM On a nveau de ref a 2m ???
        rhino(i, k) = (z(i,k)-zref)*rg*(zthvd-zthvu)/(zdu2*0.5*(zthvd+zthvu))

        IF (rhino(i,k)>=ricr) THEN
          pblh(i) = z(i, k-1) + (z(i,k-1)-z(i,k))*(ricr-rhino(i,k-1))/(rhino( &
            i,k-1)-rhino(i,k))
          ! test04
          pblh(i) = pblh(i) + 100.
          pblt(i) = t(i, k-1) + (t(i,k)-t(i,k-1))*(pblh(i)-z(i,k-1))/(z(i,k)- &
            z(i,k-1))
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
      phiminv(i) = (1.-binm*pblh(i)*unsobklen(i))**onet
      ! ***************************************************
      ! Wm ? et W* ? c'est la formule pour z/h < .1
      ! ;Calcul de w* ;;
      ! ;;;;;;;;;;;;;;;;
      ! w_star=((g/tcls)*fcsv*z(ind))^(1/3.) [ou prendre la premiere approx
      ! de h)
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
      ! <Tv'^2> = (1+1.2q).100.T* + 1.2Tv.sqrt(20*100).T*.q* +
      ! (.608*Tv)^2*20.q*^2;
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
      ! ======================================================================

      ! HBTM        therm(i) = heatv(i)*fak/wm(i)
      ! forme Mathieu :
      q_star = kqfs(i)/wm(i)
      t_star = khfs(i)/wm(i)
      ! IM 091204 BEG
      IF (1==0) THEN
        IF (t_star<0. .OR. q_star<0.) THEN
          PRINT *, 'i t_star q_star khfs kqfs wm', i, t_star, q_star, &
            khfs(i), kqfs(i), wm(i)
        END IF
      END IF
      ! IM 091204 END
      ! AM Nveau cde ref 2m =>
      ! AM        therm(i) = sqrt( b1*(1.+2.*RETV*q(i,1))*t_star**2
      ! AM     +             + (RETV*T(i,1))**2*b2*q_star**2
      ! AM     +             + 2.*RETV*T(i,1)*b212*q_star*t_star
      ! AM     +                 )
      ! IM 091204 BEG
      a1 = b1*(1.+2.*retv*qt_th(i))*t_star**2
      a2 = (retv*th_th(i))**2*b2*q_star*q_star
      a3 = 2.*retv*th_th(i)*b212*q_star*t_star
      aa = a1 + a2 + a3
      IF (1==0) THEN
        IF (aa<0.) THEN
          PRINT *, 'i a1 a2 a3 aa', i, a1, a2, a3, aa
          PRINT *, 'i qT_th Th_th t_star q_star RETV b1 b2 b212', i, &
            qt_th(i), th_th(i), t_star, q_star, retv, b1, b2, b212
        END IF
      END IF
      ! IM 091204 END
      therm(i) = sqrt(b1*(1.+2.*retv*qt_th(i))*t_star**2+(retv*th_th( &
        i))**2*b2*q_star*q_star &  ! IM 101204  +             +
                                   ! 2.*RETV*Th_th(i)*b212*q_star*t_star
        +max(0.,2.*retv*th_th(i)*b212*q_star*t_star))

      ! Theta et qT du thermique (forme H&B) avec exces
      ! (attention, on ajoute therm(i) qui est virtuelle ...)
      ! pourquoi pas sqrt(b1)*t_star ?
      ! dqs = b2sr*kqfs(i)/wm(i)
      qt_th(i) = qt_th(i) + b2sr*q_star
      ! new on differre le calcul de Theta_e
      ! The_th(i) = The_th(i) + therm(i) + RLvCp*qT_th(i)
      ! ou:    The_th(i) = The_th(i) + sqrt(b1)*khfs(i)/wm(i) +
      ! RLvCp*qT_th(i)
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
        ! test     zdu2 =
        ! (u(i,k)-u(i,1))**2+(v(i,k)-v(i,1))**2+fac*ustar(i)**2
        zdu2 = u(i, k)**2 + v(i, k)**2
        zdu2 = max(zdu2, 1.0E-20)
        ! Theta_v environnement
        zthvd = t(i, k)/s(i, k)*(1.+retv*q(i,k))

        ! et therm Theta_v (avec hypothese de constance de H&B,
        ! zthvu=(t(i,1)+therm(i))/s(i,1)*(1.+RETV*q(i,1))
        zthvu = th_th(i)*(1.+retv*qt_th(i)) + therm(i)


        ! Le Ri par Theta_v
        ! AM Niveau de ref 2m
        ! AM         rhino(i,k) = (z(i,k)-z(i,1))*RG*(zthvd-zthvu)
        ! AM     .               /(zdu2*0.5*(zthvd+zthvu))
        rhino(i, k) = (z(i,k)-zref)*rg*(zthvd-zthvu)/(zdu2*0.5*(zthvd+zthvu))


        IF (rhino(i,k)>=ricr) THEN
          pblh(i) = z(i, k-1) + (z(i,k-1)-z(i,k))*(ricr-rhino(i,k-1))/(rhino( &
            i,k-1)-rhino(i,k))
          ! test04
          pblh(i) = pblh(i) + 100.
          pblt(i) = t(i, k-1) + (t(i,k)-t(i,k-1))*(pblh(i)-z(i,k-1))/(z(i,k)- &
            z(i,k-1))
          check(i) = .FALSE.
          ! IM 170305 BEG
          IF (1==0) THEN
            ! debug print -120;34       -34-        58 et    0;26 wamp
            IF (i==950 .OR. i==192 .OR. i==624 .OR. i==118) THEN
              PRINT *, ' i,Th_th,Therm,qT :', i, th_th(i), therm(i), qt_th(i)
              q_star = kqfs(i)/wm(i)
              t_star = khfs(i)/wm(i)
              PRINT *, 'q* t*, b1,b2,b212 ', q_star, t_star, &
                b1*(1.+2.*retv*qt_th(i))*t_star**2, &
                (retv*th_th(i))**2*b2*q_star**2, 2.*retv*th_th(i)*b212*q_star &
                *t_star
              PRINT *, 'zdu2 ,100.*ustar(i)**2', zdu2, fac*ustar(i)**2
            END IF
          END IF !(1.EQ.0) THEN
          ! IM 170305 END
          ! q_star = kqfs(i)/wm(i)
          ! t_star = khfs(i)/wm(i)
          ! trmb1(i) = b1*(1.+2.*RETV*q(i,1))*t_star**2
          ! trmb2(i) = (RETV*T(i,1))**2*b2*q_star**2
          ! Omega now   trmb3(i) = 2.*RETV*T(i,1)*b212*q_star*t_star
        END IF
      END IF
    END DO
  END DO

  ! Set pbl height to maximum value where computation exceeds number of
  ! layers allowed

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
  ! latitude value for f so that c = 0.07/f = 700.

  DO i = 1, knon
    pblmin = 700.0*ustar(i)
    pblh(i) = max(pblh(i), pblmin)
    ! par exemple :
    pblt(i) = t(i, 2) + (t(i,3)-t(i,2))*(pblh(i)-z(i,2))/(z(i,3)-z(i,2))
  END DO

  ! ********************************************************************
  ! pblh is now available; do preparation for diffusivity calculation :
  ! ********************************************************************
  DO i = 1, knon
    check(i) = .TRUE.
    zsat(i) = .FALSE.
    ! omegafl utilise pour prolongement CAPE
    omegafl(i) = .FALSE.
    cape(i) = 0.
    kape(i) = 0.
    eauliq(i) = 0.
    ctei(i) = 0.
    pblk(i) = 0.0
    fak1(i) = ustar(i)*pblh(i)*vk

    ! Do additional preparation for unstable cases only, set temperature
    ! and moisture perturbations depending on stability.
    ! *** Rq: les formule sont prises dans leur forme CS ***
    IF (unstbl(i)) THEN
      ! AM Niveau de ref du thermique
      ! AM          zxt=(t(i,1)-z(i,1)*0.5*RG/RCPD/(1.+RVTMP2*q(i,1)))
      ! AM     .         *(1.+RETV*q(i,1))
      zxt = (th_th(i)-zref*0.5*rg/rcpd/(1.+rvtmp2*qt_th(i)))* &
        (1.+retv*qt_th(i))
      phiminv(i) = (1.-binm*pblh(i)*unsobklen(i))**onet
      phihinv(i) = sqrt(1.-binh*pblh(i)*unsobklen(i))
      wm(i) = ustar(i)*phiminv(i)
      fak2(i) = wm(i)*pblh(i)*vk
      wstar(i) = (heatv(i)*rg*pblh(i)/zxt)**onet
      fak3(i) = fakn*wstar(i)/wm(i)
    ELSE
      wstar(i) = 0.
    END IF
    ! Computes Theta_e for thermal (all cases : to be modified)
    ! attention ajout therm(i) = virtuelle
    the_th(i) = th_th(i) + therm(i) + rlvcp*qt_th(i)
    ! ou:    The_th(i) = Th_th(i) + sqrt(b1)*khfs(i)/wm(i) + RLvCp*qT_th(i)
  END DO

  ! Main level loop to compute the diffusivities and
  ! counter-gradient terms:

  DO k = 2, isommet

    ! Find levels within boundary layer:

    DO i = 1, knon
      unslev(i) = .FALSE.
      stblev(i) = .FALSE.
      zm(i) = z(i, k-1)
      zp(i) = z(i, k)
      IF (zkmin==0.0 .AND. zp(i)>pblh(i)) zp(i) = pblh(i)
      IF (zm(i)<pblh(i)) THEN
        zmzp = 0.5*(zm(i)+zp(i))
        ! debug
        ! if (i.EQ.1864) then
        ! print*,'i,pblh(1864),obklen(1864)',i,pblh(i),obklen(i)
        ! endif

        zh(i) = zmzp/pblh(i)
        zl(i) = zmzp*unsobklen(i)
        zzh(i) = 0.
        IF (zh(i)<=1.0) zzh(i) = (1.-zh(i))**2

        ! stblev for points zm < plbh and stable and neutral
        ! unslev for points zm < plbh and unstable

        IF (unstbl(i)) THEN
          unslev(i) = .TRUE.
        ELSE
          stblev(i) = .TRUE.
        END IF
      END IF
    END DO
    ! print*,'fin calcul niveaux'

    ! Stable and neutral points; set diffusivities; counter-gradient
    ! terms zero for stable case:

    DO i = 1, knon
      IF (stblev(i)) THEN
        IF (zl(i)<=1.) THEN
          pblk(i) = fak1(i)*zh(i)*zzh(i)/(1.+betas*zl(i))
        ELSE
          pblk(i) = fak1(i)*zh(i)*zzh(i)/(betas+zl(i))
        END IF
        ! pcfm(i,k) = pblk(i)
        ! pcfh(i,k) = pcfm(i,k)
      END IF
    END DO

    ! unssrf, unstable within surface layer of pbl
    ! unsout, unstable within outer   layer of pbl

    DO i = 1, knon
      unssrf(i) = .FALSE.
      unsout(i) = .FALSE.
      IF (unslev(i)) THEN
        IF (zh(i)<sffrac) THEN
          unssrf(i) = .TRUE.
        ELSE
          unsout(i) = .TRUE.
        END IF
      END IF
    END DO

    ! Unstable for surface layer; counter-gradient terms zero

    DO i = 1, knon
      IF (unssrf(i)) THEN
        term = (1.-betam*zl(i))**onet
        pblk(i) = fak1(i)*zh(i)*zzh(i)*term
        pr(i) = term/sqrt(1.-betah*zl(i))
      END IF
    END DO
    ! print*,'fin counter-gradient terms zero'

    ! Unstable for outer layer; counter-gradient terms non-zero:

    DO i = 1, knon
      IF (unsout(i)) THEN
        pblk(i) = fak2(i)*zh(i)*zzh(i)
        ! cgs(i,k) = fak3(i)/(pblh(i)*wm(i))
        ! cgh(i,k) = khfs(i)*cgs(i,k)
        pr(i) = phiminv(i)/phihinv(i) + ccon*fak3(i)/fak
        ! cgq(i,k) = kqfs(i)*cgs(i,k)
      END IF
    END DO
    ! print*,'fin counter-gradient terms non zero'

    ! For all unstable layers, compute diffusivities and ctrgrad ter m

    ! DO i = 1, knon
    ! IF (unslev(i)) THEN
    ! pcfm(i,k) = pblk(i)
    ! pcfh(i,k) = pblk(i)/pr(i)
    ! etc cf original
    ! ENDIF
    ! ENDDO

    ! For all layers, compute integral info and CTEI

    DO i = 1, knon
      IF (check(i) .OR. omegafl(i)) THEN
        IF (.NOT. zsat(i)) THEN
          ! Th2 = The_th(i) - RLvCp*qT_th(i)
          th2 = th_th(i)
          t2 = th2*s(i, k)
          ! thermodyn functions
          zdelta = max(0., sign(1.,rtt-t2))
          qqsat = r2es*foeew(t2, zdelta)/pplay(i, k)
          qqsat = min(0.5, qqsat)
          zcor = 1./(1.-retv*qqsat)
          qqsat = qqsat*zcor

          IF (qqsat<qt_th(i)) THEN
            ! on calcule lcl
            IF (k==2) THEN
              plcl(i) = z(i, k)
            ELSE
              plcl(i) = z(i, k-1) + (z(i,k-1)-z(i,k))*(qt_th(i)-qsatbef(i))/( &
                qsatbef(i)-qqsat)
            END IF
            zsat(i) = .TRUE.
            tbef(i) = t2
          END IF

          qsatbef(i) = qqsat ! bug dans la version orig ???
        END IF
        ! amn ???? cette ligne a deja ete faite normalement ?
      END IF
      ! print*,'hbtm2 i,k=',i,k
    END DO
  END DO ! end of level loop
  ! IM 170305 BEG
  IF (1==0) THEN
    PRINT *, 'hbtm2  ok'
  END IF !(1.EQ.0) THEN
  ! IM 170305 END
  RETURN
END SUBROUTINE hbtm
