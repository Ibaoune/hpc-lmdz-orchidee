SUBROUTINE thermcell_2002(ngrid, nlay, ptimestep, iflag_thermals, pplay, &
    pplev, pphi, pu, pv, pt, po, pduadj, pdvadj, pdtadj, pdoadj, fm0, entr0, &
    fraca, wa_moy, r_aspect, l_mix, w2di, tho)

  USE dimphy
  USE write_field_phy
  IMPLICIT NONE

  ! =======================================================================

  ! Calcul du transport verticale dans la couche limite en presence
  ! de "thermiques" explicitement representes

  ! Réécriture à partir d'un listing papier à Habas, le 14/02/00

  ! le thermique est supposé homogène et dissipé par mélange avec
  ! son environnement. la longueur l_mix contrôle l'efficacité du
  ! mélange

  ! Le calcul du transport des différentes espèces se fait en prenant
  ! en compte:
  ! 1. un flux de masse montant
  ! 2. un flux de masse descendant
  ! 3. un entrainement
  ! 4. un detrainement

  ! =======================================================================

  ! -----------------------------------------------------------------------
  ! declarations:
  ! -------------

  include "YOMCST.h"

  ! arguments:
  ! ----------

  INTEGER ngrid, nlay, w2di, iflag_thermals
  REAL tho
  REAL ptimestep, l_mix, r_aspect
  REAL pt(ngrid, nlay), pdtadj(ngrid, nlay)
  REAL pu(ngrid, nlay), pduadj(ngrid, nlay)
  REAL pv(ngrid, nlay), pdvadj(ngrid, nlay)
  REAL po(ngrid, nlay), pdoadj(ngrid, nlay)
  REAL pplay(ngrid, nlay), pplev(ngrid, nlay+1)
  REAL pphi(ngrid, nlay)
  REAL fraca(ngrid, nlay+1), zw2(ngrid, nlay+1)

  INTEGER, SAVE :: idetr = 3, lev_out = 1
  !$OMP THREADPRIVATE(idetr,lev_out)

  ! local:
  ! ------

  INTEGER, SAVE :: dvdq = 0, flagdq = 0, dqimpl = 1
  LOGICAL, SAVE :: debut = .TRUE.
  !$OMP THREADPRIVATE(dvdq,flagdq,debut,dqimpl)

  INTEGER ig, k, l, lmax(klon, klev+1), lmaxa(klon), lmix(klon)
  REAL zmax(klon), zw, zz, ztva(klon, klev), zzz

  REAL zlev(klon, klev+1), zlay(klon, klev)
  REAL zh(klon, klev), zdhadj(klon, klev)
  REAL ztv(klon, klev)
  REAL zu(klon, klev), zv(klon, klev), zo(klon, klev)
  REAL wh(klon, klev+1)
  REAL wu(klon, klev+1), wv(klon, klev+1), wo(klon, klev+1)
  REAL zla(klon, klev+1)
  REAL zwa(klon, klev+1)
  REAL zld(klon, klev+1)
  REAL zwd(klon, klev+1)
  REAL zsortie(klon, klev)
  REAL zva(klon, klev)
  REAL zua(klon, klev)
  REAL zoa(klon, klev)

  REAL zha(klon, klev)
  REAL wa_moy(klon, klev+1)
  REAL fracc(klon, klev+1)
  REAL zf, zf2
  REAL thetath2(klon, klev), wth2(klon, klev)
  ! common/comtherm/thetath2,wth2

  REAL count_time

  LOGICAL sorties
  REAL rho(klon, klev), rhobarz(klon, klev+1), masse(klon, klev)
  REAL zpspsk(klon, klev)

  REAL wmax(klon, klev), wmaxa(klon)

  REAL wa(klon, klev, klev+1)
  REAL wd(klon, klev+1)
  REAL larg_part(klon, klev, klev+1)
  REAL fracd(klon, klev+1)
  REAL xxx(klon, klev+1)
  REAL larg_cons(klon, klev+1)
  REAL larg_detr(klon, klev+1)
  REAL fm0(klon, klev+1), entr0(klon, klev), detr(klon, klev)
  REAL pu_therm(klon, klev), pv_therm(klon, klev)
  REAL fm(klon, klev+1), entr(klon, klev)
  REAL fmc(klon, klev+1)

  CHARACTER (LEN=2) :: str2
  CHARACTER (LEN=10) :: str10

  CHARACTER (LEN=20) :: modname = 'thermcell2002'
  CHARACTER (LEN=80) :: abort_message

  LOGICAL vtest(klon), down

  EXTERNAL scopy

  INTEGER ncorrec, ll
  SAVE ncorrec
  DATA ncorrec/0/
  !$OMP THREADPRIVATE(ncorrec)


  ! -----------------------------------------------------------------------
  ! initialisation:
  ! ---------------

  sorties = .TRUE.
  IF (ngrid/=klon) THEN
    PRINT *
    PRINT *, 'STOP dans convadj'
    PRINT *, 'ngrid    =', ngrid
    PRINT *, 'klon  =', klon
  END IF

  ! -----------------------------------------------------------------------
  ! incrementation eventuelle de tendances precedentes:
  ! ---------------------------------------------------

  ! print*,'0 OK convect8'

  DO l = 1, nlay
    DO ig = 1, ngrid
      zpspsk(ig, l) = (pplay(ig,l)/pplev(ig,1))**rkappa
      zh(ig, l) = pt(ig, l)/zpspsk(ig, l)
      zu(ig, l) = pu(ig, l)
      zv(ig, l) = pv(ig, l)
      zo(ig, l) = po(ig, l)
      ztv(ig, l) = zh(ig, l)*(1.+0.61*zo(ig,l))
    END DO
  END DO

  ! print*,'1 OK convect8'
  ! --------------------


  ! + + + + + + + + + + +


  ! wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
  ! wh,wt,wo ...

  ! + + + + + + + + + + +  zh,zu,zv,zo,rho


  ! --------------------   zlev(1)
  ! \\\\\\\\\\\\\\\\\\\\



  ! -----------------------------------------------------------------------
  ! Calcul des altitudes des couches
  ! -----------------------------------------------------------------------

  IF (debut) THEN
    flagdq = (iflag_thermals-1000)/100
    dvdq = (iflag_thermals-(1000+flagdq*100))/10
    IF (flagdq==2) dqimpl = -1
    IF (flagdq==3) dqimpl = 1
    debut = .FALSE.
  END IF
  PRINT *, 'TH flag th ', iflag_thermals, flagdq, dvdq, dqimpl

  DO l = 2, nlay
    DO ig = 1, ngrid
      zlev(ig, l) = 0.5*(pphi(ig,l)+pphi(ig,l-1))/rg
    END DO
  END DO
  DO ig = 1, ngrid
    zlev(ig, 1) = 0.
    zlev(ig, nlay+1) = (2.*pphi(ig,klev)-pphi(ig,klev-1))/rg
  END DO
  DO l = 1, nlay
    DO ig = 1, ngrid
      zlay(ig, l) = pphi(ig, l)/rg
    END DO
  END DO

  ! print*,'2 OK convect8'
  ! -----------------------------------------------------------------------
  ! Calcul des densites
  ! -----------------------------------------------------------------------

  DO l = 1, nlay
    DO ig = 1, ngrid
      rho(ig, l) = pplay(ig, l)/(zpspsk(ig,l)*rd*zh(ig,l))
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      rhobarz(ig, l) = 0.5*(rho(ig,l)+rho(ig,l-1))
    END DO
  END DO

  DO k = 1, nlay
    DO l = 1, nlay + 1
      DO ig = 1, ngrid
        wa(ig, k, l) = 0.
      END DO
    END DO
  END DO

  ! print*,'3 OK convect8'
  ! ------------------------------------------------------------------
  ! Calcul de w2, quarre de w a partir de la cape
  ! a partir de w2, on calcule wa, vitesse de l'ascendance

  ! ATTENTION: Dans cette version, pour cause d'economie de memoire,
  ! w2 est stoke dans wa

  ! ATTENTION: dans convect8, on n'utilise le calcule des wa
  ! independants par couches que pour calculer l'entrainement
  ! a la base et la hauteur max de l'ascendance.

  ! Indicages:
  ! l'ascendance provenant du niveau k traverse l'interface l avec
  ! une vitesse wa(k,l).

  ! --------------------

  ! + + + + + + + + + +

  ! wa(k,l)   ----       --------------------    l
  ! /\
  ! /||\       + + + + + + + + + +
  ! ||
  ! ||        --------------------
  ! ||
  ! ||        + + + + + + + + + +
  ! ||
  ! ||        --------------------
  ! ||__
  ! |___      + + + + + + + + + +     k

  ! --------------------



  ! ------------------------------------------------------------------


  DO k = 1, nlay - 1
    DO ig = 1, ngrid
      wa(ig, k, k) = 0.
      wa(ig, k, k+1) = 2.*rg*(ztv(ig,k)-ztv(ig,k+1))/ztv(ig, k+1)* &
        (zlev(ig,k+1)-zlev(ig,k))
    END DO
    DO l = k + 1, nlay - 1
      DO ig = 1, ngrid
        wa(ig, k, l+1) = wa(ig, k, l) + 2.*rg*(ztv(ig,k)-ztv(ig,l))/ztv(ig, l &
          )*(zlev(ig,l+1)-zlev(ig,l))
      END DO
    END DO
    DO ig = 1, ngrid
      wa(ig, k, nlay+1) = 0.
    END DO
  END DO

  ! print*,'4 OK convect8'
  ! Calcul de la couche correspondant a la hauteur du thermique
  DO k = 1, nlay - 1
    DO ig = 1, ngrid
      lmax(ig, k) = k
    END DO
    DO l = nlay, k + 1, -1
      DO ig = 1, ngrid
        IF (wa(ig,k,l)<=1.E-10) lmax(ig, k) = l - 1
      END DO
    END DO
  END DO

  ! print*,'5 OK convect8'
  ! Calcule du w max du thermique
  DO k = 1, nlay
    DO ig = 1, ngrid
      wmax(ig, k) = 0.
    END DO
  END DO

  DO k = 1, nlay - 1
    DO l = k, nlay
      DO ig = 1, ngrid
        IF (l<=lmax(ig,k)) THEN
          wa(ig, k, l) = sqrt(wa(ig,k,l))
          wmax(ig, k) = max(wmax(ig,k), wa(ig,k,l))
        ELSE
          wa(ig, k, l) = 0.
        END IF
      END DO
    END DO
  END DO

  DO k = 1, nlay - 1
    DO ig = 1, ngrid
      pu_therm(ig, k) = sqrt(wmax(ig,k))
      pv_therm(ig, k) = sqrt(wmax(ig,k))
    END DO
  END DO

  ! print*,'6 OK convect8'
  ! Longueur caracteristique correspondant a la hauteur des thermiques.
  DO ig = 1, ngrid
    zmax(ig) = 500.
  END DO
  ! print*,'LMAX LMAX LMAX '
  DO k = 1, nlay - 1
    DO ig = 1, ngrid
      zmax(ig) = max(zmax(ig), zlev(ig,lmax(ig,k))-zlev(ig,k))
    END DO
    ! print*,k,lmax(1,k)
  END DO
  ! print*,'ZMAX ZMAX ZMAX ',zmax
  ! call dump2d(iim,jjm-1,zmax(2:ngrid-1),'ZMAX      ')

  ! print*,'OKl336'
  ! Calcul de l'entrainement.
  ! Le rapport d'aspect relie la largeur de l'ascendance a l'epaisseur
  ! de la couche d'alimentation en partant du principe que la vitesse
  ! maximum dans l'ascendance est la vitesse d'entrainement horizontale.
  DO k = 1, nlay
    DO ig = 1, ngrid
      zzz = rho(ig, k)*wmax(ig, k)*(zlev(ig,k+1)-zlev(ig,k))/ &
        (zmax(ig)*r_aspect)
      IF (w2di==2) THEN
        entr(ig, k) = entr(ig, k) + ptimestep*(zzz-entr(ig,k))/tho
      ELSE
        entr(ig, k) = zzz
      END IF
      ztva(ig, k) = ztv(ig, k)
    END DO
  END DO


  ! print*,'7 OK convect8'
  DO k = 1, klev + 1
    DO ig = 1, ngrid
      zw2(ig, k) = 0.
      fmc(ig, k) = 0.
      larg_cons(ig, k) = 0.
      larg_detr(ig, k) = 0.
      wa_moy(ig, k) = 0.
    END DO
  END DO

  ! print*,'8 OK convect8'
  DO ig = 1, ngrid
    lmaxa(ig) = 1
    lmix(ig) = 1
    wmaxa(ig) = 0.
  END DO


  ! print*,'OKl372'
  DO l = 1, nlay - 2
    DO ig = 1, ngrid
      ! if (zw2(ig,l).lt.1.e-10.and.ztv(ig,l).gt.ztv(ig,l+1)) then
      ! print*,'COUCOU ',l,zw2(ig,l),ztv(ig,l),ztv(ig,l+1)
      IF (zw2(ig,l)<1.E-10 .AND. ztv(ig,l)>ztv(ig,l+1) .AND. &
          entr(ig,l)>1.E-10) THEN
        ! print*,'COUCOU cas 1'
        ! Initialisation de l'ascendance
        ! lmix(ig)=1
        ztva(ig, l) = ztv(ig, l)
        fmc(ig, l) = 0.
        fmc(ig, l+1) = entr(ig, l)
        zw2(ig, l) = 0.
        ! if (.not.ztv(ig,l+1).gt.150.) then
        ! print*,'ig,l+1,ztv(ig,l+1)'
        ! print*, ig,l+1,ztv(ig,l+1)
        ! endif
        zw2(ig, l+1) = 2.*rg*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig, l+1)* &
          (zlev(ig,l+1)-zlev(ig,l))
        larg_detr(ig, l) = 0.
      ELSE IF (zw2(ig,l)>=1.E-10 .AND. fmc(ig,l)+entr(ig,l)>1.E-10) THEN
        ! Incrementation...
        fmc(ig, l+1) = fmc(ig, l) + entr(ig, l)
        ! if (.not.fmc(ig,l+1).gt.1.e-15) then
        ! print*,'ig,l+1,fmc(ig,l+1)'
        ! print*, ig,l+1,fmc(ig,l+1)
        ! print*,'Fmc ',(fmc(ig,ll),ll=1,klev+1)
        ! print*,'W2 ',(zw2(ig,ll),ll=1,klev+1)
        ! print*,'Tv ',(ztv(ig,ll),ll=1,klev)
        ! print*,'Entr ',(entr(ig,ll),ll=1,klev)
        ! endif
        ztva(ig, l) = (fmc(ig,l)*ztva(ig,l-1)+entr(ig,l)*ztv(ig,l))/ &
          fmc(ig, l+1)
        ! mise a jour de la vitesse ascendante (l'air entraine de la couche
        ! consideree commence avec une vitesse nulle).
        zw2(ig, l+1) = zw2(ig, l)*(fmc(ig,l)/fmc(ig,l+1))**2 + &
          2.*rg*(ztva(ig,l)-ztv(ig,l))/ztv(ig, l)*(zlev(ig,l+1)-zlev(ig,l))
      END IF
      IF (zw2(ig,l+1)<0.) THEN
        zw2(ig, l+1) = 0.
        lmaxa(ig) = l
      ELSE
        wa_moy(ig, l+1) = sqrt(zw2(ig,l+1))
      END IF
      IF (wa_moy(ig,l+1)>wmaxa(ig)) THEN
        ! lmix est le niveau de la couche ou w (wa_moy) est maximum
        lmix(ig) = l + 1
        wmaxa(ig) = wa_moy(ig, l+1)
      END IF
      ! print*,'COUCOU cas 2 LMIX=',lmix(ig),wa_moy(ig,l+1),wmaxa(ig)
    END DO
  END DO

  ! print*,'9 OK convect8'
  ! print*,'WA1 ',wa_moy

  ! determination de l'indice du debut de la mixed layer ou w decroit

  ! calcul de la largeur de chaque ascendance dans le cas conservatif.
  ! dans ce cas simple, on suppose que la largeur de l'ascendance provenant
  ! d'une couche est égale à la hauteur de la couche alimentante.
  ! La vitesse maximale dans l'ascendance est aussi prise comme estimation
  ! de la vitesse d'entrainement horizontal dans la couche alimentante.

  ! print*,'OKl439'
  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (l<=lmaxa(ig)) THEN
        zw = max(wa_moy(ig,l), 1.E-10)
        larg_cons(ig, l) = zmax(ig)*r_aspect*fmc(ig, l)/(rhobarz(ig,l)*zw)
      END IF
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (l<=lmaxa(ig)) THEN
        ! if (idetr.eq.0) then
        ! cette option est finalement en dur.
        larg_detr(ig, l) = sqrt(l_mix*zlev(ig,l))
        ! else if (idetr.eq.1) then
        ! larg_detr(ig,l)=larg_cons(ig,l)
        ! s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
        ! else if (idetr.eq.2) then
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        ! s            *sqrt(wa_moy(ig,l))
        ! else if (idetr.eq.4) then
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        ! s            *wa_moy(ig,l)
        ! endif
      END IF
    END DO
  END DO

  ! print*,'10 OK convect8'
  ! print*,'WA2 ',wa_moy
  ! calcul de la fraction de la maille concernée par l'ascendance en tenant
  ! compte de l'epluchage du thermique.

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (larg_cons(ig,l)>1.) THEN
        ! print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
        fraca(ig, l) = (larg_cons(ig,l)-larg_detr(ig,l))/(r_aspect*zmax(ig))
        IF (l>lmix(ig)) THEN
          xxx(ig, l) = (lmaxa(ig)+1.-l)/(lmaxa(ig)+1.-lmix(ig))
          IF (idetr==0) THEN
            fraca(ig, l) = fraca(ig, lmix(ig))
          ELSE IF (idetr==1) THEN
            fraca(ig, l) = fraca(ig, lmix(ig))*xxx(ig, l)
          ELSE IF (idetr==2) THEN
            fraca(ig, l) = fraca(ig, lmix(ig))*(1.-(1.-xxx(ig,l))**2)
          ELSE
            fraca(ig, l) = fraca(ig, lmix(ig))*xxx(ig, l)**2
          END IF
        END IF
        ! print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
        fraca(ig, l) = max(fraca(ig,l), 0.)
        fraca(ig, l) = min(fraca(ig,l), 0.5)
        fracd(ig, l) = 1. - fraca(ig, l)
        fracc(ig, l) = larg_cons(ig, l)/(r_aspect*zmax(ig))
      ELSE
        ! wa_moy(ig,l)=0.
        fraca(ig, l) = 0.
        fracc(ig, l) = 0.
        fracd(ig, l) = 1.
      END IF
    END DO
  END DO

  ! print*,'11 OK convect8'
  ! print*,'Ea3 ',wa_moy
  ! ------------------------------------------------------------------
  ! Calcul de fracd, wd
  ! somme wa - wd = 0
  ! ------------------------------------------------------------------


  DO ig = 1, ngrid
    fm(ig, 1) = 0.
    fm(ig, nlay+1) = 0.
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      fm(ig, l) = fraca(ig, l)*wa_moy(ig, l)*rhobarz(ig, l)
    END DO
    DO ig = 1, ngrid
      IF (fracd(ig,l)<0.1) THEN
        abort_message = 'fracd trop petit'
        CALL abort_physic(modname, abort_message, 1)
      ELSE
        ! vitesse descendante "diagnostique"
        wd(ig, l) = fm(ig, l)/(fracd(ig,l)*rhobarz(ig,l))
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      ! masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
      masse(ig, l) = (pplev(ig,l)-pplev(ig,l+1))/rg
    END DO
  END DO

  ! print*,'12 OK convect8'
  ! print*,'WA4 ',wa_moy
  ! c------------------------------------------------------------------
  ! calcul du transport vertical
  ! ------------------------------------------------------------------

  GO TO 4444
  ! print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
  DO l = 2, nlay - 1
    DO ig = 1, ngrid
      IF (fm(ig,l+1)*ptimestep>masse(ig,l) .AND. fm(ig,l+1)*ptimestep>masse( &
          ig,l+1)) THEN
        ! print*,'WARN!!! FM>M ig=',ig,' l=',l,'  FM='
        ! s         ,fm(ig,l+1)*ptimestep
        ! s         ,'   M=',masse(ig,l),masse(ig,l+1)
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (entr(ig,l)*ptimestep>masse(ig,l)) THEN
        ! print*,'WARN!!! E>M ig=',ig,' l=',l,'  E=='
        ! s         ,entr(ig,l)*ptimestep
        ! s         ,'   M=',masse(ig,l)
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (.NOT. fm(ig,l)>=0. .OR. .NOT. fm(ig,l)<=10.) THEN
        ! print*,'WARN!!! fm exagere ig=',ig,'   l=',l
        ! s         ,'   FM=',fm(ig,l)
      END IF
      IF (.NOT. masse(ig,l)>=1.E-10 .OR. .NOT. masse(ig,l)<=1.E4) THEN
        ! print*,'WARN!!! masse exagere ig=',ig,'   l=',l
        ! s         ,'   M=',masse(ig,l)
        ! print*,'rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)',
        ! s                 rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)
        ! print*,'zlev(ig,l+1),zlev(ig,l)'
        ! s                ,zlev(ig,l+1),zlev(ig,l)
        ! print*,'pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)'
        ! s                ,pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)
      END IF
      IF (.NOT. entr(ig,l)>=0. .OR. .NOT. entr(ig,l)<=10.) THEN
        ! print*,'WARN!!! entr exagere ig=',ig,'   l=',l
        ! s         ,'   E=',entr(ig,l)
      END IF
    END DO
  END DO

4444 CONTINUE
  ! print*,'OK 444 '

  IF (w2di==1) THEN
    fm0 = fm0 + ptimestep*(fm-fm0)/tho
    entr0 = entr0 + ptimestep*(entr-entr0)/tho
  ELSE
    fm0 = fm
    entr0 = entr
  END IF

  IF (flagdq==0) THEN
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zh, zdhadj, &
      zha)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zo, pdoadj, &
      zoa)
    PRINT *, 'THERMALS OPT 1'
  ELSE IF (flagdq==1) THEN
    CALL dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zh, &
      zdhadj, zha)
    CALL dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zo, &
      pdoadj, zoa)
    PRINT *, 'THERMALS OPT 2'
  ELSE
    CALL thermcell_dq(ngrid, nlay, dqimpl, ptimestep, fm0, entr0, masse, zh, &
      zdhadj, zha, lev_out)
    CALL thermcell_dq(ngrid, nlay, dqimpl, ptimestep, fm0, entr0, masse, zo, &
      pdoadj, zoa, lev_out)
    PRINT *, 'THERMALS OPT 3', dqimpl
  END IF

  PRINT *, 'TH VENT ', dvdq
  IF (dvdq==0) THEN
    ! print*,'TH VENT OK ',dvdq
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zu, pduadj, &
      zua)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zv, pdvadj, &
      zva)
  ELSE IF (dvdq==1) THEN
    CALL dvthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zmax, &
      zu, zv, pduadj, pdvadj, zua, zva)
  ELSE IF (dvdq==2) THEN
    CALL thermcell_dv2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, &
      zmax, zu, zv, pduadj, pdvadj, zua, zva, lev_out)
  ELSE IF (dvdq==3) THEN
    CALL thermcell_dq(ngrid, nlay, dqimpl, ptimestep, fm0, entr0, masse, zu, &
      pduadj, zua, lev_out)
    CALL thermcell_dq(ngrid, nlay, dqimpl, ptimestep, fm0, entr0, masse, zv, &
      pdvadj, zva, lev_out)
  END IF

  ! CALL writefield_phy('duadj',pduadj,klev)

  DO l = 1, nlay
    DO ig = 1, ngrid
      zf = 0.5*(fracc(ig,l)+fracc(ig,l+1))
      zf2 = zf/(1.-zf)
      thetath2(ig, l) = zf2*(zha(ig,l)-zh(ig,l))**2
      wth2(ig, l) = zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
    END DO
  END DO



  ! print*,'13 OK convect8'
  ! print*,'WA5 ',wa_moy
  DO l = 1, nlay
    DO ig = 1, ngrid
      pdtadj(ig, l) = zdhadj(ig, l)*zpspsk(ig, l)
    END DO
  END DO


  ! do l=1,nlay
  ! do ig=1,ngrid
  ! if(abs(pdtadj(ig,l))*86400..gt.500.) then
  ! print*,'WARN!!! ig=',ig,'  l=',l
  ! s         ,'   pdtadj=',pdtadj(ig,l)
  ! endif
  ! if(abs(pdoadj(ig,l))*86400..gt.1.) then
  ! print*,'WARN!!! ig=',ig,'  l=',l
  ! s         ,'   pdoadj=',pdoadj(ig,l)
  ! endif
  ! enddo
  ! enddo

  ! print*,'14 OK convect8'
  ! ------------------------------------------------------------------
  ! Calculs pour les sorties
  ! ------------------------------------------------------------------

  IF (sorties) THEN
    DO l = 1, nlay
      DO ig = 1, ngrid
        zla(ig, l) = (1.-fracd(ig,l))*zmax(ig)
        zld(ig, l) = fracd(ig, l)*zmax(ig)
        IF (1.-fracd(ig,l)>1.E-10) zwa(ig, l) = wd(ig, l)*fracd(ig, l)/ &
          (1.-fracd(ig,l))
      END DO
    END DO

    DO l = 1, nlay
      DO ig = 1, ngrid
        detr(ig, l) = fm(ig, l) + entr(ig, l) - fm(ig, l+1)
        IF (detr(ig,l)<0.) THEN
          entr(ig, l) = entr(ig, l) - detr(ig, l)
          detr(ig, l) = 0.
          ! print*,'WARNING !!! detrainement negatif ',ig,l
        END IF
      END DO
    END DO
  END IF

  ! print*,'15 OK convect8'


  ! if(wa_moy(1,4).gt.1.e-10) stop

  ! print*,'19 OK convect8'
  RETURN
END SUBROUTINE thermcell_2002

SUBROUTINE thermcell_cld(ngrid, nlay, ptimestep, pplay, pplev, pphi, zlev, &
    debut, pu, pv, pt, po, pduadj, pdvadj, pdtadj, pdoadj, fm0, entr0, zqla, &
    lmax, zmax_sec, wmax_sec, zw_sec, lmix_sec, ratqscth, ratqsdiff & ! s
                                                                      ! ,pu_therm,pv_therm
    , r_aspect, l_mix, w2di, tho)

  USE dimphy
  IMPLICIT NONE

  ! =======================================================================

  ! Calcul du transport verticale dans la couche limite en presence
  ! de "thermiques" explicitement representes

  ! Réécriture à partir d'un listing papier à Habas, le 14/02/00

  ! le thermique est supposé homogène et dissipé par mélange avec
  ! son environnement. la longueur l_mix contrôle l'efficacité du
  ! mélange

  ! Le calcul du transport des différentes espèces se fait en prenant
  ! en compte:
  ! 1. un flux de masse montant
  ! 2. un flux de masse descendant
  ! 3. un entrainement
  ! 4. un detrainement

  ! =======================================================================

  ! -----------------------------------------------------------------------
  ! declarations:
  ! -------------

  include "YOMCST.h"
  include "YOETHF.h"
  include "FCTTRE.h"

  ! arguments:
  ! ----------

  INTEGER ngrid, nlay, w2di
  REAL tho
  REAL ptimestep, l_mix, r_aspect
  REAL pt(ngrid, nlay), pdtadj(ngrid, nlay)
  REAL pu(ngrid, nlay), pduadj(ngrid, nlay)
  REAL pv(ngrid, nlay), pdvadj(ngrid, nlay)
  REAL po(ngrid, nlay), pdoadj(ngrid, nlay)
  REAL pplay(ngrid, nlay), pplev(ngrid, nlay+1)
  REAL pphi(ngrid, nlay)

  INTEGER idetr
  SAVE idetr
  DATA idetr/3/
  !$OMP THREADPRIVATE(idetr)

  ! local:
  ! ------

  INTEGER ig, k, l, lmaxa(klon), lmix(klon)
  REAL zsortie1d(klon)
  ! CR: on remplace lmax(klon,klev+1)
  INTEGER lmax(klon), lmin(klon), lentr(klon)
  REAL linter(klon)
  REAL zmix(klon), fracazmix(klon)
  REAL alpha
  SAVE alpha
  DATA alpha/1./
  !$OMP THREADPRIVATE(alpha)

  ! RC
  REAL zmax(klon), zw, zz, zw2(klon, klev+1), ztva(klon, klev), zzz
  REAL zmax_sec(klon)
  REAL zmax_sec2(klon)
  REAL zw_sec(klon, klev+1)
  INTEGER lmix_sec(klon)
  REAL w_est(klon, klev+1)
  ! on garde le zmax du pas de temps precedent
  ! real zmax0(klon)
  ! save zmax0
  ! real zmix0(klon)
  ! save zmix0
  REAL, SAVE, ALLOCATABLE :: zmax0(:), zmix0(:)
  !$OMP THREADPRIVATE(zmax0, zmix0)

  REAL zlev(klon, klev+1), zlay(klon, klev)
  REAL deltaz(klon, klev)
  REAL zh(klon, klev), zdhadj(klon, klev)
  REAL zthl(klon, klev), zdthladj(klon, klev)
  REAL ztv(klon, klev)
  REAL zu(klon, klev), zv(klon, klev), zo(klon, klev)
  REAL zl(klon, klev)
  REAL wh(klon, klev+1)
  REAL wu(klon, klev+1), wv(klon, klev+1), wo(klon, klev+1)
  REAL zla(klon, klev+1)
  REAL zwa(klon, klev+1)
  REAL zld(klon, klev+1)
  REAL zwd(klon, klev+1)
  REAL zsortie(klon, klev)
  REAL zva(klon, klev)
  REAL zua(klon, klev)
  REAL zoa(klon, klev)

  REAL zta(klon, klev)
  REAL zha(klon, klev)
  REAL wa_moy(klon, klev+1)
  REAL fraca(klon, klev+1)
  REAL fracc(klon, klev+1)
  REAL zf, zf2
  REAL thetath2(klon, klev), wth2(klon, klev), wth3(klon, klev)
  REAL q2(klon, klev)
  REAL dtheta(klon, klev)
  ! common/comtherm/thetath2,wth2

  REAL ratqscth(klon, klev)
  REAL sum
  REAL sumdiff
  REAL ratqsdiff(klon, klev)
  REAL count_time
  INTEGER ialt

  LOGICAL sorties
  REAL rho(klon, klev), rhobarz(klon, klev+1), masse(klon, klev)
  REAL zpspsk(klon, klev)

  ! real wmax(klon,klev),wmaxa(klon)
  REAL wmax(klon), wmaxa(klon)
  REAL wmax_sec(klon)
  REAL wmax_sec2(klon)
  REAL wa(klon, klev, klev+1)
  REAL wd(klon, klev+1)
  REAL larg_part(klon, klev, klev+1)
  REAL fracd(klon, klev+1)
  REAL xxx(klon, klev+1)
  REAL larg_cons(klon, klev+1)
  REAL larg_detr(klon, klev+1)
  REAL fm0(klon, klev+1), entr0(klon, klev), detr(klon, klev)
  REAL massetot(klon, klev)
  REAL detr0(klon, klev)
  REAL alim0(klon, klev)
  REAL pu_therm(klon, klev), pv_therm(klon, klev)
  REAL fm(klon, klev+1), entr(klon, klev)
  REAL fmc(klon, klev+1)

  REAL zcor, zdelta, zcvm5, qlbef
  REAL tbef(klon), qsatbef(klon)
  REAL dqsat_dt, dt, num, denom
  REAL reps, rlvcp, ddt0
  REAL ztla(klon, klev), zqla(klon, klev), zqta(klon, klev)
  ! CR niveau de condensation
  REAL nivcon(klon)
  REAL zcon(klon)
  REAL zqsat(klon, klev)
  REAL zqsatth(klon, klev)
  PARAMETER (ddt0=.01)


  ! CR:nouvelles variables
  REAL f_star(klon, klev+1), entr_star(klon, klev)
  REAL detr_star(klon, klev)
  REAL alim_star_tot(klon), alim_star2(klon)
  REAL entr_star_tot(klon)
  REAL detr_star_tot(klon)
  REAL alim_star(klon, klev)
  REAL alim(klon, klev)
  REAL nu(klon, klev)
  REAL nu_e(klon, klev)
  REAL nu_min
  REAL nu_max
  REAL nu_r
  REAL f(klon)
  ! real f(klon), f0(klon)
  ! save f0
  REAL, SAVE, ALLOCATABLE :: f0(:)
  !$OMP THREADPRIVATE(f0)

  REAL f_old
  REAL zlevinter(klon)
  LOGICAL, SAVE :: first = .TRUE.
  !$OMP THREADPRIVATE(first)
  ! data first /.false./
  ! save first
  LOGICAL nuage
  ! save nuage
  LOGICAL boucle
  LOGICAL therm
  LOGICAL debut
  LOGICAL rale
  INTEGER test(klon)
  INTEGER signe_zw2
  ! RC

  CHARACTER *2 str2
  CHARACTER *10 str10

  CHARACTER (LEN=20) :: modname = 'thermcell_cld'
  CHARACTER (LEN=80) :: abort_message

  LOGICAL vtest(klon), down
  LOGICAL zsat(klon)

  EXTERNAL scopy

  INTEGER ncorrec, ll
  SAVE ncorrec
  DATA ncorrec/0/
  !$OMP THREADPRIVATE(ncorrec)



  ! -----------------------------------------------------------------------
  ! initialisation:
  ! ---------------

  IF (first) THEN
    ALLOCATE (zmix0(klon))
    ALLOCATE (zmax0(klon))
    ALLOCATE (f0(klon))
    first = .FALSE.
  END IF

  sorties = .FALSE.
  ! print*,'NOUVEAU DETR PLUIE '
  IF (ngrid/=klon) THEN
    PRINT *
    PRINT *, 'STOP dans convadj'
    PRINT *, 'ngrid    =', ngrid
    PRINT *, 'klon  =', klon
  END IF

  ! Initialisation
  rlvcp = rlvtt/rcpd
  reps = rd/rv
  ! initialisations de zqsat
  DO ll = 1, nlay
    DO ig = 1, ngrid
      zqsat(ig, ll) = 0.
      zqsatth(ig, ll) = 0.
    END DO
  END DO

  ! on met le first a true pour le premier passage de la journée
  DO ig = 1, klon
    test(ig) = 0
  END DO
  IF (debut) THEN
    DO ig = 1, klon
      test(ig) = 1
      f0(ig) = 0.
      zmax0(ig) = 0.
    END DO
  END IF
  DO ig = 1, klon
    IF ((.NOT. debut) .AND. (f0(ig)<1.E-10)) THEN
      test(ig) = 1
    END IF
  END DO
  ! do ig=1,klon
  ! print*,'test(ig)',test(ig),zmax0(ig)
  ! enddo
  nuage = .FALSE.
  ! -----------------------------------------------------------------------
  ! AM Calcul de T,q,ql a partir de Tl et qT
  ! ---------------------------------------------------

  ! Pr Tprec=Tl calcul de qsat
  ! Si qsat>qT T=Tl, q=qT
  ! Sinon DDT=(-Tprec+Tl+RLVCP (qT-qsat(T')) / (1+RLVCP dqsat/dt)
  ! On cherche DDT < DDT0

  ! defaut
  DO ll = 1, nlay
    DO ig = 1, ngrid
      zo(ig, ll) = po(ig, ll)
      zl(ig, ll) = 0.
      zh(ig, ll) = pt(ig, ll)
    END DO
  END DO
  DO ig = 1, ngrid
    zsat(ig) = .FALSE.
  END DO


  DO ll = 1, nlay
    ! les points insatures sont definitifs
    DO ig = 1, ngrid
      tbef(ig) = pt(ig, ll)
      zdelta = max(0., sign(1.,rtt-tbef(ig)))
      qsatbef(ig) = r2es*foeew(tbef(ig), zdelta)/pplev(ig, ll)
      qsatbef(ig) = min(0.5, qsatbef(ig))
      zcor = 1./(1.-retv*qsatbef(ig))
      qsatbef(ig) = qsatbef(ig)*zcor
      zsat(ig) = (max(0.,po(ig,ll)-qsatbef(ig))>1.E-10)
    END DO

    DO ig = 1, ngrid
      IF (zsat(ig) .AND. (1==1)) THEN
        qlbef = max(0., po(ig,ll)-qsatbef(ig))
        ! si sature: ql est surestime, d'ou la sous-relax
        dt = 0.5*rlvcp*qlbef
        ! write(18,*),'DT0=',DT
        ! on pourra enchainer 2 ou 3 calculs sans Do while
        DO WHILE (abs(dt)>ddt0)
          ! il faut verifier si c,a conserve quand on repasse en insature ...
          tbef(ig) = tbef(ig) + dt
          zdelta = max(0., sign(1.,rtt-tbef(ig)))
          qsatbef(ig) = r2es*foeew(tbef(ig), zdelta)/pplev(ig, ll)
          qsatbef(ig) = min(0.5, qsatbef(ig))
          zcor = 1./(1.-retv*qsatbef(ig))
          qsatbef(ig) = qsatbef(ig)*zcor
          ! on veut le signe de qlbef
          qlbef = po(ig, ll) - qsatbef(ig)
          zdelta = max(0., sign(1.,rtt-tbef(ig)))
          zcvm5 = r5les*(1.-zdelta) + r5ies*zdelta
          zcor = 1./(1.-retv*qsatbef(ig))
          dqsat_dt = foede(tbef(ig), zdelta, zcvm5, qsatbef(ig), zcor)
          num = -tbef(ig) + pt(ig, ll) + rlvcp*qlbef
          denom = 1. + rlvcp*dqsat_dt
          IF (denom<1.E-10) THEN
            PRINT *, 'pb denom'
          END IF
          dt = num/denom
        END DO
        ! on ecrit de maniere conservative (sat ou non)
        zl(ig, ll) = max(0., qlbef)
        ! T = Tl +Lv/Cp ql
        zh(ig, ll) = pt(ig, ll) + rlvcp*zl(ig, ll)
        zo(ig, ll) = po(ig, ll) - zl(ig, ll)
      END IF
      ! on ecrit zqsat
      zqsat(ig, ll) = qsatbef(ig)
    END DO
  END DO
  ! AM fin

  ! -----------------------------------------------------------------------
  ! incrementation eventuelle de tendances precedentes:
  ! ---------------------------------------------------

  ! print*,'0 OK convect8'

  DO l = 1, nlay
    DO ig = 1, ngrid
      zpspsk(ig, l) = (pplay(ig,l)/100000.)**rkappa
      ! zpspsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**RKAPPA
      ! zh(ig,l)=pt(ig,l)/zpspsk(ig,l)
      zu(ig, l) = pu(ig, l)
      zv(ig, l) = pv(ig, l)
      ! zo(ig,l)=po(ig,l)
      ! ztv(ig,l)=zh(ig,l)*(1.+0.61*zo(ig,l))
      ! AM attention zh est maintenant le profil de T et plus le profil de
      ! theta !

      ! T-> Theta
      ztv(ig, l) = zh(ig, l)/zpspsk(ig, l)
      ! AM Theta_v
      ztv(ig, l) = ztv(ig, l)*(1.+retv*(zo(ig,l))-zl(ig,l))
      ! AM Thetal
      zthl(ig, l) = pt(ig, l)/zpspsk(ig, l)

    END DO
  END DO

  ! print*,'1 OK convect8'
  ! --------------------


  ! + + + + + + + + + + +


  ! wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
  ! wh,wt,wo ...

  ! + + + + + + + + + + +  zh,zu,zv,zo,rho


  ! --------------------   zlev(1)
  ! \\\\\\\\\\\\\\\\\\\\



  ! -----------------------------------------------------------------------
  ! Calcul des altitudes des couches
  ! -----------------------------------------------------------------------

  DO l = 2, nlay
    DO ig = 1, ngrid
      zlev(ig, l) = 0.5*(pphi(ig,l)+pphi(ig,l-1))/rg
    END DO
  END DO
  DO ig = 1, ngrid
    zlev(ig, 1) = 0.
    zlev(ig, nlay+1) = (2.*pphi(ig,klev)-pphi(ig,klev-1))/rg
  END DO
  DO l = 1, nlay
    DO ig = 1, ngrid
      zlay(ig, l) = pphi(ig, l)/rg
    END DO
  END DO
  ! calcul de deltaz
  DO l = 1, nlay
    DO ig = 1, ngrid
      deltaz(ig, l) = zlev(ig, l+1) - zlev(ig, l)
    END DO
  END DO

  ! print*,'2 OK convect8'
  ! -----------------------------------------------------------------------
  ! Calcul des densites
  ! -----------------------------------------------------------------------

  DO l = 1, nlay
    DO ig = 1, ngrid
      ! rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*zh(ig,l))
      rho(ig, l) = pplay(ig, l)/(zpspsk(ig,l)*rd*ztv(ig,l))
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      rhobarz(ig, l) = 0.5*(rho(ig,l)+rho(ig,l-1))
    END DO
  END DO

  DO k = 1, nlay
    DO l = 1, nlay + 1
      DO ig = 1, ngrid
        wa(ig, k, l) = 0.
      END DO
    END DO
  END DO
  ! Cr:ajout:calcul de la masse
  DO l = 1, nlay
    DO ig = 1, ngrid
      ! masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
      masse(ig, l) = (pplev(ig,l)-pplev(ig,l+1))/rg
    END DO
  END DO
  ! print*,'3 OK convect8'
  ! ------------------------------------------------------------------
  ! Calcul de w2, quarre de w a partir de la cape
  ! a partir de w2, on calcule wa, vitesse de l'ascendance

  ! ATTENTION: Dans cette version, pour cause d'economie de memoire,
  ! w2 est stoke dans wa

  ! ATTENTION: dans convect8, on n'utilise le calcule des wa
  ! independants par couches que pour calculer l'entrainement
  ! a la base et la hauteur max de l'ascendance.

  ! Indicages:
  ! l'ascendance provenant du niveau k traverse l'interface l avec
  ! une vitesse wa(k,l).

  ! --------------------

  ! + + + + + + + + + +

  ! wa(k,l)   ----       --------------------    l
  ! /\
  ! /||\       + + + + + + + + + +
  ! ||
  ! ||        --------------------
  ! ||
  ! ||        + + + + + + + + + +
  ! ||
  ! ||        --------------------
  ! ||__
  ! |___      + + + + + + + + + +     k

  ! --------------------



  ! ------------------------------------------------------------------

  ! CR: ponderation entrainement des couches instables
  ! def des alim_star tels que alim=f*alim_star
  DO l = 1, klev
    DO ig = 1, ngrid
      alim_star(ig, l) = 0.
      alim(ig, l) = 0.
    END DO
  END DO
  ! determination de la longueur de la couche d entrainement
  DO ig = 1, ngrid
    lentr(ig) = 1
  END DO

  ! on ne considere que les premieres couches instables
  therm = .FALSE.
  DO k = nlay - 2, 1, -1
    DO ig = 1, ngrid
      IF (ztv(ig,k)>ztv(ig,k+1) .AND. ztv(ig,k+1)<=ztv(ig,k+2)) THEN
        lentr(ig) = k + 1
        therm = .TRUE.
      END IF
    END DO
  END DO

  ! determination du lmin: couche d ou provient le thermique
  DO ig = 1, ngrid
    lmin(ig) = 1
  END DO
  DO ig = 1, ngrid
    DO l = nlay, 2, -1
      IF (ztv(ig,l-1)>ztv(ig,l)) THEN
        lmin(ig) = l - 1
      END IF
    END DO
  END DO

  ! definition de l'entrainement des couches
  DO l = 1, klev - 1
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. l>=lmin(ig) .AND. l<lentr(ig)) THEN
        ! def possibles pour alim_star: zdthetadz, dthetadz, zdtheta
        alim_star(ig, l) = max((ztv(ig,l)-ztv(ig,l+1)), 0.) & ! s
                                                              ! *(zlev(ig,l+1)-zlev(ig,l))
          *sqrt(zlev(ig,l+1))
        ! alim_star(ig,l)=zlev(ig,l+1)*(1.-(zlev(ig,l+1)
        ! s                         /zlev(ig,lentr(ig)+2)))**(3./2.)
      END IF
    END DO
  END DO

  ! pas de thermique si couche 1 stable
  DO ig = 1, ngrid
    ! if (lmin(ig).gt.1) then
    ! CRnouveau test
    IF (alim_star(ig,1)<1.E-10) THEN
      DO l = 1, klev
        alim_star(ig, l) = 0.
      END DO
    END IF
  END DO
  ! calcul de l entrainement total
  DO ig = 1, ngrid
    alim_star_tot(ig) = 0.
    entr_star_tot(ig) = 0.
    detr_star_tot(ig) = 0.
  END DO
  DO ig = 1, ngrid
    DO k = 1, klev
      alim_star_tot(ig) = alim_star_tot(ig) + alim_star(ig, k)
    END DO
  END DO

  ! Calcul entrainement normalise
  DO ig = 1, ngrid
    IF (alim_star_tot(ig)>1.E-10) THEN
      ! do l=1,lentr(ig)
      DO l = 1, klev
        ! def possibles pour entr_star: zdthetadz, dthetadz, zdtheta
        alim_star(ig, l) = alim_star(ig, l)/alim_star_tot(ig)
      END DO
    END IF
  END DO

  ! print*,'fin calcul alim_star'

  ! AM:initialisations
  DO k = 1, nlay
    DO ig = 1, ngrid
      ztva(ig, k) = ztv(ig, k)
      ztla(ig, k) = zthl(ig, k)
      zqla(ig, k) = 0.
      zqta(ig, k) = po(ig, k)
      zsat(ig) = .FALSE.
    END DO
  END DO
  DO k = 1, klev
    DO ig = 1, ngrid
      detr_star(ig, k) = 0.
      entr_star(ig, k) = 0.
      detr(ig, k) = 0.
      entr(ig, k) = 0.
    END DO
  END DO
  ! print*,'7 OK convect8'
  DO k = 1, klev + 1
    DO ig = 1, ngrid
      zw2(ig, k) = 0.
      fmc(ig, k) = 0.
      ! CR
      f_star(ig, k) = 0.
      ! RC
      larg_cons(ig, k) = 0.
      larg_detr(ig, k) = 0.
      wa_moy(ig, k) = 0.
    END DO
  END DO

  ! n     print*,'8 OK convect8'
  DO ig = 1, ngrid
    linter(ig) = 1.
    lmaxa(ig) = 1
    lmix(ig) = 1
    wmaxa(ig) = 0.
  END DO

  nu_min = l_mix
  nu_max = 1000.
  ! do ig=1,ngrid
  ! nu_max=wmax_sec(ig)
  ! enddo
  DO ig = 1, ngrid
    DO k = 1, klev
      nu(ig, k) = 0.
      nu_e(ig, k) = 0.
    END DO
  END DO
  ! Calcul de l'excès de température du à la diffusion turbulente
  DO ig = 1, ngrid
    DO l = 1, klev
      dtheta(ig, l) = 0.
    END DO
  END DO
  DO ig = 1, ngrid
    DO l = 1, lentr(ig) - 1
      dtheta(ig, l) = sqrt(10.*0.4*zlev(ig,l+1)**2*1.*((ztv(ig,l+1)- &
        ztv(ig,l))/(zlev(ig,l+1)-zlev(ig,l)))**2)
    END DO
  END DO
  ! do l=1,nlay-2
  DO l = 1, klev - 1
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. alim_star(ig,l)>1.E-10 .AND. &
          zw2(ig,l)<1E-10) THEN
        ! AM
        ! test:on rajoute un excès de T dans couche alim
        ! ztla(ig,l)=zthl(ig,l)+dtheta(ig,l)
        ztla(ig, l) = zthl(ig, l)
        ! test: on rajoute un excès de q dans la couche alim
        ! zqta(ig,l)=po(ig,l)+0.001
        zqta(ig, l) = po(ig, l)
        zqla(ig, l) = zl(ig, l)
        ! AM
        f_star(ig, l+1) = alim_star(ig, l)
        ! test:calcul de dteta
        zw2(ig, l+1) = 2.*rg*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig, l+1)* &
          (zlev(ig,l+1)-zlev(ig,l))*0.4*pphi(ig, l)/(pphi(ig,l+1)-pphi(ig,l))
        w_est(ig, l+1) = zw2(ig, l+1)
        larg_detr(ig, l) = 0.
        ! print*,'coucou boucle 1'
      ELSE IF ((zw2(ig,l)>=1E-10) .AND. (f_star(ig,l)+alim_star(ig, &
          l))>1.E-10) THEN
        ! print*,'coucou boucle 2'
        ! estimation du detrainement a partir de la geometrie du pas
        ! precedent
        IF ((test(ig)==1) .OR. ((.NOT. debut) .AND. (f0(ig)<1.E-10))) THEN
          detr_star(ig, l) = 0.
          entr_star(ig, l) = 0.
          ! print*,'coucou test(ig)',test(ig),f0(ig),zmax0(ig)
        ELSE
          ! print*,'coucou debut detr'
          ! tests sur la definition du detr
          IF (zqla(ig,l-1)>1.E-10) THEN
            nuage = .TRUE.
          END IF

          w_est(ig, l+1) = zw2(ig, l)*((f_star(ig,l))**2)/(f_star(ig,l)+ &
            alim_star(ig,l))**2 + 2.*rg*(ztva(ig,l-1)-ztv(ig,l))/ztv(ig, l)*( &
            zlev(ig,l+1)-zlev(ig,l))
          IF (w_est(ig,l+1)<0.) THEN
            w_est(ig, l+1) = zw2(ig, l)
          END IF
          IF (l>2) THEN
            IF ((w_est(ig,l+1)>w_est(ig,l)) .AND. (zlev(ig, &
                l+1)<zmax_sec(ig)) .AND. (zqla(ig,l-1)<1.E-10)) THEN
              detr_star(ig, l) = max(0., (rhobarz(ig, &
                l+1)*sqrt(w_est(ig,l+1))*sqrt(nu(ig,l)* &
                zlev(ig,l+1))-rhobarz(ig,l)*sqrt(w_est(ig,l))*sqrt(nu(ig,l)* &
                zlev(ig,l)))/(r_aspect*zmax_sec(ig)))
            ELSE IF ((zlev(ig,l+1)<zmax_sec(ig)) .AND. (zqla(ig, &
                l-1)<1.E-10)) THEN
              detr_star(ig, l) = -f0(ig)*f_star(ig, lmix(ig))/(rhobarz(ig, &
                lmix(ig))*wmaxa(ig))*(rhobarz(ig,l+1)*sqrt(w_est(ig, &
                l+1))*((zmax_sec(ig)-zlev(ig,l+1))/((zmax_sec(ig)-zlev(ig, &
                lmix(ig)))))**2.-rhobarz(ig,l)*sqrt(w_est(ig, &
                l))*((zmax_sec(ig)-zlev(ig,l))/((zmax_sec(ig)-zlev(ig,lmix(ig &
                )))))**2.)
            ELSE
              detr_star(ig, l) = 0.002*f0(ig)*f_star(ig, l)* &
                (zlev(ig,l+1)-zlev(ig,l))

            END IF
          ELSE
            detr_star(ig, l) = 0.
          END IF

          detr_star(ig, l) = detr_star(ig, l)/f0(ig)
          IF (nuage) THEN
            entr_star(ig, l) = 0.4*detr_star(ig, l)
          ELSE
            entr_star(ig, l) = 0.4*detr_star(ig, l)
          END IF

          IF ((detr_star(ig,l))>f_star(ig,l)) THEN
            detr_star(ig, l) = f_star(ig, l)
            ! entr_star(ig,l)=0.
          END IF

          IF ((l<lentr(ig))) THEN
            entr_star(ig, l) = 0.
            ! detr_star(ig,l)=0.
          END IF

          ! print*,'ok detr_star'
        END IF
        ! prise en compte du detrainement dans le calcul du flux
        f_star(ig, l+1) = f_star(ig, l) + alim_star(ig, l) + &
          entr_star(ig, l) - detr_star(ig, l)
        ! test
        ! if (f_star(ig,l+1).lt.0.) then
        ! f_star(ig,l+1)=0.
        ! entr_star(ig,l)=0.
        ! detr_star(ig,l)=f_star(ig,l)+alim_star(ig,l)
        ! endif
        ! test sur le signe de f_star
        IF (f_star(ig,l+1)>1.E-10) THEN
          ! then
          ! test
          ! if (((f_star(ig,l+1)+detr_star(ig,l)).gt.1.e-10)) then
          ! AM on melange Tl et qt du thermique
          ! on rajoute un excès de T dans la couche alim
          ! if (l.lt.lentr(ig)) then
          ! ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+
          ! s
          ! (alim_star(ig,l)+entr_star(ig,l))*(zthl(ig,l)+dtheta(ig,l)))
          ! s     /(f_star(ig,l+1)+detr_star(ig,l))
          ! else
          ztla(ig, l) = (f_star(ig,l)*ztla(ig,l-1)+(alim_star(ig, &
            l)+entr_star(ig,l))*zthl(ig,l))/(f_star(ig,l+1)+detr_star(ig,l))
          ! s                    /(f_star(ig,l+1))
          ! endif
          ! on rajoute un excès de q dans la couche alim
          ! if (l.lt.lentr(ig)) then
          ! zqta(ig,l)=(f_star(ig,l)*zqta(ig,l-1)+
          ! s           (alim_star(ig,l)+entr_star(ig,l))*(po(ig,l)+0.001))
          ! s                 /(f_star(ig,l+1)+detr_star(ig,l))
          ! else
          zqta(ig, l) = (f_star(ig,l)*zqta(ig,l-1)+(alim_star(ig, &
            l)+entr_star(ig,l))*po(ig,l))/(f_star(ig,l+1)+detr_star(ig,l))
          ! s                   /(f_star(ig,l+1))
          ! endif
          ! AM on en deduit thetav et ql du thermique
          ! CR test
          ! Tbef(ig)=ztla(ig,l)*zpspsk(ig,l)
          tbef(ig) = ztla(ig, l)*zpspsk(ig, l)
          zdelta = max(0., sign(1.,rtt-tbef(ig)))
          qsatbef(ig) = r2es*foeew(tbef(ig), zdelta)/pplev(ig, l)
          qsatbef(ig) = min(0.5, qsatbef(ig))
          zcor = 1./(1.-retv*qsatbef(ig))
          qsatbef(ig) = qsatbef(ig)*zcor
          zsat(ig) = (max(0.,zqta(ig,l)-qsatbef(ig))>1.E-10)

          IF (zsat(ig) .AND. (1==1)) THEN
            qlbef = max(0., zqta(ig,l)-qsatbef(ig))
            dt = 0.5*rlvcp*qlbef
            ! write(17,*)'DT0=',DT
            DO WHILE (abs(dt)>ddt0)
              ! print*,'aie'
              tbef(ig) = tbef(ig) + dt
              zdelta = max(0., sign(1.,rtt-tbef(ig)))
              qsatbef(ig) = r2es*foeew(tbef(ig), zdelta)/pplev(ig, l)
              qsatbef(ig) = min(0.5, qsatbef(ig))
              zcor = 1./(1.-retv*qsatbef(ig))
              qsatbef(ig) = qsatbef(ig)*zcor
              qlbef = zqta(ig, l) - qsatbef(ig)

              zdelta = max(0., sign(1.,rtt-tbef(ig)))
              zcvm5 = r5les*(1.-zdelta) + r5ies*zdelta
              zcor = 1./(1.-retv*qsatbef(ig))
              dqsat_dt = foede(tbef(ig), zdelta, zcvm5, qsatbef(ig), zcor)
              num = -tbef(ig) + ztla(ig, l)*zpspsk(ig, l) + rlvcp*qlbef
              denom = 1. + rlvcp*dqsat_dt
              IF (denom<1.E-10) THEN
                PRINT *, 'pb denom'
              END IF
              dt = num/denom
              ! write(17,*)'DT=',DT
            END DO
            zqla(ig, l) = max(0., zqta(ig,l)-qsatbef(ig))
            zqla(ig, l) = max(0., qlbef)
            ! zqla(ig,l)=0.
          END IF
          ! zqla(ig,l) = max(0.,zqta(ig,l)-qsatbef(ig))

          ! on ecrit de maniere conservative (sat ou non)
          ! T = Tl +Lv/Cp ql
          ! CR rq utilisation de humidite specifique ou rapport de melange?
          ztva(ig, l) = ztla(ig, l)*zpspsk(ig, l) + rlvcp*zqla(ig, l)
          ztva(ig, l) = ztva(ig, l)/zpspsk(ig, l)
          ! on rajoute le calcul de zha pour diagnostiques (temp potentielle)
          zha(ig, l) = ztva(ig, l)
          ! if (l.lt.lentr(ig)) then
          ! ztva(ig,l) = ztva(ig,l)*(1.+RETV*(zqta(ig,l)
          ! s              -zqla(ig,l))-zqla(ig,l)) + 0.1
          ! else
          ztva(ig, l) = ztva(ig, l)*(1.+retv*(zqta(ig,l)-zqla(ig, &
            l))-zqla(ig,l))
          ! endif
          ! ztva(ig,l) = ztla(ig,l)*zpspsk(ig,l)+RLvCp*zqla(ig,l)
          ! s                 /(1.-retv*zqla(ig,l))
          ! ztva(ig,l) = ztva(ig,l)/zpspsk(ig,l)
          ! ztva(ig,l) = ztva(ig,l)*(1.+RETV*(zqta(ig,l)
          ! s                 /(1.-retv*zqta(ig,l))
          ! s              -zqla(ig,l)/(1.-retv*zqla(ig,l)))
          ! s              -zqla(ig,l)/(1.-retv*zqla(ig,l)))
          ! write(13,*)zqla(ig,l),zqla(ig,l)/(1.-retv*zqla(ig,l))
          ! on ecrit zqsat
          zqsatth(ig, l) = qsatbef(ig)
          ! enddo
          ! DO ig=1,ngrid
          ! if (zw2(ig,l).ge.1.e-10.and.
          ! s               f_star(ig,l)+entr_star(ig,l).gt.1.e-10) then
          ! mise a jour de la vitesse ascendante (l'air entraine de la couche
          ! consideree commence avec une vitesse nulle).

          ! if (f_star(ig,l+1).gt.1.e-10) then
          zw2(ig, l+1) = zw2(ig, l)* & ! s
                                       ! ((f_star(ig,l)-detr_star(ig,l))**2)
          ! s                  /f_star(ig,l+1)**2+
            ((f_star(ig,l))**2)/(f_star(ig,l+1)+detr_star(ig,l))**2 + & ! s
                                                                        ! /(f_star(ig,l+1))**2+
            2.*rg*(ztva(ig,l)-ztv(ig,l))/ztv(ig, l)*(zlev(ig,l+1)-zlev(ig,l))
          ! s                   *(f_star(ig,l)/f_star(ig,l+1))**2

        END IF
      END IF

      IF (zw2(ig,l+1)<0.) THEN
        linter(ig) = (l*(zw2(ig,l+1)-zw2(ig,l))-zw2(ig,l))/(zw2(ig,l+1)-zw2( &
          ig,l))
        zw2(ig, l+1) = 0.
        ! print*,'linter=',linter(ig)
        ! else if ((zw2(ig,l+1).lt.1.e-10).and.(zw2(ig,l+1).ge.0.)) then
        ! linter(ig)=l+1
        ! print*,'linter=l',zw2(ig,l),zw2(ig,l+1)
      ELSE
        wa_moy(ig, l+1) = sqrt(zw2(ig,l+1))
        ! wa_moy(ig,l+1)=zw2(ig,l+1)
      END IF
      IF (wa_moy(ig,l+1)>wmaxa(ig)) THEN
        ! lmix est le niveau de la couche ou w (wa_moy) est maximum
        lmix(ig) = l + 1
        wmaxa(ig) = wa_moy(ig, l+1)
      END IF
    END DO
  END DO
  PRINT *, 'fin calcul zw2'

  ! Calcul de la couche correspondant a la hauteur du thermique
  DO ig = 1, ngrid
    lmax(ig) = lentr(ig)
  END DO
  DO ig = 1, ngrid
    DO l = nlay, lentr(ig) + 1, -1
      IF (zw2(ig,l)<=1.E-10) THEN
        lmax(ig) = l - 1
      END IF
    END DO
  END DO
  ! pas de thermique si couche 1 stable
  DO ig = 1, ngrid
    IF (lmin(ig)>1) THEN
      lmax(ig) = 1
      lmin(ig) = 1
      lentr(ig) = 1
    END IF
  END DO

  ! Determination de zw2 max
  DO ig = 1, ngrid
    wmax(ig) = 0.
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (l<=lmax(ig)) THEN
        IF (zw2(ig,l)<0.) THEN
          PRINT *, 'pb2 zw2<0'
        END IF
        zw2(ig, l) = sqrt(zw2(ig,l))
        wmax(ig) = max(wmax(ig), zw2(ig,l))
      ELSE
        zw2(ig, l) = 0.
      END IF
    END DO
  END DO

  ! Longueur caracteristique correspondant a la hauteur des thermiques.
  DO ig = 1, ngrid
    zmax(ig) = 0.
    zlevinter(ig) = zlev(ig, 1)
  END DO
  DO ig = 1, ngrid
    ! calcul de zlevinter
    zlevinter(ig) = (zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*linter(ig) + &
      zlev(ig, lmax(ig)) - lmax(ig)*(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))
    ! pour le cas ou on prend tjs lmin=1
    ! zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
    zmax(ig) = max(zmax(ig), zlevinter(ig)-zlev(ig,1))
    zmax0(ig) = zmax(ig)
    WRITE (11, *) 'ig,lmax,linter', ig, lmax(ig), linter(ig)
    WRITE (12, *) 'ig,zlevinter,zmax', ig, zmax(ig), zlevinter(ig)
  END DO

  ! Calcul de zmax_sec et wmax_sec
  CALL fermeture_seche(ngrid, nlay, pplay, pplev, pphi, zlev, rhobarz, f0, &
    zpspsk, alim, zh, zo, lentr, lmin, nu_min, nu_max, r_aspect, zmax_sec2, &
    wmax_sec2)

  PRINT *, 'avant fermeture'
  ! Fermeture,determination de f
  ! en lmax f=d-e
  DO ig = 1, ngrid
    ! entr_star(ig,lmax(ig))=0.
    ! f_star(ig,lmax(ig)+1)=0.
    ! detr_star(ig,lmax(ig))=f_star(ig,lmax(ig))+entr_star(ig,lmax(ig))
    ! s                       +alim_star(ig,lmax(ig))
  END DO

  DO ig = 1, ngrid
    alim_star2(ig) = 0.
  END DO
  ! calcul de entr_star_tot
  DO ig = 1, ngrid
    DO k = 1, lmix(ig)
      entr_star_tot(ig) = entr_star_tot(ig) & ! s
                                              ! +entr_star(ig,k)
        +alim_star(ig, k)
      ! s                        -detr_star(ig,k)
      detr_star_tot(ig) = detr_star_tot(ig) & ! s
                                              ! +alim_star(ig,k)
        -detr_star(ig, k) + entr_star(ig, k)
    END DO
  END DO

  DO ig = 1, ngrid
    IF (alim_star_tot(ig)<1.E-10) THEN
      f(ig) = 0.
    ELSE
      ! do k=lmin(ig),lentr(ig)
      DO k = 1, lentr(ig)
        alim_star2(ig) = alim_star2(ig) + alim_star(ig, k)**2/(rho(ig,k)*( &
          zlev(ig,k+1)-zlev(ig,k)))
      END DO
      IF ((zmax_sec(ig)>1.E-10) .AND. (1==1)) THEN
        f(ig) = wmax_sec(ig)/(max(500.,zmax_sec(ig))*r_aspect*alim_star2(ig))
        f(ig) = f(ig) + (f0(ig)-f(ig))*exp((-ptimestep/zmax_sec(ig))*wmax_sec &
          (ig))
      ELSE
        f(ig) = wmax(ig)/(max(500.,zmax(ig))*r_aspect*alim_star2(ig))
        f(ig) = f(ig) + (f0(ig)-f(ig))*exp((-ptimestep/zmax(ig))*wmax(ig))
      END IF
    END IF
    f0(ig) = f(ig)
  END DO
  PRINT *, 'apres fermeture'
  ! Calcul de l'entrainement
  DO ig = 1, ngrid
    DO k = 1, klev
      alim(ig, k) = f(ig)*alim_star(ig, k)
    END DO
  END DO
  ! CR:test pour entrainer moins que la masse
  ! do ig=1,ngrid
  ! do l=1,lentr(ig)
  ! if ((alim(ig,l)*ptimestep).gt.(0.9*masse(ig,l))) then
  ! alim(ig,l+1)=alim(ig,l+1)+alim(ig,l)
  ! s                       -0.9*masse(ig,l)/ptimestep
  ! alim(ig,l)=0.9*masse(ig,l)/ptimestep
  ! endif
  ! enddo
  ! enddo
  ! calcul du détrainement
  DO ig = 1, klon
    DO k = 1, klev
      detr(ig, k) = f(ig)*detr_star(ig, k)
      IF (detr(ig,k)<0.) THEN
        ! print*,'detr1<0!!!'
      END IF
    END DO
    DO k = 1, klev
      entr(ig, k) = f(ig)*entr_star(ig, k)
      IF (entr(ig,k)<0.) THEN
        ! print*,'entr1<0!!!'
      END IF
    END DO
  END DO

  ! do ig=1,ngrid
  ! do l=1,klev
  ! if (((detr(ig,l)+entr(ig,l)+alim(ig,l))*ptimestep).gt.
  ! s          (masse(ig,l))) then
  ! print*,'d2+e2+a2>m2','ig=',ig,'l=',l,'lmax(ig)=',lmax(ig),'d+e+a='
  ! s,(detr(ig,l)+entr(ig,l)+alim(ig,l))*ptimestep,'m=',masse(ig,l)
  ! endif
  ! enddo
  ! enddo
  ! Calcul des flux

  DO ig = 1, ngrid
    DO l = 1, lmax(ig)
      ! do l=1,klev
      ! fmc(ig,l+1)=f(ig)*f_star(ig,l+1)
      fmc(ig, l+1) = fmc(ig, l) + alim(ig, l) + entr(ig, l) - detr(ig, l)
      ! print*,'??!!','ig=',ig,'l=',l,'lmax=',lmax(ig),'lmix=',lmix(ig),
      ! s  'e=',entr(ig,l),'d=',detr(ig,l),'a=',alim(ig,l),'f=',fmc(ig,l),
      ! s  'f+1=',fmc(ig,l+1)
      IF (fmc(ig,l+1)<0.) THEN
        PRINT *, 'fmc1<0', l + 1, lmax(ig), fmc(ig, l+1)
        fmc(ig, l+1) = fmc(ig, l)
        detr(ig, l) = alim(ig, l) + entr(ig, l)
        ! fmc(ig,l+1)=0.
        ! print*,'fmc1<0',l+1,lmax(ig),fmc(ig,l+1)
      END IF
      ! if ((fmc(ig,l+1).gt.fmc(ig,l)).and.(l.gt.lentr(ig))) then
      ! f_old=fmc(ig,l+1)
      ! fmc(ig,l+1)=fmc(ig,l)
      ! detr(ig,l)=detr(ig,l)+f_old-fmc(ig,l+1)
      ! endif

      ! if ((fmc(ig,l+1).gt.fmc(ig,l)).and.(l.gt.lentr(ig))) then
      ! f_old=fmc(ig,l+1)
      ! fmc(ig,l+1)=fmc(ig,l)
      ! detr(ig,l)=detr(ig,l)+f_old-fmc(ig,l)
      ! endif
      ! rajout du test sur alpha croissant
      ! if test
      ! if (1.eq.0) then

      IF (l==klev) THEN
        PRINT *, 'THERMCELL PB ig=', ig, '   l=', l
        abort_message = 'THERMCELL PB'
        CALL abort_physic(modname, abort_message, 1)
      END IF
      ! if ((zw2(ig,l+1).gt.1.e-10).and.(zw2(ig,l).gt.1.e-10).and.
      ! s     (l.ge.lentr(ig)).and.
      IF ((zw2(ig,l+1)>1.E-10) .AND. (zw2(ig,l)>1.E-10) .AND. (l>=lentr(ig))) &
          THEN
        IF (((fmc(ig,l+1)/(rhobarz(ig,l+1)*zw2(ig,l+1)))>(fmc(ig,l)/ &
            (rhobarz(ig,l)*zw2(ig,l))))) THEN
          f_old = fmc(ig, l+1)
          fmc(ig, l+1) = fmc(ig, l)*rhobarz(ig, l+1)*zw2(ig, l+1)/ &
            (rhobarz(ig,l)*zw2(ig,l))
          detr(ig, l) = detr(ig, l) + f_old - fmc(ig, l+1)
          ! detr(ig,l)=(fmc(ig,l+1)-fmc(ig,l))/(0.4-1.)
          ! entr(ig,l)=0.4*detr(ig,l)
          ! entr(ig,l)=fmc(ig,l+1)-fmc(ig,l)+detr(ig,l)
        END IF
      END IF
      IF ((fmc(ig,l+1)>fmc(ig,l)) .AND. (l>lentr(ig))) THEN
        f_old = fmc(ig, l+1)
        fmc(ig, l+1) = fmc(ig, l)
        detr(ig, l) = detr(ig, l) + f_old - fmc(ig, l+1)
      END IF
      IF (detr(ig,l)>fmc(ig,l)) THEN
        detr(ig, l) = fmc(ig, l)
        entr(ig, l) = fmc(ig, l+1) - alim(ig, l)
      END IF
      IF (fmc(ig,l+1)<0.) THEN
        detr(ig, l) = detr(ig, l) + fmc(ig, l+1)
        fmc(ig, l+1) = 0.
        PRINT *, 'fmc2<0', l + 1, lmax(ig)
      END IF

      ! test pour ne pas avoir f=0 et d=e/=0
      ! if (fmc(ig,l+1).lt.1.e-10) then
      ! detr(ig,l+1)=0.
      ! entr(ig,l+1)=0.
      ! zqla(ig,l+1)=0.
      ! zw2(ig,l+1)=0.
      ! lmax(ig)=l+1
      ! zmax(ig)=zlev(ig,lmax(ig))
      ! endif
      IF (zw2(ig,l+1)>1.E-10) THEN
        IF ((((fmc(ig,l+1))/(rhobarz(ig,l+1)*zw2(ig,l+1)))>1.)) THEN
          f_old = fmc(ig, l+1)
          fmc(ig, l+1) = rhobarz(ig, l+1)*zw2(ig, l+1)
          zw2(ig, l+1) = 0.
          zqla(ig, l+1) = 0.
          detr(ig, l) = detr(ig, l) + f_old - fmc(ig, l+1)
          lmax(ig) = l + 1
          zmax(ig) = zlev(ig, lmax(ig))
          PRINT *, 'alpha>1', l + 1, lmax(ig)
        END IF
      END IF
      ! write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
      ! endif test
      ! endif
    END DO
  END DO
  DO ig = 1, ngrid
    ! if (fmc(ig,lmax(ig)+1).ne.0.) then
    fmc(ig, lmax(ig)+1) = 0.
    entr(ig, lmax(ig)) = 0.
    detr(ig, lmax(ig)) = fmc(ig, lmax(ig)) + entr(ig, lmax(ig)) + &
      alim(ig, lmax(ig))
    ! endif
  END DO
  ! test sur le signe de fmc
  DO ig = 1, ngrid
    DO l = 1, klev + 1
      IF (fmc(ig,l)<0.) THEN
        PRINT *, 'fm1<0!!!', 'ig=', ig, 'l=', l, 'a=', alim(ig, l-1), 'e=', &
          entr(ig, l-1), 'f=', fmc(ig, l-1), 'd=', detr(ig, l-1), 'f+1=', &
          fmc(ig, l)
      END IF
    END DO
  END DO
  ! test de verification
  DO ig = 1, ngrid
    DO l = 1, lmax(ig)
      IF ((abs(fmc(ig,l+1)-fmc(ig,l)-alim(ig,l)-entr(ig,l)+ &
          detr(ig,l)))>1.E-4) THEN
        ! print*,'pbcm!!','ig=',ig,'l=',l,'lmax=',lmax(ig),'lmix=',lmix(ig),
        ! s  'e=',entr(ig,l),'d=',detr(ig,l),'a=',alim(ig,l),'f=',fmc(ig,l),
        ! s  'f+1=',fmc(ig,l+1)
      END IF
      IF (detr(ig,l)<0.) THEN
        PRINT *, 'detrdemi<0!!!'
      END IF
    END DO
  END DO

  ! RC
  ! CR def de  zmix continu (profil parabolique des vitesses)
  DO ig = 1, ngrid
    IF (lmix(ig)>1.) THEN
      ! test
      IF (((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig)))- &
          (zlev(ig,lmix(ig)+1)))-(zw2(ig,lmix(ig))- &
          zw2(ig,lmix(ig)+1))*((zlev(ig,lmix(ig)-1))- &
          (zlev(ig,lmix(ig)))))>1E-10) THEN

        zmix(ig) = ((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig)) &
          )**2-(zlev(ig,lmix(ig)+1))**2)-(zw2(ig,lmix(ig))-zw2(ig, &
          lmix(ig)+1))*((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))/ &
          (2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig)))- &
          (zlev(ig,lmix(ig)+1)))-(zw2(ig,lmix(ig))- &
          zw2(ig,lmix(ig)+1))*((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
      ELSE
        zmix(ig) = zlev(ig, lmix(ig))
        PRINT *, 'pb zmix'
      END IF
    ELSE
      zmix(ig) = 0.
    END IF
    ! test
    IF ((zmax(ig)-zmix(ig))<=0.) THEN
      zmix(ig) = 0.9*zmax(ig)
      ! print*,'pb zmix>zmax'
    END IF
  END DO
  DO ig = 1, klon
    zmix0(ig) = zmix(ig)
  END DO

  ! calcul du nouveau lmix correspondant
  DO ig = 1, ngrid
    DO l = 1, klev
      IF (zmix(ig)>=zlev(ig,l) .AND. zmix(ig)<zlev(ig,l+1)) THEN
        lmix(ig) = l
      END IF
    END DO
  END DO

  ! ne devrait pas arriver!!!!!
  DO ig = 1, ngrid
    DO l = 1, klev
      IF (detr(ig,l)>(fmc(ig,l)+alim(ig,l))+entr(ig,l)) THEN
        PRINT *, 'detr2>fmc2!!!', 'ig=', ig, 'l=', l, 'd=', detr(ig, l), &
          'f=', fmc(ig, l), 'lmax=', lmax(ig)
        ! detr(ig,l)=fmc(ig,l)+alim(ig,l)+entr(ig,l)
        ! entr(ig,l)=0.
        ! fmc(ig,l+1)=0.
        ! zw2(ig,l+1)=0.
        ! zqla(ig,l+1)=0.
        PRINT *, 'pb!fm=0 et f_star>0', l, lmax(ig)
        ! lmax(ig)=l
      END IF
    END DO
  END DO
  DO ig = 1, ngrid
    DO l = lmax(ig) + 1, klev + 1
      ! fmc(ig,l)=0.
      ! detr(ig,l)=0.
      ! entr(ig,l)=0.
      ! zw2(ig,l)=0.
      ! zqla(ig,l)=0.
    END DO
  END DO

  ! Calcul du detrainement lors du premier passage
  ! print*,'9 OK convect8'
  ! print*,'WA1 ',wa_moy

  ! determination de l'indice du debut de la mixed layer ou w decroit

  ! calcul de la largeur de chaque ascendance dans le cas conservatif.
  ! dans ce cas simple, on suppose que la largeur de l'ascendance provenant
  ! d'une couche est égale à la hauteur de la couche alimentante.
  ! La vitesse maximale dans l'ascendance est aussi prise comme estimation
  ! de la vitesse d'entrainement horizontal dans la couche alimentante.

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (l<=lmax(ig) .AND. (test(ig)==1)) THEN
        zw = max(wa_moy(ig,l), 1.E-10)
        larg_cons(ig, l) = zmax(ig)*r_aspect*fmc(ig, l)/(rhobarz(ig,l)*zw)
      END IF
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (l<=lmax(ig) .AND. (test(ig)==1)) THEN
        ! if (idetr.eq.0) then
        ! cette option est finalement en dur.
        IF ((l_mix*zlev(ig,l))<0.) THEN
          PRINT *, 'pb l_mix*zlev<0'
        END IF
        ! CR: test: nouvelle def de lambda
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        IF (zw2(ig,l)>1.E-10) THEN
          larg_detr(ig, l) = sqrt((l_mix/zw2(ig,l))*zlev(ig,l))
        ELSE
          larg_detr(ig, l) = sqrt(l_mix*zlev(ig,l))
        END IF
        ! else if (idetr.eq.1) then
        ! larg_detr(ig,l)=larg_cons(ig,l)
        ! s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
        ! else if (idetr.eq.2) then
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        ! s            *sqrt(wa_moy(ig,l))
        ! else if (idetr.eq.4) then
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        ! s            *wa_moy(ig,l)
        ! endif
      END IF
    END DO
  END DO

  ! print*,'10 OK convect8'
  ! print*,'WA2 ',wa_moy
  ! cal1cul de la fraction de la maille concernée par l'ascendance en tenant
  ! compte de l'epluchage du thermique.


  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (larg_cons(ig,l)>1. .AND. (test(ig)==1)) THEN
        ! print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
        fraca(ig, l) = (larg_cons(ig,l)-larg_detr(ig,l))/(r_aspect*zmax(ig))
        ! test
        fraca(ig, l) = max(fraca(ig,l), 0.)
        fraca(ig, l) = min(fraca(ig,l), 0.5)
        fracd(ig, l) = 1. - fraca(ig, l)
        fracc(ig, l) = larg_cons(ig, l)/(r_aspect*zmax(ig))
      ELSE
        ! wa_moy(ig,l)=0.
        fraca(ig, l) = 0.
        fracc(ig, l) = 0.
        fracd(ig, l) = 1.
      END IF
    END DO
  END DO
  ! CR: calcul de fracazmix
  DO ig = 1, ngrid
    IF (test(ig)==1) THEN
      fracazmix(ig) = (fraca(ig,lmix(ig)+1)-fraca(ig,lmix(ig)))/ &
        (zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))*zmix(ig) + &
        fraca(ig, lmix(ig)) - zlev(ig, lmix(ig))*(fraca(ig,lmix(ig)+1)-fraca( &
        ig,lmix(ig)))/(zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))
    END IF
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (larg_cons(ig,l)>1. .AND. (test(ig)==1)) THEN
        IF (l>lmix(ig)) THEN
          ! test
          IF (zmax(ig)-zmix(ig)<1.E-10) THEN
            ! print*,'pb xxx'
            xxx(ig, l) = (lmax(ig)+1.-l)/(lmax(ig)+1.-lmix(ig))
          ELSE
            xxx(ig, l) = (zmax(ig)-zlev(ig,l))/(zmax(ig)-zmix(ig))
          END IF
          IF (idetr==0) THEN
            fraca(ig, l) = fracazmix(ig)
          ELSE IF (idetr==1) THEN
            fraca(ig, l) = fracazmix(ig)*xxx(ig, l)
          ELSE IF (idetr==2) THEN
            fraca(ig, l) = fracazmix(ig)*(1.-(1.-xxx(ig,l))**2)
          ELSE
            fraca(ig, l) = fracazmix(ig)*xxx(ig, l)**2
          END IF
          ! print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
          fraca(ig, l) = max(fraca(ig,l), 0.)
          fraca(ig, l) = min(fraca(ig,l), 0.5)
          fracd(ig, l) = 1. - fraca(ig, l)
          fracc(ig, l) = larg_cons(ig, l)/(r_aspect*zmax(ig))
        END IF
      END IF
    END DO
  END DO

  PRINT *, 'fin calcul fraca'
  ! print*,'11 OK convect8'
  ! print*,'Ea3 ',wa_moy
  ! ------------------------------------------------------------------
  ! Calcul de fracd, wd
  ! somme wa - wd = 0
  ! ------------------------------------------------------------------


  DO ig = 1, ngrid
    fm(ig, 1) = 0.
    fm(ig, nlay+1) = 0.
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (test(ig)==1) THEN
        fm(ig, l) = fraca(ig, l)*wa_moy(ig, l)*rhobarz(ig, l)
        ! CR:test
        IF (alim(ig,l-1)<1E-10 .AND. fm(ig,l)>fm(ig,l-1) .AND. l>lmix(ig)) &
            THEN
          fm(ig, l) = fm(ig, l-1)
          ! write(1,*)'ajustement fm, l',l
        END IF
        ! write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
        ! RC
      END IF
    END DO
    DO ig = 1, ngrid
      IF (fracd(ig,l)<0.1 .AND. (test(ig)==1)) THEN
        abort_message = 'fracd trop petit'
        CALL abort_physic(modname, abort_message, 1)
      ELSE
        ! vitesse descendante "diagnostique"
        wd(ig, l) = fm(ig, l)/(fracd(ig,l)*rhobarz(ig,l))
      END IF
    END DO
  END DO

  DO l = 1, nlay + 1
    DO ig = 1, ngrid
      IF (test(ig)==0) THEN
        fm(ig, l) = fmc(ig, l)
      END IF
    END DO
  END DO

  ! fin du first
  DO l = 1, nlay
    DO ig = 1, ngrid
      ! masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
      masse(ig, l) = (pplev(ig,l)-pplev(ig,l+1))/rg
    END DO
  END DO

  ! print*,'12 OK convect8'
  ! print*,'WA4 ',wa_moy
  ! c------------------------------------------------------------------
  ! calcul du transport vertical
  ! ------------------------------------------------------------------

  GO TO 4444
  ! print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
  DO l = 2, nlay - 1
    DO ig = 1, ngrid
      IF (fm(ig,l+1)*ptimestep>masse(ig,l) .AND. fm(ig,l+1)*ptimestep>masse( &
          ig,l+1)) THEN
        PRINT *, 'WARN!!! FM>M ig=', ig, ' l=', l, '  FM=', &
          fm(ig, l+1)*ptimestep, '   M=', masse(ig, l), masse(ig, l+1)
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF ((alim(ig,l)+entr(ig,l))*ptimestep>masse(ig,l)) THEN
        PRINT *, 'WARN!!! E>M ig=', ig, ' l=', l, '  E==', &
          (entr(ig,l)+alim(ig,l))*ptimestep, '   M=', masse(ig, l)
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (.NOT. fm(ig,l)>=0. .OR. .NOT. fm(ig,l)<=10.) THEN
        ! print*,'WARN!!! fm exagere ig=',ig,'   l=',l
        ! s         ,'   FM=',fm(ig,l)
      END IF
      IF (.NOT. masse(ig,l)>=1.E-10 .OR. .NOT. masse(ig,l)<=1.E4) THEN
        ! print*,'WARN!!! masse exagere ig=',ig,'   l=',l
        ! s         ,'   M=',masse(ig,l)
        ! print*,'rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)',
        ! s                 rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)
        ! print*,'zlev(ig,l+1),zlev(ig,l)'
        ! s                ,zlev(ig,l+1),zlev(ig,l)
        ! print*,'pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)'
        ! s                ,pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)
      END IF
      IF (.NOT. alim(ig,l)>=0. .OR. .NOT. alim(ig,l)<=10.) THEN
        ! print*,'WARN!!! entr exagere ig=',ig,'   l=',l
        ! s         ,'   E=',entr(ig,l)
      END IF
    END DO
  END DO

4444 CONTINUE

  ! CR:redefinition du entr
  ! CR:test:on ne change pas la def du entr mais la def du fm
  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (test(ig)==1) THEN
        detr(ig, l) = fm(ig, l) + alim(ig, l) - fm(ig, l+1)
        IF (detr(ig,l)<0.) THEN
          ! entr(ig,l)=entr(ig,l)-detr(ig,l)
          fm(ig, l+1) = fm(ig, l) + alim(ig, l)
          detr(ig, l) = 0.
          ! write(11,*)'l,ig,entr',l,ig,entr(ig,l)
          ! print*,'WARNING !!! detrainement negatif ',ig,l
        END IF
      END IF
    END DO
  END DO
  ! RC

  IF (w2di==1) THEN
    fm0 = fm0 + ptimestep*(fm-fm0)/tho
    entr0 = entr0 + ptimestep*(alim+entr-entr0)/tho
  ELSE
    fm0 = fm
    entr0 = alim + entr
    detr0 = detr
    alim0 = alim
    ! zoa=zqta
    ! entr0=alim
  END IF

  IF (1==1) THEN
    ! call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
    ! .    ,zh,zdhadj,zha)
    ! call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
    ! .    ,zo,pdoadj,zoa)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zthl, &
      zdthladj, zta)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, po, pdoadj, &
      zoa)
  ELSE
    CALL dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zh, &
      zdhadj, zha)
    CALL dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zo, &
      pdoadj, zoa)
  END IF

  IF (1==0) THEN
    CALL dvthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zmax, &
      zu, zv, pduadj, pdvadj, zua, zva)
  ELSE
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zu, pduadj, &
      zua)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zv, pdvadj, &
      zva)
  END IF

  ! Calcul des moments
  ! do l=1,nlay
  ! do ig=1,ngrid
  ! zf=0.5*(fracc(ig,l)+fracc(ig,l+1))
  ! zf2=zf/(1.-zf)
  ! thetath2(ig,l)=zf2*(zha(ig,l)-zh(ig,l))**2
  ! wth2(ig,l)=zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
  ! enddo
  ! enddo






  ! print*,'13 OK convect8'
  ! print*,'WA5 ',wa_moy
  DO l = 1, nlay
    DO ig = 1, ngrid
      ! pdtadj(ig,l)=zdhadj(ig,l)*zpspsk(ig,l)
      pdtadj(ig, l) = zdthladj(ig, l)*zpspsk(ig, l)
    END DO
  END DO


  ! do l=1,nlay
  ! do ig=1,ngrid
  ! if(abs(pdtadj(ig,l))*86400..gt.500.) then
  ! print*,'WARN!!! ig=',ig,'  l=',l
  ! s         ,'   pdtadj=',pdtadj(ig,l)
  ! endif
  ! if(abs(pdoadj(ig,l))*86400..gt.1.) then
  ! print*,'WARN!!! ig=',ig,'  l=',l
  ! s         ,'   pdoadj=',pdoadj(ig,l)
  ! endif
  ! enddo
  ! enddo

  ! print*,'14 OK convect8'
  ! ------------------------------------------------------------------
  ! Calculs pour les sorties
  ! ------------------------------------------------------------------
  ! calcul de fraca pour les sorties
  DO l = 2, klev
    DO ig = 1, klon
      IF (zw2(ig,l)>1.E-10) THEN
        fraca(ig, l) = fm(ig, l)/(rhobarz(ig,l)*zw2(ig,l))
      ELSE
        fraca(ig, l) = 0.
      END IF
    END DO
  END DO
  IF (sorties) THEN
    DO l = 1, nlay
      DO ig = 1, ngrid
        zla(ig, l) = (1.-fracd(ig,l))*zmax(ig)
        zld(ig, l) = fracd(ig, l)*zmax(ig)
        IF (1.-fracd(ig,l)>1.E-10) zwa(ig, l) = wd(ig, l)*fracd(ig, l)/ &
          (1.-fracd(ig,l))
      END DO
    END DO
    ! CR calcul du niveau de condensation
    ! initialisation
    DO ig = 1, ngrid
      nivcon(ig) = 0.
      zcon(ig) = 0.
    END DO
    DO k = nlay, 1, -1
      DO ig = 1, ngrid
        IF (zqla(ig,k)>1E-10) THEN
          nivcon(ig) = k
          zcon(ig) = zlev(ig, k)
        END IF
        ! if (zcon(ig).gt.1.e-10) then
        ! nuage=.true.
        ! else
        ! nuage=.false.
        ! endif
      END DO
    END DO

    DO l = 1, nlay
      DO ig = 1, ngrid
        zf = fraca(ig, l)
        zf2 = zf/(1.-zf)
        thetath2(ig, l) = zf2*(zha(ig,l)-zh(ig,l)/zpspsk(ig,l))**2
        wth2(ig, l) = zf2*(zw2(ig,l))**2
        ! print*,'wth2=',wth2(ig,l)
        wth3(ig, l) = zf2*(1-2.*fraca(ig,l))/(1-fraca(ig,l))*zw2(ig, l)* &
          zw2(ig, l)*zw2(ig, l)
        q2(ig, l) = zf2*(zqta(ig,l)*1000.-po(ig,l)*1000.)**2
        ! test: on calcul q2/po=ratqsc
        ! if (nuage) then
        ratqscth(ig, l) = sqrt(q2(ig,l))/(po(ig,l)*1000.)
        ! else
        ! ratqscth(ig,l)=0.
        ! endif
      END DO
    END DO
    ! calcul du ratqscdiff
    sum = 0.
    sumdiff = 0.
    ratqsdiff(:, :) = 0.
    DO ig = 1, ngrid
      DO l = 1, lentr(ig)
        sum = sum + alim_star(ig, l)*zqta(ig, l)*1000.
      END DO
    END DO
    DO ig = 1, ngrid
      DO l = 1, lentr(ig)
        zf = fraca(ig, l)
        zf2 = zf/(1.-zf)
        sumdiff = sumdiff + alim_star(ig, l)*(zqta(ig,l)*1000.-sum)**2
        ! ratqsdiff=ratqsdiff+alim_star(ig,l)*
        ! s          (zqta(ig,l)*1000.-po(ig,l)*1000.)**2
      END DO
    END DO
    DO l = 1, klev
      DO ig = 1, ngrid
        ratqsdiff(ig, l) = sqrt(sumdiff)/(po(ig,l)*1000.)
        ! write(11,*)'ratqsdiff=',ratqsdiff(ig,l)
      END DO
    END DO

  END IF

  ! print*,'19 OK convect8'
  RETURN
END SUBROUTINE thermcell_cld

SUBROUTINE thermcell_eau(ngrid, nlay, ptimestep, pplay, pplev, pphi, pu, pv, &
    pt, po, pduadj, pdvadj, pdtadj, pdoadj, fm0, entr0 & ! s
                                                         ! ,pu_therm,pv_therm
    , r_aspect, l_mix, w2di, tho)

  USE dimphy
  IMPLICIT NONE

  ! =======================================================================

  ! Calcul du transport verticale dans la couche limite en presence
  ! de "thermiques" explicitement representes

  ! Réécriture à partir d'un listing papier à Habas, le 14/02/00

  ! le thermique est supposé homogène et dissipé par mélange avec
  ! son environnement. la longueur l_mix contrôle l'efficacité du
  ! mélange

  ! Le calcul du transport des différentes espèces se fait en prenant
  ! en compte:
  ! 1. un flux de masse montant
  ! 2. un flux de masse descendant
  ! 3. un entrainement
  ! 4. un detrainement

  ! =======================================================================

  ! -----------------------------------------------------------------------
  ! declarations:
  ! -------------

  include "YOMCST.h"
  include "YOETHF.h"
  include "FCTTRE.h"

  ! arguments:
  ! ----------

  INTEGER ngrid, nlay, w2di
  REAL tho
  REAL ptimestep, l_mix, r_aspect
  REAL pt(ngrid, nlay), pdtadj(ngrid, nlay)
  REAL pu(ngrid, nlay), pduadj(ngrid, nlay)
  REAL pv(ngrid, nlay), pdvadj(ngrid, nlay)
  REAL po(ngrid, nlay), pdoadj(ngrid, nlay)
  REAL pplay(ngrid, nlay), pplev(ngrid, nlay+1)
  REAL pphi(ngrid, nlay)

  INTEGER idetr
  SAVE idetr
  DATA idetr/3/
  !$OMP THREADPRIVATE(idetr)

  ! local:
  ! ------

  INTEGER ig, k, l, lmaxa(klon), lmix(klon)
  REAL zsortie1d(klon)
  ! CR: on remplace lmax(klon,klev+1)
  INTEGER lmax(klon), lmin(klon), lentr(klon)
  REAL linter(klon)
  REAL zmix(klon), fracazmix(klon)
  ! RC
  REAL zmax(klon), zw, zz, zw2(klon, klev+1), ztva(klon, klev), zzz

  REAL zlev(klon, klev+1), zlay(klon, klev)
  REAL zh(klon, klev), zdhadj(klon, klev)
  REAL zthl(klon, klev), zdthladj(klon, klev)
  REAL ztv(klon, klev)
  REAL zu(klon, klev), zv(klon, klev), zo(klon, klev)
  REAL zl(klon, klev)
  REAL wh(klon, klev+1)
  REAL wu(klon, klev+1), wv(klon, klev+1), wo(klon, klev+1)
  REAL zla(klon, klev+1)
  REAL zwa(klon, klev+1)
  REAL zld(klon, klev+1)
  REAL zwd(klon, klev+1)
  REAL zsortie(klon, klev)
  REAL zva(klon, klev)
  REAL zua(klon, klev)
  REAL zoa(klon, klev)

  REAL zta(klon, klev)
  REAL zha(klon, klev)
  REAL wa_moy(klon, klev+1)
  REAL fraca(klon, klev+1)
  REAL fracc(klon, klev+1)
  REAL zf, zf2
  REAL thetath2(klon, klev), wth2(klon, klev)
  ! common/comtherm/thetath2,wth2

  REAL count_time
  INTEGER ialt

  LOGICAL sorties
  REAL rho(klon, klev), rhobarz(klon, klev+1), masse(klon, klev)
  REAL zpspsk(klon, klev)

  ! real wmax(klon,klev),wmaxa(klon)
  REAL wmax(klon), wmaxa(klon)
  REAL wa(klon, klev, klev+1)
  REAL wd(klon, klev+1)
  REAL larg_part(klon, klev, klev+1)
  REAL fracd(klon, klev+1)
  REAL xxx(klon, klev+1)
  REAL larg_cons(klon, klev+1)
  REAL larg_detr(klon, klev+1)
  REAL fm0(klon, klev+1), entr0(klon, klev), detr(klon, klev)
  REAL pu_therm(klon, klev), pv_therm(klon, klev)
  REAL fm(klon, klev+1), entr(klon, klev)
  REAL fmc(klon, klev+1)

  REAL zcor, zdelta, zcvm5, qlbef
  REAL tbef(klon), qsatbef(klon)
  REAL dqsat_dt, dt, num, denom
  REAL reps, rlvcp, ddt0
  REAL ztla(klon, klev), zqla(klon, klev), zqta(klon, klev)

  PARAMETER (ddt0=.01)

  ! CR:nouvelles variables
  REAL f_star(klon, klev+1), entr_star(klon, klev)
  REAL entr_star_tot(klon), entr_star2(klon)
  REAL f(klon), f0(klon)
  REAL zlevinter(klon)
  LOGICAL first
  DATA first/.FALSE./
  SAVE first
  !$OMP THREADPRIVATE(first)

  ! RC

  CHARACTER *2 str2
  CHARACTER *10 str10

  CHARACTER (LEN=20) :: modname = 'thermcell_eau'
  CHARACTER (LEN=80) :: abort_message

  LOGICAL vtest(klon), down
  LOGICAL zsat(klon)

  EXTERNAL scopy

  INTEGER ncorrec, ll
  SAVE ncorrec
  DATA ncorrec/0/
  !$OMP THREADPRIVATE(ncorrec)



  ! -----------------------------------------------------------------------
  ! initialisation:
  ! ---------------

  sorties = .TRUE.
  IF (ngrid/=klon) THEN
    PRINT *
    PRINT *, 'STOP dans convadj'
    PRINT *, 'ngrid    =', ngrid
    PRINT *, 'klon  =', klon
  END IF

  ! Initialisation
  rlvcp = rlvtt/rcpd
  reps = rd/rv

  ! -----------------------------------------------------------------------
  ! AM Calcul de T,q,ql a partir de Tl et qT
  ! ---------------------------------------------------

  ! Pr Tprec=Tl calcul de qsat
  ! Si qsat>qT T=Tl, q=qT
  ! Sinon DDT=(-Tprec+Tl+RLVCP (qT-qsat(T')) / (1+RLVCP dqsat/dt)
  ! On cherche DDT < DDT0

  ! defaut
  DO ll = 1, nlay
    DO ig = 1, ngrid
      zo(ig, ll) = po(ig, ll)
      zl(ig, ll) = 0.
      zh(ig, ll) = pt(ig, ll)
    END DO
  END DO
  DO ig = 1, ngrid
    zsat(ig) = .FALSE.
  END DO


  DO ll = 1, nlay
    ! les points insatures sont definitifs
    DO ig = 1, ngrid
      tbef(ig) = pt(ig, ll)
      zdelta = max(0., sign(1.,rtt-tbef(ig)))
      qsatbef(ig) = r2es*foeew(tbef(ig), zdelta)/pplev(ig, ll)
      qsatbef(ig) = min(0.5, qsatbef(ig))
      zcor = 1./(1.-retv*qsatbef(ig))
      qsatbef(ig) = qsatbef(ig)*zcor
      zsat(ig) = (max(0.,po(ig,ll)-qsatbef(ig))>0.00001)
    END DO

    DO ig = 1, ngrid
      IF (zsat(ig)) THEN
        qlbef = max(0., po(ig,ll)-qsatbef(ig))
        ! si sature: ql est surestime, d'ou la sous-relax
        dt = 0.5*rlvcp*qlbef
        ! on pourra enchainer 2 ou 3 calculs sans Do while
        DO WHILE (dt>ddt0)
          ! il faut verifier si c,a conserve quand on repasse en insature ...
          tbef(ig) = tbef(ig) + dt
          zdelta = max(0., sign(1.,rtt-tbef(ig)))
          qsatbef(ig) = r2es*foeew(tbef(ig), zdelta)/pplev(ig, ll)
          qsatbef(ig) = min(0.5, qsatbef(ig))
          zcor = 1./(1.-retv*qsatbef(ig))
          qsatbef(ig) = qsatbef(ig)*zcor
          ! on veut le signe de qlbef
          qlbef = po(ig, ll) - qsatbef(ig)
          ! dqsat_dT
          zdelta = max(0., sign(1.,rtt-tbef(ig)))
          zcvm5 = r5les*(1.-zdelta) + r5ies*zdelta
          zcor = 1./(1.-retv*qsatbef(ig))
          dqsat_dt = foede(tbef(ig), zdelta, zcvm5, qsatbef(ig), zcor)
          num = -tbef(ig) + pt(ig, ll) + rlvcp*qlbef
          denom = 1. + rlvcp*dqsat_dt
          dt = num/denom
        END DO
        ! on ecrit de maniere conservative (sat ou non)
        zl(ig, ll) = max(0., qlbef)
        ! T = Tl +Lv/Cp ql
        zh(ig, ll) = pt(ig, ll) + rlvcp*zl(ig, ll)
        zo(ig, ll) = po(ig, ll) - zl(ig, ll)
      END IF
    END DO
  END DO
  ! AM fin

  ! -----------------------------------------------------------------------
  ! incrementation eventuelle de tendances precedentes:
  ! ---------------------------------------------------

  ! print*,'0 OK convect8'

  DO l = 1, nlay
    DO ig = 1, ngrid
      zpspsk(ig, l) = (pplay(ig,l)/pplev(ig,1))**rkappa
      ! zh(ig,l)=pt(ig,l)/zpspsk(ig,l)
      zu(ig, l) = pu(ig, l)
      zv(ig, l) = pv(ig, l)
      ! zo(ig,l)=po(ig,l)
      ! ztv(ig,l)=zh(ig,l)*(1.+0.61*zo(ig,l))
      ! AM attention zh est maintenant le profil de T et plus le profil de
      ! theta !

      ! T-> Theta
      ztv(ig, l) = zh(ig, l)/zpspsk(ig, l)
      ! AM Theta_v
      ztv(ig, l) = ztv(ig, l)*(1.+retv*(zo(ig,l))-zl(ig,l))
      ! AM Thetal
      zthl(ig, l) = pt(ig, l)/zpspsk(ig, l)

    END DO
  END DO

  ! print*,'1 OK convect8'
  ! --------------------


  ! + + + + + + + + + + +


  ! wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
  ! wh,wt,wo ...

  ! + + + + + + + + + + +  zh,zu,zv,zo,rho


  ! --------------------   zlev(1)
  ! \\\\\\\\\\\\\\\\\\\\



  ! -----------------------------------------------------------------------
  ! Calcul des altitudes des couches
  ! -----------------------------------------------------------------------

  DO l = 2, nlay
    DO ig = 1, ngrid
      zlev(ig, l) = 0.5*(pphi(ig,l)+pphi(ig,l-1))/rg
    END DO
  END DO
  DO ig = 1, ngrid
    zlev(ig, 1) = 0.
    zlev(ig, nlay+1) = (2.*pphi(ig,klev)-pphi(ig,klev-1))/rg
  END DO
  DO l = 1, nlay
    DO ig = 1, ngrid
      zlay(ig, l) = pphi(ig, l)/rg
    END DO
  END DO

  ! print*,'2 OK convect8'
  ! -----------------------------------------------------------------------
  ! Calcul des densites
  ! -----------------------------------------------------------------------

  DO l = 1, nlay
    DO ig = 1, ngrid
      ! rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*zh(ig,l))
      rho(ig, l) = pplay(ig, l)/(zpspsk(ig,l)*rd*ztv(ig,l))
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      rhobarz(ig, l) = 0.5*(rho(ig,l)+rho(ig,l-1))
    END DO
  END DO

  DO k = 1, nlay
    DO l = 1, nlay + 1
      DO ig = 1, ngrid
        wa(ig, k, l) = 0.
      END DO
    END DO
  END DO

  ! print*,'3 OK convect8'
  ! ------------------------------------------------------------------
  ! Calcul de w2, quarre de w a partir de la cape
  ! a partir de w2, on calcule wa, vitesse de l'ascendance

  ! ATTENTION: Dans cette version, pour cause d'economie de memoire,
  ! w2 est stoke dans wa

  ! ATTENTION: dans convect8, on n'utilise le calcule des wa
  ! independants par couches que pour calculer l'entrainement
  ! a la base et la hauteur max de l'ascendance.

  ! Indicages:
  ! l'ascendance provenant du niveau k traverse l'interface l avec
  ! une vitesse wa(k,l).

  ! --------------------

  ! + + + + + + + + + +

  ! wa(k,l)   ----       --------------------    l
  ! /\
  ! /||\       + + + + + + + + + +
  ! ||
  ! ||        --------------------
  ! ||
  ! ||        + + + + + + + + + +
  ! ||
  ! ||        --------------------
  ! ||__
  ! |___      + + + + + + + + + +     k

  ! --------------------



  ! ------------------------------------------------------------------

  ! CR: ponderation entrainement des couches instables
  ! def des entr_star tels que entr=f*entr_star
  DO l = 1, klev
    DO ig = 1, ngrid
      entr_star(ig, l) = 0.
    END DO
  END DO
  ! determination de la longueur de la couche d entrainement
  DO ig = 1, ngrid
    lentr(ig) = 1
  END DO

  ! on ne considere que les premieres couches instables
  DO k = nlay - 1, 1, -1
    DO ig = 1, ngrid
      IF (ztv(ig,k)>ztv(ig,k+1) .AND. ztv(ig,k+1)<ztv(ig,k+2)) THEN
        lentr(ig) = k
      END IF
    END DO
  END DO

  ! determination du lmin: couche d ou provient le thermique
  DO ig = 1, ngrid
    lmin(ig) = 1
  END DO
  DO ig = 1, ngrid
    DO l = nlay, 2, -1
      IF (ztv(ig,l-1)>ztv(ig,l)) THEN
        lmin(ig) = l - 1
      END IF
    END DO
  END DO

  ! definition de l'entrainement des couches
  DO l = 1, klev - 1
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. l>=lmin(ig) .AND. l<=lentr(ig)) THEN
        entr_star(ig, l) = (ztv(ig,l)-ztv(ig,l+1))*(zlev(ig,l+1)-zlev(ig,l))
      END IF
    END DO
  END DO
  ! pas de thermique si couche 1 stable
  DO ig = 1, ngrid
    IF (lmin(ig)>1) THEN
      DO l = 1, klev
        entr_star(ig, l) = 0.
      END DO
    END IF
  END DO
  ! calcul de l entrainement total
  DO ig = 1, ngrid
    entr_star_tot(ig) = 0.
  END DO
  DO ig = 1, ngrid
    DO k = 1, klev
      entr_star_tot(ig) = entr_star_tot(ig) + entr_star(ig, k)
    END DO
  END DO

  DO k = 1, klev
    DO ig = 1, ngrid
      ztva(ig, k) = ztv(ig, k)
    END DO
  END DO
  ! RC
  ! AM:initialisations
  DO k = 1, nlay
    DO ig = 1, ngrid
      ztva(ig, k) = ztv(ig, k)
      ztla(ig, k) = zthl(ig, k)
      zqla(ig, k) = 0.
      zqta(ig, k) = po(ig, k)
      zsat(ig) = .FALSE.
    END DO
  END DO

  ! print*,'7 OK convect8'
  DO k = 1, klev + 1
    DO ig = 1, ngrid
      zw2(ig, k) = 0.
      fmc(ig, k) = 0.
      ! CR
      f_star(ig, k) = 0.
      ! RC
      larg_cons(ig, k) = 0.
      larg_detr(ig, k) = 0.
      wa_moy(ig, k) = 0.
    END DO
  END DO

  ! print*,'8 OK convect8'
  DO ig = 1, ngrid
    linter(ig) = 1.
    lmaxa(ig) = 1
    lmix(ig) = 1
    wmaxa(ig) = 0.
  END DO

  ! CR:
  DO l = 1, nlay - 2
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. entr_star(ig,l)>1.E-10 .AND. &
          zw2(ig,l)<1E-10) THEN
        ! AM
        ztla(ig, l) = zthl(ig, l)
        zqta(ig, l) = po(ig, l)
        zqla(ig, l) = zl(ig, l)
        ! AM
        f_star(ig, l+1) = entr_star(ig, l)
        ! test:calcul de dteta
        zw2(ig, l+1) = 2.*rg*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig, l+1)* &
          (zlev(ig,l+1)-zlev(ig,l))*0.4*pphi(ig, l)/(pphi(ig,l+1)-pphi(ig,l))
        larg_detr(ig, l) = 0.
      ELSE IF ((zw2(ig,l)>=1E-10) .AND. (f_star(ig,l)+entr_star(ig, &
          l)>1.E-10)) THEN
        f_star(ig, l+1) = f_star(ig, l) + entr_star(ig, l)

        ! AM on melange Tl et qt du thermique
        ztla(ig, l) = (f_star(ig,l)*ztla(ig,l-1)+entr_star(ig,l)*zthl(ig,l))/ &
          f_star(ig, l+1)
        zqta(ig, l) = (f_star(ig,l)*zqta(ig,l-1)+entr_star(ig,l)*po(ig,l))/ &
          f_star(ig, l+1)

        ! ztva(ig,l)=(f_star(ig,l)*ztva(ig,l-1)+entr_star(ig,l)
        ! s                    *ztv(ig,l))/f_star(ig,l+1)

        ! AM on en deduit thetav et ql du thermique
        tbef(ig) = ztla(ig, l)*zpspsk(ig, l)
        zdelta = max(0., sign(1.,rtt-tbef(ig)))
        qsatbef(ig) = r2es*foeew(tbef(ig), zdelta)/pplev(ig, l)
        qsatbef(ig) = min(0.5, qsatbef(ig))
        zcor = 1./(1.-retv*qsatbef(ig))
        qsatbef(ig) = qsatbef(ig)*zcor
        zsat(ig) = (max(0.,zqta(ig,l)-qsatbef(ig))>0.00001)
      END IF
    END DO
    DO ig = 1, ngrid
      IF (zsat(ig)) THEN
        qlbef = max(0., zqta(ig,l)-qsatbef(ig))
        dt = 0.5*rlvcp*qlbef
        DO WHILE (dt>ddt0)
          tbef(ig) = tbef(ig) + dt
          zdelta = max(0., sign(1.,rtt-tbef(ig)))
          qsatbef(ig) = r2es*foeew(tbef(ig), zdelta)/pplev(ig, l)
          qsatbef(ig) = min(0.5, qsatbef(ig))
          zcor = 1./(1.-retv*qsatbef(ig))
          qsatbef(ig) = qsatbef(ig)*zcor
          qlbef = zqta(ig, l) - qsatbef(ig)

          zdelta = max(0., sign(1.,rtt-tbef(ig)))
          zcvm5 = r5les*(1.-zdelta) + r5ies*zdelta
          zcor = 1./(1.-retv*qsatbef(ig))
          dqsat_dt = foede(tbef(ig), zdelta, zcvm5, qsatbef(ig), zcor)
          num = -tbef(ig) + ztla(ig, l)*zpspsk(ig, l) + rlvcp*qlbef
          denom = 1. + rlvcp*dqsat_dt
          dt = num/denom
        END DO
        zqla(ig, l) = max(0., zqta(ig,l)-qsatbef(ig))
      END IF
      ! on ecrit de maniere conservative (sat ou non)
      ! T = Tl +Lv/Cp ql
      ztva(ig, l) = ztla(ig, l)*zpspsk(ig, l) + rlvcp*zqla(ig, l)
      ztva(ig, l) = ztva(ig, l)/zpspsk(ig, l)
      ztva(ig, l) = ztva(ig, l)*(1.+retv*(zqta(ig,l)-zqla(ig,l))-zqla(ig,l))

    END DO
    DO ig = 1, ngrid
      IF (zw2(ig,l)>=1.E-10 .AND. f_star(ig,l)+entr_star(ig,l)>1.E-10) THEN
        ! mise a jour de la vitesse ascendante (l'air entraine de la couche
        ! consideree commence avec une vitesse nulle).

        zw2(ig, l+1) = zw2(ig, l)*(f_star(ig,l)/f_star(ig,l+1))**2 + &
          2.*rg*(ztva(ig,l)-ztv(ig,l))/ztv(ig, l)*(zlev(ig,l+1)-zlev(ig,l))
      END IF
      ! determination de zmax continu par interpolation lineaire
      IF (zw2(ig,l+1)<0.) THEN
        linter(ig) = (l*(zw2(ig,l+1)-zw2(ig,l))-zw2(ig,l))/(zw2(ig,l+1)-zw2( &
          ig,l))
        zw2(ig, l+1) = 0.
        lmaxa(ig) = l
      ELSE
        wa_moy(ig, l+1) = sqrt(zw2(ig,l+1))
      END IF
      IF (wa_moy(ig,l+1)>wmaxa(ig)) THEN
        ! lmix est le niveau de la couche ou w (wa_moy) est maximum
        lmix(ig) = l + 1
        wmaxa(ig) = wa_moy(ig, l+1)
      END IF
    END DO
  END DO

  ! Calcul de la couche correspondant a la hauteur du thermique
  DO ig = 1, ngrid
    lmax(ig) = lentr(ig)
  END DO
  DO ig = 1, ngrid
    DO l = nlay, lentr(ig) + 1, -1
      IF (zw2(ig,l)<=1.E-10) THEN
        lmax(ig) = l - 1
      END IF
    END DO
  END DO
  ! pas de thermique si couche 1 stable
  DO ig = 1, ngrid
    IF (lmin(ig)>1) THEN
      lmax(ig) = 1
      lmin(ig) = 1
    END IF
  END DO

  ! Determination de zw2 max
  DO ig = 1, ngrid
    wmax(ig) = 0.
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (l<=lmax(ig)) THEN
        zw2(ig, l) = sqrt(zw2(ig,l))
        wmax(ig) = max(wmax(ig), zw2(ig,l))
      ELSE
        zw2(ig, l) = 0.
      END IF
    END DO
  END DO

  ! Longueur caracteristique correspondant a la hauteur des thermiques.
  DO ig = 1, ngrid
    zmax(ig) = 500.
    zlevinter(ig) = zlev(ig, 1)
  END DO
  DO ig = 1, ngrid
    ! calcul de zlevinter
    zlevinter(ig) = (zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*linter(ig) + &
      zlev(ig, lmax(ig)) - lmax(ig)*(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))
    zmax(ig) = max(zmax(ig), zlevinter(ig)-zlev(ig,lmin(ig)))
  END DO

  ! Fermeture,determination de f
  DO ig = 1, ngrid
    entr_star2(ig) = 0.
  END DO
  DO ig = 1, ngrid
    IF (entr_star_tot(ig)<1.E-10) THEN
      f(ig) = 0.
    ELSE
      DO k = lmin(ig), lentr(ig)
        entr_star2(ig) = entr_star2(ig) + entr_star(ig, k)**2/(rho(ig,k)*( &
          zlev(ig,k+1)-zlev(ig,k)))
      END DO
      ! Nouvelle fermeture
      f(ig) = wmax(ig)/(zmax(ig)*r_aspect*entr_star2(ig))*entr_star_tot(ig)
      ! test
      IF (first) THEN
        f(ig) = f(ig) + (f0(ig)-f(ig))*exp(-ptimestep/zmax(ig)*wmax(ig))
      END IF
    END IF
    f0(ig) = f(ig)
    first = .TRUE.
  END DO

  ! Calcul de l'entrainement
  DO k = 1, klev
    DO ig = 1, ngrid
      entr(ig, k) = f(ig)*entr_star(ig, k)
    END DO
  END DO
  ! Calcul des flux
  DO ig = 1, ngrid
    DO l = 1, lmax(ig) - 1
      fmc(ig, l+1) = fmc(ig, l) + entr(ig, l)
    END DO
  END DO

  ! RC


  ! print*,'9 OK convect8'
  ! print*,'WA1 ',wa_moy

  ! determination de l'indice du debut de la mixed layer ou w decroit

  ! calcul de la largeur de chaque ascendance dans le cas conservatif.
  ! dans ce cas simple, on suppose que la largeur de l'ascendance provenant
  ! d'une couche est égale à la hauteur de la couche alimentante.
  ! La vitesse maximale dans l'ascendance est aussi prise comme estimation
  ! de la vitesse d'entrainement horizontal dans la couche alimentante.

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (l<=lmaxa(ig)) THEN
        zw = max(wa_moy(ig,l), 1.E-10)
        larg_cons(ig, l) = zmax(ig)*r_aspect*fmc(ig, l)/(rhobarz(ig,l)*zw)
      END IF
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (l<=lmaxa(ig)) THEN
        ! if (idetr.eq.0) then
        ! cette option est finalement en dur.
        larg_detr(ig, l) = sqrt(l_mix*zlev(ig,l))
        ! else if (idetr.eq.1) then
        ! larg_detr(ig,l)=larg_cons(ig,l)
        ! s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
        ! else if (idetr.eq.2) then
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        ! s            *sqrt(wa_moy(ig,l))
        ! else if (idetr.eq.4) then
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        ! s            *wa_moy(ig,l)
        ! endif
      END IF
    END DO
  END DO

  ! print*,'10 OK convect8'
  ! print*,'WA2 ',wa_moy
  ! calcul de la fraction de la maille concernée par l'ascendance en tenant
  ! compte de l'epluchage du thermique.

  ! CR def de  zmix continu (profil parabolique des vitesses)
  DO ig = 1, ngrid
    IF (lmix(ig)>1.) THEN
      zmix(ig) = ((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig))) &
        **2-(zlev(ig,lmix(ig)+1))**2)-(zw2(ig,lmix(ig))-zw2(ig, &
        lmix(ig)+1))*((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))/ &
        (2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig)))- &
        (zlev(ig,lmix(ig)+1)))-(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))*((zlev( &
        ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
    ELSE
      zmix(ig) = 0.
    END IF
  END DO

  ! calcul du nouveau lmix correspondant
  DO ig = 1, ngrid
    DO l = 1, klev
      IF (zmix(ig)>=zlev(ig,l) .AND. zmix(ig)<zlev(ig,l+1)) THEN
        lmix(ig) = l
      END IF
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (larg_cons(ig,l)>1.) THEN
        ! print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
        fraca(ig, l) = (larg_cons(ig,l)-larg_detr(ig,l))/(r_aspect*zmax(ig))
        ! test
        fraca(ig, l) = max(fraca(ig,l), 0.)
        fraca(ig, l) = min(fraca(ig,l), 0.5)
        fracd(ig, l) = 1. - fraca(ig, l)
        fracc(ig, l) = larg_cons(ig, l)/(r_aspect*zmax(ig))
      ELSE
        ! wa_moy(ig,l)=0.
        fraca(ig, l) = 0.
        fracc(ig, l) = 0.
        fracd(ig, l) = 1.
      END IF
    END DO
  END DO
  ! CR: calcul de fracazmix
  DO ig = 1, ngrid
    fracazmix(ig) = (fraca(ig,lmix(ig)+1)-fraca(ig,lmix(ig)))/ &
      (zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))*zmix(ig) + &
      fraca(ig, lmix(ig)) - zlev(ig, lmix(ig))*(fraca(ig,lmix(ig)+1)-fraca(ig &
      ,lmix(ig)))/(zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (larg_cons(ig,l)>1.) THEN
        IF (l>lmix(ig)) THEN
          xxx(ig, l) = (zmax(ig)-zlev(ig,l))/(zmax(ig)-zmix(ig))
          IF (idetr==0) THEN
            fraca(ig, l) = fracazmix(ig)
          ELSE IF (idetr==1) THEN
            fraca(ig, l) = fracazmix(ig)*xxx(ig, l)
          ELSE IF (idetr==2) THEN
            fraca(ig, l) = fracazmix(ig)*(1.-(1.-xxx(ig,l))**2)
          ELSE
            fraca(ig, l) = fracazmix(ig)*xxx(ig, l)**2
          END IF
          ! print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
          fraca(ig, l) = max(fraca(ig,l), 0.)
          fraca(ig, l) = min(fraca(ig,l), 0.5)
          fracd(ig, l) = 1. - fraca(ig, l)
          fracc(ig, l) = larg_cons(ig, l)/(r_aspect*zmax(ig))
        END IF
      END IF
    END DO
  END DO

  ! print*,'11 OK convect8'
  ! print*,'Ea3 ',wa_moy
  ! ------------------------------------------------------------------
  ! Calcul de fracd, wd
  ! somme wa - wd = 0
  ! ------------------------------------------------------------------


  DO ig = 1, ngrid
    fm(ig, 1) = 0.
    fm(ig, nlay+1) = 0.
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      fm(ig, l) = fraca(ig, l)*wa_moy(ig, l)*rhobarz(ig, l)
      ! CR:test
      IF (entr(ig,l-1)<1E-10 .AND. fm(ig,l)>fm(ig,l-1) .AND. l>lmix(ig)) THEN
        fm(ig, l) = fm(ig, l-1)
        ! write(1,*)'ajustement fm, l',l
      END IF
      ! write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
      ! RC
    END DO
    DO ig = 1, ngrid
      IF (fracd(ig,l)<0.1) THEN
        abort_message = 'fracd trop petit'
        CALL abort_physic(modname, abort_message, 1)
      ELSE
        ! vitesse descendante "diagnostique"
        wd(ig, l) = fm(ig, l)/(fracd(ig,l)*rhobarz(ig,l))
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      ! masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
      masse(ig, l) = (pplev(ig,l)-pplev(ig,l+1))/rg
    END DO
  END DO

  ! print*,'12 OK convect8'
  ! print*,'WA4 ',wa_moy
  ! c------------------------------------------------------------------
  ! calcul du transport vertical
  ! ------------------------------------------------------------------

  GO TO 4444
  ! print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
  DO l = 2, nlay - 1
    DO ig = 1, ngrid
      IF (fm(ig,l+1)*ptimestep>masse(ig,l) .AND. fm(ig,l+1)*ptimestep>masse( &
          ig,l+1)) THEN
        ! print*,'WARN!!! FM>M ig=',ig,' l=',l,'  FM='
        ! s         ,fm(ig,l+1)*ptimestep
        ! s         ,'   M=',masse(ig,l),masse(ig,l+1)
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (entr(ig,l)*ptimestep>masse(ig,l)) THEN
        ! print*,'WARN!!! E>M ig=',ig,' l=',l,'  E=='
        ! s         ,entr(ig,l)*ptimestep
        ! s         ,'   M=',masse(ig,l)
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (.NOT. fm(ig,l)>=0. .OR. .NOT. fm(ig,l)<=10.) THEN
        ! print*,'WARN!!! fm exagere ig=',ig,'   l=',l
        ! s         ,'   FM=',fm(ig,l)
      END IF
      IF (.NOT. masse(ig,l)>=1.E-10 .OR. .NOT. masse(ig,l)<=1.E4) THEN
        ! print*,'WARN!!! masse exagere ig=',ig,'   l=',l
        ! s         ,'   M=',masse(ig,l)
        ! print*,'rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)',
        ! s                 rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)
        ! print*,'zlev(ig,l+1),zlev(ig,l)'
        ! s                ,zlev(ig,l+1),zlev(ig,l)
        ! print*,'pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)'
        ! s                ,pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)
      END IF
      IF (.NOT. entr(ig,l)>=0. .OR. .NOT. entr(ig,l)<=10.) THEN
        ! print*,'WARN!!! entr exagere ig=',ig,'   l=',l
        ! s         ,'   E=',entr(ig,l)
      END IF
    END DO
  END DO

4444 CONTINUE

  IF (w2di==1) THEN
    fm0 = fm0 + ptimestep*(fm-fm0)/tho
    entr0 = entr0 + ptimestep*(entr-entr0)/tho
  ELSE
    fm0 = fm
    entr0 = entr
  END IF

  IF (1==1) THEN
    ! call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
    ! .    ,zh,zdhadj,zha)
    ! call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
    ! .    ,zo,pdoadj,zoa)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zthl, &
      zdthladj, zta)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, po, pdoadj, &
      zoa)
  ELSE
    CALL dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zh, &
      zdhadj, zha)
    CALL dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zo, &
      pdoadj, zoa)
  END IF

  IF (1==0) THEN
    CALL dvthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zmax, &
      zu, zv, pduadj, pdvadj, zua, zva)
  ELSE
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zu, pduadj, &
      zua)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zv, pdvadj, &
      zva)
  END IF

  DO l = 1, nlay
    DO ig = 1, ngrid
      zf = 0.5*(fracc(ig,l)+fracc(ig,l+1))
      zf2 = zf/(1.-zf)
      thetath2(ig, l) = zf2*(zha(ig,l)-zh(ig,l))**2
      wth2(ig, l) = zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
    END DO
  END DO



  ! print*,'13 OK convect8'
  ! print*,'WA5 ',wa_moy
  DO l = 1, nlay
    DO ig = 1, ngrid
      ! pdtadj(ig,l)=zdhadj(ig,l)*zpspsk(ig,l)
      pdtadj(ig, l) = zdthladj(ig, l)*zpspsk(ig, l)
    END DO
  END DO


  ! do l=1,nlay
  ! do ig=1,ngrid
  ! if(abs(pdtadj(ig,l))*86400..gt.500.) then
  ! print*,'WARN!!! ig=',ig,'  l=',l
  ! s         ,'   pdtadj=',pdtadj(ig,l)
  ! endif
  ! if(abs(pdoadj(ig,l))*86400..gt.1.) then
  ! print*,'WARN!!! ig=',ig,'  l=',l
  ! s         ,'   pdoadj=',pdoadj(ig,l)
  ! endif
  ! enddo
  ! enddo

  ! print*,'14 OK convect8'
  ! ------------------------------------------------------------------
  ! Calculs pour les sorties
  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE thermcell_eau

SUBROUTINE thermcell(ngrid, nlay, ptimestep, pplay, pplev, pphi, pu, pv, pt, &
    po, pduadj, pdvadj, pdtadj, pdoadj, fm0, entr0 & ! s
                                                     ! ,pu_therm,pv_therm
    , r_aspect, l_mix, w2di, tho)

  USE dimphy
  IMPLICIT NONE

  ! =======================================================================

  ! Calcul du transport verticale dans la couche limite en presence
  ! de "thermiques" explicitement representes

  ! Réécriture à partir d'un listing papier à Habas, le 14/02/00

  ! le thermique est supposé homogène et dissipé par mélange avec
  ! son environnement. la longueur l_mix contrôle l'efficacité du
  ! mélange

  ! Le calcul du transport des différentes espèces se fait en prenant
  ! en compte:
  ! 1. un flux de masse montant
  ! 2. un flux de masse descendant
  ! 3. un entrainement
  ! 4. un detrainement

  ! =======================================================================

  ! -----------------------------------------------------------------------
  ! declarations:
  ! -------------

  include "YOMCST.h"

  ! arguments:
  ! ----------

  INTEGER ngrid, nlay, w2di
  REAL tho
  REAL ptimestep, l_mix, r_aspect
  REAL pt(ngrid, nlay), pdtadj(ngrid, nlay)
  REAL pu(ngrid, nlay), pduadj(ngrid, nlay)
  REAL pv(ngrid, nlay), pdvadj(ngrid, nlay)
  REAL po(ngrid, nlay), pdoadj(ngrid, nlay)
  REAL pplay(ngrid, nlay), pplev(ngrid, nlay+1)
  REAL pphi(ngrid, nlay)

  INTEGER idetr
  SAVE idetr
  DATA idetr/3/
  !$OMP THREADPRIVATE(idetr)

  ! local:
  ! ------

  INTEGER ig, k, l, lmaxa(klon), lmix(klon)
  REAL zsortie1d(klon)
  ! CR: on remplace lmax(klon,klev+1)
  INTEGER lmax(klon), lmin(klon), lentr(klon)
  REAL linter(klon)
  REAL zmix(klon), fracazmix(klon)
  ! RC
  REAL zmax(klon), zw, zz, zw2(klon, klev+1), ztva(klon, klev), zzz

  REAL zlev(klon, klev+1), zlay(klon, klev)
  REAL zh(klon, klev), zdhadj(klon, klev)
  REAL ztv(klon, klev)
  REAL zu(klon, klev), zv(klon, klev), zo(klon, klev)
  REAL wh(klon, klev+1)
  REAL wu(klon, klev+1), wv(klon, klev+1), wo(klon, klev+1)
  REAL zla(klon, klev+1)
  REAL zwa(klon, klev+1)
  REAL zld(klon, klev+1)
  REAL zwd(klon, klev+1)
  REAL zsortie(klon, klev)
  REAL zva(klon, klev)
  REAL zua(klon, klev)
  REAL zoa(klon, klev)

  REAL zha(klon, klev)
  REAL wa_moy(klon, klev+1)
  REAL fraca(klon, klev+1)
  REAL fracc(klon, klev+1)
  REAL zf, zf2
  REAL thetath2(klon, klev), wth2(klon, klev)
  ! common/comtherm/thetath2,wth2

  REAL count_time
  INTEGER ialt

  LOGICAL sorties
  REAL rho(klon, klev), rhobarz(klon, klev+1), masse(klon, klev)
  REAL zpspsk(klon, klev)

  ! real wmax(klon,klev),wmaxa(klon)
  REAL wmax(klon), wmaxa(klon)
  REAL wa(klon, klev, klev+1)
  REAL wd(klon, klev+1)
  REAL larg_part(klon, klev, klev+1)
  REAL fracd(klon, klev+1)
  REAL xxx(klon, klev+1)
  REAL larg_cons(klon, klev+1)
  REAL larg_detr(klon, klev+1)
  REAL fm0(klon, klev+1), entr0(klon, klev), detr(klon, klev)
  REAL pu_therm(klon, klev), pv_therm(klon, klev)
  REAL fm(klon, klev+1), entr(klon, klev)
  REAL fmc(klon, klev+1)

  ! CR:nouvelles variables
  REAL f_star(klon, klev+1), entr_star(klon, klev)
  REAL entr_star_tot(klon), entr_star2(klon)
  REAL f(klon), f0(klon)
  REAL zlevinter(klon)
  LOGICAL first
  DATA first/.FALSE./
  SAVE first
  !$OMP THREADPRIVATE(first)
  ! RC

  CHARACTER *2 str2
  CHARACTER *10 str10

  CHARACTER (LEN=20) :: modname = 'thermcell'
  CHARACTER (LEN=80) :: abort_message

  LOGICAL vtest(klon), down

  EXTERNAL scopy

  INTEGER ncorrec, ll
  SAVE ncorrec
  DATA ncorrec/0/
  !$OMP THREADPRIVATE(ncorrec)


  ! -----------------------------------------------------------------------
  ! initialisation:
  ! ---------------

  sorties = .TRUE.
  IF (ngrid/=klon) THEN
    PRINT *
    PRINT *, 'STOP dans convadj'
    PRINT *, 'ngrid    =', ngrid
    PRINT *, 'klon  =', klon
  END IF

  ! -----------------------------------------------------------------------
  ! incrementation eventuelle de tendances precedentes:
  ! ---------------------------------------------------

  ! print*,'0 OK convect8'

  DO l = 1, nlay
    DO ig = 1, ngrid
      zpspsk(ig, l) = (pplay(ig,l)/pplev(ig,1))**rkappa
      zh(ig, l) = pt(ig, l)/zpspsk(ig, l)
      zu(ig, l) = pu(ig, l)
      zv(ig, l) = pv(ig, l)
      zo(ig, l) = po(ig, l)
      ztv(ig, l) = zh(ig, l)*(1.+0.61*zo(ig,l))
    END DO
  END DO

  ! print*,'1 OK convect8'
  ! --------------------


  ! + + + + + + + + + + +


  ! wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
  ! wh,wt,wo ...

  ! + + + + + + + + + + +  zh,zu,zv,zo,rho


  ! --------------------   zlev(1)
  ! \\\\\\\\\\\\\\\\\\\\



  ! -----------------------------------------------------------------------
  ! Calcul des altitudes des couches
  ! -----------------------------------------------------------------------

  DO l = 2, nlay
    DO ig = 1, ngrid
      zlev(ig, l) = 0.5*(pphi(ig,l)+pphi(ig,l-1))/rg
    END DO
  END DO
  DO ig = 1, ngrid
    zlev(ig, 1) = 0.
    zlev(ig, nlay+1) = (2.*pphi(ig,klev)-pphi(ig,klev-1))/rg
  END DO
  DO l = 1, nlay
    DO ig = 1, ngrid
      zlay(ig, l) = pphi(ig, l)/rg
    END DO
  END DO

  ! print*,'2 OK convect8'
  ! -----------------------------------------------------------------------
  ! Calcul des densites
  ! -----------------------------------------------------------------------

  DO l = 1, nlay
    DO ig = 1, ngrid
      rho(ig, l) = pplay(ig, l)/(zpspsk(ig,l)*rd*zh(ig,l))
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      rhobarz(ig, l) = 0.5*(rho(ig,l)+rho(ig,l-1))
    END DO
  END DO

  DO k = 1, nlay
    DO l = 1, nlay + 1
      DO ig = 1, ngrid
        wa(ig, k, l) = 0.
      END DO
    END DO
  END DO

  ! print*,'3 OK convect8'
  ! ------------------------------------------------------------------
  ! Calcul de w2, quarre de w a partir de la cape
  ! a partir de w2, on calcule wa, vitesse de l'ascendance

  ! ATTENTION: Dans cette version, pour cause d'economie de memoire,
  ! w2 est stoke dans wa

  ! ATTENTION: dans convect8, on n'utilise le calcule des wa
  ! independants par couches que pour calculer l'entrainement
  ! a la base et la hauteur max de l'ascendance.

  ! Indicages:
  ! l'ascendance provenant du niveau k traverse l'interface l avec
  ! une vitesse wa(k,l).

  ! --------------------

  ! + + + + + + + + + +

  ! wa(k,l)   ----       --------------------    l
  ! /\
  ! /||\       + + + + + + + + + +
  ! ||
  ! ||        --------------------
  ! ||
  ! ||        + + + + + + + + + +
  ! ||
  ! ||        --------------------
  ! ||__
  ! |___      + + + + + + + + + +     k

  ! --------------------



  ! ------------------------------------------------------------------

  ! CR: ponderation entrainement des couches instables
  ! def des entr_star tels que entr=f*entr_star
  DO l = 1, klev
    DO ig = 1, ngrid
      entr_star(ig, l) = 0.
    END DO
  END DO
  ! determination de la longueur de la couche d entrainement
  DO ig = 1, ngrid
    lentr(ig) = 1
  END DO

  ! on ne considere que les premieres couches instables
  DO k = nlay - 2, 1, -1
    DO ig = 1, ngrid
      IF (ztv(ig,k)>ztv(ig,k+1) .AND. ztv(ig,k+1)<=ztv(ig,k+2)) THEN
        lentr(ig) = k
      END IF
    END DO
  END DO

  ! determination du lmin: couche d ou provient le thermique
  DO ig = 1, ngrid
    lmin(ig) = 1
  END DO
  DO ig = 1, ngrid
    DO l = nlay, 2, -1
      IF (ztv(ig,l-1)>ztv(ig,l)) THEN
        lmin(ig) = l - 1
      END IF
    END DO
  END DO

  ! definition de l'entrainement des couches
  DO l = 1, klev - 1
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. l>=lmin(ig) .AND. l<=lentr(ig)) THEN
        entr_star(ig, l) = (ztv(ig,l)-ztv(ig,l+1))*(zlev(ig,l+1)-zlev(ig,l))
      END IF
    END DO
  END DO
  ! pas de thermique si couches 1->5 stables
  DO ig = 1, ngrid
    IF (lmin(ig)>5) THEN
      DO l = 1, klev
        entr_star(ig, l) = 0.
      END DO
    END IF
  END DO
  ! calcul de l entrainement total
  DO ig = 1, ngrid
    entr_star_tot(ig) = 0.
  END DO
  DO ig = 1, ngrid
    DO k = 1, klev
      entr_star_tot(ig) = entr_star_tot(ig) + entr_star(ig, k)
    END DO
  END DO

  PRINT *, 'fin calcul entr_star'
  DO k = 1, klev
    DO ig = 1, ngrid
      ztva(ig, k) = ztv(ig, k)
    END DO
  END DO
  ! RC
  ! print*,'7 OK convect8'
  DO k = 1, klev + 1
    DO ig = 1, ngrid
      zw2(ig, k) = 0.
      fmc(ig, k) = 0.
      ! CR
      f_star(ig, k) = 0.
      ! RC
      larg_cons(ig, k) = 0.
      larg_detr(ig, k) = 0.
      wa_moy(ig, k) = 0.
    END DO
  END DO

  ! print*,'8 OK convect8'
  DO ig = 1, ngrid
    linter(ig) = 1.
    lmaxa(ig) = 1
    lmix(ig) = 1
    wmaxa(ig) = 0.
  END DO

  ! CR:
  DO l = 1, nlay - 2
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. entr_star(ig,l)>1.E-10 .AND. &
          zw2(ig,l)<1E-10) THEN
        f_star(ig, l+1) = entr_star(ig, l)
        ! test:calcul de dteta
        zw2(ig, l+1) = 2.*rg*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig, l+1)* &
          (zlev(ig,l+1)-zlev(ig,l))*0.4*pphi(ig, l)/(pphi(ig,l+1)-pphi(ig,l))
        larg_detr(ig, l) = 0.
      ELSE IF ((zw2(ig,l)>=1E-10) .AND. (f_star(ig,l)+entr_star(ig, &
          l)>1.E-10)) THEN
        f_star(ig, l+1) = f_star(ig, l) + entr_star(ig, l)
        ztva(ig, l) = (f_star(ig,l)*ztva(ig,l-1)+entr_star(ig,l)*ztv(ig,l))/ &
          f_star(ig, l+1)
        zw2(ig, l+1) = zw2(ig, l)*(f_star(ig,l)/f_star(ig,l+1))**2 + &
          2.*rg*(ztva(ig,l)-ztv(ig,l))/ztv(ig, l)*(zlev(ig,l+1)-zlev(ig,l))
      END IF
      ! determination de zmax continu par interpolation lineaire
      IF (zw2(ig,l+1)<0.) THEN
        ! test
        IF (abs(zw2(ig,l+1)-zw2(ig,l))<1E-10) THEN
          PRINT *, 'pb linter'
        END IF
        linter(ig) = (l*(zw2(ig,l+1)-zw2(ig,l))-zw2(ig,l))/(zw2(ig,l+1)-zw2( &
          ig,l))
        zw2(ig, l+1) = 0.
        lmaxa(ig) = l
      ELSE
        IF (zw2(ig,l+1)<0.) THEN
          PRINT *, 'pb1 zw2<0'
        END IF
        wa_moy(ig, l+1) = sqrt(zw2(ig,l+1))
      END IF
      IF (wa_moy(ig,l+1)>wmaxa(ig)) THEN
        ! lmix est le niveau de la couche ou w (wa_moy) est maximum
        lmix(ig) = l + 1
        wmaxa(ig) = wa_moy(ig, l+1)
      END IF
    END DO
  END DO
  PRINT *, 'fin calcul zw2'

  ! Calcul de la couche correspondant a la hauteur du thermique
  DO ig = 1, ngrid
    lmax(ig) = lentr(ig)
  END DO
  DO ig = 1, ngrid
    DO l = nlay, lentr(ig) + 1, -1
      IF (zw2(ig,l)<=1.E-10) THEN
        lmax(ig) = l - 1
      END IF
    END DO
  END DO
  ! pas de thermique si couches 1->5 stables
  DO ig = 1, ngrid
    IF (lmin(ig)>5) THEN
      lmax(ig) = 1
      lmin(ig) = 1
    END IF
  END DO

  ! Determination de zw2 max
  DO ig = 1, ngrid
    wmax(ig) = 0.
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (l<=lmax(ig)) THEN
        IF (zw2(ig,l)<0.) THEN
          PRINT *, 'pb2 zw2<0'
        END IF
        zw2(ig, l) = sqrt(zw2(ig,l))
        wmax(ig) = max(wmax(ig), zw2(ig,l))
      ELSE
        zw2(ig, l) = 0.
      END IF
    END DO
  END DO

  ! Longueur caracteristique correspondant a la hauteur des thermiques.
  DO ig = 1, ngrid
    zmax(ig) = 0.
    zlevinter(ig) = zlev(ig, 1)
  END DO
  DO ig = 1, ngrid
    ! calcul de zlevinter
    zlevinter(ig) = (zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*linter(ig) + &
      zlev(ig, lmax(ig)) - lmax(ig)*(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))
    zmax(ig) = max(zmax(ig), zlevinter(ig)-zlev(ig,lmin(ig)))
  END DO

  PRINT *, 'avant fermeture'
  ! Fermeture,determination de f
  DO ig = 1, ngrid
    entr_star2(ig) = 0.
  END DO
  DO ig = 1, ngrid
    IF (entr_star_tot(ig)<1.E-10) THEN
      f(ig) = 0.
    ELSE
      DO k = lmin(ig), lentr(ig)
        entr_star2(ig) = entr_star2(ig) + entr_star(ig, k)**2/(rho(ig,k)*( &
          zlev(ig,k+1)-zlev(ig,k)))
      END DO
      ! Nouvelle fermeture
      f(ig) = wmax(ig)/(max(500.,zmax(ig))*r_aspect*entr_star2(ig))* &
        entr_star_tot(ig)
      ! test
      ! if (first) then
      ! f(ig)=f(ig)+(f0(ig)-f(ig))*exp(-ptimestep/zmax(ig)
      ! s             *wmax(ig))
      ! endif
    END IF
    ! f0(ig)=f(ig)
    ! first=.true.
  END DO
  PRINT *, 'apres fermeture'

  ! Calcul de l'entrainement
  DO k = 1, klev
    DO ig = 1, ngrid
      entr(ig, k) = f(ig)*entr_star(ig, k)
    END DO
  END DO
  ! Calcul des flux
  DO ig = 1, ngrid
    DO l = 1, lmax(ig) - 1
      fmc(ig, l+1) = fmc(ig, l) + entr(ig, l)
    END DO
  END DO

  ! RC


  ! print*,'9 OK convect8'
  ! print*,'WA1 ',wa_moy

  ! determination de l'indice du debut de la mixed layer ou w decroit

  ! calcul de la largeur de chaque ascendance dans le cas conservatif.
  ! dans ce cas simple, on suppose que la largeur de l'ascendance provenant
  ! d'une couche est égale à la hauteur de la couche alimentante.
  ! La vitesse maximale dans l'ascendance est aussi prise comme estimation
  ! de la vitesse d'entrainement horizontal dans la couche alimentante.

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (l<=lmaxa(ig)) THEN
        zw = max(wa_moy(ig,l), 1.E-10)
        larg_cons(ig, l) = zmax(ig)*r_aspect*fmc(ig, l)/(rhobarz(ig,l)*zw)
      END IF
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (l<=lmaxa(ig)) THEN
        ! if (idetr.eq.0) then
        ! cette option est finalement en dur.
        IF ((l_mix*zlev(ig,l))<0.) THEN
          PRINT *, 'pb l_mix*zlev<0'
        END IF
        larg_detr(ig, l) = sqrt(l_mix*zlev(ig,l))
        ! else if (idetr.eq.1) then
        ! larg_detr(ig,l)=larg_cons(ig,l)
        ! s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
        ! else if (idetr.eq.2) then
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        ! s            *sqrt(wa_moy(ig,l))
        ! else if (idetr.eq.4) then
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        ! s            *wa_moy(ig,l)
        ! endif
      END IF
    END DO
  END DO

  ! print*,'10 OK convect8'
  ! print*,'WA2 ',wa_moy
  ! calcul de la fraction de la maille concernée par l'ascendance en tenant
  ! compte de l'epluchage du thermique.

  ! CR def de  zmix continu (profil parabolique des vitesses)
  DO ig = 1, ngrid
    IF (lmix(ig)>1.) THEN
      ! test
      IF (((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig)))- &
          (zlev(ig,lmix(ig)+1)))-(zw2(ig,lmix(ig))- &
          zw2(ig,lmix(ig)+1))*((zlev(ig,lmix(ig)-1))- &
          (zlev(ig,lmix(ig)))))>1E-10) THEN

        zmix(ig) = ((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig)) &
          )**2-(zlev(ig,lmix(ig)+1))**2)-(zw2(ig,lmix(ig))-zw2(ig, &
          lmix(ig)+1))*((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))/ &
          (2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig)))- &
          (zlev(ig,lmix(ig)+1)))-(zw2(ig,lmix(ig))- &
          zw2(ig,lmix(ig)+1))*((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
      ELSE
        zmix(ig) = zlev(ig, lmix(ig))
        PRINT *, 'pb zmix'
      END IF
    ELSE
      zmix(ig) = 0.
    END IF
    ! test
    IF ((zmax(ig)-zmix(ig))<0.) THEN
      zmix(ig) = 0.99*zmax(ig)
      ! print*,'pb zmix>zmax'
    END IF
  END DO

  ! calcul du nouveau lmix correspondant
  DO ig = 1, ngrid
    DO l = 1, klev
      IF (zmix(ig)>=zlev(ig,l) .AND. zmix(ig)<zlev(ig,l+1)) THEN
        lmix(ig) = l
      END IF
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (larg_cons(ig,l)>1.) THEN
        ! print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
        fraca(ig, l) = (larg_cons(ig,l)-larg_detr(ig,l))/(r_aspect*zmax(ig))
        ! test
        fraca(ig, l) = max(fraca(ig,l), 0.)
        fraca(ig, l) = min(fraca(ig,l), 0.5)
        fracd(ig, l) = 1. - fraca(ig, l)
        fracc(ig, l) = larg_cons(ig, l)/(r_aspect*zmax(ig))
      ELSE
        ! wa_moy(ig,l)=0.
        fraca(ig, l) = 0.
        fracc(ig, l) = 0.
        fracd(ig, l) = 1.
      END IF
    END DO
  END DO
  ! CR: calcul de fracazmix
  DO ig = 1, ngrid
    fracazmix(ig) = (fraca(ig,lmix(ig)+1)-fraca(ig,lmix(ig)))/ &
      (zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))*zmix(ig) + &
      fraca(ig, lmix(ig)) - zlev(ig, lmix(ig))*(fraca(ig,lmix(ig)+1)-fraca(ig &
      ,lmix(ig)))/(zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (larg_cons(ig,l)>1.) THEN
        IF (l>lmix(ig)) THEN
          ! test
          IF (zmax(ig)-zmix(ig)<1.E-10) THEN
            ! print*,'pb xxx'
            xxx(ig, l) = (lmaxa(ig)+1.-l)/(lmaxa(ig)+1.-lmix(ig))
          ELSE
            xxx(ig, l) = (zmax(ig)-zlev(ig,l))/(zmax(ig)-zmix(ig))
          END IF
          IF (idetr==0) THEN
            fraca(ig, l) = fracazmix(ig)
          ELSE IF (idetr==1) THEN
            fraca(ig, l) = fracazmix(ig)*xxx(ig, l)
          ELSE IF (idetr==2) THEN
            fraca(ig, l) = fracazmix(ig)*(1.-(1.-xxx(ig,l))**2)
          ELSE
            fraca(ig, l) = fracazmix(ig)*xxx(ig, l)**2
          END IF
          ! print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
          fraca(ig, l) = max(fraca(ig,l), 0.)
          fraca(ig, l) = min(fraca(ig,l), 0.5)
          fracd(ig, l) = 1. - fraca(ig, l)
          fracc(ig, l) = larg_cons(ig, l)/(r_aspect*zmax(ig))
        END IF
      END IF
    END DO
  END DO

  PRINT *, 'fin calcul fraca'
  ! print*,'11 OK convect8'
  ! print*,'Ea3 ',wa_moy
  ! ------------------------------------------------------------------
  ! Calcul de fracd, wd
  ! somme wa - wd = 0
  ! ------------------------------------------------------------------


  DO ig = 1, ngrid
    fm(ig, 1) = 0.
    fm(ig, nlay+1) = 0.
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      fm(ig, l) = fraca(ig, l)*wa_moy(ig, l)*rhobarz(ig, l)
      ! CR:test
      IF (entr(ig,l-1)<1E-10 .AND. fm(ig,l)>fm(ig,l-1) .AND. l>lmix(ig)) THEN
        fm(ig, l) = fm(ig, l-1)
        ! write(1,*)'ajustement fm, l',l
      END IF
      ! write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
      ! RC
    END DO
    DO ig = 1, ngrid
      IF (fracd(ig,l)<0.1) THEN
        abort_message = 'fracd trop petit'
        CALL abort_physic(modname, abort_message, 1)
      ELSE
        ! vitesse descendante "diagnostique"
        wd(ig, l) = fm(ig, l)/(fracd(ig,l)*rhobarz(ig,l))
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      ! masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
      masse(ig, l) = (pplev(ig,l)-pplev(ig,l+1))/rg
    END DO
  END DO

  ! print*,'12 OK convect8'
  ! print*,'WA4 ',wa_moy
  ! c------------------------------------------------------------------
  ! calcul du transport vertical
  ! ------------------------------------------------------------------

  GO TO 4444
  ! print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
  DO l = 2, nlay - 1
    DO ig = 1, ngrid
      IF (fm(ig,l+1)*ptimestep>masse(ig,l) .AND. fm(ig,l+1)*ptimestep>masse( &
          ig,l+1)) THEN
        ! print*,'WARN!!! FM>M ig=',ig,' l=',l,'  FM='
        ! s         ,fm(ig,l+1)*ptimestep
        ! s         ,'   M=',masse(ig,l),masse(ig,l+1)
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (entr(ig,l)*ptimestep>masse(ig,l)) THEN
        ! print*,'WARN!!! E>M ig=',ig,' l=',l,'  E=='
        ! s         ,entr(ig,l)*ptimestep
        ! s         ,'   M=',masse(ig,l)
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (.NOT. fm(ig,l)>=0. .OR. .NOT. fm(ig,l)<=10.) THEN
        ! print*,'WARN!!! fm exagere ig=',ig,'   l=',l
        ! s         ,'   FM=',fm(ig,l)
      END IF
      IF (.NOT. masse(ig,l)>=1.E-10 .OR. .NOT. masse(ig,l)<=1.E4) THEN
        ! print*,'WARN!!! masse exagere ig=',ig,'   l=',l
        ! s         ,'   M=',masse(ig,l)
        ! print*,'rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)',
        ! s                 rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)
        ! print*,'zlev(ig,l+1),zlev(ig,l)'
        ! s                ,zlev(ig,l+1),zlev(ig,l)
        ! print*,'pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)'
        ! s                ,pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)
      END IF
      IF (.NOT. entr(ig,l)>=0. .OR. .NOT. entr(ig,l)<=10.) THEN
        ! print*,'WARN!!! entr exagere ig=',ig,'   l=',l
        ! s         ,'   E=',entr(ig,l)
      END IF
    END DO
  END DO

4444 CONTINUE

  ! CR:redefinition du entr
  DO l = 1, nlay
    DO ig = 1, ngrid
      detr(ig, l) = fm(ig, l) + entr(ig, l) - fm(ig, l+1)
      IF (detr(ig,l)<0.) THEN
        entr(ig, l) = entr(ig, l) - detr(ig, l)
        detr(ig, l) = 0.
        ! print*,'WARNING !!! detrainement negatif ',ig,l
      END IF
    END DO
  END DO
  ! RC
  IF (w2di==1) THEN
    fm0 = fm0 + ptimestep*(fm-fm0)/tho
    entr0 = entr0 + ptimestep*(entr-entr0)/tho
  ELSE
    fm0 = fm
    entr0 = entr
  END IF

  IF (1==1) THEN
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zh, zdhadj, &
      zha)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zo, pdoadj, &
      zoa)
  ELSE
    CALL dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zh, &
      zdhadj, zha)
    CALL dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zo, &
      pdoadj, zoa)
  END IF

  IF (1==0) THEN
    CALL dvthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zmax, &
      zu, zv, pduadj, pdvadj, zua, zva)
  ELSE
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zu, pduadj, &
      zua)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zv, pdvadj, &
      zva)
  END IF

  DO l = 1, nlay
    DO ig = 1, ngrid
      zf = 0.5*(fracc(ig,l)+fracc(ig,l+1))
      zf2 = zf/(1.-zf)
      thetath2(ig, l) = zf2*(zha(ig,l)-zh(ig,l))**2
      wth2(ig, l) = zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
    END DO
  END DO



  ! print*,'13 OK convect8'
  ! print*,'WA5 ',wa_moy
  DO l = 1, nlay
    DO ig = 1, ngrid
      pdtadj(ig, l) = zdhadj(ig, l)*zpspsk(ig, l)
    END DO
  END DO


  ! do l=1,nlay
  ! do ig=1,ngrid
  ! if(abs(pdtadj(ig,l))*86400..gt.500.) then
  ! print*,'WARN!!! ig=',ig,'  l=',l
  ! s         ,'   pdtadj=',pdtadj(ig,l)
  ! endif
  ! if(abs(pdoadj(ig,l))*86400..gt.1.) then
  ! print*,'WARN!!! ig=',ig,'  l=',l
  ! s         ,'   pdoadj=',pdoadj(ig,l)
  ! endif
  ! enddo
  ! enddo

  ! print*,'14 OK convect8'
  ! ------------------------------------------------------------------
  ! Calculs pour les sorties
  ! ------------------------------------------------------------------

  IF (sorties) THEN
    DO l = 1, nlay
      DO ig = 1, ngrid
        zla(ig, l) = (1.-fracd(ig,l))*zmax(ig)
        zld(ig, l) = fracd(ig, l)*zmax(ig)
        IF (1.-fracd(ig,l)>1.E-10) zwa(ig, l) = wd(ig, l)*fracd(ig, l)/ &
          (1.-fracd(ig,l))
      END DO
    END DO

    ! deja fait
    ! do l=1,nlay
    ! do ig=1,ngrid
    ! detr(ig,l)=fm(ig,l)+entr(ig,l)-fm(ig,l+1)
    ! if (detr(ig,l).lt.0.) then
    ! entr(ig,l)=entr(ig,l)-detr(ig,l)
    ! detr(ig,l)=0.
    ! print*,'WARNING !!! detrainement negatif ',ig,l
    ! endif
    ! enddo
    ! enddo

    ! print*,'15 OK convect8'


    ! #define und
    GO TO 123
#ifdef und
    CALL writeg1d(1, nlay, wd, 'wd      ', 'wd      ')
    CALL writeg1d(1, nlay, zwa, 'wa      ', 'wa      ')
    CALL writeg1d(1, nlay, fracd, 'fracd      ', 'fracd      ')
    CALL writeg1d(1, nlay, fraca, 'fraca      ', 'fraca      ')
    CALL writeg1d(1, nlay, wa_moy, 'wam         ', 'wam         ')
    CALL writeg1d(1, nlay, zla, 'la      ', 'la      ')
    CALL writeg1d(1, nlay, zld, 'ld      ', 'ld      ')
    CALL writeg1d(1, nlay, pt, 'pt      ', 'pt      ')
    CALL writeg1d(1, nlay, zh, 'zh      ', 'zh      ')
    CALL writeg1d(1, nlay, zha, 'zha      ', 'zha      ')
    CALL writeg1d(1, nlay, zu, 'zu      ', 'zu      ')
    CALL writeg1d(1, nlay, zv, 'zv      ', 'zv      ')
    CALL writeg1d(1, nlay, zo, 'zo      ', 'zo      ')
    CALL writeg1d(1, nlay, wh, 'wh      ', 'wh      ')
    CALL writeg1d(1, nlay, wu, 'wu      ', 'wu      ')
    CALL writeg1d(1, nlay, wv, 'wv      ', 'wv      ')
    CALL writeg1d(1, nlay, wo, 'w15uo     ', 'wXo     ')
    CALL writeg1d(1, nlay, zdhadj, 'zdhadj      ', 'zdhadj      ')
    CALL writeg1d(1, nlay, pduadj, 'pduadj      ', 'pduadj      ')
    CALL writeg1d(1, nlay, pdvadj, 'pdvadj      ', 'pdvadj      ')
    CALL writeg1d(1, nlay, pdoadj, 'pdoadj      ', 'pdoadj      ')
    CALL writeg1d(1, nlay, entr, 'entr        ', 'entr        ')
    CALL writeg1d(1, nlay, detr, 'detr        ', 'detr        ')
    CALL writeg1d(1, nlay, fm, 'fm          ', 'fm          ')

    CALL writeg1d(1, nlay, pdtadj, 'pdtadj    ', 'pdtadj    ')
    CALL writeg1d(1, nlay, pplay, 'pplay     ', 'pplay     ')
    CALL writeg1d(1, nlay, pplev, 'pplev     ', 'pplev     ')

    ! recalcul des flux en diagnostique...
    ! print*,'PAS DE TEMPS ',ptimestep
    CALL dt2f(pplev, pplay, pt, pdtadj, wh)
    CALL writeg1d(1, nlay, wh, 'wh2     ', 'wh2     ')
#endif
123 CONTINUE

  END IF

  ! if(wa_moy(1,4).gt.1.e-10) stop

  ! print*,'19 OK convect8'
  RETURN
END SUBROUTINE thermcell

SUBROUTINE dqthermcell(ngrid, nlay, ptimestep, fm, entr, masse, q, dq, qa)
  USE dimphy
  IMPLICIT NONE

  ! =======================================================================

  ! Calcul du transport verticale dans la couche limite en presence
  ! de "thermiques" explicitement representes
  ! calcul du dq/dt une fois qu'on connait les ascendances

  ! =======================================================================

  INTEGER ngrid, nlay

  REAL ptimestep
  REAL masse(ngrid, nlay), fm(ngrid, nlay+1)
  REAL entr(ngrid, nlay)
  REAL q(ngrid, nlay)
  REAL dq(ngrid, nlay)

  REAL qa(klon, klev), detr(klon, klev), wqd(klon, klev+1)

  INTEGER ig, k

  ! calcul du detrainement

  DO k = 1, nlay
    DO ig = 1, ngrid
      detr(ig, k) = fm(ig, k) - fm(ig, k+1) + entr(ig, k)
      ! test
      IF (detr(ig,k)<0.) THEN
        entr(ig, k) = entr(ig, k) - detr(ig, k)
        detr(ig, k) = 0.
        ! print*,'detr2<0!!!','ig=',ig,'k=',k,'f=',fm(ig,k),
        ! s         'f+1=',fm(ig,k+1),'e=',entr(ig,k),'d=',detr(ig,k)
      END IF
      IF (fm(ig,k+1)<0.) THEN
        ! print*,'fm2<0!!!'
      END IF
      IF (entr(ig,k)<0.) THEN
        ! print*,'entr2<0!!!'
      END IF
    END DO
  END DO

  ! calcul de la valeur dans les ascendances
  DO ig = 1, ngrid
    qa(ig, 1) = q(ig, 1)
  END DO

  DO k = 2, nlay
    DO ig = 1, ngrid
      IF ((fm(ig,k+1)+detr(ig,k))*ptimestep>1.E-5*masse(ig,k)) THEN
        qa(ig, k) = (fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))/ &
          (fm(ig,k+1)+detr(ig,k))
      ELSE
        qa(ig, k) = q(ig, k)
      END IF
      IF (qa(ig,k)<0.) THEN
        ! print*,'qa<0!!!'
      END IF
      IF (q(ig,k)<0.) THEN
        ! print*,'q<0!!!'
      END IF
    END DO
  END DO

  DO k = 2, nlay
    DO ig = 1, ngrid
      ! wqd(ig,k)=fm(ig,k)*0.5*(q(ig,k-1)+q(ig,k))
      wqd(ig, k) = fm(ig, k)*q(ig, k)
      IF (wqd(ig,k)<0.) THEN
        ! print*,'wqd<0!!!'
      END IF
    END DO
  END DO
  DO ig = 1, ngrid
    wqd(ig, 1) = 0.
    wqd(ig, nlay+1) = 0.
  END DO

  DO k = 1, nlay
    DO ig = 1, ngrid
      dq(ig, k) = (detr(ig,k)*qa(ig,k)-entr(ig,k)*q(ig,k)-wqd(ig,k)+wqd(ig,k+ &
        1))/masse(ig, k)
      ! if (dq(ig,k).lt.0.) then
      ! print*,'dq<0!!!'
      ! endif
    END DO
  END DO

  RETURN
END SUBROUTINE dqthermcell
SUBROUTINE dvthermcell(ngrid, nlay, ptimestep, fm, entr, masse, fraca, larga, &
    u, v, du, dv, ua, va)
  USE dimphy
  IMPLICIT NONE

  ! =======================================================================

  ! Calcul du transport verticale dans la couche limite en presence
  ! de "thermiques" explicitement representes
  ! calcul du dq/dt une fois qu'on connait les ascendances

  ! =======================================================================

  INTEGER ngrid, nlay

  REAL ptimestep
  REAL masse(ngrid, nlay), fm(ngrid, nlay+1)
  REAL fraca(ngrid, nlay+1)
  REAL larga(ngrid)
  REAL entr(ngrid, nlay)
  REAL u(ngrid, nlay)
  REAL ua(ngrid, nlay)
  REAL du(ngrid, nlay)
  REAL v(ngrid, nlay)
  REAL va(ngrid, nlay)
  REAL dv(ngrid, nlay)

  REAL qa(klon, klev), detr(klon, klev)
  REAL wvd(klon, klev+1), wud(klon, klev+1)
  REAL gamma0, gamma(klon, klev+1)
  REAL dua, dva
  INTEGER iter

  INTEGER ig, k

  ! calcul du detrainement

  DO k = 1, nlay
    DO ig = 1, ngrid
      detr(ig, k) = fm(ig, k) - fm(ig, k+1) + entr(ig, k)
    END DO
  END DO

  ! calcul de la valeur dans les ascendances
  DO ig = 1, ngrid
    ua(ig, 1) = u(ig, 1)
    va(ig, 1) = v(ig, 1)
  END DO

  DO k = 2, nlay
    DO ig = 1, ngrid
      IF ((fm(ig,k+1)+detr(ig,k))*ptimestep>1.E-5*masse(ig,k)) THEN
        ! On itère sur la valeur du coeff de freinage.
        ! gamma0=rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k))
        gamma0 = masse(ig, k)*sqrt(0.5*(fraca(ig,k+1)+fraca(ig, &
          k)))*0.5/larga(ig)
        ! gamma0=0.
        ! la première fois on multiplie le coefficient de freinage
        ! par le module du vent dans la couche en dessous.
        dua = ua(ig, k-1) - u(ig, k-1)
        dva = va(ig, k-1) - v(ig, k-1)
        DO iter = 1, 5
          gamma(ig, k) = gamma0*sqrt(dua**2+dva**2)
          ua(ig, k) = (fm(ig,k)*ua(ig,k-1)+(entr(ig,k)+gamma(ig, &
            k))*u(ig,k))/(fm(ig,k+1)+detr(ig,k)+gamma(ig,k))
          va(ig, k) = (fm(ig,k)*va(ig,k-1)+(entr(ig,k)+gamma(ig, &
            k))*v(ig,k))/(fm(ig,k+1)+detr(ig,k)+gamma(ig,k))
          ! print*,k,ua(ig,k),va(ig,k),u(ig,k),v(ig,k),dua,dva
          dua = ua(ig, k) - u(ig, k)
          dva = va(ig, k) - v(ig, k)
        END DO
      ELSE
        ua(ig, k) = u(ig, k)
        va(ig, k) = v(ig, k)
        gamma(ig, k) = 0.
      END IF
    END DO
  END DO

  DO k = 2, nlay
    DO ig = 1, ngrid
      wud(ig, k) = fm(ig, k)*u(ig, k)
      wvd(ig, k) = fm(ig, k)*v(ig, k)
    END DO
  END DO
  DO ig = 1, ngrid
    wud(ig, 1) = 0.
    wud(ig, nlay+1) = 0.
    wvd(ig, 1) = 0.
    wvd(ig, nlay+1) = 0.
  END DO

  DO k = 1, nlay
    DO ig = 1, ngrid
      du(ig, k) = ((detr(ig,k)+gamma(ig,k))*ua(ig,k)-(entr(ig,k)+gamma(ig, &
        k))*u(ig,k)-wud(ig,k)+wud(ig,k+1))/masse(ig, k)
      dv(ig, k) = ((detr(ig,k)+gamma(ig,k))*va(ig,k)-(entr(ig,k)+gamma(ig, &
        k))*v(ig,k)-wvd(ig,k)+wvd(ig,k+1))/masse(ig, k)
    END DO
  END DO

  RETURN
END SUBROUTINE dvthermcell
SUBROUTINE dqthermcell2(ngrid, nlay, ptimestep, fm, entr, masse, frac, q, dq, &
    qa)
  USE dimphy
  IMPLICIT NONE

  ! =======================================================================

  ! Calcul du transport verticale dans la couche limite en presence
  ! de "thermiques" explicitement representes
  ! calcul du dq/dt une fois qu'on connait les ascendances

  ! =======================================================================

  INTEGER ngrid, nlay

  REAL ptimestep
  REAL masse(ngrid, nlay), fm(ngrid, nlay+1)
  REAL entr(ngrid, nlay), frac(ngrid, nlay)
  REAL q(ngrid, nlay)
  REAL dq(ngrid, nlay)

  REAL qa(klon, klev), detr(klon, klev), wqd(klon, klev+1)
  REAL qe(klon, klev), zf, zf2

  INTEGER ig, k

  ! calcul du detrainement

  DO k = 1, nlay
    DO ig = 1, ngrid
      detr(ig, k) = fm(ig, k) - fm(ig, k+1) + entr(ig, k)
    END DO
  END DO

  ! calcul de la valeur dans les ascendances
  DO ig = 1, ngrid
    qa(ig, 1) = q(ig, 1)
    qe(ig, 1) = q(ig, 1)
  END DO

  DO k = 2, nlay
    DO ig = 1, ngrid
      IF ((fm(ig,k+1)+detr(ig,k))*ptimestep>1.E-5*masse(ig,k)) THEN
        zf = 0.5*(frac(ig,k)+frac(ig,k+1))
        zf2 = 1./(1.-zf)
        qa(ig, k) = (fm(ig,k)*qa(ig,k-1)+zf2*entr(ig,k)*q(ig,k))/ &
          (fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2)
        qe(ig, k) = (q(ig,k)-zf*qa(ig,k))*zf2
      ELSE
        qa(ig, k) = q(ig, k)
        qe(ig, k) = q(ig, k)
      END IF
    END DO
  END DO

  DO k = 2, nlay
    DO ig = 1, ngrid
      ! wqd(ig,k)=fm(ig,k)*0.5*(q(ig,k-1)+q(ig,k))
      wqd(ig, k) = fm(ig, k)*qe(ig, k)
    END DO
  END DO
  DO ig = 1, ngrid
    wqd(ig, 1) = 0.
    wqd(ig, nlay+1) = 0.
  END DO

  DO k = 1, nlay
    DO ig = 1, ngrid
      dq(ig, k) = (detr(ig,k)*qa(ig,k)-entr(ig,k)*qe(ig,k)-wqd(ig,k)+wqd(ig,k &
        +1))/masse(ig, k)
    END DO
  END DO

  RETURN
END SUBROUTINE dqthermcell2
SUBROUTINE dvthermcell2(ngrid, nlay, ptimestep, fm, entr, masse, fraca, &
    larga, u, v, du, dv, ua, va)
  USE dimphy
  IMPLICIT NONE

  ! =======================================================================

  ! Calcul du transport verticale dans la couche limite en presence
  ! de "thermiques" explicitement representes
  ! calcul du dq/dt une fois qu'on connait les ascendances

  ! =======================================================================

  INTEGER ngrid, nlay

  REAL ptimestep
  REAL masse(ngrid, nlay), fm(ngrid, nlay+1)
  REAL fraca(ngrid, nlay+1)
  REAL larga(ngrid)
  REAL entr(ngrid, nlay)
  REAL u(ngrid, nlay)
  REAL ua(ngrid, nlay)
  REAL du(ngrid, nlay)
  REAL v(ngrid, nlay)
  REAL va(ngrid, nlay)
  REAL dv(ngrid, nlay)

  REAL qa(klon, klev), detr(klon, klev), zf, zf2
  REAL wvd(klon, klev+1), wud(klon, klev+1)
  REAL gamma0, gamma(klon, klev+1)
  REAL ue(klon, klev), ve(klon, klev)
  REAL dua, dva
  INTEGER iter

  INTEGER ig, k

  ! calcul du detrainement

  DO k = 1, nlay
    DO ig = 1, ngrid
      detr(ig, k) = fm(ig, k) - fm(ig, k+1) + entr(ig, k)
    END DO
  END DO

  ! calcul de la valeur dans les ascendances
  DO ig = 1, ngrid
    ua(ig, 1) = u(ig, 1)
    va(ig, 1) = v(ig, 1)
    ue(ig, 1) = u(ig, 1)
    ve(ig, 1) = v(ig, 1)
  END DO

  DO k = 2, nlay
    DO ig = 1, ngrid
      IF ((fm(ig,k+1)+detr(ig,k))*ptimestep>1.E-5*masse(ig,k)) THEN
        ! On itère sur la valeur du coeff de freinage.
        ! gamma0=rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k))
        gamma0 = masse(ig, k)*sqrt(0.5*(fraca(ig,k+1)+fraca(ig, &
          k)))*0.5/larga(ig)*1.
        ! s         *0.5
        ! gamma0=0.
        zf = 0.5*(fraca(ig,k)+fraca(ig,k+1))
        zf = 0.
        zf2 = 1./(1.-zf)
        ! la première fois on multiplie le coefficient de freinage
        ! par le module du vent dans la couche en dessous.
        dua = ua(ig, k-1) - u(ig, k-1)
        dva = va(ig, k-1) - v(ig, k-1)
        DO iter = 1, 5
          ! On choisit une relaxation lineaire.
          gamma(ig, k) = gamma0
          ! On choisit une relaxation quadratique.
          gamma(ig, k) = gamma0*sqrt(dua**2+dva**2)
          ua(ig, k) = (fm(ig,k)*ua(ig,k-1)+(zf2*entr(ig,k)+gamma(ig, &
            k))*u(ig,k))/(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2+gamma(ig,k) &
            )
          va(ig, k) = (fm(ig,k)*va(ig,k-1)+(zf2*entr(ig,k)+gamma(ig, &
            k))*v(ig,k))/(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2+gamma(ig,k) &
            )
          ! print*,k,ua(ig,k),va(ig,k),u(ig,k),v(ig,k),dua,dva
          dua = ua(ig, k) - u(ig, k)
          dva = va(ig, k) - v(ig, k)
          ue(ig, k) = (u(ig,k)-zf*ua(ig,k))*zf2
          ve(ig, k) = (v(ig,k)-zf*va(ig,k))*zf2
        END DO
      ELSE
        ua(ig, k) = u(ig, k)
        va(ig, k) = v(ig, k)
        ue(ig, k) = u(ig, k)
        ve(ig, k) = v(ig, k)
        gamma(ig, k) = 0.
      END IF
    END DO
  END DO

  DO k = 2, nlay
    DO ig = 1, ngrid
      wud(ig, k) = fm(ig, k)*ue(ig, k)
      wvd(ig, k) = fm(ig, k)*ve(ig, k)
    END DO
  END DO
  DO ig = 1, ngrid
    wud(ig, 1) = 0.
    wud(ig, nlay+1) = 0.
    wvd(ig, 1) = 0.
    wvd(ig, nlay+1) = 0.
  END DO

  DO k = 1, nlay
    DO ig = 1, ngrid
      du(ig, k) = ((detr(ig,k)+gamma(ig,k))*ua(ig,k)-(entr(ig,k)+gamma(ig, &
        k))*ue(ig,k)-wud(ig,k)+wud(ig,k+1))/masse(ig, k)
      dv(ig, k) = ((detr(ig,k)+gamma(ig,k))*va(ig,k)-(entr(ig,k)+gamma(ig, &
        k))*ve(ig,k)-wvd(ig,k)+wvd(ig,k+1))/masse(ig, k)
    END DO
  END DO

  RETURN
END SUBROUTINE dvthermcell2
SUBROUTINE thermcell_sec(ngrid, nlay, ptimestep, pplay, pplev, pphi, zlev, &
    pu, pv, pt, po, pduadj, pdvadj, pdtadj, pdoadj, fm0, entr0 & ! s
                                                                 ! ,pu_therm,pv_therm
    , r_aspect, l_mix, w2di, tho)

  USE dimphy
  IMPLICIT NONE

  ! =======================================================================

  ! Calcul du transport verticale dans la couche limite en presence
  ! de "thermiques" explicitement representes

  ! Réécriture à partir d'un listing papier à Habas, le 14/02/00

  ! le thermique est supposé homogène et dissipé par mélange avec
  ! son environnement. la longueur l_mix contrôle l'efficacité du
  ! mélange

  ! Le calcul du transport des différentes espèces se fait en prenant
  ! en compte:
  ! 1. un flux de masse montant
  ! 2. un flux de masse descendant
  ! 3. un entrainement
  ! 4. un detrainement

  ! =======================================================================

  ! -----------------------------------------------------------------------
  ! declarations:
  ! -------------

  include "YOMCST.h"

  ! arguments:
  ! ----------

  INTEGER ngrid, nlay, w2di
  REAL tho
  REAL ptimestep, l_mix, r_aspect
  REAL pt(ngrid, nlay), pdtadj(ngrid, nlay)
  REAL pu(ngrid, nlay), pduadj(ngrid, nlay)
  REAL pv(ngrid, nlay), pdvadj(ngrid, nlay)
  REAL po(ngrid, nlay), pdoadj(ngrid, nlay)
  REAL pplay(ngrid, nlay), pplev(ngrid, nlay+1)
  REAL pphi(ngrid, nlay)

  INTEGER idetr
  SAVE idetr
  DATA idetr/3/
  !$OMP THREADPRIVATE(idetr)

  ! local:
  ! ------

  INTEGER ig, k, l, lmaxa(klon), lmix(klon)
  REAL zsortie1d(klon)
  ! CR: on remplace lmax(klon,klev+1)
  INTEGER lmax(klon), lmin(klon), lentr(klon)
  REAL linter(klon)
  REAL zmix(klon), fracazmix(klon)
  ! RC
  REAL zmax(klon), zw, zz, zw2(klon, klev+1), ztva(klon, klev), zzz

  REAL zlev(klon, klev+1), zlay(klon, klev)
  REAL zh(klon, klev), zdhadj(klon, klev)
  REAL ztv(klon, klev)
  REAL zu(klon, klev), zv(klon, klev), zo(klon, klev)
  REAL wh(klon, klev+1)
  REAL wu(klon, klev+1), wv(klon, klev+1), wo(klon, klev+1)
  REAL zla(klon, klev+1)
  REAL zwa(klon, klev+1)
  REAL zld(klon, klev+1)
  REAL zwd(klon, klev+1)
  REAL zsortie(klon, klev)
  REAL zva(klon, klev)
  REAL zua(klon, klev)
  REAL zoa(klon, klev)

  REAL zha(klon, klev)
  REAL wa_moy(klon, klev+1)
  REAL fraca(klon, klev+1)
  REAL fracc(klon, klev+1)
  REAL zf, zf2
  REAL thetath2(klon, klev), wth2(klon, klev)
  ! common/comtherm/thetath2,wth2

  REAL count_time
  INTEGER ialt

  LOGICAL sorties
  REAL rho(klon, klev), rhobarz(klon, klev+1), masse(klon, klev)
  REAL zpspsk(klon, klev)

  ! real wmax(klon,klev),wmaxa(klon)
  REAL wmax(klon), wmaxa(klon)
  REAL wa(klon, klev, klev+1)
  REAL wd(klon, klev+1)
  REAL larg_part(klon, klev, klev+1)
  REAL fracd(klon, klev+1)
  REAL xxx(klon, klev+1)
  REAL larg_cons(klon, klev+1)
  REAL larg_detr(klon, klev+1)
  REAL fm0(klon, klev+1), entr0(klon, klev), detr(klon, klev)
  REAL pu_therm(klon, klev), pv_therm(klon, klev)
  REAL fm(klon, klev+1), entr(klon, klev)
  REAL fmc(klon, klev+1)

  ! CR:nouvelles variables
  REAL f_star(klon, klev+1), entr_star(klon, klev)
  REAL entr_star_tot(klon), entr_star2(klon)
  REAL f(klon), f0(klon)
  REAL zlevinter(klon)
  LOGICAL first
  DATA first/.FALSE./
  SAVE first
  !$OMP THREADPRIVATE(first)
  ! RC

  CHARACTER *2 str2
  CHARACTER *10 str10

  CHARACTER (LEN=20) :: modname = 'thermcell_sec'
  CHARACTER (LEN=80) :: abort_message

  LOGICAL vtest(klon), down

  EXTERNAL scopy

  INTEGER ncorrec, ll
  SAVE ncorrec
  DATA ncorrec/0/
  !$OMP THREADPRIVATE(ncorrec)


  ! -----------------------------------------------------------------------
  ! initialisation:
  ! ---------------

  sorties = .TRUE.
  IF (ngrid/=klon) THEN
    PRINT *
    PRINT *, 'STOP dans convadj'
    PRINT *, 'ngrid    =', ngrid
    PRINT *, 'klon  =', klon
  END IF

  ! -----------------------------------------------------------------------
  ! incrementation eventuelle de tendances precedentes:
  ! ---------------------------------------------------

  ! print*,'0 OK convect8'

  DO l = 1, nlay
    DO ig = 1, ngrid
      zpspsk(ig, l) = (pplay(ig,l)/pplev(ig,1))**rkappa
      zh(ig, l) = pt(ig, l)/zpspsk(ig, l)
      zu(ig, l) = pu(ig, l)
      zv(ig, l) = pv(ig, l)
      zo(ig, l) = po(ig, l)
      ztv(ig, l) = zh(ig, l)*(1.+0.61*zo(ig,l))
    END DO
  END DO

  ! print*,'1 OK convect8'
  ! --------------------


  ! + + + + + + + + + + +


  ! wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
  ! wh,wt,wo ...

  ! + + + + + + + + + + +  zh,zu,zv,zo,rho


  ! --------------------   zlev(1)
  ! \\\\\\\\\\\\\\\\\\\\



  ! -----------------------------------------------------------------------
  ! Calcul des altitudes des couches
  ! -----------------------------------------------------------------------

  DO l = 2, nlay
    DO ig = 1, ngrid
      zlev(ig, l) = 0.5*(pphi(ig,l)+pphi(ig,l-1))/rg
    END DO
  END DO
  DO ig = 1, ngrid
    zlev(ig, 1) = 0.
    zlev(ig, nlay+1) = (2.*pphi(ig,klev)-pphi(ig,klev-1))/rg
  END DO
  DO l = 1, nlay
    DO ig = 1, ngrid
      zlay(ig, l) = pphi(ig, l)/rg
    END DO
  END DO

  ! print*,'2 OK convect8'
  ! -----------------------------------------------------------------------
  ! Calcul des densites
  ! -----------------------------------------------------------------------

  DO l = 1, nlay
    DO ig = 1, ngrid
      rho(ig, l) = pplay(ig, l)/(zpspsk(ig,l)*rd*zh(ig,l))
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      rhobarz(ig, l) = 0.5*(rho(ig,l)+rho(ig,l-1))
    END DO
  END DO

  DO k = 1, nlay
    DO l = 1, nlay + 1
      DO ig = 1, ngrid
        wa(ig, k, l) = 0.
      END DO
    END DO
  END DO

  ! print*,'3 OK convect8'
  ! ------------------------------------------------------------------
  ! Calcul de w2, quarre de w a partir de la cape
  ! a partir de w2, on calcule wa, vitesse de l'ascendance

  ! ATTENTION: Dans cette version, pour cause d'economie de memoire,
  ! w2 est stoke dans wa

  ! ATTENTION: dans convect8, on n'utilise le calcule des wa
  ! independants par couches que pour calculer l'entrainement
  ! a la base et la hauteur max de l'ascendance.

  ! Indicages:
  ! l'ascendance provenant du niveau k traverse l'interface l avec
  ! une vitesse wa(k,l).

  ! --------------------

  ! + + + + + + + + + +

  ! wa(k,l)   ----       --------------------    l
  ! /\
  ! /||\       + + + + + + + + + +
  ! ||
  ! ||        --------------------
  ! ||
  ! ||        + + + + + + + + + +
  ! ||
  ! ||        --------------------
  ! ||__
  ! |___      + + + + + + + + + +     k

  ! --------------------



  ! ------------------------------------------------------------------

  ! CR: ponderation entrainement des couches instables
  ! def des entr_star tels que entr=f*entr_star
  DO l = 1, klev
    DO ig = 1, ngrid
      entr_star(ig, l) = 0.
    END DO
  END DO
  ! determination de la longueur de la couche d entrainement
  DO ig = 1, ngrid
    lentr(ig) = 1
  END DO

  ! on ne considere que les premieres couches instables
  DO k = nlay - 2, 1, -1
    DO ig = 1, ngrid
      IF (ztv(ig,k)>ztv(ig,k+1) .AND. ztv(ig,k+1)<=ztv(ig,k+2)) THEN
        lentr(ig) = k
      END IF
    END DO
  END DO

  ! determination du lmin: couche d ou provient le thermique
  DO ig = 1, ngrid
    lmin(ig) = 1
  END DO
  DO ig = 1, ngrid
    DO l = nlay, 2, -1
      IF (ztv(ig,l-1)>ztv(ig,l)) THEN
        lmin(ig) = l - 1
      END IF
    END DO
  END DO

  ! definition de l'entrainement des couches
  DO l = 1, klev - 1
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. l>=lmin(ig) .AND. l<=lentr(ig)) THEN
        entr_star(ig, l) = (ztv(ig,l)-ztv(ig,l+1))** & ! s
                                                       ! (zlev(ig,l+1)-zlev(ig,l))
          sqrt(zlev(ig,l+1))
      END IF
    END DO
  END DO
  ! pas de thermique si couche 1 stable
  DO ig = 1, ngrid
    IF (lmin(ig)>1) THEN
      DO l = 1, klev
        entr_star(ig, l) = 0.
      END DO
    END IF
  END DO
  ! calcul de l entrainement total
  DO ig = 1, ngrid
    entr_star_tot(ig) = 0.
  END DO
  DO ig = 1, ngrid
    DO k = 1, klev
      entr_star_tot(ig) = entr_star_tot(ig) + entr_star(ig, k)
    END DO
  END DO

  ! print*,'fin calcul entr_star'
  DO k = 1, klev
    DO ig = 1, ngrid
      ztva(ig, k) = ztv(ig, k)
    END DO
  END DO
  ! RC
  ! print*,'7 OK convect8'
  DO k = 1, klev + 1
    DO ig = 1, ngrid
      zw2(ig, k) = 0.
      fmc(ig, k) = 0.
      ! CR
      f_star(ig, k) = 0.
      ! RC
      larg_cons(ig, k) = 0.
      larg_detr(ig, k) = 0.
      wa_moy(ig, k) = 0.
    END DO
  END DO

  ! print*,'8 OK convect8'
  DO ig = 1, ngrid
    linter(ig) = 1.
    lmaxa(ig) = 1
    lmix(ig) = 1
    wmaxa(ig) = 0.
  END DO

  ! CR:
  DO l = 1, nlay - 2
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. entr_star(ig,l)>1.E-10 .AND. &
          zw2(ig,l)<1E-10) THEN
        f_star(ig, l+1) = entr_star(ig, l)
        ! test:calcul de dteta
        zw2(ig, l+1) = 2.*rg*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig, l+1)* &
          (zlev(ig,l+1)-zlev(ig,l))*0.4*pphi(ig, l)/(pphi(ig,l+1)-pphi(ig,l))
        larg_detr(ig, l) = 0.
      ELSE IF ((zw2(ig,l)>=1E-10) .AND. (f_star(ig,l)+entr_star(ig, &
          l)>1.E-10)) THEN
        f_star(ig, l+1) = f_star(ig, l) + entr_star(ig, l)
        ztva(ig, l) = (f_star(ig,l)*ztva(ig,l-1)+entr_star(ig,l)*ztv(ig,l))/ &
          f_star(ig, l+1)
        zw2(ig, l+1) = zw2(ig, l)*(f_star(ig,l)/f_star(ig,l+1))**2 + &
          2.*rg*(ztva(ig,l)-ztv(ig,l))/ztv(ig, l)*(zlev(ig,l+1)-zlev(ig,l))
      END IF
      ! determination de zmax continu par interpolation lineaire
      IF (zw2(ig,l+1)<0.) THEN
        ! test
        IF (abs(zw2(ig,l+1)-zw2(ig,l))<1E-10) THEN
          ! print*,'pb linter'
        END IF
        linter(ig) = (l*(zw2(ig,l+1)-zw2(ig,l))-zw2(ig,l))/(zw2(ig,l+1)-zw2( &
          ig,l))
        zw2(ig, l+1) = 0.
        lmaxa(ig) = l
      ELSE
        IF (zw2(ig,l+1)<0.) THEN
          ! print*,'pb1 zw2<0'
        END IF
        wa_moy(ig, l+1) = sqrt(zw2(ig,l+1))
      END IF
      IF (wa_moy(ig,l+1)>wmaxa(ig)) THEN
        ! lmix est le niveau de la couche ou w (wa_moy) est maximum
        lmix(ig) = l + 1
        wmaxa(ig) = wa_moy(ig, l+1)
      END IF
    END DO
  END DO
  ! print*,'fin calcul zw2'

  ! Calcul de la couche correspondant a la hauteur du thermique
  DO ig = 1, ngrid
    lmax(ig) = lentr(ig)
  END DO
  DO ig = 1, ngrid
    DO l = nlay, lentr(ig) + 1, -1
      IF (zw2(ig,l)<=1.E-10) THEN
        lmax(ig) = l - 1
      END IF
    END DO
  END DO
  ! pas de thermique si couche 1 stable
  DO ig = 1, ngrid
    IF (lmin(ig)>1) THEN
      lmax(ig) = 1
      lmin(ig) = 1
    END IF
  END DO

  ! Determination de zw2 max
  DO ig = 1, ngrid
    wmax(ig) = 0.
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (l<=lmax(ig)) THEN
        IF (zw2(ig,l)<0.) THEN
          ! print*,'pb2 zw2<0'
        END IF
        zw2(ig, l) = sqrt(zw2(ig,l))
        wmax(ig) = max(wmax(ig), zw2(ig,l))
      ELSE
        zw2(ig, l) = 0.
      END IF
    END DO
  END DO

  ! Longueur caracteristique correspondant a la hauteur des thermiques.
  DO ig = 1, ngrid
    zmax(ig) = 0.
    zlevinter(ig) = zlev(ig, 1)
  END DO
  DO ig = 1, ngrid
    ! calcul de zlevinter
    zlevinter(ig) = (zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*linter(ig) + &
      zlev(ig, lmax(ig)) - lmax(ig)*(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))
    zmax(ig) = max(zmax(ig), zlevinter(ig)-zlev(ig,lmin(ig)))
  END DO

  ! print*,'avant fermeture'
  ! Fermeture,determination de f
  DO ig = 1, ngrid
    entr_star2(ig) = 0.
  END DO
  DO ig = 1, ngrid
    IF (entr_star_tot(ig)<1.E-10) THEN
      f(ig) = 0.
    ELSE
      DO k = lmin(ig), lentr(ig)
        entr_star2(ig) = entr_star2(ig) + entr_star(ig, k)**2/(rho(ig,k)*( &
          zlev(ig,k+1)-zlev(ig,k)))
      END DO
      ! Nouvelle fermeture
      f(ig) = wmax(ig)/(max(500.,zmax(ig))*r_aspect*entr_star2(ig))* &
        entr_star_tot(ig)
      ! test
      ! if (first) then
      ! f(ig)=f(ig)+(f0(ig)-f(ig))*exp(-ptimestep/zmax(ig)
      ! s             *wmax(ig))
      ! endif
    END IF
    ! f0(ig)=f(ig)
    ! first=.true.
  END DO
  ! print*,'apres fermeture'

  ! Calcul de l'entrainement
  DO k = 1, klev
    DO ig = 1, ngrid
      entr(ig, k) = f(ig)*entr_star(ig, k)
    END DO
  END DO
  ! CR:test pour entrainer moins que la masse
  DO ig = 1, ngrid
    DO l = 1, lentr(ig)
      IF ((entr(ig,l)*ptimestep)>(0.9*masse(ig,l))) THEN
        entr(ig, l+1) = entr(ig, l+1) + entr(ig, l) - &
          0.9*masse(ig, l)/ptimestep
        entr(ig, l) = 0.9*masse(ig, l)/ptimestep
      END IF
    END DO
  END DO
  ! CR: fin test
  ! Calcul des flux
  DO ig = 1, ngrid
    DO l = 1, lmax(ig) - 1
      fmc(ig, l+1) = fmc(ig, l) + entr(ig, l)
    END DO
  END DO

  ! RC


  ! print*,'9 OK convect8'
  ! print*,'WA1 ',wa_moy

  ! determination de l'indice du debut de la mixed layer ou w decroit

  ! calcul de la largeur de chaque ascendance dans le cas conservatif.
  ! dans ce cas simple, on suppose que la largeur de l'ascendance provenant
  ! d'une couche est égale à la hauteur de la couche alimentante.
  ! La vitesse maximale dans l'ascendance est aussi prise comme estimation
  ! de la vitesse d'entrainement horizontal dans la couche alimentante.

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (l<=lmaxa(ig)) THEN
        zw = max(wa_moy(ig,l), 1.E-10)
        larg_cons(ig, l) = zmax(ig)*r_aspect*fmc(ig, l)/(rhobarz(ig,l)*zw)
      END IF
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (l<=lmaxa(ig)) THEN
        ! if (idetr.eq.0) then
        ! cette option est finalement en dur.
        IF ((l_mix*zlev(ig,l))<0.) THEN
          ! print*,'pb l_mix*zlev<0'
        END IF
        ! CR: test: nouvelle def de lambda
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        IF (zw2(ig,l)>1.E-10) THEN
          larg_detr(ig, l) = sqrt((l_mix/zw2(ig,l))*zlev(ig,l))
        ELSE
          larg_detr(ig, l) = sqrt(l_mix*zlev(ig,l))
        END IF
        ! RC
        ! else if (idetr.eq.1) then
        ! larg_detr(ig,l)=larg_cons(ig,l)
        ! s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
        ! else if (idetr.eq.2) then
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        ! s            *sqrt(wa_moy(ig,l))
        ! else if (idetr.eq.4) then
        ! larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
        ! s            *wa_moy(ig,l)
        ! endif
      END IF
    END DO
  END DO

  ! print*,'10 OK convect8'
  ! print*,'WA2 ',wa_moy
  ! calcul de la fraction de la maille concernée par l'ascendance en tenant
  ! compte de l'epluchage du thermique.

  ! CR def de  zmix continu (profil parabolique des vitesses)
  DO ig = 1, ngrid
    IF (lmix(ig)>1.) THEN
      ! test
      IF (((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig)))- &
          (zlev(ig,lmix(ig)+1)))-(zw2(ig,lmix(ig))- &
          zw2(ig,lmix(ig)+1))*((zlev(ig,lmix(ig)-1))- &
          (zlev(ig,lmix(ig)))))>1E-10) THEN

        zmix(ig) = ((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig)) &
          )**2-(zlev(ig,lmix(ig)+1))**2)-(zw2(ig,lmix(ig))-zw2(ig, &
          lmix(ig)+1))*((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))/ &
          (2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))*((zlev(ig,lmix(ig)))- &
          (zlev(ig,lmix(ig)+1)))-(zw2(ig,lmix(ig))- &
          zw2(ig,lmix(ig)+1))*((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
      ELSE
        zmix(ig) = zlev(ig, lmix(ig))
        ! print*,'pb zmix'
      END IF
    ELSE
      zmix(ig) = 0.
    END IF
    ! test
    IF ((zmax(ig)-zmix(ig))<0.) THEN
      zmix(ig) = 0.99*zmax(ig)
      ! print*,'pb zmix>zmax'
    END IF
  END DO

  ! calcul du nouveau lmix correspondant
  DO ig = 1, ngrid
    DO l = 1, klev
      IF (zmix(ig)>=zlev(ig,l) .AND. zmix(ig)<zlev(ig,l+1)) THEN
        lmix(ig) = l
      END IF
    END DO
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (larg_cons(ig,l)>1.) THEN
        ! print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
        fraca(ig, l) = (larg_cons(ig,l)-larg_detr(ig,l))/(r_aspect*zmax(ig))
        ! test
        fraca(ig, l) = max(fraca(ig,l), 0.)
        fraca(ig, l) = min(fraca(ig,l), 0.5)
        fracd(ig, l) = 1. - fraca(ig, l)
        fracc(ig, l) = larg_cons(ig, l)/(r_aspect*zmax(ig))
      ELSE
        ! wa_moy(ig,l)=0.
        fraca(ig, l) = 0.
        fracc(ig, l) = 0.
        fracd(ig, l) = 1.
      END IF
    END DO
  END DO
  ! CR: calcul de fracazmix
  DO ig = 1, ngrid
    fracazmix(ig) = (fraca(ig,lmix(ig)+1)-fraca(ig,lmix(ig)))/ &
      (zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))*zmix(ig) + &
      fraca(ig, lmix(ig)) - zlev(ig, lmix(ig))*(fraca(ig,lmix(ig)+1)-fraca(ig &
      ,lmix(ig)))/(zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      IF (larg_cons(ig,l)>1.) THEN
        IF (l>lmix(ig)) THEN
          ! test
          IF (zmax(ig)-zmix(ig)<1.E-10) THEN
            ! print*,'pb xxx'
            xxx(ig, l) = (lmaxa(ig)+1.-l)/(lmaxa(ig)+1.-lmix(ig))
          ELSE
            xxx(ig, l) = (zmax(ig)-zlev(ig,l))/(zmax(ig)-zmix(ig))
          END IF
          IF (idetr==0) THEN
            fraca(ig, l) = fracazmix(ig)
          ELSE IF (idetr==1) THEN
            fraca(ig, l) = fracazmix(ig)*xxx(ig, l)
          ELSE IF (idetr==2) THEN
            fraca(ig, l) = fracazmix(ig)*(1.-(1.-xxx(ig,l))**2)
          ELSE
            fraca(ig, l) = fracazmix(ig)*xxx(ig, l)**2
          END IF
          ! print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
          fraca(ig, l) = max(fraca(ig,l), 0.)
          fraca(ig, l) = min(fraca(ig,l), 0.5)
          fracd(ig, l) = 1. - fraca(ig, l)
          fracc(ig, l) = larg_cons(ig, l)/(r_aspect*zmax(ig))
        END IF
      END IF
    END DO
  END DO

  ! print*,'fin calcul fraca'
  ! print*,'11 OK convect8'
  ! print*,'Ea3 ',wa_moy
  ! ------------------------------------------------------------------
  ! Calcul de fracd, wd
  ! somme wa - wd = 0
  ! ------------------------------------------------------------------


  DO ig = 1, ngrid
    fm(ig, 1) = 0.
    fm(ig, nlay+1) = 0.
  END DO

  DO l = 2, nlay
    DO ig = 1, ngrid
      fm(ig, l) = fraca(ig, l)*wa_moy(ig, l)*rhobarz(ig, l)
      ! CR:test
      IF (entr(ig,l-1)<1E-10 .AND. fm(ig,l)>fm(ig,l-1) .AND. l>lmix(ig)) THEN
        fm(ig, l) = fm(ig, l-1)
        ! write(1,*)'ajustement fm, l',l
      END IF
      ! write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
      ! RC
    END DO
    DO ig = 1, ngrid
      IF (fracd(ig,l)<0.1) THEN
        abort_message = 'fracd trop petit'
        CALL abort_physic(modname, abort_message, 1)
      ELSE
        ! vitesse descendante "diagnostique"
        wd(ig, l) = fm(ig, l)/(fracd(ig,l)*rhobarz(ig,l))
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      ! masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
      masse(ig, l) = (pplev(ig,l)-pplev(ig,l+1))/rg
    END DO
  END DO

  ! print*,'12 OK convect8'
  ! print*,'WA4 ',wa_moy
  ! c------------------------------------------------------------------
  ! calcul du transport vertical
  ! ------------------------------------------------------------------

  GO TO 4444
  ! print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
  DO l = 2, nlay - 1
    DO ig = 1, ngrid
      IF (fm(ig,l+1)*ptimestep>masse(ig,l) .AND. fm(ig,l+1)*ptimestep>masse( &
          ig,l+1)) THEN
        ! print*,'WARN!!! FM>M ig=',ig,' l=',l,'  FM='
        ! s         ,fm(ig,l+1)*ptimestep
        ! s         ,'   M=',masse(ig,l),masse(ig,l+1)
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (entr(ig,l)*ptimestep>masse(ig,l)) THEN
        ! print*,'WARN!!! E>M ig=',ig,' l=',l,'  E=='
        ! s         ,entr(ig,l)*ptimestep
        ! s         ,'   M=',masse(ig,l)
      END IF
    END DO
  END DO

  DO l = 1, nlay
    DO ig = 1, ngrid
      IF (.NOT. fm(ig,l)>=0. .OR. .NOT. fm(ig,l)<=10.) THEN
        ! print*,'WARN!!! fm exagere ig=',ig,'   l=',l
        ! s         ,'   FM=',fm(ig,l)
      END IF
      IF (.NOT. masse(ig,l)>=1.E-10 .OR. .NOT. masse(ig,l)<=1.E4) THEN
        ! print*,'WARN!!! masse exagere ig=',ig,'   l=',l
        ! s         ,'   M=',masse(ig,l)
        ! print*,'rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)',
        ! s                 rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)
        ! print*,'zlev(ig,l+1),zlev(ig,l)'
        ! s                ,zlev(ig,l+1),zlev(ig,l)
        ! print*,'pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)'
        ! s                ,pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)
      END IF
      IF (.NOT. entr(ig,l)>=0. .OR. .NOT. entr(ig,l)<=10.) THEN
        ! print*,'WARN!!! entr exagere ig=',ig,'   l=',l
        ! s         ,'   E=',entr(ig,l)
      END IF
    END DO
  END DO

4444 CONTINUE

  ! CR:redefinition du entr
  DO l = 1, nlay
    DO ig = 1, ngrid
      detr(ig, l) = fm(ig, l) + entr(ig, l) - fm(ig, l+1)
      IF (detr(ig,l)<0.) THEN
        entr(ig, l) = entr(ig, l) - detr(ig, l)
        detr(ig, l) = 0.
        ! print*,'WARNING !!! detrainement negatif ',ig,l
      END IF
    END DO
  END DO
  ! RC
  IF (w2di==1) THEN
    fm0 = fm0 + ptimestep*(fm-fm0)/tho
    entr0 = entr0 + ptimestep*(entr-entr0)/tho
  ELSE
    fm0 = fm
    entr0 = entr
  END IF

  IF (1==1) THEN
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zh, zdhadj, &
      zha)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zo, pdoadj, &
      zoa)
  ELSE
    CALL dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zh, &
      zdhadj, zha)
    CALL dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zo, &
      pdoadj, zoa)
  END IF

  IF (1==0) THEN
    CALL dvthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca, zmax, &
      zu, zv, pduadj, pdvadj, zua, zva)
  ELSE
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zu, pduadj, &
      zua)
    CALL dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse, zv, pdvadj, &
      zva)
  END IF

  DO l = 1, nlay
    DO ig = 1, ngrid
      zf = 0.5*(fracc(ig,l)+fracc(ig,l+1))
      zf2 = zf/(1.-zf)
      thetath2(ig, l) = zf2*(zha(ig,l)-zh(ig,l))**2
      wth2(ig, l) = zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
    END DO
  END DO



  ! print*,'13 OK convect8'
  ! print*,'WA5 ',wa_moy
  DO l = 1, nlay
    DO ig = 1, ngrid
      pdtadj(ig, l) = zdhadj(ig, l)*zpspsk(ig, l)
    END DO
  END DO


  ! do l=1,nlay
  ! do ig=1,ngrid
  ! if(abs(pdtadj(ig,l))*86400..gt.500.) then
  ! print*,'WARN!!! ig=',ig,'  l=',l
  ! s         ,'   pdtadj=',pdtadj(ig,l)
  ! endif
  ! if(abs(pdoadj(ig,l))*86400..gt.1.) then
  ! print*,'WARN!!! ig=',ig,'  l=',l
  ! s         ,'   pdoadj=',pdoadj(ig,l)
  ! endif
  ! enddo
  ! enddo

  ! print*,'14 OK convect8'
  ! ------------------------------------------------------------------
  ! Calculs pour les sorties
  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE thermcell_sec

