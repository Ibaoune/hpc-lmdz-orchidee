
! $Id: thermcell.F90 3102 2017-12-03 20:27:42Z oboucher $

SUBROUTINE calcul_sec(ngrid, nlay, ptimestep, pplay, pplev, pphi, zlev, pu, &
    pv, pt, po, zmax, wmax, zw2, lmix & ! s
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
  REAL zmax(klon), zw, zw2(klon, klev+1), ztva(klon, klev)

  REAL zlev(klon, klev+1), zlay(klon, klev)
  REAL zh(klon, klev), zdhadj(klon, klev)
  REAL ztv(klon, klev)
  REAL zu(klon, klev), zv(klon, klev), zo(klon, klev)
  REAL wh(klon, klev+1)
  REAL wu(klon, klev+1), wv(klon, klev+1), wo(klon, klev+1)
  REAL zla(klon, klev+1)
  REAL zwa(klon, klev+1)
  REAL zld(klon, klev+1)
  ! real zwd(klon,klev+1)
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
  ! integer isplit,nsplit
  INTEGER isplit, nsplit, ialt
  PARAMETER (nsplit=10)
  DATA isplit/0/
  SAVE isplit
  !$OMP THREADPRIVATE(isplit)

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
  REAL zalim(klon)
  INTEGER lalim(klon)
  REAL norme(klon)
  REAL f(klon), f0(klon)
  REAL zlevinter(klon)
  LOGICAL therm
  LOGICAL first
  DATA first/.FALSE./
  SAVE first
  !$OMP THREADPRIVATE(first)
  ! RC

  CHARACTER *2 str2
  CHARACTER *10 str10

  CHARACTER (LEN=20) :: modname = 'calcul_sec'
  CHARACTER (LEN=80) :: abort_message


  ! LOGICAL vtest(klon),down

  EXTERNAL scopy

  INTEGER ncorrec
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
  therm = .FALSE.
  DO k = nlay - 2, 1, -1
    DO ig = 1, ngrid
      IF (ztv(ig,k)>ztv(ig,k+1) .AND. ztv(ig,k+1)<=ztv(ig,k+2)) THEN
        lentr(ig) = k + 1
        therm = .TRUE.
      END IF
    END DO
  END DO
  ! limitation de la valeur du lentr
  ! do ig=1,ngrid
  ! lentr(ig)=min(5,lentr(ig))
  ! enddo
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
  ! initialisations
  DO ig = 1, ngrid
    zalim(ig) = 0.
    norme(ig) = 0.
    lalim(ig) = 1
  END DO
  DO k = 1, klev - 1
    DO ig = 1, ngrid
      zalim(ig) = zalim(ig) + zlev(ig, k)*max(0., (ztv(ig,k)-ztv(ig, &
        k+1))/(zlev(ig,k+1)-zlev(ig,k)))
      ! s         *(zlev(ig,k+1)-zlev(ig,k))
      norme(ig) = norme(ig) + max(0., (ztv(ig,k)-ztv(ig,k+1))/(zlev(ig, &
        k+1)-zlev(ig,k)))
      ! s          *(zlev(ig,k+1)-zlev(ig,k))
    END DO
  END DO
  DO ig = 1, ngrid
    IF (norme(ig)>1.E-10) THEN
      zalim(ig) = max(10.*zalim(ig)/norme(ig), zlev(ig,2))
      ! zalim(ig)=min(zalim(ig),zlev(ig,lentr(ig)))
    END IF
  END DO
  ! détermination du lalim correspondant
  DO k = 1, klev - 1
    DO ig = 1, ngrid
      IF ((zalim(ig)>zlev(ig,k)) .AND. (zalim(ig)<=zlev(ig,k+1))) THEN
        lalim(ig) = k
      END IF
    END DO
  END DO

  ! definition de l'entrainement des couches
  DO l = 1, klev - 1
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. l>=lmin(ig) .AND. l<lentr(ig)) THEN
        entr_star(ig, l) = max((ztv(ig,l)-ztv(ig,l+1)), 0.) & ! s
                                                              ! *(zlev(ig,l+1)-zlev(ig,l))
          *sqrt(zlev(ig,l+1))
        ! autre def
        ! entr_star(ig,l)=zlev(ig,l+1)*(1.-(zlev(ig,l+1)
        ! s                         /zlev(ig,lentr(ig)+2)))**(3./2.)
      END IF
    END DO
  END DO
  ! nouveau test
  ! if (therm) then
  DO l = 1, klev - 1
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. l>=lmin(ig) .AND. l<=lalim(ig) .AND. &
          zalim(ig)>1.E-10) THEN
        ! if (l.le.lentr(ig)) then
        ! entr_star(ig,l)=zlev(ig,l+1)*(1.-(zlev(ig,l+1)
        ! s                         /zalim(ig)))**(3./2.)
        ! write(10,*)zlev(ig,l),entr_star(ig,l)
      END IF
    END DO
  END DO
  ! endif
  ! pas de thermique si couche 1 stable
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
  ! Calcul entrainement normalise
  DO ig = 1, ngrid
    IF (entr_star_tot(ig)>1.E-10) THEN
      ! do l=1,lentr(ig)
      DO l = 1, klev
        ! def possibles pour entr_star: zdthetadz, dthetadz, zdtheta
        entr_star(ig, l) = entr_star(ig, l)/entr_star_tot(ig)
      END DO
    END IF
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
    ! lmax(ig)=lalim(ig)
  END DO
  DO ig = 1, ngrid
    DO l = nlay, lentr(ig) + 1, -1
      ! do l=nlay,lalim(ig)+1,-1
      IF (zw2(ig,l)<=1.E-10) THEN
        lmax(ig) = l - 1
      END IF
    END DO
  END DO
  ! pas de thermique si couche 1 stable
  DO ig = 1, ngrid
    IF (lmin(ig)>5) THEN
      lmax(ig) = 1
      lmin(ig) = 1
      lentr(ig) = 1
      lalim(ig) = 1
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
  DO ig = 1, ngrid
    ! write(8,*)zmax(ig),lmax(ig),lentr(ig),lmin(ig)
  END DO
  ! on stope après les calculs de zmax et wmax
  RETURN

  ! print*,'avant fermeture'
  ! Fermeture,determination de f
  ! Attention! entrainement normalisé ou pas?
  DO ig = 1, ngrid
    entr_star2(ig) = 0.
  END DO
  DO ig = 1, ngrid
    IF (entr_star_tot(ig)<1.E-10) THEN
      f(ig) = 0.
    ELSE
      DO k = lmin(ig), lentr(ig)
        ! do k=lmin(ig),lalim(ig)
        entr_star2(ig) = entr_star2(ig) + entr_star(ig, k)**2/(rho(ig,k)*( &
          zlev(ig,k+1)-zlev(ig,k)))
      END DO
      ! Nouvelle fermeture
      f(ig) = wmax(ig)/(max(500.,zmax(ig))*r_aspect*entr_star2(ig))
      ! s            *entr_star_tot(ig)
      ! test
      ! if (first) then
      f(ig) = f(ig) + (f0(ig)-f(ig))*exp(-ptimestep/zmax(ig)*wmax(ig))
      ! endif
    END IF
    f0(ig) = f(ig)
    ! first=.true.
  END DO
  ! print*,'apres fermeture'
  ! on stoppe après la fermeture
  RETURN
  ! Calcul de l'entrainement
  DO k = 1, klev
    DO ig = 1, ngrid
      entr(ig, k) = f(ig)*entr_star(ig, k)
    END DO
  END DO
  ! on stoppe après le calcul de entr
  ! RETURN
  ! CR:test pour entrainer moins que la masse
  ! do ig=1,ngrid
  ! do l=1,lentr(ig)
  ! if ((entr(ig,l)*ptimestep).gt.(0.9*masse(ig,l))) then
  ! entr(ig,l+1)=entr(ig,l+1)+entr(ig,l)
  ! s                       -0.9*masse(ig,l)/ptimestep
  ! entr(ig,l)=0.9*masse(ig,l)/ptimestep
  ! endif
  ! enddo
  ! enddo
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
        ! entr(ig,l)=entr(ig,l)-detr(ig,l)
        fm(ig, l+1) = fm(ig, l) + entr(ig, l)
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

    isplit = isplit + 1


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
END SUBROUTINE calcul_sec

SUBROUTINE fermeture_seche(ngrid, nlay, pplay, pplev, pphi, zlev, rhobarz, &
    f0, zpspsk, alim_star, zh, zo, lentr, lmin, nu_min, nu_max, r_aspect, &
    zmax, wmax)

  USE dimphy
  IMPLICIT NONE

  include "YOMCST.h"

  INTEGER ngrid, nlay
  REAL pplay(ngrid, nlay), pplev(ngrid, nlay+1)
  REAL pphi(ngrid, nlay)
  REAL zlev(klon, klev+1)
  REAL alim_star(klon, klev)
  REAL f0(klon)
  INTEGER lentr(klon)
  INTEGER lmin(klon)
  REAL zmax(klon)
  REAL wmax(klon)
  REAL nu_min
  REAL nu_max
  REAL r_aspect
  REAL rhobarz(klon, klev+1)
  REAL zh(klon, klev)
  REAL zo(klon, klev)
  REAL zpspsk(klon, klev)

  INTEGER ig, l

  REAL f_star(klon, klev+1)
  REAL detr_star(klon, klev)
  REAL entr_star(klon, klev)
  REAL zw2(klon, klev+1)
  REAL linter(klon)
  INTEGER lmix(klon)
  INTEGER lmax(klon)
  REAL zlevinter(klon)
  REAL wa_moy(klon, klev+1)
  REAL wmaxa(klon)
  REAL ztv(klon, klev)
  REAL ztva(klon, klev)
  REAL nu(klon, klev)
  ! real zmax0_sec(klon)
  ! save zmax0_sec
  REAL, SAVE, ALLOCATABLE :: zmax0_sec(:)
  !$OMP THREADPRIVATE(zmax0_sec)
  LOGICAL, SAVE :: first = .TRUE.
  !$OMP THREADPRIVATE(first)

  IF (first) THEN
    ALLOCATE (zmax0_sec(klon))
    first = .FALSE.
  END IF

  DO l = 1, nlay
    DO ig = 1, ngrid
      ztv(ig, l) = zh(ig, l)/zpspsk(ig, l)
      ztv(ig, l) = ztv(ig, l)*(1.+retv*zo(ig,l))
    END DO
  END DO
  DO l = 1, nlay - 2
    DO ig = 1, ngrid
      IF (ztv(ig,l)>ztv(ig,l+1) .AND. alim_star(ig,l)>1.E-10 .AND. &
          zw2(ig,l)<1E-10) THEN
        f_star(ig, l+1) = alim_star(ig, l)
        ! test:calcul de dteta
        zw2(ig, l+1) = 2.*rg*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig, l+1)* &
          (zlev(ig,l+1)-zlev(ig,l))*0.4*pphi(ig, l)/(pphi(ig,l+1)-pphi(ig,l))
      ELSE IF ((zw2(ig,l)>=1E-10) .AND. (f_star(ig,l)+alim_star(ig, &
          l))>1.E-10) THEN
        ! estimation du detrainement a partir de la geometrie du pas
        ! precedent
        ! tests sur la definition du detr
        nu(ig, l) = (nu_min+nu_max)/2.*(1.-(nu_max-nu_min)/(nu_max+nu_min)* &
          tanh((((ztva(ig,l-1)-ztv(ig,l))/ztv(ig,l))/0.0005)))

        detr_star(ig, l) = rhobarz(ig, l)*sqrt(zw2(ig,l))/ &
          (r_aspect*zmax0_sec(ig))* & ! s
                                      ! /(r_aspect*zmax0(ig))*
          (sqrt(nu(ig,l)*zlev(ig,l+1)/sqrt(zw2(ig,l)))-sqrt(nu(ig,l)*zlev(ig, &
          l)/sqrt(zw2(ig,l))))
        detr_star(ig, l) = detr_star(ig, l)/f0(ig)
        IF ((detr_star(ig,l))>f_star(ig,l)) THEN
          detr_star(ig, l) = f_star(ig, l)
        END IF
        entr_star(ig, l) = 0.9*detr_star(ig, l)
        IF ((l<lentr(ig))) THEN
          entr_star(ig, l) = 0.
          ! detr_star(ig,l)=0.
        END IF
        ! print*,'ok detr_star'
        ! prise en compte du detrainement dans le calcul du flux
        f_star(ig, l+1) = f_star(ig, l) + alim_star(ig, l) + &
          entr_star(ig, l) - detr_star(ig, l)
        ! test sur le signe de f_star
        IF ((f_star(ig,l+1)+detr_star(ig,l))>1.E-10) THEN
          ! AM on melange Tl et qt du thermique
          ztva(ig, l) = (f_star(ig,l)*ztva(ig,l-1)+(entr_star(ig, &
            l)+alim_star(ig,l))*ztv(ig,l))/(f_star(ig,l+1)+detr_star(ig,l))
          zw2(ig, l+1) = zw2(ig, l)*(f_star(ig,l)/(f_star(ig, &
            l+1)+detr_star(ig,l)))**2 + 2.*rg*(ztva(ig,l)-ztv(ig,l))/ztv(ig, &
            l)*(zlev(ig,l+1)-zlev(ig,l))
        END IF
      END IF

      IF (zw2(ig,l+1)<0.) THEN
        linter(ig) = (l*(zw2(ig,l+1)-zw2(ig,l))-zw2(ig,l))/(zw2(ig,l+1)-zw2( &
          ig,l))
        zw2(ig, l+1) = 0.
        ! print*,'linter=',linter(ig)
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
    ! pour le cas ou on prend tjs lmin=1
    ! zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
    zmax(ig) = max(zmax(ig), zlevinter(ig)-zlev(ig,1))
    zmax0_sec(ig) = zmax(ig)
  END DO

  RETURN
END SUBROUTINE fermeture_seche
