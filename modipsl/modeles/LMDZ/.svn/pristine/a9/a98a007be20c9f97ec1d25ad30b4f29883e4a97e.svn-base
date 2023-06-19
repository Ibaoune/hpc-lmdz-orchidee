
! $Header$

SUBROUTINE vdif_kcay(ngrid, dt, g, rconst, plev, temp, zlev, zlay, u, v, &
    teta, cd, q2, q2diag, km, kn, ustar, l_mix)
  USE dimphy
  IMPLICIT NONE

  ! dt : pas de temps
  ! g  : g
  ! zlev : altitude a chaque niveau (interface inferieure de la couche
  ! de meme indice)
  ! zlay : altitude au centre de chaque couche
  ! u,v : vitesse au centre de chaque couche
  ! (en entree : la valeur au debut du pas de temps)
  ! teta : temperature potentielle au centre de chaque couche
  ! (en entree : la valeur au debut du pas de temps)
  ! cd : cdrag
  ! (en entree : la valeur au debut du pas de temps)
  ! q2 : $q^2$ au bas de chaque couche
  ! (en entree : la valeur au debut du pas de temps)
  ! (en sortie : la valeur a la fin du pas de temps)
  ! km : diffusivite turbulente de quantite de mouvement (au bas de chaque
  ! couche)
  ! (en sortie : la valeur a la fin du pas de temps)
  ! kn : diffusivite turbulente des scalaires (au bas de chaque couche)
  ! (en sortie : la valeur a la fin du pas de temps)

  ! .......................................................................
  REAL dt, g, rconst
  REAL plev(klon, klev+1), temp(klon, klev)
  REAL ustar(klon), snstable
  REAL zlev(klon, klev+1)
  REAL zlay(klon, klev)
  REAL u(klon, klev)
  REAL v(klon, klev)
  REAL teta(klon, klev)
  REAL cd(klon)
  REAL q2(klon, klev+1), q2s(klon, klev+1)
  REAL q2diag(klon, klev+1)
  REAL km(klon, klev+1)
  REAL kn(klon, klev+1)
  REAL sq(klon), sqz(klon), zz(klon, klev+1), zq, long0(klon)

  INTEGER l_mix, iii
  ! .......................................................................

  ! nlay : nombre de couches
  ! nlev : nombre de niveaux
  ! ngrid : nombre de points de grille
  ! unsdz : 1 sur l'epaisseur de couche
  ! unsdzdec : 1 sur la distance entre le centre de la couche et le
  ! centre de la couche inferieure
  ! q : echelle de vitesse au bas de chaque couche
  ! (valeur a la fin du pas de temps)

  ! .......................................................................
  INTEGER nlay, nlev, ngrid
  REAL unsdz(klon, klev)
  REAL unsdzdec(klon, klev+1)
  REAL q(klon, klev+1)

  ! .......................................................................

  ! kmpre : km au debut du pas de temps
  ! qcstat : q : solution stationnaire du probleme couple
  ! (valeur a la fin du pas de temps)
  ! q2cstat : q2 : solution stationnaire du probleme couple
  ! (valeur a la fin du pas de temps)

  ! .......................................................................
  REAL kmpre(klon, klev+1)
  REAL qcstat
  REAL q2cstat
  REAL sss, sssq
  ! .......................................................................

  ! long : longueur de melange calculee selon Blackadar

  ! .......................................................................
  REAL long(klon, klev+1)
  ! .......................................................................

  ! kmq3 : terme en q^3 dans le developpement de km
  ! (valeur au debut du pas de temps)
  ! kmcstat : valeur de km solution stationnaire du systeme {q2 ; du/dz}
  ! (valeur a la fin du pas de temps)
  ! knq3 : terme en q^3 dans le developpement de kn
  ! mcstat : valeur de m solution stationnaire du systeme {q2 ; du/dz}
  ! (valeur a la fin du pas de temps)
  ! m2cstat : valeur de m2 solution stationnaire du systeme {q2 ; du/dz}
  ! (valeur a la fin du pas de temps)
  ! m : valeur a la fin du pas de temps
  ! mpre : valeur au debut du pas de temps
  ! m2 : valeur a la fin du pas de temps
  ! n2 : valeur a la fin du pas de temps

  ! .......................................................................
  REAL kmq3
  REAL kmcstat
  REAL knq3
  REAL mcstat
  REAL m2cstat
  REAL m(klon, klev+1)
  REAL mpre(klon, klev+1)
  REAL m2(klon, klev+1)
  REAL n2(klon, klev+1)
  ! .......................................................................

  ! gn : intermediaire pour les coefficients de stabilite
  ! gnmin : borne inferieure de gn (-0.23 ou -0.28)
  ! gnmax : borne superieure de gn (0.0233)
  ! gninf : vrai si gn est en dessous de sa borne inferieure
  ! gnsup : vrai si gn est en dessus de sa borne superieure
  ! gm : drole d'objet bien utile
  ! ri : nombre de Richardson
  ! sn : coefficient de stabilite pour n
  ! snq2 : premier terme du developement limite de sn en q2
  ! sm : coefficient de stabilite pour m
  ! smq2 : premier terme du developement limite de sm en q2

  ! .......................................................................
  REAL gn
  REAL gnmin
  REAL gnmax
  LOGICAL gninf
  LOGICAL gnsup
  REAL gm
  ! REAL ri(klon,klev+1)
  REAL sn(klon, klev+1)
  REAL snq2(klon, klev+1)
  REAL sm(klon, klev+1)
  REAL smq2(klon, klev+1)
  ! .......................................................................

  ! kappa : consatnte de Von Karman (0.4)
  ! long00 : longueur de reference pour le calcul de long (160)
  ! a1,a2,b1,b2,c1 : constantes d'origine pour les  coefficients
  ! de stabilite (0.92/0.74/16.6/10.1/0.08)
  ! cn1,cn2 : constantes pour sn
  ! cm1,cm2,cm3,cm4 : constantes pour sm

  ! .......................................................................
  REAL kappa
  REAL long00
  REAL a1, a2, b1, b2, c1
  REAL cn1, cn2
  REAL cm1, cm2, cm3, cm4
  ! .......................................................................

  ! termq : termes en $q$ dans l'equation de q2
  ! termq3 : termes en $q^3$ dans l'equation de q2
  ! termqm2 : termes en $q*m^2$ dans l'equation de q2
  ! termq3m2 : termes en $q^3*m^2$ dans l'equation de q2

  ! .......................................................................
  REAL termq
  REAL termq3
  REAL termqm2
  REAL termq3m2
  ! .......................................................................

  ! q2min : borne inferieure de q2
  ! q2max : borne superieure de q2

  ! .......................................................................
  REAL q2min
  REAL q2max
  ! .......................................................................
  ! knmin : borne inferieure de kn
  ! kmmin : borne inferieure de km
  ! .......................................................................
  REAL knmin
  REAL kmmin
  ! .......................................................................
  INTEGER ilay, ilev, igrid
  REAL tmp1, tmp2
  ! .......................................................................
  PARAMETER (kappa=0.4E+0)
  PARAMETER (long00=160.E+0)
  ! PARAMETER (gnmin=-10.E+0)
  PARAMETER (gnmin=-0.28)
  PARAMETER (gnmax=0.0233E+0)
  PARAMETER (a1=0.92E+0)
  PARAMETER (a2=0.74E+0)
  PARAMETER (b1=16.6E+0)
  PARAMETER (b2=10.1E+0)
  PARAMETER (c1=0.08E+0)
  PARAMETER (knmin=1.E-5)
  PARAMETER (kmmin=1.E-5)
  PARAMETER (q2min=1.E-5)
  PARAMETER (q2max=1.E+2)
  ! ym      PARAMETER (nlay=klev)
  ! ym      PARAMETER (nlev=klev+1)

  PARAMETER (cn1=a2*(1.E+0-6.E+0*a1/b1))
  PARAMETER (cn2=-3.E+0*a2*(6.E+0*a1+b2))
  PARAMETER (cm1=a1*(1.E+0-3.E+0*c1-6.E+0*a1/b1))
  PARAMETER (cm2=a1*(-3.E+0*a2*((b2-3.E+0*a2)*(1.E+0-6.E+0*a1/b1)- &
    3.E+0*c1*(b2+6.E+0*a1))))
  PARAMETER (cm3=-3.E+0*a2*(6.E+0*a1+b2))
  PARAMETER (cm4=-9.E+0*a1*a2)

  LOGICAL first
  SAVE first
  DATA first/.TRUE./
  !$OMP THREADPRIVATE(first)
  ! .......................................................................
  ! traitment des valeur de q2 en entree
  ! .......................................................................

  ! Initialisation de q2
  nlay = klev
  nlev = klev + 1

  CALL yamada(ngrid, dt, g, rconst, plev, temp, zlev, zlay, u, v, teta, cd, &
    q2diag, km, kn, ustar, l_mix)
  IF (first .AND. 1==1) THEN
    first = .FALSE.
    q2 = q2diag
  END IF

  DO ilev = 1, nlev
    DO igrid = 1, ngrid
      q2(igrid, ilev) = amax1(q2(igrid,ilev), q2min)
      q(igrid, ilev) = sqrt(q2(igrid,ilev))
    END DO
  END DO

  DO igrid = 1, ngrid
    tmp1 = cd(igrid)*(u(igrid,1)**2+v(igrid,1)**2)
    q2(igrid, 1) = b1**(2.E+0/3.E+0)*tmp1
    q2(igrid, 1) = amax1(q2(igrid,1), q2min)
    q(igrid, 1) = sqrt(q2(igrid,1))
  END DO

  ! .......................................................................
  ! les increments verticaux
  ! .......................................................................

  ! !!!!! allerte !!!!!c
  ! !!!!! zlev n'est pas declare a nlev !!!!!c
  ! !!!!! ---->
  DO igrid = 1, ngrid
    zlev(igrid, nlev) = zlay(igrid, nlay) + (zlay(igrid,nlay)-zlev(igrid,nlev &
      -1))
  END DO
  ! !!!!! <----
  ! !!!!! allerte !!!!!c

  DO ilay = 1, nlay
    DO igrid = 1, ngrid
      unsdz(igrid, ilay) = 1.E+0/(zlev(igrid,ilay+1)-zlev(igrid,ilay))
    END DO
  END DO
  DO igrid = 1, ngrid
    unsdzdec(igrid, 1) = 1.E+0/(zlay(igrid,1)-zlev(igrid,1))
  END DO
  DO ilay = 2, nlay
    DO igrid = 1, ngrid
      unsdzdec(igrid, ilay) = 1.E+0/(zlay(igrid,ilay)-zlay(igrid,ilay-1))
    END DO
  END DO
  DO igrid = 1, ngrid
    unsdzdec(igrid, nlay+1) = 1.E+0/(zlev(igrid,nlay+1)-zlay(igrid,nlay))
  END DO

  ! .......................................................................
  ! le cisaillement et le gradient de temperature
  ! .......................................................................

  DO igrid = 1, ngrid
    m2(igrid, 1) = (unsdzdec(igrid,1)*u(igrid,1))**2 + &
      (unsdzdec(igrid,1)*v(igrid,1))**2
    m(igrid, 1) = sqrt(m2(igrid,1))
    mpre(igrid, 1) = m(igrid, 1)
  END DO

  ! -----------------------------------------------------------------------
  DO ilev = 2, nlev - 1
    DO igrid = 1, ngrid
      ! -----------------------------------------------------------------------

      n2(igrid, ilev) = g*unsdzdec(igrid, ilev)*(teta(igrid,ilev)-teta(igrid, &
        ilev-1))/(teta(igrid,ilev)+teta(igrid,ilev-1))*2.E+0
      ! n2(igrid,ilev)=0.

      ! --->
      ! on ne sais traiter que les cas stratifies. et l'ajustement
      ! convectif est cense faire en sorte que seul des configurations
      ! stratifiees soient rencontrees en entree de cette routine.
      ! mais, bon ... on sait jamais (meme on sait que n2 prends
      ! quelques valeurs negatives ... parfois) alors :
      ! <---

      IF (n2(igrid,ilev)<0.E+0) THEN
        n2(igrid, ilev) = 0.E+0
      END IF

      m2(igrid, ilev) = (unsdzdec(igrid,ilev)*(u(igrid,ilev)-u(igrid, &
        ilev-1)))**2 + (unsdzdec(igrid,ilev)*(v(igrid,ilev)-v(igrid, &
        ilev-1)))**2
      m(igrid, ilev) = sqrt(m2(igrid,ilev))
      mpre(igrid, ilev) = m(igrid, ilev)

      ! -----------------------------------------------------------------------
    END DO
  END DO
  ! -----------------------------------------------------------------------

  DO igrid = 1, ngrid
    m2(igrid, nlev) = m2(igrid, nlev-1)
    m(igrid, nlev) = m(igrid, nlev-1)
    mpre(igrid, nlev) = m(igrid, nlev)
  END DO

  ! .......................................................................
  ! calcul des fonctions de stabilite
  ! .......................................................................

  IF (l_mix==4) THEN
    DO igrid = 1, ngrid
      sqz(igrid) = 1.E-10
      sq(igrid) = 1.E-10
    END DO
    DO ilev = 2, nlev - 1
      DO igrid = 1, ngrid
        zq = sqrt(q2(igrid,ilev))
        sqz(igrid) = sqz(igrid) + zq*zlev(igrid, ilev)*(zlay(igrid,ilev)-zlay &
          (igrid,ilev-1))
        sq(igrid) = sq(igrid) + zq*(zlay(igrid,ilev)-zlay(igrid,ilev-1))
      END DO
    END DO
    DO igrid = 1, ngrid
      long0(igrid) = 0.2*sqz(igrid)/sq(igrid)
    END DO
  ELSE IF (l_mix==3) THEN
    long0(igrid) = long00
  END IF

  ! (abd 5 2)      print*,'LONG0=',long0

  ! -----------------------------------------------------------------------
  DO ilev = 2, nlev - 1
    DO igrid = 1, ngrid
      ! -----------------------------------------------------------------------

      tmp1 = kappa*(zlev(igrid,ilev)-zlev(igrid,1))
      IF (l_mix>=10) THEN
        long(igrid, ilev) = l_mix
      ELSE
        long(igrid, ilev) = tmp1/(1.E+0+tmp1/long0(igrid))
      END IF
      long(igrid, ilev) = max(min(long(igrid,ilev),0.5*sqrt(q2(igrid,ilev))/ &
        sqrt(max(n2(igrid,ilev),1.E-10))), 5.)

      gn = -long(igrid, ilev)**2/q2(igrid, ilev)*n2(igrid, ilev)
      gm = long(igrid, ilev)**2/q2(igrid, ilev)*m2(igrid, ilev)

      gninf = .FALSE.
      gnsup = .FALSE.
      long(igrid, ilev) = long(igrid, ilev)
      long(igrid, ilev) = long(igrid, ilev)

      IF (gn<gnmin) THEN
        gninf = .TRUE.
        gn = gnmin
      END IF

      IF (gn>gnmax) THEN
        gnsup = .TRUE.
        gn = gnmax
      END IF

      sn(igrid, ilev) = cn1/(1.E+0+cn2*gn)
      sm(igrid, ilev) = (cm1+cm2*gn)/((1.E+0+cm3*gn)*(1.E+0+cm4*gn))

      IF ((gninf) .OR. (gnsup)) THEN
        snq2(igrid, ilev) = 0.E+0
        smq2(igrid, ilev) = 0.E+0
      ELSE
        snq2(igrid, ilev) = -gn*(-cn1*cn2/(1.E+0+cn2*gn)**2)
        smq2(igrid, ilev) = -gn*(cm2*(1.E+0+cm3*gn)*(1.E+0+cm4*gn)-(cm3*( &
          1.E+0+cm4*gn)+cm4*(1.E+0+cm3*gn))*(cm1+cm2*gn))/((1.E+0+cm3*gn)*( &
          1.E+0+cm4*gn))**2
      END IF

      ! abd
      ! if(ilev.le.57.and.ilev.ge.37) then
      ! print*,'L=',ilev,'   GN=',gn,'  SM=',sm(igrid,ilev)
      ! endif
      ! --->
      ! la decomposition de Taylor en q2 n'a de sens que
      ! dans les cas stratifies ou sn et sm sont quasi
      ! proportionnels a q2. ailleurs on laisse le meme
      ! algorithme car l'ajustement convectif fait le travail.
      ! mais c'est delirant quand sn et snq2 n'ont pas le meme
      ! signe : dans ces cas, on ne fait pas la decomposition.
      ! <---

      IF (snq2(igrid,ilev)*sn(igrid,ilev)<=0.E+0) snq2(igrid, ilev) = 0.E+0
      IF (smq2(igrid,ilev)*sm(igrid,ilev)<=0.E+0) smq2(igrid, ilev) = 0.E+0

      ! Correction pour les couches stables.
      ! Schema repris de JHoltzlag Boville, lui meme venant de...

      IF (1==1) THEN
        snstable = 1. - zlev(igrid, ilev)/(700.*max(ustar(igrid),0.0001))
        snstable = 1. - zlev(igrid, ilev)/400.
        snstable = max(snstable, 0.)
        snstable = snstable*snstable

        ! abde       print*,'SN ',ilev,sn(1,ilev),snstable
        IF (sn(igrid,ilev)<snstable) THEN
          sn(igrid, ilev) = snstable
          snq2(igrid, ilev) = 0.
        END IF

        IF (sm(igrid,ilev)<snstable) THEN
          sm(igrid, ilev) = snstable
          smq2(igrid, ilev) = 0.
        END IF

      END IF

      ! sn : coefficient de stabilite pour n
      ! snq2 : premier terme du developement limite de sn en q2
      ! -----------------------------------------------------------------------
    END DO
  END DO
  ! -----------------------------------------------------------------------

  ! .......................................................................
  ! calcul de km et kn au debut du pas de temps
  ! .......................................................................

  DO igrid = 1, ngrid
    kn(igrid, 1) = knmin
    km(igrid, 1) = kmmin
    kmpre(igrid, 1) = km(igrid, 1)
  END DO

  ! -----------------------------------------------------------------------
  DO ilev = 2, nlev - 1
    DO igrid = 1, ngrid
      ! -----------------------------------------------------------------------

      kn(igrid, ilev) = long(igrid, ilev)*q(igrid, ilev)*sn(igrid, ilev)
      km(igrid, ilev) = long(igrid, ilev)*q(igrid, ilev)*sm(igrid, ilev)
      kmpre(igrid, ilev) = km(igrid, ilev)

      ! -----------------------------------------------------------------------
    END DO
  END DO
  ! -----------------------------------------------------------------------

  DO igrid = 1, ngrid
    kn(igrid, nlev) = kn(igrid, nlev-1)
    km(igrid, nlev) = km(igrid, nlev-1)
    kmpre(igrid, nlev) = km(igrid, nlev)
  END DO

  ! .......................................................................
  ! boucle sur les niveaux 2 a nlev-1
  ! .......................................................................

  ! ---->
  DO ilev = 2, nlev - 1
    ! ---->
    DO igrid = 1, ngrid

      ! .......................................................................

      ! calcul des termes sources et puits de l'equation de q2
      ! ------------------------------------------------------

      knq3 = kn(igrid, ilev)*snq2(igrid, ilev)/sn(igrid, ilev)
      kmq3 = km(igrid, ilev)*smq2(igrid, ilev)/sm(igrid, ilev)

      termq = 0.E+0
      termq3 = 0.E+0
      termqm2 = 0.E+0
      termq3m2 = 0.E+0

      tmp1 = dt*2.E+0*km(igrid, ilev)*m2(igrid, ilev)
      tmp2 = dt*2.E+0*kmq3*m2(igrid, ilev)
      termqm2 = termqm2 + dt*2.E+0*km(igrid, ilev)*m2(igrid, ilev) - &
        dt*2.E+0*kmq3*m2(igrid, ilev)
      termq3m2 = termq3m2 + dt*2.E+0*kmq3*m2(igrid, ilev)

      termq = termq - dt*2.E+0*kn(igrid, ilev)*n2(igrid, ilev) + &
        dt*2.E+0*knq3*n2(igrid, ilev)
      termq3 = termq3 - dt*2.E+0*knq3*n2(igrid, ilev)

      termq3 = termq3 - dt*2.E+0*q(igrid, ilev)**3/(b1*long(igrid,ilev))

      ! .......................................................................

      ! resolution stationnaire couplee avec le gradient de vitesse local
      ! -----------------------------------------------------------------

      ! -----{on cherche le cisaillement qui annule l'equation de q^2
      ! supposee en q3}

      tmp1 = termq + termq3
      tmp2 = termqm2 + termq3m2
      m2cstat = m2(igrid, ilev) - (tmp1+tmp2)/(dt*2.E+0*km(igrid,ilev))
      mcstat = sqrt(m2cstat)

      ! abde      print*,'M2 L=',ilev,mpre(igrid,ilev),mcstat

      ! -----{puis on ecrit la valeur de q qui annule l'equation de m
      ! supposee en q3}

      IF (ilev==2) THEN
        kmcstat = 1.E+0/mcstat*(unsdz(igrid,ilev)*kmpre(igrid,ilev+1)*mpre( &
          igrid,ilev+1)+unsdz(igrid,ilev-1)*cd(igrid)*(sqrt(u(igrid,3)**2+ &
          v(igrid,3)**2)-mcstat/unsdzdec(igrid,ilev)-mpre(igrid, &
          ilev+1)/unsdzdec(igrid,ilev+1))**2)/(unsdz(igrid,ilev)+unsdz(igrid, &
          ilev-1))
      ELSE
        kmcstat = 1.E+0/mcstat*(unsdz(igrid,ilev)*kmpre(igrid,ilev+1)*mpre( &
          igrid,ilev+1)+unsdz(igrid,ilev-1)*kmpre(igrid,ilev-1)*mpre(igrid, &
          ilev-1))/(unsdz(igrid,ilev)+unsdz(igrid,ilev-1))
      END IF
      tmp2 = kmcstat/(sm(igrid,ilev)/q2(igrid,ilev))/long(igrid, ilev)
      qcstat = tmp2**(1.E+0/3.E+0)
      q2cstat = qcstat**2

      ! .......................................................................

      ! choix de la solution finale
      ! ---------------------------

      q(igrid, ilev) = qcstat
      q2(igrid, ilev) = q2cstat
      m(igrid, ilev) = mcstat
      ! abd       if(ilev.le.57.and.ilev.ge.37) then
      ! print*,'L=',ilev,'   M2=',m2(igrid,ilev),m2cstat,
      ! s     'N2=',n2(igrid,ilev)
      ! abd       endif
      m2(igrid, ilev) = m2cstat

      ! --->
      ! pour des raisons simples q2 est minore
      ! <---

      IF (q2(igrid,ilev)<q2min) THEN
        q2(igrid, ilev) = q2min
        q(igrid, ilev) = sqrt(q2min)
      END IF

      ! .......................................................................

      ! calcul final de kn et km
      ! ------------------------

      gn = -long(igrid, ilev)**2/q2(igrid, ilev)*n2(igrid, ilev)
      IF (gn<gnmin) gn = gnmin
      IF (gn>gnmax) gn = gnmax
      sn(igrid, ilev) = cn1/(1.E+0+cn2*gn)
      sm(igrid, ilev) = (cm1+cm2*gn)/((1.E+0+cm3*gn)*(1.E+0+cm4*gn))
      kn(igrid, ilev) = long(igrid, ilev)*q(igrid, ilev)*sn(igrid, ilev)
      km(igrid, ilev) = long(igrid, ilev)*q(igrid, ilev)*sm(igrid, ilev)
      ! abd
      ! if(ilev.le.57.and.ilev.ge.37) then
      ! print*,'L=',ilev,'   GN=',gn,'  SM=',sm(igrid,ilev)
      ! endif

      ! .......................................................................

    END DO

  END DO

  ! .......................................................................


  DO igrid = 1, ngrid
    kn(igrid, 1) = knmin
    km(igrid, 1) = kmmin
    ! kn(igrid,1)=cd(igrid)
    ! km(igrid,1)=cd(igrid)
    q2(igrid, nlev) = q2(igrid, nlev-1)
    q(igrid, nlev) = q(igrid, nlev-1)
    kn(igrid, nlev) = kn(igrid, nlev-1)
    km(igrid, nlev) = km(igrid, nlev-1)
  END DO

  ! CALCUL DE LA DIFFUSION VERTICALE DE Q2
  IF (1==1) THEN

    DO ilev = 2, klev - 1
      sss = sss + plev(1, ilev-1) - plev(1, ilev+1)
      sssq = sssq + (plev(1,ilev-1)-plev(1,ilev+1))*q2(1, ilev)
    END DO
    ! print*,'Q2moy avant',sssq/sss
    ! print*,'Q2q20 ',(q2(1,ilev),ilev=1,10)
    ! print*,'Q2km0 ',(km(1,ilev),ilev=1,10)
    ! ! C'est quoi ca qu'etait dans l'original???
    ! do igrid=1,ngrid
    ! q2(igrid,1)=10.
    ! enddo
    ! q2s=q2
    ! do iii=1,10
    ! call vdif_q2(dt,g,rconst,plev,temp,km,q2)
    ! do ilev=1,klev+1
    ! write(iii+49,*) q2(1,ilev),zlev(1,ilev)
    ! enddo
    ! enddo
    ! stop
    ! do ilev=1,klev
    ! print*,zlev(1,ilev),q2s(1,ilev),q2(1,ilev)
    ! enddo
    ! q2s=q2-q2s
    ! do ilev=1,klev
    ! print*,q2s(1,ilev),zlev(1,ilev)
    ! enddo
    DO ilev = 2, klev - 1
      sss = sss + plev(1, ilev-1) - plev(1, ilev+1)
      sssq = sssq + (plev(1,ilev-1)-plev(1,ilev+1))*q2(1, ilev)
    END DO
    PRINT *, 'Q2moy apres', sssq/sss


    DO ilev = 1, nlev
      DO igrid = 1, ngrid
        q2(igrid, ilev) = max(q2(igrid,ilev), q2min)
        q(igrid, ilev) = sqrt(q2(igrid,ilev))

        ! .......................................................................

        ! calcul final de kn et km
        ! ------------------------

        gn = -long(igrid, ilev)**2/q2(igrid, ilev)*n2(igrid, ilev)
        IF (gn<gnmin) gn = gnmin
        IF (gn>gnmax) gn = gnmax
        sn(igrid, ilev) = cn1/(1.E+0+cn2*gn)
        sm(igrid, ilev) = (cm1+cm2*gn)/((1.E+0+cm3*gn)*(1.E+0+cm4*gn))
        ! Correction pour les couches stables.
        ! Schema repris de JHoltzlag Boville, lui meme venant de...

        IF (1==1) THEN
          snstable = 1. - zlev(igrid, ilev)/(700.*max(ustar(igrid),0.0001))
          snstable = 1. - zlev(igrid, ilev)/400.
          snstable = max(snstable, 0.)
          snstable = snstable*snstable

          ! abde      print*,'SN ',ilev,sn(1,ilev),snstable
          IF (sn(igrid,ilev)<snstable) THEN
            sn(igrid, ilev) = snstable
            snq2(igrid, ilev) = 0.
          END IF

          IF (sm(igrid,ilev)<snstable) THEN
            sm(igrid, ilev) = snstable
            smq2(igrid, ilev) = 0.
          END IF

        END IF

        ! sn : coefficient de stabilite pour n
        kn(igrid, ilev) = long(igrid, ilev)*q(igrid, ilev)*sn(igrid, ilev)
        km(igrid, ilev) = long(igrid, ilev)*q(igrid, ilev)

      END DO
    END DO
    ! print*,'Q2km1 ',(km(1,ilev),ilev=1,10)

  END IF

  RETURN
END SUBROUTINE vdif_kcay
