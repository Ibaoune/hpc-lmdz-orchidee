!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE yamada4(ni, nsrf, ngrid, dt, g, rconst, plev, temp, zlev, zlay, u, v, teta, &
    cd, tke, km, kn, kq, ustar, iflag_pbl, drgpro)

  USE dimphy
  USE ioipsl_getin_p_mod, ONLY : getin_p

  IMPLICIT NONE
  include "iniprint.h"
  ! .......................................................................
  ! ym#include "dimensions.h"
  ! ym#include "dimphy.h"
  ! ************************************************************************************************
  !
  ! yamada4: subroutine qui calcule le transfert turbulent avec une fermeture d'ordre 2 ou 2.5
  ! 
  ! Reference: Simulation of nocturnal drainage flows by a q2l Turbulence Closure Model
  !            Yamada T.
  !            J. Atmos. Sci, 40, 91-106, 1983
  !
  !************************************************************************************************
  ! Input :
  !'======
  ! ni: indice horizontal sur la grille de base, non restreinte
  ! nsrf: type de surface
  ! ngrid: nombre de mailles concern??es sur l'horizontal
  ! dt : pas de temps
  ! g  : g
  ! rconst: constante de l'air sec
  ! zlev : altitude a chaque niveau (interface inferieure de la couche
  ! de meme indice)
  ! zlay : altitude au centre de chaque couche
  ! u,v : vitesse au centre de chaque couche
  ! (en entree : la valeur au debut du pas de temps)
  ! teta : temperature potentielle virtuelle au centre de chaque couche
  ! (en entree : la valeur au debut du pas de temps)
  ! cd : cdrag pour la quantit?? de mouvement
  ! (en entree : la valeur au debut du pas de temps)
  ! ustar: vitesse de friction calcul??e par une formule diagnostique 
  ! iflag_pbl: flag pour choisir des options du sch??ma de turbulence
  !
  !             iflag_pbl doit valoir entre 6 et 9
  !             l=6, on prend  systematiquement une longueur d'equilibre
  !             iflag_pbl=6 : MY 2.0
  !             iflag_pbl=7 : MY 2.0.Fournier
  !             iflag_pbl=8/9 : MY 2.5
  !             iflag_pbl=8 with special obsolete treatments for convergence
  !             with Cmpi5 NPv3.1 simulations
  !             iflag_pbl=10/11 :  New scheme M2 and N2 explicit and dissiptation exact
  !             iflag_pbl=12 = 11 with vertical diffusion off q2
  !
  !             2013/04/01 (FH hourdin@lmd.jussieu.fr)
  !             Correction for very stable PBLs (iflag_pbl=10 and 11)
  !             iflag_pbl=8 converges numerically with NPv3.1
  !             iflag_pbl=11 -> the model starts with NP from start files created by ce0l
  !                          -> the model can run with longer time-steps.
  !             2016/11/30 (EV etienne.vignon@univ-grenoble-alpes.fr)
  !               On met tke (=q2/2) en entr??e plut??t que q2
  !               On corrige l'update de la tke
  !
  ! Inpout/Output :
  !==============
  ! tke : tke au bas de chaque couche
  ! (en entree : la valeur au debut du pas de temps)
  ! (en sortie : la valeur a la fin du pas de temps)
 
  ! Outputs:
  !==========
  ! km : diffusivite turbulente de quantite de mouvement (au bas de chaque
  ! couche)
  ! (en sortie : la valeur a la fin du pas de temps)
  ! kn : diffusivite turbulente des scalaires (au bas de chaque couche)
  ! (en sortie : la valeur a la fin du pas de temps)
  !
  !.......................................................................

  !=======================================================================
  ! Declarations:
  !=======================================================================


  ! Inputs/Outputs
  !----------------
  REAL dt, g, rconst
  REAL plev(klon, klev+1), temp(klon, klev)
  REAL ustar(klon)
  REAL kmin, qmin, pblhmin(klon), coriol(klon)
  REAL zlev(klon, klev+1)
  REAL zlay(klon, klev)
  REAL u(klon, klev)
  REAL v(klon, klev)
  REAL teta(klon, klev)
  REAL cd(klon)
  REAL tke(klon, klev+1)
  REAL unsdz(klon, klev)
  REAL unsdzdec(klon, klev+1)
  REAL kn(klon, klev+1)
  REAL km(klon, klev+1)
  INTEGER iflag_pbl, ngrid, nsrf
  INTEGER ni(klon)

!FC
  REAL drgpro(klon,klev)
  REAL winds(klon,klev)

  ! Local
  !-------

  INCLUDE "clesphys.h"

  REAL q2(klon, klev+1)
  REAL kmpre(klon, klev+1), tmp2, qpre
  REAL mpre(klon, klev+1)
  REAL kq(klon, klev+1)
  REAL ff(klon, klev+1), delta(klon, klev+1)
  REAL aa(klon, klev+1), aa0, aa1
  INTEGER nlay, nlev

  LOGICAL,SAVE :: hboville=.TRUE.
  REAL,SAVE :: viscom,viscoh
  !$OMP THREADPRIVATE( hboville,viscom,viscoh)
  INTEGER ig, k
  REAL ri, zrif, zalpha, zsm, zsn
  REAL rif(klon, klev+1), sm(klon, klev+1), alpha(klon, klev)
  REAL m2(klon, klev+1), dz(klon, klev+1), zq, n2(klon, klev+1)
  REAL dtetadz(klon, klev+1)
  REAL m2cstat, mcstat, kmcstat
  REAL l(klon, klev+1)
  REAL zz(klon, klev+1)
  INTEGER iter
  REAL dissip(klon,klev), tkeprov,tkeexp, shear(klon,klev), buoy(klon,klev)
  REAL :: disseff

  REAL,SAVE :: ric0,ric,rifc, b1, kap
  !$OMP THREADPRIVATE(ric0,ric,rifc,b1,kap)
  DATA b1, kap/16.6, 0.4/
  REAL,SAVE :: seuilsm, seuilalpha
  !$OMP THREADPRIVATE(seuilsm, seuilalpha)
  REAL,SAVE :: lmixmin
  !$OMP THREADPRIVATE(lmixmin)
  LOGICAL, SAVE :: new_yamada4
  INTEGER, SAVE :: yamada4_num
  !$OMP THREADPRIVATE(new_yamada4,yamada4_num)
  REAL, SAVE :: yun,ydeux
  !$OMP THREADPRIVATE(yun,ydeux)

  REAL frif, falpha, fsm
  REAL rino(klon, klev+1), smyam(klon, klev), styam(klon, klev), &
    lyam(klon, klev), knyam(klon, klev), w2yam(klon, klev), t2yam(klon, klev)
  LOGICAL, SAVE :: firstcall = .TRUE.
  !$OMP THREADPRIVATE(firstcall)

  CHARACTER (len = 20) :: modname = 'yamada4'
  CHARACTER (len = 80) :: abort_message



  ! Fonctions utilis??es
  !--------------------

  frif(ri) = 0.6588*(ri+0.1776-sqrt(ri*ri-0.3221*ri+0.03156))
  falpha(ri) = 1.318*(0.2231-ri)/(0.2341-ri)
  fsm(ri) = 1.96*(0.1912-ri)*(0.2341-ri)/((1.-ri)*(0.2231-ri))
  

  IF (firstcall) THEN
! Seuil dans le code de turbulence 
    new_yamada4=.false.
    CALL getin_p('new_yamada4',new_yamada4)

    IF (new_yamada4) THEN
! Corrections et reglages issus du travail de these d'Etienne Vignon.
       ric=0.143 ! qui donne des valeurs proches des seuils proposes
                 ! dans YAMADA 1983 : sm=0.0845 (0.085 dans Y83)
                 !                    sm=1.1213 (1.12  dans Y83)
       CALL getin_p('yamada4_ric',ric)
       ric0=0.19489      ! ric=0.195 originalement, mais produisait sm<0
       ric=min(ric,ric0) ! Au dela de ric0, sm devient n??gatif
       rifc=frif(ric)
       seuilsm=fsm(frif(ric))
       seuilalpha=falpha(frif(ric))
       yun=1.
       ydeux=2.
       hboville=.FALSE.
       viscom=1.46E-5
       viscoh=2.06E-5
       lmixmin=0.
       yamada4_num=5
    ELSE
       ric=0.195
       rifc=0.191
       seuilalpha=1.12
       seuilsm=0.085
       yun=2.
       ydeux=1.
       hboville=.TRUE.
       viscom=0.
       viscoh=0.
       lmixmin=1.
       yamada4_num=0
    ENDIF

    WRITE(lunout,*)'YAMADA4 RIc, RIfc, Sm_min, Alpha_min',ric,rifc,seuilsm,seuilalpha
    firstcall = .FALSE.
    CALL getin_p('lmixmin',lmixmin)
    CALL getin_p('yamada4_hboville',hboville)
    CALL getin_p('yamada4_num',yamada4_num)
  END IF



!===============================================================================
! Flags, tests et ??valuations de constantes
!===============================================================================

! On utilise ou non la routine de Holstalg Boville pour les cas tres stables


  IF (.NOT. (iflag_pbl>=6 .AND. iflag_pbl<=12)) THEN
    abort_message='probleme de coherence dans appel a MY'
    CALL abort_physic(modname,abort_message,1)
  END IF


  nlay = klev
  nlev = klev + 1


!========================================================================
! Calcul des increments verticaux
!=========================================================================

  
! Attention: zlev n'est pas declare a nlev
  DO ig = 1, ngrid
    zlev(ig, nlev) = zlay(ig, nlay) + (zlay(ig,nlay)-zlev(ig,nlev-1))
  END DO


  DO k = 1, nlay
    DO ig = 1, ngrid
      unsdz(ig, k) = 1.E+0/(zlev(ig,k+1)-zlev(ig,k))
    END DO
  END DO
  DO ig = 1, ngrid
    unsdzdec(ig, 1) = 1.E+0/(zlay(ig,1)-zlev(ig,1))
  END DO
  DO k = 2, nlay
    DO ig = 1, ngrid
      unsdzdec(ig, k) = 1.E+0/(zlay(ig,k)-zlay(ig,k-1))
    END DO
  END DO
  DO ig = 1, ngrid
    unsdzdec(ig, nlay+1) = 1.E+0/(zlev(ig,nlay+1)-zlay(ig,nlay))
  END DO

!=========================================================================
! Richardson number and stability functions
!=========================================================================
 
! initialize arrays:

  m2(1:ngrid, :) = 0.0
  sm(1:ngrid, :) = 0.0
  rif(1:ngrid, :) = 0.0

!------------------------------------------------------------
  DO k = 2, klev

    DO ig = 1, ngrid
      dz(ig, k) = zlay(ig, k) - zlay(ig, k-1)
      m2(ig, k) = ((u(ig,k)-u(ig,k-1))**2+(v(ig,k)-v(ig, &
        k-1))**2)/(dz(ig,k)*dz(ig,k))
      dtetadz(ig, k) = (teta(ig,k)-teta(ig,k-1))/dz(ig, k)
      n2(ig, k) = g*2.*dtetadz(ig, k)/(teta(ig,k-1)+teta(ig,k))
      ri = n2(ig, k)/max(m2(ig,k), 1.E-10)
      IF (ri<ric) THEN
        rif(ig, k) = frif(ri)
      ELSE
        rif(ig, k) = rifc
      END IF
if (new_yamada4) then
        alpha(ig, k) = max(falpha(rif(ig,k)),seuilalpha)
        sm(ig, k) = max(fsm(rif(ig,k)),seuilsm)
else 
      IF (rif(ig,k)<0.16) THEN
        alpha(ig, k) = falpha(rif(ig,k))
        sm(ig, k) = fsm(rif(ig,k))
      ELSE
        alpha(ig, k) = seuilalpha
        sm(ig, k) = seuilsm
      END IF

end if
      zz(ig, k) = b1*m2(ig, k)*(1.-rif(ig,k))*sm(ig, k)
    END DO
  END DO





  !=======================================================================
  !     DIFFERENT TYPES  DE SCHEMA de  YAMADA
  !=======================================================================

  ! On commence par calculer q2 a partir de la tke

  IF (new_yamada4) THEN
      DO k=1,klev+1
         q2(1:ngrid,k)=tke(1:ngrid,k)*ydeux
      ENDDO
  ELSE
      DO k=1,klev+1
         q2(1:ngrid,k)=tke(1:ngrid,k)
      ENDDO
  ENDIF

! ====================================================================
! Computing the mixing length
! ====================================================================

 
  CALL mixinglength(ni,nsrf,ngrid,iflag_pbl,pbl_lmixmin_alpha,lmixmin,zlay,zlev,u,v,q2,n2, l)


  !--------------
  ! Yamada 2.0
  !--------------
  IF (iflag_pbl==6) THEN
 
    DO k = 2, klev
      q2(1:ngrid, k) = l(1:ngrid, k)**2*zz(1:ngrid, k)
    END DO


  !------------------
  ! Yamada 2.Fournier
  !------------------

  ELSE IF (iflag_pbl==7) THEN


    ! Calcul de l,  km, au pas precedent
    !....................................
    DO k = 2, klev
      DO ig = 1, ngrid
        delta(ig, k) = q2(ig, k)/(l(ig,k)**2*sm(ig,k))
        kmpre(ig, k) = l(ig, k)*sqrt(q2(ig,k))*sm(ig, k)
        mpre(ig, k) = sqrt(m2(ig,k))
      END DO
    END DO

    DO k = 2, klev - 1
      DO ig = 1, ngrid
        m2cstat = max(alpha(ig,k)*n2(ig,k)+delta(ig,k)/b1, 1.E-12)
        mcstat = sqrt(m2cstat)

     ! Puis on ecrit la valeur de q qui annule l'equation de m supposee en q3
     !.........................................................................

        IF (k==2) THEN
          kmcstat = 1.E+0/mcstat*(unsdz(ig,k)*kmpre(ig,k+1)*mpre(ig,k+1)+ &
            unsdz(ig,k-1)*cd(ig)*(sqrt(u(ig,3)**2+v(ig,3)**2)-mcstat/unsdzdec &
            (ig,k)-mpre(ig,k+1)/unsdzdec(ig,k+1))**2)/(unsdz(ig,k)+unsdz(ig,k &
            -1))
        ELSE
          kmcstat = 1.E+0/mcstat*(unsdz(ig,k)*kmpre(ig,k+1)*mpre(ig,k+1)+ &
            unsdz(ig,k-1)*kmpre(ig,k-1)*mpre(ig,k-1))/ &
            (unsdz(ig,k)+unsdz(ig,k-1))
        END IF

        tmp2 = kmcstat/(sm(ig,k)/q2(ig,k))/l(ig, k)
        q2(ig, k) = max(tmp2, 1.E-12)**(2./3.)

      END DO
    END DO


    ! ------------------------
    ! Yamada 2.5 a la Didi
    !-------------------------

  ELSE IF (iflag_pbl==8 .OR. iflag_pbl==9) THEN

    ! Calcul de l, km, au pas precedent
    !....................................
    DO k = 2, klev
      DO ig = 1, ngrid
        delta(ig, k) = q2(ig, k)/(l(ig,k)**2*sm(ig,k))
        IF (delta(ig,k)<1.E-20) THEN
          delta(ig, k) = 1.E-20
        END IF
        km(ig, k) = l(ig, k)*sqrt(q2(ig,k))*sm(ig, k)
        aa0 = (m2(ig,k)-alpha(ig,k)*n2(ig,k)-delta(ig,k)/b1)
        aa1 = (m2(ig,k)*(1.-rif(ig,k))-delta(ig,k)/b1)
        aa(ig, k) = aa1*dt/(delta(ig,k)*l(ig,k))
        qpre = sqrt(q2(ig,k))
        IF (aa(ig,k)>0.) THEN
          q2(ig, k) = (qpre+aa(ig,k)*qpre*qpre)**2
        ELSE
          q2(ig, k) = (qpre/(1.-aa(ig,k)*qpre))**2
        END IF
        ! else ! iflag_pbl=9
        ! if (aa(ig,k)*qpre.gt.0.9) then
        ! q2(ig,k)=(qpre*10.)**2
        ! else
        ! q2(ig,k)=(qpre/(1.-aa(ig,k)*qpre))**2
        ! endif
        ! endif
        q2(ig, k) = min(max(q2(ig,k),1.E-10), 1.E4)
      END DO
    END DO

  ELSE IF (iflag_pbl>=10) THEN

    IF (yamada4_num>=1) THEN
 
    DO k = 2, klev - 1
      DO ig=1,ngrid
      q2(ig, k) = min(max(q2(ig,k),1.E-10), 1.E4)
      km(ig, k) = l(ig, k)*sqrt(q2(ig,k))*sm(ig, k)
      shear(ig,k)=km(ig, k)*m2(ig, k)
      buoy(ig,k)=km(ig, k)*m2(ig, k)*(-1.*rif(ig,k))
      dissip(ig,k)=((sqrt(q2(ig,k)))**3)/(b1*l(ig,k))
     ENDDO
    ENDDO

    IF (yamada4_num==1) THEN ! Schema du MAR tel quel
       DO k = 2, klev - 1
         DO ig=1,ngrid
         tkeprov=q2(ig,k)/ydeux
         tkeprov= tkeprov*                           &
           &  (tkeprov+dt*(shear(ig,k)+max(0.,buoy(ig,k))))/ &
           &  (tkeprov+dt*((-1.)*min(0.,buoy(ig,k))+dissip(ig,k)))
         q2(ig,k)=tkeprov*ydeux
        ENDDO
       ENDDO
    ELSE IF (yamada4_num==2) THEN ! version modifiee avec integration exacte pour la dissipation
       DO k = 2, klev - 1
         DO ig=1,ngrid
         tkeprov=q2(ig,k)/ydeux
         disseff=dissip(ig,k)-min(0.,buoy(ig,k))
         tkeprov = tkeprov/(1.+dt*disseff/(2.*tkeprov))**2
         tkeprov= tkeprov+dt*(shear(ig,k)+max(0.,buoy(ig,k)))
         q2(ig,k)=tkeprov*ydeux
         ! En cas stable, on traite la flotabilite comme la
         ! dissipation, en supposant que buoy/q2^3 est constant.
         ! Puis on prend la solution exacte
        ENDDO
       ENDDO
    ELSE IF (yamada4_num==3) THEN ! version modifiee avec integration exacte pour la dissipation
       DO k = 2, klev - 1
         DO ig=1,ngrid
         tkeprov=q2(ig,k)/ydeux
         disseff=dissip(ig,k)-min(0.,buoy(ig,k))
         tkeprov=tkeprov*exp(-dt*disseff/tkeprov)
         tkeprov= tkeprov+dt*(shear(ig,k)+max(0.,buoy(ig,k)))
         q2(ig,k)=tkeprov*ydeux
         ! En cas stable, on traite la flotabilite comme la
         ! dissipation, en supposant que buoy/q2^3 est constant.
         ! Puis on prend la solution exacte
        ENDDO
       ENDDO
    ELSE IF (yamada4_num==4) THEN ! version modifiee avec integration exacte pour la dissipation
       DO k = 2, klev - 1
         DO ig=1,ngrid
         tkeprov=q2(ig,k)/ydeux
         tkeprov= tkeprov+dt*(shear(ig,k)+max(0.,buoy(ig,k)))
         tkeprov= tkeprov*                           &
           &  tkeprov/ &
           &  (tkeprov+dt*((-1.)*min(0.,buoy(ig,k))+dissip(ig,k)))
         q2(ig,k)=tkeprov*ydeux
         ! En cas stable, on traite la flotabilite comme la
         ! dissipation, en supposant que buoy/q2^3 est constant.
         ! Puis on prend la solution exacte
        ENDDO
       ENDDO
    ELSE IF (yamada4_num==5) THEN ! version modifiee avec integration exacte pour la dissipation
       DO k = 2, klev - 1
         DO ig=1,ngrid
         tkeprov=q2(ig,k)/ydeux

!             if(ifl_pbltree .eq. 0) then
!         disseff=dissip(ig,k)-min(0.,buoy(ig,k))
!         tkeexp=exp(-dt*disseff/tkeprov)
!         tkeprov= shear(ig,k)*tkeprov/disseff*(1.-tkeexp)+tkeprov*tkeexp
!           else
!FC on ajoute la dissipation due aux arbres
         disseff=dissip(ig,k)-min(0.,buoy(ig,k)) + drgpro(ig,k)*tkeprov
         tkeexp=exp(-dt*disseff/tkeprov)
! on prend en compte la tke cree par les arbres
         winds(ig,k)=sqrt(u(ig,k)**2+v(ig,k)**2)
         tkeprov= (shear(ig,k)+ &
          & drgpro(ig,k)*(winds(ig,k))**3)*tkeprov/disseff*(1.-tkeexp)+tkeprov*tkeexp
!               endif

         q2(ig,k)=tkeprov*ydeux

         ! En cas stable, on traite la flotabilite comme la
         ! dissipation, en supposant que buoy/q2^3 est constant.
         ! Puis on prend la solution exacte
        ENDDO
       ENDDO
    ELSE IF (yamada4_num==6) THEN ! version modifiee avec integration exacte pour la dissipation
       DO k = 2, klev - 1
         DO ig=1,ngrid
         tkeprov=q2(ig,k)/ydeux
         tkeprov=tkeprov+max(buoy(ig,k)+shear(ig,k),0.)*dt
         disseff=dissip(ig,k)-min(0.,buoy(ig,k)+shear(ig,k))
         tkeexp=exp(-dt*disseff/tkeprov)
         tkeprov= tkeprov*tkeexp
         q2(ig,k)=tkeprov*ydeux
         ! En cas stable, on traite la flotabilite comme la
         ! dissipation, en supposant que buoy/q2^3 est constant.
         ! Puis on prend la solution exacte
        ENDDO
       ENDDO
    ENDIF

    DO k = 2, klev - 1
      DO ig=1,ngrid
      q2(ig, k) = min(max(q2(ig,k),1.E-10), 1.E4)
      ENDDO
    ENDDO

   ELSE 

    DO k = 2, klev - 1
      km(1:ngrid, k) = l(1:ngrid, k)*sqrt(q2(1:ngrid,k))*sm(1:ngrid, k)
      q2(1:ngrid, k) = q2(1:ngrid, k) + ydeux*dt*km(1:ngrid, k)*m2(1:ngrid, k)*(1.-rif(1:ngrid,k))
!     q2(1:ngrid, k) = q2(1:ngrid, k) + dt*km(1:ngrid, k)*m2(1:ngrid, k)*(1.-rif(1:ngrid,k))
      q2(1:ngrid, k) = min(max(q2(1:ngrid,k),1.E-10), 1.E4)
       q2(1:ngrid, k) = 1./(1./sqrt(q2(1:ngrid,k))+dt/(yun*l(1:ngrid,k)*b1))
!     q2(1:ngrid, k) = 1./(1./sqrt(q2(1:ngrid,k))+dt/(2*l(1:ngrid,k)*b1))
      q2(1:ngrid, k) = q2(1:ngrid, k)*q2(1:ngrid, k)
    END DO

  ENDIF

  ELSE
     abort_message='Cas nom prevu dans yamada4'
     CALL abort_physic(modname,abort_message,1)

  END IF ! Fin du cas 8


  ! ====================================================================
  ! Calcul des coefficients de melange
  ! ====================================================================

  DO k = 2, klev
    DO ig = 1, ngrid
      zq = sqrt(q2(ig,k))
      km(ig, k) = l(ig, k)*zq*sm(ig, k)     ! For momentum
      kn(ig, k) = km(ig, k)*alpha(ig, k)    ! For scalars
      kq(ig, k) = l(ig, k)*zq*0.2           ! For TKE
    END DO
  END DO


  !====================================================================
  ! Transport diffusif vertical de la TKE par la TKE
  !====================================================================


    ! initialize near-surface and top-layer mixing coefficients
    !...........................................................

  kq(1:ngrid, 1) = kq(1:ngrid, 2)    ! constant (ie no gradient) near the surface
  kq(1:ngrid, klev+1) = 0            ! zero at the top

    ! Transport diffusif vertical de la TKE.
    !.......................................

  IF (iflag_pbl>=12) THEN
    q2(1:ngrid, 1) = q2(1:ngrid, 2)
    CALL vdif_q2(dt, g, rconst, ngrid, plev, temp, kq, q2)
  END IF


  !====================================================================
  ! Traitement particulier pour les cas tres stables, introduction d'une
  ! longueur de m??lange minimale
  !====================================================================
  !
  ! Reference: Local versus Nonlocal boundary-layer diffusion in a global climate model
  !            Holtslag A.A.M. and Boville B.A.
  !            J. Clim., 6, 1825-1842, 1993


 IF (hboville) THEN


  IF (prt_level>1) THEN
    WRITE(lunout,*) 'YAMADA4 0'
  END IF 

  DO ig = 1, ngrid
    coriol(ig) = 1.E-4
    pblhmin(ig) = 0.07*ustar(ig)/max(abs(coriol(ig)), 2.546E-5)
  END DO

  IF (1==1) THEN
    IF (iflag_pbl==8 .OR. iflag_pbl==10) THEN

      DO k = 2, klev
        DO ig = 1, ngrid
          IF (teta(ig,2)>teta(ig,1)) THEN
            qmin = ustar(ig)*(max(1.-zlev(ig,k)/pblhmin(ig),0.))**2
            kmin = kap*zlev(ig, k)*qmin
          ELSE
            kmin = -1. ! kmin n'est utilise que pour les SL stables.
          END IF
          IF (kn(ig,k)<kmin .OR. km(ig,k)<kmin) THEN

            kn(ig, k) = kmin
            km(ig, k) = kmin
            kq(ig, k) = kmin

 ! la longueur de melange est suposee etre l= kap z
 ! K=l q Sm d'ou q2=(K/l Sm)**2

            q2(ig, k) = (qmin/sm(ig,k))**2
          END IF
        END DO
      END DO

    ELSE
      DO k = 2, klev
        DO ig = 1, ngrid
          IF (teta(ig,2)>teta(ig,1)) THEN
            qmin = ustar(ig)*(max(1.-zlev(ig,k)/pblhmin(ig),0.))**2
            kmin = kap*zlev(ig, k)*qmin
          ELSE
            kmin = -1. ! kmin n'est utilise que pour les SL stables.
          END IF
          IF (kn(ig,k)<kmin .OR. km(ig,k)<kmin) THEN
            kn(ig, k) = kmin
            km(ig, k) = kmin
            kq(ig, k) = kmin
 ! la longueur de melange est suposee etre l= kap z
 ! K=l q Sm d'ou q2=(K/l Sm)**2
            sm(ig, k) = 1.
            alpha(ig, k) = 1.
            q2(ig, k) = min((qmin/sm(ig,k))**2, 10.)
            zq = sqrt(q2(ig,k))
            km(ig, k) = l(ig, k)*zq*sm(ig, k)
            kn(ig, k) = km(ig, k)*alpha(ig, k)
            kq(ig, k) = l(ig, k)*zq*0.2
          END IF
        END DO
      END DO
    END IF

  END IF

 END IF ! hboville

! Ajout d'une viscosite moleculaire
   km(1:ngrid,2:klev)=km(1:ngrid,2:klev)+viscom
   kn(1:ngrid,2:klev)=kn(1:ngrid,2:klev)+viscoh
   kq(1:ngrid,2:klev)=kq(1:ngrid,2:klev)+viscoh

  IF (prt_level>1) THEN
    WRITE(lunout,*)'YAMADA4 1'
  END IF !(prt_level>1) THEN


 !======================================================
 ! Estimations de w'2 et T'2 d'apres Abdela et McFarlane
 !======================================================
 !
 ! Reference: A New Second-Order Turbulence Closure Scheme for the Planetary Boundary Layer
 !            Abdella K and McFarlane N
 !            J. Atmos. Sci., 54, 1850-1867, 1997

  ! Diagnostique pour stokage
  !..........................

  IF (1==0) THEN
    rino = rif
    smyam(1:ngrid, 1) = 0.
    styam(1:ngrid, 1) = 0.
    lyam(1:ngrid, 1) = 0.
    knyam(1:ngrid, 1) = 0.
    w2yam(1:ngrid, 1) = 0.
    t2yam(1:ngrid, 1) = 0.

    smyam(1:ngrid, 2:klev) = sm(1:ngrid, 2:klev)
    styam(1:ngrid, 2:klev) = sm(1:ngrid, 2:klev)*alpha(1:ngrid, 2:klev)
    lyam(1:ngrid, 2:klev) = l(1:ngrid, 2:klev)
    knyam(1:ngrid, 2:klev) = kn(1:ngrid, 2:klev)


  ! Calcul de w'2 et T'2
  !.......................

    w2yam(1:ngrid, 2:klev) = q2(1:ngrid, 2:klev)*0.24 + &
      lyam(1:ngrid, 2:klev)*5.17*kn(1:ngrid, 2:klev)*n2(1:ngrid, 2:klev)/ &
      sqrt(q2(1:ngrid,2:klev))

    t2yam(1:ngrid, 2:klev) = 9.1*kn(1:ngrid, 2:klev)* &
      dtetadz(1:ngrid, 2:klev)**2/sqrt(q2(1:ngrid,2:klev))* &
      lyam(1:ngrid, 2:klev)
  END IF



!============================================================================
! Mise a jour de la tke
!============================================================================

  IF (new_yamada4) THEN
     DO k=1,klev+1
        tke(1:ngrid,k)=q2(1:ngrid,k)/ydeux
     ENDDO
  ELSE
     DO k=1,klev+1
        tke(1:ngrid,k)=q2(1:ngrid,k)
     ENDDO
  ENDIF


!============================================================================

  RETURN


END SUBROUTINE yamada4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE vdif_q2(timestep, gravity, rconst, ngrid, plev, temp, kmy, q2)

  USE dimphy
  IMPLICIT NONE
 
!    vdif_q2: subroutine qui calcule la diffusion de la TKE par la TKE
!             avec un schema implicite en temps avec 
!             inversion d'un syst??me tridiagonal
! 
!     Reference: Description of the interface with the surface and 
!                the computation of the turbulet diffusion in LMDZ
!                Technical note on LMDZ
!                Dufresne, J-L, Ghattas, J. and Grandpeix, J-Y
!
!============================================================================
! Declarations
!============================================================================

  REAL plev(klon, klev+1)
  REAL temp(klon, klev)
  REAL timestep
  REAL gravity, rconst
  REAL kstar(klon, klev+1), zz
  REAL kmy(klon, klev+1)
  REAL q2(klon, klev+1)
  REAL deltap(klon, klev+1)
  REAL denom(klon, klev+1), alpha(klon, klev+1), beta(klon, klev+1)
  INTEGER ngrid

  INTEGER i, k


!=========================================================================
! Calcul
!=========================================================================

  DO k = 1, klev
    DO i = 1, ngrid
      zz = (plev(i,k)+plev(i,k+1))*gravity/(rconst*temp(i,k))
      kstar(i, k) = 0.125*(kmy(i,k+1)+kmy(i,k))*zz*zz/ &
        (plev(i,k)-plev(i,k+1))*timestep
    END DO
  END DO

  DO k = 2, klev
    DO i = 1, ngrid
      deltap(i, k) = 0.5*(plev(i,k-1)-plev(i,k+1))
    END DO
  END DO
  DO i = 1, ngrid
    deltap(i, 1) = 0.5*(plev(i,1)-plev(i,2))
    deltap(i, klev+1) = 0.5*(plev(i,klev)-plev(i,klev+1))
    denom(i, klev+1) = deltap(i, klev+1) + kstar(i, klev)
    alpha(i, klev+1) = deltap(i, klev+1)*q2(i, klev+1)/denom(i, klev+1)
    beta(i, klev+1) = kstar(i, klev)/denom(i, klev+1)
  END DO

  DO k = klev, 2, -1
    DO i = 1, ngrid
      denom(i, k) = deltap(i, k) + (1.-beta(i,k+1))*kstar(i, k) + &
        kstar(i, k-1)
      alpha(i, k) = (q2(i,k)*deltap(i,k)+kstar(i,k)*alpha(i,k+1))/denom(i, k)
      beta(i, k) = kstar(i, k-1)/denom(i, k)
    END DO
  END DO

  ! Si on recalcule q2(1)
  !.......................
  IF (1==0) THEN
    DO i = 1, ngrid
      denom(i, 1) = deltap(i, 1) + (1-beta(i,2))*kstar(i, 1)
      q2(i, 1) = (q2(i,1)*deltap(i,1)+kstar(i,1)*alpha(i,2))/denom(i, 1)
    END DO
  END IF


  DO k = 2, klev + 1
    DO i = 1, ngrid
      q2(i, k) = alpha(i, k) + beta(i, k)*q2(i, k-1)
    END DO
  END DO

  RETURN
END SUBROUTINE vdif_q2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 SUBROUTINE vdif_q2e(timestep, gravity, rconst, ngrid, plev, temp, kmy, q2)
  
   USE dimphy
  IMPLICIT NONE

! vdif_q2e: subroutine qui calcule la diffusion de TKE par la TKE
!           avec un schema explicite en temps


!====================================================
! Declarations
!====================================================

  REAL plev(klon, klev+1)
  REAL temp(klon, klev)
  REAL timestep
  REAL gravity, rconst
  REAL kstar(klon, klev+1), zz
  REAL kmy(klon, klev+1)
  REAL q2(klon, klev+1)
  REAL deltap(klon, klev+1)
  REAL denom(klon, klev+1), alpha(klon, klev+1), beta(klon, klev+1)
  INTEGER ngrid
  INTEGER i, k


!==================================================
! Calcul
!==================================================

  DO k = 1, klev
    DO i = 1, ngrid
      zz = (plev(i,k)+plev(i,k+1))*gravity/(rconst*temp(i,k))
      kstar(i, k) = 0.125*(kmy(i,k+1)+kmy(i,k))*zz*zz/ &
        (plev(i,k)-plev(i,k+1))*timestep
    END DO
  END DO

  DO k = 2, klev
    DO i = 1, ngrid
      deltap(i, k) = 0.5*(plev(i,k-1)-plev(i,k+1))
    END DO
  END DO
  DO i = 1, ngrid
    deltap(i, 1) = 0.5*(plev(i,1)-plev(i,2))
    deltap(i, klev+1) = 0.5*(plev(i,klev)-plev(i,klev+1))
  END DO

  DO k = klev, 2, -1
    DO i = 1, ngrid
      q2(i, k) = q2(i, k) + (kstar(i,k)*(q2(i,k+1)-q2(i, &
        k))-kstar(i,k-1)*(q2(i,k)-q2(i,k-1)))/deltap(i, k)
    END DO
  END DO

  DO i = 1, ngrid
    q2(i, 1) = q2(i, 1) + (kstar(i,1)*(q2(i,2)-q2(i,1)))/deltap(i, 1)
    q2(i, klev+1) = q2(i, klev+1) + (-kstar(i,klev)*(q2(i,klev+1)-q2(i, &
      klev)))/deltap(i, klev+1)
  END DO

  RETURN
END SUBROUTINE vdif_q2e

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE mixinglength(ni, nsrf, ngrid,iflag_pbl,pbl_lmixmin_alpha,lmixmin,zlay,zlev,u,v,q2,n2, lmix)



  USE dimphy
  USE phys_state_var_mod, only: zstd, zsig, zmea
  USE phys_local_var_mod, only: l_mixmin, l_mix

 ! zstd: ecart type de la'altitud e sous-maille
 ! zmea: altitude moyenne sous maille
 ! zsig: pente moyenne de le maille

  USE geometry_mod, only: cell_area
  ! aire_cell: aire de la maille

  IMPLICIT NONE
!*************************************************************************
! Subrourine qui calcule la longueur de m??lange dans le sch??ma de turbulence
! avec la formule de Blackadar 
! Calcul d'un  minimum en fonction de l'orographie sous-maille:
! L'id??e est la suivante: plus il y a de relief, plus il y a du m??lange
! induit par les circulations meso et submeso ??chelles.
!
! References: * The vertical distribution of wind and turbulent exchange in a neutral atmosphere
!               Blackadar A.K.
!               J. Geophys. Res., 64, No 8, 1962
!
!             * An evaluation of neutral and convective planetary boundary-layer parametrisations relative 
!               to large eddy simulations
!               Ayotte K et al 
!               Boundary Layer Meteorology, 79, 131-175, 1996
!
!
!             * Local Similarity in the Stable Boundary Layer and Mixing length Approaches: consistency of concepts
!               Van de Wiel B.J.H et al
!               Boundary-Lay Meteorol, 128, 103-166, 2008
!
!
! Histoire:
!----------
! * premi??re r??daction, Etienne et Frederic, 09/06/2016
!
! ***********************************************************************

!==================================================================
! Declarations
!==================================================================

! Inputs
!-------
 INTEGER            ni(klon)           ! indice sur la grille original (non restreinte)
 INTEGER            nsrf               ! Type de surface
 INTEGER            ngrid              ! Nombre de points concern??s sur l'horizontal
 INTEGER            iflag_pbl          ! Choix du sch??ma de turbulence
 REAL            pbl_lmixmin_alpha  ! on active ou non le calcul de la longueur de melange minimum
 REAL               lmixmin            ! Minimum absolu de la longueur de m??lange
 REAL               zlay(klon, klev)   ! altitude du centre de la couche
 REAL               zlev(klon, klev+1) ! atitude de l'interface inf??rieure de la couche
 REAL               u(klon, klev)      ! vitesse du vent zonal
 REAL               v(klon, klev)      ! vitesse du vent meridional
 REAL               q2(klon, klev+1)   ! energie cin??tique turbulente 
 REAL               n2(klon, klev+1)   ! frequence de Brunt-Vaisala

!In/out
!-------

  LOGICAL, SAVE :: firstcall = .TRUE.
  !$OMP THREADPRIVATE(firstcall)

! Outputs
!---------

 REAL               lmix(klon, klev+1)    ! Longueur de melange  


! Local
!-------
  
 INTEGER  ig,jg, k
 REAL     h_oro(klon)
 REAL     hlim(klon)
 REAL, SAVE :: kap=0.4,kapb=0.4
  !$OMP THREADPRIVATE(kap,kapb)
 REAL zq
 REAL sq(klon), sqz(klon)
 REAL, ALLOCATABLE, SAVE :: l0(:)
  !$OMP THREADPRIVATE(l0)
 REAL fl, zzz, zl0, zq2, zn2
 REAL famorti, zzzz, zh_oro, zhlim
 REAL l1(klon, klev+1), l2(klon,klev+1)
 REAL winds(klon, klev)
 REAL xcell
 REAL zstdslope(klon)  
 REAL lmax
 REAL l2strat, l2neutre, extent  
 REAL l2limit(klon)
!===============================================================
! Fonctions utiles
!===============================================================

! Calcul de l suivant la formule de Blackadar 1962 adapt??e par Ayotte 1996
!..........................................................................

 fl(zzz, zl0, zq2, zn2) = max(min(l0(ig)*kap*zlev(ig, &
    k)/(kap*zlev(ig,k)+l0(ig)),0.5*sqrt(q2(ig,k))/sqrt( &
    max(n2(ig,k),1.E-10))), 1.E-5)
 
! Fonction d'amortissement de la turbulence au dessus de la montagne
! On augmente l'amortissement en diminuant la valeur de hlim (extent) dans le code
!.....................................................................

 famorti(zzzz, zh_oro, zhlim)=(-1.)*ATAN((zzzz-zh_oro)/(zhlim-zh_oro))*2./3.1416+1.    

  IF (ngrid==0) RETURN

  IF (firstcall) THEN
    ALLOCATE (l0(klon))
    firstcall = .FALSE.
  END IF


!=====================================================================
!         CALCUL de la LONGUEUR de m??lange suivant BLACKADAR: l1
!=====================================================================


  IF (iflag_pbl==8 .OR. iflag_pbl==10) THEN

    
    ! Iterative computation of l0
    ! This version is kept for iflag_pbl only for convergence
    ! with NPv3.1 Cmip5 simulations
    !...................................................................

    DO ig = 1, ngrid
      sq(ig) = 1.E-10
      sqz(ig) = 1.E-10
    END DO
    DO k = 2, klev - 1
      DO ig = 1, ngrid
        zq = sqrt(q2(ig,k))
        sqz(ig) = sqz(ig) + zq*zlev(ig, k)*(zlay(ig,k)-zlay(ig,k-1))
        sq(ig) = sq(ig) + zq*(zlay(ig,k)-zlay(ig,k-1))
      END DO
    END DO
    DO ig = 1, ngrid
      l0(ig) = 0.2*sqz(ig)/sq(ig)
    END DO
    DO k = 2, klev
      DO ig = 1, ngrid
        l1(ig, k) = fl(zlev(ig,k), l0(ig), q2(ig,k), n2(ig,k))
      END DO
    END DO

  ELSE

    
    ! In all other case, the assymptotic mixing length l0 is imposed (150m)
    !......................................................................

    l0(1:ngrid) = 150.
    DO k = 2, klev
      DO ig = 1, ngrid
        l1(ig, k) = fl(zlev(ig,k), l0(ig), q2(ig,k), n2(ig,k))
      END DO
    END DO

  END IF

!=================================================================================
!  CALCUL d'une longueur de melange en fonctions de la topographie sous maille: l2
! si plb_lmixmin_alpha=TRUE et si on se trouve sur de la terre ( pas actif sur les 
! glacier, la glace de mer et les oc??ans)
!=================================================================================

   l2(1:ngrid,:)=0.0
   l_mixmin(1:ngrid,:,nsrf)=0.
   l_mix(1:ngrid,:,nsrf)=0.

   IF (nsrf .EQ. 1) THEN

! coefficients 
!--------------

     extent=2.                                                         ! On ??tend l'impact du relief jusqu'?? extent*h, extent >1.  
     lmax=150.                                                         ! Longueur de m??lange max dans l'absolu

! calculs
!---------

     DO ig=1,ngrid

      ! On calcule la hauteur du relief
      !.................................
      ! On ne peut pas prendre zstd seulement pour caracteriser le relief sous maille
      ! car sur un terrain pentu mais sans relief, zstd est non nul (comme en Antarctique, C. Genthon)
      ! On corrige donc zstd avec l'ecart type de l'altitude dans la maille sans relief 
      ! (en gros, une maille de taille xcell avec une pente constante zstdslope)
      jg=ni(ig) 
!     IF (zsig(jg) .EQ. 0.) THEN
!          zstdslope(ig)=0.         
!     ELSE
!     xcell=sqrt(cell_area(jg))
!     zstdslope(ig)=max((xcell*zsig(jg)-zmea(jg))**3 /(3.*zsig(jg)),0.)
!     zstdslope(ig)=sqrt(zstdslope(ig))
!     END IF
      
!     h_oro(ig)=max(zstd(jg)-zstdslope(ig),0.)   ! Hauteur du relief
      h_oro(ig)=zstd(jg)
      hlim(ig)=extent*h_oro(ig)     
     ENDDO

     l2limit(1:ngrid)=0.

     DO k=2,klev
        DO ig=1,ngrid
           winds(ig,k)=sqrt(u(ig,k)**2+v(ig,k)**2)
           IF (zlev(ig,k) .LE. h_oro(ig)) THEN  ! sous l'orographie
              l2strat= kapb*pbl_lmixmin_alpha*winds(ig,k)/sqrt(max(n2(ig,k),1.E-10))  ! si stratifi??, amplitude d'oscillation * kappab (voir Van de Wiel et al 2008)
              l2neutre=kap*zlev(ig,k)*h_oro(ig)/(kap*zlev(ig,k)+h_oro(ig))            ! Dans le cas neutre, formule de blackadar. tend asymptotiquement vers h
              l2neutre=MIN(l2neutre,lmax)                                             ! On majore par lmax 
              l2limit(ig)=MIN(l2neutre,l2strat)                                       ! Calcule de l2 (minimum de la longueur en cas neutre et celle en situation stratifi??e)
              l2(ig,k)=l2limit(ig)
                                      
           ELSE IF (zlev(ig,k) .LE. hlim(ig)) THEN ! Si on est au dessus des montagnes, mais affect?? encore par elles

      ! Au dessus des montagnes, on prend la l2limit au sommet des montagnes 
      ! (la derni??re calcul??e dans la boucle k, vu que k est un indice croissant avec z)
      ! et on multiplie l2limit par une fonction qui d??croit entre h et hlim
              l2(ig,k)=l2limit(ig)*famorti(zlev(ig,k),h_oro(ig), hlim(ig))
           ELSE                                                                    ! Au dessus de extent*h, on prend l2=l0 
              l2(ig,k)=0.
           END IF
        ENDDO
     ENDDO
   ENDIF                                                                        ! pbl_lmixmin_alpha

!==================================================================================
! On prend le max entre la longueur de melange de blackadar et celle calcul??e
! en fonction de la topographie 
!===================================================================================


 DO k=2,klev
    DO ig=1,ngrid
       lmix(ig,k)=MAX(MAX(l1(ig,k), l2(ig,k)),lmixmin)
   ENDDO
 ENDDO

! Diagnostics

 DO k=2,klev
    DO ig=1,ngrid
       jg=ni(ig)
       l_mix(jg,k,nsrf)=lmix(ig,k)
       l_mixmin(jg,k,nsrf)=l2(ig,k)
    ENDDO
 ENDDO
 DO ig=1,ngrid
    jg=ni(ig)
    l_mix(jg,1,nsrf)=hlim(ig)
 ENDDO



END SUBROUTINE mixinglength
