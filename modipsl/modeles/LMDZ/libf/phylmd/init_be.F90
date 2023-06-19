!$Id $

SUBROUTINE init_be(pctsrf,pplay,masktr,tautr,vdeptr,scavtr,srcbe)
!!!SUBROUTINE init_be(pctsrf,masktr,tautr,vdeptr,scavtr,srcbe)

  USE dimphy
  USE infotrac_phy, ONLY : nbtr
  USE indice_sol_mod
  USE geometry_mod, ONLY : longitude, latitude
    
  IMPLICIT NONE 
!=====================================================================
! Objet : prescription d'une source de Beryllium 7 
!         pour 19 niveaux verticaux
!        (d'apres le diagramme de Lal and Peters, 1967)
!
!
! written by : O. Coindreau (CEA/LDG) 05/2005
! last modified by : A. Jamelot (LMD/CEA)  04/03/2009 
!=====================================================================

  INCLUDE "YOMCST.h"
  INCLUDE "YOECUMF.h" 

!
! Input Arguments
!
  REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: pctsrf !Pourcentage de sol (f(nature du sol))
  REAL,DIMENSION(klon,klev), INTENT(IN) :: pplay  ! Pressions en milieu de couches
!
! Output Arguments
!
  REAL,DIMENSION(klon),INTENT(OUT)      :: masktr ! Masque de l'echange avec la surface (possible => 1 )
  REAL,INTENT(OUT)                      :: tautr  ! Constante de decroissance radioactive
  REAL,INTENT(OUT)                      :: vdeptr ! Vitesse de depot sec dans la couche Brownienne
  REAL,INTENT(OUT)                      :: scavtr ! Coefficient de lessivage
  REAL,DIMENSION(klon,klev),INTENT(OUT) :: srcbe  ! source volumique de 7Be      
!
! Local Variables
!
!!!  INTEGER              :: iref      ! numero d'un point oceanique donnant la grille de pression de reference
  REAL,DIMENSION(klon) :: rlatgeo   ! latitudes geomagnetiques de la grille
  REAL                 :: glt       ! latitude du pole geomagnetique
  REAL                 :: glg       ! longitude du pole geomagnetique
  REAL                 :: latgeo,qcos
  INTEGER              :: k,i, kref, k2
  INTEGER              :: nref
  PARAMETER (nref=39)
  REAL,DIMENSION(nref), SAVE :: pref      ! grille de pression de reference (bas des couches)
  DATA pref  /   &
      101249.99999999994, 100387.17261011522, 99447.35334189111,  98357.43412194174,   &
      97046.47707771382,  95447.1116450629,   93496.85259615642,  91139.46548240296,   &
      88326.55568744117,  85019.60710580258,  81192.7404556645,   76836.48366938648,   &
      71962.81275769137,  66611.56331321516,  60857.914829743604, 54819.84484441629,   &
      48663.06257114699,  42598.95465845692,  36869.104365898806, 31709.927925633147,  &
      27296.757208636915, 23682.282929080895, 20766.025578936627, 18336.105961406534,  &
      16178.04816768436,  14168.286905562818, 12275.719926478887, 10507.798835225762,  &
      8876.585404909414,  7391.283929569539,  6057.514475749798,  4877.165909157005,   &
      3848.34936408203,   2965.444753540027,  2219.2391544640013, 1597.15366044666,    &
      1083.5531161631498, 660.1311067852655,  306.36072267002805 / 
!$OMP THREADPRIVATE(pref)

  WRITE(*,*)'PASSAGE init_be ...'

! la source est maintenant définie independemment de la valeur de klev.
!!! Source actuellement definie pour klev = 19 et klev >= 39
!!  IF (klev /= 19 .AND. klev<39) CALL abort_physic("init_be","Source du be7 necessite klev=19 ou klev>=39",1)
!!!
! Definition des constantes
! -------------------------
  tautr = 6645000.
  vdeptr = 1.E-3 
  scavtr = 0.5 
!!!!!jyg le 13/03/2013; puis 20/03/2013 : pref est maintenant une table.
!!!
!!!   Recherche d'un point rlat=0., rlon=180.
!!      iref=(klon+1)/2
!!      DO i = 1,klon
!!        IF (abs(rlatd(i)) .LT. 0.15 .AND. cos(rlond(i)) .LT. -0.85) iref=i
!!      ENDDO
!!!
!!!   Grille de pression de reference (= approx de sommets de couches)
!!      pref(1) = pplay(iref,1)+0.5*(pplay(iref,1)-pplay(iref,2))
!!      DO k = 2,klev
!!        pref(k) = 0.5*(pplay(iref,k-1)+pplay(iref,k))
!!      ENDDO
!!!

  WRITE(*,*) '-------------- SOURCE DE BERYLLIUM ------------------- '
  WRITE(*,*)'Decroissance (s): ', tautr
  WRITE(*,*)'Vitesse de depot sec: ',vdeptr
  WRITE(*,*)'Facteur de lessivage: ',scavtr

  DO i = 1,klon
     masktr(i) = 0.
     IF ( NINT(pctsrf(i,1)) .EQ. 1 ) masktr(i) = 1.
  END DO

! Premiers niveaux: source nulle
! ------------------------------
  DO k = 1,6
     DO i = 1,klon
        srcbe(i,k) = 0.
     END DO
  END DO
!
! Pour les autres niveaux:
! 1-passer des coordonnees geographiques a la latitude geomagnetique
! 2-prescrire la source de Be (en 10exp5 at/g/s) dans ce repere
! 3-mettre la source de Be ds la bonne unite (en at/kgA/s)
!
  glt =  78.5*rpi/180.
  glg = -69.0*rpi/180.

  DO i = 1,klon
     qcos=sin(glt)*sin(latitude(i))
!!jyg
!!     qcos=qcos+cos(glt)*cos(latitude(i))*cos(longitude(i)+glg)
     qcos=qcos+cos(glt)*cos(latitude(i))*cos(longitude(i)-glg)
!!jyg end
     IF ( qcos .LT. -1.) qcos = -1.
     IF ( qcos .GT. 1.)  qcos = 1.
     rlatgeo(i)=rpi/2.-acos(qcos)
  ENDDO

!!!===========================
!!!  Cas 19 niveaux verticaux
!!!===========================
!!  IF (klev.eq.19) then
!!     DO k = 1,klev
!!        DO i = 1,klon
!!!!!jyg le 13/03/2013
!!!
!!! k est le niveau dans la grille locale
!!! Determination du niveau kref dans la grille de refernce
!!      kref = 1
!!      DO k2 = 1,klev
!!        IF (pref(k2) .GT. pplay(i,k)) kref=k2
!!      ENDDO
!!!!!
!!           latgeo=(180./rpi)*abs(rlatgeo(i))
!!           IF ( kref .EQ. 1 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=0.1
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.09
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.07
!!           END IF
!!           IF ( kref .EQ. 2 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=0.12
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.1
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.09
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.07
!!           END IF
!!           IF ( kref .EQ. 3 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=0.14
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.12
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.1
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.09
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.07
!!           END IF
!!           IF ( kref .EQ. 4 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=0.175
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.16
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.14
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.12
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.1
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.09
!!           END IF
!!           IF ( kref .EQ. 5 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=0.28
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.26
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.23
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.175
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.14
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.12
!!           END IF
!!           IF ( kref .EQ. 6 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=0.56
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.49
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.42
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.28
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.26
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.245
!!           END IF
!!           IF ( kref .EQ. 7 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=1.05
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.875
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.7
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.52
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.44
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.385
!!           END IF
!!           IF ( kref .EQ. 8 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=2.
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=1.8
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=1.5
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=1.
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.8
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.75
!!           END IF
!!           IF ( kref .EQ. 9 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=4.
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=3.5
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=3.
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=2.5
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=1.8
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=1.4
!!           END IF
!!           IF ( kref .EQ. 10 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=8.5
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=8.
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=7.
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=4.5
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=3.5
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=3.
!!           END IF
!!           IF ( kref .EQ. 11 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=17.
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=15.
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=11.
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=8.
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=5.
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=4.
!!           END IF
!!           IF ( kref .EQ. 12 ) THEN
!!              IF (latgeo.GE.50.0) srcbe(i,k)=25.
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=22.
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=17.
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=11.
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=7.5
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=7.
!!           END IF
!!           IF ( kref .EQ. 13 ) THEN
!!              IF (latgeo.GE.60.0) srcbe(i,k)=33.
!!              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=32.
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=30.
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=22.
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=15.
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=11.
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=8.
!!           END IF
!!           IF ( kref .EQ. 14 ) THEN
!!              IF (latgeo.GE.60.0) srcbe(i,k)=48.
!!              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=45.
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=36.
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=26.
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=17.5
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=12.5
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=10.
!!           END IF
!!           IF ( kref .EQ. 15 ) THEN
!!              IF (latgeo.GE.70.0) srcbe(i,k)=58.
!!              IF (latgeo.GE.60.0 .and. latgeo.LT.70.0) srcbe(i,k)=57.
!!              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=50.
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=38.
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=25.
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=15.
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=12.5
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=10.
!!           END IF
!!           IF ( kref .EQ. 16 ) THEN
!!              IF (latgeo.GE.70.0) srcbe(i,k)=70.
!!              IF (latgeo.GE.60.0 .and. latgeo.LT.70.0) srcbe(i,k)=65.
!!              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=50.
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=32.
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=20.
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=13.
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=9.
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=7.5
!!           END IF
!!           IF ( kref .GE. 17 ) THEN
!!              IF (latgeo.GE.70.0) srcbe(i,k)=80.
!!              IF (latgeo.GE.60.0 .and. latgeo.LT.70.0) srcbe(i,k)=70.
!!              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=45.
!!              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=27.
!!              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=17.5
!!              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=12.
!!              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=8.
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=7.
!!           END IF
!!        END DO
!!     END DO
!!  END IF ! fin de 19 niveaux verticaux
!!
!!!!!!  IF (klev .ge. 39) then
     DO k = 1,klev
        DO i = 1,klon
!!!jyg le 13/03/2013
!
! k est le niveau dans la grille locale
! Determination du niveau kref dans la grille de refernce
      kref = 1
      DO k2 = 1,nref
        IF (pref(k2) .GT. pplay(i,k)) kref=k2
      ENDDO
!!!
           latgeo=(180./rpi)*abs(rlatgeo(i))
           IF ( kref .LE. 4 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=0.07
           END IF
           IF ( kref .EQ. 5 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=0.1
              IF (latgeo.GE.20.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.09
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.07
           END IF
           IF ( kref .EQ. 6 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=0.14
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.12
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.1
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.09
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.07
           END IF
           IF ( kref .EQ. 7 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=0.16
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.16
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.14
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.12
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.1
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.09
           END IF
           IF ( kref .EQ. 8 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=0.175
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.16
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.14
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.12
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.1
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.1
           END IF
           IF ( kref .EQ. 9 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=0.245
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.21
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.175
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.14
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.12
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.12
           END IF
           IF ( kref .EQ. 10 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=0.31
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.28
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.245
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.21
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.16
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.14
           END IF
           IF ( kref .EQ. 11 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=0.35
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.3
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.3
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.2
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.18
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.16
           END IF
           IF ( kref .EQ. 12 ) THEN
              IF (latgeo.GE.40.0) srcbe(i,k)=0.5
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.4
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.35
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.3
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.25
           END IF
           IF ( kref .EQ. 13 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=0.8
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=0.7
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.6
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.5
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.4
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.35
           END IF
           IF ( kref .EQ. 14 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=1.2
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=1.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=0.75
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.6
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.5
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.4
           END IF
           IF ( kref .EQ. 15 ) THEN
              IF (latgeo.GE.60.0) srcbe(i,k)=1.75
              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=1.8 
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=1.6
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=1.4
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=0.9
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=0.75
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.65
           END IF
           IF ( kref .EQ. 16 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=3.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=2.5
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=1.8
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=1.5
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=1.2
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=0.9
           END IF
           IF ( kref .EQ. 17 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=4.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=3.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=2.5
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=2.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=1.6
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=1.4
           END IF
           IF ( kref .EQ. 18 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=7.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=6.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=4.5
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=3.5
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=3.
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=2.
           END IF
           IF ( kref .EQ. 19 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=8.5
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=8.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=7.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=4.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=3.5
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=3.
           END IF
           IF ( kref .EQ. 20 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=12.5
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=12.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=8.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=6.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=4.
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=3.5
           END IF
           IF ( kref .EQ. 21 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=16.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=13.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=10.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=7.5
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=4.5
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=4.
           END IF
           IF ( kref .EQ. 22 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=20.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=17.5
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=12.5
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=9.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=6.
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=4.5
           END IF
           IF ( kref .EQ. 23 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=25.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=22.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=15.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=10.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=7.5
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=6.
           END IF
           IF ( kref .EQ. 24 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=28.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=26.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=18.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=12.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=8.5
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=7.
           END IF
           IF ( kref .EQ. 25 ) THEN
              IF (latgeo.GE.50.0) srcbe(i,k)=33.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=28.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=20.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=14.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=10.
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=8.5
           END IF
           IF ( kref .EQ. 26 ) THEN
              IF (latgeo.GE.60.0) srcbe(i,k)=38.
              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=36.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=32.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=24.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=15.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=11.5
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=10.
!!jyg
!!              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=6.
!!jyg end
           END IF
           IF ( kref .EQ. 27 ) THEN
              IF (latgeo.GE.60.0) srcbe(i,k)=46.
              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=44.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=35.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=25.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=16.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=12.5
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=10.
           END IF
           IF ( kref .EQ. 28 ) THEN
              IF (latgeo.GE.60.0) srcbe(i,k)=53.
              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=48.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=37.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=25.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=16.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=12.5
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=10.
           END IF
           IF ( kref .EQ. 29 ) THEN
              IF (latgeo.GE.70.0) srcbe(i,k)=58.
              IF (latgeo.GE.60.0 .and. latgeo.LT.70.0) srcbe(i,k)=56.
              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=50.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=36.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=24.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=15.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=11.5
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=10.
           END IF
           IF ( kref .EQ. 30 ) THEN
              IF (latgeo.GE.70.0) srcbe(i,k)=65.
              IF (latgeo.GE.60.0 .and. latgeo.LT.70.0) srcbe(i,k)=60.
              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=50.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=35.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=22.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=14.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=10.
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=9.
           END IF
           IF ( kref .EQ. 31 ) THEN
              IF (latgeo.GE.70.0) srcbe(i,k)=70.
              IF (latgeo.GE.60.0 .and. latgeo.LT.70.0) srcbe(i,k)=62.
              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=48.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=32.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=21.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=13.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=9.
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=7.6
           END IF
           IF ( kref .EQ. 32 ) THEN
              IF (latgeo.GE.70.0) srcbe(i,k)=80.
              IF (latgeo.GE.60.0 .and. latgeo.LT.70.0) srcbe(i,k)=60.
              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=46.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=30.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=17.5
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=11.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=8.
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=7.4
           END IF
           IF ( kref .GE. 33 ) THEN
              IF (latgeo.GE.70.0) srcbe(i,k)=80.
              IF (latgeo.GE.60.0 .and. latgeo.LT.70.0) srcbe(i,k)=70.
              IF (latgeo.GE.50.0 .and. latgeo.LT.60.0) srcbe(i,k)=45.
              IF (latgeo.GE.40.0 .and. latgeo.LT.50.0) srcbe(i,k)=27.
              IF (latgeo.GE.30.0 .and. latgeo.LT.40.0) srcbe(i,k)=15.
              IF (latgeo.GE.20.0 .and. latgeo.LT.30.0) srcbe(i,k)=10.
              IF (latgeo.GE.10.0 .and. latgeo.LT.20.0) srcbe(i,k)=7.6
              IF (latgeo.GE.0.0 .and. latgeo.LT.10.0) srcbe(i,k)=7.
           END IF
        END DO
     END DO
!!!!!!  END IF ! fin de 39 niveaux verticaux


!====================================
! Conversion de la source en U/s/kgA
!====================================
  DO k = 1,klev
     DO i = 1,klon
       ! La source est  at/min/m3 -> at/s/kgA
       ! avec une masse volumique de l'air = 1.295 kg/m3
       ! 1/(60*1.295) = 0.01287
       srcbe(i,k)=srcbe(i,k)*0.01287
!!       print *,'  k, srcbe(i,k) ',   &
!!                  k, srcbe(i,k)
       ! La source est  at/min/m3 -> at/s/m3
       ! srcbe(i,k)=srcbe(i,k)*0.0166667
    END DO
 END DO

END SUBROUTINE init_be
