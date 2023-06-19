!
! $Id $
!
SUBROUTINE cvltr_noscav(it,pdtime,da, phi, mp,wght_cvfd,paprs,pplay,x,upd,dnd,dx)
  USE dimphy
  USE infotrac_phy, ONLY : nbtr
  IMPLICIT NONE 
!=====================================================================
! Objet : convection des traceurs / KE
! Auteurs: M-A Filiberti and J-Y Grandpeix
!=====================================================================
  include "YOMCST.h"
  include "YOECUMF.h" 

! Entree
  REAL,INTENT(IN)                           :: pdtime
  INTEGER, INTENT(IN)                       :: it
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: da
  REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: phi
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: mp
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: wght_cvfd  ! weights of the layers feeding convection
  REAL,DIMENSION(klon,klev+1),INTENT(IN)    :: paprs ! pression aux 1/2 couches (bas en haut)
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: pplay ! pression pour le milieu de chaque couche
  REAL,DIMENSION(klon,klev,nbtr),INTENT(IN)      :: x     ! q de traceur (bas en haut) 
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: upd   ! saturated updraft mass flux
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: dnd   ! saturated downdraft mass flux

! Sortie
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT) :: dx ! tendance de traceur  (bas en haut)

! Variables locales     
! REAL,DIMENSION(klon,klev)       :: zed
  REAL,DIMENSION(klon,klev,klev)  :: zmd
  REAL,DIMENSION(klon,klev,klev)  :: za
  REAL,DIMENSION(klon,klev)       :: zmfd,zmfa
  REAL,DIMENSION(klon,klev)       :: zmfp,zmfu
  REAL,DIMENSION(klon,nbtr)       :: qfeed     ! tracer concentration feeding convection
  REAL,DIMENSION(klon,klev)       :: deltap
  INTEGER                         :: i,k,j 
  REAL                            :: pdtimeRG
  REAL                            :: smallest_mp
  real conserv
  real smfd
  real smfu
  real smfa
  real smfp
! =========================================
! calcul des tendances liees au downdraft
! =========================================
!
  smallest_mp = tiny(mp(1,1))
!cdir collapse
  qfeed(:,it) = 0.
  DO j=1,klev
  DO i=1,klon
!   zed(i,j)=0.
    zmfd(i,j)=0.
    zmfa(i,j)=0.
    zmfu(i,j)=0.
    zmfp(i,j)=0.
  END DO
  END DO
!cdir collapse
  DO k=1,klev
  DO j=1,klev
  DO i=1,klon
    zmd(i,j,k)=0.
    za (i,j,k)=0.
  END DO
  END DO
  END DO
! entrainement
! DO k=1,klev-1
!    DO i=1,klon
!       zed(i,k)=max(0.,mp(i,k)-mp(i,k+1))
!    END DO
! END DO

! calcul de la matrice d echange
! matrice de distribution de la masse entrainee en k

  DO k=1,klev-1
     DO i=1,klon
        zmd(i,k,k)=max(0.,mp(i,k)-mp(i,k+1))
     END DO
  END DO
  DO k=2,klev
     DO j=k-1,1,-1
        DO i=1,klon
!!           if(mp(i,j+1).ne.0) then
!!              zmd(i,j,k)=zmd(i,j+1,k)*min(1.,mp(i,j)/mp(i,j+1))
!!           ENDif
           zmd(i,j,k)=zmd(i,j+1,k)*mp(i,j)/max(mp(i,j),mp(i,j+1),smallest_mp)
        END DO
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev-1
        DO i=1,klon
           za(i,j,k)=max(0.,zmd(i,j+1,k)-zmd(i,j,k))
        END DO
     END DO
  END DO
!
! rajout du terme lie a l ascendance induite
!
  DO j=2,klev
     DO i=1,klon
        za(i,j,j-1)=za(i,j,j-1)+mp(i,j)
     END DO
  END DO
!
! tendances
!            
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfd(i,j)=zmfd(i,j)+za(i,j,k)*(x(i,k,it)-x(i,j,it))
        END DO
     END DO
  END DO
!
! =========================================
! calcul des tendances liees aux flux satures
! =========================================
!RL
!  Feeding concentrations
  DO j=1,klev
     DO i=1,klon
        qfeed(i,it)=qfeed(i,it)+wght_cvfd(i,j)*x(i,j,it)
     END DO
  END DO
!RL
!
  DO j=1,klev
     DO i=1,klon
!RL
!!        zmfa(i,j,it)=da(i,j)*(x(i,1,it)-x(i,j,it))                     ! da
        zmfa(i,j)=da(i,j)*(qfeed(i,it)-x(i,j,it))                     ! da
!RL
     END DO
  END DO
!
!!  print *,'it, qfeed(1,it), x(1,1,it) ', it, qfeed(1,it), x(1,1,it)  !jyg
!!  print *,'wght_cvfd ', (j, wght_cvfd(1,j), j=1,5)                     !jyg
!
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfp(i,j)=zmfp(i,j)+phi(i,j,k)*(x(i,k,it)-x(i,j,it))
        END DO
     END DO
  END DO
  DO j=1,klev-1
     DO i=1,klon
        zmfu(i,j)=max(0.,upd(i,j+1)+dnd(i,j+1))*(x(i,j+1,it)-x(i,j,it))
     END DO
  END DO
  DO j=2,klev
     DO i=1,klon
        zmfu(i,j)=zmfu(i,j)+min(0.,upd(i,j)+dnd(i,j))*(x(i,j,it)-x(i,j-1,it))
     END DO
  END DO

! =========================================
! calcul final des tendances
! =========================================
  DO k=1, klev
     DO i=1, klon
        deltap(i,k)=paprs(i,k)-paprs(i,k+1)
     ENDDO
  ENDDO
  pdtimeRG=pdtime*RG
!cdir collapse
  DO k=1, klev
     DO i=1, klon
        dx(i,k,it)=(zmfd(i,k)+zmfu(i,k)       &
                +zmfa(i,k)+zmfp(i,k))*pdtimeRG/deltap(i,k)
     ENDDO
  ENDDO

!! test de conservation du traceur
      conserv=0.
      smfd = 0.
      smfu = 0.
      smfa = 0.
      smfp = 0.
      DO k=1, klev
        DO i=1, klon
         conserv=conserv+dx(i,k,it)*   &
          deltap(i,k)/RG
         smfd = smfd + zmfd(i,k)*pdtime
         smfu = smfu + zmfu(i,k)*pdtime
         smfa = smfa + zmfa(i,k)*pdtime
         smfp = smfp + zmfp(i,k)*pdtime
        ENDDO
      ENDDO
!!      print *,'it',it,'cvltr_noscav conserv, smfd, smfu, smfa, smfp ',conserv,  &
!!               smfd, smfu, smfa, smfp
     
END SUBROUTINE cvltr_noscav
