SUBROUTINE alpale ( debut, itap, dtime, paprs, omega, t_seri,   &
                    alp_offset, it_wape_prescr,  wape_prescr, fip_prescr, &
                    ale_bl_prescr, alp_bl_prescr, &
                    wake_pe, wake_fip,  &
                    Ale_bl, Ale_bl_trig, Alp_bl, &
                    Ale, Alp, Ale_wake, Alp_wake )

! **************************************************************
! *
! ALPALE                                                       *
! *
! *
! written by   : Jean-Yves Grandpeix, 12/05/2016              *
! modified by :                                               *
! **************************************************************

  USE dimphy
  USE ioipsl_getin_p_mod, ONLY : getin_p
  USE print_control_mod, ONLY: mydebug=>debug , lunout, prt_level
  USE phys_local_var_mod, ONLY: zw2       ! Variables internes non sauvegardees de la physique
!
  IMPLICIT NONE

!================================================================
! Auteur(s)   : Jean-Yves Grandpeix, 12/05/2016
! Objet : Sums up all contributions to Ale and Alp
!================================================================

! Input arguments
!----------------
  LOGICAL, INTENT(IN)                                        :: debut
  INTEGER, INTENT(IN)                                        :: itap
  REAL, INTENT(IN)                                           :: dtime
  INTEGER, INTENT(IN)                                        :: it_wape_prescr
  REAL, INTENT(IN)                                           :: wape_prescr, fip_prescr
  REAL, INTENT(IN)                                           :: Ale_bl_prescr, Alp_bl_prescr
  REAL, INTENT(IN)                                           :: alp_offset
  REAL, DIMENSION(klon,klev+1), INTENT(IN)                   :: paprs
  REAL, DIMENSION(klon,klev), INTENT(IN)                     :: t_seri
  REAL, DIMENSION(klon,klev), INTENT(IN)                     :: omega
  REAL, DIMENSION(klon), INTENT(IN)                          :: wake_pe, wake_fip
  REAL, DIMENSION(klon), INTENT(IN)                          :: Ale_bl, Ale_bl_trig, Alp_bl


! Output arguments
!----------------
  REAL, DIMENSION(klon), INTENT(OUT)                         :: Ale, Alp
  REAL, DIMENSION(klon), INTENT(OUT)                         :: Ale_wake, Alp_wake

  include "thermcell.h"
  include "YOMCST.h"
  include "YOETHF.h"

! Local variables
!----------------
  INTEGER                                                    :: i, k
  REAL, DIMENSION(klon)                                      :: www
  REAL, SAVE                                                 :: ale_max=1000.
  REAL, SAVE                                                 :: alp_max=2.
  CHARACTER*20 modname
  CHARACTER*80 abort_message


    !$OMP THREADPRIVATE(ale_max,alp_max)

       ! Calcul de l'energie disponible ALE (J/kg) et de la puissance
       ! disponible ALP (W/m2) pour le soulevement des particules dans
       ! le modele convectif
       !
       do i = 1,klon
          ALE(i) = 0.
          ALP(i) = 0.
       enddo
       !
       !calcul de ale_wake et alp_wake
       if (iflag_wake>=1) then
          if (itap .le. it_wape_prescr) then
             do i = 1,klon
                ale_wake(i) = wape_prescr
                alp_wake(i) = fip_prescr
             enddo
          else
             do i = 1,klon
                !jyg  ALE=WAPE au lieu de ALE = 1/2 Cstar**2
                !cc           ale_wake(i) = 0.5*wake_cstar(i)**2
                ale_wake(i) = wake_pe(i)
                alp_wake(i) = wake_fip(i)
             enddo
          endif
       else
          do i = 1,klon
             ale_wake(i) = 0.
             alp_wake(i) = 0.
          enddo
       endif
       !combinaison avec ale et alp de couche limite: constantes si pas
       !de couplage, valeurs calculees dans le thermique sinon
       if (iflag_coupl.eq.0) then
          if (debut.and.prt_level.gt.9) &
               WRITE(lunout,*)'ALE et ALP imposes'
          do i = 1,klon
             !on ne couple que ale
             !           ALE(i) = max(ale_wake(i),Ale_bl(i))
             ALE(i) = max(ale_wake(i),ale_bl_prescr)
             !on ne couple que alp
             !           ALP(i) = alp_wake(i) + Alp_bl(i)
             ALP(i) = alp_wake(i) + alp_bl_prescr
          enddo
       else
          IF(prt_level>9)WRITE(lunout,*)'ALE et ALP couples au thermique'
          !         do i = 1,klon
          !             ALE(i) = max(ale_wake(i),Ale_bl(i))
          ! avant        ALP(i) = alp_wake(i) + Alp_bl(i)
          !             ALP(i) = alp_wake(i) + Alp_bl(i) + alp_offset ! modif sb
          !         write(20,*)'ALE',ALE(i),Ale_bl(i),ale_wake(i)
          !         write(21,*)'ALP',ALP(i),Alp_bl(i),alp_wake(i)
          !         enddo

          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Modif FH 2010/04/27. Sans doute temporaire.
          ! Deux options pour le alp_offset : constant si >?? 0 ou
          ! proportionnel ??a w si <0
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Estimation d'une vitesse verticale effective pour ALP
          if (1==0) THEN
             www(1:klon)=0.
             do k=2,klev-1
                do i=1,klon
                   www(i)=max(www(i),-omega(i,k)*RD*t_seri(i,k) &
                        /(RG*paprs(i,k)) *zw2(i,k)*zw2(i,k))
                   ! if (paprs(i,k)>pbase(i)) then
                   ! calcul approche de la vitesse verticale en m/s
                   !  www(i)=max(www(i),-omega(i,k)*RD*temp(i,k)/(RG*paprs(i,k))
                   !             endif
                   !   Le 0.1 est en gros H / ps = 1e4 / 1e5
                enddo
             enddo
             do i=1,klon
                if (www(i)>0. .and. ale_bl(i)>0. ) www(i)=www(i)/ale_bl(i)
             enddo
          ENDIF


          do i = 1,klon
             ALE(i) = max(ale_wake(i),Ale_bl(i))
             !cc nrlmd le 10/04/2012----------Stochastic triggering------------
             if (iflag_trig_bl.ge.1) then
                ALE(i) = max(ale_wake(i),Ale_bl_trig(i))
             endif
             !cc fin nrlmd le 10/04/2012
             if (alp_offset>=0.) then
                ALP(i) = alp_wake(i) + Alp_bl(i) + alp_offset ! modif sb
             else
                abort_message ='Ne pas passer la car www non calcule'
                CALL abort_physic (modname,abort_message,1)

                ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !                                _                  _
                ! Ajout d'une composante 3 * A * w w'2 a w'3 avec
                ! w=www : w max sous pbase ou A est la fraction
                ! couverte par les ascendances w' on utilise le fait
                ! que A * w'3 = ALP et donc A * w'2 ~ ALP / sqrt(ALE)
                ! (on ajoute 0.1 pour les singularites)
                ALP(i)=alp_wake(i)*(1.+3.*www(i)/( sqrt(ale_wake(i))+0.1) ) &
                     +alp_bl(i)  *(1.+3.*www(i)/( sqrt(ale_bl(i))  +0.1) )
                !    ALP(i)=alp_wake(i)+Alp_bl(i)+alp_offset*min(omega(i,6),0.)
                !             if (alp(i)<0.) then
                !                print*,'ALP ',alp(i),alp_wake(i) &
                !                     ,Alp_bl(i),alp_offset*min(omega(i,6),0.)
                !             endif
             endif
          enddo
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       endif
       do i=1,klon
          if (alp(i)>alp_max) then
             IF(prt_level>9)WRITE(lunout,*)                             &
                  'WARNING SUPER ALP (seuil=',alp_max, &
                  '): i, alp, alp_wake,ale',i,alp(i),alp_wake(i),ale(i)
             alp(i)=alp_max
          endif
          if (ale(i)>ale_max) then
             IF(prt_level>9)WRITE(lunout,*)                             &
                  'WARNING SUPER ALE (seuil=',ale_max, &
                  '): i, alp, alp_wake,ale',i,ale(i),ale_wake(i),alp(i)
             ale(i)=ale_max
          endif
       enddo

       !fin calcul ale et alp
       !=======================================================================


  RETURN
  END

