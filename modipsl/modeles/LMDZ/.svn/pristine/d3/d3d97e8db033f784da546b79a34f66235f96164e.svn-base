!
! $Id$
!
  SUBROUTINE freinage(knon, uu, vv,  &
       tt,veget,lai, height,ypaprs,ypplay,drag_pro,d_u,d_v)

    !ONLINE:
    use dimphy, only: klon, klev
!    USE control, ONLY: nvm
!    USE indice_sol_mod, only : nvm_orch

    IMPLICIT NONE


    include "YOMCST.h"
    include "clesphys.h"
    include "YOEGWD.h"
!FC 
    include "dimpft.h"
    include "compbl.h"

    ! 0. DECLARATIONS:

    ! 0.1 INPUTS

    REAL, DIMENSION(klon,klev), INTENT(IN)         :: ypplay
    REAL, DIMENSION(klon,klev+1), INTENT(IN)       :: ypaprs


     REAL, DIMENSION(klon, klev), INTENT(IN)     :: uu
     REAL, DIMENSION(klon, klev), INTENT(IN)     :: vv
     REAL, DIMENSION(klon, klev), INTENT(IN)     :: tt
     REAL, DIMENSION(klon,nvm_lmdz), INTENT(IN)          :: veget,lai
     REAL, DIMENSION(klon,nvm_lmdz), INTENT(IN)          :: height

     REAL, DIMENSION(klon,klev)         :: wind
     REAL, DIMENSION(klon, klev)        :: yzlay
     INTEGER knon

    ! 0.2 OUTPUTS

      REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_v        ! change in v
      REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_u        ! change in v
    !knon nombre de points concernes 
      REAL, DIMENSION(klon,klev)         :: sumveg        ! change in v
    
     REAL,  DIMENSION(klon,klev), INTENT(OUT)          :: drag_pro
    ! (KLON, KLEV) tendencies on winds


    INTEGER k,jv,i


!FCCCC    REAL Cd_frein

    ! 0.3.1 LOCAL VARIABLE


    !-----------------------------------------------------------------

    ! 1. INITIALISATIONS

    
!    Cd_frein = 7.5E-2 ! (0.075) ! Drag from MASSON 2009
!FC ESSAI
!    Cd_frein = 1.5E-2 ! (0.075) ! Drag from MASSON 2009
!    Cd_frein = 0.005 ! (0.075) ! Drag from MASSON 2009

! initialisation 
      d_u(:,:) =0.
      d_v(:,:) =0.
      drag_pro(:,:) =0.
      sumveg(:,:) =0.
!!        print*, "Cd_frein" , Cd_frein
      
       wind(:,:)= sqrt(uu(:,:)*uu(:,:)+vv(:,:)*vv(:,:))

       yzlay(1:knon,1)= &
            RD*tt(1:knon,1)/(0.5*(ypaprs(1:knon,1)+ypplay(1:knon,1))) &
            *(ypaprs(1:knon,1)-ypplay(1:knon,1))/RG
       DO k=2,klev
             yzlay(1:knon,k)= &
                  yzlay(1:knon,k-1)+RD*0.5*(tt(1:knon,k-1)+tt(1:knon,k)) &
                  /ypaprs(1:knon,k)*(ypplay(1:knon,k-1)-ypplay(1:knon,k))/RG
       END DO

!    verifier les indexes ..... 
!!       print*, " calcul de drag_pro FC "
   
      do k= 1,klev

      do jv=2,nvm_lmdz   !   (on peut faire 9 ?)

      do i=1,knon

      sumveg(i,k)= sumveg(i,k)+ veget(i,jv)

!      if  ( (height(i,jv) .gt. yzlay(i,k)) .AND. (height(i,jv) .gt. 0.1) .and. LAI(i,jv).gt.0. ) then                     
      if  ( (height(i,jv) .gt. yzlay(i,k)) .AND. (height(i,jv) .gt. 0.1) ) then                     
!FC attention veut on le test sur le LAI ?
         if (ifl_pbltree.eq.1) then
      drag_pro(i,k)= drag_pro(i,k)+ &
      veget(i,jv)
          elseif (ifl_pbltree.eq.2) then
      drag_pro(i,k)= drag_pro(i,k)+ &
      6*LAI(i,jv)*veget(i,jv)*( yzlay(i,k)*(height(i,jv)-yzlay(i,k))/(height(i,jv)*height(i,jv)+ 0.01))
          elseif (ifl_pbltree.eq.3) then
      drag_pro(i,k)= drag_pro(i,k)+ &
      veget(i,jv)*( yzlay(i,k)*(height(i,jv)-yzlay(i,k))/(height(i,jv)*height(i,jv)+ 0.01))
          elseif (ifl_pbltree.eq.0) then
          drag_pro(i,k)=0.0
           endif
      else
      drag_pro(i,k)= drag_pro(i,k)
      endif


      enddo
      enddo
     enddo 
      do k=1,klev
        where (sumveg(1:knon,k) > 0.05 ) 
!        drag_pro(1:knon,k)=Cd_frein*drag_pro(1:knon,k)/sumveg(1:knon,k)
        drag_pro(1:knon,k)=Cd_frein*drag_pro(1:knon,k)
        elsewhere
        drag_pro(1:knon,k)=0.0
       endwhere
        d_u(1:knon,k) =(-1)*drag_pro(1:knon,k)*uu(1:knon,k)*wind(1:knon,k)
        d_v(1:knon,k) =(-1)*drag_pro(1:knon,k)*vv(1:knon,k)*wind(1:knon,k)
      enddo
      return

 END SUBROUTINE freinage

