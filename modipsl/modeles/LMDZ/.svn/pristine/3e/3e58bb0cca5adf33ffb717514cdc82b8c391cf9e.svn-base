!
! $Id: thermcell_plume.F90 2311 2015-06-25 07:45:24Z emillour $
!
      SUBROUTINE thermcell_alim(flag,ngrid,klev,ztv,d_temp,zlev,alim_star,lalim)
IMPLICIT NONE

!--------------------------------------------------------------------------
! FH : 2015/11/06
! thermcell_alim: calcule la distribution verticale de l'alimentation 
! laterale a la base des panaches thermiques
!--------------------------------------------------------------------------

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "thermcell.h"

!      fort(10) ptimestep,ztv,zthl,po,zl,rhobarz,zlev,pplev,pphi,zpspsk,f0
      INTEGER, INTENT(IN) :: ngrid,klev
      REAL, INTENT(IN) :: ztv(ngrid,klev)
      REAL, INTENT(IN) :: d_temp(ngrid)
      REAL, INTENT(IN) :: zlev(ngrid,klev+1)
      REAL, INTENT(OUT) :: alim_star(ngrid,klev)
      INTEGER, INTENT(OUT) :: lalim(ngrid)
      INTEGER, INTENT(IN) :: flag

      REAL :: alim_star_tot(ngrid),zi(ngrid),zh(ngrid)
      REAL :: zlay(ngrid,klev)
      REAL ztv_parcel

      INTEGER ig,l

      REAL h,z,falim
      falim(h,z)=0.2*((z-h)**5+h**5)


!===================================================================

   lalim(:)=1
   alim_star_tot(:)=0.

!-------------------------------------------------------------------------
! Definition de l'alimentation a l'origine dans thermcell_init
!-------------------------------------------------------------------------
   IF (flag==0) THEN ! CMIP5 version
      do l=1,klev-1
         do ig=1,ngrid
            if (ztv(ig,l)> ztv(ig,l+1) .and. ztv(ig,1)>=ztv(ig,l) ) then
               alim_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.)  &
     &                       *sqrt(zlev(ig,l+1)) 
               lalim(ig)=l+1
               alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
            endif
         enddo
      enddo
      do l=1,klev
         do ig=1,ngrid 
            if (alim_star_tot(ig) > 1.e-10 ) then
               alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
            endif
         enddo
      enddo
      alim_star_tot(:)=1.

!-------------------------------------------------------------------------
! Nouvelle definition avec possibilite d'introduire un DT en surface
! On suppose que la forme du profile d'alimentation scale avec la hauteur
! d'inversion calculée avec une particule partant de la premieere couche

! Fonction  f(z) = z ( h - z ) , avec h = zi/3
! On utilise l'integralle
! Int_0^z f(z') dz' = z^2 ( h/2 - z/3 ) = falim(h,z)
! Pour calculer l'alimentation des couches
!-------------------------------------------------------------------------
   ELSE
! Computing inversion height zi and zh=zi/3.
      zi(:)=0.
! Il faut recalculer zlay qui n'est pas dispo dans thermcell_plume
! A changer eventuellement.
      do l=1,klev
         zlay(:,l)=0.5*(zlev(:,l)+zlev(:,l+1))
      enddo

      do l=klev-1,1,-1
         do ig=1,ngrid
            ztv_parcel=ztv(ig,1)+d_temp(ig)
            if (ztv_parcel<ztv(ig,l+1)) lalim(ig)=l
         enddo
      enddo

      do ig=1,ngrid
         l=lalim(ig)
         IF (l==1) THEN
            zi(ig)=0.
         ELSE
            ztv_parcel=ztv(ig,1)+d_temp(ig)
            zi(ig)=zlay(ig,l)+(zlay(ig,l+1)-zlay(ig,l))/(ztv(ig,l+1)-ztv(ig,l))*(ztv_parcel-ztv(ig,l))
         ENDIF
      enddo

      zh(:)=zi(:)/2.
      alim_star_tot(:)=0.
      alim_star(:,:)=0.
      lalim(:)=0
      do l=1,klev-1
         do ig=1,ngrid
            IF (zh(ig)==0.) THEN
               alim_star(ig,l)=0.
               lalim(ig)=1
            ELSE IF (zlev(ig,l+1)<=zh(ig)) THEN
               alim_star(ig,l)=(falim(zh(ig),zlev(ig,l+1))-falim(zh(ig),zlev(ig,l)))/falim(zh(ig),zh(ig))
               lalim(ig)=l
            ELSE IF (zlev(ig,l)<=zh(ig)) THEN
               alim_star(ig,l)=(falim(zh(ig),zh(ig))-falim(zh(ig),zlev(ig,l)))/falim(zh(ig),zh(ig))
               lalim(ig)=l
            ELSE
               alim_star(ig,l)=0.
            ENDIF
         ENDDO
         alim_star_tot(:)=alim_star_tot(:)+alim_star(:,l)
      ENDDO
      IF (ngrid==1) print*,'NEW ALIM CALCUL DE ZI ',alim_star_tot,lalim,zi,zh
      alim_star_tot(:)=1.

   ENDIF


RETURN
END
