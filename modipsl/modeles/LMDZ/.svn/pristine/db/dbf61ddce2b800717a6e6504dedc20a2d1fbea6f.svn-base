      SUBROUTINE bcscav_spl(pdtime,flxr,flxs,alpha_r,alpha_s,x,dx) 

      USE dimphy
      IMPLICIT NONE 
c=====================================================================
c Objet : below-cloud scavenging of tracers
c Date : september 1999
c Auteur: O. Boucher (LOA) 
c=====================================================================
c
#include "dimensions.h"
#include "chem.h"
#include "YOMCST.h"
#include "YOECUMF.h"
c
      REAL pdtime, alpha_r, alpha_s, R_r, R_s 
      PARAMETER (R_r=0.001)          !--mean raindrop radius (m)
      PARAMETER (R_s=0.001)          !--mean snow crystal radius (m)
      REAL flxr(klon,klev)         ! liquid precipitation rate (kg/m2/s)
      REAL flxs(klon,klev)         ! solid  precipitation rate (kg/m2/s)
      REAL flxr_aux(klon,klev+1)
      REAL flxs_aux(klon,klev+1)
      REAL x(klon,klev)              ! q de traceur  
      REAL dx(klon,klev)             ! tendance de traceur
c
c--variables locales      
      INTEGER i, k 
      REAL pr, ps, ice, water
c
c------------------------------------------
c
! NHL
! Auxiliary variables defined to deal with the fact that precipitation
! fluxes are defined on klev levels only.
! NHL
!
      flxr_aux(:,klev+1)=0.0
      flxs_aux(:,klev+1)=0.0
      flxr_aux(:,1:klev)=flxr(:,:)
      flxs_aux(:,1:klev)=flxs(:,:)
!
      DO k=1, klev
      DO i=1, klon 
       pr=0.5*(flxr_aux(i,k)+flxr_aux(i,k+1))
       ps=0.5*(flxs_aux(i,k)+flxs_aux(i,k+1))
       water=pr*alpha_r/R_r/rho_water
       ice=ps*alpha_s/R_s/rho_ice
       dx(i,k)=-3./4.*x(i,k)*pdtime*(water+ice)
ctmp       dx(i,k)=-3./4.*x(i,k)*pdtime* 
ctmp     .         (pr*alpha_r/R_r/rho_water+ps*alpha_s/R_s/rho_ice)
      ENDDO 
      ENDDO
c
      RETURN 
      END
