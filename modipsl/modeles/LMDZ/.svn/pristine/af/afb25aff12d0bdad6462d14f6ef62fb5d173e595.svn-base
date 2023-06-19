      SUBROUTINE inscav_spl(pdtime,it,masse,henry,kk,qliq,
     .                   flxr,flxs,zrho,zdz,t,x,
     .                   his_dh) 
      USE dimphy
      IMPLICIT NONE 
c=====================================================================
c Objet : depot humide de traceurs
c Date : mars 1998
c Auteur: O. Boucher (LOA) 
c=====================================================================
c
#include "dimensions.h"
#include "chem.h"
#include "YOMCST.h"
#include "YOECUMF.h" 
c
      INTEGER it
      REAL pdtime              ! pas de temps (s)
      REAL masse               ! molar mass (except for BC/OM/IF/DUST=Nav)
      REAL henry               ! constante de Henry en mol/l/atm
      REAL kk                  ! coefficient de dependence en T (K)
      REAL qliq                ! contenu en eau liquide dans le nuage (kg/kg)
!      REAL flxr(klon,klev+1)   ! flux precipitant de pluie
!      REAL flxs(klon,klev+1)   ! flux precipitant de neige
      REAL flxr(klon,klev)   ! flux precipitant de pluie   ! Titane
      REAL flxs(klon,klev)   ! flux precipitant de neige   ! Titane
      REAL flxr_aux(klon,klev+1)
      REAL flxs_aux(klon,klev+1)
      REAL zrho(klon,klev)
      REAL zdz(klon,klev)
      REAL t(klon,klev)
      REAL x(klon,klev)        ! q de traceur  
      REAL his_dh(klon)        ! tendance de traceur integre verticalement
c
c--variables locales      
      INTEGER i, k
c
      REAL dx      ! tendance de traceur
      REAL f_a     !--rapport de la phase aqueuse a la phase gazeuse
      REAL beta    !--taux de conversion de l'eau en pluie
      REAL henry_t !--constante de Henry a T t  (mol/l/atm)
      REAL scav(klon,klev)    !--fraction aqueuse du constituant
      REAL K1, K2, ph, frac 
      REAL frac_gas, frac_aer !-cste pour la reevaporation 
      PARAMETER (ph=5., frac_gas=1.0, frac_aer=0.5)
c---cste de dissolution pour le depot humide
      REAL frac_fine_scav,frac_coar_scav
c---added by nhl
      REAL aux_cte

      PARAMETER (frac_fine_scav=0.7)
      PARAMETER (frac_coar_scav=0.7)

c--101.325  m3/l x Pa/atm
c--R        Pa.m3/mol/K
c
c------------------------------------------
c
cnhl      IF (it.EQ.2.OR.it.EQ.3) THEN !--aerosol  ! AS IT WAS FIRST
      IF (it.EQ.2.OR.it.EQ.3.OR.it.EQ.4) THEN !--aerosol
        frac=frac_aer
      ELSE                                                !--gas
        frac=frac_gas
      ENDIF
c
      IF (it.EQ.1) THEN 
      DO k=1, klev
      DO i=1, klon 
        henry_t=henry*exp(-kk*(1./298.-1./t(i,k)))    !--mol/l/atm
        K1=1.2e-2*exp(-2010*(1/298.-1/t(i,k)))
        K2=6.6e-8*exp(-1510*(1/298.-1/t(i,k)))
        henry_t=henry_t*(1 + K1/10.**(-ph) + K1*K2/(10.**(-ph))**2)
        f_a=henry_t/101.325*R*t(i,k)*qliq*zrho(i,k)/rho_water
        scav(i,k)=f_a/(1.+f_a)
      ENDDO 
      ENDDO
      ELSEIF (it.EQ.2) THEN
      DO k=1, klev
      DO i=1, klon 
        scav(i,k)=frac_fine_scav
      ENDDO 
      ENDDO
      ELSEIF (it.EQ.3) THEN
      DO k=1, klev
      DO i=1, klon 
        scav(i,k)=frac_coar_scav
      ENDDO 
      ENDDO 
      ELSEIF (it.EQ.4) THEN
      DO k=1, klev
      DO i=1, klon 
        scav(i,k)=frac_coar_scav
      ENDDO 
      ENDDO 
      ELSE 
        PRINT *,'it non pris en compte'
        STOP
      ENDIF
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
      DO k=klev, 1, -1 
      DO i=1, klon 
c--scavenging
        beta=flxr_aux(i,k)-flxr_aux(i,k+1)+flxs_aux(i,k)-flxs_aux(i,k+1)
        beta=beta/zdz(i,k)/qliq/zrho(i,k)
        beta=MAX(0.0,beta)
        dx=x(i,k)*(exp(-scav(i,k)*beta*pdtime)-1.)
        x(i,k)=x(i,k)+dx
        his_dh(i)=his_dh(i)-dx/RNAVO*
     .            masse*1.e3*1.e6*zdz(i,k)/pdtime !--mgS/m2/s
c--reevaporation
        beta=flxr_aux(i,k)-flxr_aux(i,k+1)+flxs_aux(i,k)-flxs_aux(i,k+1)
        IF (beta.LT.0.) beta=beta/(flxr_aux(i,k+1)+flxs_aux(i,k+1))
        IF (flxr_aux(i,k)+flxs_aux(i,k).EQ.0) THEN  !--reevaporation totale
          beta=MIN(MAX(0.0,-beta),1.0)
        ELSE                          !--reevaporation non totale pour aerosols
          !print *,'FRAC USED IN INSCAV_SPL'
          beta=MIN(MAX(0.0,-beta)*frac,1.0)
        ENDIF
        dx=beta*his_dh(i)*RNAVO/masse/1.e3/1.e6/zdz(i,k)*pdtime !ORIG LINE
! funny line for TL/AD
! AD test works without (x) and for xd = dxd*1.e5 : 2.79051851638 times the 0.
! AD test does not work with the line : 754592404.083 times the 0.
! problem seems to be linked to the largest dx wrt x
!       x(i, k) = x(i, k) + dx
!        x(i, k) = x(i, k) + dx         ! THIS LINE WAS COMMENTED OUT ORIGINALY !nhl
        his_dh(i)=(1.-beta)*his_dh(i)
      ENDDO 
      ENDDO
c
      RETURN 
      END
