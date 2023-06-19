!
! $Id $
!
      SUBROUTINE sw_case_williamson91_6_loc(vcov,ucov,teta,masse,ps)

c=======================================================================
c
c   Author:    Thomas Dubos      original: 26/01/2010
c   -------
c
c   Subject:
c   ------
c   Realise le cas-test 6 de Williamson et al. (1991) : onde de Rossby-Haurwitz
c
c   Method:
c   --------
c
c   Interface:
c   ----------
c
c      Input:
c      ------
c
c      Output:
c      -------
c
c=======================================================================
      USE parallel_lmdz
      USE comconst_mod, ONLY: cpp, omeg, rad
      USE comvert_mod, ONLY: ap, bp, preff

      IMPLICIT NONE
c-----------------------------------------------------------------------
c   Declararations:
c   ---------------

      include "dimensions.h"
      include "paramet.h"
      include "comgeom.h"
      include "iniprint.h"

c   Arguments:
c   ----------

c   variables dynamiques
      REAL vcov(ijb_v:ije_v,llm),ucov(ijb_u:ije_u,llm) ! vents covariants
      REAL teta(ijb_u:ije_u,llm)                 ! temperature potentielle
      REAL ps(ijb_u:ije_u)                       ! pression  au sol
      REAL masse(ijb_u:ije_u,llm)                ! masse d'air
      REAL phis(ijb_u:ije_u)                     ! geopotentiel au sol

c   Local:
c   ------

      real,allocatable :: ucov_glo(:,:)
      real,allocatable :: vcov_glo(:,:)
      real,allocatable :: teta_glo(:,:)
      real,allocatable :: masse_glo(:,:)
      real,allocatable :: ps_glo(:)

!      REAL p (ip1jmp1,llmp1  )               ! pression aux interfac.des couches
!      REAL pks(ip1jmp1)                      ! exner au  sol
!      REAL pk(ip1jmp1,llm)                   ! exner au milieu des couches
!      REAL pkf(ip1jmp1,llm)                  ! exner filt.au milieu des couches
!      REAL alpha(ip1jmp1,llm),beta(ip1jmp1,llm)

      real,allocatable :: p(:,:)
      real,allocatable :: pks(:)
      real,allocatable :: pk(:,:)
      real,allocatable :: pkf(:,:)
      real,allocatable :: alpha(:,:),beta(:,:)

      REAL :: sinth,costh,costh2, Ath,Bth,Cth, lon,dps
      INTEGER i,j,ij

      REAL, PARAMETER    :: rho=1 ! masse volumique de l'air (arbitraire)
      REAL, PARAMETER    :: K    = 7.848e-6  ! K = \omega
      REAL, PARAMETER    :: gh0  = 9.80616 * 8e3 
      INTEGER, PARAMETER :: R0=4, R1=R0+1, R2=R0+2         ! mode 4
c NB : rad = 6371220 dans W91 (6371229 dans LMDZ)
c      omeg = 7.292e-5 dans W91 (7.2722e-5 dans LMDZ)


       ! allocate (global) arrays
       allocate(vcov_glo(ip1jm,llm))
       allocate(ucov_glo(ip1jmp1,llm))
       allocate(teta_glo(ip1jmp1,llm))
       allocate(ps_glo(ip1jmp1))
       allocate(masse_glo(ip1jmp1,llm))

       allocate(p(ip1jmp1,llmp1))
       allocate(pks(ip1jmp1))
       allocate(pk(ip1jmp1,llm))
       allocate(pkf(ip1jmp1,llm))
       allocate(alpha(ip1jmp1,llm))
       allocate(beta(ip1jmp1,llm))
 
      IF(0==0) THEN
!c Williamson et al. (1991) : onde de Rossby-Haurwitz
         teta_glo(:,:) = preff/rho/cpp
!c geopotentiel (pression de surface)
         do j=1,jjp1
            costh2 = cos(rlatu(j))**2
            Ath = (R0+1)*(costh2**2) + (2*R0*R0-R0-2)*costh2 - 2*R0*R0
            Ath = .25*(K**2)*(costh2**(R0-1))*Ath
            Ath = .5*K*(2*omeg+K)*costh2 + Ath 
            Bth = (R1*R1+1)-R1*R1*costh2
            Bth = 2*(omeg+K)*K/(R1*R2) * (costh2**(R0/2))*Bth
            Cth = R1*costh2 - R2
            Cth = .25*K*K*(costh2**R0)*Cth
            do i=1,iip1
               ij=(j-1)*iip1+i
               lon = rlonv(i)
               dps = Ath + Bth*cos(R0*lon) + Cth*cos(2*R0*lon)
               ps_glo(ij) = rho*(gh0 + (rad**2)*dps)
            enddo
         enddo
!         write(lunout,*) 'W91 ps', MAXVAL(ps), MINVAL(ps)
c vitesse zonale ucov
         do j=1,jjp1
            costh  = cos(rlatu(j))
            costh2 = costh**2
            Ath = rad*K*costh
            Bth = R0*(1-costh2)-costh2
            Bth = rad*K*Bth*(costh**(R0-1))
            do i=1,iip1
               ij=(j-1)*iip1+i
               lon = rlonu(i)
               ucov_glo(ij,1) = (Ath + Bth*cos(R0*lon))
            enddo
         enddo
!         write(lunout,*) 'W91 u', MAXVAL(ucov(:,1)), MINVAL(ucov(:,1))
         ucov_glo(:,1)=ucov_glo(:,1)*cu
c vitesse meridienne vcov
         do j=1,jjm
            sinth  = sin(rlatv(j))
            costh  = cos(rlatv(j))
            Ath = -rad*K*R0*sinth*(costh**(R0-1))
            do i=1,iip1
               ij=(j-1)*iip1+i
               lon = rlonv(i)
               vcov_glo(ij,1) = Ath*sin(R0*lon)
            enddo
         enddo
         write(lunout,*) 'W91 v', MAXVAL(vcov(:,1)), MINVAL(vcov(:,1))
         vcov_glo(:,1)=vcov_glo(:,1)*cv
        
c         ucov_glo=0
c         vcov_glo=0
      ELSE
c test non-tournant, onde se propageant en latitude
         do j=1,jjp1
            do i=1,iip1
               ij=(j-1)*iip1+i
               ps_glo(ij) = 1e5*(1 + .1*exp(-100*(1+sin(rlatu(j)))**2))
            enddo
         enddo
         
c     rho = preff/(cpp*teta)
         teta_glo(:,:) = .01*preff/cpp   ! rho = 100 ; phi = ps/rho = 1e3 ; c=30 m/s = 2600 km/j = 23 degres / j
         ucov_glo(:,:)=0.
         vcov_glo(:,:)=0.
      END IF      
      
      CALL pression ( ip1jmp1, ap, bp, ps_glo, p       )
      CALL massdair(p,masse_glo)

      ! copy data from global array to local array:
      teta(ijb_u:ije_u,:)=teta_glo(ijb_u:ije_u,:)
      ucov(ijb_u:ije_u,:)=ucov_glo(ijb_u:ije_u,:)
      vcov(ijb_v:ije_v,:)=vcov_glo(ijb_v:ije_v,:)
      masse(ijb_u:ije_u,:)=masse_glo(ijb_u:ije_u,:)
      ps(ijb_u:ije_u)=ps_glo(ijb_u:ije_u)

      ! cleanup
      deallocate(teta_glo)
      deallocate(ucov_glo)
      deallocate(vcov_glo)
      deallocate(masse_glo)
      deallocate(ps_glo)
      deallocate(p)
      deallocate(pks)
      deallocate(pk)
      deallocate(pkf)
      deallocate(alpha)
      deallocate(beta)

      END
c-----------------------------------------------------------------------
