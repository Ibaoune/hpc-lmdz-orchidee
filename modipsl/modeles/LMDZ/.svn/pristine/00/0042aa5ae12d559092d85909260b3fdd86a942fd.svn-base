!***************************************************************************************
! tend_to_tke.F90
!*************
!
! Subroutine that adds a tendency on the TKE created by the 
! fluxes of momentum retrieved from the wind speed tendencies
! of the physics.
! 
! The basic concept is the following:
! the TKE equation writes  de/dt = -u'w' du/dz -v'w' dv/dz +g/theta dtheta/dz +......
!
!
! We expect contributions to the term u'w' and v'w' that do not come from the Yamada
! scheme, for instance: gravity waves, drag from high vegetation..... These contributions
! need to be accounted for.
! we explicitely calculate the fluxes, integrating the wind speed 
!                        tendency from the top of the atmospher
!
!
!
! contacts: Frederic Hourdin, Etienne Vignon
!
! History:
!---------
! - 1st redaction, Etienne, 15/10/2016
! Ajout des 4 sous surfaces pour la tke
! on sort l'ajout des tendances du if sur les deux cas, pour ne pas
! dupliuqer les lignes
! on enleve le pas de temps qui disprait dans les calculs
!
!
!**************************************************************************************

 SUBROUTINE tend_to_tke(dt,plev,exner,temp,windu,windv,dt_a,du_a,dv_a,pctsrf,tke)

 USE dimphy, ONLY: klon, klev
 USE indice_sol_mod, ONLY: nbsrf

IMPLICIT NONE
#include "YOMCST.h"

! Declarations
!==============


! Inputs
!-------
  REAL dt                   ! Time step [s]
  REAL plev(klon,klev+1)    ! inter-layer pressure [Pa]
  REAL temp(klon,klev)      ! temperature [K], grid-cell average or for a one subsurface
  REAL windu(klon,klev)     ! zonal wind [m/s], grid-cell average or for a one subsurface
  REAL windv(klon,klev)     ! meridonal wind [m/s], grid-cell average or for a one subsurface
  REAL exner(klon,klev)     ! Fonction d'Exner = T/theta
  REAL dt_a(klon,klev)      ! Temperature tendency [K], grid-cell average or for a one subsurface
  REAL du_a(klon,klev)      ! Zonal wind speed tendency [m/s], grid-cell average or for a one subsurface
  REAL dv_a(klon,klev)      ! Meridional wind speed tendency [m/s], grid-cell average or for a one subsurface
  REAL pctsrf(klon,nbsrf+1)       ! Turbulent Kinetic energy [m2/s2], grid-cell average or for a subsurface

! Inputs/Outputs
!---------------
  REAL tke(klon,klev+1,nbsrf+1)       ! Turbulent Kinetic energy [m2/s2], grid-cell average or for a subsurface


! Local
!-------


  INTEGER i,k,isrf                 ! indices
  REAL    masse(klon,klev)          ! mass in the layers [kg/m2]
  REAL    unsmasse(klon,klev+1)     ! linear mass in the layers [kg/m2]
  REAL    flux_rhotw(klon,klev+1)   ! flux massique de tempe. pot. rho*u'*theta'
  REAL    flux_rhouw(klon,klev+1)   ! flux massique de quantit?? de mouvement rho*u'*w' [kg/m/s2]
  REAL    flux_rhovw(klon,klev+1)   ! flux massique de quantit?? de mouvement rho*v'*w' [kg/m/s2]
  REAL    tendt(klon,klev)        ! new temperature tke tendency [m2/s2/s]
  REAL    tendu(klon,klev)        ! new zonal tke tendency [m2/s2/s]
  REAL    tendv(klon,klev)        ! new meridonal tke tendency [m2/s2/s]
  



! First calculations:
!=====================

      unsmasse(:,:)=0.
      DO k=1,klev
         masse(:,k)=(plev(:,k)-plev(:,k+1))/RG
         unsmasse(:,k)=unsmasse(:,k)+0.5/masse(:,k)
         unsmasse(:,k+1)=unsmasse(:,k+1)+0.5/masse(:,k)
      END DO

      tendu(:,:)=0.0
      tendv(:,:)=0.0

! Method 1: Calculation of fluxes using a downward integration
!============================================================


 
! Flux calculation

 flux_rhotw(:,klev+1)=0.
 flux_rhouw(:,klev+1)=0.
 flux_rhovw(:,klev+1)=0.

   DO k=klev,1,-1
      flux_rhotw(:,k)=flux_rhotw(:,k+1)+masse(:,k)*dt_a(:,k)/exner(:,k)
      flux_rhouw(:,k)=flux_rhouw(:,k+1)+masse(:,k)*du_a(:,k)
      flux_rhovw(:,k)=flux_rhovw(:,k+1)+masse(:,k)*dv_a(:,k)
   ENDDO


! TKE update:

   DO k=2,klev
      tendt(:,k)=-flux_rhotw(:,k)*(exner(:,k)-exner(:,k-1))*unsmasse(:,k)*RCPD
      tendu(:,k)=-flux_rhouw(:,k)*(windu(:,k)-windu(:,k-1))*unsmasse(:,k)
      tendv(:,k)=-flux_rhovw(:,k)*(windv(:,k)-windv(:,k-1))*unsmasse(:,k)
   ENDDO
   tendt(:,1)=-flux_rhotw(:,1)*(exner(:,1)-1.)*unsmasse(:,1)*RCPD
   tendu(:,1)=-1.*flux_rhouw(:,1)*windu(:,1)*unsmasse(:,1)
   tendv(:,1)=-1.*flux_rhovw(:,1)*windv(:,1)*unsmasse(:,1)


 DO isrf=1,nbsrf
    DO k=1,klev
       DO i=1,klon
          IF (pctsrf(i,isrf)>0.) THEN
            tke(i,k,isrf)= tke(i,k,isrf)+tendu(i,k)+tendv(i,k)+tendt(i,k)
            tke(i,k,isrf)= max(tke(i,k,isrf),1.e-10)
          ENDIF
       ENDDO
    ENDDO
 ENDDO


!  IF (klon==1) THEN
!  CALL iophys_ecrit('u',klev,'u','',windu)
!  CALL iophys_ecrit('v',klev,'v','',windu)
!  CALL iophys_ecrit('t',klev,'t','',temp)
!  CALL iophys_ecrit('tke1',klev,'tke1','',tke(:,1:klev,1))
!  CALL iophys_ecrit('tke2',klev,'tke2','',tke(:,1:klev,2))
!  CALL iophys_ecrit('tke3',klev,'tke3','',tke(:,1:klev,3))
!  CALL iophys_ecrit('tke4',klev,'tke4','',tke(:,1:klev,4))
!  CALL iophys_ecrit('theta',klev,'theta','',temp/exner)
!  CALL iophys_ecrit('Duv',klev,'Duv','',tendu(:,1:klev)+tendv(:,1:klev))
!  CALL iophys_ecrit('Dt',klev,'Dt','',tendt(:,1:klev))
!  ENDIF

 END SUBROUTINE tend_to_tke
