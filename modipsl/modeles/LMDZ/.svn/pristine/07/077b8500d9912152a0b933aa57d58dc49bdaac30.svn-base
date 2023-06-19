!
! $Header$
!
SUBROUTINE soil(ptimestep, indice, knon, snow, ptsrf, &
     ptsoil, pcapcal, pfluxgrd)
  
  USE dimphy
  USE mod_phys_lmdz_para
  USE indice_sol_mod
  USE print_control_mod, ONLY: lunout

  IMPLICIT NONE

!=======================================================================
!
!   Auteur:  Frederic Hourdin     30/01/92
!   -------
!
!   Object:  Computation of : the soil temperature evolution
!   -------                   the surfacic heat capacity "Capcal"
!                            the surface conduction flux pcapcal
!
!
!   Method: Implicit time integration
!   -------
!   Consecutive ground temperatures are related by:
!           T(k+1) = C(k) + D(k)*T(k)  (*)
!   The coefficients C and D are computed at the t-dt time-step.
!   Routine structure:
!   1) C and D coefficients are computed from the old temperature
!   2) new temperatures are computed using (*)
!   3) C and D coefficients are computed from the new temperature
!      profile for the t+dt time-step
!   4) the coefficients A and B are computed where the diffusive
!      fluxes at the t+dt time-step is given by
!             Fdiff = A + B Ts(t+dt)
!      or     Fdiff = F0 + Capcal (Ts(t+dt)-Ts(t))/dt
!             with F0 = A + B (Ts(t))
!                 Capcal = B*dt
!           
!   Interface:
!   ----------
!
!   Arguments:
!   ----------
!   ptimestep            physical timestep (s)
!   indice               sub-surface index
!   snow(klon)           snow
!   ptsrf(klon)          surface temperature at time-step t (K)
!   ptsoil(klon,nsoilmx) temperature inside the ground (K)
!   pcapcal(klon)        surfacic specific heat (W*m-2*s*K-1)
!   pfluxgrd(klon)       surface diffusive flux from ground (Wm-2)
!   
!=======================================================================
  INCLUDE "YOMCST.h"
  INCLUDE "dimsoil.h"
  INCLUDE "comsoil.h"
!-----------------------------------------------------------------------
! Arguments
! ---------
  REAL, INTENT(IN)                     :: ptimestep
  INTEGER, INTENT(IN)                  :: indice, knon
  REAL, DIMENSION(klon), INTENT(IN)    :: snow
  REAL, DIMENSION(klon), INTENT(IN)    :: ptsrf
  
  REAL, DIMENSION(klon,nsoilmx), INTENT(INOUT) :: ptsoil
  REAL, DIMENSION(klon), INTENT(OUT)           :: pcapcal
  REAL, DIMENSION(klon), INTENT(OUT)           :: pfluxgrd

!-----------------------------------------------------------------------
! Local variables
! ---------------
  INTEGER                             :: ig, jk, ierr
  REAL                                :: min_period,dalph_soil
  REAL, DIMENSION(nsoilmx)            :: zdz2
  REAL                                :: z1s
  REAL, DIMENSION(klon)               :: ztherm_i
  REAL, DIMENSION(klon,nsoilmx,nbsrf) :: C_coef, D_coef

! Local saved variables
! ---------------------
  REAL, SAVE                     :: lambda
!$OMP THREADPRIVATE(lambda)
  REAL, DIMENSION(nsoilmx), SAVE :: dz1, dz2
!$OMP THREADPRIVATE(dz1,dz2)
  LOGICAL, SAVE                  :: firstcall=.TRUE.
!$OMP THREADPRIVATE(firstcall)
    
!-----------------------------------------------------------------------
!   Depthts:
!   --------
  REAL fz,rk,fz1,rk1,rk2
  fz(rk)=fz1*(dalph_soil**rk-1.)/(dalph_soil-1.)


!-----------------------------------------------------------------------
! Calculation of some constants
! NB! These constants do not depend on the sub-surfaces
!-----------------------------------------------------------------------

  IF (firstcall) THEN
!-----------------------------------------------------------------------
!   ground levels 
!   grnd=z/l where l is the skin depth of the diurnal cycle:
!-----------------------------------------------------------------------

     min_period=1800. ! en secondes
     dalph_soil=2.    ! rapport entre les epaisseurs de 2 couches succ.
!$OMP MASTER
     IF (is_mpi_root) THEN
        OPEN(99,file='soil.def',status='old',form='formatted',iostat=ierr)
        IF (ierr == 0) THEN ! Read file only if it exists
           READ(99,*) min_period
           READ(99,*) dalph_soil
           WRITE(lunout,*)'Discretization for the soil model'
           WRITE(lunout,*)'First level e-folding depth',min_period, &
                '   dalph',dalph_soil
           CLOSE(99)
        END IF
     ENDIF
!$OMP END MASTER
     CALL bcast(min_period)
     CALL bcast(dalph_soil)

!   la premiere couche represente un dixieme de cycle diurne
     fz1=SQRT(min_period/3.14)
     
     DO jk=1,nsoilmx
        rk1=jk
        rk2=jk-1
        dz2(jk)=fz(rk1)-fz(rk2)
     ENDDO
     DO jk=1,nsoilmx-1
        rk1=jk+.5
        rk2=jk-.5
        dz1(jk)=1./(fz(rk1)-fz(rk2))
     ENDDO
     lambda=fz(.5)*dz1(1)
     WRITE(lunout,*)'full layers, intermediate layers (seconds)'
     DO jk=1,nsoilmx
        rk=jk
        rk1=jk+.5
        rk2=jk-.5
        WRITE(lunout,*)'fz=', &
             fz(rk1)*fz(rk2)*3.14,fz(rk)*fz(rk)*3.14
     ENDDO

     firstcall =.FALSE.
  END IF


!-----------------------------------------------------------------------
!   Calcul de l'inertie thermique a partir de la variable rnat.
!   on initialise a inertie_sic meme au-dessus d'un point de mer au cas 
!   ou le point de mer devienne point de glace au pas suivant
!   on corrige si on a un point de terre avec ou sans glace
!
!   iophys can be used to write the ztherm_i variable in a phys.nc file
!   and check the results; to do so, add "CALL iophys_ini" in physiq_mod
!   and add knindex to the list of inputs in all the calls to soil.F90
!   (and to soil.F90 itself !)
!-----------------------------------------------------------------------

  IF (indice == is_sic) THEN
     DO ig = 1, knon
        ztherm_i(ig)   = inertie_sic
     ENDDO
     IF (iflag_sic == 0) THEN
       DO ig = 1, knon
         IF (snow(ig) > 0.0) ztherm_i(ig)   = inertie_sno
       ENDDO
!      Otherwise sea-ice keeps the same inertia, even when covered by snow
     ENDIF
!    CALL iophys_ecrit_index('ztherm_sic', 1, 'ztherm_sic', 'USI', &
!      knon, knindex, ztherm_i)
  ELSE IF (indice == is_lic) THEN
     DO ig = 1, knon
        ztherm_i(ig)   = inertie_lic
        IF (snow(ig) > 0.0) ztherm_i(ig)   = inertie_sno
     ENDDO
!    CALL iophys_ecrit_index('ztherm_lic', 1, 'ztherm_lic', 'USI', &
!      knon, knindex, ztherm_i)
  ELSE IF (indice == is_ter) THEN
     DO ig = 1, knon
        ztherm_i(ig)   = inertie_sol
        IF (snow(ig) > 0.0) ztherm_i(ig)   = inertie_sno
     ENDDO
!    CALL iophys_ecrit_index('ztherm_ter', 1, 'ztherm_ter', 'USI', &
!      knon, knindex, ztherm_i)
  ELSE IF (indice == is_oce) THEN
     DO ig = 1, knon
!       This is just in case, but SST should be used by the model anyway
        ztherm_i(ig)   = inertie_sic
     ENDDO
!    CALL iophys_ecrit_index('ztherm_oce', 1, 'ztherm_oce', 'USI', &
!      knon, knindex, ztherm_i)
  ELSE
     WRITE(lunout,*) "valeur d indice non prevue", indice
     call abort_physic("soil", "", 1)
  ENDIF


!-----------------------------------------------------------------------
! 1)
! Calculation of Cgrf and Dgrd coefficients using soil temperature from 
! previous time step.
!
! These variables are recalculated on the local compressed grid instead 
! of saved in restart file.
!-----------------------------------------------------------------------
  DO jk=1,nsoilmx
     zdz2(jk)=dz2(jk)/ptimestep
  ENDDO
  
  DO ig=1,knon
     z1s = zdz2(nsoilmx)+dz1(nsoilmx-1)
     C_coef(ig,nsoilmx-1,indice)= &
          zdz2(nsoilmx)*ptsoil(ig,nsoilmx)/z1s
     D_coef(ig,nsoilmx-1,indice)=dz1(nsoilmx-1)/z1s
  ENDDO
  
  DO jk=nsoilmx-1,2,-1
     DO ig=1,knon
        z1s = 1./(zdz2(jk)+dz1(jk-1)+dz1(jk) &
             *(1.-D_coef(ig,jk,indice)))
        C_coef(ig,jk-1,indice)= &
             (ptsoil(ig,jk)*zdz2(jk)+dz1(jk)*C_coef(ig,jk,indice)) * z1s
        D_coef(ig,jk-1,indice)=dz1(jk-1)*z1s
     ENDDO
  ENDDO

!-----------------------------------------------------------------------
! 2)
! Computation of the soil temperatures using the Cgrd and Dgrd
! coefficient computed above
!
!-----------------------------------------------------------------------

!    Surface temperature
  DO ig=1,knon
     ptsoil(ig,1)=(lambda*C_coef(ig,1,indice)+ptsrf(ig))/  &
          (lambda*(1.-D_coef(ig,1,indice))+1.)
  ENDDO
  
!   Other temperatures
  DO jk=1,nsoilmx-1
     DO ig=1,knon
        ptsoil(ig,jk+1)=C_coef(ig,jk,indice)+D_coef(ig,jk,indice) &
             *ptsoil(ig,jk)
     ENDDO
  ENDDO

  IF (indice == is_sic) THEN
     DO ig = 1 , knon
        ptsoil(ig,nsoilmx) = RTT - 1.8
     END DO
  ENDIF

!-----------------------------------------------------------------------
! 3)
! Calculate the Cgrd and Dgrd coefficient corresponding to actual soil 
! temperature
!-----------------------------------------------------------------------
  DO ig=1,knon
     z1s = zdz2(nsoilmx)+dz1(nsoilmx-1)
     C_coef(ig,nsoilmx-1,indice) = zdz2(nsoilmx)*ptsoil(ig,nsoilmx)/z1s
     D_coef(ig,nsoilmx-1,indice) = dz1(nsoilmx-1)/z1s
  ENDDO
  
  DO jk=nsoilmx-1,2,-1
     DO ig=1,knon
        z1s = 1./(zdz2(jk)+dz1(jk-1)+dz1(jk) &
             *(1.-D_coef(ig,jk,indice)))
        C_coef(ig,jk-1,indice) = &
             (ptsoil(ig,jk)*zdz2(jk)+dz1(jk)*C_coef(ig,jk,indice)) * z1s
        D_coef(ig,jk-1,indice) = dz1(jk-1)*z1s
     ENDDO
  ENDDO

!-----------------------------------------------------------------------
! 4)
! Computation of the surface diffusive flux from ground and
! calorific capacity of the ground
!-----------------------------------------------------------------------
  DO ig=1,knon
     pfluxgrd(ig) = ztherm_i(ig)*dz1(1)* &
          (C_coef(ig,1,indice)+(D_coef(ig,1,indice)-1.)*ptsoil(ig,1))
     pcapcal(ig)  = ztherm_i(ig)* &
          (dz2(1)+ptimestep*(1.-D_coef(ig,1,indice))*dz1(1))
     z1s = lambda*(1.-D_coef(ig,1,indice))+1.
     pcapcal(ig)  = pcapcal(ig)/z1s
     pfluxgrd(ig) = pfluxgrd(ig) &
          + pcapcal(ig) * (ptsoil(ig,1) * z1s &
          - lambda * C_coef(ig,1,indice) &
          - ptsrf(ig)) &
          /ptimestep
  ENDDO
    
END SUBROUTINE soil
