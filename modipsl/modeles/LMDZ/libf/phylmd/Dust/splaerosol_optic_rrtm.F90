! $Id: splaerosol_optic_rrtm.F90 2644 2016-10-02 16:55:08Z oboucher $
!
SUBROUTINE splaerosol_optic_rrtm( ok_alw, pplay, paprs, t_seri, rhcl, &
     tr_seri, mass_solu_aero, mass_solu_aero_pi, &
     tau_aero, piz_aero, cg_aero, &
     tausum_aero, tau3d_aero )

  ! This routine will :
  ! 1) recevie the aerosols(already read and interpolated) corresponding to flag_aerosol
  ! 2) calculate the optical properties for the aerosols
  !

  USE dimphy
  USE aero_mod
  USE infotrac_phy
  USE YOMCST, ONLY: RD, RG

  IMPLICIT NONE

  INCLUDE "clesphys.h"


  ! Input arguments
  !****************************************************************************************
  LOGICAL, INTENT(IN) :: ok_alw                      ! Apply aerosol LW effect or not
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: pplay
  REAL, DIMENSION(klon,klev+1), INTENT(IN) :: paprs
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: t_seri
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: rhcl   ! humidite relative ciel clair
  REAL, DIMENSION(klon,klev,nbtr), INTENT(IN) :: tr_seri ! concentration tracer

  ! Output arguments
  !****************************************************************************************
  REAL, DIMENSION(klon,klev), INTENT(OUT)     :: mass_solu_aero    ! Total mass for all soluble aerosols
  REAL, DIMENSION(klon,klev), INTENT(OUT)     :: mass_solu_aero_pi !     -"-     preindustrial values
  REAL, DIMENSION(klon,klev,2,NSW), INTENT(OUT) :: tau_aero    ! Aerosol optical thickness
  REAL, DIMENSION(klon,klev,2,NSW), INTENT(OUT) :: piz_aero    ! Single scattering albedo aerosol
  REAL, DIMENSION(klon,klev,2,NSW), INTENT(OUT) :: cg_aero     ! asymmetry parameter aerosol
  REAL, DIMENSION(klon,nwave,naero_tot), INTENT(OUT)       :: tausum_aero
  REAL, DIMENSION(klon,klev,nwave,naero_tot), INTENT(OUT)  :: tau3d_aero

  INTEGER i, k, itr
  REAL, DIMENSION(klon,klev) :: zdm, zdh
  REAL zrho, pdel
  !
  ! Calculate the total mass of all soluble accumulation mode aerosols
  ! to be revisited for AR6
  !
  mass_solu_aero(:,:)    = 0.0
  mass_solu_aero_pi(:,:) = 0.0
  !
  DO itr=1,nbtr
    IF (tname(itr+nqo)=='FINE') THEN
      mass_solu_aero(:,:)    = tr_seri(:,:,itr)
      mass_solu_aero_pi(:,:) = tr_seri(:,:,itr)
    ENDIF
  ENDDO

  ! Calculate layer thickness
  DO k = 1, klev
    DO i = 1, klon
      pdel=paprs(i,k)-paprs(i,k+1)
      zrho=pplay(i,k)/t_seri(i,k)/RD             ! kg/m3
      zdh(i,k)=pdel/(RG*zrho)                    ! m
      zdm(i,k)=pdel/RG                           ! kg/m2
    ENDDO
  ENDDO

!--new aerosol properties
! aeropt_6bands for rrtm
  CALL splaeropt_6bands_rrtm(       &
       zdm, tr_seri, rhcl,          & 
       tau_aero, piz_aero, cg_aero )

! aeropt_5wv only for validation and diagnostics
  CALL splaeropt_5wv_rrtm(          &
       zdm, zdh, tr_seri, rhcl,     & 
       tausum_aero, tau3d_aero )

! LW optical properties for tropospheric aerosols
  CALL splaeropt_lw_rrtm(ok_alw,zdm,tr_seri)

END SUBROUTINE splaerosol_optic_rrtm
