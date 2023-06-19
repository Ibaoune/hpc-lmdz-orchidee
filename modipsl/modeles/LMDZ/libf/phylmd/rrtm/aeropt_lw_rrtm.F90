!
! aeropt_lw_rrtm.F90 2014-05-13 C. Kleinschmitt
!                    2016-05-03 O. Boucher 
!                    2016-12-17 O. Boucher
!
! This routine feeds aerosol LW properties to RRTM
! we only consider absorption (not scattering)
! we only consider dust for now

SUBROUTINE AEROPT_LW_RRTM(ok_alw, pdel, zrho, flag_aerosol, m_allaer, m_allaer_pi)

  USE dimphy
  USE aero_mod
  USE phys_state_var_mod, ONLY: tau_aero_lw_rrtm
  USE YOERAD, ONLY: NLW
  USE YOMCST, ONLY: RG

  IMPLICIT NONE

  INCLUDE "clesphys.h"
  !
  ! Input arguments:
  !
  LOGICAL, INTENT(IN)                              :: ok_alw
  INTEGER, INTENT(IN)                              :: flag_aerosol
  REAL, DIMENSION(klon,klev), INTENT(IN)           :: pdel, zrho
  REAL, DIMENSION(klon,klev,naero_tot), INTENT(IN) :: m_allaer, m_allaer_pi
  !
  INTEGER inu, i, k
  REAL :: zdh(klon,klev)
  REAL :: tmp_var, tmp_var_pi
  CHARACTER*20 modname
  !
  !--absorption coefficient for CIDUST
  REAL:: alpha_abs_CIDUST_16bands(nbands_lw_rrtm)   !--unit m2/g 
  DATA alpha_abs_CIDUST_16bands /                         &
  0.001, 0.003, 0.005, 0.006, 0.012, 0.030, 0.148, 0.098, &
  0.017, 0.053, 0.031, 0.008, 0.010, 0.011, 0.013, 0.015  /
  !
  modname='aeropt_lw_rrtm'
  !
  IF (NLW.NE.nbands_lw_rrtm) THEN
    CALL abort_physic(modname,'Erreur NLW doit etre egal a 16 pour cette routine',1)
  ENDIF
  ! 
  IF (ok_alw) THEN                                   !--aerosol LW effects
   !
   IF (flag_aerosol.EQ.5.OR.flag_aerosol.EQ.6.OR.flag_aerosol.EQ.7) THEN  !-Dust
    !
    zdh(:,:)=pdel(:,:)/(RG*zrho(:,:))      ! m
    !
    DO k=1, klev
      DO i=1, klon
         !
         tmp_var   =m_allaer(i,k,id_CIDUSTM_phy)   /1.e6*zdh(i,k)  !--g/m2
         tmp_var_pi=m_allaer_pi(i,k,id_CIDUSTM_phy)/1.e6*zdh(i,k)  !--g/m2
         !
         DO inu=1, NLW
           !
           !--total aerosol
           tau_aero_lw_rrtm(i,k,2,inu) = MAX(1.e-15,tmp_var*alpha_abs_CIDUST_16bands(inu))
           !--natural aerosol 
!           tau_aero_lw_rrtm(:,:,1,inu) = MAX(1.e-15,tmp_var_pi*alpha_abs_CIDUST_16bands(inu))
           tau_aero_lw_rrtm(i,k,1,inu) = 1.e-15  !--test
           !
         ENDDO
      ENDDO
      !
    ENDDO
    ! 
   ENDIF
   !
  ELSE !--no aerosol LW effects
    !
    tau_aero_lw_rrtm = 1.e-15
  ENDIF
  !
END SUBROUTINE AEROPT_LW_RRTM
