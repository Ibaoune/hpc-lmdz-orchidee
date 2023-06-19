!
! splaeropt_lw_rrtm.F90 2014-05-13 C. Kleinschmitt
!                      2016-05-03 O. Boucher 
!
! This routine feeds aerosol LW properties to RRTM
! we only consider absorption (not scattering)

SUBROUTINE SPLAEROPT_LW_RRTM(ok_alw,zdm,tr_seri)

  USE dimphy
  USE aero_mod
  USE infotrac_phy
  USE phys_state_var_mod, ONLY : tau_aero_lw_rrtm
  USE YOERAD, ONLY : NLW

  IMPLICIT NONE

  INCLUDE "clesphys.h"
  !
  ! Input arguments:
  !
  LOGICAL, INTENT(IN) :: ok_alw
  REAL, DIMENSION(klon,klev), INTENT(IN)      :: zdm
  REAL, DIMENSION(klon,klev,nbtr), INTENT(IN) :: tr_seri
  !
  ! Local arguments : 
  !
  INTEGER, PARAMETER :: naero_soluble=2    ! 1- accumulation soluble; 2- coarse soluble
  INTEGER, PARAMETER :: naero_insoluble=2  ! 1- coarse dust; 2- supercoarse dust
  INTEGER, PARAMETER :: naero=naero_soluble+naero_insoluble
  !
  INTEGER inu, itr, spinsol
  CHARACTER*20 modname
  !
  !--absorption coefficient for coarse and super-coarse DUST
  REAL:: alpha_abs_CIDUST_16bands(nbands_lw_rrtm,naero_insoluble)   !--unit m2/g 
  DATA alpha_abs_CIDUST_16bands /                         &
   ! Dust CO insoluble
  0.001, 0.003, 0.005, 0.006, 0.011, 0.031, 0.157, 0.102, &
  0.017, 0.056, 0.032, 0.008, 0.010, 0.011, 0.013, 0.016, & 
   ! Dust SC insoluble
  0.002, 0.004, 0.007, 0.010, 0.018, 0.043, 0.099, 0.071, &
  0.021, 0.056, 0.033, 0.011, 0.013, 0.014, 0.016, 0.018 /

  modname='splaeropt_lw_rrtm'
  ! 
  IF (NLW.NE.nbands_lw_rrtm) THEN
    CALL abort_physic(modname,'Erreur NLW doit etre egal a 16 pour cette routine',1)
  ENDIF
  !
  IF (ok_alw) THEN 
    !
    !--initialisation
    tau_aero_lw_rrtm = 0.0
    !
    DO itr=1,nbtr
      !
      IF (tname(itr+nqo)=='PREC') THEN       !--precursor
        CYCLE
      ELSE IF (tname(itr+nqo)=='FINE') THEN  !--fine mode accumulation mode
        CYCLE
      ELSE IF (tname(itr+nqo)=='COSS') THEN  !--coarse mode sea salt
        CYCLE
      ELSE IF (tname(itr+nqo)=='CODU') THEN  !--coarse mode dust
        spinsol=1
      ELSE IF (tname(itr+nqo)=='SCDU') THEN  !--super coarse mode dust
        spinsol=2
      ELSE
         CALL abort_physic(modname,'I cannot do aerosol optics for '//tname(itr+nqo),1)
      ENDIF
      !
      DO inu=1,NLW
        !
        !--total aerosol
        tau_aero_lw_rrtm(:,:,2,inu) = tau_aero_lw_rrtm(:,:,2,inu) + tr_seri(:,:,itr)*zdm(:,:)*alpha_abs_CIDUST_16bands(inu,spinsol)
        !--no aerosol at all
        tau_aero_lw_rrtm(:,:,1,inu) = tau_aero_lw_rrtm(:,:,1,inu) + 0.0
        !
      ENDDO
    !
    ENDDO
    !
    !--avoid very small values
    tau_aero_lw_rrtm = MAX(tau_aero_lw_rrtm,1.e-15)
    ! 
  ELSE
    !--default value
    tau_aero_lw_rrtm = 1.e-15
  ENDIF
  !
END SUBROUTINE SPLAEROPT_LW_RRTM
