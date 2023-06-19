  SUBROUTINE reevap (klon,klev,iflag_ice_thermo,t_seri,q_seri,ql_seri,qs_seri, &
   &         d_t_eva,d_q_eva,d_ql_eva,d_qs_eva)

    ! flag to include modifications to ensure energy conservation (if flag >0)
    USE add_phys_tend_mod, only : fl_cor_ebil 
    
    IMPLICIT none
    !>======================================================================

    INTEGER klon,klev,iflag_ice_thermo
    REAL, DIMENSION(klon,klev), INTENT(in) :: t_seri,q_seri,ql_seri,qs_seri
    REAL, DIMENSION(klon,klev), INTENT(out) :: d_t_eva,d_q_eva,d_ql_eva,d_qs_eva

    REAL za,zb,zdelta,zlvdcp,zlsdcp
    INTEGER i,k

    !--------Stochastic Boundary Layer Triggering: ALE_BL--------
    !---Propri\'et\'es du thermiques au LCL 
    include "YOMCST.h"
    include "YOETHF.h"
    include "FCTTRE.h"
    !IM 100106 BEG : pouvoir sortir les ctes de la physique
    !
    ! Re-evaporer l'eau liquide nuageuse
    !
!print *,'rrevap ; fl_cor_ebil:',fl_cor_ebil,' iflag_ice_thermo:',iflag_ice_thermo,' RVTMP2',RVTMP2
    DO k = 1, klev  ! re-evaporation de l'eau liquide nuageuse
       DO i = 1, klon
        if (fl_cor_ebil .GT. 0) then
          zlvdcp=RLVTT/RCPD/(1.0+RVTMP2*(q_seri(i,k)+ql_seri(i,k)+qs_seri(i,k)))
          zlsdcp=RLSTT/RCPD/(1.0+RVTMP2*(q_seri(i,k)+ql_seri(i,k)+qs_seri(i,k)))
        else
          zlvdcp=RLVTT/RCPD/(1.0+RVTMP2*q_seri(i,k))
          !jyg<
          !  Attention : Arnaud a propose des formules completement differentes
          !                  A verifier !!!
          zlsdcp=RLSTT/RCPD/(1.0+RVTMP2*q_seri(i,k))
        end if
          IF (iflag_ice_thermo .EQ. 0) THEN
             zlsdcp=zlvdcp
          ENDIF
          !>jyg

          IF (iflag_ice_thermo.eq.0) THEN   
             !pas necessaire a priori

             zdelta = MAX(0.,SIGN(1.,RTT-t_seri(i,k)))
  zdelta = 0.
             zb = MAX(0.0,ql_seri(i,k))
             za = - MAX(0.0,ql_seri(i,k)) &
                  * (zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
             d_t_eva(i,k) = za
             d_q_eva(i,k) = zb
             d_ql_eva(i,k) = -ql_seri(i,k)
             d_qs_eva(i,k) = 0.

          ELSE

             !CR: on r\'e-\'evapore eau liquide et glace

             !        zdelta = MAX(0.,SIGN(1.,RTT-t_seri(i,k)))
             !        zb = MAX(0.0,ql_seri(i,k))
             !        za = - MAX(0.0,ql_seri(i,k)) &
             !             * (zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
             zb = MAX(0.0,ql_seri(i,k)+qs_seri(i,k))
             za = - MAX(0.0,ql_seri(i,k))*zlvdcp & 
                  - MAX(0.0,qs_seri(i,k))*zlsdcp
             d_t_eva(i,k) = za
             d_q_eva(i,k) = zb
             d_ql_eva(i,k) = -ql_seri(i,k)
             d_qs_eva(i,k) = -qs_seri(i,k)
          ENDIF

       ENDDO
    ENDDO

RETURN

END SUBROUTINE reevap
