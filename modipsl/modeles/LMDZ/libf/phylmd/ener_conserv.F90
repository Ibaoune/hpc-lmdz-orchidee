subroutine ener_conserv(klon,klev,pdtphys, &
 &                      puo,pvo,pto,pqo,pql0,pqs0, &
 &                      pun,pvn,ptn,pqn,pqln,pqsn,dtke,masse,exner,d_t_ec)

!=============================================================
! Energy conservation
! Based on the TKE equation
! The M2 and N2 terms at the origin of TKE production are
! concerted into heating in the d_t_ec term
! Option 1 is the standard
!        101 is for M2 term only
!        101 for N2 term only
!         -1 is a previours treatment for kinetic energy only
!  FH (hourdin@lmd.jussieu.fr), 2013/04/25
!=============================================================

!=============================================================
! Declarations
!=============================================================

! From module
USE phys_local_var_mod, ONLY : d_u_vdf,d_v_vdf,d_t_vdf,d_u_ajs,d_v_ajs,d_t_ajs, &
 &                             d_u_con,d_v_con,d_t_con,d_t_diss
USE phys_local_var_mod, ONLY : d_t_eva,d_t_lsc,d_q_eva,d_q_lsc
USE phys_local_var_mod, ONLY : d_u_oro,d_v_oro,d_u_lif,d_v_lif
USE phys_local_var_mod, ONLY : du_gwd_hines,dv_gwd_hines,dv_gwd_front,dv_gwd_rando
USE phys_state_var_mod, ONLY : du_gwd_front,du_gwd_rando
USE phys_output_var_mod, ONLY : bils_ec,bils_ech,bils_tke,bils_kinetic,bils_enthalp,bils_latent,bils_diss
USE add_phys_tend_mod, ONLY : fl_cor_ebil 


IMPLICIT none
#include "YOMCST.h"
#include "YOETHF.h"
#include "clesphys.h"
#include "compbl.h"

! Arguments
INTEGER, INTENT(IN) :: klon,klev
REAL, INTENT(IN) :: pdtphys
REAL, DIMENSION(klon,klev), INTENT(IN)      :: puo,pvo,pto,pqo,pql0,pqs0
REAL, DIMENSION(klon,klev), INTENT(IN)      :: pun,pvn,ptn,pqn,pqln,pqsn
REAL, DIMENSION(klon,klev), INTENT(IN)      :: masse,exner
REAL, DIMENSION(klon,klev+1), INTENT(IN)    :: dtke
!
REAL, DIMENSION(klon,klev), INTENT(OUT)     :: d_t_ec

! Local
      integer k,i
REAL, DIMENSION(klon,klev+1) :: fluxu,fluxv,fluxt
REAL, DIMENSION(klon,klev+1) :: dddu,dddv,dddt
REAL, DIMENSION(klon,klev) :: d_u,d_v,d_t,zv,zu,d_t_ech
REAL ZRCPD

character*80 abort_message
character*20 :: modname


modname='ener_conser'
d_t_ec(:,:)=0.

IF (iflag_ener_conserv==-1) THEN
!+jld ec_conser
   DO k = 1, klev
   DO i = 1, klon
     IF (fl_cor_ebil .GT. 0) then
       ZRCPD = RCPD*(1.0+RVTMP2*(pqn(i,k)+pqln(i,k)+pqsn(i,k)))
     ELSE
       ZRCPD = RCPD*(1.0+RVTMP2*pqn(i,k))
     ENDIF
     d_t_ec(i,k)=0.5/ZRCPD &
 &     *(puo(i,k)**2+pvo(i,k)**2-pun(i,k)**2-pvn(i,k)**2)
   ENDDO
   ENDDO
!-jld ec_conser



ELSEIF (iflag_ener_conserv>=1) THEN

   IF (iflag_ener_conserv<=2) THEN
!     print*,'ener_conserv pbl=',iflag_pbl
      IF (iflag_pbl>=20 .AND. iflag_pbl<=27) THEN !d_t_diss accounts for conserv
         d_t(:,:)=d_t_ajs(:,:)   ! d_t_ajs = adjust + thermals
         d_u(:,:)=d_u_ajs(:,:)+d_u_con(:,:)
         d_v(:,:)=d_v_ajs(:,:)+d_v_con(:,:)
      ELSE
         d_t(:,:)=d_t_vdf(:,:)+d_t_ajs(:,:)   ! d_t_ajs = adjust + thermals
         d_u(:,:)=d_u_vdf(:,:)+d_u_ajs(:,:)+d_u_con(:,:)
         d_v(:,:)=d_v_vdf(:,:)+d_v_ajs(:,:)+d_v_con(:,:)
      ENDIF
   ELSEIF (iflag_ener_conserv==101) THEN
      d_t(:,:)=0.
      d_u(:,:)=d_u_vdf(:,:)+d_u_ajs(:,:)+d_u_con(:,:)
      d_v(:,:)=d_v_vdf(:,:)+d_v_ajs(:,:)+d_v_con(:,:)
   ELSEIF (iflag_ener_conserv==110) THEN
      d_t(:,:)=d_t_vdf(:,:)+d_t_ajs(:,:)
      d_u(:,:)=0.
      d_v(:,:)=0.

   ELSEIF (iflag_ener_conserv==3) THEN
      d_t(:,:)=0.
      d_u(:,:)=0.
      d_v(:,:)=0.
   ELSEIF (iflag_ener_conserv==4) THEN
      d_t(:,:)=0.
      d_u(:,:)=d_u_vdf(:,:)
      d_v(:,:)=d_v_vdf(:,:)
   ELSEIF (iflag_ener_conserv==5) THEN
      d_t(:,:)=d_t_vdf(:,:)
      d_u(:,:)=d_u_vdf(:,:)
      d_v(:,:)=d_v_vdf(:,:)
   ELSEIF (iflag_ener_conserv==6) THEN
      d_t(:,:)=d_t_vdf(:,:)
      d_u(:,:)=d_u_vdf(:,:)+d_u_ajs(:,:)
      d_v(:,:)=d_v_vdf(:,:)+d_v_ajs(:,:)
   ELSEIF (iflag_ener_conserv==7) THEN
      d_t(:,:)=d_t_vdf(:,:)+d_t_ajs(:,:)
      d_u(:,:)=d_u_vdf(:,:)+d_u_ajs(:,:)
      d_v(:,:)=d_v_vdf(:,:)+d_v_ajs(:,:)
   ELSEIF (iflag_ener_conserv==8) THEN
      d_t(:,:)=d_t_vdf(:,:)
      d_u(:,:)=d_u_vdf(:,:)+d_u_ajs(:,:)+d_u_con(:,:)
      d_v(:,:)=d_v_vdf(:,:)+d_v_ajs(:,:)+d_v_con(:,:)
   ELSEIF (iflag_ener_conserv==9) THEN
      d_t(:,:)=d_t_vdf(:,:)
      d_u(:,:)=d_u_vdf(:,:)+d_u_ajs(:,:)+d_u_con(:,:)+d_u_oro(:,:)
      d_v(:,:)=d_v_vdf(:,:)+d_v_ajs(:,:)+d_v_con(:,:)+d_v_oro(:,:)
   ELSEIF (iflag_ener_conserv==10) THEN
      d_t(:,:)=d_t_vdf(:,:)
      d_u(:,:)=d_u_vdf(:,:)+d_u_ajs(:,:)+d_u_con(:,:)+d_u_oro(:,:)+d_u_lif(:,:)
      d_v(:,:)=d_v_vdf(:,:)+d_v_ajs(:,:)+d_v_con(:,:)+d_v_oro(:,:)+d_v_lif(:,:)
   ELSEIF (iflag_ener_conserv==11) THEN
      d_t(:,:)=d_t_vdf(:,:)
      d_u(:,:)=d_u_vdf(:,:)+d_u_ajs(:,:)+d_u_con(:,:)+d_u_oro(:,:)+d_u_lif(:,:)
      d_v(:,:)=d_v_vdf(:,:)+d_v_ajs(:,:)+d_v_con(:,:)+d_v_oro(:,:)+d_v_lif(:,:) 
      IF (ok_hines) THEN
         d_u_vdf(:,:)=d_u_vdf(:,:)+du_gwd_hines(:,:)
         d_v_vdf(:,:)=d_v_vdf(:,:)+dv_gwd_hines(:,:)
      ENDIF
      IF (.not. ok_hines .and. ok_gwd_rando) THEN
         d_u_vdf(:,:)=d_u_vdf(:,:)+du_gwd_front(:,:)
         d_v_vdf(:,:)=d_v_vdf(:,:)+dv_gwd_front(:,:)
      ENDIF
      IF (ok_gwd_rando) THEN
         d_u_vdf(:,:)=d_u_vdf(:,:)+du_gwd_rando(:,:)
         d_v_vdf(:,:)=d_v_vdf(:,:)+dv_gwd_rando(:,:)
      ENDIF
   ELSE
      abort_message = 'iflag_ener_conserv non prevu'
      CALL abort_physic (modname,abort_message,1)
   ENDIF

!----------------------------------------------------------------------------
! Two options wether we consider time integration in the energy conservation
!----------------------------------------------------------------------------

   if (iflag_ener_conserv==2) then
      zu(:,:)=puo(:,:)
      zv(:,:)=pvo(:,:)
   else
      IF (iflag_pbl>=20 .AND. iflag_pbl<=27) THEN
         zu(:,:)=puo(:,:)+d_u_vdf(:,:)+0.5*d_u(:,:)
         zv(:,:)=pvo(:,:)+d_v_vdf(:,:)+0.5*d_v(:,:)
      ELSE
         zu(:,:)=puo(:,:)+0.5*d_u(:,:)
         zv(:,:)=pvo(:,:)+0.5*d_v(:,:)
      ENDIF
   endif

   fluxu(:,klev+1)=0.
   fluxv(:,klev+1)=0.
   fluxt(:,klev+1)=0.

   do k=klev,1,-1
      fluxu(:,k)=fluxu(:,k+1)+masse(:,k)*d_u(:,k)
      fluxv(:,k)=fluxv(:,k+1)+masse(:,k)*d_v(:,k)
      fluxt(:,k)=fluxt(:,k+1)+masse(:,k)*d_t(:,k)/exner(:,k)
   enddo

   dddu(:,1)=2*zu(:,1)*fluxu(:,1)
   dddv(:,1)=2*zv(:,1)*fluxv(:,1)
   dddt(:,1)=(exner(:,1)-1.)*fluxt(:,1)

   do k=2,klev
      dddu(:,k)=(zu(:,k)-zu(:,k-1))*fluxu(:,k)
      dddv(:,k)=(zv(:,k)-zv(:,k-1))*fluxv(:,k)
      dddt(:,k)=(exner(:,k)-exner(:,k-1))*fluxt(:,k)
   enddo
   dddu(:,klev+1)=0.
   dddv(:,klev+1)=0.
   dddt(:,klev+1)=0.

   do k=1,klev
      d_t_ech(:,k)=-(rcpd*(dddt(:,k)+dddt(:,k+1)))/(2.*rcpd*masse(:,k))
      d_t_ec(:,k)=-(dddu(:,k)+dddu(:,k+1)+dddv(:,k)+dddv(:,k+1))/(2.*rcpd*masse(:,k))+d_t_ech(:,k)
   enddo

ENDIF

!================================================================
!  Computation of integrated enthalpie and kinetic energy variation
!  FH (hourdin@lmd.jussieu.fr), 2013/04/25
!  bils_ec : energie conservation term
!  bils_ech : part of this term linked to temperature
!  bils_tke : change of TKE
!  bils_diss : dissipation of TKE (when activated)
!  bils_kinetic : change of kinetic energie of the column
!  bils_enthalp : change of enthalpie
!  bils_latent  : change of latent heat. Computed between
!          after reevaporation (at the beginning of the physics)
!          and before large scale condensation (fisrtilp)
!================================================================

      bils_ec(:)=0.
      bils_ech(:)=0.
      bils_tke(:)=0.
      bils_diss(:)=0.
      bils_kinetic(:)=0.
      bils_enthalp(:)=0.
      bils_latent(:)=0.
      DO k=1,klev
        bils_ec(:)=bils_ec(:)-d_t_ec(:,k)*masse(:,k)
        bils_diss(:)=bils_diss(:)-d_t_diss(:,k)*masse(:,k)
        bils_kinetic(:)=bils_kinetic(:)+masse(:,k)* &
     &           (pun(:,k)*pun(:,k)+pvn(:,k)*pvn(:,k) &
     &            -puo(:,k)*puo(:,k)-pvo(:,k)*pvo(:,k))
        bils_enthalp(:)= &
     &  bils_enthalp(:)+masse(:,k)*(ptn(:,k)-pto(:,k)+d_t_ec(:,k)-d_t_eva(:,k)-d_t_lsc(:,k))
!    &  bils_enthalp(:)+masse(:,k)*(ptn(:,k)-pto(:,k)+d_t_ec(:,k))
        bils_latent(:)=bils_latent(:)+masse(:,k)* &
!    &             (pqn(:,k)-pqo(:,k))
     &             (pqn(:,k)-pqo(:,k)-d_q_eva(:,k)-d_q_lsc(:,k))
      ENDDO
      bils_ec(:)=rcpd*bils_ec(:)/pdtphys
      bils_diss(:)=rcpd*bils_diss(:)/pdtphys
      bils_kinetic(:)= 0.5*bils_kinetic(:)/pdtphys
      bils_enthalp(:)=rcpd*bils_enthalp(:)/pdtphys
      bils_latent(:)=rlvtt*bils_latent(:)/pdtphys
!jyg<
      IF (iflag_pbl > 1) THEN
        DO k=1,klev
          bils_tke(:)=bils_tke(:)+0.5*(dtke(:,k)+dtke(:,k+1))*masse(:,k)
        ENDDO
        bils_tke(:)=bils_tke(:)/pdtphys
      ENDIF  ! (iflag_pbl > 1)
!>jyg

IF (iflag_ener_conserv>=1) THEN
      bils_ech(:)=0.
      DO k=1,klev
        bils_ech(:)=bils_ech(:)-d_t_ech(:,k)*masse(:,k)
      ENDDO
      bils_ech(:)=rcpd*bils_ech(:)/pdtphys
ENDIF

RETURN

END
