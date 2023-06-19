C Subroutine that calculates the emission of aerosols precursors
      SUBROUTINE precuremission(ftsol,u10m_ec,v10m_ec,
     .                          pctsrf,u_seri,v_seri,paprs,pplay,cdragh,
     .                          cdragm,t_seri,q_seri,tsol,fracso2emis,
     .                          frach2sofso2,bateau,zdz,zalt,
     .                          kminbc,kmaxbc,pdtphys,scale_param_bb,
     .                          scale_param_ind,iregion_ind,iregion_bb,
     .                          nbreg_ind, nbreg_bb,
     .                          lmt_so2ff_l,lmt_so2ff_h,lmt_so2nff,
     .                          lmt_so2ba,lmt_so2bb_l,lmt_so2bb_h,
     .                          lmt_so2volc_cont,lmt_altvolc_cont,
     .                          lmt_so2volc_expl,lmt_altvolc_expl,
     .                          lmt_dmsbio,lmt_h2sbio, lmt_dmsconc, 
     .                          lmt_dms,id_prec,id_fine,
     .                                 flux_sparam_ind,flux_sparam_bb,
     .                                 source_tr,flux_tr,tr_seri)

      USE dimphy
      USE indice_sol_mod
      USE infotrac
!      USE phytracr_spl_mod, ONLY : nbreg_dust, nbreg_ind, nbreg_bb
      IMPLICIT NONE

#include "dimensions.h"
#include "chem.h"
#include "chem_spla.h"
#include "YOMCST.h"
#include "paramet.h"

c============================= INPUT ===================================
      INTEGER kminbc, kmaxbc
      REAL ftsol(klon,nbsrf)  ! temperature du sol par type
      REAL tsol(klon)         ! temperature du sol moyenne
      REAL t_seri(klon,klev)  ! temperature
      REAL u_seri(klon,klev)  ! vent
      REAL v_seri(klon,klev)  ! vent
      REAL q_seri(klon,klev)  ! vapeur d eau kg/kg
      REAL u10m_ec(klon), v10m_ec(klon)  ! vent a 10 metres
      REAL pctsrf(klon,nbsrf)
      REAL pdtphys  ! pas d'integration pour la physique (seconde)
      REAL paprs(klon,klev+1)  ! pression pour chaque inter-couche (en Pa)
      REAL pplay(klon,klev)  ! pression pour le mileu de chaque couche (en Pa)
      REAL cdragh(klon), cdragm(klon)      
      REAL fracso2emis        !--fraction so2 emis en so2
      REAL frach2sofso2       !--fraction h2s from so2
      REAL zdz(klon,klev)
      LOGICAL edgar, bateau
      INTEGER id_prec,id_fine
c
c------------------------- Scaling Parameters --------------------------
c
      INTEGER nbreg_ind, nbreg_bb
      INTEGER iregion_ind(klon)  !Defines regions for SO2, BC & OM
      INTEGER iregion_bb(klon)  !Defines regions for SO2, BC & OM
      REAL scale_param_bb(nbreg_bb)  !Scaling parameter for biomas burning
      REAL scale_param_ind(nbreg_ind) !Scaling parameter for industrial emissions
c
c============================= OUTPUT ==================================
c
      REAL source_tr(klon,nbtr)
      REAL flux_tr(klon,nbtr)
      REAL tr_seri(klon,klev,nbtr) ! traceur
      REAL flux_sparam_ind(klon), flux_sparam_bb(klon)
c========================= LOCAL VARIABLES =============================
      INTEGER i, k, kkk_cont(klon), kkk_expl(klon)
      REAL zalt(klon,klev), zaltmid(klon,klev)
      REAL zzdz
c------------------------- SULFUR emissions ----------------------------
      REAL lmt_so2volc_cont(klon)  ! emissions so2 volcan (continuous)
      REAL lmt_altvolc_cont(klon)  ! altitude  so2 volcan (continuous)
      REAL lmt_so2volc_expl(klon)  ! emissions so2 volcan (explosive)
      REAL lmt_altvolc_expl(klon)  ! altitude  so2 volcan (explosive)
      REAL lmt_so2ff_l(klon)       ! emissions so2 fossil fuel (low)
      REAL lmt_so2ff_h(klon)       ! emissions so2 fossil fuel (high)
      REAL lmt_so2nff(klon)        ! emissions so2 non-fossil fuel
      REAL lmt_so2bb_l(klon)       ! emissions de so2 biomass burning (low)
      REAL lmt_so2bb_h(klon)       ! emissions de so2 biomass burning (high)
      REAL lmt_so2ba(klon)         ! emissions de so2 bateau
      REAL lmt_dms(klon)           ! emissions de dms
      REAL lmt_dmsconc(klon)       ! concentration de dms oceanique
      REAL lmt_dmsbio(klon)        ! emissions de dms bio
      REAL lmt_h2sbio(klon)        ! emissions de h2s bio
                        
      EXTERNAL condsurfs, liss, nightingale
c=========================================================================
c Modifications introduced by NHL
c -Variables to save fluxes were introduced
c -lmt_so2ba was multiplied by fracso2emis in line 117
c -scale_param_bb was introduced in line 105
c The last two modifications were errors existing in the original version 
c=========================================================================
c=========================================================================
c                        LOW LEVEL EMISSIONS
c=========================================================================
                        
         CALL nightingale(u_seri, v_seri, u10m_ec, v10m_ec, paprs,
     .                 pplay, cdragh, cdragm, t_seri, q_seri, ftsol,
     .                 tsol, pctsrf, lmt_dmsconc, lmt_dms)

      IF (.not.bateau) THEN
        DO i=1, klon      
          lmt_so2ba(i)=0.0
        ENDDO
      ENDIF

      DO i=1, klon
         IF (iregion_ind(i).GT.0) THEN
       IF(id_prec>0) source_tr(i,id_prec)=source_tr(i,id_prec) 
     .             + fracso2emis
     .             *scale_param_ind(iregion_ind(i))*lmt_so2ff_l(i)*1.e4
     .             +scale_param_ind(iregion_ind(i))*lmt_so2ff_l(i)*1.e4
     .                   *frach2sofso2            ! molec/m2/s
c
      IF(id_fine>0) source_tr(i,id_fine)=
     .                source_tr(i,id_fine)+(1-fracso2emis)
     .                *scale_param_ind(iregion_ind(i))*lmt_so2ff_l(i)
     .                *1.e4*masse_ammsulfate/RNAVO  ! g/m2/s
c
       IF(id_prec>0)   flux_tr(i,id_prec)=flux_tr(i,id_prec) + (
     .               scale_param_ind(iregion_ind(i))*(lmt_so2ff_l(i)+
     .                lmt_so2ff_h(i))
     .                *frach2sofso2
     .               +scale_param_ind(iregion_ind(i))*(lmt_so2ff_l(i)+
     .                lmt_so2ff_h(i))
     .                *fracso2emis
     .                )*1.e4/RNAVO*masse_s*1.e3          ! mgS/m2/s
c
      IF(id_fine>0)  flux_tr(i,id_fine)=
     . flux_tr(i,id_fine)+(1-fracso2emis)
     .               *scale_param_ind(iregion_ind(i))*(lmt_so2ff_l(i)+
     .                lmt_so2ff_h(i))
     .               *1.e4/RNAVO*masse_ammsulfate*1.e3    ! mgS/m2/s
c
      flux_sparam_ind(i)=flux_sparam_ind(i)+ (1-fracso2emis)
     .               *scale_param_ind(iregion_ind(i))*(lmt_so2ff_l(i)+
     .                lmt_so2ff_h(i))
     .               *1.e4/RNAVO*masse_ammsulfate*1.e3    ! mgS/m2/s
         ENDIF
         IF (iregion_bb(i).GT.0) THEN
      IF(id_prec>0) source_tr(i,id_prec)= 
     .                  source_tr(i,id_prec) + fracso2emis
     .                 *scale_param_bb(iregion_bb(i))*lmt_so2bb_l(i)
     .                   *(1.-pctsrf(i,is_oce))*1.e4
c
      IF(id_fine>0)     source_tr(i,id_fine)=
     .                   source_tr(i,id_fine)+(1-fracso2emis)
     .                  *scale_param_bb(iregion_bb(i))*lmt_so2bb_l(i)*
     .                   (1.-pctsrf(i,is_oce))*1.e4*
     .                   masse_ammsulfate/RNAVO  ! g/m2/s
c
      IF(id_prec>0)     flux_tr(i,id_prec)=flux_tr(i,id_prec) +
     .               (scale_param_bb(iregion_bb(i))*lmt_so2bb_l(i)
     .                 +scale_param_bb(iregion_bb(i))*lmt_so2bb_h(i))
     .                      * (1.-pctsrf(i,is_oce))*fracso2emis
     .                 *1.e4/RNAVO*masse_s*1.e3          ! mgS/m2/s
c
      IF(id_fine>0) flux_tr(i,id_fine)=
     .                flux_tr(i,id_fine)+(1-fracso2emis)
     .               *(scale_param_bb(iregion_bb(i))*lmt_so2bb_l(i)
     .                +scale_param_bb(iregion_bb(i))*lmt_so2bb_h(i))
     .                *(1.-pctsrf(i,is_oce))
     .                *1.e4/RNAVO*masse_ammsulfate*1.e3    ! mgS/m2/s
c
           flux_sparam_bb(i)=
     .                scale_param_bb(iregion_bb(i))*(lmt_so2bb_l(i)+
     .                                        lmt_so2bb_h(i))
     .                      * (1.-pctsrf(i,is_oce))*fracso2emis
     .              *1.e4/RNAVO*masse_s*1.e3          ! mgS/m2/s
           flux_sparam_bb(i)= flux_sparam_bb(i) + (1-fracso2emis) *
     .               (scale_param_bb(iregion_bb(i))*lmt_so2bb_l(i)+
     .                scale_param_bb(iregion_bb(i))*lmt_so2bb_h(i))
     .                *(1.-pctsrf(i,is_oce))
     .               *1.e4/RNAVO*masse_ammsulfate*1.e3    ! mgS/m2/s
         ENDIF
      IF(id_prec>0)   source_tr(i,id_prec)=source_tr(i,id_prec) 
     .                 + fracso2emis
     .                 *(lmt_so2ba(i)+lmt_so2nff(i))*1.e4
     .                 +(lmt_h2sbio(i)
     .                   +lmt_dms(i)+lmt_dmsbio(i))*1.e4            ! molec/m2/s
c
      IF(id_fine>0)   source_tr(i,id_fine)=source_tr(i,id_fine)
     .                +(1-fracso2emis)
     .                *(lmt_so2ba(i)+lmt_so2nff(i))*1.e4*
     .                   masse_ammsulfate/RNAVO  ! g/m2/s
c
      IF(id_prec>0)   flux_tr(i,id_prec)=flux_tr(i,id_prec) 
     .               + (lmt_h2sbio(i)
     .               +lmt_so2volc_cont(i)+lmt_so2volc_expl(i)
     .                 +(lmt_so2ba(i)+lmt_so2nff(i))*fracso2emis
     .                 +lmt_dms(i)+lmt_dmsbio(i) )
     .              *1.e4/RNAVO*masse_s*1.e3          ! mgS/m2/s
c
      IF(id_fine>0)   flux_tr(i,id_fine)=flux_tr(i,id_fine)
     .               +(1-fracso2emis)
     .               *(lmt_so2ba(i) + lmt_so2nff(i))
     .               *1.e4/RNAVO*masse_ammsulfate*1.e3    ! mgS/m2/s
c
         flux_sparam_ind(i)=flux_sparam_ind(i)+ (1-fracso2emis)
     .               *lmt_so2nff(i)
     .               *1.e4/RNAVO*masse_ammsulfate*1.e3    ! mgS/m2/s
c
      ENDDO

c========================================================================
c                        HIGH LEVEL EMISSIONS
c========================================================================
c  Source de SO2 volcaniques
      DO i = 1, klon
        kkk_cont(i)=1
        kkk_expl(i)=1
      ENDDO
      DO k=1, klev-1
      DO i = 1, klon
        zaltmid(i,k)=zalt(i,k)+zdz(i,k)/2.
        IF (zalt(i,k+1).LT.lmt_altvolc_cont(i)) kkk_cont(i)=k+1
        IF (zalt(i,k+1).LT.lmt_altvolc_expl(i)) kkk_expl(i)=k+1
      ENDDO
      ENDDO
      IF(id_prec>0) THEN
      DO i = 1, klon
        tr_seri(i,kkk_cont(i),id_prec)=tr_seri(i,kkk_cont(i),id_prec) +
     .               lmt_so2volc_cont(i)/zdz(i,kkk_cont(i))/100.*pdtphys
        tr_seri(i,kkk_expl(i),id_prec)=tr_seri(i,kkk_expl(i),id_prec) +
     .               lmt_so2volc_expl(i)/zdz(i,kkk_expl(i))/100.*pdtphys
      ENDDO
      ENDIF                                           
c  Sources hautes de SO2      
      
c
c--only GEIA SO2 emissions has high emissions
c--unit: molec/cm2/s divided by layer height (in cm) multiplied by timestep
c
      k=2                             !introducing emissions in level 2
      DO i = 1, klon
c
         IF (iregion_bb(i).GT.0) THEN
      IF(id_prec>0)   tr_seri(i,k,id_prec)=
     .              tr_seri(i,k,id_prec) + fracso2emis
     .              *scale_param_bb(iregion_bb(i))*lmt_so2bb_h(i)
     .              /zdz(i,k)/100.*pdtphys
c
      IF(id_fine>0)     tr_seri(i,k,id_fine)=tr_seri(i,k,id_fine) 
     .              + (1.-fracso2emis)
     .              *scale_param_bb(iregion_bb(i))*lmt_so2bb_h(i)
     .              *masse_ammsulfate/RNAVO/zdz(i,k)/100.*pdtphys   !g/cm3
         ENDIF
         IF (iregion_ind(i).GT.0) THEN
       IF(id_prec>0)  tr_seri(i,k,id_prec)=
     .              tr_seri(i,k,id_prec) + (fracso2emis
     .              *scale_param_ind(iregion_ind(i))*lmt_so2ff_h(i)
     .              + frach2sofso2
     .              *scale_param_ind(iregion_ind(i))*lmt_so2ff_h(i))
     .              /zdz(i,k)/100.*pdtphys
c
       IF(id_fine>0)    tr_seri(i,k,id_fine)=tr_seri(i,k,id_fine) 
     .               + (1.-fracso2emis)
     .              *scale_param_ind(iregion_ind(i))*lmt_so2ff_h(i)
     .              *masse_ammsulfate/RNAVO/zdz(i,k)/100.*pdtphys   !g/cm3
         ENDIF
c
      ENDDO

      END
