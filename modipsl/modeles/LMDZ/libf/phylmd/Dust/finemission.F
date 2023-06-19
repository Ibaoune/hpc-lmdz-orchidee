C This subroutine calculates the emissions of BLACK CARBON and ORGANIC
C MATTER
      SUBROUTINE finemission(zdz,pdtphys,zalt,kminbc,kmaxbc,
     .                       scale_param_bb,scale_param_ff,
     .                       iregion_ind,iregion_bb,
     .                       nbreg_ind,nbreg_bb,
     .                       lmt_bcff,lmt_bcnff,lmt_bcbb_l,lmt_bcbb_h,
     .                       lmt_bcba,lmt_omff,lmt_omnff,lmt_ombb_l,
     .                       lmt_ombb_h,lmt_omnat,lmt_omba,id_fine,
     .                                    flux_sparam_bb,flux_sparam_ff,
     .                                        source_tr,flux_tr,tr_seri)

      USE dimphy
      USE indice_sol_mod
      USE infotrac
!      USE phytracr_spl_mod, ONLY : nbreg_dust, nbreg_ind, nbreg_bb
      IMPLICIT NONE

#include "dimensions.h"
#include "chem.h"
#include "YOMCST.h"
#include "paramet.h"

      INTEGER i, k, kminbc, kmaxbc
c============================= INPUT ===================================
      REAL pdtphys  ! pas d'integration pour la physique (seconde)
      REAL zalt(klon,klev)
      REAL zdz(klon,klev)
c
c------------------------- Scaling Parameters --------------------------
c   
      INTEGER nbreg_ind,nbreg_bb
      INTEGER iregion_ind(klon)  !Defines regions for SO2, BC & OM
      INTEGER iregion_bb(klon)   !Defines regions for SO2, BC & OM
      REAL scale_param_bb(nbreg_bb) !Scaling parameter for biomas burning
      REAL scale_param_ff(nbreg_ind) !Scaling parameter for industrial emissions (fossil fuel)
      INTEGER id_fine
c============================= OUTPUT ==================================
      REAL source_tr(klon,nbtr)
      REAL flux_tr(klon,nbtr)
      REAL tr_seri(klon,klev,nbtr) ! traceur
      REAL flux_sparam_bb(klon), flux_sparam_ff(klon)
c========================= LOCAL VARIABLES =============================
      REAL zzdz
c------------------------- BLACK CARBON emissions ----------------------
      REAL lmt_bcff(klon)       ! emissions de BC fossil fuels
      REAL lmt_bcnff(klon)      ! emissions de BC non-fossil fuels
      REAL lmt_bcbb_l(klon)     ! emissions de BC biomass basses
      REAL lmt_bcbb_h(klon)     ! emissions de BC biomass hautes
      REAL lmt_bcba(klon)       ! emissions de BC bateau
c------------------------ ORGANIC MATTER emissions ---------------------
      REAL lmt_omff(klon)     ! emissions de OM fossil fuels
      REAL lmt_omnff(klon)    ! emissions de OM non-fossil fuels
      REAL lmt_ombb_l(klon)   ! emissions de OM biomass basses
      REAL lmt_ombb_h(klon)   ! emissions de OM biomass hautes
      REAL lmt_omnat(klon)    ! emissions de OM Natural
      REAL lmt_omba(klon)     ! emissions de OM bateau
                                
      EXTERNAL condsurfc                                
c========================================================================
c                        LOW LEVEL EMISSIONS
c========================================================================
            
c corresponds to bc_source.EQ.3

      DO i=1,klon
         IF (iregion_ind(i).GT.0) THEN
       IF(id_fine>0)    source_tr(i,id_fine)=source_tr(i,id_fine)+
     .                (scale_param_ff(iregion_ind(i))*lmt_bcff(i)+     !g/m2/s
     .                 scale_param_ff(iregion_ind(i))*lmt_omff(i)
     .                 )     * 1.e4                        !g/m2/s
c
      IF(id_fine>0)     flux_tr(i,id_fine)=flux_tr(i,id_fine)+
     .                (scale_param_ff(iregion_ind(i))*lmt_bcff(i)+     !mg/m2/s
     .                 scale_param_ff(iregion_ind(i))*lmt_omff(i)
     .                 )     * 1.e4 *1.e3                  !mg/m2/s
c
           flux_sparam_ff(i)= flux_sparam_ff(i) +
     .                     scale_param_ff(iregion_ind(i))*
     .                     ( lmt_bcff(i)+lmt_omff(i))
     .                     *1.e4*1.e3
         ENDIF
         IF (iregion_bb(i).GT.0) THEN
       IF(id_fine>0)    source_tr(i,id_fine)=source_tr(i,id_fine)+
     .                (scale_param_bb(iregion_bb(i))*lmt_bcbb_l(i)+   !g/m2/s
     .                 scale_param_bb(iregion_bb(i))*lmt_ombb_l(i)    !g/m2/s
     .                 )     * 1.e4                        !g/m2/s
c
       IF(id_fine>0)    flux_tr(i,id_fine)=flux_tr(i,id_fine)+
     .                (scale_param_bb(iregion_bb(i))*lmt_bcbb_l(i)+   !mg/m2/s
     .                 scale_param_bb(iregion_bb(i))*lmt_ombb_l(i)+   !mg/m2/s
     .                 scale_param_bb(iregion_bb(i))*lmt_bcbb_h(i)+   !mg/m2/s
     .                 scale_param_bb(iregion_bb(i))*lmt_ombb_h(i)    !mg/m2/s
     .                 )     * 1.e4 *1.e3                  !mg/m2/s
c
           flux_sparam_bb(i)=flux_sparam_bb(i) +
     .                   scale_param_bb(iregion_bb(i))*(lmt_bcbb_l(i) +
     .                   lmt_bcbb_h(i) + lmt_ombb_l(i) + lmt_ombb_h(i))
     .                   *1.e4*1.e3
         ENDIF
       IF(id_fine>0)  source_tr(i,id_fine)=source_tr(i,id_fine)+
     .                (lmt_bcnff(i)+lmt_bcba(i)+lmt_omnff(i)+
     .                 lmt_omnat(i)+lmt_omba(i))     * 1.e4           !g/m2/s
c
       IF(id_fine>0)  flux_tr(i,id_fine)=flux_tr(i,id_fine)+
     .                (lmt_bcnff(i)+lmt_omnff(i)+lmt_omnat(i)+
     .                 lmt_omba(i)+lmt_bcba(i))     * 1.e4 *1.e3      !mg/m2/s
c
         flux_sparam_ff(i)= flux_sparam_ff(i) + 
     .                      (lmt_omba(i)+lmt_bcba(i))*1.e4*1.e3
      ENDDO

c========================================================================
c                        HIGH LEVEL EMISSIONS
c========================================================================
                             
c  Sources hautes de BC/OM

c
c HIGH LEVEL EMISSIONS OF SO2 ARE IN PRECUREMISSION.F
c
      k=2                             !introducing emissions in level 2
cnhl      DO i = 1, klon
c
cnhl      tr_seri(i,k,id_fine)=tr_seri(i,k,id_fine)+scale_param_ff(iregion_ind(i))*
cnhl     .               (lmt_bcff_h(i)+lmt_omff_h(i))/zdz(i,k)/100.*pdtphys
c
cnhl      ENDDO

      DO k=kminbc, kmaxbc
      DO i = 1, klon
          zzdz=zalt(i,kmaxbc+1)-zalt(i,kminbc)
c
         IF (iregion_bb(i) .GT.0) THEN
        IF(id_fine>0)   tr_seri(i,k,id_fine)=tr_seri(i,k,id_fine)+
     .              (scale_param_bb(iregion_bb(i))*lmt_bcbb_h(i)+
     .               scale_param_bb(iregion_bb(i))*lmt_ombb_h(i))
     .                              /zzdz/100.*pdtphys
         ENDIF
c
      ENDDO
      ENDDO
c
      END
