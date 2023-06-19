c Subroutine that estimates the Deposition velocities and the depostion
C for the different tracers
      subroutine deposition(vdep_oce,vdep_sic,vdep_ter,vdep_lic,pctsrf,
     .                      zrho,zdz,pdtphys,RHcl,masse,t_seri,pplay,
     .                      paprs,lminmax,qmin,qmax,
     .                        his_ds,source_tr,tr_seri)

      USE dimphy
      USE infotrac
      USE indice_sol_mod

      IMPLICIT NONE

#include "dimensions.h"
#include "chem.h"
#include "YOMCST.h"
#include "paramet.h"

c----------------------------- INPUT -----------------------------------
      LOGICAL lminmax
      REAL qmin, qmax
      REAL vdep_oce(nbtr), vdep_sic(nbtr)
      REAL vdep_ter(nbtr), vdep_lic(nbtr)      
      REAL pctsrf(klon,nbsrf)
      REAL zrho(klon,klev)        !Density of air at mid points of Z (kg/m3)
      REAL zdz(klon,klev)        
      REAL pdtphys  ! pas d'integration pour la physique (seconde)
      REAL RHcl(klon,klev)  ! humidite relativen ciel clair
      REAL t_seri(klon,klev)  ! temperature
      REAL pplay(klon,klev)  ! pression pour le mileu de chaque couche (en Pa)
      REAL paprs(klon, klev+1)    !pressure at interface of layers Z (Pa)
      REAL masse(nbtr)
                                          
c----------------------------- OUTPUT ----------------------------------
      REAL his_ds(klon,nbtr)                                          
      REAL source_tr(klon,nbtr)
      REAL tr_seri(klon, klev,nbtr) !conc of tracers
c--------------------- INTERNAL VARIABLES ------------------------------      
      INTEGER i, it
      REAL vdep        !sed. velocity

      DO it=1, nbtr
      DO i=1, klon
          vdep=vdep_oce(it)*pctsrf(i,is_oce)+
     .         vdep_sic(it)*pctsrf(i,is_sic)+
     .         vdep_ter(it)*pctsrf(i,is_ter)+
     .         vdep_lic(it)*pctsrf(i,is_lic)
c--Unit: molec/m2/s for it=1 to nbtr-3, mg/m2/s for it=nbtr-2 to nbtr
          source_tr(i,it)=source_tr(i,it)
     .                    -vdep*tr_seri(i,1,it)*zrho(i,1)/1.e2
          his_ds(i,it)=vdep*tr_seri(i,1,it)*zrho(i,1)/1.e2
     .                 /RNAVO*masse(it)*1.e3               ! mg/m2/s
      ENDDO
      ENDDO
c
      END
