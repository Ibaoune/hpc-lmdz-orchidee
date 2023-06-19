      SUBROUTINE gastoparticle(pdtphys,zdz,zrho,xlat,pplay,t_seri,
     .         id_prec,id_fine,
     .         tr_seri,his_g2pgas ,his_g2paer )
cnhl     .                         fluxso4chem, flux_sparam_sulf,

      USE dimphy
      USE infotrac
c      USE indice_sol_mod

      IMPLICIT NONE 
c
#include "dimensions.h"
#include "chem.h"
#include "chem_spla.h"
#include "YOMCST.h"
#include "YOECUMF.h"
c
      REAL pdtphys
      REAL zrho(klon,klev)
      REAL zdz(klon,klev)
      REAL tr_seri(klon,klev,nbtr)   ! traceurs
      REAL tend                 ! tendance par espece
      REAL xlat(klon)       ! latitudes pour chaque point
      REAL pi
c   JE: 20140120
      REAL his_g2pgas(klon)
      REAL his_g2paer(klon)
      REAL tendincm3(klon,klev)
      REAL tempvar(klon,klev)      
      REAL pplay(klon,klev)
      REAL t_seri(klon,klev)
      REAL tend2d(klon,klev)
      INTEGER id_prec,id_fine
c
c------------------------- Scaling Parameter --------------------------
c
c      REAL scale_param_so4(klon)  !Scaling parameter for sulfate

      INTEGER i, k
      REAL tau_chem     !---chemical lifetime in s
c
c------------------------- Variables to save --------------------------
c
cnhl      REAL fluxso4chem(klon,klev)
cnhl      REAL flux_sparam_sulf(klon,klev)

c======================================================================
      pi=atan(1.)*4.
c
      IF (id_prec>0 .AND. id_fine>0) THEN
      DO k = 1, klev
      DO i = 1, klon
c
c        tau_chem=scale_param_so4(i)*86400.*(8.-5.*cos(xlat(i)*pi/180.))    !tchemfctn2
cnhl        tau_chem=86400.*(8.-5.*cos(xlat(i)*pi/180.))    !tchemfctn2
        tau_chem=86400.*(5.-4.*cos(xlat(i)*pi/180.))    !
        tend=tr_seri(i,k,id_prec)*(1.-exp(-pdtphys/tau_chem)) ! Sulfate production
cnhl        tend=(1.-exp(-pdtphys/tau_chem))
cnhl        tend=scale_param_so4(i) !as this it works
c      
        tr_seri(i,k,id_prec) =tr_seri(i,k,id_prec) - tend
        tr_seri(i,k,id_fine) =tr_seri(i,k,id_fine) + 
     .                      tend/RNAVO*masse_ammsulfate  !--gAER/KgAir
        tend2d(i,k)=tend
c
cnhl        fluxso4chem(i,k) = tend/RNAVO*masse_ammsulfate
cnhl        flux_sparam_sulf(i,k) = tend/RNAVO*masse_ammsulfate
      ENDDO
      ENDDO
   


        tempvar=tend2d
         CALL kg_to_cm3(pplay,t_seri,tempvar)
        tendincm3=tempvar

      DO k = 1, klev
      DO i = 1, klon

c        his_g2pgas(i) = his_g2pgas(i) + tendincm3(i,k)*1e6*zdz(i,k)/pdtphys
        his_g2paer(i) = his_g2paer(i) +      
     .    tendincm3(i,k)/RNAVO*masse_ammsulfate*1.e3*
     .    1.e6*zdz(i,k)/pdtphys    ! mg/m2/s
        his_g2pgas(i) = his_g2paer(i)*masse_s/masse_ammsulfate ! mg-S/m2/s

      ENDDO
      ENDDO
      ENDIF

c
      RETURN 
      END
