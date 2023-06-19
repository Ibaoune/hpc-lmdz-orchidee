c Subroutine that calculates the effect of precipitation in scavenging
c BELOW the cloud, for large scale as well as convective precipitation
      SUBROUTINE blcloud_scav(lminmax,qmin,qmax,pdtphys,prfl,psfl,
     .                        pmflxr,pmflxs,zdz,alpha_r,alpha_s,masse,
     .                                  his_dhbclsc,his_dhbccon,tr_seri)

      USE dimphy
      USE indice_sol_mod
      USE infotrac
      IMPLICIT NONE

#include "dimensions.h"
#include "chem.h"
#include "YOMCST.h"
#include "paramet.h"

c============================= INPUT ===================================
      REAL qmin,qmax
      REAL pdtphys  ! pas d'integration pour la physique (seconde)
!      REAL prfl(klon,klev),   psfl(klon,klev)     !--large-scale
!      REAL pmflxr(klon,klev), pmflxs(klon,klev)   !--convection
      REAL alpha_r(nbtr)!--coefficient d'impaction pour la pluie
      REAL alpha_s(nbtr)!--coefficient d'impaction pour la neige      
      REAL masse(nbtr)
      LOGICAL lminmax
      REAL zdz(klon,klev)
      REAL prfl(klon,klev+1),   psfl(klon,klev+1)     !--large-scale  ! Titane
      REAL pmflxr(klon,klev+1), pmflxs(klon,klev+1)   !--convection   ! Titane
c============================= OUTPUT ==================================
      REAL tr_seri(klon,klev,nbtr) ! traceur
      REAL aux_var1(klon,klev) ! traceur
      REAL aux_var2(klon,klev) ! traceur
      REAL his_dhbclsc(klon,nbtr), his_dhbccon(klon,nbtr)
c========================= LOCAL VARIABLES =============================      
      INTEGER it, k, i, j
      REAL d_tr(klon,klev,nbtr)

      EXTERNAL minmaxqfi, bcscav_spl
      
      DO it=1, nbtr
c
      DO j=1,klev
      DO i=1,klon
        aux_var1(i,j)=tr_seri(i,j,it)
        aux_var2(i,j)=d_tr(i,j,it)
      ENDDO
      ENDDO
c
cnhl      CALL bcscav_spl(pdtphys,prfl,psfl,alpha_r(it),alpha_s(it),
cnhl     .                tr_seri(1,1,it),d_tr(1,1,it))
      CALL bcscav_spl(pdtphys,prfl,psfl,alpha_r(it),alpha_s(it),
     .                aux_var1,aux_var2)
c
      DO j=1,klev
      DO i=1,klon
        tr_seri(i,j,it)=aux_var1(i,j)
        d_tr(i,j,it)=aux_var2(i,j)
      ENDDO
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         tr_seri(i,k,it) = tr_seri(i,k,it) + d_tr(i,k,it)
         his_dhbclsc(i,it)=his_dhbclsc(i,it)-d_tr(i,k,it)/RNAVO*
     .                masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys  !--mgS/m2/s
              
      ENDDO
      ENDDO
c
      DO i=1,klon
      DO j=1,klev
        aux_var1(i,j)=tr_seri(i,j,it)
        aux_var2(i,j)=d_tr(i,j,it)
      ENDDO
      ENDDO
c
      IF (lminmax) THEN
        CALL minmaxqfi(aux_var1,qmin,qmax,'depot humide bc lsc')
cnhl      CALL minmaxqfi(tr_seri(1,1,it),qmin,qmax,'depot humide bc lsc')
      ENDIF
c
c-scheme for convective scavenging
c
cnhl      CALL bcscav_spl(pdtphys,pmflxr,pmflxs,alpha_r(it),alpha_s(it),
cnhl     .                tr_seri(1,1,it),d_tr(1,1,it))


      CALL bcscav_spl(pdtphys,pmflxr,pmflxs,alpha_r(it),alpha_s(it),
     .                aux_var1,aux_var2)


c
      DO i=1,klon
      DO j=1,klev
        tr_seri(i,j,it)=aux_var1(i,j)
        d_tr(i,j,it)=aux_var2(i,j)
      ENDDO
      ENDDO
c
      DO k = 1, klev
      DO i = 1, klon
         tr_seri(i,k,it) = tr_seri(i,k,it) + d_tr(i,k,it)
         his_dhbccon(i,it)=his_dhbccon(i,it)-d_tr(i,k,it)/RNAVO*
     .                masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys    !--mgS/m2/s
      ENDDO
      ENDDO
c
      IF (lminmax) THEN
        DO j=1,klev
        DO i=1,klon
          aux_var1(i,j)=tr_seri(i,j,it)
        ENDDO
        ENDDO
        CALL minmaxqfi(aux_var1,qmin,qmax,'depot humide bc con')
cnhl      CALL minmaxqfi(tr_seri(1,1,it),qmin,qmax,'depot humide bc con')
        DO j=1,klev
        DO i=1,klon
          tr_seri(i,j,it)=aux_var1(i,j)
        ENDDO
        ENDDO
      ENDIF
c
c
      ENDDO !--boucle sur it
c
      END
