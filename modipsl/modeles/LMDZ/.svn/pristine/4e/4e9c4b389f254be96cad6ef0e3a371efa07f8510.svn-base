c Subroutine that calculates the effect of precipitation in scavenging 
c WITHIN the cloud, for large scale as well as convective precipitation
      SUBROUTINE incloud_scav_lsc(lminmax,qmin,qmax,masse,henry,kk,prfl,
     .                      psfl,pmflxr,pmflxs,zrho,zdz,t_seri,pdtphys,
     .                                   his_dhlsc,his_dhcon1,tr_seri)

      USE dimphy
      USE infotrac
      USE indice_sol_mod

      IMPLICIT NONE

#include "dimensions.h"
#include "chem.h"
#include "YOMCST.h"
#include "paramet.h"

c============================= INPUT ===================================
      REAL qmin, qmax
      REAL masse(nbtr)
      REAL henry(nbtr)         !--cste de Henry  mol/l/atm
      REAL kk(nbtr)            !--coefficient de var avec T (K)
      REAL prfl(klon,klev+1),   psfl(klon,klev+1)     !--large-scale
!      REAL prfl(klon,klev),   psfl(klon,klev)     !--large-scale
      REAL pmflxr(klon,klev+1), pmflxs(klon,klev+1)   !--convection
!      REAL pmflxr(klon,klev), pmflxs(klon,klev)   !--convection
      REAL zrho(klon,klev), zdz(klon,klev)
      REAL t_seri(klon,klev)
      LOGICAL lminmax
      REAL pdtphys
!      REAL prfl(klon,klev+1),   psfl(klon,klev+1)     !--large-scale  ! Titane
!      REAL pmflxr(klon,klev+1), pmflxs(klon,klev+1)   !--convection   ! Titane
c============================= OUTPUT ==================================
      REAL tr_seri(klon,klev,nbtr) ! traceur
      REAL aux_var1(klon,klev) ! traceur
      REAL aux_var2(klon) ! traceur
      REAL aux_var3(klon) ! traceur
      REAL his_dhlsc(klon,nbtr)        ! in-cloud scavenging lsc
      REAL his_dhcon1(klon,nbtr)       ! in-cloud scavenging con
c========================= LOCAL VARIABLES =============================      
      INTEGER it, i, j
      
      EXTERNAL minmaxqfi, inscav_spl
      DO it=1, nbtr
c
      DO i=1,klon
        aux_var2(i)=his_dhlsc(i,it)
        aux_var3(i)=his_dhcon1(i,it)
      ENDDO
      DO j=1,klev
      DO i=1,klon
        aux_var1(i,j)=tr_seri(i,j,it)
      ENDDO
      ENDDO
c      
      IF (lminmax) THEN
        CALL minmaxqfi(aux_var1,qmin,qmax,'avt inscav')
cnhl      CALL minmaxqfi(tr_seri(1,1,it),qmin,qmax,'avt inscav')
      ENDIF
c
cnhl      CALL inscav_spl(pdtphys,it,masse(it),henry(it),kk(it),0.5e-3,
cnhl     .             prfl,psfl,zrho,zdz,t_seri,tr_seri(1,1,it),
cnhl     .             his_dhlsc(1,it))
      CALL inscav_spl(pdtphys,it,masse(it),henry(it),kk(it),0.5e-3,
     .             prfl,psfl,zrho,zdz,t_seri,aux_var1,aux_var2)
c
      IF (lminmax) THEN
        CALL minmaxqfi(aux_var1,qmin,qmax,'depot humide lsc')
cnhl      CALL minmaxqfi(tr_seri(1,1,it),qmin,qmax,'depot humide lsc')
      ENDIF
c
c
c-scheme for convective in-cloud scavenging
c
cnhl      CALL inscav_spl(pdtphys,it,masse(it),henry(it),kk(it),1.e-3,
cnhl     .             pmflxr,pmflxs,zrho,zdz,t_seri,tr_seri(1,1,it),
cnhl     .             his_dhcon1(1,it))

c      print *,'JE inscav0'
c      IF (iflag_con.LT.3) THEN
c
c      print *,'JE inscav1'
c      print *,'iflag_con',iflag_con
c      CALL inscav_spl(pdtphys,it,masse(it),henry(it),kk(it),1.e-3,
c     .             pmflxr,pmflxs,zrho,zdz,t_seri,aux_var1,aux_var3)
c
cc
c      IF (lminmax) THEN
c        CALL minmaxqfi(aux_var1,qmin,qmax,'depot humide con')
ccnhl      CALL minmaxqfi(tr_seri(1,1,it),qmin,qmax,'depot humide con')
c
c      ENDIF
c
c      ENDIF ! iflag_con

c
c      print *,'JE inscav2'
      DO j=1,klev
      DO i=1,klon
        tr_seri(i,j,it)=aux_var1(i,j)
      ENDDO
      ENDDO
      DO i=1,klon
        his_dhlsc(i,it)=aux_var2(i)
        his_dhcon1(i,it)=aux_var3(i)
      ENDDO

c
      ENDDO !--boucle sur it

c      print *,'JE inscav3'
      END
