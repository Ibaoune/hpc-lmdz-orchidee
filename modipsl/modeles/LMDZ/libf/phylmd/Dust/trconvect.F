c Subroutine that computes the convective mixing and transport
      SUBROUTINE trconvect(pplay,t_seri,pdtphys,pmfu,pmfd,pen_u,pde_u,
     .          pen_d,pde_d,paprs,zdz,xconv,qmin,qmax,lminmax,masse,
     .                                                dtrconv,tr_seri)

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
      REAL xconv(nbtr), masse(nbtr)
      REAL pplay(klon,klev)    ! pression pour le mileu de chaque couche (en Pa)
      REAL t_seri(klon,klev)   ! temperature       
      REAL zdz(klon,klev)      ! zdz
      REAL paprs(klon,klev+1)  ! pression pour chaque inter-couche (en Pa)
      REAL pmfu(klon,klev)     ! flux de masse dans le panache montant
      REAL pmfd(klon,klev)     ! flux de masse dans le panache descendant
      REAL pen_u(klon,klev)    ! flux entraine dans le panache montant
      REAL pde_u(klon,klev)    ! flux detraine dans le panache montant
      REAL pen_d(klon,klev)    ! flux entraine dans le panache descendant
      REAL pde_d(klon,klev)    ! flux detraine dans le panache descendant  
      LOGICAL lminmax
      REAL pdtphys
c============================= OUTPUT ==================================
      REAL aux_var1(klon,klev)
      REAL aux_var2(klon,klev)
      REAL tr_seri(klon,klev,nbtr) ! traceur
      REAL dtrconv(klon,nbtr) ! traceur
c========================= LOCAL VARIABLES =============================
      INTEGER it, k, i, j
      REAL d_tr(klon,klev,nbtr)

      EXTERNAL nflxtr, tiedqneg, minmaxqfi
      
      DO it=1, nbtr
c
      DO i=1, klon
        dtrconv(i,it)=0.0
      ENDDO
      DO i=1,klon
      DO j=1,klev
        aux_var1(i,j)=tr_seri(i,j,it)
        aux_var2(i,j)=d_tr(i,j,it)
      ENDDO
      ENDDO
                                                 
c
cnhl      CALL nflxtr(pdtphys, pmfu, pmfd, pen_u, pde_u, pen_d, pde_d,
cnhl     .            pplay, paprs, tr_seri(1,1,it), d_tr(1,1,it) )
      CALL nflxtr(pdtphys, pmfu, pmfd, pen_u, pde_u, pen_d, pde_d,
     .            pplay, paprs, aux_var1, aux_var2 )
c
      CALL tiedqneg(paprs,aux_var1, aux_var2)
cnhl      CALL tiedqneg(paprs,tr_seri(1,1,it), d_tr(1,1,it))
      DO i=1,klon
      DO j=1,klev
        tr_seri(i,j,it)=aux_var1(i,j)
        d_tr(i,j,it)=aux_var2(i,j)
      ENDDO
      ENDDO
c
      DO k = 1, klev
      DO i = 1, klon
        IF (d_tr(i,k,it).LT.0.) THEN
          tr_seri(i,k,it)=tr_seri(i,k,it)+d_tr(i,k,it)
        ELSE
          tr_seri(i,k,it)=tr_seri(i,k,it)+d_tr(i,k,it)*xconv(it)
        ENDIF
      ENDDO
      ENDDO
c
cnhl      CALL kg_to_cm3(pplay,t_seri,d_tr(1,1,it))
      CALL kg_to_cm3(pplay,t_seri,aux_var2)
      DO i=1,klon
      DO j=1,klev
        d_tr(i,j,it)=aux_var2(i,j)
      ENDDO
      ENDDO

      DO k = 1, klev
      DO i = 1, klon
        IF (d_tr(i,k,it).GE.0.) THEN
        dtrconv(i,it)=dtrconv(i,it)+(1.-xconv(it))*d_tr(i,k,it)
     .                /RNAVO*masse(it)*1.e3*1.e6*zdz(i,k)/pdtphys
        ENDIF
      ENDDO
      ENDDO

      IF (lminmax) THEN
        DO i=1,klon
        DO j=1,klev
          aux_var1(i,j)=tr_seri(i,j,it)
        ENDDO
        ENDDO
        CALL minmaxqfi(aux_var1,qmin,qmax,'apr convection')
cnhl      CALL minmaxqfi(tr_seri(1,1,it),qmin,qmax,'apr convection')
        DO i=1,klon
        DO j=1,klev
          tr_seri(i,j,it)=aux_var1(i,j)
        ENDDO
        ENDDO
      ENDIF
c
      ENDDO

      END
