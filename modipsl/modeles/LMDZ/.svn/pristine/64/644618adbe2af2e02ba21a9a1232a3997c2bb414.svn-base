      SUBROUTINE tiedqneg (pres_h,q,d_q)
c
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): CG (LGGE/CNRS) date: 19950201
c            O. Boucher (LOA/CNRS) date 19961125
c Objet:  Correction eventuelle des valeurs negatives d'humidite
c induites par le schema de convection de Tiedke 
c======================================================================
c Arguments:
c pres_h--input-R-la valeur de la pression aux interfaces
c q-------input-R-quantite de traceur
c d_q-----input-output-R-increment du traceur
c======================================================================
c
#include "dimensions.h"
c #include "dimphy.h"
      REAL pres_h(klon,klev+1)
      REAL q(klon,klev)
      REAL d_q(klon,klev)
      INTEGER nb_neg
      INTEGER i, l
c
      REAL qmin
      PARAMETER (qmin=0.0)
c
      DO l = klev,2,-1
        nb_neg = 0
        DO i = 1,klon
          IF (q(i,l)+d_q(i,l).LT.qmin) THEN 
          nb_neg = nb_neg + 1
          d_q(i,l-1) = d_q(i,l-1) + (q(i,l)+d_q(i,l)-qmin)
     .       *(pres_h(i,l)-pres_h(i,l+1))/(pres_h(i,l-1)-pres_h(i,l))
            d_q(i,l) = qmin - q(i,l)
          ENDIF
        ENDDO
c        IF (nb_neg.NE.0) THEN 
c        PRINT *,'niveau ', l,' ' , nb_neg, ' valeurs negatives'
c        ENDIF
      ENDDO
c
      DO l = 1, klev-1
        nb_neg = 0
        DO i = 1,klon
          IF (q(i,l)+d_q(i,l).LT.qmin) THEN
          nb_neg = nb_neg + 1
          d_q(i,l+1) = d_q(i,l+1) + (q(i,l)+d_q(i,l)-qmin)
     .      *(pres_h(i,l)-pres_h(i,l+1))/(pres_h(i,l+1)-pres_h(i,l+2))
          d_q(i,l) = qmin - q(i,l)
          ENDIF
        ENDDO
c        IF (nb_neg.NE.0) THEN 
c        PRINT *,'niveau ', l,' ' , nb_neg, ' valeurs negatives'
c        ENDIF
      ENDDO
c
      l = klev
      DO i = 1,klon
        IF (q(i,l)+d_q(i,l).LT.qmin) THEN 
          d_q(i,l) = qmin - q(i,l)
        ENDIF
      ENDDO
c
      RETURN
      END
