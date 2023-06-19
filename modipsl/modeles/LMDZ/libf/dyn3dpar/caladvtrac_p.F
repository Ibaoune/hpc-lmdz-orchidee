!
! $Id: caladvtrac_p.F 2597 2016-07-22 06:44:47Z emillour $
!
c
c
            SUBROUTINE caladvtrac_p(q,pbaru,pbarv ,
     *                   p ,masse, dq ,  teta,
     *                   flxw, pk, iapptrac)
      USE parallel_lmdz
      USE infotrac, ONLY : nqtot
      USE control_mod, ONLY : iapp_tracvl,planet_type
      USE comconst_mod, ONLY: dtvr
c
      IMPLICIT NONE
c
c     Auteurs:   F.Hourdin , P.Le Van, F.Forget, F.Codron  
c
c     F.Codron (10/99) : ajout humidite specifique pour eau vapeur
c=======================================================================
c
c       Shema de  Van Leer
c
c=======================================================================


#include "dimensions.h"
#include "paramet.h"

c   Arguments:
c   ----------
      REAL pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm),masse(ip1jmp1,llm)
      REAL p( ip1jmp1,llmp1),q( ip1jmp1,llm,nqtot)
      real :: dq( ip1jmp1,llm,nqtot)
      REAL teta( ip1jmp1,llm),pk( ip1jmp1,llm)
      REAL               :: flxw(ip1jmp1,llm)

      integer ijb,ije,jjb,jje

c  ..................................................................
c
c  .. dq n'est utilise et dimensionne que pour l'eau  vapeur et liqu.
c
c  ..................................................................
c
c   Local:
c   ------

      INTEGER ij,l, iq, iapptrac
      REAL finmasse(ip1jmp1,llm), dtvrtrac
      
cc
c
C initialisation
cym      ijb=ij_begin
cym      ije=ij_end

      
cym      dq(ijb:ije,1:llm,1:2)=q(ijb:ije,1:llm,1:2)

c  test des valeurs minmax
cc        CALL minmaxq(q(1,1,1),1.e33,-1.e33,'Eau vapeur (a) ')
cc        CALL minmaxq(q(1,1,2),1.e33,-1.e33,'Eau liquide(a) ')

c   advection
c      print *,'appel a advtrac'

      CALL advtrac_p( pbaru,pbarv, 
     *             p,  masse,q,iapptrac, teta,
     .             flxw, pk)

         goto 9999
         IF( iapptrac.EQ.iapp_tracvl ) THEN
c
cc          CALL minmaxq(q(1,1,1),1.e33,-1.e33,'Eau vapeur     ')
cc          CALL minmaxq(q(1,1,2),1.e33,-1.e33,'Eau liquide    ')

cc     ....  Calcul  de deltap  qu'on stocke dans finmasse   ...
c
          DO l = 1, llm
           DO ij = ijb, ije
             finmasse(ij,l) =  p(ij,l) - p(ij,l+1) 
           ENDDO
          ENDDO

	  if (planet_type.eq."earth") then
! Earth-specific treatment of first 2 tracers (water)
            CALL qminimum_p( q, 2, finmasse )
	  endif


cym   --> le reste ne set a rien
          goto 9999
	  
c          CALL SCOPY   ( ip1jmp1*llm, masse, 1, finmasse,       1 )
          finmasse(ijb:ije,:)=masse(ijb:ije,:)         
          
          jjb=jj_begin
          jje=jj_end
          CALL filtreg_p ( finmasse ,jjb,jje,  jjp1,  llm, 
     *                     -2, 2, .TRUE., 1 )
c
c   *****  Calcul de dq pour l'eau , pour le passer a la physique ******
c   ********************************************************************
c
          dtvrtrac = iapp_tracvl * dtvr
c
           DO iq = 1 , 2
            DO l = 1 , llm
             DO ij = ijb,ije
             dq(ij,l,iq) = ( q(ij,l,iq) - dq(ij,l,iq) ) * finmasse(ij,l)
     *                               /  dtvrtrac
             ENDDO
            ENDDO
           ENDDO
c
         ELSE
cym   --> le reste ne set a rien
          goto 9999
	  
           DO iq = 1 , 2
           DO l  = 1, llm
             DO ij = ijb,ije
              dq(ij,l,iq)  = 0.
             ENDDO
           ENDDO
           ENDDO

         ENDIF
c


c  ... On appelle  qminimum uniquement  pour l'eau vapeur et liquide  ..

 
 9999 RETURN
      END


