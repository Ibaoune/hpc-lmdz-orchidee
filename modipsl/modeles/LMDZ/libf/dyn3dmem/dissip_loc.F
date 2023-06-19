!
! $Id: $
!
      SUBROUTINE dissip_loc( vcov,ucov,teta,p, dv,du,dh )
c
      USE parallel_lmdz
      USE write_field_loc
      USE dissip_mod, ONLY: dissip_allocate
      USE comconst_mod, ONLY: dtdiss
      IMPLICIT NONE


c ..  Avec nouveaux operateurs star :  gradiv2 , divgrad2, nxgraro2  ...
c                                 (  10/01/98  )

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c   Dissipation horizontale
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

      include "dimensions.h"
      include "paramet.h"
      include "comgeom.h"
      include "comdissnew.h"
      include "comdissipn.h"

c   Arguments:
c   ----------

      REAL,INTENT(IN) :: vcov(ijb_v:ije_v,llm) ! covariant meridional wind
      REAL,INTENT(IN) :: ucov(ijb_u:ije_u,llm) ! covariant zonal wind
      REAL,INTENT(IN) :: teta(ijb_u:ije_u,llm) ! potential temperature
      REAL,INTENT(IN) :: p(ijb_u:ije_u,llmp1) ! interlayer pressure
      ! tendencies (.../s) on covariant winds and potential temperature
      REAL,INTENT(OUT) :: dv(ijb_v:ije_v,llm)
      REAL,INTENT(OUT) :: du(ijb_u:ije_u,llm)
      REAL,INTENT(OUT) :: dh(ijb_u:ije_u,llm)

c   Local:
c   ------

      REAL gdx(ijb_u:ije_u,llm),gdy(ijb_v:ije_v,llm)
      REAL grx(ijb_u:ije_u,llm),gry(ijb_v:ije_v,llm)
      REAL te1dt(llm),te2dt(llm),te3dt(llm)
      REAL deltapres(ijb_u:ije_u,llm)

      INTEGER l,ij

      REAL  SSUM
      integer :: ijb,ije
      
      LOGICAl,SAVE :: first=.TRUE.
!$OMP THREADPRIVATE(first)

      IF (first) THEN
        CALL dissip_allocate
        first=.FALSE.
      ENDIF
c-----------------------------------------------------------------------
c   initialisations:
c   ----------------

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
         te1dt(l) = tetaudiv(l) * dtdiss
         te2dt(l) = tetaurot(l) * dtdiss
         te3dt(l) = tetah(l)    * dtdiss
      ENDDO
c$OMP END DO NOWAIT
c      CALL initial0( ijp1llm, du )
c      CALL initial0( ijmllm , dv )
c      CALL initial0( ijp1llm, dh )
     
      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
        du(ijb:ije,l)=0
        dh(ijb:ije,l)=0
      ENDDO
c$OMP END DO NOWAIT
      
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
        dv(ijb:ije,l)=0
      ENDDO
c$OMP END DO NOWAIT
     
c-----------------------------------------------------------------------
c   Calcul de la dissipation:
c   -------------------------

c   Calcul de la partie   grad  ( div ) :
c   -------------------------------------
      
     
      
      IF(lstardis) THEN
c      IF (.FALSE.) THEN
         CALL gradiv2_loc( llm,ucov,vcov,nitergdiv,gdx,gdy )
      ELSE
!         CALL gradiv_p ( llm,ucov,vcov,nitergdiv,gdx,gdy )
      ENDIF

#ifdef DEBUG_IO    
      call WriteField_u('gdx',gdx)
      call WriteField_v('gdy',gdy)
#endif

      ijb=ij_begin
      ije=ij_end
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
         if (pole_nord) then
           DO ij = 1, iip1
              gdx(     ij ,l) = 0.
           ENDDO
         endif
         
         if (pole_sud) then
           DO ij = 1, iip1
              gdx(ij+ip1jm,l) = 0.
           ENDDO
         endif
         
         if (pole_nord) ijb=ij_begin+iip1
         DO ij = ijb,ije
            du(ij,l) = du(ij,l) - te1dt(l) *gdx(ij,l)
         ENDDO

         if (pole_nord) ijb=ij_begin
         DO ij = ijb,ije
            dv(ij,l) = dv(ij,l) - te1dt(l) *gdy(ij,l)
         ENDDO

       ENDDO
c$OMP END DO NOWAIT
c   calcul de la partie   n X grad ( rot ):
c   ---------------------------------------

      IF(lstardis) THEN
c      IF (.FALSE.) THEN
         CALL nxgraro2_loc( llm,ucov, vcov, nitergrot,grx,gry )
      ELSE
!         CALL nxgrarot_p( llm,ucov, vcov, nitergrot,grx,gry )
      ENDIF

#ifdef DEBUG_IO    
      call WriteField_u('grx',grx)
      call WriteField_v('gry',gry)
#endif


      ijb=ij_begin
      ije=ij_end
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l=1,llm
         
         if (pole_nord) then
           DO ij = 1, iip1
              grx(ij,l) = 0.
           ENDDO
         endif
         
         if (pole_nord) ijb=ij_begin+iip1
         DO ij = ijb,ije
            du(ij,l) = du(ij,l) - te2dt(l) * grx(ij,l)
         ENDDO
         
         if (pole_nord) ijb=ij_begin
         DO ij =  ijb, ije
            dv(ij,l) = dv(ij,l) - te2dt(l) * gry(ij,l)
         ENDDO
      
      ENDDO
c$OMP END DO NOWAIT

c   calcul de la partie   div ( grad ):
c   -----------------------------------

        
      IF(lstardis) THEN
c      IF (.FALSE.) THEN
    
      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       DO l = 1, llm
          DO ij = ijb, ije
            deltapres(ij,l) = AMAX1( 0.,  p(ij,l) - p(ij,l+1) )
          ENDDO
       ENDDO
c$OMP END DO NOWAIT
         CALL divgrad2_loc( llm,teta, deltapres  ,niterh, gdx )
      ELSE
!         CALL divgrad_p ( llm,teta, niterh, gdx        )
      ENDIF

#ifdef DEBUG_IO    
      call WriteField_u('gdx',gdx)
#endif


      ijb=ij_begin
      ije=ij_end
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1,llm
         DO ij = ijb,ije
            dh( ij,l ) = dh( ij,l ) - te3dt(l) * gdx( ij,l )
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

      RETURN
      END
