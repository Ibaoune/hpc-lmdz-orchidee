!
! $Id: writedynav_p.F 1279 2009-12-10 09:02:56Z fairhead $
!
      subroutine writehist_loc( time, vcov, ucov,teta,ppk,phi,q,
     .                          masse,ps,phis)

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
      USE ioipsl
#endif
      USE parallel_lmdz
      USE misc_mod
      USE infotrac, ONLY : nqtot, ttext
      use com_io_dyn_mod, only : histid,histvid,histuid
      USE comconst_mod, ONLY: cpp
      USE temps_mod, ONLY: itau_dyn
      
      implicit none

C
C   Ecriture du fichier histoire au format IOIPSL
C
C   Appels succesifs des routines: histwrite
C
C   Entree:
C      histid: ID du fichier histoire
C      time: temps de l'ecriture
C      vcov: vents v covariants
C      ucov: vents u covariants
C      teta: temperature potentielle
C      phi : geopotentiel instantane
C      q   : traceurs
C      masse: masse
C      ps   :pression au sol
C      phis : geopotentiel au sol
C      
C
C   Sortie:
C      fileid: ID du fichier netcdf cree
C
C   L. Fairhead, LMD, 03/99
C
C =====================================================================
C
C   Declarations
      include "dimensions.h"
      include "paramet.h"
      include "comgeom.h"
      include "description.h"
      include "iniprint.h"

C
C   Arguments
C

      REAL vcov(ijb_v:ije_v,llm),ucov(ijb_u:ije_u,llm) 
      REAL teta(ijb_u:ije_u,llm),phi(ijb_u:ije_u,llm)
      REAL ppk(ijb_u:ije_u,llm)                  
      REAL ps(ijb_u:ije_u),masse(ijb_u:ije_u,llm)                   
      REAL phis(ijb_u:ije_u)                  
      REAL q(ijb_u:ije_u,llm,nqtot)
      integer time


#ifdef CPP_IOIPSL
! This routine needs IOIPSL
C   Variables locales
C
      INTEGER,SAVE,ALLOCATABLE :: ndex2d(:),ndexu(:),ndexv(:)
      INTEGER :: iq, ii, ll
      REAL,SAVE,ALLOCATABLE :: tm(:,:)
      REAL,SAVE,ALLOCATABLE :: vnat(:,:),unat(:,:) 
      logical ok_sync
      integer itau_w
      integer :: ijb,ije,jjn
      LOGICAL,SAVE :: first=.TRUE.
!$OMP THREADPRIVATE(first)

C
C  Initialisations
C
      if (adjust) return
      
      IF (first) THEN
!$OMP BARRIER
!$OMP MASTER
        ALLOCATE(unat(ijb_u:ije_u,llm))
        ALLOCATE(vnat(ijb_v:ije_v,llm)) 
        ALLOCATE(tm(ijb_u:ije_u,llm))
        ALLOCATE(ndex2d(ijnb_u*llm))
        ALLOCATE(ndexu(ijnb_u*llm))
        ALLOCATE(ndexv(ijnb_v*llm))
        ndex2d = 0
        ndexu = 0
        ndexv = 0
!$OMP END MASTER
!$OMP BARRIER
        first=.FALSE.
      ENDIF
      
      ok_sync = .TRUE.
      itau_w = itau_dyn + time

C Passage aux composantes naturelles du vent
      call covnat_loc(llm, ucov, vcov, unat, vnat)

C
C  Appels a histwrite pour l'ecriture des variables a sauvegarder
C
C  Vents U
C

!$OMP BARRIER      
!$OMP MASTER
      ijb=ij_begin
      ije=ij_end
      jjn=jj_nb
      
      call histwrite(histuid, 'u', itau_w, unat(ijb:ije,:), 
     .               iip1*jjn*llm, ndexu)
!$OMP END MASTER      

C
C  Vents V
C
      ije=ij_end
      if (pole_sud) jjn=jj_nb-1
      if (pole_sud) ije=ij_end-iip1
!$OMP BARRIER
!$OMP MASTER      
      call histwrite(histvid, 'v', itau_w, vnat(ijb:ije,:), 
     .               iip1*jjn*llm, ndexv)
!$OMP END MASTER      


C
C  Temperature potentielle
C
      ijb=ij_begin
      ije=ij_end
      jjn=jj_nb
!$OMP MASTER      
      call histwrite(histid, 'theta', itau_w, teta(ijb:ije,:), 
     .                iip1*jjn*llm, ndexu)
!$OMP END MASTER      

C
C  Temperature
C

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      do ll=1,llm
        do ii = ijb, ije
          tm(ii,ll) = teta(ii,ll) * ppk(ii,ll)/cpp
        enddo
      enddo
!$OMP ENDDO

!$OMP MASTER      
      call histwrite(histid, 'temp', itau_w, tm(ijb:ije,:), 
     .                iip1*jjn*llm, ndexu)
!$OMP END MASTER


C
C  Geopotentiel
C
!$OMP MASTER      
      call histwrite(histid, 'phi', itau_w, phi(ijb:ije,:), 
     .                iip1*jjn*llm, ndexu)
!$OMP END MASTER


C
C  Traceurs
C
!!$OMP MASTER      
!        DO iq=1,nqtot
!          call histwrite(histid, ttext(iq), itau_w, q(ijb:ije,:,iq), 
!     .                   iip1*jjn*llm, ndexu)
!        enddo
!!$OMP END MASTER


C
C  Masse
C
!$OMP MASTER      
       call histwrite(histid, 'masse', itau_w, masse(ijb:ije,:),
     .                iip1*jjn*llm, ndexu)
!$OMP END MASTER


C
C  Pression au sol
C
!$OMP MASTER      
       call histwrite(histid, 'ps', itau_w, ps(ijb:ije), 
     .                 iip1*jjn, ndex2d)
!$OMP END MASTER

C
C  Geopotentiel au sol
C
!$OMP MASTER      
!       call histwrite(histid, 'phis', itau_w, phis(ijb:ije),
!     .                 iip1*jjn, ndex2d)
!$OMP END MASTER

C
C  Fin
C
!$OMP MASTER      
      if (ok_sync) then
        call histsync(histid)
        call histsync(histvid)
        call histsync(histuid)
      endif
!$OMP END MASTER
#else
      write(lunout,*)'writehist_loc: Needs IOIPSL to function'
#endif
! #endif of #ifdef CPP_IOIPSL
      end
