!
! $Id: initdynav_p.F 1279 2009-12-10 09:02:56Z fairhead $
!
      subroutine initdynav_loc(day0,anne0,tstep,t_ops,t_wrt)

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
       USE IOIPSL
#endif
       USE parallel_lmdz
       use Write_field
       use misc_mod
       USE infotrac
       use com_io_dyn_mod, only : histaveid,histvaveid,histuaveid,       &
     &        dynhistave_file,dynhistvave_file,dynhistuave_file
       USE comconst_mod, ONLY: pi
       USE comvert_mod, ONLY: presnivs
       USE temps_mod, ONLY: itau_dyn
       
       implicit none

C
C   Routine d'initialisation des ecritures des fichiers histoires LMDZ
C   au format IOIPSL. Initialisation du fichier histoire moyenne.
C
C   Appels succesifs des routines: histbeg
C                                  histhori
C                                  histver
C                                  histdef
C                                  histend
C
C   Entree:
C
C      day0,anne0: date de reference
C      tstep : frequence d'ecriture
C      t_ops: frequence de l'operation pour IOIPSL
C      t_wrt: frequence d'ecriture sur le fichier
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

C   Arguments
C
      integer*4 day0, anne0
      real tstep, t_ops, t_wrt

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
C   Variables locales
C
      integer tau0
      real zjulian
      integer iq
      real rlong(iip1,jjp1), rlat(iip1,jjp1)
      integer uhoriid, vhoriid, thoriid
      integer zvertiid,zvertiidv,zvertiidu
      integer ii,jj
      integer zan, dayref
      integer :: jjb,jje,jjn

! definition du domaine d'ecriture pour le rebuild

      INTEGER,DIMENSION(2) :: ddid
      INTEGER,DIMENSION(2) :: dsg
      INTEGER,DIMENSION(2) :: dsl
      INTEGER,DIMENSION(2) :: dpf
      INTEGER,DIMENSION(2) :: dpl
      INTEGER,DIMENSION(2) :: dhs
      INTEGER,DIMENSION(2) :: dhe 
      
      INTEGER :: dynhistave_domain_id
      INTEGER :: dynhistvave_domain_id
      INTEGER :: dynhistuave_domain_id
      
      if (adjust) return

C
C  Initialisations
C
      pi = 4. * atan (1.)
C
C  Appel a histbeg: creation du fichier netcdf et initialisations diverses
C         

      zan = anne0
      dayref = day0
      CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)
      tau0 = itau_dyn
      
      do jj = 1, jjp1
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj)  = rlatu(jj) * 180. / pi
        enddo
      enddo


! Creation de 3 fichiers pour les differentes grilles horizontales
! Restriction de IOIPSL: seulement 2 coordonnees dans le meme fichier
! Grille Scalaire       

      jjb=jj_begin
      jje=jj_end
      jjn=jj_nb

      ddid=(/ 1,2 /)
      dsg=(/ iip1,jjp1 /)
      dsl=(/ iip1,jjn /)
      dpf=(/ 1,jjb /)
      dpl=(/ iip1,jje /)
      dhs=(/ 0,0 /)
      dhe=(/ 0,0 /)


      call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, 
     .                 'box',dynhistave_domain_id)
             
      call histbeg(dynhistave_file,iip1, rlong(:,1), jjn,
     .             rlat(1,jjb:jje), 1, iip1, 1, jjn,tau0,
     .             zjulian, tstep, thoriid,
     .             histaveid,dynhistave_domain_id)


C  Creation du fichier histoire pour les grilles en V et U (oblige pour l'instant,
C  IOIPSL ne permet pas de grilles avec des nombres de point differents dans
C  un meme fichier)
! Grille V

      jjb=jj_begin
      jje=jj_end
      jjn=jj_nb
      IF (pole_sud) jjn=jjn-1
      IF (pole_sud) jje=jje-1
      
      do jj = jjb, jje
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj) = rlatv(jj) * 180. / pi
        enddo
      enddo

      ddid=(/ 1,2 /)
      dsg=(/ iip1,jjm /)
      dsl=(/ iip1,jjn /)
      dpf=(/ 1,jjb /)
      dpl=(/ iip1,jje /)
      dhs=(/ 0,0 /)
      dhe=(/ 0,0 /)


      call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, 
     .                 'box',dynhistvave_domain_id)

      call histbeg(dynhistvave_file,iip1, rlong(:,1), jjn,
     .             rlat(1,jjb:jje), 1, iip1, 1, jjn,tau0,
     .             zjulian, tstep, vhoriid,
     .             histvaveid,dynhistvave_domain_id)
     
! Grille U

      do jj = 1, jjp1
        do ii = 1, iip1
          rlong(ii,jj) = rlonu(ii) * 180. / pi
          rlat(ii,jj) = rlatu(jj) * 180. / pi
        enddo
      enddo

      jjb=jj_begin
      jje=jj_end
      jjn=jj_nb

      ddid=(/ 1,2 /)
      dsg=(/ iip1,jjp1 /)
      dsl=(/ iip1,jjn /)
      dpf=(/ 1,jjb /)
      dpl=(/ iip1,jje /)
      dhs=(/ 0,0 /)
      dhe=(/ 0,0 /)


      call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, 
     .                 'box',dynhistuave_domain_id)
             
      call histbeg(dynhistuave_file,iip1, rlong(:,1), jjn,
     .             rlat(1,jjb:jje), 1, iip1, 1, jjn,tau0,
     .             zjulian, tstep, uhoriid,
     .             histuaveid,dynhistuave_domain_id)
     
     
C
C  Appel a histvert pour la grille verticale
C
      call histvert(histaveid,'presnivs','Niveaux Pression
     &     approximatifs','mb',llm, presnivs/100., zvertiid,'down')
      call histvert(histuaveid,'presnivs','Niveaux Pression
     &     approximatifs','mb',llm, presnivs/100., zvertiidv,'down')
      call histvert(histvaveid,'presnivs','Niveaux Pression
     &     approximatifs','mb',llm, presnivs/100., zvertiidu,'down')

C
C  Appels a histdef pour la definition des variables a sauvegarder
C
C  Vents U
C
      jjn=jj_nb
      call histdef(histuaveid, 'u', 'vent u moyen ',
     .             'm/s', iip1, jjn, uhoriid, llm, 1, llm, zvertiidu,
     .             32, 'ave(X)', t_ops, t_wrt)

C
C  Vents V
C
      if (pole_sud) jjn=jj_nb-1
      call histdef(histvaveid, 'v', 'vent v moyen',
     .             'm/s', iip1, jjn, vhoriid, llm, 1, llm, zvertiidv,
     .             32, 'ave(X)', t_ops, t_wrt)

C
C  Temperature
C
      jjn=jj_nb
      call histdef(histaveid, 'temp', 'temperature moyenne', 'K',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', t_ops, t_wrt)
C
C  Temperature potentielle
C
      call histdef(histaveid, 'theta', 'temperature potentielle', 'K',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', t_ops, t_wrt)


C
C  Geopotentiel
C
      call histdef(histaveid, 'phi', 'geopotentiel moyen', '-',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', t_ops, t_wrt)
C
C  Traceurs
C
!        DO iq=1,nqtot
!          call histdef(histaveid, ttext(iq), ttext(iq), '-',
!     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
!     .             32, 'ave(X)', t_ops, t_wrt)
!        enddo
C
C  Masse
C
      call histdef(histaveid, 'masse', 'masse moyenne', 'kg',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', t_ops, t_wrt)
C
C  Pression au sol
C
      call histdef(histaveid, 'ps', 'pression naturelle au sol', 'Pa',
     .             iip1, jjn, thoriid, 1, 1, 1, -99,
     .             32, 'ave(X)', t_ops, t_wrt)
C
C  Geopotentiel au sol
C
!      call histdef(histaveid, 'phis', 'geopotentiel au sol', '-',
!     .             iip1, jjn, thoriid, 1, 1, 1, -99,
!     .             32, 'ave(X)', t_ops, t_wrt)
C
C  Fin
C
      call histend(histaveid)
      call histend(histuaveid)
      call histend(histvaveid)
#else
      write(lunout,*)'initdynav_loc: Needs IOIPSL to function'
#endif
! #endif of #ifdef CPP_IOIPSL
      end
