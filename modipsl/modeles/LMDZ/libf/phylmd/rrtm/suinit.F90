SUBROUTINE SUINIT(klon,klev)
#ifdef DOC

!     **** *SUINIT* - SCM initialization.

!     Purpose.
!     --------

!     **   Interface.
!     ----------

!     Explicit arguments :    None.
!     --------------------

!     Implicit arguments :    None.
!     --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!     Eric Bazile, Francois Bouyssel et Jean-Marcel Piriou

!     Modifications.
!     --------------
!     Original :97-02-01
!               Jozef Vivoda, SHMI: calling sequence as in 3D model
!                                   and ECMWF setup
!     2001-11-27 P. Marquet : several printout on listing (NULOUT=15)

!     ------------------------------------------------------------------
#endif

USE PARKIND1  ,ONLY : JPIM     ,JPRB
!#include "tsmbkind.h"

USE PARDIM, ONLY : JPMXLE
USE YOMCT0B  , ONLY : LECMWF
USE YOMRIP   , ONLY : NINDAT   ,NSSSSS
USE YOMDIM  
USE YOMDPHY  
! MPL 29042010: NDLNPR,RHYDR0 non initialises et pour ne pas mettre tout sudyn.F90
USE YOMDYN  , ONLY : TSTEP , NDLNPR , RHYDR0        ! MPL 29042010
!USE YOMEVOL  , ONLY : TECH     ,FREQFS   ,FREQFE   , FREQDDH
!USE YOMCT0   , ONLY : LFROG
! quelques ajouts qui viennent de suallo
USE YOMGEM   , ONLY : VDELA    , VDELB   ,VC       ,NLOEN    ,NLOENG  ,NGPTOT
USE YOMSTA   , ONLY : STZ      ,STPREH   ,STPRE    ,STPHI    ,STTEM   ,STDEN
USE YOEAERD  , ONLY : CVDAES   ,CVDAEL   ,CVDAEU   ,CVDAED
USE YOEOVLP  , ONLY : RA1OVLP
USE YOECLD   , ONLY : CETA
USE YOECND   , ONLY : CEVAPCU
USE YOMTOPH  , ONLY : RMESOU   ,RMESOT   ,RMESOQ
USE YOMGC    , ONLY : GEMU     ,GELAM    ,GELAT    ,GECLO    ,GESLO    ,GM       ,GAW


IMPLICIT NONE
LOGICAL LLTRACE, LLDEBUG
integer klon,klev
CHARACTER*200 CFICP
CHARACTER*200 CFLUX
CHARACTER*200 CLIST
CHARACTER*200 CFDDH
CHARACTER*80 CNMEXP


LLTRACE=.TRUE.
LLDEBUG=.TRUE.

!     ------------------------
!     *    READ NAMELISTS.
!     ------------------------

!----------------------------------------------------------------
! Elements indispensables de SUNAM pour faire tourner RRTM dans LMDZ
!-------------------------------------------------------------------
CFICP='Profile'
CFLUX='Output'
CLIST='Listing'
CFDDH='DHFDL'
CNMEXP='SCM'
TSTEP=450
! MPL 29042010 - RHYDR0 - upper boundary contition for hydrostatic
RHYDR0=1._JPRB
! MPL 29042010
! NDLNPR : NDLNPR=0: conventional formulation of delta, i.e. ln(P(l)/P(l-1)).
!          NDLNPR=1: formulation of delta used in non hydrostatic model,
NDLNPR=0
print *,'SUINIT: RHYDR0 NDLNPR',RHYDR0,NDLNPR

!----------------------------------------------------------------
! Elements indispensables de SUDIM pour faire tourner RRTM dans LMDZ
!-------------------------------------------------------------------
NDLON=klon
NFLEVG=klev
NPROMA=klon

!-------------------------------------------------------------------
!JV    Initialize constants
!     ---------------------
!JV
IF (LLTRACE)  WRITE(*,*) " coucou SUINIT : avant SUCST"
WRITE(*,FMT='('' ---------------- '')')
WRITE(*,FMT='(''     SUCST : '')')
WRITE(*,FMT='('' ---------------- '')')
NINDAT=20090408      !!!!! A REVOIR (MPL)
NSSSSS=0  ! LMDZ demarre tjrs a 00h -- MPL 15.04.09
CALL SUCST(6,NINDAT,NSSSSS,1)
print *,'SUINIT: NINDAT, NSSSSS',NINDAT, NSSSSS

IF (LLDEBUG) THEN
WRITE(*,FMT='(''  SUINIT / apres : SUCST '')')
ENDIF


!     ------------------------
!     *    ALLOCATES RECUPERES DE SUALLO
!     ------------------------
ALLOCATE(VDELA  (MAX(JPMXLE,NFLEVG)))
ALLOCATE(VDELB  (MAX(JPMXLE,NFLEVG)))
VDELB = 0  !ym missing init
ALLOCATE( VC      (NFLEVG) )
VC = 0    !ym missing init
ALLOCATE( NLOEN   (NPROMA) )
ALLOCATE( NLOENG   (NPROMA) )
ALLOCATE( STZ     (NFLEVG) )
ALLOCATE( CVDAES  (NFLEVG+1))
ALLOCATE( CVDAEL  (NFLEVG+1))
ALLOCATE( CVDAEU  (NFLEVG+1))
ALLOCATE( CVDAED  (NFLEVG+1))
ALLOCATE(RA1OVLP(NFLEVG))

ALLOCATE(STPREH(0:NFLEVG)) ! Nouvel ajout MPL 22062010
ALLOCATE(STPRE(NFLEVG))
ALLOCATE(STPHI(NFLEVG))
ALLOCATE(STTEM(NFLEVG))
ALLOCATE(STDEN(NFLEVG))

ALLOCATE(CETA(NFLEVG))    ! Nouvel ajout MPL 28062010
ALLOCATE(CEVAPCU(NFLEVG))
ALLOCATE(RMESOU(NFLEVG))
ALLOCATE(RMESOT(NFLEVG))
ALLOCATE(RMESOQ(NFLEVG))

!     ------------------------
!     *    ALLOCATES RECUPERES DE SUGEM2
!     ------------------------

ALLOCATE(GEMU   (NGPTOT)) ! Nouvel ajout MPL 28062010
ALLOCATE(GELAM  (NGPTOT))
ALLOCATE(GELAT  (NGPTOT))
ALLOCATE(GECLO  (NGPTOT))
ALLOCATE(GESLO  (NGPTOT))
ALLOCATE(GM     (NGPTOT))
ALLOCATE(GAW    (NGPTOT))
!  
!     ------------------------------------------------------------------

END SUBROUTINE SUINIT
