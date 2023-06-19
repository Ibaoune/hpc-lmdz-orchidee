MODULE YOMDPHY

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!     DIMENSION DES TABLEAUX POINT DE GRILLE PHYSIQUE

!     NVXP : number of variables in the generic EXTRP.
!     NVXP2: number of variables in the generic XTRP2.

!     NCXP : number of levels in EXTRP
!     NCSI : number of sea-ice levels
!     NCSNEC: number of snow levels in EC physics (specific for ECMWF) 
!     NTILES: number of surface tiles

!     NVEXTR : number of variables in the generic VEXTRA (extra-fields)
!     NVEXTRDYN : number of extra-fields comming from the dynamics
!     NVXTR2 : number of variables in the generic VEXTR2

!     NCEXTR : number of levels in the generic VEXTRA

!     NTSL : nombre de types de sols nus.
!     NTSV : nombre de types de vegetations.
!     NTOZ1D: 1 si representation 1D des NVCLIS variables , 0 sinon
!     NTOZ2D: 1 si representation 2D des NVCLIS variables , 0 sinon
!     NTOZ3D: 1 si representation 3D des NVCLIS variables , 0 sinon
!     NLOA : nombre de longueurs d'ondes pour le spectre d'albedo.
!     NLOE : nombre de longueurs d'ondes pour le spectre d'emissivite.
!     NTSSG : number of surface temperatures for subgrid diagnostics

!     NCHAC: nombre de champs a accumuler.
!     NCHIN: nombre de champs instantanes.
!     NSIRA: nombre d'intervales spectraux pour les diagnostiques de
!            flux radiatif a p=0 et p=ps.

!     NVTEND: number of tendencies used in 2.order scheme

!     LTPROF: .T. if more than 1 vertical layer in deep soil

INTEGER(KIND=JPIM) :: NVXP
INTEGER(KIND=JPIM) :: NVXP2
INTEGER(KIND=JPIM) :: NCXP
INTEGER(KIND=JPIM) :: NCSI
INTEGER(KIND=JPIM) :: NCSNEC
INTEGER(KIND=JPIM) :: NTILES
INTEGER(KIND=JPIM) :: NTSL
INTEGER(KIND=JPIM) :: NVEXTR
INTEGER(KIND=JPIM) :: NVEXTRDYN
INTEGER(KIND=JPIM) :: NVXTR2
INTEGER(KIND=JPIM) :: NCEXTR
INTEGER(KIND=JPIM) :: NVCLIS
INTEGER(KIND=JPIM) :: NTOZ1D
INTEGER(KIND=JPIM) :: NTOZ2D
INTEGER(KIND=JPIM) :: NTOZ3D
INTEGER(KIND=JPIM) :: NTSSG
INTEGER(KIND=JPIM) :: NTVG
INTEGER(KIND=JPIM) :: NLOA
INTEGER(KIND=JPIM) :: NLOE
INTEGER(KIND=JPIM) :: NCHAC
INTEGER(KIND=JPIM) :: NCHIN
INTEGER(KIND=JPIM) :: NSIRA
INTEGER(KIND=JPIM) :: NVTEND
LOGICAL            :: LTPROF
!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(ltprof,ncextr,nchac,nchin,ncsi,ncsnec,ncxp,nloa,nloe,nsira,ntiles,ntoz1d,ntoz2d,ntoz3d,ntsl,ntssg)
!$OMP THREADPRIVATE(ntvg,nvclis,nvextr,nvextrdyn,nvtend,nvxp,nvxp2,nvxtr2)
END MODULE YOMDPHY
