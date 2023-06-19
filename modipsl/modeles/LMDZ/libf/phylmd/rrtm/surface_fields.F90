MODULE SURFACE_FIELDS

!     Purpose.
!     --------

!      SURFACE_FIELDS contains data structures and manipulation routines
!      for the surface (physics) fields in the IFS            

!      This module is a mix of declarations, type definitions and 
!       subroutines linked with surface fields. There are four parts:
!       1/ Declaration of dimensions (including some parameter variables).
!       2/ Definition of types.
!       3/ Declarations:
!          Declaration of variables SP_[group], YSP_[group]D, YSP_[group]
!           (prognostic surface fields).
!          Declaration of variables SD_[group], YSD_[group]D, YSD_[group]
!           (diagnostic surface fields).
!       4/ Some routines linked to the surface data flow:
!          * INI_SFLP3: Initialize 3-D surface field group
!          * SETUP_SFLP3: Setup 3-D surface field
!          * INI_SFLP2: Initialize 2-D surface field group
!          * SETUP_SFLP2: Setup 2-D surface field
!          * GPPOPER: Operations on prognostic surface fields
!          * GPOPER: Operations on ALL surface groups
!          * GPOPER_2: Operations on 2-D surface groups
!          * GPOPER_3: Operations on 3-D surface groups
!          * SURF_STORE: Store all surface fields
!          * SURF_RESTORE: Restore all surface fields
!          * ALLO_SURF: Allocate surface field arrays
!          * DEALLO_SURF: Deallocate surface field arrays

!     Author.
!     -------
!     Mats Hamrud(ECMWF)

!     Modifications.
!     --------------
!        Original : 2006-07-01
!        Modifications:
!        K. Yessad (25 Oct 2006): rephase ALARO0 contribution.
!        K. Yessad (26 Oct 2006): add missing comments.

!-------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMDIM    ,ONLY : NPROMA, NGPBLKS, NUNDEFLD
USE YOMLUN    ,ONLY : NULOUT, NULERR
USE YOMCT0    ,ONLY : LTWOTL
USE YOMDYN   , ONLY : REPSP1
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
IMPLICIT NONE
SAVE

!#include "abor1.intfb.h"
!     -------------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPMAXSFLDS=100 ! Max number of fields in individual group
INTEGER(KIND=JPIM), PARAMETER :: JPMAXSTRAJ=100 ! Dimension of NSTRAJGRIB
INTEGER(KIND=JPIM) :: NSURF=0              ! Number of surf var.
INTEGER(KIND=JPIM) :: NSURFL=0             ! Number of surf flds (fields*levels)
INTEGER(KIND=JPIM) :: NDIMSURF=0           ! Total of surf var (includes timelevels etc)
INTEGER(KIND=JPIM) :: NDIMSURFL=0          ! Total dimension of all surface variables
INTEGER(KIND=JPIM) :: NPROGSURF=0          ! Number of prognostic surf var.
INTEGER(KIND=JPIM) :: NPROGSURFL=0         ! Number of prognostic surf flds (fields*levels)
INTEGER(KIND=JPIM) :: NOFFTRAJ             ! Offset in surf trajectory
INTEGER(KIND=JPIM) :: NOFFTRAJ_CST         ! Offset in "constant" surf trajectory
INTEGER(KIND=JPIM) :: NPTRSURF             ! Used by routine GPOPER
INTEGER(KIND=JPIM) :: NSTRAJGRIB(JPMAXSTRAJ) ! Used in trajectory setup

REAL(KIND=JPRB),ALLOCATABLE   :: SURF_STORE_ARRAY(:,:,:) ! Backup array for surf (see routineSURF_STORE )
! General type defintions

! 2D surface field structure
TYPE TYPE_SURF_MTL_2D
INTEGER(KIND=JPIM) :: MP                   ! Basic field pointer
INTEGER(KIND=JPIM) :: MP0                  ! Field pointer timelevel  0 (prognostic fields)
INTEGER(KIND=JPIM) :: MP9                  ! Field pointer timelevel -1 (prognostic fields)
INTEGER(KIND=JPIM) :: MP1                  ! Field pointer timelevel +1 (prognostic fields)
INTEGER(KIND=JPIM) :: MP5                  ! Field pointer trajectory
INTEGER(KIND=JPIM) :: IGRBCODE             ! GRIB parameter code (default: -999)
CHARACTER(LEN=16)  :: CNAME                ! ARPEGE field name   (default: all spaces)
REAL(KIND=JPRB)    :: REFVALI              ! Default value       (default: 0.0)
INTEGER(KIND=JPIM) :: NREQIN               ! -1 - initial value from default (default)
                                           ! +1 - initial value from reading file
                                           !  0 - no initial value
INTEGER(KIND=JPIM) :: ITRAJ                !  0 not in trajectory (default)
                                           !  1 in trajectory
                                           !  2 in "constant" trajectory
END TYPE TYPE_SURF_MTL_2D

! 3D surface field structure
TYPE TYPE_SURF_MTL_3D
INTEGER(KIND=JPIM) :: MP   ! Basic field pointer
INTEGER(KIND=JPIM) :: MP0  ! Field pointer timelevel  0 (prognostic fields)
INTEGER(KIND=JPIM) :: MP9  ! Field pointer timelevel -1 (prognostic fields)
INTEGER(KIND=JPIM) :: MP1  ! Field pointer timelevel +1 (prognostic fields)
INTEGER(KIND=JPIM) :: MP5  ! Field pointer trajectory
INTEGER(KIND=JPIM),POINTER :: IGRBCODE(:)  ! GRIB parameter code (default: -999)
CHARACTER(LEN=16) ,POINTER :: CNAME(:)     ! ARPEGE field name   (default: all spaces)
REAL(KIND=JPRB)   ,POINTER :: REFVALI(:)   ! Default value       (default: 0.0)
INTEGER(KIND=JPIM),POINTER :: NREQIN(:)    ! -1 - initial value from default (default)
                                           ! +1 - initial value from reading file
                                           !  0 - no initial value
INTEGER(KIND=JPIM) :: ITRAJ                !  0 not in trajectory (default)
                                           !  1 in trajectory
                                           !  2 in "constant" trajectory
END TYPE TYPE_SURF_MTL_3D

! Descriptor pertaining to group
TYPE TYPE_SURF_GEN
INTEGER(KIND=JPIM) :: NUMFLDS         ! Number of field in group
INTEGER(KIND=JPIM) :: NDIM            ! Field dimenion
INTEGER(KIND=JPIM) :: NLEVS           ! Number of levels (for multi level groups)
INTEGER(KIND=JPIM) :: IPTR            ! Internal use
INTEGER(KIND=JPIM) :: IPTR5           ! Internal use
INTEGER(KIND=JPIM) :: NDIM5           ! Dimension of trajectory array
INTEGER(KIND=JPIM) :: NOFFTRAJ        ! Internal use
INTEGER(KIND=JPIM) :: NOFFTRAJ_CST    ! Internal use
CHARACTER(LEN=16)  :: CGRPNAME        ! Name of group (for prints)
LOGICAL            :: L3D             ! TRUE if multi-level field (3-D)
LOGICAL            :: LMTL            ! TRUE if prognostic field (multi time level)
END TYPE TYPE_SURF_GEN

! Type descriptor for derived type for communicating with GPOPER (see below)
TYPE TYPE_SFL_COMM
INTEGER(KIND=JPIM) :: IGRBCODE
LOGICAL            :: L_OK 
CHARACTER(LEN=16)  :: CNAME
INTEGER(KIND=JPIM) :: IFLDNUM
REAL(KIND=JPRB)    :: VALUE
INTEGER(KIND=JPIM) :: IPTRSURF
INTEGER(KIND=JPIM) :: ICODES(JPMAXSFLDS)
INTEGER(KIND=JPIM) :: ICOUNT
END TYPE TYPE_SFL_COMM

! Group specific type definitions

! * Group SB=SOILB: soil prognostic quantities for the different reservoirs
!    (four reservoirs at ECMWF, deep reservoir at METEO-FRANCE):
TYPE TYPE_SFL_SOILB
TYPE(TYPE_SURF_MTL_3D),POINTER :: YT    ! temperature
TYPE(TYPE_SURF_MTL_3D),POINTER :: YQ    ! liquid water content
TYPE(TYPE_SURF_MTL_3D),POINTER :: YTL   ! ice water content (for MF)
TYPE(TYPE_SURF_MTL_3D),POINTER :: YSB(:)
END TYPE TYPE_SFL_SOILB

! * Group SG=SNOWG: surface snow prognostic quantities:
TYPE TYPE_SFL_SNOWG
TYPE(TYPE_SURF_MTL_2D),POINTER :: YF    ! content of surface snow
TYPE(TYPE_SURF_MTL_2D),POINTER :: YA    ! snow albedo
TYPE(TYPE_SURF_MTL_2D),POINTER :: YR    ! snow density
TYPE(TYPE_SURF_MTL_2D),POINTER :: YT    ! total albedo (diagnostic for MF for LVGSN)
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSG(:)
END TYPE TYPE_SFL_SNOWG

! * Group RR=RESVR: surface prognostic quantities (ECMWF) or
!   surface + superficial reservoir prognostic quantities (MF):
!   Remark:
!    at ECMWF there are 4 soil reservoirs and there is a
!    clear distinction between the soil reservoirs (group SOILB)
!    and the surface (group RESVR);
!    at METEO-FRANCE there is a deep reservoir (group SOILB) and a
!    superficial reservoir (group RESVR):
!    - there is a skin surface temperature (Ts) which is the temperature at the
!      interface surface/superficial reservoir (and not two separate quantities
!      for superficial reservoir and surface)
!    - there is a skin surface water content (denoted by Wl) and a superficial
!      reservoir water content (denoted by Ws).
!    - there is a superficial reservoir ice content but no surface ice content.
!    (remark k.y.: it would have been more logical to use group name
!    RESVR for internal reservoirs and group name SOILB for surface!).
TYPE TYPE_SFL_RESVR
TYPE(TYPE_SURF_MTL_2D),POINTER :: YT    ! skin temperature (Ts)
TYPE(TYPE_SURF_MTL_2D),POINTER :: YW    ! skin water content (Wskin) at ECMWF
                                        ! superficial reservoir water content (Ws) at MF
TYPE(TYPE_SURF_MTL_2D),POINTER :: YFC   ! skin water content (Wl) at MF
TYPE(TYPE_SURF_MTL_2D),POINTER :: YIC   ! superficial reservoir ice
TYPE(TYPE_SURF_MTL_2D),POINTER :: YFP1  ! interpolated Ts for 2nd part of 927-FULLPOS
TYPE(TYPE_SURF_MTL_2D),POINTER :: YRR(:)
END TYPE TYPE_SFL_RESVR

! * Group WS=WAVES: surface prognostic quantities over sea:
TYPE TYPE_SFL_WAVES
TYPE(TYPE_SURF_MTL_2D),POINTER :: YWS(:)
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCHAR ! Charnock constant
END TYPE TYPE_SFL_WAVES

! * Group EP=EXTRP: extra 3-d prognostic fields:
TYPE TYPE_SFL_EXTRP
TYPE(TYPE_SURF_MTL_3D),POINTER :: YEP(:)
END TYPE TYPE_SFL_EXTRP

! * Group X2=XTRP2: extra 2-d prognostic fields:
!   (is used for precipitation fields in CANARI)
TYPE TYPE_SFL_XTRP2
TYPE(TYPE_SURF_MTL_2D),POINTER :: YX2(:)
END TYPE TYPE_SFL_XTRP2

! * Group CI=CANRI: 2-d prognostic fields for CANARI:
TYPE TYPE_SFL_CANRI
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCI(:)
END TYPE TYPE_SFL_CANRI

! * Group VF=VARSF: climatological/geographical diagnostic fields:
TYPE TYPE_SFL_VARSF
TYPE(TYPE_SURF_MTL_2D),POINTER :: YZ0F    ! gravity * surface roughness length
TYPE(TYPE_SURF_MTL_2D),POINTER :: YALBF   ! surface shortwave albedo
TYPE(TYPE_SURF_MTL_2D),POINTER :: YEMISF  ! surface longwave emissivity
TYPE(TYPE_SURF_MTL_2D),POINTER :: YGETRL  ! standard deviation of orography
TYPE(TYPE_SURF_MTL_2D),POINTER :: YITM    ! land-sea mask
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVEG    ! vegetation cover
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVRLAN  ! anisotropy of the sub-grid scale orography
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVRLDI  ! angle of the direction of orography with the x axis
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSIG    ! characteristic orographic slope
TYPE(TYPE_SURF_MTL_2D),POINTER :: YALBSF  ! soil shortwave albedo
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCONT   ! fraction of land
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSST    ! (open) sea surface temperature
TYPE(TYPE_SURF_MTL_2D),POINTER :: YLZ0H   ! logarithm of roughness length for heat
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCVL    ! low vegetation cover
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCVH    ! high vegetation cover
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTVL    ! low vegetation type
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTVH    ! high vegetation type
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCI     ! sea ice fraction
TYPE(TYPE_SURF_MTL_2D),POINTER :: YUCUR   ! U-component of the ocean current
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVCUR   ! V-component of the ocean current
TYPE(TYPE_SURF_MTL_2D),POINTER :: YZ0RLF  ! gravity * vegetation roughness length
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCO2O   ! oceanic CO2 flux
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCO2B   ! biosphere CO2 flux
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCO2A   ! anthropogenic CO2 flux
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSDFOR  ! SD filtered orography
TYPE(TYPE_SURF_MTL_2D),POINTER :: YALUVP  ! MODIS-derived parallel albedo for shortwave radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YALUVD  ! MODIS-derived diffuse albedo for shortwave radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YALNIP  ! MODIS-derived parallel albedo for longwave radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YALNID  ! MODIS-derived diffuse albedo for longwave radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSF6    ! anthropogenic SF6 flux
TYPE(TYPE_SURF_MTL_2D),POINTER :: YFP1    ! surface orography in the 2nd part of FULLPOS-927
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVF(:)
END TYPE TYPE_SFL_VARSF

! * Group VP=VCLIP: deep soil diagnostic fields:
TYPE TYPE_SFL_VCLIP
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTPC    ! climatological deep layer temperature
TYPE(TYPE_SURF_MTL_2D),POINTER :: YWPC    ! climatological deep layer moisture
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVP(:)
END TYPE TYPE_SFL_VCLIP

! * Group VV=VCLIV: vegetation diagnostic fields:
TYPE TYPE_SFL_VCLIV
TYPE(TYPE_SURF_MTL_2D),POINTER :: YARG    ! silt percentage within soil
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSAB    ! percentage of sand within the soil
TYPE(TYPE_SURF_MTL_2D),POINTER :: YD2     ! soil depth
TYPE(TYPE_SURF_MTL_2D),POINTER :: YIVEG   ! type of vegetation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YRSMIN  ! stomatal minimum resistance
TYPE(TYPE_SURF_MTL_2D),POINTER :: YLAI    ! leaf area index
TYPE(TYPE_SURF_MTL_2D),POINTER :: YHV     ! resistance to evapotranspiration
TYPE(TYPE_SURF_MTL_2D),POINTER :: YZ0H    ! gravity * roughness length for heat
TYPE(TYPE_SURF_MTL_2D),POINTER :: YALS    ! albedo of bare ground
TYPE(TYPE_SURF_MTL_2D),POINTER :: YALV    ! albedo of vegetation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVV(:)
END TYPE TYPE_SFL_VCLIV

! * Group VN=VCLIN: cloudiness diagnostic predictors:
TYPE TYPE_SFL_VCLIN
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTOP    ! index of convective cloud top
TYPE(TYPE_SURF_MTL_2D),POINTER :: YBAS    ! index of convective cloud base
TYPE(TYPE_SURF_MTL_2D),POINTER :: YACPR   ! averaged convective precipitaion rate
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVN(:)
END TYPE TYPE_SFL_VCLIN

! * Group VH=VCLIH: convective cloud diagnostic fields:
TYPE TYPE_SFL_VCLIH
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTCCH   ! total convective cloudiness
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSCCH   ! convective cloud summit
TYPE(TYPE_SURF_MTL_2D),POINTER :: YBCCH   ! convective cloud base
TYPE(TYPE_SURF_MTL_2D),POINTER :: YPBLH   ! PBL height
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSPSH   ! variable for prognostic convection scheme (ALARO)
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVH(:)
END TYPE TYPE_SFL_VCLIH

! * Group VA=VCLIA: aerosol diagnostic fields:
TYPE TYPE_SFL_VCLIA
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSEA    ! aerosol: sea
TYPE(TYPE_SURF_MTL_2D),POINTER :: YLAN    ! aerosol: land
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSOO    ! aerosol: soot
TYPE(TYPE_SURF_MTL_2D),POINTER :: YDES    ! aerosol: desert
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSUL    ! aerosol: sulfate
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVOL    ! aerosol: volcano
TYPE(TYPE_SURF_MTL_2D),POINTER :: YNUD    ! aerosol: nudging
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVA(:)
END TYPE TYPE_SFL_VCLIA

! * Group VG=VCLIG: ice-coupler diagnostic fields:
TYPE TYPE_SFL_VCLIG
TYPE(TYPE_SURF_MTL_2D),POINTER :: YICFR   ! sea-ice fraction
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSOUP   ! upward solar flux over sea-ice
TYPE(TYPE_SURF_MTL_2D),POINTER :: YIRUP   ! upward IR flux over sea-ice
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCHSS   ! sensible heat over sea-ice
TYPE(TYPE_SURF_MTL_2D),POINTER :: YEVAP   ! evaporation over sea-ice
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTAUX   ! U-component of stress over sea-ice
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTAUY   ! V-component of stress over sea-ice
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVG(:)
END TYPE TYPE_SFL_VCLIG

! * Group VC=VO3ABC: A,B and C (Climatological ozone profiles) diagnostic fields:
TYPE TYPE_SFL_VO3ABC
TYPE(TYPE_SURF_MTL_2D),POINTER :: YA      ! A climatological ozone profile
TYPE(TYPE_SURF_MTL_2D),POINTER :: YB      ! B climatological ozone profile
TYPE(TYPE_SURF_MTL_2D),POINTER :: YC      ! C climatological ozone profile
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVC(:)
END TYPE TYPE_SFL_VO3ABC

! * Group VD=VDIAG: (ECMWF) diagnostic fields:
TYPE TYPE_SFL_VDIAG
TYPE(TYPE_SURF_MTL_2D),POINTER :: YLSP   !Large scale precipitation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCP    !Convective precipitation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSF    !Snowfall
TYPE(TYPE_SURF_MTL_2D),POINTER :: YBLD   !Boundary layer dissipation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSSHF  !Surface sensible heat flux
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSLHF  !Surface latent heat flux
TYPE(TYPE_SURF_MTL_2D),POINTER :: YMSL   !Mean sea level pressure
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTCC   !Total cloud cover
TYPE(TYPE_SURF_MTL_2D),POINTER :: Y10U   !U-wind at 10 m
TYPE(TYPE_SURF_MTL_2D),POINTER :: Y10V   !V-wind at 10 m
TYPE(TYPE_SURF_MTL_2D),POINTER :: Y2T    !Temperature at 2 m
TYPE(TYPE_SURF_MTL_2D),POINTER :: Y2D    !Dewpoint temperature at 2 m
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSSR   !Surface solar radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSTR   !Surface thermal radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTSR   !Top solar radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTTR   !Top thermal radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YEWSS  !Instantaneous surface U-wind stress
TYPE(TYPE_SURF_MTL_2D),POINTER :: YNSSS  !Instantaneous surface V-wind stress
TYPE(TYPE_SURF_MTL_2D),POINTER :: YE     !Water evaporation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCCC   !Convective cloud cover
TYPE(TYPE_SURF_MTL_2D),POINTER :: YLCC   !Low cloud cover
TYPE(TYPE_SURF_MTL_2D),POINTER :: YMCC   !Medium cloud cover
TYPE(TYPE_SURF_MTL_2D),POINTER :: YHCC   !High cloud cover
TYPE(TYPE_SURF_MTL_2D),POINTER :: YLGWS  !Zonal gravity wave stress
TYPE(TYPE_SURF_MTL_2D),POINTER :: YMGWS  !Meridian gravity wave stress
TYPE(TYPE_SURF_MTL_2D),POINTER :: YGWD   !Gravity wave dissipation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YMX2T  !Maximum temperature at 2 m
TYPE(TYPE_SURF_MTL_2D),POINTER :: YMN2T  !Minimum temperature at 2 m
TYPE(TYPE_SURF_MTL_2D),POINTER :: YRO    !Runoff
TYPE(TYPE_SURF_MTL_2D),POINTER :: YALB   !(surface shortwave) albedo
TYPE(TYPE_SURF_MTL_2D),POINTER :: YIEWSS !Instantaneous surface zonal component of stress
TYPE(TYPE_SURF_MTL_2D),POINTER :: YINSSS !Instantaneous surface meridian component of stress
TYPE(TYPE_SURF_MTL_2D),POINTER :: YISSHF !Instantaneous surface heat flux
TYPE(TYPE_SURF_MTL_2D),POINTER :: YIE    !Instantaneous surface moisture flux
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCSF   !Convective snow fall
TYPE(TYPE_SURF_MTL_2D),POINTER :: YLSSF  !Large scale snowfall
TYPE(TYPE_SURF_MTL_2D),POINTER :: YZ0F   !Gravity * surface roughness length
TYPE(TYPE_SURF_MTL_2D),POINTER :: YLZ0H  !Logarithm of z0 times heat flux
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTCW   !Total water content in a vertical column
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTCWV  !Total water vapor content in a vertical column
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTCLW  !Total liquid water content in a vertical column
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTCIW  !Total ice water content in a vertical column
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSSRD  !Downward surface solar radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSTRD  !Downward surface thermic radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YBLH   !Height of boundary layer
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSUND  !Sunshine duration
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSPAR  !Surface downward PARadiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSUVB  !Surface downward UV-B radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YCAPE  !Conv.avail.potential energy (CAPE)
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTSRC  !Top solar radiation clear sky
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTTRC  !Top thermal radiation clear sky
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSSRC  !Surface solar radiation clear sky
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSTRC  !Surface thermal radiation clear sky
TYPE(TYPE_SURF_MTL_2D),POINTER :: YES    !Evaporation of snow
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSMLT  !Snow melt
TYPE(TYPE_SURF_MTL_2D),POINTER :: Y10FG  !Wind gust at 10 m (max since previous pp)
TYPE(TYPE_SURF_MTL_2D),POINTER :: YLSPF  !Large scale precipitation fraction
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTCO3  !Total ozone content in a vertical column
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVIMD  !Vertically integrated mass divergence
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSPARC !Surface clear-sky parallel radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSTINC !TOA (top of atmosph?) incident solar radiation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTCGHG(:)  !Total column greenhouse gases
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTCGRG(:)  !Total column reactive gases
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTCTRAC(:) !Total column tracers
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVD(:)
END TYPE TYPE_SFL_VDIAG

! * Group VX=VCLIX: auxilary climatological diagnostic fields:
TYPE TYPE_SFL_VCLIX
TYPE(TYPE_SURF_MTL_2D),POINTER :: YORO    ! climatological surface geopotential
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTSC    ! climatological surface temperature
TYPE(TYPE_SURF_MTL_2D),POINTER :: YPWS    ! climatological surface max. prop. moisture
TYPE(TYPE_SURF_MTL_2D),POINTER :: YPWP    ! climatological deep soil max. prop. moisture
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSNO    ! climatological snow cover
TYPE(TYPE_SURF_MTL_2D),POINTER :: YTPC    ! climatological deep soil temperature
TYPE(TYPE_SURF_MTL_2D),POINTER :: YSAB    ! climatologic percentage of sand within the soil
TYPE(TYPE_SURF_MTL_2D),POINTER :: YXD2    ! climatologic soil depth
TYPE(TYPE_SURF_MTL_2D),POINTER :: YLSM    ! climatologic land sea mask
TYPE(TYPE_SURF_MTL_2D),POINTER :: YIVEG   ! climatologic type of vegetation
TYPE(TYPE_SURF_MTL_2D),POINTER :: YVX(:)
END TYPE TYPE_SFL_VCLIX

! * Group XA=VEXTRA: extra 3-d diagnostic fields:
TYPE TYPE_SFL_VEXTRA
TYPE(TYPE_SURF_MTL_3D),POINTER :: YXA(:)
END TYPE TYPE_SFL_VEXTRA

! * Group X2=VEXTR2: extra 2-d diagnostic fields:
TYPE TYPE_SFL_VEXTR2
TYPE(TYPE_SURF_MTL_2D),POINTER :: YX2(:)
END TYPE TYPE_SFL_VEXTR2

! End of type definitions

! Data structures

! Prognostic (multi time level) fields

! Soilb
REAL(KIND=JPRB),ALLOCATABLE :: SP_SB (:,:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSP_SBD
TYPE(TYPE_SFL_SOILB)   :: YSP_SB

! Snowg
REAL(KIND=JPRB),ALLOCATABLE :: SP_SG (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSP_SGD
TYPE(TYPE_SFL_SNOWG)   :: YSP_SG

! Resvr
REAL(KIND=JPRB),ALLOCATABLE :: SP_RR (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSP_RRD
TYPE(TYPE_SFL_RESVR)   :: YSP_RR


! Extrp
REAL(KIND=JPRB),ALLOCATABLE :: SP_EP (:,:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSP_EPD
TYPE(TYPE_SFL_EXTRP)   :: YSP_EP 

! Xtrp2
REAL(KIND=JPRB),ALLOCATABLE :: SP_X2 (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSP_X2D
TYPE(TYPE_SFL_XTRP2)   :: YSP_X2

! Canri
REAL(KIND=JPRB),ALLOCATABLE :: SP_CI (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSP_CID
TYPE(TYPE_SFL_CANRI)   :: YSP_CI

! One time level fields

! Varsf
REAL(KIND=JPRB),ALLOCATABLE :: SD_VF (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_VFD
TYPE(TYPE_SFL_VARSF)   :: YSD_VF

! Vclip
REAL(KIND=JPRB),ALLOCATABLE :: SD_VP (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_VPD
TYPE(TYPE_SFL_VCLIP)   :: YSD_VP

! Vcliv
REAL(KIND=JPRB),ALLOCATABLE :: SD_VV (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_VVD
TYPE(TYPE_SFL_VCLIV)   :: YSD_VV

! Vclin
REAL(KIND=JPRB),ALLOCATABLE :: SD_VN (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_VND
TYPE(TYPE_SFL_VCLIN)   :: YSD_VN

! Vclih
REAL(KIND=JPRB),ALLOCATABLE :: SD_VH (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_VHD
TYPE(TYPE_SFL_VCLIH)   :: YSD_VH

! Vclia
REAL(KIND=JPRB),ALLOCATABLE :: SD_VA (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_VAD
TYPE(TYPE_SFL_VCLIA)   :: YSD_VA

! Vo3abc
REAL(KIND=JPRB),ALLOCATABLE :: SD_VC (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_VCD
TYPE(TYPE_SFL_VO3ABC)  :: YSD_VC

! Vdiag
REAL(KIND=JPRB),ALLOCATABLE :: SD_VD (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_VDD
TYPE(TYPE_SFL_VDIAG)   :: YSD_VD

! Waves
REAL(KIND=JPRB),ALLOCATABLE :: SD_WS (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_WSD
TYPE(TYPE_SFL_WAVES)   :: YSD_WS

! Vclix
REAL(KIND=JPRB),ALLOCATABLE :: SD_VX (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_VXD
TYPE(TYPE_SFL_VCLIX)   :: YSD_VX

! Vextra

REAL(KIND=JPRB),ALLOCATABLE :: SD_XA (:,:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_XAD
TYPE(TYPE_SFL_VEXTRA)  :: YSD_XA

! Vextr2

REAL(KIND=JPRB),ALLOCATABLE :: SD_X2 (:,:,:)
TYPE(TYPE_SURF_GEN)    :: YSD_X2D
TYPE(TYPE_SFL_VEXTR2)  :: YSD_X2

!$OMP THREADPRIVATE(ndimsurf,ndimsurfl,nofftraj,nofftraj_cst,nprogsurf)
!$OMP THREADPRIVATE(nprogsurfl,nptrsurf,nstrajgrib,nsurf,nsurfl,ysd_va,ysd_vad)
!$OMP THREADPRIVATE(ysd_vc,ysd_vcd,ysd_vd,ysd_vdd,ysd_vf,ysd_vfd,ysd_vh,ysd_vhd)
!$OMP THREADPRIVATE(ysd_vn,ysd_vnd,ysd_vp,ysd_vpd,ysd_vv,ysd_vvd,ysd_vx,ysd_vxd)
!$OMP THREADPRIVATE(ysd_ws,ysd_wsd,ysd_x2,ysd_x2d,ysd_xa,ysd_xad,ysp_ci,ysp_cid)
!$OMP THREADPRIVATE(ysp_ep,ysp_epd,ysp_rr,ysp_rrd,ysp_sb,ysp_sbd,ysp_sg,ysp_sgd)
!$OMP THREADPRIVATE(ysp_x2,ysp_x2d)

!$OMP THREADPRIVATE(sd_va,sd_vc,sd_vd,sd_vf,sd_vh,sd_vn,sd_vp,sd_vv,sd_vx,sd_ws)
!$OMP THREADPRIVATE(sd_x2,sd_xa,sp_ci,sp_ep,sp_rr,sp_sb,sp_sg,sp_x2,surf_store_array)


!-------------------------------------------------------------------------

CONTAINS

!=========================================================================

SUBROUTINE INI_SFLP3(YDSC,YD,KFLDS,KLEVS,LDMTL,CDGRPNAME)
! Initialize 3-D surface field group
TYPE(TYPE_SURF_GEN),INTENT(INOUT)    :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(INOUT) :: YD(:)
INTEGER(KIND=JPIM),INTENT(IN)        :: KFLDS
INTEGER(KIND=JPIM),INTENT(IN)        :: KLEVS
LOGICAL,INTENT(IN)                   :: LDMTL
CHARACTER(LEN=*),INTENT(IN)          :: CDGRPNAME

INTEGER(KIND=JPIM) :: JFLD, IMAXF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:INI_SFLP3',0,ZHOOK_HANDLE)

IMAXF = SIZE(YD)
YDSC%NUMFLDS = KFLDS
YDSC%NLEVS = KLEVS
YDSC%IPTR  = 1
YDSC%LMTL  = LDMTL
YDSC%CGRPNAME = CDGRPNAME
YDSC%NDIM5 = 0
YDSC%NOFFTRAJ = NOFFTRAJ
YDSC%NOFFTRAJ_CST = NOFFTRAJ_CST

NSURF = NSURF+YDSC%NUMFLDS
NSURFL = NSURFL+YDSC%NUMFLDS*YDSC%NLEVS
IF(LDMTL) THEN
  NPROGSURF = NPROGSURF+YDSC%NUMFLDS
  NPROGSURFL = NPROGSURFL+YDSC%NUMFLDS*YDSC%NLEVS
ENDIF

IF(LDMTL) THEN
  IF (LTWOTL) THEN
    YDSC%NDIM = 2*YDSC%NUMFLDS
  ELSE
    YDSC%NDIM = 3*YDSC%NUMFLDS
  ENDIF
ELSE
  YDSC%NDIM = YDSC%NUMFLDS
ENDIF  
NDIMSURF = NDIMSURF + YDSC%NDIM    
NDIMSURFL = NDIMSURFL + YDSC%NDIM*YDSC%NLEVS
    
DO JFLD=1,KFLDS
  ALLOCATE(YD(JFLD)%IGRBCODE(KLEVS))
  ALLOCATE(YD(JFLD)%CNAME(KLEVS))
  ALLOCATE(YD(JFLD)%REFVALI(KLEVS))
  ALLOCATE(YD(JFLD)%NREQIN(KLEVS))
  YD(JFLD)%IGRBCODE(:) = -999
  YD(JFLD)%CNAME(:) = ''
  YD(JFLD)%REFVALI(:) = 0.0_JPRB
  YD(JFLD)%NREQIN(:) = -1
  YD(JFLD)%MP = JFLD
  IF (YDSC%LMTL) THEN
    YD(JFLD)%MP0 = YD(JFLD)%MP
    IF(LTWOTL) THEN
      YD(JFLD)%MP9 = YD(JFLD)%MP0
      YD(JFLD)%MP1 = YD(JFLD)%MP0+YDSC%NUMFLDS
    ELSE
      YD(JFLD)%MP9 = YD(JFLD)%MP0+YDSC%NUMFLDS
      YD(JFLD)%MP1 = YD(JFLD)%MP0+2*YDSC%NUMFLDS
    ENDIF
  ELSE
    YD(JFLD)%MP0 = NUNDEFLD
    YD(JFLD)%MP9 = NUNDEFLD
    YD(JFLD)%MP1 = NUNDEFLD
  ENDIF
  YD(JFLD)%MP5 = NUNDEFLD
  YD(JFLD)%ITRAJ = 0
ENDDO

DO JFLD=KFLDS+1,IMAXF
  YD(JFLD)%MP  = NUNDEFLD
  YD(JFLD)%MP0 = NUNDEFLD
  YD(JFLD)%MP9 = NUNDEFLD
  YD(JFLD)%MP1 = NUNDEFLD
  YD(JFLD)%MP5 = NUNDEFLD
  YD(JFLD)%ITRAJ = 0
ENDDO

WRITE(NULOUT,*) 'INITIALIZING 3-D SURFACE FIELD GROUP ', YDSC%CGRPNAME
WRITE(NULOUT,*) 'NUMFLDS=',YDSC%NUMFLDS,' NLEVS=',YDSC%NLEVS,' LMTL=',YDSC%LMTL 

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:INI_SFLP3',1,ZHOOK_HANDLE)
END SUBROUTINE INI_SFLP3

!=========================================================================

SUBROUTINE SETUP_SFLP3(YDSC,YD,KGRIB,CDNAME,PDEFAULT,KTRAJ,KREQIN)
! Setup 3-D surface field 
TYPE(TYPE_SURF_GEN),INTENT(INOUT)      :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(INOUT)   :: YD
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KGRIB(:)
CHARACTER(LEN=16) ,OPTIONAL,INTENT(IN) :: CDNAME(:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN) :: PDEFAULT(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KTRAJ
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KREQIN(:)

INTEGER(KIND=JPIM) :: IPTR,JLEV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:SETUP_SFLP3',0,ZHOOK_HANDLE)
IPTR = YDSC%IPTR
IF(IPTR > YDSC%NUMFLDS) THEN
  WRITE(NULERR,*) 'SURFACE FIELDS UNDER-DIMENSINED - GROUP ',&
   & YDSC%CGRPNAME,YDSC%NUMFLDS,KGRIB(1),CDNAME(1)
  CALL ABOR1('IPTR > YDSC%NUMFLDS')
ENDIF
IF(PRESENT(KGRIB)) THEN
  YD%IGRBCODE(:) = KGRIB(:)
ENDIF
IF(PRESENT(KREQIN)) THEN
  YD%NREQIN(:) = KREQIN(:)
ENDIF
IF(PRESENT(CDNAME)) THEN
  YD%CNAME(:)    = CDNAME(:)
ENDIF
IF(PRESENT(PDEFAULT)) THEN
  YD%REFVALI(:) = PDEFAULT
ENDIF
IF(PRESENT(KTRAJ)) THEN
  IF(KTRAJ == 1) THEN
    DO JLEV=1,YDSC%NLEVS
      NSTRAJGRIB(NOFFTRAJ+JLEV) = YD%IGRBCODE(JLEV)
    ENDDO
    NOFFTRAJ = NOFFTRAJ+YDSC%NLEVS
  ELSEIF(KTRAJ == 2) THEN
    NOFFTRAJ_CST = NOFFTRAJ_CST+YDSC%NLEVS
  ELSEIF(KTRAJ /= 0) THEN
    CALL ABOR1('SURFACE_FIELDS:SETUP_SFLP3 - UNKNOWN KTRAJ')
  ENDIF
  YD%ITRAJ = KTRAJ
  YDSC%NDIM5 = YDSC%NDIM5+1
  YD%MP5 = YDSC%NDIM5
ENDIF
DO JLEV=1,YDSC%NLEVS
  IF(YDSC%LMTL) THEN
    WRITE(NULOUT,'(1X,A,2I4,1X,A,6I4)') &
     & YDSC%CGRPNAME(1:6),YDSC%IPTR,JLEV,YD%CNAME(JLEV),YD%IGRBCODE(JLEV),&
     & YD%MP0,YD%MP9,YD%MP1,YD%ITRAJ,YD%NREQIN(JLEV)
  ELSE
    WRITE(NULOUT,'(1X,A,2I4,1X,A,4I4)') &
     & YDSC%CGRPNAME(1:6),YDSC%IPTR,JLEV,YD%CNAME(JLEV),YD%IGRBCODE(JLEV),&
     & YD%MP,YD%ITRAJ,YD%NREQIN(JLEV)
  ENDIF
ENDDO
YDSC%IPTR = YDSC%IPTR+1
  
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:SETUP_SFLP3',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_SFLP3

!=========================================================================

SUBROUTINE INI_SFLP2(YDSC,YD,KFLDS,LDMTL,CDGRPNAME)
! Initialize 2-D surface field group
TYPE(TYPE_SURF_GEN),INTENT(INOUT)    :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(INOUT) :: YD(:)
INTEGER(KIND=JPIM),INTENT(IN)        :: KFLDS
LOGICAL,INTENT(IN)                   :: LDMTL
CHARACTER(LEN=*),INTENT(IN)          :: CDGRPNAME

INTEGER(KIND=JPIM) :: JFLD, IMAXF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:INI_SFLP2',0,ZHOOK_HANDLE)

IMAXF = SIZE(YD)
YDSC%NUMFLDS = KFLDS
YDSC%NLEVS = -1
YDSC%IPTR  = 1
YDSC%LMTL  = LDMTL
YDSC%CGRPNAME = CDGRPNAME
YDSC%NDIM5 = 0
YDSC%NOFFTRAJ = NOFFTRAJ 
YDSC%NOFFTRAJ_CST = NOFFTRAJ_CST 

NSURF = NSURF+YDSC%NUMFLDS
NSURFL = NSURFL+YDSC%NUMFLDS
IF(LDMTL) THEN
  NPROGSURF = NPROGSURF+YDSC%NUMFLDS
  NPROGSURFL = NPROGSURFL+YDSC%NUMFLDS
ENDIF

IF(LDMTL) THEN
  IF (LTWOTL) THEN
    YDSC%NDIM = 2*YDSC%NUMFLDS
  ELSE
    YDSC%NDIM = 3*YDSC%NUMFLDS
  ENDIF
ELSE
  YDSC%NDIM = YDSC%NUMFLDS
ENDIF  
NDIMSURF = NDIMSURF + YDSC%NDIM    
NDIMSURFL = NDIMSURFL + YDSC%NDIM    
DO JFLD=1,KFLDS
  YD(JFLD)%IGRBCODE = -999
  YD(JFLD)%CNAME = ''
  YD(JFLD)%REFVALI = 0.0_JPRB
  YD(JFLD)%NREQIN = -1
  YD(JFLD)%MP = JFLD
  IF (YDSC%LMTL) THEN
    YD(JFLD)%MP0 = YD(JFLD)%MP
    IF(LTWOTL) THEN
      YD(JFLD)%MP9 = YD(JFLD)%MP0
      YD(JFLD)%MP1 = YD(JFLD)%MP0+YDSC%NUMFLDS
    ELSE
      YD(JFLD)%MP9 = YD(JFLD)%MP0+YDSC%NUMFLDS
      YD(JFLD)%MP1 = YD(JFLD)%MP0+2*YDSC%NUMFLDS
    ENDIF
  ELSE
    YD(JFLD)%MP0 = NUNDEFLD
    YD(JFLD)%MP9 = NUNDEFLD
    YD(JFLD)%MP1 = NUNDEFLD
  ENDIF
  YD(JFLD)%MP5 = NUNDEFLD
  YD(JFLD)%ITRAJ = 0
ENDDO
  
DO JFLD=KFLDS+1,IMAXF
  YD(JFLD)%MP  = NUNDEFLD
  YD(JFLD)%MP0 = NUNDEFLD
  YD(JFLD)%MP9 = NUNDEFLD
  YD(JFLD)%MP1 = NUNDEFLD
  YD(JFLD)%MP5 = NUNDEFLD
  YD(JFLD)%ITRAJ = 0
ENDDO

WRITE(NULOUT,*) 'INITIALIZING 2-D SURFACE FIELD GROUP ', YDSC%CGRPNAME
WRITE(NULOUT,*) 'NUMFLDS=',YDSC%NUMFLDS,' LMTL=',YDSC%LMTL 

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:INI_SFLP2',1,ZHOOK_HANDLE)
END SUBROUTINE INI_SFLP2

!=========================================================================

SUBROUTINE SETUP_SFLP2(YDSC,YD,KGRIB,CDNAME,PDEFAULT,KTRAJ,KREQIN)
! Setup 2-D surface field
TYPE(TYPE_SURF_GEN),INTENT(INOUT)      :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(INOUT)   :: YD
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KGRIB
CHARACTER(LEN=16) ,OPTIONAL,INTENT(IN) :: CDNAME
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN) :: PDEFAULT
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KTRAJ
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KREQIN

INTEGER(KIND=JPIM) :: IPTR
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:SETUP_SFLP2',0,ZHOOK_HANDLE)
IPTR = YDSC%IPTR
IF(IPTR > YDSC%NUMFLDS) THEN
  WRITE(NULERR,*) 'SURFACE FIELDS UNDER-DIMENSINED - GROUP ',YDSC%CGRPNAME,YDSC%NUMFLDS,KGRIB,CDNAME
  CALL ABOR1('IPTR > YDSC%NUMFLDS')
ENDIF
IF(PRESENT(KGRIB)) THEN
  YD%IGRBCODE = KGRIB
ENDIF
IF(PRESENT(KREQIN)) THEN
  YD%NREQIN = KREQIN
ENDIF
IF(PRESENT(CDNAME)) THEN
  YD%CNAME    = CDNAME
ENDIF
IF(PRESENT(PDEFAULT)) THEN
  YD%REFVALI = PDEFAULT
ENDIF
IF(PRESENT(KTRAJ)) THEN
  IF(KTRAJ == 1) THEN
    NSTRAJGRIB(NOFFTRAJ+1) = YD%IGRBCODE
    NOFFTRAJ = NOFFTRAJ+1
  ELSEIF(KTRAJ == 2) THEN
    NOFFTRAJ_CST = NOFFTRAJ_CST+1
  ELSEIF(KTRAJ /= 0) THEN
    CALL ABOR1('SURFACE_FIELDS:SETUP_SFLP2 - UNKNOWN KTRAJ')
  ENDIF
  YD%ITRAJ = KTRAJ
  YDSC%NDIM5 = YDSC%NDIM5+1
  YD%MP5 = YDSC%NDIM5
ENDIF
IF(YDSC%LMTL) THEN
  WRITE(NULOUT,'(1X,A,I4,1X,A,6I4)') &
   & YDSC%CGRPNAME(1:6),YDSC%IPTR,YD%CNAME,YD%IGRBCODE,&
   & YD%MP0,YD%MP9,YD%MP1,YD%ITRAJ,YD%NREQIN
ELSE
  WRITE(NULOUT,'(1X,A,I4,1X,A,4I4)') &
   & YDSC%CGRPNAME(1:6),YDSC%IPTR,YD%CNAME,YD%IGRBCODE,YD%MP,YD%ITRAJ,YD%NREQIN
ENDIF
  
YDSC%IPTR = YDSC%IPTR+1
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:SETUP_SFLP2',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_SFLP2

!=========================================================================

SUBROUTINE GPPOPER(CDACT,KBL,PSP_SB,PSP_SG,PSP_RR,PSP_EP,PSP_X2,YDCOM)
! Operations on prognostic surface fields
CHARACTER(LEN=*),INTENT(IN)            :: CDACT
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KBL
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PSP_SB(:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PSP_SG(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PSP_RR(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PSP_EP(:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PSP_X2(:,:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT)   :: YDCOM


REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:GPPOPER',0,ZHOOK_HANDLE)
IF(PRESENT(KBL)) THEN
  CALL GPOPER_3(CDACT,SP_SB(:,:,:,KBL),YSP_SBD,YSP_SB%YSB,YDCOM)
  CALL GPOPER_2(CDACT,SP_SG(:,:,KBL)  ,YSP_SGD,YSP_SG%YSG,YDCOM)
  CALL GPOPER_2(CDACT,SP_RR(:,:,KBL)  ,YSP_RRD,YSP_RR%YRR,YDCOM)
  CALL GPOPER_3(CDACT,SP_EP(:,:,:,KBL),YSP_EPD,YSP_EP%YEP,YDCOM)
  CALL GPOPER_2(CDACT,SP_X2(:,:,KBL)  ,YSP_X2D,YSP_X2%YX2,YDCOM)
ELSE
  IF(PRESENT(PSP_SB)) CALL GPOPER_3(CDACT,PSP_SB(:,:,:),YSP_SBD,YSP_SB%YSB,YDCOM)
  IF(PRESENT(PSP_SG)) CALL GPOPER_2(CDACT,PSP_SG(:,:)  ,YSP_SGD,YSP_SG%YSG,YDCOM)
  IF(PRESENT(PSP_RR)) CALL GPOPER_2(CDACT,PSP_RR(:,:)  ,YSP_RRD,YSP_RR%YRR,YDCOM)
  IF(PRESENT(PSP_EP)) CALL GPOPER_3(CDACT,PSP_EP(:,:,:),YSP_EPD,YSP_EP%YEP,YDCOM)
  IF(PRESENT(PSP_X2)) CALL GPOPER_2(CDACT,PSP_X2(:,:)  ,YSP_X2D,YSP_X2%YX2,YDCOM)
ENDIF
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:GPPOPER',1,ZHOOK_HANDLE)
END SUBROUTINE GPPOPER

!=========================================================================

SUBROUTINE GPOPER(CDACT,KBL,PSP_SB,PSP_SG,PSP_RR,PSD_VF,PSD_VV,YDCOM,PFIELD,PFIELD2)
!Operations on ALL surface groups
CHARACTER(LEN=*),INTENT(IN)            :: CDACT
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KBL
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PSP_SB(:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PSP_SG(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PSP_RR(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PSD_VF(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PSD_VV(:,:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT)   :: YDCOM
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PFIELD(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PFIELD2(:,:)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:GPOPER',0,ZHOOK_HANDLE)
IF(CDACT == 'PUTALLFLDS' .OR. CDACT == 'GETALLFLDS'   .OR.&
 & CDACT == 'TRAJSTORE'  .OR. CDACT == 'TRAJSTORECST' .OR. &
 & CDACT == 'SET0TOTRAJ' .OR. CDACT == 'GETTRAJ'          ) THEN
  IF(.NOT.PRESENT(PFIELD)) CALL ABOR1('SURFACE_FIELDS:GPOPER - PFIELD MISSING')
  IF(SIZE(PFIELD,1) < NPROMA)  CALL ABOR1('SURFACE_FIELDS:GPOPER - SIZE(PFIELD,1) < NPROMA)')
ENDIF
IF(CDACT == 'PUTALLFLDS' .OR. CDACT == 'GETALLFLDS') THEN
  IF(SIZE(PFIELD,2) < NPROGSURFL) CALL ABOR1('SURFACE_FIELDS:GPOPER - SIZE(PFIELD,2) < NPROGSURFL)')
ENDIF
IF(CDACT == 'GETTRAJ') THEN
  IF(.NOT.PRESENT(PFIELD2)) CALL ABOR1('SURFACE_FIELDS:GPOPER - PFIELD2 MISSING')
  IF(SIZE(PFIELD2,1) < NPROMA)  CALL ABOR1('SURFACE_FIELDS:GPOPER - SIZE(PFIELD2,1) < NPROMA)')
ENDIF
IF(PRESENT(YDCOM)) THEN
  YDCOM%L_OK = .FALSE.
  YDCOM%IPTRSURF = 0
  YDCOM%ICOUNT = 0
ENDIF

NPTRSURF = 0
IF(PRESENT(KBL)) THEN
  IF(YSP_SBD%NDIM > 0) THEN
    CALL GPOPER_3(CDACT,SP_SB(:,:,:,KBL),YSP_SBD,YSP_SB%YSB,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSP_SGD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SP_SG(:,:,KBL)  ,YSP_SGD,YSP_SG%YSG,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSP_RRD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SP_RR(:,:,KBL)  ,YSP_RRD,YSP_RR%YRR,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSP_EPD%NDIM > 0) THEN
    CALL GPOPER_3(CDACT,SP_EP(:,:,:,KBL),YSP_EPD,YSP_EP%YEP,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSP_X2D%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SP_X2(:,:,KBL)  ,YSP_X2D,YSP_X2%YX2,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VFD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_VF(:,:,KBL)  ,YSD_VFD,YSD_VF%YVF,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VPD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_VP(:,:,KBL)  ,YSD_VPD,YSD_VP%YVP,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VVD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_VV(:,:,KBL)  ,YSD_VVD,YSD_VV%YVV,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VND%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_VN(:,:,KBL)  ,YSD_VND,YSD_VN%YVN,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VHD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_VH(:,:,KBL)  ,YSD_VHD,YSD_VH%YVH,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VAD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_VA(:,:,KBL)  ,YSD_VAD,YSD_VA%YVA,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VCD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_VC(:,:,KBL)  ,YSD_VCD,YSD_VC%YVC,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VDD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_VD(:,:,KBL)  ,YSD_VDD,YSD_VD%YVD,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_WSD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_WS(:,:,KBL)  ,YSD_WSD,YSD_WS%YWS,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_XAD%NDIM > 0) THEN
    CALL GPOPER_3(CDACT,SD_XA(:,:,:,KBL),YSD_XAD,YSD_XA%YXA,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_X2D%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_X2(:,:,KBL)  ,YSD_X2D,YSD_X2%YX2,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VXD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,SD_VX(:,:,KBL)  ,YSD_VXD,YSD_VX%YVX,YDCOM,PFIELD,PFIELD2)
  ENDIF
ELSE
  IF(YSP_SBD%NDIM > 0) THEN
    IF(PRESENT(PSP_SB)) &
     & CALL GPOPER_3(CDACT,PSP_SB,YSP_SBD,YSP_SB%YSB,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSP_SGD%NDIM > 0) THEN
    IF(PRESENT(PSP_SG)) &
     & CALL GPOPER_2(CDACT,PSP_SG,YSP_SGD,YSP_SG%YSG,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSP_RRD%NDIM > 0) THEN
    IF(PRESENT(PSP_RR)) &
     & CALL GPOPER_2(CDACT,PSP_RR,YSP_RRD,YSP_RR%YRR,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VFD%NDIM > 0) THEN
    IF(PRESENT(PSD_VF)) &
     & CALL GPOPER_2(CDACT,PSD_VF,YSD_VFD,YSD_VF%YVF,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YSD_VVD%NDIM > 0) THEN
    IF(PRESENT(PSD_VV)) &
     & CALL GPOPER_2(CDACT,PSD_VV,YSD_VVD,YSD_VV%YVV,YDCOM,PFIELD,PFIELD2)
  ENDIF
ENDIF
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:GPOPER',1,ZHOOK_HANDLE)
END SUBROUTINE GPOPER

!=========================================================================

SUBROUTINE GPOPER_2(CDACT,PFLD,YDSC,YD,YDCOM,PFIELD,PFIELD2)
! Operations on 2-D surface groups
CHARACTER(LEN=*),INTENT(IN)                :: CDACT
REAL(KIND=JPRB),INTENT(INOUT)              :: PFLD(:,:)
TYPE(TYPE_SURF_GEN),INTENT(IN)             :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(IN)          :: YD(:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT) :: YDCOM
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD2(:,:)

INTEGER(KIND=JPIM) :: J,IPTR,IPTR2
REAL(KIND=JPRB) :: ZZPHY
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:GPOPER_2',0,ZHOOK_HANDLE)
IF(CDACT == 'SET9TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = PFLD(:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO9') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP9)
  ENDDO
ELSEIF(CDACT == 'SET1TO9AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP9)+PFLD(:,YD(J)%MP1)
    PFLD(:,YD(J)%MP1) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'SET0TO1AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP1)+PFLD(:,YD(J)%MP0)
    PFLD(:,YD(J)%MP0) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET9TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = PFLD(:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILT') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = REPSP1*PFLD(:,YD(J)%MP1)+ZZPHY*PFLD(:,YD(J)%MP0)
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILTAD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP1)+PFLD(:,YD(J)%MP0)
    PFLD(:,YD(J)%MP0) = 0.0_JPRB
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP1)+REPSP1*PFLD(:,YD(J)%MP9)
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP0)+ZZPHY *PFLD(:,YD(J)%MP9)
    PFLD(:,YD(J)%MP9) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP0) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET9TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET1TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETALLTOVAL') THEN
  DO J=1,YDSC%NDIM
    PFLD(:,J) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETDEFAULT') THEN
  DO J=1,YDSC%NUMFLDS
    IF(YD(J)%NREQIN == -1) THEN
      PFLD(:,YD(J)%MP) = YD(J)%REFVALI
    ENDIF
  ENDDO
ELSEIF(CDACT == 'TRAJSTORE') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        IPTR = IPTR+1
        PFIELD(:,IPTR) = PFLD(:,YD(J)%MP)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'TRAJSTORECST') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 2) THEN
        IPTR2 = IPTR2+1
        PFIELD(:,IPTR2) = PFLD(:,YD(J)%MP)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'SET0TOTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        IPTR = IPTR+1
        PFLD(:,YD(J)%MP) = PFIELD(:,IPTR)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        IPTR = IPTR+1
        PFLD(:,YD(J)%MP5) = PFIELD(:,IPTR)
      ELSEIF(YD(J)%ITRAJ == 2) THEN
        IPTR2 = IPTR2+1
        PFLD(:,YD(J)%MP5) = PFIELD2(:,IPTR2)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETALLFLDS') THEN
  DO J=1,YDSC%NDIM
    NPTRSURF = NPTRSURF+1
    PFIELD(:,NPTRSURF) = PFLD(:,J)
  ENDDO
ELSEIF(CDACT == 'PUTALLFLDS') THEN
  DO J=1,YDSC%NDIM
    NPTRSURF = NPTRSURF+1
    PFLD(:,J) = PFIELD(:,NPTRSURF)
  ENDDO
ELSEIF(CDACT == 'GETGRIBPOS') THEN
  DO J=1,YDSC%NUMFLDS  
    YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
    IF(YD(J)%IGRBCODE == YDCOM%IGRBCODE) THEN
      YDCOM%IFLDNUM  = YDCOM%IPTRSURF
      YDCOM%L_OK = .TRUE.
    ENDIF
  ENDDO
ELSEIF(CDACT == 'GETFIELD') THEN
  DO J=1,YDSC%NUMFLDS  
    YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
    IF(YDCOM%IPTRSURF == YDCOM%IFLDNUM) THEN
      PFIELD(:,1) = PFLD(:,J)
      YDCOM%L_OK = .TRUE.
    ENDIF
  ENDDO
ELSEIF(CDACT == 'GRIBIN') THEN
  DO J=1,YDSC%NUMFLDS  
    YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
    IF(YD(J)%NREQIN == 1) THEN
      YDCOM%ICOUNT = YDCOM%ICOUNT+1
      YDCOM%ICODES(YDCOM%ICOUNT) = YD(J)%IGRBCODE
    ENDIF
  ENDDO
ELSE
  WRITE(NULOUT,*) 'SURFACE_FIELD:GPPOPER UNKNOWN ACTION - ',CDACT
  CALL ABOR1('SURFACE_FIELD:GPPOPER - UNKNOWN ACTION')
ENDIF
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:GPOPER_2',1,ZHOOK_HANDLE)
END SUBROUTINE GPOPER_2

!=========================================================================

SUBROUTINE GPOPER_3(CDACT,PFLD,YDSC,YD,YDCOM,PFIELD,PFIELD2)
! Operations on 3-D surface groups
CHARACTER(LEN=*),INTENT(IN)                :: CDACT
REAL(KIND=JPRB),INTENT(INOUT)              :: PFLD(:,:,:)
TYPE(TYPE_SURF_GEN),INTENT(IN)             :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(IN)          :: YD(:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT) :: YDCOM
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD2(:,:)

INTEGER(KIND=JPIM) :: J,JLEV,IPTR,IPTR2
REAL(KIND=JPRB) :: ZZPHY
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:GPOPER_3',0,ZHOOK_HANDLE)
IF(CDACT == 'SET9TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = PFLD(:,:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO9') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP9)
  ENDDO
ELSEIF(CDACT == 'SET1TO9AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = PFLD(:,:,YD(J)%MP9)+PFLD(:,:,YD(J)%MP1)
    PFLD(:,:,YD(J)%MP1) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP0) = PFLD(:,:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'SET0TO1AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP1)+PFLD(:,:,YD(J)%MP0)
    PFLD(:,:,YD(J)%MP0) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET9TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = PFLD(:,:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILT') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = REPSP1*PFLD(:,:,YD(J)%MP1)+ZZPHY*PFLD(:,:,YD(J)%MP0)
    PFLD(:,:,YD(J)%MP0) = PFLD(:,:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILTAD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP1)+PFLD(:,:,YD(J)%MP0)
    PFLD(:,:,YD(J)%MP0) = 0.0_JPRB
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP1)+REPSP1*PFLD(:,:,YD(J)%MP9)
    PFLD(:,:,YD(J)%MP0) = PFLD(:,:,YD(J)%MP0)+ZZPHY *PFLD(:,:,YD(J)%MP9)
    PFLD(:,:,YD(J)%MP9) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP0) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET9TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET1TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETALLTOVAL') THEN
  DO J=1,YDSC%NDIM
    PFLD(:,:,J) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETDEFAULT') THEN
  DO J=1,YDSC%NUMFLDS
    DO JLEV=1,YDSC%NLEVS
      IF(YD(J)%NREQIN(JLEV) == -1) THEN
        PFLD(:,JLEV,YD(J)%MP) = YD(J)%REFVALI(JLEV)
      ENDIF
    ENDDO
  ENDDO
ELSEIF(CDACT == 'TRAJSTORE') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR = IPTR+1
          PFIELD(:,IPTR) = PFLD(:,JLEV,YD(J)%MP)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'TRAJSTORECST') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 2) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR2 = IPTR2+1
          PFIELD(:,IPTR2) = PFLD(:,JLEV,YD(J)%MP)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'SET0TOTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR = IPTR+1
          PFLD(:,JLEV,YD(J)%MP) = PFIELD(:,IPTR)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR = IPTR+1
          PFLD(:,JLEV,YD(J)%MP5) = PFIELD(:,IPTR)
        ENDDO
      ELSEIF(YD(J)%ITRAJ == 2) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR2 = IPTR2+1
          PFLD(:,JLEV,YD(J)%MP5) = PFIELD2(:,IPTR2)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETALLFLDS') THEN
  DO J=1,YDSC%NDIM
    DO JLEV=1,YDSC%NLEVS
      NPTRSURF = NPTRSURF+1
      PFIELD(:,NPTRSURF) = PFLD(:,JLEV,J)
    ENDDO
  ENDDO
ELSEIF(CDACT == 'PUTALLFLDS') THEN
  DO J=1,YDSC%NDIM
    DO JLEV=1,YDSC%NLEVS
      NPTRSURF = NPTRSURF+1
      PFLD(:,JLEV,J) = PFIELD(:,NPTRSURF)
    ENDDO
  ENDDO
ELSEIF(CDACT == 'GETGRIBPOS') THEN
  DO J=1,YDSC%NUMFLDS  
    DO JLEV=1,YDSC%NLEVS
      YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
      IF(YD(J)%IGRBCODE(JLEV) == YDCOM%IGRBCODE) THEN
        YDCOM%IFLDNUM  = YDCOM%IPTRSURF
        YDCOM%L_OK = .TRUE.
      ENDIF
    ENDDO
  ENDDO
ELSEIF(CDACT == 'GETFIELD') THEN
  DO J=1,YDSC%NUMFLDS  
    DO JLEV=1,YDSC%NLEVS
      YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
      IF(YDCOM%IPTRSURF == YDCOM%IFLDNUM) THEN
        PFIELD(:,1) = PFLD(:,JLEV,J)
        YDCOM%L_OK = .TRUE.
      ENDIF
    ENDDO
  ENDDO
ELSEIF(CDACT == 'GRIBIN') THEN
  DO J=1,YDSC%NUMFLDS  
    DO JLEV=1,YDSC%NLEVS
      YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
      IF(YD(J)%NREQIN(JLEV) == 1) THEN
        YDCOM%ICOUNT = YDCOM%ICOUNT+1
        YDCOM%ICODES(YDCOM%ICOUNT) = YD(J)%IGRBCODE(JLEV)
      ENDIF
    ENDDO
  ENDDO
ELSE
  WRITE(NULOUT,*) 'SURFACE_FIELD:GPPOPER UNKNOWN ACTION - ',CDACT
  CALL ABOR1('SURFACE_FIELD:GPPOPER - UNKNOWN ACTION')
ENDIF
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:GPOPER_3',1,ZHOOK_HANDLE)
END SUBROUTINE GPOPER_3

!=========================================================================

SUBROUTINE SURF_STORE
! Store all surface fields
INTEGER(KIND=JPIM) :: JBL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:SURF_STORE',0,ZHOOK_HANDLE)
ALLOCATE(SURF_STORE_ARRAY(NPROMA,NDIMSURFL,NGPBLKS))
DO JBL=1,NGPBLKS
  CALL GPOPER('GETALLFLDS',KBL=JBL,PFIELD=SURF_STORE_ARRAY(:,:,JBL))
ENDDO
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:SURF_STORE',1,ZHOOK_HANDLE)
END SUBROUTINE SURF_STORE

!=========================================================================

SUBROUTINE SURF_RESTORE
! Restore all surface fields
INTEGER(KIND=JPIM) :: JBL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:SURF_RESTORE',0,ZHOOK_HANDLE)
IF(.NOT. ALLOCATED(SURF_STORE_ARRAY)) &
 & CALL ABOR1('SURFACE_FIELDS:SURF_RESTORE - SURF_STORE NOT ALLOCATED')
DO JBL=1,NGPBLKS
  CALL GPOPER('PUTALLFLDS',KBL=JBL,PFIELD=SURF_STORE_ARRAY(:,:,JBL))
ENDDO
DEALLOCATE(SURF_STORE_ARRAY)
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:SURF_RESTORE',1,ZHOOK_HANDLE)

END SUBROUTINE SURF_RESTORE

!=========================================================================

SUBROUTINE ALLO_SURF
! Allocate surface field arrays
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:ALLO_SURF',0,ZHOOK_HANDLE)
ALLOCATE(SP_SB(NPROMA,YSP_SBD%NLEVS,YSP_SBD%NDIM,NGPBLKS))
ALLOCATE(SP_SG(NPROMA,YSP_SGD%NDIM,NGPBLKS))
ALLOCATE(SP_RR(NPROMA,YSP_RRD%NDIM,NGPBLKS))
ALLOCATE(SP_EP(NPROMA,YSP_EPD%NLEVS,YSP_EPD%NDIM,NGPBLKS))
ALLOCATE(SP_X2(NPROMA,YSP_X2D%NDIM,NGPBLKS))
ALLOCATE(SP_CI(NPROMA,YSP_CID%NDIM,NGPBLKS))
ALLOCATE(SD_VF(NPROMA,YSD_VFD%NDIM,NGPBLKS))
ALLOCATE(SD_VP(NPROMA,YSD_VPD%NDIM,NGPBLKS))
ALLOCATE(SD_VV(NPROMA,YSD_VVD%NDIM,NGPBLKS))
ALLOCATE(SD_VN(NPROMA,YSD_VND%NDIM,NGPBLKS))
ALLOCATE(SD_VH(NPROMA,YSD_VHD%NDIM,NGPBLKS))
ALLOCATE(SD_VA(NPROMA,YSD_VAD%NDIM,NGPBLKS))
ALLOCATE(SD_VC(NPROMA,YSD_VCD%NDIM,NGPBLKS))
ALLOCATE(SD_VD(NPROMA,YSD_VDD%NDIM,NGPBLKS))
ALLOCATE(SD_WS(NPROMA,YSD_WSD%NDIM,NGPBLKS))
ALLOCATE(SD_XA(NPROMA,YSD_XAD%NLEVS,YSD_XAD%NDIM,NGPBLKS))
ALLOCATE(SD_X2(NPROMA,YSD_X2D%NDIM,NGPBLKS))
ALLOCATE(SD_VX(NPROMA,YSD_VXD%NDIM,NGPBLKS))

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:ALLO_SURF',1,ZHOOK_HANDLE)
END SUBROUTINE ALLO_SURF

!=========================================================================

SUBROUTINE DEALLO_SURF
! Deallocate surface field arrays
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:DEALLO_SURF',0,ZHOOK_HANDLE)
IF(ALLOCATED(SP_SB)) DEALLOCATE(SP_SB)
IF(ALLOCATED(SP_SG)) DEALLOCATE(SP_SG)
IF(ALLOCATED(SP_RR)) DEALLOCATE(SP_RR)
IF(ALLOCATED(SP_EP)) DEALLOCATE(SP_EP)
IF(ALLOCATED(SP_X2)) DEALLOCATE(SP_X2)
IF(ALLOCATED(SP_CI)) DEALLOCATE(SP_CI)
IF(ALLOCATED(SD_VF)) DEALLOCATE(SD_VF)
IF(ALLOCATED(SD_VP)) DEALLOCATE(SD_VP)
IF(ALLOCATED(SD_VV)) DEALLOCATE(SD_VV)
IF(ALLOCATED(SD_VN)) DEALLOCATE(SD_VN)
IF(ALLOCATED(SD_VH)) DEALLOCATE(SD_VH)
IF(ALLOCATED(SD_VA)) DEALLOCATE(SD_VA)
IF(ALLOCATED(SD_VC)) DEALLOCATE(SD_VC)
IF(ALLOCATED(SD_VD)) DEALLOCATE(SD_VD)
IF(ALLOCATED(SD_WS)) DEALLOCATE(SD_WS)
IF(ALLOCATED(SD_XA)) DEALLOCATE(SD_XA)
IF(ALLOCATED(SD_X2)) DEALLOCATE(SD_X2)
IF(ALLOCATED(SD_VX)) DEALLOCATE(SD_VX)
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS:DEALLO_SURF',1,ZHOOK_HANDLE)
END SUBROUTINE DEALLO_SURF

!=========================================================================

END MODULE SURFACE_FIELDS
