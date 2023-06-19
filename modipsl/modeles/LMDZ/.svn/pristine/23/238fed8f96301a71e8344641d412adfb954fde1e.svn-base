MODULE YOMFA

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TYPE_FADS, ONLY : FAD

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    GRIB packing options / Options du compactage GRIB

!     NVGRIB  :  Level of GRIB packing (0:no packing, 1:v0 GRIB, 2:v0mod GRIB)
!                Niveau de codage GRIB (0:pas de cod., 1:v0 GRIB, 2:v0mod GRIB)
!     NBITPG  :  Number of bits to code each grid point
!                Nombre de bits pour coder les points de grille
!     NBITCS  :  Number of bits to code each spectral coefficient
!                Nombre de bits pour coder les coefficients spectraux
!     NSTRON  :  Non-packed sub-truncation
!                Niveau de sous-troncature non compactee
!     NPULAP  :  Laplacian power
!                Puissance de Laplacien

INTEGER(KIND=JPIM) :: NVGRIB
INTEGER(KIND=JPIM) :: NBITPG
INTEGER(KIND=JPIM) :: NBITCS
INTEGER(KIND=JPIM) :: NSTRON
INTEGER(KIND=JPIM) :: NPULAP

!     Fields descriptors

TYPE(FAD) :: YFAOROG ! Surface geopotential
TYPE(FAD) :: YFASP   ! Surface pressure
TYPE(FAD) :: YFAPSI  ! Velocity potential
TYPE(FAD) :: YFAKHI  ! Stream function
TYPE(FAD) :: YFAUGEO ! U-Geographical wind
TYPE(FAD) :: YFAVGEO ! V-Geographical wind
TYPE(FAD) :: YFAT    ! Temperature
TYPE(FAD) :: YFAPD   ! Pressure departure (NH)
TYPE(FAD) :: YFAVD   ! Vertical divergence (NH)
TYPE(FAD) :: YFAQ    ! Specific humidity
TYPE(FAD) :: YFAL    ! Liquid water
TYPE(FAD) :: YFAI    ! Ice
TYPE(FAD) :: YFAS    ! Snow
TYPE(FAD) :: YFAR    ! Rain
TYPE(FAD) :: YFAG    ! Graupels
TYPE(FAD) :: YFATKE  ! Turbulent Kinetic Energy
TYPE(FAD) :: YFAO3   ! Ozone mixing ratio
TYPE(FAD) :: YFACLF  ! Cloud fraction
TYPE(FAD) :: YFACPF  ! Convective precipitation flux
TYPE(FAD) :: YFASPF  ! Stratiform precipitation flux
TYPE(FAD) :: YFACVGQ ! CVGQ for French physics
TYPE(FAD) :: YFASDSAT ! Standard deviation of the Saturation Depression
TYPE(FAD) :: YFACVV  ! Convective Vertical Velocity

! Pronostic convection variables
TYPE(FAD) :: YFAUOM  ! Updraught vertic velocity
TYPE(FAD) :: YFAUAL  ! Updraught mesh fraction
TYPE(FAD) :: YFADOM  ! Downdraught vertic velocity
TYPE(FAD) :: YFADAL  ! Downdraught mesh fraction
TYPE(FAD) :: YFAUEN  ! Updraught entrainment
TYPE(FAD) :: YFAUNEBH! Pseudo Hist Conv cloud fraction

! Filtered surface ln Ps for monitoring coupling updates:
TYPE(FAD) :: YFAFSP1  ,YFAFSP2  ,YFAFSP3  ,YFAFSP4  ,YFAFSP5
TYPE(FAD) :: YFASRC  ! Second order flux for Arome Turb. Scheme
TYPE(FAD) :: YFAQVA  ! Total humidity amplitude variation of Q+L+I
!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(nbitcs,nbitpg,npulap,nstron,nvgrib,yfaclf,yfacpf,yfacvgq,yfacvv,yfadal,yfadom,yfafsp1)
!$OMP THREADPRIVATE(yfafsp2,yfafsp3,yfafsp4,yfafsp5,yfag,yfai,yfakhi,yfal,yfao3,yfaorog,yfapd,yfapsi,yfaq)
!$OMP THREADPRIVATE(yfaqva,yfar,yfas,yfasdsat,yfasp,yfaspf,yfasrc,yfat,yfatke,yfaual,yfauen,yfaugeo)
!$OMP THREADPRIVATE(yfaunebh,yfauom,yfavd,yfavgeo)
END MODULE YOMFA
