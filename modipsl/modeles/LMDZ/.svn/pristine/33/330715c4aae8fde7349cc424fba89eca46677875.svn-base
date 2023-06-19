MODULE YOMGRB

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOM_YGFL  , ONLY : JPGHG, JPGRG , JPTRAC

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    GRIB CODING DESCRIPTORS

! NLOCGRB - ECMWF LOCAL USAGE IDENTIFIER
! NTOTENS - TOTAL NUMBER OF FORCASTS IN ENSEMBLE
! NENSFNB - ENSAMBLE FORECAST NUMBER
! NCALVAL - cal/val number
! NSTREAM - EXPLICIT STREAM NUMBER (OR ZERO FOR DEFAULT STREAMS TO BE USED)
! NSYSTEM - FOR USE IN SEASONAL STREAM (DIFFERENT OPERATIONAL SYSTEMS)
! NMETHOD - FOR USE IN SEASONAL STREAM (DIFFERENT ENSEMBLES)
! NREFERENCE - FOR USE IN HINDCAST STREAM 
! NCONSENSUS - FOR MULTI_ANALYSIS STREAM (1 if consensus)
! NDWD    - FOR MULTI_ANALYSIS STREAM (1 if DWD analysis used)
! NMFR    - FOR MULTI_ANALYSIS STREAM
! NNCEP   - FOR MULTI_ANALYSIS STREAM
! NUKM   - FOR MULTI_ANALYSIS STREAM
! NSTEPLPP(II) - STEP FOR PREVIOUS P.P. OF VARIABLE WITH GRIB CODE II
! NSMAXNP - SUBTRUNCTATION FOR COMPLEX PACKING OF SPECTRAL COEFF
! NBITSSH - NUMBER OF BITS FOR PACKING SPECTRAL COEFF.
! NBITSGG - NUMBER OF BITS FOR PACKING GRID-POINT DATA
! NJDIAG  - SENSITIVITY DIAGNOSTIC NUMBER (1=>J1, 2=>J2, 3=>J3, 4=>J4)
! NJDOMAI - SENSITIVITY DIAGNOSTIC MASK REGION (0=>global, 1=>europe,
!                                               2=>N.H., 3=>S.H.)
! NLEG    - current VAREPS leg number (eg 1(2) for the T399(T255) part of a T399-T255 VAREPS)

! NGRB...  - GRIB codes according to ECMWF local Code Table 2

! NGRBTH   -   3 Potential Temperature
! NGRBALUVP-  15 MODIS albedo UV-vis parallel radiation
! NGRBALUVD-  16 MODIS albedo UV-vis diffuse radiation
! NGRBALNIP-  17 MODIS albedo Near-IR parallel radiation
! NGRBALNID-  18 MODIS albedo Near-IR diffuse radiation
! NGRBSPARC-  20 surface clear-sky PARadiation

! NGRBCVL  -  27 Low vegetation cover
! NGRBCVH  -  28 High vegetation cover
! NGRBTVL  -  29 Low vegetation type
! NGRBTVH  -  30 High vegetation type
! NGRBCI   -  31 Sea ice cover
! NGRBASN  -  32 Snow albedo
! NGRBRSN  -  33 Snow density
! NGRBSST  -  34 Sea surface temperature
! NGRBISTL1-  35 Ice surface temperature layer 1
! NGRBISTL2-  36 Ice surface temperature layer 2
! NGRBISTL3-  37 Ice surface temperature layer 3
! NGRBISTL4-  38 Ice surface temperature layer 4
! NGRBSWL1 -  39 Volumetric soil water content layer 1
! NGRBSWL2 -  40 Volumetric soil water content layer 2
! NGRBSWL3 -  41 Volumetric soil water content layer 3      
! NGRBSWL4 -  42 Volumetric soil water content layer 4 
! NGRBES   -  44 Evaporation of snow
! NGRBSMLT -  45 Snow melt
! NGRB10FG -  49 gust at 10 m level
! NGRBLSPF -  50 large scale precipitation fraction 

! NGRBMONT - 053 Montgomery Geopotential
! NGRBPTHPV- 054 Pressure on Theta and PV surfaces

! NGRBSPAR -  57 surface PARadiation
! NGRBSUVB -  58 surface UV-B radiation
! NGRBCAPE -  59 convect.avail.potential energy
! NGRBPV   -  60 Potential Vorticity

! NGRBSDFOR - 74 standard deviation of a filtered orography
! NGRBSF6  -  210185 SF6 - anthropogenic emissions

! NGRBTCLW -  78 Total column liquid water
! NGRBTCIW -  79 Total column ice water
! NGRBSPD  -  80 !! 80 and 81 extra grib code introduced to 
! NGRBSVD  -  81 !! introduce extra fields for NH

! NGRB082 to NGRB117 reserved for extra fields. Do not use for permanent post-processed fields

! NGRBTCTRAC(JPTRAC)  - 210184 Total column TRAC2: Sulfure hexaFluoride
!                     - 210183 Total column TRAC1: Radon

! Codes for 2D and 3D extra fields
! NGRBMINXTRA  to NGRBMAXXTRA

! NGRBEMIS - 124 Surface Longwave emissivity # has replaced NGRB212

! NGRBAT   - 127 Atmospheric tide
! NGRBBV   - 128 Budget values
! NGRBZ    - 129 Geopotential (at the surface orography)
! NGRBT    - 130 Temperature
! NGRBU    - 131 U-velocity
! NGRBV    - 132 V-velocity
! NGRBQ    - 133 Specific humidity
! NGRBSP   - 134 Surface pressure
! NGRBW    - 135 Vertical velocity
! NGRBTCW  - 136 Total column water
! NGRBTCWV - 137 Total column water vapour
! NGRBVO   - 138 Vorticity (relative)
! NGRBSTL1 - 139 Surface temperature level 1
! NGRBSD   - 141 Snow depth
! NGRBLSP  - 142 Large scale precipitation
! NGRBCP   - 143 Convective precipitation
! NGRBSF   - 144 Snow fall
! NGRBBLD  - 145 Boundary layer dissipation
! NGRBSSHF - 146 Surface sensible heat flux
! NGRBSLHF - 147 Surface latent heat flux
! NGRBCHAR - 148 Charnock parameter 
! NGRB149  - 149 Not used
! NGRB150  - 150 Not used
! NGRBMSL  - 151 Mean sea level pressure
! NGRBLNSP - 152 Log surface pressure
! NGRB153  - 153 Not used
! NGRB154  - 154 Not used
! NGRBD    - 155 Divergence
! NGRBGH   - 156 Height (geopotential)
! NGRBR    - 157 Relative humidity
! NGRBTSP  - 158 Tendency of surface pressure
! NGRBBLH  - 159 Boundary layer height
! NGRBSDOR - 160 Standard deviation of orography
! NGRBISOR - 161 Anisotropy of subgrid scale orography
! NGRBANOR - 162 Angle of subgrid scale orography
! NGRBSLOR - 163 Slope of subgrid scale orography
! NGRBTCC  - 164 Total cloud cover
! NGRB10U  - 165 10 metre u wind
! NGRB10V  - 166 10 metre v wind
! NGRB2T   - 167 2 metre temperature
! NGRB2D   - 168 2 metre dewpoint temperature
! NGRBSSRD - 169 Surface solar radiation downwards
! NGRBSTL2 - 170 Soil temperature level 2
! NGRBLSM  - 172 Land/sea mask
! NGRBSR   - 173 Surface roughness
! NGRBAL   - 174 Albedo
! NGRBSTRD - 175 Surface thermal radiation downwards
! NGRBSSR  - 176 Surface solar radiation
! NGRBSTR  - 177 Surface thermal radiation
! NGRBTSR  - 178 Top solar radiation
! NGRBTTR  - 179 Top thermal radiation
! NGRBEWSS - 180 U-stress
! NGRBNSSS - 181 V-stress
! NGRBE    - 182 Evaporation
! NGRBSTL3 - 183 Soil temperature level 3
! NGRBCCC  - 185 Convective cloud cocer
! NGRBLCC  - 186 Low cloud cover
! NGRBMCC  - 187 Medium cloud cover
! NGRBHCC  - 188 High cloud cover
! NGRBSUND - 189 Sunshine duration
! NGRBEWOV - 190 EW component of sub-grid scale orographic variance
! NGRBNSOV - 191 NS component of sub-grid scale orographic variance
! NGRBNWOV - 192 NWSE component of sub-grid scale orographic variance
! NGRBNEOV - 193 NESW component of sub-grid scale orographic variance
! NGRBTBT  - 194 Brightness temperature (K)
! NGRBLGWS - 195 Latitudinal component of gravity wave stress
! NGRBMGWS - 196 Meridional component of gravity wave stress
! NGRBGWD  - 197 Gravity wave dissipation
! NGRBSRC  - 198 Skin reservoir content
! NGRBVEG  - 199 Percentage of vegetation
! NGRBVSO  - 200 variance of sub-grid scale orogrophy
! NGRBMX2T - 201 Maximum temperature at 2m since last post-processing
! NGRBMN2T - 202 Minimum temperature at 2m since last post-processing
! NGRBO3   - 203 Ozone mixing ratio (EC prognostic ozone)
! NGRBPAW  - 204 Precipitation analysis weights
! NGRBRO   - 205 Runoff
! NGRBTCO3 - 206 Total column ozone 
! NGRB207  - 207 Not used
! NGRBTSRC - 208 Top solar radiation clear sky 
! NGRBTTRC - 209 Top thermal radiation clear sky
! NGRBSSRC - 210 Surface solar radiation clear sky
! NGRBSTRC - 211 Surface thermal radiation clear sky

! NGRBSTINC- 212 TOA incident solar radiation

!-- bunch of codes, confusing in their use ...
! NGRB214  - 214 PROFPROP.RMAX.EA
! NGRB215  - 215 SURFPROP.RMAX.EA
! NGRB216  - 216 RELAPROP.RMAX.EA
! NGRB217  - 217 PROFRESERV.EAU
! NGRB218  - 218 INTSURFTEMPERATU
! NGRB219  - 219 PROFTEMPERATURE

! NGRB221  - 221 Not used
! NGRB222  - 222 Not used
! NGRB223  - 223 Not used
! NGRB224  - 224 used in suafn1
! NGRB225  - 225 used in suafn1
! NGRB226  - 226 Not used
! NGRB227  - 227 used in suafn1
! NGRBTP   - 228 Total precipitation
! NGRBIEWS - 229 Intantaneous X-surface stress
! NGRBINSS - 230 Intantaneous Y-surface stress
! NGRBISHF - 231 Intantaneous surface heat flux
! NGRBIE   - 232 Intantaneous moisture flux (evaporation)
! NGRBLSRH - 234 Logarithm of surface roughness length for heat
! NGRBSKT  - 235 Skin temperature
! NGRBSTL4 - 236 Soil temperature level 4
! NGRBTSN  - 238 Temperature of snow layer
! NGRBCSF  - 239 Convective snow-fall
! NGRBLSF  - 240 Large scale snow-fall
! NGRB241  - 241 Not used
! NGRB242  - 242 Not used
! NGRBFAL  - 243 Forecast albedo
! NGRBFSR  - 244 Forecast surface roughness
! NGRBFLSR - 245 Forecast logarithm of surface roughness for heat
! NGRBCLWC - 246 Cloud liquid water content
! NGRBCIWC - 247 Cloud ice water content
! NGRBCC   - 248 Cloud cover



!-- aerosols in Table 210 -------------------------
! NGRBAERMR01 - 001 aerosol mixing ratio 1
! NGRBAERMR02 - 002 aerosol mixing ratio 2
! NGRBAERMR03 - 003 aerosol mixing ratio 3
! NGRBAERMR04 - 004 aerosol mixing ratio 4
! NGRBAERMR05 - 005 aerosol mixing ratio 5
! NGRBAERMR06 - 006 aerosol mixing ratio 6
! NGRBAERMR07 - 007 aerosol mixing ratio 7
! NGRBAERMR08 - 008 aerosol mixing ratio 8
! NGRBAERMR09 - 009 aerosol mixing ratio 9
! NGRBAERMR10 - 010 aerosol mixing ratio 10
! NGRBAERMR11 - 011 aerosol mixing ratio 11
! NGRBAERMR12 - 012 aerosol mixing ratio 12
! NGRBAERMR13 - 013 aerosol mixing ratio 13
! NGRBAERMR14 - 014 aerosol mixing ratio 14
! NGRBAERMR15 - 015 aerosol mixing ratio 15
! NGRBAERGN01 - 016 aerosol gain acc. 1		2D
! NGRBAERGN02 - 017 aerosol gain acc. 2		2D
! NGRBAERGN03 - 018 aerosol gain acc. 3		2D
! NGRBAERGN04 - 019 aerosol gain acc. 4		2D
! NGRBAERGN05 - 020 aerosol gain acc. 5		2D
! NGRBAERGN06 - 021 aerosol gain acc. 6		2D
! NGRBAERGN07 - 022 aerosol gain acc. 7		2D
! NGRBAERGN08 - 023 aerosol gain acc. 8		2D
! NGRBAERGN09 - 024 aerosol gain acc. 9		2D
! NGRBAERGN10 - 025 aerosol gain acc. 10	2D
! NGRBAERGN11 - 026 aerosol gain acc. 11	2D
! NGRBAERGN12 - 027 aerosol gain acc. 12	2D
! NGRBAERGN13 - 028 aerosol gain acc. 13	2D
! NGRBAERGN14 - 029 aerosol gain acc. 14	2D
! NGRBAERGN15 - 030 aerosol gain acc. 15	2D
! NGRBAERLS01 - 031 aerosol loss acc. 1		2D
! NGRBAERLS02 - 032 aerosol loss acc. 2		2D
! NGRBAERLS03 - 033 aerosol loss acc. 3		2D
! NGRBAERLS04 - 034 aerosol loss acc. 4		2D
! NGRBAERLS05 - 035 aerosol loss acc. 5		2D
! NGRBAERLS06 - 036 aerosol loss acc. 6		2D
! NGRBAERLS07 - 037 aerosol loss acc. 7		2D
! NGRBAERLS08 - 038 aerosol loss acc. 8		2D
! NGRBAERLS09 - 039 aerosol loss acc. 9		2D
! NGRBAERLS10 - 040 aerosol loss acc. 10	2D
! NGRBAERLS11 - 041 aerosol loss acc. 11	2D
! NGRBAERLS12 - 042 aerosol loss acc. 12	2D
! NGRBAERLS13 - 043 aerosol loss acc. 13	2D
! NGRBAERLS14 - 044 aerosol loss acc. 14	2D
! NGRBAERLS15 - 045 aerosol loss acc. 15	2D
! NGRBAERPR   - 046 aerosol precursor mixing ratio 
! NGRBAERSM   - 047 small aerosols mixing ratio 
! NGRBAERLG   - 048 large aerosols mixing ratio 
! NGRBAODPR   - 049 aerosol precursor opt.depth 2D
! NGRBAODSM   - 050 small aerosols opt. depth	2D
! NGRBAODLG   - 051 large aerosols opt. depth	2D
! NGRBAERDEP  - 052 dust emission potential clim2D 
! NGRBAERLTS  - 053 lifting threshold speed clim2D
! NGRBAERSCC  - 054 soli clay content       clim2D

!---------------------------------------------------

! NGRBGHG(JPGHG)   - 210061 GHG1: Carbon dioxide
!                  - 210062 GHG2: Methane
!                  - 210063 GHG3: Nitrous oxide
! NGRBTCGHG(JPGHG) - 210064 Total column GHG1: Carbon Dioxide
!                  - 210065 Total column GHG2: Methane
!                  - 210066 Total column GHG3: Nitrous Oxide

!---------------------------------------------------
! NGRBCO2O         - 210067 CO2 - ocean flux
! NGRBCO2B         - 210068 CO2 - biosphere flux
! NGRBCO2A         - 210069 CO2 - anthropogenic emissions
!
!---------------------------------------------------
! NGRBGRG(JPGRG)   - 210121 GRG1: Nitrogen dioxide
!                  - 210122 GRG2: Sulphur dioxide
!                  - 210123 GRG3: Carbon monoxide
!                  - 210124 GRG4: Formaldehyde

! NGRBTCGRG(JPGRG) - 210125 Total column GRG1: Nitrogen dioxide
!                  - 210126 Total column GRG2: Sulphur dioxide
!                  - 210127 Total column GRG3: Carbon monoxide
!                  - 210128 Total column GRG4: Formaldehyde

!                  - 210203 GRG5: Ozone
!                  - 210206 Total column GRG5: GEMS Ozone

! NGRBTRAC(JPTRAC)  - 210182 TRAC2: SF6
!                   - 210181 TRAC1: Radon


! NGRBSP2, NGRBSP3 - Grib codes for fields in SPA2 and SPA3
! NGRBGP2, NGRBGP3 - Grib codes for grid point fields

INTEGER(KIND=JPIM), PARAMETER :: JPN3SP=10
INTEGER(KIND=JPIM), PARAMETER :: JPN2SP=4

INTEGER(KIND=JPIM) :: NSEC0(2)
INTEGER(KIND=JPIM) :: NSEC1(2048)
INTEGER(KIND=JPIM) :: NSEC2SPP(22)
INTEGER(KIND=JPIM) :: NSEC2SPM(22)
INTEGER(KIND=JPIM) :: NSEC3(2)
INTEGER(KIND=JPIM) :: NSEC4(42)
INTEGER(KIND=JPIM) :: NGRBS3(0:JPN3SP)
INTEGER(KIND=JPIM) :: NGRBS2(0:JPN2SP)
INTEGER(KIND=JPIM) :: NSTEPLPP(255)
REAL(KIND=JPRB) :: RSEC3(2)

INTEGER(KIND=JPIM),ALLOCATABLE:: NSEC2GG(:)
REAL(KIND=JPRB),ALLOCATABLE:: RSEC2(:)

INTEGER(KIND=JPIM) :: MSEC0(2)
INTEGER(KIND=JPIM) :: MSEC1(2048)
INTEGER(KIND=JPIM) :: MSEC2SPP(22)
INTEGER(KIND=JPIM) :: MSEC2SPM(22)
INTEGER(KIND=JPIM) :: MSEC3(2)
INTEGER(KIND=JPIM) :: MSEC4(42)
INTEGER(KIND=JPIM) :: MGRBS3(0:JPN3SP)
INTEGER(KIND=JPIM) :: MGRBS2(0:JPN2SP)
REAL(KIND=JPRB) :: SSEC3(2)

INTEGER(KIND=JPIM),ALLOCATABLE:: MSEC2GG(:)
REAL(KIND=JPRB),ALLOCATABLE:: SSEC2(:)

INTEGER(KIND=JPIM), ALLOCATABLE :: NGRBSP2(:), NGRBSP3(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NGRBGP2(:), NGRBGP3(:)

INTEGER(KIND=JPIM) :: NLOCGRB
INTEGER(KIND=JPIM) :: NTOTENS
INTEGER(KIND=JPIM) :: NENSFNB
INTEGER(KIND=JPIM) :: NCALVAL
INTEGER(KIND=JPIM) :: NSTREAM
INTEGER(KIND=JPIM) :: NSYSTEM
INTEGER(KIND=JPIM) :: NMETHOD
INTEGER(KIND=JPIM) :: NREFERENCE
INTEGER(KIND=JPIM) :: NCONSENSUS
INTEGER(KIND=JPIM) :: NDWD
INTEGER(KIND=JPIM) :: NMFR
INTEGER(KIND=JPIM) :: NNCEP
INTEGER(KIND=JPIM) :: NUKM
INTEGER(KIND=JPIM) :: NSMAXNP
INTEGER(KIND=JPIM) :: NBITSSH
INTEGER(KIND=JPIM) :: NBITSGG
INTEGER(KIND=JPIM) :: NJDIAG
INTEGER(KIND=JPIM) :: NJDOMAI
INTEGER(KIND=JPIM) :: NJITER
INTEGER(KIND=JPIM) :: MLOCGRB
INTEGER(KIND=JPIM) :: MTOTENS
INTEGER(KIND=JPIM) :: MENSFNB
INTEGER(KIND=JPIM) :: MSMAXNP
INTEGER(KIND=JPIM) :: MBITSSH
INTEGER(KIND=JPIM) :: MBITSGG
INTEGER(KIND=JPIM) :: MJDIAG
INTEGER(KIND=JPIM) :: MJDOMAI
INTEGER(KIND=JPIM) :: MJITER

INTEGER(KIND=JPIM) :: NLEG

INTEGER(KIND=JPIM) :: NGRBTH   =  3
INTEGER(KIND=JPIM) :: NGRBALUVP= 15
INTEGER(KIND=JPIM) :: NGRBALUVD= 16
INTEGER(KIND=JPIM) :: NGRBALNIP= 17
INTEGER(KIND=JPIM) :: NGRBALNID= 18
INTEGER(KIND=JPIM) :: NGRBSPARC= 20
INTEGER(KIND=JPIM) :: NGRB21   = 21
INTEGER(KIND=JPIM) :: NGRB22   = 22
INTEGER(KIND=JPIM) :: NGRB23   = 23
INTEGER(KIND=JPIM) :: NGRBCVL  = 27
INTEGER(KIND=JPIM) :: NGRBCVH  = 28
INTEGER(KIND=JPIM) :: NGRBTVL  = 29
INTEGER(KIND=JPIM) :: NGRBTVH  = 30
INTEGER(KIND=JPIM) :: NGRBCI   = 31
INTEGER(KIND=JPIM) :: NGRBASN  = 32
INTEGER(KIND=JPIM) :: NGRBRSN  = 33
INTEGER(KIND=JPIM) :: NGRBSST  = 34
INTEGER(KIND=JPIM) :: NGRBISTL1= 35
INTEGER(KIND=JPIM) :: NGRBISTL2= 36
INTEGER(KIND=JPIM) :: NGRBISTL3= 37
INTEGER(KIND=JPIM) :: NGRBISTL4= 38
INTEGER(KIND=JPIM) :: NGRBSWL1 = 39
INTEGER(KIND=JPIM) :: NGRBSWL2 = 40
INTEGER(KIND=JPIM) :: NGRBSWL3 = 41
INTEGER(KIND=JPIM) :: NGRBSWL4 = 42
INTEGER(KIND=JPIM) :: NGRBES   = 44
INTEGER(KIND=JPIM) :: NGRBSMLT = 45
INTEGER(KIND=JPIM) :: NGRB10FG = 49
INTEGER(KIND=JPIM) :: NGRBLSPF = 50

INTEGER(KIND=JPIM) :: NGRBMONT = 53
INTEGER(KIND=JPIM) :: NGRBPTHPV= 54

INTEGER(KIND=JPIM) :: NGRBSPAR = 57
INTEGER(KIND=JPIM) :: NGRBSUVB = 58
INTEGER(KIND=JPIM) :: NGRBCAPE = 59
INTEGER(KIND=JPIM) :: NGRBPV   = 60

INTEGER(KIND=JPIM) :: NGRBSDFOR= 74
INTEGER(KIND=JPIM) :: NGRBTCLW = 78
INTEGER(KIND=JPIM) :: NGRBTCIW = 79
! LARPEGE
INTEGER(KIND=JPIM) :: NGRBSPD = 80
INTEGER(KIND=JPIM) :: NGRBSVD = 81
! LECMWF
INTEGER(KIND=JPIM) :: NGRB080 = 80
INTEGER(KIND=JPIM) :: NGRB081 = 81

INTEGER(KIND=JPIM) :: NGRBSF6  = 210185
INTEGER(KIND=JPIM) :: NGRBMINXTRA  = 082
INTEGER(KIND=JPIM) :: NGRBMAXXTRA  = 117

INTEGER(KIND=JPIM) :: NGRB082  = 082
INTEGER(KIND=JPIM) :: NGRB083  = 083
INTEGER(KIND=JPIM) :: NGRB084  = 084
INTEGER(KIND=JPIM) :: NGRB085  = 085
INTEGER(KIND=JPIM) :: NGRB086  = 086
INTEGER(KIND=JPIM) :: NGRB087  = 087
INTEGER(KIND=JPIM) :: NGRB088  = 088
INTEGER(KIND=JPIM) :: NGRB089  = 089
INTEGER(KIND=JPIM) :: NGRB090  = 090
INTEGER(KIND=JPIM) :: NGRB091  = 091
INTEGER(KIND=JPIM) :: NGRB092  = 092
INTEGER(KIND=JPIM) :: NGRB093  = 093
INTEGER(KIND=JPIM) :: NGRB094  = 094
INTEGER(KIND=JPIM) :: NGRB095  = 095
INTEGER(KIND=JPIM) :: NGRB096  = 096
INTEGER(KIND=JPIM) :: NGRB097  = 097
INTEGER(KIND=JPIM) :: NGRB098  = 098
INTEGER(KIND=JPIM) :: NGRB099  = 099
INTEGER(KIND=JPIM) :: NGRB100  = 100
INTEGER(KIND=JPIM) :: NGRB101  = 101
INTEGER(KIND=JPIM) :: NGRB102  = 102
INTEGER(KIND=JPIM) :: NGRB103  = 103
INTEGER(KIND=JPIM) :: NGRB104  = 104
INTEGER(KIND=JPIM) :: NGRB105  = 105
INTEGER(KIND=JPIM) :: NGRB106  = 106
INTEGER(KIND=JPIM) :: NGRB107  = 107
INTEGER(KIND=JPIM) :: NGRB108  = 108
INTEGER(KIND=JPIM) :: NGRB109  = 109
INTEGER(KIND=JPIM) :: NGRB110  = 110
INTEGER(KIND=JPIM) :: NGRB111  = 111
INTEGER(KIND=JPIM) :: NGRB112  = 112
INTEGER(KIND=JPIM) :: NGRB113  = 113
INTEGER(KIND=JPIM) :: NGRB114  = 114
INTEGER(KIND=JPIM) :: NGRB115  = 115
INTEGER(KIND=JPIM) :: NGRB116  = 116
INTEGER(KIND=JPIM) :: NGRB117  = 117

INTEGER(KIND=JPIM) :: NGRB118  = 118
INTEGER(KIND=JPIM) :: NGRB119  = 119
INTEGER(KIND=JPIM) :: NGRB120  = 120

INTEGER(KIND=JPIM) :: NGRBEMIS = 124

INTEGER(KIND=JPIM) :: NGRBAT   = 127
INTEGER(KIND=JPIM) :: NGRBBV   = 128
INTEGER(KIND=JPIM) :: NGRBZ    = 129
INTEGER(KIND=JPIM) :: NGRBT    = 130
INTEGER(KIND=JPIM) :: NGRBU    = 131
INTEGER(KIND=JPIM) :: NGRBV    = 132
INTEGER(KIND=JPIM) :: NGRBQ    = 133
INTEGER(KIND=JPIM) :: NGRBSP   = 134
INTEGER(KIND=JPIM) :: NGRBW    = 135
INTEGER(KIND=JPIM) :: NGRBTCW  = 136
INTEGER(KIND=JPIM) :: NGRBTCWV = 137
INTEGER(KIND=JPIM) :: NGRBVO   = 138
INTEGER(KIND=JPIM) :: NGRBSTL1 = 139
INTEGER(KIND=JPIM) :: NGRBSD   = 141
INTEGER(KIND=JPIM) :: NGRBLSP  = 142
INTEGER(KIND=JPIM) :: NGRBCP   = 143
INTEGER(KIND=JPIM) :: NGRBSF   = 144
INTEGER(KIND=JPIM) :: NGRBBLD  = 145
INTEGER(KIND=JPIM) :: NGRBSSHF = 146
INTEGER(KIND=JPIM) :: NGRBSLHF = 147
INTEGER(KIND=JPIM) :: NGRBCHAR = 148
INTEGER(KIND=JPIM) :: NGRB149  = 149
INTEGER(KIND=JPIM) :: NGRB150  = 150
INTEGER(KIND=JPIM) :: NGRBMSL  = 151
INTEGER(KIND=JPIM) :: NGRBLNSP = 152
INTEGER(KIND=JPIM) :: NGRB153  = 153
INTEGER(KIND=JPIM) :: NGRB154  = 154
INTEGER(KIND=JPIM) :: NGRBD    = 155
INTEGER(KIND=JPIM) :: NGRBGH   = 156
INTEGER(KIND=JPIM) :: NGRBR    = 157
INTEGER(KIND=JPIM) :: NGRBTSP  = 158
INTEGER(KIND=JPIM) :: NGRBBLH  = 159
INTEGER(KIND=JPIM) :: NGRBSDOR = 160
INTEGER(KIND=JPIM) :: NGRBISOR = 161
INTEGER(KIND=JPIM) :: NGRBANOR = 162
INTEGER(KIND=JPIM) :: NGRBSLOR = 163
INTEGER(KIND=JPIM) :: NGRBTCC  = 164
INTEGER(KIND=JPIM) :: NGRB10U  = 165
INTEGER(KIND=JPIM) :: NGRB10V  = 166
INTEGER(KIND=JPIM) :: NGRB2T   = 167
INTEGER(KIND=JPIM) :: NGRB2D   = 168
INTEGER(KIND=JPIM) :: NGRBSSRD = 169
INTEGER(KIND=JPIM) :: NGRBSTL2 = 170
INTEGER(KIND=JPIM) :: NGRBLSM  = 172
INTEGER(KIND=JPIM) :: NGRBSR   = 173
INTEGER(KIND=JPIM) :: NGRBAL   = 174
INTEGER(KIND=JPIM) :: NGRBSTRD = 175
INTEGER(KIND=JPIM) :: NGRBSSR  = 176
INTEGER(KIND=JPIM) :: NGRBSTR  = 177
INTEGER(KIND=JPIM) :: NGRBTSR  = 178
INTEGER(KIND=JPIM) :: NGRBTTR  = 179
INTEGER(KIND=JPIM) :: NGRBEWSS = 180
INTEGER(KIND=JPIM) :: NGRBNSSS = 181
INTEGER(KIND=JPIM) :: NGRBE    = 182
INTEGER(KIND=JPIM) :: NGRBSTL3 = 183
INTEGER(KIND=JPIM) :: NGRBCCC  = 185
INTEGER(KIND=JPIM) :: NGRBLCC  = 186
INTEGER(KIND=JPIM) :: NGRBMCC  = 187
INTEGER(KIND=JPIM) :: NGRBHCC  = 188
INTEGER(KIND=JPIM) :: NGRBSUND = 189
INTEGER(KIND=JPIM) :: NGRBEWOV = 190
INTEGER(KIND=JPIM) :: NGRBNSOV = 191
INTEGER(KIND=JPIM) :: NGRBNWOV = 192
INTEGER(KIND=JPIM) :: NGRBNEOV = 193
INTEGER(KIND=JPIM) :: NGRBTBT  = 194
INTEGER(KIND=JPIM) :: NGRBLGWS = 195
INTEGER(KIND=JPIM) :: NGRBMGWS = 196
INTEGER(KIND=JPIM) :: NGRBGWD  = 197
INTEGER(KIND=JPIM) :: NGRBSRC  = 198
INTEGER(KIND=JPIM) :: NGRBVEG  = 199
INTEGER(KIND=JPIM) :: NGRBVSO  = 200
INTEGER(KIND=JPIM) :: NGRBMX2T = 201
INTEGER(KIND=JPIM) :: NGRBMN2T = 202
INTEGER(KIND=JPIM) :: NGRBO3   = 203
INTEGER(KIND=JPIM) :: NGRBPAW  = 204
INTEGER(KIND=JPIM) :: NGRBRO   = 205
INTEGER(KIND=JPIM) :: NGRBTCO3 = 206
INTEGER(KIND=JPIM) :: NGRB207  = 207
INTEGER(KIND=JPIM) :: NGRBTSRC = 208
INTEGER(KIND=JPIM) :: NGRBTTRC = 209
INTEGER(KIND=JPIM) :: NGRBSSRC = 210
INTEGER(KIND=JPIM) :: NGRBSTRC = 211

INTEGER(KIND=JPIM) :: NGRB214  = 214
INTEGER(KIND=JPIM) :: NGRB215  = 215
INTEGER(KIND=JPIM) :: NGRB216  = 216
INTEGER(KIND=JPIM) :: NGRB217  = 217
INTEGER(KIND=JPIM) :: NGRB218  = 218
INTEGER(KIND=JPIM) :: NGRB219  = 219

INTEGER(KIND=JPIM) :: NGRBSTINC= 212
INTEGER(KIND=JPIM) :: NGRBVIMD = 213

INTEGER(KIND=JPIM) :: NGRB222  = 222
INTEGER(KIND=JPIM) :: NGRB223  = 223
INTEGER(KIND=JPIM) :: NGRB224  = 224
INTEGER(KIND=JPIM) :: NGRB225  = 225
INTEGER(KIND=JPIM) :: NGRB226  = 226
INTEGER(KIND=JPIM) :: NGRB227  = 227
INTEGER(KIND=JPIM) :: NGRBTP   = 228
INTEGER(KIND=JPIM) :: NGRBIEWS = 229
INTEGER(KIND=JPIM) :: NGRBINSS = 230
INTEGER(KIND=JPIM) :: NGRBISHF = 231
INTEGER(KIND=JPIM) :: NGRBIE   = 232
INTEGER(KIND=JPIM) :: NGRBLSRH = 234
INTEGER(KIND=JPIM) :: NGRBSKT  = 235
INTEGER(KIND=JPIM) :: NGRBSTL4 = 236
INTEGER(KIND=JPIM) :: NGRBTSN  = 238
INTEGER(KIND=JPIM) :: NGRBCSF  = 239
INTEGER(KIND=JPIM) :: NGRBLSF  = 240
INTEGER(KIND=JPIM) :: NGRB241  = 241
INTEGER(KIND=JPIM) :: NGRB242  = 242
INTEGER(KIND=JPIM) :: NGRBFAL  = 243
INTEGER(KIND=JPIM) :: NGRBFSR  = 244
INTEGER(KIND=JPIM) :: NGRBFLSR = 245
INTEGER(KIND=JPIM) :: NGRBCLWC = 246
INTEGER(KIND=JPIM) :: NGRBCIWC = 247
INTEGER(KIND=JPIM) :: NGRBCC   = 248
INTEGER(KIND=JPIM) :: NGRB249  = 249
INTEGER(KIND=JPIM) :: NGRB250  = 250
INTEGER(KIND=JPIM) :: NGRB251  = 251
INTEGER(KIND=JPIM) :: NGRB252  = 252
INTEGER(KIND=JPIM) :: NGRB253  = 253
INTEGER(KIND=JPIM) :: NGRB254  = 254
INTEGER(KIND=JPIM) :: NGRB255  = 255


!-- aerosols -- Table 210 --------------------------
INTEGER(KIND=JPIM) :: NGRBAERMR01=210001
INTEGER(KIND=JPIM) :: NGRBAERMR02=210002
INTEGER(KIND=JPIM) :: NGRBAERMR03=210003
INTEGER(KIND=JPIM) :: NGRBAERMR04=210004
INTEGER(KIND=JPIM) :: NGRBAERMR05=210005
INTEGER(KIND=JPIM) :: NGRBAERMR06=210006
INTEGER(KIND=JPIM) :: NGRBAERMR07=210007
INTEGER(KIND=JPIM) :: NGRBAERMR08=210008
INTEGER(KIND=JPIM) :: NGRBAERMR09=210009
INTEGER(KIND=JPIM) :: NGRBAERMR10=210010
INTEGER(KIND=JPIM) :: NGRBAERMR11=210011
INTEGER(KIND=JPIM) :: NGRBAERMR12=210012
INTEGER(KIND=JPIM) :: NGRBAERMR13=210013
INTEGER(KIND=JPIM) :: NGRBAERMR14=210014
INTEGER(KIND=JPIM) :: NGRBAERMR15=210015

INTEGER(KIND=JPIM) :: NGRBAERGN01=210016
INTEGER(KIND=JPIM) :: NGRBAERGN02=210017
INTEGER(KIND=JPIM) :: NGRBAERGN03=210018
INTEGER(KIND=JPIM) :: NGRBAERGN04=210019
INTEGER(KIND=JPIM) :: NGRBAERGN05=210020
INTEGER(KIND=JPIM) :: NGRBAERGN06=210021
INTEGER(KIND=JPIM) :: NGRBAERGN07=210022
INTEGER(KIND=JPIM) :: NGRBAERGN08=210023
INTEGER(KIND=JPIM) :: NGRBAERGN09=210024
INTEGER(KIND=JPIM) :: NGRBAERGN10=210025
INTEGER(KIND=JPIM) :: NGRBAERGN11=210026
INTEGER(KIND=JPIM) :: NGRBAERGN12=210027
INTEGER(KIND=JPIM) :: NGRBAERGN13=210028
INTEGER(KIND=JPIM) :: NGRBAERGN14=210029
INTEGER(KIND=JPIM) :: NGRBAERGN15=210030

INTEGER(KIND=JPIM) :: NGRBAERLS01=210031
INTEGER(KIND=JPIM) :: NGRBAERLS02=210032
INTEGER(KIND=JPIM) :: NGRBAERLS03=210033
INTEGER(KIND=JPIM) :: NGRBAERLS04=210034
INTEGER(KIND=JPIM) :: NGRBAERLS05=210035
INTEGER(KIND=JPIM) :: NGRBAERLS06=210036
INTEGER(KIND=JPIM) :: NGRBAERLS07=210037
INTEGER(KIND=JPIM) :: NGRBAERLS08=210038
INTEGER(KIND=JPIM) :: NGRBAERLS09=210039
INTEGER(KIND=JPIM) :: NGRBAERLS10=210040
INTEGER(KIND=JPIM) :: NGRBAERLS11=210041
INTEGER(KIND=JPIM) :: NGRBAERLS12=210042
INTEGER(KIND=JPIM) :: NGRBAERLS13=210043
INTEGER(KIND=JPIM) :: NGRBAERLS14=210044
INTEGER(KIND=JPIM) :: NGRBAERLS15=210045

INTEGER(KIND=JPIM) :: NGRBAERPR  =210046
INTEGER(KIND=JPIM) :: NGRBAERSM  =210047
INTEGER(KIND=JPIM) :: NGRBAERLG  =210048
INTEGER(KIND=JPIM) :: NGRBAODPR  =210049
INTEGER(KIND=JPIM) :: NGRBAODSM  =210050
INTEGER(KIND=JPIM) :: NGRBAODLG  =210051
INTEGER(KIND=JPIM) :: NGRBAERDEP =210052
INTEGER(KIND=JPIM) :: NGRBAERLTS =210053
INTEGER(KIND=JPIM) :: NGRBAERSCC =210054
!---------------------------------------------------
INTEGER(KIND=JPIM) :: NGRBCO2O = 210067
INTEGER(KIND=JPIM) :: NGRBCO2B = 210068
INTEGER(KIND=JPIM) :: NGRBCO2A = 210069

!---------------------------------------------------

INTEGER(KIND=JPIM), DIMENSION(JPGHG) :: NGRBGHG = (/&
 & 210061, 210062, 210063/)
INTEGER(KIND=JPIM), DIMENSION(JPTRAC) :: NGRBTRAC = (/&
 & 210181, 210182/)

INTEGER(KIND=JPIM), DIMENSION(JPGHG) :: NGRBTCGHG = (/&
 & 210064, 210065, 210066/)
INTEGER(KIND=JPIM), DIMENSION(JPTRAC) :: NGRBTCTRAC = (/&
 & 210183, 210184/)

INTEGER(KIND=JPIM), DIMENSION(JPGRG) :: NGRBGRG = (/&
 & 210121, 210122, 210123, 210124, 210203/)

INTEGER(KIND=JPIM), DIMENSION(JPGRG) :: NGRBTCGRG = (/&
 & 210125, 210126, 210127, 210128, 210206/)

!     ------------------------------------------------------------------
!$OMP THREADPRIVATE(mbitsgg,mbitssh,mensfnb,mgrbs2,mgrbs3,mjdiag,mjdomai,mjiter,mlocgrb,msec0,msec1,msec2spm,msec2spp)
!$OMP THREADPRIVATE(msec3,msec4,msmaxnp,mtotens,nbitsgg,nbitssh,ncalval,nconsensus,ndwd,nensfnb,ngrb080,ngrb081)
!$OMP THREADPRIVATE(ngrb082,ngrb083,ngrb084,ngrb085,ngrb086,ngrb087,ngrb088,ngrb089,ngrb090,ngrb091,ngrb092,ngrb093)
!$OMP THREADPRIVATE(ngrb094,ngrb095,ngrb096,ngrb097,ngrb098,ngrb099,ngrb100,ngrb101,ngrb102,ngrb103,ngrb104,ngrb105)
!$OMP THREADPRIVATE(ngrb106,ngrb107,ngrb108,ngrb109,ngrb10fg,ngrb10u,ngrb10v,ngrb110,ngrb111,ngrb112,ngrb113,ngrb114)
!$OMP THREADPRIVATE(ngrb115,ngrb116,ngrb117,ngrb118,ngrb119,ngrb120,ngrb149,ngrb150,ngrb153,ngrb154,ngrb207,ngrb21)
!$OMP THREADPRIVATE(ngrb214,ngrb215,ngrb216,ngrb217,ngrb218,ngrb219,ngrb22,ngrb222,ngrb223,ngrb224,ngrb225,ngrb226)
!$OMP THREADPRIVATE(ngrb227,ngrb23,ngrb241,ngrb242,ngrb249,ngrb250,ngrb251,ngrb252,ngrb253,ngrb254,ngrb255,ngrb2d)
!$OMP THREADPRIVATE(ngrb2t,ngrbaerdep,ngrbaergn01,ngrbaergn02,ngrbaergn03,ngrbaergn04,ngrbaergn05,ngrbaergn06)
!$OMP THREADPRIVATE(ngrbaergn07,ngrbaergn08,ngrbaergn09,ngrbaergn10,ngrbaergn11,ngrbaergn12,ngrbaergn13,ngrbaergn14)
!$OMP THREADPRIVATE(ngrbaergn15,ngrbaerlg,ngrbaerls01,ngrbaerls02,ngrbaerls03,ngrbaerls04,ngrbaerls05,ngrbaerls06)
!$OMP THREADPRIVATE(ngrbaerls07,ngrbaerls08,ngrbaerls09,ngrbaerls10,ngrbaerls11,ngrbaerls12,ngrbaerls13,ngrbaerls14)
!$OMP THREADPRIVATE(ngrbaerls15,ngrbaerlts,ngrbaermr01,ngrbaermr02,ngrbaermr03,ngrbaermr04,ngrbaermr05,ngrbaermr06)
!$OMP THREADPRIVATE(ngrbaermr07,ngrbaermr08,ngrbaermr09,ngrbaermr10,ngrbaermr11,ngrbaermr12,ngrbaermr13,ngrbaermr14)
!$OMP THREADPRIVATE(ngrbaermr15,ngrbaerpr,ngrbaerscc,ngrbaersm,ngrbal,ngrbalnid,ngrbalnip,ngrbaluvd,ngrbaluvp,ngrbanor)
!$OMP THREADPRIVATE(ngrbaodlg,ngrbaodpr,ngrbaodsm,ngrbasn,ngrbat,ngrbbld,ngrbblh,ngrbbv,ngrbcape,ngrbcc,ngrbccc)
!$OMP THREADPRIVATE(ngrbchar,ngrbci,ngrbciwc,ngrbclwc,ngrbco2a,ngrbco2b,ngrbco2o,ngrbcp,ngrbcsf,ngrbcvh,ngrbcvl,ngrbd)
!$OMP THREADPRIVATE(ngrbe,ngrbemis,ngrbes,ngrbewov,ngrbewss,ngrbfal,ngrbflsr,ngrbfsr,ngrbgh,ngrbghg,ngrbgrg,ngrbgwd)
!$OMP THREADPRIVATE(ngrbhcc,ngrbie,ngrbiews,ngrbinss,ngrbishf,ngrbisor,ngrbistl1,ngrbistl2,ngrbistl3,ngrbistl4,ngrblcc)
!$OMP THREADPRIVATE(ngrblgws,ngrblnsp,ngrblsf,ngrblsm,ngrblsp,ngrblspf,ngrblsrh,ngrbmaxxtra,ngrbmcc,ngrbmgws)
!$OMP THREADPRIVATE(ngrbminxtra,ngrbmn2t,ngrbmont,ngrbmsl,ngrbmx2t,ngrbneov,ngrbnsov,ngrbnsss,ngrbnwov,ngrbo3,ngrbpaw)
!$OMP THREADPRIVATE(ngrbpthpv,ngrbpv,ngrbq,ngrbr,ngrbro,ngrbrsn,ngrbs2,ngrbs3,ngrbsd,ngrbsdfor,ngrbsdor,ngrbsf,ngrbsf6)
!$OMP THREADPRIVATE(ngrbskt,ngrbslhf,ngrbslor,ngrbsmlt,ngrbsp,ngrbspar,ngrbsparc,ngrbspd,ngrbsr,ngrbsrc,ngrbsshf,ngrbssr)
!$OMP THREADPRIVATE(ngrbssrc,ngrbssrd,ngrbsst,ngrbstinc,ngrbstl1,ngrbstl2,ngrbstl3,ngrbstl4,ngrbstr,ngrbstrc,ngrbstrd)
!$OMP THREADPRIVATE(ngrbsund,ngrbsuvb,ngrbsvd,ngrbswl1,ngrbswl2,ngrbswl3,ngrbswl4,ngrbt,ngrbtbt,ngrbtcc,ngrbtcghg)
!$OMP THREADPRIVATE(ngrbtcgrg,ngrbtciw,ngrbtclw,ngrbtco3,ngrbtctrac,ngrbtcw,ngrbtcwv,ngrbth,ngrbtp,ngrbtrac,ngrbtsn)
!$OMP THREADPRIVATE(ngrbtsp,ngrbtsr,ngrbtsrc,ngrbttr,ngrbttrc,ngrbtvh,ngrbtvl,ngrbu,ngrbv,ngrbveg,ngrbvimd,ngrbvo)
!$OMP THREADPRIVATE(ngrbvso,ngrbw,ngrbz,njdiag,njdomai,njiter,nleg,nlocgrb,nmethod,nmfr,nncep,nreference,nsec0,nsec1)
!$OMP THREADPRIVATE(nsec2spm,nsec2spp,nsec3,nsec4,nsmaxnp,nsteplpp,nstream,nsystem,ntotens,nukm,rsec3,ssec3)
!$OMP THREADPRIVATE(msec2gg,ngrbgp2,ngrbgp3,ngrbsp2,ngrbsp3,nsec2gg,rsec2,ssec2)
END MODULE YOMGRB
