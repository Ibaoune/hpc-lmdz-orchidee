SUBROUTINE RRTM_RTRN1A_140GP (KLEV,K_ISTART,K_IEND,K_ICLDLYR,P_CLDFRAC,P_TAUCLD,P_ABSS1,&
 & P_OD,P_TAUSF1,P_CLFNET,P_CLHTR,P_FNET,P_HTR,P_TOTDFLUC,P_TOTDFLUX,P_TOTUFLUC,P_TOTUFLUX,&
 & P_TAVEL,PZ,P_TZ,P_TBOUND,PFRAC,P_SEMISS,P_SEMISLW,K_IREFLECT)  

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714
!     Speed-up by D.Salmond, ECMWF, 9907
!     Bug-fix by M.J. Iacono, AER, Inc., 9911
!     Bug-fix by JJMorcrette, ECMWF, 991209 (RAT1, RAT2 initialization)
!     Speed-up by D. Salmond, ECMWF, 9912
!     Bug-fix by JJMorcrette, ECMWF, 0005 (extrapolation T<160K)
!     Speed-up by D. Salmond, ECMWF, 000515

!-* This program calculates the upward fluxes, downward fluxes,
!   and heating rates for an arbitrary atmosphere.  The input to
!   this program is the atmospheric profile and all Planck function
!   information.  First-order "numerical" quadrature is used for the 
!   angle integration, i.e. only one exponential is computed per layer
!   per g-value per band.  Cloud overlap is treated with a generalized
!   maximum/random method in which adjacent cloud layers are treated
!   with maximum overlap, and non-adjacent cloud groups are treated
!   with random overlap.  For adjacent cloud layers, cloud information
!   is carried from the previous two layers.

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARRRTM  , ONLY : JPBAND   ,JPGPT   ,JPLAY
USE YOERRTAB , ONLY : BPADE
USE YOERRTWN , ONLY : TOTPLNK  ,DELWAVE
USE YOERRTFTR, ONLY : NGB

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_ISTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_ICLDLYR(JPLAY) ! Cloud indicator
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_CLDFRAC(JPLAY) ! Cloud fraction
REAL(KIND=JPRB)                  :: Z_CLDFRAC(JPLAY) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUCLD(JPLAY,JPBAND) ! Spectral optical thickness
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ABSS1(JPGPT*JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_OD(JPGPT,JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUSF1(JPGPT*JPLAY) 
REAL(KIND=JPRB)                  :: P_CLFNET(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_CLHTR(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_FNET(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: P_HTR(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TOTDFLUC(0:JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TOTDFLUX(0:JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TOTUFLUC(0:JPLAY) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TOTUFLUX(0:JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAVEL(JPLAY) 
REAL(KIND=JPRB)                  :: PZ(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TZ(0:JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TBOUND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRAC(JPGPT,JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SEMISS(JPBAND) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SEMISLW 
INTEGER(KIND=JPIM)               :: K_IREFLECT ! Argument NOT used
!- from PROFILE             
!- from SP             
!- from SURFACE             
INTEGER(KIND=JPIM) :: INDLAY(JPLAY),INDLEV(0:JPLAY)

REAL(KIND=JPRB) :: Z_BBU1(JPGPT*JPLAY),Z_BBUTOT1(JPGPT*JPLAY)
REAL(KIND=JPRB) :: Z_TLAYFRAC(JPLAY),Z_TLEVFRAC(0:JPLAY)
REAL(KIND=JPRB) :: Z_BGLEV(JPGPT)
!-- DS_000515
REAL(KIND=JPRB) :: Z_PLVL(JPBAND+1,0:JPLAY),Z_PLAY(JPBAND+1,0:JPLAY),Z_WTNUM(3)
!-- DS_000515
REAL(KIND=JPRB) :: Z_ODCLDNW(JPGPT,JPLAY)
REAL(KIND=JPRB) :: Z_SEMIS(JPGPT),Z_RADUEMIT(JPGPT)

REAL(KIND=JPRB) :: Z_RADCLRU1(JPGPT) ,Z_RADCLRD1(JPGPT)
REAL(KIND=JPRB) :: Z_RADLU1(JPGPT)   ,Z_RADLD1(JPGPT)
!-- DS_000515
REAL(KIND=JPRB) :: Z_TRNCLD(JPLAY,JPBAND+1)
!-- DS_000515
REAL(KIND=JPRB) :: Z_ABSCLDNW(JPGPT,JPLAY)
REAL(KIND=JPRB) :: Z_ATOT1(JPGPT*JPLAY)

REAL(KIND=JPRB) :: Z_SURFEMIS(JPBAND),Z_PLNKEMIT(JPBAND)

! dimension of arrays required for cloud overlap calculations

REAL(KIND=JPRB) :: Z_CLRRADU(jpgpt),Z_CLDRADU(jpgpt),Z_OLDCLD(jpgpt)
REAL(KIND=JPRB) :: Z_OLDCLR(jpgpt),Z_RAD(jpgpt),Z_FACCLD1(jplay+1),Z_FACCLD2(jplay+1)
REAL(KIND=JPRB) :: Z_FACCLR1(jplay+1),Z_FACCLR2(jplay+1)
REAL(KIND=JPRB) :: Z_FACCMB1(jplay+1),Z_FACCMB2(jplay+1)
REAL(KIND=JPRB) :: Z_FACCLD1D(0:jplay),Z_FACCLD2D(0:jplay),Z_FACCLR1D(0:jplay)
REAL(KIND=JPRB) :: Z_FACCLR2D(0:jplay),Z_FACCMB1D(0:jplay),Z_FACCMB2D(0:jplay)
REAL(KIND=JPRB) :: Z_CLRRADD(jpgpt),Z_CLDRADD(jpgpt)
INTEGER(KIND=JPIM) :: istcld(jplay+1),istcldd(0:jplay)
!******

!REAL_B :: ZPLVL(JPGPT+1,JPLAY)  ,ZPLAY(JPGPT+1,JPLAY)
!REAL_B :: ZTRNCLD(JPGPT+1,JPLAY),ZTAUCLD(JPGPT+1,JPLAY)

INTEGER(KIND=JPIM) :: IBAND, ICLDDN, IENT, INDBOUND, INDEX, IPR, I_LAY, I_LEV, I_NBI

REAL(KIND=JPRB) :: Z_BBD, Z_BBDTOT, Z_BGLAY, Z_CLDSRC, Z_DBDTLAY, Z_DBDTLEV,&
 & Z_DELBGDN, Z_DELBGUP, Z_DRAD1, Z_DRADCL1, Z_FACTOT1, &
 & Z_FMAX, Z_FMIN, Z_GASSRC, Z_ODSM, Z_PLANKBND, Z_RADCLD, Z_RADD, Z_RADMOD, Z_RAT1, Z_RAT2, Z_SUMPL, &
 & Z_SUMPLEM, Z_TBNDFRAC, Z_TRNS, Z_TTOT, Z_URAD1, Z_URADCL1, ZEXTAU  
REAL(KIND=JPRB) :: ZHOOK_HANDLE



REAL(KIND=JPRB)                  :: CLFNET(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: CLHTR(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: FNET(0:JPLAY) ! Argument NOT used
REAL(KIND=JPRB)                  :: HTR(0:JPLAY) ! Argument NOT used



!--------------------------------------------------------------------------
! Input
!  JPLAY                 ! Maximum number of model layers
!  JPGPT                 ! Total number of g-point subintervals
!  JPBAND                ! Number of longwave spectral bands
!  SECANG                ! Diffusivity angle
!  WTNUM                 ! Weight for radiance to flux conversion
!  KLEV                  ! Number of model layers
!  PAVEL(JPLAY)          ! Mid-layer pressures (hPa)
!  PZ(0:JPLAY)           ! Interface pressures (hPa)
!  TAVEL(JPLAY)          ! Mid-layer temperatures (K)
!  TZ(0:JPLAY)           ! Interface temperatures (K)
!  TBOUND                ! Surface temperature
!  CLDFRAC(JPLAY)        ! Layer cloud fraction
!  TAUCLD(JPLAY,JPBAND)  ! Layer cloud optical thickness
!  ITR
!  PFRAC(JPGPT,JPLAY)    ! Planck function fractions
!  ICLDLYR(JPLAY)        ! Flag for cloudy layers
!  ICLD                  ! Flag for cloudy column
!  IREFLECT              ! Flag for specular reflection
!  SEMISS(JPBAND)        ! Surface spectral emissivity
!  BPADE                 ! Pade constant
!  OD                    ! Clear-sky optical thickness
!  TAUSF1                ! 
!  ABSS1                 !  

!  ABSS(JPGPT*JPLAY)     !
!  ABSCLD(JPLAY)         !
!  ATOT(JPGPT*JPLAY)     !
!  ODCLR(JPGPT,JPLAY)    ! 
!  ODCLD(JPBAND,JPLAY)   !
!  EFCLFR1(JPBAND,JPLAY) ! Effective cloud fraction
!  RADLU(JPGPT)          ! Upward radiance
!  URAD                  ! Spectrally summed upward radiance
!  RADCLRU(JPGPT)        ! Clear-sky upward radiance
!  CLRURAD               ! Spectrally summed clear-sky upward radiance
!  RADLD(JPGPT)          ! Downward radiance
!  DRAD                  ! Spectrally summed downward radiance
!  RADCLRD(JPGPT)        ! Clear-sky downward radiance
!  CLRDRAD               ! Spectrally summed clear-sky downward radiance

! Output
!  TOTUFLUX(0:JPLAY)     ! Upward longwave flux
!  TOTDFLUX(0:JPLAY)     ! Downward longwave flux
!  TOTUFLUC(0:JPLAY)     ! Clear-sky upward longwave flux
!  TOTDFLUC(0:JPLAY)     ! Clear-sky downward longwave flux

! Maximum/Random cloud overlap variables
! for upward radiaitve transfer
!  FACCLR2  fraction of clear radiance from previous layer that needs to 
!           be switched to cloudy stream
!  FACCLR1  fraction of the radiance that had been switched in the previous
!           layer from cloudy to clear that needs to be switched back to
!           cloudy in the current layer
!  FACCLD2  fraction of cloudy radiance from previous layer that needs to 
!           be switched to clear stream
!           be switched to cloudy stream
!  FACCLD1  fraction of the radiance that had been switched in the previous
!           layer from clear to cloudy that needs to be switched back to
!           clear in the current layer
! for downward radiaitve transfer
!  FACCLR2D fraction of clear radiance from previous layer that needs to 
!           be switched to cloudy stream
!  FACCLR1D fraction of the radiance that had been switched in the previous
!           layer from cloudy to clear that needs to be switched back to
!           cloudy in the current layer
!  FACCLD2D fraction of cloudy radiance from previous layer that needs to 
!           be switched to clear stream
!           be switched to cloudy stream
!  FACCLD1D fraction of the radiance that had been switched in the previous
!           layer from clear to cloudy that needs to be switched back to
!           clear in the current layer

!--------------------------------------------------------------------------

! CORRECTION PROVISOIRE BUG POTENTIEL MPLFH
! on initialise le niveau klev+1 de p_cldfrac, tableau surdimensionne
! a 100 mais apparemment non initialise en klev+1
Z_CLDFRAC(1:KLEV)=P_CLDFRAC(1:KLEV)
Z_CLDFRAC(KLEV+1)=0.0_JPRB
IF (LHOOK) CALL DR_HOOK('RRTM_RTRN1A_140GP',0,ZHOOK_HANDLE)
Z_WTNUM(1)=0.5_JPRB
Z_WTNUM(2)=0.0_JPRB
Z_WTNUM(3)=0.0_JPRB

DO I_LAY = 0, KLEV
ENDDO
!-start JJM_000511
IF (P_TBOUND < 339._JPRB .AND. P_TBOUND >= 160._JPRB ) THEN
  INDBOUND = P_TBOUND - 159._JPRB
  Z_TBNDFRAC = P_TBOUND - INT(P_TBOUND)
ELSEIF (P_TBOUND >= 339._JPRB ) THEN
  INDBOUND = 180
  Z_TBNDFRAC = P_TBOUND - 339._JPRB
ELSEIF (P_TBOUND < 160._JPRB ) THEN
  INDBOUND = 1
  Z_TBNDFRAC = P_TBOUND - 160._JPRB
ENDIF  
!-end JJM_000511
  
DO I_LAY = 0, KLEV
  P_TOTUFLUC(I_LAY) = 0.0_JPRB
  P_TOTDFLUC(I_LAY) = 0.0_JPRB
  P_TOTUFLUX(I_LAY) = 0.0_JPRB
  P_TOTDFLUX(I_LAY) = 0.0_JPRB
!-start JJM_000511
  IF (P_TZ(I_LAY) < 339._JPRB .AND. P_TZ(I_LAY) >= 160._JPRB ) THEN
    INDLEV(I_LAY) = P_TZ(I_LAY) - 159._JPRB
    Z_TLEVFRAC(I_LAY) = P_TZ(I_LAY) - INT(P_TZ(I_LAY))
  ELSEIF (P_TZ(I_LAY) >= 339._JPRB ) THEN
    INDLEV(I_LAY) = 180
    Z_TLEVFRAC(I_LAY) = P_TZ(I_LAY) - 339._JPRB
  ELSEIF (P_TZ(I_LAY) < 160._JPRB ) THEN
    INDLEV(I_LAY) = 1
    Z_TLEVFRAC(I_LAY) = P_TZ(I_LAY) - 160._JPRB
  ENDIF    
!-end JJM_000511
ENDDO

!_start_jjm 991209
DO I_LEV=0,KLEV
  Z_FACCLD1(I_LEV+1) = 0.0_JPRB
  Z_FACCLD2(I_LEV+1) = 0.0_JPRB
  Z_FACCLR1(I_LEV+1) = 0.0_JPRB
  Z_FACCLR2(I_LEV+1) = 0.0_JPRB
  Z_FACCMB1(I_LEV+1) = 0.0_JPRB
  Z_FACCMB2(I_LEV+1) = 0.0_JPRB
  Z_FACCLD1D(I_LEV) = 0.0_JPRB
  Z_FACCLD2D(I_LEV) = 0.0_JPRB
  Z_FACCLR1D(I_LEV) = 0.0_JPRB
  Z_FACCLR2D(I_LEV) = 0.0_JPRB
  Z_FACCMB1D(I_LEV) = 0.0_JPRB
  Z_FACCMB2D(I_LEV) = 0.0_JPRB
ENDDO  

Z_RAT1 = 0.0_JPRB
Z_RAT2 = 0.0_JPRB

!_end_jjm 991209

Z_SUMPL   = 0.0_JPRB
Z_SUMPLEM = 0.0_JPRB

ISTCLD(1) = 1
ISTCLDD(KLEV) = 1

DO I_LEV = 1, KLEV
!-- DS_000515
!-start JJM_000511
  IF (P_TAVEL(I_LEV) < 339._JPRB .AND. P_TAVEL(I_LEV) >= 160._JPRB ) THEN
    INDLAY(I_LEV) = P_TAVEL(I_LEV) - 159._JPRB
    Z_TLAYFRAC(I_LEV) = P_TAVEL(I_LEV) - INT(P_TAVEL(I_LEV))
  ELSEIF (P_TAVEL(I_LEV) >= 339._JPRB ) THEN
    INDLAY(I_LEV) = 180
    Z_TLAYFRAC(I_LEV) = P_TAVEL(I_LEV) - 339._JPRB
  ELSEIF (P_TAVEL(I_LEV) < 160._JPRB ) THEN
    INDLAY(I_LEV) = 1
    Z_TLAYFRAC(I_LEV) = P_TAVEL(I_LEV) - 160._JPRB
  ENDIF  
!-end JJM_000511
ENDDO
!-- DS_000515

!-- DS_000515
!OCL SCALAR

DO I_LEV = 1, KLEV
  IF (K_ICLDLYR(I_LEV) == 1) THEN

!mji    
    ISTCLD(I_LEV+1) = 0
    IF (I_LEV  ==  KLEV) THEN
      Z_FACCLD1(I_LEV+1) = 0.0_JPRB
      Z_FACCLD2(I_LEV+1) = 0.0_JPRB
      Z_FACCLR1(I_LEV+1) = 0.0_JPRB
      Z_FACCLR2(I_LEV+1) = 0.0_JPRB
!-- DS_000515      
!SB debug >>
     Z_FACCMB1(I_LEV+1) =0.0_JPRB
     Z_FACCMB2(I_LEV+1) =0.0_JPRB
!SB debug <<
!mji      ISTCLD(LEV+1) = _ZERO_
    ELSEIF (Z_CLDFRAC(I_LEV+1)  >=  Z_CLDFRAC(I_LEV)) THEN
      Z_FACCLD1(I_LEV+1) = 0.0_JPRB
      Z_FACCLD2(I_LEV+1) = 0.0_JPRB
      IF (ISTCLD(I_LEV)  ==  1) THEN
!mji        ISTCLD(LEV+1) = 0
        Z_FACCLR1(I_LEV+1) = 0.0_JPRB
!mji        
        Z_FACCLR2(I_LEV+1) = 0.0_JPRB
        IF (Z_CLDFRAC(I_LEV) < 1.0_JPRB) THEN
          Z_FACCLR2(I_LEV+1) = (Z_CLDFRAC(I_LEV+1)-Z_CLDFRAC(I_LEV))/&
           & (1.0_JPRB-Z_CLDFRAC(I_LEV))  
        ENDIF  
!SB debug >>
      Z_FACCLR2(I_LEV) = 0.0_JPRB
      Z_FACCLD2(I_LEV) = 0.0_JPRB
!SB debug <<
      ELSE
        Z_FMAX = MAX(Z_CLDFRAC(I_LEV),Z_CLDFRAC(I_LEV-1))
!mji
        IF (Z_CLDFRAC(I_LEV+1)  >  Z_FMAX) THEN
          Z_FACCLR1(I_LEV+1) = Z_RAT2
          Z_FACCLR2(I_LEV+1) = (Z_CLDFRAC(I_LEV+1)-Z_FMAX)/(1.0_JPRB-Z_FMAX)
!mji          
        ELSEIF (Z_CLDFRAC(I_LEV+1) < Z_FMAX) THEN
          Z_FACCLR1(I_LEV+1) = (Z_CLDFRAC(I_LEV+1)-Z_CLDFRAC(I_LEV))/&
           & (Z_CLDFRAC(I_LEV-1)-Z_CLDFRAC(I_LEV))  
          Z_FACCLR2(I_LEV+1) = 0.0_JPRB
!mji
        ELSE
          Z_FACCLR1(I_LEV+1) = Z_RAT2  
          Z_FACCLR2(I_LEV+1) = 0.0_JPRB
        ENDIF
      ENDIF
      IF (Z_FACCLR1(I_LEV+1) > 0.0_JPRB .OR. Z_FACCLR2(I_LEV+1) > 0.0_JPRB) THEN
        Z_RAT1 = 1.0_JPRB
        Z_RAT2 = 0.0_JPRB
!SB debug >>
!      ENDIF
      ELSE
        Z_RAT1 = 0.0_JPRB
        Z_RAT2 = 0.0_JPRB
      ENDIF
!SB debug <<
    ELSE
      Z_FACCLR1(I_LEV+1) = 0.0_JPRB
      Z_FACCLR2(I_LEV+1) = 0.0_JPRB
      IF (ISTCLD(I_LEV)  ==  1) THEN
!mji        ISTCLD(LEV+1) = 0
        Z_FACCLD1(I_LEV+1) = 0.0_JPRB
        Z_FACCLD2(I_LEV+1) = (Z_CLDFRAC(I_LEV)-Z_CLDFRAC(I_LEV+1))/Z_CLDFRAC(I_LEV)
!SB debug >>
        Z_FACCLR2(I_LEV) = 0.0_JPRB
        Z_FACCLD2(I_LEV) = 0.0_JPRB
!SB debug <<
      ELSE
        Z_FMIN = MIN(Z_CLDFRAC(I_LEV),Z_CLDFRAC(I_LEV-1))
        IF (Z_CLDFRAC(I_LEV+1)  <=  Z_FMIN) THEN
          Z_FACCLD1(I_LEV+1) = Z_RAT1
          Z_FACCLD2(I_LEV+1) = (Z_FMIN-Z_CLDFRAC(I_LEV+1))/Z_FMIN
        ELSE
          Z_FACCLD1(I_LEV+1) = (Z_CLDFRAC(I_LEV)-Z_CLDFRAC(I_LEV+1))/&
           & (Z_CLDFRAC(I_LEV)-Z_FMIN)  
          Z_FACCLD2(I_LEV+1) = 0.0_JPRB
        ENDIF
      ENDIF
      IF (Z_FACCLD1(I_LEV+1) > 0.0_JPRB .OR. Z_FACCLD2(I_LEV+1) > 0.0_JPRB) THEN
        Z_RAT1 = 0.0_JPRB
        Z_RAT2 = 1.0_JPRB
!SB debug >>
!      ENDIF
      ELSE 
        Z_RAT1 = 0.0_JPRB
        Z_RAT2 = 0.0_JPRB
      ENDIF
!SB debug <<
    ENDIF
!fcc

!SB debug >>
!    IF (I_LEV == 1) THEN
!      Z_FACCMB1(I_LEV+1) = 0.
!      Z_FACCMB2(I_LEV+1) = Z_FACCLD1(I_LEV+1) * Z_FACCLR2(I_LEV)
!    ELSE
!      Z_FACCMB1(I_LEV+1) = Z_FACCLR1(I_LEV+1) * Z_FACCLD2(I_LEV) *Z_CLDFRAC(I_LEV-1)
!      Z_FACCMB2(I_LEV+1) = Z_FACCLD1(I_LEV+1) * Z_FACCLR2(I_LEV) *&
!       & (1.0_JPRB - Z_CLDFRAC(I_LEV-1))   
!    ENDIF
     if(istcld(i_lev).ne.1.and.i_lev.ne.1) then
        z_faccmb1(i_lev+1) = max(0.,min(z_cldfrac(i_lev+1)-z_cldfrac(i_lev), &
               z_cldfrac(i_lev-1)-z_cldfrac(i_lev)))
        z_faccmb2(i_lev+1) = max(0.,min(z_cldfrac(i_lev)-z_cldfrac(i_lev+1), &
               z_cldfrac(i_lev)-z_cldfrac(i_lev-1)))
     endif
!SB debug <<
!end fcc
  ELSE
!-- DS_000515
    ISTCLD(I_LEV+1) = 1
  ENDIF
ENDDO

!_start_jjm 991209
Z_RAT1 = 0.0_JPRB
Z_RAT2 = 0.0_JPRB
!_end_jjm 991209

!-- DS_000515
!OCL SCALAR

DO I_LEV = KLEV, 1, -1
  IF (K_ICLDLYR(I_LEV) == 1) THEN
!mji
    ISTCLDD(I_LEV-1) = 0  
    IF (I_LEV  ==  1) THEN
      Z_FACCLD1D(I_LEV-1) = 0.0_JPRB
      Z_FACCLD2D(I_LEV-1) = 0.0_JPRB
      Z_FACCLR1D(I_LEV-1) = 0.0_JPRB
      Z_FACCLR2D(I_LEV-1) = 0.0_JPRB
      Z_FACCMB1D(I_LEV-1) = 0.0_JPRB
      Z_FACCMB2D(I_LEV-1) = 0.0_JPRB
!mji      ISTCLDD(LEV-1) = _ZERO_
    ELSEIF (Z_CLDFRAC(I_LEV-1)  >=  Z_CLDFRAC(I_LEV)) THEN
      Z_FACCLD1D(I_LEV-1) = 0.0_JPRB
      Z_FACCLD2D(I_LEV-1) = 0.0_JPRB
      IF (ISTCLDD(I_LEV)  ==  1) THEN
!mji        ISTCLDD(LEV-1) = 0
        Z_FACCLR1D(I_LEV-1) = 0.0_JPRB
        Z_FACCLR2D(I_LEV-1) = 0.0_JPRB
        IF (Z_CLDFRAC(I_LEV) < 1.0_JPRB) THEN
          Z_FACCLR2D(I_LEV-1) = (Z_CLDFRAC(I_LEV-1)-Z_CLDFRAC(I_LEV))/&
           & (1.0_JPRB-Z_CLDFRAC(I_LEV))  
        ENDIF
!SB debug >>
       z_facclr2d(i_lev)=0.0_JPRB        
       z_faccld2d(i_lev)=0.0_JPRB
!SB debug <<
      ELSE
        Z_FMAX = MAX(Z_CLDFRAC(I_LEV),Z_CLDFRAC(I_LEV+1))
!mji
        IF (Z_CLDFRAC(I_LEV-1)  >  Z_FMAX) THEN
          Z_FACCLR1D(I_LEV-1) = Z_RAT2
          Z_FACCLR2D(I_LEV-1) = (Z_CLDFRAC(I_LEV-1)-Z_FMAX)/(1.0_JPRB-Z_FMAX)
!mji
        ELSEIF (Z_CLDFRAC(I_LEV-1) < Z_FMAX) THEN
          Z_FACCLR1D(I_LEV-1) = (Z_CLDFRAC(I_LEV-1)-Z_CLDFRAC(I_LEV))/&
           & (Z_CLDFRAC(I_LEV+1)-Z_CLDFRAC(I_LEV))  
          Z_FACCLR2D(I_LEV-1) = 0.0_JPRB
!mji
        ELSE          
          Z_FACCLR1D(I_LEV-1) = Z_RAT2
          Z_FACCLR2D(I_LEV-1) = 0.0_JPRB
        ENDIF
      ENDIF
      IF (Z_FACCLR1D(I_LEV-1) > 0.0_JPRB .OR. Z_FACCLR2D(I_LEV-1) > 0.0_JPRB)THEN
        Z_RAT1 = 1.0_JPRB
        Z_RAT2 = 0.0_JPRB
!SB debug >>
!      ENDIF
      else 
        Z_RAT1 = 0.0_JPRB
        Z_RAT2 = 0.0_JPRB
      endif       
!SB debug <<
    ELSE
      Z_FACCLR1D(I_LEV-1) = 0.0_JPRB
      Z_FACCLR2D(I_LEV-1) = 0.0_JPRB
      IF (ISTCLDD(I_LEV)  ==  1) THEN
!mji        ISTCLDD(LEV-1) = 0
        Z_FACCLD1D(I_LEV-1) = 0.0_JPRB
        Z_FACCLD2D(I_LEV-1) = (Z_CLDFRAC(I_LEV)-Z_CLDFRAC(I_LEV-1))/Z_CLDFRAC(I_LEV)
!SB debug >>
        z_facclr2d(i_lev)=0.0_JPRB
        z_faccld2d(i_lev)=0.0_JPRB
!SB debug <<
      ELSE
        Z_FMIN = MIN(Z_CLDFRAC(I_LEV),Z_CLDFRAC(I_LEV+1))
        IF (Z_CLDFRAC(I_LEV-1)  <=  Z_FMIN) THEN
          Z_FACCLD1D(I_LEV-1) = Z_RAT1
          Z_FACCLD2D(I_LEV-1) = (Z_FMIN-Z_CLDFRAC(I_LEV-1))/Z_FMIN
        ELSE
          Z_FACCLD1D(I_LEV-1) = (Z_CLDFRAC(I_LEV)-Z_CLDFRAC(I_LEV-1))/&
           & (Z_CLDFRAC(I_LEV)-Z_FMIN)  
          Z_FACCLD2D(I_LEV-1) = 0.0_JPRB
        ENDIF
      ENDIF
      IF (Z_FACCLD1D(I_LEV-1) > 0.0_JPRB .OR. Z_FACCLD2D(I_LEV-1) > 0.0_JPRB)THEN
        Z_RAT1 = 0.0_JPRB
        Z_RAT2 = 1.0_JPRB
!SB debug >>
!      ENDIF
      ELSE
        Z_RAT1 = 0.0_JPRB
        Z_RAT2 = 0.0_JPRB
      ENDIF
!SB debug <<
    ENDIF
!SB debug >>
!    Z_FACCMB1D(I_LEV-1) = Z_FACCLR1D(I_LEV-1) * Z_FACCLD2D(I_LEV) *Z_CLDFRAC(I_LEV+1)
!    Z_FACCMB2D(I_LEV-1) = Z_FACCLD1D(I_LEV-1) * Z_FACCLR2D(I_LEV) *&
!     & (1.0_JPRB - Z_CLDFRAC(I_LEV+1))  
    if (istcldd(i_lev).ne.1.and.i_lev.ne.1) then
       z_faccmb1d(i_lev-1) = max(0.,min(z_cldfrac(i_lev+1)-z_cldfrac(i_lev), &
                            z_cldfrac(i_lev-1)-z_cldfrac(i_lev)))
       z_faccmb2d(i_lev-1) = max(0.,min(z_cldfrac(i_lev)-z_cldfrac(i_lev+1), &
                    z_cldfrac(i_lev)-z_cldfrac(i_lev-1)))
    endif
!SB debug <<
  ELSE
    ISTCLDD(I_LEV-1) = 1
  ENDIF
ENDDO

!- Loop over frequency bands.

DO IBAND = K_ISTART, K_IEND
  Z_DBDTLEV = TOTPLNK(INDBOUND+1,IBAND)-TOTPLNK(INDBOUND,IBAND)
  Z_PLANKBND = DELWAVE(IBAND) * (TOTPLNK(INDBOUND,IBAND) + Z_TBNDFRAC * Z_DBDTLEV)
  Z_DBDTLEV = TOTPLNK(INDLEV(0)+1,IBAND) -TOTPLNK(INDLEV(0),IBAND)
!-- DS_000515
  Z_PLVL(IBAND,0) = DELWAVE(IBAND)&
   & * (TOTPLNK(INDLEV(0),IBAND) + Z_TLEVFRAC(0)*Z_DBDTLEV)  

  Z_SURFEMIS(IBAND) = P_SEMISS(IBAND)
  Z_PLNKEMIT(IBAND) = Z_SURFEMIS(IBAND) * Z_PLANKBND
  Z_SUMPLEM  = Z_SUMPLEM + Z_PLNKEMIT(IBAND)
  Z_SUMPL    = Z_SUMPL   + Z_PLANKBND
!--DS
ENDDO
!---

!-- DS_000515
DO I_LEV = 1, KLEV
  DO IBAND = K_ISTART, K_IEND
! print *,'RTRN1A: I_LEV JPLAY IBAND INDLAY',I_LEV,JPLAY,IBAND,INDLAY(I_LEV)
!----              
!- Calculate the integrated Planck functions for at the
!  level and layer temperatures.
!  Compute cloud transmittance for cloudy layers.
    Z_DBDTLEV = TOTPLNK(INDLEV(I_LEV)+1,IBAND) - TOTPLNK(INDLEV(I_LEV),IBAND)
    Z_DBDTLAY = TOTPLNK(INDLAY(I_LEV)+1,IBAND) - TOTPLNK(INDLAY(I_LEV),IBAND)
!-- DS_000515
    Z_PLAY(IBAND,I_LEV) = DELWAVE(IBAND)&
     & *(TOTPLNK(INDLAY(I_LEV),IBAND)+Z_TLAYFRAC(I_LEV)*Z_DBDTLAY)  
    Z_PLVL(IBAND,I_LEV) = DELWAVE(IBAND)&
     & *(TOTPLNK(INDLEV(I_LEV),IBAND)+Z_TLEVFRAC(I_LEV)*Z_DBDTLEV)  
    IF (K_ICLDLYR(I_LEV) > 0) THEN
      ZEXTAU = MIN( P_TAUCLD(I_LEV,IBAND), 200._JPRB)
      Z_TRNCLD(I_LEV,IBAND) = EXP( -ZEXTAU )
    ENDIF
!-- DS_000515
  ENDDO

ENDDO

P_SEMISLW = Z_SUMPLEM / Z_SUMPL

!--DS
!O IPR = 1, JPGPT
! NBI = NGB(IPR)
! DO LEV =  1 , KLEV
!-- DS_000515
!   ZPLAY(IPR,LEV) = PLAY(LEV,NGB(IPR))
!   ZPLVL(IPR,LEV) = PLVL(LEV-1,NGB(IPR))
!   ZTAUCLD(IPR,LEV) = TAUCLD(LEV,NGB(IPR))
!   ZTRNCLD(IPR,LEV) = TRNCLD(LEV,NGB(IPR))
!-- DS_000515
! ENDDO
!NDDO
!----      

!- For cloudy layers, set cloud parameters for radiative transfer.
DO I_LEV = 1, KLEV
  IF (K_ICLDLYR(I_LEV) > 0) THEN
    DO IPR = 1, JPGPT
!--DS          
!            NBI = NGB(IPR)
      Z_ODCLDNW(IPR,I_LEV) = P_TAUCLD(I_LEV,NGB(IPR))
      Z_ABSCLDNW(IPR,I_LEV) = 1.0_JPRB - Z_TRNCLD(I_LEV,NGB(IPR))
!----            
!            EFCLFRNW(IPR,LEV) = ABSCLDNW(IPR,LEV) * CLDFRAC(LEV)
    ENDDO
  ENDIF
ENDDO

!- Initialize for radiative transfer.
DO IPR = 1, JPGPT
  Z_RADCLRD1(IPR) = 0.0_JPRB
  Z_RADLD1(IPR)   = 0.0_JPRB
  I_NBI = NGB(IPR)
  Z_SEMIS(IPR) = Z_SURFEMIS(I_NBI)
  Z_RADUEMIT(IPR) = PFRAC(IPR,1) * Z_PLNKEMIT(I_NBI)
!-- DS_000515
  Z_BGLEV(IPR) = PFRAC(IPR,KLEV) * Z_PLVL(I_NBI,KLEV)
ENDDO

!- Downward radiative transfer.
!  *** DRAD1 holds summed radiance for total sky stream
!  *** DRADCL1 holds summed radiance for clear sky stream

ICLDDN = 0
DO I_LEV = KLEV, 1, -1
  Z_DRAD1   = 0.0_JPRB
  Z_DRADCL1 = 0.0_JPRB

  IF (K_ICLDLYR(I_LEV) == 1) THEN

!  *** Cloudy layer
    ICLDDN = 1
    IENT = JPGPT * (I_LEV-1)
    DO IPR = 1, JPGPT
      INDEX = IENT + IPR
!--DS            
!            NBI = NGB(IPR)
      Z_BGLAY = PFRAC(IPR,I_LEV) * Z_PLAY(NGB(IPR),I_LEV)
!----            
      Z_DELBGUP     = Z_BGLEV(IPR) - Z_BGLAY
      Z_BBU1(INDEX) = Z_BGLAY + P_TAUSF1(INDEX) * Z_DELBGUP
!--DS            
      Z_BGLEV(IPR) = PFRAC(IPR,I_LEV) * Z_PLVL(NGB(IPR),I_LEV-1)
!----            
      Z_DELBGDN = Z_BGLEV(IPR) - Z_BGLAY
      Z_BBD = Z_BGLAY + P_TAUSF1(INDEX) * Z_DELBGDN
!- total-sky downward flux          
      Z_ODSM = P_OD(IPR,I_LEV) + Z_ODCLDNW(IPR,I_LEV)
      Z_FACTOT1 = Z_ODSM / (BPADE + Z_ODSM)
      Z_BBUTOT1(INDEX) = Z_BGLAY + Z_FACTOT1 * Z_DELBGUP
      Z_ATOT1(INDEX) = P_ABSS1(INDEX) + Z_ABSCLDNW(IPR,I_LEV)&
       & - P_ABSS1(INDEX) * Z_ABSCLDNW(IPR,I_LEV)  
      Z_BBDTOT = Z_BGLAY + Z_FACTOT1 * Z_DELBGDN
      Z_GASSRC = Z_BBD * P_ABSS1(INDEX)
!***
      IF (ISTCLDD(I_LEV)  ==  1) THEN
        Z_CLDRADD(IPR) = Z_CLDFRAC(I_LEV) * Z_RADLD1(IPR)
        Z_CLRRADD(IPR) = Z_RADLD1(IPR) - Z_CLDRADD(IPR)
        Z_OLDCLD(IPR) = Z_CLDRADD(IPR)
        Z_OLDCLR(IPR) = Z_CLRRADD(IPR)
        Z_RAD(IPR) = 0.0_JPRB
      ENDIF
      Z_TTOT = 1.0_JPRB - Z_ATOT1(INDEX)
      Z_CLDSRC = Z_BBDTOT * Z_ATOT1(INDEX)
      
! Separate RT equations for clear and cloudy streams      
      Z_CLDRADD(IPR) = Z_CLDRADD(IPR) * Z_TTOT + Z_CLDFRAC(I_LEV) * Z_CLDSRC
      Z_CLRRADD(IPR) = Z_CLRRADD(IPR) * (1.0_JPRB-P_ABSS1(INDEX)) +&
       & (1.0_JPRB - Z_CLDFRAC(I_LEV)) * Z_GASSRC  

!  Total sky downward radiance
      Z_RADLD1(IPR) = Z_CLDRADD(IPR) + Z_CLRRADD(IPR)
      Z_DRAD1 = Z_DRAD1 + Z_RADLD1(IPR)
      
!  Clear-sky downward radiance          
      Z_RADCLRD1(IPR) = Z_RADCLRD1(IPR)+(Z_BBD-Z_RADCLRD1(IPR))*P_ABSS1(INDEX)
      Z_DRADCL1 = Z_DRADCL1 + Z_RADCLRD1(IPR)

!* Code to account for maximum/random overlap:
!   Performs RT on the radiance most recently switched between clear and
!   cloudy streams
      Z_RADMOD = Z_RAD(IPR) * (Z_FACCLR1D(I_LEV-1) * (1.0_JPRB-P_ABSS1(INDEX)) +&
       & Z_FACCLD1D(I_LEV-1) *  Z_TTOT) - &
       & Z_FACCMB1D(I_LEV-1) * Z_GASSRC + &
       & Z_FACCMB2D(I_LEV-1) * Z_CLDSRC  
       
!   Computes what the clear and cloudy streams would have been had no
!   radiance been switched       
      Z_OLDCLD(IPR) = Z_CLDRADD(IPR) - Z_RADMOD
      Z_OLDCLR(IPR) = Z_CLRRADD(IPR) + Z_RADMOD
      
!   Computes the radiance to be switched between clear and cloudy.      
      Z_RAD(IPR) = -Z_RADMOD + Z_FACCLR2D(I_LEV-1)*Z_OLDCLR(IPR) -&
       & Z_FACCLD2D(I_LEV-1)*Z_OLDCLD(IPR)  
      Z_CLDRADD(IPR) = Z_CLDRADD(IPR) + Z_RAD(IPR)
      Z_CLRRADD(IPR) = Z_CLRRADD(IPR) - Z_RAD(IPR)
!***

    ENDDO

  ELSE

!  *** Clear layer
!  *** DRAD1 holds summed radiance for total sky stream
!  *** DRADCL1 holds summed radiance for clear sky stream

    IENT = JPGPT * (I_LEV-1)
    IF (ICLDDN == 1) THEN
      DO IPR = 1, JPGPT
        INDEX = IENT + IPR
!--DS         
!           NBI = NGB(IPR)
        Z_BGLAY = PFRAC(IPR,I_LEV) * Z_PLAY(NGB(IPR),I_LEV)
!----            
        Z_DELBGUP     = Z_BGLEV(IPR) - Z_BGLAY
        Z_BBU1(INDEX) = Z_BGLAY + P_TAUSF1(INDEX) * Z_DELBGUP
!--DS            
        Z_BGLEV(IPR) = PFRAC(IPR,I_LEV) * Z_PLVL(NGB(IPR),I_LEV-1)
!----                      
        Z_DELBGDN = Z_BGLEV(IPR) - Z_BGLAY
        Z_BBD = Z_BGLAY + P_TAUSF1(INDEX) * Z_DELBGDN
        
!- total-sky downward radiance
        Z_RADLD1(IPR) = Z_RADLD1(IPR)+(Z_BBD-Z_RADLD1(IPR))*P_ABSS1(INDEX)
        Z_DRAD1 = Z_DRAD1 + Z_RADLD1(IPR)
        
!- clear-sky downward radiance
!-  Set clear sky stream to total sky stream as long as layers
!-  remain clear.  Streams diverge when a cloud is reached.
        Z_RADCLRD1(IPR) = Z_RADCLRD1(IPR)+(Z_BBD-Z_RADCLRD1(IPR))*P_ABSS1(INDEX)
        Z_DRADCL1 = Z_DRADCL1 + Z_RADCLRD1(IPR)
      ENDDO
            
    ELSE
        
      DO IPR = 1, JPGPT
        INDEX = IENT + IPR
!--DS         
!           NBI = NGB(IPR)
        Z_BGLAY = PFRAC(IPR,I_LEV) * Z_PLAY(NGB(IPR),I_LEV)
!----            
        Z_DELBGUP     = Z_BGLEV(IPR) - Z_BGLAY
        Z_BBU1(INDEX) = Z_BGLAY + P_TAUSF1(INDEX) * Z_DELBGUP
!--DS            
        Z_BGLEV(IPR) = PFRAC(IPR,I_LEV) * Z_PLVL(NGB(IPR),I_LEV-1)
!----                      
        Z_DELBGDN = Z_BGLEV(IPR) - Z_BGLAY
        Z_BBD = Z_BGLAY + P_TAUSF1(INDEX) * Z_DELBGDN
!- total-sky downward flux          
        Z_RADLD1(IPR) = Z_RADLD1(IPR)+(Z_BBD-Z_RADLD1(IPR))*P_ABSS1(INDEX)
        Z_DRAD1 = Z_DRAD1 + Z_RADLD1(IPR)
!- clear-sky downward flux          
!-  Set clear sky stream to total sky stream as long as layers
!-  remain clear.  Streams diverge when a cloud is reached.
        Z_RADCLRD1(IPR) = Z_RADLD1(IPR)
      ENDDO
      Z_DRADCL1 = Z_DRAD1
    ENDIF
    
  ENDIF

  P_TOTDFLUC(I_LEV-1) = Z_DRADCL1 * Z_WTNUM(1)
  P_TOTDFLUX(I_LEV-1) = Z_DRAD1   * Z_WTNUM(1)

ENDDO

! Spectral reflectivity and reflectance
! Includes the contribution of spectrally varying longwave emissivity 
! and reflection from the surface to the upward radiative transfer.
! Note: Spectral and Lambertian reflections are identical for the one
! angle flux integration used here.

Z_URAD1   = 0.0_JPRB
Z_URADCL1 = 0.0_JPRB

!start JJM_000511
!IF (IREFLECT  ==  0) THEN
!- Lambertian reflection.
DO IPR = 1, JPGPT
! Clear-sky radiance
!    RADCLD = _TWO_ * (RADCLRD1(IPR) * WTNUM(1) )
  Z_RADCLD = Z_RADCLRD1(IPR)
  Z_RADCLRU1(IPR) = Z_RADUEMIT(IPR) + (1.0_JPRB - Z_SEMIS(IPR)) * Z_RADCLD
  Z_URADCL1 = Z_URADCL1 + Z_RADCLRU1(IPR)

! Total sky radiance
!    RADD = _TWO_ * (RADLD1(IPR) * WTNUM(1) )
  Z_RADD = Z_RADLD1(IPR)
  Z_RADLU1(IPR) = Z_RADUEMIT(IPR) + (1.0_JPRB - Z_SEMIS(IPR)) * Z_RADD
  Z_URAD1 = Z_URAD1 + Z_RADLU1(IPR)
ENDDO
P_TOTUFLUC(0) = Z_URADCL1 * 0.5_JPRB
P_TOTUFLUX(0) = Z_URAD1 * 0.5_JPRB
!ELSE
!!- Specular reflection.
!  DO IPR = 1, JPGPT
!    RADCLU = RADUEMIT(IPR)
!    RADCLRU1(IPR) = RADCLU + (_ONE_ - SEMIS(IPR)) * RADCLRD1(IPR)
!    URADCL1 = URADCL1 + RADCLRU1(IPR)

!    RADU = RADUEMIT(IPR)
!    RADLU1(IPR) = RADU + (_ONE_ - SEMIS(IPR)) * RADLD1(IPR)
!    URAD1 = URAD1 + RADLU1(IPR)
!  ENDDO
!  TOTUFLUC(0) = URADCL1 * WTNUM(1)
!  TOTUFLUX(0) = URAD1   * WTNUM(1)
!ENDIF

!- Upward radiative transfer.
!- *** URAD1 holds the summed radiance for total sky stream
!- *** URADCL1 holds the summed radiance for clear sky stream
DO I_LEV = 1, KLEV
  Z_URAD1   = 0.0_JPRB
  Z_URADCL1 = 0.0_JPRB

! Check flag for cloud in current layer
  IF (K_ICLDLYR(I_LEV) == 1) THEN

!- *** Cloudy layer
    IENT = JPGPT * (I_LEV-1)
    DO IPR = 1, JPGPT
      INDEX = IENT + IPR
!- total-sky upward flux          
      Z_GASSRC = Z_BBU1(INDEX) * P_ABSS1(INDEX)

!- If first cloudy layer in sequence, split up radiance into clear and
!    cloudy streams depending on cloud fraction
      IF (ISTCLD(I_LEV)  ==  1) THEN
        Z_CLDRADU(IPR) = Z_CLDFRAC(I_LEV) * Z_RADLU1(IPR)
        Z_CLRRADU(IPR) = Z_RADLU1(IPR) - Z_CLDRADU(IPR)
        Z_OLDCLD(IPR) = Z_CLDRADU(IPR)
        Z_OLDCLR(IPR) = Z_CLRRADU(IPR)
        Z_RAD(IPR) = 0.0_JPRB
      ENDIF
      Z_TTOT = 1.0_JPRB - Z_ATOT1(INDEX)
      Z_TRNS = 1.0_JPRB - P_ABSS1(INDEX)
      Z_CLDSRC = Z_BBUTOT1(INDEX) * Z_ATOT1(INDEX)

!- Separate RT equations for clear and cloudy streams      
      Z_CLDRADU(IPR) = Z_CLDRADU(IPR) * Z_TTOT + Z_CLDFRAC(I_LEV) * Z_CLDSRC
      Z_CLRRADU(IPR) = Z_CLRRADU(IPR) * Z_TRNS +(1.0_JPRB - Z_CLDFRAC(I_LEV)) * Z_GASSRC
!***

!- total sky upward flux
      Z_RADLU1(IPR) = Z_CLDRADU(IPR) + Z_CLRRADU(IPR)
      Z_URAD1 = Z_URAD1 + Z_RADLU1(IPR)
      
!- clear-sky upward flux
      Z_RADCLRU1(IPR) = Z_RADCLRU1(IPR) + (Z_BBU1(INDEX)-Z_RADCLRU1(IPR))&
       & *P_ABSS1(INDEX)  
      Z_URADCL1 = Z_URADCL1 + Z_RADCLRU1(IPR)

!* Code to account for maximum/random overlap:
!   Performs RT on the radiance most recently switched between clear and
!   cloudy streams
      Z_RADMOD = Z_RAD(IPR) * (Z_FACCLR1(I_LEV+1) * Z_TRNS +&
       & Z_FACCLD1(I_LEV+1) *  Z_TTOT) - &
       & Z_FACCMB1(I_LEV+1) * Z_GASSRC + &
       & Z_FACCMB2(I_LEV+1) * Z_CLDSRC  
       
!   Computes what the clear and cloudy streams would have been had no
!   radiance been switched       
      Z_OLDCLD(IPR) = Z_CLDRADU(IPR) - Z_RADMOD
      Z_OLDCLR(IPR) = Z_CLRRADU(IPR) + Z_RADMOD
      
!   Computes the radiance to be switched between clear and cloudy.      
      Z_RAD(IPR) = -Z_RADMOD + Z_FACCLR2(I_LEV+1)*Z_OLDCLR(IPR) -&
       & Z_FACCLD2(I_LEV+1)*Z_OLDCLD(IPR)  
      Z_CLDRADU(IPR) = Z_CLDRADU(IPR) + Z_RAD(IPR)
      Z_CLRRADU(IPR) = Z_CLRRADU(IPR) - Z_RAD(IPR)
!***
    ENDDO

  ELSE

!- *** Clear layer
    IENT = JPGPT * (I_LEV-1)
    DO IPR = 1, JPGPT
      INDEX = IENT + IPR
!- total-sky upward flux          
      Z_RADLU1(IPR) = Z_RADLU1(IPR)+(Z_BBU1(INDEX)-Z_RADLU1(IPR))*P_ABSS1(INDEX)
      Z_URAD1 = Z_URAD1 + Z_RADLU1(IPR)
!- clear-sky upward flux
!   Upward clear and total sky streams must be separate because surface
!   reflectance is different for each.
      Z_RADCLRU1(IPR) = Z_RADCLRU1(IPR)+(Z_BBU1(INDEX)-Z_RADCLRU1(IPR))*P_ABSS1(INDEX)
      Z_URADCL1 = Z_URADCL1 + Z_RADCLRU1(IPR)
    ENDDO

  ENDIF

  P_TOTUFLUC(I_LEV) = Z_URADCL1 * Z_WTNUM(1)
  P_TOTUFLUX(I_LEV) = Z_URAD1   * Z_WTNUM(1)

ENDDO

!* Convert radiances to fluxes and heating rates for total and clear sky.
! ** NB: moved to calling routine
!      TOTUFLUC(0) = TOTUFLUC(0) * FLUXFAC
!      TOTDFLUC(0) = TOTDFLUC(0) * FLUXFAC
!      TOTUFLUX(0) = TOTUFLUX(0) * FLUXFAC
!      TOTDFLUX(0) = TOTDFLUX(0) * FLUXFAC

!      CLFNET(0) = (P_TOTUFLUC(0) - P_TOTDFLUC(0))
!      FNET(0)   = (P_TOTUFLUX(0) - P_TOTDFLUX(0))
!      DO LEV = 1, KLEV
!        TOTUFLUC(LEV) = TOTUFLUC(LEV) * FLUXFAC
!        TOTDFLUC(LEV) = TOTDFLUC(LEV) * FLUXFAC
!         CLFNET(LEV) =(P_TOTUFLUC(LEV) - P_TOTDFLUC(LEV))

!        TOTUFLUX(LEV) = TOTUFLUX(LEV) * FLUXFAC
!        TOTDFLUX(LEV) = TOTDFLUX(LEV) * FLUXFAC
!        FNET(LEV) = (P_TOTUFLUX(LEV) - P_TOTDFLUX(LEV))
!        L = LEV - 1

!- Calculate Heating Rates.
!        CLHTR(L)=HEATFAC*(CLFNET(L)-CLFNET(LEV))/(PZ(L)-PZ(LEV)) 
!        HTR(L)  =HEATFAC*(FNET(L)  -FNET(LEV))  /(PZ(L)-PZ(LEV)) 
!      END DO
!      CLHTR(KLEV) = 0.0
!      HTR(KLEV)   = 0.0



IF (LHOOK) CALL DR_HOOK('RRTM_RTRN1A_140GP',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_RTRN1A_140GP
