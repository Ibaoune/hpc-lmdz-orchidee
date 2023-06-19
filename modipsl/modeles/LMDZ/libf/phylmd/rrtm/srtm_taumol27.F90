SUBROUTINE SRTM_TAUMOL27 &
 & ( KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1     , P_ONEMINUS,&
 & P_COLMOL  , P_COLO3,&
 & K_LAYTROP,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    &
 & )  

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 27:  29000-38000 cm-1 (low - O3; high - O3)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARSRTM  , ONLY : JPLAY, JPG, NG27
USE YOESRTA27, ONLY : ABSA, ABSB &
 & , SFLUXREFC, RAYLC &
 & , LAYREFFR, SCALEKUR  
USE YOESRTWN , ONLY : NSPA, NSPB

IMPLICIT NONE

!-- Output
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC00(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC01(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC10(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC11(JPLAY) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JP(JPLAY) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT(JPLAY) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT1(JPLAY) 
REAL(KIND=JPRB)                  :: P_ONEMINUS ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLO3(JPLAY) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SFLUXZEN(JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUG(JPLAY,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUR(JPLAY,JPG) 
!- from INTFAC      
!- from INTIND
!- from PRECISE             
!- from PROFDATA             
!- from SELF             
INTEGER(KIND=JPIM) :: IG, IND0, IND1, I_LAY, I_LAYSOLFR, I_NLAYERS

REAL(KIND=JPRB) ::  &
 & Z_TAURAY  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL27',0,ZHOOK_HANDLE)
I_NLAYERS = KLEV

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  
DO I_LAY = 1, K_LAYTROP
  IND0 = ((K_JP(I_LAY)-1)*5+(K_JT(I_LAY)-1))*NSPA(27) + 1
  IND1 = (K_JP(I_LAY)*5+(K_JT1(I_LAY)-1))*NSPA(27) + 1

!  DO IG = 1, NG(27)
  DO IG = 1 , NG27
    Z_TAURAY = P_COLMOL(I_LAY) * RAYLC(IG)
    P_TAUG(I_LAY,IG) = P_COLO3(I_LAY) * &
     & (P_FAC00(I_LAY) * ABSA(IND0,IG) + &
     & P_FAC10(I_LAY) * ABSA(IND0+1,IG) + &
     & P_FAC01(I_LAY) * ABSA(IND1,IG) + &
     & P_FAC11(I_LAY) * ABSA(IND1+1,IG))  
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    P_TAUR(I_LAY,IG) = Z_TAURAY
  ENDDO
ENDDO

I_LAYSOLFR = I_NLAYERS

DO I_LAY = K_LAYTROP+1, I_NLAYERS
  IF (K_JP(I_LAY-1) < LAYREFFR .AND. K_JP(I_LAY) >= LAYREFFR) &
   & I_LAYSOLFR = I_LAY  
  IND0 = ((K_JP(I_LAY)-13)*5+(K_JT(I_LAY)-1))*NSPB(27) + 1
  IND1 = ((K_JP(I_LAY)-12)*5+(K_JT1(I_LAY)-1))*NSPB(27) + 1

!  DO IG = 1, NG(27)
  DO IG = 1 , NG27
    Z_TAURAY = P_COLMOL(I_LAY) * RAYLC(IG)
    P_TAUG(I_LAY,IG) = P_COLO3(I_LAY) * &
     & (P_FAC00(I_LAY) * ABSB(IND0,IG) + &
     & P_FAC10(I_LAY) * ABSB(IND0+1,IG) + &
     & P_FAC01(I_LAY) * ABSB(IND1,IG) +  &
     & P_FAC11(I_LAY) * ABSB(IND1+1,IG))  
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (I_LAY == I_LAYSOLFR) P_SFLUXZEN(IG) = SCALEKUR * SFLUXREFC(IG) 
    P_TAUR(I_LAY,IG) = Z_TAURAY
  ENDDO
ENDDO

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL27',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_TAUMOL27

