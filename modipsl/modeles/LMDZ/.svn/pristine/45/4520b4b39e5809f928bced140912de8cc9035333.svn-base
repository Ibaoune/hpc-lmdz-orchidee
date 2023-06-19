SUBROUTINE SWUVO3 &
 & ( KIDIA,KFDIA,KLON,KNU,KABS,&
 & PU, PTR &
 & )  
  
!**** *SWUVO3* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR OZONE
!     IN THE UV and VISIBLE SPECTRAL INTERVALS.

!**   INTERFACE.
!     ----------
!          *SWUVO3* IS CALLED FROM *SW1S*.

!        EXPLICIT ARGUMENTS :
!        --------------------
! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
! KABS   :                     ; NUMBER OF ABSORBERS
! PU     : (KLON,KABS)         ; ABSORBER AMOUNT
!     ==== OUTPUTS ===
! PTR    : (KLON,KABS)         ; TRANSMISSION FUNCTION

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          TRANSMISSION FUNCTION ARE COMPUTED USING SUMS OF EXPONENTIALS

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------
!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 00-12-18
!        Modified J. HAGUE          03-01-03 MASS Vector Functions       
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
   
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOESW    , ONLY : NEXPO3, REXPO3
USE YOMJFH   , ONLY : N_VMASS
USE write_field_phy

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KABS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KABS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTR(KLON,KABS) 
!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-----------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZU(KLON)
REAL(KIND=JPRB) :: ZTMP1(KFDIA-KIDIA+1+N_VMASS)
REAL(KIND=JPRB) :: ZTMP2(KFDIA-KIDIA+1+N_VMASS)

INTEGER(KIND=JPIM) ::  JA, JL, IEXP, JX, JLEN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
LOGICAL LLDEBUG

IF (LHOOK) CALL DR_HOOK('SWUVO3',0,ZHOOK_HANDLE)
IEXP=NEXPO3(KNU)
LLDEBUG=.FALSE.

!print *,'Dans SWUVO3, N_VMASS= ',N_VMASS
IF(N_VMASS > 0) THEN
  JLEN=KFDIA-KIDIA+N_VMASS-MOD(KFDIA-KIDIA,N_VMASS)
  IF(KFDIA-KIDIA+1 /= JLEN) THEN
    ZTMP1(KFDIA-KIDIA+2:JLEN) = 0.0_JPRB
  ENDIF
ENDIF

DO JA = 1,KABS
  DO JL=KIDIA,KFDIA
    PTR(JL,JA)=0.0_JPRB
  ENDDO
  
! Ce qui concerne N_VMASS commente par MPL 20.11.08
! IF(N_VMASS <= 0) THEN ! Do not use Vector Mass

!       WRITE(*,'("---> Dans SWUVO3 ")')
    DO JX=1,IEXP
      DO JL = KIDIA,KFDIA
        ZU(JL) = PU(JL,JA)
        PTR(JL,JA) = PTR(JL,JA)+REXPO3(KNU,1,JX)*EXP(-REXPO3(KNU,2,JX)*ZU(JL))
!       WRITE(*,'("                 PTR ",E12.5)') (PTR(JL,JA))
!       WRITE(*,'("REXPO3-1 ",E12.5)') (REXPO3(KNU,1,JX))
!       WRITE(*,'("REXPO3-2 ",E12.5)') (REXPO3(KNU,2,JX))
!       WRITE(*,'("ZU ",E12.5)') (ZU(JL))
!       WRITE(*,'("KNU KABS IEXP ",3I6)') KNU,KABS,IEXP
      ENDDO
    ENDDO


! ELSE  ! Use Vector MASS

!   DO JX=1,IEXP
!     DO JL = KIDIA,KFDIA
!       ZTMP1(JL-KIDIA+1)=-REXPO3(KNU,2,JX)*PU(JL,JA)
!     ENDDO
  
!     CALL VEXP(ZTMP2,ZTMP1,JLEN)
  
!     DO JL = KIDIA,KFDIA
!       PTR(JL,JA) = PTR(JL,JA)+REXPO3(KNU,1,JX)*ZTMP2(JL-KIDIA+1)
!     ENDDO
!   ENDDO    

! ENDIF

ENDDO

IF(LLDEBUG) THEN
    call writefield_phy("swuvo3_pu",pu,kabs)
    call writefield_phy("swuvo3_ptr",ptr,kabs)
ENDIF

IF (LHOOK) CALL DR_HOOK('SWUVO3',1,ZHOOK_HANDLE)
END SUBROUTINE SWUVO3
