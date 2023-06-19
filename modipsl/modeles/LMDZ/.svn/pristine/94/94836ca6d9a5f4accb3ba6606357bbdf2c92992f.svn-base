SUBROUTINE SWTT1 ( KIDIA,KFDIA,KLON,KNU,KABS,KIND, PU, PTR )

!**** *SWTT1* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
!     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
!     INTERVALS.

!**   INTERFACE.
!     ----------
!          *SWTT1* IS CALLED FROM *SW1S*.

!        EXPLICIT ARGUMENTS :
!        --------------------
! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
! KABS   :                     ; NUMBER OF ABSORBERS
! KIND   : (KABS)              ; INDICES OF THE ABSORBERS
! PU     : (KLON,KABS)         ; ABSORBER AMOUNT
!     ==== OUTPUTS ===
! PTR    : (KLON,KABS)         ; TRANSMISSION FUNCTION

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
!     AND HORNER'S ALGORITHM.

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
!        ORIGINAL : 95-01-20
!        03-10-10 Deborah Salmond and Marta Janiskova Optimisation
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
   
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOESW    , ONLY : APAD     ,BPAD     ,D

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KABS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIND(KABS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KABS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTR(KLON,KABS) 
!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-----------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZR1(KLON), ZR2(KLON), ZU(KLON)
REAL(KIND=JPRB) :: ZRR

INTEGER(KIND=JPIM) :: IA, JA, JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

!*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION

IF (LHOOK) CALL DR_HOOK('SWTT1',0,ZHOOK_HANDLE)
DO JA = 1,KABS
  IA=KIND(JA)
! print *,'SWTT1: KNU', KNU
  DO JL = KIDIA,KFDIA
    ZU(JL) = PU(JL,JA)
    ZR1(JL) = APAD(KNU,IA,1) + ZU(JL) * (APAD(KNU,IA,2) + ZU(JL)&
     & * ( APAD(KNU,IA,3) + ZU(JL) * (APAD(KNU,IA,4) + ZU(JL)&
     & * ( APAD(KNU,IA,5) + ZU(JL) * (APAD(KNU,IA,6) + ZU(JL)&
     & * ( APAD(KNU,IA,7) ))))))  
!    print *,'SWTT1 ZU APAD',IA,ZU(JL),APAD(KNU,IA,1),APAD(KNU,IA,2),&
!    &APAD(KNU,IA,3),APAD(KNU,IA,4),APAD(KNU,IA,5),APAD(KNU,IA,6),APAD(KNU,IA,7)

    ZR2(JL) = BPAD(KNU,IA,1) + ZU(JL) * (BPAD(KNU,IA,2) + ZU(JL)&
     & * ( BPAD(KNU,IA,3) + ZU(JL) * (BPAD(KNU,IA,4) + ZU(JL)&
     & * ( BPAD(KNU,IA,5) + ZU(JL) * (BPAD(KNU,IA,6) + ZU(JL)&
     & * ( BPAD(KNU,IA,7) ))))))  
    ZRR=1.0_JPRB/ZR2(JL)

!*         2.      ADD THE BACKGROUND TRANSMISSION

    PTR(JL,JA) = (ZR1(JL)*ZRR) * (1.0_JPRB-D(KNU,IA)) + D(KNU,IA)
  ENDDO
ENDDO
!WRITE(*,'("---> Dans SWTT1, PTR : "10E12.5)') (PTR(1,JA),JA=1,KABS)

IF (LHOOK) CALL DR_HOOK('SWTT1',1,ZHOOK_HANDLE)
END SUBROUTINE SWTT1
