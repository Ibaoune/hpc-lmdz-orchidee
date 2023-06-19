MODULE YOERRTRWT

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

REAL(KIND=JPRB) :: FREFA  (JPGPT,13)
REAL(KIND=JPRB) :: FREFB  (JPGPT,6)
REAL(KIND=JPRB) :: FREFADF(JPGPT,13)
REAL(KIND=JPRB) :: FREFBDF(JPGPT,6)
REAL(KIND=JPRB) :: RWGT   (JPG*JPBAND)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! FREFA  :  REAL   :
! FREFB  :  REAL   :
! FREFADF:  REAL   :
! FREFBDF:  REAL   :
! RWT    :  REAL   : 
!    -------------------------------------------------------------------
!$OMP THREADPRIVATE(frefa,frefadf,frefb,frefbdf,rwgt)
END MODULE YOERRTRWT
