MODULE VAR_SV

IMPLICIT NONE
! +
! +   LMDZ_SV dims
! +   ===================
! +
!      INTEGER,PARAMETER :: nsol=nsoilmx-1
!      INTEGER,PARAMETER :: nsno=nsnowmx
      INTEGER,PARAMETER    :: nsol=10      !8
      INTEGER,PARAMETER    :: nsno=35
      INTEGER,PARAMETER :: nb_wri=1
      INTEGER,SAVE      :: klonv
!$OMP THREADPRIVATE(klonv)
      INTEGER,SAVE      :: knonv
!$OMP THREADPRIVATE(knonv)

END MODULE VAR_SV
