SUBROUTINE checkmass(zq,RNAVOG,mass,zdz,pplay,t_seri,iscm3,comment)
  USE dimphy
USE geometry_mod , ONLY:cell_area
  IMPLICIT NONE

#include "YOMCST.h"

! Entrees
  REAL,DIMENSION(klon,klev), INTENT(IN)   :: zq
  REAL,DIMENSION(klon,klev), INTENT(IN)   :: zdz
  REAL,INTENT(IN)                         :: mass
  REAL,INTENT(IN)                         :: RNAVOG
  CHARACTER(LEN=*),INTENT(IN)             :: comment
REAL,DIMENSION(klon,klev) ,INTENT(IN)  :: pplay
REAL,DIMENSION(klon,klev) ,INTENT(IN) :: t_seri
LOGICAL :: iscm3


! Local  
!  INTEGER,DIMENSION(klon)     :: jadrs
  INTEGER                     :: i, k
  REAL                              :: totmass
  REAL,DIMENSION(klon)              ::  burden
  REAL,DIMENSION(klon,klev)   :: zqv


DO i=1, klon
burden(i)=0.0
  DO k=1, klev
zqv(i,k)=zq(i,k)

ENDDO

ENDDO

totmass=0.0


IF (iscm3) THEN

ELSE

CALL kg_to_cm3(pplay,t_seri,zqv)

ENDIF

  DO k=1, klev
      DO i=1, klon
        burden(i)=burden(i)+(zqv(i,k)*1.e6*zdz(i,k)* &
                 mass*1.e3/RNAVOG)     !--mg S/m2
    ENDDO
  ENDDO
    DO i=1, klon
        totmass=totmass+burden(i)*cell_area(i)
ENDDO
WRITE(*,9999) comment," totmass= ",totmass
9999  format(a20,a10,e15.8)

END SUBROUTINE checkmass
