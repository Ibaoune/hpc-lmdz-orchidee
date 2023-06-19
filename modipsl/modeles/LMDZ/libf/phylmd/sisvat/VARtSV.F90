MODULE VARtSV

USE VAR_SV, only : nsol, nsno, klonv

IMPLICIT NONE
! +
      INTEGER, PARAMETER      ::   ntaver = 4
! +
      REAL,ALLOCATABLE,SAVE                ::   toicSV(:)    ! Snow to ice balance    
!$OMP THREADPRIVATE(toicSV)

      REAL,ALLOCATABLE,SAVE                ::   dz1_SV(:,:) ! "inverse" layer thicknes
!$OMP THREADPRIVATE(dz1_SV)
      REAL,ALLOCATABLE,SAVE                ::   dz2_SV(:,:) ! layer thickness
!$OMP THREADPRIVATE(dz2_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE  ::   Tsf_SV  ! surface temperature !hj12032010
!$OMP THREADPRIVATE(Tsf_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE  ::   TsfnSV  ! new surface temperature 
!$OMP THREADPRIVATE(TsfnSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE  ::   AcoHSV, BcoHSV ! coefficients for Enthalpy
!$OMP THREADPRIVATE(AcoHSV, BcoHSV)                             ! evolution,from atmosphere
      REAL, DIMENSION(:),ALLOCATABLE,SAVE  ::   AcoQSV, BcoQSV ! coefficients for Humidity
!$OMP THREADPRIVATE( AcoQSV, BcoQSV)                            ! evolution,from atmosphere
      REAL, DIMENSION(:),ALLOCATABLE,SAVE  ::   ps__SV  ! surface pressure 
!$OMP THREADPRIVATE(ps__SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE  ::   p1l_SV  ! 1st layer pressure 
!$OMP THREADPRIVATE(p1l_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE  ::   cdH_SV  ! drag coeff Energy (?)
!$OMP THREADPRIVATE(cdH_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE  ::   rsolSV  ! Radiation balance surface
!$OMP THREADPRIVATE(rsolSV)
      REAL,SAVE                            ::   lambSV  ! coefficient soil layers
!$OMP THREADPRIVATE(lambSV)

                                                                   ! #AW 
      REAL, ALLOCATABLE,SAVE               ::   V__mem(:,:)                         
!$OMP THREADPRIVATE(V__mem)                                        ! #AW
      REAL, ALLOCATABLE,SAVE               ::   VVmmem(:)                              
!$OMP THREADPRIVATE(VVmmem)                                        ! #AH 
      REAL, ALLOCATABLE,SAVE               ::   T__mem(:,:)                           
!$OMP THREADPRIVATE(T__mem)                                        ! #AH 
      REAL, ALLOCATABLE,SAVE               ::   dTmmem(:)                             
!$OMP THREADPRIVATE(dTmmem)
        
CONTAINS

  SUBROUTINE INIT_VARtSV
  IMPLICIT NONE
  
      ALLOCATE(toicSV(klonv))

      ALLOCATE(dz1_SV(klonv,-nsol:nsno))     ! "inverse" layer thicknes
      ALLOCATE(dz2_SV(klonv,-nsol:nsno))     ! layer thickness
      ALLOCATE(Tsf_SV(klonv))                ! surface temperature !hj12032010
      ALLOCATE(TsfnSV(klonv))                ! new surface temperature 
      ALLOCATE(AcoHSV(klonv), BcoHSV(klonv)) ! coefficients for Enthalpy
                                             ! evolution,from atmosphere
      ALLOCATE(AcoQSV(klonv), BcoQSV(klonv)) ! coefficients for Humidity
                                             ! evolution,from atmosphere
      ALLOCATE(ps__SV(klonv))                ! surface pressure 
      ALLOCATE(p1l_SV(klonv))                ! 1st layer pressure 
      ALLOCATE(cdH_SV(klonv))                ! drag coeff Energy (?)
      ALLOCATE(rsolSV(klonv))                ! Radiation balance surface

                                                                   ! #AW
      ALLOCATE(V__mem(klonv,ntaver))                               ! #AW
      ALLOCATE(VVmmem(klonv))                                      ! #AH
      ALLOCATE(T__mem(klonv,ntaver))                               ! #AH
      ALLOCATE(dTmmem(klonv))
  END SUBROUTINE INIT_VARtSV

END MODULE VARtSV

