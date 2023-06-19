MODULE LMDZ_KL

USE VAR_SV, only : klonv

IMPLICIT NONE
! +
    REAL, DIMENSION(:),SAVE,ALLOCATABLE :: SOsoKL    ! Abs Solar Radiation 
!$OMP THREADPRIVATE(SOsoKL)
    REAL, DIMENSION(:),SAVE,ALLOCATABLE :: IRsoKL    ! Abs IR    Radiation 
!$OMP THREADPRIVATE(IRsoKL)
    REAL, DIMENSION(:),SAVE,ALLOCATABLE :: HSsoKL    ! Abs Sensible Ht Flux 
!$OMP THREADPRIVATE(HSsoKL)
    REAL, DIMENSION(:),SAVE,ALLOCATABLE :: HLsoKL    ! Abs Latent Heat Flux 
!$OMP THREADPRIVATE(HLsoKL)
    REAL, DIMENSION(:),SAVE,ALLOCATABLE :: HLs_KL    ! Evaporation          
!$OMP THREADPRIVATE(HLs_KL)
    REAL, DIMENSION(:),SAVE,ALLOCATABLE :: HLv_KL    ! Transpiration        
!$OMP THREADPRIVATE(HLv_KL)
    
!   REAL, DIMENSION(klonv)                 :: SOsoNC    ! Abs Solar Radiation  
!   REAL, DIMENSION(klonv)                 :: IRsoNC    ! Abs IR    Radiation  
!   REAL, DIMENSION(klonv)                 :: HSsoNC    ! Abs Sensible Ht Flux
!   REAL, DIMENSION(klonv)                 :: HLsoNC    ! Abs Latent Heat Flux
!   REAL, DIMENSION(klonv)                 :: HLs_NC    ! Evaporation      
!   REAL, DIMENSION(klonv)                 :: HLv_NC    ! Transpiration       

CONTAINS
  
  SUBROUTINE INIT_LMDZ_KL
  IMPLICIT NONE
  
    ALLOCATE(SOsoKL(klonv))    ! Abs Solar Radiation 
    ALLOCATE(IRsoKL(klonv))    ! Abs IR    Radiation 
    ALLOCATE(HSsoKL(klonv))    ! Abs Sensible Ht Flux 
    ALLOCATE(HLsoKL(klonv))    ! Abs Latent Heat Flux 
    ALLOCATE(HLs_KL(klonv))    ! Evaporation          
    ALLOCATE(HLv_KL(klonv))    ! Transpiration        
       
  END SUBROUTINE INIT_LMDZ_KL
  
  
END MODULE LMDZ_KL

