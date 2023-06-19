MODULE VARySV

USE VAR_SV, only : klonv, nsol, nsno


IMPLICIT NONE

! +
      INTEGER, DIMENSION(:),SAVE,ALLOCATABLE ::   NLaysv  ! New   Snow     Layer   Switch   
!$OMP THREADPRIVATE(NLaysv)
      INTEGER, DIMENSION(:),SAVE,ALLOCATABLE ::   i_thin  ! Index of the thinest Layer      
!$OMP THREADPRIVATE(i_thin)
      INTEGER, DIMENSION(:),SAVE,ALLOCATABLE ::   LIndsv  ! Contiguous Layer relative Index
!$OMP THREADPRIVATE(LIndsv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   albisv  ! Integrated Surface Albedo
!$OMP THREADPRIVATE(albisv)
      !hj prelim                                            
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   alb1sv  ! Surface Albedo VIS      
!$OMP THREADPRIVATE(alb1sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   alb2sv  ! Surface Albedo NIR      
!$OMP THREADPRIVATE(alb2sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   alb3sv  ! Surface Albedo FIR      
!$OMP THREADPRIVATE(alb3sv)

      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   albssv  ! Soil               Albedo [-]   
!$OMP THREADPRIVATE(albssv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   SoCasv  ! Canopy  Absorbed Solar Radiat.  
!$OMP THREADPRIVATE(SoCasv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   SoSosv  !? Surface Absorbed Solar Radiat.  
!$OMP THREADPRIVATE(SoSosv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   IRv_sv  ! Vegetation IR Flux  [W/m2]      
!$OMP THREADPRIVATE(IRv_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   Evg_sv  ! Emissivity of Vegetation+Snow   
!$OMP THREADPRIVATE(Evg_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   Eso_sv  !? Emissivity of       Soil+Snow   
!$OMP THREADPRIVATE(Eso_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   tau_sv  ! Transmited Radiation Fraction   
!$OMP THREADPRIVATE(tau_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   rrMxsv  ! Canopy Maximum Intercepted Rain 
!$OMP THREADPRIVATE(rrMxsv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   LAIesv  ! effective LAI for transpirati.  
!$OMP THREADPRIVATE(LAIesv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   LAI_sv  ! corrected LAI in case of snow   
!$OMP THREADPRIVATE(LAI_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   glf_sv  ! Green  Leaf Fraction            
!$OMP THREADPRIVATE(glf_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   Sigmsv  ! Canopy Ventilation  Factor      
!$OMP THREADPRIVATE(Sigmsv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   HSv_sv  ! Sensible Heat Flux  [W/m2]      
!$OMP THREADPRIVATE(HSv_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   HLv_sv  ! Latent   Heat Flux  [W/m2]      
!$OMP THREADPRIVATE(HLv_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   HSs_sv  !? Sensible Heat Flux (t)          
!$OMP THREADPRIVATE(HSs_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   HLs_sv  !? Latent   Heat Flux (t)          
!$OMP THREADPRIVATE(HLs_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   sqrCm0  ! in Neutral Drag Coef.Moment.    
!$OMP THREADPRIVATE(sqrCm0)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   sqrCh0  ! in Neutral Drag Coef.Heat       
!$OMP THREADPRIVATE(sqrCh0)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   Lx_H2O  ! Latent Heat of Vaporiz./Sublim. 
!$OMP THREADPRIVATE(Lx_H2O)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   ram_sv  ! Aerodyn.Resistance (Moment.)    
!$OMP THREADPRIVATE(ram_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   rah_sv  ! Aerodyn.Resistance (Heat)       
!$OMP THREADPRIVATE(rah_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   Fh__sv  ! Stability Function              
!$OMP THREADPRIVATE(Fh__sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   dFh_sv  ! Stability Function (Deriv.)     
!$OMP THREADPRIVATE(dFh_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   Evp_sv  !x Evaporation        [kg/m2]      
!$OMP THREADPRIVATE(Evp_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   EvT_sv  !x Evapotranspiration [kg/m2]      
!$OMP THREADPRIVATE(EvT_sv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   LSdzsv  ! Land/Sea Vert. Discretiz. Fact. 
!$OMP THREADPRIVATE(LSdzsv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   Tsrfsv  ! Surface    Temperature          
!$OMP THREADPRIVATE(Tsrfsv)
      REAL,SAVE,ALLOCATABLE  ::   sEX_sv(:,:)  ! Verticaly Integr.Extinct.Coef.  
!$OMP THREADPRIVATE(sEX_sv)
      REAL,SAVE,ALLOCATABLE  ::   zzsnsv(:,:)    ! Snow  Pack Thickness      [m]   
!$OMP THREADPRIVATE(zzsnsv)
      REAL,SAVE,ALLOCATABLE  ::   psi_sv(:,:)    ! Soil   Water        Potential   
!$OMP THREADPRIVATE(psi_sv)
      REAL,SAVE,ALLOCATABLE  ::   Khydsv(:,:)    ! Soil   Hydraulic    Conductiv.  
!$OMP THREADPRIVATE(Khydsv)
      REAL,SAVE,ALLOCATABLE  ::   Rootsv(:,:)       ! Root Water Pump      [kg/m2/s]  
!$OMP THREADPRIVATE(Rootsv)
      REAL, DIMENSION(:),SAVE,ALLOCATABLE    ::   EExcsv  ! Energy in Excess, current       
!$OMP THREADPRIVATE(EExcsv)

CONTAINS
  SUBROUTINE INIT_VARySV
  IMPLICIT NONE
      INTEGER                ::   ikl

      ALLOCATE(NLaysv(klonv))  ! New   Snow     Layer   Switch   
      ALLOCATE(i_thin(klonv))  ! Index of the thinest Layer      
      ALLOCATE(LIndsv(klonv))  ! Contiguous Layer relative Index
      ALLOCATE(albisv(klonv))  ! Integrated Surface Albedo
      !hj prelim                                            
      ALLOCATE(alb1sv(klonv))  ! Surface Albedo VIS      
      ALLOCATE(alb2sv(klonv))  ! Surface Albedo NIR      
      ALLOCATE(alb3sv(klonv))  ! Surface Albedo FIR       

      !
      ALLOCATE(albssv(klonv))  ! Soil               Albedo [-]   
      ALLOCATE(SoCasv(klonv))  ! Canopy  Absorbed Solar Radiat.  
      ALLOCATE(SoSosv(klonv))  !? Surface Absorbed Solar Radiat.  
      ALLOCATE(IRv_sv(klonv))  ! Vegetation IR Flux  [W/m2]      
      ALLOCATE(Evg_sv(klonv))  ! Emissivity of Vegetation+Snow   
      ALLOCATE(Eso_sv(klonv))  !? Emissivity of       Soil+Snow   
      ALLOCATE(tau_sv(klonv))  ! Transmited Radiation Fraction   
      ALLOCATE(rrMxsv(klonv))  ! Canopy Maximum Intercepted Rain 
      ALLOCATE(LAIesv(klonv))  ! effective LAI for transpirati.  
      ALLOCATE(LAI_sv(klonv))  ! corrected LAI in case of snow   
      ALLOCATE(glf_sv(klonv))  ! Green  Leaf Fraction            
      ALLOCATE(Sigmsv(klonv))  ! Canopy Ventilation  Factor      
      ALLOCATE(HSv_sv(klonv))  ! Sensible Heat Flux  [W/m2]      
      ALLOCATE(HLv_sv(klonv))  ! Latent   Heat Flux  [W/m2]      
      ALLOCATE(HSs_sv(klonv))  !? Sensible Heat Flux (t)          
      ALLOCATE(HLs_sv(klonv))  !? Latent   Heat Flux (t)          
      ALLOCATE(sqrCm0(klonv))  ! in Neutral Drag Coef.Moment.    
      ALLOCATE(sqrCh0(klonv))  ! in Neutral Drag Coef.Heat       
      ALLOCATE(Lx_H2O(klonv))  ! Latent Heat of Vaporiz./Sublim. 
      ALLOCATE(ram_sv(klonv))  ! Aerodyn.Resistance (Moment.)    
      ALLOCATE(rah_sv(klonv))  ! Aerodyn.Resistance (Heat)       
      ALLOCATE(Fh__sv(klonv))  ! Stability Function              
      ALLOCATE(dFh_sv(klonv))  ! Stability Function (Deriv.)     
      ALLOCATE(Evp_sv(klonv))  !x Evaporation        [kg/m2]      
      ALLOCATE(EvT_sv(klonv))  !x Evapotranspiration [kg/m2]      
      ALLOCATE(LSdzsv(klonv))  ! Land/Sea Vert. Discretiz. Fact. 
      ALLOCATE(Tsrfsv(klonv))  ! Surface    Temperature          
      ALLOCATE(sEX_sv(klonv,-nsol:nsno+1))  ! Verticaly Integr.Extinct.Coef.  
      ALLOCATE(zzsnsv(klonv,    0:nsno))    ! Snow  Pack Thickness      [m]   
      ALLOCATE(psi_sv(klonv,-nsol:0   ))    ! Soil   Water        Potential   
      ALLOCATE(Khydsv(klonv,-nsol:0   ))    ! Soil   Hydraulic    Conductiv.  
      ALLOCATE(Rootsv(klonv,-nsol:0))       ! Root Water Pump      [kg/m2/s]  
      ALLOCATE(EExcsv(klonv))  ! Energy in Excess, current         

      DO ikl=1,klonv
        NLaysv(ikl)   = 0
        i_thin(ikl)   = 0
        LIndsv(ikl)   = 0
      END DO

  END SUBROUTINE INIT_VARySV

END MODULE VARySV

