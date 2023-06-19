MODULE VARdCP

IMPLICIT NONE
! +                  
!C +--Col de Porte specific Constants                                            
!C +  ===============================                                            
                                                                                
      LOGICAL :: ColPrt                ! Col de Porte Switch               
                                                                                
                                                                                
!C +--Fractions of total solar irradiances in 3 spectral intervals               
!C +  ------------------------------------------------------------               
                                                                                
      REAL, PARAMETER :: Dr_1SN=0.59, Dr_2SN=0.31, Dr_3SN=0.10               
      REAL, PARAMETER :: Df_1SN=0.95, Df_2SN=0.05, Df_3SN=0.00              
      REAL, PARAMETER :: Dfc1SN=0.66, Dfc2SN=0.27, Dfc3SN=0.07               
      REAL :: DirSol,DifSol,TotSol,Clouds                               

!     Dr_1SN,Dr_2SN,Dr_3SN  ! Direct  Radiation                 
!     Df_1SN,Df_2SN,Df_3SN  ! Diffuse Radiation, Clear  Sky     
!     Dfc1SN,Dfc2SN,Dfc3SN  ! Diffuse Radiation, Cloudy Sky     
!C +--DATA - from common block COL DE PORTE                            
!C +  ====                                                                       
!c #CP data     ColPrt /.true./                                   
!                                                                                
!C +.           0.3--0.8micr.m   0.8--1.5micr.m   1.5--2.8micr.m                 
!C +            Fractions of total solar irradiance in 3 spectral intervals      
!C +***        (see Eric Martin Sept. 1996, CROCUS, Subroutine METEO)          
                                                                                
END MODULE VARdCP   