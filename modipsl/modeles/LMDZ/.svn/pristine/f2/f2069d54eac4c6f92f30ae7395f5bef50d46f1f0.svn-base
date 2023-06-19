MODULE PHY_SV

IMPLICIT NONE

! Constants
                                                                              
      CHARACTER(LEN=1), PARAMETER :: labnum(0:9) =                            &
              (/'0','1','2','3','4','5','6','7','8','9'/)             
!C +...                labnum: Alphanumeric Character                          
                                                                               
      INTEGER, PARAMETER ::  izr=0,iun=1                                             
      REAL, PARAMETER    ::  zer0 = 0.0e+0, half = 0.5e+0, un_1 = 1.0e+0,     &
     &                       eps6 = 1.0e-6, epsn = 1.0e-9, R_1000=1.e3             
      REAL               ::  ea_MAX,ea_MIN 
!C +...                    zer0 = 0.0e+0                                       
!C +                       half = 0.5e+0                                 
!C +                       un_1 = 1.0e+0                             
!C +                       eps6 = 1.0e-6  ! Arbirary Small Value           !    1.e-6     
!C +                       epsn = 1.0e-9  ! Arbirary Small Value           !    1.e-6 
!C +                       R_1000=1.e3     
!C +                       ea_MIN   ! MAX allowed exponential Argum. ! computed by HOST     
!C +                       ea_MIN   ! MIN allowed exponential Argum. ! computed by HOST     
                                                                              
      REAL, PARAMETER    ::  piNmbr = 3.141592653589793238462643e0
      REAL, SAVE         ::  Dg2Rad = piNmbr / 180.
!C +...                    piNmbr   ! pi    = 3.141592653589793238462643e0             
!C +...                    Dg2rad   ! pi / 180                                       
      REAL, PARAMETER    ::  Grav_F = 9.81e0, Grav_I = 1./Grav_F                       
!C +                   Grav_F:  Gravitational  Force         =  9.81    m/s2           
!C +                   Grav_I: 1 / Grav_F                    = 1 /9.81    s2/m    
                                                                              
      REAL, PARAMETER    ::  A_MolV = 1.35e-5, vonKrm = 0.40e0, A_Stab = 5.8,  &
     &                       AhStab = 5.4,     AsStab = 4.0,    r_Stab = 3.0            
!C +...                A_MolV: Air Viscosity                 = 1.35d-5 m2/s   
!C +                   vonKrm: von Karman constant           = 0.4            
!C +                   A_Stab: Stability  Coefficient Moment = 5.8             
!C +                   AhStab: Stability  Coefficient Heat   = 5.4             
!C +                   AsStab: Stability  Coefficient Blown * = 4.0            
!C +                   r_Stab: Turbulent Diffusivities Ratio K*/Km             
                                                                              
      REAL, PARAMETER    ::  R_DAir=287.05967, CpdAir=1004.708845,             &
     &                       RCp=R_DAir/CpdAir,p0_kap=3.73010e0,               &
     &                       StefBo=5.67e-8         
      REAL, PARAMETER    ::  LhfH2O=3.34e+5, LhvH2O=2.5008e+6,&
     &                       LhsH2O=2.8345e+6      
!C +...                R_DAir: perfect gas law constant for dry air (287 J/kg/K) 
!C +                   CpdAir: dry air specific heat at constant p (1004 J/kg/K) 
!C +                   RCp   : R / Cp                      
!C +                   p0_kap: 100 kPa ** (R/Cp)                                    
!C +                   LhfH2O: Latent Heat of Fusion of Snow  =    3.34d+5 J/kg
!C +                   LhvH2O: Latent Heat Vaporisation Snow  = 2500.00d+3 J/kg  
!C +                   LhsH2O: Latent Heat Sublimation  Snow  = 2833.60d+3 J/kg  
!C +                   StefBo: Stefan-Boltzman Constant       =    5.67d-8 W/m2/K4
                                                                                
      REAL, PARAMETER    ::  rhoWat=1000.00e0, hC_Wat=4186.00e0        
!C +                   rhoWat: Water Specific Mass            = 1000.00d+0 kg/m3 
!C +                   hC_Wat: Water Heat Capacity            = 4186.00d+0 J/kg/K
                                                                                
      REAL, PARAMETER    ::  rhoIce=920.e0, BSnoRo=3.30e+2                             
!C +...                rhoIce:  Density      of Pure Ice      =  920.00d+0 kg/m3
!C +                   BSnoRo : Blowed Snow Density           =  255.00d+0 kg/m3 
                                                                                
      REAL, PARAMETER    ::  Tf_Sno=273.16e+0,Tf_Sea=271.2e+0             
!C +...                Tf_Sno: Snow  Melting  Point           =  273.16    K     
!C +...                Tf_Sea: Sea   Melting  Point           =  271.2     K

                      
                                                                                

END MODULE PHY_SV
