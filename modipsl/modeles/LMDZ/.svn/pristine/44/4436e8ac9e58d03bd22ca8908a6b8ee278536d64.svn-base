MODULE VARphy

IMPLICIT NONE

! MAR Constants
                                                                              
      CHARACTER(LEN=1), PARAMETER :: labnum(0:9) =                            &
              (/'0','1','2','3','4','5','6','7','8','9'/)             
!C +...                labnum: Alphanumeric Character                          
                                                                               
      INTEGER, PARAMETER ::  izr=0,iun=1                                             
      REAL, PARAMETER    ::  zero = 0.0e+0, demi = 0.5e+0, unun = 1.0e+0,     &
     &                       epsi = 1.0e-6, eps9 = 1.0e-9, eps12= 1.0e-12,    &
     &                       third=1/3,     thous=1.e3             
      REAL               ::  argmin,argmax      !hj valeurs??
!C +...                    zero = 0.0e+0                                       
!C +                            demi = 0.5e+0                                 
!C +                                 unun = 1.0e+0                             
!C +                                      epsi = 1.0e-6                        
!C +                                           eps9 = 1.0e-9                   
!C +                                                eps12= 1.0e-12             
!C +                                                      third=1/3            
!C +                                                            thous=1.e3     
!C +                       argmin:        Function exp(x) Minimum Argument     
!C +                              argmax: Function exp(x) Maximum Argument     
!C +                                      (Min/Max are Machine Dependant)      
                                                                              
      REAL, PARAMETER    ::  pi = 3.141592653589793238462643e0,               &
     &                       degrad = pi / 180.d0, hourad = pi /  12.d0 
!C +...                pi    = 3.141592653589793238462643e0                   
!C +                   degrad= pi / 180.d0                                     
!C +                   hourad= pi /  12.d0                                    
                                                                               
      REAL, PARAMETER    ::  earthr = 6371.229e+3, earthv = 729.217e-7,       &
     &                       gravit = 9.81e0,      gravi2 = gravit**2 ,       &
     &                       grvinv = 1./gravit                       
!C +...                earthr: Earth Radius                  =6371.d+3    m    
!C +...                earthv: Earth Angular Velocity        = 7.29d-5  s-1   
!C +                   gravit: Earth Gravity Acceleration    = 9.81    m/s2    
!C +                   gravi2: idem (squared)                                  
!C +                   grvinv: idem (inverse)                                  
                                                                              
      REAL, PARAMETER    ::  akmol = 1.35e-5, vonkar = 0.40e0, A_Turb = 5.8,  &
     &                       AhTurb = 5.4,    AsTurb = 4.0,    r_Turb = 3.0            
!C +...                akmol : Air Viscosity                 = 1.35d-5 m2/s   
!C +                   vonkar: von Karman constant           = 0.4            
!C +                   A_Turb: Stability  Coefficient Moment = 5.8             
!C +                   AhTurb: Stability  Coefficient Heat   = 5.4             
!C +                   AsTurb: Stability  Coefficient BLOW * = 4.0            
!C +                   r_Turb: Turbulent Diffusivities Ratio K*/Km             
                                                                              
!      REAL    ::  RDryAi,cp,CvDrya,racv,cap,pcap,stefan,qv_MIN           
      REAL, PARAMETER    ::  RDryAi=287.05967, ra=287.05967, cp=1004.708845,  &
     &                       CvDrya=Cp-RDryAi, racv=RDryAi/CvDryA,            &
     &                       cap=0.28586e0,    pcap=3.73010e0,                &
     &                       stefan=5.67e-8,   qv_MIN=3.00e-6          
      REAL, PARAMETER    ::  RVapor=461.e0,    epsq=3.00e-6, Lv_H2O=2.5008e+6,&
     &                       Ls_H2O=2.8345e+6, r_LvCp=2490.04,r_LcCp=332.27,  &
     &                       r_LsCp=2822.31            
!C +...                RDryAi: perfect gas law constant for dry air (287 J/kg/K) 
!C +                   ra    : perfect gas law constant for dry air (287 J/kg/K) 
!C +                   cp    : dry air specific heat at constant p (1004 J/kg/K) 
!C +                   CvDryA: dry air specific heat at constant V ( 717 J/kg/K) 
!C +                   racv  = ra / cv      = ra /(cp-ra)                        
!C +                   cap   = ra / cp                                           
!C +                   pcap  = 100.[kPa]**cap                                    
!C +...                RVapor: perfect gas law constant Water Vapor (461 J/kg/K) 
!C +                   epsq  : Minimum Water Vapor Content    =    3.00d-6 kg/kg 
!C +                   Lv_H2O: Latent Heat Vaporisation/Water = 2500.00d+3 J/kg  
!C +                   Ls_H2O: Latent Heat Sublimation /Ice   = 2833.60d+3 J/kg  
!C +                   r_LvCp [Lv=2500000J/kg]/[Cp=1004J/kg/K]= 2490.04    K/[kg/
!C +                   r_LvCp [Lc= 333600J/kg]/[Cp=1004J/kg/K]=  332.27    K/[kg/
!C +                   r_LvCp [Ls=2833600J/kg]/[Cp=1004J/kg/K]= 2822.31    K/[kg/
!C +                   stefan: Stefan Constant                =    5.67d-8 W/m2/K
!C +                   qv_MIN: Minimum Specific Humidity      =    3.00d-6 kg/kg 
                                                                                
      REAL, PARAMETER    ::  tfrwat=271.20e0,  siclf=302.00e6,cdsice=2.04e0,  &
     &                       hic0=1.00e0, ro_Wat=1000.00e0, C__Wat=4186.00e0, &
     &                       fracoh=0.75e0             
!C +...                tfrwat: Sea-Water freezing Temperature =  271.20d+0 K     
!C +                   siclf : Sea-Ice Heat of Fusion         =  302.00d+6 J/m3  
!C +                   cdsice: Sea-Ice Thermal Conductivity   =    2.04d+0 W/mK  
!C +                   hic0  : Sea-Ice Initial Thickness      =    2.00d+0 m     
!C +                   ro_Wat: Density       of Water         = 1000.00d+0 kg/m3 
!C +                   C__Wat: Heat Capacity of Water         = 4186.00d+0 J/kg/K
!C +                   fracoh: 25% of cooling = Oceanic Heat Flux                
!C +***                        (Hibler 1984)    (ANTARCTIC Ocean)                
                                                                                
      REAL, PARAMETER    ::  ro_Ice=920.e0, cdice=2.10e0                                 
!C +...                ro_Ice:  Density      of Pure Ice       =  920.00d+0 kg/m3
!C +                   cdice :  Conductivity of Pure Ice       =    2.51d+0 W/m/K
                                                                                
      REAL, PARAMETER    ::  TfSnow=273.15e+0,csnow=2105.00e+0,r0sno=3.00e+1, &
     &                       blsno=3.30e+2,   Lf_H2O=3.337e+5                          
!C +...                TfSnow:        Snow melting Temperature=  273.15d+0 K     
!C +                    csnow:Heat Capacity of Snow             2105      J/kg/K 
!C +                         (Loth et al. 1993, JGR 98 D6, 2.2.2 2e para p.10453)
!C +                   r0sno : Fresh  Snow Density                50.00d+0 kg/m3 
!C +                   blsno : Blowed Snow Density               330.00d+0 kg/m3 
!C +                   Lf_H2O: Latent Heat of Fusion of Snow  =    3.34d+5 J/kg  

! MAR Constants
                      
                                                                                

END MODULE VARphy
