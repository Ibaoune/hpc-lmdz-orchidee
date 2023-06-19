MODULE VARxSV                      

USE VAR_SV, only : klonv, nsol, nsno, nb_wri
                                                          
IMPLICIT NONE                                                                                
! +--SISVAT INPUT        Variables                                              
! +  -----------------------------                                              
  
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE ::   LSmask  ! Land-Sea   Mask                 
!$OMP THREADPRIVATE(LSmask)
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE ::   isotSV  ! Soil       Type                 
!$OMP THREADPRIVATE(isotSV)
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE ::   iWaFSV  ! Soil       Drainage:(1,0)=(y,n) 
!$OMP THREADPRIVATE(iWaFSV)
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE ::   ivgtSV  ! Vegetation Type                 
!$OMP THREADPRIVATE(ivgtSV)
                                                                              
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   coszSV  ! Cosine of Sun zenithal Angle    
!$OMP THREADPRIVATE(coszSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   sol_SV  ! Downward  Solar    Radiation    
!$OMP THREADPRIVATE(sol_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   IRd_SV  ! Downward  Longwave Radiation    
!$OMP THREADPRIVATE(IRd_SV)
                                                                                
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   drr_SV  ! Rain  Intensity       [kg/m2/s] 
!$OMP THREADPRIVATE(drr_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   dsn_SV  ! Snow  Intensity       [kg/m2/s] 
!$OMP THREADPRIVATE(dsn_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   dsnbSV  ! Idem, fraction, from Drift  [-] 
!$OMP THREADPRIVATE(dsnbSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   esnbSV  ! Idem, fraction, from Drift  [-]
!$OMP THREADPRIVATE(esnbSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   dbs_SV  ! Drift Amount            [kg/m2] 
!$OMP THREADPRIVATE(dbs_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   BrosSV  ! Buffer Snow Layer Density
!$OMP THREADPRIVATE(BrosSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   BG1sSV  ! Buffer Snow Layer Dendr/Sphe[-]
!$OMP THREADPRIVATE(BG1sSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   BG2sSV  ! Buffer Snow Layer Spher/Size[-][0.0001m]
!$OMP THREADPRIVATE(BG2sSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   dz0_SV  ! dz0(Sastrugi dh)            [m] 
!$OMP THREADPRIVATE(dz0_SV)
                                                                               
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   cld_SV  ! Cloudiness (seen from SBL)      
!$OMP THREADPRIVATE(cld_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   za__SV  ! SBL Height                      
!$OMP THREADPRIVATE(za__SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   VV__SV  !(SBL Top)  Wind Velocity         
!$OMP THREADPRIVATE(VV__SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   VVs_SV  !(Sastr,V)  Relevance
!$OMP THREADPRIVATE(VVs_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   RRs_SV  !(Sastr,V)  Counter
!$OMP THREADPRIVATE(RRs_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   DDs_SV  !(Sastr,V)  Angle
!$OMP THREADPRIVATE(DDs_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   TaT_SV  ! SBL Top   Temperature           
!$OMP THREADPRIVATE(TaT_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   ExnrSV  ! Exner     Potential             
!$OMP THREADPRIVATE(ExnrSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   dSdTSV  ! Sensible Heat Flux T Derivat.   
!$OMP THREADPRIVATE(dSdTSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   dLdTSV  ! Latent   Heat Flux T Derivat.   
!$OMP THREADPRIVATE(dLdTSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   rhT_SV  ! SBL Top   Air  Density          
!$OMP THREADPRIVATE(rhT_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   QaT_SV  ! SBL Top   Specific Humidity     
!$OMP THREADPRIVATE(QaT_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   dQa_SV  ! SBL Flux  Limitation of Qa      
!$OMP THREADPRIVATE(dQa_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   qsnoSV  ! SBL Mean  Snow       Content    
!$OMP THREADPRIVATE(qsnoSV)
                                                                            
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   LAI0SV  ! Nominal Leaf Area Index         
!$OMP THREADPRIVATE(LAI0SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   glf0SV  ! Green   Leaf Fraction           
!$OMP THREADPRIVATE(glf0SV)
                                                                               
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   alb0SV  ! Soil    Albedo                  
!$OMP THREADPRIVATE(alb0SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   slopSV  ! Snow/Ice/Soil-Water Surf. Slope 
!$OMP THREADPRIVATE(slopSV)

                                                              
      REAL,SAVE                      ::   zSBLSV  ! SBL Height (Initial Value)      
!$OMP THREADPRIVATE(zSBLSV)
      REAL,SAVE                      ::   dt__SV  ! Time Step                       
!$OMP THREADPRIVATE(dt__SV)
      CHARACTER (len=18),SAVE        ::   daHost  ! Date Host Model                 
!$OMP THREADPRIVATE(daHost)
                                                                               
                                                                               
! +--SISVAT INPUT/OUTPUT Variables                                             
! +  -----------------------------                                             
                                                                               
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE ::   isnoSV  ! Nb of Ice/Snow Layers           
!$OMP THREADPRIVATE(isnoSV)
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE ::   ispiSV  ! Uppermost superimposed ice      
!$OMP THREADPRIVATE(ispiSV)
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE ::   iiceSV  ! Nb of Ice      Layers           
!$OMP THREADPRIVATE(iiceSV)
      INTEGER ,ALLOCATABLE,SAVE    ::   istoSV(:,:)  ! Snow Layer     History          
!$OMP THREADPRIVATE(istoSV)
                                                                               
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   alb_SV  ! Surface-Canopy Albedo           
!$OMP THREADPRIVATE(alb_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   emi_SV  ! Surface-Canopy Emissivity       
!$OMP THREADPRIVATE(emi_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   IRs_SV  ! Soil           IR Flux          
!$OMP THREADPRIVATE(IRs_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   LMO_SV  ! Monin-Obukhov  Scale            
!$OMP THREADPRIVATE(LMO_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   us__SV  ! Friction       Velocity         
!$OMP THREADPRIVATE(us__SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   uts_SV  ! Temperature Turbulent Scale     
!$OMP THREADPRIVATE(uts_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   cutsSV  ! Temperature Turbulent Scale C.  
!$OMP THREADPRIVATE(cutsSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   uqs_SV  ! Spec.Humid. Turbulent Scale     
!$OMP THREADPRIVATE(uqs_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   uss_SV  ! Blow.Snow   Turbulent Scale     
!$OMP THREADPRIVATE(uss_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   usthSV  ! Blowing Snow Erosion Thresh.    
!$OMP THREADPRIVATE(usthSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   rCDmSV  ! Square  Root Contribut. Drag_m  
!$OMP THREADPRIVATE(rCDmSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   rCDhSV  ! Square  Root Contribut. Drag_h  
!$OMP THREADPRIVATE(rCDhSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0m_SV  ! Momentum     Roughness Length   
!$OMP THREADPRIVATE(Z0m_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0mmSV  !  z0(Momentum,    Time Mean) [m] 
!$OMP THREADPRIVATE(Z0mmSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0mnSV  !  z0(Momentum,    instanta.) [m] 
!$OMP THREADPRIVATE(Z0mnSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0roSV  ! Subgrid Topo Roughness Length   
!$OMP THREADPRIVATE(Z0roSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0SaSV  !  z0(Sastrugi  h)            [m] 
!$OMP THREADPRIVATE(Z0SaSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0e_SV  !  z0(Snow eroded)            [m] 
!$OMP THREADPRIVATE(Z0e_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0emSV  !  z0(Snow eroded, Time Mean) [m] 
!$OMP THREADPRIVATE(Z0emSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0enSV  !  z0(Snow eroded, instanta.) [m] 
!$OMP THREADPRIVATE(Z0enSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0h_SV  ! Heat         Roughness Length   
!$OMP THREADPRIVATE(Z0h_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0hmSV  !  z0(Heat,        Time Mean) [m] 
!$OMP THREADPRIVATE(Z0hmSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   Z0hnSV  !  z0(Heat,        instanta.) [m] 
!$OMP THREADPRIVATE(Z0hnSV)
                                                                          
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   snCaSV  ! Canopy  Snow   Thickness        
!$OMP THREADPRIVATE(snCaSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   rrCaSV  ! Canopy  Water  Content          
!$OMP THREADPRIVATE(rrCaSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   psivSV  ! Leaf    Water  Potential        
!$OMP THREADPRIVATE(psivSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   TvegSV  ! Vegetation     Temperature      
!$OMP THREADPRIVATE(TvegSV)
                                                                        
      REAL ,ALLOCATABLE,SAVE   ::   TsisSV(:,:)  ! Snow/Ice/Soil-Water Temperature 
!$OMP THREADPRIVATE(TsisSV)
      REAL ,ALLOCATABLE,SAVE   ::   ro__SV(:,:)  ! Snow/Ice/Soil-Water VolumicMass 
!$OMP THREADPRIVATE(ro__SV)
      REAL,ALLOCATABLE,SAVE    ::   eta_SV(:,:)  ! Snow/Ice/Soil     Water Content 
!$OMP THREADPRIVATE(eta_SV)
      REAL,ALLOCATABLE,SAVE    ::   G1snSV(:,:)  ! Snow Dendricity/Sphericity      
!$OMP THREADPRIVATE(G1snSV)
      REAL,ALLOCATABLE,SAVE    ::   G2snSV(:,:)  ! Snow Sphericity/Size            
!$OMP THREADPRIVATE(G2snSV)
      REAL,ALLOCATABLE,SAVE    ::   dzsnSV(:,:)  ! Snow Layer  Thickness           
!$OMP THREADPRIVATE(dzsnSV)
      REAL,ALLOCATABLE,SAVE    ::   agsnSV(:,:)  ! Snow Age                        
!$OMP THREADPRIVATE(agsnSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   BufsSV  ! Snow Buffer Layer               
!$OMP THREADPRIVATE(BufsSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   rusnSV  ! Surficial   Water               
!$OMP THREADPRIVATE(rusnSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   SWf_SV  ! Normalized  Decay               
!$OMP THREADPRIVATE(SWf_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   SWS_SV  ! Surficial Water Status          
!$OMP THREADPRIVATE(SWS_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   HFraSV  ! Frazil      Thickness           
!$OMP THREADPRIVATE(HFraSV)
                                                                            
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   zWE_SV  ! Current   Snow Thickness [mmWE] 
!$OMP THREADPRIVATE(zWE_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   zWEcSV  ! Compacted Snow Thickness [mmWE] 
!$OMP THREADPRIVATE(zWEcSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   wem_SV  ! Only Melting             [mmWE] 
!$OMP THREADPRIVATE(wem_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   wer_SV  ! Refreezing               [mmWE] 
!$OMP THREADPRIVATE(wer_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   wes_SV  ! Sublimation              [mmWE] 
!$OMP THREADPRIVATE(wes_SV)
                                                                             
                                                                              
! +--SISVAT OUTPUT       Variables                                            
! +  -----------------------------                                            
                                                                               
      INTEGER,DIMENSION(nb_wri),SAVE ::   no__SV  ! OUTPUT file Unit Number         
!$OMP THREADPRIVATE(no__SV)
      INTEGER,DIMENSION(nb_wri),SAVE ::   i___SV  ! OUTPUT point   i Coordinate    
!$OMP THREADPRIVATE(i___SV)
      INTEGER,DIMENSION(nb_wri),SAVE ::   j___SV  ! OUTPUT point   j Coordinate     
!$OMP THREADPRIVATE(j___SV)
      INTEGER,DIMENSION(nb_wri),SAVE ::   n___SV  ! OUTPUT point   n Coordinate     
!$OMP THREADPRIVATE(n___SV)
      INTEGER,DIMENSION(nb_wri),SAVE ::   lwriSV  ! OUTPUT point vec Index          
!$OMP THREADPRIVATE(lwriSV)
!                                                                              
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE ::   ii__SV  ! WORK   point   i Coordinate     
!$OMP THREADPRIVATE(ii__SV)
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE ::   jj__SV  ! WORK   point   j Coordinate     
!$OMP THREADPRIVATE(jj__SV)
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE ::   nn__SV  ! WORK   point   n Coordinate     
!$OMP THREADPRIVATE(nn__SV)
                                                                               
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   IRu_SV  ! UPward    IR Flux (effective)   
!$OMP THREADPRIVATE(IRu_SV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   hSalSV  ! Saltating Layer Height          
!$OMP THREADPRIVATE(hSalSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   qSalSV  ! Saltating Snow  Concentration   
!$OMP THREADPRIVATE(qSalSV)
      REAL, DIMENSION(:),ALLOCATABLE,SAVE    ::   RnofSV  ! RunOFF    Intensity             
!$OMP THREADPRIVATE(RnofSV)
                  
CONTAINS



  SUBROUTINE INIT_VARxSV
  IMPLICIT NONE

      INTEGER                                :: ikl

      ALLOCATE(LSmask(klonv))  ! Land-Sea   Mask                 
      ALLOCATE(isotSV(klonv))  ! Soil       Type                 
      ALLOCATE(iWaFSV(klonv))  ! Soil       Drainage:(1,0)=(y,n) 
      ALLOCATE(ivgtSV(klonv))  ! Vegetation Type                 
                                                                              
      ALLOCATE(coszSV(klonv))  ! Cosine of Sun zenithal Angle    
      ALLOCATE(sol_SV(klonv))  ! Downward  Solar    Radiation    
      ALLOCATE(IRd_SV(klonv))  ! Downward  Longwave Radiation    
                                                                                
      ALLOCATE(drr_SV(klonv))  ! Rain  Intensity       [kg/m2/s] 
      ALLOCATE(dsn_SV(klonv))  ! Snow  Intensity       [kg/m2/s] 
      ALLOCATE(dsnbSV(klonv))  ! Idem, from Drift        [kg/m2] 
      ALLOCATE(esnbSV(klonv))  ! Idem, from Drift        [kg/m2] 
      ALLOCATE(dbs_SV(klonv))  ! Drift Amount            [kg/m2] 
      ALLOCATE(BrosSV(klonv))
      ALLOCATE(BG1sSV(klonv))          
      ALLOCATE(BG2sSV(klonv))
      ALLOCATE(dz0_SV(klonv))  ! dz0(Sastrugi dh)            [m] 
                                                                               
      ALLOCATE(cld_SV(klonv))  ! Cloudiness (seen from SBL)      
      ALLOCATE(za__SV(klonv))  ! SBL Height                      
      ALLOCATE(VV__SV(klonv))  !(SBL Top)  Wind Velocity         
      ALLOCATE(VVs_SV(klonv))
      ALLOCATE(RRs_SV(klonv))
      ALLOCATE(DDs_SV(klonv))
      ALLOCATE(TaT_SV(klonv))  ! SBL Top   Temperature           
      ALLOCATE(ExnrSV(klonv))  ! Exner     Potential             
      ALLOCATE(dSdTSV(klonv))  ! Sensible Heat Flux T Derivat.   
      ALLOCATE(dLdTSV(klonv))  ! Latent   Heat Flux T Derivat.   
      ALLOCATE(rhT_SV(klonv))  ! SBL Top   Air  Density          
      ALLOCATE(QaT_SV(klonv))  ! SBL Top   Specific Humidity     
      ALLOCATE(dQa_SV(klonv))  ! SBL Flux  Limitation of Qa      
      ALLOCATE(qsnoSV(klonv))  ! SBL Mean  Snow       Content    
                                                                            
      ALLOCATE(LAI0SV(klonv))  ! Nominal Leaf Area Index         
      ALLOCATE(glf0SV(klonv))  ! Green   Leaf Fraction           
                                                                               
      ALLOCATE(alb0SV(klonv))  ! Soil    Albedo                  
      ALLOCATE(slopSV(klonv))  ! Snow/Ice/Soil-Water Surf. Slope 

                                                                               
                                                                               
! +--SISVAT INPUT/OUTPUT Variables                                             
! +  -----------------------------                                             
                                                                               
      ALLOCATE(isnoSV(klonv))  ! Nb of Ice/Snow Layers           
      ALLOCATE(ispiSV(klonv))  ! Uppermost superimposed ice      
      ALLOCATE(iiceSV(klonv))  ! Nb of Ice      Layers           
      ALLOCATE(istoSV(klonv,0:nsno))  ! Snow Layer     History          
                                                                               
      ALLOCATE(alb_SV(klonv))  ! Surface-Canopy Albedo           
      ALLOCATE(emi_SV(klonv))  ! Surface-Canopy Emissivity       
      ALLOCATE(IRs_SV(klonv))  ! Soil           IR Flux          
      ALLOCATE(LMO_SV(klonv))  ! Monin-Obukhov  Scale            
      ALLOCATE(us__SV(klonv))  ! Friction       Velocity         
      ALLOCATE(uts_SV(klonv))  ! Temperature Turbulent Scale     
      ALLOCATE(cutsSV(klonv))  ! Temperature Turbulent Scale C.  
      ALLOCATE(uqs_SV(klonv))  ! Spec.Humid. Turbulent Scale     
      ALLOCATE(uss_SV(klonv))  ! Blow.Snow   Turbulent Scale     
      ALLOCATE(usthSV(klonv))  ! Blowing Snow Erosion Thresh.    
      ALLOCATE(rCDmSV(klonv))  ! Square  Root Contribut. Drag_m  
      ALLOCATE(rCDhSV(klonv))  ! Square  Root Contribut. Drag_h  
      ALLOCATE(Z0m_SV(klonv))  ! Momentum     Roughness Length   
      ALLOCATE(Z0mmSV(klonv))  !  z0(Momentum,    Time Mean) [m] 
      ALLOCATE(Z0mnSV(klonv))  !  z0(Momentum,    instanta.) [m] 
      ALLOCATE(Z0roSV(klonv))  ! Subgrid Topo Roughness Length   
      ALLOCATE(Z0SaSV(klonv))  !  z0(Sastrugi  h)            [m] 
      ALLOCATE(Z0e_SV(klonv))  !  z0(Snow eroded)            [m] 
      ALLOCATE(Z0emSV(klonv))  !  z0(Snow eroded, Time Mean) [m] 
      ALLOCATE(Z0enSV(klonv))  !  z0(Snow eroded, instanta.) [m] 
      ALLOCATE(Z0h_SV(klonv))  ! Heat         Roughness Length   
      ALLOCATE(Z0hmSV(klonv))  !  z0(Heat,        Time Mean) [m] 
      ALLOCATE(Z0hnSV(klonv))  !  z0(Heat,        instanta.) [m] 
                                                                          
      ALLOCATE(snCaSV(klonv))  ! Canopy  Snow   Thickness        
      ALLOCATE(rrCaSV(klonv))  ! Canopy  Water  Content          
      ALLOCATE(psivSV(klonv))  ! Leaf    Water  Potential        
      ALLOCATE(TvegSV(klonv)) ! Vegetation     Temperature      
                                                                        
      ALLOCATE(TsisSV(klonv,-nsol:nsno))  ! Snow/Ice/Soil-Water Temperature 
      ALLOCATE(ro__SV(klonv,-nsol:nsno))  ! Snow/Ice/Soil-Water VolumicMass 
      ALLOCATE(eta_SV(klonv,-nsol:nsno))  ! Snow/Ice/Soil     Water Content 
      ALLOCATE(G1snSV(klonv,    0:nsno))  ! Snow Dendricity/Sphericity      
      ALLOCATE(G2snSV(klonv,    0:nsno))  ! Snow Sphericity/Size            
      ALLOCATE(dzsnSV(klonv,    0:nsno))  ! Snow Layer  Thickness           
      ALLOCATE(agsnSV(klonv,    0:nsno))  ! Snow Age                        
      ALLOCATE(BufsSV(klonv))  ! Snow Buffer Layer               
      ALLOCATE(rusnSV(klonv))  ! Surficial   Water               
      ALLOCATE(SWf_SV(klonv))  ! Normalized  Decay               
      ALLOCATE(SWS_SV(klonv))  ! Surficial Water Status          
      ALLOCATE(HFraSV(klonv))  ! Frazil      Thickness           
                                                                            
      ALLOCATE(zWE_SV(klonv))  ! Current   Snow Thickness [mmWE] 
      ALLOCATE(zWEcSV(klonv))  ! Compacted Snow Thickness [mmWE] 
      ALLOCATE(wem_SV(klonv)) ! Only Melting             [mmWE] 
      ALLOCATE(wer_SV(klonv)) ! Refreezing               [mmWE] 
      ALLOCATE(wes_SV(klonv))  ! Sublimation              [mmWE] 
                                                                             
                                                                              
! +--SISVAT OUTPUT       Variables                                            
! +  -----------------------------                                            
                                                                               
!                                                                              
      ALLOCATE(ii__SV(klonv))  ! WORK   point   i Coordinate     
      ALLOCATE(jj__SV(klonv))  ! WORK   point   j Coordinate     
      ALLOCATE(nn__SV(klonv))  ! WORK   point   n Coordinate     
                                                                               
      ALLOCATE(IRu_SV(klonv))  ! UPward    IR Flux (effective)   
      ALLOCATE(hSalSV(klonv))  ! Saltating Layer Height          
      ALLOCATE(qSalSV(klonv))  ! Saltating Snow  Concentration   
      ALLOCATE(RnofSV(klonv))  ! RunOFF    Intensity      

      DO ikl=1,klonv        
        LSmask(ikl)   = 0
        isotSV(ikl)   = 0                
        iWaFSV(ikl)   = 0 
        ivgtSV(ikl)   = 0
        isnoSV(ikl)   = 0  
        ispiSV(ikl)   = 0
        iiceSV(ikl)   = 0          
        istoSV(ikl,:) = 0
        ii__SV(ikl)   = 0
        jj__SV(ikl)   = 0
        nn__SV(ikl)   = 0 
      END DO
  END SUBROUTINE INIT_VARxSV


                                               
END MODULE VARxSV
