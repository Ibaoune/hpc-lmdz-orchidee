MODULE VAR0SV

USE VAR_SV, only : nsol, nsno 
USE VARdSV, only : nsot, nvgt

IMPLICIT NONE
                                                  
      INTEGER    :: islpSV(-nsol:0)                                           
      INTEGER    :: isnpSV(      nsno)                                        
      INTEGER    :: islmSV(-nsol:0)    
                                      
      INTEGER,PARAMETER :: nkhy=50                                            
      REAL       :: Implic,Explic                                             
      REAL       :: dzmiSV(-nsol:0)          ! dz_(i-1/2)                     
      REAL       :: dzi_SV(-nsol:0)          ! dz_(i-1)/(dz_(i)+dz_(i-1))     
      REAL       :: dziiSV(-nsol:0)          ! dz_(i)  /(dz_(i)+dz_(i-1))     
      REAL       :: dtz_SV(-nsol:0)          ! dt / dz                        
      REAL       :: dz78SV(-nsol:0)          ! 7/8 (dz)                       
      REAL       :: dz34SV(-nsol:0)          ! 3/4 (dz)                       
      REAL       :: dz_8SV(-nsol:0)          ! 1/8 (dz)                       
      REAL       :: dzAvSV(-nsol:0)          ! 1/8dz_(-1)+3/4dz+1/8dz_(+1)    
      REAL       :: OcndSV                   ! Swab Ocean / Soil Ratio        
      REAL       :: RF__SV( 0:nvgt,-nsol:0)  ! Root Fraction                  
      REAL       :: rocsSV( 0:nsot)          ! Soil Contribution to (ro c)_s  
      REAL       :: etamSV( 0:nsot)          ! Soil Minimum Humidity          
      REAL       :: s1__SV( 0:nsot)          ! ... X eta**( b+2), DR97(3.36)  
      REAL       :: s2__SV( 0:nsot)          ! ... X eta**(2b+3), DR97(3.35)  
      REAL       :: aKdtSV( 0:nsot, 0:nkhy)  ! Khyd=a*eta+b: a * dt           
      REAL       :: bKdtSV( 0:nsot, 0:nkhy)  ! Khyd=a*eta+b: b * dt           

END MODULE VAR0SV
