MODULE LMDZCTR

IMPLICIT NONE
! +                                                                
                                                                                
      logical           reaVAR    ! Input INI: Prev.Dyn.Simulat. (MAR, GCM)     
      logical           reaLBC    ! Input LBC: Prev.Dyn.Simulat. (MAR, GCM)     
      logical           safVAR    ! Full Output on Saving Files   MARxxx.DAT    
      logical           sALONE    !                                             
      logical           geoNST    !                                             
      logical           conmas    ! Mass       Conserv. Constraint Init.Switch  
      logical           potvor    ! P.-Vortic. Conserv. Constraint Init.Switch  
      logical           hamfil    ! Initial Filtered Fields (Time, Hamming)     
      logical           brocam    ! Brown and Campana  Time      Scheme Switch  
      logical           LFrBAK    ! Leap-Frog Backward Advection Scheme Switch  
      logical           openLB    !                                             
      logical           sommlb    !                                             
      logical           FirstC    !                                             
      logical           qqmass    ! Mass       Conservation             Switch  
      logical           turhor    ! Horizontal Diffusion (Smagorinsky)  Switch  
      logical           SBLitr    !                                             
      logical           tur_25    !                                             
      logical           convec    ! Convective  Adjustment              Switch  
      logical           MFLX_d    ! Convective  Adjustment (deep)       Switch  
      logical           MFLX_s    ! Convective  Adjustment (shallow)    Switch  
      logical           micphy    ! Cloud       Microphysics            Switch  
      logical           fracld    ! Fractional  Cloudiness              Switch  
      logical           chimod    ! Atmospheric Chemical Model          Switch  
      logical           physic    ! Atmospheric/Surface Physics         Switch  
      logical           polmod    ! Interactive Polynya    is           Switch  
      logical           snomod    ! Interactive Snow Model is           Switch  
      logical           BloMod    !                                             
      logical           vegmod    ! Interactive SVAT                    Switch  
      logical           VSISVAT   !                                             
      logical           no_vec    ! Scalar     (NO Vectorization)       Switch  
                                                                                
                                                                                
      integer           itexpe,iterun,nterun,nbhour,itConv                      
      integer           iboucl,nboucl,nprint,ipr_nc,npr_nc                      
      integer           maptyp                                                  
      integer           log_1D                                                  
                                                                                
      real*8            Robert                                                  
      real*8            rrmin ,rrmax                                            
      real*8            rxbase,rxfact                                           
      real*8            fxlead                                                  
      real*8            tMFLXd        ! d(time) between 2 deep convection CALL  
      real*8            aMFLXd,aMFLXs ! adjustment times (deep, shallow)        
                                                                               
      character*16      fnam                                                    


END MODULE LMDZCTR

