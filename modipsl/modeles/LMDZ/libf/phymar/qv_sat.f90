      function qv_sat(TTK,ss,pstar,pt,lsf)

!--------------------------------------------------------------------------+
!   MAR PHYSICS                                         Mc 30-05-2007  MAR |
!     Function qv_sat computes the Saturation Specific Humidity    (kg/kg) |
!                                                                          |
!--------------------------------------------------------------------------+
!                                                                          |
!     INPUT :   TTK            : Air Temperature                       (K) |
!     ^^^^^^^   pstar * ss + pt: Pressure of sigma level ss          (kPa) |
!                                                                          |
!     OUTPUT :  e__sat: Saturation Vapor    Pressure                 (hPa) |
!     ^^^^^^^   qv_sat: Saturation Specific Humidity               (kg/kg) |
!                                                                          |
!--------------------------------------------------------------------------+


      IMPLICIT NONE


! Global Variables
! ================

      real     ::  qv_sat
      real     ::  TTK   
      real     ::  ss    
      real     ::  pstar 
      real     ::  pt    
      integer  ::  lsf


! Local  Variables
! ================

      real     ::  pr    
      real     ::  e__sat
      real     ::  r273p1 = 273.16
      real     ::  zer0   =   1.00
      real     ::  un_1   =   1.00
      real     ::  eps9   =   1.e-9
      real     ::  pr__75 =  75.00
      real     ::  pr_b75


! Saturation Vapor    Pressure
! ============================


      pr   = 10.d0 *(pstar *ss + pt)                                    !  pressure (hPa)

      IF (TTK.ge.273.16d0.or.lsf.eq.0)                              THEN

        e__sat =  6.1078d0 * exp (5.138d0*log(  r273p1     /TTK))      &!  saturated vapor pressure with respect to water
     &                     * exp (6827.d0*(un_1/r273p1-un_1/TTK))       !  Dudhia (1989) JAS, (B1) and (B2) p.3103
                                                                        !  See also Pielke (1984), p.234 and Stull (1988), p.276 

      ELSE
        e__sat =  6.107d0  * exp (6150.d0*(un_1/r273p1-un_1/TTK))       !  saturated vapor pressure with respect to ice
                                                                        !  Dudhia (1989) JAS, (B1) and (B2) p.3103
      END IF

        pr_b75 = max(zer0  ,           sign(un_1,pr - pr__75))



!       ******
        qv_sat = pr_b75*max(eps9  , .622d0*e__sat/(pr-.378d0*e__sat))  &!
     &    + (1.0-pr_b75)  * 1.e-3                                       !
!       ******



      return
      end      
