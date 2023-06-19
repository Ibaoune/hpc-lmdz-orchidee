!
! $Id: fcg_racmo.h 2011-02-16 17:02:56Z cheruy $
!
      real   :: a_guide 
      logical :: ok_invertp 
      integer :: forc_trb
      character*31 :: fich_racmo

      common /fcg_racmo/forc_trb,ok_invertp,a_guide,fich_racmo

!$OMP THREADPRIVATE(/fcg_racmo/)


































