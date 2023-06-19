!
! $Id $
!
subroutine phyetat0(fichnom)
! Load initial state for the physics
! and do some resulting initializations

use iostart, only : open_startphy,get_field,close_startphy
use iophy, only : init_iophy_new
use phys_state_var_mod, only : rlat,rlon

implicit none

character(len=*),intent(in) :: fichnom ! input file name

! open physics initial state file:
call open_startphy(fichnom)

! read latitudes
call get_field("latitude",rlat)

! read longitudes
call get_field("longitude",rlon)

! read in other variables here ...

! close file
call close_startphy

! do some more initializations
call init_iophy_new(rlat,rlon)

end subroutine phyetat0
