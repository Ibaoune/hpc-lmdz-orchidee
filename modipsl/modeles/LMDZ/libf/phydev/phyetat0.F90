!
! $Id $
!
SUBROUTINE phyetat0(fichnom)
! Load initial state for the physics
! and do some resulting initializations

  USE dimphy, only: klon
  USE iostart, ONLY : open_startphy,get_field,close_startphy
  USE iophy, ONLY : init_iophy_new
  USE geometry_mod, ONLY : longitude_deg, latitude_deg

  IMPLICIT NONE

  CHARACTER(len=*),INTENT(in) :: fichnom ! input file name

  REAL :: lon_startphy(klon), lat_startphy(klon)
  INTEGER :: i

  ! open physics initial state file:
  CALL open_startphy(fichnom)

  ! read latitudes and make a sanity check (because already known from dyn)
  CALL get_field("latitude",lat_startphy)
  DO i=1,klon
    IF (ABS(lat_startphy(i)-latitude_deg(i))>=1) THEN
      WRITE(*,*) "phyetat0: Error! Latitude discrepancy wrt startphy file:",&
                 " i=",i," lat_startphy(i)=",lat_startphy(i),&
                 " latitude_deg(i)=",latitude_deg(i)
      ! This is presumably serious enough to abort run
      CALL abort_physic("phyetat0","discrepancy in latitudes!",1)
    ENDIF
    IF (ABS(lat_startphy(i)-latitude_deg(i))>=0.0001) THEN
      WRITE(*,*) "phyetat0: Warning! Latitude discrepancy wrt startphy file:",&
                 " i=",i," lat_startphy(i)=",lat_startphy(i),&
                 " latitude_deg(i)=",latitude_deg(i)
    ENDIF
  ENDDO

  ! read longitudes and make a sanity check (because already known from dyn)
  CALL get_field("longitude",lon_startphy)
  DO i=1,klon
    IF (ABS(lon_startphy(i)-longitude_deg(i))>=1) THEN
      WRITE(*,*) "phyetat0: Error! Longitude discrepancy wrt startphy file:",&
                 " i=",i," lon_startphy(i)=",lon_startphy(i),&
                 " longitude_deg(i)=",longitude_deg(i)
      ! This is presumably serious enough to abort run
      CALL abort_physic("phyetat0","discrepancy in longitudes!",1)
    ENDIF
    IF (ABS(lon_startphy(i)-longitude_deg(i))>=0.0001) THEN
      WRITE(*,*) "phyetat0: Warning! Longitude discrepancy wrt startphy file:",&
                 " i=",i," lon_startphy(i)=",lon_startphy(i),&
                 " longitude_deg(i)=",longitude_deg(i)
    ENDIF
  ENDDO

  ! read in other variables here ...

  ! close file
  CALL close_startphy

  ! do some more initializations
  CALL init_iophy_new(latitude_deg,longitude_deg)

END SUBROUTINE phyetat0
