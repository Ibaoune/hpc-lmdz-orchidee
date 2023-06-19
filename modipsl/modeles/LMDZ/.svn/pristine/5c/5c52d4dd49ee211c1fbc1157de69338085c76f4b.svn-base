!
! $Id: $
!
SUBROUTINE phyredem (fichnom)

  USE geometry_mod, ONLY : longitude_deg, latitude_deg
  USE iostart, ONLY: open_restartphy, close_restartphy, enddef_restartphy, put_field, put_var


  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: fichnom

  INTEGER,PARAMETER :: tab_cntrl_len=100
  REAL :: tab_cntrl(tab_cntrl_len)
  INTEGER :: pass

  ! open file

  CALL open_restartphy(fichnom)

  ! tab_cntrl() contains run parameters

  tab_cntrl(:)=0.0

  DO pass=1,2   ! pass=1 netcdf definition ; pass=2 netcdf write

    CALL put_var(pass, "controle", "Control parameters", tab_cntrl)

  ! coordinates

    CALL put_field(pass, "longitude", "Longitudes on physics grid", longitude_deg)
     
    CALL put_field(pass, "latitude", "Latitudes on physics grid", latitude_deg)

    IF (pass==1) CALL enddef_restartphy
    IF (pass==2) CALL close_restartphy

  ENDDO
  ! close file

  CALL close_restartphy
  !$OMP BARRIER

END SUBROUTINE phyredem
