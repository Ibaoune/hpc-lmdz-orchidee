MODULE cal_tools_m

  USE ioipsl, ONLY: ioconf_calendar, ioget_mon_len, lock_calendar,             &
                     ioget_calendar, ioget_year_len

CONTAINS

!-------------------------------------------------------------------------------
!
FUNCTION year_len(y,cal_in)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER                                :: year_len
  INTEGER,                    INTENT(IN) :: y
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: cal_in
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=20) :: cal_lmdz
!-------------------------------------------------------------------------------
  !--- No specified calendar: we need lmdz days number for current year
  IF(.NOT.PRESENT(cal_in)) THEN; year_len=ioget_year_len(y); RETURN; END IF

  !--- Get the lmdz calendar to reset at the end of the function
  CALL ioget_calendar(cal_lmdz)

  !--- Unlock calendar and set it to wanted one
  CALL lock_calendar(.FALSE.); CALL ioconf_calendar(TRIM(cal_in))

  !--- Get the number of days in this year
  year_len=ioget_year_len(y)

  !--- Back to original calendar
  CALL lock_calendar(.FALSE.); CALL ioconf_calendar(TRIM(cal_lmdz))

END FUNCTION year_len
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
FUNCTION mid_month(y,cal_in)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,                    INTENT(IN) :: y           ! year
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: cal_in      ! calendar
  REAL,                    DIMENSION(14) :: mid_month   ! mid-bins times
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=20)      :: cal_lmdz          ! lmdz current calendar
  INTEGER, DIMENSION(14) :: tlen              ! months lengths (days)
  INTEGER                :: m                 ! months counter
  INTEGER                :: nd                ! number of days
!-------------------------------------------------------------------------------
  IF(PRESENT(cal_in)) THEN
    CALL ioget_calendar(cal_lmdz)             !--- Keep track of lmdz calendar
    CALL lock_calendar(.FALSE.)               !--- Unlock calendar
    CALL ioconf_calendar(TRIM(cal_in))        !--- Change calendar to "cal_in"
  END IF

  !--- Get the length of each month
  tlen(1 )=ioget_mon_len(y-1,12)
  DO m=1,12; tlen(m+1)=ioget_mon_len(y,m); END DO
  tlen(14)=ioget_mon_len(y+1, 1)

  !--- Mid-bins times
  mid_month(1)=-0.5*REAL(tlen(1))
  DO m=2,14; mid_month(m)=mid_month(m-1)+0.5*REAL(tlen(m-1)+tlen(m)); END DO

  IF(PRESENT(cal_in)) THEN
    CALL lock_calendar(.FALSE.)               !--- Unlock calendar
    CALL ioconf_calendar(TRIM(cal_lmdz))      !--- Restore original calendar
  END IF

END FUNCTION mid_month
!
!-------------------------------------------------------------------------------


END MODULE cal_tools_m
!
!*******************************************************************************

