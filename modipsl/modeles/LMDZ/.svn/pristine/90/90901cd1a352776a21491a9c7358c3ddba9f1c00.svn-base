! $Id$
MODULE open_climoz_m

  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE

CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE open_climoz(ncid, press_in_cen, press_in_edg, time_in, daily, adjust)
!
!-------------------------------------------------------------------------------
  USE netcdf95, ONLY: nf95_open, nf95_close, nf95_gw_var, nf95_inq_varid
  USE netcdf,   ONLY: nf90_nowrite
  USE mod_phys_lmdz_mpi_data,      ONLY: is_mpi_root
  USE mod_phys_lmdz_mpi_transfert, ONLY: bcast_mpi
  USE phys_cal_mod,                ONLY: calend, year_len, year_cur
!-------------------------------------------------------------------------------
! Purpose: This procedure should be called once per "gcm" run, by a single
!          thread of each MPI process.
!          The root MPI process opens "climoz_LMDZ.nc", reads the pressure
! levels and the times and broadcasts them to the other processes.
!          We assume that, in "climoz_LMDZ.nc", the pressure levels are in hPa
! and in strictly ascending order.
!-------------------------------------------------------------------------------
! Arguments (OpenMP shared):
  INTEGER, INTENT(OUT):: ncid      !--- "climoz_LMDZ.nc" identifier
  REAL, POINTER :: press_in_cen(:) !--- at cells centers
  REAL, POINTER :: press_in_edg(:) !--- at the interfaces (pressure intervals)
  REAL, POINTER :: time_in(:)      !--- records times, in days since Jan. 1st
  LOGICAL, INTENT(IN) :: daily     !--- daily files (calendar dependent days nb)
  LOGICAL, INTENT(IN) :: adjust    !--- tropopause adjustement required
!   pressure levels press_in_cen/edg are in Pa a,d strictly ascending order.
!   time_in is only used for monthly files (14 records)
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: varid                 !--- NetCDF variables identifier
  INTEGER :: nlev, ntim            !--- pressure levels and time records numbers
  CHARACTER(LEN=80)  :: sub
  CHARACTER(LEN=320) :: msg
!-------------------------------------------------------------------------------
  sub="open_climoz"
  WRITE(lunout,*)"Entering routine "//TRIM(sub)

  IF(is_mpi_root) THEN

    !--- OPEN FILE, READ PRESSURE LEVELS AND TIME VECTOR
    CALL nf95_open("climoz_LMDZ.nc", nf90_nowrite, ncid)
    CALL nf95_inq_varid(ncid, "plev", varid)
    CALL nf95_gw_var(ncid, varid, press_in_cen)
    ! Convert from hPa to Pa because "paprs" and "pplay" are in Pa:
    press_in_cen = press_in_cen * 100.
    nlev = SIZE(press_in_cen)
    CALL NF95_INQ_VARID(ncID, "time", varID)
    CALL NF95_GW_VAR(ncid, varid, time_in)
    ntim = SIZE(time_in)

    !--- BUILD EDGES OF PRESSURE INTERVALS: HALFWAY IN LOGARITHMS
    ALLOCATE(press_in_edg(nlev+1))
    press_in_edg=[0.,SQRT(press_in_cen(1:nlev-1)*press_in_cen(2:nlev)),HUGE(0.)]

    !--- CHECK RECORDS NUMBER AND DISPLAY CORRESPONDING INFORMATION
    IF(daily.AND.ntim/=year_len) THEN
      WRITE(msg,'(a,3(i4,a))')TRIM(sub)//': Expecting a daily ozone file with',&
     &year_len,' records (year ',year_cur,') ; found ',ntim,' instead'
      CALL abort_physic(sub, msg, 1)
    ELSE IF(ALL([360,14]/=ntim)) THEN
      WRITE(msg,'(a,i4,a)')TRIM(sub)//': Expecting an ozone file with 14 (mont'&
     &//'hly case) or 360 (old style files) records ; found ',ntim,' instead'
      CALL abort_physic(sub, msg, 1)
    ELSE
      IF(daily) THEN
         WRITE(msg,'(a,2(i4,a))')'daily file (',ntim,' days in ',year_cur,')'
       ELSE IF(ntim==14) THEN
         msg='14 records monthly file'
       ELSE
         msg='360 records files (old convention)'
       END IF
      WRITE(lunout,*)TRIM(sub)//': Using a '//TRIM(msg)
    END IF

    !--- MESSAGE ABOUT OPTIONAL STRETCHING FOR TROPOPAUSES MATCHING
    IF(adjust) THEN
      WRITE(lunout,*)TRIM(sub)//': Adjusting O3 field to match gcm tropopause.'
    ELSE
      WRITE(lunout,*)TRIM(sub)//': Interpolating O3 field directly on gcm levels.'
    END IF

  END IF
  CALL bcast_mpi(nlev)
  IF(.NOT.is_mpi_root) ALLOCATE(press_in_cen(nlev  )); CALL bcast_mpi(press_in_cen)
  IF(.NOT.is_mpi_root) ALLOCATE(press_in_edg(nlev+1)); CALL bcast_mpi(press_in_edg)
  CALL bcast_mpi(ntim)
  IF(.NOT.is_mpi_root) ALLOCATE(time_in(ntim));        CALL bcast_mpi(time_in)

END SUBROUTINE open_climoz
!
!-------------------------------------------------------------------------------

END MODULE open_climoz_m
