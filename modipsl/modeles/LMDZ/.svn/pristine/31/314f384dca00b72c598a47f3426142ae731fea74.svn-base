SUBROUTINE read_rsun_rrtm(debut)

!****************************************************************************************
! This routine will read the solar constant fraction per band
!
! Olivier Boucher with inputs from Marion Marchand 
!****************************************************************************************

  USE netcdf95, ONLY: nf95_close, nf95_inq_varid, nf95_open, nf95_gw_var
  USE netcdf, ONLY: nf90_get_var, nf90_noerr, nf90_nowrite

  USE phys_cal_mod, ONLY : days_elapsed, year_len

  USE mod_phys_lmdz_mpi_data, ONLY: is_mpi_root
  USE mod_phys_lmdz_omp_data, ONLY: is_omp_root
  USE mod_phys_lmdz_para

  USE YOESW, ONLY : RSUN

  IMPLICIT NONE

  INCLUDE "clesphys.h"

  ! Input arguments
  LOGICAL, INTENT(IN) :: debut

! Local variables
  INTEGER :: ncid, dimid, varid, ncerr, nbday
  REAL, POINTER :: wlen(:), time(:)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: SSI_FRAC
!$OMP THREADPRIVATE(SSI_FRAC)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: TSI(:)
!$OMP THREADPRIVATE(TSI)

  INTEGER, SAVE :: day_pre=1
!$OMP THREADPRIVATE(day_pre)

!--only one processor reads
    IF (debut) THEN 

    ALLOCATE(SSI_FRAC(NSW,year_len))
    ALLOCATE(TSI(year_len))

    IF (is_mpi_root.AND.is_omp_root) THEN

       CALL nf95_open('solarforcing.nc', NF90_NOWRITE, ncid)

       CALL nf95_inq_varid(ncid, 'wlen', varid)
       CALL nf95_gw_var(ncid, varid, wlen)

       CALL nf95_inq_varid(ncid, 'time', varid)
       CALL nf95_gw_var(ncid, varid, time)

       IF (NSW.NE.size(wlen)) THEN 
         PRINT *,'read_rsun_rrtm NSW <> wlen = ',NSW, size(wlen)
         CALL abort_physic('read_rsun_rrtm','size of SSI is different from NSW',1)
       ENDIF

!--test if time is different from year_len but allow a mismatch of 1 day
       IF (size(time).NE.year_len.AND.size(time).NE.year_len+1) THEN
         PRINT *,'read_rsun_rrtm time <> year_len = ', size(time), year_len
         CALL abort_physic('read_rsun_rrtm','time dim should be the number of days in year',1)
       ENDIF
!--warning only if forcing file has 366 days but year_len has only 365
       IF (size(time).EQ.year_len+1) THEN 
         PRINT *,'Warning read_rsun_rrtm uses a leap year rsun for a noleap year'
       ENDIF

       CALL nf95_inq_varid(ncid, 'ssi_frac', varid)
       ncerr = nf90_get_var(ncid, varid, SSI_FRAC)

       CALL nf95_inq_varid(ncid, 'tsi', varid)
       ncerr = nf90_get_var(ncid, varid, TSI)

       CALL nf95_close(ncid)

       DO nbday=1, year_len
         IF (ABS(SUM(SSI_FRAC(:,nbday))-1.).GT.1.e-6) THEN 
           PRINT *,'somme SSI_FRAC=', SUM(SSI_FRAC(:,nbday))
           CALL abort_physic('read_rsun_rrtm','somme SSI_FRAC <> 1',1)
         ENDIF
       ENDDO
     
    ENDIF ! is_mpi_root .AND. is_omp_root

!$OMP BARRIER
    CALL bcast(SSI_FRAC)
    CALL bcast(TSI)

    ENDIF

!--only read at beginning of day
!--day in year is provided as days_elapsed since the beginning of the year +1
    IF (debut.OR.days_elapsed+1.NE.day_pre) THEN

!--keep memory of previous day
      day_pre=days_elapsed+1

!--copy 
      RSUN(1:NSW)=SSI_FRAC(:,days_elapsed+1)
      solaire=TSI(days_elapsed+1)

      print *,'READ_RSUN_RRTM day=', days_elapsed+1,' solaire=', solaire, ' RSUN=', RSUN(1:NSW)

    ENDIF !--fin allocation

END SUBROUTINE read_rsun_rrtm
