! $Id: phys_cal_mod.F90 3452 2019-01-28 09:45:20Z fairhead $
MODULE phys_cal_mod
! This module contains information on the calendar at the current time step

  INTEGER,SAVE :: year_cur      ! current year
!$OMP THREADPRIVATE(year_cur)
  INTEGER,SAVE :: mth_cur       ! current month
!$OMP THREADPRIVATE(mth_cur)
  INTEGER,SAVE :: day_cur       ! current day
!$OMP THREADPRIVATE(day_cur)
  INTEGER,SAVE :: days_elapsed  ! number of whole days since start of the current year
!$OMP THREADPRIVATE(days_elapsed)
  INTEGER,SAVE :: mth_len       ! number of days in the current month
!$OMP THREADPRIVATE(mth_len)
  INTEGER,SAVE :: year_len      ! number of days in the current year
!$OMP THREADPRIVATE(year_len)
  REAL,SAVE    :: hour         ! seconds elapsed (in the current day) since midnight
!$OMP THREADPRIVATE(hour)
  REAL,SAVE    :: jD_1jan
!$OMP THREADPRIVATE(jD_1jan)
  REAL,SAVE    :: jH_1jan
!$OMP THREADPRIVATE(jH_1jan)
  REAL,SAVE    :: xjour
!$OMP THREADPRIVATE(xjour)
  REAL,SAVE    :: jD_cur  ! jour courant a l'appel de la physique (jour julien)
!$OMP THREADPRIVATE(jD_cur)
  REAL,SAVE    :: jH_cur  ! heure courante a l'appel de la physique (jour julien)
!$OMP THREADPRIVATE(jH_cur)
  REAL,SAVE    :: jD_ref  ! jour du demarage de la simulation (jour julien)
!$OMP THREADPRIVATE(jD_ref)
 CHARACTER (len=10) :: calend ! type of calendar to use
                              ! (can be earth_360d, earth_365d or earth_366d)
!$OMP THREADPRIVATE(calend)

CONTAINS
  
  SUBROUTINE phys_cal_init(annee_ref,day_ref)

    USE IOIPSL, ONLY:  ymds2ju, ioconf_calendar
    USE mod_phys_lmdz_para, ONLY:  is_master,is_omp_master
    USE ioipsl_getin_p_mod, ONLY: getin_p

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: annee_ref
    INTEGER,INTENT(IN) :: day_ref

    ! Find out which type of calendar we are using
    calend = 'earth_360d' ! default
    CALL getin_p("calend",calend)

    IF (is_omp_master) THEN
      IF (calend == 'earth_360d') THEN
        CALL ioconf_calendar('360d')
      ELSE IF (calend == 'earth_365d') THEN
        CALL ioconf_calendar('noleap')
      ELSE IF (calend == 'earth_366d' .OR. calend == 'gregorian') THEN
        CALL ioconf_calendar('gregorian')
      ELSE
        CALL abort_physic('phys_cal_init','Mauvais choix de calendrier',1)
      ENDIF
    ENDIF
!$OMP BARRIER
     
    CALL ymds2ju(annee_ref, 1, day_ref, 0., jD_ref)
    jD_ref=INT(jD_ref)
  
  END SUBROUTINE  phys_cal_init

  SUBROUTINE phys_cal_update(julian_date)
    ! This subroutine updates the module saved variables.

    USE IOIPSL, only: ju2ymds, ymds2ju, ioget_mon_len, ioget_year_len
    IMPLICIT NONE
    REAL, INTENT(IN) :: julian_date

    jD_cur=INT(julian_date)
    jH_cur=julian_date-jD_cur
    
    CALL ju2ymds(jD_cur+jH_cur, year_cur, mth_cur, day_cur, hour)
    CALL ymds2ju(year_cur, 1, 1, 0., jD_1jan)
    
    jH_1jan = jD_1jan - int (jD_1jan)
    jD_1jan = int (jD_1jan) 
    xjour = jD_cur - jD_1jan
    days_elapsed = jD_cur - jD_1jan

    ! Get lenght of current month
    mth_len = ioget_mon_len(year_cur,mth_cur)

    ! Get length of current year
    year_len = ioget_year_len(year_cur)

  END SUBROUTINE phys_cal_update

END MODULE phys_cal_mod
