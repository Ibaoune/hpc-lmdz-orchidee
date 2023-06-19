MODULE create_limit_unstruct_mod
    PRIVATE
    INTEGER, PARAMETER                             :: lmdep=12

    PUBLIC create_limit_unstruct

CONTAINS


  SUBROUTINE create_limit_unstruct
  USE dimphy
#ifdef CPP_XIOS
  USE xios
  USE ioipsl,             ONLY : ioget_year_len
  USE time_phylmdz_mod, ONLY : annee_ref
  USE indice_sol_mod
  USE phys_state_var_mod
  USE mod_phys_lmdz_para
  IMPLICIT NONE
    INCLUDE "iniprint.h"
    REAL,    DIMENSION(:,:),ALLOCATABLE            :: sic
    REAL,    DIMENSION(:,:),ALLOCATABLE            :: sst
    REAL,    DIMENSION(klon,lmdep)                 :: rugos
    REAL,    DIMENSION(klon,lmdep)                 :: albedo
    REAL,    DIMENSION(:,:),ALLOCATABLE            :: sic_mpi
    REAL,    DIMENSION(:,:),ALLOCATABLE            :: sst_mpi
    REAL,    DIMENSION(klon_mpi,lmdep)             :: rugos_mpi
    REAL,    DIMENSION(klon_mpi,lmdep)             :: albedo_mpi
    INTEGER                                        :: ndays
    REAL                                           :: fi_ice(klon)
    REAL, ALLOCATABLE                              :: sic_year(:,:)
    REAL, ALLOCATABLE                              :: sst_year(:,:)
    REAL, ALLOCATABLE                              :: rugos_year(:,:)
    REAL, ALLOCATABLE                              :: albedo_year(:,:)
    REAL, ALLOCATABLE                              :: pctsrf_t(:,:,:)
    REAL, ALLOCATABLE                              :: phy_bil(:,:)
    REAL, ALLOCATABLE                              :: sst_year_mpi(:,:)
    REAL, ALLOCATABLE                              :: rugos_year_mpi(:,:)
    REAL, ALLOCATABLE                              :: albedo_year_mpi(:,:)
    REAL, ALLOCATABLE                              :: pctsrf_t_mpi(:,:,:)
    REAL, ALLOCATABLE                              :: phy_bil_mpi(:,:)
    INTEGER :: l,k
    INTEGER :: nbad
    INTEGER :: sic_time_axis_size 
    INTEGER :: sst_time_axis_size
    CHARACTER(LEN=99)                  :: mess            ! error message
    
     
    ndays=ioget_year_len(annee_ref)
    
    IF (is_omp_master) CALL xios_get_axis_attr("time_sic",n_glo=sic_time_axis_size)
    CALL bcast_omp(sic_time_axis_size)
    ALLOCATE(sic_mpi(klon_mpi,sic_time_axis_size))
    ALLOCATE(sic(klon,sic_time_axis_size))
    
    
    IF (is_omp_master) CALL xios_get_axis_attr("time_sst",n_glo=sst_time_axis_size)
    CALL bcast_omp(sst_time_axis_size)
    ALLOCATE(sst_mpi(klon_mpi,sst_time_axis_size))
    ALLOCATE(sst(klon,sst_time_axis_size))
    
    IF (is_omp_master) THEN
      CALL xios_recv_field("sic_limit",sic_mpi)
      CALL xios_recv_field("sst_limit",sst_mpi)
      CALL xios_recv_field("rugos_limit",rugos_mpi)
      CALL xios_recv_field("albedo_limit",albedo_mpi)
    ENDIF
    CALL scatter_omp(sic_mpi,sic)
    CALL scatter_omp(sst_mpi,sst)
    CALL scatter_omp(rugos_mpi,rugos)
    CALL scatter_omp(albedo_mpi,albedo)
    
    ALLOCATE(sic_year(klon,ndays))
    ALLOCATE(sst_year(klon,ndays))
    ALLOCATE(rugos_year(klon,ndays))
    ALLOCATE(albedo_year(klon,ndays))
    ALLOCATE(pctsrf_t(klon,nbsrf,ndays))
    ALLOCATE(phy_bil(klon,ndays)); phy_bil=0.0


! sic
    IF (sic_time_axis_size==lmdep) THEN
      CALL time_interpolation(ndays,sic,'gregorian',sic_year)
    ELSE IF (sic_time_axis_size==ndays) THEN
      sic_year=sic
    ELSE
      WRITE(mess,*) 'sic time axis is nor montly, nor daily. sic time interpolation ',&
                    'is requiered but is not currently managed'
      CALL abort_physic('create_limit_unstruct',TRIM(mess),1)
    ENDIF
    
    sic_year(:,:)=sic_year(:,:)/100.  ! convert percent to fraction
    WHERE(sic_year(:,:)>1.0) sic_year(:,:)=1.0    ! Some fractions have some time large negative values
    WHERE(sic_year(:,:)<0.0) sic_year(:,:)=0.0    ! probably better to apply alse this filter before horizontal interpolation
    
! sst
    IF (sst_time_axis_size==lmdep) THEN
      CALL time_interpolation(ndays,sst,'gregorian',sst_year)
    ELSE IF (sst_time_axis_size==ndays) THEN
      sst_year=sst
    ELSE
      WRITE(mess,*)'sic time axis is nor montly, nor daily. sic time interpolation ',&
                   'is requiered but is not currently managed'
      CALL abort_physic('create_limit_unstruct',TRIM(mess),1)
    ENDIF
    WHERE(sst_year(:,:)<271.38) sst_year(:,:)=271.38


! rugos    
    DO l=1, lmdep
      WHERE(NINT(zmasq(:))/=1) rugos(:,l)=0.001
    ENDDO
    CALL time_interpolation(ndays,rugos,'360d',rugos_year)

! albedo    
    CALL time_interpolation(ndays,albedo,'360d',albedo_year)


    DO k=1,ndays
      fi_ice=sic_year(:,k)
      WHERE(fi_ice>=1.0  ) fi_ice=1.0
      WHERE(fi_ice<EPSFRA) fi_ice=0.0
      pctsrf_t(:,is_ter,k)=pctsrf(:,is_ter)       ! land soil
      pctsrf_t(:,is_lic,k)=pctsrf(:,is_lic)       ! land ice

!!     IF (icefile==trim(fcpldsic)) THEN           ! SIC=pICE*(1-LIC-TER)
!!        pctsrf_t(:,is_sic,k)=fi_ice(:)*(1.-pctsrf(:,is_lic)-pctsrf(:,is_ter))
!!     ELSE IF (icefile==trim(fhistsic)) THEN      ! SIC=pICE
!!        pctsrf_t(:,is_sic,k)=fi_ice(:)
!!     ELSE ! icefile==famipsic                    ! SIC=pICE-LIC
        pctsrf_t(:,is_sic,k)=fi_ice-pctsrf_t(:,is_lic,k)
!     END IF
      WHERE(pctsrf_t(:,is_sic,k)<=0) pctsrf_t(:,is_sic,k)=0.
      WHERE(1.0-zmasq<EPSFRA)
        pctsrf_t(:,is_sic,k)=0.0
        pctsrf_t(:,is_oce,k)=0.0
      ELSEWHERE
        WHERE(pctsrf_t(:,is_sic,k)>=1.0-zmasq)
          pctsrf_t(:,is_sic,k)=1.0-zmasq
          pctsrf_t(:,is_oce,k)=0.0
        ELSEWHERE
          pctsrf_t(:,is_oce,k)=1.0-zmasq-pctsrf_t(:,is_sic,k)
          WHERE(pctsrf_t(:,is_oce,k)<EPSFRA)
             pctsrf_t(:,is_oce,k)=0.0
             pctsrf_t(:,is_sic,k)=1.0-zmasq
          END WHERE
        END WHERE
      END WHERE
      nbad=COUNT(pctsrf_t(:,is_oce,k)<0.0)
      IF(nbad>0) WRITE(lunout,*) 'pb sous maille pour nb point = ',nbad
      nbad=COUNT(abs(sum(pctsrf_t(:,:,k),dim=2)-1.0)>EPSFRA)
      IF(nbad>0) WRITE(lunout,*) 'pb sous surface pour nb points = ',nbad
    END DO
    
    ALLOCATE(sst_year_mpi(klon_mpi,ndays))
    ALLOCATE(rugos_year_mpi(klon_mpi,ndays))
    ALLOCATE(albedo_year_mpi(klon_mpi,ndays))
    ALLOCATE(pctsrf_t_mpi(klon_mpi,nbsrf,ndays))
    ALLOCATE(phy_bil_mpi(klon_mpi,ndays))
    
    CALL gather_omp(pctsrf_t   , pctsrf_t_mpi)
    CALL gather_omp(sst_year   , sst_year_mpi)
    CALL gather_omp(phy_bil    , phy_bil_mpi)
    CALL gather_omp(albedo_year, albedo_year_mpi)
    CALL gather_omp(rugos_year , rugos_year_mpi)

    IF (is_omp_master) THEN
      CALL xios_send_field("foce_limout",pctsrf_t_mpi(:,is_oce,:))
      CALL xios_send_field("fsic_limout",pctsrf_t_mpi(:,is_sic,:))
      CALL xios_send_field("fter_limout",pctsrf_t_mpi(:,is_ter,:))
      CALL xios_send_field("flic_limout",pctsrf_t_mpi(:,is_lic,:))
      CALL xios_send_field("sst_limout", sst_year_mpi)
      CALL xios_send_field("bils_limout",phy_bil_mpi)
      CALL xios_send_field("alb_limout", albedo_year_mpi) 
      CALL xios_send_field("rug_limout", rugos_year_mpi) 
    ENDIF
#endif
  END SUBROUTINE create_limit_unstruct
  
  
  SUBROUTINE time_interpolation(ndays,field_in,calendar,field_out)
  USE pchsp_95_m, only: pchsp_95
  USE pchfe_95_m, only: pchfe_95
  USE arth_m, only: arth
  USE dimphy, ONLY : klon
  USE ioipsl,             ONLY : ioget_year_len
  USE time_phylmdz_mod, ONLY : annee_ref
  USE mod_phys_lmdz_para
  IMPLICIT NONE
   INCLUDE "iniprint.h"

   INTEGER,         INTENT(IN)  :: ndays
   REAL,            INTENT(IN)  :: field_in(klon,lmdep)
   CHARACTER(LEN=*),INTENT(IN)  :: calendar
   REAL,            INTENT(OUT) :: field_out(klon,ndays)
 
   INTEGER :: ndays_in
   REAL    :: timeyear(lmdep)   
   REAL    :: yder(lmdep)   
   INTEGER :: ij,ierr, n_extrap
   LOGICAL :: skip

   CHARACTER (len = 50)         :: modname = 'create_limit_unstruct.time_interpolation'
   CHARACTER (len = 80)         :: abort_message

  
   IF (is_omp_master) ndays_in=year_len(annee_ref, calendar)
   CALL bcast_omp(ndays_in)
   IF (is_omp_master) timeyear=mid_months(annee_ref, calendar, lmdep)
   CALL bcast_omp(timeyear)
    
   n_extrap = 0
   skip=.FALSE.
   DO ij=1,klon
     yder = pchsp_95(timeyear, field_in(ij, :), ibeg=2, iend=2, vc_beg=0., vc_end=0.)
     CALL pchfe_95(timeyear, field_in(ij, :), yder, skip, arth(0., real(ndays_in) / ndays, ndays), field_out(ij, :), ierr)
     if (ierr < 0) then
        abort_message='error in pchfe_95'
        CALL abort_physic(modname,abort_message,1)
     endif
     n_extrap = n_extrap + ierr
   END DO
   
   IF (n_extrap /= 0) then
     WRITE(lunout,*) "get_2Dfield pchfe_95: n_extrap = ", n_extrap
   ENDIF 
  
  
  END SUBROUTINE time_interpolation
  !-------------------------------------------------------------------------------
  !
  FUNCTION year_len(y,cal_in)
  !
  !-------------------------------------------------------------------------------
    USE ioipsl, ONLY : ioget_calendar,ioconf_calendar,lock_calendar,ioget_year_len
    IMPLICIT NONE
  !-------------------------------------------------------------------------------
  ! Arguments:
    INTEGER                       :: year_len
    INTEGER,           INTENT(IN) :: y
    CHARACTER(LEN=*),  INTENT(IN) :: cal_in
  !-------------------------------------------------------------------------------
  ! Local variables:
    CHARACTER(LEN=20)             :: cal_out              ! calendar (for outputs)
  !-------------------------------------------------------------------------------
  !--- Getting the input calendar to reset at the end of the function
    CALL ioget_calendar(cal_out)
  
  !--- Unlocking calendar and setting it to wanted one
    CALL lock_calendar(.FALSE.); CALL ioconf_calendar(TRIM(cal_in))
  
  !--- Getting the number of days in this year
    year_len=ioget_year_len(y)
  
  !--- Back to original calendar
    CALL lock_calendar(.FALSE.); CALL ioconf_calendar(TRIM(cal_out))
  
  END FUNCTION year_len
  !
  !-------------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------------
  !
  FUNCTION mid_months(y,cal_in,nm)
  !
  !-------------------------------------------------------------------------------
    USE ioipsl, ONLY : ioget_calendar,ioconf_calendar,lock_calendar,ioget_mon_len
    IMPLICIT NONE
  !-------------------------------------------------------------------------------
  ! Arguments:
    INTEGER,                INTENT(IN) :: y               ! year
    CHARACTER(LEN=*),       INTENT(IN) :: cal_in          ! calendar
    INTEGER,                INTENT(IN) :: nm              ! months/year number
    REAL,    DIMENSION(nm)             :: mid_months      ! mid-month times
  !-------------------------------------------------------------------------------
  ! Local variables:
    CHARACTER(LEN=99)                  :: mess            ! error message
    CHARACTER(LEN=20)                  :: cal_out         ! calendar (for outputs)
    INTEGER, DIMENSION(nm)             :: mnth            ! months lengths (days)
    INTEGER                            :: m               ! months counter
    INTEGER                            :: nd              ! number of days
    INTEGER                            :: k
  !-------------------------------------------------------------------------------
    nd=year_len(y,cal_in)
  
    IF(nm==12) THEN
  
    !--- Getting the input calendar to reset at the end of the function
      CALL ioget_calendar(cal_out)
  
    !--- Unlocking calendar and setting it to wanted one
      CALL lock_calendar(.FALSE.); CALL ioconf_calendar(TRIM(cal_in))
  
    !--- Getting the length of each month
      DO m=1,nm; mnth(m)=ioget_mon_len(y,m); END DO
  
    !--- Back to original calendar
      CALL lock_calendar(.FALSE.); CALL ioconf_calendar(TRIM(cal_out))
  
    ELSE IF(MODULO(nd,nm)/=0) THEN
      WRITE(mess,'(a,i3,a,i3,a)')'Unconsistent calendar: ',nd,' days/year, but ',&
        nm,' months/year. Months number should divide days number.'
      CALL abort_physic('mid_months',TRIM(mess),1)
  
    ELSE
      mnth=(/(m,m=1,nm,nd/nm)/)
    END IF
  
  !--- Mid-months times
    mid_months(1)=0.5*REAL(mnth(1))
    DO k=2,nm
      mid_months(k)=mid_months(k-1)+0.5*REAL(mnth(k-1)+mnth(k))
    END DO
  
  END FUNCTION mid_months
  

END MODULE create_limit_unstruct_mod
