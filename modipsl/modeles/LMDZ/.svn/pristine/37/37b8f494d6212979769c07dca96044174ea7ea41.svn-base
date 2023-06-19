! $Header$

SUBROUTINE limit_slab(itime, dtime, jour, lmt_bils, diff_sst, diff_siv)

  USE dimphy
  USE mod_grid_phy_lmdz, ONLY: klon_glo
  USE mod_phys_lmdz_para
  USE netcdf 
  USE indice_sol_mod
  USE ocean_slab_mod, ONLY: nslay

  IMPLICIT NONE

  INCLUDE "clesphys.h"

! In- and ouput arguments
!****************************************************************************************
  INTEGER, INTENT(IN) :: itime   ! numero du pas de temps courant
  INTEGER, INTENT(IN) :: jour    ! jour a lire dans l'annee
  REAL   , INTENT(IN) :: dtime   ! pas de temps de la physique (en s)
  REAL, DIMENSION(klon), INTENT(OUT) ::  diff_sst, diff_siv
  REAL, DIMENSION(klon,nslay), INTENT(OUT) :: lmt_bils

! Locals variables with attribute SAVE
!****************************************************************************************
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: diff_siv_save, diff_sst_save
  REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: bils_save
!$OMP THREADPRIVATE(bils_save, diff_sst_save, diff_siv_save)

! Locals variables
!****************************************************************************************
  INTEGER                  :: lmt_pas   
  INTEGER                  :: nvarid, nid, ierr, i
  INTEGER, DIMENSION(2)    :: start, epais 
  REAL, DIMENSION(klon_glo):: sst_l_glo, sst_lp1_glo, diff_sst_glo
  REAL, DIMENSION(klon_glo):: siv_l_glo, siv_lp1_glo, diff_siv_glo
  REAL, DIMENSION(klon_glo,nslay):: bils_glo
  CHARACTER (len = 20)     :: modname = 'limit_slab'
  CHARACTER*2 str2
  LOGICAL                  :: read_bils,read_sst,read_siv

! End declaration
!****************************************************************************************

  ! calculate number of time steps for one day
  lmt_pas = NINT(86400./dtime)
  
! Initialize saved variables
     IF (.NOT. ALLOCATED(bils_save)) THEN
        ALLOCATE(bils_save(klon,nslay), diff_sst_save(klon), diff_siv_save(klon), stat=ierr)
        IF (ierr /= 0) CALL abort_physic('limit_slab', 'pb in allocation',1)
     END IF

  ! F. Codron 5/14: add defaults for bils, diff_sst (0)
  IF (MOD(itime-1, lmt_pas) == 0) THEN   ! time to read
!$OMP MASTER  ! Only master thread
     IF (is_mpi_root) THEN ! Only master processus
        print*,'in limit_slab time to read, itime=',itime
        read_bils=.TRUE.
        read_sst=.TRUE.
        read_siv=.TRUE.
        
        ierr = NF90_OPEN ('limit_slab.nc', NF90_NOWRITE, nid)
        IF (ierr /= NF90_NOERR) THEN
            PRINT *,'LIMIT_SLAB file not found'
            read_bils=.FALSE.
            read_sst=.FALSE.
            read_siv=.FALSE.
        ELSE ! read file
        
        ! La tranche de donnees a lire:
        start(1) = 1
        start(2) = jour
        epais(1) = klon_glo
        epais(2) = 1

!****************************************************************************************
! 2) Read bils and SST/ ice volume tendency
!
!****************************************************************************************
!
! Read bils_glo
        bils_glo(:,:)=0.
        ! First read first layer
        ! try first "BILS_OCE01"
        ierr = NF90_INQ_VARID(nid, 'BILS_OCE01', nvarid)
        IF (ierr /= NF90_NOERR) THEN
            ! Else BILS_OCE
            ierr = NF90_INQ_VARID(nid, 'BILS_OCE', nvarid)
            IF (ierr /= NF90_NOERR) THEN
              read_bils=.FALSE.
            ELSE
              ierr = NF90_GET_VAR(nid,nvarid,bils_glo(:,1),start,epais)
              IF (ierr /= NF90_NOERR) read_bils=.FALSE.
            ENDIF
        ELSE
            ierr = NF90_GET_VAR(nid,nvarid,bils_glo(:,1),start,epais)
            IF (ierr /= NF90_NOERR) read_bils=.FALSE.
        END IF
        ! Try next layers if more than 1
        IF ((nslay.GT.1).AND.read_bils) THEN
          DO i=2,nslay
            WRITE(str2,'(i2.2)') i
            ierr = NF90_INQ_VARID(nid,'BILS_OCE'//str2, nvarid)
            IF (ierr.EQ.NF90_NOERR) THEN
              ierr = NF90_GET_VAR(nid,nvarid,bils_glo(:,i),start,epais)
            ENDIF
            IF (ierr /= NF90_NOERR) THEN
              print *,'WARNING : BILS_OCE not found for layer 2'
            ENDIF
          ENDDO
        ENDIF

! Read sst_glo for this day
        ierr = NF90_INQ_VARID(nid, 'SST', nvarid)
        IF (ierr /= NF90_NOERR)  THEN
            read_sst=.FALSE.
        ELSE
            ierr = NF90_GET_VAR(nid,nvarid,sst_l_glo,start,epais)
            IF (ierr /= NF90_NOERR) read_sst=.FALSE.
! Read sst_glo for one day ahead
            start(2) = jour + 1
            IF (start(2) > 360) start(2)=1
            ierr = NF90_GET_VAR(nid,nvarid,sst_lp1_glo,start,epais)
            IF (ierr /= NF90_NOERR) read_sst=.FALSE.
        END IF

! Read siv_glo for this day
        ierr = NF90_INQ_VARID(nid, 'SICV', nvarid)
        IF (ierr /= NF90_NOERR)  THEN
            read_siv=.FALSE.
        ELSE
            start(2) = jour
            ierr = NF90_GET_VAR(nid,nvarid,siv_l_glo,start,epais)
            IF (ierr /= NF90_NOERR) read_siv=.FALSE.
! Read siv_glo for one day ahead
            start(2) = jour + 1
            IF (start(2) > 360) start(2)=1
            ierr = NF90_GET_VAR(nid,nvarid,siv_lp1_glo,start,epais)
            IF (ierr /= NF90_NOERR) read_siv=.FALSE.
        END IF

!****************************************************************************************
! 5) Close file and distribute variables to all processus
!
!****************************************************************************************
        ierr = NF90_CLOSE(nid)
        IF (ierr /= NF90_NOERR) CALL abort_physic(modname,'Pb when closing file', 1)
        END IF ! Read File 
        IF (read_sst) THEN
! Calculate difference in temperature between this day and one ahead
            DO i=1, klon_glo
               diff_sst_glo(i) = sst_lp1_glo(i) - sst_l_glo(i)
            END DO
        END IF !read_sst
        IF (read_siv) THEN
! Calculate difference in temperature between this day and one ahead
            DO i=1, klon_glo
               diff_siv_glo(i) = siv_lp1_glo(i) - siv_l_glo(i)
            END DO
        END IF !read_siv
     ENDIF ! is_mpi_root

!$OMP END MASTER
!$OMP BARRIER
       
! Send fields to all processes
! Give default values if needed
     CALL bcast(read_bils)
     CALL bcast(read_sst)
     CALL bcast(read_siv)
     PRINT *,'limit_slab sst',read_sst,'siv',read_siv,'qflux',read_bils
     IF (read_bils) THEN
         CALL Scatter(bils_glo, bils_save)
     ELSE
         bils_save(:,:)=0.
     END IF
     IF (read_sst) THEN
         CALL Scatter(diff_sst_glo, diff_sst_save)
     ELSE
         diff_sst_save(:)=0.
     END IF
     IF (read_siv) THEN
         CALL Scatter(diff_siv_glo, diff_siv_save)
     ELSE
         diff_siv_save(:)=0.
     END IF
     
  ENDIF ! time to read

  lmt_bils(:,:) = bils_save(:,:)
  diff_sst(:) = diff_sst_save(:)
  diff_siv(:) = diff_siv_save(:)

END SUBROUTINE limit_slab
