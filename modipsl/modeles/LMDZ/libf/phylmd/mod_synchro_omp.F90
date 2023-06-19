!
! $Id: mod_synchro_omp.F90 3531 2019-06-06 15:08:45Z fairhead $
!
MODULE mod_synchro_omp

    INTEGER,SAVE :: exit_omp  

CONTAINS

  SUBROUTINE Init_synchro_omp
  USE mod_phys_lmdz_para 
  IMPLICIT NONE

  CALL omp_barrier
!$OMP MASTER
  exit_omp=0
!$OMP END MASTER    
  CALL omp_barrier

  END SUBROUTINE Init_Synchro_omp
  
  SUBROUTINE Synchro_omp  
  USE mod_phys_lmdz_para
  IMPLICIT NONE
  LOGICAL :: out
  CHARACTER (LEN=20) :: modname='synchro_omp'
  CHARACTER (LEN=80) :: abort_message
  
    out=.FALSE.
!$OMP BARRIER
!$OMP BARRIER
!$OMP ATOMIC  
    exit_omp=exit_omp+1
!$OMP BARRIER
!$OMP BARRIER
    IF (exit_omp==omp_size) THEN
      out=.TRUE.
    ENDIF
    
    DO WHILE (.NOT. out)
!$OMP BARRIER
      IF (exit_omp==omp_size) out=.TRUE.
!$OMP BARRIER
    ENDDO

!$OMP BARRIER
!$OMP MASTER
    exit_omp=0
!$OMP END MASTER
!$OMP BARRIER

    IF (exit_omp/=0) THEN
       abort_message='synchro_omp'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

  END SUBROUTINE Synchro_omp


END MODULE mod_synchro_omp
