MODULE etat0_limit_unstruct_mod

  LOGICAL, SAVE  :: create_etat0_limit
!$OMP THREADPRIVATE(create_etat0_limit) 




CONTAINS
  
  SUBROUTINE init_etat0_limit_unstruct
#ifdef CPP_XIOS
  USE xios, ONLY: xios_set_axis_attr, xios_set_fieldgroup_attr, &
                  xios_set_filegroup_attr, xios_set_file_attr
  USE mod_phys_lmdz_para, ONLY: is_omp_master
  USE mod_grid_phy_lmdz, ONLY: grid_type, unstructured
  USE ioipsl, ONLY : ioget_year_len
  USE ioipsl_getin_p_mod, ONLY: getin_p
  USE time_phylmdz_mod, ONLY : annee_ref
  USE create_etat0_unstruct_mod, ONLY: init_create_etat0_unstruct
  IMPLICIT NONE
  
    INTEGER :: iflag_phys,i
    INTEGER :: ndays
    REAL,ALLOCATABLE :: value(:)
    
      IF (grid_type==unstructured) THEN
        CALL getin_p("iflag_phys",iflag_phys)
        
        CALL getin_p('create_etat0_limit',create_etat0_limit)
        
        ndays=ioget_year_len(annee_ref)
        ALLOCATE(value(ndays))
        DO i=1,ndays
          value(i)=i-1
        ENDDO
        
        IF (is_omp_master) CALL xios_set_axis_attr("time_year",n_glo=ndays,value=value) 
        
        IF (create_etat0_limit) THEN
          IF (iflag_phys<100) THEN
            IF (is_omp_master) CALL xios_set_fieldgroup_attr("etat0_limit_read",read_access=.TRUE.,enabled=.TRUE.)
            IF (is_omp_master) CALL xios_set_filegroup_attr("etat0_limit_read",enabled=.TRUE.)
          ENDIF
          IF (is_omp_master) CALL xios_set_file_attr("limit_write",enabled=.TRUE.)
          CALL init_create_etat0_unstruct
        ENDIF
      
      ENDIF  

#endif
  END SUBROUTINE init_etat0_limit_unstruct
  
  SUBROUTINE create_etat0_limit_unstruct
#ifdef CPP_XIOS
  USE mod_grid_phy_lmdz, ONLY: grid_type, unstructured
  USE create_etat0_unstruct_mod, ONLY: create_etat0_unstruct
  USE create_limit_unstruct_mod, ONLY: create_limit_unstruct
  USE phyaqua_mod, ONLY: iniaqua
  USE phys_cal_mod, only: year_len
  USE mod_phys_lmdz_para, ONLY: is_omp_master
  USE ioipsl_getin_p_mod, ONLY: getin_p
  USE dimphy, ONLY: klon
  USE xios, ONLY: xios_context_finalize, xios_set_current_context, &
                  xios_finalize
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE
      INTEGER :: iflag_phys
      INTEGER :: ierr
      CHARACTER (LEN=20) :: modname='create_etat0_limit_unstruct'
      CHARACTER (LEN=80) :: abort_message
     
      IF (grid_type==unstructured) THEN
  
        CALL getin_p("iflag_phys",iflag_phys)

        IF (iflag_phys<100) THEN
          IF ( create_etat0_limit) THEN
              CALL create_etat0_unstruct
              CALL create_limit_unstruct
              IF (is_omp_master)  THEN
                CALL xios_context_finalize()
                CALL xios_set_current_context("icosagcm")   ! very bad, need to find an other solution
                CALL xios_context_finalize()
                CALL xios_finalize()
#ifdef CPP_MPI    
                CALL MPI_Finalize(ierr)
#endif
                abort_message='create_etat0_limit_unstruct, Initial state file are created, all is fine' 
                CALL abort_physic(modname,abort_message,0)
                STOP 0
              ENDIF
!$OMP BARRIER
              abort_message='create_etat0_limit_unstruct, Initial state file are created, all is fine'
              CALL abort_physic(modname,abort_message,0)
          ENDIF 
        ELSE
          IF (create_etat0_limit) THEN 
            CALL iniaqua(klon,year_len,iflag_phys)
              IF (is_omp_master)  THEN
                CALL xios_context_finalize()
                CALL xios_set_current_context("icosagcm")   ! very bad, need to find an other solution
                CALL xios_context_finalize()
                CALL xios_finalize()
#ifdef CPP_MPI
                CALL MPI_Finalize(ierr)
#endif
              ENDIF
!$OMP BARRIER
              abort_message='create_etat0_limit_unstruct, Initial state file are created, all is fine' 
              CALL abort_physic(modname,abort_message,0)
              STOP 0
          ENDIF
        ENDIF
      ENDIF
        
#endif
  END SUBROUTINE create_etat0_limit_unstruct
  
END MODULE etat0_limit_unstruct_mod

