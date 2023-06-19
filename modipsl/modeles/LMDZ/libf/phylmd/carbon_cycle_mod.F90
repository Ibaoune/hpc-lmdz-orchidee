MODULE carbon_cycle_mod
!=======================================================================
!   Authors: Patricia Cadule and Laurent Fairhead
!            base sur un travail anterieur mene par Patricia Cadule et Josefine Ghattas
!
!  Purpose and description:
!  -----------------------
! Control module for the carbon CO2 tracers :
!   - Initialisation of carbon cycle fields
!   - Definition of fluxes to be exchanged
!
! Rest of code is in tracco2i.F90
!
! Le cas online/offline est defini par le flag carbon_cycle_cpl (y/n)
! Le transport du traceur CO2 est defini par le flag carbon_cycle_tr (y/n)
! la provenance des champs (termes de puits) est denini par le flag level_coupling_esm
!
! level_coupling_esm : level of coupling of the biogeochemical fields between
! LMDZ, ORCHIDEE and NEMO
! Definitions of level_coupling_esm in physiq.def
! level_coupling_esm = 0  ! No field exchange between LMDZ and ORCHIDEE models
!                         ! No field exchange between LMDZ and NEMO
! level_coupling_esm = 1  ! Field exchange between LMDZ and ORCHIDEE models
!                         ! No field exchange between LMDZ and NEMO models
! level_coupling_esm = 2  ! No field exchange between LMDZ and ORCHIDEE models
!                         ! Field exchange between LMDZ and NEMO models
! level_coupling_esm = 3  ! Field exchange between LMDZ and ORCHIDEE models
!                         ! Field exchange between LMDZ and NEMO models
!=======================================================================

  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: carbon_cycle_init, infocfields_init

! Variables read from parmeter file physiq.def
  LOGICAL, PUBLIC :: carbon_cycle_cpl       ! Coupling of CO2 fluxes between LMDZ/ORCHIDEE and LMDZ/OCEAN(PISCES) 
!$OMP THREADPRIVATE(carbon_cycle_cpl)
  LOGICAL, PUBLIC :: carbon_cycle_tr        ! 3D transport of CO2 in the atmosphere, parameter read in conf_phys
!$OMP THREADPRIVATE(carbon_cycle_tr)
  LOGICAL, PUBLIC :: carbon_cycle_rad       ! CO2 interactive radiatively 
!$OMP THREADPRIVATE(carbon_cycle_rad)
  INTEGER, PUBLIC :: level_coupling_esm ! Level of coupling for the ESM - 0, 1, 2, 3
!$OMP THREADPRIVATE(level_coupling_esm)
  REAL, PUBLIC :: RCO2_glo
!$OMP THREADPRIVATE(RCO2_glo)
  REAL, PUBLIC :: RCO2_tot
!$OMP THREADPRIVATE(RCO2_tot)

  LOGICAL :: carbon_cycle_emis_comp_omp=.FALSE.
  LOGICAL :: carbon_cycle_emis_comp=.FALSE. ! Calculation of emission compatible
!$OMP THREADPRIVATE(carbon_cycle_emis_comp)

  LOGICAL :: RCO2_inter_omp
  LOGICAL :: RCO2_inter  ! RCO2 interactive : if true calculate new value RCO2 for the radiation scheme
!$OMP THREADPRIVATE(RCO2_inter)

! Scalare values when no transport, from physiq.def
  REAL :: fos_fuel_s_omp
  REAL :: fos_fuel_s  ! carbon_cycle_fos_fuel dans physiq.def
!$OMP THREADPRIVATE(fos_fuel_s)
  REAL :: emis_land_s ! not yet implemented
!$OMP THREADPRIVATE(emis_land_s)

  REAL :: airetot     ! Total area of the earth surface
!$OMP THREADPRIVATE(airetot)

  INTEGER :: ntr_co2  ! Number of tracers concerning the carbon cycle
!$OMP THREADPRIVATE(ntr_co2)

! fco2_ocn_day : flux CO2 from ocean for 1 day (cumulated) [gC/m2/d]. Allocation and initalization done in cpl_mod
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_ocn_day 
!$OMP THREADPRIVATE(fco2_ocn_day)

  REAL, DIMENSION(:), ALLOCATABLE :: fco2_land_day   ! flux CO2 from land for 1 day (cumulated)  [gC/m2/d]
!$OMP THREADPRIVATE(fco2_land_day)
  REAL, DIMENSION(:), ALLOCATABLE :: fco2_lu_day     ! Emission from land use change for 1 day (cumulated) [gC/m2/d]
!$OMP THREADPRIVATE(fco2_lu_day)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_ff ! Emission from fossil fuel [kgCO2/m2/s]
!$OMP THREADPRIVATE(fco2_ff)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_bb ! Emission from biomass burning [kgCO2/m2/s]
!$OMP THREADPRIVATE(fco2_bb)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_land  ! Net flux from terrestrial ecocsystems [kgCO2/m2/s]
!$OMP THREADPRIVATE(fco2_land)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_land_nbp  ! Net flux from terrestrial ecocsystems [kgCO2/m2/s]
!$OMP THREADPRIVATE(fco2_land_nbp)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_land_nep  ! Net flux from terrestrial ecocsystems [kgCO2/m2/s]
!$OMP THREADPRIVATE(fco2_land_nep)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_land_fLuc  ! Net flux from terrestrial ecocsystems [kgCO2/m2/s]
!$OMP THREADPRIVATE(fco2_land_fLuc)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_land_fwoodharvest  ! Net flux from terrestrial ecocsystems [kgCO2/m2/s]
!$OMP THREADPRIVATE(fco2_land_fwoodharvest)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_land_fHarvest  ! Net flux from terrestrial ecocsystems [kgCO2/m2/s]
!$OMP THREADPRIVATE(fco2_land_fHarvest)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_ocean ! Net flux from ocean [kgCO2/m2/s]
!$OMP THREADPRIVATE(fco2_ocean)

  REAL, DIMENSION(:,:), ALLOCATABLE :: dtr_add       ! Tracer concentration to be injected 
!$OMP THREADPRIVATE(dtr_add)

! Following 2 fields will be allocated and initialized in surf_land_orchidee
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_land_inst  ! flux CO2 from land at one time step
!$OMP THREADPRIVATE(fco2_land_inst)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_lu_inst    ! Emission from land use change at one time step
!$OMP THREADPRIVATE(fco2_lu_inst)

! Calculated co2 field to be send to the ocean via the coupler and to ORCHIDEE 
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: co2_send ! Field allocated in phyetat0
!$OMP THREADPRIVATE(co2_send)

  INTEGER, PARAMETER, PUBLIC :: id_CO2=1              !--temporaire OB -- to be changed

! nbfields : total number of fields 
  INTEGER, PUBLIC :: nbcf
!$OMP THREADPRIVATE(nbcf)

! nbcf_in : number of fields IN
  INTEGER, PUBLIC  :: nbcf_in
!$OMP THREADPRIVATE(nbcf_in)

! nbcf_in_orc : number of fields IN 
  INTEGER, PUBLIC  :: nbcf_in_orc
!$OMP THREADPRIVATE(nbcf_in_orc)

! nbcf_in_inca : number of fields IN (from INCA)
  INTEGER, PUBLIC  :: nbcf_in_inca
!$OMP THREADPRIVATE(nbcf_in_inca)

! nbcf_in_nemo : number of fields IN (from nemo)
  INTEGER, PUBLIC  :: nbcf_in_nemo
!$OMP THREADPRIVATE(nbcf_in_nemo)

! nbcf_in_ant : number of fields IN (from anthropogenic sources)
  INTEGER, PUBLIC  :: nbcf_in_ant
!$OMP THREADPRIVATE(nbcf_in_ant)

! nbcf_out : number of fields OUT
  INTEGER, PUBLIC :: nbcf_out
!$OMP THREADPRIVATE(nbcf_out)

! Name of variables
  CHARACTER(len=25), ALLOCATABLE, DIMENSION(:), PUBLIC :: cfname     ! coupling field short name for restart (?) and diagnostics
!$OMP THREADPRIVATE(cfname)

  CHARACTER(len=25), ALLOCATABLE, DIMENSION(:), PUBLIC :: cfname_in  ! coupling field short name for restart (?) and diagnostics
!$OMP THREADPRIVATE(cfname_in)

  CHARACTER(len=25), ALLOCATABLE, DIMENSION(:), PUBLIC :: cfname_out ! coupling field short name for restart (?) and diagnostics
!$OMP THREADPRIVATE(cfname_out)

  CHARACTER(len=15), ALLOCATABLE, DIMENSION(:), PUBLIC :: cfunits_in  !  coupling field units for diagnostics
!$OMP THREADPRIVATE(cfunits_in)

  CHARACTER(len=15), ALLOCATABLE, DIMENSION(:), PUBLIC :: cfunits_out !  coupling field units for diagnostics
!$OMP THREADPRIVATE(cfunits_out)

  CHARACTER(len=120), ALLOCATABLE, DIMENSION(:), PUBLIC :: cftext_in  ! coupling field long name for diagnostics
!$OMP THREADPRIVATE(cftext_in)

  CHARACTER(len=120), ALLOCATABLE, DIMENSION(:), PUBLIC :: cftext_out ! coupling field long name for diagnostics
!$OMP THREADPRIVATE(cftext_out)

  CHARACTER(len=5), ALLOCATABLE, DIMENSION(:), PUBLIC :: cfmod1 ! model 1 (rreference) : LMDz
!$OMP THREADPRIVATE(cfmod1)

  CHARACTER(len=5), ALLOCATABLE, DIMENSION(:), PUBLIC :: cfmod2 ! model 2
!$OMP THREADPRIVATE(cfmod2)

  CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:), PUBLIC :: field_out_names
!$OMP THREADPRIVATE(field_out_names)

  CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:), PUBLIC :: field_in_names
!$OMP THREADPRIVATE(field_in_names)

  REAL, ALLOCATABLE, DIMENSION(:,:), PUBLIC :: fields_in   !  klon,nbcf_in
!$OMP THREADPRIVATE(fields_in)

  REAL, ALLOCATABLE, DIMENSION(:,:), PUBLIC :: yfields_in  !  knon,nbcf_in
!$OMP THREADPRIVATE(yfields_in)

  REAL, ALLOCATABLE, DIMENSION(:,:), PUBLIC :: fields_out  !  klon,nbcf_out
!$OMP THREADPRIVATE(fields_out)

  REAL, ALLOCATABLE, DIMENSION(:,:), PUBLIC :: yfields_out !  knon,nbcf_out
!$OMP THREADPRIVATE(yfields_out)

  TYPE, PUBLIC :: co2_trac_type
     CHARACTER(len = 8) :: name       ! Tracer name in tracer.def
     INTEGER            :: id         ! Index in total tracer list, tr_seri
     CHARACTER(len=30)  :: file       ! File name
     LOGICAL            :: cpl        ! True if this tracers is coupled from ORCHIDEE or PISCES. 
                                      ! False if read from file.
     INTEGER            :: updatefreq ! Frequence to inject in second
     INTEGER            :: readstep   ! Actual time step to read in file
     LOGICAL            :: updatenow  ! True if this tracer should be updated this time step
  END TYPE co2_trac_type
  INTEGER,PARAMETER :: maxco2trac=5  ! Maximum number of different CO2 fluxes
  TYPE(co2_trac_type), DIMENSION(maxco2trac) :: co2trac

CONTAINS
  
  SUBROUTINE carbon_cycle_init()
    ! This subroutine is called from tracco2i_init, which is called from phytrac_init only at first timestep.
    ! - Allocate variables. These variables must be allocated before first call to phys_output_write in physiq.

    USE dimphy
    USE IOIPSL
    USE print_control_mod, ONLY: lunout

    IMPLICIT NONE
    INCLUDE "clesphys.h"
 
! Local variables
    INTEGER               :: ierr

    IF (carbon_cycle_cpl) THEN

       ierr=0

       IF (.NOT.ALLOCATED(fco2_land)) ALLOCATE(fco2_land(klon), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('carbon_cycle_init', 'pb in allocation fco2_land',1)
       fco2_land(1:klon) = 0.

       IF (.NOT.ALLOCATED(fco2_land_nbp)) ALLOCATE(fco2_land_nbp(klon), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('carbon_cycle_init', 'pb in allocation fco2_land_nbp',1)
       fco2_land_nbp(1:klon) = 0.

       IF (.NOT.ALLOCATED(fco2_land_nep)) ALLOCATE(fco2_land_nep(klon), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('carbon_cycle_init', 'pb in allocation fco2_land_nep',1)
       fco2_land_nep(1:klon) = 0.

       IF (.NOT.ALLOCATED(fco2_land_fLuc)) ALLOCATE(fco2_land_fLuc(klon), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('carbon_cycle_init', 'pb in allocation fco2_land_fLuc',1)
       fco2_land_fLuc(1:klon) = 0.

       IF (.NOT.ALLOCATED(fco2_land_fwoodharvest)) ALLOCATE(fco2_land_fwoodharvest(klon), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('carbon_cycle_init', 'pb in allocation fco2_land_fwoodharvest',1)
       fco2_land_fwoodharvest(1:klon) = 0.

       IF (.NOT.ALLOCATED(fco2_land_fHarvest)) ALLOCATE(fco2_land_fHarvest(klon), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('carbon_cycle_init', 'pb in allocation fco2_land_fHarvest',1)
       fco2_land_fHarvest(1:klon) = 0.

       IF (.NOT.ALLOCATED(fco2_ff)) ALLOCATE(fco2_ff(klon), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('carbon_cycle_init', 'pb in allocation fco2_ff',1)
       fco2_ff(1:klon) = 0.

       IF (.NOT.ALLOCATED(fco2_bb)) ALLOCATE(fco2_bb(klon), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('carbon_cycle_init', 'pb in allocation fco2_bb',1)
       fco2_bb(1:klon) = 0.

       IF (.NOT.ALLOCATED(fco2_ocean)) ALLOCATE(fco2_ocean(klon), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('carbon_cycle_init', 'pb in allocation fco2_ocean',1)
       fco2_bb(1:klon) = 0.
    ENDIF

  END SUBROUTINE carbon_cycle_init

  SUBROUTINE infocfields_init

!    USE control_mod, ONLY: planet_type
    USE phys_cal_mod, ONLY : mth_cur
    USE mod_synchro_omp
    USE mod_phys_lmdz_para, ONLY: is_mpi_root, is_omp_root
    USE mod_phys_lmdz_transfert_para
    USE mod_phys_lmdz_omp_transfert
    USE dimphy, ONLY: klon

    IMPLICIT NONE

!=======================================================================
!
!   Authors: Patricia Cadule and Laurent Fairhead  
!   -------
!
!  Purpose and description:
!  -----------------------
!
! Infofields
! this routine enables to define the field exchanges in both directions between
! the atmospheric circulation model (LMDZ) and ORCHIDEE. In the future this
! routing might apply to other models (e.g., NEMO, INCA, ...).
! Therefore, currently with this routine, it is possible to define the coupling
! fields only between LMDZ and ORCHIDEE.
! The coupling_fields.def file enables to define the name of the exchanged
! fields at the coupling interface.
! field_in_names : the set of names of the exchanged fields in input to ORCHIDEE
! (LMDZ to ORCHIDEE)
! field_out_names : the set of names of the exchanged fields in output of
! ORCHIDEE (ORCHIDEE to LMDZ)
! n : the number of exchanged fields at th coupling interface
! nb_fields_in : number of inputs fields to ORCHIDEE (LMDZ to ORCHIDEE)
! nb_fields_out : number of ouput fields of ORCHIDEE (ORCHIDEE to LMDZ)
!
! The syntax for coupling_fields.def is as follows:
! IMPORTANT: each column entry must be separated from the previous one by 3
! spaces and only that
! field name         coupling          model 1         model 2         long_name
!                    direction
!   10char  -3spaces-  3char  -3spaces- 4char -3spaces- 4char -3spaces-  30char 
!
! n
! FIELD1 IN LMDZ ORC
! ....
! FIELD(j) IN LMDZ ORC
! FIELD(j+1) OUT LMDZ ORC
! ...
! FIELDn OUT LMDZ ORC   
! 
!=======================================================================
!   ... 22/12/2017 ....
!-----------------------------------------------------------------------
! Declarations

  INCLUDE "clesphys.h"
  INCLUDE "dimensions.h"
  INCLUDE "iniprint.h"

! Local variables

  INTEGER :: iq,  ierr, stat, error

  CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:), SAVE  :: cfname_root
  CHARACTER(LEN=120), ALLOCATABLE, DIMENSION(:), SAVE :: cftext_root
  CHARACTER(LEN=15), ALLOCATABLE, DIMENSION(:), SAVE  :: cfunits_root

  CHARACTER(len=3), ALLOCATABLE, DIMENSION(:) :: cfintent_root
  CHARACTER(len=5), ALLOCATABLE, DIMENSION(:) :: cfmod1_root
  CHARACTER(len=5), ALLOCATABLE, DIMENSION(:) :: cfmod2_root

  LOGICAL, ALLOCATABLE, DIMENSION(:), SAVE :: mask_in_root
  LOGICAL, ALLOCATABLE, DIMENSION(:), SAVE :: mask_out_root

  CHARACTER(len=*),parameter :: modname="infocfields"

  CHARACTER(len=10),SAVE :: planet_type="earth"

!-----------------------------------------------------------------------

nbcf=0
nbcf_in=0
nbcf_out=0

IF (planet_type=='earth') THEN

    IF (is_mpi_root .AND. is_omp_root) THEN

       IF (level_coupling_esm.GT.0) THEN

          OPEN(200,file='coupling_fields.def',form='formatted',status='old', iostat=ierr)

          IF (ierr.EQ.0) THEN

             WRITE(lunout,*) trim(modname),': Open coupling_fields.def : ok'
             READ(200,*) nbcf
             WRITE(lunout,*) 'infocfields_mod.F90 --- nbcf=',nbcf
             ALLOCATE(cfname_root(nbcf))
             ALLOCATE(cfintent_root(nbcf))
             ALLOCATE(cfmod1_root(nbcf))
             ALLOCATE(cfmod2_root(nbcf))
             ALLOCATE(cftext_root(nbcf))
             ALLOCATE(cfunits_root(nbcf))
             ALLOCATE(mask_in_root(nbcf))
             ALLOCATE(mask_out_root(nbcf))

             nbcf_in=0
             nbcf_out=0

             DO iq=1,nbcf
                WRITE(lunout,*) 'infofields : field=',iq
                READ(200,'(A15,3X,A3,3X,A5,3X,A5,3X,A120,3X,A15)',IOSTAT=ierr) &
                   cfname_root(iq),cfintent_root(iq),cfmod1_root(iq),cfmod2_root(iq),cftext_root(iq),cfunits_root(iq)
                cfname_root(iq)=TRIM(cfname_root(iq))
                cfintent_root(iq)=TRIM(cfintent_root(iq))
                cfmod1_root(iq)=TRIM(cfmod1_root(iq))
                cfmod2_root(iq)=TRIM(cfmod2_root(iq))
                cftext_root(iq)=TRIM(cftext_root(iq))
                cfunits_root(iq)=TRIM(cfunits_root(iq))
                WRITE(lunout,*) 'coupling field: ',cfname_root(iq), &  
                               ', number: ',iq,', INTENT: ',cfintent_root(iq)
                WRITE(lunout,*) 'coupling field: ',cfname_root(iq), & 
                               ', number: ',iq,', model 1 (ref): ',cfmod1_root(iq),', model 2: ',cfmod2_root(iq)
                WRITE(lunout,*) 'coupling field: ',cfname_root(iq), & 
                               ', number: ',iq,', long name: ',cftext_root(iq),', units ',cfunits_root(iq)
                IF (nbcf_in+nbcf_out.LT.nbcf) THEN
                  IF (cfintent_root(iq).NE.'OUT') THEN
                    nbcf_in=nbcf_in+1
                    mask_in_root(iq)=.TRUE.
                    mask_out_root(iq)=.FALSE.
                  ELSE IF (cfintent_root(iq).EQ.'OUT') THEN
                    nbcf_out=nbcf_out+1
                    mask_in_root(iq)=.FALSE.
                    mask_out_root(iq)=.TRUE.
                  ENDIF
                ELSE
                  WRITE(lunout,*) 'abort_gcm --- nbcf    : ',nbcf
                  WRITE(lunout,*) 'abort_gcm --- nbcf_in : ',nbcf_in
                  WRITE(lunout,*) 'abort_gcm --- nbcf_out: ',nbcf_out
                  CALL abort_physic('infocfields_init','Problem in the definition of the coupling fields',1)
               ENDIF
             ENDDO !DO iq=1,nbcf
          ELSE
             WRITE(lunout,*) trim(modname),': infocfields_mod.F90 --- Problem in opening coupling_fields.def'
             WRITE(lunout,*) trim(modname),': infocfields_mod.F90 --- WARNING using defaut values'
          ENDIF ! ierr
          CLOSE(200)

       ENDIF ! level_coupling_esm

    ENDIF !   (is_mpi_root .AND. is_omp_root)
!$OMP BARRIER

    CALL bcast(nbcf)
    CALL bcast(nbcf_in)
    CALL bcast(nbcf_out)

    WRITE(lunout,*) 'infocfields_mod.F90 --- nbcf    =',nbcf
    WRITE(lunout,*) 'infocfields_mod.F90 --- nbcf_in =',nbcf_in
    WRITE(lunout,*) 'infocfields_mod.F90 --- nbcf_out=',nbcf_out

    ALLOCATE(cfname(nbcf))
    ALLOCATE(cfname_in(nbcf_in))
    ALLOCATE(cftext_in(nbcf_in))
    ALLOCATE(cfname_out(nbcf_out))
    ALLOCATE(cftext_out(nbcf_out))
    ALLOCATE(cfmod1(nbcf))
    ALLOCATE(cfmod2(nbcf))
    ALLOCATE(cfunits_in(nbcf_in))
    ALLOCATE(cfunits_out(nbcf_out))
       
    IF (is_mpi_root .AND. is_omp_root) THEN

        IF (nbcf.GT.0)     cfname=cfname_root
        IF (nbcf_in.GT.0)  cfname_in=PACK(cfname_root,mask_in_root)
        IF (nbcf_out.GT.0) cfname_out=PACK(cfname_root,mask_out_root)
        IF (nbcf_in.GT.0)  cftext_in=PACK(cftext_root,mask_in_root)
        IF (nbcf_out.GT.0) cftext_out=PACK(cftext_root,mask_out_root)
        IF (nbcf.GT.0)     cfmod1=cfmod1_root
        IF (nbcf.GT.0)     cfmod2=cfmod2_root
        IF (nbcf_in.GT.0)  cfunits_in=PACK(cfunits_root,mask_in_root)
        IF (nbcf_out.GT.0) cfunits_out=PACK(cfunits_root,mask_out_root)

        nbcf_in_orc=0
        nbcf_in_nemo=0
        nbcf_in_inca=0
        nbcf_in_ant=0

        DO iq=1,nbcf
            IF (cfmod1(iq) == "ORC")  nbcf_in_orc  = nbcf_in_orc  + 1  
            IF (cfmod1(iq) == "NEMO") nbcf_in_nemo = nbcf_in_nemo + 1  
            IF (cfmod1(iq) == "INCA") nbcf_in_inca = nbcf_in_inca + 1
            IF (cfmod1(iq) == "ALL")  nbcf_in_orc  = nbcf_in_orc  + 1  ! ALL = ORC/NEMO/INCA
            IF (cfmod1(iq) == "ALL")  nbcf_in_nemo = nbcf_in_nemo + 1  ! ALL = ORC/NEMO/INCA 
            IF (cfmod1(iq) == "ALL")  nbcf_in_inca = nbcf_in_inca + 1  ! ALL = ORC/NEMO/INCA 
            IF (cfmod1(iq) == "ANT")  nbcf_in_ant  = nbcf_in_ant  + 1  
        ENDDO

    ENDIF !   (is_mpi_root .AND. is_omp_root)
!$OMP BARRIER

    CALL bcast(nbcf_in_orc)
    CALL bcast(nbcf_in_nemo)
    CALL bcast(nbcf_in_inca)
    CALL bcast(nbcf_in_ant)

    WRITE(lunout,*) 'nbcf_in_orc  =',nbcf_in_orc
    WRITE(lunout,*) 'nbcf_in_nemo =',nbcf_in_nemo
    WRITE(lunout,*) 'nbcf_in_inca =',nbcf_in_inca
    WRITE(lunout,*) 'nbcf_in_ant  =',nbcf_in_ant

    IF (nbcf_in.GT.0) THEN 
        DO iq=1,nbcf_in
          CALL bcast(cfname_in(iq))
          CALL bcast(cftext_in(iq))
          CALL bcast(cfunits_in(iq))
        ENDDO
    ENDIF

    IF (nbcf_out.GT.0) THEN
        DO iq=1,nbcf_out
          CALL bcast(cfname_out(iq))
          CALL bcast(cftext_out(iq))
          CALL bcast(cfunits_out(iq))
        ENDDO
    ENDIF

    IF (nbcf.GT.0) THEN
        DO iq=1,nbcf
          CALL bcast(cfmod1(iq))
          CALL bcast(cfmod2(iq))
        ENDDO
    ENDIF

    IF (nbcf_in.GT.0)  WRITE(lunout,*)'infocfields_mod --- cfname_in: ',cfname_in
    IF (nbcf_out.GT.0) WRITE(lunout,*)'infocfields_mod --- cfname_out: ',cfname_out

    IF (nbcf_in.GT.0)  WRITE(lunout,*)'infocfields_mod --- cftext_in: ',cftext_in
    IF (nbcf_out.GT.0) WRITE(lunout,*)'infocfields_mod --- cftext_out: ',cftext_out

    IF (nbcf.GT.0) WRITE(lunout,*)'infocfields_mod --- cfmod1: ',cfmod1 
    IF (nbcf.GT.0) WRITE(lunout,*)'infocfields_mod --- cfmod2: ',cfmod2 

    IF (nbcf_in.GT.0)  WRITE(lunout,*)'infocfunits_mod --- cfunits_in: ',cfunits_in
    IF (nbcf_out.GT.0) WRITE(lunout,*)'infocfunits_mod --- cfunits_out: ',cfunits_out

    IF (nbcf_in.GT.0)  WRITE(*,*)'infocfields_init --- number of fields in to LMDZ: ',nbcf_in
    IF (nbcf_out.GT.0) WRITE(*,*)'infocfields_init --- number of fields out of LMDZ: ',nbcf_out

 ELSE
 ! Default values for other planets
    nbcf=0 
    nbcf_in=0
    nbcf_out=0
 ENDIF ! planet_type

 ALLOCATE(fields_in(klon,nbcf_in),stat=error)
 IF (error /= 0)  CALL abort_physic(modname,'Pb in allocation fields_in',1)
 ALLOCATE(yfields_in(klon,nbcf_in),stat=error)
 IF (error /= 0)  CALL abort_physic(modname,'Pb in allocation yfields_in',1)
 ALLOCATE(fields_out(klon,nbcf_out),stat=error)
 IF (error /= 0)  CALL abort_physic(modname,'Pb in allocation fields_out',1)
 ALLOCATE(yfields_out(klon,nbcf_out),stat=error)
 IF (error /= 0)  CALL abort_physic(modname,'Pb in allocation yfields_out',1)

END SUBROUTINE infocfields_init

END MODULE carbon_cycle_mod
