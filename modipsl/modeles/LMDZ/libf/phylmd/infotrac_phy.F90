
! $Id: $

MODULE infotrac_phy

! Infotrac for physics; for now contains the same information as infotrac for
! the dynamics (could be further cleaned) and is initialized using values
! provided by the dynamics

! nqtot : total number of tracers and higher order of moment, water vapor and liquid included
  INTEGER, SAVE :: nqtot
!$OMP THREADPRIVATE(nqtot)

!CR: on ajoute le nombre de traceurs de l eau
  INTEGER, SAVE :: nqo
!$OMP THREADPRIVATE(nqo)

! nbtr : number of tracers not including higher order of moment or water vapor or liquid
!        number of tracers used in the physics
  INTEGER, SAVE :: nbtr
!$OMP THREADPRIVATE(nbtr)

#ifdef CPP_StratAer
! nbtr_bin: number of aerosol bins for StratAer model
! nbtr_sulgas: number of sulfur gases for StratAer model
  INTEGER, SAVE :: nbtr_bin, nbtr_sulgas
!$OMP THREADPRIVATE(nbtr_bin,nbtr_sulgas)
  INTEGER, SAVE :: id_OCS_strat, id_SO2_strat, id_H2SO4_strat, id_BIN01_strat
!$OMP THREADPRIVATE(id_OCS_strat, id_SO2_strat, id_H2SO4_strat, id_BIN01_strat)
#endif

! CRisi: nb traceurs pères= directement advectés par l'air
  INTEGER, SAVE :: nqperes
!$OMP THREADPRIVATE(nqperes)

! Name variables
  CHARACTER(len=20), ALLOCATABLE, DIMENSION(:), SAVE :: tname ! tracer short name for restart and diagnostics
  CHARACTER(len=23), ALLOCATABLE, DIMENSION(:), SAVE :: ttext ! tracer long name for diagnostics
!$OMP THREADPRIVATE(tname,ttext)

!! iadv  : index of trasport schema for each tracer
!  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: iadv

! niadv : vector keeping the coorspondance between all tracers(nqtot) treated in the 
!         dynamic part of the code and the tracers (nbtr+2) used in the physics part of the code. 
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: niadv ! equivalent dyn / physique
!$OMP THREADPRIVATE(niadv)

! CRisi: tableaux de fils
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: nqfils
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: nqdesc ! nombres de fils + nombre de tous les petits fils sur toutes les générations
  INTEGER, SAVE :: nqdesc_tot
  INTEGER, ALLOCATABLE, DIMENSION(:,:), SAVE    :: iqfils
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: iqpere
!$OMP THREADPRIVATE(nqfils,nqdesc,nqdesc_tot,iqfils,iqpere)

! conv_flg(it)=0 : convection desactivated for tracer number it 
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: conv_flg
!$OMP THREADPRIVATE(conv_flg)

! pbl_flg(it)=0  : boundary layer diffusion desactivaded for tracer number it 
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: pbl_flg
!$OMP THREADPRIVATE(pbl_flg)

  CHARACTER(len=4),SAVE :: type_trac
!$OMP THREADPRIVATE(type_trac)
  CHARACTER(len=8),DIMENSION(:),ALLOCATABLE, SAVE :: solsym
!$OMP THREADPRIVATE(solsym)
   
    ! CRisi: cas particulier des isotopes
    LOGICAL,SAVE :: ok_isotopes,ok_iso_verif,ok_isotrac,ok_init_iso
!$OMP THREADPRIVATE(ok_isotopes,ok_iso_verif,ok_isotrac,ok_init_iso)
    INTEGER :: niso_possibles   
    PARAMETER ( niso_possibles=5)
    real, DIMENSION (niso_possibles),SAVE :: tnat,alpha_ideal
!$OMP THREADPRIVATE(tnat,alpha_ideal)
    LOGICAL, DIMENSION(niso_possibles),SAVE ::  use_iso
!$OMP THREADPRIVATE(use_iso)
    INTEGER, ALLOCATABLE, DIMENSION(:,:), SAVE ::  iqiso ! donne indice iq en fn de (ixt,phase) 
!$OMP THREADPRIVATE(iqiso)
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE ::  iso_num ! donne numéro iso entre 1 et niso_possibles en fn de nqtot
!$OMP THREADPRIVATE(iso_num)
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE ::  iso_indnum ! donne numéro iso entre 1 et niso effectif en fn de nqtot
!$OMP THREADPRIVATE(iso_indnum)
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE ::  zone_num ! donne numéro de la zone de tracage en fn de nqtot
!$OMP THREADPRIVATE(zone_num)
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE ::  phase_num ! donne numéro de la zone de tracage en fn de nqtot
!$OMP THREADPRIVATE(phase_num)
    INTEGER, DIMENSION(niso_possibles), SAVE :: indnum_fn_num ! donne indice entre entre 1 et niso en fonction du numéro d isotope entre 1 et niso_possibles
!$OMP THREADPRIVATE(indnum_fn_num)
    INTEGER, ALLOCATABLE, DIMENSION(:,:), SAVE ::  index_trac ! numéro ixt en fn izone, indnum entre 1 et niso
!$OMP THREADPRIVATE(index_trac)
    INTEGER,SAVE :: niso,ntraceurs_zone,ntraciso
!$OMP THREADPRIVATE(niso,ntraceurs_zone,ntraciso)
 
CONTAINS

  SUBROUTINE init_infotrac_phy(nqtot_,nqo_,nbtr_,tname_,ttext_,type_trac_,&
                               niadv_,conv_flg_,pbl_flg_,solsym_,&
                               nqfils_,nqdesc_,nqdesc_tot_,iqfils_,iqpere_,&
                               ok_isotopes_,ok_iso_verif_,ok_isotrac_,&
                               ok_init_iso_,niso_possibles_,tnat_,&
                               alpha_ideal_,use_iso_,iqiso_,iso_num_,&
                               iso_indnum_,zone_num_,phase_num_,&
                               indnum_fn_num_,index_trac_,&
                               niso_,ntraceurs_zone_,ntraciso_&
#ifdef CPP_StratAer
                               ,nbtr_bin_,nbtr_sulgas_&
                               ,id_OCS_strat_,id_SO2_strat_,id_H2SO4_strat_,id_BIN01_strat_&
#endif
                               )

    ! transfer information on tracers from dynamics to physics
    USE print_control_mod, ONLY: prt_level, lunout
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nqtot_
    INTEGER,INTENT(IN) :: nqo_
    INTEGER,INTENT(IN) :: nbtr_
#ifdef CPP_StratAer
    INTEGER,INTENT(IN) :: nbtr_bin_
    INTEGER,INTENT(IN) :: nbtr_sulgas_
    INTEGER,INTENT(IN) :: id_OCS_strat_
    INTEGER,INTENT(IN) :: id_SO2_strat_
    INTEGER,INTENT(IN) :: id_H2SO4_strat_
    INTEGER,INTENT(IN) :: id_BIN01_strat_
#endif
    CHARACTER(len=20),INTENT(IN) :: tname_(nqtot_) ! tracer short name for restart and diagnostics
    CHARACTER(len=23),INTENT(IN) :: ttext_(nqtot_) ! tracer long name for diagnostics
    CHARACTER(len=4),INTENT(IN) :: type_trac_
    INTEGER,INTENT(IN) :: niadv_ (nqtot_) ! equivalent dyn / physique
    INTEGER,INTENT(IN) :: conv_flg_(nbtr_)
    INTEGER,INTENT(IN) :: pbl_flg_(nbtr_)
    CHARACTER(len=8),INTENT(IN) :: solsym_(nbtr_)
    ! Isotopes:
    INTEGER,INTENT(IN) :: nqfils_(nqtot_)
    INTEGER,INTENT(IN) :: nqdesc_(nqtot_)
    INTEGER,INTENT(IN) :: nqdesc_tot_
    INTEGER,INTENT(IN) :: iqfils_(nqtot_,nqtot_)
    INTEGER,INTENT(IN) :: iqpere_(nqtot_)
    LOGICAL,INTENT(IN) :: ok_isotopes_
    LOGICAL,INTENT(IN) :: ok_iso_verif_
    LOGICAL,INTENT(IN) :: ok_isotrac_
    LOGICAL,INTENT(IN) :: ok_init_iso_
    INTEGER,INTENT(IN) :: niso_possibles_
    REAL,INTENT(IN) :: tnat_(niso_possibles_)
    REAL,INTENT(IN) :: alpha_ideal_(niso_possibles_)
    LOGICAL,INTENT(IN) :: use_iso_(niso_possibles_)
    INTEGER,INTENT(IN) :: iqiso_(ntraciso_,nqo_)
    INTEGER,INTENT(IN) :: iso_num_(nqtot_)
    INTEGER,INTENT(IN) :: iso_indnum_(nqtot_)
    INTEGER,INTENT(IN) :: zone_num_(nqtot_)
    INTEGER,INTENT(IN) :: phase_num_(nqtot_)
    INTEGER,INTENT(IN) :: indnum_fn_num_(niso_possibles_)
    INTEGER,INTENT(IN) :: index_trac_(ntraceurs_zone_,niso_)
    INTEGER,INTENT(IN) :: niso_
    INTEGER,INTENT(IN) :: ntraceurs_zone_
    INTEGER,INTENT(IN) :: ntraciso_

    CHARACTER(LEN=30) :: modname="init_infotrac_phy"

    nqtot=nqtot_
    nqo=nqo_
    nbtr=nbtr_
#ifdef CPP_StratAer
    nbtr_bin=nbtr_bin_
    nbtr_sulgas=nbtr_sulgas_
    id_OCS_strat=id_OCS_strat_
    id_SO2_strat=id_SO2_strat_
    id_H2SO4_strat=id_H2SO4_strat_
    id_BIN01_strat=id_BIN01_strat_
#endif
    ALLOCATE(tname(nqtot))
    tname(:) = tname_(:)
    ALLOCATE(ttext(nqtot))
    ttext(:) = ttext_(:)
    type_trac = type_trac_
    ALLOCATE(niadv(nqtot))
    niadv(:)=niadv_(:)
    ALLOCATE(conv_flg(nbtr))
    conv_flg(:)=conv_flg_(:)
    ALLOCATE(pbl_flg(nbtr))
    pbl_flg(:)=pbl_flg_(:)
    ALLOCATE(solsym(nbtr))
    solsym(:)=solsym_(:)
  
    IF(prt_level.ge.1) THEN
      write(lunout,*) TRIM(modname)//": nqtot,nqo,nbtr",nqtot,nqo,nbtr
    ENDIF
    
    ! Isotopes:
    
    ! First check that the "niso_possibles" has the correct value
    IF (niso_possibles.ne.niso_possibles_) THEN
      CALL abort_physic(modname,&
           "wrong value for parameter niso_possibles in infotrac_phy",1)
    ENDIF
    
    ok_isotopes=ok_isotopes_
    ok_iso_verif=ok_iso_verif_
    ok_isotrac=ok_isotrac_
    ok_init_iso=ok_init_iso_
    
    niso=niso_
    ntraceurs_zone=ntraceurs_zone_
    ntraciso=ntraciso_
    
    IF (ok_isotopes) THEN
      ALLOCATE(nqfils(nqtot))
      nqfils(:)=nqfils_(:)
      ALLOCATE(nqdesc(nqtot))
      nqdesc(:)=nqdesc_(:)
      nqdesc_tot=nqdesc_tot_
      ALLOCATE(iqfils(nqtot,nqtot))
      iqfils(:,:)=iqfils_(:,:)
      ALLOCATE(iqpere(nqtot))
      iqpere(:)=iqpere_(:)
    
      tnat(:)=tnat_(:)
      alpha_ideal(:)=alpha_ideal_(:)
      use_iso(:)=use_iso_(:)
    
      ALLOCATE(iqiso(ntraciso,nqo))
      iqiso(:,:)=iqiso_(:,:)
      ALLOCATE(iso_num(nqtot))
      iso_num(:)=iso_num_(:)
      ALLOCATE(iso_indnum(nqtot))
      iso_indnum(:)=iso_indnum_(:)
      ALLOCATE(zone_num(nqtot))
      zone_num(:)=zone_num_(:)
      ALLOCATE(phase_num(nqtot))
      phase_num(:)=phase_num_(:)
      
      indnum_fn_num(:)=indnum_fn_num_(:)
      
      ALLOCATE(index_trac(ntraceurs_zone,niso))
      index_trac(:,:)=index_trac_(:,:)
    ENDIF ! of IF(ok_isotopes)
  
  END SUBROUTINE init_infotrac_phy

END MODULE infotrac_phy
