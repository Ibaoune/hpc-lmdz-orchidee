!
! $Id: iniphysiq_mod.F90 3677 2020-05-06 15:18:32Z oboucher $
!
MODULE iniphysiq_mod

CONTAINS

SUBROUTINE iniphysiq(ii,jj,nlayer, &
                     nbp, communicator, &
                     punjours, pdayref,ptimestep, &
                     rlatudyn,rlatvdyn,rlonudyn,rlonvdyn,airedyn,cudyn,cvdyn, &
                     prad,pg,pr,pcpp,iflag_phys)
  USE dimphy, ONLY: init_dimphy
  USE inigeomphy_mod, ONLY: inigeomphy
  USE mod_grid_phy_lmdz, ONLY: nbp_lon,nbp_lat,nbp_lev,klon_glo ! number of atmospheric columns (on full grid)
  USE mod_phys_lmdz_para, ONLY: klon_omp ! number of columns (on local omp grid)
  USE vertical_layers_mod, ONLY : init_vertical_layers
  USE infotrac, ONLY: nqtot,nqo,nbtr,tname,ttext,type_trac,&
                      niadv,conv_flg,pbl_flg,solsym,&
                      nqfils,nqdesc,nqdesc_tot,iqfils,iqpere,&
                      ok_isotopes,ok_iso_verif,ok_isotrac,&
                      ok_init_iso,niso_possibles,tnat,&
                      alpha_ideal,use_iso,iqiso,iso_num,&
                      iso_indnum,zone_num,phase_num,&
                      indnum_fn_num,index_trac,&
                      niso,ntraceurs_zone,ntraciso
#ifdef CPP_StratAer
  USE infotrac, ONLY: nbtr_bin, nbtr_sulgas, id_OCS_strat, &
                      id_SO2_strat, id_H2SO4_strat, id_BIN01_strat

#endif
#ifdef REPROBUS
  USE CHEM_REP, ONLY : Init_chem_rep_phys
#ifdef CPP_PARA
  USE parallel_lmdz, ONLY : mpi_size, mpi_rank
  USE bands, ONLY : distrib_phys
#endif
  USE mod_phys_lmdz_omp_data, ONLY: klon_omp
#endif
  USE control_mod, ONLY: dayref,anneeref,day_step,nday,offline, iphysiq, config_inca
  USE inifis_mod, ONLY: inifis
  USE time_phylmdz_mod, ONLY: init_time
  USE temps_mod, ONLY: annee_ref, day_ini, day_ref, start_time, calend, year_len
  USE infotrac_phy, ONLY: init_infotrac_phy
  USE phystokenc_mod, ONLY: init_phystokenc
  USE phyaqua_mod, ONLY: iniaqua
  USE comconst_mod, ONLY: omeg, rad
#ifdef INCA
  USE indice_sol_mod, ONLY: nbsrf, is_oce, is_sic, is_ter, is_lic
#ifdef CPP_PARA
  USE parallel_lmdz, ONLY : mpi_size, mpi_rank
  USE bands, ONLY : distrib_phys
#endif
  USE mod_phys_lmdz_omp_data, ONLY: klon_omp 
#endif
  USE ioipsl_getin_p_mod, ONLY: getin_p
  USE slab_heat_transp_mod, ONLY: ini_slab_transp_geom
  IMPLICIT NONE

  ! =======================================================================
  ! Initialisation of the physical constants and some positional and
  ! geometrical arrays for the physics
  ! =======================================================================

  include "dimensions.h"
  include "paramet.h"
  include "iniprint.h"
  include "tracstoke.h"
  include "comgeom.h"

  REAL, INTENT (IN) :: prad ! radius of the planet (m)
  REAL, INTENT (IN) :: pg ! gravitational acceleration (m/s2)
  REAL, INTENT (IN) :: pr ! ! reduced gas constant R/mu
  REAL, INTENT (IN) :: pcpp ! specific heat Cp
  REAL, INTENT (IN) :: punjours ! length (in s) of a standard day
  INTEGER, INTENT (IN) :: nlayer ! number of atmospheric layers
  INTEGER, INTENT (IN) :: ii ! number of atmospheric columns along longitudes
  INTEGER, INTENT (IN) :: jj ! number of atompsheric columns along latitudes
  INTEGER, INTENT(IN) :: nbp ! number of physics columns for this MPI process
  INTEGER, INTENT(IN) :: communicator ! MPI communicator
  REAL, INTENT (IN) :: rlatudyn(jj+1) ! latitudes of the physics grid
  REAL, INTENT (IN) :: rlatvdyn(jj) ! latitude boundaries of the physics grid
  REAL, INTENT (IN) :: rlonvdyn(ii+1) ! longitudes of the physics grid
  REAL, INTENT (IN) :: rlonudyn(ii+1) ! longitude boundaries of the physics grid
  REAL, INTENT (IN) :: airedyn(ii+1,jj+1) ! area of the dynamics grid (m2)
  REAL, INTENT (IN) :: cudyn((ii+1)*(jj+1)) ! cu coeff. (u_covariant = cu * u)
  REAL, INTENT (IN) :: cvdyn((ii+1)*jj) ! cv coeff. (v_covariant = cv * v)
  INTEGER, INTENT (IN) :: pdayref ! reference day of for the simulation
  REAL, INTENT (IN) :: ptimestep !physics time step (s)
  INTEGER, INTENT (IN) :: iflag_phys ! type of physics to be called

  INTEGER :: ibegin, iend, offset
  INTEGER :: i,j,k
  CHARACTER (LEN=20) :: modname = 'iniphysiq'
  CHARACTER (LEN=80) :: abort_message
  
  LOGICAL :: slab_hdiff
  INTEGER :: slab_ekman
  CHARACTER (LEN = 6) :: type_ocean 

#ifndef CPP_PARA
  INTEGER,PARAMETER :: mpi_rank=0
  INTEGER, PARAMETER :: mpi_size = 1 
  INTEGER :: distrib_phys(mpi_rank:mpi_rank)=(jjm-1)*iim+2
#endif

  ! --> initialize physics distribution, global fields and geometry
  ! (i.e. things in phy_common or dynphy_lonlat)
  CALL inigeomphy(ii,jj,nlayer, &
               nbp, communicator, &
               rlatudyn,rlatvdyn, &
               rlonudyn,rlonvdyn, &
               airedyn,cudyn,cvdyn)

  ! --> now initialize things specific to the phylmd physics package
  
!!$OMP PARALLEL DEFAULT(SHARED) COPYIN(/temps/)
!$OMP PARALLEL DEFAULT(SHARED) &
!	Copy all threadprivate variables in temps_mod
!$OMP COPYIN(annee_ref, day_ini, day_ref, start_time)

  ! Initialize physical constants in physics:
  CALL inifis(punjours,prad,pg,pr,pcpp)

  CALL init_time(annee_ref,day_ref,day_ini,start_time,nday,ptimestep)

  ! Initialize dimphy module (unless in 1D where it has already been done)
!  IF (klon_glo>1) CALL Init_dimphy(klon_omp,nlayer)

  ! Copy over "offline" settings
  CALL init_phystokenc(offline,istphy)

  ! Initialization for slab heat transport
  type_ocean="force"
  CALL getin_p('type_ocean',type_ocean)
  slab_hdiff=.FALSE.
  CALL getin_p('slab_hdiff',slab_hdiff)
  slab_ekman=0
  CALL getin_p('slab_ekman',slab_ekman)
  IF ((type_ocean=='slab').AND.(slab_hdiff.OR.(slab_ekman.GT.0))) THEN
     CALL ini_slab_transp_geom(ip1jm,ip1jmp1,unsairez,fext,unsaire,&
                                  cu,cuvsurcv,cv,cvusurcu, &
                                  aire,apoln,apols, &
                                  aireu,airev,rlatvdyn,rad,omeg) 
  END IF

  ! Initialize tracer names, numbers, etc. for physics
  CALL init_infotrac_phy(nqtot,nqo,nbtr,tname,ttext,type_trac,&
                         niadv,conv_flg,pbl_flg,solsym,&
                         nqfils,nqdesc,nqdesc_tot,iqfils,iqpere,&
                         ok_isotopes,ok_iso_verif,ok_isotrac,&
                         ok_init_iso,niso_possibles,tnat,&
                         alpha_ideal,use_iso,iqiso,iso_num,&
                         iso_indnum,zone_num,phase_num,&
                         indnum_fn_num,index_trac,&
                         niso,ntraceurs_zone,ntraciso&
#ifdef CPP_StratAer
                         ,nbtr_bin,nbtr_sulgas&
                         ,id_OCS_strat,id_SO2_strat,id_H2SO4_strat,id_BIN01_strat&
#endif
                         )

  ! Initializations for Reprobus
  IF (type_trac == 'repr') THEN
#ifdef REPROBUS
    call Init_chem_rep_phys(klon_omp,nlayer)
    call init_reprobus_para( &
          nbp_lon,nbp_lat,nbp_lev,klon_glo,mpi_size, &
          distrib_phys,communicator)
#endif
  ENDIF
!$OMP END PARALLEL

  IF (type_trac == 'inca') THEN
#ifdef INCA
     call init_const_lmdz( &
          anneeref,dayref, iphysiq,day_step,nday,  &
          nbsrf, is_oce,is_sic, is_ter,is_lic, calend, &
          config_inca)
     call init_inca_para( &
          nbp_lon,nbp_lat,nbp_lev,klon_glo,mpi_size, &
          distrib_phys,communicator)
#endif
  END IF
  IF (type_trac == 'repr') THEN
#ifdef REPROBUS
    call init_reprobus_para( &
          nbp_lon,nbp_lat,nbp_lev,klon_glo,mpi_size, &
          distrib_phys,communicator)
#endif
  ENDIF

!!$OMP PARALLEL DEFAULT(SHARED) COPYIN(/temps/)
!$OMP PARALLEL DEFAULT(SHARED)
  ! Additional initializations for aquaplanets 
  IF (iflag_phys>=100) THEN
    CALL iniaqua(klon_omp,year_len,iflag_phys)
  END IF

  IF (type_trac == 'inca') THEN
#ifdef INCA
     CALL init_inca_dim(klon_omp,nbp_lev,nbp_lon,nbp_lat - 1, &
          rlonudyn,rlatudyn,rlonvdyn,rlatvdyn)
#endif
    IF (type_trac == 'repr') THEN
#ifdef REPROBUS
       CALL Init_chem_rep_phys(klon_omp,nbp_lev)
#endif
    END IF
  END IF

!$OMP END PARALLEL

END SUBROUTINE iniphysiq

END MODULE iniphysiq_mod
