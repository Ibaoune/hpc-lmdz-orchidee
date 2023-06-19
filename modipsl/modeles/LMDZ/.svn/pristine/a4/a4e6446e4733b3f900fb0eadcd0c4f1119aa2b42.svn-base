!
! $Id: iniphysiq.F 1403 2010-07-01 09:02:53Z fairhead $
!
MODULE iniphysiq_mod

CONTAINS

SUBROUTINE iniphysiq(ii,jj,nlayer, &
                     nbp, communicator, &
                     punjours, pdayref,ptimestep, &
                     rlatu,rlatv,rlonu,rlonv,aire,cu,cv,         &
                     prad,pg,pr,pcpp,iflag_phys)
  USE dimphy, ONLY: init_dimphy
  USE inigeomphy_mod, ONLY: inigeomphy
  USE vertical_layers_mod, ONLY : init_vertical_layers
  USE infotrac, ONLY: nqtot
  USE comcstphy, ONLY: rradius, & ! planet radius (m)
                       rr, & ! recuced gas constant: R/molar mass of atm
                       rg, & ! gravity
                       rcpp  ! specific heat of the atmosphere
  USE infotrac_phy, ONLY: init_infotrac_phy
  USE nrtype, ONLY: pi
  IMPLICIT NONE
  !
  !=======================================================================
  !   Initialisation of the physical constants and some positional and 
  !   geometrical arrays for the physics
  !=======================================================================
 
 
  include "iniprint.h"

  REAL,INTENT(IN) :: prad ! radius of the planet (m)
  REAL,INTENT(IN) :: pg ! gravitational acceleration (m/s2)
  REAL,INTENT(IN) :: pr ! ! reduced gas constant R/mu
  REAL,INTENT(IN) :: pcpp ! specific heat Cp
  REAL,INTENT(IN) :: punjours ! length (in s) of a standard day
  INTEGER, INTENT (IN) :: nlayer ! number of atmospheric layers
  INTEGER, INTENT (IN) :: ii ! number of atmospheric coulumns along longitudes
  INTEGER, INTENT (IN) :: jj  ! number of atompsheric columns along latitudes
  INTEGER, INTENT(IN) :: nbp ! number of physics columns for this MPI process
  INTEGER, INTENT(IN) :: communicator ! MPI communicator
  REAL, INTENT (IN) :: rlatu(jj+1) ! latitudes of the physics grid
  REAL, INTENT (IN) :: rlatv(jj) ! latitude boundaries of the physics grid
  REAL, INTENT (IN) :: rlonv(ii+1) ! longitudes of the physics grid
  REAL, INTENT (IN) :: rlonu(ii+1) ! longitude boundaries of the physics grid
  REAL, INTENT (IN) :: aire(ii+1,jj+1) ! area of the dynamics grid (m2)
  REAL, INTENT (IN) :: cu((ii+1)*(jj+1)) ! cu coeff. (u_covariant = cu * u)
  REAL, INTENT (IN) :: cv((ii+1)*jj) ! cv coeff. (v_covariant = cv * v)
  INTEGER, INTENT (IN) :: pdayref ! reference day of for the simulation
  REAL,INTENT(IN) :: ptimestep !physics time step (s)
  INTEGER,INTENT(IN) :: iflag_phys ! type of physics to be called

  INTEGER :: ibegin,iend,offset
  INTEGER :: i,j,k
  CHARACTER (LEN=20) :: modname='iniphysiq'
  CHARACTER (LEN=80) :: abort_message


  ! --> initialize physics distribution, global fields and geometry
  ! (i.e. things in phy_common or dynphy_lonlat)
  CALL inigeomphy(ii,jj,nlayer, &
               nbp, communicator, &
               rlatu,rlatv, &
               rlonu,rlonv, &
               aire,cu,cv)

  ! --> now initialize things specific to the phymar physics package
  
!$OMP PARALLEL

  ! Initialize tracer names, numbers, etc. for physics
  CALL init_infotrac_phy(nqtot)

! copy some fundamental parameters to physics
  rradius=prad
  rg=pg
  rr=pr
  rcpp=pcpp

!$OMP END PARALLEL

!      print*,'ATTENTION !!! TRAVAILLER SUR INIPHYSIQ'
!      print*,'CONTROLE DES LATITUDES, LONGITUDES, PARAMETRES ...'

! Additional initializations for aquaplanets 
!!$OMP PARALLEL
!      if (iflag_phys>=100) then
!        call iniaqua(klon_omp,rlatd,rlond,iflag_phys)
!      endif
!!$OMP END PARALLEL

END SUBROUTINE iniphysiq

END MODULE iniphysiq_mod
