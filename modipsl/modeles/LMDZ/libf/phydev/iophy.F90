!
! $Id: iophy.F90 2588 2016-07-13 06:54:39Z emillour $
!
module iophy
  
! abd  REAL,private,allocatable,dimension(:),save :: io_lat
! abd  REAL,private,allocatable,dimension(:),save :: io_lon
  REAL,allocatable,dimension(:),save :: io_lat
  REAL,allocatable,dimension(:),save :: io_lon
  INTEGER, save :: phys_domain_id
  INTEGER, save :: npstn
  INTEGER, allocatable, dimension(:), save :: nptabij
  

#ifdef CPP_XIOS
! interfaces for both IOIPSL and XIOS
  INTERFACE histwrite_phy
    MODULE PROCEDURE histwrite2d_phy,histwrite3d_phy,histwrite2d_xios,histwrite3d_xios
  END INTERFACE
#else
! interfaces for IOIPSL
  INTERFACE histwrite_phy
    MODULE PROCEDURE histwrite2d_phy,histwrite3d_phy
  END INTERFACE
#endif

#ifdef CPP_XIOS
! interfaces for both IOIPSL and XIOS
  INTERFACE histbeg_phy_all
    MODULE PROCEDURE histbeg_phy, histbeg_phyxios
  END INTERFACE
#else
! interfaces for IOIPSL
  INTERFACE histbeg_phy_all
    MODULE PROCEDURE histbeg_phy
  END INTERFACE
#endif

contains

  subroutine init_iophy_new(rlat,rlon)
  USE dimphy, only: klon
  USE mod_phys_lmdz_para, only: gather, bcast, &
                                jj_nb, jj_begin, jj_end, ii_begin, ii_end, &
                                mpi_size, mpi_rank, klon_mpi, &
                                is_sequential, is_south_pole_dyn
  USE mod_grid_phy_lmdz, only: nbp_lon, nbp_lat, klon_glo
  USE print_control_mod, ONLY: lunout, prt_level
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
#ifdef CPP_IOIPSL
  USE ioipsl, only: flio_dom_set
#endif
#ifdef CPP_XIOS
  use wxios, only: wxios_domain_param
#endif
  implicit none
    real,dimension(klon),intent(in) :: rlon
    real,dimension(klon),intent(in) :: rlat

    REAL,dimension(klon_glo)        :: rlat_glo
    REAL,dimension(klon_glo)        :: rlon_glo
    
    INTEGER,DIMENSION(2) :: ddid
    INTEGER,DIMENSION(2) :: dsg
    INTEGER,DIMENSION(2) :: dsl
    INTEGER,DIMENSION(2) :: dpf
    INTEGER,DIMENSION(2) :: dpl
    INTEGER,DIMENSION(2) :: dhs
    INTEGER,DIMENSION(2) :: dhe 
    INTEGER :: i    
    integer :: data_ibegin,data_iend

    CALL gather(rlat,rlat_glo)
    CALL bcast(rlat_glo)
    CALL gather(rlon,rlon_glo)
    CALL bcast(rlon_glo)
    
!$OMP MASTER  
    ALLOCATE(io_lat(nbp_lat))
    io_lat(1)=rlat_glo(1)
    io_lat(nbp_lat)=rlat_glo(klon_glo)
    IF ((nbp_lon*nbp_lat) > 1) then
      DO i=2,nbp_lat-1
        io_lat(i)=rlat_glo(2+(i-2)*nbp_lon)
      ENDDO
    ENDIF

    ALLOCATE(io_lon(nbp_lon))
    IF ((nbp_lon*nbp_lat) > 1) THEN
      io_lon(:)=rlon_glo(2:nbp_lon+1)
    ELSE
      io_lon(1)=rlon_glo(1)
    ENDIF
!! (I) dtnb   : total number of domains
!! (I) dnb    : domain number
!! (I) did(:) : distributed dimensions identifiers
!!              (up to 5 dimensions are supported)
!! (I) dsg(:) : total number of points for each dimension
!! (I) dsl(:) : local number of points for each dimension
!! (I) dpf(:) : position of first local point for each dimension
!! (I) dpl(:) : position of last local point for each dimension
!! (I) dhs(:) : start halo size for each dimension
!! (I) dhe(:) : end halo size for each dimension
!! (C) cdnm   : Model domain definition name.
!!              The names actually supported are :
!!              "BOX", "APPLE", "ORANGE".
!!              These names are case insensitive.
    ddid=(/ 1,2 /)
    dsg=(/ nbp_lon, nbp_lat /)
    dsl=(/ nbp_lon, jj_nb /)
    dpf=(/ 1,jj_begin /)
    dpl=(/ nbp_lon, jj_end /)
    dhs=(/ ii_begin-1,0 /)
    if (mpi_rank==mpi_size-1) then
      dhe=(/0,0/)
    else
      dhe=(/ nbp_lon-ii_end,0 /)  
    endif
    
#ifndef CPP_IOIPSL_NO_OUTPUT
    call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, &
                      'APPLE',phys_domain_id)
#endif
#ifdef CPP_XIOS
    ! Set values for the mask:
    IF (mpi_rank == 0) THEN
        data_ibegin = 0
    ELSE 
        data_ibegin = ii_begin - 1
    END IF

    IF (mpi_rank == mpi_size-1) THEN
        data_iend = nbp_lon
    ELSE
        data_iend = ii_end + 1
    END IF

    if (prt_level>=10) then
      write(lunout,*) "init_iophy_new: mpirank=",mpi_rank," iibegin=",ii_begin , " ii_end=",ii_end," jjbegin=",jj_begin," jj_nb=",jj_nb," jj_end=",jj_end
      write(lunout,*) "init_iophy_new: mpirank=",mpi_rank," nbp_lon=",nbp_lon," nbp_lat=",nbp_lat
      write(lunout,*) "init_iophy_new: mpirank=",mpi_rank," data_ibegin=",data_ibegin," data_iend=",data_iend
      write(lunout,*) "init_iophy_new: mpirank=",mpi_rank," data_ibegin=",data_ibegin," data_iend=",data_iend
      write(lunout,*) "init_iophy_new: mpirank=",mpi_rank," is_south_pole=",is_south_pole_dyn
    endif

    ! Initialize the XIOS domain coreesponding to this process:
    CALL wxios_domain_param("dom_glo", is_sequential, nbp_lon, jj_nb, nbp_lon, nbp_lat, &
                            1, nbp_lon, ii_begin, ii_end, jj_begin, jj_end,             &
                            klon_mpi+2*(nbp_lon-1), data_ibegin, data_iend,             &
                            io_lat, io_lon,is_south_pole_dyn,mpi_rank)
#endif
!$OMP END MASTER
      
  END SUBROUTINE init_iophy_new
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine histbeg_phy(name,itau0,zjulian,dtime,nhori,nid_day)
  USE mod_phys_lmdz_para, only: is_sequential, jj_begin, jj_end, jj_nb
  use ioipsl, only: histbeg
  USE print_control_mod, ONLY: prt_level, lunout
  USE mod_grid_phy_lmdz, ONLY: nbp_lon
  implicit none
    
    character*(*), intent(IN) :: name
    integer, intent(in) :: itau0
    real,intent(in) :: zjulian
    real,intent(in) :: dtime
    integer,intent(out) :: nhori
    integer,intent(out) :: nid_day

!$OMP MASTER    
    if (is_sequential) then
      call histbeg(name,nbp_lon,io_lon, jj_nb,io_lat(jj_begin:jj_end), &
                   1,nbp_lon,1,jj_nb,itau0, zjulian, dtime, nhori, nid_day)
    else
      call histbeg(name,nbp_lon,io_lon, jj_nb,io_lat(jj_begin:jj_end), &
                   1,nbp_lon,1,jj_nb,itau0, zjulian, dtime, nhori, nid_day,phys_domain_id)
    endif
!$OMP END MASTER
  
  end subroutine histbeg_phy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CPP_XIOS

! SUBROUTINE histbeg_phyxios(name,itau0,zjulian,dtime,ffreq,lev,nhori,nid_day)
 SUBROUTINE histbeg_phyxios(name,ffreq,lev)
  USE mod_phys_lmdz_para, only: is_using_mpi, is_mpi_root
  use wxios, only: wxios_add_file
  IMPLICIT NONE
    
    character*(*), INTENT(IN) :: name
!    integer, INTENT(IN) :: itau0
!    REAL,INTENT(IN) :: zjulian
!    REAL,INTENT(IN) :: dtime
    character(LEN=*), INTENT(IN) :: ffreq
    INTEGER,INTENT(IN) :: lev
!    integer,intent(out) :: nhori
!    integer,intent(out) :: nid_day

!$OMP MASTER    

    ! ug OMP en chantier...
    IF((.NOT. is_using_mpi) .OR. is_mpi_root) THEN
        ! ug Création du fichier
        CALL wxios_add_file(name, ffreq, lev)
    END IF

!$OMP END MASTER
  
  END SUBROUTINE histbeg_phyxios

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  subroutine histwrite2d_phy(nid,lpoint,name,itau,field)
  USE dimphy, only: klon
  USE mod_phys_lmdz_para, only: Gather_omp, grid1Dto2D_mpi, &
                                is_sequential, klon_mpi_begin, klon_mpi_end, &
                                jj_nb, klon_mpi
  USE ioipsl, only: histwrite
  USE mod_grid_phy_lmdz, ONLY: nbp_lon
  implicit none
    
    integer,intent(in) :: nid
    logical,intent(in) :: lpoint 
    character*(*), intent(IN) :: name
    integer, intent(in) :: itau
    real,dimension(:),intent(in) :: field
    REAL,dimension(klon_mpi) :: buffer_omp
    INTEGER, allocatable, dimension(:) :: index2d
    REAL :: Field2d(nbp_lon,jj_nb)

    integer :: ip
    real,allocatable,dimension(:) :: fieldok

    IF (size(field)/=klon) CALL abort_physic('iophy::histwrite2d','Field first dimension not equal to klon',1)
    
    CALL Gather_omp(field,buffer_omp)    
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,Field2d)
    if(.NOT.lpoint) THEN
     ALLOCATE(index2d(nbp_lon*jj_nb))
     ALLOCATE(fieldok(nbp_lon*jj_nb))
     CALL histwrite(nid,name,itau,Field2d,nbp_lon*jj_nb,index2d)
    else
     ALLOCATE(fieldok(npstn))
     ALLOCATE(index2d(npstn))

     if(is_sequential) then
!     klon_mpi_begin=1
!     klon_mpi_end=klon
      DO ip=1, npstn
       fieldok(ip)=buffer_omp(nptabij(ip))
      ENDDO
     else
      DO ip=1, npstn
!     print*,'histwrite2d is_sequential npstn ip name nptabij',npstn,ip,name,nptabij(ip)
       IF(nptabij(ip).GE.klon_mpi_begin.AND. &
          nptabij(ip).LE.klon_mpi_end) THEN
         fieldok(ip)=buffer_omp(nptabij(ip)-klon_mpi_begin+1)
       ENDIF
      ENDDO
     endif
     CALL histwrite(nid,name,itau,fieldok,npstn,index2d)
!
    endif
    deallocate(index2d)
    deallocate(fieldok)
!$OMP END MASTER    
  end subroutine histwrite2d_phy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine histwrite3d_phy(nid,lpoint,name,itau,field)
  USE dimphy, only: klon
  USE mod_phys_lmdz_para, only: Gather_omp, grid1Dto2D_mpi, &
                                is_sequential, klon_mpi_begin, klon_mpi_end, &
                                jj_nb, klon_mpi
  USE ioipsl, only: histwrite
  USE mod_grid_phy_lmdz, ONLY: nbp_lon
  implicit none
    
    integer,intent(in) :: nid
    logical,intent(in) :: lpoint
    character*(*), intent(IN) :: name
    integer, intent(in) :: itau
    real,dimension(:,:),intent(in) :: field  ! --> field(klon,:)
    REAL,dimension(klon_mpi,size(field,2)) :: buffer_omp
    REAL :: Field3d(nbp_lon,jj_nb,size(field,2))
    INTEGER :: ip, n, nlev
    INTEGER, ALLOCATABLE, dimension(:) :: index3d
    real,allocatable, dimension(:,:) :: fieldok

    IF (size(field,1)/=klon) CALL abort_physic('iophy::histwrite3d','Field first dimension not equal to klon',1)
    nlev=size(field,2)

    CALL Gather_omp(field,buffer_omp)
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,field3d)
    if(.NOT.lpoint) THEN
     ALLOCATE(index3d(nbp_lon*jj_nb*nlev))
     ALLOCATE(fieldok(nbp_lon*jj_nb,nlev))
     CALL histwrite(nid,name,itau,Field3d,nbp_lon*jj_nb*nlev,index3d)
    else
      nlev=size(field,2)
      ALLOCATE(index3d(npstn*nlev))
      ALLOCATE(fieldok(npstn,nlev))

      if(is_sequential) then
!      klon_mpi_begin=1
!      klon_mpi_end=klon
       DO n=1, nlev
       DO ip=1, npstn
        fieldok(ip,n)=buffer_omp(nptabij(ip),n)
       ENDDO
       ENDDO
      else
       DO n=1, nlev
       DO ip=1, npstn
        IF(nptabij(ip).GE.klon_mpi_begin.AND. &
         nptabij(ip).LE.klon_mpi_end) THEN
         fieldok(ip,n)=buffer_omp(nptabij(ip)-klon_mpi_begin+1,n)
        ENDIF
       ENDDO
       ENDDO
      endif
      CALL histwrite(nid,name,itau,fieldok,npstn*nlev,index3d)
    endif 
  deallocate(index3d)
  deallocate(fieldok)
!$OMP END MASTER    
  end subroutine histwrite3d_phy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! VERSION DES HISTWRITE DEDIEES AU TOUT-XIOS-XML DEJA UTILISEE DANS PHYDEV
#ifdef CPP_XIOS
  SUBROUTINE histwrite2d_xios(field_name,field)
  USE dimphy, only: klon
  USE mod_phys_lmdz_para, only: gather_omp, grid1Dto2D_mpi, &
                                jj_nb, klon_mpi
  USE xios, only: xios_send_field
  USE print_control_mod, ONLY: prt_level, lunout
  USE mod_grid_phy_lmdz, ONLY: nbp_lon
  IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: field_name
    REAL, DIMENSION(:), INTENT(IN) :: field
      
    REAL,DIMENSION(klon_mpi) :: buffer_omp
    REAL :: Field2d(nbp_lon,jj_nb)

    IF (prt_level >= 10) WRITE(lunout,*)'Begin histrwrite2d_xios ',trim(field_name)

    IF (SIZE(field)/=klon) CALL abort_physic('iophy::histwrite2d_xios','Field first DIMENSION not equal to klon',1)
    
    CALL Gather_omp(field,buffer_omp)    
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,Field2d)
    
    CALL xios_send_field(field_name, Field2d)
!$OMP END MASTER   

    IF (prt_level >= 10) WRITE(lunout,*)'End histrwrite2d_xios ',trim(field_name)
  END SUBROUTINE histwrite2d_xios
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! VERSION DES HISTWRITE DEDIEES AU TOUT-XIOS-XML DEJA UTILISEE DANS PHYDEV
#ifdef CPP_XIOS
  SUBROUTINE histwrite3d_xios(field_name, field)
  USE dimphy, only: klon, klev
  USE mod_phys_lmdz_para, only: gather_omp, grid1Dto2D_mpi, &
                                jj_nb, klon_mpi
  USE xios, only: xios_send_field
  USE print_control_mod, ONLY: prt_level,lunout
  USE mod_grid_phy_lmdz, ONLY: nbp_lon

  IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: field_name
    REAL, DIMENSION(:,:), INTENT(IN) :: field ! --> field(klon,:)

    REAL,DIMENSION(klon_mpi,SIZE(field,2)) :: buffer_omp
    REAL :: Field3d(nbp_lon,jj_nb,SIZE(field,2))
    INTEGER :: ip, n, nlev

  IF (prt_level >= 10) write(lunout,*)'Begin histrwrite3d_xios ',trim(field_name)

    !Et on.... écrit
    IF (SIZE(field,1)/=klon) CALL abort_physic('iophy::histwrite3d','Field first DIMENSION not equal to klon',1)
    nlev=SIZE(field,2)


    CALL Gather_omp(field,buffer_omp)
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,field3d)

    CALL xios_send_field(field_name, Field3d(:,:,1:nlev))
!$OMP END MASTER   

    IF (prt_level >= 10) write(lunout,*)'End histrwrite3d_xios ',trim(field_name)
  END SUBROUTINE histwrite3d_xios
#endif

end module iophy
