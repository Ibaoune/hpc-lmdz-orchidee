!
! $Id: dimphy.F90 3442 2019-01-25 10:12:29Z fairhead $
!
  MODULE dimphy
  
  INTEGER,SAVE :: klon
  INTEGER,SAVE :: kdlon
  INTEGER,SAVE :: kfdia
  INTEGER,SAVE :: kidia
  INTEGER,SAVE :: klev
  INTEGER,SAVE :: klevp1
  INTEGER,SAVE :: klevm1
  INTEGER,SAVE :: kflev

!$OMP THREADPRIVATE(klon,kfdia,kidia,kdlon)
  REAL,save,allocatable,dimension(:) :: zmasq
!$OMP THREADPRIVATE(zmasq)   

CONTAINS
  
  SUBROUTINE Init_dimphy(klon0,klev0)
  IMPLICIT NONE
  
    INTEGER, INTENT(in) :: klon0
    INTEGER, INTENT(in) :: klev0
    
    klon=klon0
    kdlon=klon
    kidia=1
    kfdia=klon
!$OMP MASTER 
    klev=klev0
    klevp1=klev+1
    klevm1=klev-1
    kflev=klev
!$OMP END MASTER    
    ALLOCATE(zmasq(klon))    
    zmasq=0.
    
  END SUBROUTINE Init_dimphy

  SUBROUTINE Init_dimphy1D(klon0,klev0)
! 1D special version of dimphy without ALLOCATE(zmasq)
! which will be allocated in iniphysiq
  IMPLICIT NONE
  
    INTEGER, INTENT(in) :: klon0
    INTEGER, INTENT(in) :: klev0
    
    klon=klon0
    kdlon=klon
    kidia=1
    kfdia=klon
    klev=klev0
    klevp1=klev+1
    klevm1=klev-1
    kflev=klev
    
  END SUBROUTINE Init_dimphy1D

  
END MODULE dimphy
