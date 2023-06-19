MODULE allocate_field_mod

  INTERFACE allocate_u
    MODULE PROCEDURE allocate1d_u1d,allocate2d_u1d,allocate3d_u1d
  END INTERFACE allocate_u

  INTERFACE switch_u
    MODULE PROCEDURE switch1d_u1d,switch2d_u1d,switch3d_u1d
  END INTERFACE switch_u

  INTERFACE switch_v
    MODULE PROCEDURE switch1d_v1d,switch2d_v1d,switch3d_v1d
  END INTERFACE switch_v

  INTERFACE allocate_v
    MODULE PROCEDURE allocate1d_v1d,allocate2d_v1d,allocate3d_v1d
  END INTERFACE allocate_v

  INTERFACE allocate2d_u
    MODULE PROCEDURE allocate1d_u2d,allocate2d_u2d,allocate3d_u2d
  END INTERFACE allocate2d_u

  INTERFACE allocate2d_v
    MODULE PROCEDURE allocate1d_v2d,allocate2d_v2d,allocate3d_v2d
  END INTERFACE allocate2d_v

  INTERFACE switch2d_u
    MODULE PROCEDURE switch1d_u2d,switch2d_u2d,switch3d_u2d
  END INTERFACE switch2d_u

  INTERFACE switch2d_v
    MODULE PROCEDURE switch1d_v2d,switch2d_v2d,switch3d_v2d
  END INTERFACE switch2D_v

  REAL :: nan

CONTAINS

  SUBROUTINE Init_nan
  IMPLICIT NONE
    REAL*8 :: rnan
    INTEGER :: inan(2)
    EQUIVALENCE(rnan,inan)
    
    inan(1)=2147483647
    inan(2)=2147483647
    
    nan=rnan
  
  END SUBROUTINE Init_nan

  SUBROUTINE allocate1d_u1d(field,d)
  USE parallel_lmdz
  IMPLICIT NONE
  REAL,POINTER :: field(:)
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(d%ijb_u:d%ije_u))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate1d_u1d


  SUBROUTINE allocate2d_u1d(field,dim1,d)
  USE parallel_lmdz
  IMPLICIT NONE
  REAL,POINTER :: field(:,:)
  INTEGER      :: dim1
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(d%ijb_u:d%ije_u,dim1))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate2d_u1d 

  SUBROUTINE allocate3d_u1d(field,dim1,dim2,d)
  USE parallel_lmdz
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:)
  INTEGER      :: dim1,dim2
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(d%ijb_u:d%ije_u,dim1,dim2))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate3d_u1d   



  SUBROUTINE allocate1d_v1d(field,d)
  USE parallel_lmdz
  IMPLICIT NONE
  REAL,POINTER :: field(:)
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(d%ijb_v:d%ije_v))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate1d_v1d


  SUBROUTINE allocate2d_v1d(field,dim1,d)
  USE parallel_lmdz
  IMPLICIT NONE
  REAL,POINTER :: field(:,:)
  INTEGER      :: dim1
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(d%ijb_v:d%ije_v,dim1))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate2d_v1d 

  SUBROUTINE allocate3d_v1d(field,dim1,dim2,d)
  USE parallel_lmdz
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:)
  INTEGER      :: dim1,dim2
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(d%ijb_v:d%ije_v,dim1,dim2))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate3d_v1d   









  SUBROUTINE allocate1d_u2d(field,d)
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:)
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(iip1,d%jjb_u:d%jje_u))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate1d_u2d


  SUBROUTINE allocate2d_u2d(field,dim1,d)
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:)
  INTEGER      :: dim1
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(iip1,d%jjb_u:d%jje_u,dim1))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate2d_u2d

  SUBROUTINE allocate3d_u2d(field,dim1,dim2,d)
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:,:)
  INTEGER      :: dim1,dim2
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(iip1,d%jjb_u:d%jje_u,dim1,dim2))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate3d_u2d   



  SUBROUTINE allocate1d_v2d(field,d)
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:)
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(iip1,d%jjb_v:d%jje_v))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate1d_v2d


  SUBROUTINE allocate2d_v2d(field,dim1,d)
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:)
  INTEGER      :: dim1
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(iip1,d%jjb_v:d%jje_v,dim1))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate2d_v2d

  SUBROUTINE allocate3d_v2d(field,dim1,dim2,d)
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:,:)
  INTEGER      :: dim1,dim2
  TYPE(distrib),INTENT(IN) :: d

!$OMP BARRIER
!$OMP MASTER    
    IF (ASSOCIATED(field)) DEALLOCATE(field)
    ALLOCATE(field(iip1,d%jjb_v:d%jje_v,dim1,dim2))
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE allocate3d_v2d   



















  SUBROUTINE switch1d_u1d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  IMPLICIT NONE
  REAL,POINTER :: field(:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down
  
  REAL,POINTER,SAVE :: new_field(:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(new_dist%ijb_u:new_dist%ije_u))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField_u(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER
   
    CALL barrier
  END SUBROUTINE switch1d_u1d  
  
  SUBROUTINE switch2d_u1d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  IMPLICIT NONE
  REAL,POINTER :: field(:,:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:,:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(new_dist%ijb_u:new_dist%ije_u,size(field,2)))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField_u(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER
    CALL barrier

  END SUBROUTINE switch2d_u1d  

  SUBROUTINE switch3d_u1d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:,:,:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(new_dist%ijb_u:new_dist%ije_u,size(field,2),size(field,3)))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField_u(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER
    CALL barrier

  END SUBROUTINE switch3d_u1d  




  SUBROUTINE switch1d_v1d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  IMPLICIT NONE
  REAL,POINTER :: field(:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(new_dist%ijb_v:new_dist%ije_v))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField_v(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER

    CALL barrier
  END SUBROUTINE switch1d_v1d  
  
  SUBROUTINE switch2d_v1d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  IMPLICIT NONE
  REAL,POINTER :: field(:,:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:,:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(new_dist%ijb_v:new_dist%ije_v,size(field,2)))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField_v(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER
    CALL barrier

  END SUBROUTINE switch2d_v1d  

  SUBROUTINE switch3d_v1d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:,:,:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(new_dist%ijb_v:new_dist%ije_v,size(field,2),size(field,3)))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField_v(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER
    CALL barrier

  END SUBROUTINE switch3d_v1d  












  SUBROUTINE switch1d_u2d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:,:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(iip1,new_dist%jjb_u:new_dist%jje_u))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField2d_u(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER
    CALL barrier

  END SUBROUTINE switch1d_u2d  
  
  SUBROUTINE switch2d_u2d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:,:,:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(iip1,new_dist%jjb_u:new_dist%jje_u,size(field,3)))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField2d_u(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER
    CALL barrier

  END SUBROUTINE switch2d_u2d  

  SUBROUTINE switch3d_u2d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:,:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:,:,:,:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(iip1,new_dist%jjb_u:new_dist%jje_u,size(field,3),size(field,4)))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField2d_u(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER
    CALL barrier

  END SUBROUTINE switch3d_u2d  




  SUBROUTINE switch1d_v2d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:,:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(iip1,new_dist%jjb_v:new_dist%jje_v))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField2d_v(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER
    CALL barrier

  END SUBROUTINE switch1d_v2d  
  
  SUBROUTINE switch2d_v2d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:,:,:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(iip1,new_dist%jjb_v:new_dist%jje_v,size(field,3)))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField2d_v(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER

    CALL barrier
  END SUBROUTINE switch2d_v2d  

  SUBROUTINE switch3d_v2d(field,old_dist,new_dist,up,down)
  USE parallel_lmdz
  USE mod_hallo
  USE dimensions_mod
  IMPLICIT NONE
  REAL,POINTER :: field(:,:,:,:)
  TYPE(distrib),INTENT(IN) :: old_dist
  TYPE(distrib),INTENT(IN) :: new_dist
  INTEGER, OPTIONAL,INTENT(IN) :: up
  INTEGER, OPTIONAL,INTENT(IN) :: down

  REAL,POINTER,SAVE :: new_field(:,:,:,:)
  TYPE(request) :: req
  
  !$OMP BARRIER
  !$OMP MASTER    
    ALLOCATE(new_field(iip1,new_dist%jjb_v:new_dist%jje_v,size(field,3),size(field,4)))
    new_field=nan
  !$OMP END MASTER
  !$OMP BARRIER
    CALL Register_SwapField2d_v(field,new_field,new_dist,req,old_dist=old_dist,up=up,down=down)
  
    CALL SendRequest(req)

  !$OMP BARRIER
    CALL WaitRequest(req)     
  !$OMP BARRIER
    
  !$OMP MASTER
    DEALLOCATE(field)
    field=>new_field
  !$OMP END MASTER
  !$OMP BARRIER

    CALL barrier
  END SUBROUTINE switch3d_v2d 

END MODULE allocate_field_mod
  
  
  
  
