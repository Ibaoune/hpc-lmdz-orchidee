module mod_Hallo
USE parallel_lmdz
implicit none
  logical,save :: use_mpi_alloc
  integer, parameter :: MaxProc=512
  integer, parameter :: DefaultMaxBufferSize=1024*1024*100
  integer, SAVE :: MaxBufferSize=0
  integer, parameter :: ListSize=1000
  
  integer,save       :: MaxBufferSize_Used
!$OMP THREADPRIVATE( MaxBufferSize_Used)

   real,save,pointer,dimension(:) :: Buffer
!$OMP THREADPRIVATE(Buffer)

   integer,save,dimension(Listsize) :: Buffer_Pos
   integer,save :: Index_Pos
!$OMP THREADPRIVATE(Buffer_Pos,Index_pos)
   
  type Hallo
    real, dimension(:,:),pointer :: Field
    integer :: offset
    integer :: size
    integer :: NbLevel
    integer :: Stride
  end type Hallo
  
  type request_SR
    integer :: NbRequest=0
    integer :: NbRequestMax=0
    integer :: BufferSize
    integer :: Pos
    integer :: Index 
    type(Hallo), POINTER :: Hallo(:)
    integer :: MSG_Request
  end type request_SR

  type request
    type(request_SR),dimension(0:MaxProc-1) :: RequestSend
    type(request_SR),dimension(0:MaxProc-1) :: RequestRecv
    integer :: tag=1
  end type request
  
   TYPE(distrib),SAVE :: distrib_gather


  INTERFACE Register_SwapField_u
    MODULE PROCEDURE Register_SwapField1d_u,Register_SwapField2d_u1d,Register_SwapField3d_u, &
                     Register_SwapField1d_u_bis,Register_SwapField2d_u1d_bis,Register_SwapField3d_u_bis
  END INTERFACE Register_SwapField_u

  INTERFACE Register_SwapField_v
    MODULE PROCEDURE Register_SwapField1d_v,Register_SwapField2d_v1d,Register_SwapField3d_v,&
                     Register_SwapField1d_v_bis,Register_SwapField2d_v1d_bis,Register_SwapField3d_v_bis
  END INTERFACE Register_SwapField_v

  INTERFACE Register_SwapField2d_u
    MODULE PROCEDURE Register_SwapField1d_u2d,Register_SwapField2d_u2d,Register_SwapField3d_u2d, &
                     Register_SwapField1d_u2d_bis,Register_SwapField2d_u2d_bis,Register_SwapField3d_u2d_bis
  END INTERFACE Register_SwapField2d_u

  INTERFACE Register_SwapField2d_v
    MODULE PROCEDURE Register_SwapField1d_v2d,Register_SwapField2d_v2d,Register_SwapField3d_v2d, &
                     Register_SwapField1d_v2d_bis,Register_SwapField2d_v2d_bis,Register_SwapField3d_v2d_bis
  END INTERFACE Register_SwapField2d_v

  contains

  subroutine Init_mod_hallo
  USE dimensions_mod
  USE IOIPSL
    implicit none
    integer :: jj_nb_gather(0:mpi_size-1)
    
    Index_Pos=1
    Buffer_Pos(Index_Pos)=1
    MaxBufferSize_Used=0
!$OMP MASTER     
    MaxBufferSize=DefaultMaxBufferSize
    CALL getin("mpi_buffer_size",MaxBufferSize)
!$OMP END MASTER
!$OMP BARRIER
    
    IF (use_mpi_alloc .AND. using_mpi) THEN
      CALL create_global_mpi_buffer
    ELSE 
      CALL create_standard_mpi_buffer
    ENDIF
     
!$OMP MASTER     
     jj_nb_gather(:)=0
     jj_nb_gather(0)=jjp1
     
     CALL create_distrib(jj_nb_gather,distrib_gather) 
!$OMP END MASTER
!$OMP BARRIER

  end subroutine init_mod_hallo

  SUBROUTINE create_standard_mpi_buffer
  IMPLICIT NONE
    
    ALLOCATE(Buffer(MaxBufferSize))
    
  END SUBROUTINE create_standard_mpi_buffer
  
  SUBROUTINE create_global_mpi_buffer
  IMPLICIT NONE
#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif  
    POINTER (Pbuffer,MPI_Buffer(MaxBufferSize))
    REAL :: MPI_Buffer
#ifdef CPP_MPI
    INTEGER(KIND=MPI_ADDRESS_KIND) :: BS 
#else
    INTEGER(KIND=8) :: BS
#endif
    INTEGER :: i,ierr

!  Allocation du buffer MPI
      Bs=8*MaxBufferSize
!$OMP CRITICAL (MPI)
#ifdef CPP_MPI
      CALL MPI_ALLOC_MEM(BS,MPI_INFO_NULL,Pbuffer,ierr)
#endif
!$OMP END CRITICAL (MPI)
      DO i=1,MaxBufferSize
	MPI_Buffer(i)=i
      ENDDO
     
      CALL  Associate_buffer(MPI_Buffer)
      
  CONTAINS
     
     SUBROUTINE Associate_buffer(MPI_Buffer)
     IMPLICIT NONE
       REAL,DIMENSION(:),target :: MPI_Buffer  

         Buffer=>MPI_Buffer
 
      END SUBROUTINE  Associate_buffer
                                      
  END SUBROUTINE create_global_mpi_buffer
 
      
  subroutine allocate_buffer(Size,Index,Pos)
  implicit none
    integer :: Size
    integer :: Index
    integer :: Pos

    if (Buffer_pos(Index_pos)+Size>MaxBufferSize_Used) MaxBufferSize_Used=Buffer_pos(Index_pos)+Size  
    if (Buffer_pos(Index_pos)+Size>MaxBufferSize) then
      print *,'STOP :: La taille de MaxBufferSize dans mod_hallo.F90 est trop petite !!!!'
      stop
    endif
    
    if (Index_pos>=ListSize) then
      print *,'STOP :: La taille de ListSize dans mod_hallo.F90 est trop petite !!!!'
      stop
    endif
     
    Pos=Buffer_Pos(Index_Pos)
    Buffer_Pos(Index_pos+1)=Buffer_Pos(Index_Pos)+Size
    Index_Pos=Index_Pos+1
    Index=Index_Pos
    
  end subroutine allocate_buffer
     
  subroutine deallocate_buffer(Index)
  implicit none
    integer :: Index
    
    Buffer_Pos(Index)=-1
    
    do while (Buffer_Pos(Index_Pos)==-1 .and. Index_Pos>1)
      Index_Pos=Index_Pos-1
    end do

  end subroutine deallocate_buffer  
  
  subroutine SetTag(a_request,tag)
  implicit none
    type(request):: a_request
    integer :: tag
    
    a_request%tag=tag
  end subroutine SetTag
  
  
  subroutine New_Hallo(Field,Stride,NbLevel,offset,size,Ptr_request)
    integer :: Stride
    integer :: NbLevel
    integer :: size
    integer :: offset
    real, dimension(Stride,NbLevel),target :: Field
    type(request_SR),pointer :: Ptr_request
    type(Hallo),POINTER :: NewHallos(:),HalloSwitch(:), NewHallo
    
    Ptr_Request%NbRequest=Ptr_Request%NbRequest+1
    IF(Ptr_Request%NbRequestMax==0) THEN
       Ptr_Request%NbRequestMax=10
       ALLOCATE(Ptr_Request%Hallo(Ptr_Request%NbRequestMax))
    ELSE IF ( Ptr_Request%NbRequest > Ptr_Request%NbRequestMax) THEN
      Ptr_Request%NbRequestMax=INT(Ptr_Request%NbRequestMax*1.2)
      ALLOCATE(NewHallos(Ptr_Request%NbRequestMax))
      NewHallos(1:Ptr_Request%NbRequest-1)=Ptr_Request%hallo(1:Ptr_Request%NbRequest-1)
      HalloSwitch=>Ptr_Request%hallo
      Ptr_Request%hallo=>NewHallos
      DEALLOCATE(HalloSwitch)
    ENDIF
    
    NewHallo=>Ptr_Request%hallo(Ptr_Request%NbRequest)
          
    NewHallo%Field=>Field
    NewHallo%Stride=Stride
    NewHallo%NbLevel=NbLevel
    NewHallo%size=size
    NewHallo%offset=offset
    
  end subroutine New_Hallo
  
  subroutine Register_SendField(Field,ij,ll,offset,size,target,a_request)
  USE dimensions_mod
  implicit none

    
      INTEGER :: ij,ll,offset,size,target
      REAL, dimension(ij,ll) :: Field
      type(request),target :: a_request
      type(request_SR),pointer :: Ptr_request

      Ptr_Request=>a_request%RequestSend(target)

      call New_Hallo(Field,ij,ll,offset,size,Ptr_request)
      
   end subroutine Register_SendField      
      
  subroutine Register_RecvField(Field,ij,ll,offset,size,target,a_request)
  USE dimensions_mod
  implicit none

   
      INTEGER :: ij,ll,offset,size,target
      REAL, dimension(ij,ll) :: Field
      type(request),target :: a_request
      type(request_SR),pointer :: Ptr_request

      Ptr_Request=>a_request%RequestRecv(target)
            
      call New_Hallo(Field,ij,ll,offset,size,Ptr_request)

      
   end subroutine Register_RecvField      
  
  subroutine Register_SwapField(FieldS,FieldR,ij,ll,jj_Nb_New,a_request)
  USE dimensions_mod
      implicit none

    
    INTEGER :: ij,ll
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    type(request) :: a_request
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    integer ::i,jje,jjb
    
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo
    
    do i=0,MPI_Size-1
      if (i /= MPI_Rank) then
        jjb=max(jj_begin_new(i),jj_begin)
        jje=min(jj_end_new(i),jj_end)
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_SendField(FieldS,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
        jjb=max(jj_begin_new(MPI_Rank),jj_begin_Para(i))
        jje=min(jj_end_new(MPI_Rank),jj_end_Para(i))
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_RecvField(FieldR,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
      endif
    enddo
    
  end subroutine Register_SwapField    
  

  
  subroutine Register_SwapFieldHallo(FieldS,FieldR,ij,ll,jj_Nb_New,Up,Down,a_request)
  USE dimensions_mod
  
      implicit none
    
    INTEGER :: ij,ll,Up,Down
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    type(request) :: a_request
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    integer ::i,jje,jjb
    
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo
    
    do i=0,MPI_Size-1
      jj_begin_New(i)=max(1,jj_begin_New(i)-Up)
      jj_end_New(i)=min(jjp1,jj_end_new(i)+Down)
    enddo
   
    do i=0,MPI_Size-1
      if (i /= MPI_Rank) then
        jjb=max(jj_begin_new(i),jj_begin)
        jje=min(jj_end_new(i),jj_end)
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_SendField(FieldS,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
        jjb=max(jj_begin_new(MPI_Rank),jj_begin_Para(i))
        jje=min(jj_end_new(MPI_Rank),jj_end_Para(i))
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_RecvField(FieldR,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
      endif
    enddo
    
  end subroutine Register_SwapFieldHallo



  SUBROUTINE Register_SwapField1d_u(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    REAL, DIMENSION(current_dist%ijb_u:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_u:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down

    CALL  Register_SwapField_gen_u(FieldS,FieldR,1,current_dist,new_dist,halo_up,halo_down,a_request)
        
  END SUBROUTINE  Register_SwapField1d_u 

  SUBROUTINE Register_SwapField1d_u_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    TYPE(distrib),INTENT(IN)          :: old_dist
    REAL, DIMENSION(old_dist%ijb_u:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_u:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down

    CALL  Register_SwapField_gen_u(FieldS,FieldR,1,old_dist,new_dist,halo_up,halo_down,a_request)
        
  END SUBROUTINE  Register_SwapField1d_u_bis 


  SUBROUTINE Register_SwapField2d_u1d(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod
    IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    REAL, DIMENSION(current_dist%ijb_u:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_u:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,2)
    
    CALL  Register_SwapField_gen_u(FieldS,FieldR,ll,current_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField2d_u1d

  SUBROUTINE Register_SwapField2d_u1d_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod
    IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    TYPE(distrib),INTENT(IN) :: old_dist
    REAL, DIMENSION(old_dist%ijb_u:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_u:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,2)
    
    CALL  Register_SwapField_gen_u(FieldS,FieldR,ll,old_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField2d_u1d_bis
   

  SUBROUTINE Register_SwapField3d_u(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    REAL, DIMENSION(current_dist%ijb_u:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_u:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,2)*size(FieldS,3)
    
    CALL  Register_SwapField_gen_u(FieldS,FieldR,ll,current_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField3d_u 

  SUBROUTINE Register_SwapField3d_u_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    TYPE(distrib),INTENT(IN) :: old_dist
    REAL, DIMENSION(old_dist%ijb_u:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_u:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,2)*size(FieldS,3)
    
    CALL  Register_SwapField_gen_u(FieldS,FieldR,ll,old_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField3d_u_bis 
  


 SUBROUTINE Register_SwapField1d_u2d(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod

      IMPLICIT NONE

    TYPE(distrib),INTENT(IN)          :: new_dist !LF
    REAL, DIMENSION(current_dist%jjb_u:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_u:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down

    CALL  Register_SwapField_gen_u(FieldS,FieldR,1,current_dist,new_dist,halo_up,halo_down,a_request)
        
  END SUBROUTINE  Register_SwapField1d_u2d 

 SUBROUTINE Register_SwapField1d_u2d_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod

      IMPLICIT NONE

    TYPE(distrib),INTENT(IN)          :: new_dist !LF
    TYPE(distrib),INTENT(IN)          :: old_dist
    REAL, DIMENSION(old_dist%jjb_u:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_u:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down

    CALL  Register_SwapField_gen_u(FieldS,FieldR,1,old_dist,new_dist,halo_up,halo_down,a_request)
        
  END SUBROUTINE  Register_SwapField1d_u2d_bis 


  SUBROUTINE Register_SwapField2d_u2d(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod

      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    REAL, DIMENSION(current_dist%jjb_u:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_u:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,3)
    
    CALL  Register_SwapField_gen_u(FieldS,FieldR,ll,current_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField2d_u2d

  SUBROUTINE Register_SwapField2d_u2d_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod

      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    TYPE(distrib),INTENT(IN) :: old_dist
    REAL, DIMENSION(old_dist%jjb_u:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_u:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,3)
    
    CALL  Register_SwapField_gen_u(FieldS,FieldR,ll,old_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField2d_u2d_bis
   

  SUBROUTINE Register_SwapField3d_u2d(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    REAL, DIMENSION(current_dist%jjb_u:,:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_u:,:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,3)*size(FieldS,4)
    
    CALL  Register_SwapField_gen_u(FieldS,FieldR,ll,current_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField3d_u2d 

  SUBROUTINE Register_SwapField3d_u2d_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    TYPE(distrib),INTENT(IN) :: old_dist
    REAL, DIMENSION(old_dist%jjb_u:,:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_u:,:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,3)*size(FieldS,4)
    
    CALL  Register_SwapField_gen_u(FieldS,FieldR,ll,old_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField3d_u2d_bis 







  SUBROUTINE Register_SwapField1d_v(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    REAL, DIMENSION(current_dist%ijb_v:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_v:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down

    CALL  Register_SwapField_gen_v(FieldS,FieldR,1,current_dist,new_dist,halo_up,halo_down,a_request)
        
  END SUBROUTINE  Register_SwapField1d_v 

  SUBROUTINE Register_SwapField1d_v_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    TYPE(distrib),INTENT(IN) :: old_dist
    REAL, DIMENSION(old_dist%ijb_v:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_v:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down

    CALL  Register_SwapField_gen_v(FieldS,FieldR,1,old_dist,new_dist,halo_up,halo_down,a_request)
        
  END SUBROUTINE  Register_SwapField1d_v_bis 


  SUBROUTINE Register_SwapField2d_v1d(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
   
    TYPE(distrib),INTENT(IN)          :: new_dist
    REAL, DIMENSION(current_dist%ijb_v:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_v:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,2)
    
    CALL  Register_SwapField_gen_v(FieldS,FieldR,ll,current_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField2d_v1d
  
  SUBROUTINE Register_SwapField2d_v1d_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
   
    TYPE(distrib),INTENT(IN)          :: new_dist
    TYPE(distrib),INTENT(IN)          :: old_dist
    REAL, DIMENSION(old_dist%ijb_v:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_v:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,2)
    
    CALL  Register_SwapField_gen_v(FieldS,FieldR,ll,old_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField2d_v1d_bis
  
   

  SUBROUTINE Register_SwapField3d_v(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    REAL, DIMENSION(current_dist%ijb_v:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_v:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,2)*size(FieldS,3)
    
    CALL  Register_SwapField_gen_v(FieldS,FieldR,ll,current_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField3d_v 

  SUBROUTINE Register_SwapField3d_v_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    TYPE(distrib),INTENT(IN) :: old_dist
    REAL, DIMENSION(old_dist%ijb_v:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%ijb_v:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,2)*size(FieldS,3)
    
    CALL  Register_SwapField_gen_v(FieldS,FieldR,ll,old_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField3d_v_bis 




  SUBROUTINE Register_SwapField1d_v2d(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist !LF
    REAL, DIMENSION(current_dist%jjb_v:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_v:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down

    CALL  Register_SwapField_gen_v(FieldS,FieldR,1,current_dist,new_dist,halo_up,halo_down,a_request)
        
  END SUBROUTINE  Register_SwapField1d_v2d

  SUBROUTINE Register_SwapField1d_v2d_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist !LF
    TYPE(distrib),INTENT(IN) :: old_dist
    REAL, DIMENSION(old_dist%jjb_v:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_v:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down

    CALL  Register_SwapField_gen_v(FieldS,FieldR,1,old_dist,new_dist,halo_up,halo_down,a_request)
        
  END SUBROUTINE  Register_SwapField1d_v2d_bis


  SUBROUTINE Register_SwapField2d_v2d(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    REAL, DIMENSION(current_dist%jjb_v:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_v:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,3)
    
    CALL  Register_SwapField_gen_v(FieldS,FieldR,ll,current_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField2d_v2d
   
  SUBROUTINE Register_SwapField2d_v2d_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    TYPE(distrib),INTENT(IN) :: old_dist
    REAL, DIMENSION(old_dist%jjb_v:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_v:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,3)
    
    CALL  Register_SwapField_gen_v(FieldS,FieldR,ll,old_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField2d_v2d_bis
   

  SUBROUTINE Register_SwapField3d_v2d(FieldS,FieldR,new_dist,a_request,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    REAL, DIMENSION(current_dist%jjb_v:,:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_v:,:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,3)*size(FieldS,4)
    
    CALL  Register_SwapField_gen_v(FieldS,FieldR,ll,current_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField3d_v2d 
  
  SUBROUTINE Register_SwapField3d_v2d_bis(FieldS,FieldR,new_dist,a_request,old_dist,up,down)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
    
    TYPE(distrib),INTENT(IN)          :: new_dist
    TYPE(distrib),INTENT(IN) :: old_dist
    REAL, DIMENSION(old_dist%jjb_v:,:,:,:),INTENT(IN)     :: FieldS
    REAL, DIMENSION(new_dist%jjb_v:,:,:,:),INTENT(OUT)    :: FieldR
    INTEGER,OPTIONAL,INTENT(IN)       :: up
    INTEGER,OPTIONAL,INTENT(IN)       :: down      
    TYPE(request),INTENT(INOUT)         :: a_request

    INTEGER                           :: halo_up
    INTEGER                           :: halo_down
    INTEGER                           :: ll
        
    
    halo_up=0
    halo_down=0
    IF (PRESENT(up))   halo_up=up
    IF (PRESENT(down)) halo_down=down
    
    ll=size(FieldS,3)*size(FieldS,4)
    
    CALL  Register_SwapField_gen_v(FieldS,FieldR,ll,old_dist,new_dist,halo_up,halo_down,a_request)
    
  END SUBROUTINE  Register_SwapField3d_v2d_bis 
  
  

  SUBROUTINE Register_SwapField_gen_u(FieldS,FieldR,ll,old_dist,new_dist,Up,Down,a_request)
  USE parallel_lmdz
  USE dimensions_mod
      IMPLICIT NONE
   
    INTEGER :: ll,Up,Down
    TYPE(distrib)  :: old_dist
    TYPE(distrib)  :: new_dist
    REAL, DIMENSION(old_dist%ijb_u:old_dist%ije_u,ll) :: FieldS
    REAL, DIMENSION(new_dist%ijb_u:new_dist%ije_u,ll) :: FieldR
    TYPE(request) :: a_request
    INTEGER,DIMENSION(0:MPI_Size-1) :: jj_Nb_New   
    INTEGER,DIMENSION(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    INTEGER ::i,l,jje,jjb,ijb,ije
    
    DO i=0,MPI_Size-1
      jj_begin_New(i)=max(1,new_dist%jj_begin_para(i)-Up)
      jj_end_New(i)=min(jjp1,new_dist%jj_end_para(i)+Down)
    ENDDO
   
    DO i=0,MPI_Size-1
      IF (i /= MPI_Rank) THEN
        jjb=max(jj_begin_new(i),old_dist%jj_begin)
        jje=min(jj_end_new(i),old_dist%jj_end)
        
        IF (jje >= jjb) THEN
          CALL Register_SendField(FieldS,old_dist%ijnb_u,ll,jjb-old_dist%jjb_u+1,jje-jjb+1,i,a_request) 
        ENDIF
        
        jjb=max(jj_begin_new(MPI_Rank),old_dist%jj_begin_Para(i))
        jje=min(jj_end_new(MPI_Rank),old_dist%jj_end_Para(i))
        
        IF (jje >= jjb) THEN
          CALL Register_RecvField(FieldR,new_dist%ijnb_u,ll,jjb-new_dist%jjb_u+1,jje-jjb+1,i,a_request) 
        ENDIF
      ELSE
        jjb=max(jj_begin_new(i),old_dist%jj_begin)
        jje=min(jj_end_new(i),old_dist%jj_end)
        ijb=(jjb-1)*iip1+1
        ije=jje*iip1
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)           
        DO l=1,ll
          FieldR(ijb:ije,l)=FieldS(ijb:ije,l)              
        ENDDO
!$OMP END DO NOWAIT
      ENDIF
    ENDDO
    
  END SUBROUTINE Register_SwapField_gen_u



  SUBROUTINE Register_SwapField_gen_v(FieldS,FieldR,ll,old_dist,new_dist,Up,Down,a_request)
  USE parallel_lmdz
  USE dimensions_mod
    IMPLICIT NONE
    
    INTEGER :: ll,Up,Down
    TYPE(distrib)  :: old_dist
    TYPE(distrib)  :: new_dist
    REAL, DIMENSION(old_dist%ijb_v:old_dist%ije_v,ll) :: FieldS
    REAL, DIMENSION(new_dist%ijb_v:new_dist%ije_v,ll) :: FieldR
    TYPE(request) :: a_request
    INTEGER,DIMENSION(0:MPI_Size-1) :: jj_Nb_New   
    INTEGER,DIMENSION(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    INTEGER ::i,l,jje,jjb,ijb,ije
    
    DO i=0,MPI_Size-1
      jj_begin_New(i)=max(1,new_dist%jj_begin_para(i)-Up)
      jj_end_New(i)=min(jjp1,new_dist%jj_end_para(i)+Down)
    ENDDO
   
    DO i=0,MPI_Size-1
      IF (i /= MPI_Rank) THEN
        jjb=max(jj_begin_new(i),old_dist%jj_begin)
        jje=min(jj_end_new(i),old_dist%jj_end)

        IF (jje==jjp1) jje=jjm        

        IF (jje >= jjb) THEN
          CALL Register_SendField(FieldS,old_dist%ijnb_v,ll,jjb-old_dist%jjb_v+1,jje-jjb+1,i,a_request) 
        ENDIF
        
        jjb=max(jj_begin_new(MPI_Rank),old_dist%jj_begin_Para(i))
        jje=min(jj_end_new(MPI_Rank),old_dist%jj_end_Para(i))

        IF (jje==jjp1) jje=jjm
        
        IF (jje >= jjb) THEN
          CALL Register_RecvField(FieldR,new_dist%ijnb_v,ll,jjb-new_dist%jjb_v+1,jje-jjb+1,i,a_request) 
        ENDIF
      ELSE
        jjb=max(jj_begin_new(i),old_dist%jj_begin)
        jje=min(jj_end_new(i),old_dist%jj_end)
        IF (jje==jjp1) jje=jjm
        ijb=(jjb-1)*iip1+1
        ije=jje*iip1
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)           
        DO l=1,ll
          FieldR(ijb:ije,l)=FieldS(ijb:ije,l)
        ENDDO              
!$OMP END DO NOWAIT
      ENDIF
    ENDDO
    
  END SUBROUTINE Register_SwapField_gen_v


 

  
  subroutine Register_Hallo(Field,ij,ll,RUp,Rdown,SUp,SDown,a_request)
  USE dimensions_mod
      implicit none

#ifdef CPP_MPI
    include 'mpif.h'
#endif    
      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: Sup,Sdown,rup,rdown
      type(request) :: a_request
      type(Hallo),pointer :: PtrHallo
      LOGICAL :: SendUp,SendDown
      LOGICAL :: RecvUp,RecvDown
   
 
      SendUp=.TRUE.
      SendDown=.TRUE.
      RecvUp=.TRUE.
      RecvDown=.TRUE.
        
      IF (pole_nord) THEN
        SendUp=.FALSE.
        RecvUp=.FALSE.
      ENDIF
  
      IF (pole_sud) THEN
        SendDown=.FALSE.
        RecvDown=.FALSE.
      ENDIF
      
      if (Sup.eq.0) then
        SendUp=.FALSE.
       endif
      
      if (Sdown.eq.0) then
        SendDown=.FALSE.
      endif

      if (Rup.eq.0) then
        RecvUp=.FALSE.
      endif
      
      if (Rdown.eq.0) then
        RecvDown=.FALSE.
      endif
      
      IF (SendUp) THEN
        call Register_SendField(Field,ij,ll,jj_begin,SUp,MPI_Rank-1,a_request)
      ENDIF
  
      IF (SendDown) THEN
        call Register_SendField(Field,ij,ll,jj_end-SDown+1,SDown,MPI_Rank+1,a_request)
      ENDIF
    
  
      IF (RecvUp) THEN
        call Register_RecvField(Field,ij,ll,jj_begin-Rup,RUp,MPI_Rank-1,a_request)
      ENDIF
  
      IF (RecvDown) THEN
        call Register_RecvField(Field,ij,ll,jj_end+1,RDown,MPI_Rank+1,a_request)
      ENDIF
  
    end subroutine Register_Hallo


  subroutine Register_Hallo_u(Field,ll,RUp,Rdown,SUp,SDown,a_request)
  USE dimensions_mod
      implicit none
#ifdef CPP_MPI
    include 'mpif.h'
#endif    
      INTEGER :: ll
      REAL, dimension(ijb_u:ije_u,ll) :: Field
      INTEGER :: Sup,Sdown,rup,rdown
      type(request) :: a_request
      type(Hallo),pointer :: PtrHallo
      LOGICAL :: SendUp,SendDown
      LOGICAL :: RecvUp,RecvDown
   
 
      SendUp=.TRUE.
      SendDown=.TRUE.
      RecvUp=.TRUE.
      RecvDown=.TRUE.
        
      IF (pole_nord) THEN
        SendUp=.FALSE.
        RecvUp=.FALSE.
      ENDIF
  
      IF (pole_sud) THEN
        SendDown=.FALSE.
        RecvDown=.FALSE.
      ENDIF
      
      if (Sup.eq.0) then
        SendUp=.FALSE.
       endif
      
      if (Sdown.eq.0) then
        SendDown=.FALSE.
      endif

      if (Rup.eq.0) then
        RecvUp=.FALSE.
      endif
      
      if (Rdown.eq.0) then
        RecvDown=.FALSE.
      endif
      
      IF (SendUp) THEN
        call Register_SendField(Field,ijnb_u,ll,jj_begin-jjb_u+1,SUp,MPI_Rank-1,a_request)
      ENDIF
  
      IF (SendDown) THEN
        call Register_SendField(Field,ijnb_u,ll,jj_end-SDown+1-jjb_u+1,SDown,MPI_Rank+1,a_request)
      ENDIF
    
  
      IF (RecvUp) THEN
        call Register_RecvField(Field,ijnb_u,ll,jj_begin-Rup-jjb_u+1,RUp,MPI_Rank-1,a_request)
      ENDIF
  
      IF (RecvDown) THEN
        call Register_RecvField(Field,ijnb_u,ll,jj_end+1-jjb_u+1,RDown,MPI_Rank+1,a_request)
      ENDIF
  
    end subroutine Register_Hallo_u

  subroutine Register_Hallo_v(Field,ll,RUp,Rdown,SUp,SDown,a_request)
  USE dimensions_mod
      implicit none
#ifdef CPP_MPI
    include 'mpif.h'
#endif    
      INTEGER :: ll
      REAL, dimension(ijb_v:ije_v,ll) :: Field
      INTEGER :: Sup,Sdown,rup,rdown
      type(request) :: a_request
      type(Hallo),pointer :: PtrHallo
      LOGICAL :: SendUp,SendDown
      LOGICAL :: RecvUp,RecvDown
   
 
      SendUp=.TRUE.
      SendDown=.TRUE.
      RecvUp=.TRUE.
      RecvDown=.TRUE.
        
      IF (pole_nord) THEN
        SendUp=.FALSE.
        RecvUp=.FALSE.
      ENDIF
  
      IF (pole_sud) THEN
        SendDown=.FALSE.
        RecvDown=.FALSE.
      ENDIF
      
      if (Sup.eq.0) then
        SendUp=.FALSE.
       endif
      
      if (Sdown.eq.0) then
        SendDown=.FALSE.
      endif

      if (Rup.eq.0) then
        RecvUp=.FALSE.
      endif
      
      if (Rdown.eq.0) then
        RecvDown=.FALSE.
      endif
      
      IF (SendUp) THEN
        call Register_SendField(Field,ijnb_v,ll,jj_begin-jjb_v+1,SUp,MPI_Rank-1,a_request)
      ENDIF
  
      IF (SendDown) THEN
        call Register_SendField(Field,ijnb_v,ll,jj_end-SDown+1-jjb_v+1,SDown,MPI_Rank+1,a_request)
      ENDIF
    
  
      IF (RecvUp) THEN
        call Register_RecvField(Field,ijnb_v,ll,jj_begin-Rup-jjb_v+1,RUp,MPI_Rank-1,a_request)
      ENDIF
  
      IF (RecvDown) THEN
        call Register_RecvField(Field,ijnb_v,ll,jj_end+1-jjb_v+1,RDown,MPI_Rank+1,a_request)
      ENDIF
  
    end subroutine Register_Hallo_v
    
    subroutine SendRequest(a_Request)
    USE dimensions_mod
      implicit none

#ifdef CPP_MPI
      include 'mpif.h'
#endif

      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer :: SizeBuffer
      integer :: i,rank,l,ij,Pos,ierr
      integer :: offset
      real,dimension(:,:),pointer :: Field
      integer :: Nb
       
      do rank=0,MPI_SIZE-1
      
        Req=>a_Request%RequestSend(rank)
        
        SizeBuffer=0
        do i=1,Req%NbRequest
          PtrHallo=>Req%Hallo(i)
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
          DO l=1,PtrHallo%NbLevel
            SizeBuffer=SizeBuffer+PtrHallo%size*iip1
          ENDDO
!$OMP ENDDO NOWAIT          
        enddo
      
         Req%BufferSize=SizeBuffer
         if (Req%NbRequest>0) then
       
          call allocate_buffer(SizeBuffer,Req%Index,Req%pos)

          Pos=Req%Pos
          do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)
            offset=(PtrHallo%offset-1)*iip1+1
            Nb=iip1*PtrHallo%size-1
            Field=>PtrHallo%Field

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)           
            do l=1,PtrHallo%NbLevel
!cdir NODEP
              do ij=0,Nb
	        Buffer(Pos+ij)=Field(Offset+ij,l)
	      enddo
              
              Pos=Pos+Nb+1
            enddo
!$OMP END DO NOWAIT            
          enddo
    
         if (SizeBuffer>0) then
!$OMP CRITICAL (MPI)
         
#ifdef CPP_MPI
         call MPI_ISEND(Buffer(req%Pos),SizeBuffer,MPI_REAL_LMDZ,rank,a_request%tag+1000*omp_rank,     &
                         COMM_LMDZ,Req%MSG_Request,ierr)
#endif
         IF (.NOT.using_mpi) THEN
           PRINT *,'Erreur, echange MPI en mode sequentiel !!!'
           STOP
         ENDIF
!         PRINT *,"-------------------------------------------------------------------"
!         PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->"
!         PRINT *,"Requete envoye au proc :",rank,"tag :",a_request%tag+1000*omp_rank
!         PRINT *,"Taille du message :",SizeBuffer,"requete no :",Req%MSG_Request
!         PRINT *,"-------------------------------------------------------------------"
!$OMP END CRITICAL (MPI)
        endif
       endif
    enddo
   
           
      do rank=0,MPI_SIZE-1
         
          Req=>a_Request%RequestRecv(rank)
          SizeBuffer=0
          
	  do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
            DO l=1,PtrHallo%NbLevel
              SizeBuffer=SizeBuffer+PtrHallo%size*iip1
            ENDDO
!$OMP ENDDO NOWAIT          
          enddo
          
          Req%BufferSize=SizeBuffer
          
          if (Req%NbRequest>0) then
          call allocate_buffer(SizeBuffer,Req%Index,Req%Pos)
   
          if (SizeBuffer>0) then

!$OMP CRITICAL (MPI)

#ifdef CPP_MPI
             call MPI_IRECV(Buffer(Req%Pos),SizeBuffer,MPI_REAL_LMDZ,rank,a_request%tag+1000*omp_rank,     &
                           COMM_LMDZ,Req%MSG_Request,ierr)
#endif             
             IF (.NOT.using_mpi) THEN
               PRINT *,'Erreur, echange MPI en mode sequentiel !!!'
               STOP
             ENDIF

!         PRINT *,"-------------------------------------------------------------------"
!         PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->"
!         PRINT *,"Requete en attente du proc :",rank,"tag :",a_request%tag+1000*omp_rank
!         PRINT *,"Taille du message :",SizeBuffer,"requete no :",Req%MSG_Request
!         PRINT *,"-------------------------------------------------------------------"

!$OMP END CRITICAL (MPI)
          endif
        endif
      
      enddo
                        
   end subroutine SendRequest 
   
   subroutine WaitRequest(a_Request)
   USE dimensions_mod
   implicit none
   
#ifdef CPP_MPI
      include 'mpif.h'   
#endif
      
      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer, dimension(2*mpi_size) :: TabRequest
#ifdef CPP_MPI
      integer, dimension(MPI_STATUS_SIZE,2*mpi_size) :: TabStatus
#else
      integer, dimension(1,2*mpi_size) :: TabStatus
#endif
      integer :: NbRequest
      integer :: i,rank,pos,ij,l,ierr
      integer :: offset
      integer :: Nb
      
      
      NbRequest=0
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0 .AND. Req%BufferSize > 0) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0 .AND. Req%BufferSize > 0 ) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
     
      if (NbRequest>0) then
!$OMP CRITICAL (MPI)
!        PRINT *,"-------------------------------------------------------------------"
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"en attente"
!        PRINT *,"No des requetes :",TabRequest(1:NbRequest)
#ifdef CPP_MPI
        call MPI_WAITALL(NbRequest,TabRequest,TabStatus,ierr)
#endif
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"complete"
!        PRINT *,"-------------------------------------------------------------------"
!$OMP END CRITICAL (MPI)
      endif
      do rank=0,MPI_Size-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          Pos=Req%Pos
          do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)
            offset=(PtrHallo%offset-1)*iip1+1
	    Nb=iip1*PtrHallo%size-1

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
	    do l=1,PtrHallo%NbLevel
!cdir NODEP
              do ij=0,Nb
	        PtrHallo%Field(offset+ij,l)=Buffer(Pos+ij)
	      enddo

              Pos=Pos+Nb+1
	    enddo
!$OMP ENDDO NOWAIT	    
          enddo
        endif
      enddo
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
              
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
     
      a_request%tag=1
    end subroutine WaitRequest
     
   subroutine WaitSendRequest(a_Request)
   USE dimensions_mod
   implicit none
   
#ifdef CPP_MPI
      include 'mpif.h'   
#endif      
      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer, dimension(mpi_size) :: TabRequest
#ifdef CPP_MPI
      integer, dimension(MPI_STATUS_SIZE,mpi_size) :: TabStatus
#else
      integer, dimension(1,mpi_size) :: TabStatus
#endif
      integer :: NbRequest
      integer :: i,rank,pos,ij,l,ierr
      integer :: offset
      
      
      NbRequest=0
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
      

      if (NbRequest>0 .AND. Req%BufferSize > 0 ) THEN 
!$OMP CRITICAL (MPI)     
!        PRINT *,"-------------------------------------------------------------------"
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"en attente"
!        PRINT *,"No des requetes :",TabRequest(1:NbRequest)
#ifdef CPP_MPI
        call MPI_WAITALL(NbRequest,TabRequest,TabStatus,ierr)
#endif
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"complete"
!        PRINT *,"-------------------------------------------------------------------"

!$OMP END CRITICAL (MPI)
      endif      
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
              
      a_request%tag=1
    end subroutine WaitSendRequest
    
   subroutine WaitRecvRequest(a_Request)
   USE dimensions_mod
   implicit none
   
#ifdef CPP_MPI
      include 'mpif.h'   
#endif
      
      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer, dimension(mpi_size) :: TabRequest
#ifdef CPP_MPI
      integer, dimension(MPI_STATUS_SIZE,mpi_size) :: TabStatus
#else
      integer, dimension(1,mpi_size) :: TabStatus
#endif
      integer :: NbRequest
      integer :: i,rank,pos,ij,l,ierr
      integer :: offset,Nb
      
      
      NbRequest=0
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0 .AND. Req%BufferSize > 0 ) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
     
      
      if (NbRequest>0) then
!$OMP CRITICAL (MPI)     
!        PRINT *,"-------------------------------------------------------------------"
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"en attente"
!        PRINT *,"No des requetes :",TabRequest(1:NbRequest)
#ifdef CPP_MPI
        call MPI_WAITALL(NbRequest,TabRequest,TabStatus,ierr)
#endif
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"complete"
!        PRINT *,"-------------------------------------------------------------------"
!$OMP END CRITICAL (MPI)     
      endif
      
      do rank=0,MPI_Size-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          Pos=Req%Pos
          do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)
            offset=(PtrHallo%offset-1)*iip1+1
	    Nb=iip1*PtrHallo%size-1
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
	    do l=1,PtrHallo%NbLevel
!cdir NODEP
              do ij=0,Nb
	        PtrHallo%Field(offset+ij,l)=Buffer(Pos+ij)
	      enddo
                 Pos=Pos+Nb+1
            enddo
!$OMP END DO NOWAIT
          enddo
        endif
      enddo
      
           
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
     
      a_request%tag=1
    end subroutine WaitRecvRequest
    
    
    
    subroutine CopyField(FieldS,FieldR,ij,ll,jj_Nb_New)
    USE dimensions_mod
  
      implicit none
    
    INTEGER :: ij,ll,l
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    integer ::i,jje,jjb,ijb,ije
    
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo
    
    jjb=max(jj_begin,jj_begin_new(MPI_Rank))
    jje=min(jj_end,jj_end_new(MPI_Rank))
    if (ij==ip1jm) jje=min(jje,jjm)

    if (jje >= jjb) then
      ijb=(jjb-1)*iip1+1
      ije=jje*iip1

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      do l=1,ll
        FieldR(ijb:ije,l)=FieldS(ijb:ije,l)
      enddo
!$OMP ENDDO NOWAIT
    endif


  end subroutine CopyField    

  subroutine CopyFieldHallo(FieldS,FieldR,ij,ll,jj_Nb_New,Up,Down)
  USE dimensions_mod
  
      implicit none
    
    INTEGER :: ij,ll,Up,Down
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New

    integer ::i,jje,jjb,ijb,ije,l

     
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo

        
    jjb=max(jj_begin,jj_begin_new(MPI_Rank)-Up)
    jje=min(jj_end,jj_end_new(MPI_Rank)+Down)
    if (ij==ip1jm) jje=min(jje,jjm)
    
    
    if (jje >= jjb) then
      ijb=(jjb-1)*iip1+1
      ije=jje*iip1

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      do l=1,ll
        FieldR(ijb:ije,l)=FieldS(ijb:ije,l)
      enddo
!$OMP ENDDO NOWAIT

    endif
   end subroutine CopyFieldHallo        

   subroutine Gather_field_u(field_loc,field_glo,ll)
   USE dimensions_mod
   implicit none
     integer :: ll
     real :: field_loc(ijb_u:ije_u,ll)
     real :: field_glo(ip1jmp1,ll)
     type(request) :: request_gather
     integer       :: l


!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
     DO l=1,ll
       field_glo(ij_begin:ij_end,l)=field_loc(ij_begin:ij_end,l)
     ENDDO
     
     call register_SwapField(field_glo,field_glo,ip1jmp1,ll,distrib_gather%jj_nb_para,request_gather)
     call SendRequest(request_gather)
!$OMP BARRIER
     call WaitRequest(request_gather)       
!$OMP BARRIER

    end subroutine Gather_field_u
        
   subroutine Gather_field_v(field_loc,field_glo,ll)
   USE dimensions_mod
   implicit none
     integer :: ll
     real :: field_loc(ijb_v:ije_v,ll)
     real :: field_glo(ip1jm,ll)
     type(request) :: request_gather
     integer :: ijb,ije
     integer       :: l
     
   
     ijb=ij_begin
     ije=ij_end
     if (pole_sud) ije=ij_end-iip1
        
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
     DO l=1,ll
       field_glo(ijb:ije,l)=field_loc(ijb:ije,l)
     ENDDO
     
     call register_SwapField(field_glo,field_glo,ip1jm,ll,distrib_gather%jj_nb_para,request_gather)
     call SendRequest(request_gather)
!$OMP BARRIER
     call WaitRequest(request_gather)       
!$OMP BARRIER

    end subroutine Gather_field_v
     
   subroutine Scatter_field_u(field_glo,field_loc,ll)
   USE dimensions_mod
   implicit none
     integer :: ll
     real :: field_glo(ip1jmp1,ll)
     real :: field_loc(ijb_u:ije_u,ll)
     type(request) :: request_gather
     TYPE(distrib) :: distrib_swap
     integer       :: l
     
!$OMP BARRIER
!$OMP MASTER     
     call get_current_distrib(distrib_swap)
     call set_Distrib(distrib_gather)
!$OMP END MASTER
!$OMP BARRIER
 
     call register_SwapField(field_glo,field_glo,ip1jmp1,ll,distrib_swap%jj_nb_para,request_gather)
     call SendRequest(request_gather)
!$OMP BARRIER
     call WaitRequest(request_gather)       
!$OMP BARRIER
!$OMP MASTER     
     call set_Distrib(distrib_swap)
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
       DO l=1,ll
         field_loc(ij_begin:ij_end,l)=field_glo(ij_begin:ij_end,l)
       ENDDO

    end subroutine Scatter_field_u

   subroutine Scatter_field_v(field_glo,field_loc,ll)
   USE dimensions_mod
   implicit none
     integer :: ll
     real :: field_glo(ip1jmp1,ll)
     real :: field_loc(ijb_v:ije_v,ll)
     type(request) :: request_gather
     TYPE(distrib) :: distrib_swap
     integer       :: ijb,ije,l
     

!$OMP BARRIER
!$OMP MASTER     
     call get_current_distrib(distrib_swap)
     call set_Distrib(distrib_gather)
!$OMP END MASTER
!$OMP BARRIER
     call register_SwapField(field_glo,field_glo,ip1jm,ll,distrib_swap%jj_nb_para,request_gather)
     call SendRequest(request_gather)
!$OMP BARRIER
     call WaitRequest(request_gather)       
!$OMP BARRIER
!$OMP MASTER
     call set_Distrib(distrib_swap)
!$OMP END MASTER
!$OMP BARRIER
     ijb=ij_begin
     ije=ij_end
     if (pole_sud) ije=ij_end-iip1
     
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
       DO l=1,ll
         field_loc(ijb:ije,l)=field_glo(ijb:ije,l)
       ENDDO

    end subroutine Scatter_field_v
              
end module mod_Hallo 
   
