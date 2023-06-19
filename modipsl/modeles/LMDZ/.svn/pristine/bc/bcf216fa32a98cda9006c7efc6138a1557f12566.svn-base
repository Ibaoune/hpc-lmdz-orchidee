module write_field_loc
implicit none
  
  interface WriteField_u
    module procedure Write_field1d_u,Write_Field2d_u
  end interface WriteField_u

  interface WriteField_v
    module procedure Write_field1d_v,Write_Field2d_v
  end interface WriteField_v
  
  contains
  
  subroutine write_field1D_u(name,Field)
    character(len=*)   :: name
    real, dimension(:) :: Field

    CALL write_field_u_gen(name,Field,1)

  end subroutine write_field1D_u

  subroutine write_field2D_u(name,Field)
    implicit none
      
    character(len=*)   :: name
    real, dimension(:,:) :: Field
    integer :: ll
    
    ll=size(field,2)    
    CALL write_field_u_gen(name,Field,ll)
    
    end subroutine write_field2D_u


   SUBROUTINE write_field_u_gen(name,Field,ll)
    USE parallel_lmdz
    USE write_field
    USE mod_hallo
    implicit none
    include 'dimensions.h'
    include 'paramet.h'
      
    character(len=*)   :: name
    real, dimension(ijb_u:ije_u,ll) :: Field
    real, allocatable,SAVE :: New_Field(:,:,:)
    integer,dimension(0:mpi_size-1) :: jj_nb_master
    type(Request),SAVE :: Request_write
!$OMP THREADPRIVATE(Request_write)
    integer :: ll,i
    
    
    jj_nb_master(:)=0
    jj_nb_master(0)=jjp1
!$OMP BARRIER
!$OMP MASTER
    allocate(New_Field(iip1,jjp1,ll))
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    DO i=1,ll    
      New_Field(:,jj_begin:jj_end,i)=reshape(Field(ij_begin:ij_end,i),(/iip1,jj_nb/))
    ENDDO
!$OMP BARRIER    
    call Register_SwapField(new_field,new_field,ip1jmp1,ll,jj_Nb_master,Request_write)
    call SendRequest(Request_write)
!$OMP BARRIER
    call WaitRequest(Request_write)     
!$OMP BARRIER

!$OMP MASTER
    if (MPI_Rank==0) call WriteField(name,New_Field)
    DEALLOCATE(New_Field)
!$OMP END MASTER        
!$OMP BARRIER
    END SUBROUTINE write_field_u_gen


  subroutine write_field1D_v(name,Field)
    character(len=*)   :: name
    real, dimension(:) :: Field

    CALL write_field_v_gen(name,Field,1)

  end subroutine write_field1D_v

  subroutine write_field2D_v(name,Field)
    implicit none
      
    character(len=*)   :: name
    real, dimension(:,:) :: Field
    integer :: ll
    
    ll=size(field,2)    
    CALL write_field_v_gen(name,Field,ll)
    
    end subroutine write_field2D_v


   SUBROUTINE write_field_v_gen(name,Field,ll)
    USE parallel_lmdz
    USE write_field
    USE mod_hallo
    implicit none
    include 'dimensions.h'
    include 'paramet.h'
      
    character(len=*)   :: name
    real, dimension(ijb_v:ije_v,ll) :: Field
    real, allocatable,SAVE :: New_Field(:,:,:)
    integer,dimension(0:mpi_size-1) :: jj_nb_master
    type(Request),SAVE :: Request_write
!$OMP THREADPRIVATE(Request_write)    
    integer :: ll,i,jje,ije,jjn
    
    
    jj_nb_master(:)=0
    jj_nb_master(0)=jjp1

!$OMP BARRIER
!$OMP MASTER
    allocate(New_Field(iip1,jjm,ll))
!$OMP END MASTER
!$OMP BARRIER

   IF (pole_sud) THEN
     jje=jj_end-1
     ije=ij_end-iip1
     jjn=jj_nb-1
   ELSE
     jje=jj_end
     ije=ij_end
     jjn=jj_nb
   ENDIF
   
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    DO i=1,ll    
      New_Field(:,jj_begin:jje,i)=reshape(Field(ij_begin:ije,i),(/iip1,jjn/))
    ENDDO
!$OMP BARRIER    
    call Register_SwapField(new_field,new_field,ip1jm,ll,jj_Nb_master,Request_write)
    call SendRequest(Request_write)
!$OMP BARRIER
    call WaitRequest(Request_write)     
!$OMP BARRIER

!$OMP MASTER
    if (MPI_Rank==0) call WriteField(name,New_Field)
    DEALLOCATE(New_Field)
!$OMP END MASTER        
!$OMP BARRIER
    END SUBROUTINE write_field_v_gen
    
end module write_field_loc
  
