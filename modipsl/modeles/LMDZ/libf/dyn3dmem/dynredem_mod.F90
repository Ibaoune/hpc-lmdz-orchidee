MODULE dynredem_mod

  USE dimensions_mod
  USE parallel_lmdz
  USE mod_hallo
  USE netcdf
  PRIVATE
  PUBLIC :: dynredem_write_u, dynredem_write_v, dynredem_read_u, err
  PUBLIC :: cre_var, get_var1, put_var, fil, modname, msg
  CHARACTER(LEN=256), SAVE :: fil, modname
  INTEGER,            SAVE :: nvarid


CONTAINS


!===============================================================================
!
SUBROUTINE dynredem_write_u(ncid,id,var,ll)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,          INTENT(IN) :: ncid
  CHARACTER(LEN=*), INTENT(IN) :: id
  REAL,             INTENT(IN) :: var(ijb_u:ije_u,ll)
  INTEGER,          INTENT(IN) :: ll
!===============================================================================
! Local variables:
  REAL, ALLOCATABLE, SAVE :: var_tmp(:,:), var_glo(:)
  INTEGER :: start(4), count(4), l, ierr
!===============================================================================
  start(:)=[1,1,1,1]; count(:)=[iip1,jjp1,1,1]

!$OMP MASTER
  IF(mpi_rank==0) CALL err(NF90_INQ_VARID(ncid,id,nvarid),"inq",id)
!$OMP END MASTER

!$OMP MASTER
  ALLOCATE(var_tmp(ijb_u:ije_u,ll),var_glo(ip1jmp1))
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,ll; var_tmp(:,l)=var(:,l); END DO
  DO l=1,ll
    CALL gather_field_u(var_tmp(:,l),var_glo,1)
    IF(mpi_rank==0) THEN
    !$OMP MASTER
      start(3)=l
      CALL err(NF90_PUT_VAR(ncid,nvarid,var_glo,start,count),"put",id)
    !$OMP END MASTER
    END IF
  END DO
!$OMP BARRIER
!$OMP MASTER
  DEALLOCATE(var_glo,var_tmp)
!$OMP END MASTER
!$OMP BARRIER
  
END SUBROUTINE dynredem_write_u
!
!===============================================================================


!===============================================================================
!
SUBROUTINE dynredem_write_v(ncid,id,var,ll)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,          INTENT(IN) :: ncid
  CHARACTER(LEN=*), INTENT(IN) :: id
  REAL,             INTENT(IN) :: var(ijb_v:ije_v,ll)
  INTEGER,          INTENT(IN) :: ll
!===============================================================================
! Local variables:
  REAL, ALLOCATABLE, SAVE :: var_tmp(:,:), var_glo(:)
  INTEGER :: start(4), count(4), l, ierr
!===============================================================================
  start(:)=[1,1,1,1]; count(:)=[iip1,jjm,1,1]

!$OMP MASTER
  IF(mpi_rank==0) CALL err(NF90_INQ_VARID(ncid,id,nvarid),"inq",id)
!$OMP END MASTER

!$OMP MASTER
  ALLOCATE(var_tmp(ijb_v:ije_v,ll),var_glo(ip1jm))
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,ll; var_tmp(:,l)=var(:,l); END DO
  DO l=1,ll
    CALL gather_field_v(var_tmp(:,l),var_glo,1)
    IF(mpi_rank==0) THEN
    !$OMP MASTER
      start(3)=l
      CALL err(NF90_PUT_VAR(ncid,nvarid,var_glo,start,count),"put",id)
    !$OMP END MASTER
    END IF
  END DO
!$OMP BARRIER
!$OMP MASTER
  DEALLOCATE(var_glo,var_tmp)
!$OMP END MASTER
!$OMP BARRIER
  
END SUBROUTINE dynredem_write_v
!
!===============================================================================


!===============================================================================
!
SUBROUTINE dynredem_read_u(ncid,id,var,ll)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,          INTENT(IN)  :: ncid
  CHARACTER(LEN=*), INTENT(IN)  :: id
  REAL,             INTENT(OUT) :: var(ijb_u:ije_u,ll)
  INTEGER,          INTENT(IN)  :: ll
!===============================================================================
! Local variables:
  REAL, ALLOCATABLE, SAVE :: var_tmp(:,:), var_glo(:)
  INTEGER :: start(4), count(4), l, ierr
!===============================================================================
  start(:)=[1,1,1,1]; count(:)=[iip1,jjp1,1,1]

!$OMP MASTER
  IF(mpi_rank==0) CALL err(NF90_INQ_VARID(ncid,id,nvarid),'inq',id)
!$OMP END MASTER

!$OMP MASTER
  ALLOCATE(var_tmp(ijb_u:ije_u,ll),var_glo(ip1jmp1))
!$OMP END MASTER
!$OMP BARRIER

  DO l=1,ll
    IF(mpi_rank==0) THEN
    !$OMP MASTER
      start(3)=l
      CALL err(NF90_GET_VAR(ncid,nvarid,var_glo,start,count),"get",id)
    !$OMP END MASTER
    END IF
    CALL scatter_field_u(var_glo,var_tmp(:,l),1)
  END DO

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,ll; var(:,l)=var_tmp(:,l); END DO
    
!$OMP BARRIER
!$OMP MASTER
  DEALLOCATE(var_glo,var_tmp)
!$OMP END MASTER
!$OMP BARRIER
  
END SUBROUTINE dynredem_read_u    
!
!===============================================================================


!===============================================================================
!
SUBROUTINE cre_var(ncid,var,title,did,units)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,                    INTENT(IN) :: ncid
  CHARACTER(LEN=*),           INTENT(IN) :: var, title
  INTEGER,                    INTENT(IN) :: did(:)
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: units
!===============================================================================
#ifdef NC_DOUBLE
  CALL err(NF90_DEF_VAR(ncid,var,NF90_DOUBLE,did,nvarid),"inq",var)
#else
  CALL err(NF90_DEF_VAR(ncid,var,NF90_FLOAT ,did,nvarid),"inq",var)
#endif
  IF(title/="")      CALL err(NF90_PUT_ATT(ncid,nvarid,"title",title),var)
  IF(PRESENT(units)) CALL err(NF90_PUT_ATT(ncid,nvarid,"units",units),var)

END SUBROUTINE cre_var
!
!===============================================================================


!===============================================================================
!
SUBROUTINE put_var(ncid,var,title,did,v,units)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,                    INTENT(IN) :: ncid
  CHARACTER(LEN=*),           INTENT(IN) :: var, title
  INTEGER,                    INTENT(IN) :: did(:)
  REAL,                       INTENT(IN) :: v(:)
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: units
!===============================================================================
  INTEGER :: nd, k, nn(2)
  IF(     PRESENT(units)) CALL cre_var(ncid,var,title,did,units)
  IF(.NOT.PRESENT(units)) CALL cre_var(ncid,var,title,did)
  CALL err(NF90_ENDDEF(ncid))
  nd=SIZE(did)
  DO k=1,nd; CALL err(NF90_INQUIRE_DIMENSION(ncid,did(k),len=nn(k))); END DO
  IF(nd==1) CALL err(NF90_PUT_VAR(ncid,nvarid,RESHAPE(v,nn(1:1))),var)
  IF(nd==2) CALL err(NF90_PUT_VAR(ncid,nvarid,RESHAPE(v,nn(1:2))),var)
  CALL err(NF90_REDEF(ncid))
END SUBROUTINE put_var
!
!===============================================================================


!===============================================================================
!
FUNCTION msg(typ,nam)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  CHARACTER(LEN=256)                     :: msg    !--- STANDARDIZED MESSAGE
  CHARACTER(LEN=*),           INTENT(IN) :: typ    !--- TYPE OF OPERATION
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: nam    !--- FIELD NAME
!===============================================================================
  SELECT CASE(typ)
    CASE('open');  msg="Opening failed for <"//TRIM(fil)//">"
    CASE('close'); msg="Closing failed for <"//TRIM(fil)//">"
    CASE('get');   msg="Reading failed for <"//TRIM(nam)//">"
    CASE('put');   msg="Writting failed for <"//TRIM(nam)//">"
    CASE('inq');   msg="Missing field <"//TRIM(nam)//">"
    CASE('fnd');   msg="Found field <"//TRIM(nam)//">"
  END SELECT
  msg=TRIM(msg)//" in file <"//TRIM(fil)//">"

END FUNCTION msg
!
!===============================================================================


!===============================================================================
!
SUBROUTINE err(ierr,typ,nam)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,                    INTENT(IN) :: ierr   !--- NetCDF ERROR CODE
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: typ    !--- TYPE OF OPERATION
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: nam    !--- FIELD NAME
!===============================================================================
  IF(ierr==NF90_NoERR) RETURN
  IF(.NOT.PRESENT(typ)) THEN
    CALL ABORT_gcm(modname,NF90_STRERROR(ierr),ierr)
  ELSE
    CALL ABORT_gcm(modname,msg(typ,nam),ierr)
  END IF

END SUBROUTINE err
!
!===============================================================================

END MODULE dynredem_mod   

    
    
