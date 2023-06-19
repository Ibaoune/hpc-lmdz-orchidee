SUBROUTINE dynetat0_loc(fichnom,vcov,ucov,teta,q,masse,ps,phis,time)
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , L.Fairhead
!-------------------------------------------------------------------------------
! Purpose: Initial state reading.
!-------------------------------------------------------------------------------
  USE parallel_lmdz
  USE infotrac
  USE netcdf, ONLY: NF90_OPEN,  NF90_INQUIRE_DIMENSION, NF90_INQ_VARID,        &
      NF90_NOWRITE, NF90_CLOSE, NF90_INQUIRE_VARIABLE,  NF90_GET_VAR, NF90_NoErr
  USE control_mod, ONLY: planet_type
  USE assert_eq_m, ONLY: assert_eq
  USE comvert_mod, ONLY: pa,preff
  USE comconst_mod, ONLY: cpp, daysec, dtvr, g, im, jm, kappa, lllm, &
                          omeg, rad
  USE logic_mod, ONLY: fxyhypb, ysinus
  USE serre_mod, ONLY: clon, clat, grossismx, grossismy
  USE temps_mod, ONLY: annee_ref,day_ref,itau_dyn, &
                       start_time,day_ini,hour_ini
  USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0
  
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
  include "description.h"
  include "iniprint.h"
!===============================================================================
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: fichnom          !--- FILE NAME
  REAL, INTENT(OUT) ::  vcov(ijb_v:ije_v,llm)      !--- V COVARIANT WIND
  REAL, INTENT(OUT) ::  ucov(ijb_u:ije_u,llm)      !--- U COVARIANT WIND
  REAL, INTENT(OUT) ::  teta(ijb_u:ije_u,llm)      !--- POTENTIAL TEMP.
  REAL, INTENT(OUT) ::     q(ijb_u:ije_u,llm,nqtot)!--- TRACERS
  REAL, INTENT(OUT) :: masse(ijb_u:ije_u,llm)      !--- MASS PER CELL
  REAL, INTENT(OUT) ::    ps(ijb_u:ije_u)          !--- GROUND PRESSURE
  REAL, INTENT(OUT) ::  phis(ijb_u:ije_u)          !--- GEOPOTENTIAL
!===============================================================================
! Local variables:
  CHARACTER(LEN=256) :: msg, var, modname
  INTEGER, PARAMETER :: length=100
  INTEGER :: iq, fID, vID, idecal, ierr
  REAL    :: time, tab_cntrl(length)               !--- RUN PARAMS TABLE
  REAL,             ALLOCATABLE :: vcov_glo(:,:),masse_glo(:,:),   ps_glo(:)
  REAL,             ALLOCATABLE :: ucov_glo(:,:),    q_glo(:,:), phis_glo(:)
  REAL,             ALLOCATABLE :: teta_glo(:,:)
!-------------------------------------------------------------------------------
  modname="dynetat0_loc"

!--- Initial state file opening
  var=fichnom
  CALL err(NF90_OPEN(var,NF90_NOWRITE,fID),"open",var)
  CALL get_var1("controle",tab_cntrl)

!!! AS: idecal is a hack to be able to read planeto starts...
!!!     .... while keeping everything OK for LMDZ EARTH
  IF(planet_type=="generic") THEN
    WRITE(lunout,*)'NOTE NOTE NOTE : Planeto-like start files'
    idecal = 4
    annee_ref  = 2000
  ELSE
    WRITE(lunout,*)'NOTE NOTE NOTE : Earth-like start files'
    idecal = 5
    annee_ref  = tab_cntrl(5)
  END IF
  im         = tab_cntrl(1)
  jm         = tab_cntrl(2)
  lllm       = tab_cntrl(3)
  day_ref    = tab_cntrl(4)
  rad        = tab_cntrl(idecal+1)
  omeg       = tab_cntrl(idecal+2)
  g          = tab_cntrl(idecal+3)
  cpp        = tab_cntrl(idecal+4)
  kappa      = tab_cntrl(idecal+5)
  daysec     = tab_cntrl(idecal+6)
  dtvr       = tab_cntrl(idecal+7)
  etot0      = tab_cntrl(idecal+8)
  ptot0      = tab_cntrl(idecal+9)
  ztot0      = tab_cntrl(idecal+10)
  stot0      = tab_cntrl(idecal+11)
  ang0       = tab_cntrl(idecal+12)
  pa         = tab_cntrl(idecal+13)
  preff      = tab_cntrl(idecal+14)
!
  clon       = tab_cntrl(idecal+15)
  clat       = tab_cntrl(idecal+16)
  grossismx  = tab_cntrl(idecal+17)
  grossismy  = tab_cntrl(idecal+18)
!
  IF ( tab_cntrl(idecal+19)==1. )  THEN
    fxyhypb  = .TRUE.
!   dzoomx   = tab_cntrl(25)
!   dzoomy   = tab_cntrl(26)
!   taux     = tab_cntrl(28)
!   tauy     = tab_cntrl(29)
  ELSE
    fxyhypb = .FALSE.
    ysinus  = tab_cntrl(idecal+22)==1.
  END IF

  day_ini    = tab_cntrl(30)
  itau_dyn   = tab_cntrl(31)
  start_time = tab_cntrl(32)

!-------------------------------------------------------------------------------
  WRITE(lunout,*)TRIM(modname)//': rad,omeg,g,cpp,kappa',rad,omeg,g,cpp,kappa
  CALL check_dim(im,iim,'im','im')
  CALL check_dim(jm,jjm,'jm','jm')
  CALL check_dim(lllm,llm,'lm','lllm')
  CALL get_var1("rlonu",rlonu)
  CALL get_var1("rlatu",rlatu)
  CALL get_var1("rlonv",rlonv)
  CALL get_var1("rlatv",rlatv)
  CALL get_var1("cu"  ,cu)
  CALL get_var1("cv"  ,cv)
  CALL get_var1("aire",aire)

  var="temps"
  IF(NF90_INQ_VARID(fID,var,vID)/=NF90_NoErr) THEN
    WRITE(lunout,*)TRIM(modname)//": missing field <temps>"
    WRITE(lunout,*)TRIM(modname)//": trying with <Time>"; var="Time"
    CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  END IF
  CALL err(NF90_GET_VAR(fID,vID,time),"get",var)

  ALLOCATE(phis_glo(ip1jmp1))
  CALL get_var1("phisinit",phis_glo)
  phis (ijb_u:ije_u)  =phis_glo(ijb_u:ije_u);    DEALLOCATE(phis_glo)

  ALLOCATE(ucov_glo(ip1jmp1,llm))
  CALL get_var2("ucov",ucov_glo)
  ucov (ijb_u:ije_u,:)=ucov_glo(ijb_u:ije_u,:);  DEALLOCATE(ucov_glo)

  ALLOCATE(vcov_glo(ip1jm,llm))
  CALL get_var2("vcov",vcov_glo)
  vcov (ijb_v:ije_v,:)=vcov_glo(ijb_v:ije_v,:);  DEALLOCATE(vcov_glo)

  ALLOCATE(teta_glo(ip1jmp1,llm))
  CALL get_var2("teta",teta_glo)
  teta (ijb_u:ije_u,:)=teta_glo(ijb_u:ije_u,:);  DEALLOCATE(teta_glo)

  ALLOCATE(masse_glo(ip1jmp1,llm))
  CALL get_var2("masse",masse_glo)
  masse(ijb_u:ije_u,:)=masse_glo(ijb_u:ije_u,:); DEALLOCATE(masse_glo)
  
  ALLOCATE(ps_glo(ip1jmp1))
  CALL get_var1("ps",ps_glo)
  ps   (ijb_u:ije_u)  =   ps_glo(ijb_u:ije_u);   DEALLOCATE(ps_glo)

!--- Tracers
  ALLOCATE(q_glo(ip1jmp1,llm))
  DO iq=1,nqtot
    var=tname(iq)
#ifdef INCA
    IF (var .eq. "O3" ) THEN 
       IF(NF90_INQ_VARID(fID,var,vID) == NF90_NoErr) THEN
          CALL get_var2(var,q_glo); q(ijb_u:ije_u,:,iq)=q_glo(ijb_u:ije_u,:); CYCLE
       ELSE
          WRITE(lunout,*) 'Tracer O3 is missing - it is initialized to OX'
          IF(NF90_INQ_VARID(fID,"OX",vID) == NF90_NoErr) THEN
             CALL get_var2("OX",q_glo); q(ijb_u:ije_u,:,iq)=q_glo(ijb_u:ije_u,:); CYCLE
          ENDIF
       ENDIF
    ENDIF
#endif
    IF(NF90_INQ_VARID(fID,var,vID)==NF90_NoErr) THEN
      CALL get_var2(var,q_glo); q(ijb_u:ije_u,:,iq)=q_glo(ijb_u:ije_u,:); CYCLE
    END IF
    WRITE(lunout,*)TRIM(modname)//": Tracer <"//TRIM(var)//"> is missing"
    WRITE(lunout,*)"         It is hence initialized to zero"
    q(ijb_u:ije_u,:,iq)=0.
   !--- CRisi: for isotops, theoretical initialization using very simplified
   !           Rayleigh distillation las.
    IF(ok_isotopes.AND.iso_num(iq)>0) THEN
      IF(zone_num(iq)==0) q(:,:,iq)=q(:,:,iqpere(iq))*tnat(iso_num(iq))        &
     &           *(q(:,:,iqpere(iq))/30.e-3)**(alpha_ideal(iso_num(iq))-1)
      IF(zone_num(iq)==1) q(:,:,iq)=q(:,:,iqiso(iso_indnum(iq),phase_num(iq)))
    END IF
  END DO
  DEALLOCATE(q_glo)
  CALL err(NF90_CLOSE(fID),"close",fichnom)
  day_ini=day_ini+INT(time)
  time=time-INT(time)


  CONTAINS


SUBROUTINE check_dim(n1,n2,str1,str2)
  INTEGER,          INTENT(IN) :: n1, n2
  CHARACTER(LEN=*), INTENT(IN) :: str1, str2
  CHARACTER(LEN=256) :: s1, s2
  IF(n1/=n2) THEN
    s1='value of '//TRIM(str1)//' ='
    s2=' read in starting file differs from parametrized '//TRIM(str2)//' ='
    WRITE(msg,'(10x,a,i4,2x,a,i4)'),s1,n1,s2,n2
    CALL ABORT_gcm(TRIM(modname),TRIM(msg),1)
  END IF
END SUBROUTINE check_dim


SUBROUTINE get_var1(var,v)
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:)
  REAL,             ALLOCATABLE :: w2(:,:), w3(:,:,:)
  INTEGER :: nn(3), dids(3), k, nd, ntot

  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  ierr=NF90_INQUIRE_VARIABLE(fID,vID,ndims=nd)
  IF(nd==1) THEN
    CALL err(NF90_GET_VAR(fID,vID,v),"get",var); RETURN
  END IF
  ierr=NF90_INQUIRE_VARIABLE(fID,vID,dimids=dids)
  DO k=1,nd; ierr=NF90_INQUIRE_DIMENSION(fID,dids(k),len=nn(k)); END DO
  ntot=PRODUCT(nn(1:nd))
  SELECT CASE(nd)
    CASE(2); ALLOCATE(w2(nn(1),nn(2)))
      CALL err(NF90_GET_VAR(fID,vID,w2),"get",var)
      v=RESHAPE(w2,[ntot]); DEALLOCATE(w2)
    CASE(3); ALLOCATE(w3(nn(1),nn(2),nn(3)))
      CALL err(NF90_GET_VAR(fID,vID,w3),"get",var)
      v=RESHAPE(w3,[ntot]); DEALLOCATE(w3)
  END SELECT
END SUBROUTINE get_var1


SUBROUTINE get_var2(var,v)
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:,:)
  REAL,             ALLOCATABLE :: w4(:,:,:,:)
  INTEGER :: nn(4), dids(4), k, nd
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  ierr=NF90_INQUIRE_VARIABLE(fID,vID,dimids=dids,ndims=nd)
  DO k=1,nd; ierr=NF90_INQUIRE_DIMENSION(fID,dids(k),len=nn(k)); END DO
  ALLOCATE(w4(nn(1),nn(2),nn(3),nn(4)))
  CALL err(NF90_GET_VAR(fID,vID,w4),"get",var)
  v=RESHAPE(w4,[nn(1)*nn(2),nn(3)]); DEALLOCATE(w4)
END SUBROUTINE get_var2


SUBROUTINE err(ierr,typ,nam)
  INTEGER,          INTENT(IN) :: ierr   !--- NetCDF ERROR CODE
  CHARACTER(LEN=*), INTENT(IN) :: typ    !--- TYPE OF OPERATION
  CHARACTER(LEN=*), INTENT(IN) :: nam    !--- FIELD/FILE NAME
  IF(ierr==NF90_NoERR) RETURN
  SELECT CASE(typ)
    CASE('inq');   msg="Field <"//TRIM(nam)//"> is missing"
    CASE('get');   msg="Reading failed for <"//TRIM(nam)//">"
    CASE('open');  msg="File opening failed for <"//TRIM(nam)//">"
    CASE('close'); msg="File closing failed for <"//TRIM(nam)//">"
  END SELECT
  CALL ABORT_gcm(TRIM(modname),TRIM(msg),ierr)
END SUBROUTINE err

END SUBROUTINE dynetat0_loc
