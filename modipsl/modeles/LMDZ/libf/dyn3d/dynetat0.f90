SUBROUTINE dynetat0(fichnom,vcov,ucov,teta,q,masse,ps,phis,time)
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , L.Fairhead
!-------------------------------------------------------------------------------
! Purpose: Initial state reading.
!-------------------------------------------------------------------------------
  USE infotrac
  USE netcdf,      ONLY: NF90_OPEN,  NF90_NOWRITE, NF90_INQ_VARID, NF90_NoErr, &
                         NF90_CLOSE, NF90_GET_VAR
  USE control_mod, ONLY: planet_type
  USE assert_eq_m, ONLY: assert_eq
  USE comvert_mod, ONLY: pa,preff
  USE comconst_mod, ONLY: cpp, daysec, dtvr, g, im, jm, kappa, lllm, omeg, rad
  USE logic_mod, ONLY: fxyhypb, ysinus
  USE serre_mod, ONLY: clon, clat, grossismx, grossismy
  USE temps_mod, ONLY: annee_ref, day_ini, day_ref, itau_dyn, start_time
  USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0

  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom2.h"
  include "description.h"
  include "iniprint.h"
!===============================================================================
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: fichnom          !--- FILE NAME
  REAL, INTENT(OUT) ::  vcov(iip1,jjm, llm)        !--- V COVARIANT WIND
  REAL, INTENT(OUT) ::  ucov(iip1,jjp1,llm)        !--- U COVARIANT WIND
  REAL, INTENT(OUT) ::  teta(iip1,jjp1,llm)        !--- POTENTIAL TEMP.
  REAL, INTENT(OUT) ::     q(iip1,jjp1,llm,nqtot)  !--- TRACERS
  REAL, INTENT(OUT) :: masse(iip1,jjp1,llm)        !--- MASS PER CELL
  REAL, INTENT(OUT) ::    ps(iip1,jjp1)            !--- GROUND PRESSURE
  REAL, INTENT(OUT) ::  phis(iip1,jjp1)            !--- GEOPOTENTIAL
!===============================================================================
! Local variables:
  CHARACTER(LEN=256) :: msg, var, modname
  INTEGER, PARAMETER :: length=100
  INTEGER :: iq, fID, vID, idecal!, iml, jml, lml, nqt
  REAL    :: time, tab_cntrl(length)               !--- RUN PARAMS TABLE
!-------------------------------------------------------------------------------
  modname="dynetat0"

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
  CALL get_var2("cu"   ,cu)
  CALL get_var2("cv"   ,cv)
  CALL get_var2("aire" ,aire)
  var="temps"
  IF(NF90_INQ_VARID(fID,var,vID)/=NF90_NoErr) THEN
    WRITE(lunout,*)TRIM(modname)//": missing field <temps>"
    WRITE(lunout,*)TRIM(modname)//": trying with <Time>"; var="Time"
    CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  END IF
  CALL err(NF90_GET_VAR(fID,vID,time),"get",var)
  CALL get_var2("phisinit",phis)
  CALL get_var3("ucov",ucov)
  CALL get_var3("vcov",vcov)
  CALL get_var3("teta",teta)
  CALL get_var3("masse",masse)
  CALL get_var2("ps",ps)

!--- Tracers
  DO iq=1,nqtot
    var=tname(iq)
    IF(NF90_INQ_VARID(fID,var,vID)==NF90_NoErr) THEN
      CALL err(NF90_GET_VAR(fID,vID,q(:,:,:,iq)),"get",var); CYCLE
    END IF
    WRITE(lunout,*)TRIM(modname)//": Tracer <"//TRIM(var)//"> is missing"
    WRITE(lunout,*)"         It is hence initialized to zero"
    q(:,:,:,iq)=0.
   !--- CRisi: for isotops, theoretical initialization using very simplified
   !           Rayleigh distillation las.
    IF(ok_isotopes.AND.iso_num(iq)>0) THEN
      IF(zone_num(iq)==0) q(:,:,:,iq)=q(:,:,:,iqpere(iq))*tnat(iso_num(iq))    &
     &             *(q(:,:,:,iqpere(iq))/30.e-3)**(alpha_ideal(iso_num(iq))-1)
      IF(zone_num(iq)==1) q(:,:,:,iq)=q(:,:,:,iqiso(iso_indnum(iq),phase_num(iq)))
    END IF
  END DO

  CALL err(NF90_CLOSE(fID),"close",fichnom)
  day_ini=day_ini+INT(time)
  time=time-INT(time)


  CONTAINS


SUBROUTINE check_dim(n1,n2,str1,str2)
  INTEGER,          INTENT(IN) :: n1, n2
  CHARACTER(LEN=*), INTENT(IN) :: str1, str2
  CHARACTER(LEN=100) :: s1, s2
  IF(n1/=n2) THEN
    s1='value of '//TRIM(str1)//' ='
    s2=' read in starting file differs from parametrized '//TRIM(str2)//' ='
    WRITE(msg,'(10x,a,i4,2x,a,i4)'),TRIM(ADJUSTL(s1)),n1,TRIM(ADJUSTL(s2)),n2
    CALL ABORT_gcm(TRIM(modname),TRIM(msg),1)
  END IF
END SUBROUTINE check_dim


SUBROUTINE get_var1(var,v)
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:)
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v),"get",var)
END SUBROUTINE get_var1


SUBROUTINE get_var2(var,v)
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:,:)
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v),"get",var)
END SUBROUTINE get_var2


SUBROUTINE get_var3(var,v)
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:,:,:)
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v),"get",var)
END SUBROUTINE get_var3


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

END SUBROUTINE dynetat0
