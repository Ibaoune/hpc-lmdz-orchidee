SUBROUTINE dynredem0_loc(fichnom,iday_end,phis)
!
!-------------------------------------------------------------------------------
! Write the NetCDF restart file (initialization).
!-------------------------------------------------------------------------------
#ifdef CPP_IOIPSL
  USE IOIPSL
#endif
  USE parallel_lmdz
  USE mod_hallo
  USE infotrac
  USE netcdf, ONLY: NF90_CREATE, NF90_DEF_DIM, NF90_INQ_VARID, NF90_GLOBAL,    &
                    NF90_CLOSE,  NF90_PUT_ATT, NF90_UNLIMITED, NF90_CLOBBER
  USE dynredem_mod, ONLY: cre_var, put_var, err, modname, fil
  USE comvert_mod, ONLY: ap,bp,aps,bps,presnivs,pseudoalt,pa,preff, &
                         nivsig,nivsigs
  USE comconst_mod, ONLY: cpp, daysec, dtvr, g, kappa, omeg, rad
  USE logic_mod, ONLY: fxyhypb, ysinus
  USE serre_mod, ONLY: clon,clat,grossismx,grossismy,dzoomx,dzoomy, &
                       taux,tauy
  USE temps_mod, ONLY: annee_ref, day_ref, itau_dyn, itaufin, start_time
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
  INTEGER,          INTENT(IN) :: iday_end         !--- 
  REAL,             INTENT(IN) :: phis(ijb_u:ije_u)!--- GROUND GEOPOTENTIAL
!===============================================================================
! Local variables:
  INTEGER :: iq, l
  INTEGER, PARAMETER :: length=100
  REAL    :: tab_cntrl(length)                     !--- RUN PARAMETERS TABLE
  REAL    :: phis_glo(ip1jmp1)
!   For NetCDF:
  CHARACTER(LEN=30) :: unites
  INTEGER :: indexID
  INTEGER :: rlonuID, rlonvID, rlatuID, rlatvID
  INTEGER :: sID, sigID, nID, vID, timID
  INTEGER :: yyears0, jjour0, mmois0
  REAL    :: zan0, zjulian, hours
!===============================================================================
  modname='dynredem0'; fil=fichnom
  CALL Gather_field_u(phis,phis_glo,1)
  IF(mpi_rank/=0) RETURN

#ifdef CPP_IOIPSL
  CALL ymds2ju(annee_ref, 1, iday_end, 0.0, zjulian)
  CALL ju2ymds(zjulian, yyears0, mmois0, jjour0, hours)
#else
! set yyears0, mmois0, jjour0 to 0,1,1 (hours is not used)
  yyears0=0
  mmois0=1
  jjour0=1
#endif        

  tab_cntrl(:)  = 0.
  tab_cntrl(1)  = REAL(iim)
  tab_cntrl(2)  = REAL(jjm)
  tab_cntrl(3)  = REAL(llm)
  tab_cntrl(4)  = REAL(day_ref)
  tab_cntrl(5)  = REAL(annee_ref)
  tab_cntrl(6)  = rad
  tab_cntrl(7)  = omeg
  tab_cntrl(8)  = g
  tab_cntrl(9)  = cpp
  tab_cntrl(10) = kappa
  tab_cntrl(11) = daysec
  tab_cntrl(12) = dtvr
  tab_cntrl(13) = etot0
  tab_cntrl(14) = ptot0
  tab_cntrl(15) = ztot0
  tab_cntrl(16) = stot0
  tab_cntrl(17) = ang0
  tab_cntrl(18) = pa
  tab_cntrl(19) = preff

!    .....    parameters for zoom    ......   
  tab_cntrl(20) = clon
  tab_cntrl(21) = clat
  tab_cntrl(22) = grossismx
  tab_cntrl(23) = grossismy
!
  IF ( fxyhypb )   THEN
    tab_cntrl(24) = 1.
    tab_cntrl(25) = dzoomx
    tab_cntrl(26) = dzoomy
    tab_cntrl(27) = 0.
    tab_cntrl(28) = taux
    tab_cntrl(29) = tauy
  ELSE
    tab_cntrl(24) = 0.
    tab_cntrl(25) = dzoomx
    tab_cntrl(26) = dzoomy
    tab_cntrl(27) = 0.
    tab_cntrl(28) = 0.
    tab_cntrl(29) = 0.
    IF( ysinus )  tab_cntrl(27) = 1.
  END IF
  tab_cntrl(30) = REAL(iday_end)
  tab_cntrl(31) = REAL(itau_dyn + itaufin)
! start_time: start_time of simulation (not necessarily 0.)
  tab_cntrl(32) = start_time

!--- File creation
  CALL err(NF90_CREATE(fichnom,NF90_CLOBBER,nid))

!--- Some global attributes
  CALL err(NF90_PUT_ATT(nid,NF90_GLOBAL,"title","Fichier demarrage dynamique"))

!--- Dimensions
  CALL err(NF90_DEF_DIM(nid,"index", length, indexID))
  CALL err(NF90_DEF_DIM(nid,"rlonu", iip1,   rlonuID))
  CALL err(NF90_DEF_DIM(nid,"rlatu", jjp1,   rlatuID))
  CALL err(NF90_DEF_DIM(nid,"rlonv", iip1,   rlonvID))
  CALL err(NF90_DEF_DIM(nid,"rlatv", jjm,    rlatvID))
  CALL err(NF90_DEF_DIM(nid,"sigs",  llm,        sID))
  CALL err(NF90_DEF_DIM(nid,"sig",   llmp1,    sigID))
  CALL err(NF90_DEF_DIM(nid,"temps", NF90_UNLIMITED, timID))

!--- Define and save invariant fields
  CALL put_var(nid,"controle","Parametres de controle" ,[indexID],tab_cntrl)
  CALL put_var(nid,"rlonu"   ,"Longitudes des points U",[rlonuID],rlonu)
  CALL put_var(nid,"rlatu"   ,"Latitudes des points U" ,[rlatuID],rlatu)
  CALL put_var(nid,"rlonv"   ,"Longitudes des points V",[rlonvID],rlonv)
  CALL put_var(nid,"rlatv"   ,"Latitudes des points V" ,[rlatvID],rlatv)
  CALL put_var(nid,"nivsigs" ,"Numero naturel des couches s"    ,[sID]  ,nivsigs)
  CALL put_var(nid,"nivsig"  ,"Numero naturel des couches sigma",[sigID],nivsig)
  CALL put_var(nid,"ap"      ,"Coefficient A pour hybride"      ,[sigID],ap)
  CALL put_var(nid,"bp"      ,"Coefficient B pour hybride"      ,[sigID],bp)
  CALL put_var(nid,"presnivs",""                                ,[sID]  ,presnivs)
! covariant <-> contravariant <-> natural conversion coefficients
  CALL put_var(nid,"cu","Coefficient de passage pour U",[rlonuID,rlatuID],cu)
  CALL put_var(nid,"cv","Coefficient de passage pour V",[rlonvID,rlatvID],cv)
  CALL put_var(nid,"aire","Aires de chaque maille"     ,[rlonvID,rlatuID],aire)
  CALL put_var(nid,"phisinit","Geopotentiel au sol"    ,[rlonvID,rlatuID],phis_glo)

!--- Define fields saved later
  WRITE(unites,"('days since ',i4,'-',i2.2,'-',i2.2,' 00:00:00')"),&
               yyears0,mmois0,jjour0
  CALL cre_var(nid,"temps","Temps de simulation",[timID],unites)
  CALL cre_var(nid,"ucov" ,"Vitesse U"  ,[rlonuID,rlatuID,sID,timID])
  CALL cre_var(nid,"vcov" ,"Vitesse V"  ,[rlonvID,rlatvID,sID,timID])
  CALL cre_var(nid,"teta" ,"Temperature",[rlonvID,rlatuID,sID,timID])
  DO iq=1,nqtot
    CALL cre_var(nid,tname(iq),ttext(iq),[rlonvID,rlatuID,sID,timID])
  END DO
  CALL cre_var(nid,"masse","Masse d air"    ,[rlonvID,rlatuID,sID,timID])
  CALL cre_var(nid,"ps"   ,"Pression au sol",[rlonvID,rlatuID    ,timID])
  CALL err(NF90_CLOSE (nid))

  WRITE(lunout,*)TRIM(modname)//': iim,jjm,llm,iday_end',iim,jjm,llm,iday_end
  WRITE(lunout,*)TRIM(modname)//': rad,omeg,g,cpp,kappa',rad,omeg,g,cpp,kappa

END SUBROUTINE dynredem0_loc
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE dynredem1_loc(fichnom,time,vcov,ucov,teta,q,masse,ps)
!
!-------------------------------------------------------------------------------
! Purpose: Write the NetCDF restart file (append).
!-------------------------------------------------------------------------------
  USE parallel_lmdz
  USE mod_hallo
  USE infotrac
  USE control_mod
  USE netcdf,   ONLY: NF90_OPEN,  NF90_NOWRITE, NF90_GET_VAR, NF90_INQ_VARID,  &
                      NF90_CLOSE, NF90_WRITE,   NF90_PUT_VAR, NF90_NoErr
  USE dynredem_mod, ONLY: dynredem_write_u, dynredem_write_v, dynredem_read_u, &
                          err, modname, fil, msg
  USE temps_mod, ONLY: itau_dyn, itaufin
  
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "description.h"
  include "comgeom.h"
  include "iniprint.h"
!===============================================================================
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: fichnom              !-- FILE NAME
  REAL, INTENT(IN)    ::  time                         !-- TIME
  REAL, INTENT(IN)    ::  vcov(ijb_v:ije_v,llm)        !-- V COVARIANT WIND
  REAL, INTENT(IN)    ::  ucov(ijb_u:ije_u,llm)        !-- U COVARIANT WIND
  REAL, INTENT(IN)    ::  teta(ijb_u:ije_u,llm)        !-- POTENTIAL TEMPERATURE
  REAL, INTENT(INOUT) ::     q(ijb_u:ije_u,llm,nqtot)  !-- TRACERS
  REAL, INTENT(IN)    :: masse(ijb_u:ije_u,llm)        !-- MASS PER CELL
  REAL, INTENT(IN)    ::    ps(ijb_u:ije_u)            !-- GROUND PRESSURE
!===============================================================================
! Local variables:
  INTEGER :: l, iq, nid, vID, ierr, nid_trac, vID_trac
  INTEGER, SAVE :: nb=0
  INTEGER, PARAMETER :: length=100
  REAL               :: tab_cntrl(length) ! tableau des parametres du run
  CHARACTER(LEN=256) :: var, dum
  LOGICAL            :: lread_inca
!===============================================================================

!$OMP MASTER
  IF(mpi_rank==0) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  modname='dynredem1_loc'; fil=fichnom
  fil=fichnom
  CALL err(NF90_OPEN(fil,NF90_WRITE,nid),"open",fil)

!--- Write/extend time coordinate
  nb = nb + 1
  var="temps"
  CALL err(NF90_INQ_VARID(nid,var,vID),"inq",var)
  CALL err(NF90_PUT_VAR(nid,vID,[time]),"put",var)
  WRITE(lunout,*)TRIM(modname)//": Saving for ", nb, time

!--- Rewrite control table (itaufin undefined in dynredem0)
  var="controle"
  CALL err(NF90_INQ_VARID(nid,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(nid,vID,tab_cntrl),"get",var)
  tab_cntrl(31)=DBLE(itau_dyn + itaufin)
  CALL err(NF90_INQ_VARID(nid,var,vID),"inq",var)
  CALL err(NF90_PUT_VAR(nid,vID,tab_cntrl),"put",var)
  END IF               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!$OMP END MASTER

!--- Save fields
  CALL dynredem_write_u(nid,"ucov" ,ucov ,llm)
  CALL dynredem_write_v(nid,"vcov" ,vcov ,llm)
  CALL dynredem_write_u(nid,"teta" ,teta ,llm)
  CALL dynredem_write_u(nid,"masse",masse,llm)
  CALL dynredem_write_u(nid,"ps"   ,ps   ,1)

!--- Tracers in file "start_trac.nc" (added by Anne)
  lread_inca=.FALSE.
!$OMP MASTER
  fil="start_trac.nc"
  IF(type_trac=='inca') INQUIRE(FILE=fil,EXIST=lread_inca)
  IF(lread_inca) CALL err(NF90_OPEN(fil,NF90_NOWRITE,nid_trac),"open")
!$OMP END MASTER
!$OMP BARRIER

!--- Save tracers
  DO iq=1,nqtot; var=tname(iq); ierr=-1
    IF(lread_inca) THEN                  !--- Possibly read from "start_trac.nc"
!$OMP MASTER      
      fil="start_trac.nc"
      ierr=NF90_INQ_VARID(nid_trac,var,vID_trac)
      dum='inq'; IF(ierr==NF90_NoErr) dum='fnd'
      WRITE(lunout,*)msg(dum,var)
!$OMP END MASTER
!$OMP BARRIER
      IF(ierr==NF90_NoErr) CALL dynredem_read_u(nid_trac,var,q(:,:,iq),llm)
    END IF
    fil=fichnom
    CALL dynredem_write_u(nid,var,q(:,:,iq),llm)
  END DO

!$OMP MASTER
  IF(mpi_rank==0) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL err(NF90_CLOSE(nid),"close")
  fil="start_trac.nc"
  IF(lread_inca) CALL err(NF90_CLOSE(nid_trac),"close")
  END IF               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!$OMP END MASTER

END SUBROUTINE dynredem1_loc

