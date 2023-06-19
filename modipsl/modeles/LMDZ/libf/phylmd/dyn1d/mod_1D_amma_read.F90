MODULE mod_1D_amma_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Declarations specifiques au cas AMMA
        character*80 :: fich_amma
! Option du cas AMMA ou on impose la discretisation verticale (Ap,Bp)
        integer nlev_amma, nt_amma


        integer year_ini_amma, day_ini_amma, mth_ini_amma
        real heure_ini_amma
        real day_ju_ini_amma   ! Julian day of amma first day
        parameter (year_ini_amma=2006)
        parameter (mth_ini_amma=7)
        parameter (day_ini_amma=10)  ! 10 = 10Juil2006
        parameter (heure_ini_amma=0.) !0h en secondes
        real dt_amma
        parameter (dt_amma=1800.)

!profils initiaux:
        real, allocatable::  plev_amma(:)

        real, allocatable::  z_amma(:)
        real, allocatable::  th_amma(:),q_amma(:)
        real, allocatable::  u_amma(:)
        real, allocatable::  v_amma(:)

        real, allocatable::  th_ammai(:),q_ammai(:)
        real, allocatable::  u_ammai(:)
        real, allocatable::  v_ammai(:)
        real, allocatable::  vitw_ammai(:)
        real, allocatable::  ht_ammai(:)
        real, allocatable::  hq_ammai(:)
        real, allocatable::  vt_ammai(:)
        real, allocatable::  vq_ammai(:)

!forcings
        real, allocatable::  ht_amma(:,:)
        real, allocatable::  hq_amma(:,:)
        real, allocatable::  vitw_amma(:,:)
        real, allocatable::  lat_amma(:),sens_amma(:)

!champs interpoles
        real, allocatable::  vitw_profamma(:)
        real, allocatable::  ht_profamma(:)
        real, allocatable::  hq_profamma(:)
        real lat_profamma,sens_profamma
        real, allocatable::  vt_profamma(:)
        real, allocatable::  vq_profamma(:)
        real, allocatable::  th_profamma(:)
        real, allocatable::  q_profamma(:)
        real, allocatable::  u_profamma(:)
        real, allocatable::  v_profamma(:)


CONTAINS

SUBROUTINE read_1D_cases
      implicit none

#include "netcdf.inc"

      INTEGER nid,rid,ierr

      fich_amma='amma.nc'
      print*,'fich_amma ',fich_amma
      ierr = NF_OPEN(fich_amma,NF_NOWRITE,nid)
      print*,'fich_amma,NF_NOWRITE,nid ',fich_amma,NF_NOWRITE,nid
      if (ierr.NE.NF_NOERR) then
         write(*,*) 'ERROR: GROS Pb opening forcings nc file '
         write(*,*) NF_STRERROR(ierr)
         stop ""
      endif
!.......................................................................
      ierr=NF_INQ_DIMID(nid,'lev',rid)
      IF (ierr.NE.NF_NOERR) THEN
         print*, 'Oh probleme lecture dimension zz'
      ENDIF
      ierr=NF_INQ_DIMLEN(nid,rid,nlev_amma)
      print*,'OK nid,rid,nlev_amma',nid,rid,nlev_amma
!.......................................................................
      ierr=NF_INQ_DIMID(nid,'time',rid)
      print*,'nid,rid',nid,rid
      nt_amma=0
      IF (ierr.NE.NF_NOERR) THEN
        stop 'probleme lecture dimension sens'
      ENDIF
      ierr=NF_INQ_DIMLEN(nid,rid,nt_amma)
      print*,'nid,rid,nlev_amma',nid,rid,nt_amma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!profils initiaux:
        allocate(plev_amma(nlev_amma))
        
        allocate(z_amma(nlev_amma))
        allocate(th_amma(nlev_amma),q_amma(nlev_amma))
        allocate(u_amma(nlev_amma))
        allocate(v_amma(nlev_amma))

!forcings
        allocate(ht_amma(nlev_amma,nt_amma))
        allocate(hq_amma(nlev_amma,nt_amma))
        allocate(vitw_amma(nlev_amma,nt_amma))
        allocate(lat_amma(nt_amma),sens_amma(nt_amma))

!profils initiaux:
        allocate(th_ammai(nlev_amma),q_ammai(nlev_amma))
        allocate(u_ammai(nlev_amma))
        allocate(v_ammai(nlev_amma))
        allocate(vitw_ammai(nlev_amma) )
        allocate(ht_ammai(nlev_amma))
        allocate(hq_ammai(nlev_amma))
        allocate(vt_ammai(nlev_amma))
        allocate(vq_ammai(nlev_amma))

!champs interpoles
        allocate(vitw_profamma(nlev_amma))
        allocate(ht_profamma(nlev_amma))
        allocate(hq_profamma(nlev_amma))
        allocate(vt_profamma(nlev_amma))
        allocate(vq_profamma(nlev_amma))
        allocate(th_profamma(nlev_amma))
        allocate(q_profamma(nlev_amma))
        allocate(u_profamma(nlev_amma))
        allocate(v_profamma(nlev_amma))

        print*,'Allocations OK'
        call read_amma(nid,nlev_amma,nt_amma                                  &
     &     ,z_amma,plev_amma,th_amma,q_amma,u_amma,v_amma,vitw_amma         &
     &     ,ht_amma,hq_amma,sens_amma,lat_amma)

END SUBROUTINE read_1D_cases



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE deallocate_1D_cases
!profils initiaux:
        deallocate(plev_amma)
        
        deallocate(z_amma)
        deallocate(th_amma,q_amma)
        deallocate(u_amma)
        deallocate(v_amma)

        deallocate(th_ammai,q_ammai)
        deallocate(u_ammai)
        deallocate(v_ammai)
        deallocate(vitw_ammai )
        deallocate(ht_ammai)
        deallocate(hq_ammai)
        deallocate(vt_ammai)
        deallocate(vq_ammai)
        
!forcings
        deallocate(ht_amma)
        deallocate(hq_amma)
        deallocate(vitw_amma)
        deallocate(lat_amma,sens_amma)

!champs interpoles
        deallocate(vitw_profamma)
        deallocate(ht_profamma)
        deallocate(hq_profamma)
        deallocate(vt_profamma)
        deallocate(vq_profamma)
        deallocate(th_profamma)
        deallocate(q_profamma)
        deallocate(u_profamma)
        deallocate(v_profamma)
END SUBROUTINE deallocate_1D_cases


END MODULE mod_1D_amma_read
!=====================================================================
      subroutine read_amma(nid,nlevel,ntime                          &
     &     ,zz,pp,temp,qv,u,v,dw                   &
     &     ,dt,dq,sens,flat)

!program reading forcings of the AMMA case study
      implicit none
#include "netcdf.inc"

      integer ntime,nlevel

      real zz(nlevel)
      real temp(nlevel),pp(nlevel)
      real qv(nlevel),u(nlevel)
      real v(nlevel)
      real dw(nlevel,ntime)
      real dt(nlevel,ntime)
      real dq(nlevel,ntime)
      real flat(ntime),sens(ntime)


      integer nid, ierr,rid
      integer nbvar3d
      parameter(nbvar3d=30)
      integer var3didin(nbvar3d)

       ierr=NF_INQ_VARID(nid,"zz",var3didin(1)) 
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'lev'
         endif


      ierr=NF_INQ_VARID(nid,"temp",var3didin(2))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'temp'
         endif

      ierr=NF_INQ_VARID(nid,"qv",var3didin(3))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'qv'
         endif

      ierr=NF_INQ_VARID(nid,"u",var3didin(4))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'u'
         endif

      ierr=NF_INQ_VARID(nid,"v",var3didin(5))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'v'
         endif

      ierr=NF_INQ_VARID(nid,"dw",var3didin(6))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'dw'
         endif

      ierr=NF_INQ_VARID(nid,"dt",var3didin(7))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'dt'
         endif

      ierr=NF_INQ_VARID(nid,"dq",var3didin(8))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'dq'
         endif
      
      ierr=NF_INQ_VARID(nid,"sens",var3didin(9))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'sens'
         endif

      ierr=NF_INQ_VARID(nid,"flat",var3didin(10))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'flat'
         endif

      ierr=NF_INQ_VARID(nid,"pp",var3didin(11))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
      endif

!dimensions lecture
!      call catchaxis(nid,ntime,nlevel,time,z,ierr)
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(1),zz)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(1),zz)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture z ok',zz

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(2),temp)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(2),temp)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture th ok',temp

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(3),qv)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(3),qv)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture qv ok',qv
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(4),u)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(4),u)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture u ok',u

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(5),v)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(5),v)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture v ok',v

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(6),dw)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(6),dw)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture w ok',dw

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(7),dt)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(7),dt)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dt ok',dt

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(8),dq)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(8),dq)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dq ok',dq

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(9),sens)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(9),sens)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture sens ok',sens

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(10),flat)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(10),flat)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture flat ok',flat

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(11),pp)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(11),pp)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture pp ok',pp

         return 
         end subroutine read_amma
!======================================================================
        SUBROUTINE interp_amma_time(day,day1,annee_ref                     &
     &         ,year_ini_amma,day_ini_amma,nt_amma,dt_amma,nlev_amma       &
     &         ,vitw_amma,ht_amma,hq_amma,lat_amma,sens_amma               &
     &         ,vitw_prof,ht_prof,hq_prof,lat_prof,sens_prof)
        implicit none

!---------------------------------------------------------------------------------------
! Time interpolation of a 2D field to the timestep corresponding to day
!
! day: current julian day (e.g. 717538.2)
! day1: first day of the simulation
! nt_amma: total nb of data in the forcing (e.g. 48 for AMMA)
! dt_amma: total time interval (in sec) between 2 forcing data (e.g. 30min for AMMA)
!---------------------------------------------------------------------------------------

#include "compar1d.h"

! inputs:
        integer annee_ref
        integer nt_amma,nlev_amma
        integer year_ini_amma
        real day, day1,day_ini_amma,dt_amma
        real vitw_amma(nlev_amma,nt_amma)
        real ht_amma(nlev_amma,nt_amma)
        real hq_amma(nlev_amma,nt_amma)
        real lat_amma(nt_amma)
        real sens_amma(nt_amma)
! outputs:
        real vitw_prof(nlev_amma)
        real ht_prof(nlev_amma)
        real hq_prof(nlev_amma)
        real lat_prof,sens_prof
! local:
        integer it_amma1, it_amma2,k
        real timeit,time_amma1,time_amma2,frac


        if (forcing_type.eq.6) then
! Check that initial day of the simulation consistent with AMMA case:
       if (annee_ref.ne.2006) then
        print*,'Pour AMMA, annee_ref doit etre 2006'
        print*,'Changer annee_ref dans run.def'
        stop
       endif
       if (annee_ref.eq.2006 .and. day1.lt.day_ini_amma) then
        print*,'AMMA a débuté le 10 juillet 2006',day1,day_ini_amma
        print*,'Changer dayref dans run.def'
        stop
       endif
       if (annee_ref.eq.2006 .and. day1.gt.day_ini_amma+1) then
        print*,'AMMA a fini le 11 juillet'
        print*,'Changer dayref ou nday dans run.def'
        stop
       endif
       endif

! Determine timestep relative to the 1st day of AMMA:
!       timeit=(day-day1)*86400.
!       if (annee_ref.eq.1992) then
!        timeit=(day-day_ini_toga)*86400.
!       else
!        timeit=(day+61.-1.)*86400. ! 61 days between Nov01 and Dec31 1992
!       endif
      timeit=(day-day_ini_amma)*86400

! Determine the closest observation times:
!       it_amma1=INT(timeit/dt_amma)+1
!       it_amma2=it_amma1 + 1
!       time_amma1=(it_amma1-1)*dt_amma
!       time_amma2=(it_amma2-1)*dt_amma

       it_amma1=INT(timeit/dt_amma)+1
       IF (it_amma1 .EQ. nt_amma) THEN
       it_amma2=it_amma1 
       ELSE
       it_amma2=it_amma1 + 1
       ENDIF
       time_amma1=(it_amma1-1)*dt_amma
       time_amma2=(it_amma2-1)*dt_amma

       if (it_amma1 .gt. nt_amma) then
        write(*,*) 'PB-stop: day, it_amma1, it_amma2, timeit: '            &
     &        ,day,day_ini_amma,it_amma1,it_amma2,timeit/86400.
        stop
       endif

! time interpolation:
       IF (it_amma1 .EQ. it_amma2) THEN
          frac=0.
       ELSE 
          frac=(time_amma2-timeit)/(time_amma2-time_amma1)
          frac=max(frac,0.0)
       ENDIF

       lat_prof = lat_amma(it_amma2)                                       &
     &          -frac*(lat_amma(it_amma2)-lat_amma(it_amma1)) 
       sens_prof = sens_amma(it_amma2)                                     &
     &          -frac*(sens_amma(it_amma2)-sens_amma(it_amma1))

       do k=1,nlev_amma
        vitw_prof(k) = vitw_amma(k,it_amma2)                               &
     &          -frac*(vitw_amma(k,it_amma2)-vitw_amma(k,it_amma1))
        ht_prof(k) = ht_amma(k,it_amma2)                                   &
     &          -frac*(ht_amma(k,it_amma2)-ht_amma(k,it_amma1))
        hq_prof(k) = hq_amma(k,it_amma2)                                   &
     &          -frac*(hq_amma(k,it_amma2)-hq_amma(k,it_amma1))
        enddo

        return
        END

