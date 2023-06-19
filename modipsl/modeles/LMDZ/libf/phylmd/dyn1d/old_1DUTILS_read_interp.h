!======================================================================
      SUBROUTINE read_togacoare(fich_toga,nlev_toga,nt_toga                     &
     &             ,ts_toga,plev_toga,t_toga,q_toga,u_toga,v_toga,w_toga        &
     &             ,ht_toga,vt_toga,hq_toga,vq_toga)
      implicit none

!-------------------------------------------------------------------------
! Read TOGA-COARE forcing data 
!-------------------------------------------------------------------------

      integer nlev_toga,nt_toga
      real ts_toga(nt_toga),plev_toga(nlev_toga,nt_toga)
      real t_toga(nlev_toga,nt_toga),q_toga(nlev_toga,nt_toga)
      real u_toga(nlev_toga,nt_toga),v_toga(nlev_toga,nt_toga)
      real w_toga(nlev_toga,nt_toga)
      real ht_toga(nlev_toga,nt_toga),vt_toga(nlev_toga,nt_toga)
      real hq_toga(nlev_toga,nt_toga),vq_toga(nlev_toga,nt_toga)
      character*80 fich_toga

      integer k,ip
      real bid

      integer iy,im,id,ih
      
       real plev_min

       plev_min = 55.  ! pas de tendance de vap. d eau au-dessus de 55 hPa

      open(21,file=trim(fich_toga),form='formatted')
      read(21,'(a)') 
      do ip = 1, nt_toga
      read(21,'(a)') 
      read(21,'(a)') 
      read(21,223) iy, im, id, ih, bid, ts_toga(ip), bid,bid,bid,bid
      read(21,'(a)') 
      read(21,'(a)') 

       do k = 1, nlev_toga
         read(21,230) plev_toga(k,ip), t_toga(k,ip), q_toga(k,ip)          &
     &       ,u_toga(k,ip), v_toga(k,ip), w_toga(k,ip)                     &
     &       ,ht_toga(k,ip), vt_toga(k,ip), hq_toga(k,ip), vq_toga(k,ip)

! conversion in SI units:
         t_toga(k,ip)=t_toga(k,ip)+273.15     ! K
         q_toga(k,ip)=q_toga(k,ip)*0.001      ! kg/kg
         w_toga(k,ip)=w_toga(k,ip)*100./3600. ! Pa/s
! no water vapour tendency above 55 hPa
         if (plev_toga(k,ip) .lt. plev_min) then
          q_toga(k,ip) = 0.
          hq_toga(k,ip) = 0.
          vq_toga(k,ip) =0.
         endif
       enddo

         ts_toga(ip)=ts_toga(ip)+273.15       ! K
       enddo
       close(21)

  223 format(4i3,6f8.2)
  230 format(6f9.3,4e11.3)

          return
          end

!-------------------------------------------------------------------------
      SUBROUTINE read_sandu(fich_sandu,nlev_sandu,nt_sandu,ts_sandu)
      implicit none

!-------------------------------------------------------------------------
! Read I.SANDU case forcing data
!-------------------------------------------------------------------------

      integer nlev_sandu,nt_sandu
      real ts_sandu(nt_sandu)
      character*80 fich_sandu

      integer ip
      integer iy,im,id,ih

      real plev_min

      print*,'nlev_sandu',nlev_sandu
      plev_min = 55000.  ! pas de tendance de vap. d eau au-dessus de 55 hPa

      open(21,file=trim(fich_sandu),form='formatted')
      read(21,'(a)')
      do ip = 1, nt_sandu
      read(21,'(a)')
      read(21,'(a)')
      read(21,223) iy, im, id, ih, ts_sandu(ip)
      print *,'ts=',iy,im,id,ih,ip,ts_sandu(ip)
      enddo
      close(21)

  223 format(4i3,f8.2)

          return
          end

!=====================================================================
!-------------------------------------------------------------------------
      SUBROUTINE read_astex(fich_astex,nlev_astex,nt_astex,div_astex,      &
     & ts_astex,ug_astex,vg_astex,ufa_astex,vfa_astex)
      implicit none

!-------------------------------------------------------------------------
! Read Astex case forcing data
!-------------------------------------------------------------------------

      integer nlev_astex,nt_astex
      real div_astex(nt_astex),ts_astex(nt_astex),ug_astex(nt_astex)
      real vg_astex(nt_astex),ufa_astex(nt_astex),vfa_astex(nt_astex)
      character*80 fich_astex

      integer ip
      integer iy,im,id,ih

       real plev_min

      print*,'nlev_astex',nlev_astex
       plev_min = 55000.  ! pas de tendance de vap. d eau au-dessus de 55 hPa

      open(21,file=trim(fich_astex),form='formatted')
      read(21,'(a)')
      read(21,'(a)')
      do ip = 1, nt_astex
      read(21,'(a)')
      read(21,'(a)')
      read(21,223) iy, im, id, ih, div_astex(ip),ts_astex(ip),             &
     &ug_astex(ip),vg_astex(ip),ufa_astex(ip),vfa_astex(ip)
      ts_astex(ip)=ts_astex(ip)+273.15
      print *,'ts=',iy,im,id,ih,ip,div_astex(ip),ts_astex(ip),             &
     &ug_astex(ip),vg_astex(ip),ufa_astex(ip),vg_astex(ip)
      enddo
      close(21)

  223 format(4i3,e13.2,f7.2,f7.3,f7.2,f7.3,f7.2)

          return
          end
!=====================================================================
      subroutine read_twpice(fich_twpice,nlevel,ntime                       &
     &     ,T_srf,plev,T,q,u,v,omega                                       &
     &     ,T_adv_h,T_adv_v,q_adv_h,q_adv_v)

!program reading forcings of the TWP-ICE experiment

!      use netcdf

      implicit none

#include "netcdf.inc"

      integer ntime,nlevel
      integer l,k
      character*80 :: fich_twpice
      real*8 time(ntime)
      real*8 lat, lon, alt, phis
      real*8 lev(nlevel)
      real*8 plev(nlevel,ntime)

      real*8 T(nlevel,ntime)
      real*8 q(nlevel,ntime),u(nlevel,ntime)
      real*8 v(nlevel,ntime)
      real*8 omega(nlevel,ntime), div(nlevel,ntime)
      real*8 T_adv_h(nlevel,ntime)
      real*8 T_adv_v(nlevel,ntime), q_adv_h(nlevel,ntime)
      real*8 q_adv_v(nlevel,ntime)
      real*8 s(nlevel,ntime), s_adv_h(nlevel,ntime)
      real*8 s_adv_v(nlevel,ntime)
      real*8 p_srf_aver(ntime), p_srf_center(ntime)
      real*8 T_srf(ntime)

      integer nid, ierr
      integer nbvar3d
      parameter(nbvar3d=20)
      integer var3didin(nbvar3d)

      ierr = NF_OPEN(fich_twpice,NF_NOWRITE,nid)
      if (ierr.NE.NF_NOERR) then
         write(*,*) 'ERROR: Pb opening forcings cdf file '
         write(*,*) NF_STRERROR(ierr)
         stop ""
      endif

      ierr=NF_INQ_VARID(nid,"lat",var3didin(1))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'lat'
         endif
      
       ierr=NF_INQ_VARID(nid,"lon",var3didin(2))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'lon'
         endif

       ierr=NF_INQ_VARID(nid,"alt",var3didin(3))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'alt'
         endif

      ierr=NF_INQ_VARID(nid,"phis",var3didin(4))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'phis'
         endif

      ierr=NF_INQ_VARID(nid,"T",var3didin(5))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'T'
         endif

      ierr=NF_INQ_VARID(nid,"q",var3didin(6))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'q'
         endif

      ierr=NF_INQ_VARID(nid,"u",var3didin(7))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'u'
         endif

      ierr=NF_INQ_VARID(nid,"v",var3didin(8))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'v'
         endif

      ierr=NF_INQ_VARID(nid,"omega",var3didin(9))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'omega'
         endif

      ierr=NF_INQ_VARID(nid,"div",var3didin(10))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'div'
         endif

      ierr=NF_INQ_VARID(nid,"T_adv_h",var3didin(11))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'T_adv_h'
         endif

      ierr=NF_INQ_VARID(nid,"T_adv_v",var3didin(12))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'T_adv_v'
         endif

      ierr=NF_INQ_VARID(nid,"q_adv_h",var3didin(13))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'q_adv_h'
         endif

      ierr=NF_INQ_VARID(nid,"q_adv_v",var3didin(14))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'q_adv_v'
         endif

      ierr=NF_INQ_VARID(nid,"s",var3didin(15))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 's'
         endif

      ierr=NF_INQ_VARID(nid,"s_adv_h",var3didin(16))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 's_adv_h'
         endif
    
      ierr=NF_INQ_VARID(nid,"s_adv_v",var3didin(17))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 's_adv_v'
         endif

      ierr=NF_INQ_VARID(nid,"p_srf_aver",var3didin(18))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'p_srf_aver'
         endif

      ierr=NF_INQ_VARID(nid,"p_srf_center",var3didin(19))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'p_srf_center'
         endif

      ierr=NF_INQ_VARID(nid,"T_srf",var3didin(20))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'T_srf'
         endif

!dimensions lecture
      call catchaxis(nid,ntime,nlevel,time,lev,ierr)

!pressure 
       do l=1,ntime
       do k=1,nlevel
          plev(k,l)=lev(k)
       enddo
       enddo
          
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(1),lat)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(1),lat)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture lat ok',lat

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(2),lon)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(2),lon)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture lon ok',lon
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(3),alt)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(3),alt)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture alt ok',alt
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(4),phis)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(4),phis)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture phis ok',phis
          
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(5),T)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(5),T)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture T ok'

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(6),q)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(6),q)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture q ok'
!q in kg/kg
       do l=1,ntime
       do k=1,nlevel
          q(k,l)=q(k,l)/1000.
       enddo
       enddo
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(7),u)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(7),u)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture u ok'

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(8),v)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(8),v)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture v ok'

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(9),omega)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(9),omega)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture omega ok'
!omega in mb/hour
       do l=1,ntime
       do k=1,nlevel
          omega(k,l)=omega(k,l)*100./3600.
       enddo
       enddo

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(10),div)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(10),div)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture div ok'

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(11),T_adv_h)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(11),T_adv_h)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture T_adv_h ok'
!T adv in K/s
       do l=1,ntime
       do k=1,nlevel
          T_adv_h(k,l)=T_adv_h(k,l)/3600.
       enddo
       enddo


#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(12),T_adv_v)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(12),T_adv_v)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture T_adv_v ok'
!T adv in K/s
       do l=1,ntime
       do k=1,nlevel
          T_adv_v(k,l)=T_adv_v(k,l)/3600.
       enddo
       enddo

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(13),q_adv_h)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(13),q_adv_h)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture q_adv_h ok'
!q adv in kg/kg/s
       do l=1,ntime
       do k=1,nlevel
          q_adv_h(k,l)=q_adv_h(k,l)/1000./3600.
       enddo
       enddo


#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(14),q_adv_v)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(14),q_adv_v)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture q_adv_v ok'
!q adv in kg/kg/s
       do l=1,ntime
       do k=1,nlevel
          q_adv_v(k,l)=q_adv_v(k,l)/1000./3600.
       enddo
       enddo


#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(15),s)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(15),s)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(16),s_adv_h)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(16),s_adv_h)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(17),s_adv_v)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(17),s_adv_v)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(18),p_srf_aver)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(18),p_srf_aver)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(19),p_srf_center)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(19),p_srf_center)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(20),T_srf)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(20),T_srf)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!         write(*,*)'lecture T_srf ok', T_srf

         return 
         end subroutine read_twpice
!=====================================================================
         subroutine catchaxis(nid,ttm,llm,time,lev,ierr)

!         use netcdf

         implicit none
#include "netcdf.inc"
         integer nid,ttm,llm
         real*8 time(ttm)
         real*8 lev(llm)
         integer ierr

         integer timevar,levvar
         integer timelen,levlen
         integer timedimin,levdimin

! Control & lecture on dimensions
! ===============================
         ierr=NF_INQ_DIMID(nid,"time",timedimin)
         ierr=NF_INQ_VARID(nid,"time",timevar)
         if (ierr.NE.NF_NOERR) then
            write(*,*) 'ERROR: Field <time> is missing'
            stop ""  
         endif
         ierr=NF_INQ_DIMLEN(nid,timedimin,timelen)

         ierr=NF_INQ_DIMID(nid,"lev",levdimin)
         ierr=NF_INQ_VARID(nid,"lev",levvar)
         if (ierr.NE.NF_NOERR) then
             write(*,*) 'ERROR: Field <lev> is lacking'
             stop "" 
         endif
         ierr=NF_INQ_DIMLEN(nid,levdimin,levlen)

         if((timelen/=ttm).or.(levlen/=llm)) then
            write(*,*) 'ERROR: Not the good lenght for axis'
            write(*,*) 'longitude: ',timelen,ttm+1
            write(*,*) 'latitude: ',levlen,llm
            stop ""  
         endif

!#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,timevar,time)
         ierr = NF_GET_VAR_DOUBLE(nid,levvar,lev)
!#else
!        ierr = NF_GET_VAR_REAL(nid,timevar,time)
!        ierr = NF_GET_VAR_REAL(nid,levvar,lev)
!#endif

       return
       end
!=====================================================================

       SUBROUTINE interp_sandu_vertical(play,nlev_sandu,plev_prof          &
     &         ,t_prof,thl_prof,q_prof,u_prof,v_prof,w_prof                &
     &         ,omega_prof,o3mmr_prof                                      &
     &         ,t_mod,thl_mod,q_mod,u_mod,v_mod,w_mod                      &
     &         ,omega_mod,o3mmr_mod,mxcalc)

       implicit none

#include "dimensions.h"

!-------------------------------------------------------------------------
! Vertical interpolation of SANDUREF forcing data onto model levels
!-------------------------------------------------------------------------

       integer nlevmax
       parameter (nlevmax=41)
       integer nlev_sandu,mxcalc
!       real play(llm), plev_prof(nlevmax)
!       real t_prof(nlevmax),q_prof(nlevmax)
!       real u_prof(nlevmax),v_prof(nlevmax), w_prof(nlevmax)
!       real ht_prof(nlevmax),vt_prof(nlevmax)
!       real hq_prof(nlevmax),vq_prof(nlevmax)

       real play(llm), plev_prof(nlev_sandu)
       real t_prof(nlev_sandu),thl_prof(nlev_sandu),q_prof(nlev_sandu)
       real u_prof(nlev_sandu),v_prof(nlev_sandu), w_prof(nlev_sandu)
       real omega_prof(nlev_sandu),o3mmr_prof(nlev_sandu)

       real t_mod(llm),thl_mod(llm),q_mod(llm)
       real u_mod(llm),v_mod(llm), w_mod(llm)
       real omega_mod(llm),o3mmr_mod(llm)

       integer l,k,k1,k2
       real frac,frac1,frac2,fact

       do l = 1, llm

        if (play(l).ge.plev_prof(nlev_sandu)) then

        mxcalc=l
         k1=0
         k2=0

         if (play(l).le.plev_prof(1)) then

         do k = 1, nlev_sandu-1
          if (play(l).le.plev_prof(k).and. play(l).gt.plev_prof(k+1)) then
            k1=k
            k2=k+1
          endif
         enddo

         if (k1.eq.0 .or. k2.eq.0) then
          write(*,*) 'PB! k1, k2 = ',k1,k2
          write(*,*) 'l,play(l) = ',l,play(l)/100
         do k = 1, nlev_sandu-1
          write(*,*) 'k,plev_prof(k) = ',k,plev_prof(k)/100
         enddo
         endif

         frac = (plev_prof(k2)-play(l))/(plev_prof(k2)-plev_prof(k1))
         t_mod(l)= t_prof(k2) - frac*(t_prof(k2)-t_prof(k1))
         thl_mod(l)= thl_prof(k2) - frac*(thl_prof(k2)-thl_prof(k1))
         q_mod(l)= q_prof(k2) - frac*(q_prof(k2)-q_prof(k1))
         u_mod(l)= u_prof(k2) - frac*(u_prof(k2)-u_prof(k1))
         v_mod(l)= v_prof(k2) - frac*(v_prof(k2)-v_prof(k1))
         w_mod(l)= w_prof(k2) - frac*(w_prof(k2)-w_prof(k1))
         omega_mod(l)=omega_prof(k2)-frac*(omega_prof(k2)-omega_prof(k1))
         o3mmr_mod(l)=o3mmr_prof(k2)-frac*(o3mmr_prof(k2)-o3mmr_prof(k1))

         else !play>plev_prof(1)

         k1=1
         k2=2
         frac1 = (play(l)-plev_prof(k2))/(plev_prof(k1)-plev_prof(k2))
         frac2 = (play(l)-plev_prof(k1))/(plev_prof(k1)-plev_prof(k2))
         t_mod(l)= frac1*t_prof(k1) - frac2*t_prof(k2)
         thl_mod(l)= frac1*thl_prof(k1) - frac2*thl_prof(k2)
         q_mod(l)= frac1*q_prof(k1) - frac2*q_prof(k2)
         u_mod(l)= frac1*u_prof(k1) - frac2*u_prof(k2)
         v_mod(l)= frac1*v_prof(k1) - frac2*v_prof(k2)
         w_mod(l)= frac1*w_prof(k1) - frac2*w_prof(k2)
         omega_mod(l)= frac1*omega_prof(k1) - frac2*omega_prof(k2)
         o3mmr_mod(l)= frac1*o3mmr_prof(k1) - frac2*o3mmr_prof(k2)

         endif ! play.le.plev_prof(1)

        else ! above max altitude of forcing file

!jyg
         fact=20.*(plev_prof(nlev_sandu)-play(l))/plev_prof(nlev_sandu) !jyg
         fact = max(fact,0.)                                           !jyg
         fact = exp(-fact)                                             !jyg
         t_mod(l)= t_prof(nlev_sandu)                                   !jyg
         thl_mod(l)= thl_prof(nlev_sandu)                                   !jyg
         q_mod(l)= q_prof(nlev_sandu)*fact                              !jyg
         u_mod(l)= u_prof(nlev_sandu)*fact                              !jyg
         v_mod(l)= v_prof(nlev_sandu)*fact                              !jyg
         w_mod(l)= w_prof(nlev_sandu)*fact                              !jyg
         omega_mod(l)= omega_prof(nlev_sandu)*fact                      !jyg
         o3mmr_mod(l)= o3mmr_prof(nlev_sandu)*fact                      !jyg

        endif ! play

       enddo ! l

       do l = 1,llm
!      print *,'t_mod(l),thl_mod(l),q_mod(l),u_mod(l),v_mod(l) ',
!    $        l,t_mod(l),thl_mod(l),q_mod(l),u_mod(l),v_mod(l)
       enddo

          return
          end
!=====================================================================
       SUBROUTINE interp_astex_vertical(play,nlev_astex,plev_prof          &
     &         ,t_prof,thl_prof,qv_prof,ql_prof,qt_prof,u_prof,v_prof      &
     &         ,w_prof,tke_prof,o3mmr_prof                                 &
     &         ,t_mod,thl_mod,qv_mod,ql_mod,qt_mod,u_mod,v_mod,w_mod       &
     &         ,tke_mod,o3mmr_mod,mxcalc)

       implicit none

#include "dimensions.h"

!-------------------------------------------------------------------------
! Vertical interpolation of Astex forcing data onto model levels
!-------------------------------------------------------------------------

       integer nlevmax
       parameter (nlevmax=41)
       integer nlev_astex,mxcalc
!       real play(llm), plev_prof(nlevmax)
!       real t_prof(nlevmax),qv_prof(nlevmax)
!       real u_prof(nlevmax),v_prof(nlevmax), w_prof(nlevmax)
!       real ht_prof(nlevmax),vt_prof(nlevmax)
!       real hq_prof(nlevmax),vq_prof(nlevmax)

       real play(llm), plev_prof(nlev_astex)
       real t_prof(nlev_astex),thl_prof(nlev_astex),qv_prof(nlev_astex)
       real u_prof(nlev_astex),v_prof(nlev_astex), w_prof(nlev_astex)
       real o3mmr_prof(nlev_astex),ql_prof(nlev_astex)
       real qt_prof(nlev_astex),tke_prof(nlev_astex)

       real t_mod(llm),thl_mod(llm),qv_mod(llm)
       real u_mod(llm),v_mod(llm), w_mod(llm),tke_mod(llm)
       real o3mmr_mod(llm),ql_mod(llm),qt_mod(llm)

       integer l,k,k1,k2
       real frac,frac1,frac2,fact

       do l = 1, llm

        if (play(l).ge.plev_prof(nlev_astex)) then

        mxcalc=l
         k1=0
         k2=0

         if (play(l).le.plev_prof(1)) then

         do k = 1, nlev_astex-1
          if (play(l).le.plev_prof(k).and. play(l).gt.plev_prof(k+1)) then
            k1=k
            k2=k+1
          endif
         enddo

         if (k1.eq.0 .or. k2.eq.0) then
          write(*,*) 'PB! k1, k2 = ',k1,k2
          write(*,*) 'l,play(l) = ',l,play(l)/100
         do k = 1, nlev_astex-1
          write(*,*) 'k,plev_prof(k) = ',k,plev_prof(k)/100
         enddo
         endif

         frac = (plev_prof(k2)-play(l))/(plev_prof(k2)-plev_prof(k1))
         t_mod(l)= t_prof(k2) - frac*(t_prof(k2)-t_prof(k1))
         thl_mod(l)= thl_prof(k2) - frac*(thl_prof(k2)-thl_prof(k1))
         qv_mod(l)= qv_prof(k2) - frac*(qv_prof(k2)-qv_prof(k1))
         ql_mod(l)= ql_prof(k2) - frac*(ql_prof(k2)-ql_prof(k1))
         qt_mod(l)= qt_prof(k2) - frac*(qt_prof(k2)-qt_prof(k1))
         u_mod(l)= u_prof(k2) - frac*(u_prof(k2)-u_prof(k1))
         v_mod(l)= v_prof(k2) - frac*(v_prof(k2)-v_prof(k1))
         w_mod(l)= w_prof(k2) - frac*(w_prof(k2)-w_prof(k1))
         tke_mod(l)= tke_prof(k2) - frac*(tke_prof(k2)-tke_prof(k1))
         o3mmr_mod(l)=o3mmr_prof(k2)-frac*(o3mmr_prof(k2)-o3mmr_prof(k1))

         else !play>plev_prof(1)

         k1=1
         k2=2
         frac1 = (play(l)-plev_prof(k2))/(plev_prof(k1)-plev_prof(k2))
         frac2 = (play(l)-plev_prof(k1))/(plev_prof(k1)-plev_prof(k2))
         t_mod(l)= frac1*t_prof(k1) - frac2*t_prof(k2)
         thl_mod(l)= frac1*thl_prof(k1) - frac2*thl_prof(k2)
         qv_mod(l)= frac1*qv_prof(k1) - frac2*qv_prof(k2)
         ql_mod(l)= frac1*ql_prof(k1) - frac2*ql_prof(k2)
         qt_mod(l)= frac1*qt_prof(k1) - frac2*qt_prof(k2)
         u_mod(l)= frac1*u_prof(k1) - frac2*u_prof(k2)
         v_mod(l)= frac1*v_prof(k1) - frac2*v_prof(k2)
         w_mod(l)= frac1*w_prof(k1) - frac2*w_prof(k2)
         tke_mod(l)= frac1*tke_prof(k1) - frac2*tke_prof(k2)
         o3mmr_mod(l)= frac1*o3mmr_prof(k1) - frac2*o3mmr_prof(k2)

         endif ! play.le.plev_prof(1)

        else ! above max altitude of forcing file

!jyg
         fact=20.*(plev_prof(nlev_astex)-play(l))/plev_prof(nlev_astex) !jyg
         fact = max(fact,0.)                                           !jyg
         fact = exp(-fact)                                             !jyg
         t_mod(l)= t_prof(nlev_astex)                                   !jyg
         thl_mod(l)= thl_prof(nlev_astex)                                   !jyg
         qv_mod(l)= qv_prof(nlev_astex)*fact                              !jyg
         ql_mod(l)= ql_prof(nlev_astex)*fact                              !jyg
         qt_mod(l)= qt_prof(nlev_astex)*fact                              !jyg
         u_mod(l)= u_prof(nlev_astex)*fact                              !jyg
         v_mod(l)= v_prof(nlev_astex)*fact                              !jyg
         w_mod(l)= w_prof(nlev_astex)*fact                              !jyg
         tke_mod(l)= tke_prof(nlev_astex)*fact                              !jyg
         o3mmr_mod(l)= o3mmr_prof(nlev_astex)*fact                      !jyg

        endif ! play

       enddo ! l

       do l = 1,llm
!      print *,'t_mod(l),thl_mod(l),qv_mod(l),u_mod(l),v_mod(l) ',
!    $        l,t_mod(l),thl_mod(l),qv_mod(l),u_mod(l),v_mod(l)
       enddo

          return
          end

!======================================================================
      SUBROUTINE read_rico(fich_rico,nlev_rico,ps_rico,play                &
     &             ,ts_rico,t_rico,q_rico,u_rico,v_rico,w_rico             &
     &             ,dth_dyn,dqh_dyn)
      implicit none

!-------------------------------------------------------------------------
! Read RICO forcing data 
!-------------------------------------------------------------------------
#include "dimensions.h"


      integer nlev_rico
      real ts_rico,ps_rico
      real t_rico(llm),q_rico(llm)
      real u_rico(llm),v_rico(llm)
      real w_rico(llm)
      real dth_dyn(llm)
      real dqh_dyn(llm)
      

      real play(llm),zlay(llm)
     

      real prico(nlev_rico),zrico(nlev_rico)

      character*80 fich_rico

      integer k,l

      
      print*,fich_rico
      open(21,file=trim(fich_rico),form='formatted')
        do k=1,llm
      zlay(k)=0.
         enddo
      
        read(21,*) ps_rico,ts_rico
        prico(1)=ps_rico
        zrico(1)=0.0
      do l=2,nlev_rico
        read(21,*) k,prico(l),zrico(l)
      enddo
       close(21)

      do k=1,llm
        do l=1,80
          if(prico(l)>play(k)) then
              if(play(k)>prico(l+1)) then
                zlay(k)=zrico(l)+(play(k)-prico(l)) *                      &
     &              (zrico(l+1)-zrico(l))/(prico(l+1)-prico(l))
              else 
                zlay(k)=zrico(l)+(play(k)-prico(80))*                      &
     &              (zrico(81)-zrico(80))/(prico(81)-prico(80))
              endif
          endif
        enddo
        print*,k,zlay(k)
        ! U
        if(0 < zlay(k) .and. zlay(k) < 4000) then
          u_rico(k)=-9.9 + (-1.9 + 9.9)*zlay(k)/4000
        elseif(4000 < zlay(k) .and. zlay(k) < 12000) then
       u_rico(k)=  -1.9 + (30.0 + 1.9) /                                   &
     &          (12000 - 4000) * (zlay(k) - 4000)
        elseif(12000 < zlay(k) .and. zlay(k) < 13000) then
          u_rico(k)=30.0
        elseif(13000 < zlay(k) .and. zlay(k) < 20000) then
          u_rico(k)=30.0 - (30.0) /                                        &
     & (20000 - 13000) * (zlay(k) - 13000)
        else
          u_rico(k)=0.0
        endif

!Q_v
        if(0 < zlay(k) .and. zlay(k) < 740) then
          q_rico(k)=16.0 + (13.8 - 16.0) / (740) * zlay(k)
        elseif(740 < zlay(k) .and. zlay(k) < 3260) then
          q_rico(k)=13.8 + (2.4 - 13.8) /                                   &
     &          (3260 - 740) * (zlay(k) - 740)
        elseif(3260 < zlay(k) .and. zlay(k) < 4000) then
          q_rico(k)=2.4 + (1.8 - 2.4) /                                    &
     &               (4000 - 3260) * (zlay(k) - 3260)
        elseif(4000 < zlay(k) .and. zlay(k) < 9000) then
          q_rico(k)=1.8 + (0 - 1.8) /                                      &
     &             (9000 - 4000) * (zlay(k) - 4000)
        else
          q_rico(k)=0.0
        endif

!T
        if(0 < zlay(k) .and. zlay(k) < 740) then
          t_rico(k)=299.2 + (292.0 - 299.2) / (740) * zlay(k)
        elseif(740 < zlay(k) .and. zlay(k) < 4000) then
          t_rico(k)=292.0 + (278.0 - 292.0) /                              &                       
     &       (4000 - 740) * (zlay(k) - 740)
        elseif(4000 < zlay(k) .and. zlay(k) < 15000) then
          t_rico(k)=278.0 + (203.0 - 278.0) /                              &
     &       (15000 - 4000) * (zlay(k) - 4000)
        elseif(15000 < zlay(k) .and. zlay(k) < 17500) then
          t_rico(k)=203.0 + (194.0 - 203.0) /                              & 
     &       (17500 - 15000)* (zlay(k) - 15000)
        elseif(17500 < zlay(k) .and. zlay(k) < 20000) then
          t_rico(k)=194.0 + (206.0 - 194.0) /                              &
     &       (20000 - 17500)* (zlay(k) - 17500)
        elseif(20000 < zlay(k) .and. zlay(k) < 60000) then
          t_rico(k)=206.0 + (270.0 - 206.0) /                              & 
     &        (60000 - 20000)* (zlay(k) - 20000)
        endif

! W
        if(0 < zlay(k) .and. zlay(k) < 2260 ) then
          w_rico(k)=- (0.005/2260) * zlay(k)
        elseif(2260 < zlay(k) .and. zlay(k) < 4000 ) then
          w_rico(k)=- 0.005
        elseif(4000 < zlay(k) .and. zlay(k) < 5000 ) then
       w_rico(k)=- 0.005 + (0.005/ (5000 - 4000)) * (zlay(k) - 4000)
        else
          w_rico(k)=0.0
        endif

! dThrz+dTsw0+dTlw0
        if(0 < zlay(k) .and. zlay(k) < 4000) then
          dth_dyn(k)=- 2.51 / 86400 + (-2.18 + 2.51 )/                     &
     &               (86400*4000) * zlay(k)
        elseif(4000 < zlay(k) .and. zlay(k) < 5000) then
          dth_dyn(k)=- 2.18 / 86400 + ( 2.18 ) /                           &
     &           (86400*(5000 - 4000)) * (zlay(k) - 4000)
        else
          dth_dyn(k)=0.0
        endif
! dQhrz
        if(0 < zlay(k) .and. zlay(k) < 3000) then
          dqh_dyn(k)=-1.0 / 86400 + (0.345 + 1.0)/                         &
     &                    (86400*3000) * (zlay(k))
        elseif(3000 < zlay(k) .and. zlay(k) < 4000) then
          dqh_dyn(k)=0.345 / 86400
        elseif(4000 < zlay(k) .and. zlay(k) < 5000) then
          dqh_dyn(k)=0.345 / 86400 +                                       &
     &   (-0.345)/(86400 * (5000 - 4000)) * (zlay(k)-4000)
        else
          dqh_dyn(k)=0.0
        endif

!?        if(play(k)>6e4) then
!?          ratqs0(1,k)=ratqsbas*(plev(1)-play(k))/(plev(1)-6e4)
!?        elseif((play(k)>3e4).and.(play(k)<6e4)) then
!?          ratqs0(1,k)=ratqsbas+(ratqshaut-ratqsbas)&
!?                          *(6e4-play(k))/(6e4-3e4)
!?        else
!?          ratqs0(1,k)=ratqshaut
!?        endif

      enddo

      do k=1,llm
      q_rico(k)=q_rico(k)/1e3 
      dqh_dyn(k)=dqh_dyn(k)/1e3
      v_rico(k)=-3.8
      enddo

          return
          end

!======================================================================
        SUBROUTINE interp_sandu_time(day,day1,annee_ref                    &
     &             ,year_ini_sandu,day_ini_sandu,nt_sandu,dt_sandu         &
     &             ,nlev_sandu,ts_sandu,ts_prof)
        implicit none

!---------------------------------------------------------------------------------------
! Time interpolation of a 2D field to the timestep corresponding to day
!
! day: current julian day (e.g. 717538.2)
! day1: first day of the simulation
! nt_sandu: total nb of data in the forcing (e.g. 13 for Sanduref)
! dt_sandu: total time interval (in sec) between 2 forcing data (e.g. 6h for Sanduref)
!---------------------------------------------------------------------------------------
! inputs:
        integer annee_ref
        integer nt_sandu,nlev_sandu
        integer year_ini_sandu
        real day, day1,day_ini_sandu,dt_sandu
        real ts_sandu(nt_sandu)
! outputs:
        real ts_prof
! local:
        integer it_sandu1, it_sandu2
        real timeit,time_sandu1,time_sandu2,frac
! Check that initial day of the simulation consistent with SANDU period:
       if (annee_ref.ne.2006 ) then
        print*,'Pour SANDUREF, annee_ref doit etre 2006 '
        print*,'Changer annee_ref dans run.def'
        stop
       endif
!      if (annee_ref.eq.2006 .and. day1.lt.day_ini_sandu) then
!       print*,'SANDUREF debute le 15 Juillet 2006 (jour julien=196)'
!       print*,'Changer dayref dans run.def'
!       stop
!      endif

! Determine timestep relative to the 1st day of TOGA-COARE:
!       timeit=(day-day1)*86400.
!       if (annee_ref.eq.1992) then
!        timeit=(day-day_ini_sandu)*86400.
!       else
!        timeit=(day+61.-1.)*86400. ! 61 days between Nov01 and Dec31 1992
!       endif
      timeit=(day-day_ini_sandu)*86400

! Determine the closest observation times:
       it_sandu1=INT(timeit/dt_sandu)+1
       it_sandu2=it_sandu1 + 1
       time_sandu1=(it_sandu1-1)*dt_sandu
       time_sandu2=(it_sandu2-1)*dt_sandu
       print *,'timeit day day_ini_sandu',timeit,day,day_ini_sandu
       print *,'it_sandu1,it_sandu2,time_sandu1,time_sandu2',              &
     &          it_sandu1,it_sandu2,time_sandu1,time_sandu2

       if (it_sandu1 .ge. nt_sandu) then
        write(*,*) 'PB-stop: day, it_sandu1, it_sandu2, timeit: '          &
     &        ,day,it_sandu1,it_sandu2,timeit/86400.
        stop
       endif

! time interpolation:
       frac=(time_sandu2-timeit)/(time_sandu2-time_sandu1)
       frac=max(frac,0.0)

       ts_prof = ts_sandu(it_sandu2)                                       &
     &          -frac*(ts_sandu(it_sandu2)-ts_sandu(it_sandu1))

         print*,                                                           &
     &'day,annee_ref,day_ini_sandu,timeit,it_sandu1,it_sandu2,SST:',       &
     &day,annee_ref,day_ini_sandu,timeit/86400.,it_sandu1,                  &
     &it_sandu2,ts_prof

        return
        END
!=====================================================================
!-------------------------------------------------------------------------
      SUBROUTINE read_armcu(fich_armcu,nlev_armcu,nt_armcu,                &
     & sens,flat,adv_theta,rad_theta,adv_qt)
      implicit none

!-------------------------------------------------------------------------
! Read ARM_CU case forcing data
!-------------------------------------------------------------------------

      integer nlev_armcu,nt_armcu
      real sens(nt_armcu),flat(nt_armcu)
      real adv_theta(nt_armcu),rad_theta(nt_armcu),adv_qt(nt_armcu)
      character*80 fich_armcu

      integer ip

      integer iy,im,id,ih,in

      print*,'nlev_armcu',nlev_armcu

      open(21,file=trim(fich_armcu),form='formatted')
      read(21,'(a)')
      do ip = 1, nt_armcu
      read(21,'(a)')
      read(21,'(a)')
      read(21,223) iy, im, id, ih, in, sens(ip),flat(ip),                  &
     &             adv_theta(ip),rad_theta(ip),adv_qt(ip)
      print *,'forcages=',iy,im,id,ih,in, sens(ip),flat(ip),               &
     &             adv_theta(ip),rad_theta(ip),adv_qt(ip)
      enddo
      close(21)

  223 format(5i3,5f8.3)

          return
          end

!=====================================================================
       SUBROUTINE interp_toga_vertical(play,nlev_toga,plev_prof            &
     &         ,t_prof,q_prof,u_prof,v_prof,w_prof                         &
     &         ,ht_prof,vt_prof,hq_prof,vq_prof                            &
     &         ,t_mod,q_mod,u_mod,v_mod,w_mod                              &
     &         ,ht_mod,vt_mod,hq_mod,vq_mod,mxcalc)
 
       implicit none
 
#include "dimensions.h"

!-------------------------------------------------------------------------
! Vertical interpolation of TOGA-COARE forcing data onto model levels
!-------------------------------------------------------------------------
 
       integer nlevmax
       parameter (nlevmax=41)
       integer nlev_toga,mxcalc
!       real play(llm), plev_prof(nlevmax) 
!       real t_prof(nlevmax),q_prof(nlevmax)
!       real u_prof(nlevmax),v_prof(nlevmax), w_prof(nlevmax)
!       real ht_prof(nlevmax),vt_prof(nlevmax)
!       real hq_prof(nlevmax),vq_prof(nlevmax)
 
       real play(llm), plev_prof(nlev_toga) 
       real t_prof(nlev_toga),q_prof(nlev_toga)
       real u_prof(nlev_toga),v_prof(nlev_toga), w_prof(nlev_toga)
       real ht_prof(nlev_toga),vt_prof(nlev_toga)
       real hq_prof(nlev_toga),vq_prof(nlev_toga)
 
       real t_mod(llm),q_mod(llm)
       real u_mod(llm),v_mod(llm), w_mod(llm)
       real ht_mod(llm),vt_mod(llm)
       real hq_mod(llm),vq_mod(llm)
 
       integer l,k,k1,k2
       real frac,frac1,frac2,fact
 
       do l = 1, llm

        if (play(l).ge.plev_prof(nlev_toga)) then
 
        mxcalc=l
         k1=0
         k2=0

         if (play(l).le.plev_prof(1)) then

         do k = 1, nlev_toga-1
          if (play(l).le.plev_prof(k).and. play(l).gt.plev_prof(k+1)) then
            k1=k
            k2=k+1
          endif
         enddo

         if (k1.eq.0 .or. k2.eq.0) then
          write(*,*) 'PB! k1, k2 = ',k1,k2
          write(*,*) 'l,play(l) = ',l,play(l)/100
         do k = 1, nlev_toga-1
          write(*,*) 'k,plev_prof(k) = ',k,plev_prof(k)/100
         enddo
         endif

         frac = (plev_prof(k2)-play(l))/(plev_prof(k2)-plev_prof(k1))
         t_mod(l)= t_prof(k2) - frac*(t_prof(k2)-t_prof(k1))
         q_mod(l)= q_prof(k2) - frac*(q_prof(k2)-q_prof(k1))
         u_mod(l)= u_prof(k2) - frac*(u_prof(k2)-u_prof(k1))
         v_mod(l)= v_prof(k2) - frac*(v_prof(k2)-v_prof(k1))
         w_mod(l)= w_prof(k2) - frac*(w_prof(k2)-w_prof(k1))
         ht_mod(l)= ht_prof(k2) - frac*(ht_prof(k2)-ht_prof(k1))
         vt_mod(l)= vt_prof(k2) - frac*(vt_prof(k2)-vt_prof(k1))
         hq_mod(l)= hq_prof(k2) - frac*(hq_prof(k2)-hq_prof(k1))
         vq_mod(l)= vq_prof(k2) - frac*(vq_prof(k2)-vq_prof(k1))
     
         else !play>plev_prof(1)

         k1=1
         k2=2
         frac1 = (play(l)-plev_prof(k2))/(plev_prof(k1)-plev_prof(k2))
         frac2 = (play(l)-plev_prof(k1))/(plev_prof(k1)-plev_prof(k2))
         t_mod(l)= frac1*t_prof(k1) - frac2*t_prof(k2)
         q_mod(l)= frac1*q_prof(k1) - frac2*q_prof(k2)
         u_mod(l)= frac1*u_prof(k1) - frac2*u_prof(k2)
         v_mod(l)= frac1*v_prof(k1) - frac2*v_prof(k2)
         w_mod(l)= frac1*w_prof(k1) - frac2*w_prof(k2)
         ht_mod(l)= frac1*ht_prof(k1) - frac2*ht_prof(k2)
         vt_mod(l)= frac1*vt_prof(k1) - frac2*vt_prof(k2)
         hq_mod(l)= frac1*hq_prof(k1) - frac2*hq_prof(k2)
         vq_mod(l)= frac1*vq_prof(k1) - frac2*vq_prof(k2)

         endif ! play.le.plev_prof(1)

        else ! above max altitude of forcing file
 
!jyg
         fact=20.*(plev_prof(nlev_toga)-play(l))/plev_prof(nlev_toga) !jyg
         fact = max(fact,0.)                                           !jyg
         fact = exp(-fact)                                             !jyg
         t_mod(l)= t_prof(nlev_toga)                                   !jyg
         q_mod(l)= q_prof(nlev_toga)*fact                              !jyg
         u_mod(l)= u_prof(nlev_toga)*fact                              !jyg
         v_mod(l)= v_prof(nlev_toga)*fact                              !jyg
         w_mod(l)= 0.0                                                 !jyg
         ht_mod(l)= ht_prof(nlev_toga)                                 !jyg
         vt_mod(l)= vt_prof(nlev_toga)                                 !jyg
         hq_mod(l)= hq_prof(nlev_toga)*fact                            !jyg
         vq_mod(l)= vq_prof(nlev_toga)*fact                            !jyg
 
        endif ! play
 
       enddo ! l

!       do l = 1,llm
!       print *,'t_mod(l),q_mod(l),ht_mod(l),hq_mod(l) ',
!     $        l,t_mod(l),q_mod(l),ht_mod(l),hq_mod(l)
!       enddo
 
          return
          end
 
!=====================================================================
       SUBROUTINE interp_case_vertical(play,nlev_cas,plev_prof_cas            &
     &         ,t_prof_cas,q_prof_cas,u_prof_cas,v_prof_cas,ug_prof_cas,vg_prof_cas,vitw_prof_cas                         &
     &         ,du_prof_cas,hu_prof_cas,vu_prof_cas,dv_prof_cas,hv_prof_cas,vv_prof_cas           &
     &         ,dt_prof_cas,ht_prof_cas,vt_prof_cas,dtrad_prof_cas,dq_prof_cas,hq_prof_cas,vq_prof_cas                            &
     &         ,t_mod_cas,q_mod_cas,u_mod_cas,v_mod_cas,ug_mod_cas,vg_mod_cas,w_mod_cas                              &
     &         ,du_mod_cas,hu_mod_cas,vu_mod_cas,dv_mod_cas,hv_mod_cas,vv_mod_cas               &
     &         ,dt_mod_cas,ht_mod_cas,vt_mod_cas,dtrad_mod_cas,dq_mod_cas,hq_mod_cas,vq_mod_cas,mxcalc)
 
       implicit none
 
#include "dimensions.h"

!-------------------------------------------------------------------------
! Vertical interpolation of TOGA-COARE forcing data onto mod_casel levels
!-------------------------------------------------------------------------
 
       integer nlevmax
       parameter (nlevmax=41)
       integer nlev_cas,mxcalc
!       real play(llm), plev_prof(nlevmax) 
!       real t_prof(nlevmax),q_prof(nlevmax)
!       real u_prof(nlevmax),v_prof(nlevmax), w_prof(nlevmax)
!       real ht_prof(nlevmax),vt_prof(nlevmax)
!       real hq_prof(nlevmax),vq_prof(nlevmax)
 
       real play(llm), plev_prof_cas(nlev_cas) 
       real t_prof_cas(nlev_cas),q_prof_cas(nlev_cas)
       real u_prof_cas(nlev_cas),v_prof_cas(nlev_cas)
       real ug_prof_cas(nlev_cas),vg_prof_cas(nlev_cas), vitw_prof_cas(nlev_cas)
       real du_prof_cas(nlev_cas),hu_prof_cas(nlev_cas),vu_prof_cas(nlev_cas)
       real dv_prof_cas(nlev_cas),hv_prof_cas(nlev_cas),vv_prof_cas(nlev_cas)
       real dt_prof_cas(nlev_cas),ht_prof_cas(nlev_cas),vt_prof_cas(nlev_cas),dtrad_prof_cas(nlev_cas)
       real dq_prof_cas(nlev_cas),hq_prof_cas(nlev_cas),vq_prof_cas(nlev_cas)
 
       real t_mod_cas(llm),q_mod_cas(llm)
       real u_mod_cas(llm),v_mod_cas(llm)
       real ug_mod_cas(llm),vg_mod_cas(llm), w_mod_cas(llm)
       real du_mod_cas(llm),hu_mod_cas(llm),vu_mod_cas(llm)
       real dv_mod_cas(llm),hv_mod_cas(llm),vv_mod_cas(llm)
       real dt_mod_cas(llm),ht_mod_cas(llm),vt_mod_cas(llm),dtrad_mod_cas(llm)
       real dq_mod_cas(llm),hq_mod_cas(llm),vq_mod_cas(llm)
 
       integer l,k,k1,k2
       real frac,frac1,frac2,fact
 
       do l = 1, llm

        if (play(l).ge.plev_prof_cas(nlev_cas)) then
 
        mxcalc=l
         k1=0
         k2=0

         if (play(l).le.plev_prof_cas(1)) then

         do k = 1, nlev_cas-1
          if (play(l).le.plev_prof_cas(k).and. play(l).gt.plev_prof_cas(k+1)) then
            k1=k
            k2=k+1
          endif
         enddo

         if (k1.eq.0 .or. k2.eq.0) then
          write(*,*) 'PB! k1, k2 = ',k1,k2
          write(*,*) 'l,play(l) = ',l,play(l)/100
         do k = 1, nlev_cas-1
          write(*,*) 'k,plev_prof_cas(k) = ',k,plev_prof_cas(k)/100
         enddo
         endif

         frac = (plev_prof_cas(k2)-play(l))/(plev_prof_cas(k2)-plev_prof_cas(k1))
         t_mod_cas(l)= t_prof_cas(k2) - frac*(t_prof_cas(k2)-t_prof_cas(k1))
         q_mod_cas(l)= q_prof_cas(k2) - frac*(q_prof_cas(k2)-q_prof_cas(k1))
         u_mod_cas(l)= u_prof_cas(k2) - frac*(u_prof_cas(k2)-u_prof_cas(k1))
         v_mod_cas(l)= v_prof_cas(k2) - frac*(v_prof_cas(k2)-v_prof_cas(k1))
         ug_mod_cas(l)= ug_prof_cas(k2) - frac*(ug_prof_cas(k2)-ug_prof_cas(k1))
         vg_mod_cas(l)= vg_prof_cas(k2) - frac*(vg_prof_cas(k2)-vg_prof_cas(k1))
         w_mod_cas(l)= vitw_prof_cas(k2) - frac*(vitw_prof_cas(k2)-vitw_prof_cas(k1))
         du_mod_cas(l)= du_prof_cas(k2) - frac*(du_prof_cas(k2)-du_prof_cas(k1))
         hu_mod_cas(l)= hu_prof_cas(k2) - frac*(hu_prof_cas(k2)-hu_prof_cas(k1))
         vu_mod_cas(l)= vu_prof_cas(k2) - frac*(vu_prof_cas(k2)-vu_prof_cas(k1))
         dv_mod_cas(l)= dv_prof_cas(k2) - frac*(dv_prof_cas(k2)-dv_prof_cas(k1))
         hv_mod_cas(l)= hv_prof_cas(k2) - frac*(hv_prof_cas(k2)-hv_prof_cas(k1))
         vv_mod_cas(l)= vv_prof_cas(k2) - frac*(vv_prof_cas(k2)-vv_prof_cas(k1))
         dt_mod_cas(l)= dt_prof_cas(k2) - frac*(dt_prof_cas(k2)-dt_prof_cas(k1))
         ht_mod_cas(l)= ht_prof_cas(k2) - frac*(ht_prof_cas(k2)-ht_prof_cas(k1))
         vt_mod_cas(l)= vt_prof_cas(k2) - frac*(vt_prof_cas(k2)-vt_prof_cas(k1))
         dq_mod_cas(l)= dq_prof_cas(k2) - frac*(dq_prof_cas(k2)-dq_prof_cas(k1))
         hq_mod_cas(l)= hq_prof_cas(k2) - frac*(hq_prof_cas(k2)-hq_prof_cas(k1))
         vq_mod_cas(l)= vq_prof_cas(k2) - frac*(vq_prof_cas(k2)-vq_prof_cas(k1))
         dtrad_mod_cas(l)= dtrad_prof_cas(k2) - frac*(dtrad_prof_cas(k2)-dtrad_prof_cas(k1))
     
         else !play>plev_prof_cas(1)

         k1=1
         k2=2
         frac1 = (play(l)-plev_prof_cas(k2))/(plev_prof_cas(k1)-plev_prof_cas(k2))
         frac2 = (play(l)-plev_prof_cas(k1))/(plev_prof_cas(k1)-plev_prof_cas(k2))
         t_mod_cas(l)= frac1*t_prof_cas(k1) - frac2*t_prof_cas(k2)
         q_mod_cas(l)= frac1*q_prof_cas(k1) - frac2*q_prof_cas(k2)
         u_mod_cas(l)= frac1*u_prof_cas(k1) - frac2*u_prof_cas(k2)
         v_mod_cas(l)= frac1*v_prof_cas(k1) - frac2*v_prof_cas(k2)
         ug_mod_cas(l)= frac1*ug_prof_cas(k1) - frac2*ug_prof_cas(k2)
         vg_mod_cas(l)= frac1*vg_prof_cas(k1) - frac2*vg_prof_cas(k2)
         w_mod_cas(l)= frac1*vitw_prof_cas(k1) - frac2*vitw_prof_cas(k2)
         du_mod_cas(l)= frac1*du_prof_cas(k1) - frac2*du_prof_cas(k2)
         hu_mod_cas(l)= frac1*hu_prof_cas(k1) - frac2*hu_prof_cas(k2)
         vu_mod_cas(l)= frac1*vu_prof_cas(k1) - frac2*vu_prof_cas(k2)
         dv_mod_cas(l)= frac1*dv_prof_cas(k1) - frac2*dv_prof_cas(k2)
         hv_mod_cas(l)= frac1*hv_prof_cas(k1) - frac2*hv_prof_cas(k2)
         vv_mod_cas(l)= frac1*vv_prof_cas(k1) - frac2*vv_prof_cas(k2)
         dt_mod_cas(l)= frac1*dt_prof_cas(k1) - frac2*dt_prof_cas(k2)
         ht_mod_cas(l)= frac1*ht_prof_cas(k1) - frac2*ht_prof_cas(k2)
         vt_mod_cas(l)= frac1*vt_prof_cas(k1) - frac2*vt_prof_cas(k2)
         dq_mod_cas(l)= frac1*dq_prof_cas(k1) - frac2*dq_prof_cas(k2)
         hq_mod_cas(l)= frac1*hq_prof_cas(k1) - frac2*hq_prof_cas(k2)
         vq_mod_cas(l)= frac1*vq_prof_cas(k1) - frac2*vq_prof_cas(k2)
         dtrad_mod_cas(l)= frac1*dtrad_prof_cas(k1) - frac2*dtrad_prof_cas(k2)

         endif ! play.le.plev_prof_cas(1)

        else ! above max altitude of forcing file
 
!jyg
         fact=20.*(plev_prof_cas(nlev_cas)-play(l))/plev_prof_cas(nlev_cas) !jyg
         fact = max(fact,0.)                                           !jyg
         fact = exp(-fact)                                             !jyg
         t_mod_cas(l)= t_prof_cas(nlev_cas)                                   !jyg
         q_mod_cas(l)= q_prof_cas(nlev_cas)*fact                              !jyg
         u_mod_cas(l)= u_prof_cas(nlev_cas)*fact                              !jyg
         v_mod_cas(l)= v_prof_cas(nlev_cas)*fact                              !jyg
         ug_mod_cas(l)= ug_prof_cas(nlev_cas)*fact                              !jyg
         vg_mod_cas(l)= vg_prof_cas(nlev_cas)*fact                              !jyg
         w_mod_cas(l)= 0.0                                                 !jyg
         du_mod_cas(l)= du_prof_cas(nlev_cas)*fact 
         hu_mod_cas(l)= hu_prof_cas(nlev_cas)*fact                            !jyg
         vu_mod_cas(l)= vu_prof_cas(nlev_cas)*fact                            !jyg
         dv_mod_cas(l)= dv_prof_cas(nlev_cas)*fact 
         hv_mod_cas(l)= hv_prof_cas(nlev_cas)*fact                            !jyg
         vv_mod_cas(l)= vv_prof_cas(nlev_cas)*fact                            !jyg
         dt_mod_cas(l)= dt_prof_cas(nlev_cas) 
         ht_mod_cas(l)= ht_prof_cas(nlev_cas)                                 !jyg
         vt_mod_cas(l)= vt_prof_cas(nlev_cas)                                 !jyg
         dq_mod_cas(l)= dq_prof_cas(nlev_cas)*fact 
         hq_mod_cas(l)= hq_prof_cas(nlev_cas)*fact                            !jyg
         vq_mod_cas(l)= vq_prof_cas(nlev_cas)*fact                            !jyg
         dtrad_mod_cas(l)= dtrad_prof_cas(nlev_cas)*fact                      !jyg
 
        endif ! play
 
       enddo ! l

!       do l = 1,llm
!       print *,'t_mod_cas(l),q_mod_cas(l),ht_mod_cas(l),hq_mod_cas(l) ',
!     $        l,t_mod_cas(l),q_mod_cas(l),ht_mod_cas(l),hq_mod_cas(l)
!       enddo
 
          return
          end
!***************************************************************************** 
!=====================================================================
       SUBROUTINE interp_dice_vertical(play,nlev_dice,nt_dice,plev_prof   &
     &         ,th_prof,qv_prof,u_prof,v_prof,o3_prof                     &
     &         ,ht_prof,hq_prof,hu_prof,hv_prof,w_prof,omega_prof         &
     &         ,th_mod,qv_mod,u_mod,v_mod,o3_mod                          &
     &         ,ht_mod,hq_mod,hu_mod,hv_mod,w_mod,omega_mod,mxcalc)
 
       implicit none
 
#include "dimensions.h"

!-------------------------------------------------------------------------
! Vertical interpolation of Dice forcing data onto model levels
!-------------------------------------------------------------------------
 
       integer nlevmax
       parameter (nlevmax=41)
       integer nlev_dice,mxcalc,nt_dice
 
       real play(llm), plev_prof(nlev_dice) 
       real th_prof(nlev_dice),qv_prof(nlev_dice)
       real u_prof(nlev_dice),v_prof(nlev_dice) 
       real o3_prof(nlev_dice)
       real ht_prof(nlev_dice),hq_prof(nlev_dice)
       real hu_prof(nlev_dice),hv_prof(nlev_dice)
       real w_prof(nlev_dice),omega_prof(nlev_dice)
 
       real th_mod(llm),qv_mod(llm)
       real u_mod(llm),v_mod(llm), o3_mod(llm)
       real ht_mod(llm),hq_mod(llm)
       real hu_mod(llm),hv_mod(llm),w_mod(llm),omega_mod(llm)
 
       integer l,k,k1,k2,kp
       real aa,frac,frac1,frac2,fact
 
       do l = 1, llm

        if (play(l).ge.plev_prof(nlev_dice)) then
 
        mxcalc=l
         k1=0
         k2=0

         if (play(l).le.plev_prof(1)) then

         do k = 1, nlev_dice-1
          if (play(l).le.plev_prof(k) .and. play(l).gt.plev_prof(k+1)) then
            k1=k
            k2=k+1
          endif
         enddo

         if (k1.eq.0 .or. k2.eq.0) then
          write(*,*) 'PB! k1, k2 = ',k1,k2
          write(*,*) 'l,play(l) = ',l,play(l)/100
         do k = 1, nlev_dice-1
          write(*,*) 'k,plev_prof(k) = ',k,plev_prof(k)/100
         enddo
         endif

         frac = (plev_prof(k2)-play(l))/(plev_prof(k2)-plev_prof(k1))
         th_mod(l)= th_prof(k2) - frac*(th_prof(k2)-th_prof(k1))
         qv_mod(l)= qv_prof(k2) - frac*(qv_prof(k2)-qv_prof(k1))
         u_mod(l)= u_prof(k2) - frac*(u_prof(k2)-u_prof(k1))
         v_mod(l)= v_prof(k2) - frac*(v_prof(k2)-v_prof(k1))
         o3_mod(l)= o3_prof(k2) - frac*(o3_prof(k2)-o3_prof(k1))
         ht_mod(l)= ht_prof(k2) - frac*(ht_prof(k2)-ht_prof(k1))
         hq_mod(l)= hq_prof(k2) - frac*(hq_prof(k2)-hq_prof(k1))
         hu_mod(l)= hu_prof(k2) - frac*(hu_prof(k2)-hu_prof(k1))
         hv_mod(l)= hv_prof(k2) - frac*(hv_prof(k2)-hv_prof(k1))
         w_mod(l)= w_prof(k2) - frac*(w_prof(k2)-w_prof(k1))
         omega_mod(l)= omega_prof(k2) - frac*(omega_prof(k2)-omega_prof(k1))
     
         else !play>plev_prof(1)

         k1=1
         k2=2
         frac1 = (play(l)-plev_prof(k2))/(plev_prof(k1)-plev_prof(k2))
         frac2 = (play(l)-plev_prof(k1))/(plev_prof(k1)-plev_prof(k2))
         th_mod(l)= frac1*th_prof(k1) - frac2*th_prof(k2)
         qv_mod(l)= frac1*qv_prof(k1) - frac2*qv_prof(k2)
         u_mod(l)= frac1*u_prof(k1) - frac2*u_prof(k2)
         v_mod(l)= frac1*v_prof(k1) - frac2*v_prof(k2)
         o3_mod(l)= frac1*o3_prof(k1) - frac2*o3_prof(k2)
         ht_mod(l)= frac1*ht_prof(k1) - frac2*ht_prof(k2)
         hq_mod(l)= frac1*hq_prof(k1) - frac2*hq_prof(k2)
         hu_mod(l)= frac1*hu_prof(k1) - frac2*hu_prof(k2)
         hv_mod(l)= frac1*hv_prof(k1) - frac2*hv_prof(k2)
         w_mod(l)= frac1*w_prof(k1) - frac2*w_prof(k2)
         omega_mod(l)= frac1*omega_prof(k1) - frac2*omega_prof(k2)

         endif ! play.le.plev_prof(1)

        else ! above max altitude of forcing file
 
!jyg
         fact=20.*(plev_prof(nlev_dice)-play(l))/plev_prof(nlev_dice) !jyg
         fact = max(fact,0.)                                           !jyg
         fact = exp(-fact)                                             !jyg
         th_mod(l)= th_prof(nlev_dice)                                 !jyg
         qv_mod(l)= qv_prof(nlev_dice)*fact                            !jyg
         u_mod(l)= u_prof(nlev_dice)*fact                              !jyg
         v_mod(l)= v_prof(nlev_dice)*fact                              !jyg
         o3_mod(l)= o3_prof(nlev_dice)*fact                            !jyg
         ht_mod(l)= ht_prof(nlev_dice)                                 !jyg
         hq_mod(l)= hq_prof(nlev_dice)*fact                            !jyg
         hu_mod(l)= hu_prof(nlev_dice)                                 !jyg
         hv_mod(l)= hv_prof(nlev_dice)                                 !jyg
         w_mod(l)= 0.                                                  !jyg
         omega_mod(l)= 0.                                              !jyg
 
        endif ! play
 
       enddo ! l

!       do l = 1,llm
!       print *,'t_mod(l),q_mod(l),ht_mod(l),hq_mod(l) ',
!     $        l,t_mod(l),q_mod(l),ht_mod(l),hq_mod(l)
!       enddo
 
          return
          end

!======================================================================
        SUBROUTINE interp_astex_time(day,day1,annee_ref                    &
     &             ,year_ini_astex,day_ini_astex,nt_astex,dt_astex         &
     &             ,nlev_astex,div_astex,ts_astex,ug_astex,vg_astex        &
     &             ,ufa_astex,vfa_astex,div_prof,ts_prof,ug_prof,vg_prof   &
     &             ,ufa_prof,vfa_prof)
        implicit none

!---------------------------------------------------------------------------------------
! Time interpolation of a 2D field to the timestep corresponding to day
!
! day: current julian day (e.g. 717538.2)
! day1: first day of the simulation
! nt_astex: total nb of data in the forcing (e.g. 41 for Astex)
! dt_astex: total time interval (in sec) between 2 forcing data (e.g. 1h for Astex)
!---------------------------------------------------------------------------------------

! inputs:
        integer annee_ref
        integer nt_astex,nlev_astex
        integer year_ini_astex
        real day, day1,day_ini_astex,dt_astex
        real div_astex(nt_astex),ts_astex(nt_astex),ug_astex(nt_astex)
        real vg_astex(nt_astex),ufa_astex(nt_astex),vfa_astex(nt_astex)
! outputs:
        real div_prof,ts_prof,ug_prof,vg_prof,ufa_prof,vfa_prof
! local:
        integer it_astex1, it_astex2
        real timeit,time_astex1,time_astex2,frac

! Check that initial day of the simulation consistent with ASTEX period:
       if (annee_ref.ne.1992 ) then
        print*,'Pour Astex, annee_ref doit etre 1992 '
        print*,'Changer annee_ref dans run.def'
        stop
       endif
       if (annee_ref.eq.1992 .and. day1.lt.day_ini_astex) then
        print*,'Astex debute le 13 Juin 1992 (jour julien=165)'
        print*,'Changer dayref dans run.def'
        stop
       endif

! Determine timestep relative to the 1st day of TOGA-COARE:
!       timeit=(day-day1)*86400.
!       if (annee_ref.eq.1992) then
!        timeit=(day-day_ini_astex)*86400.
!       else
!        timeit=(day+2.-1.)*86400. ! 2 days between Jun13 and Jun15 1992
!       endif
      timeit=(day-day_ini_astex)*86400

! Determine the closest observation times:
       it_astex1=INT(timeit/dt_astex)+1
       it_astex2=it_astex1 + 1
       time_astex1=(it_astex1-1)*dt_astex
       time_astex2=(it_astex2-1)*dt_astex
       print *,'timeit day day_ini_astex',timeit,day,day_ini_astex
       print *,'it_astex1,it_astex2,time_astex1,time_astex2',              &
     &          it_astex1,it_astex2,time_astex1,time_astex2

       if (it_astex1 .ge. nt_astex) then
        write(*,*) 'PB-stop: day, it_astex1, it_astex2, timeit: '          &
     &        ,day,it_astex1,it_astex2,timeit/86400.
        stop
       endif

! time interpolation:
       frac=(time_astex2-timeit)/(time_astex2-time_astex1)
       frac=max(frac,0.0)

       div_prof = div_astex(it_astex2)                                     &
     &          -frac*(div_astex(it_astex2)-div_astex(it_astex1))
       ts_prof = ts_astex(it_astex2)                                        &
     &          -frac*(ts_astex(it_astex2)-ts_astex(it_astex1))
       ug_prof = ug_astex(it_astex2)                                       &
     &          -frac*(ug_astex(it_astex2)-ug_astex(it_astex1))
       vg_prof = vg_astex(it_astex2)                                       &
     &          -frac*(vg_astex(it_astex2)-vg_astex(it_astex1))
       ufa_prof = ufa_astex(it_astex2)                                     &
     &          -frac*(ufa_astex(it_astex2)-ufa_astex(it_astex1))
       vfa_prof = vfa_astex(it_astex2)                                     &
     &          -frac*(vfa_astex(it_astex2)-vfa_astex(it_astex1))

         print*,                                                           &
     &'day,annee_ref,day_ini_astex,timeit,it_astex1,it_astex2,SST:',       &
     &day,annee_ref,day_ini_astex,timeit/86400.,it_astex1,                 &
     &it_astex2,div_prof,ts_prof,ug_prof,vg_prof,ufa_prof,vfa_prof 

        return
        END

!======================================================================
        SUBROUTINE interp_toga_time(day,day1,annee_ref                     &
     &             ,year_ini_toga,day_ini_toga,nt_toga,dt_toga,nlev_toga   &
     &             ,ts_toga,plev_toga,t_toga,q_toga,u_toga,v_toga,w_toga   &
     &             ,ht_toga,vt_toga,hq_toga,vq_toga                        &
     &             ,ts_prof,plev_prof,t_prof,q_prof,u_prof,v_prof,w_prof   &
     &             ,ht_prof,vt_prof,hq_prof,vq_prof)
        implicit none

!---------------------------------------------------------------------------------------
! Time interpolation of a 2D field to the timestep corresponding to day
!
! day: current julian day (e.g. 717538.2)
! day1: first day of the simulation
! nt_toga: total nb of data in the forcing (e.g. 480 for TOGA-COARE)
! dt_toga: total time interval (in sec) between 2 forcing data (e.g. 6h for TOGA-COARE)
!---------------------------------------------------------------------------------------

#include "compar1d.h"

! inputs:
        integer annee_ref
        integer nt_toga,nlev_toga
        integer year_ini_toga
        real day, day1,day_ini_toga,dt_toga
        real ts_toga(nt_toga)
        real plev_toga(nlev_toga,nt_toga),t_toga(nlev_toga,nt_toga)
        real q_toga(nlev_toga,nt_toga),u_toga(nlev_toga,nt_toga)
        real v_toga(nlev_toga,nt_toga),w_toga(nlev_toga,nt_toga)
        real ht_toga(nlev_toga,nt_toga),vt_toga(nlev_toga,nt_toga)
        real hq_toga(nlev_toga,nt_toga),vq_toga(nlev_toga,nt_toga)
! outputs:
        real ts_prof
        real plev_prof(nlev_toga),t_prof(nlev_toga)
        real q_prof(nlev_toga),u_prof(nlev_toga)
        real v_prof(nlev_toga),w_prof(nlev_toga)
        real ht_prof(nlev_toga),vt_prof(nlev_toga)
        real hq_prof(nlev_toga),vq_prof(nlev_toga)
! local:
        integer it_toga1, it_toga2,k
        real timeit,time_toga1,time_toga2,frac


        if (forcing_type.eq.2) then
! Check that initial day of the simulation consistent with TOGA-COARE period:
       if (annee_ref.ne.1992 .and. annee_ref.ne.1993) then
        print*,'Pour TOGA-COARE, annee_ref doit etre 1992 ou 1993'
        print*,'Changer annee_ref dans run.def'
        stop
       endif
       if (annee_ref.eq.1992 .and. day1.lt.day_ini_toga) then
        print*,'TOGA-COARE a debute le 1er Nov 1992 (jour julien=306)'
        print*,'Changer dayref dans run.def'
        stop
       endif
       if (annee_ref.eq.1993 .and. day1.gt.day_ini_toga+119) then
        print*,'TOGA-COARE a fini le 28 Feb 1993 (jour julien=59)'
        print*,'Changer dayref ou nday dans run.def'
        stop
       endif

       else if (forcing_type.eq.4) then

! Check that initial day of the simulation consistent with TWP-ICE period:
       if (annee_ref.ne.2006) then
        print*,'Pour TWP-ICE, annee_ref doit etre 2006'
        print*,'Changer annee_ref dans run.def'
        stop
       endif
       if (annee_ref.eq.2006 .and. day1.lt.day_ini_toga) then
        print*,'TWP-ICE a debute le 17 Jan 2006 (jour julien=17)'
        print*,'Changer dayref dans run.def'
        stop
       endif
       if (annee_ref.eq.2006 .and. day1.gt.day_ini_toga+26) then
        print*,'TWP-ICE a fini le 12 Feb 2006 (jour julien=43)'
        print*,'Changer dayref ou nday dans run.def'
        stop
       endif

       endif

! Determine timestep relative to the 1st day of TOGA-COARE:
!       timeit=(day-day1)*86400.
!       if (annee_ref.eq.1992) then
!        timeit=(day-day_ini_toga)*86400.
!       else
!        timeit=(day+61.-1.)*86400. ! 61 days between Nov01 and Dec31 1992
!       endif
      timeit=(day-day_ini_toga)*86400

! Determine the closest observation times:
       it_toga1=INT(timeit/dt_toga)+1
       it_toga2=it_toga1 + 1
       time_toga1=(it_toga1-1)*dt_toga
       time_toga2=(it_toga2-1)*dt_toga

       if (it_toga1 .ge. nt_toga) then
        write(*,*) 'PB-stop: day, it_toga1, it_toga2, timeit: '            &
     &        ,day,it_toga1,it_toga2,timeit/86400.
        stop
       endif

! time interpolation:
       frac=(time_toga2-timeit)/(time_toga2-time_toga1)
       frac=max(frac,0.0)

       ts_prof = ts_toga(it_toga2)                                         &
     &          -frac*(ts_toga(it_toga2)-ts_toga(it_toga1))

!        print*,
!     :'day,annee_ref,day_ini_toga,timeit,it_toga1,it_toga2,SST:',
!     :day,annee_ref,day_ini_toga,timeit/86400.,it_toga1,it_toga2,ts_prof

       do k=1,nlev_toga
        plev_prof(k) = 100.*(plev_toga(k,it_toga2)                         &
     &          -frac*(plev_toga(k,it_toga2)-plev_toga(k,it_toga1)))
        t_prof(k) = t_toga(k,it_toga2)                                     &
     &          -frac*(t_toga(k,it_toga2)-t_toga(k,it_toga1))
        q_prof(k) = q_toga(k,it_toga2)                                     &
     &          -frac*(q_toga(k,it_toga2)-q_toga(k,it_toga1))
        u_prof(k) = u_toga(k,it_toga2)                                     &
     &          -frac*(u_toga(k,it_toga2)-u_toga(k,it_toga1))
        v_prof(k) = v_toga(k,it_toga2)                                     &
     &          -frac*(v_toga(k,it_toga2)-v_toga(k,it_toga1))
        w_prof(k) = w_toga(k,it_toga2)                                     &
     &          -frac*(w_toga(k,it_toga2)-w_toga(k,it_toga1))
        ht_prof(k) = ht_toga(k,it_toga2)                                   &
     &          -frac*(ht_toga(k,it_toga2)-ht_toga(k,it_toga1))
        vt_prof(k) = vt_toga(k,it_toga2)                                   &
     &          -frac*(vt_toga(k,it_toga2)-vt_toga(k,it_toga1))
        hq_prof(k) = hq_toga(k,it_toga2)                                   &
     &          -frac*(hq_toga(k,it_toga2)-hq_toga(k,it_toga1))
        vq_prof(k) = vq_toga(k,it_toga2)                                   &
     &          -frac*(vq_toga(k,it_toga2)-vq_toga(k,it_toga1))
        enddo

        return
        END

!======================================================================
        SUBROUTINE interp_dice_time(day,day1,annee_ref                    &
     &             ,year_ini_dice,day_ini_dice,nt_dice,dt_dice            &
     &             ,nlev_dice,shf_dice,lhf_dice,lwup_dice,swup_dice       &
     &             ,tg_dice,ustar_dice,psurf_dice,ug_dice,vg_dice         &
     &             ,ht_dice,hq_dice,hu_dice,hv_dice,w_dice,omega_dice     &
     &             ,shf_prof,lhf_prof,lwup_prof,swup_prof,tg_prof         &
     &             ,ustar_prof,psurf_prof,ug_prof,vg_prof                 &
     &             ,ht_prof,hq_prof,hu_prof,hv_prof,w_prof,omega_prof)
        implicit none

!---------------------------------------------------------------------------------------
! Time interpolation of a 2D field to the timestep corresponding to day
!
! day: current julian day (e.g. 717538.2)
! day1: first day of the simulation
! nt_dice: total nb of data in the forcing (e.g. 145 for Dice)
! dt_dice: total time interval (in sec) between 2 forcing data (e.g. 30min. for Dice)
!---------------------------------------------------------------------------------------

#include "compar1d.h"

! inputs:
        integer annee_ref
        integer nt_dice,nlev_dice
        integer year_ini_dice
        real day, day1,day_ini_dice,dt_dice
        real shf_dice(nt_dice),lhf_dice(nt_dice),lwup_dice(nt_dice)
        real swup_dice(nt_dice),tg_dice(nt_dice),ustar_dice(nt_dice)
        real psurf_dice(nt_dice),ug_dice(nt_dice),vg_dice(nt_dice)
        real ht_dice(nlev_dice,nt_dice),hq_dice(nlev_dice,nt_dice)
        real hu_dice(nlev_dice,nt_dice),hv_dice(nlev_dice,nt_dice)
        real w_dice(nlev_dice,nt_dice),omega_dice(nlev_dice,nt_dice)
! outputs:
        real tg_prof,shf_prof,lhf_prof,lwup_prof,swup_prof
        real ustar_prof,psurf_prof,ug_prof,vg_prof
        real ht_prof(nlev_dice),hq_prof(nlev_dice)
        real hu_prof(nlev_dice),hv_prof(nlev_dice)
        real w_prof(nlev_dice),omega_prof(nlev_dice)
! local:
        integer it_dice1, it_dice2,k
        real timeit,time_dice1,time_dice2,frac


        if (forcing_type.eq.7) then
! Check that initial day of the simulation consistent with Dice period:
       print *,'annee_ref=',annee_ref
       print *,'day1=',day1
       print *,'day_ini_dice=',day_ini_dice
       if (annee_ref.ne.1999) then
        print*,'Pour Dice, annee_ref doit etre 1999'
        print*,'Changer annee_ref dans run.def'
        stop
       endif
       if (annee_ref.eq.1999 .and. day1.gt.day_ini_dice) then
        print*,'Dice a debute le 23 Oct 1999 (jour julien=296)'
        print*,'Changer dayref dans run.def',day1,day_ini_dice
        stop
       endif
       if (annee_ref.eq.1999 .and. day1.gt.day_ini_dice+2) then
        print*,'Dice a fini le 25 Oct 1999 (jour julien=298)'
        print*,'Changer dayref ou nday dans run.def',day1,day_ini_dice
        stop
       endif

       endif

! Determine timestep relative to the 1st day of TOGA-COARE:
!       timeit=(day-day1)*86400.
!       if (annee_ref.eq.1992) then
!        timeit=(day-day_ini_dice)*86400.
!       else
!        timeit=(day+61.-1.)*86400. ! 61 days between Nov01 and Dec31 1992
!       endif
      timeit=(day-day_ini_dice)*86400

! Determine the closest observation times:
       it_dice1=INT(timeit/dt_dice)+1
       it_dice2=it_dice1 + 1
       time_dice1=(it_dice1-1)*dt_dice
       time_dice2=(it_dice2-1)*dt_dice

       if (it_dice1 .ge. nt_dice) then
        write(*,*) 'PB-stop: day, it_dice1, it_dice2, timeit: ',day,it_dice1,it_dice2,timeit/86400.
        stop
       endif

! time interpolation:
       frac=(time_dice2-timeit)/(time_dice2-time_dice1)
       frac=max(frac,0.0)

       shf_prof = shf_dice(it_dice2)-frac*(shf_dice(it_dice2)-shf_dice(it_dice1))
       lhf_prof = lhf_dice(it_dice2)-frac*(lhf_dice(it_dice2)-lhf_dice(it_dice1))
       lwup_prof = lwup_dice(it_dice2)-frac*(lwup_dice(it_dice2)-lwup_dice(it_dice1))
       swup_prof = swup_dice(it_dice2)-frac*(swup_dice(it_dice2)-swup_dice(it_dice1))
       tg_prof = tg_dice(it_dice2)-frac*(tg_dice(it_dice2)-tg_dice(it_dice1))
       ustar_prof = ustar_dice(it_dice2)-frac*(ustar_dice(it_dice2)-ustar_dice(it_dice1))
       psurf_prof = psurf_dice(it_dice2)-frac*(psurf_dice(it_dice2)-psurf_dice(it_dice1))
       ug_prof = ug_dice(it_dice2)-frac*(ug_dice(it_dice2)-ug_dice(it_dice1))
       vg_prof = vg_dice(it_dice2)-frac*(vg_dice(it_dice2)-vg_dice(it_dice1))

!        print*,
!     :'day,annee_ref,day_ini_dice,timeit,it_dice1,it_dice2,SST:',
!     :day,annee_ref,day_ini_dice,timeit/86400.,it_dice1,it_dice2,ts_prof

       do k=1,nlev_dice
        ht_prof(k) = ht_dice(k,it_dice2)-frac*(ht_dice(k,it_dice2)-ht_dice(k,it_dice1))
        hq_prof(k) = hq_dice(k,it_dice2)-frac*(hq_dice(k,it_dice2)-hq_dice(k,it_dice1))
        hu_prof(k) = hu_dice(k,it_dice2)-frac*(hu_dice(k,it_dice2)-hu_dice(k,it_dice1))
        hv_prof(k) = hv_dice(k,it_dice2)-frac*(hv_dice(k,it_dice2)-hv_dice(k,it_dice1))
        w_prof(k) = w_dice(k,it_dice2)-frac*(w_dice(k,it_dice2)-w_dice(k,it_dice1))
        omega_prof(k) = omega_dice(k,it_dice2)-frac*(omega_dice(k,it_dice2)-omega_dice(k,it_dice1))
        enddo

        return
        END

!======================================================================
        SUBROUTINE interp_gabls4_time(day,day1,annee_ref                              &
     &             ,year_ini_gabls4,day_ini_gabls4,nt_gabls4,dt_gabls4,nlev_gabls4    &
     &             ,ug_gabls4,vg_gabls4,ht_gabls4,hq_gabls4,tg_gabls4                          &
     &             ,ug_prof,vg_prof,ht_prof,hq_prof,tg_prof)
        implicit none

!---------------------------------------------------------------------------------------
! Time interpolation of a 2D field to the timestep corresponding to day
!
! day: current julian day 
! day1: first day of the simulation
! nt_gabls4: total nb of data in the forcing (e.g. 37 for gabls4)
! dt_gabls4: total time interval (in sec) between 2 forcing data (e.g. 60min. for gabls4)
!---------------------------------------------------------------------------------------

#include "compar1d.h"

! inputs:
        integer annee_ref
        integer nt_gabls4,nlev_gabls4
        integer year_ini_gabls4
        real day, day1,day_ini_gabls4,dt_gabls4
        real ug_gabls4(nlev_gabls4,nt_gabls4),vg_gabls4(nlev_gabls4,nt_gabls4)
        real ht_gabls4(nlev_gabls4,nt_gabls4),hq_gabls4(nlev_gabls4,nt_gabls4)
        real tg_gabls4(nt_gabls4), tg_prof
! outputs:
        real ug_prof(nlev_gabls4),vg_prof(nlev_gabls4)
        real ht_prof(nlev_gabls4),hq_prof(nlev_gabls4)
! local:
        integer it_gabls41, it_gabls42,k
        real timeit,time_gabls41,time_gabls42,frac



! Check that initial day of the simulation consistent with gabls4 period:
       if (forcing_type.eq.8 ) then
       print *,'annee_ref=',annee_ref
       print *,'day1=',day1
       print *,'day_ini_gabls4=',day_ini_gabls4
       if (annee_ref.ne.2009) then
        print*,'Pour gabls4, annee_ref doit etre 2009'
        print*,'Changer annee_ref dans run.def'
        stop
       endif
       if (annee_ref.eq.2009 .and. day1.gt.day_ini_gabls4) then
        print*,'gabls4 a debute le 11 dec 2009 (jour julien=345)'
        print*,'Changer dayref dans run.def',day1,day_ini_gabls4
        stop
       endif
       if (annee_ref.eq.2009 .and. day1.gt.day_ini_gabls4+2) then
        print*,'gabls4 a fini le 12 dec 2009 (jour julien=346)'
        print*,'Changer dayref ou nday dans run.def',day1,day_ini_gabls4
        stop
       endif
       endif

      timeit=(day-day_ini_gabls4)*86400
       print *,'day,day_ini_gabls4=',day,day_ini_gabls4
       print *,'nt_gabls4,dt,timeit=',nt_gabls4,dt_gabls4,timeit

! Determine the closest observation times:
       it_gabls41=INT(timeit/dt_gabls4)+1
       it_gabls42=it_gabls41 + 1
       time_gabls41=(it_gabls41-1)*dt_gabls4
       time_gabls42=(it_gabls42-1)*dt_gabls4

       if (it_gabls41 .ge. nt_gabls4) then
        write(*,*) 'PB-stop: day, it_gabls41, it_gabls42, timeit: ',day,it_gabls41,it_gabls42,timeit/86400.
        stop
       endif

! time interpolation:
       frac=(time_gabls42-timeit)/(time_gabls42-time_gabls41)
       frac=max(frac,0.0)


       do k=1,nlev_gabls4
        ug_prof(k) = ug_gabls4(k,it_gabls42)-frac*(ug_gabls4(k,it_gabls42)-ug_gabls4(k,it_gabls41))
        vg_prof(k) = vg_gabls4(k,it_gabls42)-frac*(vg_gabls4(k,it_gabls42)-vg_gabls4(k,it_gabls41))
        ht_prof(k) = ht_gabls4(k,it_gabls42)-frac*(ht_gabls4(k,it_gabls42)-ht_gabls4(k,it_gabls41))
        hq_prof(k) = hq_gabls4(k,it_gabls42)-frac*(hq_gabls4(k,it_gabls42)-hq_gabls4(k,it_gabls41))
        enddo
        tg_prof=tg_gabls4(it_gabls42)-frac*(tg_gabls4(it_gabls42)-tg_gabls4(it_gabls41))
        return
        END

!======================================================================
        SUBROUTINE interp_armcu_time(day,day1,annee_ref                    &
     &             ,year_ini_armcu,day_ini_armcu,nt_armcu,dt_armcu         &
     &             ,nlev_armcu,fs_armcu,fl_armcu,at_armcu,rt_armcu         &
     &             ,aqt_armcu,fs_prof,fl_prof,at_prof,rt_prof,aqt_prof)
        implicit none

!---------------------------------------------------------------------------------------
! Time interpolation of a 2D field to the timestep corresponding to day
!
! day: current julian day (e.g. 717538.2)
! day1: first day of the simulation
! nt_armcu: total nb of data in the forcing (e.g. 31 for armcu)
! dt_armcu: total time interval (in sec) between 2 forcing data (e.g. 1/2h for armcu)
! fs= sensible flux
! fl= latent flux
! at,rt,aqt= advective and radiative tendencies
!---------------------------------------------------------------------------------------

! inputs:
        integer annee_ref
        integer nt_armcu,nlev_armcu
        integer year_ini_armcu
        real day, day1,day_ini_armcu,dt_armcu
        real fs_armcu(nt_armcu),fl_armcu(nt_armcu),at_armcu(nt_armcu)
        real rt_armcu(nt_armcu),aqt_armcu(nt_armcu)
! outputs:
        real fs_prof,fl_prof,at_prof,rt_prof,aqt_prof
! local:
        integer it_armcu1, it_armcu2,k
        real timeit,time_armcu1,time_armcu2,frac

! Check that initial day of the simulation consistent with ARMCU period:
       if (annee_ref.ne.1997 ) then
        print*,'Pour ARMCU, annee_ref doit etre 1997 '
        print*,'Changer annee_ref dans run.def'
        stop
       endif

      timeit=(day-day_ini_armcu)*86400

! Determine the closest observation times:
       it_armcu1=INT(timeit/dt_armcu)+1
       it_armcu2=it_armcu1 + 1
       time_armcu1=(it_armcu1-1)*dt_armcu
       time_armcu2=(it_armcu2-1)*dt_armcu
       print *,'timeit day day_ini_armcu',timeit,day,day_ini_armcu
       print *,'it_armcu1,it_armcu2,time_armcu1,time_armcu2',              &
     &          it_armcu1,it_armcu2,time_armcu1,time_armcu2

       if (it_armcu1 .ge. nt_armcu) then
        write(*,*) 'PB-stop: day, it_armcu1, it_armcu2, timeit: '          &
     &        ,day,it_armcu1,it_armcu2,timeit/86400.
        stop
       endif

! time interpolation:
       frac=(time_armcu2-timeit)/(time_armcu2-time_armcu1)
       frac=max(frac,0.0)

       fs_prof = fs_armcu(it_armcu2)                                       &
     &          -frac*(fs_armcu(it_armcu2)-fs_armcu(it_armcu1))
       fl_prof = fl_armcu(it_armcu2)                                       &
     &          -frac*(fl_armcu(it_armcu2)-fl_armcu(it_armcu1))
       at_prof = at_armcu(it_armcu2)                                       &
     &          -frac*(at_armcu(it_armcu2)-at_armcu(it_armcu1))
       rt_prof = rt_armcu(it_armcu2)                                       &
     &          -frac*(rt_armcu(it_armcu2)-rt_armcu(it_armcu1))
       aqt_prof = aqt_armcu(it_armcu2)                                       &
     &          -frac*(aqt_armcu(it_armcu2)-aqt_armcu(it_armcu1))

         print*,                                                           &
     &'day,annee_ref,day_ini_armcu,timeit,it_armcu1,it_armcu2,SST:',       &
     &day,annee_ref,day_ini_armcu,timeit/86400.,it_armcu1,                 &
     &it_armcu2,fs_prof,fl_prof,at_prof,rt_prof,aqt_prof

        return
        END

!=====================================================================
      subroutine readprofiles(nlev_max,kmax,ntrac,height,                  &
     &           thlprof,qtprof,uprof,                                     &
     &           vprof,e12prof,ugprof,vgprof,                              &
     &           wfls,dqtdxls,dqtdyls,dqtdtls,                             &
     &           thlpcar,tracer,nt1,nt2)
      implicit none

        integer nlev_max,kmax,kmax2,ntrac
        logical :: llesread = .true.

        real height(nlev_max),thlprof(nlev_max),qtprof(nlev_max),          &
     &       uprof(nlev_max),vprof(nlev_max),e12prof(nlev_max),            &
     &       ugprof(nlev_max),vgprof(nlev_max),wfls(nlev_max),             &
     &       dqtdxls(nlev_max),dqtdyls(nlev_max),dqtdtls(nlev_max),        &
     &           thlpcar(nlev_max),tracer(nlev_max,ntrac)

        real height1(nlev_max)

        integer, parameter :: ilesfile=1
        integer :: ierr,k,itrac,nt1,nt2

        if(.not.(llesread)) return

       open (ilesfile,file='prof.inp.001',status='old',iostat=ierr)
        if (ierr /= 0) stop 'ERROR:Prof.inp does not exist'
        read (ilesfile,*) kmax
        do k=1,kmax
          read (ilesfile,*) height1(k),thlprof(k),qtprof (k),               &
     &                      uprof (k),vprof  (k),e12prof(k)
        enddo
        close(ilesfile)

       open(ilesfile,file='lscale.inp.001',status='old',iostat=ierr)
        if (ierr /= 0) stop 'ERROR:Lscale.inp does not exist'
        read (ilesfile,*) kmax2
        if (kmax .ne. kmax2) then
          print *, 'fichiers prof.inp et lscale.inp incompatibles :'
          print *, 'nbre de niveaux : ',kmax,' et ',kmax2
          stop 'lecture profiles'
        endif
        do k=1,kmax
          read (ilesfile,*) height(k),ugprof(k),vgprof(k),wfls(k),         &
     &                      dqtdxls(k),dqtdyls(k),dqtdtls(k),thlpcar(k)
        end do
        do k=1,kmax
          if (height(k) .ne. height1(k)) then
            print *, 'fichiers prof.inp et lscale.inp incompatibles :'
            print *, 'les niveaux different : ',k,height1(k), height(k)
            stop
          endif
        end do
        close(ilesfile)

       open(ilesfile,file='trac.inp.001',status='old',iostat=ierr)
        if (ierr /= 0) then
            print*,'WARNING : trac.inp does not exist'
        else
        read (ilesfile,*) kmax2,nt1,nt2
        if (nt2>ntrac) then
          stop 'Augmenter le nombre de traceurs dans traceur.def'
        endif
        if (kmax .ne. kmax2) then
          print *, 'fichiers prof.inp et lscale.inp incompatibles :'
          print *, 'nbre de niveaux : ',kmax,' et ',kmax2
          stop 'lecture profiles'
        endif
        do k=1,kmax
          read (ilesfile,*) height(k),(tracer(k,itrac),itrac=nt1,nt2)
        end do
        close(ilesfile)
        endif

        return
        end
!======================================================================
      subroutine readprofile_sandu(nlev_max,kmax,height,pprof,tprof,       &
     &       thlprof,qprof,uprof,vprof,wprof,omega,o3mmr)
!======================================================================
      implicit none

        integer nlev_max,kmax
        logical :: llesread = .true.

        real height(nlev_max),pprof(nlev_max),tprof(nlev_max)
        real thlprof(nlev_max)
        real qprof(nlev_max),uprof(nlev_max),vprof(nlev_max)
        real wprof(nlev_max),omega(nlev_max),o3mmr(nlev_max)

        integer, parameter :: ilesfile=1
        integer :: k,ierr

        if(.not.(llesread)) return

       open (ilesfile,file='prof.inp.001',status='old',iostat=ierr)
        if (ierr /= 0) stop 'ERROR:Prof.inp does not exist'
        read (ilesfile,*) kmax
        do k=1,kmax
          read (ilesfile,*) height(k),pprof(k),  tprof(k),thlprof(k),      &
     &                      qprof (k),uprof(k),  vprof(k),  wprof(k),      &
     &                      omega (k),o3mmr(k)
        enddo
        close(ilesfile)

        return
        end

!======================================================================
      subroutine readprofile_astex(nlev_max,kmax,height,pprof,tprof,       &
     &    thlprof,qvprof,qlprof,qtprof,uprof,vprof,wprof,tkeprof,o3mmr)
!======================================================================
      implicit none

        integer nlev_max,kmax
        logical :: llesread = .true.

        real height(nlev_max),pprof(nlev_max),tprof(nlev_max),             &
     &  thlprof(nlev_max),qlprof(nlev_max),qtprof(nlev_max),               &
     &  qvprof(nlev_max),uprof(nlev_max),vprof(nlev_max),                  &
     &  wprof(nlev_max),tkeprof(nlev_max),o3mmr(nlev_max)

        integer, parameter :: ilesfile=1
        integer :: ierr,k

        if(.not.(llesread)) return

       open (ilesfile,file='prof.inp.001',status='old',iostat=ierr)
        if (ierr /= 0) stop 'ERROR:Prof.inp does not exist'
        read (ilesfile,*) kmax
        do k=1,kmax
          read (ilesfile,*) height(k),pprof(k),  tprof(k),thlprof(k),      &
     &                qvprof (k),qlprof (k),qtprof (k),                    &
     &                uprof(k),  vprof(k),  wprof(k),tkeprof(k),o3mmr(k)
        enddo
        close(ilesfile)

        return
        end



!======================================================================
      subroutine readprofile_armcu(nlev_max,kmax,height,pprof,uprof,       &
     &       vprof,thetaprof,tprof,qvprof,rvprof,aprof,bprof)
!======================================================================
      implicit none

        integer nlev_max,kmax
        logical :: llesread = .true.

        real height(nlev_max),pprof(nlev_max),tprof(nlev_max)
        real thetaprof(nlev_max),rvprof(nlev_max)
        real qvprof(nlev_max),uprof(nlev_max),vprof(nlev_max)
        real aprof(nlev_max+1),bprof(nlev_max+1)

        integer, parameter :: ilesfile=1
        integer, parameter :: ifile=2
        integer :: ierr,jtot,k

        if(.not.(llesread)) return

! Read profiles at full levels
       IF(nlev_max.EQ.19) THEN
       open (ilesfile,file='prof.inp.19',status='old',iostat=ierr)
       print *,'On ouvre prof.inp.19'
       ELSE
       open (ilesfile,file='prof.inp.40',status='old',iostat=ierr)
       print *,'On ouvre prof.inp.40'
       ENDIF
        if (ierr /= 0) stop 'ERROR:Prof.inp does not exist'
        read (ilesfile,*) kmax
        do k=1,kmax
          read (ilesfile,*) height(k)    ,pprof(k),  uprof(k), vprof(k),   &
     &                      thetaprof(k) ,tprof(k), qvprof(k),rvprof(k)
        enddo
        close(ilesfile)

! Vertical coordinates half levels for eta-coordinates (plev = alpha + beta * psurf) 
       IF(nlev_max.EQ.19) THEN
       open (ifile,file='proh.inp.19',status='old',iostat=ierr)
       print *,'On ouvre proh.inp.19'
       if (ierr /= 0) stop 'ERROR:Proh.inp.19 does not exist'
       ELSE
       open (ifile,file='proh.inp.40',status='old',iostat=ierr)
       print *,'On ouvre proh.inp.40'
       if (ierr /= 0) stop 'ERROR:Proh.inp.40 does not exist'
       ENDIF
        read (ifile,*) kmax
        do k=1,kmax
          read (ifile,*) jtot,aprof(k),bprof(k)
        enddo
        close(ifile)

        return
        end

!=====================================================================
      subroutine read_fire(fich_fire,nlevel,ntime                          &
     &     ,zz,thl,qt,u,v,tke                                              &
     &     ,ug,vg,wls,dqtdx,dqtdy,dqtdt,thl_rad)

!program reading forcings of the FIRE case study


      implicit none

#include "netcdf.inc"

      integer ntime,nlevel
      character*80 :: fich_fire
      real*8 zz(nlevel)

      real*8 thl(nlevel)
      real*8 qt(nlevel),u(nlevel)
      real*8 v(nlevel),tke(nlevel)
      real*8 ug(nlevel,ntime),vg(nlevel,ntime),wls(nlevel,ntime)
      real*8 dqtdx(nlevel,ntime),dqtdy(nlevel,ntime)
      real*8 dqtdt(nlevel,ntime),thl_rad(nlevel,ntime)

      integer nid, ierr
      integer nbvar3d
      parameter(nbvar3d=30)
      integer var3didin(nbvar3d)

      ierr = NF_OPEN(fich_fire,NF_NOWRITE,nid)
      if (ierr.NE.NF_NOERR) then
         write(*,*) 'ERROR: Pb opening forcings nc file '
         write(*,*) NF_STRERROR(ierr)
         stop ""
      endif


       ierr=NF_INQ_VARID(nid,"zz",var3didin(1)) 
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'lev'
         endif


      ierr=NF_INQ_VARID(nid,"thetal",var3didin(2))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'temp'
         endif

      ierr=NF_INQ_VARID(nid,"qt",var3didin(3))
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

      ierr=NF_INQ_VARID(nid,"tke",var3didin(6))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'tke'
         endif

      ierr=NF_INQ_VARID(nid,"ugeo",var3didin(7))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'ug'
         endif

      ierr=NF_INQ_VARID(nid,"vgeo",var3didin(8))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'vg'
         endif
      
      ierr=NF_INQ_VARID(nid,"wls",var3didin(9))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'wls'
         endif

      ierr=NF_INQ_VARID(nid,"dqtdx",var3didin(10))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'dqtdx'
         endif

      ierr=NF_INQ_VARID(nid,"dqtdy",var3didin(11))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'dqtdy'
      endif

      ierr=NF_INQ_VARID(nid,"dqtdt",var3didin(12))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'dqtdt'
      endif

      ierr=NF_INQ_VARID(nid,"thl_rad",var3didin(13))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'thl_rad'
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
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(2),thl)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(2),thl)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture thl ok',thl

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(3),qt)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(3),qt)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture qt ok',qt
 
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
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(6),tke)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(6),tke)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture tke ok',tke

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(7),ug)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(7),ug)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture ug ok',ug

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(8),vg)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(8),vg)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture vg ok',vg

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(9),wls)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(9),wls)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture wls ok',wls

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(10),dqtdx)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(10),dqtdx)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dqtdx ok',dqtdx

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(11),dqtdy)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(11),dqtdy)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dqtdy ok',dqtdy

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(12),dqtdt)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(12),dqtdt)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture dqtdt ok',dqtdt

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(13),thl_rad)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(13),thl_rad)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture thl_rad ok',thl_rad

         return 
         end subroutine read_fire
!=====================================================================
      subroutine read_dice(fich_dice,nlevel,ntime                         &
     &     ,zz,pres,t,qv,u,v,o3                                          &
     &     ,shf,lhf,lwup,swup,tg,ustar,psurf,ug,vg                        &
     &     ,hadvt,hadvq,hadvu,hadvv,w,omega)

!program reading initial profils and forcings of the Dice case study


      implicit none

#include "netcdf.inc"
#include "YOMCST.h"

      integer ntime,nlevel
      integer l,k
      character*80 :: fich_dice
      real*8 time(ntime)
      real*8 zz(nlevel)

      real*8 th(nlevel),pres(nlevel),t(nlevel)
      real*8 qv(nlevel),u(nlevel),v(nlevel),o3(nlevel)
      real*8 shf(ntime),lhf(ntime),lwup(ntime),swup(ntime),tg(ntime)
      real*8 ustar(ntime),psurf(ntime),ug(ntime),vg(ntime)
      real*8 hadvt(nlevel,ntime),hadvq(nlevel,ntime),hadvu(nlevel,ntime)
      real*8 hadvv(nlevel,ntime),w(nlevel,ntime),omega(nlevel,ntime)
      real*8 pzero

      integer nid, ierr
      integer nbvar3d
      parameter(nbvar3d=30)
      integer var3didin(nbvar3d)

      pzero=100000.
      ierr = NF_OPEN(fich_dice,NF_NOWRITE,nid)
      if (ierr.NE.NF_NOERR) then
         write(*,*) 'ERROR: Pb opening forcings nc file '
         write(*,*) NF_STRERROR(ierr)
         stop ""
      endif


       ierr=NF_INQ_VARID(nid,"height",var3didin(1)) 
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'height'
         endif

       ierr=NF_INQ_VARID(nid,"pf",var3didin(11)) 
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'pf'
         endif

      ierr=NF_INQ_VARID(nid,"theta",var3didin(12))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'theta'
         endif

      ierr=NF_INQ_VARID(nid,"qv",var3didin(13))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'qv'
         endif

      ierr=NF_INQ_VARID(nid,"u",var3didin(14))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'u'
         endif

      ierr=NF_INQ_VARID(nid,"v",var3didin(15))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'v'
         endif

      ierr=NF_INQ_VARID(nid,"o3mmr",var3didin(16))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'o3'
         endif

      ierr=NF_INQ_VARID(nid,"shf",var3didin(2))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'shf'
         endif

      ierr=NF_INQ_VARID(nid,"lhf",var3didin(3))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'lhf'
         endif
      
      ierr=NF_INQ_VARID(nid,"lwup",var3didin(4))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'lwup'
         endif

      ierr=NF_INQ_VARID(nid,"swup",var3didin(5))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'dqtdx'
         endif

      ierr=NF_INQ_VARID(nid,"Tg",var3didin(6))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'Tg'
      endif

      ierr=NF_INQ_VARID(nid,"ustar",var3didin(7))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'ustar'
      endif

      ierr=NF_INQ_VARID(nid,"psurf",var3didin(8))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'psurf'
      endif

      ierr=NF_INQ_VARID(nid,"Ug",var3didin(9))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'Ug'
      endif

      ierr=NF_INQ_VARID(nid,"Vg",var3didin(10))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'Vg'
      endif

      ierr=NF_INQ_VARID(nid,"hadvT",var3didin(17))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hadvT'
      endif

      ierr=NF_INQ_VARID(nid,"hadvq",var3didin(18))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hadvq'
      endif

      ierr=NF_INQ_VARID(nid,"hadvu",var3didin(19))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hadvu'
      endif

      ierr=NF_INQ_VARID(nid,"hadvv",var3didin(20))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hadvv'
      endif

      ierr=NF_INQ_VARID(nid,"w",var3didin(21))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'w'
      endif

      ierr=NF_INQ_VARID(nid,"omega",var3didin(22))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'omega'
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
!          write(*,*)'lecture zz ok',zz
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(11),pres)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(11),pres)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture pres ok',pres

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(12),th)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(12),th)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture th ok',th
           do k=1,nlevel
             t(k)=th(k)*(pres(k)/pzero)**rkappa
           enddo

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(13),qv)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(13),qv)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture qv ok',qv
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(14),u)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(14),u)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture u ok',u

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(15),v)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(15),v)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture v ok',v

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(16),o3)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(16),o3)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture o3 ok',o3

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(2),shf)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(2),shf)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture shf ok',shf

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(3),lhf)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(3),lhf)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture lhf ok',lhf

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(4),lwup)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(4),lwup)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture lwup ok',lwup

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(5),swup)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(5),swup)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture swup ok',swup

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(6),tg)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(6),tg)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture tg ok',tg

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(7),ustar)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(7),ustar)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture ustar ok',ustar

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(8),psurf)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(8),psurf)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture psurf ok',psurf

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(9),ug)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(9),ug)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture ug ok',ug

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(10),vg)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(10),vg)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture vg ok',vg

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(17),hadvt)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(17),hadvt)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture hadvt ok',hadvt

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(18),hadvq)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(18),hadvq)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture hadvq ok',hadvq

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(19),hadvu)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(19),hadvu)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture hadvu ok',hadvu

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(20),hadvv)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(20),hadvv)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture hadvv ok',hadvv

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(21),w)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(21),w)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture w ok',w

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(22),omega)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(22),omega)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
!          write(*,*)'lecture omega ok',omega

         return 
         end subroutine read_dice
!=====================================================================
      subroutine read_gabls4(fich_gabls4,nlevel,ntime,nsol                    &
     &     ,zz,depth_sn,ug,vg,pf,th,t,qv,u,v,hadvt,hadvq,tg,tsnow,snow_dens)

!program reading initial profils and forcings of the Gabls4 case study


      implicit none

#include "netcdf.inc"

      integer ntime,nlevel,nsol
      integer l,k
      character*80 :: fich_gabls4
      real*8 time(ntime)

!  ATTENTION: visiblement quand on lit gabls4_driver.nc on recupere les donnees 
! dans un ordre inverse par rapport a la convention LMDZ
! ==> il faut tout inverser  (MPL 20141024)
! les variables indexees "_i" sont celles qui sont lues dans gabls4_driver.nc
      real*8 zz_i(nlevel),th_i(nlevel),pf_i(nlevel),t_i(nlevel)
      real*8 qv_i(nlevel),u_i(nlevel),v_i(nlevel),ug_i(nlevel,ntime),vg_i(nlevel,ntime)
      real*8 hadvt_i(nlevel,ntime),hadvq_i(nlevel,ntime)

      real*8 zz(nlevel),th(nlevel),pf(nlevel),t(nlevel)
      real*8 qv(nlevel),u(nlevel),v(nlevel),ug(nlevel,ntime),vg(nlevel,ntime)
      real*8 hadvt(nlevel,ntime),hadvq(nlevel,ntime)

      real*8 depth_sn(nsol),tsnow(nsol),snow_dens(nsol)
      real*8 tg(ntime)
      integer nid, ierr
      integer nbvar3d
      parameter(nbvar3d=30)
      integer var3didin(nbvar3d)

      ierr = NF_OPEN(fich_gabls4,NF_NOWRITE,nid)
      if (ierr.NE.NF_NOERR) then
         write(*,*) 'ERROR: Pb opening forcings nc file '
         write(*,*) NF_STRERROR(ierr)
         stop ""
      endif


       ierr=NF_INQ_VARID(nid,"height",var3didin(1)) 
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'height'
         endif

      ierr=NF_INQ_VARID(nid,"depth_sn",var3didin(2))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'depth_sn'
      endif

      ierr=NF_INQ_VARID(nid,"Ug",var3didin(3))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'Ug'
      endif

      ierr=NF_INQ_VARID(nid,"Vg",var3didin(4))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'Vg'
      endif
       ierr=NF_INQ_VARID(nid,"pf",var3didin(5)) 
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'pf'
         endif

      ierr=NF_INQ_VARID(nid,"theta",var3didin(6))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'theta'
         endif

      ierr=NF_INQ_VARID(nid,"tempe",var3didin(7))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'tempe'
         endif

      ierr=NF_INQ_VARID(nid,"qv",var3didin(8))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'qv'
         endif

      ierr=NF_INQ_VARID(nid,"u",var3didin(9))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'u'
         endif

      ierr=NF_INQ_VARID(nid,"v",var3didin(10))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'v'
         endif

      ierr=NF_INQ_VARID(nid,"hadvT",var3didin(11))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hadvt'
         endif

      ierr=NF_INQ_VARID(nid,"hadvQ",var3didin(12))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'hadvq'
      endif

      ierr=NF_INQ_VARID(nid,"Tsnow",var3didin(14))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'tsnow'
      endif

      ierr=NF_INQ_VARID(nid,"snow_density",var3didin(15))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'snow_density'
      endif

      ierr=NF_INQ_VARID(nid,"Tg",var3didin(16))
         if(ierr/=NF_NOERR) then
           write(*,*) NF_STRERROR(ierr)
           stop 'Tg'
      endif


!dimensions lecture
!      call catchaxis(nid,ntime,nlevel,time,z,ierr)
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(1),zz_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(1),zz_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(2),depth_sn)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(2),depth_sn)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(3),ug_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(3),ug_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(4),vg_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(4),vg_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(5),pf_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(5),pf_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(6),th_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(6),th_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(7),t_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(7),t_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(8),qv_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(8),qv_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(9),u_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(9),u_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(10),v_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(10),v_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(11),hadvt_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(11),hadvt_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(12),hadvq_i)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(12),hadvq_i)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(14),tsnow)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(14),tsnow)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif
 
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(15),snow_dens)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(15),snow_dens)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif

#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid,var3didin(16),tg)
#else
         ierr = NF_GET_VAR_REAL(nid,var3didin(16),tg)
#endif
         if(ierr/=NF_NOERR) then
            write(*,*) NF_STRERROR(ierr)
            stop "getvarup"
         endif

! On remet les variables lues dans le bon ordre des niveaux (MPL 20141024)
         do k=1,nlevel
           zz(k)=zz_i(nlevel+1-k)
           ug(k,:)=ug_i(nlevel+1-k,:)
           vg(k,:)=vg_i(nlevel+1-k,:)
           pf(k)=pf_i(nlevel+1-k)
           print *,'pf=',pf(k)
           th(k)=th_i(nlevel+1-k)
           t(k)=t_i(nlevel+1-k)
           qv(k)=qv_i(nlevel+1-k)
           u(k)=u_i(nlevel+1-k)
           v(k)=v_i(nlevel+1-k)
           hadvt(k,:)=hadvt_i(nlevel+1-k,:)
           hadvq(k,:)=hadvq_i(nlevel+1-k,:)
         enddo
         return 
 end subroutine read_gabls4
!=====================================================================

!     Reads CIRC input files      

      SUBROUTINE read_circ(nlev_circ,cf,lwp,iwp,reliq,reice,t,z,p,pm,h2o,o3,sza)
      
      parameter (ncm_1=49180)
#include "YOMCST.h"

      real albsfc(ncm_1), albsfc_w(ncm_1)
      real cf(nlev_circ), icefra(nlev_circ), deice(nlev_circ), &
           reliq(nlev_circ), reice(nlev_circ), lwp(nlev_circ), iwp(nlev_circ)
      real t(nlev_circ+1), z(nlev_circ+1), dz(nlev_circ), p(nlev_circ+1)
      real aer_beta(nlev_circ), waer(nlev_circ), gaer(nlev_circ)
      real pm(nlev_circ), tm(nlev_circ), h2o(nlev_circ), o3(nlev_circ)
      real co2(nlev_circ), n2o(nlev_circ), co(nlev_circ), ch4(nlev_circ), &
           o2(nlev_circ), ccl4(nlev_circ), f11(nlev_circ), f12(nlev_circ)
!     za= zenital angle
!     sza= cosinus angle zenital
      real wavn(ncm_1), ssf(ncm_1),za,sza
      integer nlev


!     Open the files

      open (11, file='Tsfc_sza_nlev_case.txt', status='old')
      open (12, file='level_input_case.txt', status='old')
      open (13, file='layer_input_case.txt', status='old')
      open (14, file='aerosol_input_case.txt', status='old')
      open (15, file='cloud_input_case.txt', status='old')
      open (16, file='sfcalbedo_input_case.txt', status='old')
      
!     Read scalar information
      do iskip=1,5
         read (11, *)
      enddo
      read (11, '(i8)') nlev
      read (11, '(f10.2)') tsfc
      read (11, '(f10.2)') za
      read (11, '(f10.4)') sw_dn_toa
      sza=cos(za/180.*RPI)
      print *,'nlev,tsfc,sza,sw_dn_toa,RPI',nlev,tsfc,sza,sw_dn_toa,RPI
      close(11)

!     Read level information
      read (12, *)
      do il=1,nlev
         read (12, 302) ilev, z(il), p(il), t(il)
         z(il)=z(il)*1000.    ! z donne en km
         p(il)=p(il)*100.     ! p donne en mb
      enddo
302   format (i8, f8.3, 2f9.2)
      close(12)

!     Read layer information (midpoint values)
      do iskip=1,3
         read (13, *)
      enddo
      do il=1,nlev-1
         read (13, 303) ilev,pm(il),tm(il),h2o(il),co2(il),o3(il), &
                        n2o(il),co(il),ch4(il),o2(il),ccl4(il), &
                        f11(il),f12(il)
         pm(il)=pm(il)*100.
      enddo
303   format (i8, 2f9.2, 10(2x,e13.7))      
      close(13)
      
!     Read aerosol layer information
      do iskip=1,3
         read (14, *)
      enddo
      read (14, '(f10.2)') aer_alpha
      read (14, *)
      read (14, *)
      do il=1,nlev-1
         read (14, 304) ilev, aer_beta(il), waer(il), gaer(il)
      enddo
304   format (i8, f9.5, 2f8.3)
      close(14)
      
!     Read cloud information
      do iskip=1,3
         read (15, *)
      enddo
      do il=1,nlev-1
         read (15, 305) ilev, cf(il), lwp(il), iwp(il), reliq(il), reice(il)
         lwp(il)=lwp(il)/1000.          ! lwp donne en g/kg
         iwp(il)=iwp(il)/1000.          ! iwp donne en g/kg
         reliq(il)=reliq(il)/1000000.   ! reliq donne en microns
         reice(il)=reice(il)/1000000.   ! reice donne en microns
      enddo
305   format (i8, f8.3, 4f9.2)
      close(15)

!     Read surface albedo (weighted & unweighted) and spectral solar irradiance
      do iskip=1,6
         read (16, *)
      enddo
      do icm_1=1,ncm_1
         read (16, 306) wavn(icm_1), albsfc(icm_1), albsfc_w(icm_1), ssf(icm_1)
      enddo
306   format(f10.1, 2f12.5, f14.8)
      close(16)
 
      return 
      end subroutine read_circ
!=====================================================================
!     Reads RTMIP input files      

      SUBROUTINE read_rtmip(nlev_rtmip,play,plev,t,h2o,o3)
      
#include "YOMCST.h"

      real t(nlev_rtmip), pt(nlev_rtmip),pb(nlev_rtmip),h2o(nlev_rtmip), o3(nlev_rtmip)
      real temp(nlev_rtmip), play(nlev_rtmip),ovap(nlev_rtmip), oz(nlev_rtmip),plev(nlev_rtmip+1)
      integer nlev


!     Open the files

      open (11, file='low_resolution_profile.txt', status='old')
      
!     Read level information
      read (11, *)
      do il=1,nlev_rtmip
         read (11, 302) pt(il), pb(il), t(il),h2o(il),o3(il)
      enddo
      do il=1,nlev_rtmip
         play(il)=pt(nlev_rtmip-il+1)*100.     ! p donne en mb
         temp(il)=t(nlev_rtmip-il+1)
         ovap(il)=h2o(nlev_rtmip-il+1)
         oz(il)=o3(nlev_rtmip-il+1)
      enddo
      do il=1,39
         plev(il)=play(il)+(play(il+1)-play(il))/2.
         print *,'il p t ovap oz=',il,plev(il),temp(il),ovap(il),oz(il)
      enddo
      plev(41)=101300.
302   format (e16.10,3x,e16.10,3x,e16.10,3x,e12.6,3x,e12.6)
      close(12)
 
      return 
      end subroutine read_rtmip
!=====================================================================
