      SUBROUTINE read_vent(debutphy, step, nbjour, u10m_ec, v10m_ec)
      USE dimphy
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
!      USE write_field_phy
      IMPLICIT NONE 
#include "dimensions.h"
c #include "dimphy.h"
#include "paramet.h"
#include "netcdf.inc"
c
      INTEGER step, nbjour
      LOGICAL debutphy
      real u10m_ec(klon), v10m_ec(klon)
      real u10m_ec_glo(klon_glo), v10m_ec_glo(klon_glo)
c
!      real u10m_nc(iip1,jjp1) !, v10m_nc(iip1,jjm) ! dim 97x72
!      real v10m_nc(iip1,jjp1)  ! dim 97x73
      real u10m_nc_glo(nbp_lon+1,nbp_lat) !, v10m_nc(iip1,jjm) ! dim 97x72
      real v10m_nc_glo(nbp_lon+1,nbp_lat)  ! dim 97x73
      real rcode
      integer ncidu1, varidu1, ncidv1, varidv1
      save ncidu1, varidu1, ncidv1, varidv1
!$OMP THREADPRIVATE(ncidu1, varidu1, ncidv1, varidv1)
      integer start(4),count(4), status
      integer i, j, ig


c
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      if (debutphy) then
c
         ncidu1=NCOPN('u10m.nc',NCNOWRIT,rcode)
         varidu1=NCVID(ncidu1,'UWND',rcode)
         ncidv1=NCOPN('v10m.nc',NCNOWRIT,rcode)
         varidv1=NCVID(ncidv1,'VWND',rcode)
c
      endif
c
      start(1)=1
      start(2)=1
      start(4)=0

!      count(1)=iip1
      count(1)=nbp_lon+1
!      count(2)=jjp1
      count(2)=nbp_lat
      count(3)=1
      count(4)=0
c
      start(3)=step
c
#ifdef NC_DOUBLE
!      status=NF_GET_VARA_DOUBLE(ncidu1,varidu1,start,count,u10m_nc)
      status=NF_GET_VARA_DOUBLE(ncidu1,varidu1,start,count,u10m_nc_glo)
#else
!      status=NF_GET_VARA_REAL(ncidu1,varidu1,start,count,u10m_nc)
      status=NF_GET_VARA_REAL(ncidu1,varidu1,start,count,u10m_nc_glo)
#endif
!       print *,status
c
#ifdef NC_DOUBLE
!      status=NF_GET_VARA_DOUBLE(ncidv1,varidv1,start,count,v10m_nc)
      status=NF_GET_VARA_DOUBLE(ncidv1,varidv1,start,count,v10m_nc_glo)
#else
!      status=NF_GET_VARA_REAL(ncidv1,varidv1,start,count,v10m_nc)
      status=NF_GET_VARA_REAL(ncidv1,varidv1,start,count,v10m_nc_glo)
#endif
c

!      print *,'beforebidcor u10m_nc', u10m_nc(1,jjp1)
!      print *,'beforebidcor v10m_nc', v10m_nc(1,jjp1)

!       print *,status
!      call correctbid(iim,jjp1,u10m_nc)
!      call correctbid(iim,jjp1,v10m_nc)
      call correctbid(nbp_lon,nbp_lat,u10m_nc_glo)
      call correctbid(nbp_lon,nbp_lat,v10m_nc_glo)

!      print *,'afterbidcor u10m_nc', u10m_nc(1,jjp1)
!      print *,'afterbidcor v10m_nc', v10m_nc(1,jjp1)
c
c--upside down + physical grid
c
!      u10m_ec(1)=u10m_nc(1,jjp1)
!      v10m_ec(1)=v10m_nc(1,jjp1)
      u10m_ec_glo(1)=u10m_nc_glo(1,nbp_lat)
      v10m_ec_glo(1)=v10m_nc_glo(1,nbp_lat)
      ig=2
!      DO j=2,jjm
!         DO i = 1, iim
      DO j=2,nbp_lat-1
         DO i = 1, nbp_lon
!           u10m_ec(ig)=u10m_nc(i,jjp1+1-j)
!           v10m_ec(ig)=v10m_nc(i,jjp1+1-j)
           u10m_ec_glo(ig)=u10m_nc_glo(i,nbp_lat+1-j)
           v10m_ec_glo(ig)=v10m_nc_glo(i,nbp_lat+1-j)
           ig=ig+1
!         print *,u10m_ec(ig) ,v10m_ec(ig) 
         ENDDO
      ENDDO
      u10m_ec_glo(ig)=u10m_nc_glo(1,1)
      v10m_ec_glo(ig)=v10m_nc_glo(1,1)


!      end if master
      ENDIF 
!$OMP END MASTER
!$OMP BARRIER
      CALL scatter(u10m_ec_glo,u10m_ec)
      CALL scatter(v10m_ec_glo,v10m_ec)

!      print *,'JE  tamagno viento ig= ', ig
!      print *,'READ_VENT U = ',SUM(u10m_ec),MINVAL(u10m_ec),
!     .                                      MAXVAL(u10m_ec)
!      print *,'READ_VENT V = ',SUM(v10m_ec),MINVAL(v10m_ec),
!     .                                      MAXVAL(v10m_ec)
!       print *,'u v 1 ', u10m_ec(1),v10m_ec(1)
!       print *,'u v klon ', u10m_ec(klon),v10m_ec(klon)
      RETURN 
      END

c added by JE from the nh SPLA, dyn3d/read_reanalyse.F which is not available any more
      subroutine correctbid(iim,nl,x)
      integer iim,nl
      real x(iim+1,nl)
      integer i,l
      real zz

      do l=1,nl
         do i=2,iim-1
            if(abs(x(i,l)).gt.1.e10) then
               zz=0.5*(x(i-1,l)+x(i+1,l))
c              print*,'correction ',i,l,x(i,l),zz
               x(i,l)=zz
            endif
         enddo
      enddo

      return
      end



