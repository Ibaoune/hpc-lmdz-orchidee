      SUBROUTINE read_dust(debutphy, step, nbjour, dust_ec)
      USE dimphy
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      IMPLICIT NONE 
c
#include "dimensions.h"
#include "paramet.h"
#include "netcdf.inc"
c
      INTEGER step, nbjour
      LOGICAL debutphy
      real dust_ec(klon)
      real dust_ec_glo(klon_glo)
c
      real dust_nc(iip1,jjp1)
      real dust_nc_glo(nbp_lon+1,nbp_lat)
      real rcode
      integer ncid1, varid1, ncid2, varid2 

      save ncid1, varid1, ncid2, varid2
!$OMP THREADPRIVATE(ncid1, varid1, ncid2, varid2)
      integer start(4),count(4), status
      integer i, j, ig
c
!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
      if (debutphy) then
c
         ncid1=NCOPN('dust.nc',NCNOWRIT,rcode)
         varid1=NCVID(ncid1,'EMISSION',rcode)
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
!      status=NF_GET_VARA_DOUBLE(ncid1,varid1,start,count,dust_nc)
      status=NF_GET_VARA_DOUBLE(ncid1,varid1,start,count,dust_nc_glo)
#else
!      status=NF_GET_VARA_REAL(ncid1,varid1,start,count,dust_nc)
      status=NF_GET_VARA_REAL(ncid1,varid1,start,count,dust_nc_glo)
#endif
c
!      call correctbid(iim,jjp1,dust_nc)
      call correctbid(nbp_lon,nbp_lat,dust_nc_glo)
c
c--upside down + physical grid
c
c--OB=change jjp1 to 1 here
!      dust_ec(1)=MAX(dust_nc(1,jjp1),0.0)
      dust_ec(1)=MAX(dust_nc(1,nbp_lat),0.0)
      ig=2
!      DO j=2,jjm
      DO j=2,nbp_lat-1
!         DO i = 1, iim
         DO i = 1, nbp_lon
c--OB=change jjp1+1-j to j here
!           dust_ec(ig)=MAX(dust_nc(i,jjp1+1-j),0.0)
           dust_ec_glo(ig)=MAX(dust_nc_glo(i,nbp_lat+1-j),0.0)
           ig=ig+1 
         ENDDO
      ENDDO
c--OB=change second 1 to jjp1 here
      dust_ec_glo(ig)=MAX(dust_nc_glo(1,1),0.0)
!      end if master
      ENDIF
!$OMP END MASTER
!$OMP BARRIER
      CALL scatter(dust_ec_glo,dust_ec)
c
      RETURN 
      END
