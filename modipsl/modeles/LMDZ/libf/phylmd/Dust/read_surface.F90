       subroutine read_surface(name,surfa)

     
! common
! ------
       USE ioipsl
!       USE comgeomphy
       USE dimphy
       USE mod_grid_phy_lmdz
       USE mod_phys_lmdz_para
       USE iophy
!       USE netcdf
       IMPLICIT NONE

       include "netcdf.inc"
#include "dimensions.h"
#include "paramet.h"

       character*10 name
       character*10 varname
!
       real tmp_dyn(iip1,jjp1)
       real tmp_dyn_glo(nbp_lon+1,nbp_lat)
!       real tmp_dyn_glo(nbp_lon,nbp_lat)
       REAL tmp_dyn_invers(iip1,jjp1)
       real tmp_dyn_invers_glo(nbp_lon+1,nbp_lat)
!       real tmp_dyn_invers_glo(nbp_lon,nbp_lat)
       real tmp_fi(klon)
       real tmp_fi_glo(klon_glo)
       real surfa(klon,5)
       real surfa_glo(klon_glo,5)
!
       integer ncid
       integer varid
       real rcode
       integer start(2),count(2),status
       integer i,j,l,ig
       character*1 str1

!JE20140526<<
      character*4 ::  latstr,aux4s
      logical :: outcycle, isinversed
      real, dimension(jjp1) :: lats
      real, dimension(nbp_lat) :: lats_glo
      real :: rcode2
      integer, dimension(1) :: startj,endj
!JE20140526>>
!$OMP MASTER
       IF (is_mpi_root .AND. is_omp_root) THEN

       print*,'Lecture du fichier donnees_lisa.nc' 
       ncid=NCOPN('donnees_lisa.nc',NCNOWRIT,rcode)

!JE20140526<<: check if are inversed or not the latitude grid in donnes_lisa
      outcycle=.false.
      latstr='null'
      isinversed=.false.
      do i=1,5
       if (i==1) aux4s='latu'
       if (i==2) aux4s='LATU'
       if (i==3) aux4s='LatU'
       if (i==4) aux4s='Latu'
       if (i==5) aux4s='latU'
       status = NF_INQ_VARID (ncid, aux4s, rcode)
!       print *,'stat,i',status,i,outcycle,aux4s
!       print *,'ifclause',status.NE. NF_NOERR ,outcycle == .false.
       IF ((.not.(status.NE. NF_NOERR) ).and.( .not. outcycle )) THEN
         outcycle=.true.
         latstr=aux4s
       ENDIF
      enddo ! check if it inversed lat
      startj(1)=1
!      endj(1)=jjp1
      endj(1)=nbp_lat
      varid=NCVID(ncid,latstr,rcode)

#ifdef NC_DOUBLE
          status=NF_GET_VARA_DOUBLE(ncid,varid,startj,endj,lats_glo)
#else
          status=NF_GET_VARA_REAL(ncid,varid,startj,endj,lats_glo)
#endif
!      print *,latstr,varid,status,jjp1,rcode
!      IF (status .NE. NF_NOERR) print*,'NOOOOOOO'
!      print *,lats
!stop

! check if netcdf is latitude inversed or not. 
      if (lats_glo(1)<lats_glo(2)) isinversed=.true.
! JE20140526>>


       DO i=1,5
          write(str1,'(i1)') i
          varname=trim(name)//str1
       print*,'lecture variable:',varname
          varid=NCVID(ncid,trim(varname),rcode)
!          varid=NCVID(ncid,varname,rcode)

!  dimensions pour les champs scalaires et le vent zonal
!  -----------------------------------------------------

          start(1)=1
          start(2)=1      
          count(1)=nbp_lon+1
!          count(1)=iip1
          count(2)=nbp_lat
!          count(2)=jjp1

! mise a zero des tableaux 
! ------------------------
          tmp_dyn(:,:)=0.0
          tmp_fi(:)=0.0
! Lecture
! -----------------------
#ifdef NC_DOUBLE
          status=NF_GET_VARA_DOUBLE(ncid,varid,start,count,tmp_dyn_glo)
#else
          status=NF_GET_VARA_REAL(ncid,varid,start,count,tmp_dyn_glo)
#endif

!      call dump2d(iip1,jjp1,tmp_dyn,'tmp_dyn   ')
       DO j=1, nbp_lat
          DO ig=1, nbp_lon+1
             tmp_dyn_invers_glo(ig,j)=tmp_dyn_glo(ig,nbp_lat-j+1)
          ENDDO
       ENDDO

       
!JE20140522!          call gr_dyn_fi_p(1, iip1, jjp1, klon, tmp_dyn_invers, tmp_fi)

!JE20140526<<
!              call gr_dyn_fi(1, iip1, jjp1, klon, tmp_dyn_invers, tmp_fi)
           if (isinversed) then
                        call gr_dyn_fi(1, nbp_lon+1, nbp_lat, klon_glo, &
     & tmp_dyn_invers_glo, tmp_fi_glo)
!              call gr_dyn_fi(1, iip1, jjp1, klon, tmp_dyn_invers, tmp_fi)
!              call gr_dyn_fi_p(1, iip1, jjp1, klon, tmp_dyn_invers, tmp_fi)
           else      
                        call gr_dyn_fi(1, nbp_lon+1, nbp_lat, klon_glo, &
     &   tmp_dyn_glo, tmp_fi_glo)
!              call gr_dyn_fi(1, iip1, jjp1, klon, tmp_dyn, tmp_fi)
!              call gr_dyn_fi_p(1, iip1, jjp1, klon, tmp_dyn, tmp_fi)
           endif
!JE20140526>>
!      call dump2d(iim,jjm-1,tmp_fi(2),'tmp_fi   ')
!
          DO j=1,klon_glo

                surfa_glo(j,i)=tmp_fi_glo(j)

          ENDDO ! Fin de recopie du tableau
!
       ENDDO ! Fin boucle 1 a 5
       print*,'Passage Grille Dyn -> Phys'


      ENDIF !mpi 
!$OMP END MASTER
!$OMP BARRIER
      call scatter(surfa_glo,surfa)


       return
       end subroutine read_surface
