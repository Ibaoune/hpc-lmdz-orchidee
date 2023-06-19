MODULE cfmip_point_locations
  IMPLICIT NONE

CONTAINS

 SUBROUTINE read_CFMIP_point_locations(npCFMIP, tab, lonCFMIP, latCFMIP)
  IMPLICIT none
  INTEGER :: npCFMIP
  REAL, DIMENSION(npCFMIP) :: lonCFMIP, latCFMIP
  INTEGER :: tab(npCFMIP), np

  WRITE(*,*) 'npCFMIP=',npCFMIP
! OPEN(20, file="pointlocations.txt",status='old')
  OPEN(20, file="pointlocations.txt",status='old',err=999)
  OPEN(21, file="pointlocations_lon180.txt")
  np=1
10 READ(20,*) tab(np), lonCFMIP(np), latCFMIP(np)
!!! passage de 0-360 a -180/180
   IF (lonCFMIP(np).GT.180.) THEN
    lonCFMIP(np)=lonCFMIP(np)-360.
   ENDIF 
   WRITE(21,*) np, lonCFMIP(np), latCFMIP(np)
   np=np+1
   IF(np.LE.npCFMIP) THEN 
    GOTO 10
   ENDIF
   CLOSE(20)
   CLOSE(21)
999 RETURN
 END SUBROUTINE read_CFMIP_point_locations

 SUBROUTINE LMDZ_CFMIP_point_locations(npCFMIP, lonCFMIP, latCFMIP, &
  tabijGCM, lonGCM, latGCM, ipt, jpt)
  USE dimphy
  USE iophy
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, klon_glo

  IMPLICIT none 
  INTEGER :: npCFMIP
  REAL, DIMENSION(npCFMIP) :: lonCFMIP, latCFMIP
  INTEGER :: i, j, np, ip
  INTEGER, DIMENSION(npCFMIP) :: ipt, jpt 
  REAL :: dlon1, dlon2
  REAL :: dlat1, dlat2
  REAL, DIMENSION(nbp_lon+1) :: lon
  INTEGER, DIMENSION(npCFMIP) :: tabijGCM
  REAL, DIMENSION(npCFMIP) :: lonGCM, latGCM

  lon(1:nbp_lon)=io_lon(:)
  lon(nbp_lon+1)=-1*lon(1)
  OPEN(22, file="LMDZ_pointsCFMIP.txt")
  DO np=1, npCFMIP
  DO i=1, nbp_lon
!
! PRINT*,'IM np i lonCF lonGCM lonGCM+1',np,i,lonCFMIP(np),lon(i), &
!  lon(i+1)
!
   IF(lonCFMIP(np).GE.lon(i).AND.lonCFMIP(np).LT.lon(i+1)) THEN
    dlon1 = abs (lonCFMIP(np) - lon(i))
    dlon2 = abs (lonCFMIP(np) - lon(i+1))
    IF (dlon1.LE.dlon2) THEN
     ipt(np)=i
    ELSE
     ipt(np)=i+1
    ENDIF
   ENDIF
  END DO
  END DO
!
   np=1
30 j=1
40 IF(latCFMIP(np).LE.io_lat(j).AND.latCFMIP(np).GE.io_lat(j+1)) THEN
    dlat1 = abs (latCFMIP(np) - io_lat(j))
    dlat2 = abs (latCFMIP(np) - io_lat(j+1))
    IF (dlat1.LE.dlat2) THEN
     jpt(np)=j
    ELSE
     jpt(np)=j+1
    ENDIF
    np=np+1
    IF(np.LE.npCFMIP) THEN
     GOTO 30
    ENDIF 
   ELSE
    j=j+1
    IF(j.LE.nbp_lat-1) THEN 
     GOTO 40
    ENDIF
   ENDIF

  DO np=1, npCFMIP
   WRITE(22,*) lon(ipt(np)), io_lat(jpt(np))
  ENDDO
  CLOSE(22)

  OPEN(23, file="pointsCFMIPvsLMDZ.txt")
    DO ip=1, npCFMIP
     lonGCM(ip)=lon(ipt(ip))
     latGCM(ip)=io_lat(jpt(ip))
     if(jpt(ip).GE.2.AND.jpt(ip).LE.nbp_lat-1) THEN     
      tabijGCM(ip)=1+(jpt(ip)-2)*nbp_lon+ipt(ip)
     else if(jpt(ip).EQ.1) THEN
      tabijGCM(ip)=1
     else if(jpt(ip).EQ.nbp_lat) THEN
      tabijGCM(ip)=klon_glo
     else 
      print*,'ip jpt tabijGCM',ip,jpt(ip),tabijGCM(ip)
     endif
!    PRINT*,'CFMIP ip lon lat tabijGCM',ip,lonGCM(ip),latGCM(ip),tabijGCM(ip)
    ENDDO
    DO ip=1, npCFMIP
     if(lonGCM(ip).EQ.io_lon(1)) lonGCM(ip)=360.+lonGCM(ip)
    ENDDO
   DO i=1, npCFMIP
    WRITE(23,*) i, lonCFMIP(i), latCFMIP(i), lonGCM(i), latGCM(i), tabijGCM(i)
   ENDDO
   CLOSE(23)
 END SUBROUTINE LMDZ_CFMIP_point_locations 

END MODULE CFMIP_point_locations
