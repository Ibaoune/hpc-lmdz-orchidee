      SUBROUTINE condsurfs(jour, edgar, flag_dms,
     .                     lmt_so2h, lmt_so2b, lmt_so2bb, lmt_so2ba,
     .                     lmt_so2volc, lmt_altvolc,  
     .                     lmt_dmsbio, lmt_h2sbio, lmt_dms, lmt_dmsconc)
       USE dimphy
      IMPLICIT none
c
c Lire les conditions aux limites du modele pour la chimie.
c --------------------------------------------------------
c
#include "dimensions.h"
#include "netcdf.inc"
c
      REAL lmt_so2h(klon), lmt_so2b(klon), lmt_so2bb(klon)
      REAL lmt_dmsbio(klon), lmt_h2sbio(klon), lmt_so2ba(klon)
      REAL lmt_so2volc(klon), lmt_altvolc(klon)
      REAL lmt_dms(klon), lmt_dmsconc(klon)
      LOGICAL edgar
      INTEGER flag_dms
c
      INTEGER jour, i
      INTEGER ierr
      INTEGER nid,nvarid
      INTEGER debut(2),epais(2)
c
      IF (jour.LT.0 .OR. jour.GT.(360-1)) THEN
         IF ((jour.GT.(360-1)) .AND. (jour.LE.367)) THEN
           jour=360-1
           print *,'JE: jour changed to jour= ',jour
         ELSE
           PRINT*,'Le jour demande n est pas correcte:', jour
           CALL ABORT
         ENDIF
      ENDIF
c
      ierr = NF_OPEN ("limitsoufre.nc", NF_NOWRITE, nid)
      if (ierr.ne.NF_NOERR) then
        write(6,*)' Pb d''ouverture du fichier limitsoufre.nc'
        write(6,*)' ierr = ', ierr
        call exit(1)
      endif
c
c Tranche a lire:
      debut(1) = 1
      debut(2) = jour+1
      epais(1) = klon
      epais(2) = 1
c
      ierr = NF_INQ_VARID (nid, "VOLC", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_so2volc)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_so2volc)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources so2 volcan'
         CALL exit(1)
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "ALTI", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_altvolc)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_altvolc)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les altitudes volcan'
         CALL exit(1)
      ENDIF
c
      IF (edgar) THEN   !--EDGAR w/o ship and biomass burning
c
      ierr = NF_INQ_VARID (nid, "SO2ED95L", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_so2b)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_so2b)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources so2 edgar low'
         CALL exit(1)
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "SO2ED95H", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_so2h)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_so2h)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources so2 edgar high'
         CALL exit(1)
      ENDIF
c
      ELSE  !--GEIA 
c
      ierr = NF_INQ_VARID (nid, "SO2H", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_so2h)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_so2h)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources so2 haut'
         CALL exit(1)
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "SO2B", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_so2b)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_so2b)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources so2 bas'
         CALL exit(1)
      ENDIF
c
      ENDIF  !--edgar
c
      ierr = NF_INQ_VARID (nid, "SO2BB", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_so2bb)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_so2bb)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources so2 bb'
         CALL exit(1)
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "SO2BA", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_so2ba)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_so2ba)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources so2 bateau'
         CALL exit(1)
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "DMSB", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_dmsbio)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dmsbio)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms bio'
         CALL exit(1)
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "H2SB", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_h2sbio)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_h2sbio)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources h2s bio'
         CALL exit(1)
      ENDIF
c
      IF (flag_dms.EQ.1) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSL", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_dms)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dms)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms liss'
         CALL exit(1)
      ENDIF
c
      ELSEIF (flag_dms.EQ.2) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSW", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, lmt_dms)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dms)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms wann'
         CALL exit(1)
      ENDIF
c
      ELSEIF (flag_dms.EQ.3) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSC1", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut,epais,lmt_dmsconc)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dmsconc)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dmsconc old'
         CALL exit(1)
      ENDIF
c
      ELSEIF (flag_dms.EQ.4) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSC2", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut,epais,lmt_dmsconc)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dmsconc)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms conc 2'
         CALL exit(1)
      ENDIF
c
      ELSEIF (flag_dms.EQ.5) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSC3", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut,epais,lmt_dmsconc)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dmsconc)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms conc 3'
         CALL exit(1)
      ENDIF
c
      ELSEIF (flag_dms.EQ.6) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSC4", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut,epais,lmt_dmsconc)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dmsconc)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms conc 4'
         CALL exit(1)
      ENDIF
c
      ELSEIF (flag_dms.EQ.7) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSC5", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut,epais,lmt_dmsconc)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dmsconc)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms conc 5'
         CALL exit(1)
      ENDIF
c
      ELSEIF (flag_dms.EQ.8) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSC6", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut,epais,lmt_dmsconc)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dmsconc)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms conc 6'
         CALL exit(1)
      ENDIF
c
      ELSEIF (flag_dms.EQ.9) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSC7", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut,epais,lmt_dmsconc)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dmsconc)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms conc 7'
         CALL exit(1)
      ENDIF
c
      ELSEIF (flag_dms.EQ.10) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSC8", nvarid)
cnhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut,epais,lmt_dmsconc)
cnhl #else
cnhl       ierr = NF_GET_VARA_REAL (nid, nvarid, debut, epais, lmt_dmsconc)
cnhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms conc 8'
         CALL exit(1)
      ENDIF
c
      ELSE 
c
         PRINT *,'choix non possible pour flag_dms'
         STOP
c
      ENDIF
c
      ierr = NF_CLOSE(nid)
c
      IF (flag_dms.LE.2) THEN 
      DO i=1, klon 
         lmt_dmsconc(i)=0.0
      ENDDO
      ELSE 
      DO i=1, klon 
         lmt_dms(i)=0.0
      ENDDO
      ENDIF
c
      PRINT*, 'Sources SOUFRE lues pour jour: ', jour
c
      RETURN
      END
