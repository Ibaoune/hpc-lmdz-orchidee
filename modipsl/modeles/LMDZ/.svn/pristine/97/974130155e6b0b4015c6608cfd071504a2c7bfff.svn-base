      SUBROUTINE condsurfc(jour,lmt_bcff,lmt_bcbb,
     .                     lmt_bcbbl,lmt_bcbbh,lmt_bc_penner,
     .                     lmt_omff,lmt_ombb,lmt_ombbl,lmt_ombbh,
     .                     lmt_omnat)
      USE dimphy
      IMPLICIT none
!
! Lire les conditions aux limites du modele pour la chimie.
! --------------------------------------------------------
!
#include "dimensions.h"
#include "netcdf.inc"
      REAL lmt_bcff(klon), lmt_bcbb(klon),lmt_bc_penner(klon)
      REAL lmt_omff(klon), lmt_ombb(klon)
      REAL lmt_bcbbl(klon), lmt_bcbbh(klon)
      REAL lmt_ombbl(klon), lmt_ombbh(klon)
      REAL lmt_omnat(klon)
      REAL lmt_terp(klon)
!
      INTEGER jour, i
      INTEGER ierr
      INTEGER nid1,nvarid
      INTEGER debut(2),epais(2)
!
      IF (jour.LT.0 .OR. jour.GT.(360-1)) THEN
         IF (jour.GT.(360-1).AND.jour.LE.367) THEN
           jour=360-1
           print *,'JE: jour changed to jour= ',jour
         ELSE
           PRINT*,'Le jour demande n est pas correcte:', jour
           CALL ABORT
         ENDIF
      ENDIF
!
      ierr = NF_OPEN ("limitcarbon.nc", NF_NOWRITE, nid1)
      if (ierr.ne.NF_NOERR) then
        write(6,*)' Pb d''ouverture du fichier limitbc.nc'
        write(6,*)' ierr = ', ierr
        call exit(1)
      endif
!
! Tranche a lire:
      debut(1) = 1
      debut(2) = jour+1
      epais(1) = klon
      epais(2) = 1
!
!
      ierr = NF_INQ_VARID (nid1, "BCFF", nvarid)
!nhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais, lmt_bcff)
!      print *,'IERR = ',ierr
!      print *,'NF_NOERR = ',NF_NOERR
!      print *,'debut = ',debut
!      print *,'epais = ',epais
!nhl #else
!nhl       ierr = NF_GET_VARA_REAL (nid1, nvarid, debut, epais, lmt_bcff)
!nhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources BC'
         CALL exit(1)
      ENDIF
!
!
      ierr = NF_INQ_VARID (nid1, "BCBB", nvarid)
!nhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais, lmt_bcbb)
!nhl #else
!nhl       ierr = NF_GET_VARA_REAL (nid1, nvarid, debut, epais, lmt_bcbb)
!nhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources BC-biomass'
         CALL exit(1)
      ENDIF
!
!
      ierr = NF_INQ_VARID (nid1, "BCBL", nvarid)
!nhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais, lmt_bcbbl)
!nhl #else
!nhl       ierr = NF_GET_VARA_REAL (nid1, nvarid, debut, epais, lmt_bcbbl)
!nhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources BC low'
         CALL exit(1)
      ENDIF
!
!
      ierr = NF_INQ_VARID (nid1, "BCBH", nvarid)
!nhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais, lmt_bcbbh)
!nhl #else
!nhl       ierr = NF_GET_VARA_REAL (nid1, nvarid, debut, epais, lmt_bcbbh)
!nhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources BC high'
         CALL exit(1)
      ENDIF
!
      ierr = NF_INQ_VARID (nid1, "TERP", nvarid)
!nhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais, lmt_terp)
!nhl #else
!nhl       ierr = NF_GET_VARA_REAL (nid1, nvarid, debut, epais, lmt_terp)
!nhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources Terpene'
         CALL exit(1)
      ENDIF
!
!
      ierr = NF_INQ_VARID (nid1, "BC_penner", nvarid)
!nhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, 
     .       epais, lmt_bc_penner)
!nhl #else
!nhl       ierr = NF_GET_VARA_REAL (nid1, nvarid, debut, epais, 
!nhl      .       lmt_bc_penner)
!nhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources BC Penner'
         CALL exit(1)
      ENDIF
!
!
      ierr = NF_INQ_VARID (nid1, "OMFF", nvarid)
!nhl #ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais, lmt_omff)
!nhl #else
!nhl       ierr = NF_GET_VARA_REAL (nid1, nvarid, debut, epais, lmt_omff)
!nhl #endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources om-ifossil'
         CALL exit(1)
      ENDIF
!
      DO i=1,klon
        lmt_ombb(i)  = lmt_bcbb(i)*7.0*1.6      !OC/BC=7.0;OM/OC=1.6
        lmt_ombbl(i) = lmt_bcbbl(i)*7.0*1.6
        lmt_ombbh(i) = lmt_bcbbh(i)*7.0*1.6
        lmt_omff(i)  = lmt_omff(i)*1.4          !--OM/OC=1.4
        lmt_omnat(i)  = lmt_terp(i)*0.11*1.4 !-- 11% Terpene is OC
      ENDDO
!
      ierr = NF_CLOSE(nid1)
      PRINT*, 'Carbon sources lues pour jour: ', jour
!
      RETURN
      END
