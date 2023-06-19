      SUBROUTINE condsurfc_new(jour,lmt_bcff, lmt_bcnff,
     .                         lmt_bcbbl,lmt_bcbbh, lmt_bcba,
     .                         lmt_omff,lmt_omnff,lmt_ombbl,lmt_ombbh,
     .                                             lmt_omnat, lmt_omba)
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      USE dimphy
      IMPLICIT none
c
c Lire les conditions aux limites du modele pour la chimie.
c --------------------------------------------------------
c
#include "dimensions.h"
#include "netcdf.inc"
      REAL lmt_bcff(klon), lmt_bcnff(klon), lmt_bcba(klon)
      REAL lmt_omff(klon), lmt_omnff(klon), lmt_ombb(klon)
      REAL lmt_bcbbl(klon), lmt_bcbbh(klon)
      REAL lmt_ombbl(klon), lmt_ombbh(klon)
      REAL lmt_omnat(klon), lmt_omba(klon)
      REAL lmt_terp(klon)
c
      REAL lmt_bcff_glo(klon_glo), lmt_bcnff_glo(klon_glo)
      REAL lmt_bcba_glo(klon_glo)
      REAL lmt_omff_glo(klon_glo), lmt_omnff_glo(klon_glo)
      REAL lmt_ombb_glo(klon_glo)
      REAL lmt_bcbbl_glo(klon_glo), lmt_bcbbh_glo(klon_glo)
      REAL lmt_ombbl_glo(klon_glo), lmt_ombbh_glo(klon_glo)
      REAL lmt_omnat_glo(klon_glo), lmt_omba_glo(klon_glo)
      REAL lmt_terp_glo(klon_glo)
!
      INTEGER jour, i
      INTEGER ierr
      INTEGER nid1,nvarid
      INTEGER debut(2),epais(2)
c
!      IF (jour.LT.0 .OR. jour.GT.(366-1)) THEN
      IF (jour.LT.0 .OR. jour.GT.366) THEN
         PRINT*,'Le jour demande n est pas correcte:', jour
         print *,'JE: FORCED TO CONTINUE (emissions have
     . to be longer than 1 year!!!! )'
!JE         CALL ABORT
      ENDIF

!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN
!
! Tranche a lire:
      debut(1) = 1
      debut(2) = jour
      epais(1) = klon_glo
!      epais(1) = klon
      epais(2) = 1
!
!=======================================================================
!                        BC EMISSIONS
!=======================================================================
!
      ierr = NF_OPEN ("carbon_emissions.nc", NF_NOWRITE, nid1)
      if (ierr.ne.NF_NOERR) then
        write(6,*)' Pb d''ouverture du fichier limitbc.nc'
        write(6,*)' ierr = ', ierr
        call exit(1)
      endif
!
! BC emissions from fossil fuel combustion
!
      ierr = NF_INQ_VARID (nid1, "BCFF", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais,
     .  lmt_bcff_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources BC'
         CALL exit(1)
      ENDIF
      !print *,'lmt_bcff = ',lmt_bcff
      !stop
!
! BC emissions from non fossil fuel combustion
!
      ierr = NF_INQ_VARID (nid1, "BCNFF", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais,
     .    lmt_bcnff_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources BC'
         CALL exit(1)
      ENDIF
!
! Low BC emissions from biomass burning
!
      ierr = NF_INQ_VARID (nid1, "BCBBL", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais,
     .  lmt_bcbbl_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources BC low'
         CALL exit(1)
      ENDIF
!
! High BC emissions from biomass burning
!
      ierr = NF_INQ_VARID (nid1, "BCBBH", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais,
     .      lmt_bcbbh_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources BC high'
         CALL exit(1)
      ENDIF
!
! BC emissions from ship transport
!
      ierr = NF_INQ_VARID (nid1, "BCBA", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais, 
     .   lmt_bcba_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources BC'
         CALL exit(1)
      ENDIF
!
!=======================================================================
!                        OM EMISSIONS 
!=======================================================================
!

!
! OM emissions from fossil fuel combustion
!
      ierr = NF_INQ_VARID (nid1, "OMFF", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais,
     .  lmt_omff_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources OM'
         CALL exit(1)
      ENDIF
!
! OM emissions from non fossil fuel combustion
!
      ierr = NF_INQ_VARID (nid1, "OMNFF", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais,
     .   lmt_omnff_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources OM'
         CALL exit(1)
      ENDIF
!
! Low OM emissions from biomass burning - low
!
      ierr = NF_INQ_VARID (nid1, "OMBBL", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais, 
     .  lmt_ombbl_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources OM low'
         CALL exit(1)
      ENDIF
!
! High OM emissions from biomass burning - high
!
      ierr = NF_INQ_VARID (nid1, "OMBBH", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais,
     .  lmt_ombbh_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources OM high'
         CALL exit(1)
      ENDIF
!
! High OM emissions from ship
!
      ierr = NF_INQ_VARID (nid1, "OMBA", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais,
     .   lmt_omba_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources OM ship'
         CALL exit(1)
      ENDIF
!
! Natural Terpene emissions => Natural OM emissions
!
      ierr = NF_INQ_VARID (nid1, "TERP", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid1, nvarid, debut, epais,
     .  lmt_terp_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources Terpene'
         CALL exit(1)
      ENDIF
!
      DO i=1,klon_glo
        lmt_omnat_glo(i)  = lmt_terp_glo(i)*0.11*1.4 !-- 11% Terpene is OC
      ENDDO

      ierr = NF_CLOSE(nid1)
!
      PRINT*, 'Carbon sources lues pour jour: ', jour
!      lmt_bcff(klon)=0.0
!      lmt_bcnff(klon)=0.0
!      lmt_omff(klon)=0.0
!      lmt_omnff(klon)=0.0
!      lmt_ombb(klon)=0.0
!      lmt_bcbbl(klon)=0.0
!      lmt_bcbbh(klon)=0.0
!      lmt_ombbl(klon)=0.0
!      lmt_ombbh(klon)=0.0
!      lmt_omnat(klon)=0.0
!      lmt_omba(klon)=0.0
!      lmt_terp(klon)=0.0


      ENDIF
!$OMP END MASTER
!$OMP BARRIER
      call scatter( lmt_bcff_glo   , lmt_bcff )   
      call scatter( lmt_bcnff_glo  , lmt_bcnff )
      call scatter( lmt_bcbbl_glo  , lmt_bcbbl )
      call scatter( lmt_bcbbh_glo  , lmt_bcbbh )
      call scatter( lmt_bcba_glo   , lmt_bcba )
      call scatter( lmt_omff_glo   , lmt_omff )
      call scatter( lmt_omnff_glo  , lmt_omnff )
      call scatter( lmt_ombbl_glo  , lmt_ombbl )
      call scatter( lmt_ombbh_glo  , lmt_ombbh )
      call scatter( lmt_omba_glo   , lmt_omba )
      call scatter( lmt_terp_glo   , lmt_terp )
      call scatter( lmt_omnat_glo  , lmt_omnat )





      RETURN
      END
