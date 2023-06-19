      SUBROUTINE condsurfs_new(jour, edgar, flag_dms,
     .                         lmt_so2b, lmt_so2h, lmt_so2nff,
     .                         lmt_so2bb_l, lmt_so2bb_h, lmt_so2ba,
     .                         lmt_so2volc_cont, lmt_altvolc_cont,  
     .                         lmt_so2volc_expl, lmt_altvolc_expl,  
     .                         lmt_dmsbio, lmt_h2sbio, lmt_dms, 
     .                                                      lmt_dmsconc)
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
c
      REAL lmt_so2b(klon), lmt_so2h(klon), lmt_so2nff(klon)
      REAL lmt_so2bb_l(klon), lmt_so2bb_h(klon)
      REAL lmt_dmsbio(klon), lmt_h2sbio(klon), lmt_so2ba(klon)
      REAL lmt_so2volc_cont(klon), lmt_altvolc_cont(klon)
      REAL lmt_so2volc_expl(klon), lmt_altvolc_expl(klon)
      REAL lmt_dms(klon), lmt_dmsconc(klon)

      REAL lmt_so2b_glo(klon_glo), lmt_so2h_glo(klon_glo)
      REAL lmt_so2nff_glo(klon_glo)
      REAL lmt_so2bb_l_glo(klon_glo), lmt_so2bb_h_glo(klon_glo)
      REAL lmt_dmsbio_glo(klon_glo), lmt_h2sbio_glo(klon_glo)
      REAL lmt_so2ba_glo(klon_glo)
      REAL lmt_so2volc_cont_glo(klon_glo),lmt_altvolc_cont_glo(klon_glo)
      REAL lmt_so2volc_expl_glo(klon_glo),lmt_altvolc_expl_glo(klon_glo)
      REAL lmt_dms_glo(klon_glo), lmt_dmsconc_glo(klon_glo)
      LOGICAL edgar
      INTEGER flag_dms
c
      INTEGER jour, i
      INTEGER ierr
      INTEGER nid,nvarid
      INTEGER debut(2),epais(2)
c
      IF (jour.LT.0 .OR. jour.GT.(366-1)) THEN
         PRINT*,'Le jour demande n est pas correcte:', jour
         print *,'JE: FORCED TO CONTINUE (emissions have
     . to be longer than 1 year!!!! )'
!         CALL ABORT
      ENDIF
!

!$OMP MASTER
      IF (is_mpi_root .AND. is_omp_root) THEN

c Tranche a lire:
      debut(1) = 1
      debut(2) = jour
!      epais(1) = klon
      epais(1) = klon_glo
      epais(2) = 1
!=======================================================================
!                 READING NEW EMISSIONS FROM RCP
!=======================================================================
!
      ierr = NF_OPEN ("sulphur_emissions_antro.nc", NF_NOWRITE, nid)
      if (ierr.ne.NF_NOERR) then
        write(6,*)' Pb d''ouverture du fichier sulphur_emissions_antro'
        write(6,*)' ierr = ', ierr
        call exit(1)
      endif

!
! SO2 Low level emissions
!
      ierr = NF_INQ_VARID (nid, "SO2FF_LOW", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais,lmt_so2b_glo)
      IF (ierr .NE. NF_NOERR) THEN
        PRINT*, 'Pb de lecture pour les sources so2 low'
        print *,'JE klon, jour, debut ,epais ',klon_glo,jour,debut,epais
        CALL HANDLE_ERR(ierr)
        print *,'error ierr= ',ierr
        CALL exit(1)
      ENDIF
!
! SO2 High level emissions
!
      ierr = NF_INQ_VARID (nid, "SO2FF_HIGH", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais,lmt_so2h_glo)
      IF (ierr .NE. NF_NOERR) THEN
        PRINT*, 'Pb de lecture pour les sources so2 high'
        CALL exit(1)
      ENDIF
!
! SO2 Biomass burning High level emissions
!
      ierr = NF_INQ_VARID (nid, "SO2BBH", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, 
     . epais, lmt_so2bb_h_glo)
      IF (ierr .NE. NF_NOERR) THEN
        PRINT*, 'Pb de lecture pour les sources so2 BB high'
        CALL exit(1)
      ENDIF
!
! SO2 biomass burning low level emissions
!
      ierr = NF_INQ_VARID (nid, "SO2BBL", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, 
     . epais, lmt_so2bb_l_glo)
      IF (ierr .NE. NF_NOERR) THEN
        PRINT*, 'Pb de lecture pour les sources so2 BB low'
        CALL exit(1)
      ENDIF
!
! SO2 ship emissions
!
      ierr = NF_INQ_VARID (nid, "SO2BA", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut,epais,lmt_so2ba_glo)
      IF (ierr .NE. NF_NOERR) THEN
        PRINT*, 'Pb de lecture pour les sources so2 ship'
        CALL exit(1)
      ENDIF
!
! SO2 Non Fossil Fuel Emissions
!
      ierr = NF_INQ_VARID (nid, "SO2NFF", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais,
     .  lmt_so2nff_glo)
      IF (ierr .NE. NF_NOERR) THEN
        PRINT*, 'Pb de lecture pour les sources so2 non FF'
        CALL exit(1)
      ENDIF
!
      ierr = NF_CLOSE(nid)
!
!=======================================================================
!                      READING NATURAL EMISSIONS
!=======================================================================
      ierr = NF_OPEN ("sulphur_emissions_nat.nc", NF_NOWRITE, nid)
      if (ierr.ne.NF_NOERR) then
        write(6,*)' Pb d''ouverture du fichier sulphur_emissions_nat'
        write(6,*)' ierr = ', ierr
        call exit(1)
      endif
c
c Biologenic source of DMS
c
      ierr = NF_INQ_VARID (nid, "DMSB", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid,debut,epais,lmt_dmsbio_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms bio'
         CALL exit(1)
      ENDIF
c
c Biologenic source of H2S
c
      ierr = NF_INQ_VARID (nid, "H2SB", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid,debut,epais,lmt_h2sbio_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources h2s bio'
         CALL exit(1)
      ENDIF
c
c Ocean surface concentration of dms (emissions are computed later)
c
      IF (flag_dms.EQ.4) THEN 
c
      ierr = NF_INQ_VARID (nid, "DMSC2", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid,nvarid,debut,epais,lmt_dmsconc_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources dms conc 2'
         CALL exit(1)
      ENDIF
c
      DO i=1, klon 
!         lmt_dms(i)=0.0
         lmt_dms_glo(i)=0.0
      ENDDO
c
      ELSE 
c
         PRINT *,'choix non possible pour flag_dms'
         STOP

      ENDIF
c
      ierr = NF_CLOSE(nid)
c
!=======================================================================
!                      READING VOLCANIC EMISSIONS
!=======================================================================
      print *,'   ***      READING VOLCANIC EMISSIONS   ***   '
      print *,' Jour = ',jour
      ierr = NF_OPEN ("sulphur_emissions_volc.nc", NF_NOWRITE, nid)
      if (ierr.ne.NF_NOERR) then
        write(6,*)' Pb d''ouverture du fichier sulphur_emissions_volc'
        write(6,*)' ierr = ', ierr
        call exit(1)
      endif
c
c Continuous Volcanic emissions
c
!      ierr = NF_INQ_VARID (nid, "VOLC", nvarid)
      ierr = NF_INQ_VARID (nid, "flx_volc_cont", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, 
     .                           lmt_so2volc_cont_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources so2 volcan (cont)'
         CALL exit(1)
      ENDIF
      print *,'SO2 volc cont (in read) = ',SUM(lmt_so2volc_cont_glo),
     +      MINVAL(lmt_so2volc_cont_glo),MAXVAL(lmt_so2volc_cont_glo)
!      lmt_so2volc(:)=0.0
c
c Altitud of continuous volcanic emissions
c
!      ierr = NF_INQ_VARID (nid, "ALTI", nvarid)
      ierr = NF_INQ_VARID (nid, "flx_volc_altcont", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, 
     .                           lmt_altvolc_cont_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les altitudes volcan (cont)'
         CALL exit(1)
      ENDIF
c
c Explosive Volcanic emissions
c
      ierr = NF_INQ_VARID (nid, "flx_volc_expl", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, 
     .                           lmt_so2volc_expl_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les sources so2 volcan (expl)'
         CALL exit(1)
      ENDIF
!      lmt_so2volc_expl(:)=0.0
      print *,'SO2 volc expl (in read) = ',SUM(lmt_so2volc_expl_glo),
     +      MINVAL(lmt_so2volc_expl_glo),MAXVAL(lmt_so2volc_expl_glo)
c
c Altitud of explosive volcanic emissions
c
      ierr = NF_INQ_VARID (nid, "flx_volc_altexpl", nvarid)
      ierr = NF_GET_VARA_DOUBLE (nid, nvarid, debut, epais, 
     .                           lmt_altvolc_expl_glo)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, 'Pb de lecture pour les altitudes volcan'
         CALL exit(1)
      ENDIF
!      lmt_altvolc_expl(:)=0.0

      ierr = NF_CLOSE(nid)
c
      PRINT*, 'Sources SOUFRE lues pour jour: ', jour
c


      ENDIF
!$OMP END MASTER
!$OMP BARRIER
      call scatter( lmt_so2b_glo        , lmt_so2b )
      call scatter(lmt_so2h_glo         , lmt_so2h )  
      call scatter(lmt_so2bb_h_glo      , lmt_so2bb_h )
      call scatter(lmt_so2bb_l_glo      , lmt_so2bb_l)
      call scatter(lmt_so2ba_glo        , lmt_so2ba)
      call scatter(lmt_so2nff_glo       , lmt_so2nff)
      call scatter(lmt_dmsbio_glo       , lmt_dmsbio)
      call scatter(lmt_h2sbio_glo       , lmt_h2sbio)
      call scatter(lmt_dmsconc_glo      , lmt_dmsconc)
      call scatter(lmt_dms_glo          , lmt_dms)
      call scatter(lmt_so2volc_cont_glo , lmt_so2volc_cont)
      call scatter(lmt_altvolc_cont_glo , lmt_altvolc_cont)
      call scatter(lmt_so2volc_expl_glo , lmt_so2volc_expl)
      call scatter(lmt_altvolc_expl_glo , lmt_altvolc_expl)


      RETURN
      END
