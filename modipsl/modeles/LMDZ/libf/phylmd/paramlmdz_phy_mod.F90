MODULE paramLMDZ_phy_mod 

! Olivier 13/07/2016
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

  SUBROUTINE ini_paramLMDZ_phy(dtime,nid_ctesGCM)

    USE iophy 
    USE dimphy
    USE ioipsl, only: histbeg, histvert, histdef, histend, ymds2ju
    USE mod_phys_lmdz_mpi_data, ONLY: is_mpi_root
    USE geometry_mod, ONLY: longitude_deg, latitude_deg
    USE mod_grid_phy_lmdz, only: nbp_lon, nbp_lat, klon_glo
    USE time_phylmdz_mod, ONLY: annee_ref, day_ref, itau_phy, pdtphys
    USE mod_phys_lmdz_transfert_para, ONLY: gather, bcast

    IMPLICIT NONE

    include "clesphys.h"
    include "YOMCST.h"

    REAL, INTENT(OUT)    :: dtime
    INTEGER, INTENT(OUT) :: nid_ctesGCM

    REAL,DIMENSION(klon_glo)        :: rlat_glo
    REAL,DIMENSION(klon_glo)        :: rlon_glo
    INTEGER i, idayref, ISW, itau_w
    REAL zstophy, zout
    REAL zx_lon(nbp_lon,nbp_lat)
    REAL zx_lat(nbp_lon,nbp_lat)

    CHARACTER*1 ch1
    INTEGER nhori
    INTEGER, PARAMETER :: np=1

    REAL zjulian
    SAVE zjulian
!$OMP THREADPRIVATE(zjulian)

!IM    Implemente en modes sequentiel et parallele

       CALL gather(latitude_deg,rlat_glo)
       CALL bcast(rlat_glo)
       CALL gather(longitude_deg,rlon_glo)
       CALL bcast(rlon_glo)

!$OMP MASTER
      IF (is_mpi_root) THEN
!
!       zstophy = pdtphys
!       zout = -1
!--OB modified for daily output
       zstophy = 86400.
       zout = 86400.
!
       idayref = day_ref
       CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
!
       CALL gr_fi_ecrit(1,klon,nbp_lon,nbp_lat,rlon_glo,zx_lon)
       IF (nbp_lon.GT.1) THEN
       DO i = 1, nbp_lon
         zx_lon(i,1) = rlon_glo(i+1)
         zx_lon(i,nbp_lat) = rlon_glo(i+1)
       ENDDO
       ENDIF
!
       CALL gr_fi_ecrit(1,klon,nbp_lon,nbp_lat,rlat_glo,zx_lat)
!
       CALL histbeg("paramLMDZ_phy.nc",  &
                       np,zx_lon(np:np,1), np,zx_lat(1,np:np), &
                       1,1,1,1, &
                       itau_phy, zjulian, dtime, &
                       nhori, nid_ctesGCM)
!
       CALL histdef(nid_ctesGCM, "R_ecc",  &
                      "Excentricite","-", &
                      1,1,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
       CALL histdef(nid_ctesGCM, "R_peri",  &
                      "Equinoxe","-", &
                      1,1,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
       CALL histdef(nid_ctesGCM, "R_incl",  &
                      "Inclinaison","deg", &
                      1,1,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
       CALL histdef(nid_ctesGCM, "solaire",  &
                      "Constante solaire","W/m2", &
                      1,1,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
#ifdef CPP_RRTM
       IF (iflag_rrtm.EQ.1) THEN
         DO ISW=1, NSW
           WRITE(ch1,'(i1)') ISW
           CALL histdef(nid_ctesGCM, "rsun"//ch1,  &
                      "Fraction constante solaire bande "//ch1,"W/m2", &
                      1,1,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
         ENDDO
       ENDIF
#endif
!
       CALL histdef(nid_ctesGCM, "co2_ppm",  &
                      "Concentration du CO2", "ppm", &
                      1,1,nhori, 1,1,1, -99, 32,  &
                      "ave(X)", zstophy,zout)
!
       CALL histdef(nid_ctesGCM, "CH4_ppb",  &
                      "Concentration du CH4", "ppb", &
                      1,1,nhori, 1,1,1, -99, 32,  &
                      "ave(X)", zstophy,zout)
!
       CALL histdef(nid_ctesGCM, "N2O_ppb", &
                      "Concentration du N2O", "ppb", &
                      1,1,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
       CALL histdef(nid_ctesGCM, "CFC11_ppt", &
                      "Concentration du CFC11", "ppt", &
                      1,1,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
       CALL histdef(nid_ctesGCM, "CFC12_ppt", &
                      "Concentration du CFC12", "ppt", &
                      1,1,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
       CALL histend(nid_ctesGCM)
       
       ENDIF !(is_mpi_root)
!$OMP END MASTER

  END SUBROUTINE ini_paramLMDZ_phy

!=================================================================

  SUBROUTINE write_paramLMDZ_phy(itap,nid_ctesGCM,ok_sync)

    USE mod_phys_lmdz_mpi_data, ONLY: is_mpi_root
    USE time_phylmdz_mod, ONLY: day_step_phy, annee_ref, itau_phy, start_time
    USE mod_grid_phy_lmdz, only: nbp_lon, nbp_lat, klon_glo

    USE iophy
    USE ioipsl, ONLY: histwrite, histsync

#ifdef CPP_RRTM
    USE YOESW, ONLY : RSUN
#endif

    IMPLICIT NONE

    include "clesphys.h"
    include "YOMCST.h"

    INTEGER, INTENT(IN) :: itap, nid_ctesGCM
    LOGICAL, INTENT(IN) :: ok_sync

    INTEGER itau_w, ISW
    INTEGER ndex2d(nbp_lon*nbp_lat)
    REAL :: zx_tmp_0d(1,1)
    INTEGER, PARAMETER :: np=1

    CHARACTER*1 ch1

!$OMP MASTER
      IF (is_mpi_root) THEN      
!
      ndex2d = 0
      itau_w = itau_phy + itap + int(start_time * day_step_phy)
!
! Variables globales
!
      zx_tmp_0d=R_ecc
      CALL histwrite(nid_ctesGCM,"R_ecc",itau_w, &
                     zx_tmp_0d,np,ndex2d)
!
      zx_tmp_0d=R_peri
      CALL histwrite(nid_ctesGCM,"R_peri",itau_w, &
                     zx_tmp_0d,np,ndex2d)
!
      zx_tmp_0d=R_incl
      CALL histwrite(nid_ctesGCM,"R_incl",itau_w, &
                     zx_tmp_0d,np,ndex2d)
!
      zx_tmp_0d=solaire
      CALL histwrite(nid_ctesGCM,"solaire",itau_w, &
                     zx_tmp_0d,np,ndex2d)
!
#ifdef CPP_RRTM
      IF (iflag_rrtm.EQ.1) THEN
        DO ISW=1, NSW
          WRITE(ch1,'(i1)') ISW
          zx_tmp_0d=RSUN(ISW)
          CALL histwrite(nid_ctesGCM,"rsun"//ch1,itau_w, &
                         zx_tmp_0d,np,ndex2d)
        ENDDO 
      ENDIF
#endif
!
      zx_tmp_0d=co2_ppm
      CALL histwrite(nid_ctesGCM,"co2_ppm",itau_w, &
                     zx_tmp_0d,np,ndex2d)
!
      zx_tmp_0d=CH4_ppb
      CALL histwrite(nid_ctesGCM,"CH4_ppb",itau_w, &
                     zx_tmp_0d,np,ndex2d)
!
      zx_tmp_0d=N2O_ppb
      CALL histwrite(nid_ctesGCM,"N2O_ppb",itau_w, &
                     zx_tmp_0d,np,ndex2d)
!
      zx_tmp_0d=CFC11_ppt
      CALL histwrite(nid_ctesGCM,"CFC11_ppt",itau_w, &
                     zx_tmp_0d,np,ndex2d)
!
      zx_tmp_0d=CFC12_ppt
      CALL histwrite(nid_ctesGCM,"CFC12_ppt",itau_w, &
                     zx_tmp_0d,np,ndex2d)
!
!=================================================================
!
      IF (ok_sync) THEN
        call histsync(nid_ctesGCM)
      ENDIF

      ENDIF !(is_mpi_root) then      
!$OMP END MASTER

  END SUBROUTINE write_paramLMDZ_phy

END MODULE paramLMDZ_phy_mod
