c This subroutine calculates the emissions of SEA SALT and DUST, part of
C which goes to tracer 2 and other part to tracer 3.
      SUBROUTINE coarsemission(pctsrf,pdtphys,
     .                         t_seri,pmflxr,pmflxs,prfl,psfl,
     .                         xlat,xlon,debutphy,
     .                         zu10m,zv10m,wstar,ale_bl,ale_wake,
     .                         scale_param_ssacc,scale_param_sscoa,
     .                         scale_param_dustacc,scale_param_dustcoa,
     .                         scale_param_dustsco,
     .                         nbreg_dust,
     .                         iregion_dust,dust_ec,
     .                  param_wstarBLperregion,param_wstarWAKEperregion,
     .                  nbreg_wstardust,
     .                  iregion_wstardust,
     .                         lmt_sea_salt,qmin,qmax,
     .                             flux_sparam_ddfine,flux_sparam_ddcoa,
     .                             flux_sparam_ddsco,
     .                             flux_sparam_ssfine,flux_sparam_sscoa,
     .                          id_prec,id_fine,id_coss,id_codu,id_scdu,
     .                          ok_chimeredust,
     .                                                source_tr,flux_tr)
!     .                         wth,cly,zprecipinsoil,lmt_sea_salt,

!      CALL dustemission( debutphy, xlat, xlon, pctsrf, 
!     .               zu10m     zv10m,wstar,ale_bl,ale_wake)
!

      USE dimphy
      USE indice_sol_mod
      USE infotrac
      USE dustemission_mod,  ONLY : dustemission 
!      USE phytracr_spl_mod, ONLY : nbreg_dust, nbreg_ind, nbreg_bb
      IMPLICIT NONE

#include "dimensions.h"
#include "chem.h"
#include "chem_spla.h"
#include "YOMCST.h"
#include "paramet.h"
      
c============================== INPUT ==================================
      INTEGER nbjour
      LOGICAL ok_chimeredust
      REAL pdtphys  ! pas d'integration pour la physique (seconde)
      REAL t_seri(klon,klev)  ! temperature
      REAL pctsrf(klon,nbsrf)
      REAL pmflxr(klon,klev+1), pmflxs(klon,klev+1)   !--convection
!      REAL pmflxr(klon,klev), pmflxs(klon,klev)   !--convection
      REAL prfl(klon,klev+1),   psfl(klon,klev+1)     !--large-scale
!      REAL prfl(klon,klev),   psfl(klon,klev)     !--large-scale
      LOGICAL debutphy, lafinphy
      REAL, intent(in) ::  xlat(klon)    ! latitudes pour chaque point 
      REAL, intent(in) ::  xlon(klon)    ! longitudes pour chaque point 
      REAL,DIMENSION(klon),INTENT(IN)    :: zu10m
      REAL,DIMENSION(klon),INTENT(IN)    :: zv10m
      REAL,DIMENSION(klon),INTENT(IN)    :: wstar,Ale_bl,ale_wake

c
c------------------------- Scaling Parameters --------------------------
c
      INTEGER iregion_dust(klon) !Defines  dust regions
      REAL scale_param_ssacc  !Scaling parameter for Fine Sea Salt
      REAL scale_param_sscoa  !Scaling parameter for Coarse Sea Salt
      REAL scale_param_dustacc(nbreg_dust)  !Scaling parameter for Fine Dust
      REAL scale_param_dustcoa(nbreg_dust)  !Scaling parameter for Coarse Dust
      REAL scale_param_dustsco(nbreg_dust)  !Scaling parameter for SCoarse Dust
!JE20141124<<
      INTEGER iregion_wstardust(klon) !Defines dust regions in terms of wstar
      REAL param_wstarBLperregion(nbreg_wstardust)  !
      REAL param_wstarWAKEperregion(nbreg_wstardust)  !
      REAL param_wstarBL(klon)  !parameter for surface wind correction..
      REAL param_wstarWAKE(klon)  !parameter for surface wind correction..
      INTEGER  nbreg_wstardust
!JE20141124>>
      INTEGER  nbreg_dust
      INTEGER, INTENT(IN) :: id_prec,id_fine,id_coss,id_codu,id_scdu
c============================== OUTPUT =================================
      REAL source_tr(klon,nbtr)
      REAL flux_tr(klon,nbtr)
      REAL flux_sparam_ddfine(klon), flux_sparam_ddcoa(klon)
      REAL flux_sparam_ddsco(klon)
      REAL flux_sparam_ssfine(klon), flux_sparam_sscoa(klon)
c=========================== LOCAL VARIABLES ===========================            
      INTEGER i, j
      REAL pct_ocean(klon)
!      REAL zprecipinsoil(klon)
!      REAL cly(klon), wth(klon)
      REAL clyfac, avgdryrate, drying
            
c---------------------------- SEA SALT emissions ------------------------
      REAL lmt_sea_salt(klon,ss_bins) !Sea salt 0.03-8.0 um                         
c
c--------vent 10 m CEPMMT
c
      REAL dust_ec(klon)

      real tmp_var2(klon,nbtr) ! auxiliary variable to replace source
      REAL qmin, qmax
!----------------------DUST Sahara ---------------
      REAL, DIMENSION(klon) :: dustsourceacc,dustsourcecoa,dustsourcesco
      INTEGER, DIMENSION(klon) :: maskd
C*********************** DUST EMMISSIONS *******************************
c
      
!     avgdryrate=300./365.*pdtphys/86400.
c
!     DO i=1, klon
c
!       IF (cly(i).LT.9990..AND.wth(i).LT.9990.) THEN
!        zprecipinsoil(i)=zprecipinsoil(i) +
!    .        (pmflxr(i,1)+pmflxs(i,1)+prfl(i,1)+psfl(i,1))*pdtphys
c
!        clyfac=MIN(16., cly(i)*0.4+8.) ![mm] max amount of water hold in top soil
!        drying=avgdryrate*exp(0.03905491*
!    .                    exp(0.17446*(t_seri(i,1)-273.15))) ! [mm]
!        zprecipinsoil(i)=min(max(0.,zprecipinsoil(i)-drying),clyfac) ! [mm]
c         
!       ENDIF
c
!     ENDDO
c                
c ==================== CALCULATING DUST EMISSIONS ======================
c
!      IF (lminmax) THEN
      DO j=1,nbtr
      DO i=1,klon
         tmp_var2(i,j)=source_tr(i,j)
      ENDDO
      ENDDO
      CALL minmaxsource(tmp_var2,qmin,qmax,'src: before DD emiss')
!      print *,'Source = ',SUM(source_tr),MINVAL(source_tr),
!     .                                     MAXVAL(source_tr)
!      ENDIF

c
      IF (.NOT. ok_chimeredust)  THEN  
      DO i=1, klon
!!     IF (cly(i).GE.9990..OR.wth(i).GE.9990..OR.
!!    .    t_seri(i,1).LE.273.15.OR.zprecipinsoil(i).GT.1.e-8) THEN
!!          dust_ec(i)=0.0
!!     ENDIF
!c Corresponds to dust_emission.EQ.3        
!!!!!!!****************AQUIIIIIIIIIIIIIIIIIIIIIIIIIIII
!! Original line (4 tracers)
!JE<<  old 4 tracer(nhl scheme)        source_tr(i,id_fine)=scale_param_dustacc(iregion_dust(i))*
!     .                  dust_ec(i)*1.e3*0.093   ! g/m2/s
!         source_tr(i,id_codu)=scale_param_dustcoa(iregion_dust(i))*
!     .                  dust_ec(i)*1.e3*0.905   ! g/m2/s   bin 0.5-10um
!! Original line (4 tracers)
!         flux_tr(i,id_fine)=scale_param_dustacc(iregion_dust(i))*
!     .                  dust_ec(i)*1.e3*0.093*1.e3  !mg/m2/s
!         flux_tr(i,id_codu)=scale_param_dustcoa(iregion_dust(i))*
!     .                  dust_ec(i)*1.e3*0.905*1.e3  !mg/m2/s bin 0.5-10um
!         flux_sparam_ddfine(i)=scale_param_dustacc(iregion_dust(i)) *
!     .                            dust_ec(i)*1.e3*0.093*1.e3
!         flux_sparam_ddcoa(i)=scale_param_dustcoa(iregion_dust(i)) *
!     .                            dust_ec(i)*1.e3*0.905*1.e3
      IF(id_fine>0)     source_tr(i,id_fine)=
     . scale_param_dustacc(iregion_dust(i))*
     .                  dust_ec(i)*1.e3*0.093   ! g/m2/s 
      IF(id_codu>0)   source_tr(i,id_codu)=
     . scale_param_dustcoa(iregion_dust(i))*
     .                  dust_ec(i)*1.e3*0.905   ! g/m2/s   bin 0.5-10um
      IF(id_scdu>0)  source_tr(i,id_scdu)=0.   ! no supercoarse
! Original line (4 tracers)
       IF(id_fine>0)   flux_tr(i,id_fine)=
     .  scale_param_dustacc(iregion_dust(i))*
     .                  dust_ec(i)*1.e3*0.093*1.e3  !mg/m2/s
       IF(id_codu>0)  flux_tr(i,id_codu)=
     . scale_param_dustcoa(iregion_dust(i))*
     .                  dust_ec(i)*1.e3*0.905*1.e3  !mg/m2/s bin 0.5-10um
       IF(id_scdu>0) flux_tr(i,id_scdu)=0.

         flux_sparam_ddfine(i)=scale_param_dustacc(iregion_dust(i)) *
     .                            dust_ec(i)*1.e3*0.093*1.e3
         flux_sparam_ddcoa(i)=scale_param_dustcoa(iregion_dust(i)) *
     .                            dust_ec(i)*1.e3*0.905*1.e3
         flux_sparam_ddsco(i)=0.
      ENDDO
      ENDIF
!*****************NEW CHIMERE DUST EMISSION Sahara*****
! je  20140522
      IF(ok_chimeredust) THEN
      print *,'MIX- NEW SAHARA DUST SOURCE SCHEME...'

      DO i=1,klon
      param_wstarBL(i)  =param_wstarBLperregion(iregion_wstardust(i))
      param_wstarWAKE(i)=param_wstarWAKEperregion(iregion_wstardust(i))
      ENDDO


      CALL dustemission( debutphy, xlat, xlon, pctsrf,
     .                  zu10m,zv10m,wstar,ale_bl,ale_wake,
     .                  param_wstarBL, param_wstarWAKE,
     .                  dustsourceacc,dustsourcecoa,
     .                  dustsourcesco,maskd)
      
      DO i=1,klon   
         if (maskd(i).gt.0) then 
      IF(id_fine>0)    source_tr(i,id_fine)=
     . scale_param_dustacc(iregion_dust(i))*
     .                  dustsourceacc(i)*1.e3   ! g/m2/s  bin 0.03-0.5
      IF(id_codu>0)    source_tr(i,id_codu)=
     . scale_param_dustcoa(iregion_dust(i))*
     .                  dustsourcecoa(i)*1.e3   ! g/m2/s   bin 0.5-3um
      IF(id_scdu>0)   source_tr(i,id_scdu)=
     . scale_param_dustsco(iregion_dust(i))*
     .                  dustsourcesco(i)*1.e3   ! g/m2/s   bin 3-15um
! Original line (4 tracers)
       IF(id_fine>0)  flux_tr(i,id_fine)=
     .  scale_param_dustacc(iregion_dust(i))*
     .                  dustsourceacc(i)*1.e3*1.e3  !mg/m2/s
       IF(id_codu>0)  flux_tr(i,id_codu)=
     . scale_param_dustcoa(iregion_dust(i))*
     .                  dustsourcecoa(i)*1.e3*1.e3  !mg/m2/s bin 0.5-3um
       IF(id_scdu>0)  flux_tr(i,id_scdu)=
     . scale_param_dustsco(iregion_dust(i))*
     .                  dustsourcesco(i)*1.e3*1.e3  !mg/m2/s bin 3-15um
         flux_sparam_ddfine(i)=scale_param_dustacc(iregion_dust(i)) *
     .                            dustsourceacc(i)*1.e3*1.e3
         flux_sparam_ddcoa(i)=scale_param_dustcoa(iregion_dust(i)) *
     .                            dustsourcecoa(i)*1.e3*1.e3
         flux_sparam_ddsco(i)=scale_param_dustsco(iregion_dust(i)) *
     .                            dustsourcesco(i)*1.e3*1.e3
         else
        IF(id_fine>0) source_tr(i,id_fine)=
     .  scale_param_dustacc(iregion_dust(i))*
     .                  dust_ec(i)*1.e3*0.114   ! g/m2/s
        IF(id_codu>0) source_tr(i,id_codu)=
     .  scale_param_dustcoa(iregion_dust(i))*
     .                  dust_ec(i)*1.e3*0.108   ! g/m2/s   bin 0.5-3um
        IF(id_scdu>0) source_tr(i,id_scdu)=
     .  scale_param_dustsco(iregion_dust(i))*
     .                  dust_ec(i)*1.e3*0.778   ! g/m2/s   bin 3-15um
! Original line (4 tracers)
        IF(id_fine>0) flux_tr(i,id_fine)=
     . scale_param_dustacc(iregion_dust(i))*
     .                  dust_ec(i)*1.e3*0.114*1.e3  !mg/m2/s
        IF(id_codu>0) flux_tr(i,id_codu)=
     . scale_param_dustcoa(iregion_dust(i))*
     .                  dust_ec(i)*1.e3*0.108*1.e3  !mg/m2/s bin 0.5-3um
        IF(id_scdu>0) flux_tr(i,id_scdu)=
     . scale_param_dustsco(iregion_dust(i))*
     .                  dust_ec(i)*1.e3*0.778*1.e3  !mg/m2/s bin 0.5-3um

         flux_sparam_ddfine(i)=scale_param_dustacc(iregion_dust(i)) *
     .                            dust_ec(i)*1.e3*0.114*1.e3
         flux_sparam_ddcoa(i)=scale_param_dustcoa(iregion_dust(i)) *
     .                            dust_ec(i)*1.e3*0.108*1.e3
         flux_sparam_ddsco(i)=scale_param_dustsco(iregion_dust(i)) *
     .                            dust_ec(i)*1.e3*0.778*1.e3

         endif
      ENDDO





      ENDIF
!*****************************************************                                   
C******************* SEA SALT EMMISSIONS *******************************
      DO i=1,klon
         pct_ocean(i)=pctsrf(i,is_oce)
      ENDDO
c
!      IF (lminmax) THEN
      DO j=1,nbtr
      DO i=1,klon
         tmp_var2(i,j)=source_tr(i,j)
      ENDDO
      ENDDO
      CALL minmaxsource(tmp_var2,qmin,qmax,'src: before SS emiss')
      IF(id_coss>0) then
      print *,'Source = ',SUM(source_tr(:,id_coss)),
     .     MINVAL(source_tr(:,id_coss)), MAXVAL(source_tr(:,id_coss))
      ENDIF
!
      DO i=1,klon
! Original line (4 tracers)
         IF(id_fine>0) source_tr(i,id_fine)=
     . source_tr(i,id_fine)+scale_param_ssacc*
     .                                 lmt_sea_salt(i,1)*1.e4       !g/m2/s

! Original line (4 tracers)
       IF(id_fine>0)  flux_tr(i,id_fine)=
     . flux_tr(i,id_fine)+scale_param_ssacc
     .                            *lmt_sea_salt(i,1)*1.e4*1.e3      !mg/m2/s
!
      IF(id_coss>0)  source_tr(i,id_coss)=
     . scale_param_sscoa*lmt_sea_salt(i,2)*1.e4    !g/m2/s
      IF(id_coss>0)  flux_tr(i,id_coss)=
     . scale_param_sscoa*lmt_sea_salt(i,2)*1.e4*1.e3 !mg/m2/s
c
         flux_sparam_ssfine(i)=scale_param_ssacc *
     .                                  lmt_sea_salt(i,1)*1.e4*1.e3
         flux_sparam_sscoa(i)=scale_param_sscoa *
     .                                  lmt_sea_salt(i,2)*1.e4*1.e3
      ENDDO
!      IF (lminmax) THEN
      DO j=1,nbtr
      DO i=1,klon
         tmp_var2(i,j)=source_tr(i,j)
      ENDDO
      ENDDO
      CALL minmaxsource(tmp_var2,qmin,qmax,'src: after SS emiss')
      IF(id_coss>0) then
      print *,'Source = ',SUM(source_tr(:,id_coss)),
     .  MINVAL(source_tr(:,id_coss)), MAXVAL(source_tr(:,id_coss))
      ENDIF
c     

      END
