c----- This subroutine calculates the sedimentation flux of Tracers 
c
      SUBROUTINE sediment_mod(t_seri,pplay,zrho,paprs,time_step,RHcl,
     .                                       id_coss,id_codu,id_scdu,
     .                                        ok_chimeredust,
     .                           sed_ss,sed_dust,sed_dustsco,
     .                        sed_ss3D,sed_dust3D,sed_dustsco3D,tr_seri)
cnhl     .                                       xlon,xlat,
c
       USE dimphy
       USE infotrac
      IMPLICIT NONE 
c
#include "dimensions.h"
#include "chem.h"
c #include "dimphy.h"
#include "YOMCST.h"
#include "YOECUMF.h"
c
       REAL RHcl(klon,klev)     ! humidite relative ciel clair
       REAL tr_seri(klon, klev,nbtr) !conc of tracers
       REAL sed_ss(klon) !sedimentation flux of Sea Salt (g/m2/s)
       REAL sed_dust(klon) !sedimentation flux of dust (g/m2/s)
       REAL sed_dustsco(klon) !sedimentation flux of scoarse  dust (g/m2/s)
       REAL sed_ss3D(klon,klev) !sedimentation flux of Sea Salt (g/m2/s)
       REAL sed_dust3D(klon,klev) !sedimentation flux of dust (g/m2/s)
       REAL sed_dustsco3D(klon,klev) !sedimentation flux of scoarse  dust (g/m2/s)
       REAL t_seri(klon, klev)   !Temperature at mid points of Z (K)
       REAL v_dep_ss(klon,klev)  ! sed. velocity for SS m/s
       REAL v_dep_dust(klon,klev)  ! sed. velocity for dust m/s
       REAL v_dep_dustsco(klon,klev)  ! sed. velocity for dust m/s
       REAL pplay(klon, klev)    !pressure at mid points of Z (Pa)
       REAL zrho(klon, klev)     !Density of air at mid points of Z (kg/m3)
       REAL paprs(klon, klev+1)    !pressure at interface of layers Z (Pa)
       REAL time_step            !time step (sec)
       LOGICAL ok_chimeredust
       REAL xlat(klon)       ! latitudes pour chaque point 
       REAL xlon(klon)       ! longitudes pour chaque point 
       INTEGER id_coss,id_codu,id_scdu
c
c------local variables
c 
       INTEGER i, k, nbre_RH
       PARAMETER(nbre_RH=12)
c
       REAL lambda, ss_g            
       REAL mmd_ss      !mass median diameter of SS (um) 
       REAL mmd_dust          !mass median diameter of dust (um) 
       REAL mmd_dustsco          !mass median diameter of scoarse dust (um) 
       REAL rho_ss(nbre_RH),rho_ss1 !density of sea salt (kg/m3)
       REAL rho_dust          !density of dust(kg/m3)
       REAL v_stokes, CC, v_sed, ss_growth_f(nbre_RH)
       REAL sed_flux(klon,klev)  ! sedimentation flux g/m2/s
       REAL air_visco(klon,klev)
       REAL zdz(klon,klev)       ! layers height (m)
       REAL temp                 ! temperature in degree Celius 
c
       INTEGER RH_num
       REAL RH_MAX, DELTA, rh, RH_tab(nbre_RH)
       PARAMETER (RH_MAX=95.)
c
       DATA RH_tab/0.,10.,20.,30.,40.,50.,60.,70.,80.,85.,90.,95./
c
c
       DATA rho_ss/2160. ,2160. ,2160.,  2160,  1451.6, 1367.9,
     .             1302.9,1243.2,1182.7, 1149.5,1111.6, 1063.1/
c
       DATA ss_growth_f/0.503, 0.503, 0.503, 0.503, 0.724, 0.782, 
     .                  0.838, 0.905, 1.000, 1.072, 1.188, 1.447/
c
c
       mmd_ss=12.7   !dia -um at 80% for bin 0.5-20 um but 90% of real mmd
!               obsolete      mmd_dust=2.8  !micrometer for bin 0.5-20 and 0.5-10 um
! 4tracer SPLA:       mmd_dust=11.0  !micrometer for bin 0.5-20 and 0.5-10 um
!3days       mmd_dust=3.333464  !micrometer for bin 0.5-20 and 0.5-10 um
!3days       mmd_dustsco=12.91315  !micrometer for bin 0.5-20 and 0.5-10 um
!JE20140911       mmd_dust=3.002283  !micrometer for bin 0.5-20 and 0.5-10 um
!JE20140911       mmd_dustsco=13.09771  !micrometer for bin 0.5-20 and 0.5-10 um
!JE20140911        mmd_dust=5.156346  !micrometer for bin 0.5-20 and 0.5-10 um
!JE20140911        mmd_dustsco=15.56554  !micrometer for bin 0.5-20 and 0.5-10 um
        IF (ok_chimeredust) THEN 
!JE20150212<< : changes in ustar in dustmod changes emission distribution
!        mmd_dust=3.761212  !micrometer for bin 0.5-3 and 0.5-10 um
!        mmd_dustsco=15.06167  !micrometer for bin 3-20 and 0.5-10 um
!JE20150212>>
!JE20150618: Change in div3 of dustmod changes distribution. now is div3=6
!div=3        mmd_dust=3.983763
!div=3        mmd_dustsco=15.10854 
        mmd_dust=3.898047
        mmd_dustsco=15.06167 
        ELSE 
        mmd_dust=11.0  !micrometer for bin 0.5-20 and 0.5-10 um
        mmd_dustsco=100. ! absurd value, bin not used in this scheme
        ENDIF


       rho_dust=2600. !kg/m3
c
c--------- Air viscosity (poise=0.1 kg/m-sec)-----------
c
       DO k=1, klev
       DO i=1, klon
c
       zdz(i,k)=(paprs(i,k)-paprs(i,k+1))/zrho(i,k)/RG
c
       temp=t_seri(i,k)-RTT
c
       IF (temp.LT.0.) THEN
         air_visco(i,k)=(1.718+0.0049*temp-1.2e-5*temp*temp)*1.e-4 
       ELSE
         air_visco(i,k)=(1.718+0.0049*temp)*1.e-4 
       ENDIF
c
       ENDDO
       ENDDO
c
c--------- for Sea Salt -------------------
c
c
c
       IF(id_coss>0) THEN
       DO k=1, klev
       DO i=1,klon
c
c---cal. correction factor hygroscopic growth of aerosols
c
        rh=MIN(RHcl(i,k)*100.,RH_MAX)
        RH_num = INT( rh/10. + 1.)
        IF (rh.gt.85.) RH_num=10
        IF (rh.gt.90.) RH_num=11
        DELTA=(rh-RH_tab(RH_num))/(RH_tab(RH_num+1)-RH_tab(RH_num))
c
        ss_g=ss_growth_f(rh_num) +
     .       DELTA*(ss_growth_f(RH_num+1)-ss_growth_f(RH_num))

        rho_ss1=rho_ss(rh_num) +
     .       DELTA*(rho_ss(RH_num+1)-rho_ss(RH_num))              
c
        v_stokes=RG*(rho_ss1-zrho(i,k))*      !m/sec
     .           (mmd_ss*ss_g)*(mmd_ss*ss_g)*
     .           1.e-12/(18.0*air_visco(i,k)/10.)
c
       lambda=6.6*1.e-8*(103125/pplay(i,k))*(t_seri(i,k)/293.15)
c
       CC=1.0+1.257*lambda/(mmd_ss*ss_g)/1.e6  ! C-correction factor 
c
       v_sed=v_stokes*CC                       ! m/sec !orig
c
c---------check for v_sed*dt<zdz
c
       IF (v_sed*time_step.GT.zdz(i,k)) THEN
         v_sed=zdz(i,k)/time_step      
       ENDIF
c
       v_dep_ss(i,k)= v_sed
       sed_flux(i,k)= tr_seri(i,k,id_coss)*v_sed !g/cm3*m/sec
       !sed_ss3D(i,k)= -sed_flux(i,k)/zdz(i,k)      !g/cm3*sec !!!!!!!
      ! conc_sed_ss3D(i,k)=sed_flux(i,k)*1.e6      !g/m3*sec !!!!!!!
c
       ENDDO          !klon
       ENDDO          !klev
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       sed_ss3D(:,:)=0.0  ! initialisation 
      
       DO k=1, klev
       DO i=1, klon
       sed_ss3D(i,k)=sed_ss3D(i,k)-
     .        sed_flux(i,k)/zdz(i,k) !!!!!!!!!!!!!!!!!!!!!!
       ENDDO          !klon
       ENDDO          !klev
c
       DO k=1, klev-1
       DO i=1, klon
        sed_ss3D(i,k)=sed_ss3D(i,k)+                   
     .                  sed_flux(i,k+1)/zdz(i,k) !!!!!!!!

       ENDDO          !klon
       ENDDO          !klev

      DO k = 1, klev
      DO i = 1, klon
          tr_seri(i,k,id_coss)=tr_seri(i,k,id_coss)+
     s   sed_ss3D(i,k)*time_step
      ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
       DO i=1, klon
         sed_ss(i)=sed_flux(i,1)*1.e6*1.e3    !--unit mg/m2/s
       ENDDO          !klon
       ELSE
        DO i=1, klon
          sed_ss(i)=0.
        ENDDO
       ENDIF
c
c

c--------- For dust ------------------
c
c
       IF(id_codu>0) THEN
       DO k=1, klev
       DO i=1,klon
c
        v_stokes=RG*(rho_dust-zrho(i,k))*      !m/sec
     .           mmd_dust*mmd_dust*
     .           1.e-12/(18.0*air_visco(i,k)/10.)
c
       lambda=6.6*1.e-8*(103125/pplay(i,k))*(t_seri(i,k)/293.15)
       CC=1.0+1.257*lambda/(mmd_dust)/1.e6        !dimensionless
       v_sed=v_stokes*CC                       !m/sec
c
c---------check for v_sed*dt<zdz
c
       IF (v_sed*time_step.GT.zdz(i,k)) THEN
         v_sed=zdz(i,k)/time_step      
       ENDIF

c
       v_dep_dust(i,k)= v_sed
       sed_flux(i,k)  = tr_seri(i,k,id_codu)*v_sed !g/cm3.m/sec
       !sed_dust3D(i,k)= -sed_flux(i,k)/zdz(i,k)      !g/cm3*sec !!!!!!!
c
       ENDDO          !klon
       ENDDO          !klev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       sed_dust3D(:,:)=0.0  ! initialisation 

       DO k=1, klev
       DO i=1, klon
       sed_dust3D(i,k)=sed_dust3D(i,k)-
     .                  sed_flux(i,k)/zdz(i,k)
       ENDDO          !klon
       ENDDO          !klev

c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       DO k=1, klev-1
       DO i=1, klon
        sed_dust3D(i,k)=sed_dust3D(i,k) +
     .                  sed_flux(i,k+1)/zdz(i,k)
       ENDDO          !klon
       ENDDO          !klev
c
      DO k = 1, klev
      DO i = 1, klon
         tr_seri(i,k,id_codu)=tr_seri(i,k,id_codu)+
     s    sed_dust3D(i,k)*time_step
      ENDDO
      ENDDO


       DO i=1, klon
         sed_dust(i)=sed_flux(i,1)*1.e6*1.e3    !--unit mg/m2/s
       ENDDO          !klon
       ELSE
        DO i=1, klon
          sed_dust(i)=0.
        ENDDO
       ENDIF
c


c--------- For scoarse  dust ------------------
c
c
       IF(id_scdu>0) THEN
       DO k=1, klev
       DO i=1,klon
c
        v_stokes=RG*(rho_dust-zrho(i,k))*      !m/sec
     .           mmd_dustsco*mmd_dustsco*
     .           1.e-12/(18.0*air_visco(i,k)/10.)
c
       lambda=6.6*1.e-8*(103125/pplay(i,k))*(t_seri(i,k)/293.15)
       CC=1.0+1.257*lambda/(mmd_dustsco)/1.e6        !dimensionless
       v_sed=v_stokes*CC                       !m/sec
c
c---------check for v_sed*dt<zdz


       IF (v_sed*time_step.GT.zdz(i,k)) THEN
         v_sed=zdz(i,k)/time_step
       ENDIF

c
       v_dep_dustsco(i,k)= v_sed
       sed_flux(i,k)     = tr_seri(i,k,id_scdu)*v_sed !g/cm3.m/sec
       !sed_dustsco3D(i,k)= -sed_flux(i,k)/zdz(i,k)      !g/cm3*sec !!!!!!!
c
       ENDDO          !klon
       ENDDO          !klev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       sed_dustsco3D(:,:)=0.0  ! initialisation 

       DO k=1, klev
       DO i=1, klon
       sed_dustsco3D(i,k)=sed_dustsco3D(i,k)-
     .                  sed_flux(i,k)/zdz(i,k)
       ENDDO          !klon
       ENDDO          !klev
c
       DO k=1, klev-1
       DO i=1, klon
        sed_dustsco3D(i,k)=sed_dustsco3D(i,k) +
     .                  sed_flux(i,k+1)/zdz(i,k)
       ENDDO          !klon
       ENDDO          !klev

      DO k = 1, klev
      DO i = 1, klon
       tr_seri(i,k,id_scdu)=tr_seri(i,k,id_scdu)+
     s  sed_dustsco3D(i,k)*time_step
      ENDDO
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c
       DO i=1, klon
         sed_dustsco(i)=sed_flux(i,1)*1.e6*1.e3    !--unit mg/m2/s
       ENDDO          !klon
       ELSE
        DO i=1, klon
          sed_dustsco(i)=0.
        ENDDO
       ENDIF
c




c
       RETURN
       END
