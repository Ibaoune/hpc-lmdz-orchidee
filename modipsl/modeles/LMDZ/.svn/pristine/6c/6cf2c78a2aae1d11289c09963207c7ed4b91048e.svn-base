C Routine to read the emissions of the different species
C
      subroutine read_newemissions(julien, jH_emi ,edgar, flag_dms,
     I                             debutphy, 
     I                             pdtphys,lafinphy, nbjour, pctsrf,
     I                             t_seri, xlat, xlon,
     I                             pmflxr, pmflxs, prfl, psfl,
     O                             u10m_ec, v10m_ec, dust_ec, 
     O                             lmt_sea_salt, lmt_so2ff_l, 
     O                             lmt_so2ff_h, lmt_so2nff, lmt_so2ba,
     O                             lmt_so2bb_l, lmt_so2bb_h,
     O                             lmt_so2volc_cont, lmt_altvolc_cont,
     O                             lmt_so2volc_expl, lmt_altvolc_expl,
     O                             lmt_dmsbio, lmt_h2sbio, lmt_dmsconc, 
     O                             lmt_bcff, lmt_bcnff, lmt_bcbb_l, 
     O                             lmt_bcbb_h, lmt_bcba, lmt_omff, 
     O                             lmt_omnff, lmt_ombb_l, lmt_ombb_h, 
     O                             lmt_omnat, lmt_omba)
      
      USE dimphy
      USE indice_sol_mod
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para

      IMPLICIT NONE


#include "dimensions.h"
c      INCLUDE 'dimphy.h'      
      INCLUDE 'paramet.h'      
      INCLUDE 'chem.h'      
      INCLUDE 'chem_spla.h'
c      INCLUDE 'indicesol.h'

      logical debutphy, lafinphy, edgar
      INTEGER test_vent, test_day, step_vent, flag_dms, nbjour
      INTEGER julien, i, iday
      SAVE step_vent, test_vent, test_day, iday
!$OMP THREADPRIVATE(step_vent, test_vent, test_day, iday)
      REAL pct_ocean(klon), pctsrf(klon,nbsrf)
      REAL pdtphys  ! pas d'integration pour la physique (seconde)      
      REAL t_seri(klon,klev)  ! temperature

      REAL xlat(klon)       ! latitudes pour chaque point 
      REAL xlon(klon)       ! longitudes pour chaque point 
      
c
c   Emissions:
c   ---------
c
c---------------------------- SEA SALT & DUST emissions ------------------------
      REAL lmt_sea_salt(klon,ss_bins) !Sea salt 0.03-8.0 um !NOT SAVED OK
      REAL clyfac, avgdryrate, drying
c je      REAL u10m_ec1(klon), v10m_ec1(klon), dust_ec1(klon)
c je      REAL u10m_ec2(klon), v10m_ec2(klon), dust_ec2(klon)

      REAL, SAVE, ALLOCATABLE :: u10m_ec1(:), v10m_ec1(:), dust_ec1(:)
      REAL, SAVE, ALLOCATABLE :: u10m_ec2(:), v10m_ec2(:), dust_ec2(:)
!$OMP THREADPRIVATE(u10m_ec1, v10m_ec1, dust_ec1)
!$OMP THREADPRIVATE(u10m_ec2, v10m_ec2, dust_ec2)
      REAL u10m_nc(iip1,jjp1), v10m_nc(iip1,jjp1)
      REAL u10m_ec(klon), v10m_ec(klon), dust_ec(klon)
c      REAL cly(klon), wth(klon), zprecipinsoil(klon)
      REAL, SAVE, ALLOCATABLE :: cly(:), wth(:), zprecipinsoil(:)
      REAL :: cly_glo(klon_glo), wth_glo(klon_glo)
      REAL :: zprecipinsoil_glo(klon_glo)
!$OMP THREADPRIVATE(cly,wth,zprecipinsoil)


c je     SAVE u10m_ec2, v10m_ec2, dust_ec2
c je      SAVE u10m_ec1, v10m_ec1, dust_ec1   ! Added on titane
c je      SAVE cly, wth, zprecipinsoil        ! Added on titane
!     SAVE cly, wth, zprecipinsoil, u10m_ec2, v10m_ec2, dust_ec2
c------------------------- BLACK CARBON emissions ----------------------
      REAL lmt_bcff(klon)       ! emissions de BC fossil fuels
      REAL lmt_bcnff(klon)      ! emissions de BC non-fossil fuels
      REAL lmt_bcbb_l(klon)     ! emissions de BC biomass basses
      REAL lmt_bcbb_h(klon)     ! emissions de BC biomass hautes
      REAL lmt_bcba(klon)       ! emissions de BC bateau
c------------------------ ORGANIC MATTER emissions ---------------------
      REAL lmt_omff(klon)     ! emissions de OM fossil fuels
      REAL lmt_omnff(klon)    ! emissions de OM non-fossil fuels
      REAL lmt_ombb_l(klon)   ! emissions de OM biomass basses
      REAL lmt_ombb_h(klon)   ! emissions de OM biomass hautes
      REAL lmt_omnat(klon)    ! emissions de OM Natural
      REAL lmt_omba(klon)     ! emissions de OM bateau
c------------------------- SULFUR emissions ----------------------------
      REAL lmt_so2ff_l(klon)       ! emissions so2 fossil fuels (low)
      REAL lmt_so2ff_h(klon)       ! emissions so2 fossil fuels (high)
      REAL lmt_so2nff(klon)        ! emissions so2 non-fossil fuels
      REAL lmt_so2bb_l(klon)       ! emissions de so2 biomass burning basse
      REAL lmt_so2bb_h(klon)       ! emissions de so2 biomass burning hautes
      REAL lmt_so2ba(klon)         ! emissions de so2 bateau
      REAL lmt_so2volc_cont(klon)  ! emissions so2 volcan continuous
      REAL lmt_altvolc_cont(klon)  ! altitude  so2 volcan continuous
      REAL lmt_so2volc_expl(klon)  ! emissions so2 volcan explosive
      REAL lmt_altvolc_expl(klon)  ! altitude  so2 volcan explosive
      REAL lmt_dmsconc(klon)       ! concentration de dms oceanique
      REAL lmt_dmsbio(klon)        ! emissions de dms bio
      REAL lmt_h2sbio(klon)        ! emissions de h2s bio

      REAL,SAVE,ALLOCATABLE ::  lmt_dms(:)           ! emissions de dms
!$OMP THREADPRIVATE(lmt_dms)
c
c  Lessivage
c  ---------
c
      REAL pmflxr(klon,klev+1), pmflxs(klon,klev+1) !--convection
      REAL prfl(klon,klev+1),   psfl(klon,klev+1)   !--large-scale
!      REAL pmflxr(klon,klev), pmflxs(klon,klev) !--convection
!      REAL prfl(klon,klev),   psfl(klon,klev)   !--large-scale
c
c  Variable interne
c  ----------------
c      
      INTEGER icount
      REAL tau_1, tau_2                  
      REAL max_flux, min_flux                                          
      INTRINSIC MIN, MAX
c
c JE: Changes due to new pdtphys in new physics.
c      REAL windintime ! time in hours of the wind input files resolution 
c      REAL dayemintime ! time in hours of the other emissions input files resolution 
      REAL jH_init ! shift in the hour (count as days) respecto to
!                  ! realhour = (pdtphys*i)/3600/24 -days_elapsed
      REAL jH_emi,jH_vent,jH_day
      SAVE jH_init,jH_vent,jH_day
!$OMP THREADPRIVATE(jH_init,jH_vent,jH_day)
      REAL,PARAMETER :: vent_resol = 6. ! resolution of winds in hours
      REAL,PARAMETER :: day_resol = 24. ! resolution of daily emmis. in hours
!      INTEGER   test_day1
!      SAVE test_day1
!      REAL tau_1j,tau_2j
c je
c allocate if necessary
c

      IF (.NOT. ALLOCATED(u10m_ec1)) ALLOCATE(u10m_ec1(klon))
      IF (.NOT. ALLOCATED(v10m_ec1)) ALLOCATE(v10m_ec1(klon))
      IF (.NOT. ALLOCATED(dust_ec1)) ALLOCATE(dust_ec1(klon))
      IF (.NOT. ALLOCATED(u10m_ec2)) ALLOCATE(u10m_ec2(klon))
      IF (.NOT. ALLOCATED(v10m_ec2)) ALLOCATE(v10m_ec2(klon))
      IF (.NOT. ALLOCATED(dust_ec2)) ALLOCATE(dust_ec2(klon))
      IF (.NOT. ALLOCATED(cly)) ALLOCATE(cly(klon))
      IF (.NOT. ALLOCATED(wth)) ALLOCATE(wth(klon))
      IF (.NOT. ALLOCATED(zprecipinsoil)) ALLOCATE(zprecipinsoil(klon))
      IF (.NOT. ALLOCATED(lmt_dms)) ALLOCATE(lmt_dms(klon))
c end je nov2013
c
C***********************************************************************
C DUST EMISSIONS
C***********************************************************************
c
      IF (debutphy) THEN 
C---Fields are read only at the beginning of the period
c--reading wind and dust
        iday=julien
        step_vent=1
        test_vent=0
        test_day=0
        CALL read_vent(.true.,step_vent,nbjour,u10m_ec2,v10m_ec2)
        print *,'Read (debut) dust emissions: step_vent,julien,nbjour', 
     .                                   step_vent,julien,nbjour
        CALL read_dust(.true.,step_vent,nbjour,dust_ec2)
C Threshold velocity map
!$OMP MASTER
       IF (is_mpi_root .AND. is_omp_root) THEN
        zprecipinsoil_glo(:)=0.0
        OPEN(51,file='wth.dat',status='unknown',form='formatted')
        READ(51,'(G18.10)') (wth_glo(i),i=1,klon_glo)
        CLOSE(51)
c Clay content
        OPEN(52,file='cly.dat',status='unknown',form='formatted')
        READ(52,'(G18.10)') (cly_glo(i),i=1,klon_glo)
        CLOSE(52)
        OPEN(53,file='precipinsoil.dat',
     .        status='old',form='formatted',err=999)
        READ(53,'(G18.10)') (zprecipinsoil_glo(i),i=1,klon_glo)
        PRINT *,'lecture precipinsoil.dat'
 999    CONTINUE
        CLOSE(53)
       ENDIF
!$OMP END MASTER
!$OMP BARRIER
       call scatter(wth_glo,wth)
       call scatter(cly_glo,cly)
       call scatter(zprecipinsoil_glo,zprecipinsoil)

!JE20140908<<        GOTO 1000
!        DO i=1, klon
!          zprecipinsoil(i)=0.0
!        ENDDO
! 1000   CLOSE(53)
!JE20140908>>
        jH_init=jH_emi
        jH_vent=jH_emi
        jH_day=jH_emi
!        test_day1=0
!JE end
c
      
      ENDIF !--- debutphy
        
      print *,'READ_EMISSION: test_vent & test_day = ',test_vent, 
     +                                                 test_day
      IF (test_vent.EQ.0) THEN    !--on lit toutes les 6 h                                                                      
        CALL SCOPY(klon, u10m_ec2, 1, u10m_ec1, 1)
        CALL SCOPY(klon, v10m_ec2, 1, v10m_ec1, 1)
        CALL SCOPY(klon, dust_ec2, 1, dust_ec1, 1)
        step_vent=step_vent+1
        !PRINT *,'step_vent=', step_vent
        CALL read_vent(.false.,step_vent,nbjour,u10m_ec2,v10m_ec2)
        print *,'Reading dust emissions: step_vent, julien, nbjour ', 
     .                                   step_vent, julien, nbjour
        !print *,'test_vent, julien = ',test_vent, julien
        CALL read_dust(.false.,step_vent,nbjour,dust_ec2)
      
      ENDIF !--test_vent

c     ubicacion original
c      test_vent=test_vent+1
c      IF (test_vent.EQ.(6*2)) test_vent=0 !on remet a zero ttes les 6 h
            
!JE      tau_2=FLOAT(test_vent)/12.
!JE      tau_1=1.-tau_2
      tau_2=(jH_vent-jH_init)*24./(vent_resol)
      tau_1=1.-tau_2
!      print*,'JEdec jHv,JHi,ventres',jH_vent,jH_init,vent_resol
!      print*,'JEdec tau2,tau1',tau_2,tau_1
!      print*,'JEdec step_vent',step_vent
      DO i=1, klon
!      PRINT*,'JE tau_2,tau_2j',tau_2,tau_2j
        u10m_ec(i)=tau_1*u10m_ec1(i)+tau_2*u10m_ec2(i)
        v10m_ec(i)=tau_1*v10m_ec1(i)+tau_2*v10m_ec2(i)
        dust_ec(i)=tau_1*dust_ec1(i)+tau_2*dust_ec2(i)
      ENDDO
c
cJE      IF (test_vent.EQ.(6*2)) THEN
cJE        PRINT *,'6 hrs interval reached'
cJE        print *,'day in read_emission, test_vent = ',julien, test_vent
cJE      ENDIF
cJE
!JE      test_vent=test_vent+1
!JE      IF (test_vent.EQ.(6*2)) test_vent=0 !on remet a zero ttes les 6 h
c JE  
      jH_vent=jH_vent+pdtphys/(24.*3600.)
      test_vent=test_vent+1
      IF (jH_vent.GT.(vent_resol)/24.) THEN
          test_vent=0
          jH_vent=jH_init
      ENDIF
!      PRINT*,'JE test_vent,test_vent1,jH_vent ', test_vent,test_vent1
!     .     ,jH_vent
c endJEi
c
      avgdryrate=300./365.*pdtphys/86400.
c
      DO i=1, klon
c
        IF (cly(i).LT.9990..AND.wth(i).LT.9990.) THEN
          zprecipinsoil(i)=zprecipinsoil(i) +
     .           (pmflxr(i,1)+pmflxs(i,1)+prfl(i,1)+psfl(i,1))*pdtphys
c
          clyfac=MIN(16., cly(i)*0.4+8.) ![mm] max amount of water hold in top soil
          drying=avgdryrate*exp(0.03905491*
     .                    exp(0.17446*(t_seri(i,1)-273.15))) ! [mm]
          zprecipinsoil(i)=min(max(0.,zprecipinsoil(i)-drying),clyfac) ! [mm]
        ENDIF
!        zprecipinsoil(i)=0.0 ! Temporarely introduced to reproduce obelix result
      ENDDO

!      print *,'cly = ',sum(cly),maxval(cly),minval(cly)
!      print *,'wth = ',sum(wth),maxval(wth),minval(wth)
!      print *,'t_seri = ',sum(t_seri),maxval(t_seri),minval(t_seri)
!      print *,'precipinsoil = ',sum(zprecipinsoil),maxval(zprecipinsoil)
!     .                      ,minval(zprecipinsoil)
      icount=0
      DO i=1, klon
        IF (cly(i).GE.9990..OR.wth(i).GE.9990..OR.
     .     t_seri(i,1).LE.273.15.OR.zprecipinsoil(i).GT.1.e-8) THEN
             dust_ec(i)=0.0 ! commented out for test dustemtest
!             print *,'Dust emissions surpressed at grid = ',i
!             icount=icount+1
        ENDIF
      ENDDO                            
c
      print *,'Total N of grids with surpressed emission = ',icount
      print *,'dust_ec = ',SUM(dust_ec),MINVAL(dust_ec),
     .                                  MAXVAL(dust_ec)
cnhl Transitory scaling of desert dust emissions
      
cnhl      DO i=1, klon
cnhl         dust_ec(i)=dust_ec(i)/2.
cnhl      ENDDO                            

C-saving precipitation field to be read in next simulation
      
      IF (lafinphy) THEN
c
        CALL gather(zprecipinsoil,zprecipinsoil_glo)
!$OMP MASTER
        IF (is_mpi_root .AND. is_omp_root) THEN

        OPEN(53,file='newprecipinsoil.dat',
     .          status='unknown',form='formatted')
        WRITE(53,'(G18.10)') (zprecipinsoil_glo(i),i=1,klon_glo)
        CLOSE(53)
        ENDIF
!$OMP END MASTER
!$OMP BARRIER
c
      ENDIF
c
C***********************************************************************
C SEA SALT EMISSIONS
C***********************************************************************
c
      DO i=1,klon
        pct_ocean(i)=pctsrf(i,is_oce)
      ENDDO

      print *,'IS_OCE = ',is_oce
      CALL seasalt(v10m_ec, u10m_ec, pct_ocean, lmt_sea_salt) !mgSeaSalt/cm2/s
!      print *,'SUM, MAX & MIN Sea Salt = ',SUM(lmt_sea_salt),
!     .               MAXVAL(lmt_sea_salt),MINVAL(lmt_sea_salt)
c      
C***********************************************************************
C SULFUR & CARBON EMISSIONS
C***********************************************************************
c
      
      IF (test_day.EQ.0) THEN
        print *,'Computing SULFATE emissions for day : ',iday,julien, 
     .                                                   step_vent
        CALL condsurfs_new(iday, edgar, flag_dms,
     O                      lmt_so2ff_l, lmt_so2ff_h, lmt_so2nff, 
     O                      lmt_so2bb_l, lmt_so2bb_h, lmt_so2ba, 
     O                      lmt_so2volc_cont, lmt_altvolc_cont, 
     O                      lmt_so2volc_expl, lmt_altvolc_expl, 
     O                      lmt_dmsbio, lmt_h2sbio, lmt_dms,lmt_dmsconc)
        print *,'Computing CARBON emissions for day : ',iday,julien,
     .                                                   step_vent
        CALL condsurfc_new(iday,
     O                       lmt_bcff,lmt_bcnff,lmt_bcbb_l,lmt_bcbb_h,
     O                       lmt_bcba,lmt_omff,lmt_omnff,lmt_ombb_l,
     O                       lmt_ombb_h, lmt_omnat, lmt_omba)
        print *,'IDAY = ',iday
        iday=iday+1
        print *,'BCBB_L emissions :',SUM(lmt_bcbb_l), MAXVAL(lmt_bcbb_l)
     .                              ,MINVAL(lmt_bcbb_l)
        print *,'BCBB_H emissions :',SUM(lmt_bcbb_h), MAXVAL(lmt_bcbb_h)
     .                              ,MINVAL(lmt_bcbb_h)
      ENDIF
      
!JE      test_day=test_day+1
!JE      IF (test_day.EQ.(24*2.)) THEN
!JE        test_day=0 !on remet a zero ttes les 24 h
!JE        print *,'LAST TIME STEP OF DAY ',julien
!JE      ENDIF


      jH_day=jH_day+pdtphys/(24.*3600.)
      test_day=test_day+1
      IF (jH_day.GT.(day_resol)/24.) THEN
          print *,'LAST TIME STEP OF DAY ',julien
          test_day=0
          jH_day=jH_init
      ENDIF
!      PRINT*,'test_day,test_day1',test_day,test_day1

      END
