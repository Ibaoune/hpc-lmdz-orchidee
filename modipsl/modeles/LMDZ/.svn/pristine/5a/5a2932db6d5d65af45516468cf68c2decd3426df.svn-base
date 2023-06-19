
SUBROUTINE satellite_out_spla(jD_cur,jH_cur,pdtphys,rlat,rlon, &
                              masque_aqua, masque_terra )

  USE dimphy
  USE IOIPSL
  USE write_field_phy

  IMPLICIT NONE

  REAL :: pdtphys, jD_cur,jH_cur, hour
  INTEGER :: year_cur, mth_cur, day_cur
  INTEGER :: masque_polder(klon)  ! masque polder
  INTEGER :: masque_aqua(klon)  ! masque polder
  INTEGER :: masque_terra(klon)  ! masque polder
  INTEGER :: i
  REAL :: overpassaqua, overpassterra
  REAL,dimension(klon) :: rlat,rlon


  masque_polder(:) = 0.

  CALL ju2ymds(jD_cur+jH_cur, year_cur, mth_cur, day_cur, hour)
!  print*,'JDcur=',jD_cur,'JHcur=',jH_cur,'year_cur' ,year_cur,'mth_cur' ,mth_cur, 'day_cur',day_cur,'hour' ,hour

!  IF ( (year_cur*100.+mth_cur .GE. 199611 ) .AND. (year_cur*100.+mth_cur .LE. 199706)) THEN 
!          CALL swathpolder(year_cur,mth_cur,day_cur,hour/86400., &
!            pdtphys,rlon,rlat,masque_polder)
!  ENDIF
!
!  DO i=1,klon
!    IF ( masque_polder(i) .EQ. 1 ) THEN
!        print *,'polder output point, lon:', rlon(i),', lat: ',rlat(i)
!    ENDIF
!  ENDDO
!  CALL writefield_phy("masque_polder",float(masque_polder),1)

! Aqua
  masque_aqua(:) = 0.
  overpassaqua=48600.  ! 13.30 p.m. local time
  CALL swathpolarsat(year_cur,mth_cur,day_cur,hour, &
            pdtphys,rlon,rlat,overpassaqua,masque_aqua)

!  DO i=1,klon
!    IF ( masque_aqua(i) .EQ. 1 ) THEN
!        print *,'aqua output point, lon:', rlon(i),', lat: ',rlat(i)
!    ENDIF
!  ENDDO
!  CALL writefield_phy("masque_aqua",float(masque_aqua),1)

  masque_terra(:) = 0.
  overpassterra=37800.  ! 10.30 a.m. local time
  CALL swathpolarsat(year_cur,mth_cur,day_cur,hour, &
            pdtphys,rlon,rlat,overpassterra,masque_terra)

!  DO i=1,klon
!    IF ( masque_terra(i) .EQ. 1 ) THEN
!        print *,'terra output point, lon:', rlon(i),', lat: ',rlat(i)
!    ENDIF
!  ENDDO
!  CALL writefield_phy("masque_terra",float(masque_terra),1)

  RETURN
END SUBROUTINE satellite_out_spla




SUBROUTINE swathpolarsat(annee,mois,jour,heure,pdtphys, &
     rlon,rlat,overpasstime,masque)
  ! Adaptation from the simple satellite simulator of AeroCom working group Indirect
  ! forcing ( Johannes Quaas, MPI for Meteorology, Hamburg )
  ! http://wiki.esipfed.org/index.php/Indirect_forcing
  
  USE dimphy
  IMPLICIT NONE

  INTEGER :: annee, mois, jour, i
  REAL :: heure                      !--heure en jour
  REAL :: pdtphys                    !--pas de temps en seconde
  REAL :: rlon(1:klon), rlat(1:klon) !--longitude et latitude
  INTEGER :: masque(1:klon)
  REAL :: localtime(1:klon)
  REAL :: overpasstime, utctime

  masque(:) = 0
  utctime=heure
  DO i=1,klon
        localtime(i) = utctime + 240. * rlon(i) ! for each degree of longitude east,4 min earlier local time
!        IF ( localtime(i) > 86400. ) THEN ! this is still the previous day
!                localtime(i) = localtime(i) - 86400.
!        ENDIF
        ! Select 10.30 a.m. ± dt/2
!       IF ( ABS( localtime(i) - overpasstime ) <= pdtphys/2. ) THEN 
        IF ( ABS(MOD(localtime(i)+86400.*100,86400.) - overpasstime ) <= pdtphys/2. ) THEN 
                masque(i) = 1
        ENDIF
   ENDDO


END SUBROUTINE swathpolarsat

SUBROUTINE swathpolder(annee,mois,jour,heure,pdtphys, &
     rlon,rlat,masque)
!  Adapted from INCA

!  USE inca_dim
  USE dimphy
  IMPLICIT NONE

  !--Auteurs : Francois-Marie Breon + Olivier Boucher

  !-- adapted to be used in INCA aerosol module Michael Schulz

  ! not needed?

  INTEGER :: annee, mois, jour 
  REAL :: heure                      !--heure en jour
  REAL :: pdtphys                    !--pas de temps en seconde
  REAL :: rlon(1:klon), rlat(1:klon) !--longitude et latitude
  INTEGER :: masque(1:klon)

  REAL :: J0              !--origine des temps pour les orbites ADEOS 
  PARAMETER (J0=183.91267)
  REAL :: secinday
  PARAMETER (secinday=86400.)
  REAL :: duree_orb       !--Duree d une orbite ADEOS en jour 
  PARAMETER (duree_orb=6055.3715/secinday)
  REAL :: deltalon        !--Decalage en longitude entre 2 orbites successives
  PARAMETER (deltalon=41./585.*360.)
  REAL :: demi_larg_eq    !--demi-largeur d'une orbite a l equateur
  PARAMETER (demi_larg_eq=11.)
  REAL :: incli, inclideg
  PARAMETER (inclideg=98.59)
  REAL :: RADEG, DTOR, RPI
  REAL :: demi_periode
  INTEGER :: jacum(1:12)
  DATA jacum/0,31,59,90,120,151,181,212,243,273,304,334/
  INTEGER :: an, orb, i, j
  REAL :: timepolder, lon0, posnorm, lim_nord, lim_sud
  REAL :: tempo, lon_cen, demi_larg
  REAL :: lat_debut, lat_fin, lon_west, lon_east
  REAL :: zlon            !--rlon mais remis entre 0 et 360
  REAL :: deltat          !--plage de temps a considerer en jours  


  RPI = 4 * atan (1.0)
  deltat=pdtphys/secinday
  demi_periode=deltat/2./duree_orb
  RADEG=180./RPI
  DTOR=RPI/180.
  incli=inclideg*DTOR

  an = MOD(annee, 100) 
  timepolder=FLOAT((an-96)*365+jacum(mois)+jour)+heure
  orb=INT((timepolder-J0)/duree_orb+0.5) 
  lon0=360.-MOD(168.02+FLOAT(orb)*deltalon,360.)
  posnorm=(timepolder-j0)/duree_orb-FLOAT(orb)
  j=jacum(mois)+jour
  lim_nord=60.5 + 25.*(1.-COS((FLOAT(J)+10.)/365.*2.*RPI))
  lim_sud=-53.0 - 20.*(1.+COS(FLOAT(J)/365.*2*RPI))
  !--lat de debut
  lat_debut=MIN(MIN(90.,(-posnorm+demi_periode)*360.),lim_nord) 
  !--lat de fin
  lat_fin  =MAX(MAX(-90.,(-posnorm-demi_periode)*360.),lim_sud) 

  DO i=1, klon
     masque(i)=0
     tempo=ASIN( MAX(-1.,MIN(1., -SIN(rlat(i)*DTOR)/SIN(incli))))
     lon_cen=(ATAN(TAN(tempo)*COS(incli))-duree_orb*tempo)*RADEG
     demi_larg=demi_larg_eq/COS(rlat(i)*DTOR)
     IF (ABS(SIN(rlat(i)*DTOR)/SIN(incli)).GE.1.0) demi_larg=200.0
     IF (rlat(i).GE.lat_fin.AND.rlat(i).LE.lat_debut) THEN 
        IF (demi_larg.GE. 180.) THEN 
           masque(i)=1
        ELSE  
           lon_west = MOD(lon0+lon_cen-demi_larg+720., 360.)
           lon_east = MOD(lon0+lon_cen+demi_larg,      360.)
           zlon     = MOD(rlon(i)+360.,                360.)
           IF (lon_west.LE.lon_east) THEN 
              IF (zlon.GE.lon_west.AND.zlon.LE.lon_east) masque(i)=1
           ELSE               
              IF (zlon.GE.lon_west.OR.zlon.LE.lon_east) masque(i)=1
           ENDIF
        ENDIF
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE swathpolder


