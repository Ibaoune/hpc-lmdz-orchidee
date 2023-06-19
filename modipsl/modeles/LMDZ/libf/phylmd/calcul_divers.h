!
! $Id: calcul_divers.h 3435 2019-01-22 15:21:59Z fairhead $
!
!
! Initialisations diverses au tout debut
      IF(itap.EQ.1) THEN
         itapm1=0
      ENDIF

! Initialisation debut de mois
      IF(itap.EQ.itapm1+1) THEN
        nday_rain(:)=0.
!       print*,'initialisation mois suivants day_rain itap',itap
      ENDIF

! Calcul fin de journee : total_rain, nday_rain
      IF(MOD(itap,NINT(un_jour/phys_tstep)).EQ.0) THEN
!        print*,'calcul nday_rain itap ',itap
         DO i = 1, klon
            total_rain(i)=rain_fall(i)+snow_fall(i)  
            IF(total_rain(i).GT.0.) nday_rain(i)=nday_rain(i)+1.
         ENDDO
      ENDIF

! Initialisation fin de mois
      IF(MOD(itap-itapm1,NINT(mth_len*un_jour/phys_tstep)).EQ.0) THEN
        itapm1=itapm1+NINT(mth_len*un_jour/phys_tstep)
!       print*,'initialisation itapm1 ',itapm1
      ENDIF
!
! calcul temperatures minimale et maximale moyennees sur le mois 
!
!initialisation debut de mois ou de journee pour les fichiers mensuels
  IF(itap.EQ.itapm1+1) THEN
     t2m_min_mon=0.
     t2m_max_mon=0.
  ENDIF
  IF(MOD(itap,NINT(un_jour/phys_tstep)).EQ.1) THEN
     zt2m_min_mon=zt2m
     zt2m_max_mon=zt2m
  ENDIF
!calcul a chaque pas de temps pour les fichiers mensuels
     DO i = 1, klon
        zt2m_min_mon(i)=MIN(zt2m(i),zt2m_min_mon(i))
        zt2m_max_mon(i)=MAX(zt2m(i),zt2m_max_mon(i))
     ENDDO
!fin de journee
  IF(MOD(itap,NINT(un_jour/phys_tstep)).EQ.0) THEN
   t2m_min_mon=t2m_min_mon+zt2m_min_mon 
   t2m_max_mon=t2m_max_mon+zt2m_max_mon 
  ENDIF
!fin mois
  IF(itap==itapm1) THEN
   t2m_min_mon=t2m_min_mon/mth_len
   t2m_max_mon=t2m_max_mon/mth_len
  ENDIF
!
