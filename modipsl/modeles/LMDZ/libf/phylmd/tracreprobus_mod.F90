MODULE tracreprobus_mod
!
! This module prepares and calls the Reprobus main subroutine 
!

CONTAINS

  SUBROUTINE tracreprobus(pdtphys, gmtime, debutphy, julien, &
       presnivs, xlat, xlon, pphis, pphi, &
       t_seri, pplay, paprs, sh , &
       tr_seri)

    USE dimphy
    USE infotrac_phy
#ifdef REPROBUS 
    USE CHEM_REP, ONLY : pdt_rep, &  ! pas de temps reprobus
         daynum, iter, &             ! jourjulien, iteration chimie
         pdel,&
         d_q_rep,d_ql_rep,d_qi_rep
#endif
    IMPLICIT NONE

! Input argument
!---------------
    REAL,INTENT(IN)    :: pdtphys    ! Pas d'integration pour la physique (seconde)
    REAL,INTENT(IN)    :: gmtime     ! Heure courante
    LOGICAL,INTENT(IN) :: debutphy   ! le flag de l'initialisation de la physique
    INTEGER,INTENT(IN) :: julien     ! Jour julien

    REAL,DIMENSION(klev),INTENT(IN)        :: presnivs! pressions approximat. des milieux couches (en PA)
    REAL,DIMENSION(klon),INTENT(IN)        :: xlat    ! latitudes pour chaque point 
    REAL,DIMENSION(klon),INTENT(IN)        :: xlon    ! longitudes pour chaque point 
    REAL,DIMENSION(klon),INTENT(IN)        :: pphis   ! geopotentiel du sol
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pphi    ! geopotentiel de chaque couche

    REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
    REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs   ! pression pour chaque inter-couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: sh      ! humidite specifique   


! Output argument
!----------------
    REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT)  :: tr_seri ! Concentration Traceur [U/KgA]  
  

! Local variables
!----------------
    INTEGER :: it, k, niter

#ifdef REPROBUS
    !   -- CHIMIE REPROBUS --
!    pdt_rep=pdtphys/2.
    niter=pdtphys/pdt_rep
    write(*,*)'nb d appel de REPROBUS',niter
    
    DO k = 1, klev
       pdel(:,k) = paprs(:,k) - paprs (:,k+1)
    END DO
    
    ! initialisation de ozone passif a ozone en debut d hiver HN et HS
    IF (julien == 341 .OR. julien == 181) THEN
       tr_seri(:,:,11)=tr_seri(:,:,8)
    END IF

    d_q_rep(:,:)  =0.
    d_ql_rep(:,:) =0.
    d_qi_rep(:,:) =0.
    
    DO  iter = 1,niter
       daynum = FLOAT(julien) + gmtime + (iter-1)*pdt_rep/86400.
       
!       DO it=1, nbtr
!     WRITE(lunout,*)it,' ',minval(tr_seri(:,:,it)),maxval(tr_seri(:,:,it))
! seulement pour les especes chimiques (pas l'age de l'air)
! verif valeurs extremes
! correction: a 1.e-30 quand =0 ou negatif et
! call abort si >ou= 1.e10
!          WRITE(*,*)it,'nqtot',nqtot,'nbtr',nbtr
!          IF (it < nqtot) THEN
!             WRITE(*,*)'iciav',it,nqtot
!#ifdef REPROBUS
!             CALL minmaxqfi_chimie(it,tr_seri(1,1,it),0.,1.e10,'avant chimie ')
!#endif
!             WRITE(*,*)iter,'avpres'
!          ENDIF
!       ENDDO
       
#ifdef REPROBUS
       CALL chemmain_rlong_1401( &
            tr_seri, & !argument phytrac (change de nom apres: vmr)
            xlon,    & !argument phytrac (change de nom apres: lon)
            xlat,    & !argument phytrac (change de nom apres: lat)
            t_seri,  & !argument phytrac (meme nom)
            pplay,   & !argument phytrac (meme nom)
            pphi,    & !argument phytrac (meme nom)
            pphis,   & !argument phytrac (meme nom)
            presnivs, & !argument phytrac (meme nom)
            sh,      & !argument phytrac (meme nom)
            debutphy) !argument phytrac (change de nom apres: debut)
       ! pdel, pdt_rep, daynum : definit dans phytrac et utilise dans chemmain
       !                 et transporte par CHEM_REP

!       DO it=1, nbtr
!     WRITE(lunout,*)it,' ',minval(tr_seri(:,:,it)),maxval(tr_seri(:,:,it))
! seulement pour les especes chimiques (pas l'age de l'air)
! verif valeurs extremes
! correction: a 1.e-30 quand =0 ou negatif et
! call abort si >ou= 1.e10
!          WRITE(*,*)it,'nqtot',nqtot,'nbtr',nbtr
!          IF (it < nqtot) THEN
!             WRITE(*,*)'iciap',it,nqtot
!             CALL minmaxqfi_chimie(it,tr_seri(1,1,it),0.,1.e10,'apres chemmain')
!             WRITE(*,*)iter,'appres'
!          ENDIF
!       ENDDO

#endif        
       
    END DO
    

    ! 
    DO it=1,nbtr
       WRITE(solsym(it),'(i2)') it
    END DO
#endif
  END SUBROUTINE tracreprobus

END MODULE tracreprobus_mod
