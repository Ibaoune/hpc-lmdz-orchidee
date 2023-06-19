!
! $Id: splaeropt_5wv_rrtm.F90 2644 2016-10-02 16:55:08Z oboucher $
!

SUBROUTINE SPLAEROPT_5WV_RRTM(  &
   zdm, zdh, tr_seri, RHcl,     &
   tausum, tau )

  USE DIMPHY
  USE aero_mod
  USE infotrac_phy
  USE phys_local_var_mod, ONLY: od550aer,od865aer,ec550aer,od550lt1aer
  !
  ! Olivier Boucher Jan 2017
  ! Based on Mie routines on ciclad CMIP6
  !
  IMPLICIT NONE
  !
  ! Input arguments:
  !
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: zdh      !--m
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: zdm      !--kg/m2
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: RHcl     ! humidite relative ciel clair
  REAL, DIMENSION(klon,klev,nbtr), INTENT(IN) :: tr_seri
  !
  ! Output arguments:
  !
  REAL, DIMENSION(klon,nwave,naero_tot), INTENT(OUT)      :: tausum
  REAL, DIMENSION(klon,klev,nwave,naero_tot), INTENT(OUT) :: tau
  !
  ! Local
  !
  INTEGER, PARAMETER :: las = nwave_sw
  LOGICAL :: soluble
  
  INTEGER :: i, k, m, itr, irh, aerindex
  INTEGER :: spsol, spinsol, la
  INTEGER :: RH_num(klon,klev)
  INTEGER, PARAMETER :: la443 = 1
  INTEGER, PARAMETER :: la550 = 2
  INTEGER, PARAMETER :: la670 = 3
  INTEGER, PARAMETER :: la765 = 4
  INTEGER, PARAMETER :: la865 = 5
  INTEGER, PARAMETER :: nbre_RH=12
  INTEGER, PARAMETER :: naero_soluble=2
  INTEGER, PARAMETER :: naero_insoluble=2
  INTEGER, PARAMETER :: naero=naero_soluble+naero_insoluble

  REAL, PARAMETER :: RH_tab(nbre_RH)=(/0.,10.,20.,30.,40.,50.,60.,70.,80.,85.,90.,95./)
  REAL, PARAMETER :: RH_MAX=95.
  REAL :: delta(klon,klev), rh(klon,klev)
  REAL :: tau_ae5wv_int   ! Intermediate computation of epaisseur optique aerosol
  REAL :: zrho
  CHARACTER*20 modname
  
  ! Soluble components 1-accumulation mode soluble; 2- seasalt coarse
  REAL :: alpha_aers_5wv(nbre_RH,las,naero_soluble)   ! Ext. coeff. ** m2/g 
  ! Insoluble components 1- Dust: 2- BC; 3- POM
  REAL :: alpha_aeri_5wv(las,naero_insoluble)         ! Ext. coeff. ** m2/g 
  !
  ! Proprietes optiques
  !
  REAL :: fact_RH(nbre_RH)

! From here on we look at the optical parameters at 5 wavelengths:  
! 443nm, 550, 670, 765 and 865 nm 
!                                   le 12 AVRIL 2006 
!  
 DATA alpha_aers_5wv/ & 
                           ! accumulation mode (sulfate+2% bc) soluble 
       4.632, 4.632, 4.632, 4.632, 6.206, 6.827, 7.616, 8.716,10.514,12.025,14.688,21.539, &
       3.981, 3.981, 3.981, 3.981, 5.346, 5.923, 6.662, 7.704, 9.437,10.914,13.562,20.591, &
       3.265, 3.265, 3.265, 3.265, 4.400, 4.909, 5.565, 6.500, 8.081, 9.449,11.943,18.791, &
       2.761, 2.761, 2.761, 2.761, 3.731, 4.182, 4.767, 5.606, 7.041, 8.294,10.610,17.118, &
       2.307, 2.307, 2.307, 2.307, 3.129, 3.522, 4.034, 4.774, 6.052, 7.180, 9.286,15.340, &

                           ! seasalt seasalt Coarse Soluble (CS)      
       0.576, 0.690, 0.738, 0.789, 0.855, 0.935, 1.046, 1.212, 1.512, 1.785, 2.258, 3.449, &
       0.595, 0.713, 0.763, 0.814, 0.880, 0.963, 1.079, 1.248, 1.550, 1.826, 2.306, 3.507, &
       0.617, 0.738, 0.789, 0.842, 0.911, 0.996, 1.113, 1.286, 1.592, 1.871, 2.369, 3.562, &
       0.632, 0.755, 0.808, 0.862, 0.931, 1.018, 1.140, 1.316, 1.626, 1.909, 2.409, 3.622, &
       0.645, 0.771, 0.825, 0.880, 0.951, 1.039, 1.164, 1.344, 1.661, 1.948, 2.455, 3.682  /

  DATA alpha_aeri_5wv/ &
                                 ! coarse dust insoluble 
       0.605, 0.611, 0.661, 0.714, 0.760, &
                                 ! super coarse insoluble 
       0.153, 0.156, 0.158, 0.157, 0.161 /
  ! 
  ! Initialisations
  tausum(:,:,:)=0.
  tau(:,:,:,:)=0.

  modname='splaeropt_5wv_rrtm'

  IF (naero.GT.naero_tot) THEN 
    CALL abort_physic(modname,'Too many aerosol types',1)
  ENDIF

  DO irh=1,nbre_RH-1
    fact_RH(irh)=1./(RH_tab(irh+1)-RH_tab(irh))
  ENDDO
   
  DO k=1, klev
    DO i=1, klon
      rh(i,k)=MIN(RHcl(i,k)*100.,RH_MAX)
      RH_num(i,k) = INT( rh(i,k)/10. + 1.)
      IF (rh(i,k).GT.85.) RH_num(i,k)=10
      IF (rh(i,k).GT.90.) RH_num(i,k)=11
      delta(i,k)=(rh(i,k)-RH_tab(RH_num(i,k)))*fact_RH(RH_num(i,k))
    ENDDO
  ENDDO

  DO itr=1,nbtr    !--loop over tracers  

    IF (tname(itr+nqo)=='PREC') THEN       !--fine mode accumulation mode
      CYCLE
    ELSE IF (tname(itr+nqo)=='FINE') THEN  !--fine mode accumulation mode
      soluble=.TRUE.
      spsol=1
      aerindex=1
    ELSE IF (tname(itr+nqo)=='COSS') THEN  !--coarse mode sea salt
      soluble=.TRUE.
      spsol=2
      aerindex=2
    ELSE IF (tname(itr+nqo)=='CODU') THEN  !--coarse mode dust
      soluble=.FALSE.
      spinsol=1
      aerindex=3
    ELSE IF (tname(itr+nqo)=='SCDU') THEN  !--super coarse mode dust
      soluble=.FALSE.
      spinsol=2
      aerindex=4
    ELSE
       CALL abort_physic(modname,'I cannot do aerosol optics for '//tname(itr+nqo),1)
    ENDIF

    DO la=1,las

      !--only 550 and 865 nm are used
      IF (la.NE.la550.AND.la.NE.la865) CYCLE

      IF (soluble) THEN  !--soluble aerosol with RH dependence

        DO k=1, klev
          DO i=1, klon
            tau_ae5wv_int = alpha_aers_5wv(RH_num(i,k),la,spsol)+DELTA(i,k)* &
                           (alpha_aers_5wv(RH_num(i,k)+1,la,spsol) - & 
                            alpha_aers_5wv(RH_num(i,k),la,spsol))
            tau(i,k,la,aerindex) = tr_seri(i,k,itr)*zdm(i,k)*tau_ae5wv_int
            tausum(i,la,aerindex)=tausum(i,la,aerindex)+tau(i,k,la,aerindex)
          ENDDO
        ENDDO

      ELSE               !--cases of insoluble aerosol

        DO k=1, klev
          DO i=1, klon
            tau_ae5wv_int = alpha_aeri_5wv(la,spinsol)
            tau(i,k,la,aerindex) = tr_seri(i,k,itr)*zdm(i,k)*tau_ae5wv_int
            tausum(i,la,aerindex)= tausum(i,la,aerindex)+tau(i,k,la,aerindex)
          ENDDO
        ENDDO

      ENDIF

    ENDDO   ! Boucle sur les longueurs d'onde
  ENDDO     ! Boucle sur les masses de traceurs

!--AOD calculations for diagnostics
  od550aer(:)=SUM(tausum(:,la550,1:naero),dim=2)
  od865aer(:)=SUM(tausum(:,la865,1:naero),dim=2)

!--extinction coefficient for diagnostic
  ec550aer(:,:)=SUM(tau(:,:,la550,1:naero),dim=3)/zdh(:,:)

!--aod for particles lower than 1 micron
  od550lt1aer(:)=tausum(:,la550,1)+tausum(:,la550,2)*0.3+tausum(:,la550,3)*0.2

END SUBROUTINE SPLAEROPT_5WV_RRTM
