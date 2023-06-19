!
! $Id: splaeropt_6bands_rrtm.F90 2644 2016-10-02 16:55:08Z oboucher $
!
SUBROUTINE SPLAEROPT_6BANDS_RRTM ( &
     zdm, tr_seri, RHcl,           &
     tau_allaer, piz_allaer, cg_allaer )

  USE dimphy
  USE aero_mod
  USE infotrac_phy
  USE phys_local_var_mod, ONLY: absvisaer

  ! Olivier Boucher Jan 2017
  ! based on Mie routines on ciclad CMIP6
  !
  IMPLICIT NONE

  INCLUDE "clesphys.h"
  !
  ! Input arguments:
  !
  REAL, DIMENSION(klon,klev),     INTENT(IN)  :: zdm        !--hauteur des couches en kg/m2
  REAL, DIMENSION(klon,klev),     INTENT(IN)  :: RHcl       ! humidite relative ciel clair
  REAL, DIMENSION(klon,klev,nbtr),INTENT(IN)  :: tr_seri
  !
  ! Output arguments:
  ! 2= total aerosols 
  ! 1= natural aerosols
  !
  REAL, DIMENSION(klon,klev,2,nbands_sw_rrtm), INTENT(OUT) :: tau_allaer ! epaisseur optique aerosol
  REAL, DIMENSION(klon,klev,2,nbands_sw_rrtm), INTENT(OUT) :: piz_allaer ! single scattering albedo aerosol
  REAL, DIMENSION(klon,klev,2,nbands_sw_rrtm), INTENT(OUT) :: cg_allaer  ! asymmetry parameter aerosol
  !
  ! Local
  !
  LOGICAL :: soluble
  INTEGER :: i, k, irh, itr, inu
  INTEGER :: aerindex, spsol, spinsol
  INTEGER :: RH_num(klon,klev)

  INTEGER, PARAMETER :: naero_soluble=2    ! 1- accumulation soluble; 2- coarse soluble
  INTEGER, PARAMETER :: naero_insoluble=2  ! 1- coarse dust; 2- supercoarse dust
  INTEGER, PARAMETER :: naero=naero_soluble+naero_insoluble

  INTEGER, PARAMETER :: nbre_RH=12
  REAL,PARAMETER :: RH_tab(nbre_RH)=(/0.,10.,20.,30.,40.,50.,60.,70.,80.,85.,90.,95./)
  REAL, PARAMETER :: RH_MAX=95.
  REAL :: delta(klon,klev), rh(klon,klev)
  REAL :: tau_ae2b_int   ! Intermediate computation of epaisseur optique aerosol
  REAL :: piz_ae2b_int   ! Intermediate computation of Single scattering albedo
  REAL :: cg_ae2b_int    ! Intermediate computation of Assymetry parameter
  REAL :: fact_RH(nbre_RH), tmp_var
  !
  REAL, DIMENSION(klon,klev,naero_tot,nbands_sw_rrtm) ::  tau_ae
  REAL, DIMENSION(klon,klev,naero_tot,nbands_sw_rrtm) ::  piz_ae
  REAL, DIMENSION(klon,klev,naero_tot,nbands_sw_rrtm) ::  cg_ae
  !
  ! Proprietes optiques
  !
  REAL:: alpha_aers_6bands(nbre_RH,nbands_sw_rrtm,naero_soluble)   !--unit m2/g aer
  REAL:: alpha_aeri_6bands(nbands_sw_rrtm,naero_insoluble)         !--unit m2/g aer
  REAL:: cg_aers_6bands(nbre_RH,nbands_sw_rrtm,naero_soluble)      !--unitless 
  REAL:: cg_aeri_6bands(nbands_sw_rrtm,naero_insoluble)            !--unitless
  REAL:: piz_aers_6bands(nbre_RH,nbands_sw_rrtm,naero_soluble)     !--unitless
  REAL:: piz_aeri_6bands(nbands_sw_rrtm,naero_insoluble)           !--unitless
  !
  REAL, PARAMETER :: tau_min = 1.e-8
  CHARACTER*20 modname

!***************************************************************************
!--the order of the soluble   species has to follow the spsol   index below
!--the order of the insoluble species has to follow the spinsol index below

  DATA alpha_aers_6bands/  & 
       ! accumulation (sulfate+2% bc) mode soluble
  5.212, 5.212, 5.212, 5.212, 6.973, 7.581, 8.349, 9.400,11.078,12.463,14.857,20.837, &
  4.906, 4.906, 4.906, 4.906, 6.568, 7.195, 7.989, 9.088,10.869,12.354,14.951,21.545, &
  3.940, 3.940, 3.940, 3.940, 5.291, 5.861, 6.591, 7.620, 9.332,10.791,13.410,20.370, &
  2.292, 2.292, 2.292, 2.292, 3.105, 3.493, 3.996, 4.724, 5.978, 7.084, 9.146,15.067, &
  0.762, 0.762, 0.762, 0.762, 1.050, 1.201, 1.401, 1.699, 2.232, 2.721, 3.675, 6.678, &
  0.090, 0.090, 0.090, 0.090, 0.122, 0.141, 0.166, 0.204, 0.275, 0.344, 0.484, 0.973, &
        ! coarse seasalt
  0.547, 0.657, 0.705, 0.754, 0.817, 0.896, 1.008, 1.169, 1.456, 1.724, 2.199, 3.358, &
  0.566, 0.679, 0.727, 0.776, 0.840, 0.920, 1.032, 1.196, 1.492, 1.760, 2.238, 3.416, &
  0.596, 0.714, 0.764, 0.816, 0.882, 0.965, 1.081, 1.250, 1.552, 1.828, 2.310, 3.509, &
  0.644, 0.771, 0.825, 0.880, 0.951, 1.040, 1.164, 1.345, 1.666, 1.957, 2.462, 3.700, &
  0.640, 0.772, 0.829, 0.887, 0.965, 1.061, 1.198, 1.398, 1.758, 2.085, 2.658, 4.031, &
  0.452, 0.562, 0.609, 0.659, 0.728, 0.813, 0.938, 1.125, 1.471, 1.797, 2.384, 3.855  /

  DATA alpha_aeri_6bands/  & 
       ! coarse dust insoluble
  0.594, 0.610, 0.619, 0.762, 0.791, 0.495, &
       ! supercoarse dust insoluble 
  0.151, 0.152, 0.155, 0.159, 0.168, 0.167  /

  DATA cg_aers_6bands/ &
       ! accumulation (sulfate+2% bc) mode soluble 
  0.692, 0.692, 0.692, 0.692, 0.735, 0.739, 0.744, 0.749, 0.755, 0.759, 0.765, 0.772, &
  0.690, 0.690, 0.690, 0.690, 0.736, 0.740, 0.746, 0.752, 0.760, 0.765, 0.771, 0.779, &
  0.678, 0.678, 0.678, 0.678, 0.727, 0.733, 0.740, 0.748, 0.759, 0.766, 0.775, 0.787, &
  0.641, 0.641, 0.641, 0.641, 0.692, 0.700, 0.710, 0.721, 0.736, 0.746, 0.760, 0.781, &
  0.553, 0.553, 0.553, 0.553, 0.603, 0.615, 0.627, 0.643, 0.664, 0.678, 0.699, 0.735, &
  0.343, 0.343, 0.343, 0.343, 0.386, 0.399, 0.414, 0.433, 0.460, 0.480, 0.510, 0.569, &
       ! seasalt coarse Soluble
  0.754, 0.770, 0.776, 0.781, 0.784, 0.791, 0.797, 0.805, 0.815, 0.822, 0.828, 0.840, &
  0.736, 0.753, 0.759, 0.765, 0.771, 0.778, 0.785, 0.793, 0.804, 0.811, 0.820, 0.831, &
  0.716, 0.735, 0.742, 0.748, 0.754, 0.762, 0.769, 0.778, 0.789, 0.796, 0.807, 0.819, &
  0.704, 0.725, 0.733, 0.739, 0.745, 0.752, 0.759, 0.768, 0.778, 0.784, 0.792, 0.803, &
  0.716, 0.737, 0.744, 0.751, 0.756, 0.763, 0.770, 0.777, 0.786, 0.790, 0.795, 0.800, &
  0.688, 0.730, 0.741, 0.751, 0.761, 0.771, 0.782, 0.795, 0.810, 0.820, 0.833, 0.849  /

  DATA cg_aeri_6bands/ &
       ! coarse dust insoluble
  0.801, 0.779, 0.709, 0.698, 0.710, 0.687, &
       ! super coarse dust insoluble
  0.862, 0.871, 0.852, 0.799, 0.758, 0.651  /

  DATA piz_aers_6bands/&
       ! accumulation (sulfate+2% bc) mode soluble
  0.941, 0.941, 0.941, 0.941, 0.958, 0.961, 0.965, 0.969, 0.974, 0.977, 0.981, 0.987, &
  0.941, 0.941, 0.941, 0.941, 0.959, 0.963, 0.967, 0.971, 0.976, 0.979, 0.983, 0.988, &
  0.953, 0.953, 0.953, 0.953, 0.967, 0.971, 0.974, 0.978, 0.982, 0.984, 0.988, 0.992, &
  0.955, 0.955, 0.955, 0.955, 0.969, 0.972, 0.976, 0.980, 0.984, 0.986, 0.989, 0.994, &
  0.936, 0.936, 0.936, 0.936, 0.955, 0.961, 0.966, 0.972, 0.978, 0.982, 0.987, 0.993, &
  0.792, 0.792, 0.792, 0.792, 0.848, 0.867, 0.887, 0.907, 0.931, 0.944, 0.960, 0.980, &
       ! seasalt coarse soluble
  1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.001, 1.000, &
  1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
  1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
  1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
  0.994, 0.994, 0.995, 0.995, 0.995, 0.995, 0.996, 0.996, 0.996, 0.996, 0.996, 0.996, &
  0.976, 0.867, 0.837, 0.814, 0.796, 0.774, 0.754, 0.735, 0.713, 0.702, 0.690, 0.675  /

  DATA piz_aeri_6bands/ &
       ! coarse dust insoluble
  0.866, 0.875, 0.915, 0.977, 0.993, 0.971, &
       ! super coarse dust insoluble
  0.789, 0.749, 0.791, 0.918, 0.970, 0.909  /
!
  spsol = 0
  spinsol = 0 

  modname='splaeropt_6bands_rrt'

  IF (NSW.NE.nbands_sw_rrtm) THEN
     CALL abort_physic(modname,'Erreur NSW doit etre egal a 6 pour cette routine',1)
  ENDIF

  DO irh=1,nbre_RH-1
    fact_RH(irh)=1./(RH_tab(irh+1)-RH_tab(irh))
  ENDDO
   
  DO k=1, klev
    DO i=1, klon
      rh(i,k)=MIN(RHcl(i,k)*100.,RH_MAX)
      RH_num(i,k) = INT(rh(i,k)/10. + 1.)
      IF (rh(i,k).GT.85.) RH_num(i,k)=10
      IF (rh(i,k).GT.90.) RH_num(i,k)=11
      delta(i,k)=(rh(i,k)-RH_tab(RH_num(i,k)))*fact_RH(RH_num(i,k))
    ENDDO
  ENDDO

  tau_ae(:,:,:,:)=0.
  piz_ae(:,:,:,:)=0.
  cg_ae(:,:,:,:)=0.
    
  DO itr=1, nbtr

    IF (tname(itr+nqo)=='PREC') THEN       !--precursor
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

    IF (soluble) THEN ! For aerosol soluble components

         DO k=1, klev
           DO i=1, klon

             tmp_var=tr_seri(i,k,itr)*zdm(i,k) !-- g/m2

             DO inu=1,NSW

               tau_ae2b_int= alpha_aers_6bands(RH_num(i,k),inu,spsol)+ & 
                             delta(i,k)* (alpha_aers_6bands(RH_num(i,k)+1,inu,spsol) - & 
                             alpha_aers_6bands(RH_num(i,k),inu,spsol))
                   
               piz_ae2b_int = piz_aers_6bands(RH_num(i,k),inu,spsol) + & 
                              delta(i,k)* (piz_aers_6bands(RH_num(i,k)+1,inu,spsol) - & 
                              piz_aers_6bands(RH_num(i,k),inu,spsol))
                   
               cg_ae2b_int = cg_aers_6bands(RH_num(i,k),inu,spsol) + & 
                             delta(i,k)* (cg_aers_6bands(RH_num(i,k)+1,inu,spsol) - & 
                             cg_aers_6bands(RH_num(i,k),inu,spsol))

               tau_ae(i,k,aerindex,inu)    = tmp_var*tau_ae2b_int
               piz_ae(i,k,aerindex,inu)    = piz_ae2b_int
               cg_ae(i,k,aerindex,inu)     = cg_ae2b_int

             ENDDO
           ENDDO
         ENDDO

    ELSE    ! For all aerosol insoluble components

       DO k=1, klev
         DO i=1, klon

           tmp_var=tr_seri(i,k,itr)*zdm(i,k) !-- g/m2

           DO inu=1,NSW
             tau_ae2b_int = alpha_aeri_6bands(inu,spinsol)
             piz_ae2b_int = piz_aeri_6bands(inu,spinsol)
             cg_ae2b_int = cg_aeri_6bands(inu,spinsol) 

             tau_ae(i,k,aerindex,inu) = tmp_var*tau_ae2b_int
             piz_ae(i,k,aerindex,inu) = piz_ae2b_int
             cg_ae(i,k,aerindex,inu)= cg_ae2b_int
           ENDDO
         ENDDO
       ENDDO

     ENDIF ! soluble / insoluble

  ENDDO  ! nbtr

!--all (natural + anthropogenic) aerosol
  tau_allaer(:,:,2,:)=SUM(tau_ae(:,:,1:naero,:),dim=3) 
  tau_allaer(:,:,2,:)=MAX(tau_allaer(:,:,2,:),tau_min)

  piz_allaer(:,:,2,:)=SUM(tau_ae(:,:,1:naero,:)*piz_ae(:,:,1:naero,:),dim=3)/tau_allaer(:,:,2,:)
  piz_allaer(:,:,2,:)=MIN(MAX(piz_allaer(:,:,2,:),0.01),1.0)
  WHERE (tau_allaer(:,:,2,:).LE.tau_min) piz_allaer(:,:,2,:)=1.0

  cg_allaer(:,:,2,:)=SUM(tau_ae(:,:,1:naero,:)*piz_ae(:,:,1:naero,:)*cg_ae(:,:,1:naero,:),dim=3)/  & 
                     (tau_allaer(:,:,2,:)*piz_allaer(:,:,2,:))
  cg_allaer(:,:,2,:)=MIN(MAX(cg_allaer(:,:,2,:),0.0),1.0)

!--no aerosol
  tau_allaer(:,:,1,:)=tau_min 
  piz_allaer(:,:,1,:)=1.0 
  cg_allaer(:,:,1,:)=0.0

!--waveband 2 and all aerosol (third index = 2)
  inu=2
  absvisaer(:)=SUM((1-piz_allaer(:,:,2,inu))*tau_allaer(:,:,2,inu),dim=2)

END SUBROUTINE SPLAEROPT_6BANDS_RRTM
