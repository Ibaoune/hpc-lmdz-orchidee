      module Mod_CMiPhy_loc

!--------------------------------------------------------------------------+
!                                                     Mon 17-Jun-2013  MAR |
!     module Mod_CMiPhy_loc contains local variables of CMiPhy             |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Mon 17-Jun-2013      |
!           Last Modification by H. Gallee,           Mon 17-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


! Global  Variables
! =================

      use Mod_Real



! Internal Variables
! ==================

      IMPLICIT NONE


      real(kind=real8), SAVE,ALLOCATABLE,dimension(:)    :: qw_io0            ! Droplets   Concentration entering CMiPhy (for io)       [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:)    :: qi_io0            ! Ice  Part. Concentration entering CMiPhy (for io)       [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: qs___0            ! Snow Part. Concentration entering CMiPhy                [kg/kg]
! #qg real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: qg___0            ! Graupels   Concentration entering CMiPhy                [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: qr___0            ! Rain Drops Concentration entering CMiPhy                [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: Ta_dgC            ! Air   Temperature                                         [dgC]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: sqrrro            ! sqrt(roa(mzp)/roa(k))                                       [-]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: qsiEFF            ! EFFective Saturation Specific Humidity over Ice         [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: Fletch            ! Monodisperse Nb of hexagonal Plates, Fletcher (1962)        [-]

      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: lamdaS            ! Marshall-Palmer distribution parameter for Snow Particl.
! #qg real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: lamdaG            ! Marshall-Palmer distribution parameter for Graupels     
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: lamdaR            ! Marshall-Palmer distribution parameter for Rain Derops  
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: ps_ACW            ! Accretion of Cloud Drop. by Snow Particl.          Rate [kg/kg/s]

      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: ps_ACR            ! Accretion of Snow        by Rain                   Rate [kg/kg/s]

      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: FallVw            !     Sedimentation Velocity   of Droplets
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: FallVi            !     Sedimentation Velocity   of Ice  Particles
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: FallVs            !     Sedimentation Velocity   of Snow Particles
! #qg real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: FallVg            !     Sedimentation Velocity   of Snow Particles
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: FallVr            !     Sedimentation Velocity   of Rain Drops 
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:)    :: qwLoss            ! Mass Loss related to Sedimentation of Rain Droplets
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:)    :: qiLoss            ! Mass Loss related to Sedimentation of Ice  Crystals
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:)    :: qsLoss            ! Mass Loss related to Sedimentation of Snow Particles
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:)    :: qrLoss            ! Mass Loss related to Sedimentation of Rain Drops


!  Debug Variables
!  ---------------

! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wihm1             ! Cloud Droplets Freezing
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wihm2             ! Ice   Crystals Homogeneous Sublimation
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wicnd             ! Ice   Crystals Nucleation              (Emde & Kahlig)
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: widep             ! Ice   Crystals Growth Bergeron Process (Emde & Kahlig)
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wisub             ! Ice   Crystals             Sublimation (Levkov)
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wimlt             ! Ice   Crystals Melting  
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wwevp             ! Water Vapor Condensation / Evaporation (Fractional Cloudiness)
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wraut             ! Cloud Droplets AUTO-Conversion
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wsaut             ! Ice   Crystals AUTO-Conversion
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wracw             ! Accretion of Cloud Droplets by Rain, Ta > 0, --> Rain
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wsacw             ! Accretion of Cloud Droplets by Rain, Ta < 0, --> Snow
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wsaci             ! Accretion of Ice   Crystals by Snow          --> Snow
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wraci             ! Accretion of Ice   Crystals by Rain          --> Snow
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wiacr             ! Accretion of Rain by Ice   Crystals          --> Snow
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wsacr             ! Accretion of Rain by Snow                    --> Snow
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wracs             ! Accretion of Snow by Rain                    --> Snow, Rain
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wrevp             ! Rain  Drops     Evaporation  
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wssub             ! Snow  Particles Sublimation
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wsmlt             ! Snow  Particles Melting 
! #WH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: wsfre             ! Rain  Drops     Freezing

! #wH real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  :: debugV            ! Debug Variable (of 16 microphysical processes)


      end module Mod_CMiPhy_loc
