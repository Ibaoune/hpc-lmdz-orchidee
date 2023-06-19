      module Mod_SISVATLTSo

!--------------------------------------------------------------------------+
!                                                     Wed 26-Jun-2013  MAR |
!     module Mod_SISVATLTSo contains local variables of SISVAT_TSo         |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!     version 3.p.4.1 created by H. Gallee,           Mon 17-Jun-2013      |
!           Last Modification by H. Gallee,           Wed 26-Jun-2013      |
!                                                                          |
!--------------------------------------------------------------------------+


! Global  Variables
! =================

      use Mod_Real



! Internal Variables
! ==================

      IMPLICIT NONE

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Tsisva
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  Fsisva       
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  dza__1       

      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  mu_sno        !     Snow thermal Conductivity
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  mu__dz        ! mu_(eta,sno)   / dz
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  dtC_sv        ! dt      / C
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  IRs__D        ! UpwardIR Previous Iter.Contr.
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  dIRsdT        ! UpwardIR           T Derivat.
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  f_HSHL        ! Factor common to HS and HL
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  dRidTs        ! d(Rib)/d(Ts)
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  HS___D        ! Sensible Heat Flux Atm.Contr.
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  f___HL        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  HL___D        ! Latent   Heat Flux Atm.Contr.
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  TSurf0        ! Previous Surface Temperature
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  qsatsg        ! Soil   Saturat. Spec. Humidity
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  dqs_dT        ! d(qsatsg)/dTv
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  Psi           ! 1st Soil Layer Water Potential
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  RHuSol        ! Soil Surface Relative Humidity
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  RHu_av        ! Soil Surface Relative Humidity
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Diag_A        ! A Diagonal
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Diag_B        ! B Diagonal
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Diag_C        ! C Diagonal
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Term_D        !   Independant Term
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Aux__P        ! P Auxiliary Variable
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:,:)  ::  Aux__Q        ! Q Auxiliary Variable
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  etaBAK        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  etaNEW        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  etEuBk        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  fac_dt        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  faceta        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  PsiArg        !
      real(kind=real8), SAVE, ALLOCATABLE, dimension(:,:)    ::  SHuSol        !


      end module Mod_SISVATLTSo
