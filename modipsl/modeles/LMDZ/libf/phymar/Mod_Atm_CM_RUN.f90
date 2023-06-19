      module Mod_Atm_CM_RUN

!--------------------------------------------------------------------------+
!                                                     Mon 17-Jun-2013  MAR |
!     module Mod_Atm_CM_RUN contains local variables of PHY_Atm_CM_RUN     |
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


      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  ::  qv__00           !  Air   Specific Humidity  before Cloud Microphys.  [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  ::  qw__00           !  Cloud Droplets Concentr. before Cloud Microphys.  [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  ::  qi__00           !  Cloud Crystals Concentr. before Cloud Microphys.  [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  ::  qs__00           !  Snow  Particl. Concentr. before Cloud Microphys.  [kg/kg]
! #qg real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  ::  qg__00           !  Graupels       Concentr. before Cloud Microphys.  [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  ::  qr__00           !  Rain  Drops    Concentr. before Cloud Microphys.  [kg/kg]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  ::  CFra00           !  Cloud Fraction           before Cloud Microphys.      [-]
! #cw real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  ::  CCNw00           !  Cloud Condens. Nuclei    before Cloud Microphys.      [-]
      real(kind=real8), SAVE,ALLOCATABLE,dimension(:,:)  ::  CCNi00           !  Cloud Ice      Nuclei    before Cloud Microphys.      [-]


      end module Mod_Atm_CM_RUN
