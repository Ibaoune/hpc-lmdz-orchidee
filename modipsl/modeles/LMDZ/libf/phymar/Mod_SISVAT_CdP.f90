     module Mod_SISVAT_CdP

!--------------------------------------------------------------------------+
!    module Mod_SISVAT_CdP                            Fri  1-Feb-2013  MAR |
!    module Mod_SISVAT_CdP  contains specific (Col de Porte) constants of  |
!                Soil/Ice Snow Vegetation Atmosphere Transfer Scheme       |
!                                                                          |
!--------------------------------------------------------------------------+



! General Variables
! =================

      use Mod_Real


      IMPLICIT NONE


! Col de Porte specific Constants
! ===============================

! Fractions of total solar irradiances in 3 spectral intervals
! ------------------------------------------------------------

      real(kind=real8), SAVE  ::  Dr_1SN = 0.59,   Dr_2SN = 0.31,   Dr_3SN = 0.10  ! Direct  Radiation
      real(kind=real8), SAVE  ::  Df_1SN = 0.95,   Df_2SN = 0.05,   Df_3SN = 0.00  ! Diffuse Radiation, Clear  Sky
      real(kind=real8), SAVE  ::  Dfc1SN = 0.66,   Dfc2SN = 0.27,   Dfc3SN = 0.07  ! Diffuse Radiation, Cloudy Sky
!                           0.3--0.8micr.m   0.8--1.5micr.m   1.5--2.8micr.m
!                           Fractions of total solar irradiance in 3 spectral intervals
!                          (see Eric Martin Sept. 1996, CROCUS, Subroutine METEO)

      real(kind=real8), SAVE  ::  DirSol,DifSol,TotSol,Clouds



     end module Mod_SISVAT_CdP
