      subroutine PHY_SISVAT_OUT

!------------------------------------------------------------------------------+
!                                                         Wed 26-Jun-2013  MAR |
!     SubRoutine PHY_SISVAT_OUT writes main variables from                     |
!                    SISVAT    (Soil                                           |
!                               Ice                                            |
!                               Snow                                           |
!                               Vegetation                                     |
!                               Atmosphere                                     |
!                               Transfer   Scheme)                             |
!                                                                              |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Mon 24-Jun-2013      |
!           Last Modification by H. Gallee,               Wed 26-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


! Global  Variables
! =================

      use Mod_Real
      use Mod_PHY____grd
      use Mod_SISVAT_grd
      use Mod_SISVAT_dat
      use Mod_SISVAT_kkl
      use Mod_SISVAT_flx
      use Mod_SISVAT_cdf




! Local   Variables
! =================

! Dummy
! -----

      integer  ::  i     ,j     ,nm    ,ikp
      integer  ::  isl




! From 1D to 2D Horizontal Grid
! =============================

        DO ikp=1,kcolp
           i   = ii__AP(ikp)
           j   = jj__AP(ikp)
        DO nm=1,mwp

! Energy Fluxes                                             (OUTPUT/NetCDF)
! ^^^^^^^^^^^^^                                              ^^^^^^^^^^^^^
           SOsoNC_xyn(i,j,nm)      = SOsoKL(ikp,nm)                     ! Absorb.Sol.Rad.            [W/m2]
           IRsoNC_xyn(i,j,nm)      = IRsoKL(ikp,nm)                     ! Absorb.IR  Rad.            [W/m2]
           HSsoNC_xyn(i,j,nm)      = HSsoKL(ikp,nm)                     ! HS                         [W/m2]
           HLsoNC_xyn(i,j,nm)      = HLsoKL(ikp,nm)                     ! HL                         [W/m2]
           HLs_NC_xyn(i,j,nm)      = HLs_KL(ikp,nm)                     ! Evaporation           [mm w.e./s]
           HLv_NC_xyn(i,j,nm)      = HLv_KL(ikp,nm)                     ! Transpiration              [W/m2]

           eta_NC_xyn(i,j,nm)      = 0.                                 !
         DO isl = -nsoil,0                                              !
           eta_NC_xyn(i,j,nm)      = eta_NC_xyn(i,j,nm)                &! Soil Moisture
     &    +eta_SV(ikp,nm,isl)      * dz_dSV(    isl)                    !
         END DO                                                         !

        END DO
        END DO




      end subroutine PHY_SISVAT_OUT
