!
! $Id: $
!
      SUBROUTINE phyredem (fichnom)

      use iostart, only : open_restartphy,close_restartphy,
     &                    put_var,put_field
      use phys_state_var_mod, only : rlon,rlat

      implicit none

      character(len=*),intent(in) :: fichnom

      integer,parameter :: tab_cntrl_len=100
      real :: tab_cntrl(tab_cntrl_len)

! open file

      CALL open_restartphy(fichnom)

! tab_cntrl() contains run parameters

      tab_cntrl(:)=0.0
 

      CALL put_var("controle","Parametres de controle",tab_cntrl)

! coordinates

      CALL put_field("longitude",
     .               "Longitudes de la grille physique",rlon)
     
      CALL put_field("latitude","Latitudes de la grille physique",rlat)

! close file

      CALL close_restartphy
!$OMP BARRIER

      END
