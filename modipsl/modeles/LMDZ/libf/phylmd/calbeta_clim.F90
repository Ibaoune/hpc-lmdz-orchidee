!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/calbeta.F90,v 1.2 2007/06/22 12:49:51
! fairhead Exp $
!

SUBROUTINE calbeta_clim(klon,time,lat_radian,beta)

  !======================================================================
  ! Auteur(s): A.K. TRAORE
  !======================================================================

  !USE phys_local_var_mod, ONLY : ideal_beta !pour faire la variable dans le
  ! physiq.f pour des sorties directes de beta

  USE phys_cal_mod, only: year_len
  USE print_control_mod, ONLY: prt_level

  implicit none
  integer klon,nt,j,it
  real logbeta(klon),pi
  real lat(klon),lat_radian(klon)
  integer time
  real time_radian
  real lat_sahel,beta(klon)
  real lat_nord,lat_sud

  !==============================================

  pi=2.*asin(1.)
  beta=0.

  !calcul des cordonnees

  ! print*,'LATITUDES BETA ',lat_radian
  time_radian=(time+15.)*2.*pi / year_len

  if (prt_level >= 1) print *, 'time_radian time', time_radian, time

  lat(:)=180.*lat_radian(:)/pi !lat(:)=lat_radian(:)

  lat_sahel=-5*sin(time_radian)+13
  lat_nord=lat_sahel+25.
  lat_sud=lat_sahel-25.
  do j=1,klon
     !===========
     if (lat(j) < 5. ) then

        logbeta(j)=0.2*(lat(j)-lat_sud)-1.6
        beta(j)=10**(logbeta(j))
        beta(j)=max(beta(j),0.03)
        beta(j)=min(beta(j),0.22)
        ! print*,'j,lat,lat_radian,beta',j,lat(j),lat_radian(j),beta(j)
        !===========
     elseif (lat(j) < 22.) then !lat(j)<22.

        logbeta(j)=-0.25*(lat(j)-lat_sahel)-1.6
        beta(j)=10**(logbeta(j))
        beta(j)=max(beta(j),1.e-2)
        beta(j)=min(beta(j),0.22)
        ! print*,'j,lat,lat_radian,beta',j,lat(j),lat_radian(j),beta(j)
        !===========
     else
        logbeta(j)=0.25*(lat(j)-lat_nord)-1.
        beta(j)=10**(logbeta(j))
        beta(j)=max(beta(j),1.e-2)
        beta(j)=min(beta(j),0.25)
        ! print*,'j,lat,lat_radian,beta',j,lat(j),lat_radian(j),beta(j)
     endif
     !===========
  enddo

end SUBROUTINE calbeta_clim
