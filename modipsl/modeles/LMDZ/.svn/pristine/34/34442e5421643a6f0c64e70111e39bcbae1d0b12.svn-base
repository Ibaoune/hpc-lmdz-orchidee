      subroutine PHY________OUT(MOTIVATION)

!------------------------------------------------------------------------------+
!                                                         Sat 22-Jun-2013  MAR |
!   MAR          PHY________OUT                                                |
!     subroutine PHY________OUT OUTPUTs    MAR PHYsical parameterizations      |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 12-Mar-2013      |
!           Last Modification by H. Gallee,               Sat 22-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_PHY_DY_kkl
      use Mod_PHY_CM_kkl
      use Mod_PHY_AT_kkl
      use Mod_SISVAT_gpt


      IMPLICIT NONE


      character(len=50)                     ::  MOTIVATION
      integer                               ::       i,     j,   ikl            !
      integer                               ::       k,    kz                   !
      real(kind=real8)                      ::  RelHum




! OUTPUT
! ======


! txt file
! --------

      i      = i_x0
      j      = j_y0
      ikl    = ikl0

      write(4,50)       MOTIVATION
 50   format(3(/,1x),'  OUTPUT for VERIFICATION: ',a50,               &
     &         /,1x ,'  ***********************',/,1x)

      write(4,55) Day_TU,Mon_TU,YearTU,HourTU,MinuTU,Sec_TU,it_EXP
 55   format(3x,2(i2,'-'),i4,4x,3(i2,'-'),'  Simulation Iteration No ',i6)

      write(4,62)
      write(4,60)
 60   format(7x,'|','  sigma  ',' |','     Z   ',' |','    T    ',' |' &!
     &     ,'   U  ',' |','   V  ',' |',' RH   ',' |','  Qv   ',' |'   &!
     &                                ,'  Qw   ',' |','  Qi   ',' |'   &!
     &                                ,'  TKE  ',' |','  eps  ',' |'   &!
     &                                ,'  Kzh L',' |','  Kzh  ',' |'   &!
     &,/,'       |','         ',' |','   [km]  ',' |','   [K]   ',' |' &!
     &     ,' [m/s]',' |',' [m/s]',' |',' [-]  ',' |',' [g/kg]',' |'   &!
     &                                ,' [g/kg]',' |',' [g/kg]',' |'   &!
     &                                ,'[m2/s2]',' |','[m3/s2]',' |'   &!
     &                                ,' [m2/s]',' |',' [m2/s]',' |'   &!
     &,/,' ------+','---------','-+','---------','-+','---------','-+' &!
     &     ,'------','-+','------','-+','------','-+','-------','-+'   &!
     &                                ,'-------','-+','-------','-+'   &!
     &                                ,'-------','-+','-------','-+'   &!
     &                                ,'-------','-+','-------','-+')   !
      DO k = 1,mzp
      kz=min(k,mzp)
      RelHum =      min(un_1,qv__DY(ikl,k)   /max(qvswCM(ikl,k),eps6))
      write(4,61) k,sigma(k),Z___DY(ikl,k)*1.e-3 ,Ta__DY(ikl,k)        &!
     &                      ,Ua__DY(ikl,k)       ,Va__DY(ikl,k)        &!
     &                      ,RelHum       ,  1.e3*qv__DY(ikl,k)        &!
     &                 ,1.e3*qw__CM(ikl,kz), 1.e3*qi__CM(ikl,kz)       &!
     &                 ,     TKE_AT(ikl,kz),      eps_AT(ikl,kz)       &!
     &                 ,     Kzh_AT(ikl,kz),      Kzh0AT(ikl,kz)
 61   format(i6,' |',f9.6,' |',  f9.4,' |' ,  f9.3,' |' ,2(f6.1,' |')  &!
     &              ,f6.2,' |',  f7.3,' |' ,2(f7.3,' |')               &!
     &                        ,2(f7.2,' |'),2(f7.3,' |')              )
      IF (mod(k,15).EQ.0)                                           THEN
      write(4,62)
 62   format(                                                          &!
     &   ' ------+','---------','-+','---------','-+','---------','-+' &!
     &     ,'------','-+','------','-+','------','-+','-------','-+'   &!
     &                                ,'-------','-+','-------','-+'   &!
     &                                ,'-------','-+','-------','-+'   &!
     &                                ,'-------','-+','-------','-+')   !

      write(4,60)
      END IF
      ENDDO

      write(4,62)
      write(4,63) uts_SV_gpt(ikl),1.e3*uqs_SV_gpt(ikl),us__SV_gpt(ikl) &
     &          , HsenSV_gpt(ikl),     HLatSV_gpt(ikl)
 63   format(                                                          &!
     &       6x,' |', 9x ,' |', ( 9x ,' |'),' u*T*     |'              &!
     &                                     ,2( 6x ,' |'), ' u*q*  |'   &!
     &                        ,   7x ,' |' ,2( 7x ,' |'),'  u*    |'   &!
     &                        ,   7x ,' |' ,2( 7x ,' |')               &!
     &    ,/,6x,' |', 9x ,' |', ( 9x ,' |'),  f9.6,' |'                &!
     &                                     ,2( 6x ,' |'),f6.3,' |'     &!
     &                        ,   7x ,' |' ,2( 7x ,' |'),f7.2,' |'     &!
     &                        ,   7x ,' |' ,2( 7x ,' |')               &!
     &    ,/,6x,' |', 9x ,' |', ( 9x ,' |'),  f9.3,' |'                &!
     &                                     ,2( 6x ,' |'),f6.0,' |'     &!
     &                        ,   7x ,' |' ,2( 7x ,' |'), 7x ,' |'     &!
     &                        ,   7x ,' |' ,2( 7x ,' |')               &!
     &    ,/,6x,' |', 9x ,' |', ( 9x ,' |'),' [W/m2]   |'              &!
     &                                     ,2( 6x ,' |'), ' [W/m2]|'   &!
     &                        ,   7x ,' |' ,2( 7x ,' |'),'        |'   &!
     &                        ,   7x ,' |' ,2( 7x ,' |')            )

      k =      mzpp
      write(4,62)
      RelHum =      min(un_1,qv__DY(ikl,k)   /max(qvswCM(ikl,k),eps6))
      write(4,64) k,sigma(k),Z___DY(ikl,k)*1.e-3 ,Ta__DY(ikl,k)        &!
     &                      ,zer0                ,zer0                 &!
     &                      ,RelHum       ,  1.e3*qv__DY(ikl,k)
 64   format(i6,' |',f9.6,' |',  f9.4,' |' ,  f9.3,' |' ,2(f6.1,' |')  &!
     &              ,f6.2,' |',  f7.3,' |' ,2( 7x ,' |')               &!
     &                        ,2( 7x ,' |'),2( 7x ,' |')              )
      write(4,62)



! cdf file
! --------


      end subroutine PHY________OUT
