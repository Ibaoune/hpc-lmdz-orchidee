      SUBROUTINE aeropt_spl(zdz, tr_seri, RHcl,
     .              id_prec, id_fine, id_coss, id_codu, id_scdu,
     .              ok_chimeredust,
     .                      ztaue550,ztaue670,ztaue865,
     .                      taue550_tr2,taue670_tr2,taue865_tr2,
     .                      taue550_ss,taue670_ss,taue865_ss,
     .                      taue550_dust,taue670_dust,taue865_dust,
     .             taue550_dustsco,taue670_dustsco,taue865_dustsco)
c
      USE dimphy
      USE infotrac
      IMPLICIT none
c
#include "chem.h"
#include "dimensions.h"
c #include "dimphy.h"
#include "YOMCST.h"
c
c Arguments:
c
c======================== INPUT ==================================      
      REAL zdz(klon,klev)
      REAL tr_seri(klon,klev,nbtr) ! masse of tracer 
      REAL RHcl(klon,klev)     ! humidite relativen ciel clair
      INTEGER  id_prec, id_fine, id_coss, id_codu, id_scdu
      LOGICAL ok_chimeredust
c============================== OUTPUT =================================
      REAL ztaue550(klon) ! epaisseur optique aerosol 550 nm
      REAL ztaue670(klon) ! epaisseur optique aerosol 670 nm
      REAL ztaue865(klon) ! epaisseur optique aerosol 865 nm
      REAL taue550_tr2(klon) ! epaisseur optique aerosol 550 nm, diagnostic
      REAL taue670_tr2(klon) ! epaisseur optique aerosol 670 nm, diagnostic
      REAL taue865_tr2(klon) ! epaisseur optique aerosol 865 nm, diagnostic
      REAL taue550_ss(klon) ! epaisseur optique aerosol 550 nm, diagnostic
      REAL taue670_ss(klon) ! epaisseur optique aerosol 670 nm, diagnostic
      REAL taue865_ss(klon) ! epaisseur optique aerosol 865 nm, diagnostic
      REAL taue550_dust(klon) ! epaisseur optique aerosol 550 nm, diagnostic
      REAL taue670_dust(klon) ! epaisseur optique aerosol 670 nm, diagnostic
      REAL taue865_dust(klon) ! epaisseur optique aerosol 865 nm, diagnostic
      REAL taue550_dustsco(klon) ! epaisseur optique aerosol 550 nm, diagnostic
      REAL taue670_dustsco(klon) ! epaisseur optique aerosol 670 nm, diagnostic
      REAL taue865_dustsco(klon) ! epaisseur optique aerosol 865 nm, diagnostic
c===================== LOCAL VARIABLES ===========================      
      INTEGER nb_lambda,nbre_RH
      PARAMETER (nb_lambda=3,nbre_RH=12)
      INTEGER i, k, RH_num
      REAL rh, RH_MAX, DELTA, RH_tab(nbre_RH)
      PARAMETER (RH_MAX=95.)
      INTEGER rh_int
      PARAMETER (rh_int=12)
      REAL auxreal
c      REAL ss_a(nb_lambda,int,nbtr-1)
c      DATA ss_a/72*1./
      REAL ss_dust(nb_lambda), ss_acc550(rh_int), alpha_acc
      REAL ss_dustsco(nb_lambda)
      REAL ss_acc670(rh_int), ss_acc865(rh_int)
      REAL ss_ssalt550(rh_int)
      REAL ss_ssalt670(rh_int), ss_ssalt865(rh_int)
      REAL burden_ss(klon)
      DATA ss_acc550 /3.135,3.135,3.135, 3.135, 4.260, 4.807,
     .                5.546,6.651,8.641,10.335,13.534,22.979/
      DATA ss_acc670 /2.220,2.220,2.220, 2.220, 3.048, 3.460,
     .                4.023,4.873,6.426, 7.761,10.322,18.079/
      DATA ss_acc865 /1.329,1.329,1.329, 1.329, 1.855, 2.124,
     .                2.494,3.060,4.114, 5.033, 6.831,12.457/  
!old4tracers      DATA ss_dust/0.564, 0.614, 0.700/ !for bin 0.5-10um
!      DATA ss_dust/0.553117, 0.610185, 0.7053460 / !for bin 0.5-3um radius
!      DATA ss_dustsco/0.1014, 0.102156, 0.1035538 / !for bin 3-15um radius
!20140902      DATA ss_dust/0.5345737, 0.5878828, 0.6772957/ !for bin 0.5-3um radius
!20140902      DATA ss_dustsco/0.1009634, 0.1018700, 0.1031178/ !for bin 3-15um radius
!3days      DATA ss_dust/0.4564216, 0.4906738, 0.5476248/ !for bin 0.5-3um radius
!3days      DATA ss_dustsco/0.1015022, 0.1024051, 0.1036622/ !for bin 3-15um radius
!JE20140911      DATA ss_dust/0.5167768,0.5684330,0.6531643/ !for bin 0.5-3um radius
!JE20140911      DATA ss_dustsco/0.1003391,0.1012288,0.1024651/ !for bin 3-15um radius
!JE20140915      DATA ss_dust/0.3188754,0.3430106,0.3829019/ !for bin 0.5-5um radius
!JE20140915      DATA ss_dustsco/8.0582686E-02,8.1255026E-02,8.1861295E-02/ !for bin 5-15um radius

!      DATA ss_dust/0.5167768,0.5684330,0.6531643/ !for bin 0.5-3um radius
!      DATA ss_dustsco/0.1003391,0.1012288,0.1024651/ !for bin 3-15um radius


      DATA ss_ssalt550/0.182,0.182,0.182,0.182,0.366,0.430,
     .                 0.484,0.551,0.648,0.724,0.847,1.218/ !for bin 0.5-20 um, fit_v2
      DATA ss_ssalt670/0.193,0.193,0.193,0.193,0.377,0.431,
     .                 0.496,0.587,0.693,0.784,0.925,1.257/ !for bin 0.5-20 um
      DATA ss_ssalt865/0.188,0.188,0.188,0.188,0.384,0.443,
     .                 0.502,0.580,0.699,0.799,0.979,1.404/ !for bin 0.5-20 um

      DATA RH_tab/0.,10.,20.,30.,40.,50.,60.,70.,80.,85.,90.,95./
c
      IF (ok_chimeredust) THEN 
!JE20150212<< : changes in ustar in dustmod changes emission distribution
!      ss_dust=(/0.5167768,0.5684330,0.6531643/)
!      ss_dustsco=(/0.1003391,0.1012288,0.1024651/)
! JE20150618: Change in dustmodule, div3 is now =6: change distributions
! div3=3      ss_dust   =(/0.4670522 , 0.5077308 , 0.5745184/)
! div3=3      ss_dustsco=(/0.099858  , 0.1007395 , 0.1019673/)
      ss_dust   =(/0.4851232 , 0.5292494 , 0.5935509/)
      ss_dustsco=(/0.1001981 , 0.1011043 , 0.1023113/)

!JE20150212>>

      ELSE
      ss_dust=(/0.564, 0.614, 0.700/)
      ss_dustsco=(/0.,0.,0./)
      ENDIF

      DO i=1, klon
        ztaue550(i)=0.0
        ztaue670(i)=0.0
        ztaue865(i)=0.0
        taue550_tr2(i)=0.0
        taue670_tr2(i)=0.0
        taue865_tr2(i)=0.0
        taue550_ss(i)=0.0
        taue670_ss(i)=0.0
        taue865_ss(i)=0.0
        taue550_dust(i)=0.0
        taue670_dust(i)=0.0
        taue865_dust(i)=0.0
        taue550_dustsco(i)=0.0
        taue670_dustsco(i)=0.0
        taue865_dustsco(i)=0.0
        burden_ss(i)=0.0
      ENDDO

      DO k=1, klev
      DO i=1, klon
c     
      rh=MIN(RHcl(i,k)*100.,RH_MAX)
      RH_num = INT( rh/10. + 1.)
      IF (rh.gt.85.) RH_num=10
      IF (rh.gt.90.) RH_num=11
c      IF (rh.gt.40.) THEN 
c          RH_num=5   ! Added by NHL temporarily
c          print *,'TEMPORARY CASE'
c      ENDIF
      DELTA=(rh-RH_tab(RH_num))/(RH_tab(RH_num+1)-RH_tab(RH_num))


c*******************************************************************
c                       AOD at 550 NM
c*******************************************************************
        alpha_acc=ss_acc550(RH_num) + DELTA*(ss_acc550(RH_num+1)-
     .            ss_acc550(RH_num))              !--m2/g
cnhl_test TOTAL AOD
       auxreal=0.
      IF(id_fine>0) auxreal=auxreal+alpha_acc*tr_seri(i,k,id_fine)
      IF(id_coss>0) auxreal=auxreal+ss_ssalt550(RH_num)*
     .                       tr_seri(i,k,id_coss)
      IF(id_codu>0) auxreal=auxreal+ss_dust(1)*tr_seri(i,k,id_codu)
      IF(id_scdu>0) auxreal=auxreal+ss_dustsco(1)*tr_seri(i,k,id_scdu)
      ztaue550(i)=ztaue550(i)+auxreal*zdz(i,k)*1.e6
    
!JE20150128        ztaue550(i)=ztaue550(i)+(alpha_acc*tr_seri(i,k,id_fine)+
!     .                 ss_ssalt550(RH_num)*tr_seri(i,k,id_coss)+
!     .                 ss_dust(1)*tr_seri(i,k,id_codu)+
!     .              ss_dustsco(1)*tr_seri(i,k,id_scdu)  )*zdz(i,k)*1.e6

cnhl_test TOTAL AOD IS NOW AOD COARSE MODE ONLY
cnhl_test        ztaue550(i)=ztaue550(i)+(                     
cnhl_test     .                 ss_ssalt550(RH_num)*tr_seri(i,k,3)+
cnhl_test     .                 ss_dust(1)*tr_seri(i,k,4))*zdz(i,k)*1.e6

        IF(id_fine>0) taue550_tr2(i)=taue550_tr2(i)
     .               + alpha_acc*tr_seri(i,k,id_fine)*zdz(i,k)*1.e6
        IF(id_coss>0) taue550_ss(i)=taue550_ss(i)+
     .                ss_ssalt550(RH_num)*tr_seri(i,k,id_coss)*
     .                zdz(i,k)*1.e6
        IF(id_codu>0) taue550_dust(i)=taue550_dust(i)+
     .                ss_dust(1)*tr_seri(i,k,id_codu)*
     .                zdz(i,k)*1.e6
        IF(id_scdu>0) taue550_dustsco(i)=taue550_dustsco(i)+
     .                ss_dustsco(1)*tr_seri(i,k,id_scdu)*
     .                zdz(i,k)*1.e6
!        print *,'taue550_ss = ',SUM(taue550_ss),MINVAL(taue550_ss),
!     .                                          MAXVAL(taue550_ss)

c*******************************************************************
c                       AOD at 670 NM
c*******************************************************************
        alpha_acc=ss_acc670(RH_num) + DELTA*(ss_acc670(RH_num+1)-
     .            ss_acc670(RH_num))              !--m2/g
      auxreal=0.
      IF(id_fine>0) auxreal=auxreal+alpha_acc*tr_seri(i,k,id_fine)
      IF(id_coss>0) auxreal=auxreal+ss_ssalt670(RH_num)
     .                      *tr_seri(i,k,id_coss)
      IF(id_codu>0) auxreal=auxreal+ss_dust(2)*tr_seri(i,k,id_codu)
      IF(id_scdu>0) auxreal=auxreal+ss_dustsco(2)*tr_seri(i,k,id_scdu)
      ztaue670(i)=ztaue670(i)+auxreal*zdz(i,k)*1.e6

!JE20150128        ztaue670(i)=ztaue670(i)+(alpha_acc*tr_seri(i,k,id_fine)+
!     .                 ss_ssalt670(RH_num)*tr_seri(i,k,id_coss)+
!     .                 ss_dust(2)*tr_seri(i,k,id_codu)+
!     .               ss_dustsco(2)*tr_seri(i,k,id_scdu))*zdz(i,k)*1.e6

      IF(id_fine>0)  taue670_tr2(i)=taue670_tr2(i)+
     .                alpha_acc*tr_seri(i,k,id_fine)*
     .                 zdz(i,k)*1.e6
      IF(id_coss>0)  taue670_ss(i)=taue670_ss(i)+
     .                ss_ssalt670(RH_num)*tr_seri(i,k,id_coss)*
     .                zdz(i,k)*1.e6
      IF(id_codu>0)  taue670_dust(i)=taue670_dust(i)
     .                +ss_dust(2)*tr_seri(i,k,id_codu)*
     .                zdz(i,k)*1.e6
      IF(id_scdu>0)  taue670_dustsco(i)=taue670_dustsco(i)+
     .                ss_dustsco(2)*tr_seri(i,k,id_scdu)*
     .                zdz(i,k)*1.e6

c*******************************************************************
c                       AOD at 865 NM
c*******************************************************************
        alpha_acc=ss_acc865(RH_num) + DELTA*(ss_acc865(RH_num+1)-
     .            ss_acc865(RH_num))              !--m2/g
        auxreal=0.
      IF(id_fine>0) auxreal=auxreal+alpha_acc*tr_seri(i,k,id_fine)
      IF(id_coss>0) auxreal=auxreal
     .                     +ss_ssalt865(RH_num)*tr_seri(i,k,id_coss)
      IF(id_codu>0) auxreal=auxreal+ss_dust(3)*tr_seri(i,k,id_codu)
      IF(id_scdu>0) auxreal=auxreal+ss_dustsco(3)*tr_seri(i,k,id_scdu)
        ztaue865(i)=ztaue865(i)+auxreal*zdz(i,k)*1.e6
!JE20150128        ztaue865(i)=ztaue865(i)+(alpha_acc*tr_seri(i,k,id_fine)+
!     .                 ss_ssalt865(RH_num)*tr_seri(i,k,id_coss)+
!     .                 ss_dust(3)*tr_seri(i,k,id_codu)+
!     .               ss_dustsco(3)*tr_seri(i,k,id_scdu))*zdz(i,k)*1.e6
      IF(id_fine>0) taue865_tr2(i)=taue865_tr2(i)
     .                +alpha_acc*tr_seri(i,k,id_fine)*
     .                 zdz(i,k)*1.e6
      IF(id_coss>0) taue865_ss(i)=taue865_ss(i)+
     .                ss_ssalt865(RH_num)*tr_seri(i,k,id_coss)*
     .                zdz(i,k)*1.e6
      IF(id_codu>0)  taue865_dust(i)=taue865_dust(i)
     .                +ss_dust(3)*tr_seri(i,k,id_codu)*
     .                zdz(i,k)*1.e6
      IF(id_scdu>0)  taue865_dustsco(i)=taue865_dustsco(i)+
     .                ss_dustsco(3)*tr_seri(i,k,id_scdu)*
     .                zdz(i,k)*1.e6


c
      IF(id_coss>0)  burden_ss(i)=burden_ss(i)
     .                +tr_seri(i,k,id_coss)*1.e6*1.e3*zdz(i,k)
      ENDDO            !-loop on klev
      ENDDO            !-loop on klon
!      print *,'tr_seri = ',SUM(tr_seri(:,:,3)),MINVAL(tr_seri(:,:,3)),
!     .                                          MAXVAL(tr_seri(:,:,3))
!      print *,'taue550_ss = ',SUM(taue550_ss),MINVAL(taue550_ss),
!     .                                        MAXVAL(taue550_ss)
c
      RETURN
      END
