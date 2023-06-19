
! $Id: aeropt.F90 2346 2015-08-21 15:13:46Z emillour $

SUBROUTINE aeropt(pplay, paprs, t_seri, msulfate, rhcl, tau_ae, piz_ae, &
    cg_ae, ai)

  USE dimphy
  IMPLICIT NONE



  include "YOMCST.h"

  ! Arguments:

  REAL, INTENT (IN) :: paprs(klon, klev+1)
  REAL, INTENT (IN) :: pplay(klon, klev), t_seri(klon, klev)
  REAL, INTENT (IN) :: msulfate(klon, klev) ! masse sulfate ug SO4/m3  [ug/m^3]
  REAL, INTENT (IN) :: rhcl(klon, klev) ! humidite relative ciel clair
  REAL, INTENT (OUT) :: tau_ae(klon, klev, 2) ! epaisseur optique aerosol
  REAL, INTENT (OUT) :: piz_ae(klon, klev, 2) ! single scattering albedo aerosol
  REAL, INTENT (OUT) :: cg_ae(klon, klev, 2) ! asymmetry parameter aerosol
  REAL, INTENT (OUT) :: ai(klon) ! POLDER aerosol index

  ! Local

  INTEGER i, k, inu
  INTEGER rh_num, nbre_rh
  PARAMETER (nbre_rh=12)
  REAL rh_tab(nbre_rh)
  REAL rh_max, delta, rh
  PARAMETER (rh_max=95.)
  DATA rh_tab/0., 10., 20., 30., 40., 50., 60., 70., 80., 85., 90., 95./
  REAL zrho, zdz
  REAL taue670(klon) ! epaisseur optique aerosol absorption 550 nm
  REAL taue865(klon) ! epaisseur optique aerosol extinction 865 nm
  REAL alpha_aer_sulfate(nbre_rh, 5) !--unit m2/g SO4
  REAL alphasulfate

  CHARACTER (LEN=20) :: modname = 'aeropt'
  CHARACTER (LEN=80) :: abort_message


  ! Proprietes optiques

  REAL alpha_aer(nbre_rh, 2) !--unit m2/g SO4
  REAL cg_aer(nbre_rh, 2)
  DATA alpha_aer/.500130E+01, .500130E+01, .500130E+01, .500130E+01, &
    .500130E+01, .616710E+01, .826850E+01, .107687E+02, .136976E+02, &
    .162972E+02, .211690E+02, .354833E+02, .139460E+01, .139460E+01, &
    .139460E+01, .139460E+01, .139460E+01, .173910E+01, .244380E+01, &
    .332320E+01, .440120E+01, .539570E+01, .734580E+01, .136038E+02/
  DATA cg_aer/.619800E+00, .619800E+00, .619800E+00, .619800E+00, &
    .619800E+00, .662700E+00, .682100E+00, .698500E+00, .712500E+00, &
    .721800E+00, .734600E+00, .755800E+00, .545600E+00, .545600E+00, &
    .545600E+00, .545600E+00, .545600E+00, .583700E+00, .607100E+00, &
    .627700E+00, .645800E+00, .658400E+00, .676500E+00, .708500E+00/
  DATA alpha_aer_sulfate/4.910, 4.910, 4.910, 4.910, 6.547, 7.373, 8.373, &
    9.788, 12.167, 14.256, 17.924, 28.433, 1.453, 1.453, 1.453, 1.453, 2.003, &
    2.321, 2.711, 3.282, 4.287, 5.210, 6.914, 12.305, 4.308, 4.308, 4.308, &
    4.308, 5.753, 6.521, 7.449, 8.772, 11.014, 12.999, 16.518, 26.772, 3.265, &
    3.265, 3.265, 3.265, 4.388, 5.016, 5.775, 6.868, 8.745, 10.429, 13.457, &
    22.538, 2.116, 2.116, 2.116, 2.116, 2.882, 3.330, 3.876, 4.670, 6.059, &
    7.327, 9.650, 16.883/

  DO i = 1, klon
    taue670(i) = 0.0
    taue865(i) = 0.0
  END DO

  DO k = 1, klev
    DO i = 1, klon
      IF (t_seri(i,k)==0) WRITE (*, *) 'aeropt T ', i, k, t_seri(i, k)
      IF (pplay(i,k)==0) WRITE (*, *) 'aeropt p ', i, k, pplay(i, k)
      zrho = pplay(i, k)/t_seri(i, k)/rd ! kg/m3
      zdz = (paprs(i,k)-paprs(i,k+1))/zrho/rg ! m
      rh = min(rhcl(i,k)*100., rh_max)
      rh_num = int(rh/10.+1.)
      IF (rh<0.) THEN
        abort_message = 'aeropt: RH < 0 not possible'
        CALL abort_physic(modname, abort_message, 1)
      END IF
      IF (rh>85.) rh_num = 10
      IF (rh>90.) rh_num = 11
      delta = (rh-rh_tab(rh_num))/(rh_tab(rh_num+1)-rh_tab(rh_num))

      inu = 1
      tau_ae(i, k, inu) = alpha_aer(rh_num, inu) + delta*(alpha_aer(rh_num+1, &
        inu)-alpha_aer(rh_num,inu))
      tau_ae(i, k, inu) = tau_ae(i, k, inu)*msulfate(i, k)*zdz*1.E-6
      piz_ae(i, k, inu) = 1.0
      cg_ae(i, k, inu) = cg_aer(rh_num, inu) + delta*(cg_aer(rh_num+1,inu)- &
        cg_aer(rh_num,inu))

      inu = 2
      tau_ae(i, k, inu) = alpha_aer(rh_num, inu) + delta*(alpha_aer(rh_num+1, &
        inu)-alpha_aer(rh_num,inu))
      tau_ae(i, k, inu) = tau_ae(i, k, inu)*msulfate(i, k)*zdz*1.E-6
      piz_ae(i, k, inu) = 1.0
      cg_ae(i, k, inu) = cg_aer(rh_num, inu) + delta*(cg_aer(rh_num+1,inu)- &
        cg_aer(rh_num,inu))
      ! jq
      ! jq for aerosol index

      alphasulfate = alpha_aer_sulfate(rh_num, 4) + &
        delta*(alpha_aer_sulfate(rh_num+1,4)-alpha_aer_sulfate(rh_num,4)) !--m2/g

      taue670(i) = taue670(i) + alphasulfate*msulfate(i, k)*zdz*1.E-6

      alphasulfate = alpha_aer_sulfate(rh_num, 5) + &
        delta*(alpha_aer_sulfate(rh_num+1,5)-alpha_aer_sulfate(rh_num,5)) !--m2/g

      taue865(i) = taue865(i) + alphasulfate*msulfate(i, k)*zdz*1.E-6

    END DO
  END DO

  DO i = 1, klon
    ai(i) = (-log(max(taue670(i),0.0001)/max(taue865(i), &
      0.0001))/log(670./865.))*taue865(i)
  END DO

  RETURN
END SUBROUTINE aeropt
