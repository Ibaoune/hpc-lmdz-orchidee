SUBROUTINE cv3a_compress(len, nloc, ncum, nd, ntra, compress, &
                         iflag1, nk1, icb1, icbs1, &
                         plcl1, tnk1, qnk1, gznk1, hnk1, unk1, vnk1, &
                         wghti1, pbase1, buoybase1, &
                         t1, q1, qs1, t1_wake, q1_wake, qs1_wake, s1_wake, &
                         u1, v1, gz1, th1, th1_wake, &
                         tra1, &
                         h1, lv1, lf1, cpn1, p1, ph1, tv1, tp1, tvp1, clw1, &
                         h1_wake, lv1_wake, lf1_wake, cpn1_wake, tv1_wake, &
                         sig1, w01, ptop21, &
                         Ale1, Alp1, omega1, &
                         iflag, nk, icb, icbs, &
                         plcl, tnk, qnk, gznk, hnk, unk, vnk, &
                         wghti, pbase, buoybase, &
                         t, q, qs, t_wake, q_wake, qs_wake, s_wake, &
                         u, v, gz, th, th_wake, &
                         tra, &
                         h, lv, lf, cpn, p, ph, tv, tp, tvp, clw, &
                         h_wake, lv_wake, lf_wake, cpn_wake, tv_wake, &
                         sig, w0, ptop2, &
                         Ale, Alp, omega)
  ! **************************************************************
  ! *
  ! CV3A_COMPRESS                                               *
  ! *
  ! *
  ! written by   : Sandrine Bony-Lena , 17/05/2003, 11.22.15    *
  ! modified by  : Jean-Yves Grandpeix, 23/06/2003, 10.28.09    *
  ! **************************************************************

  IMPLICIT NONE

  include "cv3param.h"

  ! inputs:
  INTEGER, INTENT (IN)                               :: len, nloc, nd, ntra
!jyg<
  LOGICAL, INTENT (IN)                               :: compress  ! compression is performed if compress is true
!>jyg
  INTEGER, DIMENSION (len), INTENT (IN)              :: iflag1, nk1, icb1, icbs1
  REAL, DIMENSION (len), INTENT (IN)                 :: plcl1, tnk1, qnk1, gznk1
  REAL, DIMENSION (len), INTENT (IN)                 :: hnk1, unk1, vnk1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: wghti1(len, nd)
  REAL, DIMENSION (len), INTENT (IN)                 :: pbase1, buoybase1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: t1, q1, qs1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: t1_wake, q1_wake, qs1_wake
  REAL, DIMENSION (len), INTENT (IN)                 :: s1_wake
  REAL, DIMENSION (len, nd), INTENT (IN)             :: u1, v1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: gz1, th1, th1_wake
  REAL, DIMENSION (len, nd,ntra), INTENT (IN)        :: tra1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: h1, lv1, lf1, cpn1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: p1
  REAL, DIMENSION (len, nd+1), INTENT (IN)           :: ph1(len, nd+1)
  REAL, DIMENSION (len, nd), INTENT (IN)             :: tv1, tp1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: tvp1, clw1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: h1_wake, lv1_wake, cpn1_wake
  REAL, DIMENSION (len, nd), INTENT (IN)             :: tv1_wake, lf1_wake
  REAL, DIMENSION (len, nd), INTENT (IN)             :: sig1, w01
  REAL, DIMENSION (len), INTENT (IN)                 :: ptop21
  REAL, DIMENSION (len), INTENT (IN)                 :: Ale1, Alp1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: omega1
!
  ! in/out 
  INTEGER, INTENT (INOUT)                            :: ncum
!
  ! outputs:
  ! en fait, on a nloc=len pour l'instant (cf cv_driver)
  INTEGER, DIMENSION (nloc), INTENT (OUT)            ::  iflag, nk, icb, icbs
  REAL, DIMENSION (nloc), INTENT (OUT)               ::  plcl, tnk, qnk, gznk
  REAL, DIMENSION (nloc), INTENT (OUT)               ::  hnk, unk, vnk
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  wghti
  REAL, DIMENSION (nloc), INTENT (OUT)               ::  pbase, buoybase
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  t, q, qs
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  t_wake, q_wake, qs_wake
  REAL, DIMENSION (nloc), INTENT (OUT)               ::  s_wake
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  u, v
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  gz, th, th_wake
  REAL, DIMENSION (nloc, nd,ntra), INTENT (OUT)      ::  tra
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  h, lv, lf, cpn
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  p
  REAL, DIMENSION (nloc, nd+1), INTENT (OUT)         ::  ph
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  tv, tp
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  tvp, clw
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  h_wake, lv_wake, cpn_wake
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  tv_wake, lf_wake
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  sig, w0
  REAL, DIMENSION (nloc), INTENT (OUT)               ::  ptop2
  REAL, DIMENSION (nloc), INTENT (OUT)               ::  Ale, Alp
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           ::  omega

  ! local variables:
  INTEGER i, k, nn, j

  CHARACTER (LEN=20) :: modname = 'cv3a_compress'
  CHARACTER (LEN=80) :: abort_message

!jyg<
  IF (compress) THEN
!>jyg

  DO k = 1, nl + 1
    nn = 0
    DO i = 1, len
      IF (iflag1(i)==0) THEN
        nn = nn + 1
        wghti(nn, k) = wghti1(i, k)
        t(nn, k) = t1(i, k)
        q(nn, k) = q1(i, k)
        qs(nn, k) = qs1(i, k)
        t_wake(nn, k) = t1_wake(i, k)
        q_wake(nn, k) = q1_wake(i, k)
        qs_wake(nn, k) = qs1_wake(i, k)
        u(nn, k) = u1(i, k)
        v(nn, k) = v1(i, k)
        gz(nn, k) = gz1(i, k)
        th(nn, k) = th1(i, k)
        th_wake(nn, k) = th1_wake(i, k)
        h(nn, k) = h1(i, k)
        lv(nn, k) = lv1(i, k)
        lf(nn, k) = lf1(i, k)
        cpn(nn, k) = cpn1(i, k)
        p(nn, k) = p1(i, k)
        ph(nn, k) = ph1(i, k)
        tv(nn, k) = tv1(i, k)
        tp(nn, k) = tp1(i, k)
        tvp(nn, k) = tvp1(i, k)
        clw(nn, k) = clw1(i, k)
        h_wake(nn, k) = h1_wake(i, k)
        lv_wake(nn, k) = lv1_wake(i, k)
        lf_wake(nn, k) = lf1_wake(i, k)
        cpn_wake(nn, k) = cpn1_wake(i, k)
        tv_wake(nn, k) = tv1_wake(i, k)
        sig(nn, k) = sig1(i, k)
        w0(nn, k) = w01(i, k)
        omega(nn, k) = omega1(i, k)
      END IF
    END DO
  END DO
!
  ! AC!      do 121 j=1,ntra
  ! AC!ccccc      do 111 k=1,nl+1
  ! AC!      do 111 k=1,nd
  ! AC!       nn=0
  ! AC!      do 101 i=1,len
  ! AC!      if(iflag1(i).eq.0)then
  ! AC!       nn=nn+1
  ! AC!       tra(nn,k,j)=tra1(i,k,j)
  ! AC!      endif
  ! AC! 101  continue
  ! AC! 111  continue
  ! AC! 121  continue

  IF (nn/=ncum) THEN
    PRINT *, 'WARNING nn not equal to ncum: ', nn, ncum
    abort_message = ''
    CALL abort_physic(modname, abort_message, 1)
  END IF

  nn = 0
  DO i = 1, len
    IF (iflag1(i)==0) THEN
      nn = nn + 1
      s_wake(nn) = s1_wake(i)
      iflag(nn) = iflag1(i)
      nk(nn) = nk1(i)
      icb(nn) = icb1(i)
      icbs(nn) = icbs1(i)
      plcl(nn) = plcl1(i)
      tnk(nn) = tnk1(i)
      qnk(nn) = qnk1(i)
      gznk(nn) = gznk1(i)
      hnk(nn) = hnk1(i)
      unk(nn) = unk1(i)
      vnk(nn) = vnk1(i)
      pbase(nn) = pbase1(i)
      buoybase(nn) = buoybase1(i)
      sig(nn, nd) = sig1(i, nd)
      ptop2(nn) = ptop2(i)
      Ale(nn) = Ale1(i)
      Alp(nn) = Alp1(i)
    END IF
  END DO

  IF (nn/=ncum) THEN
    PRINT *, 'WARNING nn not equal to ncum: ', nn, ncum
    abort_message = ''
    CALL abort_physic(modname, abort_message, 1)
  END IF
!
!jyg<
  ELSE  !(compress)
!
      ncum = len
!
      wghti(:,1:nl+1) = wghti1(:,1:nl+1)
      t(:,1:nl+1) = t1(:,1:nl+1)
      q(:,1:nl+1) = q1(:,1:nl+1)
      qs(:,1:nl+1) = qs1(:,1:nl+1)
      t_wake(:,1:nl+1) = t1_wake(:,1:nl+1)
      q_wake(:,1:nl+1) = q1_wake(:,1:nl+1)
      qs_wake(:,1:nl+1) = qs1_wake(:,1:nl+1)
      u(:,1:nl+1) = u1(:,1:nl+1)
      v(:,1:nl+1) = v1(:,1:nl+1)
      gz(:,1:nl+1) = gz1(:,1:nl+1)
      th(:,1:nl+1) = th1(:,1:nl+1)
      th_wake(:,1:nl+1) = th1_wake(:,1:nl+1)
      h(:,1:nl+1) = h1(:,1:nl+1)
      lv(:,1:nl+1) = lv1(:,1:nl+1)
      lf(:,1:nl+1) = lf1(:,1:nl+1)
      cpn(:,1:nl+1) = cpn1(:,1:nl+1)
      p(:,1:nl+1) = p1(:,1:nl+1)
      ph(:,1:nl+1) = ph1(:,1:nl+1)
      tv(:,1:nl+1) = tv1(:,1:nl+1)
      tp(:,1:nl+1) = tp1(:,1:nl+1)
      tvp(:,1:nl+1) = tvp1(:,1:nl+1)
      clw(:,1:nl+1) = clw1(:,1:nl+1)
      h_wake(:,1:nl+1) = h1_wake(:,1:nl+1)
      lv_wake(:,1:nl+1) = lv1_wake(:,1:nl+1)
      lf_wake(:,1:nl+1) = lf1_wake(:,1:nl+1)
      cpn_wake(:,1:nl+1) = cpn1_wake(:,1:nl+1)
      tv_wake(:,1:nl+1) = tv1_wake(:,1:nl+1)
      sig(:,1:nl+1) = sig1(:,1:nl+1)
      w0(:,1:nl+1) = w01(:,1:nl+1)
      omega(:,1:nl+1) = omega1(:,1:nl+1)
!
      s_wake(:) = s1_wake(:)
      iflag(:) = iflag1(:)
      nk(:) = nk1(:)
      icb(:) = icb1(:)
      icbs(:) = icbs1(:)
      plcl(:) = plcl1(:)
      tnk(:) = tnk1(:)
      qnk(:) = qnk1(:)
      gznk(:) = gznk1(:)
      hnk(:) = hnk1(:)
      unk(:) = unk1(:)
      vnk(:) = vnk1(:)
      pbase(:) = pbase1(:)
      buoybase(:) = buoybase1(:)
      sig(:, nd) = sig1(:, nd)
      ptop2(:) = ptop2(:)
      Ale(:) = Ale1(:)
      Alp(:) = Alp1(:)
!
  ENDIF !(compress)
!>jyg

  RETURN
END SUBROUTINE cv3a_compress
