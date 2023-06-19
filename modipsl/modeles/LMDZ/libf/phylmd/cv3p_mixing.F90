SUBROUTINE cv3p_mixing(nloc, ncum, nd, na, ntra, icb, nk, inb, &
                       ph, t, rr, rs, u, v, tra, h, lv, lf, frac, qta, &
                       unk, vnk, hp, tv, tvp, ep, clw, sig, &
                       Ment, Qent, hent, uent, vent, nent, &
                       Sigij, elij, supmax, Ments, Qents, traent)
! **************************************************************
! *
! CV3P_MIXING : compute mixed draught properties and,         *
! within a scaling factor, mixed draught        *
! mass fluxes.                                  *
! written by  : VTJ Philips,JY Grandpeix, 21/05/2003, 09.14.15*
! modified by :                                               *
! **************************************************************

  USE print_control_mod, ONLY: mydebug=>debug , lunout, prt_level
  USE ioipsl_getin_p_mod, ONLY: getin_p
  USE add_phys_tend_mod, ONLY: fl_cor_ebil

  IMPLICIT NONE

  include "cvthermo.h"
  include "cv3param.h"
  include "YOMCST2.h"
  include "cvflag.h"

!inputs:
  INTEGER, INTENT (IN)                               :: ncum, nd, na
  INTEGER, INTENT (IN)                               :: ntra, nloc
  INTEGER, DIMENSION (nloc), INTENT (IN)             :: icb, inb, nk
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: sig
  REAL, DIMENSION (nloc), INTENT (IN)                :: unk, vnk
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: qta
  REAL, DIMENSION (nloc, nd+1), INTENT (IN)          :: ph
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: t, rr, rs
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: u, v
  REAL, DIMENSION (nloc, nd, ntra), INTENT (IN)      :: tra ! input of convect3
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: lv
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: lf
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: frac !ice fraction in condensate
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: h  !liquid water static energy of environment
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: hp !liquid water static energy of air shed from adiab. asc.
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: tv, tvp
  REAL, DIMENSION (nloc, na), INTENT (IN)            :: ep, clw

!outputs:
  REAL, DIMENSION (nloc, na, na), INTENT (OUT)       :: Ment, Qent
  REAL, DIMENSION (nloc, na, na), INTENT (OUT)       :: uent, vent
  REAL, DIMENSION (nloc, na, na), INTENT (OUT)       :: Sigij, elij
  REAL, DIMENSION (nloc, na), INTENT (OUT)           :: supmax           ! Highest mixing fraction of mixed
                                                                         ! updraughts with the sign of (h-hp)
  REAL, DIMENSION (nloc, nd, nd, ntra), INTENT (OUT) :: traent
  REAL, DIMENSION (nloc, nd, nd), INTENT (OUT)       :: Ments, Qents
  REAL, DIMENSION (nloc, nd, nd), INTENT (OUT)       :: hent
  INTEGER, DIMENSION (nloc, nd), INTENT (OUT)        :: nent

!local variables:
  INTEGER i, j, k, il, im, jm
  INTEGER num1, num2
  REAL                               :: rti, bf2, anum, denom, dei, altem, cwat, stemp
  REAL                               :: alt, delp, delm
  REAL, DIMENSION (nloc)             :: Qmixmax, Rmixmax, sqmrmax
  REAL, DIMENSION (nloc)             :: Qmixmin, Rmixmin, sqmrmin
  REAL, DIMENSION (nloc)             :: signhpmh
  REAL, DIMENSION (nloc)             :: Sx
  REAL                               :: Scrit2
  REAL, DIMENSION (nloc)             :: Smid, Sjmin, Sjmax
  REAL, DIMENSION (nloc)             :: Sbef, sup, smin
  REAL, DIMENSION (nloc)             :: ASij, ASij_inv, smax, Scrit
  REAL, DIMENSION (nloc, nd, nd)     :: Sij
  REAL, DIMENSION (nloc, nd)         :: csum
  REAL                               :: awat
  REAL                               :: cpm        !Mixed draught heat capacity
  REAL                               :: Tm         !Mixed draught temperature
  LOGICAL, DIMENSION (nloc)          :: lwork

  REAL amxupcrit, df, ff
  INTEGER nstep

  INTEGER,SAVE                                       :: igout=1
!$OMP THREADPRIVATE(igout)

! --   Mixing probability distribution functions

  REAL Qcoef1, Qcoef2, QFF, QFFF, Qmix, Rmix, Qmix1, Rmix1, Qmix2, Rmix2, F

  Qcoef1(F) = tanh(F/gammas)
  Qcoef2(F) = (tanh(F/gammas)+gammas*log(cosh((1.-F)/gammas)/cosh(F/gammas)))
  QFF(F) = max(min(F,1.), 0.)
  QFFf(F) = min(QFF(F), scut)
  Qmix1(F) = (tanh((QFF(F)-Fmax)/gammas)+Qcoef1max)/Qcoef2max
  Rmix1(F) = (gammas*log(cosh((QFF(F)-Fmax)/gammas))+QFF(F)*Qcoef1max)/Qcoef2max
  Qmix2(F) = -log(1.-QFFf(F))/scut
  Rmix2(F) = (QFFf(F)+(1.-QFF(F))*log(1.-QFFf(F)))/scut
  Qmix(F) = qqa1*Qmix1(F) + qqa2*Qmix2(F)
  Rmix(F) = qqa1*Rmix1(F) + qqa2*Rmix2(F)

  INTEGER, SAVE :: ifrst
  DATA ifrst/0/
!$OMP THREADPRIVATE(ifrst)


! =====================================================================
! --- INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
! =====================================================================

! -- Initialize mixing PDF coefficients
  IF (ifrst==0) THEN
    ifrst = 1
    Qcoef1max = Qcoef1(Fmax)
    Qcoef2max = Qcoef2(Fmax)
!<jyg
   print*, 'fmax, gammas, qqa1, qqa2, Qcoef1max, Qcoef2max ', &
            fmax, gammas, qqa1, qqa2, Qcoef1max, Qcoef2max
!>jyg
!
  END IF


! ori        do 360 i=1,ncum*nlp
  DO j = 1, nl
    DO i = 1, ncum
      nent(i, j) = 0
! in convect3, m is computed in cv3_closure
! ori          m(i,1)=0.0
    END DO
  END DO

! ori      do 400 k=1,nlp
! ori       do 390 j=1,nlp
  DO j = 1, nl
    DO k = 1, nl
      DO i = 1, ncum
        Qent(i, k, j) = rr(i, j)
        uent(i, k, j) = u(i, j)
        vent(i, k, j) = v(i, j) 
        elij(i, k, j) = 0.0
        hent(i, k, j) = 0.0
!AC!            Ment(i,k,j)=0.0
!AC!            Sij(i,k,j)=0.0
      END DO
    END DO
  END DO

!AC!
  Ment(1:ncum, 1:nd, 1:nd) = 0.0
  Sij(1:ncum, 1:nd, 1:nd) = 0.0
!AC!
!ym
  Sigij(1:ncum, 1:nd, 1:nd) = 0.0
!ym 

!jyg!  DO k = 1, ntra
!jyg!    DO j = 1, nd ! instead nlp
!jyg!      DO i = 1, nd ! instead nlp
!jyg!        DO il = 1, ncum
!jyg!          traent(il, i, j, k) = tra(il, j, k)
!jyg!        END DO
!jyg!      END DO
!jyg!    END DO
!jyg!  END DO

! =====================================================================
! --- CALCULATE ENTRAINED AIR MASS FLUX (Ment), TOTAL WATER MIXING
! --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
! --- FRACTION (Sij)
! =====================================================================

  DO i = minorig + 1, nl

    IF (ok_entrain) THEN
      DO j = minorig, nl
        DO il = 1, ncum
          IF ((i>=icb(il)) .AND. (i<=inb(il)) .AND. (j>=(icb(il)-1)) &
                           .AND. (j<=inb(il))) THEN

!!            rti = qnk(il) - ep(il, i)*clw(il, i)
            rti = qta(il,i-1) - ep(il, i)*clw(il, i)
            bf2 = 1. + lv(il, j)*lv(il, j)*rs(il, j)/(rrv*t(il,j)*t(il,j)*cpd)
!jyg(from aj)<
            IF (cvflag_ice) THEN
! print*,cvflag_ice,'cvflag_ice dans do 700'
              IF (t(il,j)<=263.15) THEN
                bf2 = 1. + (lf(il,j)+lv(il,j))*(lv(il,j)+frac(il,j)* &
                     lf(il,j))*rs(il, j)/(rrv*t(il,j)*t(il,j)*cpd)
              END IF
            END IF
!>jyg
            anum = h(il, j) - hp(il, i) + (cpv-cpd)*t(il, j)*(rti-rr(il,j))
            denom = h(il, i) - hp(il, i) + (cpd-cpv)*(rr(il,i)-rti)*t(il, j)
            dei = denom
            IF (abs(dei)<0.01) dei = 0.01
            Sij(il, i, j) = anum/dei
            Sij(il, i, i) = 1.0
            altem = Sij(il, i, j)*rr(il, i) + (1.-Sij(il,i,j))*rti - rs(il, j)
            altem = altem/bf2
            cwat = clw(il, j)*(1.-ep(il,j))
            stemp = Sij(il, i, j)
            IF ((stemp<0.0 .OR. stemp>1.0 .OR. altem>cwat) .AND. j>i) THEN
!jyg(from aj)<
              IF (cvflag_ice) THEN
                anum = anum - (lv(il,j)+frac(il,j)*lf(il,j))*(rti-rs(il,j)-cwat*bf2)
                denom = denom + (lv(il,j)+frac(il,j)*lf(il,j))*(rr(il,i)-rti)
              ELSE
                anum = anum - lv(il, j)*(rti-rs(il,j)-cwat*bf2)
                denom = denom + lv(il, j)*(rr(il,i)-rti)
              END IF
!>jyg
              IF (abs(denom)<0.01) denom = 0.01
              Sij(il, i, j) = anum/denom
              altem = Sij(il, i, j)*rr(il, i) + (1.-Sij(il,i,j))*rti - rs(il, j)
              altem = altem - (bf2-1.)*cwat
            END IF
            IF (Sij(il,i,j)>0.0) THEN
!!!                 Ment(il,i,j)=m(il,i)
              Ment(il, i, j) = 1.
              elij(il, i, j) = altem
              elij(il, i, j) = amax1(0.0, elij(il,i,j))
              nent(il, i) = nent(il, i) + 1
            END IF

            Sij(il, i, j) = amax1(0.0, Sij(il,i,j))
            Sij(il, i, j) = amin1(1.0, Sij(il,i,j))
          ELSE IF (j > i) THEN
            IF (prt_level >= 10) THEN
              print *,'cv3p_mixing i, j, Sij given by the no-precip eq. ', i, j, Sij(il,i,j)
            ENDIF
          END IF ! new
        END DO
      END DO
    ELSE  ! (ok_entrain)
      DO il = 1,ncum
        nent(il,i) = 0
      ENDDO
    ENDIF ! (ok_entrain)

!jygdebug<
    IF (prt_level >= 10) THEN
      print *,'cv3p_mixing i, nent(i), icb, inb ',i, nent(igout,i), icb(igout), inb(igout)
      IF (nent(igout,i) .gt. 0) THEN
        print *,'i,(j,Sij(i,j),j=icb-1,inb) ',i,(j,Sij(igout,i,j),j=icb(igout)-1,inb(igout))
      ENDIF
    ENDIF
!>jygdebug

! ***   if no air can entrain at level i assume that updraft detrains  ***
! ***   at that level and calculate detrained air flux and properties  ***


! @      do 170 i=icb(il),inb(il)

    DO il = 1, ncum
      IF ((i>=icb(il)) .AND. (i<=inb(il)) .AND. (nent(il,i)==0)) THEN
! @      if(nent(il,i).eq.0)then
!!!       Ment(il,i,i)=m(il,i)
        Ment(il, i, i) = 1.
!!        Qent(il, i, i) = qnk(il) - ep(il, i)*clw(il, i)
        Qent(il, i, i) = qta(il,i-1) - ep(il, i)*clw(il, i)
        uent(il, i, i) = unk(il)
        vent(il, i, i) = vnk(il)
        IF (fl_cor_ebil .GE. 2) THEN
          hent(il, i, i) = hp(il,i)
        ENDIF
        elij(il, i, i) = clw(il, i)*(1.-ep(il,i))
        Sij(il, i, i) = 0.0
      END IF
    END DO
  END DO ! i = minorig + 1, nl

!jyg!  DO j = 1, ntra
!jyg!    DO i = minorig + 1, nl
!jyg!      DO il = 1, ncum
!jyg!        IF (i>=icb(il) .AND. i<=inb(il) .AND. nent(il,i)==0) THEN
!jyg!          traent(il, i, i, j) = tra(il, nk(il), j)
!jyg!        END IF
!jyg!      END DO
!jyg!    END DO
!jyg!  END DO

  DO j = minorig, nl
    DO i = minorig, nl
      DO il = 1, ncum
        IF ((j>=(icb(il)-1)) .AND. (j<=inb(il)) .AND. &
            (i>=icb(il)) .AND. (i<=inb(il))) THEN
          Sigij(il, i, j) = Sij(il, i, j)
        END IF
      END DO
    END DO
  END DO
! @      enddo

! @170   continue

! =====================================================================
! ---  NORMALIZE ENTRAINED AIR MASS FLUXES
! ---  TO REPRESENT EQUAL PROBABILITIES OF MIXING
! =====================================================================

  CALL zilch(csum, nloc*nd)

  DO il = 1, ncum
    lwork(il) = .FALSE.
  END DO

! ---------------------------------------------------------------
  DO i = minorig + 1, nl      !Loop on origin level "i"
! ---------------------------------------------------------------

    num1 = 0
    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il)) num1 = num1 + 1
    END DO
    IF (num1<=0) GO TO 789


!JYG1    Find maximum of SIJ for J>I, if any.

    Sx(:) = 0.

    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il)) THEN
        signhpmh(il) = sign(1., hp(il,i)-h(il,i))
        Sbef(il) = max(0., signhpmh(il))
      END IF
    END DO

    DO j = i + 1, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. j<=inb(il)) THEN
          IF (Sbef(il)<Sij(il,i,j)) THEN
            Sx(il) = max(Sij(il,i,j), Sx(il))
          END IF
          Sbef(il) = Sij(il, i, j)
        END IF
      END DO
    END DO


    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il)) THEN
        lwork(il) = (nent(il,i)/=0)
!!        rti = qnk(il) - ep(il, i)*clw(il, i)
        rti = qta(il,i-1) - ep(il, i)*clw(il, i)
!jyg<
        IF (cvflag_ice) THEN

          anum = h(il, i) - hp(il, i) - (lv(il,i)+frac(il,i)*lf(il,i))* &
                       (rti-rs(il,i)) + (cpv-cpd)*t(il, i)*(rti-rr(il,i))
          denom = h(il, i) - hp(il, i) + (lv(il,i)+frac(il,i)*lf(il,i))* &
                       (rr(il,i)-rti) + (cpd-cpv)*t(il, i)*(rr(il,i)-rti)
        ELSE

          anum = h(il, i) - hp(il, i) - lv(il, i)*(rti-rs(il,i)) + &
                       (cpv-cpd)*t(il, i)*(rti-rr(il,i))
          denom = h(il, i) - hp(il, i) + lv(il, i)*(rr(il,i)-rti) + &
                       (cpd-cpv)*t(il, i)*(rr(il,i)-rti)
        END IF
!>jyg
        IF (abs(denom)<0.01) denom = 0.01
        Scrit(il) = min(anum/denom, 1.)
        alt = rti - rs(il, i) + Scrit(il)*(rr(il,i)-rti)

!JYG1    Find new critical value Scrit2
!         such that : Sij > Scrit2  => mixed draught will detrain at J<I
!                     Sij < Scrit2  => mixed draught will detrain at J>I

        Scrit2 = min(Scrit(il), Sx(il))*max(0., -signhpmh(il)) + &
                 Scrit(il)*max(0., signhpmh(il))

        Scrit(il) = Scrit2

!JYG    Correction pour la nouvelle logique; la correction pour ALT
! est un peu au hazard
        IF (Scrit(il)<=0.0) Scrit(il) = 0.0
        IF (alt<=0.0) Scrit(il) = 1.0

        smax(il) = 0.0
        ASij(il) = 0.0
        sup(il) = 0.      ! upper S-value reached by descending draughts
      END IF
    END DO

! ---------------------------------------------------------------
    DO j = minorig, nl         !Loop on destination level "j"
! ---------------------------------------------------------------

      num2 = 0
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. &
            j>=(icb(il)-1) .AND. j<=inb(il) .AND. &
            lwork(il)) num2 = num2 + 1
      END DO
      IF (num2<=0) GO TO 175

! -----------------------------------------------
      IF (j>i) THEN
! -----------------------------------------------
        DO il = 1, ncum
          IF (i>=icb(il) .AND. i<=inb(il) .AND. &
              j>=(icb(il)-1) .AND. j<=inb(il) .AND. &
              lwork(il)) THEN
            IF (Sij(il,i,j)>0.0) THEN
              Smid(il) = min(Sij(il,i,j), Scrit(il))
              Sjmax(il) = Smid(il)
              Sjmin(il) = Smid(il)
              IF (Smid(il)<smin(il) .AND. Sij(il,i,j+1)<Smid(il)) THEN
                smin(il) = Smid(il)
                Sjmax(il) = min((Sij(il,i,j+1)+Sij(il,i,j))/2., Sij(il,i,j), Scrit(il))
                Sjmin(il) = max((Sbef(il)+Sij(il,i,j))/2., Sij(il,i,j))
                Sjmin(il) = min(Sjmin(il), Scrit(il))
                Sbef(il) = Sij(il, i, j)
              END IF
            END IF
          END IF
        END DO
! -----------------------------------------------
      ELSE IF (j==i) THEN
! -----------------------------------------------
        DO il = 1, ncum
          IF (i>=icb(il) .AND. i<=inb(il) .AND. &
              j>=(icb(il)-1) .AND. j<=inb(il) .AND. &
              lwork(il)) THEN
            IF (Sij(il,i,j)>0.0) THEN
              Smid(il) = 1.
              Sjmin(il) = max((Sij(il,i,j-1)+Smid(il))/2., Scrit(il))*max(0., -signhpmh(il)) + &
                          min((Sij(il,i,j+1)+Smid(il))/2., Scrit(il))*max(0., signhpmh(il))
              Sjmin(il) = max(Sjmin(il), sup(il))
              Sjmax(il) = 1.

! -             preparation des variables Scrit, Smin et Sbef pour la partie j>i
              Scrit(il) = min(Sjmin(il), Sjmax(il), Scrit(il))

              smin(il) = 1.
              Sbef(il) = max(0., signhpmh(il))
              supmax(il, i) = sign(Scrit(il), -signhpmh(il))
            END IF
          END IF
        END DO
! -----------------------------------------------
      ELSE IF (j<i) THEN
! -----------------------------------------------
        DO il = 1, ncum
          IF (i>=icb(il) .AND. i<=inb(il) .AND. &
              j>=(icb(il)-1) .AND. j<=inb(il) .AND. &
              lwork(il)) THEN
            IF (Sij(il,i,j)>0.0) THEN
              Smid(il) = max(Sij(il,i,j), Scrit(il))
              Sjmax(il) = Smid(il)
              Sjmin(il) = Smid(il)
              IF (Smid(il)>smax(il) .AND. Sij(il,i,j+1)>Smid(il)) THEN
                smax(il) = Smid(il)
                Sjmax(il) = max((Sij(il,i,j+1)+Sij(il,i,j))/2., Sij(il,i,j))
                Sjmax(il) = max(Sjmax(il), Scrit(il))
                Sjmin(il) = min((Sbef(il)+Sij(il,i,j))/2., Sij(il,i,j))
                Sjmin(il) = max(Sjmin(il), Scrit(il))
                Sbef(il) = Sij(il, i, j)
              END IF
              IF (abs(Sjmin(il)-Sjmax(il))>1.E-10) &
                             sup(il) = max(Sjmin(il), Sjmax(il), sup(il))
            END IF
          END IF
        END DO
! -----------------------------------------------
      END IF
! -----------------------------------------------


      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. &
            j>=(icb(il)-1) .AND. j<=inb(il) .AND. &
            lwork(il)) THEN
          IF (Sij(il,i,j)>0.0) THEN
!!            rti = qnk(il) - ep(il, i)*clw(il, i)
            rti = qta(il,i-1) - ep(il, i)*clw(il, i)
            Qmixmax(il) = Qmix(Sjmax(il))
            Qmixmin(il) = Qmix(Sjmin(il))
            Rmixmax(il) = Rmix(Sjmax(il))
            Rmixmin(il) = Rmix(Sjmin(il))
            sqmrmax(il) = Sjmax(il)*Qmix(Sjmax(il)) - Rmix(Sjmax(il))
            sqmrmin(il) = Sjmin(il)*Qmix(Sjmin(il)) - Rmix(Sjmin(il))

            Ment(il, i, j) = abs(Qmixmax(il)-Qmixmin(il))*Ment(il, i, j)

! Sigij(i,j) is the 'true' mixing fraction of mixture Ment(i,j)
            IF (abs(Qmixmax(il)-Qmixmin(il))>1.E-10) THEN
              Sigij(il, i, j) = (sqmrmax(il)-sqmrmin(il))/(Qmixmax(il)-Qmixmin(il))
            ELSE
              Sigij(il, i, j) = 0.
            END IF

! --    Compute Qent, uent, vent according to the true mixing fraction
            Qent(il, i, j) = (1.-Sigij(il,i,j))*rti     + Sigij(il, i, j)*rr(il, i)
            uent(il, i, j) = (1.-Sigij(il,i,j))*unk(il) + Sigij(il, i, j)*u(il, i)
            vent(il, i, j) = (1.-Sigij(il,i,j))*vnk(il) + Sigij(il, i, j)*v(il, i)

! --     Compute liquid water static energy of mixed draughts
!    IF (j .GT. i) THEN
!      awat=elij(il,i,j)-(1.-ep(il,j))*clw(il,j)
!      awat=amax1(awat,0.0)
!    ELSE
!      awat = 0.
!    ENDIF
!    Hent(il,i,j) = (1.-Sigij(il,i,j))*HP(il,i)
!    :         + Sigij(il,i,j)*H(il,i)
!    :         + (LV(il,j)+(cpd-cpv)*t(il,j))*awat
!IM 301008 beg
            hent(il, i, j) = (1.-Sigij(il,i,j))*hp(il, i) + Sigij(il, i, j)*h(il, i)

!jyg<
!            elij(il, i, j) = Qent(il, i, j) - rs(il, j)
!            elij(il, i, j) = elij(il, i, j) + &
!                             ((h(il,j)-hent(il,i,j))*rs(il,j)*lv(il,j) / &
!                              ((cpd*(1.-Qent(il,i,j))+Qent(il,i,j)*cpv)*rrv*t(il,j)*t(il,j)))
!            elij(il, i, j) = elij(il, i, j) / &
!                             (1.+lv(il,j)*lv(il,j)*rs(il,j) / &
!                              ((cpd*(1.-Qent(il,i,j))+Qent(il,i,j)*cpv)*rrv*t(il,j)*t(il,j)))
!
!       Computation of condensate amount Elij, taking into account the ice fraction frac
!       Warning : the same saturation humidity rs is used over both liquid water and ice; this
!                 should be corrected.
!
!  Heat capacity of mixed draught
    cpm = cpd+Qent(il,i,j)*(cpv-cpd)
!
    IF (cvflag_ice .and. frac(il,j) .gt. 0.) THEN
            elij(il, i, j) = Qent(il, i, j) - rs(il, j)
            elij(il, i, j) = elij(il, i, j) + &
                             (h(il,j)-hent(il,i,j)+(cpv-cpd)*(Qent(il,i,j)-rr(il,j))*t(il,j))* &
                             rs(il,j)*lv(il,j) / (cpm*rrv*t(il,j)*t(il,j))
            elij(il, i, j) = elij(il, i, j) / &
                             (1.+(lv(il,j)+frac(il,j)*lf(il,j))*lv(il,j)*rs(il,j) / &
                              (cpm*rrv*t(il,j)*t(il,j)))
    ELSE
            elij(il, i, j) = Qent(il, i, j) - rs(il, j)
            elij(il, i, j) = elij(il, i, j) + &
                             (h(il,j)-hent(il,i,j)+(cpv-cpd)*(Qent(il,i,j)-rr(il,j))*t(il,j))* &
                             rs(il,j)*lv(il,j) / (cpm*rrv*t(il,j)*t(il,j))
            elij(il, i, j) = elij(il, i, j) / &
                             (1.+lv(il,j)*lv(il,j)*rs(il,j) / &
                              (cpm*rrv*t(il,j)*t(il,j)))
    ENDIF
!>jyg
            elij(il, i, j) = max(elij(il,i,j), 0.)

            elij(il, i, j) = min(elij(il,i,j), Qent(il,i,j))

            IF (j>i) THEN
              awat = elij(il, i, j) - (1.-ep(il,j))*clw(il, j)
              awat = amax1(awat, 0.0)
            ELSE
              awat = 0.
            END IF

! print *,h(il,j)-hent(il,i,j),LV(il,j)*rs(il,j)/(cpd*rrv*t(il,j)*
! :         t(il,j))

!jyg<
!            hent(il, i, j) = hent(il, i, j) + (lv(il,j)+(cpd-cpv)*t(il,j))*awat
! Mixed draught temperature at level j
    IF (cvflag_ice .and. frac(il,j) .gt. 0.) THEN
          Tm = t(il,j) + (Qent(il,i,j)-elij(il,i,j)-rs(il,j))*rrv*t(il,j)*t(il,j)/(lv(il,j)*rs(il,j))
          hent(il, i, j) = hent(il, i, j) + (lv(il,j)+frac(il,j)*lf(il,j)+(cpd-cpv)*Tm)*awat
    ELSE
          Tm = t(il,j) + (Qent(il,i,j)-elij(il,i,j)-rs(il,j))*rrv*t(il,j)*t(il,j)/(lv(il,j)*rs(il,j))
          hent(il, i, j) = hent(il, i, j) + (lv(il,j)+(cpd-cpv)*Tm)*awat
    ENDIF
!>jyg

!IM 301008 end

! print *,'mix : i,j,hent(il,i,j),Sigij(il,i,j) ',
! :               i,j,hent(il,i,j),Sigij(il,i,j)

! --      ASij is the integral of P(F) over the relevant F interval
            ASij(il) = ASij(il) + abs(Qmixmax(il)*(1.-Sjmax(il))+Rmixmax(il) - &
                                      Qmixmin(il)*(1.-Sjmin(il))-Rmixmin(il))

          END IF
        END IF
      END DO
!jyg!      DO k = 1, ntra
!jyg!        DO il = 1, ncum
!jyg!          IF ((i>=icb(il)) .AND. (i<=inb(il)) .AND. &
!jyg!              (j>=(icb(il)-1)) .AND. (j<=inb(il)) .AND. &
!jyg!              lwork(il)) THEN
!jyg!            IF (Sij(il,i,j)>0.0) THEN
!jyg!              traent(il, i, j, k) = Sigij(il, i, j)*tra(il, i, k) + &
!jyg!                                    (1.-Sigij(il,i,j))*tra(il, nk(il), k)
!jyg!            END IF
!jyg!          END IF
!jyg!        END DO
!jyg!      END DO

! --    If I=J (detrainement and entrainement at the same level), then only the
! --    adiabatic ascent part of the mixture is considered
      IF (i==j) THEN
        DO il = 1, ncum
          IF (i>=icb(il) .AND. i<=inb(il) .AND. &
              j>=(icb(il)-1) .AND. j<=inb(il) .AND. &
              lwork(il)) THEN
            IF (Sij(il,i,j)>0.0) THEN
!!              rti = qnk(il) - ep(il, i)*clw(il, i)
              rti = qta(il,i-1) - ep(il, i)*clw(il, i)
!!!             Ment(il,i,i) = m(il,i)*abs(Qmixmax(il)*(1.-Sjmax(il))
              Ment(il, i, i) = abs(Qmixmax(il)*(1.-Sjmax(il))+Rmixmax(il) - &
                                   Qmixmin(il)*(1.-Sjmin(il))-Rmixmin(il))
              Qent(il, i, i) = rti
              uent(il, i, i) = unk(il)
              vent(il, i, i) = vnk(il)
              hent(il, i, i) = hp(il, i)
              elij(il, i, i) = clw(il, i)*(1.-ep(il,i))
              Sigij(il, i, i) = 0.
            END IF
          END IF
        END DO
!jyg!        DO k = 1, ntra
!jyg!          DO il = 1, ncum
!jyg!            IF ((i>=icb(il)) .AND. (i<=inb(il)) .AND. &
!jyg!                (j>=(icb(il)-1)) .AND. (j<=inb(il)) .AND. &
!jyg!                lwork(il)) THEN
!jyg!              IF (Sij(il,i,j)>0.0) THEN
!jyg!                traent(il, i, i, k) = tra(il, nk(il), k)
!jyg!              END IF
!jyg!            END IF
!jyg!          END DO
!jyg!        END DO

      END IF

! ---------------------------------------------------------------
175 END DO        ! End loop on destination level "j"
! ---------------------------------------------------------------

    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il)) THEN
        ASij(il) = amax1(1.0E-16, ASij(il))
!jyg+lluis<
!!        ASij(il) = 1.0/ASij(il)
        ASij_inv(il) = 1.0/ASij(il)
!   IF the F-interval spanned by possible mixtures is less than 0.01, no mixing occurs
        IF (ASij_inv(il) > 100.)  ASij_inv(il) = 0.
!>jyg+lluis
        csum(il, i) = 0.0
      END IF
    END DO

    DO j = minorig, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. &
            j>=(icb(il)-1) .AND. j<=inb(il)) THEN
!jyg          Ment(il, i, j) = Ment(il, i, j)*ASij(il)
          Ment(il, i, j) = Ment(il, i, j)*ASij_inv(il)
        END IF
      END DO
    END DO

    DO j = minorig, nl
      DO il = 1, ncum
        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. &
            j>=(icb(il)-1) .AND. j<=inb(il)) THEN
          csum(il, i) = csum(il, i) + Ment(il, i, j)
        END IF
      END DO
    END DO

    DO il = 1, ncum
      IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. csum(il,i)<1.) THEN
! cc     :     .and. csum(il,i).lt.m(il,i) ) then
        nent(il, i) = 0
! cc        Ment(il,i,i)=m(il,i)
        Ment(il, i, i) = 1.
!!        Qent(il, i, i) = qnk(il) - ep(il, i)*clw(il, i)
        Qent(il, i, i) = qta(il,i-1) - ep(il, i)*clw(il, i)
        uent(il, i, i) = unk(il)
        vent(il, i, i) = vnk(il)
        elij(il, i, i) = clw(il, i)*(1.-ep(il,i))
        IF (fl_cor_ebil .GE. 2) THEN
          hent(il, i, i) = hp(il,i)
          Sigij(il, i, i) = 0.0
        ELSE
          Sij(il, i, i) = 0.0
        ENDIF
      END IF
    END DO ! il

!jyg!    DO j = 1, ntra
!jyg!      DO il = 1, ncum
!jyg!        IF (i>=icb(il) .AND. i<=inb(il) .AND. lwork(il) .AND. csum(il,i)<1.) THEN
!jyg!! cc     :     .and. csum(il,i).lt.m(il,i) ) then
!jyg!          traent(il, i, i, j) = tra(il, nk(il), j)
!jyg!        END IF
!jyg!      END DO
!jyg!    END DO

! ---------------------------------------------------------------
789 END DO              ! End loop on origin level "i"
! ---------------------------------------------------------------


  RETURN
END SUBROUTINE cv3p_mixing

