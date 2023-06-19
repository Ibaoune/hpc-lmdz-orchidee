
! $Header$

SUBROUTINE tlift(p, t, rr, rs, gz, plcl, icb, nk, tvp, tpk, clw, nd, nl, &
    dtvpdt1, dtvpdq1)
  IMPLICIT NONE
  ! Argument NK ajoute (jyg) = Niveau de depart de la
  ! convection
  INTEGER icb, nk, nd, nl
  INTEGER,PARAMETER :: na=60
  REAL gz(nd), tpk(nd), clw(nd), plcl
  REAL t(nd), rr(nd), rs(nd), tvp(nd), p(nd)
  REAL dtvpdt1(nd), dtvpdq1(nd) ! Derivatives of parcel virtual
  ! temperature wrt T1 and Q1

  REAL clw_new(na), qi(na)
  REAL dtpdt1(na), dtpdq1(na) ! Derivatives of parcel temperature
  ! wrt T1 and Q1
  REAL gravity, cpd, cpv, cl, ci, cpvmcl, clmci, eps, alv0, alf0
  REAL cpp, cpinv, ah0, alf, tg, s, ahg, tc, denom, alv, es, esi
  REAL qsat_new, snew
  INTEGER icbl, i, imin, j, icb1

  LOGICAL ice_conv

  ! ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***

  ! sb        CPD=1005.7
  ! sb      CPV=1870.0
  ! sb        CL=4190.0
  ! sb        CPVMCL=2320.0
  ! sb        RV=461.5
  ! sb        RD=287.04
  ! sb        EPS=RD/RV
  ! sb        ALV0=2.501E6
  ! cccccccccccccccccccccc
  ! constantes coherentes avec le modele du Centre Europeen
  ! sb      RD = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 28.9644
  ! sb      RV = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 18.0153
  ! sb      CPD = 3.5 * RD
  ! sb      CPV = 4.0 * RV
  ! sb      CL = 4218.0
  ! sb      CI=2090.0
  ! sb      CPVMCL=CL-CPV
  ! sb      CLMCI=CL-CI
  ! sb      EPS=RD/RV
  ! sb      ALV0=2.5008E+06
  ! sb      ALF0=3.34E+05

  ! ccccccccccc
  ! on utilise les constantes thermo du Centre Europeen: (SB)

  include "YOMCST.h"
  gravity = rg !sb: Pr que gravite ne devienne pas humidite!

  cpd = rcpd
  cpv = rcpv
  cl = rcw
  ci = rcs
  cpvmcl = cl - cpv
  clmci = cl - ci
  eps = rd/rv
  alv0 = rlvtt
  alf0 = rlmlt ! (ALF0 = RLSTT-RLVTT)

  ! ccccccccccccccccccccc

  ! ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***

  icb1 = max(icb, 2)
  icb1 = min(icb, nl)

  ! jyg1
  ! C      CPP=CPD*(1.-RR(1))+RR(1)*CPV
  cpp = cpd*(1.-rr(nk)) + rr(nk)*cpv
  ! jyg2
  cpinv = 1./cpp
  ! jyg1
  ! ICB may be below condensation level
  ! CC        DO 100 I=1,ICB1-1
  ! CC         TPK(I)=T(1)-GZ(I)*CPINV
  ! CC         TVP(I)=TPK(I)*(1.+RR(1)/EPS)
  DO i = 1, icb1
    clw(i) = 0.0
  END DO

  DO i = nk, icb1
    tpk(i) = t(nk) - (gz(i)-gz(nk))*cpinv
    ! jyg1
    ! CC         TVP(I)=TPK(I)*(1.+RR(NK)/EPS)
    tvp(i) = tpk(i)*(1.+rr(nk)/eps-rr(nk))
    ! jyg2
    dtvpdt1(i) = 1. + rr(nk)/eps - rr(nk)
    dtvpdq1(i) = tpk(i)*(1./eps-1.)

    ! jyg2

  END DO


  ! ***  FIND LIFTED PARCEL TEMPERATURE AND MIXING RATIO    ***

  ! jyg1
  ! C        AH0=(CPD*(1.-RR(1))+CL*RR(1))*T(1)
  ! C     $     +RR(1)*(ALV0-CPVMCL*(T(1)-273.15))
  ah0 = (cpd*(1.-rr(nk))+cl*rr(nk))*t(nk) + rr(nk)*(alv0-cpvmcl*(t(nk)-273.15 &
    )) + gz(nk)
  ! jyg2

  ! jyg1
  imin = icb1
  ! If ICB is below LCL, start loop at ICB+1
  IF (plcl<p(icb1)) imin = min(imin+1, nl)

  ! CC        DO 300 I=ICB1,NL
  DO i = imin, nl
    ! jyg2
    alv = alv0 - cpvmcl*(t(i)-273.15)
    alf = alf0 + clmci*(t(i)-273.15)

    rg = rs(i)
    tg = t(i)
    ! S=CPD+ALV*ALV*RG/(RV*T(I)*T(I))
    ! jyg1
    ! C        S=CPD*(1.-RR(1))+CL*RR(1)+ALV*ALV*RG/(RV*T(I)*T(I))
    s = cpd*(1.-rr(nk)) + cl*rr(nk) + alv*alv*rg/(rv*t(i)*t(i))
    ! jyg2
    s = 1./s

    DO j = 1, 2
      ! jyg1
      ! C         AHG=CPD*TG+(CL-CPD)*RR(1)*TG+ALV*RG+GZ(I)
      ahg = cpd*tg + (cl-cpd)*rr(nk)*tg + alv*rg + gz(i)
      ! jyg2
      tg = tg + s*(ah0-ahg)
      tc = tg - 273.15
      denom = 243.5 + tc
      denom = max(denom, 1.0)

      ! FORMULE DE BOLTON POUR PSAT

      es = 6.112*exp(17.67*tc/denom)
      rg = eps*es/(p(i)-es*(1.-eps))


    END DO

    ! jyg1
    ! C        TPK(I)=(AH0-GZ(I)-ALV*RG)/(CPD+(CL-CPD)*RR(1))
    tpk(i) = (ah0-gz(i)-alv*rg)/(cpd+(cl-cpd)*rr(nk))
    ! jyg2
    ! TPK(I)=(AH0-GZ(I)-ALV*RG-(CL-CPD)*T(I)*RR(1))/CPD

    ! jyg1
    ! C        CLW(I)=RR(1)-RG
    clw(i) = rr(nk) - rg
    ! jyg2
    clw(i) = max(0.0, clw(i))
    ! jyg1
    ! CC        TVP(I)=TPK(I)*(1.+RG/EPS)
    tvp(i) = tpk(i)*(1.+rg/eps-rr(nk))
    ! jyg2

    ! jyg1       Derivatives

    dtpdt1(i) = cpd*s
    dtpdq1(i) = alv*s

    dtvpdt1(i) = dtpdt1(i)*(1.+rg/eps-rr(nk)+alv*rg/(rd*tpk(i)))
    dtvpdq1(i) = dtpdq1(i)*(1.+rg/eps-rr(nk)+alv*rg/(rd*tpk(i))) - tpk(i)

    ! jyg2

  END DO

  ice_conv = .FALSE.

  IF (ice_conv) THEN

    ! JAM
    ! RAJOUT DE LA PROCEDURE ICEFRAC

    ! sb        CALL ICEFRAC(T,CLW,CLW_NEW,QI,ND,NL)

    DO i = icb1, nl
      IF (t(i)<263.15) THEN
        tg = tpk(i)
        tc = tpk(i) - 273.15
        denom = 243.5 + tc
        es = 6.112*exp(17.67*tc/denom)
        alv = alv0 - cpvmcl*(t(i)-273.15)
        alf = alf0 + clmci*(t(i)-273.15)

        DO j = 1, 4
          esi = exp(23.33086-(6111.72784/tpk(i))+0.15215*log(tpk(i)))
          qsat_new = eps*esi/(p(i)-esi*(1.-eps))
          ! CC        SNEW=
          ! CPD*(1.-RR(1))+CL*RR(1)+ALV*ALV*QSAT_NEW/(RV*TPK(I)*TPK(I))
          snew = cpd*(1.-rr(nk)) + cl*rr(nk) + alv*alv*qsat_new/(rv*tpk(i)* &
            tpk(i))

          snew = 1./snew
          tpk(i) = tg + (alf*qi(i)+alv*rg*(1.-(esi/es)))*snew
          ! @$$        PRINT*,'################################'
          ! @$$        PRINT*,TPK(I)
          ! @$$        PRINT*,(ALF*QI(I)+ALV*RG*(1.-(ESI/ES)))*SNEW
        END DO
        ! CC        CLW(I)=RR(1)-QSAT_NEW
        clw(i) = rr(nk) - qsat_new
        clw(i) = max(0.0, clw(i))
        ! jyg1
        ! CC        TVP(I)=TPK(I)*(1.+QSAT_NEW/EPS)
        tvp(i) = tpk(i)*(1.+qsat_new/eps-rr(nk))
        ! jyg2
      ELSE
        CONTINUE
      END IF

    END DO

  END IF


  ! *****************************************************
  ! * BK :  RAJOUT DE LA TEMPERATURE DES ASCENDANCES
  ! *   NON DILUES AU  NIVEAU KLEV = ND
  ! *   POSONS LE ENVIRON EGAL A CELUI DE KLEV-1
  ! *******************************************************

  tpk(nl+1) = tpk(nl)

  ! ******************************************************

  rg = gravity ! RG redevient la gravite de YOMCST (sb)


  RETURN
END SUBROUTINE tlift








