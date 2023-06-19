
! $Header$


! ================================================================================

SUBROUTINE clouds_gno(klon, nd, r, rs, qsub, ptconv, ratqsc, cldf)
  IMPLICIT NONE

  ! --------------------------------------------------------------------------------

  ! Inputs:

  ! ND----------: Number of vertical levels
  ! R--------ND-: Domain-averaged mixing ratio of total water
  ! RS-------ND-: Mean saturation humidity mixing ratio within the gridbox
  ! QSUB-----ND-: Mixing ratio of condensed water within clouds associated
  ! with SUBGRID-SCALE condensation processes (here, it is
  ! predicted by the convection scheme)
  ! Outputs:

  ! PTCONV-----ND-: Point convectif = TRUE
  ! RATQSC-----ND-: Largeur normalisee de la distribution
  ! CLDF-----ND-: Fraction nuageuse

  ! --------------------------------------------------------------------------------


  INTEGER klon, nd
  REAL r(klon, nd), rs(klon, nd), qsub(klon, nd)
  LOGICAL ptconv(klon, nd)
  REAL ratqsc(klon, nd)
  REAL cldf(klon, nd)

  ! -- parameters controlling the iteration:
  ! --    nmax    : maximum nb of iterations (hopefully never reached)
  ! --    epsilon : accuracy of the numerical resolution
  ! --    vmax    : v-value above which we use an asymptotic expression for
  ! ERF(v)

  INTEGER nmax
  PARAMETER (nmax=10)
  REAL epsilon, vmax0, vmax(klon)
  PARAMETER (epsilon=0.02, vmax0=2.0)

  REAL min_mu, min_q
  PARAMETER (min_mu=1.E-12, min_q=1.E-12)

  INTEGER i, k, n, m
  REAL mu(klon), qsat, delta(klon), beta(klon)
  REAL zu2, zv2
  REAL xx(klon), aux(klon), coeff, block
  REAL dist, fprime, det
  REAL pi, u, v, erfcu, erfcv
  REAL xx1, xx2
  REAL erf, hsqrtlog_2, v2
  REAL sqrtpi, sqrt2, zx1, zx2, exdel
  ! lconv = true si le calcul a converge (entre autre si qsub < min_q)
  LOGICAL lconv(klon)

  ! cdir arraycomb
  cldf(1:klon, 1:nd) = 0.0 ! cym
  ratqsc(1:klon, 1:nd) = 0.0
  ptconv(1:klon, 1:nd) = .FALSE.
  ! cdir end arraycomb

  pi = acos(-1.)
  sqrtpi = sqrt(pi)
  sqrt2 = sqrt(2.)
  hsqrtlog_2 = 0.5*sqrt(log(2.))

  DO k = 1, nd

    DO i = 1, klon ! vector
      mu(i) = r(i, k)
      mu(i) = max(mu(i), min_mu)
      qsat = rs(i, k)
      qsat = max(qsat, min_mu)
      delta(i) = log(mu(i)/qsat)
      ! enddo ! vector


      ! ***          There is no subgrid-scale condensation;        ***
      ! ***   the scheme becomes equivalent to an "all-or-nothing"  ***
      ! ***             large-scale condensation scheme.            ***



      ! ***     Some condensation is produced at the subgrid-scale       ***
      ! ***                                                              ***
      ! ***       PDF = generalized log-normal distribution (GNO)        ***
      ! ***   (k<0 because a lower bound is considered for the PDF)      ***
      ! ***                                                              ***
      ! ***  -> Determine x (the parameter k of the GNO PDF) such        ***
      ! ***  that the contribution of subgrid-scale processes to         ***
      ! ***  the in-cloud water content is equal to QSUB(K)              ***
      ! ***  (equations (13), (14), (15) + Appendix B of the paper)      ***
      ! ***                                                              ***
      ! ***    Here, an iterative method is used for this purpose        ***
      ! ***    (other numerical methods might be more efficient)         ***
      ! ***                                                              ***
      ! ***          NB: the "error function" is called ERF              ***
      ! ***                 (ERF in double precision)                   ***


      ! On commence par eliminer les cas pour lesquels on n'a pas
      ! suffisamment d'eau nuageuse.

      ! do i=1,klon ! vector

      IF (qsub(i,k)<min_q) THEN
        ptconv(i, k) = .FALSE.
        ratqsc(i, k) = 0.
        lconv(i) = .TRUE.

        ! Rien on a deja initialise

      ELSE

        lconv(i) = .FALSE.
        vmax(i) = vmax0

        beta(i) = qsub(i, k)/mu(i) + exp(-min(0.0,delta(i)))

        ! --  roots of equation v > vmax:

        det = delta(i) + vmax(i)*vmax(i)
        IF (det<=0.0) vmax(i) = vmax0 + 1.0
        det = delta(i) + vmax(i)*vmax(i)

        IF (det<=0.) THEN
          xx(i) = -0.0001
        ELSE
          zx1 = -sqrt2*vmax(i)
          zx2 = sqrt(1.0+delta(i)/(vmax(i)*vmax(i)))
          xx1 = zx1*(1.0-zx2)
          xx2 = zx1*(1.0+zx2)
          xx(i) = 1.01*xx1
          IF (xx1>=0.0) xx(i) = 0.5*xx2
        END IF
        IF (delta(i)<0.) xx(i) = -hsqrtlog_2

      END IF

    END DO ! vector

    ! ----------------------------------------------------------------------
    ! Debut des nmax iterations pour trouver la solution.
    ! ----------------------------------------------------------------------

    DO n = 1, nmax

      DO i = 1, klon ! vector
        IF (.NOT. lconv(i)) THEN

          u = delta(i)/(xx(i)*sqrt2) + xx(i)/(2.*sqrt2)
          v = delta(i)/(xx(i)*sqrt2) - xx(i)/(2.*sqrt2)
          v2 = v*v

          IF (v>vmax(i)) THEN

            IF (abs(u)>vmax(i) .AND. delta(i)<0.) THEN

              ! -- use asymptotic expression of erf for u and v large:
              ! ( -> analytic solution for xx )
              exdel = beta(i)*exp(delta(i))
              aux(i) = 2.0*delta(i)*(1.-exdel)/(1.+exdel)
              IF (aux(i)<0.) THEN
                ! print*,'AUX(',i,',',k,')<0',aux(i),delta(i),beta(i)
                aux(i) = 0.
              END IF
              xx(i) = -sqrt(aux(i))
              block = exp(-v*v)/v/sqrtpi
              dist = 0.0
              fprime = 1.0

            ELSE

              ! -- erfv -> 1.0, use an asymptotic expression of erfv for v
              ! large:

              erfcu = 1.0 - erf(u)
              ! !!! ATTENTION : rajout d'un seuil pour l'exponentiel
              aux(i) = sqrtpi*erfcu*exp(min(v2,100.))
              coeff = 1.0 - 0.5/(v2) + 0.75/(v2*v2)
              block = coeff*exp(-v2)/v/sqrtpi
              dist = v*aux(i)/coeff - beta(i)
              fprime = 2.0/xx(i)*(v2)*(exp(-delta(i))-u*aux(i)/coeff)/coeff

            END IF ! ABS(u)

          ELSE

            ! -- general case:

            erfcu = 1.0 - erf(u)
            erfcv = 1.0 - erf(v)
            block = erfcv
            dist = erfcu/erfcv - beta(i)
            zu2 = u*u
            zv2 = v2
            IF (zu2>20. .OR. zv2>20.) THEN
              ! print*,'ATTENTION !!! xx(',i,') =', xx(i)
              ! print*,'ATTENTION !!! klon,ND,R,RS,QSUB,PTCONV,RATQSC,CLDF',
              ! .klon,ND,R(i,k),RS(i,k),QSUB(i,k),PTCONV(i,k),RATQSC(i,k),
              ! .CLDF(i,k)
              ! print*,'ATTENTION !!! zu2 zv2 =',zu2(i),zv2(i)
              zu2 = 20.
              zv2 = 20.
              fprime = 0.
            ELSE
              fprime = 2./sqrtpi/xx(i)/(erfcv*erfcv)* &
                (erfcv*v*exp(-zu2)-erfcu*u*exp(-zv2))
            END IF
          END IF ! x

          ! -- test numerical convergence:

          ! if (beta(i).lt.1.e-10) then
          ! print*,'avant test ',i,k,lconv(i),u(i),v(i),beta(i)
          ! stop
          ! endif
          IF (abs(fprime)<1.E-11) THEN
            ! print*,'avant test fprime<.e-11 '
            ! s        ,i,k,lconv(i),u(i),v(i),beta(i),fprime(i)
            ! print*,'klon,ND,R,RS,QSUB',
            ! s        klon,ND,R(i,k),rs(i,k),qsub(i,k)
            fprime = sign(1.E-11, fprime)
          END IF


          IF (abs(dist/beta(i))<epsilon) THEN
            ! print*,'v-u **2',(v(i)-u(i))**2
            ! print*,'exp v-u **2',exp((v(i)-u(i))**2)
            ptconv(i, k) = .TRUE.
            lconv(i) = .TRUE.
            ! borne pour l'exponentielle
            ratqsc(i, k) = min(2.*(v-u)*(v-u), 20.)
            ratqsc(i, k) = sqrt(exp(ratqsc(i,k))-1.)
            cldf(i, k) = 0.5*block
          ELSE
            xx(i) = xx(i) - dist/fprime
          END IF
          ! print*,'apres test ',i,k,lconv(i)

        END IF ! lconv
      END DO ! vector

      ! ----------------------------------------------------------------------
      ! Fin des nmax iterations pour trouver la solution.
    END DO ! n
    ! ----------------------------------------------------------------------


  END DO
  ! K
  RETURN
END SUBROUTINE clouds_gno



