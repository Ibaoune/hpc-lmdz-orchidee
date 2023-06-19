
! $Header$


! ================================================================================

SUBROUTINE clouds_bigauss(klon, nd, r, rs, qtc, sigt, ptconv, ratqsc, cldf)
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
  REAL r(klon, nd), rs(klon, nd), qtc(klon, nd), sigt(klon, nd)
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
  REAL mu, qsat, delta
  REAL sigma1, sigma2, alpha, qconv
  REAL xconv, xenv
  REAL cconv, cenv
  REAL pi, u, v
  REAL erf
  REAL sqrtpi, sqrt2
  ! lconv = true si le calcul a converge (entre autre si qsub < min_q)
  LOGICAL lconv(klon)


  cldf(1:klon, 1:nd) = 0.0 ! cym
  ratqsc(1:klon, 1:nd) = 0.0
  ptconv(1:klon, 1:nd) = .FALSE.
  ! cdir end arraycomb

  pi = acos(-1.)
  sqrtpi = sqrt(pi)
  sqrt2 = sqrt(2.)


  DO k = 1, nd

  DO i = 1, klon ! vector

      mu = r(i, k)
      mu = max(mu, min_mu)
      qsat = rs(i, k)
      qsat = max(qsat, min_mu)
      delta = log(mu/qsat)
      qconv=qtc(i,k)
      alpha=sigt(i,k)

     IF (qconv<min_q) THEN
        ptconv(i, k) = .FALSE.
        ratqsc(i, k) = 0.

        ! Rien on a deja initialise

      ELSE
    
      sigma1=0.1*((qconv-mu)**2)**0.5+0.002*mu
      sigma2=0.1*((qconv-mu)**2)**0.5+0.002*qconv 

!      sigma2=0.09*((qconv-mu)**2)**0.5/(alpha+0.01)**0.5+0.002*qconv 
!-----------------------------------------------------------------------------------------------------------------
! Calcul de la couverture nuageuse et de ratqs
!-----------------------------------------------------------------------------------------------------------------

      xconv=(qsat-qconv)/(sqrt(2.)*sigma2)
      xenv=(qsat-mu)/(sqrt(2.)*sigma1)

      cconv=0.5*(1.-1.*erf(xconv))
      cenv=0.5*(1.-1.*erf(xenv)) 
      cldf(i,k)=alpha*cconv+(1.-1.*alpha)*cenv
      ratqsc(i,k)= alpha*sigma1+(1.-1.*alpha)*sigma2
      ptconv(i,k)= .TRUE.

     END IF

  END DO ! vector


  END DO  ! K

  RETURN
  END SUBROUTINE clouds_bigauss



