SUBROUTINE cv3_inicp()

  ! **************************************************************
  ! *
  ! CV3_INIP Lecture des choix de lois de probabilité de mélange*
  ! et calcul de leurs coefficients normalisés.        *
  ! *
  ! written by   : Jean-Yves Grandpeix, 06/06/2006, 19.39.27    *
  ! modified by :                                               *
  ! **************************************************************
  IMPLICIT NONE
  include "YOMCST2.h"

  INTEGER iflag_clos
  CHARACTER (LEN=20) :: modname = 'cv3_inicp'
  CHARACTER (LEN=80) :: abort_message

  ! --   Mixing probability distribution functions

  REAL qcoef1, qcoef2, qff, qfff, qmix, rmix, qmix1, rmix1, qmix2, rmix2, f
  REAL :: sumcoef,sigma,aire,pdf,mu,df,ff

  qcoef1(f) = tanh(f/gammas)
  qcoef2(f) = (tanh(f/gammas)+gammas*log(cosh((1.-f)/gammas)/cosh(f/gammas)))
  qff(f) = max(min(f,1.), 0.)
  qfff(f) = min(qff(f), scut)
  qmix1(f) = (tanh((qff(f)-fmax)/gammas)+qcoef1max)/qcoef2max
  rmix1(f) = (gammas*log(cosh((qff(f)-fmax)/gammas))+qff(f)*qcoef1max)/ &
    qcoef2max
  qmix2(f) = -log(1.-qfff(f))/scut
  rmix2(f) = (qfff(f)+(1.-qff(f))*log(1.-qfff(f)))/scut
  qmix(f) = qqa1*qmix1(f) + qqa2*qmix2(f)
  rmix(f) = qqa1*rmix1(f) + qqa2*rmix2(f)

  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  ! ===========================================================================
  ! READ IN PARAMETERS FOR THE MIXING DISTRIBUTION
  ! AND PASS THESE THROUGH A COMMON BLOCK TO SUBROUTINE CONVECT etc.
  ! (Written by V.T.J. Phillips, 20-30/Jan/99)
  ! ===========================================================================

  ! line 1:  a flag (0 or 1) to decide whether P(F) = 1 or the general P(F)
  ! is to be
  ! used, followed by SCUT, which is the cut-off value of F in CONVECT
  ! line 2:  blank
  ! line 3:  the coefficients for the linear combination of P(F)s to
  ! make the general P(F)
  ! line 4:  blank
  ! line 5:  gammas, Fmax for the cosh^2 component of P(F)
  ! line 6:  blank
  ! line 7:  alphas for the 1st irrational P(F)
  ! line 8:  blank
  ! line 9:  betas  for the 2nd irrational P(F)


  ! open(57,file='parameter_mix.data')

  ! read(57,*) iflag_clos
  ! read(57,*) iflag_mix, scut
  ! read(57,*)
  ! if(iflag_mix .gt. 0) then
  ! read(57,*) qqa1, qqa2
  ! read(57,*)
  ! read(57,*) gammas, Fmax
  ! read(57,*)
  ! read(57,*) alphas
  ! endif
  ! close(57)

  IF (iflag_mix>0) THEN

    ! --      Normalize Pdf weights

    sumcoef = qqa1 + qqa2
    qqa1 = qqa1/sumcoef
    qqa2 = qqa2/sumcoef

    qcoef1max = qcoef1(fmax)
    qcoef2max = qcoef2(fmax)

    sigma = 0.
    aire = 0.0
    pdf = 0.0
    mu = 0.0
    df = 0.0001

    ! do ff = 0.0 + df, 1.0 - 2.*df, df
    ff = df
    DO WHILE (ff<=1.0-2.*df)
      pdf = (qmix(ff+df)-qmix(ff))*(1.-ff)/df
      aire = aire + (qmix(ff+df)-qmix(ff))*(1.-ff)
      mu = mu + pdf*ff*df
      ! c              write(*,*) pdf,  Qmix(ff), aire, ff
      ff = ff + df
    END DO

    ! do ff=0.0+df,1.0 - 2.*df,df
    ff = df
    DO WHILE (ff<=1.0-2.*df)
      pdf = (qmix(ff+df)-qmix(ff))*(1.-ff)/df
      sigma = sigma + pdf*(ff-mu)*(ff-mu)*df
      ff = ff + df
    END DO
    sigma = sqrt(sigma)

    IF (abs(aire-1.0)>0.02) THEN
      PRINT *, 'WARNING:: AREA OF MIXING PDF IS::', aire
      abort_message = ''
      CALL abort_physic(modname, abort_message, 1)
    ELSE
      PRINT *, 'Area, mean & std deviation are ::', aire, mu, sigma
    END IF
  END IF !  (iflag_mix .gt. 0)

  RETURN
END SUBROUTINE cv3_inicp
