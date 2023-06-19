SUBROUTINE cv3_inip()
  ! *******************************************************************
  ! *                                                                  *
  ! CV3_INIP Input = choice of mixing probability laws                 *
  !          Output = normalized coefficients of the probability laws. *
  ! *                                                                  *
  ! written by   : Jean-Yves Grandpeix, 06/06/2006, 19.39.27           *
  ! modified by :                                                      *
  ! *******************************************************************
!
!----------------------------------------------
!  INPUT (from Common YOMCST2 in "YOMCST2.h") : 
! iflag_mix
! gammas    
! alphas
! betas
! Fmax
! scut
!
!----------------------------------------------
!  INPUT/OUTPUT (from and to Common YOMCST2 in "YOMCST2.h") : 
! qqa1
! qqa2
!
!----------------------------------------------
!  OUTPUT (to Common YOMCST2 in "YOMCST2.h") :
! Qcoef1max
! Qcoef2max
!
!----------------------------------------------

  USE print_control_mod, ONLY: prt_level, lunout
  IMPLICIT NONE

  include "YOMCST2.h"

!----------------------------------------------
! Local variables :
   CHARACTER (LEN=20) :: modname = 'cv3_inip'
   CHARACTER (LEN=80) :: abort_message

   REAL               :: sumcoef
   REAL               :: sigma, aire, pdf, mu, df
   REAL               :: ff


  ! --   Mixing probability distribution functions

  REAL qcoef1, qcoef2, qff, qfff, qmix, rmix, qmix1, rmix1, qmix2, rmix2, f

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


  ! c$$$        open(57,file='parameter_mix.data')
  ! c$$$
  ! c$$$        read(57,*) iflag_mix, scut
  ! c$$$        read(57,*)
  ! c$$$        if(iflag_mix .gt. 0) then
  ! c$$$	      read(57,*) qqa1, qqa2
  ! c$$$              read(57,*)
  ! c$$$              read(57,*) gammas, Fmax
  ! c$$$              read(57,*)
  ! c$$$              read(57,*) alphas
  ! c$$$         endif
  ! c$$$	 close(57)


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
      IF (prt_level>9) WRITE (lunout, *) pdf, qmix(ff), aire, ff
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
      WRITE (lunout, *) 'WARNING:: AREA OF MIXING PDF IS::', aire
      abort_message = ''
      CALL abort_physic(modname, abort_message, 1)
    ELSE
      PRINT *, 'Area, mean & std deviation are ::', aire, mu, sigma
    END IF
  END IF !  (iflag_mix .gt. 0)

  RETURN
END SUBROUTINE cv3_inip
