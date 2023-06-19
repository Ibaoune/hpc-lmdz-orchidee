!
! $Header$
!
!------------------------------------------------------------
! Parameters for convectL:
! (includes - microphysical parameters, 
!			- parameters that control the rate of approach 
!               to quasi-equilibrium)
!			- noff & minorig (previously in input of convect1)
!------------------------------------------------------------

      integer noff, minorig, nl, nlp, nlm
      real elcrit, tlcrit
      real entp
      real sigs, sigd
      real omtrain, omtsnow, coeffr, coeffs
      real dtmax
      real cu
      real betad
      real alpha, damp
      real delta

      COMMON /cvparam/ noff, minorig, nl, nlp, nlm &
                      ,elcrit, tlcrit &
                      ,entp, sigs, sigd &
                      ,omtrain, omtsnow, coeffr, coeffs &
                      ,dtmax, cu, betad, alpha, damp, delta

!$OMP THREADPRIVATE(/cvparam/)
