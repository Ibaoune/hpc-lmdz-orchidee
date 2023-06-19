!
! $Header$
!
!------------------------------------------------------------
! Parameters for convectL, iflag_con=30:
! (includes - microphysical parameters, 
!			- parameters that control the rate of approach 
!               to quasi-equilibrium)
!			- noff & minorig (previously in input of convect1)
!------------------------------------------------------------

      integer noff, minorig, nl, nlp, nlm
      real sigd, spfac
!IM cf. FH : pour compatibilite avec conema3 TEMPORAIRE   real pbcrit, ptcrit, epmax
      real pbcrit, ptcrit
      real omtrain
      real dtovsh, dpbase, dttrig
      real dtcrit, tau, beta, alpha
      real delta
      real betad

      COMMON /cv30param/  noff, minorig, nl, nlp, nlm &
                      ,  sigd, spfac &
!IM cf. FH : pour compatibilite avec conema3 TEMPORAIRE  :                ,pbcrit, ptcrit, epmax
                      ,pbcrit, ptcrit &
                      ,omtrain &
                      ,dtovsh, dpbase, dttrig &
                      ,dtcrit, tau, beta, alpha, delta, betad

!$OMP THREADPRIVATE(/cv30param/)
