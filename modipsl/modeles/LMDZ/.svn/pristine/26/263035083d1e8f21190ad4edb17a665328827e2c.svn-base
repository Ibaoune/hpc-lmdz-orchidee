!------------------------------------------------------------
! Parameters for convectL, iflag_con=3:
! (includes - microphysical parameters,
!			- parameters that control the rate of approach
!               to quasi-equilibrium)
!			- noff & minorig (previously in input of convect1)
!------------------------------------------------------------

      integer flag_epKEorig
      real flag_wb
      integer cv_flag_feed
      integer noff, minorig, nl, nlp, nlm
      real sigdz, spfac
      real pbcrit, ptcrit
      real elcrit, tlcrit
      real coef_peel
      real omtrain
      real dtovsh, dpbase, dttrig
      real dtcrit, tau, beta, alpha, alpha1
      real T_top_max
      real tau_stop, noconv_stop
      real wbmax
      real delta
      real betad
      real ejectliq
      real ejectice

      COMMON /cv3param/ sigdz, spfac &
                      ,pbcrit, ptcrit &
                      ,elcrit, tlcrit &
                      ,coef_peel &
                      ,omtrain &
                      ,dtovsh, dpbase, dttrig &
                      ,dtcrit, tau, beta, alpha, alpha1 &
                      ,T_top_max &
                      ,tau_stop, noconv_stop &
                      ,wbmax &
                      ,delta, betad  &
                      ,ejectliq, ejectice &
                      ,flag_epKEorig &
                      ,flag_wb, cv_flag_feed &
                      ,noff, minorig, nl, nlp, nlm
!$OMP THREADPRIVATE(/cv3param/)

