!
! $Header$
!
!  ATTENTION : ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!     -----------------------------------------------------------------
!*    *COMMON* *YOEGWD* - PARAMETERS FOR GRAVITY WAVE DRAG CALCULATIONS
!     -----------------------------------------------------------------
!
      integer NKTOPG,NSTRA
      real GFRCRIT,GKWAKE,GRCRIT,GVCRIT,GKDRAG,GKLIFT
      real GHMAX,GRAHILO,GSIGCR,GSSEC,GTSEC,GVSEC

      REAL GWD_RANDO_RUWMAX
!     Maximum Eliassen-Palm flux at launch level, in "FLOTT_GWD_rando"

      REAL GWD_RANDO_SAT ! saturation parameter in "FLOTT_GWD_rando"
!     S_c in equation (12) of Lott (JGR, vol 118, page 8897, 2013)

      REAL GWD_FRONT_RUWMAX,GWD_FRONT_SAT 
! Same as GWD_RANDO params but for fronal GWs


      COMMON/YOEGWD/ GFRCRIT,GKWAKE,GRCRIT,GVCRIT,GKDRAG,GKLIFT,        &
     &     GHMAX,GRAHILO,GSIGCR,NKTOPG,NSTRA,GSSEC,GTSEC,GVSEC,         &
     &     GWD_RANDO_RUWMAX, gwd_rando_sat,                             &
     &     GWD_FRONT_RUWMAX, gwd_front_sat

      save /YOEGWD/
!$OMP THREADPRIVATE(/YOEGWD/)
