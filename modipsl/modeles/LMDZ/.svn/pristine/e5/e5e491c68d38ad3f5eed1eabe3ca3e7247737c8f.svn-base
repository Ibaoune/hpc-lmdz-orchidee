
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/plevel.F,v 1.1.1.1.10.1 2006/08/17
! 15:41:51 fairhead Exp $

! ================================================================
! ================================================================
SUBROUTINE plevel_new(ilon, ilev, klevstd, lnew, pgcm, pres, qgcm, qpres)
  ! ================================================================
  ! ================================================================
  USE netcdf
  USE dimphy
#ifdef CPP_IOIPSL 
  USE phys_state_var_mod, ONLY: missing_val_nf90
#endif
#ifdef CPP_XIOS
  USE wxios, ONLY: missing_val
#endif

  IMPLICIT NONE

  ! ================================================================

  ! Interpoler des champs 3-D u, v et g du modele a un niveau de
  ! pression donnee (pres)

  ! INPUT:  ilon ----- nombre de points
  ! ilev ----- nombre de couches
  ! lnew ----- true si on doit reinitialiser les poids
  ! pgcm ----- pressions modeles
  ! pres ----- pression vers laquelle on interpolle
  ! Qgcm ----- champ GCM
  ! Qpres ---- champ interpolle au niveau pres

  ! ================================================================

  ! arguments :
  ! -----------

  INTEGER ilon, ilev, klevstd
  LOGICAL lnew

  REAL pgcm(ilon, ilev)
  REAL qgcm(ilon, ilev)
  REAL pres(klevstd)
  REAL qpres(ilon, klevstd)

  ! local :
  ! -------

  ! ym      INTEGER lt(klon), lb(klon)
  ! ym      REAL ptop, pbot, aist(klon), aisb(klon)

  ! ym      save lt,lb,ptop,pbot,aist,aisb
  INTEGER, ALLOCATABLE, SAVE, DIMENSION (:, :) :: lt, lb
  REAL, ALLOCATABLE, SAVE, DIMENSION (:, :) :: aist, aisb
  !$OMP THREADPRIVATE(lt,lb,aist,aisb)
  REAL, SAVE :: ptop, pbot
  !$OMP THREADPRIVATE(ptop, pbot)
  LOGICAL, SAVE :: first = .TRUE.
  INTEGER :: nlev
  !$OMP THREADPRIVATE(first)
  INTEGER i, k

! REAL missing_val
#ifndef CPP_XIOS
  REAL :: missing_val
#endif

! missing_val = nf90_fill_real

#ifndef CPP_XIOS
      missing_val=missing_val_nf90
#endif

  IF (first) THEN
    ALLOCATE (lt(klon,klevstd), lb(klon,klevstd))
    ALLOCATE (aist(klon,klevstd), aisb(klon,klevstd))
    first = .FALSE.
  END IF

  ! =====================================================================
  IF (lnew) THEN
    ! on reinitialise les reindicages et les poids
    ! =====================================================================


    ! Chercher les 2 couches les plus proches du niveau a obtenir

    ! Eventuellement, faire l'extrapolation a partir des deux couches
    ! les plus basses ou les deux couches les plus hautes:


    DO nlev = 1, klevstd
      DO i = 1, klon
        IF (abs(pres(nlev)-pgcm(i,ilev))<abs(pres(nlev)-pgcm(i,1))) THEN
          lt(i, nlev) = ilev ! 2
          lb(i, nlev) = ilev - 1 ! 1
        ELSE
          lt(i, nlev) = 2
          lb(i, nlev) = 1
        END IF
      END DO
      DO k = 1, ilev - 1
        DO i = 1, klon
          pbot = pgcm(i, k)
          ptop = pgcm(i, k+1)
          IF (ptop<=pres(nlev) .AND. pbot>=pres(nlev)) THEN
            lt(i, nlev) = k + 1
            lb(i, nlev) = k
          END IF
        END DO
      END DO

      ! Interpolation lineaire:
      DO i = 1, klon
        ! interpolation en logarithme de pression:

        ! ...   Modif . P. Le Van    ( 20/01/98) ....
        ! Modif Frederic Hourdin (3/01/02)

        aist(i, nlev) = log(pgcm(i,lb(i,nlev))/pres(nlev))/log(pgcm(i,lb(i, &
          nlev))/pgcm(i,lt(i,nlev)))
        aisb(i, nlev) = log(pres(nlev)/pgcm(i,lt(i,nlev)))/log(pgcm(i,lb(i, &
          nlev))/pgcm(i,lt(i,nlev)))
      END DO
    END DO

  END IF ! lnew

  ! ======================================================================
  ! inteprollation
  ! ET je mets les vents a zero quand je rencontre une montagne
  ! ======================================================================

  DO nlev = 1, klevstd
    DO i = 1, klon
      IF (pgcm(i,1)<pres(nlev)) THEN
        qpres(i, nlev) = missing_val
      ELSE
        qpres(i, nlev) = qgcm(i, lb(i,nlev))*aisb(i, nlev) + &
          qgcm(i, lt(i,nlev))*aist(i, nlev)
      END IF
    END DO
  END DO


  RETURN
END SUBROUTINE plevel_new
