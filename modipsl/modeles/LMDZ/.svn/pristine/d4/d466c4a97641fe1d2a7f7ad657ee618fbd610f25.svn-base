
! $Header$

! ================================================================
! ================================================================
SUBROUTINE plevel(ilon, ilev, lnew, pgcm, pres, qgcm, qpres)
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

  INTEGER ilon, ilev
  LOGICAL lnew

  REAL pgcm(ilon, ilev)
  REAL qgcm(ilon, ilev)
  REAL pres
  REAL qpres(ilon)

  ! local :
  ! -------

  ! ym      INTEGER lt(klon), lb(klon)
  ! ym      REAL ptop, pbot, aist(klon), aisb(klon)

  ! ym      save lt,lb,ptop,pbot,aist,aisb
  INTEGER, ALLOCATABLE, SAVE, DIMENSION (:) :: lt, lb
  REAL, ALLOCATABLE, SAVE, DIMENSION (:) :: aist, aisb
  !$OMP THREADPRIVATE(lt,lb,aist,aisb)
  REAL, SAVE :: ptop, pbot
  !$OMP THREADPRIVATE(ptop, pbot)
  LOGICAL, SAVE :: first = .TRUE.
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
    ALLOCATE (lt(klon), lb(klon), aist(klon), aisb(klon))
    first = .FALSE.
  END IF

  ! =====================================================================
  IF (lnew) THEN
    ! on r�nitialise les r�ndicages et les poids
    ! =====================================================================


    ! Chercher les 2 couches les plus proches du niveau a obtenir

    ! Eventuellement, faire l'extrapolation a partir des deux couches
    ! les plus basses ou les deux couches les plus hautes:
    DO i = 1, klon
      IF (abs(pres-pgcm(i,ilev))<abs(pres-pgcm(i,1))) THEN
        lt(i) = ilev ! 2
        lb(i) = ilev - 1 ! 1
      ELSE
        lt(i) = 2
        lb(i) = 1
      END IF
    END DO
    DO k = 1, ilev - 1
      DO i = 1, klon
        pbot = pgcm(i, k)
        ptop = pgcm(i, k+1)
        IF (ptop<=pres .AND. pbot>=pres) THEN
          lt(i) = k + 1
          lb(i) = k
        END IF
      END DO
    END DO

    ! Interpolation lineaire:

    DO i = 1, klon
      ! interpolation en logarithme de pression:

      ! ...   Modif . P. Le Van    ( 20/01/98) ....
      ! Modif Fr��ic Hourdin (3/01/02)

      aist(i) = log(pgcm(i,lb(i))/pres)/log(pgcm(i,lb(i))/pgcm(i,lt(i)))
      aisb(i) = log(pres/pgcm(i,lt(i)))/log(pgcm(i,lb(i))/pgcm(i,lt(i)))
    END DO


  END IF ! lnew

  ! ======================================================================
  ! inteprollation
  ! ======================================================================

  DO i = 1, klon
    qpres(i) = qgcm(i, lb(i))*aisb(i) + qgcm(i, lt(i))*aist(i)
  END DO

  ! Je mets les vents a zero quand je rencontre une montagne
  DO i = 1, klon
    IF (pgcm(i,1)<pres) THEN
      qpres(i) = missing_val
    END IF
  END DO


  RETURN
END SUBROUTINE plevel
