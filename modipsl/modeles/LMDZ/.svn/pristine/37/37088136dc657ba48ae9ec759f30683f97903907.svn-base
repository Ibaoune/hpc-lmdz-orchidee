
! $Header$

SUBROUTINE histo_o500_pctau(nbreg, pct_ocean, w, histo, histow, nhisto)
  USE dimphy
  IMPLICIT NONE

  INTEGER :: ij, k, l, nw
  INTEGER :: nreg, nbreg
  INTEGER, PARAMETER :: kmax = 8, lmax = 8
  INTEGER, PARAMETER :: kmaxm1 = kmax - 1, lmaxm1 = lmax - 1
  INTEGER, PARAMETER :: iwmax = 40

  INTEGER, DIMENSION (klon) :: iw
  REAL, DIMENSION (klon) :: w
  REAL, PARAMETER :: wmin = -200., pas_w = 10.
  REAL, DIMENSION (kmaxm1, lmaxm1, iwmax, nbreg) :: histow, nhisto
  REAL, DIMENSION (klon, kmaxm1, lmaxm1) :: histo

  ! LOGICAL, dimension(klon,nbreg) :: pct_ocean
  INTEGER, DIMENSION (klon, nbreg) :: pct_ocean

  ! initialisation
  histow(:, :, :, :) = 0.
  nhisto(:, :, :, :) = 0.

  ! calcul de l'histogramme de chaque regime dynamique
  DO nreg = 1, nbreg
    DO ij = 1, klon
      iw(ij) = int((w(ij)-wmin)/pas_w) + 1
      ! IF(pct_ocean(ij,nreg)) THEN
      ! IF(pct_ocean(ij,nreg).EQ.1) THEN
      IF (iw(ij)>=1 .AND. iw(ij)<=iwmax) THEN
        DO l = 1, lmaxm1
          DO k = 1, kmaxm1
            IF (histo(ij,k,l)>0.) THEN
              histow(k, l, iw(ij), nreg) = histow(k, l, iw(ij), nreg) + &
                histo(ij, k, l)*pct_ocean(ij, nreg)
              nhisto(k, l, iw(ij), nreg) = nhisto(k, l, iw(ij), nreg) + &
                pct_ocean(ij, nreg)
            END IF
          END DO !k
        END DO !l
        ! ELSE IF (iw(ij).LE.0.OR.iw(ij).GT.iwmax) THEN !iw
        ! PRINT*,'ij,iw=',ij,iw(ij)
      END IF !iw
      ! ENDIF !pct_ocean
    END DO !ij
    ! normalisation
    DO nw = 1, iwmax
      DO l = 1, lmaxm1
        DO k = 1, kmaxm1
          IF (nhisto(k,l,nw,nreg)/=0.) THEN
            histow(k, l, nw, nreg) = 100.*histow(k, l, nw, nreg)/ &
              nhisto(k, l, nw, nreg)
            ! PRINT*,'k,l,nw,nreg,histoW',k,l,nw,nreg,
            ! &     histoW(k,l,nw,nreg)
          END IF
        END DO !k
      END DO !l
    END DO !nw
  END DO !nreg

  RETURN
END SUBROUTINE histo_o500_pctau
