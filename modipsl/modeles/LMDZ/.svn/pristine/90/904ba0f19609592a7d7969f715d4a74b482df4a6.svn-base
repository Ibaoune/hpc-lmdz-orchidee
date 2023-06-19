
! $Header$

SUBROUTINE ajsec(paprs, pplay, t, q, limbas, d_t, d_q)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: ajustement sec (adaptation du GCM du LMD)
  ! ======================================================================
  ! Arguments:
  ! t-------input-R- Temperature

  ! d_t-----output-R-Incrementation de la temperature
  ! ======================================================================
  include "YOMCST.h"
  REAL paprs(klon, klev+1), pplay(klon, klev)
  REAL t(klon, klev), q(klon, klev)
  REAL d_t(klon, klev), d_q(klon, klev)

  INTEGER limbas(klon), limhau ! les couches a ajuster

  LOGICAL mixq
  ! cc      PARAMETER (mixq=.TRUE.)
  PARAMETER (mixq=.FALSE.)

  REAL zh(klon, klev)
  REAL zho(klon, klev)
  REAL zq(klon, klev)
  REAL zpk(klon, klev)
  REAL zpkdp(klon, klev)
  REAL hm, sm, qm
  LOGICAL modif(klon), down
  INTEGER i, k, k1, k2

  ! Initialisation:

  ! ym
  limhau = klev

  DO k = 1, klev
    DO i = 1, klon
      d_t(i, k) = 0.0
      d_q(i, k) = 0.0
    END DO
  END DO
  ! ------------------------------------- detection des profils a modifier
  DO k = 1, limhau
    DO i = 1, klon
      zpk(i, k) = pplay(i, k)**rkappa
      zh(i, k) = rcpd*t(i, k)/zpk(i, k)
      zho(i, k) = zh(i, k)
      zq(i, k) = q(i, k)
    END DO
  END DO

  DO k = 1, limhau
    DO i = 1, klon
      zpkdp(i, k) = zpk(i, k)*(paprs(i,k)-paprs(i,k+1))
    END DO
  END DO

  DO i = 1, klon
    modif(i) = .FALSE.
  END DO
  DO k = 2, limhau
    DO i = 1, klon
      IF (.NOT. modif(i) .AND. k-1>limbas(i)) THEN
        IF (zh(i,k)<zh(i,k-1)) modif(i) = .TRUE.
      END IF
    END DO
  END DO
  ! ------------------------------------- correction des profils instables
  DO i = 1, klon
    IF (modif(i)) THEN
      k2 = limbas(i)
8000  CONTINUE
      k2 = k2 + 1
      IF (k2>limhau) GO TO 8001
      IF (zh(i,k2)<zh(i,k2-1)) THEN
        k1 = k2 - 1
        k = k1
        sm = zpkdp(i, k2)
        hm = zh(i, k2)
        qm = zq(i, k2)
8020    CONTINUE
        sm = sm + zpkdp(i, k)
        hm = hm + zpkdp(i, k)*(zh(i,k)-hm)/sm
        qm = qm + zpkdp(i, k)*(zq(i,k)-qm)/sm
        down = .FALSE.
        IF (k1/=limbas(i)) THEN
          IF (hm<zh(i,k1-1)) down = .TRUE.
        END IF
        IF (down) THEN
          k1 = k1 - 1
          k = k1
        ELSE
          IF ((k2==limhau)) GO TO 8021
          IF ((zh(i,k2+1)>=hm)) GO TO 8021
          k2 = k2 + 1
          k = k2
        END IF
        GO TO 8020
8021    CONTINUE
        ! ------------ nouveau profil : constant (valeur moyenne)
        DO k = k1, k2
          zh(i, k) = hm
          zq(i, k) = qm
        END DO
        k2 = k2 + 1
      END IF
      GO TO 8000
8001  CONTINUE
    END IF
  END DO

  DO k = 1, limhau
    DO i = 1, klon
      d_t(i, k) = (zh(i,k)-zho(i,k))*zpk(i, k)/rcpd
      d_q(i, k) = zq(i, k) - q(i, k)
    END DO
  END DO

  ! FH : les d_q et d_t sont maintenant calcules de facon a valoir
  ! effectivement 0. si on ne fait rien.

  ! IF (limbas.GT.1) THEN
  ! DO k = 1, limbas-1
  ! DO i = 1, klon
  ! d_t(i,k) = 0.0
  ! d_q(i,k) = 0.0
  ! ENDDO
  ! ENDDO
  ! ENDIF

  ! IF (limhau.LT.klev) THEN
  ! DO k = limhau+1, klev
  ! DO i = 1, klon
  ! d_t(i,k) = 0.0
  ! d_q(i,k) = 0.0
  ! ENDDO
  ! ENDDO
  ! ENDIF

  IF (.NOT. mixq) THEN
    DO k = 1, klev
      DO i = 1, klon
        d_q(i, k) = 0.0
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE ajsec

SUBROUTINE ajsec_convv2(paprs, pplay, t, q, d_t, d_q)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: ajustement sec (adaptation du GCM du LMD)
  ! ======================================================================
  ! Arguments:
  ! t-------input-R- Temperature

  ! d_t-----output-R-Incrementation de la temperature
  ! ======================================================================
  include "YOMCST.h"
  REAL paprs(klon, klev+1), pplay(klon, klev)
  REAL t(klon, klev), q(klon, klev)
  REAL d_t(klon, klev), d_q(klon, klev)

  INTEGER limbas, limhau ! les couches a ajuster
  ! cc      PARAMETER (limbas=klev-3, limhau=klev)
  ! ym      PARAMETER (limbas=1, limhau=klev)

  LOGICAL mixq
  ! cc      PARAMETER (mixq=.TRUE.)
  PARAMETER (mixq=.FALSE.)

  REAL zh(klon, klev)
  REAL zq(klon, klev)
  REAL zpk(klon, klev)
  REAL zpkdp(klon, klev)
  REAL hm, sm, qm
  LOGICAL modif(klon), down
  INTEGER i, k, k1, k2

  ! Initialisation:

  ! ym
  limbas = 1
  limhau = klev

  DO k = 1, klev
    DO i = 1, klon
      d_t(i, k) = 0.0
      d_q(i, k) = 0.0
    END DO
  END DO
  ! ------------------------------------- detection des profils a modifier
  DO k = limbas, limhau
    DO i = 1, klon
      zpk(i, k) = pplay(i, k)**rkappa
      zh(i, k) = rcpd*t(i, k)/zpk(i, k)
      zq(i, k) = q(i, k)
    END DO
  END DO

  DO k = limbas, limhau
    DO i = 1, klon
      zpkdp(i, k) = zpk(i, k)*(paprs(i,k)-paprs(i,k+1))
    END DO
  END DO

  DO i = 1, klon
    modif(i) = .FALSE.
  END DO
  DO k = limbas + 1, limhau
    DO i = 1, klon
      IF (.NOT. modif(i)) THEN
        IF (zh(i,k)<zh(i,k-1)) modif(i) = .TRUE.
      END IF
    END DO
  END DO
  ! ------------------------------------- correction des profils instables
  DO i = 1, klon
    IF (modif(i)) THEN
      k2 = limbas
8000  CONTINUE
      k2 = k2 + 1
      IF (k2>limhau) GO TO 8001
      IF (zh(i,k2)<zh(i,k2-1)) THEN
        k1 = k2 - 1
        k = k1
        sm = zpkdp(i, k2)
        hm = zh(i, k2)
        qm = zq(i, k2)
8020    CONTINUE
        sm = sm + zpkdp(i, k)
        hm = hm + zpkdp(i, k)*(zh(i,k)-hm)/sm
        qm = qm + zpkdp(i, k)*(zq(i,k)-qm)/sm
        down = .FALSE.
        IF (k1/=limbas) THEN
          IF (hm<zh(i,k1-1)) down = .TRUE.
        END IF
        IF (down) THEN
          k1 = k1 - 1
          k = k1
        ELSE
          IF ((k2==limhau)) GO TO 8021
          IF ((zh(i,k2+1)>=hm)) GO TO 8021
          k2 = k2 + 1
          k = k2
        END IF
        GO TO 8020
8021    CONTINUE
        ! ------------ nouveau profil : constant (valeur moyenne)
        DO k = k1, k2
          zh(i, k) = hm
          zq(i, k) = qm
        END DO
        k2 = k2 + 1
      END IF
      GO TO 8000
8001  CONTINUE
    END IF
  END DO

  DO k = limbas, limhau
    DO i = 1, klon
      d_t(i, k) = zh(i, k)*zpk(i, k)/rcpd - t(i, k)
      d_q(i, k) = zq(i, k) - q(i, k)
    END DO
  END DO

  IF (limbas>1) THEN
    DO k = 1, limbas - 1
      DO i = 1, klon
        d_t(i, k) = 0.0
        d_q(i, k) = 0.0
      END DO
    END DO
  END IF

  IF (limhau<klev) THEN
    DO k = limhau + 1, klev
      DO i = 1, klon
        d_t(i, k) = 0.0
        d_q(i, k) = 0.0
      END DO
    END DO
  END IF

  IF (.NOT. mixq) THEN
    DO k = 1, klev
      DO i = 1, klon
        d_q(i, k) = 0.0
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE ajsec_convv2
SUBROUTINE ajsec_old(paprs, pplay, t, d_t)
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: ajustement sec (adaptation du GCM du LMD)
  ! ======================================================================
  ! Arguments:
  ! t-------input-R- Temperature

  ! d_t-----output-R-Incrementation de la temperature
  ! ======================================================================
  include "YOMCST.h"
  REAL paprs(klon, klev+1), pplay(klon, klev)
  REAL t(klon, klev)
  REAL d_t(klon, klev)

  REAL local_h(klon, klev)
  REAL hm, sm
  LOGICAL modif(klon), down
  INTEGER i, l, l1, l2
  ! ------------------------------------- detection des profils a modifier
  DO i = 1, klon
    modif(i) = .FALSE.
  END DO

  DO l = 1, klev
    DO i = 1, klon
      local_h(i, l) = rcpd*t(i, l)/(pplay(i,l)**rkappa)
    END DO
  END DO

  DO l = 2, klev
    DO i = 1, klon
      IF (local_h(i,l)<local_h(i,l-1)) THEN
        modif(i) = .TRUE.
      ELSE
        modif(i) = modif(i)
      END IF
    END DO
  END DO
  ! ------------------------------------- correction des profils instables
  DO i = 1, klon
    IF (modif(i)) THEN
      l2 = 1
8000  CONTINUE
      l2 = l2 + 1
      IF (l2>klev) GO TO 8001
      IF (local_h(i,l2)<local_h(i,l2-1)) THEN
        l1 = l2 - 1
        l = l1
        sm = pplay(i, l2)**rkappa*(paprs(i,l2)-paprs(i,l2+1))
        hm = local_h(i, l2)
8020    CONTINUE
        sm = sm + pplay(i, l)**rkappa*(paprs(i,l)-paprs(i,l+1))
        hm = hm + pplay(i, l)**rkappa*(paprs(i,l)-paprs(i,l+1))*(local_h(i,l) &
          -hm)/sm
        down = .FALSE.
        IF (l1/=1) THEN
          IF (hm<local_h(i,l1-1)) THEN
            down = .TRUE.
          END IF
        END IF
        IF (down) THEN
          l1 = l1 - 1
          l = l1
        ELSE
          IF ((l2==klev)) GO TO 8021
          IF ((local_h(i,l2+1)>=hm)) GO TO 8021
          l2 = l2 + 1
          l = l2
        END IF
        GO TO 8020
8021    CONTINUE
        ! ------------ nouveau profil : constant (valeur moyenne)
        DO l = l1, l2
          local_h(i, l) = hm
        END DO
        l2 = l2 + 1
      END IF
      GO TO 8000
8001  CONTINUE
    END IF
  END DO

  DO l = 1, klev
    DO i = 1, klon
      d_t(i, l) = local_h(i, l)*(pplay(i,l)**rkappa)/rcpd - t(i, l)
    END DO
  END DO

  RETURN
END SUBROUTINE ajsec_old
