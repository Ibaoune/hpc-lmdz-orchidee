!
MODULE stdlevvar_mod
!
! This module contains main procedures for calculation
! of temperature, specific humidity and wind at a reference level
!
  USE cdrag_mod
  USE screenp_mod
  USE screenc_mod
  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!
!r original routine svn3623
!
      SUBROUTINE stdlevvar(klon, knon, nsrf, zxli, &
                           u1, v1, t1, q1, z1, &
                           ts1, qsurf, z0m, z0h, psol, pat1, &
                           t_2m, q_2m, t_10m, q_10m, u_10m, ustar)
      IMPLICIT NONE
!-------------------------------------------------------------------------
!
! Objet : calcul de la temperature et l'humidite relative a 2m et du 
!         module du vent a 10m a partir des relations de Dyer-Businger et
!         des equations de Louis.
!
! Reference : Hess, Colman et McAvaney (1995)        
!
! I. Musat, 01.07.2002
!
!AM On rajoute en sortie t et q a 10m pr le calcule d'hbtm2 dans clmain
!
!-------------------------------------------------------------------------
!
! klon----input-I- dimension de la grille physique (= nb_pts_latitude X nb_pts_longitude)
! knon----input-I- nombre de points pour un type de surface
! nsrf----input-I- indice pour le type de surface; voir indice_sol_mod.F90
! zxli----input-L- TRUE si calcul des cdrags selon Laurent Li
! u1------input-R- vent zonal au 1er niveau du modele
! v1------input-R- vent meridien au 1er niveau du modele
! t1------input-R- temperature de l'air au 1er niveau du modele
! q1------input-R- humidite relative au 1er niveau du modele
! z1------input-R- geopotentiel au 1er niveau du modele
! ts1-----input-R- temperature de l'air a la surface
! qsurf---input-R- humidite relative a la surface
! z0m, z0h---input-R- rugosite
! psol----input-R- pression au sol
! pat1----input-R- pression au 1er niveau du modele
!
! t_2m---output-R- temperature de l'air a 2m
! q_2m---output-R- humidite relative a 2m
! u_10m--output-R- vitesse du vent a 10m
!AM
! t_10m--output-R- temperature de l'air a 10m
! q_10m--output-R- humidite specifique a 10m
! ustar--output-R- u*
!
      INTEGER, intent(in) :: klon, knon, nsrf
      LOGICAL, intent(in) :: zxli
      REAL, dimension(klon), intent(in) :: u1, v1, t1, q1, z1, ts1
      REAL, dimension(klon), intent(in) :: qsurf, z0m, z0h
      REAL, dimension(klon), intent(in) :: psol, pat1
!
      REAL, dimension(klon), intent(out) :: t_2m, q_2m, ustar
      REAL, dimension(klon), intent(out) :: u_10m, t_10m, q_10m
!-------------------------------------------------------------------------
      include "flux_arp.h"
      include "YOMCST.h"
!IM PLUS
      include "YOETHF.h"
!
! Quelques constantes et options:
!
! RKAR : constante de von Karman
      REAL, PARAMETER :: RKAR=0.40
! niter : nombre iterations calcul "corrector"
!     INTEGER, parameter :: niter=6, ncon=niter-1
      INTEGER, parameter :: niter=2, ncon=niter-1
!
! Variables locales
      INTEGER :: i, n
      REAL :: zref
      REAL, dimension(klon) :: speed
! tpot : temperature potentielle
      REAL, dimension(klon) :: tpot
      REAL, dimension(klon) :: zri1, cdran
      REAL, dimension(klon) :: cdram, cdrah
! ri1 : nb. de Richardson entre la surface --> la 1ere couche
      REAL, dimension(klon) :: ri1 
      REAL, dimension(klon) :: testar, qstar
      REAL, dimension(klon) :: zdte, zdq   
! lmon : longueur de Monin-Obukhov selon Hess, Colman and McAvaney 
      DOUBLE PRECISION, dimension(klon) :: lmon
      DOUBLE PRECISION, parameter :: eps=1.0D-20
      REAL, dimension(klon) :: delu, delte, delq
      REAL, dimension(klon) :: u_zref, te_zref, q_zref  
      REAL, dimension(klon) :: temp, pref
      LOGICAL :: okri
      REAL, dimension(klon) :: u_zref_p, te_zref_p, temp_p, q_zref_p
!convertgence
      REAL, dimension(klon) :: te_zref_con, q_zref_con
      REAL, dimension(klon) :: u_zref_c, te_zref_c, temp_c, q_zref_c
      REAL, dimension(klon) :: ok_pred, ok_corr
!     REAL, dimension(klon) :: conv_te, conv_q
!------------------------------------------------------------------------- 
      DO i=1, knon
       speed(i)=SQRT(u1(i)**2+v1(i)**2)
       ri1(i) = 0.0
      ENDDO
!
      okri=.FALSE.
!      CALL coefcdrag(klon, knon, nsrf, zxli, &
! &                   speed, t1, q1, z1, psol, &
! &                   ts1, qsurf, rugos, okri, ri1,  &         
! &                   cdram, cdrah, cdran, zri1, pref)            
! Fuxing WANG, 04/03/2015, replace the coefcdrag by the merged version: cdrag

      CALL cdrag(knon, nsrf, &
 &                   speed, t1, q1, z1, &
 &                   psol, ts1, qsurf, z0m, z0h, &
 &                   cdram, cdrah, zri1, pref)

! --- special Dice: on force cdragm ( a defaut de forcer ustar) MPL 05082013
     IF (ok_prescr_ust) then
      DO i = 1, knon
       print *,'cdram avant=',cdram(i)
       cdram(i) = ust*ust/speed(i)/speed(i)
       print *,'cdram ust speed apres=',cdram(i),ust,speed
      ENDDO
     ENDIF
!
!---------Star variables----------------------------------------------------
!
      DO i = 1, knon
        ri1(i) = zri1(i)
        tpot(i) = t1(i)* (psol(i)/pat1(i))**RKAPPA
        ustar(i) = sqrt(cdram(i) * speed(i) * speed(i))
        zdte(i) = tpot(i) - ts1(i)
        zdq(i) = max(q1(i),0.0) - max(qsurf(i),0.0)
!
!
!IM BUG BUG BUG       zdte(i) = max(zdte(i),1.e-10)
        zdte(i) = sign(max(abs(zdte(i)),1.e-10),zdte(i))
!
        testar(i) = (cdrah(i) * zdte(i) * speed(i))/ustar(i)
        qstar(i) = (cdrah(i) * zdq(i) * speed(i))/ustar(i)
        lmon(i) = (ustar(i) * ustar(i) * tpot(i))/ &
 &                (RKAR * RG * testar(i))
      ENDDO
!
!----------First aproximation of variables at zref --------------------------
      zref = 2.0
      CALL screenp(klon, knon, nsrf, speed, tpot, q1, &
 &                 ts1, qsurf, z0m, lmon, &
 &                 ustar, testar, qstar, zref, &
 &                 delu, delte, delq)
!
      DO i = 1, knon
        u_zref(i) = delu(i)
        q_zref(i) = max(qsurf(i),0.0) + delq(i)
        te_zref(i) = ts1(i) + delte(i)
        temp(i) = te_zref(i) * (psol(i)/pat1(i))**(-RKAPPA)
        q_zref_p(i) = q_zref(i)
!       te_zref_p(i) = te_zref(i)
        temp_p(i) = temp(i)
      ENDDO
!
! Iteration of the variables at the reference level zref : corrector calculation ; see Hess & McAvaney, 1995 
!
      DO n = 1, niter
!
        okri=.TRUE.
        CALL screenc(klon, knon, nsrf, zxli, &
 &                   u_zref, temp, q_zref, zref, &
 &                   ts1, qsurf, z0m, z0h, psol, &           
 &                   ustar, testar, qstar, okri, ri1, &
 &                   pref, delu, delte, delq) 
!
        DO i = 1, knon
          u_zref(i) = delu(i)
          q_zref(i) = delq(i) + max(qsurf(i),0.0)
          te_zref(i) = delte(i) + ts1(i) 
!
! return to normal temperature
!
          temp(i) = te_zref(i) * (psol(i)/pref(i))**(-RKAPPA)
!         temp(i) = te_zref(i) - (zref* RG)/RCPD/ &
!                 (1 + RVTMP2 * max(q_zref(i),0.0))
!
!IM +++
!         IF(temp(i).GT.350.) THEN
!           WRITE(*,*) 'temp(i) GT 350 K !!',i,nsrf,temp(i)
!         ENDIF
!IM ---
!
        IF(n.EQ.ncon) THEN
          te_zref_con(i) = te_zref(i)
          q_zref_con(i) = q_zref(i)
        ENDIF 
!
        ENDDO 
!
      ENDDO 
!
! verifier le critere de convergence : 0.25% pour te_zref et 5% pour qe_zref
!
!       DO i = 1, knon
!         conv_te(i) = (te_zref(i) - te_zref_con(i))/te_zref_con(i)
!         conv_q(i) = (q_zref(i) - q_zref_con(i))/q_zref_con(i)
!IM +++
!         IF(abs(conv_te(i)).GE.0.0025.AND.abs(conv_q(i)).GE.0.05) THEN
!           PRINT*,'DIV','i=',i,te_zref_con(i),te_zref(i),conv_te(i), &
!           q_zref_con(i),q_zref(i),conv_q(i)
!         ENDIF
!IM ---
!       ENDDO
!
      DO i = 1, knon
        q_zref_c(i) = q_zref(i)
        temp_c(i) = temp(i)
!
!       IF(zri1(i).LT.0.) THEN
!         IF(nsrf.EQ.1) THEN
!           ok_pred(i)=1.
!           ok_corr(i)=0.
!         ELSE
!           ok_pred(i)=0.
!           ok_corr(i)=1.
!         ENDIF
!       ELSE
!         ok_pred(i)=0.
!         ok_corr(i)=1.
!       ENDIF
!
        ok_pred(i)=0.
        ok_corr(i)=1.
!
        t_2m(i) = temp_p(i) * ok_pred(i) + temp_c(i) * ok_corr(i)
        q_2m(i) = q_zref_p(i) * ok_pred(i) + q_zref_c(i) * ok_corr(i)
!IM +++
!       IF(n.EQ.niter) THEN
!       IF(t_2m(i).LT.t1(i).AND.t_2m(i).LT.ts1(i)) THEN
!         PRINT*,' BAD t2m LT ',i,nsrf,t_2m(i),t1(i),ts1(i) 
!       ELSEIF(t_2m(i).GT.t1(i).AND.t_2m(i).GT.ts1(i)) THEN
!         PRINT*,' BAD t2m GT ',i,nsrf,t_2m(i),t1(i),ts1(i) 
!       ENDIF
!       ENDIF
!IM ---
      ENDDO
!
!
!----------First aproximation of variables at zref --------------------------
!
      zref = 10.0
      CALL screenp(klon, knon, nsrf, speed, tpot, q1, &
 &                 ts1, qsurf, z0m, lmon, &
 &                 ustar, testar, qstar, zref, &
 &                 delu, delte, delq)
!
      DO i = 1, knon
        u_zref(i) = delu(i)
        q_zref(i) = max(qsurf(i),0.0) + delq(i)
        te_zref(i) = ts1(i) + delte(i)
        temp(i) = te_zref(i) * (psol(i)/pat1(i))**(-RKAPPA)
!       temp(i) = te_zref(i) - (zref* RG)/RCPD/ &
!                 (1 + RVTMP2 * max(q_zref(i),0.0))
        u_zref_p(i) = u_zref(i)
      ENDDO
!
! Iteration of the variables at the reference level zref : corrector ; see Hess & McAvaney, 1995 
!
      DO n = 1, niter
!
        okri=.TRUE.
        CALL screenc(klon, knon, nsrf, zxli, &
 &                   u_zref, temp, q_zref, zref, &
 &                   ts1, qsurf, z0m, z0h, psol, &
 &                   ustar, testar, qstar, okri, ri1, &
 &                   pref, delu, delte, delq)
!
        DO i = 1, knon
          u_zref(i) = delu(i)
          q_zref(i) = delq(i) + max(qsurf(i),0.0)
          te_zref(i) = delte(i) + ts1(i)
          temp(i) = te_zref(i) * (psol(i)/pref(i))**(-RKAPPA)
!         temp(i) = te_zref(i) - (zref* RG)/RCPD/ &
!                   (1 + RVTMP2 * max(q_zref(i),0.0))
        ENDDO 
!
      ENDDO
!
      DO i = 1, knon
        u_zref_c(i) = u_zref(i)
!
        u_10m(i) = u_zref_p(i) * ok_pred(i) + u_zref_c(i) * ok_corr(i)
!
!AM
        q_zref_c(i) = q_zref(i)
        temp_c(i) = temp(i)
        t_10m(i) = temp_p(i) * ok_pred(i) + temp_c(i) * ok_corr(i)
        q_10m(i) = q_zref_p(i) * ok_pred(i) + q_zref_c(i) * ok_corr(i)
!MA
      ENDDO
! 
      RETURN
      END subroutine stdlevvar
!
      SUBROUTINE stdlevvarn(klon, knon, nsrf, zxli, &
                           u1, v1, t1, q1, z1, &
                           ts1, qsurf, z0m, z0h, psol, pat1, &
                           t_2m, q_2m, t_10m, q_10m, u_10m, ustar, &
                           n2mout)
!
      USE ioipsl_getin_p_mod, ONLY : getin_p
      IMPLICIT NONE
!-------------------------------------------------------------------------
!
! Objet : calcul de la temperature et l'humidite relative a 2m et du 
!         module du vent a 10m a partir des relations de Dyer-Businger et
!         des equations de Louis.
!
! Reference : Hess, Colman et McAvaney (1995)        
!
! I. Musat, 01.07.2002
!
!AM On rajoute en sortie t et q a 10m pr le calcule d'hbtm2 dans clmain
!
!-------------------------------------------------------------------------
!
! klon----input-I- dimension de la grille physique (= nb_pts_latitude X nb_pts_longitude)
! knon----input-I- nombre de points pour un type de surface
! nsrf----input-I- indice pour le type de surface; voir indice_sol_mod.F90
! zxli----input-L- TRUE si calcul des cdrags selon Laurent Li
! u1------input-R- vent zonal au 1er niveau du modele
! v1------input-R- vent meridien au 1er niveau du modele
! t1------input-R- temperature de l'air au 1er niveau du modele
! q1------input-R- humidite relative au 1er niveau du modele
! z1------input-R- geopotentiel au 1er niveau du modele
! ts1-----input-R- temperature de l'air a la surface
! qsurf---input-R- humidite relative a la surface
! z0m, z0h---input-R- rugosite
! psol----input-R- pression au sol
! pat1----input-R- pression au 1er niveau du modele
!
! t_2m---output-R- temperature de l'air a 2m
! q_2m---output-R- humidite relative a 2m
! u_2m--output-R- vitesse du vent a 2m
! u_10m--output-R- vitesse du vent a 10m
! ustar--output-R- u*
!AM
! t_10m--output-R- temperature de l'air a 10m
! q_10m--output-R- humidite specifique a 10m
!
      INTEGER, intent(in) :: klon, knon, nsrf
      LOGICAL, intent(in) :: zxli
      REAL, dimension(klon), intent(in) :: u1, v1, t1, q1, z1, ts1
      REAL, dimension(klon), intent(in) :: qsurf, z0m, z0h
      REAL, dimension(klon), intent(in) :: psol, pat1
!
      REAL, dimension(klon), intent(out) :: t_2m, q_2m, ustar
      REAL, dimension(klon), intent(out) :: u_10m, t_10m, q_10m
      INTEGER, dimension(klon, 6), intent(out) :: n2mout
!
      REAL, dimension(klon) :: u_2m
      REAL, dimension(klon) :: cdrm2m, cdrh2m, ri2m
      REAL, dimension(klon) :: cdram, cdrah, zri1
      REAL, dimension(klon) :: cdmn1, cdhn1, fm1, fh1
      REAL, dimension(klon) :: cdmn2m, cdhn2m, fm2m, fh2m
      REAL, dimension(klon) :: ri2m_new
!-------------------------------------------------------------------------
      include "flux_arp.h"
      include "YOMCST.h"
!IM PLUS
      include "YOETHF.h"
!
! Quelques constantes et options:
!
! RKAR : constante de von Karman
      REAL, PARAMETER :: RKAR=0.40
! niter : nombre iterations calcul "corrector"
!     INTEGER, parameter :: niter=6, ncon=niter-1
!IM 071020     INTEGER, parameter :: niter=2, ncon=niter-1
      INTEGER, parameter :: niter=2, ncon=niter
!     INTEGER, parameter :: niter=6, ncon=niter
!
! Variables locales
      INTEGER :: i, n
      REAL :: zref
      REAL, dimension(klon) :: speed
! tpot : temperature potentielle
      REAL, dimension(klon) :: tpot
      REAL, dimension(klon) :: cdran
! ri1 : nb. de Richardson entre la surface --> la 1ere couche
      REAL, dimension(klon) :: ri1 
      DOUBLE PRECISION, parameter :: eps=1.0D-20
      REAL, dimension(klon) :: delu, delte, delq
      REAL, dimension(klon) :: delh, delm
      REAL, dimension(klon) :: delh_new, delm_new
      REAL, dimension(klon) :: u_zref, te_zref, q_zref  
      REAL, dimension(klon) :: u_zref_pnew, te_zref_pnew, q_zref_pnew
      REAL, dimension(klon) :: temp, pref
      REAL, dimension(klon) :: temp_new, pref_new
      LOGICAL :: okri
      REAL, dimension(klon) :: u_zref_p, te_zref_p, temp_p, q_zref_p
      REAL, dimension(klon) :: u_zref_p_new, te_zref_p_new, temp_p_new, q_zref_p_new
!convergence
      REAL, dimension(klon) :: te_zref_con, q_zref_con
      REAL, dimension(klon) :: u_zref_c, te_zref_c, temp_c, q_zref_c
      REAL, dimension(klon) :: ok_pred, ok_corr
!
      REAL, dimension(klon) :: cdrm10m, cdrh10m, ri10m
      REAL, dimension(klon) :: cdmn10m, cdhn10m, fm10m, fh10m
      REAL, dimension(klon) :: cdn2m, cdn1
      REAL :: CEPDUE,zdu2 
      INTEGER :: nzref, nz1
      LOGICAL, dimension(klon) :: ok_t2m_toosmall, ok_t2m_toobig
      LOGICAL, dimension(klon) :: ok_q2m_toosmall, ok_q2m_toobig
      LOGICAL, dimension(klon) :: ok_u2m_toobig
      LOGICAL, dimension(klon) :: ok_t10m_toosmall, ok_t10m_toobig
      LOGICAL, dimension(klon) :: ok_q10m_toosmall, ok_q10m_toobig
      LOGICAL, dimension(klon) :: ok_u10m_toobig
      INTEGER, dimension(klon, 6) :: n10mout

!------------------------------------------------------------------------- 
      CEPDUE=0.1
!
! n2mout : compteur des pas de temps ou t2m,q2m ou u2m sont en dehors des intervalles
!          [tsurf, temp], [qsurf, q1] ou [0, speed]
! n10mout : compteur des pas de temps ou t10m,q10m ou u10m sont en dehors des intervalles
!          [tsurf, temp], [qsurf, q1] ou [0, speed]
!
      n2mout(:,:)=0
      n10mout(:,:)=0
      
      DO i=1, knon
       speed(i)=MAX(SQRT(u1(i)**2+v1(i)**2),CEPDUE)
       ri1(i) = 0.0
      ENDDO
!
      okri=.FALSE.
      CALL cdrag(knon, nsrf, &
 &                   speed, t1, q1, z1, &
 &                   psol, ts1, qsurf, z0m, z0h, &
 &                   cdram, cdrah, zri1, pref)

!
      DO i = 1, knon
        ri1(i) = zri1(i)
        tpot(i) = t1(i)* (psol(i)/pat1(i))**RKAPPA
        zdu2 = MAX(CEPDUE*CEPDUE, speed(i)**2)
        ustar(i) = sqrt(cdram(i) * zdu2)
!
      ENDDO
!
!----------First aproximation of variables at zref --------------------------
      zref = 2.0
!
! calcul first-guess en utilisant dans les calculs à 2m
! le Richardson de la premiere couche atmospherique 
!
       CALL screencn(klon, knon, nsrf, zxli, &
 &                   speed, tpot, q1, zref, &
 &                   ts1, qsurf, z0m, z0h, psol, &           
 &                   cdram, cdrah,  okri, &
 &                   ri1, 1, &
 &                   pref_new, delm_new, delh_new, ri2m)
!
       DO i = 1, knon
         u_zref(i) = delm_new(i)*speed(i)
         u_zref_p(i) = u_zref(i)
         q_zref(i) = delh_new(i)*max(q1(i),0.0) + &
         &           max(qsurf(i),0.0)*(1-delh_new(i))
         q_zref_p(i) = q_zref(i)
         te_zref(i) = delh_new(i)*tpot(i) + ts1(i)*(1-delh_new(i))
         te_zref_p(i) = te_zref(i)
!
! return to normal temperature
         temp(i) = te_zref(i) * (psol(i)/pref_new(i))**(-RKAPPA)
         temp_p(i) = temp(i)
!
! compteurs ici
!
         ok_t2m_toosmall(i)=te_zref(i).LT.tpot(i).AND. &
         & te_zref(i).LT.ts1(i)
         ok_t2m_toobig(i)=te_zref(i).GT.tpot(i).AND. &
         & te_zref(i).GT.ts1(i)
         ok_q2m_toosmall(i)=q_zref(i).LT.q1(i).AND. &
         & q_zref(i).LT.qsurf(i)
         ok_q2m_toobig(i)=q_zref(i).GT.q1(i).AND. &
         & q_zref(i).GT.qsurf(i)
         ok_u2m_toobig(i)=u_zref(i).GT.speed(i)
!
         IF(ok_t2m_toosmall(i).OR.ok_t2m_toobig(i)) THEN
             n2mout(i,1)=n2mout(i,1)+1
         ENDIF
         IF(ok_q2m_toosmall(i).OR.ok_q2m_toobig(i)) THEN
             n2mout(i,3)=n2mout(i,3)+1
         ENDIF
         IF(ok_u2m_toobig(i)) THEN
             n2mout(i,5)=n2mout(i,5)+1
         ENDIF
!
         IF(ok_t2m_toosmall(i).OR.ok_t2m_toobig(i).OR. &
          & ok_q2m_toosmall(i).OR.ok_q2m_toobig(i).OR. &
          & ok_u2m_toobig(i)) THEN 
             delm_new(i)=min(max(delm_new(i),0.),1.)
             delh_new(i)=min(max(delh_new(i),0.),1.)
             u_zref(i) = delm_new(i)*speed(i)
             u_zref_p(i) = u_zref(i)
             q_zref(i) = delh_new(i)*max(q1(i),0.0) + &
         &               max(qsurf(i),0.0)*(1-delh_new(i))
             q_zref_p(i) = q_zref(i)
             te_zref(i) = delh_new(i)*tpot(i) + ts1(i)*(1-delh_new(i))
             te_zref_p(i) = te_zref(i)
!
! return to normal temperature
             temp(i) = te_zref(i) * (psol(i)/pref_new(i))**(-RKAPPA)
             temp_p(i) = temp(i)
         ENDIF 
!
       ENDDO
!
! Iteration of the variables at the reference level zref : corrector calculation ; see Hess & McAvaney, 1995 
!
      DO n = 1, niter
!
        okri=.TRUE.
        CALL screencn(klon, knon, nsrf, zxli, &
 &                   u_zref, temp, q_zref, zref, &
 &                   ts1, qsurf, z0m, z0h, psol, &
 &                   cdram, cdrah,  okri, &
 &                   ri1, 0, &
 &                   pref, delm, delh, ri2m)
!
        DO i = 1, knon
          u_zref(i) = delm(i)*speed(i)
          q_zref(i) = delh(i)*max(q1(i),0.0) + &
          &           max(qsurf(i),0.0)*(1-delh(i))
          te_zref(i) = delh(i)*tpot(i) + ts1(i)*(1-delh(i))
!
! return to normal temperature
          temp(i) = te_zref(i) * (psol(i)/pref(i))**(-RKAPPA)
!
! compteurs ici
!
          ok_t2m_toosmall(i)=te_zref(i).LT.tpot(i).AND. &
          & te_zref(i).LT.ts1(i)
          ok_t2m_toobig(i)=te_zref(i).GT.tpot(i).AND. &
          & te_zref(i).GT.ts1(i)
          ok_q2m_toosmall(i)=q_zref(i).LT.q1(i).AND. &
          & q_zref(i).LT.qsurf(i)
          ok_q2m_toobig(i)=q_zref(i).GT.q1(i).AND. &
          & q_zref(i).GT.qsurf(i)
          ok_u2m_toobig(i)=u_zref(i).GT.speed(i)
!
          IF(ok_t2m_toosmall(i).OR.ok_t2m_toobig(i)) THEN
              n2mout(i,2)=n2mout(i,2)+1
          ENDIF
          IF(ok_q2m_toosmall(i).OR.ok_q2m_toobig(i)) THEN
              n2mout(i,4)=n2mout(i,4)+1
          ENDIF
          IF(ok_u2m_toobig(i)) THEN
              n2mout(i,6)=n2mout(i,6)+1
          ENDIF
!
          IF(ok_t2m_toosmall(i).OR.ok_t2m_toobig(i).OR. &
           & ok_q2m_toosmall(i).OR.ok_q2m_toobig(i).OR. &
           & ok_u2m_toobig(i)) THEN 
              delm(i)=min(max(delm(i),0.),1.)
              delh(i)=min(max(delh(i),0.),1.)
              u_zref(i) = delm(i)*speed(i)
              q_zref(i) = delh(i)*max(q1(i),0.0) + &
          &           max(qsurf(i),0.0)*(1-delh(i))
              te_zref(i) = delh(i)*tpot(i) + ts1(i)*(1-delh(i))
              temp(i) = te_zref(i) * (psol(i)/pref(i))**(-RKAPPA)
          ENDIF
!
!
          IF(n.EQ.ncon) THEN
            te_zref_con(i) = te_zref(i)
            q_zref_con(i) = q_zref(i)
          ENDIF 
!
        ENDDO 
!
      ENDDO 
!
      DO i = 1, knon
        q_zref_c(i) = q_zref(i)
        temp_c(i) = temp(i)
!
        ok_pred(i)=0.
        ok_corr(i)=1.
!
        t_2m(i) = temp_p(i) * ok_pred(i) + temp_c(i) * ok_corr(i)
        q_2m(i) = q_zref_p(i) * ok_pred(i) + q_zref_c(i) * ok_corr(i)
!
        u_zref_c(i) = u_zref(i)
        u_2m(i) = u_zref_p(i) * ok_pred(i) + u_zref_c(i) * ok_corr(i)
      ENDDO
!
!
!----------First aproximation of variables at zref --------------------------
!
      zref = 10.0
!
       CALL screencn(klon, knon, nsrf, zxli, &
 &                   speed, tpot, q1, zref, &
 &                   ts1, qsurf, z0m, z0h, psol, &           
 &                   cdram, cdrah,  okri, &
 &                   ri1, 1, &
 &                   pref_new, delm_new, delh_new, ri10m)
!
       DO i = 1, knon
         u_zref(i) = delm_new(i)*speed(i)
         q_zref(i) = delh_new(i)*max(q1(i),0.0) + &
         &           max(qsurf(i),0.0)*(1-delh_new(i))
         te_zref(i) = delh_new(i)*tpot(i) + ts1(i)*(1-delh_new(i))
         temp(i) = te_zref(i) * (psol(i)/pref_new(i))**(-RKAPPA)
         u_zref_p(i) = u_zref(i)
!
! compteurs ici
!
         ok_t10m_toosmall(i)=te_zref(i).LT.tpot(i).AND. &
         & te_zref(i).LT.ts1(i)
         ok_t10m_toobig(i)=te_zref(i).GT.tpot(i).AND. &
         & te_zref(i).GT.ts1(i)
         ok_q10m_toosmall(i)=q_zref(i).LT.q1(i).AND. &
         & q_zref(i).LT.qsurf(i)
         ok_q10m_toobig(i)=q_zref(i).GT.q1(i).AND. &
         & q_zref(i).GT.qsurf(i)
         ok_u10m_toobig(i)=u_zref(i).GT.speed(i)
!
         IF(ok_t10m_toosmall(i).OR.ok_t10m_toobig(i)) THEN
             n10mout(i,1)=n10mout(i,1)+1
         ENDIF
         IF(ok_q10m_toosmall(i).OR.ok_q10m_toobig(i)) THEN
             n10mout(i,3)=n10mout(i,3)+1
         ENDIF
         IF(ok_u10m_toobig(i)) THEN
             n10mout(i,5)=n10mout(i,5)+1
         ENDIF
!
         IF(ok_t10m_toosmall(i).OR.ok_t10m_toobig(i).OR. &
          & ok_q10m_toosmall(i).OR.ok_q10m_toobig(i).OR. &
          & ok_u10m_toobig(i)) THEN 
             delm_new(i)=min(max(delm_new(i),0.),1.)
             delh_new(i)=min(max(delh_new(i),0.),1.)
             u_zref(i) = delm_new(i)*speed(i)
             u_zref_p(i) = u_zref(i)
             q_zref(i) = delh_new(i)*max(q1(i),0.0) + &
         &               max(qsurf(i),0.0)*(1-delh_new(i))
             te_zref(i) = delh_new(i)*tpot(i) + ts1(i)*(1-delh_new(i))
             temp(i) = te_zref(i) * (psol(i)/pref_new(i))**(-RKAPPA)
         ENDIF 
!
       ENDDO
!
! Iteration of the variables at the reference level zref : corrector calculation ; see Hess & McAvaney, 1995 
!
      DO n = 1, niter
!
        okri=.TRUE.
        CALL screencn(klon, knon, nsrf, zxli, &
 &                   u_zref, temp, q_zref, zref, &
 &                   ts1, qsurf, z0m, z0h, psol, &
 &                   cdram, cdrah,  okri, &
 &                   ri1, 0, &
 &                   pref, delm, delh, ri10m)
!
        DO i = 1, knon
          u_zref(i) = delm(i)*speed(i)
          q_zref(i) = delh(i)*max(q1(i),0.0) + &
          &           max(qsurf(i),0.0)*(1-delh(i))
          te_zref(i) = delh(i)*tpot(i) + ts1(i)*(1-delh(i))
!
! return to normal temperature
          temp(i) = te_zref(i) * (psol(i)/pref(i))**(-RKAPPA)
!
! compteurs ici
!
          ok_t10m_toosmall(i)=te_zref(i).LT.tpot(i).AND. &
          & te_zref(i).LT.ts1(i)
          ok_t10m_toobig(i)=te_zref(i).GT.tpot(i).AND. &
          & te_zref(i).GT.ts1(i)
          ok_q10m_toosmall(i)=q_zref(i).LT.q1(i).AND. &
          & q_zref(i).LT.qsurf(i)
          ok_q10m_toobig(i)=q_zref(i).GT.q1(i).AND. &
          & q_zref(i).GT.qsurf(i)
          ok_u10m_toobig(i)=u_zref(i).GT.speed(i)
!
          IF(ok_t10m_toosmall(i).OR.ok_t10m_toobig(i)) THEN
              n10mout(i,2)=n10mout(i,2)+1
          ENDIF
          IF(ok_q10m_toosmall(i).OR.ok_q10m_toobig(i)) THEN
              n10mout(i,4)=n10mout(i,4)+1
          ENDIF
          IF(ok_u10m_toobig(i)) THEN
              n10mout(i,6)=n10mout(i,6)+1
          ENDIF
!
          IF(ok_t10m_toosmall(i).OR.ok_t10m_toobig(i).OR. &
           & ok_q10m_toosmall(i).OR.ok_q10m_toobig(i).OR. &
           & ok_u10m_toobig(i)) THEN 
              delm(i)=min(max(delm(i),0.),1.)
              delh(i)=min(max(delh(i),0.),1.)
              u_zref(i) = delm(i)*speed(i)
              q_zref(i) = delh(i)*max(q1(i),0.0) + &
          &           max(qsurf(i),0.0)*(1-delh(i))
              te_zref(i) = delh(i)*tpot(i) + ts1(i)*(1-delh(i))
              temp(i) = te_zref(i) * (psol(i)/pref(i))**(-RKAPPA)
          ENDIF
!
!
          IF(n.EQ.ncon) THEN
            te_zref_con(i) = te_zref(i)
            q_zref_con(i) = q_zref(i)
          ENDIF 
!
        ENDDO 
!
      ENDDO 
!
      DO i = 1, knon
        q_zref_c(i) = q_zref(i)
        temp_c(i) = temp(i)
!
        ok_pred(i)=0.
        ok_corr(i)=1.
!
        t_10m(i) = temp_p(i) * ok_pred(i) + temp_c(i) * ok_corr(i)
        q_10m(i) = q_zref_p(i) * ok_pred(i) + q_zref_c(i) * ok_corr(i)
!
        u_zref_c(i) = u_zref(i)
        u_10m(i) = u_zref_p(i) * ok_pred(i) + u_zref_c(i) * ok_corr(i)
      ENDDO
!
      RETURN
      END subroutine stdlevvarn

END MODULE stdlevvar_mod
