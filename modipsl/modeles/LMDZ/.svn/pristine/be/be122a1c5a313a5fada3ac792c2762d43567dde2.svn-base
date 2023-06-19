SUBROUTINE SRTM_CLDPROP &
 & ( KLEV, K_ICLDATM, K_INFLAG, K_ICEFLAG, K_LIQFLAG, K_NSTR,&
 & P_CLDFRAC, P_CLDDAT1, P_CLDDAT2, P_CLDDAT3, P_CLDDAT4, P_CLDDATMOM,&
 & P_TAUCLDORIG, P_TAUCLOUD, P_SSACLOUD, P_XMOM &
 & )  

!  path:      $Source: /storm/rc1/cvsroot/rc/rrtm_sw/src/cldprop_sw.f,v $
!  author:    $Author: jdelamer $
!  revision:  $Revision: 2.6 $
!  created:   $Date: 2002/04/04 18:29:47 $

!  PURPOSE:  COMPUTE THE CLOUD OPTICAL DEPTH(S) FOR EACH CLOUDY
!  LAYER.  NOTE:  ONLY INFLAG = 0 AND INFLAG=2/LIQFLAG=1,
!  ICEFLAG=3
!  (HU & STAMNES, Q. FU) ARE IMPLEMENTED.

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARSRTM ,  ONLY : JPLAY, JPBAND, JPB1, JPB2
USE YOESRTOP,  ONLY : EXTLIQ1, SSALIQ1, ASYLIQ1 &
 & , EXTICE3, SSAICE3, ASYICE3, FDLICE3 &
 & , FDELTA , EXTCOICE, SSACOICE, GICE, FORWICE &
 & , EXTCOLIQ, SSACOLIQ, GLIQ, FORWLIQ  

!     ------------------------------------------------------------------

IMPLICIT NONE

!-- real arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_ICLDATM 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INFLAG 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_ICEFLAG 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LIQFLAG 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_NSTR 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_CLDFRAC(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_CLDDAT1(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_CLDDAT2(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_CLDDAT3(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_CLDDAT4(JPLAY) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_CLDDATMOM(0:16,JPLAY) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUCLDORIG(JPLAY,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUCLOUD(JPLAY,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SSACLOUD(JPLAY,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_XMOM(0:16,JPLAY,JPBAND) 
!-- integer arguments

!-- real locals
REAL(KIND=JPRB) :: Z_EPS
REAL(KIND=JPRB) :: Z_TAUCLDORIG_A, Z_FFP, Z_FFP1, Z_FFPSSA, Z_SSACLOUD_A, Z_TAUCLOUD_A
REAL(KIND=JPRB) :: Z_CWP, Z_FICE, Z_RADICE, Z_FACTOR, Z_FINT, Z_FLIQ, Z_RADLIQ
REAL(KIND=JPRB) :: Z_TAUICEORIG, Z_SCATICE, Z_SSAICE, Z_TAUICE &
 & , Z_TAULIQORIG, Z_SCATLIQ, Z_SSALIQ, Z_TAULIQ  

!-- integer locals
INTEGER(KIND=JPIM) :: I_NCBANDS, I_NLAYERS
INTEGER(KIND=JPIM) :: IB, IB1, IB2, I_LAY, ISTR , INDEX
INTEGER(KIND=JPIM) :: I_NDBUG
REAL(KIND=JPRB) :: ZHOOK_HANDLE
 
!     INCLUDE 'param.f'

!      COMMON /CONTROL/   IAER, NSTR, IOUT, ISTART, IEND, ICLD 

!      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
!     &                   PZ(0:MXLAY),TZ(0:MXLAY)

!      COMMON /CLOUDIN/   INFLAG,CLDDAT1(MXLAY),CLDDAT2(MXLAY),
!     &                   ICEFLAG,LIQFLAG,CLDDAT3(MXLAY),CLDDAT4(MXLAY),
!     &                   CLDDATMOM(0:16,MXLAY)

!      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),
!     &     TAUCLOUD(MXLAY,NBANDS),SSACLOUD(MXLAY,NBANDS),
!     &     XMOM(0:16,MXLAY,NBANDS),TAUCLDORIG(MXLAY,NBANDS)

!      COMMON /HVERSN/    HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
!     *                   HVDUM1(4),HVRUTL,HVREXT,
!     *                   HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
!     *                   HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR

!      CHARACTER*15 HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
!     *            HVDUM1,HVRUTL,HVREXT,
!     *            HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
!     *            HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR

!      DIMENSION ABSICE0(2), ABSICE1(2,5), ABSICE2(40,16)
!      DIMENSION ABSLIQ1(58,16)

!      DIMENSION ABSCOICE(NBANDS), ABSCOLIQ(NBANDS)
!      DIMENSION EXTCOLIQ(NBANDS),SSACOLIQ(NBANDS),GLIQ(NBANDS)
!      DIMENSION EXTCOICE(NBANDS),SSACOICE(NBANDS),GICE(NBANDS)
!      DIMENSION FORWLIQ(NBANDS),FORWICE(NBANDS)

! Added for SW
!      DIMENSION EXTLIQ1(58,NBANDS), SSALIQ1(58,NBANDS), 
!     &     ASYLIQ1(58,NBANDS)
!      DIMENSION EXTICE3(27,NBANDS), SSAICE3(27,NBANDS),
!     &     ASYICE3(27,NBANDS),FDLICE3(28,NBANDS)
!      DIMENSION FDELTA(NBANDS)

!DATA EPS /1.E-6/

! HVRCLD = '$Revision: 2.6 $'

IF (LHOOK) CALL DR_HOOK('SRTM_CLDPROP',0,ZHOOK_HANDLE)
Z_EPS = 1.E-06_JPRB
I_NDBUG = 3

K_ICLDATM = 0
I_NCBANDS = 29
I_NLAYERS = KLEV
IB1     = JPB1
IB2     = JPB2

IF (I_NDBUG <= 2) THEN
  print *,'cldprop before loop K_INFLAG, K_ICEFLAG, K_LIQFLAG:',K_INFLAG,K_ICEFLAG,K_LIQFLAG,IB1,IB2
ENDIF

DO I_LAY = 1, I_NLAYERS

  IF (P_CLDFRAC(I_LAY) >= Z_EPS) THEN
    K_ICLDATM = 1

    IF (I_NDBUG <= 2) THEN
      print 9101,I_LAY,K_ICLDATM,P_CLDFRAC(I_LAY),P_CLDDAT1(I_LAY),P_CLDDAT2(I_LAY),P_CLDDAT3(I_LAY)&
       & ,P_CLDDAT4(I_LAY),(P_CLDDATMOM(ISTR,I_LAY),ISTR=0,K_NSTR)  
      9101  format(1x,'Cld :',2I3,f7.4,7E12.5)
    ENDIF

!   Ice clouds and water clouds combined.
    IF (K_INFLAG == 0) THEN
      Z_TAUCLDORIG_A = P_CLDDAT1(I_LAY)
      Z_FFP = P_CLDDATMOM(K_NSTR,I_LAY)
      Z_FFP1 = 1.0 - Z_FFP
      Z_FFPSSA = 1.0 - Z_FFP * P_CLDDAT2(I_LAY)
      Z_SSACLOUD_A = Z_FFP1*P_CLDDAT2(I_LAY)/Z_FFPSSA
      Z_TAUCLOUD_A = Z_FFPSSA*Z_TAUCLDORIG_A

!       DO IB = 16,NBANDS
      DO IB = IB1 , IB2
        P_TAUCLDORIG(I_LAY,IB) = Z_TAUCLDORIG_A
        P_SSACLOUD(I_LAY,IB) = Z_SSACLOUD_A
        P_TAUCLOUD(I_LAY,IB) = Z_TAUCLOUD_A

        DO ISTR = 0,K_NSTR
          P_XMOM(ISTR,I_LAY,IB) = (P_CLDDATMOM(ISTR,I_LAY) - Z_FFP)/ &
           & (Z_FFP1)  
        ENDDO
      ENDDO

!   Separate treatement of ice clouds and water clouds.
    ELSEIF(K_INFLAG == 2) THEN       
      Z_CWP = P_CLDDAT1(I_LAY)
      Z_FICE = P_CLDDAT2(I_LAY)
      Z_RADICE = P_CLDDAT3(I_LAY)

      IF (I_NDBUG <= 1) THEN
        print 9102,I_LAY,Z_CWP,Z_FICE,Z_RADICE
        9102     format(1x,'A',I3,3E13.6)
      ENDIF

!   Calculation of absorption coefficients due to ice clouds.
      IF (Z_FICE == 0.0) THEN
!         NCBANDS = 29
!         DO IB = 16, NCBANDS
        DO IB = IB1 , IB2
          EXTCOICE(IB) = 0.0_JPRB
          SSACOICE(IB) = 1.0_JPRB
          GICE(IB)     = 1.0_JPRB
          FORWICE(IB)  = 0.0_JPRB

          IF (I_NDBUG <= 1) THEN
            print 9103,I_LAY,Z_FICE,Z_CWP,Z_RADICE,IB,EXTCOICE(IB),SSACOICE(IB),GICE(IB),FORWICE(IB)
            9103         format(1x,'B',I3,F6.3,2E13.6,I3,4E12.5)
          ENDIF

        ENDDO

      ELSEIF (K_ICEFLAG == 3) THEN
        IF (Z_RADICE < 10.0 .OR. Z_RADICE > 140.0) STOP 'ICE EFFECTIVE SIZE OUT OF BOUNDS'
!         NCBANDS = 29
        Z_FACTOR = (Z_RADICE - 5._JPRB)/5._JPRB
        INDEX = INT(Z_FACTOR)
        IF (INDEX == 27) INDEX = 26
        Z_FINT = Z_FACTOR - REAL(INDEX)

!         DO IB=16,NCBANDS
        DO IB = IB1 , IB2
          EXTCOICE(IB) = Z_FICE * (EXTICE3(INDEX,IB) + Z_FINT * &
           & (EXTICE3(INDEX+1,IB) - EXTICE3(INDEX,IB)))  
          SSACOICE(IB) = SSAICE3(INDEX,IB) + Z_FINT * &
           & (SSAICE3(INDEX+1,IB) - SSAICE3(INDEX,IB))  
          GICE(IB) = ASYICE3(INDEX,IB) + Z_FINT * &
           & (ASYICE3(INDEX+1,IB) - ASYICE3(INDEX,IB))  
          FDELTA(IB) = FDLICE3(INDEX,IB) + Z_FINT * &
           & (FDLICE3(INDEX+1,IB) - FDLICE3(INDEX,IB))  
          if (fdelta(ib) < 0.0) STOP 'FDELTA LESS THAN 0.0'
          if (fdelta(ib) > 1.0) STOP 'FDELTA GT THAN 1.0'                     
          FORWICE(IB) = FDELTA(IB) + 0.5 / SSACOICE(IB)
! See Fu 1996 p. 2067 
          IF (FORWICE(IB) > GICE(IB)) FORWICE(IB) = GICE(IB)
! Check to ensure all calculated quantities are within physical limits.
          if (extcoice(ib) < 0.0_JPRB) STOP 'ICE EXTINCTION LESS THAN 0.0'
          if (ssacoice(ib) > 1.0_JPRB) STOP 'ICE SSA GRTR THAN 1.0'
          if (ssacoice(ib) < 0.0_JPRB) STOP 'ICE SSA LESS THAN 0.0'
          if (gice(ib) > 1.0_JPRB) STOP 'ICE ASYM GRTR THAN 1.0'
          if (gice(ib) < 0.0_JPRB) STOP 'ICE ASYM LESS THAN 0.0'

          IF (I_NDBUG <= 1) THEN
            print 9104,I_LAY,Z_FICE,Z_CWP,Z_RADICE,IB,EXTCOICE(IB),SSACOICE(IB),GICE(IB),FORWICE(IB),FDELTA(IB)
            9104         format(1x,'C',I3,F5.3,2E13.6,I3,5E12.5)
          ENDIF

        ENDDO
      ENDIF
      print *,'end of ice computations for I_LAY=',I_LAY
                  
!  Calculation of absorption coefficients due to water clouds.
      Z_FLIQ = 1. - Z_FICE
      IF (Z_FLIQ == 0.0) THEN
!         NCBANDS = 29
!         DO IB = 16, NCBANDS
        DO IB = IB1 , IB2
          EXTCOLIQ(IB) = 0.0
          SSACOLIQ(IB) = 1.0
          GLIQ(IB) = 1.0
          FORWLIQ(IB) = 0.0

          IF (I_NDBUG <= 1) THEN
            print 9105,I_LAY,Z_FLIQ,Z_CWP,IB,EXTCOLIQ(IB),SSACOLIQ(IB),GLIQ(IB),FORWLIQ(IB)
            9105         format(1x,'D',I3,F5.3,1E13.6,I3,4E12.5)
          ENDIF

        ENDDO

      ELSEIF (K_LIQFLAG == 1) THEN
        Z_RADLIQ = P_CLDDAT4(I_LAY)
        IF (Z_RADLIQ < 1.5 .OR. Z_RADLIQ > 60.) STOP 'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
        INDEX = INT(Z_RADLIQ - 1.5)
        IF (INDEX == 0) INDEX = 1
        IF (INDEX == 58) INDEX = 57
        Z_FINT = Z_RADLIQ - 1.5 - REAL(INDEX)
!         NCBANDS = 29

!         DO IB = 16, NCBANDS
        DO IB = IB1 , IB2
          EXTCOLIQ(IB) = Z_FLIQ * (EXTLIQ1(INDEX,IB) + Z_FINT * &
           & (EXTLIQ1(INDEX+1,IB) - (EXTLIQ1(INDEX,IB))))  
          SSACOLIQ(IB) = SSALIQ1(INDEX,IB) + Z_FINT * &
           & (SSALIQ1(INDEX+1,IB) - SSALIQ1(INDEX,IB))  
          GLIQ(IB) = ASYLIQ1(INDEX,IB) + Z_FINT * &
           & (ASYLIQ1(INDEX+1,IB) - ASYLIQ1(INDEX,IB))  
          FORWLIQ(IB) = GLIQ(IB)**K_NSTR
! Check to ensure all calculated quantities are within physical limits.
          if (extcoliq(ib) < 0.0_JPRB) STOP 'LIQUID EXTINCTION LESS THAN 0.0'
          if (ssacoliq(ib) > 1.0_JPRB) STOP 'LIQUID SSA GRTR THAN 1.0'
          if (ssacoliq(ib) < 0.0_JPRB) STOP 'LIQUID SSA LESS THAN 0.0'
          if (gliq(ib) > 1.0_JPRB) STOP 'LIQUID ASYM GRTR THAN 1.0'
          if (gliq(ib) < 0.0_JPRB) STOP 'LIQUID ASYM LESS THAN 0.0'

          IF (I_NDBUG <= 1) THEN
            print 9106,I_LAY,Z_FLIQ,Z_CWP,Z_RADLIQ,IB,EXTCOLIQ(IB),SSACOLIQ(IB),GLIQ(IB),FORWLIQ(IB)
            9106         format(1x,'E',I3,F5.3,2E13.6,I3,5E12.5)
          ENDIF

        ENDDO
      ENDIF

      IF (I_NDBUG <= 1) THEN
        print *,'end of liquid water computations for I_LAY=',I_LAY
      ENDIF

!       DO IB = 16, NCBANDS
      DO IB = IB1 , IB2
        Z_TAULIQORIG = Z_CWP * EXTCOLIQ(IB)
        Z_TAUICEORIG = Z_CWP * EXTCOICE(IB)
        P_TAUCLDORIG(I_LAY,IB) = Z_TAULIQORIG + Z_TAUICEORIG

        IF (I_NDBUG <= 1) THEN
          print 9107,IB,Z_TAULIQORIG,Z_TAUICEORIG,P_TAUCLDORIG(I_LAY,IB),Z_CWP &
           & ,EXTCOLIQ(IB),EXTCOICE(IB),SSACOLIQ(IB),SSACOICE(IB) &
           & ,FORWLIQ(IB),FORWICE(IB)  
          9107       format(1x,'F',I3,10E12.5)
        ENDIF

        Z_SSALIQ = SSACOLIQ(IB) * (1. - FORWLIQ(IB)) / &
         & (1. - FORWLIQ(IB) * SSACOLIQ(IB))  
        Z_TAULIQ =  (1. - FORWLIQ(IB) * SSACOLIQ(IB)) * &
         & Z_TAULIQORIG  
        Z_SSAICE = SSACOICE(IB) * (1. - FORWICE(IB)) / &
         & (1. - FORWICE(IB) * SSACOICE(IB))  
        Z_TAUICE =  (1. - FORWICE(IB) * SSACOICE(IB)) * &
         & Z_TAUICEORIG  
        Z_SCATLIQ = Z_SSALIQ * Z_TAULIQ
        Z_SCATICE = Z_SSAICE * Z_TAUICE
        P_TAUCLOUD(I_LAY,IB) = Z_TAULIQ + Z_TAUICE
        P_SSACLOUD(I_LAY,IB) = (Z_SCATLIQ + Z_SCATICE) / &
         & P_TAUCLOUD(I_LAY,IB)  
        P_XMOM(0,I_LAY,IB) = 1.0

        IF (I_NDBUG <= 1) THEN
          print 9108,IB,Z_TAULIQORIG,Z_TAUICEORIG,Z_SSALIQ,Z_TAULIQ,Z_SCATLIQ,Z_SSAICE,Z_TAUICE,Z_SCATICE
          9108       format(1x,'G',I3,8E13.6)
        ENDIF

        DO ISTR = 1, K_NSTR
!This commented code is the standard method for delta-m scaling. In accordance
!  with the 1996 Fu paper, equation A.3, the moments for ice were calculated
!  as in the uncommented code.
!                     XMOM(ISTR,LAY,IB) = (SCATLIQ *  &
!     &                    (GLIQ(IB)**ISTR - FORWLIQ(IB)) / &
!     &                    (1. - FORWLIQ(IB)) &
!     &                    + SCATICE * &
!     &                    (GICE(IB)**ISTR - FORWICE(IB)) /  &
!     &                    (1. - FORWICE(IB)))/(SCATLIQ + SCATICE)

          P_XMOM(ISTR,I_LAY,IB) = (1.0/(Z_SCATLIQ+Z_SCATICE))* &
           & (Z_SCATLIQ*(GLIQ(IB)**ISTR - FORWLIQ(IB)) / &
           & (1. - FORWLIQ(IB)) &
           & + Z_SCATICE * &
           & ((gice(ib)-forwice(ib))/(1.0-forwice(ib)))**ISTR)  
        ENDDO
      ENDDO

    ENDIF

  ENDIF

ENDDO

IF (I_NDBUG <= 1) THEN
  print *,'about to leave SRTM_CLDPROP'
ENDIF

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CLDPROP',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_CLDPROP

