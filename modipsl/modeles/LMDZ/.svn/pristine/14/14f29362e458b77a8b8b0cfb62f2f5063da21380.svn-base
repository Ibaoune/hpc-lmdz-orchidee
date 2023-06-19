!
MODULE cdrag_mod
!
! This module contains some procedures for calculation of the cdrag
! coefficients for turbulent diffusion at surface
!
  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!
!r original routine svn3623
!
 SUBROUTINE cdrag(knon,  nsrf,   &
     speed, t1,    q1,    zgeop1, &
     psol,  tsurf, qsurf, z0m, z0h,  &
     cdm,  cdh,  zri,   pref)

  USE dimphy
  USE indice_sol_mod
  USE print_control_mod, ONLY: lunout, prt_level
  USE ioipsl_getin_p_mod, ONLY : getin_p

  IMPLICIT NONE
! ================================================================= c
!
! Objet : calcul des cdrags pour le moment (pcfm) et 
!         les flux de chaleur sensible et latente (pcfh) d'apr??s
!         Louis 1982, Louis 1979, King et al 2001
!         ou Zilitinkevich et al 2002  pour les cas stables, Louis 1979
!         et 1982 pour les cas instables 
!
! Modified history: 
!  writting on the 20/05/2016
!  modified on the 13/12/2016 to be adapted to LMDZ6
!
! References:
!   Louis, J. F., 1979: A parametric model of vertical eddy fluxes in the
!     atmosphere. Boundary-Layer Meteorology. 01/1979; 17(2):187-202. 
!   Louis, J. F., Tiedtke, M. and Geleyn, J. F., 1982: `A short history of the
!     operational PBL parametrization at ECMWF'. Workshop on boundary layer
!     parametrization, November 1981, ECMWF, Reading, England. 
!     Page: 19. Equations in Table 1.
!   Mascart P, Noilhan J, Giordani H 1995.A MODIFIED PARAMETERIZATION OF FLUX-PROFILE RELATIONSHIPS 
!    IN THE SURFACE LAYER USING DIFFERENT ROUGHNESS LENGTH VALUES FOR HEAT AND MOMENTUM
!    Boundary-Layer Meteorology 72: 331-344
!   Anton Beljaars. May 1992. The parametrization of the planetary boundary layer.  
!     European Centre for Medium-Range Weather Forecasts.
!     Equations: 110-113. Page 40.
!   Miller,M.J., A.C.M.Beljaars, T.N.Palmer. 1992. The sensitivity of the ECMWF
!     model to the parameterization of evaporation from the tropical oceans. J.
!     Climate, 5:418-434.
!   King J.C, Connolley, W.M ad Derbyshire S.H. 2001, Sensitivity of Modelled Antarctic climate
!   to surface and boundary-layer flux parametrizations
!   QJRMS, 127, pp 779-794
!
! ================================================================= c
! ================================================================= c
! On choisit le couple de fonctions de correction avec deux flags:
! Un pour les cas instables, un autre pour les cas stables
!
! iflag_corr_insta:
!                   1: Louis 1979 avec les modifications de Mascart 1995 (z0/= z0h)
!                   2: Louis 1982 
!                   3: Laurent Li
!
! iflag_corr_sta:
!                   1: Louis 1979 avec les modifications de Mascart 1995 (z0/= z0h)
!                   2: Louis 1982 
!                   3: Laurent Li
!                   4: King  2001 (SHARP)
!                   5: MO 1st order theory (allow collapse of turbulence)
!            
!
!*****************************************************************
! Parametres d'entree
!***************************************************************** 
  
  INTEGER, INTENT(IN)                      :: knon, nsrf ! nombre de points de grille sur l'horizontal + type de surface
  REAL, DIMENSION(klon), INTENT(IN)        :: speed ! module du vent au 1er niveau du modele
  REAL, DIMENSION(klon), INTENT(IN)        :: zgeop1! geopotentiel au 1er niveau du modele
  REAL, DIMENSION(klon), INTENT(IN)        :: tsurf ! Surface temperature (K)
  REAL, DIMENSION(klon), INTENT(IN)        :: qsurf ! Surface humidity (Kg/Kg)
  REAL, DIMENSION(klon), INTENT(IN)        :: z0m, z0h ! Rugosity at surface (m) 
  REAL, DIMENSION(klon), INTENT(IN)        :: t1  ! Temperature au premier niveau (K) 
  REAL, DIMENSION(klon), INTENT(IN)        :: q1  ! humidite specifique au premier niveau (kg/kg)
  REAL, DIMENSION(klon), INTENT(IN)        :: psol ! pression au sol



! Parametres de sortie
!******************************************************************
  REAL, DIMENSION(klon), INTENT(OUT)       :: cdm   ! Drag coefficient for momentum
  REAL, DIMENSION(klon), INTENT(OUT)       :: cdh   ! Drag coefficient for heat flux
  REAL, DIMENSION(klon), INTENT(OUT)       :: zri   ! Richardson number
  REAL, DIMENSION(klon), INTENT(OUT)       :: pref  ! Pression au niveau zgeop/RG

! Variables Locales
!******************************************************************
 

  INCLUDE "YOMCST.h"
  INCLUDE "YOETHF.h"
  INCLUDE "clesphys.h"


  REAL, PARAMETER       :: CKAP=0.40, CKAPT=0.42
  REAL                   CEPDU2
  REAL                   ALPHA
  REAL                   CB,CC,CD,C2,C3
  REAL                   MU, CM, CH, B, CMstar, CHstar
  REAL                   PM, PH, BPRIME
  REAL                   C
  INTEGER                ng_q1                      ! Number of grids that q1 < 0.0
  INTEGER                ng_qsurf                   ! Number of grids that qsurf < 0.0
  INTEGER                i
  REAL                   zdu2, ztsolv
  REAL                   ztvd, zscf
  REAL                   zucf, zcr
  REAL                   friv, frih
  REAL, DIMENSION(klon) :: FM, FH                   ! stability functions
  REAL, DIMENSION(klon) :: cdmn, cdhn               ! Drag coefficient in neutral conditions
  REAL zzzcd

  LOGICAL, SAVE :: firstcall = .TRUE.
  !$OMP THREADPRIVATE(firstcall)
  INTEGER, SAVE :: iflag_corr_sta
  !$OMP THREADPRIVATE(iflag_corr_sta)
  INTEGER, SAVE :: iflag_corr_insta
  !$OMP THREADPRIVATE(iflag_corr_insta)

!===================================================================c
! Valeurs numeriques des constantes
!===================================================================c


! Minimum du carre du vent

 CEPDU2 = (0.1)**2

! Louis 1982

  CB=5.0
  CC=5.0
  CD=5.0


! King 2001

  C2=0.25
  C3=0.0625

  
! Louis 1979

  BPRIME=4.7
  B=9.4
  

!MO

  ALPHA=5.0
  

! ================================================================= c
! Tests avant de commencer
! Fuxing WANG, 04/03/2015
! To check if there are negative q1, qsurf values.
!====================================================================c
  ng_q1 = 0     ! Initialization
  ng_qsurf = 0  ! Initialization
  DO i = 1, knon
     IF (q1(i).LT.0.0)     ng_q1 = ng_q1 + 1
     IF (qsurf(i).LT.0.0)  ng_qsurf = ng_qsurf + 1
  ENDDO
  IF (ng_q1.GT.0 .and. prt_level > 5) THEN
      WRITE(lunout,*)" *** Warning: Negative q1(humidity at 1st level) values in cdrag.F90 !"
      WRITE(lunout,*)" The total number of the grids is: ", ng_q1
      WRITE(lunout,*)" The negative q1 is set to zero "
!      abort_message="voir ci-dessus"
!      CALL abort_physic(modname,abort_message,1)
  ENDIF
  IF (ng_qsurf.GT.0 .and. prt_level > 5) THEN
      WRITE(lunout,*)" *** Warning: Negative qsurf(humidity at surface) values in cdrag.F90 !"
      WRITE(lunout,*)" The total number of the grids is: ", ng_qsurf
      WRITE(lunout,*)" The negative qsurf is set to zero "
!      abort_message="voir ci-dessus"
!      CALL abort_physic(modname,abort_message,1)
  ENDIF



!=============================================================================c
! Calcul du cdrag
!=============================================================================c

! On choisit les fonctions de stabilite utilisees au premier appel
!**************************************************************************
  IF (firstcall) THEN
   iflag_corr_sta=2
   iflag_corr_insta=2
 
   CALL getin_p('iflag_corr_sta',iflag_corr_sta)
   CALL getin_p('iflag_corr_insta',iflag_corr_insta)

   firstcall = .FALSE.
 ENDIF

!xxxxxxxxxxxxxxxxxxxxxxx
  DO i = 1, knon        ! Boucle sur l'horizontal
!xxxxxxxxxxxxxxxxxxxxxxx


! calculs preliminaires:
!***********************
     

     zdu2 = MAX(CEPDU2, speed(i)**2)
     pref(i) = EXP(LOG(psol(i)) - zgeop1(i)/(RD*t1(i)* &
                 (1.+ RETV * max(q1(i),0.0))))           ! negative q1 set to zero 
     ztsolv = tsurf(i) * (1.0+RETV*max(qsurf(i),0.0))    ! negative qsurf set to zero
     ztvd = (t1(i)+zgeop1(i)/RCPD/(1.+RVTMP2*max(q1(i),0.0))) &
          *(1.+RETV*max(q1(i),0.0))                      ! negative q1 set to zero
     zri(i) = zgeop1(i)*(ztvd-ztsolv)/(zdu2*ztvd)


! Coefficients CD neutres : k^2/ln(z/z0) et k^2/(ln(z/z0)*ln(z/z0h)):
!********************************************************************

     zzzcd=CKAP/LOG(1.+zgeop1(i)/(RG*z0m(i)))
     cdmn(i) = zzzcd*zzzcd
     cdhn(i) = zzzcd*(CKAP/LOG(1.+zgeop1(i)/(RG*z0h(i))))


! Calcul des fonctions de stabilit?? FMs, FHs, FMi, FHi :
!*******************************************************

!''''''''''''''
! Cas instables
!''''''''''''''

 IF (zri(i) .LT. 0.) THEN    


        SELECT CASE (iflag_corr_insta)
    
        CASE (1) ! Louis 1979 + Mascart 1995

           MU=LOG(MAX(z0m(i)/z0h(i),0.01))
           CMstar=6.8741+2.6933*MU-0.3601*(MU**2)+0.0154*(MU**3)
           PM=0.5233-0.0815*MU+0.0135*(MU**2)-0.001*(MU**3)
           CHstar=3.2165+4.3431*MU+0.536*(MU**2)-0.0781*(MU**3)
           PH=0.5802-0.1571*MU+0.0327*(MU**2)-0.0026*(MU**3)
           CH=CHstar*B*CKAP/LOG(z0m(i)+zgeop1(i)/(RG*z0m(i))) &
            & * CKAPT/LOG(z0h(i)+zgeop1(i)/(RG*z0h(i)))       &
            & * ((zgeop1(i)/(RG*z0h(i)))**PH)
           CM=CMstar*B*CKAP/LOG(z0m(i)+zgeop1(i)/(RG*z0m(i))) &
            & *CKAP/LOG(z0m(i)+zgeop1(i)/(RG*z0m(i)))         &
            & * ((zgeop1(i)/(RG*z0m(i)))**PM)
         



           FM(i)=1.-B*zri(i)/(1.+CM*SQRT(ABS(zri(i))))
           FH(i)=1.-B*zri(i)/(1.+CH*SQRT(ABS(zri(i))))

        CASE (2) ! Louis 1982

           zucf = 1./(1.+3.0*CB*CC*cdmn(i)*SQRT(ABS(zri(i)) &
                *(1.0+zgeop1(i)/(RG*z0m(i)))))
           FM(i) = AMAX1((1.-2.0*CB*zri(i)*zucf),f_ri_cd_min)
           FH(i) = AMAX1((1.-3.0*CB*zri(i)*zucf),f_ri_cd_min)


        CASE (3) ! Laurent Li

           
           FM(i) = MAX(SQRT(1.0-18.0*zri(i)),f_ri_cd_min)
           FH(i) = MAX(SQRT(1.0-18.0*zri(i)),f_ri_cd_min)



         CASE default ! Louis 1982

           zucf = 1./(1.+3.0*CB*CC*cdmn(i)*SQRT(ABS(zri(i)) &
                *(1.0+zgeop1(i)/(RG*z0m(i)))))
           FM(i) = AMAX1((1.-2.0*CB*zri(i)*zucf),f_ri_cd_min)
           FH(i) = AMAX1((1.-3.0*CB*zri(i)*zucf),f_ri_cd_min)


         END SELECT



! Calcul des drags


       cdm(i)=cdmn(i)*FM(i)
       cdh(i)=f_cdrag_ter*cdhn(i)*FH(i)


! Traitement particulier des cas oceaniques
! on applique Miller et al 1992 en l'absence de gustiness 

  IF (nsrf == is_oce) THEN
!        cdh(i)=f_cdrag_oce*cdhn(i)*FH(i)

        IF(iflag_gusts==0) THEN
           zcr = (0.0016/(cdmn(i)*SQRT(zdu2)))*ABS(ztvd-ztsolv)**(1./3.)
           cdh(i) =f_cdrag_oce* cdhn(i)*(1.0+zcr**1.25)**(1./1.25)
        ENDIF


        cdm(i)=MIN(cdm(i),cdmmax)
        cdh(i)=MIN(cdh(i),cdhmax)

  END IF



 ELSE                          

!'''''''''''''''
! Cas stables :
!'''''''''''''''
        zri(i) = MIN(20.,zri(i))

       SELECT CASE (iflag_corr_sta)
    
        CASE (1) ! Louis 1979 + Mascart 1995
   
           FM(i)=MAX(1./((1+BPRIME*zri(i))**2),f_ri_cd_min)
           FH(i)=FM(i)


        CASE (2) ! Louis 1982
           
           zscf = SQRT(1.+CD*ABS(zri(i)))
           FM(i)= AMAX1(1. / (1.+2.*CB*zri(i)/zscf), f_ri_cd_min)
           FH(i)= AMAX1(1./ (1.+3.*CB*zri(i)*zscf), f_ri_cd_min )


        CASE (3) ! Laurent Li
   
           FM(i)=MAX(1.0 / (1.0+10.0*zri(i)*(1+8.0*zri(i))),f_ri_cd_min)
           FH(i)=FM(i)


        CASE (4)  ! King 2001
          
           if (zri(i) .LT. C2/2.) then
           FM(i)=MAX((1.-zri(i)/C2)**2,f_ri_cd_min)
           FH(i)=  FM(i)


           else
           FM(i)=MAX(C3*((C2/zri(i))**2),f_ri_cd_min)
           FH(i)= FM(i)
           endif


        CASE (5) ! MO
   
          if (zri(i) .LT. 1./alpha) then

           FM(i)=MAX((1.-alpha*zri(i))**2,f_ri_cd_min)
           FH(i)=FM(i)

           else 


           FM(i)=MAX(1E-7,f_ri_cd_min)
           FH(i)=FM(i)

          endif





        CASE default ! Louis 1982

           zscf = SQRT(1.+CD*ABS(zri(i)))
           FM(i)= AMAX1(1. / (1.+2.*CB*zri(i)/zscf), f_ri_cd_min)
           FH(i)= AMAX1(1./ (1.+3.*CB*zri(i)*zscf), f_ri_cd_min )



   END SELECT

! Calcul des drags


       cdm(i)=cdmn(i)*FM(i)
       cdh(i)=f_cdrag_ter*cdhn(i)*FH(i)

       IF(nsrf.EQ.is_oce) THEN

        cdh(i)=f_cdrag_oce*cdhn(i)*FH(i)
        cdm(i)=MIN(cdm(i),cdmmax)
        cdh(i)=MIN(cdh(i),cdhmax)

      ENDIF


 ENDIF

!xxxxxxxxxxx
  END DO   !  Fin de la boucle sur l'horizontal
!xxxxxxxxxxx
! ================================================================= c
     
END SUBROUTINE cdrag

!
 SUBROUTINE cdragn_ri1(knon,  nsrf,   &
     speed, t1,    q1,    zgeop1, &
     psol,  tsurf, qsurf, z0m, z0h,  &
     ri1, iri1, &
     cdm,  cdh,  zri,   pref)

  USE dimphy
  USE indice_sol_mod
  USE print_control_mod, ONLY: lunout, prt_level
  USE ioipsl_getin_p_mod, ONLY : getin_p

  IMPLICIT NONE
! ================================================================= c
!
! Objet : calcul des cdrags pour le moment (pcfm) et 
!         les flux de chaleur sensible et latente (pcfh) d'apr??s
!         Louis 1982, Louis 1979, King et al 2001
!         ou Zilitinkevich et al 2002  pour les cas stables, Louis 1979
!         et 1982 pour les cas instables 
!
! Modified history: 
!  writting on the 20/05/2016
!  modified on the 13/12/2016 to be adapted to LMDZ6
!
! References:
!   Louis, J. F., 1979: A parametric model of vertical eddy fluxes in the
!     atmosphere. Boundary-Layer Meteorology. 01/1979; 17(2):187-202. 
!   Louis, J. F., Tiedtke, M. and Geleyn, J. F., 1982: `A short history of the
!     operational PBL parametrization at ECMWF'. Workshop on boundary layer
!     parametrization, November 1981, ECMWF, Reading, England. 
!     Page: 19. Equations in Table 1.
!   Mascart P, Noilhan J, Giordani H 1995.A MODIFIED PARAMETERIZATION OF FLUX-PROFILE RELATIONSHIPS 
!    IN THE SURFACE LAYER USING DIFFERENT ROUGHNESS LENGTH VALUES FOR HEAT AND MOMENTUM
!    Boundary-Layer Meteorology 72: 331-344
!   Anton Beljaars. May 1992. The parametrization of the planetary boundary layer.  
!     European Centre for Medium-Range Weather Forecasts.
!     Equations: 110-113. Page 40.
!   Miller,M.J., A.C.M.Beljaars, T.N.Palmer. 1992. The sensitivity of the ECMWF
!     model to the parameterization of evaporation from the tropical oceans. J.
!     Climate, 5:418-434.
!   King J.C, Connolley, W.M ad Derbyshire S.H. 2001, Sensitivity of Modelled Antarctic climate
!   to surface and boundary-layer flux parametrizations
!   QJRMS, 127, pp 779-794
!
! ================================================================= c
! ================================================================= c
! On choisit le couple de fonctions de correction avec deux flags:
! Un pour les cas instables, un autre pour les cas stables
!
! iflag_corr_insta:
!                   1: Louis 1979 avec les modifications de Mascart 1995 (z0/= z0h)
!                   2: Louis 1982 
!                   3: Laurent Li
!
! iflag_corr_sta:
!                   1: Louis 1979 avec les modifications de Mascart 1995 (z0/= z0h)
!                   2: Louis 1982 
!                   3: Laurent Li
!                   4: King  2001 (SHARP)
!                   5: MO 1st order theory (allow collapse of turbulence)
!            
!
!*****************************************************************
! Parametres d'entree
!***************************************************************** 
  
  INTEGER, INTENT(IN)                      :: knon, nsrf ! nombre de points de grille sur l'horizontal + type de surface
  REAL, DIMENSION(klon), INTENT(IN)        :: speed ! module du vent au 1er niveau du modele
  REAL, DIMENSION(klon), INTENT(IN)        :: zgeop1! geopotentiel au 1er niveau du modele
  REAL, DIMENSION(klon), INTENT(IN)        :: tsurf ! Surface temperature (K)
  REAL, DIMENSION(klon), INTENT(IN)        :: qsurf ! Surface humidity (Kg/Kg)
  REAL, DIMENSION(klon), INTENT(IN)        :: z0m, z0h ! Rugosity at surface (m) 
  REAL, DIMENSION(klon), INTENT(IN)        :: ri1   ! Richardson 1ere couche 
  INTEGER, INTENT(IN)                      :: iri1   !
  REAL, DIMENSION(klon), INTENT(IN)        :: t1  ! Temperature au premier niveau (K) 
  REAL, DIMENSION(klon), INTENT(IN)        :: q1  ! humidite specifique au premier niveau (kg/kg)
  REAL, DIMENSION(klon), INTENT(IN)        :: psol ! pression au sol



! Parametres de sortie
!******************************************************************
  REAL, DIMENSION(klon), INTENT(OUT)       :: cdm  ! Drag coefficient for heat flux
  REAL, DIMENSION(klon), INTENT(OUT)       :: cdh  ! Drag coefficient for momentum
  REAL, DIMENSION(klon), INTENT(OUT)       :: zri   ! Richardson number
  REAL, DIMENSION(klon), INTENT(OUT)       :: pref  ! Pression au niveau zgeop/RG
! Variables Locales
!******************************************************************
 

  INCLUDE "YOMCST.h"
  INCLUDE "YOETHF.h"
  INCLUDE "clesphys.h"


  REAL, PARAMETER       :: CKAP=0.40, CKAPT=0.42
  REAL                   CEPDU2
  REAL                   ALPHA
  REAL                   CB,CC,CD,C2,C3
  REAL                   MU, CM, CH, B, CMstar, CHstar
  REAL                   PM, PH, BPRIME
  REAL                   C
  INTEGER                ng_q1                      ! Number of grids that q1 < 0.0
  INTEGER                ng_qsurf                   ! Number of grids that qsurf < 0.0
  INTEGER                i
  REAL                   zdu2, ztsolv
  REAL                   ztvd, zscf
  REAL                   zucf, zcr
  REAL                   friv, frih
  REAL, DIMENSION(klon) :: FM, FH                   ! stability functions
  REAL, DIMENSION(klon) :: cdmn, cdhn               ! Drag coefficient in neutral conditions
  REAL zzzcd

  LOGICAL, SAVE :: firstcall = .TRUE.
  !$OMP THREADPRIVATE(firstcall)
  INTEGER, SAVE :: iflag_corr_sta
  !$OMP THREADPRIVATE(iflag_corr_sta)
  INTEGER, SAVE :: iflag_corr_insta
  !$OMP THREADPRIVATE(iflag_corr_insta)

!===================================================================c
! Valeurs numeriques des constantes
!===================================================================c


! Minimum du carre du vent

 CEPDU2 = (0.1)**2
! Louis 1982

  CB=5.0
  CC=5.0
  CD=5.0


! King 2001

  C2=0.25
  C3=0.0625

  
! Louis 1979

  BPRIME=4.7
  B=9.4
  

!MO

  ALPHA=5.0
  

! ================================================================= c
! Tests avant de commencer
! Fuxing WANG, 04/03/2015
! To check if there are negative q1, qsurf values.
!====================================================================c
  ng_q1 = 0     ! Initialization
  ng_qsurf = 0  ! Initialization
  DO i = 1, knon
     IF (q1(i).LT.0.0)     ng_q1 = ng_q1 + 1
     IF (qsurf(i).LT.0.0)  ng_qsurf = ng_qsurf + 1
  ENDDO
  IF (ng_q1.GT.0 .and. prt_level > 5) THEN
      WRITE(lunout,*)" *** Warning: Negative q1(humidity at 1st level) values in cdrag.F90 !"
      WRITE(lunout,*)" The total number of the grids is: ", ng_q1
      WRITE(lunout,*)" The negative q1 is set to zero "
!      abort_message="voir ci-dessus"
!      CALL abort_physic(modname,abort_message,1)
  ENDIF
  IF (ng_qsurf.GT.0 .and. prt_level > 5) THEN
      WRITE(lunout,*)" *** Warning: Negative qsurf(humidity at surface) values in cdrag.F90 !"
      WRITE(lunout,*)" The total number of the grids is: ", ng_qsurf
      WRITE(lunout,*)" The negative qsurf is set to zero "
!      abort_message="voir ci-dessus"
!      CALL abort_physic(modname,abort_message,1)
  ENDIF



!=============================================================================c
! Calcul du cdrag
!=============================================================================c

! On choisit les fonctions de stabilite utilisees au premier appel
!**************************************************************************
  IF (firstcall) THEN
   iflag_corr_sta=2
   iflag_corr_insta=2
 
   CALL getin_p('iflag_corr_sta',iflag_corr_sta)
   CALL getin_p('iflag_corr_insta',iflag_corr_insta)

   firstcall = .FALSE.
 ENDIF

!xxxxxxxxxxxxxxxxxxxxxxx
  DO i = 1, knon        ! Boucle sur l'horizontal
!xxxxxxxxxxxxxxxxxxxxxxx


! calculs preliminaires:
!***********************
     

     zdu2 = MAX(CEPDU2, speed(i)**2)
     pref(i) = EXP(LOG(psol(i)) - zgeop1(i)/(RD*t1(i)* &
                 (1.+ RETV * max(q1(i),0.0))))           ! negative q1 set to zero 
     ztsolv = tsurf(i) * (1.0+RETV*max(qsurf(i),0.0))    ! negative qsurf set to zero
     ztvd = (t1(i)+zgeop1(i)/RCPD/(1.+RVTMP2*max(q1(i),0.0))) &
          *(1.+RETV*max(q1(i),0.0))                      ! negative q1 set to zero
     zri(i) = zgeop1(i)*(ztvd-ztsolv)/(zdu2*ztvd)

     IF (iri1.EQ.1) THEN
       zri(i) = ri1(i)
     ENDIF

! Coefficients CD neutres : k^2/ln(z/z0) et k^2/(ln(z/z0)*ln(z/z0h)):
!********************************************************************

     zzzcd=CKAP/LOG(1.+zgeop1(i)/(RG*z0m(i)))
     cdmn(i) = zzzcd*zzzcd
     cdhn(i) = zzzcd*(CKAP/LOG(1.+zgeop1(i)/(RG*z0h(i))))


! Calcul des fonctions de stabilit?? FMs, FHs, FMi, FHi :
!*******************************************************

!''''''''''''''
! Cas instables
!''''''''''''''

 IF (zri(i) .LT. 0.) THEN    


        SELECT CASE (iflag_corr_insta)
    
        CASE (1) ! Louis 1979 + Mascart 1995

           MU=LOG(MAX(z0m(i)/z0h(i),0.01))
           CMstar=6.8741+2.6933*MU-0.3601*(MU**2)+0.0154*(MU**3)
           PM=0.5233-0.0815*MU+0.0135*(MU**2)-0.001*(MU**3)
           CHstar=3.2165+4.3431*MU+0.536*(MU**2)-0.0781*(MU**3)
           PH=0.5802-0.1571*MU+0.0327*(MU**2)-0.0026*(MU**3)
           CH=CHstar*B*CKAP/LOG(z0m(i)+zgeop1(i)/(RG*z0m(i))) &
            & * CKAPT/LOG(z0h(i)+zgeop1(i)/(RG*z0h(i)))       &
            & * ((zgeop1(i)/(RG*z0h(i)))**PH)
           CM=CMstar*B*CKAP/LOG(z0m(i)+zgeop1(i)/(RG*z0m(i))) &
            & *CKAP/LOG(z0m(i)+zgeop1(i)/(RG*z0m(i)))         &
            & * ((zgeop1(i)/(RG*z0m(i)))**PM)
         



           FM(i)=1.-B*zri(i)/(1.+CM*SQRT(ABS(zri(i))))
           FH(i)=1.-B*zri(i)/(1.+CH*SQRT(ABS(zri(i))))

        CASE (2) ! Louis 1982

           zucf = 1./(1.+3.0*CB*CC*cdmn(i)*SQRT(ABS(zri(i)) &
                *(1.0+zgeop1(i)/(RG*z0m(i)))))
           FM(i) = AMAX1((1.-2.0*CB*zri(i)*zucf),f_ri_cd_min)
           FH(i) = AMAX1((1.-3.0*CB*zri(i)*zucf),f_ri_cd_min)


        CASE (3) ! Laurent Li

           
           FM(i) = MAX(SQRT(1.0-18.0*zri(i)),f_ri_cd_min)
           FH(i) = MAX(SQRT(1.0-18.0*zri(i)),f_ri_cd_min)



         CASE default ! Louis 1982

           zucf = 1./(1.+3.0*CB*CC*cdmn(i)*SQRT(ABS(zri(i)) &
                *(1.0+zgeop1(i)/(RG*z0m(i)))))
           FM(i) = AMAX1((1.-2.0*CB*zri(i)*zucf),f_ri_cd_min)
           FH(i) = AMAX1((1.-3.0*CB*zri(i)*zucf),f_ri_cd_min)


         END SELECT



! Calcul des drags


       cdm(i)=cdmn(i)*FM(i)
       cdh(i)=f_cdrag_ter*cdhn(i)*FH(i)


! Traitement particulier des cas oceaniques
! on applique Miller et al 1992 en l'absence de gustiness 

  IF (nsrf == is_oce) THEN
!        cdh(i)=f_cdrag_oce*cdhn(i)*FH(i)

        IF(iflag_gusts==0) THEN
           zcr = (0.0016/(cdmn(i)*SQRT(zdu2)))*ABS(ztvd-ztsolv)**(1./3.)
           cdh(i) =f_cdrag_oce* cdhn(i)*(1.0+zcr**1.25)**(1./1.25)
        ENDIF


        cdm(i)=MIN(cdm(i),cdmmax)
        cdh(i)=MIN(cdh(i),cdhmax)

  END IF



 ELSE                          

!'''''''''''''''
! Cas stables :
!'''''''''''''''
        zri(i) = MIN(20.,zri(i))

       SELECT CASE (iflag_corr_sta)
    
        CASE (1) ! Louis 1979 + Mascart 1995
   
           FM(i)=MAX(1./((1+BPRIME*zri(i))**2),f_ri_cd_min)
           FH(i)=FM(i)


        CASE (2) ! Louis 1982
           
           zscf = SQRT(1.+CD*ABS(zri(i)))
           FM(i)= AMAX1(1. / (1.+2.*CB*zri(i)/zscf), f_ri_cd_min)
           FH(i)= AMAX1(1./ (1.+3.*CB*zri(i)*zscf), f_ri_cd_min )


        CASE (3) ! Laurent Li
   
           FM(i)=MAX(1.0 / (1.0+10.0*zri(i)*(1+8.0*zri(i))),f_ri_cd_min)
           FH(i)=FM(i)


        CASE (4)  ! King 2001
          
           if (zri(i) .LT. C2/2.) then
           FM(i)=MAX((1.-zri(i)/C2)**2,f_ri_cd_min)
           FH(i)=  FM(i)


           else
           FM(i)=MAX(C3*((C2/zri(i))**2),f_ri_cd_min)
           FH(i)= FM(i)
           endif


        CASE (5) ! MO
   
          if (zri(i) .LT. 1./alpha) then

           FM(i)=MAX((1.-alpha*zri(i))**2,f_ri_cd_min)
           FH(i)=FM(i)

           else 


           FM(i)=MAX(1E-7,f_ri_cd_min)
           FH(i)=FM(i)

          endif





        CASE default ! Louis 1982

           zscf = SQRT(1.+CD*ABS(zri(i)))
           FM(i)= AMAX1(1. / (1.+2.*CB*zri(i)/zscf), f_ri_cd_min)
           FH(i)= AMAX1(1./ (1.+3.*CB*zri(i)*zscf), f_ri_cd_min )



   END SELECT

! Calcul des drags


       cdm(i)=cdmn(i)*FM(i)
       cdh(i)=f_cdrag_ter*cdhn(i)*FH(i)

       IF(nsrf.EQ.is_oce) THEN

        cdh(i)=f_cdrag_oce*cdhn(i)*FH(i)
        cdm(i)=MIN(cdm(i),cdmmax)
        cdh(i)=MIN(cdh(i),cdhmax)

      ENDIF


 ENDIF

!xxxxxxxxxxx
  END DO   !  Fin de la boucle sur l'horizontal
!xxxxxxxxxxx
! ================================================================= c
     
END SUBROUTINE cdragn_ri1

END MODULE cdrag_mod
