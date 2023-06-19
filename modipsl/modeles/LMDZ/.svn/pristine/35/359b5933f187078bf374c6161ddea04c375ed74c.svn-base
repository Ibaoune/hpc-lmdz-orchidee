!
! $Id$
!
!
MODULE add_phys_tend_mod

  IMPLICIT NONE  
  ! flag to compute diagnostics to check energy conservation. If  fl_ebil==0, no check
  INTEGER, SAVE ::   fl_ebil
!$OMP THREADPRIVATE(fl_ebil)
  ! flag to include modifcations to ensure energy conservation. If  fl_cor_ebil==0, no corrections
  ! Note that with time, all these modifications should be included by default
  INTEGER, SAVE ::   fl_cor_ebil
!$OMP THREADPRIVATE(fl_cor_ebil)

CONTAINS

SUBROUTINE add_pbl_tend(zdu, zdv, zdt, zdq, zdql, zdqi, paprs, text,abortphy,flag_inhib_tend, itap)
  ! ======================================================================
  ! Ajoute les tendances de couche limite, soit determinees par la
  ! parametrisation
  ! physique, soit forcees,  aux variables d etat de la dynamique t_seri,
  ! q_seri ...
  ! ======================================================================


  ! ======================================================================
  ! Declarations
  ! ======================================================================

  USE dimphy, ONLY: klon, klev
!  USE dimphy
  USE phys_local_var_mod
  USE phys_state_var_mod
  USE mod_grid_phy_lmdz, ONLY: nbp_lev
  IMPLICIT NONE
  REAL,SAVE,ALLOCATABLE :: hthturb_gcssold(:)
  REAL,SAVE,ALLOCATABLE :: hqturb_gcssold(:)
!$OMP THREADPRIVATE(hthturb_gcssold,hqturb_gcssold)
  REAL,SAVE :: dtime_frcg
  LOGICAL,SAVE :: turb_fcg_gcssold
  LOGICAL,SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall,turb_fcg_gcssold,dtime_frcg)
  INTEGER abortphy
!  COMMON /turb_forcing/dtime_frcg, hthturb_gcssold, hqturb_gcssold, &
!    turb_fcg_gcssold

  ! Arguments :
  ! ------------
  REAL zdu(klon, klev), zdv(klon, klev)
  REAL zdt(klon, klev), zdq(klon, klev), zdql(klon, klev), zdqi(klon, klev)
  CHARACTER *(*) text
  REAL paprs(klon,klev+1)
  INTEGER flag_inhib_tend ! if flag_inhib_tend != 0, tendencies are not added
  INTEGER itap ! time step number

  ! Local :
  ! --------
  REAL zzdt(klon, klev), zzdq(klon, klev)
  INTEGER i, k

  IF (firstcall) THEN
    ALLOCATE(hthturb_gcssold(nbp_lev))
    ALLOCATE(hqturb_gcssold(nbp_lev))
    firstcall=.false.
  ENDIF

  IF (turb_fcg_gcssold) THEN
    DO k = 1, klev
      DO i = 1, klon
        zzdt(i, k) = hthturb_gcssold(k)*dtime_frcg
        zzdq(i, k) = hqturb_gcssold(k)*dtime_frcg
      END DO
    END DO
    PRINT *, ' add_pbl_tend, dtime_frcg ', dtime_frcg
    PRINT *, ' add_pbl_tend, zzdt ', zzdt
    PRINT *, ' add_pbl_tend, zzdq ', zzdq
    CALL add_phys_tend(zdu, zdv, zzdt, zzdq, zdql, zdqi, paprs, text,abortphy,flag_inhib_tend, itap, 0)
  ELSE
    CALL add_phys_tend(zdu, zdv, zdt, zdq, zdql, zdqi, paprs, text,abortphy,flag_inhib_tend, itap, 0)
  END IF


  RETURN
END SUBROUTINE add_pbl_tend
!
! $Id: add_phys_tend.F90 2611 2016-08-03 15:41:26Z jyg $
!
SUBROUTINE add_phys_tend (zdu,zdv,zdt,zdq,zdql,zdqi,paprs,text, &
                          abortphy,flag_inhib_tend, itap, diag_mode)
!======================================================================
! Ajoute les tendances des variables physiques aux variables 
! d'etat de la dynamique t_seri, q_seri ...
! On en profite pour faire des tests sur les tendances en question.
!======================================================================


!======================================================================
! Declarations
!======================================================================

USE dimphy, ONLY: klon, klev
USE phys_state_var_mod, ONLY : phys_tstep
USE phys_local_var_mod, ONLY: u_seri, v_seri, ql_seri, qs_seri, q_seri, t_seri
USE phys_state_var_mod, ONLY: ftsol
USE geometry_mod, ONLY: longitude_deg, latitude_deg
USE print_control_mod, ONLY: prt_level
USE cmp_seri_mod
USE phys_output_var_mod, ONLY : d_qw_col, d_ql_col, d_qs_col, d_qt_col, d_ek_col, d_h_dair_col &
  &           , d_h_qw_col, d_h_ql_col, d_h_qs_col, d_h_col
IMPLICIT none
  include "YOMCST.h"
  include "clesphys.h"

! Arguments :
!------------
REAL, DIMENSION(klon,klev),     INTENT(IN)    :: zdu, zdv
REAL, DIMENSION(klon,klev),     INTENT(IN)    :: zdt, zdql, zdqi
REAL, DIMENSION(klon,klev+1),   INTENT(IN)    :: paprs
CHARACTER*(*),                  INTENT(IN)    :: text
INTEGER,                        INTENT(IN)    :: abortphy
INTEGER,                        INTENT(IN)    :: flag_inhib_tend ! if not 0, tendencies are not added
INTEGER,                        INTENT(IN)    :: itap            ! time step number
INTEGER,                        INTENT(IN)    :: diag_mode       ! 0 -> normal effective mode
                                                                 ! 1 -> only conservation stats are computed
!
REAL, DIMENSION(klon,klev),     INTENT(INOUT) :: zdq

! Local :
!--------
REAL zt,zq
REAL zq_int, zqp_int, zq_new

REAL zqp(klev)

! Save variables, used in diagnostic mode (diag_mode=1).
REAL, DIMENSION(klon,klev)   :: sav_u_seri, sav_v_seri
REAL, DIMENSION(klon,klev)   :: sav_ql_seri, sav_qs_seri, sav_q_seri
REAL, DIMENSION(klon,klev)   :: sav_t_seri
REAL, DIMENSION(klon,klev)   :: sav_zdq
!
INTEGER i, k,j, n
INTEGER jadrs(klon*klev), jbad
INTEGER jqadrs(klon*klev), jqbad
INTEGER kadrs(klon*klev)
INTEGER kqadrs(klon*klev)

LOGICAL done(klon)

integer debug_level
logical, save :: first=.true.
!$OMP THREADPRIVATE(first)
!
!======================================================================
! Variables for energy conservation tests
!======================================================================
!

! zh_col-------  total enthalpy of vertical air column
! (air with watter vapour, liquid and solid) (J/m2)
! zh_dair_col--- total enthalpy of dry air (J/m2)
! zh_qw_col----  total enthalpy of watter vapour (J/m2)
! zh_ql_col----  total enthalpy of liquid watter (J/m2)
! zh_qs_col----  total enthalpy of solid watter  (J/m2)
! zqw_col------  total mass of watter vapour (kg/m2)
! zql_col------  total mass of liquid watter (kg/m2)
! zqs_col------  total mass of solid watter (kg/m2)
! zek_col------  total kinetic energy (kg/m2)
!
REAL zairm(klon, klev) ! layer air mass (kg/m2)
REAL zqw_col(klon,2)
REAL zql_col(klon,2)
REAL zqs_col(klon,2)
REAL zek_col(klon,2)
REAL zh_dair_col(klon,2)
REAL zh_qw_col(klon,2), zh_ql_col(klon,2), zh_qs_col(klon,2)
REAL zh_col(klon,2)

REAL zcpvap, zcwat, zcice

!======================================================================
! Initialisations

     IF (prt_level >= 5) then
        write (*,*) "In add_phys_tend, after ",text
        call flush
     end if

     ! if flag_inhib_tend != 0, tendencies are not added
     IF (flag_inhib_tend /= 0) then
        ! If requiered, diagnostics are shown
        IF (flag_inhib_tend > 0) then
           ! print some diagnostics if xxx_seri have changed
           call cmp_seri(flag_inhib_tend,text)
        END IF
        RETURN ! on n ajoute pas les tendance
     END IF

     IF (abortphy==1) RETURN ! on n ajoute pas les tendance si le modele
                              ! a deja plante.

     debug_level=10
     if (first) then
        print *,"TestJLD rcpv, rcw, rcs",rcpv, rcw, rcs
        first=.false.
     endif
!
!  print *,'add_phys_tend: paprs ',paprs
! When in diagnostic mode, save initial values of out variables
  IF (diag_mode == 1) THEN
    sav_u_seri(:,:)  = u_seri(:,:)
    sav_v_seri(:,:)  = v_seri(:,:)
    sav_ql_seri(:,:) = ql_seri(:,:)
    sav_qs_seri(:,:) = qs_seri(:,:)
    sav_q_seri(:,:)  = q_seri(:,:)
    sav_t_seri(:,:)  = t_seri(:,:)
    sav_zdq(:,:)     = zdq(:,:)   
  ENDIF ! (diag_mode == 1)
!======================================================================
! Diagnostics for energy conservation tests
!======================================================================
  DO k = 1, klev
    ! layer air mass
    zairm(:, k) = (paprs(:,k)-paprs(:,k+1))/rg
  END DO

  if (fl_ebil .GT. 0) then
    ! ------------------------------------------------
    ! Compute vertical sum for each atmospheric column
    ! ------------------------------------------------
    n=1   ! begining of time step

    zcpvap = rcpv
    zcwat = rcw
    zcice = rcs

    CALL integr_v(klon, klev, zcpvap, &
                  t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, zairm, &
                    zqw_col(:,n), zql_col(:,n), zqs_col(:,n), zek_col(:,n), zh_dair_col(:,n), &
                    zh_qw_col(:,n), zh_ql_col(:,n), zh_qs_col(:,n), zh_col(:,n))

  end if ! end if (fl_ebil .GT. 0)

!======================================================================
! Ajout des tendances sur le vent et l'eau liquide
!======================================================================

     u_seri(:,:)=u_seri(:,:)+zdu(:,:)
     v_seri(:,:)=v_seri(:,:)+zdv(:,:)
     ql_seri(:,:)=ql_seri(:,:)+zdql(:,:)
     qs_seri(:,:)=qs_seri(:,:)+zdqi(:,:)

!======================================================================
! On ajoute les tendances de la temperature et de la vapeur d'eau
! en verifiant que ca ne part pas dans les choux
!======================================================================

      jbad=0
      jqbad=0
      DO k = 1, klev
         DO i = 1, klon
            zt=t_seri(i,k)+zdt(i,k)
            zq=q_seri(i,k)+zdq(i,k)
            IF ( zt>370. .or. zt<130. .or. abs(zdt(i,k))>50. ) then
            jbad = jbad + 1
            jadrs(jbad) = i
            kadrs(jbad) = k
            ENDIF
            IF ( zq<0. .or. zq>0.1 .or. abs(zdq(i,k))>1.e-2 ) then
            jqbad = jqbad + 1
            jqadrs(jqbad) = i
            kqadrs(jqbad) = k
            ENDIF
            t_seri(i,k)=zt
            q_seri(i,k)=zq
         ENDDO
      ENDDO

!=====================================================================================
! Impression et stop en cas de probleme important
!=====================================================================================

IF (jbad .GT. 0) THEN
      DO j = 1, jbad
         i=jadrs(j)
         if(prt_level.ge.debug_level) THEN
          print*,'PLANTAGE POUR LE POINT i lon lat =',&
                 i,longitude_deg(i),latitude_deg(i),text
          print*,'l    T     dT       Q     dQ    '
          DO k = 1, klev
             write(*,'(i3,2f14.4,2e14.2)') k,t_seri(i,k),zdt(i,k),q_seri(i,k),zdq(i,k)
          ENDDO
          call print_debug_phys(i,debug_level,text)
         endif
      ENDDO
ENDIF
!
!=====================================================================================
! Impression, warning et correction en cas de probleme moins important
!=====================================================================================
IF (jqbad .GT. 0) THEN
      done(:) = .false.                         !jyg
      DO j = 1, jqbad
        i=jqadrs(j)
          if(prt_level.ge.debug_level) THEN
           print*,'WARNING  : EAU POUR LE POINT i lon lat =',&
                  i,longitude_deg(i),latitude_deg(i),text
           print*,'l    T     dT       Q     dQ    '
           DO k = 1, klev
              write(*,'(i3,2f14.4,2e14.2)') k,t_seri(i,k),zdt(i,k),q_seri(i,k),zdq(i,k)
           ENDDO
          endif
          IF (ok_conserv_q) THEN
!jyg<20140228 Corrections pour conservation de l'eau
            IF (.NOT.done(i)) THEN                  !jyg
              DO k = 1, klev
                zqp(k) = max(q_seri(i,k),1.e-15)
              ENDDO
              zq_int  = 0.
              zqp_int = 0.
              DO k = 1, klev
                zq_int  = zq_int  + q_seri(i,k)*(paprs(i,k)-paprs(i,k+1))/Rg
                zqp_int = zqp_int + zqp(k)     *(paprs(i,k)-paprs(i,k+1))/Rg
              ENDDO
              if(prt_level.ge.debug_level) THEN
               print*,' cas q_seri<1.e-15 i k zq_int zqp_int zq_int/zqp_int :', &
                                    i, kqadrs(j), zq_int, zqp_int, zq_int/zqp_int
              endif
              DO k = 1, klev
                zq_new = zqp(k)*zq_int/zqp_int
                zdq(i,k) = zdq(i,k) + zq_new - q_seri(i,k)
                q_seri(i,k) = zq_new
              ENDDO
              done(i) = .true.
            ENDIF !(.NOT.done(i))
          ELSE 
!jyg>
            DO k = 1, klev
              zq=q_seri(i,k)+zdq(i,k)
              if (zq.lt.1.e-15) then
                 if (q_seri(i,k).lt.1.e-15) then
                  if(prt_level.ge.debug_level) THEN
                   print*,' cas q_seri<1.e-15 i k q_seri zq zdq :',i,k,q_seri(i,k),zq,zdq(i,k)
                  endif
                  q_seri(i,k)=1.e-15
                  zdq(i,k)=(1.e-15-q_seri(i,k))
                 endif
              endif
!              zq=q_seri(i,k)+zdq(i,k)
!              if (zq.lt.1.e-15) then
!                 zdq(i,k)=(1.e-15-q_seri(i,k))
!              endif
            ENDDO
!jyg<20140228
          ENDIF   !  (ok_conserv_q)
!jyg>
      ENDDO ! j = 1, jqbad
ENDIF
!

!IM ajout memes tests pour reverifier les jbad, jqbad beg
      jbad=0
      jqbad=0
      DO k = 1, klev
         DO i = 1, klon
            IF ( t_seri(i,k)>370. .or. t_seri(i,k)<130. .or. abs(zdt(i,k))>50. ) then
            jbad = jbad + 1
            jadrs(jbad) = i
!            if(prt_level.ge.debug_level) THEN
!             print*,'cas2 i k t_seri zdt',i,k,t_seri(i,k),zdt(i,k)
!            endif
            ENDIF
            IF ( q_seri(i,k)<0. .or. q_seri(i,k)>0.1 .or. abs(zdq(i,k))>1.e-2 ) then
            jqbad = jqbad + 1
            jqadrs(jqbad) = i
            kqadrs(jqbad) = k
!            if(prt_level.ge.debug_level) THEN
!             print*,'cas2 i k q_seri zdq',i,k,q_seri(i,k),zdq(i,k)
!            endif
            ENDIF
         ENDDO
      ENDDO
IF (jbad .GT. 0) THEN
      DO j = 1, jbad
         i=jadrs(j)
         k=kadrs(j)
         if(prt_level.ge.debug_level) THEN
          print*,'PLANTAGE2 POUR LE POINT i itap lon lat txt jbad zdt t',&
                 i,itap,longitude_deg(i),latitude_deg(i),text,jbad, &
       &        zdt(i,k),t_seri(i,k)-zdt(i,k)
!!!       if(prt_level.ge.10.and.itap.GE.229.and.i.EQ.3027) THEN 
          print*,'l    T     dT       Q     dQ    '
          DO k = 1, klev
             write(*,'(i3,2f14.4,2e14.2)') k,t_seri(i,k),zdt(i,k),q_seri(i,k),zdq(i,k)
          ENDDO
          call print_debug_phys(i,debug_level,text)
         endif
      ENDDO
ENDIF 
!
IF (jqbad .GT. 0) THEN
      DO j = 1, jqbad
         i=jqadrs(j)
         k=kqadrs(j)
         if(prt_level.ge.debug_level) THEN
          print*,'WARNING  : EAU2 POUR LE POINT i itap lon lat txt jqbad zdq q zdql ql',&
                 i,itap,longitude_deg(i),latitude_deg(i),text,jqbad,&
       &        zdq(i,k), q_seri(i,k)-zdq(i,k), zdql(i,k), ql_seri(i,k)-zdql(i,k)
!!!       if(prt_level.ge.10.and.itap.GE.229.and.i.EQ.3027) THEN 
          print*,'l    T     dT       Q     dQ    '
          DO k = 1, klev
            write(*,'(i3,2f14.4,2e14.2)') k,t_seri(i,k),zdt(i,k),q_seri(i,k),zdq(i,k)
          ENDDO
          call print_debug_phys(i,debug_level,text)
         endif
      ENDDO
ENDIF

!======================================================================
! Contrôle des min/max pour arrêt du modèle
! Si le modele est en mode abortphy, on retire les tendances qu'on
! vient d'ajouter. Pas exactement parce qu'on ne tient pas compte des
! seuils.
!======================================================================

      CALL hgardfou(t_seri,ftsol,text,abortphy) 
      IF (abortphy==1) THEN
        Print*,'ERROR ABORT hgardfou dans ',text
! JLD pourquoi on ne modifie pas de meme t_seri et q_seri ?
        u_seri(:,:)=u_seri(:,:)-zdu(:,:)
        v_seri(:,:)=v_seri(:,:)-zdv(:,:)
        ql_seri(:,:)=ql_seri(:,:)-zdql(:,:)
        qs_seri(:,:)=qs_seri(:,:)-zdqi(:,:)
      ENDIF

!======================================================================
! Diagnostics for energy conservation tests
!======================================================================

  if (fl_ebil .GT. 0) then
  
    ! ------------------------------------------------
    ! Compute vertical sum for each atmospheric column
    ! ------------------------------------------------
    n=2   ! end of time step

    CALL integr_v(klon, klev, zcpvap, &
                  t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, zairm, &
                    zqw_col(:,n), zql_col(:,n), zqs_col(:,n), zek_col(:,n), zh_dair_col(:,n), &
                    zh_qw_col(:,n), zh_ql_col(:,n), zh_qs_col(:,n), zh_col(:,n))

    ! ------------------------------------------------
    ! Compute the changes by unit of time
    ! ------------------------------------------------

    d_qw_col(:) = (zqw_col(:,2)-zqw_col(:,1))/phys_tstep
    d_ql_col(:) = (zql_col(:,2)-zql_col(:,1))/phys_tstep
    d_qs_col(:) = (zqs_col(:,2)-zqs_col(:,1))/phys_tstep
    d_qt_col(:) = d_qw_col(:) + d_ql_col(:) + d_qs_col(:)

    d_ek_col(:) = (zek_col(:,2)-zek_col(:,1))/phys_tstep

    d_h_dair_col(:) = (zh_dair_col(:,2)-zh_dair_col(:,1))/phys_tstep
    d_h_qw_col(:) = (zh_qw_col(:,2)-zh_qw_col(:,1))/phys_tstep
    d_h_ql_col(:) = (zh_ql_col(:,2)-zh_ql_col(:,1))/phys_tstep
    d_h_qs_col(:) = (zh_qs_col(:,2)-zh_qs_col(:,1))/phys_tstep

    d_h_col = (zh_col(:,2)-zh_col(:,1))/phys_tstep 

  end if ! end if (fl_ebil .GT. 0)
!
! When in diagnostic mode, restore "out" variables to initial values.
  IF (diag_mode == 1) THEN
    u_seri(:,:)  = sav_u_seri(:,:)
    v_seri(:,:)  = sav_v_seri(:,:)
    ql_seri(:,:) = sav_ql_seri(:,:)
    qs_seri(:,:) = sav_qs_seri(:,:)
    q_seri(:,:)  = sav_q_seri(:,:)
    t_seri(:,:)  = sav_t_seri(:,:)
    zdq(:,:)     = sav_zdq(:,:)   
  ENDIF ! (mode == 1)

  RETURN
END SUBROUTINE add_phys_tend

SUBROUTINE diag_phys_tend (nlon, nlev, uu, vv, temp, qv, ql, qs, &
                          zdu,zdv,zdt,zdq,zdql,zdqs,paprs,text)
!======================================================================
! Ajoute les tendances des variables physiques aux variables 
! d'etat de la dynamique t_seri, q_seri ...
! On en profite pour faire des tests sur les tendances en question.
!======================================================================


!======================================================================
! Declarations
!======================================================================

USE phys_state_var_mod, ONLY : phys_tstep, ftsol
USE geometry_mod, ONLY: longitude_deg, latitude_deg
USE print_control_mod, ONLY: prt_level
USE cmp_seri_mod
USE phys_output_var_mod, ONLY : d_qw_col, d_ql_col, d_qs_col, d_qt_col, d_ek_col, d_h_dair_col &
  &           , d_h_qw_col, d_h_ql_col, d_h_qs_col, d_h_col
IMPLICIT none
  include "YOMCST.h"
  include "clesphys.h"

! Arguments :
!------------
INTEGER, INTENT(IN)                           :: nlon, nlev
REAL, DIMENSION(nlon,nlev),     INTENT(IN)    :: uu, vv
REAL, DIMENSION(nlon,nlev),     INTENT(IN)    :: temp, qv, ql, qs
REAL, DIMENSION(nlon,nlev),     INTENT(IN)    :: zdu, zdv
REAL, DIMENSION(nlon,nlev),     INTENT(IN)    :: zdt, zdq, zdql, zdqs
REAL, DIMENSION(nlon,nlev+1),   INTENT(IN)    :: paprs
CHARACTER*(*),                  INTENT(IN)    :: text

! Local :
!--------
REAL, DIMENSION(nlon,nlev)      :: uu_n, vv_n
REAL, DIMENSION(nlon,nlev)      :: temp_n, qv_n, ql_n, qs_n


!
INTEGER k, n

integer debug_level
logical, save :: first=.true.
!$OMP THREADPRIVATE(first)
!
!======================================================================
! Variables for energy conservation tests
!======================================================================
!

! zh_col-------  total enthalpy of vertical air column
! (air with watter vapour, liquid and solid) (J/m2)
! zh_dair_col--- total enthalpy of dry air (J/m2)
! zh_qw_col----  total enthalpy of watter vapour (J/m2)
! zh_ql_col----  total enthalpy of liquid watter (J/m2)
! zh_qs_col----  total enthalpy of solid watter  (J/m2)
! zqw_col------  total mass of watter vapour (kg/m2)
! zql_col------  total mass of liquid watter (kg/m2)
! zqs_col------  total mass of solid watter (kg/m2)
! zek_col------  total kinetic energy (kg/m2)
!
REAL zairm(nlon, nlev) ! layer air mass (kg/m2)
REAL zqw_col(nlon,2)
REAL zql_col(nlon,2)
REAL zqs_col(nlon,2)
REAL zek_col(nlon,2)
REAL zh_dair_col(nlon,2)
REAL zh_qw_col(nlon,2), zh_ql_col(nlon,2), zh_qs_col(nlon,2)
REAL zh_col(nlon,2)

!======================================================================
! Initialisations

     IF (prt_level >= 5) then
        write (*,*) "In diag_phys_tend, after ",text
        call flush
     end if

     debug_level=10
     if (first) then
        print *,"TestJLD rcpv, rcw, rcs",rcpv, rcw, rcs
        first=.false.
     endif
!
!  print *,'add_phys_tend: paprs ',paprs
!======================================================================
! Diagnostics for energy conservation tests
!======================================================================
  DO k = 1, nlev
    ! layer air mass
    zairm(:, k) = (paprs(:,k)-paprs(:,k+1))/rg
  END DO

  if (fl_ebil .GT. 0) then
    ! ------------------------------------------------
    ! Compute vertical sum for each atmospheric column
    ! ------------------------------------------------
    n=1   ! begining of time step

    CALL integr_v(nlon, nlev, rcpv, &
                  temp, qv, ql, qs, uu, vv, zairm, &
                    zqw_col(:,n), zql_col(:,n), zqs_col(:,n), zek_col(:,n), zh_dair_col(:,n), &
                    zh_qw_col(:,n), zh_ql_col(:,n), zh_qs_col(:,n), zh_col(:,n))

  end if ! end if (fl_ebil .GT. 0)

!======================================================================
! Ajout des tendances sur le vent, la temperature et les diverses phases de l'eau
!======================================================================

     uu_n(:,:)=uu(:,:)+zdu(:,:)
     vv_n(:,:)=vv(:,:)+zdv(:,:)
     qv_n(:,:)=qv(:,:)+zdq(:,:)
     ql_n(:,:)=ql(:,:)+zdql(:,:)
     qs_n(:,:)=qs(:,:)+zdqs(:,:)
     temp_n(:,:)=temp(:,:)+zdt(:,:)



!======================================================================
! Diagnostics for energy conservation tests
!======================================================================

  if (fl_ebil .GT. 0) then
  
    ! ------------------------------------------------
    ! Compute vertical sum for each atmospheric column
    ! ------------------------------------------------
    n=2   ! end of time step

    CALL integr_v(nlon, nlev, rcpv, &
                  temp_n, qv_n, ql_n, qs_n, uu_n, vv_n, zairm, &
                    zqw_col(:,n), zql_col(:,n), zqs_col(:,n), zek_col(:,n), zh_dair_col(:,n), &
                    zh_qw_col(:,n), zh_ql_col(:,n), zh_qs_col(:,n), zh_col(:,n))

    ! ------------------------------------------------
    ! Compute the changes by unit of time
    ! ------------------------------------------------

    d_qw_col(:) = (zqw_col(:,2)-zqw_col(:,1))/phys_tstep
    d_ql_col(:) = (zql_col(:,2)-zql_col(:,1))/phys_tstep
    d_qs_col(:) = (zqs_col(:,2)-zqs_col(:,1))/phys_tstep
    d_qt_col(:) = d_qw_col(:) + d_ql_col(:) + d_qs_col(:)

    d_ek_col(:) = (zek_col(:,2)-zek_col(:,1))/phys_tstep

   print *,'zdu ', zdu
   print *,'zdv ', zdv
   print *,'d_ek_col, zek_col(2), zek_col(1) ',d_ek_col(1), zek_col(1,2), zek_col(1,1)

    d_h_dair_col(:) = (zh_dair_col(:,2)-zh_dair_col(:,1))/phys_tstep
    d_h_qw_col(:) = (zh_qw_col(:,2)-zh_qw_col(:,1))/phys_tstep
    d_h_ql_col(:) = (zh_ql_col(:,2)-zh_ql_col(:,1))/phys_tstep
    d_h_qs_col(:) = (zh_qs_col(:,2)-zh_qs_col(:,1))/phys_tstep

    d_h_col = (zh_col(:,2)-zh_col(:,1))/phys_tstep 

  end if ! end if (fl_ebil .GT. 0)
!

  RETURN
END SUBROUTINE diag_phys_tend

SUBROUTINE integr_v(nlon, nlev, zcpvap, &
                    temp, qv, ql, qs, uu, vv, zairm,  &
                    zqw_col, zql_col, zqs_col, zek_col, zh_dair_col, &
                    zh_qw_col, zh_ql_col, zh_qs_col, zh_col)

IMPLICIT none
  include "YOMCST.h"

INTEGER,                    INTENT(IN)    :: nlon,nlev
REAL,                       INTENT(IN)    :: zcpvap
REAL, DIMENSION(nlon,nlev), INTENT(IN)    :: temp, qv, ql, qs, uu, vv
REAL, DIMENSION(nlon,nlev), INTENT(IN)    :: zairm
REAL, DIMENSION(nlon),      INTENT(OUT)   :: zqw_col
REAL, DIMENSION(nlon),      INTENT(OUT)   :: zql_col
REAL, DIMENSION(nlon),      INTENT(OUT)   :: zqs_col
REAL, DIMENSION(nlon),      INTENT(OUT)   :: zek_col
REAL, DIMENSION(nlon),      INTENT(OUT)   :: zh_dair_col
REAL, DIMENSION(nlon),      INTENT(OUT)   :: zh_qw_col
REAL, DIMENSION(nlon),      INTENT(OUT)   :: zh_ql_col
REAL, DIMENSION(nlon),      INTENT(OUT)   :: zh_qs_col
REAL, DIMENSION(nlon),      INTENT(OUT)   :: zh_col

INTEGER          :: i, k


  ! Reset variables
    zqw_col(:) = 0.
    zql_col(:) = 0.
    zqs_col(:) = 0.
    zek_col(:) = 0.
    zh_dair_col(:) = 0.
    zh_qw_col(:) = 0.
    zh_ql_col(:) = 0.
    zh_qs_col(:) = 0.

!JLD    write (*,*) "rcpd, zcpvap, zcwat, zcice ",rcpd, zcpvap, zcwat, zcice
  
    ! ------------------------------------------------
    ! Compute vertical sum for each atmospheric column
    ! ------------------------------------------------
    DO k = 1, nlev
      DO i = 1, nlon
        ! Watter mass
        zqw_col(i) = zqw_col(i) + qv(i, k)*zairm(i, k)
        zql_col(i) = zql_col(i) + ql(i, k)*zairm(i, k)
        zqs_col(i) = zqs_col(i) + qs(i, k)*zairm(i, k)
        ! Kinetic Energy
        zek_col(i) = zek_col(i) + 0.5*(uu(i,k)**2+vv(i,k)**2)*zairm(i, k)
        ! Air enthalpy : dry air, water vapour, liquid, solid
        zh_dair_col(i) = zh_dair_col(i) + rcpd*(1.-qv(i,k)-ql(i,k)-qs(i,k))* &
                                               zairm(i, k)*temp(i, k)
        zh_qw_col(i) = zh_qw_col(i) +  zcpvap*temp(i, k)         *qv(i, k)*zairm(i, k)    !jyg
        zh_ql_col(i) = zh_ql_col(i) + (zcpvap*temp(i, k) - rlvtt)*ql(i, k)*zairm(i, k)   !jyg
        zh_qs_col(i) = zh_qs_col(i) + (zcpvap*temp(i, k) - rlstt)*qs(i, k)*zairm(i, k)   !jyg
      END DO
    END DO
    ! compute total air enthalpy
    zh_col(:) = zh_dair_col(:) + zh_qw_col(:) + zh_ql_col(:) + zh_qs_col(:)

END SUBROUTINE integr_v

SUBROUTINE prt_enerbil (text, itap)
!======================================================================
! Print enenrgy budget diagnotics for the 1D case
!======================================================================

!======================================================================
! Declarations
!======================================================================

USE dimphy, ONLY: klon, klev
USE phys_state_var_mod, ONLY : phys_tstep
USE phys_state_var_mod, ONLY : topsw, toplw, solsw, sollw, rain_con, snow_con
USE geometry_mod, ONLY: longitude_deg, latitude_deg
USE print_control_mod, ONLY: prt_level
USE cmp_seri_mod
USE phys_output_var_mod, ONLY : d_qw_col, d_ql_col, d_qs_col, d_qt_col, d_ek_col, d_h_dair_col &
  &           , d_h_qw_col, d_h_ql_col, d_h_qs_col, d_h_col
USE phys_local_var_mod, ONLY: evap, sens
USE phys_local_var_mod, ONLY: u_seri, v_seri, ql_seri, qs_seri, q_seri, t_seri &
   &  , rain_lsc, snow_lsc
USE climb_hq_mod, ONLY : d_h_col_vdf, f_h_bnd
IMPLICIT none
include "YOMCST.h"

! Arguments :
!------------
CHARACTER*(*) text ! text specifing the involved parametrization 
integer itap        ! time step number
! local variables
! ---------------
real bilq_seuil,  bilh_seuil ! thresold on error in Q and H budget
real bilq_error,  bilh_error ! erros in Q and H budget
real bilq_bnd,  bilh_bnd     ! Q and H budget due to exchange with boundaries
integer bilq_ok,  bilh_ok
CHARACTER*(12) status

bilq_seuil = 1.E-10
bilh_seuil = 1.E-1
bilq_ok=0
bilh_ok=0

!!print *,'prt_level:',prt_level,' fl_ebil:',fl_ebil,' fl_cor_ebil:',fl_cor_ebil
if ( (fl_ebil .GT. 0) .and. (klon .EQ. 1)) then

  bilq_bnd = 0.
  bilh_bnd = 0.

  param: SELECT CASE (text)
  CASE("vdf") param
      bilq_bnd = evap(1)
      bilh_bnd = sens(1)+(rcpv-rcpd)*evap(1)*t_seri(1,1)
  CASE("lsc") param
      bilq_bnd = - rain_lsc(1) - snow_lsc(1)
      bilh_bnd = (-(rcw-rcpd)*t_seri(1,1) + rlvtt) * rain_lsc(1) &
    &         + (-(rcs-rcpd)*t_seri(1,1) + rlstt) * snow_lsc(1)
  CASE("convection") param
      bilq_bnd = - rain_con(1) - snow_con(1)
      bilh_bnd = (-(rcw-rcpd)*t_seri(1,1) + rlvtt) * rain_con(1) &
    &         + (-(rcs-rcpd)*t_seri(1,1) + rlstt) * snow_con(1)
  CASE("SW") param
      bilh_bnd = topsw(1) - solsw(1)
  CASE("LW") param
      bilh_bnd = -(toplw(1) + sollw(1))
  CASE DEFAULT param
      bilq_bnd = 0.
      bilh_bnd = 0.
  END SELECT param

  bilq_error = d_qt_col(1) - bilq_bnd
  bilh_error = d_h_col(1) - bilh_bnd
! are the errors too large?
  if ( abs(bilq_error) .gt. bilq_seuil) bilq_ok=1
  if ( abs(bilh_error) .gt. bilh_seuil) bilh_ok=1
!
! Print diagnostics
! =================
  if ( (bilq_ok .eq. 0).and.(bilh_ok .eq. 0) ) then
    status="enerbil-OK"
  else
    status="enerbil-PB"
  end if

  if ( prt_level .GE. 3) then
    write(*,9010) text,status," itap:",itap,"enerbilERROR: Q", bilq_error,"  H", bilh_error
9010  format (1x,A8,2x,A12,A6,I4,A18,E15.6,A5,E15.6)
  end if
  if ( prt_level .GE. 3) then
    write(*,9000) text,"enerbil: Q,H,KE budget", d_qt_col(1), d_h_col(1),d_ek_col(1)
  end if
  if ( prt_level .GE. 5) then
    write(*,9000) text,"enerbil at boundaries: Q, H",bilq_bnd, bilh_bnd
    write(*,9000) text,"enerbil: water budget",d_qt_col(1),d_qw_col(1),d_ql_col(1),d_qs_col(1)
    write(*,9000) text,"enerbil: enthalpy budget",d_h_col(1),d_h_dair_col(1),d_h_qw_col(1),d_h_ql_col(1),d_h_qs_col(1)
  end if

  specific_diag: SELECT CASE (text)
  CASE("vdf") specific_diag
    if ( prt_level .GE. 5) then
      write(*,9000) text,"enerbil: d_h, bilh, sens,t_seri", d_h_col(1), bilh_bnd, sens(1), t_seri(1,1)
      write(*,9000) text,"enerbil: d_h_col_vdf, f_h, diff",d_h_col_vdf, f_h_bnd, bilh_bnd-sens(1)
    end if
  CASE("lsc") specific_diag
    if ( prt_level .GE. 5) then
      write(*,9000) text,"enerbil: rain, bil_lat, bil_sens", rain_lsc(1), rlvtt * rain_lsc(1), -(rcw-rcpd)*t_seri(1,1) * rain_lsc(1)
      write(*,9000) text,"enerbil: snow, bil_lat, bil_sens", snow_lsc(1), rlstt * snow_lsc(1), -(rcs-rcpd)*t_seri(1,1) * snow_lsc(1)
    end if
  CASE("convection") specific_diag
    if ( prt_level .GE. 5) then
      write(*,9000) text,"enerbil: rain, bil_lat, bil_sens", rain_con(1), rlvtt * rain_con(1), -(rcw-rcpd)*t_seri(1,1) * rain_con(1)
      write(*,9000) text,"enerbil: snow, bil_lat, bil_sens", snow_con(1), rlstt * snow_con(1), -(rcs-rcpd)*t_seri(1,1) * snow_con(1)
    end if
  END SELECT specific_diag

9000 format (1x,A8,2x,A35,10E15.6)

end if ! end if (fl_ebil .GT. 0)

END SUBROUTINE prt_enerbil

END MODULE add_phys_tend_mod
