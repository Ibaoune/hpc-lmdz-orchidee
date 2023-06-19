MODULE tropopause_m

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_tropopause

CONTAINS

!-------------------------------------------------------------------------------
!
FUNCTION dyn_tropopause(t, ts, paprs, pplay, rot, itrop, thet0, pvor0)
!
!-------------------------------------------------------------------------------
  USE assert_m,     ONLY: assert
  USE assert_eq_m,  ONLY: assert_eq
  USE dimphy,       ONLY: klon, klev
  USE geometry_mod, ONLY: latitude_deg, longitude_deg
  USE vertical_layers_mod, ONLY: aps, bps, preff
#ifdef REPROBUS
  USE chem_rep, ONLY: itroprep
#endif

!-------------------------------------------------------------------------------
! Arguments:
  REAL ::     dyn_tropopause(klon) !--- Pressure at tropopause
  REAL, INTENT(IN)  ::      t(:,:) !--- Cells-centers temperature
  REAL, INTENT(IN)  ::     ts(:)   !--- Surface       temperature
  REAL, INTENT(IN)  ::  paprs(:,:) !--- Cells-edges   pressure
  REAL, INTENT(IN)  ::  pplay(:,:) !--- Cells-centers pressure
  REAL, INTENT(IN)  ::    rot(:,:) !--- Cells-centers relative vorticity
  INTEGER, INTENT(OUT), OPTIONAL :: itrop(klon) !--- Last tropospheric layer idx
  REAL,    INTENT(IN),  OPTIONAL :: thet0, pvor0
!-------------------------------------------------------------------------------
! Local variables:
  include "YOMCST.h"
  REAL, PARAMETER :: DynPTrMin =8.E+3 !--- Thresholds for minimum and maximum
  REAL, PARAMETER :: DynPTrMax =4.E+4 !    dynamical tropopause pressure (Pa).
  CHARACTER(LEN=80)  :: sub
  INTEGER :: i, k, kb, kt, kp, ib, ie, nw
  REAL    :: al, th0, pv0
  REAL,    DIMENSION(klon,klev) :: tpot_cen, tpot_edg, pvor_cen
  REAL,    PARAMETER :: sg0=0.75  !--- Start level for PV=cte search loop
  INTEGER, PARAMETER :: nadj=3    !--- Adjacent levs nb for thresholds detection
  REAL,    PARAMETER :: w(5)=[0.1,0.25,0.3,0.25,0.1] !--- Vertical smoothing
  INTEGER, SAVE :: k0
  LOGICAL, SAVE :: first=.TRUE.
!$OMP THREADPRIVATE(k0,first)
!-------------------------------------------------------------------------------
  sub='dyn_tropopause'
  CALL assert(SIZE(t ,1)==klon, TRIM(sub)//" t klon")
  CALL assert(SIZE(t ,2)==klev, TRIM(sub)//" t klev")
  CALL assert(SIZE(ts,1)==klon, TRIM(sub)//" ts klon")
  CALL assert(SHAPE(paprs)==[klon,klev+1],TRIM(sub)//" paprs shape")
  CALL assert(SHAPE(pplay)==[klon,klev  ],TRIM(sub)//" pplay shape")
  CALL assert(SHAPE(rot)  ==[klon,klev  ],TRIM(sub)//" rot shape")

  !--- DEFAULT THRESHOLDS
  th0=380.; IF(PRESENT(thet0)) th0=thet0   !--- In kelvins
  pv0=  2.; IF(PRESENT(pvor0)) pv0=pvor0   !--- In PVU
  IF(first) THEN
    DO k0=1,klev; IF(aps(k0)/preff+bps(k0)<sg0) EXIT; END DO; first=.FALSE.
  END IF

  !--- POTENTIAL TEMPERATURE AT CELLS CENTERS AND INTERFACES
  DO i = 1,klon
    tpot_cen(i,1) = t(i,1)*(preff/pplay(i,1))**RKAPPA
    tpot_edg(i,1) = ts(i) *(preff/paprs(i,1))**RKAPPA
    DO k=2,klev
      al = LOG(pplay(i,k-1)/paprs(i,k))/LOG(pplay(i,k-1)/pplay(i,k))
      tpot_cen(i,k) =  t(i,k)                        *(preff/pplay(i,k))**RKAPPA
      tpot_edg(i,k) = (t(i,k-1)+al*(t(i,k)-t(i,k-1)))*(preff/paprs(i,k))**RKAPPA
      !--- FORCE QUANTITIES TO BE GROWING
      IF(tpot_edg(i,k)<tpot_edg(i,k-1)) tpot_edg(i,k)=tpot_edg(i,k-1)+1.E-5
      IF(tpot_cen(i,k)<tpot_cen(i,k-1)) tpot_cen(i,k)=tpot_cen(i,k-1)+1.E-5
    END DO
    !--- VERTICAL SMOOTHING
    tpot_cen(i,:)=smooth(tpot_cen(i,:),w)
    tpot_edg(i,:)=smooth(tpot_edg(i,:),w)
  END DO

  !--- ERTEL POTENTIAL VORTICITY AT CELLS CENTERS (except in top layer)
  DO i = 1, klon
    DO k= 1, klev-1
      pvor_cen(i,k)=-1.E6*RG*(rot(i,k)+2.*ROMEGA*SIN(latitude_deg(i)*RPI/180.))&
                   * (tpot_edg(i,k+1)-tpot_edg(i,k)) / (paprs(i,k+1)-paprs(i,k))
    END DO
    !--- VERTICAL SMOOTHING
    pvor_cen(i,1:klev-1)=smooth(pvor_cen(i,1:klev-1),w)
  END DO

  !--- LOCATE TROPOPAUSE: LOWEST POINT BETWEEN THETA=380K AND PV=2PVU SURFACES.
  DO i = 1, klon
    !--- UPPER TROPOPAUSE: |PV|=2PVU POINT STARTING FROM TOP
    DO kt=klev-1,1,-1; IF(ALL(ABS(pvor_cen(i,kt-nadj:kt))<=pv0)) EXIT; END DO
    !--- LOWER TROPOPAUSE: |PV|=2PVU POINT STARTING FROM BOTTOM
    DO kb=k0,klev-1;   IF(ALL(ABS(pvor_cen(i,kb:kb+nadj))> pv0)) EXIT; END DO; kb=kb-1
    !--- ISO-THETA POINT: THETA=380K       STARTING FROM TOP
    DO kp=klev-1,1,-1; IF(ALL(ABS(tpot_cen(i,kp-nadj:kp))<=th0)) EXIT; END DO
    !--- CHOOSE BETWEEN LOWER AND UPPER TROPOPAUSE
    IF(2*COUNT(ABS(pvor_cen(i,kb:kt))>pv0)>kt-kb+1) kt=kb
    !--- PV-DEFINED TROPOPAUSE
    al = (ABS(pvor_cen(i,kt+1))-pv0)/ABS(pvor_cen(i,kt+1)-pvor_cen(i,kt))
    dyn_tropopause(i) = pplay(i,kt+1)*(pplay(i,kt)/pplay(i,kt+1))**al
    !--- THETA=380K IN THE TROPICAL REGION
    al = (tpot_cen(i,kp+1)-th0)/(tpot_cen(i,kp+1)-tpot_cen(i,kp))
    dyn_tropopause(i) = MAX( pplay(i,kp+1)*(pplay(i,kp)/pplay(i,kp+1))**al,    &
                            dyn_tropopause(i) )
    !--- UNREALISTIC VALUES DETECTION
    IF(dyn_tropopause(i)<DynPTrMin.OR.dyn_tropopause(i)>DynPTrMax) THEN
      dyn_tropopause(i)=MIN(MAX(dyn_tropopause(i),DynPTrMax),DynPTrMin)
      DO kt=1,klev-1; IF(pplay(i,kt+1)>dyn_tropopause(i)) EXIT; END DO; kp=kt
    END IF
#ifdef REPROBUS
    itroprep(i)=MAX(kt,kp)
#endif
    !--- LAST TROPOSPHERIC LAYER INDEX NEEDED
    IF(PRESENT(itrop)) itrop(i)=MAX(kt,kp)
  END DO

END FUNCTION dyn_tropopause


!-------------------------------------------------------------------------------
!
FUNCTION smooth(v,w)
!
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN)         :: v(:), w(:)
  REAL, DIMENSION(SIZE(v)) :: smooth
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: nv, nw, k, kb, ke, lb, le
!-------------------------------------------------------------------------------
  nv=SIZE(v); nw=(SIZE(w)-1)/2
  DO k=1,nv
    kb=MAX(k-nw,1 ); lb=MAX(2+nw   -k,1)
    ke=MIN(k+nw,nv); le=MIN(1+nw+nv-k,1+2*nw)
    smooth(k)=SUM(v(kb:ke)*w(lb:le))/SUM(w(lb:le))
  END DO

END FUNCTION smooth
!
!-------------------------------------------------------------------------------

END MODULE tropopause_m
