!
! $Id: inidissip.F90 2603 2016-07-25 09:31:56Z emillour $
!
SUBROUTINE inidissip ( lstardis,nitergdiv,nitergrot,niterh  , &
     tetagdiv,tetagrot,tetatemp, vert_prof_dissip)
  !=======================================================================
  !   initialisation de la dissipation horizontale
  !=======================================================================
  !-----------------------------------------------------------------------
  !   declarations:
  !   -------------

  USE control_mod, only : dissip_period,iperiod
  USE comconst_mod, ONLY: dissip_deltaz, dissip_factz, dissip_zref, &
                          dtdiss, dtvr, rad
  USE comvert_mod, ONLY: preff, presnivs

  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comdissipn.h"
  include "iniprint.h"

  LOGICAL,INTENT(in) :: lstardis
  INTEGER,INTENT(in) :: nitergdiv,nitergrot,niterh
  REAL,INTENT(in) :: tetagdiv,tetagrot,tetatemp

  integer, INTENT(in):: vert_prof_dissip
  ! Vertical profile of horizontal dissipation
  ! Allowed values:
  ! 0: rational fraction, function of pressure
  ! 1: tanh of altitude

! Local variables:
  REAL fact,zvert(llm),zz
  REAL zh(ip1jmp1),zu(ip1jmp1), gx(ip1jmp1), divgra(ip1jmp1)
  real zv(ip1jm), gy(ip1jm), deltap(ip1jmp1,llm)
  REAL ullm,vllm,umin,vmin,zhmin,zhmax
  REAL zllm

  INTEGER l,ij,idum,ii
  REAL tetamin
  REAL pseudoz
  character (len=80) :: abort_message

  REAL ran1


  !-----------------------------------------------------------------------
  !
  !   calcul des valeurs propres des operateurs par methode iterrative:
  !   -----------------------------------------------------------------

  crot     = -1.
  cdivu    = -1.
  cdivh    = -1.

  !   calcul de la valeur propre de divgrad:
  !   --------------------------------------
  idum = 0
  DO l = 1, llm
     DO ij = 1, ip1jmp1
        deltap(ij,l) = 1.
     ENDDO
  ENDDO

  idum  = -1
  zh(1) = RAN1(idum)-.5
  idum  = 0
  DO ij = 2, ip1jmp1
     zh(ij) = RAN1(idum) -.5
  ENDDO

  CALL filtreg (zh,jjp1,1,2,1,.TRUE.,1)

  CALL minmax(iip1*jjp1,zh,zhmin,zhmax )

  IF ( zhmin .GE. zhmax  )     THEN
     write(lunout,*)'  Inidissip  zh min max  ',zhmin,zhmax
     abort_message='probleme generateur alleatoire dans inidissip'
     call abort_gcm('inidissip',abort_message,1)
  ENDIF

  zllm = ABS( zhmax )
  DO l = 1,50
     IF(lstardis) THEN
        CALL divgrad2(1,zh,deltap,niterh,divgra)
     ELSE
        CALL divgrad (1,zh,niterh,divgra)
     ENDIF

     zllm  = ABS(maxval(divgra))
     zh = divgra / zllm
  ENDDO

  IF(lstardis) THEN
     cdivh = 1./ zllm
  ELSE
     cdivh = zllm ** ( -1./niterh )
  ENDIF

  !   calcul des valeurs propres de gradiv (ii =1) et  nxgrarot(ii=2)
  !   -----------------------------------------------------------------
  write(lunout,*)'inidissip: calcul des valeurs propres'

  DO    ii = 1, 2
     !
     DO ij = 1, ip1jmp1
        zu(ij)  = RAN1(idum) -.5
     ENDDO
     CALL filtreg (zu,jjp1,1,2,1,.TRUE.,1)
     DO ij = 1, ip1jm
        zv(ij) = RAN1(idum) -.5
     ENDDO
     CALL filtreg (zv,jjm,1,2,1,.FALSE.,1)

     CALL minmax(iip1*jjp1,zu,umin,ullm )
     CALL minmax(iip1*jjm, zv,vmin,vllm )

     ullm = ABS ( ullm )
     vllm = ABS ( vllm )

     DO    l = 1, 50
        IF(ii.EQ.1) THEN
           !cccc             CALL covcont( 1,zu,zv,zu,zv )
           IF(lstardis) THEN
              CALL gradiv2( 1,zu,zv,nitergdiv,gx,gy )
           ELSE
              CALL gradiv ( 1,zu,zv,nitergdiv,gx,gy )
           ENDIF
        ELSE
           IF(lstardis) THEN
              CALL nxgraro2( 1,zu,zv,nitergrot,gx,gy )
           ELSE
              CALL nxgrarot( 1,zu,zv,nitergrot,gx,gy )
           ENDIF
        ENDIF

        zllm = max(abs(maxval(gx)), abs(maxval(gy)))
        zu = gx / zllm
        zv = gy / zllm
     end DO

     IF ( ii.EQ.1 ) THEN
        IF(lstardis) THEN
           cdivu  = 1./zllm
        ELSE
           cdivu  = zllm **( -1./nitergdiv )
        ENDIF
     ELSE
        IF(lstardis) THEN
           crot   = 1./ zllm
        ELSE
           crot   = zllm **( -1./nitergrot )
        ENDIF
     ENDIF

  end DO

  !   petit test pour les operateurs non star:
  !   ----------------------------------------

  !     IF(.NOT.lstardis) THEN
  fact    = rad*24./REAL(jjm)
  fact    = fact*fact
  write(lunout,*)'inidissip: coef u ', fact/cdivu, 1./cdivu
  write(lunout,*)'inidissip: coef r ', fact/crot , 1./crot
  write(lunout,*)'inidissip: coef h ', fact/cdivh, 1./cdivh
  !     ENDIF

  !-----------------------------------------------------------------------
  !   variation verticale du coefficient de dissipation:
  !   --------------------------------------------------

  if (vert_prof_dissip == 1) then
     do l=1,llm
        pseudoz=8.*log(preff/presnivs(l))
        zvert(l)=1+ &
             (tanh((pseudoz-dissip_zref)/dissip_deltaz)+1.)/2. &
             *(dissip_factz-1.)
     enddo
  else
     DO l=1,llm
        zvert(l)=1.
     ENDDO
     fact=2.
     DO l = 1, llm
        zz      = 1. - preff/presnivs(l)
        zvert(l)= fact -( fact-1.)/( 1.+zz*zz )
     ENDDO
  endif


  write(lunout,*)'inidissip: Constantes de temps de la diffusion horizontale'

  tetamin =  1.e+6

  DO l=1,llm
     tetaudiv(l)   = zvert(l)/tetagdiv
     tetaurot(l)   = zvert(l)/tetagrot
     tetah(l)      = zvert(l)/tetatemp

     IF( tetamin.GT. (1./tetaudiv(l)) ) tetamin = 1./ tetaudiv(l)
     IF( tetamin.GT. (1./tetaurot(l)) ) tetamin = 1./ tetaurot(l)
     IF( tetamin.GT. (1./   tetah(l)) ) tetamin = 1./    tetah(l)
  ENDDO

  ! If dissip_period=0 calculate value for dissipation period, else keep value read from gcm.def
  IF (dissip_period == 0) THEN
     dissip_period = INT( tetamin/( 2.*dtvr*iperiod) ) * iperiod
     write(lunout,*)'inidissip: tetamin dtvr iperiod dissip_period(intermed) ',tetamin,dtvr,iperiod,dissip_period
     dissip_period = MAX(iperiod,dissip_period)
  END IF

  dtdiss  = dissip_period * dtvr
  write(lunout,*)'inidissip: dissip_period=',dissip_period,' dtdiss=',dtdiss,' dtvr=',dtvr

  DO l = 1,llm
     write(lunout,*)zvert(l),dtdiss*tetaudiv(l),dtdiss*tetaurot(l), &
          dtdiss*tetah(l)
  ENDDO

END SUBROUTINE inidissip
