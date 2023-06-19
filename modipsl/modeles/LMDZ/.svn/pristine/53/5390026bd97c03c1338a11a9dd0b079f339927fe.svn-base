!
! $Id $
!
      SUBROUTINE dynetat0(fichnom,vcov,ucov,
     .                    teta,q,masse,ps,phis,time)

      USE infotrac
      use netcdf, only: nf90_get_var

      use control_mod, only : planet_type
      USE comvert_mod, ONLY: pa,preff
      USE comconst_mod, ONLY: cpp, daysec, dtvr, g, im, jm, kappa,
     &                        lllm, omeg, rad
      USE logic_mod, ONLY: fxyhypb, ysinus
      USE serre_mod, ONLY: clon,clat,grossismx,grossismy
      USE temps_mod, ONLY: annee_ref,day_ref,itau_dyn,
     &                     start_time,day_ini,hour_ini
      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0

      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van / L.Fairhead
c   -------
c
c   objet:
c   ------
c
c   Lecture de l'etat initial
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"
#include "netcdf.inc"
#include "description.h"
#include "iniprint.h"

c   Arguments:
c   ----------

      CHARACTER*(*) fichnom
      REAL vcov(iip1, jjm,llm),ucov(iip1, jjp1,llm),teta(iip1, jjp1,llm)
      REAL q(iip1,jjp1,llm,nqtot),masse(iip1, jjp1,llm)
      REAL ps(iip1, jjp1),phis(iip1, jjp1)

      REAL time

c   Variables 
c
      INTEGER length,iq
      PARAMETER (length = 100)
      REAL tab_cntrl(length) ! tableau des parametres du run
      INTEGER ierr, nid, nvarid

      INTEGER idecal

c-----------------------------------------------------------------------

c  Ouverture NetCDF du fichier etat initial

      ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
      IF (ierr.NE.NF_NOERR) THEN
        write(lunout,*)'dynetat0: Pb d''ouverture du fichier start.nc'
        write(lunout,*)' ierr = ', ierr
        CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

c
      ierr = NF_INQ_VARID (nid, "controle", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <controle> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, tab_cntrl)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echoue pour <controle>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      !!! AS: idecal is a hack to be able to read planeto starts...
      !!!     .... while keeping everything OK for LMDZ EARTH
      if (planet_type.eq."generic") then
          print*,'NOTE NOTE NOTE : Planeto-like start files'
          idecal = 4
          annee_ref  = 2000
      else
          print*,'NOTE NOTE NOTE : Earth-like start files'
          idecal = 5
          annee_ref  = tab_cntrl(5)
      endif


      im         = tab_cntrl(1)
      jm         = tab_cntrl(2)
      lllm       = tab_cntrl(3)
      day_ref    = tab_cntrl(4)
      rad        = tab_cntrl(idecal+1)
      omeg       = tab_cntrl(idecal+2)
      g          = tab_cntrl(idecal+3)
      cpp        = tab_cntrl(idecal+4)
      kappa      = tab_cntrl(idecal+5)
      daysec     = tab_cntrl(idecal+6)
      dtvr       = tab_cntrl(idecal+7)
      etot0      = tab_cntrl(idecal+8)
      ptot0      = tab_cntrl(idecal+9)
      ztot0      = tab_cntrl(idecal+10)
      stot0      = tab_cntrl(idecal+11)
      ang0       = tab_cntrl(idecal+12)
      pa         = tab_cntrl(idecal+13)
      preff      = tab_cntrl(idecal+14)
c
      clon       = tab_cntrl(idecal+15)
      clat       = tab_cntrl(idecal+16)
      grossismx  = tab_cntrl(idecal+17)
      grossismy  = tab_cntrl(idecal+18)
c
      IF ( tab_cntrl(idecal+19).EQ.1. )  THEN
        fxyhypb  = . TRUE .
c        dzoomx   = tab_cntrl(25)
c        dzoomy   = tab_cntrl(26)
c        taux     = tab_cntrl(28)
c        tauy     = tab_cntrl(29)
      ELSE
        fxyhypb = . FALSE .
        ysinus  = . FALSE .
        IF( tab_cntrl(idecal+22).EQ.1. ) ysinus = . TRUE. 
      ENDIF

      day_ini = tab_cntrl(30)
      itau_dyn = tab_cntrl(31)
      start_time = tab_cntrl(32)
c   .................................................................
c
c
      write(lunout,*)'dynetat0: rad,omeg,g,cpp,kappa',
     &               rad,omeg,g,cpp,kappa

      IF(   im.ne.iim           )  THEN
          PRINT 1,im,iim
          STOP
      ELSE  IF( jm.ne.jjm       )  THEN
          PRINT 2,jm,jjm
          STOP
      ELSE  IF( lllm.ne.llm     )  THEN
          PRINT 3,lllm,llm
          STOP
      ENDIF

      ierr = NF_INQ_VARID (nid, "rlonu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <rlonu> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, rlonu)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <rlonu>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "rlatu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <rlatu> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, rlatu)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <rlatu>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "rlonv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <rlonv> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, rlonv)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <rlonv>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "rlatv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <rlatv> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, rlatv)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour rlatv"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "cu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <cu> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, cu)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <cu>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "cv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <cv> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, cv)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <cv>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "aire", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <aire> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, aire)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <aire>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "phisinit", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <phisinit> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, phis)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <phisinit>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "temps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <temps> est absent"
         write(lunout,*)"dynetat0: J essaie <Time>"
         ierr = NF_INQ_VARID (nid, "Time", nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            write(lunout,*)"dynetat0: Le champ <Time> est absent"
            CALL ABORT_gcm("dynetat0", "", 1)
         ENDIF
      ENDIF
      ierr = nf90_get_var(nid, nvarid, time)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee <temps>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "ucov", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <ucov> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, ucov)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <ucov>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
 
      ierr = NF_INQ_VARID (nid, "vcov", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <vcov> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, vcov)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <vcov>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "teta", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <teta> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, teta)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <teta>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF


      IF(nqtot.GE.1) THEN
      DO iq=1,nqtot
        ierr =  NF_INQ_VARID (nid, tname(iq), nvarid)
        IF (ierr .NE. NF_NOERR) THEN
           write(lunout,*)"dynetat0: Le traceur <"//trim(tname(iq))//
     &                    "> est absent"
           write(lunout,*)"          Il est donc initialise a zero"
           q(:,:,:,iq)=0.
        ELSE
           ierr = NF90_GET_VAR(nid, nvarid, q(:,:,:,iq))
          IF (ierr .NE. NF_NOERR) THEN
            write(lunout,*)"dynetat0: Lecture echouee pour "//tname(iq)
            CALL ABORT_gcm("dynetat0", "", 1)
          ENDIF
        ENDIF
      ENDDO
      ENDIF

      ierr = NF_INQ_VARID (nid, "masse", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <masse> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, masse)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <masse>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_INQ_VARID (nid, "ps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <ps> est absent"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF
      ierr = nf90_get_var(nid, nvarid, ps)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <ps>"
         CALL ABORT_gcm("dynetat0", "", 1)
      ENDIF

      ierr = NF_CLOSE(nid)

       day_ini=day_ini+INT(time)
       time=time-INT(time)

  1   FORMAT(//10x,'la valeur de im =',i4,2x,'lue sur le fichier de dem
     *arrage est differente de la valeur parametree iim =',i4//)
   2  FORMAT(//10x,'la valeur de jm =',i4,2x,'lue sur le fichier de dem
     *arrage est differente de la valeur parametree jjm =',i4//)
   3  FORMAT(//10x,'la valeur de lmax =',i4,2x,'lue sur le fichier dema
     *rrage est differente de la valeur parametree llm =',i4//)
   4  FORMAT(//10x,'la valeur de dtrv =',i4,2x,'lue sur le fichier dema
     *rrage est differente de la valeur  dtinteg =',i4//)

      RETURN
      END
