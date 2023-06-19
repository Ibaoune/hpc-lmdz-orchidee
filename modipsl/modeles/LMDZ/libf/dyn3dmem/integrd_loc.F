!
! $Id: integrd_p.F 1299 2010-01-20 14:27:21Z fairhead $
!
      SUBROUTINE integrd_loc
     $  (  nq,vcovm1,ucovm1,tetam1,psm1,massem1,
     $     dv,du,dteta,dq,dp,vcov,ucov,teta,q,ps0,masse,phis) !,finvmaold)
      USE parallel_lmdz
      USE control_mod
      USE mod_filtreg_p
      USE write_field_loc
      USE write_field
      USE integrd_mod
      USE infotrac, ONLY: ok_iso_verif ! ajout CRisi
      USE comconst_mod, ONLY: pi
      USE logic_mod, ONLY: leapf
      USE comvert_mod, ONLY: ap, bp
      USE temps_mod, ONLY: dt
      
      IMPLICIT NONE


c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   objet:
c   ------
c
c   Incrementation des tendances dynamiques
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

      include "dimensions.h"
      include "paramet.h"
      include "comgeom.h"
      include "iniprint.h"

c   Arguments:
c   ----------

      INTEGER,intent(in) :: nq ! number of tracers to handle in this routine

      REAL,INTENT(INOUT) :: vcov(ijb_v:ije_v,llm) ! covariant meridional wind
      REAL,INTENT(INOUT) :: ucov(ijb_u:ije_u,llm) ! covariant zonal wind
      REAL,INTENT(INOUT) :: teta(ijb_u:ije_u,llm) ! potential temperature
      REAL,INTENT(INOUT) :: q(ijb_u:ije_u,llm,nq) ! advected tracers
      REAL,INTENT(INOUT) :: ps0(ijb_u:ije_u) ! surface pressure
      REAL,INTENT(INOUT) :: masse(ijb_u:ije_u,llm) ! atmospheric mass
      REAL,INTENT(INOUT) :: phis(ijb_u:ije_u) ! ground geopotential !!! unused
      ! values at previous time step
      REAL,INTENT(INOUT) :: vcovm1(ijb_v:ije_v,llm)
      REAL,INTENT(INOUT) :: ucovm1(ijb_u:ije_u,llm)
      REAL,INTENT(INOUT) :: tetam1(ijb_u:ije_u,llm)
      REAL,INTENT(INOUT) :: psm1(ijb_u:ije_u)
      REAL,INTENT(INOUT) :: massem1(ijb_u:ije_u,llm)
      ! the tendencies to add
      REAL,INTENT(INOUT) :: dv(ijb_v:ije_v,llm)
      REAL,INTENT(INOUT) :: du(ijb_u:ije_u,llm)
      REAL,INTENT(INOUT) :: dteta(ijb_u:ije_u,llm)
      REAL,INTENT(INOUT) :: dp(ijb_u:ije_u)
      REAL,INTENT(INOUT) :: dq(ijb_u:ije_u,llm,nq) !!! unused
!      REAL,INTENT(INOUT) ::finvmaold(ijb_u:ije_u,llm) !!! unused

c   Local:
c   ------

      REAL vscr( ijb_v:ije_v ),uscr( ijb_u:ije_u )
      REAL hscr( ijb_u:ije_u ),pscr(ijb_u:ije_u)
      REAL massescr( ijb_u:ije_u,llm )
!      REAL finvmasse(ijb_u:ije_u,llm)
      REAL tpn,tps,tppn(iim),tpps(iim)
      REAL qpn,qps,qppn(iim),qpps(iim)

      INTEGER  l,ij,iq,i,j

      REAL SSUM
      EXTERNAL SSUM
      INTEGER ijb,ije,jjb,jje
      LOGICAL :: checksum
      LOGICAL,SAVE :: checksum_all=.TRUE.
      INTEGER :: stop_it
      INTEGER :: ierr

      !write(*,*) 'integrd 88: entree, nq=',nq
c-----------------------------------------------------------------------

c$OMP BARRIER     
      if (pole_nord) THEN
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO  l = 1,llm
          DO  ij = 1,iip1
           ucov(    ij    , l) = 0.
           uscr(     ij      ) = 0.
           ENDDO
        ENDDO
c$OMP END DO NOWAIT        
      ENDIF

      if (pole_sud) THEN
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
        DO  l = 1,llm
          DO  ij = 1,iip1
           ucov( ij +ip1jm, l) = 0.
           uscr( ij +ip1jm   ) = 0.
          ENDDO
        ENDDO
c$OMP END DO NOWAIT      
      ENDIF

c    ............    integration  de       ps         ..............

c      CALL SCOPY(ip1jmp1*llm, masse, 1, massescr, 1)

      ijb=ij_begin
      ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO  l = 1,llm
        massescr(ijb:ije,l)=masse(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT 

c$OMP DO SCHEDULE(STATIC)
      DO 2 ij = ijb,ije
       pscr (ij)    = ps0(ij)
       ps (ij)      = psm1(ij) + dt * dp(ij)     

   2  CONTINUE

c$OMP END DO  
c$OMP BARRIER
c --> ici synchro OPENMP pour ps
       
      checksum=.TRUE.
      stop_it=0

c$OMP MASTER
!c$OMP DO SCHEDULE(STATIC)
      DO ij = ijb,ije
         IF( ps(ij).LT.0. ) THEN
           IF (checksum) stop_it=ij
           checksum=.FALSE.
         ENDIF
       ENDDO
!c$OMP END DO NOWAIT 
       
!      CALL MPI_ALLREDUCE(checksum,checksum_all,1,
!     &                   MPI_LOGICAL,MPI_LOR,COMM_LMDZ,ierr)
      IF( .NOT. checksum ) THEN
         write(lunout,*) "integrd: ps = ", ps(stop_it)
         write(lunout,*) " at node ij =", stop_it
         ! since ij=j+(i-1)*jjp1 , we have
          j=modulo(stop_it,jjp1)
          i=1+(stop_it-j)/jjp1
          write(lunout,*) " lon = ",rlonv(i)*180./pi, " deg",
     &                    " lat = ",rlatu(j)*180./pi, " deg"
         call abort_gcm("integrd_loc", "negative surface pressure", 1)
      ENDIF

c$OMP END MASTER
c$OMP BARRIER
        !write(*,*) 'integrd 170'
      IF (.NOT. Checksum_all) THEN
        call WriteField_v('int_vcov',vcov)
        call WriteField_u('int_ucov',ucov)
        call WriteField_u('int_teta',teta)
        call WriteField_u('int_ps0',ps0)
        call WriteField_u('int_masse',masse)
        call WriteField_u('int_phis',phis)
        call WriteField_v('int_vcovm1',vcovm1)
        call WriteField_u('int_ucovm1',ucovm1)
        call WriteField_u('int_tetam1',tetam1)
        call WriteField_u('int_psm1',psm1)
        call WriteField_u('int_massem1',massem1)

        call WriteField_v('int_dv',dv)
        call WriteField_u('int_du',du)
        call WriteField_u('int_dteta',dteta)
        call WriteField_u('int_dp',dp)
!        call WriteField_u('int_finvmaold',finvmaold)
        do j=1,nq
          call WriteField_u('int_q'//trim(int2str(j)),
     .                q(:,:,j))
          call WriteField_u('int_dq'//trim(int2str(j)),
     .                dq(:,:,j))
        enddo
        call abort_gcm("integrd_loc", "", 1)
      ENDIF
    
       
c
        !write(*,*) 'integrd 200'
C$OMP MASTER
      if (pole_nord) THEN
      
        DO  ij    = 1, iim
         tppn(ij) = aire(   ij   ) * ps(  ij    )
        ENDDO
         tpn      = SSUM(iim,tppn,1)/apoln
        DO ij   = 1, iip1
         ps(   ij   )  = tpn
        ENDDO
      
      ENDIF
      
      if (pole_sud) THEN
      
        DO  ij    = 1, iim
         tpps(ij) = aire(ij+ip1jm) * ps(ij+ip1jm)
        ENDDO
         tps      = SSUM(iim,tpps,1)/apols
        DO ij   = 1, iip1
         ps(ij+ip1jm)  = tps
        ENDDO
      
      ENDIF
c$OMP END MASTER
c$OMP BARRIER
      !write(*,*) 'integrd 217'  
c
c  ... Calcul  de la nouvelle masse d'air au dernier temps integre t+1 ...
c

      CALL pression_loc ( ip1jmp1, ap, bp, ps, p )

c$OMP BARRIER
      CALL massdair_loc (     p  , masse         )

! Ehouarn : we don't use/need finvmaold and finvmasse,
!           so might as well not compute them
!c      CALL   SCOPY( ijp1llm  , masse, 1, finvmasse,  1      )
!      ijb=ij_begin
!      ije=ij_end
!      
!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
!      DO  l = 1,llm
!        finvmasse(ijb:ije,l)=masse(ijb:ije,l)
!      ENDDO
!c$OMP END DO NOWAIT

!      jjb=jj_begin
!      jje=jj_end
!      CALL filtreg_p( finvmasse,jjb_u,jje_u,jjb,jje, jjp1, llm,
!     &                -2, 2, .TRUE., 1  )
c

c    ............   integration  de  ucov, vcov,  h     ..............

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO 10 l = 1,llm
      
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
      DO 4 ij = ijb,ije
      uscr( ij )   =  ucov( ij,l )
      ucov( ij,l ) = ucovm1( ij,l ) + dt * du( ij,l )
   4  CONTINUE

      ijb=ij_begin
      ije=ij_end
      if (pole_sud)  ije=ij_end-iip1
      
      DO 5 ij = ijb,ije
      vscr( ij )   =  vcov( ij,l )
      vcov( ij,l ) = vcovm1( ij,l ) + dt * dv( ij,l )
   5  CONTINUE
      
      ijb=ij_begin
      ije=ij_end
      
      DO 6 ij = ijb,ije
      hscr( ij )    =  teta(ij,l)
      teta ( ij,l ) = tetam1(ij,l) *  massem1(ij,l) / masse(ij,l) 
     $                + dt * dteta(ij,l) / masse(ij,l)
   6  CONTINUE

c   ....  Calcul de la valeur moyenne, unique  aux poles pour  teta    ......
c
c
        !write(*,*) 'integrd 291'
      IF (pole_nord) THEN
       
        DO  ij   = 1, iim
          tppn(ij) = aire(   ij   ) * teta(  ij    ,l)
        ENDDO
          tpn      = SSUM(iim,tppn,1)/apoln

        DO ij   = 1, iip1
          teta(   ij   ,l)  = tpn
        ENDDO
      
      ENDIF
      
      IF (pole_sud) THEN
       
        DO  ij   = 1, iim
          tpps(ij) = aire(ij+ip1jm) * teta(ij+ip1jm,l)
        ENDDO
          tps      = SSUM(iim,tpps,1)/apols

        DO ij   = 1, iip1
          teta(ij+ip1jm,l)  = tps
        ENDDO
      
      ENDIF
c

      IF(leapf)  THEN
c         CALL SCOPY ( ip1jmp1, uscr(1), 1, ucovm1(1, l), 1 )
c         CALL SCOPY (   ip1jm, vscr(1), 1, vcovm1(1, l), 1 )
c         CALL SCOPY ( ip1jmp1, hscr(1), 1, tetam1(1, l), 1 )
        ijb=ij_begin
        ije=ij_end
        ucovm1(ijb:ije,l)=uscr(ijb:ije)
        tetam1(ijb:ije,l)=hscr(ijb:ije)
        if (pole_sud) ije=ij_end-iip1
        vcovm1(ijb:ije,l)=vscr(ijb:ije)
      
      END IF

  10  CONTINUE
c$OMP END DO NOWAIT

c
c   .......  integration de   q   ......
c
      ijb=ij_begin
      ije=ij_end

         if (planet_type.eq."earth") then
! Earth-specific treatment of first 2 tracers (water)
c$OMP BARRIER
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
          DO l = 1, llm
           DO ij = ijb, ije
            deltap(ij,l) =  p(ij,l) - p(ij,l+1) 
           ENDDO
          ENDDO
          
c$OMP END DO NOWAIT
c$OMP BARRIER

        if (ok_iso_verif) then
           call check_isotopes(q,ijb,ije,'integrd 342')
        endif !if (ok_iso_verif) then

        !write(*,*) 'integrd 341'
        CALL qminimum_loc( q, nq, deltap )
        !write(*,*) 'integrd 343'

        if (ok_iso_verif) then
           call check_isotopes(q,ijb,ije,'integrd 346')
        endif !if (ok_iso_verif) then
c
c    .....  Calcul de la valeur moyenne, unique  aux poles pour  q .....
c
c$OMP BARRIER
      IF (pole_nord) THEN 
      
        DO iq = 1, nq
        
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l = 1, llm
  
             DO ij = 1, iim
               qppn(ij) = aire(   ij   ) * q(   ij   ,l,iq)
             ENDDO
               qpn  =  SSUM(iim,qppn,1)/apoln
      
             DO ij = 1, iip1
               q(   ij   ,l,iq)  = qpn
             ENDDO    
  
          ENDDO
c$OMP END DO NOWAIT

        ENDDO
      
      ENDIF

      IF (pole_sud) THEN 
      
        DO iq = 1, nq

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l = 1, llm
  
             DO ij = 1, iim
               qpps(ij) = aire(ij+ip1jm) * q(ij+ip1jm,l,iq)
             ENDDO
               qps  =  SSUM(iim,qpps,1)/apols 
  
             DO ij = 1, iip1
               q(ij+ip1jm,l,iq)  = qps
             ENDDO    
  
          ENDDO
c$OMP END DO NOWAIT

        ENDDO
      
      ENDIF

        if (ok_iso_verif) then
           call check_isotopes(q,ijb,ije,'integrd 409')
        endif !if (ok_iso_verif) then
      
! Ehouarn: forget about finvmaold
!c         CALL  SCOPY( ijp1llm , finvmasse, 1, finvmaold, 1 )

!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
!      DO l = 1, llm      
!        finvmaold(ijb:ije,l)=finvmasse(ijb:ije,l)        
!      ENDDO
!c$OMP END DO NOWAIT

      endif ! of if (planet_type.eq."earth")

c
c
c     .....   FIN  de l'integration  de   q    .......

15    continue
          !write(*,*) 'integrd 410'

c$OMP DO SCHEDULE(STATIC)
      DO ij=ijb,ije  
        ps0(ij)=ps(ij)
      ENDDO
c$OMP END DO NOWAIT 

c    .................................................................


      IF( leapf )  THEN
c       CALL SCOPY (    ip1jmp1 ,  pscr   , 1,   psm1  , 1 )
c       CALL SCOPY ( ip1jmp1*llm, massescr, 1,  massem1, 1 )
c$OMP DO SCHEDULE(STATIC)
      DO ij=ijb,ije  
        psm1(ij)=pscr(ij)
      ENDDO
c$OMP END DO NOWAIT 

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l = 1, llm
            massem1(ijb:ije,l)=massescr(ijb:ije,l)
          ENDDO
c$OMP END DO NOWAIT          
      END IF
c$OMP BARRIER
      RETURN
      END
