!
! $Id$
!
      SUBROUTINE addfi_loc(pdt, leapf, forward,
     S          pucov, pvcov, pteta, pq   , pps ,
     S          pdufi, pdvfi, pdhfi,pdqfi, pdpfi  )
      USE parallel_lmdz
      USE infotrac, ONLY : nqtot
      USE control_mod, ONLY : planet_type
      IMPLICIT NONE
c
c=======================================================================
c
c    Addition of the physical tendencies
c
c    Interface :
c    -----------
c
c      Input :
c      -------
c      pdt                    time step of integration
c      leapf                  logical
c      forward                logical
c      pucov(ip1jmp1,llm)     first component of the covariant velocity
c      pvcov(ip1ip1jm,llm)    second component of the covariant velocity
c      pteta(ip1jmp1,llm)     potential temperature
c      pts(ip1jmp1,llm)       surface temperature
c      pdufi(ip1jmp1,llm)     |
c      pdvfi(ip1jm,llm)       |   respective
c      pdhfi(ip1jmp1)         |      tendencies
c      pdtsfi(ip1jmp1)        |
c
c      Output :
c      --------
c      pucov
c      pvcov
c      ph
c      pts
c
c
c=======================================================================
c
c-----------------------------------------------------------------------
c
c    0.  Declarations :
c    ------------------
c
      include "dimensions.h"
      include "paramet.h"
      include "comgeom.h"
c
c    Arguments :
c    -----------
c
      REAL,INTENT(IN) :: pdt ! time step for the integration (s)
c
      REAL,INTENT(INOUT) :: pvcov(ijb_v:ije_v,llm) ! covariant meridional wind
      REAL,INTENT(INOUT) :: pucov(ijb_u:ije_u,llm) ! covariant zonal wind
      REAL,INTENT(INOUT) :: pteta(ijb_u:ije_u,llm) ! potential temperature
      REAL,INTENT(INOUT) :: pq(ijb_u:ije_u,llm,nqtot) ! tracers
      REAL,INTENT(INOUT) :: pps(ijb_u:ije_u) ! surface pressure (Pa)
c respective tendencies (.../s) to add
      REAL,INTENT(IN) :: pdvfi(ijb_v:ije_v,llm)
      REAL,INTENT(IN) :: pdufi(ijb_u:ije_u,llm)
      REAL,INTENT(IN) :: pdqfi(ijb_u:ije_u,llm,nqtot)
      REAL,INTENT(IN) :: pdhfi(ijb_u:ije_u,llm)
      REAL,INTENT(IN) :: pdpfi(ijb_u:ije_u)
c
      LOGICAL,INTENT(IN) :: leapf,forward ! not used
c
c
c    Local variables :
c    -----------------
c
      REAL xpn(iim),xps(iim),tpn,tps
      INTEGER j,k,iq,ij
      REAL,PARAMETER :: qtestw = 1.0e-15
      REAL,PARAMETER :: qtestt = 1.0e-40

      REAL SSUM
      EXTERNAL SSUM
      
      INTEGER :: ijb,ije
c
c-----------------------------------------------------------------------
      
      ijb=ij_begin
      ije=ij_end
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO k = 1,llm
         DO j = ijb,ije
            pteta(j,k)= pteta(j,k) + pdhfi(j,k) * pdt
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

      if (pole_nord) then
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO  k    = 1, llm
         DO  ij   = 1, iim
           xpn(ij) = aire(   ij   ) * pteta(  ij    ,k)
         ENDDO
         tpn      = SSUM(iim,xpn,1)/ apoln

         DO ij   = 1, iip1
           pteta(   ij   ,k)  = tpn
         ENDDO
       ENDDO
c$OMP END DO NOWAIT
      endif

      if (pole_sud) then
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO  k    = 1, llm
         DO  ij   = 1, iim
           xps(ij) = aire(ij+ip1jm) * pteta(ij+ip1jm,k)
         ENDDO
         tps      = SSUM(iim,xps,1)/ apols

         DO ij   = 1, iip1
           pteta(ij+ip1jm,k)  = tps
         ENDDO
       ENDDO
c$OMP END DO NOWAIT
      endif
c

      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO k = 1,llm
         DO j = ijb,ije
            pucov(j,k)= pucov(j,k) + pdufi(j,k) * pdt
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

      if (pole_nord) ijb=ij_begin

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO k = 1,llm
         DO j = ijb,ije
            pvcov(j,k)= pvcov(j,k) + pdvfi(j,k) * pdt
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

c
      if (pole_sud)  ije=ij_end
c$OMP MASTER
      DO j = ijb,ije
         pps(j) = pps(j) + pdpfi(j) * pdt
      ENDDO
c$OMP END MASTER
 
      if (planet_type=="earth") then
      ! earth case, special treatment for first 2 tracers (water)
      DO iq = 1, 2
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO k = 1,llm
            DO j = ijb,ije
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestw )
            ENDDO
         ENDDO
c$OMP END DO NOWAIT
      ENDDO

      DO iq = 3, nqtot
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO k = 1,llm
            DO j = ijb,ije
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestt )
            ENDDO
         ENDDO
c$OMP END DO NOWAIT
      ENDDO
      else
      ! general case, treat all tracers equally)
       DO iq = 1, nqtot
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO k = 1,llm
            DO j = ijb,ije
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestt )
            ENDDO
         ENDDO
c$OMP END DO NOWAIT
       ENDDO
      endif ! of if (planet_type=="earth")

c$OMP MASTER
      if (pole_nord) then
      
        DO  ij   = 1, iim
          xpn(ij) = aire(   ij   ) * pps(  ij     )
        ENDDO

        tpn      = SSUM(iim,xpn,1)/apoln

        DO ij   = 1, iip1
          pps (   ij     )  = tpn
        ENDDO
      
      endif

      if (pole_sud) then
      
        DO  ij   = 1, iim
          xps(ij) = aire(ij+ip1jm) * pps(ij+ip1jm )
        ENDDO

        tps      = SSUM(iim,xps,1)/apols

        DO ij   = 1, iip1
          pps ( ij+ip1jm )  = tps
        ENDDO
      
      endif
c$OMP END MASTER

      if (pole_nord) then
        DO iq = 1, nqtot
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO  k    = 1, llm
            DO  ij   = 1, iim
              xpn(ij) = aire(   ij   ) * pq(  ij    ,k,iq)
            ENDDO
            tpn      = SSUM(iim,xpn,1)/apoln
  
            DO ij   = 1, iip1
              pq (   ij   ,k,iq)  = tpn
            ENDDO
          ENDDO
c$OMP END DO NOWAIT	  
        ENDDO
      endif
      
      if (pole_sud) then
        DO iq = 1, nqtot
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO  k    = 1, llm
            DO  ij   = 1, iim
              xps(ij) = aire(ij+ip1jm) * pq(ij+ip1jm,k,iq)
            ENDDO
            tps      = SSUM(iim,xps,1)/apols
  
            DO ij   = 1, iip1
              pq (ij+ip1jm,k,iq)  = tps
            ENDDO
          ENDDO
c$OMP END DO NOWAIT	  
        ENDDO
      endif
      
      
      RETURN
      END
