      subroutine PHY_SISVAT_RUN(                                        &
!------------------------------------------------------------------------------+
!                                                                              |
!     SubRoutine PHY_SISVAT_RUN interfaces 3d-grid        Wed 26-Jun-2013 MAR  |
!               with SISVAT                2d-grid                             |
!                    SISVAT as  Soil                                           |
!                               Ice                                            |
!                               Snow                                           |
!                               Vegetation                                     |
!                               Atmosphere                                     |
!                               Transfer   Scheme                              |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Tue 26-Feb-2013      |
!           Last Modification by H. Gallee,               Wed 26-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+
! #TC&                        qxTC,uqTC  ,qsTC  ,                       &      ! Aerosols: Atm.Conc., Surf.Flux
     &                     )

!------------------------------------------------------------------------------+
!                                                                              |
!     INPUT:    it_RUN: Run Iterations Counter                                 |
!     ^^^^^                                                                    |
!                                                                              |
!     INPUT    (via module Mod_SISVAT_ctr)                                     |
!     ^^^^^     VegMod: SISVAT    is set up when .T.                           |
!               SnoMod: Snow Pack is set up when .T.                           |
!               inpSBC: Update Bound.Condit.when .T.                           |
!                                                                              |
!     INPUT    (via common block)                                              |
!     ^^^^^     xxxxTV: SISVAT/MAR interfacing variables                       |
!                                                                              |
!     Preprocessing  Option: SISVAT PHYSICS                                    |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^                                    |
!   #                        #SN: Snow         Model                           |
!   #                        #GP  LAI and GLF  Variations not specified        |
!   #                        #OP  SST       is interactive                     |
!                                                                              |
!   #                        #DS: diffuse radiation differing from direct      |
!                                (variable RADsod must still be included)      |
!                                                                              |
!     Preprocessing  Option: SISVAT PHYSICS: Col de Porte                      |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^                      |
!   #                        #CP: SBL,                       Col de Porte      |
!   #                        #cp  Solar Radiation,           Col de Porte      |
!   #                        #AG: Snow Ageing,               Col de Porte      |
!                                                                              |
!                                                                              |
!     Preprocessing  Option: STANDARD Possibility                              |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^                              |
!     #HY: Explicit Cloud MICROPHYSICS              may be turned ON           |
!     #BS: Explicit Cloud MICROPHYSICS: Blow. *(Snow)         Model)           |
!     #SI: SISVAT: Sea-Ice Fraction calculated from prescribed SST             |
!     #IP: SISVAT: Sea-Ice Fraction prescribed from SMMR and SSM/I             |
!     #SN: SNOW Model                               may be turned ON           |
!     #AO: COUPLING with  NEMO Ocean-Sea-Ice Model using OASIS                 |
!     #AE: TURBULENCE: Aerosols Erosion / Turbulent Diffusion Coeff.           |
!     #BD: TraCer   Aeolian Erosion  Submodel           is turned ON           |
!     #OR: SBL: Orography Roughness included from z0__SV (MARdom)              |
!     #SZ: SBL: Mom.: Roughn.Length= F(u*) Andreas &al.(1995)  Snow            |
!     #ZM: SBL: M/H   Roughn.Length: Box Moving Average (in Time)              |
!     #IB: OUTPUT: Ice-Sheet Surface Mass Balance  (on NetCDF File )           |
!     #vL: PORTABILITY: Vectorization enhanced                                 |
!     #vS: PORTABILITY: Vectorization enhanced: Snow    Model *                |
!     #vD: PORTABILITY: Vectorization enhanced: Blowing  Dust .                |
!     #vZ: PORTABILITY: Vectorization enhanced: Av.Roughness Length            |
!     #AA  TURBULENCE: SBL  Time Mean (BOX Moving Average)                     |
!     #AW  TURBULENCE: Wind Time Mean (BOX Moving Average)                     |
!     #AH  TURBULENCE: Ta-T Time Mean (BOX Moving Average)                     |
!                                                                              |
!                                                                              |
!     Preprocessing  Option: OBSOLETE                                          |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^                                          |
!     #AM  TURBULENCE: u*   Time Mean (BOX Moving Average)                     |
!     #AT  TURBULENCE: u*T* Time Mean (BOX Moving Average)                     |
!     #AS  TURBULENCE: u*s* Time Mean (BOX Moving Average)                     |
!                                                                              |
!                                                                              |
!     Preprocessing  Option:                                                   |
!     ^^^^^^^^^^^^^^^^^^^^^                                                    |
!     #TC: TraCer   Advection-Diffusion Equation        is turned ON           |
!     #CM: Update Green Leaf Fraction may be performed (in routine INIglf)     |
!     #VM: TURBULENCE: Ua   minimum is prescribed                              |
!     #GP: Soil /Vegetation Model: LAI, GLF Variations NOT prescrib.           |
!     #OP: SISVAT: Interactive Sea Surface Temperature     turned ON           |
!     #op: SISVAT: SST Nudging -->   prescribed values     turned ON           |
!     #RE: SISVAT: SST Condition     modified   by         epsilon             |
!     #PO: POLYNYA Model                            may be turned ON           |
!                                                                              |
!                                                                              |
!     Preprocessing  Option: SISVAT IO (not always a standard preprocess.)     |
!     ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^                                         |
!     FILE                 |      CONTENT                                      |
!     ~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     |
!   # SISVAT_iii_jjj_n     | #e0: OUTPUT on ASCII  File (SISVAT Variables)     |
!   #                      |(#e0  MUST BE PREPROCESSED BEFORE #e1 & #e2 !)     |
!   # SISVAT_iii_jjj_n     | #e1: OUTPUT/Verification: Energy Conservation     |
!   # SISVAT_iii_jjj_n     | #e2: OUTPUT/Verification: Energy Consrv.2e pt.    |
!                          |                           (no premature stop)     |
!   # SISVAT_iii_jjj_n     | #e3: OUTPUT/Verification: Energy Consrv. (HL)     |
!   # SISVAT_iii_jjj_n     | #e4: OUTPUT/Verification: Energy Consrv. (HLS)    |
!   # SISVAT_iii_jjj_n     | #el: OUTPUT Format      : Addition.Frame Line     |
!                          |                                                   |
!   # SISVAT_iii_jjj_n     | #m0: OUTPUT/Verification: H2O    Conservation     |
!   # SISVAT_iii_jjj_n     | #m1: OUTPUT/Verification: * Mass Conservation     |
!   # SISVAT_iii_jjj_n     | #m2: OUTPUT/Verification: SeaIce Conservation     |
!   # stdout               | #mw: OUTPUT/Verification: H2O    Conservation     |
!                          |      unit  6, SubRoutine  SISVAT_qSo **ONLY**     |
!                          |                                                   |
!   # SISVAT_iii_jjj_n     | #mb: OUTPUT/Verification: Blowing Snow            |
!                          |               SubRoutine  PHY_SISVAT **ONLY**     |
!                          |                                                   |
!   # SISVAT_zSn.vz        | #vz: OUTPUT/Verification: Snow Layers Agrega.     |
!                          |      unit 41, SubRoutine  SISVAT_zSn **ONLY**     |
!   # SISVAT_qSo.vw        | #vw: OUTPUT/Verif+Detail: H2O    Conservation     |
!                          |      unit 42, SubRoutine  SISVAT_qSo **ONLY**     |
!   # SISVAT_qSn.vm        | #vm: OUTPUT/Verification: Energy/Water Budget     |
!                          |      unit 43, SubRoutine  SISVAT_qSn **ONLY**     |
!   # SISVAT_qSn.vu        | #vu: OUTPUT/Verification: Slush  Parameteriz.     |
!                          |      unit 44, SubRoutine  SISVAT_qSn **ONLY**     |
!   # SISVAT_wEq.ve        | #ve: OUTPUT/Verification: Snow/Ice Water Eqv.     |
!                          |      unit 45, SubRoutine  SISVAT_wEq **ONLY**     |
!   # SnOptP____.va        | #va: OUTPUT/Verification: Albedo Parameteriz.     |
!                          |      unit 46, SubRoutine  SnOptP     **ONLY**     |
!   # SISVAT_GSn.vp        | #vp: OUTPUT/Verification: Snow   Properties       |
!                          |      unit 47, SubRoutines SISVAT_zSn, _GSn        |
!   # PHY_SISVAT.v0        | #v0: OUTPUT/Verification: DUMP                    |
!                          |      unit 50, SubRoutine  PHY_SISVAT **ONLY**     |
!                          |                                                   |
!   # stdout               | #b0: OUTPUT of Snow Erosion Variables             |
!                          |      unit  6, SubRoutine  SISVAT_BSn **ONLY**     |
!   # stdout               | #b1: OUTPUT of Snow Erosion Turbulence            |
!                          |      unit  6, SubRoutine  SISVATeSBL **ONLY**     |
!   # stdout               | #b2: OUTPUT of Snow Erosion Turbulence (2)        |
!                          |      unit  6, SubRoutine  SISVATeSBL **ONLY**     |
!                          |                                                   |
!   # stdout               | #s0: OUTPUT of Snow Buffer Layer                  |
!                          |      unit  6, SubRoutine  SISVAT     **ONLY**     |
!   # stdout               | #s1: OUTPUT of Snow Layers Agregation             |
!                          |      unit  6, SubRoutine  SISVAT_zSn, _zAg        |
!   # stdout               | #s2: OUTPUT of SnowFall, Snow Buffer              |
!                          |      unit  6, SubRoutine  SISVAT_BSn, _qSn        |
!   # stdout               | #sb: OUTPUT/Verification: SISVAT_SBL              |
!                          |      unit  6, SubRoutine  SISVAT_SBL **ONLY**     |
!   # stdout               | #sd: OUTPUT/Verification: Snow Thinest Layer      |
!                          |      unit  6, SubRoutine  SISVAT_zSn **ONLY**     |
!   # stdout               | #sf: OUTPUT of SnowFall, Z0 and Drag Coeff.       |
!                          |      unit  6, SubRoutines PHY_SISVAT, SISVAT      |
!   # stdout               | #sg: OUTPUT/Verification: Gravitational Front     |
!                          |      unit  6, SubRoutine  SISVAT_qSo **ONLY**     |
!   # stdout               | #sh: OUTPUT/Verification: SBL Turbulent Fluxes    |
!                          |      unit  6, SubRoutine  PHY_SISVAT **ONLY**     |
!   # stdout               | #si: OUTPUT of Sea Ice Fraction, Temperature      |
!                          |      unit  6, SubRoutine  PHY_SISVAT **ONLY**     |
!   # stdout               | #sm: OUTPUT/Verification: Mosaic Fractions Sum    |
!                          |      unit  6, SubRoutine  PHY_SISVAT **ONLY**     |
!   # stdout               | #so: OUTPUT/Verification: Soil Vertic.Discret.    |
!                          |      unit  6, SubRoutine  SISVAT_ini **ONLY**     |
!   # stdout               | #ss: OUTPUT of Blowing Snow Variables             |
!                          |      unit  6, SubRoutine  SISVATeSBL **ONLY**     |
!   # stdout               | #sz: OUTPUT of Roughness Length & Drag Coeff.     |
!                          |      unit  6, SubRoutine  SISVAT     **ONLY**     |
!                          |                                                   |
!------------------------------------------------------------------------------+


! Global Variables
! ================

      use Mod_Real
      use Mod_PHY____dat
      use Mod_PHY____grd
      use Mod_SISVAT_dim
      use Mod_SISVAT_grd
      use Mod_SISVAT_ctr
      use Mod_SISVAT_gpt

      use Mod_PHY_S0_dat
      use Mod_PHY_S0_kkl
      use Mod_PHY_RT_kkl
      use Mod_PHY_DY_kkl
      use Mod_PHY_AT_kkl
      use Mod_PHY_CM_kkl


! INTERFACE Variables
! ===================

      use Mod_SISVAT_dat
      use Mod_SISVAT_dzS
      use Mod_SISVAT_kkl
      use Mod_SISVAT_loc
! #AW use Mod_SISVAT_xAW
! #AH use Mod_SISVAT_xAH
! #ZM use Mod_SISVAT_xZM
! #WL use Mod_SISVAT_xWL



      IMPLICIT NONE



      INTEGER    ntrac
      PARAMETER (ntrac=1)


! Tracers
! ^^^^^^^
! #TC real       uqTC  ( mx, my,    ntrac)  ! q*     (Tracer)
! #TC real       qsTC  ( mx, my,    ntrac)  ! qSurf  (Tracer)
! #TC real       qxTC  ( mx, my, mz,ntrac)  ! q      (Tracer)

! Auxiliary Variables
! ^^^^^^^^^^^^^^^^^^^
! #CM integer                                          ::  newglfSIS    !
! #IP integer                                          ::  newsicSI     !




! Internal  Variables
! ===================

      logical            glfFIX

      character(len=1)   cha  ,chb
      integer            jjtime
      integer            iwr
! #SI integer            l
      integer            isl  ,isn
      integer            i    ,j    ,k     ,ikp
      integer            n    ,nm

      real(kind=real8)   rhAir                ! Air     Densitity
      real(kind=real8)   Ua_min               ! Minimum Air Velocity
      real(kind=real8)   d_snow,SnowOK        ! Snow Precip.: Total
! #BS real(kind=real8)   dbsnow               ! Snow Precip.: from Drift
! #WL real(kind=real8)   uqstar               ! u*q*

! #SI real(kind=real8)   SrfSIC,SIc0OK        ! Oceanic Fraction: previous
! #SI real(kind=real8)   FraOcn,SIceOK        ! Oceanic Fraction

! #Za real(kind=real8)   S_Eros               ! Snow Erosion (status) = (1,0)
! #Za real(kind=real8)   SnEros               ! Snow Erosion (status) = (1,0)

! OUTPUT/Verification:   Blowing Snow
! #mb integer            noUNIT
! #mb real(kind=real8)   BlowST,SnowSB


! DATA
! ====

      data     glfFIX   /  .false./
      data     cha      /'-'/
      data     chb      /':'/



! SISVAT Time
! ===========

          jjtime = HourTU*3600+minuTU*60+Sec_TU




! SISVAT Variables:      (x,y,z) Domain
! =====================================

      DO ikp=1,kcolp
      DO nm =1,mwp


! SISVAT Variables: Interpolation of V at 10 m a.g.l.
! ---------------------------------------------------

        IF (k_SL.EQ.mzp)                                            THEN
          VV10SV(ikp,nm) = r_SL10 * WindSV(ikp,nm,k_SL)                 !
        ELSE
          VV10SV(ikp,nm) =          WindSV(ikp,nm,k_SL)                &!
     &                   + r_SL10 *(WindSV(ikp,nm,k_SL-1)              &!
     &                             -WindSV(ikp,nm,k_SL))                !
        END IF


! Sastrugi Height decreased by Precipitation if V < 6 m/s (Kotlyakov, 1961)
! -------------------------------------------------------------------------

! #ZS   dz0_SV(ikp,nm)  =      max(zer0,(SnowCM(ikp)-Sno0CM(ikp))      &!
! #ZS&    /max(0.05,0.104*sqrt(max(zer0, WindSV(ikp,nm,mzp)-6.00))))    !
! #ZS   dz0_SV(ikp,nm)  =                dz0_SV(ikp,nm) *0.01*max(2-n,0)! dz0(Sastrugi dh)


! Influence of the Angle(Wind,Sastrugi) (Andreas, 1995, CCREL report 95-16)
! -------------------------------------------------------------------------

! #Za   S_Eros     = max(zer0,sign(un_1,-uss_CM(ikp) -eps9))
! #Za   SnEros     = max(zer0,sign(un_1, uss_CM(ikp) +eps9))
! #Za   VVs_SV(ikp,nm)    = VVs_SV(ikp,nm) * RRsxSV(ikp,nm)
! #Za   VVs_SV(ikp,nm)    =                                            &
! #Za&           SnEros* VVs_SV(ikp,nm)                                &
! #Za&         + S_Eros*(VVs_SV(ikp,nm) * Adz0dt     + VV10SV(ikp,nm))
! #Za   RRsxSV(ikp,nm)    =                                            &
! #Za&           SnEros* RRsxSV(ikp,nm)                                &
! #Za&         + S_Eros*(RRsxSV(ikp,nm) * Adz0dt     + 1.0 )
! #Za   DDsxSV(ikp,nm)    =                                            &
! #Za&           SnEros* DDsxSV(ikp,nm)                                &
! #Za&         + S_Eros* DDsxSV(ikp,nm) * Adz0dt                       &
! #Za&         +        (Va__SV(ikp,nm) *(Ua__SV(ikp,nm) -Ua0_SV(ikp,nm))  &
! #Za&                  -Ua__SV(ikp,nm) *(Va__SV(ikp,nm) -Va0_SV(ikp,nm))) &
! #Za&   /(degrad*max(.3,WindSV(ikp,nm,mzp) * WindSV(ikp,nm,mzp)))
! #Za    IF (DDsxSV(ikp,nm)   .GT.360.) DDsxSV(ikp,nm)  = DDsxSV(ikp,nm) - 360.
! #Za    IF (DDsxSV(ikp,nm)   .LT.  0.) DDsxSV(ikp,nm)  = DDsxSV(ikp,nm) + 360.
      END DO
      END DO


! Water Vapor Flux Limitor
! ------------------------

! #WL DO ikp=1,kcolp
! #WL DO nm =1,mwp

! #WL  DO n=1,nLimi-1
! #WL    WL_mem(ikp,nm,n) = WL_mem(ikp,nm,n+1)
! #WL  END DO

! #WL    uqstar  =  max(abs(uqs_SV_gpt(ikp),eps6)*sign(1.,uqs_SV_gpt(ikp))
! #WL    WL_mem(ikp,nm,n)                                              &
! #WL& = Kzh_AT(ikp,nm,k1m(mzp))
! #WL&    *(qv__SV(ikp,nm,mzpp-k2m(mzp))-qv__SV(ikp,nm,mzpp-k1m(mzp))) &
! #WL&    /(uqstar    *(hsigma(k2m(mzp))            -hsigma(k1m(mzp))))

! #WL     WLmmem(ikp,nm) = 0.
! #WL  DO n=1,nLimi
! #WL     WLmmem(ikp,nm) = WLmmem(ikp,nm)+WL_mem(ikp,nm,n)
! #WL  ENDDO
! #WL     WLmmem(ikp,nm) = WLmmem(ikp,nm)/nLimit

! #WL END DO
! #WL END DO



! SWD / SWA INPUT
! ---------------

      DO ikp = 1,kcolp

      IF (.NOT.InpSWD)                                              THEN
        SWDsRT(ikp) = SWAsRT(ikp)/(1.d0-Alb_SV_gpt(ikp))                ! => downward Solar
      END IF


! Grid Averages
! -------------

        Alb_SV_gpt(ikp)  = 0.
        EmisSV_gpt(ikp)  = 0.
        LWUsRT    (ikp)  = 0.
        LMO_SV_gpt(ikp)  = 0.
        us__SV_gpt(ikp)  = 0.
        uts_SV_gpt(ikp)  = 0.
        uqs_SV_gpt(ikp)  = 0.
        uss_CM(ikp)      = 0.
        qs__CM(ikp,mzpp) = 0.
        hFraSV_gpt(ikp)  = 0.
        Tas_SV_gpt(ikp)  = 0.
        qv__DY(ikp,mzpp) = 0.
! #TC     uqTC(i,j,1)    = 0.
! #TC     qsTC(i,j,1)    = 0.
      END DO



! SISVAT Variables: from (x,y,z) Domain to (vector,mosaic) Domain
! ===============================================================

! SISVAT Variables Update
! ^^^^^^^^^^^^^^^^^^^^^^^
      DO      ikp         = 1,kcolp
      DO      nm          = 1,mwp
              i           = ii__AP(ikp)
              j           = jj__AP(ikp)

! OUTPUT of SnowFall, Roughness Length and Drag Coefficients
! #sf         IF (ikp.EQ.1 .AND. Sec_TU.EQ.0) 
! #sf&        write(6,6659) 
! #sf 6659    format(20x,'   dsn_SV   us__SV   Z0SaSi   Z0Sa_N'        &
! #sf&                  ,'   Z0SaSV   Z0m_Sn   Z0m_SV')

! Atmospheric Forcing                                       (INPUT)
! ^^^^^^^^^^^^^^^^^^^                                        ^^^^^
            DO isn=1,mzp
              isl =  mzp + 1  - isn                                     !
              Kz__SV(ikp,nm,isn) = Kzh_AT(ikp,isl)                      !
              pktaSV(ikp,nm,isn) = pkt0SV(ikp,nm,isn)                   !
            ENDDO
              za__SV(ikp,nm) =    zza_SV(ikp,nm,1)                      ! [m]
                                                                        !
              VV__SV(ikp,nm) =    WindSV(ikp,nm,mzp)                    !

              Ua_min         = eps6                                     !
! #VM         Ua_min         = 0.2 * sqrt(za__SV(ikp,nm) )              !
              VV__SV(ikp,nm) = max(Ua_min,VV__SV(ikp,nm) )              !
                                                                        !
! #Za         VVs_SV(ikp,nm) = VVs_SV(ikp,nm) / RRsxSV(ikp,nm)          !
! #Za         DDs_SV(ikp,nm) =         max(zer0,DDsxSV(ikp,nm)-180.)   &
! #Za&        +180.*min(un_1,zer0-min(zer0,DDsxSV(ikp,nm)     -180.))  &
! #Za&                           +min(zer0,DDsxSV(ikp,nm)     -180.)
              rhT_SV(ikp,nm) = pkPaSV(ikp,nm)   *1.e3                  &! [kg/m3]
     &                       /(TaT_SV(ikp,nm)   *R_DAir)                !
              QaT_SV(ikp,nm) = qv__SV(ikp,nm,1)                         !
! #WL         dQa_SV(ikp,nm) =  max(0.,1.-WLmmem(ikp,nm))              &! Water  Vapor
! #WL&                        * dt__SV/hsigma(mzp)                      ! Flux Limitor
              qsnoSV(ikp,nm) =      0.                                 &!
     &                       + qs__CM(ikp,mzp)                         &!
! #TC&                       +   qxTC(i,j,mzp,1)                       &!
     &                       +   0.

! Energy Fluxes                                             (INPUT)
! ^^^^^^^^^^^^^                                              ^^^^^
              coszSV(ikp,nm) = max(cszEPS,csz_S0(ikp))                  ! cos(zenith.Dist.)
              sol_SV(ikp,nm) =            SWDsRT(ikp)                   !    downward Solar
              IRd_SV(ikp,nm) =            LWDsRT(ikp)                   !    downward IR

! Water  Fluxes                                             (INPUT)
! ^^^^^^^^^^^^^                                              ^^^^^
              drr_SV(ikp,nm) =(RainCM(ikp)-Rai0CM(ikp)) *1.e3 /dt__SV   ! [m/s] -> [mm/s] .OR. [kg/m2/s]
              d_snow         = SnowCM(ikp)-SnobCM(ikp)                  ! Only SnowFall
              dsn_SV(ikp,nm) = d_snow                   *1.e3 /dt__SV   ! Erosion NOT incl.
              SnowOK         =                                         &! Correction
     &         max(zer0,sign(un_1,       qs__CM(ikp,mzp)-eps6))        &!
     &        *max(zer0,sign(un_1,Tf_Sno-TaT_SV(ikp,nm) -eps6))         !
              dsn_SV(ikp,nm) = dsn_SV(ikp,nm)+drr_SV(ikp,nm)*SnowOK     !
              drr_SV(ikp,nm) = drr_SV(ikp,nm)           *(1.-SnowOK)    !
! #BS         dbsnow         =-ussxSV(ikp,nm)                          &! Erosion
! #BS&                        *dt__SV     *rhT_SV(ikp,nm)               !
! #bS         dsnbSV(ikp,nm) =min(max(zer0,dbsnow)                     &!
! #bS&                    /       max(eps6,d_snow),un_1)                !
! #BS         dsnbSV(ikp,nm) =1.0-min(qs__CM(ikp,k_zb)                 &! k_zb level ~ 25 magl
! #BS&                           /max(qs__CM(ikp,mzp ),epsn),un_1)      !(BS negligib.at k_zb)
!             dsnbSV is used and modified in SISVAT_BSn,
!                  then used for Buffer Layer Update
! #BS         dbs_SV(ikp,nm) =        ussbSV(ikp,nm)                    !          [kg/m2]

! Soil/Canopy                                               (INPUT)
! ^^^^^^^^^^^                                                ^^^^^
              slorSV(ikp,nm) = atan(slopSV(ikp,nm))                     ! Fall Line Slope       [radian]

! Energy Fluxes                                             (INPUT/OUTPUT)
! ^^^^^^^^^^^^^                                              ^^^^^^^^^^^^
              cld_SV(ikp,nm) =      ClouRT(ikp)                         ! Cloudiness

! Water  Fluxes                                             (INPUT/OUTPUT)
! ^^^^^^^^^^^^^                                              ^^^^^^^^^^^^
! #BS         uss_SV(ikp,nm) =      ussxSV(ikp,nm)                      ! u*s*

! #vL ENDDO
! #vL ENDDO

! Snow Roughness                                            (INPUT/OUTPUT)
! ^^^^^^^^^^^^^^                                             ^^^^^^^^^^^^
! #vL DO      ikp = 1,kcolp
! #vL DO      nm  = 1,mwp

! #ZM         Z0mmSV(ikp,nm) =     0.                                   !
! #ZM         Z0emSV(ikp,nm) =     0.                                   !
! #ZM         Z0hmSV(ikp,nm) =     0.                                   !
! #ZM       DO k  =   1,ntavz                                           !
! #ZM         Z0mmSV(ikp,nm) =     Z0mmSV(ikp,nm)                      &!
! #ZM&                       +     z0_mem(ikp,nm,k)                     !
! #ZM         Z0emSV(ikp,nm) =     Z0emSV(ikp,nm)                      &!
! #ZM&                       +     b0_mem(ikp,nm,k)                     !
! #ZM         Z0hmSV(ikp,nm) =     Z0hmSV(ikp,nm)                       !
! #ZM&                       +     r0_mem(ikp,nm,k)                     !
! #ZM       ENDDO                                                       !
! #ZM         Z0mmSV(ikp,nm) = min(Z0mmSV(ikp,nm)  /ntavz              &!  z0(Mom., Box Av.)
! #ZM&                            ,hsigma(mzp)     /    3.)             !
! #ZM         Z0emSV(ikp,nm) =     Z0emSV(ikp,nm)  /ntavz               !  b0(Eros, Box Av.)
! #ZM         Z0hmSV(ikp,nm) =     Z0hmSV(ikp,nm)  /ntavz               !  r0(Heat, Box Av.)

! V,  dT(a-s)    Time Moving Averages
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #AA       DO k =  1,ntaver-1
! #AW         V__mem(ikp,nm,k    ) = V__mem(ikp,nm,k+1)
! #AH         T__mem(ikp,nm,k    ) = T__mem(ikp,nm,k+1)
! #AA       ENDDO
! #AW         V__mem(ikp,nm,ntave) = VV__SV(ikp,nm)
! #AH         T__mem(ikp,nm,ntave) = TaT_SV(ikp,nm)-Ts__SV(ikp,nm)

! #AW         VVmmem(ikp,nm)       = 0.0
! #AH         dTmmem(ikp,nm)       = 0.0
! #AA       DO k=1,ntaver
! #AW         VVmmem(ikp,nm)       = VVmmem(ikp,nm)+V__mem(ikp,nm,k)
! #AH         dTmmem(ikp,nm)       = dTmmem(ikp,nm)+T__mem(ikp,nm,k)
! #AA       ENDDO
! #AW         VVmmem(ikp,nm)        = VVmmem(ikp,nm)/ntaver
! #AH         dTmmem(ikp,nm)        = dTmmem(ikp,nm)/ntaver

! #vL ENDDO
! #vL ENDDO

! Snow Pack                                                 (INPUT/OUTPUT)
! ^^^^^^^^^                                                  ^^^^^^^^^^^^
! #vS DO      ikp = 1,kcolp
! #vS DO      nm  = 1,mwp
              dsn_SV(ikp,nm) =     dsn_SV(ikp,nm)                      &!
     &                        +max(BufsSV(ikp,nm)-SMndSV,0.)           &!
     &                        /    dt__SV                               !
              BufsSV(ikp,nm) = min(BufsSV(ikp,nm),SMndSV   )            !
! #vS ENDDO
! #vS ENDDO

      DO      isn = 1,nsnow
! #vS DO      ikp = 1,kcolp
! #vS DO      nm  = 1,mwp
              G1snSV(ikp,nm,isn) = max(-G1_dSV,min(G1_dSV,             &!
     &                                  G1snSV(ikp,nm,isn)))            ! [-]        [-]
              G2snSV(ikp,nm,isn) = max(-G1_dSV,min(G1_dSV,             &!
     &                                  G2snSV(ikp,nm,isn)))            ! [-] [0.0001 m]
! #vS END DO
! #vS END DO
      END DO

! #vL DO      ikp = 1,kcolp
! #vL DO      nm  = 1,mwp
! #vL         i   = ii__AP(ikp)   
! #vL         j   = jj__AP(ikp)   

! #SI         HFraSV(ikp,nm)     = 0.                         ! Frazil Thickness

! RunOFF Intensity                                          (INPUT/OUTPUT)
! ^^^^^^^^^^^^^^^^                                           ^^^^^^^^^^^^
              RnofSV(ikp,nm)     = 0.                         ! RunOFF Intensity

! Blown Snow OUTPUT for a specific Grid Point
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #mb    IF  (lwriSV(ikp,nm).ne.0.AND.it_RUN.gt.0)                     THEN
! #mb         noUNIT =   no__SV(lwriSV(ikp,nm))
! #mb         write(noUNIT,5012)
! #mb 5012    format(/,1x)
! #mb         write(noUNIT,5013)
! #mb 5013    format(' -----+--------+--------+--------+--------+',    &
! #mb&             '--------+--------+--------+--------+--------+',    &
! #mb&             '--------+')
! #mb         write(noUNIT,5014)
! #mb 5014    format('    n |     z  |     qs |      V |        |',    &
! #mb&             '     T  | TKE^0.5|        |        |        |',    &
! #mb&             '        |',                                        &
! #mb&             /,'      |    [m] | [g/kg] |  [m/s] |        |',    &
! #mb&             '    [K] |  [m/s] |        |        |        |',    &
! #mb&             '        |')
! #mb           BlowST=0.
! #mb           k=0
! #mb 5011    CONTINUE
! #mb           k=k+1
! #mb         IF (               k             .gt.mzp )      GO TO 5010
! #mb         IF (Z___DY    (ikp,k)-sh__AP(ikp).lt.100.)            THEN
! #mb           BlowST=BlowST+WindSV(ikp,nm,k)*qs__CM    (ikp,k)       &
! #mb&                 *1.e3*(pkPaSV(ikp,nm)  -pt__DY)*d1Sigm(k) *Grav_I
! #mb           write(noUNIT,5015) mzpp-k                              &
! #mb&           ,     Z___DY(ikp,k)   -sh__AP(ikp)                    &
! #mb&           ,     qs__CM(ikp,k)   *1.e3                           &
! #mb&           ,     WindSV(ikp,nm,k),Ta__DY(ikp,k)                  &
! #mb&           ,sqrt(TKE_AT(ikp,k))          
! #mb 5015      format(i5,' |',f7.2,' |',f7.3,' |',f7.2,' |',          &
! #mb&                  8x,'|',f7.2,' |',f7.3,' |', 4(8x,'|'))
! #mb         END IF
! #mb         GO TO 5011
! #mb 5010    CONTINUE

! #mb           SnowSB =   BufsSV(ikp,nm)      
! #mb         IF        (  isnoSV(ikp,nm)   .GT.     0      )       THEN
! #mb         DO isn=max(0,isnoSV(ikp,nm)  ) ,isnoSV(ikp,nm)  
! #mb           SnowSB =   SnowSB                                      &
! #mb&                +    dzsnSV(ikp,nm,isn)*ro__SV(ikp,nm,isn)
! #mb         END DO
! #mb         END IF

! #mb         write(noUNIT,5016) BlowST,SnowSB
! #mb 5016    format(' * TRANSPORT = ',e12.3,' kg/m/s'  ,  8x,'|',     &
! #mb&               ' * BUDGET    = ',f12.6,' mm w.e.|',2(8x,'|'))
! #mb         write(noUNIT,5013)
! #mb    END IF

! OUTPUT, for Stand-Alone VERIFICATION
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! #v0    IF  (ikp .le. kcolp .and. nm .le. mxpp)                    THEN
! #v0         write(50,5001) it_RUN,ikp,nm                             &
! #v0&                   za__SV(ikp,nm),VV__SV(ikp,nm),TaT_SV(ikp,nm), &
! #v0&                   rhT_SV(ikp,nm),QaT_SV(ikp,nm),qsnoSV(ikp,nm), &
! #v0&                   coszSV(ikp,nm),sol_SV(ikp,nm),IRd_SV(ikp,nm), &
! #v0&                   drr_SV(ikp,nm),dsn_SV(ikp,nm),dbs_SV(ikp,nm), &
! #v0&                   LSmask(ikp,nm),isotSV(ikp,nm),alb0SV(ikp,nm), &
! #v0&                   IRs_SV(ikp,nm),                               &
! #v0&                   ivgtSV(ikp,nm),LAI0SV(ikp,nm),glf0SV(ikp,nm), &
! #v0&                   TvegSV(ikp,nm),LMO_SV(ikp,nm),us__SV(ikp,nm), &
! #v0&                   uqs_SV(ikp,nm),uts_SV(ikp,nm),uss_SV(ikp,nm), &
! #v0&                   snCaSV(ikp,nm),rrCaSV(ikp,nm),psivSV(ikp,nm)
! #v0 5001    format(/,'c #INFO   it_RUN             = ',i15  ,        &
! #v0&               /,'c #INFO   ikp                = ',i15  ,        &
! #v0&               /,'c #INFO   nm                 = ',i15  ,        &
! #v0&               /,'          za__SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          VV__SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          TaT_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          rhT_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          QaT_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          qsnoSV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          coszSV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          sol_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          IRd_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          drr_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          dsn_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          dbs_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          LSmask(ikp,nm)     = ',i15  ,        &
! #v0&               /,'          isotSV(ikp,nm)     = ',i15  ,        &
! #v0&               /,'          alb0SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          IRs_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          ivgtSV(ikp,nm)     = ',i15  ,        &
! #v0&               /,'          LAI0SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          glf0SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          TvegSV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          LMO_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          us__SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          uqs_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          uts_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          uss_SV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          snCaSV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          rrCaSV(ikp,nm)     = ',e15.6,        &
! #v0&               /,'          psivSV(ikp,nm)     = ',e15.6)
! #v0         DO isl = -nsoil,0
! #v0           write(50,5002)    isl,TsisSV(ikp,nm,isl),              &
! #v0&                            isl,eta_SV(ikp,nm,isl)
! #v0 5002      format('          TsisSV(ikp,nm,',i2,')  = ',e15.6,    &
! #v0&                 '          eta_SV(ikp,nm,',i2,')  = ',e15.6)
! #v0         END DO
! #v0         DO isl =       1,nsnow
! #v0           write(50,5003)    isl,TsisSV(ikp,nm,isl),              &
! #v0&                            isl,dzsnSV(ikp,nm,isl)
! #v0 5003      format('          TsisSV(ikp,nm,',i2,')  = ',e15.6,    &
! #v0&                 '          dzsnSV(ikp,nm,',i2,')  = ',e15.6)
! #v0         END DO
! #v0    END IF
      END DO
      END DO

! SISVAT Execution
! ^^^^^^^^^^^^^^^^
          write(daHost,'(i2,a3,i4,i3,2(a1,i2))')                        &
     &          Day_TU,LabMon(Mon_TU),YearTU,                           &
     &          HourTU,chb,minuTU,chb,Sec_TU

! OUTPUT of SnowFall, Roughness Length and Drag Coefficients
! #sf     write(6,666) Day_TU,Mon_TU,YearTU,HourTU,minuTU,Sec_TU
! #sf 666 format(2(i2,'-'),2i4,2(':',i2),3x,$)



!           ******
      call  SISVAT(jjtime,kcolv)
!           ******



! MAR    Variables Update
! ^^^^^^^^^^^^^^^^^^^^^^^
      DO      ikp = 1,kcolp
      DO      k   =    1,mzp 
              dpktSV_gpt(ikp,k) = 0.
      ENDDO
      ENDDO

      DO      ikp = 1,kcolp
      DO      nm  = 1,mwp

! Atmospheric Forcing                                       (OUTPUT)
! ^^^^^^^^^^^^^^^^^^^                                        ^^^^^^
            DO isn=1,mzp
              k   =  mzp + 1  - isn                                     !
              dpktSV_gpt(ikp,k)     =                                  &!
      &       dpktSV_gpt(ikp,k)     +                                  &!
      &       FracSV(ikp,nm)        *                                  &!
      &      (pktaSV(ikp,nm,isn)-pkt0SV(ikp,nm,isn)) / dt__SV           !
            ENDDO

! Water  Fluxes                                             (INPUT/OUTPUT)
! ^^^^^^^^^^^^^                                              ^^^^^^^^^^^^
! #BS     ussxSV(ikp,nm)         = 0.                                   ! 
! #BS     ussxSV(ikp,nm)         =                                     &! MAX        Erosion 
! #BS&   (ussbSV(ikp,nm)  -        dbs_SV(ikp,nm))                     &!-unconsumed Erosion
! #BS&                    /(dt__SV*rhT_SV(ikp,nm))                      ! ==> actual u*s*
! #vL ENDDO
! #vL ENDDO

! #vL DO      ikp = 1,kcolp
! #vL DO      nm  = 1,mwp  
! #BS     ussbSV(ikp,nm)  = dt__SV*uss_SV(ikp,nm)                      &! NEW  MAX   Erosion
! #BS&                            *rhT_SV(ikp,nm)                       ! rho u*s* dt[kg/m2]

! Soil/Canopy                                               (INPUT/OUTPUT)
! ^^^^^^^^^^^                                                ^^^^^^^^^^^^
          Z0m_SV(ikp,nm)  =    min(Z0m_SV(ikp,nm)                      &! Moment.Roughn.L.
      &                           ,hsigma(mzp)/3.)                      !
! #vL END DO
! #vL END DO

! Snow Roughness                                            (INPUT/OUTPUT)
! ^^^^^^^^^^^^^^                                             ^^^^^^^^^^^^
! #vZ DO      ikp = 1,kcolp
! #vZ DO      nm  = 1,mwp  
! #ZM       DO k=1,ntavz-1                                              !
! #ZM         z0_mem(ikp,nm,k)       = z0_mem(ikp,nm,k+1)               !
! #ZM         b0_mem(ikp,nm,k)       = b0_mem(ikp,nm,k+1)               !
! #ZM         r0_mem(ikp,nm,k)       = r0_mem(ikp,nm,k+1)               !
! #ZM       ENDDO                                                       !
! #vZ ENDDO
! #vZ ENDDO

! #vL DO      ikp = 1,kcolp
! #vL DO      nm  = 1,mwp  
! #ZM         z0_mem(ikp,nm,ntavz)   = Z0mnSV(ikp,nm)                   ! z0(Momentum)
! #ZM         b0_mem(ikp,nm,ntavz)   = Z0enSV(ikp,nm)                   ! z0(Mom., Erosion)
! #ZM         r0_mem(ikp,nm,ntavz)   = Z0hnSV(ikp,nm)                   ! z0(Heat)
! #vL END DO
! #vL END DO

! Dust   Fluxes                                             (INPUT/OUTPUT)
! ^^^^^^^^^^^^^
! #vD DO      ikp = 1,kcolp
! #vD DO      nm  = 1,mwp  
! #BD         uds_SV(ikp,nm)         =(1-min(1,isnoSV(ikp,nm)))        &! Snow Free  Surface 
! #BD&                                        *uss_SV(ikp,nm)          &! DUST       Erosion
! #BD&                                  *max(1,(2-my)*3)                ! Tuning Factor (2D)
! #vD END DO
! #vD END DO

! Snow Pack                                                 (INPUT/OUTPUT)
! ^^^^^^^^^                                                  ^^^^^^^^^^^^
! #vL DO      ikp = 1,kcolp
! #vL DO      nm  = 1,mwp  
! #IB         wes0SV(ikp,nm)         = dwesSV(ikp,nm)                   ! Depo. / Subli.
! #IB         wem0SV(ikp,nm)         = dwemSV(ikp,nm)                   ! Melting
! #IB         wer0SV(ikp,nm)         = dwerSV(ikp,nm)                   ! Refreezing

! Radiative Properties                                            (OUTPUT)
! ^^^^^^^^^^^^^^^^^^^^                                             ^^^^^^
          albcSV(ikp,nm)             = alb_SV(ikp,nm)                  &! Mosaic Albedo
! #AO&                             *   LSmask(ikp,nm)                  &!
! #AO&                             +   Alb_AO(ikp,nm)                  &! Mosaic AlbedoNEMO
! #AO&                             *(1-LSmask(ikp,nm))                 &!
     &                             + 0.
          ROF_SV(ikp,nm)             = RnofSV(ikp,nm)   *dt__SV        &!
     &                             +   ROF_SV(ikp,nm)                   ! Integrated Run OFF

      END DO
      END DO



! Surface Temperature: Prescription of relevant Medium (Snow, precribed SST)
! ==========================================================================

! Control Martin
!PRINT*, 'nsoil=',nsoil
!PRINT*, 'FixSST=',FixSST
!PRINT*, 'LSMask=',LSMask
! Control Martin

      DO ikp=1,kcolp
            DO isl =   -nsoil,0                                          !


! Open Ocean
! ----------

              eta_SV(ikp,1,isl) =                                      &!
     &        eta_SV(ikp,1,isl) *   LSMask(ikp,1)                      &!
     &                          +(1-LSMask(ikp,1)    )                  ! Sea: Humidity:=1

              TsisSV(ikp,1,isl) =                                      &!
     &       (TsisSV(ikp,1,isl) *   LSMask(ikp,1)                      &! Soil Temperature
     &       +sst_SB(ikp)       *(1-LSMask(ikp,1)    )) * FixSST       &! Prescribed   SST
! #OP&                                                  * FixSST       &!
! #OP&      +(TsisSV(ikp,1,isl)                                        &!
! #op&      +(sst_SB(ikp)                                              &!~Prescribed   SST
! #op&       -TsisSV(ikp,1,isl))*(1-LSMask(ikp,1)    )  * SSTnud       &! (Nudging)
! #OP&                                               )  * VarSST       &! Interactive  SST
     &      + 0.


! Sea Ice
! -------

! #AO         eta_SV(ikp,2,isl) =                                      &!
! #AO&        eta_SV(ikp,2,isl) *   LSMask(ikp,2)                       ! Sea: Humidity:=0

! #AO         TsisSV(ikp,2,isl) =                                      &!
! #AO&       (TsisSV(ikp,2,isl) =   LSMask(ikp,2)                      &! Soil Temperature
! #AO&       +  271.2           *(1-LSMask(ikp,2)    ))                 ! Prescribed    ST
            END DO

! #AO       DO isl =         0,nsnow                                    !
! #AO         TsisSV(ikp,2,isl) =                                      &!
! #AO&        s_T_AO(ikp,2)     *(1-LSMask(ikp,2)    )                  ! Prescribed   SIT
! #AO       END DO
      ENDDO


      DO      ikp = 1,kcolp
      DO      nm  = 1,mwp  
          Ts__SV(ikp,nm)         =                                     &! Surf.Temperature
     &    TsisSV(ikp,nm,0)                                             &!
     &                    *(1-min(1,isnoSV    (ikp,nm)  ))             &!
     &   +TsisSV(ikp,nm,      max(1,isnoSV    (ikp,nm)  ))             &!
     &                    *   min(1,isnoSV    (ikp,nm)  )              &!
     &       +0.
      ENDDO
      ENDDO



! Mosaic Cleaning
! ===============

      DO      ikp = 1,kcolp
        IF (LSMask(ikp,1)    .EQ.   0)                              THEN
            isnoSV    (ikp,1)     = 0
            ispiSV    (ikp,1)     = 0
            iiceSV    (ikp,1)     = 0
          DO ISL=1,nsnow
            TsisSV    (ikp,1,isl) = 0.
            dzsnSV    (ikp,1,isl) = 0.
            ro__SV    (ikp,1,isl) = 0.
            eta_SV    (ikp,1,isl) = 0.
            G1snSV    (ikp,1,isl) = 0.
            G2snSV    (ikp,1,isl) = 0.
            agsnSV    (ikp,1,isl) = 0.
            istoSV    (ikp,1,isl) = 0.
          ENDDO
        ENDIF
        IF (FracSV(ikp,n2) .LT. epsn .AND. n2 .GT. 1)               THEN
            Ts__SV(ikp,n2)         = Ts__SV(ikp,1)  
          DO isl=-nsoil,0 
            TsisSV    (ikp,n2,isl) = TsisSV(ikp,1,isl)
          ENDDO    ! #n2  
            isnoSV    (ikp,n2    ) = isnoSV(ikp,1    ) * LSMask(ikp,n2)
            ispiSV    (ikp,n2    ) = ispiSV(ikp,1    ) * LSMask(ikp,n2)
            iiceSV    (ikp,n2    ) = iiceSV(ikp,1    ) * LSMask(ikp,n2)
          DO isl=1,nsnow      !   
            TsisSV    (ikp,n2,isl) = TsisSV(ikp,1,isl) * LSMask(ikp,n2)
            dzsnSV    (ikp,n2,isl) = dzsnSV(ikp,1,isl) * LSMask(ikp,n2)
            ro__SV    (ikp,n2,isl) = ro__SV(ikp,1,isl) * LSMask(ikp,n2)
            eta_SV    (ikp,n2,isl) = eta_SV(ikp,1,isl) * LSMask(ikp,n2)
            G1snSV    (ikp,n2,isl) = G1snSV(ikp,1,isl) * LSMask(ikp,n2)
            G2snSV    (ikp,n2,isl) = G2snSV(ikp,1,isl) * LSMask(ikp,n2)
            agsnSV    (ikp,n2,isl) = agsnSV(ikp,1,isl) * LSMask(ikp,n2)
            istoSV    (ikp,n2,isl) = istoSV(ikp,1,isl) * LSMask(ikp,n2)
          ENDDO    ! #n2  
        ENDIF      ! #n2  
      ENDDO



! Grid Averages / Diagnostics
! ===========================

! Grid Averages                                                   (OUTPUT)
! ^^^^^^^^^^^^^                                                    ^^^^^^
      DO ikp = 1,kcolp
      DO nm  = 1,mwp
! #IB         wes_SV    (ikp,nm)= - wes0SV(ikp,nm) + wes_SV(ikp,nm)     ! Depo. / Subli.
! #IB         wem_SV    (ikp,nm)=   wem0SV(ikp,nm) + wem_SV(ikp,nm)     ! Melting
! #IB         wer_SV    (ikp,nm)=   wer0SV(ikp,nm) + wer_SV(ikp,nm)     ! Refreezing
! #IB         wee_SV    (ikp,nm)= - uqs_SV(ikp,nm) + wee_SV(ikp,nm)    &! Evapotranspiration
! #IB&                            * dt__SV         * roa_SV(ikp,nm,1)   !

              Alb_SV_gpt(ikp)   =   Alb_SV_gpt(ikp)                    &! Grid   Albedo
     &      + FracSV    (ikp,nm)  * albcSV    (ikp,nm)                  ! Mosaic Albedo
              EmisSV_gpt(ikp)   =   EmisSV_gpt(ikp)                    &! Grid   Emissivity
     &      + FracSV    (ikp,nm)  * emi_SV    (ikp,nm)                  ! Mosaic Emissivity
              LWUsRT    (ikp)   =   LWUsRT    (ikp)                    &!
     &      + FracSV    (ikp,nm)  * IRu_SV    (ikp,nm)                  ! Mosaic Upw.IR
              LMO_SV_gpt(ikp)   =   LMO_SV_gpt(ikp)                    &!
     &      + FracSV    (ikp,nm)  * LMO_SV    (ikp,nm)                  ! Mosaic Mon.Ob.
              us__SV_gpt(ikp)   =   us__SV_gpt(ikp)                    &! Grid   u*
     &      + FracSV    (ikp,nm)  * us__SV    (ikp,nm)                  ! Mosaic u*
              uts_SV_gpt(ikp)   =   uts_SV_gpt(ikp)                    &! Grid   u*T*
     &      + FracSV    (ikp,nm)  * uts_SV    (ikp,nm)                  ! Mosaic u*T*
              uqs_SV_gpt(ikp)   =   uqs_SV_gpt(ikp)                    &! Grid   u*q*
     &      + FracSV    (ikp,nm)  * uqs_SV    (ikp,nm)                  ! Mosaic u*q*
! #BS         uss_CM    (ikp)   =   uss_CM    (ikp)                    &! Grid   u*s*
! #BS&      + FracSV    (ikp,nm)  * ussxSV    (ikp,nm)                  ! Mosaic u*s*
!     NO !    ussxSV    (ikp,nm)  = uss_SV    (ikp,nm)                  !        u*s*
!        !    Upper Update = wrong Source of Atmospher.Snow             !
! #BS         qs__CM    (ikp,mzpp) =                                   &! Salt.Part.Concent.
! #BS&        qs__CM    (ikp,mzpp)                                     &!
! #BS&      + FracSV    (ikp,nm)  * qSalSV    (ikp,nm)                 &! 
! #BS&                       *min(1,isnoSV    (ikp,nm)  )               !
! #PO         hFraSV_gpt(ikp)     = hFraSV_gpt(ikp)                    &! Frazil  Thickness
! #PO&      + FracSV    (ikp,nm)  * HFraSV    (ikp,nm)                  ! 
              Tas_SV_gpt(ikp)     = Tas_SV_gpt(ikp)                    &! Surface Air
     &      + FracSV    (ikp,nm)  * Ts__SV    (ikp,nm)                  !         Temperatur
              qv__DY    (ikp,mzpp)= qv__DY    (ikp,mzpp)               &!
     &      + FracSV    (ikp,nm)  * SHumSV    (ikp,nm)                  ! Surf.Specif.Humid.
                                                                        ! to adapt over soil
! #TC         uqTC      (i,j,1)   = uqTC      (i,j,1)                  &! Grid   u*d*
! #TC&      + FracSV    (ikp,nm)  * uds_SV    (ikp,nm)                  ! Mosaic u*d*
! #TC         qsTC      (i,j,1)   = qsTC      (i,j,1)                  &! Salt.Part.Concent.
! #TC&      + FracSV    (ikp,nm)  * qSalSV    (ikp,nm)                 &! 
! #TC&                    *(1-min(1,isnoSV    (ikp,nm)  ))              !
      ENDDO
      ENDDO



! De la grille "SISVAT" vers la grille "AtmPHY"
! =============================================

      DO ikp = 1,kcolp
        Sno0CM    (ikp) = SnowCM    (ikp)                               !
         rhAir          = roa_SV    (ikp,1,1)                           ! Air    Densitity
        Ta__DY    (ikp,mzpp) =                                         &!
     &                    Tas_SV_gpt(ikp)                               !
        pkt_DY    (ikp,mzpp) =                                         &!
     &                    Tas_SV_gpt(ikp) /pkPaSV(ikp,1) ** RCp         !
        HsenSV_gpt(ikp) =-uts_SV_gpt(ikp) *rhAir    *CpdAir             ! Sensible Heat Flux
        HLatSV_gpt(ikp) =-uqs_SV_gpt(ikp) *rhAir    *LhvH2O             ! Latent   Heat Flux
        WE2aSV_gpt(ikp) = WE2aSV_gpt(ikp)                              &! Total    Evaporat.
     &                   -uqs_SV_gpt(ikp) *rhAir    *dt__SV             ! [mm w.e.]
      END DO



! Blown Snow/Dust Accumulation
! ============================

! #AE DO ikp=1,kcolp
! #AE DO nm =1,mwp

! #BS    SnowCM(ikp) = SnowCM(ikp)                                     &!
! #BS&               + dt__SV *1.e-3 *roa_SV(ikp,1) *uss_CM(ikp)

! #AE END DO
! #AE END DO



! Formation of Lakes                                                                 \VER
! ==================



! Sea-Ice Ice Floe Size
! =====================

! Prescription from SST
! ---------------------

! #SI IF (VarSST.le.eps6)                                           THEN
! #SI   DO ikp=1,kcolp
! #SI        FraOcn          =(TsisSV    (ikp,1,0)-TOF_SV)/TSIdSV       ! Prescribed
!                                                                       ! from SST
! #IP        FraOcn          =  1.-sif_SB(ikp)                          ! Prescribed
!                                                                       ! from SSM/I
! #SI        FraOcn          = min(  un_1,FraOcn)                       ! UpperLimit
! #SI        FraOcn          = max(OcnMin,FraOcn)                       ! LowerLimit
! #SI        FracSV(ikp,1)   =     LSMask(ikp,1)     * FracSV(ikp,1)   &! New Ocean
! #SI&                        + (1-LSMask(ikp,1)    )* FraOcn           !
! #SI        SrfSIC          =                         FracSV(ikp,2)    ! Old Sea Ice
! #SI        SIc0OK          = max(zer0, sign(un_1,SrfSIC-eps6))        !
! #SI        FracSV(ikp,2)     *   LSMask(ikp,2)     * FracSV(ikp,2)   &! New Sea Ice
! #SI&                        + (1-LSMask(ikp,2)    )*(1.-FraOcn)       !
! #SI        SIceOK          = max(zer0, sign(un_1,FracSV(ikp,2)       &!
! #SI&                                                   -eps6))        !

! Sea-Ice Vertical Discretization
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #SI        isnoSV(ikp,2)     =                                       &
! #SI&       isnoSV(ikp,2)                     *   LSMask(ikp,2)       &
! #SI&     +(max(1                                                     &
! #SI&      ,isnoSV(ikp,2)  )  *SIc0OK                                 &
! #SI&      +     3        *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2))

! #SI        iiceSV(ikp,2)     =                                       &
! #SI&       iiceSV(ikp,2)                     *   LSMask(ikp,2)       &
! #SI&     +(max(1                                                     &
! #SI&      ,iiceSV(ikp,2)  )  *SIc0OK                                 &
! #SI&      +     3        *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2))
! #SI        ispiSV(ikp,2)     =                   iiceSV(ikp,2)  

! #SI     DO l=1,nsnow
! #SI        dzsnSV(ikp,2,l) =                                         &
! #SI&       dzsnSV(ikp,2,l)                   *   LSMask(ikp,2)       &
! #SI&     +(max                                                       &
! #SI&      (SIc_OK(min(2,l))*    SIcMIN                               &
! #SI&      ,dzsnSV(ikp,2,l))*SIc0OK                                   &
! #SI&      +dzSIce(min(4,l))*(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2))

! #SI        TsisSV(ikp,2,l) =                                         &
! #SI&       TsisSV(ikp,2,l)                 *   LSMask(ikp,2)         &
! #SI&     +(TsisSV(ikp,2,l) *    SIc0OK                               &
! #SI&      +TsisSV(ikp,1,0) *(1.-SIc0OK)   )*(1-LSMask(ikp,2))

! #SI        ro__SV(ikp,2,l) =                                         &
! #SI&       ro__SV(ikp,2,l)                 *   LSMask(ikp,2)         &
! #SI&     +(max                                                       &
! #SI&      (SIc_OK(min(2,l))*rhoIce                                   &
! #SI&      ,ro__SV(ikp,2,l))*SIc0OK                                   &
! #SI&      +rhoIce      *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2))

! #SI        G1snSV(ikp,2,l) =                                         &
! #SI&       G1snSV(ikp,2,l)                 *   LSMask(ikp,2)         &
! #SI&     +(G1snSV(ikp,2,l) *SIc0OK                                   &
! #SI&      +G1_dSV      *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2))

! #SI        G2snSV(ikp,2,l) =                                         &
! #SI&       G2snSV(ikp,2,l)                 *   LSMask(ikp,2)         &
! #SI&     +(G2snSV(ikp,2,l) *SIc0OK                                   &
! #SI&      +30.         *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2))

! #SI        istoSV(ikp,2,l) =                                         &
! #SI&       istoSV(ikp,2,l)                 *   LSMask(ikp,2)         &
!!#SI&new  +(istoSV(ikp,2,l) *SIc0OK                                   &
!!#SI&new   +istdSV(2)   *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2))        &
! #SI&     + istdSV(2)   *                    (1-LSMask(ikp,2))
! #SI     END DO
! #SI     DO l=-nsoil,0
! #SI        TsisSV(ikp,2,l) =                                         &
! #SI&       TsisSV(ikp,2,l)                 *   LSMask(ikp,2)         &
! #SI&     +(TsisSV(ikp,2,l) *    SIc0OK                               &
! #SI&      +TsisSV(ikp,1,l) *(1.-SIc0OK)   )*(1-LSMask(ikp,2))

! #SI        eta_SV(ikp,2,l) =                                         &
! #SI&       eta_SV(ikp,2,l)                 *   LSMask(ikp,2)         &
! #SI&     + eta_SV(ikp,2,l) *    SIc0OK     *(1-LSMask(ikp,2))
!                                 No Pore in Ice => No Water
! #SI     END DO

! OUTPUT of Sea Ice Fraction, Temperature, Discretization
! #si        write(6,6001) Day_TU,LabMon(Mon_TU),YearTU                &
! #si&                    ,HourTU,minuTU,Sec_TU ,TsisSV(ikp,1,0)       &
! #si&                    ,FraOcn,FracSV(ikp,1) ,TsisSV(ikp,2,0)       &
! #si&                    ,       iiceSV(ikp,2) ,isnoSV(ikp,2)  

! #SI   END DO
! #SI END IF

! Otherwise SST and FrLead have been computed in the Sea-Ice Polynya Model
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! Rainfall, Snowfall Time Integral at previous Time Step
! ------------------------------------------------------

        DO ikp = 1,kcolp
           Rai0CM(ikp) = RainCM(ikp)                                    ! Rainfall Time Integral
           SnobCM(ikp) = SnowCM(ikp)                                    ! Snowfall Time Integral
                                                                        ! Erosion  skipped


!  Wind Horizontal Components         at previous Time Step
!  --------------------------------------------------------

        DO nm=1,mwp
! #Za      Ua0_SV(ikp,nm) = Ua__SV(ikp,nm)                              !
! #Za      Va0_SV(ikp,nm) = Va__SV(ikp,nm)                              !
        END DO
        END DO


! Additional OUTPUT for VERIFICATION
! ----------------------------------

! OUTPUT/Verification: SBL Turbulent Fluxes
! #sh i   = i_x0 + 10.*111.111e3/dxHOST
! #sh j   = j_y0
! #sh nm  = 1
! #sh ikp = ikl_AP(i,j)
! #sh write(6,6060) it_EXP,Day_TU,LabMon(Mon_TU),YearTU                    &
! #sh&                    ,HourTU,minuTU        ,lat__r(ikp)/Dg2Rad        &
! #sh&         ,TaT_SV    (ikp,nm)              ,roa_SV(ikp,nm,1)          &
! #sh&         ,HsenSV_gpt(ikp),HLatSV_gpt(ikp),-86400.*uqs_SV_gpt(ikp)    &
! #sh&        ,1.e3*RainCM(ikp),WE2aSV_gpt(ikp),        ROF_SV    (ikp,nm)  
! #sh 6060 format(i6,i3,'-',a3,'-',i4,':',i2,':',i2,f6.2,' N',             &
! #sh&            f9.3,' K'     ,f6.3,' kg/m3',2(f6.1,' W/m2'),            &
! #sh&            f6.3,' mm/day',3(f9.3,' mm'))



! SISVAT OUTPUTs
! ==============

!     IF (WRIsbc)                                                   THEN

!            **************
!       CALL PHY_SISVAT_IOs
!            **************

!     END IF



      return
      end subroutine PHY_SISVAT_RUN
