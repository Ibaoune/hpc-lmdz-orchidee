      subroutine PHY_SISVAT_INI

!------------------------------------------------------------------------------+
!                                                         Sat 29-Jun-2013  MAR |
!     SubRoutine PHY_SISVAT_INI initializes                                    |
!                    SISVAT    (Soil                                           |
!                               Ice                                            |
!                               Snow                                           |
!                               Vegetation                                     |
!                               Atmosphere                                     |
!                               Transfer   Scheme)                             |
!                                                                              |
!     version 3.p.4.1 created by H. Gallee,               Mon  4-Feb-2013      |
!           Last Modification by H. Gallee,               Sat 29-Jun-2013      |
!                                                                              |
!------------------------------------------------------------------------------+


! Global  Variables
! =================

      use Mod_Real
      use Mod_PHY____grd
      use Mod_PHY____dat
      use Mod_PHY____kkl
      use Mod_PHY_DY_kkl
      use Mod_SISVAT_ctr
      use Mod_SISVAT_grd
      use Mod_SISVAT_dat
      use Mod_SISVAT_loc
      use Mod_SISVAT_kkl
      use Mod_SISVAT_gpt



      IMPLICIT NONE


      integer                 ::  isl   ,n     ,k     ,l     ,jo    ,no
      integer                 ::  i     ,j     ,ikp   ,nm



! Normalized Decay of the Surficial Water Content: data
! ===============================================

      real(kind=real8)  ::  c1_zuo = 12.960e+4      !  Run Off Parameters
      real(kind=real8)  ::  c2_zuo =  2.160e+6      !  86400*1.5 day     ...*25 days 
                                                    !  86400*0.3 day    (Modif. ETH Camp)
      real(kind=real8)  ::  c3_zuo =  1.400e+2      !  (Zuo and Oerlemans 1996, J.Glacio. 42, 305--317)

!     real(kind=real8)  ::  c1_zuo =  2.796e+4      !  Run Off Parameters (Tuning)



  
! INITIALISATION:  BEGIN  ++++++++++++++++++++++++++++++++++++++++++++  


! Initialisation of Mod_SISVAT_dat (Ocean Surface Parameters)
! ================================

! #SI VarSST =   0.
! #OP VarSST =   1.                  ! Variable (0.) / Fixed    (1.) SSTs
      FixSST =   1.-VarSST           ! Fixed    (1.) / Variable (0.) SSTs 
      SSTnud = exp(-dt__SV/2.592e6)  ! SST    Nudging:
!                                    ! e-folding time: 30 Days

! #SI TOF_SV = TocnSI                ! Ocn Grid Cell Freez.Temperature  [K]
! #RE TOF_SV = 271.35 + eps6         !


  

! Initialization of Mod_SISVAT_grd 
! ================================

! Parameters used in the Interpolation of V(10 m)
! -----------------------------------------------

        IF   (hsigma(1).GT.10.          )                           THEN
          k = 0
  301     CONTINUE
          k = k + 1
            !if (hsigma(k).LT.10.OR.k.gt.mzp)                   GO TO 300
            if (hsigma(k).LT.10.OR.k.ge.mzp)                   GO TO 300
                                                               GO TO 301
  300     CONTINUE
          k_SL = k

          IF (k_SL.EQ.mzp)                                          THEN
              r_SL10 = log(10.           / 0.002)   &! 0.002: typical Z0
     &                /log(hsigma(k_SL)  / 0.002)    !

          ELSE
              r_SL10 =    (10.           - hsigma(k_SL))               &
     &                /   (hsigma(k_SL-1)- hsigma(k_SL))
          END IF
        ELSE
              k_SL   =     mzp
              r_SL10 =      1.
        END IF



! Level of negligible blown Particles Concentration ( z_zb ~ 25 magl)      SNOW
! -------------------------------------------------                   .OR. DUST

! #AE        k_zb  =mzp
! #AE 11 CONTINUE
! #AE    IF (hsigma(k_zb  ).GT.z_zb.OR.k_zb  .LE.1)             GO TO 10
! #AE        k_zb  =k_zb  -1
! #AE                                                           GO TO 11
! #AE 10 CONTINUE
! #AE        write(6,1000)             k_zb
! #AE 1000   format(/,' BS : Level of negligible '                     &
! #AE&                     ,'blown Particles Concentration is',i4      &
! #AE&                     ,' (i.e., ~  25. magl)',/)




! Initialization of SISVAT constants and parameters
! =================================================

!           **********
      CALL  SISVAT_ini
!           **********




! Initialization of Mod_SISVAT_kkl 
! ================================

! Surface Fall Line Slope
! -----------------------

        IF (kcolp .EQ. 1)                                           THEN
            slopAP    (1)    =      slop1d
          DO n   = 1,mwp
            slopSV(1  ,n)    =      slop1d
          ENDDO
        ELSE
          DO ikp=1,kcolp
          DO nm =1,mwp
            slopSV(ikp,nm)   = slopAP(ikp)
          END DO
          END DO
        END IF

        IF (SnoMod)                                                 THEN
          DO ikp=1,kcolp
          DO nm =1,mwp
            SWf_SV(ikp,nm)      =                            &! Normalized Decay of the
     &       exp(-dt__SV                                  &! Surficial Water Content
     &          /(c1_zuo                                  &! Zuo and Oerlemans 1996,
     &           +c2_zuo*exp(-c3_zuo*slopSV(ikp,nm) )))       ! J.Glacio. 42, 305--317
          END DO
          END DO
        END IF





! OUTPUT point (i,j,n) coordinates
! --------------------------------

! stdout
! ~~~~~~
! Martin CONTROL
          iwr_SV = 23
          jwr_SV = 17
          nwr_SV = 1
          !iwr_SV = 1
          !jwr_SV = 1
          !nwr_SV = 1
! Martin CONTROL

! Particular txt files (for a regular spacing in 1D arrays)
! ~~~~~~~~~~~~~~~~~~~~
        DO  ikp=1,kcolp
        DO  nm =1,mwp
          lwriSV(ikp,nm) = 0
        END DO
        END DO

          jo             = kcolp*mwp/nbwri
          jo             = max(1,jo)
          no             =               0
          ikp            =      -jo /    2
 11     CONTINUE
          nm             =       1
          no             = no  + 1
          ikp            = ikp + jo
          ikp            = max(    1,ikp)
          ikp            = min(kcolp,ikp)

 101    CONTINUE
          i___SV(no)     = ii__AP(ikp)
          j___SV(no)     = jj__AP(ikp)
          n___SV(no)     = nm
          lwriSV(ikp,nm) = 1
          nm             = nm  + 1
          no             = no  + 1
        IF (nm.LE.mwp .AND. no.LE.nbwri)                       GO TO 101
          no          = no  - 1

        IF (no.LT.nbwri)                                       GO TO 11




! +++++  ++++++++++++++  +++++  ++++++++++++++++++++++++++++++++++++++++
! FIRST  INITIALISATION: BEGIN
! +++++  ++++++++++++++  +++++

!                 ===========
        IF       (it_EXP.EQ.1)                                      THEN
!                 ===========


! Ocean                          1st Initialization (Grid Cells)
! -------------------------------------------------

! #SI     DO ikp=1,kcolp
! #SI       IF  (MaskSV_gpt(ikp) .EQ. 0 .AND. FracSV(ikp,1) .lt. 1.)  THEN
! #SI            i   = ii__AP(ikp)
! #SI            j   = jj__AP(ikp)
! #SI          DO n=1,mwp
! #SI            write(6,6000)i,j,n,FracSV(ikp,nm)
! #SI 6000       format(' WARNING on Grid Point',2i4,' Mosaic',i3      &
! #SI&                 ,' Fraction = f4.3,' : ISLAND excluded')
! #SI            FracSV    (ikp,nm)  =   0.
! #SI            iVgTSV    (ikp,nm)  =   0
! #SI          END DO
! #SI            FracSV    (ikp,1)   =   1.
! #SI       END IF
! #SI     END DO


! Prescription from SST
! ---------------------

! #SI     DO ikp=1,kcolp
! #SI        FraOcn          =(TsisSV(ikp,1,0)-TOF_SV)/TSIdSV           ! Open Ocean
! #IP        FraOcn          = 1.0000         -sif_SB(ikp)              ! Prescribed
! #SI        FraOcn          = min(  un_1,FraOcn)                       !      Fract.
! #SI        FraOcn          = max(OcnMin,FraOcn)                       !
! #SI        FracSV(ikp,1)   = LSMask(ikp,1)     *    FracSV(ikp,1)    &! New  Ocean
! #SI&                     +(1-LSMask(ikp,1)    )*    FraOcn            !
! #SI        SrfSIC          =                        FracSV(ikp,2)     ! Old  Sea Ice
! #SI        SIc0OK          = max(zer0, sign(un_1,   SrfSIC-eps6))     !
! #SI        FracSV(ikp,2)   = LSMask(ikp,1)     *    FracSV(ikp,2)    &! New  Sea Ice
! #SI&                     +(1-LSMask(ikp,1)    )*(1.-FraOcn)           !
! #SI        SIceOK          = max(zer0, sign(un_1,   FracSV(ikp,2)    &!
! #SI&                                                      -eps6))     !

! Sea-Ice Vertical Discretization
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #SI        isnoSV    (ikp,2)     =                                   &
! #SI&       isnoSV    (ikp,2)                   *   LSMask(ikp,2)     &
! #SI&     +(isnoSV    (ikp,2)  * SIc0OK                               &
! #SI&      +     3          *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2) ) 

! #SI        iiceSV    (ikp,2)     =                                   &
! #SI&       iiceSV    (ikp,2)                   *   LSMask(ikp,2)     &
! #SI&     +(iiceSV    (ikp,2)     *SIc0OK                             &
! #SI&      +     3          *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2) )
! #SI        ispiSV    (ikp,2)     =                 iiceSV(ikp,2)  

! #SI       DO l=1,nsnow
! #SI        dzsnSV    (ikp,2,l) =                                     &
! #SI&       dzsnSV    (ikp,2,l)                 *   LSMask(ikp,2)     &
! #SI&     +(dzsnSV    (ikp,2,l) *SIc0OK                               &
! #SI&      +dzSIce(min(4,l))*(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2) )

! #SI        TsisSV    (ikp,2,l) =                                     &
! #SI&       TsisSV    (ikp,2,l)                 *   LSMask(ikp,2)     &
! #SI&     +(TsisSV    (ikp,2,l) *    SIc0OK                           &
! #SI&      +TsisSV    (ikp,1,0) *(1.-SIc0OK)   )*(1-LSMask(ikp,2) )

! #SI        ro__SV    (ikp,2,l) =                                     &
! #SI&       ro__SV    (ikp,2,l)                 *   LSMask(ikp,2)     &
! #SI&     +(ro__SV    (ikp,2,l) *SIc0OK                               &
! #SI&      +rhoIce          *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2) )

! #SI        G1snSV    (ikp,2,l) =                                     &
! #SI&       G1snSV    (ikp,2,l)                 *   LSMask(ikp,2)     &
! #SI&     +(G1snSV    (ikp,2,l) *SIc0OK                               &
! #SI&      +G1_dSV          *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2) )

! #SI        G2snSV    (ikp,2,l) =                                     &
! #SI&       G2snSV    (ikp,2,l)                 *   LSMask(ikp,2)     &
! #SI&     +(G2snSV    (ikp,2,l) *SIc0OK                               &
! #SI&      +30.             *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2) )

! #SI        istoSV    (ikp,2,l) =                                     &
! #SI&       istoSV    (ikp,2,l)                 *   LSMask(ikp,2)     &
!!#SI&new  +(istoSV    (ikp,2,l) *SIc0OK                               &
!!#SI&new   +istdSV(2)       *(1.-SIc0OK)*SIceOK)*(1-LSMask(ikp,2) )   &
! #SI&     + istdSV(2)           *                (1-LSMask(ikp,2) )
! #SI       END DO
! #SI       DO l=-nsoil,0
! #SI        TsisSV    (ikp,2,l) =                                     &
! #SI&       TsisSV    (ikp,2,l)                 *   LSMask(ikp,2)     &
! #SI&     +(TsisSV    (ikp,2,l) *    SIc0OK                           &
! #SI&      +TsisSV    (ikp,1,l) *(1.-SIc0OK)   )*(1-LSMask(ikp,2) )

! #SI        eta_SV    (ikp,2,l) =                                     &
! #SI&       eta_SV    (ikp,2,l)                 *   LSMask(ikp,2)     &
! #SI&     + eta_SV    (ikp,2,l) *    SIc0OK     *(1-LSMask(ikp,2) )
!                                 No Pore in Ice => No Water
! #SI       END DO

! OUTPUT of Sea Ice Fraction, Temperature, Discretization
! #si        write(6,6001) Day_TU,LabMon(Mon_TU),YearTU                &
! #si&                    ,HourTU,minuTU,Sec_TU ,TsisSV(ikp,1,0)       &
! #si&                    ,FraOcn,FracSV(ikp,1) ,TsisSV(ikp,2,0)       &
! #si&                    ,iiceSV(ikp,2)        ,isnoSV(ikp,2)  
! #si 6001        format(/,98('_')                                     &
! #si&         ,/, i3,'-',a3,'-',i4,3(':',i2)                          &
! #si&         ,   2x,'T OCN = ',f7.3,4x,'% OCN = ',f7.3,'(',f4.3,')'  &
! #si&         ,   2x,'T ICE = ',f7.3                                  &
! #si&         ,/,43x,'NbIce = ',i3, 11x,'NbSno = ',i3)      

! #SI     END DO


! SBL                            1st Initialization
! -------------------------------------------------

! Influence of the Angle(Wind,Sastrugi) (Andreas, 1995, CCREL report 95-16)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #Za     DO ikp=1,kcolp
! #Za     DO nm =1,mwp
! #Za       ua0_SV(ikp,nm)  = ua__SV(ikp,nm)
! #Za       va0_SV(ikp,nm)  = va__SV(ikp,nm)
! #Za     END DO
! #Za     END DO

! H2O  Upward IR Flux                Initialization
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          DO ikp=1,kcolp
            WE2aSV_gpt(ikp) = 0.0
          END DO


! GUESS, eventually from DATA   (1st Initialization)
! --------------------------------------------------

          DO ikp = 1,kcolp
          DO nm  = 1,mwp

            IRs_SV(ikp,nm)    =-StefBo*Ta__DY(ikp,mzpp) *Ta__DY(ikp,mzpp)     &! Upward IR Flux
     &                                *Ta__DY(ikp,mzpp) *Ta__DY(ikp,mzpp)      !


! SBL                            1st Initialization (Mosaics)
! -------------------------------------------------

! Drag Coefficient               1st Initialization
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rCDmSV(ikp,nm)    =    0.04
            rCDhSV(ikp,nm)    =    0.04

! Turbulent Scales               1st Initialization
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            us__SV(ikp,nm)    = rCDmSV(ikp,nm) * WindSV(ikp,nm,mzp)
            uts_SV(ikp,nm)    = rCDhSV(ikp,nm) *(Ta__DY(ikp,mzp)-Ta__DY(ikp,mzpp) ) &
     &                        * us__SV(ikp,nm)
            uqs_SV(ikp,nm)    =    0.0
            uqs_SV_gpt(ikp)   =    0.0                                   !  redondance sur n, sans consequences
            uss_SV(ikp,nm)    =    0.0
! #BS       ussxSV(ikp,nm)    =    0.0

! Orography Roughness Lengths
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #ZO       Z0roSV(ikp,nm)    = min(z0__SV(ikp,nm),hsigma(mzp)/3.)

! Roughness Lengths              1st Initialization
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #ZS       Z0SaSV(ikp,nm)    =    0.0                                   !  z0(Sastrugi  h)
! #ZM      DO k=1,ntavz
! #ZM       z0_mem(ikp,nm,k)  =    0.0001
! #ZM       r0_mem(ikp,nm,k)  =    0.0001
! #ZM       b0_mem(ikp,nm,k)  =    0.0001
! #ZM      END DO

! SBL Wind Speed and Vertical Temperature Gradient
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #AA      DO k=1,ntave
! #AW       V__mem(ikp,nm,k)  =    WindSV(ikp,nm,mzp)
! #AH       T__mem(ikp,nm,k)  =    0.0000
! #AA      END DO

! SBL Water Vapor Flux Limitor
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #WL      DO n =1,nLimi
! #WL        WL_mem(ikp,nm,n)  =   1.
! #WL      END DO


! Frazil                         1st Initialization
! -------------------------------------------------

             HFraSV(ikp,nm)    =   0.


! Roughness Length               1st Initialization
! -------------------------------------------------

             Z0roSV(ikp,nm)    =   z0__SV(ikp,nm)  


! Vegetation                     1st Initialization
! -------------------------------------------------

             TvegSV(ikp,nm)    = Ta__DY(ikp,mzpp)     ! Vegetation skin Temperature
             snCaSV(ikp,nm)    =          0.          ! Canopy intercepted Snow
             rrCaSV(ikp,nm)    =          0.          ! Canopy intercepted Raiw
             psivSV(ikp,nm)    =          0.          ! Leaf   Water    Potential     [m]


! Blowing Snow                   1st Initialization
! -------------------------------------------------

! #Za        VVs_SV(ikp,nm)    =         10.          ! Wind Speed,         (Sastrugi), Relevance  [m/s]
! #Za        RRsxSV(ikp,nm)    =          1.          ! Wind Speed Counter, (Sastrugi), Relevance    [-]
! #Za        DDsxSV(ikp,nm)    =          0.          ! Wind Direction    , (Sastrugi), Relevance    [-]

             BufsSV(ikp,nm)    =          0.          ! Fallen Snow Buffer          (Mosaic)   [mm w.e.]
             ussbSV(ikp,nm)    =          0.          ! Eroded Snow Buffer          (Mosaic)   [mm w.e.]



! Snow Pack                      1st Initialization
! -------------------------------------------------

             zWE_SV(ikp,nm)    =            0.        ! Thickness 
             qSalSV(ikp,nm)    =            0.        ! Saltating Particles Concentration
             BrosSV(ikp,nm)    =          300.        ! Buffer Snow Layer: initial density: Polar Snow
             BG1sSV(ikp,nm)    =       G1_dSV         ! Buffer Snow Layer: initial G1     : Polar Snow
             BG2sSV(ikp,nm)    =       ADSdSV         ! Buffer Snow Layer: initial G2     : Polar Snow
                                                   ! Buffer Snow Layer  initial characteristics for 0-mass

             wes_SV(ikp,nm)    =            0.        ! Depo. / Subli.
             wem_SV(ikp,nm)    =            0.        ! Melting
             wer_SV(ikp,nm)    =            0.        ! Refreezing
             wee_SV(ikp,nm)    =            0.        ! Evapotranspiration


! Surficial  Water                   Initialization
! -------------------------------------------------

             SWf_SV(ikp,nm)    =            0.        ! Normalized Decay of Surficial Water Content  [-]
             rusnSV(ikp,nm)    =            0.        ! Surficial Water Mass        (Mosaic)   [mm w.e.]
             SWS_SV(ikp,nm)    =            1.        ! Surficial Water Status      (Mosaic)   [0,1=m,f]
                                                   ! Freezng Initial Status   is assumed
 

! Cumulative Run-Off                 Initialization
! -------------------------------------------------

             ROF_SV(ikp,nm)       =  0.


! Soil Roots                     1st Initialization
! -------------------------------------------------

           DO isl = -nsoil,0
             Rootsv(ikp,nm,isl) =    0.
           END DO

          END DO
          END DO

!                 ===========
        END IF ! (it_EXP.EQ.0)
!                 ===========
!
! +++++  ++++++++++++++  +++
! FIRST  INITIALISATION: END
! +++++  ++++++++++++++  +++  ++++++++++++++++++++++++++++++++++++++++++




! OUTPUT Files Set-Up
! =======================

! #v0   open(unit=50,status='unknown',file='PHY_SISVAT.v0')
! #v0   rewind    50




! OUTPUT
! ======

      IF (kcolp.EQ.1)                                               THEN
         write(6,6) (iVgTSV    (1,nm),nm=1,mwp)
 6       format(/ ,' Vegetation Type: ',6i9)
      END IF



! +++ INITIALISATION:  END  ++++++++++++++++++++++++++++++++++++++++++  

      return
      end subroutine PHY_SISVAT_INI
