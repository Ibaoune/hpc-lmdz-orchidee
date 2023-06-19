SUBROUTINE SU_AERP

!**** *SU_AERP*   - INITIALIZE MODULES YOEAERSRC, YOEAERSNK

!     PURPOSE.
!     --------
!           INITIALIZE YOEAERSRC AND YOEAERSNK, THE MODULES THAT CONTAINS 
!           COEFFICIENTS NEEDED TO RUN THE PROGNOSTIC AEROSOLS

!**   INTERFACE.
!     ----------
!        *CALL* *SU_AERP

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        YOEAERSRC, YOEAERSNK, YOEAERATM

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE *ECMWF*
!        from O.BOUCHER (LOA, 1998-03) 

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2004-05-10

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOEAERSRC , ONLY : RSSFLX

USE YOEAERSNK , ONLY : R_R, R_S, RALPHAR, RALPHAS, RFRAER, RFRGAS, &
  &  RRHMAX, RRHTAB, RRHO_SS, RSSGROW, RMMD_SS, RMMD_DD, RRHO_DD, &
  &  RFRBC , RFRIF , RFROM  , RFRSO4 , RFRDD  , RFRSS , RHO_WAT, RHO_ICE, &
  &  RVDPOCE, RVDPSIC, RVDPLND, RVDPLIC, RVSEDOCE, RVSEDSIC, RVSEDLND, RVSEDLIC, &
  &  NBRH

USE YOEAERATM , ONLY : RMASSE, RMFMIN, NDD1, NSS1

IMPLICIT NONE

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!     ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SU_AERP',0,ZHOOK_HANDLE)

!-- For the ECMWF model, the following tables when dimensioned to 12 
!   can refer to 12 values of RH
!      (RHTAB, RHHO_SS, RSSGROW)
!   or to 12 types/bins of aerosols with the following mapping: NTYPAER
!   1- 3  sea-salt  0.03 - 0.5 -  5  - 20 microns		   1
!   4- 6  dust      0.03 - 0.5 - 0.9 - 20 microns		   2
!   7- 8  POM	    hydrophilic, hydrophobic			   3
!   9-10  BC	    hydrophilic, hydrophobic			   4
!  11     sulfate						   5
!  12     fly ash						   6
!  13     stratospheric aerosols				   7
!  14     volcanic aerosols					   8
!  15    							   9
!      (RVDPOCE, RVDSIC, RVDPLND, RVDPLIC)
!      (RVSEDOCE,RVSEDSIC,RVSEDLND,RVSEDLIC)



!*      1.    PARAMETERS RELATED TO SOURCES
!             -----------------------------

!-- OB's original 12 types and 24 different values had the following mapping
! DMS SO2 SO4 H2S DMSO MSA H2O2
! BC(2) OM(2) FlyAsh DU(2) SS(10)

!-- OB's original SS 10 bins     
! bin sizes: 0.03-0.06-0.13-0.25-0.5-1.0-2.0-5.0-10.-15.-20
!RSSFLX = (/ &
!  &  0.20526E-09_JPRB, 0.49292E-09_JPRB, 0.97079E-09_JPRB, 0.31938E-08_JPRB &
!  &, 0.16245E-07_JPRB, 0.86292E-07_JPRB, 0.31326E-06_JPRB, 0.24671E-06_JPRB &
!  &, 0.14109E-06_JPRB, 0.11784E-06_JPRB /)

! maximum possible number of aerosol types
!NMAXTAER=9   already defined in SU_AERW

!N.B. Fluxes of sea salt for each size bin are given in mg m-2 s-1 at wind 
!     speed of 1 m s-1 at 10m height (at 80% RH) in OB's seasalt.F
! RSSFLX also in mg m-2 s-1       
!-- OB's ECMWF 3 bins of sea salt: 0.03, 0.5, 5, 20 microns  
RSSFLX = (/ 4.85963536E-09_JPRB, 4.15358556E-07_JPRB, 5.04905813E-07_JPRB /)

! OB's original vdep_oce, vdep_sic, vdep_ter, vdep_lic were given over 24 values 

! following 12 values for 10 SS and 2 DU in cm s-1
!RVDPOCE = (/   0.1_JPRB,   1.2_JPRB,   0.1_JPRB,   0.1_JPRB,   0.1_JPRB,   0.1_JPRB &
!          &,   0.1_JPRB,   1.2_JPRB,   1.2_JPRB,   1.2_JPRB,   1.5_JPRB,   1.5_JPRB /)
!
!RVDPSIC = (/   0.1_JPRB,   1.2_JPRB,   0.1_JPRB,   0.1_JPRB,   0.1_JPRB,   0.1_JPRB &
!          &,   0.1_JPRB,   1.2_JPRB,   1.2_JPRB,   1.2_JPRB,   1.5_JPRB,   1.5_JPRB /)
!
!RVDPLND = (/   0.1_JPRB,   1.2_JPRB,   0.1_JPRB,   0.1_JPRB,   0.1_JPRB,   0.1_JPRB &
!          &,   0.1_JPRB,   1.2_JPRB,   1.2_JPRB,   1.2_JPRB,   1.5_JPRB,   1.5_JPRB /)
!
!RVDPLIC = (/   0.1_JPRB,   1.2_JPRB,   0.1_JPRB,   0.1_JPRB,   0.1_JPRB,   0.1_JPRB &
!          &,   0.1_JPRB,   1.2_JPRB,   1.2_JPRB,   1.2_JPRB,   1.5_JPRB,   1.5_JPRB /)



!*      2.    PARAMETERS RELATED TO SINKS
!             ---------------------------

R_R = 0.001_JPRB
R_S = 0.001_JPRB

RFRAER = 0.5_JPRB
RFRGAS = 1.0_JPRB

!*      2.1   SEA SALT
!             -------- 
!-- parameters related to SEA SALT: 12 relates to 12 values of relative humidity

NBRH=12
RRHMAX = 95._JPRB
RRHTAB = (/ 0._JPRB, 10._JPRB, 20._JPRB, 30._JPRB, 40._JPRB, 50._JPRB &
       & , 60._JPRB, 70._JPRB, 80._JPRB, 85._JPRB, 90._JPRB, 95._JPRB /)
RRHO_SS = (/ 2160._JPRB, 2160._JPRB, 2160._JPRB, 2160._JPRB, 1451.6_JPRB &
     & , 1367.9_JPRB, 1302.9_JPRB, 1243.2_JPRB, 1182.7_JPRB, 1149.5_JPRB &
     & , 1111.6_JPRB, 1063.1_JPRB /)
RSSGROW = (/ 0.503_JPRB, 0.503_JPRB, 0.503_JPRB, 0.503_JPRB, 0.724_JPRB &
     & , 0.782_JPRB, 0.838_JPRB, 0.905_JPRB, 1.000_JPRB, 1.072_JPRB &
     & , 1.188_JPRB, 1.447_JPRB /)

!-- OB's original 10 bins !RMMD_SS = (/ 0.09_JPRB, 0.19_JPRB, 0.38_JPRB, 0.75_JPRB, 1.50_JPRB &
!         & , 3.00_JPRB, 7.00_JPRB, 15.0_JPRB, 25.0_JPRB, 35.0_JPRB /)

!-- OB's ECMWF 3 bins of sea salt
!  bins are 0.03 - 0.5 -  5.0  - 20 microns

RMMD_SS = (/ 0.30_JPRB, 3.00_JPRB, 10.00_JPRB /)
RFRSS   = (/  0.7_JPRB,  0.7_JPRB,   0.7_JPRB /)
RHO_WAT = 1000._JPRB
RHO_ICE = 500._JPRB

!-  computed off-line by gems_ss.f  (m s-1)

RVSEDOCE(1:3) = (/   0.24E-04_JPRB,   0.20E-02_JPRB,   0.20E-01_JPRB /)
RVSEDSIC(1:3) = (/   0.24E-04_JPRB,   0.20E-02_JPRB,   0.20E-01_JPRB /)
RVSEDLND(1:3) = (/   0.24E-04_JPRB,   0.20E-02_JPRB,   0.20E-01_JPRB /)
RVSEDLIC(1:3) = (/   0.24E-04_JPRB,   0.20E-02_JPRB,   0.20E-01_JPRB /)

! adapted from LMDZ   (m s-1)

RVDPOCE(1:3) = (/ 0.100E-02_JPRB, 0.110E-01_JPRB,   0.145E-01_JPRB /)
RVDPSIC(1:3) = (/ 0.100E-02_JPRB, 0.110E-01_JPRB,   0.145E-01_JPRB /)
RVDPLND(1:3) = (/ 0.100E-02_JPRB, 0.110E-01_JPRB,   0.145E-01_JPRB /)
RVDPLIC(1:3) = (/ 0.100E-02_JPRB, 0.110E-01_JPRB,   0.145E-01_JPRB /)


!*      2.2   DESERT DUST
!             ----------- 
!- parameters related to DESERT DUST  (OB's ECMWF 3 bins)
!  bins are 0.03 - 0.55 -  0.9  - 20 microns

RMMD_DD = (/  0.32_JPRB,  0.75_JPRB,   9.0_JPRB /)
RRHO_DD = (/ 2600._JPRB, 2600._JPRB, 2600._JPRB /)
RFRDD   = (/   0.7_JPRB,   0.7_JPRB,   0.7_JPRB /)

!-  computed off-line by gems_dust.f

RVSEDOCE(4:6) = (/ 0.70E-04_JPRB, 0.20E-03_JPRB, 0.20E-02_JPRB /)
RVSEDSIC(4:6) = (/ 0.70E-04_JPRB, 0.20E-03_JPRB, 0.20E-02_JPRB /)
RVSEDLND(4:6) = (/ 0.70E-04_JPRB, 0.20E-03_JPRB, 0.20E-02_JPRB /)
RVSEDLIC(4:6) = (/ 0.70E-04_JPRB, 0.20E-03_JPRB, 0.20E-02_JPRB /)

! adapted from LMDZ   (m s-1)

RVDPOCE(4:6)  = (/ 0.100E-02_JPRB, 0.500E-02_JPRB, 0.110E-01_JPRB /)
RVDPSIC(4:6)  = (/ 0.100E-02_JPRB, 0.500E-02_JPRB, 0.110E-01_JPRB /)
RVDPLND(4:6)  = (/ 0.100E-02_JPRB, 0.500E-02_JPRB, 0.110E-01_JPRB /)
RVDPLIC(4:6)  = (/ 0.100E-02_JPRB, 0.500E-02_JPRB, 0.110E-01_JPRB /)


!*      2.3   OTHER AEROSOLS (to be improved later!)
!             --------------
!- parameters related to other aerosol types
!- particulate organic matter POM
RFROM = (/  0.0_JPRB,  0.7_JPRB /)

RVSEDOCE(7:8) = (/   0.10E+00_JPRB,   0.10E+00_JPRB /)
RVSEDSIC(7:8) = (/   0.10E+00_JPRB,   0.10E+00_JPRB /)
RVSEDLND(7:8) = (/   0.10E+00_JPRB,   0.10E+00_JPRB /)
RVSEDLIC(7:8) = (/   0.10E+00_JPRB,   0.10E+00_JPRB /)

! adapted from LMDZ   (m s-1)

RVDPOCE(7:8) = (/   0.10E-02_JPRB,   0.10E-02_JPRB /)
RVDPSIC(7:8) = (/   0.10E-02_JPRB,   0.10E-02_JPRB /)
RVDPLND(7:8) = (/   0.10E-02_JPRB,   0.10E-02_JPRB /)
RVDPLIC(7:8) = (/   0.10E-02_JPRB,   0.10E-02_JPRB /)


!- black carbon
RFRBC = (/  0.0_JPRB,  0.7_JPRB /)

RVSEDOCE(9:10) = (/   0.10E+00_JPRB,   0.10E+00_JPRB /)
RVSEDSIC(9:10) = (/   0.10E+00_JPRB,   0.10E+00_JPRB /)
RVSEDLND(9:10) = (/   0.10E+00_JPRB,   0.10E+00_JPRB /)
RVSEDLIC(9:10) = (/   0.10E+00_JPRB,   0.10E+00_JPRB /)

! adapted from LMDZ   (m s-1)

RVDPOCE(9:10) = (/   0.10E-02_JPRB,   0.10E-02_JPRB /)
RVDPSIC(9:10) = (/   0.10E-02_JPRB,   0.10E-02_JPRB /)
RVDPLND(9:10) = (/   0.10E-02_JPRB,   0.10E-02_JPRB /)
RVDPLIC(9:10) = (/   0.10E-02_JPRB,   0.10E-02_JPRB /)

!- sulfate
RFRSO4 = 0.7_JPRB

RVSEDOCE(11) = 0.05_JPRB
RVSEDSIC(11) = 0.25_JPRB
RVSEDLND(11) = 0.25_JPRB
RVSEDLIC(11) = 0.25_JPRB

! adapted from LMDZ   (m s-1)

RVDPOCE(11) = 0.05E-02_JPRB
RVDPSIC(11) = 0.25E-02_JPRB
RVDPLND(11) = 0.25E-02_JPRB
RVDPLIC(11) = 0.25E-02_JPRB

!- fly ash
RFRIF  = 0.7_JPRB

RVSEDOCE(12) = 0.20E+00_JPRB
RVSEDSIC(12) = 0.20E+00_JPRB
RVSEDLND(12) = 0.20E+00_JPRB
RVSEDLIC(12) = 0.20E+00_JPRB

! adapted from LMDZ   (m s-1)

RVDPOCE(12) = 0.20E-02_JPRB
RVDPSIC(12) = 0.20E-02_JPRB
RVDPLND(12) = 0.20E-02_JPRB
RVDPLIC(12) = 0.20E-02_JPRB




!- NB: 15 values for all possible types of ECMWF aerosols
RALPHAR = (/ &
 & 0.001_JPRB, 0.001_JPRB, 0.001_JPRB, 0.001_JPRB, 0.001_JPRB, 0.001_JPRB, 0.001_JPRB, 0.001_JPRB, &
 & 0.001_JPRB, 0.001_JPRB, 0.001_JPRB, 0.001_JPRB, 0.001_JPRB, 0.001_JPRB, 0.001_JPRB /)
RALPHAS = (/ &
 &  0.01_JPRB,  0.01_JPRB,  0.01_JPRB,  0.01_JPRB,  0.01_JPRB,  0.01_JPRB,  0.01_JPRB,  0.01_JPRB, &
 &  0.01_JPRB,  0.01_JPRB,  0.01_JPRB,  0.01_JPRB,  0.01_JPRB,  0.01_JPRB,  0.01_JPRB /)

										 
!*      3.    PARAMETERS RELATED TO TRANSPORT WITHIN THE FREE ATMOSPHERE
!             ----------------------------------------------------------

NDD1=4
NSS1=1

RMFMIN = 1.E-10_JPRB

RMASSE = (/ &
  &  6.02E+23_JPRB, 6.02E+23_JPRB, 6.02E+23_JPRB, 6.02E+23_JPRB, 6.02E+23_JPRB, 6.02E+23_JPRB &
  &, 6.02E+23_JPRB, 6.02E+23_JPRB, 6.02E+23_JPRB, 6.02E+23_JPRB, 6.02E+23_JPRB, 6.02E+23_JPRB &
  &, 6.02E+23_JPRB, 6.02E+23_JPRB, 6.02E+23_JPRB /)


!     ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SU_AERP',1,ZHOOK_HANDLE)
END SUBROUTINE SU_AERP

