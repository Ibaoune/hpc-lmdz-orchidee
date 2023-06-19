SUBROUTINE SUAER15

!**** *SUAER15*   - INITIALIZE COMMON YOMAER15

!     PURPOSE.
!     --------
!           INITIALIZE YOMAER15, THE COMMON THAT CONTAINS THE
!           RADIATIVE CHARACTERISTICS OF THE AEROSOLS

!**   INTERFACE.
!     ----------
!              -----        -----

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOMAER15

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     04-06: F. Bouyssel. Meteo-France
!     05-09: A. Alias - PAESOD (black_carbon) is added to PAELAN (oraganic)
!                     - the sulfate is put in a separate type 
!                     (see Hu Ron Ming work) - P.Marquet

!=======================================================================
!-- The (old) five aerosol types were respectively:

!  1/ continental average (+desert)       2/ maritime
!  3/ urban                               4/ volcanic active
!  5/ stratospheric background

!=======================================================================

!-- The (new) six aerosol types are respectively:

!  1/ continental average                 2/ maritime
!  3/ desert                              4/ volcanic active
!  5/ stratospheric background            6/ sulfate

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMRAD15 , ONLY : LNEWAER15
USE YOMAER15 , ONLY : RAER15    ,TAUA15   ,RPIZA15  ,RCGA15

IMPLICIT NONE

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SUAER15',0,ZHOOK_HANDLE)
IF (.NOT.LNEWAER15) THEN
  
  RAER15=RESHAPE((/&
   & .038520_JPRB, .037196_JPRB, .040532_JPRB, .054934_JPRB, .038520_JPRB ,&
   & .12613_JPRB , .18313_JPRB , .10357_JPRB , .064106_JPRB, .126130_JPRB ,&
   & .012579_JPRB, .013649_JPRB, .018652_JPRB, .025181_JPRB, .012579_JPRB ,&
   & .011890_JPRB, .016142_JPRB, .021105_JPRB, .028908_JPRB, .011890_JPRB ,&
   & .013792_JPRB, .026810_JPRB, .052203_JPRB, .066338_JPRB, .013792_JPRB ,&
   & .012579_JPRB, .013649_JPRB, .018652_JPRB, .025181_JPRB, .012579_JPRB /)&
   & ,SHAPE=(/5,6/))  

  TAUA15=RESHAPE((/&
   & .730719_JPRB,.730719_JPRB,&
   & .912819_JPRB,.912819_JPRB,&
   & .725059_JPRB,.725059_JPRB,&
   & .745405_JPRB,.745405_JPRB,&
   & .682188_JPRB,.682188_JPRB,&
   & 1.35059_JPRB,.725059_JPRB /)&
   & ,SHAPE=(/2,6/))  

  RPIZA15=RESHAPE((/&
   & .872212_JPRB,.872212_JPRB,&
   & .982545_JPRB,.982545_JPRB,&
   & .623143_JPRB,.623143_JPRB,&
   & .944887_JPRB,.944887_JPRB,&
   & .997975_JPRB,.997975_JPRB,&
   & .999999_JPRB,.999999_JPRB /)&
   & ,SHAPE=(/2,6/))  

  RCGA15=RESHAPE((/&
   & .647596_JPRB,.647596_JPRB,&
   & .739002_JPRB,.739002_JPRB,&
   & .580845_JPRB,.580845_JPRB,&
   & .662657_JPRB,.662657_JPRB,&
   & .624246_JPRB,.624246_JPRB,&
   & .680845_JPRB,.680845_JPRB /)&
   & ,SHAPE=(/2,6/))  

ELSE

  RAER15( :, 1)= (/&
   & .036271_JPRB, .030153_JPRB, .017343_JPRB, .015002_JPRB, .008806_JPRB /)  
  RAER15( :, 2)= (/&
   & .026561_JPRB, .032657_JPRB, .017977_JPRB, .014210_JPRB, .016775_JPRB /)  
  RAER15( :, 3)= (/&
   & .014897_JPRB, .016359_JPRB, .019789_JPRB, .030777_JPRB, .013341_JPRB /)  
  RAER15( :, 4)= (/&
   & .011890_JPRB, .016142_JPRB, .021105_JPRB, .028908_JPRB, .011890_JPRB /)  
  RAER15( :, 5)= (/&
   & .013792_JPRB, .026810_JPRB, .052203_JPRB, .066338_JPRB, .013792_JPRB /)  
  RAER15( :, 6)= (/&
   & .012579_JPRB, .013649_JPRB, .018652_JPRB, .025181_JPRB, .012579_JPRB /)

  TAUA15(1, :)= (/&
   & 1.69446_JPRB, 1.11855_JPRB, 1.09212_JPRB, 1.03858_JPRB, 1.12044_JPRB, 1.35059_JPRB  /)  
  TAUA15(2, :)= (/&
   & 0.40174_JPRB, 0.89383_JPRB, 0.89546_JPRB, 0.51143_JPRB, 0.32646_JPRB, .725059_JPRB /)  
 
  RPIZA15(1, :)= (/&
   & .9148907_JPRB, .9956173_JPRB, .7504584_JPRB, .9401905_JPRB, .9999999_JPRB, .999999_JPRB/)  
  RPIZA15(2, :)= (/&
   & .8814597_JPRB, .9920407_JPRB, .9239428_JPRB, .9515548_JPRB, .9938563_JPRB, .999999_JPRB/)  
 
  RCGA15(1, :)= (/&
   & 0.729019_JPRB, 0.803129_JPRB, 0.784592_JPRB, .7008249_JPRB, .7270548_JPRB, .680845_JPRB/)  
  RCGA15(2, :)= (/&
   & 0.663224_JPRB, 0.793746_JPRB, 0.696315_JPRB, .6608509_JPRB, .6318786_JPRB, .680845_JPRB/)  

ENDIF

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUAER15',1,ZHOOK_HANDLE)
END SUBROUTINE SUAER15
