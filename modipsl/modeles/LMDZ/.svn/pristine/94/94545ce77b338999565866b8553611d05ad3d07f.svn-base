
! $Header$

SUBROUTINE printflag(tabcntr0, radpas, ok_journe, ok_instan, ok_region)



  ! Auteur :  P. Le Van

  IMPLICIT NONE

  REAL tabcntr0(100)
  LOGICAL cycle_diurn0, soil_model0, new_oliq0, ok_orodr0
  LOGICAL ok_orolf0, ok_limitvr0
  LOGICAL ok_journe, ok_instan, ok_region
  INTEGER radpas, radpas0

  include "clesphys.h"


  PRINT 100
  PRINT *, ' ******************************************************* &
    &                                                         &
    &  ************'
  PRINT *, ' ********   Choix  des principales  cles de la physique &
    &                                                         &
    &      *********'
  PRINT *, ' ******************************************************* &
    &                                                         &
    &  ************'
  PRINT 100
  PRINT 10, iflag_cycle_diurne.GE.1, soil_model
  PRINT 100

  IF (iflag_con==1) THEN
    PRINT *, ' *****           Shema  convection   LMD        &
      &                                                       &
      &                   ******'
  ELSE IF (iflag_con==2) THEN
    PRINT *, ' *****           Shema  convection  Tiedtke     &
      &                                                       &
      &                   ******'
  ELSE IF (iflag_con>=3) THEN
    PRINT *, ' *****           Shema  convection    Emanuel   &
      &                                                       &
      &                   ******'
  END IF
  PRINT 100

  PRINT 11, new_oliq, ok_orodr, ok_orolf
  PRINT 100

  PRINT 7, ok_limitvrai
  PRINT 100

  PRINT 12, nbapp_rad
  PRINT 100

  PRINT 8, radpas
  PRINT 100

  PRINT 4, ok_journe, ok_instan, ok_region
  PRINT 100
  PRINT 100


  cycle_diurn0 = .FALSE.
  soil_model0 = .FALSE.
  new_oliq0 = .FALSE.
  ok_orodr0 = .FALSE.
  ok_orolf0 = .FALSE.
  ok_limitvr0 = .FALSE.

  IF (tabcntr0(7)==1.) cycle_diurn0 = .TRUE.
  IF (tabcntr0(8)==1.) soil_model0 = .TRUE.
  IF (tabcntr0(9)==1.) new_oliq0 = .TRUE.
  IF (tabcntr0(10)==1.) ok_orodr0 = .TRUE.
  IF (tabcntr0(11)==1.) ok_orolf0 = .TRUE.
  IF (tabcntr0(12)==1.) ok_limitvr0 = .TRUE.

  PRINT *, ' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ &
    &                                                         &
    & $$$$$$$$$$$$$'
  PRINT 100

  IF (int(tabcntr0(5))/=iflag_con) THEN
    PRINT 20, int(tabcntr0(5)), iflag_con
    PRINT 100
  END IF

  IF (int(tabcntr0(6))/=nbapp_rad) THEN
    PRINT 21, int(tabcntr0(6)), nbapp_rad
    ! radpas0  = NINT( 86400./tabcntr0(1)/INT( tabcntr0(6) ) )
    PRINT 100
    ! PRINT 22, radpas0, radpas
    PRINT 100
  END IF

  IF (cycle_diurn0 .AND. .NOT. (iflag_cycle_diurne.GE.1) .OR. .NOT. cycle_diurn0 .AND. &
      (iflag_cycle_diurne.GE.1) ) THEN
    PRINT 13, cycle_diurn0, iflag_cycle_diurne
    PRINT 100
  END IF

  IF (soil_model0 .AND. .NOT. soil_model .OR. .NOT. soil_model0 .AND. &
      soil_model) THEN
    PRINT 14, soil_model0, soil_model
    PRINT 100
  END IF

  IF (new_oliq0 .AND. .NOT. new_oliq .OR. .NOT. new_oliq0 .AND. new_oliq) &
      THEN
    PRINT 16, new_oliq0, new_oliq
    PRINT 100
  END IF

  IF (ok_orodr0 .AND. .NOT. ok_orodr .OR. .NOT. ok_orodr0 .AND. ok_orodr) &
      THEN
    PRINT 15, ok_orodr0, ok_orodr
    PRINT 100
  END IF

  IF (ok_orolf0 .AND. .NOT. ok_orolf .OR. .NOT. ok_orolf0 .AND. ok_orolf) &
      THEN
    PRINT 17, ok_orolf0, ok_orolf
    PRINT 100
  END IF

  IF (ok_limitvr0 .AND. .NOT. ok_limitvrai .OR. .NOT. ok_limitvr0 .AND. &
      ok_limitvrai) THEN
    PRINT 18, ok_limitvr0, ok_limitvrai
    PRINT 100
  END IF

  PRINT 100
  PRINT *, ' ******************************************************* &
    &                                                         &
    &  ************'
  PRINT 100

4 FORMAT (2X, 5('*'), '  ok_journe= ', L3, 3X, ',ok_instan = ', L3, 3X, &
    ',ok_region = ', L3, 3X, 5('*'))

7 FORMAT (2X, 5('*'), 15X, '      ok_limitvrai   = ', L3, 16X, 5('*'))

8 FORMAT (2X, '*****             radpas    =                      ', I4, 6X, &
    ' *****')

10 FORMAT (2X, 5('*'), '    Cycle_diurne = ', L3, 4X, ', Soil_model = ', L3, &
    12X, 6('*'))


11 FORMAT (2X, 5('*'), '  new_oliq = ', L3, 3X, ', Ok_orodr = ', L3, 3X, &
    ', Ok_orolf = ', L3, 3X, 5('*'))


12 FORMAT (2X, '*****  Nb d appels /jour des routines de rayonn. = ', I4, 6X, &
    ' *****')

13 FORMAT (2X, '$$$$$$$$   Attention !!  cycle_diurne  different  sur', /1X, &
    10X, ' startphy = ', L3, 2X, ' et  run.def = ', L3)

14 FORMAT (2X, '$$$$$$$$   Attention !!    soil_model  different  sur', /1X, &
    10X, ' startphy = ', L3, 2X, ' et  run.def = ', L3)

15 FORMAT (2X, '$$$$$$$$   Attention !!      ok_orodr  different  sur', /1X, &
    10X, ' startphy = ', L3, 2X, ' et  run.def = ', L3)

16 FORMAT (2X, '$$$$$$$$   Attention !!      new_oliq  different  sur', /1X, &
    10X, ' startphy = ', L3, 2X, ' et  run.def = ', L3)

17 FORMAT (2X, '$$$$$$$$   Attention !!      ok_orolf  different  sur', /1X, &
    10X, ' startphy = ', L3, 2X, ' et  run.def = ', L3)

18 FORMAT (2X, '$$$$$$$$   Attention !!  ok_limitvrai  different  sur', /1X, &
    10X, ' startphy = ', L3, 2X, ' et  run.def = ', L3)

20 FORMAT (/2X, '$$$$$$$$   Attention !!    iflag_con  different  sur', /1X, &
    10X, ' startphy = ', I3, 2X, ' et  run.def = ', I3)

21 FORMAT (2X, '$$$$$$$$   Attention !!     nbapp_rad  different  sur', /1X, &
    10X, ' startphy = ', I3, 2X, ' et  run.def = ', I3)

22 FORMAT (2X, '$$$$$$$$   Attention !!        radpas  different  sur', /1X, &
    10X, ' startphy = ', I3, 2X, ' et  run.def = ', I3)

100 FORMAT (/)

  RETURN
END SUBROUTINE printflag
