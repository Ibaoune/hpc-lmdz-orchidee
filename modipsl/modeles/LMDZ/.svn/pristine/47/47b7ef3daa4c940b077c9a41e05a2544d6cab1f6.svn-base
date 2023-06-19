
! $Id $

SUBROUTINE isccp_cloud_types(debug, debugcol, npoints, sunlit, nlev, ncol, &
    seed, pfull, phalf, qv, cc, conv, dtau_s, dtau_c, top_height, overlap, &
    tautab, invtau, skt, emsfc_lw, at, dem_s, dem_c, fq_isccp, totalcldarea, &
    meanptop, meantaucld, boxtau, boxptop)


  ! Copyright Steve Klein and Mark Webb 2002 - all rights reserved.

  ! This code is available without charge with the following conditions:

  ! 1. The code is available for scientific purposes and is not for
  ! commercial use.
  ! 2. Any improvements you make to the code should be made available
  ! to the to the authors for incorporation into a future release.
  ! 3. The code should not be used in any way that brings the authors
  ! or their employers into disrepute.

  IMPLICIT NONE

  ! NOTE:   the maximum number of levels and columns is set by
  ! the following parameter statement

  INTEGER ncolprint

  ! -----
  ! Input
  ! -----

  INTEGER npoints !  number of model points in the horizontal
  ! PARAMETER(npoints=6722)
  INTEGER nlev !  number of model levels in column
  INTEGER ncol !  number of subcolumns

  INTEGER sunlit(npoints) !  1 for day points, 0 for night time

  INTEGER seed(npoints) !  seed value for random number generator
  ! !  ( see Numerical Recipes Chapter 7)
  ! !  It is recommended that the seed is set
  ! !  to a different value for each model
  ! !  gridbox it is called on, as it is
  ! !  possible that the choice of the samec
  ! !  seed value every time may introduce some
  ! !  statistical bias in the results, particularly
  ! !  for low values of NCOL.

  REAL pfull(npoints, nlev) !  pressure of full model levels (Pascals)
  ! !  pfull(npoints,1)    is    top level of model
  ! !  pfull(npoints,nlev) is bottom level of model

  REAL phalf(npoints, nlev+1) !  pressure of half model levels (Pascals)
  ! !  phalf(npoints,1)    is    top       of model
  ! !  phalf(npoints,nlev+1) is the surface pressure

  REAL qv(npoints, nlev) !  water vapor specific humidity (kg vapor/ kg air)
  ! !         on full model levels

  REAL cc(npoints, nlev) !  input cloud cover in each model level (fraction)
  ! !  NOTE:  This is the HORIZONTAL area of each
  ! !         grid box covered by clouds

  REAL conv(npoints, nlev) !  input convective cloud cover in each model level (fraction)
  ! !  NOTE:  This is the HORIZONTAL area of each
  ! !         grid box covered by convective clouds

  REAL dtau_s(npoints, nlev) !  mean 0.67 micron optical depth of stratiform
  ! !  clouds in each model level
  ! !  NOTE:  this the cloud optical depth of only the
  ! !         cloudy part of the grid box, it is not weighted
  ! !         with the 0 cloud optical depth of the clear
  ! !         part of the grid box

  REAL dtau_c(npoints, nlev) !  mean 0.67 micron optical depth of convective
  ! !  clouds in each
  ! !  model level.  Same note applies as in dtau_s.

  INTEGER overlap !  overlap type

  ! 1=max

  ! 2=rand
  ! 3=max/rand

  INTEGER top_height !  1 = adjust top height using both a computed
  ! !  infrared brightness temperature and the visible
  ! !  optical depth to adjust cloud top pressure. Note
  ! !  that this calculation is most appropriate to compare
  ! !  to ISCCP data during sunlit hours.
  ! !  2 = do not adjust top height, that is cloud top
  ! !  pressure is the actual cloud top pressure
  ! !  in the model
  ! !  3 = adjust top height using only the computed
  ! !  infrared brightness temperature. Note that this
  ! !  calculation is most appropriate to compare to ISCCP
  ! !  IR only algortihm (i.e. you can compare to nighttime
  ! !  ISCCP data with this option)

  REAL tautab(0:255) !  ISCCP table for converting count value to
  ! !  optical thickness

  INTEGER invtau(-20:45000) !  ISCCP table for converting optical thickness
  ! !  to count value

  ! The following input variables are used only if top_height = 1 or
  ! top_height = 3

  REAL skt(npoints) !  skin Temperature (K)
  REAL emsfc_lw !  10.5 micron emissivity of surface (fraction)
  REAL at(npoints, nlev) !  temperature in each model level (K)
  REAL dem_s(npoints, nlev) !  10.5 micron longwave emissivity of stratiform
  ! !  clouds in each
  ! !  model level.  Same note applies as in dtau_s.
  REAL dem_c(npoints, nlev) !  10.5 micron longwave emissivity of convective
  ! !  clouds in each
  ! !  model level.  Same note applies as in dtau_s.
  ! IM reg.dyn BEG
  REAL t1, t2
  ! REAL w(npoints)                   !vertical wind at 500 hPa
  ! LOGICAL pct_ocean(npoints)        !TRUE if oceanic point, FALSE otherway
  ! INTEGER iw(npoints) , nw
  ! REAL wmin, pas_w
  ! INTEGER k, l, iwmx
  ! PARAMETER(wmin=-100.,pas_w=10.,iwmx=30)
  ! REAL fq_dynreg(7,7,iwmx)
  ! REAL nfq_dynreg(7,7,iwmx)
  ! LOGICAL pctj(7,7,iwmx)
  ! IM reg.dyn END
  ! ------
  ! Output
  ! ------

  REAL fq_isccp(npoints, 7, 7) !  the fraction of the model grid box covered by
  ! !  each of the 49 ISCCP D level cloud types

  REAL totalcldarea(npoints) !  the fraction of model grid box columns
  ! !  with cloud somewhere in them.  This should
  ! !  equal the sum over all entries of fq_isccp


  ! ! The following three means are averages over the cloudy areas only.  If
  ! no
  ! ! clouds are in grid box all three quantities should equal zero.

  REAL meanptop(npoints) !  mean cloud top pressure (mb) - linear averaging
  ! !  in cloud top pressure.

  REAL meantaucld(npoints) !  mean optical thickness
  ! !  linear averaging in albedo performed.

  REAL boxtau(npoints, ncol) !  optical thickness in each column

  REAL boxptop(npoints, ncol) !  cloud top pressure (mb) in each column



  ! ------
  ! Working variables added when program updated to mimic Mark Webb's PV-Wave
  ! code
  ! ------

  REAL frac_out(npoints, ncol, nlev) ! boxes gridbox divided up into
  ! ! Equivalent of BOX in original version, but
  ! ! indexed by column then row, rather than
  ! ! by row then column

  REAL tca(npoints, 0:nlev) ! total cloud cover in each model level (fraction)
  ! ! with extra layer of zeroes on top
  ! ! in this version this just contains the values input
  ! ! from cc but with an extra level
  REAL cca(npoints, nlev) ! convective cloud cover in each model level (fraction)
  ! ! from conv

  REAL threshold(npoints, ncol) ! pointer to position in gridbox
  REAL maxocc(npoints, ncol) ! Flag for max overlapped conv cld
  REAL maxosc(npoints, ncol) ! Flag for max overlapped strat cld

  REAL boxpos(npoints, ncol) ! ordered pointer to position in gridbox

  REAL threshold_min(npoints, ncol) ! minimum value to define range in with new threshold
  ! ! is chosen

  REAL dem(npoints, ncol), bb(npoints) !  working variables for 10.5 micron longwave
  ! !  emissivity in part of
  ! !  gridbox under consideration

  REAL ran(npoints) ! vector of random numbers
  REAL ptrop(npoints)
  REAL attrop(npoints)
  REAL attropmin(npoints)
  REAL atmax(npoints)
  REAL atmin(npoints)
  REAL btcmin(npoints)
  REAL transmax(npoints)

  INTEGER i, j, ilev, ibox, itrop(npoints)
  INTEGER ipres(npoints)
  INTEGER itau(npoints), ilev2
  INTEGER acc(nlev, ncol)
  INTEGER match(npoints, nlev-1)
  INTEGER nmatch(npoints)
  INTEGER levmatch(npoints, ncol)

  ! !variables needed for water vapor continuum absorption
  REAL fluxtop_clrsky(npoints), trans_layers_above_clrsky(npoints)
  REAL taumin(npoints)
  REAL dem_wv(npoints, nlev), wtmair, wtmh20, navo, grav, pstd, t0
  REAL press(npoints), dpress(npoints), atmden(npoints)
  REAL rvh20(npoints), wk(npoints), rhoave(npoints)
  REAL rh20s(npoints), rfrgn(npoints)
  REAL tmpexp(npoints), tauwv(npoints)

  CHARACTER *1 cchar(6), cchar_realtops(6)
  INTEGER icycle
  REAL tau(npoints, ncol)
  LOGICAL box_cloudy(npoints, ncol)
  REAL tb(npoints, ncol)
  REAL ptop(npoints, ncol)
  REAL emcld(npoints, ncol)
  REAL fluxtop(npoints, ncol)
  REAL trans_layers_above(npoints, ncol)
  REAL isccp_taumin, fluxtopinit(npoints), tauir(npoints)
  REAL meanalbedocld(npoints)
  REAL albedocld(npoints, ncol)
  REAL boxarea
  INTEGER debug ! set to non-zero value to print out inputs
  ! ! with step debug
  INTEGER debugcol ! set to non-zero value to print out column
  ! ! decomposition with step debugcol

  INTEGER index1(npoints), num1, jj
  REAL rec2p13, tauchk

  CHARACTER *10 ftn09

  DATA isccp_taumin/0.3/
  DATA cchar/' ', '-', '1', '+', 'I', '+'/
  DATA cchar_realtops/' ', ' ', '1', '1', 'I', 'I'/

  tauchk = -1.*log(0.9999999)
  rec2p13 = 1./2.13

  ncolprint = 0

  ! IM
  ! PRINT*,' isccp_cloud_types npoints=',npoints

  ! if ( debug.ne.0 ) then
  ! j=1
  ! write(6,'(a10)') 'j='
  ! write(6,'(8I10)') j
  ! write(6,'(a10)') 'debug='
  ! write(6,'(8I10)') debug
  ! write(6,'(a10)') 'debugcol='
  ! write(6,'(8I10)') debugcol
  ! write(6,'(a10)') 'npoints='
  ! write(6,'(8I10)') npoints
  ! write(6,'(a10)') 'nlev='
  ! write(6,'(8I10)') nlev
  ! write(6,'(a10)') 'ncol='
  ! write(6,'(8I10)') ncol
  ! write(6,'(a10)') 'top_height='
  ! write(6,'(8I10)') top_height
  ! write(6,'(a10)') 'overlap='
  ! write(6,'(8I10)') overlap
  ! write(6,'(a10)') 'emsfc_lw='
  ! write(6,'(8f10.2)') emsfc_lw
  ! write(6,'(a10)') 'tautab='
  ! write(6,'(8f10.2)') tautab
  ! write(6,'(a10)') 'invtau(1:100)='
  ! write(6,'(8i10)') (invtau(i),i=1,100)
  ! do j=1,npoints,debug
  ! write(6,'(a10)') 'j='
  ! write(6,'(8I10)') j
  ! write(6,'(a10)') 'sunlit='
  ! write(6,'(8I10)') sunlit(j)
  ! write(6,'(a10)') 'seed='
  ! write(6,'(8I10)') seed(j)
  ! write(6,'(a10)') 'pfull='
  ! write(6,'(8f10.2)') (pfull(j,i),i=1,nlev)
  ! write(6,'(a10)') 'phalf='
  ! write(6,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
  ! write(6,'(a10)') 'qv='
  ! write(6,'(8f10.3)') (qv(j,i),i=1,nlev)
  ! write(6,'(a10)') 'cc='
  ! write(6,'(8f10.3)') (cc(j,i),i=1,nlev)
  ! write(6,'(a10)') 'conv='
  ! write(6,'(8f10.2)') (conv(j,i),i=1,nlev)
  ! write(6,'(a10)') 'dtau_s='
  ! write(6,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
  ! write(6,'(a10)') 'dtau_c='
  ! write(6,'(8f10.2)') (dtau_c(j,i),i=1,nlev)
  ! write(6,'(a10)') 'skt='
  ! write(6,'(8f10.2)') skt(j)
  ! write(6,'(a10)') 'at='
  ! write(6,'(8f10.2)') (at(j,i),i=1,nlev)
  ! write(6,'(a10)') 'dem_s='
  ! write(6,'(8f10.3)') (dem_s(j,i),i=1,nlev)
  ! write(6,'(a10)') 'dem_c='
  ! write(6,'(8f10.2)') (dem_c(j,i),i=1,nlev)
  ! enddo
  ! endif

  ! ---------------------------------------------------!

  ! assign 2d tca array using 1d input array cc

  DO j = 1, npoints
    tca(j, 0) = 0
  END DO

  DO ilev = 1, nlev
    DO j = 1, npoints
      tca(j, ilev) = cc(j, ilev)
    END DO
  END DO

  ! assign 2d cca array using 1d input array conv

  DO ilev = 1, nlev
    ! IM pas besoin        do ibox=1,ncol
    DO j = 1, npoints
      cca(j, ilev) = conv(j, ilev)
    END DO
    ! IM        enddo
  END DO

  ! IM
  ! do j=1, iwmx
  ! do l=1, 7
  ! do k=1, 7
  ! fq_dynreg(k,l,j) =0.
  ! nfq_dynreg(k,l,j) =0.
  ! enddo !k
  ! enddo !l
  ! enddo !j
  ! IM
  ! IM
  ! if (ncolprint.ne.0) then
  ! do j=1,npoints,1000
  ! write(6,'(a10)') 'j='
  ! write(6,'(8I10)') j
  ! write (6,'(a)') 'seed:'
  ! write (6,'(I3.2)') seed(j)

  ! write (6,'(a)') 'tca_pp_rev:'
  ! write (6,'(8f5.2)')
  ! &   ((tca(j,ilev)),
  ! &      ilev=1,nlev)

  ! write (6,'(a)') 'cca_pp_rev:'
  ! write (6,'(8f5.2)')
  ! &   ((cca(j,ilev),ibox=1,ncolprint),ilev=1,nlev)
  ! enddo
  ! endif

  IF (top_height==1 .OR. top_height==3) THEN

    DO j = 1, npoints
      ptrop(j) = 5000.
      atmin(j) = 400.
      attropmin(j) = 400.
      atmax(j) = 0.
      attrop(j) = 120.
      itrop(j) = 1
    END DO

    DO ilev = 1, nlev
      DO j = 1, npoints
        IF (pfull(j,ilev)<40000. .AND. pfull(j,ilev)>5000. .AND. &
            at(j,ilev)<attropmin(j)) THEN
          ptrop(j) = pfull(j, ilev)
          attropmin(j) = at(j, ilev)
          attrop(j) = attropmin(j)
          itrop(j) = ilev
        END IF
        IF (at(j,ilev)>atmax(j)) atmax(j) = at(j, ilev)
        IF (at(j,ilev)<atmin(j)) atmin(j) = at(j, ilev)
      END DO
    END DO

  END IF

  ! -----------------------------------------------------!

  ! ---------------------------------------------------!

  ! IM
  ! do 13 ilev=1,nlev
  ! num1=0
  ! do j=1,npoints
  ! if (cc(j,ilev) .lt. 0. .or. cc(j,ilev) .gt. 1.) then
  ! num1=num1+1
  ! index1(num1)=j
  ! end if
  ! enddo
  ! do jj=1,num1
  ! j=index1(jj)
  ! write(6,*)  ' error = cloud fraction less than zero'
  ! write(6,*) ' or '
  ! write(6,*)  ' error = cloud fraction greater than 1'
  ! write(6,*) 'value at point ',j,' is ',cc(j,ilev)
  ! write(6,*) 'level ',ilev
  ! STOP
  ! enddo
  ! num1=0
  ! do j=1,npoints
  ! if (conv(j,ilev) .lt. 0. .or. conv(j,ilev) .gt. 1.) then
  ! num1=num1+1
  ! index1(num1)=j
  ! end if
  ! enddo
  ! do jj=1,num1
  ! j=index1(jj)
  ! write(6,*)
  ! &           ' error = convective cloud fraction less than zero'
  ! write(6,*) ' or '
  ! write(6,*)
  ! &           ' error = convective cloud fraction greater than 1'
  ! write(6,*) 'value at point ',j,' is ',conv(j,ilev)
  ! write(6,*) 'level ',ilev
  ! STOP
  ! enddo

  ! num1=0
  ! do j=1,npoints
  ! if (dtau_s(j,ilev) .lt. 0.) then
  ! num1=num1+1
  ! index1(num1)=j
  ! end if
  ! enddo
  ! do jj=1,num1
  ! j=index1(jj)
  ! write(6,*)
  ! &           ' error = stratiform cloud opt. depth less than zero'
  ! write(6,*) 'value at point ',j,' is ',dtau_s(j,ilev)
  ! write(6,*) 'level ',ilev
  ! STOP
  ! enddo
  ! num1=0
  ! do j=1,npoints
  ! if (dtau_c(j,ilev) .lt. 0.) then
  ! num1=num1+1
  ! index1(num1)=j
  ! end if
  ! enddo
  ! do jj=1,num1
  ! j=index1(jj)
  ! write(6,*)
  ! &           ' error = convective cloud opt. depth less than zero'
  ! write(6,*) 'value at point ',j,' is ',dtau_c(j,ilev)
  ! write(6,*) 'level ',ilev
  ! STOP
  ! enddo

  ! num1=0
  ! do j=1,npoints
  ! if (dem_s(j,ilev) .lt. 0. .or. dem_s(j,ilev) .gt. 1.) then
  ! num1=num1+1
  ! index1(num1)=j
  ! end if
  ! enddo
  ! do jj=1,num1
  ! j=index1(jj)
  ! write(6,*)
  ! &           ' error = stratiform cloud emissivity less than zero'
  ! write(6,*)'or'
  ! write(6,*)
  ! &           ' error = stratiform cloud emissivity greater than 1'
  ! write(6,*) 'value at point ',j,' is ',dem_s(j,ilev)
  ! write(6,*) 'level ',ilev
  ! STOP
  ! enddo

  ! num1=0
  ! do j=1,npoints
  ! if (dem_c(j,ilev) .lt. 0. .or. dem_c(j,ilev) .gt. 1.) then
  ! num1=num1+1
  ! index1(num1)=j
  ! end if
  ! enddo
  ! do jj=1,num1
  ! j=index1(jj)
  ! write(6,*)
  ! &           ' error = convective cloud emissivity less than zero'
  ! write(6,*)'or'
  ! write(6,*)
  ! &           ' error = convective cloud emissivity greater than 1'
  ! write (6,*)
  ! &          'j=',j,'ilev=',ilev,'dem_c(j,ilev) =',dem_c(j,ilev)
  ! STOP
  ! enddo
  ! 13    continue


  DO ibox = 1, ncol
    DO j = 1, npoints
      boxpos(j, ibox) = (ibox-.5)/ncol
    END DO
  END DO

  ! ---------------------------------------------------!
  ! Initialise working variables
  ! ---------------------------------------------------!

  ! Initialised frac_out to zero

  DO ilev = 1, nlev
    DO ibox = 1, ncol
      DO j = 1, npoints
        frac_out(j, ibox, ilev) = 0.0
      END DO
    END DO
  END DO

  ! IM
  ! if (ncolprint.ne.0) then
  ! write (6,'(a)') 'frac_out_pp_rev:'
  ! do j=1,npoints,1000
  ! write(6,'(a10)') 'j='
  ! write(6,'(8I10)') j
  ! write (6,'(8f5.2)')
  ! &     ((frac_out(j,ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)

  ! enddo
  ! write (6,'(a)') 'ncol:'
  ! write (6,'(I3)') ncol
  ! endif
  ! if (ncolprint.ne.0) then
  ! write (6,'(a)') 'last_frac_pp:'
  ! do j=1,npoints,1000
  ! write(6,'(a10)') 'j='
  ! write(6,'(8I10)') j
  ! write (6,'(8f5.2)') (tca(j,0))
  ! enddo
  ! endif

  ! ---------------------------------------------------!
  ! ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
  ! frac_out is the array that contains the information
  ! where 0 is no cloud, 1 is a stratiform cloud and 2 is a
  ! convective cloud

    !loop over vertical levels
  DO ilev = 1, nlev

    ! Initialise threshold

    IF (ilev==1) THEN
        ! If max overlap
      IF (overlap==1) THEN
          ! select pixels spread evenly
          ! across the gridbox
        DO ibox = 1, ncol
          DO j = 1, npoints
            threshold(j, ibox) = boxpos(j, ibox)
          END DO
        END DO
      ELSE
        DO ibox = 1, ncol
          CALL ran0_vec(npoints, seed, ran)
            ! select random pixels from the non-convective
            ! part the gridbox ( some will be converted into
            ! convective pixels below )
          DO j = 1, npoints
            threshold(j, ibox) = cca(j, ilev) + (1-cca(j,ilev))*ran(j)
          END DO
        END DO
      END IF
      ! IM
      ! IF (ncolprint.ne.0) then
      ! write (6,'(a)') 'threshold_nsf2:'
      ! do j=1,npoints,1000
      ! write(6,'(a10)') 'j='
      ! write(6,'(8I10)') j
      ! write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
      ! enddo
      ! ENDIF
    END IF

    ! IF (ncolprint.ne.0) then
    ! write (6,'(a)') 'ilev:'
    ! write (6,'(I2)') ilev
    ! ENDIF

    DO ibox = 1, ncol

        ! All versions
      DO j = 1, npoints
        IF (boxpos(j,ibox)<=cca(j,ilev)) THEN
          ! IM REAL           maxocc(j,ibox) = 1
          maxocc(j, ibox) = 1.0
        ELSE
          ! IM REAL           maxocc(j,ibox) = 0
          maxocc(j, ibox) = 0.0
        END IF
      END DO

        ! Max overlap
      IF (overlap==1) THEN
        DO j = 1, npoints
          threshold_min(j, ibox) = cca(j, ilev)
          ! IM REAL           maxosc(j,ibox)=1
          maxosc(j, ibox) = 1.0
        END DO
      END IF

        ! Random overlap
      IF (overlap==2) THEN
        DO j = 1, npoints
          threshold_min(j, ibox) = cca(j, ilev)
          ! IM REAL           maxosc(j,ibox)=0
          maxosc(j, ibox) = 0.0
        END DO
      END IF

        ! Max/Random overlap
      IF (overlap==3) THEN
        DO j = 1, npoints
          threshold_min(j, ibox) = max(cca(j,ilev), min(tca(j,ilev-1),tca(j, &
            ilev)))
          IF (threshold(j,ibox)<min(tca(j,ilev-1),tca(j, &
              ilev)) .AND. (threshold(j,ibox)>cca(j,ilev))) THEN
            ! IM REAL                maxosc(j,ibox)= 1
            maxosc(j, ibox) = 1.0
          ELSE
            ! IM REAL                 maxosc(j,ibox)= 0
            maxosc(j, ibox) = 0.0
          END IF
        END DO
      END IF

        ! Reset threshold
      CALL ran0_vec(npoints, seed, ran)

      DO j = 1, npoints
        threshold(j, ibox) = &       !if max overlapped conv cloud
          maxocc(j, ibox)*(boxpos(j,ibox)) + &   !else
          (1-maxocc(j,ibox))*( &     !if max overlapped strat cloud
          (maxosc(j,ibox))*( &       !threshold=boxpos
          threshold(j,ibox))+ &      !else
          (1-maxosc(j,ibox))*( &     !threshold_min=random[thrmin,1]
          threshold_min(j,ibox)+(1-threshold_min(j,ibox))*ran(j)))
      END DO

    END DO ! ibox

    ! Fill frac_out with 1's where tca is greater than the threshold

    DO ibox = 1, ncol
      DO j = 1, npoints
        IF (tca(j,ilev)>threshold(j,ibox)) THEN
          ! IM REAL             frac_out(j,ibox,ilev)=1
          frac_out(j, ibox, ilev) = 1.0
        ELSE
          ! IM REAL             frac_out(j,ibox,ilev)=0
          frac_out(j, ibox, ilev) = 0.0
        END IF
      END DO
    END DO

    ! Code to partition boxes into startiform and convective parts
    ! goes here

    DO ibox = 1, ncol
      DO j = 1, npoints
        IF (threshold(j,ibox)<=cca(j,ilev)) THEN
            ! = 2 IF threshold le cca(j)
          ! IM REAL                  frac_out(j,ibox,ilev) = 2
          frac_out(j, ibox, ilev) = 2.0
        ELSE
            ! = the same IF NOT threshold le cca(j)
          frac_out(j, ibox, ilev) = frac_out(j, ibox, ilev)
        END IF
      END DO
    END DO

    ! Set last_frac to tca at this level, so as to be tca
    ! from last level next time round

    ! IM
    ! if (ncolprint.ne.0) then

    ! do j=1,npoints ,1000
    ! write(6,'(a10)') 'j='
    ! write(6,'(8I10)') j
    ! write (6,'(a)') 'last_frac:'
    ! write (6,'(8f5.2)') (tca(j,ilev-1))

    ! write (6,'(a)') 'cca:'
    ! write (6,'(8f5.2)') (cca(j,ilev),ibox=1,ncolprint)

    ! write (6,'(a)') 'max_overlap_cc:'
    ! write (6,'(8f5.2)') (maxocc(j,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'max_overlap_sc:'
    ! write (6,'(8f5.2)') (maxosc(j,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'threshold_min_nsf2:'
    ! write (6,'(8f5.2)') (threshold_min(j,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'threshold_nsf2:'
    ! write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'frac_out_pp_rev:'
    ! write (6,'(8f5.2)')
    ! &       ((frac_out(j,ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
    ! enddo
    ! endif


  END DO

  ! ---------------------------------------------------!



  ! ---------------------------------------------------!
  ! COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
  ! put into vector tau

    !initialize tau and albedocld to zero
  ! loop over nlev
  DO ibox = 1, ncol
    DO j = 1, npoints
      tau(j, ibox) = 0.
      albedocld(j, ibox) = 0.
      boxtau(j, ibox) = 0.
      boxptop(j, ibox) = 0.
      box_cloudy(j, ibox) = .FALSE.
    END DO
  END DO

    !compute total cloud optical depth for each column
  DO ilev = 1, nlev
      !increment tau for each of the boxes
    DO ibox = 1, ncol
      DO j = 1, npoints
        ! IM REAL              if (frac_out(j,ibox,ilev).eq.1) then
        IF (frac_out(j,ibox,ilev)==1.0) THEN
          tau(j, ibox) = tau(j, ibox) + dtau_s(j, ilev)
        END IF
        ! IM REAL              if (frac_out(j,ibox,ilev).eq.2) then
        IF (frac_out(j,ibox,ilev)==2.0) THEN
          tau(j, ibox) = tau(j, ibox) + dtau_c(j, ilev)
        END IF
      END DO
    END DO ! ibox
  END DO ! ilev
  ! IM
  ! if (ncolprint.ne.0) then

  ! do j=1,npoints ,1000
  ! write(6,'(a10)') 'j='
  ! write(6,'(8I10)') j
  ! write(6,'(i2,1X,8(f7.2,1X))')
  ! &          ilev,
  ! &          (tau(j,ibox),ibox=1,ncolprint)
  ! enddo
  ! endif

  ! ---------------------------------------------------!




  ! ---------------------------------------------------!
  ! COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
  ! AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE

  ! again this is only done if top_height = 1 or 3

  ! fluxtop is the 10.5 micron radiance at the top of the
  ! atmosphere
  ! trans_layers_above is the total transmissivity in the layers
  ! above the current layer
  ! fluxtop_clrsky(j) and trans_layers_above_clrsky(j) are the clear
  ! sky versions of these quantities.

  IF (top_height==1 .OR. top_height==3) THEN


      !----------------------------------------------------------------------
      !
      !             DO CLEAR SKY RADIANCE CALCULATION FIRST
      !
      !compute water vapor continuum emissivity
      !this treatment follows Schwarkzopf and Ramasamy
      !JGR 1999,vol 104, pages 9467-9499.
      !the emissivity is calculated at a wavenumber of 955 cm-1,
      !or 10.47 microns
    wtmair = 28.9644
    wtmh20 = 18.01534
    navo = 6.023E+23
    grav = 9.806650E+02
    pstd = 1.013250E+06
    t0 = 296.
    ! IM
    ! if (ncolprint .ne. 0)
    ! &         write(6,*)  'ilev   pw (kg/m2)   tauwv(j)      dem_wv'
    DO ilev = 1, nlev
      DO j = 1, npoints
          !press and dpress are dyne/cm2 = Pascals *10
        press(j) = pfull(j, ilev)*10.
        dpress(j) = (phalf(j,ilev+1)-phalf(j,ilev))*10
          !atmden = g/cm2 = kg/m2 / 10
        atmden(j) = dpress(j)/grav
        rvh20(j) = qv(j, ilev)*wtmair/wtmh20
        wk(j) = rvh20(j)*navo*atmden(j)/wtmair
        rhoave(j) = (press(j)/pstd)*(t0/at(j,ilev))
        rh20s(j) = rvh20(j)*rhoave(j)
        rfrgn(j) = rhoave(j) - rh20s(j)
        tmpexp(j) = exp(-0.02*(at(j,ilev)-t0))
        tauwv(j) = wk(j)*1.E-20*((0.0224697*rh20s(j)*tmpexp(j))+(3.41817E-7* &
          rfrgn(j)))*0.98
        dem_wv(j, ilev) = 1. - exp(-1.*tauwv(j))
      END DO
      ! IM
      ! if (ncolprint .ne. 0) then
      ! do j=1,npoints ,1000
      ! write(6,'(a10)') 'j='
      ! write(6,'(8I10)') j
      ! write(6,'(i2,1X,3(f8.3,3X))') ilev,
      ! &           qv(j,ilev)*(phalf(j,ilev+1)-phalf(j,ilev))/(grav/100.),
      ! &           tauwv(j),dem_wv(j,ilev)
      ! enddo
      ! endif
    END DO

      !initialize variables
    DO j = 1, npoints
      fluxtop_clrsky(j) = 0.
      trans_layers_above_clrsky(j) = 1.
    END DO

    DO ilev = 1, nlev
      DO j = 1, npoints

          ! Black body emission at temperature of the layer

        bb(j) = 1/(exp(1307.27/at(j,ilev))-1.)
          !bb(j)= 5.67e-8*at(j,ilev)**4

          ! increase TOA flux by flux emitted from layer
          ! times total transmittance in layers above

        fluxtop_clrsky(j) = fluxtop_clrsky(j) + dem_wv(j, ilev)*bb(j)* &
          trans_layers_above_clrsky(j)

          ! update trans_layers_above with transmissivity
          ! from this layer for next time around loop

        trans_layers_above_clrsky(j) = trans_layers_above_clrsky(j)* &
          (1.-dem_wv(j,ilev))


      END DO
      ! IM
      ! if (ncolprint.ne.0) then
      ! do j=1,npoints ,1000
      ! write(6,'(a10)') 'j='
      ! write(6,'(8I10)') j
      ! write (6,'(a)') 'ilev:'
      ! write (6,'(I2)') ilev

      ! write (6,'(a)')
      ! &        'emiss_layer,100.*bb(j),100.*f,total_trans:'
      ! write (6,'(4(f7.2,1X))') dem_wv(j,ilev),100.*bb(j),
      ! &             100.*fluxtop_clrsky(j),trans_layers_above_clrsky(j)
      ! enddo
      ! endif

    END DO !loop over level

    DO j = 1, npoints
        !add in surface emission
      bb(j) = 1/(exp(1307.27/skt(j))-1.)
        !bb(j)=5.67e-8*skt(j)**4

      fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw*bb(j)* &
        trans_layers_above_clrsky(j)
    END DO

    ! IM
    ! if (ncolprint.ne.0) then
    ! do j=1,npoints ,1000
    ! write(6,'(a10)') 'j='
    ! write(6,'(8I10)') j
    ! write (6,'(a)') 'id:'
    ! write (6,'(a)') 'surface'

    ! write (6,'(a)') 'emsfc,100.*bb(j),100.*f,total_trans:'
    ! write (6,'(4(f7.2,1X))') emsfc_lw,100.*bb(j),
    ! &      100.*fluxtop_clrsky(j),
    ! &       trans_layers_above_clrsky(j)
    ! enddo
    ! endif


      !
      !           END OF CLEAR SKY CALCULATION
      !
      !----------------------------------------------------------------


    ! IM
    ! if (ncolprint.ne.0) then

    ! do j=1,npoints ,1000
    ! write(6,'(a10)') 'j='
    ! write(6,'(8I10)') j
    ! write (6,'(a)') 'ts:'
    ! write (6,'(8f7.2)') (skt(j),ibox=1,ncolprint)

    ! write (6,'(a)') 'ta_rev:'
    ! write (6,'(8f7.2)')
    ! &       ((at(j,ilev2),ibox=1,ncolprint),ilev2=1,nlev)

    ! enddo
    ! endif
      !loop over columns
    DO ibox = 1, ncol
      DO j = 1, npoints
        fluxtop(j, ibox) = 0.
        trans_layers_above(j, ibox) = 1.
      END DO
    END DO

    DO ilev = 1, nlev
      DO j = 1, npoints
          ! Black body emission at temperature of the layer

        bb(j) = 1/(exp(1307.27/at(j,ilev))-1.)
          !bb(j)= 5.67e-8*at(j,ilev)**4
      END DO

      DO ibox = 1, ncol
        DO j = 1, npoints

            ! emissivity for point in this layer
          ! IM REAL             if (frac_out(j,ibox,ilev).eq.1) then
          IF (frac_out(j,ibox,ilev)==1.0) THEN
            dem(j, ibox) = 1. - ((1.-dem_wv(j,ilev))*(1.-dem_s(j,ilev)))
            ! IM REAL             else if (frac_out(j,ibox,ilev).eq.2) then
          ELSE IF (frac_out(j,ibox,ilev)==2.0) THEN
            dem(j, ibox) = 1. - ((1.-dem_wv(j,ilev))*(1.-dem_c(j,ilev)))
          ELSE
            dem(j, ibox) = dem_wv(j, ilev)
          END IF


            ! increase TOA flux by flux emitted from layer
            ! times total transmittance in layers above

          fluxtop(j, ibox) = fluxtop(j, ibox) + dem(j, ibox)*bb(j)* &
            trans_layers_above(j, ibox)

            ! update trans_layers_above with transmissivity
            ! from this layer for next time around loop

          trans_layers_above(j, ibox) = trans_layers_above(j, ibox)* &
            (1.-dem(j,ibox))

        END DO ! j
      END DO ! ibox

      ! IM
      ! if (ncolprint.ne.0) then
      ! do j=1,npoints,1000
      ! write (6,'(a)') 'ilev:'
      ! write (6,'(I2)') ilev

      ! write(6,'(a10)') 'j='
      ! write(6,'(8I10)') j
      ! write (6,'(a)') 'emiss_layer:'
      ! write (6,'(8f7.2)') (dem(j,ibox),ibox=1,ncolprint)

      ! write (6,'(a)') '100.*bb(j):'
      ! write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)

      ! write (6,'(a)') '100.*f:'
      ! write (6,'(8f7.2)')
      ! &         (100.*fluxtop(j,ibox),ibox=1,ncolprint)

      ! write (6,'(a)') 'total_trans:'
      ! write (6,'(8f7.2)')
      ! &          (trans_layers_above(j,ibox),ibox=1,ncolprint)
      ! enddo
      ! endif

    END DO ! ilev


    DO j = 1, npoints
        !add in surface emission
      bb(j) = 1/(exp(1307.27/skt(j))-1.)
        !bb(j)=5.67e-8*skt(j)**4
    END DO

    DO ibox = 1, ncol
      DO j = 1, npoints

          !add in surface emission

        fluxtop(j, ibox) = fluxtop(j, ibox) + emsfc_lw*bb(j)* &
          trans_layers_above(j, ibox)

      END DO
    END DO

    ! IM
    ! if (ncolprint.ne.0) then

    ! do j=1,npoints ,1000
    ! write(6,'(a10)') 'j='
    ! write(6,'(8I10)') j
    ! write (6,'(a)') 'id:'
    ! write (6,'(a)') 'surface'

    ! write (6,'(a)') 'emiss_layer:'
    ! write (6,'(8f7.2)') (dem(1,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') '100.*bb(j):'
    ! write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)

    ! write (6,'(a)') '100.*f:'
    ! write (6,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)
    ! end do
    ! endif

      !now that you have the top of atmosphere radiance account
      !for ISCCP procedures to determine cloud top temperature

      !account for partially transmitting cloud recompute flux
      !ISCCP would see assuming a single layer cloud
      !note choice here of 2.13, as it is primarily ice
      !clouds which have partial emissivity and need the
      !adjustment performed in this section
      !
      !If it turns out that the cloud brightness temperature
      !is greater than 260K, then the liquid cloud conversion
      !factor of 2.56 is used.
      !
      !Note that this is discussed on pages 85-87 of
      !the ISCCP D level documentation (Rossow et al. 1996)

    DO j = 1, npoints
        !compute minimum brightness temperature and optical depth
      btcmin(j) = 1./(exp(1307.27/(attrop(j)-5.))-1.)
    END DO
    DO ibox = 1, ncol
      DO j = 1, npoints
        transmax(j) = (fluxtop(j,ibox)-btcmin(j))/(fluxtop_clrsky(j)-btcmin(j &
          ))
          !note that the initial setting of tauir(j) is needed so that
          !tauir(j) has a realistic value should the next if block be
          !bypassed
        tauir(j) = tau(j, ibox)*rec2p13
        taumin(j) = -1.*log(max(min(transmax(j),0.9999999),0.001))

      END DO

      IF (top_height==1) THEN
        DO j = 1, npoints
          IF (transmax(j)>0.001 .AND. transmax(j)<=0.9999999) THEN
            fluxtopinit(j) = fluxtop(j, ibox)
            tauir(j) = tau(j, ibox)*rec2p13
          END IF
        END DO
        DO icycle = 1, 2
          DO j = 1, npoints
            IF (tau(j,ibox)>(tauchk)) THEN
              IF (transmax(j)>0.001 .AND. transmax(j)<=0.9999999) THEN
                emcld(j, ibox) = 1. - exp(-1.*tauir(j))
                fluxtop(j, ibox) = fluxtopinit(j) - ((1.-emcld(j, &
                  ibox))*fluxtop_clrsky(j))
                fluxtop(j, ibox) = max(1.E-06, (fluxtop(j,ibox)/emcld(j, &
                  ibox)))
                tb(j, ibox) = 1307.27/(log(1.+(1./fluxtop(j,ibox))))
                IF (tb(j,ibox)>260.) THEN
                  tauir(j) = tau(j, ibox)/2.56
                END IF
              END IF
            END IF
          END DO
        END DO

      END IF

      DO j = 1, npoints
        IF (tau(j,ibox)>(tauchk)) THEN
            !cloudy box
          tb(j, ibox) = 1307.27/(log(1.+(1./fluxtop(j,ibox))))
          IF (top_height==1 .AND. tauir(j)<taumin(j)) THEN
            tb(j, ibox) = attrop(j) - 5.
            tau(j, ibox) = 2.13*taumin(j)
          END IF
        ELSE
            !clear sky brightness temperature
          tb(j, ibox) = 1307.27/(log(1.+(1./fluxtop_clrsky(j))))
        END IF
      END DO ! j
    END DO ! ibox

    ! IM
    ! if (ncolprint.ne.0) then

    ! do j=1,npoints,1000
    ! write(6,'(a10)') 'j='
    ! write(6,'(8I10)') j

    ! write (6,'(a)') 'attrop:'
    ! write (6,'(8f7.2)') (attrop(j))

    ! write (6,'(a)') 'btcmin:'
    ! write (6,'(8f7.2)') (btcmin(j))

    ! write (6,'(a)') 'fluxtop_clrsky*100:'
    ! write (6,'(8f7.2)')
    ! &      (100.*fluxtop_clrsky(j))

    ! write (6,'(a)') '100.*f_adj:'
    ! write (6,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'transmax:'
    ! write (6,'(8f7.2)') (transmax(ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'tau:'
    ! write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'emcld:'
    ! write (6,'(8f7.2)') (emcld(j,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'total_trans:'
    ! write (6,'(8f7.2)')
    ! &          (trans_layers_above(j,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'total_emiss:'
    ! write (6,'(8f7.2)')
    ! &          (1.0-trans_layers_above(j,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'total_trans:'
    ! write (6,'(8f7.2)')
    ! &          (trans_layers_above(j,ibox),ibox=1,ncolprint)

    ! write (6,'(a)') 'ppout:'
    ! write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
    ! enddo ! j
    ! endif

  END IF

  ! ---------------------------------------------------!


  ! ---------------------------------------------------!
  ! DETERMINE CLOUD TOP PRESSURE

  ! again the 2 methods differ according to whether
  ! or not you use the physical cloud top pressure (top_height = 2)
  ! or the radiatively determined cloud top pressure (top_height = 1 or 3)


    !compute cloud top pressure
  DO ibox = 1, ncol
      !segregate according to optical thickness
    IF (top_height==1 .OR. top_height==3) THEN
        !find level whose temperature
        !most closely matches brightness temperature
      DO j = 1, npoints
        nmatch(j) = 0
      END DO
      DO ilev = 1, nlev - 1
        ! cdir nodep
        DO j = 1, npoints
          IF ((at(j,ilev)>=tb(j,ibox) .AND. at(j,ilev+1)<tb(j, &
              ibox)) .OR. (at(j,ilev)<=tb(j,ibox) .AND. at(j,ilev+1)>tb(j, &
              ibox))) THEN

            nmatch(j) = nmatch(j) + 1
            IF (abs(at(j,ilev)-tb(j,ibox))<abs(at(j,ilev+1)-tb(j,ibox))) THEN
              match(j, nmatch(j)) = ilev
            ELSE
              match(j, nmatch(j)) = ilev + 1
            END IF
          END IF
        END DO
      END DO

      DO j = 1, npoints
        IF (nmatch(j)>=1) THEN
          ptop(j, ibox) = pfull(j, match(j,nmatch(j)))
          levmatch(j, ibox) = match(j, nmatch(j))
        ELSE
          IF (tb(j,ibox)<atmin(j)) THEN
            ptop(j, ibox) = ptrop(j)
            levmatch(j, ibox) = itrop(j)
          END IF
          IF (tb(j,ibox)>atmax(j)) THEN
            ptop(j, ibox) = pfull(j, nlev)
            levmatch(j, ibox) = nlev
          END IF
        END IF
      END DO ! j

    ELSE ! if (top_height .eq. 1 .or. top_height .eq. 3)

      DO j = 1, npoints
        ptop(j, ibox) = 0.
      END DO
      DO ilev = 1, nlev
        DO j = 1, npoints
          IF ((ptop(j,ibox)==0.) & ! IM  &
                                   ! .and.(frac_out(j,ibox,ilev) .ne. 0))
                                   ! then
              .AND. (frac_out(j,ibox,ilev)/=0.0)) THEN
            ptop(j, ibox) = pfull(j, ilev)
            levmatch(j, ibox) = ilev
          END IF
        END DO
      END DO
    END IF

    DO j = 1, npoints
      IF (tau(j,ibox)<=(tauchk)) THEN
        ptop(j, ibox) = 0.
        levmatch(j, ibox) = 0
      END IF
    END DO

  END DO



  ! ---------------------------------------------------!



  ! ---------------------------------------------------!
  ! DETERMINE ISCCP CLOUD TYPE FREQUENCIES

  ! Now that ptop and tau have been determined,
  ! determine amount of each of the 49 ISCCP cloud
  ! types

  ! Also compute grid box mean cloud top pressure and
  ! optical thickness.  The mean cloud top pressure and
  ! optical thickness are averages over the cloudy
  ! area only. The mean cloud top pressure is a linear
  ! average of the cloud top pressures.  The mean cloud
  ! optical thickness is computed by converting optical
  ! thickness to an albedo, averaging in albedo units,
  ! then converting the average albedo back to a mean
  ! optical thickness.


    !compute isccp frequencies

    !reset frequencies
  DO ilev = 1, 7
    DO ilev2 = 1, 7
      DO j = 1, npoints !
        fq_isccp(j, ilev, ilev2) = 0.
      END DO
    END DO
  END DO

    !reset variables need for averaging cloud properties
  DO j = 1, npoints
    totalcldarea(j) = 0.
    meanalbedocld(j) = 0.
    meanptop(j) = 0.
    meantaucld(j) = 0.
  END DO ! j

  boxarea = 1./real(ncol)

    !determine optical depth category
  ! IM       do 39 j=1,npoints
  ! IM       do ibox=1,ncol
  DO ibox = 1, ncol
    DO j = 1, npoints

      ! IM
      ! CALL CPU_time(t1)
      ! IM

      IF (tau(j,ibox)>(tauchk) .AND. ptop(j,ibox)>0.) THEN
        box_cloudy(j, ibox) = .TRUE.
      END IF

      ! IM
      ! CALL CPU_time(t2)
      ! print*,'IF tau t2 - t1',t2 - t1

      ! CALL CPU_time(t1)
      ! IM

      IF (box_cloudy(j,ibox)) THEN

          ! totalcldarea always diagnosed day or night
        totalcldarea(j) = totalcldarea(j) + boxarea

        IF (sunlit(j)==1) THEN

            ! tau diagnostics only with sunlight

          boxtau(j, ibox) = tau(j, ibox)

            !convert optical thickness to albedo
          albedocld(j, ibox) = real(invtau(min(nint(100.*tau(j,ibox)), &
            45000)))

            !contribute to averaging
          meanalbedocld(j) = meanalbedocld(j) + albedocld(j, ibox)*boxarea

        END IF

      END IF

      ! IM
      ! CALL CPU_time(t2)
      ! print*,'IF box_cloudy t2 - t1',t2 - t1

      ! CALL CPU_time(t1)
      ! IM BEG
      ! IM           !convert ptop to millibars
      ptop(j, ibox) = ptop(j, ibox)/100.

      ! IM           !save for output cloud top pressure and optical
      ! thickness
      boxptop(j, ibox) = ptop(j, ibox)
      ! IM END

      ! IM BEG
        !reset itau(j), ipres(j)
      itau(j) = 0
      ipres(j) = 0

      IF (tau(j,ibox)<isccp_taumin) THEN
        itau(j) = 1
      ELSE IF (tau(j,ibox)>=isccp_taumin .AND. tau(j,ibox)<1.3) THEN
        itau(j) = 2
      ELSE IF (tau(j,ibox)>=1.3 .AND. tau(j,ibox)<3.6) THEN
        itau(j) = 3
      ELSE IF (tau(j,ibox)>=3.6 .AND. tau(j,ibox)<9.4) THEN
        itau(j) = 4
      ELSE IF (tau(j,ibox)>=9.4 .AND. tau(j,ibox)<23.) THEN
        itau(j) = 5
      ELSE IF (tau(j,ibox)>=23. .AND. tau(j,ibox)<60.) THEN
        itau(j) = 6
      ELSE IF (tau(j,ibox)>=60.) THEN
        itau(j) = 7
      END IF

        !determine cloud top pressure category
      IF (ptop(j,ibox)>0. .AND. ptop(j,ibox)<180.) THEN
        ipres(j) = 1
      ELSE IF (ptop(j,ibox)>=180. .AND. ptop(j,ibox)<310.) THEN
        ipres(j) = 2
      ELSE IF (ptop(j,ibox)>=310. .AND. ptop(j,ibox)<440.) THEN
        ipres(j) = 3
      ELSE IF (ptop(j,ibox)>=440. .AND. ptop(j,ibox)<560.) THEN
        ipres(j) = 4
      ELSE IF (ptop(j,ibox)>=560. .AND. ptop(j,ibox)<680.) THEN
        ipres(j) = 5
      ELSE IF (ptop(j,ibox)>=680. .AND. ptop(j,ibox)<800.) THEN
        ipres(j) = 6
      ELSE IF (ptop(j,ibox)>=800.) THEN
        ipres(j) = 7
      END IF
      ! IM END

      IF (sunlit(j)==1 .OR. top_height==3) THEN

        ! IM         !convert ptop to millibars
        ! IM           ptop(j,ibox)=ptop(j,ibox) / 100.

        ! IM         !save for output cloud top pressure and optical
        ! thickness
        ! IM             boxptop(j,ibox) = ptop(j,ibox)

        IF (box_cloudy(j,ibox)) THEN

          meanptop(j) = meanptop(j) + ptop(j, ibox)*boxarea

          ! IM             !reset itau(j), ipres(j)
          ! IM           itau(j) = 0
          ! IM           ipres(j) = 0

          ! if (tau(j,ibox) .lt. isccp_taumin) then
          ! itau(j)=1
          ! else if (tau(j,ibox) .ge. isccp_taumin
          ! &
          ! &          .and. tau(j,ibox) .lt. 1.3) then
          ! itau(j)=2
          ! else if (tau(j,ibox) .ge. 1.3
          ! &          .and. tau(j,ibox) .lt. 3.6) then
          ! itau(j)=3
          ! else if (tau(j,ibox) .ge. 3.6
          ! &          .and. tau(j,ibox) .lt. 9.4) then
          ! itau(j)=4
          ! else if (tau(j,ibox) .ge. 9.4
          ! &          .and. tau(j,ibox) .lt. 23.) then
          ! itau(j)=5
          ! else if (tau(j,ibox) .ge. 23.
          ! &          .and. tau(j,ibox) .lt. 60.) then
          ! itau(j)=6
          ! else if (tau(j,ibox) .ge. 60.) then
          ! itau(j)=7
          ! end if

          ! !determine cloud top pressure category
          ! if (    ptop(j,ibox) .gt. 0.
          ! &          .and.ptop(j,ibox) .lt. 180.) then
          ! ipres(j)=1
          ! else if(ptop(j,ibox) .ge. 180.
          ! &          .and.ptop(j,ibox) .lt. 310.) then
          ! ipres(j)=2
          ! else if(ptop(j,ibox) .ge. 310.
          ! &          .and.ptop(j,ibox) .lt. 440.) then
          ! ipres(j)=3
          ! else if(ptop(j,ibox) .ge. 440.
          ! &          .and.ptop(j,ibox) .lt. 560.) then
          ! ipres(j)=4
          ! else if(ptop(j,ibox) .ge. 560.
          ! &          .and.ptop(j,ibox) .lt. 680.) then
          ! ipres(j)=5
          ! else if(ptop(j,ibox) .ge. 680.
          ! &          .and.ptop(j,ibox) .lt. 800.) then
          ! ipres(j)=6
          ! else if(ptop(j,ibox) .ge. 800.) then
          ! ipres(j)=7
          ! end if

            !update frequencies
          IF (ipres(j)>0 .AND. itau(j)>0) THEN
            fq_isccp(j, itau(j), ipres(j)) = fq_isccp(j, itau(j), ipres(j)) + &
              boxarea
          END IF

          ! IM calcul stats regime dynamique BEG
          ! iw(j) = int((w(j)-wmin)/pas_w) +1
          ! pctj(itau(j),ipres(j),iw(j))=.FALSE.
          ! !update frequencies W500
          ! if (pct_ocean(j)) then
          ! if (ipres(j) .gt. 0.and.itau(j) .gt. 0) then
          ! if (iw(j) .gt. int(wmin).and.iw(j) .le. iwmx) then
          ! print*,' ISCCP iw=',iw(j),j
          ! fq_dynreg(itau(j),ipres(j),iw(j))=
          ! &          fq_dynreg(itau(j),ipres(j),iw(j))+
          ! &          boxarea
          ! &          fq_isccp(j,itau(j),ipres(j))
          ! pctj(itau(j),ipres(j),iw(j))=.TRUE.
          ! nfq_dynreg(itau(j),ipres(j),iw(j))=
          ! &          nfq_dynreg(itau(j),ipres(j),iw(j))+1.
          ! end if
          ! end if
          ! end if
          ! IM calcul stats regime dynamique END
        END IF !IM boxcloudy

      END IF !IM sunlit

      ! IM
      ! CALL CPU_time(t2)
      ! print*,'IF sunlit boxcloudy t2 - t1',t2 - t1
      ! IM
    END DO !IM ibox/j


    ! IM ajout stats s/ W500 BEG
    ! IM ajout stats s/ W500 END

    ! if(j.EQ.6722) then
    ! print*,' ISCCP',w(j),iw(j),ipres(j),itau(j)
    ! endif

    ! if (pct_ocean(j)) then
    ! if (ipres(j) .gt. 0.and.itau(j) .gt. 0) then
    ! if (iw(j) .gt. int(wmin).and.iw(j) .le. iwmx) then
    ! if(pctj(itau(j),ipres(j),iw(j))) THEN
    ! nfq_dynreg(itau(j),ipres(j),iw(j))=
    ! &    nfq_dynreg(itau(j),ipres(j),iw(j))+1.
    ! if(itau(j).EQ.4.AND.ipres(j).EQ.2.AND.
    ! &    iw(j).EQ.10) then
    ! PRINT*,' isccp AVANT',
    ! &    nfq_dynreg(itau(j),ipres(j),iw(j)),
    ! &    fq_dynreg(itau(j),ipres(j),iw(j))
    ! endif
    ! endif
    ! endif
    ! endif
    ! endif

  END DO
    !compute mean cloud properties
  ! IM j/ibox
  DO j = 1, npoints
    IF (totalcldarea(j)>0.) THEN
      meanptop(j) = meanptop(j)/totalcldarea(j)
      IF (sunlit(j)==1) THEN
        meanalbedocld(j) = meanalbedocld(j)/totalcldarea(j)
        meantaucld(j) = tautab(min(255,max(1,nint(meanalbedocld(j)))))
      END IF
    END IF
  END DO ! j

  ! IM ajout stats s/ W500 BEG
  ! do nw = 1, iwmx
  ! do l = 1, 7
  ! do k = 1, 7
  ! if (nfq_dynreg(k,l,nw).GT.0.) then
  ! fq_dynreg(k,l,nw) = fq_dynreg(k,l,nw)/nfq_dynreg(k,l,nw)
  ! if(k.EQ.4.AND.l.EQ.2.AND.nw.EQ.10) then
  ! print*,' isccp APRES',nfq_dynreg(k,l,nw),
  ! &   fq_dynreg(k,l,nw)
  ! endif
  ! else
  ! if(fq_dynreg(k,l,nw).NE.0.) then
  ! print*,'nfq_dynreg = 0 tau,pc,nw',k,l,nw,fq_dynreg(k,l,nw)
  ! endif
  ! fq_dynreg(k,l,nw) = -1.E+20
  ! nfq_dynreg(k,l,nw) = 1.E+20
  ! end if
  ! enddo !k
  ! enddo !l
  ! enddo !nw
  ! IM ajout stats s/ W500 END
  ! ---------------------------------------------------!

  ! ---------------------------------------------------!
  ! OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM

  ! cIM
  ! if (debugcol.ne.0) then

  ! do j=1,npoints,debugcol

  ! !produce character output
  ! do ilev=1,nlev
  ! do ibox=1,ncol
  ! acc(ilev,ibox)=0
  ! enddo
  ! enddo

  ! do ilev=1,nlev
  ! do ibox=1,ncol
  ! acc(ilev,ibox)=frac_out(j,ibox,ilev)*2
  ! if (levmatch(j,ibox) .eq. ilev)
  ! &                 acc(ilev,ibox)=acc(ilev,ibox)+1
  ! enddo
  ! enddo

    !print test

  ! write(ftn09,11) j
  ! 11        format('ftn09.',i4.4)
  ! open(9, FILE=ftn09, FORM='FORMATTED')

  ! write(9,'(a1)') ' '
  ! write(9,'(10i5)')
  ! &                  (ilev,ilev=5,nlev,5)
  ! write(9,'(a1)') ' '

  ! do ibox=1,ncol
  ! write(9,'(40(a1),1x,40(a1))')
  ! &           (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev)
  ! &           ,(cchar(acc(ilev,ibox)+1),ilev=1,nlev)
  ! end do
  ! close(9)

  ! IM
  ! if (ncolprint.ne.0) then
  ! write(6,'(a1)') ' '
  ! write(6,'(a2,1X,5(a7,1X),a50)')
  ! &                  'ilev',
  ! &                  'pfull','at',
  ! &                  'cc*100','dem_s','dtau_s',
  ! &                  'cchar'

  ! do 4012 ilev=1,nlev
  ! write(6,'(60i2)') (box(i,ilev),i=1,ncolprint)
  ! write(6,'(i2,1X,5(f7.2,1X),50(a1))')
  ! &                  ilev,
  ! &                  pfull(j,ilev)/100.,at(j,ilev),
  ! &                  cc(j,ilev)*100.0,dem_s(j,ilev),dtau_s(j,ilev)
  ! &                  ,(cchar(acc(ilev,ibox)+1),ibox=1,ncolprint)
  ! 4012           continue
  ! write (6,'(a)') 'skt(j):'
  ! write (6,'(8f7.2)') skt(j)

  ! write (6,'(8I7)') (ibox,ibox=1,ncolprint)

  ! write (6,'(a)') 'tau:'
  ! write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)

  ! write (6,'(a)') 'tb:'
  ! write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)

  ! write (6,'(a)') 'ptop:'
  ! write (6,'(8f7.2)') (ptop(j,ibox),ibox=1,ncolprint)
  ! endif

  ! enddo

  ! end if

  RETURN
END SUBROUTINE isccp_cloud_types
