!
! $Id:  ??? $
!
MODULE cmp_seri_mod
  !======================================================================
  ! to check if the x_seri vaariables have changed
  !
  !======================================================================

  REAL,SAVE,ALLOCATABLE ::  u_prev(:,:), v_prev(:,:), ql_prev(:,:), &
       qs_prev(:,:), q_prev(:,:), t_prev(:,:)
  !$OMP THREADPRIVATE(u_prev, v_prev, ql_prev, qs_prev, q_prev, t_prev)

CONTAINS

  SUBROUTINE init_cmp_seri
    !======================================================================
    !
    !======================================================================

    !======================================================================
    ! Declarations
    !======================================================================
    USE dimphy, ONLY: klon, klev
    USE phys_local_var_mod, ONLY: u_seri, v_seri, ql_seri, qs_seri, q_seri, &
         t_seri
    !
    IMPLICIT none
    ! Local :
    !--------
    logical, save :: first=.true.
    !$OMP THREADPRIVATE(first)

    !
    !======================================================================
    !
    if (first) then
       ALLOCATE(u_prev(klon,klev))
       ALLOCATE(v_prev(klon,klev))
       ALLOCATE(ql_prev(klon,klev))
       ALLOCATE(qs_prev(klon,klev))
       ALLOCATE(q_prev(klon,klev))
       ALLOCATE(t_prev(klon,klev))
       first=.false.
    endif
    ! store current values of x_seri variables
    u_prev(:,:) =u_seri(:,:)
    v_prev(:,:) =v_seri(:,:)
    ql_prev(:,:)=ql_seri(:,:)
    qs_prev(:,:)=qs_seri(:,:)
    t_prev(:,:) =t_seri(:,:)
    q_prev(:,:) =q_seri(:,:)

  END SUBROUTINE init_cmp_seri

  !======================================================================
  !
  !======================================================================
  SUBROUTINE cmp_seri(iflag,text)
    !
    USE dimphy, ONLY: klon, klev
    USE phys_local_var_mod, ONLY: u_seri, v_seri, ql_seri, qs_seri, q_seri, &
         t_seri
    USE print_control_mod, ONLY: prt_level
    IMPLICIT none
    ! Arguments :
    !------------
    integer iflag
    CHARACTER*(*) text
    ! Local :
    !--------
    integer i, k
    real du_seri(klon,klev), dv_seri(klon,klev), dql_seri(klon,klev), &
         dqs_seri(klon,klev), dq_seri(klon,klev), dt_seri(klon,klev)
    logical t_differ, u_differ, v_differ, ql_differ, qs_differ, q_differ
    logical any_differ
    !
    any_differ = .false.
    !
    dt_seri(:,:)=t_seri(:,:)-t_prev(:,:)
    IF (any(t_seri/=t_prev)) then
       if ( .not. any_differ ) write (*,*) "In cmp_seri, after ",text
       write (*,*) '*** t_seri has changed',minval(dt_seri),maxval(dt_seri)
       t_differ=.true.
       any_differ = .true.
    ELSE
       t_differ=.false.
    END IF
    !
    dq_seri(:,:)=q_seri(:,:)-q_prev(:,:)
    IF (any(q_seri/=q_prev)) then
       if ( .not. any_differ ) write (*,*) "In cmp_seri, after ",text
       write (*,*) '*** q_seri has changed',minval(dq_seri),maxval(dq_seri)
       q_differ=.true.
       any_differ = .true.
    ELSE
       q_differ=.false.
    END IF
    !
    dql_seri(:,:)=ql_seri(:,:)-ql_prev(:,:)
    IF (any(ql_seri/=ql_prev)) then
       if ( .not. any_differ ) write (*,*) "In cmp_seri, after ",text
       write (*,*) '*** ql_seri has changed',minval(dql_seri),maxval(dql_seri)
       ql_differ=.true.
       any_differ = .true.
    ELSE
       ql_differ=.false.
    END IF
    !
    dqs_seri(:,:)=qs_seri(:,:)-qs_prev(:,:)
    IF (any(qs_seri/=qs_prev)) then
       if ( .not. any_differ ) write (*,*) "In cmp_seri, after ",text
       write (*,*) '*** qs_seri has changed',minval(dqs_seri),maxval(dqs_seri)
       qs_differ=.true.
       any_differ = .true.
    ELSE
       qs_differ=.false.
    END IF
    !
    du_seri(:,:)=u_seri(:,:)-u_prev(:,:)
    IF (any(u_seri/=u_prev)) then
       if ( .not. any_differ ) write (*,*) "In cmp_seri, after ",text
       write (*,*) '*** u_seri has changed',minval(du_seri),maxval(du_seri)
       u_differ=.true.
       any_differ = .true.
    ELSE
       u_differ=.false.
    END IF
    !
    dv_seri(:,:)=v_seri(:,:)-v_prev(:,:)
    IF (any(v_seri/=v_prev)) then
       if ( .not. any_differ ) write (*,*) "In cmp_seri, after ",text
       write (*,*) '*** v_seri has changed',minval(dv_seri),maxval(dv_seri)
       v_differ=.true.
       any_differ = .true.
    ELSE
       v_differ=.false.
    END IF
    !
    if (any_differ) then
       IF (iflag >= 2 .and. prt_level >= 5) then
          i=1
          DO k = 1, klev
             write(*,'(i3,6e14.6)')k, t_seri(i,k)-t_prev(i,k), u_seri(i,k)-u_prev(i,k), &
                  v_seri(i,k)-v_prev(i,k), ql_seri(i,k)-ql_prev(i,k),&
                  qs_seri(i,k)-qs_prev(i,k), q_seri(i,k)-q_prev(i,k)
          ENDDO
       END IF
    else
       write (*,*) 'No change for any xxx_seri after ',text
    end if

  END SUBROUTINE cmp_seri

END MODULE cmp_seri_mod

