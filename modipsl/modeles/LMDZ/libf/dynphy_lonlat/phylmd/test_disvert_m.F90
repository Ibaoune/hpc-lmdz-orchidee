module test_disvert_m

  implicit none

contains

  subroutine test_disvert

    ! Author: Lionel GUEZ

    ! This procedure tests the order of pressure values at half-levels
    ! and full levels. We arbitrarily choose to test ngrid values of
    ! the surface pressure, which sample possible values on Earth.

    use exner_hyb_m, only: exner_hyb
    use vertical_layers_mod, only: ap,bp,preff
    use comconst_mod, only: kappa, cpp

    ! For llm:
    include "dimensions.h"

    ! Local:
    integer l, i
    integer, parameter:: ngrid = 7
    real p(ngrid, llm + 1) ! pressure at half-level, in Pa
    real pks(ngrid) ! exner function at the surface, in J K-1 kg-1 
    real pk(ngrid, llm) ! exner function at full level, in J K-1 kg-1 
    real ps(ngrid) ! surface pressure, in Pa
    real p_lay(ngrid, llm) ! pressure at full level, in Pa
    real delta_ps ! in Pa

    !---------------------

    print *, "Call sequence information: test_disvert"

    delta_ps = 6e4 / (ngrid - 1)
    ps = (/(5e4 + delta_ps * i, i = 0, ngrid - 1)/)
    forall (l = 1: llm + 1) p(:, l) = ap(l) + bp(l) * ps
    call exner_hyb(ngrid, ps, p, pks, pk)
    p_lay = preff * (pk / cpp)**(1. / kappa)

    ! Are pressure values in the right order?
    if (any(p(:, :llm) <= p_lay .or. p_lay <= p(:, 2:))) then
       ! List details and stop:
       do l = 1, llm
          do i = 1, ngrid
             if (p(i, l) <= p_lay(i, l)) then
                print 1000, "ps = ", ps(i) / 100., "hPa, p(level ",  l, &
                     ") = ", p(i, l) / 100., " hPa <= p(layer ", l, ") = ", &
                     p_lay(i, l) / 100., " hPa"
             end if
             if (p_lay(i, l) <= p(i, l + 1)) then
                print 1000, "ps = ", ps(i) / 100., "hPa, p(layer ", l, ") = ", &
                     p_lay(i, l) / 100., " hPa <= p(level ", l + 1, ") = ", &
                     p(i, l + 1) / 100., " hPa"
             end if
          end do
       end do
       call abort_physic("test_disvert", "bad order of pressure values", 1)
    end if

1000 format (3(a, g10.4, a, i0))

  end subroutine test_disvert

end module test_disvert_m
