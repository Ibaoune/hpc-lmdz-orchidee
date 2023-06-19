module ifirstloc_m

  implicit none

contains

  INTEGER FUNCTION ifirstloc(mask)

    ! Location of first true value in a logical array, returned as an
    ! integer. Returns size(mask)+1 if mask has zero element or all
    ! elements of mask are false. So the result is always >= 1.

    ! See notes on programming choices.

    LOGICAL, INTENT(IN):: mask(:)

    ! Local:
    integer n

    !-------------------------------------------------------

    n = size(mask)
    ifirstloc = 1

    if (n >= 1) then
       do while (ifirstloc <= n - 1 .and. .not. mask(ifirstloc))
          ifirstloc = ifirstloc + 1
       end do
       ! {1 <= ifirstloc <= n}
       if (.not. mask(ifirstloc)) ifirstloc = n + 1
    end if

  END FUNCTION ifirstloc

end module ifirstloc_m
